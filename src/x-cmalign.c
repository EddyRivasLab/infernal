/* cmalign: align sequences to a CM.
 * 
 * EPN, Fri Dec 30 10:13:34 2011 [Updated for v1.1]
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"	

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "easel.h"		/* general seq analysis library   */
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile2.h"
#include "esl_msaweight.h"
#include "esl_random.h"		
#include "esl_sq.h"		
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_sse.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif /*HAVE_MPI*/

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define OUTALPHOPTS  "--rna,--dna"                               /* Exclusive choice for output alphabet */

#define CMALIGN_MAX_NSEQ_ILEAVED  100000     /* 100k sequences, most we allow to be aligned and output in interleaved format */
#define CMALIGN_MAX_NRES_ILEAVED  100000000  /* 100Mb, most we allow to be aligned and output in interleaved format */
/* What I will eventually set these as: */
///#define CMALIGN_MAX_NRES 10000000  /* 10 Mb, average parsetree is 25 bytes/position this means ~250Mb for all parsetrees */

#define CMALIGN_MAX_NRES 1000

#define DEBUGSERIAL 0
#define DEBUGMPI    0

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  CM_t             *cm;     /* a covariance model */
  CM_ALNDATA      **dataA;  /* array of CM_ALNDATA objects with ptrs to sqs, parsetrees, scores */
  int               n;      /* size of outdataA   */
  float             mxsize; /* max size in Mb of allowable DP mx */
  ESL_STOPWATCH    *w;      /* stopwatch for timing stages (band calc, alignment) */
  ESL_STOPWATCH    *w_tot;  /* stopwatch for timing total time for processing 1 seq */
} WORKER_INFO;

#define ALGOPTS      "--cyk,--optacc,--sample"         /* Exclusive choice for algorithm */
#define BIGALGOPTS   "--cyk,--optacc,--sample,--small" /* Incompatible with --optacc,--sample (except their selves) */
#if defined (HMMER_THREADS) && defined (HAVE_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

static ESL_OPTIONS options[] = {
  /* name                  type   default   env        range      toggles      reqs       incomp  help  docgroup*/
  { "-h",           eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "show brief help on version and usage",   1 },
  { "-o",        eslARG_OUTFILE,     NULL,  NULL,      NULL,      NULL,        NULL,        NULL, "output the alignment to file <f>, not stdout", 1 },
  { "-g",           eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "configure CM for glocal alignment [default: local]", 1 },
  { "-i",           eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "output in interleaved format (WARNING: memory intensive for large inputs)", 1 },
#ifdef HMMER_THREADS 
  { "--cpu",        eslARG_INT,      NULL,"HMMER_NCPU", "n>=0",   NULL,        NULL,     CPUOPTS, "number of parallel CPU workers to use for multithreads", 1 },
#endif
#ifdef HAVE_MPI
  { "--mpi",        eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,     MPIOPTS, "run as an MPI parallel program", 1 },  
#endif
  { "--devhelp",    eslARG_NONE,     NULL,  NULL,      NULL,      NULL,        NULL,        NULL, "show list of undocumented developer options", 1 },
  /* options controlling the alignment algorithm */
  { "--optacc",     eslARG_NONE,"default",  NULL,      NULL,"--optacc",        NULL,  BIGALGOPTS, "align with the Holmes/Durbin optimal accuracy algorithm", 2 },
  { "--cyk",        eslARG_NONE,    FALSE,  NULL,      NULL,"--optacc",        NULL,     ALGOPTS, "align with the CYK algorithm", 2 },
  { "--sample",     eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,     ALGOPTS, "sample alignment of each seq from posterior distribution", 2 },
  ///  { "--viterbi", eslARG_NONE,   FALSE, NULL, NULL,"--optacc","--no-prob",   ALGOPTS, "align to a CM Plan 9 HMM with the Viterbi algorithm",2 },
  { "--sub",        eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,"-g,--notrunc",      NULL, "build sub CM for columns b/t HMM predicted start/end points", 2 },
  { "--small",      eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,"--cyk,--noprob,--nonbanded",NULL, "use small memory divide and conquer (d&c) algorithm", 2 },
  { "--notrunc",    eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "do not use truncated alignment algorithm", 2 },
  { "--nonbanded",  eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "do not use bands to accelerate aln algorithm", 2 },
  { "--tau",        eslARG_REAL,   "1E-7",  NULL,"1E-18<x<1",     NULL,        NULL,"--nonbanded","set tail loss prob for HMM bands to <x>", 2 },
  { "--mxsize",     eslARG_REAL,  "256.0",  NULL,   "x>0.",       NULL,        NULL,   "--small", "set maximum allowable DP matrix size to <x> Mb", 2 },
  { "--seed",        eslARG_INT,    "181",  NULL,   "n>=0",       NULL,  "--sample",        NULL, "w/--sample, set RNG seed to <n> (if 0: one-time arbitrary seed)", 2 },
  /* options for including a preset alignment */
  { "--mapali",   eslARG_INFILE,     NULL,  NULL,      NULL,      NULL,      NULL,          NULL, "include alignment in file <f> (same ali that CM came from)",       2 },
  { "--mapstr",     eslARG_NONE,     NULL,  NULL,      NULL,      NULL, "--mapali",         NULL, "include structure (w/pknots) from <f> from --mapali <f>", 4 },
  /* options that modify how the output alignment is created */
  { "--dna",        eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 3 },
  { "--oneline",    eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "output in non-interleaved (1 line/seq) format ", 3 },
  { "--noannot",    eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "do not add cmalign execution annotation to the alignment", 3 },
  { "--noprob",     eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "do not include posterior probabilities in the alignment", 3 },
  { "--matchonly",  eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "include only match columns in output alignment", 3 },
  /* options controlling optional output */
  { "--tfile",   eslARG_OUTFILE,     NULL,  NULL,      NULL,      NULL,        NULL,        NULL, "dump individual sequence parsetrees to file <f>", 4 },
  { "--ifile",   eslARG_OUTFILE,     NULL,  NULL,      NULL,      NULL,        NULL,        NULL, "dump information on per-sequence inserts to file <f>", 4 },
  { "--elfile",  eslARG_OUTFILE,     NULL,  NULL,      NULL,      NULL,        NULL,        NULL, "dump information on per-sequence EL inserts to file <f>", 4 },
  { "--sfile",   eslARG_OUTFILE,     NULL,  NULL,      NULL,      NULL,        NULL,        NULL, "dump alignment score information to file <f>", 4 },
  /* developer options the average user doesn't need to know about (only shown if --devhelp) */
  { "--regress", eslARG_OUTFILE,     NULL,  NULL,      NULL,      NULL,        "-i",  "--mapali", "save regression test data to file <f>", 101}, 
#ifdef HAVE_MPI
  { "--stall",      eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "arrest after start: for debugging MPI under gdb", 101 },  
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char            *cmfile;	       /* name of input CM file  */ 
  char            *sqfile;	       /* name of sequence file  */ 
  CM_FILE         *cmfp;	       /* open input CM file stream       */
  ESL_SQFILE      *sqfp;               /* open sequence input file stream */
  ESL_ALPHABET    *abc;                /* alphabet for input */
  ESL_ALPHABET    *abc_out;            /* alphabet for output */

  /* mpi */
  int              do_mpi;
  int              my_rank;
  int              nproc;
  int              do_stall;             /* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  FILE            *tmpfp;		 /* the temporary output file where alignments are initially written */
  FILE            *ofp;		         /* output file where alignments are ultimately written (default is stdout) */
  FILE            *tfp;       	         /* optional output for parsetrees  */
  FILE            *ifp;	                 /* optional output for insert info */
  FILE            *efp;	                 /* optional output for EL insert info */
  FILE            *sfp;                  /* optional output for alignment scores */
  FILE            *rfp;       	         /* optional output for --regress alignment */
};

static char usage[]  = "[-options] <cmfile> <sequence file>";
static char banner[] = "align sequences to a CM";

static void serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, char *errbuf, ESL_SQ_BLOCK **sq_blockA, int nblocks);

/* Functions to avoid code duplication for common tasks */
static int  init_master_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  init_shared_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  output_alignment(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, FILE *ofp, CM_ALNDATA **dataA, int ndata, char *map_sscons);
static void output_info_file_header(FILE *fp, char *firstline, char *elstring);
static int  output_scores(FILE *ofp, CM_t *cm, char *errbuf, CM_ALNDATA **dataA, int ndata, int first_idx);
static int  output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile, char *sqfile, CM_t *cm);
static void process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile, char **ret_sqfile);
static int  map_alignment(const char *msafile, CM_t *cm, char *errbuf, CM_ALNDATA ***ret_dataA, int *ret_ndata, char **ret_ss);

/* Functions that enable memory efficiency by storing only a fraction
 * of the seqs/parsetrees from target file in memory at once.
 */
static int  create_and_output_final_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nali, char *tmpfile);
static void update_maxins_and_maxel(ESL_MSA *msa, int clen, int64_t alen, int *maxins, int *maxel);
static int  determine_gap_columns_to_add(ESL_MSA *msa, int *maxins, int *maxel, int clen, int **ret_ngap_insA, int **ret_ngap_elA, int **ret_ngap_eitherA, char *errbuf);
static void inflate_gc_with_gaps_and_els(FILE *ofp, ESL_MSA *msa, int *ngap_insA, int *ngap_elA, char **ret_ss_cons2print, char **ret_rf2print);

#ifdef HMMER_THREADS
static int  thread_loop(WORKER_INFO *info, char *errbuf, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ_BLOCK **sq_blockA, int nblocks);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

#ifdef HAVE_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif /*HAVE_MPI*/

static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile, char **ret_sqfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     cm_Fail("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_GetBoolean(go, "--devhelp")) { 
    cm_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);
    puts("\nBasic options:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 100); /* 1= group; 2 = indentation; 100=textwidth*/
    puts("\nOptions controlling alignment algorithm:");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 100); 
    puts("\nOptional output files:");
    esl_opt_DisplayHelp(stdout, go, 3, 2, 100); 
    if(esl_opt_GetBoolean(go, "--devhelp")) { 
      puts("\nUndocumented developer options:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 100); 
    }   
    exit(0);
  } 

  if (esl_opt_ArgNumber(go)                 != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_cmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <cmfile> argument on command line");  goto ERROR; }
  if ((*ret_sqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  
  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}

/* output_header(): 
 *
 * Differently from other Infernal applications, which output header to
 * stdout, we output the header to stdout only if the user has
 * specified a non-stdout output file for the alignment. Otherwise,
 * the alignment will be printed to stdout, without a header because
 * we want it to be a valid Stockholm format alignment.
 */

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile, char *sqfile, CM_t *cm)
{
  cm_banner(ofp, go->argv[0], banner);
                                               fprintf(ofp, "# CM file:                                     %s\n", cmfile);
					       fprintf(ofp, "# sequence file:                               %s\n", sqfile);
                                               fprintf(ofp, "# CM name:                                     %s\n", cm->name);
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu"))       {     fprintf(ofp, "# number of worker threads:                    %d\n", esl_opt_GetInteger(go, "--cpu")); }
#endif
#ifdef HAVE_MPI
  if (esl_opt_IsUsed(go, "--mpi"))       {     fprintf(ofp, "# MPI:                                         on\n"); }
#endif
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  return eslOK;
}

int
main(int argc, char **argv)
{
  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
  struct cfg_s     cfg;           /* configuration data                      */
  
  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet)
   */
  cfg.cmfile       = NULL;
  cfg.sqfile       = NULL;
  cfg.cmfp         = NULL; 
  cfg.sqfp         = NULL; 
  cfg.do_mpi       = FALSE;               /* this gets reset below, if we init MPI */
  cfg.nproc        = 0;                   /* this gets reset below, if we init MPI */
  cfg.my_rank      = 0;                   /* this gets reset below, if we init MPI */
  cfg.abc          = NULL; 

  cfg.tmpfp        = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ofp          = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tfp          = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ifp          = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.efp          = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.sfp          = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.rfp          = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */

  /* Initializations */
  init_ilogsum();
  FLogsumInit();
  process_commandline(argc, argv, &go, &(cfg.cmfile), &(cfg.sqfile));

  /* update cfg now that we have go */
  cfg.abc_out = esl_opt_GetBoolean(go, "--dna") ? esl_alphabet_Create(eslDNA) : esl_alphabet_Create(eslRNA);

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI

  /* TEMP */
  pid_t pid;
  /* get the process id */
  pid = getpid();
  printf("The process id is %d\n", pid);
  fflush(stdout);
  /* TEMP */

  /* pause the execution of the programs execution until the user has a
   * chance to attach with a debugger and send a signal to resume execution
   * i.e. (gdb) signal SIGCONT
   */
  if (esl_opt_GetBoolean(go, "--stall")) pause();

  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you execute mpirun?)");

      if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
      else 		    status = mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      serial_master(go, &cfg);
    }

  /* Close output files */
  if(cfg.ofp   != NULL && esl_opt_IsUsed(go, "-o")) fclose(cfg.ofp); 
  if(cfg.tmpfp != NULL) fclose(cfg.tmpfp); 
  if(cfg.tfp   != NULL) fclose(cfg.tfp); 
  if(cfg.ifp   != NULL) fclose(cfg.ifp); 
  if(cfg.efp   != NULL) fclose(cfg.efp); 
  if(cfg.sfp   != NULL) fclose(cfg.sfp); 
  if(cfg.cmfp  != NULL) cm_file_Close(cfg.cmfp);
  if(cfg.sqfp  != NULL) esl_sqfile_Close(cfg.sqfp);

  esl_getopts_Destroy(go);

  return status;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 *
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errbuf) by the ESL_FAIL mech.
 */
static int
init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  
  /* initialize cfg variables used by masters and workers */
  if((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) return status;

  /* open output files */
  cfg->ofp = stdout;
  if (esl_opt_IsUsed(go, "-o")) { 
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)      ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } 
  if (esl_opt_IsUsed(go, "--tfile")) { 
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
  }
  if (esl_opt_IsUsed(go, "--ifile")) { 
    if ((cfg->ifp = fopen(esl_opt_GetString(go, "--ifile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --ifile output file %s\n", esl_opt_GetString(go, "--ifile"));
    output_info_file_header(cfg->ifp, "Insert information file created by cmalign.", "");
  }
  if (esl_opt_IsUsed(go, "--elfile")) { 
    if ((cfg->efp = fopen(esl_opt_GetString(go, "--efile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --elfile output file %s\n", esl_opt_GetString(go, "--elfile"));
    output_info_file_header(cfg->efp, "EL state (local end) insert information file created by cmalign.", "EL ");
  }
  if (esl_opt_IsUsed(go, "--sfile")) { 
    if ((cfg->sfp = fopen(esl_opt_GetString(go, "--sfile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --sfile output file %s\n", esl_opt_GetString(go, "--sfile"));
  }
  if (esl_opt_IsUsed(go, "--regress")) { 
    if ((cfg->rfp = fopen(esl_opt_GetString(go, "--regress"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
  }
  return eslOK;
}


/* init_shared_cfg()
 * Called by serial masters and mpi workers and masters.
 *
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errbuf) by the ESL_FAIL mech.
 */
static int
init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;

  /* open CM file */
  status = cm_file_Open(cfg->cmfile, NULL, FALSE, &(cfg->cmfp), errbuf);
  if      (status == eslENOTFOUND) return status;
  else if (status == eslEFORMAT)   return status;
  else if (status != eslOK)        return status;

  /* open sequence file */
  status = esl_sqfile_Open(cfg->sqfile, eslSQFILE_UNKNOWN, p7_SEQDBENV, &(cfg->sqfp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open sequence file %s for reading\n",          cfg->sqfile);
  else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Sequence file %s is empty or misformatted\n",            cfg->sqfile);
  else if (status == eslEINVAL)    ESL_FAIL(status, errbuf, "Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d opening sequence file %s\n", status, cfg->sqfile);  

  return eslOK;
}

/* serial_master()
 * The serial version of cmalign.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;                   /* Easel status */
  char     errbuf[eslERRBUFSIZE];    /* for printing error messages */
  CM_t    *cm = NULL;                /* a CM */
  int      i, j, k;                  /* counter over blocks, parsetrees and workers */
  int      nali;                     /* index of the (possibly temporary) alignment we are working on */
  int      nseq_cur;                 /* number of sequences in current alignment */
  int      nseq_tot;                 /* number of sequences in all alignments */

  /* variables related to output, we may use a tmpfile if seqfile is large */
  int      do_ileaved;               /* TRUE to output interleaved alignment */ 
  int      use_tmpfile;              /* print out current alignment to tmpfile? */
  int      created_tmpfile = FALSE;  /* TRUE if we've created a tmp file for current CM */
  char tmpfile[32] = "esltmpXXXXXX"; /* name of the tmpfile */
  CM_ALNDATA **merged_dataA = NULL; /* array of all CM_ALNDATA pointers for current alignment */

  /* variables related to reading sequence blocks */
  int            sstatus = eslOK;  /* status from esl_sq_ReadBlock() */
  ESL_SQ_BLOCK  *big_sq_block;     /* a large sequence block, this gets split into sq_blockA */
  ESL_SQ_BLOCK  *nxt_big_sq_block; /* a large sequence block */
  ESL_SQ_BLOCK **sq_blockA;        /* array of seq blocks, these point at seqs in big_sq_block */
  int            target_nblocks;   /* target number of blocks we'll split big_sq_block into */
  int            nblocks;          /* actual number of blocks we split big_sq_block into */
  int            block_max_nres;   /* maximum number of residues we'll allow in a sq_blockA[] block */
  int            block_max_nseq;   /* maximum number of sequences we'll allow in a sq_blockA[] block */
  int            nres_big_block;   /* number of residues in big_sq_block */
  int            nres_this_block;  /* number of residues in sq_blockA[] block */
  int            nseq_this_block;  /* number of sequences in sq_blockA[] block */
  int     target_nres_per_block;   /* target number of residues per sq_blockA[] block */
  int            reached_eof;      /* TRUE if we've reached EOF in target sequence file */
  int            x;                /* counter over sequences in big_sq_block */
  ESL_SQ        *sqp = NULL;       /* pointer to a sequence in a block */

  /* variables related to --mapali */
  char         *map_file   = NULL; /* name of alignment file from --mapali */
  CM_ALNDATA **map_dataA  = NULL; /* array of CM_ALNDATA pointers for mapali alignment */
  int           nmap_data  = 0;    /* number of CM_ALNDATA ptrs in map_dataA */
  int           nmap_cur   = 0;    /* number of CM_ALNDATA ptrs to include in current iteration, 0 unless nali==0 */
  char         *map_sscons = NULL; /* SS_cons from mapali, only used if --mapstr */

  /* variables related to threaded implementation */
  int              ncpus = 0;        /* number of CPUs working */
  WORKER_INFO     *info      = NULL; /* the worker info */
  int              infocnt   = 0;    /* number of worker infos */
#ifdef HMMER_THREADS
  ESL_THREADS     *threadObj = NULL;
  ESL_WORK_QUEUE  *queue     = NULL;
#endif
  
  /* General notes on serial_master()'s strategy: 
   * 
   * Ideally, we'd read in all sequences, align them all, and output
   * the alignment. But we're worried about running out of memory, so
   * we read in sequence blocks (big_sq_block) of at most
   * CMALIGN_MAX_NRES (10 Mb) at a time, and process each in turn. (A
   * parsetree is about 25 bytes per residue, so that should be about
   * 250 Mb). If there's more than one such block, we output each to a
   * tmp file and free the parsetrees afterwards so we don't require
   * too much memory. Once finished, we go through the tmp file and
   * merge all the alignments within it into a single one (without
   * ever storing all of them simultaneously) and output it to the
   * standard output file. In this case we need to use Pfam (1
   * line/seq) format so we don't have to store the full alignment/set
   * of parsetrees at once. If there's only one block we just output
   * it to the standard output file in interleaved format (which is
   * what previous versions of cmalign did).
   * 
   * The code could be simplified if we read the sequences into 
   * a array of sequence blocks of preset size, instead of the two-step
   * process of into a big sequence block and then dividing it up.
   * The reason we do the latter is so we can optimize the distribution
   * of work amongst all available workers, regardless of the remaining 
   * amount of sequence data in the file. One undesirable side effect
   * is that we point our small blocks at data in the big blocks (to 
   * avoid copying the data) and so we have to be careful about how
   * we free the data in the blocks.
   * 
   * Exceptions: 
   *
   * if -i: the user wants interleaved output, which we'll do but 
   *        only if we can read the full sequence file into a single
   *        big_sq_block that doesn't exceed CMALIGN_MAX_NSEQ_ILEAVED
   *        (100K) sequences or CMALIGN_MAX_NRES_ILEAVED (100 Mb)
   *        residues. If it does, we die immediately once we realize
   *        it's too big, and tell the user. (I considered allowing
   *        arbitrarily large alignments with -i, but it would
   *        significantly complicate the code, and the user can still
   *        create such an alignment with esl-alimerge and/or
   *        esl-reformat).
   *
   * if --oneline: alignment will always be output in 1 line/seq format.
   */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  do_ileaved = esl_opt_GetBoolean(go, "-i") ? TRUE : FALSE;

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn  (go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                             esl_threads_CPUCount(&ncpus);
  printf("NCPUS: %d\n", ncpus);
  if (ncpus > 0) {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(WORKER_INFO) * infocnt);

  /* Read one CM, and make sure there's only one. This fills cfg->abc. */
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);
  if(status != eslOK) cm_Fail(cfg->cmfp->errbuf);
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), NULL);
  if(status != eslEOF) cm_Fail("HMM file %s does not contain just one CM\n", cfg->cmfp->fname);

  for (k = 0; k < infocnt; ++k)    {
    info[k].cm     = NULL;
    info[k].dataA  = NULL;
    info[k].n      = 0;
    info[k].mxsize = esl_opt_GetReal(go, "--mxsize");
    info[k].w      = esl_stopwatch_Create();
    info[k].w_tot   = esl_stopwatch_Create();
#ifdef HMMER_THREADS
    info[k].queue  = queue;
#endif
  }
  
#ifdef HMMER_THREADS    
  for (k = 0; k < ncpus * 2; ++k) {
    if((big_sq_block = esl_sq_CreateDigitalBlock(1, cfg->abc)) == NULL)   cm_Fail("Failed to allocate sequence block");
    if((status       = esl_workqueue_Init(queue, big_sq_block)) != eslOK) cm_Fail("Failed to add block to work queue");
  }
#endif

  /* initialization */
  nali = nseq_cur = nseq_tot = 0;
  if((status = initialize_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
  for (k = 0; k < infocnt; ++k) {
    if((status = cm_Clone(cm, errbuf, &(info[k].cm))) != eslOK) cm_Fail(errbuf);
  }

  /* include the mapali, if nec */
  if((map_file = esl_opt_GetString(go, "--mapali")) != NULL) { 
    if((status = map_alignment(map_file, cm, errbuf, &map_dataA, &nmap_data, &map_sscons)) != eslOK) cm_Fail(errbuf);
    if(esl_opt_GetBoolean(go, "--mapstr") && map_sscons == NULL) cm_Fail("Failed to read SS_cons for --mapstr from %s", map_file);
  }

  /* Our main loop will loop over reading a single large block
   * (<big_sq_block>) of sequences, up to CMALIGN_MAX_NRES residues or
   * <block_max_nseq> sequences, but potentially less if we reach the
   * end of the sequence file first. Then, based on <big_sq_block>'s
   * size, we divide that block up into smaller blocks. Each of these
   * smaller blocks will be a workunit for a worker.  If (do_ileaved)
   * we require a single big block of a different maximum size.
   */
  if(do_ileaved) { 
    block_max_nres = CMALIGN_MAX_NRES_ILEAVED;
    block_max_nseq = CMALIGN_MAX_NSEQ_ILEAVED;
  }
  else { 
    block_max_nres = CMALIGN_MAX_NRES;
    block_max_nseq = CMALIGN_MAX_NRES / ESL_MIN(cm->clen/2, 100); 
    /* We guess that most sequences will be 100 residues or larger, or
     * half of clen if clen < 200. In most cases they'll probably be
     * roughly clen, but we are specifically worried about aligning many
     * short truncated SSU reads. If sequences are larger,
     * CMALIGN_MAX_NRES will kick in before <block_max_nseq>.
     */
  }

  /* Read the first big block */
  big_sq_block = esl_sq_CreateDigitalBlock(block_max_nseq, cfg->abc);
  sstatus = esl_sqio_ReadBlock(cfg->sqfp, big_sq_block, block_max_nres, FALSE); /* FALSE says: read complete sequences */
  nxt_big_sq_block = big_sq_block; /* special case of first block read */

  while(sstatus == eslOK) { 
#ifdef HMMER_THREADS
    if (ncpus > 0) { 
      for (k = 0; k < infocnt; ++k) esl_threads_AddThread(threadObj, &info[k]);
    }
#endif
    /* First, read the next big sequence block, so we can determine if
     * we've reached the end of the seqfile. We need to know this for
     * two reasons:
     *
     * (1) if the first big block read above included all sequences,
     * we don't need to go into memory-saving mode and output to a
     * tmpfile, we can output (in interleaved mode) to the final 
     * output file. 
     *
     * (2) if -i used (do_ileaved = TRUE) we need to fail if we still
     * have sequences left, because the sequence file exceeded the
     * size limits. And we want to fail *before* we align all the 
     * sequences, so the user doesn't get pissed when the job fails
     * after seemingly going along fine for a while.
     */
    big_sq_block               = nxt_big_sq_block; /* our current big_sq_block becomes the one we read on the previous iteration */
    big_sq_block->first_seqidx = nseq_tot;
    nxt_big_sq_block = esl_sq_CreateDigitalBlock(block_max_nseq, cfg->abc);
    sstatus = esl_sqio_ReadBlock(cfg->sqfp, nxt_big_sq_block, block_max_nres, FALSE); /* FALSE says: read complete sequences */
    if(sstatus == eslEOF) reached_eof = TRUE; /* nxt_big_sq_block will be NULL */
    if(! reached_eof && do_ileaved) esl_fatal("Error: to use -i the sequence file must have less than %d sequences and %d residues", CMALIGN_MAX_NSEQ_ILEAVED, CMALIGN_MAX_NRES_ILEAVED);

    /* if we have workers, divide up big_sq_block into smaller blocks */
    if(infocnt == 1) { /* one worker, no need to divide up the block: sq_blockA will have 1 element: big_sq_block */
      target_nblocks = nblocks = 1;
      ESL_ALLOC(sq_blockA, sizeof(ESL_SQ_BLOCK *) * 1);
      sq_blockA[0] = big_sq_block;
    }
    else { /* we have multiple workers, divide up big_sq_block into smaller blocks */
      target_nblocks = ESL_MIN(10 * ncpus, big_sq_block->count); /* create 10 times as many blocks as workers, in case some blocks happen to take longer */ 
      /* determine how many residues each small block should have */
      nres_big_block = 0;
      for(i = 0; i < big_sq_block->count; i++) { 
	sqp = big_sq_block->list + i; 
	nres_big_block += sqp->n;
      }
      /* if there's fewer sequences than CPUs set target_nres to 1, this forces 1 seq/block */
      target_nres_per_block = (big_sq_block->count > ncpus) ? nres_big_block / target_nblocks : 1;

      /* create the small blocks and point them at the ESL_SQ elements in big_sq_block */
      ESL_ALLOC(sq_blockA, sizeof(ESL_SQ_BLOCK *) * target_nblocks);
      for(i = 0; i < target_nblocks; i++) sq_blockA[i] = NULL;
      nblocks = 0;
      x = 0; /* counter over sequences in big_sq_block */
      for(i = 0; i < target_nblocks; i++) { 
	nres_this_block = 0;
	nseq_this_block = 0;
	while(nres_this_block < target_nres_per_block && x < big_sq_block->count) { 
	  sqp = big_sq_block->list + x;
	  nres_this_block += sqp->n;
	  nseq_this_block++;
	  x++;
	}
	if(nseq_this_block > 0) { 
	  sq_blockA[i]               = esl_sq_CreateDigitalBlock(nseq_this_block, cfg->abc);
	  sq_blockA[i]->list         = big_sq_block->list + (x - nseq_this_block);
	  sq_blockA[i]->first_seqidx = nseq_tot + x - nseq_this_block;
	  sq_blockA[i]->count        = nseq_this_block;
	  sq_blockA[i]->complete     = TRUE;
	  nblocks++;
	  printf("created block: %2d/%2d with %5d sequences %8d/%8d residues\n", nblocks, target_nblocks, nseq_this_block, nres_this_block, target_nres_per_block);
	  /* be careful not to free big_sq_block->list, sq_blockA[i] is pointing at elements within it! */
	}
      }
    }
    nseq_cur  = big_sq_block->count;
    nseq_tot += big_sq_block->count;

    /* now we have an array of <nblocks> valid seq blocks, each is a workunit, align them all */
#ifdef HMMER_THREADS
    if (ncpus > 0)  status = thread_loop(info, errbuf, threadObj, queue, sq_blockA, nblocks);
    else            status = serial_loop(info, errbuf, sq_blockA, nblocks);
#else
    status = serial_loop(info, errbuf, sq_blockA, nblocks);
#endif
    if(status != eslOK) cm_Fail(errbuf);

    /* create a single array of all CM_ALNDATA objects, in original (input) order */
    nmap_cur = (nali == 0) ? nmap_data : 0;
    ESL_ALLOC(merged_dataA, sizeof(CM_ALNDATA *) * (nseq_cur + nmap_cur));
    /* prepend mapali data if nec */
    if(nmap_cur > 0) {
      for(j = 0; j < nmap_cur; j++) merged_dataA[j] = map_dataA[j];
      free(map_dataA); /* don't free the CM_ALNDATA objects, merged_dataA is pointing at them */
      map_dataA = NULL;
    }
    for(k = 0; k < infocnt; ++k) { 
      for(j = 0; j < info[k].n; j++) merged_dataA[info[k].dataA[j]->idx + nmap_cur] = info[k].dataA[j];
      free(info[k].dataA); /* don't free the CM_ALNDATA objects, merged_dataA is pointing at them */
      info[k].dataA = NULL;
      info[k].n     = 0;
    }

    /* output alignment (if do_ileaved we died above if we didn't reach EOF yet) */
    use_tmpfile = reached_eof ? FALSE : TRUE;
    if(use_tmpfile && (! created_tmpfile)) { 
      /* first aln for temporary output file, open the file */	
      if ((status = esl_tmpfile_named(tmpfile, &(cfg->tmpfp))) != eslOK) cm_Fail("Failed to open temporary output file (status %d)", status);
      created_tmpfile = TRUE;
    }
    if((status   = output_alignment(go, cfg, errbuf, cm, (use_tmpfile ? cfg->tmpfp : cfg->ofp), merged_dataA, nseq_cur + nmap_cur, map_sscons)) != eslOK) cm_Fail(errbuf);
    /* optionally output same alignment to regress file */
    if(cfg->rfp != NULL) { 
      if((status = output_alignment(go, cfg, errbuf, cm, cfg->rfp,                              merged_dataA, nseq_cur + nmap_cur, map_sscons)) != eslOK) cm_Fail(errbuf);
    }    
    nali++;

    /* output scores to stdout, if -o used */
    if(cfg->ofp != stdout) { 
      if(nali == 1) output_header(stdout, go, cfg->cmfile, cfg->sqfile, cm);
      if((status =  output_scores(stdout,   cm, errbuf, merged_dataA, nseq_tot, nmap_cur)) != eslOK) cm_Fail(errbuf);
    }
    /* output scores to scores file, if --sfp used */
    if(cfg->sfp != NULL) { 
      if(nali == 1) output_header(stdout, go, cfg->cmfile, cfg->sqfile, cm);
      if((status = output_scores(cfg->sfp, cm, errbuf, merged_dataA, nseq_tot, nmap_cur)) != eslOK) cm_Fail(errbuf);
    }

    /* free big block and worker data */
    esl_sq_DestroyBlock(big_sq_block);
    big_sq_block = NULL;
    if(infocnt > 1) { 
      for(i = 0; i < nblocks; i++) { 
	if(sq_blockA[i] != NULL) { 
	  /* remember: sequences inside sq_blockA[] were just pointing at big_sq_block's sequences */
	  free(sq_blockA[i]);
	}      
      }
      if(sq_blockA != NULL) free(sq_blockA);
      sq_blockA = NULL;
    }
  } /* end of outer while loop 'while(sstatus == eslOK)' */
  if     (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n", cfg->sqfp->filename, esl_sqfile_GetErrorBuf(cfg->sqfp));
  else if(sstatus == eslEMEM)    cm_Fail("Out of memory");
  else if(sstatus != eslEOF)     cm_Fail("Unexpected error while reading sequence file");
    
  /* if nec, close tmpfile then merge all alignments in it */
  if(created_tmpfile) { 
    fclose(cfg->tmpfp); /* we're done writing to tmpfp */
    cfg->tmpfp = NULL;
    /* merge all temporary alignments now in cfg->tmpfp, and output merged alignment */
    if((status = create_and_output_final_msa(go, cfg, errbuf, cm, nali, tmpfile)) != eslOK) cm_Fail(errbuf);
    remove(tmpfile); 
  }
    
  /* finish insert and el files */
  if(cfg->ifp != NULL) { fprintf(cfg->ifp, "//\n"); }
  if(cfg->efp != NULL) { fprintf(cfg->efp, "//\n"); }
      
  FreeCM(cm);

  return;
  
  ERROR:
  cm_Fail("Memory allocation error.");
  return;
}

/* serial_loop(): 
 * 
 * Align all sequences in array of sequence blocks (sq_blockA) 
 * and store parsetrees.
 */
static int
serial_loop(WORKER_INFO *info, char *errbuf, ESL_SQ_BLOCK **sq_blockA, int nblocks)
{
  int status;
  int i, j;                    /* counter over blocks, parsetrees */
  int cur_n = 0;               /* current number of parsetrees valid in info->trA */
  CM_ALNDATA **bdataA = NULL; /* array of CM_ALNDATA for current block */

  /* allocate dataA */
  info->n = 0;
  for(i = 0; i < nblocks; i++) { 
    if(sq_blockA[i] != NULL) info->n += sq_blockA[i]->count; 
  }
  ESL_ALLOC(info->dataA, sizeof(CM_ALNDATA *) * info->n);
  for(j = 0; j < info->n; j++) info->dataA[j] = NULL;

  for(i = 0; i < nblocks; i++) { 
    if(sq_blockA[i] != NULL && sq_blockA[i]->count > 0) { 

      if((status = ProcessAlignmentWorkunit(info->cm, errbuf, sq_blockA[i], info->mxsize, info->w, info->w_tot, &bdataA)) != eslOK) cm_Fail(errbuf);
      /* point info->dataA at newly collected data, including parsetrees */
      for(j = 0; j < sq_blockA[i]->count; j++) info->dataA[cur_n + j] = bdataA[j];
      cur_n += sq_blockA[i]->count;
      free(bdataA);  /* don't free CM_ALNDATA objects, info->dataA is pointing at them */
      bdataA = NULL;
    }
  }
      
  return eslOK;
  
 ERROR: 
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVERREACHED */
}
 
#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, char *errbuf, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ_BLOCK **sq_blockA, int nblocks)
{
  int           status = eslOK;
  int           i, k;            /* counter over blocks, workers */
  ESL_SQ_BLOCK *sq_block;
  void         *new_sq_block;
  ESL_SQ_BLOCK *empty_sq_block;
  int           nworkers = esl_threads_GetWorkerCount(obj);

  esl_workqueue_Reset(queue);
#if DEBUGSERIAL
  printf("master threads reset\n");
#endif

  esl_threads_WaitForStart(obj);

#if DEBUGSERIAL
  printf("master threads started\n");
#endif

  status = esl_workqueue_ReaderUpdate(queue, NULL, &new_sq_block);
  if (status != eslOK) cm_Fail("Work queue reader failed");

#if DEBUGSERIAL
  printf("master initial update (blocks: %d)\n", nblocks);
#endif 

  /* main loop: */
  for(i = 0; i < nblocks; i++) { 
    if(sq_blockA[i] != NULL && sq_blockA[i]->count > 0) { 
      sq_block = (ESL_SQ_BLOCK *) new_sq_block;
      sq_block = sq_blockA[i];
      status = esl_workqueue_ReaderUpdate(queue, sq_block, &new_sq_block);
      if (status != eslOK) cm_Fail("Work queue reader failed");

#if DEBUGSERIAL
      printf("master internal update\n");
#endif
    }
  }

  /* now send a empty sq_block to all workers signaling them to stop */
  empty_sq_block = esl_sq_CreateBlock(1);
  for(k = 0; k < nworkers; k++) { 
    status = esl_workqueue_ReaderUpdate(queue, empty_sq_block, &new_sq_block);
    if (status != eslOK) cm_Fail("Work queue reader failed");
#if DEBUGSERIAL
    printf("master termination update\n");
#endif
  }

  status = esl_workqueue_ReaderUpdate(queue, sq_block, NULL);
#if DEBUGSERIAL
  printf("master final update\n");
#endif

  /* wait for all the threads to complete */
  esl_threads_WaitForFinish(obj);
#if DEBUGSERIAL
  printf("master got finish\n");
#endif

  esl_workqueue_Complete(queue);  
#if DEBUGSERIAL
  printf("master completed\n");
#endif

  esl_sq_DestroyBlock(empty_sq_block);
  return status;
}

/* pipeline_thread()
 * 
 * Receive a block of sequences from the master, 
 * align them and store their parsetrees.
 */

static void 
pipeline_thread(void *arg)
{
  int           status;
  int           j;             /* counter over parsetrees */
  CM_ALNDATA **bdataA = NULL; /* CM_ALNDATA array for the current block */
  int           workeridx;
  WORKER_INFO  *info;
  ESL_THREADS  *obj;
  ESL_SQ_BLOCK *sq_block = NULL;
  void         *new_sq_block;
  char          errbuf[eslERRBUFSIZE];
  
#ifdef HAVE_FLUSH_ZERO_MODE
  /* In order to avoid the performance penalty dealing with sub-normal
   * values in the floating point calculations, set the processor flag
   * so sub-normals are "flushed" immediately to zero.
   * On OS X, need to reset this flag for each thread
   * (see TW notes 05/08/10 for details)
   */
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

#if DEBUGSERIAL
  printf("started thread %d\n", workeridx);
#endif

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

#if DEBUGSERIAL
  printf("got data %d\n", workeridx);
#endif

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &new_sq_block);
  if (status != eslOK) cm_Fail("Work queue worker failed\n");

#if DEBUGSERIAL
  printf("initial update %d\n", workeridx);
#endif

  /* loop until all sequences have been processed */
  sq_block = (ESL_SQ_BLOCK *) new_sq_block;
  while (sq_block->count > 0) { 
    if((status = ProcessAlignmentWorkunit(info->cm, errbuf, sq_block, info->mxsize, info->w, info->w_tot, &bdataA)) != eslOK) cm_Fail(errbuf);
    /* point info->dataA at newly collected data, including parsetrees */
    ESL_REALLOC(info->dataA, sizeof(CM_ALNDATA *) * (info->n + sq_block->count));
    for(j = 0; j < sq_block->count; j++) info->dataA[info->n + j] = bdataA[j];
    info->n += sq_block->count;
    free(bdataA);  /* don't free actual CM_ALNDATA objects, info->trA is pointing at them */
    bdataA = NULL;

    status = esl_workqueue_WorkerUpdate(info->queue, sq_block, &new_sq_block);
    if (status != eslOK) cm_Fail("Work queue worker failed");
    sq_block = (ESL_SQ_BLOCK *) new_sq_block;

#if DEBUGSERIAL
    printf("internal update %d\n", workeridx);
#endif
  }

  status = esl_workqueue_WorkerUpdate(info->queue, sq_block, NULL);
  if (status != eslOK) cm_Fail("Work queue worker failed");
  
#if DEBUGSERIAL
  printf("final update %d\n", workeridx);
#endif

  esl_threads_Finished(obj, workeridx);
  return;

 ERROR: 
  cm_Fail("out of memory");
  return;  /* NEVERREACHED */
}
#endif   /* HMMER_THREADS */


/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * set flags and a few parameters. cm_Configure configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;

  ///if(! esl_opt_GetBoolean(go, "--nonull3")) cm->search_opts |= CM_SEARCH_NULL3;

  /* set up alignment options in cm->align_opts */
  if(  esl_opt_GetBoolean(go, "--nonbanded")) cm->align_opts |= CM_ALIGN_NONBANDED;
  else                                        cm->align_opts |= CM_ALIGN_HBANDED;
  if(  esl_opt_GetBoolean(go, "--optacc"))    cm->align_opts |= CM_ALIGN_OPTACC;
  if(  esl_opt_GetBoolean(go, "--cyk"))       cm->align_opts |= CM_ALIGN_CYK;
  if(! esl_opt_GetBoolean(go, "--noprob"))    cm->align_opts |= CM_ALIGN_POST;
  if(  esl_opt_GetBoolean(go, "--sample"))    cm->align_opts |= CM_ALIGN_SAMPLE;
  if(! esl_opt_GetBoolean(go, "--notrunc"))   cm->align_opts |= CM_ALIGN_TRUNC;
  if(  esl_opt_GetBoolean(go, "--sub"))       cm->align_opts |= CM_ALIGN_SUB;   /* --sub requires --notrunc and -g */
  if(  esl_opt_GetBoolean(go, "--small"))     cm->align_opts |= CM_ALIGN_SMALL; /* --small requires --noprob --nonbanded --cyk */

  /* set up configuration options in cm->config_opts */
  if(! esl_opt_GetBoolean(go, "--notrunc"))   cm->config_opts |= CM_CONFIG_TRUNC;
  if(  esl_opt_GetBoolean(go, "--sub"))       cm->config_opts |= CM_CONFIG_SUB;   /* --sub requires --notrunc and -g */
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }

  /* configure */
  if((status = cm_Configure(cm, errbuf, -1)) != eslOK) return status; 

  return eslOK;
}

static int
output_alignment(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, FILE *ofp, CM_ALNDATA **dataA, int ndata, char *map_sscons)
{
  int           status;
  ESL_MSA      *msa = NULL;
  int           j;
  ESL_SQ      **sqpA = NULL;   /*  array of sequence pointers,  only nec for Parsetrees2Alignment() */
  Parsetree_t **trA  = NULL;   /*  array of Parsetree pointers, only nec for Parsetrees2Alignment() */
  char        **ppstrA = NULL; /*  array of PP string pointers, only nec for Parsetrees2Alignment() */
  float         sc;
  float         struct_sc;
  int           first_ali = (dataA[0]->idx == 0) ? TRUE : FALSE;
  int           afmt;         /* format to output in */
  int           cpos, apos;   /* counters over consensus positions, alignment positions */

  /* contract check */
  if(esl_opt_GetBoolean(go, "-i") && ofp == cfg->tmpfp)        ESL_FAIL(eslEINVAL, errbuf, "-i enabled, but trying to output to temporary alignment file. This shouldn't happen.");
  if(esl_opt_GetBoolean(go, "--mapstr") && map_sscons == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mapstr enabled, but SS_cons not read from the --mapali alignment.");

  /* output the parsetrees, if nec */
  if(cfg->tfp != NULL) { 
    for (j = 0; j < ndata; j++) { 
      if((status = ParsetreeScore(cm, NULL, errbuf, dataA[j]->tr, dataA[j]->sqp->dsq, FALSE, &sc, &struct_sc, NULL, NULL, NULL)) != eslOK) return status;
      fprintf(cfg->tfp, ">%s\n", dataA[j]->sqp->name);
      fprintf(cfg->tfp, "  %16s %.2f bits\n", "SCORE:", sc);
      fprintf(cfg->tfp, "  %16s %.2f bits\n", "STRUCTURE SCORE:", struct_sc);
      ParsetreeDump(cfg->tfp, dataA[j]->tr, cm, dataA[j]->sqp->dsq, NULL, NULL); /* NULLs are dmin, dmax */
      fprintf(cfg->tfp, "//\n");
    }
  }

  /* print per-CM info to insertfp and elfp, if nec */
  if(first_ali == 0 && cfg->ifp != NULL) { fprintf(cfg->ifp, "%s %d\n", cm->name, cm->clen); } 
  if(first_ali == 0 && cfg->efp != NULL) { fprintf(cfg->efp, "%s %d\n", cm->name, cm->clen); } 

  /* create the alignment, Parsetrees2Alignment frees the parsetrees as it is creating the alignment */
  ESL_ALLOC(sqpA,   sizeof(ESL_SQ *)      * ndata); for(j = 0; j < ndata; j++) sqpA[j]   = dataA[j]->sqp;
  ESL_ALLOC(trA,    sizeof(Parsetree_t *) * ndata); for(j = 0; j < ndata; j++) trA[j]    = dataA[j]->tr;
  ESL_ALLOC(ppstrA, sizeof(char *)        * ndata); for(j = 0; j < ndata; j++) ppstrA[j] = dataA[j]->ppstr;
  if((status = Parsetrees2Alignment(cm, errbuf, cfg->abc_out, sqpA, NULL, trA, ppstrA, ndata, cfg->ifp, cfg->efp, TRUE, esl_opt_GetBoolean(go, "--matchonly"), &msa)) != eslOK) return status;

  if(ofp == cfg->rfp) { /* --regress file, remove GF author annotation */
    free(msa->au);
    msa->au = NULL;
  }

  /* rewrite SS_cons if --mapstr used */
  if(esl_opt_GetBoolean(go, "--mapstr")) { 
    /* step along the existing SS_cons, overwriting consensus positions in place */
    cpos = 0; /* span 0..clen-1 */
    for(apos = 0; apos < msa->alen; apos++) { /* span 0..alen-1 */
      if((! esl_abc_CIsGap    (cm->abc, msa->rf[apos])) && 
	 (! esl_abc_CIsMissing(cm->abc, msa->rf[apos]))) { 
	msa->ss_cons[apos] = map_sscons[cpos++];    
      }
    }
  }

  /* Determine format: we print in interleaved Stockholm if we're not
   * outputting to the tmpfile and --oneline was not used. This will
   * happen if -i was used or if we predicted we could fit all target
   * sequence parsetrees in a reasonable amount of memory. If this
   * wasn't the case, we're outputting to a tmpfile because we're
   * concerned we may run out of memory, and we'll go back and merge
   * all the alignments in the tmpfile once we're finished with 
   * all target sequences.
   */
  afmt   = (ofp == cfg->tmpfp || esl_opt_GetBoolean(go, "--oneline")) ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM;
  status = eslx_msafile_Write(ofp, msa, afmt);
  /* the contract asserted that if -i then ofp != cfg->tmpfp */
  if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
  else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);

  if(msa    != NULL) esl_msa_Destroy(msa);
  if(sqpA   != NULL) free(sqpA);
  if(trA    != NULL) free(trA);
  if(ppstrA != NULL) free(ppstrA);

  return eslOK;

 ERROR:
  if(msa    != NULL) esl_msa_Destroy(msa);
  if(sqpA   != NULL) free(sqpA);
  if(trA    != NULL) free(trA);
  if(ppstrA != NULL) free(ppstrA);
  return status;
}


/* Function: output_info_file_header
 * Date:     EPN, Fri Dec  4 08:15:31 2009
 *
 * Purpose:  Print the header section of an insert or EL insert
 *           (--ifile, --elfile) information file.
 *
 * Returns:  void
 */
void
output_info_file_header(FILE *fp, char *firstline, char *elstring)
{
  fprintf(fp, "# %s\n", firstline);
  fprintf(fp, "# This file includes 2+<nseq> non-'#' pre-fixed lines per model used for alignment,\n");
  fprintf(fp, "# where <nseq> is the number of sequences in the target file.\n");
  fprintf(fp, "# The first non-'#' prefixed line per model includes 2 tokens, separated by a single space (' '):\n");
  fprintf(fp, "# The first token is the model name and the second is the consensus length of the model (<clen>).\n");
  fprintf(fp, "# The following <nseq> lines include (4+3*<n>) whitespace delimited tokens per line.\n");
  fprintf(fp, "# The format for these <nseq> lines is:\n");
  fprintf(fp, "#   <seqname> <seqlen> <spos> <epos> <c_1> <u_1> <i_1> <c_2> <u_2> <i_2> .... <c_x> <u_x> <i_x> .... <c_n> <u_n> <i_n>\n");
  fprintf(fp, "#   indicating <seqname> has >= 1 %sinserted residues after <n> different consensus positions,\n", elstring);
  fprintf(fp, "#   <seqname> is the name of the sequence\n");
  fprintf(fp, "#   <seqlen>  is the unaligned length of the sequence\n");
  fprintf(fp, "#   <spos>    is the first (5'-most) consensus position filled by a nongap for this sequence (-1 if 0 nongap consensus posns)\n");
  fprintf(fp, "#   <epos>    is the final (3'-most) consensus position filled by a nongap for this sequence (-1 if 0 nongap consensus posns)\n");
  fprintf(fp, "#   <c_x> is a consensus position (between 0 and <clen>; if 0: inserts before 1st consensus posn)\n");
  fprintf(fp, "#   <u_x> is the *unaligned* position (b/t 1 and <seqlen>) in <seqname> of the first %sinserted residue after <c_x>.\n", elstring);
  fprintf(fp, "#   <i_x> is the number of %sinserted residues after position <c_x> for <seqname>.\n", elstring);
  fprintf(fp, "# Lines for sequences with 0 %sinserted residues will include only <seqname> <seqlen> <spos> <epos>.\n", elstring);
  fprintf(fp, "# The final non-'#' prefixed line per model includes only '//', indicating the end of info for a model.\n");
  fprintf(fp, "#\n");

  return;
}

/* Function: output_scores()
 * Date:     EPN, Tue Jan  3 14:49:56 2012
 *
 * Purpose:  Print scores and other information to a scores file.
 *
 * Returns:  eslOK on success.
 *           eslEMEM if out of memory.
 */
int
output_scores(FILE *ofp, CM_t *cm, char *errbuf, CM_ALNDATA **dataA, int ndata, int first_idx)
{
  int   status;               /* easel status */
  int   i;                    /* counter */
  int   namewidth = 8;        /* length of 'seq name' */
  char *namedashes = NULL;    /* namewidth-long string of dashes */
  
  /* alignment options */
  int do_nonbanded = (cm->align_opts & CM_ALIGN_NONBANDED) ? TRUE : FALSE;
  int do_post      = (cm->align_opts & CM_ALIGN_POST)      ? TRUE : FALSE;

  for(i = first_idx; i < ndata; i++) namewidth = ESL_MAX(namewidth, strlen(dataA[i]->sqp->name));

  ESL_ALLOC(namedashes, sizeof(char) * (namewidth+1));
  namedashes[namewidth] = '\0';
  for(i = 0; i < namewidth; i++) namedashes[i] = '-';

  fprintf(ofp, "# %9s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %-30s\n",    "",          namewidth,         "",      " ",        "",        "",      "",         "",       "", "         running time");
  fprintf(ofp, "# %9s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %30s\n",     "",          namewidth,         "",      " ",        "",        "",      "",         "",       "", "-------------------------------");
  fprintf(ofp, "# %9s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %9s  %9s  %9s\n", "seq idx",   namewidth, "seq name", "length", "cm from",   "cm to", "trunc",   "bit sc", "avg pp", "band calc", "alignment", "total");
  fprintf(ofp, "# %9s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %9s  %9s  %9s\n", "---------", namewidth, namedashes, "------", "-------", "-------", "-----", "--------", "------", "---------", "---------", "---------");

  for(i = first_idx; i < ndata; i++) { 
    fprintf(ofp, "  %9" PRId64 "  %-*s  %6" PRId64 "  %7d  %7d", dataA[i]->idx+1, namewidth, dataA[i]->sqp->name, dataA[i]->sqp->n, dataA[i]->cm_from, dataA[i]->cm_to);
    if     (dataA[i]->tr->mode[0] == TRMODE_T) fprintf(ofp, "  %5s", "5'&3'");
    else if(dataA[i]->tr->mode[0] == TRMODE_L) fprintf(ofp, "  %5s", "3'");
    else if(dataA[i]->tr->mode[0] == TRMODE_R) fprintf(ofp, "  %5s", "5'");
    else                                       fprintf(ofp, "  %5s", "no");
    fprintf(ofp, "  %8.2f", dataA[i]->sc);
    if(do_post)        fprintf(ofp, "  %6.3f", dataA[i]->pp);
    else               fprintf(ofp, "  %6s",   "-");
    if(! do_nonbanded) fprintf(ofp, "  %9.2f", dataA[i]->secs_bands);
    else               fprintf(ofp, "  %9s",   "-");
    fprintf(ofp, "  %9.2f  %9.2f\n", dataA[i]->secs_aln, dataA[i]->secs_tot);
  }

  if(namedashes != NULL) free(namedashes);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVER REACHED */
}

/* Function: create_and_output_final_msa
 * Incept:   EPN, Mon Dec 14 05:35:51 2009
 *
 * Purpose:  Read the >=1 MSAs that were written to a temporary file,
 *           merge them and output the merged MSA to a file without
 *           storing any of the full MSAs (incl. the final one) in
 *           memory.  To accomplish this a first pass of reading is
 *           done to determine how many gap columns must be added to
 *           each MSA to create the merged MSA during which only non
 *           per-sequence information is stored. After this pass, with
 *           the size of the merged alignment known, a second pass
 *           occurs during which only GS annotation is regurgitated
 *           (if any exists in at least 1 aln). Then a final pass
 *           occurs during which all other per-sequence data (PPs,
 *           aligned seqs) are regurgitated, taking care to add gap
 *           columns as necessary to make each input alignment the
 *           correct width of the merged alignment.
 *
 * Args:     go      - options
 *           cfg     - cmalign config
 *           errbuf  - for error messages
 *           cm      - CM used for alignment, useful for cm->clen
 *           tmpfile - name of temporary file with alignments to merge
 * 
 * Returns:   <eslOK> on success. 
 *            Returns <eslEOF> if there are no more alignments in <afp>.
 *            <eslEFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem. <eslEMEM> on allocation
 *            error.
 *
 * Xref:      /groups/eddy/home/nawrockie/notebook/9_1211_inf_cmalign_memeff/
 */
int 
create_and_output_final_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nali, char *tmpfile) 
{
  int           status;
  int           ai;                            /* counters over alignments */
  int           nseq_tot;                      /* number of sequences in all alignments */
  int           nseq_cur;                      /* number of sequences in current alignment */
  int64_t       alen_cur;                      /* length of current alignment */
  int64_t      *alenA = NULL;                  /* [0..nali_tot-1] alignment length of input msas (even after 
						* potentially removingeinserts (--rfonly)) */
  ESL_MSA     **msaA = NULL;                   /* [0..nali_tot-1] all msas read from all files */
  int          *maxins = NULL;                 /* [0..cpos..cm->clen+1] max number of inserts 
						* before each consensus position in all alignments */
  int          *maxel = NULL;                  /* [0..cpos..cm->clen+1] max number of EL inserts ('~' missing data symbols) 
						* before each consensus position in all alignments */
  int           cur_clen;                      /* consensus length (non-gap #=GC RF length) of current alignment */
  int           apos;                          /* alignment position */
  ESL_MSA      *fmsa = NULL;                   /* the merged alignment created by merging all alignments in msaA */
  int           alen_fmsa;                     /* number of columns in merged MSA */
  int          *ngap_insA = NULL;               /* [0..alen] number of insert gap columns to add after each alignment column when merging */
  int          *ngap_elA = NULL;                /* [0..alen] number of missing data ('~') gap columns to add after each alignment column when merging */
  int          *ngap_eitherA = NULL;            /* [0..apos..alen] = ngap_insA[apos] + ngap_elA[apos] */
  char         *rf2print = NULL;                /* #=GC RF annotation for final alignment */
  char         *ss_cons2print = NULL;           /* #=GC SS_cons annotation for final alignment */

  /* variables only used in small mode */
  int           ngs_cur;                       /* number of GS lines in current alignment (only used if do_small) */
  int           gs_exists = FALSE;             /* set to TRUE if do_small and any input aln has >= 1 GS line */
  int           maxname, maxgf, maxgc, maxgr;  /* max length of seqname, GF tag, GC tag, GR tag in all input alignments */
  int           maxname_cur, maxgf_cur, maxgc_cur, maxgr_cur; /* max length of seqname, GF tag, GC tag, GR tag in current input alignment */
  int           margin = 0;                    /* total margin length for output msa */
  int           regurg_header = FALSE;         /* set to TRUE if we're printing out header */
  int           regurg_gf     = FALSE;         /* set to TRUE if we're printing out GF */
  ESL_MSAFILE2 *afp;

  /* Allocate and initialize */
  ESL_ALLOC(msaA,   sizeof(ESL_MSA *) * nali);
  ESL_ALLOC(alenA,  sizeof(int64_t) * nali);

  /****************************************************************************
   * Read alignments one at a time, storing all non-sequence info, separately *
   ****************************************************************************/
  if((status = esl_msafile2_Open(tmpfile, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading", tmpfile);

  ai = 0;
  nseq_tot = 0;
  maxname = maxgf = maxgc = maxgr = 0;

  /* allocate maxins */
  ESL_ALLOC(maxins, sizeof(int) * (cm->clen+1)); 
  esl_vec_ISet(maxins, (cm->clen+1), 0);
  /* allocate maxel */
  ESL_ALLOC(maxel, sizeof(int) * (cm->clen+1)); 
  esl_vec_ISet(maxel, (cm->clen+1), 0); /* these will all stay 0 unless we see '~' in the alignments */

  /* read all alignments, there should be nali of them */
  for(ai = 0; ai < nali; ai++) { 
    status = esl_msafile2_ReadInfoPfam(afp, NULL, cfg->abc, -1, NULL, NULL, &(msaA[ai]), &nseq_cur, &alen_cur, &ngs_cur, &maxname_cur, &maxgf_cur, &maxgc_cur, &maxgr_cur, NULL, NULL, NULL, NULL, NULL);
    if      (status == eslEFORMAT) cm_Fail("Rereading alignment %d for merging, parse error:\n%s\n", ai+1, afp->errbuf);
    else if (status == eslEINVAL)  cm_Fail("Rereading alignment %d for merging, parse error:\n%s\n", ai+1, afp->errbuf);
    else if (status != eslOK)      cm_Fail("Rereading alignment %d for merging, parse error:\n%s\n", ai+1, afp->errbuf);

    msaA[ai]->abc = cfg->abc; 
    if(msaA[ai]->rf == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "When rereading alignment %d for merging, no RF annotation found.", ai+1);
    cur_clen = 0;
    for(apos = 0; apos < (int) alen_cur; apos++) { 
      if((! esl_abc_CIsGap(msaA[ai]->abc, msaA[ai]->rf[apos])) && (! esl_abc_CIsMissing(msaA[ai]->abc, msaA[ai]->rf[apos]))) cur_clen++;
    }
    if(cur_clen != cm->clen) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "When rereading alignment %d for merging, consensus length wrong (%d, when %d was expected)", ai, cur_clen, cm->clen);
    maxname = ESL_MAX(maxname, maxname_cur); 
    maxgf   = ESL_MAX(maxgf, maxgf_cur); 
    maxgc   = ESL_MAX(maxgc, maxgc_cur); 
    maxgr   = ESL_MAX(maxgr, maxgr_cur); 
    msaA[ai]->alen = alen_cur;
    alenA[ai]      = alen_cur; /* to remember total width of aln to expect in second pass */
    nseq_tot += nseq_cur;
    if(ngs_cur > 0) gs_exists = TRUE; 
      
    /* determine max number inserts and ELs between each position */
    update_maxins_and_maxel(msaA[ai], cm->clen, msaA[ai]->alen, maxins, maxel);
  }
  /* final check, make sure we've read all msas from the file, we should have, we only printed nali */
  status = esl_msafile2_ReadInfoPfam(afp, NULL, cfg->abc, -1, NULL, NULL, NULL, 
				     NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
				     NULL, NULL, NULL, NULL, NULL);
  if(status != eslEOF) ESL_FAIL(status, errbuf, "More alignments in temp file than expected.");
  esl_msafile2_Close(afp);
  
  /*********************************************
   *  Merge all alignments into the merged MSA *
   *********************************************/

  /* We allocate space for all sequences, but leave sequences as NULL
   * (nseq = -1).  We didn't store the sequences on the first pass
   * through the alignment files, and we'll never allocate space for
   * the sequences in fmsa, we'll just output them as we reread them
   * on another pass through the individual alignments. If we read >=
   * 1 GS line in any of the temporary alignments, we need to do an
   * additional pass through them, outputting only GS data. Then, in a
   * final (3rd) pass we'll output aligned data.
   */     
  fmsa = esl_msa_Create(nseq_tot, -1); 
  alen_fmsa = cm->clen + esl_vec_ISum(maxins, (cm->clen+1)); 

  /* if there was any GS annotation in any of the individual alignments,
   * do second pass through alignment files, outputting GS annotation as we go. */
  if(gs_exists) { 
    if((status = esl_msafile2_Open(tmpfile, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading on second pass", tmpfile);
    for(ai = 0; ai < nali; ai++) { 
      regurg_header = (ai == 0) ? TRUE : FALSE;
      regurg_gf     = (ai == 0) ? TRUE : FALSE;
      status = esl_msafile2_RegurgitatePfam(afp, cfg->ofp, 
					    maxname, maxgf, maxgc, maxgr, /* max width of a seq name, gf tag, gc tag, gr tag */
					    regurg_header, /* regurgitate stockholm header ? */
					    FALSE,         /* regurgitate // trailer ? */
					    regurg_header, /* regurgitate blank lines */
					    regurg_header, /* regurgitate comments */
					    regurg_gf,     /* regurgitate GF ? */
					    TRUE,          /* regurgitate GS ? */
					    FALSE,         /* regurgitate GC ? */
					    FALSE,         /* regurgitate GR ? */
					    FALSE,         /* regurgitate aseq ? */
					    NULL,          /* output all seqs, not just those stored in a keyhash */
					    NULL,          /* output all seqs, don't skip those listed in a keyhash */
					    NULL,          /* useme,  irrelevant, we're only outputting GS */
					    NULL,          /* add2me, irrelevant, we're only outputting GS */
					    alenA[ai], /* alignment length, as we read it in first pass (inserts may have been removed since then) */
					    '.', NULL, NULL);
      if(status == eslEOF) cm_Fail("Second pass, error out of temp alignments too soon, when trying to read aln %d", ai);
      if(status != eslOK)  cm_Fail("Second pass, error reading temp alignment %d %s", ai, afp->errbuf); 
      fflush(cfg->ofp);
    }
    esl_msafile2_Close(afp);
    fprintf(cfg->ofp, "\n"); /* a single blank line to separate GS annotation from aligned data */
  }
  /* do another (either second or third) pass through alignment files, outputting aligned sequence data (and GR) as we go */

  if((status = esl_msafile2_Open(tmpfile, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading on second (or third) pass", tmpfile);

  for(ai = 0; ai < nali; ai++) { 
    /* determine how many all gap columns to insert after each alignment position
     * of the temporary msa when copying it to the merged msa */
    if((status = determine_gap_columns_to_add(msaA[ai], maxins, maxel, cm->clen, &(ngap_insA), &(ngap_elA), &(ngap_eitherA), errbuf)) != eslOK) 
      cm_Fail("error determining number of all gap columns to add to temp alignment %d", ai);
    regurg_header = ((! gs_exists) && (ai == 0)) ? TRUE : FALSE;
    regurg_gf     = ((! gs_exists) && (ai == 0)) ? TRUE : FALSE;

    status = esl_msafile2_RegurgitatePfam(afp, cfg->ofp,
					  maxname, maxgf, maxgc, maxgr, /* max width of a seq name, gf tag, gc tag, gr tag */
					  regurg_header,  /* regurgitate stockholm header ? */
					  FALSE,          /* regurgitate // trailer ? */
					  regurg_header,  /* regurgitate blank lines */
					  regurg_header,  /* regurgitate comments */
					  regurg_gf,      /* regurgitate GF ? */
					  FALSE,          /* regurgitate GS ? */
					  FALSE,          /* regurgitate GC ? */
					  TRUE,           /* regurgitate GR ? */
					  TRUE,           /* regurgitate aseq ? */
					  NULL,           /* output all seqs, not just those stored in a keyhash */
					  NULL,           /* output all seqs, don't skip those stored in a keyhash */
					  NULL,           /* useme, not nec b/c we want to keep all columns */
					  ngap_eitherA,   /* number of all gap columns to add after each apos */
					  alenA[ai],      /* alignment length, as we read it in first pass, not strictly necessary */
					  '.', NULL, NULL);
    if(status == eslEOF) cm_Fail("Second pass, error out of alignments too soon, when trying to read temp aln %d", ai);
    if(status != eslOK)  cm_Fail("Second pass, error reading temp alignment %d: %s", ai, afp->errbuf); 
    if(ai == 0) { 
      /* create the GC SS_cons and GC RF to print from the first alignment,
       * we use the first alignment b/c this is the one potentially with rewritten pknots
       * from --withpknots.
       */
      inflate_gc_with_gaps_and_els(cfg->ofp, msaA[ai], ngap_insA, ngap_elA, &ss_cons2print, &rf2print);
    }
    free(ngap_insA);
    free(ngap_elA);
    free(ngap_eitherA);
    
    esl_msa_Destroy(msaA[ai]);
    msaA[ai] = NULL;
    fflush(cfg->ofp);
  }
  /* output SS_cons and RF */
  margin = maxname+1;
  if (maxgc > 0 && maxgc+6         > margin) margin = maxgc+6;
  if (maxgr > 0 && maxname+maxgr+7 > margin) margin = maxname+maxgr+7; 
  fprintf(cfg->ofp, "#=GC %-*s %s\n", margin-6, "SS_cons", ss_cons2print);
  fprintf(cfg->ofp, "#=GC %-*s %s\n", margin-6, "RF", rf2print);
  fprintf(cfg->ofp, "//\n");

  esl_msafile2_Close(afp);

  if(ss_cons2print != NULL) free(ss_cons2print);
  if(rf2print != NULL) free(rf2print);
  if(alenA != NULL)  free(alenA);
  if(msaA != NULL)   free(msaA);
  if(maxins != NULL) free(maxins);
  if(maxel != NULL)  free(maxel);
  if(fmsa != NULL)   esl_msa_Destroy(fmsa);
  return eslOK;

 ERROR: 
  esl_fatal("Out of memory. Reformat to Pfam with esl-reformat and try esl-alimerge --savemem.");
  return eslEMEM; /*NEVERREACHED*/
}

/* Function: update_maxins_and_maxel
 * Date:     EPN, Sun Nov 22 09:40:48 2009
 * 
 * Update maxins[] and maxel[], arrays that keeps track of the
 * max number of inserted ('.' gap #=GC RF) columns and inserted EL
 * emissions ('~' gap #=GC RF) columns before each cpos (consensus
 * (non-gap #=GC RF) column)
 *
 * Consensus columns are index [0..cpos..clen].
 * 
 * max{ins,el}[0]      is number of {IL/IR inserts, EL inserts} before 1st cpos.
 * max{ins,el}[clen-1] is number of {IL/IR inserts, EL inserts} before final cpos.
 * max{ins,el}[clen]   is number of {IL/IR inserts, EL inserts} after  final cpos.
 * 
 * Caller has already checked that msa->rf != NULL
 * and its non-gap length is clen. If we find either
 * of these is not true, we die (but this shouldn't happen).
 * 
 * Returns: void.
 */
void
update_maxins_and_maxel(ESL_MSA *msa, int clen, int64_t alen, int *maxins, int *maxel) 
{
  int apos;
  int cpos = 0;
  int nins = 0;
  int nel = 0;

  for(apos = 0; apos < alen; apos++) { 
    if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      nins++;
    }
    else if (esl_abc_CIsMissing(msa->abc, msa->rf[apos])) { 
      nel++;
    }
    else {
      maxins[cpos] = ESL_MAX(maxins[cpos], nins);
      maxel[cpos]  = ESL_MAX(maxel[cpos], nel);
      cpos++;
      nins = 0;
      nel = 0;
    }
  }
      
  /* update final value, max{ins,el}[clen+1], the number of inserts
   * after the final consensus position */
  maxins[cpos] = ESL_MAX(maxins[cpos], nins);
  maxel[cpos]  = ESL_MAX(maxel[cpos], nel);
  if(cpos != clen) cm_Fail("Unexpected error in update_maxins_and_maxel(), expected clen (%d) not equal to actual clen (%d).\n", clen, cpos);

  return;
}

/* determine_gap_columns_to_add
 *                   
 * Given <maxins> and <maxel>, two arrays of the number of gap RF
 * (inserts) positions and '~' RF (EL inserts) after each non-gap RF 
 * (consensus) position in the eventual final merged alignment, 
 * calculate how many inserts and missing data inserts
 * we need to add at each position of <msa> to expand it out to the 
 * appropriate size of the eventual merged alignment.
 * 
 * max{ins,el}[0]      is number of inserts,ELs before 1st cpos in merged aln
 * max{ins,el}[cpos]   is number of inserts,ELs after  final cpos in merged aln
 *                             for cpos = 1..clen 
 * clen is the number of non-gap RF positions in msa (and in eventual merged msa).             
 * 
 * We allocate fill and return ret_ngap_insA[0..msa->alen], ret_ngap_elA[0..msa->alen], 
 * and ret_ngap_eitherA[0..msa->alen] here.
 *
 * ret_n{ins,el}gapA[0]      is number of inserts,ELs to add before 1st position of msa 
 * ret_n{ins,el}gapA[apos]   is number of inserts,ELs to add after alignment position apos
 *                             for apos = 1..msa->alen
 * 
 * ret_ngap_eitherA[apos] = ngap_insA[apos] + ngap_elA[apos]
 * 
 * This is similar to the esl_msa.c helper function of the same name,
 * but that function does not bother with missing data '~'.
 * 
 * Returns eslOK on success.
 *         eslEMEM on memory alloaction error 
 *         eslERANGE if a value exceeds what we expected (based on earlier
 *                   checks before this function was entered).
 *         if !eslOK, errbuf if filled.
 */
int
determine_gap_columns_to_add(ESL_MSA *msa, int *maxins, int *maxel, int clen, int **ret_ngap_insA, int **ret_ngap_elA, int **ret_ngap_eitherA, char *errbuf)
{
  int status;
  int apos;
  int prv_cpos = 0;  /* alignment position corresponding to consensus position cpos-1 */
  int cpos = 0;
  int nins = 0;
  int nel = 0;
  int *ngap_insA = NULL;
  int *ngap_elA = NULL;
  int *ngap_eitherA = NULL;
  int insert_before_el_flag = FALSE; /* set to TRUE for a cpos if we observe an insert column 5' of a missing data column between cpos-1 and cpos */
  int el_before_insert_flag = FALSE; /* set to TRUE for a cpos if we observe a missing data column 5' of an insert column between cpos-1 and cpos */

  /* contract check */
  if(maxel[0] != 0) ESL_FAIL(eslEINVAL, errbuf, "missing characters exist prior to first cpos, this shouldn't happen.\n");
  if(msa->ss_cons == NULL) ESL_FAIL(eslEINVAL, errbuf, "MSA's SS_cons is null in determine_gap_columns_to_add.\n");

  ESL_ALLOC(ngap_insA, sizeof(int) * (msa->alen+1));
  ESL_ALLOC(ngap_elA, sizeof(int) * (msa->alen+1));
  ESL_ALLOC(ngap_eitherA, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(ngap_insA, (msa->alen+1), 0);
  esl_vec_ISet(ngap_elA, (msa->alen+1), 0);
  esl_vec_ISet(ngap_eitherA, (msa->alen+1), 0);
  
  for(apos = 0; apos < msa->alen; apos++) { 
    if(esl_abc_CIsMissing(msa->abc, msa->rf[apos])) { 
      nel++;
      if(nins > 0) insert_before_el_flag = TRUE;
    }
    else if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      nins++;
      if(nel > 0) el_before_insert_flag = TRUE;
    }
    else { /* a consensus position */
      /* a few sanity checks */
      if(nins  > maxins[cpos])  ESL_FAIL(eslEINCONCEIVABLE, errbuf, "%d inserts before cpos %d greater than max expected (%d).\n", nins, cpos, maxins[cpos]); 
      if(nel > maxel[cpos]) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "%d EL inserts before cpos %d greater than max expected (%d).\n", nel, cpos, maxel[cpos]); 
      if(insert_before_el_flag && el_before_insert_flag) ESL_XFAIL(eslEINVAL, errbuf, "found inserts then missing chars '~' then inserts in RF before cpos %d.\n", cpos);

      if (cpos == 0) { 
	if(nel != 0) ESL_FAIL(eslEINVAL, errbuf, "found missing chars prior to first cpos, shouldn't happen\n");
	ngap_insA[prv_cpos]  = maxins[cpos] - nins; /* inserts before first position: flush right (so add all-gap columns after leftmost column) */
	/* we already checked that maxel[0] is 0 during contract check above */
      }
      else {
	/* Determine where to place inserts and/or missing data.
	 * Handle each of 4 possibilities separately, note that if 
	 * maxins[cpos] == 0 then nins == 0, and if maxel[cpos] == 0 then nel == 0 (see sanity check above). 
	 */
	if(maxins[cpos] >  0 && maxel[cpos] == 0) { /* most common case */
	  ngap_insA[prv_cpos + 1 + (nins/2)] = maxins[cpos] - nins; /* internal cpos: split */
	}
	else if(maxins[cpos] == 0 && maxel[cpos] > 0) { 
	  ngap_elA[prv_cpos + 1 + (nel/2)] = maxel[cpos] - nel; /* internal cpos: split */
	}
	else if(maxins[cpos] >  0 && maxel[cpos] > 0) { 
	  /* We have to determine 5'->3' order: it could be inserts 5'
	   * (before) missing data or missing data 5' (before)
	   * inserts.  Rule for this is convoluted but we can
	   * determine which order based on the consensus structure:
	   * ff cpos-1 and cpos form a base pair (SS_cons[cpos-1] ==
	   * '<' && SS_cons[cpos] == '>'), then we do inserts 5' of
	   * missing data, else we do missing data 5' of inserts.  
	   * 
	   * The reason is b/c all ELs (missing data) occur at the end
	   * of a stem (prior to END_nd's left consensus position).
	   * Inserts emitted by states immediately prior to an END_nd
	   * are involved in the lone case of ambiguity in the CM
	   * grammar because it indicates an insert position to which
	   * two insert states (one IL and one IR) emit. This
	   * ambiguity is dealt with by 'detaching' (see
	   * cm_modelmaker.c:cm_find_and_detach_dual_inserts()) insert
	   * states with index equal to the state index of an END_E
	   * state *minus 1*. These detached states are always IL
	   * states that would emit inserts 5' of any ELs (missing
	   * data) (e.g. a MATL_IL just before an END_E) *except* when
	   * a MATP_nd is immediately prior to an END_E node, in which
	   * case the detached state is an IR state (always MATP_IR
	   * just before END_E). Because these states are detached we
	   * know that the *other* state that inserts at that position
	   * must be used.  Either an IR (ROOT_IR, MATR_IR or MATP_IR)
	   * in the former case or a MATP_IL in the latter. The IR
	   * will insert 3' of an EL, and the MATP_IL will insert 5'
	   * of an EL. We can determine the MATP_IL cases by asking if
	   * cpos-1 and cpos form a basepair in the SS_cons, b/c
	   * that's the only case in which a MATP node exists
	   * immediately before an END_nd. Note that this basepair
	   * will always be a '<' '>' basepair (as opposed to a '('
	   * ')' or '{' '}' for example) b/c it is the lowest nesting
	   * order (enclosing basepair of stem. (Told you it was
	   * convoluted.)
	   */
	  if((msa->ss_cons[prv_cpos] == '<') && (msa->ss_cons[apos] == '>')) { /* special case */
	    /* inserts should be 5' of missing data (ELs) */
	    ngap_insA[prv_cpos + 1        + (nins/2)] = maxins[cpos] - nins; /* internal cpos: split */
	    ngap_elA[prv_cpos  + 1 + nins + (nel/2)] = maxel[cpos] - nel; /* internal cpos: split */
	  }
	  else { 
	    /* missing data (ELs) should be 5' of inserts */
	    ngap_elA[prv_cpos + 1 + (nel/2)] = maxel[cpos] - nel; /* internal cpos: split */
	    ngap_insA[prv_cpos + 1 + nel + (nins/2)] = maxins[cpos] - nins; /* internal cpos: split */
	  }
	}
	/* final case is if (maxins[cpos] == 0 && maxel[cpos] == 0) 
	 * in this case we do nothing. 
	 */
      }

      cpos++;
      prv_cpos = apos;
      nins = 0;
      nel = 0;
      el_before_insert_flag = FALSE;
      insert_before_el_flag = FALSE;
    }
  }
  cpos = clen;
  apos = msa->alen;
  /* deal with inserts after final consensus position, we could have both inserts and missing data, but only in case where 
   * CM has 0 bps (model is ROOT, MATL, MATL, ..., MATL, END) so we know that MATL_IL has been detached, ROOT_IR emits inserts
   * after final cpos, and thus missing data (ELs) should be 5' of inserts 
   */
  if(maxins[cpos] > 0 && maxel[cpos] == 0) { /* most common case */
    ngap_insA[apos] = maxins[cpos] - nins; /* flush left inserts (no ELs) */
  }
  else if(maxins[cpos] == 0 && maxel[cpos] > 0) { 
    ngap_elA[apos] = maxel[cpos] - nel; /* flush left ELs (no inserts) */
  }
  else if(maxins[cpos] > 0 && maxel[cpos] > 0) { 
    if((msa->ss_cons[prv_cpos] == '<') && (msa->ss_cons[apos] == '>')) ESL_FAIL(eslEINVAL, errbuf, "final position needs ELs and inserts after it, and is right half of a base pair, this shouldn't happen");
    /* missing data (ELs) should be 5' of inserts */
    ngap_elA[apos]      = maxel[cpos] - nel; /* flush left */
    ngap_insA[apos+nel] = maxins[cpos] - nins; /* flush left, after ELs */
  }

  /* determine ngap_eitherA[], the number of gaps due to either inserts or missing data after each apos */
  for(apos = 0; apos <= msa->alen; apos++) { 
    ngap_eitherA[apos] = ngap_insA[apos] + ngap_elA[apos];
  }

  /* validate that clen is what it should be */
  if(cpos != clen) { 
    if(ngap_insA != NULL) free(ngap_insA);
    if(ngap_elA != NULL) free(ngap_elA);
    if(ngap_eitherA != NULL) free(ngap_eitherA);
    ESL_FAIL(eslEINCONCEIVABLE, errbuf, "consensus length (%d) is not the expected length (%d).", cpos, clen);
  }

  *ret_ngap_insA  = ngap_insA;
  *ret_ngap_elA = ngap_elA;
  *ret_ngap_eitherA = ngap_eitherA;

  return eslOK;

 ERROR: 
  if(ngap_insA  != NULL) free(ngap_insA);
  if(ngap_elA != NULL) free(ngap_elA);
  if(ngap_eitherA != NULL) free(ngap_eitherA);
  ESL_FAIL(status, errbuf, "Memory allocation error.");
  return status; /*NEVERREACHED*/
}

/* inflate_gc_with_gaps_and_els
 *                   
 * Given an MSA and two arrays specifying the number of inserts '.'
 * and EL inserts '~' to add after each position, create the 
 * SS_cons and RF strings to output for the merged alignment
 * after adding the gaps and missing data symbols ('~') and
 * return them in ret_ss_cons2print and ret_rf2print. Caller
 * is responsible for freeing them.
 *
 * Returns void. If something unexpected occurs, including an 
 * allocation error, we die here and print error message.
 */
void
inflate_gc_with_gaps_and_els(FILE *ofp, ESL_MSA *msa, int *ngap_insA, int *ngap_elA, char **ret_ss_cons2print, char **ret_rf2print) 
{
  int status;
  int apos  = 0;
  int apos2print  = 0;
  int i, j;
  int el_before_ins = FALSE;
  int prv_cpos = -1;
  int alen2print = 0;
  char *rf2print;
  char *ss_cons2print;

  alen2print = msa->alen + esl_vec_ISum(ngap_insA, msa->alen+1) + esl_vec_ISum(ngap_elA, msa->alen+1);
  ESL_ALLOC(rf2print,      sizeof(char) * (alen2print+1));
  ESL_ALLOC(ss_cons2print, sizeof(char) * (alen2print+1));
  rf2print[alen2print] = '\0';
  ss_cons2print[alen2print] = '\0';

  if(msa->ss_cons == NULL) cm_Fail("Error trying to add inserts to SS_cons, but unexpectedly it doesn't exist.");
  if(msa->rf      == NULL) cm_Fail("Error trying to add inserts to SS_cons, but unexpectedly it doesn't exist.");
  for(apos = 0; apos <= msa->alen; apos++) { 
    el_before_ins = TRUE; /* until proven otherwise */
    if(ngap_insA[apos] > 0 && ngap_elA[apos] > 0) { 
      /* Rare case: we need to add at least one '.' and '~'. We need
       * to determine which to add first. The only possible way we do
       * '.'s and then '~'s is if we're inserting between the two
       * match positions of a MATP node and the following node is an
       * END (search for 'convoluted' in comments in
       * determine_gap_columns_to_add() for explanation of why).
       * This will only occur if previously seen consensus SS_cons 
       * char was a '<' and next consensus SS_cons char is a '>'.
       * First, check if prev consensus SS_cons char is '<',
       * if it is check if following one is a '>' by searching for
       * it. This is inefficient but should be rare so it's okay 
       * (we'll only do this at most once per END node in model).
       */
      if((prv_cpos == -1) || (msa->ss_cons[prv_cpos] != '<')) el_before_ins = TRUE; /* normal case */
      else { /* prv_cpos is '<', check if next one is a '>' */
	j = i+1; 
	while(j < msa->alen && msa->ss_cons[j] == '.') j++;
	if((j < msa->alen) && (msa->ss_cons[j] == '>')) el_before_ins = FALSE;
      }
    }
    if(el_before_ins) { 
      for(i = 0; i < ngap_elA[apos]; i++) { 
	rf2print[apos2print] = '~';
	ss_cons2print[apos2print++] = '~';
      }
      for(i = 0; i < ngap_insA[apos]; i++) { 
	rf2print[apos2print] = '.';
	ss_cons2print[apos2print++] = '.';
      }
    } 
    else { /* insert before els, rare, only occurs if el_before_ins set to FALSE above */
      for(i = 0; i < ngap_insA[apos]; i++) { 
	rf2print[apos2print] = '.';
	ss_cons2print[apos2print++] = '.';
      }
      for(i = 0; i < ngap_elA[apos]; i++) { 
	rf2print[apos2print] = '~';
	ss_cons2print[apos2print++] = '~';
      }
    }
    if(apos < msa->alen) { 
      rf2print[apos2print]        = msa->rf[apos];
      ss_cons2print[apos2print++] = msa->ss_cons[apos];
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) prv_cpos = apos;
    }	
  }    

  *ret_ss_cons2print = ss_cons2print;
  *ret_rf2print = rf2print;

  return;

 ERROR:
  cm_Fail("Allocation error when creating final alignment RF and SS_cons.");
  return; /* NEVERREACHED */
}

/* map_alignment()
 *                   
 * Called if the --mapali <f> option is used. Open and read a 
 * MSA from <f>, confirm it was the same alignment used to 
 * build the CM, and convert its aligned sequences to 
 * parsetrees. Return data (sequences and parsetrees) gets populated into <ret_dataA>.
 * 
 * Also, return the dealigned (cm-> clen length) SS_cons
 * from the alignment in <ret_ss>. This will be used to 
 * overwrite the output alignment's SS_cons if --mapstr 
 * was used. 
 */
static int
map_alignment(const char *msafile, CM_t *cm, char *errbuf, CM_ALNDATA ***ret_dataA, int *ret_ndata, char **ret_ss)
{
  int            status;
  ESLX_MSAFILE  *afp       = NULL;
  ESL_MSA       *msa       = NULL;
  ESL_ALPHABET  *abc       = (ESL_ALPHABET *) cm->abc; /* removing const'ness to make compiler happy. Safe. */
  uint32_t       chksum    = 0;
  int            i, x;              /* counters */
  int            apos, uapos, cpos; /* counter over aligned, unaligned, consensus positions */
  int           *a2u_map   = NULL;  /* map from aligned to unaligned positions */
  Parsetree_t   *mtr       = NULL;  /* the guide tree for mapali */
  CM_ALNDATA  **dataA     = NULL;  /* includes ptrs to sq and parsetrees */
  char          *aseq      = NULL;  /* aligned sequence, req'd by Transmogrify() */
  char          *ss        = NULL;  /* msa's SS_cons, if there is one, dealigned to length cm->clen */
  
  status = eslx_msafile_Open(&abc, msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);

  status = eslx_msafile_Read(afp, &msa);
  if (status != eslOK) eslx_msafile_ReadFailure(afp, status);

  if (! (cm->flags & CMH_CHKSUM))  cm_Fail("CM has no checksum. --mapali unreliable without it.");
  if (! (cm->flags & CMH_MAP))     cm_Fail("CM has no map. --mapali can't work without it.");
  esl_msa_Checksum(msa, &chksum);
  if (cm->checksum != chksum)      cm_Fail("--mapali MSA %s isn't same as the one CM came from (checksum mismatch)", msafile);

  /* allocate and initialize dataA */
  ESL_ALLOC(dataA, sizeof(CM_ALNDATA *) * msa->nseq);
  for(i = 0; i < msa->nseq; i++) dataA[i] = cm_alndata_Create();

  /* get SS_cons from the msa possibly for --mapstr, important to do it here, before it is potentially deknotted in HandModelMaker() */
  if(msa->ss_cons != NULL) { 
    ESL_ALLOC(ss, sizeof(char) * (cm->clen+1));
    ss[cm->clen] = '\0';
    for(cpos = 1; cpos <= cm->clen; cpos++) ss[cpos-1] = msa->ss_cons[cm->map[cpos]-1];
  }

  /* add RF annotation to the msa, so we can use it in HandModelMaker() */
  if(msa->rf != NULL) free(msa->rf);
  ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
  /* init to all inserts, then set match states based on cm->map */
  for (apos = 0; apos <  msa->alen; apos++) msa->rf[apos] = '.';
  for (cpos = 1; cpos <= cm->clen;  cpos++) msa->rf[cm->map[cpos]-1] = 'x'; /* note off by one */

  /* create a guide tree, which we'll need to convert aligned sequences to parsetrees */
  status = HandModelmaker(msa, errbuf, 
			  TRUE,  /* use_rf */
			  FALSE, /* use_wts, irrelevant */
			  0.5,   /* gapthresh, irrelevant */
			  NULL,  /* returned CM, irrelevant */
			  &mtr); /* guide tree */
  if(status != eslOK) return status;

  /* create a parsetree from each aligned sequence */
  ESL_ALLOC(aseq,    sizeof(char) * (msa->alen+1));
  ESL_ALLOC(a2u_map, sizeof(int)  * (msa->alen+1));
  a2u_map[0] = -1; /* invalid */
  for (i = 0; i < msa->nseq; i++) { 
    /* we need a text sequence in addition to the digitized sequence for Transmogrify() */
    esl_abc_Textize(msa->abc, msa->ax[i], msa->alen, aseq);
    dataA[i]->tr = Transmogrify(cm, mtr, msa->ax[i], aseq, msa->alen);
    /* dataA[i]->tr is in alignment coords, convert it to unaligned coords.
     * First we construct a map of aligned to unaligned coords, then
     * we use it to convert. 
     */
    uapos = 1;
    for(apos = 1; apos <= msa->alen; apos++) { 
      a2u_map[apos] = (esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) ? -1 : uapos++; 
    }
    for(x = 0; x < dataA[i]->tr->n; x++) { 
      if(dataA[i]->tr->emitl[x] != -1) dataA[i]->tr->emitl[x] = a2u_map[dataA[i]->tr->emitl[x]];
      if(dataA[i]->tr->emitr[x] != -1) dataA[i]->tr->emitr[x] = a2u_map[dataA[i]->tr->emitr[x]];
    }
  }   

  /* get sequences */
  for (i = 0; i < msa->nseq; i++) esl_sq_FetchFromMSA(msa, i, &(dataA[i]->sqp));

  *ret_dataA = dataA;
  *ret_ndata = msa->nseq;
  *ret_ss    = ss;

  eslx_msafile_Close(afp);
  esl_msa_Destroy(msa);
  FreeParsetree(mtr);
  free(a2u_map);
  free(aseq);

  return eslOK;

 ERROR:
  *ret_ndata = 0;
  *ret_dataA = NULL;
  if (dataA     != NULL) { 
    for(i = 0; i < msa->nseq; i++) cm_alndata_Destroy(dataA[i], TRUE); 
    dataA = NULL;
  }
  if (afp       != NULL) eslx_msafile_Close(afp);
  if (msa       != NULL) esl_msa_Destroy(msa);
  if (a2u_map   != NULL) free(a2u_map);
  if (aseq      != NULL) free(aseq);  
  ESL_FAIL(status, errbuf, "out of memory");
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

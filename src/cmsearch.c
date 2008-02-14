/* cmsearch.c
 * SRE, Fri May  3 13:58:18 2002
 * SVN $Id$
 * 
 * Search sequences with a CM.
 * 
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"              /* better general sequence analysis library */
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define STRATOPTS1  "--cyk,--inside,--viterbi,--forward"               /* incompatible with --cyk, --inside (besides themselves) */
#define STRATOPTS2  "--cyk,--inside,--viterbi,--forward,--fil-hmm"     /* incompatible with --viterbi,--forward  (besides themselves) */
#define ALPHOPTS    "--rna,--dna"                                      /* exclusive choice for output alphabet */

#define CUTOPTS1    "-E,-T,--ga,--tc,--nc"                             /* incompatible with -E, -T (besides themselves) */
#define CUTOPTS2    "-E,-T,--ga,--tc,--nc,--viterbi,--forward"         /* incompatible with --ga, --tc, --nc (besides themselves) */
#define HMMONLYOPTS "--viterbi,--forward"                              /* with these options, there are no filters, use only HMM */
								     
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  /* basic options */
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-o",        eslARG_OUTFILE,NULL,  NULL, NULL,      NULL,      NULL,        NULL, "direct output to file <f>, not stdout", 1 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "configure CM/HMM for glocal alignment [default: local]", 1 },
  { "--informat",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--toponly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the top strand", 1 },
  { "--bottomonly", eslARG_NONE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the bottom strand", 1 },
  { "--null2",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, HMMONLYOPTS, "turn on the post hoc second null model", 1 },
  /* 4 --p* options below are hopefully temporary b/c if we have E-values for the CM using a certain cm->pbegin, cm->pend,
   * changing those values in cmsearch invalidates the E-values, so we should pick hard-coded values for cm->pbegin cm->pend */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-g,--pbegin","set all local begins as equiprobable", 1 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      NULL, "-g,--pend",  "set all local end probs to <x>", 1 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        "-g", "set aggregate local begin prob to <x>", 1 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        "-g", "set aggregate local end prob to <x>", 1 },
  /* options for algorithm for final round of search */
  { "--inside",  eslARG_NONE,  "default",NULL,NULL,    NULL,       NULL,    STRATOPTS1, "use scanning CM Inside algorithm", 2 },
  { "--cyk",     eslARG_NONE,  FALSE, NULL, NULL, "--inside",      NULL,    STRATOPTS1, "use scanning CM CYK algorithm", 2 },
  { "--viterbi", eslARG_NONE,  FALSE, NULL, NULL, "--fil-hmm,--fil-qdb",     NULL,    STRATOPTS2, "use scanning HMM Viterbi algorithm", 2 },
  { "--forward", eslARG_NONE,  FALSE, NULL, NULL, "--fil-hmm,--fil-qdb",     NULL,    STRATOPTS2, "use scanning HMM Forward algorithm", 2 },
  /* CM cutoff options */
  { "-E",        eslARG_REAL,   "0.1", NULL, "x>0.",    NULL,      NULL,    CUTOPTS1, "use cutoff E-value of <x> for final round of search", 3 },
  { "-T",        eslARG_REAL,   "0.0", NULL, NULL,      NULL,      NULL,    CUTOPTS1, "use cutoff bit score of <x> for final round of search", 3 },
  { "--ga",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    CUTOPTS2, "use CM Rfam GA gathering threshold as cutoff bit score", 3 },
  { "--tc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    CUTOPTS2, "use CM Rfam TC trusted cutoff as cutoff bit score", 3 },
  { "--nc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    CUTOPTS2, "use CM Rfam NC noise cutoff as cutoff bit score", 3 },
  /* banded options (for final round of searching) */
  { "--no-qdb",  eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,  HMMONLYOPTS, "do not use QDBs in final round of searching (after >= 0 filters)", 4 },
  { "--beta",    eslARG_REAL,  "1e-13",NULL, "0<x<1",   NULL,      NULL,  HMMONLYOPTS, "set tail loss prob for QDB calculation to <x>", 4 },
  { "--hbanded", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  HMMONLYOPTS, "calculate and use HMM bands in final round of CM search", 4 },
  { "--tau",     eslARG_REAL,   "1e-7",NULL, "0<x<1",   NULL,"--hbanded", HMMONLYOPTS, "set tail loss prob for --hbanded to <x>", 4 },
  { "--aln2bands",eslARG_NONE, FALSE, NULL, NULL,      NULL, "--hbanded", HMMONLYOPTS, "w/--hbanded derive HMM bands w/o scanning Forward/Backward", 4 },
  /* filtering options, by default do HMM, then CYK filter */
  { "--fil-qdb",   eslARG_NONE, "default", NULL, NULL,  NULL,      NULL,"--fil-no-qdb", "filter with CM QDB (banded) CYK algorithm", 5 },
  { "--fil-beta",  eslARG_REAL, NULL,      NULL, "x>0", NULL,      NULL,"--fil-no-qdb", "set tail loss prob for filter QDB and W calculation to <x>", 4 },
  { "--fil-no-qdb",eslARG_NONE, FALSE,     NULL, NULL,  "--fil-qdb",NULL,         NULL, "do not filter with CM banded CYK", 5 },
  { "--fil-hmm",   eslARG_NONE, "default", NULL, NULL,  NULL,      NULL,"--fil-no-hmm", "filter with HMM forward algorithm", 5 },
  { "--fil-no-hmm",eslARG_NONE, FALSE,     NULL, NULL,  "--fil-hmm",NULL,         NULL, "do not filter with HMM forward algorithm", 5 },
  /* filter cutoff options */
  {"--fil-S-qdb",eslARG_REAL,  "0.02",NULL, "0<x<1.",    NULL,      NULL, "--fil-T-qdb", "set QDB CM filter cutoff to achieve survival fraction <x>", 6 },
  ///{ "--fil-E-qdb",eslARG_REAL,  "0.02",NULL, "x>0.",    NULL,      NULL, "--fil-T-qdb", "use E-value of <x> QDB CM filter", 6 },
  { "--fil-T-qdb",eslARG_REAL,  "0.0", NULL, NULL,      NULL,      NULL, "--fil-S-qdb", "use cutoff bit score of <x> for QDB CM filter", 6 },
  { "--fil-S-hmm",eslARG_REAL,  "0.02",NULL, "0<x<1",    NULL,      NULL, "--fil-T-hmm", "set HMM filter cutoff to achieve survival fraction <x>", 6 },
  ///  { "--fil-E-hmm",eslARG_REAL,  "0.02",NULL, "x>0.",    NULL,      NULL, "--fil-T-hmm", "use E-value cut of <x> for HMM filter", 6 }, 
  { "--fil-T-hmm",eslARG_REAL,  "0.0", NULL, NULL,      NULL,      NULL,"--fil-S-hmm", "use cutoff bit score of <x> for HMM filter", 6 },
  { "--fil-Smax-hmm",eslARG_REAL,NULL, NULL, "0<x<1",    NULL,      NULL,"--fil-T-hmm,--fil-S-hmm", "set maximum HMM survival fraction (predicted) as <x>", 6 },
  /* alignment options */
  { "-p",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "append posterior probabilities to hit alignments", 7 },
  { "--noalign", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,       NULL, "find start/stop/score only; don't do alignments", 7 },
  { "--optacc",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "align hits with the Holmes/Durbin optimal accuracy algorithm", 7 },
  { "--addx",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "add line to output alnments marking non-compensatory bps with 'x'", 7 },
  /* verbose output files */
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump parsetrees for each hit to file <f>", 8 },
  { "--gcfile",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save GC content stats of target sequence file to <f>", 8 },
  { "--bfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save query-dependent bands (QDBs) for each state to file <f>", 8 },
  /* Setting output alphabet */
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as RNA sequence data", 9 },
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 9 },
  /* Other options */
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb", 10 },  
  { "--mxsize",  eslARG_REAL, "256.0", NULL, "x>0.",    NULL,      NULL,        NULL, "set maximum allowable HMM banded DP matrix size to <x> Mb", 10 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program", 10 },  
#endif
  /* Development options */
#ifdef HAVE_DEVOPTS
  { "--rtrans",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--viterbi,--forward", "replace CM transition scores from <cmfile> with RSEARCH scores", 11 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 11 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 * NOTE: MPI not yet implemented.
 */
struct cfg_s {
  char         *cmfile;	        /* name of input CM file  */ 
  char         *sqfile;	        /* name of sequence file  */ 
  ESL_SQFILE   *sqfp;           /* open sequence input file stream */
  FILE         *ofp;            /* output file (default is stdout) */
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  long          dbsize;         /* database size in nucleotides (doubled if doing rev comp) */
  int           ncm;            /* number CM we're at in file */
  int           do_rc;          /* should we search reverse complement? (for convenience */
  int           init_rci;       /* initial strand to search 0 for top, 1 for bottom (only 1 if --bottomonly enabled) */
  float         avg_hit_len;    /* average CM hit length, calc'ed using QDB calculation algorithm */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (mainly i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *tfp;	        /* optional output for parsetrees  */
  FILE         *bfp;	        /* optional output for qdbs */
  FILE         *filhfp;	        /* optional output for filter thr calc histgram */
  FILE         *filrfp;	        /* optional output for filter thr calc R info file */
  ESL_ALPHABET *abc_out; 	/* digital alphabet for writing */
};

static char usage[]  = "[-options] <cmfile> <sequence file>";
static char banner[] = "align sequences to an RNA CM";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
/* static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); */

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int set_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int print_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, FILE *fp, CM_t *cm, long N, char *errbuf);
static int read_next_search_seq(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int do_revcomp, dbseq_t **ret_dbseq);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  struct cfg_s     cfg;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  /*********************************************** 
   * Parse command line
   ***********************************************/

  /* Process command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nalgorithm for final round of search (after >= 0 filters): [default: --inside]");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\ncutoff options for final round of search:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\noptions for banded DP in final round of search:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\nfiltering options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nfilter cutoff options (survival fractions are predicted, not guaranteed):");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\noptions for returning alignments of search hits:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\nverbose output files:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\nother options:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
#ifdef HAVE_DEVOPTS
      puts("\nnon-standard options used only in development:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80);
#endif 
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 2) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\n  where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.sqfile     = esl_opt_GetArg(go, 2); 
  cfg.ofp        = NULL;
  cfg.sqfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  if   (esl_opt_IsDefault(go, "--informat")) cfg.fmt = eslSQFILE_UNKNOWN; /* autodetect sequence file format by default. */ 
  else { 
    cfg.fmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
    if(cfg.fmt == eslSQFILE_UNKNOWN) cm_Fail("Can't recognize sequence file format: %s. valid options are: fasta, embl, genbank, ddbj, uniprot, stockholm, or pfam\n", esl_opt_GetString(go, "--informat"));
  }
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if      (esl_opt_GetBoolean(go, "--rna")) cfg.abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna")) cfg.abc_out = esl_alphabet_Create(eslDNA);
  else    cm_Fail("Can't determine output alphabet");
  cfg.dbsize     = 0;                      /* db size  */
  cfg.ncm        = 0;                      /* what number CM we're on, updated in masters, stays 0 (irrelevant) for workers */
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tfp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.bfp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.filhfp     = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.filrfp     = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.do_rc      = (! esl_opt_GetBoolean(go, "--toponly")); 
  cfg.init_rci   = esl_opt_GetBoolean(go, "--bottomonly") ? 1 : 0; 
  cfg.avg_hit_len= 0.;

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");
  
  /* This is our stall point, if we need to wait until we get a
   * debugger attached to this process for debugging (especially
   * useful for MPI):
   */
  while (cfg.do_stall); 

  /* Start timing. */
  esl_stopwatch_Start(w);

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you execute mpirun?)");

      if (cfg.my_rank > 0)  mpi_worker(go, &cfg);
      else 		    mpi_master(go, &cfg);

      esl_stopwatch_Stop(w);
      esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }
  if (cfg.my_rank == 0) esl_stopwatch_Display(stdout, w, "# CPU time: ");

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (! esl_opt_IsDefault(go, "-o")) { fclose(cfg.ofp); }
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.tfp       != NULL) fclose(cfg.tfp);
    if (cfg.bfp       != NULL) fclose(cfg.bfp);
    if (cfg.filhfp    != NULL) fclose(cfg.filhfp);
    if (cfg.filrfp    != NULL) fclose(cfg.filrfp);
  }
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.abc_out   != NULL) esl_alphabet_Destroy(cfg.abc_out);
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->cmfile      - command line arg 1
 *    cfg->sqfile      - command line arg 2
 *    cfg->fmt         - format of output file
 * Allocates/Sets: 
 *    cfg->sqfp        - open sequence file                
 *    cfg->cmfp        - open CM file                
 *    cfg->tfp         - optional output file
 *    cfg->bfp         - optional output file
 *    cfg->filhfp      - optional output file
 *    cfg->filrfp      - optional output file
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

  /* open output file, or set to stdout if none */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } else cfg->ofp = stdout;

  /* open input sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(cfg->sqfp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", cfg->sqfile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", cfg->sqfile);
  else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
  cfg->fmt = cfg->sqfp->format;

  /* GetDBInfo() reads all sequences, rewinds seq file and returns db size */
  GetDBInfo(NULL, cfg->sqfp, &(cfg->dbsize), NULL);  
  if ((! esl_opt_GetBoolean(go, "--toponly")) && (! esl_opt_GetBoolean(go, "--bottomonly"))) cfg->dbsize *= 2;

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open bands file */
  if (esl_opt_GetString(go, "--bfile") != NULL) {
    if ((cfg->bfp = fopen(esl_opt_GetString(go, "--bfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --bfile output file %s\n", esl_opt_GetString(go, "--bfile"));
    }

  return eslOK;
}

/* serial_master()
 * The serial version of cmsearch.
 * 
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int            status;
  char           errbuf[cmERRBUFSIZE];
  CM_t          *cm = NULL;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */
  int            using_e_cutoff;
  int            rci;
  dbseq_t       *dbseq = NULL;
  int            do_top;


  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/
  do_top = (cfg->init_rci == 0) ? TRUE : FALSE; 

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      /* initialize the flags/options/params and configuration of the CM */
      if((  status = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((  status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) cm_Fail(errbuf);
      if((  status = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = UpdateGumbelsForDBSize(cm, errbuf, cfg->dbsize))            != eslOK) cm_Fail(errbuf);
      if((status = set_searchinfo(go, cfg, errbuf, cm))                     != eslOK) cm_Fail(errbuf);
      print_searchinfo(go, cfg, stdout, cm, cfg->dbsize, errbuf);
      using_e_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF) ? TRUE : FALSE;
	 
      while ((status = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK)
	{
	  for(rci = cfg->init_rci; rci <= cfg->do_rc; rci++) {
	    /*printf("SEARCHING >%s %d\n", dbseq->sq[reversed]->name, reversed);*/
	    if ((status = ProcessSearchWorkunit(cm, errbuf, dbseq->sq[rci]->dsq, dbseq->sq[rci]->n, &dbseq->results[rci], esl_opt_GetReal(go, "--mxsize"), cfg->my_rank)) != eslOK) cm_Fail(errbuf);
	    RemoveOverlappingHits(dbseq->results[rci], 1, dbseq->sq[rci]->n);
	    if(using_e_cutoff) RemoveHitsOverECutoff(cm, cm->si, dbseq->results[rci], dbseq->sq[rci]); 
	  }
	  PrintResults (cm, cfg->ofp, cm->si, cfg->abc_out, cons, dbseq, do_top, cfg->do_rc, esl_opt_GetBoolean(go, "--addx"));
	  for(rci = 0; rci <= cfg->do_rc; rci++) { /* we can free results for top strand even if cfg->init_rci is 1, due to --bottomonly */
	    FreeResults(dbseq->results[rci]);
	    esl_sq_Destroy(dbseq->sq[rci]);
	  }
	  free(dbseq);
	}
      if (status != eslEOF) cm_Fail("Parse failed, line %d, file %s:\n%s", 
				    cfg->sqfp->linenumber, cfg->sqfp->filename, cfg->sqfp->errbuf);
      FreeCM(cm);
      FreeCMConsensus(cons);
      esl_sqio_Rewind(cfg->sqfp); /* we may be searching this file again with another CM */
    }
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmsearch
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include (hopefully) all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * cm_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static void
mpi_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      xstatus       = eslOK;	/* changes from OK on recoverable error */
  int      status;
  int      have_work     = TRUE;	/* TRUE while work remains  */
  int      nproc_working = 0;	        /* number of worker processes working, up to nproc-1 */
  int      wi;          	        /* rank of next worker to get an alignment to work on */
  char    *buf           = NULL;	/* input/output buffer, for packed MPI messages */
  int      bn            = 0;
  int      pos = 1;
  int      using_e_cutoff; 
  int      wi_error = 0;                /* worker index that sent back an error message, if an error occurs */

  CM_t *cm;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */

  int si      = 0;        /* sequence index */
  int si_recv = 1;        /* sequence index of the sequence we've just received results for from a worker */
  /* properties of the workers, indexed 1..wi..nproc-1 */
  int *silist = NULL;     /* [0..wi..nproc-1], the sequence index worker wi is working on */
  int in_rc = FALSE;      /* are we currently on the reverse complement? */
  int *rclist = NULL;     /* [0..wi..nproc-1] 0 if worker wi is searching top strand, 1 if wi is searching bottom strand */
  int rci;                /* index that ranges from 0 to 1 */
  int seqpos = 1;         /* sequence position in the current sequence */
  int *seqposlist = NULL; /* [0..wi..nproc-1] the first position of the sequence that worker wi is searching */
  int len;                /* length of chunk */
  int *lenlist = NULL;    /* [0..wi..nproc-1] length of chunk worker wi is searching */
  /* properties of the sequences currently being worked on, we can have at most 1 per worker, so these are of size 
   * cfg->nproc, but indexed by si = 0..nproc-2, cfg->nproc-1 is never used. */
  int *sentlist = NULL;   /* [0..si..nproc-1] TRUE if all chunks for sequence index si have been sent, FALSE otherwise */
  int ndbseq = 0;         /* ndbseq is the number of currently active sequences, we can read a new seq IFF ndbseq < (cfg->nproc-1) */
  dbseq_t **dbseqlist= NULL; /* pointers to the dbseq_t objects that hold the actual sequence data, and the results data */
  dbseq_t  *dbseq = NULL;  /* a database sequence */
  
  char     errbuf[cmERRBUFSIZE];
  MPI_Status mpistatus; 
  int      n;

  int need_seq = TRUE;
  int chunksize;
  search_results_t *worker_results;

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  /*if (xstatus == eslOK) { if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }*/
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((silist     = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((rclist     = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqposlist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((lenlist    = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((dbseqlist  = malloc(sizeof(dbseq_t *) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((sentlist   = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) cm_Fail(errbuf);
  ESL_DPRINTF1(("MPI master is initialized\n"));

  for (wi = 0; wi < cfg->nproc; wi++) 
  { 
    silist[wi] = rclist[wi] = seqposlist[wi] = lenlist[wi] = -1;
    dbseqlist[wi] = NULL;
    sentlist[wi] = FALSE;
  }
  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  MPI_Bcast(&(cfg->dbsize), 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) cm_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));

  /* Main loop: combining load workers, send/receive, clear workers loops;
   * also, catch error states and die later, after clean shutdown of workers.
   * 
   * When a recoverable error occurs, have_work = FALSE, xstatus !=
   * eslOK, and errbuf is set to an informative message. No more
   * errbuf's can be received after the first one. We wait for all the
   * workers to clear their work units, then send them shutdown signals,
   * then finally print our errbuf and exit.
   * 
   * Unrecoverable errors just crash us out with cm_Fail().
   */

  while (xstatus == eslOK && CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      cfg->ncm++;  
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((status   = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) cm_Fail(errbuf);
      if((status   = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = UpdateGumbelsForDBSize(cm, errbuf, cfg->dbsize))            != eslOK) cm_Fail(errbuf);
      if((status = set_searchinfo(go, cfg, errbuf, cm))                     != eslOK) cm_Fail(errbuf);

      print_searchinfo(go, cfg, stdout, cm, cfg->dbsize, errbuf);
      using_e_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF) ? TRUE : FALSE;

      /* reset vars for searching with current CM */
      wi = 1;
      ndbseq = 0;
      need_seq = TRUE;
      have_work = TRUE;	/* TRUE while work remains  */
      seqpos = 1;
      in_rc = FALSE;
      while (have_work || nproc_working)
	{
	  if (need_seq) 
	    {
	      need_seq = FALSE;
	      /* read a new seq */
	      if((status = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK) 
		{
		  ndbseq++;
		  ESL_DASSERT1((ndbseq < cfg->nproc));

		  dbseq->chunks_sent = 0;
		  dbseq->alignments_sent = -1;     /* None sent yet */
		  for(rci = 0; rci <= cfg->do_rc; rci++) {
		    dbseq->results[rci] = CreateResults(INIT_RESULTS);
		  }
		  in_rc = (cfg->init_rci == 0) ? FALSE : TRUE; /* if --bottomonly --> cfg->init_rci = 1, and we only search bottom strand */
		  seqpos = 1;
		  
		  si = 0;
		  while(dbseqlist[si] != NULL) si++;
		  ESL_DASSERT1((si < cfg->nproc));
		  dbseqlist[si] = dbseq;
		  sentlist[si]  = FALSE;
		  have_work = TRUE;
		  chunksize = DetermineSeqChunksize(cfg->nproc, dbseq->sq[0]->n, cm->W);
		  ESL_DPRINTF1(("L: %d chunksize: %d\n", dbseq->sq[0]->n, chunksize));
		}
	      else if(status == eslEOF) have_work = FALSE;
	      else goto ERROR;
	    }
	
	  if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0))
	    {
	      /* we're waiting to receive */
	      if (MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi probe failed");
	      if (MPI_Get_count(&mpistatus, MPI_PACKED, &n)                != 0) cm_Fail("mpi get count failed");
	      wi = mpistatus.MPI_SOURCE;
	      ESL_DPRINTF1(("MPI master sees a result of %d bytes from worker %d\n", n, wi));
	      
	      if (n > bn) {
		if ((buf = realloc(buf, sizeof(char) * n)) == NULL) cm_Fail("reallocation failed");
		bn = n; 
	      }
	      if (MPI_Recv(buf, bn, MPI_PACKED, wi, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi recv failed");
	      ESL_DPRINTF1(("MPI master has received the buffer\n"));
	      
	      /* If we're in a recoverable error state, we're only clearing worker results;
	       * just receive them, don't unpack them or print them.
	       * But if our xstatus is OK, go ahead and process the result buffer.
	       */
	      if (xstatus == eslOK) /* worker reported success. Get the result. */
		{
		  pos = 0;
		  if (MPI_Unpack(buf, bn, &pos, &xstatus, 1, MPI_INT, MPI_COMM_WORLD)     != 0)     cm_Fail("mpi unpack failed");
		  if (xstatus == eslOK) /* worker reported success. Get the results. */
		    {
		      si_recv = silist[wi];
		      ESL_DPRINTF1(("MPI master sees that the result buffer contains search results (si_recv:%d)\n", si_recv));
		      if ((status = cm_search_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_results)) != eslOK) cm_Fail("search results unpack failed");
		      ESL_DPRINTF1(("MPI master has unpacked search results\n"));
		      
		      /* worker_results will be NULL if 0 results (hits) sent back */
		      int x;
		      if(worker_results != NULL) { 
			/* add results to dbseqlist[si_recv]->results[rclist[wi]] */
			if(! esl_opt_GetBoolean(go, "--noalign")) { 
			  for(x = 0; x < worker_results->num_results; x++) {
			    assert(worker_results->data[x].tr != NULL);
			    assert(worker_results->data[x].tr->n > 0);
			  }
			}
			AppendResults(worker_results, dbseqlist[si_recv]->results[rclist[wi]], seqposlist[wi]);
			/* careful, dbseqlist[si_recv]->results[rclist[wi]] now points to the nodes in worker_results->data,
			 * don't free those (don't use FreeResults(worker_results)) */
			free(worker_results);
			worker_results = NULL;
		      }
		      dbseqlist[si_recv]->chunks_sent--;
		      if(sentlist[si_recv] && dbseqlist[si_recv]->chunks_sent == 0)
			{
			  for(rci = 0; rci <= cfg->do_rc; rci++) {
			    RemoveOverlappingHits(dbseqlist[si_recv]->results[rci], 1, dbseqlist[si_recv]->sq[rci]->n);
			    if(using_e_cutoff) RemoveHitsOverECutoff(cm, cm->si, dbseqlist[si_recv]->results[rci], dbseqlist[si_recv]->sq[rci]);
			  }					      
			  PrintResults(cm, cfg->ofp, cm->si, cfg->abc_out, cons, dbseqlist[si_recv], TRUE, cfg->do_rc, esl_opt_GetBoolean(go, "--addx"));
			  for(rci = 0; rci <= cfg->do_rc; rci++) {
			    esl_sq_Destroy(dbseqlist[si_recv]->sq[rci]);
			    FreeResults(dbseqlist[si_recv]->results[rci]);
			  }
			  free(dbseqlist[si_recv]);
			  dbseqlist[si_recv] = NULL;
			  ndbseq--;
			}
		    }
		  else	/* worker reported an error. Get the errbuf. */
		    {
		      if (MPI_Unpack(buf, bn, &pos, errbuf, cmERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack of errbuf failed");
		      ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
		      have_work = FALSE;
		      wi_error  = wi;
		    }
		}
	      nproc_working--;
	    }
	  
	  if (have_work)
	    {   
	      /* send new search job */
	      len = (chunksize < (dbseqlist[si]->sq[0]->n - seqpos + 1)) ? chunksize : (dbseqlist[si]->sq[0]->n - seqpos + 1);
	      ESL_DPRINTF1(("MPI master is sending sequence i0..j0 %d..%d to search to worker %d\n", seqpos, seqpos+len-1, wi));
	      assert(seqpos > 0);
	      if ((status = cm_dsq_MPISend(dbseqlist[si]->sq[in_rc]->dsq+seqpos-1, len, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
	      
	      silist[wi]      = si;
	      seqposlist[wi]  = seqpos;
	      lenlist[wi]     = len;
	      rclist[wi]      = in_rc;
	      dbseqlist[si]->chunks_sent++;
	      
	      wi++;
	      nproc_working++;
	      
	      if(len == chunksize) seqpos += len - cm->W + 1;
	      else if(cfg->do_rc && !in_rc) {
		in_rc = TRUE;
		seqpos = 1; 
	      }
	      else {
		need_seq     = TRUE;
		sentlist[si] = TRUE; /* we've sent all chunks from this seq */
	      }
	    }
	}
      ESL_DPRINTF1(("MPI master: done with this CM. Telling all workers\n"));
      for (wi = 1; wi < cfg->nproc; wi++) 
	if ((status = cm_dsq_MPISend(NULL, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");
      FreeCM(cm);
      FreeCMConsensus(cons);
      esl_sqio_Rewind(cfg->sqfp); /* we may be searching this file again with another CM */
    }
  
  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((status = cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); cm_Fail(errbuf); }
  else                  return;

 ERROR: 
  cm_Fail("memory allocation error.");
  return; /* NOTREACHED */
}


static void
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int           xstatus = eslOK;
  int           status;
  /*int           type;*/
  CM_t         *cm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
  int           pos;
  char          errbuf[cmERRBUFSIZE];
  /*float         Smin;*/
  /*MPI_Status  mpistatus;*/
  ESL_DSQ      *dsq = NULL;
  int           L;
  search_results_t *results = NULL;

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
  
  /* Master now broadcasts worker initialization information (db size N) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  MPI_Bcast(&(cfg->dbsize), 1, MPI_LONG, 0, MPI_COMM_WORLD);
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized N: %ld\n", cfg->my_rank, cfg->dbsize));
  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      if((status   = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) goto ERROR;
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = UpdateGumbelsForDBSize(cm, errbuf, cfg->dbsize))       != eslOK) goto ERROR;
      if((status = set_searchinfo(go, cfg, errbuf, cm))                     != eslOK) goto ERROR;
      
      /* print_searchinfo(go, cfg, stdout, cm, cm_mode, cp9_mode, cfg->dbsize, errbuf); */
      
      while((status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &dsq, &L)) == eslOK)
	{
	  ESL_DPRINTF1(("worker %d: has received search job, length: %d\n", cfg->my_rank, L));
	  if ((status = ProcessSearchWorkunit(cm, errbuf, dsq, L, &results, esl_opt_GetReal(go, "--mxsize"), cfg->my_rank)) != eslOK) goto ERROR;
	  ESL_DPRINTF1(("worker %d: has gathered search results\n", cfg->my_rank));
	  
	  n = 0;
	  if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	  n += sz;
	  if (cm_search_results_MPIPackSize(results, MPI_COMM_WORLD, &sz) != eslOK)
	    ESL_XFAIL(eslFAIL, errbuf, "cm_serch_results_MPIPackSize() call failed"); 
	  n += sz;  

	  if (n > wn) {
	    void *tmp;
	    ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	    wn = n;
	  }
	  ESL_DPRINTF1(("worker %d: has calculated the search results will pack into %d bytes\n", cfg->my_rank, n));
	  status = eslOK;

	  pos = 0;
	  if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	  if (cm_search_results_MPIPack(results, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
	    ESL_XFAIL(eslFAIL, errbuf, "cm_search_results_MPIPack() call failed"); 
	  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	  ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));

	  FreeResults(results);
	  free(dsq);
	}
      if(status == eslEOD)ESL_DPRINTF1(("worker %d: has seen message to stop with this CM.\n", cfg->my_rank));
      else ESL_XFAIL(eslFAIL, errbuf, "within CM loop, unexpected status code: %d received from cm_dsq_MPIRecv()\n", status);

      FreeCM(cm);
      cm = NULL;
    }
  if (status == eslEOD) ESL_DPRINTF1(("Worker %d told CMs are done.\n", cfg->my_rank));
  else ESL_XFAIL(eslFAIL, errbuf, "outside CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);
  
  if (wbuf != NULL) free(wbuf);
  return;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

  /* if we get here this worker failed and sent an error message, now the master knows a worker
   * failed but it has to send the message to all other workers (besides this one) to abort so they 
   * can be shut down cleanly. As currently implemented, this means we have to wait here for that 
   * signal which comes in the form of a special 'empty' work packet that tells us we're done with
   * the current CM, and then a 'empty' CM broadcast that tells us we're done with all CMs in the file.
   */
  status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &dsq, &L);
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);
  /* status after each of the above calls should be eslEOD, but if it isn't we can't really do anything 
   * about it b/c we've already sent our error message, so in that scenario the MPI will break uncleanly 
   */

  return;
}
#endif /*HAVE_MPI*/

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM. 
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int use_hmmonly;
  int nstarts, nexits, nd;
  
  /* set up CM parameters that are option-changeable */
  cm->beta_qdb = esl_opt_GetReal(go, "--beta");
  cm->tau      = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  use_hmmonly = (esl_opt_GetBoolean(go, "--viterbi") || esl_opt_GetBoolean(go, "--forward")) ? TRUE : FALSE;

  /* Update cm->config_opts and cm->align_opts based on command line options */

  /* config_opts */
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  /* config QDB for final round of search? yes, unless --no-qdb, otherwise we only use QDB to filter */
  if(! esl_opt_GetBoolean(go, "--no-qdb")) cm->config_opts |= CM_CONFIG_QDB;

  /* search_opts */
  if(! use_hmmonly) 
    if(  esl_opt_GetBoolean(go, "--inside"))      cm->search_opts |= CM_SEARCH_INSIDE;
  if(  esl_opt_GetBoolean(go, "--noalign"))     cm->search_opts |= CM_SEARCH_NOALIGN;
  if(  esl_opt_GetBoolean(go, "--null2"))       cm->search_opts |= CM_SEARCH_NULL2;
  if(  esl_opt_GetBoolean(go, "--no-qdb"))      cm->search_opts |= CM_SEARCH_NOQDB;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts |= CM_SEARCH_HBANDED;
  if(  esl_opt_GetBoolean(go, "--aln2bands"))   cm->search_opts |= CM_SEARCH_HMMALNBANDS;
  if(  esl_opt_GetBoolean(go, "--viterbi"))  { 
    cm->search_opts |= CM_SEARCH_HMMVITERBI;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }
  if(  esl_opt_GetBoolean(go, "--forward"))  { 
    cm->search_opts |= CM_SEARCH_HMMFORWARD;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }

  /* align_opts, by default, align with HMM bands */
  cm->align_opts |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--optacc"))       cm->align_opts |= CM_ALIGN_OPTACC;
  if(esl_opt_GetBoolean(go, "-p"))             cm->align_opts |= CM_ALIGN_POST;

  /* handle special developer's options, not recommend for normal users */
#ifdef HAVE_DEVOPTS
  if(  esl_opt_GetBoolean(go, "--rtrans"))      cm->flags       |= CM_RSEARCHTRANS;
  if(  esl_opt_GetBoolean(go, "--sums"))        cm->search_opts |= CM_SEARCH_SUMS;
#endif

  /* TEMPORARY BLOCK */
  /* set aggregate local begin/end probs, set with --pbegin, --pend, defaults are DEFAULT_PBEGIN, DEFAULT_PEND */
  cm->pbegin = esl_opt_GetReal(go, "--pbegin");
  cm->pend   = esl_opt_GetReal(go, "--pend");
  /* possibly overwrite local begin probs such that all begin points are equiprobable (--pebegin) */
  if(esl_opt_GetBoolean(go, "--pebegin")) {
    nstarts = 0;
    for (nd = 2; nd < cm->nodes; nd++) 
      if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
	nstarts++;
    cm->pbegin = 1.- (1./(1+nstarts));
  }
  /* possibly overwrite cm->pend so that local end prob from all legal states is fixed,
   * this is strange in that cm->pend may be placed as a number greater than 1., this number
   * is then divided by nexits in ConfigLocalEnds() to get the prob for each v --> EL transition,
   * this is guaranteed by the way we calculate it to be < 1.,  it's the argument from --pfend */
  if(! esl_opt_IsDefault(go, "--pfend")) {
    nexits = 0;
    for (nd = 1; nd < cm->nodes; nd++) {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	nexits++;
    }
    cm->pend = nexits * esl_opt_GetReal(go, "--pfend");
  }
  /* END TEMPORARY BLOCK */

  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  ConfigCM(cm, TRUE);  /* TRUE says: calculate W */

  /* Setup ScanMatrix for CYK/Inside scanning functions, we can't 
   * do it in initialize_cm(), b/c it's W dependent; W was just set.
   * We don't need it if we're only using an HMM though.
   */
  if(use_hmmonly) cm->smx = NULL;
  else { 
    int do_float = TRUE;
    int do_int   = FALSE;
    if(cm->search_opts & CM_SEARCH_INSIDE) { do_float = FALSE; do_int = TRUE; }
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
    if(cm->smx == NULL) cm_Fail("initialize_cm(), use_hmmonly is FALSE, CreateScanMatrixForCM() call failed, mx is NULL.");
  }

  ESL_DPRINTF1(("cm->pbegin: %.3f\n", cm->pbegin));
  ESL_DPRINTF1(("cm->pend: %.3f\n", cm->pend));
  /* print qdbs to file if nec */
  if(! esl_opt_IsDefault(go, "--bfile")) {
    fprintf(cfg->bfp, "beta:%f\n", cm->beta_qdb);
    debug_print_bands(cfg->bfp, cm, cm->dmin, cm->dmax);
  }

  if(cfg->my_rank == 0) fprintf(cfg->ofp, "CM %d: %s\n", (cfg->ncm), cm->name);
  return eslOK;
}

/* Function: set_searchinfo()
 * Date:     EPN, Mon Jan 21 08:56:04 2008 (updated)
 * 
 * Purpose:  Determine how many rounds of searching we will do (all rounds 
 *           but last round are filters), and set the relevant info in the
 *           SearchInfo_t <cm->si> object, including cutoffs.
 *
 * Filters:
 *
 * User can specify 0 or 1 round of HMM filtering with Forward, followed by 0 or 
 * 1 round of CM filtering with QDB CYK. This is determined as described below. 
 * The default behavior depends on whether or not the CM file has gumbel and 
 * filter threshold information from cmcalibrate.
 * 
 * HMM filtering:
 * ('none' below means none of --fil-no-hmm, --fil-T-hmm, --fil-S-hmm, --fil-Smax-hmm) 
 * 
 *                    gumbel &
 *                  filter thr 
 *   options        in cmfile?    filter/don't filter and cutoff determination
 *   -------       ----------    ------------------------------------------------
 * 1. none                 yes    filter (usually). automatically determine appropriate HMM 
 *                                filter cutoff for final round cutoff using cmfile's filter 
 *                                threshold info from cmcalibrate. Sometimes this info may
 *                                tell us it's inappropriate to use an HMM filter, in which
 *                                case, don't filter.
 * 
 * 2. none                  no    filter. set filter bit score cutoff to default value <x> for 
 *                                --fil-T-hmm (3.0 bits)
 *                               
 * 3. --fil-no-hmm      yes/no    don't filter.
 * 
 * 4. --fil-T-hmm <x>   yes/no    filter. set filter bit score cutoff as <x>
 * 
 * 5. --fil-S-hmm <x>      yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 6. --fil-S-hmm <x>       no    die. we can't deal with this combo, tell the user.
 * 
 * 7. --fil-Smax-hmm <x>   yes    filter, same as in case 1 above, but if the 'appropriate'
 *                                CM E-value cutoff gives a predicted survival fraction > <x>,
 *                                reset it as the E-value cutoff that gives <x> (even if 
 *                                cmfile told us to turn filtering off).
 * 
 * 8. --fil-Smax-hmm <x>    no    die. we can't deal with this combo, tell the user.
 *
 * 
 * CM filtering: 
 * ('none' below means none of --fil-no-qdb, --fil-T-qdb, --fil-S-qdb) 
 *                               
 *                                final 
 *                    gumbel &    round 
 *                  filter thr    cutoff
 *   options        in cmfile?    type    filter/don't filter and cutoff determination
 *   -------       -----------    ------  -----------------------------------------
 * 1. none                 yes    E val   filter. set E-value filter cutoff to E-value that gives
 *                                        predicted survival fraction of <x>. NOTE: larger models,
 *                                        with larger Ws, will have smaller E-value cutoffs for
 *                                        same <x> (a somewhat unsettling behavior).
 * 
 *HERE HERE HERE 
 *                                final 
 *                    gumbel &    round 
 *                  filter thr    cutoff
 *   options        in cmfile?    type    filter/don't filter and cutoff determination
 *   -------       -----------    ------  -----------------------------------------
 * 1. none                 yes    E val   filter. set E-value filter cutoff to E-value that gives
 *                                        predicted survival fraction of <x>. NOTE: larger models,
 *                                        with larger Ws, will have smaller E-value cutoffs for
 *                                        same <x> (a somewhat unsettling behavior).

 * 2. none                  no    filter. set filter bit score cutoff to default value <x> for 
 *                                --fil-T-qdb (0.0 bits)
 *                               
 * 3. --fil-no-qdb      yes/no    don't filter. 
 * 
 * 4. --fil-T-qdb <x>   yes/no    filter. set filter bit score cutoff as <x>
 * 
 * 5. --fil-S-qdb <x>      yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 6. --fil-S-qdb <x>       no    die. we can't deal with this combo, tell the user.
 *
 * ('none' below means none of --fil-no-qdb, --fil-T-qdb, --fil-S-qdb) 
 * 
 *                    gumbel &
 *                  filter thr 
 *   options        in cmfile?    filter/don't filter and cutoff determination
 *   -------       -----------    ------------------------------------------------
 * 1. none                 yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 2. none                  no    filter. set filter bit score cutoff to default value <x> for 
 *                                --fil-T-qdb (0.0 bits)
 *                               
 * 3. --fil-no-qdb      yes/no    don't filter. 
 * 
 * 4. --fil-T-qdb <x>   yes/no    filter. set filter bit score cutoff as <x>
 * 
 * 5. --fil-S-qdb <x>      yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 6. --fil-S-qdb <x>       no    die. we can't deal with this combo, tell the user.
 *
 * 
 * beta, the tail loss probability for the QDB calculation is determined as follows:
 * If --fil-beta <x> is not enabled, cm->beta_qdb is set as cm->beta_W as read in the cmfile
 *               and initially set by cmbuild.
 * If --fil-beta <x> is enabled, <x> is used as cm->beta_qdb
 *
 * Final round related options (after all filtering is complete):
 *
 * --inside:     search with CM inside (TRUE by default)
 * --cyk:        search with CM CYK 
 * --forward:    search with HMM Forward
 * --viterbi:    search with HMM Viterbi
 * -T:           CM/HMM bit score threshold
 * -E:           CM/HMM E-value threshold (requires gumbel info in CM file)
 * --ga:         use Rfam gathering threshold (bit sc) from CM file 
 * --tc:         use Rfam trusted cutoff      (bit sc) from CM file
 * --nc:         use Rfam noise cutoff        (bit sc) from CM file
 *
 * NOTE: --ga, --nc, --tc are incompatible with --forward and --viterbi
 *
 */
static int
set_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;                 /* easel status code */
  int n;                      /* counter over rounds */
  int stype;                  /* type of filter */
  int search_opts;            /* search_opts for filter */
  int use_hmmonly;            /* TRUE if --viterbi or --forward */
  int we_have_gumbels;        /* TRUE if cm->flags & CMH_GUMBEL_STATS */
  int we_have_fthr;           /* TRUE if cm->flags & CMH_FILTER_STATS */
  int fthr_mode;              /* filter threshold mode */
  float final_E  = -1.;       /* final round E-value cutoff, stays -1 if !CMH_GUMBEL_STATS */
  float fqdb_E   = -1.;       /* QDB filter round E-value cutoff, stays -1 if !CMH_GUMBEL_STATS */
  float fhmm_E   = -1.;       /* HMM filter round E-value cutoff, stays -1 if !CMH_GUMBEL_STATS */
  float final_sc = -1.;       /* final round bit score cutoff */
  float fqdb_sc  = -1.;       /* QDB filter round bit score cutoff */
  float fhmm_sc  = -1.;       /* HMM filter round bit score cutoff */
  int   final_ctype;          /* final round cutoff type, SCORE_CUTOFF or E_CUTOFF */
  int   fqdb_ctype;           /* QDB filter round cutoff type, SCORE_CUTOFF or E_CUTOFF */
  int   fhmm_ctype;           /* HMM filter round cutoff type, SCORE_CUTOFF or E_CUTOFF */
  float final_S = -1;         /* predicted survival fraction from final round */
  float fqdb_S = -1;          /* predicted survival fraction from qdb filter round */
  float fhmm_S = -1;          /* predicted survival fraction from HMM filter round */
  int   do_qdb_filter = TRUE; /* TRUE to add QDB filter, FALSE not to */
  int   do_hmm_filter = TRUE; /* TRUE to add HMM filter, FALSE not to */
  double fqdb_beta_qdb;       /* beta for QDBs in QDB filter round */
  double fqdb_beta_W;         /* beta for W in QDB filter round */
  int    fqdb_W;              /* W for QDB filter round */
  int   *fqdb_dmin, *fqdb_dmax; /* d bands (QDBs) for QDB filter round */
  int safe_windowlen;         /* used to get QDBs */
  ScanMatrix_t *fqdb_smx = NULL;/* the scan matrix for the QDB filter round */
  int cm_mode, hmm_mode;      /* CM gumbel mode and CP9 HMM gumbel mode for gumbel statistics */
  int qdb_mode;               /* CM gumbel mode during QDB filter round */
  int cut_point;              /* HMM forward E-value cut point from filter threshold stats */

  if(cm->si != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo(), cm->si is not NULL, shouldn't happen.\n");

  /* Create SearchInfo, specifying no filtering, we change the threshold below */
  CreateSearchInfo(cm, SCORE_CUTOFF, 0., -1.);
  if(cm->si == NULL) cm_Fail("set_searchinfo(), CreateSearchInfo() call failed.");
  SearchInfo_t *si = cm->si; 
  if(si->nrounds > 0) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_search_info(), si->nrounds (%d) > 0\n", si->nrounds);
  
  /* First, set up cutoff for final round, this will be round si->nrounds == 0, b/c no filters have been added yet */
  n           = si->nrounds;
  stype       = si->stype[n];
  search_opts = si->search_opts[n];
  we_have_gumbels = (cm->flags & CMH_GUMBEL_STATS);
  we_have_fthr    = (cm->flags & CMH_FILTER_STATS);
  use_hmmonly = ((search_opts & CM_SEARCH_HMMVITERBI) || (search_opts & CM_SEARCH_HMMFORWARD));
  if(use_hmmonly) do_hmm_filter = do_qdb_filter = FALSE; /* don't filter if we're searching only with the HMM */
  if((! use_hmmonly) && (stype != SEARCH_WITH_CM)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), search_opts for final round of search does not have HMMVITERBI or HMMFORWARD flags raised, but is not of type SEARCH_WITH_CM.");
  /* determine configuration of CM and CP9 HMM based on cm->flags & cm->search_opts */
  CM2Gumbel_mode(cm, search_opts, &cm_mode, &hmm_mode); 

  final_S = final_E = final_sc = -1.;
  /* set up final round cutoff, either 0 or 1 of 5 options is enabled. 
   * esl_opt_IsDefault() returns FALSE even if option is enabled with default value.
   */
  if(esl_opt_IsDefault(go, "-E") && 
     esl_opt_IsDefault(go, "-T") && 
     esl_opt_IsDefault(go, "--ga") && 
     esl_opt_IsDefault(go, "--tc") && 
     esl_opt_IsDefault(go, "--nc")) { 
    /* Choose from, in order of priority:
     * 1. default E value if CM file has Gumbel stats
     * 3. default bit score
     */
    if(we_have_gumbels) { /* use default E-value cutoff */
      final_ctype = E_CUTOFF;
      final_E     = esl_opt_GetReal(go, "-E");
      if((status = E2Score(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E, &final_sc)) != eslOK) return status;
    }
    else { /* no Gumbel stats in CM file, use default bit score cutoff */
      final_ctype = SCORE_CUTOFF;
      final_sc    = esl_opt_GetReal(go, "-T");
      final_E     = -1.; /* invalid, we'll never use it */  
    }
  }
  else if(! esl_opt_IsDefault(go, "-E")) { /* -E enabled, use that */
    if(!we_have_gumbels) ESL_FAIL(eslEINVAL, errbuf, "-E requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
    final_ctype = E_CUTOFF;
    final_E     = esl_opt_GetReal(go, "-E");
    if((status = E2Score(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E, &final_sc)) != eslOK) return status;
  }
  else if(! esl_opt_IsDefault(go, "-T")) { /* -T enabled, use that */
    final_ctype = SCORE_CUTOFF;
    final_sc    = esl_opt_GetReal(go, "-T");
    final_E     = -1.; /* invalid, we'll never use it */  
  }
  else if(! esl_opt_IsDefault(go, "--ga")) { /* --ga enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--ga is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_GA)) ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->ga;
    final_E     = -1.; /* we'll never use it */
  }
  else if(! esl_opt_IsDefault(go, "--tc")) { /* --tc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--tc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_TC)) ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->tc;
    final_E     = -1.; /* we'll never use it */
  }
  else if(! esl_opt_IsDefault(go, "--nc")) { /* --nc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--nc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_NC)) ESL_FAIL(eslEINVAL, errbuf, "No NC noise cutoff in CM file, can't use --nc.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->nc;
    final_E     = -1.; /* we'll never use it */
  }

  else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No final round cutoff selected. This shouldn't happen.");

  if(we_have_gumbels) { 
    if(final_ctype == SCORE_CUTOFF) { /* determine max E-value that corresponds to final_sc bit sc cutoff  across all partitions */
      if((status = Score2E(cm, errbuf, hmm_mode, final_sc, &final_E)) != eslOK) return status;
    }
    else if(final_ctype == E_CUTOFF) { /* determine min bit sc that corresponds to final_E E-val cutoff across all partitions */
      if((status  = E2Score(cm, errbuf, hmm_mode, final_E, &final_sc)) != eslOK) return status;
    }
    final_S = E2SurvFract(final_E, cm->W, cfg->avg_hit_len, cfg->dbsize);
  }
  /* update the search info, which holds the thresholds for final round */
  UpdateSearchInfoCutoff(cm, cm->si->nrounds, final_ctype, final_sc, final_E);   
  ValidateSearchInfo(cm, cm->si);
  /* DumpSearchInfo(cm->si); */

  if(we_have_gumbels) { 
    final_S = E2SurvFract(final_E, cm->W, cfg->avg_hit_len, cfg->dbsize);
    if(final_S >= 0.09) do_qdb_filter = FALSE; /* we want QDB filter to let through 10X survival fraction of final round, so if final round S > 0.09, qdb filter would only filter out 10% of database, so we turn it off */
  }
  else final_S = final_E = -1.;

  /* done with threshold for final round */

  /* Set up the filters and their thresholds 
   * A. determine thresholds/stats for HMM filter to add
   * B. determine thresholds/stats for CM QDB filter to add
   * C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first)
   * D. add HMM filter, if necessary (after QDB filter)
   */

  /* 1. determine thresholds/stats for HMM filter */
  do_hmm_filter = esl_opt_GetBoolean(go, "--fil-hmm");
  fhmm_S = fhmm_E = fhmm_sc = -1.;
  if(do_hmm_filter) { /* determine thresholds for HMM forward filter */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), --fil-hmm enabled, along with --viterbi or --forward, shouldn't happen.");
    if(esl_opt_IsDefault(go, "--fil-S-hmm") && esl_opt_IsDefault(go, "--fil-T-hmm")) { /* default: use HMM filter threshold stats, if they exist in cmfile, else use default bit score cutoff */
      /* No relevant options selected. Choose from, in order of priority:
       * 1. appropriate HMM filter E value cutoff from cmfile (if cmfile has filter threshold stats)
       * 2. default HMM filter bit score
       */
      if(we_have_fthr) { /* we have filter stats in cmfile, use appropriate HMM E-value cutoff */
	/* determine filter threshold mode, the mode of final stage of searching, either FTHR_CM_LC,
	 * FTHR_CM_LI, FTHR_CM_GC, FTHR_CM_GI (can't be an HMM mode b/c --viterbi and --forward toggle --fil-hmm off)
	 */
	if((status = CM2FthrMode(cm, errbuf, cm->search_opts, &fthr_mode)) != eslOK) return status;
	HMMFilterInfo_t *hfi_ptr = cm->stats->hfiA[fthr_mode]; /* for convenience */
	if(hfi_ptr->is_valid == FALSE) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), --fgiven enabled, cm's CMH_FILTER_STATS is raised, but best filter info for fthr_mode %d is invalid.", fthr_mode);
	fhmm_ctype = E_CUTOFF;
	/* determine the appropriate filter cut point <cut_point> to use, for details
	 * see comments in function: searchinfo.c:GetHMMFilterFwdECutGivenCME() for details */
	if(final_ctype == SCORE_CUTOFF) { /* final round has bit score cutoff */
	  if((status = GetHMMFilterFwdECutGivenCMBitScore(hfi_ptr, errbuf, final_sc, cfg->dbsize, &cut_point, cm, cm_mode)) != eslOK) return status; 
	}
	else { /* final round has E-value cutoff */
	  if((status = GetHMMFilterFwdECutGivenCME(hfi_ptr, errbuf, final_E, cfg->dbsize, &cut_point)) != eslOK) return status; 
	}
	if(cut_point != -1) { /* it's worth it to filter */
	  fhmm_E = hfi_ptr->fwd_E_cut[cut_point] * ((double) cfg->dbsize / (double) hfi_ptr->dbsize); 
	  fhmm_S = E2SurvFract(fhmm_E, cm->W, cfg->avg_hit_len, cfg->dbsize);
	  /* check if --fil-Smax-hmm applies */
	  if(! esl_opt_IsDefault(go, "--fil-Smax-hmm")) { 
	    if(fhmm_S > esl_opt_GetReal(go, "--fil-Smax-hmm")) { /* predicted survival fraction exceeds maximum allowed, set E cutoff as E value that gives max allowed survival fraction */
	      fhmm_E = SurvFract2E(esl_opt_GetReal(go, "--fil-Smax-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
	    }
	  }
	  if((status  = E2Score(cm, errbuf, cm_mode, fhmm_E, &fhmm_sc)) != eslOK) return status; /* note: use cm_mode, not fthr_mode */
	}
	else { 
	  do_hmm_filter = FALSE;
	  /* it's not worth it to filter, our HMM filter cutoff would be so low, 
	   * letting so much of the db survive, the filter is a waste of time */
	  ESL_DPRINTF1(("cut_point -1, always_better FALSE\n"));
	}
      }
      else { /* no HMM filter relevant options selected, and cmfile does not have HMM filter stats, use default bit score */
	fhmm_ctype = SCORE_CUTOFF;
	fhmm_sc    = esl_opt_GetReal(go, "--fil-T-hmm");
      }
    } /* end of if(esl_opt_IsDefault(go, "--fil-S-hmm") && esl_opt_IsDefault(go, "--fil-T-hmm")),
       * if we get here, we have no HMM filter relevant options selected. */
    else if(! esl_opt_IsDefault(go, "--fil-S-hmm")) {
      if(!we_have_gumbels) ESL_FAIL(eslEINVAL, errbuf, "--fil-S-hmm requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      fhmm_ctype = E_CUTOFF;
      fhmm_E     = SurvFract2E(esl_opt_GetReal(go, "--fil-S-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
      fhmm_E     = ESL_MIN(fhmm_E, 1.0); /* we don't allow filter E cutoffs below 1. */
    }
    else if(! esl_opt_IsDefault(go, "--fil-T-hmm")) {
      fhmm_ctype = SCORE_CUTOFF;
      fhmm_sc    = esl_opt_GetReal(go, "--fil-T-hmm");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No HMM filter cutoff selected. This shouldn't happen.");

    if(we_have_gumbels) { 
      if(fhmm_ctype == SCORE_CUTOFF) { /* determine max E-value that corresponds to fhmm_sc bit sc cutoff  across all partitions */
	if((status = Score2E(cm, errbuf, hmm_mode, fhmm_sc, &fhmm_E)) != eslOK) return status;
      }
      else if(fhmm_ctype == E_CUTOFF) { /* determine min bit sc that corresponds to fhmm_E E-val cutoff across all partitions */
	if((status  = E2Score(cm, errbuf, hmm_mode, fhmm_E, &fhmm_sc)) != eslOK) return status;
      }
      fhmm_S = E2SurvFract(fhmm_E, cm->W, cfg->avg_hit_len, cfg->dbsize);
    }
  }

  /* B. determine thresholds/stats for CM QDB filter to add (do this after HMM stats b/c we use fhmm_S when setting fqdb_E */
  qdb_mode = (GumModeIsLocal(cm_mode)) ? GUM_CM_LC : GUM_CM_GC; /* always do CYK with QDB filter, only question is local or glocal? */
  fqdb_S = fqdb_E = fqdb_sc = -1.;
  
  if(do_qdb_filter && esl_opt_GetBoolean(go, "--fil-qdb")) { /* determine thresholds, beta for qdb cyk filter */
    /* build the ScanMatrix_t for the QDB filter round, requires calcing dmin, dmax */
    if(! esl_opt_IsDefault(go, "--fil-beta")) fqdb_beta_qdb = esl_opt_GetReal(go, "--fil-beta");
    else fqdb_beta_qdb = cm->beta_W; /* use beta used to calc W in CM file by default */
    safe_windowlen = cm->W * 3;
    while(!(BandCalculationEngine(cm, safe_windowlen, fqdb_beta_qdb, FALSE, &fqdb_dmin, &fqdb_dmax, NULL, NULL))) {
      free(fqdb_dmin);
      free(fqdb_dmax);
      fqdb_dmin = NULL;
      fqdb_dmax = NULL;
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), band calculation safe_windowlen big: %d\n", safe_windowlen);
    }
    /* tricky step here, W is calculated for filter using maximum of fqdb_beta_qdb and cm->beta_W, this is b/c 
     * model was (possibly) calibrated with W as calc'ed with cm->beta_W and we don't want a bigger hit
     * to be possible then was during calibration (could overestimate significance of scores), 
     * W can be less than cm->beta_W though because that would only lead to possible underestimation 
     * of significance of scores (the good direction).
     */
    if(cm->beta_W > fqdb_beta_qdb) { 
      fqdb_beta_W = cm->beta_W;
      fqdb_W  = cm->W;
    }
    else { 
      fqdb_beta_W = fqdb_beta_qdb;
      fqdb_W  = fqdb_dmax[0];
    }
    fqdb_smx = cm_CreateScanMatrix(cm, fqdb_W, fqdb_dmin, fqdb_dmax, fqdb_beta_W, fqdb_beta_qdb, TRUE, TRUE, FALSE);
    
    if(esl_opt_IsDefault(go, "--fil-S-qdb") && esl_opt_IsDefault(go, "--fil-T-qdb")) {
      
      /* No relevant options selected. Choose from, in order of priority:
       * 1. CM filter E value that gives default survival fraction 
       * 2. default CM filter bit score
       */
      if(we_have_gumbels) { /* set CM E-value cutoff as that which gives default predicted survival fraction */
	fqdb_ctype = E_CUTOFF;
	fqdb_S     = esl_opt_GetReal(go, "--fil-S-qdb");
	fqdb_S     = ESL_MAX(fqdb_S, final_S * 10.); /* final_S must be < 0.09 (enforced above), or else do_qdb_filter was set to FALSE */
	if(do_hmm_filter && fqdb_S > fhmm_S) { 
	  /* predicted survival fraction from QDB filter is higher than predicted survival fraction from the HMM filter,
	   * this means HMM should be a better filter, turn OFF QDB filter */
	  do_qdb_filter = FALSE; 
	}
	else { 
	  fqdb_E     = SurvFract2E(fqdb_S, fqdb_smx->W, cfg->avg_hit_len, cfg->dbsize);
	  fqdb_E     = ESL_MAX(fqdb_E, 1.0); /* we don't allow filter E cutoffs below 1. */
	}
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	fqdb_ctype = SCORE_CUTOFF;
	fqdb_sc    = esl_opt_GetReal(go, "--fil-T-qdb");
      }
    }
    else if(! esl_opt_IsDefault(go, "--fil-S-qdb")) {
      if(!we_have_gumbels) ESL_FAIL(eslEINVAL, errbuf, "--fil-S-qdb requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      fqdb_ctype = E_CUTOFF;
      fqdb_S     = esl_opt_GetReal(go, "--fil-S-qdb");
      fqdb_S     = ESL_MAX(fqdb_S, final_S * 10.); /* final_S must be < 0.09 (enforced above), or else do_qdb_filter was set to FALSE */
      if(do_hmm_filter) fqdb_S = ESL_MAX(fqdb_S, fhmm_S); /* predicted survival fraction from QDB filter can't be more than predicted survival fraction from the HMM filter */
      fqdb_E     = SurvFract2E(fqdb_S, fqdb_smx->W, cfg->avg_hit_len, cfg->dbsize);
      fqdb_E     = ESL_MAX(fqdb_E, 1.0); /* we don't allow filter E cutoffs below 1. */
    }
    else if(! esl_opt_IsDefault(go, "--fil-T-qdb")) {
      if(we_have_gumbels) ESL_FAIL(eslEINVAL, errbuf, "--fil-T-qdb is not allowed when <cmfile> has Gumbel statistics, use --fil-S-qdb instead.");
      fqdb_ctype = SCORE_CUTOFF;
      fqdb_sc    = esl_opt_GetReal(go, "--fil-T-qdb");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM filter cutoff selected. This shouldn't happen.");

    if(we_have_gumbels) { 
      if(fqdb_ctype == SCORE_CUTOFF) { /* determine max E-value that corresponds to fqdb_sc bit sc cutoff  across all partitions */
	if((status = Score2E(cm, errbuf, qdb_mode, fqdb_sc, &fqdb_E)) != eslOK) return status;
      }
      else if(fqdb_ctype == E_CUTOFF) { /* determine min bit sc that corresponds to fqdb_E E-val cutoff across all partitions */
	if((status  = E2Score(cm, errbuf, qdb_mode, fqdb_E, &fqdb_sc)) != eslOK) return status;
      }
      fqdb_S = E2SurvFract(fqdb_E, cm->W, cfg->avg_hit_len, cfg->dbsize);

    }
    if(fqdb_S > 0.9) do_qdb_filter = FALSE; /* turn off QDB filter if predicted survival fraction is crappy */
  }
  else do_qdb_filter = FALSE;

  /* C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first) */
  if(do_qdb_filter) { 
    AddFilterToSearchInfo(cm, TRUE, FALSE, FALSE, FALSE, FALSE, fqdb_smx, NULL, fqdb_ctype, fqdb_sc, fqdb_E);
    /* DumpSearchInfo(cm->si); */
  }
  else if (fqdb_smx != NULL) cm_FreeScanMatrix(cm, fqdb_smx); 
  /* D. add HMM filter, if necessary (after QDB filter, filters added like a stack, HMM filter is added last but used first) */
  if (do_hmm_filter) { 
    AddFilterToSearchInfo(cm, FALSE, FALSE, FALSE, TRUE, FALSE, NULL, NULL, fhmm_ctype, fhmm_sc, fhmm_E);
    /* DumpSearchInfo(cm->si); */
  }
  ValidateSearchInfo(cm, cm->si);
  return eslOK;
}

/*
 * Function: print_searchinfo
 * Date:     EPN, Thu May 17 14:47:36 2007
 * Purpose:  Print info about search (cutoffs, algorithm, etc.) to file or stdout 
 */
int print_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, FILE *fp, CM_t *cm, long N, char *errbuf)
{
  int status;
  int p, n;
  float surv_fract = 1.;
  float prv_surv_fract = 1.;
  int cutoff_type;
  float sc_cutoff;
  float e_cutoff;
  int using_filters;
  int cm_mode;
  int cp9_mode;
  int gum_mode;             /* index to use for gumbel stats in cm->stats->gumAA[], gum_mode = (stype == SEARCH_WITH_CM) cm_mode : cp9_mode; */
  int stype;
  int search_opts;
  ScanMatrix_t *smx;
  HybridScanInfo_t *hsi;
  float avg_hit_len;        /* average length of a hit, from qdb calculation */
  float hmm_ncalcs_per_res; /* millions of dp calcs per residue for HMM */
  float cm_ncalcs_per_res;  /* millions of dp calcs per residue for CM, changes per round due to QDBs */
  int   W;                  /* W for current round, can change per round due to QDBs */
  float Mc_per_res;         /* either cm_ncalcs_per_res or hmm_ncalcs_per_res */
  float seconds;            /* predicted number of seconds per round */
  int   use_qdb;            /* are we using qdb for current round? */

  /* contract check */
  if(cm->si == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo(), cm->si is NULL, shouldn't happen.\n");
  SearchInfo_t *si = cm->si;

  /* Could use ESL_GETOPTS here, but using the CM flags assures we're reporting
   * on how the CM is actually config'ed, not how we want it to be
   */
    
  using_filters = (si->nrounds > 0) ? TRUE : FALSE;

  fprintf(cfg->ofp, "#\n");
  fprintf(cfg->ofp, "# %3s  %3s  %3s  %3s  %5s  %16s  %18s\n",               ""    , "",    "",    "",    "",      "    cutoffs     ",   "   predictions    ");
  fprintf(cfg->ofp, "# %3s  %3s  %3s  %3s  %5s  %16s  %18s\n",               "",     "",    "",    "",    "",      "----------------",   "------------------");
  fprintf(cfg->ofp, "# %3s  %3s  %3s  %3s  %5s  %8s  %6s  %6s  %10s\n", "rnd",  "mod", "cfg", "alg", "beta",  "E value",  "bit sc", "surv",   "time (s)");
  fprintf(cfg->ofp, "# %3s  %3s  %3s  %3s  %5s  %8s  %6s  %6s  %10s\n", "---",  "---", "---", "---", "-----", "--------", "------", "------", "----------");

  if((status = cm_GetAvgHitLen        (cm,      errbuf, &avg_hit_len))        != eslOK) return status;
  if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &hmm_ncalcs_per_res)) != eslOK) return status;

  for(n = 0; n <= si->nrounds; n++) {
    stype       = si->stype[n];
    search_opts = si->search_opts[n];
    cutoff_type = si->cutoff_type[n];
    sc_cutoff   = si->sc_cutoff[n];
    e_cutoff    = si->e_cutoff[n];
    smx         = si->smx[n];
    hsi         = si->hsi[n];

    /* Determine configuration of CM and CP9 based on cm->flags & cm->search_opts */
    CM2Gumbel_mode(cm, search_opts, &cm_mode, &cp9_mode); 
    gum_mode = (stype == SEARCH_WITH_CM) ? cm_mode : cp9_mode; 

    use_qdb     = (smx == NULL || (smx->dmin == NULL && smx->dmax == NULL)) ? FALSE : TRUE;
    if(use_qdb) { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, smx->beta_qdb, &cm_ncalcs_per_res, &W))  != eslOK) return status; }
    else        { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, TRUE,  cm->beta_W,    &cm_ncalcs_per_res, &W))  != eslOK) return status; }

    if(cm->flags & CMH_GUMBEL_STATS) { 
      prv_surv_fract = surv_fract;
      surv_fract = E2SurvFract(e_cutoff, W, avg_hit_len, cfg->dbsize);
      Mc_per_res = (stype == SEARCH_WITH_CM) ? cm_ncalcs_per_res : hmm_ncalcs_per_res;
      seconds    = prv_surv_fract * cfg->dbsize * Mc_per_res;
      if(stype == SEARCH_WITH_CM) seconds = (search_opts & CM_SEARCH_INSIDE) ?     (seconds /  75.) : (seconds / 275.);  /*  75 Mc/S inside;  275 Mc/S CYK */
      else                        seconds = (search_opts & CM_SEARCH_HMMFORWARD) ? (seconds / 175.) : (seconds / 380.);  /* 175 Mc/S forward; 380 Mc/S viterbi */

      fprintf(cfg->ofp, "  %3d", (n+1));
      if(stype == SEARCH_WITH_CM) { 
	fprintf(cfg->ofp, "  %3s  %3s  %3s  ", "cm", ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"), ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"));
	if(use_qdb) fprintf(cfg->ofp, "%5g", smx->beta_qdb);
	else        fprintf(cfg->ofp, "%5s", "-");
      }
      else { 
	fprintf(cfg->ofp, "  %3s  %3s  %3s  %5s", "hmm", ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), "-");
      }
      if(e_cutoff < -0.1)  if((status = Score2E(cm, errbuf, gum_mode, sc_cutoff, &e_cutoff)) != eslOK) return status;
      if(e_cutoff < 0.01)  fprintf(cfg->ofp, "  %4.2e", e_cutoff);
      else                 fprintf(cfg->ofp, "  %8.3f", e_cutoff);

      fprintf(cfg->ofp, "  %6.2f  %6.4f  %10.2f\n", sc_cutoff, surv_fract, seconds);
    }
    //else cm_Fail("write code for print_searchinfo without E-values\n");
  }
  fprintf(fp, "\n");
  fflush(fp);
  return eslOK;
}

/*
 * Function: read_next_search_seq
 *
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easeled: EPN, Fri Dec  8 11:40:20 2006
 *
 * Purpose:  Given a dbfp and whether or not to take the reverse complement,
 *           reads in the next sequence and prepares reverse complement.
 *
 * Returns:  eslOK on success; eslEOF if end of file, 
 *           some other status code from esl_sqio_Read() if an error occurs.
 */
int read_next_search_seq (const ESL_ALPHABET *abc, ESL_SQFILE *dbfp, int do_revcomp, dbseq_t **ret_dbseq) 
{
  int status;
  dbseq_t *dbseq = NULL;

  ESL_ALLOC(dbseq, sizeof(dbseq_t));
  dbseq->sq[0] = NULL;
  dbseq->sq[1] = NULL;

  dbseq->sq[0] = esl_sq_CreateDigital(abc);

  while((status = esl_sqio_Read(dbfp, dbseq->sq[0])) == eslOK && (dbseq->sq[0]->n == 0)) /* skip zero length seqs */
    esl_sq_Reuse(dbseq->sq[0]);

  if(status != eslOK) goto ERROR;
  if (do_revcomp) {
    /* make a new ESL_SQ object, to store the reverse complement */
    if((dbseq->sq[1] = esl_sq_CreateDigitalFrom(abc, dbseq->sq[0]->name, dbseq->sq[0]->dsq, 
						dbseq->sq[0]->n, dbseq->sq[0]->desc, 
						dbseq->sq[0]->acc, dbseq->sq[0]->ss)) == NULL) goto ERROR;
    /* reverse complement it in place */
    revcomp(dbseq->sq[1]->abc, dbseq->sq[1], dbseq->sq[1]);
  }
  dbseq->results[0] = NULL;
  dbseq->results[1] = NULL;

  *ret_dbseq = dbseq;
  return eslOK;

 ERROR:
  if(dbseq->sq[0] != NULL) esl_sq_Destroy(dbseq->sq[0]);
  if(dbseq->sq[1] != NULL) esl_sq_Destroy(dbseq->sq[1]);
  if(dbseq != NULL) free(dbseq);
  return status;
}

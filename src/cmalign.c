/* cmalign.c
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 * SVN $Id$
 * 
 * Align sequences to a CM.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */
#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"		/* general seq analysis library   */
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_sqio.h"		
#include "esl_stopwatch.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define ALGOPTS  "--cyk,--inside,--outside,--hmmonly"        /* Exclusive choice for scoring algorithms */
#define MEMOPTS  "--small,--nosmall"                         /* Exclusive choice for memory choice */
#define ACCOPTS  "--nonbanded,--hbanded,--qdb"               /* Exclusive choice for acceleration strategies */
#define ALPHOPTS "--rna,--dna"                               /* Exclusive choice for output alphabet */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "align locally w.r.t. the model",         1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "output the alignment to file <f>, not stdout", 1 },
  { "-q",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "quiet; suppress verbose banner",         1 },
  { "-f",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "include all  match columns in output alignment", 1 },
  { "-m",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "include only match columns in output alignment", 1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use tabular output summary format, 1 line per sequence", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program",                    1 },  
#endif
  /* Miscellaneous expert options */
  { "--informat",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "specify input alignment is in format <s>, not Stockholm",  2 },
  { "--time",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print timings for alignment, band calculation, etc.", 2 },
  { "--regress", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regresion test data to file <f>", 2 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "allow informative insert emissions, do not zero them", 1 },
  /* Algorithm options */
  { "--cyk",     eslARG_NONE,"default",NULL, NULL,   ALGOPTS,      NULL,        NULL, "align with the CYK algorithm", 3 },
  { "--hmmonly", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "align to a CM Plan 9 HMM with the Viterbi algorithm",3 },
  { "--inside",  eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "don't align; return scores from the Inside algorithm", 3 },
  { "--outside", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "don't align; return scores from the Outside algorithm", 3 },
  { "--post",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "align with CYK and append posterior probabilities", 3 },
  { "--checkpost",eslARG_NONE,  FALSE, NULL, NULL,      NULL,  "--post",        NULL, "check that posteriors are correctly calc'ed", 3 },
  { "--sub",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "build sub CM for columns b/t HMM predicted start/end points", 3 },
  /* Memory options */
  { "--small",   eslARG_NONE,"default",  NULL, NULL,  MEMOPTS,      NULL, "--nosmall", "use divide and conquer (d&c) alignment algorithm", 4 },
  { "--nosmall", eslARG_NONE,   FALSE,  NULL, NULL,   MEMOPTS,      NULL,   "--small", "use normal alignment algorithm, not d&c", 4 },
  /* Banded alignment */
  { "--nonbanded",eslARG_NONE,"default",NULL, NULL,  ACCOPTS,      NULL,        NULL, "do not use bands to accelerate aln algorithm", 5 },
  { "--hbanded", eslARG_NONE,   FALSE,  NULL, NULL,  ACCOPTS,"--nosmall",       NULL, "accelerate using CM plan 9 HMM banded CYK aln algorithm", 5 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 5 },
  { "--hsafe",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "realign (non-banded) seqs with HMM banded CYK score < 0 bits", 5 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 5 },
  { "--qdb",     eslARG_NONE,   FALSE, NULL, NULL,   ACCOPTS,      NULL,        NULL, "use query dependent banded CYK alignment algorithm", 5 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,   "--qdb",        NULL, "set tail loss prob for --qdb to <x>", 5 },
  /* Including a preset alignment */
  { "--withali", eslARG_STRING, NULL,  NULL, NULL,      NULL,    "--cyk",       NULL, "incl. alignment in <f> (must be aln <cm file> was built from)", 6 },
  { "--rf",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--withali",       NULL, "--rf was originally used with cmbuild", 6 },
  { "--gapthresh",eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,"--withali",       NULL, "--gapthresh <x> was originally used with cmbuild", 6 },
  /* Enforcing a subsequence */
  { "--enfstart",eslARG_INT,    FALSE, NULL, "n>0",     NULL,"--enfseq",        NULL, "enforce MATL stretch starting at consensus position <n>", 7 },
  { "--enfseq",  eslARG_STRING, NULL,  NULL, NULL,      NULL,"--enfstart",      NULL, "enforce MATL stretch starting at --enfstart <n> emits seq <s>", 7 },
  /* Verbose output files/debugging */
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump individual sequence tracebacks to file <f>", 8 },
  { "--banddump",eslARG_INT,    "0",   NULL, "0<=n<=3", NULL,      NULL,        NULL, "set verbosity of band info print statements to <n>", 8 },
  { "--dlev",    eslARG_INT,    "0",   NULL, "0<=n<=3", NULL,      NULL,        NULL, "set verbosity of debugging print statements to <n>", 8 },
/* Setting output alphabet */
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as RNA sequence data", 9},
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 9},
/* Other options */
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "arrest after start: for debugging MPI under gdb",   10 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char         *cmfile;	        /* name of input CM file  */ 
  char         *sqfile;	        /* name of sequence file  */ 
  ESL_SQFILE   *sqfp;           /* open sequence input file stream */
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  int           ncm;            /* number cm we're on */
  int           be_verbose;	

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *ofp;		/* output file (default is stdout) */
  FILE         *tracefp;	/* optional output for parsetrees  */
  FILE         *regressfp;	/* optional output for regression test  */
  ESL_MSAFILE  *withalifp;	/* optional input alignment to include */
  ESL_MSA      *withmsa;	/* MSA from withalifp to include */
  Parsetree_t  *withali_mtr;	/* guide tree for MSA from withalifp */
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

static int process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);

static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int check_withali(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr);
static int include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, int *ret_nseq, char *errbuf);
static int compare_cm_guide_trees(CM_t *cm1, CM_t *cm2);
static int make_aligned_string(char *aseq, char *gapstring, int alen, char *ss, char **ret_s);
static int determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker);
static int add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  struct cfg_s     cfg;

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
      puts("\nexpert miscellaneous options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions specifying alignment algorithm:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\nmemory related options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nbanded dynamic programming acceleration options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\noptions for including a fixed alignment within output alignment:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\noptions for enforcing alignment of a single-stranded subsequence:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\nverbose output files and debugging:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\n  other options:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
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
  cfg.sqfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.fmt        = eslSQFILE_UNKNOWN;      /* autodetect sequence file format by default. */ 
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if (esl_opt_GetBoolean(go, "-1")) cfg.be_verbose = FALSE;        
  else                              cfg.be_verbose = TRUE;        
  if      (esl_opt_GetBoolean(go, "--rna")) cfg.abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna")) cfg.abc_out = esl_alphabet_Create(eslDNA);
  else    esl_fatal("Can't determine output alphabet");
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ofp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tracefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.regressfp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withalifp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withmsa    = NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.withali_mtr= NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.ncm        = 0;

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
      cfg.be_verbose = FALSE;

      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("ERROR, MPI mode, but only 1 processor running...");

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
  if (cfg.my_rank == 0) esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (! esl_opt_IsDefault(go, "-o")) { fclose(cfg.ofp); }
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.tracefp   != NULL) fclose(cfg.tracefp);
    if (cfg.regressfp != NULL) fclose(cfg.regressfp);
    if (cfg.withalifp != NULL) esl_msafile_Close(cfg.withalifp);
    if (cfg.withmsa   != NULL) esl_msa_Destroy(cfg.withmsa);
    if (cfg.withali_mtr != NULL) FreeParsetree(cfg.withali_mtr);
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
 * Sets: 
 *    cfg->sqfp        - open sequence file                
 *    cfg->ofp         - output file (stdout by default)
 *    cfg->cmfp        - open CM file                
 *    cfg->abc         - digital input alphabet
 *    cfg->tracefp     - optional output file
 *    cfg->regressfp   - optional output file
 *    cfg->withalifp   - optional input alignment file to include
 *    cfg->withmsa     - MSA from --withali file 
 *    cfg->withali_mtr - guide tree for MSA from --withali file 
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

  /* open input sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(cfg->sqfp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", cfg->sqfile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", cfg->sqfile);
  else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
  cfg->fmt = cfg->sqfp->format;

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, NULL, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
    } else cfg->ofp = stdout;

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tracefp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open regression file */
  if (esl_opt_GetString(go, "--regress") != NULL) {
    if ((cfg->regressfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
    }

  /* optionally, open withali file for reading */
  if(esl_opt_GetString(go, "--withali") != NULL)
    {
      status = esl_msafile_OpenDigital(cfg->abc, esl_opt_GetString(go, "--withali"), eslMSAFILE_UNKNOWN, NULL, &(cfg->withalifp));
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --withali alignment %s\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
    }

  return eslOK;
}

/* serial_master()
 * The serial version of cmalign.
 * Align each sequence to the CM.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[eslERRBUFSIZE];
  CM_t     *cm;
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, holds seqs, parsetrees, CP9 traces, postcodes */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/
  
  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;
      
      /* initialize the flags/options/params and configuration of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK)    cm_Fail(errbuf);
      /* read in all sequences, this is wasteful, but Parsetrees2Alignment() requires all seqs in memory */
      seqs_to_aln = CreateSeqsToAln(100, FALSE);
      if((status = ReadSeqsToAln(cfg->abc, cfg->sqfp, 0, TRUE, seqs_to_aln, FALSE)) != eslEOF) cm_Fail("Error reading sqfile: %s\n", cfg->sqfile);
      /* align all sequences */
      if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
      if ((status = output_result   (go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
      
      /* clean up */
      FreeSeqsToAln(seqs_to_aln);
      FreeCM(cm);
    }   
  return;
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmalign
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * EPN: GOAL OF IMPLEMENTATION FOLLOWS IN LOWERCASE.
 * IT IS NOT YET ACHIEVED.
 * TO ACHIEVE WE'LL NEED ALL FUNCS CALLED BY MPI TO
 * RETURN CLEANLY ALWAYS - BIG TASK TO REWRITE THOSE.
 * CURRENTLY NEARLY ALL ERRORS ARE UNRECOVERABLE, BUT THESE
 * ARE NOT LIMITED TO MPI COMMUNICATION ERRORS.
 *
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include all worker-side errors, and any
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

  CM_t *cm;
  int nseq_per_worker;
  int nseq_this_worker;
  int nseq_prev;

  seqs_to_aln_t  *all_seqs_to_aln    = NULL;
  seqs_to_aln_t  *worker_seqs_to_aln = NULL;
  int            *seqidx         = NULL;
  
  char     errbuf[eslERRBUFSIZE];
  MPI_Status mpistatus; 
  int      n;

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  /*if (xstatus == eslOK) { if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }*/
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqidx  = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  for (wi = 0; wi < cfg->nproc; wi++) seqidx[wi] = 0;
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) cm_Fail(errbuf);
  ESL_DPRINTF1(("MPI master is initialized\n"));

  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
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

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      cfg->ncm++;  
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      determine_nseq_per_worker(go, cfg, cm, &nseq_per_worker); /* this func dies internally if there's some error */
      printf("nseq_per_worker: %d\n", nseq_per_worker);

      wi = 1;
      all_seqs_to_aln = CreateSeqsToAln(100, TRUE);
      while (have_work || nproc_working)
	{
	  if (have_work) 
	    {
	      nseq_prev = all_seqs_to_aln->nseq;
	      if((status = ReadSeqsToAln(cfg->abc, cfg->sqfp, nseq_per_worker, FALSE, all_seqs_to_aln, TRUE)) == eslOK)
		{
		  nseq_this_worker = all_seqs_to_aln->nseq - nseq_prev;
		  ESL_DPRINTF1(("MPI master read %d seqs\n", all_seqs_to_aln->nseq));
		}
	      else 
		{
		  have_work = FALSE;
		  if (status != eslEOF) cm_Fail("Sequence file read unexpectedly failed with code %d\n", status); 
		  ESL_DPRINTF1(("MPI master has run out of sequences to read (having read %d)\n", 0));
		} 
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
		  ESL_DPRINTF1(("MPI master sees that the result buffer contains aligned sequences (seqidx: %d)\n", seqidx[wi]));
		  if ((status = cm_seqs_to_aln_MPIUnpack(cfg->abc, buf, bn, &pos, MPI_COMM_WORLD, &worker_seqs_to_aln)) != eslOK) cm_Fail("search results unpack failed");
		  ESL_DPRINTF1(("MPI master has unpacked search results\n"));
		  if ((status = add_worker_seqs_to_master(all_seqs_to_aln, worker_seqs_to_aln, seqidx[wi])) != eslOK) cm_Fail("adding worker results to master results failed");
		  /* careful not to free data from worker_seqs_to_aln we've
		   * just added to all_seqs_to_aln. we didn't copy it, we just
		   * had pointers in all_seqs_to_aln point to it. We can 
		   * free the worker's pointers to those pointers though */
		  if(worker_seqs_to_aln->sq       != NULL) free(worker_seqs_to_aln->sq);
		  if(worker_seqs_to_aln->tr       != NULL) free(worker_seqs_to_aln->tr);
		  if(worker_seqs_to_aln->cp9_tr   != NULL) free(worker_seqs_to_aln->cp9_tr);
		  if(worker_seqs_to_aln->postcode != NULL) free(worker_seqs_to_aln->postcode);
		  if(worker_seqs_to_aln->sc       != NULL) free(worker_seqs_to_aln->sc);
		  free(worker_seqs_to_aln);
								
		}
	      else	/* worker reported an error. Get the errbuf. */
		{
		  if (MPI_Unpack(buf, bn, &pos, errbuf, eslERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack of errbuf failed");
		  ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
		}
	      nproc_working--;
	    }
	  
	  if (have_work)
	    {   
	      /* send new alignment job */
	      ESL_DPRINTF1(("MPI master is sending sequence to search to worker %d\n", wi));
	      if ((status = cm_seqs_to_aln_MPISend(all_seqs_to_aln, all_seqs_to_aln->nseq - nseq_this_worker, nseq_this_worker, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
	      seqidx[wi] = all_seqs_to_aln->nseq - nseq_this_worker;
	      wi++;
	      nproc_working++;
	    }
	}
      /* we have all worker's results, output the alignment */
      if ((status = output_result(go, cfg, errbuf, cm, all_seqs_to_aln)) != eslOK) cm_Fail(errbuf);
      ESL_DPRINTF1(("MPI master: done with this CM. Telling all workers\n"));
      /* send workers the message that we're done with this stage */
      for (wi = 1; wi < cfg->nproc; wi++) 
	if ((status = cm_seqs_to_aln_MPISend(NULL, 0, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");
      FreeCM(cm);
      esl_sqio_Rewind(cfg->sqfp); /* we may be searching this file again with another CM */
    }
  
  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((status = cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if (xstatus != eslOK) cm_Fail(errbuf);
  else                  return;

}



static void
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int           xstatus = eslOK;
  int           status;
  CM_t         *cm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
  int           pos;
  char          errbuf[eslERRBUFSIZE];
  seqs_to_aln_t *seqs_to_aln = NULL;
  int           i;
  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
  
  /* Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized\n", cfg->my_rank));
  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      
      while((status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln)) == eslOK)
	{
	  ESL_DPRINTF1(("worker %d: has received alignment job, nseq: %d\n", cfg->my_rank, seqs_to_aln->nseq));
	  /* align all sequences */
	  if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  ESL_DPRINTF1(("worker %d: has gathered alignment results\n", cfg->my_rank));

	  /* free the sequences, master already has a copy so we don't send them back */
	  for(i = 0; i < seqs_to_aln->nseq; i++) esl_sq_Destroy(seqs_to_aln->sq[i]);
	  free(seqs_to_aln->sq);
	  seqs_to_aln->sq = NULL;
	  n = 0;
	  if (cm_seqs_to_aln_MPIPackSize(seqs_to_aln, 0, seqs_to_aln->nseq, MPI_COMM_WORLD, &sz) != eslOK) goto ERROR; n += sz;  
	  if (n > wn) {
	    void *tmp;
	    ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	    wn = n;
	  }
	  ESL_DPRINTF1(("worker %d: has calculated the alignment results will pack into %d bytes\n", cfg->my_rank, n));
	  status = eslOK;

	  pos = 0;
	  if (cm_seqs_to_aln_MPIPack(seqs_to_aln, 0, seqs_to_aln->nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) goto ERROR;
	  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	  ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));

	  FreeSeqsToAln(seqs_to_aln);
	}
      if(status == eslEOD)ESL_DPRINTF1(("worker %d: has seen message to stop with this CM.\n", cfg->my_rank));
      else goto ERROR;
      FreeCM(cm);
      cm = NULL;
    }
  if (status == eslEOD) ESL_DPRINTF1(("worker %d told CMs are done.\n", cfg->my_rank));
  else goto ERROR;
  
  if (wbuf != NULL) free(wbuf);
  return;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  eslERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
  return;
}
#endif /*HAVE_MPI*/

static int
output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln)
{
  int status;
  ESL_MSA *msa = NULL;
  int i, imax;

  /* create a new MSA, if we didn't do either --inside or --outside */
  if(esl_opt_GetBoolean(go, "--cyk") || esl_opt_GetBoolean(go, "--hmmonly"))
    {
      /* optionally include a fixed alignment provided with --withali,
       * this has already been checked to see it matches the CM structure */
      if(esl_opt_GetString(go, "--withali") != NULL)
	{
	  if((status = include_withali(go, cfg, cm, &(seqs_to_aln->sq), &(seqs_to_aln->tr), &(seqs_to_aln->nseq), errbuf)) != eslOK)
	    ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have a SS_cons compatible with the CM\n", esl_opt_GetString(go, "--withali"));
	}
      
      if(esl_opt_GetBoolean(go, "--hmmonly"))
	{
	  assert(seqs_to_aln->cp9_tr != NULL);
	  if((status = CP9Traces2Alignment(cm, cfg->abc_out, seqs_to_aln->sq, NULL, seqs_to_aln->nseq, seqs_to_aln->cp9_tr, 
					   esl_opt_GetBoolean(go, "-f"), esl_opt_GetBoolean(go, "-m"), &msa)) != eslOK)
	    goto ERROR;
	}
      else
	{
	  assert(seqs_to_aln->tr != NULL);
	  if((status = Parsetrees2Alignment(cm, cfg->abc_out, seqs_to_aln->sq, NULL, seqs_to_aln->tr, seqs_to_aln->nseq, 
					    esl_opt_GetBoolean(go, "-f"), esl_opt_GetBoolean(go, "-m"), &msa)) != eslOK)
	    goto ERROR;
	}
      if(esl_opt_GetBoolean(go, "--post")) 
	{                                                                              
	  char *apostcode;   /* aligned posterior decode array */
	  if(seqs_to_aln->postcode == NULL) 
	    cm_Fail("ERROR --post enabled, but {serial,parallel}_align_targets() did not return post codes.\n");

	  imax = seqs_to_aln->nseq - 1;
	  if(cfg->withmsa != NULL) imax -= cfg->withmsa->nseq;
	  for (i = 0; i <= imax; i++)                                                   
	    {                                                                          
	      if((status =make_aligned_string(msa->aseq[i], "-_.", msa->alen, seqs_to_aln->postcode[i], &apostcode)) != eslOK)
		ESL_FAIL(status, errbuf, "error creating posterior string\n");
	      /*esl_msa_AppendGR(msa, "POST", i, apostcode);*/                                  
	      free(apostcode);                                                         
	    }                                                                          
	}                                                                              
      status = esl_msa_Write(cfg->ofp, msa, eslMSAFILE_STOCKHOLM);
      if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
      else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);
      
      /* Detailed traces for debugging training set. */
      if(cfg->tracefp != NULL)
	{
	  if(esl_opt_GetBoolean(go,"--hmmonly")) { printf("%-40s ... ", "Saving CP9 HMM traces"); fflush(stdout); }
	  else                                   { printf("%-40s ... ", "Saving CM parsetrees");  fflush(stdout); }
	  for (i = 0; i < msa->nseq; i++) 
	    {
	      fprintf(cfg->tracefp, "> %s\n", seqs_to_aln->sq[i]->name);
	      if(esl_opt_GetBoolean(go,"--hmmonly")) 
		{
		  fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", CP9TraceScore(cm->cp9, seqs_to_aln->sq[i]->dsq, seqs_to_aln->cp9_tr[i]));
		  CP9PrintTrace(cfg->tracefp, seqs_to_aln->cp9_tr[i], cm->cp9, seqs_to_aln->sq[i]->dsq);
		}
	      else
		{
		  fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE));
		  ParsetreeDump(cfg->tracefp, seqs_to_aln->tr[i], cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); /* NULLs are dmin, dmax */
		}
	      fprintf(cfg->tracefp, "//\n");
	    }
	  printf("done. [%s]\n", esl_opt_GetString(go, "--tfile"));
	}
      if (cfg->regressfp != NULL)
	{
	  /* Must delete author info from msa, because it contains version
	   * and won't diff clean in regression tests. */
	  if(msa->au != NULL) free(msa->au); msa->au = NULL;
	  status = esl_msa_Write(cfg->ofp, msa, eslMSAFILE_STOCKHOLM);
	  if (status == eslEMEM)    ESL_FAIL(status, errbuf, "Memory error when outputting regression file\n");
	  else if (status != eslOK) ESL_FAIL(status, errbuf, "Writing regression file failed with error %d\n", status);
	}
    }
  if(msa != NULL) esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if(msa != NULL) esl_msa_Destroy(msa);
  return status;
}

/* An alignment work unit consists a seqs_to_aln_t object which contains sequences to align, 
 * and space for their parsetrees, or CP9 traces, and possibly postal codes.
 * The job is to align the sequences and create parsetrees or cp9 traces and maybe postal codes.
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, 
		 seqs_to_aln_t *seqs_to_aln)
{
  actually_align_targets(cm, seqs_to_aln,
			 NULL, NULL,   /* we're not aligning search hits */
			 esl_opt_GetInteger(go, "--banddump"),
			 esl_opt_GetInteger(go, "--dlev"), esl_opt_GetBoolean(go, "-q"), NULL);
  return eslOK;
  
  /* ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_search_workunit\n", cfg->my_rank));
  FreeCM(cm);
  FreeResults(results);
  return status;*/
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults.
 * Configures the CM with a ConfigCM() call at end.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;

  /* set up params/flags/options of the CM */
  cm->beta   = esl_opt_GetReal(go, "--beta"); /* this will be DEFAULT_BETA unless changed at command line */
  cm->tau    = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  /* update cm->config->opts */
  if(esl_opt_GetBoolean(go, "-l"))
    {
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cm->config_opts |= CM_CONFIG_HMMEL;
    }
  if(! esl_opt_GetBoolean(go, "--iins"))      cm->config_opts |= CM_CONFIG_ZEROINSERTS;

  /* update cm->align->opts */
  if(esl_opt_GetBoolean(go, "--hbanded"))     cm->align_opts  |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--sums"))        cm->align_opts  |= CM_ALIGN_SUMS;
  if(esl_opt_GetBoolean(go, "--sub"))         cm->align_opts  |= CM_ALIGN_SUB;
  if(esl_opt_GetBoolean(go, "--hmmonly"))     cm->align_opts  |= CM_ALIGN_HMMONLY;
  if(esl_opt_GetBoolean(go, "--inside"))      cm->align_opts  |= CM_ALIGN_INSIDE;
  if(esl_opt_GetBoolean(go, "--outside"))     cm->align_opts  |= CM_ALIGN_OUTSIDE;
  if(esl_opt_GetBoolean(go, "--nosmall"))     cm->align_opts  |= CM_ALIGN_NOSMALL;
  if(esl_opt_GetBoolean(go, "--post"))        cm->align_opts  |= CM_ALIGN_POST;
  if(esl_opt_GetBoolean(go, "--time"))        cm->align_opts  |= CM_ALIGN_TIME;
  if(esl_opt_GetBoolean(go, "--checkpost"))   cm->align_opts  |= CM_ALIGN_CHECKINOUT;
  if(esl_opt_GetBoolean(go, "--hsafe"))       cm->align_opts  |= CM_ALIGN_HMMSAFE;
  if(esl_opt_GetString (go, "--enfseq") != NULL)
    {
      cm->config_opts |= CM_CONFIG_ENFORCE;
      cm->enf_start    = EnforceFindEnfStart(cm, esl_opt_GetInteger(go, "--enfstart"));
      cm->enf_seq      = esl_opt_GetString(go, "--enfseq");
    }

  /* config QDB? */
  if(esl_opt_GetBoolean(go, "--qdb"))          
    { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }

  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  ConfigCM(cm, NULL, NULL); 
  if(cm->config_opts & CM_CONFIG_ENFORCE) ConfigCMEnforce(cm);

  if(cfg->my_rank == 0) printf("CM %d: %s\n", (cfg->ncm), cm->name);

  /* if we're master and we're trying to include an alignment, make sure it is consistent with CM structure */
  if(cfg->my_rank == 0 && (! esl_opt_IsDefault(go, "--withali"))) { 
    if((status = check_withali(go, cfg, cm, &(cfg->withmsa), &(cfg->withali_mtr))) != eslOK)
      ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have a SS_cons compatible with the CM\n", esl_opt_GetString(go, "--withali"));
    }
  return eslOK;
}

/* Function: compare_cm_guide_trees()
 * EPN, Tue Mar  6 08:32:12 2007
 *
 * Purpose:  Given two CMs, cm1 and cm2, compare them, returning TRUE 
 *           iff they have the same guide tree (same node architecture).
 *
 * Args:     cm1          - covariance model number 1
 *           cm2          - covariance model number 2
 * 
 * Returns:  TRUE if CMs have same guide tree, FALSE otherwise
 */
static int compare_cm_guide_trees(CM_t *cm1, CM_t *cm2)
{
  int          nd; 
  if(cm1->nodes != cm2->nodes) return FALSE;
  for(nd = 0; nd < cm1->nodes; nd++)
    if(cm1->ndtype[nd] != cm2->ndtype[nd]) return FALSE;
  return TRUE;
}

/* Function: check_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Ensure that the alignment to include has a secondary
 *           structure that matches our CM. Pass the alignment back
 *           as *ret_msa.
 *
 * Returns:  <eslOK> on success.
 *           <eslEINCOMPAT> if alignment doesn't match the CM 
 */
static int check_withali(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr)
{
  int           status;
  ESL_MSA      *msa      = NULL; /* alignment we're including  */
  CM_t         *new_cm   = NULL; /* CM built from MSA, we check it has same guide tree as 'cm' */
  CM_t         *newer_cm = NULL; /* used briefly if we rebalance new_cm */
  Parsetree_t  *mtr      = NULL; /* master structure tree from the alignment*/
  char          errbuf[eslERRBUFSIZE];

  /* cfg->withalifp is open */
  status = esl_msa_Read(cfg->withalifp, &msa);
  if (status == eslEFORMAT)  cm_Fail("--withali alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->withalifp->linenumber, cfg->withalifp->fname, cfg->withalifp->errbuf, cfg->withalifp->buf);
  else if (status != eslOK)       cm_Fail("--withali alignment file read unexpectedly failed with code %d\n", status);

  /* Some input data cleaning. */
  if (esl_opt_GetBoolean(go, "--rf") && msa->rf == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--rf invoked but --withali alignment has no reference coord annotation.\n");
  if (msa->ss_cons == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--withali alignment did not contain consensus structure annotation.\n");
  if (! clean_cs(msa->ss_cons, msa->alen))
    ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation for --withali alignment\n");

  /* Build a CM from a master guide tree built from the msa, 
   * then check to make sure this CM has same emit map as the CM
   * we've had passed in. This is fragile and hopefully temporary. 
   * Another solution would be to use a checksum, but CM files don't 
   * have checksums yet.
   */
  HandModelmaker(msa, esl_opt_GetBoolean(go, "--rf"), esl_opt_GetReal(go, "--gapthresh"), &new_cm, &mtr);
  if(!(compare_cm_guide_trees(cm, new_cm)))
    {
      CM_t *newer_cm;
      newer_cm = CMRebalance(new_cm);
      FreeCM(new_cm);
      new_cm = NULL;
      if(!(compare_cm_guide_trees(cm, newer_cm)))
	{
	  status = eslEINCOMPAT;
	  goto ERROR;
	}
    }

  /* if we get here, the CM guide trees match */
  if(new_cm   != NULL) FreeCM(new_cm);
  if(newer_cm != NULL) FreeCM(newer_cm);
  *ret_mtr = mtr;
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if(msa != NULL)      esl_msa_Destroy(msa);
  if(new_cm   != NULL) FreeCM(new_cm);
  if(newer_cm != NULL) FreeCM(newer_cm);
  if(mtr != NULL)      FreeParsetree(mtr);
  return eslEINCOMPAT;
}

/* Function: include_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Infer the implicit parses of sequences in an
 *           MSA to a CM and append to passed in data structures
 *           We've already checked to make sure the MSA's consensus
 *           structure matches the CM.
 *
 * Args:     go           - command line options
 *           cfg          - cmalign configuration, includes msa to add
 *           cm           - CM we're aligning to 
 *           ret_sq       - pre-existing sequences, to append msa seqs to
 *           ret_tr       - pre-existing parsetrees, to append to
 *           ret_nseq     - number of exisiting seqs, updated at end of function
 *           errbuf       - easel error message
 * 
 * Returns:  eslOK on success, eslEMEM on memory error
 *           Also new, realloc'ed arrays for sq, tr in ret_seq, ret_tr; 
 *           *ret_nseq is increased by number of seqs in cfg->withmsa.
*/
static int include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, int *ret_nseq, char *errbuf)
{
  int           status;
  void         *tmp;      /* for ESL_RALLOC() */
  int           i;	  /* counter over aseqs       */
  int           ip;	  /* offset counter over aseqs */
  char        **uaseq;    /* unaligned seqs, dealigned from the MSA */
  char         **aseq;    /*   aligned text seqs */
  int           apos;     /*   aligned position index */
  int           uapos;    /* unaligned position index */
  int           x;        /* counter of parsetree nodes */
  int         **map;      /* [0..msa->nseq-1][0..msa->alen] map from aligned
			   * positions to unaligned (non-gap) positions */

  /* Contract check */
  if(cfg->withmsa == NULL) esl_fatal("ERROR in include_withali() withmsa is NULL.\n");
  if(! (cfg->withmsa->flags & eslMSA_DIGITAL)) esl_fatal("ERROR in include_withali() withmsa is not digitized.\n");

  /* For each seq in the MSA, map the aligned sequences coords to 
   * the unaligned coords, we stay in digitized seq coords (1..alen),
   * we need this for converting parsetrees from Transmogrify (which
   * have emitl and emitr in aligned coords) to unaligned coords, so 
   * we can call Parsetrees2Alignment() with them. */
  ESL_ALLOC(map,   sizeof(int *)  * cfg->withmsa->nseq);
  ESL_ALLOC(uaseq, sizeof(char *) * cfg->withmsa->nseq);
  ESL_ALLOC(aseq,  sizeof(char *) * cfg->withmsa->nseq);
  for (i = 0; i < cfg->withmsa->nseq; i++)
    {
      ESL_ALLOC(map[i],   sizeof(int)  * (cfg->withmsa->alen+1));
      ESL_ALLOC(aseq[i],  sizeof(char) * (cfg->withmsa->alen+1));
      map[i][0] = -1; /* invalid */
      uapos = 1;
      for(apos = 1; apos <= cfg->withmsa->alen; apos++)
	{
	  if (!esl_abc_XIsGap(cfg->withmsa->abc, cfg->withmsa->ax[i][apos]))
	    map[i][apos] = uapos++;
	  else
	    map[i][apos] = -1;
	}
      /* we need digitized AND text seqs for Transmogrify */
      esl_abc_Textize(cfg->withmsa->abc, cfg->withmsa->ax[i], cfg->withmsa->alen, aseq[i]);
      esl_strdup(aseq[i], -1, &(uaseq[i]));
      esl_sq_Dealign(uaseq[i], uaseq[i], "-_.", cfg->withmsa->alen);
    }
  ESL_RALLOC((*ret_tr),  tmp, (sizeof(Parsetree_t *)  * (*ret_nseq + cfg->withmsa->nseq)));
  ESL_RALLOC((*ret_sq),  tmp, (sizeof(ESL_SQ *)       * (*ret_nseq + cfg->withmsa->nseq)));

  /* Transmogrify each aligned seq to get a parsetree */
  /*for (i = 0; i < cfg->withmsa->nseq; i++)*/
  for (i = *ret_nseq; i < (*ret_nseq + cfg->withmsa->nseq); i++)
    {
      ip = i - *ret_nseq;
      (*ret_tr)[i] = Transmogrify(cm, cfg->withali_mtr, cfg->withmsa->ax[ip], aseq[ip], cfg->withmsa->alen);
      /* ret_tr[i] is in alignment coords, convert it to unaligned coords, */
      for(x = 0; x < (*ret_tr)[i]->n; x++)
	{
	  if((*ret_tr)[i]->emitl[x] != -1)
	    (*ret_tr)[i]->emitl[x] = map[ip][(*ret_tr)[i]->emitl[x]];
	  if((*ret_tr)[i]->emitr[x] != -1)
	    (*ret_tr)[i]->emitr[x] = map[ip][(*ret_tr)[i]->emitr[x]];
	}
      (*ret_sq)[i]      = esl_sq_CreateFrom(cfg->withmsa->sqname[ip], uaseq[ip], NULL, NULL, NULL);
      esl_sq_Digitize(cm->abc, (*ret_sq)[i]);
    }
  *ret_nseq += cfg->withmsa->nseq; /* we added these seqs to sq, tr */

  /* Clean up and exit. */
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  esl_Free2D((void **) aseq,  cfg->withmsa->nseq);
  return eslOK;

 ERROR:
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  return status;
}


/* Function: make_aligned_string() 
 * Incept:   EPN, Thu Aug  2 18:24:49 2007
 *           stolen from Squid:alignio.c:MakeAlignedString()
 * 
 * Purpose:  Given a raw string of some type (secondary structure, say),
 *           align it to a given aseq by putting gaps wherever the
 *           aseq has gaps. 
 *           
 * Args:     aseq:  template for alignment
 *           gapstring: defines all gap chars ex: "-_."
 *           alen:  length of aseq
 *           ss:    raw string to align to aseq
 *           ret_s: RETURN: aligned ss
 *           
 * Return:   eslOK on success, eslEMEM on memory allocation error,
 *           eslEINCONCEIVABLE on strange error,
 *           ret_ss is alloc'ed here and must be free'd by caller.
 */
int
make_aligned_string(char *aseq, char *gapstring, int alen, char *ss, char **ret_s)
{
  int status;
  char *new; 
  int   apos, rpos;

  ESL_ALLOC(new, (sizeof(char) * (alen+1)));
  for (apos = rpos = 0; apos < alen; apos++)
    {
      if (strchr(gapstring, aseq[apos]) != NULL)
	new[apos] = aseq[apos];
      else
	new[apos] = ss[rpos++];
    }
  new[apos] = '\0';

  if (rpos != strlen(ss))
    {
      if(new != NULL) free(new);
      return eslEINCONCEIVABLE;
    }
  *ret_s = new;
  return eslOK;

 ERROR:
  if(new != NULL) free(new);
  return status;
}

/* determine_nseq_per_worker()
 * Given a CM, return the number of sequences we think we should send
 * to each worker (we don't know the number of sequences in the file).
 * The calculation is based on trying to get a worker to spend 
 * a specific amount of time MPI_WORKER_ALIGN_TARGET_SEC, a constant
 * from structs.h. 
 */
static int
determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker)
{
  *ret_nseq_worker = 5;
  return eslOK;
}

/* add_worker_seqs_to_master
 * Add results (parstrees or CP9 traces, and possibly postcodes) from a
 * worker's seqs_to_aln object to a master seqs_to_aln object.
 */
static int
add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset)
{
  int x;

  if(worker_seqs->sq != NULL) cm_Fail("add_worker_seqs_to_master(), worker_seqs->sq non-NULL.");
  if(master_seqs->nseq < (offset + worker_seqs->nseq)) cm_Fail("add_worker_seqs_to_master(), master->nseq: %d, offset %d, worker->nseq: %d\n", master_seqs->nseq, offset, worker_seqs->nseq);

  if(worker_seqs->tr != NULL) {
    if(master_seqs->tr == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned parsetrees, master->tr is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->tr[x] == NULL); 
      master_seqs->tr[x] = worker_seqs->tr[(x-offset)];
    }
  }

  if(worker_seqs->cp9_tr != NULL) {
    if(master_seqs->cp9_tr == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned cp9 traces, master->cp9_tr is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->cp9_tr[x] == NULL); 
      master_seqs->cp9_tr[x] = worker_seqs->cp9_tr[(x-offset)];
    }
  }

  if(worker_seqs->postcode != NULL) {
    if(master_seqs->postcode == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned postcodes, master->postcode is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->postcode[x] == NULL); 
      master_seqs->postcode[x] = worker_seqs->postcode[(x-offset)];
    }
  }

  if(worker_seqs->sc != NULL) {
    if(master_seqs->sc == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned scores, master->sc is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(!(NOT_IMPOSSIBLE(master_seqs->sc[x])));
      master_seqs->sc[x] = worker_seqs->sc[(x-offset)];
    }
  }

  return eslOK;
}

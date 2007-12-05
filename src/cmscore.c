/* SRE, Thu Aug  3 17:08:44 2000 [StL]
 * SVN $Id$
 * 
 * Score a CM against unaligned sequence examples.
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
#include <math.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_sqio.h"		
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define ALGOPTS  "--std,--qdb,--qdbsmall,--qdbboth,--hbanded,--hmmviterbi"  /* Exclusive choice for scoring algorithms */
#define SEQOPTS  "--emit,--random,--infile"                                 /* Exclusive choice for sequence input */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,     "--sub", "align locally w.r.t. the model",         1 },
  { "-i",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print individual timings & scores, not just summary", 1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use tabular output summary format, 1 line per sequence", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program",                    1 },  
#endif
  /* Miscellaneous expert options */
  { "--stringent",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "require the two parsetrees to be identical", 2 },
  { "--sub",      eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        "-l", "build sub CM for columns b/t HMM predicted start/end points", 2 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "allow informative insert emissions, do not zero them", 2 },
  /* options for generating/reading sequences to score */
  { "--emit",    eslARG_INT,    "100", NULL, "n>0",     NULL,      NULL,     SEQOPTS, "emit <n> sequences from each CM [default: 100]", 3 },
  { "--random",  eslARG_INT,    "100", NULL, "n>0",     NULL,      NULL,     SEQOPTS, "emit <n> random seq from cm->null model (length = CM consensus)", 3},
  { "--infile",  eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL,     SEQOPTS, "read sequences to align from file <s>", 3 },
  { "--outfile", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,  "--infile", "save seqs to file <s>", 3 },
  { "--seed",    eslARG_INT,     NULL, NULL, "n>0",     NULL,      NULL,    "--seed", "set random number seed to <n>", 3 },
  /* Stage 2 algorithm options */
  { "--std",     eslARG_NONE,"default",NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare divide and conquer versus standard CYK", 4 },
  { "--qdb",     eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare non-banded d&c versus QDB standard CYK", 4 },
  { "--qdbsmall",eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare non-banded d&c versus QDB d&c", 4 },
  { "--qdbboth", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare        QDB d&c versus QDB standard CYK", 4 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      NULL,        NULL, "set tail loss prob for QDB to <x>", 4 },
  { "--hbanded", eslARG_NONE,   FALSE,  NULL, NULL,  ALGOPTS,      NULL,        NULL, "accelerate using CM plan 9 HMM banded CYK aln algorithm", 4 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 4 },
  { "--hsafe",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "realign (non-banded) seqs with HMM banded CYK score < 0 bits", 4 },
  { "--hmmviterbi",eslARG_NONE, FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "align to a CM Plan 9 HMM with the Viterbi algorithm", 4 },
  { "--scoreonly",eslARG_NONE,  FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "for standard CYK stage, do only score, save memory", 4 },
  /* Options for testing multiple rounds of banded alignment, stage 2->N alignment */
  { "--betas",   eslARG_INT,  NULL,    NULL, "0<n<50",   NULL, "--betae",       NULL, "set initial (stage 2) tail loss prob to 10E-<x> for QDB", 5 },
  { "--betae",   eslARG_INT,  NULL,    NULL, "0<n<50",   NULL, "--betas",       NULL, "set final   (stage N) tail loss prob to 10E-<x> for QDB", 5 },
  { "--taus",    eslARG_INT,  NULL,    NULL, "0<n<50",   NULL,"--hbanded,--taue",NULL,"set initial (stage 2) tail loss prob to 10E-<x> for HMM banding", 5 },
  { "--taue",    eslARG_INT,  NULL,    NULL, "0<n<50",   NULL,"--hbanded,--taus",NULL,"set final   (stage N) tail loss prob to 10E-<x> for HMM banding", 5 },
  /* Output options */
  { "--regress", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test data to file <f>", 6 },
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump parsetrees to file <f>",  6 },
  /* Other options */
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "arrest after start: for debugging MPI under gdb",   7 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  int           be_verbose;	/* one-line-per-seq summary */
  int           nseq;           /* which number sequence this is in file (only valid in serial mode) */
  int           nstages;        /* number of stages of alignment we'll do */
  int           s;              /* which stage we're on [0..nstages-1], 0 = stage 1, 1 = stage 2 etc. */
  double       *beta;           /* [0..nstages-1] beta for each stage, NULL if not doing QDB */
  double       *tau;            /* [0..nstages-1] tau  for each stage, NULL if not doing HMM banding */
  float        *s1_sc;          /* [0..cfg->nseq] scores for stage 1 parses, filled in 1st output_result() call*/
  ESL_STOPWATCH *s1_w;          /* stopwatch for timing stage 1 */
  ESL_STOPWATCH *w;             /* stopwatch for timing stages 2+ */
  ESL_RANDOMNESS *r;            /* source of randomness for generating sequences */
  int           infmt;		/* format code for input seqfile */
  int           ncm;            /* number CM we're on */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  char         *cmfile;	        /* name of input CM file  */ 
  CMFILE       *cmfp;		/* open input CM file stream       */
  ESL_SQFILE   *sqfp;           /* open sequence input file stream */
  FILE         *tracefp;	/* optional output for parsetrees  */
  FILE         *regressfp;	/* optional output for regression test  */
  FILE         *ofp;	        /* name of output sequence file  */ 

};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "score RNA covariance model against sequences";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); 

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int summarize_align_options(const struct cfg_s *cfg, CM_t *cm);
static int get_sequences(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int i_am_mpi_master, seqs_to_aln_t **ret_seqs_to_aln);
static int determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker);
static int add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
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
      puts("\nexpert miscellaneous options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for source of input sequences:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\nstage 2 alignment options, to compare to stage 1 (D&C non-banded), are:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for testing multiple tau/beta values for --hbanded/--qdb:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      puts("\noutput options are:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 
      puts("\n  other options:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\n  where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  /* Check for incompatible option combinations I don't know how to disallow with esl_getopts */
  /* Can't have --betas and --betae without a --qdb* option */
  if ((! esl_opt_IsDefault(go, "--betas")) && (! esl_opt_IsDefault(go, "--betae")))
    if(! ((esl_opt_GetBoolean(go, "--qdb")) || (esl_opt_GetBoolean(go, "--qdbsmall"))))
    {
	printf("Error parsing options, --betas and --betae combination requires either --qdb or --qdbsmall.\n");
	exit(1);
      }
  if ((! esl_opt_IsDefault(go, "--betas")) && (! esl_opt_IsDefault(go, "--betae")))
    if((esl_opt_GetInteger(go, "--betas")) > (esl_opt_GetInteger(go, "--betae")))
      {
	printf("Error parsing options, --betas <n> argument must be less than --betae <n> argument.\n");
	exit(1);
      }
  if ((! esl_opt_IsDefault(go, "--taus")) && (! esl_opt_IsDefault(go, "--taue")))
    if((esl_opt_GetInteger(go, "--taus")) > (esl_opt_GetInteger(go, "--taue")))
      {
	printf("Error parsing options, --taus <n> argument must be less than --taue <n> argument.\n");
	exit(1);
      }


  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.infmt      = eslSQFILE_UNKNOWN;      /* autodetect sequence file format by default. */ 
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if (esl_opt_GetBoolean(go, "-1")) cfg.be_verbose = FALSE;        
  else                              cfg.be_verbose = TRUE;        
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.sqfp       = NULL;                   /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ofp        = NULL;                   /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tracefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.regressfp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.nseq       = 0;		           /* this counter is incremented in masters */
  cfg.nstages    = 0;                      /* set in init_master_cfg() in masters, stays 0 for workers */
  cfg.s          = 0;                      /* initialized to 0 in init_master_cfg() in masters, stays 0 for workers */
  cfg.beta       = NULL;                   /* alloc'ed and filled in init_master_cfg() in masters, stays NULL in workers */
  cfg.tau        = NULL;                   /* alloc'ed and filled in init_master_cfg() in masters, stays NULL in workers */
  cfg.s1_sc      = NULL;                   /* alloc'ed and filled in first call to output_result() in masters, stays NULL in workers */
  cfg.s1_w       = NULL;                   /* created in init_master_cfg in masters, stays NULL in workers */
  cfg.w          = NULL;                   /* created in init_master_cfg in masters, stays NULL in workers */
  cfg.r          = NULL;                   /* created in init_master_cfg in masters, stays NULL in workers */

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");

  /* This is our stall point, if we need to wait until we get a
   * debugger attached to this process for debugging (especially
   * useful for MPI):
   */
  while (cfg.do_stall); 

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

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you run mpirun?)");

      if (cfg.my_rank > 0)  mpi_worker(go, &cfg);
      else 		    mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      serial_master(go, &cfg);
    }

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.tracefp   != NULL) {
      printf("Parsetrees saved in file %s.\n", esl_opt_GetString(go, "--tfile"));
      fclose(cfg.tracefp);
    }
    if (cfg.regressfp   != NULL) {
      printf("Regression data (parsetrees) saved in file %s.\n", esl_opt_GetString(go, "--regress"));
      fclose(cfg.regressfp);
    }
    if (cfg.ofp       != NULL) { 
      printf("Sequences scored against the CM saved in file %s.\n", esl_opt_GetString(go, "--outfile"));
      fclose(cfg.ofp);
    }
    if (cfg.s1_sc     != NULL) free(cfg.s1_sc);
  }
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.beta      != NULL) free(cfg.beta);
  if (cfg.tau       != NULL) free(cfg.tau);
  if (cfg.s1_w      != NULL) esl_stopwatch_Destroy(cfg.s1_w);
  if (cfg.w         != NULL) esl_stopwatch_Destroy(cfg.w);
  if (cfg.r         != NULL) esl_randomness_Destroy(cfg.r);
  esl_getopts_Destroy(go);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->cmfile      - command line arg 1
 *    cfg->infmt       - format of input file
 * Sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->abc         - digital input alphabet
 *    cfg->ofp         - optional output sequence file
 *    cfg->sqfp        - optional input sequence file
 *    cfg->tracefp     - optional output file
 *    cfg->regressfp   - optional output file
 *    cfg->r           - source of randomness 
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

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* optionally, open the input sequence file */
  if(! esl_opt_IsDefault(go, "--infile")) { 
    status = esl_sqfile_Open(esl_opt_GetString(go, "--infile"), cfg->infmt, NULL, &(cfg->sqfp));
    if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--infile"));
    else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--infile"));
    else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can$.1Žòùt autodetect stdin or .gz."); 
    else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
    cfg->infmt = cfg->sqfp->format;
  }
  /* optionally, open output file */
  if (esl_opt_GetString(go, "--outfile") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "--outfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --outfile output file %s\n", esl_opt_GetString(go, "--outfile"));
  } 

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tracefp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	cm_Fail("Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open regression file */
  if (esl_opt_GetString(go, "--regress") != NULL) {
    if ((cfg->regressfp = fopen(esl_opt_GetString(go, "--regress"), "w")) == NULL) 
	cm_Fail("Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
    }

  /* create RNG */
  if (! esl_opt_IsDefault(go, "--seed")) 
    cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "--seed"));
  else cfg->r = esl_randomness_CreateTimeseeded();

  if (cfg->r       == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create random number generator: probably out of memory");
  return eslOK;
}

/* init_shared_cfg() 
 * Shared initialization of cfg
 *
 * Sets/creates:
 *    cfg->ncm
 *    cfg->nstages
 *    cfg->beta
 *    cfg->tau
 *    cfg->w
 *    cfg->s1_w
 *    
 * Because this is called from an MPI worker, it cannot print; 
 * it must return error messages, not print them.
 */
static int
init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int s;

  cfg->ncm = 0;
  /* determine the number of stages and beta and tau values for each stage */
  cfg->beta = NULL;
  cfg->tau  = NULL;

  if((! (esl_opt_IsDefault(go, "--betas"))) && (! (esl_opt_IsDefault(go, "--betae"))))
    {
      cfg->nstages = esl_opt_GetInteger(go, "--betae") - esl_opt_GetInteger(go, "--betas") + 2;
      ESL_ALLOC(cfg->beta, sizeof(double) * cfg->nstages);
      cfg->beta[0] = 0.; /* this won't matter b/c stage 1 is non-QDB */
      cfg->beta[1] = pow(10., (-1. * esl_opt_GetInteger(go, "--betas")));
      for(s = 2; s < cfg->nstages; s++)
	cfg->beta[s] = cfg->beta[(s-1)] / 10.;
    }
  else if((! (esl_opt_IsDefault(go, "--taus"))) && (! (esl_opt_IsDefault(go, "--taue"))))
    {
      cfg->nstages = esl_opt_GetInteger(go, "--taue") - esl_opt_GetInteger(go, "--taus") + 2;
      ESL_ALLOC(cfg->tau, sizeof(double) * cfg->nstages);
      cfg->tau[1] = pow(10., (-1. * esl_opt_GetInteger(go, "--taus")));
      for(s = 2; s < cfg->nstages; s++)
	cfg->tau[s] = cfg->tau[(s-1)] / 10.;
    }
  else
    {
      cfg->nstages = 2;
      if(esl_opt_GetBoolean(go, "--qdb") || esl_opt_GetBoolean(go, "--qdbsmall"))
	{
	  ESL_ALLOC(cfg->beta, sizeof(double) * cfg->nstages);
	  cfg->beta[0] = 0.; /* this won't matter b/c stage 1 is non-QDB */
	  cfg->beta[1] = esl_opt_GetReal(go, "--beta");
	}
      else if(esl_opt_GetBoolean(go, "--hbanded"))
	{
	  ESL_ALLOC(cfg->tau, sizeof(double) * cfg->nstages);
	  cfg->tau[0] = 0.; /* this won't matter b/c stage 1 is non-banded */
	  cfg->tau[1] = esl_opt_GetReal(go, "--tau");
	}
    }  
  cfg->s1_w  = esl_stopwatch_Create();
  cfg->w     = esl_stopwatch_Create();

  if (cfg->s1_w    == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create stopwatch: probably out of memory");
  if (cfg->w       == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create stopwatch: probably out of memory");

  return eslOK;

 ERROR:
  return status;
}

/* serial_master()
 * The serial version of cmscore.
 * Score each sequence against the CM with specified
 * scoring algorithms.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[cmERRBUFSIZE];
  CM_t    *cm;
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, holds seqs, parsetrees, CP9 traces, postcodes */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      /* get sequences, either generate them (--emit (default) or --random) or read them (--infile) */
      if((status = get_sequences(go, cfg, errbuf, cm, FALSE, &seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  
      /* align sequences cfg->nstages times */
      for(cfg->s = 0; cfg->s < cfg->nstages; cfg->s++) 
	{
	  /* start timing */
	  if(cfg->s == 0) esl_stopwatch_Start(cfg->s1_w);
	  else            esl_stopwatch_Start(cfg->w);
	  
	  /* initialize the flags/options/params of the CM for the current stage */
	  if((status = initialize_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);

	  /* align all sequences, keep scores in sc */
	  if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  
	  /* stop timing, and output result */
	  if(cfg->s == 0) esl_stopwatch_Stop(cfg->s1_w);
	  else            esl_stopwatch_Stop(cfg->w);
	  if ((status = output_result(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);

	  /* clean up, free everything in seqs_to_aln but the sqs, which we'll reuse for each stage */
	  FreePartialSeqsToAln(seqs_to_aln, FALSE, TRUE, TRUE, TRUE, TRUE);
	                                 /* sq,    tr,  cp9_tr,post, sc  */ 
	}
      FreeSeqsToAln(seqs_to_aln); 
      FreeCM(cm);
      printf("//\n");
    }
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmscore.
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
  int nseq_sent;

  seqs_to_aln_t  *all_seqs_to_aln    = NULL;
  seqs_to_aln_t  *worker_seqs_to_aln = NULL;
  int            *seqidx         = NULL;

  char     errbuf[cmERRBUFSIZE];
  MPI_Status mpistatus; 
  int      n;
  int      i;

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
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
      
      /* align sequences cfg->nstages times */
      for(cfg->s = 0; cfg->s < cfg->nstages; cfg->s++) 
	{
	  have_work     = TRUE;	/* TRUE while work remains  */
	  /* Start timing. */
	  if(cfg->s == 0) esl_stopwatch_Start(cfg->s1_w);
	  else            esl_stopwatch_Start(cfg->w);
	  
	  /* initialize the flags/options/params of the CM for the current stage */
	  if((status = initialize_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);

	  if(cfg->s == 0) {
	    /* get sequences, either generate them (--emit (default) or --random) or read them (--infile) */
	    if((status = get_sequences(go, cfg, errbuf, cm, TRUE, &all_seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  }

	  /* determine number of sequences per worker, depends on the stage (b/c stage 1 is slow, stage 2+ may be fast) */
	  determine_nseq_per_worker(go, cfg, cm, &nseq_per_worker); /* this func dies internally if there's some error */
	  ESL_DPRINTF1(("cfg->s: %d nseq_per_worker: %d\n", cfg->s, nseq_per_worker));
	  
	  wi = 1;
	  nseq_sent = 0;
	  while (have_work || nproc_working)
	    {
	      if (have_work) 
		{
		  if(nseq_sent < all_seqs_to_aln->nseq) {
		      nseq_this_worker = (nseq_sent + nseq_per_worker <= all_seqs_to_aln->nseq) ? 
			nseq_per_worker : (all_seqs_to_aln->nseq - nseq_sent);
		  }
		  else {
		    have_work = FALSE;
		    ESL_DPRINTF1(("MPI master has run out of sequences to read (having read %d)\n", all_seqs_to_aln->nseq));
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
		      if(worker_seqs_to_aln->postcode1!= NULL) free(worker_seqs_to_aln->postcode1);
		      if(worker_seqs_to_aln->postcode2!= NULL) free(worker_seqs_to_aln->postcode2);
		      if(worker_seqs_to_aln->sc       != NULL) free(worker_seqs_to_aln->sc);
		      free(worker_seqs_to_aln);
		    }
		  else	/* worker reported an error. Get the errbuf. */
		    {
		      if (MPI_Unpack(buf, bn, &pos, errbuf, cmERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack of errbuf failed");
		      ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
		    }
		  nproc_working--;
		}
	      
	      if (have_work)
		{   
		  /* send new alignment job */
		  ESL_DPRINTF1(("MPI master is sending sequences to align to worker %d\n", wi));
		  if ((status = cm_seqs_to_aln_MPISend(all_seqs_to_aln, nseq_sent, nseq_this_worker, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
		  seqidx[wi] = nseq_sent;
		  nseq_sent += nseq_this_worker;
		  wi++;
		  nproc_working++;
		}
	    }
	  /* send workers the message that we're done with this stage */
	  for (wi = 1; wi < cfg->nproc; wi++) 
	    if ((status = cm_seqs_to_aln_MPISend(NULL, 0, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");

	  /* stop timing, and output result */
	  if(cfg->s == 0) {
	      esl_stopwatch_Stop(cfg->s1_w);
	      esl_stopwatch_MPIReduce(cfg->s1_w, 0, MPI_COMM_WORLD);
	  }
	  else {
	      esl_stopwatch_Stop(cfg->w);
	      esl_stopwatch_MPIReduce(cfg->w, 0, MPI_COMM_WORLD);
	  }
	  /* we have all worker's results for this stage output the results */
	  if ((status = output_result(go, cfg, errbuf, cm, all_seqs_to_aln)) != eslOK) cm_Fail(errbuf);

	  /* clean up, free everything in all_seqs_to_aln but the sqs, which we'll reuse for each stage */
	  if(all_seqs_to_aln->tr       != NULL) { for (i=0; i < all_seqs_to_aln->nseq; i++) if(all_seqs_to_aln->tr[i] != NULL)       { FreeParsetree(all_seqs_to_aln->tr[i]);  all_seqs_to_aln->tr[i] = NULL; } }
	  if(all_seqs_to_aln->cp9_tr   != NULL) { for (i=0; i < all_seqs_to_aln->nseq; i++) if(all_seqs_to_aln->cp9_tr[i] != NULL)   { CP9FreeTrace(all_seqs_to_aln->cp9_tr[i]); all_seqs_to_aln->tr[i] = NULL; } }
	  if(all_seqs_to_aln->postcode1!= NULL) { for (i=0; i < all_seqs_to_aln->nseq; i++) if(all_seqs_to_aln->postcode1[i] != NULL) { free(all_seqs_to_aln->postcode1[i]); all_seqs_to_aln->tr[i] = NULL; } }
	  if(all_seqs_to_aln->postcode2!= NULL) { for (i=0; i < all_seqs_to_aln->nseq; i++) if(all_seqs_to_aln->postcode2[i] != NULL) { free(all_seqs_to_aln->postcode2[i]); all_seqs_to_aln->tr[i] = NULL; } }
	  for (i=0; i < all_seqs_to_aln->nseq; i++) all_seqs_to_aln->sc[i] = IMPOSSIBLE;
	}
      ESL_DPRINTF1(("MPI master: done with this CM.\n"));

      FreeSeqsToAln(all_seqs_to_aln); 
      FreeCM(cm);
      printf("//\n");
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
  char          errbuf[cmERRBUFSIZE];
  seqs_to_aln_t *seqs_to_aln = NULL;
  int           do_free_tr = TRUE;
  int           do_free_cp9_tr = TRUE;

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
  if (xstatus == eslOK) { if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK)   xstatus = status;  }
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
      cfg->ncm++;  
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));

      /* align sequences cfg->nstages times */
      for(cfg->s = 0; cfg->s < cfg->nstages; cfg->s++)
	{
	  /* Start timing. */
	  if(cfg->s == 0) esl_stopwatch_Start(cfg->s1_w);
	  else            esl_stopwatch_Start(cfg->w);

	  /* initialize the flags/options/params of the CM for current stage */
	  if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      
	  while((status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln)) == eslOK)
	    {
	      ESL_DPRINTF1(("worker %d: has received alignment job, nseq: %d\n", cfg->my_rank, seqs_to_aln->nseq));
	      /* align all sequences */
	      if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	      ESL_DPRINTF1(("worker %d: has gathered alignment results\n", cfg->my_rank));
	      
	      /* clean up, free everything in seqs_to_aln but the scores, and maybe the parsetrees or 
	       * cp9_traces (only if --regress or --tfile enabled though), which we'll pass back to the master */
	      do_free_tr = do_free_cp9_tr = TRUE;
	      if((! esl_opt_IsDefault(go, "--regress")) || (! esl_opt_IsDefault(go, "--tfile"))) do_free_tr = do_free_cp9_tr = FALSE;
	      FreePartialSeqsToAln(seqs_to_aln, TRUE, do_free_tr, do_free_cp9_tr, TRUE, FALSE);
                                             /* sq,   tr,         cp9_tr,         post, sc   */ 

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
	      
	    }
	  if(status == eslEOD) ESL_DPRINTF1(("worker %d: has seen message to stop for this stage with this CM.\n", cfg->my_rank));
	  else goto ERROR;

	  /* stop timing */
	  if(cfg->s == 0) {
	      esl_stopwatch_Stop(cfg->s1_w);
	      esl_stopwatch_MPIReduce(cfg->s1_w, 0, MPI_COMM_WORLD);
	  }
	  else {
	      esl_stopwatch_Stop(cfg->w);
	      esl_stopwatch_MPIReduce(cfg->w, 0, MPI_COMM_WORLD);
	  }
	}
      FreeCM(cm);
    }
  if (status == eslEOD) ESL_DPRINTF1(("worker %d told CMs are done.\n", cfg->my_rank));
  else goto ERROR;
  
  if (wbuf != NULL) free(wbuf);
  return;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
  return;
}
#endif /*HAVE_MPI*/

static int
output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln)
{
  int status;
  int i;

  /* print the parsetrees to regression file or parse file */
  for(i = 0; i < seqs_to_aln->nseq; i++)
    {
      if (cfg->regressfp != NULL) 
	{
	  fprintf(cfg->regressfp, "> %s\n", seqs_to_aln->sq[i]->name);
	  if(esl_opt_GetBoolean(go,"--hmmviterbi")) 
	    {
	      ESL_DASSERT1((seqs_to_aln->cp9_tr != NULL));
	      fprintf(cfg->regressfp, "  SCORE : %.2f bits\n", CP9TraceScore(cm->cp9, seqs_to_aln->sq[i]->dsq, seqs_to_aln->cp9_tr[i]));
	      CP9PrintTrace(cfg->regressfp, seqs_to_aln->cp9_tr[i], cm->cp9, seqs_to_aln->sq[i]->dsq);
	    }
	  else
	    {
	      ESL_DASSERT1((seqs_to_aln->tr != NULL));
	      fprintf(cfg->regressfp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE));
	      ParsetreeDump(cfg->regressfp, seqs_to_aln->tr[i], cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); /* NULLs are dmin, dmax */
	    }
	  fprintf(cfg->regressfp, "//\n");
	}
      if (cfg->tracefp != NULL) 
	{
	  fprintf(cfg->tracefp, "> %s\n", seqs_to_aln->sq[i]->name);
	  if(esl_opt_GetBoolean(go,"--hmmviterbi")) 
	    {
	      ESL_DASSERT1((seqs_to_aln->cp9_tr != NULL));
	      fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", CP9TraceScore(cm->cp9, seqs_to_aln->sq[i]->dsq, seqs_to_aln->cp9_tr[i]));
	      CP9PrintTrace(cfg->tracefp, seqs_to_aln->cp9_tr[i], cm->cp9, seqs_to_aln->sq[i]->dsq);
	    }
	  else
	    {
	      ESL_DASSERT1((seqs_to_aln->tr != NULL));
	      fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE));
	      ParsetreeDump(cfg->tracefp, seqs_to_aln->tr[i], cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); /* NULLs are dmin, dmax */
	    }
	  fprintf(cfg->tracefp, "//\n");
	}
    }

  /* print info about scores of parsetrees */
  if(cfg->s == 0) /* store the scores, only */
    {
      ESL_ALLOC(cfg->s1_sc, sizeof(float) * seqs_to_aln->nseq);
      for(i = 0; i < seqs_to_aln->nseq; i++) {
	cfg->s1_sc[i] = seqs_to_aln->sc[i];
	/*cfg->s1_sc[i] = ParsetreeScore(cm, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE);
	  ESL_DASSERT1(((fabs(cfg->s1_sc[i] - seqs_to_aln->sc[i])) < 1e-4));*/
      }
    }
  else /* if(cfg->s > 0) we don't do the comparison test for stage 0 */
    {
      /* Compare parsetrees from stage 1 and stage s (current stage) and collect stats */
      double diff_sc = 0.; /* difference in summed parse scores for this stage versus stage 1 */
      int    diff_ct = 0.; /* number of parses different between this stage and stage 1 */
      double spdup;        /* speed-up versus stage 1 */

      for(i = 0; i < seqs_to_aln->nseq; i++)
	{
	  /* TO DO: write function that inside actually_align_targets() takes
	   * a CP9 parse, and converts it to a CM parsetree */
	  if(esl_opt_GetBoolean(go, "-i"))
	    printf("%-12s S1: %.3f S%d: %.3f diff: %.3f\n", seqs_to_aln->sq[i]->name, cfg->s1_sc[i], (cfg->s+1), seqs_to_aln->sc[i], (fabs(cfg->s1_sc[i] - seqs_to_aln->sc[i])));
	  if(fabs(cfg->s1_sc[i] -  seqs_to_aln->sc[i]) > 0.0001) {
	    diff_ct++;
	    diff_sc += fabs(cfg->s1_sc[i] - seqs_to_aln->sc[i]);
	  }
	}
      /* Print summary for this stage versus stage 1 */ 
      printf("Results summary for stage %d:\n", (cfg->s+1));
      printf("---------------------------------\n");
      printf("Number seqs aligned:     %d\n", seqs_to_aln->nseq);
      esl_stopwatch_Display(stdout, cfg->s1_w, "Stage  1 time:           ");
      printf("Stage %2d time:           ", (cfg->s+1));
      esl_stopwatch_Display(stdout, cfg->w, "");
      spdup = cfg->s1_w->user / cfg->w->user;
      printf("Speedup (user):          %.2f\n", spdup);

      if(! esl_opt_GetBoolean(go, "--hmmviterbi"))
	{
	  printf("Avg bit score diff:      %.2f\n", (diff_sc / ((float) seqs_to_aln->nseq)));
	  if(diff_ct == 0)
	    printf("Avg sc diff(>1e-4):      %.2f\n", 0.);
	  else
	    printf("Avg sc diff(>1e-4):      %.2f\n", (diff_sc / ((float) diff_ct)));
	  printf("Num   diff (>1e-4):      %d\n", (diff_ct));
	  printf("Fract diff (>1e-4):      %.5f\n", (((float) diff_ct) / ((float) seqs_to_aln->nseq)));
	  printf("---------------------------------\n");
	  printf("\n");
	}
    }
  return eslOK;

 ERROR:
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
  int status;

  if((status = ActuallyAlignTargets(cm, errbuf, seqs_to_aln,
				    NULL, NULL, 0,  /* we're not aligning search hits */
				    FALSE, 0, TRUE, NULL)) != eslOK) goto ERROR;

  return eslOK;
  
  ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_search_workunit\n", cfg->my_rank));
  FreeCM(cm);
  return status;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults
 * for the specified stage alignment. We only set flags and 
 * a few parameters. ConfigCM() configures the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  /* Some stuff we do no matter what stage we're on */

  /* set up params/flags/options of the CM */
  if(cfg->beta != NULL) cm->beta = cfg->beta[cfg->s];
  if(cfg->tau  != NULL) cm->tau  = cfg->tau[cfg->s];

  /* enable option to check parsetree score against the alignment score */
  cm->align_opts  |= CM_ALIGN_CHECKPARSESC;

  /* Update cm->config_opts and cm->align_opts based on command line options */
  if(esl_opt_GetBoolean(go, "-l"))
    {
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cm->config_opts |= CM_CONFIG_HMMEL;
    }
  if(esl_opt_GetBoolean(go, "--sub"))         
    { 
      cm->align_opts  |=  CM_ALIGN_SUB;
      cm->align_opts  &= ~CM_ALIGN_CHECKPARSESC; /* parsetree score won't match aln score */
    }
  /*if(  esl_opt_GetBoolean(go, "-i"))         cm->align_opts  |= CM_ALIGN_TIME;*/
  if(! esl_opt_GetBoolean(go, "--iins"))     cm->config_opts |= CM_CONFIG_ZEROINSERTS;
    
  /* do stage 1 specific stuff */
  if(cfg->s == 0) { /* set up stage 1 alignment we'll compare all other stages to */
    cm->align_opts |= CM_ALIGN_SMALL;
    /* only one option allows cmscore NOT to do standard CYK as first stage aln */
    if(esl_opt_GetBoolean(go, "--qdbboth")) { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }
    /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
     * set local mode, make cp9 HMM, calculate QD bands etc. 
     */
    ConfigCM(cm, NULL, NULL);
  }
  else { /* cfg->s > 0, we're at least on stage 2, 
	    don't call ConfigCM() again, only info that may change is QDBs, and align_opts */
    /* Clear QDBs if they exist */
    if(cm->flags & CMH_QDB) {
      free(cm->dmin);
      free(cm->dmax);
      cm->dmin = NULL;
      cm->dmax = NULL;
      cm->flags &= ~CMH_QDB;
    }
    cm->align_opts  = 0;  /* clear alignment options from previous stage */
    cm->config_opts = 0; /* clear configure options from previous stage */

    if(esl_opt_GetBoolean(go, "--hbanded"))     cm->align_opts  |= CM_ALIGN_HBANDED;
    if(esl_opt_GetBoolean(go, "--hmmviterbi"))  cm->align_opts  |= CM_ALIGN_HMMVITERBI;
    if(esl_opt_GetBoolean(go, "--hsafe"))       cm->align_opts  |= CM_ALIGN_HMMSAFE;
    if(esl_opt_GetBoolean(go, "--scoreonly"))   cm->align_opts  |= CM_ALIGN_SCOREONLY;
    if(esl_opt_GetBoolean(go, "--qdb") || esl_opt_GetBoolean(go, "--qdbsmall")) {                    
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
      /* calc QDBs for this stage */
      ConfigQDB(cm);
    }
    /* only one way stage 2+ alignment will be D&C, if --qdbsmall was enabled */
    if(esl_opt_GetBoolean(go, "--qdbsmall"))  cm->align_opts  |= CM_ALIGN_SMALL;
  }
  if(cfg->my_rank == 0) 
    {
      if(cfg->s == 0 && cfg->ncm == 1) {
	if     (! esl_opt_IsDefault(go, "--infile")) printf("Sequence mode: file (%s)\n", esl_opt_GetString(go, "--infile"));
	else if(! esl_opt_IsDefault(go, "--random")) printf("Sequence mode: random (seed: %ld)\n", esl_randomness_GetSeed(cfg->r));
	else                                         printf("Sequence mode: CM emitted (seed: %ld)\n", esl_randomness_GetSeed(cfg->r));
      }
      if(cfg->s == 0) { 
	printf("=================================\n");
	printf("CM: ");
	if(cm->name == NULL) printf("%d\n", (cfg->ncm+1));
	else printf("%s\n", cm->name);
	printf("=================================\n");
      }
      printf("Stage %2d alignment:\n", (cfg->s+1));
      summarize_align_options(cfg, cm);
    }
  return eslOK;
}

/* Function: summarize_align_options
 * Date:     EPN, Wed Jan 17 09:08:18 2007
 * Purpose:  Print out alignment options in pretty format. 
 */
int summarize_align_options(const struct cfg_s *cfg, CM_t *cm)
{
  printf("---------------------------------\n");
  /* Algorithm */
  if(cm->align_opts & CM_ALIGN_INSIDE)
    printf("Algorithm:               Inside\n");
  else if(cm->align_opts & CM_ALIGN_HMMVITERBI) 
    printf("Algorithm:               CP9 HMM Viterbi\n");
  else if(cm->align_opts & CM_ALIGN_SCOREONLY)
    printf("Algorithm:               CYK Standard (score only)\n");
  else if(cm->align_opts & CM_ALIGN_SMALL)
    printf("Algorithm:               CYK D&C\n");
  else 
    printf("Algorithm:               CYK Standard\n");

  /* Bands */
  if(cm->align_opts & CM_ALIGN_HBANDED)
    {
      if(cm->align_opts & CM_ALIGN_SUMS)
	printf("Bands:                   CP9 HMM (sums)\n");
      else
	printf("Bands:                   CP9 HMM\n"); 
      printf("Tail loss:               %g\n", cm->tau);
    }
  else if(cm->align_opts & CM_ALIGN_QDB)
    {
      printf("Bands:                   QDB\n"); 
      printf("Tail loss:               %g\n", cm->beta);
    }
  else
    {
      printf("Bands:                   None\n"); 
      printf("Tail loss:               0.0\n");
    }

  /* Locality */
  if(cm->config_opts & CM_CONFIG_LOCAL)
    printf("Local:                   Yes\n");
  else
    printf("Local:                   No\n");

  /* Sub mode? */
  if(cm->align_opts & CM_ALIGN_SUB)
    printf("Sub mode.\n");

  if(cfg->ncm > 1) printf("---------------------------------\n");
  return eslOK;
}

/* get_sequences()
 * Get sequences to score by either generating them
 * or reading them from a file depending on the input options.
 *
 * Sequences are allocated slightly different if the MPI master
 * calls this function, to allow us to store them after receiving
 * them back from workers in any order.
 */
static int
get_sequences(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int i_am_mpi_master, seqs_to_aln_t **ret_seqs_to_aln)
{
  int status = eslOK;
  int do_emit   = FALSE;
  int do_random = FALSE;
  int do_infile = FALSE;
  seqs_to_aln_t *seqs_to_aln = NULL;
  double *dnull = NULL;
  int i;

  if  (  esl_opt_IsDefault(go, "--emit") && esl_opt_IsDefault(go, "--random") && esl_opt_IsDefault(go, "--infile")) do_emit = TRUE;
  if  (! esl_opt_IsDefault(go, "--emit"))   do_emit = TRUE;
  if  (! esl_opt_IsDefault(go, "--random")) do_random = TRUE;
  if  (! esl_opt_IsDefault(go, "--infile")) do_infile = TRUE;

  ESL_DASSERT1(((do_emit + do_random + do_infile) == 1));

  if(do_emit) {
    seqs_to_aln = CMEmitSeqsToAln(cfg->r, cm, cfg->ncm, esl_opt_GetInteger(go, "--emit"), i_am_mpi_master);
  }
  else if(do_random) {
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);
    seqs_to_aln = RandomEmitSeqsToAln(cfg->r, cm->abc, dnull, cfg->ncm, esl_opt_GetInteger(go, "--random"), cm->clen, i_am_mpi_master);
    free(dnull);
  }
  else if(do_infile)
    {
      seqs_to_aln = CreateSeqsToAln(100, i_am_mpi_master);
      if((status = ReadSeqsToAln(cfg->abc, cfg->sqfp, 0, TRUE, seqs_to_aln, i_am_mpi_master)) != eslEOF) 
	cm_Fail("Error reading sqfile: %s\n", esl_opt_GetString(go, "--infile"));
      /* rewind the sqfile so we can read the seqs again */
      esl_sqio_Rewind(cfg->sqfp); /* we may be searching this file again with another CM */
    }    
  else cm_Fail("get_sequences() error, !do_emit, !do_random and !do_infile.");

  /* optionally, print out the sequences to outfile */
  if(cfg->ofp != NULL) {
    for(i = 0; i < seqs_to_aln->nseq; i++)
      if((esl_sqio_Write(cfg->ofp, seqs_to_aln->sq[i], eslSQFILE_FASTA)) != eslOK) cm_Fail("Error writing unaligned sequences to %s.", esl_opt_GetString(go, "--outfile"));
    fprintf(cfg->ofp, "//\n");
  }

  *ret_seqs_to_aln = seqs_to_aln;
  return eslOK;

 ERROR:
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
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

  if(worker_seqs->postcode1 != NULL) {
    if(master_seqs->postcode1 == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned postcodes, master->postcode1 is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->postcode1[x] == NULL); 
      master_seqs->postcode1[x] = worker_seqs->postcode1[(x-offset)];
    }
  }

  if(worker_seqs->postcode2 != NULL) {
    if(master_seqs->postcode2 == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned postcodes, master->postcode2 is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->postcode2[x] == NULL); 
      master_seqs->postcode2[x] = worker_seqs->postcode2[(x-offset)];
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

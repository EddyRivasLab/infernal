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

#ifdef HAVEDEVOPTS
#define ALGOPTS  "--nonbanded,--qdb,--qdbsmall,--qdbboth,--hbanded,--hmmviterbi,--hmmforward"  /* Exclusive choice for scoring algorithms ifdef #HAVE_DEVOPTS*/
#else
#define ALGOPTS  "--nonbanded,--hbanded,--hmmviterbi,--hmmforward"                             /* Exclusive choice for scoring algorithms */
#endif
#define SEQOPTS  "--emit,--random,--infile"                                 /* Exclusive choice for sequence input */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-n",        eslARG_INT,     "10", NULL, "n>0",     NULL,      NULL,  "--infile", "generate <n> sequences",  1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,     "--sub", "align locally w.r.t. the model",         1 },
  { "-s",        eslARG_INT,     NULL, NULL, "n>0",     NULL,      NULL,  "--infile", "set random number seed to <n>", 1 },
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print individual timings & scores, not just a summary", 1 },
  { "--sub",      eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL, "-l,--search", "build sub CM for columns b/t HMM predicted start/end points", 1 },
  { "--mxsize",  eslARG_REAL, "512.0", NULL, "x>0.",    NULL,      NULL,        NULL, "set maximum allowable DP matrix size to <x> Mb", 1 },
  /* 4 --p* options below are hopefully temporary b/c if we have E-values for the CM using a certain cm->pbegin, cm->pend,
   * changing those values in cmsearch invalidates the E-values, so we should pick hard-coded values for cm->pbegin cm->pend */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      "-l",  "--pbegin", "set all local begins as equiprobable", 1 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      "-l",    "--pend", "set all local end probs to <x>", 1 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      "-l",        NULL, "set aggregate local begin prob to <x>", 1 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      "-l",        NULL, "set aggregate local end prob to <x>", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program",                    1 },  
#endif
  /* options for generating/reading sequences to score */
  { "--emit",    eslARG_NONE,"default",NULL,"n>0",   SEQOPTS,      NULL,     SEQOPTS, "emit <n> sequences from each CM", 2 },
  { "--random",  eslARG_NONE,   FALSE, NULL,"n>0",   SEQOPTS,      NULL,     SEQOPTS, "emit <n> random seq from cm->null model", 2},
  { "--infile",  eslARG_INFILE,  NULL, NULL, NULL,   SEQOPTS,      NULL,     SEQOPTS, "read sequences to align from file <s>", 2 },
  { "--outfile", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,  "--infile", "save seqs to file <s> in FASTA format", 2 },
  { "--Lmin",    eslARG_INT,    FALSE, NULL,"0<n<=1000000",NULL,"--random,--Lmax", NULL, "with --random, specify minimum length of random sequences as <n>", 2},
  { "--Lmax",    eslARG_INT,    FALSE, NULL,"0<n<=1000000",NULL,"--random,--Lmin", NULL, "with --random, specify maximum length of random sequences as <n>", 2},
  { "--pad",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--emit,--search", NULL, "with --emit, pad (W-L) residues on each side of emitted sequence", 2},
  /* Stage 2 algorithm options */
  { "--hbanded", eslARG_NONE,"default",NULL, NULL,  ALGOPTS,      NULL,        NULL, "compare d&c optimal CYK versus HMM banded CYK", 3 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,"--hbanded",      NULL, "set tail loss prob for --hbanded to <x>", 3 },
  { "--aln2bands",eslARG_NONE, FALSE, NULL, NULL,      NULL,"--hbanded,--search",NULL, "w/--hbanded derive HMM bands w/o scanning Forward/Backward", 3 },
  { "--hsafe",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded","--search", "realign (non-banded) seqs with HMM banded CYK score < 0 bits", 3 },
  { "--nonbanded",eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,"--search", "compare divide and conquer (d&c) versus standard non-banded CYK", 3 },
  { "--scoreonly",eslARG_NONE,  FALSE, NULL, NULL,      NULL,"--nonbanded","--tfile,--search", "with --nonbanded, do only score, save memory", 3 },
  { "--hmmviterbi",eslARG_NONE, FALSE, NULL, NULL,   ALGOPTS,      NULL,       NULL, "align to a CM Plan 9 HMM with the Viterbi algorithm", 3 },
#ifdef HAVE_DEVOPTS
  { "--qdb",     eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare non-banded d&c versus QDB standard CYK", 3 },
  { "--qdbsmall",eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,  "--search", "compare non-banded d&c versus QDB d&c", 3 },
  { "--qdbboth", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,  "--search", "compare        QDB d&c versus QDB standard CYK", 3 },
  { "--beta",    eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      NULL,        NULL, "set tail loss prob for QDB to <x>", 3 },
#endif
  /* Options for testing multiple rounds of banded alignment, stage 2->N alignment */
  { "--taus",    eslARG_INT,  NULL,    NULL, "0<n<=30", NULL,"--hbanded,--taue",NULL,"set initial (stage 2) tail loss prob to 1E-<x> for HMM banding", 4 },
  { "--taue",    eslARG_INT,  NULL,    NULL, "0<n<=30", NULL,"--hbanded,--taus",NULL,"set final   (stage N) tail loss prob to 1E-<x> for HMM banding", 4 },
#ifdef HAVE_DEVOPTS
  { "--betas",   eslARG_INT,  NULL,    NULL, "0<n<=30", NULL, "--betae",       NULL, "set initial (stage 2) tail loss prob to 1E-<x> for QDB", 4 },
  { "--betae",   eslARG_INT,  NULL,    NULL, "0<n<=30", NULL, "--betas",       NULL, "set final   (stage N) tail loss prob to 1E-<x> for QDB", 4 },
#endif
  /* Verbose output files/debugging */
  { "--search",  eslARG_NONE,  FALSE,  NULL, NULL,      NULL,      NULL,        NULL, "run algorithms in scanning search mode", 5 },
  { "--hmmforward",eslARG_NONE, FALSE, NULL, NULL,   ALGOPTS,"--search",        NULL, "align to a CM Plan 9 HMM with the Forward algorithm", 5 },
  { "--inside",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--search",        NULL, "with --search use inside instead of CYK", 5 },
  { "--old",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use old hmm to cm band calculation algorithm", 5},
  { "--regress", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test data to file <f>", 5 },
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,  "--search", "dump parsetrees to file <f>",  5 },
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb",   5 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  int           nseq;           /* which number sequence this is in file (only valid in serial mode) */
  int           nstages;        /* number of stages of alignment we'll do */
  int           s;              /* which stage we're on [0..nstages-1], 0 = stage 1, 1 = stage 2 etc. */
  double       *beta;           /* [0..nstages-1] beta for each stage, NULL if not doing QDB */
  double       *tau;            /* [0..nstages-1] tau  for each stage, NULL if not doing HMM banding */
  float        *s1_sc;          /* [0..cfg->nseq] scores for stage 1 parses, filled in 1st output_result() call*/
  ESL_STOPWATCH *s1_w;          /* stopwatch for timing stage 1 */
  ESL_STOPWATCH *s_w;           /* stopwatch for timing stages 2+ */
  ESL_RANDOMNESS *r;            /* source of randomness for generating sequences */
  int           infmt;		/* format code for input seqfile */
  int           ncm;            /* number CM we're on */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* info for the comlog we'll add to the cmfiles */
  char            *ccom;        /* command line used in this execution of cmscore */
  char            *cdate;       /* date of this execution of cmscore */

  /* Masters only (i/o streams) */
  char         *cmfile;	        /* name of input CM file  */ 
  CMFILE       *cmfp;		/* open input CM file stream       */
  ESL_SQFILE   *sqfp;           /* open sequence input file stream */
  FILE         *ofp;            /* output file (default is stdout) */
  FILE         *tracefp;	/* optional output for parsetrees  */
  FILE         *regressfp;	/* optional output for regression test  */
  FILE         *sfp;	        /* optional output for sequence file  */ 
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

static int process_align_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int process_cmscore_search_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int initialize_cm_for_align(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int initialize_cm_for_search(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int get_sequences(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int i_am_mpi_master, seqs_to_aln_t **ret_seqs_to_aln);
static int dispatch_search_for_cmscore(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, float *ret_sc);
extern int get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static void print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq);
static void print_stage_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg);
static void print_seq_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg);
static int print_align_options(const struct cfg_s *cfg, CM_t *cm);
static int print_search_options(const struct cfg_s *cfg, CM_t *cm);
#ifdef HAVE_MPI
static int determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker);
static int add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset);
#endif

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory error, stopwatch not created.\n");
  esl_stopwatch_Start(w);
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
      puts("\noptions for source of input sequences:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nstage 2 alignment options, to compare to stage 1 (D&C non-banded), are:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
#ifdef HAVE_DEVOPTS
      puts("\noptions for testing multiple tau/beta values for --hbanded/--qdb:");
#else
      puts("\noptions for testing multiple tau values for --hbanded:");
#endif
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\nverbose output files and debugging:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
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
  if ((! esl_opt_IsDefault(go, "--taus")) && (! esl_opt_IsDefault(go, "--taue")))
    if((esl_opt_GetInteger(go, "--taus")) > (esl_opt_GetInteger(go, "--taue")))
      {
	printf("Error parsing options, --taus <n> argument must be less than --taue <n> argument.\n");
	exit(1);
      }
#ifdef HAVE_DEVOPTS
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
#endif

  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.infmt      = eslSQFILE_UNKNOWN;      /* autodetect sequence file format by default. */ 
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.sqfp       = NULL;                   /* opened in init_master_cfg() in masters, stays NULL for workers */
  stdout        = NULL;                   /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.sfp        = NULL;                   /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tracefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.regressfp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.nseq       = 0;		           /* this counter is incremented in masters */
  cfg.nstages    = 0;                      /* set in init_master_cfg() in masters, stays 0 for workers */
  cfg.s          = 0;                      /* initialized to 0 in init_master_cfg() in masters, stays 0 for workers */
  cfg.beta       = NULL;                   /* alloc'ed and filled in init_master_cfg() in masters, stays NULL in workers */
  cfg.tau        = NULL;                   /* alloc'ed and filled in init_master_cfg() in masters, stays NULL in workers */
  cfg.s1_sc      = NULL;                   /* alloc'ed and filled in first call to output_result() in masters, stays NULL in workers */
  cfg.s1_w       = NULL;                   /* created in init_master_cfg in masters, stays NULL in workers */
  cfg.s_w        = NULL;                   /* created in init_master_cfg in masters, stays NULL in workers */
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

      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you run mpirun?)");

      if (cfg.my_rank > 0)  mpi_worker(go, &cfg);
      else {
	cm_banner(stdout, argv[0], banner);
	mpi_master(go, &cfg);
      }
      esl_stopwatch_Stop(w);
      esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      cm_banner(stdout, argv[0], banner);
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.tracefp   != NULL) {
      printf("# Parsetrees saved in file %s.\n", esl_opt_GetString(go, "--tfile"));
      fclose(cfg.tracefp);
    }
    if (cfg.regressfp   != NULL) {
      printf("# Regression data (parsetrees) saved in file %s.\n", esl_opt_GetString(go, "--regress"));
      fclose(cfg.regressfp);
    }
    if (cfg.sfp       != NULL) { 
      printf("# Sequences scored against the CM saved in file %s.\n", esl_opt_GetString(go, "--outfile"));
      fclose(cfg.sfp);
    }
    if (cfg.s1_sc     != NULL) free(cfg.s1_sc);
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
    esl_stopwatch_Destroy(w);
    if (stdout       != NULL) fclose(stdout);
  }
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.beta      != NULL) free(cfg.beta);
  if (cfg.tau       != NULL) free(cfg.tau);
  if (cfg.s1_w      != NULL) esl_stopwatch_Destroy(cfg.s1_w);
  if (cfg.s_w       != NULL) esl_stopwatch_Destroy(cfg.s_w);
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
    if ((cfg->sfp = fopen(esl_opt_GetString(go, "--outfile"), "w")) == NULL) 
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
  if (! esl_opt_IsDefault(go, "-s")) 
    cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
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
 *    cfg->s_w
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

  if((! (esl_opt_IsDefault(go, "--taus"))) && (! (esl_opt_IsDefault(go, "--taue"))))
    {
      cfg->nstages = esl_opt_GetInteger(go, "--taue") - esl_opt_GetInteger(go, "--taus") + 2;
      ESL_ALLOC(cfg->tau, sizeof(double) * cfg->nstages);
      cfg->tau[1] = pow(10., (-1. * esl_opt_GetInteger(go, "--taus")));
      for(s = 2; s < cfg->nstages; s++)
	cfg->tau[s] = cfg->tau[(s-1)] / 10.;
    }
#ifdef HAVE_DEVOPTS
  else if((! (esl_opt_IsDefault(go, "--betas"))) && (! (esl_opt_IsDefault(go, "--betae"))))
    {
      cfg->nstages = esl_opt_GetInteger(go, "--betae") - esl_opt_GetInteger(go, "--betas") + 2;
      ESL_ALLOC(cfg->beta, sizeof(double) * cfg->nstages);
      cfg->beta[0] = 0.; /* this won't matter b/c stage 1 is non-QDB */
      cfg->beta[1] = pow(10., (-1. * esl_opt_GetInteger(go, "--betas")));
      for(s = 2; s < cfg->nstages; s++)
	cfg->beta[s] = cfg->beta[(s-1)] / 10.;
    }
#endif
  else
    {
      cfg->nstages = 2;
      if(esl_opt_GetBoolean(go, "--hbanded"))
	{
	  ESL_ALLOC(cfg->tau, sizeof(double) * cfg->nstages);
	  cfg->tau[0] = 0.; /* this won't matter b/c stage 1 is non-banded */
	  cfg->tau[1] = esl_opt_GetReal(go, "--tau");
	}
#ifdef HAVE_DEVOPTS
      else if(esl_opt_GetBoolean(go, "--qdb") || esl_opt_GetBoolean(go, "--qdbsmall"))
	{
	  ESL_ALLOC(cfg->beta, sizeof(double) * cfg->nstages);
	  cfg->beta[0] = 0.; /* this won't matter b/c stage 1 is non-QDB */
	  if(! esl_opt_IsDefault(go, "--beta")) cm->beta_qdb = esl_opt_GetReal(go, "--beta");
	  /* else cm->beta_qdb will be equal to beta from CM file */
	  cfg->beta[1] = cm->beta_qdb;
	}
#endif
    }  
  cfg->s1_w  = esl_stopwatch_Create();
  cfg->s_w   = esl_stopwatch_Create();

  if (cfg->s1_w    == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create stopwatch: probably out of memory");
  if (cfg->s_w     == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create stopwatch: probably out of memory");

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
  if ((status = print_run_info (go, cfg, errbuf))  != eslOK) cm_Fail(errbuf);

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      /* align sequences cfg->nstages times */
      for(cfg->s = 0; cfg->s < cfg->nstages; cfg->s++) 
	{
	  /* start timing */
	  if(cfg->s == 0) esl_stopwatch_Start(cfg->s1_w);
	  else            esl_stopwatch_Start(cfg->s_w);
	  
	  if(esl_opt_GetBoolean(go, "--search")) {
	    /* initialize the flags/options/params of the CM for the current stage */
	    if((status = initialize_cm_for_search(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
	  }
	  else { 
	    /* initialize the flags/options/params of the CM for the current stage */
	    if((status = initialize_cm_for_align(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
	  }

	  if(cfg->s == 0) {
	    /* get sequences, either generate them (--emit (default) or --random) or read them (--infile) */
	    if((status = get_sequences(go, cfg, errbuf, cm, FALSE, &seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	    /* print the per-cm info */
	    print_cm_info (go, cfg, errbuf, cm, seqs_to_aln->nseq);
	  }

	  if(esl_opt_GetBoolean(go, "--search")) {
	    /* align all sequences, keep scores in seqs_to_aln->sc */
	    if ((status = process_cmscore_search_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  }
	  else {
	    /* align all sequences, keep scores in seqs_to_aln->sc */
	    if ((status = process_align_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  }

	  /* stop timing, and output result */
	  if(cfg->s == 0) esl_stopwatch_Stop(cfg->s1_w);
	  else            esl_stopwatch_Stop(cfg->s_w);
	  if ((status = output_result(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);

	  /* clean up, free everything in seqs_to_aln but the sqs, which we'll reuse for each stage */
	  FreePartialSeqsToAln(seqs_to_aln, FALSE, TRUE, TRUE, TRUE, TRUE);
	                                 /* sq,    tr,  cp9_tr,post, sc  */ 
	}
      FreeSeqsToAln(seqs_to_aln); 
      FreeCM(cm);
    }
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmscore.
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
  int      wi_error = 0;                /* worker index that sent back an error message, if an error occurs */

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
  if (xstatus == eslOK) { if ((seqidx = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((status = print_run_info(go, cfg, errbuf))  != eslOK) xstatus = status; }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) cm_Fail(errbuf);
  ESL_DPRINTF1(("MPI master is initialized\n"));

  for (wi = 0; wi < cfg->nproc; wi++) seqidx[wi] = 0;

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

  while (xstatus == eslOK && CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      cfg->ncm++;  
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* align sequences cfg->nstages times */
      for(cfg->s = 0; cfg->s < cfg->nstages; cfg->s++) 
	{
	  if(xstatus == eslOK) have_work = TRUE;	
	  /* have_work stays TRUE while work remains, if a worker has seen an error and sent us back an error message xstatus != eslOK,
	   * and have_work will be FALSE for all remaining stages, so no more work is performed, then once we leave this loop, we shut down
	   */

	  /* Start timing. */
	  if(cfg->s == 0) esl_stopwatch_Start(cfg->s1_w);
	  else            esl_stopwatch_Start(cfg->s_w);
	  
	  /* initialize the flags/options/params of the CM for current stage */
	  if(esl_opt_GetBoolean(go, "--search")) { 
	    if((status = initialize_cm_for_search(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
	  }
	  else {
	    if((status = initialize_cm_for_align(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
	  }

	  if(cfg->s == 0) {
	    /* get sequences, either generate them (--emit (default) or --random) or read them (--infile) */
	    if((status = get_sequences(go, cfg, errbuf, cm, TRUE, &all_seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	    /* print the per-cm info */
	    print_cm_info (go, cfg, errbuf, cm, all_seqs_to_aln->nseq);
	  }

	  /* determine number of sequences per worker */
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
		      if (MPI_Unpack(buf, bn, &pos, &xstatus, 1, MPI_INT, MPI_COMM_WORLD)     != 0)     cm_Fail("mpi unpack failed");
		      if (xstatus == eslOK) /* worker reported success. Get the results. */
			{
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
			  have_work = FALSE;
			  wi_error  = wi;
			}
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
	    esl_stopwatch_Stop(cfg->s_w);
	    esl_stopwatch_MPIReduce(cfg->s_w, 0, MPI_COMM_WORLD);
	  }
	  /* if we've got valid results for this stage, output them */
	  if (xstatus == eslOK) {
	    if ((status = output_result(go, cfg, errbuf, cm, all_seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  }
	  ESL_DPRINTF1(("MPI master: done with this stage for this CM. Telling all workers\n"));

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
      fprintf(stdout, "//\n");
    }

  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((status = cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); cm_Fail(errbuf); }
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
	  else            esl_stopwatch_Start(cfg->s_w);

	  /* initialize the flags/options/params of the CM for current stage */
	  if(esl_opt_GetBoolean(go, "--search")) { 
	    if((status = initialize_cm_for_search(go, cfg, errbuf, cm)) != eslOK) goto ERROR;
	  }
	  else {
	    if((status = initialize_cm_for_align(go, cfg, errbuf, cm)) != eslOK) goto ERROR;
	  }
      
	  while((status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln)) == eslOK)
	    {
	      ESL_DPRINTF1(("worker %d: has received alignment job, nseq: %d\n", cfg->my_rank, seqs_to_aln->nseq));

	      if(esl_opt_GetBoolean(go, "--search")) { 
		/* align all sequences, keep scores in seqs_to_aln->sc */
		if ((status = process_cmscore_search_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
		ESL_DPRINTF1(("worker %d: has gathered search results\n", cfg->my_rank));
	      }
	      else {
		/* align all sequences */
		if ((status = process_align_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) goto ERROR;
		ESL_DPRINTF1(("worker %d: has gathered alignment results\n", cfg->my_rank));
	      }

	      /* clean up, free everything in seqs_to_aln but the scores, and maybe the parsetrees or 
	       * cp9_traces (only if --regress or --tfile enabled though), which we'll pass back to the master */
	      do_free_tr = do_free_cp9_tr = TRUE;
	      if((! esl_opt_IsDefault(go, "--regress")) || (! esl_opt_IsDefault(go, "--tfile"))) do_free_tr = do_free_cp9_tr = FALSE;
	      FreePartialSeqsToAln(seqs_to_aln, TRUE, do_free_tr, do_free_cp9_tr, TRUE, FALSE);
                                             /* sq,   tr,         cp9_tr,         post, sc   */ 

	      n = 0;
	      if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
		ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	      n += sz;
	      if (cm_seqs_to_aln_MPIPackSize(seqs_to_aln, 0, seqs_to_aln->nseq, MPI_COMM_WORLD, &sz) != eslOK) 
		ESL_XFAIL(eslFAIL, errbuf, "cm_seqs_to_aln_MPIPackSize() call failed"); 
	      n += sz;
	      if (n > wn) {
		void *tmp;
		ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
		wn = n;
	      }
	      ESL_DPRINTF1(("worker %d: has calculated the alignment results will pack into %d bytes\n", cfg->my_rank, n));
	      status = eslOK;
	      
	      pos = 0;
	      if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
		ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	      if (cm_seqs_to_aln_MPIPack(seqs_to_aln, 0, seqs_to_aln->nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) 
		ESL_XFAIL(eslFAIL, errbuf, "cm_seqs_to_aln_MPIPack() call failed"); 
	      MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	      ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));
	      
	    }
	  if(status == eslEOD) ESL_DPRINTF1(("worker %d: has seen message to stop for this stage with this CM.\n", cfg->my_rank));
	  else ESL_XFAIL(eslFAIL, errbuf, "within CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);

	  /* stop timing */
	  if(cfg->s == 0) {
	      esl_stopwatch_Stop(cfg->s1_w);
	      esl_stopwatch_MPIReduce(cfg->s1_w, 0, MPI_COMM_WORLD);
	  }
	  else {
	      esl_stopwatch_Stop(cfg->s_w);
	      esl_stopwatch_MPIReduce(cfg->s_w, 0, MPI_COMM_WORLD);
	  }
	}
      FreeCM(cm);
    }
  if (status == eslEOD) ESL_DPRINTF1(("worker %d told CMs are done.\n", cfg->my_rank));
  else ESL_XFAIL(eslFAIL, errbuf, "outside CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);
  
  if (wbuf != NULL) free(wbuf);
  return;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

  /* if we get here, this worker failed and sent an error message, now the master knows a worker
   * failed but it has to send the message to all other workers (besides this one) to abort so they 
   * can be shut down cleanly. As currently implemented, this means we have to wait here for that 
   * signal which comes in the form of a special 'empty' work packet that tells us we're done with
   * the current CM, and then a 'empty' CM broadcast that tells us we're done with all CMs in the file.
   */
  for(; cfg->s < cfg->nstages; cfg->s++) {
    status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln);
    /* stop timing (stupid, but necessary to avoid massive reorg of current implementation (poor design)) */
    if(cfg->s == 0) {
      esl_stopwatch_Stop(cfg->s1_w);
      esl_stopwatch_MPIReduce(cfg->s1_w, 0, MPI_COMM_WORLD);
    }
    else {
      esl_stopwatch_Stop(cfg->s_w);
      esl_stopwatch_MPIReduce(cfg->s_w, 0, MPI_COMM_WORLD);
    }
  }
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);
  /* status after each of the above calls should be eslEOD, but if it isn't we can't really do anything 
   * about it b/c we've already sent our error message, so in that scenario the MPI will break uncleanly 
   */
  return;
}
#endif /*HAVE_MPI*/

static int
output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln)
{
  int status;
  int i;
  char time_buf[128];	  /* another string for printing elapsed time */

  /* print the parsetrees to regression file or parse file */
  for(i = 0; i < seqs_to_aln->nseq; i++)
    {
      if (cfg->regressfp != NULL) 
	{
	  fprintf(cfg->regressfp, "> %s\n", seqs_to_aln->sq[i]->name);
	  if(cfg->s > 0 && esl_opt_GetBoolean(go,"--hmmviterbi")) 
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
	  if(cfg->s > 0 && esl_opt_GetBoolean(go,"--hmmviterbi")) 
	    {
	      ESL_DASSERT1((seqs_to_aln->cp9_tr != NULL));
	      fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", CP9TraceScore(cm->cp9, seqs_to_aln->sq[i]->dsq, seqs_to_aln->cp9_tr[i]));
	      CP9PrintTrace(cfg->tracefp, seqs_to_aln->cp9_tr[i], cm->cp9, seqs_to_aln->sq[i]->dsq);
	    }
	  else
	    {
	      ESL_DASSERT1((seqs_to_aln->tr != NULL));
	      if(esl_opt_GetBoolean(go, "--sub")) { 
		fprintf(cfg->tracefp, "  SUB CM PARSE SCORE                               : %.2f bits\n", seqs_to_aln->sc[i]);
		fprintf(cfg->tracefp, "  SUB CM ALIGNMENT MAPPED ONTO ORIG CM PARSE SCORE : %.2f bits\n", ParsetreeScore(cm, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE));
	      }
	      else { 
		fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE));
	      }
	      ParsetreeDump(cfg->tracefp, seqs_to_aln->tr[i], cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); /* NULLs are dmin, dmax */
	    }
	  fprintf(cfg->tracefp, "//\n");
	}
    }

  /* print info about scores of parsetrees */
  if(cfg->s == 0) /* store the scores, only */
    {
      FormatTimeString(time_buf, cfg->s1_w->user, TRUE);
      ESL_ALLOC(cfg->s1_sc, sizeof(float) * seqs_to_aln->nseq);
      for(i = 0; i < seqs_to_aln->nseq; i++) {
	cfg->s1_sc[i] = seqs_to_aln->sc[i];
	/*cfg->s1_sc[i] = ParsetreeScore(cm, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE);
	  ESL_DASSERT1(((fabs(cfg->s1_sc[i] - seqs_to_aln->sc[i])) < 1e-4));*/
      }
      /* Print summary for stage 1 */ 
      print_stage_column_headings(go, cfg);
      fprintf(stdout, "  %5d", (cfg->s+1)); /* stage number */
      if(esl_opt_GetBoolean(go, "--search")) print_search_options(cfg, cm);
      else                                   print_align_options(cfg, cm);
      fprintf(stdout, "  %11s  %6s  %7s  %7s  %6s\n", 
	      time_buf,                  /* time */
	      "-", "-", "-", "-");  /* comparisons with stage 1 are all N/A */
    }
  else /* if(cfg->s > 0) we don't do the comparison test for stage 0 */
    {
      /* Compare parsetrees from stage 1 and stage s (current stage) and collect stats */
      double diff_sc = 0.; /* difference in summed parse scores for this stage versus stage 1 */
      int    diff_ct = 0.; /* number of parses different between this stage and stage 1 */

      FormatTimeString(time_buf, cfg->s_w->user, TRUE);
      if(esl_opt_GetBoolean(go, "-a")) print_seq_column_headings(go, cfg);
      for(i = 0; i < seqs_to_aln->nseq; i++)
	{
	  /* TO DO: write function that inside DispatchAlignments() takes
	   * a CP9 parse, and converts it to a CM parsetree */
	  if(esl_opt_GetBoolean(go, "-a")) 
	    fprintf(stdout, "  %-25.25s  %6d  %11.4f  %11.4f  %10.4f\n", seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, cfg->s1_sc[i], seqs_to_aln->sc[i], cfg->s1_sc[i] - seqs_to_aln->sc[i]);
	  if(fabs(cfg->s1_sc[i] - seqs_to_aln->sc[i]) > 0.01) {
	    diff_ct++;
	    diff_sc += cfg->s1_sc[i] - seqs_to_aln->sc[i]; /* don't take absolute value in case cur stage sc > stage 1 sc, for example with -l --hmmviterbi */
	  }
	  /*////*/if(seqs_to_aln->sc[i] < -10000.) { 
	    ///cm_Fail("cm: %s seq %d got impossible score: %f\n", cm->name, i, seqs_to_aln->sc[i]); 
	    printf("cm: %s seq %d got impossible score: %f\n", cm->name, i, seqs_to_aln->sc[i]); 
	  }
	  if(esl_opt_GetBoolean(go, "--nonbanded") && (fabs(cfg->s1_sc[i] -  seqs_to_aln->sc[i]) >= 0.01))
	    cm_Fail("Non-banded standard CYK score (%.3f bits) differs from D&C CYK score (%.3f bits)", cfg->s1_sc[i], seqs_to_aln->sc[i]);
	}
      if(esl_opt_GetBoolean(go, "-a")) print_stage_column_headings(go, cfg);

      /* Print summary for this stage versus stage 1 */ 
      fprintf(stdout, "  %5d", (cfg->s+1)); /* stage number */
      if(esl_opt_GetBoolean(go, "--search")) print_search_options(cfg, cm);
      else                                   print_align_options(cfg, cm);
      fprintf(stdout, "  %11s  %6.2f  %7d  %7.5f  %6.2f\n", 
	      time_buf,                                            /* time */
	      cfg->s1_w->user / cfg->s_w->user,                    /* speedup versus stage 1 */
	      diff_ct,                                             /* number of seqs with different scores */
	      ((float) diff_ct) / ((float) seqs_to_aln->nseq),     /* fraction that are different */
	      (diff_ct == 0) ? 0. : (diff_sc / ((float) diff_ct)));/* avg score diff for those that are diff */
    }
  return eslOK;

 ERROR:
  return status;
}

/* An alignment work unit consists a seqs_to_aln_t object which contains sequences to align, 
 * and space for their parsetrees, or CP9 traces, and postal codes.
 * The job is to align the sequences and collect alignment scores and possibly
 * create parsetrees or cp9 traces.
 */
static int
process_align_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, 
		 seqs_to_aln_t *seqs_to_aln)
{
  int status;

  if((status = DispatchAlignments(cm, errbuf, seqs_to_aln,
				  NULL, NULL, 0,  /* we're not aligning search hits */
				  FALSE, 0, TRUE, NULL, 
				  esl_opt_GetReal(go, "--mxsize"), stdout)) != eslOK) goto ERROR;

  return eslOK;
  
  ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_align_workunit\n", cfg->my_rank));
  FreeCM(cm);
  return status;
}

/* A search work unit consists a seqs_to_aln_t object which contains sequences to search, 
 * and space for their parsetrees, or CP9 traces, and postal codes, but all we care about
 * is the score of the best hit in each sequence.
 * The job is to search the sequences and return the scores.
 */
static int
process_cmscore_search_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, 
				seqs_to_aln_t *seqs_to_aln)
{
  int status;
  int i;
  float *scA; /* will store scores of all sequences */

  ESL_ALLOC(scA, sizeof(float) * seqs_to_aln->nseq);
  for(i = 0; i < seqs_to_aln->nseq; i++) { 
    if((status = dispatch_search_for_cmscore(cm, errbuf, seqs_to_aln->sq[i]->dsq, 1, seqs_to_aln->sq[i]->n,
					     esl_opt_GetReal(go, "--mxsize"), &(scA[i]))) != eslOK) goto ERROR;
  }
  seqs_to_aln->sc = scA;
  return eslOK;
  
  ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_cmscore_search_workunit\n", cfg->my_rank));
  FreeCM(cm);
  return status;
}

/* initialize_cm_for_align()
 * Setup the CM based on the command-line options/defaults
 * for the specified stage alignment. We only set flags and 
 * a few parameters. ConfigCM() configures the CM.
 */
static int
initialize_cm_for_align(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  /* Some stuff we do no matter what stage we're on */
  cm->align_opts  = 0;  /* clear alignment options from previous stage */
  cm->config_opts = 0;  /* clear configure options from previous stage */

  /* set up params/flags/options of the CM */
  if(cfg->beta != NULL) cm->beta_qdb = cfg->beta[cfg->s];
  if(cfg->tau  != NULL) cm->tau  = cfg->tau[cfg->s];

  /* enable option to check parsetree score against the alignment score */
  cm->align_opts  |= CM_ALIGN_CHECKPARSESC;

  /* Update cm->config_opts and cm->align_opts based on command line options */
  if(esl_opt_GetBoolean(go, "-l")) {
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  /* BEGIN (POTENTIALLY) TEMPORARY BLOCK */
  int nstarts, nexits, nd;
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
  /* END (POTENTIALLY) TEMPORARY BLOCK */

  if(esl_opt_GetBoolean(go, "--sub")) {        
    cm->align_opts  |=  CM_ALIGN_SUB;
    cm->align_opts  &= ~CM_ALIGN_CHECKPARSESC; /* parsetree score won't match aln score */
  }
    
  /* do stage 1 specific stuff */
  if(cfg->s == 0) { /* set up stage 1 alignment we'll compare all other stages to */
    cm->align_opts |= CM_ALIGN_SMALL;
#ifdef HAVE_DEVOPTS
    /* only one option allows cmscore NOT to do standard CYK as first stage aln */
    if(esl_opt_GetBoolean(go, "--qdbboth")) { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }
#endif
    /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
     * set local mode, make cp9 HMM, calculate QD bands etc. 
     */
    ConfigCM(cm, FALSE); /* FALSE says don't bother calc'ing W, we won't need it */
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

    if(esl_opt_GetBoolean(go, "--hbanded"))     cm->align_opts  |= CM_ALIGN_HBANDED;
    if(esl_opt_GetBoolean(go, "--old"))         cm->align_opts  |= CM_ALIGN_HMM2IJOLD;
    if(esl_opt_GetBoolean(go, "--hmmviterbi"))  cm->align_opts  |= CM_ALIGN_HMMVITERBI;
    if(esl_opt_GetBoolean(go, "--hsafe"))       cm->align_opts  |= CM_ALIGN_HMMSAFE;
    if(esl_opt_GetBoolean(go, "--scoreonly"))   cm->align_opts  |= CM_ALIGN_SCOREONLY;
#ifdef HAVE_DEVOPTS
    if(esl_opt_GetBoolean(go, "--qdb") || esl_opt_GetBoolean(go, "--qdbsmall") || esl_opt_GetBoolean(go, "--qdbboth")) {                    
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
      /* calc QDBs for this stage */
      ConfigQDB(cm);
    }
    /* only one way stage 2+ alignment will be D&C, if --qdbsmall was enabled */
    if(esl_opt_GetBoolean(go, "--qdbsmall"))  cm->align_opts  |= CM_ALIGN_SMALL;
#endif
  }
  return eslOK;
}


/* initialize_cm_for_search()
 * Setup the CM based on the command-line options/defaults
 * for the specified stage alignment. We only set flags and 
 * a few parameters. ConfigCM() configures the CM.
 */
static int
initialize_cm_for_search(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  /* Some stuff we do no matter what stage we're on */
    cm->search_opts  = 0;  /* clear alignment options from previous stage */
    cm->config_opts = 0;  /* clear configure options from previous stage */

  /* set up params/flags/options of the CM */
  if(cfg->beta != NULL) cm->beta_qdb = cfg->beta[cfg->s];
  if(cfg->tau  != NULL) cm->tau  = cfg->tau[cfg->s];

  /* Update cm->config_opts and cm->align_opts based on command line options */
  if(esl_opt_GetBoolean(go, "-l")) {
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  /* BEGIN (POTENTIALLY) TEMPORARY BLOCK */
  int nstarts, nexits, nd;
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
  /* END (POTENTIALLY) TEMPORARY BLOCK */
  if(esl_opt_GetBoolean(go, "--inside"))      cm->search_opts  |= CM_SEARCH_INSIDE;
  cm->search_opts |= CM_SEARCH_NOQDB;
  cm->search_opts |= CM_SEARCH_NOALIGN;
      
  /* do stage 1 specific stuff */
  if(cfg->s == 0) { /* set up stage 1 alignment we'll compare all other stages to */
    cm->search_opts |= CM_SEARCH_NOQDB;
    /* configure the CM for search based on cm->config_opts and cm->align_opts.
     * set local mode, make cp9 HMM, calculate QD bands etc. 
     */
    ConfigCM(cm, TRUE); /* TRUE says calculate W */
  }
  else { /* cfg->s > 0, we're at least on stage 2, 
	    don't call ConfigCM() again, only info that may change is QDBs, and search_opts */
    /* Clear QDBs if they exist */
    if(cm->flags & CMH_QDB) {
      free(cm->dmin);
      free(cm->dmax);
      cm->dmin = NULL;
      cm->dmax = NULL;
      cm->flags &= ~CMH_QDB;
    }

    if(esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts  |= CM_SEARCH_HBANDED;
    if(esl_opt_GetBoolean(go, "--aln2bands"))   cm->search_opts  |= CM_SEARCH_HMMALNBANDS;
    if(esl_opt_GetBoolean(go, "--old"))         cm->search_opts  |= CM_SEARCH_HMM2IJOLD;
    if(esl_opt_GetBoolean(go, "--hmmviterbi"))  { 
      cm->search_opts  |= CM_SEARCH_HMMVITERBI;
      cm->search_opts  &= ~CM_SEARCH_INSIDE;
    }
    if(esl_opt_GetBoolean(go, "--hmmforward")) {
      cm->search_opts  |= CM_SEARCH_HMMFORWARD;
      cm->search_opts  &= ~CM_SEARCH_INSIDE;
    }
#ifdef HAVE_DEVOPTS    
    if(esl_opt_GetBoolean(go, "--qdb")) {
      cm->search_opts &= ~CM_CONFIG_NOQDB;
      cm->config_opts |= CM_CONFIG_QDB;
      /* calc QDBs for this stage */
      ConfigQDB(cm);
    }
#endif
  }
  /* create scan matrix (for all rounds, including first round, round 0) */
  if(cm->flags & CMH_SCANMATRIX) { cm_FreeScanMatrixForCM(cm); cm->smx = NULL; }
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) ESL_FAIL(eslFAIL, errbuf, "initialize_cm_for_search(), CreateScanMatrixForCM() call failed.");
    
  /* create the search info (for all rounds, including first round, round 0) */
  if(cm->si != NULL) { FreeSearchInfo(cm->si, cm); cm->si = NULL; }

  CreateSearchInfo(cm, SCORE_CUTOFF, 0.0, -1.); /* 0.0 is score threshold, it's irrelevant we find best score per seq regardless */
  ValidateSearchInfo(cm, cm->si);

  return eslOK;
}

/* Function: print_align_options
 * Date:     EPN, Wed Jan 17 09:08:18 2007
 * Purpose:  Print out alignment options in pretty format. 
 */
int print_align_options(const struct cfg_s *cfg, CM_t *cm)
{
  /* algorithm */
  if     (cm->align_opts & CM_ALIGN_HMMVITERBI) fprintf(stdout, "  %7s", "hmm-vit");
  else if(cm->align_opts & CM_ALIGN_SMALL)      fprintf(stdout, "  %7s", "cyk-d&c");
  else                                          fprintf(stdout, "  %7s", "cyk-std");
  /* bands and beta/tau*/
  if     (cm->align_opts & CM_ALIGN_HBANDED)    fprintf(stdout, "  %5s  %6.0e", "hmm", cm->tau);
  else if(cm->align_opts & CM_ALIGN_QDB)        fprintf(stdout, "  %5s  %6.0e", "qdb", cm->beta_qdb);
  else                                          fprintf(stdout, "  %5s  %6s", "-", "-");

  return eslOK;
}

/* Function: print_search_options
 * Date:     EPN, Fri Jan 25 09:26:32 2008
 * Purpose:  Print out search options in pretty format. 
 */
int print_search_options(const struct cfg_s *cfg, CM_t *cm)
{
  /* algorithm */
  if     (cm->search_opts & CM_SEARCH_HMMVITERBI) fprintf(stdout, "  %7s", "hmm-vit");
  else if(cm->search_opts & CM_SEARCH_HMMVITERBI) fprintf(stdout, "  %7s", "hmm-fwd");
  else if(cm->search_opts & CM_SEARCH_INSIDE)     fprintf(stdout, "  %7s", "inside");
  else                                            fprintf(stdout, "  %7s", "cyk");
  /* bands and beta/tau*/
  if     (cm->search_opts & CM_SEARCH_HBANDED)    fprintf(stdout, "  %5s  %6.0e", "hmm", cm->tau);
  else if(! (cm->search_opts & CM_SEARCH_NOQDB))  fprintf(stdout, "  %5s  %6.0e", "qdb", cm->beta_qdb);
  else                                            fprintf(stdout, "  %5s  %6s", "-", "-");

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
  int do_emit   =    esl_opt_GetBoolean(go, "--emit");
  int do_random =    esl_opt_GetBoolean(go, "--random");
  int do_infile = (! esl_opt_IsDefault (go, "--infile"));
  int nseq      =    esl_opt_GetInteger(go, "-n");
  seqs_to_aln_t *seqs_to_aln = NULL;
  double        *dnull = NULL;
  int            i;
  int            safe_windowlen = cm->clen * 2;
  double       **gamma = NULL;
  double        *Ldistro = NULL;
  int            lengths_specified = FALSE;
  int            L, Lmin, Lmax;

  assert((do_emit + do_random + do_infile) == 1);
  ESL_DASSERT1(((do_emit + do_random + do_infile) == 1));
  ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
  for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
  esl_vec_DNorm(dnull, cm->abc->K);

  if(do_emit) {
    seqs_to_aln = CMEmitSeqsToAln(cfg->r, cm, cfg->ncm, nseq, esl_opt_GetBoolean(go, "--pad"), dnull, i_am_mpi_master);
  }
  else if(do_random) {
    lengths_specified = (esl_opt_IsDefault(go, "--Lmin") && esl_opt_IsDefault(go, "--Lmax")) ? FALSE : TRUE;
    if(!lengths_specified) { /* set random sequence length distribution as length distribution of generative CM, obtained from QDB calc */
      while(!(BandCalculationEngine(cm, safe_windowlen, DEFAULT_HS_BETA, TRUE, NULL, NULL, &(gamma), NULL))) {
	safe_windowlen *= 2;
	if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "Error trying to get gamma[0], safe_windowlen big: %d\n", safe_windowlen);
	FreeBandDensities(cm, gamma);
      }
      Ldistro = gamma[0];
      Lmax    = safe_windowlen;
    }
    else { /* --Lmin <n1> and --Lmax <n2> enabled, set length distribution as uniform from <n1>..<n2> inclusive */
      Lmin = esl_opt_GetInteger(go, "--Lmin"); 
      Lmax = esl_opt_GetInteger(go, "--Lmax"); 
      if(Lmin > Lmax) ESL_FAIL(eslEINCOMPAT, errbuf, "with --Lmin <n1> and --Lmax <n2>, <n2> must be >= <n1>.");
      ESL_ALLOC(Ldistro, sizeof(double) * (Lmax + 1));
      for(L = 0;    L < Lmin;  L++)  Ldistro[L] = 0.;
      for(L = Lmin; L <= Lmax; L++)  Ldistro[L] = 1. / (float) (Lmax - Lmin + 1);
    }
    seqs_to_aln = RandomEmitSeqsToAln(cfg->r, cm->abc, dnull, cfg->ncm, nseq, Ldistro, Lmax, i_am_mpi_master);

    if(gamma != NULL) FreeBandDensities(cm, gamma);
    else              free(Ldistro);
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
  if(cfg->sfp != NULL) {
    for(i = 0; i < seqs_to_aln->nseq; i++)
      if((esl_sqio_Write(cfg->sfp, seqs_to_aln->sq[i], eslSQFILE_FASTA)) != eslOK) cm_Fail("Error writing unaligned sequences to %s.", esl_opt_GetString(go, "--outfile"));
  }

  free(dnull);
  *ret_seqs_to_aln = seqs_to_aln;
  return eslOK;

 ERROR:
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: dispatch_search_for_cmscore()
 * Incept:   EPN, Fri Jan 25 09:43:52 2008
 *            
 * Purpose:  Given a CM and a sequence, call the correct search algorithm
 *           based on search_info and return the score of the best
 *           scoring sequence in the hit in <ret_sc>. Based on dispatch.c:DispatchSearch(),
 *           but simpler, no filtering is allowed and we don't care about storing
 *           results, all we want is the score of the best hit.
 * 
 * Args:     cm              - the covariance model
 *           errbuf          - char buffer for reporting errors
 *           dsq             - the target sequence (digitized)
 *           i0              - start of target subsequence (often 1, beginning of dsq)
 *           j0              - end of target subsequence (often L, end of dsq)
 *           size_limit      - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           ret_sc          - RETURN: Highest scoring hit from search (even if below cutoff).
 *
 * Returns: eslOK on success. eslERANGE if we're doing HMM banded alignment and requested matrix is too big.
 */
int dispatch_search_for_cmscore(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, float *ret_sc)
{
  int               status;          /* easel status code */
  float             sc;              /* score of best hit in seq */
  int               i, j;            /* subseq start/end points */
  int               sround;          /* round 0 */
  SearchInfo_t     *si = cm->si;     /* the SearchInfo */

  /* convenience pointers to cm->si for this 'filter round' of searching */
  float             cutoff;          /* cutoff for this round, HMM or CM, whichever is relevant for this round */
  int               stype;           /* search type for this round SEARCH_WITH_HMM, SEARCH_WITH_HYBRID, or SEARCH_WITH_CM */
  ScanMatrix_t     *smx;             /* scan matrix for this round, != NULL only if SEARCH_WITH_CM, and must == cm->smx if we're in the final round */
  HybridScanInfo_t *hsi;             /* hybrid scan info for this round, NULL unless stype is SEARCH_WITH_HYBRID */

  /* Contract checks */
  if(!(cm->flags & CMH_BITS))          ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(), CMH_BITS flag down.\n");
  if(si == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(): search info cm->si is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(): dsq is NULL.");
  if(!(cm->flags & CMH_BITS))          ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(): CMH_BITS flag down.\n");
  if(si->nrounds != 0)                 ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(): si->nrounds != 0\n");
  if(si->stype[0] == SEARCH_WITH_HYBRID) ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(): hybrid filtering not yet implemented.\n");

  /* copy info for this round from SearchInfo fi */
  sround = 0;
  cm->search_opts = si->search_opts[sround]; 
  cutoff          = si->sc_cutoff[sround]; /* this will be a bit score regardless of whether the cutoff_type == E_CUTOFF */
  stype           = si->stype[sround];
  smx             = si->smx[sround]; /* may be NULL */
  hsi             = si->hsi[sround]; /* may be NULL */

  /* SEARCH_WITH_HMM section */
  if(stype == SEARCH_WITH_HMM) { 
    /* some SEARCH_WITH_HMM specific contract checks */
    if(cm->cp9 == NULL)                    ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(), trying to use CP9 HMM that is NULL.\n");
    if(!(cm->cp9->flags & CPLAN9_HASBITS)) ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(), trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");
    if(hsi != NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(), SEARCH_WITH_HMM but hsi != NULL.\n");
    if(! ((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD)))
      ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(), search type for this round is SEARCH_WITH_HMM, but CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags are both down.");

    /* Scan the (sub)seq in forward direction w/Viterbi or Forward, find score of best hit and it's endpoint j, then if CM_SEARCH_HMMFORWARD, go backwards to get it's start point, and a better guess at it's score */
    if(cm->search_opts & CM_SEARCH_HMMVITERBI) { 
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, i0, j0, cm->W, cutoff, 
			       NULL,   /* don't store hits */
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       NULL, NULL, NULL,  /* don't return best score at each posn, best scoring posn, or traces */
			       &sc)) != eslOK) return status;
    }
    else if(cm->search_opts & CM_SEARCH_HMMFORWARD) { 
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, i0, j0, cm->W, cutoff, 
			       NULL,   /* don't store hits */
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       NULL,   /* don't return best score at each posn */
			       &j,     /* return end point j of best scoring hit */
			       &sc)) != eslOK) return status;
      if((status = cp9_Backward(cm, errbuf, cm->cp9_mx, dsq, i0, j, cm->W, cutoff, 
				NULL,   /* don't report hits */
				TRUE,   /* we're scanning */
				FALSE,  /* we're not ultimately aligning */
				TRUE,   /* be memory efficient */
				NULL,   /* don't return best score at each posn */
				&i,     /* return start point i of best scoring hit, not used actually */
			       &sc)) != eslOK) return status;
      /* now sc is score of Forward hit from i..j, this is the *probably* the best hit in the sequence i0..j0,
       * but we can't be sure, it's possible that a hit in the cp9_Forward() run that ended at j' != j, had a 
       * lower cumulative score from i0..j' then did i0..j, and then when we found i' that maximized the score
       * from i'..j' the score of i'..j' > sc. This is possible, but I'm not sure how we could test for it, or
       * if we even care.
       */
    }
  }
  else { /* stype == SEARCH_WITH_CM */
    ESL_DASSERT1((stype == SEARCH_WITH_CM));
    if(smx == NULL)                             ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(), SEARCH_WITH_CM but smx == NULL.\n");
    if(hsi != NULL)                             ESL_FAIL(eslEINCOMPAT, errbuf, "dispatch_search_for_cmscore(): SEARCH_WITH_CM, but hsi is NULL\n");

    if(cm->search_opts & CM_SEARCH_HBANDED) {
      if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, i0, j0, cm->cp9b, TRUE, 0)) != eslOK) return status; 
      if(cm->search_opts & CM_SEARCH_INSIDE) { if((status = FastFInsideScanHB(cm, errbuf, dsq, i0, j0, cutoff, NULL, cm->hbmx, size_limit, &sc)) != eslOK) return status; }
      else                                   { if((status = FastCYKScanHB    (cm, errbuf, dsq, i0, j0, cutoff, NULL, cm->hbmx, size_limit, &sc)) != eslOK) return status; }
    }
    else { /* don't do HMM banded search */
      if(cm->search_opts & CM_SEARCH_INSIDE) { if((status = FastIInsideScan(cm, errbuf, smx, dsq, i0, j0, cutoff, NULL, NULL, &sc)) != eslOK) return status; }
      else                                   { if((status = FastCYKScan    (cm, errbuf, smx, dsq, i0, j0, cutoff, NULL, NULL, &sc)) != eslOK) return status; }
    }    
    /* now sc is score of best hit found by the relevant CM scanning algorithm */
  }
  /* don't do alignments */

  if(ret_sc != NULL) *ret_sc = sc;
  return eslOK;
}  

/* Function: print_cm_info
 * Date:     EPN, Fri Jan 25 13:43:28 2008
 *
 * Purpose:  Print per-CM info to output file (stdout unless -o). 
 */
static void
print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq)
{
  fprintf(stdout, "# %4s  %-25s  %6s  %6s  %3s  %6s\n", "idx",  "cm name",                   "strat",  "config", "sub", "nseq"  ); 
  fprintf(stdout, "# %4s  %-25s  %6s  %6s  %3s  %6s\n", "----", "-------------------------", "------", "------", "---", "------"); 
  fprintf(stdout, "# %4d  %-25.25s  %6s  %6s  %3s  %6d\n", 
	  cfg->ncm, 
	  cm->name,
	  (esl_opt_GetBoolean(go, "--search")) ? "search" : "align", 
	  (esl_opt_GetBoolean(go, "-l")) ? "local" : "glocal",
	  (esl_opt_GetBoolean(go, "--sub")) ? "yes" : "no",
	  nseq);
  return;
}

/* Function: print_stage_column_headings()
 * Date:     EPN, Fri Jan 25 15:17:09 2008
 *
 * Purpose:  Print per stage column headings to output file (stdout unless -o). 
 *
 * Returns:  eslOK on success
 */
static void
print_stage_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg)
{
  int do_qdb = FALSE;
#ifdef HAVE_DEVOPTS 
  if((esl_opt_GetBoolean(go, "--qdb")) || (esl_opt_GetBoolean(go, "--qdbsmall")) || (esl_opt_GetBoolean(go, "--qdbboth")) || (! esl_opt_IsDefault(go, "--betas"))) do_qdb = TRUE;
#endif
  fprintf(stdout, "#\n");
  fprintf(stdout, "# %5s  %7s  %5s  %6s  %11s  %32s\n",               "",      "",        "",      "",                         "",            "    comparison with stage 1    ");
  fprintf(stdout, "# %5s  %7s  %5s  %6s  %11s  %32s\n",               "",      "",        "",      "",                         "",            "--------------------------------");
  fprintf(stdout, "# %5s  %7s  %5s  %6s  %11s  %6s  %7s  %7s  %6s\n", "stage", "alg",     "bands", (do_qdb) ? "beta" : "tau",  "run time",    "spdup",  "num dif", "frc dif", "sc dif");
  fprintf(stdout, "# %5s  %7s  %5s  %6s  %11s  %6s  %7s  %7s  %6s\n", "-----", "-------", "-----", "------",                   "-----------", "------", "-------", "-------", "------");
  return;
}


/* Function: print_seq_column_headings()
 * Date:     EPN, Fri Jan 25 15:17:09 2008
 *
 * Purpose:  Print sequence column headings to output file (stdout unless -o). 
 *           These are printed only if -a enabled. 
 *
 * Returns:  eslOK on success
 */
static void
print_seq_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg)
{
  fprintf(stdout, "#\n");
  fprintf(stdout, "# %-25s  %6s  %11s  %5s %2d %2s  %10s\n", "seq name",                  "length", "stage 1 sc",  "stage", cfg->s+1, "sc", "sc dif");
  fprintf(stdout, "# %25s  %6s  %11s  %11s  %10s\n",        "-------------------------",  "------", "-----------", "-----------",           "----------");
  return;
}

/* Function: print_run_info
 * Date:     EPN, Fri Jan 25 13:43:28 2008
 *
 * Purpose:  Print information on this run of cmscore.
 *           Command used to run it, and execution date.
 *
 * Returns:  eslOK on success
 */
static int
print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  int status;
  char *command;
  char *date;

  if((status = get_command(go, errbuf, &command)) != eslOK) return status;
  if((status = GetDate    (errbuf, &date))    != eslOK) return status;

  fprintf(stdout, "%-10s %s\n",  "# command:", command);
  fprintf(stdout, "%-10s %s\n",  "# date:",    date);
  fprintf(stdout, "%-10s %ld\n", "# seed:",    esl_randomness_GetSeed(cfg->r));
  if(cfg->nproc > 1) fprintf(stdout, "# %-8s %d\n", "nproc:", cfg->nproc);
  if     (! esl_opt_IsDefault(go, "--infile")) fprintf(stdout, "%-10s input file (%s)\n", "# mode:", esl_opt_GetString(go, "--infile"));
  else if( esl_opt_GetBoolean(go, "--random")) fprintf(stdout, "%-10s random\n", "# mode:");
  else                                         fprintf(stdout, "%-10s cm emitted\n", "# mode:");

  fprintf(stdout, "#\n");
  free(command);
  free(date);
  return eslOK;
}


/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call cmscore
 *           in <ret_command>.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
int 
get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command)
{
  int status;
  int i;
  char *command = NULL;

  for (i = 0; i < go->argc; i++) { /* copy all command line options and args */
    if((status = esl_strcat(&(command),  -1, go->argv[i], -1)) != eslOK) goto ERROR;
    if(i < (go->argc-1)) if((status = esl_strcat(&(command), -1, " ", 1)) != eslOK) goto ERROR;
  }
  *ret_command = command;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_command(): memory allocation error.");
  return status;
}

#ifdef HAVE_MPI
/* determine_nseq_per_worker()
 * Given a CM, return the number of sequences we think we should send
 * to each worker (we don't know the number of sequences in the file).
 */
static int
determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker)
{
  if     (cm->clen <= 200) *ret_nseq_worker = 5;
  else if(cm->clen <= 400) *ret_nseq_worker = 4;
  else if(cm->clen <= 600) *ret_nseq_worker = 3;
  else if(cm->clen <= 800) *ret_nseq_worker = 2;
  else                     *ret_nseq_worker = 1;
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
#endif /* #ifdef HAVE_MPI */

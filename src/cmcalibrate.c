/* cmcalibrate.c
 * Score a CM against random sequence data to set the statistical 
 * parameters for E-value determination.
 * 
 * EPN, Wed May  2 07:02:52 2007 [Updated for v1.1 (removed CP9 and fthr)
 * EPN, Wed May  2 07:02:52 2007
 * based on HMMER-2.3.2's hmmcalibrate.c from SRE
 *
 ************************************************************
 * @LICENSE@
 ************************************************************
 */

#include "esl_config.h"
#include "config.h"	

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_mpi.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_ratematrix.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#define MPI_FINISHED_EXP    -1 /* message to send to workers */
#define MPI_FINISHED_FILTER -2 /* message to send to workers */
#define REALLYSMALLX        1e-20
#define EXPTAIL_CHUNKLEN    10000 /* sequence chunk length for random sequence searches */

#include "funcs.h"		/* external functions                   */
#include "structs.h"

static ESL_OPTIONS options[] = {
  /* name                type           default env   range     toggles      reqs       incomp  help  docgroup*/
  { "-h",                eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "show brief help on version and usage",   1 },
  { "-s",                eslARG_INT,    "181",  NULL, "n>=0",   NULL,        NULL,        NULL, "set RNG seed to <n> (if 0: one-time arbitrary seed)", 1 },
  { "--forecast",        eslARG_INT,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "don't do calibration, forecast running time with <n> processors", 1 },
  { "--devhelp",         eslARG_NONE,   NULL,   NULL, NULL,     NULL,        NULL,        NULL, "show list of undocumented developer options", 1 },
#ifdef HAVE_MPI
  { "--mpi",            eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "run as an MPI parallel program", 1 },  
  { "--Wchunk",         eslARG_INT,     NULL,   NULL, "n>4",    NULL,     "--mpi",        NULL, "set W multiplier for sequence chunk size to <n>", 1 },
#endif
  /* options for exp tail fitting */
  { "-L",           eslARG_REAL,    "1.5",  NULL, "0.1<=x<=1000.", NULL,NULL,         NULL, "set random seq length to search in Mb to <x>", 2 },
  { "--gtailn",     eslARG_INT,    "250",  NULL, "n>=100",  NULL,        NULL,   "--tailp", "fit the top <n> hits/Mb in histogram for  CM local modes", 2 },
  { "--ltailn",     eslARG_INT,    "750",  NULL, "n>=100",  NULL,        NULL,   "--tailp", "fit the top <n> hits/Mb in histogram for  CM glocal modes", 2 },
  { "--tailp",      eslARG_REAL,    NULL,   NULL, "0.0<x<0.6",NULL,      NULL,        NULL, "set fraction of histogram tail to fit to exp tail to <x>", 2 },
  { "--tailxn",     eslARG_INT,     "1000", NULL, "n>=50",  NULL,        NULL,   "--tailp", "w/--tailp, set max num hits in tail to fit as <n>", 2 },
  { "--beta",       eslARG_REAL,    "1E-15",NULL, "x>0",    NULL,        NULL,  "--no-qdb", "set tail loss prob for QDB to <x>", 2 },
  { "--no-qdb",     eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "do not use QDBs for calibrating CM search modes", 2 },
  { "--hfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save fitted score histogram(s) to file <f>", 2 },
  { "--sfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save survival plot to file <f>", 2 },
  { "--qqfile",     eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save Q-Q plot for score histogram(s) to file <f>", 2 },
  { "--ffile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save lambdas for different tail fit probs to file <f>", 2 },
  /* All options below are developer options, only shown if --devhelp invoked */
  /* Developer option, print extra info */
  { "-v",                eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL, "--forecast", "print arguably interesting info",  101},
#ifdef HAVE_MPI
  /* Developer option, for debugging */
  { "--stall",          eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "arrest after start: for debugging MPI under gdb", 101 },  
#endif
  /* Developer options the average user doesn't need to know about */
  { "-T",           eslARG_REAL,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "set bit sc cutoff for exp tail fitting to <x> [df: -INFTY]", 102 },
  { "--random",     eslARG_NONE,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "use GC content of random null background model of CM",  102},
  { "--gc",         eslARG_INFILE,  NULL,   NULL, NULL,     NULL,        NULL,        NULL, "use GC content distribution from file <f>",  102},
  { "--no-null3", eslARG_NONE,   FALSE,  NULL, NULL,      NULL,        NULL,        NULL, "turn OFF the NULL3 post hoc additional null model", 102 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char              *cmfile;	         /* name of input CM file  */ 
  ESL_RANDOMNESS    *r;
  ESL_ALPHABET      *abc;
  ESL_STOPWATCH     *w_stage;            /* stopwatch for each stage */
  double            *gc_freq;
  ExpInfo_t       ***expAA;              /* the exponential tail info, 1st dim: 1 for each CM, 2nd dim: EXP_NMODES */
  int                ncm;                /* what number CM we're on */
  int                be_verbose;	 /* print extra info? */
  int                cmalloc;            /* number of expAA we have allocated (1st dim) */
  char              *tmpfile;            /* tmp file we're writing to */
  char              *mode;               /* write mode, "w" or "wb"                     */
  double            *dnull;              /* double version of cm->null, for generating random seqs */
  float              sc_cutoff;          /* minimum score of a hit we'll consider (-eslINFINITY by default) */

  /* the HMM that generates sequences for exponential tail fitting */
  int               ghmm_nstates;        /* number of states in the HMM */
  double           *ghmm_sA;             /* start probabilities [0..ghmm_nstates-1] */
  double          **ghmm_tAA;            /* transition probabilities [0..nstates-1][0..nstates-1] */
  double          **ghmm_eAA;            /* emission probabilities   [0..nstates-1][0..abc->K-1] */

  /* number of sequences and the length of each seq for exp tail
   * fitting, set such that: exp_cmN is the number of 10 Kb
   * seqs we'll search for CM local/glocal exponential tail fitting:
   *
   * expN = (esl_opt_GetBoolean(go, "-L")  * 10^6) / EXPTAIL_CHUNKLEN(10000); 
   *
   * We don't search just 1 long sequence (i.e. 1.5 Mb) b/c using
   * sequence lengths above 10 Kb for exp tail calibration can yield
   * millions of hits (for CM searches) before overlaps are removed,
   * which requires a lot of memory.
   */
  int              expN;        /* number of 10 Kb seqs for  local CM exp tail fitting */
  int              expL;        /* the size of seq chunks to search, set as 10,000 (10 Kb) */

  /* info for the comlog we'll add to the cmfiles */
  char            *ccom;               /* command line used in this execution of cmcalibrate */
  char            *cdate;              /* date of this execution of cmcalibrate */

  /* mpi */
  int              do_mpi;
  int              my_rank;
  int              nproc;
  int              do_stall;          /* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CM_FILE         *cmfp;	      /* open input CM file stream       */
  FILE            *exphfp;            /* optional output for exp tail histograms */
  FILE            *expsfp;            /* optional output for exp tail survival plot */
  FILE            *expqfp;            /* optional output for exp tail QQ file */
  FILE            *exptfitfp;         /* optional output for exp tail fit file */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit exponential tails for CM E-values";

static int init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static int   mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int   mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int  initialize_cm   (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int do_local);
static int  initialize_stats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static int  fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, int exp_mode, double *ret_mu, double *ret_lambda, int *ret_nrandhits, float *ret_tailp);
static int  get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int L, ESL_DSQ **ret_dsq);
static int  get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm);
static int  set_dnull(struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int  print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int  print_per_cm_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  print_per_cm_summary        (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, double psec, double asec);
static int  print_exp_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, int expN, int expL, double psec);
static int  print_post_calibration_info (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, double *exp_psecA, double *exp_asecA);
static int  estimate_time_for_exp_round (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_res);
static int  process_exp_search_workunit(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float cutoff, CM_TOPHITS **ret_th);

int
main(int argc, char **argv)
{
  int              status;
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory allocation error, stopwatch could not be created.");
  esl_stopwatch_Start(w);
  struct cfg_s     cfg;
  int total_nt;
  int i;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

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
  if (esl_opt_GetBoolean(go, "--devhelp") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nexponential tail distribution fitting options :");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nundocumented developer options for debugging:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\nundocumented exp tail related developer options:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80);
      puts("\nother undocumented developer options:");
      esl_opt_DisplayHelp(stdout, go, 103, 2, 80);
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nexponential tail distribution fitting options :");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  /* Initialize configuration shared across all kinds of masters
   * and workers in this .c file.
   */
  cfg.cmfile  = esl_opt_GetArg(go, 1);
  if (cfg.cmfile == NULL) cm_Fail("Failed to read <cmfile> argument from command line.");
  cfg.r           = NULL; 
  cfg.abc         = NULL; 
  cfg.w_stage     = NULL; 
  cfg.be_verbose  = FALSE;
  cfg.cmalloc     = 128;
  cfg.tmpfile     = NULL;
  cfg.mode        = NULL;
  cfg.expL        = EXPTAIL_CHUNKLEN; /* 10 Kb chunks are searched */
  
  cfg.dnull        = NULL;
  cfg.ghmm_nstates = 0;
  cfg.ghmm_sA      = NULL;
  cfg.ghmm_tAA     = NULL;
  cfg.ghmm_eAA     = NULL;

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cfg.cmalloc  = 128;
  ESL_ALLOC(cfg.expAA, sizeof(ExpInfo_t **)       * cfg.cmalloc);
  cfg.ncm      = 0;

  cfg.ccom      = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */
  cfg.cdate     = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */

  cfg.do_mpi   = FALSE;
  cfg.my_rank  = 0;
  cfg.nproc    = 0;
  cfg.do_stall = FALSE;
#ifdef HAVE_MPI
  cfg.do_stall = esl_opt_GetBoolean(go, "--stall");
#endif

  /* calculate sequence lengths and quantities for exp tail fitting, */
  total_nt  = (int) (1000000. * esl_opt_GetReal(go, "-L"));

  /* determine the number of 10 Kb chunks (cfg.expL = 10000) to search to reach the totals */
  if(total_nt < (cfg.expL + 1)) cm_Fail("with -L <x>, <x> must be at least %.3f.", cfg.expL / 1000000.);
  cfg.expN = (int) (((float) total_nt  / (float) cfg.expL) + 0.999999); 

  cfg.cmfp     = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.exphfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.expsfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.expqfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.exptfitfp= NULL; /* ALWAYS remains NULL for mpi workers */

  ESL_DASSERT1((EXP_CM_GC  == 0));
  ESL_DASSERT1((EXP_CM_GI  == 1));
  ESL_DASSERT1((EXP_CM_LC  == 2));
  ESL_DASSERT1((EXP_CM_LI  == 3));
  ESL_DASSERT1((EXP_NMODES == 4));

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
      char             errbuf[cmERRBUFSIZE]; /* for error messages in mpi_master() */
      if(esl_opt_IsOn(go, "--forecast")) cm_Fail("--forecast is incompatible with --mpi.");
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you run mpirun?)");

      if (cfg.my_rank > 0)  { status = mpi_worker(go, &cfg); }
      else { 
	cm_banner(stdout, argv[0], banner);
	status = mpi_master(go, &cfg, errbuf);
      }
      /* check status, if eslOK, we continue, else we exit. either way we call MPI_Finalize() */
      if(status == eslOK) { 
	esl_stopwatch_Stop(w);
	esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
	MPI_Finalize();
      }
      else { /* status != eslOK, master has error message in errbuf, worker does not */
	MPI_Finalize();
	if(cfg.my_rank == 0) cm_Fail(errbuf); /* master */
	else                 return 0;        /* worker */
      }
    }
  else
#endif /*HAVE_MPI*/
    {
      cm_banner(stdout, argv[0], banner);
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }

  if(! esl_opt_IsOn(go, "--forecast") && cfg.my_rank == 0) { /* master, serial or mpi */
    /* Rewind the CM file for a second pass.
     * Write a temporary CM file with new stats information in it
     */
    int   status;
    int   cmi;
    CM_t *cm;
    FILE *outfp;
    sigset_t blocksigs;  /* list of signals to protect from             */
    char     errbuf[cmERRBUFSIZE];

    cm_file_Position(cfg.cmfp, (off_t) 0);

    if (esl_FileExists(cfg.tmpfile))                    cm_Fail("Ouch. Temporary file %s appeared during the run.", cfg.tmpfile);
    if ((outfp = fopen(cfg.tmpfile, cfg.mode)) == NULL) cm_Fail("Ouch. Temporary file %s couldn't be opened for writing.", cfg.tmpfile); 
    
    for (cmi = 0; cmi < cfg.ncm; cmi++) {
      if ((status = cm_file_Read(cfg.cmfp, TRUE, &(cfg.abc), &cm)) != eslOK) cm_Fail("Ran out of CMs too early in pass 2");
      if (cm == NULL)                                                        cm_Fail("CM file %s was corrupted? Parse failed in pass 2", cfg.cmfile);

      /* update the cm->comlog info */
      if((status = update_comlog(go, errbuf, cfg.ccom, cfg.cdate, cm)) != eslOK) cm_Fail(errbuf);
	
      if(cm->expA != NULL) { 
	for(i = 0; i < EXP_NMODES;  i++)  free(cm->expA[i]); free(cm->expA);
      }
      ESL_ALLOC(cm->expA, sizeof(ExpInfo_t *) * EXP_NMODES);

      cm->expA   = cfg.expAA[cmi];
      cm->flags |= CMH_EXPTAIL_STATS; 
      if(cfg.cmfp->is_binary) { 
	if ((status = cm_file_WriteBinary(outfp, -1, cm, NULL)) != eslOK) ESL_FAIL(status, errbuf, "binary CM save failed");
      }
      else { 
	if ((status = cm_file_WriteASCII(outfp, -1, cm)) != eslOK) ESL_FAIL(status, errbuf, "CM save failed");
      }
      FreeCM(cm);
    } /* end of from idx = 0 to ncm */
    
    /* Now, carefully remove original file and replace it
     * with the tmpfile. Note the protection from signals;
     * we wouldn't want a user to ctrl-C just as we've deleted
     * their CM file but before the new one is moved.
     */
    cm_file_Close(cfg.cmfp);
    if (fclose(outfp)   != 0)                            cm_Fail("system error during rewrite of CM file");
    if (sigemptyset(&blocksigs) != 0)                    cm_Fail("system error during rewrite of CM file.");;
    if (sigaddset(&blocksigs, SIGINT) != 0)              cm_Fail("system error during rewrite of CM file.");;
    if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   cm_Fail("system error during rewrite of CM file.");;
    if (remove(cfg.cmfile) != 0)                         cm_Fail("system error during rewrite of CM file.");;
    if (rename(cfg.tmpfile, cfg.cmfile) != 0)            cm_Fail("system error during rewrite of CM file.");;
    if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) cm_Fail("system error during rewrite of CM file.");;
    free(cfg.tmpfile);
    
    /* master specific cleaning */
    if (cfg.exphfp   != NULL) { 
      fclose(cfg.exphfp);
      printf("# Histogram of high scoring hits in random seqs saved to file %s.\n", esl_opt_GetString(go, "--hfile"));
    }
    if (cfg.expsfp   != NULL) { 
      fclose(cfg.expsfp);
      printf("# Survival plot for exponential tails saved to file %s.\n", esl_opt_GetString(go, "--sfile"));
    }
    if (cfg.expqfp   != NULL) { 
      fclose(cfg.expqfp);
      printf("# Exponential tail QQ plots saved to file %s.\n", esl_opt_GetString(go, "--qqfile"));
    }
    if (cfg.exptfitfp   != NULL) { 
      fclose(cfg.exptfitfp);
      printf("# Exponential tail fit points saved to file %s.\n", esl_opt_GetString(go, "--ffile"));
    }
    if (cfg.ccom  != NULL) free(cfg.ccom);
    if (cfg.cdate != NULL) free(cfg.cdate);
    if (cfg.expAA != NULL) free(cfg.expAA);
  }

  /* clean up */
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.w_stage   != NULL) esl_stopwatch_Destroy(cfg.w_stage);
  if (cfg.r         != NULL) esl_randomness_Destroy(cfg.r);
  if (cfg.my_rank == 0) { 
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
  }
  if(cfg.ghmm_eAA != NULL) { 
    for(i = 0; i < cfg.ghmm_nstates; i++) free(cfg.ghmm_eAA[i]); 
    free(cfg.ghmm_eAA);
  }
  if(cfg.ghmm_tAA != NULL) { 
    for(i = 0; i < cfg.ghmm_nstates; i++) free(cfg.ghmm_tAA[i]); 
    free(cfg.ghmm_tAA);
  }
  if(cfg.ghmm_sA != NULL) free(cfg.ghmm_sA);

  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;

 ERROR: 
  cm_Fail("Memory allocation error.");
  return 1; /* NEVERREACHED */
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Allocates/sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->exphfp      - optional output file
 *    cfg->expsfp      - optional output file
 *    cfg->expqfp      - optional output file
 *    cfg->gc_freq     - observed GC freqs (if --gc invoked)
 *    cfg->expAA       - the exp tail stats, allocated only
 *    cfg->r           - source of randomness
 *    cfg->tmpfile     - temp file for rewriting cm file
 *    cfg->be_verbose  - print extra info? 
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
  status = cm_file_Open(cfg->cmfile, NULL, FALSE, &(cfg->cmfp), errbuf);
  if      (status == eslENOTFOUND) return status;
  else if (status == eslEFORMAT)   return status;
  else if (status != eslOK)        return status;

  /* optionally, open exp tail histogram file */
  if (esl_opt_GetString(go, "--hfile") != NULL) 
    {
      if ((cfg->exphfp = fopen(esl_opt_GetString(go, "--hfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail histogram save file %s for writing\n", esl_opt_GetString(go, "--hfile"));
    }

  /* optionally, open survival plot */
  if (esl_opt_GetString(go, "--sfile") != NULL) 
    {
      if ((cfg->expsfp = fopen(esl_opt_GetString(go, "--sfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open survival plot save file %s for writing\n", esl_opt_GetString(go, "--sfile"));
    }

  /* optionally, open exp tail QQ plot file */
  if (esl_opt_GetString(go, "--qqfile") != NULL) 
    {
      if ((cfg->expqfp = fopen(esl_opt_GetString(go, "--qqfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail QQ plot save file %s for writing\n", esl_opt_GetString(go, "--qqfile"));
    }

  /* optionally, open exp tail tail fit prob file */
  if (esl_opt_GetString(go, "--ffile") != NULL) 
    {
      if ((cfg->exptfitfp = fopen(esl_opt_GetString(go, "--ffile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail save file %s for writing\n", esl_opt_GetString(go, "--ffile"));
    }

  /* optionally, get distribution of GC content from an input database (default is use cm->null for GC distro) */
  if(esl_opt_GetString(go, "--gc") != NULL) {
    ESL_ALPHABET *tmp_abc = NULL;
    tmp_abc = esl_alphabet_Create(eslRNA);
    ESL_SQFILE      *dbfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--gc"), eslSQFILE_UNKNOWN, NULL, &dbfp);
    if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "No such file: %s.", esl_opt_GetString(go, "--gc")); 
    else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "file: %s format unrecognized.", esl_opt_GetString(go, "--gc")); 
    else if (status != eslOK)      ESL_FAIL(status, errbuf, "Failed to open sequence database file %s, code %d.", esl_opt_GetString(go, "--gc"), status); 
    if((status = GetDBInfo(tmp_abc, dbfp, errbuf, NULL, NULL, &(cfg->gc_freq))) != eslOK) return status; 
    esl_vec_DNorm(cfg->gc_freq, GC_SEGMENTS);
    esl_alphabet_Destroy(tmp_abc);
    esl_sqfile_Close(dbfp); 
  }

  cfg->be_verbose = FALSE;
  if (esl_opt_GetBoolean(go, "-v")) cfg->be_verbose = TRUE;        
  
  /* seed master's RNG */
  cfg->r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  /* create the stopwatch */
  cfg->w_stage = esl_stopwatch_Create();

  /* From HMMER 2.4X hmmcalibrate.c:
   * Generate calibrated CM(s) in a tmp file in the current
   * directory. When we're finished, we delete the original
   * CM file and rename() this one. That way, the worst
   * effect of a catastrophic failure should be that we
   * leave a tmp file lying around, but the original CM
   * file remains uncorrupted. tmpnam() doesn't work portably here,
   * because it'll put the file in /tmp and we won't
   * necessarily be able to rename() it from there.
   */
  ESL_ALLOC(cfg->tmpfile, (sizeof(char) * (strlen(cfg->cmfile) + 5)));
  strcpy(cfg->tmpfile, cfg->cmfile);
  strcat(cfg->tmpfile, ".xxx");	/* could be more inventive here... */
  if (esl_FileExists(cfg->tmpfile)) ESL_FAIL(eslFAIL, errbuf, "temporary file %s already exists; please delete it first", cfg->tmpfile);
  if (cfg->cmfp->is_binary) cfg->mode = "wb";
  else                      cfg->mode = "w"; 

  if(cfg->r        == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create master RNG.");
  if(cfg->w_stage  == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create stopwatch.");

  /* fill cfg->ccom, and cfg->cdate */
  if((status = get_cmcalibrate_comlog_info(go, cfg, errbuf)) != eslOK) return status;

  if(esl_opt_IsOn(go, "-T")) cfg->sc_cutoff = esl_opt_GetReal(go, "-T");
  else                       cfg->sc_cutoff = -eslINFINITY;

  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "init_master_cfg(), memory allocation error."); 
  return eslEMEM; /* NEVER REACHED */
}

/* serial_master()
 * The serial version of cmcalibrate.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;                /* Easel status */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm    = NULL;          /* the CM */
  CM_t    *nc_cm = NULL;          /* a non-configured copy of the CM */
  int      cmi;                   /* CM index, which model we're working on */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  uint32_t seed;                  /* RNG seed */
  int      exp_mode;              /* ctr over exp tail modes */
  double   cm_psec;               /* predicted number of seconds for calibrating current CM */
  double   cm_asec;               /* predicted number of seconds for calibrating current CM */
  double   total_psec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   total_asec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double  *exp_asecA;             /* stores actual timings for each exp tail fit stage, for each CM */
  double  *exp_psecA;             /* stores predicted timings for each exp tail fit stage, for each CM */
  double   psec;                  /* predicted seconds */
  /* exptail related vars */
  CM_TOPHITS       *th = NULL;            /* list of hits */
  int               exp_scN = 0;          /* number of hits reported thus far, for all seqs */
  float            *exp_scA = NULL;       /* [0..exp_scN-1] hit scores for all seqs */
  ESL_DSQ          *dsq = NULL;           /* digitized sequence to search */
  double           *dnull = NULL;         /* double version of cm->null, for generating random seqs */
  int               i;                    /* counter over sequences */
  int               h;                    /* counter over hits */
  double            tmp_mu, tmp_lambda;   /* temporary mu and lambda used for setting exp tails */
  int               tmp_nrandhits;        /* temporary number of rand hits found */
  float             tmp_tailp;            /* temporary tail mass probability fit to an exponential */
  
  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = print_run_info (go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  
  while ((status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm)) == eslOK) 
    {
      if      (status == eslEOD)  cm_Fail("read failed, CM file %s may be truncated?", cfg->cmfile);
      else if (status != eslOK)   cm_Fail(cfg->cmfp->errbuf);

      /* clone the non-configured CM we just read, we'll come back to it when we switch from global to local */
      if((status = cm_Clone(cm, errbuf, &nc_cm)) != eslOK) cm_Fail("unable to clone CM");

      if(cfg->ncm == 0) { /* first CM read, initialize our abc-dependent cfg values */
	if((status = CreateGenomicHMM(cm->abc, errbuf, &(cfg->ghmm_sA), &(cfg->ghmm_tAA), &(cfg->ghmm_eAA), &cfg->ghmm_nstates)) != eslOK) cm_Fail("unable to create generative HMM\n%s", errbuf);
	if((status = set_dnull(cfg, cm, errbuf)) != eslOK) cm_Fail("unable to create set_dnull\n%s\n", errbuf);
      }

      cfg->ncm++;
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->expAA, tmp, sizeof(ExpInfo_t **) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;

      if((status = initialize_cm(go, cfg, errbuf, cm, FALSE))          != eslOK) cm_Fail(errbuf);
      if((status = initialize_stats(go, cfg, errbuf))                  != eslOK) cm_Fail(errbuf);
      if((status = print_per_cm_column_headings(go, cfg, errbuf, cm))  != eslOK) cm_Fail(errbuf);
      
      /* allocate the exp_{a,p}secAA and fil_{a,p}secA arrays that hold {actual,predicted} times */
      ESL_ALLOC(exp_asecA, sizeof(double) * (EXP_NMODES));
      ESL_ALLOC(exp_psecA, sizeof(double) * (EXP_NMODES));
      esl_vec_DSet(exp_asecA, EXP_NMODES, 0.);
      esl_vec_DSet(exp_psecA, EXP_NMODES, 0.);
      cm_psec = cm_asec = 0.;

      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
	/* do we need to switch from global configuration to local? */
	if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	  /* switch from global to local by copying the current exptail stats from <cm>
	   * into <nc_cm> and then configure <nc_cm> for local mode. We do it this
	   * way because as a rule we don't allow reconfiguration of CMs (to limit
	   * execution paths through configuration functions)
	   */
	  FreeCM(cm);
	  cm = nc_cm;
	  if((status = initialize_cm(go, cfg, errbuf, cm, TRUE)) != eslOK) cm_Fail(errbuf);
	}
	/* update search info for round 0 (final round) for exp tail mode */
	if(ExpModeIsInside(exp_mode)) cm->search_opts |= CM_SEARCH_INSIDE;
	else                          cm->search_opts &= ~CM_SEARCH_INSIDE;

	/* We want to use the same seqs for exp tail fittings of all modes,
	 * so we free RNG, then create a new one and reseed it with the initial seed.
	 */
	seed = esl_randomness_GetSeed(cfg->r);
	esl_randomness_Destroy(cfg->r);
	cfg->r = esl_randomness_Create(seed);
	
	/************************************/
	/* exponential tail fitting section */
	/************************************/
	/* calculate exp tails for this exp mode */
	/* determine length of seqs to search for exp tail fitting */
	/* estimate time for this round */
	if((status = estimate_time_for_exp_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf); 
	psec *= cfg->expN * cfg->expL; /* psec was per residue */
	/* with --forecast, take into account parallelization */
	if((esl_opt_IsOn(go, "--forecast")) && (esl_opt_GetInteger(go, "--forecast") > 1)) psec /= (esl_opt_GetInteger(go, "--forecast") - 1);
	exp_psecA[exp_mode] = psec;
	cm_psec    += psec;
	total_psec += psec;
	print_exp_line(go, cfg, errbuf, exp_mode, cfg->expN, cfg->expL, psec);
	if(esl_opt_IsUsed(go, "--forecast")) continue; /* special mode, we don't do the calibration, just print the predicting timings */

	esl_stopwatch_Start(cfg->w_stage);
	fflush(stdout);
	
	ESL_DPRINTF1(("\n\ncalling process_exp_search_workunit to fit exp tail EXP mode: %d\n", exp_mode));
	
	exp_scN  = 0;
	for(i = 0; i < cfg->expN; i++) { 
	  /* do the work, fit the histogram, update exp tail info in expAA */
	  /* generate sequence, either randomly from background null or from hard-wired 5 state HMM that emits genome like sequence */
	  if(esl_opt_GetBoolean(go, "--random")) {
	    if((status = get_random_dsq(cfg, errbuf, cm, cfg->expL, &dsq)) != eslOK) cm_Fail(errbuf); 
	  }
	  else { 
	    if((status = SampleGenomicSequenceFromHMM(cfg->r, cm->abc, errbuf, cfg->ghmm_sA, cfg->ghmm_tAA, cfg->ghmm_eAA, cfg->ghmm_nstates, cfg->expL, &dsq)) != eslOK) cm_Fail(errbuf);
	  }
	  
	  /* to print seqs to stdout, uncomment this block */
	  /* ESL_SQ *mytmpdsq;
	     mytmpdsq = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, cfg->expL, NULL, NULL, NULL);
	     esl_sq_Textize(mytmpdsq);
	     printf(">seq%d\n%s\n", i, mytmpdsq->seq);
	     esl_sq_Destroy(mytmpdsq);
	     fflush(stdout);
	  */

	  if((status = process_exp_search_workunit(cm, errbuf, dsq, cfg->expL, cfg->sc_cutoff, &th)) != eslOK) cm_Fail(errbuf);
	  
	  printf("# mode: %12s\n", DescribeExpMode(exp_mode));
	  printf("Tophits dump before:\n");
	  cm_tophits_Dump(stdout, th);

	  cm_tophits_SortByPosition(th);
	  cm_tophits_RemoveOverlaps(th);
	  if(th->N > 0) { 
	    if(i == 0) ESL_ALLOC  (exp_scA, sizeof(float) * (exp_scN + th->N));
	    else       ESL_REALLOC(exp_scA, sizeof(float) * (exp_scN + th->N));
	    for(h = 0; h < th->N; h++) exp_scA[(exp_scN+h)] = th->hit[h]->score;
	    exp_scN += th->N;
	  }
	  /* TEMP  printf("after nresults: %8d\n", results->num_results); */
	  /* TEMP for(zz = 0; zz < results->num_results; zz++) 
	     printf("%5d  %5d  %10.3f\n", results->data[zz].start, results->data[zz].stop, results->data[zz].score);
	     fflush(stdout); */

	  printf("Tophits dump after:\n");
	  cm_tophits_Dump(stdout, th);
	  
	  cm_tophits_Destroy(th);
	  free(dsq);
	}
	if(cfg->exptfitfp != NULL) { 
	  fprintf(cfg->exptfitfp, "# CM: %s\n", cm->name);
	  fprintf(cfg->exptfitfp, "# mode: %12s\n", DescribeExpMode(exp_mode));
	}
	if((status = fit_histogram(go, cfg, errbuf, exp_scA, exp_scN, exp_mode, &tmp_mu, &tmp_lambda, &tmp_nrandhits, &tmp_tailp)) != eslOK) cm_Fail(errbuf);
	SetExpInfo(cfg->expAA[cmi][exp_mode], tmp_lambda, tmp_mu, (long) (cfg->expL * cfg->expN), tmp_nrandhits, tmp_tailp);
	
	esl_stopwatch_Stop(cfg->w_stage);
	exp_asecA[exp_mode] = cfg->w_stage->elapsed;
	cm_asec    += cfg->w_stage->elapsed;
	total_asec += cfg->w_stage->elapsed;
	FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
	printf("  %10s\n", time_buf);
	fflush(stdout);
	free(exp_scA);
      } /* end of for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) */

      if(cfg->be_verbose) if((status = debug_print_expinfo_array(cm, errbuf, cfg->expAA[cmi])) != eslOK) cm_Fail(errbuf);
      print_per_cm_summary(go, cfg, errbuf, cm, cm_psec, cm_asec);
      if(!esl_opt_IsOn(go, "--forecast")) { if((status = print_post_calibration_info(go, cfg, errbuf, stdout, cm, exp_psecA, exp_asecA)) != eslOK) cm_Fail(errbuf); }
      free(dnull);
      free(exp_asecA);
      free(exp_psecA);
      FreeCM(cm);

      printf("//\n");
      fflush(stdout);
    } /* end of while(cm_file_Read()) */
  if(status != eslEOF) cm_Fail(cfg->cmfp->errbuf);
  
  if(cfg->ncm > 1 && (esl_opt_IsOn(go, "--forecast"))) { 
    fprintf(stdout, "#\n");
    FormatTimeString(time_buf, total_psec, FALSE);
    fprintf(stdout, "# total predicted time for all %d CMs: %s\n", cfg->ncm, time_buf);
    fprintf(stdout, "#\n");
  }
  return;
      
 ERROR:
  cm_Fail("Memory allocation error.");
  return; /*NEVERREACHED*/
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmcalibrate
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * A master returns eslOK if it's successful.  Errors in an MPI master
 * come in two classes: recoverable and nonrecoverable.  If a
 * recoverable error occurs, errbuf is filled with an error message
 * from the master or a worker, and it's sent back while returning a
 * non-eslOK error code.
 * 
 * Recoverable errors include most worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers. The 
 * 
 * Some worker side errors (such as ESL_ALLOCs) are likely to be 
 * unrecoverable and will almost certainly cause MPI to crash
 * uncleanly, they're only here because I couldn't find a way around
 * them without massive reimplementation. Hopefully they rarely occur.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * cm_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int 
mpi_master(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int      xstatus       = eslOK; /* changes from OK on recoverable error */
  int      status;                /* Easel status */
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  int      n, i;                  /* counters */
  int      have_work     = TRUE;  /* TRUE while work remains  */
  int      nproc_working = 0;	  /* number of worker processes working, up to nproc-1 */
  int      wi;          	  /* rank of next worker to get a job to work on */
  int      wi_error = 0;          /* worker index that sent back an error message, if an error occurs */
  char    *buf           = NULL;  /* input/output buffer, for packed MPI messages */
  int      bn            = 0;     /* size of buf */
  int      pos = 1;               /* posn in buf */
  int      nseq_sent        = 0;  /* number of seqs we've told workers to work on */
  int      nseq_this_worker = 0;  /* number of seqs to tell current worker to work on */
  int      nseq_just_recv   = 0;  /* number of seqs we just received scores for from a worker */
  int      nseq_recv        = 0;  /* number of seqs we've received thus far this round from workers */
  int      seq_offset       = 0;  /* offset in master's list of sequences, rec'd then sent back, unchanged */
  MPI_Status mpistatus;           /* MPI status... */
  int      msg;                   /* holds integer telling workers we've finished current stage */
  int      exp_mode;              /* ctr over exp tail modes */
  int      h;                     /* ctr over hits */
  uint32_t seed;                  /* for seeding the master's RNG */
  /* variables for predicted and actual timings */
  double   psec;                  /* predicted number of seconds */
  double   cm_psec;               /* predicted number of seconds for calibrating current CM */
  double   cm_asec;               /* predicted number of seconds for calibrating current CM */
  double   total_asec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   total_psec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double  *exp_asecA;             /* stores actual timings for each exp tail fit stage, for each CM */
  double  *fil_asecA;             /* stores actual timings for each filter stage for each CM */
  double  *exp_psecA;             /* stores predicted timings for each exp tail fit stage, for each CM */
  double  *fil_psecA;             /* stores predicted timings for each filter stage for each CM */
  /* exponential tail related vars */
  int      exp_scN = 0;                                 /* number of hits reported thus far, for all seqs */
  float   *exp_scA = NULL;                              /* [0..exp_scN-1] hit scores for all seqs */
  int      exp_cm_cyk_mode;                             /* CYK    exp mode CM is in EXP_CM_LC or EXP_CM_GC */
  int      exp_cm_ins_mode;                             /* Inside exp mode CM is in EXP_CM_LI or EXP_CM_GI */
  int      si;                                          /* sequence index, the index of the last sequence generated */
  int      si_recv;                                     /* sequence index of the sequence we've just received results for from a worker */
  int      seqpos = 1;                                  /* sequence position in the current sequence */
  int      len;                                         /* length of a sequence chunk */
  int      chunksize;                                   /* size of chunks for each worker */
  double   tmp_mu, tmp_lambda;                          /* temporary mu and lambda used for setting exp tails */
  int      tmp_nrandhits;                               /* temporary number of rand hits found */
  float    tmp_tailp;                                   /* temporary tail mass probability fit to an exponential */
  int      need_seq;                                    /* TRUE if we are ready to generate a new seq */
  int      z;                                           /* counter */
  search_results_t *worker_results;                     /* results for seq si_recv we've just rec'd from a worker, we copy it to results_slist[si_recv] */
  /* *_slist variables: lists of data that are specific to each sequence 0..si..cfg->expN-1 */
  ESL_DSQ          **dsq_slist = NULL;                  /* [0..si..cfg->expN-1], the digitized sequences to search, when finished, they're freed */
  search_results_t **results_slist = NULL;              /* [0..si..cfg->expN-1], the compiled results from searching each seq si, when finished, copied to exp_scA and freed */
  int               *chunks_slist = NULL;               /* [0..si..cfg->expN-1], number of chunks of seq si currently being searched by workers */
  int               *sent_slist = NULL;                 /* [0..si..cfg->expN-1], TRUE if all chunks of seq si have been sent to workers, FALSE if not */
  /* *_wlist variables: lists of data that are specific to each worker 1..wi..nproc-1 */
  int *si_wlist = NULL;                                 /* [0..wi..nproc-1], the sequence index worker wi is working on */
  int *seqpos_wlist= NULL;                              /* [0..wi..nproc-1] the first position of the sequence that worker wi is searching */
  int *len_wlist = NULL;                                /* [0..wi..nproc-1] length of chunk worker wi is searching */

  /* filter threshold related vars */
  int       filN = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences to search for filter threshold calculation */
  int       fthr_mode = 0;         /* CM mode for filter threshold calculation, FTHR_CM_GC, FTHR_CM_GI, FTHR_CM_LC, FTHR_CM_LI */
  int       fil_cm_cyk_mode;       /* CYK    fthr mode CM is in FTHR_CM_LC or FTHR_CM_GC */
  int       fil_cm_ins_mode;       /* Inside fthr mode CM is in FTHR_CM_LI or FTHR_CM_GI */
  int       fil_nseq_per_worker  = (filN / (cfg->nproc-1)); /* when calcing filters, number of seqs to tell each worker to work on */
  /* full arrays of CYK, Inside, Fwd scores, [0..filN-1] */
  float    *fil_cyk_scA = NULL;    /* [0..filN-1] best cm cyk score for each emitted seq */
  float    *fil_ins_scA = NULL;    /* [0..filN-1] best cm insidei score for each emitted seq */
  float    *fil_fwd_scA = NULL;    /* [0..filN-1] best cp9 Forward score for each emitted seq */
  ESL_DSQ **fil_dsqA    = NULL;    /* [0..filN-1] CM emitted seq */
  int      *fil_L_A     = NULL;    /* [0..filN-1] length of emitted seq */
  /* worker's arrays of CYK, Inside, Fwd scores, [0..nseq_per_worker-1], rec'd from workers, copied to full arrays (ex: fil_cyk_scA) */
  float    *wkr_fil_cyk_scA = NULL;/* rec'd from worker: best cm cyk score for each emitted seq */
  float    *wkr_fil_ins_scA = NULL;/* rec'd from worker: best cm insidei score for each emitted seq */
  float    *wkr_fil_fwd_scA = NULL;/* rec'd from worker: best cp9 Forward score for each emitted seq */
  
  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { if ((status = print_run_info (go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn))             == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((si_wlist       = malloc(sizeof(int)  * cfg->nproc))     == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqpos_wlist   = malloc(sizeof(int)  * cfg->nproc))     == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((len_wlist      = malloc(sizeof(int)  * cfg->nproc))     == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* errbuf was filled above */
  ESL_DPRINTF1(("MPI master is initialized\n"));
  
  for (wi = 0; wi < cfg->nproc; wi++) {
    si_wlist[wi] = seqpos_wlist[wi] = len_wlist[wi] = -1;
  }
  
  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) cm_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));
  
  /* if we've used the --gc option, we read in a seq file 
   * to fill cfg->gc_freq, and we need to broadcast that info to workers
   */
  if( esl_opt_IsOn(go, "--gc")) { /* receive gc_freq info from master */
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

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
  
  while ((xstatus == eslOK) && ((status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm)) == eslOK)) 
    {
      cfg->ncm++;  
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->expAA, tmp, sizeof(ExpInfo_t **) * cfg->cmalloc);
	ESL_RALLOC(cfg->hfiAA, tmp, sizeof(HMMFilterInfo_t **) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;

      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm   (go, cfg, errbuf, cm))                                != eslOK) cm_Fail(errbuf);
      if((status = initialize_stats(go, cfg, errbuf, cm))                                != eslOK) cm_Fail(errbuf);
      if((status = print_per_cm_column_headings(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      
      /* allocate the exp_{a,p}secA and fil_{a,p}secA arrays that hold {actual,predicted} times */
      ESL_ALLOC(exp_asecA, sizeof(double) * (EXP_NMODES));
      ESL_ALLOC(fil_asecA, sizeof(double) * (EXP_NMODES));
      ESL_ALLOC(exp_psecA, sizeof(double) * (EXP_NMODES));
      ESL_ALLOC(fil_psecA, sizeof(double) * (EXP_NMODES));
      esl_vec_DSet(fil_asecA, EXP_NMODES, 0.);
      esl_vec_DSet(fil_psecA, EXP_NMODES, 0.);
      esl_vec_DSet(exp_asecA, EXP_NMODES, 0.);
      esl_vec_DSet(exp_psecA, EXP_NMODES, 0.);
      cm_psec = cm_asec = 0.;

      ESL_ALLOC(fil_cyk_scA, sizeof(float) * filN);
      ESL_ALLOC(fil_ins_scA, sizeof(float) * filN);
      ESL_ALLOC(fil_fwd_scA, sizeof(float) * filN);
      ESL_ALLOC(fil_L_A,     sizeof(int) *   filN);
      ESL_ALLOC(fil_dsqA,    sizeof(ESL_DSQ *) *   filN);
      
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
	/* allocate and initialize sequence lists */
	ESL_ALLOC(dsq_slist,     sizeof(ESL_DSQ *) * cfg->expN);
	ESL_ALLOC(results_slist, sizeof(search_results_t *) * cfg->expN);
	ESL_ALLOC(chunks_slist,  sizeof(int) * cfg->expN);
	ESL_ALLOC(sent_slist,    sizeof(int) * cfg->expN);
	for(z = 0; z < cfg->expN; z++) { 
	  dsq_slist[z]    = NULL; 
	  results_slist[z]= NULL;
	  chunks_slist[z] = 0;
	  sent_slist[z]   = FALSE;
	}
	
	/* do we need to switch from glocal configuration to local? */
	if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf))      != eslOK) cm_Fail(errbuf);
	}
	chunksize = (esl_opt_IsOn(go, "--Wchunk")) ? (cm->W * esl_opt_GetInteger(go, "--Wchunk")) : DetermineSeqChunksize(cfg->nproc, cfg->expL, cm->W);
	
	/* update search info for round 0 (final round) for exp tail mode */
	UpdateSearchInfoForExpMode(cm, 0, exp_mode);
	
	/* We want to use the same seqs for exp tail fittings of all CM modes and HMM modes, 
	 * so we free RNG, then create a new one and reseed it with the initial seed,
	 * The following pairs of modes will have identical sequences used for each member of the pair:
	 * 1. EXP_CP9_GV and EXP_CP9_GF
	 * 2. EXP_CM_GC  and EXP_CM_GI
	 * 3. EXP_CP9_LV and EXP_CP9_LF
	 * 4. EXP_CM_LC  and EXP_CM_LI
	 * Also the first min(--cmN-{loc,glc} <n>, --hmmN-{loc-glc} <n>) sequences between 1 and 2, and between 3 and 4,
	 * will also be identical.
	 */
	seed = esl_randomness_GetSeed(cfg->r);
	esl_randomness_Destroy(cfg->r);
	cfg->r = esl_randomness_Create(seed);
    
	/************************************/
	/* exponential tail fitting section */
	/************************************/
	/* fit exponential tails for this exp mode */
	ESL_DPRINTF1(("MPI master: CM: %d exp tail mode: %d\n", cfg->ncm, exp_mode));
	/* estimate time for this round, assuming all workers have same processor speed as master */
	if((status = estimate_time_for_exp_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf); 
	psec *= cfg->expN * cfg->expL; /* psec was per residue */
	if(cfg->nproc > 1) psec /= (cfg->nproc-1); /* parallelization will speed us up */
	exp_psecA[exp_mode] = psec;
	cm_psec    += psec;
	total_psec += psec;
	print_exp_line(go, cfg, errbuf, exp_mode, cfg->expN, cfg->expL, psec);

	esl_stopwatch_Start(cfg->w_stage);
	exp_scN  = 0;
	
	if(xstatus == eslOK) have_work     = TRUE;	/* TRUE while work remains  */
	
	wi = 1;
	si   = -1;
	need_seq = TRUE;
	have_work = TRUE;	/* TRUE while work remains  */
	seqpos = 1;
	
	while (have_work || nproc_working)
	  {
	    if (need_seq) 
	      {
		need_seq = FALSE;
		/* generate a new seq */
		si++;
		if(si < cfg->expN)
		  {
		    /* generate sequence, either randomly from background null or from hard-wired 5 state HMM that emits genome like sequence */
		    if(esl_opt_GetBoolean(go, "--random")) { 
		      if((status = get_random_dsq(cfg, errbuf, cm, cfg->expL, &(dsq_slist[si]))) != eslOK) goto ERROR; 
		    }
		    else { 
		      if((status = SampleGenomicSequenceFromHMM(cfg->r, cm->abc, errbuf, cfg->ghmm_sA, cfg->ghmm_tAA, cfg->ghmm_eAA, cfg->ghmm_nstates, cfg->expL, &(dsq_slist[si]))) != eslOK) goto ERROR;
		    }
		    /* TEMP */
		    /* ESL_SQ *tmp;
		       tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq_slist[si], cfg->expL, NULL, NULL, NULL);
		       esl_sq_Textize(tmp);
		       printf(">seq%d\n%s\n", si, tmp->seq);
		       esl_sq_Destroy(tmp);
		       fflush(stdout); */
		    /* TEMP */
		    
		    results_slist[si] = CreateResults(INIT_RESULTS);
		    sent_slist[si]    = FALSE;
		    chunks_slist[si]  = 0;
		    seqpos = 1;
		    have_work = TRUE;
		  }
		else if(si == cfg->expN) have_work = FALSE; 
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
			si_recv = si_wlist[wi];
			ESL_DPRINTF1(("MPI master sees that the result buffer contains search results\n"));
			if ((status = cm_search_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_results)) != eslOK) cm_Fail("search results unpack failed");
			ESL_DPRINTF1(("MPI master has unpacked search results\n"));
			
			if(worker_results != NULL) { 
			  /* add results to seqlist[si_recv]->results[rclist[wi]] */
			  AppendResults(worker_results, results_slist[si_recv], seqpos_wlist[wi]);
			  /* careful, dbseqlist[si_recv]->results[rclist[wi]] now points to the traces and postal codes in worker_results->data,
			   * don't free those (don't use FreeResults(worker_results)) */
			  free(worker_results->data);
			  free(worker_results);
			  worker_results = NULL;
			}
			chunks_slist[si_recv]--;
			if(sent_slist[si_recv] && chunks_slist[si_recv] == 0) 
			  { /* we're done with sequence si_recv; remove overlapping hits, copy scores of remaining hits, then free data */
			    if(results_slist[si_recv]->num_results > 0) 
			      { 
				/* TEMP printf("# mode: %12s\n", DescribeExpMode(exp_mode)); */
				/* TEMP printf("si_recv: %4d expL: %d before nresults: %8d\n", si_recv, cfg->expL, results_slist[si_recv]->num_results); */
				/* TEMP int zz;
				   for(zz = 0; zz < results_slist[si_recv]->num_results; zz++) 
				   printf("%5d  %5d  %10.3f\n", results_slist[si_recv]->data[zz].start, results_slist[si_recv]->data[zz].stop, results_slist[si_recv]->data[zz].score);
				*/
				RemoveOverlappingHits(results_slist[si_recv], 1, cfg->expL);
				if(exp_scA == NULL) ESL_ALLOC (exp_scA,      sizeof(float) * (exp_scN + results_slist[si_recv]->num_results));
				else                ESL_RALLOC(exp_scA, tmp, sizeof(float) * (exp_scN + results_slist[si_recv]->num_results));
				for(h = 0; h < results_slist[si_recv]->num_results; h++) exp_scA[(exp_scN+h)] = results_slist[si_recv]->data[h].score;
				exp_scN += results_slist[si_recv]->num_results;
				/* TEMP printf("after nresults: %8d\n", results_slist[si_recv]->num_results); */
				/* TEMP for(zz = 0; zz < results_slist[si_recv]->num_results; zz++) 
				   printf("%5d  %5d  %10.3f\n", results_slist[si_recv]->data[zz].start, results_slist[si_recv]->data[zz].stop, results_slist[si_recv]->data[zz].score);
				   fflush(stdout); */
			      }
			    free(dsq_slist[si_recv]);
			    FreeResults(results_slist[si_recv]);
			    dsq_slist[si_recv] = NULL;
			    results_slist[si_recv] = NULL;
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
		len = (chunksize < (cfg->expL - seqpos + 1)) ? chunksize : (cfg->expL - seqpos + 1);
		ESL_DPRINTF1(("MPI master is sending sequence i0..j0 %d..%d to search to worker %d\n", seqpos, seqpos+len-1, wi));
		assert(seqpos > 0);
		
		if ((status = cm_dsq_MPISend(dsq_slist[si]+seqpos-1, len, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
		
		si_wlist[wi]     = si;
		seqpos_wlist[wi] = seqpos;
		len_wlist[wi]    = len;
		chunks_slist[si]++;
		
		/* TEMP printf("sending chunk of seq: %4d of len: %5d to wi: %4d seqpos: %4d\n", si, len, wi, seqpos); */
		
		wi++;
		nproc_working++;
		
		if(len == chunksize) seqpos += len - cm->W + 1;
		else {
		  need_seq       = TRUE;
		  sent_slist[si] = TRUE; /* we've sent all chunks from this seq */
		}
	      }
	  }
	ESL_DPRINTF1(("MPI master: done with this exp tail mode. Telling all workers\n"));
	for (wi = 1; wi < cfg->nproc; wi++) 
	  if ((status = cm_dsq_MPISend(NULL, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");
	
	/* fill a histogram with the exp_scN scores in exp_scA and fit it to an exponential tail */
	if(cfg->exptfitfp != NULL) { 
	  fprintf(cfg->exptfitfp, "# CM: %s\n", cm->name);
	  fprintf(cfg->exptfitfp, "# mode: %12s\n", DescribeExpMode(exp_mode));
	}
	if((status = fit_histogram(go, cfg, errbuf, exp_scA, exp_scN, exp_mode, &tmp_mu, &tmp_lambda, &tmp_nrandhits, &tmp_tailp)) != eslOK) cm_Fail(errbuf);
	SetExpInfo(cfg->expAA[cmi][exp_mode], tmp_lambda, tmp_mu, (long) (cfg->expL * cfg->expN), tmp_nrandhits, tmp_tailp);
	
	for(si = 0; si < cfg->expN; si++) {
	  ESL_DASSERT1((dsq_slist[si] == NULL));
	  ESL_DASSERT1((results_slist[si] == NULL));
	  ESL_DASSERT1((chunks_slist[si] == 0));
	  ESL_DASSERT1((sent_slist[si] == TRUE));
	}
	esl_stopwatch_Stop(cfg->w_stage);
	exp_asecA[exp_mode] = cfg->w_stage->elapsed;
	cm_asec += cfg->w_stage->elapsed;
	total_asec += cfg->w_stage->elapsed;
	FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
	printf("  %10s\n", time_buf);
	fflush(stdout);
	
	free(exp_scA); 
	exp_scA = NULL;
	free(dsq_slist);
	free(results_slist);
	free(chunks_slist);
	free(sent_slist);

	/****************************/
	/* filter threshold section */
	/****************************/
	if(exp_mode == EXP_CM_GI || exp_mode == EXP_CM_LI) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  if((status = estimate_time_for_fil_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf);
	  psec *= filN;
	  if(cfg->nproc > 1) psec /= (cfg->nproc-1); /* parallelization will speed us up */
	  fil_psecA[exp_mode] = psec;
	  cm_psec    += psec;
	  total_psec += psec;
	  print_fil_line(go, cfg, errbuf, exp_mode, psec);
	  
	  esl_stopwatch_Start(cfg->w_stage);
	  fthr_mode = ExpModeToFthrMode(exp_mode);
	  ESL_DPRINTF1(("MPI master: CM: %d fthr mode: %d\n", cfg->ncm, fthr_mode));
	  
	  /* generate all filN sequences before sending any */
	  ESL_ALLOC(fil_dsqA,  sizeof(ESL_DSQ *) * filN);
	  ESL_ALLOC(fil_L_A,   sizeof(int) * filN);
	  for(i = 0; i < filN; i++) { 
	    if((status = get_cmemit_dsq(cfg, errbuf, cm, &fil_L_A[i], &fil_dsqA[i])) != eslOK) cm_Fail(errbuf);
	  }
	  for(i = 0; i < filN; i++) { 
	    if((status = get_cmemit_dsq(cfg, errbuf, cm, &fil_L_A[i], &fil_dsqA[i])) != eslOK) cm_Fail(errbuf);
	    /* TEMP to print seqs to stdout uncomment this block */
	    /* TEMP ESL_SQ *tmp;
	       tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", fil_dsqA[i], fil_L_A[i], NULL, NULL, NULL);
	       esl_sq_Textize(tmp);
	       printf(">seq%d\n%s\n", i, tmp->seq);
	       esl_sq_Destroy(tmp);
	       fflush(stdout); */
	  }
	  
	  if(xstatus == eslOK) have_work = TRUE;  /* TRUE while work remains  */
	  else                 have_work = FALSE; /* we've seen an error and are trying to finish cleanly */
	  
	  wi = 1;
	  
	  nseq_sent = 0;
	  nseq_recv = 0;
	  while (have_work || nproc_working)
	    {
	      if(have_work) { 
		if(nseq_sent < filN) {
		  nseq_this_worker = ((nseq_sent + fil_nseq_per_worker) <= filN) ? 
		    fil_nseq_per_worker : (filN - nseq_sent);
		}
		else { 
		  have_work = FALSE;
		  ESL_DPRINTF1(("MPI master has run out of numbers of sequences to dole out (%d doled)\n", nseq_sent));
		}
	      }
	      if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0)) {
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
			ESL_DPRINTF1(("MPI master sees that the result buffer contains HMM filter results\n"));
			if ((status = cmcalibrate_filter_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &wkr_fil_cyk_scA, &wkr_fil_ins_scA, &wkr_fil_fwd_scA, &nseq_just_recv, &seq_offset)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
			ESL_DPRINTF1(("MPI master has unpacked HMM filter results\n"));
			for(i = 0; i < nseq_just_recv; i++) {
			  fil_cyk_scA[seq_offset+i] = wkr_fil_cyk_scA[i];
			  fil_ins_scA[seq_offset+i] = wkr_fil_ins_scA[i];
			  fil_fwd_scA[seq_offset+i] = wkr_fil_fwd_scA[i];
			}
			free(wkr_fil_cyk_scA);
			free(wkr_fil_ins_scA);
			free(wkr_fil_fwd_scA);
			nseq_recv += nseq_just_recv;
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
		  /* send number of sequences, then seq offset, then sequences */
		  MPI_Send(&(nseq_this_worker), 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
		  MPI_Send(&(nseq_sent),        1, MPI_INT, wi, 0, MPI_COMM_WORLD);
		  for(i = nseq_sent; i < (nseq_sent + nseq_this_worker); i++) { 
		    ESL_DPRINTF1(("MPI master is sending HMM filter seq %d to worker %d\n", i, wi));
		    if ((status = cm_dsq_MPISend(fil_dsqA[i], fil_L_A[i], wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI filter job send failed");
		    free(fil_dsqA[i]);
		  }
		  
		  wi++;
		  nproc_working++;
		  nseq_sent += nseq_this_worker;
		}
	    }
	  
	  if(xstatus == eslOK) { 
	    exp_cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? EXP_CM_LC  : EXP_CM_GC;
	    exp_cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? EXP_CM_LI  : EXP_CM_GI;
	    fil_cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? FTHR_CM_LC : FTHR_CM_GC;
	    fil_cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? FTHR_CM_LI : FTHR_CM_GI;
	    /* set cutoffs for forward HMM filters, first for CYK, then for Inside */
	    if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_cyk_scA, fil_fwd_scA, exp_cm_cyk_mode, cfg->hfiAA[cmi][fil_cm_cyk_mode])) != eslOK) cm_Fail(errbuf);
	    if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_ins_scA, fil_fwd_scA, exp_cm_ins_mode, cfg->hfiAA[cmi][fil_cm_ins_mode])) != eslOK) cm_Fail(errbuf);
	  }
	  ESL_DPRINTF1(("MPI master: done with HMM filter calc for fthr mode %d for this CM.\n", fthr_mode));
	  
	  for (wi = 1; wi < cfg->nproc; wi++) { 
	    msg = MPI_FINISHED_FILTER;
	    MPI_Send(&msg, 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	  }
	  esl_stopwatch_Stop(cfg->w_stage);
	  fil_asecA[exp_mode] = cfg->w_stage->elapsed;
	  cm_asec +=  cfg->w_stage->elapsed;
	  total_asec += cfg->w_stage->elapsed;
	  FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
	  printf("  %10s\n", time_buf);
	  fflush(stdout);
	}
	ESL_DPRINTF1(("MPI master: done with exp tail mode %d for this CM.\n", exp_mode));
      }
      ESL_DPRINTF1(("MPI master: done with this CM.\n"));
      if(xstatus == eslOK) if(cfg->be_verbose) if((status = debug_print_expinfo_and_filterinfo_arrays(cm, errbuf, cfg->expAA[cmi], cfg->hfiAA[cmi])) != eslOK) cm_Fail(errbuf);
      print_per_cm_summary(go, cfg, errbuf, cm, cm_psec, cm_asec);
      if(! esl_opt_IsOn(go, "--forecast")) { if((status = print_post_calibration_info(go, cfg, errbuf, stdout, cm, exp_psecA, fil_psecA, exp_asecA, fil_asecA)) != eslOK) cm_Fail(errbuf); }
      free(fil_cyk_scA);
      free(fil_ins_scA);
      free(fil_fwd_scA);
      free(fil_dsqA);
      free(fil_L_A);
      free(exp_asecA);
      free(exp_psecA);
      free(fil_asecA); 
      free(fil_psecA); 
      FreeCM(cm);
      
      printf("//\n");
      fflush(stdout);
    }
  free(si_wlist);
  free(seqpos_wlist);
  free(len_wlist);
  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if     (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); return xstatus; }
  else if(status != eslEOF) { strcpy(errbuf, cfg->cmfp->errbuf); return status; }
  else                      return eslOK;
  
 ERROR: 
  ESL_FAIL(status, errbuf, "mpi_master() memory allocation error.");
  return eslOK; /* NOTREACHED */
}

static int
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;                /* Easel status */
  int      xstatus = eslOK;       /* changes from OK on recoverable error */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  int      i;                     /* counter over seqs */
  char    *wbuf = NULL;	          /* packed send/recv buffer  */
  int      wn   = 0;	          /* allocation size for wbuf */
  int      sz, n;		  /* size of a packed message */
  int      pos;                   /* posn in wbuf */
  int      nseq;                  /* number of seqs to emit/search for current job */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  MPI_Status mpistatus;           /* MPI status... */

  /* exponential tail related vars */
  ESL_DSQ          *exp_dsq = NULL;     /* dsq chunk received from master */
  int               expL;               /* length of exp_dsq received from master */
  search_results_t *exp_results = NULL; /* hits found in exp_dsq, sent back to master */
  int               exp_mode  = 0;

  /* filter threshold related vars */
  int       fthr_mode = 0;         /* CM mode for filter threshold calculation, FTHR_CM_GC, FTHR_CM_GI, FTHR_CM_LC, FTHR_CM_LI */
  float    *fil_cyk_scA = NULL;    /* [0..nseq-1] best cm cyk score for each emitted seq */
  float    *fil_ins_scA = NULL;    /* [0..nseq-1] best cm insidei score for each emitted seq */
  float    *fil_fwd_scA = NULL;    /* [0..nseq-1] best cp9 Forward score for each emitted seq */
  ESL_DSQ **fil_dsqA    = NULL;    /* [0..nseq-1] CM emitted seqs, rec'd from master */
  int      *fil_L_A     = NULL;    /* [0..nseq-1] lengths of seqs, rec'd from master */
  int       seq_offset;            /* sent by master, we send back so master knows where to store worker's results */
  int       in_fil_section_flag = FALSE;/* set to TRUE while we're in the filter threshold calculation
					 * section, we need to know this when we goto ERROR, b/c we have
					 * to know how many more MPI_Recv() calls to make to match up
					 * with the Master's sends before we can shut down.
					 */

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* master saw an error code; workers do an immediate normal shutdown. */
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
    return xstatus; /* shutdown; we passed the error back for the master to deal with. */
  }

  /* if we've used the --gc option, we read in a seq file to fill
   * cfg->gc_freq, and we need that info here for the worker, so we receive
   * it's broadcast from the master
   */
  if( esl_opt_IsOn(go, "--gc")) { /* receive gc_freq info from master */
    ESL_DASSERT1((cfg->gc_freq == NULL));
    ESL_ALLOC(cfg->gc_freq,  sizeof(double) * GC_SEGMENTS);
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else cfg->gc_freq = NULL; /* default */
  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      cfg->ncm++;  
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->expAA, tmp, sizeof(ExpInfo_t **) * cfg->cmalloc);
	ESL_RALLOC(cfg->hfiAA, tmp, sizeof(HMMFilterInfo_t **) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm(go, cfg, errbuf, cm))                                   != eslOK) goto ERROR;
      if((status = initialize_stats(go, cfg, errbuf, cm))                                != eslOK) goto ERROR;
      
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {

	/* do we need to switch from glocal configuration to local? */
	if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf))      != eslOK) goto ERROR;
	}
	/* update search info for round 0 (final round) for exp tail mode */
	UpdateSearchInfoForExpMode(cm, 0, exp_mode);

	/************************************/
	/* exponential tail fitting section */
	/************************************/
	ESL_DPRINTF1(("worker %d exp_mode: %d\n", cfg->my_rank, exp_mode));
	
	while((status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &exp_dsq, &expL)) == eslOK)
	  {
	    ESL_DPRINTF1(("worker %d: has received dsq chunk of length L: %d\n", cfg->my_rank, expL));
	    if ((status = ProcessSearchWorkunit(cm, errbuf, exp_dsq, expL, 
						ESL_MAX(1, ((int) cfg->avg_hit_len)), /* guess at average hit len for HMM scanning functions */
						&exp_results, esl_opt_GetReal(go, "--mxsize"), cfg->my_rank)) != eslOK) goto ERROR;
	    
	    /*RemoveOverlappingHits(exp_results, 1, expL);*/
	    ESL_DPRINTF1(("worker %d: has gathered search results\n", cfg->my_rank));
	    
	    n = 0;
	    if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	      ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	    n += sz;
	    if (cm_search_results_MPIPackSize(exp_results, MPI_COMM_WORLD, &sz) != eslOK)
	      ESL_XFAIL(eslFAIL, errbuf, "cm_search_results_MPIPackSize() call failed"); 
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
	    if (cm_search_results_MPIPack(exp_results, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
	      ESL_XFAIL(eslFAIL, errbuf, "cm_search_results_MPIPack() call failed"); 
	    MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	    ESL_DPRINTF1(("worker %d: has sent results (%d hits) to master in message of %d bytes\n", cfg->my_rank, exp_results->num_results, pos));
	    /* TEMP printf("worker %d: has sent results (%d hits) to master in message of %d bytes\n", cfg->my_rank, exp_results->num_results, pos); */
	    
	    FreeResults(exp_results);
	    free(exp_dsq);
	  }
	ESL_DPRINTF1(("worker %d finished exp_mode: %d\n", cfg->my_rank, exp_mode));

	
	/****************************/
	/* filter threshold section */
	/****************************/
	if(exp_mode == EXP_CM_GI || exp_mode == EXP_CM_LI) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  in_fil_section_flag = TRUE;
	  fthr_mode = ExpModeToFthrMode(exp_mode);
	  ESL_DPRINTF1(("worker %d fthr_mode: %d\n", cfg->my_rank, fthr_mode));

	  if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	  while(nseq != MPI_FINISHED_FILTER) {
	    ESL_DPRINTF1(("worker %d: has received hmm filter nseq: %d\n", cfg->my_rank, nseq));
	    
	    if(MPI_Recv(&seq_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");

	    /* allocate for the sequences */
	    ESL_ALLOC(fil_dsqA,  sizeof(ESL_DSQ *) * nseq);
	    ESL_ALLOC(fil_L_A,   sizeof(int) * nseq);
	    for(i = 0; i < nseq; i++) { 
	      if((status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &fil_dsqA[i], &fil_L_A[i])) != eslOK) return status;
	    }
	      
	    if((status = process_filter_workunit (go, cfg, errbuf, cm, fil_dsqA, fil_L_A, nseq, &fil_cyk_scA, &fil_ins_scA, &fil_fwd_scA)) != eslOK) return status;
	    ESL_DPRINTF1(("worker %d: has gathered HMM filter results\n", cfg->my_rank));
	    n = 0;

	    if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	      ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	    n += sz;
	    if(cmcalibrate_filter_results_MPIPackSize(nseq, MPI_COMM_WORLD, &sz) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPackSize() call failed"); 
	    n += sz;  
	    if (n > wn) {
	      void *tmp;
	      ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	      wn = n;
	    }
	    ESL_DPRINTF1(("worker %d: has calculated the HMM filter results will pack into %d bytes\n", cfg->my_rank, n));
	    status = eslOK;
	    pos = 0;

	    if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	      ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	    if (cmcalibrate_filter_results_MPIPack(fil_cyk_scA, fil_ins_scA, fil_fwd_scA, nseq, seq_offset, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPack() call failed"); 
	    free(fil_cyk_scA);
	    free(fil_ins_scA);
	    free(fil_fwd_scA);
	    for(i = 0; i < nseq; i++) free(fil_dsqA[i]);
	    free(fil_dsqA);
	    free(fil_L_A);

	    MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	    ESL_DPRINTF1(("worker %d: has sent CP9 filter results to master in message of %d bytes\n", cfg->my_rank, pos));
	    /* receive next number of sequences, if MPI_FINISHED_EXPTAIL, we'll stop */
	    if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	  }
	  in_fil_section_flag = FALSE;
	}
      } /* end of for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) */

      FreeCM(cm);
      cm = NULL;
      ESL_DPRINTF1(("worker %d finished all exp_modes for this cm.\n", cfg->my_rank));
    }
  if (status == eslEOD) ESL_DPRINTF1(("Worker %d told CMs are done.\n", cfg->my_rank));
  else goto ERROR;
  
  if (wbuf != NULL) free(wbuf);
  return eslOK;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  if(status == eslEMEM) sprintf(errbuf, "Memory allocation error.");
  MPI_Pack(&status, 1,               MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

  /* if we get here this worker failed and sent an error message, now the master knows a worker
   * failed but it has to continue through the mpi_master() code, sending the messages that
   * the workers expect, telling them to continue to move through the loops in those functions.
   * Minimal work will be done, but this is necessary so that we shut down cleanly. 
   * Because the master is sending messages to us still, we have to receive them. We can't
   * check that they're the expected messages though (codes telling us to keep moving through
   * the loops) because even if they were the wrong messages we couldn't do anything about it,
   * we've already entered error mode.
   */
  if(in_fil_section_flag) MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  for(; exp_mode < EXP_NMODES; exp_mode++) {
    MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
    if(ExpModeIsForCM(exp_mode)) {
      MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
    }
  }
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);

  return eslFAIL; /* recoverable error, master has error message and will print it */
}
#endif /*HAVE_MPI*/

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int do_local)
{
  int status;

  /* config QDB? */
  if(esl_opt_GetBoolean(go, "--no-qdb")) { 
    cm->search_opts |= CM_SEARCH_NONBANDED; /* don't use QDB to search */
    /* no need to recalculate QDBs, don't raise CM_CONFIG_QDB */
  }
  else {
    cm->search_opts |= CM_SEARCH_QDB; /* use QDB to search */
    /* check if cm->qdbinfo->beta2 == <x> from --beta, if so we don't need to recalculate QDBs */
    if(CheckCMQDBInfo(cm->qdbinfo, 0., FALSE, esl_opt_GetReal(go, "--beta"), TRUE) != eslOK) {
      /* we'll use beta2 for calibration, setting them both as equal makes it slightly more efficient */
      cm->config_opts |= CM_CONFIG_QDB;   /* configure QDB */
      cm->qdbinfo->beta1 = esl_opt_GetReal(go, "--beta"); 
      cm->qdbinfo->beta2 = esl_opt_GetReal(go, "--beta"); 
    }
  }

  cm->search_opts |= CM_SEARCH_NOALIGN;

  if(! esl_opt_GetBoolean(go, "--no-null3")) cm->search_opts |= CM_SEARCH_NULL3;

  /* ALWAYS use the greedy overlap resolution algorithm to return hits for exp calculation
   * it's irrelevant for filter threshold stats, we return best score per seq for that */
  /* don't turn on CM_SEARCH_CMNOTGREEDY */

  if(do_local) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  /* we'll need a scan matrix too */
  cm->config_opts |= CM_CONFIG_SCANMX;

  /* configure */
  if((status = cm_Configure(cm, errbuf)) != eslOK) return status; 

  if(cm->smx == NULL) ESL_FAIL(eslEINVAL, errbuf, "unable to create scan matrix for CM");
  return eslOK;
}

/* initialize_stats()
 * Allocate and initialize cfg->expAA */
static int
initialize_stats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int i;
  int cmi = cfg->ncm-1;

  ESL_DPRINTF1(("initializing cmstats\n"));

  ESL_ALLOC(cfg->expAA[cmi], sizeof(ExpInfo_t *) * EXP_NMODES);
  for(i = 0; i < EXP_NMODES; i++) cfg->expAA[cmi][i] = CreateExpInfo();
  
  return eslOK;

 ERROR: 
  return status;
}


/* fit_histogram()
 * Create, fill and fit the tail of a histogram to an exponential tail. Data to fill the histogram
 * is given as <scores>.
 */
static int
fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, int exp_mode, double *ret_mu, double *ret_lambda, int *ret_nrandhits, float *ret_tailp)
{
  int status;
  double mu;
  double lambda;
  int i;
  double *xv;         /* raw data from histogram */
  int     n,z;  
  float tailp;
  double  params[2];
  int     nrandhits; 
  float   a;
  float   nhits_to_fit;

  ESL_HISTOGRAM *h = NULL;       /* histogram of scores */

  /* Initialize histogram; these numbers are guesses */
  if((h = esl_histogram_CreateFull(-100., 100., .1)) == NULL) return eslEMEM;    

  /* fill histogram */
  for(i = 0; i < nscores; i++) {
    if((status = esl_histogram_Add(h, scores[i])) != eslOK) ESL_FAIL(status, errbuf, "fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);
    /* printf("%4d %.3f\n", i, scores[i]); */
  }

  /* fit scores to an exponential tail */
  if(cfg->exptfitfp != NULL) { 
    /* fit to 41 different tailp values and print lambda, mu to a save file*/
    fprintf(cfg->exptfitfp, "# %11s  %10s  %10s  %12s\n", "tail pmass",  "lambda",     "mu",         "nhits");
    fprintf(cfg->exptfitfp, "# %11s  %10s  %10s  %12s\n", "-----------", "----------", "----------", "------------");
    for(a = 0.; a >= -4.; a -= 0.1) { 
      tailp = pow(10., a);
      esl_histogram_GetTailByMass(h, tailp, &xv, &n, &z); 
      if(n > 1) { 
	esl_exp_FitComplete(xv, n, &(params[0]), &(params[1]));
	esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params);
	fprintf(cfg->exptfitfp, "  %.9f  %10.6f  %10.4f  %12d\n", tailp, params[1], params[0], n);
      }
      else { 
	fprintf(cfg->exptfitfp, "  %.9f  %10s  %10s  %12d\n", tailp, "N/A", "N/A", n);
      }
    }
    fprintf(cfg->exptfitfp, "//\n");
  }
  /* end of if cfg->exptfitfp != NULL) */

  /* determine the fraction of the tail to fit, if --tail-p, it's easy */
  if(esl_opt_IsOn(go, "--tailp")) { 
    tailp = esl_opt_GetReal(go, "--tailp");
    tailp = ESL_MIN(tailp, ((float) esl_opt_GetInteger(go, "--tailxn") / (float) h->n)); /* ensure we don't exceed our max nhits in tail */
  }
  else { /* number of hits is per Mb and specific to local or glocal fits */
    if(ExpModeIsLocal(exp_mode)) { /* local mode */
      nhits_to_fit = (float) esl_opt_GetInteger(go, "--ltailn") * ((cfg->expN * cfg->expL) / 1000000.);
      tailp = nhits_to_fit / (float) h->n;
      if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--ltailn <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --tailp.", esl_opt_GetInteger(go, "--ltailn"), (h->n / ((float) cfg->expN * ((float) cfg->expL) / 1000000.)));
    }
    else { /* glocal mode */
      nhits_to_fit = (float) esl_opt_GetInteger(go, "--gtailn") * ((cfg->expN * cfg->expL) / 1000000.);
      tailp = nhits_to_fit / (float) h->n;
      if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--gtailn <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --tailp.", esl_opt_GetInteger(go, "--gtailn"), (h->n / ((float) cfg->expN * ((float) cfg->expL) / 1000000.)));
    }
  }

  esl_histogram_GetTailByMass(h, tailp, &xv, &n, &z); /* fit to right 'tailfit' fraction, 0.01 by default */
  if(n <= 1) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), too few points in right tailfit: %f fraction of histogram. Increase <x> with -L <x>.", tailp);
  esl_exp_FitComplete(xv, n, &(params[0]), &(params[1]));
  esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params);

  /* printf("# Exponential fit to %.7f%% tail: lambda = %f\n", tailp*100.0, params[1]); */
  mu = params[0];
  lambda = params[1];
  if(isnan(lambda)) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), exp tail fit lambda is NaN, too few hits in histogram. Increase <x> with -L <x>.");
  if(isinf(lambda)) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), exp tail fit lambda is inf, too few hits in histogram. Increase <x> with -L <x>.");
  nrandhits = h->n; /* total number of hits in the histogram */

  /* print to output files if nec */
  if(cfg->exphfp != NULL) esl_histogram_Plot(cfg->exphfp, h);
  if(cfg->expqfp != NULL) esl_histogram_PlotQQ(cfg->expqfp, h, &esl_exp_generic_invcdf, params);
  if (cfg->expsfp != NULL) {
    esl_histogram_PlotSurvival(cfg->expsfp, h);
    esl_exp_Plot(cfg->expsfp, (params[0] - log(1./tailp) / params[1]), 0.693147, esl_exp_surv, h->xmin - 5., h->xmax + 5., 0.1); /* extrapolate mu */
  }
  esl_histogram_Destroy(h);

  *ret_mu     = mu;
  *ret_lambda = lambda;
  *ret_nrandhits = nrandhits;
  *ret_tailp = tailp;
  return eslOK;
}

/* Function: get_random_dsq()
 * Date:     EPN, Tue Sep 11 08:31:47 2007
 * 
 * Purpose:  Generate a random digitized seq and return it.
 *           Two possible modes:
 *           1. if(cfg->gc_freq == NULL && dnull != NULL) 
 *              use dnull disto (a double version of cm->null) to generate
 *           2. if(cfg->gc_freq != NULL && dnull == NULL) 
 *              use choose a GC frequency from cfg->gc_freq
 *              and generate with that
 *
 * Returns:  eslOK on success, ret_dsq filled with newly alloc'ed ESL_DSQ *,
 *           some other status code on failure.
 */
int
get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int L, ESL_DSQ **ret_dsq)
{
  int      status;
  double   gc_comp;
  double  *distro = NULL;
  int      do_free_distro = FALSE;
  ESL_DSQ *dsq = NULL;

  /* contract check, make sure we're in a valid mode */
  if(cfg->gc_freq == NULL && cfg->dnull == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->gc_freq == NULL and cfg->dnull == NULL");
  if(cfg->gc_freq != NULL && cfg->dnull != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->gc_freq != NULL and cfg->dnull != NULL");

  /* determine mode */ /* generate sequence */
  if      (cfg->gc_freq == NULL && cfg->dnull != NULL) distro = cfg->dnull;
  else if (cfg->gc_freq != NULL && cfg->dnull == NULL) {
    assert(cm->abc->K == 4);
    ESL_ALLOC(distro, sizeof(double) * cm->abc->K);
    do_free_distro = TRUE;
    gc_comp = 0.01 * esl_rnd_DChoose(cfg->r, cfg->gc_freq, GC_SEGMENTS);
    distro[1] = distro[2] = 0.5 * gc_comp;
    distro[0] = distro[3] = 0.5 * (1. - gc_comp);
  }
  /* generate sequence */
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if ((status = esl_rsq_xIID(cfg->r, distro, cm->abc->K, L, dsq) != eslOK)) return status;

  if (do_free_distro) free(distro);
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  return status;
}

/* Function: get_cmcalibrate_comlog_info
 * Date:     EPN, Mon Dec 31 14:59:52 2007
 *
 * Purpose:  Create the cmcalibrate command info and creation date info 
 *           to eventually be set in the CM's ComLog_t data structure.
 *
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation.
 */
static int
get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int i;
  uint32_t seed;
  int  temp;
  int  seedlen;
  char *seedstr;
  time_t date = time(NULL);

  if(cfg->ccom  != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->ccom  is non-NULL.");
  if(cfg->cdate != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->cdate is non-NULL.");
  
  /* Set the cmcalibrate command info, the cfg->ccom string */
  for (i = 0; i < go->optind; i++) { /* copy all command line options, but not the command line args yet, we may need to append '-s ' before the args */
    esl_strcat(&(cfg->ccom),  -1, go->argv[i], -1);
    esl_strcat(&(cfg->ccom),  -1, " ", 1);
  }
  /* if -s NOT enabled, we need to append the seed info also */
  seed = esl_randomness_GetSeed(cfg->r);
  if(esl_opt_IsUsed(go, "-s")) { /* -s was used on command line, we'll do a sanity check */
    if(seed != esl_opt_GetInteger(go, "-s")) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "get_cmcalibrate_comlog_info(), cfg->r's seed is %" PRIu32 ", but -s was enabled with argument: %" PRIu32 "!, this shouldn't happen.", seed, esl_opt_GetInteger(go, "-s"));
  } else {
    temp = seed; 
    seedlen = 1; 
    while(temp > 0) { temp/=10; seedlen++; } /* determine length of stringized version of seed */
    seedlen += 4; /* strlen(' -s ') */
    ESL_ALLOC(seedstr, sizeof(char) * (seedlen+1));
    sprintf(seedstr, " -s %" PRIu32 " ", seed);
    esl_strcat((&cfg->ccom), -1, seedstr, seedlen);
    free(seedstr);
  }

  for (i = go->optind; i < go->argc; i++) { /* copy command line args yet */
    esl_strcat(&(cfg->ccom), -1, go->argv[i], -1);
    if(i < (go->argc-1)) esl_strcat(&(cfg->ccom), -1, " ", 1);
  }
  
  /* Set the cmcalibrate call date, the cfg->cdate string */
  if((status = esl_strdup(ctime(&date), -1, &(cfg->cdate))) != eslOK) goto ERROR;
  esl_strchop(cfg->cdate, -1); /* doesn't return anything but eslOK */

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_cmcalibrate_comlog_info() error status: %d, probably out of memory.", status);
  return status; 
}


/* Function: update_comlog
 * Date:     EPN, Mon Dec 31 15:14:26 2007
 *
 * Purpose:  Update the CM's comlog info to reflect the current
 *           cmcalibrate call.
 *
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation.
 */
static int
update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm)
{
  int status;
  if(ccom  == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), ccom  is non-NULL.");
  if(cdate == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), cdate is non-NULL.");

  /* free all cmcalibrate comlog info, we're about to overwrite any information that any previous cmcalibrate
   * call could have written to the cm file.
   */
  if(cm->comlog->ccom  != NULL)  { free(cm->comlog->ccom);  cm->comlog->ccom = NULL;  }
  if(cm->comlog->cdate != NULL)  { free(cm->comlog->cdate); cm->comlog->cdate = NULL; }
  
  if((status = esl_strdup(ccom, -1, &(cm->comlog->ccom)))  != eslOK) goto ERROR; 
  if((status = esl_strdup(cdate,-1, &(cm->comlog->cdate))) != eslOK) goto ERROR; 
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "update_comlog() error status: %d, probably out of memory.", status);
  return status; 
}


/* Function: set_dnull
 * Date:     EPN, Thu Jan 24 09:48:54 2008
 *
 * Purpose:  Allocate, fill and return dnull, a double version of cm->null used
 *           for generating random seqs.
 *
 * Returns:  eslOK on success
 */
static int
set_dnull(struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  int status;
  int i;

  if(cfg->dnull != NULL) { free(cfg->dnull); }
  ESL_ALLOC(cfg->dnull, sizeof(double) * cm->abc->K);
  for(i = 0; i < cm->abc->K; i++) cfg->dnull[i] = (double) cm->null[i];
  esl_vec_DNorm(cfg->dnull, cm->abc->K);    

  return eslOK;

 ERROR:
  ESL_FAIL(eslEINCOMPAT, errbuf, "set_dnull(), memory allocation error.");
}


/* Function: print_run_info
 * Date:     EPN, Sun Mar  2 16:57:25 2008
 *
 * Purpose:  Print information on this run of cmcalibrate.
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
  if((status = GetDate    (errbuf, &date))        != eslOK) return status;

  fprintf(stdout, "%-10s %s\n",  "# command:", command);
  fprintf(stdout, "%-10s %s\n",  "# date:",    date);
  fprintf(stdout, "%-10s %" PRIu32 "\n", "# seed:", esl_randomness_GetSeed(cfg->r));
  if(cfg->nproc > 1) fprintf(stdout, "# %-8s %d\n", "nproc:", cfg->nproc);

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

/* Function: print_post_calibration_info
 * Date:     EPN, Wed Mar  5 05:25:02 2008
 * Purpose:  Print info about calibration for a CM we just calibrated including
 *           timings. 
 */
int print_post_calibration_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, double *exp_psecA, double *exp_asecA)
{
  char  time_buf[128];	      /* for printing run time */
  int   exp_mode;             /* counter over exp tail modes */
  double total_psec = 0.;     /* predicted number of seconds for cm, all stages */
  double total_asec = 0.;     /* actual number of seconds for cm, all stages */
  float      L_Mb;            /* total seq length we'll calibrate exp tails on in Mb */
  ExpInfo_t *exp;             /* pointer to current exp tail info, for convenience */

  if(exp_psecA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "print_post_calibration_info, exp_psecA is NULL");
  if(exp_asecA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "print_post_calibration_info, exp_asecA is NULL");

  fprintf(fp, "#\n");
  fprintf(fp, "# Post-calibration info for CM %d: %s\n", cfg->ncm, cm->name);
  fprintf(fp, "#\n");
  fprintf(fp, "# Exponential tail fitting:\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",  "",   "",       "",        "",            "    running time     ");
  fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",  "",   "",       "",        "",            "---------------------");
  fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %10s %10s\n","mod", "cfg", "alg", "L (Mb)",  "mu",     "lambda",   "nhits", "predicted",     "actual");
  fprintf(fp, "# %3s  %3s  %3s %7s %6s %6s %7s %10s %10s\n", "---", "---", "---", "-------", "------", "------", "-------", "----------", "----------");

  for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
    L_Mb = ((float) cfg->expN * (float) cfg->expL) / 1000000.;
    total_psec += exp_psecA[exp_mode];
    total_asec += exp_asecA[exp_mode];
    FormatTimeString(time_buf, exp_psecA[exp_mode], FALSE);
    exp = cfg->expAA[cfg->ncm-1][exp_mode];
    fprintf(fp, "  %-12s %7.2f %6.2f %6.3f %7d %10s", DescribeExpMode(exp_mode), L_Mb, exp->mu_orig, exp->lambda, exp->nrandhits, time_buf);
    FormatTimeString(time_buf, exp_asecA[exp_mode], FALSE);
    fprintf(fp, " %10s\n", time_buf);
  }

  fprintf(fp, "#\n");
  FormatTimeString(time_buf, total_psec, FALSE);
  fprintf(fp, "# total predicted time for CM: %s\n", time_buf);
  FormatTimeString(time_buf, total_asec, FALSE);
  fprintf(fp, "# total actual    time for CM: %s\n", time_buf);
  
  fflush(fp);
  return eslOK;
}
  
/* Function: estimate_time_for_exp_round
 * Date:     EPN, Wed Mar  5 05:46:45 2008
 * Purpose:  Estimate search time for round of exp tail fitting
 *           of exp mode <exp_mode>. This is done by actually 
 *           searching a sequence with the appropriate algorithm. 
 *           The length of the sequence to search is set such 
 *           that it should take about <targ_sec> seconds.
 */
int estimate_time_for_exp_round(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_res)
{
  int    status;
  int    L;                /* length of sequence we'll generate and search to get time estimate */
  double psec_per_Mc;      /* rough prediction at seconds per Mc based on empirical run times I've witnessed, 
                            * doesn't need to be very accurate as we just use it to set length of seq to search to get real prediction */
  float  Mc;               /* millions of DP calculations we're going to do */
  float  Mc_per_res_corr;  /* millions of dp calcs per residue, if searching with CM, corrects         for first W residues requiring less dp calcs */
  float  Mc_per_res;       /* millions of dp calcs per residue, if searching with CM, does not correct for first W residues requiring less dp calcs */
  int    orig_search_opts; /* cm->search_opts when function was entered */
  float  sec_per_res;      /* seconds per residue */
  float  targ_sec = 0.1;   /* target number of seconds our timing expt will take */
  int    Lmin = 100;       /* minimum number of residues to search to get timing */
  int    Lmax = 10000;     /* maximum number of residues to search to get timing */
  int    use_qdb;          /* TRUE if we're using QDB, FALSE if not */
  ESL_DSQ *dsq;            /* the random seq we'll create and search to get predicted time */
  ESL_STOPWATCH *w  = esl_stopwatch_Create(); /* for timings */
  CM_TOPHITS *th    = NULL;

  if(w == NULL)               ESL_FAIL(status, errbuf, "estimate_time_for_exp_round(): memory error, stopwatch not created.\n");
  if(ret_sec_per_res == NULL) ESL_FAIL(status, errbuf, "estimate_time_for_exp_round(): ret_sec_per_res is NULL");

  orig_search_opts = cm->search_opts; /* we'll restore this at end of func, just to be safe */
  if(ExpModeIsInside(exp_mode)) cm->search_opts |= CM_SEARCH_INSIDE;
  else                          cm->search_opts &= ~CM_SEARCH_INSIDE;

  if(cm->smx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "estimate_time_for_exp_round(), cm->smx is NULL");
  use_qdb  = (cm->search_opts & CM_SEARCH_QDB) ? TRUE : FALSE;

  if(use_qdb) { if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*cm->W, NULL,               NULL,               cm->W, TRUE,  NULL, &Mc_per_res_corr)) != eslOK) return status; }
  else        { if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*cm->W, cm->qdbinfo->dmin2, cm->qdbinfo->dmax2, cm->W, TRUE,  NULL, &Mc_per_res_corr)) != eslOK) return status; }

  psec_per_Mc = (cm->search_opts & CM_SEARCH_INSIDE) ? (1. /  75.) : (1. / 275.);  /*  75 Mc/S inside;  275 Mc/S CYK */
  /* determine L that will take about <targ_sec> seconds */
  L = targ_sec / (psec_per_Mc * Mc_per_res_corr);
  L = ESL_MAX(L, Lmin); /* we have to search at least <Lmin> residues */
  L = ESL_MIN(L, Lmax); /* we want to search at most  <Lmax> residues */
  /* now determine exactly how many dp calculations we'd do if we search L residues, 
   * this won't be the same as Mc_per_res * L b/c Mc_per_res from cm_GetNCalcsPerResidueGivenBeta
   * b/c that was calculated after correcting for the fact that the first W residues have fewer
   * DP calcs than all other residues, b/c d < W for first W residues.
   */  
  if(use_qdb) { if((status = cm_CountSearchDPCalcs(cm, errbuf, L, NULL,               NULL,               cm->W, FALSE,  NULL, &Mc_per_res)) != eslOK) return status; }
  else        { if((status = cm_CountSearchDPCalcs(cm, errbuf, L, cm->qdbinfo->dmin2, cm->qdbinfo->dmax2, cm->W, FALSE,  NULL, &Mc_per_res)) != eslOK) return status; }

  /* FALSE says don't correct for fewer dp calcs for first W residues, we want to know how many total DP calcs
   * there will be in L residues */
  Mc = Mc_per_res * L;

  esl_stopwatch_Start(w);
  /* simulate a workunit, generate a sequence, search it, and remove overlaps */
  /*printf("exptL: %d\n", L);*/
  if(esl_opt_GetBoolean(go, "--random")) { 
    if((status = get_random_dsq(cfg, errbuf, cm, L, &dsq)) != eslOK) return status;
  }
  else { 
    if((status = SampleGenomicSequenceFromHMM(cfg->r, cm->abc, errbuf, cfg->ghmm_sA, cfg->ghmm_tAA, cfg->ghmm_eAA, cfg->ghmm_nstates, cfg->expL, &dsq)) != eslOK) cm_Fail(errbuf);
  }

  if((status = process_exp_search_workunit(cm, errbuf, dsq, L, cfg->sc_cutoff, &th)) != eslOK) cm_Fail(errbuf);
  cm_tophits_SortByPosition(th);
  cm_tophits_RemoveOverlaps(th);

  esl_stopwatch_Stop(w);

  cm->search_opts = orig_search_opts;
  sec_per_res = w->user * (Mc_per_res_corr / Mc);
  ESL_DPRINTF1(("L: %d\n", L));
  ESL_DPRINTF1(("w->user: %f\n", w->user));
  ESL_DPRINTF1(("sec_per_res: %f\n", sec_per_res));
  ESL_DPRINTF1(("Mc_per_res: %f\n", Mc_per_res));
  ESL_DPRINTF1(("Mc: %f\n", Mc));

  /*printf("L: %d\n", L);
    printf("w->user: %f\n", w->user);
    printf("sec_per_res: %f\n", sec_per_res);
    printf("Mc_per_res: %f\n", Mc_per_res);
    printf("Mc: %f\n", Mc);*/
  
  if(w != NULL) esl_stopwatch_Destroy(w);
  cm_tophits_Destroy(th);

  free(dsq);
  
  *ret_sec_per_res = sec_per_res;
  return eslOK;
}


/* Function: print_per_cm_column_headings
 * Date:     EPN, Tue Jan  8 05:51:47 2008
 *
 * Purpose:  Print per-CM info to stdout. 
 *
 * Returns:  eslOK on success
 */
static int
print_per_cm_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  printf("#\n");
  if(! esl_opt_IsOn(go, "--forecast")) { 
    printf("# Calibrating CM %d: %s\n", cfg->ncm, cm->name);
    printf("#\n");
    printf("# %8s  %3s  %3s  %3s  %9s %22s\n",           "",      "",    "",    "",    "", "     running time    ");
    printf("# %8s  %3s  %3s  %3s  %9s %22s\n",           "",      "",    "",    "",    "", "----------------------");
    printf("# %-8s  %3s  %3s  %3s  %9s %10s  %10s\n", "stage",    "mod", "cfg", "alg", "expL (Mb)", "predicted", "actual");
    printf("# %8s  %3s  %3s  %3s  %9s %10s  %10s\n",  "--------", "---", "---", "---", "---------", "----------", "----------");
  }
  else { 
    printf("# Forecasting time for %d processor(s) to calibrate CM %d: %s\n", esl_opt_GetInteger(go, "--forecast"), cfg->ncm, cm->name);
    printf("#\n");
    printf("# %-8s  %3s  %3s  %3s  %9s %14s\n", "stage",    "mod", "cfg", "alg", "expL (Mb)", "predicted time");
    printf("# %8s  %3s  %3s  %3s  %9s %14s\n",  "--------", "---", "---", "---", "---------", "--------------");
  }
  return eslOK;
}


/* Function: print_per_cm_summary
 * Date:     EPN, Thu Mar  6 10:35:49 2008
 *
 * Purpose:  Print per-CM summary to stdout. 
 *
 * Returns:  eslOK on success
 */
static int
print_per_cm_summary(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, double psec, double asec)
{
  char  time_buf[128];	      /* for printing run time */
  if(! esl_opt_IsOn(go, "--forecast")) { 
    FormatTimeString(time_buf, psec, FALSE);
    printf("# %8s  %3s  %3s  %3s  %9s %10s  %10s\n", "--------", "---", "---", "---", "---------", "----------", "----------");
    printf("# %-8s  %3s  %3s  %3s  %9s %10s",       "all",        "-",   "-",   "-",         "-", time_buf);
    FormatTimeString(time_buf, asec, FALSE);
    printf("  %10s  (hr:min:sec)\n", time_buf);
  }
  else { 
    FormatTimeString(time_buf, psec, FALSE);
    printf("# %8s  %3s  %3s  %3s  %9s %14s\n", "--------", "---", "---", "---",  "---------", "--------------");
    printf("# %-8s  %3s  %3s  %3s  %9s %14s  (hr:min:sec)\n", "all",        "-",   "-",   "-",         "-",     time_buf);
  }
  return eslOK;
}

/* Function: print_exp_line
 * Date:     EPN, Thu Mar  6 13:38:01 2008
 *
 * Purpose:  Print a line describing exp tail fitting for a given mode.
 */
static int
print_exp_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, int expN, int expL, double psec)
{
  char  time_buf[128];	      /* for printing run time */
  float expL_Mb;              /* exp seq length in Mb */

  FormatTimeString(time_buf, psec, FALSE);

  expL_Mb =  (float) expN * (float) expL; 
  expL_Mb /= 1000000.;

  if(esl_opt_IsOn(go, "--forecast")) { 
    printf("  %-8s  %-12s  %9.2f %14s\n",              "exp tail", DescribeExpMode(exp_mode), expL_Mb, time_buf);
  }
  else { 
    printf("  %-8s  %-12s  %9.2f %10s",              "exp tail", DescribeExpMode(exp_mode), expL_Mb, time_buf);
  }
  fflush(stdout);
  return eslOK;
}


/* Function: process_exp_search_workunit()
 * Date:     EPN, Thu Dec  8 13:48:02 2011
 *
 * Purpose:  Perform search workunit, which consists of a CM, digitized sequence
 *           and indices i and j. The job is to search dsq from i..j and return 
 *           search results in <*ret_results>. Called by cmsearch and cmcalibrate,
 *           which is why it's here and not local in cmsearch.c.
 *
 * Args:     cm              - the covariance model, must have valid searchinfo (si).
 *           errbuf          - char buffer for reporting errors
 *           dsq             - the digitized sequence
 *           L               - length of target sequence 
 *           cutoff          - minimum bit score cutoff to report
 *           ret_th          - search_results_t to create and fill
 *
 * Returns:  eslOK on succes;
 *           <ret_th> is filled with a newly alloc'ed and filled CM_TOPHITS structure, must be freed by caller
 */
int
process_exp_search_workunit(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float cutoff, CM_TOPHITS **ret_th)
{
  int status;
  CM_TOPHITS *th = NULL;
  int use_qdbs = cm->search_opts & CM_SEARCH_QDB ? TRUE : FALSE;

  th = cm_tophits_Create();
  if(th == NULL) ESL_FAIL(eslEMEM, errbuf, "out of memory");
  
  if(cm->search_opts & CM_SEARCH_INSIDE) { 
    if((status = FastIInsideScan(cm, errbuf, cm->smx,                   
				 use_qdbs ? SMX_NOQDB : SMX_QDB2_LOOSE, /* qdbidx, indicates which QDBs to use */
				 dsq, 1, L,                             /* sequence, bounds */
				 cutoff,                                /* minimum score to report */
				 th,                                    /* hitlist to add to */
				 cm->search_opts & CM_SEARCH_NULL3,     /* do the NULL3 correction? */
				 0., NULL, NULL,                        /* vars for redefining envelopes, which we won't do */
				 NULL, NULL))                           /* ret_vsc, ret_sc, irrelevant here */
       != eslOK) return status;
  }
  else { 
   if((status = FastCYKScan(cm, errbuf, cm->smx, 
			    use_qdbs ? SMX_NOQDB : SMX_QDB2_LOOSE, /* qdbidx, indicates which QDBs to use */
			    dsq, 1, L,                             /* sequence, bounds */
			    cutoff,                                /* minimum score to report */
			    th,                                    /* hitlist to add to */
			    cm->search_opts & CM_SEARCH_NULL3,     /* do the NULL3 correction? */
			    0., NULL, NULL,                        /* vars for redefining envelopes, which we won't do */
			    NULL, NULL))                           /* ret_vsc, ret_sc, irrelevant here */
       != eslOK) return status;
  }
  
  *ret_th = th;
  return eslOK;
}

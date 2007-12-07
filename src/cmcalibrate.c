/* cmcalibrate.c
 * Score a CM and a CM Plan 9 HMM against random sequence 
 * data to set the statistical parameters for E-value determination,
 * and CP9 HMM filtering thresholds. 
 * 
 * EPN, Wed May  2 07:02:52 2007
 * based on HMMER-2.3.2's hmmcalibrate.c from SRE
 *
 * MPI example:  
 * qsub -N testrun -j y -R y -b y -cwd -V -pe lam-mpi-tight 32 'mpirun -l C ./mpi-cmcalibrate foo.cm > foo.out'
 * -l forces line buffered output
 *  
 ************************************************************
 * @LICENSE@
 ************************************************************
 */

#include "esl_config.h"
#include "config.h"	

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_histogram.h"
#include "esl_mpi.h"
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#define MPI_NEXT_PARTITION -1 /* message to send to workers */

#include "funcs.h"		/* external functions                   */
#include "structs.h"

#define CUTOPTS  "--eval,--ga,--nc,--tc,--all"  /* Exclusive choice for filter threshold score cutoff */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-s",        eslARG_INT,      "0", NULL, "n>=0",    NULL,      NULL,        NULL, "set random number seed to <n>",   1 },
  { "-t",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print timings for Gumbel fitting and CP9 filter calculation",  1},
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "allow informative insert emissions, do not zero them", 1 },
  /* options for gumbel estimation */
  { "--gumonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--filonly", "only estimate Gumbels, don't calculate filter thresholds", 2},
  { "--cmN",     eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CM gumbel estimation",    2 },
  { "--hmmN",    eslARG_INT,   "5000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CP9 HMM gumbel estimation",    2 },
  { "--dbfile",  eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--filonly", "use GC content distribution from file <s>",  2},
  { "--pfile",   eslARG_STRING,  NULL, NULL, NULL,      NULL,"--dbfile",        NULL, "read partition info for Gumbels from file <s>", 2},
  { "--gumhfile",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save fitted Gumbel histogram(s) to file <s>", 2 },
  { "--gumqqfile",eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save Q-Q plot for Gumbel histogram(s) to file <s>", 2 },
  { "--beta",    eslARG_REAL,  "1e-7", NULL, "x>0",     NULL,      NULL,    "--noqdb", "set tail loss prob for Gumbel calculation to <x>", 5 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "DO NOT use query dependent banding (QDB) Gumbel searches", 5 },
  /* options for filter threshold calculation */
  { "--filonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "only calculate filter thresholds, don't estimate Gumbels", 3},
  { "--filN",    eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--gumonly", "number of emitted sequences for HMM filter threshold calc",    3 },
  { "--fbeta",   eslARG_REAL,  "1e-3",NULL, "x>0",     NULL,      NULL, "--gumonly", "set tail loss prob for filtering sub-CMs QDB to <x>", 5 },
  { "--F",       eslARG_REAL,  "0.95", NULL, "0<x<=1",  NULL,      NULL, "--gumonly", "required fraction of seqs that survive HMM filter", 3},
  { "--fstep",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "step from F to 1.0 while S < Starg", 3},
  { "--starg",   eslARG_REAL,  "0.01", NULL, "0<x<=1",  NULL,      NULL, "--gumonly", "target filter survival fraction", 3},
  { "--spad",    eslARG_REAL,  "1.0",  NULL, "0<=x<=1", NULL,      NULL, "--gumonly", "fraction of (sc(S) - sc(Starg)) to add to sc(S)", 3},
  { "--fast",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "calculate filter thr quickly, assume parsetree sc is optimal", 3},
  { "--gemit",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "when calc'ing filter thresholds, always emit globally from CM",  3},
  { "--filhfile",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--gumonly", "save CP9 filter threshold histogram(s) to file <s>", 3},
  { "--filrfile",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--gumonly", "save CP9 filter threshold information in R format to file <s>", 3},
  /* exclusive choice of filter threshold cutoff */
  { "--eval",    eslARG_REAL,   "0.1", NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "min CM E-val (for a 1MB db) to consider for filter thr calc", 4}, 
  { "--ga",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM gathering threshold as minimum sc for filter thr calc", 4}, 
  { "--nc",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM noise cutoff as minimum sc for filter thr calc", 4}, 
  { "--tc",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM trusted cutoff as minimum sc for filter thr calc", 4},   
  { "--all",     eslARG_NONE,   FALSE, NULL, NULL,   CUTOPTS,      NULL, "--gumonly", "accept all CM hits for filter calc, DO NOT use cutoff", 4}, 
/* Other options */
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb", 5 },  
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program", 5 },  
#endif

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char            *cmfile;	      /* name of input CM file  */ 
  ESL_RANDOMNESS  *r;
  ESL_ALPHABET    *abc;
  double          *gc_freq;;
  double          *pgc_freq;
  int              be_verbose;	
  ESL_STOPWATCH   *w;
  CMStats_t      **cmstatsA;          /* the CM stats data structures, 1 for each CM */
  HybridScanInfo_t *hsi;              /* information for a hybrid scan */ 
  int              ncm;                /* what number CM we're on */
  int              cmalloc;            /* number of cmstats we have allocated */
  char            *tmpfile;            /* tmp file we're writing to */
  char            *mode;               /* write mode, "w" or "wb"                     */
  long             dbsize;             /* size of DB for gumbel stats (impt for E-value cutoffs for filters) */ 
  int              np;                 /* number of partitions, 1 unless --pfile invoked */
  int             *pbegin;             /* [0..p..np-1], begin points for partitions, end pts are implicit */

  /* the following data is modified for each CM, and in some cases for each stat mode for each CM,
   * it is assumed to be 'current' in many functions.
   */
  float           *cutoffA;            /* bit score cutoff for each partition, changes to reflect
				        * current mode CM is in, on masters and workers */
  float           *full_vcalcs;        /* [0..v..cm->M-1] millions of calcs for each subtree to scan 1 residue with --beta  */
  //int             *dmin_gumbel;        /* [0..v..cm->M-1] dmin for current CM and CM configuration, updated in update_qdbs() */
  //int             *dmax_gumbel;        /* [0..v..cm->M-1] dmax for current CM and CM configuration, updated in update_qdbs() */
  double         **vmuAA;              /* [0..np-1][0..cm->M-1], mu for each partition, each state, 
				        * if vmuAA[p][v] == -1 : we're not fitting state v to a gumbel */
  double         **vlambdaAA;          /* same as vmuAA, but lambda */

  /* mpi */
  int              do_mpi;
  int              my_rank;
  int              nproc;
  int              do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *gumhfp;        /* optional output for gumbel histograms */
  FILE         *gumqfp;        /* optional output for gumbel QQ file */
  FILE         *filhfp;        /* optional output for filter histograms */
  FILE         *filrfp;        /* optional output for filter info for R */

};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit Gumbels for E-value stats and calculate HMM filter threshold stats";

static int init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int init_worker_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif
static int process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			    int emit_from_cm, float ***ret_vscAA, float **ret_cp9scA);
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
//static int initialize_sub_filter_info(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm);

static int update_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int fthr_mode);
//static int update_qdbs(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int doing_filter);
static int set_partition_gc_freq(struct cfg_s *cfg, int p);
static int fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, double *ret_mu, double *ret_lambda);
static int cm_fit_histograms(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **vscA, int nscores, int p);
static ESL_DSQ * get_random_dsq(const struct cfg_s *cfg, CM_t *cm, double *dnull, int L);
static ESL_DSQ * get_cmemit_dsq(const struct cfg_s *cfg, CM_t *cm, int *ret_L, int *ret_p, Parsetree_t **ret_tr);
//static float search_target_cm_calibration(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float **ret_vsc);
static int cm_find_hit_above_cutoff(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq, 
				    Parsetree_t *tr, int L, float cutoff);
static int calc_best_filter(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, float **fil_vscAA, float *fil_cp9scA);
static void estimate_workunit_time(const ESL_GETOPTS *go, const struct cfg_s *cfg, int nseq, int L, int gum_mode);
static int read_partition_file(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

int
main(int argc, char **argv)
{

  ESL_GETOPTS     *go	   = NULL;     /* command line processing                     */
  struct cfg_s     cfg;

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
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nGumbel distribution fitting options :");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\ngeneral CP9 HMM filter threshold calculation options :");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\noptions for CM score cutoff to to use for filter threshold calculation:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (! esl_opt_IsDefault(go, "--ga"))
    cm_Fail("--ga not yet implemented, implement it.");
  if (! esl_opt_IsDefault(go, "--nc"))
    cm_Fail("--nc not yet implemented, implement it.");
  if (! esl_opt_IsDefault(go, "--tc"))
    cm_Fail("--tc not yet implemented, implement it.");

  /* Initialize configuration shared across all kinds of masters
   * and workers in this .c file.
   */
  cfg.cmfile  = esl_opt_GetArg(go, 1);
  if (cfg.cmfile == NULL) cm_Fail("Failed to read <cmfile> argument from command line.");
  cfg.cmfp     = NULL;
  cfg.gc_freq  = NULL; 
  cfg.pgc_freq = NULL; 
  cfg.r        = NULL; 
  cfg.w        = NULL; 
  cfg.ncm      = 0;
  cfg.cmstatsA = NULL;
  cfg.hsi      = NULL;
  cfg.tmpfile  = NULL;
  cfg.mode     = NULL;
  cfg.dbsize   = 1000000; /* default DB size = 1MB, changed ONLY if --dbfile enabled */
  cfg.cutoffA  = NULL; 
  cfg.full_vcalcs = NULL;
  //cfg.dmin_gumbel = NULL;
  //cfg.dmax_gumbel = NULL;
  cfg.vmuAA     = NULL;
  cfg.vlambdaAA = NULL;
  cfg.np        = 1;     /* default number of partitions is 1, changed if --pfile */
  cfg.pbegin    = NULL;  /* allocated (by default to size 1) in init_master_cfg() */

  cfg.gumhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gumqfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.filhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.filrfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.abc      = NULL; 

  cfg.do_mpi   = FALSE;
  cfg.my_rank  = 0;
  cfg.nproc    = 0;
  cfg.do_stall = esl_opt_GetBoolean(go, "--stall");

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
      else 		    mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      serial_master(go, &cfg);
    }

  /* Rewind the CM file for a second pass.
   * Write a temporary CM file with new stats information in it
   */
  int   cmi;
  CM_t *cm;
  FILE *outfp;
  sigset_t blocksigs;  /* list of signals to protect from             */
  CMFileRewind(cfg.cmfp);
  if (esl_FileExists(cfg.tmpfile))                    cm_Fail("Ouch. Temporary file %s appeared during the run.", cfg.tmpfile);
  if ((outfp = fopen(cfg.tmpfile, cfg.mode)) == NULL) cm_Fail("Ouch. Temporary file %s couldn't be opened for writing.", cfg.tmpfile); 

  for (cmi = 0; cmi < cfg.ncm; cmi++) {
    if (!CMFileRead(cfg.cmfp, &(cfg.abc), &cm)) cm_Fail("Ran out of CMs too early in pass 2");
    if (cm == NULL)                               cm_Fail("CM file %s was corrupted? Parse failed in pass 2", cfg.cmfile);
    cm->stats = cfg.cmstatsA[cmi];
    if(!(esl_opt_GetBoolean(go, "--filonly"))) cm->flags |= CMH_GUMBEL_STATS; 
    /*if(!(esl_opt_GetBoolean(go, "--gumonly"))) cm->flags |= CMH_FTHR_STATS; */
    CMFileWrite(outfp, cm, cfg.cmfp->is_binary);
    FreeCM(cm);
  } /* end of from idx = 0 to ncm */
      
  /* Now, carefully remove original file and replace it
   * with the tmpfile. Note the protection from signals;
   * we wouldn't want a user to ctrl-C just as we've deleted
   * their CM file but before the new one is moved.
   */
  CMFileClose(cfg.cmfp);
  if (fclose(outfp)   != 0)                            cm_Fail("system error during rewrite of CM file");
  if (sigemptyset(&blocksigs) != 0)                    cm_Fail("system error during rewrite of CM file.");;
  if (sigaddset(&blocksigs, SIGINT) != 0)              cm_Fail("system error during rewrite of CM file.");;
  if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   cm_Fail("system error during rewrite of CM file.");;
  if (remove(cfg.cmfile) != 0)                         cm_Fail("system error during rewrite of CM file.");;
  if (rename(cfg.tmpfile, cfg.cmfile) != 0)            cm_Fail("system error during rewrite of CM file.");;
  if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) cm_Fail("system error during rewrite of CM file.");;
  free(cfg.tmpfile);

  /* Clean up */
  if (cfg.my_rank == 0) {
    if (cfg.gumhfp   != NULL) fclose(cfg.gumhfp);
    if (cfg.gumqfp   != NULL) fclose(cfg.gumqfp);
    if (cfg.filhfp   != NULL) fclose(cfg.filhfp);
    if (cfg.filrfp   != NULL) fclose(cfg.filrfp);
  }
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->gumhfp      - optional output file
 *    cfg->gumqfp      - optional output file
 *    cfg->filhfp      - optional output file
 *    cfg->filrfp      - optional output file
 *    cfg->pbegin      - array of partition starts 
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

  /* optionally, open gumbel histogram file */
  if (esl_opt_GetString(go, "--gumhfile") != NULL) 
    {
      if ((cfg->gumhfp = fopen(esl_opt_GetString(go, "--gumhfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open gumbel histogram save file %s for writing\n", esl_opt_GetString(go, "--gumhfile"));
    }

  /* optionally, open gumbel QQ plot file */
  if (esl_opt_GetString(go, "--gumqqfile") != NULL) 
    {
      if ((cfg->gumqfp = fopen(esl_opt_GetString(go, "--gumqqfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open gumbel QQ plot save file %s for writing\n", esl_opt_GetString(go, "--gumqqfile"));
    }

  /* optionally, open filter threshold calc histogram file */
  if (esl_opt_GetString(go, "--filhfile") != NULL) {
    if ((cfg->filhfp = fopen(esl_opt_GetString(go, "--filhfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --filhfile output file %s\n", esl_opt_GetString(go, "--filhfile"));
    }

  /* optionally, open filter threshold calc info file */
  if (esl_opt_GetString(go, "--filrfile") != NULL) {
    if ((cfg->filrfp = fopen(esl_opt_GetString(go, "--filrfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --filrfile output file %s\n", esl_opt_GetString(go, "--filrfile"));
    }

  /* optionally, get distribution of GC content from an input database (default is use cm->null for GC distro) */
  if(esl_opt_GetString(go, "--dbfile") != NULL) {
    ESL_ALPHABET *tmp_abc = NULL;
    tmp_abc = esl_alphabet_Create(eslRNA);
    ESL_SQFILE      *dbfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--dbfile"), eslSQFILE_UNKNOWN, NULL, &dbfp);
    if (status == eslENOTFOUND)    cm_Fail("No such file."); 
    else if (status == eslEFORMAT) cm_Fail("Format unrecognized."); 
    else if (status == eslEINVAL)  cm_Fail("Can’t autodetect stdin or .gz."); 
    else if (status != eslOK)      cm_Fail("Failed to open sequence database file, code %d.", status); 
    GetDBInfo(tmp_abc, dbfp, &(cfg->dbsize), &(cfg->gc_freq));
    esl_vec_DNorm(cfg->gc_freq, GC_SEGMENTS);
    esl_alphabet_Destroy(cfg->abc);
    esl_sqfile_Close(dbfp); 
   /* allocate pgc_freq, the gc freqs per partition, used to sample seqs for different partitions */
    ESL_ALLOC(cfg->pgc_freq, sizeof(double) * GC_SEGMENTS);
  }

  /* set up the partition data that's used for all CMs */
  if(esl_opt_IsDefault(go, "--pfile")) { /* by default we have 1 partition 0..100 */
    ESL_ALLOC(cfg->pbegin, sizeof(int) * 1);
    cfg->np        = 1;
    cfg->pbegin[0] = 0;
  }
  else { /* setup cfg->np and cfg->pbegin in read_partition_file() */
    if((status = read_partition_file(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  }

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cfg->cmalloc  = 128;
  ESL_ALLOC(cfg->cmstatsA, sizeof(CMStats_t *) * cfg->cmalloc);
  cfg->ncm      = 0;

  /* seed master's RNG */
  if(esl_opt_GetInteger(go, "-s") > 0) {
      if ((cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"))) == NULL)
	cm_Fail("Failed to create random number generator: probably out of memory");
  }
  else {
    if ((cfg->r = esl_randomness_CreateTimeseeded()) == NULL)
      cm_Fail("Failed to create random number generator: probably out of memory");
  }
  printf("Random seed: %ld\n", esl_randomness_GetSeed(cfg->r));

  /* create the stopwatch */
  cfg->w     = esl_stopwatch_Create();

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
  if (esl_FileExists(cfg->tmpfile))
    cm_Fail("temporary file %s already exists; please delete it first", cfg->tmpfile);
  if (cfg->cmfp->is_binary) cfg->mode = "wb";
  else                      cfg->mode = "w"; 

  if(cfg->w == NULL) cm_Fail("Failed to create stopwatch.");
  if(cfg->r == NULL) cm_Fail("Failed to create master RNG.");

  return eslOK;

 ERROR:
  return status;
}

/* init_worker_cfg() 
 * Worker initialization of cfg, worker
 * will get all the info it needs sent to it
 * by the master, so we initialize worker's cfg
 * pointers to NULL, and other values to default.
 * 
 * Because this is called from an MPI worker, it cannot print; 
 * it must return error messages, not print them.
 */
static int
init_worker_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  cfg->cmfile = NULL;
  cfg->abc      = NULL;
  cfg->gc_freq  = NULL;
  cfg->pgc_freq = NULL;
  cfg->be_verbose = FALSE;
  cfg->w = NULL;
  cfg->cmstatsA = NULL;
  cfg->hsi      = NULL;
  cfg->ncm      = 0;
  cfg->cmalloc  = 0;
  cfg->tmpfile  = NULL;
  cfg->mode     = NULL;
  cfg->dbsize  = 0;
  cfg->cutoffA = NULL;
  cfg->full_vcalcs = NULL;
  cfg->vmuAA  = NULL;
  cfg->vlambdaAA = NULL;
  
  cfg->cmfp = NULL;
  cfg->gumhfp = NULL;
  cfg->gumqfp = NULL;
  cfg->filhfp = NULL;
  cfg->filrfp = NULL;
  
  /* we may receive info pertaining to these from master 
   * inside mpi_worker(), at which time we'll update them
   */
  cfg->np       = 0;
  cfg->pbegin   = NULL;
  cfg->r        = NULL; /* we'll create this when seed is received from master */

  return eslOK;
}

/* serial_master()
 * The serial version of cmcalibrate.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int            status;
  char           errbuf[cmERRBUFSIZE];
  CM_t          *cm = NULL;
  int            cmN  = esl_opt_GetInteger(go, "--cmN");
  int            hmmN = esl_opt_GetInteger(go, "--hmmN");
  int            filN = esl_opt_GetInteger(go, "--filN");
  /*int           *vwin = NULL;
    int            nwin = 0;*/
  int        gum_mode = 0;
  int               p;
  float  **gum_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each random seq */
  float  **fil_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each emitted seq */
  float   *gum_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each random seq */
  float   *fil_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each emitted seq */
  int            cmi;
  int      do_time_estimates = TRUE;
  /* int      do_time_estimates = esl_opt_GetBoolean(go, "--etime"); */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;
      cmi = cfg->ncm-1;

      if((status = initialize_cm(go, cfg, errbuf, cm))          != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf))         != eslOK) cm_Fail(errbuf);
      //if((status = initialize_sub_filter_info(go, cfg, cm))     != eslOK) cm_Fail("initialize_sub_info failed unexpectedly (status: %d)", status);

      for(gum_mode = 0; gum_mode < NGUMBELMODES; gum_mode++) {
	//for(gum_mode = 0; gum_mode < 5; gum_mode += 4) { /* will fit CM_LC and CP9_LV only */
	//for(gum_mode = 1; gum_mode < 2; gum_mode++) {
	  ConfigForGumbelMode(cm, gum_mode);
	  /* recalc QDBs for fitting Gumbels */
	  //update_qdbs(go, cfg, cm, FALSE);
	  for (p = 0; p < cfg->cmstatsA[cmi]->np; p++) 
	    {
	      if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	      if(gum_mode < NCMMODES) 
		{ /* CM mode */
		  /* search random sequences to get gumbels for each candidate sub CM root state (including 0, the root of the full model) */
		  if(do_time_estimates) estimate_workunit_time(go, cfg, cmN, cm->W*2, gum_mode);
		  ESL_DPRINTF1(("\n\ncalling process_workunit to fit gumbel for p: %d CM mode: %d\n", p, gum_mode));
		  if((status = process_workunit (go, cfg, errbuf, cm, cmN,   FALSE, &gum_vscAA, NULL)) != eslOK) cm_Fail(errbuf);
		  if((status = cm_fit_histograms(go, cfg, errbuf, cm, gum_vscAA, cmN, p))              != eslOK) cm_Fail(errbuf);
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->mu     = cfg->vmuAA[p][0];
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->lambda = cfg->vlambdaAA[p][0];
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->L      = cm->W*2;
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->N      = cmN;
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->isvalid = TRUE;
		}
	      else 
		{ /* CP9 mode, fit gumbel for full HMM */
		  if(do_time_estimates) estimate_workunit_time(go, cfg, cmN, cm->W*2, gum_mode);
		  ESL_DPRINTF1(("\n\ncalling process_workunit to fit gumbel for p: %d CP9 mode: %d\n", p, gum_mode));
		  if((status = process_workunit (go, cfg, errbuf, cm, hmmN, FALSE, NULL, &gum_cp9scA))                 != eslOK) cm_Fail(errbuf);
		  if((status = fit_histogram(go, cfg, errbuf, gum_cp9scA, hmmN, &(cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->mu), 
					     &(cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->lambda)))     != eslOK) cm_Fail(errbuf);
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->L      = cm->W*2;
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->N      = hmmN;
		  cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->isvalid = TRUE;
		}
	    }
	  continue; /* temporarily skip filter thr section */
	  /* filter threshold calculation section */
	  if((! esl_opt_GetBoolean(go, "--gumonly")) && gum_mode < NCMMODES) {
	    if((status = update_cutoffs(go, cfg, errbuf, cm, gum_mode)) != eslOK) cm_Fail(errbuf);
	    //update_qdbs(go, cfg, cm, TRUE); /* recalc QDBs for filtering sub CMs */
	    printf("\n\ncutoffsA[0]: %10.4f\n\n", cfg->cutoffA[0]);

	    /* search emitted sequences to get filter thresholds for HMM and each candidate sub CM root state */
	    ESL_DPRINTF1(("\n\ncalling process_workunit to get filter threshold for p: %d mode: %d\n", p, gum_mode));
	    if((status = process_workunit (go, cfg, errbuf, cm, filN, TRUE, &fil_vscAA, &fil_cp9scA)) != eslOK) cm_Fail(errbuf);

	    /* determine the 'best' way to filter using a heuristic:
	     * 1. predict speedup for HMM-only filter
	     * 2. add 'best' sub CM root independent of those already chosen, and predict speedup for a hybrid scan
	     * 3. repeat 2 until predicted speedup drops
	     * best filter is either HMM only or fastest combo of sub CM roots
	     */
	    if((status = calc_best_filter(go, cfg, cm, fil_vscAA, fil_cp9scA)) != eslOK) cm_Fail("unexpected error code: %d in calc_best_filter.");

	    /* determine best n-state sets of sub-CM filter states */
	    //if((status = calc_sub_filter_sets(go, cfg, cm, fil_vscAA)) != eslOK) cm_Fail("unexpected error code: %d in calc_sub_filter_sets.");
	       
	    /* determine HMM filter threshold */
	    /* if((status = calc_cp9_filter_thr(go, cfg, errbuf, cm, fil_vscAA)) != eslOK) cm_Fail("unexpected error code: %d in calc_cp9_filter_thr.");*/
	  }
	  /* free muA and lambdaA */
	} /* end of for(gum_mode = 0; gum_mode < NGUMBELMODES1; gum_mode++) */
      debug_print_cmstats(cfg->cmstatsA[cmi], (! esl_opt_GetBoolean(go, "--gumonly")));
      if(cfg->hsi != NULL) cm_FreeHybridScanInfo(cfg->hsi, cm);
      cfg->hsi = NULL;
      FreeCM(cm);
    }
  return;
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmcalibrate
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

  CM_t          *cm = NULL;
  int            cmN  = esl_opt_GetInteger(go, "--cmN");
  int            hmmN = esl_opt_GetInteger(go, "--hmmN");
  int            filN = esl_opt_GetInteger(go, "--filN");

  int        gum_mode = 0;
  int               p;
  float  **gum_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each random seq */
  float  **fil_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each emitted seq */
  float   *gum_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each random seq */
  float   *fil_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each emitted seq */
  float  **worker_vscAA = NULL;
  float   *worker_cp9scA = NULL;

  long *seedlist = NULL;
  char  errbuf[cmERRBUFSIZE];
  MPI_Status mpistatus; 
  int   n, v, i;
  int working_on_cm;        /* TRUE when gum_mode is for CM gumbel */
  int working_on_cp9;       /* TRUE when gum_mode is for CP9 gumbel */
  int nseq_sent        = 0; /* number of seqs we've told workers to work on */
  int nseq_per_worker  = 0; /* number of seqs to tell each worker to work on */
  int nseq_this_worker = 0; /* number of seqs to tell current worker to work on */
  int nseq_this_round  = 0; /* number of seqs for current round */
  int nseq_just_recv   = 0; /* number of seqs we just received scores for from a worker */
  int nseq_recv        = 0; /* number of seqs we've received thus far this round from workers */
  int msg;
  int cmi;                  /* CM index, which number CM we're working on */

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL)    { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seedlist  = malloc(sizeof(long) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }

  ESL_ALLOC(seedlist, sizeof(long) * cfg->nproc);
  for (wi = 0; wi < cfg->nproc; wi++) 
    {
      /* not sure what to do here */
      seedlist[wi] = wi;
    }
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) cm_Fail(errbuf);
  ESL_DPRINTF1(("MPI master is initialized\n"));

  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  for (wi = 1; wi < cfg->nproc; wi++)
    MPI_Send(&(seedlist[wi]), 1, MPI_LONG, wi, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) cm_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));

  /* 2 special (annoying) case:
   * case 1: if we've used the --dbfile option, we read in a seq file to fill
   * cfg->gc_freq, and we need to broadcast that info to workers
   *
   * case 2: if we are calculating stats for more than 1 partition, 
   * (--pfile invoked), we need to broadcast that information to 
   * the workers. 
   */
  if(! (esl_opt_IsDefault(go, "--dbfile"))) { /* receive gc_freq info from master */
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  if(! (esl_opt_IsDefault(go, "--pfile"))) { /* broadcast partition info to workers */
    MPI_Bcast(&(cfg->np),  1,       MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pbegin != NULL));
    MPI_Bcast(cfg->pbegin, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
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

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      cfg->ncm++;  
      cmi = cfg->ncm-1;
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm(go, cfg, errbuf, cm))          != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf))         != eslOK) cm_Fail(errbuf);

      ESL_ALLOC(gum_vscAA, sizeof(float *) * cm->M);
      ESL_ALLOC(gum_cp9scA, sizeof(float)  * hmmN);
      for(v = 0; v < cm->M; v++) ESL_ALLOC(gum_vscAA[v], sizeof(float) * cmN);

      for(gum_mode = 0; gum_mode < NGUMBELMODES; gum_mode++) {
	ConfigForGumbelMode(cm, gum_mode);
	working_on_cm   = (gum_mode < NCMMODES) ? TRUE  : FALSE;
	working_on_cp9  = (gum_mode < NCMMODES) ? FALSE : TRUE;
	nseq_per_worker = working_on_cm ? (int) (cmN / (cfg->nproc-1)) : (int) (hmmN / (cfg->nproc-1));
	nseq_this_round = working_on_cm ? cmN : hmmN;

	for (p = 0; p < cfg->np; p++) {

	  ESL_DPRINTF1(("MPI master: CM: %d gumbel mode: %d partition: %d\n", cfg->ncm, gum_mode, p));

	  have_work     = TRUE;	/* TRUE while work remains  */
	  esl_stopwatch_Start(cfg->w);  /* Start timing. */
	  
	  wi = 1;
	  nseq_sent = 0;
	  nseq_recv = 0;
	  while (have_work || nproc_working)
	    {
	      if(have_work) { 
		if(nseq_sent < nseq_this_round) {
		  nseq_this_worker = (nseq_sent + nseq_per_worker <= nseq_this_round) ? 
		    nseq_per_worker : (nseq_this_round - nseq_sent);
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
		    ESL_DPRINTF1(("MPI master sees that the result buffer contains calibration results\n"));
		    if(working_on_cm) {
		      if ((status = cmcalibrate_cm_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, cm->M, &worker_vscAA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
		      ESL_DPRINTF1(("MPI master has unpacked CM gumbel results\n"));
		      ESL_DASSERT1((nseq_just_recv > 0));
		      for(v = 0; v < cm->M; v++) {
			for(i = 0; i < nseq_just_recv; i++) {
			  ESL_DPRINTF2(("\tscore from worker v: %d i: %d sc: %f\n", i, v, worker_vscAA[v][i]));
			  gum_vscAA[v][nseq_recv+i] = worker_vscAA[v][i];
			}
			free(worker_vscAA[v]);
		      }
		      free(worker_vscAA);
		    }
		    else { /* working on cp9 */
		      if ((status = cmcalibrate_cp9_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_cp9scA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
		      ESL_DPRINTF1(("MPI master has unpacked CP9 gumbel results\n"));
		      ESL_DASSERT1((nseq_just_recv > 0));
		      for(i = 0; i < nseq_just_recv; i++) 
			gum_cp9scA[nseq_recv+i] = worker_cp9scA[i];
		      free(worker_cp9scA);
		    }
		    nseq_recv += nseq_just_recv;
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
		  /* send new search job */
		  ESL_DPRINTF1(("MPI master is sending nseq %d to worker %d\n", nseq_this_worker, wi));
		  MPI_Send(&(nseq_this_worker), 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	      
		  wi++;
		  nproc_working++;
		  nseq_sent += nseq_this_worker;
		}
	    }
	  /* fit gumbels for this partition p, this gumbel mode gum_mode */
	  if(working_on_cm) { 
	    if((status = cm_fit_histograms(go, cfg, errbuf, cm, gum_vscAA, cmN, p)) != eslOK) cm_Fail(errbuf);
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->mu     = cfg->vmuAA[p][0];
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->lambda = cfg->vlambdaAA[p][0];
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->L      = cm->W*2;
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->N      = cmN;
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->isvalid = TRUE;
	  }
	  else /* working on CP9 */ {
	    if((status = fit_histogram(go, cfg, errbuf, gum_cp9scA, hmmN, &(cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->mu), 
				       &(cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->lambda)))     != eslOK) cm_Fail(errbuf);
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->L      = cm->W*2;
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->N      = hmmN;
	    cfg->cmstatsA[cmi]->gumAA[gum_mode][p]->isvalid = TRUE;
	  }

	  ESL_DPRINTF1(("MPI master: done with partition: %d for gumbel mode: %d for this CM. Telling all workers\n", p, gum_mode));
	  for (wi = 1; wi < cfg->nproc; wi++) { 
	    msg = MPI_NEXT_PARTITION;
	    MPI_Send(&msg, 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	  }
	}
	ESL_DPRINTF1(("MPI master: done with gumbel mode %d for this CM.\n", gum_mode));
      }
      ESL_DPRINTF1(("MPI master: done with this CM.\n"));
      debug_print_cmstats(cfg->cmstatsA[cmi], (! esl_opt_GetBoolean(go, "--gumonly")));
      for(v = 0; v < cm->M; v++) free(gum_vscAA[v]);
      free(gum_vscAA);
      free(gum_cp9scA);
      FreeCM(cm);
    }
  
  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((status = cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if (xstatus != eslOK) cm_Fail(errbuf);
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
  CM_t         *cm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
  int           pos;
  char          errbuf[cmERRBUFSIZE];
  MPI_Status  mpistatus;
  float  **gum_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each random seq */
  float   *gum_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each random seq */
  long     seed;  /* seed for RNG */
  int      gum_mode;
  int working_on_cm;        /* TRUE when gum_mode is for CM gumbel */
  int working_on_cp9;       /* TRUE when gum_mode is for CP9 gumbel */
  int nseq;
  int v, p;
  int cmi;

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
	   
  /* Master now sends worker initialization information (RNG seed) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  if (MPI_Recv(&seed, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  if (xstatus == eslOK) { if ((status = init_worker_cfg(go, cfg, errbuf)) != eslOK)   xstatus = status;  }
  if (xstatus == eslOK) { if((cfg->r = esl_randomness_Create(seed)) == NULL)          xstatus = eslEMEM; }
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized seed: %ld\n", cfg->my_rank, seed));

  /* 2 special (annoying) cases: 
   * case 1: if we've used the --dbfile option, we read in a seq file to fill
   * cfg->gc_freq, and we need that info here for the worker, so we receive
   * it's broadcast from the master
   * 
   * case 2: if we are calculating stats for more than 1 
   * partition, (--pfile invoked), we need to receive that information 
   * via broadcast from master. Otherwise we need to setup the default partition info
   * (single partition, 0..100 GC content)
   */
  if(! (esl_opt_IsDefault(go, "--dbfile"))) { /* receive gc_freq info from master */
    ESL_DASSERT1((cfg->gc_freq == NULL));
    ESL_ALLOC(cfg->gc_freq,  sizeof(double) * GC_SEGMENTS);
    ESL_ALLOC(cfg->pgc_freq, sizeof(double) * GC_SEGMENTS);
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else cfg->gc_freq = NULL; /* default */
  if(! (esl_opt_IsDefault(go, "--pfile"))) { /* receive partition info from master */
    MPI_Bcast(&(cfg->np),     1, MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pbegin == NULL));
    ESL_ALLOC(cfg->pbegin, sizeof(int) * cfg->np);
    MPI_Bcast(cfg->pbegin, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else { /* no --pfile, set up default partition info */  
    cfg->np     = 1;
    ESL_ALLOC(cfg->pbegin, sizeof(int) * cfg->np);
    cfg->pbegin = 0;
  }
  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      
      for(gum_mode = 0; gum_mode < NGUMBELMODES; gum_mode++) {
	ESL_DPRINTF1(("worker: %d gum_mode: %d nparts: %d\n", cfg->my_rank, gum_mode, cfg->np));
	ConfigForGumbelMode(cm, gum_mode);
	working_on_cm   = (gum_mode < NCMMODES) ? TRUE  : FALSE;
	working_on_cp9  = (gum_mode < NCMMODES) ? FALSE : TRUE;
	
	for (p = 0; p < cfg->np; p++) { /* for each partition */
	    
	  ESL_DPRINTF1(("worker %d gum_mode: %d partition: %d\n", cfg->my_rank, gum_mode, p));
	  if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	  
	  if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
	  while(nseq != MPI_NEXT_PARTITION) {
	    ESL_DPRINTF1(("worker %d: has received nseq: %d\n", cfg->my_rank, nseq));
	    
	    if(working_on_cm) {
	      if((status = process_workunit (go, cfg, errbuf, cm, nseq, FALSE, &gum_vscAA, NULL)) != eslOK) goto CLEANERROR;
	    }
	    else  { /* working on cp9 */
	      if((status = process_workunit (go, cfg, errbuf, cm, nseq, FALSE, NULL, &gum_cp9scA)) != eslOK) goto CLEANERROR;
	    }
	    ESL_DPRINTF1(("worker %d: has gathered gumbel results\n", cfg->my_rank));
	    n = 0;
	    if(working_on_cm) { 
	      if (cmcalibrate_cm_results_MPIPackSize(gum_vscAA, nseq, cm->M, MPI_COMM_WORLD, &sz) != eslOK) goto CLEANERROR; n += sz;  
	    }
	    else { /* working on cp9 */
	      if (cmcalibrate_cp9_results_MPIPackSize(gum_cp9scA, nseq, MPI_COMM_WORLD, &sz) != eslOK) goto CLEANERROR; n += sz;  
	    }
	    if (n > wn) {
	      void *tmp;
	      ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	      wn = n;
	    }
	    ESL_DPRINTF1(("worker %d: has calculated the gumbel results will pack into %d bytes\n", cfg->my_rank, n));
	    status = eslOK;
	    pos = 0;
	    if(working_on_cm) {
	      if (cmcalibrate_cm_results_MPIPack(gum_vscAA, nseq, cm->M, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) goto ERROR;
	    }
	    else {
	      if (cmcalibrate_cp9_results_MPIPack(gum_cp9scA, nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) goto ERROR;
	    }	    
	    MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	    ESL_DPRINTF1(("worker %d: has sent gumbel results to master in message of %d bytes\n", cfg->my_rank, pos));
	    
	    if(working_on_cm) { 
	      for(v = 0; v < cm->M; v++) free(gum_vscAA[v]);
	      free(gum_vscAA);
	    }
	    else { /* working on cp9 */
	      free(gum_cp9scA);
	    }
	    /* receive next number of sequences, if MPI_NEXT_PARTITION, we'll stop */
	    if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
	  }
	  ESL_DPRINTF1(("worker %d gum_mode: %d finished partition: %d\n", cfg->my_rank, gum_mode, p));
	}
	ESL_DPRINTF1(("worker %d finished all partitions for gum_mode: %d\n", cfg->my_rank, gum_mode));
      }
      FreeCM(cm);
      cm = NULL;
      ESL_DPRINTF1(("worker %d finished all gum_modes for this cm.\n", cfg->my_rank));
    }
  if (status == eslEOD) ESL_DPRINTF1(("Worker %d told CMs are done.\n", cfg->my_rank));
  else goto ERROR;
  
  if (wbuf != NULL) free(wbuf);
  return;

 CLEANERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
  return;

 ERROR:
  cm_Fail("Allocation error in mpi_worker");
  return;
}
#endif /*HAVE_MPI*/


/* Function: process_workunit()
 * Date:     EPN, Mon Sep 10 16:55:09 2007
 *
 * Purpose:  A work unit consists of a CM, a int specifying a number of 
 *           sequences <nseq>, and a flag indicated how to generate those
 *           sequences. The job is to generate <nseq> sequences and search
 *           them with a CM and/or CP9, saving scores, which are returned.
 *
 *           This function can be run in 1 of 3 modes, determined by the
 *           status of the input variables:
 *         
 *           Mode 1. Gumbel calculation for CM. 
 *           <emit_from_cm> is FALSE, <ret_vscAA> != NULL, <ret_cp9scA> == NULL
 *           Emit randomly and search only with the CM. <ret_vscAA> is filled
 *           with the best CM score at each state for each sequence.
 *
 *           Mode 2. Gumbel calculation for CP9.
 *           <emit_from_cm> is FALSE, <ret_vscAA> == NULL, <ret_cp9scA> != NULL
 *           Emit randomly and search only with the CP9. <ret_cp9scA> is filled
 *           with the best CP9 score for each sequence.
 *
 *           Mode 3. Scores will eventually be used to determine filter thresholds.
 *           <emit_from_cm> is TRUE, <ret_vscAA> != NULL, <ret_cp9scA> != NULL
 *           Emit from the CM (which is already configured how we want it). Search
 *           with the CM and CP9. <ret_vscAA> and <ret_cp9scA> filled as explained
 *           in modes 1 and 2. Importantly, in this mode, each sequence must have
 *           a NON-BANDED CYK scan hit above a given cutoff. That cutoff is given
 *           as a bit score in cfg->cutoffA[p], where p is the partition for the
 *           sequence (p is determined in get_cmemit_dsq() called from this function). 
 *           Sequences that have no hit better than cutoff are not accepted. 
 *           The cutoff[p] value is assumed to be already set before this
 *           function is entered.
 *
 *           The ability to run 3 different modes complicates the code a bit,
 *           but I prefered it to making 2 separate functions b/c a significant
 *           part of those 2 functions would have identical code.
 *
 * Args:     go           - getopts
 *           cfg          - cmcalibrate's configuration
 *           errbuf       - for writing out error messages
 *           cm           - the CM (already configured as we want it)
 *           nseq         - number of seqs to generate
 *           emit_from_cm - TRUE to emit from CM; FALSE emit random 
 *           ret_vscAA    - RETURN: [0..v..cm->M-1][0..nseq-1] best 
 *                                  score at each state v for each seq
 *           ret_cp9scA   - RETURN: [0..nseq-1] best CP9 score for each seq
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
		 int emit_from_cm, float ***ret_vscAA, float **ret_cp9scA)
{
  int            status;
  int            use_cm      = FALSE;
  int            use_cp9     = FALSE;
  float        **vscAA       = NULL;  /* [0..v..cm->M-1][0..i..nseq-1] best CM score for each state, each seq */
  float         *cur_vscA    = NULL;  /* [0..v..cm->M-1]               best CM score for each state cur seq */
  float         *cp9scA      = NULL;  /*                [0..i..nseq-1] best CP9 score for each seq */
  double        *dnull       = NULL;
  int            p;                   /* what partition we're in, not used unless emit_from_cm = TRUE */
  int            i;
  int            v;
  int            L;
  int            nfailed = 0;
  Parsetree_t   *tr;
  ESL_DSQ       *dsq;
  float          sc1;
  /* float          sc2; */
  /*  float          sc3;*/

  /* contract check */
  if(ret_vscAA  == NULL && ret_cp9scA == NULL)                  ESL_FAIL(eslEINVAL, errbuf, "process_workunit, ret_vscAA and ret_cp9scA both NULL."); 
  if(ret_vscAA  != NULL && ret_cp9scA != NULL && !emit_from_cm) ESL_FAIL(eslEINVAL, errbuf, "process_workunit, ret_vscAA and ret_cp9scA both non-NULL, but emit_from_cm is FALSE.");

  /* determine mode, are we calc'ing CM stats, CP9 stats, or both */
  if(ret_vscAA   != NULL) use_cm  = TRUE;
  if(ret_cp9scA  != NULL) use_cp9 = TRUE;
  
  ESL_DPRINTF1(("in process_workunit use_cm: %d use_cp9: %d emit_from_cm: %d nseq: %d ", use_cm, use_cp9, emit_from_cm, nseq));
  if(emit_from_cm == FALSE) {
    if(use_cm)  ESL_DPRINTF1(("mode 1 CM gumbel\n"));
    if(use_cp9) ESL_DPRINTF1(("mode 2 CP9 gumbel\n"));
  }
  else ESL_DPRINTF1(("mode 3 filter thr\n"));

  int do_cyk     = TRUE;
  int do_inside  = FALSE;
  int do_viterbi = FALSE;
  int do_forward = FALSE;
  if(cm->search_opts & CM_SEARCH_INSIDE)     { do_cyk = FALSE; do_inside = TRUE; }
  if(cm->search_opts & CM_SEARCH_HMMVITERBI) do_viterbi = TRUE;
  if(cm->search_opts & CM_SEARCH_HMMFORWARD) do_forward = TRUE;
  if((do_viterbi + do_forward) > 1)               ESL_FAIL(eslEINVAL, errbuf, "process_workunit, cm->search_opts CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags both raised.");
  if(use_cp9 && ((do_viterbi + do_forward) != 1)) ESL_FAIL(eslEINVAL, errbuf, "process_workunit, use_cp9: %d use_cm: %d but CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags both down.", use_cp9, use_cm);
  ESL_DPRINTF1(("do_cyk:     %d\ndo_inside:  %d\ndo_viterbi: %d\ndo_forward: %d\n", do_cyk, do_inside, do_viterbi, do_forward)); 
  ESL_DPRINTF1(("cm->M: %d\n", cm->M));
  
  /* allocate the score arrays we'll pass back */
  if(use_cm) { 
    ESL_ALLOC(vscAA, sizeof(float *) * cm->M);
    for(v = 0; v < cm->M; v++) ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    ESL_ALLOC(cur_vscA, sizeof(float) * cm->M);
  }
  if(use_cp9) ESL_ALLOC(cp9scA, sizeof(float) * nseq);
  
  /* fill dnull, a double version of cm->null, but only if we're going to need it to generate random seqs */
  if(!emit_from_cm && cfg->pgc_freq == NULL) {
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);    
  }
  ESL_DPRINTF1(("cm->M: %d\n", cm->M));
  
  /* generate dsqs one at a time and collect best CM scores at each state and/or best overall CP9 score */
  for(i = 0; i < nseq; i++) {
    if(emit_from_cm) { /* if emit_from_cm == TRUE, use_cm == TRUE */
      if(nfailed > 1000 * nseq) { cm_Fail("Max number of failures (%d) reached while trying to emit %d seqs.\n", nfailed, nseq); }
      dsq = get_cmemit_dsq(cfg, cm, &L, &p, &tr);
      while(!(cm_find_hit_above_cutoff(go, cfg, errbuf, cm, dsq, tr, L, cfg->cutoffA[p]))) { 
	free(dsq); 	
	/* parsetree tr is freed in cm_find_hit_above_cutoff() */
	dsq = get_cmemit_dsq(cfg, cm, &L, &p, &tr);
	nfailed++;
      }
    }
    else { 
      dsq = get_random_dsq(cfg, cm, dnull, cm->W*2); 
      L = cm->W*2; 
    }
    /* if nec, search with CM */
    if (use_cm) {
      if(cm->search_opts & CM_SEARCH_INSIDE) { /* inside */
	if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc1)) != eslOK) return status;
      }
      else { /* cyk */
	if((status = FastCYKScan    (cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc1)) != eslOK) return status;
      }
      /* sc1 = search_target_cm_calibration(cm, dsq, cm->dmin, cm->dmax, 1, L, cm->W, &(cur_vscA)); */
      /*sc2 = actually_search_target(cm, dsq, 1, L, 0., 0., NULL, FALSE, FALSE, FALSE, NULL, FALSE); */
      /*printf("i: %4d sc1: %10.4f sc2: %10.4f\n", i, sc1, sc2);
	fflush(stdout);
	if(fabs(sc1 - sc2) > 0.01) 
	cm_Fail("i: %4d sc1: %10.4f sc2: %10.4f\n", i, sc1, sc2);*/
      /*if(fabs(sc1 - sc3) > 0.01) 
	cm_Fail("i: %4d sc1: %10.4f sc3: %10.4f\n", i, sc1, sc3);
      */
      for(v = 0; v < cm->M; v++) vscAA[v][i] = cur_vscA[v];
      free(cur_vscA);
    }
    /* if nec, search with CP9 */
    if(use_cp9) {
      if(cm->search_opts & CM_SEARCH_HMMFORWARD) { /* forward */
	if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
				 TRUE,   /* yes, we are scanning */
				 FALSE,  /* no, we are not aligning */
				 FALSE,  /* don't be memory efficient */
				 NULL,   /* don't want best score at each posn back */
				 NULL,   /* don't want the max scoring posn back */
				 &(cp9scA[i]))) != eslOK) return status;
      }
      else { /* viterbi */
	if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
				 TRUE,   /* yes, we are scanning */
				 FALSE,  /* no, we are not aligning */
				 FALSE,  /* don't be memory efficient */
				 NULL,   /* don't want best score at each posn back */
				 NULL,   /* don't want the max scoring posn back */
				 NULL,   /* don't want traces back */
				 &(cp9scA[i]))) != eslOK) return status;
      }
    }
    free(dsq);
  }
  if(dnull != NULL) free(dnull);
  if(use_cm)  *ret_vscAA  = vscAA;
  if(use_cp9) *ret_cp9scA = cp9scA;
  return eslOK;

 ERROR:
  return status;
}


/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;

  cm->beta   = esl_opt_GetReal(go, "--beta"); /* this will be 1e-7 (default beta) unless changed at command line */

  if(!(esl_opt_GetBoolean(go, "--iins"))) cm->config_opts |= CM_CONFIG_ZEROINSERTS;

  /* config QDB? Yes, unless --noqdb enabled */
  if(! (esl_opt_GetBoolean(go, "--noqdb"))) 
    cm->config_opts |= CM_CONFIG_QDB;   /* configure QDB */
  else
    cm->search_opts |= CM_SEARCH_NOQDB; /* don't use QDB to search */

  ConfigCM(cm, NULL, NULL);
  
  /* count number of DP calcs */
  if(cfg->full_vcalcs != NULL) free(cfg->full_vcalcs);
  if((status = cm_CountSearchDPCalcs(cm, errbuf, 1000, cm->dmin, cm->dmax, cm->W, &(cfg->full_vcalcs), NULL)) != eslOK) return status;

  /* create and initialize scan info for CYK/Inside scanning functions */
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("initialize_cm(), CreateScanMatrixForCM() call failed.");
  
  /* create and initialize hybrid scan info */
  ESL_DASSERT1((cfg->hsi == NULL));
  cfg->hsi = cm_CreateHybridScanInfo(cm, esl_opt_GetReal(go, "--fbeta"), cfg->full_vcalcs[0]);
  if(cfg->hsi == NULL) cm_Fail("initialize_cm(), CreateHybridScanInfo() call failed.");

  return eslOK;
}


/* initialize_cmstats()
 * Allocate and initialize a cmstats object in the cfg->cmstatsA array. 
 */
static int
initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int i;
  int p;
  int cmi = cfg->ncm-1;
  int np = cfg->np;

  ESL_DPRINTF1(("initializing cmstats for %d partitions\n", np));

  cfg->cmstatsA[cmi] = AllocCMStats(np);
  ESL_DASSERT1((cfg->pbegin[0] == 0));
  for(p = 0; p < np;     p++) cfg->cmstatsA[cmi]->ps[p] = cfg->pbegin[p];
  for(p = 0; p < (np-1); p++) cfg->cmstatsA[cmi]->pe[p] = cfg->pbegin[p+1]-1;
  cfg->cmstatsA[cmi]->pe[(np-1)] = GC_SEGMENTS-1; /* this is 100 */
  
  for(p = 0; p < np; p++)
    for(i = cfg->cmstatsA[cmi]->ps[p]; i <= cfg->cmstatsA[cmi]->pe[p]; i++)
      cfg->cmstatsA[cmi]->gc2p[i] = p; 
  
  if(cfg->cutoffA != NULL) free(cfg->cutoffA);
  ESL_ALLOC(cfg->cutoffA, sizeof(float) * np); /* number of partitions */
  
  /* allocate cfg->vmuAA, and cfg->vlambdaAA if they're not already allocated */
  if(cfg->vmuAA == NULL) {
    ESL_ALLOC(cfg->vmuAA,     sizeof(double *) * cfg->cmstatsA[cmi]->np);
    for(p = 0; p < cfg->cmstatsA[cmi]->np; p++) {
      cfg->vmuAA[p]     = NULL;
    }
  }
  if(cfg->vlambdaAA == NULL) { 
    ESL_ALLOC(cfg->vlambdaAA, sizeof(double *) * cfg->cmstatsA[cmi]->np);
    for(p = 0; p < cfg->cmstatsA[cmi]->np; p++) {
      cfg->vlambdaAA[p] = NULL;
    }
  }
  return eslOK;
    
  ERROR:
  sprintf(errbuf, "initialize_cmstats(), memory allocation error (status: %d).", status);
  return status;
}

/* update_cutoffs()
 * Update the cfg->cutoffA array to have the bit score cutoff for each partition
 * for the 'current' cm (number ncm-1).
 */
static int
update_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int fthr_mode)
{
  double         tmp_K;          /* used for recalc'ing Gumbel stats for DB size */
  double         ecutoff;       /* filled if we're using an E-value cutoff */
  int            p;
  double         mu;

  if(esl_opt_GetBoolean(go, "--all")) {
    for (p = 0; p < cfg->cmstatsA[cfg->ncm-1]->np; p++)
      cfg->cutoffA[p] = -eslINFINITY;
  }
  else if (esl_opt_GetBoolean(go, "--ga"))   cm_Fail("update_cutoffs() --ga not yet implemented.");
  else if (esl_opt_GetBoolean(go, "--nc"))   cm_Fail("update_cutoffs() --ga not yet implemented.");
  else if (esl_opt_GetBoolean(go, "--tc"))   cm_Fail("update_cutoffs() --ga not yet implemented.");
  else ecutoff = esl_opt_GetReal(go, "--eval"); /* default, use --eval cutoff */

  for (p = 0; p < cfg->cmstatsA[cfg->ncm-1]->np; p++)
    {
      /* First determine mu based on db_size */
      tmp_K = exp(cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->mu * cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda) / 
	cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->L;
      mu = log(tmp_K  * ((double) cfg->dbsize)) / cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cfg->cutoffA[p] = mu - (log(ecutoff) / cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda);
    }
  return eslOK;
}  

/* Function: set_partition_gc_freq()
 * Date:     EPN, Mon Sep 10 08:00:27 2007
 *
 * Purpose:  Set up the GC freq to sample from for the current partition. 
 *           Only used if --dbfile used to read in dbseq from which to derive
 *           GC distributions for >= 1 partition.
 *
 * Returns:  eslOK on success;
 */
int
set_partition_gc_freq(struct cfg_s *cfg, int p)
{
  int i, begin, end;
  ESL_DASSERT1((cfg->pgc_freq != NULL));
  ESL_DASSERT1((cfg->gc_freq != NULL));

  esl_vec_DSet(cfg->pgc_freq, GC_SEGMENTS, 0.);
  begin = cfg->pbegin[p];
  if(p == (cfg->np-1)) end = (GC_SEGMENTS-1); /* this is 100 */
  else end = cfg->pbegin[p+1] - 1;
  for (i = begin; i <= end; i++) 
    cfg->pgc_freq[i] = cfg->gc_freq[i];
  esl_vec_DNorm(cfg->pgc_freq, GC_SEGMENTS);

  return eslOK;
}

/* fit_histogram()
 * Create, fill and fit a histogram to a gum. Data to fill the histogram
 * is given as <data>.
 */
static int
fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores,
	      double *ret_mu, double *ret_lambda)
{
  double mu;
  double lambda;
  int i;
  double *xv;         /* raw data from histogram */
  int     n,z;  

  ESL_HISTOGRAM *h = NULL;       /* histogram of scores */


  /* Initialize histogram; these numbers are guesses */
  h = esl_histogram_CreateFull(-100., 100., .25);    

  /* fill histogram */
  for(i = 0; i < nscores; i++)
    esl_histogram_Add(h, scores[i]);

  /* fit scores to a gumbel */
  esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
  esl_gumbel_FitCensored(xv, n, z, xv[0], &mu, &lambda);

  /* print to output files if nec */
  if(cfg->gumhfp != NULL)
    esl_histogram_Plot(cfg->gumhfp, h);
  if(cfg->gumqfp != NULL) {
      double  params[2];  
      params[0] = mu;
      params[1] = lambda;
      esl_histogram_PlotQQ(cfg->gumqfp, h, &esl_exp_generic_invcdf, params);
  }

  esl_histogram_Destroy(h);

  *ret_mu     = mu;
  *ret_lambda = lambda;
  return eslOK;
}

/* cm_fit_histograms()
 * We want gumbels for each cm state we can do a legal local begin into.
 * Call fit_histogram() for each such state.
 */
static int
cm_fit_histograms(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, 
		  float **vscA, int nscores, int p)
{
  int status;
  int v;

  if(cfg->vmuAA[p]     != NULL) free(cfg->vmuAA[p]);
  if(cfg->vlambdaAA[p] != NULL) free(cfg->vlambdaAA[p]);
  
  ESL_ALLOC(cfg->vmuAA[p],     sizeof(double) * cm->M);
  ESL_ALLOC(cfg->vlambdaAA[p], sizeof(double) * cm->M);

  for(v = 0; v < cm->M; v++) {
    if(cfg->hsi->iscandA[v]) {
      /* printf("FITTING v: %d sttype: %d\n", v, cm->sttype[v]); */
      fit_histogram(go, cfg, errbuf, vscA[v], nscores, &(cfg->vmuAA[p][v]), &(cfg->vlambdaAA[p][v]));
    }
    else cfg->vmuAA[p][v] = cfg->vlambdaAA[p][v] = -1.;
  }

  return eslOK;

 ERROR:
  return status;
}

/* Function: pick_stemwinners()
 * Date:     
 *
 * Returns:  eslOK on success;
 */
int 
pick_stemwinners(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, double *muA, double *lambdaA, float *avglen, int **ret_vwin, int *ret_nwin)
{
  return eslOK;
}

/* Function: get_random_dsq()
 * Date:     EPN, Tue Sep 11 08:31:47 2007
 * 
 * Purpose:  Generate a random digitized seq and return it.
 *           Two possible modes:
 *           1. if(cfg->pgc_freq == NULL && dnull != NULL) 
 *              use dnull disto (a double version of cm->null) to generate
 *           2. if(cfg->pgc_freq != NULL && dnull == NULL) 
 *              use choose a GC frequency from cfg->pgc_freq
 *              and generate with that
 *
 * Returns:  eslOK on success;
 */
ESL_DSQ *
get_random_dsq(const struct cfg_s *cfg, CM_t *cm, double *dnull, int L)
{
  int status;
  double  gc_comp;
  double *distro = NULL;
  int do_free_distro = FALSE;
  ESL_DSQ *dsq = NULL;

  /* contract check, make sure we're in a valid mode */
  if(cfg->pgc_freq == NULL && dnull == NULL) cm_Fail("get_random_dsq(), cfg->pgc_freq == NULL and dnull == NULL");
  if(cfg->pgc_freq != NULL && dnull != NULL) cm_Fail("get_random_dsq(), cfg->pgc_freq != NULL and dnull != NULL");

  /* determine mode */ /* generate sequence */
  if      (cfg->pgc_freq == NULL && dnull != NULL) distro = dnull;
  else if (cfg->pgc_freq != NULL && dnull == NULL) {
    assert(cm->abc->K == 4);
    ESL_ALLOC(distro, sizeof(double) * cm->abc->K);
    do_free_distro = TRUE;
    gc_comp = 0.01 * esl_rnd_DChoose(cfg->r, cfg->pgc_freq, GC_SEGMENTS);
    distro[1] = distro[2] = 0.5 * gc_comp;
    distro[0] = distro[3] = 0.5 * (1. - gc_comp);
  }
  /* generate sequence */
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if (esl_rnd_xIID(cfg->r, distro, cm->abc->K, L, dsq) != eslOK) cm_Fail("get_random_dsq(): failure creating random sequence.");

  if (do_free_distro) free(distro);
  return dsq;

 ERROR:
  cm_Fail("get_random_dsq() memory allocation error.");
  return NULL; /*NEVERREACHED */
}

/* Function: get_cmemit_dsq()
 * Date:     EPN, Tue Sep 11 08:51:33 2007
 * 
 * Purpose:  Generate a dsq from a CM and return it.
 *
 * Returns:  eslOK on success;
 */
ESL_DSQ *
get_cmemit_dsq(const struct cfg_s *cfg, CM_t *cm, int *ret_L, int *ret_p, Parsetree_t **ret_tr)
{
  int p;
  int L;
  ESL_SQ *sq;
  ESL_DSQ *dsq;
  Parsetree_t *tr;

  EmitParsetree(cm, cfg->r, "irrelevant", TRUE, &tr, &sq, &L);
  while(L == 0) { 
    FreeParsetree(tr); 
    esl_sq_Destroy(sq); 
    EmitParsetree(cm, cfg->r, "irrelevant", TRUE, &tr, &sq, &L);
  }

  /* determine the partition */
  p = cfg->cmstatsA[cfg->ncm-1]->gc2p[(get_gc_comp(sq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */

  /* free everything allocated by a esl_sqio.c:esl_sq_CreateFrom() call, but the dsq */
  dsq = sq->dsq;
  free(sq->name);
  free(sq->acc);
  free(sq->desc);
  free(sq);

  *ret_L  = L;
  *ret_p  = p;
  *ret_tr = tr;
  return dsq;
}

/* Function: calc_best_filter()
 * Date:     EPN, Thu Nov  1 15:05:03 2007
 *
 * Purpose:  Given a CM and scores for a CP9 and CM scan of target seqs
 *           determine the best filter we could use, either an HMM only
 *           or a hybrid scan with >= 1 sub CM roots.
 *            
 * Returns:  eslOK on success;
 *           Dies immediately on an error.
 */
int
calc_best_filter(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, float **fil_vscAA, float *fil_cp9scA)
{
  int    status;
  int    v;
  float  *sorted_fil_cp9scA;
  float **sorted_fil_vscAA;
  float **sorted_fil_EAA;
  int    filN  = esl_opt_GetInteger(go, "--filN");
  int    Fidx;
  float  sc;
  float  E;
  float  fil_calcs;
  float  surv_calcs;
  float  nonfil_calcs;
  float  spdup;
  int    i;
  int    cmi = cfg->ncm-1;
  
  float F = esl_opt_GetReal(go, "--F");
  Fidx  = (int) (1. - F) * filN;

  if(cfg->cmstatsA[cfg->ncm-1]->np != 1) cm_Fail("calc_sub_filter_sets(), not yet implemented for multiple partitions.\nYou'll need to keep track of partition of each sequence OR\nstore E-values not scores inside process_workunit.");

  /* Determine the predicted CP9 filter speedup */
  
  ESL_ALLOC(sorted_fil_cp9scA, sizeof(float) * filN);
  esl_vec_FCopy(fil_cp9scA, filN, sorted_fil_cp9scA); 
  esl_vec_FSortIncreasing(sorted_fil_cp9scA, filN);
  sc = sorted_fil_cp9scA[Fidx];

  printf("\n\n***********************************\nvscAA[0] scores:\n");
  for(i = 0; i < filN; i++)
    printf("HMM i: %4d sc: %10.4f\n", i, sorted_fil_cp9scA[i]);
  printf("***********************************\n\n");

  /* E is expected number of hits for db of length 2 * cm->W */
  /* SHOULD THIS BE CP9_LF???? */
  E  = RJK_ExtremeValueE(sc, cfg->cmstatsA[cmi]->gumAA[CP9_LF][0]->mu, cfg->cmstatsA[cmi]->gumAA[CP9_LF][0]->lambda);
  fil_calcs  = cfg->hsi->full_cp9_ncalcs;
  surv_calcs = E * (cm->W * 2) * cfg->hsi->full_cm_ncalcs;
  printf("HMM(L)        sc: %10.4f E: %10.4f filt: %10.4f ", sc, E, fil_calcs);
  fil_calcs += surv_calcs;
  nonfil_calcs = cfg->hsi->full_cm_ncalcs;
  spdup = nonfil_calcs / fil_calcs;
  printf("surv: %10.4f sum : %10.4f full: %10.4f spdup %10.4f\n", surv_calcs, fil_calcs, nonfil_calcs, spdup);

  exit(1);
  /*****************TEMPORARY PRINTF BEGIN***************************/
  ESL_ALLOC(sorted_fil_vscAA, sizeof(float *) * cm->M);
  ESL_ALLOC(sorted_fil_EAA, sizeof(float *) * cm->M);

  /*printf("\n\n***********************************\nvscAA[0] scores:\n");
  for(i = 0; i < filN; i++)
    printf("i: %4d sc: %10.4f\n", i, fil_vscAA[0][i]);
    printf("***********************************\n\n");*/

  for(v = 0; v < cm->M; v++)
    {
      ESL_ALLOC(sorted_fil_vscAA[v], sizeof(float) * filN);
      esl_vec_FCopy(fil_vscAA[v], filN, sorted_fil_vscAA[v]); 
      esl_vec_FSortIncreasing(sorted_fil_vscAA[v], filN);
    }
  for(v = 0; v < cm->M; v++) {
    if(cfg->hsi->iscandA[v]) {
      sc = sorted_fil_vscAA[v][Fidx];
      /* E is expected number of hits for db of length 2 * cm->W */
      E  = RJK_ExtremeValueE(sc, cfg->vmuAA[0][v], cfg->vlambdaAA[0][v]);
      /* note partition = 0, this is bogus if more than 1 partition, that's why we die if there are more (see above). */
      fil_calcs  = cfg->hsi->cm_vcalcs[v];
      surv_calcs = E * (cm->W * 2) * cfg->hsi->full_cm_ncalcs;
      printf("SUB %3d sg: %2d sc: %10.4f E: %10.4f filt: %10.4f ", v, cfg->hsi->startA[v], sc, E, fil_calcs);
      fil_calcs += surv_calcs;
      nonfil_calcs = cfg->hsi->full_cm_ncalcs;
      spdup = nonfil_calcs / fil_calcs;
      printf("surv: %10.4f sum : %10.4f full: %10.4f spdup %10.4f\n", surv_calcs, fil_calcs, nonfil_calcs, spdup);
    }  
  }
  
  for(v = 0; v < cm->M; v++) 
    {
      ESL_ALLOC(sorted_fil_EAA[v], sizeof(float *) * filN);
      /* E is expected number of hits for db of length 2 * cm->W */
      /* assumes only 1 partition */
      /*sorted_fil_EAA[v] = RJK_ExtremeValueE(sorted_fil_vscAA[v], cfg->vmuAA[0][v], cfg->vlambdaAA[0][v]);*/
    }      

  /*****************TEMPORARY PRINTF END***************************/

  return eslOK;

 ERROR:
  cm_Fail("calc_sub_filter_sets(), memory allocation error.");
  return status; /* NEVERREACHED */
}


/*
 * Function: cm_find_hit_above_cutoff()
 * Date:     EPN, Wed Sep 12 04:59:08 2007
 *
 * Purpose:  Given a CM, a sequence, and a cutoff, try to 
 *           *quickly* answer the question: Does this sequence 
 *           contain a hit to the CM above the cutoff?
 *           To do this we first check the parsetree score, and
 *           then do do up to 3 iterations of search.
 *           The first 2 are performend with j and d bands 
 *           (of decreasing tightness), then default 
 *           search (with QDB unless --noqdb enabled) is done.
 *           We return TRUE if any search finds a hit above
 *           cutoff, and FALSE otherwise.
 *
 * Args:     go              - getopts
 *           cfg             - cmcalibrate's configuration
 *           errbuf          - char buffer for error message
 *           cm              - CM to emit from
 *           dsq             - the digitized sequence to search
 *           tr              - parsetree for dsq
 *           L               - length of sequence
 *           cutoff          - bit score cutoff 
 *
 * Returns:  TRUE if a hit above cutoff is found.
 *           FALSE otherwise.
 *           Dies immediately on failure with informative error message.
 */
int 
cm_find_hit_above_cutoff(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq,
			 Parsetree_t *tr, int L, float cutoff)
{
  int status;
  int init_flags       = cm->flags;
  int init_search_opts = cm->search_opts;
  int turn_qdb_back_on = FALSE;
  int turn_hbanded_back_off = FALSE;
  int turn_hmmscanbands_back_off = FALSE;
  double orig_tau = cm->tau;
  float sc;

  /* Determine if this sequence has a hit in it above the cutoff as quickly as possible. 
   * Stage 0: Check parsetree score
   * Stage 1: HMM banded search tau = 1e-2
   * Stage 2: HMM banded search with scanning bands, tau = 1e-10
   * Stage 3: QDB search (CYK or inside), beta = --beta, (THIS IS MOST LENIENT SEARCH WE'LL DO)
   *
   * The earliest stage at which we find a hit > cutoff at any stage, we return cm->flags, cm->search_opts
   * to how they were when we entered, and return TRUE.
   *
   * NOTE: We don't do a full non-banded parse to be 100% sure we don't exceed the cutoff, 
   * unless --noqdb was enabled (ScanMatrix_t *smx stores dn/dx (min/max d) for each state), 
   * because we assume the --beta value used in *this* cmcalibrate 
   * run will also be used for any cmsearch runs.
   */

  sc = ParsetreeScore(cm, tr, dsq, FALSE); 
  FreeParsetree(tr);
  if(sc > cutoff) { /* parse score exceeds cutoff */
    assert(cm->flags       == init_flags);
    assert(cm->search_opts == init_search_opts);
    /* printf("0 sc: %10.4f\n", sc); */
    return TRUE;
  } 

  if(!(cm->search_opts & CM_SEARCH_NOQDB))        turn_qdb_back_on = TRUE;
  if(!(cm->search_opts & CM_SEARCH_HBANDED))      turn_hbanded_back_off = TRUE;
  if(!(cm->search_opts & CM_SEARCH_HMMSCANBANDS)) turn_hmmscanbands_back_off = TRUE;

  cm->search_opts |= CM_SEARCH_NOQDB;

  /* stage 1 */
  /*cm->search_opts |= CM_SEARCH_HBANDED;
  cm->tau = 0.01;
  if((sc = actually_search_target(cm, dsq, 1, L, 0., 0., NULL, FALSE, FALSE, FALSE, NULL, FALSE)) > cutoff) { 
    if(turn_qdb_back_on)        cm->search_opts &= ~CM_SEARCH_NOQDB; 
    if(turn_hbanded_back_off) { cm->search_opts &= ~CM_SEARCH_HBANDED; cm->tau = orig_tau; }
    assert(cm->flags       == init_flags);
    assert(cm->search_opts == init_search_opts);
    printf("1 sc: %10.4f ", sc);
    return TRUE;
    } */

  /* stage 2 */
  /*cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
  cm->tau = 1e-10;
  if((sc = actually_search_target(cm, dsq, 1, L, 0., 0., NULL, FALSE, FALSE, FALSE, NULL, FALSE)) > cutoff) {
    if(turn_qdb_back_on)             cm->search_opts &= ~CM_SEARCH_NOQDB; 
    if(turn_hbanded_back_off)      { cm->search_opts &= ~CM_SEARCH_HBANDED;      cm->tau = orig_tau; }
    if(turn_hmmscanbands_back_off) { cm->search_opts &= ~CM_SEARCH_HMMSCANBANDS; cm->tau = orig_tau; }
    assert(cm->flags       == init_flags);
    assert(cm->search_opts == init_search_opts);
    printf("2 sc: %10.4f ", sc);
    return TRUE;
    } */

  /* stage 3, use 'default' dmin, dmax (which could be NULL) */
  cm->search_opts &= ~CM_SEARCH_HBANDED;
  cm->search_opts &= ~CM_SEARCH_HMMSCANBANDS;
  if(turn_qdb_back_on) cm->search_opts &= ~CM_SEARCH_NOQDB; 

  if((status = FastCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &sc)) != eslOK) return status;
  /* sc = search_target_cm_calibration(cm, dsq, dmin_default_beta, dmax_default_beta, 1, L, cm->W, NULL); */
  if(!turn_hbanded_back_off)      { cm->search_opts |= CM_SEARCH_HBANDED;      cm->tau = orig_tau; }
  if(!turn_hmmscanbands_back_off) { cm->search_opts |= CM_SEARCH_HMMSCANBANDS; cm->tau = orig_tau; }
  assert(cm->flags       == init_flags);
  assert(cm->search_opts == init_search_opts);

  if(sc > cutoff) { printf("3 sc: %10.4f\n", sc);  return TRUE; }
  else return FALSE;
}

/* Function: estimate_workunit_time()
 * Date:     EPN, Thu Nov  1 17:57:20 2007
 * 
 * Purpose:  Estimate time req'd for a cmcalibrate workunit
 *
 * Returns:  eslOK on success;
 */
void
estimate_workunit_time(const ESL_GETOPTS *go, const struct cfg_s *cfg, int nseq, int L, int gum_mode)
{
  /* these are ballparks for a 3 GHz machine with optimized code */
  float cyk_megacalcs_per_sec = 250.;
  float ins_megacalcs_per_sec =  40.;
  float fwd_megacalcs_per_sec = 175.;
  float vit_megacalcs_per_sec = 380.;
  
  float seconds = 0.;

  switch(gum_mode) { 
  case CM_LC: 
  case CM_GC: 
    seconds = cfg->full_vcalcs[0] * (float) L * (float) nseq / cyk_megacalcs_per_sec;
    break;
  case CM_LI:
  case CM_GI:
    seconds = cfg->full_vcalcs[0] * (float) L * (float) nseq / ins_megacalcs_per_sec;
    break;
  case CP9_LV: 
  case CP9_GV: 
    seconds = cfg->hsi->full_cp9_ncalcs * (float) L * (float) nseq / vit_megacalcs_per_sec;
    break;
  case CP9_LF: 
  case CP9_GF: 
    seconds = cfg->hsi->full_cp9_ncalcs * (float) L * (float) nseq / fwd_megacalcs_per_sec;
    break;
  }
  printf("Estimated time for this workunit: %10.2f seconds\n", seconds);

  return;
}


/* Function: read_partition_file
 * Date:     EPN, Fri Dec  7 08:38:41 2007
 * 
 * Called when --pfile is invoked. 
 * Opens and reads a partition file of 
 * with 2 * <npartitions> tokens, every odd token is
 * a partition start <pbegin>, and every even token is 
 * a parititon end <pend>. First <pbegin> must be 0,
 * other <pbegin>s must be 1 more than previous
 * <pend>. The last <pend> must be 100, other <pends>
 * must be 1 less than following <pbegin>.
 *
 * Example of file that implies 3 partitions: 
 * 0..39, 40..60, and 61.100
 * 
 * ~~~~~~~~~~~~~~~~
 * 0 39
 * 40 60
 * 61 100
 * ~~~~~~~~~~~~~~~~
 * 
 * After reading the file and checking it's legit,
 * set up the cfg->np and cfg->pbegin data.
 *
 * Returns:  eslOK on success, eslEINVAL if file is 
 *           in wrong format, or doesn't follow rules described above.
 */
int
read_partition_file(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;
  int            *begin;
  int             end=0;
  int             nread=0;
  int             p;

  printf("in read_partition_file, mp: %d gc: %d\n", MAX_PARTITIONS, GC_SEGMENTS);

  ESL_DASSERT1((MAX_PARTITIONS < GC_SEGMENTS));
  if(esl_opt_IsDefault(go, "--pfile")) ESL_FAIL(eslEINVAL, errbuf, "read_partition_file, but --pfile not invoked!\n");

  if (esl_fileparser_Open(esl_opt_GetString(go, "--pfile"), &efp) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "failed to open %s in read_mask_file\n", esl_opt_GetString(go, "--pfile"));
  esl_fileparser_SetCommentChar(efp, '#');
  

  ESL_ALLOC(begin, sizeof(int) * GC_SEGMENTS);
  begin[0] = 0;

  while((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslEOF) {
    begin[nread] = atoi(tok);
    if(nread == 0) {
      if(atoi(tok) != 0) ESL_FAIL(eslEINVAL, errbuf, "first partition begin must be 0 in %s\n", esl_opt_GetString(go, "--pfile"));
    }
    else if (begin[nread] != (end+1)) {
      if(atoi(tok) != 0) ESL_FAIL(eslEINVAL, errbuf, "partition %d begin point (%d) is not exactly 1 more than prev partition end pt %d in %s\n", (nread+1), begin[nread], end, esl_opt_GetString(go, "--pfile"));
    }      
    if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "no end point for each partition %d's begin (%d) in partition file %s\n", (nread+1), begin[nread], esl_opt_GetString(go, "--pfile"));
    end = atoi(tok);
    if(end < begin[nread]) ESL_FAIL(eslEINVAL, errbuf, "partition %d end point (%d) < begin point (%d) in %s\n", (nread+1), end, begin[nread], esl_opt_GetString(go, "--pfile"));
    nread++;
    if(nread > MAX_PARTITIONS) ESL_FAIL(eslEINVAL, errbuf, "partition file %s has at least %d partitions, but max num partitions is %d\n", esl_opt_GetString(go, "--pfile"), nread, MAX_PARTITIONS);
  }
  if(nread == 0) ESL_FAIL(eslEINVAL, errbuf, "failed to read a single token from %s\n", esl_opt_GetString(go, "--pfile"));
  if(end != 100) ESL_FAIL(eslEINVAL, errbuf, "final partitions end point must be 100, but it's %d in %s\n", end, esl_opt_GetString(go, "--pfile"));

  /* create cfg->pbegin */
  ESL_DASSERT1((cfg->pbegin == NULL));
  ESL_ALLOC(cfg->pbegin, sizeof(int) * nread);
  for(p = 0; p < nread; p++) cfg->pbegin[p] = begin[p];
  free(begin);
  cfg->np = nread;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return status;
}

#if 0
/*
 * Function: cm_emit_seqs_to_aln_above_cutoff()
 * Date:     EPN, Mon Sep 10 17:31:36 2007
 *
 * Purpose:  Create a seqs_to_aln object by generating sequences
 *           from a CM. Only accept sequences that have a CM hit
 *           within them above a bit score cutoff.
 *
 * Args:     go              - getopts
 *           cfg             - cmcalibrate's configuration
 *           cm              - CM to emit from
 *           nseq            - number of seqs to emit
 *           cutoff          - bit score cutoff 
 *
 * Returns:  Ptr to a newly allocated seqs_to_aln object with nseq sequences to align.
 *           Dies immediately on failure with informative error message.
 */
seqs_to_aln_t *cm_emit_seqs_to_aln_above_cutoff(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int nseq)
{
  int status;
  seqs_to_aln_t *seqs_to_aln = NULL;
  char *name = NULL;
  int namelen;
  int L;
  int i;
  int do_cyk = TRUE;
  Parsetree_t *tr = NULL;

  if(cm->dmin == NULL || cm->dmax == NULL) cm_Fail("cm_emit_seqs_to_aln_above_cutoff(), dmin, dmax are NULL.");
  if(cm->search_opts & CM_SEARCH_NOQDB) cm_Fail("cm_emit_seqs_to_aln_above_cutoff(), search opt NOQDB enabled.");

  seqs_to_aln = CreateSeqsToAln(nseq, FALSE);

  namelen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
  if(cm->name != NULL) namelen += strlen(cm->name) + 1;
  ESL_ALLOC(name, sizeof(char) * namelen);

  while(i < nseq)
    {
      if(cm->name != NULL) sprintf(name, "%s-%d", cm->name, i+1);
      else                 sprintf(name, "%d-%d", cfg->ncm-1, i+1);
      L = 0; 
      EmitParsetree(cm, cfg->r, name, TRUE, &tr, &(seqs_to_aln->sq[i]), &L);
      while(L == 0) { FreeParsetree(tr); esl_sq_Destroy(seqs_to_aln->sq[i]); EmitParsetree(cm, cfg->r, name, TRUE, &tr, &(seqs_to_aln->sq[i]), &L); }
      p = cfg->cmstatsA[(ncm-1)]->gc2p[(get_gc_comp(seqs_to_aln->sq[i], 1, L))]; /* in get_gc_comp() should be i and j of best hit */

      sc = ParsetreeScore(cm, tr, seqs_to_aln->sq[i]->dsq, FALSE); 
      FreeParsetree(tr);
      if(sc > cfg->cutoffA[p]) { i++; continue; }

      /* If we get here, parse score is not above cfg->cutoffA[p], we want to determine if 
       * this sequence has a hit in it above the cfg->cutoffA[p] as quickly as possible. 
       *
       * Stage 1: HMM banded search tau = 1e-2
       * Stage 2: HMM banded search with scanning bands, tau = 1e-10
       * Stage 3: QDB search (CYK or inside), beta = cm->beta (should be default beta)
       *
       * If we find a hit > cfg->cutoffA[p] at any stage, we accept the seq, increment i and move on.
       * We don't do a full non-banded parse to ensure that we don't exceed the cfg->cutoffA[p], 
       * because QDB is on in cmsearch by default.
       */

      /* stage 1 */
      cm->search_opts |= CM_SEARCH_HBANDED;
      cm->tau = 0.01;
      if((sc = actually_search_target(cm, seqs_to_aln->sq[i]->dsq, 1, L, 0., 0., NULL, FALSE, FALSE, FALSE, NULL, FALSE)) > cfg->cutoffA[p]) 
	{ i++; break; }
      s1_np++;
      /* stage 2 */
      cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
      cm->tau = 1e-10;
      if((sc = actually_search_target(cm, seqs_to_aln->sq[i]->dsq, 1, L, 0., 0., NULL, FALSE, FALSE, FALSE, NULL, FALSE)) > cfg->cutoffA[p]) 
	{ i++; break; }
      s2_np++;
      /* stage 3 */
      cm->search_opts &= ~CM_SEARCH_HBANDED;
      cm->search_opts &= ~CM_SEARCH_HBANDED;
      if((sc = search_target_cm_calibration(cm, seqs_to_aln->sq[i]->dsq, cm->dmin, cm->dmax, 1, seqs_to_aln->sq[i]->n, cm->W, NULL)) > cfg->cutoffA[p])
	{ i++; break; }
      s3_np++;
      if(s3_np > (1000 * nseq)) cm_Fail("cm_emit_seqs_to_aln_above_cutoff(), wanted %d seqs above cutoff: %d bits, reached limit of %d seqs\n", nseq, cfg->cutoffA[p], (1000 * nseq));

      /* didn't pass */
      esl_sq_Destroy(seqs_to_aln->sq[i]);
    }

  seqs_to_aln->nseq = nseq;

  free(name);
  return seqs_to_aln;

 ERROR:
  cm_Fail("memory allocation error");
  return NULL;
}


/* Function: search_target_cm_calibration() based on bandcyk.c:CYKBandedScan()
 * Date:     EPN, Sun Sep  9 19:05:07 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm. If bands are NULL, reverts to non-banded
 *           (scancyk.c:CYKScan()). 
 *
 *           Special local cmcalibrate function that only cares about
 *           collecting the best score at each state for the sequence.
 *
 * Args:     cm        - the covariance model
 *           dsq       - the digitized sequence
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           ret_vsc  - RETURN: [0..v..M-1] best score at each state v
 *
 * Returns:  score of best overall hit (vsc[0]).
 *           dies immediately if some error occurs.
 */
float 
search_target_cm_calibration(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float **ret_vsc)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int    ***ialpha;             /* Inside DP score matrix, [v][j][d] */
  float    *vsc;           /* best score for each state (float) */
  float    *ivsc;          /* best score for each state (int, only used if do_inside) */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_w;  	        /* offset j for state w */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       dn;                 /* temporary value for min d in for loops */
  int       dx;                 /* temporary value for max d in for loops */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */

  int       do_cyk    = FALSE;  /* TRUE: do cyk; FALSE: do_inside */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  float     ret_val;

  /* determine if we're doing cyk/inside banded/non-banded */
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) do_cyk    = TRUE;
  if(dmin != NULL && dmax != NULL)           do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 

  ESL_ALLOC(vsc, sizeof(float) * cm->M);
  esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);

  if(do_cyk) { 

    /*****************************************************************
     * alpha allocations.
     * The scanning matrix is indexed [v][j][d]. 
     *    v ranges from 0..M-1 over states in the model.
     *    j takes values 0 or 1: only the previous (prv) or current (cur) row
     *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
     *      deck in memory, and j ranges from 0..W, and yes it must be square
     *      because we'll use a rolling pointer trick thru it
     *    d ranges from 0..W over subsequence lengths.
     * Note that E memory is shared: all E decks point at M-1 deck.
     *****************************************************************/
    ESL_ALLOC(alpha, (sizeof(float **) * cm->M));
    for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
      if (cm->stid[v] == BEGL_S)
	{
	  ESL_ALLOC(alpha[v], (sizeof(float *) * (W+1)));
	  for (j = 0; j <= W; j++)
	    ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
	}
      else if (cm->sttype[v] == E_st && v < cm->M-1) 
	alpha[v] = alpha[cm->M-1];
      else 
	{
	  ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	  for (j = 0; j < 2; j++) 
	    ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
	}
    }

    /*****************************************************************
     * alpha initializations.
     * We initialize on d=0, subsequences of length 0; these are
     * j-independent. Any generating state (P,L,R) is impossible on d=0.
     * E=0 for d=0. B,S,D must be calculated. 
     * Also, for MP, d=1 is impossible.
     * Also, for E, all d>0 are impossible.
     *
     * and, for banding: any cell outside our bands is impossible.
     * These inits are never changed in the recursion, so even with the
     * rolling, matrix face reuse strategy, this works.
     *****************************************************************/ 
    for (v = cm->M-1; v >= 0; v--)
      {
	alpha[v][0][0] = IMPOSSIBLE;

	if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
	else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = IMPOSSIBLE;
	else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	  {
	    y = cm->cfirst[v];
	    alpha[v][0][0] = cm->endsc[v];
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	      alpha[v][0][0] = ESL_MAX(alpha[v][0][0], (alpha[y+yoffset][0][0] + cm->tsc[v][yoffset]));
	    /* ...we don't bother to look at local alignment starts here... */
	    alpha[v][0][0] = ESL_MAX(alpha[v][0][0], IMPOSSIBLE);
	  }
	else if (cm->sttype[v] == B_st) 
	  {
	    w = cm->cfirst[v];
	    y = cm->cnum[v];
	    alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	  }

	alpha[v][1][0] = alpha[v][0][0];
	if (cm->stid[v] == BEGL_S) 
	  for (j = 2; j <= W; j++) 
	    alpha[v][j][0] = alpha[v][0][0];
      }
    /* Impose the bands.
     *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
     *    the first loop will be skipped, the second initializes the E states.)
     */
    if(do_banded) { 
      for (v = 0; v < cm->M; v++) {
	jmax = (cm->stid[v] == BEGL_S) ? W : 1;
	for (d = 0; d < dmin[v] && d <=W; d++) 
	  for(j = 0; j <= jmax; j++)
	    alpha[v][j][d] = IMPOSSIBLE;
	for (d = dmax[v]+1; d <= W;      d++) 
	  for(j = 0; j <= jmax; j++)
	    alpha[v][j][d] = IMPOSSIBLE;
      }
    }

    /* The main loop: scan the sequence from position i0 to j0.
     */
    for (j = i0; j <= j0; j++) 
      {
	cur = j%2;
	prv = (j-1)%2;
	for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	  {
	    /* determine min/max d we're allowing for this state v and this position j */
	    if(do_banded) { 
	      dn = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	      dx = ESL_MIN((j-i0+1), dmax[v]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { 
	      dn = (cm->sttype[v] == MP_st) ? 2 : 1;
	      dx = ESL_MIN((j-i0+1), W); 
	    }

	    jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	    jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	    sd   = StateDelta(cm->sttype[v]);

	    if(cm->sttype[v] == B_st) {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      for (d = dn; d <= dx; d++) {
		/* k is the length of the right fragment */
		/* Careful, make sure k is consistent with bands in state w and state y. */
		if(do_banded) {
		  kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		  kmin = ESL_MAX(kmin, 0);
		  kmax = ESL_MIN(dmax[y], (d-dmin[w]));
		}
		else { kmin = 0; kmax = d; }

		alpha[v][jp_v][d] = ESL_MAX(IMPOSSIBLE, cm->endsc[v] + (cm->el_selfsc * (d - sd)));
		for (k = kmin; k <= kmax; k++) { 
		  jp_w = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      alpha[v][jp_v][d] = ESL_MAX(alpha[v][jp_v][d], (alpha[w][jp_w][d-k] + alpha[y][jp_y][k]));
		}
		vsc[v] = ESL_MAX(vsc[v], alpha[v][jp_v][d]);
	      }
	    }
	    else { /* if cm->sttype[v] != B_st */
	      for (d = dn; d <= dx; d++) {
		alpha[v][jp_v][d] = ESL_MAX (IMPOSSIBLE, (cm->endsc[v] + (cm->el_selfsc * (d - sd))));
		y = cm->cfirst[v];
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  alpha[v][jp_v][d] = ESL_MAX (alpha[v][jp_v][d], (alpha[y+yoffset][jp_y][d - sd] + cm->tsc[v][yoffset]));
		
		/* add in emission score, if any */
		i = j-d+1;
		switch (cm->sttype[v]) {
		case MP_st: 
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][jp_v][d] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  break;
		case ML_st:
		case IL_st:
		  alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  break;
		case MR_st:
		case IR_st:
		  alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  break;
		} /* end of switch */
		vsc[v] = ESL_MAX(vsc[v], alpha[v][jp_v][d]);
	      } /* end of d = dn; d <= dx; d++ */
	    } /* end of else (v != B_st) */
	  } /*loop over decks v>0 */

	/* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
	 * 
	 * If local begins are off, the hit must be rooted at v=0.
	 * With local begins on, the hit is rooted at the second state in
	 * the traceback (e.g. after 0), the internal entry point. Divide & conquer
	 * can only handle this if it's a non-insert state; this is guaranteed
	 * by the way local alignment is parameterized (other transitions are
	 * -INFTY), which is probably a little too fragile of a method. 
	 */

	/* determine min/max d we're allowing for the root state and this position j */
	if(do_banded) { 
	  dn = ESL_MAX(dmin[0], 1); 
	  dx = ESL_MIN((j-i0+1), dmax[0]); 
	  dx = ESL_MIN(dx, W);
	}
	else { 
	  dn = 1; 
	  dx = ESL_MIN((j-i0+1), W); 
	}
	jp_v = cur;
	jp_y = cur;
	for (d = dn; d <= dx; d++) {
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = ESL_MAX(IMPOSSIBLE, alpha[y][cur][d] + cm->tsc[0][0]);
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	    alpha[0][cur][d] = ESL_MAX (alpha[0][cur][d], (alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]));
	  vsc[0] = ESL_MAX(vsc[0], alpha[0][cur][d]);
	}
	
	if (cm->flags & CM_LOCAL_BEGIN) {
	  for (y = 1; y < cm->M; y++) {
	    if(do_banded) {
	      dn = (cm->sttype[y] == MP_st) ? ESL_MAX(dmin[y], 2) : ESL_MAX(dmin[y], 1); 
	      dn = ESL_MAX(dn, dmin[0]);
	      dx = ESL_MIN((j-i0+1), dmax[y]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { 
	      dn = 1; 
	      dx = ESL_MIN((j-i0+1), W); 
	    }
	    jp_y = (cm->stid[y] == BEGL_S) ? (j % (W+1)) : cur;
	    for (d = dn; d <= dx; d++) {
	      alpha[0][cur][d] = ESL_MAX(alpha[0][cur][d], alpha[y][jp_y][d] + cm->beginsc[y]);
	      vsc[0] = ESL_MAX(vsc[0], alpha[0][cur][d]);
	    }
	  }
	}
      } /* end loop over end positions j */
    /* free alpha, we only care about vsc 
     */
    for (v = 0; v < cm->M; v++) 
      {
	if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	  for (j = 0; j <= W; j++) free(alpha[v][j]);
	  free(alpha[v]);
	} else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	  continue;
	} else {
	  free(alpha[v][0]);
	  free(alpha[v][1]);
	  free(alpha[v]);
	}
      }
    free(alpha);
  }
  /*********************
   * end of if(do_cyk) *
   *********************/
  else { /* ! do_cyk, do_inside, with scaled int log odds scores instead of floats */

    ESL_ALLOC(ivsc, sizeof(int) * cm->M);
    esl_vec_FSet(ivsc, cm->M, -INFTY);
    
    /* ialpha allocations. (see comments for do_cyk section */ 
    ESL_ALLOC(ialpha, sizeof(int **) * cm->M);
    for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(ialpha[v], sizeof(int *) * (W+1));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(ialpha[v][j], sizeof(int) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      ialpha[v] = ialpha[cm->M-1];
    else 
      {
	ESL_ALLOC(ialpha[v], sizeof(int *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(ialpha[v][j], sizeof(int) * (W+1));
      }
    }
    /* ialpha initializations. (see comments for do_cyk section */
    for (v = cm->M-1; v >= 0; v--)  {
	ialpha[v][0][0] = -INFTY;
	if      (cm->sttype[v] == E_st)  ialpha[v][0][0] = 0;
	else if (cm->sttype[v] == MP_st) ialpha[v][0][1] = ialpha[v][1][1] = -INFTY;
	else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	  {
	    y = cm->cfirst[v];
	    ialpha[v][0][0] = cm->iendsc[v];
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	      ialpha[v][0][0] = ILogsum(ialpha[v][0][0], (ialpha[y+yoffset][0][0] 
							+ cm->itsc[v][yoffset]));
	    /* ...we don't bother to look at local alignment starts here... */
	    /* ! */
	    if (ialpha[v][0][0] < -INFTY) ialpha[v][0][0] = -INFTY;	
	  }
	else if (cm->sttype[v] == B_st)  {
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  ialpha[v][0][0] = ialpha[w][0][0] + ialpha[y][0][0]; 
	}
      ialpha[v][1][0] = ialpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  ialpha[v][j][0] = ialpha[v][0][0];
    }
    /* Impose the bands.
     *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
     *    the first loop will be skipped, the second initializes the E states.)
     */
    if(do_banded) {
      for (v = 0; v < cm->M; v++) {
	if(cm->stid[v] == BEGL_S) jmax = W; 
	else jmax = 1;
	
	dx = ESL_MIN(dmin[v], W);
	for (d = 0; d < dx; d++) 
	  for(j = 0; j <= jmax; j++)
	    ialpha[v][j][d] = -INFTY;
	
	for (d = dmax[v]+1; d <= W;      d++) 
	  for(j = 0; j <= jmax; j++)
	    ialpha[v][j][d] = -INFTY;
      }
    }

    /* The main loop: scan the sequence from position i0 to j0.
     */
    for (j = i0; j <= j0; j++) 
      {
	cur = j%2;
	prv = (j-1)%2;
	for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	  {
	    /* determine min/max d we're allowing for this state v and this position j */
	    if(do_banded) { 
	      dn = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	      dx = ESL_MIN((j-i0+1), dmax[v]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { 
	      dn = (cm->sttype[v] == MP_st) ? 2 : 1;
	      dx = ESL_MIN((j-i0+1), W); 
	    }

	    jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	    jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	    sd   = StateDelta(cm->sttype[v]);

	    if(cm->sttype[v] == B_st) {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      for (d = dn; d <= dx; d++) {
		/* k is the length of the right fragment */
		/* Careful, make sure k is consistent with bands in state w and state y. */
		if(do_banded) {
		  kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		  kmin = ESL_MAX(kmin, 0);
		  kmax = ESL_MIN(dmax[y], (d-dmin[w]));
		}
		else { kmin = 0; kmax = d; }

		ialpha[v][jp_v][d] = ESL_MAX(-INFTY, cm->iendsc[v] + (cm->iel_selfsc * (d - sd)));
		for (k = kmin; k <= kmax; k++) { 
		  jp_w = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      ialpha[v][jp_v][d] = ESL_MAX(ialpha[v][jp_v][d], (ialpha[w][jp_w][d-k] + ialpha[y][jp_y][k]));
		}
		ivsc[v] = ESL_MAX(ivsc[v], ialpha[v][jp_v][d]);
	      }
	    }
	    else { /* if cm->sttype[v] != B_st */
	      for (d = dn; d <= dx; d++) {
		ialpha[v][jp_v][d] = ESL_MAX (-INFTY, (cm->iendsc[v] + (cm->iel_selfsc * (d - sd))));
		y = cm->cfirst[v];
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  ialpha[v][jp_v][d] = ESL_MAX (ialpha[v][jp_v][d], (ialpha[y+yoffset][jp_y][d - sd] + cm->itsc[v][yoffset]));
		
		/* add in emission score, if any */
		i = j-d+1;
		switch (cm->sttype[v]) {
		case MP_st: 
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    ialpha[v][jp_v][d] += cm->iesc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		  else
		    ialpha[v][cur][d] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);
		  break;
		case ML_st:
		case IL_st:
		  ialpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);
		  break;
		case MR_st:
		case IR_st:
		  ialpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);
		  break;
		} /* end of switch */
		ivsc[v] = ESL_MAX(ivsc[v], ialpha[v][jp_v][d]);
	      } /* end of d = dn; d <= dx; d++ */
	    } /* end of else (v != B_st) */
	  } /*loop over decks v>0 */
	/* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
	 * 
	 * If local begins are off, the hit must be rooted at v=0.
	 * With local begins on, the hit is rooted at the second state in
	 * the traceback (e.g. after 0), the internal entry point. Divide & conquer
	 * can only handle this if it's a non-insert state; this is guaranteed
	 * by the way local alignment is parameterized (other transitions are
	 * -INFTY), which is probably a little too fragile of a method. 
	 */

	/* determine min/max d we're allowing for the root state and this position j */
	if(do_banded) { 
	  dn = ESL_MAX(dmin[0], 1); 
	  dx = ESL_MIN((j-i0+1), dmax[0]); 
	  dx = ESL_MIN(dx, W);
	}
	else { 
	  dn = 1; 
	  dx = ESL_MIN((j-i0+1), W); 
	}
	jp_v = cur;
	jp_y = cur;
	for (d = dn; d <= dx; d++) {
	  y = cm->cfirst[0];
	  ialpha[0][cur][d] = ESL_MAX(IMPOSSIBLE, ialpha[y][cur][d] + cm->itsc[0][0]);
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	    ialpha[0][cur][d] = ESL_MAX (ialpha[0][cur][d], (ialpha[y+yoffset][cur][d] + cm->itsc[0][yoffset]));
	  ivsc[0] = ESL_MAX(ivsc[0], ialpha[0][cur][d]);
	}
	
	if (cm->flags & CM_LOCAL_BEGIN) {
	  for (y = 1; y < cm->M; y++) {
	    if(do_banded) {
	      dn = ESL_MAX(1,  dmin[y]);
	      dn = ESL_MAX(dn, dmin[0]);
	      dx = ESL_MIN((j-i0+1), dmax[y]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { dn = 1; dx = W; }
	    jp_y = (cm->stid[y] == BEGL_S) ? (j % (W+1)) : cur;
	    for (d = dn; d <= dx; d++) {
	      ialpha[0][cur][d] = ESL_MAX(ialpha[0][cur][d], ialpha[y][jp_y][d] + cm->ibeginsc[y]);
	      ivsc[0] = ESL_MAX(ivsc[0], ialpha[0][cur][d]);
	    }
	  }
	}
      } /* end loop over end positions j */
    /* free ialpha, we only care about ivsc 
     */
    for (v = 0; v < cm->M; v++) 
      {
	if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	  for (j = 0; j <= W; j++) free(ialpha[v][j]);
	  free(ialpha[v]);
	} else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	  continue;
	} else {
	  free(ialpha[v][0]);
	  free(ialpha[v][1]);
	  free(ialpha[v]);
	}
      }
    free(ialpha);
    /* convert ivsc to floats in vsc */
    ESL_ALLOC(vsc, sizeof(float) * cm->M);
    for(v = 0; v < cm->M; v++)
      vsc[v] = Scorify(ivsc[v]);
    free(ivsc);
  }
  /**************************
   * end of else (do_inside)
   **************************/

  ret_val = vsc[0];
  if (ret_vsc != NULL) *ret_vsc = vsc;
  else free(vsc);
  
  return ret_val;

  ERROR:
    cm_Fail("Memory allocation error.\n");
    return 0.; /* NEVERREACHED */
}
#endif


/* cmcalibrate.c
 * Score a CM and a CM Plan 9 HMM against random sequence 
 * data to set the statistical parameters for E-value determination,
 * and CP9 HMM filtering thresholds. 
 * 
 * EPN, Wed May  2 07:02:52 2007
 * based on HMMER-2.3.2's hmmcalibrate.c from SRE
 *
 * MPI example:  
 * qsub -N testrun -j y -R y -b y -cwd -V -pe lam-mpi-tight 32 'mpirun C ./mpi-cmcalibrate foo.cm > foo.out'
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

#define MPI_FINISHED_GUMBEL     -1 /* message to send to workers */
#define MPI_FINISHED_CP9_FILTER -2 /* message to send to workers */

#include "funcs.h"		/* external functions                   */
#include "structs.h"

#define CUTOPTS  "--filE,--ga,--nc,--tc,--all"  /* Exclusive choice for filter threshold score cutoff */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-s",        eslARG_INT,    NULL,  NULL, "n>0",     NULL,      NULL,        NULL, "set random number generator seed to <n>",  1 },
  { "-t",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print timings for Gumbel fitting and CP9 filter calculation",  1},
#ifdef HAVE_DEVOPTS 
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print arguably interesting info",  1},
#endif
  { "--id",        eslARG_NONE, FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use same seqs for all Gumbel fittings", 1},
  /* 4 --p* options below are hopefully temporary b/c if we have E-values for the CM using a certain cm->pbegin, cm->pend,
   * changing those values in cmsearch invalidates the E-values, so we should pick hard-coded values for cm->pbegin cm->pend */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--pbegin","set all local begins as equiprobable", 1 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      NULL,    "--pend", "set all local end probs to <x>", 1 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        NULL, "set aggregate local begin prob to <x>", 1 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        NULL, "set aggregate local end prob to <x>", 1 },

  /* options for gumbel estimation */
  { "--gumonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--filonly", "only estimate Gumbels, don't calculate filter thresholds", 2},
  { "--gumL",    eslARG_INT,     NULL, NULL, "n>0",     NULL,      NULL, "--filonly", "set length of random seqs for Gumbel estimation to <n>", 2},
  { "--cmN",     eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CM gumbel estimation",    2 },
  { "--hmmN",    eslARG_INT,   "5000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CP9 HMM gumbel estimation",    2 },
  { "--gcfromdb",eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL, "--filonly", "use GC content distribution from file <s>",  2},
  { "--pfile",   eslARG_INFILE,  NULL, NULL, NULL,      NULL,"--gcfromdb",        NULL, "read partition info for Gumbels from file <s>", 2},
  { "--gumhfile",eslARG_OUTFILE,  NULL, NULL, NULL,     NULL,      NULL, "--filonly", "save fitted Gumbel histogram(s) to file <s>", 2 },
  { "--gumqqfile",eslARG_OUTFILE, NULL, NULL, NULL,     NULL,      NULL, "--filonly", "save Q-Q plot for Gumbel histogram(s) to file <s>", 2 },
  { "--beta",    eslARG_REAL,  "1e-7", NULL, "x>0",     NULL,      NULL,    "--noqdb", "set tail loss prob for Gumbel calculation to <x>", 5 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "DO NOT use query dependent banding (QDB) Gumbel searches", 5 },
  { "--gtail",   eslARG_REAL,    "0.5",NULL, "0.0<x<0.6",NULL,     NULL,        NULL, "set fraction of right histogram tail to fit to Gumbel to <x>", 5 },
  /* options for filter threshold calculation */
  { "--filonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "only calculate filter thresholds, don't estimate Gumbels", 3},
  { "--filN",    eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--gumonly", "number of emitted sequences for HMM filter threshold calc",    3 },
  { "--F",       eslARG_REAL,  "0.99", NULL, "0<x<=1",  NULL,      NULL, "--gumonly", "required fraction of seqs that survive HMM filter", 3},
  { "--db",      eslARG_INFILE,  NULL, NULL, NULL,      NULL,  "--filE", "--gumonly", "with --filE, set database size to size of <f>, not 1 Mb",  2},
  { "--starg",   eslARG_REAL,  "0.01", NULL, "0<x<=1",  NULL,      NULL, "--gumonly", "target filter survival fraction", 3},
  { "--gemit",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "when calc'ing filter thresholds, always emit globally from CM",  3},
  { "--fviterbi",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly,--fforward", "always choose CP9 Viterbi filter over CP9 Forward filter",  3},
  { "--fforward",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly,--fviterbi", "always choose CP9 Forward filter over CP9 Viterbi filter",  3},
  { "--filhfile",eslARG_OUTFILE, NULL,  NULL, NULL,     NULL,      NULL, "--gumonly", "save CP9 filter threshold histogram(s) to file <s>", 3},
  { "--filrfile",eslARG_OUTFILE, NULL,  NULL, NULL,     NULL,      NULL, "--gumonly", "save CP9 filter threshold information in R format to file <s>", 3},
  #ifdef HAVE_DEVOPTS
  { "--hybrid",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--filonly,--gumonly", "try hybrid filters, and keep them if better than HMM", 3},
  { "--filbeta", eslARG_REAL,  "1e-3",NULL, "x>0",     NULL,      NULL, "--gumonly", "set tail loss prob for filtering sub-CMs QDB to <x>", 5 },
#endif 
  /* exclusive choice of filter threshold cutoff */
  { "--filE",    eslARG_REAL,   "0.1", NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "min CM E-val (for a 1 Mb db) to consider for filter thr calc", 4}, 
  { "--ga",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM gathering threshold as minimum sc for filter thr calc", 4}, 
  { "--nc",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM noise cutoff as minimum sc for filter thr calc", 4}, 
  { "--tc",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM trusted cutoff as minimum sc for filter thr calc", 4},   
  { "--all",     eslARG_NONE,   FALSE, NULL, NULL,   CUTOPTS,      NULL, "--gumonly", "accept all CM hits for filter calc, DO NOT use cutoff", 4}, 
/* Other options */
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb", 5 },  
  { "--mxsize",  eslARG_REAL, "256.0", NULL, "x>0.",    NULL,      NULL,        NULL, "set maximum allowable HMM banded DP matrix size to <x> Mb", 9 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program", 5 },  
#endif

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char            *cmfile;	      /* name of input CM file  */ 
  ESL_RANDOMNESS  *r;
  ESL_ALPHABET    *abc;
  ESL_STOPWATCH   *w_stage;           /* stopwatch for each gumbel, filter stage */
  double          *gc_freq;;
  double          *pgc_freq;
  CMStats_t      **cmstatsA;          /* the CM stats data structures, 1 for each CM */
  HybridScanInfo_t *hsi;              /* information for a hybrid scan */ 
  int              ncm;                /* what number CM we're on */
  int              be_verbose;	       /* print extra info, only can be TRUE if #ifdef HAVE_DEVOPTS */
  int              cmalloc;            /* number of cmstats we have allocated */
  char            *tmpfile;            /* tmp file we're writing to */
  char            *mode;               /* write mode, "w" or "wb"                     */
  long             dbsize;             /* size of DB for gumbel stats (impt for E-value cutoffs for filters) */ 
  int              np;                 /* number of partitions, must be 1 unless --pfile or --filonly invoked,
					* once set, is never changed, once case we look out for: 
					* if --filonly invoked and CM file has CMs with diff number of partitions,
					* this is unlikely but possible, in this case we die in initialize_cmstats() 
					*/
  int             *pstart;             /* [0..p..np-1], begin points for partitions, end pts are implicit */
  float           *avglen;             /* [0..v..M-1] average hit len for subtree rooted at each state v for current CM */

  /* info for the comlog we'll add to the cmfiles */
  char            *ccom;               /* command line used in this execution of cmcalibrate */
  char            *cdate;              /* date of this execution of cmcalibrate */


  /* the following data is modified for each CM, and in some cases for each Gumbel mode for each CM,
   * it is assumed to be 'current' in many functions.
   */
  float           *cutoffA;            /* bit score cutoff for each partition, changes to reflect
				        * current mode CM is in, on masters and workers */
  float           *full_vcalcs;        /* [0..v..cm->M-1] millions of calcs for each subtree to scan 1 residue with --beta  */
  double         **vmuAA;              /* [0..np-1][0..cm->M-1], mu for each partition, each state, 
				        * if vmuAA[p][v] == -1 : we're not fitting state v to a gumbel */
  double         **vlambdaAA;          /* same as vmuAA, but lambda */
  GumbelInfo_t   **gum_hybA;           /* [0..np-1], hybrid gumbel info for each partition, rewritten 
					* for each candidate set of hybrid sub cm roots */
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
static char banner[] = "fit Gumbels for E-values and calculate HMM filter thresholds";

static int init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
#ifdef HAVE_MPI
static int init_worker_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
#endif

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
				   float ***ret_vscAA, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, float **ret_hyb_scA, int **ret_partA);
static int process_gumbel_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq, int L,
				   float ***ret_vscAA, float **ret_cp9scA, float **ret_hybscA);

static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);

static int update_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int fthr_mode);
static int set_partition_gc_freq(struct cfg_s *cfg, int p);
static int fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, double *ret_mu, double *ret_lambda);
static int cm_fit_histograms(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **vscA, int nscores, int p);
static int get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, double *dnull, int L, ESL_DSQ **ret_dsq);
static int get_cmemit_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int *ret_L, int *ret_p, Parsetree_t **ret_tr, ESL_DSQ **ret_dsq);
static int cm_find_hit_above_cutoff(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq, Parsetree_t *tr, int L, float cutoff, float *ret_sc);
static void estimate_workunit_time(const ESL_GETOPTS *go, const struct cfg_s *cfg, int nseq, int L, int gum_mode);
static int read_partition_file(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int update_avg_hit_len(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int switch_global_to_local(CM_t *cm, char *errbuf);
static int predict_cp9_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_vit_cp9scA, float *fil_fwd_cp9scA, int *fil_partA, BestFilterInfo_t *bf);
static int get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm);
static void format_time_string(char *buf, double sec, int do_frac);
static int print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
#ifdef HAVE_DEVOPTS
static int predict_best_sub_cm_roots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **fil_vscAA, int **ret_best_sub_roots);
static int predict_hybrid_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_hybscA, int *fil_partA, GumbelInfo_t **gum_hybA, BestFilterInfo_t *bf, int *ret_getting_faster);
#endif

int
main(int argc, char **argv)
{

  ESL_GETOPTS     *go	   = NULL;     /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory allocation error, stopwatch could not be created.");
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
  cfg.ncm      = 0;
  cfg.cmstatsA = NULL;
  cfg.hsi      = NULL;
  cfg.tmpfile  = NULL;
  cfg.mode     = NULL;
  cfg.dbsize   = 1000000; /* default DB size, changed if --db */
  cfg.cutoffA  = NULL; 
  cfg.full_vcalcs = NULL;
  cfg.vmuAA     = NULL;
  cfg.vlambdaAA = NULL;
  cfg.gum_hybA  = NULL;
  cfg.np        = 1;     /* default number of partitions is 1, changed if --pfile */
  cfg.pstart    = NULL;  /* allocated (by default to size 1) in init_master_cfg() */
  cfg.ccom      = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */
  cfg.cdate     = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */

  cfg.gumhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gumqfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.filhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.filrfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.abc      = NULL; 
  cfg.w_stage  = NULL; 

  cfg.do_mpi   = FALSE;
  cfg.my_rank  = 0;
  cfg.nproc    = 0;
  cfg.do_stall = esl_opt_GetBoolean(go, "--stall");
  cfg.be_verbose = FALSE;

  ESL_DASSERT1((GUM_CP9_GV == 0));
  ESL_DASSERT1((GUM_CP9_GF == 1));
  ESL_DASSERT1((GUM_CM_GC  == 2));
  ESL_DASSERT1((GUM_CM_GI  == 3));
  ESL_DASSERT1((GUM_CP9_LV == 4));
  ESL_DASSERT1((GUM_CP9_LF == 5));
  ESL_DASSERT1((GUM_CM_LC  == 6));
  ESL_DASSERT1((GUM_CM_LI  == 7));
  ESL_DASSERT1((GUM_NMODES == 8));
  ESL_DASSERT1((FTHR_CM_GC == 0));
  ESL_DASSERT1((FTHR_CM_GI == 1));
  ESL_DASSERT1((FTHR_CM_LC == 2));
  ESL_DASSERT1((FTHR_CM_LI == 3));
  ESL_DASSERT1((FTHR_NMODES== 4));

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

  if(cfg.my_rank == 0) { /* master, serial or mpi */
    /* Rewind the CM file for a second pass.
     * Write a temporary CM file with new stats information in it
     */
    int   status;
    int   cmi;
    CM_t *cm;
    FILE *outfp;
    sigset_t blocksigs;  /* list of signals to protect from             */
    char     errbuf[cmERRBUFSIZE];

    CMFileRewind(cfg.cmfp);
    if (esl_FileExists(cfg.tmpfile))                    cm_Fail("Ouch. Temporary file %s appeared during the run.", cfg.tmpfile);
    if ((outfp = fopen(cfg.tmpfile, cfg.mode)) == NULL) cm_Fail("Ouch. Temporary file %s couldn't be opened for writing.", cfg.tmpfile); 
    
    for (cmi = 0; cmi < cfg.ncm; cmi++) {
      if (!CMFileRead(cfg.cmfp, &(cfg.abc), &cm)) cm_Fail("Ran out of CMs too early in pass 2");
      if (cm == NULL)                             cm_Fail("CM file %s was corrupted? Parse failed in pass 2", cfg.cmfile);

      /* update the cm->comlog info */
      if((status = update_comlog(go, errbuf, cfg.ccom, cfg.cdate, cm)) != eslOK) cm_Fail(errbuf);
	
      cm->stats = cfg.cmstatsA[cmi];
      cm->flags &= ~CMH_FILTER_STATS; /* forget that CM may have had filter stats, if --gumonly was invoked, it won't anymore */

      if(!(esl_opt_GetBoolean(go, "--filonly"))) cm->flags |= CMH_GUMBEL_STATS; 
      if(!(esl_opt_GetBoolean(go, "--gumonly"))) cm->flags |= CMH_FILTER_STATS; 
      if((status = CMFileWrite(outfp, cm, cfg.cmfp->is_binary, errbuf)) != eslOK) cm_Fail(go->errbuf);
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
    
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
    
    /* master specific cleaning */
    if (cfg.gumhfp   != NULL) fclose(cfg.gumhfp);
    if (cfg.gumqfp   != NULL) fclose(cfg.gumqfp);
    if (cfg.filhfp   != NULL) fclose(cfg.filhfp);
    if (cfg.filrfp   != NULL) fclose(cfg.filrfp);
    if (cfg.ccom     != NULL) free(cfg.ccom);
    if (cfg.cdate    != NULL) free(cfg.cdate);

  }
  /* clean up */
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Allocates/sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->gumhfp      - optional output file
 *    cfg->gumqfp      - optional output file
 *    cfg->filhfp      - optional output file
 *    cfg->filrfp      - optional output file
 *    cfg->gc_freq     - observed GC freqs (if --gcfromdb invoked)
 *    cfg->cmstatsA    - the stats, allocated only
 *    cfg->np          - number of partitions, never changes once set 
 *    cfg->pstart      - array of partition starts 
 *    cfg->r           - source of randomness
 *    cfg->tmpfile     - temp file for rewriting cm file
 *    cfg->be_verbose  - print extra info? (only available #ifdef HAVE_DEVOPTS) 
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
  if(esl_opt_GetString(go, "--gcfromdb") != NULL) {
    ESL_ALPHABET *tmp_abc = NULL;
    tmp_abc = esl_alphabet_Create(eslRNA);
    ESL_SQFILE      *dbfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--gcfromdb"), eslSQFILE_UNKNOWN, NULL, &dbfp);
    if (status == eslENOTFOUND)    cm_Fail("No such file."); 
    else if (status == eslEFORMAT) cm_Fail("Format unrecognized."); 
    else if (status == eslEINVAL)  cm_Fail("Can’t autodetect stdin or .gz."); 
    else if (status != eslOK)      cm_Fail("Failed to open sequence database file, code %d.", status); 
    GetDBInfo(tmp_abc, dbfp, NULL, &(cfg->gc_freq)); 
    esl_vec_DNorm(cfg->gc_freq, GC_SEGMENTS);
    esl_alphabet_Destroy(cfg->abc);
    esl_sqfile_Close(dbfp); 
   /* allocate pgc_freq, the gc freqs per partition, used to sample seqs for different partitions */
    ESL_ALLOC(cfg->pgc_freq, sizeof(double) * GC_SEGMENTS);
  }

  /* optionally, set cfg->dbsize as size of db in file <f> from --db <f> */
  if(esl_opt_GetString(go, "--db") != NULL) {
    ESL_SQFILE      *dbfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--gcfromdb"), eslSQFILE_UNKNOWN, NULL, &dbfp);
    if (status == eslENOTFOUND)    cm_Fail("No such file."); 
    else if (status == eslEFORMAT) cm_Fail("Format unrecognized."); 
    else if (status == eslEINVAL)  cm_Fail("Can’t autodetect stdin or .gz."); 
    else if (status != eslOK)      cm_Fail("Failed to open sequence database file, code %d.", status); 
    GetDBInfo(NULL, dbfp, &(cfg->dbsize), NULL);  
    esl_sqfile_Close(dbfp); 
  }

  /* set up the partition data that's used for all CMs */
  if(esl_opt_IsDefault(go, "--pfile")) { /* by default we have 1 partition 0..100 */
    ESL_ALLOC(cfg->pstart, sizeof(int) * 1);
    cfg->np        = 1;
    cfg->pstart[0] = 0;
  }
  else { /* setup cfg->np and cfg->pstart in read_partition_file() */
    if((status = read_partition_file(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  }

  cfg->be_verbose = FALSE;
#ifdef HAVE_DEVOPTS
  if (esl_opt_GetBoolean(go, "-v")) cfg->be_verbose = TRUE;        
#endif

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cfg->cmalloc  = 128;
  ESL_ALLOC(cfg->cmstatsA, sizeof(CMStats_t *) * cfg->cmalloc);
  cfg->ncm      = 0;

  /* seed master's RNG */
  if (! esl_opt_IsDefault(go, "-s")) 
    cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
  else cfg->r = esl_randomness_CreateTimeseeded();

  /* create the stopwatch */
  cfg->w_stage = esl_stopwatch_Create();

  printf("Random number generator seed: %ld\n", esl_randomness_GetSeed(cfg->r));

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

  cfg->avglen = NULL; /* this will be allocated and filled inside serial_master() or mpi_master() */
  if(cfg->r       == NULL) cm_Fail("Failed to create master RNG.");
  if(cfg->w_stage == NULL) cm_Fail("Failed to create stopwatch.");

  /* fill cfg->ccom, and cfg->cdate */
  if((status = get_cmcalibrate_comlog_info(go, cfg, errbuf)) != eslOK) return status;

  return eslOK;

 ERROR:
  return status;
}

#ifdef HAVE_MPI
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
  int status; 
  cfg->cmfile = NULL;
  cfg->abc      = NULL;
  cfg->gc_freq  = NULL;
  cfg->pgc_freq = NULL;
  cfg->be_verbose = FALSE;
  cfg->tmpfile  = NULL;
  cfg->mode     = NULL;
  cfg->dbsize  = 0;
  cfg->cutoffA = NULL;
  cfg->full_vcalcs = NULL;
  cfg->vmuAA  = NULL;
  cfg->vlambdaAA = NULL;
  cfg->gum_hybA  = NULL;
  cfg->avglen = NULL;
  
  cfg->cmfp = NULL;
  cfg->gumhfp = NULL;
  cfg->gumqfp = NULL;
  cfg->filhfp = NULL;
  cfg->filrfp = NULL;
  
  /* allocate cmstats results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cfg->cmalloc  = 128;
  ESL_ALLOC(cfg->cmstatsA, sizeof(CMStats_t *) * cfg->cmalloc);
  cfg->ncm      = 0;

  /* we may receive info pertaining to these from master 
   * inside mpi_worker(), at which time we'll update them
   */
  cfg->np       = 0;  
  cfg->pstart   = NULL;
  cfg->r        = NULL; /* we'll create this when seed is received from master */

  return eslOK;

 ERROR:
  return status;
}
#endif /* #ifdef HAVE_MPI */

/* serial_master()
 * The serial version of cmcalibrate.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[cmERRBUFSIZE];
  CM_t    *cm = NULL;
  int      cmN  = esl_opt_GetInteger(go, "--cmN");
  int      hmmN = esl_opt_GetInteger(go, "--hmmN");
  int      filN = esl_opt_GetInteger(go, "--filN");
  int      gum_mode  = 0;
  int      fthr_mode = 0;
  int      p;
  int      cmi;
  float  **gum_vscAA      = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each random seq */
  float  **fil_vscAA      = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each emitted seq */
  float   *gum_cp9scA     = NULL; /*                [0..nseq-1] best cp9 score for each random seq */
  float   *fil_vit_cp9scA = NULL; /*                [0..nseq-1] best cp9 Viterbi score for each emitted seq */
  float   *fil_fwd_cp9scA = NULL; /*                [0..nseq-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA      = NULL; /*                [0..nseq-1] partition of CM emitted seq */
  double   tmp_mu, tmp_lambda;    /* temporary mu and lambda used for setting HMM gumbels */
  int      L;                     /* length of sequences to search for gumbel fitting, L==cm->W*2 unless --gumL enabled, 
				   * in which case L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL") */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  void    *tmp;
  long     seed;
#ifdef HAVE_DEVOPTS
  /* variables for --hybrid */
  float   *gum_hybscA     = NULL; /*                [0..nseq-1] best hybrid score for each random seq */
  float   *fil_hybscA     = NULL; /*                [0..nseq-1] best hybrid score for each emitted seq */
  int     *best_sub_roots = NULL; /* [0..cfg->hsi->nstarts-1], best predicted sub CM root filter for each start group */
  int      s = 0;
  int      getting_faster = TRUE;
#endif

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->cmstatsA, tmp, sizeof(CMStats_t *) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;
      if (esl_opt_GetBoolean(go, "--filonly") && (! (cm->flags & CMH_GUMBEL_STATS))) cm_Fail("--filonly invoked, but CM %s (CM number %d) does not have Gumbel stats in CM file\n", cm->name, (cmi+1));

      if((status = initialize_cm(go, cfg, errbuf, cm))      != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = print_cm_info     (go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);

      for(gum_mode = 0; gum_mode < GUM_NMODES; gum_mode++) {

#ifdef HAVE_DEVOPTS
	/* free and recalculate hybrid scan info, b/c when investigating hybrid filters we may have added sub CM roots */
	if(cfg->hsi != NULL) cm_FreeHybridScanInfo(cfg->hsi, cm);
	cfg->hsi = cm_CreateHybridScanInfo(cm, esl_opt_GetReal(go, "--filbeta"), cfg->full_vcalcs[0]);
#endif
	
	/* do we need to switch from glocal configuration to local? */
	if(gum_mode > 0 && (! GumModeIsLocal(gum_mode-1)) && GumModeIsLocal(gum_mode)) {
	  if((status = switch_global_to_local(cm, errbuf)) != eslOK)      cm_Fail(errbuf);
	  if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
	}
	/* update search opts for gumbel mode */
	GumModeToSearchOpts(cm, gum_mode);
	/* if --id, we free RNG, then create a new one and reseed it with the initial seed,
	 * so we use same seqs for Gumbel fittings of all CM modes and HMM modes 
	 * (but if --cmN <n> != --hmmN <n>) they won't be the same between CM and HMM modes. 
	 */
	if(esl_opt_GetBoolean(go, "--id")) {
	  seed = esl_randomness_GetSeed(cfg->r);
	  esl_randomness_Destroy(cfg->r);
	  cfg->r = esl_randomness_Create(seed);
	}

	/* gumbel fitting section */
	if(! (esl_opt_GetBoolean(go, "--filonly"))) {
	  /* calculate gumbels for this gum mode */
	  /* determine length of seqs to search for gumbel fitting */
	  if(esl_opt_IsDefault(go, "--gumL")) L = cm->W*2; 
	  else                                L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL")); /* minimum L we allow is 2 * cm->W, this is enforced silently (!) */

	  ESL_DASSERT1((cfg->np == cfg->cmstatsA[cmi]->np));
	  for (p = 0; p < cfg->np; p++) {
	    esl_stopwatch_Start(cfg->w_stage);
	    if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	    if(GumModeIsForCM(gum_mode)) { /* CM mode */
	      /* search random sequences to get gumbels for each candidate sub CM root state (including 0, the root of the full model) */
	      if(cfg->be_verbose) estimate_workunit_time(go, cfg, cmN, L, gum_mode);
	      ESL_DPRINTF1(("\n\ncalling process_gumbel_workunit to fit gumbel for p: %d CM mode: %d\n", p, gum_mode));
	      printf("gumbel  %-12s %5d %6d %5s %5s [", DescribeGumMode(gum_mode), cmN, L, "-", "-");
	      fflush(stdout);
	      if((status = process_gumbel_workunit (go, cfg, errbuf, cm, cmN, L, &gum_vscAA, NULL, NULL)) != eslOK) cm_Fail(errbuf);
	      if((status = cm_fit_histograms(go, cfg, errbuf, cm, gum_vscAA, cmN, p)) != eslOK) cm_Fail(errbuf);
	      SetGumbelInfo(cfg->cmstatsA[cmi]->gumAA[gum_mode][p], cfg->vmuAA[p][0], cfg->vlambdaAA[p][0], L, cmN);
	    }
	    else { /* CP9 mode, fit gumbel for full HMM */
	      if(cfg->be_verbose) estimate_workunit_time(go, cfg, hmmN, L, gum_mode);
	      ESL_DPRINTF1(("\n\ncalling process_gumbel_workunit to fit gumbel for p: %d CP9 mode: %d\n", p, gum_mode));
	      printf("gumbel  %-12s %5d %6d %5s %5s [", DescribeGumMode(gum_mode), hmmN, L, "-", "-");
	      fflush(stdout);
	      if((status = process_gumbel_workunit (go, cfg, errbuf, cm, hmmN, L, NULL, &gum_cp9scA, NULL)) != eslOK) cm_Fail(errbuf);
	      if((status = fit_histogram(go, cfg, errbuf, gum_cp9scA, hmmN, &tmp_mu, &tmp_lambda))       != eslOK) cm_Fail(errbuf);
	      SetGumbelInfo(cfg->cmstatsA[cmi]->gumAA[gum_mode][p], tmp_mu, tmp_lambda, L, hmmN);
	    }
	    esl_stopwatch_Stop(cfg->w_stage);
	    format_time_string(time_buf, cfg->w_stage->elapsed, 0);
	    printf(" %10s\n", time_buf);
	  } /* end of for loop over partitions */
	} /* end of if ! --filonly */
	
	/* filter threshold section */
	if(! (esl_opt_GetBoolean(go, "--gumonly"))) {
	  if(GumModeIsForCM(gum_mode)) { /* CM mode, we want to do filter threshold calculations */
	    esl_stopwatch_Start(cfg->w_stage);
	    fthr_mode = GumModeToFthrMode(gum_mode);
	    ESL_DASSERT1((fthr_mode != -1));
	    if((status = update_cutoffs(go, cfg, errbuf, cm, gum_mode)) != eslOK) cm_Fail(errbuf);
	  
	    /* search emitted sequences to get filter thresholds for HMM and each candidate sub CM root state */
	    ESL_DPRINTF1(("\n\ncalling process_filter_workunit to get HMM filter thresholds for p: %d mode: %d\n", p, gum_mode));
	    if(fil_partA != NULL) free(fil_partA);
	    ///printf("filter  %3s  %-8s %5s %6s %5d %5g [", "-", DescribeFthrMode(fthr_mode), "-", "-", filN, esl_opt_GetReal(go, "--filE"));
	    printf("filter  %3s  %-8s %5s %6s %5d ", "-", DescribeFthrMode(fthr_mode), "-", "-", filN);
	    if     (esl_opt_GetBoolean(go, "--all")) printf("infty [");
	    else if(esl_opt_GetBoolean(go, "--ga"))  printf("   GA [");
	    else if(esl_opt_GetBoolean(go, "--nc"))  printf("   NC [");
	    else if(esl_opt_GetBoolean(go, "--tc"))  printf("   TC [");
	    else                                     printf("%5g [", esl_opt_GetReal(go, "--filE"));
	    fflush(stdout);

#ifdef HAVE_DEVOPTS
	    if(esl_opt_GetBoolean(go, "--hybrid")) { /* we want fil_vscAA filled with hybrid scanning scores */
	      if((status = process_filter_workunit (go, cfg, errbuf, cm, filN, &fil_vscAA, &fil_vit_cp9scA, &fil_fwd_cp9scA, NULL, &fil_partA)) != eslOK) cm_Fail(errbuf);
	    }
	    else { /* we don't want fil_vscAA filled with hybrid scanning scores */
#endif
	      if((status = process_filter_workunit (go, cfg, errbuf, cm, filN, NULL, &fil_vit_cp9scA, &fil_fwd_cp9scA, NULL, &fil_partA)) != eslOK) cm_Fail(errbuf);
#ifdef HAVE_DEVOPTS
	    }	      
#endif
	    /* determine the 'best' way to filter using a heuristic:
	     * 1. predict speedup for HMM-only filter
	     * if --hybrid, continue to 2, else use HMM filter as 'best' filter.
	     * 2. add 'best' sub CM root independent of those already chosen, and predict speedup for a hybrid scan
	     * 3. repeat 2 until predicted speedup drops
	     * if --hybrid, best filter is either HMM only or fastest combo of sub CM roots
	     */

	    /* 1. predict speedup for HMM-only filter */
	    if((status = predict_cp9_filter_speedup(go, cfg, errbuf, cm, fil_vit_cp9scA, fil_fwd_cp9scA, fil_partA, cfg->cmstatsA[cmi]->bfA[fthr_mode])) != eslOK) cm_Fail(errbuf);
	    if(cfg->be_verbose) DumpBestFilterInfo(cfg->cmstatsA[cmi]->bfA[fthr_mode]);
	    
#ifdef HAVE_DEVOPTS	    
	    /* if desired, try hybrid filters to see if we can do better than the HMM */
	    if(esl_opt_GetBoolean(go, "--hybrid")) {
	      /* 2. predict best sub CM filter root state in each start group */
	      if((status = predict_best_sub_cm_roots(go, cfg, errbuf, cm, fil_vscAA, &best_sub_roots)) != eslOK) cm_Fail(errbuf);
	      
	      s = 0;
	      getting_faster = TRUE;
	      while(getting_faster && s < cfg->hsi->nstarts) { 
		/* add next fastest (predicted) sub CM root, if it's compatible, otherwise skip it (this is arguably a bad idea, but it's tough to remove roots in this case) */
		if(cm_CheckCompatibleWithHybridScanInfo(cm, cfg->hsi, best_sub_roots[s])) {
		  cm_AddRootToHybridScanInfo(cm, cfg->hsi, best_sub_roots[s++]);
		  /* fit gumbels for new hybrid scanner */
		  ESL_DASSERT1((cfg->np == cfg->cmstatsA[cmi]->np));
		  for (p = 0; p < cfg->np; p++) {
		    if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
		    if((status = process_gumbel_workunit (go, cfg, errbuf, cm, cmN, L, NULL, NULL, &gum_hybscA)) != eslOK) cm_Fail(errbuf);
		    if((status = fit_histogram(go, cfg, errbuf, gum_hybscA, cmN, &tmp_mu, &tmp_lambda))       != eslOK) cm_Fail(errbuf);
		    SetGumbelInfo(cfg->gum_hybA[p], tmp_mu, tmp_lambda, L, cmN);
		    if(cfg->be_verbose) debug_print_gumbelinfo(cfg->gum_hybA[p]);
		  }		
		  /* get hybrid scores of CM emitted seqs */
		  if(fil_partA != NULL) free(fil_partA);
		  if((status = process_filter_workunit (go, cfg, errbuf, cm, filN, NULL, NULL, NULL, &fil_hybscA, &fil_partA)) != eslOK) cm_Fail(errbuf);
		  /* determine predicted speedup of hybrid scanner */
		  if((status = predict_hybrid_filter_speedup(go, cfg, errbuf, cm, fil_hybscA, fil_partA, cfg->gum_hybA, cfg->cmstatsA[cmi]->bfA[fthr_mode], &getting_faster)) != eslOK) cm_Fail(errbuf);
		  DumpBestFilterInfo(cfg->cmstatsA[cmi]->bfA[fthr_mode]);
		}
		else { /* root we were going to add was incompatible, skip it */
		  s++;
		}
	      }
	    }
#endif
	    esl_stopwatch_Stop(cfg->w_stage);
	    format_time_string(time_buf, cfg->w_stage->elapsed, 0);
	    printf(" %10s\n", time_buf);
	  }
	}
	/* free muA and lambdaA */
      } /* end of for(gum_mode = 0; gum_mode < NCMMODES; gum_mode++) */
      if(cfg->be_verbose) debug_print_cmstats(cfg->cmstatsA[cmi], (! esl_opt_GetBoolean(go, "--gumonly")));
#ifdef HAVE_DEVOPTS
      if(cfg->hsi != NULL) cm_FreeHybridScanInfo(cfg->hsi, cm);
      cfg->hsi = NULL;
#endif
      FreeCM(cm);
      printf("//\n");
      fflush(stdout);
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
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include most worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
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
  void    *tmp;
  int      wi_error = 0;                /* worker index that sent back an error message, if an error occurs */

  CM_t           *cm  = NULL;
  int            cmN  = esl_opt_GetInteger(go, "--cmN");
  int            hmmN = esl_opt_GetInteger(go, "--hmmN");
  int            filN = esl_opt_GetInteger(go, "--filN");

  int        gum_mode = 0;
  int       fthr_mode = 0;
  int               p;
  float  **gum_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each random seq */
  float   *gum_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each random seq */

  float  **fil_vscAA      = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each emitted seq */
  float   *fil_vit_cp9scA = NULL; /*                [0..nseq-1] best viterbi cp9 score for each emitted seq */
  float   *fil_fwd_cp9scA = NULL; /*                [0..nseq-1] best forward cp9 score for each emitted seq */
  int     *fil_partA      = NULL; /*                [0..nseq-1] partition each CM emitted seq belongs to */

  /* data received from workers */
  float  **worker_vscAA = NULL;   
  float   *worker_cp9scA = NULL;
  float   *worker_vit_cp9scA = NULL;
  float   *worker_fwd_cp9scA = NULL;
  int     *worker_partA = NULL;

  long *seedlist = NULL;
  char  errbuf[cmERRBUFSIZE];
  MPI_Status mpistatus; 
  int   n, v, i;
  int working_on_cm;        /* TRUE when gum_mode is for CM gumbel */
  int working_on_cp9;       /* TRUE when gum_mode is for CP9 gumbel */
  int gum_nseq_per_worker  = 0; /* when calcing gumbels, number of seqs to tell each worker to work on */
  int gum_nseq_this_round  = 0; /* when calcing gumbels, number of seqs for current round */
  int fil_nseq_per_worker  = (filN / (cfg->nproc-1)); /* when calcing filters, number of seqs to tell each worker to work on */
  int nseq_sent        = 0; /* number of seqs we've told workers to work on */
  int nseq_this_worker = 0; /* number of seqs to tell current worker to work on */
  int nseq_just_recv   = 0; /* number of seqs we just received scores for from a worker */
  int nseq_recv        = 0; /* number of seqs we've received thus far this round from workers */
  int msg;
  int cmi;                  /* CM index, which number CM we're working on */
  double tmp_mu, tmp_lambda;/* temporary mu and lambda used for setting HMM gumbels */
  int            L;  /* length of sequences to search for gumbel fitting, L==cm->W*2 unless --gumL enabled, 
		      * in which case L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL") */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  float    update_i;

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL)    { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seedlist  = malloc(sizeof(long) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) cm_Fail(errbuf);
  ESL_DPRINTF1(("MPI master is initialized\n"));

  ESL_ALLOC(seedlist, sizeof(long) * cfg->nproc);
  for (wi = 0; wi < cfg->nproc; wi++) 
    {
      seedlist[wi] = esl_rnd_Choose(cfg->r, 1000000000); /* not sure what to use as max for seed */
      ESL_DPRINTF1(("wi %d seed: %ld\n", wi, seedlist[wi]));
    }

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

  /* 3 special (annoying) case:
   * case 1: if we've used the --gcfromdb option, we read in a seq file to fill
   * cfg->gc_freq, and we need to broadcast that info to workers
   *
   * case 2: if we are calculating stats for more than 1 partition, 
   * (--pfile invoked), we need to broadcast that information to 
   * the workers. 
   *
   * case 3: if we've changed the default dbsize for the CM E-value cutoff
   * in the filter threshold calculation (--db invoked), we need to broadcast 
   * that information to the workers.
   */
  if(! (esl_opt_IsDefault(go, "--gcfromdb"))) { /* receive gc_freq info from master */
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  if(! (esl_opt_IsDefault(go, "--pfile"))) { /* broadcast partition info to workers */
    ESL_DASSERT1((! (esl_opt_GetBoolean(go, "--filonly"))));
    MPI_Bcast(&(cfg->np),  1,       MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pstart != NULL));
    MPI_Bcast(cfg->pstart, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if(! (esl_opt_IsDefault(go, "--db"))) { /* receive cfg->dbsize info from master */
    MPI_Bcast(&(cfg->dbsize), 1, MPI_LONG, 0, MPI_COMM_WORLD);
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
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->cmstatsA, tmp, sizeof(CMStats_t *) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;
      if (esl_opt_GetBoolean(go, "--filonly") && (! (cm->flags & CMH_GUMBEL_STATS))) cm_Fail("--filonly invoked, but CM %s (CM number %d) does not have Gumbel stats in CM file\n", cm->name, (cmi+1));

      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm     (go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = print_cm_info     (go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);

      if(! (esl_opt_GetBoolean(go, "--filonly"))) { 
	ESL_ALLOC(gum_vscAA, sizeof(float *) * cm->M);
	for(v = 0; v < cm->M; v++) ESL_ALLOC(gum_vscAA[v], sizeof(float) * cmN);
	ESL_ALLOC(gum_cp9scA, sizeof(float)  * hmmN);
      }

      if(! (esl_opt_GetBoolean(go, "--gumonly"))) { 
	ESL_ALLOC(fil_vscAA, sizeof(float *) * cm->M);
	for(v = 0; v < cm->M; v++) ESL_ALLOC(fil_vscAA[v], sizeof(float) * filN);
	ESL_ALLOC(fil_vit_cp9scA, sizeof(float) * filN);
	ESL_ALLOC(fil_fwd_cp9scA, sizeof(float) * filN);
	ESL_ALLOC(fil_partA,      sizeof(int) *   filN);
      }

      for(gum_mode = 0; gum_mode < GUM_NMODES; gum_mode++) {

	if(GumModeIsForCM(gum_mode)) { working_on_cm = TRUE;  working_on_cp9 = FALSE; }
	else                         { working_on_cm = FALSE; working_on_cp9 = TRUE;  }
	gum_nseq_per_worker = working_on_cm ? (int) (cmN / (cfg->nproc-1)) : (int) (hmmN / (cfg->nproc-1));
	gum_nseq_this_round = working_on_cm ? cmN : hmmN;
	update_i = gum_nseq_this_round / 20.;

#ifdef HAVE_DEVOPTS
	/* free and recalculate hybrid scan info, b/c when investigating hybrid filters we may have added sub CM roots */
	if(cfg->hsi != NULL) cm_FreeHybridScanInfo(cfg->hsi, cm);
	cfg->hsi = cm_CreateHybridScanInfo(cm, esl_opt_GetReal(go, "--filbeta"), cfg->full_vcalcs[0]);
#endif

	/* do we need to switch from glocal configuration to local? */
	if(gum_mode > 0 && (! GumModeIsLocal(gum_mode-1)) && GumModeIsLocal(gum_mode)) {
	  if((status = switch_global_to_local(cm, errbuf)) != eslOK)      cm_Fail(errbuf);
	  if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
	}
	/* update search opts for gumbel mode */
	GumModeToSearchOpts(cm, gum_mode);

	/* gumbel fitting section */
	if(! (esl_opt_GetBoolean(go, "--filonly"))) {
	  if(esl_opt_IsDefault(go, "--gumL")) L = cm->W*2; 
	  else                                L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL")); /* minimum L we allow is 2 * cm->W, this is enforced silently (!) */

	  for (p = 0; p < cfg->np; p++) {
	    ESL_DPRINTF1(("MPI master: CM: %d gumbel mode: %d partition: %d\n", cfg->ncm, gum_mode, p));
	    
	    esl_stopwatch_Start(cfg->w_stage);
	    printf("gumbel  %-12s %5d %6d %5s %5s [", DescribeGumMode(gum_mode), gum_nseq_this_round, L, "-", "-");
	    fflush(stdout);

	    if(xstatus == eslOK) have_work     = TRUE;	/* TRUE while work remains  */
	    
	    wi = 1;
	    nseq_sent = 0;
	    nseq_recv = 0;
	    while (have_work || nproc_working)
	      {
		if(have_work) { 
		  if(nseq_sent < gum_nseq_this_round) {
		    nseq_this_worker = (nseq_sent + gum_nseq_per_worker <= gum_nseq_this_round) ? 
		      gum_nseq_per_worker : (gum_nseq_this_round - nseq_sent);
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
			  ESL_DPRINTF1(("MPI master sees that the result buffer contains calibration results\n"));
			  if(working_on_cm) {
			    if ((status = cmcalibrate_cm_gumbel_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, cm->M, &worker_vscAA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
			    ESL_DPRINTF1(("MPI master has unpacked CM gumbel results\n"));
			    ESL_DASSERT1((nseq_just_recv > 0));
			    for(v = 0; v < cm->M; v++) {
			      for(i = 0; i < nseq_just_recv; i++) {
				ESL_DPRINTF3(("\tscore from worker v: %d i: %d sc: %f\n", i, v, worker_vscAA[v][i]));
				gum_vscAA[v][nseq_recv+i] = worker_vscAA[v][i];
				if(nseq_recv+i > update_i) {
				  printf("=");
				  fflush(stdout); 
				  update_i += gum_nseq_this_round / 20.; 
				}
			      }
			      free(worker_vscAA[v]);
			    }
			    free(worker_vscAA);
			  }
			  else { /* working on cp9 */
			    if ((status = cmcalibrate_cp9_gumbel_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_cp9scA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
			    ESL_DPRINTF1(("MPI master has unpacked CP9 gumbel results\n"));
			    ESL_DASSERT1((nseq_just_recv > 0));
			    for(i = 0; i < nseq_just_recv; i++) {
			      gum_cp9scA[nseq_recv+i] = worker_cp9scA[i];
			      if(nseq_recv+i > update_i) {
				printf("=");
				fflush(stdout); 
				update_i += gum_nseq_this_round / 20.; 
			      }
			    }
			    free(worker_cp9scA);
			  }
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
		    /* send new search job */
		    ESL_DPRINTF1(("MPI master is sending nseq %d to worker %d\n", nseq_this_worker, wi));
		    MPI_Send(&(nseq_this_worker), 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	      
		    wi++;
		    nproc_working++;
		    nseq_sent += nseq_this_worker;
		  }
	      }

	    if(xstatus == eslOK) { 
	      /* fit gumbels for this partition p, this gumbel mode gum_mode */
	      if(working_on_cm) { 
		if((status = cm_fit_histograms(go, cfg, errbuf, cm, gum_vscAA, cmN, p)) != eslOK) cm_Fail(errbuf);
		SetGumbelInfo(cfg->cmstatsA[cmi]->gumAA[gum_mode][p], cfg->vmuAA[p][0], cfg->vlambdaAA[p][0], L, cmN);
	      }
	      else /* working on CP9 */ {
		if((status = fit_histogram(go, cfg, errbuf, gum_cp9scA, hmmN, &tmp_mu, &tmp_lambda))       != eslOK) cm_Fail(errbuf);
		SetGumbelInfo(cfg->cmstatsA[cmi]->gumAA[gum_mode][p], tmp_mu, tmp_lambda, L, hmmN);
	      }
	    }
	    esl_stopwatch_Stop(cfg->w_stage);
	    format_time_string(time_buf, cfg->w_stage->elapsed, 0);
	    printf("=] %10s\n", time_buf);
	  }
	  ESL_DPRINTF1(("MPI master: done with partition: %d for gumbel mode: %d for this CM. Telling all workers\n", p, gum_mode));
	  for (wi = 1; wi < cfg->nproc; wi++) { 
	    msg = MPI_FINISHED_GUMBEL;
	    MPI_Send(&msg, 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	  }
	} /* end of if ! --filonly */

	/* filter threshold section */
	if(GumModeIsForCM(gum_mode) && (! (esl_opt_GetBoolean(go, "--gumonly")))) {
	  fthr_mode = GumModeToFthrMode(gum_mode);
	  ESL_DASSERT1((fthr_mode != -1));
	  ESL_DPRINTF1(("MPI master: CM: %d fthr mode: %d\n", cfg->ncm, fthr_mode));

	  esl_stopwatch_Start(cfg->w_stage);
	  ///printf("filter  %3s  %-8s %5s %6s %5d %5g [", "-", DescribeFthrMode(fthr_mode), "-", "-", filN, esl_opt_GetReal(go, "--filE"));
	  printf("filter  %3s  %-8s %5s %6s %5d ", "-", DescribeFthrMode(fthr_mode), "-", "-", filN);
	  if     (esl_opt_GetBoolean(go, "--all")) printf("infty [");
	  else if(esl_opt_GetBoolean(go, "--ga"))  printf("   GA [");
	  else if(esl_opt_GetBoolean(go, "--nc"))  printf("   NC [");
	  else if(esl_opt_GetBoolean(go, "--tc"))  printf("   TC [");
	  else                                     printf("%5g [", esl_opt_GetReal(go, "--filE"));
	  fflush(stdout);

	  if(xstatus == eslOK) { if((status = update_cutoffs(go, cfg, errbuf, cm, gum_mode)) != eslOK) cm_Fail(errbuf); }
	  else { /* a worker has seen an error and we're trying to finish cleanly, but we still need to broadcast cutoffs, 
		  * their values don't matter */
	    for (p = 0; p < cfg->np; p++) cfg->cutoffA[p] = -eslINFINITY;
	  }
	  /* broadcast cutoffs */
	  ESL_DASSERT1((cfg->cutoffA != NULL));
	  MPI_Bcast(cfg->cutoffA, cfg->np, MPI_FLOAT, 0, MPI_COMM_WORLD);

	  if(xstatus == eslOK) have_work = TRUE;  /* TRUE while work remains  */
	  else                 have_work = FALSE; /* we've seen an error and are trying to finish cleanly */
	    
	  wi = 1;

	  nseq_sent = 0;
	  nseq_recv = 0;
	  update_i = filN / 20.;
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
#ifdef HAVE_DEVOPTS
			if (esl_opt_GetBoolean(go, "--hybrid")) { /* we receive worker_vscAA along with worker_vit_cp9scA and worker_fwd_cp9scA */
			  if ((status = cmcalibrate_cp9_filter_results_hyb_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, cm->M, &worker_vscAA, &worker_vit_cp9scA, &worker_fwd_cp9scA, &worker_partA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
			}
			else { /* only receive worker_vit_cp9scA, worker_fwd_cp9scA, no worker_vscAA, it's only used if --hybrid */
#endif
			  if ((status = cmcalibrate_cp9_filter_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_vit_cp9scA, &worker_fwd_cp9scA, &worker_partA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
#ifdef HAVE_DEVOPTS
			}
#endif
			ESL_DPRINTF1(("MPI master has unpacked HMM filter results\n"));
			ESL_DASSERT1((nseq_just_recv > 0));
			ESL_DASSERT1(((nseq_recv + nseq_just_recv) <= filN));
			for(i = 0; i < nseq_just_recv; i++) {
			  fil_vit_cp9scA[nseq_recv+i] = worker_vit_cp9scA[i];
			  fil_fwd_cp9scA[nseq_recv+i] = worker_fwd_cp9scA[i];
			  fil_partA[nseq_recv+i]      = worker_partA[i];
			  ESL_DASSERT1((fil_partA[nseq_recv+i] < cfg->np));
			  if(nseq_recv+i > update_i) {
			    printf("=");
			    fflush(stdout); 
			    update_i += filN / 20.; 
			  }
			}
#ifdef HAVE_DEVOPTS
			if(esl_opt_GetBoolean(go, "--hybrid")) {
			  for(v = 0; v < cm->M; v++) {
			    for(i = 0; i < nseq_just_recv; i++) {
			      ESL_DPRINTF3(("\tscore from worker v: %d i: %d sc: %f\n", i, v, worker_vscAA[v][i]));
			      printf("\tscore from worker v: %d i: %d sc: %f\n", i, v, worker_vscAA[v][i]);
			      fil_vscAA[v][nseq_recv+i] = worker_vscAA[v][i];
			    }
			    free(worker_vscAA[v]);
			  }
			  free(worker_vscAA);
			}
#endif
			free(worker_vit_cp9scA);
			free(worker_fwd_cp9scA);
			free(worker_partA);
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
		  /* send new search job */
		  ESL_DPRINTF1(("MPI master is sending HMM filter nseq %d to worker %d\n", nseq_this_worker, wi));
		  MPI_Send(&(nseq_this_worker), 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	      
		  wi++;
		  nproc_working++;
		  nseq_sent += nseq_this_worker;
		}
	    }

	  if(xstatus == eslOK) { 
	    /* predict speedup for HMM-only filter */
	    if((status = predict_cp9_filter_speedup(go, cfg, errbuf, cm, fil_vit_cp9scA, fil_fwd_cp9scA, fil_partA, cfg->cmstatsA[cmi]->bfA[fthr_mode])) != eslOK) cm_Fail(errbuf);
	    if(cfg->be_verbose) DumpBestFilterInfo(cfg->cmstatsA[cmi]->bfA[fthr_mode]);
	  }
	  ESL_DPRINTF1(("MPI master: done with HMM filter calc for fthr mode %d for this CM.\n", fthr_mode));
	  
	  for (wi = 1; wi < cfg->nproc; wi++) { 
	    msg = MPI_FINISHED_CP9_FILTER;
	    MPI_Send(&msg, 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	  }

	  esl_stopwatch_Stop(cfg->w_stage);
	  format_time_string(time_buf, cfg->w_stage->elapsed, 0);
	  printf("=] %10s\n", time_buf);
	}
	ESL_DPRINTF1(("MPI master: done with gumbel mode %d for this CM.\n", gum_mode));
      }
      ESL_DPRINTF1(("MPI master: done with this CM.\n"));
      if(xstatus == eslOK) if(cfg->be_verbose) { debug_print_cmstats(cfg->cmstatsA[cmi], (! esl_opt_GetBoolean(go, "--gumonly"))); }
      
      if(! (esl_opt_GetBoolean(go, "--filonly"))) { 
	for(v = 0; v < cm->M; v++) free(gum_vscAA[v]);
	free(gum_vscAA);
	free(gum_cp9scA);
      }
      if(! (esl_opt_GetBoolean(go, "--gumonly"))) { 
	for(v = 0; v < cm->M; v++) free(fil_vscAA[v]);
	free(fil_vscAA);
	free(fil_vit_cp9scA);
	free(fil_fwd_cp9scA);
	free(fil_partA);
      }
      FreeCM(cm);
      printf("//\n");
      fflush(stdout);
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
  CM_t         *cm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
  int           pos;
  char          errbuf[cmERRBUFSIZE];
  MPI_Status  mpistatus;
  float  **gum_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each random seq */
  float   *gum_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each random seq */
  float  **fil_vscAA  = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each emitted seq */
  float   *fil_vit_cp9scA = NULL; /*                [0..nseq-1] best cp9 Viterbi score for each emitted seq */
  float   *fil_fwd_cp9scA = NULL; /*                [0..nseq-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA      = NULL; /*                [0..nseq-1] partition of CM emitted seq */

  long     seed;  /* seed for RNG */
  int      gum_mode;
  int      fthr_mode;
  int working_on_cm;        /* TRUE when gum_mode is for CM gumbel */
  int working_on_cp9;       /* TRUE when gum_mode is for CP9 gumbel */
  int nseq;
  int v, p;
  void *tmp;
  int  cmi;
  int            L;  /* length of sequences to search for gumbel fitting, L==cm->W*2 unless --gumL enabled, 
		      * in which case L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL") */
  int in_fil_section_flag = FALSE; /* set to TRUE while we're in the filter threshold calculation
				    * section, we need to know this when we goto ERROR, b/c we have
				    * to know how many more MPI_Recv() calls to make to match up
				    * with the Master's sends before we can shut down.
				    */

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
   * case 1: if we've used the --gcfromdb option, we read in a seq file to fill
   * cfg->gc_freq, and we need that info here for the worker, so we receive
   * it's broadcast from the master
   * 
   * case 2: if we are calculating stats for more than 1 
   * partition, (--pfile invoked), we need to receive that information 
   * via broadcast from master. Otherwise we need to setup the default partition info
   * (single partition, 0..100 GC content)
   */
  if(! (esl_opt_IsDefault(go, "--gcfromdb"))) { /* receive gc_freq info from master */
    ESL_DASSERT1((cfg->gc_freq == NULL));
    ESL_ALLOC(cfg->gc_freq,  sizeof(double) * GC_SEGMENTS);
    ESL_ALLOC(cfg->pgc_freq, sizeof(double) * GC_SEGMENTS);
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else cfg->gc_freq = NULL; /* default */
  if(! (esl_opt_IsDefault(go, "--pfile"))) { /* receive partition info from master */
    MPI_Bcast(&(cfg->np),     1, MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pstart == NULL));
    ESL_ALLOC(cfg->pstart, sizeof(int) * cfg->np);
    MPI_Bcast(cfg->pstart, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else { /* no --pfile, set up default partition info */  
    cfg->np     = 1;
    ESL_ALLOC(cfg->pstart, sizeof(int) * cfg->np);
    cfg->pstart[0] = 0;
  }
  if(! (esl_opt_IsDefault(go, "--db"))) { /* receive dbsize for CM evalue cutoff for filter thr calc from master */
    MPI_Bcast(&(cfg->dbsize), 1, MPI_LONG, 0, MPI_COMM_WORLD);
  } /* else cfg->dbsize set to default 1 Mb when cfg was created */

  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      cfg->ncm++;  
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->cmstatsA, tmp, sizeof(CMStats_t *) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm(go, cfg, errbuf, cm))      != eslOK) goto ERROR;
      if((status = initialize_cmstats(go, cfg, errbuf, cm)) != eslOK) goto ERROR;
      if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) goto ERROR;
      
      for(gum_mode = 0; gum_mode < GUM_NMODES; gum_mode++) {

#ifdef HAVE_DEVOPTS
	/* free and recalculate hybrid scan info, b/c when investigating hybrid filters we may have added sub CM roots */
	if(cfg->hsi != NULL) cm_FreeHybridScanInfo(cfg->hsi, cm);
	cfg->hsi = cm_CreateHybridScanInfo(cm, esl_opt_GetReal(go, "--filbeta"), cfg->full_vcalcs[0]);
#endif

	ESL_DPRINTF1(("worker: %d gum_mode: %d nparts: %d\n", cfg->my_rank, gum_mode, cfg->np));
	if(GumModeIsForCM(gum_mode)) { working_on_cm = TRUE;  working_on_cp9 = FALSE; }
	else                         { working_on_cm = FALSE; working_on_cp9 = TRUE;  }

	/* do we need to switch from glocal configuration to local? */
	if(gum_mode > 0 && (! GumModeIsLocal(gum_mode-1)) && GumModeIsLocal(gum_mode)) 
	  if((status = switch_global_to_local(cm, errbuf)) != eslOK) goto ERROR;
	/* update search opts for gumbel mode */
	GumModeToSearchOpts(cm, gum_mode);
	/* if --id, we free RNG, then create a new one and reseed it with the initial seed,
	 * so we use same seqs for Gumbel fittings of all CM modes and HMM modes 
	 * (but if --cmN <n> != --hmmN <n>) they won't be the same between CM and HMM modes. 
	 */
	if(esl_opt_GetBoolean(go, "--id")) {
	  seed = esl_randomness_GetSeed(cfg->r);
	  esl_randomness_Destroy(cfg->r);
	  cfg->r = esl_randomness_Create(seed);
	}
	
	/* gumbel fitting section */
	if(! (esl_opt_GetBoolean(go, "--filonly"))) {
	  if(esl_opt_IsDefault(go, "--gumL")) L = cm->W*2; 
	  else                                L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL")); /* minimum L we allow is 2 * cm->W, this is enforced silently (!) */
	  for (p = 0; p < cfg->np; p++) { /* for each partition */
	    
	    ESL_DPRINTF1(("worker %d gum_mode: %d partition: %d\n", cfg->my_rank, gum_mode, p));
	    if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	  
	    if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	    while(nseq != MPI_FINISHED_GUMBEL) {
	      ESL_DPRINTF1(("worker %d: has received nseq: %d\n", cfg->my_rank, nseq));
	    
	      if(working_on_cm) {
		if((status = process_gumbel_workunit (go, cfg, errbuf, cm, nseq, L, &gum_vscAA, NULL, NULL)) != eslOK) goto ERROR;
		ESL_DPRINTF1(("worker %d: has gathered CM gumbel results\n", cfg->my_rank));
		n = 0;
		if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
		  ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
		n += sz;
		if (cmcalibrate_cm_gumbel_results_MPIPackSize(gum_vscAA, nseq, cm->M, MPI_COMM_WORLD, &sz) != eslOK)
		  ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cm_gumbel_results_MPIPackSize() call failed"); 
		n += sz;  
		if (n > wn) {
		  void *tmp;
		  ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
		  wn = n;
		}
		ESL_DPRINTF1(("worker %d: has calculated the CM gumbel results will pack into %d bytes\n", cfg->my_rank, n));
		status = eslOK;
		pos = 0;
		if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
		  ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
		if ((status = cmcalibrate_cm_gumbel_results_MPIPack(gum_vscAA, nseq, cm->M, wbuf, wn, &pos, MPI_COMM_WORLD)) != eslOK) 
		  ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cm_gumbel_results_MPIPack() call failed.");
		for(v = 0; v < cm->M; v++) free(gum_vscAA[v]);
	      }
	      else { /* working on cp9 */
		if((status = process_gumbel_workunit (go, cfg, errbuf, cm, nseq, L, NULL, &gum_cp9scA, NULL)) != eslOK) goto ERROR;
		ESL_DPRINTF1(("worker %d: has gathered CP9 gumbel results\n", cfg->my_rank));
		n = 0;
		if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
		  ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
		n += sz;
		if ((status = cmcalibrate_cp9_gumbel_results_MPIPackSize(gum_cp9scA, nseq, MPI_COMM_WORLD, &sz)) != eslOK)
		  ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_gumbel_results_MPIPackSize() call failed"); 
		n += sz;  

		if (n > wn) {
		  void *tmp;
		  ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
		  wn = n;
		}
		ESL_DPRINTF1(("worker %d: has calculated the CP9 gumbel results will pack into %d bytes\n", cfg->my_rank, n));
		status = eslOK;
		pos = 0;
		if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
		  ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
		if (cmcalibrate_cp9_gumbel_results_MPIPack(gum_cp9scA, nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) 
		  ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_gumbel_results_MPIPack() call failed.");
		free(gum_cp9scA);
	      }	    

	      MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	      ESL_DPRINTF1(("worker %d: has sent gumbel results to master in message of %d bytes\n", cfg->my_rank, pos));

	      /* receive next number of sequences, if MPI_FINISHED_GUMBEL, we'll stop */
	      if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	    }
	    ESL_DPRINTF1(("worker %d gum_mode: %d finished partition: %d\n", cfg->my_rank, gum_mode, p));
	  }
	  ESL_DPRINTF1(("worker %d finished all partitions for gum_mode: %d\n", cfg->my_rank, gum_mode));
	} /* end of if ! --filonly */
	
	/* filter threshold section */
	if(GumModeIsForCM(gum_mode) && (! (esl_opt_GetBoolean(go, "--gumonly")))) {
	  in_fil_section_flag = TRUE;
	  fthr_mode = GumModeToFthrMode(gum_mode);
	  ESL_DASSERT1((fthr_mode != -1));
	  ESL_DPRINTF1(("worker %d fthr_mode: %d\n", cfg->my_rank, fthr_mode));

	  /* get cutoffs for each partition from master, cfg->np never changes */
	  ESL_DASSERT1((cfg->cutoffA != NULL));
	  MPI_Bcast(cfg->cutoffA, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);

	  if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	  while(nseq != MPI_FINISHED_CP9_FILTER) {
	    ESL_DPRINTF1(("worker %d: has received hmm filter nseq: %d\n", cfg->my_rank, nseq));
	    

#ifdef HAVE_DEVOPTS
	    if(esl_opt_GetBoolean(go, "--hybrid")) { /* we want fil_vscAA filled with CM scanning scores */
	      if((status = process_filter_workunit (go, cfg, errbuf, cm, nseq, &fil_vscAA, &fil_vit_cp9scA, &fil_fwd_cp9scA, NULL, &fil_partA)) != eslOK) goto ERROR;
	    }
	    else { /* we don't want fil_vscAA filled with CM scanning scores */
#endif
	      if((status = process_filter_workunit (go, cfg, errbuf, cm, nseq, NULL, &fil_vit_cp9scA, &fil_fwd_cp9scA, NULL, &fil_partA)) != eslOK) goto ERROR;
#ifdef HAVE_DEVOPTS
	    }
#endif
	    ESL_DPRINTF1(("worker %d: has gathered HMM filter results\n", cfg->my_rank));
	    n = 0;

	    if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	      ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	    n += sz;

#ifdef HAVE_DEVOPTS
	    if (esl_opt_GetBoolean(go, "--hybrid")) { /* we send back fil_vscAA along with fil_vit_cp9scA and fil_fwd_cp9scA */
	      if(cmcalibrate_cp9_filter_results_hyb_MPIPackSize(nseq, cm->M, MPI_COMM_WORLD, &sz) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPackSize() call failed"); 
	    }
	    else { /* only send back fil_vit_cp9scA, fil_fwd_cp9scA, fil_vscAA is irrelevant, it's only used if --hybrid */
#endif
	      if(cmcalibrate_cp9_filter_results_MPIPackSize(nseq, MPI_COMM_WORLD, &sz) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPackSize() call failed"); 
#ifdef HAVE_DEVOPTS
	      }
#endif
	    n += sz;  
	    if (n > wn) {
	      void *tmp;
	      ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	      wn = n;
	    }
	    ESL_DPRINTF1(("worker %d: has calculated the HMM filter results will pack into %d bytes\n", cfg->my_rank, n));
	    status = eslOK;
	    pos = 0;
	    int i; for(i = 0; i < nseq; i++) assert(fil_partA[i] < cfg->np);

	    if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	      ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
#ifdef HAVE_DEVOPTS
	    if (esl_opt_GetBoolean(go, "--hybrid")) { /* we send back fil_vscAA along with fil_vit_cp9scA and fil_fwd_cp9scA */
	      if (cmcalibrate_cp9_filter_results_hyb_MPIPack(fil_vscAA, fil_vit_cp9scA, fil_fwd_cp9scA, fil_partA, nseq, cm->M, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPack() call failed"); 
	      for(v = 0; v < cm->M; v++) free(fil_vscAA[v]);
	      free(fil_vscAA);
	    }
	    else { /* only send back fil_vit_cp9scA, fil_fwd_cp9scA, fil_vscAA is irrelevant, it's only used if --hybrid */
#endif
	      if (cmcalibrate_cp9_filter_results_MPIPack(fil_vit_cp9scA, fil_fwd_cp9scA, fil_partA, nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPack() call failed"); 
#ifdef HAVE_DEVOPTS
	    }

#endif
	    free(fil_vit_cp9scA);
	    free(fil_fwd_cp9scA);
	    free(fil_partA);

	    MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	    ESL_DPRINTF1(("worker %d: has sent CP9 filter results to master in message of %d bytes\n", cfg->my_rank, pos));
	    /* receive next number of sequences, if MPI_FINISHED_GUMBEL, we'll stop */
	    if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	  }
	  in_fil_section_flag = FALSE;
	}
      } /* end of for(gum_mode = 0; gum_mode < GUM_NMODES; gum_mode++) */

      FreeCM(cm);
      cm = NULL;
      ESL_DPRINTF1(("worker %d finished all gum_modes for this cm.\n", cfg->my_rank));
    }
  if (status == eslEOD) ESL_DPRINTF1(("Worker %d told CMs are done.\n", cfg->my_rank));
  else goto ERROR;
  
  if (wbuf != NULL) free(wbuf);
  return;

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
  for(; gum_mode < GUM_NMODES; gum_mode++) {
    MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
    if(GumModeIsForCM(gum_mode) && (! (esl_opt_GetBoolean(go, "--gumonly")))) {
      MPI_Bcast(cfg->cutoffA, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
    }
  }
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);

  return;
}
#endif /*HAVE_MPI*/



/* Function: process_gumbel_workunit()
 * Date:     EPN, Mon Dec 10 06:09:09 2007
 *
 * Purpose:  A gumbel work unit consists of a CM, and an int specifying a 
 *           number of sequences <nseq>. The job is to randomly generate <nseq> 
 *           sequences using the cm->null background distribution, and 
 *           search them with either (a) the CM, (b) the CM's CP9 HMM, or
 *           (c) a hybrid CM/CP9 CYK/Viterbi scanning algorithm, with hybrid
 *           scanning info in cfg->hsi.
 *
 *           Thus, this function can be run in 1 of 3 modes, determined by the
 *           status of the input variables:
 *         
 *           Mode 1. Gumbel calculation for CM. 
 *           <ret_vscAA> != NULL, <ret_cp9scA> == NULL, <ret_hybscA> == NULL.
 *           Search random sequences with only the CM, either CYK or Inside
 *           (as specified by cm->search_opts>. <ret_vscAA> is filled
 *           with the best CM score at each state for each sequence.
 *
 *           Mode 2. Gumbel calculation for the CP9. 
 *           <ret_vscAA> == NULL, <ret_cp9scA> != NULL, <ret_hybscA> == NULL.
 *           Search random sequences with only the CP9, either Viterbi or Forward
 *           (as specified by cm->search_opts). <ret_cp9scA> is filled
 *           with the best CP9 score for each sequence.
 *
 *           Mode 3. Gumbel calculation for hybrid scanner.
 *           <ret_vscAA> == NULL, <ret_cp9scA> == NULL, <ret_hybscA> != NULL.
 *           Search random sequences with only a hybrid CM/CP9 scanner, 
 *           using hybrid info in cfg->hsi. <ret_hybscA> is filled
 *           with the best hybrid score for each sequence.
 *
 * Args:     go           - getopts
 *           cfg          - cmcalibrate's configuration
 *           errbuf       - for writing out error messages
 *           cm           - the CM (already configured as we want it)
 *           nseq         - number of seqs to generate
 *           L            - length of sequences to search, L==cm->W*2 unless --gumL enabled, in which case
 *                          L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL")
 *           ret_vscAA    - RETURN: [0..v..cm->M-1][0..nseq-1] best score at each state v for each seq
 *           ret_cp9scA   - RETURN: [0..nseq-1] best CP9 score for each seq
 *           ret_hybscA   - RETURN: [0..nseq-1] best hybrid score for each seq
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_gumbel_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			int L, float ***ret_vscAA, float **ret_cp9scA, float **ret_hybscA)
{
  int            status;
  int            mode; /* 1, 2, or 3, determined by status of input args, as explained in 'Purpose' above. */
  float        **vscAA        = NULL;  /* [0..v..cm->M-1][0..i..nseq-1] best CM score for each state, each seq */
  float         *cur_vscA     = NULL;  /* [0..v..cm->M-1]               best CM score for each state cur seq */
  float         *cp9scA       = NULL;  /*                [0..i..nseq-1] best CP9 score for each seq, */
  float         *hybscA       = NULL;  /*                [0..i..nseq-1] best hybrid score for each seq */
  double        *dnull        = NULL; /* double version of cm->null, for generating random seqs */
  int            i;
  int            v;
  ESL_DSQ       *dsq;
  float          sc;
  float          update_i = nseq / 20.;

  /* determine mode, and enforce mode-specific contract */
  if     (ret_vscAA != NULL && ret_cp9scA == NULL && ret_hybscA == NULL) mode = 1; /* calcing CM     gumbel stats */
  else if(ret_vscAA == NULL && ret_cp9scA != NULL && ret_hybscA == NULL) mode = 2; /* calcing CP9    gumbel stats */
  else if(ret_vscAA == NULL && ret_cp9scA == NULL && ret_hybscA != NULL) mode = 3; /* calcing hybrid gumbel stats */
  else ESL_FAIL(eslEINCOMPAT, errbuf, "can't determine mode in process_gumbel_workunit.");

  ESL_DPRINTF1(("in process_gumbel_workunit nseq: %d L: %d mode: %d\n", nseq, L, mode));

  int do_cyk     = FALSE;
  int do_inside  = FALSE;
  int do_viterbi = FALSE;
  int do_forward = FALSE;
  int do_hybrid  = FALSE;
  /* determine algs we'll use and allocate the score arrays we'll pass back */
  if(mode == 1) {
    if(cm->search_opts & CM_SEARCH_INSIDE) do_inside = TRUE;
    else                                   do_cyk    = TRUE;
    ESL_ALLOC(vscAA, sizeof(float *) * cm->M);
    for(v = 0; v < cm->M; v++) ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    ESL_ALLOC(cur_vscA, sizeof(float) * cm->M);
  }
  else if(mode == 2) {
    if(cm->search_opts & CM_SEARCH_HMMVITERBI) do_viterbi = TRUE;
    if(cm->search_opts & CM_SEARCH_HMMFORWARD) do_forward = TRUE;
    if((do_viterbi + do_forward) > 1) ESL_FAIL(eslEINVAL, errbuf, "process_gumbel_workunit, mode 2, and cm->search_opts CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags both raised.");
    ESL_ALLOC(cp9scA, sizeof(float) * nseq); /* will hold Viterbi or Forward scores */
  }
#ifdef HAVE_DEVOPTS
  else if(mode == 3) {
    do_hybrid = TRUE;
    ESL_ALLOC(hybscA,       sizeof(float) * nseq); /* will hold hybrid scores */
  }
#endif
  else if(mode == 3) { /* never entered if HAVE_DEVOPTS is defined */
    ESL_FAIL(eslEINCOMPAT, errbuf, "process_gumbel_workunit(), mode 3 unavailable (HAVE_DEVOPTS is undefined)");
  }

  ESL_DPRINTF1(("do_cyk:     %d\ndo_inside:  %d\ndo_viterbi: %d\ndo_forward: %d\ndo_hybrid: %d", do_cyk, do_inside, do_viterbi, do_forward, do_hybrid)); 
  
  /* fill dnull, a double version of cm->null, but only if we're going to need it to generate random seqs */
  if(cfg->pgc_freq == NULL) {
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);    
  }
  
  /* generate dsqs one at a time and collect best CM scores at each state and/or best overall CP9 score */
  for(i = 0; i < nseq; i++) {
    if(cfg->my_rank == 0 && i > update_i) { /* print status update to stdout */
      printf("=");
      fflush(stdout); 
      update_i += nseq / 20.; 
    }
    if((status = get_random_dsq(cfg, errbuf, cm, dnull, L, &dsq)) != eslOK) return status; 

    /* if nec, search with CM */
    if (do_cyk)    if((status = FastCYKScan    (cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc)) != eslOK) return status;
    if (do_inside) if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc)) != eslOK) return status;
    /* if nec, search with CP9 */
    if (do_viterbi) 
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       NULL,   /* don't want traces back */
			       &(cp9scA[i]))) != eslOK) return status;
    if (do_forward) {
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       &(cp9scA[i]))) != eslOK) return status;
    }
#ifdef HAVE_DEVOPTS
    if (do_hybrid) {
      if((status = cm_cp9_HybridScan(cm, errbuf, cm->cp9_mx, dsq, cfg->hsi, 1, L, cfg->hsi->W, 0., 
				     NULL, /* don't report results */
				     NULL, /* don't want best score at each posn back */
				     NULL, /* don't want the max scoring posn back */
				     &(hybscA[i]))) != eslOK) return status;
    }
#endif

    /*to print seqs to stdout uncomment this block 
    ESL_SQ *tmp;
    tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
    esl_sq_Textize(tmp);
    printf(">seq%d\n%s\n", i, tmp->seq);
    esl_sq_Destroy(tmp);
    */

    free(dsq);
    if (cur_vscA != NULL) /* will be NULL if do_cyk == do_inside == FALSE (mode 2) */
      for(v = 0; v < cm->M; v++) vscAA[v][i] = cur_vscA[v];
    free(cur_vscA);
  }
  if(cfg->my_rank == 0) { printf("=]"); }

  if(dnull != NULL) free(dnull);
  if(ret_vscAA  != NULL) *ret_vscAA  = vscAA;
  if(ret_cp9scA != NULL) *ret_cp9scA = cp9scA;
  if(ret_hybscA != NULL) *ret_hybscA = hybscA;
  return eslOK;

 ERROR:
  return status;
  }


/* Function: process_filter_workunit()
 * Date:     EPN, Mon Dec 10 05:48:35 2007
 *
 * Purpose:  A filter work unit consists of a CM, an int specifying a 
 *           number of sequences <nseq>, and a flag indicating how to search
 *           the sequences. The job is to generate <nseq> sequences from the
 *           CM and search them, the way they're searched is mode dependent
 *           (see below).  with either (a) the CM using bands from
 *           hybrid scanning info in cfg->hsi, then the CP9 HMM with Viterbi and 
 *           Forward or (b) using the hybrid CM/CP9 CYK/Viterbi algorithm
 *           with the hybrid scanning info in cfg->hsi.
 *
 *           Thus, this function can be run in 1 of 3 modes, determined by the
 *           status of the input variables. Note modes 2 and 3 are only possible
 *           if the --hybrid option is enabled, which is only even available if
 *           HAVE_DEVOPTS is defined.
 *         
 *           Mode 1. Scores will be used for calc'ing filter threshold of CP9 HMM.
 *           <ret_vscAA> == NULL, <ret_vit_cp9scA> != NULL, <ret_fwd_cp9scA> != NULL, <ret_hyb_cmscA> == NULL
 *           Emit from CM and score with CP9 Viterbi and Forward, <ret_vit_cp9scA> 
 *           and <ret_fwd_cp9scA> are filled with the best CP9 Viterbi/HMM score 
 *           for each sequence.
 *
 *           Mode 2. Scores will be used for calc'ing filter threshold of CP9 HMM
 *           and CM scores will be used to predict which sub CM roots will be good 
 *           at filtering.
 *           <ret_vscAA> != NULL, <ret_vit_cp9scA> != NULL, <ret_fwd_cp9scA> != NULL, <ret_hyb_cmscA> == NULL
 *           Emit from CM and search first with CM using QDBs from hybrid scanning
 *           info in cfg->hsi, best score from each state of the CM for each 
 *           seq is stored in >ret_vscAA>. Then search (same seq) with CM CP9 Viterbi and Forward, 
 *           <ret_vit_cp9scA> and <ret_fwd_cp9scA> are filled with the best CP9 
 *           Viterbi/HMM score for each sequence.
 *
 *           Mode 3. Scores will be used for calc'ing filter threshold of hybrid scanner.
 *           <ret_vscAA> == NULL, <ret_vit_cp9scA> == NULL, <ret_fwd_cp9scA> == NULL, <ret_hybscA> != NULL
 *           Emit from CM and score with hybrid CM/CP9 CYK/Viterbi scanner, 
 *           <ret_hybscA> are filled with the best hybrid scanner scores 
 *           for each sequence.
 *
 * Args:     go             - getopts
 *           cfg            - cmcalibrate's configuration
 *           errbuf         - for writing out error messages
 *           cm             - the CM (already configured as we want it)
 *           nseq           - number of seqs to generate
 *           ret_vscAA      - RETURN: [0..v..cm->M-1][0..nseq-1] best score at each state v for each seq
 *           ret_vit_cp9scA - RETURN: [0..nseq-1] best Viterbi CP9 score for each seq
 *           ret_fwd_cp9scA - RETURN: [0..nseq-1] best Forward CP9 score for each seq
 *           ret_hybscA     - RETURN: [0..nseq-1] best Hybrid CM/CP9 score for each seq
 *           ret_partA      - RETURN: [0..nseq-1] partition of each seq 
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			float ***ret_vscAA, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, float **ret_hybscA, int **ret_partA)
{
  int            status;
  int            mode; /* 1 or 2 determined by status of input args, as explained in 'Purpose' above. */
  float        **vscAA      = NULL;  /* [0..v..cm->M-1][0..i..nseq-1] best CM score for each state, each seq */
  float         *cur_vscA   = NULL;  /* [0..v..cm->M-1]               best CM score for each state cur seq */
  float         *vit_cp9scA = NULL;  /* [0..i..nseq-1] best CP9 Viterbi score for each seq */
  float         *fwd_cp9scA = NULL;  /* [0..i..nseq-1] best CP9 Forward score for each seq */
  float         *hybscA     = NULL;  /* [0..i..nseq-1] best hybrid CM/CP9 scanner score for each seq */
  int           *partA      = NULL;  /* [0..i..nseq-1] partitions of each seq */
  int            p;                  /* what partition we're in, not used unless emit_from_cm = TRUE */
  int            i, v;
  int            L;
  int            nfailed = 0;
  Parsetree_t   *tr;
  ESL_DSQ       *dsq;
  float          sc;
  int            inside_flag_raised = FALSE;
  float          update_i = nseq / 20.;


  /* determine mode, and enforce mode-specific contract */
  if     (ret_vscAA == NULL && ret_vit_cp9scA != NULL && ret_fwd_cp9scA != NULL && ret_hybscA == NULL) mode = 1; /* running CP9 Viterbi and Forward */
#if HAVE_DEVOPTS
  else if(ret_vscAA != NULL && ret_vit_cp9scA != NULL && ret_fwd_cp9scA != NULL && ret_hybscA == NULL) mode = 2; /* running CM CYK and CP9 Viterbi and Forward */
  else if(ret_vscAA == NULL && ret_vit_cp9scA == NULL && ret_fwd_cp9scA == NULL && ret_hybscA != NULL) mode = 3; /* running hybrid CM/CP9 scanner */
#endif
  else ESL_FAIL(eslEINCOMPAT, errbuf, "can't determine mode in process_filter_workunit.");
  if (ret_partA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "ret_partA is NULL in process_filter_workunit.");

  ESL_DPRINTF1(("in process_filter_workunit nseq: %d mode: %d\n", nseq, mode));
  /* if we get this far, if HAVE_DEVOPTS is undefined, mode MUST be 1 */

#ifndef HAVE_DEVOPTS
  if(mode != 1) ESL_FAIL(eslEINCOMPAT, errbuf, "HAVE_DEVOPTS is undefined, but mode is not 1 in process_filter_workunit(). This shoudln't happen.");
#endif

  /* determine algs we'll use and allocate the score arrays we'll pass back */
  ESL_ALLOC(partA, sizeof(int) * nseq); /* will hold partitions */

  if(mode == 1 || mode == 2) {
    ESL_ALLOC(vit_cp9scA, sizeof(float) * nseq); /* will hold Viterbi scores */
    ESL_ALLOC(fwd_cp9scA, sizeof(float) * nseq); /* will hold Forward scores */
    if(mode == 2) { 
      ESL_ALLOC(vscAA, sizeof(float *) * cm->M);
      for(v = 0; v < cm->M; v++) ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
      ESL_ALLOC(cur_vscA, sizeof(float) * cm->M);
    }
    if(cm->search_opts & CM_SEARCH_INSIDE) { inside_flag_raised = TRUE; cm->search_opts &= ~CM_SEARCH_INSIDE; }
    else inside_flag_raised = FALSE;
  }
  else  /* mode == 3 */
    ESL_ALLOC(hybscA, sizeof(float) * nseq); /* will hold hybrid scores */

  /* generate dsqs one at a time and collect best CM scores at each state and/or best overall CP9 score */
  for(i = 0; i < nseq; i++) {
    if(cfg->my_rank == 0 && i > update_i) { /* print status update to stdout */
      printf("=");
      fflush(stdout); 
      update_i += nseq / 20.; 
    }
    if((status = get_cmemit_dsq(cfg, errbuf, cm, &L, &p, &tr, &dsq)) != eslOK) return status;
    /* we only want to use emitted seqs with a sc > cutoff */
    if((status = cm_find_hit_above_cutoff(go, cfg, errbuf, cm, dsq, tr, L, cfg->cutoffA[p], &sc)) != eslOK) return status;
    while(sc < cfg->cutoffA[p]) { 
      free(dsq); 	
      /* parsetree tr is freed in cm_find_hit_above_cutoff() */
      if((status = get_cmemit_dsq(cfg, errbuf, cm, &L, &p, &tr, &dsq)) != eslOK) return status;
      nfailed++;
      if(nfailed > 1000 * nseq) ESL_FAIL(eslERANGE, errbuf, "process_filter_workunit(), max number of failures (%d) reached while trying to emit %d seqs.\n", nfailed, nseq);
      if((status = cm_find_hit_above_cutoff(go, cfg, errbuf, cm, dsq, tr, L, cfg->cutoffA[p], &sc)) != eslOK) return status;
    }

    /*to print seqs to stdout uncomment this block  
    ESL_SQ *tmp;
    tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
    esl_sq_Textize(tmp);
    printf(">seq%d\n%s\n", i, tmp->seq);
    esl_sq_Destroy(tmp);
    */

    partA[i] = p;
    assert(partA[i] < cfg->np);
    ESL_DPRINTF1(("i: %d nfailed: %d cutoff: %.3f p: %d\n", i, nfailed, cfg->cutoffA[p], p));

    /* search dsq with mode-specific search algs */
    if(mode == 1 || mode == 2) {
      /* Note: in mode 2, with FastCYKScan, we use cfg->hsi->smx scan matrix, which may have qdbs calc'ed differently than cm->smx */
      if(mode == 2) { if((status = FastCYKScan(cm, errbuf, cfg->hsi->smx, dsq, 1, L, 0., NULL, &(cur_vscA), NULL)) != eslOK) return status; }
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       NULL,   /* don't want traces back */
			       &(vit_cp9scA[i]))) != eslOK) return status;
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       &(fwd_cp9scA[i]))) != eslOK) return status;
    }
    else { /* mode == 3 */
      if((status = cm_cp9_HybridScan(cm, errbuf, cm->cp9_mx, dsq, cfg->hsi, 1, L, cfg->hsi->W, 0., 
				     NULL, /* don't report results */
				     NULL, /* don't want best score at each posn back */
				     NULL, /* don't want the max scoring posn back */
				     &(hybscA[i]))) != eslOK) return status;
    }
    free(dsq);
    if (cur_vscA != NULL) /* will be NULL if do_cyk == do_inside == FALSE (mode 3) */
      for(v = 0; v < cm->M; v++) vscAA[v][i] = cur_vscA[v];
    free(cur_vscA);
  }
  if(cfg->my_rank == 0) { printf("=]"); }
  *ret_partA = partA;
  if(ret_vscAA      != NULL)  *ret_vscAA      = vscAA;
  if(ret_vit_cp9scA != NULL)  *ret_vit_cp9scA = vit_cp9scA;
  if(ret_fwd_cp9scA != NULL)  *ret_fwd_cp9scA = fwd_cp9scA;
  if(ret_hybscA != NULL)      *ret_hybscA     = hybscA;

  if(inside_flag_raised) cm->search_opts |= CM_SEARCH_INSIDE;

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
  int nstarts, nexits, nd;

  cm->beta   = esl_opt_GetReal(go, "--beta"); /* this will be 1e-7 (default beta) unless changed at command line */

  /* config QDB? Yes, unless --noqdb enabled */
  if(! (esl_opt_GetBoolean(go, "--noqdb"))) 
    cm->config_opts |= CM_CONFIG_QDB;   /* configure QDB */
  else
    cm->search_opts |= CM_SEARCH_NOQDB; /* don't use QDB to search */

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
  /* process the --gemit option, this option forces all emitted parsetrees to be 'global'
   * in that they'll never contain a local begin or local end. */
  if(esl_opt_GetBoolean(go, "--gemit")) { 
    cm->flags |= CM_EMIT_NO_LOCAL_BEGINS; 
    cm->flags |= CM_EMIT_NO_LOCAL_ENDS;
  }

  ConfigCM(cm, NULL, NULL);
  
  if(esl_opt_GetBoolean(go, "--noqdb")) { /* setup cm->W */ 
    if(cm->dmin != NULL || cm->dmax != NULL) 
      cm_Fail("initialize_cm() --noqdb enabled, but cm->dmin and cm->dmax non-null. This shouldn't happen.");
    int *dmin;
    int *dmax;
    int safe_windowlen = cm->clen * 2;
    while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta, 0, &(dmin), &(dmax), NULL, NULL)))
      {
	free(dmin);
	free(dmax);
	safe_windowlen *= 2;
	if(safe_windowlen > (cm->clen * 1000))
	  cm_Fail("initialize_cm(), safe_windowlen big: %d\n", safe_windowlen);
      }
    cm->W = dmax[0];
    free(dmin);
    free(dmax);
    CMLogoddsify(cm); /* QDB calculation invalidates log odds scores */
  }

  /* count number of DP calcs */
  if(cfg->full_vcalcs != NULL) free(cfg->full_vcalcs);
  if((status = cm_CountSearchDPCalcs(cm, errbuf, 1000, cm->dmin, cm->dmax, cm->W, &(cfg->full_vcalcs), NULL)) != eslOK) return status;

  /* create and initialize scan info for CYK/Inside scanning functions */
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("initialize_cm(), CreateScanMatrixForCM() call failed.");
  
  return eslOK;
}


/* initialize_cmstats()
 * Allocate and initialize a cmstats object in the cfg->cmstatsA array. 
 */
static int
initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;
  int i;
  int p;
  int cmi = cfg->ncm-1;

  ESL_DPRINTF1(("initializing cmstats for %d partitions\n", cfg->np));

  cfg->cmstatsA[cmi] = AllocCMStats(cfg->np);
  
  if(esl_opt_GetBoolean(go, "--filonly")) { 
    /* set the cfg->np if this is the first CM */
    if(cfg->ncm == 1) { 
      cfg->np = cm->stats->np;
      if(cfg->cutoffA != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "--filonly invoked, and we're on the first CM, but cfg->cutoffA already allocated.");
      ESL_ALLOC(cfg->cutoffA, sizeof(float) * cfg->np);
    }
    /* Make sure the rare, rare case that would be a real pain in the ass to implement isn't happening:
     * the case when we use --filonly with a CM file with > 1 CMs, and at least 2 of those CMs have
     * gumbel stats for different numbers of partitions. If we did deal with this case then 
     * it'd be a bitch to deal with, b/c cfg->np could change (and we'd have to send that
     * info to the workers for each CM in MPI mode).
     */
    ESL_DASSERT1((esl_opt_IsDefault(go, "--gcfromdb"))); /* getopts should enforce this */
    ESL_DASSERT1((esl_opt_IsDefault(go, "--pfile")));  /* getopts should enforce this */
    if(! (cm->flags & CMH_GUMBEL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "--filonly invoked by CM has no gumbel stats in initialize_cmstats()\n");
    if(cfg->np != cm->stats->np)         ESL_FAIL(eslEINCOMPAT, errbuf, "--filonly invoked and CM file has CMs with different numbers of partitions. We can't deal. Either split CMs into different files, or recalibrate them fully (without --filonly)");    
    /* with --filonly, we're only calc'ing filter thresholds, so we copy the Gumbel stats from cm->stats. */
    esl_vec_ICopy(cm->stats->ps,   cfg->np, cfg->cmstatsA[cmi]->ps);
    esl_vec_ICopy(cm->stats->pe,   cfg->np, cfg->cmstatsA[cmi]->pe);
    esl_vec_ICopy(cm->stats->gc2p, GC_SEGMENTS, cfg->cmstatsA[cmi]->gc2p);
    for(i = 0; i < GUM_NMODES; i++) { 
      for(p = 0; p < cfg->np; p++) {
	cfg->cmstatsA[cmi]->gumAA[i][p]->N      = cm->stats->gumAA[i][p]->N;
	cfg->cmstatsA[cmi]->gumAA[i][p]->L      = cm->stats->gumAA[i][p]->L;
	cfg->cmstatsA[cmi]->gumAA[i][p]->mu     = cm->stats->gumAA[i][p]->mu;
	cfg->cmstatsA[cmi]->gumAA[i][p]->lambda = cm->stats->gumAA[i][p]->lambda;
	cfg->cmstatsA[cmi]->gumAA[i][p]->is_valid = cm->stats->gumAA[i][p]->is_valid;
      }
    }

    return eslOK;
  }

  /* if we get here --filonly was not invoked */
  ESL_DASSERT1((cfg->pstart[0] == 0));
  for(p = 0; p < cfg->np;     p++) cfg->cmstatsA[cmi]->ps[p] = cfg->pstart[p];
  for(p = 0; p < (cfg->np-1); p++) cfg->cmstatsA[cmi]->pe[p] = cfg->pstart[p+1]-1;
  cfg->cmstatsA[cmi]->pe[(cfg->np-1)] = GC_SEGMENTS-1; /* this is 100 */
  
  for(p = 0; p < cfg->np; p++)
    for(i = cfg->cmstatsA[cmi]->ps[p]; i <= cfg->cmstatsA[cmi]->pe[p]; i++)
      cfg->cmstatsA[cmi]->gc2p[i] = p; 
  
  /* master only allocations, workers don't need this */
  if(cfg->my_rank == 0) { 
    /* if they're NULL, allocate cfg->vmuAA, cfg->vlambdaAA, and cfg->gum_hybA
     * otherwise they're fine as they are because number of partitions never changes.
     */
    if(cfg->vmuAA == NULL) {
      ESL_ALLOC(cfg->vmuAA,     sizeof(double *) * cfg->np);
      ESL_ALLOC(cfg->vlambdaAA, sizeof(double *) * cfg->np);
      for(p = 0; p < cfg->np; p++) {
	cfg->vmuAA[p]     = NULL;
	cfg->vlambdaAA[p] = NULL;
      }
    }
    if(cfg->gum_hybA == NULL) {
      ESL_ALLOC(cfg->gum_hybA, sizeof(GumbelInfo_t *) * cfg->np);
      for(p = 0; p < cfg->np; p++) { 
	ESL_ALLOC(cfg->gum_hybA[p], sizeof(GumbelInfo_t));
	cfg->gum_hybA[p]->is_valid = FALSE;
      }
    }
  }
  if(cfg->cutoffA == NULL) ESL_ALLOC(cfg->cutoffA, sizeof(float) * cfg->np);

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
  double         e_cutoff;       /* E-value cutoff for each partition */
  int            sc_cutoff;      /* bit score cutoff for each partition */
  int            p;              /* partition index */
  double         mu;             /* mu for a requested db size (which is 1Mb unless --db enabled) */
  int            revert_to_default_filE; /* if --ga, --nc, or --tc enabled but CM does not have GA, NC or TC cutoff in the CM file,
					  * and we've already calibrated at least 1 CM in this CM file, pretend like --ga, --nc, or --tc,
					  * was not enabled by reverting to the default --filE value of 0.1
					  */
  /* if --ga, --nc, or --tc: 
   * cfg->filE:                set as E-value (in 1Mb db) that the GA, NC, or TC bit score corresponds to for each partition in this fthr_mode
   * cfg->cutoffA[0..p..np-1]: set as GA, NC, or TC bit score for all p for this fthr_mode
   */
  if ((esl_opt_GetBoolean(go, "--ga")) || (esl_opt_GetBoolean(go, "--nc")) || (esl_opt_GetBoolean(go, "--tc"))) {
    if(esl_opt_GetBoolean(go, "--ga")) { 
      if(! (cm->flags & CMH_GA)) {
	if(cfg->ncm > 1) revert_to_default_filE = TRUE;
	else             cm_Fail("--ga enabled but first CM in CM file does not have a Rfam GA cutoff.");
      }
      else sc_cutoff = cm->ga;
    }
    if(esl_opt_GetBoolean(go, "--nc")) { 
      if(! (cm->flags & CMH_NC)) {
	if(cfg->ncm > 1) revert_to_default_filE = TRUE;
	else             cm_Fail("--nc enabled but first CM in CM file does not have a Rfam NC cutoff.");
      }
      else sc_cutoff = cm->nc;
    }
    if(esl_opt_GetBoolean(go, "--tc")) { 
      if(! (cm->flags & CMH_TC)) {
	if(cfg->ncm > 1) revert_to_default_filE = TRUE;
	else             cm_Fail("--tc enabled but first CM in CM file does not have a Rfam TC cutoff.");
      }
      else sc_cutoff = cm->tc;
    }
    if(! revert_to_default_filE) { /* we've set sc_cutoff above, now determine e_cutoff for each partition */
      for (p = 0; p < cfg->np; p++) 
	cfg->cutoffA[p] = sc_cutoff; /* either cm->ga, cm->nc, or cm->tc as set above */
      return eslOK; /* we're done */
    }
  }

  if(esl_opt_GetBoolean(go, "--all")) {
    for (p = 0; p < cfg->np; p++)
      cfg->cutoffA[p] = -eslINFINITY;
  }
  else { /* either none of: --filE, --ga, --nc, --tc, --all were enabled, or --filE was enabled, 
	  * or --ga, --nc, --tc were enabled, but CM does not have cm->ga, cm->nc, or cm->tc and CM is not first in file */
    e_cutoff = esl_opt_GetReal(go, "--filE"); 
    for (p = 0; p < cfg->np; p++) {
      /* first determine mu based on db_size */
      tmp_K = exp(cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->mu * cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda) / 
	cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->L;
      mu = log(tmp_K  * ((double) cfg->dbsize)) / cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cfg->cutoffA[p] = mu - (log(e_cutoff) / cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda);
    }
  }
  return eslOK;
}  

/* Function: set_partition_gc_freq()
 * Date:     EPN, Mon Sep 10 08:00:27 2007
 *
 * Purpose:  Set up the GC freq to sample from for the current partition. 
 *           Only used if --gcfromdb used to read in dbseq from which to derive
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
  begin = cfg->pstart[p];
  if(p == (cfg->np-1)) end = (GC_SEGMENTS-1); /* this is 100 */
  else end = cfg->pstart[p+1] - 1;
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
  int status;
  double mu;
  double lambda;
  int i;
  double *xv;         /* raw data from histogram */
  int     n,z;  
  float tailfit;

  ESL_HISTOGRAM *h = NULL;       /* histogram of scores */


  /* Initialize histogram; these numbers are guesses */
  if((h = esl_histogram_CreateFull(-100., 100., .25)) == NULL) return eslEMEM;    

  /* fill histogram */
  for(i = 0; i < nscores; i++)
    if((status = esl_histogram_Add(h, scores[i])) != eslOK) ESL_FAIL(status, errbuf, "fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);

  /* fit scores to a gumbel */
  tailfit = esl_opt_GetReal(go, "--gtail");
  esl_histogram_GetTailByMass(h, tailfit, &xv, &n, &z); /* fit to right 'tailfit' fraction, 0.5 by default */
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
  
  if(esl_opt_GetBoolean(go, "--hybrid")) { /* fit gumbels for each candidate sub CM root for a hybrid filter */
    ESL_ALLOC(cfg->vmuAA[p],     sizeof(double) * cm->M);
    ESL_ALLOC(cfg->vlambdaAA[p], sizeof(double) * cm->M);
    
    for(v = 0; v < cm->M; v++) {
      if(cfg->hsi->iscandA[v]) {
	/* printf("FITTING v: %d sttype: %d\n", v, cm->sttype[v]); */
	if((status = fit_histogram(go, cfg, errbuf, vscA[v], nscores, &(cfg->vmuAA[p][v]), &(cfg->vlambdaAA[p][v]))) != eslOK) return status;
      }
      else cfg->vmuAA[p][v] = cfg->vlambdaAA[p][v] = -1.;
    }
  } 
  else { /* only fit root state 0, the full model */
    ESL_ALLOC(cfg->vmuAA[p],     sizeof(double) * 1);
    ESL_ALLOC(cfg->vlambdaAA[p], sizeof(double) * 1);
    if((status = fit_histogram(go, cfg, errbuf, vscA[0], nscores, &(cfg->vmuAA[p][0]), &(cfg->vlambdaAA[p][0]))) != eslOK) return status;
  }
  return eslOK;

 ERROR:
  return status;
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
 * Returns:  eslOK on success, ret_dsq filled with newly alloc'ed ESL_DSQ *,
 *           some other status code on failure.
 */
int
get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, double *dnull, int L, ESL_DSQ **ret_dsq)
{
  int status;
  double  gc_comp;
  double *distro = NULL;
  int do_free_distro = FALSE;
  ESL_DSQ *dsq = NULL;

  /* contract check, make sure we're in a valid mode */
  if(cfg->pgc_freq == NULL && dnull == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->pgc_freq == NULL and dnull == NULL");
  if(cfg->pgc_freq != NULL && dnull != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->pgc_freq != NULL and dnull != NULL");

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
  if ((status = esl_rnd_xIID(cfg->r, distro, cm->abc->K, L, dsq) != eslOK)) return status;

  if (do_free_distro) free(distro);
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  return status;
}

/* Function: get_cmemit_dsq()
 * Date:     EPN, Tue Sep 11 08:51:33 2007
 * 
 * Purpose:  Generate a dsq from a CM and return it.
 *
 * Returns:  eslOK on success, ESL_DSQ is filled with newly alloc'ed dsq; some other status code on an error, 
 */
int
get_cmemit_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int *ret_L, int *ret_p, Parsetree_t **ret_tr, ESL_DSQ **ret_dsq)
{
  int status;
  int p;
  int L;
  ESL_SQ *sq;
  ESL_DSQ *dsq;
  Parsetree_t *tr;

  if((status = EmitParsetree(cm, errbuf, cfg->r, "irrelevant", TRUE, &tr, &sq, &L)) != eslOK) return status;
  while(L == 0) { 
    FreeParsetree(tr); 
    esl_sq_Destroy(sq); 
    if((status = EmitParsetree(cm, errbuf, cfg->r, "irrelevant", TRUE, &tr, &sq, &L)) != eslOK) return status;
  }

  /* determine the partition */
  p = cfg->cmstatsA[cfg->ncm-1]->gc2p[(get_gc_comp(sq, 1, L))]; /* this is slightly wrong, 1,L for get_gc_comp() should be i and j of best hit */
  assert(p < cfg->np);
  ESL_DASSERT1((p < cfg->np));

  /* free everything allocated by a esl_sqio.c:esl_sq_CreateFrom() call, but the dsq */
  dsq = sq->dsq;
  free(sq->name);
  free(sq->acc);
  free(sq->desc);
  free(sq);

  *ret_L  = L;
  *ret_p  = p;
  *ret_tr = tr;
  *ret_dsq = dsq;
  return eslOK;
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
 *           ret_sc          - score of a hit within dsq, if < cutoff,
 *                             this is score of best hit within dsq, which
 *                             means no hit with sc > cutoff exists. If > cutoff,
 *                             not necessarily score of best hit within dsq.
 *
 * Returns:  eslOK on success. other status code upon failure, errbuf filled with error message.
 */
int 
cm_find_hit_above_cutoff(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq,
			 Parsetree_t *tr, int L, float cutoff, float *ret_sc)
{
  int status;
  int turn_qdb_back_on = FALSE;
  int turn_hbanded_back_off = FALSE;
  int turn_hmmscanbands_back_off = FALSE;
  double orig_tau = cm->tau;
  float sc;
  float size_limit = esl_opt_GetReal(go, "--mxsize");

#if eslDEBUGLEVEL >= 1
  int init_flags       = cm->flags;
  int init_search_opts = cm->search_opts;
#endif

  if(ret_sc == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_find_hit_above_cutoff(), ret_sc == NULL.\n");

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
  if(sc > cutoff || L == 0) { /* parse score exceeds cutoff, or zero length sequence (only 1 path is possible, must be parse score) */
    ESL_DASSERT1((cm->flags       == init_flags));
    ESL_DASSERT1((cm->search_opts == init_search_opts));
    /* printf("0 sc: %10.4f\n", sc); */
    *ret_sc = sc;
    return eslOK;
  } 

  if(!(cm->search_opts & CM_SEARCH_NOQDB))        turn_qdb_back_on = TRUE;
  if(!(cm->search_opts & CM_SEARCH_HBANDED))      turn_hbanded_back_off = TRUE;
  if(!(cm->search_opts & CM_SEARCH_HMMSCANBANDS)) turn_hmmscanbands_back_off = TRUE;

  cm->search_opts |= CM_SEARCH_NOQDB;

  /* stage 1 */
  cm->search_opts |= CM_SEARCH_HBANDED;
  cm->tau = 0.01;
  if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, TRUE, 0)) != eslOK) return status;
  status = FastCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, cm->hbmx, size_limit, &sc);
  if(status == eslOK) { /* FastCYKScanHB() successfully finished */
    if(sc > cutoff) { /* score exceeds cutoff, we're done, reset search_opts and return */
      if(turn_qdb_back_on)        cm->search_opts &= ~CM_SEARCH_NOQDB; 
      if(turn_hbanded_back_off) { cm->search_opts &= ~CM_SEARCH_HBANDED; cm->tau = orig_tau; }
      ESL_DASSERT1((cm->flags       == init_flags));
      ESL_DASSERT1((cm->search_opts == init_search_opts));
      *ret_sc = sc;
      return eslOK;
    }
  }
  else if (status != eslERANGE) return status; /* else if status == eslERANGE, FastCYKScanHB() couldn't grow its DP matrix big enough, move onto next stage */

  /* stage 2 */
  cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
  cm->tau = 1e-10;
  if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, TRUE, 0)) != eslOK) return status;
  status = FastCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, cm->hbmx, size_limit, &sc);
  if(status == eslOK) { /* FastCYKScanHB() successfully finished */
    if(sc > cutoff) { /* score exceeds cutoff, we're done, reset search_opts and return */
      if(turn_qdb_back_on)             cm->search_opts &= ~CM_SEARCH_NOQDB; 
      if(turn_hbanded_back_off)      { cm->search_opts &= ~CM_SEARCH_HBANDED;      cm->tau = orig_tau; }
      if(turn_hmmscanbands_back_off) { cm->search_opts &= ~CM_SEARCH_HMMSCANBANDS; cm->tau = orig_tau; }
      ESL_DASSERT1((cm->flags       == init_flags));
      ESL_DASSERT1((cm->search_opts == init_search_opts));
      *ret_sc = sc;
      return eslOK;
    }
  }
  else if (status != eslERANGE) return status; /* else if status == eslERANGE, FastCYKScanHB() couldn't grow its DP matrix big enough, move onto next stage */

  /* stage 3, use 'default' dmin, dmax (which could be NULL) CYK or Inside */
  cm->search_opts &= ~CM_SEARCH_HBANDED;
  cm->search_opts &= ~CM_SEARCH_HMMSCANBANDS;
  if(turn_qdb_back_on) cm->search_opts &= ~CM_SEARCH_NOQDB; 

  if(cm->search_opts & CM_SEARCH_INSIDE) {
    if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &sc)) != eslOK) return status;
  }
  else { 
    if((status = FastCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &sc)) != eslOK) return status;
  }
  if(!turn_hbanded_back_off)      { cm->search_opts |= CM_SEARCH_HBANDED;      cm->tau = orig_tau; }
  if(!turn_hmmscanbands_back_off) { cm->search_opts |= CM_SEARCH_HMMSCANBANDS; cm->tau = orig_tau; }
  ESL_DASSERT1((cm->flags       == init_flags));
  ESL_DASSERT1((cm->search_opts == init_search_opts));

  /*if(sc > cutoff) { printf("3 sc: %10.4f\n", sc); }*/
  *ret_sc = sc;
  return eslOK;
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
  float cyk_megacalcs_per_sec = 275.;
  float ins_megacalcs_per_sec =  75.;
  float fwd_megacalcs_per_sec = 175.;
  float vit_megacalcs_per_sec = 380.;
  
  float seconds = 0.;

  if(! esl_opt_IsDefault(go, "--gumL")) L = ESL_MAX(L, esl_opt_GetInteger(go, "--gumL")); /* minimum L we allow is 2 * cm->W (L is sent into this func as 2 * cm->W), this is enforced silently (!) */

  switch(gum_mode) { 
  case GUM_CM_LC: 
  case GUM_CM_GC: 
    seconds = cfg->full_vcalcs[0] * (float) L * (float) nseq / cyk_megacalcs_per_sec;
    break;
  case GUM_CM_LI:
  case GUM_CM_GI:
    seconds = cfg->full_vcalcs[0] * (float) L * (float) nseq / ins_megacalcs_per_sec;
    break;
  case GUM_CP9_LV: 
  case GUM_CP9_GV: 
    seconds = cfg->hsi->full_cp9_ncalcs * (float) L * (float) nseq / vit_megacalcs_per_sec;
    break;
  case GUM_CP9_LF: 
  case GUM_CP9_GF: 
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
 * a partition start <pstart>, and every even token is 
 * a parititon end <pend>. First <pstart> must be 0,
 * other <pstart>s must be 1 more than previous
 * <pend>. The last <pend> must be 100, other <pends>
 * must be 1 less than following <pstart>.
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
 * set up the cfg->np and cfg->pstart data.
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

  /* create cfg->pstart */
  ESL_DASSERT1((cfg->pstart == NULL));
  ESL_ALLOC(cfg->pstart, sizeof(int) * nread);
  for(p = 0; p < nread; p++) cfg->pstart[p] = begin[p];
  free(begin);
  cfg->np = nread;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return status;
}


/* Function: update_avg_hit_len()
 * Date:     EPN, Sun Dec  9 15:50:39 2007
 * 
 * Purpose:  Calculate the average subseq length rooted at each state
 *           using the QDB calculation.
 *
 * Returns:  eslOK on success;
 */
int
update_avg_hit_len(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int safe_windowlen;
  float *avglen = NULL;

  safe_windowlen = cm->W * 2;
  while(!(BandCalculationEngine(cm, safe_windowlen, 1E-15, TRUE, NULL, NULL, NULL, &avglen))) {
    safe_windowlen *= 2;
    if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "update_avg_hit_len(), band calculation safe_windowlen big: %d\n", safe_windowlen);
  }
  if(cfg->avglen != NULL) free(cfg->avglen);
  cfg->avglen = avglen;
  return eslOK;
}

/* Function: switch_global_to_local()
 * Incept:   EPN, Mon Dec 10 08:43:32 2007
 * 
 * Purpose:  Switch a CM and it's CP9 HMM from global configuration
 *           to local configuration. Purposefully a local static function 
 *           in cmcalibrate.c, b/c we don't check if CM is in rsearch mode
 *           or any other jazz that'll never happen in cmcalibrate.
 *
 * Args:      cm - the model
 *
 * Returns:   eslOK on succes, othewise some other easel status code and
 *            errbuf is filled with error message.
 */
int 
switch_global_to_local(CM_t *cm, char *errbuf)
{
  if(cm->flags & CMH_LOCAL_BEGIN) ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CMH_LOCAL_BEGIN flag already raised.\n");
  if(cm->flags & CMH_LOCAL_END)   ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CMH_LOCAL_END flag already raised.\n");
  if(! (cm->flags & CMH_CP9))     ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CMH_CP9 flag down.\n");
  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CPLAN9_LOCAL_BEGIN flag already raised.\n");
  if(cm->cp9->flags & CPLAN9_LOCAL_END)   ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CPLAN9_LOCAL_END flag already raised.\n");
  if(cm->cp9->flags & CPLAN9_EL)          ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CPLAN9_EL flag already raised.\n");

  /* ConfigLocal() puts CM in local mode, recalcs QDBs (if they exist), remakes cm's scan matrix, 
   * logoddsifies CM, and makes inserts equiprobable (if nec) */
  ConfigLocal(cm, cm->pbegin, cm->pend); 
  /* CPlan9SWConfig() configures CP9 for local alignment, then logoddisfies CP9 (wastefully in this case) */
  CPlan9SWConfig(cm->cp9, cm->pbegin, cm->pbegin); 
  /* CPlan9ELConfig() configures CP9 for CM EL local ends, then logoddisfies CP9 */
  CPlan9ELConfig(cm);
  return eslOK;
}

/* Function: predict_cp9_filter_speedup()
 * Date:     EPN, Mon Dec 10 11:55:24 2007
 *
 * Purpose:  Given a CM and scores for a CP9 Viterbi and Forward scan
 *           of target seqs predict the speedup with an HMM filter, Forward and
 *           Viterbi, then update a BestFilterInfo_t object to 
 *           hold info on the faster of the two.
 *            
 * Args:     go  - command line options
 *           cfg - cmcalibrate's cfg object, mucho data (probably too much)
 *           errbuf - for printing error messages
 *           cm - the model
 *           fil_vit_cp9scA - [0..i..filN-1] best Viterbi score in sequence i 
 *           fil_fwd_cp9scA - [0..i..filN-1] best Foward score in sequence i 
 *           fil_partA      - [0..i..filN-1] partition of sequence i 
 *           bf             - BestFilterInfo_t object, we'll update this to hold info on a Viterbi or Forward filter strategy 
 *
 * Returns:  Updates BestFilterInfo_t object <bf> to hold info on fastest HMM filter, Viterbi or Forward
 *           eslOK on success;
 *           Other easel status code on an error with errbuf filled with error message.
 */
int
predict_cp9_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_vit_cp9scA, float *fil_fwd_cp9scA, int *fil_partA, BestFilterInfo_t *bf)
{
  int    status;
  float *sorted_fil_vit_EA;       /* sorted Viterbi E-values, so we can easily choose a threshold */
  float *sorted_fil_fwd_EA;       /* sorted Forward E-values, so we can easily choose a threshold */
  float  vit_E, fwd_E;            /* a Viterbi and Forward E value */
  int    cp9_vit_mode, cp9_fwd_mode; /* a Viterbi, Forward Gumbel mode, respectively  */
  int    evalue_L;                /* database length used for calcing E-values in CP9 gumbels from cfg->cmstats */
  float  fil_calcs;               /* number of million dp calcs predicted for the HMM filter scan */
  float  vit_surv_calcs;          /* number of million dp calcs predicted for the CM scan of Viterbi filter survivors */
  float  fwd_surv_calcs;          /* number of million dp calcs predicted for the CM scan of Forward filter survivors */
  float  vit_fil_plus_surv_calcs; /* Viterbi filter calcs plus survival calcs */ 
  float  fwd_fil_plus_surv_calcs; /* Foward filter calcs plus survival calcs */ 
  float  nonfil_calcs;            /* number of million dp calcs predicted for a full CM scan */
  float  vit_spdup, fwd_spdup;    /* predicted speedups for Viterbi and Forward */
  int    i, p;                    /* counters */
  int    cmi = cfg->ncm-1;        /* CM index we're on */
  int    Fidx;                    /* index in sorted scores that threshold will be set at (1-F) * N */
  float  filE = esl_opt_GetReal(go, "--filE"); /* E-value cutoff for accepting CM hits for filter test */
  float  F = esl_opt_GetReal(go, "--F"); /* fraction of CM seqs we require filter to let pass */
  int    filN  = esl_opt_GetInteger(go, "--filN"); /* number of sequences we emitted from CM for filter test */
  float  Starg   = esl_opt_GetReal(go, "--starg"); /* target survival fraction */
  float  E_Starg;                 /* E-value threshold that exactly satisifies Starg survival fraction */
  float  vit_E_F1;                /* viterbi E value if F were == 1.0 */
  float  fwd_E_F1;                /* forward E value if F were == 1.0 */
  float  surv_res_per_hit;        /* expected number of residues to survive filter from DB for each hit 2*W-avglen[0], twice W minus the average lenght of a hit from QDB calc */
  float  Smin;                    /* minimally useful survival fraction, any less than this and our filter would be (predicted to be) doing more than 10X the work of the CM */
  float  E_min;                   /* E-value threshold that exactly satisifies Smin survival fraction */
  float  vit_surv_fract;          /* predicted survival fraction for Viterbi filter */
  float  fwd_surv_fract;          /* predicted survival fraction for Forward filter */
  
  cp9_vit_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? GUM_CP9_LV : GUM_CP9_GV;
  cp9_fwd_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? GUM_CP9_LF : GUM_CP9_GF;

  /* contract checks */
  if(! (cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), gumbel stats for CP9 viterbi mode: %d are not valid.\n", cp9_vit_mode);
  if(! (cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), gumbel stats for CP9 forward mode: %d are not valid.\n", cp9_fwd_mode);
  if(cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->L != cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), db length for gumbel stats for CP9 viterbi (%d) and forward (%d) differ.\n", cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->L, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->L);

  evalue_L = cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->L;

  /* contract checks specific to case when there is more than 1 partition */
  if(cfg->cmstatsA[cfg->ncm-1]->np != 1) { 
    for(p = 1; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) {
      if(evalue_L != cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), partition %d db length (%d) for Viterbi gumbel stats differ than from partition 1 Viterbi db length (%d).\n", p, cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->L, evalue_L);
      if(evalue_L != cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), partition %d db length (%d) for Forward gumbel stats differ than from partition 1 Viterbi db length (%d).\n", p, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->L, evalue_L);
    }
  }

  Fidx  = (int) ((1. - F) * (float) filN);

  /* convert bit scores to E-values and sort them */
  ESL_ALLOC(sorted_fil_vit_EA, sizeof(float) * filN);
  ESL_ALLOC(sorted_fil_fwd_EA, sizeof(float) * filN);
  for(i = 0; i < filN; i++) { 
    p = fil_partA[i];
    sorted_fil_vit_EA[i] = RJK_ExtremeValueE(fil_vit_cp9scA[i], cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->mu, cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->lambda); 
    sorted_fil_fwd_EA[i] = RJK_ExtremeValueE(fil_fwd_cp9scA[i], cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->mu, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->lambda); 
  }
  esl_vec_FSortDecreasing(sorted_fil_vit_EA, filN);
  vit_E = sorted_fil_vit_EA[Fidx];
  esl_vec_FSortDecreasing(sorted_fil_fwd_EA, filN);
  fwd_E = sorted_fil_fwd_EA[Fidx];

  /* now vit_E and fwd_E are expected number of CP9 Viterbi/Forward hits with score above threshold 
   * (Fidx'th best score) in a sequence DB of length evalue_L, convert that DB size to cfg->dbsize */
  vit_E *= cfg->dbsize / evalue_L;
  fwd_E *= cfg->dbsize / evalue_L;

  fil_calcs  = cfg->hsi->full_cp9_ncalcs; /* fil_calcs is millions of DP calcs for CP9 scan of 1 residue */
  fil_calcs *= cfg->dbsize;               /* fil_calcs is millions of DP calcs for CP9 scan of length cfg->dbsize */
  nonfil_calcs = cfg->hsi->full_cm_ncalcs;      /* total number of millions of DP calculations for full CM scan of 1 residue */
  nonfil_calcs *= cfg->dbsize;                  /* now nonfil-calcs corresponds to cfg->dbsize */

  /* determine our thresholds:
   * 1. if vit_E and/or fwd_E yield predicted survival fractions that are less than Starg:
   *    rewrite vit_E or fwd_E as the maximum E-value threshold that sees ALL hits (when F==1.0,
   *    this is {vit,fwd}_E_F1 below), and the E-value threshold that exactly satisfies Starg (E_Starg below).
   * 2. if after 1, vit_E and/or fwd_E still yield predicted survival fractions less than Smin,
   *    the survival fraction at which the number of DP calcs for the survivors of the filter is only
   *    10% the number of dp calcs for the filter, then we set vit_E or fwd_E to the E value that
   *    exactly satisfies Smin.
   */
  surv_res_per_hit = (2. * cm->W - (cfg->avglen[0])); /* avg length of surviving fraction of db from a single hit (cfg->avglen[0] is avg subseq len in subtree rooted at v==0, from QDB calculation) */
  vit_E_F1= (sorted_fil_vit_EA[0] * (cfg->dbsize / evalue_L));
  fwd_E_F1= (sorted_fil_fwd_EA[0] * (cfg->dbsize / evalue_L));
  E_Starg = (Starg * (float) cfg->dbsize) / surv_res_per_hit;
  Smin    = fil_calcs / (10. * nonfil_calcs);
  E_min   = (Smin * (float) cfg->dbsize) / surv_res_per_hit;
  if(Starg < Smin) { /* we never go less than Smin */
    Starg = Smin;
    E_Starg = E_min;
  }
  
  vit_surv_fract = (vit_E * surv_res_per_hit) / (float) cfg->dbsize;
  fwd_surv_fract = (fwd_E * surv_res_per_hit) / (float) cfg->dbsize;
  
  ESL_DPRINTF1(("vit_E:    %.5f\n", vit_E));
  ESL_DPRINTF1(("vit_surv: %.10f\n", vit_surv_fract));
  ESL_DPRINTF1(("fwd_E:    %.5f\n", fwd_E));
  ESL_DPRINTF1(("fwd_surv: %.10f\n", fwd_surv_fract));

  ESL_DPRINTF1(("vit_E_F1: %.5f\n", vit_E_F1));
  ESL_DPRINTF1(("fwd_E_F1: %.5f\n", fwd_E_F1));
  ESL_DPRINTF1(("Starg:    %.5f\n", Starg));
  ESL_DPRINTF1(("E_Starg:  %.5f\n", E_Starg));
  ESL_DPRINTF1(("Smin:     %.5f\n", Smin));
  ESL_DPRINTF1(("E_min:    %.5f\n", E_min));

  if(vit_surv_fract < Starg) { 
    if(vit_E_F1 < E_Starg) { 
      vit_E = vit_E_F1;
      if(cfg->be_verbose) printf("set vit_E as vit_E_F1: %.5f\n", vit_E);
      else { ESL_DPRINTF1(("set vit_E as vit_E_F1: %.5f\n", vit_E)); }
    }
    else { 
      vit_E = E_Starg;
      if(cfg->be_verbose) printf("set vit_E as E_Starg: %.5f\n", vit_E);
      else { ESL_DPRINTF1(("set vit_E as E_Starg: %.5f\n", vit_E)); }
    }
    if(vit_E < E_min) {
      vit_E = E_min;
      if(cfg->be_verbose) printf("set vit_E as E_min: %.5f\n", vit_E);
      else { ESL_DPRINTF1(("set vit_E as E_min: %.5f\n", vit_E)); }
    }      
  }

  if(fwd_surv_fract < Starg) { 
    if(fwd_E_F1 < E_Starg) { 
      fwd_E = fwd_E_F1;
      if(cfg->be_verbose) printf("set fwd_E as fwd_E_F1: %.5f\n", fwd_E);
      else { ESL_DPRINTF1(("set fwd_E as fwd_E_F1: %.5f\n", fwd_E)); }
    }
    else { 
      fwd_E = E_Starg;
      if(cfg->be_verbose) printf("set fwd_E as E_Starg: %.5f\n", fwd_E);
      else { ESL_DPRINTF1(("set fwd_E as E_Starg: %.5f\n", fwd_E)); }
    }
    if(fwd_E < E_min) {
      fwd_E = E_min;
      if(cfg->be_verbose) printf("set fwd_E as E_min: %.5f\n", fwd_E);
      else { ESL_DPRINTF1(("set fwd_E as E_min: %.5f\n", fwd_E)); }
    }      
  }

  for(i = 0; i < filN; i++) ESL_DPRINTF1(("HMM i: %4d vit E: %15.10f fwd E: %15.10f\n", i, sorted_fil_vit_EA[i], sorted_fil_fwd_EA[i]));

  /* calculate speedup for Viterbi */
  vit_surv_calcs = vit_E *     /* number of hits expected to survive filter */
    surv_res_per_hit *         /* avg length of surviving fraction of db from a single hit */
    cfg->hsi->full_cm_ncalcs;  /* number of calculations for full CM scan of 1 residue */

  fwd_surv_calcs = fwd_E *     /* number of hits expected to survive filter */
    surv_res_per_hit *         /* avg length of surviving fraction of db from a single hit */
    cfg->hsi->full_cm_ncalcs;  /* number of calculations for full CM scan of 1 residue */

  vit_fil_plus_surv_calcs = fil_calcs + vit_surv_calcs; /* total number of millions of DP calculations expected using the CP9 viterbi filter for scan of cfg->dbsize */
  vit_spdup = nonfil_calcs / vit_fil_plus_surv_calcs;
  fwd_fil_plus_surv_calcs = (fil_calcs * 2.) + fwd_surv_calcs; /* total number of millions of DP calculations expected using the CP9 forward filter for scan of cfg->dbsize (logsum corrected, Forward calcs *= 2.) */
  fwd_spdup = nonfil_calcs / fwd_fil_plus_surv_calcs;
  /* We multiply number of forward calculations by 2.0 to correct for the fact that Forward takes about 2X as long as Viterbi, b/c it requires logsum operations instead of ESL_MAX's,
   * so we factor this in when calc'ing the predicted speedup. */

  if(cfg->be_verbose) { 
    printf("\nHMM(vit) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", vit_E, fil_calcs, vit_surv_calcs, vit_fil_plus_surv_calcs, nonfil_calcs, vit_spdup);
    printf("HMM(fwd) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", fwd_E, fil_calcs, fwd_surv_calcs, fwd_fil_plus_surv_calcs, nonfil_calcs, fwd_spdup);
  }
  else {
    ESL_DPRINTF1(("\nHMM(vit) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", vit_E, fil_calcs, vit_surv_calcs, vit_fil_plus_surv_calcs, nonfil_calcs, vit_spdup));
    ESL_DPRINTF1(("HMM(fwd) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", fwd_E, fil_calcs, fwd_surv_calcs, fwd_fil_plus_surv_calcs, nonfil_calcs, fwd_spdup));
  }

  if(esl_opt_GetBoolean(go, "--fviterbi")) { /* user specified Viterbi */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_VITERBI, vit_E, fil_calcs, vit_fil_plus_surv_calcs)) != eslOK) return status;
  }
  else if(esl_opt_GetBoolean(go, "--fforward")) { /* user specified Forward */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_FORWARD, fwd_E, fil_calcs, fwd_fil_plus_surv_calcs)) != eslOK) return status;
  }
  else if (vit_spdup > fwd_spdup) { /* Viterbi is winner */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_VITERBI, vit_E, fil_calcs, vit_fil_plus_surv_calcs)) != eslOK) return status;
  }
  else { /* Forward is winner */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_FORWARD, fwd_E, fil_calcs, fwd_fil_plus_surv_calcs)) != eslOK) return status;
  }
  return eslOK;

 ERROR:
  return status; 
}

/* Function: predict_hybrid_filter_speedup()
 * Date:     EPN, Tue Dec 11 04:56:39 2007
 *
 * Purpose:  Given a CM and scores for a hybrid CYK/Viterbi scan
 *           of target seqs predict the speedup with a hybrid filter,
 *           then if it's faster than the existing best filter in
 *           BestFilterInfo_t object <bf>, update <bf> to hold info 
 *           on the hybrid filter.
 *            
 * Args:     go  - command line options
 *           cfg - cmcalibrate's cfg object, mucho data (probably too much)
 *           errbuf - for printing error messages
 *           cm - the model
 *           fil_hybscA     - [0..i..filN-1] best Foward score in sequence i 
 *           gum_hybA       - [0..cfg->np]   hybrid gumbels for each partition
 *           fil_partA      - [0..i..filN-1] partition of sequence i 
 *           bf             - BestFilterInfo_t object, we'll update this to hold info on a Viterbi or Forward filter strategy 
 * 
 * Returns:  possibly updates BestFilterInfo_t object <bf> to hold info on hybrid filter
 *           eslOK on success;
 *           Other easel status code on an error with errbuf filled with error message.
 *           <ret_getting_faster> set to TRUE if hybrid scanner replaced previous best filter,
 *           FALSE if not.
 */
int
predict_hybrid_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_hybscA, int *fil_partA, GumbelInfo_t **gum_hybA, BestFilterInfo_t *bf, int *ret_getting_faster)
{
  int    status;
  float  *sorted_fil_hybEA;       /* sorted hybrid E-values, so we can easily choose a threshold */
  float  E;                       /* E-value */
  int    evalue_L;                /* length used for calc'ing E values */
  float  fil_calcs;               /* number of million dp calcs predicted for the hybrid scan */
  float  surv_calcs;              /* number of million dp calcs predicted for the CM scan of filter survivors */
  float  fil_plus_surv_calcs;     /* filter calcs plus survival calcs */ 
  float  nonfil_calcs;            /* number of million dp calcs predicted for a full CM scan */
  float  spdup;                   /* predicted speedups for Viterbi and Forward, and a temporary one */
  int    i, p;                    /* counters */
  int    Fidx;                     /* index in sorted scores that threshold will be set at (1-F) * N */
  float  F = esl_opt_GetReal(go, "--F"); /* fraction of CM seqs we require filter to let pass */
  float  filE = esl_opt_GetReal(go, "--filE"); /* E-value of cutoff for accepting CM seqs */
  int    filN  = esl_opt_GetInteger(go, "--filN"); /* number of sequences we emitted from CM for filter test */

  /* contract checks */
  if(! (gum_hybA[0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_hybrid_filter_speedup(), gumbel stats for hybrid scanner are not valid.\n");
  evalue_L = gum_hybA[0]->L;
  /* contract checks specific to case when there is more than 1 partition */
  if(cfg->cmstatsA[cfg->ncm-1]->np != 1) { 
    for(p = 1; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) {
      if(evalue_L != gum_hybA[p]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_hybrid_filter_speedup(), partition %d db length (%d) for hybrid gumbel stats differ than from partition 1 hybrid db length (%d).\n", p, gum_hybA[p]->L, evalue_L);
    }
  }

  Fidx  = (int) ((1. - F) * (float) filN);

  /* convert bit scores to E-values and sort them */
  ESL_ALLOC(sorted_fil_hybEA, sizeof(float) * filN);
  for(i = 0; i < filN; i++) { 
    p = fil_partA[i];
    sorted_fil_hybEA[i] = RJK_ExtremeValueE(fil_hybscA[i], gum_hybA[p]->mu, gum_hybA[p]->lambda); 
  }
  esl_vec_FSortDecreasing(sorted_fil_hybEA, filN);
  E = sorted_fil_hybEA[Fidx];

  for(i = 0; i < filN; i++) ESL_DPRINTF1(("HYBRID i: %4d E: %10.4f\n", i, sorted_fil_hybEA[i]));
  
  /* calculate speedup */
  /* E is expected number of hybrid hits with score above threshold (Fidx'th best score) at least sc in sequence DB of length evalue_L */
  E *= cfg->dbsize / evalue_L;
  /* E is now expected number of CP9 Viterbi or Forward hits with score above threshold in sequence DB of length cfg->dbsize */
  fil_calcs  = cfg->hsi->hybrid_ncalcs; /* fil_calcs is millions of DP calcs for hybrid scan of 1 residue */
  fil_calcs *= cfg->dbsize;             /* fil_calcs is millions of DP calcs for hybrid scan of length cfg->dbsize */
  surv_calcs = E *     /* number of hits expected to survive filter */
    (2. * cm->W - (cfg->avglen[0])) * /* average length of surviving fraction of db from a single hit (cfg->avglen[0] is avg subseq len in  subtree rooted at v==0, from QDB calculation, so slightly inappropriate b/c we're concerned with hybrid hits here) */
    cfg->hsi->full_cm_ncalcs; /* number of calculations for full CM scan of 1 residue */
  fil_plus_surv_calcs = fil_calcs  + surv_calcs; /* total number of millions of DP calculations expected using the hybrid filter for scan of cfg->dbsize (logsum corrected, Forward calcs *= 2.) */
  nonfil_calcs = cfg->hsi->full_cm_ncalcs;      /* total number of millions of DP calculations for full CM scan of 1 residue */
  nonfil_calcs *= cfg->dbsize;                  /* now nonfil-calcs corresponds to cfg->dbsize */
  spdup = nonfil_calcs / fil_plus_surv_calcs;
  
  if(cfg->be_verbose) printf("HYBRID E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup);
  else { ESL_DPRINTF1(("HYBRID E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup)); }

  if(spdup > (bf->full_cm_ncalcs / bf->fil_plus_surv_ncalcs)) { /* hybrid is best filter strategy so far */
    if((status = SetBestFilterInfoHybrid(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, E, fil_calcs, fil_plus_surv_calcs, cfg->hsi, cfg->cmstatsA[cfg->ncm-1]->np, gum_hybA)) != eslOK) return status;
    *ret_getting_faster = TRUE;
  }
  else *ret_getting_faster = FALSE;

  return eslOK;

 ERROR:
  return status; 
}


/* Function: predict_best_sub_cm_roots()
 * Date:     EPN, Mon Dec 10 15:56:00 2007
 *
 * Purpose:  Given a CM and scores for a CM scan of target seqs
 *           predict the best sub CM roots we could use to 
 *           filter with.
 *            
 * Returns:  eslOK on success;
 *           Other status code on error, with error message in errbuf.
 */
int 
predict_best_sub_cm_roots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **fil_vscAA, int **ret_sorted_best_roots_v)
{
  int    status;
  float **sorted_fil_vscAA;       /* [0..v..cm->M-1][0..filN-1] best score for each state v, each target seq */
  float  fil_calcs;               /* number of million dp calcs predicted for the HMM filter scan */
  float  surv_calcs;              /* number of million dp calcs predicted for the CM scan of filter survivors */
  float  fil_plus_surv_calcs;     /* filter calcs plus survival calcs */ 
  float  nonfil_calcs;            /* number of million dp calcs predicted for a full CM scan */
  float  spdup;                   /* predicted speedup a sub CM filter */
  float  E, tmp_E;                /* E value */
  float  sc;                      /* bit score */
  int    i, p, s, v;              /* counters */
  int    Fidx;                    /* index in sorted scores that threshold will be set at (1-F) * N */
  float  F = esl_opt_GetReal(go, "--F"); /* fraction of CM seqs we require filter to let pass */
  int    filN  = esl_opt_GetInteger(go, "--filN"); /* number of sequences we emitted from CM for filter test */
  int    nstarts;                  /* # start states (and start groups) in the CM, from cfg->hsi */                                 
  int   *best_per_start_v;         /* sub CM filter state v that gives best speedup per start group */
  float *best_per_start_spdup;     /* best sub CM filter state speedup per start group */

  if(ret_sorted_best_roots_v == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_best_sub_cm_roots, ret_sorted_best_roots_v == NULL.\n");

  Fidx  = (int) ((1. - F) * (float) filN);

  ESL_ALLOC(sorted_fil_vscAA, sizeof(float *) * cm->M);
  /*ESL_ALLOC(sorted_fil_EAA, sizeof(float *) * cm->M);*/

  nstarts = cfg->hsi->nstarts;
  ESL_ALLOC(best_per_start_v,     sizeof(int)   * nstarts);
  ESL_ALLOC(best_per_start_spdup, sizeof(float) * nstarts);
  for(s = 0; s < nstarts; s++) {
    best_per_start_v[s] = -1;
    best_per_start_spdup[s] = -eslINFINITY;
  }

  for(v = 0; v < cm->M; v++) {
    ESL_ALLOC(sorted_fil_vscAA[v], sizeof(float) * filN);
    esl_vec_FCopy(fil_vscAA[v], filN, sorted_fil_vscAA[v]); 
    esl_vec_FSortIncreasing(sorted_fil_vscAA[v], filN);
  }
  ESL_DPRINTF1(("vscAA[0] scores:\n"));
  for(i = 0; i < filN; i++) ESL_DPRINTF1(("i: %4d sc: %10.4f\n", i, sorted_fil_vscAA[0][i]));

  for(v = 0; v < cm->M; v++) {
    if(cfg->hsi->iscandA[v]) {
      sc = sorted_fil_vscAA[v][Fidx];
      /* set E as E-value for sc from partition that gives sc the lowest E-value (conservative, sc will be at least as significant as E across all partitions) */
      E  = RJK_ExtremeValueE(sc, cfg->vmuAA[0][v], cfg->vlambdaAA[0][v]); 
      for(p = 1; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) {
	tmp_E = RJK_ExtremeValueE(sc, cfg->vmuAA[p][v], cfg->vlambdaAA[p][v]); 
	if(tmp_E < E) E = tmp_E;
      }
      /* E is now expected number of hits for db of cfg->length 2 * cm->W */
      E *= cfg->dbsize / (cm->W * 2.);
      /* E is now expected number of hits for db of cfg->dbsize */

      fil_calcs   = cfg->hsi->cm_vcalcs[v]; /* fil_calcs is millions of DP calcs for sub CM (root = v) scan of 1 residue */
      fil_calcs  *= cfg->dbsize;            /* fil_calcs is millions of DP calcs for sub CM (root = v) scan of length cfg->dbsize */
      surv_calcs = E *     /* number of hits expected to survive filter */
	(2. * cm->W - (cfg->avglen[v])) * /* average length of surviving fraction of db from a single hit (cfg->avglen[v] is avg subseq len in  subtree rooted at v */
	cfg->hsi->full_cm_ncalcs; /* number of calculations for full CM scan of 1 residue */
      fil_plus_surv_calcs = fil_calcs + surv_calcs;
      nonfil_calcs = cfg->hsi->full_cm_ncalcs;      /* total number of millions of DP calculations for full CM scan of 1 residue */
      nonfil_calcs *= cfg->dbsize;                  /* now nonfil-calcs corresponds to cfg->dbsize */
      spdup = nonfil_calcs / fil_plus_surv_calcs;
      if(cfg->be_verbose) { printf("SUB %3d sg: %2d sc: %10.4f E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", v, cfg->hsi->startA[v], sc, E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup); } 
      else                { ESL_DPRINTF1(("SUB %3d sg: %2d sc: %10.4f E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", v, cfg->hsi->startA[v], sc, E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup)); } 
      s = cfg->hsi->startA[v];
      if(spdup > best_per_start_spdup[s] && cm->ndidx[v] != 0) { /* can't filter with a state in node 0 */
	best_per_start_v[s]     = v;
	best_per_start_spdup[s] = spdup;
      }	
    }  
  }
  for(s = 0; s < nstarts; s++) { 
    if(cfg->be_verbose) printf("START %d v: %d spdup: %10.4f\n", s, best_per_start_v[s], best_per_start_spdup[s]);
    else {       ESL_DPRINTF1(("START %d v: %d spdup: %10.4f\n", s, best_per_start_v[s], best_per_start_spdup[s])); }
    ESL_DASSERT1((best_per_start_v[s] != -1));
  }

  /* sort the best sub CM roots (1 per start group) by their speedup,
   * this is an embarassing N^2 sorting, but biggest RNAs have ~ 100 starts, so this is okay I guess (LSU has ~140 starts) 
   */
  int *sorted_best_roots_v; 
  int *sorted_best_roots_start; 
  float *sorted_best_roots_spdup;
  int *already_chosen;
  int s1, s2;
  int best_cur_v;
  int best_cur_start;
  float best_cur_spdup;

  ESL_ALLOC(sorted_best_roots_v,     sizeof(int) * nstarts);
  ESL_ALLOC(sorted_best_roots_start, sizeof(int) * nstarts);
  ESL_ALLOC(sorted_best_roots_spdup, sizeof(int) * nstarts);
  ESL_ALLOC(already_chosen,          sizeof(int) * nstarts);
  esl_vec_ISet(already_chosen, nstarts, FALSE);
  for(s1 = 0; s1 < nstarts; s1++) {
    best_cur_v = -1;
    best_cur_start = -1;
    best_cur_spdup = -eslINFINITY;
    for(s2 = 0; s2 < nstarts; s2++) { 
      if(! already_chosen[s2]) {
	if(best_per_start_spdup[s2] > best_cur_spdup) { 
	  best_cur_v = best_per_start_v[s2];
	  best_cur_start = s2;
	  best_cur_spdup = best_per_start_spdup[s2];
	}
      }
    }
    sorted_best_roots_v[s1] = best_cur_v;
    sorted_best_roots_start[s1] = best_cur_start;
    sorted_best_roots_spdup[s1] = best_cur_spdup;
    already_chosen[best_cur_start] = TRUE;
  }
  for(s1 = 0; s1 < nstarts; s1++) {
    if(cfg->be_verbose) printf("SORTED rank: %d v: %d spdup: %.5f start: %d\n", s1, sorted_best_roots_v[s1], sorted_best_roots_spdup[s1], sorted_best_roots_start[s1]);
    else { ESL_DPRINTF1(("SORTED rank: %d v: %d spdup: %.5f start: %d\n", s1, sorted_best_roots_v[s1], sorted_best_roots_spdup[s1], sorted_best_roots_start[s1])); } 
    ESL_DASSERT1((sorted_best_roots_v[s1] != -1));
  }
  *ret_sorted_best_roots_v = sorted_best_roots_v;

  free(sorted_best_roots_start);
  free(sorted_best_roots_spdup);
  free(already_chosen);
  free(best_per_start_v);
  free(best_per_start_spdup);
  for(v = 0; v < cm->M; v++) free(sorted_fil_vscAA[v]);
  free(sorted_fil_vscAA);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "predict_best_sub_cm_roots(), memory allocation error.");
  return status; /* NEVERREACHED */
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
  long seed;
  long temp;
  int  seedlen;
  char *seedstr;

  if(cfg->ccom  != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->ccom  is non-NULL.");
  if(cfg->cdate != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->cdate is non-NULL.");
  
  
  /* Set the cmbuild command info, the cfg->clog->bcom string */
  for (i = 0; i < go->optind; i++) { /* copy all command line options, but not the command line args yet, we may need to append '-s ' before the args */
    esl_strcat(&(cfg->ccom),  -1, go->argv[i], -1);
    esl_strcat(&(cfg->ccom),  -1, " ", 1);
  }
  /* if -s NOT enabled, we need to append the seed info also */
  seed = esl_randomness_GetSeed(cfg->r);
  if(esl_opt_IsDefault(go, "-s")) {
    if(cfg->r == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "get_cmcalibrate_comlog_info(), cfg->r is NULL but --gibbs enabled, shouldn't happen.");
    temp = seed; 
    seedlen = 1; 
    while(temp > 0) { temp/=10; seedlen++; } /* determine length of stringized version of seed */
    seedlen += 4; /* strlen(' -s ') */
    ESL_ALLOC(seedstr, sizeof(char) * (seedlen+1));
    sprintf(seedstr, " -s %ld ", seed);
    esl_strcat((&cfg->ccom), -1, seedstr, seedlen);
  }
  else { /* -s was enabled, we'll do a sanity check */
    if(seed != (long) esl_opt_GetInteger(go, "-s")) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "get_cmcalibrate_comlog_info(), cfg->r's seed is %ld, but -s was enabled with argument: %ld!, this shouldn't happen.", seed, (long) esl_opt_GetInteger(go, "-s"));
  }

  for (i = go->optind; i < go->argc; i++) { /* copy all command line options, but not the command line args yet, we may need to append '-s ' before the args */
    esl_strcat(&(cfg->ccom), -1, go->argv[i], -1);
    if(i < (go->argc-1)) esl_strcat(&(cfg->ccom), -1, " ", 1);
  }

  /* Set the cmcalibrate call date, the cfg->clog->bdate string */
  time_t date = time(NULL);
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
  int which_comlog;       /* 1 or 2, 1 if we'll write to cm->ccom1 and cm->cdate1, 2 if we'll write to 
			   * cm->ccom2 and cm->cdate2. Can only be 2 if --filonly enabled, otherwise
			   * we're creating new gumbel stats and/or filter thresholds, so any previous info 
			   * in the CM file from previous cmcalibrate calls will be deleted/overwritten. 
			   */ 
  
  if(ccom  == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), ccom  is non-NULL.");
  if(cdate == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), cdate is non-NULL.");

  which_comlog = (esl_opt_GetBoolean(go, "--filonly")) ? 2 : 1; /* only way to write to comlog->ccom2, comlog->cdate2 is if --filonly was enabled */
  
  /* 2 possible cases, case 1: we're writing/overwriting cm->comlog->ccom1, cm->comlog->cdate1, 
   * case 2: we're writing/overwriting cm->comlog->ccom2, cm->comlog->cdate2 
   */

  if(which_comlog == 1) {
    /* free all cmcalibrate comlog info, we're about to overwrite any information that any previous cmcalibrate
     * call could have written to the cm file.
     */
    if(cm->comlog->ccom1  != NULL)  { free(cm->comlog->ccom1);  cm->comlog->ccom1 = NULL;  }
    if(cm->comlog->cdate1 != NULL)  { free(cm->comlog->cdate1); cm->comlog->cdate1 = NULL; }
    if(cm->comlog->ccom2  != NULL)  { free(cm->comlog->ccom2);  cm->comlog->ccom2 = NULL;  }
    if(cm->comlog->cdate2 != NULL)  { free(cm->comlog->cdate2); cm->comlog->cdate2 = NULL; }
    
    if((status = esl_strdup(ccom, -1, &(cm->comlog->ccom1)))  != eslOK) goto ERROR; 
    if((status = esl_strdup(cdate,-1, &(cm->comlog->cdate1))) != eslOK) goto ERROR; 
  }
  else { /* which_comlog == 2 */
    /* if it exists, free comlog info comlog->ccom2, comlog->cdate2,, we're about to overwrite the info
     * corresponding to that cmcalibrate call (the filter threshold information ONLY) 
     * First, assert we have comlog->ccom1, comlog->cdate1 */
    if(cm->comlog->ccom1  == NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), --filonly enabled, but cm->comlog->ccom1 is NULL.");
    if(cm->comlog->ccom1  == NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), --filonly enabled, but cm->comlog->cdate1 is NULL.");

    if(cm->comlog->ccom2  != NULL)  { free(cm->comlog->ccom2);  cm->comlog->ccom2 = NULL;  }
    if(cm->comlog->cdate2 != NULL)  { free(cm->comlog->cdate2); cm->comlog->cdate2 = NULL; }

    if((status = esl_strdup(ccom, -1, &(cm->comlog->ccom2)))  != eslOK) goto ERROR; 
    if((status = esl_strdup(cdate,-1, &(cm->comlog->cdate2))) != eslOK) goto ERROR; 
  }
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "update_comlog() error status: %d, probably out of memory.", status);
  return status; 
}


/* format_time_string()
 * Date:     SRE, Fri Nov 26 15:06:28 1999 [St. Louis]
 *
 * Purpose:  Given a number of seconds, format into
 *           hh:mm:ss.xx in a provided buffer.
 *
 * Args:     buf     - allocated space (128 is plenty!)
 *           sec     - number of seconds
 *           do_frac - TRUE (1) to include hundredths of a sec
 */
static void
format_time_string(char *buf, double sec, int do_frac)
{
  int h, m, s, hs;
  
  h  = (int) (sec / 3600.);
  m  = (int) (sec / 60.) - h * 60;
  s  = (int) (sec) - h * 3600 - m * 60;
  if (do_frac) {
    hs = (int) (sec * 100.) - h * 360000 - m * 6000 - s * 100;
    sprintf(buf, "%02d:%02d:%02d.%02d", h,m,s,hs);
  } else {
    sprintf(buf, "%02d:%02d:%02d", h,m,s);
  }
}

/* Function: print_cm_info
 * Date:     EPN, Tue Jan  8 05:51:47 2008
 *
 * Purpose:  Print per-CM info to stdout. 
 *
 * Returns:  eslOK on success
 */
static int
print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  printf("%-8s %s\n", "CM:",      cm->name); 
  if(esl_opt_GetBoolean(go, "--noqdb")) printf("%-8s %s %g\n", "beta:", "0.0 (--noqdb), for W calc: ", cm->beta); 
  else                                  printf("%-8s %g\n", "beta:", cm->beta); 
  printf("%-8s %s\n", "command:", cfg->ccom);
  printf("%-8s %s\n", "date:",    cfg->cdate);
  printf("%-8s %d\n", "nproc:",   (cfg->nproc == 0) ? 1 : cfg->nproc);
  printf("%-8s ", "fil cut:");
  if     (esl_opt_GetBoolean(go,"--gumonly")) printf("N/A (gumbel only)\n");
  else if(esl_opt_GetBoolean(go,"--nc"))      printf("NC (%.2f bits)\n", cm->nc);
  else if(esl_opt_GetBoolean(go,"--tc"))      printf("TC (%.2f bits)\n", cm->tc);
  else if(esl_opt_GetBoolean(go,"--ga"))      printf("GA (%.2f bits)\n", cm->ga);
  else if(esl_opt_GetBoolean(go,"--all"))     printf("none, all CM seqs accepted\n");
  else printf("E value (%g)\n", esl_opt_GetReal(go, "--filE"));
  printf("%-8s %.3f Mb\n", "db size:", (cfg->dbsize/1000000.));
  printf("\n");
  printf("%6s  %3s  %3s  %3s %5s %6s %5s %5s    percent complete         %10s\n", "",       "",     "",    "",     "",        "",     "",         "",  "");
  printf("%6s  %3s  %3s  %3s %5s %6s %5s %5s [5.......50.......100] %10s\n", "stage",  "mod",  "cfg", "alg",  "gumN",    "len",  "filN",     "cut", "elapsed");
  printf("%6s  %3s  %3s  %3s %5s %6s %5s %5s %22s %10s\n", "------", "---", "---", "---", "-----", "------", "-----", "-----", "----------------------", "----------");
  return eslOK;
}

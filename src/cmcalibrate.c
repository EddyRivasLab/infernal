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

#define MPI_FINISHED_GUMBEL -1 /* message to send to workers */
#define MPI_FINISHED_FILTER -2 /* message to send to workers */

#include "funcs.h"		/* external functions                   */
#include "structs.h"

static ESL_OPTIONS options[] = {
  /* name                type           default env   range     toggles      reqs       incomp  help  docgroup*/
  { "-h",                eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "show brief help on version and usage",   1 },
  { "-s",                eslARG_INT,    NULL,   NULL, "n>0",    NULL,        NULL,        NULL, "set random number generator seed to <n>",  1 },
  { "-t",                eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "print timings for gumbel fitting and HMM filter calculation",  1},
#ifdef HAVE_DEVOPTS 
  { "-v",                eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "print arguably interesting info",  1},
#endif
  /* 4 --p* options below are hopefully temporary b/c if we have E-values for the CM using a certain cm->pbegin, cm->pend,
   * changing those values in cmsearch invalidates the E-values, so we should pick hard-coded values for cm->pbegin cm->pend */
  { "--pebegin",         eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL,  "--pbegin","set all local begins as equiprobable", 1 },
  { "--pfend",           eslARG_REAL,   NULL,   NULL, "0<x<1",  NULL,        NULL,    "--pend", "set all local end probs to <x>", 1 },
  { "--pbegin",          eslARG_REAL,   "0.05", NULL, "0<x<1",  NULL,        NULL,        NULL, "set aggregate local begin prob to <x>", 1 },
  { "--pend",            eslARG_REAL,   "0.05", NULL, "0<x<1",  NULL,        NULL,        NULL, "set aggregate local end prob to <x>", 1 },

  /* options for gumbel estimation */
  { "--gum-only",       eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,"--fil-only", "only estimate gumbels, don't calculate HMM filter thresholds", 2},
  { "--gum-L",          eslARG_INT,     NULL,   NULL, "n>0",    NULL,        NULL,"--fil-only", "set length of random seqs for gumbel estimation to <n>", 2},
  { "--gum-cmN",        eslARG_INT,     "1000", NULL, "n>0",    NULL,        NULL,"--fil-only", "number of random sequences for CM gumbel estimation",    2 },
  { "--gum-hmmN",       eslARG_INT,     "5000", NULL, "n>0",    NULL,        NULL,"--fil-only", "number of random sequences for CP9 HMM gumbel estimation",    2 },
  { "--gum-gcfromdb",   eslARG_INFILE,  NULL,   NULL, NULL,     NULL,        NULL,"--fil-only", "use GC content distribution from file <s>",  2},
  { "--gum-beta",       eslARG_REAL,    NULL,   NULL, "x>0.000000000000001", NULL,NULL,   NULL, "set tail loss prob for QDBs in gumbel calc to <x>, [default: no QDBs]", 5 },
  { "--gum-gtail",      eslARG_REAL,    "0.5",  NULL, "0.0<x<0.6",NULL,      NULL,        NULL, "set fraction of right histogram tail to fit to gumbel to <x>", 5 },
  { "--gum-pfile",      eslARG_INFILE,  NULL,   NULL, NULL,     NULL,"--gcfromdb",        NULL, "read partition info for gumbels from file <s>", 2},
  { "--gum-hfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,"--fil-only", "save fitted gumbel histogram(s) to file <s>", 2 },
  { "--gum-sfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,"--fil-only", "save survival plot to file <s>", 2 },
  { "--gum-qqfile",     eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,"--fil-only", "save Q-Q plot for gumbel histogram(s) to file <s>", 2 },
  /* options for HMM filter threshold calculation */
  { "--fil-only",       eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,"--gum-only", "only calculate filter thresholds, don't estimate gumbels", 3},
  { "--fil-N",          eslARG_INT,     "1000", NULL, "10<=n<=100000",NULL,  NULL,"--gum-only", "number of emitted sequences for HMM filter threshold calc",    3 },
  { "--fil-F",          eslARG_REAL,    "0.99", NULL, "0<x<=1", NULL,        NULL,"--gum-only", "required fraction of seqs that survive HMM filter", 3},
  { "--fil-beta",       eslARG_REAL,    "1e-7", NULL, "x>0",    NULL,        NULL,        NULL, "set tail loss prob for QDBs in filter calculation to <x>", 5 },
  { "--fil-noqdb",      eslARG_NONE,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "do not use QDBs for filter calculation", 5 },
  { "--fil-xhmm",       eslARG_REAL,    "2.0",  NULL, "x>=1.1", NULL,        NULL,"--gum-only", "set target time for filtered search as <x> times HMM time", 3},
  { "--fil-hbanded",    eslARG_NONE,    NULL,   NULL, NULL,     NULL,        NULL,"--gum-only", "use HMM banded search for filter calculation", 3},
  { "--fil-tau",        eslARG_REAL,    "1e-7", NULL, "0<x<1",  NULL, "--fbanded",        NULL, "set tail loss prob for --fbanded to <x>", 6 },
  { "--fil-scan2bands", eslARG_NONE,    FALSE,  NULL, NULL,     NULL, "--fbanded",        NULL, "derive HMM bands from scanning Forward/Backward", 6 },
  { "--fil-gemit",      eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,"--gum-only", "when calc'ing filter thresholds, always emit globally from CM",  3},
  { "--fil-hfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,"--gum-only", "save CP9 filter threshold histogram(s) to file <s>", 3},
  { "--fil-rfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,"--gum-only", "save CP9 filter threshold information in R format to file <s>", 3},
/* Other options */
  { "--stall",          eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "arrest after start: for debugging MPI under gdb", 5 },  
  { "--mxsize",         eslARG_REAL,    "256.0",NULL, "x>0.",   NULL,        NULL,        NULL, "set maximum allowable HMM banded DP matrix size to <x> Mb", 9 },
#ifdef HAVE_MPI
  { "--mpi",            eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "run as an MPI parallel program", 5 },  
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
  int              ncm;                /* what number CM we're on */
  int              be_verbose;	       /* print extra info, only can be TRUE if #ifdef HAVE_DEVOPTS */
  int              cmalloc;            /* number of cmstats we have allocated */
  char            *tmpfile;            /* tmp file we're writing to */
  char            *mode;               /* write mode, "w" or "wb"                     */
  int              np;                 /* number of partitions, must be 1 unless --gum-pfile or --fil-only invoked,
					* once set, is never changed, once case we look out for: 
					* if --fil-only invoked and CM file has CMs with diff number of partitions,
					* this is unlikely but possible, in this case we die in initialize_cmstats() */
  int             *pstart;             /* [0..p..np-1], begin points for partitions, end pts are implicit */
  float           *avglen;             /* [0..v..M-1] average hit len for subtree rooted at each state v for current CM */

  /* info for the comlog we'll add to the cmfiles */
  char            *ccom;               /* command line used in this execution of cmcalibrate */
  char            *cdate;              /* date of this execution of cmcalibrate */


  /* the following data is modified for each CM, and in some cases for each Gumbel mode for each CM,
   * it is assumed to be 'current' in many functions.
   */
  float            gum_cm_ncalcs;   /* millions of calcs for each CM scan of 1 residue with NO QDBs, unless --gum-beta <x> invoked, then it's with QDBs with beta=<x> */
  float            fil_cm_ncalcs;   /* millions of calcs for full CM scan of 1 residue, with --fil-beta <x>, if --fil-noqdb, beta = 0. */
  float            cp9_ncalcs;      /* millions of calcs for CP9 HMM scan of 1 residue, updated when model is localfied */
  /* mpi */
  int              do_mpi;
  int              my_rank;
  int              nproc;
  int              do_stall;          /* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE          *cmfp;	      /* open input CM file stream       */
  FILE            *gumhfp;            /* optional output for gumbel histograms */
  FILE            *gumsfp;            /* optional output for gumbel survival plot */
  FILE            *gumqfp;            /* optional output for gumbel QQ file */
  FILE            *filhfp;            /* optional output for filter histograms */
  FILE            *filrfp;            /* optional output for filter info for R */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit Gumbels for E-values and calculate HMM filter thresholds";

static int init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
				   float **ret_cyk_scA, float **ret_ins_scA, float **ret_fwd_scA, int **ret_partA);
static int process_gumbel_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq, 
				   int L, int use_cm, float **ret_scA);
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);

static int  set_partition_gc_freq(struct cfg_s *cfg, int p);
static int  fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, double *ret_mu, double *ret_lambda);
static int  get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, double *dnull, int L, ESL_DSQ **ret_dsq);
static int  get_cmemit_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int *ret_L, int *ret_p, Parsetree_t **ret_tr, ESL_DSQ **ret_dsq);
static void estimate_workunit_time(const ESL_GETOPTS *go, const struct cfg_s *cfg, int nseq, int L, int gum_mode);
static int  read_partition_file(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  update_avg_hit_len(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  switch_global_to_local(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, char *errbuf);
extern int  update_dp_calcs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm);
static void format_time_string(char *buf, double sec, int do_frac);
static int  print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  get_hmm_filter_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_scA, float *fwd_scA, int *partA, int cm_mode, BestFilterInfo_t *bf);
static int  compare_fseq_by_cm_Eval (const void *a_void, const void *b_void);
static int  compare_fseq_by_fwd_Eval(const void *a_void, const void *b_void);

int
main(int argc, char **argv)
{
  int status;
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

  /* Initialize configuration shared across all kinds of masters
   * and workers in this .c file.
   */
  cfg.cmfile  = esl_opt_GetArg(go, 1);
  if (cfg.cmfile == NULL) cm_Fail("Failed to read <cmfile> argument from command line.");
  cfg.r        = NULL; 
  cfg.abc      = NULL; 
  cfg.w_stage  = NULL; 
  cfg.gc_freq  = NULL; 
  cfg.pgc_freq = NULL; 
  cfg.be_verbose = FALSE;
  cfg.cmalloc  = 128;
  cfg.tmpfile  = NULL;
  cfg.mode     = NULL;
  cfg.np        = 1;     /* default number of partitions is 1, changed if --gum-pfile */
  cfg.pstart    = NULL;  /* allocated (by default to size 1) in init_master_cfg() */
  cfg.avglen    = NULL; 

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cfg.cmalloc  = 128;
  ESL_ALLOC(cfg.cmstatsA, sizeof(CMStats_t *) * cfg.cmalloc);
  cfg.ncm      = 0;

  cfg.ccom      = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */
  cfg.cdate     = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */

  cfg.gum_cm_ncalcs = 0;
  cfg.fil_cm_ncalcs = 0;
  cfg.cp9_ncalcs    = 0;

  cfg.do_mpi   = FALSE;
  cfg.my_rank  = 0;
  cfg.nproc    = 0;
  cfg.do_stall = esl_opt_GetBoolean(go, "--stall");

  cfg.cmfp     = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gumhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gumsfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gumqfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.filhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.filrfp   = NULL; /* ALWAYS remains NULL for mpi workers */

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
      cm->flags &= ~CMH_FILTER_STATS; /* forget that CM may have had filter stats, if --gum-only was invoked, it won't anymore */

      if(!(esl_opt_GetBoolean(go, "--fil-only"))) cm->flags |= CMH_GUMBEL_STATS; 
      if(!(esl_opt_GetBoolean(go, "--gum-only"))) cm->flags |= CMH_FILTER_STATS; 
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
    if (cfg.gumhfp   != NULL) { 
      fclose(cfg.gumhfp);
      printf("High score for random seqs histograms saved to file %s.\n", esl_opt_GetString(go, "--gum-hfile"));
    }
    if (cfg.gumsfp   != NULL) { 
      fclose(cfg.gumsfp);
      printf("Survival plot for Gumbels saved to file %s.\n", esl_opt_GetString(go, "--gum-sfile"));
    }
    if (cfg.gumqfp   != NULL) { 
      fclose(cfg.gumqfp);
      printf("Gumbel QQ plots saved to file %s.\n", esl_opt_GetString(go, "--gum-qqfile"));
    }
    if (cfg.filhfp   != NULL) { 
      fclose(cfg.filhfp);
      printf("Filter histograms saved to file %s.\n", esl_opt_GetString(go, "--fil-hfile"));
    }
    if (cfg.filrfp   != NULL) {
      fclose(cfg.filrfp);
      printf("Filter R info saved to file %s.\n", esl_opt_GetString(go, "--fil-rfile"));
    }
    if (cfg.ccom     != NULL) free(cfg.ccom);
    if (cfg.cdate    != NULL) free(cfg.cdate);
  }
  /* clean up */
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
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
 *    cfg->gumhfp      - optional output file
 *    cfg->gumsfp      - optional output file
 *    cfg->gumqfp      - optional output file
 *    cfg->filhfp      - optional output file
 *    cfg->filrfp      - optional output file
 *    cfg->gc_freq     - observed GC freqs (if --gum-gcfromdb invoked)
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
  if (esl_opt_GetString(go, "--gum-hfile") != NULL) 
    {
      if ((cfg->gumhfp = fopen(esl_opt_GetString(go, "--gum-hfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open gumbel histogram save file %s for writing\n", esl_opt_GetString(go, "--gum-hfile"));
    }

  /* optionally, open survival plot */
  if (esl_opt_GetString(go, "--gum-sfile") != NULL) 
    {
      if ((cfg->gumsfp = fopen(esl_opt_GetString(go, "--gum-sfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open survival plot save file %s for writing\n", esl_opt_GetString(go, "--gum-sfile"));
    }

  /* optionally, open gumbel QQ plot file */
  if (esl_opt_GetString(go, "--gum-qqfile") != NULL) 
    {
      if ((cfg->gumqfp = fopen(esl_opt_GetString(go, "--gum-qqfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open gumbel QQ plot save file %s for writing\n", esl_opt_GetString(go, "--gum-qqfile"));
    }

  /* optionally, open filter threshold calc histogram file */
  if (esl_opt_GetString(go, "--fil-hfile") != NULL) {
    if ((cfg->filhfp = fopen(esl_opt_GetString(go, "--fil-hfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --fil-hfile output file %s\n", esl_opt_GetString(go, "--fil-hfile"));
    }

  /* optionally, open filter threshold calc info file */
  if (esl_opt_GetString(go, "--fil-rfile") != NULL) {
    if ((cfg->filrfp = fopen(esl_opt_GetString(go, "--fil-rfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --fil-rfile output file %s\n", esl_opt_GetString(go, "--fil-rfile"));
    }

  /* optionally, get distribution of GC content from an input database (default is use cm->null for GC distro) */
  if(esl_opt_GetString(go, "--gum-gcfromdb") != NULL) {
    ESL_ALPHABET *tmp_abc = NULL;
    tmp_abc = esl_alphabet_Create(eslRNA);
    ESL_SQFILE      *dbfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--gum-gcfromdb"), eslSQFILE_UNKNOWN, NULL, &dbfp);
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

  /* set up the partition data that's used for all CMs */
  if(esl_opt_IsDefault(go, "--gum-pfile")) { /* by default we have 1 partition 0..100 */
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
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  int      p;                     /* partition index */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  long     seed;                  /* RNG seed */
  int      working_on_cm;         /* TRUE when gum_mode is for CM gumbel */
  int      working_on_cp9;        /* TRUE when gum_mode is for CP9 gumbel */
  int      gum_mode;              /* ctr over gumbel modes */

  /* gumbel related vars */
  int      cmN  = esl_opt_GetInteger(go, "--gum-cmN");  /* number of random sequences to search for CM gumbel fitting */
  int      hmmN = esl_opt_GetInteger(go, "--gum-hmmN"); /* number of random sequences to search for CM gumbel fitting */
  int      gumN;                  /* number of sequences to search for gumbel fitting of current gumbel mode, either --gum-cmN or --gum-hmmN */
  int      gumL;                  /* length of sequences to search for gumbel fitting, L==cm->W*2 unless --gum-L enabled, 
				   * in which case L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gum-L") */  
  float   *gum_scA        = NULL; /* best cm or hmm score for each random seq */
  double   tmp_mu, tmp_lambda;    /* temporary mu and lambda used for setting gumbels */

  /* filter threshold related vars */
  int      filN = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences to search for filter threshold calculation */
  int      fthr_mode = 0;         /* CM mode for filter threshold calculation, FTHR_CM_GC, FTHR_CM_GI, FTHR_CM_LC, FTHR_CM_LI */
  float   *fil_cyk_scA = NULL;    /* [0..filN-1] best cm cyk score for each emitted seq */
  float   *fil_ins_scA = NULL;    /* [0..filN-1] best cm insidei score for each emitted seq */
  float   *fil_fwd_scA = NULL;    /* [0..filN-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA   = NULL;    /* [0..filN-1] partition of CM emitted seq */
  int      cm_cyk_mode;           /* CYK    mode CM is in GUM_CM_LC or GUM_CM_GC */
  int      cm_ins_mode;           /* Inside mode CM is in GUM_CM_LI or GUM_CM_GI */

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
      if (esl_opt_GetBoolean(go, "--fil-only") && (! (cm->flags & CMH_GUMBEL_STATS))) cm_Fail("--fil-only invoked, but CM %s (CM number %d) does not have Gumbel stats in CM file\n", cm->name, (cmi+1));

      if((status = initialize_cm(go, cfg, errbuf, cm))      != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = print_cm_info     (go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);

      for(gum_mode = 0; gum_mode < GUM_NMODES; gum_mode++) {

	if(GumModeIsForCM(gum_mode)) { working_on_cm = TRUE;  working_on_cp9 = FALSE; }
	else                         { working_on_cm = FALSE; working_on_cp9 = TRUE;  }

	/* do we need to switch from glocal configuration to local? */
	if(gum_mode > 0 && (! GumModeIsLocal(gum_mode-1)) && GumModeIsLocal(gum_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
	  if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
	}
	/* update search opts for gumbel mode */
	GumModeToSearchOpts(cm, gum_mode);

	/* We want to use the same seqs for Gumbel fittings of all CM modes and HMM modes, 
	 * so we free RNG, then create a new one and reseed it with the initial seed,
	 * The following pairs of modes will have identical sequences used for each member of the pair:
	 * 1. GUM_CP9_GV and GUM_CP9_GF
	 * 2. GUM_CM_GC  and GUM_CM_GI
	 * 3. GUM_CP9_LV and GUM_CP9_LF
	 * 4. GUM_CM_LC  and GUM_CM_LI
	 * Also the first min(--gum--cmN <n>, --gum--hmmN <n>) sequences between 1 and 2, and between 3 and 4,
	 * will also be identical.
	 */
	seed = esl_randomness_GetSeed(cfg->r);
	esl_randomness_Destroy(cfg->r);
	cfg->r = esl_randomness_Create(seed);

	/**************************/
	/* gumbel fitting section */
	/**************************/
	if(! (esl_opt_GetBoolean(go, "--fil-only"))) {
	  /* calculate gumbels for this gum mode */
	  /* determine length of seqs to search for gumbel fitting */
	  if(esl_opt_IsDefault(go, "--gum-L")) gumL = cm->W*2; 
	  else                                 gumL = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gum-L")); /* minimum L we allow is 2 * cm->W, this is enforced silently (!) */
	  ESL_DASSERT1((cfg->np == cfg->cmstatsA[cmi]->np));
	  for (p = 0; p < cfg->np; p++) {
	    esl_stopwatch_Start(cfg->w_stage);
	    if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);

	    gumN = working_on_cm ? cmN : hmmN;
	    if(cfg->be_verbose) estimate_workunit_time(go, cfg, gumN, gumL, gum_mode);
	    ESL_DPRINTF1(("\n\ncalling process_gumbel_workunit to fit gumbel for p: %d GUM mode: %d\n", p, gum_mode));
	    printf("gumbel  %-12s %5d %6d %5s %5s [", DescribeGumMode(gum_mode), gumN, gumL, "-", "-");
	    fflush(stdout);

	    /* do the work, fit the histogram, update gumbel info in cmstats */
	    if((status = process_gumbel_workunit (go, cfg, errbuf, cm, gumN, gumL, working_on_cm, &gum_scA))!= eslOK) cm_Fail(errbuf);
	    if((status = fit_histogram(go, cfg, errbuf, gum_scA, gumN, &tmp_mu, &tmp_lambda))               != eslOK) cm_Fail(errbuf);
	    SetGumbelInfo(cfg->cmstatsA[cmi]->gumAA[gum_mode][p], tmp_mu, tmp_lambda, gumL, gumN);

	    esl_stopwatch_Stop(cfg->w_stage);
	    format_time_string(time_buf, cfg->w_stage->elapsed, 0);
	    printf(" %10s\n", time_buf);
	  } /* end of for loop over partitions */
	} /* end of if ! --fil-only */
	
	/****************************/
	/* filter threshold section */
	/****************************/
	if((! (esl_opt_GetBoolean(go, "--gum-only"))) && (gum_mode == GUM_CM_GI || gum_mode == GUM_CM_LI)) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  esl_stopwatch_Start(cfg->w_stage);
	  fthr_mode = GumModeToFthrMode(gum_mode);
	  /* search emitted sequences to get filter thresholds for HMM and each candidate sub CM root state */
	  ESL_DPRINTF1(("\n\ncalling process_filter_workunit to get HMM filter thresholds for p: %d mode: %d\n", p, gum_mode));
	  if(fil_partA != NULL) free(fil_partA);
	  /* FIX ME! */ printf("filter  %3s  %-8s %5s %6s %5d %5s [", "-", DescribeFthrMode(fthr_mode), "-", "-", filN, "");
	  fflush(stdout);
	  
	  if((status = process_filter_workunit (go, cfg, errbuf, cm, filN, &fil_cyk_scA, &fil_ins_scA, &fil_fwd_scA, &fil_partA)) != eslOK) cm_Fail(errbuf);
	  
	  cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN)         ? GUM_CM_LC : GUM_CM_GC;
	  cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN)         ? GUM_CM_LI : GUM_CM_GI;
	  /* set cutoffs for forward HMM filters, first for CYK, then for Inside */
	  if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_cyk_scA, fil_fwd_scA, fil_partA, cm_cyk_mode, cfg->cmstatsA[cmi]->bfA[cm_cyk_mode])) != eslOK) cm_Fail(errbuf);
	  if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_ins_scA, fil_fwd_scA, fil_partA, cm_ins_mode, cfg->cmstatsA[cmi]->bfA[cm_ins_mode])) != eslOK) cm_Fail(errbuf);
	  if(cfg->be_verbose) DumpBestFilterInfo(cfg->cmstatsA[cmi]->bfA[fthr_mode]);
	  
	  esl_stopwatch_Stop(cfg->w_stage);
	  format_time_string(time_buf, cfg->w_stage->elapsed, 0);
	  printf(" %10s\n", time_buf);
	}
      } /* end of for(gum_mode = 0; gum_mode < NCMMODES; gum_mode++) */
      if(cfg->be_verbose) debug_print_cmstats(cfg->cmstatsA[cmi], (! esl_opt_GetBoolean(go, "--gum-only")));

      if(! (esl_opt_GetBoolean(go, "--fil-only"))) free(gum_scA);
      if(! (esl_opt_GetBoolean(go, "--gum-only"))) { 
	free(fil_cyk_scA);
	free(fil_ins_scA);
	free(fil_fwd_scA);
	free(fil_partA);
      }
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
  int      status;                /* Easel status */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  int      p;                     /* partition index */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  long     seed;                  /* RNG seed */
  int      n, v, i;               /* counters */
  float    update_i;              /* when this i is reached, print a progress update to stdout */
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
  long    *seedlist = NULL;       /* seeds for worker's RNGs, we send these to workers */
  MPI_Status mpistatus;           /* MPI status... */
  int      msg;                   /* holds integer telling workers we've finished current stage */
  int      working_on_cm;         /* TRUE when gum_mode is for CM gumbel */
  int      working_on_cp9;        /* TRUE when gum_mode is for CP9 gumbel */
  int      gum_mode;              /* ctr over gumbel modes */

  /* gumbel related vars */
  int      cmN  = esl_opt_GetInteger(go, "--gum-cmN");  /* number of random sequences to search for CM gumbel fitting */
  int      hmmN = esl_opt_GetInteger(go, "--gum-hmmN"); /* number of random sequences to search for CM gumbel fitting */
  int      gumN;                  /* number of sequences to search for gumbel fitting of current gumbel mode, either --gum-cmN or --gum-hmmN */
  int      gumL;                  /* length of sequences to search for gumbel fitting, L==cm->W*2 unless --gum-L enabled, 
				   * in which case L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gum-L") */  int      gum_mode  = 0;
  int      gum_nseq_per_worker  = 0; /* when calcing gumbels, number of seqs to tell each worker to work on */
  int      gum_nseq_this_round  = 0; /* when calcing gumbels, number of seqs for current round */
  float   *gum_scA        = NULL; /* [0..gumN-1] full list best cm or hmm score for each random seq */
  /* worker's Gumbel scores [0..nseq_per_worker-1], rec'd from workers, copied to full arrays (ex: fil_cyk_scA) */
  float   *wkr_gum_scA = NULL;    /* [0..nseq_per_worker-1] best cm or hmm score for worker's random seqs, rec'd from worker */
  double   tmp_mu, tmp_lambda;    /* temporary mu and lambda used for setting gumbels */

  /* filter threshold related vars */
  int      filN = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences to search for filter threshold calculation */
  int      fthr_mode = 0;         /* CM mode for filter threshold calculation, FTHR_CM_GC, FTHR_CM_GI, FTHR_CM_LC, FTHR_CM_LI */
  int      cm_cyk_mode;           /* CYK    mode CM is in GUM_CM_LC or GUM_CM_GC */
  int      cm_ins_mode;           /* Inside mode CM is in GUM_CM_LI or GUM_CM_GI */
  int      fil_nseq_per_worker  = (filN / (cfg->nproc-1)); /* when calcing filters, number of seqs to tell each worker to work on */
  /* full arrays of CYK, Inside, Fwd scores, [0..filN-1] */
  float   *fil_cyk_scA = NULL;    /* [0..filN-1] best cm cyk score for each emitted seq */
  float   *fil_ins_scA = NULL;    /* [0..filN-1] best cm insidei score for each emitted seq */
  float   *fil_fwd_scA = NULL;    /* [0..filN-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA   = NULL;    /* [0..filN-1] partition of CM emitted seq */
  /* worker's arrays of CYK, Inside, Fwd scores, [0..nseq_per_worker-1], rec'd from workers, copied to full arrays (ex: fil_cyk_scA) */
  float   *wkr_fil_cyk_scA = NULL;/* rec'd from worker: best cm cyk score for each emitted seq */
  float   *wkr_fil_ins_scA = NULL;/* rec'd from worker: best cm insidei score for each emitted seq */
  float   *wkr_fil_fwd_scA = NULL;/* rec'd from worker: best cp9 Forward score for each emitted seq */
  int     *wkr_fil_partA = NULL;  /* rec'd from worker: partition for seq i */


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
  for (wi = 0; wi < cfg->nproc; wi++) {
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
   * case 1: if we've used the --gum-gcfromdb option, we read in a seq file 
   * to fill cfg->gc_freq, and we need to broadcast that info to workers
   *
   * case 2: if we are calculating stats for more than 1 partition, 
   * (--gum-pfile invoked), we need to broadcast that information to 
   * the workers. 
   */
  if(! (esl_opt_IsDefault(go, "--gum-gcfromdb"))) { /* receive gc_freq info from master */
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  if(! (esl_opt_IsDefault(go, "--gum-pfile"))) { /* broadcast partition info to workers */
    ESL_DASSERT1((! (esl_opt_GetBoolean(go, "--fil-only"))));
    MPI_Bcast(&(cfg->np),  1,       MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pstart != NULL));
    MPI_Bcast(cfg->pstart, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
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
      if (esl_opt_GetBoolean(go, "--fil-only") && (! (cm->flags & CMH_GUMBEL_STATS))) cm_Fail("--fil-only invoked, but CM %s (CM number %d) does not have Gumbel stats in CM file\n", cm->name, (cmi+1));

      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm     (go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      if((status = print_cm_info     (go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);

      if(! (esl_opt_GetBoolean(go, "--gum-only"))) { 
	ESL_ALLOC(fil_cyk_scA, sizeof(float) * filN);
	ESL_ALLOC(fil_ins_scA, sizeof(float) * filN);
	ESL_ALLOC(fil_fwd_scA, sizeof(float) * filN);
	ESL_ALLOC(fil_partA,   sizeof(int) *   filN);
      }

      for(gum_mode = 0; gum_mode < GUM_NMODES; gum_mode++) {

	if(GumModeIsForCM(gum_mode)) { working_on_cm = TRUE;  working_on_cp9 = FALSE; }
	else                         { working_on_cm = FALSE; working_on_cp9 = TRUE;  }
	gum_nseq_per_worker = working_on_cm ? (int) (cmN / (cfg->nproc-1)) : (int) (hmmN / (cfg->nproc-1));
	gum_nseq_this_round = working_on_cm ? cmN : hmmN;
	update_i = gum_nseq_this_round / 20.;

	/* do we need to switch from glocal configuration to local? */
	if(gum_mode > 0 && (! GumModeIsLocal(gum_mode-1)) && GumModeIsLocal(gum_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
	  if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK)     cm_Fail(errbuf);
	}
	/* update search opts for gumbel mode */
	GumModeToSearchOpts(cm, gum_mode);

	/**************************/
	/* gumbel fitting section */
	/**************************/
	if(! (esl_opt_GetBoolean(go, "--fil-only"))) {
	  if(esl_opt_IsDefault(go, "--gum-L")) gumL = cm->W*2; 
	  else                                 gumL = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gum-L")); /* minimum L we allow is 2 * cm->W, this is enforced silently (!) */
	  gumN = (working_on_cm) ? cmN : hmmN;
	  ESL_ALLOC(gum_scA, sizeof(float) * gumN);

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
			  if ((status = cmcalibrate_gumbel_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &wkr_gum_scA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
			    ESL_DPRINTF1(("MPI master has unpacked CM gumbel results\n"));
			    ESL_DASSERT1((nseq_just_recv > 0));
			    for(i = 0; i < nseq_just_recv; i++) {
			      gum_scA[nseq_recv+i] = wkr_gum_scA[i];
			      if(nseq_recv+i > update_i) { /* print progress update to stdout */
				printf("=");
				fflush(stdout); 
				update_i += gum_nseq_this_round / 20.; 
			      }
			    }
			    free(wkr_gum_scA);
			    nseq_recv += nseq_just_recv;
			}
		      else /* worker reported an error. Get the errbuf. */
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
	      if((status = fit_histogram(go, cfg, errbuf, gum_scA, gumN, &tmp_mu, &tmp_lambda))       != eslOK) cm_Fail(errbuf);
	      SetGumbelInfo(cfg->cmstatsA[cmi]->gumAA[gum_mode][p], tmp_mu, tmp_lambda, L, hmmN);
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
	} /* end of if ! --fil-only */

	/****************************/
	/* filter threshold section */
	/****************************/
	if((! (esl_opt_GetBoolean(go, "--gum-only"))) && (gum_mode == GUM_CM_GI || gum_mode == GUM_CM_LI)) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  esl_stopwatch_Start(cfg->w_stage);
	  fthr_mode = GumModeToFthrMode(gum_mode);
	  ESL_DPRINTF1(("MPI master: CM: %d fthr mode: %d\n", cfg->ncm, fthr_mode));
	  /* FIX ME! */printf("filter  %3s  %-8s %5s %6s %5d %5s [", "-", DescribeFthrMode(fthr_mode), "-", "-", "",filN); 
	  fflush(stdout);

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
			if ((status = cmcalibrate_filter_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &wkr_fil_cyk_scA, &wkr_fil_ins_scA, &wkr_fil_fwd_scA, &wkr_fil_partA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
			ESL_DPRINTF1(("MPI master has unpacked HMM filter results\n"));
			for(i = 0; i < nseq_just_recv; i++) {
			  fil_cyk_scA[nseq_recv+i] = wkr_fil_cyk_scA[i];
			  fil_ins_scA[nseq_recv+i] = wkr_fil_ins_scA[i];
			  fil_fwd_scA[nseq_recv+i] = wkr_fil_fwd_scA[i];
			  fil_partA[nseq_recv+i]   = wkr_fil_partA[i];
			  ESL_DASSERT1((fil_partA[nseq_recv+i] < cfg->np));
			  /* print progress update to output */
			  if(nseq_recv+i > update_i) {
			    printf("=");
			    fflush(stdout); 
			    update_i += filN / 20.; 
			  }
			}
			free(wkr_fil_cyk_scA);
			free(wkr_fil_ins_scA);
			free(wkr_fil_fwd_scA);
			free(wkr_fil_partA);
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
	    cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN)         ? GUM_CM_LC : GUM_CM_GC;
	    cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN)         ? GUM_CM_LI : GUM_CM_GI;
	    /* set cutoffs for forward HMM filters, first for CYK, then for Inside */
	    if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_cyk_scA, fil_fwd_scA, fil_partA, cm_cyk_mode, cfg->cmstatsA[cmi]->bfA[cm_cyk_mode])) != eslOK) cm_Fail(errbuf);
	    if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_ins_scA, fil_fwd_scA, fil_partA, cm_ins_mode, cfg->cmstatsA[cmi]->bfA[cm_ins_mode])) != eslOK) cm_Fail(errbuf);
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
      if(xstatus == eslOK) if(cfg->be_verbose) { debug_print_cmstats(cfg->cmstatsA[cmi], (! esl_opt_GetBoolean(go, "--gum-only"))); }
      
      if(! (esl_opt_GetBoolean(go, "--fil-only"))) free(gum_scA);
      if(! (esl_opt_GetBoolean(go, "--gum-only"))) { 
	free(fil_cyk_scA);
	free(fil_ins_scA);
	free(fil_fwd_scA);
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
  int      status;                /* Easel status */
  int      xstatus = eslOK;       /* changes from OK on recoverable error */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  char    *wbuf = NULL;	          /* packed send/recv buffer  */
  int      wn   = 0;	          /* allocation size for wbuf */
  int      sz, n;		  /* size of a packed message */
  int      pos;                   /* posn in wbuf */
  long     seed;                  /* seed for RNG, rec'd from master */
  int      working_on_cm;         /* TRUE when gum_mode is for CM gumbel */
  int      working_on_cp9;        /* TRUE when gum_mode is for CP9 gumbel */
  int      nseq;                  /* number of seqs to emit/search for current job */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  MPI_Status mpistatus;           /* MPI status... */

  /* gumbel related vars */
  int      gumL;                  /* length of sequences to search for gumbel fitting, L==cm->W*2 unless --gum-L enabled, 
				   * in which case L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gum-L") */  int      gum_mode  = 0;
  float   *gum_scA        = NULL; /* [0..nseq-1] list of best cm or hmm score for each random seq */

  /* filter threshold related vars */
  float   *fil_cyk_scA = NULL;    /* [0..nseq-1] best cm cyk score for each emitted seq */
  float   *fil_ins_scA = NULL;    /* [0..nseq-1] best cm insidei score for each emitted seq */
  float   *fil_fwd_scA = NULL;    /* [0..nseq-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA   = NULL;    /* [0..nseq-1] partition of CM emitted seq */
  int      in_fil_section_flag = FALSE; /* set to TRUE while we're in the filter threshold calculation
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
   * case 1: if we've used the --gum-gcfromdb option, we read in a seq file to fill
   * cfg->gc_freq, and we need that info here for the worker, so we receive
   * it's broadcast from the master
   * 
   * case 2: if we are calculating stats for more than 1 
   * partition, (--gum-pfile invoked), we need to receive that information 
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
  if(! (esl_opt_IsDefault(go, "--gum-pfile"))) { /* receive partition info from master */
    MPI_Bcast(&(cfg->np),     1, MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pstart == NULL));
    ESL_ALLOC(cfg->pstart, sizeof(int) * cfg->np);
    MPI_Bcast(cfg->pstart, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else { /* no --gum-pfile, set up default partition info */  
    cfg->np     = 1;
    ESL_ALLOC(cfg->pstart, sizeof(int) * cfg->np);
    cfg->pstart[0] = 0;
  }
  
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

	if(GumModeIsForCM(gum_mode)) { working_on_cm = TRUE;  working_on_cp9 = FALSE; }
	else                         { working_on_cm = FALSE; working_on_cp9 = TRUE;  }
	ESL_DPRINTF1(("worker: %d gum_mode: %d nparts: %d\n", cfg->my_rank, gum_mode, cfg->np));

	/* do we need to switch from glocal configuration to local? */
	if(gum_mode > 0 && (! GumModeIsLocal(gum_mode-1)) && GumModeIsLocal(gum_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf)) != eslOK) goto ERROR;
	  if((status = update_avg_hit_len(go, cfg, errbuf, cm)) != eslOK)     goto ERROR;
	}
	/* update search opts for gumbel mode */
	GumModeToSearchOpts(cm, gum_mode);

	/* We want to use the same seqs for Gumbel fittings of all CM modes and HMM modes, 
	 * so we free RNG, then create a new one and reseed it with the initial seed,
	 * The following pairs of modes will have identical sequences used for each member of the pair:
	 * 1. GUM_CP9_GV and GUM_CP9_GF
	 * 2. GUM_CM_GC  and GUM_CM_GI
	 * 3. GUM_CP9_LV and GUM_CP9_LF
	 * 4. GUM_CM_LC  and GUM_CM_LI
	 * Also the first min(--gum--cmN <n>, --gum--hmmN <n>) sequences between 1 and 2, and between 3 and 4,
	 * will also be identical.
	 */
	seed = esl_randomness_GetSeed(cfg->r);
	esl_randomness_Destroy(cfg->r);
	cfg->r = esl_randomness_Create(seed);	seed = esl_randomness_GetSeed(cfg->r);
	
	/**************************/
	/* gumbel fitting section */
	/**************************/
	if(! (esl_opt_GetBoolean(go, "--fil-only"))) {
	  if(esl_opt_IsDefault(go, "--gum-L")) L = cm->W*2; 
	  else                                 L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gum-L")); /* minimum L we allow is 2 * cm->W, this is enforced silently (!) */
	  for (p = 0; p < cfg->np; p++) { /* for each partition */
	    ESL_DPRINTF1(("worker %d gum_mode: %d partition: %d\n", cfg->my_rank, gum_mode, p));

	    if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	  
	    if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	    while(nseq != MPI_FINISHED_GUMBEL) {
	      ESL_DPRINTF1(("worker %d: has received nseq: %d\n", cfg->my_rank, nseq));
	    
	      /* do the work */
	      if((status = process_gumbel_workunit (go, cfg, errbuf, cm, nseq, gumL, working_on_cm, &gum_scA)) != eslOK) goto ERROR;
	      ESL_DPRINTF1(("worker %d: has gathered gumbel results\n", cfg->my_rank));
	      n = 0;
	      if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
		ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	      n += sz;
	      if (cmcalibrate_gumbel_results_MPIPackSize(gum_scA, nseq, MPI_COMM_WORLD, &sz) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_gumbel_results_MPIPackSize() call failed"); 
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
	      if ((status = cmcalibrate_gumbel_results_MPIPack(gum_scA, nseq, wbuf, wn, &pos, MPI_COMM_WORLD)) != eslOK) 
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cm_gumbel_results_MPIPack() call failed.");
	      free(gum_scA);
	      MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	      ESL_DPRINTF1(("worker %d: has sent gumbel results to master in message of %d bytes\n", cfg->my_rank, pos));

	      /* receive next number of sequences, if MPI_FINISHED_GUMBEL, we'll stop */
	      if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	    }
	    ESL_DPRINTF1(("worker %d gum_mode: %d finished partition: %d\n", cfg->my_rank, gum_mode, p));
	  }
	  ESL_DPRINTF1(("worker %d finished all partitions for gum_mode: %d\n", cfg->my_rank, gum_mode));
	} /* end of if ! --fil-only */
	
	/****************************/
	/* filter threshold section */
	/****************************/
	if((! (esl_opt_GetBoolean(go, "--gum-only"))) && (gum_mode == GUM_CM_GI || gum_mode == GUM_CM_LI)) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  in_fil_section_flag = TRUE;
	  fthr_mode = GumModeToFthrMode(gum_mode);
	  ESL_DPRINTF1(("worker %d fthr_mode: %d\n", cfg->my_rank, fthr_mode));

	  if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	  while(nseq != MPI_FINISHED_FILTER) {
	    ESL_DPRINTF1(("worker %d: has received hmm filter nseq: %d\n", cfg->my_rank, nseq));
	    
	    if((status = process_filter_workunit (go, cfg, errbuf, cm, nseq, &fil_cyk_scA, &fil_ins_scA, &fil_fwd_scA, &fil_partA)) != eslOK) cm_Fail(errbuf);
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
	      if (cmcalibrate_filter_results_MPIPack(fil_cyk_scA, fil_ins_scA, fil_fwd_scA, fil_partA, nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPack() call failed"); 
	    free(fil_cyk_scA);
	    free(fil_ins_scA);
	    free(fil_fwd_scA);
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
    if(GumModeIsForCM(gum_mode) && (! (esl_opt_GetBoolean(go, "--gum-only")))) {
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
 *           search them with either (a) the CM or (b) the CM's CP9 HMM.
 *
 *           Thus, this function can be run in 1 of 2 modes, determined by the
 *           <use_cm> input variable. 

 *           If <use_cm> == TRUE, search with the CM with either CYK or Inside 
 *           (as specified by cm->search_opts>. <ret_scA> is filled with the 
 *           best CM score (at root state, v==0) for each sequence. 
 *
 *           If <use_cm> is FALSE, search with the CP9 HMM with either Viterbi or
 *           Forward (as specified by cm->search_opts). <ret_scA> is filled with
 *           the best CP9 score for each sequence.

 *
 * Args:     go           - getopts
 *           cfg          - cmcalibrate's configuration
 *           errbuf       - for writing out error messages
 *           cm           - the CM (already configured as we want it)
 *           nseq         - number of seqs to generate
 *           L            - length of sequences to search, L==cm->W*2 unless --gum-L enabled, 
 *                          in which case
 *                          L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gum-L")
 *           use_cm       - TRUE to search with CM, FALSE to search with CP9 HMM,
 *                          search algorithm determined by cm->search_opts 
 *           ret_scA      - RETURN: [0..nseq-1] best CM or CP9 score for each seq
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_gumbel_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq, int L, int use_cm, float **ret_scA)
{
  int            status;
  float         *scA = NULL;   /* [0..i..nseq-1] best CM or CP9 score for each seq, */
  double        *dnull = NULL; /* double version of cm->null, for generating random seqs */
  int            i;            /* counter */
  ESL_DSQ       *dsq;          /* the digitized sequence to search */
  float          update_i = nseq / 20.; /* when to print progress update to stdout */

  /* determine mode, and enforce mode-specific contract */
  if(ret_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_gumbel_workunit(): ret_scA is NULL.");

  ESL_DPRINTF1(("in process_gumbel_workunit nseq: %d L: %d\n", nseq, L, mode));

  /* allocate scA */
  ESL_ALLOC(scA, sizeof(float) * nseq);

  /* fill dnull, a double version of cm->null, but only if we're going to need it to generate random seqs */
  if(cfg->pgc_freq == NULL) {
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);    
  }

  int do_cyk     = FALSE;
  int do_inside  = FALSE;
  int do_viterbi = FALSE;
  int do_forward = FALSE;
  /* determine algs we'll use and allocate the score arrays we'll pass back */
  if(use_cm) { 
    if(cm->search_opts & CM_SEARCH_INSIDE) do_inside = TRUE;
    else                                   do_cyk    = TRUE;
  }
  else { /* use CP9 HMM */
    if(cm->search_opts & CM_SEARCH_HMMVITERBI) do_viterbi = TRUE;
    if(cm->search_opts & CM_SEARCH_HMMFORWARD) do_forward = TRUE;
    if((do_viterbi + do_forward) > 1) ESL_FAIL(eslEINVAL, errbuf, "process_gumbel_workunit, use_cm == FALSE and cm->search_opts CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags both raised.");
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
    if (do_cyk)    if((status = FastCYKScan    (cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &(scA[i]))) != eslOK) return status;
    if (do_inside) if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &(scA[i]))) != eslOK) return status;
    /* if nec, search with CP9 */
    if (do_viterbi) 
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       NULL,   /* don't want traces back */
			       &(scA[i]))) != eslOK) return status;
    if (do_forward) {
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       &(scA[i]))) != eslOK) return status;
    }

    /*to print seqs to stdout uncomment this block 
    ESL_SQ *tmp;
    tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
    esl_sq_Textize(tmp);
    printf(">seq%d\n%s\n", i, tmp->seq);
    esl_sq_Destroy(tmp);
    */
    free(dsq);
  }
  if(cfg->my_rank == 0) { printf("=]"); }

  if(dnull != NULL) free(dnull);
  *ret_scA = scA;
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
  int safe_windowlen;
  int *dmin, *dmax;

  /* config QDB? NO!, unless --gum-beta enabled */
  if(esl_opt_IsDefault(go, "--gum-beta")) { 
    cm->search_opts |= CM_SEARCH_NOQDB; /* don't use QDB to search */
    cm->beta = DEFAULT_BETA;
  }
  else {
    cm->config_opts |= CM_CONFIG_QDB;   /* configure QDB */
    cm->beta = esl_opt_GetReal(go, "--gum-beta"); 
  }

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
  if(esl_opt_GetBoolean(go, "--fil-gemit")) { 
    cm->flags |= CM_EMIT_NO_LOCAL_BEGINS; 
    cm->flags |= CM_EMIT_NO_LOCAL_ENDS;
  }

  ConfigCM(cm, NULL, NULL);
  
  /* we may still need to determine cm->W */
  safe_windowlen = cm->clen * 2;
  if(esl_opt_IsDefault(go, "--gum-beta")) { /* NO QDBs, we still need to setup cm->W */ 
    if(cm->dmin != NULL || cm->dmax != NULL) 
      cm_Fail("initialize_cm() --gum-beta NOT enabled, but cm->dmin and cm->dmax non-null. This shouldn't happen.");
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

  /* create and initialize scan info for CYK/Inside scanning functions */
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("initialize_cm(), CreateScanMatrixForCM() call failed.");
  
  if((status = update_dp_calcs(go, cfg, errbuf, cm)) != eslOK) return status;
  return eslOK;
}

/* initialize_cmstats()
 * Allocate and initialize a cmstats object in the cfg->cmstatsA array. 
 */
static int
initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int i;
  int p;
  int cmi = cfg->ncm-1;

  ESL_DPRINTF1(("initializing cmstats for %d partitions\n", cfg->np));

  cfg->cmstatsA[cmi] = AllocCMStats(cfg->np);
  
  if(esl_opt_GetBoolean(go, "--fil-only")) { 
    /* set the cfg->np if this is the first CM */
    if(cfg->ncm == 1) { 
      cfg->np = cm->stats->np;
    }
    /* Make sure the rare, rare case that would be a real pain in the ass to implement isn't happening:
     * the case when we use --fil-only with a CM file with > 1 CMs, and at least 2 of those CMs have
     * gumbel stats for different numbers of partitions. If we did deal with this case then 
     * it'd be a bitch to deal with, b/c cfg->np could change (and we'd have to send that
     * info to the workers for each CM in MPI mode).
     */
    ESL_DASSERT1((esl_opt_IsDefault(go, "--gum-gcfromdb"))); /* getopts should enforce this */
    ESL_DASSERT1((esl_opt_IsDefault(go, "--gum-pfile")));  /* getopts should enforce this */
    if(! (cm->flags & CMH_GUMBEL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "--fil-only invoked by CM has no gumbel stats in initialize_cmstats()\n");
    if(cfg->np != cm->stats->np)         ESL_FAIL(eslEINCOMPAT, errbuf, "--fil-only invoked and CM file has CMs with different numbers of partitions. We can't deal. Either split CMs into different files, or recalibrate them fully (without --fil-only)");    
    /* with --fil-only, we're only calc'ing filter thresholds, so we copy the Gumbel stats from cm->stats. */
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

  /* if we get here --fil-only was not invoked */
  ESL_DASSERT1((cfg->pstart[0] == 0));
  for(p = 0; p < cfg->np;     p++) cfg->cmstatsA[cmi]->ps[p] = cfg->pstart[p];
  for(p = 0; p < (cfg->np-1); p++) cfg->cmstatsA[cmi]->pe[p] = cfg->pstart[p+1]-1;
  cfg->cmstatsA[cmi]->pe[(cfg->np-1)] = GC_SEGMENTS-1; /* this is 100 */
  
  for(p = 0; p < cfg->np; p++)
    for(i = cfg->cmstatsA[cmi]->ps[p]; i <= cfg->cmstatsA[cmi]->pe[p]; i++)
      cfg->cmstatsA[cmi]->gc2p[i] = p; 
  return eslOK;
}

/* Function: set_partition_gc_freq()
 * Date:     EPN, Mon Sep 10 08:00:27 2007
 *
 * Purpose:  Set up the GC freq to sample from for the current partition. 
 *           Only used if --gum-gcfromdb used to read in dbseq from which to derive
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
  double mufix;

  ESL_HISTOGRAM *h = NULL;       /* histogram of scores */


  /* Initialize histogram; these numbers are guesses */
  if((h = esl_histogram_CreateFull(-100., 100., .25)) == NULL) return eslEMEM;    

  /* fill histogram */
  for(i = 0; i < nscores; i++)
    if((status = esl_histogram_Add(h, scores[i])) != eslOK) ESL_FAIL(status, errbuf, "fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);

  /* fit scores to a gumbel */
  tailfit = esl_opt_GetReal(go, "--gum-gtail");
  esl_histogram_GetTailByMass(h, tailfit, &xv, &n, &z); /* fit to right 'tailfit' fraction, 0.5 by default */
  esl_gumbel_FitCensored(xv, n, z, xv[0], &mu, &lambda);
  esl_gumbel_FitCensoredLoc(xv, n, z, xv[0], 0.693147, &mufix);

  /* print to output files if nec */
  if(cfg->gumhfp != NULL)
    esl_histogram_Plot(cfg->gumhfp, h);
  if(cfg->gumqfp != NULL) {
      double  params[2];  
      params[0] = mu;
      params[1] = lambda;
      esl_histogram_PlotQQ(cfg->gumqfp, h, &esl_exp_generic_invcdf, params);
  }

  if (cfg->gumsfp != NULL) {
    esl_histogram_PlotSurvival(cfg->gumsfp, h);
    esl_gumbel_Plot(cfg->gumsfp, mu,    lambda,   esl_gumbel_surv, h->xmin - 5., h->xmax + 5., 0.1);
    esl_gumbel_Plot(cfg->gumsfp, mufix, 0.693147, esl_gumbel_surv, h->xmin - 5., h->xmax + 5., 0.1);
  }

  esl_histogram_Destroy(h);

  *ret_mu     = mu;
  *ret_lambda = lambda;
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

  if(! esl_opt_IsDefault(go, "--gum-L")) L = ESL_MAX(L, esl_opt_GetInteger(go, "--gum-L")); /* minimum L we allow is 2 * cm->W (L is sent into this func as 2 * cm->W), this is enforced silently (!) */

  switch(gum_mode) { 
  case GUM_CM_LC: 
  case GUM_CM_GC: 
    seconds = cfg->gum_cm_ncalcs * (float) L * (float) nseq / cyk_megacalcs_per_sec;
    break;
  case GUM_CM_LI:
  case GUM_CM_GI:
    seconds = cfg->gum_cm_ncalcs * (float) L * (float) nseq / ins_megacalcs_per_sec;
    break;
  case GUM_CP9_LV: 
  case GUM_CP9_GV: 
    seconds = cfg->cp9_ncalcs * (float) L * (float) nseq / vit_megacalcs_per_sec;
    break;
  case GUM_CP9_LF: 
  case GUM_CP9_GF: 
    seconds = cfg->cp9_ncalcs * (float) L * (float) nseq / fwd_megacalcs_per_sec;
    break;
  }
  printf("Estimated time for this workunit: %10.2f seconds\n", seconds);

  return;
}

/* Function: read_partition_file
 * Date:     EPN, Fri Dec  7 08:38:41 2007
 * 
 * Called when --gum-pfile is invoked. 
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
  if(esl_opt_IsDefault(go, "--gum-pfile")) ESL_FAIL(eslEINVAL, errbuf, "read_partition_file, but --gum-pfile not invoked!\n");

  if (esl_fileparser_Open(esl_opt_GetString(go, "--gum-pfile"), &efp) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "failed to open %s in read_mask_file\n", esl_opt_GetString(go, "--gum-pfile"));
  esl_fileparser_SetCommentChar(efp, '#');
  
  ESL_ALLOC(begin, sizeof(int) * GC_SEGMENTS);
  begin[0] = 0;

  while((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslEOF) {
    begin[nread] = atoi(tok);
    if(nread == 0) {
      if(atoi(tok) != 0) ESL_FAIL(eslEINVAL, errbuf, "first partition begin must be 0 in %s\n", esl_opt_GetString(go, "--gum-pfile"));
    }
    else if (begin[nread] != (end+1)) {
      if(atoi(tok) != 0) ESL_FAIL(eslEINVAL, errbuf, "partition %d begin point (%d) is not exactly 1 more than prev partition end pt %d in %s\n", (nread+1), begin[nread], end, esl_opt_GetString(go, "--gum-pfile"));
    }      
    if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "no end point for each partition %d's begin (%d) in partition file %s\n", (nread+1), begin[nread], esl_opt_GetString(go, "--gum-pfile"));
    end = atoi(tok);
    if(end < begin[nread]) ESL_FAIL(eslEINVAL, errbuf, "partition %d end point (%d) < begin point (%d) in %s\n", (nread+1), end, begin[nread], esl_opt_GetString(go, "--gum-pfile"));
    nread++;
    if(nread > MAX_PARTITIONS) ESL_FAIL(eslEINVAL, errbuf, "partition file %s has at least %d partitions, but max num partitions is %d\n", esl_opt_GetString(go, "--gum-pfile"), nread, MAX_PARTITIONS);
  }
  if(nread == 0) ESL_FAIL(eslEINVAL, errbuf, "failed to read a single token from %s\n", esl_opt_GetString(go, "--gum-pfile"));
  if(end != 100) ESL_FAIL(eslEINVAL, errbuf, "final partitions end point must be 100, but it's %d in %s\n", end, esl_opt_GetString(go, "--gum-pfile"));

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
 * Args:      go     - get opts
 *            cfg    - cmcalibrate's cfg
 *            cm     - the model
 *            errbuf - for printing errors
 *
 * Returns:   eslOK on succes, othewise some other easel status code and
 *            errbuf is filled with error message.
 */
int 
switch_global_to_local(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  int status;

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

  /* update cfg->fil_cm_ncalcs, cfg->gum_cm_ncalcs and cfg->cp9_ncalcs */
  if((status = update_dp_calcs(go, cfg, errbuf, cm)) != eslOK) return status;

  return eslOK;
}

/* Function: update_dp_calcs()
 * Incept:   EPN, Tue Jan 15 15:40:40 2008
 * 
 * Purpose:  Update cfg->gum_ncalcs, cfg->fil_ncalcs and cfg->cp9_ncalcs
 *           based on configuration of CM.
 *
 * Args:      go     - get opts
 *            cfg    - cmcalibrate's cfg
 *            errbuf - for printing errors
 *            cm     - the model
 *
 * Returns:   eslOK on succes, othewise some other easel status code and
 *            errbuf is filled with error message.
 */
int 
update_dp_calcs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int  status;
  int *dmin; /* potentially used to calculate temporary d bands */
  int *dmax; /* potentially used to calculate temporary d bands */
  int  safe_windowlen; /* potentially used to calculate temporary d bands */
  int  fil_W;          /* cm->W assumed to be used in cmsearch where filter thresholds will be used*/

  /* Count number of DP calcs, these counts are used to determine target, minimum and maximum survival
   * fractions for HMM filter thresholds in get_hmm_filter_cutoffs() (see that code for details).
   * Our predicted running times for searches in get_hmm_filter_cutoffs() are calculated based
   * on cfg->fil_cm_calcs and cfg->fil_cp9_calcs. cfg->fil_cm_calcs is calculated assuming we'll
   * use a QDB filter with beta=1E-7 (by default), unless --fil-beta or --fil-noqdb is enabled.
   *
   * We want to set for global mode:
   * cfg->gum_cm_ncalcs;  millions of calcs for each CM scan of 1 residue with NO QDBs, unless --gum-beta <x> invoked, then it's with QDBs with beta=<x>
   * cfg->fil_cm_ncalcs:  millions of calcs for full CM scan of 1 residue, with --fil-beta <x>, if --fil-noqdb, beta = 0. 
   * cfg->cp9_ncalcs:     millions of calcs for CP9 HMM scan of 1 residue
   *
   * This function is called twice. Once for globa mode and then for local mode when models get localized.
   */

  /* get cfg->gum_cm_ncalcs */
  if((status = cm_CountSearchDPCalcs(cm, errbuf, 10.*cm->W, cm->dmin, cm->dmax, cm->W, FALSE, NULL, &(cfg->gum_cm_ncalcs))) != eslOK) return status;

  /* get cfg->fil_cm_ncalcs */
  if(! (esl_opt_GetBoolean(go, "--fil-noqdb"))) { /* Assume QDBs will be on in cmsearch, determine those QDBs */
    while(!(BandCalculationEngine(cm, safe_windowlen, esl_opt_GetReal(go, "--fil-beta"), 0, &(dmin), &(dmax), NULL, NULL))) {
      free(dmin);
      free(dmax);
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000))
	cm_Fail("initialize_cm(), safe_windowlen big: %d\n", safe_windowlen);
    }
    fil_W = dmax[0];
    free(dmin);
    free(dmax);
    if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*fil_W, dmin, dmax, fil_W, TRUE,  NULL, &(cfg->fil_cm_ncalcs))) != eslOK) return status;
  }
  else { /* assume QDBs will be off in cmsearch, get counts */
    if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*fil_W, NULL, NULL, fil_W, TRUE,  NULL, &(cfg->fil_cm_ncalcs))) != eslOK) return status;
  }
    
  /* get cfg->cp9_ncalcs, used to determine efficiency of CP9 filters, at first it's global mode, then
   * when switch_global_to_local() is called, cfg->full_cp9_ncalcs is updated to ncalcs in local mode */
  int cp9_ntrans = NHMMSTATETYPES * NHMMSTATETYPES; /* 3*3 = 9 transitions in global mode */
  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) cp9_ntrans++; 
  if(cm->cp9->flags & CPLAN9_LOCAL_END)   cp9_ntrans++; 
  if(cm->cp9->flags & CPLAN9_EL)          cp9_ntrans++; 
  cfg->cp9_ncalcs = (cp9_ntrans * cm->cp9->M) / 1000000.; /* convert to millions of calcs per residue */

  return eslOK;
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
			   * cm->ccom2 and cm->cdate2. Can only be 2 if --fil-only enabled, otherwise
			   * we're creating new gumbel stats and/or filter thresholds, so any previous info 
			   * in the CM file from previous cmcalibrate calls will be deleted/overwritten. 
			   */ 
  
  if(ccom  == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), ccom  is non-NULL.");
  if(cdate == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), cdate is non-NULL.");

  which_comlog = (esl_opt_GetBoolean(go, "--fil-only")) ? 2 : 1; /* only way to write to comlog->ccom2, comlog->cdate2 is if --fil-only was enabled */
  
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
    if(cm->comlog->ccom1  == NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), --fil-only enabled, but cm->comlog->ccom1 is NULL.");
    if(cm->comlog->ccom1  == NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), --fil-only enabled, but cm->comlog->cdate1 is NULL.");

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
  if(esl_opt_IsDefault(go, "--gum-beta"))  printf("%-8s %s %g\n", "gumbel beta:", "0.0 (NO qdbs), for W calc: ", cm->beta); 
  else                                     printf("%-8s %g\n",    "gumbel beta:", cm->beta); 
  if(esl_opt_IsDefault(go, "--fil-noqdb")) printf("%-8s %s %g\n", "filter beta:", "0.0 (NO qdbs), for W calc: ", cm->beta); 
  else                                     printf("%-8s %g\n",    "filter beta:", esl_opt_GetReal(go, "--fil-beta")); 
  printf("%-8s %s\n", "command:", cfg->ccom);
  printf("%-8s %s\n", "date:",    cfg->cdate);
  printf("%-8s %d\n", "nproc:",   (cfg->nproc == 0) ? 1 : cfg->nproc);
  printf("\n");
  printf("%6s  %3s  %3s  %3s %5s %6s %5s %5s    percent complete         %10s\n", "",       "",     "",    "",     "",        "",     "",         "",  "");
  printf("%6s  %3s  %3s  %3s %5s %6s %5s %5s [5.......50.......100] %10s\n", "stage",  "mod",  "cfg", "alg",  "gumN",    "len",  "filN",          "",  "elapsed");
  printf("%6s  %3s  %3s  %3s %5s %6s %5s %5s %22s %10s\n", "------", "---", "---", "---", "-----", "------", "-----", "-----", "----------------------", "----------");
  return eslOK;
}


/* Function: process_filter_workunit()
 * Date:     EPN, Fri Jan 11 11:34:46 2008
 *
 * Purpose:  A filter work unit consists of a CM, an int specifying a 
 *           number of sequences <nseq>. The job is to generate <nseq> sequences 
 *           from the CM and search them first with the CM, both CYK and Inside
 *           and then with the HMM (Forward).
 *           Scores will eventually be used for calc'ing HMM filter thresolds.
 *
 * Args:     go             - getopts
 *           cfg            - cmcalibrate's configuration
 *           errbuf         - for writing out error messages
 *           cm             - the CM (already configured as we want it)
 *           nseq           - number of seqs to generate
 *           ret_cyk_scA    - RETURN: [0..nseq-1] best CM CYK score for each seq
 *           ret_ins_scA    - RETURN: [0..nseq-1] best CM Inside score for each seq
 *           ret_fwd_scA    - RETURN: [0..nseq-1] best CP9 Forward score for each seq
 *           ret_partA      - RETURN: [0..nseq-1] partition of each seq 
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			float **ret_cyk_scA, float **ret_ins_scA, float **ret_fwd_scA, int **ret_partA)
{
  int            status;
  float         *cyk_scA  = NULL;  /* [0..i..nseq-1] best CM CYK score for each state, each seq */
  float         *ins_scA  = NULL;  /* [0..i..nseq-1] best CM Inside score for each state, each seq */
  float         *fwd_scA = NULL;   /* [0..i..nseq-1] best CP9 Viterbi score for each seq */
  int           *partA  = NULL;  /* [0..i..nseq-1] partitions of each seq */
  int            p;              /* what partition we're in */
  int            i;
  int            L;
  ESL_DSQ       *dsq;
  Parsetree_t   *tr;
  int            inside_flag_raised = FALSE;
  int            hbanded_flag_raised = FALSE;
  int            scanbands_flag_raised = FALSE;
  float          update_i = nseq / 20.;

  if(ret_cyk_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_cyk_scA != NULL.");
  if(ret_ins_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_ins_scA != NULL.");
  if(ret_fwd_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_fwd_scA != NULL.");
  if(ret_partA   == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_partA != NULL.");

  ESL_DPRINTF1(("in process_filter_workunit nseq: %d mode: %d\n", nseq, mode));

  /* determine algs we'll use and allocate the score arrays we'll pass back */
  ESL_ALLOC(partA, sizeof(int) * nseq); /* will hold partitions */

  ESL_ALLOC(cyk_scA, sizeof(float) * nseq); /* will hold CM CYK scores */
  ESL_ALLOC(ins_scA, sizeof(float) * nseq); /* will hold CM Inside scores */
  ESL_ALLOC(fwd_scA, sizeof(float) * nseq);  /* will hold HMM Forward scores */

  inside_flag_raised = (cm->search_opts & CM_SEARCH_INSIDE) ? TRUE : FALSE;
  hbanded_flag_raised = (cm->search_opts & CM_SEARCH_HBANDED) ? TRUE : FALSE;
  scanbands_flag_raised = (cm->search_opts & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;

  /* generate dsqs one at a time and collect optimal CM CYK/Inside scores and/or best CP9 Forward score */
  for(i = 0; i < nseq; i++) {
    if(cfg->my_rank == 0 && i > update_i) { /* print status update to stdout */
      printf("=");
      fflush(stdout); 
      update_i += nseq / 20.; 
    }
    if((status = get_cmemit_dsq(cfg, errbuf, cm, &L, &p, &tr, &dsq)) != eslOK) return status;
    partA[i] = p;
    /*to print seqs to stdout uncomment this block  
    ESL_SQ *tmp;
    tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
    esl_sq_Textize(tmp);
    printf(">seq%d\n%s\n", i, tmp->seq);
    esl_sq_Destroy(tmp);
    */

    /* search dsq thrice, cyk, inside, fwd */
    /* Note: in mode 2, with FastCYKScan, we use cfg->hsi->smx scan matrix, which may have qdbs calc'ed differently than cm->smx */

    if(! (esl_opt_GetBoolean(go, "--fil-hbanded"))) { 
      cm->search_opts &= ~CM_SEARCH_INSIDE;
      if((status = FastCYKScan    (cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &(cyk_scA[i]))) != eslOK) return status; 
      
      cm->search_opts |= CM_SEARCH_INSIDE; 
      if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &(ins_scA[i]))) != eslOK) return status; 
    }
    else { /* search with HMM bands */
      cm->search_opts &= ~CM_SEARCH_INSIDE;
      cm->search_opts |= CM_SEARCH_HBANDED;
      cm->tau = esl_opt_GetReal(go, "--fil-tau");
      if(esl_opt_GetBoolean(go, "--fil-scan2bands")) cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
      if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, TRUE, 0)) != eslOK) return status; 
      if((status = FastCYKScanHB  (cm, errbuf, dsq, 1, L, 0., NULL, cm->hbmx, 256., &(cyk_scA[i]))) != eslOK) return status; 

      cm->search_opts |= CM_SEARCH_INSIDE; 
      if((status = FastFInsideScanHB(cm, errbuf, dsq, 1, L, 0., NULL, cm->hbmx, 256., &(ins_scA[i]))) != eslOK) return status; 
    }
    if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			     TRUE,   /* yes, we are scanning */
			     FALSE,  /* no, we are not aligning */
			     FALSE,  /* don't be memory efficient */
			     NULL,   /* don't want best score at each posn back */
			     NULL,   /* don't want the max scoring posn back */
			     &(fwd_scA[i]))) != eslOK) return status;
    free(dsq);
  }
  if(cfg->my_rank == 0) { printf("=]"); }
  /* contract enforced these are all non-NULL */
  *ret_cyk_scA = cyk_scA;
  *ret_ins_scA = ins_scA;
  *ret_fwd_scA = fwd_scA;
  *ret_partA = partA;

  if(inside_flag_raised) cm->search_opts |= CM_SEARCH_INSIDE;
  if(! hbanded_flag_raised) cm->search_opts &= ~CM_SEARCH_HBANDED;
  if(! scanbands_flag_raised) cm->search_opts &= ~CM_SEARCH_HMMSCANBANDS;

  return eslOK;

 ERROR:
  return status;
}

typedef struct fseq_Eval_s {
  float i;
  float cm_E;
  float fwd_E;
} fseq_Eval_t;

/*
 * Function: compare_fseq_by_{cm,fwd}_Eval()
 * Date:     EPN, Fri Jan 11 12:30:59 2008
 * Purpose:  Compares two fseq_Eval_t's based on CM or HMM Forward E-value
 *           and returns -1 if first is higher E-value than second, 0 if equal, 
 *           1 if first E-value is lower.  This results in sorting by E-value, 
 *           highest first.
 */
int compare_fseq_by_cm_Eval(const void *a_void, const void *b_void) {
  fseq_Eval_t *a, *b;
  a = (fseq_Eval_t *) a_void;
  b = (fseq_Eval_t *) b_void;
  if      (a->cm_E > b->cm_E) return -1;
  else if (a->cm_E < b->cm_E) return  1;
  else                        return  0;
}

int compare_fseq_by_fwd_Eval(const void *a_void, const void *b_void) {
  fseq_Eval_t *a, *b;
  a = (fseq_Eval_t *) a_void;
  b = (fseq_Eval_t *) b_void;
  if      (a->fwd_E > b->fwd_E) return -1;
  else if (a->fwd_E < b->fwd_E) return  1;
  else                          return  0;
}

/* Function: predict_hmm_filter_speedup_ONCE()
 * Date:     EPN, Fri Jan 11 12:01:07 2008
 *
 * Purpose:  Given a CM and scores for a CM and HMM Forward scan of
 *           filN target seqs predict the HMM filter threshold for
 *           each possible CM threshold. The CM scores are either CYK
 *           or Inside; so this function is called four times, Once
 *           each for glocal CYK, glocal Inside, local CYK and glocal
 *           Inside.
 *           
 *           For the CM scores, the possible CM thresholds are the
 *           E-values for the first 75% (worst scoring 75%) observed
 *           CYK/Inside scores in a ranked list of such E-values,
 *           stored in sorted order in by_cmA[], E-value from
 *           i=0..filN-1. The first 75% are the elements i=0..imax,
 *           with imax = 0.75 * filN.
 *
 *           The HMM threshold fwd_E_cut[i] for each i=0..imax is the
 *           HMM Forward E-value that recognizes F fraction of the
 *           (filN-i+1) sequences that have a CM E-value better than
 *           by_cmA[i].  Note, that when i == imax, 25% of of the
 *           sequences (0.25 * filN) have an E-value of by_cmA[i].cm_E
 *           or lower, this means that the HMM filter for i == imax
 *           will recognize (F * filN * 0.25) sequences, so filN must
 *           be appropriately large so 0.25 * filN is a reasonable
 *           sample size.
 *
 *           To achieve this, we maintain two separately sorted lists,
 *           one by CM E-values (by_cmA) and one by Forward E-values
 *           (by_fwdA), with the cmi2fwdi array providing a map
 *           between the two lists. These two lists contain the same
 *           data (the CM and Forward scores for the filN CM sampled
 *           sequences), but in a different order. If cmi2fwdi[j] ==
 *           k, this means element with rank j in by_cmA has rank k
 *           in by_fwdA, specifically by_cmA[j].i == by_fwdA[k].i.
 *
 *           With these two sorted lists we can step through the CM
 *           list and at each different threshold point (i=0..imax),
 *           update the Forward list by setting specific elements
 *           (specifically those elements whose CM score is greater
 *           (worse) than by_cmA[i].cm_E) as 'no longer used'.  Thus
 *           at each threshold point i, we have the same curN ==
 *           filN-i+1 elements being 'used' in the by_cmA and by_fwdA
 *           lists. Concurrently, we keep track of what HMM threshold
 *           is necessary to recognize F fraction of those curN hits,
 *           specifically this is the Fidx'th ranking hit OF USED
 *           ELEMENTS in the by_fwdA array. This step-through loop is
 *           implemented in the section of code below marked "main
 *           loop"
 *           
 *           It is somewhat tricky to follow the implementation of
 *           this below, at least the way I've done it, and I've found
 *           it equally tricky to comment it in a clear way, but I
 *           tried to make it clear (without spending an unreasonable
 *           amount of time on it).
 *            
 * Args:     go -      command line options
 *           cfg -     cmcalibrate's cfg object
 *           errbuf -  for printing error messages
 *           cm -      the model
 *           cm_scA -  [0..i..filN-1] best CM score (CYK or Inside) in sequence i 
 *           fwd_scA - [0..i..filN-1] best Foward score in sequence i 
 *           partA   - [0..i..filN-1] partition of sequence i 
 *           cm_mode - CM mode that explain the configuration the CM was in when scores in cm_scA were collected
 *                     either GUM_CM_LC (local CYK), GUM_CM_LI (local Inside), GUM_CM_GC (glocal CYK), or GUM_CM_GI (glocal Inside)
 *           bf      - BestFilterInfo_t object, we'll update this to hold info on Forward filter cutoffs
 *
 * Returns:  Updates BestFilterInfo_t object <bf> to hold HMM filter info
 *           eslOK on success;
 *           Other easel status code on an error with errbuf filled with error message.
 */

int
get_hmm_filter_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_scA, float *fwd_scA, int *partA, int cm_mode, BestFilterInfo_t *bf)
{
  int    status;
  float  fil_calcs;               /* number of million dp calcs predicted for the HMM filter scan */
  float  nonfil_calcs;            /* number of million dp calcs predicted for a full CM scan */
  int    evalue_L;                /* database length used for calcing E-values in CP9 gumbels from cfg->cmstats */
  int    i, p;                    /* counters */
  int    cmi   = cfg->ncm-1;      /* CM index we're on */
  float  F     = esl_opt_GetReal   (go, "--fil-F");    /* fraction of CM seqs we require filter to let pass */
  int    filN  = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences we emitted from CM for filter test */
  float  xhmm  = esl_opt_GetReal   (go, "--fil-xhmm"); /* number of sequences we emitted from CM for filter test */
  int    Fidx;                    /* index in sorted scores that threshold will be set at (1-F) * N */
  float  surv_res_per_hit;        /* expected number of residues to survive filter from DB for each hit 2*W-avglen[0], twice W minus the average lenght of a hit from QDB calc */
  float  Smax;                    /* maximally useful survival fraction, set to 0.5 */
  float  Starg;                   /* our target survival fraction, if we achieve this fraction the search should take xhmm times longer than a HMM only search */
  float  Smin;                    /* minimally useful survival fraction, any less than this and our filter would be (predicted to be) doing more than 10X the work of the CM */
  float  fwd_Emax;                /* fwd E-value that gives Smax survival fraction */
  float  fwd_Etarg;               /* fwd E-value that gives Starg survival fraction */
  float  fwd_Emin;                /* fwd E-value that gives Smin survival fraction */
  double dbsize;                  /* database size the E-values refer to, just a scaling factor */
  int   *cmi2fwdi;                /* [0..i..filN] map between by_cmA and by_fwdA elements, cmi2fwdi[j] == k ==> by_cmA[j].i == by_fwdA[k].i */
  int    j;                       /* counter */
  
  fseq_Eval_t *by_cmA;            /* [0..i..filN] list of fseq_Eval_t for all filN observed seqs, sorted by decreasing CM      E-value */
  fseq_Eval_t *by_fwdA;           /* [0..i..filN] list of fseq_Eval_t for all filN observed seqs, sorted by decreasing Forward E-value */
  float        cm_E, fwd_E;       /* a CM E-value and Forward E-value */
  int          hmm_fwd_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? GUM_CP9_LF : GUM_CP9_GF;
  
  /* contract checks */
  if(! (cfg->cmstatsA[cmi]->gumAA[cm_mode][0]->is_valid))      ESL_FAIL(eslEINCOMPAT, errbuf, "get_hmm_filter_cutoffs(), gumbel stats for CM mode: %d are not valid.\n", cm_mode);
  if(! (cfg->cmstatsA[cmi]->gumAA[hmm_fwd_mode][0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "get_hmm_filter_cutoffs(), gumbel stats for HMM fwd mode: %d are not valid.\n", hmm_fwd_mode);

  evalue_L = cfg->cmstatsA[cmi]->gumAA[hmm_fwd_mode][0]->L;
  if(evalue_L != cfg->cmstatsA[cmi]->gumAA[cm_mode][0]->L) ESL_FAIL(eslEINCOMPAT, errbuf,     "get_hmm_filter_cutoffs(), E value L used to HMM forward different than for CM mode: %d. This shouldn't happen.", cm_mode);

  /* contract checks specific to case when there is more than 1 partition */
  if(cfg->cmstatsA[cfg->ncm-1]->np != 1) { 
    for(p = 1; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) {
      if(evalue_L != cfg->cmstatsA[cmi]->gumAA[hmm_fwd_mode][p]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "get_hmm_filter_cutoffs(), partition %d db length (%d) for Forward gumbel stats differ than from partition 1 Forward db length (%d).\n", p, cfg->cmstatsA[cmi]->gumAA[hmm_fwd_mode][p]->L, evalue_L);
    }
  }

  /* TEMPORARY */ 
#define ESL_DPRINTF1(x)  printf x

  /* Determine our target, min and max E-value cutoffs for the HMM forward scan
   * these are independent of the observed scores. Each of these is calculated
   * based on how it affects <total_calcs>. <total_calcs> is the sum of the
   * number of DP calculations required for the HMM filter <fil_calcs>, plus the predicted
   * number of DP calculations required for the CM to search the survival fraction
   * <nonfil_calcs> * survival fract. So if <total_calcs> is <xhmm> * <fil_calcs>, the 
   * required timewe predict to do a HMM filter plus CM scan of survivors is <xhmm> times the
   * time it would take to do an HMM only scan. 
   * 
   * fwd_Etarg:  the target  forward E-value cutoff. If used <total_calcs> = <xhmm> * <fil_calcs>
   *             <xhmm> is obtained from the --fil-xhmm, by default it is 2.0.
   * fwd_Emin:   the minimal forward E-value cutoff, anything less than this is overkill. 
   *             If used <total_calcs> = 1.1 * fil_calcs.
   * fwd_Emax:   the maximum forward E-value cutoff we'll accept as useful. If used 
   *             <total_calcs> = 0.5 * <nonfil_calcs> (we expect HMM filter to give us 
   *                                                   2X speedup versus full CM)
   *
   * forward filter   predicted 
   * E-value cutoff   survival fraction
   * --------------   -----------------
   * fwd_Emax         Smax  = 0.5 -  (fil_calcs / nonfil_calcs)
   * fwd_Etarg        Starg = xhmm * (fil_calcs / nonfil_calcs)
   * fwd_Emin         Smin  = 0.1  * (fil_calcs / nonfil_calcs)
   */
  dbsize = FTHR_DBSIZE_MB * 1000000.;       /* FTHR_DBSIZE_MB is 1, meaning 1 Mb */
  fil_calcs        = cfg->cp9_ncalcs;       /* fil_calcs is millions of DP calcs for HMM Forward scan of 1 residue */
  fil_calcs       *= dbsize;                /* now fil_calcs is millions of DP calcs for HMM Forward scan of length 1 Mb */
  nonfil_calcs     = cfg->fil_cm_ncalcs;    /* total number of millions of DP calculations for full CM scan of 1 residue */
  nonfil_calcs    *= dbsize;                /* now nonfil-calcs corresponds to dbsize (1 Mb)*/
  surv_res_per_hit = (2. * cm->W - (cfg->avglen[0]));   /* avg length of surviving fraction of db from a single hit (cfg->avglen[0] is avg subseq len in subtree rooted at v==0, from QDB calculation) */
  Smin             = 0.1 *  (fil_calcs / nonfil_calcs); /* if survival fraction == Smin, total number of DP calcs (filter + survivors) == (1.1) * fil_calcs, so a HMM + CM search will take only 10% longer than HMM only */
  Starg            = xhmm * (fil_calcs / nonfil_calcs); /* if survival fraction == Starg, the total number DP calcs (filter + survivors) == (1 + xhmm) * fil_calcs,
							 * so if xhmm = 1 (which is default), our target is for the HMM filter + CM search of survivors (with QDB possibly) 
							 * to take 1+1 = 2 times as long as only the HMM filter would take. */ 
  Smax             = 0.5 -  (fil_calcs / nonfil_calcs); /* if survival fraction == Smax, total number of DP calcs (filter+survivors) == 1/2 * nonfil_calcs, so our predicted speedup
							 * by using the filter is 2 fold */
  if(Smin > Smax) { /* rare case, happens only if fil_calcs >= 5/11 * nonfil_calcs, in this case we don't filter */
    cm_Fail("need to implement special case when Smin > Smax\n possible strategy: write 0 points, instead of 50?");
  }
  if(Starg < Smin) { /* another rare case min value for --fil-xhmm is 0.1, so if Starg < Smin it's just because of precision issues, nevertheless we have to deal */
    Starg = Smin;
  }
  if(Starg > Smax) { /* yet another rare case, happens only if fil_calcs >= (0.5 * nonfil_calcs) / (xhmm+1), when xhmm = 1, this is when fil_calcs <= nonfil_calcs / 4 */
    Starg = Smax;
  }
  fwd_Emax         = (Smax  * (float) dbsize) / surv_res_per_hit;
  fwd_Etarg        = (Starg * (float) dbsize) / surv_res_per_hit;
  fwd_Emin         = (Smin  * (float) dbsize) / surv_res_per_hit;

  /* Copy bit scores to three separate quicksort-able data structures, 
   * Sort one by CM E-value, one by Inside E-value and one by Forward E-value.
   * The seq idx is kept in the data structure, providing a link between 
   * the three different lists after they're sorted. 
   */
  ESL_ALLOC(by_cmA, sizeof(fseq_Eval_t) * filN);
  ESL_ALLOC(by_fwdA, sizeof(fseq_Eval_t) * filN);
  /* convert bit scores to E-values and copy them to the qsortable structures */
  for(i = 0; i < filN; i++) { 
    p = partA[i];
    cm_E  = RJK_ExtremeValueE(cm_scA[i],  cfg->cmstatsA[cmi]->gumAA[cm_mode][p]->mu,  cfg->cmstatsA[cmi]->gumAA[cm_mode][p]->lambda); 
    fwd_E = RJK_ExtremeValueE(fwd_scA[i], cfg->cmstatsA[cmi]->gumAA[hmm_fwd_mode][p]->mu, cfg->cmstatsA[cmi]->gumAA[hmm_fwd_mode][p]->lambda); 
    /* E-values correspond to expected number of hits in a seq DB of length evalue_L, convert that DB size to dbsize */
    cm_E  *= dbsize / evalue_L;
    fwd_E *= dbsize / evalue_L;
    /* copy E-values to qsortable data structures */
    by_cmA[i].i     = by_fwdA[i].i     = i;
    by_cmA[i].cm_E  = by_fwdA[i].cm_E  = cm_E;
    by_cmA[i].fwd_E = by_fwdA[i].fwd_E = fwd_E;
  }
  /* qsort */
  qsort(by_cmA,  filN, sizeof(fseq_Eval_t),  compare_fseq_by_cm_Eval);
  qsort(by_fwdA, filN, sizeof(fseq_Eval_t), compare_fseq_by_fwd_Eval);

  /* determine the mapping between the sorted CM list and sorted Fwd list.
   * This is N^2 which might be of concern with large filN's 
   * empirically for filN = 10,000 this block takes 1.25s, 
   * so for filN == 100,000 it would take about 2 minutes. But
   * for filN == 100,000 the time search 100,000 seqs will drown
   * this out. Nevertheless, max <n> for --fil-N option is 100,000.
   */
  ESL_ALLOC(cmi2fwdi, sizeof(int) * filN);
  for(i = 0; i < filN; i++) {
    for(j = 0; j < filN; j++) {
      if(by_cmA[i].i == by_fwdA[j].i) cmi2fwdi[i] = j;
    }
  }
  /* A few more preparations before the main loop:
   * we want Fidx such that in the sorted arrays, by_cmA and by_fwdA, F fraction of the elements 
   * have E value <= score of element [Fidx], including element Fidx. Some examples:
   * F      filN    Fidx
   * 0.99    100       1  we'll miss 0..0   ==   1 hit  at this threshold, but recognize   99 (0.99)
   * 0.93   2000     140  we'll miss 0..139 == 140 hits at this threshold, but recognize 1860 (0.93)
   * 0.81    243      46  we'll miss 0..45  ==  46 hits at this threshold, but recognize  197 (0.8106)
   */
  Fidx = (int) ((1. - F) * (float) filN);
  /* deal with a precision issue */
  if(Fidx < ((1. - F) * (float) filN)) { /* we rounded down */
    if (((Fidx - ((1. - F) * (float) filN)) - 1.) < eslSMALLX1) Fidx++; /* due to precision issues we unnecessarily rounded down */
  }

  /* Setup fwd_useme array: fwd_useme[i] is TRUE if element with rank i in by_fwdA array is 'in use', 
   * meaning it corresponds to a sequence with CM score better than current threshold */
  int *fwd_useme;      /* [0..i..filN] TRUE if element with rank i in by_fwd array is 'in use', meaning it corresponds to a sequence with CM score better than current threshold */
  ESL_ALLOC(fwd_useme, sizeof(int) * filN);
  esl_vec_ISet(fwd_useme, filN, TRUE);
  int fwd_all = 0;     /* the index in by_fwdA of the worst ranked  (highest) forward E-value of all elements in use (fwd_useme[i] == TRUE) */
  int fwd_F   = Fidx;  /* the index in by_fwdA of the sequence with forward E-value worse than F fraction of all elements in use (fwd_useme[i] == TRUE) */

  /* Precalculate the i values for which we have to change fwd_F, as step through the i=0..imax loop.
   * when i == 0, curN == filN-i+1 == filN, and fwd_F == Fidx because we want to find F fraction of the filN hits
   * but as i increases, curN decreases (as we step through the ranked CM E-value threshold) and the integer 
   * corresponding  to F fraction of (filN-i+1) also changes, but due to precision difference of ints and floats,
   * it's non-trivial, and worth precalculating to avoid putting the messy code that does
   * it in the loop below. 
   */
  int curN    = filN;
  int *change_fwd_F;
  int prv_Fidx;
  int cur_Fidx;
  ESL_ALLOC(change_fwd_F, sizeof(int) * filN);
  prv_Fidx = Fidx;
  for(i = filN-1; i >= 0; i--) { 
    cur_Fidx = (int) ((1. - F) * (float) i);
    if(cur_Fidx < ((1. - F) * (float) i)) { /* we rounded down */
      if (((cur_Fidx - ((1. - F) * (float) i)) - 1.) < eslSMALLX1) cur_Fidx++; /* due to precision issues we unnecessarily rounded down */
    }
    if(prv_Fidx == cur_Fidx) change_fwd_F[i+1] = TRUE;
    else                     change_fwd_F[i+1] = FALSE;
    prv_Fidx = cur_Fidx;
  }
  change_fwd_F[0] = FALSE; /* we'll never access this element */

  /* The main loop:
   * step through 3/4 of the CM E-values starting with worst (so we'll stop at the hit with rank imax = (0.75 * filN)),
   * so we step through i=0..imax, and for each i, determine fwd_E_F[i], fwd_E_all[i], and fwd_E_cut[i]
   * fwd_E_F[i]:   the forward E-value cutoff that will recognize F fraction of CM hits with E-values better than or equal to by_cmA[i].cm_E
   * fwd_E_all[i]: the forward E-value cutoff that will recognize all           CM hits with E-values better than or equal to by_cmA[i].cm_E
   * fwd_E_cut[i]: the forward E-value cutoff we would report for this i, may be different than fwd_E_F[i] and fwd_E_all[i]
   *               because we consider our fwd_Emax, fwd_Etarg, and fwd_Emin values.
   * fwd_E_S[i]:   predicted survival fraction using fwd_E_cut[i], equals (fwd_E_cut[i] * surv_res_per_hit / dbsize).
   */
  int imax = (int) ((0.750001) * (float) filN); /* .750001 is to avoid unnec rounding down (for ex, if filN == 100, and we used 0.75, we may get 74.99999999 and end up rounding down) */
  float *fwd_E_F; 
  float *fwd_E_all;
  float *fwd_E_cut;
  float *fwd_E_S;  
  ESL_ALLOC(fwd_E_F,    sizeof(float) * (imax+1));
  ESL_ALLOC(fwd_E_all,  sizeof(float) * (imax+1));
  ESL_ALLOC(fwd_E_cut,  sizeof(float) * (imax+1));
  ESL_ALLOC(fwd_E_S,    sizeof(float) * (imax+1));
  int imax_above_Emax  = -1; /* at end of loop, imax_above_Emax is max i for which fwd_E_cut[i] >  fwd_Emax,  if -1, fwd_E_cut[i] < fwd_Emax  for all i=0..imax */
  int imin_at_Etarg    = -1; /* at end of loop, imin_at_Etarg   is min i for which fwd_E_cut[i] == fwd_Etarg, if -1, fwd_E_cut[i] > fwd_Etarg for all i=0..imax */
  int imax_at_Etarg    = -1; /* at end of loop, imax_at_Etarg   is max i for which fwd_E_cut[i] == fwd_Etarg, if -1, fwd_E_cut[i] > fwd_Etarg for all i=0..imax */
  int imin_at_Emin     = -1; /* at end of loop, imin_at_Emin    is min i for which fwd_E_cut[i] == fwd_Emin,  if -1, fwd_E_cut[i] > fwd_Emin  for all i=0..imax */

  /* TEMPORARY */
  ESL_DPRINTF1(("fwd_Emax:  %f\n", fwd_Emax));
  ESL_DPRINTF1(("fwd_Etarg: %f\n", fwd_Etarg));
  ESL_DPRINTF1(("fwd_Emin: %f\n", fwd_Emin));

  for(i = 0; i <= imax; i++) {
    fwd_E_all[i] = by_fwdA[fwd_all].fwd_E; /* by_fwdA[fwd_all] is HMM threshold that recognizes ALL curN CM hits w/E value <= by_cmA[i].cm_E */
    fwd_E_F[i]   = by_fwdA[fwd_F].fwd_E;   /* by_fwdA[fwd_F]   is HMM threshold that recognizes F fraction of curN CM hits w/E value <= by_cmA[i].cm_E */

    /* based on fwd_E_all[i] and fwd_E_F[i], determine what our cutoff should be */
    if(fwd_E_F[i] < fwd_Etarg) { /* if TRUE, using E-value that achieves the target survival fraction as cutoff will recognize F fraction of CM hits */
      fwd_E_cut[i]    = ESL_MIN(fwd_Etarg,    fwd_E_all[i]); /* take minimum of: E value that exactly satisfies target survival fraction, and E value that recognizes all CM hits */
      fwd_E_cut[i]    = ESL_MAX(fwd_E_cut[i], fwd_Emin);     /* never go less than fwd_Emin */
    }
    else fwd_E_cut[i] = fwd_E_F[i];       /* we didn't achieve our target survival fraction */
    fwd_E_S[i] = (fwd_E_cut[i] * (float) surv_res_per_hit) / (float) dbsize;

    ESL_DPRINTF1(("i: %3d N: %5d cm: %g fwd_E_all: %g fwd_E_F: %g cut: %g S: %g\n", i, curN, by_cmA[i].cm_E, fwd_E_all[i], fwd_E_F[i], fwd_E_cut[i], fwd_E_S[i])); 

    if(fwd_E_cut[i]  > fwd_Emax) imax_above_Emax++; /* we didn't even achieve our maximum allowed survival fraction */
    if(imin_at_Emin  == -1 && (fabs(fwd_E_cut[i] - fwd_Emin)  < eslSMALLX1)) { ESL_DPRINTF1(("\tAchieved Smin\n"));   imin_at_Emin  = i;   }
    if(imin_at_Etarg == -1 && (fabs(fwd_E_cut[i] - fwd_Etarg) < eslSMALLX1)) { ESL_DPRINTF1(("\tAchieved Starg\n"));  imin_at_Etarg = i;   }
    if(             imin_at_Etarg != -1 && imax_at_Etarg == -1 && (fabs(fwd_E_cut[i] - fwd_Etarg) > eslSMALLX1)) { ESL_DPRINTF1(("\tDone with Starg\n")); imax_at_Etarg = i-1; }
    if(i == imax && imin_at_Etarg != -1 && imax_at_Etarg == -1)                                                  { ESL_DPRINTF1(("\tBoundary case, final point is at Starg\n")); imax_at_Etarg = imax; }

    /* Now manipulate our lists for next step */
    /* mark hit cm rank i as 'no longer used' */
    fwd_useme[cmi2fwdi[i]] = FALSE;
    /* if hit cm rank i was the worst scoring fwd hit, update fwd_all to the new worst fwd hit (now that prev worst is no longer used) */
    if(cmi2fwdi[i] == fwd_all) {
      fwd_all++;
      while(!(fwd_useme[fwd_all])) fwd_all++; /* skip all 'no longer used' hits */
    }
    /* if hit cm rank i was within the (1.-F) fraction of worst hits, update fwd_F to new idx in sorted fwd hits that will recognize F fraction of CM hits */
    if(cmi2fwdi[i] <= fwd_F) { /* should we increment fwd_F? */
      if((change_fwd_F[curN]) && ((fwd_F+1) < filN)) { 
	fwd_F++;
	while((!(fwd_useme[fwd_F])) && ((fwd_F+1) < filN)) fwd_F++; /* skip all 'no longer used' hits */
      } 
    }
    curN--;
  }
  printf("imax_above_Emax: %d\n", imax_above_Emax);
  printf("imin_at_Etarg:   %d\n", imin_at_Etarg);
  printf("imax_at_Etarg:   %d\n", imax_at_Etarg);
  printf("imin_at_Emin:    %d\n", imin_at_Emin);
  ESL_DPRINTF1(("\n"));

  /*************TEMPORARY BLOCK************************/
  /* paranoid, expensive check */
  int error_flag = FALSE;
  int nmissed = 0;
  for(i = 0; i <= imax; i++) {
    nmissed = 0;
    for(j = i; j < filN; j++) { 
      if(by_cmA[j].fwd_E > fwd_E_all[i]) { 
	error_flag = 1; 
	printf("ERROR: i: %d j: %d fwd_E_all[i]: %g < fwd_E[j]: %g\n", i, j, fwd_E_all[i], by_cmA[i].fwd_E);
      }
      if(by_cmA[j].fwd_E > fwd_E_F[i]) { 
	nmissed++;
      }
    }
    if(((float) nmissed/((float) filN-i)) > F) { 
      error_flag = 1;
      printf("ERROR: i: %d nmissed: %d Fmissed: %f F: %f\n", i, nmissed, ((float) nmissed/((float)filN-i)), F);
    }
  }
  if(error_flag) cm_Fail("Implementation error dude.");
  ESL_DPRINTF1(("Passed expensive paranoia check\n"));
  /*************END OF TEMPORARY BLOCK************************/

  int ip_min = imax_above_Emax == -1 ? 0    : imax_above_Emax+1;
  int ip_max = imin_at_Emin    == -1 ? imax : imin_at_Emin;
  int ip = ip_min;
  float max_next_E_cut;
  float prev_E_cut;
  int *saveme;
  int keep_going;
  ESL_ALLOC(saveme, sizeof(int) * (ip_max + 1));
  esl_vec_ISet(saveme, (ip_max+1), FALSE);
  i = 0;
  printf("ip_min: %d\nip_max: %d\n", ip_min, ip_max);
  while(ip <= ip_max) {
    saveme[ip] = TRUE;
    max_next_E_cut    = fwd_E_cut[ip] * 0.9; /* we'll skip all subsequenct E-value cutoffs that are within 10% of the one we've just added */
    prev_E_cut        = fwd_E_cut[ip];

    printf("i: %5d fwd_E_cut[ip: %5d]: %12g CM_cut: %12g S: %12g max_next: %12g ", i, ip, fwd_E_cut[ip], by_cmA[ip].cm_E, fwd_E_S[ip], max_next_E_cut); 
    if     (fabs(fwd_E_cut[ip] - fwd_Emin)  < eslSMALLX1) printf(" [Emin]\n");
    else if(fabs(fwd_E_cut[ip] - fwd_Etarg) < eslSMALLX1) printf(" [Etarg]\n");
    else if(fabs(fwd_E_cut[ip] - fwd_Emax)  < eslSMALLX1) printf(" [Emax]\n");
    else printf("\n");

    i++;
    ip++;

    keep_going = TRUE;
    while(keep_going) { 
      if(fwd_E_cut[ip] < max_next_E_cut) keep_going = FALSE; /* the difference between this E-value cut, and the prev one we added is > 10% */ 
      if(ip == imin_at_Emin)             keep_going = FALSE; /* min i for which we acheive Smin,  we want to add this point no matter what */
      if(ip == imin_at_Etarg)            keep_going = FALSE; /* min i for which we acheive Starg, we want to add this point no matter what */
      if(ip >  imax)                     keep_going = FALSE; /* we're done with the list */
      if(keep_going) ip++; /* skip all points that fail 4 ifs() above */
    }
  }
  printf("final i: %d\n", i);

  float *fwd_E_cut2save;
  float *cm_E_cut2save;
  int n2save = i;
  ESL_ALLOC(fwd_E_cut2save, sizeof(float) * n2save);
  ESL_ALLOC(cm_E_cut2save,  sizeof(float) * n2save);
  i = 0;
  for(ip = ip_min; ip <= ip_max; ip++) {
    if(saveme[ip]) { 
      fwd_E_cut2save[i] = fwd_E_cut[ip];
      cm_E_cut2save[i]  = by_cmA[ip].cm_E;
      i++;
    }
  }
  assert(i == n2save);

#if 0  
  /* Now we want to store the CM/HMM cutoff pairs as efficiently as possible into 50 data points.
   * These 50 points will be a subset of the imax=0.75 * filN cutoff points, but deciding
   * which ones is complex. The choice of 50 points is arbitrary, we could use more but
   * these points will be stored in the CM file and I don't want to use too many. 
   * 
   * First, we can determine the minimum and maximum possible points:
   *
   * For the minimum, we don't want to use any points i for which fwd_E_cut[i] > fwd_Emax, 
   * because these points only tell us we shouldn't use a filter for CM E cutoffs >= by_cmA[i].cm_E.
   *
   * imax_above_Emax is the maximum i for which fwd_E_cut[i] <= fwd_Emax; if imax_above_Emax == -1, then
   * fwd_E_cut[i] <= fwd_Emax for all i.
   *
   * For the maximum, we only want to use a single point for which fwd_E_cut[i] < fwd_Emin,
   * namely the minimal i for which fwd_E_cut[i] < fwd_Emin, this point is imax_below_Emin, 
   * given this point we know that for any CM E cutoff less than by_cmA[imax_below_Emin] we
   * can set HMM cutoff as fwd_Emin and achieve or best possible scenario of HMM+CM running time 
   * 1.1X HMM running time. If imax_below_Emin == -1, no such i exists, fwd_E_cut[i] > fwd_Emin for all i.
   *
   * Let's define ip_min = imax_above_Emax == -1 ? 0    : imax_above_Emin;
   *              ip_max = imin_at_Emin    == -1 ? imax : imin_at_Emin;
   * 
   * Given the min and max, the remaining task is to set the 48 points in between:
   * 
   * Easy case: ip_max - ip_min + 1 <= 48, we just use all points. 
   * Hard case: ip_max - ip_min + 1 >  48.
   *
   * A case we should deal with specially is that it's possible that many 
   * points will have fwd_E_cut[i] == fwd_E_targ, the range of such points is
   * imin_at_Etarg .. imax_at_Etarg. If both of these values is -1, then 
   * fwd_E_cut[i] > fwd_Etarg for all i. 
   * (It should never happen that one of imin_at_Etarg and imax_at_Etarg is -1 
   * but the other is not).
   * This case is special b/c we only need to store 1 point for fwd_Etarg, the
   * maximal CM E-value cutoff for which fwd_Etarg is achieved, this is
   * by_cmA[imin_at_Etarg].cm_E.
   *
   * So if imin_at_Etarg != -1, it is our third point, so we have 47 to work with.
   * ip_mid is our third point, and we're left with 47. 
   * These 47 points should cover the space of CM E-value cutoffs in between
   * our end points in some well-principled way. 
   * 
   * The space we need to cover is from:
   * 
   * if(ip_mid != -1): 
   * 
   * by_cmA[ip_min].cm_E        to by_cmA[imin_at_Etarg].cm_E and
   * by_cmA[imax_at_Etarg].cm_E to by_cmA[ip_max].cm_E 
   * 
   * else if(ip_mid == -1)
   *
   * by_cmA[ip_min].cm_E to by_cmA[ip_max].cm_E
   *
   * 
   *  
   */
  /*******************************************************************************/
  /* 
   * want 50 bins that will hold (int) ip_max-ip_min+1/50. scores each
   */
  int always_less_than_Emax;
  int achieved_Emin;

  if(imin_above_Emax == -1) { /* all fwd_E_cut[i]'s are < fwd_Emax */
    always_less_than_Emax = TRUE;
    imin_above_Emax = 0;
  }
  else always_less_than_Emax = FALSE; /* imin_above_Emax will be maximum i for which fwd_E_cut[i] > fwd_Emax */
  if(imax_below_Emin == -1) { /* fwd_E_cut[i] > fwd_Emin for all i */
    achieved_Emin = FALSE;
    imax_below_Emin = imax;
  }
  else achieved_Emin = TRUE; /* imax_below_Emin will be minimum i for which fwd_E_cut[i] == fwd_Emin */
  

  int ip_min = imin_above_Emax;
  int ip_max = imax_below_Emin;
  float *cm_E2print;
  float *fwd_E_cut2print;
  int nbins = ESL_MIN(50., ((float) ip_max - ip_min + 1.));
  ESL_ALLOC(cm_E2print,     sizeof(float) * nbins);
  ESL_ALLOC(fwd_E_cut2print, sizeof(float) * nbins);
  float fbinsize = (float) ((ip_max - ip_min + 1.)) / nbins;
  int ip = ip_min;
  float ip_float = (float) ip_min;
  i = 0;
  while(ip < ip_max) {
    if(i >= nbins) { printf("ERROR, i: %d > nbins: %d\n", i, nbins); }
    cm_E2print[i]      = by_cmA[ip].cm_E;
    fwd_E_cut2print[i] = fwd_E_cut[ip];
    printf("ip: %d ip_float: %.3f i: %d cm: %g fwd: %g S: %g\n", ip, ip_float, i, cm_E2print[i], fwd_E_cut2print[i], fwd_E_S[ip]);
    ip_float += fbinsize;
    ip = (int) ip_float;
    i++;
  }
#endif

  return eslOK;

 ERROR:
  return status; 
}

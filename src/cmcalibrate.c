/************************************************************
 * @LICENSE@
 ************************************************************/
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
 */
#include "esl_config.h"
#include "config.h"	

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <time.h>

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
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

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
  { "--dbfile",  eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--filonly", "use GC content distribution from file <s>",  2},
  { "--gumhfile",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save fitted Gumbel histogram(s) to file <s>", 2 },
  { "--gumqqfile",eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save Q-Q plot for Gumbel histogram(s) to file <s>", 2 },
  /* options for filter threshold calculation */
  { "--filonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "only calculate filter thresholds, don't estimate Gumbels", 3},
  { "--filN",    eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--gumonly", "number of emitted sequences for HMM filter threshold calc",    3 },
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
  char           *cmfile;	      /* name of input CM file  */ 
  ESL_RANDOMNESS *r;
  ESL_ALPHABET   *abc;
  double         *gc_freq;
  double         *pgc_freq;
  int             be_verbose;	
  ESL_STOPWATCH  *w;
  CMStats_t     **cmstatsA;           /* the CM stats data structures, 1 for each CM */
  int             ncm;                /* what number CM we're on */
  int             cmalloc;            /* number of cmstats we have allocated */
  char           *tmpfile;            /* tmp file we're writing to */
  char           *mode;               /* write mode, "w" or "wb"                     */
  float           minlen;             /* minimum average length of a subseq for calc'ing gumbels */
  long            dbsize;             /* size of DB for gumbel stats (impt for E-value cutoffs for filters) */ 
  float          *cutoffA;            /* bit score cutoff for each partition, changes to reflect
				       * current mode CM is in, on masters and workers */

  int             do_mpi;
  int             my_rank;
  int             nproc;
  int             do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *gumhfp;        /* optional output for gumbel histograms */
  FILE         *gumqfp;        /* optional output for gumbel QQ file */
  FILE         *filhfp;        /* optional output for filter histograms */
  FILE         *filrfp;        /* optional output for filter info for R */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit Gumbels for E-value stats and calculate HMM filter threshold stats";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
/* static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); */

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif
static int process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			    int emit_from_cm, float ***ret_vscAA, float **ret_cp9scA);
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int update_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int fthr_mode);
static int set_partition_gc_freq(struct cfg_s *cfg, int p);
static int find_possible_filter_states(struct cfg_s *cfg, CM_t *cm, int **ret_fitme, float **ret_avglen);
static int calc_avg_hit_length(CM_t *cm, float **ret_hitlen, double beta);
static int pick_stemwinners(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, double *muA, double *lambdaA, float *avglen, int **ret_vwin, int *ret_nwin);
static int fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, double *ret_mu, double *ret_lambda);
static int cm_fit_histograms(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **vscA, int nscores, int *lbegin, double **ret_muA, double **ret_lambdaA);
static ESL_DSQ * get_random_dsq(const struct cfg_s *cfg, CM_t *cm, double *dnull, int L);
static ESL_DSQ * get_cmemit_dsq(const struct cfg_s *cfg, CM_t *cm, int *ret_L, int *ret_p);
static float search_target_cm_calibration(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float **ret_vsc);
int
main(int argc, char **argv)
{

  ESL_GETOPTS     *go	   = NULL;     /* command line processing                     */
  struct cfg_s     cfg;

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
    esl_fatal("--ga not yet implemented, implement it.");
  if (! esl_opt_IsDefault(go, "--nc"))
    esl_fatal("--nc not yet implemented, implement it.");
  if (! esl_opt_IsDefault(go, "--tc"))
    esl_fatal("--tc not yet implemented, implement it.");

  /* Initialize configuration shared across all kinds of masters
   * and workers in this .c file.
   */
  cfg.cmfile  = esl_opt_GetArg(go, 1);
  if (cfg.cmfile == NULL) esl_fatal("Failed to read <cmfile> argument from command line.");
  cfg.cmfp     = NULL;
  cfg.gc_freq  = NULL; 
  cfg.pgc_freq = NULL; 
  cfg.r        = NULL; 
  cfg.w        = NULL; 
  cfg.ncm      = 0;
  cfg.cmstatsA = NULL;
  cfg.tmpfile  = NULL;
  cfg.mode     = NULL;
  cfg.minlen   = 7.;      /* not user changeable, should it be? is 8 not good? */
  cfg.dbsize   = 1000000; /* default DB size = 1MB, changed ONLY if --dbfile enabled */
  cfg.cutoffA  = NULL; 

  cfg.gumhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gumqfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.filhfp   = NULL; /* ALWAYS remains NULL for mpi workers */
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

      if(cfg.nproc == 1) cm_Fail("ERROR, MPI mode, but only 1 processor running...");

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
 *    cfg->cmfp         - open CM file                
 *    cfg->gumhfp      - optional output file
 *    cfg->gumqfp      - optional output file
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

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, NULL, "Failed to open covariance model save file %s\n", cfg->cmfile);

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
    if (status == eslENOTFOUND)    esl_fatal("No such file."); 
    else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
    else if (status == eslEINVAL)  esl_fatal("Can’t autodetect stdin or .gz."); 
    else if (status != eslOK)      esl_fatal("Failed to open sequence database file, code %d.", status); 
    GetDBInfo(tmp_abc, dbfp, &(cfg->dbsize), &(cfg->gc_freq));
    esl_vec_DNorm(cfg->gc_freq, GC_SEGMENTS);
    esl_alphabet_Destroy(cfg->abc);
    esl_sqfile_Close(dbfp);
    /* allocate pgc_freq, the gc freqs per partition, used to sample seqs for different partitions */
    ESL_ALLOC(cfg->pgc_freq, sizeof(double) * GC_SEGMENTS);
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


/* serial_master()
 * The serial version of cmcalibrate.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int            status;
  char           errbuf[eslERRBUFSIZE];
  CM_t          *cm = NULL;
  int            cmN  = esl_opt_GetInteger(go, "--cmN");
  int            hmmN = esl_opt_GetInteger(go, "--hmmN");
  double        **muA = NULL;  /* [0..np-1][0..cm->M-1], mu for each partition, each state, 
				* muA[][v] == -1. if we're not fitting state v to a gumbel */
  double    **lambdaA = NULL;  /* same as muA, but lambda */
  /*int           *vwin = NULL;
    int            nwin = 0;*/
  int          *fitme = NULL;
  float       *avglen = NULL;
  int        gum_mode = 0;
  int               p;
  float **gum_vscAA = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each random seq */
  float **fil_vscAA = NULL; /* [0..v..cm->M-1][0..nseq-1] best cm score for each state, each emitted seq */
  float   *gum_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each random seq */
  float   *fil_cp9scA = NULL; /*                [0..nseq-1] best cp9 score for each emitted seq */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      if((status = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf))                   != eslOK) cm_Fail(errbuf);
      if((status = find_possible_filter_states(cfg, cm, &fitme, &avglen)) != eslOK) cm_Fail(errbuf);

      //for(gum_mode = 0; gum_mode < NGUMBELMODES; gum_mode++) {
      for(gum_mode = 0; gum_mode < 1; gum_mode++) 
	{
	  ConfigForGumbelMode(cm, gum_mode);
	  ESL_ALLOC(muA,     sizeof(double *) * cfg->cmstatsA[cfg->ncm-1]->np);
	  ESL_ALLOC(lambdaA, sizeof(double *) * cfg->cmstatsA[cfg->ncm-1]->np);
	  for (p = 0; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) 
	    {
	      if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	      if(gum_mode < NCMMODES) 
		{ /* CM mode */
		  /* search random sequences to get gumbels */
		  if((status = process_workunit (go, cfg, errbuf, cm, cmN,   FALSE, &gum_vscAA, NULL))        != eslOK) cm_Fail(errbuf);
		  if((status = cm_fit_histograms(go, cfg, errbuf, cm, gum_vscAA, cmN, fitme, &(muA[p]), &(lambdaA[p]))) != eslOK) cm_Fail(errbuf);
		  cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->mu     = muA[p][0];
		  cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->lambda = lambdaA[p][0];
		  cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->L      = cm->W*2;
		  cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->N      = cmN;
		}
	      else 
		{ /* CP9 mode, only fit gumbels */
		  if((status = process_workunit (go, cfg, errbuf, cm, hmmN, FALSE, NULL, &gum_cp9scA))                 != eslOK) cm_Fail(errbuf);
		  if((status = fit_histogram(go, cfg, errbuf, gum_cp9scA, hmmN, &(cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->mu), 
					     &(cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->lambda)))     != eslOK) cm_Fail(errbuf);
		  cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->L      = cm->W*2;
		  cfg->cmstatsA[cfg->ncm-1]->gumAA[gum_mode][p]->N      = hmmN;
		}
	    }
	  if(gum_mode < NCMMODES) {
	    /* update cfg, so cutoff[p] is bit score for current gum_mode, and i can access it inside
	     * process_workunit, by passing it to cm_emit_seqs_to_aln_above_cutoff() */
	    if((status = update_cutoffs(go, cfg, errbuf, cm, gum_mode)) != eslOK) cm_Fail(errbuf);
	    
	    /* search emitted sequences to get filter thresholds */
	    if((status = process_workunit (go, cfg, errbuf, cm, cmN, TRUE, &fil_vscAA, &fil_cp9scA)) != eslOK) cm_Fail(errbuf);
	    cm_Fail("done.");
	  }
	  /* free muA and lambdaA */
	}
    }
  return;
 ERROR:
  cm_Fail("serial_master(), memory allocation error.");
  return; /* NEVERREACHED */
}

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
 *           <emit_from_cm> is FALSE, ret_vscAA != NULL, ret_cp9scA == NULL
 *           Emit randomly and search only with the CM. <ret_vscAA> is filled
 *           with the best CM score at each state for each sequence.
 *
 *           Mode 2. Gumbel calculation for CP9.
 *           <emit_from_cm> is FALSE, ret_vscAA == NULL, ret_cp9scA != NULL
 *           Emit randomly and search only with the CP9. <ret_cp9scA> is filled
 *           with the best CP9 score for each sequence.
 *
 *           Mode 3. Scores will eventually be used to determine filter thresholds.
 *           <emit_from_cm> is TRUE, ret_vscAA != NULL, ret_cp9scA != NULL
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
  int            passed_flag;
  ESL_DSQ       *dsq;

  /* contract check */
  if(ret_vscAA  == NULL && ret_cp9scA == NULL) { sprintf(errbuf, "process_workunit, ret_vscAA and ret_cp9scA both NULL."); return eslEINVAL; } 
  if(ret_vscAA  != NULL && ret_cp9scA != NULL && !emit_from_cm) { sprintf(errbuf, "process_workunit, ret_vscAA and ret_cp9scA both non-NULL, but emit_from_cm is FALSE."); return eslEINVAL; } 

  /* determine mode, are we calc'ing CM stats, CP9 stats, or both */
  if(ret_vscAA   != NULL) use_cm  = TRUE;
  if(ret_cp9scA  != NULL) use_cp9 = TRUE;

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

  /* generate dsqs one at a time and collect best CM scores at each state and/or best overall CP9 score */
  for(i = 0; i < nseq; i++) 
    {
      passed_flag = TRUE;
      if(emit_from_cm)   dsq = get_cmemit_dsq(cfg, cm, &L, &p);
      else             { dsq = get_random_dsq(cfg, cm, dnull, cm->W*2); L = cm->W*2; }
      /* if nec, search with CM */
      if (use_cm)
	{ 
	  search_target_cm_calibration(cm, dsq, cm->dmin, cm->dmax, 1, L, cm->W, &(cur_vscA)); 
	  for(v = 0; v < cm->M; v++) vscAA[v][i] = cur_vscA[v];
	  free(cur_vscA);

	  /* if we're emitting from CM for filter stats, we have to ensure this dsq has a hit
	   * above our cutoff, if not we don't want to count it */
	  if(emit_from_cm && vscAA[0][i] < cfg->cutoffA[p]) {
	    search_target_cm_calibration(cm, dsq, NULL, NULL, 1, L, cm->W, &(cur_vscA)); 
	    if(cur_vscA[0] < cfg->cutoffA[p]) { passed_flag = FALSE; i--; } /* set flag to FALSE, decrement i */
	    free(cur_vscA);
	  }
	}
      /* if nec, search with CP9 */
      if (use_cp9 && passed_flag) cp9scA[i] = CP9Forward(cm, dsq, 1, L, cm->W, 0., 
							 NULL,   /* don't return scores of hits */
							 NULL,   /* don't return posns of hits */
							 NULL,   /* don't keep track of hits */
							 TRUE,   /* we're scanning */
							 FALSE,  /* we're not ultimately aligning */
							 FALSE,  /* we're not rescanning */
							 TRUE,   /* be memory efficient */
							 NULL);  /* don't want the DP matrix back */
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
  if(!(esl_opt_GetBoolean(go, "--iins"))) cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  cm->config_opts |= CM_CONFIG_QDB;   /* configure QDB, only to get W */
  cm->search_opts |= CM_SEARCH_NOQDB; /* don't use QDB to search */
  ConfigCM(cm, NULL, NULL);
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

  cfg->cmstatsA[cfg->ncm-1] = AllocCMStats(1); /* Always 1 partition (TEMPORARY) */
  cfg->cmstatsA[cfg->ncm-1]->ps[0] = 0;
  cfg->cmstatsA[cfg->ncm-1]->pe[0] = 100;
  for(i = 0; i < GC_SEGMENTS; i++)
    cfg->cmstatsA[cfg->ncm-1]->gc2p[i] = 0; 

  if(cfg->cutoffA != NULL) free(cfg->cutoffA);
  ESL_ALLOC(cfg->cutoffA, sizeof(float) * 1); /* number of partitions */
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

/* Function: find_possible_filter_states()
 * Date:     EPN, Mon Sep 10 15:49:43 2007
 *
 * Purpose:  Given a CM, fill an array of length cm->M with TRUE, FALSE for each
 *           state. TRUE if we should fit a Gumbel to the state b/c it may be
 *           a good sub-CM filter root state, or FALSE if not. Criteria is that
 *           a state must be a possible local entry state AND must have an average
 *           subseq length of cfg->minlen.
 *
 * Returns:  eslOK on success;
 */
int
find_possible_filter_states(struct cfg_s *cfg, CM_t *cm, int **ret_fitme, float **ret_avglen)
{
  int    status;
  int   *fitme  = NULL;
  float *avglen = NULL;
  int nd;
  int v;

  if((status = calc_avg_hit_length(cm, &avglen, 1e-15)) != eslOK) return status;

  ESL_ALLOC(fitme, sizeof(int) * cm->M);
  esl_vec_ISet(fitme, cm->M, FALSE);

  fitme[0] = TRUE; /* ROOT_S: has to be true */
  for (nd = 1; nd < cm->nodes; nd++) /* note: nd 1 must be MATP, MATL, MATR, or BIF */
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
      {
	v = cm->nodemap[nd];
	if(avglen[v] > cfg->minlen) fitme[v] = TRUE;
      }
  *ret_avglen = avglen;
  *ret_fitme  = fitme;
  return eslOK;

 ERROR:
  return status;
}

/* Function: calc_avg_hit_length()
 * Date:     EPN, Mon Sep 10 10:04:52 2007
 *
 * Purpose:  Calculate the average hit length for each state of a CM using the 
 *           QDB calculation engine, using beta = 1e-15.
 *
 * Returns:  eslOK on success;
 */
int
calc_avg_hit_length(CM_t *cm, float **ret_avglen, double beta)
{
  float *avglen = NULL;
  int safe_windowlen;

  safe_windowlen = cm->W * 2;
  while(!(BandCalculationEngine(cm, safe_windowlen, beta, TRUE, NULL, NULL, NULL, &avglen))) {
    safe_windowlen *= 2;
    if(safe_windowlen > (cm->clen * 1000)) cm_Fail("ERROR safe_windowlen big: %d\n", safe_windowlen);
  }

  int v;
  for(v = 0; v < cm->M; v++) printf("AVG LEN v: %4d d: %10.4f\n", v, avglen[v]);

  *ret_avglen = avglen;
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
  int i;

  esl_vec_DSet(cfg->pgc_freq, GC_SEGMENTS, 0.);
  for (i = cfg->cmstatsA[cfg->ncm-1]->ps[p]; i < cfg->cmstatsA[cfg->ncm-1]->pe[p]; i++) 
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
		  float **vscA, int nscores, int *fitme, double **ret_muA, double **ret_lambdaA)
{
  int status;
  double *muA = NULL;
  double *lambdaA = NULL;
  int v;

  ESL_ALLOC(muA,     sizeof(double) * cm->M);
  ESL_ALLOC(lambdaA, sizeof(double) * cm->M);

  for(v = 0; v < cm->M; v++) {
    if(fitme[v]) {
      printf("FITTING v: %d sttype: %d\n", v, cm->sttype[v]);
      fit_histogram(go, cfg, errbuf, vscA[v], nscores, &(muA[v]), &(lambdaA[v]));
    }
    else muA[v] = lambdaA[v] = -1.;
  }

  *ret_muA     = muA;
  *ret_lambdaA = lambdaA;
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
  if (esl_rnd_xIID(cfg->r, dnull, cm->abc->K, L, dsq) != eslOK) cm_Fail("get_random_dsq(): failure creating random sequence.");

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
get_cmemit_dsq(const struct cfg_s *cfg, CM_t *cm, int *ret_L, int *ret_p)
{
  int p;
  int L;
  ESL_SQ *sq;
  ESL_DSQ *dsq;

  EmitParsetree(cm, cfg->r, "irrelevant", TRUE, NULL, &sq, &L);
  /* determine the partition */
  p = cfg->cmstatsA[cfg->ncm-1]->gc2p[(get_gc_comp(sq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */

  /* free everything allocated by a esl_sqio.c:esl_sq_CreateFrom() call, but the dsq */
  dsq = sq->dsq;
  free(sq->name);
  free(sq->acc);
  free(sq->desc);
  free(sq);

  *ret_L = L;
  *ret_p = p;
  return sq->dsq;
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
		vsc[v] = ESL_MAX(vsc[0], alpha[v][jp_v][d]);
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
		ivsc[v] = ESL_MAX(ivsc[0], ialpha[v][jp_v][d]);
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

  if (ret_vsc != NULL) *ret_vsc = vsc;
  else free(vsc);

  return vsc[0];

  ERROR:
    cm_Fail("Memory allocation error.\n");
    return 0.; /* NEVERREACHED */
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

  if(cm->name != NULL) namelen = strlen(cm->name) + 50; /* 50 digit int is considered max, sloppy. */
  else                 namelen = 100;                   /* 50 digit int is considered max, sloppy. */
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
      cm->search_opts &~ CM_SEARCH_HBANDED;
      cm->search_opts &~ CM_SEARCH_HBANDED;
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
#endif

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

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_vectorops.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"
#include "esl_exponential.h"
#include "esl_random.h"
#include "esl_gumbel.h"
#include "esl_stopwatch.h"

#include "structs.h"
#include "funcs.h"		/* external functions                   */
#include "stats.h"              /* gumbel functions */
#include "cm_dispatch.h"	
#include "mpifuncs.h"	

#define CUTOPTS  "--eval,--ga,--nc,--tc,--all"  /* Exclusive choice for filter threshold score cutoff */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-s",        eslARG_INT,      "0", NULL, "n>=0",    NULL,      NULL,        NULL, "set random number seed to <n>",   1 },
  { "-t",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print timings for Gumbel fitting and CP9 filter calculation",  1},
  { "--lins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "learn insert emission scores, do not zero them",  1},

  { "--gumonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--filonly", "only estimate Gumbels, don't calculate filter thresholds", 2},
  { "--cmN",     eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CM gumbel estimation",    2 },
  { "--hmmN",    eslARG_INT,   "5000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CP9 HMM gumbel estimation",    2 },
  { "--dbfile",  eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--filonly", "use GC content distribution from file <s>",  2},
  { "--gumhfile",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save fitted Gumbel histogram(s) to file <s>", 2 },
  { "--gumqqfile",eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save Q-Q plot for Gumbel histogram(s) to file <s>", 2 },

  { "--filonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "only calculate filter thresholds, don't estimate Gumbels", 3},
  { "--filN",    eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--gumonly", "number of emitted sequences for HMM filter threshold calc",    3 },
  { "--fract",   eslARG_REAL,  "0.95", NULL, "0<x<=1",  NULL,      NULL, "--gumonly", "required fraction of CM hits that survive HMM filter", 3},
  { "--targsurv",eslARG_REAL,  "0.01", NULL, "0<x<=1",  NULL,      NULL, "--gumonly", "target filter survival fraction", 3},
  { "--fast",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "calculate filter thr quickly, assume parsetree sc is optimal", 3},
  { "--gemit",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "when calc'ing filter thresholds, always emit globally from CM",  3},
  { "--filhfile",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--gumonly", "save CP9 filter threshold histogram(s) to file <s>", 3},
  { "--filrfile",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--gumonly", "save CP9 filter threshold information in R format to file <s>", 3},

  { "--eval",    eslARG_REAL,    "10", NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "min CM E-val (for a 1MB db) to consider for filter thr calc", 4}, 
  { "--ga",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM gathering threshold as minimum sc to consider for filter thr calc", 4}, 
  { "--nc",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM noise cutoff as minimum sc to consider for filter thr calc", 4}, 
  { "--tc",      eslARG_NONE,   FALSE, NULL, "x>0",  CUTOPTS,      NULL, "--gumonly", "use CM trusted cutoff as minimum sc to consider for filter thr calc", 4},   
  { "--all",     eslARG_NONE,   FALSE, NULL, NULL,   CUTOPTS,      NULL, "--gumonly", "accept all CM hits for filter calc, DO NOT use cutoff", 4}, 

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char           *cmfile;
  ESL_RANDOMNESS *r;
  double         *gc_freq;
  int             be_verbose;	

  int             do_mpi;
  int             my_rank;
  int             nproc;
  int             do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *gum_hfp;        /* optional output for gumbel histograms */
  FILE         *gum_qfp;        /* optional output for gumbel QQ file */
  FILE         *fil_hfp;        /* optional output for filter histograms */
  FILE         *fil_rfp;        /* optional output for filter info for R */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit Gumbels for E-value stats and calculate HMM filter threshold stats";

static int cm_fit_gumbel(CM_t *cm, ESL_GETOPTS *go, struct cfg_s *cfg, CMStats_t *cmstats, int gum_mode);

int
main(int argc, char **argv)
{
  int status;			       /* status of a function call                   */
  void            *tmp;                /* for ESL_RALLOC                              */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                     */
  struct cfg_s     cfg;
  CM_t            *cm      = NULL;     /* the covariance model                        */
  CMStats_t      **cmstats;            /* the CM stats data structures, 1 for each CM */
  char            *tmpfile;            /* temporary calibrated CM file                */
  FILE            *outfp;              /* for writing CM(s) into tmpfile              */
  char            *mode;               /* write mode, "w" or "wb"                     */
  int              idx;		       /* counter over CMs                            */
  sigset_t         blocksigs;	       /* list of signals to protect from             */
  ESL_STOPWATCH   *w;                  /* watch for timing GUM, filter thr calc times */
  int              ncm;                /* number of CMs read from cmfile              */
  int              cmalloc;            /* for alloc'ing CMStats_t objects             */
  int              i;                  /* counter over GC segments                    */
  int              continue_flag;      /* we still have a CM to get stats for?        */

  /* Gumbel related variables */
  int              gum_mode;           /* counter over gum modes                      */
  ESL_SQFILE      *dbfp;               /* open file pointer for dbfile                */
  int              format;	       /* format of dbfile                            */
  long             N;                  /* database size if --dbfile enabled           */

  /* CP9 HMM filtering threshold related variables */
  int              fthr_mode;          /* counter over fthr modes                     */
  int              db_size = 1000000;  /* we assume a 1MB db for filter thr calcs     */
  float            l_eval;             /*  local HMM E-value cutoff                   */
  float            l_F;                /* fraction of CM hits that are found w/l_eval */
  float            g_eval;             /* global HMM E-value cutoff                   */
  float            g_F;                /* fraction of CM hits that are found w/g_eval */
  int              emit_mode;          /* CM_LC or CM_GC, CM mode to emit with        */

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
  cfg.cmfile    = esl_opt_GetArg(go, 1);
  if (cfg.cmfile == NULL) esl_fatal("Failed to read <cmfile> argument from command line.");
  cfg.cmfp     = NULL;
  cfg.gc_freq  = NULL; 
  cfg.r        = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gum_hfp  = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.gum_qfp  = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.fil_hfp  = NULL; /* ALWAYS remains NULL for mpi workers */

  cfg.do_mpi   = FALSE;
  cfg.my_rank  = 0;
  cfg.nproc    = 0;


#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  /* Initialize MPI, figure out who we are, and whether we're running
   * this show (proc 0) or working in it (procs >0)
   */
  cfg.do_mpi = TRUE;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &cfg.my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &cfg.nproc);
#endif
  /* Get distribution of GC content from an input database */
  if(esl_opt_GetString(go, "--dbfile") != NULL)
    {
      ESL_ALPHABET *abc;
      abc = esl_alphabet_Create(eslRNA);
      status = esl_sqfile_Open(esl_opt_GetString(go, "--dbfile"), format, NULL, &dbfp);
      if (status == eslENOTFOUND)    esl_fatal("No such file."); 
      else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
      else if (status == eslEINVAL)  esl_fatal("Can’t autodetect stdin or .gz."); 
      else if (status != eslOK)      esl_fatal("Failed to open sequence database file, code %d.", status); 
      GetDBInfo(abc, dbfp, &N, &cfg.gc_freq);
      printf("Read DB file: %s of length %ld residues (both strands) for GC distro.\n", 
	     esl_opt_GetString(go, "--dbfile"), (2*N));
      esl_vec_DNorm(cfg.gc_freq, GC_SEGMENTS);
      esl_alphabet_Destroy(abc);
    }
  else /* use 0.25 A,C,G,U to generate random sequences */
    cfg.gc_freq = NULL;

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cmalloc  = 128;
  ESL_ALLOC(cmstats, sizeof(CMStats_t *) * cmalloc);
  ncm      = 0;
  
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  if(cfg.my_rank == 0) /* master only block 1 */
    {
#endif 
  /* Seed RNG */
  if(esl_opt_GetInteger(go, "-s") > 0)
    {
      if ((cfg.r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"))) == NULL)
	esl_fatal("Failed to create random number generator: probably out of memory");
    }
  else 
    {
      if ((cfg.r = esl_randomness_CreateTimeseeded()) == NULL)
	esl_fatal("Failed to create random number generator: probably out of memory");
    }
  printf("Random seed: %ld\n", esl_randomness_GetSeed(cfg.r));
  fflush(stdout);

  if ((cfg.cmfp = CMFileOpen(cfg.cmfile, NULL)) == NULL)
    esl_fatal("Failed to open covariance model save file %s\n", cfg.cmfile);

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
  ESL_ALLOC(tmpfile, (sizeof(char) * (strlen(cfg.cmfile) + 5)));
  strcpy(tmpfile, cfg.cmfile);
  strcat(tmpfile, ".xxx");	/* could be more inventive here... */
  if (esl_FileExists(tmpfile))
    esl_fatal("temporary file %s already exists; please delete it first", tmpfile);
  if (cfg.cmfp->is_binary) mode = "wb";
  else                     mode = "w"; 

  if(esl_opt_GetBoolean(go, "--time"))
    w = esl_stopwatch_Create(); 
  CMFileRead(cfg.cmfp, NULL, &cm);
  if(cm == NULL) esl_fatal("Failed to read a CM from %s -- file corrupt?\n", cfg.cmfile);

  /* Open output files if nec */
				/* histogram files, fitted and predicted */
  if (esl_opt_GetString(go, "--gumhfile") != NULL) 
    {
      if ((cfg.gum_hfp = fopen(esl_opt_GetString(go, "--gumhfile"), "w")) == NULL)
	esl_fatal("Failed to open gumbel histogram save file for writing\n");
    }
  if (esl_opt_GetString(go, "--gumqqfile") != NULL) 
    {
      if ((cfg.gum_qfp = fopen(esl_opt_GetString(go, "--gumqqfile"), "w")) == NULL)
	esl_fatal("Failed to open gumbel histogram save file for writing\n");
    }
  if (esl_opt_GetString(go, "--filhfile") != NULL) 
    {
      if ((cfg.fil_hfp = fopen(esl_opt_GetString(go, "--filhfile"), "w")) == NULL)
	esl_fatal("Failed to open filter histogram save file for writing\n");
    }
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
    }  /* end of master only block 1 */
#endif
  ncm = 0;
  continue_flag = 1; /* crudely used in MPI mode to make non-master MPIs go
		      * through the main loop for potentially multiple CMs */
  while (continue_flag)
    {
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      MPI_Barrier(MPI_COMM_WORLD);
      broadcast_cm(&cm, cfg.my_rank, 0);
      if(cfg.my_rank == 0) /* second master only block */
	{
#endif
      /* Allocate and initialize cmstats data structure */
      cmstats[ncm] = AllocCMStats(1); /* Always 1 partition (TEMPORARY) */
      cmstats[ncm]->ps[0] = 0;
      cmstats[ncm]->pe[0] = 100;
      for(i = 0; i < GC_SEGMENTS; i++)
	cmstats[ncm]->gc2p[i] = 0; 
      if(esl_opt_GetBoolean(go, "--time")) esl_stopwatch_Start(w);
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	} /* end of second master only block */
      else /* worker */
	cmstats[ncm] = NULL;
#endif
      if(!(esl_opt_GetBoolean(go, "--lins"))) cm->config_opts |= CM_CONFIG_ZEROINSERTS;
      cm->config_opts |= CM_CONFIG_QDB; /* always use qdb (temporary?) */
      /* Configure the CM */
      ConfigCM(cm, NULL, NULL);

      /****************************************************************
       * Fit Gumbels for all 6 modes (4 w/CM, 2 w/CP9 HMM) 
       *****************************************************************/
      if(!(esl_opt_GetBoolean(go, "--filonly")))
	{
	  for(gum_mode = 0; gum_mode < NGUMBELMODES; gum_mode++)
	    cm_fit_gumbel(cm, go, &cfg, cmstats[ncm], gum_mode);
	  cm->flags |= CM_GUMBEL_STATS;
	}
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      if(cfg.my_rank == 0) /* master */
	{
#endif
      if(esl_opt_GetBoolean(go, "--time"))
	{ 
	  esl_stopwatch_Stop(w);  
	  esl_stopwatch_Display(stdout, w, "Gumbel calculation time:");
	}
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	}
#endif
      /****************************************************************
       * Determine CP9 filtering thresholds
       *****************************************************************/
      if(!(esl_opt_GetBoolean(go, "--gumonly")))
	{
	  if(cfg.my_rank == 0) /* master, MPI or serial */
	    {
	      if(esl_opt_GetBoolean(go, "--filonly")) /* we can only calc filter thr if we had Gumbel stats in input CM file */
		{
		  if(!(cm->flags & CM_GUMBEL_STATS))
		    esl_fatal("ERROR no gumbel stats for CM %d in cmfile: %s, rerun without --filonly.\n", (ncm+1), cfg.cmfile);
		  CopyCMStatsGumbel(cm->stats, cmstats[ncm]);
		}
	      if(esl_opt_GetBoolean(go, "--time")) esl_stopwatch_Start(w);
	    }
	  for(fthr_mode = 0; fthr_mode < NFTHRMODES; fthr_mode++)
	    {
	      if(esl_opt_GetBoolean(go, "--fastfil") && (fthr_mode == CM_LI || fthr_mode == CM_GI)) continue;
	      if((fthr_mode == CM_GC || fthr_mode == CM_GI) || (esl_opt_GetBoolean(go, "--gemit")))
		emit_mode = CM_GC;
	      else
		emit_mode = CM_LC;
	      /* CP9_L, HMM in local mode */
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	      MPI_Barrier(MPI_COMM_WORLD);
#endif
	      l_eval = 
		FindCP9FilterThreshold(cm, cmstats[ncm], cfg.r, 
				       esl_opt_GetReal     (go, "--fract"),
				       esl_opt_GetReal     (go, "--minsurv"),
				       -1.0, /* FIX ME! */
				       -1.0, /* FIX ME! */
				       esl_opt_GetInteger  (go, "--filN"),
				       !(esl_opt_GetBoolean(go, "--nocmeval")),
				       esl_opt_GetReal     (go, "--cmeval"),
				       db_size, emit_mode, fthr_mode, CP9_L, 
				       esl_opt_GetBoolean  (go, "--fastfil"), 
				       FALSE, /* FIX ME! */
				       cfg.my_rank, cfg.nproc, cfg.do_mpi, 
				       NULL, NULL, &l_F);
	      /* CP9_G, HMM in global mode */
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	      MPI_Barrier(MPI_COMM_WORLD);
#endif
	      g_eval = 
		FindCP9FilterThreshold(cm, cmstats[ncm], cfg.r, 
				       esl_opt_GetReal     (go, "--fract"),
				       esl_opt_GetReal     (go, "--minsurv"),
				       -1.0, /* FIX ME! */
				       -1.0, /* FIX ME! */
				       esl_opt_GetInteger  (go, "--filN"),
				       !(esl_opt_GetBoolean(go, "--nocmeval")),
				       esl_opt_GetReal     (go, "--cmeval"),
				       db_size, emit_mode, fthr_mode, CP9_G, 
				       esl_opt_GetBoolean  (go, "--fastfil"), 
				       FALSE, /* FIX ME! */
				       cfg.my_rank, cfg.nproc, cfg.do_mpi, 
				       NULL, NULL, &g_F);
	      if(cfg.my_rank == 0)
		{
		  /* If master (MPI or serial), fill in the filter thr stats */
		  printf("fthr_mode; %d l_eval: %f l_F: %f g_eval: %f g_F: %f\n\n\n", fthr_mode, l_eval, l_F, g_eval, g_F);
		  cmstats[ncm]->fthrA[fthr_mode]->l_eval   = l_eval;
		  cmstats[ncm]->fthrA[fthr_mode]->l_F      = l_F;
		  cmstats[ncm]->fthrA[fthr_mode]->g_eval   = g_eval;
		  cmstats[ncm]->fthrA[fthr_mode]->g_F      = g_F;
		  cmstats[ncm]->fthrA[fthr_mode]->N        = esl_opt_GetInteger  (go, "--filN");
		  cmstats[ncm]->fthrA[fthr_mode]->cm_eval  = esl_opt_GetReal     (go, "--cmeval");
		  cmstats[ncm]->fthrA[fthr_mode]->db_size  = db_size;
		  cmstats[ncm]->fthrA[fthr_mode]->was_fast = esl_opt_GetBoolean  (go, "--fastfil");
		}
	    }
	  if(cfg.my_rank == 0) /* master, MPI or serial */
	    {
	      if(esl_opt_GetBoolean  (go, "--fastfil"))
		{
		  /* we skipped CM_LI and CM_GI, copy CM_LC to CM_LI and CM_GC to CM_GI */
		  CopyFThrInfo(cmstats[ncm]->fthrA[CM_LC], cmstats[ncm]->fthrA[CM_LI]);
		  CopyFThrInfo(cmstats[ncm]->fthrA[CM_GC], cmstats[ncm]->fthrA[CM_GI]);
		}
	      if(esl_opt_GetBoolean(go, "--time"))
		{ 
		  esl_stopwatch_Stop(w);  
		  if(esl_opt_GetBoolean  (go, "--fastfil"))
		    esl_stopwatch_Display(stdout, w, "Fast HMM threshold calculation time:");
		  else
		    esl_stopwatch_Display(stdout, w, "Slow HMM threshold calculation time:");
		}
	      cm->flags |= CM_FTHR_STATS; /* raise flag saying we have filter threshold stats */
	    }
	}
      FreeCM(cm);
      if(cfg.my_rank == 0) /* master, MPI or serial */
	{
	  debug_print_cmstats(cmstats[ncm], (!esl_opt_GetBoolean(go, "--gumonly")));
	  /* Reallocation, if needed. */
	  ncm++;
	  if (ncm == cmalloc) 
	    {
	      cmalloc *= 2;		/* realloc by doubling */
	      ESL_RALLOC(cmstats, tmp, sizeof(CMStats_t *) * cmalloc);
	    }
	  if(!(CMFileRead(cfg.cmfp, NULL, &cm))) continue_flag = 0;
	}
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast (&continue_flag,  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    } /* end of while(continue_flag */

  if(cfg.my_rank == 0) /* master, MPI or serial */
    {
      /*****************************************************************
       * Rewind the CM file for a second pass.
       * Write a temporary CM file with new stats information in it
       *****************************************************************/
      CMFileRewind(cfg.cmfp);
      if (esl_FileExists(tmpfile))
	esl_fatal("Ouch. Temporary file %s appeared during the run.", tmpfile);
      if ((outfp = fopen(tmpfile, mode)) == NULL)
	esl_fatal("Ouch. Temporary file %s couldn't be opened for writing.", tmpfile); 
      
      for (idx = 0; idx < ncm; idx++)
	{
	  /* Sanity checks 
	   */
	  if (!CMFileRead(cfg.cmfp, NULL, &cm))
	    esl_fatal("Ran out of CMs too early in pass 2");
	  if (cm == NULL) 
	    esl_fatal("CM file %s was corrupted? Parse failed in pass 2", cfg.cmfile);
	  
	  cm->stats = cmstats[idx];
	  if(!(esl_opt_GetBoolean(go, "--filonly"))) cm->flags |= CM_GUMBEL_STATS; 
	  if(!(esl_opt_GetBoolean(go, "--gumonly"))) cm->flags |= CM_FTHR_STATS; 
	  
	  /* Save CM to tmpfile
	   */
	  CMFileWrite(outfp, cm, cfg.cmfp->is_binary);
	  FreeCM(cm);
	} /* end of from idx = 0 to ncm */
      
      /*****************************************************************
       * Now, carefully remove original file and replace it
       * with the tmpfile. Note the protection from signals;
       * we wouldn't want a user to ctrl-C just as we've deleted
       * their CM file but before the new one is moved.
       *****************************************************************/
      
      CMFileClose(cfg.cmfp);
      if (fclose(outfp)   != 0)                            esl_fatal("Unexpected Std C/POSIX error.");
      if (sigemptyset(&blocksigs) != 0)                    esl_fatal("Unexpected Std C/POSIX error.");
      if (sigaddset(&blocksigs, SIGINT) != 0)              esl_fatal("Unexpected Std C/POSIX error.");
      if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   esl_fatal("Unexpected Std C/POSIX error.");
      if (remove(cfg.cmfile) != 0)                         esl_fatal("Unexpected Std C/POSIX error.");
      if (rename(tmpfile, cfg.cmfile) != 0)                esl_fatal("Unexpected Std C/POSIX error.");
      if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) esl_fatal("Unexpected Std C/POSIX error.");
      free(tmpfile);

      /***********************************************
       * Exit
       ***********************************************/
    } /* if(cfg.my_rank == 0)  master, MPI or serial */
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      MPI_Finalize();
#endif
  return eslOK;

 ERROR:
  return status;
}

/*
 * Function: cm_fit_gumbel()
 * Date:     EPN, Wed May  9 13:23:37 2007
 * Purpose:  Fits a gumbel distribution from a histogram of scores against 
 *           random sequences.  Fills the relevant Gumbel information in the 
 *           CMStats_t data structure (cmstats) w/mu, and lambda. Determines 
 *           desired search strategy and local/glocal CM/CP9 configuration 
 *           from gum_mode. One histogram is made for each partition of GC frequency. 
 *
 *           This function should be called 6 times, once with each gum_mode
 *           below, to fill all the Gumbel stats in the cmstats data structure:
 *
 *           0. CM_LC: !(cm->search_opts & CM_SEARCH_INSIDE)  w/  local CM
 *           1. CM_GC: !(cm->search_opts & CM_SEARCH_INSIDE)  w/ glocal CM
 *           2. CM_LI:   cm->search_opts & CM_SEARCH_INSIDE   w/  local CM
 *           3. CM_GI:   cm->search_opts & CM_SEARCH_INSIDE   w/ glocal CM
 *           4. CP9_L:   cm->search_opts & CM_SEARCH_HMMONLY  w/  local CP9 HMM
 *           5. CP9_G:   cm->search_opts & CM_SEARCH_HMMONLY  w/ glocal CP9 HMM
 *
 *           This function can be called in MPI mode or not, the flag
 *           is cfg->do_mpi.
 *
 * Args:
 *           cm       - the model
 *           go       - ESL_GETOPTS structure from main
 *           cfg      - the configuration sent from main
 *           cmstats  - data structure that will hold Gumbel stats (NULL if my_rank > 0)
 *           gum_mode - gives search strategy to use for CM or CP9
 *
 */  
static int cm_fit_gumbel(CM_t *cm, ESL_GETOPTS *go, struct cfg_s *cfg, 
			 CMStats_t *cmstats, int gum_mode)
{
  int status = 0;

  /* Configure the CM based on the stat mode */
  ConfigForGumbelMode(cm, gum_mode);
  
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  if(cfg->my_rank > 0)
    {
      mpi_worker_search_target(cm, cfg->my_rank);
      return eslOK;
    }
  else /* my_rank == 0, master (MPI or serial) */
#endif
    {
      int N;                  /* number of random seq samples to use                */
      ESL_HISTOGRAM *h;       /* histogram of scores                                */
      float   sc;             /* score returned from worker                         */
      int     i;              /* counter over sampled seqs                          */
      char   *randseq;        /* a random sequence, textized                        */
      ESL_SQ *sq;             /* random sequence, digitized, to search              */
      int     L;              /* length of random sequences (cm->W*2.)              */
      double *nt_p;	      /* nt distribution for random sequences               */
      double *part_gc_freq;   /* gc content distro to choose from for cur partition */
      double  gc_comp;        /* current gc content                                 */
      double *xv;             /* raw data from histogram                            */
      GumbelInfo_t **gum;     /* for convenience, ptr to relevant cmstats->gumAA[][]*/
      int     p;              /* counter over partitions                            */
      int     n,z;            /* for fitting Gumbels                                */
      double  params[2];      /* used if printing Q-Q file                          */
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      MPI_Status mstatus;     
      int     have_work;      /* for MPI control, do we have work left?             */
      int     nproc_working;  /* number worker processors currently working         */
      int     wi;             /* worker index                                       */
      ESL_SQ **dsqlist = NULL; /* queue of digitized seqs being worked on, 1..nproc-1*/
      ESL_ALLOC(dsqlist,      sizeof(char *) * cfg->nproc);
      for (wi = 0; wi < cfg->nproc; wi++) dsqlist[wi] = NULL;
#endif

      if(gum_mode == CP9_G || gum_mode == CP9_L) N = esl_opt_GetInteger(go, "--hmmN");
      else                                       N = esl_opt_GetInteger(go, "--cmN");

      printf("in cm_fit_gumbel, cfg.my_rank: %d gum_mode: %d\n", cfg->my_rank, gum_mode);
      fflush(stdout);

      /* Check contract */
      if(cmstats == NULL)  esl_fatal("ERROR in cm_fit_gumbel(), master's cmstats is NULL\n");
      if(cfg->r == NULL)   esl_fatal("ERROR in cm_fit_gumbel(), master's source of randomness is NULL\n");
      /* Allocate for random distribution */
      L   = cm->W * 2.;
      ESL_ALLOC(nt_p,         sizeof(double) * cm->abc->K);
      ESL_ALLOC(randseq,      sizeof(char)   * (L+1));
      ESL_ALLOC(part_gc_freq, sizeof(double) * GC_SEGMENTS);
      
      gum = cmstats->gumAA[gum_mode];
      /* The main loop, for each partition... */
      for (p = 0; p < cmstats->np; p++)
	{
	  /* Initialize histogram; these numbers are guesses */
	  h = esl_histogram_CreateFull(-100., 100., .25);    
	      
	  /* If we read in a GC content distro for a file (in this case cfg->gc_freq will not be NULL),
	   * then set up part_gc_freq for this partition */
	  esl_vec_DSet(part_gc_freq, GC_SEGMENTS, 0.);
	  if(cfg->gc_freq != NULL)
	    {
	      for (i = cmstats->ps[p]; i < cmstats->pe[p]; i++) 
		part_gc_freq[i] = cfg->gc_freq[i];
	    }
	  esl_vec_DNorm(part_gc_freq, GC_SEGMENTS);
	  /****************SERIAL BLOCK******************************/
	  if(!(cfg->do_mpi))
	    {
	      /* Take N samples */
	      for (i = 0; i < N; i++) 
		{
		  /* Get random GC content from GC distro (if --dbfile enabled) or 
		  * from CM null model. */
		  if(cfg->gc_freq != NULL)
		    {
		      gc_comp = 0.01 * esl_rnd_DChoose(cfg->r, part_gc_freq, GC_SEGMENTS);
		      nt_p[1] = nt_p[2] = 0.5 * gc_comp;
		      nt_p[0] = nt_p[3] = 0.5 * (1. - gc_comp);
		    }
		  else
		    {
		      nt_p[0] = cm->null[0];
		      nt_p[1] = cm->null[1];
		      nt_p[2] = cm->null[2];
		      nt_p[3] = cm->null[3];
		    }
		  esl_vec_DNorm(nt_p, Alphabet_size);
		  /* Get random sequence */
		  /* We have to generate a text sequence for now, b/c the digitized
		   * version els_rnd_xIID() generates a ESL_DSQ seq, which actually_search_target()
		   * can't handle.
		   */
		  if (esl_rnd_IID(cfg->r, cm->abc->sym, nt_p, cm->abc->K, L, randseq)  != eslOK) 
		    esl_fatal("Failure creating random sequence.");
		  sq = esl_sq_CreateFrom("randseq", randseq, NULL, NULL, NULL, NULL);
		  esl_sq_Digitize(cm->abc, sq);
		  /* Do the scan */
		  sc = actually_search_target(cm, sq, 1, L,
					      0.,    /* cutoff is 0 bits (actually we'll find highest
						      * negative score if it's < 0.0) */
					      0.,    /* CP9 cutoff is 0 bits */
					      NULL,  /* don't keep results */
					      FALSE, /* don't filter with a CP9 HMM */
					      (!(cm->search_opts & CM_SEARCH_HMMONLY)), /* TRUE if we're calcing CM  stats */
					      (cm->search_opts & CM_SEARCH_HMMONLY),    /* TRUE if we're calcing CP9 stats */
					      NULL);          /* filter fraction N/A */
		  esl_sq_Destroy(sq);
		  /* Add best score to histogram */
		  esl_histogram_Add(h, sc);
		}
	    }
	  /*************END OF SERIAL BLOCK*************/
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	  /*************MPI BLOCK***********************/
	  else 
	    {
	      /* Sean's design pattern for data parallelization in a master/worker model:
	       * three phases: 
	       *  1. load workers;
	       *  2. recv result/send work loop;
	       *  3. collect remaining results
	       * but implemented in a single while loop to avoid redundancy.
	       */
	      have_work     = TRUE;
	      nproc_working = 0;
	      wi            = 1;
	      i             = 0;
	      while (have_work || nproc_working)
		{
		  /* Get next work unit. */
		  if (i < N)
		    {
		      /* Get random GC content from GC distro (if --dbfile enabled) or 
		       * from CM null model. */
		      if(cfg->gc_freq != NULL)
			{
			  gc_comp = 0.01 * esl_rnd_DChoose(cfg->r, part_gc_freq, GC_SEGMENTS);
			  nt_p[1] = nt_p[2] = 0.5 * gc_comp;
			  nt_p[0] = nt_p[3] = 0.5 * (1. - gc_comp);
			}
		      else
			{
			  nt_p[0] = cm->bg->f[0];
			  nt_p[1] = cm->bg->f[1];
			  nt_p[2] = cm->bg->f[2];
			  nt_p[3] = cm->bg->f[3];
			}
		      esl_vec_DNorm(nt_p, cm->abc->K);
		      
		      /* Get random sequence */
		      /* We have to generate a text sequence for now, b/c the digitized
		       * version esl_rnd_xIID() generates a ESL_DSQ seq, which actually_search_target()
		       * can't handle. */
		      
		      if (esl_rnd_IID(cfg->r, cm->abc->sym, nt_p, cm->abc->K, L, randseq)  != eslOK) 
			esl_fatal("Failure creating random seq for GC content distro");
		      sq = esl_sq_CreateFrom("randseq", randseq, NULL, NULL, NULL, NULL);
		      esl_sq_Digitize(cm->abc, sq);
		      i++;
		    }
		  else have_work = FALSE;
		  /* If we have work but no free workers, or we have no work but workers
		   * are still working, then wait for a result to return from any worker.
		   */
		  if ( (have_work && nproc_working == cfg->nproc-1) || (! have_work && nproc_working > 0))
		    {
		      MPI_Recv(&sc, 1, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mstatus);
		      wi = mstatus.MPI_SOURCE;
		      if(esl_histogram_Add(h, (double) sc) != eslOK)
			ESL_EXCEPTION(eslENOHALT, "Not eslOK from esl_histogram_Add().");
		      nproc_working--;
		      free(dsqlist[wi]);
		      dsqlist[wi] = NULL;
		    }
		  /* If we have work, assign it to a free worker;
		   * else, terminate the free worker.
		   */
		  if (have_work) 
		    {
		      dsq_MPISend(dsq, L, wi);
		      dsqlist[wi] = dsq;
		      wi++;
		      nproc_working++;
		    }
		  else dsq_MPISend(NULL, -1, wi);	
		} 
	    }
	  /********************END OF MPI BLOCK*************************/
#endif /* if defined(USE_MPI) && defined(MPI_EXECUTABLE) */

	  /* We have all the scores for this partition, fit them to a Gumbel */
	  printf("0 mu: %f lambda: %f\n", gum[p]->mu, gum[p]->lambda);
	  esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
	  esl_gumbel_FitCensored(xv, n, z, xv[0], &(gum[p]->mu), &(gum[p]->lambda));
	  printf("1 mu: %f lambda: %f\n", gum[p]->mu, gum[p]->lambda);
	  fflush(stdout);
	  gum[p]->N = N;
	  gum[p]->L = L;
	  if(cfg->gum_hfp != NULL)
	    esl_histogram_Plot(cfg->gum_hfp, h);
	  if(cfg->gum_qfp != NULL)
	    {
	      params[0] = gum[p]->mu;
	      params[1] = gum[p]->lambda;
	      esl_histogram_PlotQQ(cfg->gum_qfp, h, &esl_exp_generic_invcdf, params);
	    }
	  esl_histogram_Destroy(h);
	}
      free(part_gc_freq);
      free(nt_p);
      free(randseq);
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      free(dsqlist);
#endif
      return eslOK;
    }
 ERROR:
  return status;
}


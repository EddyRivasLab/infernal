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
#include "config.h"	
#include "squidconf.h"
#include "esl_config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <time.h>

#include "structs.h"
#include "funcs.h"		/* external functions                   */
#include "stats.h"              /* GUM functions */
#include "cm_dispatch.h"	
#include "mpifuncs.h"	
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

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-s",        eslARG_INT,      "0", NULL, "n>=0",    NULL,      NULL,        NULL, "set random number seed to <n>",   1 },
  { "--learninserts",eslARG_NONE,FALSE,NULL, NULL,      NULL,      NULL,        NULL, "DO NOT zero insert emission scores",  1},
  { "--dbfile",  eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--filonly", "use GC content distribution from file <s>",  1},
  { "--time",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print timings for Gumbel fitting and CP9 filter calculation",  1},
  { "--gumonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--filonly", "only estimate GUMs, don't calculate filter thresholds", 2},
  { "--cmN",     eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CM GUM estimation",    2 },
  { "--hmmN",    eslARG_INT,   "5000", NULL, "n>0",     NULL,      NULL, "--filonly", "number of random sequences for CP9 HMM GUM estimation",    2 },
  { "--gumhfile",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save fitted Gumbel histogram(s) to file <s>", 2 },
  { "--gumqqfile",eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL, "--filonly", "save Q-Q plot for Gumbel histogram(s) to file <s>", 2 },
  { "--filN",    eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--gumonly", "number of emitted sequences for HMM filter threshold calc",    3 },
  { "--cmeval",  eslARG_REAL,    "10", NULL, "x>0",     NULL,      NULL, "--gumonly", "min CM E-val (for a 1MB db) to consider for filter thr calc", 3}, 
  { "--nocmeval",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--cmeval", "accept all CM hits for filter calc, DO NOT use E-val cutoff", 3}, 
  { "--fract",   eslARG_REAL,  "0.95", NULL, "0<x<=1",  NULL,      NULL, "--gumonly",  "required fraction of CM hits that survive HMM filter", 3},
  { "--filonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "only calculate filter thresholds, don't estimate GUMs", 3},
  { "--fastfil", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "calc filter thr quickly, assume parsetree sc is optimal", 3},
  { "--exp",     eslARG_REAL,  "1.0",   NULL, "x>0.",   NULL,      NULL,        NULL, "exponentiate CM prior to filter threshold calculation", 3},
  { "--gemit",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--gumonly", "when calc'ing filter thresholds, always emit globally from CM",  3},
  { "--filhfile",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL, "--filonly", "save CP9 filter threshold histogram(s) to file <s>", 3},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char           *cmfile;
  CMFILE         *cmfp;
  ESL_RANDOMNESS *r;
  double         *gc_freq;
  FILE           *gum_hfp;
  FILE           *gum_qfp;
  FILE           *fil_hfp;
  int             my_rank;
  int             nproc;
  int             do_mpi;
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit Gumbels for E-value stats and calculate HMM filter threshold stats";

static int cm_fit_gumbel(CM_t *cm, ESL_GETOPTS *go, struct cfg_s *cfg, CMStats_t *cmstats, int gum_mode);

static int CopyFThrInfo(CP9FilterThr_t *src, CP9FilterThr_t *dest);
static int CopyCMStatsGumbel(CMStats_t *src, CMStats_t *dest);

int
main(int argc, char **argv)
{
  int status;			       /* status of a function call                   */
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
  double           params[2];          /* used if printing Q-Q file                   */
  ESL_SQFILE      *dbfp;               /* open file pointer for dbfile                */
  int              format;	       /* format of dbfile                            */
  long             N;                  /* database size if --dbfile enabled           */

  /* CP9 HMM filtering threshold related variables */
  int              fthr_mode;          /* counter over fthr modes                     */
  int              db_size = 1000000;  /* we assume a 1MB database, for use w/cm_ecutoff     */

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
      puts("\ngumbel distribution fitting options :");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nCP9 HMM filter threshold calculation options :");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
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
  cfg.my_rank  = 0;
  cfg.nproc    = 0;
  cfg.do_mpi   = FALSE;

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  /* Initialize MPI, figure out who we are, and whether we're running
   * this show (proc 0) or working in it (procs >0)
   */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &cfg.my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &cfg.nproc);
#endif
  /* Get distribution of GC content from a long random sequence */
  if(esl_opt_GetString(go, "--dbfile") != NULL)
    {
      status = esl_sqfile_Open(esl_opt_GetString(go, "--dbfile"), format, NULL, &dbfp);
      if (status == eslENOTFOUND)    esl_fatal("No such file."); 
      else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
      else if (status == eslEINVAL)  esl_fatal("Can’t autodetect stdin or .gz."); 
      else if (status != eslOK)      esl_fatal("Failed to open sequence database file, code %d.", status); 
      GetDBInfo(dbfp, &N, &cfg.gc_freq);
      printf("Read DB file: %s of length %ld residues (both strands) for GC distro.\n", 
	     esl_opt_GetString(go, "--dbfile"), (2*N));
    }
  else /* use default GC distro */
    {
      /* rmark-1.gc.code from ~/notebook/7_0502_inf_cmcalibrate/gc_distros/
       * Replace with RFAMSEQ derived one */
      cfg.gc_freq = MallocOrDie(sizeof(double) * GC_SEGMENTS);
      cfg.gc_freq[ 0] = 0.; 
      cfg.gc_freq[ 1] = 0.; cfg.gc_freq[ 2] = 0.; cfg.gc_freq[ 3] = 0.; cfg.gc_freq[ 4] = 0.; cfg.gc_freq[ 5] = 0.;
      cfg.gc_freq[ 6] = 1.; cfg.gc_freq[ 7] = 0.; cfg.gc_freq[ 8] = 0.; cfg.gc_freq[ 9] = 0.; cfg.gc_freq[10] = 0.;
      cfg.gc_freq[11] = 1.; cfg.gc_freq[12] = 0.; cfg.gc_freq[13] = 0.; cfg.gc_freq[14] = 0.; cfg.gc_freq[15] = 1.;
      cfg.gc_freq[16] = 3.; cfg.gc_freq[17] = 1.; cfg.gc_freq[18] = 4.; cfg.gc_freq[19] = 0.; cfg.gc_freq[20] = 2.;
      cfg.gc_freq[21] = 2.; cfg.gc_freq[22] = 2.; cfg.gc_freq[23] = 3.; cfg.gc_freq[24] = 2.; cfg.gc_freq[25] = 11.;
      cfg.gc_freq[26] = 6.; cfg.gc_freq[27] = 3.; cfg.gc_freq[28] = 13.; cfg.gc_freq[29] = 9.; cfg.gc_freq[30] = 8.;
      cfg.gc_freq[31] = 8.; cfg.gc_freq[32] = 7.; cfg.gc_freq[33] = 18.; cfg.gc_freq[34] = 10.; cfg.gc_freq[35] = 18.;
      cfg.gc_freq[36] = 27.; cfg.gc_freq[37] = 38.; cfg.gc_freq[38] = 71.; cfg.gc_freq[39] = 86.; cfg.gc_freq[40] = 129.;
      cfg.gc_freq[41] = 198.; cfg.gc_freq[42] = 235.; cfg.gc_freq[43] = 300.; cfg.gc_freq[44] = 362.; cfg.gc_freq[45] = 470.;
      cfg.gc_freq[46] = 557.; cfg.gc_freq[47] = 631.; cfg.gc_freq[48] = 689.; cfg.gc_freq[49] = 754.; cfg.gc_freq[50] = 696.;
      cfg.gc_freq[51] = 756.; cfg.gc_freq[52] = 696.; cfg.gc_freq[53] = 633.; cfg.gc_freq[54] = 560.; cfg.gc_freq[55] = 447.;
      cfg.gc_freq[56] = 403.; cfg.gc_freq[57] = 282.; cfg.gc_freq[58] = 202.; cfg.gc_freq[59] = 182.; cfg.gc_freq[60] = 116.;
      cfg.gc_freq[61] = 96.; cfg.gc_freq[62] = 56.; cfg.gc_freq[63] = 43.; cfg.gc_freq[64] = 27.; cfg.gc_freq[65] = 13.;
      cfg.gc_freq[66] = 14.; cfg.gc_freq[67] = 13.; cfg.gc_freq[68] = 14.; cfg.gc_freq[69] = 9.; cfg.gc_freq[70] = 8.;
      cfg.gc_freq[71] = 12.; cfg.gc_freq[72] = 11.; cfg.gc_freq[73] = 6.; cfg.gc_freq[74] = 4.; cfg.gc_freq[75] = 7.;
      cfg.gc_freq[76] = 3.; cfg.gc_freq[77] = 2.; cfg.gc_freq[78] = 1.; cfg.gc_freq[79] = 4.; cfg.gc_freq[80] = 0.;
      cfg.gc_freq[81] = 2.; cfg.gc_freq[82] = 0.; cfg.gc_freq[83] = 1.; cfg.gc_freq[84] = 0.; cfg.gc_freq[85] = 1.;
      cfg.gc_freq[86] = 0.; cfg.gc_freq[87] = 0.; cfg.gc_freq[88] = 0.; cfg.gc_freq[89] = 0.; cfg.gc_freq[90] = 0.;
      cfg.gc_freq[91] = 0.; cfg.gc_freq[92] = 0.; cfg.gc_freq[93] = 0.; cfg.gc_freq[94] = 0.; cfg.gc_freq[95] = 0.;
      cfg.gc_freq[96] = 0.; cfg.gc_freq[97] = 0.; cfg.gc_freq[98] = 0.; cfg.gc_freq[99] = 0.; cfg.gc_freq[100] = 0.;
      /* distro is centered at just about 50, you can see it, thanks to fixed width font */
    }
  esl_vec_DNorm(cfg.gc_freq, GC_SEGMENTS);

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cmalloc  = 128;
  cmstats  = MallocOrDie(sizeof(CMStats_t *) * cmalloc);
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
    Die("Failed to open covariance model save file %s\n%s\n", cfg.cmfile, usage);

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
  tmpfile = MallocOrDie(strlen(cfg.cmfile) + 5);
  strcpy(tmpfile, cfg.cmfile);
  strcat(tmpfile, ".xxx");	/* could be more inventive here... */
  if (FileExists(tmpfile))
    Die("temporary file %s already exists; please delete it first", tmpfile);
  if (cfg.cmfp->is_binary) mode = "wb";
  else                     mode = "w"; 

  if(esl_opt_GetBoolean(go, "--time"))
    w = esl_stopwatch_Create(); 
  CMFileRead(cfg.cmfp, &cm);
  if(cm == NULL) Die("Failed to read a CM from %s -- file corrupt?\n", cfg.cmfile);

  /* Open output files if nec */
				/* histogram files, fitted and predicted */
  if (esl_opt_GetString(go, "--gumhfile") != NULL) 
    {
      if ((cfg.gum_hfp = fopen(esl_opt_GetString(go, "--gumhfile"), "w")) == NULL)
	Die("Failed to open GUM histogram save file for writing\n");
    }
  if (esl_opt_GetString(go, "--gumqqfile") != NULL) 
    {
      if ((cfg.gum_qfp = fopen(esl_opt_GetString(go, "--gumqqfile"), "w")) == NULL)
	Die("Failed to open GUM histogram save file for writing\n");
    }
  if (esl_opt_GetString(go, "--filhfile") != NULL) 
    {
      if ((cfg.fil_hfp = fopen(esl_opt_GetString(go, "--filhfile"), "w")) == NULL)
	Die("Failed to open filter histogram save file for writing\n");
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
      if(!(esl_opt_GetBoolean(go, "--learninserts"))) cm->config_opts |= CM_CONFIG_ZEROINSERTS;
      cm->config_opts |= CM_CONFIG_QDB; /* always use qdb (temporary?) */
      /* Configure the CM */
      ConfigCM(cm, NULL, NULL);

      /****************************************************************
       * Fit Gumbels for all 6 modes (4 w/CM, 2 w/CP9 HMM) 
       *****************************************************************/
      if(!(esl_opt_GetBoolean(go, "--filonly")))
	{
	  for(gum_mode = 0; gum_mode < NGUMBELMODES; gum_mode++)
	    cm_fit_gumbel(cm, go, &cfg, cmstats[ncm], gum_mode); /* works in MPI and serial mode */
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
	  float fraction;
	  int   filN;
	  int   use_cm_cutoff;
	  float cm_ecutoff;
	  int   do_fastfil;
	  int   emit_global;
	  float exp_factor;

	  if(cfg.my_rank == 0) /* master, MPI or serial */
	    {
	      if(esl_opt_GetBoolean(go, "--filonly")) /* we can only calc filter thr if we had Gumbel stats in input CM file */
		{
		  if(!(cm->flags & CM_GUMBEL_STATS))
		    esl_fatal("ERROR no GUM stats for CM %d in cmfile: %s, rerun without --filonly.\n", (ncm+1), cfg.cmfile);
		  CopyCMStatsGumbel(cm->stats, cmstats[ncm]);
		}
	      if(esl_opt_GetBoolean(go, "--time")) esl_stopwatch_Start(w);
	    }
	  for(fthr_mode = 0; fthr_mode < NFTHRMODES; fthr_mode++)
	    {
/* 	      /\* Skip Inside calculations for  */
	      if(esl_opt_GetBoolean(go, "--fastfil") && (fthr_mode == CM_LI || fthr_mode == CM_GI)) continue;
	      if((fthr_mode == CM_GC || fthr_mode == CM_GI) ||
		 (esl_opt_GetBoolean(go, "--gemit")))
		emit_global = TRUE;
	      else
		emit_global = FALSE;
	      
	      if(fthr_mode == CM_GC || fthr_mode == CM_GI) 
		emit_global = TRUE;
	      else if(fthr_mode == CM_LC || fthr_mode == CM_LI)
		{
		  if(esl_opt_GetBoolean(go, "--gemit")) emit_global = TRUE;
		  else emit_global = FALSE;
		}
	      fraction      =   esl_opt_GetReal   (go, "--fract");
	      filN          =   esl_opt_GetInteger(go, "--filN");
	      use_cm_cutoff = !(esl_opt_GetBoolean(go, "--nocmeval"));
	      cm_ecutoff    =   esl_opt_GetReal   (go, "--cmeval");
	      do_fastfil    =   esl_opt_GetBoolean(go, "--fastfil");
	      emit_global   =   esl_opt_GetBoolean(go, "--gemit");
	      exp_factor    =   esl_opt_GetReal   (go, "--exp");

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	      if(TRUE)
		{
		  /* Two rounds, calc l_eval then g_eval (one for LOCAL CP9 search, one for GLOBAL) */
		  cmstats[ncm]->fthrA[fthr_mode]->l_eval =  
		    mpi_FindCP9FilterThreshold(cm, cmstats[ncm], cfg.r, fraction, filN, use_cm_cutoff, 
					       cm_ecutoff, db_size, emit_global, fthr_mode, 
					       CP9_L, do_fastfil, cfg.my_rank, cfg.nproc);
		  MPI_Barrier(MPI_COMM_WORLD);
		  cmstats[ncm]->fthrA[fthr_mode]->g_eval =  
		    mpi_FindCP9FilterThreshold(cm, cmstats[ncm], cfg.r, fraction, filN, use_cm_cutoff, 
					       cm_ecutoff, db_size, emit_global, fthr_mode, CP9_G, do_fastfil,
					       cfg.my_rank, cfg.nproc);
		}
	      else /* if MPI we want to skip serial calls below */
#endif
		{
		  cmstats[ncm]->fthrA[fthr_mode]->l_eval =  /*  local */
		    serial_FindCP9FilterThreshold(cm, cmstats[ncm], cfg.r, fraction, filN, use_cm_cutoff,
						  cm_ecutoff, db_size, emit_global, fthr_mode, CP9_L, do_fastfil);
		  cmstats[ncm]->fthrA[fthr_mode]->g_eval =  /* glocal */
		    serial_FindCP9FilterThreshold(cm, cmstats[ncm], cfg.r, fraction, filN, use_cm_cutoff,
						  cm_ecutoff, db_size, emit_global, fthr_mode, CP9_G, do_fastfil);
		}
	      if(cfg.my_rank == 0)
		{
		  /* If master (MPI or serial), fill in the rest of the stats */
		  cmstats[ncm]->fthrA[fthr_mode]->N        = filN;
		  cmstats[ncm]->fthrA[fthr_mode]->fraction = fraction;
		  cmstats[ncm]->fthrA[fthr_mode]->cm_eval  = cm_ecutoff;
		  cmstats[ncm]->fthrA[fthr_mode]->db_size  = db_size;
		  cmstats[ncm]->fthrA[fthr_mode]->was_fast = FALSE;
		}
	    }
	  if(cfg.my_rank == 0) /* master, MPI or serial */
	    {
	      if(do_fastfil) 
		{
		  /* we skipped CM_LI and CM_GI, copy CM_LC to CM_LI and CM_GC to CM_GI */
		  CopyFThrInfo(cmstats[ncm]->fthrA[CM_LC], cmstats[ncm]->fthrA[CM_LI]);
		  CopyFThrInfo(cmstats[ncm]->fthrA[CM_GC], cmstats[ncm]->fthrA[CM_GI]);
		}
	      if(esl_opt_GetBoolean(go, "--time"))
		{ 
		  esl_stopwatch_Stop(w);  
		  if(do_fastfil)
		    esl_stopwatch_Display(stdout, w, "Fast HMM threshold calculation time:");
		  else
		    esl_stopwatch_Display(stdout, w, "Slow HMM threshold calculation time:");
		}
	      cm->flags &= CM_FTHR_STATS;
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
	      cmstats = ReallocOrDie(cmstats, sizeof(CMStats_t *) * cmalloc);
	    }
	  if(!(CMFileRead(cfg.cmfp, &cm))) continue_flag = 0;
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
      if (FileExists(tmpfile))
	esl_fatal("Ouch. Temporary file %s appeared during the run.", tmpfile);
      if ((outfp = fopen(tmpfile, mode)) == NULL)
	esl_fatal("Ouch. Temporary file %s couldn't be opened for writing.", tmpfile); 
      
      for (idx = 0; idx < ncm; idx++)
	{
	  /* Sanity checks 
	   */
	  if (!CMFileRead(cfg.cmfp, &cm))
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
      if (fclose(outfp)   != 0)                            PANIC;
      if (sigemptyset(&blocksigs) != 0)                    PANIC;
      if (sigaddset(&blocksigs, SIGINT) != 0)              PANIC;
      if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   PANIC;
      if (remove(cfg.cmfile) != 0)                         PANIC;
      if (rename(tmpfile, cfg.cmfile) != 0)                PANIC;
      if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) PANIC;
      free(tmpfile);
      /***********************************************
       * Exit
       ***********************************************/
    } /* if(cfg.my_rank == 0)  master, MPI or serial */
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      MPI_Finalize();
#endif
  return eslOK;
}

/*
 * Function: cm_fit_gumbel()
 * Date:     EPN, Wed May  9 13:23:37 2007
 * Purpose:  Fits a gumbel distribution from a histogram of scores against 
 *           random sequences.  Fills the relevant GUM information in the 
 *           CMStats_t data structure (cmstats) w/mu, and lambda. Determines 
 *           desired search strategy and local/glocal CM/CP9 configuration 
 *           from gum_mode. One histogram is made for each partition of GC frequency. 
 *
 *           This function should be called 6 times, once with each gum_mode
 *           below, to fill all the GUM stats in the cmstats data structure:
 *
 *           This function can be called in MPI mode or not, the flag
 *           is cfg->do_mpi.
 *
 *           0. CM_LC: !cm->search_opts & CM_SEARCH_INSIDE  w/  local CM
 *           1. CM_GC: !cm->search_opts & CM_SEARCH_INSIDE  w/ glocal CM
 *           2. CM_LI:  cm->search_opts & CM_SEARCH_INSIDE  w/  local CM
 *           3. CM_GI:  cm->search_opts & CM_SEARCH_INSIDE  w/ glocal CM
 *           4. CP9_L:  cm->search_opts & CM_SEARCH_HMMONLY w/  local CP9 HMM
 *           5. CP9_G:  cm->search_opts & CM_SEARCH_HMMONLY w/ glocal CP9 HMM
 *
 * Args:
 *           cm       - the model
 *           go       - ESL_GETOPTS structure from main
 *           cfg      - the configuration sent from main
 *           cmstats  - data structure that will hold GUM stats (NULL if my_rank > 0)
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
      mpi_worker_search_target(cm, my_rank);
      return eslOK;
    }
  else /* my_rank == 0, master */
#endif
    {
      int N;                  /* number of random seq samples to use                */
      ESL_HISTOGRAM *h;       /* histogram of scores                                */
      float   sc;             /* score returned from worker                         */
      int     i;              /* counter over sampled seqs                          */
      char   *randseq;        /* a random sequence, textized                        */
      char   *dsq;            /* randseq, digitized                                 */
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
      char  **dsqlist = NULL; /* queue of digitized seqs being worked on, 1..nproc-1*/
      ESL_ALLOC(dsqlist,      sizeof(char *) * cfg->nproc);
      for (wi = 0; wi < cfg->nproc; wi++) dsqlist[wi] = NULL;
#endif

      if(gum_mode == CP9_G || gum_mode == CP9_L) N = esl_opt_GetInteger(go, "--hmmN");
      else N = esl_opt_GetInteger(go, "--cmN");

      printf("in cm_fit_gumbel, cfg.my_rank: %d gum_mode: %d\n", cfg->my_rank, gum_mode);
      fflush(stdout);

      /* Check contract */
      if(cmstats == NULL) esl_fatal("ERROR in cm_fit_gumbel(), master's cmstats is NULL\n");
      if(cfg->r == NULL)   esl_fatal("ERROR in cm_fit_gumbel(), master's source of randomness is NULL\n");

      /* Allocate for random distribution */
      L   = cm->W * 2.;
      ESL_ALLOC(nt_p,         sizeof(double) * Alphabet_size);
      ESL_ALLOC(randseq,      sizeof(char)   * (L+1));
      ESL_ALLOC(part_gc_freq, sizeof(double) * GC_SEGMENTS);
      
      gum = cmstats->gumAA[gum_mode];
      /* The main loop, for each partition... */
      for (p = 0; p < cmstats->np; p++)
	{
	  /* Initialize histogram; these numbers are guesses */
	  h = esl_histogram_CreateFull(-100., 100., .1);    
	      
	  /* Set up part_gc_freq for this partition */
	  esl_vec_DSet(part_gc_freq, GC_SEGMENTS, 0.);
	  for (i = cmstats->ps[p]; i < cmstats->pe[p]; i++) 
	    part_gc_freq[i] = cfg->gc_freq[i];
	  esl_vec_DNorm(part_gc_freq, GC_SEGMENTS);
	  
	  /****************SERIAL BLOCK******************************/
	  if(!(cfg->do_mpi))
	    {
	      /* Take N samples */
	      for (i = 0; i < N; i++) 
		{
		  /* Get random GC content */
		  gc_comp = 0.01 * esl_rnd_DChoose(cfg->r, part_gc_freq, GC_SEGMENTS);
		  nt_p[1] = nt_p[2] = 0.5 * gc_comp;
		  nt_p[0] = nt_p[3] = 0.5 * (1. - gc_comp);
		  
		  /* Get random sequence */
		  /* We have to generate a text sequence for now, b/c the digitized
		   * version els_rnd_xIID() generates a ESL_DSQ seq, which actually_search_target()
		   * can't handle.
		   */
		  if (esl_rnd_IID(cfg->r, Alphabet, nt_p, Alphabet_size, L, randseq)  != eslOK) 
		    esl_fatal("Failure creating random seq for GC content distro");
		  /* Digitize the sequence, parse it, and add to histogram */
		  dsq = DigitizeSequence (randseq, L);
		  /* Do the scan */
		  sc = actually_search_target(cm, dsq, 1, L,
					      0.,    /* cutoff is 0 bits (actually we'll find highest
						      * negative score if it's < 0.0) */
					      0.,    /* CP9 cutoff is 0 bits */
					      NULL,  /* don't keep results */
					      FALSE, /* don't filter with a CP9 HMM */
					      (!(cm->search_opts & CM_SEARCH_HMMONLY)), /* TRUE if we're calcing CM  stats */
					      (cm->search_opts & CM_SEARCH_HMMONLY),    /* TRUE if we're calcing CP9 stats */
					      NULL);          /* filter fraction N/A */
		  /* Add best score to histogram */
		  esl_histogram_Add(h, sc);
		}
	    }
	  /*************END OF SERIAL BLOCK*************/
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	  /*************MPI BLOCK***********************/
	  else 
	    {
	      /* From Sean's design pattern in hmmsim.c */
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
		      /* Get random GC content */
		      gc_comp = 0.01 * esl_rnd_DChoose(cfg->r, part_gc_freq, GC_SEGMENTS);
		      nt_p[1] = nt_p[2] = 0.5 * gc_comp;
		      nt_p[0] = nt_p[3] = 0.5 * (1. - gc_comp);
		      esl_vec_DNorm(nt_p, Alphabet_size);
		      
		      /* Get random sequence */
		      /* We have to generate a text sequence for now, b/c the digitized
		       * version esl_rnd_xIID() generates a ESL_DSQ seq, which actually_search_target()
		       * can't handle. */
		      
		      if (esl_rnd_IID(r, Alphabet, nt_p, Alphabet_size, L, randseq)  != eslOK) 
			esl_fatal("Failure creating random seq for GC content distro");
		      /* Digitize the sequence, send it off to be searched, receive score and add to histogram */
		      dsq = DigitizeSequence (randseq, L);
		      i++;
		    }
		  else have_work = FALSE;
		  /* If we have work but no free workers, or we have no work but workers
		   * are still working, then wait for a result to return from any worker.
		   */
		  if ( (have_work && nproc_working == nproc-1) || (! have_work && nproc_working > 0))
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
	  esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
	  esl_gumbel_FitCensored(xv, n, z, xv[0], &(gum[p]->mu), &(gum[p]->lambda));
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


/* Function: CopyFThrInfo()
 * Incept:   EPN, Fri May  4 15:54:51 2007
 */
int CopyFThrInfo(CP9FilterThr_t *src, CP9FilterThr_t *dest)
{
  dest->N           = src->N;
  dest->fraction    = src->fraction;
  dest->cm_eval     = src->cm_eval;
  dest->l_eval      = src->l_eval;
  dest->g_eval      = src->g_eval;
  dest->db_size     = src->db_size;
  dest->was_fast    = src->was_fast;
  return eslOK;
}

/* Function: CopyCMStatsGumbel()
 * Incept:   EPN, Mon May  7 06:04:58 2007
 * 
 * Purpose:  Copy the GUM stats in a source CMStats_t object into
 *           a pre-alloc'ed destination CMStats_t object.
 */
int CopyCMStatsGumbel(CMStats_t *src, CMStats_t *dest)
{
  int i, p;

  /* Check contract */
  if(src->np != dest->np)
    Die("ERROR in CopyCMStatsGUM() src->np: %d not equal to alloc'ed dest->np: %d\n", src->np, dest->np);

  for(p = 0; p < src->np; p++)
    {
      dest->ps[p] = src->ps[p];
      dest->pe[p] = src->pe[p];
    }
  for(i = 0; i < GC_SEGMENTS; i++)
    dest->gc2p[i] = src->gc2p[i]; 

  for(i = 0; i < NGUMBELMODES; i++)
    {
      for(p = 0; p < src->np; p++)
	{
	  dest->gumAA[i][p]->N      = src->gumAA[i][p]->N;
	  dest->gumAA[i][p]->L      = src->gumAA[i][p]->L;
	  dest->gumAA[i][p]->mu     = src->gumAA[i][p]->mu;
	  dest->gumAA[i][p]->lambda = src->gumAA[i][p]->lambda;
	}
    }
  return eslOK;
}



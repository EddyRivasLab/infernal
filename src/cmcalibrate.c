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

#include "structs.h"
#include "funcs.h"		/* external functions                   */
#include "stats.h"              /* EVD functions */
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
#include "esl_gumbel.h"
#include "esl_stopwatch.h"

static int NEW_serial_make_histogram (CM_t *cm, CMStats_t *cmstats, double *gc_freq, 
				      int N, int L, int evd_mode);
static int CopyFThrInfo(CP9FilterThr_t *src, CP9FilterThr_t *dest);
static int CopyCMStatsEVD(CMStats_t *src, CMStats_t *dest);

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
static int mpi_make_histogram (CM_t *cm, CMStats_t *cmstats, double *gc_freq, 
			       int N, int L, int evd_mode, int my_rank, int nproc);
#endif

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "--cmevdN",  eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--onlyfil", "number of random sequences for CM EVD estimation",    1 },
  { "--hmmevdN", eslARG_INT,   "5000", NULL, "n>0",     NULL,      NULL, "--onlyfil", "number of random sequences for CP9 HMM EVD estimation",    1 },
  { "--filN",    eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL, "--onlyevd", "number of emitted sequences for HMM filter threshold calc",    1 },
  { "--cmeval",  eslARG_REAL,    "10", NULL, "x>0",     NULL,      NULL, "--onlyevd", "min CM E-val (for a 1MB db) to consider for filter thr calc", 1}, 
  { "--nocmeval",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--cmeval", "accept all CM hits for filter calc, DO NOT use E-val cutoff", 1}, 
  { "--fract",   eslARG_REAL,  "0.95", NULL, "0<x<=1",  NULL,      NULL, "--onlyevd",  "required fraction of CM hits that survive HMM filter", 1},
  { "--onlyevd", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--onlyfil", "only estimate EVDs, don't calculate filter thresholds", 1},
  { "--onlyfil", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--onlyevd", "only calculate filter thresholds, don't estimate EVDs", 1},
  { "--fastfil", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--onlyevd", "calc filter thr quickly, assume parsetree sc is optimal", 1},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
static char usage[]  = "mpi-cmcalibrate [-options]";
#else
static char usage[]  = "cmcalibrate [-options]";
#endif

int
main(int argc, char **argv)
{
  int status;			       /* status of a function call                   */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                     */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                        */
  CMFILE          *cmfp;               /* open CM file for reading                    */
  CM_t            *cm      = NULL;     /* the covariance model                        */
  CMStats_t      **cmstats;            /* the CM stats data structures, 1 for each CM */
  char            *cmfile;             /* file to read CM(s) from                     */
  char            *tmpfile;            /* temporary calibrated CM file                */
  FILE            *outfp;              /* for writing CM(s) into tmpfile              */
  char            *mode;               /* write mode, "w" or "wb"                     */
  int              idx;		       /* counter over CMs                            */
  sigset_t         blocksigs;	       /* list of signals to protect from             */
  int              do_timings = TRUE;  /* TRUE to print timings, FALSE not to         */
  ESL_STOPWATCH   *w;                  /* watch for timing EVD, filter thr calc times */
  int              ncm;                /* number of CMs read from cmfile              */
  int              cmalloc;            /* for alloc'ing CMStats_t objects             */
  int              optset;             /* boolean value for esl_get_opts()            */
  double          *gc_freq;            /* GC frequency array [0..100(GC_SEGMENTS)]    */
  int              i;                  /* counter over GC segments                    */
  int              continue_flag;      /* we still have a CM to get stats for?        */

  /* EVD related variables */
  int              do_qdb = TRUE;      /* TRUE to use QDB while searching with CM     */
  int              do_evd = TRUE;      /* TRUE to calc EVDs, set FALSE if --onlyfil   */
  int              sample_length;      /* sample len used for calc'ing stats (2*W)    */
  int              cmevdN  = 1000;     /* # samples used to calculate  CM EVDs        */
  int              hmmevdN = 5000;     /* # samples used to calculate HMM EVDs        */
  int              evdN;               /* either cmevdN or hmmevdN                    */
  int              evd_mode;           /* counter over evd modes                      */

  /* CP9 HMM filtering threshold related variables */
  int              do_fil = TRUE;      /* TRUE to calc filter thr, FALSE if --onlyevd */
  int              filN    =  500;     /* # emitted seqs to use to calc filter thr    */
  int              do_fastfil = FALSE; /* TRUE: use fast hacky filter thr calc method */
  float            fraction = 0.95;    /* fraction of CM hits req'd to find with HMM  */
  float            cm_ecutoff = 10;    /* minimum CM E-value we care about            */
  int              use_cm_cutoff = TRUE; /* TRUE to use cm_ecutoff, FALSE not to      */
  int              fthr_mode;          /* counter over fthr modes                     */
  int              db_size = 1000000;  /* we assume a 1MB database, for use w/cm_ecutoff     */
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  int my_rank;              /* My rank in MPI                           */
  int nproc;                /* number of processors */
#endif

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    puts("\ngeneral options are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 2 = indentation; 80=textwidth*/
    return eslOK;
  }

  esl_opt_GetIntegerOption(go, "--cmevdN",  &cmevdN);
  esl_opt_GetIntegerOption(go, "--hmmevdN", &hmmevdN);
  esl_opt_GetIntegerOption(go, "--filN",    &filN);
  esl_opt_GetFloatOption  (go, "--cmeval",  &cm_ecutoff);
  esl_opt_GetFloatOption  (go, "--fract",   &fraction);
  esl_opt_GetBooleanOption(go, "--onlyevd", &optset); if(optset) do_fil = FALSE;
  esl_opt_GetBooleanOption(go, "--onlyfil", &optset); if(optset) do_evd = FALSE;
  esl_opt_GetBooleanOption(go, "--fastfil", &do_fastfil);
  esl_opt_GetBooleanOption(go, "--nocmeval",&optset); 
  if(optset) { use_cm_cutoff = FALSE; cm_ecutoff = -1.; }

  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }
  cmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL); /* NULL=no range checking */

  /****************************************************************
   * Get distribution of GC content from a long random sequence
   *****************************************************************/
  /* human_chr1.gc.code from ~/notebook/7_0502_inf_cmcalibrate/gc_distros/ */
  /* Replace with RFAMSEQ derived one */
  gc_freq = MallocOrDie(sizeof(double) * GC_SEGMENTS);
  gc_freq[0] = 79.; 
  gc_freq[1] = 59.; gc_freq[2] = 56.; gc_freq[3] = 79.; gc_freq[4] = 91.; gc_freq[5] = 89.;
  gc_freq[6] = 114.; gc_freq[7] = 105.; gc_freq[8] = 152.; gc_freq[9] = 149.; gc_freq[10] = 181.;
  gc_freq[11] = 203.; gc_freq[12] = 293.; gc_freq[13] = 401.; gc_freq[14] = 586.; gc_freq[15] = 894.;
  gc_freq[16] = 1256.; gc_freq[17] = 1909.; gc_freq[18] = 2901.; gc_freq[19] = 4333.; gc_freq[20] = 6421.;
  gc_freq[21] = 9027.; gc_freq[22] = 12494.; gc_freq[23] = 16458.; gc_freq[24] = 20906.; gc_freq[25] = 26706.;
  gc_freq[26] = 32624.; gc_freq[27] = 38814.; gc_freq[28] = 44862.; gc_freq[29] = 51035.; gc_freq[30] = 57425.;
  gc_freq[31] = 62700.; gc_freq[32] = 67732.; gc_freq[33] = 72385.; gc_freq[34] = 75656.; gc_freq[35] = 78714.;
  gc_freq[36] = 80693.; gc_freq[37] = 81741.; gc_freq[38] = 81731.; gc_freq[39] = 81291.; gc_freq[40] = 79600.;
  gc_freq[41] = 76883.; gc_freq[42] = 75118.; gc_freq[43] = 71834.; gc_freq[44] = 69564.; gc_freq[45] = 68058.;
  gc_freq[46] = 66882.; gc_freq[47] = 64922.; gc_freq[48] = 63846.; gc_freq[49] = 61481.; gc_freq[50] = 57994.;
  gc_freq[51] = 54322.; gc_freq[52] = 50144.; gc_freq[53] = 46483.; gc_freq[54] = 43104.; gc_freq[55] = 40144.;
  gc_freq[56] = 36987.; gc_freq[57] = 33653.; gc_freq[58] = 29844.; gc_freq[59] = 25838.; gc_freq[60] = 21861.;
  gc_freq[61] = 17858.; gc_freq[62] = 14685.; gc_freq[63] = 12233.; gc_freq[64] = 9981.; gc_freq[65] = 8047.;
  gc_freq[66] = 6495.; gc_freq[67] = 5263.; gc_freq[68] = 4210.; gc_freq[69] = 3405.; gc_freq[70] = 2653.;
  gc_freq[71] = 2157.; gc_freq[72] = 1786.; gc_freq[73] = 1458.; gc_freq[74] = 1230.; gc_freq[75] = 1048.;
  gc_freq[76] = 915.; gc_freq[77] = 802.; gc_freq[78] = 769.; gc_freq[79] = 668.; gc_freq[80] = 544.;
  gc_freq[81] = 462.; gc_freq[82] = 395.; gc_freq[83] = 329.; gc_freq[84] = 228.; gc_freq[85] = 170.;
  gc_freq[86] = 116.; gc_freq[87] = 94.; gc_freq[88] = 39.; gc_freq[89] = 46.; gc_freq[90] = 22.;
  gc_freq[91] = 11.; gc_freq[92] = 5.; gc_freq[93] = 3.; gc_freq[94] = 0.; gc_freq[95] = 1.;
  gc_freq[96] = 0.; gc_freq[97] = 1.; gc_freq[98] = 0.; gc_freq[99] = 0.; gc_freq[100] = 0.;
  /* distro is skewed towards low GC, you can see it, thanks to fixed width font */
  esl_vec_DNorm(gc_freq, GC_SEGMENTS);

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  /* Initialize MPI, figure out who we are, and whether we're running
   * this show (proc 0) or working in it (procs >0)
   */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if(my_rank == 0) /* master only block 1 */
    {
#endif
  /*****************************************************************
   * Initializations, including opening the CM file 
   *****************************************************************/
  /*if ((r = esl_randomness_CreateTimeseeded()) == NULL)*/
  if ((r = esl_randomness_Create(33)) == NULL)
    esl_fatal("Failed to create random number generator: probably out of memory");

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cmalloc  = 128;
  cmstats  = MallocOrDie(sizeof(CMStats_t *) * cmalloc);
  ncm      = 0;
  
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);

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
  tmpfile = MallocOrDie(strlen(cmfile) + 5);
  strcpy(tmpfile, cmfile);
  strcat(tmpfile, ".xxx");	/* could be more inventive here... */
  if (FileExists(tmpfile))
    Die("temporary file %s already exists; please delete it first", tmpfile);
  if (cmfp->is_binary) mode = "wb";
  else                 mode = "w"; 

  if(do_timings)
    w = esl_stopwatch_Create(); 
  CMFileRead(cmfp, &cm);
  if(cm == NULL) Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
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
      broadcast_cm(&cm, my_rank, 0);
      if(my_rank == 0) /* second master only block */
	{
#endif
      /* Allocate and initialize cmstats data structure */
      cmstats[ncm] = AllocCMStats(1); /* Always 1 partition (TEMPORARY) */
      cmstats[ncm]->ps[0] = 0;
      cmstats[ncm]->pe[0] = 100;
      for(i = 0; i < GC_SEGMENTS; i++)
	cmstats[ncm]->gc2p[i] = 0; 
      if(do_timings) esl_stopwatch_Start(w);
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	} /* end of second master only block */
#endif
      /* Configure for QDB */
      if(do_qdb) cm->config_opts |= CM_CONFIG_QDB;
      ConfigCM(cm, NULL, NULL);
      /****************************************************************
       * Determine CM and CP9 EVDs 
       *****************************************************************/
      if(do_evd)
	{
	  sample_length = 2.0 * cm->W;
	  for(evd_mode = 0; evd_mode < NEVDMODES; evd_mode++)
	    {
	      if(evd_mode == CP9_G || evd_mode == CP9_L) evdN = hmmevdN;
	      else evdN = cmevdN;

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	      MPI_Barrier(MPI_COMM_WORLD);
	      if(my_rank == 0)
		mpi_make_histogram (cm, cmstats[ncm], gc_freq, evdN, sample_length, evd_mode, my_rank, nproc);
	      else if (my_rank > 0)
		mpi_make_histogram (cm, NULL, NULL, 0, 0, evd_mode, my_rank, nproc);
	      else /* we'll never get here */
#endif
	      NEW_serial_make_histogram (cm, cmstats[ncm], gc_freq, evdN, sample_length, evd_mode);
	    }
	  cm->flags |= CM_EVD_STATS;
	}
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      if(my_rank == 0) /* master */
	{
#endif
      if(do_timings) 
	{ 
	  esl_stopwatch_Stop(w);  
	  esl_stopwatch_Display(stdout, w, "EVD calculation time:");
	}
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	}
#endif
      /****************************************************************
       * Determine CP9 filtering thresholds
       *****************************************************************/
      if(do_fil)
	{
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	  if(my_rank == 0) /* master */
	    {
#endif
	  if(!do_evd) /* we can only calc filter thresholds if we had EVD stats in input CM file */
	    {
	      if(!(cm->flags & CM_EVD_STATS))
		Die("ERROR no EVD stats for CM %d in cmfile: %s, rerun without --onlyfil.\n", (ncm+1), cmfile);
	      CopyCMStatsEVD(cm->stats, cmstats[ncm]);
	    }
	  if(do_timings) esl_stopwatch_Start(w);
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	    }
#endif
	  for(fthr_mode = 0; fthr_mode < NFTHRMODES; fthr_mode++)
	    {
	      if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI)) continue;
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	      if(TRUE)
		{
		  MPI_Barrier(MPI_COMM_WORLD);
		  if(my_rank == 0) /* master */
		    {
		      cmstats[ncm]->fthrA[fthr_mode]->l_eval =  
			mpi_FindCP9FilterThreshold(cm, cmstats[ncm], fraction, filN, use_cm_cutoff, 
						   cm_ecutoff, db_size, fthr_mode, CP9_L, do_fastfil,
						   my_rank, nproc);
		    }
		  else /* worker */
		    mpi_FindCP9FilterThreshold(cm, NULL, fraction, filN, use_cm_cutoff, 
					     cm_ecutoff, db_size, fthr_mode, CP9_L, do_fastfil,
					     my_rank, nproc);

		    
		  MPI_Barrier(MPI_COMM_WORLD);
		  if(my_rank == 0) /* master */
		    {
		      cmstats[ncm]->fthrA[fthr_mode]->g_eval =  
		      mpi_FindCP9FilterThreshold(cm, cmstats[ncm], fraction, filN, use_cm_cutoff, 
						 cm_ecutoff, db_size, fthr_mode, CP9_G, do_fastfil,
						 my_rank, nproc);
		    }
		  else
		    mpi_FindCP9FilterThreshold(cm, NULL, fraction, filN, use_cm_cutoff, 
					       cm_ecutoff, db_size, fthr_mode, CP9_G, do_fastfil,
					       my_rank, nproc);
		}
	      else /* if MPI we want to skip serial calls below */
		{
#endif
		cmstats[ncm]->fthrA[fthr_mode]->l_eval =  /*  local */
		  FindCP9FilterThreshold(cm, cmstats[ncm],fraction, filN, use_cm_cutoff,
					 cm_ecutoff, db_size, fthr_mode, CP9_L, do_fastfil);
		cmstats[ncm]->fthrA[fthr_mode]->g_eval =  /* glocal */
		  FindCP9FilterThreshold(cm, cmstats[ncm],fraction, filN, use_cm_cutoff,
					 cm_ecutoff, db_size, fthr_mode, CP9_G, do_fastfil);
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
		}
	      if(my_rank == 0)
		{
#endif
	      /* Fill in the rest of the stats */
	      cmstats[ncm]->fthrA[fthr_mode]->N        = filN;
	      cmstats[ncm]->fthrA[fthr_mode]->fraction = fraction;
	      cmstats[ncm]->fthrA[fthr_mode]->cm_eval  = cm_ecutoff;
	      cmstats[ncm]->fthrA[fthr_mode]->db_size  = db_size;
	      cmstats[ncm]->fthrA[fthr_mode]->was_fast = FALSE;
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
		}
#endif
	    }

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	  if(my_rank == 0) /* master */
	    {
#endif
	  if(do_fastfil) 
	    {
	      /* we skipped CM_LI and CM_GI, copy CM_LC to CM_LI and CM_GC to CM_GI */
	      CopyFThrInfo(cmstats[ncm]->fthrA[CM_LC], cmstats[ncm]->fthrA[CM_LI]);
	      CopyFThrInfo(cmstats[ncm]->fthrA[CM_GC], cmstats[ncm]->fthrA[CM_GI]);
	    }
	  if(do_timings) 
	    { 
	      esl_stopwatch_Stop(w);  
	      esl_stopwatch_Display(stdout, w, "Slow HMM threshold (1/4 modes) calculation time:");
	    }
	  cm->flags &= CM_FTHR_STATS;
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	    }
#endif
	}
      FreeCM(cm);
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
      if(my_rank == 0) /* master */
	{
#endif
      debug_print_cmstats(cmstats[ncm], do_fil);
      /* Reallocation, if needed.
       */
      ncm++;
      if (ncm == cmalloc) 
	{
	  cmalloc *= 2;		/* realloc by doubling */
	  cmstats = ReallocOrDie(cmstats, sizeof(CMStats_t *) * cmalloc);
	}
      if(!(CMFileRead(cmfp, &cm))) continue_flag = 0;
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
	}
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast (&continue_flag,  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    } /* end of while(continue_flag */

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  if(my_rank > 0) /* job's done worker */
    {
      MPI_Finalize();
      return 0;
    }
#endif
  esl_stopwatch_Destroy(w);

  /*****************************************************************
   * Rewind the CM file for a second pass.
   * Write a temporary CM file with new stats information in it
   *****************************************************************/
  CMFileRewind(cmfp);
  if (FileExists(tmpfile))
    Die("Ouch. Temporary file %s appeared during the run.", tmpfile);
  if ((outfp = fopen(tmpfile, mode)) == NULL)
    Die("Ouch. Temporary file %s couldn't be opened for writing.", tmpfile); 
  
  for (idx = 0; idx < ncm; idx++)
  {
  /* Sanity checks 
   */
    if (!CMFileRead(cmfp, &cm))
      Die("Ran out of CMs too early in pass 2");
    if (cm == NULL) 
      Die("CM file %s was corrupted? Parse failed in pass 2", cmfile);
    
    cm->stats = cmstats[idx];
    if(do_evd) cm->flags |= CM_EVD_STATS; 
    if(do_fil) cm->flags |= CM_FTHR_STATS; 

    /* Save CM to tmpfile
     */
    CMFileWrite(outfp, cm, cmfp->is_binary);
    FreeCM(cm);
  } /* end of from idx = 0 to ncm */
  
  /*****************************************************************
   * Now, carefully remove original file and replace it
   * with the tmpfile. Note the protection from signals;
   * we wouldn't want a user to ctrl-C just as we've deleted
   * their CM file but before the new one is moved.
   *****************************************************************/
  
  CMFileClose(cmfp);
  if (fclose(outfp)   != 0)                            PANIC;
  if (sigemptyset(&blocksigs) != 0)                    PANIC;
  if (sigaddset(&blocksigs, SIGINT) != 0)              PANIC;
  if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   PANIC;
  if (remove(cmfile) != 0)                             PANIC;
  if (rename(tmpfile, cmfile) != 0)                    PANIC;
  if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) PANIC;
  free(tmpfile);
  /***********************************************
   * Exit
   ***********************************************/
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  if(my_rank == 0) /* master */
    MPI_Finalize();
#endif
  return eslOK;
}

/*
 * Function: NEW_serial_make_histogram()
 * Date:     Mon Apr 1 2002 [St. Louis]
 * Purpose:  Makes histogram(s) using random sequences.  Fills the relevant
 *           EVD information in the CMStats_t data structure (cmstats) w/mu,
 *           K, and lambda. Determines desired search strategy and 
 *           local/glocal CM/CP9 configuration from evd_mode.
 *           One histogram is made for each partition of GC frequency. 
 *
 *           This function should be called 6 times, once with each evd_mode
 *           below, to fill all the EVD stats in the cmstats data structure:
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
 *           cmstats  - data structure that will hold EVD stats
 *           gc_freq  - GC frequencies for each GC content
 *           N        - number of random seqs to search 
 *           L        - length of random samples
 *           evd_mode - gives search strategy to use for CM or CP9
 *
 */  
static int NEW_serial_make_histogram (CM_t *cm, CMStats_t *cmstats, double *gc_freq, 
				      int N, int L, int evd_mode)
{
  int i;
  char *randseq;
  char *dsq;
  ESL_HISTOGRAM *h;
  double  *nt_p;		     /* Distribution for random sequences */
  int status = 0;
  float score;
  double cur_gc_freq[GC_SEGMENTS];
  float gc_comp;
  double *xv;
  int z;
  int n;
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  EVDInfo_t **evd; 
  int p; /* counter over partitions */
  double tmp_mu;
  double x;
  /*printf("in serial_make_histogram, nparts: %d sample_len: %d cp9_stats: %d do_ins: %d do_enf: %d\n", num_partitions, L, doing_cp9_stats, (cm->search_opts & CM_SEARCH_INSIDE), (cm->config_opts & CM_CONFIG_ENFORCE));*/

  /* Check contract */
  if(cmstats == NULL)
    Die("ERROR in serial_make_histogram(), cmstats is NULL\n");

  /*if ((r = esl_randomness_CreateTimeseeded() == NULL */
  if ((r = esl_randomness_Create(33)) == NULL)
    esl_fatal("Failed to create random number generator: probably out of memory");
  /* Allocate for random distribution */
  ESL_ALLOC(nt_p,   sizeof(double) * Alphabet_size);
  ESL_ALLOC(randseq, sizeof(char) * (L+2));

  /* Configure the CM based on the stat mode */
  ConfigForEVDMode(cm, evd_mode);
  evd = cmstats->evdAA[evd_mode];

  /* For each partition */
  for (p = 0; p < cmstats->np; p++)
    {
      /* Initialize histogram; these numbers are guesses */
      h = esl_histogram_CreateFull(0., 100., 1.);    
      
      /* Set up cur_gc_freq */
      esl_vec_DSet(cur_gc_freq, GC_SEGMENTS, 0.);
      for (i = cmstats->ps[p]; i < cmstats->pe[p]; i++) 
	cur_gc_freq[i] = gc_freq[i];
      esl_vec_DNorm(cur_gc_freq, GC_SEGMENTS);
      
      /* Take N samples */
      for (i = 0; i < N; i++) 
	{
	  /* Get random GC content */
	  gc_comp = 0.01 * esl_rnd_DChoose(r, cur_gc_freq, GC_SEGMENTS);
	  nt_p[1] = nt_p[2] = 0.5 * gc_comp;
	  nt_p[0] = nt_p[3] = 0.5 * (1. - gc_comp);
	  
	  /* Get random sequence */
	  /* We have to generate a text sequence for now, b/c the digitized
	   * version els_rnd_xIID() generates a ESL_DSQ seq, which actually_search_target()
	   * can't handle.
	   */
	  if (esl_rnd_IID(r, Alphabet, nt_p, Alphabet_size, L, randseq)  != eslOK) 
	    esl_fatal("Failure creating random seq for GC content distro");
	  /* Digitize the sequence, parse it, and add to histogram */
	  dsq = DigitizeSequence (randseq, L);
	  /* Do the scan */
	  score = 
	    actually_search_target(cm, dsq, 1, L,
				   0.,    /* cutoff is 0 bits (actually we'll find highest
					   * negative score if it's < 0.0) */
				   0.,    /* CP9 cutoff is 0 bits */
				   NULL,  /* don't keep results */
				   FALSE, /* don't filter with a CP9 HMM */
				   (!(cm->search_opts & CM_SEARCH_HMMONLY)), /* TRUE if we're calcing CM  stats */
				   (cm->search_opts & CM_SEARCH_HMMONLY),    /* TRUE if we're calcing CP9 stats */
				   NULL);          /* filter fraction N/A */
	  /* Add best score to histogram */
	  esl_histogram_Add(h, score);
	}
      /* Fit the scores to a Gumbel */
      esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
      esl_gumbel_FitCensored(xv, n, z, xv[0], &(evd[p]->mu), &(evd[p]->lambda));
      evd[p]->N = N;
      evd[p]->L = L;
    }
  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  free(nt_p);
  free(randseq);

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

/* Function: CopyCMStatsEVD()
 * Incept:   EPN, Mon May  7 06:04:58 2007
 * 
 * Purpose:  Copy the EVD stats in a source CMStats_t object into
 *           a pre-alloc'ed destination CMStats_t object.
 */
int CopyCMStatsEVD(CMStats_t *src, CMStats_t *dest)
{
  int i, p;

  /* Check contract */
  if(src->np != dest->np)
    Die("ERROR in CopyCMStatsEVD() src->np: %d not equal to alloc'ed dest->np: %d\n", src->np, dest->np);

  for(p = 0; p < src->np; p++)
    {
      dest->ps[p] = src->ps[p];
      dest->pe[p] = src->pe[p];
    }
  for(i = 0; i < GC_SEGMENTS; i++)
    dest->gc2p[i] = src->gc2p[i]; 

  for(i = 0; i < NEVDMODES; i++)
    {
      for(p = 0; p < src->np; p++)
	{
	  dest->evdAA[i][p]->N      = src->evdAA[i][p]->N;
	  dest->evdAA[i][p]->L      = src->evdAA[i][p]->L;
	  dest->evdAA[i][p]->mu     = src->evdAA[i][p]->mu;
	  dest->evdAA[i][p]->lambda = src->evdAA[i][p]->lambda;
	}
    }
  return eslOK;
}



#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
/*
 * Function: mpi_make_histogram()
 * Date:     EPN, Wed May  9 13:23:37 2007
 * Purpose:  Parallel version of serial_make_histogram for use with MPI.
 *
 * Args:
 *           cm       - the model
 *           cmstats  - data structure that will hold EVD stats (NULL if my_rank > 0)
 *           gc_freq  - GC frequencies for each GC content      (NULL if my_rank > 0)
 *           N        - number of random seqs to search   (irrelevant if my_rank > 0)
 *           L        - length of random samples          (irrelevant if my_rank > 0)
 *           evd_mode - gives search strategy to use for CM or CP9
 *           my_rank  - MPI rank
 */  
static int mpi_make_histogram (CM_t *cm, CMStats_t *cmstats, double *gc_freq, 
					int N, int L, int evd_mode, int my_rank, int nproc)
{
  int status = 0;
  /*printf("in mpi_make_histogram, my_rank: %d\n", my_rank);*/
  /* Configure the CM based on the stat mode */
  ConfigForEVDMode(cm, evd_mode);
  
  if(my_rank > 0)
    {
      mpi_worker_search_target(cm, my_rank);
      return eslOK;
    }
  else /* my_rank == 0, master */
    {
      int i;
      char *randseq;
      char *dsq;
      ESL_HISTOGRAM *h;
      double  *nt_p;		     /* Distribution for random sequences */
      float score;
      double cur_gc_freq[GC_SEGMENTS];
      float gc_comp;
      double *xv;
      int z;
      int n;
      ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
      EVDInfo_t **evd; 
      int p; /* counter over partitions */
      double tmp_mu;
      double x;
      int              have_work;
      int              nproc_working;
      int              wi;
      MPI_Status       mstatus;
      float sc;
      char        **dsqlist = NULL;      /* queue of digitized seqs being worked on, 1..nproc-1 */

      /* Check contract */
      if(cmstats == NULL)
	Die("ERROR in mpi_make_histogram(), master's cmstats is NULL\n");
      evd = cmstats->evdAA[evd_mode];

      /*if ((r = esl_randomness_CreateTimeseeded() == NULL */
      if ((r = esl_randomness_Create(33)) == NULL)
	esl_fatal("Failed to create random number generator: probably out of memory");
      /* Allocate for random distribution */
      ESL_ALLOC(nt_p,   sizeof(double) * Alphabet_size);
      ESL_ALLOC(randseq, sizeof(char) * (L+2));
      ESL_ALLOC(dsqlist, sizeof(char *) * nproc);
      for (wi = 0; wi < nproc; wi++) dsqlist[wi] = NULL;
      
      /* For each partition */
      for (p = 0; p < cmstats->np; p++)
	{
	  /* Initialize histogram; these numbers are guesses */
	  h = esl_histogram_CreateFull(0., 100., 1.);    
	      
	  /* Set up cur_gc_freq */
	  esl_vec_DSet(cur_gc_freq, GC_SEGMENTS, 0.);
	  for (i = cmstats->ps[p]; i < cmstats->pe[p]; i++) 
	    cur_gc_freq[i] = gc_freq[i];
	  esl_vec_DNorm(cur_gc_freq, GC_SEGMENTS);
	  
	  /* From Sean's design pattern in mpi-hmmsim.c */
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
		  gc_comp = 0.01 * esl_rnd_DChoose(r, cur_gc_freq, GC_SEGMENTS);
		  nt_p[1] = nt_p[2] = 0.5 * gc_comp;
		  nt_p[0] = nt_p[3] = 0.5 * (1. - gc_comp);
		  /* Get random sequence */
		  /* We have to generate a text sequence for now, b/c the digitized
		   * version els_rnd_xIID() generates a ESL_DSQ seq, which actually_search_target()
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
		  esl_histogram_Add(h, sc);
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
	    } /* we have all the scores for this partition*/
	  /* Fit the scores to a Gumbel */
	  esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
	  esl_gumbel_FitCensored(xv, n, z, xv[0], &(evd[p]->mu), &(evd[p]->lambda));
	  evd[p]->N = N;
	  evd[p]->L = L;
	}
      esl_histogram_Destroy(h);
      esl_randomness_Destroy(r);
      free(nt_p);
      free(randseq);
      return eslOK;
    }
 ERROR:
  return status;
}


#endif


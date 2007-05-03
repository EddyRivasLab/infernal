/************************************************************
 * @LICENSE@
 ************************************************************/
/* cmcalibrate.c
 * Score a CM and a CM Plan 9 HMM against random sequence 
 * data to set the statistical parameters for E-value determination,
 * and HMM filtering thresholds. 
 * 
 * EPN, Wed May  2 07:02:52 2007
 * based on HMMER-2.3.2's hmmcalibrate.c from SRE
 *  
 */
#include "config.h"	
#include "esl_config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "funcs.h"		/* external functions                   */
#include "stats.h"              /* EVD functions */
#include "cm_dispatch.h"	
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_vectorops.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"

#define ALGORITHMS "--cyk,--inside"

static int NEW_serial_make_histogram (CM_t *cm, double *gc_freq, int N, int L,
				       int statmode);
static int  FindHMMFilterThreshold(CM_t *cm, float fraction, int N, float cm_pcutoff, 
				   int statmode);
static int debug_print_cmstats(CMStats_t *cmstats);
static int debug_print_evdinfo(EVDInfo_t *evd);
static int debug_print_filterthrinfo(CP9FilterThr_t *fthr);

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "cmcalibrate [-options]";

int
main(int argc, char **argv)
{
  //int status;			       /* status of a function call               */
  //ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  CM_t            *cm      = NULL;     /* the covariance model                    */
  //int              nseq;	       /* counter over sequences                  */
  //char             errbuf[eslERRBUFSIZE];
  char            *cmfile;             /* file to read CM(s) from                 */
  CMFILE          *cmfp;        /* open CM file for reading                 */
  //enum { DO_CYK, DO_INSIDE } algorithm_choice;

  //  int  *partitions;             /* partition each GC % point seg goes to    */
  //int   num_partitions = 1;     /* number of partitions                     */
  int   cm_num_samples;         /* # samples used to calculate  CM EVDs     */
  //int   hmm_num_samples;        /* # samples used to calculate HMM EVDs     */
  int   sample_length;          /* sample len used for calc'ing stats (2*W) */
  int   do_qdb = TRUE;
  int   do_cp9_stats = TRUE;
  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    return eslOK;
  }

  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }
  cmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL); /* NULL=no range checking */

  /*****************************************************************
   * Initializations, including opening and reading the HMM file 
   *****************************************************************/
  /*if ((r = esl_randomness_CreateTimeseeded()) == NULL)*/
  if ((r = esl_randomness_Create(33)) == NULL)
    esl_fatal("Failed to create random number generator: probably out of memory");

  /* currently set up for a single CM - temporary */
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);
  if(do_qdb)          cm->config_opts |= CM_CONFIG_QDB;
  if(do_cp9_stats)    cm->search_opts |= CM_SEARCH_CP9STATS;

  cm->config_opts |= CM_CONFIG_LOCAL;

  ConfigCM(cm, NULL, NULL);

  /****************************************************************
   * Get distribution of GC content from a long random sequence
   *****************************************************************/
  /* human_chr1.gc.code from ~/notebook/7_0502_inf_cmcalibrate/gc_distros/ */
  /* Replace with RFAMSEQ derived one */
  double *gc_freq = MallocOrDie(sizeof(double) * GC_SEGMENTS); /* this should really be doubles */
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

  /****************************************************************
   * Determine CM and CP9 EVDs 
   *****************************************************************/
  CMStats_t *cmstats;
  int N; /* N is the database size (2 * length) */
  int statmode;
  int i;

  N = 2000000;
  sample_length = 2.0 * cm->W;
  cm_num_samples       = 1000;
  cmstats = AllocCMStats(1);
  cmstats->ps[0] = 0;
  cmstats->pe[0] = 100;
  for(i = 0; i < GC_SEGMENTS; i++)
    cmstats->gc2p[i] = 0; 
  cm->stats = cmstats;

  for(statmode = 0; statmode < NSTATMODES; statmode++)
    //for(statmode = 0; statmode < 1; statmode++)
    {
      printf("%-40s ... ", "Determining EVD "); fflush(stdout);
      NEW_serial_make_histogram (cm, gc_freq, cm_num_samples, sample_length, statmode);
      printf("done.\n");
    }

  /****************************************************************
   * Determine CP9 filtering thresholds
   *****************************************************************/
  float fraction = 0.95;
  float cm_pcutoff = 0.001;
  for(statmode = 0; statmode < NFTHRMODES; statmode++)
    //for(statmode = 0; statmode < 1; statmode++)
    {
      printf("%-40s ... ", "Determining HMM filter threshold "); fflush(stdout);
      FindHMMFilterThreshold(cm, fraction, 1000, cm_pcutoff, statmode);
    }
  debug_print_cmstats(cm->stats);
  
  FreeCMStats(cmstats);

  FreeCM(cm);
  return eslOK;
}

/*
 * Function: NEW_serial_make_histogram()
 * Date:     Mon Apr 1 2002 [St. Louis]
 * Purpose:  Makes histogram(s) using random sequences.  Fills the relevant
 *           EVD information in the CMStats_t data structure of the CM with
 *           K and lambda. Determines desired search strategy and 
 *           local/glocal CM/CP9 configuration from statmode.
 *           One histogram is made for each partition of GC frequency. 
 *
 *           This function should be called 6 times, once with each statmode
 *           below, to fill all the EVD stats in the cm->stats data structure:
 *
 *           0. CM_LC: !cm->search_opts & CM_SEARCH_INSIDE  w/  local CM
 *           1. CM_GC: !cm->search_opts & CM_SEARCH_INSIDE  w/ glocal CM
 *           2. CM_LI:  cm->search_opts & CM_SEARCH_INSIDE  w/  local CM
 *           3. CM_GI:  cm->search_opts & CM_SEARCH_INSIDE  w/ glocal CM
 *           4. CP9_L:  cm->search_opts & CM_SEARCH_HMMONLY w/  local CP9 HMM
 *           5. CP9_G:  cm->search_opts & CM_SEARCH_HMMONLY w/ glocal CP9 HMM
 *
 * Args:
 *           cm       - the model, including allocated cm->stats 
 *           gc_freq  - GC frequencies for each GC content
 *           N        - number of random seqs to search 
 *           L        - length of random samples
 *           statmode - gives search strategy to use for CM or CP9
 *
 */  
static int NEW_serial_make_histogram (CM_t *cm, double *gc_freq, int N, int L, int statmode)
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
  /*printf("in serial_make_histogram, nparts: %d sample_len: %d cp9_stats: %d do_ins: %d do_enf: %d\n", num_partitions, L, doing_cp9_stats, (cm->search_opts & CM_SEARCH_INSIDE), (cm->config_opts & CM_CONFIG_ENFORCE));*/

  /* Check contract */
  if(cm->stats == NULL)
    Die("ERROR in serial_make_histogram(), cm->stats is NULL\n");

  /*if ((r = esl_randomness_CreateTimeseeded() == NULL */
  if ((r = esl_randomness_Create(33)) == NULL)
    esl_fatal("Failed to create random number generator: probably out of memory");
  /* Allocate for random distribution */
  ESL_ALLOC(nt_p,   sizeof(double) * Alphabet_size);
  ESL_ALLOC(randseq, sizeof(char) * (L+1));

  /* Configure the CM based on the stat mode */
  ConfigForStatmode(cm, statmode);
  evd = cm->stats->evdAA[statmode];

  /* For each partition */
  for (p = 0; p < cm->stats->np; p++)
    {
      /* Initialize histogram; these numbers are guesses */
      h = esl_histogram_CreateFull(0., 100., 1.);    
      
      /* Set up cur_gc_freq */
      esl_vec_DSet(cur_gc_freq, GC_SEGMENTS, 0.);
      for (i = cm->stats->ps[p]; i < cm->stats->pe[p]; i++) 
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
	  /*if(i % 100 == 0)
	    printf("(%4d) SCORE: %f\n", i, score);*/
	  /* Add best score to histogram */
	  esl_histogram_Add(h, score);
	}
      /* Fit the scores to a Gumbel */
      esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
      esl_gumbel_FitCensored(xv, n, z, xv[0], &(evd[p]->mu), &(evd[p]->lambda));
      evd[p]->K = exp(evd[p]->mu * evd[p]->lambda) / L;
    }
  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  free(nt_p);
  free(randseq);

 ERROR:
  return status;
}


/*
 * Function: FindHMMFilterThreshold
 * Incept:   EPN, Wed May  2 10:00:45 2007
 *
 * Purpose:  Sample sequences from a CM and determine the CP9 HMM bit score
 *           threshold necessary to recognize a specified fraction of those
 *           hits. Fills the relevant information in the CMStats_t data
 *           strcture of the CM (cm->stats). Determines what CM local/glocal
 *           configuration should be using mode. Each call sets 
 *           cm->stats->fthr_*[i]->gT and cm->stats->fthr_*[i]->lT for
 *           "*" equal to one of the following depending on mode:
 *
 *           fthr_lc if mode = FTHR_LC (0): score w/CM in  local CYK 
 *           fthr_gc if mode = FTHR_GC (0): score w/CM in glocal CYK 
 *           fthr_li if mode = FTHR_LI (0): score w/CM in  local inside
 *           fthr_gi if mode = FTHR_GI (0): score w/CM in glocal inside
 *
 *           This function takes a while because emitted parsetree score
 *           not necessarily (a) optimal nor (b) highest score returned
 *           from a CYK/Inside scan (a subseq of the full seq could score 
 *           higher even if parsetree was optimal). This means we have
 *           to scan each emitted seq with CYK or Inside.
 *
 * Args:
 *           cm           - the CM
 *           fraction     - target fraction of CM hits to detect with CP9
 *           N            - number of sequences to sample from CM better than cm_minsc
 *           cm_pcutoff   - minimum CM P-value to accept 
 *           statmode     - gives CM search strategy to use, and EVD to use
 */
int FindHMMFilterThreshold(CM_t *cm, float fraction, int N, float cm_pcutoff, int statmode)
{
  char            **dsq;        /* digitized sequences                    */
  char            **seq;        /* alphabetic sequences                   */
  float             sc;
  float            *hmm_local_pval;
  float            *hmm_glocal_pval;
  int               nattempts = 0;
  int               max_attempts;
  int              *L;
  int               i;
  int               gc_comp;
  CP9FilterThr_t **fthr;
  EVDInfo_t **evd; 
  int passed;
  int *p;                   /* partition each seq belongs in based on 1..L gc_comp (not entirely correct) */
  int  pctr;                /* counter over partitions */
  double pval;
  /* Contract check */
  if (cm->cp9 == NULL) 
    Die("ERROR in FindHMMFilterThreshold() cp9 is NULL\n");
  if (fraction < 0. || fraction > 1.)  
   Die("ERROR in FindHMMFilterThreshold() fraction is %f, should be [0.0..1.0]\n", fraction);

  /* Configure the CM based on the stat mode */
  ConfigForStatmode(cm, statmode);
  fthr = cm->stats->fthrAA[statmode];
  evd  = cm->stats->evdAA[statmode];

  /* Strategy: emit sequences one at a time from CM, scanning each with CM 
   * (local or glocal, CYK or Inside as specified by mode). If best scan
   * score meets cm_minsc threshold, keep it. Do this until we've collected
   * N sequences or reached max_attempts. This is wasteful in that we save
   * all N sequences, but saves us from having to reconfigure HMM from glocal
   * to local for each seq.
   */

  /* Emit seqs until we have N that meet score cutoff */
  max_attempts = 500 * N;
  seq    = MallocOrDie(sizeof(char *) * N);
  dsq    = MallocOrDie(sizeof(char *) * N);
  L      = MallocOrDie(sizeof(char *) * N);
  p      = MallocOrDie(sizeof(int *) * N);
  for (i = 0; i < N; i++)
    {
      if(nattempts++ > max_attempts) 
	Die("ERROR number of attempts exceeded 500 times number of seqs.\n");
      dsq[i] = seq[i] = NULL;
      passed = FALSE;
      while(!passed && nattempts <= max_attempts)
	{
	  if(dsq[i] != NULL) free(dsq[i]);
	  if(seq[i] != NULL) free(seq[i]);
	  EmitParsetree(cm, NULL, &seq[i], &dsq[i], &L[i]); /* we don't care about the parsetree, 
							  * we have to scan each seq to get score */
	  sc = actually_search_target(cm, dsq[i], 1, L[i],
				      0.,    /* cutoff is 0 bits (actually we'll find highest
					      * negative score if it's < 0.0) */
				      0.,    /* CP9 cutoff is 0 bits */
				      NULL,  /* don't keep results */
				      FALSE, /* don't filter with a CP9 HMM */
				      FALSE, /* we're not calcing CM  stats */
				      FALSE, /* we're not calcing CP9 stats */
				      NULL); /* filter fraction N/A */
	  /* Only accept if P-value of best hit in this seq is better than our cutoff.
	   * To do that we first determine which partition to use EVD from */
	  gc_comp = get_gc_comp(seq[i], 1, L[i]); /* should be i and j of best hit, but
						* EVD construction is wrong too, it uses 
						* GC_COMP of full sample seq, not of best hit
						* within it. */
	  p[i] = cm->stats->gc2p[gc_comp];
	  pval = esl_gumbel_surv((double) sc, evd[p[i]]->mu, evd[p[i]]->lambda);
	  //printf("i: %d pval: %f\n", i, pval);
	  if(pval < cm_pcutoff)
	     passed = TRUE;
	  nattempts++;
	}
      if(nattempts > max_attempts)
	Die("ERROR number of attempts exceeded 50 times number of seqs.\n");
    }

  /****************************************
   * Scan each seq with HMM w/LOCAL Forward
   ****************************************/
  cm->search_opts |= CM_SEARCH_HMMONLY; /* turn on HMMONLY searching */    
  CPlan9SWConfig(cm->cp9, ((cm->cp9->M)-1.)/cm->cp9->M,  /* all start pts equiprobable, including 1 */
		 ((cm->cp9->M)-1.)/cm->cp9->M);          /* all end pts equiprobable, including M */
  CP9Logoddsify(cm->cp9);

  hmm_local_pval  = MallocOrDie(sizeof(float) * N);
  for (i = 0; i < N; i++)
    {
      sc = actually_search_target(cm, dsq[i], 1, L[i],
				  cm->cutoff, cm->cp9_cutoff,
				  NULL,  /* don't report hits to a results structure */
				  FALSE, /* we're not filtering with a CP9 HMM */
				  FALSE, /* we're not building a histogram for CM stats  */
				  FALSE, /* we're not building a histogram for CP9 stats */
				  NULL); /* filter fraction, irrelevant here */
      hmm_local_pval[i] = esl_gumbel_surv((double) sc, 
					cm->stats->evdAA[CP9_L][p[i]]->mu, 
					cm->stats->evdAA[CP9_L][p[i]]->lambda);
      printf("L i: %4d hmm sc: %10.4f hmm P: %10.4f\n", i, sc, hmm_local_pval[i]);
    }

  /****************************************
   * Scan each seq with HMM w/GLOCAL Forward
   ****************************************/
  CPlan9GlobalConfig(cm->cp9);
  CP9Logoddsify(cm->cp9);
  hmm_glocal_pval = MallocOrDie(sizeof(float) * N);
  for (i = 0; i < N; i++)
    {
      sc = actually_search_target(cm, dsq[i], 1, L[i],
				  cm->cutoff, cm->cp9_cutoff,
				  NULL,  /* don't report hits to a results structure */
				  FALSE, /* we're not filtering with a CP9 HMM */
				  FALSE, /* we're not building a histogram for CM stats  */
				  FALSE, /* we're not building a histogram for CP9 stats */
				  NULL); /* filter fraction, irrelevant here */
      hmm_glocal_pval[i] = esl_gumbel_surv((double) sc, 
					   cm->stats->evdAA[CP9_G][p[i]]->mu, 
					   cm->stats->evdAA[CP9_G][p[i]]->lambda);
      printf("G i: %4d hmm sc: %10.4f hmm P: %10.4f\n", i, sc, hmm_glocal_pval[i]);
    }

  /* Sort the HMM P-values with quicksort */
  esl_vec_FSortDecreasing(hmm_local_pval, N);
  esl_vec_FSortDecreasing(hmm_glocal_pval, N);
  /* Set thresholds */
  /* TEMPORARY set thresholds for all partitions as equal
   * you have to decide if you want to set them differently, or
   * just have one filter threshold for ALL partitions, in which
   * case you'll need to change the CP9FilterThr_t data structure */
  for(pctr = 0; pctr < cm->stats->np; pctr++)
    {
      fthr[pctr]->gsc = hmm_glocal_pval[(int) ((1. - fraction) * (float) N)];
      fthr[pctr]->lsc = hmm_local_pval[(int) ((1. - fraction) * (float) N)];
      fthr[pctr]->cmsc = cm_pcutoff;
      fthr[pctr]->fraction = fraction;
      fthr[pctr]->is_pval = TRUE;
    }
  printf("\n\nnattempts: %d\n", nattempts);
  free(hmm_local_pval);
  free(hmm_glocal_pval);
  free(p);
  free(L);
  for(i = 0; i < N; i++)
    {
      if(dsq[i] != NULL) free(dsq[i]);
      if(seq[i] != NULL) free(seq[i]);
    }
  free(dsq);
  free(seq);
  return eslOK;
}

/* Function: debug_print_cmstats
 */
int debug_print_cmstats(CMStats_t *cmstats)
{
  int p;
  printf("Num partitions: %d\n", cmstats->np);
  for (p = 0; p < cmstats->np; p++)
    {
      printf("Partition %d: start: %d end: %d\n", p, cmstats->ps[p], cmstats->pe[p]);
      printf("cm_lc EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_LC][p]);
      printf("cm_gc EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_GC][p]);
      printf("cm_li EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_LI][p]);
      printf("cm_gi EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_GI][p]);
      printf("cp9_l EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CP9_L][p]);
      printf("cp9_g EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CP9_G][p]);
      printf("\n\n");
      printf("fthr lc filter threshold:\t");
      debug_print_filterthrinfo(cmstats->fthrAA[CM_LC][p]);
      printf("fthr gc filter threshold:\t");
      debug_print_filterthrinfo(cmstats->fthrAA[CM_GC][p]);
      printf("fthr li filter threshold:\t");
      debug_print_filterthrinfo(cmstats->fthrAA[CM_LI][p]);
      printf("fthr gi filter threshold:\t");
      debug_print_filterthrinfo(cmstats->fthrAA[CM_GI][p]);
      printf("\n\n");
    }
  return 0;
}

/* Function: debug_print_evdinfo
 */
int debug_print_evdinfo(EVDInfo_t *evd)
{
  printf("lambda: %.5f mu: %.5f K: %.5f\n", evd->lambda, evd->mu, evd->K);
  return 0;
}

/* Function: debug_print_filterthrinfo
 */
int debug_print_filterthrinfo(CP9FilterThr_t *fthr)
{
  printf("gsc: %.5f lsc: %.5f cmsc: %.5f fraction: %.3f is_pval: %d\n", 
	 fthr->gsc, fthr->lsc, fthr->cmsc, fthr->fraction, fthr->is_pval);
  return 0;
}


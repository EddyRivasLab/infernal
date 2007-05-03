/*/************************************************************
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
static int debug_print_evdinfo(EVDInfo_t *evd);
static int debug_print_cmstats(CMStats_t *cmstats);

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "cmcalibrate [-options]";

int
main(int argc, char **argv)
{
  int status;			       /* status of a function call               */
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  CM_t            *cm      = NULL;     /* the covariance model                    */
  int              sc;		       /* a CYK or Inside score                   */
  int              nseq;	       /* counter over sequences                  */
  char             errbuf[eslERRBUFSIZE];
  char            *cmfile;             /* file to read CM(s) from                 */
  CMFILE          *cmfp;        /* open CM file for reading                 */
  enum { DO_CYK, DO_INSIDE } algorithm_choice;

  int  *partitions;             /* partition each GC % point seg goes to    */
  int   num_partitions = 1;     /* number of partitions                     */
  int   cm_num_samples;         /* # samples used to calculate  CM EVDs     */
  int   hmm_num_samples;        /* # samples used to calculate HMM EVDs     */
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
  int i;
  int N; /* N is the database size (2 * length) */
  int statmode;
  N = 2000000;
  sample_length = 2.0 * cm->W;
  cm_num_samples       = 1000;
  cmstats = AllocCMStats(1);
  cmstats->ps[0] = 0;
  cmstats->pe[0] = 100;
  cm->stats = cmstats;

  for(statmode = 0; statmode < NSTATMODES; statmode++)
    {
      printf("%-40s ... ", "Determining EVD "); fflush(stdout);
      NEW_serial_make_histogram (cm, gc_freq, cm_num_samples, sample_length, statmode);
      printf("done.\n");
    }

  debug_print_cmstats(cm->stats);

  //SE EVD STATS
  //cm->lambda[i] = 0.509238; //comment me
  //cm->K[i] = .014367;    //comment me
  //cm->mu[i] = log(cm->K[i]*N)/cm->lambda[i];

  /****************************************************************
   * Determine CP9 filtering thresholds
   *****************************************************************/
  float globE = 50.;
  double globP = 0.01;
  double globT;
  /* determine CM bit score that gives E value of globE */
  globT = esl_gumbel_invcdf(globP, cm->mu[50], cm->lambda[50]);
  printf("bit score that gives P value of %f: %f\n", globP, globT);
  printf("E value of score: %f is %f\n", globT, RJK_ExtremeValueE(globT, cm->mu[50], cm->lambda[50]));

  /* First, glocal CP9 to glocal CM */
  float fraction = 0.95;
  float thresh;
  /*for(fthrmode = 0; fthrmode < NFTHRMODES; fthrmode++)
    {
    printf("%-40s ... ", "Determining EVD "); fflush(stdout);*/
  thresh = 
    FindHMMFilterThreshold(cm, fraction, globT, FALSE,  /* glocal CM  */
			   FALSE,      /* glocal CP9 */
			   fraction, globT, 1000);
  printf("Glocal glocal thresh: %.3f\n", thresh);
  thresh = 
    FindHMMFilterThreshold(cm, FALSE, /* glocal CM  */
			   TRUE,      /* local CP9 */
			   fraction, globT, 1000);
  printf("Glocal local thresh: %.3f\n", thresh);
  thresh = 
    FindHMMFilterThreshold(cm, TRUE,   /* local CM  */
			   FALSE,      /* glocal CP9 */
			   fraction, globT, 1000);
  printf("Local glocal thresh: %.3f\n", thresh);
  thresh = 
    FindHMMFilterThreshold(cm, TRUE,  /* local CM  */
			   TRUE,      /* local CP9 */
			   fraction, globT, 1000);
  printf("Local local thresh: %.3f\n", thresh);


  /* end of from cmemit.c */
  
  FreeCMStats(cmstats);

  FreeCM(cm);
  return eslOK;
}

/*
 * Function: NEW_serial_make_histogram()
 * Date:     Mon Apr 1 2002 [St. Louis]
 * Purpose:  Makes histogram(s) using random sequences.  Fills the relevant
 *           EVD information in the CMStats_t data structure of the CM with
 *           K and lambda. Determines search strategy using cm->search_opts
 *           and local/glocal configuration of CM using cm->flags.
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
 *           statmode - 
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
  int do_cm_local = FALSE;
  int do_cp9_local = FALSE;
  EVDInfo_t **evd; 
  int p; /* counter over partitions */
  /*printf("in serial_make_histogram, nparts: %d sample_len: %d cp9_stats: %d do_ins: %d do_enf: %d\n", num_partitions, L, doing_cp9_stats, (cm->search_opts & CM_SEARCH_INSIDE), (cm->config_opts & CM_CONFIG_ENFORCE));*/

  /* Check contract */
  if(cm->stats == NULL)
    Die("ERROR in serial_make_histogram(), cm->stats is NULL\n");

  /*if ((r = esl_randomness_CreateTimeseeded() == NULL */
  if ((r = esl_randomness_Create(33)) == NULL)
    esl_fatal("Failed to create random number generator: probably out of memory");

  /* Configure CM based on mode, and save pointer
   * to EVD params we're filling */
  switch (statmode) {
  case CM_LC: /* local CYK */
    printf("CM_LC\n");
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMONLY;
    do_cm_local  = TRUE;
    evd = cm->stats->cm_lc;
    break;
  case CM_GC: /* glocal CYK */
    printf("CM_GC\n");
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMONLY;
    do_cm_local  = FALSE;
    evd = cm->stats->cm_gc;
    break;
  case CM_LI: /* local inside */
    printf("CM_LI\n");
    cm->search_opts |= CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMONLY;
    do_cm_local  = TRUE;
    evd = cm->stats->cm_li;
    break;
  case CM_GI: /* glocal inside */
    printf("CM_GI\n");
    cm->search_opts |= CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMONLY;
    do_cm_local  = FALSE;
    evd = cm->stats->cm_gi;
    break;
  case CP9_L: /* local CP9 Forward */
    printf("CP9_L\n");
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts |= CM_SEARCH_HMMONLY;
    do_cp9_local  = TRUE;
    evd = cm->stats->cp9_l;
    break;
  case CP9_G: /* glocal CP9 Forward */
    printf("CP9_G\n");
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts |= CM_SEARCH_HMMONLY;
    do_cp9_local  = FALSE;
    evd = cm->stats->cp9_g;
    break;
  default: 
    Die("ERROR unrecognized statmode: %d in serial_make_histogram");
  }    
  /* configure CM and, if needed, CP9 */
  if(cm->search_opts & CM_SEARCH_HMMONLY)
    {
      printf("configuring HMM\n");
      if(do_cp9_local)
	CPlan9SWConfig(cm->cp9, ((cm->cp9->M)-1.)/cm->cp9->M,  /* all start pts equiprobable, including 1 */
		     ((cm->cp9->M)-1.)/cm->cp9->M);          /* all end pts equiprobable, including M */
      else
	CPlan9GlobalConfig(cm->cp9);
      CP9Logoddsify(cm->cp9);
    }
  else
    {
      if(do_cm_local)
	ConfigLocal(cm, 0.5, 0.5);
      else
	ConfigGlobal(cm);
      CMLogoddsify(cm);
    }
  /* Allocate for random distribution */
  ESL_ALLOC(nt_p,   sizeof(double) * Alphabet_size);
  ESL_ALLOC(randseq, sizeof(char) * (L+1));

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
	   * version els_rnd_xIID() generates a ESL_DSQ seq, which actually_search_target
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
      /* Fit the histogram.  */
      /* Fit the scores to a Gumbel */
      /* If the esl_histogram example for 'complete data, high scores fit as
       *  censored Gumbel' is the correct approach we do this: */
      esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
      
      /* If the esl_histogram example for 'censored data, fit as censored Gumbel
       *  is the correct approach, we would do this (BUT IT'S NOT!): */
      /*esl_histogram_GetData(h, &xv, &n);*/
      
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
      debug_print_evdinfo(cmstats->cm_lc[p]);
      printf("cm_gc EVD:\t");
      debug_print_evdinfo(cmstats->cm_gc[p]);
      printf("cm_li EVD:\t");
      debug_print_evdinfo(cmstats->cm_li[p]);
      printf("cm_gi EVD:\t");
      debug_print_evdinfo(cmstats->cm_gi[p]);
      printf("cp9_l EVD:\t");
      debug_print_evdinfo(cmstats->cp9_l[p]);
      printf("cp9_g EVD:\t");
      debug_print_evdinfo(cmstats->cp9_g[p]);
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

/* fastsearch.c
 * EPN, Wed Sep 12 16:53:32 2007
 * 
 * Fast versions of CYK and Inside search functions.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"
   
#define TSC(s,k) (tsc[(k) * cp9O_NTRANS + (s)])

/***********************************************************************
 * Function: cp9_FastViterbi()
 * 
 * Purpose:  Runs the Viterbi dynamic programming algorithm on an
 *           input subsequence (i0-j0). 
 *           Somewhat flexible based on input options.
 *    
 * Note:     IDENTICAL to CP9Forward() below with maxes replacing
 *           sums in DP recursion, and the possibility of returning a
 *           CP9 trace if doing_align == TRUE. See CP9Forward() for more info,
 *           including more verbose 'Purpose' and description of arguments.
 *
 * Returns:  if(!do_scan) log P(S,tr|M)/P(S,tr|R), as a bit score
 *           else         max log P(S,tr|M)/P(S,tr|R), for argmax subseq S of input seq i0..j0,
 */
float
cp9_FastViterbi(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
		int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
		int be_efficient, CP9_dpmatrix_t **ret_mx, CP9trace_t **ret_tr)
{
  int          status;
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /* j-W: position in the subsequence                             */
  int          ip;          /* i': relative position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  int          L;           /* j0-i0+1: subsequence length                                  */
  CP9_dpmatrix_t *mx;       /* the CP9 DP matrix                                            */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cm->cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int         *scA;         /* prob (seq from j0..jp | HMM) [0..jp..cm->cp9->M]             */
  int          tmp_sc;      /* temporary int log odds score                                 */
  float        fsc;         /* float log odds score                                         */
  float        curr_sc;     /* temporary score used for filling in gamma                    */
  float       *gamma;       /* SHMM DP matrix for optimum nonoverlap resolution [0..j0-i0+1]*/
  int         *gback;       /* traceback pointers for SHMM                      [0..j0-i0+1]*/ 
  float       *savesc;      /* saves score of hit added to best parse at j      [0..j0-i0+1]*/ 
  float        best_hmm_sc; /* Best overall score from semi-HMM to return if do_scan        */
  float        best_hmm_pos;/* residue giving best_hmm_sc                                   */
  float        best_sc;     /* Best score overall, returned if 0 hits found by HMM & do_scan*/
  float        best_pos;    /* residue giving best_sc                                       */
  float        return_sc;   /* score to return, if (!do_scan) return overall Viterbi sc,    *
			     * else return best_hmm_sc if # HMM hits>0, else return best_sc */
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */
  int          c;           /* counter for EL states */
  CP9trace_t  *tr;
  int          M;
  /*debug_print_cp9_params(stdout, cm->cp9, TRUE);*/

  /*printf("in CP9Viterbi() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in CP9Viterbi, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    esl_fatal("ERROR in CP9Viterbi, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    esl_fatal("ERROR in CP9Viterbi, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    esl_fatal("ERROR in CP9Viterbi, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
  if(dsq == NULL)
    esl_fatal("ERROR in CP9Viterbi, dsq is NULL.");
    
  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_hmm_sc = IMPOSSIBLE;
  best_hmm_pos= -1;
  L = j0-i0+1;
  M = cm->cp9->M;
  if (W > L) W = L; 

  /* Transition scores */
  int *otsc;
  ESL_ALLOC(otsc,   sizeof(int)   * (M+1)  * cp9O_NTRANS);
  for (k = 0 ; k <= M; k++) {
    int *otsc_k = otsc + k*cp9O_NTRANS;
    
    otsc_k[cp9O_MM] = cm->cp9->tsc[CTMM][k];
    otsc_k[cp9O_MI] = cm->cp9->tsc[CTMI][k];
    otsc_k[cp9O_MD] = cm->cp9->tsc[CTMD][k];
    otsc_k[cp9O_IM] = cm->cp9->tsc[CTIM][k];
    otsc_k[cp9O_II] = cm->cp9->tsc[CTII][k];
    otsc_k[cp9O_DM] = cm->cp9->tsc[CTDM][k];
    otsc_k[cp9O_DD] = cm->cp9->tsc[CTDD][k];
    otsc_k[cp9O_ID] = cm->cp9->tsc[CTID][k];
    otsc_k[cp9O_DI] = cm->cp9->tsc[CTDI][k];
    otsc_k[cp9O_BM] = cm->cp9->bsc[k];
    otsc_k[cp9O_MEL]= cm->cp9->tsc[CTME][k];
    otsc_k[cp9O_ME] = cm->cp9->esc[k];
  }
  int const *tsc = otsc;
  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0.;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 2;
  else             nrows = L+1;
  mx = AllocCPlan9Matrix(nrows, M, &mmx, &imx, &dmx, &elmx, &erow); 

  /* scA will hold P(seq up to j | Model) in int log odds form */
  ESL_ALLOC(scA, sizeof(int) * (j0-i0+2));
			
  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  elmx[0][0]= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY;   
  /*printf("mmx[jp:%d][%d]: %d %f\n", 0, 0, mmx[0][0], Score2Prob(mmx[0][0], 1.));
    printf("imx[jp:%d][%d]: %d %f\n", 0, 0, imx[0][0], Score2Prob(imx[0][0], 1.));
    printf("dmx[jp:%d][%d]: %d %f\n", 0, 0, dmx[0][0], Score2Prob(dmx[0][0], 1.));
    printf("elmx[jp:%d][%d]: %d %f\n", 0, 0, elmx[0][0], Score2Prob(elmx[0][0], 1.));*/

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  int sc;
  for (k = 1; k <= M; k++)
    {
      mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
      sc = ESL_MAX(mmx[0][k-1] + TSC(cp9O_MD,k-1),
		   imx[0][k-1] + TSC(cp9O_ID,k-1));
      sc = ESL_MAX(sc, dmx[0][k-1] + TSC(cp9O_DD,k-1));
      dmx[0][k] = sc;
      /*printf("mmx[jp:%d][%d]: %d %f\n", 0, k, mmx[0][k], Score2Prob(mmx[0][k], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", 0, k, imx[0][k], Score2Prob(imx[0][k], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));
	printf("elmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));*/
    }
  /* We can do a full parse through all delete states. */
  erow[0]  = dmx[0][M] + TSC(cp9O_DM,M);
  scA[0]   = erow[0];
  fsc      = Scorify(scA[0]);
  /*printf("jp: %d j: %d fsc: %f isc: %d\n", jp, j, fsc, isc[jp]);*/
  if(fsc > best_sc) 
    {
      best_sc = fsc;
      best_pos= i0-1;
    }
  /*printf("jp: %d j: %d fsc: %f sc: %d\n", 0, i0-1, Scorify(sc[0]), sc[0]);*/

  /*****************************************************************
   * The main loop: scan the sequence from position i0 to j0.
   *****************************************************************/
  /* Recursion. */
  int ctr = 0;
  for (j = i0; j <= j0; j++)
    {
      int const *isc = cm->cp9->isc[dsq[j]];
      int const *msc = cm->cp9->msc[dsq[j]];

      int endsc     = -INFTY;
      int el_selfsc = cm->cp9->el_selfsc;
      int sc;

      jp = j-i0+1;     /* jp is relative position in the sequence 1..L */
      cur = (j-i0+1);
      prv = (j-i0);

      if(be_efficient)
	{
	  cur %= 2;
	  prv %= 2;
	}	  
      /* The 1 difference between a Viterbi scanner and the 
       * regular Viterbi. In non-scanner parse must begin in B at
       * position 0 (i0-1), in scanner we can start at any position 
       * in the seq. */
      mmx[cur][0]  = (do_scan == TRUE) ? 0 : -INFTY;
      dmx[cur][0]  = -INFTY;  /*D_0 is non-existent*/
      elmx[cur][0] = -INFTY;  /*no EL state for node 0 */

      sc = ESL_MAX(mmx[prv][0] + TSC(cp9O_MI,0),
		   imx[prv][0] + TSC(cp9O_II,0));
      sc = ESL_MAX(sc, dmx[prv][0] + TSC(cp9O_DI,0));
      imx[cur][0] = ESL_MAX(sc + isc[0], -INFTY);

      /*printf("mmx[jp:%d][%d]: %d %f\n", jp, 0, mmx[cur][0], Score2Prob(mmx[cur][0], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", jp, 0, imx[cur][0], Score2Prob(imx[cur][0], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", jp, 0, dmx[cur][0], Score2Prob(dmx[cur][0], 1.));*/

      for (k = 1; k <= M; k++)
	{
	  /*match state*/
	  //ctr++;
	  sc = ESL_MAX(    mmx[prv][k-1] + TSC(cp9O_MM,k-1),
			   imx[prv][k-1] + TSC(cp9O_IM,k-1));
	  sc = ESL_MAX(sc, dmx[prv][k-1] + TSC(cp9O_DM,k-1));
	  sc = ESL_MAX(sc, mmx[prv][0]   + TSC(cp9O_BM,k));
	  /* check possibility we came from an EL, if they're valid */
	  for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	    sc = ESL_MAX(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
	    /* transition penalty to EL incurred when EL was entered */
	  mmx[cur][k] = ESL_MAX(sc + msc[k], -INFTY);

	  /* E state update */
	  endsc = ESL_MAX(endsc, mmx[cur][k] + TSC(cp9O_ME,k));

	  /*insert state*/
	  sc = ESL_MAX(    mmx[prv][k] + TSC(cp9O_MI,k),
		           imx[prv][k] + TSC(cp9O_II,k));
	  sc = ESL_MAX(sc, dmx[prv][k] + TSC(cp9O_DI,k));
	  imx[cur][k] = ESL_MAX(sc + isc[k], -INFTY);

	  /*delete state*/
	  sc = ESL_MAX(    mmx[cur][k-1] + TSC(cp9O_MD,k-1),
		           imx[cur][k-1] + TSC(cp9O_ID,k-1));
	  sc = ESL_MAX(sc, dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	  dmx[cur][k] = sc;

	  /*el state*/
	  sc = -INFTY;
	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								    HMM nodes that map to right half of a MATP_MP) */
	    {
	      sc = ESL_MAX(sc, mmx[cur][k]  + TSC(cp9O_MEL,k));
	      sc = ESL_MAX(sc, elmx[prv][k] + el_selfsc);
	    }
	  elmx[cur][k] = sc;

	  /*printf("mmx[jp:%d][%d]: %d %f\n", jp, k, mmx[cur][k], Score2Prob(mmx[cur][k], 1.));
	    printf("imx[jp:%d][%d]: %d %f\n", jp, k, imx[cur][k], Score2Prob(imx[cur][k], 1.));
	    printf("dmx[jp:%d][%d]: %d %f\n", jp, k, dmx[cur][k], Score2Prob(dmx[cur][k], 1.));*/
	}
      endsc = ESL_MAX(endsc, dmx[cur][M] + TSC(cp9O_DM,M)); /* transition from D_M -> end */
      endsc = ESL_MAX(endsc, imx[cur][M] + TSC(cp9O_IM,M)); /* transition from I_M -> end */
      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	/* transition penalty to EL incurred when EL was entered */
	endsc = ESL_MAX(endsc, elmx[cur][cm->cp9->el_from_idx[M+1][c]]);

      erow[cur] = endsc;
      scA[jp]   = endsc;
      fsc = Scorify(endsc);
      /*printf("jp: %d j: %d fsc: %f sc: %d\n", jp, j, fsc, scA[jp]);*/
      if(fsc > best_sc)
	{
	  best_sc = fsc;
	  best_pos= j;
	}
      if (fsc > best_hmm_sc)
	{
	  best_hmm_sc = fsc;
	  best_hmm_pos= j;
	}
      if(!(cm->search_opts & CM_SEARCH_HMMGREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   */
	  gamma[jp]  = gamma[jp-1] + 0; /* extend without adding a new hit */
	  gback[jp]  = -1;
	  savesc[jp] = IMPOSSIBLE;
	  i = ((j-W+1)> i0) ? (j-W+1) : i0;
	  ip = i-i0+1;
	  curr_sc = gamma[ip-1] + fsc + cm->cp9_sc_boost;
	  /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
	   * value is 0.0 if technique not used. */
	  if (curr_sc > gamma[jp])
	    {
	      gamma[jp]  = curr_sc;
	      gback[jp]  = i;
	      savesc[jp] = fsc;
	    }
	}
      else
	{
	  /* Resolving overlaps greedily (RSEARCH style),  
	   * Return best hit for each j, IFF it's above threshold */
	  if (fsc >= cutoff) 
	    {
	      if(results != NULL) 
		{
		  i = ((j-W+1)> i0) ? (j-W+1) : i0;
		  /*printf("VIT greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i, j, fsc);*/
		  report_hit (i, j, 0, fsc, results);
		  /* 0 is for saver, which is irrelevant for HMM hits */
		}
	    }
	}
    } /* end loop over end positions j */
      
  if((!(cm->search_opts & CM_SEARCH_HMMGREEDY)) && /* resolve overlaps optimally */
     (!doing_align || do_scan)) /* else we can save time by skipping traceback */
    {
      /*****************************************************************
       * Traceback stage.
       * Recover all hits: an (i,j,sc) triple for each one and report them.
       *****************************************************************/ 
      j           = j0;
      while (j >= i0) {
	jp = j-i0+1;
	if (gback[jp] == -1) /* no hit */
	  j--; 
	else                /* a hit, a palpable hit */
	  {
	    if(savesc[jp] > best_hmm_sc) 
	      {
		best_hmm_sc = savesc[jp];
		best_hmm_pos= j;
	      }
	    if(savesc[jp] >= cutoff)
	      {
		if(results != NULL) /* report the hit */
		  {
		    /*printf("VIT reporting hit: i: %d j: %d sc: %f\n", gback[jp], j, savesc[jp]);*/
		    report_hit(gback[jp], j, 0, savesc[jp], results); 
		    /* 0 is for saver, which is irrelevant for HMM hits */
		  }
	      }
	    j = gback[jp]-1;
	  }
      }
    }
  /* clean up and exit */
  free(gback);
  free(gamma);
  free(savesc);

  /* determine score to return: (I know, too complex) */
  if(doing_align)
    {
      return_sc = Scorify(scA[(j0-i0+1)]); /* L = j0-i0+1 */
      if(ret_bestpos != NULL) *ret_bestpos = i0;
      if(ret_tr != NULL) 
	{
	  CP9ViterbiTrace(cm->cp9, dsq, i0, j0, mx, &tr);
	  /* CP9PrintTrace(stdout, tr, cm->cp9, dsq); */
	  *ret_tr = tr;
	}
    }
  else if(best_hmm_sc <= 0.) /* scanning and there were no hits found by the 
			      * semi-HMM above 0.0 bits */
    {
      return_sc = best_sc;
      if(ret_bestpos != NULL) *ret_bestpos = best_pos;
    }
  else
    {
      return_sc = best_hmm_sc;
      if(ret_bestpos != NULL) *ret_bestpos = best_hmm_pos;
    }
  if(ret_sc != NULL) *ret_sc = scA;
  else free(scA);
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);
  /*printf("Viterbi return_sc: %f\n", return_sc);*/

  printf("ctr: %d\n", ctr);
  return return_sc;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}


/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_CP9_FASTSEARCH_BENCHMARK
/* gcc -o benchmark-cp9_fastsearch -g -O2 -I. -L. -I../easel -L../easel -DIMPL_CP9_FASTSEARCH_BENCHMARK cp9_fastsearch.c -linfernal -leasel -lm
 * ./benchmark-cp9_fastsearch <cmfile>
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_sqio.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "33", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT, "500000", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-w",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute slow Viterbi scan implementation",  0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute experimental fast Viterbi scan implementation",  0 },
  { "-f",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "calculate full matrix, not just 2 rows",         0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "configure HMM for glocal alignment [default: local]", 0 },
  { "--noel",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "turn local ends off [default: on, unless -g]", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "benchmark driver for the fast scanning CM plan 9 HMM Viterbi implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  CM_t            *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  char            *cmfile = esl_opt_GetArg(go, 1);
  CMFILE          *cmfp;	/* open input CM file stream */

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  cm->config_opts |= CM_CONFIG_QDB;
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    if(! esl_opt_GetBoolean(go, "--noel")) cm->config_opts |= CM_CONFIG_HMMEL;
  }
  ConfigCM(cm, NULL, NULL);

  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, cm->null, abc->K, L, dsq);

      esl_stopwatch_Start(w);
      sc = cp9_FastViterbi(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			   TRUE,   /* we're scanning */
			   FALSE,  /* we're not ultimately aligning */
			   (! esl_opt_GetBoolean(go, "-f")),  /* memory efficient ? */
			   NULL,   /* don't want the DP matrix back */
			   NULL);  /* don't want traces back */
      printf("%4d %-30s %10.4f bits ", (i+1), "cp9_FastViterbi(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      

      if (esl_opt_GetBoolean(go, "-w")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = CP9Viterbi(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			  TRUE,   /* we're scanning */
			  FALSE,  /* we're not ultimately aligning */
			  (! esl_opt_GetBoolean(go, "-f")),  /* memory efficient ? */
			  NULL,   /* don't want the DP matrix back */
			  NULL);  /* don't want traces back */
	  printf("%4d %-30s %10.4f bits ", (i+1), "CP9Viterbi(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}
    }
  FreeCM(cm);
  free(dsq);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IMPL_FASTSEARCH_BENCHMARK*/

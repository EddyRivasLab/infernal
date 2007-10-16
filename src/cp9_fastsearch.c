/* fastsearch.c EPN, Wed Sep 12 16:53:32 2007
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
#include <math.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"
   
#define TSC(s,k) (tsc[(k) * cp9O_NTRANS + (s)])

#define IAMI_M  (1<<0)
#define IAMI_I  (1<<1)
#define IAMI_D  (1<<2)
#define IAMI_EL (1<<3)
#define nIAMI_S   4

#define IAMI_MM  (1<<0)
#define IAMI_IM  (1<<1)
#define IAMI_DM  (1<<2)
#define IAMI_BM  (1<<3)
#define IAMI_MI  (1<<4)
#define IAMI_II  (1<<5)
#define IAMI_DI  (1<<6)
#define IAMI_MD  (1<<7)
#define IAMI_ID  (1<<8)
#define IAMI_DD  (1<<9)
#define IAMI_ME  (1<<10)
#define IAMI_MEL (1<<11)
#define nIAMI_T   12

#define TDIFF_MI_M 0
#define TDIFF_MI_I 1
#define TDIFF_MI_D 2
#define nTDIFF 3

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
    cm_Fail("in CP9Viterbi, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    cm_Fail("in CP9Viterbi, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    cm_Fail("in CP9Viterbi, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    cm_Fail("in CP9Viterbi, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
  if(dsq == NULL)
    cm_Fail("in CP9Viterbi, dsq is NULL.");
    
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
    otsc_k[cp9O_MEL]= cm->cp9->tsc[CTMEL][k];
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
	}

      endsc = ESL_MAX(endsc, dmx[cur][M] + TSC(cp9O_DM,M)); /* transition from D_M -> end */
      endsc = ESL_MAX(endsc, imx[cur][M] + TSC(cp9O_IM,M)); /* transition from I_M -> end */
      for(c = 0; c < cm->cp9->el_from_ct[M+1]; c++) /* el_from_ct[k] is >= 0 */
	/* transition penalty to EL incurred when EL was entered */
	endsc = ESL_MAX(endsc, elmx[cur][cm->cp9->el_from_idx[M+1][c]]);

      erow[cur] = endsc;
      scA[jp]   = endsc;
      fsc = Scorify(endsc);

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
  free(otsc);

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

  //printf("ctr: %d\n", ctr);
  return return_sc;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}


/***********************************************************************
 * Function: cp9_FastViterbiBackward()
 * 
 * Purpose:  Runs the Viterbi dynamic programming algorithm BACKWARDS on an
 *           input subsequence (i0-j0). 
 *           Somewhat flexible based on input options.
 *    
 * Note:     IDENTICAL to CP9Backward() below with maxes replacing
 *           sums in DP recursion, and the possibility of returning a
 *           CP9 trace if doing_align == TRUE. See CP9Backward() for more info,
 *           including more verbose 'Purpose' and description of arguments.
 *
 * Returns:  if(!do_scan) log P(S,tr|M)/P(S,tr|R), as a bit score
 *           else         max log P(S,tr|M)/P(S,tr|R), for argmax subseq S of input seq i0..j0,
 */
float
cp9_FastViterbiBackward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
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
    cm_Fail("in CP9Viterbi, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    cm_Fail("in CP9Viterbi, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    cm_Fail("in CP9Viterbi, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    cm_Fail("in CP9Viterbi, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
  if(dsq == NULL)
    cm_Fail("in CP9Viterbi, dsq is NULL.");
    
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
    otsc_k[cp9O_MEL]= cm->cp9->tsc[CTMEL][k];
    otsc_k[cp9O_ME] = cm->cp9->esc[k];
  }
  int const *tsc = otsc;
  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. We use it first for
   * the Forward scan, and reuse it for the Backward scan.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+2));
  gamma[0]     = 0.; /* this is impossible, no hit can be rooted at ip=0, 
		      * we require hits to have at least 1 emission,
		      * plus, no positive scoring hit can have 0 emissions as 
		      * all transition scores are negative */
  gamma[L]   = 0.;
  ESL_ALLOC(gback, sizeof(int)   * (L+2));
  gback[L]   = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+2));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 2;
  else             nrows = L+1;
  mx = AllocCPlan9Matrix(nrows, M, &mmx, &imx, &dmx, &elmx, &erow); 

  /* scA will hold P(seq up to j | Model) in int log odds form */
  ESL_ALLOC(scA, sizeof(int) * (j0-i0+2));
			
  /* Initialization of the L (i = j0, cur = (j0-i) = (j0-j0) %2 = 0) row. */
  if(be_efficient) cur = 0;
  else cur = j0-i0+1; /* L */
  ip = j0-i0+1;  /*L */
  i  = j0;

  /*******************************************************************
   * 0 Handle EL, looking at EL_k->E for all valid k.
   * we're going backwards so we have to work out of order, we could get 
   * around this by storing the nodes each EL goes TO in an el_to_ct[] vec. */
  /* init to -INFTY */
  for (k = 1; k <= cm->cp9->M; k++) elmx[cur][k] = -INFTY;
  if(cm->cp9->flags & CPLAN9_EL) {
    for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M+1]; c++) /* el_from_ct[cm->cp9->M+1] holds # ELs that can go to END */
      elmx[cur][cm->cp9->el_from_idx[cm->cp9->M+1][c]] = 0.; /* EL<-E, penalty incurred when we enter EL (i.e. leave going backwards) */
  }
  /*******************************************************************/

  /* elmx[cur][cm->cp9->M] is either 0 (if EL_M exists (it would nec be in el_from_idx[cm->cp9->M+1] array if it does, so
   * it would be filled with 0 in above loop), or -INFTY if it doesn't exist. We don't add possibility of EL_M -> EL_M
   * self loop b/c it's impossible to do that without emitting, and we've already seen our last res emitted. 
   * either way we don't have to modify it */

  mmx[cur][cm->cp9->M]  = ESL_MAX(elmx[cur][cm->cp9->M] + cm->cp9->tsc[CTMEL][cm->cp9->M],/* M_M<-EL_M<-E, with 0 self loops in EL_M */
				  cm->cp9->esc[cm->cp9->M]);                              /* M_M<-E ... everything ends in E (the 0; 2^0=1.0) */
  mmx[cur][cm->cp9->M] += cm->cp9->msc[dsq[i]][cm->cp9->M];  /* ... + emitted match symbol */
  imx[cur][cm->cp9->M]  = cm->cp9->tsc[CTIM][cm->cp9->M];    /* I_M<-E ... everything ends in E (the 0; 2^0=1.0) */
  imx[cur][cm->cp9->M] += cm->cp9->isc[dsq[i]][cm->cp9->M];  /* ... + emitted insert symbol */
  dmx[cur][cm->cp9->M]  = cm->cp9->tsc[CTDM][cm->cp9->M];    /* D_M<-E */

  /*******************************************************************
   * No need to look at EL_k->M_M b/c elmx[cur] with cur == L means last emitted residue was L+1 
   * and this is impossible if we've come from M_M (only would be valid if we were coming from
   * E which is handled above with the EL_k->E code). 
   *******************************************************************/
  for (k = cm->cp9->M-1; k >= 1; k--)
    {
      mmx[cur][k]  = cm->cp9->esc[k];  /*M_k<- E */
      mmx[cur][k]  = ESL_MAX(mmx[cur][k], dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]);
      if(cm->cp9->flags & CPLAN9_EL)
	mmx[cur][k]  = ESL_MAX(mmx[cur][k], elmx[cur][k] + cm->cp9->tsc[CTMEL][k]);
      mmx[cur][k] += cm->cp9->msc[dsq[i]][k];

      /*******************************************************************
       * No need to look at EL_k->M_M b/c elmx[cur] with cur == L means last emitted residue was L+1 
       * and this is impossible if we've come from M_M (only would be valid if we were coming from
       * E which is handled above with the EL_k->E code). 
       *******************************************************************/
      imx[cur][k]  = dmx[cur][k+1] + cm->cp9->tsc[CTID][k];
      imx[cur][k] += cm->cp9->isc[dsq[i]][k];

      dmx[cur][k]  = dmx[cur][k+1] + cm->cp9->tsc[CTDD][k];
      /* elmx[cur][k] was set above, out of order */
    }
  
  /* remember M_0 is special, the B state, a non-emitter */
  mmx[cur][0]  = dmx[cur][1] + cm->cp9->tsc[CTMD][0]; /* M_0(B)->D_1, no seq emitted, all deletes */
  /* above line is diff from CPBackwardOLD() which has mmx[cur][0] = -INFTY; */
  imx[cur][0]  = dmx[cur][1] + cm->cp9->tsc[CTID][0];
  imx[cur][0] += cm->cp9->isc[dsq[i]][0];

  dmx[cur][0]   = -INFTY; /*D_0 doesn't exist*/
  elmx[cur][0]  = -INFTY; /*EL_0 doesn't exist*/

  scA[ip] = mmx[cur][0]; /* all parses must start in M_0, the B state */
  fsc = Scorify(scA[ip]);
  /*printf("ip: %d i: %d fsc: %f i: %d\n", ip, i, fsc, scA[ip]);*/
  /* we can't have a hit starting here, b/c it would correspond to all deletes,
   * no seq emitted, so we skip the check for if(fsc > best_sc) */

  /*printf("scA[ip:%d]: %d\n", ip, scA[ip]);*/

  /*printf("mmx[ip:%d][%d]: %d cur: %d\n", L+1, 0, mmx[cur][0], cur);
    printf("imx[ip:%d][%d]: %d cur: %d\n", L+1, 0, imx[cur][0], cur);
    printf("dmx[ip:%d][%d]: %d cur: %d\n", L+1, 0, dmx[cur][0], cur);*/
  
  /*****************************************************************
   * The main loop: scan the sequence from position j0-1 to i0.
   *****************************************************************/
  /* Recursion */
  for (i = j0-1; i >= (i0-1); i--)
    {
      ip = i-i0+1;		/* ip is relative index in dsq (0 to L-1) */
      if(be_efficient) { cur = (j0-i) % 2;  prv = (j0-i+1) % 2; }
      else             { cur = ip;          prv = ip+1;         }
      /* init EL mx to -INFTY */
      for (k = 1; k <= cm->cp9->M; k++) elmx[cur][k] = -INFTY;
      if(ip > 0) {
	/* elmx[cur][k] is possibly of coming from self (EL_k), we 
	 * can't have come from END b/c we haven't emitted the last res of the seq yet.
	 */
	if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[cm->cp9->M]))
	  elmx[cur][cm->cp9->M] = elmx[cur][cm->cp9->M] + cm->cp9->el_selfsc;
	
	mmx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTMI][cm->cp9->M];
	mmx[cur][cm->cp9->M] += cm->cp9->msc[dsq[i]][cm->cp9->M];
	
	if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[cm->cp9->M]))
	  mmx[cur][cm->cp9->M] = ESL_MAX(mmx[cur][cm->cp9->M], 
					 elmx[cur][cm->cp9->M] + cm->cp9->tsc[CTMEL][cm->cp9->M]);
	
	imx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTII][cm->cp9->M];
	imx[cur][cm->cp9->M] += cm->cp9->isc[dsq[i]][cm->cp9->M];
      }
      else { /* ip == 0 */
	mmx[cur][cm->cp9->M] = -INFTY;  /* need seq to get here */
	imx[cur][cm->cp9->M] = -INFTY;  /* need seq to get here */
	elmx[cur][cm->cp9->M]= -INFTY;  /* first emitted res can't be from an EL, need to see >= 1 matches */
      }
      dmx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTDI][cm->cp9->M]; 

      /*******************************************************************
       * 1b Handle EL, looking at EL_k->M_M for all valid k.
       * EL_k->M_M transition, which has no transition penalty */
      if(cm->cp9->flags & CPLAN9_EL)
	{
	  for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M]; c++) /* el_from_ct[cm->cp9->M] holds # ELs that can go to M_M */
	    elmx[cur][cm->cp9->el_from_idx[cm->cp9->M][c]] = ESL_MAX(elmx[cur][cm->cp9->el_from_idx[cm->cp9->M][c]], mmx[prv][cm->cp9->M]);
	}

      /* A main difference between a Backward scanner and 
       * regular Backward: a scanner can end at the END 
       * state at any position, regular can only end at
       * the final position j0. */
      if(do_scan) {	
	if(ip > 0) {
	  /*******************************************************************
	   * 2 Handle EL, looking at EL_k->E for all valid k.
	   * EL_k->M_M transition, which has no transition penalty */
	  if(cm->cp9->flags & CPLAN9_EL) {
	    for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M+1]; c++) /* el_from_ct[cm->cp9->M] holds # ELs that can go to END */
	      elmx[cur][cm->cp9->el_from_idx[cm->cp9->M+1][c]] = 0.; /* EL<-E, penalty incurred when we enter EL (i.e. leave going backwards) */
	  }
	  /*******************************************************************/
	  /* elmx[cur][cm->cp9->M] is either 0 (if EL_M exists (it would nec be in el_from_idx[cm->cp9->M+1] array if it does, so
	   * it would be filled with 0 in above loop), or -INFTY if it doesn't exist. We don't add possibility of EL_M -> EL_M
	   * self loop b/c it's impossible to do that without emitting, and we've already seen our last res emitted,
	   * either way we don't have to modify it */
	  
	  mmx[cur][cm->cp9->M] = ESL_MAX(mmx[cur][cm->cp9->M], elmx[cur][cm->cp9->M] + cm->cp9->tsc[CTMEL][cm->cp9->M]);/* M_M<-EL_M<-E, with 0 selfs in EL_M */
	  mmx[cur][cm->cp9->M] = ESL_MAX(mmx[cur][cm->cp9->M], cm->cp9->esc[cm->cp9->M]);                              /* M_M<-E ... */
	  ///mmx[cur][cm->cp9->M] += cm->cp9->msc[dsq[i]][cm->cp9->M]; /* ... + emitted match symbol */
	  /* DO NOT add contribution of emitting i from M, it's been added above */
	  
	  imx[cur][cm->cp9->M] = ESL_MAX(imx[cur][cm->cp9->M],
					 (cm->cp9->tsc[CTIM][cm->cp9->M] +            /* I_M<-E + (only in scanner)     */
					  0));                                        /* all parses end in E, 2^0 = 1.0;*/
	  ///imx[cur][cm->cp9->M] += cm->cp9->isc[dsq[i]][cm->cp9->M]; /* ... + emitted insert symbol */  
	  /* DO NOT add contribution of emitting i from M, it's been added above */
	}
	dmx[cur][cm->cp9->M] =  ESL_MAX(dmx[cur][cm->cp9->M], 
					(cm->cp9->tsc[CTDM][cm->cp9->M] +            /* D_M<-E + (only in scanner)     */
					 0));                                        /* all parses end in E, 2^0 = 1.0;*/
      }
      /*printf("mmx[ip:%d][%d]: %d cur: %d\n", ip, cm->cp9->M, mmx[cur][cm->cp9->M], cur);
	printf("imx[ip:%d][%d]: %d cur: %d\n", ip, cm->cp9->M, imx[cur][cm->cp9->M], cur);
	printf("dmx[ip:%d][%d]: %d cur: %d\n", ip, cm->cp9->M, dmx[cur][cm->cp9->M], cur);*/
      
      for (k = cm->cp9->M-1; k >= 1; k--) {
	if(ip > 0) {
	  /*******************************************************************
	   * 3 Handle EL, looking at EL_k->M_k for all valid k and EL_k->EL_k
	   * we're going backwards so we have to work out of order
	   * we could get around this by storing the nodes each EL goes TO
	   * in an el_to_ct[] vector. */
	  if(cm->cp9->flags & CPLAN9_EL) {
	    for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] holds # ELs that can go to M_k */
	      elmx[cur][cm->cp9->el_from_idx[k][c]] = ESL_MAX(elmx[cur][cm->cp9->el_from_idx[k][c]], mmx[prv][k]);
	    /* EL<-M, penalty incurred when we enter EL (i.e. leave going backwards) */
	  }
	  /*******************************************************************/
	  
	  /* Finish off elmx[cur][k] with possibility of coming from self (EL_k), 
	   * elmx[cur][k] will have been filled by block above for ks > current k,
	   * no M_k -> EL_k' with k' > k */
	  if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[k]))
	    elmx[cur][k] = ESL_MAX(elmx[cur][k], elmx[prv][k] + cm->cp9->el_selfsc);
	  
	  mmx[cur][k]  = ESL_MAX(ESL_MAX((mmx[prv][k+1] + cm->cp9->tsc[CTMM][k]),  
					 (imx[prv][k]   + cm->cp9->tsc[CTMI][k])),
				 (dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]));
	  if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[k]))
		mmx[cur][k] = ESL_MAX(mmx[cur][k], elmx[cur][k] + cm->cp9->tsc[CTMEL][k]); /* penalty for entering EL */
	      mmx[cur][k] += cm->cp9->msc[dsq[i]][k];

	      imx[cur][k]  = ESL_MAX(ESL_MAX((mmx[prv][k+1] + cm->cp9->tsc[CTIM][k]),
					     (imx[prv][k]   + cm->cp9->tsc[CTII][k])),
				     (dmx[cur][k+1] + cm->cp9->tsc[CTID][k]));
	      imx[cur][k] += cm->cp9->isc[dsq[i]][k];

	}
	else { /* ip == 0 */
	  mmx[cur][k] = -INFTY; /* need seq to get here, unless we come from E in a scanner (below) */
	  imx[cur][k] = -INFTY; /* need seq to get here */
	  elmx[cur][k]= -INFTY;  /* first emitted res can't be from an EL, need to see >= 1 matches */
	}
	if(do_scan && ip > 0) { /* add possibility of ending at this position from this state */
	  mmx[cur][k] = 
	    ESL_MAX(mmx[cur][k], 
		    (cm->cp9->esc[k] +                    /* M_k<-E + (only in scanner)     */ 
		     0));                                 /* all parses end in E, 2^0 = 1.0;*/
	  /* DO NOT add contribution of emitting i from M, it's been added above */
	  /* No EL contribution here b/c we'd be looking for M_k<-EL_k<-E, but EL_k<-E is impossible 
	   * for k != cm->cp9->M; */
	}	      
	dmx[cur][k]  = ESL_MAX(ESL_MAX((mmx[prv][k+1] + cm->cp9->tsc[CTDM][k]),
				       (imx[prv][k]   + cm->cp9->tsc[CTDI][k])),
			       (dmx[cur][k+1] + cm->cp9->tsc[CTDD][k]));
      } /* end of for (k = cm->cp9->M-1; k >= 1; k--) { */
      /* Case when k == 0 */
      /* imx[cur][0] is filled same as imx[cur][1..k] in the loop above */
      if(ip > 0) {
	imx[cur][0] = ESL_MAX(ESL_MAX((mmx[prv][1] + cm->cp9->tsc[CTIM][0]),
				      (imx[prv][0] + cm->cp9->tsc[CTII][0])),
			      (dmx[cur][1] + cm->cp9->tsc[CTID][k]));
	imx[cur][0] += cm->cp9->isc[dsq[i]][k];
      }
      else /* ip == 0 */
	imx[cur][0] = -INFTY; /* need seq to get here */
      dmx[cur][0]   = -INFTY; /* D_0 does not exist */
      elmx[cur][0]  = -INFTY; /* EL_0 does not exist */

      /*M_0 is the B state, it doesn't emit, and can be reached from any match via a begin transition */
      mmx[cur][0] = -INFTY;
      for (k = cm->cp9->M; k >= 1; k--) mmx[cur][0] = ESL_MAX(mmx[cur][0], (mmx[prv][k] + cm->cp9->bsc[k]));
      mmx[cur][0] = ESL_MAX(mmx[cur][0], (imx[prv][0] + cm->cp9->tsc[CTMI][0]));
      mmx[cur][0] = ESL_MAX(mmx[cur][0], (dmx[cur][1] + cm->cp9->tsc[CTMD][0]));     /* B->D_1 */
      /* No EL contribution here, can't go B->EL_* */
      
      /* determine isc, the int score of all possible parses starting at the current
       * position (i) of the target sequence. */
      scA[ip] = mmx[cur][0]; /* all parses must start in M_0, the B state */
      fsc = Scorify(scA[ip]);
      /*printf("ip: %d i: %d fsc: %f i: %d\n", ip, i, fsc, sc[ip]);*/
      if(fsc > best_sc) { best_sc = fsc; best_pos= i+1; } /* *off-by-one* (see *off-by-one* below) */
      if(!(cm->search_opts & CM_SEARCH_HMMGREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   * *off-by-one*:
	   * There's an off-by-one issue here: all Backward hits are rooted in 
	   * M_O, the B (begin) state, which is a non-emitter.
	   * let i = ip+i0-1 => ip = i-i0+1;
	   * so sc[ip] = backward->mmx[ip][0] = summed log prob of all parses that end at j0, 
	   * and start at position i+1 of the sequence (because i+1 is the last residue
	   * whose emission has been accounted for). As a result, gamma indexing is off-by-one
	   * with respect to sequence position, hence the i+1 or i-1 in the following
	   * code blocks, each marked by "*off-by-one*" comment below. 
	   * for example: let i0 = 2 gamma[ip=4], normally this means ip=4 corresponds to i=5 
	   *              but due to this off-by-one sc[ip=4] corresponds to hits that start at i=6
	   */
	  gamma[ip]  = gamma[ip+1] + 0; /* extend without adding a new hit */
	  /*printf("ip: %d | gamma[ip]: %f | gamma[ip+1]: %f\n", ip, gamma[ip], gamma[ip+1]);*/
	  gback[ip]  = -1;
	  savesc[ip] = IMPOSSIBLE;
	  j = (((i+1)+W-1) < j0) ? ((i+1)+W-1) : j0; /* *off-by-one* */
	  jp = j-i0+1;
	  curr_sc = gamma[jp+1-1] + fsc + cm->cp9_sc_boost; /* *off-by-one* */
	  /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
	   * value is 0.0 if technique not used. */
	  if (curr_sc > gamma[ip])
	    {
	      gamma[ip]  = curr_sc;
	      /*printf("\ti: %d | gamma[i]: %f\n", i+1, gamma[ip]);*/
	      gback[ip]  = j;
	      savesc[ip] = fsc;
	    }
	  /*printf("i: %d ip: %d gamma[ip]: %f\n", i, ip, gamma[ip]);*/
	}
      else
	{
	  /* Resolving overlaps greedily (RSEARCH style),  
	   * Return best hit for each j, IFF it's above threshold */
	  if (fsc >= cutoff) 
	    {
	      if(results != NULL) 
		{
		  j = (((i+1)+W-1) < j0) ? ((i+1)+W-1) : j0; /* *off-by-one* */
		  /*printf("BWD greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i+1, j, fsc);*/ /* *off-by-one* */
		  report_hit (i+1, j, 0, fsc, results); 
		  /* 0 is for saver, which is irrelevant for HMM hits */
		}
	    }
	}
      
      if (fsc > best_hmm_sc) { best_hmm_sc = fsc; best_hmm_pos= i; }
    }
  /* End of Backward recursion */
  
  if((!(cm->search_opts & CM_SEARCH_HMMGREEDY)) && /* resolve overlaps optimally */
     (!doing_align || do_scan)) /* else we can save time by skipping traceback */
    {
      /*****************************************************************
       * Traceback stage for Backward.
       * Recover all hits: an (i,j,sc) triple for each one.
       * Remember the off-by-one b/t seq index and gamma (see *off-by-one* above)
       *****************************************************************/ 
      i     = i0; 
      while (i <= j0) 
	{
	  ip    = (i-1)-i0+1; /* *off-by-one*, i-1, ip corresponds to i+1 
			       * (yes: i *-* 1 =>   ip corresponds to i *+* 1) */
	  if (gback[ip] == -1) /* no hit */
	    i++; 
	  else                /* a hit, a palpable hit */
	    {
	      if(savesc[ip] > best_hmm_sc) 
		{
		  best_hmm_sc = savesc[ip];
		  best_hmm_pos= i;
		}
	      if(savesc[ip] >= cutoff)
		{
		  if(results != NULL) /* report the hit */
		    {
		      /*printf("BWD optimal reporting hit i: %d j: %d sc: %f\n", i, gback[ip], savesc[ip]);*/
		      report_hit(i, gback[ip], 0, savesc[ip], results); 
		      /* 0 is for saver, which is irrelevant for HMM hits */
		    }
		}
	      i = gback[ip]+1;
	    }
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  /*printf("returning from CP9Backward()\n");*/
  if(ret_sc != NULL) *ret_sc = scA;
  else free(scA);
  /* determine score to return: (I know, too complex) */
  if(doing_align)
    {
      return_sc = Scorify(mmx[cur][0]);
      if(ret_bestpos != NULL) *ret_bestpos = i0;
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
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);
  /*printf("Backward return_sc: %f\n", return_sc);*/
  return return_sc;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}


/***********************************************************************
 * Function: cp9_FastForward()
 * 
 * Purpose:  Runs the Forward dynamic programming algorithm on an
 *           input subsequence (i0-j0). Complements cp9_FastBackward().  
 *           Somewhat flexible based on input options as follows:
 *
 *           if(be_efficient): only allocates 2 rows of the Forward
 *           matrix, else allocates full L+1 matrix.
 *
 *           if(do_scan): allows parses to start at any position i
 *           i0..j0, changing meaning of DP matrix cells as discussed
 *           below.
 *
 *           Reference for algorithm (when do_scan is FALSE): 
 *           Durbin et. al. Biological Sequence Analysis; p. 58.
 *
 *           The meaning of the Forward (F) matrix DP cells for
 *           matches (M) inserts (I) and deletes (D):
 *
 *           For relative subsequence positions ip = 0..L:
 *           For HMM nodes 1..M: 
 *           F->M[ip][k] : sum of all parses emitting seq
 *                         from i0..ip that visit node k's match 
 *                         state, which emits posn ip
 *           F->I[ip][k] : sum of all parses emitting seq from 
 *                         i0..ip that visit node k's insert
 *                         state, which emits posn ip 
 *           F->D[ip][k] : sum of all parses emitting seq from 
 *                         i0..ip that visit node k's delete
 *                         delete state, last emitted (leftmost)
 *                         posn was ip
 *
 *           For *special* HMM node 0:
 *           F->M[ip][0] : M_0 is the Begin state, which does not 
 *                         emit, so this is the sum of all parses 
 *                         emitting seq from i0..ip that start
 *                         in the begin state, the last emitted 
 *                         (leftmost) posn was ip.
 *
 *           Note: if ip=0, only D_k and M_0 states can have 
 *                 non-IMPOSSIBLE values. 
 *
 *           if(do_scan) the 'i0..ip' in the above definitions is
 *           changed to iE..ip such that i0 <= iE <= ip. Meaning
 *           any residue can be the first residue emitted in the
 *           parse. This means F->M[ip][0] is the sum of all parses
 *           emitting a subseq starting anywhere from i0..ip and 
 *           ending at ip. 
 *
 *
 * Args
 *           cm        - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq       - sequence in digitized form
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - the maximum size of a hit (often cm->W)
 *           cutoff    - minimum score to report
 *           ret_sc    - RETURN: int log odds Forward score for each end point [0..(j0-i0+1)]
 *           ret_maxres- RETURN: start position that gives maximum score max argmax_i sc[i]
 *           results   - search_results_t to add to; if NULL, don't keep results
 *           do_scan   - TRUE if we're scanning, HMM can start to emit anywhere i0..j0,
 *                       FALSE if we're not, HMM must start emitting at i0, end emitting at j0
 *           doing_align  - TRUE if reason we've called this function is to help get posteriors
 *                          for CM alignment, in which case we can skip the traceback. 
 *           doing_rescan - TRUE if we've called this function recursively, and we don't
 *                          want to again
 *           be_efficient- TRUE to keep only 2 rows of DP matrix in memory, FALSE keep whole thing
 *           ret_mx    - RETURN: CP9 Forward DP matrix, NULL if not wanted
 *
 * Returns:  if(!do_scan) log P(S|M)/P(S|R), as a bit score
 *           else         max log P(S|M)/P(S|R), for argmax subseq S of input seq i0..j0,
 */
float
cp9_FastForward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
		int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
		int doing_rescan, int be_efficient, CP9_dpmatrix_t **ret_mx)
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
  int          nrows = 2;   /* number of rows for the dp matrix */
  int          c;           /* counter for EL states */
  float        temp_sc;     /* temporary score                                              */
  int          M;
  int          accept;
  /*debug_print_cp9_params(stdout, cm->cp9, TRUE);*/
  /*printf("in cp9_FastForward() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    cm_Fail("in cp9_FastForward, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    cm_Fail("in cp9_FastForward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(doing_rescan && !do_scan) 
    cm_Fail("in cp9_FastForward, doing_rescan but not do_scan");
  if(results != NULL && !do_scan)
    cm_Fail("in cp9_FastForward, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    cm_Fail("in cp9_FastForward, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
  if(dsq == NULL)
    cm_Fail("in cp9_FastForward, dsq is NULL.");
    
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
    otsc_k[cp9O_MEL]= cm->cp9->tsc[CTMEL][k];
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

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  int sc;
  for (k = 1; k <= M; k++)
    {
      mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
      sc = ILogsum(ILogsum(mmx[0][k-1] + TSC(cp9O_MD,k-1),
			   imx[0][k-1] + TSC(cp9O_ID,k-1)),
		   dmx[0][k-1] + TSC(cp9O_DD,k-1));
      dmx[0][k] = sc;
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

  /* int ctr = 0; */
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

      sc = ILogsum(ILogsum(mmx[prv][0] + TSC(cp9O_MI,0),
			   imx[prv][0] + TSC(cp9O_II,0)),
		   dmx[prv][0] + TSC(cp9O_DI,0));
      imx[cur][0] = sc + isc[0];

      for (k = 1; k <= M; k++)
	{
	  /*match state*/
	  //ctr++;
	  sc = ILogsum(ILogsum(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
			       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
		       ILogsum(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
			       mmx[prv][0]   + TSC(cp9O_BM,k  )));
	  /* check possibility we came from an EL, if they're valid */
	  for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	    sc = ILogsum(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
	    /* transition penalty to EL incurred when EL was entered */
	  mmx[cur][k] = sc + msc[k];

	  /* E state update */
	  endsc = ILogsum(endsc, mmx[cur][k] + TSC(cp9O_ME,k));

	  /*insert state*/
	  sc = ILogsum(ILogsum(mmx[prv][k] + TSC(cp9O_MI,k),
			       imx[prv][k] + TSC(cp9O_II,k)),
		       dmx[prv][k] + TSC(cp9O_DI,k));
	  imx[cur][k] = sc + isc[k];

	  /*delete state*/
	  sc = ILogsum(ILogsum(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
			       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
		       dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	  dmx[cur][k] = sc;

	  /*el state*/
	  sc = -INFTY;
	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								    HMM nodes that map to right half of a MATP_MP) */
	    {
	      sc = ILogsum(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
			   elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	    }
	  elmx[cur][k] = sc;
	  /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
	    printf("imx [jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
	    printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);
	    printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[cur][k]);*/
	}
      endsc = ILogsum(ILogsum(endsc, dmx[cur][M] + TSC(cp9O_DM,M)), /* transition from D_M -> end */
		      imx[cur][M] + TSC(cp9O_IM,M)); /* transition from I_M -> end */
      for(c = 0; c < cm->cp9->el_from_ct[M+1]; c++) /* el_from_ct[k] is >= 0 */
	endsc = ILogsum(endsc, elmx[cur][cm->cp9->el_from_idx[M+1][c]]);
	/* transition penalty to EL incurred when EL was entered */
      /*printf("endsc: %d\n", endsc);*/

      erow[cur] = endsc;
      scA[jp]   = endsc;
      fsc = Scorify(endsc);
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
		accept = TRUE;
		/* Potentially rescan just the subseq that is the hit we're about to report.
		 * Implemented to deal with fact that --enfseq option was enforcing the subseq
		 * to have hit pass filter b/c this Forward scanning function is 'infinite' length
		 * (Weinberg). Sometimes the subseq we're about to report has a really
		 * crappy score, even though the cumulative Forward score (starting at i0) is good. */
		if(do_scan && cm->search_opts & CM_SEARCH_HMMRESCAN && doing_rescan == FALSE)
		  {
		    /*printf("rechecking hit from %d to %d\n", gback[jp], j);*/
		    temp_sc = CP9Forward(cm, dsq, gback[jp], j, cm->W, cutoff, 
					 NULL,  /* don't care about scores of each pos */
					 NULL,  /* don't care about best scoring position */
					 NULL,  /* don't report hits to results data structure */
					 TRUE,  /* we're scanning */
					 FALSE, /* we're not ultimately aligning */
					 TRUE,  /* set the doing_rescan arg to TRUE, 
						   so we don't potentially infinitely recurse */
					 TRUE,  /* be memory efficient */
					 NULL); /* don't want the DP matrix back */
		    /*printf("new score: %f old score %f\n", temp_sc, savesc[jp]);*/
		    if(temp_sc >= cutoff) 
		      { 
			accept = TRUE; 
			/*printf("rechecked hit from %d to %d\n", gback[jp], j);
			  printf("new score: %f old score %f\n", temp_sc, savesc[jp]);*/
		      }
		    else accept = FALSE;
		    savesc[jp] = temp_sc;
		  }
		if(accept)
		  {
		    if(results != NULL) /* report the hit */
		      {
			/*printf("FWD reporting hit: i: %d j: %d sc: %f\n", gback[jp], j, savesc[jp]);*/
			report_hit(gback[jp], j, 0, savesc[jp], results); 
			/* 0 is for saver, which is irrelevant for HMM hits */
		      }
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
  free(otsc);

  /* determine score to return: (I know, too complex) */
  if(doing_align)
    {
      return_sc = Scorify(scA[(j0-i0+1)]); /* L = j0-i0+1 */
      if(ret_bestpos != NULL) *ret_bestpos = i0;
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

  //printf("ctr: %d\n", ctr);
  return return_sc;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}


/***********************************************************************
 * Function: cp9_EXTPLFastForward()
 * 
 * Purpose:  Runs the Forward dynamic programming algorithm on an
 *           input subsequence (i0-j0). Complements cp9_FastBackward().  
 *           Somewhat flexible based on input options as follows:
 *
 *           if(be_efficient): only allocates 2 rows of the Forward
 *           matrix, else allocates full L+1 matrix.
 *
 *           if(do_scan): allows parses to start at any position i
 *           i0..j0, changing meaning of DP matrix cells as discussed
 *           below.
 *
 *           Reference for algorithm (when do_scan is FALSE): 
 *           Durbin et. al. Biological Sequence Analysis; p. 58.
 *
 *           The meaning of the Forward (F) matrix DP cells for
 *           matches (M) inserts (I) and deletes (D):
 *
 *           For relative subsequence positions ip = 0..L:
 *           For HMM nodes 1..M: 
 *           F->M[ip][k] : sum of all parses emitting seq
 *                         from i0..ip that visit node k's match 
 *                         state, which emits posn ip
 *           F->I[ip][k] : sum of all parses emitting seq from 
 *                         i0..ip that visit node k's insert
 *                         state, which emits posn ip 
 *           F->D[ip][k] : sum of all parses emitting seq from 
 *                         i0..ip that visit node k's delete
 *                         delete state, last emitted (leftmost)
 *                         posn was ip
 *
 *           For *special* HMM node 0:
 *           F->M[ip][0] : M_0 is the Begin state, which does not 
 *                         emit, so this is the sum of all parses 
 *                         emitting seq from i0..ip that start
 *                         in the begin state, the last emitted 
 *                         (leftmost) posn was ip.
 *
 *           Note: if ip=0, only D_k and M_0 states can have 
 *                 non-IMPOSSIBLE values. 
 *
 *           if(do_scan) the 'i0..ip' in the above definitions is
 *           changed to iE..ip such that i0 <= iE <= ip. Meaning
 *           any residue can be the first residue emitted in the
 *           parse. This means F->M[ip][0] is the sum of all parses
 *           emitting a subseq starting anywhere from i0..ip and 
 *           ending at ip. 
 *
 *
 * Args
 *           cm        - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq       - sequence in digitized form
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - the maximum size of a hit (often cm->W)
 *           cutoff    - minimum score to report
 *           ret_sc    - RETURN: int log odds Forward score for each end point [0..(j0-i0+1)]
 *           ret_maxres- RETURN: start position that gives maximum score max argmax_i sc[i]
 *           results   - search_results_t to add to; if NULL, don't keep results
 *           do_scan   - TRUE if we're scanning, HMM can start to emit anywhere i0..j0,
 *                       FALSE if we're not, HMM must start emitting at i0, end emitting at j0
 *           doing_align  - TRUE if reason we've called this function is to help get posteriors
 *                          for CM alignment, in which case we can skip the traceback. 
 *           doing_rescan - TRUE if we've called this function recursively, and we don't
 *                          want to again
 *           be_efficient- TRUE to keep only 2 rows of DP matrix in memory, FALSE keep whole thing
 *           be_safe     - TRUE to never use ILogsumNI(sc1,sc2) (which assumes sc1, sc2 != -INFTY)
 *           ret_mx    - RETURN: CP9 Forward DP matrix, NULL if not wanted
 *        
 *
 * Returns:  if(!do_scan) log P(S|M)/P(S|R), as a bit score
 *           else         max log P(S|M)/P(S|R), for argmax subseq S of input seq i0..j0,
 */
float
cp9_EXPTLFastForward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
		     int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, int doing_rescan,
		     int be_efficient, int be_safe, CP9_dpmatrix_t **ret_mx)
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
  int          accept;      /* flag used if cm->search_opts & CM_SEARCH_HMMRESCAN           */
  float        temp_sc;     /* temporary score                                              */
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */
  int          c;           /* counter for EL states                                        */
  int          M;           /* cm->cp9->M, number of nodes in the model                     */
  int          locality_mode; /* 1 of 4 locality modes, we can optimize (a bit) specifically*
			       * for each of the 4 modes                                    */
  /*debug_print_cp9_params(stdout, cm->cp9, TRUE);*/

  printf("in cp9_EXPTLFastForward() i0: %d j0: %d\n", i0, j0);  
  /* Contract checks */
  if(cm->cp9 == NULL)
    cm_Fail("in cp9_FastForward, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    cm_Fail("in cp9_FastForward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(doing_rescan && !do_scan) 
    cm_Fail("in cp9_FastForward, doing_rescan but not do_scan");
  if(results != NULL && !do_scan)
    cm_Fail("in cp9_FastForward, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    cm_Fail("in cp9_FastForward, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
  if((cm->cp9->flags & CPLAN9_LOCAL_BEGIN)     && (! (cm->cp9->flags & CPLAN9_LOCAL_BEGIN)))
    cm_Fail("ERROR in cp9_FastForward, CPLAN9_LOCAL_BEGIN flag up, but CPLAN9_LOCAL_END flag down, this shouldn't happen.");
  if((! (cm->cp9->flags & CPLAN9_LOCAL_BEGIN)) && (cm->cp9->flags & CPLAN9_LOCAL_BEGIN))
    cm_Fail("ERROR in cp9_FastForward, CPLAN9_LOCAL_BEGIN flag down, but CPLAN9_LOCAL_END flag up, this shouldn't happen.");
  if(dsq == NULL)
    cm_Fail("in cp9_FastForward, dsq is NULL.");
    
  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_hmm_sc = IMPOSSIBLE;
  best_hmm_pos= -1;
  L = j0-i0+1;
  M = cm->cp9->M;
  if (W > L) W = L; 

  /* Transition scores, we get a little speed win by reordering them here */
  int *otsc;
  ESL_ALLOC(otsc,      sizeof(int)   * (M+1)  * cp9O_NTRANS);
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
    otsc_k[cp9O_MEL]= cm->cp9->tsc[CTMEL][k];
    otsc_k[cp9O_ME] = cm->cp9->esc[k];
  }
  int const *tsc = otsc;

  /* determine and set locality mode, and check that transition guarantees hold */
  locality_mode = cp9_GetLocalityMode(cm->cp9);
  /* printf("locality mode: %d\n", locality_mode); */
  /* printf("0 be_safe: %d\n", be_safe); */
  if((status = cp9_CheckTransitionGuarantees(cm->cp9)) != eslOK) be_safe = TRUE;
  /* printf("1 be_safe: %d\n", be_safe); */

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

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  int sc;
  for (k = 1; k <= M; k++)
    {
      mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
      sc = ILogsum(ILogsum(mmx[0][k-1] + TSC(cp9O_MD,k-1),
			   imx[0][k-1] + TSC(cp9O_ID,k-1)),
		   dmx[0][k-1] + TSC(cp9O_DD,k-1));
      dmx[0][k] = sc;
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

  /* int ctr = 0; */

  /*****************************************************************************
   * special position: j = i0                                                  */
  j  = i0;
  int const *isc = cm->cp9->isc[dsq[j]];
  int const *msc = cm->cp9->msc[dsq[j]];
  int endsc     = -INFTY;
  int el_selfsc = cm->cp9->el_selfsc;
  jp = j-i0+1;     /* jp is relative position in the sequence 1..L */
  cur = (j-i0+1);
  prv = (j-i0);
  if(be_efficient)  { cur %= 2; prv %= 2; }	  
  mmx[cur][0]        = (do_scan == TRUE) ? 0     : -INFTY;
  dmx[cur][0]        = -INFTY;  /*D_0 is non-existent*/
  elmx[cur][0]       = -INFTY;  /*no EL state for node 0 */
  sc = ILogsum(ILogsum(mmx[prv][0] + TSC(cp9O_MI,0),
		       imx[prv][0] + TSC(cp9O_II,0)),
	       dmx[prv][0] + TSC(cp9O_DI,0));
  imx[cur][0] = sc + isc[0];
  for (k = 1; k <= M; k++)
    {
      /*match state*/
      sc = ILogsum(ILogsum(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
			   imx[prv][k-1] + TSC(cp9O_IM,k-1)),
		   ILogsum(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
			   mmx[prv][0]   + TSC(cp9O_BM,k  )));
      mmx[cur][k] = sc + msc[k];

      /* E state update, put me inside mmx[cur][k] calc */
      endsc = ILogsum(endsc, mmx[cur][k] + TSC(cp9O_ME,k));

      /*insert state*/
      sc = ILogsum(ILogsum(mmx[prv][k] + TSC(cp9O_MI,k),
			   imx[prv][k] + TSC(cp9O_II,k)),
		   dmx[prv][k] + TSC(cp9O_DI,k));
      imx[cur][k] = sc + isc[k];

      /*delete state*/
      sc = ILogsum(ILogsum(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
			   imx[cur][k-1] + TSC(cp9O_ID,k-1)),
		   dmx[cur][k-1] + TSC(cp9O_DD,k-1));
      dmx[cur][k] = sc;

      /*el state*/
      sc = -INFTY;
      if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) { /* not all HMM nodes have an EL state (for ex: 
								  HMM nodes that map to right half of a MATP_MP) */
	sc = ILogsum(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
		     elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
      }
      elmx[cur][k] = sc;
      /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
	printf("imx [jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
	printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);
	printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[cur][k]);*/
    }
  endsc = ILogsum(ILogsum(endsc, dmx[cur][M] + TSC(cp9O_DM,M)), /* transition from D_M -> end */
		  imx[cur][M] + TSC(cp9O_IM,M)); /* transition from I_M -> end */
  for(c = 0; c < cm->cp9->el_from_ct[M+1]; c++) /* el_from_ct[k] is >= 0 */
    endsc = ILogsum(endsc, elmx[cur][cm->cp9->el_from_idx[M+1][c]]);
  /* transition penalty to EL incurred when EL was entered */
  /*printf("endsc: %d\n", endsc);*/
  erow[cur] = endsc;
  scA[jp]   = endsc;
  fsc = Scorify(endsc);
  /*printf("jp: %d j: %d fsc: %f sc: %d\n", jp, j, fsc, scA[jp]);*/
  if (fsc > best_sc)      { best_sc     = fsc; best_pos    = j; }
  if (fsc > best_hmm_sc)  { best_hmm_sc = fsc; best_hmm_pos= j; }
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
  /* end of special case position j == i0 */

  /*****************************************************************
   * The main loop: scan the sequence from position i0+1 to j0.
   *****************************************************************/
  /* Recursion. */
  for (j = i0+1; j <= j0; j++)
    {
      int const *isc = cm->cp9->isc[dsq[j]];
      int const *msc = cm->cp9->msc[dsq[j]];
      int el_selfsc  = cm->cp9->el_selfsc;
      int endsc      = -INFTY;
      int sc;

      jp = j-i0+1;     /* jp is relative position in the sequence 1..L */
      cur = (j-i0+1);
      prv = (j-i0);

      if(be_efficient) { cur %= 2; prv %= 2; }
      /* The 1 difference between a Viterbi scanner and the 
       * regular Viterbi. In non-scanner parse must begin in B at
       * position 0 (i0-1), in scanner we can start at any position 
       * in the seq. */
      mmx[cur][0]        = (do_scan == TRUE) ? 0     : -INFTY;
      dmx[cur][0]        = -INFTY;  /*D_0 is non-existent*/
      elmx[cur][0]       = -INFTY;  /*no EL state for node 0 */

      sc = ILogsum(ILogsum(mmx[prv][0] + TSC(cp9O_MI,0),
			   imx[prv][0] + TSC(cp9O_II,0)),
		   dmx[prv][0] + TSC(cp9O_DI,0));
      imx[cur][0] = sc + isc[0];

      /*********************************************************/
      /* special case node: k == 1 */
      k = 1; 
      sc = ILogsum(ILogsum(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
			   imx[prv][k-1] + TSC(cp9O_IM,k-1)),
		   ILogsum(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
			   mmx[prv][k-1] + TSC(cp9O_BM,k  )));
      /* check possibility we came from an EL, if they're valid */
      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	sc = ILogsum(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
      /* transition penalty to EL incurred when EL was entered */
      mmx[cur][k] = sc + msc[k];			   

      /* E state update, put me inside mmx[cur][k] calc */
      endsc = ILogsum(endsc, mmx[cur][k] + TSC(cp9O_ME,k));
      
      /*insert state*/
      sc = ILogsum(ILogsum(mmx[prv][k] + TSC(cp9O_MI,k),
			   imx[prv][k] + TSC(cp9O_II,k)),
		     dmx[prv][k] + TSC(cp9O_DI,k));
      imx[cur][k] = sc + isc[k];
      
      /*delete state*/
      sc = ILogsum(ILogsum(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
			   imx[cur][k-1] + TSC(cp9O_ID,k-1)),
		   dmx[cur][k-1] + TSC(cp9O_DD,k-1));
      dmx[cur][k] = sc;
      
      /*el state*/
      sc = -INFTY;
      if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) { /* not all HMM nodes have an EL state (for ex: 
								  HMM nodes that map to right half of a MATP_MP) */
	  sc = ILogsum(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
		       elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
      }
      elmx[cur][k] = sc;
      /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
	printf("imx [jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
	printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);
	printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[cur][k]);*/
      /* 
       * end of special case node: k == 1 
       ********************************************************/

      if(be_safe) { /* never use ILogsumNI() */
	/********************************************************
	 * Switch statement for the 4 possible locality modes wraps the
	 * main chunk of the dp recursion. be_safe is TRUE, so we
	 * don't use the faster ILogsumNI(sc1,sc2) calls, which assume
	 * sc1 != -INFTY and sc2 != -INFTY. be_safe is usually TRUE
	 * for a reason, meaning we can't trust sc1 != -INFTY and 
	 * sc2 != -INFTY for all cases.
	 *
	 * Each locality mode has a different block b/c some 
	 * transitions are ALWAYS -INFTY for certain locality modes
	 * so we skip them. If said transitions are not -INFTY,
	 * there's something terribly wrong, and we've had already
	 * died in the cp9_GetLocalityMode() call above.
	 *
	 */
	switch(locality_mode) {
	  
	case CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF: 
	  /* M_0->M_K (begin) transitions for K > 1 and M_K->E (end) 
	   * transitions for K < M are -INFTY for k > 1 (so they have 
	   * no contribution and we can skip them).
	   * Also ELs are illegal, M->EL transitions are -INFTY 
	   * and ELMX cells will always be -INFTY 
	   */
	  for (k = 2; k < M; k++) /* note we only go to M-1, we have to handle M 
				   * special b/c it's the only node we can END from */
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     dmx[prv][k-1] + TSC(cp9O_DM,k-1));
	      mmx[cur][k] = sc + msc[k];
	      
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];
	      
	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      elmx[cur][k] = -INFTY;
	    }
	  break; /* end of case CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF: */

	case CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON:
	  /* M_0->M_K (begin) transitions for K > 1 and M_K->E (end) 
	   * transitions for K < M are -INFTY for k > 1 (so they have 
	   * no contribution and we can skip them).
	   */
	  for (k = 2; k < M; k++)
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     dmx[prv][k-1] + TSC(cp9O_DM,k-1));
	      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
		sc = ILogsumNI(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
	      mmx[cur][k] = sc + msc[k];
	    
	    
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];

	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      sc = -INFTY;
	      if(cm->cp9->has_el[k]) { /* not all HMM nodes have an EL state (for ex: 
					  HMM nodes that map to right half of a MATP_MP) */
		sc = ILogsumNI(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
			       elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	      }
	      elmx[cur][k] = sc;
	    }
	  break;

	case CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF:
	  /* ELs are illegal M->EL transitions are -INFTY (so they have no
	   * contribution and we can skip them)
	   * and ELMX cells will always be -INFTY */
	  for (k = 2; k < M; k++)
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     ILogsumNI(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
				       mmx[prv][0]   + TSC(cp9O_BM,k  )));
	      mmx[cur][k] = sc + msc[k];
	    
	      /* E state update */
	      endsc = ILogsumNI(endsc, mmx[cur][k] + TSC(cp9O_ME,k));
	    
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];

	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      elmx[cur][k] = -INFTY;
	    }
	  break; /* end of case CP9_LOCAL_BEGIN_END_ON_AND_EL_ON: */

	case CP9_LOCAL_BEGIN_END_ON_AND_EL_ON:
	  for (k = 2; k < M; k++)
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     ILogsumNI(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
				       mmx[prv][0]   + TSC(cp9O_BM,k  )));
	      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
		sc = ILogsumNI(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
	      mmx[cur][k] = sc + msc[k];
	    
	      /* E state update */
	      endsc = ILogsumNI(endsc, mmx[cur][k] + TSC(cp9O_ME,k));
	    
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];

	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      sc = -INFTY;
	      if(cm->cp9->has_el[k]) { /* not all HMM nodes have an EL state (for ex: 
					  HMM nodes that map to right half of a MATP_MP) */
		sc = ILogsumNI(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
			       elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	      }
	      elmx[cur][k] = sc;
	      /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
		printf("imx [jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
		printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);
		printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[cur][k]);*/
	    }
	  break; /* end of case CP9_LOCAL_BEGIN_END_ON_AND_EL_ON: */
	}
	/* end of switch(locality_mode) */
	/**********************************************************/
      } /* end of if(be_safe) */
	/**********************************************************/
      else { /* don't be safe */
	/********************************************************
	 * Switch statement for the 4 possible locality modes wraps the
	 * main chunk of the dp recursion. This chunk is because it
	 * loops over nodes k = 2..M-1, and we know that certain
	 * transition scores and dp cells for these nodes are NEVER
	 * -INFTY.  We can exploit this knowledge by calling a special
	 * ILogsumNI() function that assumes neither of the scores
	 * passed in are -INFTY. The following is true for all 4
	 * locality modes:
	 *
	 * mmx[prv][k],   imx[prv][k],   dmx[prv][k]   != -INFTY
	 * mmx[cur][k],   mmx[cur][k],   dmx[cur][k]   != -INFTY
	 * mmx[prv][k-1], imx[prv][k-1], dmx[prv][k-1] != -INFTY
	 * mmx[cur][k-1], mmx[cur][k-1], dmx[cur][k-1] != -INFTY
	 *
	 * tsc[MM][k],    tsc[IM][k],    tsc[DM][k]    != -INFTY
	 * tsc[MI][k],    tsc[II][k],    tsc[DI][k]    != -INFTY
	 * tsc[MD][k],    tsc[ID][k],    tsc[DD][k]    != -INFTY
	 * tsc[MM][k-1],  tsc[IM][k-1],  tsc[DM][k-1]  != -INFTY
	 * tsc[MI][k-1],  tsc[II][k-1],  tsc[DI][k-1]  != -INFTY
	 * tsc[MD][k-1],  tsc[ID][k-1],  tsc[DD][k-1]  != -INFTY
	 *
	 * Further, each locality mode (except the last:
	 * CP9_LOCAL_BEGIN_END_ON_AND_EL_ON) has more guarantees as
	 * commented below.
	 *
	 * All of the transition score guarantees are checked 
	 * for at the beginning of this function (at least currently)
	 * to make sure they're valid. The matrix score guarantees
	 * are checked (within ILogsumNI()) if the code is 
	 * configured/compiled in debugging mode (Although the
	 * transition score guarantees imply the matrix score
	 * guarantees are correct). Compiling in debugging mode:
	 * $ sh ./configure --enable-debugging;
	 */
	switch(locality_mode) {

	case CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF: 
	  /* M_0->M_K (begin) transitions for K > 1 and M_K->E (end) 
	   * transitions for K < M are -INFTY for k > 1 (so they have 
	   * no contribution and we can skip them).
	   * Also ELs are illegal, M->EL transitions are -INFTY 
	   * and ELMX cells will always be -INFTY 
	   */
	  for (k = 2; k < M; k++) /* note we only go to M-1, we have to handle M 
				   * special b/c it's the only node we can END from */
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     dmx[prv][k-1] + TSC(cp9O_DM,k-1));
	      mmx[cur][k] = sc + msc[k];
	    
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];

	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      elmx[cur][k] = -INFTY;
	    }
	  break; /* end of case CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF: */

	case CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON:
	  /* M_0->M_K (begin) transitions for K > 1 and M_K->E (end) 
	   * transitions for K < M are -INFTY for k > 1 (so they have 
	   * no contribution and we can skip them).
	   */
	  for (k = 2; k < M; k++)
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     dmx[prv][k-1] + TSC(cp9O_DM,k-1));
	      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
		sc = ILogsumNI(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
	      mmx[cur][k] = sc + msc[k];
	    
	    
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];

	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      sc = -INFTY;
	      if(cm->cp9->has_el[k]) { /* not all HMM nodes have an EL state (for ex: 
					  HMM nodes that map to right half of a MATP_MP) */
		sc = ILogsumNI(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
			       elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	      }
	      elmx[cur][k] = sc;
	    }
	  break;

	case CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF:
	  /* ELs are illegal M->EL transitions are -INFTY (so they have no
	   * contribution and we can skip them)
	   * and ELMX cells will always be -INFTY */
	  for (k = 2; k < M; k++)
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     ILogsumNI(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
				       mmx[prv][0]   + TSC(cp9O_BM,k  )));
	      mmx[cur][k] = sc + msc[k];
	    
	      /* E state update */
	      endsc = ILogsumNI(endsc, mmx[cur][k] + TSC(cp9O_ME,k));
	    
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];

	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      elmx[cur][k] = -INFTY;
	    }
	  break; /* end of case CP9_LOCAL_BEGIN_END_ON_AND_EL_ON: */

	case CP9_LOCAL_BEGIN_END_ON_AND_EL_ON:
	  for (k = 2; k < M; k++)
	    {
	      /*match state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
				       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
			     ILogsumNI(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
				       mmx[prv][0]   + TSC(cp9O_BM,k  )));
	      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
		sc = ILogsumNI(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
	      mmx[cur][k] = sc + msc[k];
	    
	      /* E state update */
	      endsc = ILogsumNI(endsc, mmx[cur][k] + TSC(cp9O_ME,k));
	    
	      /*insert state*/
	      sc = ILogsumNI(ILogsumNI(mmx[prv][k] + TSC(cp9O_MI,k),
				       imx[prv][k] + TSC(cp9O_II,k)),
			     dmx[prv][k] + TSC(cp9O_DI,k));
	      imx[cur][k] = sc + isc[k];

	      /*delete state*/
	      sc = ILogsumNI(ILogsumNI(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
				       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
			     dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	      dmx[cur][k] = sc;

	      /*el state*/
	      sc = -INFTY;
	      if(cm->cp9->has_el[k]) { /* not all HMM nodes have an EL state (for ex: 
					  HMM nodes that map to right half of a MATP_MP) */
		sc = ILogsumNI(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
			       elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	      }
	      elmx[cur][k] = sc;
	      /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
		printf("imx [jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
		printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);
		printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[cur][k]);*/
	    }
	  break; /* end of case CP9_LOCAL_BEGIN_END_ON_AND_EL_ON: */
	}
	/* end of switch(locality_mode) */
	/**********************************************************/
      } /* end of else { } (which came after if(be_safe) { } */
      /**********************************************************/

      /* finish up with special case k = M, we don't have the zero
       * -INFTY guarantees anymore for begins/ends, and since we're no
       * longer in a locality mode specific code block we check for
       * ELs. Since we're outside of the if(be_safe {} else {} block
       * we don't even know if any of our guarantees hold (we wouldn't
       * gain much from them here anyway, k==M is only one node).
       */
      k = M;
      /*match state*/
      sc = ILogsum(ILogsum(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
			   imx[prv][k-1] + TSC(cp9O_IM,k-1)),
		   ILogsum(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
			   mmx[prv][0]   + TSC(cp9O_BM,k  )));
      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	sc = ILogsum(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
      mmx[cur][k] = sc + msc[k];
      
      /* E state update, put me inside mmx[cur][k] calc */
      endsc = ILogsum(endsc, mmx[cur][k] + TSC(cp9O_ME,k));
      
      /*insert state*/
      sc = ILogsum(ILogsum(mmx[prv][k] + TSC(cp9O_MI,k),
			   imx[prv][k] + TSC(cp9O_II,k)),
		   dmx[prv][k] + TSC(cp9O_DI,k));
      imx[cur][k] = sc + isc[k];
      
      /*delete state*/
      sc = ILogsum(ILogsum(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
			   imx[cur][k-1] + TSC(cp9O_ID,k-1)),
		   dmx[cur][k-1] + TSC(cp9O_DD,k-1));
      dmx[cur][k] = sc;
      
      /*el state*/
      sc = -INFTY;
      if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								HMM nodes that map to right half of a MATP_MP) */
	{
	  sc = ILogsum(mmx [cur][k] + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
		       elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	}
      elmx[cur][k] = sc;
      /* end of k == M special case. */
      /******************************************************************/
      endsc = ILogsum(ILogsum(endsc, dmx[cur][M] + TSC(cp9O_DM,M)), 
		      imx[cur][M] + TSC(cp9O_IM,M)); /* transition from D_M -> end and I_M -> end*/
      for(c = 0; c < cm->cp9->el_from_ct[M+1]; c++)  /* el_from_ct[k] is >= 0 */
	endsc = ILogsum(endsc, elmx[cur][cm->cp9->el_from_idx[M+1][c]]);
	/* transition penalty to EL incurred when EL was entered */
      
      /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
	printf("imx [jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
	printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);
	printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[cur][k]);
	printf("endsc: %d\n", endsc);*/

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
		accept = TRUE;
		/* Potentially rescan just the subseq that is the hit we're about to report.
		 * Implemented to deal with fact that --enfseq option was enforcing the subseq
		 * to have hit pass filter b/c this Forward scanning function is 'infinite' length
		 * (Weinberg). Sometimes the subseq we're about to report has a really
		 * crappy score, even though the cumulative Forward score (starting at i0) is good. */
		if(do_scan && cm->search_opts & CM_SEARCH_HMMRESCAN && doing_rescan == FALSE)
		  {
		    /*printf("rechecking hit from %d to %d\n", gback[jp], j);*/
		    temp_sc = CP9Forward(cm, dsq, gback[jp], j, cm->W, cutoff, 
					 NULL,  /* don't care about scores of each pos */
					 NULL,  /* don't care about best scoring position */
					 NULL,  /* don't report hits to results data structure */
					 TRUE,  /* we're scanning */
					 FALSE, /* we're not ultimately aligning */
					 TRUE,  /* set the doing_rescan arg to TRUE, 
						   so we don't potentially infinitely recurse */
					 TRUE,  /* be memory efficient */
					 NULL); /* don't want the DP matrix back */
		    /*printf("new score: %f old score %f\n", temp_sc, savesc[jp]);*/
		    if(temp_sc >= cutoff) 
		      { 
			accept = TRUE; 
			/*printf("rechecked hit from %d to %d\n", gback[jp], j);
			  printf("new score: %f old score %f\n", temp_sc, savesc[jp]);*/
		      }
		    else accept = FALSE;
		    savesc[jp] = temp_sc;
		  }
		if(accept)
		  {
		    if(results != NULL) /* report the hit */
		      {
			/*printf("FWD reporting hit: i: %d j: %d sc: %f\n", gback[jp], j, savesc[jp]);*/
			report_hit(gback[jp], j, 0, savesc[jp], results); 
			/* 0 is for saver, which is irrelevant for HMM hits */
		      }
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
  free(otsc);

  /* determine score to return: (I know, too complex) */
  if(doing_align)
    {
      return_sc = Scorify(scA[(j0-i0+1)]); /* L = j0-i0+1 */
      if(ret_bestpos != NULL) *ret_bestpos = i0;
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

  /* printf("ctr: %d\n", ctr); */
  return return_sc;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}

/***********************************************************************
 * Function: cp9_EXPTLFastBackward()
 * 
 * Purpose:  Runs the Backward dynamic programming algorithm on an
 *           input subsequence (i0-j0). Complements CP9Forward().  
 *           Somewhat flexible based on input options as follows:
 *
 *           if(be_efficient): only allocates 2 rows of the Backward
 *           matrix, else allocates full L+1 matrix.
 *
 *           if(do_scan): allows parses to end at any position j
 *           i0..j0, changing meaning of DP matrix cells as discussed
 *           below.
 *
 *           Reference for algorithm (when do_scan is FALSE): 
 *           Durbin et. al. Biological Sequence Analysis; p. 59.
 *           With 1 IMPORTANT difference, emission scores for
 *           residue at posn j are part of the sum in DP cells
 *           for position j, but in Durbin, emission scores for 
 *           residue at posn j+1 are part of the sum in DP cells
 *           for position j. The Durbin method makes it more 
 *           straightforward to combine Backward and Forward
 *           cells to get posteriors, but it causes precision issues
 *           (overall Backward Score P(seq|X) != overall Forward Score
 *            P(seq|X) solely due to integer log odds scaling precision
 *            issues as investigated in 
 *            ~nawrockie/notebook/7_0410_inf_hmmfb_hbanded_scan/00LOG)
 *           So I've resorted to the algorithm implemented here, the
 *           meaning of the Backward DP cells is given below. This
 *           implementation requires the subtraction of an emission
 *           score when combining corresponding Forward and Backward
 *           DP cells to get posteriors, b/c they've been double
 *           counted. 
 *
 *           I've wasted a lot of time rewrapping my head around this
 *           function when I revisit it, so I'll be verbose about the
 *           the meaning of the Backward (B) matrix DP cells for
 *           matches (M) inserts (I) and deletes (D) here:
 *
 *           For relative subsequence positions jp = 0..L:
 *           For HMM nodes 1..M: 
 *           B->M[jp][k] : sum of all parses emitting seq
 *                         from jp+1..j0 that visit node k's match 
 *                         state, which emitted posn jp
 *           B->I[jp][k] : sum of all parses emitting seq from 
 *                         jp+1..j0 that visit node k's insert
 *                         state, which emitted posn jp 
 *           B->D[jp][k] : sum of all parses emitting seq from 
 *                         jp+1..j0 that visit node k's delete
 *                         delete state, last emitted (rightmost)
 *                         posn was jp+1
 *           B->EL[jp][k]: sum of all parses emitting seq from 
 *                         jp+1..j0 that visit node k's EL
 *                         state, which MAY OR MAY NOT have
 *                         emitted any posns >= jp+1, last
 *                         emitted (rightmost) posn was jp+1.
 *                         Some nodes k do not have an EL state 
 *                         (if cp9->has_el[k] == FALSE)
 *                         NOTE: EL can act as non-emitter if 0
 *                         self loops taken or emitter if >= 1
 *                         self loops taken, this is why we
 *                         treat jp as having NOT YET BEEN EMITTED
 *                         which is different than the M and I 
 *                         convention.
 *                          
 *           For *special* HMM node 0:
 *           B->M[jp][0] : M_0 is the Begin state, which does not 
 *                         emit, so this is the sum of all parses 
 *                         emitting seq from jp+1..j0 that start
 *                         in the begin state, the last emitted 
 *                         (rightmost) posn was jp+1.
 *
 *           Note: if jp=0, only D and M_0 states can have 
 *                 non-IMPOSSIBLE values. 
 *
 *           if(do_scan) the 'jp+1..j0' in the above definitions is
 *           changed to jp+1..jE such that jp+1 <= jE <= j0. Meaning
 *           any residue can be the final residue emitted in the
 *           parse. This means B->M[jp][0] is the sum of all parses
 *           emitting a subseq ending anywhere from jp+1..j0 and 
 *           starting at jp+1.
 *
 *           The *will emit* in the above definitions refers to 
 *           the fact that emission scores from a state x are not 
 *           counted in the matrix score for state x, and are only
 *           added when we calculate a matrix score for state y,
 *           after transitioning (backwards) from y to x. See
 *           code. This is done to facilitate combining Forward
 *           and Backward cells to get posterior probabilities
 *           in CP9Posterior().
 *
 * Args:     
 *           cm        - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq       - sequence in digitized form
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - the maximum size of a hit (often cm->W)
 *           cutoff    - minimum score to report
 *           ret_sc    - RETURN: int log odds Backward score for each end point [0..(j0-i0+1)]
 *           ret_bestpos- RETURN: start position that gives maximum score max argmax_i sc[i]
 *           results   - search_results_t to add to; if NULL, don't keep results
 *           do_scan   - TRUE if we're scanning, HMM can start to emit anywhere i0..j0,
 *                       FALSE if we're not, HMM must start emitting at i0, end emitting at j0
 *           doing_align  - TRUE if reason we've called this function is to help get posteriors
 *                          for CM alignment, in which case we can skip the traceback. 
 *           doing_rescan - TRUE if we've called this function recursively, and we don't
 *                          want to again
 *           be_efficient- TRUE to keep only 2 rows of DP matrix in memory, FALSE keep whole thing
 *           ret_mx    - RETURN: CP9 Forward DP matrix, NULL if not wanted
 *
 * Returns:  if(!do_scan) log P(S|M)/P(S|R), as a bit score, this is B->M[0][0];
 *           else         max log P(S|M)/P(S|R), for argmax subseq S of input seq i0..j0,
 *                            this is max_jp B->M[jp][0]
 */
float
cp9_EXPTLFastBackward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
	    int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
	    int doing_rescan, int be_efficient, CP9_dpmatrix_t **ret_mx)
{
  int          status;
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /*     j-W: position in the subsequence                         */
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
  int         *sc;          /* prob (seq from j0..jp | HMM) [0..jp..cm->cp9->M]             */
  float        fsc;         /* float log odds score                                         */
  float        curr_sc;     /* temporary score used for filling in gamma                    */
  float       *gamma;       /* SHMM DP matrix for optimum nonoverlap resolution [0..j0-i0+1]*/
  int         *gback;       /* traceback pointers for SHMM                      [0..j0-i0+1]*/ 
  float       *savesc;      /* saves score of hit added to best parse at j      [0..j0-i0+1]*/ 
  float        best_hmm_sc; /* Best overall score from semi-HMM to return if do_scan        */
  float        best_hmm_pos;/* residue giving best_hmm_sc                                   */
  float        best_sc;     /* Best score overall, returned if 0 hits found by HMM & do_scan*/
  float        best_pos;    /* residue giving best_sc                                       */
  float        return_sc;   /* score to return, if (!do_scan) return overall Backward sc,   *
			     * else return best_hmm_sc if # HMM hits>0, else return best_sc */
  int          accept;      /* Flag used if cm->search_opts & CM_SEARCH_HMMRESCAN           */
  float        temp_sc;     /* temporary score                                              */
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */
  int          c;           /* counter for EL states */

  printf("in cp9_EXPTLFastBackward() i0: %d j0: %d do_scan: %d \n", i0, j0, do_scan);  
  /* Contract checks */
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in cp9_EXPTLFastBackward, cm->cp9 is NULL.\n");
  if(doing_rescan && !do_scan) 
    esl_fatal("ERROR in cp9_EXPTLFastBackward, doing_rescan but not do_scan");
  if(be_efficient && (ret_mx != NULL))
    esl_fatal("ERROR in cp9_EXPTLFastBackward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    esl_fatal("ERROR in cp9_EXPTLFastBackward, passing in results data structure, but not in scanning mode.a\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    esl_fatal("ERROR in cp9_EXPTLFastBackward, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
  if(dsq == NULL)
    esl_fatal("ERROR in cp9_EXPTLFastBackward, dsq is NULL.");
    
  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_hmm_sc = IMPOSSIBLE;
  best_hmm_pos= -1;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. We use it first for
   * the Forward scan, and reuse it for the Backward scan.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+2));
  gamma[0]     = 0.; /* this is impossible, no hit can be rooted at ip=0, 
		      * we require hits to have at least 1 emission,
		      * plus, no positive scoring hit can have 0 emissions as 
		      * all transition scores are negative */
  gamma[L]   = 0.;
  ESL_ALLOC(gback, sizeof(int)   * (L+2));
  gback[L]   = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+2));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 2;
  else             nrows = L+1; 
  mx = AllocCPlan9Matrix(nrows, cm->cp9->M, &mmx, &imx, &dmx, &elmx, &erow);

  /* sc will hold P(seq from i..j0 | Model) for each i in int log odds form */
  ESL_ALLOC(sc, sizeof(int) * (j0-i0+3));

  /* Initialization of the L (i = j0, cur = (j0-i) = (j0-j0) %2 = 0) row. */
  if(be_efficient) cur = 0;
  else cur = j0-i0+1; /* L */
  ip = j0-i0+1;  /*L */
  i  = j0;

  /*******************************************************************
   * 0 Handle EL, looking at EL_k->E for all valid k.
   * we're going backwards so we have to work out of order, we could get 
   * around this by storing the nodes each EL goes TO in an el_to_ct[] vec. */
  /* init to -INFTY */
  for (k = 1; k <= cm->cp9->M; k++)
    elmx[cur][k] = -INFTY;
  if(cm->cp9->flags & CPLAN9_EL)
    {
      for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M+1]; c++) /* el_from_ct[cm->cp9->M+1] holds # ELs that can go to END */
	elmx[cur][cm->cp9->el_from_idx[cm->cp9->M+1][c]] = 0.; /* EL<-E, penalty incurred when we enter EL (i.e. leave going backwards) */
    }
  /*******************************************************************/

  /* elmx[cur][cm->cp9->M] is either 0 (if EL_M exists (it would nec be in el_from_idx[cm->cp9->M+1] array if it does, so
   * it would be filled with 0 in above loop), or -INFTY if it doesn't exist. We don't add possibility of EL_M -> EL_M
   * self loop b/c it's impossible to do that without emitting, and we've already seen our last res emitted. 
   * either way we don't have to modify it */

  mmx[cur][cm->cp9->M]  = 0. + 
    ILogsum(elmx[cur][cm->cp9->M] + cm->cp9->tsc[CTMEL][cm->cp9->M],/* M_M<-EL_M<-E, with 0 self loops in EL_M */
	    cm->cp9->esc[cm->cp9->M]);                             /* M_M<-E ... everything ends in E (the 0; 2^0=1.0) */
  mmx[cur][cm->cp9->M] += cm->cp9->msc[dsq[i]][cm->cp9->M];  /* ... + emitted match symbol */
  imx[cur][cm->cp9->M]  = 0. + cm->cp9->tsc[CTIM][cm->cp9->M];     /* I_M<-E ... everything ends in E (the 0; 2^0=1.0) */
  imx[cur][cm->cp9->M] += cm->cp9->isc[dsq[i]][cm->cp9->M];  /* ... + emitted insert symbol */
  dmx[cur][cm->cp9->M]  = cm->cp9->tsc[CTDM][cm->cp9->M];          /* D_M<-E */

  /*******************************************************************
   * No need to look at EL_k->M_M b/c elmx[cur] with cur == L means last emitted residue was L+1 
   * and this is impossible if we've come from M_M (only would be valid if we were coming from
   * E which is handled above with the EL_k->E code). 
   *******************************************************************/

  for (k = cm->cp9->M-1; k >= 1; k--)
    {
      mmx[cur][k]  = 0 + cm->cp9->esc[k];  /*M_k<- E */
      mmx[cur][k]  = ILogsum(mmx[cur][k], dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]);
      if(cm->cp9->flags & CPLAN9_EL)
	mmx[cur][k]  = ILogsum(mmx[cur][k], elmx[cur][k] + cm->cp9->tsc[CTMEL][k]);
      mmx[cur][k] += cm->cp9->msc[dsq[i]][k];

      /*******************************************************************
       * No need to look at EL_k->M_M b/c elmx[cur] with cur == L means last emitted residue was L+1 
       * and this is impossible if we've come from M_M (only would be valid if we were coming from
       * E which is handled above with the EL_k->E code). 
       *******************************************************************/
      imx[cur][k]  = dmx[cur][k+1] + cm->cp9->tsc[CTID][k];
      imx[cur][k] += cm->cp9->isc[dsq[i]][k];

      dmx[cur][k]  = dmx[cur][k+1] + cm->cp9->tsc[CTDD][k];
      /* elmx[cur][k] was set above, out of order */
    }
  
  /* remember M_0 is special, the B state, a non-emitter */
  mmx[cur][0]  = dmx[cur][1] + cm->cp9->tsc[CTMD][0]; /* M_0(B)->D_1, no seq emitted, all deletes */
  /* above line is diff from CPBackwardOLD() which has mmx[cur][0] = -INFTY; */
  imx[cur][0]  = dmx[cur][1] + cm->cp9->tsc[CTID][0];
  imx[cur][0] += cm->cp9->isc[dsq[i]][0];

  dmx[cur][0]   = -INFTY; /*D_0 doesn't exist*/
  elmx[cur][0]  = -INFTY; /*EL_0 doesn't exist*/

  sc[ip] = mmx[cur][0]; /* all parses must start in M_0, the B state */
  fsc = Scorify(sc[ip]);
  /*printf("ip: %d i: %d fsc: %f i: %d\n", ip, i, fsc, sc[ip]);*/
  /* we can't have a hit starting here, b/c it would correspond to all deletes,
   * no seq emitted, so we skip the check for if(fsc > best_sc) */

  /*printf("sc[ip:%d]: %d\n", ip, sc[ip]);*/

  /*printf("mmx[ip:%d][%d]: %d cur: %d\n", L+1, 0, mmx[cur][0], cur);
    printf("imx[ip:%d][%d]: %d cur: %d\n", L+1, 0, imx[cur][0], cur);
    printf("dmx[ip:%d][%d]: %d cur: %d\n", L+1, 0, dmx[cur][0], cur);*/
  
  /*****************************************************************
   * The main loop: scan the sequence from position j0-1 to i0.
   *****************************************************************/
  /* Recursion */
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  for (i = j0-1; i >= i0; i--) 
    {
      ip = i-i0+1;		/* ip is relative index in dsq (0 to L-1) */
      if(be_efficient) { cur = (j0-i)  %2; prv = (j0-i+1)%2; }	  
      else { cur = ip; prv = ip+1; }

      /* init EL mx to -INFTY */
      for (k = 1; k <= cm->cp9->M; k++) elmx[cur][k] = -INFTY;
      
      /* elmx[cur][k] is possibly of coming from self (EL_k), we 
       * can't have come from END b/c we haven't emitted the last res of the seq yet.
       */
      if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[cm->cp9->M]))
	    elmx[cur][cm->cp9->M] = elmx[cur][cm->cp9->M] + cm->cp9->el_selfsc;

      mmx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTMI][cm->cp9->M];
      mmx[cur][cm->cp9->M] += cm->cp9->msc[dsq[i]][cm->cp9->M];
      
      if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[cm->cp9->M]))
	mmx[cur][cm->cp9->M] = ILogsum(mmx[cur][cm->cp9->M], 
				       elmx[cur][cm->cp9->M] + cm->cp9->tsc[CTMEL][cm->cp9->M]);
      
      imx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTII][cm->cp9->M];
      imx[cur][cm->cp9->M] += cm->cp9->isc[dsq[i]][cm->cp9->M];
      
      dmx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTDI][cm->cp9->M]; 
      
      /*******************************************************************
       * 1b Handle EL, looking at EL_k->M_M for all valid k.
       * EL_k->M_M transition, which has no transition penalty */
      if(cm->cp9->flags & CPLAN9_EL)
	{
	  for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M]; c++) /* el_from_ct[cm->cp9->M] holds # ELs that can go to M_M */
	    elmx[cur][cm->cp9->el_from_idx[cm->cp9->M][c]] = ILogsum(elmx[cur][cm->cp9->el_from_idx[cm->cp9->M][c]], mmx[prv][cm->cp9->M]);
	}
      
      /* A main difference between a Backward scanner and 
       * regular Backward: a scanner can end at the END 
       * state at any position, regular can only end at
       * the final position j0. */
      if(do_scan)
	{	
	  /*******************************************************************
	   * 2 Handle EL, looking at EL_k->E for all valid k.
	   * EL_k->M_M transition, which has no transition penalty */
	  if(cm->cp9->flags & CPLAN9_EL)
	    {
	      for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M+1]; c++) /* el_from_ct[cm->cp9->M] holds # ELs that can go to END */
		elmx[cur][cm->cp9->el_from_idx[cm->cp9->M+1][c]] = 0.; /* EL<-E, penalty incurred when we enter EL (i.e. leave going backwards) */
	    }
	  /*******************************************************************/
	  /* elmx[cur][cm->cp9->M] is either 0 (if EL_M exists (it would nec be in el_from_idx[cm->cp9->M+1] array if it does, so
	   * it would be filled with 0 in above loop), or -INFTY if it doesn't exist. We don't add possibility of EL_M -> EL_M
	   * self loop b/c it's impossible to do that without emitting, and we've already seen our last res emitted,
	   * either way we don't have to modify it */
	  
	  mmx[cur][cm->cp9->M]  =  
	    ILogsum(mmx[cur][cm->cp9->M], 
		    ILogsum(elmx[cur][cm->cp9->M] + cm->cp9->tsc[CTMEL][cm->cp9->M],/* M_M<-EL_M<-E, with 0 selfs in EL_M */
			    cm->cp9->esc[cm->cp9->M]));                             /* M_M<-E ... */
	  ///mmx[cur][cm->cp9->M] += cm->cp9->msc[dsq[i]][cm->cp9->M]; /* ... + emitted match symbol */
	  /* DO NOT add contribution of emitting i from M, it's been added above */
	  
	  imx[cur][cm->cp9->M]  =
	    ILogsum(imx[cur][cm->cp9->M],
		    (cm->cp9->tsc[CTIM][cm->cp9->M] +            /* I_M<-E + (only in scanner)     */
		     0));                                        /* all parses end in E, 2^0 = 1.0;*/
	  ///imx[cur][cm->cp9->M] += cm->cp9->isc[dsq[i]][cm->cp9->M]; /* ... + emitted insert symbol */  
	  /* DO NOT add contribution of emitting i from M, it's been added above */
	  dmx[cur][cm->cp9->M] =  
	    ILogsum(dmx[cur][cm->cp9->M],
		    (cm->cp9->tsc[CTDM][cm->cp9->M] +            /* D_M<-E + (only in scanner)     */
		     0));                                        /* all parses end in E, 2^0 = 1.0;*/
	}
      /*printf("mmx[ip:%d][%d]: %d cur: %d\n", ip, cm->cp9->M, mmx[cur][cm->cp9->M], cur);
	printf("imx[ip:%d][%d]: %d cur: %d\n", ip, cm->cp9->M, imx[cur][cm->cp9->M], cur);
	printf("dmx[ip:%d][%d]: %d cur: %d\n", ip, cm->cp9->M, dmx[cur][cm->cp9->M], cur);*/
      
      for (k = cm->cp9->M-1; k >= 1; k--)
	{
	  /*******************************************************************
	   * 3 Handle EL, looking at EL_k->M_k for all valid k and EL_k->EL_k
	   * we're going backwards so we have to work out of order
	   * we could get around this by storing the nodes each EL goes TO
	   * in an el_to_ct[] vector. */
	  if(cm->cp9->flags & CPLAN9_EL) {
	    for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] holds # ELs that can go to M_k */
	      elmx[cur][cm->cp9->el_from_idx[k][c]] = ILogsum(elmx[cur][cm->cp9->el_from_idx[k][c]], mmx[prv][k]);
	    /* EL<-M, penalty incurred when we enter EL (i.e. leave going backwards) */
	  }
	  /*******************************************************************/
	  
	  /* Finish off elmx[cur][k] with possibility of coming from self (EL_k), 
	   * elmx[cur][k] will have been filled by block above for ks > current k,
	   * no M_k -> EL_k' with k' > k */
	  if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[k]))
	    elmx[cur][k] = ILogsum(elmx[cur][k], elmx[prv][k] + cm->cp9->el_selfsc);

	  mmx[cur][k]  = ILogsum(ILogsum((mmx[prv][k+1] + cm->cp9->tsc[CTMM][k]),  
					 (imx[prv][k]   + cm->cp9->tsc[CTMI][k])),
				 (dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]));
	  if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[k]))
	    mmx[cur][k] = ILogsum(mmx[cur][k], elmx[cur][k] + cm->cp9->tsc[CTMEL][k]); /* penalty for entering EL */
	  mmx[cur][k] += cm->cp9->msc[dsq[i]][k];
	  
	  imx[cur][k]  = ILogsum(ILogsum((mmx[prv][k+1] + cm->cp9->tsc[CTIM][k]),
					 (imx[prv][k]   + cm->cp9->tsc[CTII][k])),
				 (dmx[cur][k+1] + cm->cp9->tsc[CTID][k]));
	  imx[cur][k] += cm->cp9->isc[dsq[i]][k];
	  
	  if(do_scan) { /* add possibility of ending at this position from this state */
	    mmx[cur][k] = 
	      ILogsum(mmx[cur][k], 
		      (cm->cp9->esc[k] +                  /* M_k<-E + (only in scanner)     */ 
		       0));                               /* all parses end in E, 2^0 = 1.0;*/
	    /* DO NOT add contribution of emitting i from M, it's been added above */
	    /* No EL contribution here b/c we'd be looking for M_k<-EL_k<-E, but EL_k<-E is impossible 
	     * for k != cm->cp9->M; */
	  }	      
	  dmx[cur][k]  = ILogsum(ILogsum((mmx[prv][k+1] + cm->cp9->tsc[CTDM][k]),
					 (imx[prv][k]   + cm->cp9->tsc[CTDI][k])),
				 (dmx[cur][k+1] + cm->cp9->tsc[CTDD][k]));
	}
      /* Case when k == 0 */
      /* imx[cur][0] is filled same as imx[cur][1..k] in the loop above */
      imx[cur][0] = ILogsum(ILogsum((mmx[prv][1] + cm->cp9->tsc[CTIM][0]),
				    (imx[prv][0] + cm->cp9->tsc[CTII][0])),
			    (dmx[cur][1] + cm->cp9->tsc[CTID][k]));
      imx[cur][0] += cm->cp9->isc[dsq[i]][k];
      dmx[cur][0]   = -INFTY; /* D_0 does not exist */
      elmx[cur][0]  = -INFTY; /* EL_0 does not exist */

      /*M_0 is the B state, it doesn't emit, and can be reached from any match via a begin transition */
      mmx[cur][0] = -INFTY;
      for (k = cm->cp9->M; k >= 1; k--) 
	mmx[cur][0] = ILogsum(mmx[cur][0], (mmx[prv][k] + cm->cp9->bsc[k]));
      mmx[cur][0] = ILogsum(mmx[cur][0], (imx[prv][0] + cm->cp9->tsc[CTMI][0]));
      mmx[cur][0] = ILogsum(mmx[cur][0], (dmx[cur][1] + cm->cp9->tsc[CTMD][0]));     /* B->D_1 */
      /* No EL contribution here, can't go B->EL_* */
      
      /* determine isc, the int score of all possible parses starting at the current
       * position (i) of the target sequence. */
      sc[ip] = mmx[cur][0]; /* all parses must start in M_0, the B state */
      fsc = Scorify(sc[ip]);
      /*printf("ip: %d i: %d fsc: %f i: %d\n", ip, i, fsc, sc[ip]);*/
      if(fsc > best_sc) { best_sc = fsc; best_pos= i+1; } /* *off-by-one* (see *off-by-one* below) */
      if(!(cm->search_opts & CM_SEARCH_HMMGREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   * *off-by-one*:
	   * There's an off-by-one issue here: all Backward hits are rooted in 
	   * M_O, the B (begin) state, which is a non-emitter.
	   * let i = ip+i0-1 => ip = i-i0+1;
	   * so sc[ip] = backward->mmx[ip][0] = summed log prob of all parses that end at j0, 
	   * and start at position i+1 of the sequence (because i+1 is the last residue
	   * whose emission has been accounted for). As a result, gamma indexing is off-by-one
	   * with respect to sequence position, hence the i+1 or i-1 in the following
	   * code blocks, each marked by "*off-by-one*" comment below. 
	   * for example: let i0 = 2 gamma[ip=4], normally this means ip=4 corresponds to i=5 
	   *              but due to this off-by-one sc[ip=4] corresponds to hits that start at i=6
	   */
	  gamma[ip]  = gamma[ip+1] + 0; /* extend without adding a new hit */
	  /*printf("ip: %d | gamma[ip]: %f | gamma[ip+1]: %f\n", ip, gamma[ip], gamma[ip+1]);*/
	  gback[ip]  = -1;
	  savesc[ip] = IMPOSSIBLE;
	  j = (((i+1)+W-1) < j0) ? ((i+1)+W-1) : j0; /* *off-by-one* */
	  jp = j-i0+1;
	  curr_sc = gamma[jp+1-1] + fsc + cm->cp9_sc_boost; /* *off-by-one* */
	  /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
	   * value is 0.0 if technique not used. */
	  if (curr_sc > gamma[ip]) { 
	    gamma[ip]  = curr_sc;
	    /*printf("\ti: %d | gamma[i]: %f\n", i+1, gamma[ip]);*/
	    gback[ip]  = j;
	    savesc[ip] = fsc;
	  }
	  /*printf("i: %d ip: %d gamma[ip]: %f\n", i, ip, gamma[ip]);*/
	}
      else {
	/* Resolving overlaps greedily (RSEARCH style),  
	 * Return best hit for each j, IFF it's above threshold */
	if (fsc >= cutoff) {
	  if(results != NULL) {
	    j = (((i+1)+W-1) < j0) ? ((i+1)+W-1) : j0; /* *off-by-one* */
	    /*printf("BWD greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i+1, j, fsc);*/ /* *off-by-one* */
	    report_hit (i+1, j, 0, fsc, results); 
	    /* 0 is for saver, which is irrelevant for HMM hits */
	  }
	}
      }
      if (fsc > best_hmm_sc) {
	best_hmm_sc = fsc;
	best_hmm_pos= i;
      }
    }
  esl_stopwatch_Stop(w);
  printf("\n\n");
  esl_stopwatch_Display(stdout, w, "backward bulk CPU time: ");
  printf("\n\n");
  /*******************************************************************/
  /* Special case: ip == 0, i = i0-1; */
  ip = i-i0+1;		/* ip is relative index in dsq (0 to L-1) */
  if(be_efficient)
    {
      cur = (j0-i)  %2;
      prv = (j0-i+1)%2;
    }	  
  else
    {
      cur = ip;
      prv = ip+1;
    }
  /* init EL mx to -INFTY */
  for (k = 1; k <= cm->cp9->M; k++)
    elmx[cur][k] = -INFTY;

  mmx[cur][cm->cp9->M] = -INFTY;  /* need seq to get here */
  imx[cur][cm->cp9->M] = -INFTY;  /* need seq to get here */
  elmx[cur][cm->cp9->M]= -INFTY;  /* first emitted res can't be from an EL, need to see >= 1 matches */
  dmx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTDI][cm->cp9->M]; 
  /* A main difference between a Backward scanner and 
   * regular Backward: a scanner can end at the END 
   * state at any position, regular can only end at
   * the final position j0. */
  if(do_scan)
    {	
      dmx[cur][cm->cp9->M] =  
	ILogsum(dmx[cur][cm->cp9->M],
		(cm->cp9->tsc[CTDM][cm->cp9->M] +            /* D_M<-E + (only in scanner)     */
		 0));                                        /* all parses end in E, 2^0 = 1.0;*/
    }
  for (k = cm->cp9->M-1; k >= 1; k--)
    {
      mmx[cur][k] = -INFTY; /* need seq to get here, unless we come from E in a scanner (below) */
      imx[cur][k] = -INFTY; /* need seq to get here */
      elmx[cur][k]= -INFTY;  /* first emitted res can't be from an EL, need to see >= 1 matches */
      dmx[cur][k]  = ILogsum(ILogsum((mmx[prv][k+1] + cm->cp9->tsc[CTDM][k]),
				     (imx[prv][k]   + cm->cp9->tsc[CTDI][k])),
			     (dmx[cur][k+1] + cm->cp9->tsc[CTDD][k]));
    }

  /* Case when k == 0 */
  /* imx[cur][0] is filled same as imx[cur][1..k] in the loop above */
  imx[cur][0] = -INFTY; /* need seq to get here */
  dmx[cur][0]   = -INFTY; /* D_0 does not exist */
  elmx[cur][0]  = -INFTY; /* EL_0 does not exist */

  /*M_0 is the B state, it doesn't emit, and can be reached from any match via a begin transition */
  mmx[cur][0] = -INFTY;
  for (k = cm->cp9->M; k >= 1; k--) 
    mmx[cur][0] = ILogsum(mmx[cur][0], (mmx[prv][k] + cm->cp9->bsc[k]));
  mmx[cur][0] = ILogsum(mmx[cur][0], (imx[prv][0] + cm->cp9->tsc[CTMI][0]));
  mmx[cur][0] = ILogsum(mmx[cur][0], (dmx[cur][1] + cm->cp9->tsc[CTMD][0]));     /* B->D_1 */
  /* No EL contribution here, can't go B->EL_* */
      
  /* determine isc, the int score of all possible parses starting at the current
   * position (i) of the target sequence. */
  sc[ip] = mmx[cur][0]; /* all parses must start in M_0, the B state */
  fsc = Scorify(sc[ip]);
  /*printf("ip: %d i: %d fsc: %f i: %d\n", ip, i, fsc, sc[ip]);*/
  if(fsc > best_sc) { best_sc = fsc; best_pos= i+1; } /* *off-by-one* (see *off-by-one* below) */
  if(!(cm->search_opts & CM_SEARCH_HMMGREEDY)) { /* resolve overlaps optimally */
    gamma[ip]  = gamma[ip+1] + 0; /* extend without adding a new hit */
    /*printf("ip: %d | gamma[ip]: %f | gamma[ip+1]: %f\n", ip, gamma[ip], gamma[ip+1]);*/
    gback[ip]  = -1;
    savesc[ip] = IMPOSSIBLE;
    j = (((i+1)+W-1) < j0) ? ((i+1)+W-1) : j0; /* *off-by-one* */
    jp = j-i0+1;
    curr_sc = gamma[jp+1-1] + fsc + cm->cp9_sc_boost; /* *off-by-one* */
    /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
     * value is 0.0 if technique not used. */
    if (curr_sc > gamma[ip]) {
      gamma[ip]  = curr_sc;
      /*printf("\ti: %d | gamma[i]: %f\n", i+1, gamma[ip]);*/
      gback[ip]  = j;
      savesc[ip] = fsc;
    }
  }
  else {
    /* Resolving overlaps greedily (RSEARCH style),  
     * Return best hit for each j, IFF it's above threshold */
    if (fsc >= cutoff) {
      if(results != NULL) {
	j = (((i+1)+W-1) < j0) ? ((i+1)+W-1) : j0; /* *off-by-one* */
	/*printf("BWD greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i+1, j, fsc);*/ /* *off-by-one* */
	report_hit (i+1, j, 0, fsc, results); 
	/* 0 is for saver, which is irrelevant for HMM hits */
      }
    }
  }
  if (fsc > best_hmm_sc) {
    best_hmm_sc = fsc;
    best_hmm_pos= i;
  }
  /* end of special case, ip == 0 */
  /**********************************************************************************/
  /* End of Backward recursion */
  
  if((!(cm->search_opts & CM_SEARCH_HMMGREEDY)) && /* resolve overlaps optimally */
     (!doing_align || do_scan)) /* else we can save time by skipping traceback */
    {
      /*****************************************************************
       * Traceback stage for Backward.
       * Recover all hits: an (i,j,sc) triple for each one.
       * Remember the off-by-one b/t seq index and gamma (see *off-by-one* above)
       *****************************************************************/ 
      i     = i0; 
      while (i <= j0) 
	{
	  ip    = (i-1)-i0+1; /* *off-by-one*, i-1, ip corresponds to i+1 
			       * (yes: i *-* 1 =>   ip corresponds to i *+* 1) */
	  if (gback[ip] == -1) /* no hit */
	    i++; 
	  else                /* a hit, a palpable hit */
	    {
	      if(savesc[ip] > best_hmm_sc) 
		{
		  best_hmm_sc = savesc[ip];
		  best_hmm_pos= i;
		}
	      if(savesc[ip] >= cutoff)
		{
		  accept = TRUE;
		  /* Potentially rescan just the subseq that is the hit we're about to report.
		   * Implemented to deal with fact that --enfseq option was enforcing the subseq
		   * to have hit pass filter b/c this Forward scanning function is 'infinite' length
		   * (Weinberg). Sometimes the subseq we're about to report has a really
		   * crappy score, even though the cumulative Backward score (ending at j0) is good. */
		  if(do_scan && cm->search_opts & CM_SEARCH_HMMRESCAN && doing_rescan == FALSE)
		    {
		      /*printf("rechecking hit from %d to %d\n", i, gback[ip]);*/
		      temp_sc = cp9_EXPTLFastBackward(cm, dsq, i, gback[ip], cm->W, cutoff,  /* *off-by-one* i+1 */
					    NULL,  /* don't care about scores of each pos */
					    NULL,  /* don't care about best scoring position */
					    NULL,  /* don't report hits to results data structure */
					    TRUE,  /* we're scanning */
					    FALSE, /* we're not ultimately aligning */
					    TRUE,  /* set the doing_rescan arg to TRUE, 
						      so we don't potentially infinitely recurse */
					    TRUE,  /* be memory efficient */
					    NULL); /* don't want the DP matrix back */
		      /*printf("new score: %f old score %f\n", temp_sc, savesc[ip]);*/
		      if(temp_sc >= cutoff) 
			{ 
			  accept = TRUE; 
			  /*printf("rechecked hit from %d to %d\n", i, gback[ip]);
			    printf("new score: %f old score %f\n", temp_sc, savesc[ip]);*/
			}
		      else accept = FALSE;
		      savesc[ip] = temp_sc;
		    }
		  if(accept)
		    {
		      if(results != NULL) /* report the hit */
			{
			  /*printf("BWD optimal reporting hit i: %d j: %d sc: %f\n", i, gback[ip], savesc[ip]);*/
			  report_hit(i, gback[ip], 0, savesc[ip], results); 
			  /* 0 is for saver, which is irrelevant for HMM hits */
			}
		    }
		}
	      i = gback[ip]+1;
	    }
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  /*printf("returning from cp9_EXPTLFastBackward()\n");*/
  if(ret_sc != NULL) *ret_sc = sc;
  else free(sc);
  /* determine score to return: (I know, too complex) */
  if(doing_align)
    {
      return_sc = Scorify(mmx[cur][0]);
      if(ret_bestpos != NULL) *ret_bestpos = i0;
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
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);
  /*printf("Backward return_sc: %f\n", return_sc);*/
  return return_sc;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}

/***********************************************************************
 * Function: cp9_WorstForward()
 * 
 * Purpose:  Finds the minimum length L of any possible sequence for
 *           which the Forward score <= <thresh>. Used to determine
 *           the maximum length of sequences we can safely use the
 *           optimized Forward implementation that takes advantage 
 *           of guarantees that dp mx cells are > -INFTY 
 *           (in cp9_FastForward()).
 *
 * Args
 *           cm        - the covariance model, includes cm->cp9: a CP9 HMM
 *           thresh    - score threshold, we seek seq length that can have a sc <= thresh
 *           doing_scan- TRUE if we're scanning, HMM can start to emit at any position
 *                       FALSE if we're not, HMM must start emitting at res 1
 *           doing_align  - TRUE if reason we've called this function is to help get posteriors
 *                          for CM alignment, in which case we can skip the traceback. 
 *
 * Returns:  length L at which worst scoring seq of length L has a dp cell with 1 <= k <= M
 *           and i0 <= j <= j0 with sc <= thresh. -1 if no such length exists (this 
 *           is almost always the case if do_scan == TRUE 
 */
float
cp9_WorstForward(CM_t *cm, int thresh, int doing_scan, int doing_align)
{
  int          status;
  int          j;           /*     actual   position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  CP9_dpmatrix_t *mx;       /* the CP9 DP matrix                                            */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cm->cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */
  int          c;           /* counter for EL states */
  int          M;
  int          minsc = 0;
  int          a;
  int          keep_going;
  int          retL;

  /*debug_print_cp9_params(stdout, cm->cp9, TRUE);*/
  /*printf("in cp9_FastForward() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    cm_Fail("in cp9_FastForward, cm->cp9 is NULL.\n");
    
  M = cm->cp9->M;

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
    otsc_k[cp9O_MEL]= cm->cp9->tsc[CTMEL][k];
    otsc_k[cp9O_ME] = cm->cp9->esc[k];
  }
  int const *tsc = otsc;

  /* get minimal scoring emission residue for each insert, match state */
  int *min_isc;
  int *min_msc;
  int *isc_vec;
  int *msc_vec;
  ESL_ALLOC(min_isc, sizeof(int) * (M+1));
  ESL_ALLOC(min_msc, sizeof(int *) * (M+1));
  ESL_ALLOC(isc_vec, sizeof(int) * cm->abc->K);
  ESL_ALLOC(msc_vec, sizeof(int) * cm->abc->K);

  min_msc[0] = -1; /* invalid, M_0 (BEGIN state) doesn't emit */
  for (a = 0; a < cm->abc->K; a++)
    isc_vec[a] = cm->cp9->isc[a][0];
  min_isc[0] = esl_vec_IArgMin(isc_vec, cm->abc->K);

  for(k = 1; k <= M; k++) {
    for (a = 0; a < cm->abc->K; a++) {
      isc_vec[a] = cm->cp9->isc[a][k];
      msc_vec[a] = cm->cp9->msc[a][k];
    }
    min_isc[k] = esl_vec_IMin(isc_vec, cm->abc->K);
    min_msc[k] = esl_vec_IMin(msc_vec, cm->abc->K);
  }
  free(msc_vec);
  free(isc_vec);

  /* Allocate DP matrix, 2 rows */
  nrows = 2;
  mx = AllocCPlan9Matrix(nrows, M, &mmx, &imx, &dmx, &elmx, &erow); 

  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  elmx[0][0]= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY;   

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  int sc;
  for (k = 1; k <= M; k++)
    {
      mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
      sc = ILogsum(ILogsum(mmx[0][k-1] + TSC(cp9O_MD,k-1),
			   imx[0][k-1] + TSC(cp9O_ID,k-1)),
		   dmx[0][k-1] + TSC(cp9O_DD,k-1));
      dmx[0][k] = sc;
    }
  /* We can do a full parse through all delete states. */
  erow[0]  = dmx[0][M] + TSC(cp9O_DM,M);
  /*****************************************************************
   * The main loop: scan the 'worst possible scoring seq', until we see a DP cell
   * with a score that falls below our thresh (if doing_align), or until
   * we know we can't ever fall below that thresh (if doing_scan)
   *****************************************************************/
  /* Recursion. */

  j = 0;
  keep_going = TRUE;
  while(keep_going)
    {
      int endsc     = -INFTY;
      int el_selfsc = cm->cp9->el_selfsc;
      int sc;
      j++;
      cur = j     % 2;
      prv = (j-1) % 2;
      /* The 1 difference between a Forward scanner and the 
       * regular Forward. In non-scanner parse must begin in B at
       * position 0 (i0-1), in scanner we can start at any position 
       * in the seq. */
      mmx[cur][0]  = (doing_scan == TRUE) ? 0 : -INFTY;
      dmx[cur][0]  = -INFTY;  /*D_0 is non-existent*/
      elmx[cur][0] = -INFTY;  /*no EL state for node 0 */

      sc = ILogsum(ILogsum(mmx[prv][0] + TSC(cp9O_MI,0),
			   imx[prv][0] + TSC(cp9O_II,0)),
		   dmx[prv][0] + TSC(cp9O_DI,0));
      imx[cur][0] = sc + min_isc[0];

      for (k = 1; k <= M; k++)
	{
	  /*match state*/
	  //ctr++;
	  sc = ILogsum(ILogsum(mmx[prv][k-1] + TSC(cp9O_MM,k-1),
			       imx[prv][k-1] + TSC(cp9O_IM,k-1)),
		       ILogsum(dmx[prv][k-1] + TSC(cp9O_DM,k-1),
			       mmx[prv][0]   + TSC(cp9O_BM,k  )));
	  /* check possibility we came from an EL, if they're valid */
	  for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	    sc = ILogsum(sc, elmx[prv][cm->cp9->el_from_idx[k][c]]);
	    /* transition penalty to EL incurred when EL was entered */
	  mmx[cur][k] = sc + min_msc[k];

	  /* E state update */
	  endsc = ILogsum(endsc, mmx[cur][k] + TSC(cp9O_ME,k));

	  /*insert state*/
	  sc = ILogsum(ILogsum(mmx[prv][k] + TSC(cp9O_MI,k),
			       imx[prv][k] + TSC(cp9O_II,k)),
		       dmx[prv][k] + TSC(cp9O_DI,k));
	  imx[cur][k] = sc + min_isc[k];

	  /*delete state*/
	  sc = ILogsum(ILogsum(mmx[cur][k-1] + TSC(cp9O_MD,k-1),
			       imx[cur][k-1] + TSC(cp9O_ID,k-1)),
		       dmx[cur][k-1] + TSC(cp9O_DD,k-1));
	  dmx[cur][k] = sc;

	  /*el state*/
	  sc = -INFTY;
	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								    HMM nodes that map to right half of a MATP_MP) */
	    {
	      sc = ILogsum(mmx[cur][k]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
			   elmx[prv][k] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	    }
	  elmx[cur][k] = sc;

	  minsc = ESL_MIN(minsc, mmx[cur][k]);
	  minsc = ESL_MIN(minsc, imx[cur][k]);
	  minsc = ESL_MIN(minsc, dmx[cur][k]);
	  if((cm->cp9->flags & CPLAN9_EL) && (cm->cp9->has_el[k])) 
	    minsc = ESL_MIN(minsc, elmx[cur][k]); 
	  if(cm->cp9->flags & CM_LOCAL_BEGIN)
	    minsc = ESL_MIN(minsc, endsc);

	  /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
	    printf("imx [jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
	    printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);
	    printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[cur][k]);*/
	}
      endsc = ILogsum(ILogsum(endsc, dmx[cur][M] + TSC(cp9O_DM,M)), /* transition from D_M -> end */
		      imx[cur][M] + TSC(cp9O_IM,M)); /* transition from I_M -> end */
      for(c = 0; c < cm->cp9->el_from_ct[M+1]; c++) /* el_from_ct[k] is >= 0 */
	endsc = ILogsum(endsc, elmx[cur][cm->cp9->el_from_idx[M+1][c]]);
      if(cm->cp9->flags & CM_LOCAL_BEGIN)
	minsc = ESL_MIN(minsc, endsc);
      
      /* check if we should break the loop now */
      if(doing_scan && j == 2) {
	keep_going = FALSE;
	ESL_DPRINTF1(("cp9_WorstForward() SCAN final minsc: %d j: %d\n", minsc, j));
	if  (minsc > thresh) retL = -1; /* no seq of any length will give minsc < thresh */
	else                 retL = j;
      }
      else if((doing_align) && (minsc <= thresh || j >= (100. * cm->clen))) 
	{
	  ESL_DPRINTF1(("cp9_WorstForward() ALIGN final minsc: %d thresh: %d clen: %d j: %d\n", minsc, thresh, cm->clen, j));
	  keep_going = FALSE;
	  retL = j;
	}
    } /* end of while(keep_going) */
  free(min_isc);
  free(min_msc);
  return retL;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}

/* cp9_GetLocalityMode()
 *
 * Infer a CM plan 9 HMM's locality mode based on it's
 * flagsl
 */
int
cp9_GetLocalityMode(CP9_t *cp9)
{
  if((cp9->flags & CPLAN9_LOCAL_BEGIN)     && (! (cp9->flags & CPLAN9_LOCAL_END)))
    cm_Fail("ERROR in cp9_GetLocalityMode, CPLAN9_LOCAL_BEGIN flag up, but CPLAN9_LOCAL_END flag down, this shouldn't happen.");
  if((! (cp9->flags & CPLAN9_LOCAL_BEGIN)) && (cp9->flags & CPLAN9_LOCAL_END))
    cm_Fail("ERROR in cp9_GetLocalityMode, CPLAN9_LOCAL_BEGIN flag down, but CPLAN9_LOCAL_END flag up, this shouldn't happen.");

  if((cp9->flags & CPLAN9_LOCAL_BEGIN) && (cp9->flags & CPLAN9_LOCAL_END))
    if(cp9->flags & CPLAN9_EL) return CP9_LOCAL_BEGIN_END_ON_AND_EL_ON;
    else                       return CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF;
  else
    if(cp9->flags & CPLAN9_EL) return CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON;
    else                       return CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF;
}

/* cp9_CheckTransitionGuarantees()
 *
 * Given a cm with a CP9 HMM and it's locality mode, check to make
 * sure it's transition guarantees hold up.
 * 
 * Return eslEINVAL if they don't (which we can usually deal with
 * without dying), or die if some of the local transitions don't 
 * hold up, which should absolutely never, ever happen.
 */
int
cp9_CheckTransitionGuarantees(CP9_t *cp9)
{
  int status = eslEINVAL;
  int locality_mode = cp9_GetLocalityMode(cp9);
  int k;

  /* these should be non -INFTY for all locality modes */
    for(k = 1; k < cp9->M; k++)
      {
	if(cp9->tsc[CTMM][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTMI][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTMD][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTIM][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTII][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTID][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTDM][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTDI][k] == -INFTY) goto ERROR;
	if(cp9->tsc[CTDD][k] == -INFTY) goto ERROR;
      }
    switch(locality_mode) {
    case CP9_LOCAL_BEGIN_END_ON_AND_EL_ON:
      for(k = 2; k < cp9->M; k++)
	{
	  if(cp9->bsc[k] == -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_ON_AND_EL_ON, but bsc[%d] is -INFTY!", k);
	  if(cp9->esc[k] == -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_ON_AND_EL_ON, but esc[%d] is -INFTY!", k);
	  if(cp9->has_el[k] && cp9->tsc[CTMEL][k] == -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_ON_AND_EL_ON, k:%d has an EL state, but tsc[CTMEL] is -INFTY!", k);
	}
      break;
    case CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF:
      for(k = 2; k < cp9->M; k++)
	{
	  if(cp9->bsc[k] == -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF, but bsc[%d] is -INFTY!", k);
	  if(cp9->esc[k] == -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF, but esc[%d] is -INFTY!", k);
	  if(cp9->tsc[CTMEL][k] != -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF, but tsc[%d][CTMEL] is ! -INFTY!", k);
	}
      break;
      break;
    case CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON:
      for(k = 2; k < cp9->M; k++)
	{
	  if(cp9->bsc[k] != -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON, but bsc[%d] is ! -INFTY!", k);
	  if(cp9->esc[k] != -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON, but esc[%d] is ! -INFTY!", k);
	  if(cp9->has_el[k] && cp9->tsc[CTMEL][k] == -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON, k:%d has an EL state, but tsc[CTMEL] is -INFTY!", k);
	}
      break;
    case CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF:
      for(k = 2; k < cp9->M; k++)
	{
	  if(cp9->bsc[k] != -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF, but bsc[%d] is ! -INFTY!", k);
	  if(cp9->esc[k] != -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF, but esc[%d] is ! -INFTY!", k);
	  if(cp9->tsc[CTMEL][k] != -INFTY) cm_Fail("CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF, but tsc[%d][CTMEL] is ! -INFTY!", k);
	}
      break;
    }
    return eslOK;

 ERROR:
    return status;
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
  { "-f",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute optimized Forward scan implementation",  0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute experimental Forward scan implementation", 0},
  { "-o",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute old (version 0.8) Forward scan implementation", 0},
  { "-w",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute slow Viterbi scan implementation",  0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "configure HMM for glocal alignment [default: local]", 0 },
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "do alignment, don't scan", 0 },
  { "--noel",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "turn local ends off [default: on, unless -g]", 0 },
  { "--full",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "calculate full matrix, not just 2 rows",         0 },
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
  int             do_scan;
  int             do_align;
  int             minL = 0;
  int             be_safe;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  //cm->config_opts |= CM_CONFIG_QDB;
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
  }
  if(! esl_opt_GetBoolean(go, "--noel")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  ConfigCM(cm, NULL, NULL);
  init_ilogsum();

  if (esl_opt_GetBoolean(go, "-a"))  { do_scan = FALSE; do_align = TRUE;  }
  else                               { do_scan = TRUE;  do_align = FALSE; }

  if (esl_opt_GetBoolean(go, "-x")) { 
    /* determine the minimum length we can search safely with the optimized forward version. */
    minL = cp9_WorstForward(cm, -INFTY, do_scan, do_align);
    /* minL = 100000; */
  }

  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, cm->null, abc->K, L, dsq);

      esl_stopwatch_Start(w);
      sc = cp9_FastViterbi(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			   do_scan,   /* are we scanning? */
			   do_align,  /* are we aligning? */
			   (! esl_opt_GetBoolean(go, "--full")),  /* memory efficient ? */
			   NULL,   /* don't want the DP matrix back */
			   NULL);  /* don't want traces back */
      printf("%4d %-30s %10.4f bits ", (i+1), "cp9_FastViterbi(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      
      if (esl_opt_GetBoolean(go, "-f")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = cp9_FastForward(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			       do_scan,   /* are we scanning? */
			       do_align,  /* are we aligning? */
			       FALSE,  /* we're not rescanning */
			       (! esl_opt_GetBoolean(go, "--full")),  /* memory efficient ? */
			       NULL);  /* don't want the DP matrix back */
	  printf("%4d %-30s %10.4f bits ", (i+1), "cp9_FastForward(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "-o")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = CP9Forward(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			  do_scan,   /* are we scanning? */
			  do_align,  /* are we aligning? */
			  FALSE,  /* we're not rescanning */
			  (! esl_opt_GetBoolean(go, "--full")),  /* memory efficient ? */
			  NULL);  /* don't want the DP matrix back */
	  printf("%4d %-30s %10.4f bits ", (i+1), "CP9Forward(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}
      if (esl_opt_GetBoolean(go, "-x")) 
	{ 
	  be_safe = FALSE;
	  ESL_DPRINTF1(("minL: %d L: %d\n", minL, L));
	  if(minL != -1 && minL <= L) be_safe = TRUE;
	  esl_stopwatch_Start(w);
	  sc = cp9_EXPTLFastForward(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
				    do_scan,   /* are we scanning? */
				    do_align,  /* are we aligning? */
				    FALSE,  /* we're not rescanning */
				    (! esl_opt_GetBoolean(go, "--full")),  /* memory efficient ? */
				    be_safe,
				    NULL);  /* don't want the DP matrix back */
	  printf("%4d %-30s %10.4f bits ", (i+1), "cp9_EXPTLFastForward(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "-w")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = CP9Viterbi(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			  do_scan,   /* are we scanning? */
			  do_align,  /* are we aligning? */
			  (! esl_opt_GetBoolean(go, "--full")),  /* memory efficient ? */
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

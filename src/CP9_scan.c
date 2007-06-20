/* CP9_scan.c 
 * EPN
 * 
 * Scanning algorithms for CM Plan 9 HMMs.  These algorithms align
 * subsequences of the target sequence to the model (e.g. glocal or
 * local alignment) Global alignment algorithms are in hmmband.c.
 *
 * These functions are still under development and are very
 * experimental. No guarantees here.
 *
 *################################################################
 * CP9Viterbi()     - Viterbi algorithm, in scan mode: scan input 
 *                    sequence for high scoring Viterbi hits to 
 *                    the model.
 * CP9Forward()     - Forward algorithm, in scan mode: scan input 
 *                    sequence for high scoring Forward hits to 
 *                    the model.
 * CP9Backward()    - Backward algorithm, in scan mode: scan input
 *                    sequence for high scoring Backward hits to 
 *                    the model.
 * CP9FilteredScan()- use CP9Forward() and potentially CP9Backward()
 *                    (if CM_SEARCH_HMMFB) to filter DB, 
 *                    and pass promising subseqs to CM search,
 *                    can be run in different modes.
 *################################################################
 * 
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

#include "stopwatch.h"          /* squid's process timing module        */
#include "hmmband.h"
#include "cm_dispatch.h"
#include "mpifuncs.h"
#include "stats.h"
#include "easel.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"
#include "esl_histogram.h"
#include "esl_normal.h"
#include "esl_stats.h"
#include "esl_gev.h"
#include "esl_weibull.h"
#include "esl_stretchexp.h"
#include "esl_gamma.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"
  
/* Function: CP9Viterbi()
 * based on  P7Viterbi() <-- this function's comments below  
 *           from HMMER 2.4devl core_algorithms.c
 *
 * Purpose:  Runs the Viterbi dynamic programming algorithm on an
 *           input subsequence (i0-j0). 
 *           Somewhat flexible based on input options as follows:
 *
 *           if(be_efficient): only allocates 2 rows of the Viterbi
 *           matrix, else allocates full L+1 matrix.
 *
 *           if(do_scan): allows parses to start at any position i
 *           i0..j0, changing meaning of DP matrix cells as discussed
 *           below.
 *
 *           Reference for algorithm (when do_scan is FALSE): 
 *           Durbin et. al. Biological Sequence Analysis; p. 58.
 *
 *           The meaning of the Viterbi (V) matrix DP cells for
 *           matches (M) inserts (I) and deletes (D):
 *
 *           For relative subsequence positions ip = 0..L:
 *           For HMM nodes 1..M: 
 *           V->M[ip][k] : score of most likely parse emitting seq
 *                         from i0..ip that visit node k's match 
 *                         state, which emits posn ip
 *           F->I[ip][k] : score of most likely parse emitting seq from 
 *                         i0..ip that visit node k's insert
 *                         state, which emits posn ip 
 *           F->D[ip][k] : score of most likely parse emitting seq from 
 *                         i0..ip that visit node k's delete
 *                         delete state, last emitted (leftmost)
 *                         posn was ip
 *
 *           For *special* HMM node 0:
 *           F->M[ip][0] : M_0 is the Begin state, which does not 
 *                         emit, so this is the score of the most likely
 *                         parse emitting seq from i0..ip that starts
 *                         in the begin state, the last emitted 
 *                         (leftmost) posn was ip.
 *
 *           Note: if ip=0, only D_k and M_0 states can have 
 *                 non-IMPOSSIBLE values. 
 *
 *           if(do_scan) the 'i0..ip' in the above definitions is
 *           changed to iE..ip such that i0 <= iE <= ip. Meaning
 *           any residue can be the first residue emitted in the
 *           parse. This means F->M[ip][0] is the score of the best
 *           parse emitting a subseq starting anywhere from i0..ip and 
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
 *           ret_sc    - RETURN: int log odds Viterbi score for each end point [0..(j0-i0+1)]
 *           ret_maxres- RETURN: start position that gives maximum score max argmax_i sc[i]
 *           results   - scan_results_t to add to; if NULL, don't keep results
 *           do_scan   - TRUE if we're scanning, HMM can start to emit anywhere i0..j0,
 *                       FALSE if we're not, HMM must start emitting at i0, end emitting at j0
 *           be_efficient- TRUE to keep only 2 rows of DP matrix in memory, FALSE keep whole thing
 *           ret_mx    - RETURN: CP9 Viterbi DP matrix, NULL if not wanted
 *           ret_tr    - RETURN: CP9 traces, NULL if not wanted
 *
 * Returns:  if(!do_scan) log P(S|M)/P(S|R), as a bit score
 *           else         max log P(S|M)/P(S|R), for argmax subseq S of input seq i0..j0,
 */
float
CP9Viterbi(CM_t *cm, char *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
	   int *ret_bestpos, scan_results_t *results, int do_scan,
	   int be_efficient, CP9_dpmatrix_t **ret_mx, CP9trace_t **ret_tr)
{
  CP9trace_t  *tr;
  CP9_dpmatrix_t *mx;       /* the CP9 DP matrix                                            */
  int **mmx;
  int **imx;
  int **dmx;
  int **emx;
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int   i,k;
  int         sc;
  int   ip;		/* i': relative position in the subsequence  */
  int          L;           /* j0-i0+1: subsequence length                                  */
  float        fsc;         /* float log odds score                                         */
  float        curr_sc;     /* temporary score used for filling in gamma                    */
  float       *gamma;       /* SHMM DP matrix for optimum nonoverlap resolution [0..j0-i0+1]*/
  int         *gback;       /* traceback pointers for SHMM                      [0..j0-i0+1]*/ 
  float       *savesc;      /* saves score of hit added to best parse at j      [0..j0-i0+1]*/ 
  float        best_hmm_sc; /* Best overall score from semi-HMM to return if do_scan        */
  float        best_hmm_pos;/* residue giving best_hmm_sc                                   */
  float        best_sc;     /* Best score overall, returned if 0 hits found by HMM & do_scan*/
  float        best_pos;    /* residue giving best_sc                                       */
  float        return_sc;   /* score to return, if (!do_scan) return overall Forward sc,    *
			     * else return best_hmm_sc if # HMM hits>0, else return best_sc */
  float        temp_sc;     /* temporary score                                              */
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */

  W  = j0-i0+1;		/* the length of the subsequence */

  /*printf("in CP9Forward() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    Die("ERROR in CP9Forward, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    Die("ERROR in CP9Forward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    Die("ERROR in CP9Forward, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    Die("ERROR in CP9Forward, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
    
  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_hmm_sc = IMPOSSIBLE;
  best_hmm_pos= -1;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+1));
  gamma[0] = 0.;
  gback    = MallocOrDie(sizeof(int)   * (L+1));
  gback[0] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+1));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 2;
  else             nrows = L+1;
  mx = AllocCPlan9Matrix(nrows, cm->cp9->M, &mmx, &imx, &dmx, NULL,
			 NULL);
  /* sc will hold P(seq up to j | Model) in int log odds form */
  /*sc   = MallocOrDie(sizeof(int) * (j0-i0+2));*/

  /* ResizeCPlan9Matrix(mx, W, hmm->M, &mmx, &imx, &dmx, &emx);*/
  /* Initialization of the zero row.
   */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= cm->cp9->M; k++)
    mmx[0][k] = imx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= cm->cp9->M; k++)
    {
      dmx[0][k]  = -INFTY;
      if((sc = mmx[0][k-1] + cm->cp9->tsc[CTMD][k-1]) > dmx[0][k])
	dmx[0][k] = sc;
      if((sc = dmx[0][k-1] + cm->cp9->tsc[CTDD][k-1]) > dmx[0][k])
	dmx[0][k] = sc;
    }

  /* We can do a full parse through all delete states. */
  emx[0][0] = dmx[0][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]; 
  fsc = Scorify(emx[0][0]);
  /*printf("jp: %d j: %d fsc: %f isc: %d\n", jp, j, fsc, isc[jp]);*/
  if(fsc > best_sc) 
    {
      best_sc = fsc;
      best_pos= i0-1;
    }
  
  /*****************************************************************
   * The main loop: scan the sequence from position i0 to j0.
   *****************************************************************/

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   */
  for (j = i0; j <= j0; j++)
    {
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
      if(do_scan)
	mmx[cur][0] = 0;
      else
	mmx[cur][0] = -INFTY;

      dmx[cur][0] = -INFTY;  /*D_0 is non-existent*/
      imx[cur][0] = -INFTY;
      if((sc = mmx[prv][0] + cm->cp9->tsc[CTMI][0]) > imx[cur][0])
	imx[cur][0] = sc;
      if((sc = imx[prv][0] + cm->cp9->tsc[CTII][0]) > imx[cur][0])
	imx[cur][0] = sc;
      if((sc = dmx[prv][0] + cm->cp9->tsc[CTDI][0]) > imx[cur][0])
	imx[cur][0] = sc;
      if(imx[cur][0] != -INFTY)
	imx[cur][0] += cm->cp9->isc[(int) dsq[j]][0];
      else 
	imx[cur][0] = -INFTY;

      for (k = 1; k <= cm->cp9->M; k++)
	{
	  /*match state*/
	  mmx[cur][k] = -INFTY;
	  if((sc = mmx[prv][k-1] + cm->cp9->tsc[CTMM][k-1]) > mmx[cur][k])
	    mmx[cur][k] = sc;
	  if((sc = imx[prv][k-1] + cm->cp9->tsc[CTIM][k-1]) > mmx[cur][k])
	    mmx[cur][k] = sc;
	  if((sc = mmx[prv][0] + cm->cp9->bsc[k]) > mmx[cur][k])
	    mmx[cur][k] = sc;
	  if((sc = dmx[prv][k-1] + cm->cp9->tsc[CTDM][k-1]) > mmx[cur][k])
	    mmx[cur][k] = sc;
	  if(mmx[cur][k] != -INFTY)
	    mmx[cur][k] += cm->cp9->msc[(int) dsq[j]][k];
	  else 
	    mmx[cur][k] = -INFTY;

	  /*insert state*/
	  imx[cur][k] = -INFTY;
	  if((sc = mmx[prv][k] + cm->cp9->tsc[CTMI][k]) > imx[cur][k])
	    imx[cur][k] = sc;
	  if((sc = imx[prv][k] + cm->cp9->tsc[CTII][k]) > imx[cur][k])
	    imx[cur][k] = sc;
	  if((sc = dmx[prv][k] + cm->cp9->tsc[CTDI][k]) > imx[cur][k])
	    imx[cur][k] = sc;
	  if(imx[cur][k] != -INFTY)
	    imx[cur][k] += cm->cp9->isc[(int) dsq[j]][k];
	  else 
	    imx[cur][k] = -INFTY;

	  /*delete state*/
	  dmx[cur][k] = -INFTY;
	  if((sc = mmx[cur][k-1] + cm->cp9->tsc[CTMD][k-1]) > dmx[cur][k])
	    dmx[cur][k] = sc;
	  if((sc = imx[cur][k-1] + cm->cp9->tsc[CTID][k-1]) > dmx[cur][k])
	    dmx[cur][k] = sc;
	  if((sc = dmx[cur][k-1] + cm->cp9->tsc[CTDD][k-1]) > dmx[cur][k])
	    dmx[cur][k] = sc;
	}
      emx[0][cur] = -INFTY;
      for (k = 1; k <= cm->cp9->M; k++)
	if ((sc = mmx[cur][k] + cm->cp9->esc[k]) > emx[0][cur])
	  emx[0][cur] = sc;
      if ((sc =  dmx[cur][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]) > emx[0][cur])
	emx[0][cur] = sc;
      /* transition from D_M -> end */
      fsc = Scorify(emx[0][cur]);
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
	  curr_sc = gamma[prv] + fsc + cm->cp9_sc_boost;
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
		  /*printf("FWD greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i, j, fsc);*/
		  report_hit (i, j, 0, fsc, results);
		  /* 0 is for saver, which is irrelevant for HMM hits */
		}
	    }
	}
    } /* end loop over end positions j */
  sc = emx[0][W];
  /*printf("returing sc: %d from CPViterbi()\n", sc);*/
  if((!(cm->search_opts & CM_SEARCH_HMMGREEDY)) && /* resolve overlaps optimally */
     do_scan) /* else we can save time by skipping traceback */
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
		    printf("VITERBI reporting hit: i: %d j: %d sc: %f\n", gback[jp], j, savesc[jp]);
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
  if(!do_scan)
    {
      return_sc = Scorify(emx[0][cur]); /* probably wrong */
      if(ret_bestpos != NULL) *ret_bestpos = i0;
    }
  
  if (!do_scan && ret_tr != NULL) {
    CP9ViterbiTrace(cm->cp9, dsq, i0, j0, mx, &tr);
    /* CP9PrintTrace(stdout, tr, cm->cp9, dsq); */
    *ret_tr = tr;
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
  /*printf("Forward return_sc: %f\n", return_sc);*/
  return return_sc;
}

/***********************************************************************
 * Function: CP9Forward()
 * 
 * Purpose:  Runs the Forward dynamic programming algorithm on an
 *           input subsequence (i0-j0). Complements CP9Backward().  
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
 *           results   - scan_results_t to add to; if NULL, don't keep results
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
CP9Forward(CM_t *cm, char *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
	   int *ret_bestpos, scan_results_t *results, int do_scan, int doing_align, int doing_rescan, 
	   int be_efficient, CP9_dpmatrix_t **ret_mx)
{
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
  float        return_sc;   /* score to return, if (!do_scan) return overall Forward sc,    *
			     * else return best_hmm_sc if # HMM hits>0, else return best_sc */
  int          accept;      /* Flag used if cm->search_opts & CM_SEARCH_HMMRESCAN           */
  float        temp_sc;     /* temporary score                                              */
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */

  /*debug_print_cp9_params(cm->cp9);*/

  /*printf("in CP9Forward() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    Die("ERROR in CP9Forward, cm->cp9 is NULL.\n");
  if(doing_rescan && !do_scan) 
    Die("ERROR in CP9Forward, doing_rescan but not do_scan");
  if(be_efficient && (ret_mx != NULL))
    Die("ERROR in CP9Forward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    Die("ERROR in CP9Forward, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    Die("ERROR in CP9Forward, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
    
  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_hmm_sc = IMPOSSIBLE;
  best_hmm_pos= -1;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+1));
  gamma[0] = 0.;
  gback    = MallocOrDie(sizeof(int)   * (L+1));
  gback[0] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+1));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 2;
  else             nrows = L+1;
  mx = AllocCPlan9Matrix(nrows, cm->cp9->M, &mmx, &imx, &dmx, NULL, /*FIX ME! */
			 NULL); /* we don't use emx any more, isc replaces it,
				 * 1. b/c we want isc for each position, but
				 *    we may only be alloc'ing 2 rows of mx,
				 * 2. b/c emx is 2D for dumb memory error reasons,
				 *    which needs to be fixed. */
  /* sc will hold P(seq up to j | Model) in int log odds form */
  sc   = MallocOrDie(sizeof(int) * (j0-i0+2));
			
  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  /*printf("mmx[jp:%d][%d]: %d %f\n", 0, 0, mmx[0][0], Score2Prob(mmx[0][0], 1.));
    printf("imx[jp:%d][%d]: %d %f\n", 0, 0, imx[0][0], Score2Prob(imx[0][0], 1.));
    printf("dmx[jp:%d][%d]: %d %f\n", 0, 0, dmx[0][0], Score2Prob(dmx[0][0], 1.));*/

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= cm->cp9->M; k++)
    {
      mmx[0][k] = imx[0][k] = -INFTY;      /* need seq to get here */
      dmx[0][k] = ILogsum(ILogsum(mmx[0][k-1] + cm->cp9->tsc[CTMD][k-1],
				  imx[0][k-1] + cm->cp9->tsc[CTID][k-1]),
			  dmx[0][k-1] + cm->cp9->tsc[CTDD][k-1]);
      /*printf("mmx[jp:%d][%d]: %d %f\n", 0, k, mmx[0][k], Score2Prob(mmx[0][k], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", 0, k, imx[0][k], Score2Prob(imx[0][k], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));*/
    }
  /* We can do a full parse through all delete states. */
  sc[0] = dmx[0][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]; 
  fsc = Scorify(sc[0]);
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
      jp = j-i0+1;     /* jp is relative position in the sequence 1..L */
      cur = (j-i0+1);
      prv = (j-i0);
      if(be_efficient)
	{
	  cur %= 2;
	  prv %= 2;
	}	  
      /* The 1 difference between a Forward scanner and the 
       * regular Forward. In non-scanner parse must begin in B at
       * position 0 (i0-1), in scanner we can start at any position 
       * in the seq. */
      if(do_scan)
	mmx[cur][0] = 0;
      else
	mmx[cur][0] = -INFTY;

      dmx[cur][0] = -INFTY;  /*D_0 is non-existent*/
      imx[cur][0]  = ILogsum(ILogsum(mmx[prv][0] + cm->cp9->tsc[CTMI][0],
				     imx[prv][0] + cm->cp9->tsc[CTII][0]),
			     dmx[prv][0] + cm->cp9->tsc[CTDI][0]);
      imx[cur][0] += cm->cp9->isc[(int) dsq[j]][0];
      /*printf("mmx[jp:%d][%d]: %d %f\n", jp, 0, mmx[cur][0], Score2Prob(mmx[cur][0], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", jp, 0, imx[cur][0], Score2Prob(imx[cur][0], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", jp, 0, dmx[cur][0], Score2Prob(dmx[cur][0], 1.));*/

      for (k = 1; k <= cm->cp9->M; k++)
	{
	  mmx[cur][k]  = ILogsum(ILogsum(mmx[prv][k-1] + cm->cp9->tsc[CTMM][k-1],
				       imx[prv][k-1] + cm->cp9->tsc[CTIM][k-1]),
				 ILogsum(mmx[prv][0] + cm->cp9->bsc[k],
					 dmx[prv][k-1] + cm->cp9->tsc[CTDM][k-1]));
	  mmx[cur][k] += cm->cp9->msc[(int) dsq[j]][k];
	  dmx[cur][k]  = ILogsum(ILogsum(mmx[cur][k-1] + cm->cp9->tsc[CTMD][k-1],
					imx[cur][k-1] + cm->cp9->tsc[CTID][k-1]),
				dmx[cur][k-1] + cm->cp9->tsc[CTDD][k-1]);
	  
	  imx[cur][k]  = ILogsum(ILogsum(mmx[prv][k] + cm->cp9->tsc[CTMI][k],
				       imx[prv][k] + cm->cp9->tsc[CTII][k]),
			       dmx[prv][k] + cm->cp9->tsc[CTDI][k]);
	  imx[cur][k] += cm->cp9->isc[(int) dsq[j]][k];
	  /*printf("mmx[jp:%d][%d]: %d %f\n", jp, k, mmx[cur][k], Score2Prob(mmx[cur][k], 1.));
	    printf("imx[jp:%d][%d]: %d %f\n", jp, k, imx[cur][k], Score2Prob(imx[cur][k], 1.));
	    printf("dmx[jp:%d][%d]: %d %f\n", jp, k, dmx[cur][k], Score2Prob(dmx[cur][k], 1.));*/
	}
      /* determine sc[jp], the int score of all possible parses ending at the current
       * position (j) of the target sequence. */
      sc[jp] = -INFTY;
      for (k = 1; k <= cm->cp9->M; k++)
	sc[jp] = ILogsum(sc[jp], mmx[cur][k] + cm->cp9->esc[k]);
      /* 04.17.07 Arent' I double counting here! (at least for scanner?) */
      sc[jp] = ILogsum(sc[jp], dmx[cur][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]); 
      sc[jp] = ILogsum(sc[jp], imx[cur][cm->cp9->M] + cm->cp9->tsc[CTIM][cm->cp9->M]); 
      /* transition from D_M -> end */
      fsc = Scorify(sc[jp]);
      /*printf("jp: %d j: %d fsc: %f sc: %d\n", jp, j, fsc, sc[jp]);*/
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
	  curr_sc = gamma[prv] + fsc + cm->cp9_sc_boost;
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
		  /*printf("FWD greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i, j, fsc);*/
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

  /* determine score to return: (I know, too complex) */
  if(doing_align)
    {
      return_sc = Scorify(sc[(j0-i0+1)]); /* L = j0-i0+1 */
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
  if(ret_sc != NULL) *ret_sc = sc;
  else free(sc);
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);
  /*printf("Forward return_sc: %f\n", return_sc);*/
  return return_sc;
}

/***********************************************************************
 * Function: CP9Backward()
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
 *           results   - scan_results_t to add to; if NULL, don't keep results
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
CP9Backward(CM_t *cm, char *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
	    int *ret_bestpos, scan_results_t *results, int do_scan, int doing_align, 
	    int doing_rescan, int be_efficient, CP9_dpmatrix_t **ret_mx)
{
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

  /*printf("in CP9Backward() i0: %d j0: %d do_scan: %d \n", i0, j0, do_scan);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    Die("ERROR in CP9Backward, cm->cp9 is NULL.\n");
  if(doing_rescan && !do_scan) 
    Die("ERROR in CP9Backward, doing_rescan but not do_scan");
  if(be_efficient && (ret_mx != NULL))
    Die("ERROR in CP9Backward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    Die("ERROR in CP9Backward, passing in results data structure, but not in scanning mode.a\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    Die("ERROR in CP9Backward, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
    
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
  gamma        = MallocOrDie(sizeof(float) * (L+2));
  gamma[0]     = 0.; /* this is impossible, no hit can be rooted at ip=0, 
		      * we require hits to have at least 1 emission,
		      * plus, no positive scoring hit can have 0 emissions as 
		      * all transition scores are negative */
  gamma[L]   = 0.;
  gback        = MallocOrDie(sizeof(int)   * (L+2));
  gback[L]   = -1;
  savesc       = MallocOrDie(sizeof(float) * (L+2));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 2;
  else             nrows = L+1; 
  mx = AllocCPlan9Matrix(nrows, cm->cp9->M, &mmx, &imx, &dmx, NULL, /*FIX ME!*/
			 NULL); /* we don't use emx any more, sc replaces it,
				 * 1. b/c we want isc for each position, but
				 * we may only be alloc'ing 2 rows of mx,
				 * 2. b/c emx is 2D for dumb memory error reasons,
				 * which needs to be fixed. */
  /* sc will hold P(seq from i..j0 | Model) for each i in int log odds form */
  sc   = MallocOrDie(sizeof(int) * (j0-i0+3));

  /* Initialization of the L (i = j0, cur = (j0-i) = (j0-j0) %2 = 0) row. */
  if(be_efficient) cur = 0;
  else cur = j0-i0+1; /* L */
  ip = j0-i0+1;  /*L */
  i  = j0;

  mmx[cur][cm->cp9->M]  = 0. + cm->cp9->esc[cm->cp9->M]; /* M_M<-E ... everything ends in E (the 0; 2^0=1.0) */
  mmx[cur][cm->cp9->M] += cm->cp9->msc[(int) dsq[i]][cm->cp9->M]; /* ... + emitted match symbol */
  imx[cur][cm->cp9->M]  = 0. + cm->cp9->tsc[CTIM][cm->cp9->M]; /* I_M<-E ... everything ends in E (the 0; 2^0=1.0) */
  imx[cur][cm->cp9->M] += cm->cp9->isc[(int) dsq[i]][cm->cp9->M]; /* ... + emitted insert symbol */
  dmx[cur][cm->cp9->M]  = cm->cp9->tsc[CTDM][cm->cp9->M];    /* D_M<-E everything ends in E (the 0; 2^0=1.0) */
  for (k = cm->cp9->M-1; k >= 1; k--)
    {
      mmx[cur][k]  = 0 + cm->cp9->esc[k];  /*M_k<- E */
      mmx[cur][k]  = ILogsum(mmx[cur][k], dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]);
      mmx[cur][k] += cm->cp9->msc[(int) dsq[i]][k];

      imx[cur][k]  = dmx[cur][k+1] + cm->cp9->tsc[CTID][k];
      imx[cur][k] += cm->cp9->isc[(int) dsq[i]][k];

      dmx[cur][k]  = dmx[cur][k+1] + cm->cp9->tsc[CTDD][k];
    }
  
  /* remember M_0 is special, the B state, a non-emitter */
  mmx[cur][0]  = dmx[cur][1] + cm->cp9->tsc[CTMD][0]; /* M_0(B)->D_1, no seq emitted, all deletes */
  /* above line is diff from CPBackwardOLD() which has mmx[cur][0] = -INFTY; */
  imx[cur][0]  = dmx[cur][1] + cm->cp9->tsc[CTID][0];
  imx[cur][0] += cm->cp9->isc[(int) dsq[i]][0];

  dmx[cur][0]  = -INFTY; /*D_0 doesn't exist*/

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
  for (i = j0-1; i >= (i0-1); i--)
    {
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

      /* Now the main states. Note the boundary conditions at M. */
      if(ip > 0)
	{
	  mmx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTMI][cm->cp9->M];
	  mmx[cur][cm->cp9->M] += cm->cp9->msc[(int) dsq[i]][cm->cp9->M];
	  imx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTII][cm->cp9->M];
	  imx[cur][cm->cp9->M] += cm->cp9->isc[(int) dsq[i]][cm->cp9->M];
	}
      else /* ip == 0 */
	{
	  mmx[cur][cm->cp9->M] = -INFTY; /* need seq to get here */
	  imx[cur][cm->cp9->M] = -INFTY; /* need seq to get here */
	}
      dmx[cur][cm->cp9->M]  = imx[prv][cm->cp9->M] + cm->cp9->tsc[CTDI][cm->cp9->M]; 
	      
      /* A main difference between a Backward scanner and 
       * regular Backward: a scanner can end at the END 
       * state at any position, regular can only end at
       * the final position j0. */
      if(do_scan) 
	{	
	  if(ip > 0)
	    {
	      mmx[cur][cm->cp9->M]  =  
		ILogsum(mmx[cur][cm->cp9->M],
			(cm->cp9->esc[cm->cp9->M] +                  /* M_M<-E + (only in scanner)     */ 
			 0));                                        /* all parses end in E, 2^0 = 1.0;*/
	      mmx[cur][cm->cp9->M] += cm->cp9->msc[(int) dsq[i]][cm->cp9->M]; /* ... + emitted match symbol */
	      
	      imx[cur][cm->cp9->M]  =
		ILogsum(imx[cur][cm->cp9->M],
			(cm->cp9->tsc[CTIM][cm->cp9->M] +            /* I_M<-E + (only in scanner)     */
			 0));                                        /* all parses end in E, 2^0 = 1.0;*/
	      imx[cur][cm->cp9->M] += cm->cp9->isc[(int) dsq[i]][cm->cp9->M]; /* ... + emitted insert symbol */  
	    }
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
	  if(ip > 0) 
	    {
	      mmx[cur][k]  = ILogsum(ILogsum((mmx[prv][k+1] + cm->cp9->tsc[CTMM][k]),  
					     (imx[prv][k]   + cm->cp9->tsc[CTMI][k])),
				     (dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]));
	      mmx[cur][k] += cm->cp9->msc[(int) dsq[i]][k];

	      imx[cur][k]  = ILogsum(ILogsum((mmx[prv][k+1] + cm->cp9->tsc[CTIM][k]),
					     (imx[prv][k]   + cm->cp9->tsc[CTII][k])),
				     (dmx[cur][k+1] + cm->cp9->tsc[CTID][k]));
	      imx[cur][k] += cm->cp9->isc[(int) dsq[i]][k];
	    }
	  else
	    {
	      mmx[cur][k] = -INFTY; /* need seq to get here, unless we come from E in a scanner (below) */
	      imx[cur][k] = -INFTY; /* need seq to get here */
	    }
	  if(do_scan && ip > 0) /* add possibility of ending at his position from this state */
	    {
	      mmx[cur][k] = 
		ILogsum(mmx[cur][k], 
			(cm->cp9->esc[k] +                    /* M_k<-E + (only in scanner)     */ 
			 0));                                 /* all parses end in E, 2^0 = 1.0;*/
	      /* DO NOT add contribution of emitting i from M, it's been added above */
	    }	      
	  dmx[cur][k]  = ILogsum(ILogsum((mmx[prv][k+1] + cm->cp9->tsc[CTDM][k]),
					 (imx[prv][k]   + cm->cp9->tsc[CTDI][k])),
				 (dmx[cur][k+1] + cm->cp9->tsc[CTDD][k]));
	}

      /* Case when k == 0 */
      /* imx[cur][0] is filled same as imx[cur][1..k] in the loop above */
      if(ip > 0)
	{
	  imx[cur][0] = ILogsum(ILogsum((mmx[prv][1] + cm->cp9->tsc[CTIM][0]),
					(imx[prv][0] + cm->cp9->tsc[CTII][0])),
				(dmx[cur][1] + cm->cp9->tsc[CTID][k]));
	  imx[cur][0] += cm->cp9->isc[(int) dsq[i]][k];
	}
      else /* ip == 0 */
	imx[cur][0] = -INFTY; /* need seq to get here */
      dmx[cur][0] = -INFTY; /* D_0 does not exist */
	  
      /*M_0 is the B state, it doesn't emit, and can be reached from any match via a begin transition */
      mmx[cur][0] = -INFTY;
      for (k = cm->cp9->M; k >= 1; k--) 
	mmx[cur][0] = ILogsum(mmx[cur][0], (mmx[prv][k] + cm->cp9->bsc[k]));
      mmx[cur][0] = ILogsum(mmx[cur][0], (imx[prv][0] + cm->cp9->tsc[CTMI][0]));
      mmx[cur][0] = ILogsum(mmx[cur][0], (dmx[cur][1] + cm->cp9->tsc[CTMD][0]));     /* B->D_1 */
      
      /* determine isc, the int score of all possible parses starting at the current
       * position (i) of the target sequence. */
      sc[ip] = mmx[cur][0]; /* all parses must start in M_0, the B state */
      fsc = Scorify(sc[ip]);
      /*printf("ip: %d i: %d fsc: %f i: %d\n", ip, i, fsc, sc[ip]);*/
      if(fsc > best_sc)
	{
	  best_sc = fsc;
	  best_pos= i+1; /* *off-by-one* (see *off-by-one* below) */
	}
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
      if (fsc > best_hmm_sc)
	{
	  best_hmm_sc = fsc;
	  best_hmm_pos= i;
	}
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
		  accept = TRUE;
		  /* Potentially rescan just the subseq that is the hit we're about to report.
		   * Implemented to deal with fact that --enfseq option was enforcing the subseq
		   * to have hit pass filter b/c this Forward scanning function is 'infinite' length
		   * (Weinberg). Sometimes the subseq we're about to report has a really
		   * crappy score, even though the cumulative Backward score (ending at j0) is good. */
		  if(do_scan && cm->search_opts & CM_SEARCH_HMMRESCAN && doing_rescan == FALSE)
		    {
		      /*printf("rechecking hit from %d to %d\n", i, gback[ip]);*/
		      temp_sc = CP9Backward(cm, dsq, i, gback[ip], cm->W, cutoff,  /* *off-by-one* i+1 */
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

  /*printf("returning from CP9Backward()\n");*/
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
}

/***********************************************************************
 * Function: CP9Scan_dispatch()
 * Incept:   EPN, Tue Jan  9 06:28:49 2007
 * 
 * Purpose:  Scan a sequence with a CP9, potentially rescan CP9 hits with CM.
 *
 *           3 possible modes:
 *
 *           Mode 1: Filter mode with default pad 
 *                   (IF cm->search_opts & CM_SEARCH_HMMFILTER and 
 *                      !cm->search_opts & CM_SEARCH_HMMPAD)
 *                   Scan with CP9Forward() to get likely endpoints (j) of 
 *                   hits, for each j do a CP9Backward() from j-W+1..j to get
 *                   the most likely start point i for this j. 
 *                   Set i' = j-W+1 and
 *                       j' = i+W-1. 
 *                   Each i'..j' subsequence is rescanned with the CM.
 * 
 *           Mode 2: Filter mode with user-defined pad 
 *                   (IF cm->search_opts & CM_SEARCH_HMMFILTER and 
 *                       cm->search_opts & CM_SEARCH_HMMPAD)
 *                   Same as mode 1, but i' and j' defined differently:
 *                   Set i' = i - (cm->hmmpad) and
 *                       j' = j + (cm->hmmpad) 
 *                   Each i'..j' subsequence is rescanned with the CM.
 * 
 *           Mode 3: HMM only mode (IF cm->search_opts & CM_SEARCH_HMMONLY)
 *                   Hit boundaries i and j are defined the same as in mode 1, but
 *                   no rescanning with CM takes place. i..j hits are reported 
 *                   (note i' and j' are never calculated).
 * 
 * Args:     
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report (or keep if filtering)
 *           results    - scan_results_t to add to, only passed to 
 *                        actually_search_target()
 *           doing_cp9_stats- TRUE if we're calc'ing stats for the CP9, in this 
 *                            case we never rescan with CM
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc, score of maximally scoring end position j 
 */
float
CP9Scan_dispatch(CM_t *cm, char *dsq, int i0, int j0, int W, float cm_cutoff, 
		 float cp9_cutoff, scan_results_t *results, int doing_cp9_stats,
		 int *ret_flen)
{
  int h;
  int i;
  int min_i;
  float best_hmm_sc;
  float best_hmm_fsc;
  float cur_best_hmm_bsc;
  float best_cm_sc;
  int   flen;
  float ffrac;
  int do_collapse;
  int ipad;
  int jpad;
  int padmode;
  scan_results_t *fwd_results;
  scan_results_t *bwd_results;

  /* check contract */
  if(cm->cp9 == NULL)
    Die("ERROR in CP9Scan_dispatch(), cm->cp9 is NULL\n");
  if((cm->search_opts & CM_SEARCH_HMMPAD) &&
     (!(cm->search_opts & CM_SEARCH_HMMFILTER)))
     Die("ERROR in CP9Scan_dispatch(), CM_SEARCH_HMMPAD flag up, but CM_SEARCH_HMMFILTER flag down.\n");
  if(!doing_cp9_stats && (!((cm->search_opts & CM_SEARCH_HMMFILTER) || 
			    (cm->search_opts & CM_SEARCH_HMMONLY))))
    Die("ERROR in CP9Scan_dispatch(), not doing CP9 stats and neither CM_SEARCH_HMMFILTER nor CM_SEARCH_HMMONLY flag is up.\n");

  /*printf("in CP9Scan_dispatch(), i0: %d j0: %d\n", i0, j0);
    printf("cp9_cutoff: %f\n", cp9_cutoff);*/

  best_cm_sc = best_hmm_sc = IMPOSSIBLE;
  /* set up options for RescanFilterSurvivors() if we're filtering */
  if(cm->search_opts & CM_SEARCH_HMMFILTER)
    {
      if(cm->search_opts & CM_SEARCH_HMMPAD) /* mode 2 */
	{
	  padmode = PAD_SUBI_ADDJ;
	  ipad = jpad = cm->hmmpad; /* subtract hmmpad from i, add hmmpad to j */
	}
      else /* mode 1 */
	{
	  padmode = PAD_ADDI_SUBJ;
	  ipad = jpad = W-1; /* subtract W-1 from j, add W-1 to i */
	}
      if(cm->search_opts && CM_SEARCH_HBANDED)
	do_collapse = FALSE;
      else
	do_collapse = TRUE;
    }
  
  /* Scan the (sub)seq w/Forward, getting j end points of hits above cutoff */
  fwd_results = CreateResults(INIT_RESULTS);
  best_hmm_fsc = CP9Forward(cm, dsq, i0, j0, W, cp9_cutoff, NULL, NULL, fwd_results,
			    TRUE,   /* we're scanning */
			    FALSE,  /* we're not ultimately aligning */
			    FALSE,  /* we're not rescanning */
			    TRUE,   /* be memory efficient */
			    NULL);  /* don't want the DP matrix back */
  best_hmm_sc = best_hmm_fsc;

  /* Remove overlapping hits, if we're being greedy */
  if(cm->search_opts & CM_SEARCH_HMMGREEDY) /* resolve overlaps by being greedy */
    {
      assert(i0 == 1); 
      remove_overlapping_hits (fwd_results, i0, j0);
    }

  /* Determine start points (i) of the hits based on Backward scan starting at j,
   * report hits IFF CM_SEARCH_HMMONLY */
  bwd_results = CreateResults(INIT_RESULTS);
  for(h = 0; h < fwd_results->num_results; h++) 
    {
      min_i = (fwd_results->data[h].stop - W + 1) >= 1 ? (fwd_results->data[h].stop - W + 1) : 1;
      cur_best_hmm_bsc = CP9Backward(cm, dsq, min_i, fwd_results->data[h].stop, W, cp9_cutoff, 
				     NULL, /* don't care about score of each posn */
				     &i,   /* set i as the best scoring start point from j-W..j */
				     ((cm->search_opts & CM_SEARCH_HMMONLY) ? results : bwd_results),  
				     TRUE,  /* we're scanning */
				     /*FALSE,*/  /* we're not scanning */
				     FALSE, /* we're not ultimately aligning */
				     FALSE, /* don't rescan */
				     TRUE,  /* be memory efficient */
				     NULL); /* don't want the DP matrix back */
      //FALSE,  /* don't be memory efficient */
      //&bmx); /* give the DP matrix back */
      /* this only works if we've saved the matrices, and didn't do scan mode
       * for both Forward and Backward:
       * debug_check_CP9_FB(fmx, bmx, cm->cp9, cur_best_hmm_bsc, i0, j0, dsq); */
      
      if(cur_best_hmm_bsc > best_hmm_sc) best_hmm_sc = cur_best_hmm_bsc;
      /*printf("cur_best_hmm_bsc: %f\n", cur_best_hmm_bsc);*/
    }	  
  /* Rescan with CM if we're filtering and not doing cp9 stats */
  if(!doing_cp9_stats && (cm->search_opts & CM_SEARCH_HMMFILTER))
    {
      /* Remove overlapping hits, if we're being greedy */
      if(cm->search_opts & CM_SEARCH_HMMGREEDY) 
	{
	  assert(i0 == 1); 
	  remove_overlapping_hits (bwd_results, i0, j0);
	}
      best_cm_sc = RescanFilterSurvivors(cm, dsq, bwd_results, i0, j0, W, 
					 padmode, ipad, jpad, 
					 do_collapse, cm_cutoff, cp9_cutoff, 
					 results, &flen);
      if(flen == 0) ffrac = 100.;
      else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
      /*if(!(cm->search_opts & CM_SEARCH_HMMONLY))
	printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);*/
    }
  FreeResults (fwd_results);
  FreeResults (bwd_results);

  /*printf("in CP9Scan_dispatch, returning best_hmm_sc: %f\n", best_hmm_sc);*/
  if(doing_cp9_stats || cm->search_opts & CM_SEARCH_HMMONLY)
    return best_hmm_sc;
  else
    return best_cm_sc;
}

/***********************************************************************
 * Function: RescanFilterSurvivors()
 * Incept:   EPN, Wed Apr 11 05:51:55 2007
 * 
 * Purpose:  Given start and end points of hits that have survived
 *           a CP9 filter, pad a given amount of residues on 
 *           on each side, and rescan with a CM. Optionally,
 *           collapse overlapping subseqs into one larger subseq before
 *           rescanning (we don't always do this b/c we may want to
 *           derive HMM bands for a subseq from a Forward/Backward scan).
 * 
 *           Can be run in 2 modes, depending on input variable padmode:
 *           Mode 1: padmode = PAD_SUBI_ADDJ
 *                   For each i,j pair in hiti, hitj: 
 *                   set i' = i - ipad; and j' = j + jpad, 
 *                   Rescan subseq from i' to j'.
 *           Mode 2: padmode = PAD_ADDI_SUBJ
 *                   For each i,j pair in hiti, hitj: 
 *                   set j' = i + ipad; and i' = j - jpad, 
 *                   ensure j' >= j and i' <= i. 
 *                   Rescan subseq from i' to j'.
 *
 * Args:     
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           hmm_results- info on HMM hits that survived filter
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           padmode    - PAD_SUBI_ADDJ or PAD_ADDI_SUBJ (see above)
 *           ipad       - number of residues to subtract/add from each i 
 *           jpad       - number of residues to add/subtract from each j 
 *           do_collapse- TRUE: collapse overlapping hits (after padding) into 1
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report 
 *           results    - scan_results_t to add to, only passed to 
 *                        actually_search_target()
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc found when rescanning with CM 
 */
float
RescanFilterSurvivors(CM_t *cm, char *dsq, scan_results_t *hmm_results, int i0, int j0,
		      int W, int padmode, int ipad, int jpad, int do_collapse,
		      float cm_cutoff, float cp9_cutoff, scan_results_t *results, int *ret_flen)
{
  int h;
  int i, j;
  float best_cm_sc;
  float cm_sc;
  int   flen;
  int   prev_j;
  int   next_j;
  int   nhits;

  best_cm_sc = IMPOSSIBLE;
  flen = 0;

  /* check contract */
  if(padmode != PAD_SUBI_ADDJ && padmode != PAD_ADDI_SUBJ)
    ESL_EXCEPTION(eslEINCOMPAT, "can't determine mode.");

  /*if(padmode == PAD_SUBI_ADDJ)
    printf("in RescanFilterSurvivors(), mode: PAD_SUBI_ADDJ\n");
    else
    printf("in RescanFilterSurvivors(), mode: PAD_ADDI_SUBJ\n");
    printf("\tipad: %d, jpad: %d collapse: %d\n", ipad, jpad, do_collapse);*/

  /* For each hit, add pad according to mode and rescan by calling actually_search_target(). 
   * If do_collapse, collapse multiple overlapping hits into 1 before rescanning */
  /* hits should always be sorted by decreasing j, if this is violated - die. */
  nhits = hmm_results->num_results;
  for(h = 0; h < nhits; h++) 
    {
      if(h != 0 && hmm_results->data[h].stop > prev_j) 
	ESL_EXCEPTION(eslEINCOMPAT, "j's not in descending order");
      prev_j = hmm_results->data[h].stop;

      /* add pad */
      if(padmode == PAD_SUBI_ADDJ)
	{
	  i = ((hmm_results->data[h].start - ipad) >= 1)    ? (hmm_results->data[h].start - ipad) : 1;
	  j = ((hmm_results->data[h].stop  + jpad) <= j0)   ? (hmm_results->data[h].stop  + jpad) : j0;
	  if((h+1) < nhits)
	    next_j = ((hmm_results->data[h+1].stop + jpad) <= j0)   ? (hmm_results->data[h+1].stop + jpad) : j0;
	  else
	    next_j = -1;
	}
      else if(padmode == PAD_ADDI_SUBJ)
	{
	  i = ((hmm_results->data[h].stop  - jpad) >= 1)    ? (hmm_results->data[h].stop  - jpad) : 1;
	  j = ((hmm_results->data[h].start + ipad) <= j0)   ? (hmm_results->data[h].start + ipad) : j0;
	  if((h+1) < nhits)
	    next_j = ((hmm_results->data[h+1].start + ipad) <= j0)   ? (hmm_results->data[h+1].start + ipad) : j0;
	  else
	    next_j = -1;
	}
      /*printf("subseq: hit %d i: %d (%d) j: %d (%d)\n", h, i, hmm_results->data[h].start[h], j, hmm_results->data[h].stop[h]);*/

      if(do_collapse) /* collapse multiple hits that overlap after padding on both sides into a single hit */
	{
	  while(((h+1) < nhits) && (next_j >= i))
	    {
	      /* suck in hit */
	      h++;
	      if(padmode == PAD_SUBI_ADDJ)
		{
		  i = ((hmm_results->data[h].start - ipad) >= 1)    ? (hmm_results->data[h].start - ipad) : 1;
		  if((h+1) < nhits)
		    next_j = ((hmm_results->data[h+1].stop + jpad) <= j0)   ? (hmm_results->data[h+1].stop + jpad) : j0;
		  else
		    next_j = -1;
		}
	      else if(padmode == PAD_ADDI_SUBJ)
		{
		  i = ((hmm_results->data[h].stop - jpad) >= 1)    ? (hmm_results->data[h].stop - jpad) : 1;
		  if((h+1) < nhits)
		    next_j = ((hmm_results->data[h+1].start + ipad) <= j0)   ? (hmm_results->data[h+1].start + ipad) : j0;
		  else
		    next_j = -1;
		}
	      /*printf("\tsucked in subseq: hit %d new_i: %d j (still): %d\n", h, i, j);*/
	    }
	}
      /*printf("in RescanFilterSurvivors(): calling actually_search_target: %d %d h: %d nhits: %d\n", i, j, h, nhits);*/
      cm_sc =
	actually_search_target(cm, dsq, i, j, cm_cutoff, cp9_cutoff,
			       results, /* keep results                                 */
			       FALSE,   /* don't filter, we already have                */
			       FALSE,   /* we're not building a histogram for CM stats  */
			       FALSE,   /* we're not building a histogram for CP9 stats */
			       NULL);   /* filter fraction N/A                          */
      flen += (j - i + 1);
      if(cm_sc > best_cm_sc) best_cm_sc = cm_sc;
    }

  //if(flen == 0) ffrac = 100.;
  //else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
  if(ret_flen != NULL) *ret_flen = flen;
  return best_cm_sc;
}

/* Function: CP9ScanPosterior()
 * based on Ian Holmes' hmmer/src/postprob.c::P7EmitterPosterior()
 *
 * Purpose:  Combines Forward and Backward scanning matrices into a posterior
 *           probability matrix. 
 *
 *           The main difference between this function and CP9Posterior()
 *           in hmmband.c is that this func takes matrices from CP9ForwardBackwardScan()
 *           in which parses are allowed to start and end in any residue.
 *           In CP9Posterior(), the matrices are calc'ed in CP9Forward()
 *           and CP9Backward() which force all parses considered to start at posn
 *           1 and end at L. This means here we have to calculate probability
 *           that each residue from 1 to L is contained in any parse prior
 *           to determining the posterior probability it was emitted from
 *           each state.
 * 
 *           For emitters (match and inserts) the 
 *           entries in row i of this matrix are the logs of the posterior 
 *           probabilities of each state emitting symbol i of the sequence. 
 *           For non-emitters the entries in row i of this matrix are the 
 *           logs of the posterior probabilities of each state being 'visited' 
 *           when the last emitted residue in the parse was symbol i of the
 *           sequence. 
 *           The last point distinguishes this function from P7EmitterPosterior() 
 *           which set all posterior values for for non-emitting states to -INFTY.
 *           The caller must allocate space for the matrix, although the
 *           backward matrix can be used instead (overwriting it will not
 *           compromise the algorithm).
 *           
 * Args:     dsq      - sequence in digitized form
 *           i0       - start of target subsequence
 *           j0       - end of target subsequence 
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           mx       - pre-allocated dynamic programming matrix
 *           
 * Return:   void
 */
void
CP9ScanPosterior(char *dsq, int i0, int j0,
		     CP9_t *hmm,
		     CP9_dpmatrix_t *fmx,
		     CP9_dpmatrix_t *bmx,
		     CP9_dpmatrix_t *mx)
{
  int i;
  int ip;
  int k;
  int fb_sum; /* tmp value, the probability that the current residue (i) was
	       * visited in any parse */
  /*printf("\n\nin CP9ScanPosterior() i0: %d, j0: %d\n", i0, j0);*/
  fb_sum = -INFTY;
  for (i = i0-1; i <= j0; i++) 
    {
      ip = i-i0+1; /* ip is relative position in seq, 0..L */
      /*printf("bmx->mmx[i:%d][0]: %d\n", i, bmx->mmx[ip][0]); */
      fb_sum = ILogsum(fb_sum, (bmx->mmx[ip][0])); 
    }
  /*printf("fb_sc: %f\n", Scorify(fb_sum));*/
    /* fb_sum is the probability of all parses */

    /*for(k = 1; k <= hmm->M; k++)*/
  /*{*/
      /*fbsum_ = ILogsum(fmx->mmx[0][k] + bmx->mmx[0][k]))*/; /* residue 0 can't be emitted
								    * but we can start in BEGIN,
								    * before any residues */
      /*fb_sum = ILogsum(fb_sum, (fmx->imx[0][k] + bmx->imx[0][k]))*/; /* these will be all -INFTY */
  /*}*/      

  /* note boundary conditions, case by case by case... */
  mx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - fb_sum;
  mx->imx[0][0] = -INFTY; /*need seq to get here*/
  mx->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  for (k = 1; k <= hmm->M; k++) 
    {
      mx->mmx[0][k] = -INFTY; /*need seq to get here*/
      mx->imx[0][k] = -INFTY; /*need seq to get here*/
      mx->dmx[0][k] = fmx->dmx[0][k] + bmx->dmx[0][k] - fb_sum;
    }
      
  for (i = i0; i <= j0; i++)
    {
      ip = i-i0+1; /* ip is relative position in seq, 0..L */
      /*fb_sum = -INFTY;*/ /* this will be probability of seeing residue i in any parse */
      /*for (k = 0; k <= hmm->M; k++) 
	{
	fb_sum = ILogsum(fb_sum, (fmx->mmx[i][k] + bmx->mmx[i][k] - hmm->msc[dsq[i]][k]));*/
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
      /*fb_sum = ILogsum(fb_sum, (fmx->imx[i][k] + bmx->imx[i][k] - hmm->isc[dsq[i]][k]));*/
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
      /*}*/
      mx->mmx[ip][0] = -INFTY; /*M_0 does not emit*/
      mx->imx[ip][0] = fmx->imx[ip][0] + bmx->imx[ip][0] - hmm->isc[(int) dsq[i]][0] - fb_sum;
      /*hmm->isc[dsq[i]][0] will have been counted in both fmx->imx and bmx->imx*/
      mx->dmx[ip][0] = -INFTY; /*D_0 does not exist*/
      for (k = 1; k <= hmm->M; k++) 
	{
	  mx->mmx[ip][k] = fmx->mmx[ip][k] + bmx->mmx[ip][k] - hmm->msc[(int) dsq[i]][k] - fb_sum;
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  mx->imx[ip][k] = fmx->imx[ip][k] + bmx->imx[ip][k] - hmm->isc[(int) dsq[i]][k] - fb_sum;
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
	  mx->dmx[ip][k] = fmx->dmx[ip][k] + bmx->dmx[ip][k] - fb_sum;
	}	  
    }
  /*
  float temp_sc;
  for (i = i0; i <= j0; i++)
    {
      ip = i-i0+1; 
      for(k = 0; k <= hmm->M; k++)
	{
	  temp_sc = Score2Prob(mx->mmx[ip][k], 1.);
	  if(temp_sc > .0001)
	    printf("mx->mmx[%3d][%3d]: %9d | %8f\n", i, k, mx->mmx[ip][k], temp_sc);
	  temp_sc = Score2Prob(mx->imx[ip][k], 1.);
	  if(temp_sc > .0001)
	    printf("mx->imx[%3d][%3d]: %9d | %8f\n", i, k, mx->imx[ip][k], temp_sc);
	  temp_sc = Score2Prob(mx->dmx[ip][k], 1.);
	  if(temp_sc > .0001)
	    printf("mx->dmx[%3d][%3d]: %9d | %8f\n", i, k, mx->dmx[ip][k], temp_sc);
	}
    }
  */
}

/* Function: CP9ForwardScanDemands()
 * Date:     EPN, Fri Jun 15 10:08:40 2007
 *
 * Purpose:  Determine the number of calculations for 
 *           a CP9 Forward scanner and return it.
 *
 * Args:     cp9    - the CP9 HMM
 *           L      - length of sequence.
 * 
 * Returns: (float) the total number of DP calculations
 */
float
CP9ForwardScanDemands(CP9_t *cp9, int L)
{
  float dpcalcs;	/* # of inner loops executed for non-bif calculations */

  dpcalcs = L * cp9->M * 3.; /* 3 states per node, M nodes in the model */
  return dpcalcs;
}

/*
 * Function: FindCP9FilterThreshold()
 * Incept:   EPN, Wed May  2 10:00:45 2007
 *
 * Purpose:  Sample sequences from a CM and determine the CP9 HMM E-value
 *           threshold necessary to recognize a specified fraction of those
 *           hits. Sequences are sampled from the CM until N with a E-value
 *           better than cm_ecutoff are sampled (those with worse E-values
 *           are rejected). CP9 scans are carried out in either local or
 *           glocal mode depending on hmm_gum_mode. CM is configured in 
 *           local/glocal and sampled seqs are scored in CYK/inside depending
 *           on fthr_mode (4 possibilities). E-values are determined using
 *           lambda from cm->stats, and a recalc'ed mu using database size
 *           of 'db_size'.
 *
 *           If do_fastfil and fthr_mode is local or glocal CYK, parsetree 
 *           scores of emitted sequences are assumed to an optimal CYK scores 
 *           (not nec true). This saves a lot of time b/c no need to
 *           scan emitted seqs, but it's statistically wrong. 
 *           If do_fastfil and fthr_mode is local or glocal inside, 
 *           contract is violated and we Die. Current strategy *outside*
 *           of this function is to copy HMM filtering thresholds from
 *           CYK for Inside cases.
 *
 *           If !do_fastfil this function takes much longer b/c emitted 
 *           parsetree score is not necessarily (a) optimal nor (b) highest 
 *           score returned from a CYK scan (a subseq of the full seq could
 *           score higher even if parsetree was optimal). For Inside, Inside
 *           scan will always be higher than parstree score. This means we have
 *           to scan each emitted seq with CYK or Inside.
 *
 * Args:
 *           cm           - the CM
 *           cmstats      - CM stats object we'll get Gumbel stats from
 *           r            - source of randomness for EmitParsetree()
 *           Fmin         - minimum target fraction of CM hits to detect with CP9
 *           Smin         - minimum filter survival fraction, for lower fractions, we'll
 *                          increase F. 
 *           Starget      - target filter survival fraction, for lower fractions that 
 *                          satisfy Fmin, we'll increase F until Starget is reached
 *           Spad         - fraction of (sc(S) - sc(max(Starget, Smin))) to subtract from sc(S)
 *                          in case 2. 0.0 = fast, 1.0 = safe.
 *           N            - number of sequences to sample from CM better than cm_minsc
 *           use_cm_cutoff- TRUE to only accept CM parses w/E-values above cm_ecutoff
 *           cm_ecutoff   - minimum CM E-value to accept 
 *           db_size      - DB size (nt) to use w/cm_ecutoff to calc CM bit score cutoff 
 *           emit_mode    - CM_LC or CM_GC, CM mode to emit with
 *           fthr_mode    - gives CM search strategy to use, and Gumbel to use
 *           hmm_gum_mode - CP9_L to search with  local HMM (we're filling a fthr->l_eval)
 *                          CP9_G to search with glocal HMM (we're filling a fthr->g_eval)
 *           do_fastfil   - TRUE to use fast method: assume parsetree score
 *                          is optimal CYK score
 *           do_Fstep     - TRUE to step from F towards 1.0 while S < Starget in case 2.
 *           my_rank      - MPI rank, 0 if serial
 *           nproc        - number of processors in MPI rank, 1 if serial
 *           do_mpi       - TRUE if we're doing MPI, FALSE if not
 *           histfile     - root string for histogram files, we'll make nexps * 4 of them
 *           ret_F        - the fraction of observed CM hits we've scored with the HMM better
 *                          than return value
 * 
 * Returns: HMM E-value cutoff above which the HMM scores (ret_F * N) CM 
 *          hits with CM E-values better than cm_ecutoff 
 * 
 */
float FindCP9FilterThreshold(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
			     float Fmin, float Smin, float Starget, float Spad, int N, 
			     int use_cm_cutoff, float cm_ecutoff, int db_size, 
			     int emit_mode, int fthr_mode, int hmm_gum_mode, 
			     int do_fastfil, int do_Fstep, int my_rank, int nproc, int do_mpi, 
			     char *histfile, float *ret_F)
{

  /* Contract checks */
  if (!(cm->flags & CM_CP9) || cm->cp9 == NULL) 
    esl_fatal("ERROR in FindCP9FilterThreshold() CP9 does not exist\n");
  if (Fmin < 0. || Fmin > 1.)  
    esl_fatal("ERROR in FindCP9FilterThreshold() Fmin is %f, should be [0.0..1.0]\n", Fmin);
  if (Smin < 0. || Smin > 1.)  
    esl_fatal("ERROR in FindCP9FilterThreshold() Smin is %f, should be [0.0..1.0]\n", Smin);
  if (Starget < 0. || Starget > 1.)  
    esl_fatal("ERROR in FindCP9FilterThreshold() Starget is %f, should be [0.0..1.0]\n", Starget);
  if((fthr_mode != CM_LI) && (fthr_mode != CM_GI) && (fthr_mode != CM_LC) && (fthr_mode != CM_GC))
    esl_fatal("ERROR in FindCP9FilterThreshold() fthr_mode not CM_LI, CM_GI, CM_LC, or CM_GC\n");
  if(hmm_gum_mode != CP9_L && hmm_gum_mode != CP9_G)
    esl_fatal("ERROR in FindCP9FilterThreshold() hmm_gum_mode not CP9_L or CP9_G\n");
  if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI))
    esl_fatal("ERROR in FindCP9FilterThreshold() do_fastfil TRUE, but fthr_mode CM_GI or CM_LI\n");
  if(my_rank > 0 && !do_mpi)
    esl_fatal("ERROR in FindCP9FilterThreshold() my_rank is not 0, but do_mpi is FALSE\n");
  if(emit_mode != CM_GC && emit_mode != CM_LC)
    esl_fatal("ERROR in FindCP9FilterThreshold() emit_mode not CM_LC or CM_GC\n");
  if(emit_mode == CM_LC && (fthr_mode == CM_GC || fthr_mode == CM_GI))
    esl_fatal("ERROR in FindCP9FilterThreshold() emit_mode CM_LC but score mode CM_GC or CM_GI.\n");
  if(Spad < 0 || Spad > 1.0)
    esl_fatal("ERROR in FindCP9FilterThreshold() Spad %f not between 0.0 and 1.0\n");

#if defined(USE_MPI)  && defined(NOTDEFINED)
  /* If a worker in MPI mode, we go to worker function mpi_worker_cm_and_cp9_search */
  if(my_rank > 0) 
    {
      /* Configure the CM based on the fthr mode COULD BE DIFFERENT FROM MASTER! */
      ConfigForGumbelMode(cm, fthr_mode);
      /* Configure the HMM based on the hmm_gum_mode */
      if(hmm_gum_mode == CP9_L)
	CPlan9SWConfig(cm->cp9, 0.5, 0.5);
      else /* hmm_gum_mode == CP9_G (it's in the contract) */
	CPlan9GlobalConfig(cm->cp9);
      CP9Logoddsify(cm->cp9);

      //mpi_worker_cm_and_cp9_search(cm, do_fastfil, my_rank);
      mpi_worker_cm_and_cp9_search_maxsc(cm, do_fastfil, do_minmax, my_rank);

      *ret_F = 0.0; /* this return value is irrelevant */
      return 0.0;   /* this return value is irrelevant */
    }
#endif 
  int            status;         /* Easel status */
  CM_t          *cm_for_scoring; /* used to score parsetrees, nec b/c  *
				  * emitting mode may != scoring mode   */
  Parsetree_t   *tr = NULL;      /* parsetree                           */
  char          *dsq = NULL;     /* digitized sequence                  */
  char          *seq = NULL;     /* alphabetic sequence                 */
  float         *tr_sc;          /* scores of all parsetrees sampled    */
  float         *cm_sc_p;        /* CM scores of samples above thresh   */
  float         *hmm_sc_p;       /* HMM scores of samples above thresh  */
  float         *hmm_eval_p;     /* HMM E-values of samples above thresh*/  
  int            i  = 0;         /* counter over samples                */
  int            ip = 0;         /* counter over samples above thresh   */
  int            imax = 500 * N; /* max num samples                     */
  int            p;              /* counter over partitions             */
  int            L;              /* length of a sample                  */
  float          F;              /* fraction of hits found by HMM >= Fmin*/
  float          E;              /* HMM CP9 E-value cutoff to return    */
  float          Emin;           /* E-value that corresponds to Smin */
  float          Etarget;        /* E-value that corresponds to Starget */
  double        *cm_mu;          /* mu for each partition's CM Gumbel   */
  double        *hmm_mu;         /* mu for each partition's HMM Gumbel  */
  float         *cm_minbitsc = NULL; /* minimum CM bit score to allow to pass for each partition */
  float         *cm_maxbitsc = NULL; /* maximum CM bit score to allow to pass for each partition */
  double         tmp_K;          /* used for recalc'ing Gumbel stats for DB size */
  int            was_hmmonly;    /* flag for if HMM only search was set */
  int            nalloc;         /* for cm_sc, hmm_sc, hmm_eval, num alloc'ed */
  int            chunksize;      /* allocation chunk size               */
  float         *scores;         /* CM and HMM score returned from worker*/
  void          *tmp;            /* temp void pointer for ESL_RALLOC() */
  int            clen;           /* consensus length of CM             */
  float          avg_hit_len;    /* crude estimate of average hit len  */
  int            Fidx;           /* index within hmm_eval              */
  float          S;              /* predicted survival fraction        */
  int            init_flag;      /* used for finding F                 */ 
  int            tr_npassed = 0; /* number of traces that passed cutoffs
				  * we had to search before getting N  */
  int            passed_flag = FALSE;
  int            tr_passed_flag = FALSE;
  float          cm_sc;
  float          hmm_sc;
  float          orig_tau;
  float          hb_sc= 0.; 
  float          S_sc = 0.;
  float          Starget_sc = 0.;

#if defined(USE_MPI)  && defined(NOTDEFINED)
  int            have_work;      /* MPI: do we still have work to send out?*/
  int            nproc_working;  /* MPI: number procs currently doing work */
  int            wi;             /* MPI: worker index                      */
  char         **dsqlist = NULL; /* MPI: queue of digitized seqs being worked on, 1..nproc-1 */
  int           *plist = NULL;   /* MPI: queue of partition indices of seqs being worked on, 1..nproc-1 */
  Parsetree_t  **trlist = NULL;  /* MPI: queue of traces of seqs being worked on, 1..nproc-1 */
  MPI_Status     mstatus;        /* useful info from MPI Gods         */
#endif

  /* TEMPORARY! */
  do_mpi = FALSE;
  printf("in FindCP9FilterThreshold fthr_mode: %d hmm_gum_mode: %d emit_mode: %d\n", fthr_mode, 
	 hmm_gum_mode, emit_mode);

  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  cm->search_opts &= ~CM_SEARCH_HMMONLY;

  chunksize = 5 * N;
  nalloc    = chunksize;
  ESL_ALLOC(tr_sc,       sizeof(float) * nalloc);
  ESL_ALLOC(hmm_eval_p,  sizeof(float) * N);
  ESL_ALLOC(hmm_sc_p,    sizeof(float) * N);
  ESL_ALLOC(cm_sc_p,     sizeof(float) * N);
  ESL_ALLOC(cm_minbitsc, sizeof(float) * cmstats->np);
  ESL_ALLOC(cm_mu,       sizeof(double)* cmstats->np);
  ESL_ALLOC(hmm_mu,      sizeof(double)* cmstats->np);
  ESL_ALLOC(scores,      sizeof(float) * 2);
  
#if defined(USE_MPI) && defined(NOTDEFINED)
  if(do_mpi)
    {
      ESL_ALLOC(dsqlist,     sizeof(char *)        * nproc);
      ESL_ALLOC(plist,       sizeof(int)           * nproc);
      ESL_ALLOC(trlist,      sizeof(Parsetree_t *) * nproc);
    }
#endif

  if(use_cm_cutoff) printf("CM E cutoff: %f\n", cm_ecutoff);
  else              printf("Not using CM cutoff\n");

  /* Configure the CM based on the emit mode COULD DIFFERENT FROM WORKERS! */
  ConfigForGumbelMode(cm, emit_mode);
  /* Copy the CM into cm_for_scoring, and reconfigure it if nec.,
   * We do this, so we can exponentiate or change emission modes of
   * the original CM */
  cm_for_scoring = DuplicateCM(cm); 
  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
    ConfigForGumbelMode(cm_for_scoring, fthr_mode);

  /* Configure the HMM based on the hmm_gum_mode */
  if(hmm_gum_mode == CP9_L)
    CPlan9SWConfig(cm->cp9, 0.5, 0.5);
  else /* hmm_gum_mode == CP9_G (it's in the contract) */
    CPlan9GlobalConfig(cm->cp9);
  CP9Logoddsify(cm->cp9);

  /* Determine bit cutoff for each partition, calc'ed from cm_ecutoff */
  for (p = 0; p < cmstats->np; p++)
    {
      /* First determine mu based on db_size */
      tmp_K      = exp(cmstats->gumAA[hmm_gum_mode][p]->mu * cmstats->gumAA[hmm_gum_mode][p]->lambda) / 
	cmstats->gumAA[hmm_gum_mode][p]->L;
      hmm_mu[p]  = log(tmp_K * ((double) db_size)) / cmstats->gumAA[hmm_gum_mode][p]->lambda;
      tmp_K      = exp(cmstats->gumAA[fthr_mode][p]->mu * cmstats->gumAA[fthr_mode][p]->lambda) / 
	cmstats->gumAA[fthr_mode][p]->L;
      cm_mu[p]   = log(tmp_K  * ((double) db_size)) / cmstats->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cm_minbitsc[p] = cm_mu[p] - (log(cm_ecutoff) / cmstats->gumAA[fthr_mode][p]->lambda);
      if(use_cm_cutoff)
	printf("E: %f p: %d %d--%d bit score: %f\n", cm_ecutoff, p, 
	       cmstats->ps[p], cmstats->pe[p], cm_minbitsc[p]);
    }
  
  /* Do the work, emit parsetrees and collect the scores 
   * either in serial or MPI depending on do_mpi flag.
   */
  /*********************SERIAL BLOCK*************************************/
  int nleft = 0; /* number of seqs with scores < min CM score */
  int tr_np, tr_na, s1_np, s1_na, s2_np, s2_na, s3_np, s3_na;
  tr_np = tr_na = s1_np = s1_na = s2_np = s2_na = s3_np = s3_na = 0;
  printf("06.11.07 Min np: %5d Smin: %12f Starget: %f Spad: %.3f do_Fstep: %d\n", N, Smin, Starget, Spad, do_Fstep);
  if(!(do_mpi))
    {
      orig_tau = cm_for_scoring->tau;

      /* Serial strategy: 
       * Emit sequences one at a time and search them with the CM and HMM,
       * keeping track of scores. If do_fastfil we don't have to search 
       * with the CM we just keep track of the parsetree score. 
       *
       * If do_fastfil && emit_mode is different than scoring mode we 
       * score with 'cm_for_scoring', instead of with 'cm'.
       */

      SQINFO            sqinfo;     /* info about sequence (name/len)         */
      while(ip < N) /* while number seqs passed CM score threshold (ip) < N */
	{
	  EmitParsetree(cm, r, &tr, &seq, &dsq, &L);
	  while(L == 0) { FreeParsetree(tr); free(seq); free(dsq); EmitParsetree(cm, r, &tr, &seq, &dsq, &L); }

	  tr_na++;
	  passed_flag = FALSE;

	  p = cmstats->gc2p[(get_gc_comp(seq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */
	  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
	    tr_sc[i] = ParsetreeScore_Global2Local(cm_for_scoring, tr, dsq, FALSE);
	  else
	    tr_sc[i] = ParsetreeScore(cm_for_scoring, tr, dsq, FALSE); 
	  /*sprintf(sqinfo.name, "seq%d", i+1);
	  sqinfo.len   = L;
	  sqinfo.flags = SQINFO_NAME | SQINFO_LEN;
	  WriteSeq(stdout, SQFILE_FASTA, seq, &sqinfo);
	  ParsetreeDump(stdout, tr, cm, dsq);
	  printf("%d Parsetree Score: %f\n\n", i, tr_sc[i]);*/
	  //FreeParsetree(tr);
	  //free(seq);

	  /* If do_minmax, check if the parsetree score less than maximum allowed */
	  if(tr_sc[i] > cm_minbitsc[p]) /* we know we've passed */
	    {
	      tr_np++;
	      passed_flag = TRUE;
	      cm_sc_p[ip] = tr_sc[i];
	      //printf("TR P (P: %5d L: %5d)\n", ip, nleft);
	    }
	  else
	    {
	      /* we're not sure if our optimal score exceeds cm_minbitsc */
	      /* STAGE 1 */
	      
	      /* For speed first see if a strict (high tau) HMM banded search finds a 
	       * conditional optimal parse with score > min score */
	      s1_na++;
	      cm_for_scoring->search_opts |= CM_SEARCH_HBANDED;
	      //cm_for_scoring->tau = 0.1;
	      cm_for_scoring->tau = 0.01;
	      //cm_for_scoring->tau = 0.001;
	      
	      hb_sc = actually_search_target(cm_for_scoring, dsq, 1, L,
					     0.,    /* cutoff is 0 bits (actually we'll find highest
						     * negative score if it's < 0.0) */
					     0.,    /* CP9 cutoff is 0 bits */
					     NULL,  /* don't keep results */
					     FALSE, /* don't filter with a CP9 HMM */
					     FALSE, /* we're not calcing CM  stats */
					     FALSE, /* we're not calcing CP9 stats */
					     NULL); /* filter fraction N/A */
	      //if(!do_fastfil) printf("%4d %5d %d T: %10.4f BC: %10.4f ", ip, i, passed_flag, tr_sc[i], hb_sc);
	      if(hb_sc > cm_minbitsc[p])
		{
		  s1_np++;
		  passed_flag = TRUE;
		  cm_sc_p[ip] = hb_sc;
		  //printf("S1 P (P: %5d L: %5d)\n", ip, nleft);
		}
	      else /* Stage 2, search with another, less strict (lower tau)  HMM banded parse */
		{
		  s2_na++;
		  cm_for_scoring->search_opts |= CM_SEARCH_HMMSCANBANDS;
		  cm_for_scoring->tau = 1e-15;
		  cm_sc = actually_search_target(cm_for_scoring, dsq, 1, L,
						 0.,    /* cutoff is 0 bits (actually we'll find highest
							 * negative score if it's < 0.0) */
						 0.,    /* CP9 cutoff is 0 bits */
						 NULL,  /* don't keep results */
						 FALSE, /* don't filter with a CP9 HMM */
						 FALSE, /* we're not calcing CM  stats */
						 FALSE, /* we're not calcing CP9 stats */
						 NULL); /* filter fraction N/A */
		  if(cm_sc > cm_minbitsc[p])
		    {
		      s2_np++;
		      passed_flag = TRUE;
		      cm_sc_p[ip] = cm_sc;
		      //printf("S2 P (P: %5d L: %5d)\n", ip, nleft);
		    }
		  else /* search for the optimal parse */
		    {
		      s3_na++;
		      /* Stage 3 do QDB CYK */
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HBANDED;
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HMMSCANBANDS;
		      cm_sc = actually_search_target(cm_for_scoring, dsq, 1, L,
						     0.,    /* cutoff is 0 bits (actually we'll find highest
							     * negative score if it's < 0.0) */
						     0.,    /* CP9 cutoff is 0 bits */
						     NULL,  /* don't keep results */
						     FALSE, /* don't filter with a CP9 HMM */
						     FALSE, /* we're not calcing CM  stats */
						     FALSE, /* we're not calcing CP9 stats */
						     NULL); /* filter fraction N/A */
		      if(cm_sc > cm_minbitsc[p])
			{
			  s3_np++;
			  passed_flag = TRUE;
			  cm_sc_p[ip] = cm_sc;
			  //printf("S3 P (P: %5d L: %5d)\n", ip, nleft);
			}
		    }
		}
	    }

	  if(!passed_flag) 
	    {
	      nleft++;
	      //printf("LEFT (P: %5d L: %5d)\n", ip, nleft);
	    }
	  else if(passed_flag)
	    {
	      /* Scan seq with HMM */
	      /* DO NOT CALL actually_search_target b/c that will run Forward then 
	       * Backward to get score of best hit, but we'll be detecting by a
	       * Forward scan (then running Backward only on hits above our threshold).
	       */
	      hmm_sc_p[ip] = CP9Forward(cm, dsq, 1, L, cm->W, 0., 
					NULL,   /* don't return scores of hits */
				       NULL,   /* don't return posns of hits */
				       NULL,   /* don't keep track of hits */
				       TRUE,   /* we're scanning */
				       FALSE,  /* we're not ultimately aligning */
				       FALSE,  /* we're not rescanning */
				       TRUE,   /* be memory efficient */
				       NULL);  /* don't want the DP matrix back */
	      hmm_eval_p[ip] = RJK_ExtremeValueE(hmm_sc_p[ip], hmm_mu[p], cmstats->gumAA[hmm_gum_mode][p]->lambda);
	      ip++; /* increase counter of seqs passing threshold */
	    }
	  free(dsq);
	  
	  /* Check if we need to reallocate */
	  i++;
	  if(i > imax) esl_fatal("ERROR number of attempts exceeded 500 times number of seqs.\n");
	  if (i == nalloc) 
	    {
	      nalloc += chunksize;
	      ESL_RALLOC(tr_sc,    tmp, nalloc * sizeof(float));
	    }
	  if(ip % N == 0)
	    {
	      ESL_RALLOC(hmm_eval_p, tmp, (ip + N) * sizeof(float));
	      ESL_RALLOC(hmm_sc_p,   tmp, (ip + N) * sizeof(float));
	      ESL_RALLOC(cm_sc_p,    tmp, (ip + N) * sizeof(float));
	    }
	}
    }
  N = ip; /* update N based on number of seqs sampled */
  printf("06.11.07 np: %5d nleft: %5d\n", ip, nleft);
  printf("06.11.07 tr_a: %5d tr_p: %5d\n06.11.07 s1_a: %5d s1_p: %5d\n06.11.07 s2_a: %5d s2_p: %5d\n06.11.07 s3_a: %5d s3_p: %5d\n", tr_na, tr_np, s1_na, s1_np, s2_na, s2_np, s3_na, s3_np);
  
  /*************************END OF SERIAL BLOCK****************************/
#if defined(USE_MPI)  && defined(NOTDEFINED)
  /*************************MPI BLOCK****************************/
  if(do_mpi)
    {
      /* MPI Strategy: 
       * Emit seqs one at a time from the CM, if the CM parse tree score is greater
       * than our max E-value (if do_minmax, otherwise there is no maximum),
       * send it to a worker. The worker then tries to quickly determine if
       * the sequence is within the acceptable E-value range by doing HMM
       * banded searches. If the sequence is within the E-value range,
       * it is searched with an HMM. Both the optimal CM parse scores and HMM
       * scores are passed back to the master, whose keeping track of how
       * many seqs have been sampled within the E-value range.
       *
       * If do_fastfil, the worker skips the CM search and returns 0. bits
       * as CM score, which master replaces with parsetree score. (This is 
       * an old strategy I'm probably about to deprecate (06.07.07))
       *
       * Sean's design pattern for data parallelization in a master/worker model:
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
	  if(ip < N) /* if number seqs passed CM score threshold (ip) < N */
	    {
	      tr_passed_flag = FALSE;
	      while(!tr_passed_flag)
		{
		  EmitParsetree(cm, r, &tr, &seq, &dsq, &L);
		  while(L == 0) { FreeParsetree(tr); free(seq); free(dsq); EmitParsetree(cm, r, &tr, &seq, &dsq, &L); }

		  p = cmstats->gc2p[(get_gc_comp(seq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */
		  free(seq);
		  /*ParsetreeDump(stdout, tr, cm, dsq);
		    printf("%d Parsetree Score: %f\n\n", (nattempts), ParsetreeScore(cm, tr, dsq, FALSE)); */
		  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
		    tr_sc[i] = ParsetreeScore_Global2Local(cm_for_scoring, trlist[wi], dsqlist[wi], FALSE);
		  else
		    tr_sc[i] = ParsetreeScore(cm_for_scoring, tr, dsq, FALSE); 
		  FreeParsetree(tr);
		  /* If do_minmax, check if the parsetree score less than maximum allowed */
		  if((!do_minmax || (do_minmax && (tr_sc[i] <= cm_maxbitsc[p]))))
		    tr_passed_flag = TRUE;
		}		    
	      if(!do_fastfil)
		FreeParsetree(tr);
	    }
	  else have_work = FALSE;
	  /* If we have work but no free workers, or we have no work but workers
	   * Are still working, then wait for a result to return from any worker.
	   */
	  if ( (have_work && nproc_working == nproc-1) || (! have_work && nproc_working > 0))
	    {
	      MPI_Recv(scores,  2, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mstatus);
	      cm_sc    = scores[0];
	      hmm_sc   = scores[1];
	      wi = mstatus.MPI_SOURCE;

	      if(do_fastfil)
		{
		  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
		    cm_sc = ParsetreeScore_Global2Local(cm_for_scoring, trlist[wi], dsqlist[wi], FALSE);
		  else
		    cm_sc = ParsetreeScore(cm_for_scoring, trlist[wi], dsqlist[wi], FALSE); 
		  FreeParsetree(trlist[wi]);
		  trlist[wi] = NULL;
		}
	      if(ip < N && 
		 ( do_minmax && (cm_sc >= cm_minbitsc[plist[wi]] && cm_sc <= cm_maxbitsc[plist[wi]])) ||
		 (!do_minmax && (cm_sc >= cm_minbitsc[plist[wi]])))
		{
		  cm_sc_p[ip]    = cm_sc;
		  hmm_sc_p[ip]   = hmm_sc;
		  hmm_eval_p[ip] = RJK_ExtremeValueE(hmm_sc, hmm_mu[plist[wi]], cmstats->gumAA[hmm_gum_mode][plist[wi]]->lambda);
		  ip++; /* increase counter of seqs passing threshold */
		}
	      nproc_working--;
	      free(dsqlist[wi]);
	      dsqlist[wi] = NULL;
	      
	      /* Check if we need to reallocate */
	      i++;
	      if(i > imax) esl_fatal("ERROR number of attempts exceeded 500 times number of seqs.\n");
	      if (i == nalloc) 
		{
		  nalloc += chunksize;
		  ESL_RALLOC(tr_sc,    tmp, nalloc * sizeof(float));
		}
	    }
	  /* If we have work, assign it to a free worker;
	   * else, terminate the free worker.
	   */
	  if (have_work) 
	    {
	      dsq_maxsc_MPISend(dsq, L, cm_maxbitsc[p], wi);
	      dsqlist[wi] = dsq;
	      plist[wi]   = p;
	      if(do_fastfil)
		trlist[wi] = tr;
	      wi++;
	      nproc_working++;
	    }
	  else 
	    dsq_maxsc_MPISend(NULL, -1, -1, wi);	
	}
    }
  /*********************END OF MPI BLOCK*****************************/
#endif /* if defined(USE_MPI) */
  /* Sort the HMM E-values with quicksort */
  esl_vec_FSortIncreasing(hmm_eval_p, N);
  esl_vec_FSortDecreasing(hmm_sc_p, N); /* TEMPORARY FOR PRINTING */

  /* Determine E-value to return based on predicted filter survival fraction */
  clen = (2*CMCountStatetype(cm, MATP_MP) + 
	  CMCountStatetype(cm, MATL_ML) + 
	  CMCountStatetype(cm, MATR_MR));
  avg_hit_len = clen; /* currently, we don't correct for local hits probably being shorter */

  /* Strategy 1, enforce minimum F, but look for highest that gives survival fraction of at least
   * Starget */
  init_flag  = TRUE;
  S          = IMPOSSIBLE;
  Emin    = ((double) db_size * Smin)    / ((2. * cm_for_scoring->W) - avg_hit_len);
  Etarget = ((double) db_size * Starget) / ((2. * cm_for_scoring->W) - avg_hit_len);
  printf("Smin: %12f\nStarget: %f\n", Smin, Starget);
  printf("Emin: %12f\nEtarget: %f\n", Emin, Etarget);

  Fidx  = (int) (Fmin * (float) N) - 1; /* off by one, ex: 95% cutoff in 100 len array is index 94 */
  E     = hmm_eval_p[Fidx];

  /* Are we case 1 or case 2? 
   * Case 1: our E is greater than our Etarget, so we cannot satisfy Starget AND capture F fraction of
   *         the true CM hits, it's more impt we capture F fraction, so return S > Starget.
   * Case 2: our E is less than our Etarget so we can satisfy both criteria, we have to choose an S
   *         now. One option is If (do_Fstep) increase F toward 1.0 while E < Etarget. Then add
   *         Spad fraction of (bitscore(S) - bitscore(max(Starget, Smin))) to bitscore(S), and calculate
   *         the new S. If Spad == 0., speed trumps filter safety, if S == 1.0, we want the filter
   *         as safe as possible, we're returning S = Starget. 
   */
  if(E > Etarget) 
    {
      F = ((float) (Fidx + 1.)) / (float) N;
      S = (E * ((2. * cm_for_scoring->W) - avg_hit_len)) / ((double) db_size); 
      if(E < Emin)  /* no need for E < Emin */
	  printf("Case 1A rare case: init Emin %12f > E %f > Etarget %f F: %f S: %.12f Spad: %.3f\n", Emin, E, Etarget, F, S, Spad);
      else
	  printf("Case 1B bad  case: init E %f > Etarget %f && E > Emin %12f F: %f S: %.12f Spad: %.3f\n", E, Etarget, Emin, F, S, Spad);
    }
  else
    {
      /* If(do_Fstep): increase F until we're about to go over Etarget (which gives our target 
	 survival fraction Starget). */
      if(do_Fstep)
	while(((Fidx+1) < N) && hmm_eval_p[(Fidx+1)] < Etarget)
	  Fidx++;
      F = ((float) (Fidx + 1.)) / (float) N;

      /* Subtract Spad * (score(S) - score(max(Starget, Smin))) from score(S) and recalculate E&S */
      /* Use partition that covers GC = 50 */
      p = cm->stats->gc2p[50];
      E = hmm_eval_p[Fidx];
      printf("Before Spad E: %f\n", E);
      S_sc =    (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
		 (log(E) / cm->stats->gumAA[hmm_gum_mode][p]->lambda)); 
      if(Starget > Smin)
	Starget_sc = (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
		   (log(Etarget) / cm->stats->gumAA[hmm_gum_mode][p]->lambda));
      else
	Starget_sc = (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
		      (log(Emin) / cm->stats->gumAA[hmm_gum_mode][p]->lambda));

      printf("S_sc 0: %.3f - %.3f * %.3f = ", S_sc, Spad, S_sc - Starget_sc);
      S_sc -= Spad * (S_sc - Starget_sc); /* Spad may be 0. */
      printf(" S1: %.3f\n", S_sc);
      /* now recalculate what E and S should be based on S_sc */
      E = RJK_ExtremeValueE(S_sc, cm->stats->gumAA[hmm_gum_mode][p]->mu, 
			    cmstats->gumAA[hmm_gum_mode][p]->lambda);
      S = (E * ((2. * cm_for_scoring->W) - avg_hit_len)) / ((double) db_size); 

      if(fabs(E - Emin) < 0.01 || E < Emin)  /* no need for E < Emin */
	printf("Case 2A: best case Emin %12f > E %f < Etarget %f F: %f S: %.12f Spad: %.3f\n", Emin, E, Etarget, F, S, Spad);
       else
	printf("Case 2B: good case Emin %12f < E %f < Etarget %f F: %f S: %.12f\n Spad: %.3f", Emin, E, Etarget, F, S, Spad);
      
     }
  /* Make sure our E is less than the DB size and greater than Emin */
  if(E > ((float) db_size)) /* E-val > db_size is useless */
    {
      printf("Case 3 : worst case E > db_size\n", E, db_size);
      E = (float) db_size;
    }  
  if(E < Emin) E = Emin; /* No reason to have an E below Emin */
  
  /* Informative, temporary print statements */
  for (i = ((int) (Fmin  * (float) N) -1); i < N; i++)
    printf("%d i: %4d hmm sc: %10.4f hmm E: %10.4f\n", hmm_gum_mode, i, hmm_sc_p[i], hmm_eval_p[i]);
  printf("\nSummary: %d %d %d %d %f %f\n", fthr_mode, hmm_gum_mode, do_fastfil, emit_mode,
	 (hmm_sc_p[Fidx]), (hmm_eval_p[Fidx]));
  printf("05.21.07 %d %d %f %f\n", fthr_mode, hmm_gum_mode, hmm_sc_p[Fidx], E);

  /* Reset CM_SEARCH_HMMONLY search option as it was when function was entered */
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;

  /* Clean up and exit */
  free(tr_sc);
  free(hmm_eval_p);
  free(hmm_sc_p);
  free(cm_sc_p);
  free(hmm_mu);
  free(cm_mu);
  free(cm_minbitsc);
  free(cm_maxbitsc);
  free(scores);
  cm_for_scoring->tau = orig_tau;
  FreeCM(cm_for_scoring);
  /* Return threshold */
  *ret_F = F;
  printf("F: %f\n", *ret_F);
  return E;

 ERROR:
  esl_fatal("Reached ERROR in FindCP9FilterThreshold()\n");
  return 0.;
}

/*
 * Function: FindExpFactor()
 * Incept:   EPN, Mon Jun  4 08:14:49 2007
 *
 * Purpose:  Find the exponentiation factor 'X', that when a CM's emission
 *           and transition probabilities are exponentiated by it a specified
 *           fraction range of sample parsetrees have scores greater than
 *           a specified minimum E-value. Often used prior to calling 
 *           FindCP9FilterThreshold().
 *
 * Args:
 *           cm           - the CM
 *           cmstats      - CM stats object we'll get Gumbel stats from
 *           r            - source of randomness for EmitParsetree()
 *           use_cm_cutoff- TRUE to only accept CM parses w/E-values above cm_ecutoff
 *           cm_ecutoff   - minimum CM E-value to accept 
 *           db_size      - DB size (nt) to use w/cm_ecutoff to calc CM bit score cutoff 
 *           emit_mode    - CM_LC or CM_GC, CM mode to emit with
 *           fthr_mode    - gives CM search strategy to use, and Gumbel to use
 *           do_fastfil   - TRUE to use fast method: assume parsetree score
 *                          is optimal CYK score
 *           ntrials      - number of parsetrees to sample 
 *           fp_min       - min fraction of ntrials we require > cm_ecutoff
 *           fp_max       - max fraction of ntrials we require > cm_ecutoff
 *           
 * Returns: Exponentiation factor 'X' for which CM^X emits between (ntrails * fp_min) 
 *          and (ntrials * fp_max) parsetrees with scores better than cm_ecutoff.
 *          Dies if no such 'X' can be found.
 */
float FindExpFactor(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
		    int use_cm_cutoff, float cm_ecutoff, int db_size, 
		    int emit_mode, int fthr_mode, int do_fastfil,
		    int ntrials, float fp_min, float fp_max)
{
  /* Contract checks */
  if((fthr_mode != CM_LI) && (fthr_mode != CM_GI) && (fthr_mode != CM_LC) && (fthr_mode != CM_GC))
    Die("ERROR in FindExpFactor() fthr_mode not CM_LI, CM_GI, CM_LC, or CM_GC\n");
  if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI))
    Die("ERROR in FindExpFactor() do_fastfil TRUE, but fthr_mode CM_GI or CM_LI\n");
  if(emit_mode != CM_GC && emit_mode != CM_LC)
    Die("ERROR in FindExpFactor() emit_mode not CM_LC or CM_GC\n");
  if(emit_mode == CM_LC && (fthr_mode == CM_GC || fthr_mode == CM_GI))
    Die("ERROR in FindExpFactor() emit_mode CM_LC but score mode CM_GC or CM_GI.\n");
  if(fp_min > fp_max)
    Die("ERROR in FindExpFactor() fp_min > fp_max\n");
  if(fp_min < 0. || fp_min > 1.)
    Die("ERROR in FindExpFactor() fp_min not between and 0. and 1.\n");
  if(fp_max < 0. || fp_max > 1.)
    Die("ERROR in FindExpFactor() fp_max not between and 0. and 1.\n");

  int          keep_going = TRUE;
  float        X = 1.0;
  float        init_high = 2.0;
  float        init_low  = 0.5;
  float        high;
  float        low;
  int          X_ctr = 0;
  int          np_min = ntrials * fp_min;
  int          np_max = ntrials * fp_max;
  float        sc;
  int          np = 0;
  float        old_X;
  float       *cm_minbitsc;
  double        *cm_mu;          /* mu for each partition's CM Gumbel   */
  CM_t        *cm_for_scoring;
  Parsetree_t *tr = NULL;      /* parsetree                           */
  char        *dsq = NULL;     /* digitized sequence                  */
  char        *seq = NULL;     /* alphabetic sequence                 */
  int          status;
  int          p;
  float        tmp_K;
  int            i  = 0;         /* counter over samples                */
  int            L;              /* length of a sample                  */

  high = init_high;
  low  = init_low;
  if(!do_fastfil) esl_fatal("ERROR in FindExpFactor(), do_fastfil not TRUE, this is not yet implemented, should it be?\n");

  ESL_ALLOC(cm_minbitsc, sizeof(float) * cmstats->np);
  ESL_ALLOC(cm_mu,       sizeof(double)* cmstats->np);
  /* Configure the CM based on the emit mode COULD DIFFERENT FROM WORKERS! */
  ConfigForGumbelMode(cm, emit_mode);
  /* Copy the CM into cm_for_scoring, and reconfigure it if nec.,
   * We do this, so we can exponentiate or change emission modes of
   * the original CM */
  cm_for_scoring = DuplicateCM(cm); 
  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
    ConfigForGumbelMode(cm_for_scoring, fthr_mode);

  /* Determine bit cutoff for each partition, calc'ed from cm_ecutoff */
  for (p = 0; p < cmstats->np; p++)
    {
      /* First determine mu based on db_size */
      tmp_K      = exp(cmstats->gumAA[fthr_mode][p]->mu * cmstats->gumAA[fthr_mode][p]->lambda) / 
	cmstats->gumAA[fthr_mode][p]->L;
      cm_mu[p]   = log(tmp_K  * ((double) db_size)) / cmstats->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cm_minbitsc[p] = cm_mu[p] - (log(cm_ecutoff) / cmstats->gumAA[fthr_mode][p]->lambda);
      if(use_cm_cutoff)
	printf("FindExpFactor() E: %f p: %d %d--%d bit score: %f\n", cm_ecutoff, p, 
	       cmstats->ps[p], cmstats->pe[p], cm_minbitsc[p]);
    }

  /* Do a binary search for 'X' that gives between np_min and np_max parsetrees with sc better
   * than cm_minsc */
  while(keep_going)
    {
      ExponentiateCM(cm, X);
      old_X = X;
      np = 0;
      for(i = 0; i < ntrials; i++)
	{
	  EmitParsetree(cm, r, &tr, &seq, &dsq, &L);
	  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
	    sc = ParsetreeScore_Global2Local(cm_for_scoring, tr, dsq, FALSE);
	  else
	    sc = ParsetreeScore(cm_for_scoring, tr, dsq, FALSE); 
	  FreeParsetree(tr);
	  free(dsq);
	  p = cmstats->gc2p[(get_gc_comp(seq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */
	  free(seq);
	  if(sc > cm_minbitsc[p]) np++;
	}
      printf("06.03.07 %4d X: %.5f np: %d high: %.5f low: %.5f\n", X_ctr, X, np, high, low);
      X_ctr++;
      if(np > np_min && np < np_max)
	{
	  keep_going = FALSE;
	  printf("06.03.07 success! X: %f\n", X);
	}
      else
	{
	  if      (np < np_min) { low  = X;  X += (high - X)/2.; }
	  else if (np > np_max) { high = X;  X -= (X  - low)/2.;  }
	  if(fabs(high-low) < 0.0001) 
	    if(fabs(high - init_high) < 0.00001) /* init_high (2.0) not high enough, return init_high (2.0) */
	      return init_high;
	    else if (fabs(low - init_low) < 0.00001) /* init_low (0.5) not low enough, return init_low (0.5) */
	      return init_low;
	}
      /* Put the CM back the way it was */
      ExponentiateCM(cm, (1./old_X));
    }
  FreeCM(cm_for_scoring);
  free(cm_minbitsc);

  printf("06.03.07 done expt returning %f\n", X);
  return X;

 ERROR:
  esl_fatal("ERROR in FindExpFactor(), death imminent.\n");
  return 1.0;
}

/*
 * Function: FindExpFactor_minmax()
 * Incept:   EPN, Wed Jun  6 10:14:04 2007
 *
 * Purpose:  Find the exponentiation factor 'X', that when a CM's emission
 *           and transition probabilities are exponentiated by it a specified
 *           fraction range of sample parsetrees have scores greater than
 *           a specified minimum E-value AND less than a specified maximum
 *           E-value. This last point distinguishes this function from 
 *           FindExpFactor() which only uses a minimum E-value. 
 *           Often used prior to calling FindCP9FilterThreshold().
 *
 * Args:
 *           cm           - the CM
 *           cmstats      - CM stats object we'll get Gumbel stats from
 *           r            - source of randomness for EmitParsetree()
 *           cm_min_ecutoff - minimum CM E-value to accept 
 *           cm_max_ecutoff - maximum CM E-value to accept 
 *           db_size      - DB size (nt) to use w/cm_ecutoff to calc CM bit score cutoff 
 *           emit_mode    - CM_LC or CM_GC, CM mode to emit with
 *           fthr_mode    - gives CM search strategy to use, and Gumbel to use
 *           do_fastfil   - TRUE to use fast method: assume parsetree score
 *                          is optimal CYK score
 *           ntrials      - number of parsetrees to sample 
 *           fp_min       - min fraction of ntrials we require b/t cm_min_ecutoff and cm_max_ecutoff
 *           fp_max       - max fraction of ntrials we require b/t cm_min_ecutoff and cm_max_ecutoff
 *           
 * Returns: Exponentiation factor 'X' for which CM^X emits between (ntrails * fp_min) 
 *          and (ntrials * fp_max) parsetrees with scores better than cm_min_ecutoff 
 *          and less than cm_max_ecutoff.
 *          Dies if no such 'X' can be found.
 */
float FindExpFactor_minmax(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
			   float cm_min_ecutoff, float cm_max_ecutoff, int db_size, int emit_mode, 
			   int fthr_mode, int do_fastfil, int ntrials, float fp_min, float fp_max)
{
  /* Contract checks */
  if((fthr_mode != CM_LI) && (fthr_mode != CM_GI) && (fthr_mode != CM_LC) && (fthr_mode != CM_GC))
    Die("ERROR in FindExpFactor_minmax() fthr_mode not CM_LI, CM_GI, CM_LC, or CM_GC\n");
  if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI))
    Die("ERROR in FindExpFactor_minmax() do_fastfil TRUE, but fthr_mode CM_GI or CM_LI\n");
  if(emit_mode != CM_GC && emit_mode != CM_LC)
    Die("ERROR in FindExpFactor_minmax() emit_mode not CM_LC or CM_GC\n");
  if(emit_mode == CM_LC && (fthr_mode == CM_GC || fthr_mode == CM_GI))
    Die("ERROR in FindExpFactor_minmax() emit_mode CM_LC but score mode CM_GC or CM_GI.\n");
  if(fp_min > fp_max)
    Die("ERROR in FindExpFactor_minmax() fp_min > fp_max\n");
  if(fp_min < 0. || fp_min > 1.)
    Die("ERROR in FindExpFactor_minmax() fp_min not between and 0. and 1.\n");
  if(fp_max < 0. || fp_max > 1.)
    Die("ERROR in FindExpFactor_minmax() fp_max not between and 0. and 1.\n");

  int          keep_going = TRUE;
  float        X = 1.0;
  float        init_high = 1.0;
  float        init_low  = 0.8;
  float        high;
  float        low;
  int          X_ctr = 0;
  int          np_min = ntrials * fp_min;
  int          np_max = ntrials * fp_max;
  int          out_half = (ntrials - ((np_min + np_max) / 2.)) / 2.;
  float        sc;
  int          np = 0;
  int          nlow = 0;
  int          nhigh = 0;
  float        old_X;
  float       *cm_minbitsc;
  float       *cm_maxbitsc;
  double        *cm_mu;          /* mu for each partition's CM Gumbel   */
  CM_t        *cm_for_scoring;
  Parsetree_t *tr = NULL;      /* parsetree                           */
  Parsetree_t **hb_trA = NULL;      /* HMM banded parsetrees          */
  char        *dsq = NULL;     /* digitized sequence                  */
  char        *seq = NULL;     /* alphabetic sequence                 */
  int          status;
  int          p;
  float        tmp_K;
  int            i  = 0;         /* counter over samples                */
  int            L;              /* length of a sample                  */
  float          orig_tau;
  int          seen_passed = FALSE;

  high = init_high;
  low  = init_low;
  ESL_ALLOC(cm_minbitsc, sizeof(float) * cmstats->np);
  ESL_ALLOC(cm_maxbitsc, sizeof(float) * cmstats->np);
  ESL_ALLOC(cm_mu,       sizeof(double)* cmstats->np);
  /* Configure the CM based on the emit mode COULD DIFFERENT FROM WORKERS! */
  ConfigForGumbelMode(cm, emit_mode);
  /* Copy the CM into cm_for_scoring, and reconfigure it if nec.,
   * We do this, so we can exponentiate or change emission modes of
   * the original CM */
  cm_for_scoring = DuplicateCM(cm); 
  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
    ConfigForGumbelMode(cm_for_scoring, fthr_mode);

  /* Determine bit cutoff for each partition, calc'ed from cm_min_ecutoff, cm_max_ecutoff */
  for (p = 0; p < cmstats->np; p++)
    {
      /* First determine mu based on db_size */
      tmp_K      = exp(cmstats->gumAA[fthr_mode][p]->mu * cmstats->gumAA[fthr_mode][p]->lambda) / 
	cmstats->gumAA[fthr_mode][p]->L;
      cm_mu[p]   = log(tmp_K  * ((double) db_size)) / cmstats->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cm_minbitsc[p] = cm_mu[p] - (log(cm_min_ecutoff) / cmstats->gumAA[fthr_mode][p]->lambda);
      cm_maxbitsc[p] = cm_mu[p] - (log(cm_max_ecutoff) / cmstats->gumAA[fthr_mode][p]->lambda);
      printf("FindExpFactor() E: %f--%f p: %d %d--%d bit score: %f--%f\n", cm_min_ecutoff, cm_max_ecutoff, p, 
	     cmstats->ps[p], cmstats->pe[p], cm_minbitsc[p], cm_maxbitsc[p]);
      printf("\t ntrials: %d min to pass: %d\n", ntrials, np_min);
    }
  
  /* Do a binary search for maximum 'X' <= 1.0 that gives at least np_min optimal parses with sc better
   * than cm_minsc and less than cm_maxsc. If (do_fastfil), assume parsetree score is optimal. */
  float low_X_that_passes = 0.;
  float step_down_size = 0.05;
  orig_tau = cm_for_scoring->tau;
  while(keep_going)
    {
      printf("06.03.07 Search loop (%d): X: %f\n", X_ctr, X);
      ExponentiateCM(cm, X);
      old_X = X;
      np = nlow = nhigh = 0;
      for(i = 0; i < ntrials; i++)
	{
	  EmitParsetree(cm, r, &tr, &seq, &dsq, &L);
	  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
	    sc = ParsetreeScore_Global2Local(cm_for_scoring, tr, dsq, FALSE);
	  else
	    sc = ParsetreeScore(cm_for_scoring, tr, dsq, FALSE); 
	  p = cmstats->gc2p[(get_gc_comp(seq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */

 	  if(!do_fastfil && (sc < cm_maxbitsc[p]))
	    {
	      printf("T: %10.4f ", sc);

	      /* Accelerate via HMM banded SEARCH */
	      /* For speed, first see if a HMM banded search with big tau gives a bit score above our max */

	      /* STAGE 1 */

	      cm_for_scoring->search_opts |= CM_SEARCH_HBANDED;
	      //cm_for_scoring->tau = 0.1;
	      cm_for_scoring->tau = 0.01;
	      //cm_for_scoring->tau = 0.001;
	      sc = actually_search_target(cm_for_scoring, dsq, 1, L,
					  0.,    /* cutoff is 0 bits (actually we'll find highest
						  * negative score if it's < 0.0) */
					  0.,    /* CP9 cutoff is 0 bits */
					  NULL,  /* don't keep results */
					  FALSE, /* don't filter with a CP9 HMM */
					  FALSE, /* we're not calcing CM  stats */
					  FALSE, /* we're not calcing CP9 stats */
					  NULL); /* filter fraction N/A */
	      printf("SBC: %10.4f %d ", sc, i, np);

	      /* Accelerate via HMM banded ALIGNMENT */
	      /* For speed, first see if an HMM banded alignment gives a bit score above our max */
	      //cm_for_scoring->align_opts |= CM_ALIGN_HBANDED;
	      //cm_for_scoring->align_opts |= CM_ALIGN_NOSMALL;
	      /* Abusing the ESL_SQ structure */
	      //sq[0] = esl_sq_CreateFrom("seq", seq, NULL, NULL, NULL);
	      //asq[0]->dsq = DigitizeSequence (sq[0]->seq, L);
	      //actually_align_targets(cm_for_scoring, sq, 1, &hb_trA, 
	      //NULL, NULL, 0, 0, FALSE);
	      //sc = ParsetreeScore(cm_for_scoring, hb_trA[0], dsq, FALSE);
	      //printf("ABC: %10.4f %d ", sc, i, np);
	      //cm_for_scoring->align_opts &= ~CM_ALIGN_HBANDED;
	      //cm_for_scoring->align_opts &= ~CM_ALIGN_NOSMALL;

	      if(sc < cm_maxbitsc[p])
		{
		  /* do HMM banded search with really small tau */
		  
		  /* Stage 2 */
		  
		  cm_for_scoring->search_opts |= CM_SEARCH_HMMSCANBANDS;
		  cm_for_scoring->tau = 1e-15;
		  //if(FALSE)
		    sc = actually_search_target(cm_for_scoring, dsq, 1, L,
						0.,    /* cutoff is 0 bits (actually we'll find highest
							* negative score if it's < 0.0) */
						0.,    /* CP9 cutoff is 0 bits */
						NULL,  /* don't keep results */
						FALSE, /* don't filter with a CP9 HMM */
						FALSE, /* we're not calcing CM  stats */
						FALSE, /* we're not calcing CP9 stats */
						NULL); /* filter fraction N/A */
		  printf("C: %10.4f i: %d np: %d\n", sc, i, np);
		  
		  if(sc < cm_maxbitsc[p])
		    {
		      /* Stage 3 do QDB CYK */
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HBANDED;
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HMMSCANBANDS;
		      //if(FALSE)
			sc = actually_search_target(cm_for_scoring, dsq, 1, L,
						    0.,    /* cutoff is 0 bits (actually we'll find highest
							    * negative score if it's < 0.0) */
						    0.,    /* CP9 cutoff is 0 bits */
						    NULL,  /* don't keep results */
						    FALSE, /* don't filter with a CP9 HMM */
						    FALSE, /* we're not calcing CM  stats */
						    FALSE, /* we're not calcing CP9 stats */
						    NULL); /* filter fraction N/A */
		    }
	      }
	      else printf("\n");
	    }
	  FreeParsetree(tr);
	  free(seq);
	  free(dsq);

	  if(sc > cm_minbitsc[p] && sc < cm_maxbitsc[p]) np++;
	  if(sc < cm_minbitsc[p]) nlow++;
	  if(sc > cm_maxbitsc[p]) nhigh++;

	  /* If the # passed is at our minimum number and we've seen some high or low, 
	     break loop to save time. */
	  if(np >= np_min && (nlow > 10 || nhigh > 10))
	    i = ntrials;
	}
      printf("06.03.07 %4d X: %.5f nl: %5d np: %5d nh: %5d high: %.5f low: %.5f\n", X_ctr, X, nlow, np, nhigh, high, low);
      X_ctr++;
      if(np >= np_min)
	{
	  seen_passed = TRUE;
	  low_X_that_passes = X;
	  low  = X;  
	  X += (high - X)/2.;  
	  printf("high: %f low: %f (%f)\n", high, low, fabs(high-low));
	}
      else if(nhigh > nlow) 
	{ 
	  high = X;  
	  if(seen_passed) X -= (X  - low)/2.;  
	  else { X -= step_down_size; }
	}
      else
	{ 
	  esl_fatal("ERROR in FindExpFactor_minmax(), nhigh less than nlow\n");
	}
      if(fabs(high-low) < 0.01) 
	{
	  if(seen_passed) /* we've seen an X that passes, set X to that */
	    {
	      X = low_X_that_passes;
	      keep_going = FALSE;
	      printf("06.03.07 success! X: %f\n", X);
	    }
	  else if(fabs(high - init_high) < 0.001) /* init_high (1.0) not high enough, return init_high (2.0) */
	    return init_high;
	  else if (fabs(low - init_low) < 0.001) /* init_low (0.5) not low enough, return init_low (0.5) */
	    return init_low;
	  else
	    {
	      printf("ERROR in FindExpFactor_minmax(), convergence at high: %f low: %f\n", high, low);
	      esl_fatal("done.\n");
	    }
	}
      /* Put the CM back the way it was */
      ExponentiateCM(cm, (1./old_X));
    }
  cm_for_scoring->tau = orig_tau;
  FreeCM(cm_for_scoring);
  free(cm_minbitsc);
  free(cm_maxbitsc);
  free(cm_mu);

  printf("06.03.07 done expt returning %f\n", X);
  return X;

 ERROR:
  esl_fatal("ERROR in FindExpFactor_minmax(), death imminent.\n");
  return 1.0;
}

/*
 * Function: FindExpFactor_min()
 * Incept:   EPN, Mon Jun 11 16:26:47 2007
 *
 * Purpose:  Find the largest exponentiation factor 'X' <= 1.0 that when a CM's 
 *           emission and transition probabilities are exponentiated by it a 
 *           specified fraction range of sample parsetrees have scores WORSE than
 *           a specified minimum E-value. 
 *           Often used prior to calling FindCP9FilterThreshold().
 *
 * Args:
 *           cm           - the CM
 *           cmstats      - CM stats object we'll get Gumbel stats from
 *           r            - source of randomness for EmitParsetree()
 *           cm_min_ecutoff - minimum CM E-value to accept 
 *           db_size      - DB size (nt) to use w/cm_ecutoff to calc CM bit score cutoff 
 *           emit_mode    - CM_LC or CM_GC, CM mode to emit with
 *           fthr_mode    - gives CM search strategy to use, and Gumbel to use
 *           do_fastfil   - TRUE to use fast method: assume parsetree score
 *                          is optimal CYK score (DANGEROUS! DO NOT USE!)
 *           ntrials      - number of parsetrees to sample 
 *           fp_min       - min fraction of ntrials we require b/t cm_min_ecutoff and cm_max_ecutoff
 *           
 * Returns: Exponentiation factor 'X' we're looking for (see Purpose).
 *          Dies if no such 'X' can be found.
 */
float FindExpFactor_min(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
			float cm_min_ecutoff, int db_size, int emit_mode, 
			int fthr_mode, int do_fastfil, int ntrials, float fp_min)
{
  /* Contract checks */
  if((fthr_mode != CM_LI) && (fthr_mode != CM_GI) && (fthr_mode != CM_LC) && (fthr_mode != CM_GC))
    Die("ERROR in FindExpFactor_min() fthr_mode not CM_LI, CM_GI, CM_LC, or CM_GC\n");
  if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI))
    Die("ERROR in FindExpFactor_min() do_fastfil TRUE, but fthr_mode CM_GI or CM_LI\n");
  if(emit_mode != CM_GC && emit_mode != CM_LC)
    Die("ERROR in FindExpFactor_min() emit_mode not CM_LC or CM_GC\n");
  if(emit_mode == CM_LC && (fthr_mode == CM_GC || fthr_mode == CM_GI))
    Die("ERROR in FindExpFactor_min() emit_mode CM_LC but score mode CM_GC or CM_GI.\n");
  if(fp_min < 0. || fp_min > 1.)
    Die("ERROR in FindExpFactor_min() fp_min not between and 0. and 1.\n");

  if(do_fastfil)
    Die("ERROR in FindExpFactor_min() do_fastfil TRUE (is this restriction temporary?)\n");

  int          keep_going = TRUE;
  float        X = 1.0;
  float        init_high = 1.0;
  float        init_low  = 0.5;
  float        high;
  float        low;
  int          X_ctr = 0;
  int          np_min = ntrials * fp_min;
  float        sc;
  int          np = 0;
  int          na = 0;
  float        old_X;
  float       *cm_minbitsc;
  float       *cm_maxbitsc;
  double        *cm_mu;          /* mu for each partition's CM Gumbel   */
  CM_t        *cm_for_scoring;
  Parsetree_t *tr = NULL;      /* parsetree                           */
  Parsetree_t **hb_trA = NULL;      /* HMM banded parsetrees          */
  char        *dsq = NULL;     /* digitized sequence                  */
  char        *seq = NULL;     /* alphabetic sequence                 */
  int          status;
  int          p;
  float        tmp_K;
  int            i  = 0;         /* counter over samples                */
  int            L;              /* length of a sample                  */
  float          orig_tau;
  int          seen_passed = FALSE;
  float low_X_that_passes = 0.;
  float step_down_size = 0.05;
  
  high = init_high;
  low  = init_low;
  ESL_ALLOC(cm_minbitsc, sizeof(float) * cmstats->np);
  ESL_ALLOC(cm_mu,       sizeof(double)* cmstats->np);
  /* Configure the CM based on the emit mode COULD DIFFERENT FROM WORKERS! */
  ConfigForGumbelMode(cm, emit_mode);
  /* Copy the CM into cm_for_scoring, and reconfigure it if nec.,
   * We do this, so we can exponentiate or change emission modes of
   * the original CM */
  cm_for_scoring = DuplicateCM(cm); 
  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
    ConfigForGumbelMode(cm_for_scoring, fthr_mode);

  /* Determine bit cutoff for each partition, calc'ed from cm_min_ecutoff */
  for (p = 0; p < cmstats->np; p++)
    {
      /* First determine mu based on db_size */
      tmp_K      = exp(cmstats->gumAA[fthr_mode][p]->mu * cmstats->gumAA[fthr_mode][p]->lambda) / 
	cmstats->gumAA[fthr_mode][p]->L;
      cm_mu[p]   = log(tmp_K  * ((double) db_size)) / cmstats->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cm_minbitsc[p] = cm_mu[p] - (log(cm_min_ecutoff) / cmstats->gumAA[fthr_mode][p]->lambda);
      printf("06.03.07 FindExpFactor_min() min E: %f p: %d %d--%d bit score: %f", cm_min_ecutoff, p, 
	     cmstats->ps[p], cmstats->pe[p], cm_minbitsc[p]);
      printf("\t ntrials: %d min to pass: %d\n", ntrials, np_min);
    }
  
  /* Do a binary search for maximum 'X' <= 1.0 that gives at least np_min optimal parses with sc WORSE
   * than cm_minsc. If (do_fastfil), assume parsetree score is optimal (WHICH IS DANGEROUS!). */
  orig_tau = cm_for_scoring->tau;
  while(keep_going)
    {
      printf("06.03.07 Search loop (%d): X: %f\n", X_ctr, X);
      ExponentiateCM(cm, X);
      old_X = X;
      np = na = 0;
      for(i = 0; i < ntrials; i++)
	{
	  EmitParsetree(cm, r, &tr, &seq, &dsq, &L);
	  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
	    sc = ParsetreeScore_Global2Local(cm_for_scoring, tr, dsq, FALSE);
	  else
	    sc = ParsetreeScore(cm_for_scoring, tr, dsq, FALSE); 
	  p = cmstats->gc2p[(get_gc_comp(seq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */

 	  if(!do_fastfil && (sc < cm_minbitsc[p]))
	    {
	      printf("T: %10.4f ", sc);

	      /* Accelerate via HMM banded SEARCH */
	      /* For speed, first see if a HMM banded search with big tau gives a bit score above our max */

	      /* STAGE 1 */

	      cm_for_scoring->search_opts |= CM_SEARCH_HBANDED;
	      //cm_for_scoring->tau = 0.1;
	      cm_for_scoring->tau = 0.01;
	      //cm_for_scoring->tau = 0.001;
	      sc = actually_search_target(cm_for_scoring, dsq, 1, L,
					  0.,    /* cutoff is 0 bits (actually we'll find highest
						  * negative score if it's < 0.0) */
					  0.,    /* CP9 cutoff is 0 bits */
					  NULL,  /* don't keep results */
					  FALSE, /* don't filter with a CP9 HMM */
					  FALSE, /* we're not calcing CM  stats */
					  FALSE, /* we're not calcing CP9 stats */
					  NULL); /* filter fraction N/A */
	      printf("SBC: %10.4f %d ", sc, i, np);

	      /* Accelerate via HMM banded ALIGNMENT */
	      /* For speed, first see if an HMM banded alignment gives a bit score above our max */
	      //cm_for_scoring->align_opts |= CM_ALIGN_HBANDED;
	      //cm_for_scoring->align_opts |= CM_ALIGN_NOSMALL;
	      /* Abusing the ESL_SQ structure */
	      //sq[0] = esl_sq_CreateFrom("seq", seq, NULL, NULL, NULL);
	      //asq[0]->dsq = DigitizeSequence (sq[0]->seq, L);
	      //actually_align_targets(cm_for_scoring, sq, 1, &hb_trA, 
	      //NULL, NULL, 0, 0, FALSE);
	      //sc = ParsetreeScore(cm_for_scoring, hb_trA[0], dsq, FALSE);
	      //printf("ABC: %10.4f %d ", sc, i, np);
	      //cm_for_scoring->align_opts &= ~CM_ALIGN_HBANDED;
	      //cm_for_scoring->align_opts &= ~CM_ALIGN_NOSMALL;

	      if(sc < cm_minbitsc[p])
		{
		  /* do HMM banded search with really small tau */
		  
		  /* Stage 2 */
		  
		  cm_for_scoring->search_opts |= CM_SEARCH_HMMSCANBANDS;
		  cm_for_scoring->tau = 1e-15;
		  //if(FALSE)
		    sc = actually_search_target(cm_for_scoring, dsq, 1, L,
						0.,    /* cutoff is 0 bits (actually we'll find highest
							* negative score if it's < 0.0) */
						0.,    /* CP9 cutoff is 0 bits */
						NULL,  /* don't keep results */
						FALSE, /* don't filter with a CP9 HMM */
						FALSE, /* we're not calcing CM  stats */
						FALSE, /* we're not calcing CP9 stats */
						NULL); /* filter fraction N/A */
		    printf("C: %10.4f i: %d np: %d\n", sc, i, np);
		  
		  if(sc < cm_minbitsc[p])
		    {
		      /* Stage 3 do QDB CYK */
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HBANDED;
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HMMSCANBANDS;
		      //if(FALSE)
			sc = actually_search_target(cm_for_scoring, dsq, 1, L,
						    0.,    /* cutoff is 0 bits (actually we'll find highest
							    * negative score if it's < 0.0) */
						    0.,    /* CP9 cutoff is 0 bits */
						    NULL,  /* don't keep results */
						    FALSE, /* don't filter with a CP9 HMM */
						    FALSE, /* we're not calcing CM  stats */
						    FALSE, /* we're not calcing CP9 stats */
						    NULL); /* filter fraction N/A */
		    }
	      }
	      printf("\n");
	    }
	  FreeParsetree(tr);
	  free(seq);
	  free(dsq);
	  if(sc < cm_minbitsc[p]) np++;
	  na++;


	  /* If we've reached our minimum number, break loop to save time. */
	  if(np == np_min)
	    i = ntrials;
	}
      printf("06.03.07 %4d X: %.5f np: %5d na: %5d high: %.5f low: %.5f\n", X_ctr, X, np, na, high, low);
      X_ctr++;
      if(np >= np_min)
	{
	  seen_passed = TRUE;
	  low_X_that_passes = X;
	  low  = X;  
	  X += (high - X)/2.;  
	  printf("high: %f low: %f (%f)\n", high, low, fabs(high-low));
	}
      else
	{ 
	  high = X;  
	  if(seen_passed) X -= (X  - low)/2.;  
	  else { X -= step_down_size; }
	}
      if(fabs(high-low) < 0.005) 
	{
	  if(seen_passed) /* we've seen an X that passes, set X to that */
	    {
	      X = low_X_that_passes;
	      keep_going = FALSE;
	      printf("06.03.07 success! X: %f\n", X);
	    }
	  else if(fabs(high - init_high) < 0.001) /* init_high (1.0) not high enough, return init_high (1.0) */
	    return init_high;
	  else if (fabs(low - init_low) < 0.001) /* init_low (0.5) not low enough, return init_low (0.8) */
	    return init_low;
	  else
	    {
	      printf("ERROR in FindExpFactor_min(), convergence at high: %f low: %f\n", high, low);
	      esl_fatal("done.\n");
	    }
	}
      /* Put the CM back the way it was */
      ExponentiateCM(cm, (1./old_X));
    }
  cm_for_scoring->tau = orig_tau;
  FreeCM(cm_for_scoring);
  free(cm_minbitsc);
  free(cm_mu);

  printf("06.03.07 done expt returning %f\n", X);
  return X;

 ERROR:
  esl_fatal("ERROR in FindExpFactor_min(), death imminent.\n");
  return 1.0;
}


/*
 * Function: Filter_XFTableLookup()
 * Incept:   EPN, Tue Jun  5 07:07:39 2007
 *
 * Purpose:  Given a fraction 'F' of CM hits the HMM filter should find and
 *           the exponentiation factor 'X' used to make the CM emit parsetrees
 *           with suitable scores, lookup and return the corrected 'F' using
 *           a hardcoded map of empirically determined values, and extrapolating
 *           where necessary.      
 *
 * Args:
 *           X            - the exponentiation factor
 *           F            - desired fraction of CM hits HMM should find (0.95<=x<=0.99)
 *           emit_mode    - CM_LC or CM_GC, CM mode to emit with
 *           fthr_mode    - gives CM search strategy to use, and Gumbel to use
 *           
 * Returns: Corrected 'F' from table after extrapolation as explained above.
 */
float Filter_XFTableLookup(float X, float F, int emit_mode, int fthr_mode)
{
  printf("in Filter_XFTableLookup() X: %f F: %f emit_mode: %d fthr_mode: %d\n", X, F, emit_mode, fthr_mode);
  /* Contract checks (some are temporary) */
  if(emit_mode != CM_GC) 
    esl_fatal("ERROR in Filter_XFTableLookup(), emit mode not CM_GC.\n");
  if(fthr_mode != CM_LC && fthr_mode != CM_LI) 
    esl_fatal("ERROR in Filter_XFTableLookup(), fthr mode not CM_LC or CM_LI.\n");
  if(F > 0.990001 || F < 0.94999)
    esl_fatal("ERROR in Filter_XFTableLookup(), F must be in range [0.95..0.99]\n");
  if(X < 0.75)
    esl_fatal("in Filter_XTableLookup(), X is below 0.75\n");

  if(X > 1.249) 
    { 
      printf("in Filter_XTableLookup(), exceeds 1.249, returning 1.0.\n");
      return 1.0;
    }

  int status;
  double **XF_AA;
  int i;
  ESL_ALLOC(XF_AA, sizeof(double *) * 5);
  for(i = 0; i < 5; i++)
    {
      ESL_ALLOC(XF_AA[i], sizeof(double) * 50);
    }

  /* from ~/notebook/7_0520_inf_tail_emit/expts_0603/rmark-1/X_close_to_1/rmark_close21.table.code */
  XF_AA[0][0] = 0.604515423192598;
  XF_AA[0][1] = 0.614831485532718;
  XF_AA[0][2] = 0.640483722291663;
  XF_AA[0][3] = 0.653194799326666;
  XF_AA[0][4] = 0.664721613516714;
  XF_AA[0][5] = 0.684792103312101;
  XF_AA[0][6] = 0.696334432588007;
  XF_AA[0][7] = 0.717301102410688;
  XF_AA[0][8] = 0.738809779743004;
  XF_AA[0][9] = 0.748424707090713;
  XF_AA[0][10] = 0.764510665501093;
  XF_AA[0][11] = 0.776107759845839;
  XF_AA[0][12] = 0.793102525868823;
  XF_AA[0][13] = 0.802242539634606;
  XF_AA[0][14] = 0.817396549001533;
  XF_AA[0][15] = 0.833510018143928;
  XF_AA[0][16] = 0.845054834049552;
  XF_AA[0][17] = 0.859990635768138;
  XF_AA[0][18] = 0.873048725966369;
  XF_AA[0][19] = 0.88312357348995;
  XF_AA[0][20] = 0.893851543442886;
  XF_AA[0][21] = 0.899095896617135;
  XF_AA[0][22] = 0.91290439255526;
  XF_AA[0][23] = 0.919998809374404;
  XF_AA[0][24] = 0.929454784999846;
  XF_AA[0][25] = 0.939837938984853;
  XF_AA[0][26] = 0.945430601430853;
  XF_AA[0][27] = 0.952149230210086;
  XF_AA[0][28] = 0.961300521748108;
  XF_AA[0][29] = 0.966066709815535;
  XF_AA[0][30] = 0.970570392014799;
  XF_AA[0][31] = 0.975461211728663;
  XF_AA[0][32] = 0.977603094876106;
  XF_AA[0][33] = 0.981002374018307;
  XF_AA[0][34] = 0.984292259859661;
  XF_AA[0][35] = 0.985812772516189;
  XF_AA[0][36] = 0.988604697815925;
  XF_AA[0][37] = 0.991247813621203;
  XF_AA[0][38] = 0.993506975192497;
  XF_AA[0][39] = 0.994064191924624;
  XF_AA[0][40] = 0.995988455979326;
  XF_AA[0][41] = 0.996120992797291;
  XF_AA[0][42] = 0.996524405956255;
  XF_AA[0][43] = 0.997289275412066;
  XF_AA[0][44] = 0.998140930863923;
  XF_AA[0][45] = 0.99811346203176;
  XF_AA[0][46] = 0.999321588371466;
  XF_AA[0][47] = 0.998690299422014;
  XF_AA[0][48] = 0.998872079106245;
  XF_AA[0][49] = 0.99928186466236;

  XF_AA[1][0] = 0.651507805191574;
  XF_AA[1][1] = 0.660560549957533;
  XF_AA[1][2] = 0.686186803189465;
  XF_AA[1][3] = 0.696858774402236;
  XF_AA[1][4] = 0.706405154190228;
  XF_AA[1][5] = 0.725072729101922;
  XF_AA[1][6] = 0.735330802994684;
  XF_AA[1][7] = 0.75508335639012;
  XF_AA[1][8] = 0.775144236530571;
  XF_AA[1][9] = 0.784491364559886;
  XF_AA[1][10] = 0.79852241987608;
  XF_AA[1][11] = 0.807965324841711;
  XF_AA[1][12] = 0.824217530792584;
  XF_AA[1][13] = 0.832024779137886;
  XF_AA[1][14] = 0.844790310291407;
  XF_AA[1][15] = 0.859582957980509;
  XF_AA[1][16] = 0.86894693564019;
  XF_AA[1][17] = 0.882236426601724;
  XF_AA[1][18] = 0.893328487923062;
  XF_AA[1][19] = 0.901919407618695;
  XF_AA[1][20] = 0.912516845647004;
  XF_AA[1][21] = 0.915903778251477;
  XF_AA[1][22] = 0.928925101621968;
  XF_AA[1][23] = 0.933522379705056;
  XF_AA[1][24] = 0.941902240181435;
  XF_AA[1][25] = 0.950084243608121;
  XF_AA[1][26] = 0.955767251863576;
  XF_AA[1][27] = 0.961124663003683;
  XF_AA[1][28] = 0.968520148816714;
  XF_AA[1][29] = 0.97275666968008;
  XF_AA[1][30] = 0.976555605591042;
  XF_AA[1][31] = 0.980210559693913;
  XF_AA[1][32] = 0.982099840185598;
  XF_AA[1][33] = 0.984546921628583;
  XF_AA[1][34] = 0.987191158213102;
  XF_AA[1][35] = 0.988905332488261;
  XF_AA[1][36] = 0.991035950680023;
  XF_AA[1][37] = 0.993250800528251;
  XF_AA[1][38] = 0.995293396435803;
  XF_AA[1][39] = 0.995605041935139;
  XF_AA[1][40] = 0.996507055729084;
  XF_AA[1][41] = 0.996684719302554;
  XF_AA[1][42] = 0.997785194884525;
  XF_AA[1][43] = 0.997813869807061;
  XF_AA[1][44] = 0.998851091704242;
  XF_AA[1][45] = 0.998786788290816;
  XF_AA[1][46] = 0.999631834629549;
  XF_AA[1][47] = 0.99782431444757;
  XF_AA[1][48] = 0.999083449547723;
  XF_AA[1][49] = 0.998939543116936;

  XF_AA[2][0] = 0.698068831381895;
  XF_AA[2][1] = 0.706196770587634;
  XF_AA[2][2] = 0.730030559690639;
  XF_AA[2][3] = 0.740250612982444;
  XF_AA[2][4] = 0.749767542725801;
  XF_AA[2][5] = 0.766840329536004;
  XF_AA[2][6] = 0.775252201727982;
  XF_AA[2][7] = 0.791592319198176;
  XF_AA[2][8] = 0.81069967397328;
  XF_AA[2][9] = 0.819699794556972;
  XF_AA[2][10] = 0.831920771972702;
  XF_AA[2][11] = 0.839637857470845;
  XF_AA[2][12] = 0.853415720626429;
  XF_AA[2][13] = 0.858613091377097;
  XF_AA[2][14] = 0.87233921186971;
  XF_AA[2][15] = 0.884409633343452;
  XF_AA[2][16] = 0.893122940535288;
  XF_AA[2][17] = 0.90361957928424;
  XF_AA[2][18] = 0.912921961208646;
  XF_AA[2][19] = 0.920095554050192;
  XF_AA[2][20] = 0.929290051807059;
  XF_AA[2][21] = 0.931576686768534;
  XF_AA[2][22] = 0.942707955664224;
  XF_AA[2][23] = 0.946248893484383;
  XF_AA[2][24] = 0.952945812285395;
  XF_AA[2][25] = 0.960442605946639;
  XF_AA[2][26] = 0.965279267778482;
  XF_AA[2][27] = 0.96882712601592;
  XF_AA[2][28] = 0.97614984985738;
  XF_AA[2][29] = 0.97952055184705;
  XF_AA[2][30] = 0.981374975272281;
  XF_AA[2][31] = 0.984755218232072;
  XF_AA[2][32] = 0.986447558547585;
  XF_AA[2][33] = 0.988059530318827;
  XF_AA[2][34] = 0.989895930683065;
  XF_AA[2][35] = 0.99173856432554;
  XF_AA[2][36] = 0.992792458352072;
  XF_AA[2][37] = 0.994776057675342;
  XF_AA[2][38] = 0.996653380435322;
  XF_AA[2][39] = 0.996358177774394;
  XF_AA[2][40] = 0.99767166216276;
  XF_AA[2][41] = 0.99807797054367;
  XF_AA[2][42] = 0.998490407655929;
  XF_AA[2][43] = 0.998625900055842;
  XF_AA[2][44] = 0.999032032449554;
  XF_AA[2][45] = 0.998368634332335;
  XF_AA[2][46] = 0.999493938736146;
  XF_AA[2][47] = 0.997961464856902;
  XF_AA[2][48] = 0.99926757254451;
  XF_AA[2][49] = 0.99919760044068;

  XF_AA[3][0] = 0.757724305104733;
  XF_AA[3][1] = 0.766461466775073;
  XF_AA[3][2] = 0.784979752178282;
  XF_AA[3][3] = 0.79434397402588;
  XF_AA[3][4] = 0.802017008556587;
  XF_AA[3][5] = 0.817887452755621;
  XF_AA[3][6] = 0.824458070452807;
  XF_AA[3][7] = 0.837354545237931;
  XF_AA[3][8] = 0.854769767428774;
  XF_AA[3][9] = 0.862213006604172;
  XF_AA[3][10] = 0.873201595298575;
  XF_AA[3][11] = 0.877962153736902;
  XF_AA[3][12] = 0.888494760978758;
  XF_AA[3][13] = 0.8929426934511;
  XF_AA[3][14] = 0.904267090291344;
  XF_AA[3][15] = 0.913804159254234;
  XF_AA[3][16] = 0.919535282202327;
  XF_AA[3][17] = 0.928230503690814;
  XF_AA[3][18] = 0.935307335667513;
  XF_AA[3][19] = 0.941256475828225;
  XF_AA[3][20] = 0.949030568488772;
  XF_AA[3][21] = 0.949945934709077;
  XF_AA[3][22] = 0.957970528251316;
  XF_AA[3][23] = 0.961090198134555;
  XF_AA[3][24] = 0.966164195654882;
  XF_AA[3][25] = 0.971684199746855;
  XF_AA[3][26] = 0.976030330408388;
  XF_AA[3][27] = 0.977793072539259;
  XF_AA[3][28] = 0.98379319029689;
  XF_AA[3][29] = 0.985342821244426;
  XF_AA[3][30] = 0.988342683165911;
  XF_AA[3][31] = 0.989854392548513;
  XF_AA[3][32] = 0.990699980929024;
  XF_AA[3][33] = 0.99195007325228;
  XF_AA[3][34] = 0.993298367742526;
  XF_AA[3][35] = 0.994398012140472;
  XF_AA[3][36] = 0.995394140596778;
  XF_AA[3][37] = 0.996807787452819;
  XF_AA[3][38] = 0.99797568427179;
  XF_AA[3][39] = 0.997979619235294;
  XF_AA[3][40] = 0.998438792467191;
  XF_AA[3][41] = 0.998581294266464;
  XF_AA[3][42] = 0.998972976932633;
  XF_AA[3][43] = 0.99912961663165;
  XF_AA[3][44] = 0.999416928771567;
  XF_AA[3][45] = 0.998799776518357;
  XF_AA[3][46] = 0.999490598632985;
  XF_AA[3][47] = 0.998556001944389;
  XF_AA[3][48] = 0.999248894097374;
  XF_AA[3][49] = 0.999438980162221;

  XF_AA[4][0] = 0.835926947200063;
  XF_AA[4][1] = 0.840096285366306;
  XF_AA[4][2] = 0.854913111501865;
  XF_AA[4][3] = 0.861427817213593;
  XF_AA[4][4] = 0.868913860855596;
  XF_AA[4][5] = 0.880673045371533;
  XF_AA[4][6] = 0.885221733032393;
  XF_AA[4][7] = 0.892556750391326;
  XF_AA[4][8] = 0.906645038425893;
  XF_AA[4][9] = 0.912782761250072;
  XF_AA[4][10] = 0.920695351777993;
  XF_AA[4][11] = 0.923606158847454;
  XF_AA[4][12] = 0.929654991028677;
  XF_AA[4][13] = 0.931456550407444;
  XF_AA[4][14] = 0.941009532409196;
  XF_AA[4][15] = 0.947452390809934;
  XF_AA[4][16] = 0.951219782862265;
  XF_AA[4][17] = 0.955998843984201;
  XF_AA[4][18] = 0.960581845202378;
  XF_AA[4][19] = 0.964765994127025;
  XF_AA[4][20] = 0.969891036795612;
  XF_AA[4][21] = 0.969916152140906;
  XF_AA[4][22] = 0.976101734678514;
  XF_AA[4][23] = 0.976799430520982;
  XF_AA[4][24] = 0.979214706351899;
  XF_AA[4][25] = 0.982674464196227;
  XF_AA[4][26] = 0.986092249453835;
  XF_AA[4][27] = 0.986582538247042;
  XF_AA[4][28] = 0.991199061348433;
  XF_AA[4][29] = 0.992459929940093;
  XF_AA[4][30] = 0.993635344813326;
  XF_AA[4][31] = 0.994438724023198;
  XF_AA[4][32] = 0.994970564163626;
  XF_AA[4][33] = 0.995279388113399;
  XF_AA[4][34] = 0.996209951784272;
  XF_AA[4][35] = 0.996779356160901;
  XF_AA[4][36] = 0.997499207959605;
  XF_AA[4][37] = 0.998394104286069;
  XF_AA[4][38] = 0.999293502799076;
  XF_AA[4][39] = 0.998760821386919;
  XF_AA[4][40] = 0.999366130969548;
  XF_AA[4][41] = 0.999565664979706;
  XF_AA[4][42] = 0.999219295380022;
  XF_AA[4][43] = 0.999533070445127;
  XF_AA[4][44] = 0.999348779659682;
  XF_AA[4][45] = 0.998307864255877;
  XF_AA[4][46] = 0.999540316530718;
  XF_AA[4][47] = 0.999045950019794;
  XF_AA[4][48] = 0.998854972879331;
  XF_AA[4][49] = 0.999457304554798;

  /* Find appropriate values that bracket X above and below and smallest F
   * above F, and extrapolate  */
  float F_below = 0.95;
  float F_above;
  int i_x = 0;
  int i_f = 0;
  int keep_going = TRUE;
  int i_f_below, i_f_above;
  float X_below;
  float return_F;

  while(keep_going && i_f < 5)
    {
      if(fabs(F-F_below) > 0.0001) { i_f++; F_below += 0.01; }
      else 
	{ 
	  keep_going = FALSE;
	  i_f_below = i_f-1;
	  if(i_x == 0) i_f_below = i_f;
	  i_f_above = i_f;
	  F_above = F_below;
	  F_below -= 0.01;
	}
    }
  if(keep_going) 
    {
      printf("couldn't determine F_below F_above for F: %f\n", F);
      esl_fatal("done.\n");
    }
  printf("i_f: %d F_above: %f F_below: %f\n", i_f, F_above, F_below);

  i_x = 0;
  X_below = 0.75;
  keep_going = TRUE;
  while(keep_going && i_x < 50)
    {
      if(X > X_below) { i_x++; X_below += 0.01; }
      else keep_going = FALSE;
    }
  if(keep_going) 
    {
      printf("couldn't determine X_below for X: %f\n", X);
      esl_fatal("done.\n");
    }
  printf("i_x: %d X_below: %f\n", i_x, X_below);

  return_F = XF_AA[i_f_below][i_x] + (((F - F_below) / (F_above-F_below)) * 
				      (XF_AA[i_f_above] - XF_AA[i_f_below]));
  printf("in Filter_XFTableLookup() F: %f X: %f returning: %f\n", F, X, return_F);
  return return_F;

 ERROR:
  esl_fatal("EASEL ERROR\n");
  return 0.;
}
  

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
  mx = AllocCPlan9Matrix(nrows, cm->cp9->M, &mmx, &imx, &dmx, 
			 NULL); /* we don't use emx any more, isc replaces it,
				 * 1. b/c we want isc for each position, but
				 * we may only be alloc'ing 2 rows of mx,
				 * 2. b/c emx is 2D for dumb memory error reasons,
				 * which needs to be fixed. */
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
		  /*printf("FWD greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i, j, fsc);*/
		  report_hit (i, j, 0, fsc, results);
		  /* 0 is for saver, which is irrelevant for HMM hits */
		}
	    }
	}
      if (fsc > best_hmm_sc)
	{
	  best_hmm_sc = fsc;
	  best_hmm_pos= j;
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
      return_sc = Scorify(sc[(j0-i0+1)]); /* L = i0-j0+1 */
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
  mx = AllocCPlan9Matrix(nrows, cm->cp9->M, &mmx, &imx, &dmx, 
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
      if(!(cm->search_opts & CM_SEARCH_HMMONLY))
	printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);
    }
  FreeResults (fwd_results);
  FreeResults (bwd_results);
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
  float ffrac;
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

/*
 * Function: FindCP9FilterThreshold()
 * Incept:   EPN, Wed May  2 10:00:45 2007
 *
 * Purpose:  Sample sequences from a CM and determine the CP9 HMM E-value
 *           threshold necessary to recognize a specified fraction of those
 *           hits. Sequences are sampled from the CM until N with a E-value
 *           better than cm_ecutoff are sampled (those with worse E-values
 *           are rejected). CP9 scans are carried out in either local or
 *           glocal mode depending on hmm_evd_mode. CM is configured in 
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
 *           cmstats      - CM stats object we'll get EVD stats from
 *           fraction     - target fraction of CM hits to detect with CP9
 *           N            - number of sequences to sample from CM better than cm_minsc
 *           use_cm_cutoff- TRUE to only accept CM parses w/E-values above cm_ecutoff
 *           cm_ecutoff   - minimum CM E-value to accept 
 *           db_size      - DB size (nt) to use w/cm_ecutoff to calc CM bit score cutoff 
 *           fthr_mode    - gives CM search strategy to use, and EVD to use
 *           hmm_evd_mode - CP9_L to search with  local HMM (we're filling a fthr->l_eval)
 *                          CP9_G to search with glocal HMM (we're filling a fthr->g_eval)
 *           do_fastfil   - TRUE to use fast method: assume parsetree score
 *                          is optimal CYK score
 * 
 * Returns: HMM E-value cutoff above which the HMM scores (fraction * N) CM 
 *          hits with CM E-values better than cm_ecutoff 
 * 
 */
float FindCP9FilterThreshold(CM_t *cm, CMStats_t *cmstats, float fraction, int N, 
			     int use_cm_cutoff, float cm_ecutoff, int db_size, 
			     int fthr_mode, int hmm_evd_mode, int do_fastfil)
{
  Parsetree_t      *tr;         /* parsetree (TEMPORARY NOT REALLY NEEDED)*/
  char             *dsq;        /* digitized sequence                     */
  char             *seq;        /* alphabetic sequence                    */
  float             sc;
  float            *hmm_eval;
  float            *hmm_sc;
  int               nattempts = 0;
  int               max_attempts;
  int               L;
  int               i;
  int               gc_comp;
  EVDInfo_t **cm_evd; 
  EVDInfo_t **hmm_evd; 
  double *cm_mu;
  double *hmm_mu;
  int passed;
  int p;                    /* counter over partitions */
  double pval;
  double x;     
  float *cm_minbitsc = NULL;/* minimum CM bit score to allow to pass for each partition */
  float return_eval;
  double tmp_K;

  /* Contract check */
  if (!(cm->flags & CM_CP9) || cm->cp9 == NULL) 
    Die("ERROR in FindCP9FilterThreshold() CP9 does not exist\n");
  if (fraction < 0. || fraction > 1.)  
    Die("ERROR in FindCP9FilterThreshold() fraction is %f, should be [0.0..1.0]\n", fraction);
  if((fthr_mode != CM_LI) && (fthr_mode != CM_GI) && (fthr_mode != CM_LC) && (fthr_mode != CM_GC))
    Die("ERROR in FindCP9FilterThreshold() fthr_mode not CM_LI, CM_GI, CM_LC, or CM_GC\n");
  if(hmm_evd_mode != CP9_L && hmm_evd_mode != CP9_G)
    Die("ERROR in FindCP9FilterThreshold() hmm_evd_mode not CP9_L or CP9_G\n");
  if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI))
    Die("ERROR in FindCP9FilterThreshold() do_fastfil TRUE, but fthr_mode CM_GI or CM_LI\n");

  /* Configure the CM based on the stat mode */
  ConfigForEVDMode(cm, fthr_mode);
  cm_evd   = cmstats->evdAA[fthr_mode];
  hmm_evd  = cmstats->evdAA[hmm_evd_mode];

  hmm_eval  = MallocOrDie(sizeof(float) * N);
  hmm_sc    = MallocOrDie(sizeof(float) * N);
  max_attempts = 500 * N;
  seq    = NULL;
  dsq    = NULL;
  tr     = NULL;
  /* Configure the HMM based on the hmm_evd_mode */
  if(hmm_evd_mode == CP9_L)
    CPlan9SWConfig(cm->cp9, ((cm->cp9->M)-1.)/cm->cp9->M,  /* all start pts equiprobable, including 1 */
		   ((cm->cp9->M)-1.)/cm->cp9->M);          /* all end pts equiprobable, including M */
  else /* hmm_evd_mode == CP9_G (it's in the contract) */
    CPlan9GlobalConfig(cm->cp9);
  CP9Logoddsify(cm->cp9);

  if(use_cm_cutoff) printf("CM E cutoff: %f\n", cm_ecutoff);
  else              printf("Not using CM cutoff\n");

  /* Determine bit cutoff for each partition, calc'ed from cm_ecutoff */
  cm_minbitsc = MallocOrDie(sizeof(float)  * cmstats->np);
  cm_mu       = MallocOrDie(sizeof(double) * cmstats->np);
  hmm_mu      = MallocOrDie(sizeof(double) * cmstats->np);
  for (p = 0; p < cmstats->np; p++)
    {
      /* First determine mu based on db_size */
      tmp_K      = exp(hmm_evd[p]->mu * hmm_evd[p]->lambda) / hmm_evd[p]->L;
      printf("HMM mu: %f lambda: %f K: %f\n", hmm_evd[p]->mu, hmm_evd[p]->lambda, tmp_K);
      hmm_mu[p]  = log(tmp_K * ((double) db_size)) / hmm_evd[p]->lambda;
      printf("HMM new mu: %f\n", hmm_mu[p]);
      tmp_K      = exp(cm_evd[p]->mu * cm_evd[p]->lambda) / cm_evd[p]->L;
      printf("CM mu: %f lambda: %f K: %f\n", cm_evd[p]->mu, cm_evd[p]->lambda, tmp_K);
      cm_mu[p]   = log(tmp_K  * ((double) db_size)) / cm_evd[p]->lambda;
      printf("CM new mu: %f\n", cm_mu[p]);
      /* Now determine bit score */
      cm_minbitsc[p] = cm_mu[p] - (log(cm_ecutoff) / cm_evd[p]->lambda);
      if(use_cm_cutoff)
	printf("E: %f p: %d %d--%d bit score: %f\n", cm_ecutoff, p, 
	       cmstats->ps[p], cmstats->pe[p], cm_minbitsc[p],
	       RJK_ExtremeValueE(cm_minbitsc[p], cm_evd[p]->mu, cm_evd[p]->lambda));
    }

  /* Strategy: 
   * Emit sequences one at a time from CM, scanning each with CM 
   * (local or glocal, CYK or Inside as specified by fthr_mode). If best scan
   * score meets cm_ecutoff threshold, search it with CP9 in local or glocal 
   * mode (specified by hmm_evd_mode) and save top CP9 score. Repeat
   * until N seqs have been searched with CP9.
   */

  /* Emit seqs until we have N that meet score cutoff */
  for (i = 0; i < N; i++)
    {
      if(nattempts++ > max_attempts) 
	Die("ERROR number of attempts exceeded 500 times number of seqs.\n");
      passed = FALSE;
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      while(!passed && nattempts <= max_attempts)
	{
	  if(dsq != NULL) free(dsq);
	  if(seq != NULL) free(seq);
	  if(tr != NULL) FreeParsetree(tr);
	  EmitParsetree(cm, &tr, &seq, &dsq, &L); /* we don't care about the parsetree, 
						   * we have to scan each seq to get score */
	  if(do_fastfil) 
	    sc = ParsetreeScore(cm, tr, dsq, FALSE); 
	  else
	    sc = actually_search_target(cm, dsq, 1, L,
					0.,    /* cutoff is 0 bits (actually we'll find highest
						* negative score if it's < 0.0) */
					0.,    /* CP9 cutoff is 0 bits */
					NULL,  /* don't keep results */
					FALSE, /* don't filter with a CP9 HMM */
					FALSE, /* we're not calcing CM  stats */
					FALSE, /* we're not calcing CP9 stats */
					NULL); /* filter fraction N/A */
	  /* Only accept if E-value of best hit in this seq is better than our cutoff.
	   * To do that we first determine which partition to use EVD from */
	  gc_comp = get_gc_comp(seq, 1, L); /* should be i and j of best hit, but
					     * EVD construction is wrong too, it uses 
					     * GC_COMP of full sample seq, not of best hit
					     * within it. */
	  p    = cmstats->gc2p[gc_comp];
	  if(sc >= cm_minbitsc[p])
	    passed = TRUE;
	  nattempts++;
	}
      if(nattempts > max_attempts)
	Die("ERROR number of attempts exceeded 50 times number of seqs.\n");
      /* Scan seq with HMM */
      cm->search_opts |= CM_SEARCH_HMMONLY;
      hmm_sc[i] = actually_search_target(cm, dsq, 1, L,
					 cm->cutoff, cm->cp9_cutoff,
					 NULL,  /* don't report hits to a results structure */
					 FALSE, /* we're not filtering with a CP9 HMM */
					 FALSE, /* we're not building a histogram for CM stats  */
					 FALSE, /* we're not building a histogram for CP9 stats */
					 NULL); /* filter fraction, irrelevant here */
      /*printf("hmm_mu[p]: %f hmm_evd[p]->mu: %f\n", hmm_mu[p], hmm_evd[p]->mu);*/
      hmm_eval[i] = RJK_ExtremeValueE(hmm_sc[i], hmm_mu[p], hmm_evd[p]->lambda);
      printf("sc: %f hmm_eval[i]: %f orig eval: %f ", hmm_sc[i], hmm_eval[i], 
	RJK_ExtremeValueE(hmm_sc[i], hmm_evd[p]->mu, hmm_evd[p]->lambda));
      printf("hmm P: %f\n", esl_gumbel_surv((double) hmm_sc[i], hmm_mu[p], hmm_evd[p]->lambda));
   }
  /* Sort the HMM E-values with quicksort */
  esl_vec_FSortIncreasing(hmm_eval, N);
  esl_vec_FSortDecreasing(hmm_sc, N);

  for (i = 0; i < N; i++)
    printf("%d i: %4d hmm sc: %10.4f hmm E: %10.4f\n", hmm_evd_mode, i, hmm_sc[i], hmm_eval[i]);
  printf("\n\nnattempts: %d\n", nattempts);

  /* Clean up and exit */
  return_eval = hmm_eval[(int) ((fraction) * (float) N)];
  if(return_eval > ((float) db_size)) /* E-val > db_size is useless */
    return_eval = (float) db_size;
  free(hmm_eval);
  free(hmm_sc);
  free(hmm_mu);
  free(cm_mu);
  free(cm_minbitsc);
  if(dsq != NULL) free(dsq);
  if(seq != NULL) free(seq);
  /* Return threshold */
  return return_eval;
}

#if USE_MPI
/*
 * Function: mpi_FindCP9FilterThreshold()
 * Incept:   EPN, Thu May 10 09:03:25 2007
 *
 * Purpose:  MPI version of FindCP9FilterThreshold(). 
 *
 * Args:
 *           cm           - the CM
 *           cmstats      - CM stats object we'll get EVD stats from
 *           fraction     - target fraction of CM hits to detect with CP9
 *           N            - number of sequences to sample from CM better than cm_minsc
 *           use_cm_cutoff- TRUE to only accept CM parses w/E-values above cm_ecutoff
 *           cm_ecutoff   - minimum CM E-value to accept 
 *           db_size      - DB size (nt) to use w/cm_ecutoff to calc CM bit score cutoff 
 *           fthr_mode    - gives CM search strategy to use, and EVD to use
 *           hmm_evd_mode - CP9_L to search with  local HMM (we're filling a fthr->l_eval)
 *                          CP9_G to search with glocal HMM (we're filling a fthr->g_eval)
 *           do_fastfil   - TRUE to use fast method: assume parsetree score
 *                          is optimal CYK score
 *           my_rank      - my MPI rank, 0 is master, >0 is worker
 *           nproc       - number of MPI processors
 *
 * Returns: HMM E-value cutoff above which the HMM scores (fraction * N) CM 
 *          hits with CM E-values better than cm_ecutoff 
 * 
 */
float mpi_FindCP9FilterThreshold(CM_t *cm, CMStats_t *cmstats, float fraction, int N, 
				 int use_cm_cutoff, float cm_ecutoff, int db_size, 
				 int fthr_mode, int hmm_evd_mode, int do_fastfil,
				 int my_rank, int nproc)
{
  /* Contract checks */
  if (!(cm->flags & CM_CP9) || cm->cp9 == NULL) 
    Die("ERROR in mpi_FindCP9FilterThreshold() CP9 does not exist\n");
  if (fraction < 0. || fraction > 1.)  
    Die("ERROR in mpi_FindCP9FilterThreshold() fraction is %f, should be [0.0..1.0]\n", fraction);
  if((fthr_mode != CM_LI) && (fthr_mode != CM_GI) && (fthr_mode != CM_LC) && (fthr_mode != CM_GC))
    Die("ERROR in mpi_FindCP9FilterThreshold() fthr_mode not CM_LI, CM_GI, CM_LC, or CM_GC\n");
  if(hmm_evd_mode != CP9_L && hmm_evd_mode != CP9_G)
    Die("ERROR in mpi_FindCP9FilterThreshold() hmm_evd_mode not CP9_L or CP9_G\n");
  if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI))
    Die("ERROR in mpi_FindCP9FilterThreshold() do_fastfil TRUE, but fthr_mode CM_GI or CM_LI\n");

  /* Configure the CM based on the stat mode */
  ConfigForEVDMode(cm, fthr_mode);
  /* Configure the HMM based on the hmm_evd_mode */
  if(hmm_evd_mode == CP9_L)
    CPlan9SWConfig(cm->cp9, ((cm->cp9->M)-1.)/cm->cp9->M,  /* all start pts equiprobable, including 1 */
		   ((cm->cp9->M)-1.)/cm->cp9->M);          /* all end pts equiprobable, including M */
  else /* hmm_evd_mode == CP9_G (it's in the contract) */
    CPlan9GlobalConfig(cm->cp9);
  CP9Logoddsify(cm->cp9);
      
  if(my_rank > 0)
    {
      mpi_worker_cm_and_cp9_search(cm, do_fastfil, my_rank);
      return eslOK;
    }
  else /* my_rank == 0, master */
    {
      /* GETS HERE IN MASTER! */
      int status;
      Parsetree_t      *tr;         /* parsetree */
      char             *dsq;        /* digitized sequence                     */
      char             *seq;        /* alphabetic sequence                    */
      float             sc;
      float            *hmm_eval;
      float            *hmm_sc;
      int               nattempts = 0;
      int               max_attempts;
      int               L;
      int               i;
      int               gc_comp;
      EVDInfo_t **cm_evd; 
      EVDInfo_t **hmm_evd; 
      double *cm_mu;
      double *hmm_mu;
      double  tmp_K;
      int passed;
      int p;                    /* counter over partitions */
      double pval;
      double x;     
      float *cm_minbitsc = NULL;/* minimum CM bit score to allow to pass for each partition */
      float return_eval;
      int              have_work;
      int              nproc_working;
      int              wi;
      char        **dsqlist = NULL;      /* queue of digitized seqs being worked on, 1..nproc-1 */
      int            *plist = NULL;      /* queue of partition indices of seqs being worked on, 1..nproc-1 */
      Parsetree_t  **trlist = NULL;      /* queue of traces of seqs being worked on, 1..nproc-1 */
      float *scores;
      float cur_cm_sc;
      float cur_hmm_sc;

      MPI_Status       mstatus;
      
      cm_evd   = cmstats->evdAA[fthr_mode];
      hmm_evd  = cmstats->evdAA[hmm_evd_mode];

      ESL_ALLOC(scores,      sizeof(float)         * 2);
      ESL_ALLOC(dsqlist,     sizeof(char *)        * nproc);
      ESL_ALLOC(plist,       sizeof(int)           * nproc);
      ESL_ALLOC(trlist,      sizeof(Parsetree_t *) * nproc);
      ESL_ALLOC(hmm_eval,    sizeof(float)         * N);
      ESL_ALLOC(hmm_sc,      sizeof(float)         * N);
      ESL_ALLOC(cm_minbitsc, sizeof(float)         * cmstats->np);
      ESL_ALLOC(cm_mu,       sizeof(float)         * cmstats->np);
      ESL_ALLOC(hmm_mu,      sizeof(float)         * cmstats->np);

      max_attempts = 500 * N;
      seq    = NULL;
      dsq    = NULL;
      tr     = NULL;

      if(use_cm_cutoff) printf("CM E cutoff: %f\n", cm_ecutoff);
      else              printf("Not using CM cutoff\n");
      
      /* Determine bit cutoff for each partition, calc'ed from cm_ecutoff */
      for (p = 0; p < cmstats->np; p++)
	{
	  /* First determine mu based on db_size */
	  tmp_K      = exp(hmm_evd[p]->mu * hmm_evd[p]->lambda) / hmm_evd[p]->L;
	  hmm_mu[p]  = log(tmp_K * ((double) db_size)) / hmm_evd[p]->lambda;
	  tmp_K      = exp(cm_evd[p]->mu * cm_evd[p]->lambda) / cm_evd[p]->L;
	  cm_mu[p]   = log(tmp_K  * ((double) db_size)) / cm_evd[p]->lambda;
	  /* Now determine bit score */
	  cm_minbitsc[p] = cm_mu[p] - (log(cm_ecutoff) / cm_evd[p]->lambda);
	  if(use_cm_cutoff)
	    printf("E: %f p: %d %d--%d bit score: %f\n", cm_ecutoff, p, 
		   cmstats->ps[p], cmstats->pe[p], cm_minbitsc[p]);
	}
      fflush(stdout);
      /* Strategy: 
       * Emit seqs from the CM and send them to workers to search with 
       * CM and with CP9 HMM and return scores of both searches.
       * Master collects these scores and determines if they are better 
       * than cm_ecutoff, and continues to do so until N seqs with E-values 
       * better than cm_ecutoff have been collected.
       *
       * If do_fastfil, the worker skips the CM search and returns 0. bits
       * as CM score, which master replaces with parsetree score.
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
	  if (i < N)
	    {
	      EmitParsetree(cm, &tr, &seq, &dsq, &L);
	      gc_comp = get_gc_comp(seq, 1, L); /* should be i and j of best hit, but
						 * EVD construction is wrong too, it uses 
						 * GC_COMP of full sample seq, not of best hit
						 * within it. */
	      fflush(stdout);
	      p    = cmstats->gc2p[gc_comp];
	      free(seq);
	      if(!do_fastfil)
		FreeParsetree(tr);
	    }
	  else have_work = FALSE;
	  /* If we have work but no free workers, or we have no work but workers
	   * are still working, then wait for a result to return from any worker.
	   */
	  if ( (have_work && nproc_working == nproc-1) || (! have_work && nproc_working > 0))
	    {
	      MPI_Recv(scores,  2, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mstatus);
	      cur_cm_sc  = scores[0];
	      cur_hmm_sc = scores[1];
	      wi = mstatus.MPI_SOURCE;
	      /* could change this to keep ALL CP9 and CM bit scores */
	      if(do_fastfil)
		{
		  cur_cm_sc = ParsetreeScore(cm, trlist[wi], dsqlist[wi], FALSE); 
		  FreeParsetree(tr);
		  trlist[wi] = NULL;
		}
	      /* careful not go let i get bigger than N, which can happen due to send/recv lag */
	      if(cur_cm_sc >= cm_minbitsc[plist[wi]] && i < N)
		{
		  hmm_eval[i] = RJK_ExtremeValueE(cur_hmm_sc, hmm_mu[plist[wi]], hmm_evd[plist[wi]]->lambda);
		  hmm_sc[i]   = cur_hmm_sc;
		  printf("i: %d sc: %f E: %f P: %f\n", hmm_sc[i], hmm_eval[i], esl_gumbel_surv((double) hmm_sc[i], hmm_mu[plist[wi]], hmm_evd[plist[wi]]->lambda));
		  fflush(stdout);
		  i++;
		}		  
	      nattempts++;
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
	      plist[wi]   = p;
	      if(do_fastfil)
		trlist[wi] = tr;
	      wi++;
	      nproc_working++;
	    }
	  else 
	    dsq_MPISend(NULL, -1, wi);	
	}
      /* Sort the HMM E-values with quicksort */
      esl_vec_FSortIncreasing(hmm_eval, N);
      esl_vec_FSortDecreasing(hmm_sc, N);
      
      printf("Dumping HMM scores fthr_mode: %d hmm_mode: %d\n", fthr_mode, hmm_evd_mode);
      for (i = 0; i < N; i++)
	printf("i: %4d hmm sc: %10.4f hmm E: %10.4f\n", i, hmm_sc[i], hmm_eval[i]);
      printf("\n\nnattempts: %d\n", nattempts);
      fflush(stdout);
      
      /* Clean up and exit */
      return_eval = hmm_eval[(int) ((fraction) * (float) N)];
      if(return_eval > ((float) db_size)) /* E-val > db_size is useless */
	return_eval = (float) db_size;
      free(scores);
      free(dsqlist);
      free(plist);
      free(trlist);
      free(hmm_eval);
      free(hmm_sc);
      free(cm_minbitsc);
      free(cm_mu);
      free(hmm_mu);
      fflush(stdout);
      /* Return threshold */
      return return_eval;
    }    

 ERROR:
  return;
}
#endif

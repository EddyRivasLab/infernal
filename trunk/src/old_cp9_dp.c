/* old_cp9_dp.c: formerly CP9_scan.c 
 * EPN, Wed Dec  5 10:42:01 2007
 * 
 * Scanning algorithms for CM Plan 9 HMMs.  These algorithms align
 * subsequences of the target sequence to the model (e.g. glocal or
 * local alignment) Global alignment algorithms are in hmmband.c.
 *
 * These functions are old, used in version 0.81 of Infernal.
 * They are currently ONLY compiled for the cp9_dp.c benchmark
 * driver, in which case they are used in comparison with their
 * new version 1.0 counterparts in cp9_dp.c.
 *
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
 *################################################################
 * 
 */
#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dirichlet.h"
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_fileparser.h"
#include "esl_gamma.h"
#include "esl_getopts.h"
#include "esl_gev.h"
#include "esl_gumbel.h"
#include "esl_histogram.h"
#include "esl_hyperexp.h"
#include "esl_keyhash.h"
#include "esl_minimizer.h"
#include "esl_mixgev.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msaweight.h"
#include "esl_normal.h"
#include "esl_paml.h"
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_regexp.h"
#include "esl_rootfinder.h"
#include "esl_scorematrix.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_stopwatch.h"
#include "esl_stretchexp.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_weibull.h"
#include "esl_wuss.h"

#include "funcs.h"
#include "old_funcs.h"
#include "structs.h"

/***********************************************************************
 * Function: CP9Viterbi()
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
CP9Viterbi(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
	   int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
	   int be_efficient, CP9_MX **ret_mx, CP9trace_t **ret_tr)
{
  int          status;
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /* j-W: position in the subsequence                             */
  int          ip;          /* i': relative position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  int          L;           /* j0-i0+1: subsequence length                                  */
  CP9_MX *mx;       /* the CP9 DP matrix                                            */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cm->cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int         *sc;          /* prob (seq from j0..jp | HMM) [0..jp..cm->cp9->M]             */
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
  /*debug_print_cp9_params(stdout, cm->cp9, TRUE);*/

  /*printf("in CP9Viterbi() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    cm_Fail("ERROR in CP9Viterbi, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    cm_Fail("ERROR in CP9Viterbi, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    cm_Fail("ERROR in CP9Viterbi, passing in results data structure, but not in scanning mode.\n");
  if(dsq == NULL)
    cm_Fail("ERROR in CP9Viterbi, dsq is NULL.");
    
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
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0.;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 1;
  else             nrows = L;
  mx = CreateCP9Matrix(nrows, cm->cp9->M);
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;
  elmx = mx->elmx;
  erow = mx->erow;

  /* sc will hold P(seq up to j | Model) in int log odds form */
  ESL_ALLOC(sc, sizeof(int) * (j0-i0+2));
			
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
  for (k = 1; k <= cm->cp9->M; k++)
    {
      mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
      dmx[0][k]  = -INFTY;
      if((tmp_sc = mmx[0][k-1] + cm->cp9->tsc[CTMD][k-1]) > dmx[0][k])
	dmx[0][k] = tmp_sc;
      if((tmp_sc = imx[0][k-1] + cm->cp9->tsc[CTID][k-1]) > dmx[0][k])
	dmx[0][k] = tmp_sc;
      if((tmp_sc = dmx[0][k-1] + cm->cp9->tsc[CTDD][k-1]) > dmx[0][k])
	dmx[0][k] = tmp_sc;
      /*printf("mmx[jp:%d][%d]: %d %f\n", 0, k, mmx[0][k], Score2Prob(mmx[0][k], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", 0, k, imx[0][k], Score2Prob(imx[0][k], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));
	printf("elmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));*/
    }
  /* We can do a full parse through all delete states. */
  erow[0] = dmx[0][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]; 
  sc[0] = erow[0];
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
      /* The 1 difference between a Viterbi scanner and the 
       * regular Viterbi. In non-scanner parse must begin in B at
       * position 0 (i0-1), in scanner we can start at any position 
       * in the seq. */
      if(do_scan)
	mmx[cur][0] = 0;
      else
	mmx[cur][0] = -INFTY;

      dmx[cur][0] = -INFTY;  /*D_0 is non-existent*/
      elmx[cur][0]= -INFTY;  /*no EL state for node 0 */

      imx[cur][0] = -INFTY; /* initialization */
      if((tmp_sc = mmx[prv][0] + cm->cp9->tsc[CTMI][0]) > imx[cur][0])
	imx[cur][0] = tmp_sc;
      if((tmp_sc = imx[prv][0] + cm->cp9->tsc[CTII][0]) > imx[cur][0])
	imx[cur][0] = tmp_sc;
      if((tmp_sc = dmx[prv][0] + cm->cp9->tsc[CTDI][0]) > imx[cur][0])
	imx[cur][0] = tmp_sc;
      if(imx[cur][0] != -INFTY)
	imx[cur][0] += cm->cp9->isc[dsq[j]][0];
      /*printf("mmx[jp:%d][%d]: %d %f\n", jp, 0, mmx[cur][0], Score2Prob(mmx[cur][0], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", jp, 0, imx[cur][0], Score2Prob(imx[cur][0], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", jp, 0, dmx[cur][0], Score2Prob(dmx[cur][0], 1.));*/

      for (k = 1; k <= cm->cp9->M; k++)
	{
	  /*match state*/
	  mmx[cur][k] = -INFTY;
	  if((tmp_sc = mmx[prv][k-1] + cm->cp9->tsc[CTMM][k-1]) > mmx[cur][k])
	    mmx[cur][k] = tmp_sc;
	  if((tmp_sc = imx[prv][k-1] + cm->cp9->tsc[CTIM][k-1]) > mmx[cur][k])
	    mmx[cur][k] = tmp_sc;
	  if((tmp_sc = mmx[prv][0] + cm->cp9->bsc[k]) > mmx[cur][k])
	    mmx[cur][k] = tmp_sc;
	  if((tmp_sc = dmx[prv][k-1] + cm->cp9->tsc[CTDM][k-1]) > mmx[cur][k])
	    mmx[cur][k] = tmp_sc;

	  /* check possibility we came from an EL, if they're valid */
	  if(cm->cp9->flags & CPLAN9_EL) 
	    for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	      /* transition penalty to EL incurred when EL was entered */
	      if((tmp_sc = elmx[prv][cm->cp9->el_from_idx[k][c]]) > mmx[cur][k])
		mmx[cur][k] = tmp_sc;
	  
	  if(mmx[cur][k] != -INFTY)
	    mmx[cur][k] += cm->cp9->msc[dsq[j]][k];

	  /*insert state*/
	  imx[cur][k] = -INFTY;
	  if((tmp_sc = mmx[prv][k] + cm->cp9->tsc[CTMI][k]) > imx[cur][k])
	    imx[cur][k] = tmp_sc;
	  if((tmp_sc = imx[prv][k] + cm->cp9->tsc[CTII][k]) > imx[cur][k])
	    imx[cur][k] = tmp_sc;
	  if((tmp_sc = dmx[prv][k] + cm->cp9->tsc[CTDI][k]) > imx[cur][k])
	    imx[cur][k] = tmp_sc;
	  if(imx[cur][k] != -INFTY)
	    imx[cur][k] += cm->cp9->isc[dsq[j]][k];
	  else 
	    imx[cur][k] = -INFTY;

	  /*delete state*/
	  dmx[cur][k] = -INFTY;
	  if((tmp_sc = mmx[cur][k-1] + cm->cp9->tsc[CTMD][k-1]) > dmx[cur][k])
	    dmx[cur][k] = tmp_sc;
	  if((tmp_sc = imx[cur][k-1] + cm->cp9->tsc[CTID][k-1]) > dmx[cur][k])
	    dmx[cur][k] = tmp_sc;
	  if((tmp_sc = dmx[cur][k-1] + cm->cp9->tsc[CTDD][k-1]) > dmx[cur][k])
	    dmx[cur][k] = tmp_sc;

	  elmx[cur][k] = -INFTY;
	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								    HMM nodes that map to right half of a MATP_MP) */
	    {
	      if((tmp_sc = mmx[cur][k] + cm->cp9->tsc[CTMEL][k]) > elmx[cur][k]) /* transitioned from cur node's match state */

		elmx[cur][k] = tmp_sc;
	      if((tmp_sc = elmx[prv][k] + cm->cp9->el_selfsc) > elmx[cur][k]) /* transitioned from cur node's EL state emitted ip on transition */
		elmx[cur][k] = tmp_sc;
	    }
	  else elmx[cur][k] = -INFTY;
	  /*printf("mmx[jp:%d][%d]: %d %f\n", jp, k, mmx[cur][k], Score2Prob(mmx[cur][k], 1.));
	    printf("imx[jp:%d][%d]: %d %f\n", jp, k, imx[cur][k], Score2Prob(imx[cur][k], 1.));
	    printf("dmx[jp:%d][%d]: %d %f\n", jp, k, dmx[cur][k], Score2Prob(dmx[cur][k], 1.));*/
	}
      /* determine erow[cur] == sc[jp], the int score of all possible parses ending at the current
       * position (j) of the target sequence. */
      erow[cur] = -INFTY;

      for (k = 1; k <= cm->cp9->M; k++)
	if ((tmp_sc = mmx[cur][k] + cm->cp9->esc[k]) > erow[cur])
	  erow[cur] = tmp_sc;
      if ((tmp_sc =  dmx[cur][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]) > erow[cur])
	erow[cur] = tmp_sc;
      /* transition from D_M -> end */
      if ((tmp_sc =  imx[cur][cm->cp9->M] + cm->cp9->tsc[CTIM][cm->cp9->M]) > erow[cur])
	erow[cur] = tmp_sc;
      /* transition from I_M -> end */
      if(cm->cp9->flags & CPLAN9_EL) /* no need to waste time otherwise */
	{
	  /* check if we came from an EL */
	  for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M+1]; c++) /* el_from_ct[cm->cp9->M+1] holds # ELs that can go to END */
	    if((tmp_sc = elmx[cur][cm->cp9->el_from_idx[cm->cp9->M+1][c]]) > erow[cur])
	      erow[cur] = tmp_sc;
	  /* transition penalty to EL incurred when EL was entered */
	}
      sc[jp] = erow[cur];
      fsc = Scorify(erow[cur]);
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
	  curr_sc = gamma[ip-1] + fsc;
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
      return_sc = Scorify(sc[(j0-i0+1)]); /* L = j0-i0+1 */
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
  if(ret_sc != NULL) *ret_sc = sc;
  else free(sc);
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCP9Matrix(mx);
  /*printf("Viterbi return_sc: %f\n", return_sc);*/

  return return_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
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
 *           results   - search_results_t to add to; if NULL, don't keep results
 *           do_scan   - TRUE if we're scanning, HMM can start to emit anywhere i0..j0,
 *                       FALSE if we're not, HMM must start emitting at i0, end emitting at j0
 *           doing_align  - TRUE if reason we've called this function is to help get posteriors
 *                          for CM alignment, in which case we can skip the traceback. 
 *           be_efficient- TRUE to keep only 2 rows of DP matrix in memory, FALSE keep whole thing
 *           ret_mx    - RETURN: CP9 Forward DP matrix, NULL if not wanted
 *
 * Returns:  if(!do_scan) log P(S|M)/P(S|R), as a bit score
 *           else         max log P(S|M)/P(S|R), for argmax subseq S of input seq i0..j0,
 */
float
CP9Forward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
	   int *ret_bestpos, search_results_t *results, int do_scan, int doing_align,
	   int be_efficient, CP9_MX **ret_mx)
{
  int          status;
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /* j-W: position in the subsequence                             */
  int          ip;          /* i': relative position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  int          L;           /* j0-i0+1: subsequence length                                  */
  CP9_MX *mx;       /* the CP9 DP matrix                                            */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cm->cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int         *sc;          /* prob (seq from j0..jp | HMM) [0..jp..L]                      */
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
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */
  int          c;           /* counter for EL states */
  /*debug_print_cp9_params(stdout, cm->cp9, TRUE);*/

  /*printf("in CP9Forward() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    cm_Fail("ERROR in CP9Forward, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    cm_Fail("ERROR in CP9Forward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    cm_Fail("ERROR in CP9Forward, passing in results data structure, but not in scanning mode.\n");
  if(dsq == NULL)
    cm_Fail("ERROR in CP9Forward, dsq is NULL.");
    
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
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0.;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 1;
  else             nrows = L;
  mx = CreateCP9Matrix(nrows, cm->cp9->M);
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;
  elmx = mx->elmx;
  erow = mx->erow;

  /* sc will hold P(seq up to j | Model) in int log odds form */
  ESL_ALLOC(sc, sizeof(int) * (j0-i0+2));
			
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
  for (k = 1; k <= cm->cp9->M; k++)
    {
      mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
      dmx[0][k] = ILogsum(ILogsum(mmx[0][k-1] + cm->cp9->tsc[CTMD][k-1],
				  imx[0][k-1] + cm->cp9->tsc[CTID][k-1]),
			  dmx[0][k-1] + cm->cp9->tsc[CTDD][k-1]);
      /*printf("mmx[jp:%d][%d]: %d %f\n", 0, k, mmx[0][k], Score2Prob(mmx[0][k], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", 0, k, imx[0][k], Score2Prob(imx[0][k], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));
	printf("elmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));*/
    }
  /* We can do a full parse through all delete states. */
  erow[0] = dmx[0][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]; 
  sc[0] = erow[0];
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
      elmx[cur][0]= -INFTY;  /*no EL state for node 0 */
      imx[cur][0]  = ILogsum(ILogsum(mmx[prv][0] + cm->cp9->tsc[CTMI][0],
				     imx[prv][0] + cm->cp9->tsc[CTII][0]),
			     dmx[prv][0] + cm->cp9->tsc[CTDI][0]);
      imx[cur][0] += cm->cp9->isc[dsq[j]][0];
      /*printf("mmx[jp:%d][%d]: %d %f\n", jp, 0, mmx[cur][0], Score2Prob(mmx[cur][0], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", jp, 0, imx[cur][0], Score2Prob(imx[cur][0], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", jp, 0, dmx[cur][0], Score2Prob(dmx[cur][0], 1.));*/

      for (k = 1; k <= cm->cp9->M; k++)
	{
	  mmx[cur][k]  = ILogsum(ILogsum(mmx[prv][k-1] + cm->cp9->tsc[CTMM][k-1],
				       imx[prv][k-1] + cm->cp9->tsc[CTIM][k-1]),
				 ILogsum(mmx[prv][0] + cm->cp9->bsc[k],
					 dmx[prv][k-1] + cm->cp9->tsc[CTDM][k-1]));
	  /* Add contribution of ELs if they're valid */
	  if(cm->cp9->flags & CPLAN9_EL) /* no need to waste time */
	    for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	      /* transition penalty to EL incurred when EL was entered */
	      mmx[cur][k] = ILogsum(mmx[cur][k], (elmx[prv][cm->cp9->el_from_idx[k][c]]));
	  mmx[cur][k] += cm->cp9->msc[dsq[j]][k];

	  dmx[cur][k]  = ILogsum(ILogsum(mmx[cur][k-1] + cm->cp9->tsc[CTMD][k-1],
					imx[cur][k-1] + cm->cp9->tsc[CTID][k-1]),
				dmx[cur][k-1] + cm->cp9->tsc[CTDD][k-1]);
	  
	  imx[cur][k]  = ILogsum(ILogsum(mmx[prv][k] + cm->cp9->tsc[CTMI][k],
				       imx[prv][k] + cm->cp9->tsc[CTII][k]),
			       dmx[prv][k] + cm->cp9->tsc[CTDI][k]);
	  imx[cur][k] += cm->cp9->isc[dsq[j]][k];

	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
							    HMM nodes that map to right half of a MATP_MP) */
	    {
	      elmx[cur][k] = ILogsum(mmx[cur][k] + cm->cp9->tsc[CTMEL][k], 
				     /* transitioned from cur node's match state */
				     elmx[prv][k] + cm->cp9->el_selfsc);
	      /* transitioned from cur node's EL state emitted ip on transition */
	    }
	  else elmx[cur][k] = -INFTY;
	  /*printf("mmx[jp:%d][%d]: %d %f\n", jp, k, mmx[cur][k], Score2Prob(mmx[cur][k], 1.));
	    printf("imx[jp:%d][%d]: %d %f\n", jp, k, imx[cur][k], Score2Prob(imx[cur][k], 1.));
	    printf("dmx[jp:%d][%d]: %d %f\n", jp, k, dmx[cur][k], Score2Prob(dmx[cur][k], 1.));*/
	}
      /* determine erow[cur] == sc[jp], the int score of all possible parses ending at the current
       * position (j) of the target sequence. */
      erow[cur] = -INFTY;
      for (k = 1; k <= cm->cp9->M; k++)
	erow[cur] = ILogsum(erow[cur], mmx[cur][k] + cm->cp9->esc[k]);
      /* 04.17.07 Arent' I double counting here! (at least for scanner?) */
      erow[cur] = ILogsum(erow[cur], dmx[cur][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]); 
      erow[cur] = ILogsum(erow[cur], imx[cur][cm->cp9->M] + cm->cp9->tsc[CTIM][cm->cp9->M]); 
      if(cm->cp9->flags & CPLAN9_EL) /* no need to waste time otherwise */
	{
	  /* check if we came from an EL */
	  for(c = 0; c < cm->cp9->el_from_ct[cm->cp9->M+1]; c++) /* el_from_ct[cm->cp9->M+1] holds # ELs that can go to END */
	    erow[cur] = ILogsum(erow[cur], elmx[cur][cm->cp9->el_from_idx[cm->cp9->M+1][c]]);
	    /* transition penalty to EL incurred when EL was entered */
	}
      sc[jp] = erow[cur];
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
	  curr_sc = gamma[ip-1] + fsc;
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
		if(results != NULL) /* report the hit */
		{
		  /*printf("FWD reporting hit: i: %d j: %d sc: %f\n", gback[jp], j, savesc[jp]);*/
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
  else                FreeCP9Matrix(mx);
  /*printf("Forward return_sc: %f\n", return_sc);*/

  return return_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
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
 *           be_efficient- TRUE to keep only 2 rows of DP matrix in memory, FALSE keep whole thing
 *           ret_mx    - RETURN: CP9 Forward DP matrix, NULL if not wanted
 *
 * Returns:  if(!do_scan) log P(S|M)/P(S|R), as a bit score, this is B->M[0][0];
 *           else         max log P(S|M)/P(S|R), for argmax subseq S of input seq i0..j0,
 *                            this is max_jp B->M[jp][0]
 */
float
CP9Backward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
	    int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
	    int be_efficient, CP9_MX **ret_mx)
{
  int          status;
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /*     j-W: position in the subsequence                         */
  int          ip;          /* i': relative position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  int          L;           /* j0-i0+1: subsequence length                                  */
  CP9_MX *mx;       /* the CP9 DP matrix                                            */
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
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */
  int          c;           /* counter for EL states */

  /*printf("in CP9Backward() i0: %d j0: %d do_scan: %d \n", i0, j0, do_scan);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    cm_Fail("ERROR in CP9Backward, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    cm_Fail("ERROR in CP9Backward, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    cm_Fail("ERROR in CP9Backward, passing in results data structure, but not in scanning mode.a\n");
  if(dsq == NULL)
    cm_Fail("ERROR in CP9Backward, dsq is NULL.");
    
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
  if(be_efficient) nrows = 1;
  else             nrows = L; 
  mx = CreateCP9Matrix(nrows, cm->cp9->M);
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;
  elmx = mx->elmx;
  erow = mx->erow;

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
      /* init EL mx to -INFTY */
      for (k = 1; k <= cm->cp9->M; k++)
	elmx[cur][k] = -INFTY;
      
      if(ip > 0)
	{
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
	}
      else /* ip == 0 */
	{
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
	    elmx[cur][cm->cp9->el_from_idx[cm->cp9->M][c]] = ILogsum(elmx[cur][cm->cp9->el_from_idx[cm->cp9->M][c]], mmx[prv][cm->cp9->M]);
	}

      
      /* A main difference between a Backward scanner and 
       * regular Backward: a scanner can end at the END 
       * state at any position, regular can only end at
       * the final position j0. */
      if(do_scan)
	{	
	  if(ip > 0)
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
	      /*******************************************************************
	       * 3 Handle EL, looking at EL_k->M_k for all valid k and EL_k->EL_k
	       * we're going backwards so we have to work out of order
	       * we could get around this by storing the nodes each EL goes TO
	       * in an el_to_ct[] vector. */
	      if(cm->cp9->flags & CPLAN9_EL)
		{
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

	    }
	  else
	    {
	      mmx[cur][k] = -INFTY; /* need seq to get here, unless we come from E in a scanner (below) */
	      imx[cur][k] = -INFTY; /* need seq to get here */
	      elmx[cur][k]= -INFTY;  /* first emitted res can't be from an EL, need to see >= 1 matches */
	    }
	  if(do_scan && ip > 0) /* add possibility of ending at this position from this state */
	    {
	      mmx[cur][k] = 
		ILogsum(mmx[cur][k], 
			(cm->cp9->esc[k] +                    /* M_k<-E + (only in scanner)     */ 
			 0));                                 /* all parses end in E, 2^0 = 1.0;*/
	      /* DO NOT add contribution of emitting i from M, it's been added above */
	      /* No EL contribution here b/c we'd be looking for M_k<-EL_k<-E, but EL_k<-E is impossible 
	       * for k != cm->cp9->M; */
	      /* HERE */
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
	  imx[cur][0] += cm->cp9->isc[dsq[i]][k];
	}
      else /* ip == 0 */
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
	  curr_sc = gamma[jp+1-1] + fsc; /* *off-by-one* */
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
  else                FreeCP9Matrix(mx);
  /*printf("Backward return_sc: %f\n", return_sc);*/
  return return_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
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
CP9ScanPosterior(ESL_DSQ *dsq, int i0, int j0,
		     CP9_t *hmm,
		     CP9_MX *fmx,
		     CP9_MX *bmx,
		     CP9_MX *mx)
{
  int i;
  int ip;
  int k;
  int fb_sum; /* tmp value, the probability that the current residue (i) was
	       * visited in any parse */
  /* contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in CP9ScanPosterior(), dsq is NULL.");

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
      mx->imx[ip][0] = fmx->imx[ip][0] + bmx->imx[ip][0] - hmm->isc[dsq[i]][0] - fb_sum;
      /*hmm->isc[dsq[i]][0] will have been counted in both fmx->imx and bmx->imx*/
      mx->dmx[ip][0] = -INFTY; /*D_0 does not exist*/
      for (k = 1; k <= hmm->M; k++) 
	{
	  mx->mmx[ip][k] = fmx->mmx[ip][k] + bmx->mmx[ip][k] - hmm->msc[dsq[i]][k] - fb_sum;
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  mx->imx[ip][k] = fmx->imx[ip][k] + bmx->imx[ip][k] - hmm->isc[dsq[i]][k] - fb_sum;
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


/***********************************************************************
 * Function: CP9ForwardAlign
 * based on  P7Forward() <-- this function's comments below  
 *           from HMMER 2.4devl core_algorithms.c
 *
 * Purpose:  The Forward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Note:     Differs from CP9_scan.c:CP9Forward() in that this function
 *           does not work in 'scan' mode. CP9_scan.c:CP9Forward() can
 *           do scanning OR alignment, but is somewhat more complex than
 *           this function.
 *
 * Args:     dsq    - sequence in digitized form
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           hmm    - the model
 *           mx     - the dp matrix
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
CP9ForwardAlign(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx)
{
  int status;
  int **mmx;
  int **imx;
  int **dmx;
  int **elmx;
  int  *erow;
  int   i,k;
  int   sc;
  int   L;		/* subsequence length */
  int   ip;		/* i': relative position in the subsequence  */

  if(dsq == NULL) cm_Fail("In CP9ForwardAlign() dsq is NULL.");
  L  = j0-i0+1;		/* the length of the subsequence */

  /* Grow the DP matrix to 0..L rows */
  if((status = GrowCP9Matrix(mx, NULL, L, hmm->M, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) cm_Fail("CP9ForwardAlign() error growing dp mx.");

  /* Initialization of the zero row.
   */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= hmm->M; k++)
    dmx[0][k] = ILogsum(ILogsum(mmx[0][k-1] + hmm->tsc[CTMD][k-1],
				imx[0][k-1] + hmm->tsc[CTID][k-1]),
			dmx[0][k-1] + hmm->tsc[CTDD][k-1]);
  
  erow[0] = dmx[0][hmm->M] + hmm->tsc[CTDM][hmm->M]; 
  /*printf("i: %d score(erow[ip]): %f erow[ip]: %d\n", 0, Scorify(erow[0]), erow[0]);*/
  /* Recursion. Done as a pull.
   */
  for (ip = 1; ip <= L; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */

      mmx[ip][0] = dmx[ip][0] = -INFTY;  /*M_0 (B) and D_0 (non-existent)
					 don't emit.
				       */
      imx[ip][0]  = ILogsum(ILogsum(mmx[ip-1][0] + hmm->tsc[CTMI][0],
				    imx[ip-1][0] + hmm->tsc[CTII][0]),
			    dmx[ip-1][0] + hmm->tsc[CTDI][0]);
      imx[ip][0] += hmm->isc[dsq[i]][0];
      for (k = 1; k <= hmm->M; k++)
	{
	  mmx[ip][k]  = ILogsum(ILogsum(mmx[ip-1][k-1] + hmm->tsc[CTMM][k-1],
				       imx[ip-1][k-1] + hmm->tsc[CTIM][k-1]),
			       ILogsum(mmx[ip-1][0] + hmm->bsc[k],
				       dmx[ip-1][k-1] + hmm->tsc[CTDM][k-1]));
	  mmx[ip][k] += hmm->msc[dsq[i]][k];

	  dmx[ip][k]  = ILogsum(ILogsum(mmx[ip][k-1] + hmm->tsc[CTMD][k-1],
				       imx[ip][k-1] + hmm->tsc[CTID][k-1]),
			       dmx[ip][k-1] + hmm->tsc[CTDD][k-1]);

	  imx[ip][k]  = ILogsum(ILogsum(mmx[ip-1][k] + hmm->tsc[CTMI][k],
				       imx[ip-1][k] + hmm->tsc[CTII][k]),
			       dmx[ip-1][k] + hmm->tsc[CTDI][k]);
	  imx[ip][k] += hmm->isc[dsq[i]][k];
	}
      erow[ip] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	erow[ip] = ILogsum(erow[ip], mmx[ip][k] + hmm->esc[k]);
      erow[ip] = ILogsum(erow[ip], dmx[ip][hmm->M] + hmm->tsc[CTDM][hmm->M]); 
      erow[ip] = ILogsum(erow[ip], imx[ip][hmm->M] + hmm->tsc[CTIM][hmm->M]); 
		       /* transition from D_M -> end */
      /*printf("i: %d score(erow[ip]): %f erow[ip]: %d\n", i, Scorify(erow[ip]), erow[ip]);*/
    }		
  sc = erow[L];
  /*printf("F erow[%d]: %d\n", i, erow[L]);*/

  return Scorify(sc);		/* the total Forward score. */
}

/* Function: CP9ViterbiAlign()
 * based on  P7Viterbi() <-- this function's comments below  
 *           from HMMER 2.4devl core_algorithms.c
 *
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's. 
 *           
 * Args:     dsq    - sequence in digitized form
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           hmm    - the model
 *           mx     - reused DP matrix
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
CP9ViterbiAlign(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx,
	      CP9trace_t **ret_tr)
{
  int status;
  CP9trace_t  *tr;
  int **mmx;
  int **imx;
  int **dmx;
  int **elmx;
  int  *erow;
  int   i,k,c;
  int   sc;
  int   W;		/* subsequence length */
  int   ip;		/* i': relative position in the subsequence  */

  if(dsq == NULL)  cm_Fail("in CP9ViterbiAlign(), dsq is NULL.");

  W  = j0-i0+1;		/* the length of the subsequence */

  /* Grow the DP matrix to 0..W rows */
  if((status = GrowCP9Matrix(mx, NULL, W, hmm->M, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) cm_Fail("CP9ViterbiAlign() error growing dp mx.");
  /* Initialization of the zero row.
   */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  elmx[0][0]= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY; 

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= hmm->M; k++)
    {
      dmx[0][k]  = -INFTY;
      if((sc = mmx[0][k-1] + hmm->tsc[CTMD][k-1]) > dmx[0][k])
	dmx[0][k] = sc;
      if((sc = dmx[0][k-1] + hmm->tsc[CTDD][k-1]) > dmx[0][k])
	dmx[0][k] = sc;
    }

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   */
  for (ip = 1; ip <= W; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      mmx[ip][0] = dmx[ip][0] = -INFTY;  /* M_0 (B) and D_0 (non-existent)
					  * don't emit. */
      elmx[ip][0] = -INFTY; /* there is no EL state for node 0 */
      imx[ip][0] = -INFTY;
      if((sc = mmx[ip-1][0] + hmm->tsc[CTMI][0]) > imx[ip][0])
	imx[ip][0] = sc;
      if((sc = imx[ip-1][0] + hmm->tsc[CTII][0]) > imx[ip][0])
	imx[ip][0] = sc;
      if((sc = dmx[ip-1][0] + hmm->tsc[CTDI][0]) > imx[ip][0])
	imx[ip][0] = sc;
      if(imx[ip][0] != -INFTY)
	imx[ip][0] += hmm->isc[dsq[i]][0];
      else 
	imx[ip][0] = -INFTY;

      for (k = 1; k <= hmm->M; k++)
	{
	  /*match state*/
	  mmx[ip][k] = -INFTY;
	  if((sc = mmx[ip-1][k-1] + hmm->tsc[CTMM][k-1]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  if((sc = imx[ip-1][k-1] + hmm->tsc[CTIM][k-1]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  if((sc = mmx[ip-1][0] + hmm->bsc[k]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  if((sc = dmx[ip-1][k-1] + hmm->tsc[CTDM][k-1]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  /* Check if we came from an EL state */
	  if(hmm->flags & CPLAN9_EL) /* no need to waste time */
	    {
	      for(c = 0; c < hmm->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
		{
		  /* transition penalty to EL incurred when EL was entered */
		  if((sc = elmx[ip-1][hmm->el_from_idx[k][c]]) > mmx[ip][k])
		    mmx[ip][k] = sc;
		}
	    }
	  if(mmx[ip][k] != -INFTY)
	    mmx[ip][k] += hmm->msc[dsq[i]][k];
	  else 
	    mmx[ip][k] = -INFTY;

	  /*insert state*/
	  imx[ip][k] = -INFTY;
	  if((sc = mmx[ip-1][k] + hmm->tsc[CTMI][k]) > imx[ip][k])
	    imx[ip][k] = sc;
	  if((sc = imx[ip-1][k] + hmm->tsc[CTII][k]) > imx[ip][k])
	    imx[ip][k] = sc;
	  if((sc = dmx[ip-1][k] + hmm->tsc[CTDI][k]) > imx[ip][k])
	    imx[ip][k] = sc;
	  if(imx[ip][k] != -INFTY)
	    imx[ip][k] += hmm->isc[dsq[i]][k];
	  else 
	    imx[ip][k] = -INFTY;

	  /*delete state*/
	  dmx[ip][k] = -INFTY;
	  if((sc = mmx[ip][k-1] + hmm->tsc[CTMD][k-1]) > dmx[ip][k])
	    dmx[ip][k] = sc;
	  if((sc = imx[ip][k-1] + hmm->tsc[CTID][k-1]) > dmx[ip][k])
	    dmx[ip][k] = sc;
	  if((sc = dmx[ip][k-1] + hmm->tsc[CTDD][k-1]) > dmx[ip][k])
	    dmx[ip][k] = sc;

	  /*EL (end-local) state*/
	  elmx[ip][k] = -INFTY;
	  if((hmm->flags & CPLAN9_EL) && hmm->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
							    HMM nodes that map to right half of a MATP_MP) */
	    {
	      if((sc = mmx[ip][k] + hmm->tsc[CTMEL][k]) > elmx[ip][k])
		elmx[ip][k] = sc; /* transitioned from cur node's match state */
	      if((sc = elmx[ip-1][k] + hmm->el_selfsc) > elmx[ip][k])
		elmx[ip][k] = sc; /* transitioned from cur node's EL state emitted ip on transition */
	    }
	}
      erow[ip] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	if ((sc = mmx[ip][k] + hmm->esc[k]) > erow[ip])
	  erow[ip] = sc;
      if ((sc =  dmx[ip][hmm->M] + hmm->tsc[CTDM][hmm->M]) > erow[ip]) /* transition from D_M -> end */
	erow[ip] = sc;
      /* check if we came from an EL */
      if(hmm->flags & CPLAN9_EL) /* no need to waste time */
	{
	  for(c = 0; c < hmm->el_from_ct[hmm->M+1]; c++) /* el_from_ct[hmm->M+1] holds # ELs that can go to END */
	    {
	      /* transition penalty to EL incurred when EL was entered */
	      if((sc = elmx[ip][hmm->el_from_idx[hmm->M+1][c]]) > erow[ip])
		erow[ip] = sc;
	    }
	}
    } 
  sc = erow[W];
  /* printf("returning sc: %d from CP9ViterbiAlign()\n", sc); */
  
  if (ret_tr != NULL) {
    CP9ViterbiTrace(hmm, dsq, i0, j0, mx, &tr);
    /* CP9PrintTrace(stdout, tr, hmm, dsq); */
    *ret_tr = tr;
  }
  return Scorify(sc);		/* the total Viterbi score. */
}

/*********************************************************************
 * Function: CP9BackwardAlign
 * based on  P7Backward() <-- this function's comments below
 *           from HMMER 2.4devl core_algorithms.c
 * 
 * Note:     Differs from CP9_scan.c:CP9Backward() in that this function
 *           does not work in 'scan' mode. CP9_scan.c:CP9Backward() can
 *           do scanning OR alignment, but is somewhat more complex than
 *           this function.
 * 
 * Purpose:  The Backward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           hmm    - the model
 *           mx     - the dp matrix
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
CP9BackwardAlign(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx)
{
  int status;
  int **mmx;
  int **imx;
  int **dmx;
  int **elmx;
  int  *erow;
  int   i,k;
  int   sc;

  int   W;		/* subsequence length */
  int   ip;		/* i': relative position in the subsequence  */

  if(dsq == NULL) cm_Fail("in CP9BackwardAlign(), dsq is NULL.");

  W  = j0-i0+1;		/* the length of the subsequence */

  /* Grow the DP matrix to 0..W rows */
  if((status = GrowCP9Matrix(mx, NULL, W, hmm->M, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) cm_Fail("CP9ForwardAlign() error growing dp mx.");

  /* Initialization of the W row.
   */
  i = j0;

  erow[W] = 0; /*have to end in E*/

  mmx[W][hmm->M] = erow[W] + hmm->esc[hmm->M]; /* M<-E ...                   */
  mmx[W][hmm->M] += hmm->msc[dsq[i]][hmm->M]; /* ... + emitted match symbol */
  /* can't come from I_M b/c we've emitted a single residue, L from M_M */
  imx[W][hmm->M] = erow[W] + hmm->tsc[CTIM][hmm->M];   /* I_M(C)<-E ... */
  imx[W][hmm->M] += hmm->isc[dsq[i]][hmm->M];           /* ... + emitted match symbol */
  dmx[W][hmm->M] = erow[W] + hmm->tsc[CTDM][hmm->M];    /* D_M<-E */
  for (k = hmm->M-1; k >= 1; k--)
    {
      mmx[W][k]  = hmm->esc[k] + erow[W];
      mmx[W][k]  = ILogsum(mmx[W][k], dmx[W][k+1] + hmm->tsc[CTMD][k]);
      mmx[W][k] += hmm->msc[dsq[i]][k];

      imx[W][k] = dmx[W][k+1] + hmm->tsc[CTID][k];
      imx[W][k] += hmm->isc[dsq[i]][k];

      dmx[W][k] = dmx[W][k+1] + hmm->tsc[CTDD][k];
    }
  
  mmx[W][0] = -INFTY; /*M_0 doesn't emit*/
  imx[W][0] = dmx[W][1] + hmm->tsc[CTID][0];
  imx[W][0] += hmm->isc[dsq[i]][hmm->M];    
  dmx[W][0] = -INFTY; /*D_0 doesn't exist*/

  /* Recursion. Done as a pull.
   */

  for (ip = W-1; ip >= 1; ip--) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from j0 down to i0 */
      erow[ip] = -INFTY;
      
      /* Now the main states. Note the boundary conditions at M.
       */
      mmx[ip][hmm->M] = imx[ip+1][hmm->M] + hmm->tsc[CTMI][hmm->M];
      mmx[ip][hmm->M] += hmm->msc[dsq[i]][hmm->M];
      imx[ip][hmm->M] = imx[ip+1][hmm->M] + hmm->tsc[CTII][hmm->M];
      imx[ip][hmm->M] += hmm->isc[dsq[i]][hmm->M];
      dmx[ip][hmm->M] = imx[ip+1][hmm->M] + hmm->tsc[CTDI][hmm->M];  /* * */
      for (k = hmm->M-1; k >= 1; k--)
	{
	  mmx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTMM][k],
				       imx[ip+1][k] + hmm->tsc[CTMI][k]),
				dmx[ip][k+1] + hmm->tsc[CTMD][k]);
	  mmx[ip][k] += hmm->msc[dsq[i]][k];
	  
	  imx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTIM][k],
				       imx[ip+1][k] + hmm->tsc[CTII][k]),
			       dmx[ip][k+1] + hmm->tsc[CTID][k]);
	  imx[ip][k] += hmm->isc[dsq[i]][k];
	  
	  dmx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTDM][k],
				       imx[ip+1][k] + hmm->tsc[CTDI][k]),
			       dmx[ip][k+1] + hmm->tsc[CTDD][k]);
	}

      imx[ip][0]  = ILogsum(ILogsum(mmx[ip+1][1] + hmm->tsc[CTIM][0],
				    imx[ip+1][0] + hmm->tsc[CTII][0]),
			    dmx[ip][1] + hmm->tsc[CTID][0]);
      imx[ip][0] += hmm->isc[dsq[i]][0];
      mmx[ip][0] = -INFTY;
      /*for (k = hmm->M-1; k >= 1; k--)*/ /*M_0 is the B state, it doesn't emit*/
      for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	mmx[ip][0] = ILogsum(mmx[ip][0], mmx[ip+1][k] + hmm->bsc[k]);
      /* 04.17.07 NOTE: ERROR! mmx[ip][k==0] should takes mmx[ip+0(not 1)][k] b/c
       * M_0 does not emit, and ip was accounted for in mmx[ip][k] */
      mmx[ip][0] = ILogsum(mmx[ip][0], imx[ip+1][k] + hmm->tsc[CTMI][k]);
      mmx[ip][0] = ILogsum(mmx[ip][0], dmx[ip][k+1] + hmm->tsc[CTMD][k]);

      dmx[ip][0] = -INFTY; /* D_0 does not exist */
    }
  /* case when ip = 0 */
  erow[0] = -INFTY;
  mmx[0][hmm->M] = -INFTY; /* need seq to get here */
  imx[0][hmm->M] = -INFTY; /* need seq to get here */
  dmx[0][hmm->M] = imx[1][hmm->M] + hmm->tsc[CTDI][hmm->M];  /* * */
  for (k = hmm->M-1; k >= 1; k--)
    {
      mmx[0][k] = -INFTY; /* need seq to get here */
      imx[0][k] = -INFTY; /* need seq to get here */
      dmx[0][k]  = ILogsum(ILogsum(mmx[1][k+1] + hmm->tsc[CTDM][k],
				   imx[1][k] + hmm->tsc[CTDI][k]),
			   dmx[0][k+1] + hmm->tsc[CTDD][k]);
    }
  imx[0][0] = -INFTY; /*need seq to get here*/
  dmx[0][0] = -INFTY; /*D_0 doesn't exist*/

  mmx[0][0] = -INFTY;
  for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
    {
      mmx[0][0] = ILogsum(mmx[0][0], mmx[1][k] + hmm->bsc[k]);
    }
  mmx[0][0] = ILogsum(mmx[0][0], imx[1][0] + hmm->tsc[CTMI][0]);
  mmx[0][0] = ILogsum(mmx[0][0], dmx[0][1] + hmm->tsc[CTMD][0]);

  sc = mmx[0][0];
  /*printf("B final score: %f (i: %d)\n", Scorify(sc), sc);*/

  return Scorify(sc);		/* the total Backward score. */
}

/* old_cm_dpsearch.c
 * EPN, Tue Dec  4 17:41:14 2007 
 * 
 * note: 
 * formerly scaninside.c  EPN 05.08.06 (through v0.81, svn revision 2243)
 * absorbed scancyk.c when renamed to old_cm_dpsearch.c
 *
 * scancyk.c notes
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * scancyk.c 
 * SRE, Thu May  2 11:50:48 2002 [AA 3050 SFO->STL]
 * SVN $Id: scancyk.c 2243 2007-12-04 22:22:57Z nawrockie $
 * 
 * CYK alignment: multihit, local, database scanning mode.
 * [xref STL6 p47]
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * CM CYK/Inside alignment: multihit, local, database scanning mode.
 ***************************************************************** 
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 ***************************************************************** 
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

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
#include "old_funcs.h"		/* old external functions               */
#include "structs.h"

/**************************************************************
 * Function: InsideScan()
 * EPN 05.08.06
 * 
 * Based on CYKScan() by SRE
 *
 * Purpose:  Scan a (sub)sequence for matches to a covariance model.
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Log sums are performed (slowly) on floats by the LogSum2()
 *           function.
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; 1..L
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results    - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float 
InsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, 
	   float cutoff, search_results_t *results)

{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("in InsideScan, dsq is NULL\n");

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   *****************************************************************/

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*printf("in InsideScan i0: %d j0: %d\n", i0, j0);*/

  ESL_ALLOC(alpha, sizeof(float **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * (W+1));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], sizeof(float) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], sizeof(float) * (W+1));
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *****************************************************************/ 
  for (v = cm->M-1; v >= 0; v--)
    {
      alpha[v][0][0] = IMPOSSIBLE;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->endsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    alpha[v][0][0] = LogSum2(alpha[v][0][0], (alpha[y+yoffset][0][0] 
						     + cm->tsc[v][yoffset]));
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  if (alpha[v][0][0] < IMPOSSIBLE) alpha[v][0][0] = IMPOSSIBLE;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  for (d = 1; d <= W; d++)
    alpha[cm->M-1][0][d] = alpha[cm->M-1][1][d] = IMPOSSIBLE;

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0; 
  ESL_ALLOC(gback,  sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offeset j index in gamma* data structures */
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][jp][d] = LogSum2(alpha[v][jp][d], (alpha[y+yoffset][cur][d] 
							       + cm->tsc[v][yoffset]));
		  if (alpha[v][jp][d] < IMPOSSIBLE) alpha[v][jp][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = 2; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][prv][d-2] 
								 + cm->tsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v]))); 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][cur][d-1] 
								 + cm->tsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][prv][d-1] 
								 + cm->tsc[v][yoffset]));
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));

		  for (k = 0; k <= d; k++) /* k is length of right fragment */
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[w][jp][d-k] 
								   + alpha[y][cur][k]));
		    }
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      for (d = 1; d <= W && d <= gamma_j; d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    alpha[0][cur][d] = LogSum2(alpha[0][cur][d], (alpha[y+yoffset][cur][d] 
							 + cm->tsc[0][yoffset]));
	  if (cm->flags & CMH_LOCAL_BEGIN) {
	    /* EPN 05.08.06 Left local the same as it was in CYKScan(), I think this 
	     * is okay, we're saying that the best root of the alignment is the state that
	     * gives the highest Inside score. 
	     */
	    for (y = 1; y < cm->M; y++) {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	    }
	  }
	  if (alpha[0][cur][d] < IMPOSSIBLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	}

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; 
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = 1; d <= W && d <= gamma_j; d++) 
	{
	  i = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  sc = gamma[gamma_i-1] + alpha[0][cur][d];
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = alpha[0][cur][d]; 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	continue;
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  return best_score;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}


/* Function: CYKBandedScan()
 * Date:     SRE, Fri Sep 26 09:53:42 2003 [AA 886, returning from Salk Institute]
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm.
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *
 *           dmin,dmax set the bounds. Both are arrays, 0..v..cm->M. 
 *           The band for v is dmin[v]..dmax[v], inclusive; that is,
 *           dmin[v] is the minimum allowed d for state v (inclusive);
 *           dmax[v] is the maximum; that is, v can only account for
 *           subsequences of length >= dmin[v] and <= dmax[v].
 *
 * Args:     cm        - the covariance model
 *           dsq       - the digitized sequence
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results   - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 *
 * Note:     Dies (cm_Fail()) from memory allocation error, without
 *           cleanup.
 */
float
CYKBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
	      float cutoff, search_results_t *results)
{
  /* Contract check */
  if(j0 < i0)
    cm_Fail("ERROR in CYKBandedScan, i0: %d j0: %d\n", i0, j0);
  if(!(cm->flags & CMH_QDB))
    cm_Fail("ERROR in CYKBandedScan, QDBs invalid\n");
  if(dsq == NULL)
    cm_Fail("ERROR in CYKBandedScan, dsq is NULL.\n");
  if(dmin == NULL || dmax == NULL)
    cm_Fail("ERROR in CYKBandedScan, dmin and/or dmax are NULL.\n");

  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmax;               /* for B_st's, maximum k value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  int       curr_dmin;          /* temporary value for min d in for loops */
  int       curr_dmax;          /* temporary value for max d in for loops */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */

  /*PrintDPCellsSaved(cm, dmin, dmax, W);*/
  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*printf("in CYKBandedScan i0: %d j0: %d\n", i0, j0);*/

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   *****************************************************************/
  ESL_ALLOC(alpha, (sizeof(float **) * cm->M));
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], (sizeof(float *) * (W+1)));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
  }
  ESL_ALLOC(bestr, (sizeof(int) * (W+1)));

  /*****************************************************************
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   *****************************************************************/ 
  for (v = cm->M-1; v >= 0; v--)
    {
      alpha[v][0][0] = IMPOSSIBLE;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->endsc[v];
	  /* treat EL as emitting only on self transition */
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    if ((sc = alpha[y+yoffset][0][0] + cm->tsc[v][yoffset]) > alpha[v][0][0]) 
	      alpha[v][0][0] = sc;
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  if (alpha[v][0][0] < IMPOSSIBLE) alpha[v][0][0] = IMPOSSIBLE;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  /* Impose the bands.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  for (v = 0; v < cm->M; v++)
    {
      if(cm->stid[v] == BEGL_S) jmax = W; 
      else jmax = 1;

      for (d = 0; d < dmin[v] && d <=W; d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
      
      for (d = dmax[v]+1; d <= W;      d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
    }

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma,  sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback,  sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offset j index in gamma* data structures */
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* determine max d we're allowing for this state v and this position j */
	  curr_dmax = (gamma_j < dmax[v]) ? gamma_j : dmax[v];
	  if(curr_dmax > W) curr_dmax = W;
	  curr_dmin = (dmin[v] > 0) ? dmin[v] : 1;
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[v][yoffset]) > alpha[v][jp][d]) 
		      alpha[v][jp][d] = sc;
		  if (alpha[v][jp][d] < IMPROBABLE) alpha[v][jp][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  /* k is the length of the right fragment */
		  /* Careful, make sure k is consistent with bands in state w and state y. */
		  if(dmin[y] > (d-dmax[w])) k = dmin[y];
		  else k = d-dmax[w];
		  if(k < 0) k = 0;

		  if(dmax[y] < (d-dmin[w])) kmax = dmax[y];
		  else kmax = d-dmin[w];
		  
		  for (; k <= kmax; k++)
		    { 
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      if ((sc = alpha[w][jp][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			alpha[v][cur][d] = sc;
		    }
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      /* determine max d we're allowing for the root state and this position j */
      curr_dmax = (gamma_j < dmax[0]) ? gamma_j : dmax[0];
      if(curr_dmax > W) curr_dmax = W;
      curr_dmin = (dmin[0] > 0) ? dmin[0] : 1;

      for (d = curr_dmin; d <= curr_dmax; d++)
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
	      alpha[0][cur][d] = sc;
	  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	}
      
      /* EPN 11.09.05 
       * The following loop that deals with local begins was modified
       * to enforce bands on all states y that are possible internal entry
       * points. Old code block went from [0..d] in the d dimension
       * for state y.
       * ref: ~nawrocki/notebook/5_1109_inf_local_banded_spd/00LOG
       */

      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  d = (dmin[y] > dmin[0]) ? dmin[y]:dmin[0];
	  for (; (d <= dmax[y] && d <= gamma_j) && d <= W; d++)
	    {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	      if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	      if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	    }
	}
      }
      
      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  for (d = dmin[0]; d <= curr_dmax; d++) 
	    {
	      i = j-d+1;
	      gamma_i = j-d+1-i0+1;
	      sc = gamma[gamma_i-1] + alpha[0][cur][d];
	      if (sc > gamma[gamma_j])
		{
		  gamma[gamma_j]  = sc;
		  gback[gamma_j]  = i;
		  savesc[gamma_j] = alpha[0][cur][d]; 
		  saver[gamma_j]  = bestr[d];
		}
	    }
	}
      else
	{
	  /* Resolving overlaps greedily (RSEARCH style),  
	   * At least one hit is sent back for each j here.
	   * However, some hits can already be removed for the greedy overlap
	   * resolution algorithm.  Specifically, at the given j, any hit with a
	   * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
	  /* First, report hit with d of 1 if > cutoff */
	  if (alpha[0][cur][dmin[v]] >= cutoff) 
	    if(results != NULL) 
	      report_hit (j, j, bestr[1], alpha[0][cur][dmin[v]], results);
	  bestd = 1;
	  if (alpha[0][cur][dmin[v]] > best_score) 
	    best_score = alpha[0][cur][dmin[v]];
	  
	  /* Now, if current score is greater than maximum seen previous, report
	     it if >= cutoff and set new max */
	  for (d = dmin[v]; d <= curr_dmax; d++) 
	    {
	      if (alpha[0][cur][d] > best_score) 
		best_score = alpha[0][cur][d];
	      if (alpha[0][cur][d] > alpha[0][cur][bestd]) 
		{
		  if (alpha[0][cur][d] >= cutoff)
		    if(results != NULL) 
		      report_hit (j-d+1, j, bestr[d], alpha[0][cur][d], results);
		  bestd = d;
		}
	    }
	}
    } /* end loop over end positions j */
      
  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	continue;
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
    {
      j     = j0;
      while (j >= i0) 
	{
	  gamma_j = j-i0+1;
	  if (gback[gamma_j] == -1) /* no hit */
	    j--; 
	  else                /* a hit, a palpable hit */
	    {
	      if(savesc[gamma_j] > best_score) 
		best_score = savesc[gamma_j];
	      if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
		report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	      j = gback[gamma_j]-1;
	    }
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);
  
  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;
  /* printf("CYKBandedScan() return score: %10.4f\n", best_score); */
  return best_score;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

/******************************************************************
 * Function: InsideBandedScan()
 * EPN 05.09.06
 * 
 * Derived from CYKBandedScan() by SRE
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm.
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *           Log sums are performed (slowly) on floats by the LogSum2()
 *           function.
 *
 *           dmin,dmax set the bounds. Both are arrays, 0..v..cm->M. 
 *           The band for v is dmin[v]..dmax[v], inclusive; that is,
 *           dmin[v] is the minimum allowed d for state v (inclusive);
 *           dmax[v] is the maximum; that is, v can only account for
 *           subsequences of length >= dmin[v] and <= dmax[v].
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; 1..L
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results    - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float
InsideBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
		 float cutoff, search_results_t *results)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmax;               /* for B_st's, maximum k value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  int       curr_dmin;          /* temporary value for min d in for loops */
  int       curr_dmax;          /* temporary value for max d in for loops */

  /* Contract check */
  if(j0 < i0)
    cm_Fail("in InsideBandedScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)
    cm_Fail("ERROR in InsideBandedScan, dsq is NULL.");
  if(!(cm->flags & CMH_QDB))
    cm_Fail("in InsideBandedScan, QDBs invalid\n");

  /*printf("in InsideBandedScan i0: %d j0: %d\n", i0, j0);*/
  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*PrintDPCellsSaved(cm, dmin, dmax, W);*/

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   *****************************************************************/
  ESL_ALLOC(alpha, sizeof(float **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * (W+1));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], sizeof(float) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], sizeof(float) * (W+1));
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   *****************************************************************/ 
  for (v = cm->M-1; v >= 0; v--)
    {
      alpha[v][0][0] = IMPOSSIBLE;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->endsc[v];
	  /* treat EL as emitting only on self transition */
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    alpha[v][0][0] = LogSum2(alpha[v][0][0], (alpha[y+yoffset][0][0] 
						     + cm->tsc[v][yoffset]));
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  if (alpha[v][0][0] < IMPOSSIBLE) alpha[v][0][0] = IMPOSSIBLE;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  /* Impose the bands.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  for (v = 0; v < cm->M; v++)
    {
      if(cm->stid[v] == BEGL_S) jmax = W; 
      else jmax = 1;

      for (d = 0; d < dmin[v] && d <=W; d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
      
      for (d = dmax[v]+1; d <= W;      d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
    }

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offeset j index in gamma* data structures */
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* determine max d we're allowing for this state v and this position j */
	  curr_dmax = (gamma_j < dmax[v]) ? gamma_j : dmax[v];
	  if(curr_dmax > W) curr_dmax = W;
	  curr_dmin = (dmin[v] > 0) ? dmin[v] : 1;

	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][jp][d] = LogSum2(alpha[v][jp][d], (alpha[y+yoffset][cur][d] 
							       + cm->tsc[v][yoffset]));
		  if (alpha[v][jp][d] < IMPROBABLE) alpha[v][jp][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][prv][d-2] 
								 + cm->tsc[v][yoffset]));
		  
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][cur][d-1] 
								 + cm->tsc[v][yoffset]));
		  
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][prv][d-1] 
								 + cm->tsc[v][yoffset]));
		  
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));

		  /*EPN : Make sure k is consistent with bands in state w and state y. */

		  /* k is the length of the right fragment */
		  if(dmin[y] > (d-dmax[w])) k = dmin[y];
		  else k = d-dmax[w];
		  if(k < 0) k = 0;

		  if(dmax[y] < (d-dmin[w])) kmax = dmax[y];
		  else kmax = d-dmin[w];
		  
		  for (; k <= kmax; k++)
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[w][jp][d-k] 
								   + alpha[y][cur][k]));
		    }
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      /* determine max d we're allowing for the root state and this position j */
      curr_dmax = (gamma_j < dmax[0]) ? gamma_j : dmax[0];
      if(curr_dmax > W) curr_dmax = W;
      curr_dmin = (dmin[0] > 0) ? dmin[0] : 1;

      for (d = curr_dmin; d <= curr_dmax; d++)
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    alpha[0][cur][d] = LogSum2(alpha[0][cur][d], (alpha[y+yoffset][cur][d] 
							 + cm->tsc[0][yoffset]));

	  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	}
      
      /* EPN 11.09.05 
       * The following loop that deals with local begins was modified
       * to enforce bands on all states y that are possible internal entry
       * points. Old code block went from [0..d] in the d dimension
       * for state y.
       * ref: ~nawrocki/notebook/5_1109_inf_local_banded_spd/00LOG
       */
      /* EPN 05.08.06 
       * For Inside version, left local the same as it was in CYKScan(), I think this 
       * is okay, we're saying that the best root of the alignment is the state that
       * gives the highest Inside score. 
       */

      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  d = (dmin[y] > dmin[0]) ? dmin[y]:dmin[0];
	  /*if (dmin[y] > dmin[0]) d = dmin[y];
	    else d = dmin[0];*/
	  for (; (d <= dmax[y] && d <= gamma_j) && d <= W; d++)
	    {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      /*alpha[0][cur][d] = LogSum2(alpha[0][cur][d], sc);*/ /* EPN, Tue Nov  6 05:39:01 2007 to make it match cm_fastsearch.c */
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	      if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	      if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	    }
	}
      }
      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = dmin[0]; d <= curr_dmax; d++) 
	{
	  i = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  sc = gamma[gamma_i-1] + alpha[0][cur][d];
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = alpha[0][cur][d]; 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	continue;
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }

  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  return best_score;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* Function: InsideBandedScan_jd()
 * Date:     EPN, 07.20.06 Benasque, Spain
 * based on: hbandcyk.c::CYKBandedScan_jd() which is based on
 *           CYKBandedScan() by SRE in bandcyk.c
 *
 * Purpose:  Scan a (sub)sequence for matches to a covariance model, using HMM
 *           derived bands on the j and d dimensions. Intended for use on 
 *           subsequences with endpoints i and j, where i and j were determined
 *           using a HMM scan of the full sequence. This function then refines
 *           the positions of i and j, as well as deriving a CYK score that is
 *           more informative than an HMM based score. 
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Stores (i,j) segments that that have high inside scores: the sum
 *           over all possible alignments with start point [i,j] and end 
 *           points [i,j]. Derived from scancyk.c.
 *           Log sums are performed (slowly) on floats by the LogSum2()
 *           function.
 *
 *           jmin, jmax set the state specific bounds on the j dimension. 0..v..cm->M-1.
 *           hdmin, hdmax set the state and j specific bounds on the d dimension, indexed
 *           [0..v..cm-M-1][0..(jmax[v]-jmin[v]+1)].
 *           
 *           The j band for v is jmin[v]..jmax[v], inclusive; that is,
 *           jmin[v] is the minimum allowed j for state v (inclusive);
 *           jmax[v] is the maximum; 
 * 
 *           The d bands are v and j specific (in contrast to the a priori d bands
 *           which are only v specific), the d band for v and j is 
 *           hdmin[v][j-jmin[v]]..hdmax[v][j-jmin[v]] inclusive;
 *           hdmin[v][j-jmin[v]] is the minimum allowed d for state v and end point j
 *           hdmax[v][j-jmin[v]] is the maximum; 
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; i0..j0
 *           jmin      - minimum bound on j for state v; 0..M
 *           jmax      - maximum bound on j for state v; 0..M
 *           hdmin     - minimum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           hdmax     - maximum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - max d: max size of a hit
 *           ret_nhits - RETURN: number of hits (found only within current function call from i0 to j0)
 *           ret_hitr  - RETURN: start states of hits, 0..nhits-1
 *           ret_hiti  - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj  - RETURN: end positions of hits, 0..nhits-1
 *           ret_hitsc - RETURN: scores of hits, 0..nhits-1            
 *           min_thresh- minimum score to report (EPN via Alex Coventry 03.11.06)
 *
 * Returns:  
 *           hiti, hitj, hitsc are allocated here if nec; caller free's w/ free().
 */
void
InsideBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
		 int j0, int W, int *ret_nhits, int **ret_hitr, int **ret_hiti, 
		 int **ret_hitj, float **ret_hitsc, float min_thresh)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int      gamma_j;             /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int      gamma_i;             /* i index in the gamma matrix */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp_roll;   	        /* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       tmp_dmin, tmp_dmax; /* temp variables for ensuring we stay within d bands within loops */
  int       tmp_kmin, tmp_kmax; /* temp vars for B_st's, min/max k values consistent with bands*/

  int      jp_v, jp_y, jp_w;    /* mem eff banded j index in states v, y, and z 
				 * jp_x = j-jmin[x] */
  int      L;                   /* length of subsequence (j0-i0+1) */
  int      x;
  int      tmp_y;
  float    tmp_alpha_w, tmp_alpha_y;
  int       nhits;		/* # of hits in optimal parse */
  int      *hitr;		/* initial state indices of hits in optimal parse */
  int      *hiti;               /* start positions of hits in optimal parse */
  int      *hitj;               /* end positions of hits in optimal parse */
  float    *hitsc;              /* scores of hits in optimal parse */
  int       alloc_nhits;	/* used to grow the hit arrays */
  void     *tmp;                /* for ESL_RALLOC() */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in InsideBandedScan_jd, dsq is NULL.");

  /* EPN 08.11.05 Next line prevents wasteful computations when imposing
   * bands before the main recursion.  There is no need to worry about
   * alpha cells corresponding to subsequence distances within the windowlen
   * (W) but LONGER than the full sequence (L).  Saves a significant amount 
   * of time only if W is much larger than necessary, and the search sequences 
   * are short.
   */
  L = j0-i0+1;
  if (W > L) W = L; 

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));*/

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that unlike the other CYK scan functions, E memory is not shared: 
   * this is because the E deck will be different for different j values
   * due to the j bands. 
   * 
   *****************************************************************/
  ESL_ALLOC(alpha, sizeof(float **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * (W+1));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], sizeof(float) * (W+1));
      }
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], sizeof(float) * (W+1));
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {

      gamma_j = j-i0+1;
      cur = j%2;
      prv = (j-1)%2;

      /*****************************************************************
       * alpha initializations.
       * For the jd (HMM) banded strategy, we initialize inside the j loop,
       * because no cells are j-independent: for j's outside
       * the bands for a state v, should have ALL cells = IMPOSSIBLE.
       *****************************************************************/ 
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* Check to see if we're outside the bounds on j */
	  if(j < jmin[v] || j > jmax[v])
	    {
	      if (cm->stid[v] == BEGL_S) jp_roll = j % (W+1); else jp_roll = cur;
	      for (d = 0; d <= W; d++) 
		alpha[v][jp_roll][d] = IMPOSSIBLE;
	      /* Special boundary case: have to initialize alpha[v][prv] also */
	      if (j == i0)
		for (d = 0; d <= W; d++) 
		  alpha[v][prv][d] = IMPOSSIBLE;
	      continue;
	    }
	  
	  /* else we initialize on d = 0 */
	  alpha[v][cur][0] = IMPOSSIBLE;

	  if      (cm->sttype[v] == E_st)  alpha[v][cur][0] = 0;
	  else if (cm->sttype[v] == MP_st) alpha[v][cur][1] = alpha[v][prv][1] = IMPOSSIBLE;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[v][cur][0] = cm->endsc[v];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[v][cur][0] = LogSum2(alpha[v][cur][0], (alpha[y+yoffset][cur][0] 
							     + cm->tsc[v][yoffset]));
	      /* ...we don't bother to look at local alignment starts here... */
	      bestr[cur] = -1;
	      if (alpha[v][cur][0] < IMPOSSIBLE) alpha[v][cur][0] = IMPOSSIBLE;	
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* w is BEGL_S */
	      y = cm->cnum[v];   /* y is BEGR_S */
	      /* original line: 
	       * alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	       * we can't use that because alpha[w][0][0] and alpha[y][0][0] 
	       * may have been reset to IMPOSSIBLE, so we recalculate what they
	       * should be (this is wasteful):
	       */
	      tmp_y = cm->cfirst[w];
	      tmp_alpha_w = cm->endsc[w];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[w]; yoffset++)
		tmp_alpha_w = LogSum2(tmp_alpha_w, (alpha[tmp_y+yoffset][cur][0] 
						    + cm->tsc[w][yoffset]));
	      tmp_y = cm->cfirst[y];
	      tmp_alpha_y = cm->endsc[y];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[y]; yoffset++)
		tmp_alpha_y = LogSum2(tmp_alpha_y, (alpha[tmp_y+yoffset][cur][0] 
						    + cm->tsc[y][yoffset]));
	      alpha[v][cur][0] = tmp_alpha_w + tmp_alpha_y;
	      if (alpha[v][cur][0] < IMPOSSIBLE) alpha[v][cur][0] = IMPOSSIBLE;	
	    }
	  
	  /* Special boundary case: have to initialize alpha[v][prv] also */
	  if(j == i0) alpha[v][prv][0] = alpha[v][cur][0];

	  if (cm->stid[v] == BEGL_S) 
	    {
	      alpha[v][prv][0] = alpha[v][cur][0];
	      for (x = 2; x <= W; x++) 
		alpha[v][x][0] = alpha[v][cur][0];
	    }
	  /* done initialization */
	  
	  jp_v = j - jmin[v];
	  /* Impose the bands.
	   *   We have to do this inside the main loop because d bands are
	   *   dependent on v AND j. 
	   */
	  if (cm->stid[v] == BEGL_S) jp_roll = j % (W+1); else jp_roll = cur;
	  for (d =0; d < hdmin[v][jp_v] && d <=W; d++) 
	    alpha[v][jp_roll][d] = IMPOSSIBLE;
	  for (d = hdmax[v][jp_v]+1; d <= W;      d++) 
	    alpha[v][jp_roll][d] = IMPOSSIBLE;
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp_roll][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][jp_roll][d] = LogSum2(alpha[v][jp_roll][d], (alpha[y+yoffset][cur][d] 
									  + cm->tsc[v][yoffset]));
		  if (alpha[v][jp_roll][d] < IMPROBABLE) alpha[v][jp_roll][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][prv][d-2] 
								  + cm->tsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][cur][d-1] 
								  + cm->tsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[y+yoffset][prv][d-1] 
								  + cm->tsc[v][yoffset]));
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      /* Five inequalities must be satisfied to ensure that j and k 
	       * and k combinations correspond with alpha cells within the bands 
	       * on states y and w. 
	       * Below: jp_y = j - jmin[y] & jp_w = j - jmin[w]
	       *
	       * (1) j   >= jmin[y]          && j   <= jmax[y]
	       * (2) j-k >= jmin[w]          && j-k <= jmax[w]
	       * (3) k   >= hdmin[y][jp_y]   && k   <= hdmax[y][jp_y]
	       * (4) d-k >= hdmin[w][jp_w-k] && d-k <= hdmax[w][jp_w-k]
	       * (5) d   >= hdmin[v][jp_v]   && d   <= hdmax[v][jp_v]
	       */

	      /* initialize with endsc for all valid d for B state v */
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		}
	      /* Following code is careful, and not 'efficient' */
	      if(j >= jmin[y] && j <= jmax[y]) /* ensures (1): that j is valid for state y */
		{
		  jp_y = j - jmin[y];
		  jp_w = j - jmin[w]; 
		  i = j-d+1;
		  /*
		    printf("valid j: %d | jp_y: %d | jp_w: %d\n", j, jp_y, jp_w);
		    printf("hdmin[v][jp_v]: %d | hdmin[y][jp_y]: %d\n", hdmin[v][jp_v], hdmin[y][jp_y]);
		    printf("hdmax[v][jp_v]: %d | hdmax[y][jp_y]: %d\n", hdmax[v][jp_v], hdmax[y][jp_y]);
		  */
		  for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		    {
		      /* k is the length of the right fragment */
		      tmp_kmin = ((j-jmax[w]) > hdmin[y][jp_y]) ? (j-jmax[w]) : hdmin[y][jp_y];
		      if(tmp_kmin < 0) tmp_kmin = 0;
		      if(tmp_kmin < d-hdmax[w][jp_w-tmp_kmin]) tmp_kmin = d-hdmax[w][jp_w-tmp_kmin];
		      /* tmp_kmin is now smallest k that satisfies (2), (3), and (4) */

		      tmp_kmax = ((j-jmin[w]) < hdmax[y][jp_y]) ? (j-jmin[w]) : hdmax[y][jp_y];
		      if(tmp_kmax > d-hdmin[w][jp_w-tmp_kmax]) tmp_kmax = d-hdmin[w][jp_w-tmp_kmax];
		      /* tmp_kmax is now largest k that satisfies (2), (3), and (4) */
		      /*printf("tmp_kmin: %d | tmp_kmax: %d\n", tmp_kmin, tmp_kmax);*/
		      for (k = tmp_kmin; k <= tmp_kmax; k++)
			{
			  jp_roll = (j-k)%(W+1); /* jp_roll is rolling index into BEGL_S (state w) 
						  * deck j dimension */
			  alpha[v][cur][d] = LogSum2(alpha[v][cur][d], (alpha[w][jp_roll][d-k] 
									+ alpha[y][cur][k]));
			}
		      if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		      /*printf("B alpha[%d][%d][%d]: %f\n", v, cur, d, alpha[v][cur][d]);*/
		    }
		}
	    }
	} /* end loop over decks v>0 */
	  
      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      /* Check to see if we're within bounds on j */
      if(j < jmin[0] || j > jmax[0])
	{
	  for (d = 0; d <= W; d++) 
	    alpha[0][jp_roll][d] = IMPOSSIBLE;
	  /* Inform the little semi-Markov model that deals with multihit parsing
	   * that a hit is impossible, j is outside root band on j:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  continue;
	}
      /* if we get here, j is within ROOT_S state 0's band */

      /* first initialize on d = 0 */
      alpha[0][0][0] = IMPOSSIBLE;
      y = cm->cfirst[v];
      alpha[0][0][0] = cm->endsc[v];
      /* treat EL as emitting only on self transition */
      for (yoffset = 0; yoffset < cm->cnum[0]; yoffset++)
	alpha[0][0][0] = LogSum2(alpha[0][0][0], (alpha[y+yoffset][0][0] 
						  + cm->tsc[0][yoffset]));
      /* ...we don't bother to look at local alignment starts here... */
      bestr[0] = -1;
      if (alpha[0][0][0] < IMPOSSIBLE) alpha[0][0][0] = IMPOSSIBLE;	
      alpha[0][1][0] = alpha[0][0][0];
      /* done initialization on d = 0 */

      jp_v = j - jmin[0];
      /* Impose the bands.
       *   We have to do this here because d bands are
       *   dependent on v AND j. 
       */
      jp_roll = cur;
      for (d =0; d < hdmin[0][jp_v] && d <=W; d++) 
	alpha[0][jp_roll][d] = IMPOSSIBLE;
      for (d = hdmax[0][jp_v]+1; d <= W;      d++) 
	alpha[0][jp_roll][d] = IMPOSSIBLE;
      
      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    alpha[0][cur][d] = LogSum2(alpha[0][cur][d], (alpha[y+yoffset][cur][d] 
							  + cm->tsc[0][yoffset]));

	  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	}

      if (cm->flags & CMH_LOCAL_BEGIN) {
	/* (comment from scaninside.c::CYKScan()) 
	 * EPN 05.08.06 Left local the same as it was in CYKScan(), I think this 
	 * is okay, we're saying that the best root of the alignment is the state that
	 * gives the highest Inside score. 
	 */	
	for (y = 1; y < cm->M; y++) {
	  if(j >= jmin[y] && j <= jmax[y]) 
	    {
	      jp_y = j - jmin[y];
	      tmp_dmin = (hdmin[y][jp_y] > hdmin[0][jp_v]) ? hdmin[y][jp_y] : hdmin[0][jp_v];
	      tmp_dmax = (hdmax[y][jp_y] < hdmax[0][jp_v]) ? hdmax[y][jp_y] : hdmax[0][jp_v];
	      if(tmp_dmax > j) tmp_dmax = j;
	      for (d = tmp_dmin; d <= tmp_dmax; d++)
		{
		  if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
		  else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
		  if (sc > alpha[0][cur][d]) {
		    alpha[0][cur][d] = sc;
		    bestr[d]         = y;
		  }
		  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
		}
	    }
	}
      }
      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = hdmin[0][jp_v]; ((d <= hdmax[0][jp_v] && d <= gamma_j) && d <= W); d++) 
	{
	  i       = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  sc = gamma[gamma_i-1] + alpha[0][cur][d]  - min_thresh; 
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = alpha[0][cur][d]; 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  alloc_nhits = 10;
  ESL_ALLOC(hitr, sizeof(int)   * alloc_nhits);
  ESL_ALLOC(hitj, sizeof(int)   * alloc_nhits);
  ESL_ALLOC(hiti, sizeof(int)   * alloc_nhits);
  ESL_ALLOC(hitsc,sizeof(float) * alloc_nhits);

  j     = j0;
  nhits = 0;
  while (j >= i0) {
    gamma_j = j-i0+1;
    if (gback[gamma_j] == -1) /* no hit */
      j--; 
    else                /* a hit, a palpable hit */
      {
	hitr[nhits]   = saver[gamma_j];
	hitj[nhits]   = j;
	hiti[nhits]   = gback[gamma_j];
	hitsc[nhits]  = savesc[gamma_j];
	//printf("new hit num %d | i: %d | j: %d | sc: %f\n", nhits, hiti[nhits], hitj[nhits], hitsc[nhits]); 
	nhits++;
	j = gback[gamma_j]-1;
	
	if (nhits == alloc_nhits) {
	  ESL_RALLOC(hitr, tmp, sizeof(int)   * (alloc_nhits + 10));
	  ESL_RALLOC(hitj, tmp, sizeof(int)   * (alloc_nhits + 10));
	  ESL_RALLOC(hiti, tmp, sizeof(int)   * (alloc_nhits + 10));
	  ESL_RALLOC(hitsc,tmp, sizeof(float) * (alloc_nhits + 10));
	  alloc_nhits += 10;
	}
      }
  }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  //printf("end nhis: %d\n", nhits);
  *ret_nhits = nhits;
  *ret_hitr  = hitr;
  *ret_hiti  = hiti;
  *ret_hitj  = hitj;
  *ret_hitsc = hitsc;

  /*for (i = 0; i < nhits; i++)
    {
      printf("j0: %d end hit %-4d: %6d %6d %8.2f bits\n", j0, i, 
	     hiti[i], 
	     hitj[i],
	     hitsc[i]);
      
    }*/
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
}

/**************************************************************
 * Function: iInsideScan()
 * EPN, Tue Dec 19 11:54:43 2006
 * 
 * Based on CYKScan() by SRE
 *
 * Purpose:  Scan a (sub)sequence for matches to a covariance model.
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Log sums are performed using scaled ints with ILogsum().
 *           -INFTY is -987654321, but treated as valid score! This
 *           seems okay.
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; 1..L
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results    - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score (float - not scaled int) of best overall hit
 */
float 
iInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, 
	    float cutoff, search_results_t *results)

{
  int       status;
  int    ***alpha;              /* inside DP score matrix (scaled ints), [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in iInsideScan, dsq is NULL.");

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   *****************************************************************/
  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*printf("in iInsideScan i0: %d j0: %d\n", i0, j0);*/
  ESL_ALLOC(alpha, sizeof(int **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], sizeof(int *) * (W+1));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], sizeof(int) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(int *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], sizeof(int) * (W+1));
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *****************************************************************/ 
  for (v = cm->M-1; v >= 0; v--)
    {
      alpha[v][0][0] = -INFTY;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = -INFTY;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->iendsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    alpha[v][0][0] = ILogsum(alpha[v][0][0], (alpha[y+yoffset][0][0] 
						      + cm->itsc[v][yoffset]));
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  /* ! */
	  if (alpha[v][0][0] < -INFTY) alpha[v][0][0] = -INFTY;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  for (d = 1; d <= W; d++)
    alpha[cm->M-1][0][d] = alpha[cm->M-1][1][d] = -INFTY;

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0; 
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc,sizeof(float) * (L+1));
  ESL_ALLOC(saver, sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offeset j index in gamma* data structures */
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][jp][d] = ILogsum(alpha[v][jp][d], (alpha[y+yoffset][cur][d] 
								+ cm->itsc[v][yoffset]));
		  if (alpha[v][jp][d] < -INFTY) alpha[v][jp][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = 2; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-2] 
								  + cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v]))); 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][cur][d-1] 
								  + cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);

		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-1] 
								 + cm->itsc[v][yoffset]));
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);

		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));

		  for (k = 0; k <= d; k++) /* k is length of right fragment */
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[w][jp][d-k] 
								    + alpha[y][cur][k]));
		    }
		  /* ! */
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      for (d = 1; d <= W && d <= gamma_j; d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->itsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    alpha[0][cur][d] = ILogsum(alpha[0][cur][d], (alpha[y+yoffset][cur][d] 
							  + cm->itsc[0][yoffset]));
	  if (cm->flags & CMH_LOCAL_BEGIN) {
	    /* EPN 05.08.06 Left local the same as it was in CYKScan(), I think this 
	     * is okay, we're saying that the best root of the alignment is the state that
	     * gives the highest Inside score. 
	     */
	    for (y = 1; y < cm->M; y++) {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->ibeginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->ibeginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	    }
	  }
	  if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
	  if (Scorify(alpha[0][cur][d]) > best_neg_score) best_neg_score = Scorify(alpha[0][cur][d]);
	}

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; 
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = 1; d <= W && d <= gamma_j; d++) 
	{
	  i = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  sc = gamma[gamma_i-1] + Scorify(alpha[0][cur][d]);
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = Scorify(alpha[0][cur][d]); 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	continue;
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  /*printf("iInsideScan() returning best_score: %f best_neg_score: %f\n", best_score, best_neg_score);*/
  return best_score;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/******************************************************************
 * Function: iInsideBandedScan()
 * EPN, Tue Dec 19 09:43:14 2006
 * 
 * Derived from CYKBandedScan() by SRE
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm, and using scaled int log odd scores
 *           (normal InsideBandedScan() uses float log odd scores, which
 *            makes it much slower).
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *           Log sums are performed using scaled ints with ILogsum().
 *
 *           dmin,dmax set the bounds. Both are arrays, 0..v..cm->M. 
 *           The band for v is dmin[v]..dmax[v], inclusive; that is,
 *           dmin[v] is the minimum allowed d for state v (inclusive);
 *           dmax[v] is the maximum; that is, v can only account for
 *           subsequences of length >= dmin[v] and <= dmax[v].
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; 1..L
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results    - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score (float - not scaled int) of best overall hit
 */
float
iInsideBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
		  float cutoff, search_results_t *results)
{
  int       status;
  int    ***alpha;              /* inside DP score matrix (scaled ints), [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmax;               /* for B_st's, maximum k value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  int       curr_dmin;          /* temporary value for min d in for loops */
  int       curr_dmax;          /* temporary value for max d in for loops */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */

  /* Contract check */
  if(j0 < i0)
    cm_Fail("in iInsideBandedScan, i0: %d j0: %d\n", i0, j0);
  if(!(cm->flags & CMH_QDB))
    cm_Fail("in iInsideBandedScan, QDBs invalid\n");
  if(dsq == NULL)
    cm_Fail("ERROR in iInsideBandedScan, dsq is NULL.");

  /*printf("in iInsideBandedScan i0: %d j0: %d\n", i0, j0);*/

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*PrintDPCellsSaved(cm, dmin, dmax, W);*/

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   *****************************************************************/
  ESL_ALLOC(alpha, sizeof(int **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], sizeof(int *) * (W+1));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], sizeof(int) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(int *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], sizeof(int) * (W+1));
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));
  
  /*****************************************************************
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   *****************************************************************/ 
  for (v = cm->M-1; v >= 0; v--)
    {
      alpha[v][0][0] = -INFTY;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = -INFTY;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->iendsc[v];
	  /* treat EL as emitting only on self transition */
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    alpha[v][0][0] = ILogsum(alpha[v][0][0], (alpha[y+yoffset][0][0] 
						     + cm->itsc[v][yoffset]));
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  /* ! */
	  if (alpha[v][0][0] < -INFTY) alpha[v][0][0] = -INFTY;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  /* Impose the bands.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  for (v = 0; v < cm->M; v++)
    {
      if(cm->stid[v] == BEGL_S) jmax = W; 
      else jmax = 1;

      for (d = 0; d < dmin[v] && d <=W; d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = -INFTY;
      
      for (d = dmax[v]+1; d <= W;      d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = -INFTY;
    }

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   * gamma contains real (float) scores, not scaled ints.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc,sizeof(float) * (L+1));
  ESL_ALLOC(saver, sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offeset j index in gamma* data structures */
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* determine max d we're allowing for this state v and this position j */
	  curr_dmax = (gamma_j < dmax[v]) ? gamma_j : dmax[v];
	  if(curr_dmax > W) curr_dmax = W;
	  curr_dmin = (dmin[v] > 0) ? dmin[v] : 1;

	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		      alpha[v][jp][d] = ILogsum(alpha[v][jp][d], (alpha[y+yoffset][cur][d] 
								  + cm->itsc[v][yoffset]));
		  if (alpha[v][jp][d] < -INFTY) alpha[v][jp][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));

		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		      alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-2] 
								    + cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);

		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][cur][d-1] 
								  + cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-1] 
								  + cm->itsc[v][yoffset]));
		  
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);

		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));

		  /*EPN : Make sure k is consistent with bands in state w and state y. */

		  /* k is the length of the right fragment */
		  if(dmin[y] > (d-dmax[w])) k = dmin[y];
		  else k = d-dmax[w];
		  if(k < 0) k = 0;

		  if(dmax[y] < (d-dmin[w])) kmax = dmax[y];
		  else kmax = d-dmin[w];
		  
		  for (; k <= kmax; k++)
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[w][jp][d-k] 
								    + alpha[y][cur][k]));
		    }
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      /* determine max d we're allowing for the root state and this position j */
      curr_dmax = (gamma_j < dmax[0]) ? gamma_j : dmax[0];
      if(curr_dmax > W) curr_dmax = W;
      curr_dmin = (dmin[0] > 0) ? dmin[0] : 1;

      for (d = curr_dmin; d <= curr_dmax; d++)
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->itsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    alpha[0][cur][d] = ILogsum(alpha[0][cur][d], (alpha[y+yoffset][cur][d] 
							  + cm->itsc[0][yoffset]));

	  if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
	  if (Scorify(alpha[0][cur][d]) > best_neg_score) best_neg_score = Scorify(alpha[0][cur][d]);
	}
      
      /* EPN 11.09.05 
       * The following loop that deals with local begins was modified
       * to enforce bands on all states y that are possible internal entry
       * points. Old code block went from [0..d] in the d dimension
       * for state y.
       * ref: ~nawrocki/notebook/5_1109_inf_local_banded_spd/00LOG
       */
      /* EPN 05.08.06 
       * For Inside version, left local the same as it was in CYKScan(), I think this 
       * is okay, we're saying that the best root of the alignment is the state that
       * gives the highest Inside score. 
       */


      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  d = (dmin[y] > dmin[0]) ? dmin[y]:dmin[0];
	  /*if (dmin[y] > dmin[0]) d = dmin[y];
	    else d = dmin[0];*/
	  for (; (d <= dmax[y] && d <= gamma_j) && d <= W; d++)
	    {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->ibeginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->ibeginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	      if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
	      if (Scorify(alpha[0][cur][d]) > best_neg_score) best_neg_score = Scorify(alpha[0][cur][d]);
	    }
	}
      }
      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = dmin[0]; d <= curr_dmax; d++)
	{
	  i = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  sc = gamma[gamma_i-1] + Scorify(alpha[0][cur][d]);
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = Scorify(alpha[0][cur][d]); 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	continue;
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  /*printf("iInsideBandedScan() returning best_score: %f best_neg_score: %f\n", best_score, best_neg_score);*/
  return best_score;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}


/* Function: CYKScan()
 * Date:     SRE, Thu May  2 11:56:11 2002 [AA 3050 SFO->STL]
 *
 * Purpose:  Scan a sequence for matches to a covariance model.
 *           Multiple nonoverlapping hits and local alignment.
 *
 * Args:     cm          - the covariance model
 *           dsq         - digitized sequence to search; i0..j0
 *           i0          - start of target subsequence (1 for full seq)
 *           j0          - end of target subsequence (L for full seq)
 *           W           - max d: max size of a hit
 *           cutoff      - minimum score to report 
 *           results     - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float 
CYKScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, 
	float cutoff, search_results_t *results)
{
  /* Contract check */
  if(j0 < i0)
    cm_Fail("ERROR in CYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)
    cm_Fail("ERROR in CYKScan, dsq is NULL\n");

  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* index into BEGL_S decks: jp=j%(W+1) */
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   *****************************************************************/

  ESL_ALLOC(alpha, (sizeof(float **) * cm->M));
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], (sizeof(float *) * (W+1)));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
  }
  ESL_ALLOC(bestr, (sizeof(int) * (W+1)));

  /*****************************************************************
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *****************************************************************/ 
  for (v = cm->M-1; v >= 0; v--)
    {
      alpha[v][0][0] = IMPOSSIBLE;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->endsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    if ((sc = alpha[y+yoffset][0][0] + cm->tsc[v][yoffset]) > alpha[v][0][0]) 
	      alpha[v][0][0] = sc;
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  if (alpha[v][0][0] < IMPOSSIBLE) alpha[v][0][0] = IMPOSSIBLE;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  for (d = 1; d <= W; d++)
    alpha[cm->M-1][0][d] = alpha[cm->M-1][1][d] = IMPOSSIBLE;

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma,  sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback,  sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offset j index in gamma* data structures */
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[v][yoffset]) > alpha[v][jp][d]) 
		      alpha[v][jp][d] = sc;
		  if (alpha[v][jp][d] < IMPOSSIBLE) alpha[v][jp][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = 2; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v]))); 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  
		  for (k = 0; k <= d; k++) /* k is length of right fragment */
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      if ((sc = alpha[w][jp][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			alpha[v][cur][d] = sc; 
		    }
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      for (d = 1; d <= W && d <= gamma_j; d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
	      alpha[0][cur][d] = sc;
	  if (cm->flags & CMH_LOCAL_BEGIN) {
	    for (y = 1; y < cm->M; y++) {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	    }
	  }
	  if (alpha[0][cur][d] < IMPOSSIBLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	}

      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; 
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  for (d = 1; d <= W && d <= gamma_j; d++) 
	    {
	      i = j-d+1;
	      gamma_i = j-d+1-i0+1;
	      sc = gamma[gamma_i-1] + alpha[0][cur][d];
	      if (sc > gamma[gamma_j])
		{
		  gamma[gamma_j]  = sc;
		  gback[gamma_j]  = i;
		  savesc[gamma_j] = alpha[0][cur][d]; 
		  saver[gamma_j]  = bestr[d];
		}
	    }
	}
      else
	{
	  /* Resolving overlaps greedily (RSEARCH style),  
	   * At least one hit is sent back for each j here.
	   * However, some hits can already be removed for the greedy overlap
	   * resolution algorithm.  Specifically, at the given j, any hit with a
	   * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
	  /* First, report hit with d of 1 if > cutoff */
	  if (alpha[0][cur][1] >= cutoff) 
	    if(results != NULL) 
	      {
		report_hit (j, j, bestr[1], alpha[0][cur][1], results);
	      }
	  bestd = 1;
	  if (alpha[0][cur][1] > best_score) 
	    best_score = alpha[0][cur][1];
	  
	  /* Now, if current score is greater than maximum seen previous, report
	     it if >= cutoff and set new max */
	  for (d = 2; d <= W && d <= gamma_j; d++) 
	    {
	      if (alpha[0][cur][d] > best_score) 
		best_score = alpha[0][cur][d];
	      if (alpha[0][cur][d] > alpha[0][cur][bestd]) 
		{
		  if (alpha[0][cur][d] >= cutoff)
		    if(results != NULL) 
		      {
			report_hit (j-d+1, j, bestr[d], alpha[0][cur][d], results);
		      }
		  bestd = d;
		}
	    }
	}
    } /* end loop over end positions j */
  
  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	continue;
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
    {
      j     = j0;
      while (j >= i0) 
	{
	  gamma_j = j-i0+1;
	  if (gback[gamma_j] == -1) /* no hit */
	    j--; 
	  else                /* a hit, a palpable hit */
	    {
	      if(savesc[gamma_j] > best_score) 
		best_score = savesc[gamma_j];
	      if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
		report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	      j = gback[gamma_j]-1;
	    }
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);
  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;
  /* printf("CYKScan() return score: %10.4f\n", best_score); */
  return best_score;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

/* Function: CYKScanRequires()
 * Date:     SRE, Mon May  6 18:48:18 2002 [St. Louis]
 *
 * Purpose:  Return the memory required by CYKScan(), in megabytes.
 */
float
CYKScanRequires(CM_t *cm, int L, int W)
{
  float ram;
  int   v,j;

				/* alpha allocations */
  ram = (float) (sizeof(float **) * cm->M);
  for (v = cm->M-1; v >= 0; v--)
    {
      if (cm->stid[v] == BEGL_S)
	{
	  ram += (float) (sizeof(float *) * (W+1));
	  for (j = 0; j <= W; j++)
	    ram += (float) (sizeof(float) * (W+1));
	}
      else if (cm->sttype[v] == E_st && v < cm->M-1) 
	continue;
      else 
	{
	  ram += (float) (sizeof(float *) * 2);
	  for (j = 0; j < 2; j++) 
	    ram += (float) (sizeof(float) * (W+1));
	}
    }
  ram += (float) (sizeof(float) * (L+1)); /* gamma allocation */
  ram += (float) (sizeof(int)   * (L+1)); /* gback allocation */
  ram += (float) (sizeof(float) * (L+1)); /* savesc allocation */
  return (ram / 1000000.);
}


/* Function: CYKBandedScan_jd() 
 * Date:     EPN, 06.22.06
 * based on: CYKBandedScan() by SRE in bandcyk.c
 *
 * Purpose:  Scan a (sub)sequence for matches to a covariance model, using HMM
 *           derived bands on the j and d dimensions. Intended for use on 
 *           subsequences with endpoints i and j, where i and j were determined
 *           using a HMM scan of the full sequence. This function then refines
 *           the positions of i and j, as well as deriving a CYK score that is
 *           more informative than an HMM based score. 
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *
 *           jmin, jmax set the state specific bounds on the j dimension. 0..v..cm->M-1.
 *           hdmin, hdmax set the state and j specific bounds on the d dimension, indexed
 *           [0..v..cm-M-1][0..(jmax[v]-jmin[v]+1)].
 *           
 *           The j band for v is jmin[v]..jmax[v], inclusive; that is,
 *           jmin[v] is the minimum allowed j for state v;
 *           jmax[v] is the maximum; 
 * 
 *           The d bands are v and j specific (in contrast to the a priori d bands
 *           which are only v specific), the d band for v and j is 
 *           hdmin[v][j-jmin[v]]..hdmax[v][j-jmin[v]] inclusive;
 *           hdmin[v][j-jmin[v]] is the minimum allowed d for state v and end point j
 *           hdmax[v][j-jmin[v]] is the maximum; 
 *
 * Args:     cm        - the covariance model
 *           sq        - digitized sequence to search; i0..j0
 *           jmin      - minimum bound on j for state v; 0..M
 *           jmax      - maximum bound on j for state v; 0..M
 *           hdmin     - minimum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           hdmax     - maximum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results   - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float
CYKBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
		 int j0, int W, float cutoff, search_results_t *results)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  float  ***alpha_mem;          /* pointers to original alpha memory */
  float    *imp_row;           /* an IMPOSSIBLE deck (full of IMPOSSIBLE scores), 
				 * pointed to when j is outside j band for v */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int      gamma_j;             /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int      gamma_i;             /* i index in the gamma matrix */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp_roll;   	        /* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       tmp_dmin, tmp_dmax; /* temp variables for ensuring we stay within d bands within loops */
  int       tmp_kmin, tmp_kmax; /* temp vars for B_st's, min/max k values consistent with bands*/

  int      jp_v, jp_y, jp_w;    /* mem eff banded j index in states v, y, and z 
				 * jp_x = j-jmin[x] */
  int      L;                   /* length of subsequence (j0-i0+1) */
  int      x;
  int      tmp_y;
  float    tmp_alpha_w, tmp_alpha_y;
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  
  /* Contract checks */
  if((!(cm->search_opts & CM_SEARCH_NOQDB)) && (cm->dmin == NULL || cm->dmax == NULL))
    cm_Fail("ERROR in CYKBandedScan_jd(), trying to use QDB, but cm->dmin, cm->dmax are NULL.\n");
  if(dsq == NULL)
    cm_Fail("ERROR in CYKBandedScan_jd(), dsq is NULL.");

  if(!(cm->search_opts & CM_SEARCH_NOQDB)) /* we're doing qdb */
    combine_qdb_hmm_d_bands(cm, jmin, jmax, hdmin, hdmax);

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  /*printf("in CYKBandedScan_jd i0: %5d | j0: %5d | L: %5d | W: %5d\n", i0, j0, L, W);*/
  if (W > L) W = L; /* shouldn't look longer than seq length L */

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));*/
  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that unlike the other CYK scan functions, E memory is not shared: 
   * this is because the E deck will be different for different j values
   * due to the j bands. 
   * 
   *****************************************************************/
  ESL_ALLOC(alpha, sizeof(float **) * cm->M);
  ESL_ALLOC(alpha_mem, sizeof(float **) * cm->M);
  /* we use alpha_mem to remember where each alpha row (alpha[v][j]) is
   * in case we've set alpha[v][cur] to imp_row (the precalc'ed IMPOSSIBLE row)
   * in a prior iteration, and we are about to overwrite it, and we don't
   * want to overwrite our special IMPOSSIBLE row.
   */
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], sizeof(float *)  * (W+1));
	ESL_ALLOC(alpha_mem[v], sizeof(float *)  * (W+1));
	for (j = 0; j <= W; j++)
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(float) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	ESL_ALLOC(alpha_mem[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(float) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

  /* Initialize the impossible deck, which we'll point to for 
   * j positions that are outside of the j band on v */
  ESL_ALLOC(imp_row, sizeof(float) * (W+1));
  for (d = 0; d <= W; d++) imp_row[d] = IMPOSSIBLE;
    
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {

      gamma_j = j-i0+1;
      cur = (j-i0+1)%2; /* cur == 1 when j == i0 */
      prv = (j-i0)  %2; /* prv == 0 when j == i0 */

      /*****************************************************************
       * alpha initializations.
       * For the jd (HMM) banded strategy, we initialize inside the j loop,
       * because no cells are j-independent: for j's outside
       * the bands for a state v, should have ALL cells = IMPOSSIBLE.
       *****************************************************************/ 
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* Check to see if we're outside the bounds on j */
	  if(j < jmin[v] || j > jmax[v])
	    {
	      if (cm->stid[v] == BEGL_S) 
		{
		  jp_roll = j % (W+1); 
		  for (d = 0; d <= W; d++) 
		    alpha[v][jp_roll][d] = IMPOSSIBLE;
		}
	      else
		{
		  jp_roll = cur;
		  alpha[v][jp_roll] = imp_row;
		}
	      /* Special boundary case: have to initialize alpha[v][prv] also */
	      if (j == i0)
		alpha[v][prv] = imp_row;
	      /*for (d = 0; d <= W; d++) 
		alpha[v][prv][d] = IMPOSSIBLE;*/

	      continue;
	    }
	  if(alpha[v][cur] == imp_row)
	    {
	      /* we don't want to overwrite imp_row */
	      alpha[v][cur] = alpha_mem[v][cur];
	      for (d = 0; d <= W; d++) 
		alpha[v][cur][d] = IMPOSSIBLE;
	    }
	  
	  /* else we initialize on d = 0 */
	  alpha[v][cur][0] = IMPOSSIBLE;
	  
	  if      (cm->sttype[v] == E_st)  alpha[v][cur][0] = 0;
	  else if (cm->sttype[v] == MP_st) alpha[v][cur][1] = alpha[v][prv][1] = IMPOSSIBLE;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[v][cur][0] = cm->endsc[v];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		if ((sc = alpha[y+yoffset][cur][0] + cm->tsc[v][yoffset]) > alpha[v][cur][0]) 
		  alpha[v][cur][0] = sc;
	      /* ...we don't bother to look at local alignment starts here... */
	      bestr[cur] = -1;
	      if (alpha[v][cur][0] < IMPOSSIBLE) alpha[v][cur][0] = IMPOSSIBLE;	
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* w is BEGL_S */
	      y = cm->cnum[v];   /* y is BEGR_S */
	      /* original line: 
	       * alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	       * we can't use that because alpha[w][0][0] and alpha[y][0][0] 
	       * may have been reset to IMPOSSIBLE, so we recalculate what they
	       * should be (this is wasteful):
	       */
	      tmp_y = cm->cfirst[w];
	      tmp_alpha_w = cm->endsc[w];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[w]; yoffset++)
		{
		  if ((sc = alpha[tmp_y+yoffset][cur][0] + cm->tsc[w][yoffset]) > tmp_alpha_w)
		    tmp_alpha_w = sc;
		}
	      tmp_y = cm->cfirst[y];
	      tmp_alpha_y = cm->endsc[y];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[y]; yoffset++)
		if ((sc = alpha[tmp_y+yoffset][cur][0] + cm->tsc[y][yoffset]) > tmp_alpha_y)
		  tmp_alpha_y = sc;
	      alpha[v][cur][0] = tmp_alpha_w + tmp_alpha_y;
	      if (alpha[v][cur][0] < IMPOSSIBLE) alpha[v][cur][0] = IMPOSSIBLE;	
	    }
	  
	  /* Special boundary case: have to initialize alpha[v][prv] also */
	  if(j == i0) alpha[v][prv][0] = alpha[v][cur][0];

	  if (cm->stid[v] == BEGL_S) 
	    {
	      alpha[v][prv][0] = alpha[v][cur][0];
	      for (x = 2; x <= W; x++) 
		alpha[v][x][0] = alpha[v][cur][0];
	    }
	  /* done initialization */
	  
	  jp_v = j - jmin[v];
	  /* Impose the bands.
	   *   We have to do this inside the main loop because d bands are
	   *   dependent on v AND j. 
	   */
	  if (cm->stid[v] == BEGL_S) jp_roll = j % (W+1); else jp_roll = cur;

	  for (d =0; d < hdmin[v][jp_v] && d <=W; d++) 
	    alpha[v][jp_roll][d] = IMPOSSIBLE;
	  for (d = hdmax[v][jp_v]+1; d <= W;      d++) 
	    alpha[v][jp_roll][d] = IMPOSSIBLE;
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp_roll][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[v][yoffset]) > alpha[v][jp_roll][d]) 
		      alpha[v][jp_roll][d] = sc;
		  if (alpha[v][jp_roll][d] < IMPROBABLE) alpha[v][jp_roll][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      /* Five inequalities must be satisfied to ensure that j and k 
	       * and k combinations correspond with alpha cells within the bands 
	       * on states y and w. 
	       * Below: jp_y = j - jmin[y] & jp_w = j - jmin[w]
	       *
	       * (1) j   >= jmin[y]          && j   <= jmax[y]
	       * (2) j-k >= jmin[w]          && j-k <= jmax[w]
	       * (3) k   >= hdmin[y][jp_y]   && k   <= hdmax[y][jp_y]
	       * (4) d-k >= hdmin[w][jp_w-k] && d-k <= hdmax[w][jp_w-k]
	       * (5) d   >= hdmin[v][jp_v]   && d   <= hdmax[v][jp_v]
	       */

	      /* initialize with endsc for all valid d for B state v */
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		}
	      /* Following code is careful, and not 'efficient' */
	      if(j >= jmin[y] && j <= jmax[y]) /* ensures (1): that j is valid for state y */
		{
		  jp_y = j - jmin[y];
		  jp_w = j - jmin[w]; 
		  i = j-d+1;
		  /*
		    printf("valid j: %d | jp_y: %d | jp_w: %d\n", j, jp_y, jp_w);
		    printf("hdmin[v][jp_v]: %d | hdmin[y][jp_y]: %d\n", hdmin[v][jp_v], hdmin[y][jp_y]);
		    printf("hdmax[v][jp_v]: %d | hdmax[y][jp_y]: %d\n", hdmax[v][jp_v], hdmax[y][jp_y]);
		  */
		  for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		    {
		      /* k is the length of the right fragment */
		      tmp_kmin = ((j-jmax[w]) > hdmin[y][jp_y]) ? (j-jmax[w]) : hdmin[y][jp_y];
		      if(tmp_kmin < 0) tmp_kmin = 0;
		      /* HEREHEREHEREHEREHEREHEREHERE, ensure that hdmax[w][jp_w-tmp_kmin] is valid
		       * before accessing it! */
		      if(tmp_kmin < d-hdmax[w][jp_w-tmp_kmin]) tmp_kmin = d-hdmax[w][jp_w-tmp_kmin];
		      /* tmp_kmin is now smallest k that satisfies (2), (3), and (4) */

		      tmp_kmax = ((j-jmin[w]) < hdmax[y][jp_y]) ? (j-jmin[w]) : hdmax[y][jp_y];
		      /* HEREHEREHEREHEREHEREHEREHERE, ensure that hdmin[w][jp_w-tmp_kmax] is valid
		       * before accessing it! */
		      if(tmp_kmax > d-hdmin[w][jp_w-tmp_kmax]) tmp_kmax = d-hdmin[w][jp_w-tmp_kmax];
		      /* tmp_kmax is now largest k that satisfies (2), (3), and (4) */
		      /*printf("tmp_kmin: %d | tmp_kmax: %d\n", tmp_kmin, tmp_kmax);*/
		      for (k = tmp_kmin; k <= tmp_kmax; k++)
			{
			  jp_roll = (j-k)%(W+1); /* jp_roll is rolling index into BEGL_S (state w) 
						  * deck j dimension */
			  if ((sc = alpha[w][jp_roll][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			    alpha[v][cur][d] = sc;
			}
		      if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		      /*printf("B alpha[%d][%d][%d]: %f\n", v, cur, d, alpha[v][cur][d]);*/
		    }
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      /* Check to see if we're within bounds on j */
      if(j < jmin[0] || j > jmax[0])
	{
	  /*printf("j: %d (gamma_j: %d) IMPOSSIBLE BABY! min: %d max: %d\n", j, gamma_j, jmin[0], jmax[0]);*/
	  for (d = 0; d <= W; d++) 
	    alpha[0][jp_roll][d] = IMPOSSIBLE; 
	  /* Inform the little semi-Markov model that deals with multihit parsing
	   * that a hit is impossible, j is outside root band on j:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  continue;
	}
      /* if we get here, j is within ROOT_S state 0's band */

      /* first initialize on d = 0 */
      alpha[0][0][0] = IMPOSSIBLE;
      y = cm->cfirst[v];
      alpha[0][0][0] = cm->endsc[v];
      /* treat EL as emitting only on self transition */
      for (yoffset = 0; yoffset < cm->cnum[0]; yoffset++)
	if ((sc = alpha[y+yoffset][0][0] + cm->tsc[0][yoffset]) > alpha[0][0][0]) 
	  alpha[0][0][0] = sc;
      /* ...we don't bother to look at local alignment starts here... */
      bestr[0] = -1;
      if (alpha[0][0][0] < IMPOSSIBLE) alpha[0][0][0] = IMPOSSIBLE;	
      alpha[0][1][0] = alpha[0][0][0];
      /* done initialization on d = 0 */

      jp_v = j - jmin[0];
      /* Impose the bands.
       *   We have to do this here because d bands are
       *   dependent on v AND j. 
       */
      for (d =0; d < hdmin[0][jp_v] && d <=W; d++) 
	alpha[0][cur][d] = IMPOSSIBLE;
      for (d = hdmax[0][jp_v]+1; d <= W;      d++) 
	alpha[0][cur][d] = IMPOSSIBLE;
      
      for (d = hdmin[0][jp_v]; ((d <= hdmax[0][jp_v] && d <= gamma_j) && d <= W); d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    {
	      if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
		{
		  alpha[0][cur][d] = sc;
		}
	    }
	  /*printf("j: %d | alpha[0][cur][%d]: %f\n", j, d, alpha[0][cur][d]);*/
	  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	}

      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(j >= jmin[y] && j <= jmax[y]) 
	    {
	      jp_y = j - jmin[y];
	      tmp_dmin = (hdmin[y][jp_y] > hdmin[0][jp_v]) ? hdmin[y][jp_y] : hdmin[0][jp_v];
	      tmp_dmax = (hdmax[y][jp_y] < hdmax[0][jp_v]) ? hdmax[y][jp_y] : hdmax[0][jp_v];
	      if(tmp_dmax > j) tmp_dmax = j;
	      for (d = tmp_dmin; d <= tmp_dmax; d++)
		{
		  if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
		  else                       sc = alpha[y][cur][d]     + cm->beginsc[y]; /* BUG! alpha[y][cur][d] outside bands */
		  if (sc > alpha[0][cur][d]) {
		    alpha[0][cur][d] = sc;
		    bestr[d]         = y;
		  }
		  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
		  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
		}
	    }
	}
      }
	      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = hdmin[0][jp_v]; (d <= hdmax[0][jp_v] && d <= gamma_j) && d <= W; d++) 
	{
	  i       = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  assert(i > 0);
	  /*printf("v: %d gamma_i: %d d: %d\n", v, gamma_i, d);*/
	  /*printf("alpha[0][j:%3d][d:%3d]: %f\n", j, d, alpha[0][cur][d]);*/
	  sc = gamma[gamma_i-1] + alpha[0][cur][d];
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = alpha[0][cur][d]; 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   * be careful about only freeing our impossible deck, imp_row, once
   *****************************************************************/ 
for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) 
	  if(alpha_mem[v][j] != imp_row) free(alpha_mem[v][j]); 
	free(alpha[v]);
	free(alpha_mem[v]);
      } else {
	if(alpha_mem[v][0] != imp_row) free(alpha_mem[v][0]);
	if(alpha_mem[v][1] != imp_row) free(alpha_mem[v][1]); 
	free(alpha[v]);
	free(alpha_mem[v]);
      }
    }
  free(imp_row);
  free(alpha);
  free(alpha_mem);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  return best_score;
 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}


/* Function: iInsideBandedScan_jd() 
 * Date    : EPN, Fri Apr 27 09:32:38 2007
 *           
 *           
 * Purpose:  Identical to CYKBandedScan_jd(), but sums replaces maxes.
 *           Scan a (sub)sequence for matches to a covariance model, using HMM
 *           derived bands on the j and d dimensions. Intended for use on 
 *           subsequences with endpoints i and j, where i and j were determined
 *           using a HMM scan of the full sequence. This function then refines
 *           the positions of i and j, as well as deriving a CYK score that is
 *           more informative than an HMM based score. 
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *           Log sums are performed using scaled ints with ILogsum().
 *
 *           jmin, jmax set the state specific bounds on the j dimension. 0..v..cm->M-1.
 *           hdmin, hdmax set the state and j specific bounds on the d dimension, indexed
 *           [0..v..cm-M-1][0..(jmax[v]-jmin[v]+1)].
 *           
 *           The j band for v is jmin[v]..jmax[v], inclusive; that is,
 *           jmin[v] is the minimum allowed j for state v;
 *           jmax[v] is the maximum; 
 * 
 *           The d bands are v and j specific (in contrast to the a priori d bands
 *           which are only v specific), the d band for v and j is 
 *           hdmin[v][j-jmin[v]]..hdmax[v][j-jmin[v]] inclusive;
 *           hdmin[v][j-jmin[v]] is the minimum allowed d for state v and end point j
 *           hdmax[v][j-jmin[v]] is the maximum; 
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; i0..j0
 *           jmin      - minimum bound on j for state v; 0..M
 *           jmax      - maximum bound on j for state v; 0..M
 *           hdmin     - minimum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           hdmax     - maximum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results   - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float
iInsideBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
		    int j0, int W, float cutoff, search_results_t *results)
{
  int        status;
  int     ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int     ***alpha_mem;          /* pointers to original alpha memory */
  int       *imp_row;           /* an IMPOSSIBLE deck (full of -INFTY scores), 
				 * pointed to when j is outside j band for v */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int      gamma_j;             /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int      gamma_i;             /* i index in the gamma matrix */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp_roll;   	        /* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       tmp_dmin, tmp_dmax; /* temp variables for ensuring we stay within d bands within loops */
  int       tmp_kmin, tmp_kmax; /* temp vars for B_st's, min/max k values consistent with bands*/

  int      jp_v, jp_y, jp_w;    /* mem eff banded j index in states v, y, and z 
				 * jp_x = j-jmin[x] */
  int      L;                   /* length of subsequence (j0-i0+1) */
  int      x;
  int      tmp_y;
  int      tmp_alpha_w, tmp_alpha_y;
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  
  /* Contract checks */
  if((!(cm->search_opts & CM_SEARCH_NOQDB)) && (cm->dmin == NULL || cm->dmax == NULL))
    cm_Fail("ERROR in iInsideBandedScan_jd(), trying to use QDB, but cm->dmin, cm->dmax are NULL.\n");
  if(dsq == NULL)
    cm_Fail("ERROR in iInsideBandedScan_jd(), dsq is NULL.");

  if(!(cm->search_opts & CM_SEARCH_NOQDB)) /* we're doing qdb */
    combine_qdb_hmm_d_bands(cm, jmin, jmax, hdmin, hdmax);

  best_score     = -INFTY;
  best_neg_score = -INFTY;
  L = j0-i0+1;
  /*printf("in iInsideBandedScan_jd i0: %5d | j0: %5d | L: %5d | W: %5d\n", i0, j0, L, W);*/
  if (W > L) W = L; /* shouldn't look longer than seq length L */

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));*/
  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that unlike the other CYK scan functions, E memory is not shared: 
   * this is because the E deck will be different for different j values
   * due to the j bands. 
   * 
   *****************************************************************/
  ESL_ALLOC(alpha,     sizeof(int **) * cm->M);
  ESL_ALLOC(alpha_mem, sizeof(int **) * cm->M);
  /* we use alpha_mem to remember where each alpha row (alpha[v][j]) is
   * in case we've set alpha[v][cur] to imp_row (the precalc'ed -INFTY row)
   * in a prior iteration, and we are about to overwrite it, and we don't
   * want to overwrite our special -INFTY row.
   */
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v],     sizeof(int *)  * (W+1));
	ESL_ALLOC(alpha_mem[v], sizeof(int *)  * (W+1));
	for (j = 0; j <= W; j++)
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(int) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
    else 
      {
	ESL_ALLOC(alpha[v],     sizeof(int *) * 2);
	ESL_ALLOC(alpha_mem[v], sizeof(int *) * 2);
	for (j = 0; j < 2; j++) 
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(int) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc,sizeof(float) * (L+1));
  ESL_ALLOC(saver, sizeof(int)   * (L+1));

  /* Initialize the impossible deck, which we'll point to for 
   * j positions that are outside of the j band on v */
  ESL_ALLOC(imp_row, sizeof(int) * (W+1));
  for (d = 0; d <= W; d++) imp_row[d] = -INFTY;
    
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {

      gamma_j = j-i0+1;
      cur = (j-i0+1)%2; /* cur == 1 when j == i0 */
      prv = (j-i0)  %2; /* prv == 0 when j == i0 */

      /*****************************************************************
       * alpha initializations.
       * For the jd (HMM) banded strategy, we initialize inside the j loop,
       * because no cells are j-independent: for j's outside
       * the bands for a state v, should have ALL cells = -INFTY.
       *****************************************************************/ 
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* Check to see if we're outside the bounds on j */
	  if(j < jmin[v] || j > jmax[v])
	    {
	      if (cm->stid[v] == BEGL_S) 
		{
		  jp_roll = j % (W+1); 
		  for (d = 0; d <= W; d++) 
		    alpha[v][jp_roll][d] = -INFTY;
		}
	      else
		{
		  jp_roll = cur;
		  alpha[v][jp_roll] = imp_row;
		}
	      /* Special boundary case: have to initialize alpha[v][prv] also */
	      if (j == i0)
		alpha[v][prv] = imp_row;
	      /*for (d = 0; d <= W; d++) 
		alpha[v][prv][d] = -INFTY;*/

	      continue;
	    }
	  if(alpha[v][cur] == imp_row)
	    {
	      /* we don't want to overwrite imp_row */
	      alpha[v][cur] = alpha_mem[v][cur];
	      for (d = 0; d <= W; d++) 
		alpha[v][cur][d] = -INFTY;
	    }
	  
	  /* else we initialize on d = 0 */
	  alpha[v][cur][0] = -INFTY;
	  
	  if      (cm->sttype[v] == E_st)  alpha[v][cur][0] = 0;
	  else if (cm->sttype[v] == MP_st) alpha[v][cur][1] = alpha[v][prv][1] = -INFTY;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[v][cur][0] = cm->iendsc[v];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[v][cur][0] = ILogsum(alpha[v][cur][0], (alpha[y+yoffset][cur][0] 
							      + cm->itsc[v][yoffset]));
	      /* ...we don't bother to look at local alignment starts here... */
	      bestr[cur] = -1;
	      if (alpha[v][cur][0] < -INFTY) alpha[v][cur][0] = -INFTY;	
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* w is BEGL_S */
	      y = cm->cnum[v];   /* y is BEGR_S */
	      /* original line: 
	       * alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	       * we can't use that because alpha[w][0][0] and alpha[y][0][0] 
	       * may have been reset to -INFTY, so we recalculate what they
	       * should be (this is wasteful):
	       */
	      tmp_y = cm->cfirst[w];
	      tmp_alpha_w = cm->iendsc[w];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[w]; yoffset++)
		{
		  tmp_alpha_w = ILogsum(tmp_alpha_w, (alpha[tmp_y+yoffset][cur][0] 
						      + cm->itsc[w][yoffset]));
		}
	      tmp_y = cm->cfirst[y];
	      tmp_alpha_y = cm->iendsc[y];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[y]; yoffset++)
		tmp_alpha_y = ILogsum(tmp_alpha_y, (alpha[tmp_y+yoffset][cur][0]
						    + cm->itsc[y][yoffset]));
	      alpha[v][cur][0] = tmp_alpha_w + tmp_alpha_y;
	      /*printf("! alpha[v][j:%3d][d:%3d]: %f\n", v, j, 0, Scorify(alpha[v][cur][0]));*/
	      if (alpha[v][cur][0] < -INFTY) alpha[v][cur][0] = -INFTY;	
	    }
	  
	  /* Special boundary case: have to initialize alpha[v][prv] also */
	  if(j == i0) alpha[v][prv][0] = alpha[v][cur][0];

	  if (cm->stid[v] == BEGL_S) 
	    {
	      alpha[v][prv][0] = alpha[v][cur][0];
	      for (x = 2; x <= W; x++) 
		alpha[v][x][0] = alpha[v][cur][0];
	    }
	  /* done initialization */
	  
	  jp_v = j - jmin[v];
	  /* Impose the bands.
	   *   We have to do this inside the main loop because d bands are
	   *   dependent on v AND j. 
	   */
	  if (cm->stid[v] == BEGL_S) jp_roll = j % (W+1); else jp_roll = cur;

	  for (d =0; d < hdmin[v][jp_v] && d <=W; d++) 
	    alpha[v][jp_roll][d] = -INFTY;
	  for (d = hdmax[v][jp_v]+1; d <= W;      d++) 
	    alpha[v][jp_roll][d] = -INFTY;
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp_roll][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][jp_roll][d] = ILogsum(alpha[v][jp_roll][d], (alpha[y+yoffset][cur][d] + 
									  cm->itsc[v][yoffset]));
		  if (alpha[v][jp_roll][d] < -INFTY) alpha[v][jp_roll][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-2] + 
								  cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][cur][d-1] + 
								  cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		  /*printf("alpha[v][j:%3d][d:%3d]: %f\n", v, j, d, Scorify(alpha[v][cur][d]));*/

		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-1] + 
								  cm->itsc[v][yoffset]));
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      /* Five inequalities must be satisfied to ensure that j and k 
	       * and k combinations correspond with alpha cells within the bands 
	       * on states y and w. 
	       * Below: jp_y = j - jmin[y] & jp_w = j - jmin[w]
	       *
	       * (1) j   >= jmin[y]          && j   <= jmax[y]
	       * (2) j-k >= jmin[w]          && j-k <= jmax[w]
	       * (3) k   >= hdmin[y][jp_y]   && k   <= hdmax[y][jp_y]
	       * (4) d-k >= hdmin[w][jp_w-k] && d-k <= hdmax[w][jp_w-k]
	       * (5) d   >= hdmin[v][jp_v]   && d   <= hdmax[v][jp_v]
	       */

	      /* initialize with endsc for all valid d for B state v */
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		{
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		}
	      /* Following code is careful, and not 'efficient' */
	      if(j >= jmin[y] && j <= jmax[y]) /* ensures (1): that j is valid for state y */
		{
		  jp_y = j - jmin[y];
		  jp_w = j - jmin[w]; 
		  i = j-d+1;
		  /*
		    printf("valid j: %d | jp_y: %d | jp_w: %d\n", j, jp_y, jp_w);
		    printf("hdmin[v][jp_v]: %d | hdmin[y][jp_y]: %d\n", hdmin[v][jp_v], hdmin[y][jp_y]);
		    printf("hdmax[v][jp_v]: %d | hdmax[y][jp_y]: %d\n", hdmax[v][jp_v], hdmax[y][jp_y]);
		  */
		  for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		    {
		      /* k is the length of the right fragment */
		      tmp_kmin = ((j-jmax[w]) > hdmin[y][jp_y]) ? (j-jmax[w]) : hdmin[y][jp_y];
		      if(tmp_kmin < 0) tmp_kmin = 0;
		      if(tmp_kmin < d-hdmax[w][jp_w-tmp_kmin]) tmp_kmin = d-hdmax[w][jp_w-tmp_kmin];
		      /* tmp_kmin is now smallest k that satisfies (2), (3), and (4) */

		      tmp_kmax = ((j-jmin[w]) < hdmax[y][jp_y]) ? (j-jmin[w]) : hdmax[y][jp_y];
		      if(tmp_kmax > d-hdmin[w][jp_w-tmp_kmax]) tmp_kmax = d-hdmin[w][jp_w-tmp_kmax];
		      /* tmp_kmax is now largest k that satisfies (2), (3), and (4) */
		      /*printf("tmp_kmin: %d | tmp_kmax: %d\n", tmp_kmin, tmp_kmax);*/
		      for (k = tmp_kmin; k <= tmp_kmax; k++)
			{
			  jp_roll = (j-k)%(W+1); /* jp_roll is rolling index into BEGL_S (state w) 
						  * deck j dimension */
			  alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[w][jp_roll][d-k] + 
									alpha[y][cur][k]));
			}
		      if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		      /*printf("B alpha[%d][%d][%d]: %f\n", v, cur, d, alpha[v][cur][d]);*/
		    }
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * IMPOSSIBLE), which is probably a little too fragile of a method. 
       */
      /* Check to see if we're within bounds on j */
      if(j < jmin[0] || j > jmax[0])
	{
	  /*printf("j: %d (gamma_j: %d) -INFTY BABY! min: %d max: %d\n", j, gamma_j, jmin[0], jmax[0]);*/
	  for (d = 0; d <= W; d++) 
	    alpha[0][jp_roll][d] = -INFTY; 
	  /* Inform the little semi-Markov model that deals with multihit parsing
	   * that a hit is impossible, j is outside root band on j:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  continue;
	}
      /* if we get here, j is within ROOT_S state 0's band */

      /* first initialize on d = 0 */
      alpha[0][0][0] = -INFTY;
      y = cm->cfirst[v];
      alpha[0][0][0] = cm->iendsc[v];
      /* treat EL as emitting only on self transition */
      for (yoffset = 0; yoffset < cm->cnum[0]; yoffset++)
	alpha[0][0][0] = ILogsum(alpha[0][0][0], (alpha[y+yoffset][0][0]
						  + cm->itsc[0][yoffset]));
      /* ...we don't bother to look at local alignment starts here... */
      bestr[0] = -1;
      if (alpha[0][0][0] < -INFTY) alpha[0][0][0] = -INFTY;	
      alpha[0][1][0] = alpha[0][0][0];
      /* done initialization on d = 0 */

      jp_v = j - jmin[0];
      /* Impose the bands.
       *   We have to do this here because d bands are
       *   dependent on v AND j. 
       */
      for (d =0; d < hdmin[0][jp_v] && d <=W; d++) 
	alpha[0][cur][d] = -INFTY;
      for (d = hdmax[0][jp_v]+1; d <= W;      d++) 
	alpha[0][cur][d] = -INFTY;
      
      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->itsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    {
	      alpha[0][cur][d] = ILogsum(alpha[0][cur][d], (alpha[y+yoffset][cur][d]
							    + cm->itsc[0][yoffset]));
	    }
	  /*printf("j: %d | alpha[0][cur][%d]: %f\n", j, d, alpha[0][cur][d]);*/
	  if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
	  if (Scorify(alpha[0][cur][d]) > best_neg_score) best_neg_score = Scorify(alpha[0][cur][d]);
	}

      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(j >= jmin[y] && j <= jmax[y]) 
	    {
	      jp_y = j - jmin[y];
	      tmp_dmin = (hdmin[y][jp_y] > hdmin[0][jp_v]) ? hdmin[y][jp_y] : hdmin[0][jp_v];
	      tmp_dmax = (hdmax[y][jp_y] < hdmax[0][jp_v]) ? hdmax[y][jp_y] : hdmax[0][jp_v];
	      if(tmp_dmax > j) tmp_dmax = j;
	      for (d = tmp_dmin; d <= tmp_dmax; d++)
		{
		  if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->ibeginsc[y];
		  else                       sc = alpha[y][cur][d]     + cm->ibeginsc[y];
		  if (sc > alpha[0][cur][d]) {
		    alpha[0][cur][d] = sc;
		    bestr[d]         = y;
		  }
		  if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
		  if (Scorify(alpha[0][cur][d]) > best_neg_score) best_neg_score = Scorify(alpha[0][cur][d]);
		}
	    }
	}
      }
	      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = hdmin[0][jp_v]; (d <= hdmax[0][jp_v] && d <= gamma_j) && d <= W; d++) 
	{
	  i       = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  assert(i > 0);
	  /*printf("v: %d gamma_i: %d d: %d\n", v, gamma_i, d);*/
	  /*printf("alpha[0][j:%3d][d:%3d]: %f\n", j, d, Scorify(alpha[0][cur][d]));*/
	  sc = gamma[gamma_i-1] + Scorify(alpha[0][cur][d]);
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = Scorify(alpha[0][cur][d]); 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */
  
  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   * be careful about only freeing our impossible deck, imp_row, once
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) 
	  if(alpha[v][j] != imp_row) free(alpha[v][j]); 
	free(alpha[v]);
      } else {
	if(alpha[v][0] != imp_row) free(alpha[v][0]);
	if(alpha[v][1] != imp_row) free(alpha[v][1]); 
	free(alpha[v]);
      }
    }
  free(imp_row);
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  return best_score;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* Below is OldActuallySearchTarget(), the 0.81 version of cm_dpsearch.c:ActuallySearchTarget(),
 * it's deprecated, and can't be compiled because it uses old flags for the CM. Kept here for
 * reference.
 */
#if 0
/* 
 * Function: OldActuallySearchTarget()
 * Incept:   EPN, Mon Jan  8 06:42:59 2007
 *
 * Purpose:  Given a CM and a sequence, call the correct search algorithm
 *           based on cm->search_opts. Uses v0.81 (old) DP functions, as
 *           opposed to the fast v1.0 (newer) DP functions called by 
 *           ActuallySearchTarget().
 * 
 * Args:     cm              - the covariance model
 *           dsq             - the target sequence (digitized)
 *           i0              - start of target subsequence (often 1, beginning of dsq)
 *           j0              - end of target subsequence (often L, end of dsq)
 *           cm_cutoff       - minimum CM  score to report 
o *           cp9_cutoff      - minimum CP9 score to report (or keep if filtering)
 *           results         - search_results_t to keep results in, must be empty; if NULL, don't add to it
 *           do_filter       - TRUE if we should filter, but only if cm->search_opts tells us to 
 *           doing_cm_stats  - TRUE if the reason we're scanning this seq is to build
 *                             a histogram to calculate Gumbels for the CM, in this
 *                             case we don't filter regardless of what cm->search_opts says.
 *           doing_cp9_stats - TRUE if we're calc'ing stats for the CP9, in this 
 *                             case we always run CP9Forward()
 *           ret_flen        - RETURN: subseq len that survived filter (NULL if not filtering)
 *           do_align_hits   - TRUE to do alignments and return  parsetrees in results
 *
 * Returns: Highest scoring hit from search (even if below cutoff).
 */
float OldActuallySearchTarget(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, float cm_cutoff, 
			     float cp9_cutoff, search_results_t *results, int do_filter, 
			     int doing_cm_stats, int doing_cp9_stats, int *ret_flen,
			     int do_align_hits)
{
  int status;
  float sc;
  int flen;
  int use_cp9;    

  /*printf("in OldActuallySearchTarget: i0: %d j0: %d do_filter: %d doing_cm_stats: %d doing_cp9_stats: %d\n", i0, j0, do_filter, doing_cm_stats, doing_cp9_stats);
    printf("\ti0: %d j0: %d filter: %d\n", i0, j0, do_filter);*/

  /* Contract checks */
  if(dsq == NULL)                                 cm_Fail("OldActuallySearchTarget(): dsq is NULL.");
  if(!(cm->flags & CMH_BITS))                     cm_Fail("OldActuallySearchTarget(): CMH_BITS flag down.\n");
  if(doing_cm_stats && doing_cp9_stats)           cm_Fail("OldActuallySearchTarget(): doing_cm_stats and doing_cp9_stats both TRUE.\n");
  if(results != NULL && results->num_results > 0) cm_Fail("OldActuallySearchTarget(): there's already hits in results.\n");

  flen = (j0-i0+1);

  /* Check if we need the CP9 */
  use_cp9 = FALSE;
  /* use the CP9 b/c we're calcing CP9 Gumbel stats */
  if(doing_cp9_stats) use_cp9 = TRUE;                     
  /* use the CP9 b/c we're searching only with the CP9 HMM */
  if((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD)) use_cp9 = TRUE; 
  /* The third way we use the CP9 is if we're filtering, AND we haven't 
   * called this function recursively from AFTER filtering (the do_filter flag)
   * AND we're not determining CM Gumbel stats. */
  if((cm->search_opts & CM_SEARCH_HMMFILTER) && (do_filter && !doing_cm_stats)) use_cp9 = TRUE;

  /* Check if we have a valid CP9 (if we need it) */
  if(use_cp9) {
    if(cm->cp9 == NULL)                    cm_Fail("OldActuallySearchTarget(), trying to use CP9 HMM that is NULL\n");
    if(!(cm->cp9->flags & CPLAN9_HASBITS)) cm_Fail("OldActuallySearchTarget(), trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");
    if((cm->search_opts & CM_SEARCH_HBANDED) && cm->cp9b == NULL) cm_Fail("OldActuallySearchTarget(), trying to use CP9 HMM for HMM banded search, but cm->cp9b is NULL.\n");
  }      

  if(use_cp9)
    sc = CP9Scan_dispatch(cm, dsq, i0, j0, cm->W, cm_cutoff, cp9_cutoff, results, doing_cp9_stats, ret_flen);
  else {
      if(cm->search_opts & CM_SEARCH_HBANDED) {
	if((status = cp9_Seq2Bands(cm, NULL, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, i0, j0, cm->cp9b, TRUE, 0)) != eslOK) cm_Fail("OldActuallySearchTarget(): unrecoverable error in cp9_Seq2Bands().");

	/*debug_print_hmm_bands(stdout, (j0-i0+1), cm->cp9b, cm->tau, 3);*/
	if(cm->search_opts & CM_SEARCH_INSIDE)
	  sc = iInsideBandedScan_jd(cm, dsq, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, 
				    i0, j0, cm->W, cm_cutoff, results);
	else /* don't do inside */
	  sc = CYKBandedScan_jd(cm, dsq, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, 
				i0, j0, cm->W, cm_cutoff, results);
      }
      else if(cm->search_opts & CM_SEARCH_NOQDB) {
	if(cm->search_opts & CM_SEARCH_INSIDE)
	  sc = iInsideScan(cm, dsq, i0, j0, cm->W, cm_cutoff, results);
	else /* don't do inside */
	  sc = CYKScan (cm, dsq, i0, j0, cm->W, cm_cutoff, results);
	/* sc = FastCYKScan(cm, dsq, i0, j0, cm->W, cm_cutoff, results, NULL); */
      }
      else { /* use QDB */
	if(cm->search_opts & CM_SEARCH_INSIDE)
	  sc = iInsideBandedScan(cm, dsq, cm->dmin, cm->dmax, i0, j0, cm->W, cm_cutoff, results);
	else /* don't do inside */
	  sc = CYKBandedScan (cm, dsq, cm->dmin, cm->dmax, i0, j0, cm->W, cm_cutoff, results);
	/* sc = FastCYKScan(cm, dsq, i0, j0, cm->W, cm_cutoff, results, NULL);*/
      }
  }    
  if((results != NULL && results->num_results > 0) && do_align_hits) {
    if(cm->align_opts & CM_ALIGN_OLDDP) { 
      OldActuallyAlignTargets(cm, NULL, 
			      dsq, results,   /* put function into dsq_mode, designed for aligning search hits */
			      0, 0, 0, NULL);
    }
    else {
      ActuallyAlignTargets(cm, NULL, NULL, 
			   dsq, results,      /* put function into dsq_mode, designed for aligning search hits */
			   0, 0, 0, NULL);
    }
   }
  return sc;
}
#endif

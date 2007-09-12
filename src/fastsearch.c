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
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"



/* Function: FastCYKScan()
 * Date:     EPN, Wed Sep 12 16:55:28 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm. If bands are NULL, reverts to non-banded
 *           (scancyk.c:CYKScan()). 
 *
 * Args:     cm              - the covariance model
 *           dsq             - the digitized sequence
 *           dmin            - minimum bound on d for state v; 0..M
 *           dmax            - maximum bound on d for state v; 0..M          
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *           ret_best_hit_sc - RETURN score of best hit (reported to results) NULL if not-wanted
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
FastCYKScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float cutoff, 
	    search_results_t *results, float **ret_vsc, float *ret_best_hit_sc)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  float    *vsc;                /* best score for each state (float) */
  float     vsc_root;           /* best overall score (score at ROOT_S) */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_w;  	        /* offset j for state w */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       ip_g;               /* offset i for gamma (i-i0+1) */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       dn;                 /* temporary value for min d in for loops */
  int       dx;                 /* temporary value for max d in for loops */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  float     best_hit_sc;        /* best hit score found */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  float     sc;                 /* temp score for reporting hits */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in FastCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in FastCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in FastCYKScan, CM_SEARCH_INSIDE flag raised");

  /* determine if we're doing banded/non-banded */
  if(dmin != NULL && dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 

  vsc = NULL;
  if(ret_vsc != NULL) { 
    ESL_ALLOC(vsc, sizeof(float) * cm->M);
    esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);
  }
  best_hit_sc = IMPOSSIBLE;
  vsc_root    = IMPOSSIBLE;

  /*
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   */
  ESL_ALLOC(alpha, (sizeof(float **) * cm->M));
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S) {
	ESL_ALLOC(alpha[v], (sizeof(float *) * (W+1)));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) alpha[v] = alpha[cm->M-1];
    else {
      ESL_ALLOC(alpha[v], sizeof(float *) * 2);
      for (j = 0; j < 2; j++) ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
    }
  }
  ESL_ALLOC(bestr, (sizeof(int) * (W+1)));

  /*
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
   */
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
	    alpha[v][0][0] = ESL_MAX(alpha[v][0][0], (alpha[y+yoffset][0][0] + cm->tsc[v][yoffset]));
	  /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  alpha[v][0][0] = ESL_MAX(alpha[v][0][0], IMPOSSIBLE);
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

  /*
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   */
  if(results != NULL) { 
    ESL_ALLOC(gamma,  sizeof(float) * (L+1));
    gamma[0] = 0;
    ESL_ALLOC(gback,  sizeof(int)   * (L+1));
    gback[0] = -1;
    ESL_ALLOC(savesc, sizeof(float) * (L+1));
    ESL_ALLOC(saver,  sizeof(int)   * (L+1));
  }
  /*
   * query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) {
      jmax = (cm->stid[v] == BEGL_S) ? W : 1;
      for (d = 0; d < dmin[v] && d <=W; d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
      for (d = dmax[v]+1; d <= W;      d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
    }
  }

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* determine min/max d we're allowing for this state v and this position j */
	  if(do_banded) { 
	    dn = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	    dx = ESL_MIN(jp_g, dmax[v]); 
	    dx = ESL_MIN(dx, W);
	  }
	  else { 
	    dn = (cm->sttype[v] == MP_st) ? 2 : 1;
	    dx = ESL_MIN(jp_g, W); 
	  }

	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  sd   = StateDelta(cm->sttype[v]);

	  if(cm->sttype[v] == B_st) {
	    w = cm->cfirst[v];
	    y = cm->cnum[v];
	    for (d = dn; d <= dx; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	      }
	      else { kmin = 0; kmax = d; }

	      alpha[v][jp_v][d] = ESL_MAX(IMPOSSIBLE, cm->endsc[v] + (cm->el_selfsc * (d - sd)));
	      for (k = kmin; k <= kmax; k++) { 
		jp_w = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		alpha[v][jp_v][d] = ESL_MAX(alpha[v][jp_v][d], (alpha[w][jp_w][d-k] + alpha[y][jp_y][k]));
	      }
	    }
	  }
	  else { /* if cm->sttype[v] != B_st */
	    for (d = dn; d <= dx; d++) {
	      alpha[v][jp_v][d] = ESL_MAX (IMPOSSIBLE, (cm->endsc[v] + (cm->el_selfsc * (d - sd))));
	      y = cm->cfirst[v];
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[v][jp_v][d] = ESL_MAX (alpha[v][jp_v][d], (alpha[y+yoffset][jp_y][d - sd] + cm->tsc[v][yoffset]));
		
	      /* add in emission score, if any */
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st: 
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][jp_v][d] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		break;
	      case ML_st:
	      case IL_st:
		alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		break;
	      case MR_st:
	      case IR_st:
		alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		break;
	      } /* end of switch */
	    } /* end of d = dn; d <= dx; d++ for B_st */
	  }
	  /* end of else (v != B_st) */
	  if(vsc != NULL) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[v][jp_v][d]);
	} /*loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      /* determine min/max d we're allowing for the root state and this position j */
      if(do_banded) { 
	dn = ESL_MAX(dmin[0], 1); 
	dx = ESL_MIN(jp_g, dmax[0]); 
	dx = ESL_MIN(dx, W);
      }
      else { 
	dn = 1; 
	dx = ESL_MIN(jp_g, W); 
      }
      jp_v = cur;
      jp_y = cur;
      for (d = dn; d <= dx; d++) {
	bestr[d] = 0;	/* root of the traceback = root state 0 */
	y = cm->cfirst[0];
	alpha[0][jp_v][d] = ESL_MAX(IMPOSSIBLE, alpha[y][cur][d] + cm->tsc[0][0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[0][jp_v][d] = ESL_MAX (alpha[0][jp_v][d], (alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]));
      }
	
      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(do_banded) {
	    dn = (cm->sttype[y] == MP_st) ? ESL_MAX(dmin[y], 2) : ESL_MAX(dmin[y], 1); 
	    dn = ESL_MAX(dn, dmin[0]);
	    dx = ESL_MIN(jp_g, dmax[y]); 
	    dx = ESL_MIN(dx, W);
	  }
	  else { 
	    dn = 1; 
	    dx = ESL_MIN(jp_g, W); 
	  }
	  jp_y = (cm->stid[y] == BEGL_S) ? (j % (W+1)) : cur;
	  for (d = dn; d <= dx; d++) {
	    {
	      /* Is this more efficient:? 
		bestr[d]          = (alpha[0][jp_v][d] > (alpha[y][jp_y][d] + cm->beginsc[y])) ? bestr[d] : y;
		alpha[0][jp_v][d] = ESL_MAX(alpha[0][jp_v][d], alpha[y][jp_y][d] + cm->beginsc[y]); */
	      if(alpha[0][jp_v][d] < (alpha[y][jp_y][d] + cm->beginsc[y])) {
		alpha[0][jp_v][d] = alpha[y][jp_y][d] + cm->beginsc[y];
		bestr[d] = y;
	      }
	    }
	  }
	}
      }
      /* find the best score (we need dn and dx for state 0 again) */
      if(do_banded) { 
	dn = ESL_MAX(dmin[0], 1); 
	dx = ESL_MIN(jp_g, dmax[0]); 
	dx = ESL_MIN(dx, W);
      }
      else { 
	dn = 1; 
	dx = ESL_MIN(jp_g, W); 
      }
      for (d = dn; d <= dx; d++) 
	vsc_root = ESL_MAX(vsc_root, alpha[0][jp_v][d]);


      if(results != NULL) 
	{
	  /* get information on hits */
	  if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	    {
	      /* The little semi-Markov model that deals with multihit parsing:
	       */
	      gamma[jp_g]  = gamma[jp_g-1] + 0; 
	      gback[jp_g]  = -1;
	      savesc[jp_g] = IMPOSSIBLE;
	      saver[jp_g]  = -1;
	      for (d = dn; d <= dx; d++) {
		i    = j-d+1;
	      ip_g = i-i0+1;
	      sc = gamma[ip_g-1] + alpha[0][jp_v][d] + cm->sc_boost; 
	      /* sc_boost is experimental technique for finding hits < 0 bits. value is 0.0 if technique not used. */
	      if (sc > gamma[jp_g]) {
		gamma[jp_g]  = sc;
		gback[jp_g]  = i;
		savesc[jp_g] = alpha[0][jp_v][d]; 
		saver[jp_g]  = bestr[d];
	      }
	      }
	    }
	  else {
	    /* Resolving overlaps greedily (RSEARCH style),  
	     * At least one hit is sent back for each j here.
	     * However, some hits can already be removed for the greedy overlap
	     * resolution algorithm.  Specifically, at the given j, any hit with a
	     * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
	    /* First, report hit with d of 1 if > cutoff */
	    if (alpha[0][jp_v][1] >= cutoff) 
	      if(results != NULL) 
		report_hit (j, j, bestr[1], alpha[0][jp_v][1], results);
	    bestd = 1;
	    if (alpha[0][jp_v][1] > best_hit_sc)
	      best_hit_sc = alpha[0][jp_v][1];
	    
	    /* Now, if current score is greater than maximum seen previous, report
	       it if >= cutoff and set new max */
	    for (d = 2; d <= W && d <= jp_g; d++) 
	      {
		if (alpha[0][jp_v][d] > best_hit_sc) best_hit_sc = alpha[0][jp_v][d];
		if (alpha[0][jp_v][d] > alpha[0][jp_v][bestd]) {
		  if (alpha[0][jp_v][d] >= cutoff)
		    if(results != NULL) 
		      report_hit (j-d+1, j, bestr[d], alpha[0][jp_v][d], results);
		  bestd = d;
		}
	      }
	  }
	}
      } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* free alpha, everything we need is in gamma.
   */
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

  if(results != NULL) 
    {
      /*
       * Traceback stage.
       * Recover all hits: an (i,j,sc) triple for each one.
       */
      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	{
	  j     = j0;
	  while (j >= i0) 
	    {
	      jp_g = j-i0+1;
	      if (gback[jp_g] == -1) /* no hit */
		j--; 
	      else                /* a hit, a palpable hit */
		{
		  if(savesc[jp_g] > best_hit_sc) best_hit_sc = savesc[jp_g];
		  if(savesc[jp_g] >= cutoff && results != NULL) /* report the hit */
		    report_hit(gback[jp_g], j, saver[jp_g], savesc[jp_g], results);
		  j = gback[jp_g]-1;
		}
	    }
	}
      free(gback);
      free(gamma);
      free(savesc);
      free(saver);
    }
  
  if (ret_best_hit_sc != NULL) *ret_best_hit_sc = best_hit_sc;
  if (ret_vsc         != NULL) *ret_vsc         = vsc;
  
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

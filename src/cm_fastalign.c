/* cm_fastalign.c
 * EPN, Wed Oct 10 07:20:48 2007
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
#include <assert.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

#define TSC(s,k) (tsc[(v) * MAXCONNECT + (s)])
#define AMX(j,v,d) (alphap[(j * cm->M * (W+1)) + ((v) * (W+1) + d)])



 
/* EPN 03.29.06
 * Function: Fast_inside_b_jd_me()
 * based on inside_b_me() which was ...
 * based on inside()
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phase of a CYK alignment algorithm
 *           using bands in the j and d dimension from obtained
 *           from an HMM forwards-backwards run. This function
 *           is memory efficient in the j AND d dimension.
 * 
 *           To be able to consistently handle end states, the
 *           original SRE behavior of reusing the end deck was
 *           abandoned. Now each end state has its own deck, which
 *           makes this implementation easier because each state
 *           has its own bands on j, and thus has a state specific
 *           offset with alpha[end][jp][dp] in the banded mem eff 
 *           matrix corresponding to alpha[end][jp+jmin[end]][dp+hdmin[v][jp_v]]
 *           in the platonic matrix.
 *           
 *           The deck re-use strategy in general does not work with
 *           this implementation b/c each state has it's own j-specific
 *           bands. 
 *
 *           Notes from inside():
 *           A note on the loop conventions. We're going to keep the
 *           sequence (dsq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: only vroot..vend
 *           in the decks, i0-1..j in the rows, and up to j0-i0+1 in
 *           the columns (d dimension). Where this is handled the most
 *           is in two variables: W, which is the length of the subsequence
 *           (j0-i0+1), and is oft used in place of L in the usual CYK;
 *           and jp (read: j'), which is the *relative* j w.r.t. the
 *           subsequence, ranging from 0..W, and then d ranges from 
 *           0 to jp, and j is calculated from jp (i0-1+jp).
 *           In this banded version, there are more offset issues,
 *           these are detailed with comments in the code.
 *
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in an
 *           alpha, we expect that alpha[vroot..vend] are all NULL
 *           decks already; any other decks <vroot and >vend will
 *           be preserved. If we pass in a dpool, the decks *must* be
 *           sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *           We also deal with local begins, by keeping track of the optimal
 *           state that we could enter and account for the whole target 
 *           sequence: b = argmax_v  alpha_v(i0,j0) + log t_0(v),
 *           and bsc is the score for that. 
 *
 *           If vroot==0, i0==1, and j0==L (e.g. a complete alignment),
 *           the optimal alignment might use a local begin transition, 0->b,
 *           and we'd have to be able to trace that back. For any
 *           problem where the caller sets allow_begin, we return a valid b 
 *           (the optimal 0->b choice) and bsc (the score if 0->b is used).
 *           If a local begin is part of the optimal parse tree, the optimal
 *           alignment score returned by inside() will be bsc and yshad[0][L][L] 
 *           will be USE_LOCAL_BEGIN, telling insideT() to check b and
 *           start with a local 0->b entry transition. When inside()
 *           is called on smaller subproblems (v != 0 || i0 > 1 || j0
 *           < L), we're using inside() as an engine in divide &
 *           conquer, and we don't use the overall return score nor
 *           shadow matrices, but we do need allow_begin, b, and bsc for
 *           divide&conquer to sort out where a local begin might be used.
 *
 * Args:     cm        - the model    [0..M-1]
 *           sq        - the sequence [1..L]   
 *                     - length of the dsq
 *           vroot     - first start state of subtree (0, for whole model)
 *           vend      - last end state of subtree (cm->M-1, for whole model)
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will 
 *                       be valid. 
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           cp9b      - HMM bands for <dsq> and <cm>
 *                       
 * Returns: Score of the optimal alignment.  
 */
float 
Fast_inside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
		    CM_FHB_MX *mx, void ****ret_shadow, 
		    int allow_begin, int *ret_b, float *ret_bsc, CP9Bands_t *cp9b)
{
  /* Contract check */
  if(dsq == NULL) cm_Fail("inside_b_jd_me(), dsq is NULL.\n");
  if (mx == NULL) cm_Fail("inside_b_jd_me(), mx is NULL.\n");

  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  void  ***shadow;      /* shadow matrix for tracebacks */
  int    **kshad;       /* a shadow deck for bifurcations */
  char   **yshad;       /* a shadow deck for every other kind of state */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */

  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;

  int       *jmin = cp9b->jmin;
  int       *jmax = cp9b->jmax;
  int     **hdmin = cp9b->hdmin;
  int     **hdmax = cp9b->hdmax;

  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */

  /* grow the matrix based on the current sequence and bands */
  cm_fhb_mx_GrowTo(mx, cp9b);
  float    ***alpha = mx->dp;

  /* The shadow matrix, if caller wants a traceback.
   * We do some pointer tricks here to save memory. The shadow matrix
   * is a void ***. Decks may either be char ** (usually) or
   * int ** (for bifurcation decks). Watch out for the casts.
   * For most states we only need
   * to keep y as traceback info, and y <= 6. For bifurcations,
   * we need to keep k, and k <= L, and L might be fairly big.
   * (We could probably limit k to an unsigned short ... anyone
   * aligning an RNA > 65536 would need a big computer... but
   * we'll hold off on that for now. We could also pack more
   * traceback pointers into a smaller space since we only really
   * need 3 bits, not 8.)
   */
  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  float *el_scA;
  float **init_scAA;
  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* precalculate the initial score for all alpha[v][j][d] cells, it's independent
   * of j, so we do it here, outside the for(j...) loop */
  ESL_ALLOC(init_scAA, sizeof(float *) * (cm->M));
  for (v = 0; v < cm->M; v++) {
      ESL_ALLOC(init_scAA[v], sizeof(float) * (W+1));

      if(NOT_IMPOSSIBLE(cm->endsc[v]))
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = el_scA[d] + cm->endsc[v];
      else
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = IMPOSSIBLE;
  }

  int *yvalidA; 
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, 0);

  /* Main recursion
   */
  for (v = vend; v >= vroot; v--) 
    {
      /* First we need a shadow deck to fill in. 
       */
      if (cm->sttype[v] != E_st) {
	if (ret_shadow != NULL) {
	  if (cm->sttype[v] == B_st) {
	    kshad     = alloc_jdbanded_vjd_kshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	    shadow[v] = (void **) kshad;
	  } else {
	    yshad     = alloc_jdbanded_vjd_yshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	    shadow[v] = (void **) yshad;
	  }
	}
      }

      /* We've only allocated alpha cells that are within the bands
       * on the j and d dimensions. This means we have to deal
       * with all sorts of offset issues, but we don't have to 
       * waste time setting cells outside the bands to IMPOSSIBLE.
       */
      int sd   = StateDelta(cm->sttype[v]);
      int sdr  = StateRightDelta(cm->sttype[v]);
      int jn   = jmin[v];
      int jx   = jmax[v];
      int dn;
      int dx;
      float tsc;
      int cnum = cm->cnum[v];
      int jp_y_sdr;
      int dp_y_sd;
      int dpn;
      int dpx;
      int jpn;
      int jpx;
      int yvalid_idx;
      int yvalid_ct;
      int j_sdr;
      if(cm->sttype[v] == E_st) { 
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v = j-jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	    ESL_DASSERT1((hdmin[v][jp_v] == 0));
	    ESL_DASSERT1((hdmax[v][jp_v] == 0));
	    alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	  }
	}
	continue;
      }
      //else if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      else if(cm->sttype[v] == 100) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	    alpha[v][jp_v][dp_v]  = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
	    /* treat EL as emitting only on self transition */
	    if (ret_shadow != NULL) yshad[jp_v][dp_v]  = USED_EL; 
	    for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
		yoffset = y - cm->cfirst[v];

		if((j-sdr) >= jmin[y] && ((j-sdr) <= jmax[y])) /* Enforces j-sdr is valid for state y */
		  {
		    jp_y = j - sdr - jmin[y];

		    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
		    ESL_DASSERT1((jp_y >= 0 && jp_y <= (jmax[y]-jmin[y])));
		    /* if(cm->sttype[v] == MP_st) printf("v: %5d j: %5d jp_y_sdr: %5d\n", v, j, jp_y); */
		
		    if((d-sd) >= hdmin[y][jp_y] && (d-sd) <= hdmax[y][jp_y])
		      {
			dp_y = d - sd - hdmin[y][jp_y];  /* d index for state y 
							    in alpha w/mem eff bands */
			/* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			 * corresponding to alpha[y][j][d] in the platonic matrix.
			 */
			ESL_DASSERT1((dp_v >= 0 && dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v])));
			ESL_DASSERT1((dp_y >= 0 && dp_y <= (hdmax[y][jp_y] - hdmin[y][jp_y])));
			if ((sc = alpha[y][jp_y][dp_y] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
			  {
			    alpha[v][jp_v][dp_v] = sc; 
			    if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
			  }
		      }
		  }
	    }
	    switch(cm->sttype[v]) { 
	    case ML_st:
	    case IL_st:
	      i = j-d+1;
	      if (dsq[i] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[i]];
	      else                     alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	      break;
	    case MR_st:
	    case IR_st:
	      if (dsq[j] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[j]];
	      else             	       alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	      break;
	    case MP_st:
	      i = j-d+1;
	      if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
	      else  		                              alpha[v][jp_v][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
	      break;
	    }
	    if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	  }
	}
      }
      //else if(cm->sttype[v] == 100) { /* IL_st || cm->sttype[v] == IR_st) {*/ 
      else if(cm->sttype[v] == IL_st) {
	/* initialize all cells within v's j band */
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	    alpha[v][jp_v][dp_v] = init_scAA[v][d-sd];
	    yshad[jp_v][dp_v] = USED_EL; 
	  }
	}
	/* update alpha[v][jp_v][dp_v] cells, for IL, IR states, loop nesting order is:
	 * for j { for d { for y { } } } because they can self transit, and a 
	 * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
	 * before can start calc'ing for alpha[v][j][d+1] */
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  j_sdr = j - sdr;

	  /* determine which children y we can legally transit to for v, j */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */

	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    i = j - d + 1;
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      jp_y_sdr = j - jmin[y] - sdr;

	      if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		  if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
		    {
		      alpha[v][jp_v][dp_v] = sc; 
		      yshad[jp_v][dp_v]    = yoffset;
		    }
	      }
	    }
	    if (dsq[i] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[i]];
	    else                     alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	    i--;
	  }
	}
      }
      else if(cm->sttype[v] == IR_st) { 
	/* initialize all cells within v's j band */
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	    alpha[v][jp_v][dp_v] = init_scAA[v][d-sd];
	    yshad[jp_v][dp_v] = USED_EL; 
	  }
	}
	/* update alpha[v][jp_v][dp_v] cells, for IL, IR states, loop nesting order is:
	 * for j { for d { for y { } } } because they can self transit, and a 
	 * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
	 * before can start calc'ing for alpha[v][j][d+1] */
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  j_sdr = j - sdr;

	  /* determine which children y we can legally transit to for v, j */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */

	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      jp_y_sdr = j - jmin[y] - sdr;

	      if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		  if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
		    {
		      alpha[v][jp_v][dp_v] = sc; 
		      yshad[jp_v][dp_v]    = yoffset;
		    }
	      }
	    }
	    if (dsq[j] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[j]];
	    else          	       alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	  }
	}
      }
      else if(cm->sttype[v] != B_st) {
	/* initialize all cells within v's j band */
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	    alpha[v][jp_v][dp_v] = init_scAA[v][d-sd];
	    yshad[jp_v][dp_v] = USED_EL; 
	  }
	}
	/* for each child state y of v, update alpha[v][jp_v][dp_v] cells */
	for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	  yoffset = y - cm->cfirst[v];
	  tsc = cm->tsc[v][yoffset];

	  jn = ESL_MAX(jmin[v], ESL_MIN(jmin[y] + sdr, jmax[y]));
	  jx = ESL_MIN(jmax[v], ESL_MAX(jmax[y] + sdr, jmin[y]));
	  jx = ESL_MIN(jx, jmax[y] + sdr);
	  jpn = jn - jmin[v];
	  jpx = jx - jmin[v];
	  jp_y_sdr = jn - jmin[y] - sdr;

	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));

	    dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	    dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	    dx = ESL_MIN(dx, hdmax[y][jp_y_sdr] + sd);
	    dpn     = dn - hdmin[v][jp_v];
	    dpx     = dx - hdmin[v][jp_v];
	    dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;

	    for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	      if((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc) > alpha[v][jp_v][dp_v]) {
		alpha[v][jp_v][dp_v] = sc;
		yshad[jp_v][dp_v]    = yoffset;
	      }
	    }
	  }
	}
	/* add in emission score, if any */
	switch(cm->sttype[v]) { 
	case ML_st:
	case IL_st:
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    i     = j - hdmin[v][jp_v] + 1;
	    for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++, i--) {
	      if (dsq[i] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[i]];
	      else                     alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	    }
	  }
	  break;
	case MR_st:
	case IR_st:
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	      if (dsq[j] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[j]];
	      else          	     alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	    }
	  }
	  break;
	case MP_st:
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    i     = j - hdmin[v][jp_v] + 1;
	    for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++, i--) {
	      if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K) alpha[v][jp_v][dp_v] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
	      else  		                              alpha[v][jp_v][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
	    }
	  }
	}
	/*for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++, i--) {
	    printf("alpha[v: %4d][jp_v: %4d][dp_v: %4d]: %.4f\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v]);

	  }
	  printf("\n");
	}
	printf("\n\n");*/
      }
      else { /* B_st */
	y = cm->cfirst[v];
	z = cm->cnum[v];
	/* Any valid j must be within both state v and state z's j band 
	 * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
	 * constructed, but we'll check anyway. 
	 */
	tmp_jmin = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	tmp_jmax = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];

	/* For any values of j within v's j band but outside of z's j band,
	 * we have to set the corresponding alpha cells to IMPOSSIBLE.
	 * This is done be the following two ugly for loops, 
	 * which will only be looked at once for each B state, and
	 * even then only *very* rarely entered. This
	 * is why they're here, seemingly out of place before the 
	 * main j loop below, where similar performing code would be 
	 * looked at on the order of j times, instead of just once.
	 */
	for(j = jmin[v]; j < tmp_jmin; j++)
	  {
	    jp_v = j-jmin[v];
	    for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d-hdmin[v][jp_v];
		alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
	      }
	  }
	if(tmp_jmax < jmax[v])
	  for(j = (tmp_jmax+1); j <= jmax[v]; j++)
	    {
	      jp_v = j-jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d-hdmin[v][jp_v];
		  alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
		}
	    }
	/* the main j loop */
	for (j = tmp_jmin; j <= tmp_jmax; j++)
	  {
	    jp_v = j - jmin[v];
	    jp_y = j - jmin[y];
	    jp_z = j - jmin[z];
	    for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */

		/* Find the first k value that implies a valid cell in the y and z decks.
		 * This k must satisfy the following 6 inequalities (some may be redundant):
		 * (1) k >= j-jmax[y];
		 * (2) k <= j-jmin[y]; 
		 *     1 and 2 guarantee (j-k) is within state y's j band
		 *
		 * (3) k >= hdmin[z][j-jmin[z]];
		 * (4) k <= hdmax[z][j-jmin[z]]; 
		 *     3 and 4 guarantee k is within z's j=(j), d band
		 *
		 * (5) k >= d-hdmax[y][j-jmin[y]-k];
		 * (6) k <= d-hdmin[y][j-jmin[y]-k]; 
		 *     5 and 6 guarantee (d-k) is within state y's j=(j-k) d band
		 */
		kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
		/* kmin satisfies inequalities (1) and (3) */
		kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
		/* kmax satisfies inequalities (2) and (4) */
		/* RHS of inequalities 5 and 6 are dependent on k, so we check
		 * for these within the next for loop.
		 */
		alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* initialize */
		for(k = kmin; k <= kmax; k++)
		  {
		    if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
		      {
			/* for current k, all 6 inequalities have been satisified 
			 * so we know the cells corresponding to the platonic 
			 * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
			 * alpha[z][j][k] are all within the bands. These
			 * cells correspond to alpha[v][jp_v][dp_v], 
			 * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
			 * and alpha[z][jp_z][k-hdmin[jp_z]];
			 */
			kp_z = k-hdmin[z][jp_z];
			dp_y = d-hdmin[y][jp_y-k];

			if ((sc = alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]) 
			    > alpha[v][jp_v][dp_v])
			  {
			    alpha[v][jp_v][dp_v] = sc;
			    if (ret_shadow != NULL) kshad[jp_v][dp_v] = kp_z;
			  }
		      }
		  }
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	  }
      }				/* finished calculating deck v. */
      
      /* The following loops originally access alpha[v][j0][W] but the index W will be
	 in different positions due to the bands */

      if(j0 >= jmin[v] && j0 <= jmax[v])
	{
	  jp_v = j0 - jmin[v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v])
	    {
	      Wp = W - hdmin[v][jp_v];
	      /* If we get here alpha[v][jp_v][Wp] is a valid cell
	       * in the banded alpha matrix, corresponding to 
	       * alpha[v][j0][W] in the platonic matrix.
	       */
	      /* Check for local begin getting us to the root.
	       * This is "off-shadow": if/when we trace back, we'll handle this
	       * case separately (and we'll know to do it because we'll immediately
	       * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	       * to jump right to state b; see below)
	       */
	      if (allow_begin && alpha[v][jp_v][Wp] + cm->beginsc[v] > bsc) 
		{
		  b   = v;
		  bsc = alpha[v][jp_v][Wp] + cm->beginsc[v];
		}
	    }
	}
      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (v == 0)
	{
	  if(j0 >= jmin[0] && j0 <= jmax[0])
	    {
	      jp_v = j0 - jmin[v];
	      if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v])
		{
		  if (allow_begin && v == 0 && bsc > alpha[0][jp_v][Wp]) {
		    alpha[0][jp_v][Wp] = bsc;
		    if (ret_shadow != NULL) yshad[jp_v][Wp] = USED_LOCAL_BEGIN;
		  }
		}
	    }
	}
    } /* end loop over all v */
  /* cm_fhb_mx_Dump(stdout, mx); */
  /*debug_print_alpha_banded_jd(alpha, cm, L, jmin, jmax, hdmin, hdmax);*/

  Wp = W - hdmin[vroot][j0-jmin[vroot]];
  sc =     alpha[vroot][j0-jmin[vroot]][Wp];

  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */

  if (ret_shadow != NULL) *ret_shadow = shadow;
  /*printf("inside jd me returning sc: %f\n", sc);*/

  for (v = 0; v < cm->M; v++) free(init_scAA[v]);
  free(init_scAA);
  free(el_scA);
  free(yvalidA);

  return sc;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}


/* Function: Fast_CYKInside_b_jd()
 * Incept:   EPN, Fri Oct 26 09:31:43 2007
 * 
 * Note:     based on CYKInside_b_jd() [11.04.05] which was based on CYKInside_b() 
 *           which was based on CYKInside() [SRE, Sun Jun  3 19:48:33 2001 [St. Louis]]
 *
 * Purpose:  Wrapper for the insideT_b_jd_me() routine - solve
 *           a full alignment problem, return the traceback
 *           and the score, without dividing & conquering, using bands.
 *           
 *           Analogous to CYKDivideAndConquer() in many respects;
 *           see the more extensive comments in that function for
 *           more details on shared aspects.
 *           
 * Args:     cm     - the covariance model
 *           sq     - the sequence, 1..L
 *           r      - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *           dmin   - minimum d bound for each state v; [0..v..M-1]
 *           dmax   - maximum d bound for each state v; [0..v..M-1]
 *
 * Returns:  score of the alignment in bits.
 */
float
Fast_CYKInside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, 
		    CP9Bands_t *cp9b, CM_FHB_MX *mx)
{
  /* Contract check */
  if(dsq == NULL) cm_Fail("ERROR in Fast_CYKInside_b_jd(), dsq is NULL.\n");

  Parsetree_t *tr;
  int          z;
  float        sc;

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));
    printf("alignment strategy:CYKInside_b_jd:b:nosmall\n"); 
    printf("L: %d\n", L);*/

  /* Trust, but verify.
   * Check out input parameters.
   */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CM_LOCAL_BEGIN)) cm_Fail("internal error: we're not in local mode, but r is not root");
    if (cm->stid[r] != MATP_MP && cm->stid[r] != MATL_ML &&
	cm->stid[r] != MATR_MR && cm->stid[r] != BIF_B)
      cm_Fail("internal error: trying to do a local begin at a non-mainline start");
  }

  /* Create the parse tree, and initialize.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */
  z  = cm->M-1;
  sc = 0.;

  /* Deal with case where we already know a local entry transition 0->r
   */
  if (r != 0)
    {
      InsertTraceNode(tr, 0,  TRACE_LEFT_CHILD, i0, j0, r);
      z  =  CMSubtreeFindEnd(cm, r);
      sc =  cm->beginsc[r];
    }

  /* Solve the whole thing with one call to insideT_b_jd.  This calls
     a memory efficient insideT function, which only allocates cells
     in alpha within the bands. 
   */
  sc += Fast_insideT_b_jd_me(cm, dsq, L, tr, r, z, i0, j0, (r==0), cp9b, mx);

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  /*printf("returning from CYKInside_b_jd() sc : %f\n", sc); */

  return sc;
}

/* Function: Fast_insideT_b_jd_me()
 *           EPN 03.29.06
 * Note:     based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *           only difference is memory efficient bands on the j and d dimensions
 *
 * Purpose:  Call inside, get vjd shadow matrix;
 *           then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 */
float
Fast_insideT_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
		     int r, int z, int i0, int j0, 
		     int allow_begin, CP9Bands_t *cp9b, CM_FHB_MX *mx)
{
  /* Contract check */
  if(dsq == NULL) cm_Fail("Fast_insideT_b_jd_me(), dsq is NULL.");

  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;
  int       jp_v;               /* j-jmin[v] for current j, and current v */
  int       dp_v;               /* d-hdmin[v][jp_v] for current j, current v, current d*/
  int       jp_z;               /* j-jmin[z] for current j, and current z */
  int       kp_z;               /* the k value (d dim) from the shadow matrix
				 * giving the len of right fragment offset in deck z,
				 * k = kp_z + hdmin[z][jp_z]*/
  int       *jmin = cp9b->jmin;
  int       *jmax = cp9b->jmax;
  int     **hdmin = cp9b->hdmin;
  int     **hdmax = cp9b->hdmax;

  sc = Fast_inside_b_jd_me(cm, dsq, L, r, z, i0, j0, 
			   mx, &shadow,	     /* return a shadow matrix to me. */
			   allow_begin,      /* TRUE to allow local begins */
			   &b, &bsc,	     /* if allow_begin is TRUE, gives info on optimal b */
			   cp9b);

  pda = esl_stack_ICreate();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];

  while (1) {
    if(cm->sttype[v] != EL_st && d > hdmax[v][jp_v]) cm_Fail("insideT_b_jd(). d : %d > hdmax[%d] (%d)\n", d, v, hdmax[v]);
    if(cm->sttype[v] != EL_st && d < hdmin[v][jp_v]) cm_Fail("insideT_b_jd(). d : %d < hdmin[%d] (%d)\n", d, v, hdmin[v]);
    
    if (cm->sttype[v] == B_st) {
      kp_z = ((int **) shadow[v])[jp_v][dp_v];   /* kp = offset len of right fragment */
      z = cm->cnum[v];
      jp_z = j-jmin[z];
      k = kp_z + hdmin[z][jp_z];  /* k = offset len of right fragment */
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */
      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else {
      yoffset = ((char **) shadow[v])[jp_v][dp_v];
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	  jp_v = j;
	  dp_v = d;
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  return sc;
}


/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_FASTALIGN_BENCHMARK
/* gcc -o benchmark-fastalign -g -O2 -I. -L. -I../easel -L../easel -DIMPL_FASTALIGN_BENCHMARK cm_fastalign.c -linfernal -leasel -lm
 * ./benchmark-fastalign <cmfile> <seqfile>
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
  { "-e",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "emit sequences from CM, don't randomly create them", 0 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "configure CM/HMM for local alignment", 0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of target seqs",                          0 },
  { "-L",        eslARG_INT,     NULL, NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs, default: consensus length", 0 },
  { "-o",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute original CYK HMM banded alignment implementation", 0 },
  { "--scan",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "run in scan mode, not alignment mode", 0 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use posterior sums during HMM band calculation (widens bands)", 0 },
  { "--dlev",    eslARG_INT,    "0",   NULL, "0<=n<=3",NULL,NULL,NULL, "set verbosity of debugging print statements to <n>", 0 },
  { "--check",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "check that HMM posteriors are correctly calc'ed", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "benchmark driver for fast HMM banded CYK alignment and scanning algorithm";

int 
main(int argc, char **argv)
{
  int status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  CM_t           *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  int             i;
  float           sc, bsc;
  int             b;
  char           *cmfile = esl_opt_GetArg(go, 1);
  CMFILE         *cmfp;	    /* open input CM file stream */
  int             L;        /* length of sequence */
  int             do_random;
  int             N = esl_opt_GetInteger(go, "-N");
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, either randomly created, or emitted from CM (if -e) */
  int             do_align;
  int             v;
  CM_FHB_MX      *mx;
  Parsetree_t    *slowtr, *fasttr;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  do_align = TRUE;
  if(esl_opt_GetBoolean(go, "--scan")) do_align = FALSE; 

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  /* determine sequence length */
  if(esl_opt_IsDefault(go, "-L")) L = cm->clen;      
  else                            L = esl_opt_GetInteger(go, "-L");

  /* configure CM for HMM banded alignment */
  //cm->config_opts |= CM_CONFIG_QDB;
  
  /* cm->config_opts |= CM_CONFIG_ZEROINSERTS; */
  cm->align_opts  |= CM_ALIGN_TIME;
  if(do_align) { 
    cm->align_opts  |= CM_ALIGN_NOSMALL;
    cm->align_opts  |= CM_ALIGN_HBANDED;
    if(esl_opt_GetBoolean(go, "--sums")) cm->align_opts |= CM_ALIGN_SUMS;
  }
  else /* don't align, scan */ {
    cm->search_opts  |= CM_SEARCH_HBANDED;
    cm->search_opts  |= CM_SEARCH_HMMSCANBANDS;
    if(esl_opt_GetBoolean(go, "--sums")) cm->search_opts |= CM_SEARCH_SUMS;
  }
  if(esl_opt_GetBoolean(go, "-l")) { 
    cm->config_opts  |= CM_CONFIG_LOCAL;
    cm->config_opts  |= CM_CONFIG_HMMLOCAL;
    cm->config_opts  |= CM_CONFIG_HMMEL;
  }
  if(esl_opt_GetBoolean(go, "--check")) cm->align_opts |= CM_ALIGN_CHECKFB;

  ConfigCM(cm, NULL, NULL);
  CP9Bands_t *cp9b;
  cp9b = AllocCP9Bands(cm, cm->cp9);

  /* get sequences */
  if(do_random) {
    double *dnull;
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);
    seqs_to_aln = RandomEmitSeqsToAln(r, cm->abc, dnull, 1, N, L, FALSE);
    free(dnull);
  }
  else /* don't randomly generate seqs, emit them from the CM */
    seqs_to_aln = CMEmitSeqsToAln(r, cm, 1, N, FALSE);

  /* create the matrix, it'll be empty initially */
  mx = cm_fhb_mx_Create(cm->M);

  FILE *fpfast;
  FILE *fpslow;
  for (i = 0; i < N; i++)
    {
      L = seqs_to_aln->sq[i]->n;

      esl_stopwatch_Start(w);
      cp9_Seq2Bands(cm, seqs_to_aln->sq[i]->dsq, 1, L, cp9b, 0);
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "Exptl Band calc", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      esl_stopwatch_Start(w);
      sc = Fast_CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, NULL, cp9b, mx);
      /* sc = Fast_CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, &fasttr, cp9b, mx); */
      printf("%4d %-30s %10.4f bits ", (i+1), "CYKInside_b_jd FAST (): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /* fpfast = fopen("tempfast", "w");
	 ParsetreeDump(fpfast, fasttr, cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); */

      if(esl_opt_GetBoolean(go, "-o")) {
	esl_stopwatch_Start(w);
	/* sc = CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, &slowtr, cp9b->jmin, */
	sc = CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, NULL, cp9b->jmin, 
			    cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);
	printf("%4d %-30s %10.4f bits ", (i+1), "CYKInside_b_jd SLOW (): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*fpslow = fopen("tempslow", "w");
	  ParsetreeDump(fpslow, slowtr, cm, seqs_to_aln->sq[i]->dsq, NULL, NULL);*/
      }
    }

  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  FreeSeqsToAln(seqs_to_aln);
  FreeCP9Bands(cp9b);

  return 0;

 ERROR:
  cm_Fail("memory allocation error");
}



#endif /*IMPL_FASTALIGN_BENCHMARK*/

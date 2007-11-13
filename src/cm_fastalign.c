/* cm_fastalign.c
 * Optimized DP functions for HMM banded CM alignment.
 * 
 * Functions, and their non-optimized analogs: 
 * optimized version                slow, old, reference version      ~speedup
 * ----------------------------     --------------------------------- --------
 * fast_cyk_inside_align_hb()   --> hbandcyk.c:inside_b_jd_me()           1.5 
 * fast_cyk_insideT_align_hb()  --> hbandcyk.c:insideT_b_jd_me()          N/A
 * FastCYKInsideAlignHB()       --> hbandcyk.c:CYKInside_b_jd()           N/A
 * FastFInsideAlignHB()         --> cm_postprob.c:FInside_b_jd_me()       1.4
 * FastFOutsideAlignHB()        --> cm_postprob.c:FOutside_b_jd_me()      1.4
 * FastIInsideAlignHB()         --> cm_postprob.c:IInside_b_jd_me()       1.25
 * FastIOutsideAlignHB()        --> cm_postprob.c:IOutside_b_jd_me()      1.35
 *
 * In general the float versions of inside and outside are about 1.3 -
 * 1.4X slower than the integer versions (for optimized and non-optimized).
 * 
 * Speedups are approximate, and based on tests with 2 models, an SSU
 * model and a RNaseP model. Tests were performed with
 * benchmark-fastalign, a standalone executable included at the end of
 * this file that must be separately compiled (version tested was
 * subversion revision 2204, gcc -O2 compiled
 *  EPN, Mon Nov 12 13:41:57 2007).
 * 
 * All functions use a specialized DP matrix, a CM_FHB_MX or 
 * CM_IHB_MX data structure which only allocates cells within
 * bands derived from a HMM Forward/Backward alignment of the
 * target sequence. The bands are stored in a CP9Bands_t object
 * which must be passed into each function. 
 * 
 * All but 2 functions use the float log odds scores of the
 * CM and a CM_FHB_MX matrix. The other two (FastIInsideAlignHB() 
 * and FastIOutsideAlignHB()) use integer log odds scores
 * and a CM_IHB_MX matrix. 
 * 
 * EPN, Wed Oct 10 07:20:48 2007
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


/* 
 * Function: fast_cyk_inside_align_hb()
 * based on inside_b_me() which was ...
 * based on inside()
 * Date:     EPN 03.29.06 [EPN started] 
 *           SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phase of a CYK alignment using bands 
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 *
 *           A CM_FHB_MX DP matrix must be passed in. Only
 *           cells valid within the bands given in the CP9Bands_t <cp9b>
 *           will be valid. 
 *
 *           We deal with local begins  by keeping track of the optimal
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
 *           dsq       - the digitaized sequence [i0..j0]   
 *           L         - length of the dsq
 *           vroot     - first start state of subtree (0, for whole model)
 *           vend      - last end state of subtree (cm->M-1, for whole model)
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           allow_begin-TRUE to allow 0->b local alignment begin transitions. 
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           cp9b      - HMM bands for <dsq> and <cm>
 *           mx        - the dp matrix, only cells within bands in cp9b will 
 *                       be valid. 
 *                       
 * Returns: Score of the optimal alignment.  
 */
float 
fast_cyk_inside_align_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, void ****ret_shadow,  
			 int allow_begin, int *ret_b, float *ret_bsc, CP9Bands_t *cp9b, CM_FHB_MX *mx)
{
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
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrix */
  float ***alpha = mx->dp; /* pointer to the alpha DP matrix */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      j_sdr;              /* j - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y;         /* d index for state v/y in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Wp;                 /* W oalso changes depending on state */
  float    tsc;                /* a transition score */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */

  /* Contract check */
  if(dsq == NULL) cm_Fail("fast_cyk_inside_align_hb(), dsq is NULL.\n");
  if (mx == NULL) cm_Fail("fast_cyk_inside_align_hb(), mx is NULL.\n");

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  /* grow the matrix based on the current sequence and bands */
  cm_fhb_mx_GrowTo(mx, cp9b);

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* The shadow matrix, we always allocate it, so we don't have to 
   * check if it's null in the depths of the DP recursion.
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
  ESL_ALLOC(shadow, sizeof(void **) * cm->M);
  for (v = 0; v < cm->M; v++) shadow[v] = NULL;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to IMPOSSIBLE */
  //esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* Main recursion */
  for (v = vend; v >= vroot; v--) {
    float const *esc_v = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];
    /* Get a shadow deck to fill in and initialize all valid cells for state v */
    if (cm->sttype[v] != E_st) {
      if (cm->sttype[v] == B_st) {
	kshad     = alloc_jdbanded_vjd_kshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	shadow[v] = (void **) kshad;
	/* initialize all valid cells for state v to IMPOSSIBLE (local ends are impossible for B states) */
	assert(! (NOT_IMPOSSIBLE(cm->endsc[v])));
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	    alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	    kshad[jp_v][dp_v] = USED_EL; 
	  }
	}
      } else { /* ! B_st && ! E_st */
	yshad     = alloc_jdbanded_vjd_yshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	shadow[v] = (void **) yshad;
	/* initialize all valid cells for state v */
	if(NOT_IMPOSSIBLE(cm->endsc[v])) {
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	      alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
	      yshad[jp_v][dp_v] = USED_EL; 
	    }
	  }
	}
	else { /* cm->endsc[v] == IMPOSSIBLE */
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	      alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      yshad[jp_v][dp_v] = USED_EL; 
	    }
	  }
	}
      }
    }

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr valid for state y? */
	
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
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > alpha[v][jp_v][dp_v])
		{
		  alpha[v][jp_v][dp_v] = sc; 
		  yshad[jp_v][dp_v]    = yoffset;
		}
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == IR_st) { 
      /* update alpha[v][jp_v][dp_v] cells, for IR states, loop nesting order is:
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
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > alpha[v][jp_v][dp_v])
		{
		  alpha[v][jp_v][dp_v] = sc; 
		  yshad[jp_v][dp_v]    = yoffset;
		}
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in beta[v] are independent of each other, only depending on beta[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
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
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      /* Any valid j must be within both state v and state z's j band 
       * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
       * constructed, but we'll check anyway. 
       */
      jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
      jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
      /* the main j loop */
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	kn = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
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
	   *
	   * kn and kx were set above (outside (for (dp_v...) loop) that
	   * satisfy 1-4 (b/c 1-4 are d-independent and k-independent)
	   * RHS of inequalities 5 and 6 are dependent on k, so we check
	   * for these within the next for loop.
	   */
	  for(k = kn; k <= kx; k++) { 
	    if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k]) {
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
		  > alpha[v][jp_v][dp_v]) { 
		alpha[v][jp_v][dp_v] = sc;
		kshad[jp_v][dp_v] = kp_z;
	      }
	    }
	  }
	}
      }
    }				/* finished calculating deck v. */
         
    /* The following loops originally access alpha[v][j0][W] but the index W will be
       in different positions due to the bands */
    
    if(allow_begin) { 
      if(j0 >= jmin[v] && j0 <= jmax[v]) { 
	jp_v = j0 - jmin[v];
	Wp = W - hdmin[v][jp_v];
	if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
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
	  if (alpha[v][jp_v][Wp] + cm->beginsc[v] > bsc) { 
	    b   = v;
	    bsc = alpha[v][jp_v][Wp] + cm->beginsc[v];
	  }
	}
      }
      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (v == 0) { 
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp   = W - hdmin[v][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    if (bsc > alpha[0][jp_v][Wp]) {
	      alpha[0][jp_v][Wp] = bsc;
	      if (ret_shadow != NULL) yshad[jp_v][Wp] = USED_LOCAL_BEGIN;
	    }
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp; fp = fopen("cyk.mx", "w"); cm_fhb_mx_Dump(fp, mx); fclose(fp);*/
  
  Wp = W - hdmin[vroot][j0-jmin[vroot]];
  sc =     alpha[vroot][j0-jmin[vroot]][Wp];

  if (ret_b != NULL)      *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL)    *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */
  if (ret_shadow != NULL) *ret_shadow = shadow;
  else free_vjd_shadow_matrix(shadow, cm, i0, j0);
  free(el_scA);
  free(yvalidA);

  ESL_DPRINTF1(("fast_cyk_inside_align_hb return sc: %f\n", sc));
  return sc;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

/* Function: fast_cyk_insideT_align_hb()
 * Date:     EPN 03.29.06
 * 
 * Note:     based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *           only difference is memory efficient bands on the j and d 
 *           dimensions from CP9Bands_t <cp9b> are used. 
 *
 * Purpose:  Call fast_cyk_inside_align_hb(), get vjd shadow matrix;
 *           then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 */
float
fast_cyk_insideT_align_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			  int r, int z, int i0, int j0, 
			  int allow_begin, CP9Bands_t *cp9b, CM_FHB_MX *mx)
{
  /* Contract check */
  if(dsq == NULL) cm_Fail("fast_cyk_inside_align_hbT(), dsq is NULL.");

  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;               /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			/* subseq len for bifurcs */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  int       b;                  /* local begin state */
  float     bsc;                /* local begin score */
  int       jp_v;               /* j-jmin[v] for current j, and current v */
  int       dp_v;               /* d-hdmin[v][jp_v] for current j, current v, current d*/
  int       jp_z;               /* j-jmin[z] for current j, and current z */
  int       kp_z;               /* the k value (d dim) from the shadow matrix
				 * giving the len of right fragment offset in deck z,
				 * k = kp_z + hdmin[z][jp_z]*/
  /* pointers to cp9b data for convenience */
  int       *jmin = cp9b->jmin;
  int     **hdmin = cp9b->hdmin;

  sc = fast_cyk_inside_align_hb(cm, dsq, L, r, z, i0, j0, 
				&shadow,	     /* return a shadow matrix to me. */
				allow_begin,         /* TRUE to allow local begins */
				&b, &bsc,	     /* if allow_begin is TRUE, gives info on optimal b */
				cp9b, mx);           /* the bands and the banded mx */
  
  pda = esl_stack_ICreate();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;
  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];

  while (1) {
    ESL_DASSERT1((!(cm->sttype[v] != EL_st && d > hdmax[v][jp_v])));
    ESL_DASSERT1((!(cm->sttype[v] != EL_st && d < hdmin[v][jp_v])));
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

/* Function: FastCYKInsideAlignHB()
 * Incept:   EPN, Fri Oct 26 09:31:43 2007
 * 
 * Note:     based on CYKInside_b_jd() [11.04.05] which was based on CYKInside_b() 
 *           which was based on CYKInside() [SRE, Sun Jun  3 19:48:33 2001 [St. Louis]]
 *
 * Purpose:  Wrapper for the fast_cyk_insideT_hb() routine - solve
 *           a full alignment problem, return the traceback
 *           and the score, without dividing & conquering, but by
 *           using bands on the j and d dimensions of the DP matrix.
 *           Bands derived by HMM Forward/Backward runs.
 *           
 *           Analogous to CYKDivideAndConquer() in many respects;
 *           see the more extensive comments in that function for
 *           more details (at least on those aspects that are shared).
 *           
 * Args:     cm     - the covariance model
 *           dsq    - the digitized sequence, 1..L
 *           r      - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *           cp9b      - HMM bands for <dsq> and <cm>
 *           mx        - the dp matrix, only cells within bands in cp9b will 
 *                       be valid. 
 *
 * Returns:  score of the alignment in bits.
 */
float
FastCYKInsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, 
		     CP9Bands_t *cp9b, CM_FHB_MX *mx)
{
  /* Contract check */
  if(dsq == NULL) cm_Fail("ERROR in Fast_CYKInside_b_jd(), dsq is NULL.\n");

  Parsetree_t *tr;
  int          z;
  float        sc;

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));
   *printf("alignment strategy:CYKInside_b_jd:b:nosmall\n"); 
  */

  /* Trust, but verify.
   * Check out input parameters.
   */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CMH_LOCAL_BEGIN)) cm_Fail("internal error: we're not in local mode, but r is not root");
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
  if (r != 0) {
    InsertTraceNode(tr, 0,  TRACE_LEFT_CHILD, i0, j0, r);
    z  =  CMSubtreeFindEnd(cm, r);
    sc =  cm->beginsc[r];
  }

  /* Solve the whole thing with one call to fast_cyk_insideT_align_hb.  
   */
  sc += fast_cyk_insideT_align_hb(cm, dsq, L, tr, r, z, i0, j0, (r==0), cp9b, mx);

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  ESL_DPRINTF1(("returning from FastCYKInsideAlignHB() sc : %f\n", sc)); 
  return sc;
}

/*
 * Function: FastIInsideAlignHB()
 * Date:     EPN, Wed Nov  7 15:17:47 2007 [updated]
 *
 * Purpose:  Run the inside algorithm on a target sequence using bands 
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses int log odds scores.
 * 
 *           Very similar with fast_cyk_inside_align_hb(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *           - uses ints, not floats (FastFInsideAlignHB() uses floats).
 *
 *           This function complements FastIOutsideAlignHB().
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           cp9b      - HMM bands for <dsq> and <cm>
 *                       
 * Returns:  log P(S|M)/P(S|R), as a bit score
 */
float 
FastIInsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, CM_IHB_MX *mx)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    fsc;		/* the final score, floatized */
  int      tsc;         /* a temporary variable holding a transition score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  int      bsc;		/* summed score for using all local begin states */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  int     *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrix */
  int   ***alpha = mx->dp; /* pointer to the alpha DP matrix */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      j_sdr;              /* j - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y;         /* d index for state v/y in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Wp;                 /* W oalso changes depending on state */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */

  /* Contract check */
  if(dsq == NULL) cm_Fail("fast_cyk_inside_align_hb(), dsq is NULL.\n");
  if (mx == NULL) cm_Fail("fast_cyk_inside_align_hb(), mx is NULL.\n");

  /* Allocations and initializations */
  bsc = -INFTY;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  /* grow the matrix based on the current sequence and bands */
  cm_ihb_mx_GrowTo(mx, cp9b);
  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(int) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->iel_selfsc * d;
  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands)
   */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to -INFTY */
  esl_vec_ISet(alpha[0][0], mx->ncells_valid, -INFTY);

  /* Main recursion  */
  for (v = cm->M-1; v >= 0; v--) {
    int const *esc_v = cm->ioesc[v]; 
    int const *tsc_v = cm->itsc[v];
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];
    
    /* initialize all valid cells for state v to the local end prob, if they're allowed  */
    if(cm->endsc[v] != -INFTY) {
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) 
	  alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->iendsc[v];
      }
    }

    /* E_st: easy, no children, and d must be 0 for all valid j */
    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
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
	      alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], -INFTY);
	}
      }
    }
    else if(cm->sttype[v] == IR_st) { 
      /* update alpha[v][jp_v][dp_v] cells, for IR states, loop nesting order is:
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
	      alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], -INFTY);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in beta[v] are independent of each other, only depending on beta[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
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
	  dpn     = dn - hdmin[v][jp_v];
	  dpx     = dx - hdmin[v][jp_v];
	  dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y_sdr][dp_y_sd] + tsc));;
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	  }
      default: /* no emission */
	break;
      }
      /* ensure all cells are >= -INFTY */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], -INFTY);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      /* Any valid j must be within both state v and state z's j band 
       * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
       * constructed, but we'll check anyway. 
       */
      jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
      jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
      /* the main j loop */
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	kn = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
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
	   *
	   * kn and kx were set above (outside (for (dp_v...) loop) that
	   * satisfy 1-4 (b/c 1-4 are d-independent and k-independent)
	   * RHS of inequalities 5 and 6 are dependent on k, so we check
	   * for these within the next for loop.
	   */
	  for(k = kn; k <= kx; k++) { 
	    if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k]) {
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

	      alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]); 
	    }
	  }
	}
      }
    }				/* finished calculating deck v. */
      
    if(cm->flags & CMH_LOCAL_BEGIN) { /* if local begins are on */
      if(j0 >= jmin[v] && j0 <= jmax[v]) { 
	jp_v = j0 - jmin[v];
	Wp = W - hdmin[v][jp_v];
	if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	  /* If we get here alpha[v][jp_v][Wp] is a valid cell
	   * in the banded alpha matrix, corresponding to 
	   * alpha[v][j0][W] in the platonic matrix.
	   */
	  /* Check for local begin getting us to the root.
	   */
	  bsc = ILogsum(bsc, (alpha[v][jp_v][Wp] + cm->ibeginsc[v]));
	}
      }
      /* include the bsc as part of alpha[0][jp_v][Wp] */
      if (v == 0) { 
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp = W - hdmin[0][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    alpha[0][jp_v][Wp] = ILogsum(alpha[0][jp_v][Wp], bsc);
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp; fp = fopen("ins.mx", "w"); cm_ihb_mx_Dump(fp, mx); fclose(fp);*/

  Wp  = W - hdmin[0][j0-jmin[0]];
  fsc = Scorify(alpha[0][j0-jmin[0]][Wp]);

  free(el_scA);
  free(yvalidA);

  ESL_DPRINTF1(("FastIInsideAlignHB() return sc: %f\n", fsc));
  return fsc;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

/*
 * Function: FastIOutsideAlignHB()
 * Date:     EPN, Wed Nov  7 15:24:02 2007
 *
 * Purpose:  Run the outside algorithm using bands
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses integer log odds scores.
 *
 *           A CM_IHB_MX DP matrix must be passed in. Only
 *           cells valid within the bands given in the CP9Bands_t <cp9b>
 *           will be valid. 
 *
 *           The DP recursion has been 'optimized' for all state types
 *           except IL, IR, BEGL_S, BEGR_S. The main optimization
 *           is a change in nesting order of the for loops:
 *           optimized order:     for v { for y { for j { for d {}}}}
 *           non-optimized order: for v { for j { for d { for y {}}}}
 * 
 *           ILs and IRs are not optmized because they can self transit
 *           so mx[v][j][d] must be fully calc'ed before mx[v][j][d+1] can 
 *           be calced. BEGL_S and BEGR_S are not optimized b/c 
 *           they require searching for optimal d and k, which complicates
 *           the enforcement of the bands and makes this optimization strategy
 *           impossible.
 *
 *           If <do_check> is TRUE (the CM must not be in local mode) 
 *           we check that the outside calculations are consistent 
 *           with the inside calculations (in ins_mx). 
 *           This check is described in comments towards the end of 
 *           the function. 
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ins_mx    - the dp matrix from the Inside run calculation (required)
 *           do_check  - TRUE to attempt to check 
 *                       
 * Returns:  log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *           ends are on (see *** comment towards end of function).
 */
float 
FastIOutsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b,
		    CM_IHB_MX *mx, CM_IHB_MX *ins_mx, int do_check)
{
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  int      isc;     	       /* a temporary variable holding an int score */
  float    fsc;     	       /* a temporary variable holding a float score */
  int    **esc_vAA;            /* ptr to cm->ioesc, optimized emission scores */
  int      iescore;	       /* an emission score, tmp variable */
  int      W;		       /* subsequence length */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      jp;		       /* j': relative position in the subsequence  */
  int      bsc;		       /* total score for using local begin states */
  int      ireturn_sc;         /* P(S|M)/P(S|R), a scaled int*/
  float    freturn_sc;         /* P(S|M)/P(S|R), a float (Scorified ireturn_sc) */
  /* band related variables */
  int      dp_v;               /* d index for state v in alpha w/mem eff bands */
  int      dp_y;               /* d index for state y in alpha w/mem eff bands */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;                 /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      kmin, kmax;         /* temporary minimum/maximum allowed k */
  int      fail_flag = FALSE;  /* set to TRUE if do_check and we see a problem */
  /* variables used only if do_check */
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */
  /* indices used in the depths of the DP recursion */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      dn, dx;             /* current minimum/maximum d allowed */
  /* DP matrix variables */
  int   ***beta  = mx->dp;     /* pointer to the Oustide DP mx */
  int   ***alpha = ins_mx->dp; /* pointer to the Inside DP mx (already calc'ed and passed in) */
  int    **beta_el = NULL;     /* the cm->M EL deck, it's special b/c it has no bands */
  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  /* Contract check */
  if (dsq == NULL)                                     cm_Fail("FastIOutsideAlignHB(), dsq is NULL.\n");
  if (mx == NULL)                                      cm_Fail("FastIOutsideAlignHB(), mx is NULL.\n");
  if (ins_mx == NULL)                                  cm_Fail("FastIOutsideAlignHB(), ins_mx is NULL.\n");
  if ((cm->flags & CMH_LOCAL_END) && do_check == TRUE) cm_Fail("FastIOutsideAlignHB(), do_check = TRUE, but CMH_LOCAL_END flag raised, we can't check with local ends.n");

  /* Allocations and initializations
   */
  bsc = -INFTY;                 /* the summed prob of all local begins */
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  esc_vAA = cm->ioesc;          /* a ptr to the optimized emission scores */
  cm_ihb_mx_GrowTo(mx, cp9b);   /* grow the matrix based on the current sequence and bands */

  /* initialize all cells of the matrix to -INFTY */
  esl_vec_ISet(beta[0][0], mx->ncells_valid, -INFTY);
  /* now set beta[0][j0][W] to 0., all parses must end there */
  jp_v = j0 - jmin[0];
  Wp = W - hdmin[0][jp_v];
  assert(W >= hdmin[0][jp_v]);
  assert(W <= hdmax[0][jp_v]);
  beta[0][jp_v][Wp] = 0.;

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   * EL deck has no bands as currently implemented. Set all cells to -INFTY;
   */
  beta_el = NULL;
  if (cm->flags & CMH_LOCAL_END) {
    beta_el = Ialloc_vjd_deck(W, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++) beta_el[j][d] = -INFTY;
    }
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }
  /* assign the ptr for EL in our DP matrix to the EL deck, it should be NULL */
  assert(beta[cm->M] == NULL);
  beta[cm->M] = beta_el;
  /* If we can do a local begin into v, overwrite -INFY with the local begin score. 
   * By definition, beta[0][j0][W] == 0.
   */ 
  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (v = 1; v < cm->M; v++) {
      if(cm->ibeginsc[v] != -INFTY) { 
	if((j0 >= jmin[v]) && (j0 <= jmax[v])) {
	  jp_v = j0 - jmin[v];
	  if((W >= hdmin[v][jp_v]) && W <= hdmax[v][jp_v]) {
	    Wp = W - hdmin[v][jp_v];
	    beta[v][jp_v][Wp] = cm->ibeginsc[v];
	  }
	}
      }
    }
  }
  /* done allocation/initialization */

  /* Recursion: main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      y = cm->plast[v];	/* the parent bifurcation    */
      z = cm->cnum[y];	/* the other (right) S state */
      for (j = jmax[v]; j >= jmin[v]; j--) {
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;
	for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--) {
	  dp_v = d - hdmin[v][jp_v];
	  /* Find the first k value that implies a valid cell in the y and z decks.
	   * This k must satisfy the following 8 inequalities (some may be redundant):
	   * NOTE: these are different from those in Inside() (for one thing, v and y
	   *       (BEGL_S and BIF_B here respectively) are switched relative to Inside.
	   *
	   * (1) k <= jmax[y] - j;
	   * (2) k >= jmin[y] - j;
	   * (3) k <= jmax[z] - j;
	   * (4) k >= jmin[z] - j;
	   *     1 and 2 guarantee (j+k) is within state y's j band
	   *     3 and 4 guarantee (j+k) is within state z's j band
	   *
	   * (5) k >= hdmin[y][j-jmin[y]+k] - d;
	   * (6) k <= hdmax[y][j-jmin[y]+k] - d; 
	   *     5 and 6 guarantee k+d is within y's j=(j+k), d band
	   *
	   * (7) k >= hdmin[z][j-jmin[z]+k];
	   * (8) k <= hdmax[z][j-jmin[z]+k]; 
	   *     5 and 6 guarantee k is within state z's j=(j+k) d band
	   */
	  kmin = ESL_MAX(jmin[y], jmin[z]) - j;
	  kmax = ESL_MIN(jmax[y], jmax[z]) - j;
	  /* kmin and kmax satisfy inequalities (1-4) */
	  /* RHS of inequalities 5-8 are dependent on k, so we check
	   * for these within the next for loop. */
	  for(k = kmin; k <= kmax; k++) {
	    if(k < (hdmin[y][jp_y+k] - d) || k > (hdmax[y][jp_y+k] - d)) continue; 
	    /* above line continues if inequality 5 or 6 is violated */
	    if(k < (hdmin[z][jp_z+k])     || k > (hdmax[z][jp_z+k]))     continue; 
	    /* above line continues if inequality 7 or 8 is violated */
		  
	    /* if we get here for current k, all 8 inequalities have been satisified 
	     * so we know the cells corresponding to the platonic 
	     * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
	     * alpha[z][j+k][k] are all within the bands. These
	     * cells correspond to beta[v][jp_v][dp_v], 
	     * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
	     * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
	     */
	    kp_z = k-hdmin[z][jp_z+k];
	    dp_y = d-hdmin[y][jp_y+k];
	    beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] 
								+ alpha[z][jp_z+k][kp_z]));
	  }
	}
      }
    } /* end of 'if (cm->stid[v] == BEGL_S */
    else if (cm->stid[v] == BEGR_S) {
      y = cm->plast[v];	  /* the parent bifurcation    */
      z = cm->cfirst[y];  /* the other (left) S state  */
      jn = ESL_MAX(jmin[v], jmin[y]);
      jx = ESL_MIN(jmax[v], jmax[y]);
      for (j = jx; j >= jn; j--) {
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;

	dn = ESL_MAX(hdmin[v][jp_v], j-jmax[z]);
	dx = ESL_MIN(hdmax[v][jp_v], jp_z);
	/* above makes sure that j,d are valid for state z: (jmin[z] + d) >= j >= (jmax[z] + d) */
	for (d = dx; d >= dn; d--) {
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  /* Find the first k value that implies a valid cell in the y and z decks.
	   * This k must satisfy the following 4 inequalities (some may be redundant):
	   * NOTE: these are different from those in Inside() (for one thing, v and y
	   *       (BEGR_S and BIF_B here respectively) are switched relative to Inside.
	   *
	   * (1) k >= hdmin[y][j-jmin[y]] - d;
	   * (2) k <= hdmax[y][j-jmin[y]] - d;
	   *     1 and 2 guarantee (d+k) is within state y's j=(j) d band
	   *
	   * (3) k >= hdmin[z][j-jmin[z]-d];
	   * (4) k <= hdmax[z][j-jmin[z]-d];
	   *     3 and 4 guarantee k is within z's j=(j-d) d band
	   *
	   */
	  kmin = ESL_MAX((hdmin[y][jp_y]-d), (hdmin[z][jp_z-d]));
	  kmax = ESL_MIN((hdmax[y][jp_y]-d), (hdmax[z][jp_z-d]));
	  /* kmin and kmax satisfy inequalities (1-4) */
	  for(k = kmin; k <= kmax; k++) { 
	    /* for current k, all 4 inequalities have been satisified 
	     * so we know the cells corresponding to the platonic 
	     * matrix cells beta[v][j][d], beta[y][j][d+k], and
	     * alpha[z][j-d][k] are all within the bands. These
	     * cells correspond to beta[v][jp_v][dp_v], 
	     * beta[y][jp_y+k][d-hdmin[y][jp_y]+k],
	     * and alpha[z][jp_z-d][k-hdmin[z][jp_z-d]];
	     */
	    kp_z = k-hdmin[z][jp_z-d];
	    dp_y = d-hdmin[y][jp_y];
	    beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
								+ alpha[z][jp_z-d][kp_z]));
	  }
	}
      }
    } /* end of 'else if (cm->stid[v] == BEGR_S */
    else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
      /* ILs and IRs can self transit, this means that beta[v][j][d] must be fully calculated
       * before beta[v][j][d+1] can be started to be calculated, forcing the following nesting order:
       * for j { for d { for y { } } } 
       * for non-self-transitioners, we can do a more efficient nesting order (see below)  
       */
      for (j = jmax[v]; j >= jmin[v]; j--) {
	jp_v = j - jmin[v];
	for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--) {
	  i = j-d+1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    
	    /* Note: this looks like it can be optimized, I tried but my 'optimization' slowed the code, so I reverted [EPN] */
	    switch(cm->sttype[y]) {
	    case MP_st: 
	      if (j == j0 || d == j) continue; /* boundary condition */
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+2) < hdmin[y][(jp_y+1)] || (d+2) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
		   */
	      dp_y = d - hdmin[y][jp_y+1];  /* d index for state y */
	      iescore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
								  + cm->itsc[y][voffset] + iescore));
	      break;
	      
	    case ML_st:
	    case IL_st: 
	      if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
	      if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hdmin[y][jp_y] || (d+1) > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j][d+1] in the platonic matrix.
	       */
	      dp_y = d - hdmin[y][jp_y];  /* d index for state y */
	      iescore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
								  + cm->itsc[y][voffset] + iescore));
	      break;
	      
	    case MR_st:
	    case IR_st:
	      if (j == j0) continue;
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hdmin[y][(jp_y+1)] || (d+1) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
	       */
	      dp_y = d - hdmin[y][(jp_y+1)];  /* d index for state y */
	      iescore = esc_vAA[y][dsq[j+1]];
	      /*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
								  + cm->itsc[y][voffset] + iescore));
	      break;
	      
	    case S_st:
	    case E_st:
	    case D_st:
	      if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if (d < hdmin[y][jp_y] || d > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
	       * corresponding to alpha[y][j][d] in the platonic matrix.
	       */
	      dp_y = d - hdmin[y][jp_y];  /* d index for state y */
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->itsc[y][voffset])); 
	      break;
	    } /* end of switch(cm->sttype[y] */  
	  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	  if (beta[v][jp_v][dp_v] < -INFTY) beta[v][jp_v][dp_v] = -INFTY;
	} /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
      } /* end loop over jp. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] == IL_st || cm->sttype[v] == IR_st' */
    else { /* state v is not BEGL_S, BEGL_R IL nor IR (must be ML, MP, MR, D, S or E */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in beta[v] are independent of each other, only depending on beta[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	sdr = StateRightDelta(cm->sttype[y]);
	sd  = StateDelta(cm->sttype[y]);
	emitmode = Emitmode(cm->sttype[y]);
	/* determine min j (jn) and max j (jx) that are valid for v and y */
	jn = ESL_MAX(jmin[v], ESL_MIN(jmin[y] - sdr, jmax[y]));
	jx = ESL_MIN(jmax[v], ESL_MAX(jmax[y] - sdr, jmin[y]));
	jx = ESL_MIN(jx, jmax[y] - sdr);
	for (j = jx; j >= jn; j--) {
	  jp_v = j - jmin[v];
	  jp_y = j - jmin[y];
	  ESL_DASSERT1((j+sdr >= jmin[y] && j+sdr <= jmax[y]));
	  
	  /* determine min d (dn) and max d (dx) that are valid for v and y and j */
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y + sdr] - sd);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y + sdr] - sd);
	  dp_v = dx - hdmin[v][jp_v];
	  dp_y = dx - hdmin[y][jp_y + sdr];
	  i    = j-dx+1;
	  
	  /* for each emit mode, update beta[v][jp_v][dp_v] for all valid d = dp_v */
	  switch(emitmode) { 
	  case EMITPAIR:  /* MP_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      iescore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->itsc[y][voffset] + iescore));
	    }
	    break;
	  case EMITLEFT:  /* ML_st, IL_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      iescore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->itsc[y][voffset] + iescore));
	    }
	    break;
	  case EMITRIGHT:  /* MR_st, IR_st */
	    iescore = esc_vAA[y][dsq[j+1]]; /* not dependent on i */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->itsc[y][voffset] + iescore));
	    }
	    break;
	  case EMITNONE:  /* D_st, S_st, E_st*/
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->itsc[y][voffset]));
	    }
	    break;
	  } /* end of switch(emitmode) */
	} /* end of for j = jx; j >= jn; j-- */
      } /* end of for y = plast[v]... */
    } /* ends else entered for non-BEGL_S/BEGR_S/IL/IR states*/	
    /* we're done calculating deck v for everything but local begins */

    /* deal with local alignment end transitions v->EL (EL = deck at M.) */
    if (cm->iendsc[v] != -INFTY) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);   /* note emitmode is for state v */

      /* determine min j (jn) and max j (jx) that are valid for v and y */
      jn = ESL_MIN(jmin[v] - sdr, jmax[v]);
      jx = ESL_MAX(jmax[v] - sdr, jmin[v]);
      for (j = jx; j >= jn; j--) {
	jp_v = j - jmin[v];
	ESL_DASSERT1((j+sdr >= jmin[y] && j+sdr <= jmax[v]));
	  
	/* determine min d (dn) and max d (dx) that are valid for v and j */
	dn   = hdmin[v][jp_v + sdr] - sd;
	dx   = hdmax[v][jp_v + sdr] - sd;
	dp_v = dx - hdmin[v][jp_v];
	i    = j-dx+1;

	  /* for each emit mode, update beta[v][jp_v][dp_v] for all valid d = dp_v */
	switch(emitmode) { 
	case EMITPAIR:  /* MP_st */
	  for (d = dx; d >= dn; d--, dp_v--, i++) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    iescore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->iendsc[v] 
							    + iescore));
	  }
	  break;

	case EMITLEFT:  /* ML_st, IL_st */
	  for (d = dx; d >= dn; d--, dp_v--, i++) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    iescore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v + sdr][dp_y + sd] + cm->iendsc[v] 
							    + iescore));
	  }
	  break;

	case EMITRIGHT:  /* MR_st, IR_st */
	  iescore = esc_vAA[v][dsq[j+1]]; /* not dependent on i */
	  for (d = dx; d >= dn; d--, dp_v--) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v + sdr][dp_y + sd] + cm->iendsc[v] 
							    + iescore));
	  }
	  break;

	case EMITNONE:  /* D_st, S_st, E_st*/
	  for (d = dx; d >= dn; d--, dp_v--) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v + sdr][dp_y + sd] + cm->iendsc[v]));
	    break;
	  }
	} /* end of switch over emitmodes */
      } /* end of loop over j */
    } /* end conditional section for dealing w/ v->EL local end transitions */
  } /* end loop over decks v. */

  /*FILE *fp; fp = fopen("out.mx", "w"); cm_ihb_mx_Dump(fp, mx); fclose(fp);*/

  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, EL->EL transitions. (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->iel_selfsc));
    }
  }

  Wp = W - hdmin[0][j0-jmin[0]];
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    
    /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
     *                                                = P(S|M) / P(S|R)  
     * in v,j,d coordinates this is:
     * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                = P(S|M) / P(S|R)
     */
    
    for(n = 0; n < cm->nodes; n++) {
      isc = -INFTY;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	    /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
	      printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	    isc = ILogsum(isc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
	  }
	}
      }
      fsc = Scorify(isc);
      /*printf("checking node: %d | sc: %.6f\n", n, fsc);*/
      diff = fsc - (Scorify(alpha[0][j0-jmin[0]][Wp]));
      if(diff > 0.1 || diff < -0.1) { /* we need to allow .1 bit difference slack b/c of precision
				       * with int log odds scores, the float version of this function is more stringent */
	fail_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, fsc, Scorify(alpha[0][(j0-jmin[0])][Wp]), diff);
      }
    }
  }

  /* If not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to W (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to W (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * *** If local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END)) { 
    ireturn_sc = -INFTY;
    v = cm->M-1;
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      assert(hdmin[v][jp_v] == 0);
      /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
      ireturn_sc = ILogsum(ireturn_sc, (beta[v][jp_v][0]));
    }
    freturn_sc = Scorify(ireturn_sc);
  }
  else { /* return_sc = P(S|M) / P(S|R) from Inside() */
    freturn_sc = Scorify(alpha[0][(j0-jmin[0])][Wp]);
  }

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tFast_IOutside_b_jd_me() sc : %f\n", freturn_sc));
  else
    ESL_DPRINTF1(("\tFast_IOutside_b_jd_me() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  return freturn_sc;
}

/*
 * Function: FastFInsideAlignHB()
 * Date:     EPN, Thu Nov  8 18:24:41 2007
 *
 * Purpose:  Run the inside algorithm on a target sequence using bands 
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 * 
 *           Very similar with fast_cyk_inside_align_hb(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *
 *           This function complements FastFOutsideAlignHB().
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           cp9b      - HMM bands for <dsq> and <cm>
 *                       
 * Returns:  log P(S|M)/P(S|R), as a bit score
 */
float 
FastFInsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, CM_FHB_MX *mx)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    fsc;		/* the final score */
  float    tsc;         /* a temporary variable holding a transition score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  float    bsc;		/* summed score for using all local begins */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrix */
  float ***alpha = mx->dp;     /* pointer to the alpha DP matrix */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      j_sdr;              /* j - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y;         /* d index for state v/y in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Wp;                 /* W oalso changes depending on state */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */

  /* Contract check */
  if(dsq == NULL) cm_Fail("fast_cyk_inside_align_hb(), dsq is NULL.\n");
  if (mx == NULL) cm_Fail("fast_cyk_inside_align_hb(), mx is NULL.\n");

  /* Allocations and initializations */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  /* grow the matrix based on the current sequence and bands */
  cm_fhb_mx_GrowTo(mx, cp9b);

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(int) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;
  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands)
   */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* Main recursion  */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; 
    float const *tsc_v = cm->tsc[v];
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];
    
    /* initialize all valid cells for state v to the local end prob, if they're allowed  */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) 
	  alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
      }
    }

    /* E_st: easy, no children, and d must be 0 for all valid j */
    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
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
	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == IR_st) { 
      /* update alpha[v][jp_v][dp_v] cells, for IR states, loop nesting order is:
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
	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in beta[v] are independent of each other, only depending on beta[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
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
	  dpn     = dn - hdmin[v][jp_v];
	  dpx     = dx - hdmin[v][jp_v];
	  dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y_sdr][dp_y_sd] + tsc));;
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	  }
      default: /* no emission */
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      /* Any valid j must be within both state v and state z's j band 
       * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
       * constructed, but we'll check anyway. 
       */
      jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
      jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
      /* the main j loop */
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	kn = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
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
	   *
	   * kn and kx were set above (outside (for (dp_v...) loop) that
	   * satisfy 1-4 (b/c 1-4 are d-independent and k-independent)
	   * RHS of inequalities 5 and 6 are dependent on k, so we check
	   * for these within the next for loop.
	   */
	  for(k = kn; k <= kx; k++) { 
	    if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k]) {
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

	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]); 
	    }
	  }
	}
      }
    }				/* finished calculating deck v. */
      
    if(cm->flags & CMH_LOCAL_BEGIN) { /* if local begins are on */
      if(j0 >= jmin[v] && j0 <= jmax[v]) { 
	jp_v = j0 - jmin[v];
	Wp = W - hdmin[v][jp_v];
	if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	  /* If we get here alpha[v][jp_v][Wp] is a valid cell
	   * in the banded alpha matrix, corresponding to 
	   * alpha[v][j0][W] in the platonic matrix.
	   */
	  /* Check for local begin getting us to the root.
	   */
	  bsc = FLogsum(bsc, (alpha[v][jp_v][Wp] + cm->beginsc[v]));
	}
      }
      /* include the bsc as part of alpha[0][jp_v][Wp] */
      if (v == 0) { 
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp = W - hdmin[0][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    alpha[0][jp_v][Wp] = FLogsum(alpha[0][jp_v][Wp], bsc);
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp; fp = fopen("ins.mx", "w"); cm_ihb_mx_Dump(fp, mx); fclose(fp);*/

  Wp  = W - hdmin[0][j0-jmin[0]];
  fsc = alpha[0][j0-jmin[0]][Wp];

  free(el_scA);
  free(yvalidA);

  ESL_DPRINTF1(("FastFInsideAlignHB() return sc: %f\n", fsc));
  return fsc;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

/*
 * Function: FastFOutsideAlignHB()
 * Date:     EPN, Thu Nov  8 18:40:05 2007
 *
 * Purpose:  Run the outside algorithm using bands
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 *
 *           Identical to FastIOutsideAlignHB() with except that
 *           float scores are used instead of ints. See the
 *           'Purpose' section of that function for more details.
 *           A CM_FHB_MX DP matrix must be passed in. Only
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ins_mx    - the dp matrix from the Inside run calculation (required)
 *           do_check  - TRUE to attempt to check 
 *                       
 * Returns:  log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *           ends are on (see *** comment towards end of function).
 */
float 
FastFOutsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b,
		    CM_FHB_MX *mx, CM_FHB_MX *ins_mx, int do_check)
{
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float    fsc;     	       /* a temporary variable holding a float score */
  float  **esc_vAA;            /* ptr to cm->oesc, optimized emission scores */
  float    escore;	       /* an emission score, tmp variable */
  int      W;		       /* subsequence length */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      jp;		       /* j': relative position in the subsequence  */
  float    bsc;		       /* total score for using local begin states */
  float    freturn_sc;         /* P(S|M)/P(S|R), a float (Scorified ireturn_sc) */
  /* band related variables */
  int      dp_v;               /* d index for state v in alpha w/mem eff bands */
  int      dp_y;               /* d index for state y in alpha w/mem eff bands */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;                 /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      kmin, kmax;         /* temporary minimum/maximum allowed k */
  int      fail_flag = FALSE;  /* set to TRUE if do_check and we see a problem */
  /* variables used only if do_check */
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */
  /* indices used in the depths of the DP recursion */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      dn, dx;             /* current minimum/maximum d allowed */
  /* DP matrix variables */
  float ***beta  = mx->dp;     /* pointer to the Oustide DP mx */
  float ***alpha = ins_mx->dp; /* pointer to the Inside DP mx (already calc'ed and passed in) */
  float  **beta_el = NULL;     /* the cm->M EL deck, it's special b/c it has no bands */
  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  /* Contract check */
  if (dsq == NULL)                                     cm_Fail("FastIOutsideAlignHB(), dsq is NULL.\n");
  if (mx == NULL)                                      cm_Fail("FastIOutsideAlignHB(), mx is NULL.\n");
  if (ins_mx == NULL)                                  cm_Fail("FastIOutsideAlignHB(), ins_mx is NULL.\n");
  if ((cm->flags & CMH_LOCAL_END) && do_check == TRUE) cm_Fail("FastIOutsideAlignHB(), do_check = TRUE, but CMH_LOCAL_END flag raised, we can't check with local ends.n");

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;                 /* the summed prob of all local begins */
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  esc_vAA = cm->oesc;          /* a ptr to the optimized emission scores */
  cm_fhb_mx_GrowTo(mx, cp9b);   /* grow the matrix based on the current sequence and bands */

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(beta[0][0], mx->ncells_valid, IMPOSSIBLE);
  /* now set beta[0][j0][W] to 0., all parses must end there */
  jp_v = j0 - jmin[0];
  Wp = W - hdmin[0][jp_v];
  assert(W >= hdmin[0][jp_v]);
  assert(W <= hdmax[0][jp_v]);
  beta[0][jp_v][Wp] = 0.;

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   * EL deck has no bands as currently implemented. Set all cells to IMPOSSIBLE;
   */
  beta_el = NULL;
  if (cm->flags & CMH_LOCAL_END) {
    beta_el = alloc_vjd_deck(W, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++) beta_el[j][d] = IMPOSSIBLE;
    }
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }
  /* assign the ptr for EL in our DP matrix to the EL deck, it should be NULL */
  assert(beta[cm->M] == NULL);
  beta[cm->M] = beta_el;
  /* If we can do a local begin into v, overwrite -INFY with the local begin score. 
   * By definition, beta[0][j0][W] == 0.
   */ 
  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (v = 1; v < cm->M; v++) {
      if(NOT_IMPOSSIBLE(cm->beginsc[v])) {
	if((j0 >= jmin[v]) && (j0 <= jmax[v])) {
	  jp_v = j0 - jmin[v];
	  if((W >= hdmin[v][jp_v]) && W <= hdmax[v][jp_v]) {
	    Wp = W - hdmin[v][jp_v];
	    beta[v][jp_v][Wp] = cm->beginsc[v];
	  }
	}
      }
    }
  }
  /* done allocation/initialization */

  /* Recursion: main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      y = cm->plast[v];	/* the parent bifurcation    */
      z = cm->cnum[y];	/* the other (right) S state */
      for (j = jmax[v]; j >= jmin[v]; j--) {
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;
	for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--) {
	  dp_v = d - hdmin[v][jp_v];
	  /* Find the first k value that implies a valid cell in the y and z decks.
	   * This k must satisfy the following 8 inequalities (some may be redundant):
	   * NOTE: these are different from those in Inside() (for one thing, v and y
	   *       (BEGL_S and BIF_B here respectively) are switched relative to Inside.
	   *
	   * (1) k <= jmax[y] - j;
	   * (2) k >= jmin[y] - j;
	   * (3) k <= jmax[z] - j;
	   * (4) k >= jmin[z] - j;
	   *     1 and 2 guarantee (j+k) is within state y's j band
	   *     3 and 4 guarantee (j+k) is within state z's j band
	   *
	   * (5) k >= hdmin[y][j-jmin[y]+k] - d;
	   * (6) k <= hdmax[y][j-jmin[y]+k] - d; 
	   *     5 and 6 guarantee k+d is within y's j=(j+k), d band
	   *
	   * (7) k >= hdmin[z][j-jmin[z]+k];
	   * (8) k <= hdmax[z][j-jmin[z]+k]; 
	   *     5 and 6 guarantee k is within state z's j=(j+k) d band
	   */
	  kmin = ESL_MAX(jmin[y], jmin[z]) - j;
	  kmax = ESL_MIN(jmax[y], jmax[z]) - j;
	  /* kmin and kmax satisfy inequalities (1-4) */
	  /* RHS of inequalities 5-8 are dependent on k, so we check
	   * for these within the next for loop. */
	  for(k = kmin; k <= kmax; k++) {
	    if(k < (hdmin[y][jp_y+k] - d) || k > (hdmax[y][jp_y+k] - d)) continue; 
	    /* above line continues if inequality 5 or 6 is violated */
	    if(k < (hdmin[z][jp_z+k])     || k > (hdmax[z][jp_z+k]))     continue; 
	    /* above line continues if inequality 7 or 8 is violated */
		  
	    /* if we get here for current k, all 8 inequalities have been satisified 
	     * so we know the cells corresponding to the platonic 
	     * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
	     * alpha[z][j+k][k] are all within the bands. These
	     * cells correspond to beta[v][jp_v][dp_v], 
	     * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
	     * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
	     */
	    kp_z = k-hdmin[z][jp_z+k];
	    dp_y = d-hdmin[y][jp_y+k];
	    beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] 
								+ alpha[z][jp_z+k][kp_z]));
	  }
	}
      }
    } /* end of 'if (cm->stid[v] == BEGL_S */
    else if (cm->stid[v] == BEGR_S) {
      y = cm->plast[v];	  /* the parent bifurcation    */
      z = cm->cfirst[y];  /* the other (left) S state  */
      jn = ESL_MAX(jmin[v], jmin[y]);
      jx = ESL_MIN(jmax[v], jmax[y]);
      for (j = jx; j >= jn; j--) {
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;

	dn = ESL_MAX(hdmin[v][jp_v], j-jmax[z]);
	dx = ESL_MIN(hdmax[v][jp_v], jp_z);
	/* above makes sure that j,d are valid for state z: (jmin[z] + d) >= j >= (jmax[z] + d) */
	for (d = dx; d >= dn; d--) {
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  /* Find the first k value that implies a valid cell in the y and z decks.
	   * This k must satisfy the following 4 inequalities (some may be redundant):
	   * NOTE: these are different from those in Inside() (for one thing, v and y
	   *       (BEGR_S and BIF_B here respectively) are switched relative to Inside.
	   *
	   * (1) k >= hdmin[y][j-jmin[y]] - d;
	   * (2) k <= hdmax[y][j-jmin[y]] - d;
	   *     1 and 2 guarantee (d+k) is within state y's j=(j) d band
	   *
	   * (3) k >= hdmin[z][j-jmin[z]-d];
	   * (4) k <= hdmax[z][j-jmin[z]-d];
	   *     3 and 4 guarantee k is within z's j=(j-d) d band
	   *
	   */
	  kmin = ESL_MAX((hdmin[y][jp_y]-d), (hdmin[z][jp_z-d]));
	  kmax = ESL_MIN((hdmax[y][jp_y]-d), (hdmax[z][jp_z-d]));
	  /* kmin and kmax satisfy inequalities (1-4) */
	  for(k = kmin; k <= kmax; k++) { 
	    /* for current k, all 4 inequalities have been satisified 
	     * so we know the cells corresponding to the platonic 
	     * matrix cells beta[v][j][d], beta[y][j][d+k], and
	     * alpha[z][j-d][k] are all within the bands. These
	     * cells correspond to beta[v][jp_v][dp_v], 
	     * beta[y][jp_y+k][d-hdmin[y][jp_y]+k],
	     * and alpha[z][jp_z-d][k-hdmin[z][jp_z-d]];
	     */
	    kp_z = k-hdmin[z][jp_z-d];
	    dp_y = d-hdmin[y][jp_y];
	    beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
								+ alpha[z][jp_z-d][kp_z]));
	  }
	}
      }
    } /* end of 'else if (cm->stid[v] == BEGR_S */
    else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
      /* ILs and IRs can self transit, this means that beta[v][j][d] must be fully calculated
       * before beta[v][j][d+1] can be started to be calculated, forcing the following nesting order:
       * for j { for d { for y { } } } 
       * for non-self-transitioners, we can do a more efficient nesting order (see below)  
       */
      for (j = jmax[v]; j >= jmin[v]; j--) {
	jp_v = j - jmin[v];
	for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--) {
	  i = j-d+1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    
	    /* Note: this looks like it can be optimized, I tried but my 'optimization' slowed the code, so I reverted [EPN] */
	    switch(cm->sttype[y]) {
	    case MP_st: 
	      if (j == j0 || d == j) continue; /* boundary condition */
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+2) < hdmin[y][(jp_y+1)] || (d+2) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
		   */
	      dp_y = d - hdmin[y][jp_y+1];  /* d index for state y */
	      escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case ML_st:
	    case IL_st: 
	      if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
	      if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hdmin[y][jp_y] || (d+1) > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j][d+1] in the platonic matrix.
	       */
	      dp_y = d - hdmin[y][jp_y];  /* d index for state y */
	      escore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case MR_st:
	    case IR_st:
	      if (j == j0) continue;
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hdmin[y][(jp_y+1)] || (d+1) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
	       */
	      dp_y = d - hdmin[y][(jp_y+1)];  /* d index for state y */
	      escore = esc_vAA[y][dsq[j+1]];
	      /*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case S_st:
	    case E_st:
	    case D_st:
	      if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if (d < hdmin[y][jp_y] || d > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
	       * corresponding to alpha[y][j][d] in the platonic matrix.
	       */
	      dp_y = d - hdmin[y][jp_y];  /* d index for state y */
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->tsc[y][voffset])); 
	      break;
	    } /* end of switch(cm->sttype[y] */  
	  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	  if (beta[v][jp_v][dp_v] < IMPOSSIBLE) beta[v][jp_v][dp_v] = IMPOSSIBLE;
	} /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
      } /* end loop over jp. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] == IL_st || cm->sttype[v] == IR_st' */
    else { /* state v is not BEGL_S, BEGL_R IL nor IR (must be ML, MP, MR, D, S or E */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in beta[v] are independent of each other, only depending on beta[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	sdr = StateRightDelta(cm->sttype[y]);
	sd  = StateDelta(cm->sttype[y]);
	emitmode = Emitmode(cm->sttype[y]);
	/* determine min j (jn) and max j (jx) that are valid for v and y */
	jn = ESL_MAX(jmin[v], ESL_MIN(jmin[y] - sdr, jmax[y]));
	jx = ESL_MIN(jmax[v], ESL_MAX(jmax[y] - sdr, jmin[y]));
	jx = ESL_MIN(jx, jmax[y] - sdr);
	for (j = jx; j >= jn; j--) {
	  jp_v = j - jmin[v];
	  jp_y = j - jmin[y];
	  ESL_DASSERT1((j+sdr >= jmin[y] && j+sdr <= jmax[y]));
	  
	  /* determine min d (dn) and max d (dx) that are valid for v and y and j */
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y + sdr] - sd);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y + sdr] - sd);
	  dp_v = dx - hdmin[v][jp_v];
	  dp_y = dx - hdmin[y][jp_y + sdr];
	  i    = j-dx+1;
	  
	  /* for each emit mode, update beta[v][jp_v][dp_v] for all valid d = dp_v */
	  switch(emitmode) { 
	  case EMITPAIR:  /* MP_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITLEFT:  /* ML_st, IL_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      escore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITRIGHT:  /* MR_st, IR_st */
	    escore = esc_vAA[y][dsq[j+1]]; /* not dependent on i */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITNONE:  /* D_st, S_st, E_st*/
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hdmin[v][jp_v]        &&   d       <= hdmax[v][jp_v]));
	      ESL_DASSERT1((((d + sd) >= hdmin[y][jp_y + sdr]) && ((d + sd) <= hdmax[y][jp_y + sdr])));
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset]));
	    }
	    break;
	  } /* end of switch(emitmode) */
	} /* end of for j = jx; j >= jn; j-- */
      } /* end of for y = plast[v]... */
    } /* ends else entered for non-BEGL_S/BEGR_S/IL/IR states*/	
    /* we're done calculating deck v for everything but local begins */

    /* deal with local alignment end transitions v->EL (EL = deck at M.) */
    if (NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);   /* note emitmode is for state v */

      /* determine min j (jn) and max j (jx) that are valid for v and y */
      jn = ESL_MIN(jmin[v] - sdr, jmax[v]);
      jx = ESL_MAX(jmax[v] - sdr, jmin[v]);
      for (j = jx; j >= jn; j--) {
	jp_v = j - jmin[v];
	ESL_DASSERT1((j+sdr >= jmin[y] && j+sdr <= jmax[v]));
	  
	/* determine min d (dn) and max d (dx) that are valid for v and j */
	dn   = hdmin[v][jp_v + sdr] - sd;
	dx   = hdmax[v][jp_v + sdr] - sd;
	dp_v = dx - hdmin[v][jp_v];
	i    = j-dx+1;

	  /* for each emit mode, update beta[v][jp_v][dp_v] for all valid d = dp_v */
	switch(emitmode) { 
	case EMITPAIR:  /* MP_st */
	  for (d = dx; d >= dn; d--, dp_v--, i++) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v] 
							    + escore));
	  }
	  break;

	case EMITLEFT:  /* ML_st, IL_st */
	  for (d = dx; d >= dn; d--, dp_v--, i++) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v + sdr][dp_y + sd] + cm->endsc[v] 
							    + escore));
	  }
	  break;

	case EMITRIGHT:  /* MR_st, IR_st */
	  escore = esc_vAA[v][dsq[j+1]]; /* not dependent on i */
	  for (d = dx; d >= dn; d--, dp_v--) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v + sdr][dp_y + sd] + cm->endsc[v] 
							    + escore));
	  }
	  break;

	case EMITNONE:  /* D_st, S_st, E_st*/
	  for (d = dx; d >= dn; d--, dp_v--) { 
	    ESL_DASSERT1((((d + sd) >= hdmin[v][jp_v + sdr]) && ((d + sd) <= hdmax[v][jp_v + sdr])));
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v + sdr][dp_y + sd] + cm->endsc[v]));
	    break;
	  }
	} /* end of switch over emitmodes */
      } /* end of loop over j */
    } /* end conditional section for dealing w/ v->EL local end transitions */
  } /* end loop over decks v. */

  /*FILE *fp; fp = fopen("out.mx", "w"); cm_ihb_mx_Dump(fp, mx); fclose(fp);*/

  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, EL->EL transitions. (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->el_selfsc));
    }
  }

  Wp = W - hdmin[0][j0-jmin[0]];
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    
    /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
     *                                                = P(S|M) / P(S|R)  
     * in v,j,d coordinates this is:
     * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                = P(S|M) / P(S|R)
     */
    
    for(n = 0; n < cm->nodes; n++) {
      fsc = IMPOSSIBLE;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	    /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
	      printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	    fsc = FLogsum(fsc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
	  }
	}
      }
      /*printf("checking node: %d | sc: %.6f\n", n, fsc);*/
      diff = fsc - alpha[0][j0-jmin[0]][Wp];
      if(diff > 0.01 || diff < -0.01) { /* we're using floats, this should be pretty precise */
	fail_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, fsc, alpha[0][(j0-jmin[0])][Wp], diff);
      }
    }
  }

  /* If not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to W (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to W (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * *** If local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END)) { 
    freturn_sc = IMPOSSIBLE;
    v = cm->M-1;
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      assert(hdmin[v][jp_v] == 0);
      /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
      freturn_sc = FLogsum(freturn_sc, (beta[v][jp_v][0]));
    }
  }
  else { /* return_sc = P(S|M) / P(S|R) from Inside() */
    freturn_sc = alpha[0][(j0-jmin[0])][Wp];
  }

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tFast_FOutside_b_jd_me() sc : %f\n", freturn_sc));
  else
    ESL_DPRINTF1(("\tFast_FOutside_b_jd_me() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  return freturn_sc;
}

/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_FASTALIGN_BENCHMARK
/* gcc -o benchmark-fastalign -g -O2 -I. -L. -I../easel -L../easel -DIMPL_FASTALIGN_BENCHMARK cm_fastalign.c -linfernal -leasel -lm
 * ./benchmark-fastalign <cmfile>
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
  //{ "--scan",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "run in scan mode, not alignment mode", 0 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use posterior sums during HMM band calculation (widens bands)", 0 },
  { "--dlev",    eslARG_INT,    "0",   NULL, "0<=n<=3",NULL,NULL,NULL, "set verbosity of debugging print statements to <n>", 0 },
  { "--hmmcheck",eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "check that HMM posteriors are correctly calc'ed", 0 },
  { "--cmcheck", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "check that HMM posteriors are correctly calc'ed", 0 },
  { "--fpost",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute fast float HMM banded Inside/Outside alignment algs", 0 },
  { "--ipost",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute fast int   HMM banded Inside/Outside alignment algs", 0 },
  { "--ofpost",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute slow float HMM banded Inside/Outside alignment algs", 0 },
  { "--oipost",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute slow int   HMM banded Inside/Outside alignment algs", 0 },
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
  float           sc;
  /*float           bsc;
    int             v;
    int             b;*/
  char           *cmfile = esl_opt_GetArg(go, 1);
  CMFILE         *cmfp;	    /* open input CM file stream */
  int             L;        /* length of sequence */
  int             do_random;
  int             N = esl_opt_GetInteger(go, "-N");
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, either randomly created, or emitted from CM (if -e) */
  int             do_align;
  CM_FHB_MX      *mx;
  CM_IHB_MX      *iout_mx, *iins_mx;
  CM_FHB_MX      *fout_mx, *fins_mx;
  /* FILE *fpfast;
     FILE *fpslow; */
  /* Parsetree_t    *slowtr, *fasttr; */

  int             ***oialpha;    
  int             ***oibeta;     
  float           ***ofalpha;    
  float           ***ofbeta;     

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  do_align = TRUE;
  /* if(esl_opt_GetBoolean(go, "--scan")) do_align = FALSE; */

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  /* determine sequence length */
  if(esl_opt_IsDefault(go, "-L")) L = cm->clen;      
  else                            L = esl_opt_GetInteger(go, "-L");

  /* configure CM for HMM banded alignment */
  //cm->config_opts |= CM_CONFIG_QDB;
  
  cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  cm->align_opts  |= CM_ALIGN_TIME;
  if(do_align) { 
    cm->align_opts  |= CM_ALIGN_NOSMALL;
    cm->align_opts  |= CM_ALIGN_HBANDED;
    if(esl_opt_GetBoolean(go, "--sums")) cm->align_opts |= CM_ALIGN_SUMS;
  }
  /*  else {
    cm->search_opts  |= CM_SEARCH_HBANDED;
    cm->search_opts  |= CM_SEARCH_HMMSCANBANDS;
    if(esl_opt_GetBoolean(go, "--sums")) cm->search_opts |= CM_SEARCH_SUMS;
    }*/
  if(esl_opt_GetBoolean(go, "-l")) { 
    cm->config_opts  |= CM_CONFIG_LOCAL;
    cm->config_opts  |= CM_CONFIG_HMMLOCAL;
    cm->config_opts  |= CM_CONFIG_HMMEL;
  }
  if(esl_opt_GetBoolean(go, "--hmmcheck")) cm->align_opts |= CM_ALIGN_CHECKFB;

  ConfigCM(cm, NULL, NULL);
  CP9Bands_t *cp9b;
  cp9b = AllocCP9Bands(cm, cm->cp9);

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

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
  if(esl_opt_GetBoolean(go, "--ipost")) { 
    iins_mx = cm_ihb_mx_Create(cm->M);
    iout_mx = cm_ihb_mx_Create(cm->M);
  }
  if(esl_opt_GetBoolean(go, "--fpost")) { 
    fins_mx = cm_fhb_mx_Create(cm->M);
    fout_mx = cm_fhb_mx_Create(cm->M);
  }
  int do_check = esl_opt_GetBoolean(go, "--cmcheck");
  for (i = 0; i < N; i++)
    {
      L = seqs_to_aln->sq[i]->n;

      esl_stopwatch_Start(w);
      cp9_Seq2Bands(cm, seqs_to_aln->sq[i]->dsq, 1, L, cp9b, 0);
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "Exptl Band calc", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      esl_stopwatch_Start(w);
      sc = FastCYKInsideAlignHB(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, NULL, cp9b, mx);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKInsideAlignHB (): ", sc);
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

      if(esl_opt_GetBoolean(go, "--ipost")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FastIInsideAlignHB(cm, seqs_to_aln->sq[i]->dsq, 1, L, cp9b, iins_mx);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastIInsideAlignHB (): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FastIOutsideAlignHB(cm, seqs_to_aln->sq[i]->dsq, 1, L, cp9b, iout_mx, iins_mx, do_check);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastIOutsideAlignHB (): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      if(esl_opt_GetBoolean(go, "--fpost")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FastFInsideAlignHB(cm, seqs_to_aln->sq[i]->dsq, 1, L, cp9b, fins_mx);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastFInsideAlignHB (): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FastFOutsideAlignHB(cm, seqs_to_aln->sq[i]->dsq, 1, L, cp9b, fout_mx, fins_mx, do_check);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastFOutsideAlignHB (): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      /* do old Inside/Outside if requested */
      if(esl_opt_GetBoolean(go, "--oipost")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = IInside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			     TRUE,	    /* save full alpha so we can run outside */
			     NULL, &oialpha, /* fill alpha, and return it, needed for IOutside() */
			     NULL, NULL,    /* manage your own deckpool, I don't want it */
			     esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			     cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "IInside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = IOutside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			      TRUE,	        /* save full beta */
			      NULL, &oibeta,	/* fill beta, and return it, needed for ICMPosterior() */
			      NULL, NULL,	/* manage your own deckpool, I don't want it */
			      esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			      oialpha, &oialpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
			      do_check,      /* TRUE to check Outside probs agree with Inside */
			      cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "IOutside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      /* do old Inside/Outside if requested */
      if(esl_opt_GetBoolean(go, "--oipost")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FInside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			     TRUE,	    /* save full alpha so we can run outside */
			     NULL, &ofalpha, /* fill alpha, and return it, needed for IOutside() */
			     NULL, NULL,    /* manage your own deckpool, I don't want it */
			     esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			     cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "FInside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FOutside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			      TRUE,	        /* save full beta */
			      NULL, &ofbeta,	/* fill beta, and return it, needed for ICMPosterior() */
			      NULL, NULL,	/* manage your own deckpool, I don't want it */
			      esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			      ofalpha, &ofalpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
			      do_check,      /* TRUE to check Outside probs agree with Inside */
			      cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "FOutside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }
      printf("\n");
    }

  cm_fhb_mx_Destroy(mx);
  if(esl_opt_GetBoolean(go, "--fpost")) { 
    cm_fhb_mx_Destroy(fins_mx);
    cm_fhb_mx_Destroy(fout_mx);
  }
  if(esl_opt_GetBoolean(go, "--ipost")) { 
    cm_ihb_mx_Destroy(iins_mx);
    cm_ihb_mx_Destroy(iout_mx);
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
  return 0; /* NEVERREACHED */
}
#endif /*IMPL_FASTALIGN_BENCHMARK*/

#if 0	
/* old versions of funcs */
/*
 * Function: OLDFastIOutsideAlignHB()
 * Date:     EPN, Wed Nov  7 15:24:02 2007 [updated]
 *
 * Purpose:  Run the outside algorithm using bands in the j and d 
 *           dimension from obtained from an HMM forwards-backwards 
 *           run. This function is memory efficient in the j AND d 
 *           dimension.
 *
 *           The full sequence is aligned against the full model.
 *           This function complements FastIInsideAlignHB().
 *
 *           If <do_check> is TRUE and the CM is not in local mode, 
 *           this function checks that the outside calculations are
 *           consistent with the inside calculations (in ins_mx).
 *           This check is described in comments towards the end of 
 *           the function. 
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           cp9b      - HMM bands for <dsq> and <cm>
 *                       
 * Returns:  log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *           ends are on (see *** comment towards end of function).
 */
float 
OLDFastIOutsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b,
		       CM_IHB_MX *mx, CM_IHB_MX *ins_mx, int do_check)
{
  int      status;
  int      v,y,z;	/* indices for states */
  int      j,d,i,k;	/* indices in sequence dimensions */
  int      isc;     	/* a temporary variable holding an int score */
  float    fsc;     	/* a temporary variable holding a float score */
  int      iescore;	/* an emission score, tmp variable */
  int      W;		/* subsequence length */
  int      voffset;	/* index of v in t_v(y) transition scores */
  int      jp;		/* j': relative position in the subsequence  */
  int      bsc;		/* total score for using local begin states */
  int      ireturn_sc;  /* P(S|M)/P(S|R), a scaled int*/
  float    freturn_sc;  /* P(S|M)/P(S|R), a float (Scorified ireturn_sc) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;
  int      fail_flag = FALSE; /* set to TRUE if do_check and we see a problem */
  int     **esc_vAA = cm->ioesc;

  int       *jmin = cp9b->jmin;
  int       *jmax = cp9b->jmax;
  int     **hdmin = cp9b->hdmin;
  int     **hdmax = cp9b->hdmax;

  /* Contract check */
  if(dsq == NULL) cm_Fail("ERROR in IOutside_b_jd_me(), dsq is NULL.\n");
  if(i0 != 1)     cm_Fail("ERROR: IOutside_b_jd_me requires that i0 be 1. This function is not set up for subsequence alignment\n");
  if (mx == NULL)     cm_Fail("IOutside_b_jd_me(), mx is NULL.\n");
  if (ins_mx == NULL) cm_Fail("IOutside_b_jd_me(), ins_mx is NULL.\n");
  if (cm->flags & CMH_LOCAL_END && do_check == FALSE) cm_Fail("IOutside_b_jd_me(), do_check = TRUE, but CMH_LOCAL_END flag raised, we can't check with local ends.n");
  /* Code for checking doesn't apply in local mode. See below. */

  /* Allocations and initializations
   */
  bsc = -INFTY;
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */

  /* grow the matrix based on the current sequence and bands */
  cm_ihb_mx_GrowTo(mx, cp9b);
  int ***beta  = mx->dp;
  int ***alpha = ins_mx->dp;

  /* precalculate data that we can */
  int *el_scA;
  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(int) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->iel_selfsc * d;

  int *yvalidA; 
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, 0);

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   * EL deck has no bands as currently implemented. Set all cells to -INFTY;
   */
  int **beta_el = NULL;
  if (cm->flags & CMH_LOCAL_END) {
    beta_el = Ialloc_vjd_deck(W, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++) beta_el[j][d] = -INFTY;
    }
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }
  /* assign the ptr for EL in our DP matrix to the EL deck, it should be NULL */
  assert(beta[cm->M] == NULL);
  beta[cm->M] = beta_el;

  /* initialize all cells of the matrix to -INFTY */
  esl_vec_ISet(beta[0][0], mx->ncells_valid, -INFTY);
  /* now set beta[0][j0][W] to 0., all parses must end there */
  jp_v = j0 - jmin[0];
  Wp = W - hdmin[0][jp_v];
  assert(W >= hdmin[0][jp_v]);
  assert(W <= hdmax[0][jp_v]);
  beta[0][jp_v][Wp] = 0.;

  /* If we can do a local begin into v, overwrite -INFY with the local begin score. 
   * By definition, beta[0][j0][W] == 0.
   */ 
  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (v = 1; v < cm->M; v++) {
      if(cm->ibeginsc[v] != -INFTY) { 
	if((j0 >= jmin[v]) && (j0 <= jmax[v])) {
	  jp_v = j0 - jmin[v];
	  if((W >= hdmin[v][jp_v]) && W <= hdmax[v][jp_v]) {
	    Wp = W - hdmin[v][jp_v];
	    beta[v][jp_v][Wp] = cm->ibeginsc[v];
	  }
	}
      }
    }
  }
  
  /* Main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    /* main recursion: reorganized relative to FOutside() for simplification of
     * band-related issues.
     */
    for (j = jmax[v]; j >= jmin[v]; j--)
      {
	jp_v = j - jmin[v];
	for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--)
	  {
	    i = j-d+1;
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	    
	    if (cm->stid[v] == BEGL_S) 
	      {
		y = cm->plast[v];	/* the parent bifurcation    */
		z = cm->cnum[y];	/* the other (right) S state */
		jp_y = j - jmin[y];
		jp_z = j - jmin[z];
		
		/* Find the first k value that implies a valid cell in the y and z decks.
		 * This k must satisfy the following 8 inequalities (some may be redundant):
		 * NOTE: these are different from those in Inside() (for one thing, v and y
		 *       (BEGL_S and BIF_B here respectively) are switched relative to Inside.
		 *
		 * (1) k <= jmax[y] - j;
		 * (2) k >= jmin[y] - j;
		 * (3) k <= jmax[z] - j;
		 * (4) k >= jmin[z] - j;
		 *     1 and 2 guarantee (j+k) is within state y's j band
		 *     3 and 4 guarantee (j+k) is within state z's j band
		 *
		 * (5) k >= hdmin[y][j-jmin[y]+k] - d;
		 * (6) k <= hdmax[y][j-jmin[y]+k] - d; 
		 *     5 and 6 guarantee k+d is within y's j=(j+k), d band
		 *
		 * (7) k >= hdmin[z][j-jmin[z]+k];
		 * (8) k <= hdmax[z][j-jmin[z]+k]; 
		 *     5 and 6 guarantee k is within state z's j=(j+k) d band
		 * 
		 * Note, below:
		 * tmp_jmin = MAX(jmin[y], jmin[z];
		 * tmp_jmax = MIN(jmax[y], jmax[z];
		 */
		tmp_jmin = (jmin[y] > jmin[z]) ? jmin[y] : jmin[z];
		tmp_jmax = (jmax[y] < jmax[z]) ? jmax[y] : jmax[z];
		
		kmin = tmp_jmin - j;
		kmax = tmp_jmax - j;
		/* kmin and kmax satisfy inequalities (1-4) */
		/* RHS of inequalities 5-8 are dependent on k, so we check
		 * for these within the next for loop.
		 */
		for(k = kmin; k <= kmax; k++) {
		  if(k < (hdmin[y][jp_y+k] - d) || k > (hdmax[y][jp_y+k] - d)) continue; 
		  /* above line continues if inequality 5 or 6 is violated */
		  if(k < (hdmin[z][jp_z+k]) || k > (hdmax[z][jp_z+k])) continue; 
		  /* above line continues if inequality 7 or 8 is violated */
		  
		  /* if we get here for current k, all 8 inequalities have been satisified 
		   * so we know the cells corresponding to the platonic 
		   * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
		   * alpha[z][j+k][k] are all within the bands. These
		   * cells correspond to beta[v][jp_v][dp_v], 
		   * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
		   * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
		   */
		  kp_z = k-hdmin[z][jp_z+k];
		  dp_y = d-hdmin[y][jp_y+k];
		  beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] 
								      + alpha[z][jp_z+k][kp_z]));
		}
	      }
	    else if (cm->stid[v] == BEGR_S) {
	      y = cm->plast[v];	/* the parent bifurcation    */
	      z = cm->cfirst[y];	/* the other (left) S state */
	      
	      jp_y = j - jmin[y];
	      jp_z = j - jmin[z];
	      
	      /* For j to be valid for state y: 
	       * jmin[y] >= j >= jmax[y]
	       * These are independent of k so we check outside of k loop below
	       * For j to be valid for state z: *
	       * (jmin[z] + d) >= j >= (jmax[z] + d)
	       */
	      if(j < jmin[y] || j > jmax[y]) continue;
	      if((j < (jmin[z] + d)) || (j > (jmax[z]+d))) continue;
	      
	      /* Find the first k value that implies a valid cell in the y and z decks.
	       * This k must satisfy the following 4 inequalities (some may be redundant):
	       * NOTE: these are different from those in Inside() (for one thing, v and y
	       *       (BEGR_S and BIF_B here respectively) are switched relative to Inside.
	       *
	       * (1) k >= hdmin[y][j-jmin[y]] - d;
	       * (2) k <= hdmax[y][j-jmin[y]] - d;
	       *     1 and 2 guarantee (d+k) is within state y's j=(j) d band
	       *
	       * (3) k >= hdmin[z][j-jmin[z]-d];
	       * (4) k <= hdmax[z][j-jmin[z]-d];
	       *     3 and 4 guarantee k is within z's j=(j-d), d band
	       *
	       */
	      kmin = ((hdmin[y][jp_y]-d) > (hdmin[z][jp_z-d])) ? (hdmin[y][jp_y]-d) : (hdmin[z][jp_z-d]);
	      /* kmin satisfies inequalities (1) and (3) */
	      kmax = ((hdmax[y][jp_y]-d) < (hdmax[z][jp_z-d])) ? (hdmax[y][jp_y]-d) : (hdmax[z][jp_z-d]);
	      /* kmax satisfies inequalities (2) and (4) */
	      
	      for(k = kmin; k <= kmax; k++) { 
		/* for current k, all 4 inequalities have been satisified 
		 * so we know the cells corresponding to the platonic 
		 * matrix cells beta[v][j][d], beta[y][j][d+k], and
		 * alpha[z][j-d][k] are all within the bands. These
		 * cells correspond to beta[v][jp_v][dp_v], 
		 * beta[y][jp_y+k][d-hdmin[y][jp_y]+k],
		 * and alpha[z][jp_z-d][k-hdmin[z][jp_z-d]];
		 */
		kp_z = k-hdmin[z][jp_z-d];
		dp_y = d-hdmin[y][jp_y];
		beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
								    + alpha[z][jp_z-d][kp_z]));
	      }
	    }
	    else { 
	      for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
		voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

		switch(cm->sttype[y]) {
		case MP_st: 
		  if (j == j0 || d == j) continue; /* boundary condition */
		  if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if ((d+2) < hdmin[y][(jp_y+1)] || (d+2) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
		   * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][jp_y+1];  /* d index for state y */
		  iescore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		  beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
								      + cm->itsc[y][voffset] + iescore));
		  break;
		  
		case ML_st:
		case IL_st: 
		  if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
		  if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if ((d+1) < hdmin[y][jp_y] || (d+1) > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
		   * corresponding to alpha[y][j][d+1] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][jp_y];  /* d index for state y */
		  iescore = esc_vAA[y][dsq[i-1]];
		  beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
								      + cm->itsc[y][voffset] + iescore));
		  break;
		  
		case MR_st:
		case IR_st:
		  if (j == j0) continue;
		  if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if ((d+1) < hdmin[y][(jp_y+1)] || (d+1) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
		   * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][(jp_y+1)];  /* d index for state y */
		  iescore = esc_vAA[y][dsq[j+1]];
		  /*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
		  beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
								      + cm->itsc[y][voffset] + iescore));
		  break;
		  
		case S_st:
		case E_st:
		case D_st:
		  if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if (d < hdmin[y][jp_y] || d > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
		   * corresponding to alpha[y][j][d] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][jp_y];  /* d index for state y */
		  beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->itsc[y][voffset])); 
		  break;
		  
		default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
		}/* end switch over states*/
	      } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	      if (beta[v][jp_v][dp_v] < -INFTY) beta[v][jp_v][dp_v] = -INFTY;
	    }	/* ends else entered for non-BEGL/BEGR states*/	
	  } /* ends loop over d. We know all beta[v][j][d] in this row j*/
      } /* end loop over jp. We know the beta's for the whole deck.*/
    /* Deal with local alignment end transitions v->EL
     * (EL = deck at M.)
     */
    if (cm->iendsc[v] != -INFTY) {
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j - jmin[v];
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { 
	  i = j-d+1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    if (j == j0 || d == j) continue; /* boundary condition */
	    if (((j+1) > jmax[v]) || ((d+2) > hdmax[v][jp_v])) continue; /*boundary condition*/
	    iescore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+2] + cm->iendsc[v] 
							    + iescore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == j) continue;	
	    if ((d+1) > hdmax[v][jp_v]) continue; /*boundary condition*/
	    iescore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v+1] + cm->iendsc[v] 
							    + iescore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == j0) continue;
	    if (((j+1) > jmax[v]) || ((d+1) > hdmax[v][jp_v])) continue; /*boundary condition*/
	    iescore = esc_vAA[y][dsq[j+1]];
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+1] + cm->iendsc[v] 
							    + iescore));
	    break;
	  case S_st:
	  case D_st:
	  case E_st:
	    beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v] + cm->iendsc[v]));
							    //+ iescore));
	    break;
	  case B_st:  
	  default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	    /* note that although B is a valid vend for a segment we'd do
	       outside on, B->EL is set to be impossible, by the local alignment
	       config. There's no point in having a B->EL because B is a nonemitter
	       (indeed, it would introduce an alignment ambiguity). The same
	       alignment case is handled by the X->EL transition where X is the
	       parent consensus state (S, MP, ML, or MR) above the B. Thus,
	       this code is relying on the (cm->iendsc[v] != -INFTY) test, above,
	       to make sure the sttype[vend]=B case gets into this switch.
	    */
	  } /* end switch over parent state type v */
	} /* end inner loop over d */
      } /* end outer loop over jp */
    } /* end conditional section for dealing w/ v->EL local end transitions */
  } /* end loop over decks v. */
  /*FILE *fp
    fp = fopen("out.mx", "w");
    cm_ihb_mx_Dump(fp, mx);
    fclose(fp);*/


  /* EPN left the below SRE block out */
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->iel_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE JD, printing beta\n");
    debug_print_alpha_banded_jd(beta, cm, L, jmin, jmax, hdmin, hdmax);
    printf("OUTSIDE JD, done printing beta\n");
  */

  Wp = W - hdmin[0][j0-jmin[0]];
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    
    /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
     *                                                = P(S|M) / P(S|R)  
     * in v,j,d coordinates this is:
     * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                = P(S|M) / P(S|R)
     */
    
    for(n = 0; n < cm->nodes; n++) {
      isc = -INFTY;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	    /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
	      printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	    isc = ILogsum(isc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
	  }
	}
      }
      fsc = Scorify(isc);
      /*printf("checking node: %d | sc: %.6f\n", n, fsc);*/
      diff = fsc - (Scorify(alpha[0][j0-jmin[0]][Wp]));
      if(diff > 0.01 || diff < -0.01) {
	fail_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, fsc, Scorify(alpha[0][(j0-jmin[0])][Wp]), diff);
      }
    }
  }

  /* If not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to W (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to W (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * *** If local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END)) { 
    ireturn_sc = -INFTY;
      v = cm->M-1;
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	assert(hdmin[v][jp_v] == 0);
	/* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	ireturn_sc = ILogsum(ireturn_sc, (beta[v][jp_v][0]));
      }
      freturn_sc = Scorify(ireturn_sc);
  }
  else { /* return_sc = P(S|M) / P(S|R) from Inside() */
    freturn_sc = Scorify(alpha[0][(j0-jmin[0])][Wp]);
  }

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tFast_IOutside_b_jd_me() sc : %f\n", freturn_sc));
  else
    ESL_DPRINTF1(("\tFast_IOutside_b_jd_me() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  return freturn_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/*
 * Function: FastFOutsideAlignHB()
 * Date:     EPN, Thu Nov  8 18:40:05 2007
 *
 * Purpose:  Run the outside algorithm using bands in the j and d 
 *           dimension from obtained from an HMM forwards-backwards 
 *           run. This function is memory efficient in the j AND d 
 *           dimension. Uses Float scores, not ints. 
 *
 *           The full sequence is aligned against the full model.
 *           This function complements FastIInsideAlignHB().
 *
 *           If <do_check> is TRUE and the CM is not in local mode, 
 *           this function checks that the outside calculations are
 *           consistent with the inside calculations (in ins_mx).
 *           This check is described in comments towards the end of 
 *           the function. 
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           cp9b      - HMM bands for <dsq> and <cm>
 *                       
 * Returns:  log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *           ends are on (see *** comment towards end of function).
 */
float 
FastFOutsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b,
		    CM_FHB_MX *mx, CM_FHB_MX *ins_mx, int do_check)
{
  int      status;
  int      v,y,z;	/* indices for states */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    fsc;     	/* a temporary variable holding a float score */
  float    escore;	/* an emission score, tmp variable */
  int      W;		/* subsequence length */
  int      voffset;	/* index of v in t_v(y) transition scores */
  int      jp;		/* j': relative position in the subsequence  */
  float    bsc;		/* total score for using local begin states */
  float    freturn_sc;  /* P(S|M)/P(S|R), a float (Scorified ireturn_sc) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;
  int      fail_flag = FALSE; /* set to TRUE if do_check and we see a problem */
  float   **esc_vAA = cm->oesc;
  
  int       *jmin = cp9b->jmin;
  int       *jmax = cp9b->jmax;
  int     **hdmin = cp9b->hdmin;
  int     **hdmax = cp9b->hdmax;

  /* Contract check */
  if(dsq == NULL) cm_Fail("ERROR in IOutside_b_jd_me(), dsq is NULL.\n");
  if(i0 != 1)     cm_Fail("ERROR: IOutside_b_jd_me requires that i0 be 1. This function is not set up for subsequence alignment\n");
  if (mx == NULL)     cm_Fail("IOutside_b_jd_me(), mx is NULL.\n");
  if (ins_mx == NULL) cm_Fail("IOutside_b_jd_me(), ins_mx is NULL.\n");
  if (cm->flags & CMH_LOCAL_END && do_check == FALSE) cm_Fail("IOutside_b_jd_me(), do_check = TRUE, but CMH_LOCAL_END flag raised, we can't check with local ends.n");
  /* Code for checking doesn't apply in local mode. See below. */

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */

  /* grow the matrix based on the current sequence and bands */
  cm_fhb_mx_GrowTo(mx, cp9b);
  float ***beta  = mx->dp;
  float ***alpha = ins_mx->dp;

  /* precalculate data that we can */
  float *el_scA;
  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  int *yvalidA; 
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, 0);

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   * EL deck has no bands as currently implemented. Set all cells to IMPOSSIBLE;
   */
  float **beta_el = NULL;
  if (cm->flags & CMH_LOCAL_END) {
    beta_el = alloc_vjd_deck(W, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++) beta_el[j][d] = IMPOSSIBLE;
    }
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }
  /* assign the ptr for EL in our DP matrix to the EL deck, it should be NULL */
  assert(beta[cm->M] == NULL);
  beta[cm->M] = beta_el;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(beta[0][0], mx->ncells_valid, IMPOSSIBLE);
  /* now set beta[0][j0][W] to 0., all parses must end there */
  jp_v = j0 - jmin[0];
  Wp = W - hdmin[0][jp_v];
  assert(W >= hdmin[0][jp_v]);
  assert(W <= hdmax[0][jp_v]);
  beta[0][jp_v][Wp] = 0.;

  /* If we can do a local begin into v, overwrite -INFY with the local begin score. 
   * By definition, beta[0][j0][W] == 0.
   */ 
  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (v = 1; v < cm->M; v++) {
      if(NOT_IMPOSSIBLE(cm->beginsc[v])) { 
	if((j0 >= jmin[v]) && (j0 <= jmax[v])) {
	  jp_v = j0 - jmin[v];
	  if((W >= hdmin[v][jp_v]) && W <= hdmax[v][jp_v]) {
	    Wp = W - hdmin[v][jp_v];
	    beta[v][jp_v][Wp] = cm->beginsc[v];
	  }
	}
      }
    }
  }
  
  /* Main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    /* main recursion: reorganized relative to FOutside() for simplification of
     * band-related issues.
     */
    for (j = jmax[v]; j >= jmin[v]; j--)
      {
	jp_v = j - jmin[v];
	for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--)
	  {
	    i = j-d+1;
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	    
	    if (cm->stid[v] == BEGL_S) 
	      {
		y = cm->plast[v];	/* the parent bifurcation    */
		z = cm->cnum[y];	/* the other (right) S state */
		jp_y = j - jmin[y];
		jp_z = j - jmin[z];
		
		/* Find the first k value that implies a valid cell in the y and z decks.
		 * This k must satisfy the following 8 inequalities (some may be redundant):
		 * NOTE: these are different from those in Inside() (for one thing, v and y
		 *       (BEGL_S and BIF_B here respectively) are switched relative to Inside.
		 *
		 * (1) k <= jmax[y] - j;
		 * (2) k >= jmin[y] - j;
		 * (3) k <= jmax[z] - j;
		 * (4) k >= jmin[z] - j;
		 *     1 and 2 guarantee (j+k) is within state y's j band
		 *     3 and 4 guarantee (j+k) is within state z's j band
		 *
		 * (5) k >= hdmin[y][j-jmin[y]+k] - d;
		 * (6) k <= hdmax[y][j-jmin[y]+k] - d; 
		 *     5 and 6 guarantee k+d is within y's j=(j+k), d band
		 *
		 * (7) k >= hdmin[z][j-jmin[z]+k];
		 * (8) k <= hdmax[z][j-jmin[z]+k]; 
		 *     5 and 6 guarantee k is within state z's j=(j+k) d band
		 * 
		 * Note, below:
		 * tmp_jmin = MAX(jmin[y], jmin[z];
		 * tmp_jmax = MIN(jmax[y], jmax[z];
		 */
		tmp_jmin = (jmin[y] > jmin[z]) ? jmin[y] : jmin[z];
		tmp_jmax = (jmax[y] < jmax[z]) ? jmax[y] : jmax[z];
		
		kmin = tmp_jmin - j;
		kmax = tmp_jmax - j;
		/* kmin and kmax satisfy inequalities (1-4) */
		/* RHS of inequalities 5-8 are dependent on k, so we check
		 * for these within the next for loop.
		 */
		for(k = kmin; k <= kmax; k++) {
		  if(k < (hdmin[y][jp_y+k] - d) || k > (hdmax[y][jp_y+k] - d)) continue; 
		  /* above line continues if inequality 5 or 6 is violated */
		  if(k < (hdmin[z][jp_z+k]) || k > (hdmax[z][jp_z+k])) continue; 
		  /* above line continues if inequality 7 or 8 is violated */
		  
		  /* if we get here for current k, all 8 inequalities have been satisified 
		   * so we know the cells corresponding to the platonic 
		   * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
		   * alpha[z][j+k][k] are all within the bands. These
		   * cells correspond to beta[v][jp_v][dp_v], 
		   * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
		   * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
		   */
		  kp_z = k-hdmin[z][jp_z+k];
		  dp_y = d-hdmin[y][jp_y+k];
		  beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] 
								      + alpha[z][jp_z+k][kp_z]));
		}
	      }
	    else if (cm->stid[v] == BEGR_S) {
	      y = cm->plast[v];	/* the parent bifurcation    */
	      z = cm->cfirst[y];	/* the other (left) S state */
	      
	      jp_y = j - jmin[y];
	      jp_z = j - jmin[z];
	      
	      /* For j to be valid for state y: 
	       * jmin[y] >= j >= jmax[y]
	       * These are independent of k so we check outside of k loop below
	       * For j to be valid for state z: *
	       * (jmin[z] + d) >= j >= (jmax[z] + d)
	       */
	      if(j < jmin[y] || j > jmax[y]) continue;
	      if((j < (jmin[z] + d)) || (j > (jmax[z]+d))) continue;
	      
	      /* Find the first k value that implies a valid cell in the y and z decks.
	       * This k must satisfy the following 4 inequalities (some may be redundant):
	       * NOTE: these are different from those in Inside() (for one thing, v and y
	       *       (BEGR_S and BIF_B here respectively) are switched relative to Inside.
	       *
	       * (1) k >= hdmin[y][j-jmin[y]] - d;
	       * (2) k <= hdmax[y][j-jmin[y]] - d;
	       *     1 and 2 guarantee (d+k) is within state y's j=(j) d band
	       *
	       * (3) k >= hdmin[z][j-jmin[z]-d];
	       * (4) k <= hdmax[z][j-jmin[z]-d];
	       *     3 and 4 guarantee k is within z's j=(j-d), d band
	       *
	       */
	      kmin = ((hdmin[y][jp_y]-d) > (hdmin[z][jp_z-d])) ? (hdmin[y][jp_y]-d) : (hdmin[z][jp_z-d]);
	      /* kmin satisfies inequalities (1) and (3) */
	      kmax = ((hdmax[y][jp_y]-d) < (hdmax[z][jp_z-d])) ? (hdmax[y][jp_y]-d) : (hdmax[z][jp_z-d]);
	      /* kmax satisfies inequalities (2) and (4) */
	      
	      for(k = kmin; k <= kmax; k++) { 
		/* for current k, all 4 inequalities have been satisified 
		 * so we know the cells corresponding to the platonic 
		 * matrix cells beta[v][j][d], beta[y][j][d+k], and
		 * alpha[z][j-d][k] are all within the bands. These
		 * cells correspond to beta[v][jp_v][dp_v], 
		 * beta[y][jp_y+k][d-hdmin[y][jp_y]+k],
		 * and alpha[z][jp_z-d][k-hdmin[z][jp_z-d]];
		 */
		kp_z = k-hdmin[z][jp_z-d];
		dp_y = d-hdmin[y][jp_y];
		beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
								    + alpha[z][jp_z-d][kp_z]));
	      }
	    }
	    else { 
	      for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
		voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

		switch(cm->sttype[y]) {
		case MP_st: 
		  if (j == j0 || d == j) continue; /* boundary condition */
		  if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if ((d+2) < hdmin[y][(jp_y+1)] || (d+2) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
		   * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][jp_y+1];  /* d index for state y */
		  escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		  beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
								      + cm->tsc[y][voffset] + escore));
		  break;
		  
		case ML_st:
		case IL_st: 
		  if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
		  if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if ((d+1) < hdmin[y][jp_y] || (d+1) > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
		   * corresponding to alpha[y][j][d+1] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][jp_y];  /* d index for state y */
		  escore = esc_vAA[y][dsq[i-1]];
		  beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
								      + cm->tsc[y][voffset] + escore));
		  break;
		  
		case MR_st:
		case IR_st:
		  if (j == j0) continue;
		  if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if ((d+1) < hdmin[y][(jp_y+1)] || (d+1) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
		   * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][(jp_y+1)];  /* d index for state y */
		  escore = esc_vAA[y][dsq[j+1]];
		  /*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
		  beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
								      + cm->tsc[y][voffset] + escore));
		  break;
		  
		case S_st:
		case E_st:
		case D_st:
		  if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if (d < hdmin[y][jp_y] || d > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
		  /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
		   * corresponding to alpha[y][j][d] in the platonic matrix.
		   */
		  dp_y = d - hdmin[y][jp_y];  /* d index for state y */
		  beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->tsc[y][voffset])); 
		  break;
		  
		default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
		}/* end switch over states*/
	      } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	      if (beta[v][jp_v][dp_v] < IMPOSSIBLE) beta[v][jp_v][dp_v] = IMPOSSIBLE;
	    }	/* ends else entered for non-BEGL/BEGR states*/	
	  } /* ends loop over d. We know all beta[v][j][d] in this row j*/
      } /* end loop over jp. We know the beta's for the whole deck.*/
    /* Deal with local alignment end transitions v->EL
     * (EL = deck at M.)
     */
    if (NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j - jmin[v];
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { 
	  i = j-d+1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    if (j == j0 || d == j) continue; /* boundary condition */
	    if (((j+1) > jmax[v]) || ((d+2) > hdmax[v][jp_v])) continue; /*boundary condition*/
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+2] + cm->endsc[v] 
							    + escore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == j) continue;	
	    if ((d+1) > hdmax[v][jp_v]) continue; /*boundary condition*/
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v+1] + cm->endsc[v] 
							    + escore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == j0) continue;
	    if (((j+1) > jmax[v]) || ((d+1) > hdmax[v][jp_v])) continue; /*boundary condition*/
	    escore = esc_vAA[y][dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+1] + cm->endsc[v] 
							    + escore));
	    break;
	  case S_st:
	  case D_st:
	  case E_st:
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v] + cm->endsc[v]));
							    //+ escore));
	    break;
	  case B_st:  
	  default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	    /* note that although B is a valid vend for a segment we'd do
	       outside on, B->EL is set to be impossible, by the local alignment
	       config. There's no point in having a B->EL because B is a nonemitter
	       (indeed, it would introduce an alignment ambiguity). The same
	       alignment case is handled by the X->EL transition where X is the
	       parent consensus state (S, MP, ML, or MR) above the B. Thus,
	       this code is relying on the NOT_IMPOSSIBLE(cm->endsc[v]) test, above,
	       to make sure the sttype[vend]=B case gets into this switch.
	    */
	  } /* end switch over parent state type v */
	} /* end inner loop over d */
      } /* end outer loop over jp */
    } /* end conditional section for dealing w/ v->EL local end transitions */
  } /* end loop over decks v. */
  /*FILE *fp
    fp = fopen("out.mx", "w");
    cm_ihb_mx_Dump(fp, mx);
    fclose(fp);*/


  /* EPN left the below SRE block out */
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->el_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE JD, printing beta\n");
    debug_print_alpha_banded_jd(beta, cm, L, jmin, jmax, hdmin, hdmax);
    printf("OUTSIDE JD, done printing beta\n");
  */

  Wp = W - hdmin[0][j0-jmin[0]];
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    
    /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
     *                                                = P(S|M) / P(S|R)  
     * in v,j,d coordinates this is:
     * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                = P(S|M) / P(S|R)
     */
    
    for(n = 0; n < cm->nodes; n++) {
      fsc = IMPOSSIBLE;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	    /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
	      printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	    fsc = FLogsum(fsc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
	  }
	}
      }
      /*printf("checking node: %d | sc: %.6f\n", n, fsc);*/
      diff = fsc - alpha[0][j0-jmin[0]][Wp];
      if(diff > 0.01 || diff < -0.01) {
	fail_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, fsc, alpha[0][(j0-jmin[0])][Wp], diff);
      }
    }
  }

  /* If not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to W (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to W (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * *** If local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END)) { 
    freturn_sc = IMPOSSIBLE;
      v = cm->M-1;
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	assert(hdmin[v][jp_v] == 0);
	/* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	freturn_sc = FLogsum(freturn_sc, (beta[v][jp_v][0]));
      }
  }
  else { /* return_sc = P(S|M) / P(S|R) from Inside() */
    freturn_sc = alpha[0][(j0-jmin[0])][Wp];
  }

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tFast_IOutside_b_jd_me() sc : %f\n", freturn_sc));
  else
    ESL_DPRINTF1(("\tFast_IOutside_b_jd_me() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  return freturn_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/*
 * Function: FastFInsideAlignHB()
 * Date:     EPN, Thu Nov  8 18:24:41 2007
 *
 * Purpose:  Run the inside algorithm using bands in the j and d 
 *           dimension from obtained from an HMM forwards-backwards 
 *           run. This function is memory efficient in the j AND d 
 *           dimension. Uses Float scores. 
 * 
 *           Synchronized with fast_cyk_inside_align_hb(), see notes
 *           in that function for more details. One major difference
 *           with that function is no shadow matrix is used. Another
 *           is that we use integer log odds score, not floats (for
 *           speed).
 *
 *           The full sequence is aligned against the full model.
 *
 *           This function complements FastFOutsideAlignHB().
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           cp9b      - HMM bands for <dsq> and <cm>
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *                       
 * Returns:  log P(S|M)/P(S|R), as a bit score
 */
float 
FastFInsideAlignHB(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, CM_FHB_MX *mx)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  float    bsc;		/* total score for using local begin states */

  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kn, kx;
  int      dp_y_k;
  int      tmp_jmin, tmp_jmax, kmin, kmax;

  int       *jmin = cp9b->jmin;
  int       *jmax = cp9b->jmax;
  int     **hdmin = cp9b->hdmin;
  int     **hdmax = cp9b->hdmax;
  int       allow_begin = (cm->flags & CMH_LOCAL_BEGIN);

  /* Contract check */
  if(dsq == NULL) cm_Fail("Fast_IInside_b_jd_me(), dsq is NULL.\n");
  if (mx == NULL) cm_Fail("Fast_IInside_b_jd_me(), mx is NULL.\n");

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */

  /* grow the matrix based on the current sequence and bands */
  cm_fhb_mx_GrowTo(mx, cp9b);
  float ***alpha = mx->dp;
  
  float *el_scA;
  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;
  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands)
   */
  int *yvalidA; 
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* Main recursion
   */
  for (v = cm->M-1; v >= 0; v--) {
    /* We've only allocated alpha cells that are within the bands
     * on the j and d dimensions. This means we have to deal
     * with all sorts of offset issues, but we don't have to 
     * waste time setting cells outside the bands to -INFTY;
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
    float const *esc_v = cm->oesc[v]; 
    float const *tsc_v = cm->tsc[v];

    /* initialize all valid cells for state v */
    if (cm->sttype[v] != E_st) {
      if (cm->sttype[v] == B_st) {
	/* initialize all valid cells for state v to IMPOSSIBLE (local ends are impossible for B states) */
	assert(! (NOT_IMPOSSIBLE(cm->endsc[v])));
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) 
	    alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	}
      } else { /* ! B_st && ! E_st */
	/* initialize all valid cells for state v */
	if(NOT_IMPOSSIBLE(cm->endsc[v])) {
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) 
	      alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
	  }
	}
	else { /* cm->endsc[v] == IMPOSSIBLE */
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	      alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	  }
	}
      }
    }

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
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
	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == IR_st) { 
      /* update alpha[v][jp_v][dp_v] cells, for IR states, loop nesting order is:
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
	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) {
      /* for each child state y of v, update alpha[v][jp_v][dp_v] cells */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
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
	  dpn     = dn - hdmin[v][jp_v];
	  dpx     = dx - hdmin[v][jp_v];
	  dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	  
	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y_sdr][dp_y_sd] + tsc));;
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	  }
      default: /* no emission */
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
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
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      /* Any valid j must be within both state v and state z's j band 
       * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
       * constructed, but we'll check anyway. 
       */
      jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
      jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
      /* the main j loop */
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	kn = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
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
	   *
	   * kn and kx were set above (outside (for (dp_v...) loop) that
	   * satisfy 1-4 (b/c 1-4 are d-independent and k-independent)
	   * RHS of inequalities 5 and 6 are dependent on k, so we check
	   * for these within the next for loop.
	   */
	  for(k = kn; k <= kx; k++) { 
	    if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k]) {
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

	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]); 
	    }
	  }
	}
      }
    }				/* finished calculating deck v. */
      
    /* The following loops originally access alpha[v][j0][W] but the index W will be
       in different positions due to the bands */
    
    if(allow_begin) { 
      if(j0 >= jmin[v] && j0 <= jmax[v]) { 
	jp_v = j0 - jmin[v];
	Wp = W - hdmin[v][jp_v];
	if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
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
	  bsc = FLogsum(bsc, (alpha[v][jp_v][Wp] + cm->beginsc[v]));
	}
      }
      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (v == 0) { 
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp = W - hdmin[0][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    alpha[0][jp_v][Wp] = FLogsum(alpha[0][jp_v][Wp], bsc);
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp;
    fp = fopen("ins.mx", "w");
    cm_fhb_mx_Dump(fp, mx);
    fclose(fp);*/

  Wp = W - hdmin[0][j0-jmin[0]];
  sc =     alpha[0][j0-jmin[0]][Wp];

  free(el_scA);
  free(yvalidA);

  return sc;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

#endif 

/* cm_dpalign.c
 * Optimized DP functions for HMM banded and non-banded, non-D&C CM alignment.
 * 
 * HMM banded functions, and their non-optimized analogs: 
 * optimized version                slow, old, reference version      ~speedup
 * ----------------------------     --------------------------------- --------
 * fast_cyk_align_hb()          --> hbandcyk.c:inside_b_jd_me()           1.5 
 * optimal_accuracy_align_hb()  --> NONE
 * fast_alignT_hb()             --> hbandcyk.c:insideT_b_jd_me()          N/A
 * FastAlignHB()                --> hbandcyk.c:CYKInside_b_jd()           N/A
 * FastInsideAlignHB()          --> cm_postprob.c:FInside_b_jd_me()       1.4
 * FastOutsideAlignHB()         --> cm_postprob.c:FOutside_b_jd_me()      1.4
 *
 * FastAlignHB() and fast_alignT_hb() have extra functionality missing
 * from their 0.81 analogs. Specifically they can do either CYK or
 * Holmes/Durbin optimally accurate alignment, the latter of which 
 * was not possible in 0.81.
 * 
 * Speedups are approximate, and based on tests with 2 models, an SSU
 * model and a RNaseP model. Tests were performed with
 * benchmark-fastalign, a standalone executable included at the end of
 * this file that must be separately compiled (version tested was
 * subversion revision 2204, gcc -O2 compiled
 *  EPN, Mon Nov 12 13:41:57 2007).
 * 
 * All functions use a specialized DP matrix, a CM_HB_MX 
 * data structure which only allocates cells within bands 
 * derived from a HMM Forward/Backward alignment of the
 * target sequence. The bands are stored in a CP9Bands_t object,
 * a pointer to which must exist in the cm (CM_t object).
 * 
 * Non-banded, non-D&C alignment functions were rewritten for completeness.
 * These are consistent with their HB counterparts, but require non-banded
 * float matrices. 
 *
 * new version                      old v0.81 version      
 * ----------------------------     ---------------------------------
 * fast_cyk_align()             --> smallcyk.c:inside()    
 * optimal_accuracy_align()     --> NONE
 * fast_alignT()                --> smallcyk.c:insideT()
 * FastAlign()                  --> smallcyk.c:CYKInside()
 * FastInsideAlign()            --> cm_postprob.c:FInside()
 * FastOutsideAlign()           --> cm_postprob.c:FOutside()
 *  
 * Note: the smallcyk.c v0.81 functions are still used for D&C 
 * (divide and conquer) small-memory alignment. Though smallcyk.c
 * has been renamed cm_dpsmall.c
 *
 * Four additional functions exist in this file:
 * CMPosteriorHB()      : combine Inside/Outside matrices --> posterior matrix 
 * CMPosterior()        : non-banded version of CMPosteriorHB()
 * SampleFromInsideHB() : sample a parsetree from HMM banded inside matrix
 * SampleFromInside()   : sample a parsetree from a full inside matrix
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
#include <math.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

#define TSC(s,k) (tsc[(v) * MAXCONNECT + (s)])
#define AMX(j,v,d) (alphap[(j * cm->M * (W+1)) + ((v) * (W+1) + d)])

static float get_femission_score(CM_t *cm, ESL_DSQ *dsq, int v, int i, int j);

/* 
 * Function: fast_cyk_align_hb()
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
 *           A CM_HB_MX DP matrix must be passed in. Only
 *           cells valid within the bands given in the CP9Bands_t <cm->cp9b>
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
 *           will be USE_LOCAL_BEGIN, telling fast_alignT() to check b and
 *           start with a local 0->b entry transition. When inside()
 *           is called on smaller subproblems (v != 0 || i0 > 1 || j0
 *           < L), we're using inside() as an engine in divide &
 *           conquer, and we don't use the overall return score nor
 *           shadow matrices, but we do need allow_begin, b, and bsc for
 *           divide&conquer to sort out where a local begin might be used.
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
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
 *           mx        - the dp matrix, only cells within bands in cm->cp9b will 
 *                       be valid. 
 *           ret_sc    - score of optimal, CYK parsetree 
 *                       
 * Returns: <ret_sc>, <ret_b>, <ret_bsc>, <ret_shadow>, see 'Args'
 * 
 * Throws:  <eslOK> on success.
 *          <eslERANGE> if required CM_HB_MX size exceeds CM_HB_MBLIMIT set in structs.h, in 
 *                      this case, alignment has been aborted, ret_* variables are not valid
 */
int
fast_cyk_align_hb(CM_t *cm, char *errbuf,  ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, void ****ret_shadow,  
		  int allow_begin, int *ret_b, float *ret_bsc, CM_HB_MX *mx, float *ret_sc)
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
  if(dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_cyk_inside_align_hb(), dsq is NULL.\n");
  if (mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_cyk_inside_align_hb(), mx is NULL.\n");
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_cyk_inside_align_hb(), cm->cp9b is NULL.\n");

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrix */
  float ***alpha = mx->dp; /* pointer to the alpha DP matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  /* grow the matrix based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo(cm, mx, errbuf, cp9b, W)) != eslOK) return status;

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
  /* esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE); */

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
	  ESL_DASSERT1((j >= i0 && j <= j0));
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
	    ESL_DASSERT1((j >= i0 && j <= j0));
	    jp_v  = j - jmin[v];
	    for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	      alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
	      yshad[jp_v][dp_v] = USED_EL; 
	    }
	  }
	}
	else { /* cm->endsc[v] == IMPOSSIBLE */
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    ESL_DASSERT1((j >= i0 && j <= j0));
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
       * in alpha[v] are independent of each other, only depending on alpha[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	/* EPN, Mon Nov 26 14:48:07 2007
	 * jn = ESL_MAX(jmin[v], ESL_MIN(jmin[y] + sdr, jmax[y]));
	 * jn = ESL_MAX(jn, jmin[y] + sdr);
	 * jx = ESL_MIN(jmax[v], ESL_MAX(jmax[y] + sdr, jmin[y]));
	 * jx = ESL_MIN(jx, jmax[y] + sdr);
	 */
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);

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
    
    if((cm->flags & CMH_LOCAL_BEGIN) && allow_begin) {
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
       * in the shadow matrix telling fast_alignT() to use the b we return.
       */
      if (v == 0) { 
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp   = W - hdmin[v][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    if (bsc > alpha[0][jp_v][Wp]) {
	      alpha[0][jp_v][Wp] = bsc;
	      yshad[jp_v][Wp] = USED_LOCAL_BEGIN;
	    }
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp; fp = fopen("cyk.mx", "w"); cm_hb_mx_Dump(fp, mx); fclose(fp);*/
  
  Wp = W - hdmin[vroot][j0-jmin[vroot]];
  sc =     alpha[vroot][j0-jmin[vroot]][Wp];

  if (ret_b != NULL)      *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL)    *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */
  if (ret_shadow != NULL) *ret_shadow = shadow;
  if(ret_sc != NULL)      *ret_sc = sc;
  else free_vjd_shadow_matrix(shadow, cm, i0, j0);
  free(el_scA);
  free(yvalidA);

  ESL_DPRINTF1(("fast_cyk_align_hb return sc: %f\n", sc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* 
 * Function: fast_cyk_align()
 * Date:     EPN, Sun Nov 18 19:37:39 2007
 *           
 * Note:     Very similar to inside(), but slightly more efficient.
 *           Identical to fast_cyk_align_hb() but HMM bands are NOT
 *           used.a
 * 
 * Purpose:  Run the inside phase of a CYK alignment algorithm, on a 
 *           subsequence from i0..j0, using a subtree of a model
 *           anchored at a start state vroot, and ending at an end
 *           state vend. (It is a feature of the model layout in
 *           a CM structure that all subtrees are contiguous in the
 *           model.)
 *
 *           We deal with local begins by keeping track of the optimal
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
 *           will be USE_LOCAL_BEGIN, telling fast_alignT() to check b and
 *           start with a local 0->b entry transition. When inside()
 *           is called on smaller subproblems (v != 0 || i0 > 1 || j0
 *           < L), we're using inside() as an engine in divide &
 *           conquer, and we don't use the overall return score nor
 *           shadow matrices, but we do need allow_begin, b, and bsc for
 *           divide&conquer to sort out where a local begin might be used.
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
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
 *           ret_mx    - the dp matrix, we'll allocate it here, NULL if not wanted
 *           ret_sc    - score of optimal, CYK parsetree 
 *                       
 * Returns: <ret_sc>, <ret_b>, <ret_bsc>, <ret_mx>, <ret_shadow>, see 'Args'.
 * 
 * Throws:  <eslOK> on success.
 *          <eslERANGE> if required DP matrix size exceeds CM_NB_MBLIMIT set in structs.h, in 
 *                      this case, alignment has been aborted, ret_* variables are not valid
 */
int
fast_cyk_align(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, void ****ret_shadow,  
	       int allow_begin, int *ret_b, float *ret_bsc, float ****ret_mx, float *ret_sc)
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
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      jp;          /* offset j, j = i0-1+jp */
  int      j_sdr;       /* j - sdr */
  int      d_sd;        /* d - sd */
  float    tsc;         /* a transition score */
  float ***alpha;       /* the DP matrix, we allocate here */
  float    Mb_for_alpha;/* megabytes needed for alpha matrix */

  /* Contract check */
  if(dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_cyk_inside_align(), dsq is NULL.\n");


  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */

  /* allocate alpha (if it's small enough), allocate all decks, no deck reuse */
  Mb_for_alpha = ((float) size_vjd_deck(W, 1, W) * ((float) (cm->M)));
  if(Mb_for_alpha > CM_NB_MX_MB_LIMIT) 
    ESL_FAIL(eslERANGE, errbuf, "fast_cyk_align(), requested size of non-banded DP matrix %.2f Mb > %.2f Mb limit (CM_NB_MX_MB_LIMIT from structs.h).", Mb_for_alpha, (float) CM_NB_MX_MB_LIMIT);
  ESL_DPRINTF1(("Size of alpha matrix: %.2f\n", Mb_for_alpha));

  ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
  for (v = 0; v <= cm->M; v++) alpha[v] = alloc_vjd_deck(L, i0, j0);

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

  /* Main recursion */
  for (v = vend; v >= vroot; v--) {
    float const *esc_v = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    /* Get a shadow deck to fill in and initialize all valid cells for state v */
    if (cm->sttype[v] != E_st) {
      if (cm->sttype[v] == B_st) {
	kshad     = alloc_vjd_kshadow_deck(L, i0, j0); 
	shadow[v] = (void **) kshad;
	/* initialize all valid cells for state v to IMPOSSIBLE (local ends are impossible for B states) */
	assert(! (NOT_IMPOSSIBLE(cm->endsc[v])));
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  for (d = 0; d <= jp; d++) {
	    alpha[v][j][d] = IMPOSSIBLE;
	    kshad[j][d] = USED_EL; 
	  }
	}
      } else { /* ! B_st && ! E_st */
	yshad     = alloc_vjd_yshadow_deck(L, i0, j0);
	shadow[v] = (void **) yshad;
	/* initialize all valid cells for state v */
	if(NOT_IMPOSSIBLE(cm->endsc[v])) {
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d < sd && d <= jp; d++) { 
	      alpha[v][j][d] = IMPOSSIBLE;
	      yshad[j][d] = USED_EL; 
	    }
	    for (d = sd; d <= jp; d++) {
	      alpha[v][j][d] = el_scA[d-sd] + cm->endsc[v];
	      yshad[j][d] = USED_EL; 
	    }
	  }
	}
	else { /* cm->endsc[v] == IMPOSSIBLE */
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++) {
	      alpha[v][j][d] = IMPOSSIBLE;
	      yshad[j][d] = USED_EL; 
	    }
	  }
	}
      }
    }
    
    if(cm->sttype[v] == E_st) { 
      for (jp = 0; jp <= W; jp++) {
	j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
	alpha[v][j][0] = 0.;
	for (d = 1; d <= jp; d++) alpha[v][j][d] = IMPOSSIBLE;
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][j][d] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i    = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = alpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc; 
	      yshad[j][d]    = yoffset;
	    }
	  }
	  alpha[v][j][d] += esc_v[dsq[i--]];
	  alpha[v][j][d]  = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == IR_st) { 
      /* update alpha[v][j][d] cells, for IR states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = alpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc; 
	      yshad[j][d]    = yoffset;
	    }
	  }
	  alpha[v][j][d] += esc_v[dsq[j]];
	  alpha[v][j][d]  = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in alpha[v] are independent of each other, only depending on alpha[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];

	for (jp = sdr; jp <= W; jp++) {
	  j = i0-1+jp;
	  j_sdr = j - sdr;

	  for (d = sd; d <= jp; d++) {
	    if((sc = alpha[y][j_sdr][d - sd] + tsc) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc;
	      yshad[j][d]    = yoffset;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  i = j - 1;
	  for (d = sd; d <= jp; d++) 
	    alpha[v][j][d] += esc_v[dsq[j-d+1]];
	}
	break;
      case MR_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  for (d = sd; d <= jp; d++)
	    alpha[v][j][d] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  i = j - 1;
	  for (d = sd; d <= jp; d++)
	    alpha[v][j][d] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++)
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      for (jp = 0; jp <= W; jp++) { 
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  for (k = 0; k <= d; k++) {
	    if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) { 
	      alpha[v][j][d] = sc;
	      kshad[j][d]    = k;
	    }
	  }
	}
      }
    }
				/* finished calculating deck v. */
      
    if (allow_begin && alpha[v][j0][W] + cm->beginsc[v] > bsc) {
      b   = v;
      bsc = alpha[v][j0][W] + cm->beginsc[v];
    }
    /* Check for whether we need to store an optimal local begin score
     * as the optimal overall score, and if we need to put a flag
     * in the shadow matrix telling fast_alignT() to use the b we return.
     */
    if (allow_begin && v == 0 && bsc > alpha[0][j0][W]) {
      alpha[0][j0][W] = bsc;
      yshad[j0][W] = USED_LOCAL_BEGIN;
    }
  } /* end loop over all v */
  
  sc =     alpha[vroot][j0][W];
  if (ret_b != NULL)      *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL)    *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */
  if (ret_shadow != NULL) *ret_shadow = shadow;
  if (ret_sc     != NULL) *ret_sc = sc;
  else free_vjd_shadow_matrix(shadow, cm, i0, j0);
  if (ret_mx     != NULL) *ret_mx = alpha;
  else free_vjd_matrix(alpha, cm->M, 1, L);

  free(el_scA);

  ESL_DPRINTF1(("fast_cyk_align return sc: %f\n", sc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: fast_alignT_hb()
 * Date:     EPN 03.29.06
 * 
 * Note:     based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *
 * Purpose:  Call either fast_cyk_align_hb() (if !<do_optacc>), 
 *           or optimal_accuracy_align_hb()  (if  <do_optacc>),
 *           get vjd shadow matrix; then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 *        
 *           If (<do_optacc>) then post_mx must != NULL.
 *
 * Returns:  <ret_sc>: if(!do_optacc): score of appended parsetree
 *                     if( do_optacc): avg posterior probability of all i0..j0 residues
 *                                     in optimally accurate alignment in tr 
 *           
 * Throws:  <eslOK>     on success
 *          <eslERANGE> if required CM_HB_MX exceeds CM_HB_MBLIMIT set in structs.h, in 
 *                      this case, alignment has been aborted, tr has not been appended to
 *
 */
int
fast_alignT_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
	       int r, int z, int i0, int j0, 
	       int allow_begin, CM_HB_MX *mx, int do_optacc, CM_HB_MX *post_mx, float *ret_sc)
{
  int status;
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
  /* contract check */
  if(dsq == NULL)      ESL_FAIL(eslEINCOMPAT, errbuf, "fast_alignT_hb(), dsq == NULL.\n");
  if(cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_alignT_hb(), cm->cp9b == NULL.\n");
  if(do_optacc && post_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_alignT_hb(), do_optacc == TRUE but post_mx == NULL.\n");
			 
  /* pointers to cp9b data for convenience */
  CP9Bands_t *cp9b= cm->cp9b;
  int       *jmin = cp9b->jmin;
  int     **hdmin = cp9b->hdmin;
#if eslDEBUGLEVEL >= 1	
  int     **hdmax = cp9b->hdmax;
#endif

  if(do_optacc) {
    status = optimal_accuracy_align_hb(cm, errbuf, dsq, L, i0, j0, 
				       &shadow,	     /* return a shadow matrix to me. */
				       &b, &bsc,     /* if allow_begin is TRUE, gives info on optimal b */
				       mx,           /* the HMM banded mx to fill-in */
				       post_mx,      /* pre-calc'ed posterior matrix */
				       &sc);         /* avg post prob of all emissions in optimally accurate parsetree */
  }
  else {
    status = fast_cyk_align_hb(cm, errbuf, dsq, L, r, z, i0, j0, 
			       &shadow,	     /* return a shadow matrix to me. */
			       allow_begin,  /* TRUE to allow local begins */
			       &b, &bsc,     /* if allow_begin is TRUE, gives info on optimal b */
			       mx,           /* the HMM banded mx */
			       &sc);         /* score of CYK parsetree */
  }
  if(status != eslOK) return status;

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
      /*printf("     mx[v:%4d][jp_v:%4d][dp_v:%4d]: %10.5f\n", v, jp_v, dp_v, mx->dp[v][jp_v][dp_v]);
	if(post_mx != NULL) printf("post_mx[v:%4d][jp_v:%4d][dp_v:%4d]: %10.5f prob: %.5f\n", v, jp_v, dp_v, post_mx->dp[v][jp_v][dp_v], FScore2Prob(post_mx->dp[v][jp_v][dp_v], 1.));*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    ESL_FAIL(eslEINCONCEIVABLE, errbuf, "'Inconceivable!'\n'You keep using that word...'");
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
  if(ret_sc != NULL) *ret_sc = sc;
  return eslOK;
}


/* Function: fast_alignT()
 * Date:     EPN, Sun Nov 18 19:21:30 2007
 * 
 * Note:     based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *
 * Purpose:  Call either fast_cyk_align() (if !<do_optacc>), 
 *           or optimal_accuracy_align()  (if  <do_optacc>),
 *           get vjd shadow matrix; then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 *        
 *           If (<do_optacc>) then post_mx must != NULL.
 *
 *           Very similar to cm_dpsmall.c:insideT() in case of 
 *           CYK alignment, but uses more efficient implementation
 *           of CYK alignment (fast_cyk_align()) as opposed to
 *           inside(). 
 *
 * Returns:  <ret_sc>: if(!do_optacc): score of appended parsetree
 *                     if( do_optacc): avg posterior probability of all i0..j0 residues
 *                                     in optimally accurate alignment in tr 
 *           <ret_mx>: the DP matrix filled, NULL if not wanted 
 * 
 * Throws:  <eslOK>     on success
 *          <eslERANGE> if required DP matrix size exceeds CM_NB_MBLIMIT set in structs.h, in 
 *                      this case, alignment has been aborted, ret_* variables are not valid
 */
int
fast_alignT(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
	    int r, int z, int i0, int j0, 
	    int allow_begin, float ****ret_mx, int do_optacc, float ***post_mx, float *ret_sc)
{
  int       status;
  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;               /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			/* subseq len for bifurcs */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  int       b;                  /* local begin state */
  float     bsc;                /* local begin score */

  /* contract check */
  if(dsq == NULL)      ESL_FAIL(eslEINCOMPAT, errbuf, "fast_alignT(), dsq == NULL.\n");
  if(cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_alignT(), cm->cp9b == NULL.\n");
  if(do_optacc && post_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_alignT(), do_optacc == TRUE but post_mx == NULL.\n");
			 
  if(do_optacc) {
    status = optimal_accuracy_align(cm, errbuf, dsq, L, i0, j0, 
				    &shadow,	     /* return a shadow matrix to me. */
				    &b, &bsc,	     /* if allow_begin is TRUE, gives info on optimal b */
				    ret_mx,           /* the DP mx to fill-in */
				    post_mx,          /* pre-calc'ed posterior matrix */
				    &sc);             /* avg post prob of all emissions in optimally accurate parsetree */
  }
  else {
    status = fast_cyk_align(cm, errbuf, dsq, L, r, z, i0, j0, 
			    &shadow,	    /* return a shadow matrix to me. */
			    allow_begin,    /* TRUE to allow local begins */
			    &b, &bsc,	    /* if allow_begin is TRUE, gives info on optimal b */
			    ret_mx,         /* the DP mx to fill in */
			    &sc);           /* score of CYK parsetree */
  }
  if(status != eslOK) return status;

  pda = esl_stack_ICreate();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  while (1) {
    if (cm->sttype[v] == B_st) {
      k = ((int **) shadow[v])[j][d];   /* k = len of right fragment */

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
    } else {
      yoffset = ((char **) shadow[v])[j][d];

      /*printf("v : %d | r : %d | z : %d | i0 : %d | \n", v, r, z, i0);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
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
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  if(ret_sc != NULL) *ret_sc = sc;
  return eslOK;
}

/* Function: FastAlignHB()
 * Incept:   EPN, Fri Oct 26 09:31:43 2007
 * 
 * Note:     based on CYKInside_b_jd() [11.04.05] which was based on CYKInside_b() 
 *           which was based on CYKInside() [SRE, Sun Jun  3 19:48:33 2001 [St. Louis]]
 *
 * Purpose:  Wrapper for the fast_alignT_hb() routine - solve a full
 *           alignment problem using CYK or optimal accuracy and
 *           return the traceback and the score, without dividing &
 *           conquering, but by using bands on the j and d dimensions
 *           of the DP matrix.  Bands derived by HMM Forward/Backward
 *           runs. Optionally return a postal code.
 *           
 *           Identical to FastAlign() but HMM bands are used here.
 * 
 *           Input arguments allow this function to be run in 4 'modes':
 *
 *           mode      returns                 arguments
 *           ----  ----------------  -----------------------------------
 *                 tr        pcodes  do_optacc  post_mx   ret_pcode{1,2}
 *                 ----------------  -----------------------------------
 *              1. CYK       no      FALSE       NULL      NULL
 *              2. CYK       yes     FALSE      !NULL     !NULL
 *              3. Opt acc   no      TRUE       !NULL      NULL
 *              4. Opt acc   yes     TRUE       !NULL     !NULL
 *
 *           CYK parsetrees are most likely parsetree, 'Opt acc' parsetrees
 *           are Holmes/Durbin optimally accurate parsetrees, the parse the 
 *           maximizes the expected accuracy of all aligned residues.
 *
 *           Note: if ret_tr is NULL, parsetree is not returned.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence, 1..L
 *           L         - length of sequence 
 *           i0        - start of target subsequence (often 1, beginning of dsq)
 *           j0        - end of target subsequence (often L, end of dsq)
 *           mx        - the main dp matrix, only cells within bands in cm->cp9b will be valid. 
 *           do_optacc - TRUE to not do CYK alignment, determine the Holmes/Durbin optimally 
 *                       accurate parsetree in ret_tr, requires post_mx != NULL
 *           post_mx   - dp matrix for posterior calculation, can be NULL only if !do_optacc
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_pcode1- RETURN: postal code 1, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_pcode2- RETURN: postal code 2, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_sc    - if(!do_optacc): score of the alignment in bits.
 *                       if( do_optacc): average posterior probability of all L aligned residues 
 *                       in optimally accurate alignment
 * 
 * Returns: <ret_tr>, <ret_pcode1>, <ret_pcode2>, <ret_sc>, see 'Args' section
 * 
 * Throws:  <eslOK> on success; 
 *          <eslERANGE> if required CM_HB_MX for FastInsideAlignHB(), FastOutsideAlignHB() or
 *                      fast_cyk_align_hb() exceeds CM_HB_MBLIMIT set in structs.h, in this 
 *                      case, alignment has been aborted, ret_* variables are not valid 
 */

int
FastAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, CM_HB_MX *mx, int do_optacc,
	    CM_HB_MX *post_mx, Parsetree_t **ret_tr, char **ret_pcode1, char **ret_pcode2, float *ret_sc)
{
  int          status;
  Parsetree_t *tr;
  float        sc;
  int          do_post;
  char        *pcode1;
  char        *pcode2;
  int          have_pcodes;
  have_pcodes = (ret_pcode1 != NULL && ret_pcode2 != NULL) ? TRUE : FALSE;

  /* Contract check */
  if(dsq == NULL)                    ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlignHB(), dsq is NULL.\n");
  if(mx == NULL)                     ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlignHB(), mx is NULL.\n");
  if(post_mx == NULL && have_pcodes) ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlignHB(), post_mx == NULL but ret_pcode{1|2} != NULL.\n");
  if(do_optacc && post_mx == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlignHB(), do_optacc is TRUE, but post_mx == NULL.\n");

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));*/
  do_post = (do_optacc || have_pcodes) ? TRUE : FALSE;

  /* if doing post, fill Inside, Outside, Posterior matrices, in that order */
  if(do_post) { 
    if((status = FastInsideAlignHB (cm, errbuf, dsq, i0, j0, mx, NULL)) != eslOK) return status;
    if((status = FastOutsideAlignHB(cm, errbuf, dsq, i0, j0, post_mx, mx, ((cm->align_opts & CM_ALIGN_CHECKINOUT) && (! cm->flags & CMH_LOCAL_END)), NULL)) != eslOK) return status;
    /* Note: we can only check the posteriors in FastOutsideAlignHB() if local begin/ends are off */
    if((status = CMPosteriorHB(cm, errbuf, i0, j0, mx, post_mx, post_mx)) != eslOK) return status;   
    if(cm->align_opts & CM_ALIGN_CHECKINOUT) { 
      if((status = CMCheckPosteriorHB(cm, errbuf, i0, j0, post_mx)) != eslOK) return status;
      printf("\nHMM banded posteriors checked.\n\n");
    }
  }

  /* Create the parse tree, and initialize. */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */
  /* Fill in the parsetree (either CYK or optimally accurate (if do_optacc)), 
   * this will overwrite mx if (do_post) caused it to be filled in FastInsideAlignHB 
   */
  if((status = fast_alignT_hb(cm, errbuf, dsq, L, tr, 0, cm->M-1, i0, j0, TRUE, mx, do_optacc, post_mx, &sc)) != eslOK) return status;

  if(have_pcodes) {
    CMPostalCodeHB(cm, L, post_mx, tr, &pcode1, &pcode2);
    *ret_pcode1 = pcode1;
    *ret_pcode2 = pcode2;
  }

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  if (ret_sc != NULL) *ret_sc = sc;
  ESL_DPRINTF1(("returning from FastAlignHB() sc : %f\n", sc)); 
  return eslOK;
}


/* Function: FastAlign()
 * Date:     EPN, Sun Nov 18 19:26:45 2007
 *
 * Note:     Very similar to cm_dpsmall.c:CYKInside() for case
 *           of CYK alignment, but uses slightly more efficient
 *           implementation (fast_cyk_align() instead of inside()).
 *
 * Purpose:  Wrapper for the fast_alignT() routine - solve a full
 *           alignment problem either by CYK or using optimal
 *           accuracy, return the traceback and the score, without
 *           dividing & conquering. Optionally return a postal code.
 *           
 *           Identical to FastAlignHB() but HMM bands are NOT used here.
 * 
 *           Input arguments allow this function to be run in 4 'modes':
 *
 *           mode      returns                 arguments
 *           ----  ----------------  -----------------------------------
 *                 tr        pcodes   do_optacc  post_mx   ret_pcode{1,2}
 *                 ----------------  -----------------------------------
 *              1. CYK       no      FALSE       NULL      NULL
 *              2. CYK       yes     FALSE      !NULL     !NULL
 *              3. Opt acc   no      TRUE       !NULL      NULL
 *              4. Opt acc   yes     TRUE       !NULL     !NULL
 *
 *           CYK parsetrees are most likely parsetree, 'Opt acc' parsetrees
 *           are Holmes/Durbin optimally accurate parsetrees, the parse the 
 *           maximizes the expected accuracy of all aligned residues.
 *
 *           QDB: if dmin and dmax are non-NULL, the CYK alignment conditional
 *           on the bands in dmin dmax can be calculated, but only mode 1
 *           can be run with QDB. The reason is that there is no implementation
 *           of FInside() and FOutside() that use QDBs yet. This is not a big
 *           deal b/c it's better to use HMM bands than query-dependent bands 
 *           (QDB) for alignment anyway.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence, 1..L
 *           L         - length of sequence 
 *           i0        - start of target subsequence (often 1, beginning of dsq)
 *           j0        - end of target subsequence (often L, end of dsq)
 *           ret_mx    - RETURN: the main dp matrix, we'll allocate and fill it, must be non-NULL
 *           do_optacc - TRUE to not do CYK alignment, determine the Holmes/Durbin optimally 
 *                       accurate parsetree in ret_tr, requires post_mx != NULL
 *           ret_post_mx- dp matrix for posterior calculation, we'll allocate and fill it if nec,
 *                        can be NULL only if !do_optacc
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_pcode1- RETURN: postal code 1, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_pcode2- RETURN: postal code 2, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_sc    - if(!do_optacc): score of the alignment in bits.
 *                       if( do_optacc): average posterior probability of all L aligned residues 
 *                       in optimally accurate alignment
 * 
 * Returns: <ret_tr>, <ret_pcode1>, <ret_pcode2>, <ret_sc>, see 'Args' section
 * 
 * Throws:  <eslOK> on success; 
 */
int
FastAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, float ****ret_mx, int do_optacc,
	  float ****ret_post_mx, Parsetree_t **ret_tr, char **ret_pcode1, char **ret_pcode2, float *ret_sc)
{
  int          status;
  Parsetree_t *tr;
  float        sc;
  int          do_post;
  char        *pcode1;
  char        *pcode2;
  int          have_pcodes;
  have_pcodes = (ret_pcode1 != NULL && ret_pcode2 != NULL) ? TRUE : FALSE;

  /* Contract check */
  if(ret_mx == NULL)                     ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlign(), ret_mx == NULL.\n");
  if(dsq == NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlign(), dsq is NULL.\n");
  if(ret_post_mx == NULL && have_pcodes) ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlign(), post_mx == NULL but ret_pcode{1|2} != NULL.\n");
  if(do_optacc && ret_post_mx == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlign(), do_optacc is TRUE, but post_mx == NULL.\n");
  if(do_optacc && ret_mx == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "FastAlign(), do_optacc is TRUE, but ret_mx == NULL.\n");

  do_post = (do_optacc || have_pcodes) ? TRUE : FALSE;

  /* if doing post, fill Inside, Outside, Posterior matrices, in that order */
  if(do_post) { 
    if((status = FastInsideAlign (cm, errbuf, dsq, i0, j0, ret_mx,  NULL)) != eslOK) return status;
    if((status = FastOutsideAlign(cm, errbuf, dsq, i0, j0, ret_post_mx, *ret_mx, ((cm->align_opts & CM_ALIGN_CHECKINOUT) && (! cm->flags & CMH_LOCAL_END)), NULL)) != eslOK) return status;
    /* Note: we can only check the posteriors in FastOutsideAlign() if local begin/ends are off */
    if((status = CMPosterior(cm, errbuf, i0, j0, *ret_mx, *ret_post_mx, *ret_post_mx)) != eslOK) return status;   
    if(cm->align_opts & CM_ALIGN_CHECKINOUT) { 
      if((status = CMCheckPosterior(cm, errbuf, i0, j0, *ret_post_mx)) != eslOK) return status;
      printf("\nPosteriors checked.\n\n");
    }
    /* we have to free ret_mx, so we can check it's size, reallocate and refill it in fast_alignT()
     * this is wasteful, but if we were being efficient we'd be using HMM bands anyway... */
    free_vjd_matrix(*ret_mx, cm->M, 1, L);
  }
  /* Create the parse tree, and initialize. */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */
  /* Fill in the parsetree (either CYK or optimally accurate (if do_optacc)) */
  if(do_optacc) { /* we have to send the filled *ret_post_mx */
    if((status = fast_alignT(cm, errbuf, dsq, L, tr, 0, cm->M-1, i0, j0, TRUE, ret_mx, do_optacc, *ret_post_mx, &sc)) != eslOK) return status;
  }
  else { /* don't need to send *ret_post_mx (in fact, it could be NULL) */
    if((status = fast_alignT(cm, errbuf, dsq, L, tr, 0, cm->M-1, i0, j0, TRUE, ret_mx, do_optacc, NULL, &sc)) != eslOK) return status;
  }

  if(have_pcodes) {
    CMPostalCode(cm, L, *ret_post_mx, tr, &pcode1, &pcode2);
    *ret_pcode1 = pcode1;
    *ret_pcode2 = pcode2;
  }

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  if (ret_sc != NULL) *ret_sc = sc;
  ESL_DPRINTF1(("returning from FastAlign() sc : %f\n", sc)); 
  return eslOK;
}

/*
 * Function: FastInsideAlignHB()
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
 *           This function complements FastOutsideAlignHB().
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score
 * 
 * Returns:  <ret_sc>
 *
 * Throws:  <eslOK> on success
 *          <eslERANGE> if required CM_HB_MX for exceeds CM_HB_MBLIMIT set in structs.h, 
 *                      in this case, alignment has been aborted, ret_sc is not valid
 */
int
FastInsideAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, CM_HB_MX *mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,jp,d,i,k;	/* indices in sequence dimensions */
  float    fsc;		/* the final score */
  float    tsc;         /* a temporary variable holding a transition score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  float    bsc;		/* summed score for using all local begins */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
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
  int      have_el;            /* TRUE if local ends are on */

  /* Contract check */
  if(dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "FastInsideAlignHB(), dsq is NULL.\n");
  if (mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "FastInsideAlignHB(), mx is NULL.\n");
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "FastInsideAlignHB(), cm->cp9b is NULL.\n");

  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrix */
  float ***alpha = mx->dp;     /* pointer to the alpha DP matrix */

  /* Allocations and initializations */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  /* grow the matrix based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo(cm, mx, errbuf, cp9b, W)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;
  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands)
   */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* if local ends are on, replace the EL deck IMPOSSIBLEs with EL scores */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el) { 
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0;  d <= jp; d++) alpha[cm->M][j][d] = el_scA[d];
    }
  }

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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
       * in alpha[v] are independent of each other, only depending on alpha[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
      
    if(cm->flags & CMH_LOCAL_BEGIN && NOT_IMPOSSIBLE(cm->beginsc[v])) { /* if local begins are on */
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
    }
    /* include the bsc as part of alpha[0][jp_v][Wp] */
    if ((cm->flags & CMH_LOCAL_BEGIN) && v == 0) { 
      if(j0 >= jmin[0] && j0 <= jmax[0]) {
	jp_v = j0 - jmin[v];
	Wp = W - hdmin[0][jp_v];
	if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	  alpha[0][jp_v][Wp] = FLogsum(alpha[0][jp_v][Wp], bsc);
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp; fp = fopen("ins.mx", "w"); cm_ihb_mx_Dump(fp, mx); fclose(fp);*/

  Wp  = W - hdmin[0][j0-jmin[0]];
  fsc = alpha[0][j0-jmin[0]][Wp];

  free(el_scA);
  free(yvalidA);

  if(ret_sc != NULL) *ret_sc = fsc;
  ESL_DPRINTF1(("FastInsideAlignHB() return sc: %f\n", fsc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/*
 * Function: FastInsideAlign()
 * Date:     EPN, Mon Nov 19 06:21:51 2007
 *
 * Purpose:  Run the inside algorithm on a target sequence 
 *           without using bands. 
 *
 *           Identical to FastInsideAlignHB() but no bands
 *           are used.
 * 
 *           Very similar with fast_cyk_inside_align(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *
 *           This function complements FastOutsideAlign().
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           ret_mx    - RETURN: the dp matrix, we'll allocate and fill it here 
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score
 *                       
 * Returns:  <ret_sc>, <ret_mx>
 *
 * Throws:   <eslOK> on success.
 *           <eslERANGE> if required size of DP matrix exceeds CM_NB_MBLIMIT set in structs.h, 
 *                       in this case, alignment has been aborted, ret_sc, ret_mx are not valid
 */
int
FastInsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float ****ret_mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* the final score */
  float    tsc;         /* a temporary variable holding a transition score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  float    bsc;		/* summed score for using all local begins */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp;          /* offset j, j = i0-1+jp */
  int      j_sdr;       /* j - sdr */
  int      d_sd;        /* d - sd */
  float ***alpha;       /* the DP matrix, we allocate here */
  float    Mb_for_alpha;/* megabytes needed for alpha matrix */
  int      have_el;     /* TRUE if local ends are on */

  /* Contract check */
  if(dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "FastInsideAlign(), dsq is NULL.\n");

  /* Allocations and initializations */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  
  /* allocate alpha (if it's small enough), allocate all decks, no deck reuse */
  Mb_for_alpha = ((float) size_vjd_deck(W, 1, W) * ((float) (cm->M+1)));
  if(Mb_for_alpha > CM_NB_MX_MB_LIMIT) 
    ESL_FAIL(eslERANGE, errbuf, "FastInsideAlign(), requested size of non-banded DP matrix %.2f Mb > %.2f Mb limit (CM_NB_MX_MB_LIMIT from structs.h).", Mb_for_alpha, (float) CM_NB_MX_MB_LIMIT);
  ESL_DPRINTF1(("Size of alpha matrix: %.2f\n", Mb_for_alpha));

  ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
  for (v = 0; v <= cm->M; v++) alpha[v] = alloc_vjd_deck(W, 1, W);

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* Fill the EL cm->M deck, with EL scores if local ends are on, IMPOSSIBLE if not (they'll never be used anyway) */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el) { 
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0;  d <= jp; d++) alpha[cm->M][j][d] = el_scA[d];
    }
  }
  else { /* local ends are off */
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0;  d <= jp; d++) alpha[cm->M][j][d] = IMPOSSIBLE;
    }
  }
  
  /* Main recursion  */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; 
    float const *tsc_v = cm->tsc[v];
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    
    /* initialize all valid cells for state v to the local end prob, if they're allowed,
     * and IMPOSSIBLE if not  */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0;  d < sd && d <= jp; d++) alpha[v][j][d] = IMPOSSIBLE;
	for (d = sd;           d <= jp; d++) alpha[v][j][d]   = el_scA[d-sd] + cm->endsc[v];
      }
    }
    else {
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) 
	  alpha[v][j][d] = IMPOSSIBLE;
      }
    }
    /* E_st: easy, no children, and d must be 0 for all valid j */
    if(cm->sttype[v] == E_st) { 
      for (jp = 0; jp <= W; jp++) {
	j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
	alpha[v][j][0] = 0.;
	for (d = 1; d <= jp; d++) alpha[v][j][d] = IMPOSSIBLE;
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i    = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], alpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	  }
	  alpha[v][j][d] += esc_v[dsq[i--]];
	  alpha[v][j][d]  = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == IR_st) { 
      /* update alpha[v][jp_v][dp_v] cells, for IR states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i    = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], alpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	  }
	  alpha[v][j][d] += esc_v[dsq[j]];
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in alpha[v] are independent of each other, only depending on alpha[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];

	for (jp = sdr; jp <= W; jp++) {
	  j = i0-1+jp;
	  j_sdr = j - sdr;

	  for (d = sd; d <= jp; d++) {
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], (alpha[y][j_sdr][d-sd] + tsc));;
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  i = j - 1;
	  for (d = sd; d <= jp; d++) 
	    alpha[v][j][d] += esc_v[dsq[j-d+1]];
	}
	break;
      case MR_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  for (d = sd; d <= jp; d++)
	    alpha[v][j][d] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  i = j - 1;
	  for (d = sd; d <= jp; d++)
	    alpha[v][j][d] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++)
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      for (jp = 0; jp <= W; jp++) { 
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  for (k = 0; k <= d; k++) {
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], alpha[y][j-k][d-k] + alpha[z][j][k]); 
	  }
	}
      }
    }				/* finished calculating deck v. */
      
    if (cm->flags & CMH_LOCAL_BEGIN && NOT_IMPOSSIBLE(cm->beginsc[v])) { 
      /* add in score for local begin getting us to the root. */
      bsc = FLogsum(bsc, alpha[v][j0][W] + cm->beginsc[v]);
    }
    /* include the bsc as part of alpha[0][jp_v][Wp] */
    if(v == 0)  
      alpha[0][j0][W] = FLogsum(alpha[v][j0][W], bsc);
  } /* end loop over all v */

  sc =     alpha[0][j0][W];
  free(el_scA);
  if(ret_sc != NULL) *ret_sc = sc;
  if (ret_mx     != NULL) *ret_mx = alpha;
  else free_vjd_matrix(alpha, cm->M, 1, W);

  ESL_DPRINTF1(("FastInsideAlign() return sc: %f\n", sc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/*
 * Function: FastOutsideAlignHB()
 * Date:     EPN, Thu Nov  8 18:40:05 2007
 *
 * Purpose:  Run the outside algorithm using bands
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 *
 *           A CM_FHB_MX DP matrix must be passed in. Only
 *           cells valid within the bands given in the CP9Bands_t <cm->cp9b>
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
 *           If <do_check> is TRUE (and the CM is not in local mode) 
 *           we check that the outside calculations are consistent 
 *           with the inside calculations (in ins_mx). 
 *           This check is described in comments towards the end of 
 *           the function. 
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ins_mx    - the dp matrix from the Inside run calculation (required)
 *           do_check  - TRUE to attempt to check 
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *                       ends are on (see *** comment towards end of function).
 *
 * Returns:  <ret_sc>
 *
 * Throws:  <eslOK> on success
 *          <eslERANGE> if required CM_HB_MX for exceeds CM_HB_MBLIMIT set in structs.h, 
 *                      in this case, alignment has been aborted, ret_sc is not valid *                       
 */
int
FastOutsideAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, CM_HB_MX *mx, 
		    CM_HB_MX *ins_mx, int do_check, float *ret_sc)
{
  int      status;
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
  int      jp_el;              /* offset j in EL deck, jp_el = j - i0 */
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

  /* Contract check */
  if (dsq == NULL)                                     ESL_FAIL(eslEINCOMPAT, errbuf, "FastOutsideAlignHB(), dsq is NULL.\n");
  if (mx == NULL)                                      ESL_FAIL(eslEINCOMPAT, errbuf, "FastOutsideAlignHB(), mx is NULL.\n");
  if (cm->cp9b == NULL)                                ESL_FAIL(eslEINCOMPAT, errbuf, "FastOutsideAlignHB(), cm->cp9b is NULL.\n");
  if (ins_mx == NULL)                                  ESL_FAIL(eslEINCOMPAT, errbuf, "FastOutsideAlignHB(), ins_mx is NULL.\n");
  if (cm->flags & CMH_LOCAL_END) do_check = FALSE; /* Code for checking doesn't apply in local mode. See below. */

  /* DP matrix variables */
  float ***beta  = mx->dp;     /* pointer to the Oustide DP mx */
  float ***alpha = ins_mx->dp; /* pointer to the Inside DP mx (already calc'ed and passed in) */

  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;              /* the summed prob of all local begins */
  W   = j0-i0+1;		 /* the length of the subsequence -- used in many loops  */
				 /* if caller didn't give us a deck pool, make one */
  esc_vAA = cm->oesc;            /* a ptr to the optimized emission scores */

  /* grow the matrix based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo(cm, mx, errbuf, cp9b, W)) != eslOK) return status;

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
  if (cm->flags & CMH_LOCAL_END) {
    for(jp = 0; jp <= W; jp++) { 
      for (d = 0; d <= jp; d++) beta[cm->M][jp][d] = IMPOSSIBLE;
    }
    /* We don't have to worry about vroot -> EL transitions the way 
     * cm_dpsmall.c::outside() does, because vroot = 0.
     */
  }
  /* If we can do a local begin into v, overwrite IMPOSSIBLE with the local begin score. 
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
	ESL_DASSERT1((j >= i0 && j <= j0));
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
	jn = ESL_MAX(jmin[v], jmin[y]-sdr);
	jx = ESL_MIN(jmax[v], jmax[y]-sdr);
	for (j = jx; j >= jn; j--) {
	  ESL_DASSERT1((j >= i0 && j <= j0));
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

    /* deal with local alignment end transitions v->EL, 
     * these values get their own matrix.
     */
    /* deal with local alignment end transitions v->EL (EL = deck at M.) */
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);   /* note emitmode is for state v */
      
      jn = jmin[v] - sdr;
      jx = jmax[v] - sdr;
      for (j = jn; j <= jx; j++) {
	jp_v =  j - jmin[v];
	jp_el = j-i0+1;     /* offset j in cm->M deck */
	dn   = hdmin[v][jp_v + sdr] - sd;
	dx   = hdmax[v][jp_v + sdr] - sd;
	i    = j-dn+1;                     /* we'll decrement this in for (d... loops inside switch below */
	dp_v = dn - hdmin[v][jp_v + sdr];  /* we'll increment this in for (d... loops inside switch below */

	switch (emitmode) {
	case EMITPAIR:
	  for (d = dn; d <= dx; d++, dp_v++, i--) {
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][jp_el][d] = FLogsum(beta[cm->M][jp_el][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v] 
								    + escore));
	  }
	  break;
	case EMITLEFT:
	  for (d = dn; d <= dx; d++, dp_v++, i--) {
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][jp_el][d] = FLogsum(beta[cm->M][jp_el][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v] 
								    + escore));
	  }
	  break;
	  
	case EMITRIGHT:
	  escore = esc_vAA[v][dsq[j+1]];
	  for (d = dn; d <= dx; d++, dp_v++) {
	    beta[cm->M][jp_el][d] = FLogsum(beta[cm->M][jp_el][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v]
								    + escore));
	  }
	  break;
	  
	case EMITNONE:
	  for (d = dn; d <= dx; d++, dp_v++) {
	    beta[cm->M][jp_el][d] = FLogsum(beta[cm->M][jp_el][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v]));
	  }
	  break;
	}
      }
    }
  } /* end loop over decks v. */
  /*FILE *fp; fp = fopen("out.mx", "w"); cm_ihb_mx_Dump(fp, mx); fclose(fp);*/

  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, EL->EL transitions. (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1] + cm->el_selfsc));
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

  if(fail_flag) ESL_FAIL(eslFAIL, errbuf, "Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END)) ESL_DPRINTF1(("\tFastOutsideAlignHB() sc : %f\n", freturn_sc));
  else                             ESL_DPRINTF1(("\tFastOutsideAlignHB() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  if (ret_sc != NULL) *ret_sc = freturn_sc;
  return eslOK;
}  


/*
 * Function: FastOutsideAlign()
 * Date:     EPN, Mon Nov 19 07:00:37 2007
 *
 * Purpose:  Run the outside algorithm on a target sequence
 *           without using bands.
 *
 *           Very similar to FastInsideAlignHB() but no bands
 *           are used and recursion nesting order for all non - BEGL_S 
 *           and non BEGR_S is: for v { for j { for d { for y {}}}}.
 *           This is slower, but corrects some precision issues I was
 *           having during testing, if you want fast alignment, you should
 *           be using HMM banded alignment anyway.
 *
 *           A float, non-banded DP matrix must be passed in. Only
 *           cells valid within the bands given in the CP9Bands_t <cm->cp9b>
 *           will be valid. 
 *
 *           If <do_check> is TRUE (and the CM is not in local mode) 
 *           we check that the outside calculations are consistent 
 *           with the inside calculations (in ins_mx). 
 *           This check is described in comments towards the end of 
 *           the function. 
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           ret_mx    - the dp matrix, we'll allocate it here, NULL if not wanted
 *           ins_mx    - the pre-filled dp matrix from the Inside run calculation (required)
 *           do_check  - TRUE to attempt to check 
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *                       ends are on (see *** comment towards end of function).
 *
 * Returns:  <ret_sc>, <ret_mx>
 *
 * Throws:   <eslOK> on success
 *           <eslERANGE> if required DP matrix size exceeds CM_NB_MBLIMIT set in structs.h, in 
 *                       this case, alignment has been aborted, ret_* variables are not valid
 */
int 
FastOutsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float ****ret_mx, 
		 float ***ins_mx, int do_check, float *ret_sc)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float    fsc;     	       /* a temporary variable holding a float score */
  float  **esc_vAA;            /* ptr to cm->oesc, optimized emission scores */
  float    escore;	       /* an emission score, tmp variable */
  int      W;		       /* subsequence length */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      jp;		       /* j': relative position in the subsequence  */
  float    bsc;		       /* total score for using local begin states */
  float    freturn_sc;         /* P(S|M)/P(S|R) */
  /* variables used only if do_check */
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */
  int      fail_flag = FALSE;  /* set to TRUE if do_check and we see a problem */
  /* indices used in the depths of the DP recursion */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */
  /* DP matrix */
  float ***beta;        /* the DP matrix, we allocate here */
  float    Mb_for_beta; /* megabytes needed for alpha matrix */

  /* Contract check */
  if (dsq == NULL)                                     ESL_FAIL(eslEINCOMPAT, errbuf, "FastOutsideAlign(), dsq is NULL.\n");
  if (ins_mx == NULL)                                  ESL_FAIL(eslEINCOMPAT, errbuf, "FastOutsideAlign(), ins_mx is NULL.\n");
  if (cm->flags & CMH_LOCAL_END) do_check = FALSE; /* Code for checking doesn't apply in local mode. See below. */

  /* inside DP matrix */
  float ***alpha = ins_mx; /* pointer to the Inside DP mx (already calc'ed and passed in) */

  /* Allocations and initializations */
  bsc = IMPOSSIBLE;              /* the summed prob of all local begins */
  W   = j0-i0+1;		 /* the length of the subsequence -- used in many loops  */
				 /* if caller didn't give us a deck pool, make one */
  esc_vAA = cm->oesc;            /* a ptr to the optimized emission scores */

  /* allocate beta (if it's small enough), allocate all decks, no deck reuse */
  Mb_for_beta = ((float) size_vjd_deck(W, 1, W) * ((float) (cm->M)));
  if(Mb_for_beta > CM_NB_MX_MB_LIMIT) 
    ESL_FAIL(eslERANGE, errbuf, "FastOutsideAlign(), requested size of non-banded DP matrix %.2f Mb > %.2f Mb limit (CM_NB_MX_MB_LIMIT from structs.h).", Mb_for_beta, (float) CM_NB_MX_MB_LIMIT);
  ESL_DPRINTF1(("Size of beta matrix: %.2f\n", Mb_for_beta));

  ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
  for (v = 0; v <= cm->M; v++) beta[v] = alloc_vjd_deck(W, 1, W);

  /* Init whole matrix to IMPOSSIBLE. */
  for (v = 0; v <= cm->M; v++) {
    for(jp = 0; jp <= W; jp++) { 
      j = i0-1+jp;
      for (d = 0; d <= jp; d++) 
	beta[v][j][d] = IMPOSSIBLE;
    }
  }
  /* set beta[0][j0][W] to 0., all parses must end there */
  beta[0][j0][W] = 0.;

  /* init local begin cells for emitting full seq (j==j0 && d == W) */
  if (cm->flags & CMH_LOCAL_BEGIN) { 
    for (v = 1; v < cm->M; v++) 
      beta[v][j0][W] = cm->beginsc[v];
  }
  /* done allocation/initialization */

  /* Recursion: main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    sd  = StateDelta(cm->sttype[v]);
    sdr = StateRightDelta(cm->sttype[v]);

    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      y = cm->plast[v];	/* the parent bifurcation    */
      z = cm->cnum[y];	/* the other (right) S state */
      for(jp = 0; jp <= W; jp++) { 
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  for (k = 0; k <= (W-j); k++) {
	    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+k][d+k] + alpha[z][j+k][k]));
	  }
	}
      }
    } /* end of 'if (cm->stid[v] == BEGL_S */
    else if (cm->stid[v] == BEGR_S) {
      y = cm->plast[v];	  /* the parent bifurcation    */
      z = cm->cfirst[y];  /* the other (left) S state  */
      for(jp = 0; jp <= W; jp++) { 
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  for (k = 0; k <= (j-d); k++) {
	    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j][d+k] + alpha[z][j-d][k]));
	  }
	}
      }
    } /* end of 'else if (cm->stid[v] == BEGR_S */
    else { /* (cm->sttype[v] != BEGL_S && cm->sttype[v] != BEGR_S */ 
      for (jp = W; jp >= 0; jp--) {
	j = i0-1+jp;
	i = j-jp+1;
	for (d = jp; d >= 0; d--, i++) {
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    sd  = StateDelta(cm->sttype[y]);
	    sdr = StateRightDelta(cm->sttype[y]);
	    switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == j0 || d == jp) continue; /* boundary condition */
		escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;

	      case ML_st:
	      case IL_st: 
		if (d == jp) continue;	/* boundary condition (note when j=0, d=0*/
		escore = esc_vAA[y][dsq[i-1]];
		beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		escore = esc_vAA[y][dsq[j+1]];
		beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset]));
		break;
	    } /* end of switch(cm->sttype[y] */  
	  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	  if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;
	} /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
      } /* end loop over jp. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] != BEGL_S, BEGR_S */
    /* we're done calculating deck v for everything but local begins */

    /* deal with local alignment end transitions v->EL (EL = deck at M.) */
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);   /* note emitmode is for state v */
      
      for (jp = 0; jp <= W; jp++) { 
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  i = j-d+1;
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    if (j == j0 || d == jp) continue; /* boundary condition */
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v] 
							    + escore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == jp) continue;	
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v] 
							    + escore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == j0) continue;
	    escore = esc_vAA[v][dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v]
							    + escore));
	    break;
	  case S_st:
	  case D_st:
	  case E_st:
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v]));
	    break;
	  }
	}
      }
    }
  }
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, EL->EL transitions. (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1] + cm->el_selfsc));
    }
  }

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
	for (jp = 0; jp <= W-sdr; jp++) {
	  j = i0-1+jp;
	  for (d = 0; d <= jp-sd; d++) {
	    fsc = FLogsum(fsc, (alpha[v][j][d] + beta[v][j][d]));
	  }
	}
      }
      /*printf("checking node: %d | sc: %.6f\n", n, fsc);*/
      diff = fsc - alpha[0][j0][W];
      if(diff > 0.01 || diff < -0.01) { /* we're using floats, this should be pretty precise */
	fail_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, fsc, alpha[0][j0][W], diff);
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
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
      freturn_sc = FLogsum(freturn_sc, (beta[v][j][0]));
    }
  }
  else { /* return_sc = P(S|M) / P(S|R) from Inside() */
    freturn_sc = alpha[0][j0][W];
  }

  if(fail_flag) ESL_FAIL(eslFAIL, errbuf, "Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END)) ESL_DPRINTF1(("\tFastOutsideAlign() sc : %f\n", freturn_sc));
  else ESL_DPRINTF1(("\tFastOutsideAlign() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));
  if(ret_sc != NULL) *ret_sc = freturn_sc;
  if (ret_mx     != NULL) *ret_mx = beta;
  else free_vjd_matrix(beta, cm->M, 1, W);

  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}  


/*
 * Function: CMPosteriorHB()
 * Date:     EPN 05.27.06 
 * Note:     based on Ian Holmes' P7EmitterPosterior() from HMMER's 2.x postprob.c
 *
 * Purpose:  Combines HMM banded Inside and Outside matrices into a 
 *           posterior probability matrix. Any cells outside of
 *           HMM bands do not exist in memory. The value in post[v][jp_v][dp_v] 
 *           is the log of the posterior probability of a parse subtree rooted at v 
 *           emitting the subsequence i..j (i=j-d+1). Where j = jp_v + jmin[v],
 *           and d = dp_v + hdmin[v][jp_v]. The caller must provide a <post> CM_HB_MX
 *           matrix, but this matrix may be the same matrix as that provided as
 *           Outside <out_mx>, (overwriting it will not compromise the algorithm).
 *           
 * Args:     cm       - the model
 *           errbuf   - char buffer for reporting errors
 *           i0       - first position of target seq we're aligning, usually 1 
 *           j0       - final position of target seq we're aligning, usually L (length of seq) 
 *           ins_mx   - pre-calculated Inside matrix 
 *           out_mx   - pre-calculated Outside matrix
 *           post_mx  - pre-allocated matrix for Posteriors 
 *
 * Return:   eslOK on succes, eslEINCOMPAT on contract violation, eslEMEM on memory allocation error
 */
int
CMPosteriorHB(CM_t *cm, char *errbuf, int i0, int j0, CM_HB_MX *ins_mx, CM_HB_MX *out_mx, CM_HB_MX *post_mx)
{
  int      status;
  int      v, j, d;
  float    sc;   /* total score, the log probability of the current seq  */
  int      jp_v; /* j index for state v in alpha/beta with HMM bands */
  int      dp_v; /* d index for state v in alpha/beta with HMM bands */
  int      L;    /* length of sequence */
  int      Lp;   /* offset length */
  int      have_el; /* TRUE if we have local ends */
  L  = j0-i0+1;
  
  /* Contract check */
  if (ins_mx == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), ins_mx is NULL.\n");
  if (out_mx == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), out_mx is NULL.\n");
  if (post_mx == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), post_mx is NULL.\n");
  if (cm->cp9b == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), cm->cp9b is NULL.\n");
  if (ins_mx->L != L)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), ins_mx->L != L passed in.\n");
  if (out_mx->L != L)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), out_mx->L != L passed in.\n");
  if (ins_mx->M != cm->M) ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), ins_mx->M != cm->M.\n");
  if (out_mx->M != cm->M) ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosteriorHB(), out_mx->M != cm->M.\n");

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrices */
  float ***alpha = ins_mx->dp; /* pointer to the alpha DP matrix */
  float ***beta  = out_mx->dp; /* pointer to the beta DP matrix */
  float ***post  = post_mx->dp; /* pointer to the post DP matrix */

  /* grow our post matrix, unless it's just the outside matrix, which is the correct size already */
  if(out_mx != post_mx) if((status = cm_hb_mx_GrowTo(cm, post_mx, errbuf, cp9b, L)) != eslOK) return status; 

  Lp = L - hdmin[0][L-jmin[0]];
  sc = alpha[0][L-jmin[0]][Lp];
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;  

  /* If local ends are on, start with the EL state (cm->M), otherwise
   * M deck is not valid. Note: there are no bands on the EL state */
  if (have_el) { 
    /* fill in the cm->M deck */
    for(j = 0; j <= L; j++) {
      for (d = 0; d <= j; d++) { 
	/*printf("M: j: %4d d: %4d ins: %10.6f out: %10.6f ", j, d, alpha[cm->M][j][d], beta[cm->M][j][d]);*/
	post[cm->M][j][d] = alpha[cm->M][j][d] + beta[cm->M][j][d] - sc;
	/*printf("post: %10.6f\n", post[cm->M][j][d]);*/
	/* convention is that alpha[cm->M] from Inside is always invalid, 
	 * because we can just as easily calculate it on the fly and save memory:
	 * alpha[cm->M][j][d] = cm->el_selfsc * d; for all j, d. */
      }
    }
  }
  
  for (v = (cm->M-1); v >= 0; v--) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      ESL_DASSERT1((j >= i0 && j <= j0));
      jp_v = j - jmin[v];
      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	dp_v = d - hdmin[v][jp_v];
	/*printf("v: %3d | jp_v: %3d | dp_v: %3d | alpha: %5.2f | beta: %5.2f\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	post[v][jp_v][dp_v] = alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v] - sc;
      }  
    }
  }
  return eslOK;
}


/*
 * Function: CMPosterior() 
 * Date:     EPN, Mon Nov 19 09:02:12 2007
 * Note:     based on Ian Holmes' P7EmitterPosterior() from HMMER's 2.x postprob.c
 *
 * Purpose:  Combines non-banded Inside and Outside matrices into a 
 *           posterior probability matrix. The value in post[v][j][d] 
 *           is the log of the posterior probability of a parse subtree rooted at v 
 *           emitting the subsequence i..j (i=j-d+1). 
 *           The caller must provide a <post> float matrix, but this matrix may 
 *           be the same matrix as that provided as Outside <out_mx>, 
 *           (overwriting it will not compromise the algorithm).
 *           
 * Args:     cm       - the model
 *           errbuf   - char buffer for reporting errors
 *           i0       - first position of target seq we're aligning, must be 1 
 *           j0       - final position of target seq we're aligning, L (length of seq) 
 *           ins_mx   - pre-calculated Inside matrix 
 *           out_mx   - pre-calculated Outside matrix
 *           post_mx  - pre-allocated matrix for Posteriors 
 *
 * Return:   eslOK on succes, eslEINCOMPAT on contract violation
 */
int
CMPosterior(CM_t *cm, char *errbuf, int i0, int j0, float ***ins_mx, float ***out_mx, float ***post_mx)
{
  int   v, j, d;
  float sc;
  int   vmax;
  int      L;    /* length of sequence */
  
  /* Contract check */
  if (ins_mx == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosterior(), ins_mx is NULL.\n");
  if (out_mx == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosterior(), out_mx is NULL.\n");
  if (post_mx == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosterior(), post_mx is NULL.\n");
  if (i0 != 1)            ESL_FAIL(eslEINCOMPAT, errbuf, "CMPosterior(), i0 != 1.\n");

  L = j0-i0+1;
  sc = ins_mx[0][L][L];

  /* If local ends are on, start with the EL state (cm->M), otherwise
   * its not a valid deck. */
  vmax = (cm->flags & CMH_LOCAL_END) ? cm->M : cm->M-1;
  for (v = vmax; v >= 0; v--) 
    for (j = 0; j <= L; j++) 
      for (d = 0; d <= j; d++)
	{ 
	  /*if(v == cm->M) printf("M: j: %4d d: %4d ins: %10.6f out: %10.6f ", j, d, ins_mx[v][j][d], out_mx[v][j][d]);*/
	  post_mx[v][j][d] = ins_mx[v][j][d] + out_mx[v][j][d] - sc;
	  /*if(v == cm->M) printf("post: %10.6f\n", post_mx[v][j][d]);*/
	}
  return eslOK;
}

/* 
 * Function: optimal_accuracy_align_hb()
 * Date:     EPN, Thu Nov 15 10:48:37 2007
 *
 * Purpose:  Run the Holmes/Durbin optimal accuracy algorithm 
 *           using bands in the j and d dimensions of the DP matrix. 
 *           Bands were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 *
 *           Two CM_HB_MX DP matrices must be passed in. The first
 *           <post_mx> must be pre-filled, containing posterior values
 *           from Inside/Outside runs on the target sequence. The
 *           second <mx> will be filled with the optimal accuracy
 *           scores, where:
 *
 *           mx[v][j][d] is the log of the sum of the posterior
 *                       probabilities of the residues i=j-d+1..j
 *                       in the subtree rooted at v (in the platonic,
 *                       non-banded matrix).
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitaized sequence [i0..j0]   
 *           L         - length of the dsq
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           mx        - the dp matrix to fill in, only cells within bands in cm->cp9b will 
 *                       be valid. 
 *           post_mx   - the pre-filled posterior matrix 
 *           ret_pp    - average posterior probability of aligned res i0..j0 in optimally accurate
 *                       parse.
 *
 * Returns: <ret_sc>, <ret_b>, <ret_bsc>, <ret_shadow>, see 'Args'
 * 
 * Throws:  <eslOK> on success.
 *          <eslERANGE> if required CM_HB_MX size exceeds CM_HB_MBLIMIT set in structs.h, in 
 *                      this case, alignment has been aborted, ret_* variables are not valid
 */
int
optimal_accuracy_align_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, void ****ret_shadow,  
			  int *ret_b, float *ret_bsc, CM_HB_MX *mx, CM_HB_MX *post_mx, float *ret_pp)
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
  int      have_el;     /* TRUE if we have local ends */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp_el;              /* offset j for non-banded EL state */
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
  float    pp;                 /* avg posterior probability of aligned res i0..j0 in optimally accurate parse */

  /* Contract check */
  if(dsq == NULL)                ESL_FAIL(eslEINCOMPAT, errbuf, "optimal_accuracy_align_hb(), dsq is NULL.\n");
  if (mx == NULL)                ESL_FAIL(eslEINCOMPAT, errbuf, "optimal_accuracy_align_hb(), mx is NULL.\n");
  if (cm->cp9b == NULL)          ESL_FAIL(eslEINCOMPAT, errbuf, "optimal_accuracy_align_hb(), cm->cp9b is NULL.\n");
  if (post_mx == NULL)           ESL_FAIL(eslEINCOMPAT, errbuf, "optimal_accuracy_align_hb(), cm->cp9b is NULL.\n");

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrices */
  float ***alpha = mx->dp; /* pointer to the alpha DP matrix, we'll store optimal parse in  */
  float ***post  = post_mx->dp; /* pointer to the alpha DP matrix, prefilled posterior values */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  /* grow the matrix based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo(cm, mx, errbuf, cp9b, W)) != eslOK) return status;

  /* initialize the EL deck, if it's valid */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;  
  if(have_el) { 
    for (jp_el = 0; jp_el <= W; jp_el++) {
      j = i0-1+jp_el;
      alpha[cm->M][j][0] = IMPOSSIBLE;
      for (d = 1; d <= jp_el; d++) {
	alpha[cm->M][j][d] = FLogsum(alpha[cm->M][j][d-1], post[cm->M][j][d]); /* optimal (and only) parse for EL is to emit all d residues */
      }
    }
  }

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

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
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
	/* initialize all valid cells for state v to IMPOSSIBLE */
	ESL_DASSERT1((! (NOT_IMPOSSIBLE(cm->endsc[v]))));
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
	if(have_el && NOT_IMPOSSIBLE(cm->endsc[v])) { 
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    /* special case: if hdmin[v][jp_v] > 0, we have to sum up posterior prob of emitting 
	     * hdmin[v][jp_v] residues from the EL state, this is then stored in alpha[v][jp_v][dp_v==0]
	     * note: when d = hdmin[v][jp_v], dp_v is 0
	     */
	    alpha[v][jp_v][0] = IMPOSSIBLE;
	    yshad[jp_v][0] = USED_EL; 
	    for(d = sd+1; d <= hdmin[v][jp_v]; d++) {
	      alpha[v][jp_v][0] = FLogsum(alpha[v][jp_v][0], post[cm->M][j-sdr][d-sd]);
	    }
	    /* now finish initializing the remaining valid d values */
	    for(d = hdmin[v][jp_v]+1; d <= hdmax[v][jp_v]; d++) {
	      dp_v = d - hdmin[v][jp_v];
	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v-1], post[cm->M][j-sdr][d-sd]); /* careful, we'll emit sd residues from v and d-sd from EL (i=((j-sdr)-d+1)..(j-sdr)) */
	      yshad[jp_v][dp_v] = USED_EL; 
	    }
	  }
	}
	else { /* cm->endsc[v] == IMPOSSIBLE, init all cells to IMPOSSIBLE */
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];

	    /* Check for special initialization case, specific to
	     * optimal_accuracy alignment, normally (with CYK for
	     * example) we init shadow matrix to USED_EL for all cells
	     * b/c we know that will be overwritten for the most
	     * likely transition, but with optimal accuracy, only
	     * emissions add to the score, so when d == sd, we know
	     * we'll emit sd residues from v, so the initialization
	     * will NOT be overwritten. We get around this by
	     * specially initializing cells for d == sd and v is a
	     * state with an END_E as a child, as the v->end
	     * transition. (v->end is always the yoffset == cnum[v]-1
	     * transition)
	     */
	    dp_v = 0;
	    for (d = hdmin[v][jp_v]; d <= sd; d++, dp_v++) { 
	      alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      yshad[jp_v][dp_v] = (cm->ndtype[cm->ndidx[v]+1] == END_nd) ? cm->cnum[v] - 1 : USED_EL;
	    }
	    for (d = ESL_MAX(hdmin[v][jp_v], sd+1); d <= hdmax[v][jp_v]; d++, dp_v++) {
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
	alpha[v][jp_v][0] = IMPOSSIBLE; 
      }
    }
    else if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL/IR states, loop nesting order is:
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
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd]) > alpha[v][jp_v][dp_v])
		{
		  alpha[v][jp_v][dp_v] = sc; 
		  yshad[jp_v][dp_v]    = yoffset;
		}
	    }
	  }
	  alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], post[v][jp_v][dp_v]);
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in alpha[v] are independent of each other, only depending on alpha[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
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
	    if((sc = alpha[y][jp_y_sdr][dp_y_sd]) > alpha[v][jp_v][dp_v]) {
	      alpha[v][jp_v][dp_v] = sc;
	      yshad[jp_v][dp_v]    = yoffset;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], post[v][jp_v][dp_v]);
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) 
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], FLogsum(post[v][jp_v][dp_v], post[v][jp_v][dp_v]));
	  /* note: for MP states, we're emitting 2 residues, include 2 * the posterior probability */
	}
	break;
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
	      
	      if ((sc = FLogsum(alpha[y][jp_y-k][dp_y - k], alpha[z][jp_z][kp_z])) 
		  > alpha[v][jp_v][dp_v]) { 
		alpha[v][jp_v][dp_v] = sc;
		kshad[jp_v][dp_v] = kp_z;
		/* note: we take the logsum here, because we're keeping track of the
		 * log of the summed probability of emitting all residues up to this
		 * point, (from i..j) from left subtree (i=j-d+1..j-k) and from the 
		 * right subtree. (j-k+1..j)
		 */
	      }
	    }
	  }
	}
      }
    }				/* finished calculating deck v. */

    /* Deal with possible local begins, these will most likely be off because local
     * ends are off (b/c we can't deal with local end emissions), but we can deal
     * with local begins.
     * We do this even though we can't deal with local ends, chances are that
     * local begins will be off also (since local ends must be), but this code
     * will work with local begins, so it's left here.
     */
    if(cm->flags & CMH_LOCAL_BEGIN) {
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
	  if (alpha[v][jp_v][Wp] > bsc) { 
	    b   = v;
	    bsc = alpha[v][jp_v][Wp];
	  }
	}
      }
      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling fast_alignT_hb() to use the b we return.
       */
      if (v == 0) { 
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp   = W - hdmin[v][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    if (bsc > alpha[0][jp_v][Wp]) {
	      alpha[0][jp_v][Wp] = bsc;
	      yshad[jp_v][Wp] = USED_LOCAL_BEGIN;
	    }
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp; fp = fopen("cyk.mx", "w"); cm_hb_mx_Dump(fp, mx); fclose(fp);*/
  
  Wp = W - hdmin[0][j0-jmin[0]];
  sc =     alpha[0][j0-jmin[0]][Wp];

  if (ret_b != NULL)      *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL)    *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */
  if (ret_shadow != NULL) *ret_shadow = shadow;
  else free_vjd_shadow_matrix(shadow, cm, i0, j0);
  free(yvalidA);

  /* convert score, a log probability, into the average posterior probability of all W aligned residues */
  pp = sreEXP2(sc) / (float) W;
  if(ret_pp != NULL) *ret_pp = pp;
  ESL_DPRINTF1(("optimal_accuracy_align_hb return pp: %f\n", pp));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
  return status; /* never reached */
}


/* 
 * Function: optimal_accuracy_align()
 * Date:     EPN, Sun Nov 18 20:45:22 2007
 *           
 * Purpose:  Run the Holmes/Durbin optimal accuracy algorithm 
 *           using bands in the j and d dimensions of the DP matrix. 
 *           Bands were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 * 
 *           Two float DP matrices must be passed in. The first
 *           <post_mx> must be pre-filled, containing posterior values
 *           from Inside/Outside runs on the target sequence. The
 *           second <mx> will be filled with the optimal accuracy
 *           scores, where:
 *
 *           mx[v][j][d] is the log of the sum of the posterior
 *                       probabilities of the residues i=j-d+1..j
 *                       in the subtree rooted at v. 
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitaized sequence [i0..j0]   
 *           L         - length of the dsq
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           ret_mx    - the main dp matrix, we'll allocate and fill it
 *           post_mx   - the pre-filled posterior matrix
 *           ret_pp    - average posterior probability of aligned res i0..j0 in optimally accurate
 *                       parse.
 *
 * Returns: <ret_sc>, <ret_b>, <ret_bsc>, <ret_mx>, <ret_shadow>, see 'Args'
 * 
 * Throws:  <eslOK> on success.                        
 *          <eslERANGE> if required size of DP matrix exceeds CM_NB_MBLIMIT set in structs.h, 
 *                      in this case, alignment has been aborted, ret_sc, ret_mx are not valid
 *                      Note: this shouldn't happen as post_mx (which is the same size as ret_mx
 *                      will be) was allocated after surviving the same size check in a previous function.
 */
int
optimal_accuracy_align(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, void ****ret_shadow,  
		       int *ret_b, float *ret_bsc, float ****ret_mx, float ***post_mx, float *ret_pp)
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
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      jp;          /* offset j, j = i0-1+jp */
  int      j_sdr;       /* j - sdr */
  int      d_sd;        /* d - sd */
  float    tsc;         /* a transition score */
  float    pp;          /* avg posterior probability of aligned res i0..j0 in optimally accurate parse */
  float ***alpha;       /* the DP matrix, we allocate here */
  float    Mb_for_alpha;/* megabytes needed for alpha matrix */
  int      have_el;     /* TRUE if we have local ends */

  /* Contract check */
  if (dsq == NULL)         ESL_FAIL(eslEINCOMPAT, errbuf, "optimal_accuracy_align(), dsq is NULL.\n");
  if (cm->cp9b == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "optimal_accuracy_align(), cm->cp9b is NULL.\n");
  if (post_mx == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "optimal_accuracy_align(), post_mx is NULL.\n");

  /* the pre-filled post matrix */
  float ***post  = post_mx;

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */

  /* allocate alpha (if it's small enough), allocate all decks, no deck reuse */
  Mb_for_alpha = ((float) size_vjd_deck(W, 1, W) * ((float) (cm->M-1)));
  if(Mb_for_alpha > CM_NB_MX_MB_LIMIT) 
    ESL_FAIL(eslERANGE, errbuf, "optimal_accuracy_align(), requested size of non-banded DP matrix %.2f Mb > %.2f Mb limit (CM_NB_MX_MB_LIMIT from structs.h).", Mb_for_alpha, (float) CM_NB_MX_MB_LIMIT);
  ESL_DPRINTF1(("Size of alpha matrix: %.2f\n", Mb_for_alpha));

  ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
  for (v = 0; v <= cm->M; v++) alpha[v] = alloc_vjd_deck(L, i0, j0);

  /* initialize the EL deck, if it's valid */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;  
  if(have_el) { 
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      alpha[cm->M][j][0] = IMPOSSIBLE;
      for (d = 1; d <= jp; d++) {
	alpha[cm->M][j][d] = FLogsum(alpha[cm->M][j][d-1], post[cm->M][j][d]); /* optimal (and only) parse for EL is to emit all d residues */
      }
    }
  }

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

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);

    /* Get a shadow deck to fill in and initialize all valid cells for state v */
    if (cm->sttype[v] != E_st) {
      if (cm->sttype[v] == B_st) {
	kshad     = alloc_vjd_kshadow_deck(L, i0, j0); 
	shadow[v] = (void **) kshad;
	/* initialize all valid cells for state v to IMPOSSIBLE (local ends are impossible for B states) */
	ESL_DASSERT1((! (NOT_IMPOSSIBLE(cm->endsc[v]))));
	ESL_DASSERT1((cm->ndtype[cm->ndidx[v]+1] != END_nd));
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  for (d = 0; d <= jp; d++) {
	    alpha[v][j][d] = IMPOSSIBLE;
	    kshad[j][d] = USED_EL; 
	  }
	}
      } else { /* ! B_st && ! E_st */
	yshad     = alloc_vjd_yshadow_deck(L, i0, j0);
	shadow[v] = (void **) yshad;
	/* initialize all valid cells for state v */
	if(have_el && NOT_IMPOSSIBLE(cm->endsc[v])) {
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= sd && d <= jp; d++) { 
	      alpha[v][j][d] = IMPOSSIBLE;
	      yshad[j][d] = USED_EL; 
	    }
	    for (d = sd+1; d <= jp; d++) {
	      alpha[v][j][d] = FLogsum(alpha[v][j][d-1], post[cm->M][j-sdr][d-sd]); /* careful, we'll emit sd residues from v and d-sd from EL (i=((j-sdr)-d+1)..(j-sdr)) */
	      yshad[j][d] = USED_EL; 
	      /* printf("I alpha[v: %d][j: %d][k: %d]: %.4f\n", v, j, d, alpha[v][j][d]); */
	    }
	  }
	}
	else { /* cm->endsc[v] == IMPOSSIBLE */
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;

	    /* Handle a special initialization case, specific to
	     * optimal_accuracy alignment. Normally (with CYK for
	     * example) we init shadow matrix to USED_EL for all cells
	     * b/c we know that USED_EL will be overwritten for the most
	     * likely transition, but with optimal accuracy, only
	     * emissions add to the score, so when d <= sd, we know
	     * we'll emit sd residues from v, so the initialization
	     * will NOT be overwritten. We get around this by
	     * specially initializing cells for d == sd and v is a
	     * state with an END_E as a child, as the v->end
	     * transition. (v->end is always the yoffset == cnum[v]-1
	     * transition)
	     */
	    for (d = 0; d <= sd; d++) { 
	      alpha[v][j][d] = IMPOSSIBLE;
	      yshad[j][d] = (cm->ndtype[cm->ndidx[v]+1] == END_nd) ? cm->cnum[v] - 1 : USED_EL;
	    }
	    for (d = sd+1; d <= jp; d++) {
	      alpha[v][j][d] = IMPOSSIBLE;
	      yshad[j][d] = USED_EL; 
	    }
	  }
	}
      }
    }
    
    if(cm->sttype[v] == E_st) { 
      for (jp = 0; jp <= W; jp++) {
	j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
	alpha[v][j][0] = IMPOSSIBLE; /* we haven't seen any residues yet, score is IMPOSSIBLE */
	for (d = 1; d <= jp; d++) alpha[v][j][d] = IMPOSSIBLE;
      }
    }
    else if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      /* update alpha[v][j][d] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i    = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = alpha[y][j_sdr][d_sd]) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc; 
	      yshad[j][d]    = yoffset;
	    }
	  }
	  alpha[v][j][d] = FLogsum(alpha[v][j][d], post[v][j][d]);
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in alpha[v] are independent of each other, only depending on alpha[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];

	for (jp = sdr; jp <= W; jp++) {
	  j = i0-1+jp;
	  j_sdr = j - sdr;

	  for (d = sd; d <= jp; d++) {
	    if((sc = alpha[y][j_sdr][d - sd]) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc;
	      yshad[j][d]    = yoffset;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
      case MR_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  i = j - 1;
	  for (d = sd; d <= jp; d++) 
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], post[v][j][d]);
	}
	break;
      case MP_st:
	for (jp = 0; jp <= W; jp++) {
	  j = i0-1+jp;
	  i = j - 1;
	  for (d = sd; d <= jp; d++)
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], FLogsum(post[v][j][d], post[v][j][d]));
	  /* note: for MP states, we're emitting 2 residues, include 2 * the posterior probability */
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++)
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      for (jp = 0; jp <= W; jp++) { 
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  for (k = 0; k <= d; k++) {
	    if ((sc = FLogsum(alpha[y][j-k][d-k], alpha[z][j][k])) > alpha[v][j][d]) { 
	      alpha[v][j][d] = sc;
	      kshad[j][d]    = k;
	      /* note: we take the logsum here, because we're keeping track of the
	       * log of the summed probability of emitting all residues up to this
	       * point, (from i..j) from left subtree (i=j-d+1..j-k) and from the 
	       * right subtree. (j-k+1..j)
	       */
	    }
	  }
	}
      }
    }
				/* finished calculating deck v. */
      
    if(cm->flags & CMH_LOCAL_BEGIN) {
      if (alpha[v][j0][W] > bsc) {
	b   = v;
	bsc = alpha[v][j0][W];
      }
    }
    /* Check for whether we need to store an optimal local begin score
     * as the optimal overall score, and if we need to put a flag
     * in the shadow matrix telling fast_alignT() to use the b we return.
     * We do this even though we can't deal with local ends, chances are that
     * local begins will be off also (since local ends must be), but this code
     * will work with local begins, so it's left here.
     */
    if(v == 0 && bsc > alpha[0][j0][W]) {
      alpha[0][j0][W] = bsc;
      yshad[j0][W] = USED_LOCAL_BEGIN;
    }
  } /* end loop over all v */
  
  sc =     alpha[0][j0][W];
  if (ret_b != NULL)      *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL)    *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */
  if (ret_shadow != NULL) *ret_shadow = shadow;
  else free_vjd_shadow_matrix(shadow, cm, i0, j0);
  if (ret_mx     != NULL) *ret_mx = alpha;
  else free_vjd_matrix(alpha, cm->M, 1, L);

  /* convert score, a log probability, into the average posterior probability of all W aligned residues */
  pp = sreEXP2(sc) / (float) W;

  ESL_DPRINTF1(("optimal_accuracy_align return pp: %f\n", pp));
  if(ret_pp != NULL) *ret_pp = pp;


  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
  return status; /* never reached */
}


/*
 * Function: SampleFromInside()
 * Incept:   EPN, Thu Nov 15 16:45:32 2007
 *          
 * Purpose:  Sample a parsetree from a non-banded float Inside matrix.
 *           
 * Args:     r        - source of randomness
 *           cm       - the model
 *           errbuf   - char buffer for reporting errors
 *           dsq      - digitized sequence
 *           L        - length of dsq, alpha *must* go from 1..L
 *           mx       - pre-calculated Inside matrix (floats)
 *           ret_tr   - ptr to parsetree we'll return (*must* be non-NULL)
 *           ret_sc   - score of sampled parsetree; dies immediately with cm_Fail() if an error occurs.
 * 
 * Return:   <ret_tr>, <ret_sc> see 'Args'.
 *           eslOK on succes, eslEINCOMPAT on contract violation, eslEMEM on memory allocation error
 *           eslFAIL on unexpeced error
 */
int
SampleFromInside(ESL_RANDOMNESS *r, CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float ***mx, Parsetree_t **ret_tr, float *ret_sc)
{
  int          status;             /* easel status code */
  int          v, y, z, b;         /* state indices */
  int          yoffset;            /* transition offset in a states transition vector */
  int          i, j;               /* sequence position indices */
  int          d;                  /* j - i + 1; the current subseq length */
  int          k;                  /* right subseq fragment length for bifurcs */
  int          nd;                 /* node index */
  int          bifparent;          /* for connecting bifurcs */
  Parsetree_t *tr;                 /* trace we're building */
  ESL_STACK   *pda;                /* the stack */
  float        pvec[MAXCONNECT+1]; /* prob vector of possible paths to take, (max num children + 1 for possibility of EL) */
  float       *bifvec;             /* pvec for choosing transition out of BIF_B states */
  float       *rootvec;            /* pvec for choosing transition out of ROOT_S if local begins are on */
  float        maxsc;              /* max score in our vector of scores of possible subparses */
  int          el_is_possible;     /* TRUE if we can jump to EL from current state (and we're in local mode) FALSE if not */
  int          ntrans;             /* number of transitions for current state */
  float        fsc = 0.;           /* score of the parsetree we're sampling */

  /* contract check */
  if(ret_tr == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SampleFromInside(), ret_tr is NULL.");
  if(r      == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SampleFromInside(), source of randomness r is NULL.");
  if(mx     == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SampleFromInside(), source of randomness r is NULL.");
  
  float ***alpha = mx;

  /* initialize pvec */
  esl_vec_FSet(pvec, (MAXCONNECT+1), 0.);

  /* Create a parse tree structure and initialize it by adding the root state.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Stochastically traceback through the Inside matrix 
   * this section of code is stolen and adapted from cm_dpsmall.c:insideT() 
   */
  pda = esl_stack_ICreate();
  v = 0;

  j = d = L;
  i = 1;
  fsc = 0.;
  while (1) {
    if (cm->sttype[v] == B_st) {
      y = cm->cfirst[v];
      z = cm->cnum[v];

      ESL_ALLOC(bifvec, sizeof(float) * (d+1));
      /* set bifvec[] as (float-ized) log odds scores for each valid left fragment length */
      for(k = 0; k <= d; k++) 
	bifvec[k] = alpha[y][j-k][d-k] + alpha[z][j][k];
      maxsc = esl_vec_FMax(bifvec, (d+1));
      esl_vec_FIncrement(bifvec, (d+1), (-1. * maxsc));
      esl_vec_FScale(bifvec, (d+1), log(2));
      esl_vec_FLogNorm(bifvec, (d+1));
      k = esl_rnd_FChoose(r, bifvec, (d+1));
      free(bifvec);

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
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
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
    } else {
      if((v > 0) || (! (cm->flags & CMH_LOCAL_BEGIN))) /* ROOT_S with local begins on is a special case that we handle below */
	{ 
	  /* choose which transition we take */
	  esl_vec_FSet(pvec, (MAXCONNECT+1), IMPOSSIBLE); /* not really necessary */
	  fsc += get_femission_score(cm, dsq, v, i, j); 
	  
	  /* set pvec[] as (float-ized) log odds scores for each child we can transit to, 
	   * plus a local end (if possible) */
	  ntrans = cm->cnum[v];
	  el_is_possible = FALSE;
	  if((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) { 
	    el_is_possible = TRUE; 
	    ntrans++; 
	  }
	  for(yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = yoffset + cm->cfirst[v];
	    pvec[yoffset] = cm->tsc[v][yoffset] + 
	      alpha[y][j - StateRightDelta(cm->sttype[v])][d - StateDelta(cm->sttype[v])];
	  }
	  if(el_is_possible) pvec[cm->cnum[v]] = cm->endsc[v] + 
			       alpha[cm->M][j][d]; /* EL is silent when we transition into it from non-EL */
	  /* note: we can treat the log odds scores as log probs, because
	   * the log probability of the null model is the same for each,
	   * so essentially we've divided each score by the same constant, so 
	   * the *relative* proportion of the log odds scores is the
	   * same as the relative proportion of the log probabilities (seq | model) */
	  
	  maxsc = esl_vec_FMax(pvec, ntrans);
	  esl_vec_FIncrement(pvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (pvec, ntrans, log(2));
	  esl_vec_FLogNorm(pvec, ntrans);
	  yoffset = esl_rnd_FChoose(r, pvec, ntrans);
	  if(yoffset < cm->cnum[v]) fsc += cm->tsc[v][yoffset]; 
	  else {
	    fsc += cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
	    yoffset = USED_EL; /* we chose EL */
	  }
	}
      else /* v == 0 && (cm->flags && CMH_LOCAL_BEGIN) ( local begins are on )*/
	{
	  ntrans = cm->M; /* pretend all states are possible to begin into, but they're not as some will remain IMPOSSIBLE */
	  ESL_ALLOC(rootvec, sizeof(float) * (ntrans));
	  esl_vec_FSet(rootvec, ntrans, IMPOSSIBLE);
	  rootvec[cm->nodemap[1]] = cm->beginsc[cm->nodemap[1]] + alpha[cm->nodemap[1]][j][d]; /* ROOT_S is silent */
	  for (nd = 2; nd < cm->nodes; nd++) {
	    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
		cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
	      {
		rootvec[cm->nodemap[nd]] = cm->beginsc[cm->nodemap[nd]] + alpha[cm->nodemap[nd]][j][d]; /* ROOT_S is silent */
	      }
	  }
	  /* this block is shared with v > 0 block, but we repeat it here so we don't need another if statement */
	  maxsc = esl_vec_FMax(rootvec, ntrans);
	  esl_vec_FIncrement(rootvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (rootvec, ntrans, log(2));
	  esl_vec_FLogNorm(rootvec, ntrans);
	  b = esl_rnd_FChoose(r, rootvec, ntrans);
	  /* end of similar block with v > 0 */
	  fsc += cm->beginsc[b];
	  yoffset = USED_LOCAL_BEGIN; 
	  free(rootvec); /* we will not need this again */
	}

      /*printf("v : %d | r : %d | z : %d | 1 : %d | \n", v, r, z, 1);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    ESL_FAIL(eslEINCONCEIVABLE, errbuf, "'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */

  *ret_tr = tr; /* contract checked ret_tr was non-NULL */
  if(ret_sc != NULL) *ret_sc = fsc;
  ESL_DPRINTF1(("SampleFromInside() return sc: %f\n", fsc));
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "memory error.");
}

/*
 * Function: SampleFromInsideHB()
 * Incept:   EPN, Fri Sep  7 11:02:15 2007
 *          
 * Purpose:  Sample a parsetree from a HMM banded Inside matrix.
 *           
 * Args:     r        - source of randomness
 *           cm       - the model
 *           errbuf   - char buffer for reporting errors
 *           dsq      - digitized sequence
 *           L        - length of dsq, alpha *must* go from 1..L
 *           mx       - pre-calculated Inside matrix
 *           ret_tr   - ptr to parsetree we'll return (*must* be non-NULL)
 *           ret_sc   - score of sampled parsetree; dies immediately with cm_Fail() if an error occurs.
 * 
 * Return:   <ret_tr>, <ret_sc> see 'Args'.
 *           eslOK on succes, eslEINCOMPAT on contract violation, eslEMEM on memory allocation error
 *           eslFAIL on unexpeced error
 */
int
SampleFromInsideHB(ESL_RANDOMNESS *r, CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, CM_HB_MX *mx, Parsetree_t **ret_tr, float *ret_sc)
{
  int          status;             /* easel status code */
  int          v, y, z, b;         /* state indices */
  int          yoffset;            /* transition offset in a states transition vector */
  int          i, j;               /* sequence position indices */
  int          jp_v, jp_y, jp_z;   /* positions, offset inside j band */
  int          kmin, kmax;         /* min/max k in current d band */
  int          d;                  /* j - i + 1; the current subseq length */
  int          dp_v, dp_y;         /* length, offset inside a d band */
  int          k;                  /* right subseq fragment length for bifurcs */
  int          kp_z;               /* right fragment length, offset inside a d band */
  int          nd;                 /* node index */
  int          bifparent;          /* for connecting bifurcs */
  Parsetree_t *tr;                 /* trace we're building */
  ESL_STACK   *pda;                /* the stack */
  float        pvec[MAXCONNECT+1]; /* prob vector of possible paths to take, (max num children + 1 for possibility of EL) */
  float       *bifvec;             /* pvec for choosing transition out of BIF_B states */
  float       *rootvec;            /* pvec for choosing transition out of ROOT_S if local begins are on */
  float        maxsc;              /* max score in our vector of scores of possible subparses */
  int          el_is_possible;     /* TRUE if we can jump to EL from current state (and we're in local mode) FALSE if not */
  int          ntrans;             /* number of transitions for current state */
  float        fsc = 0.;           /* score of the parsetree we're sampling */
  int          seen_valid;         /* for checking we have at least one valid path to take  */
  int          sd;                 /* state delta for current state, residues emitted left + residues emitted right */
  int          sdr;                /* state right delta for current state, residues emitted right */

  /* contract check */
  if(ret_tr == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SampleFromInsideHB(), ret_tr is NULL.");
  if(r      == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SampleFromInsideHB(), source of randomness r is NULL.");
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SampleFromInsideHB(), cm->cp9b is NULL.\n");
  if (mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SampleFromInsideHB(), mx is NULL.\n");

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrix */
  float ***alpha = mx->dp; /* pointer to the alpha DP matrix */
  
  /* initialize pvec */
  esl_vec_FSet(pvec, (MAXCONNECT+1), 0.);

  /* Create a parse tree structure and initialize it by adding the root state.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Stochastically traceback through the Inside matrix 
   * this section of code is stolen and adapted from hbandcyk.c:insideTHB() 
   */
  pda = esl_stack_ICreate();
  v = 0;

  j = d = L;
  i = 1;
  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];
  fsc  = 0.;
  while (1) {
    if(cm->sttype[v] != EL_st && d > hdmax[v][jp_v]) ESL_FAIL(eslFAIL, errbuf, "ERROR in SampleFromInsideHB(). d : %d > hdmax[%d] (%d)\n", d, v, hdmax[v][jp_v]);
    if(cm->sttype[v] != EL_st && d < hdmin[v][jp_v]) ESL_FAIL(eslFAIL, errbuf, "ERROR in SampleFromInsideHB(). d : %d < hdmin[%d] (%d)\n", d, v, hdmin[v][jp_v]);

    if (cm->sttype[v] == B_st) {
      y = cm->cfirst[v];
      z = cm->cnum[v];
      jp_z = j-jmin[z];
      k = kp_z + hdmin[z][jp_z];  /* k = offset len of right fragment */

      ESL_ALLOC(bifvec, sizeof(float) * (d+1));
      /* set bifvec[] as (float-ized) log odds scores for each valid left fragment length,
       * we have to be careful to check that the corresponding alpha cell for each length is valid  */
      esl_vec_FSet(bifvec, (d+1), IMPOSSIBLE); /* only valid d's will be reset to a non-IMPOSSIBLE score */

      /* This search for valid k's is complex, and uncommented. It was taken from
       * cm_dpalign.c:fast_cyk_align_hb(), the B_st case. The code there is commented somewhat
       * extensively. I'm pretty sure this is the most efficient (or at least close to it) 
       * way to find the valid cells in the DP matrix we're looking for. 
       */
      jp_v = j - jmin[v];
      jp_y = j - jmin[y];
      jp_z = j - jmin[z];
      if(j < jmin[v] || j > jmax[v])               ESL_FAIL(eslFAIL, errbuf, "SampleFromInsideHB() B_st v: %d j: %d outside band jmin: %d jmax: %d\n", v, j, jmin[v], jmax[v]);
      if(d < hdmin[v][jp_v] || d > hdmax[v][jp_v]) ESL_FAIL(eslFAIL, errbuf, "SampleFromInsideHB() B_st v: %d j: %d d: %d outside band dmin: %d dmax: %d\n", v, j, d, hdmin[v][jp_v], hdmax[v][jp_v]);
      seen_valid = FALSE;
      kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
      kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
      for(k = kmin; k <= kmax; k++)
	{
	  if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
	    {
	      kp_z = k-hdmin[z][jp_z];
	      dp_y = d-hdmin[y][jp_y-k];
	      bifvec[k] = alpha[y][jp_y-k][dp_y-k] + alpha[z][jp_z][kp_z]; 
	      seen_valid = TRUE;
	    }
	}
      if(!seen_valid) ESL_FAIL(eslFAIL, errbuf, "SampleFromInsideHB() number of valid transitions (for a B_st) is 0. You thought this was impossible.");
      maxsc = esl_vec_FMax(bifvec, (d+1));
      esl_vec_FIncrement(bifvec, (d+1), (-1. * maxsc));
      esl_vec_FScale(bifvec, (d+1), log(2));
      esl_vec_FLogNorm(bifvec, (d+1));
      k = esl_rnd_FChoose(r, bifvec, (d+1));
      free(bifvec);

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
      if((v > 0) || (! (cm->flags & CMH_LOCAL_BEGIN))) /* ROOT_S with local begins on is a special case that we handle below */
	{ 
	  /* Choose which transition we take.
	   * Set pvec[] as (float-ized) log odds scores for each child we can transit to, 
	   * plus a local end (if possible). We only want to look at valid transitions, that
	   * is, those that do not violate the bands (correspond to accessing cells that actually
	   * exist in the DP matrix). 
	   */
	  seen_valid = FALSE;
	  esl_vec_FSet(pvec, (MAXCONNECT+1), IMPOSSIBLE); /* only transitions that correspond to valid cells will be reset to a non-IMPOSSIBLE score */
	  fsc += get_femission_score(cm, dsq, v, i, j); 
	  sdr = StateRightDelta(cm->sttype[v]);
	  sd  = StateDelta(cm->sttype[v]);
	  for(yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
	    {
	      y = yoffset + cm->cfirst[v];
	      if((j - sdr) >= jmin[y] && (j - sdr) <= jmax[y]) 
		{ /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if((d - sd) >= hdmin[y][jp_y-sdr] && (d - sd) <= hdmax[y][jp_y-sdr])
		    {
		      dp_y = d - hdmin[y][(jp_y - sdr)];  /* d index for state y 
							     in alpha w/mem eff bands */
		      /* if we get here alpha[y][jp_y-sdr][dp_y-sd] is a valid alpha cell
		       * corresponding to alpha[y][j-sdr][d-sd] in the platonic matrix.
		       */
		      pvec[yoffset] = cm->tsc[v][yoffset] + alpha[y][jp_y - sdr][dp_y - sd];
		      seen_valid = TRUE;
		    }
		}		
	    }
	  if(!seen_valid) {
	    ESL_FAIL(eslFAIL, errbuf, "SampleFromInsideHB() number of valid transitions is 0. You thought this was impossible.");
	  }
	  if((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) 
	    el_is_possible = TRUE; 
	  else 
	    el_is_possible = FALSE;
	  if(el_is_possible) pvec[cm->cnum[v]] = cm->endsc[v] + alpha[cm->M][j][d]; /* EL is silent when we transition into it from non-EL */
	  ntrans = cm->cnum[v] + el_is_possible;
	  maxsc = esl_vec_FMax(pvec, ntrans);
	  esl_vec_FIncrement(pvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (pvec, ntrans, log(2));
	  esl_vec_FLogNorm(pvec, ntrans);
	  yoffset = esl_rnd_FChoose(r, pvec, ntrans);
	  if(yoffset < cm->cnum[v]) fsc += cm->tsc[v][yoffset]; 
	  else {
	    fsc += cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
	    yoffset = USED_EL; /* we chose EL */
	  }
	}
      else /* v == 0 && (cm->flags && CMH_LOCAL_BEGIN) ( local begins are on )*/
	{
	  seen_valid = FALSE;
	  ntrans = cm->M; /* pretend all states are possible to begin into, but they're not as some will remain IMPOSSIBLE */
	  ESL_ALLOC(rootvec, sizeof(float) * (ntrans));
	  esl_vec_FSet(rootvec, ntrans, IMPOSSIBLE);

	  /* Set all the legal states that we can local begin into to appropriate scores.
	   * Only states y that have a non-zero cm->beginsc[y] AND have alpha[y][j][d]
	   * within their bands are legal.
	   */
	  for (nd = 1; nd < cm->nodes; nd++) {
	    if ((nd == 1) || /* we can transit into node 1 no matter what */
		(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
		 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd))
	      {
		y = cm->nodemap[nd];
		if(j >= jmin[y] && j <= jmax[y]) 
		  { /* enforces j is valid for state y */
		    jp_y = j - jmin[y];
		    if(d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y])
		      {
			dp_y = d - hdmin[y][jp_y];
			rootvec[y] = cm->beginsc[y] + alpha[y][jp_y][dp_y]; /* ROOT_S is silent */
			seen_valid = TRUE;
		      }
		  }
	      }
	  }
	  if(!seen_valid) ESL_FAIL(eslFAIL, errbuf, "SampleFromInsideHB() number of valid transitions (from ROOT_S!) is 0. You thought this was impossible.");
	  /* this block is shared with v > 0 block, but we repeat it here so we don't need another if statement */
	  maxsc = esl_vec_FMax(rootvec, ntrans);
	  esl_vec_FIncrement(rootvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (rootvec, ntrans, log(2));
	  esl_vec_FLogNorm(rootvec, ntrans);
	  b = esl_rnd_FChoose(r, rootvec, ntrans);
	  /* end of similar block with v > 0 */
	  fsc += cm->beginsc[b];
	  yoffset = USED_LOCAL_BEGIN; 
	  free(rootvec); /* we will not need this again */
	}

      /*printf("v : %d | r : %d | z : %d | 1 : %d | \n", v, r, z, 1);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    ESL_FAIL(eslEINCONCEIVABLE, errbuf, "'Inconceivable!'\n'You keep using that word...'");
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
  if(ret_sc != NULL) *ret_sc = fsc;
  ESL_DPRINTF1(("SampleFromInsideHB() return sc: %f\n", fsc));
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "memory error.");
}



/*
 * Function: get_femission_score()
 * Incept:   EPN, Thu Nov 15 16:48:56 2007
 *          
 * Purpose:  Given a CM, dsq, state index and coordinates return the float emission
 *           score.
 *           
 * Args:     cm       - the model
 *           dsq      - digitized sequence
 *           v        - state index
 *           i        - dsq index for first position of subseq for subtree at v
 *           j        - dsq index for last position of subseq for subtree at v
 *
 * Return:   float emission score, 0 if state is non-emitter.
 */
float
get_femission_score(CM_t *cm, ESL_DSQ *dsq, int v, int i, int j)
{
  if     (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) return cm->oesc[v][dsq[i]];
  else if(cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) return cm->oesc[v][dsq[j]];
  else if(cm->sttype[v] == MP_st)                           return cm->oesc[v][dsq[i]*cm->abc->Kp+dsq[j]];
  else return 0.;
}


/*
 * Function: CMCheckPosteriorHB()
 * Date:     EPN, Fri Nov 16 14:35:43 2007      
 *
 * Purpose:  Given a HMM banded posterior probability cube, 
 *           check to make sure that for each residue k of the 
 *           sequence: \sum_v p(v | k emitted from v) = 1.0
 *           To check this, we have to allow possibility that 
 *           the res at posn k was emitted from a left 
 *           emitter or a right emitter.
 *           
 * Args:     cm       - the model
 *           errbuf   - char buffer for returning error messages with ESL_FAIL
 *           i0       - first residue to check
 *           j0       - last residue to check
 *           post     - pre-filled dynamic programming cube
 *           
 * Return:   eslOK on success, eslFAIL if any residue fails check
 */
int
CMCheckPosteriorHB(CM_t *cm, char *errbuf, int i0, int j0, CM_HB_MX *post)
{
  int   v, j, d, k;
  int   jp_v, dp_v, kp_v;
  float sc;

  /* Contract check */
  if (post == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMCheckPosteriorHB(), post is NULL.\n");
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "CMCheckPosteriorHB(), cm->cp9b is NULL.\n");

  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  for (k = i0; k <= j0; k++) {
    sc = IMPOSSIBLE;
    for (v = (cm->M - 1); v >= 0; v--) {
      if((cm->sttype[v] == MP_st) ||
	 (cm->sttype[v] == ML_st) ||
	 (cm->sttype[v] == IL_st)) {
	for (j = k; j <= j0; j++) {
	  if(j >= jmin[v] && j <= jmax[v]) { 
	    jp_v = j - jmin[v]; 
	    d    = j-k+1;
	    if(d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]) {
	      dp_v = d - hdmin[v][jp_v];
	      sc = FLogsum(sc, (post->dp[v][jp_v][dp_v]));
	      /*printf("sc: %.4f added L v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", sc, v, (j-d+1), j, d, post->dp[v][jp_v][dp_v]); */
	    }
	  }
	}
      }
      if(k >= jmin[v] && k <= jmax[v]) {
	kp_v = k - jmin[v];
	if((cm->sttype[v] == MP_st) ||
	   (cm->sttype[v] == MR_st) ||
	   (cm->sttype[v] == IR_st)) {
	  for (d = 1; d <= k; d++) { 
	    if(d >= hdmin[v][kp_v] && d <= hdmax[v][kp_v]) {
	      dp_v = d - hdmin[v][kp_v];
	      sc = FLogsum(sc, (post->dp[v][kp_v][dp_v]));
	      /*printf("sc: %.4f added R v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", sc, v, (k-d+1), k, d, post->dp[v][kp_v][dp_v]); */
	    }
	  }
	}
      }
    }
    /* Finally factor in possibility of a local end, i.e. that the EL state
     * may have "emitted" this residue.
     */
    if (cm->flags & CMH_LOCAL_END) {
      for (j = k; j <= j0; j++) {
	d = j-k+1;
	sc = FLogsum(sc, (post->dp[cm->M][j][d]));
	/*printf("sc: %.4f added EL v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", sc, cm->M, (j-d+1), j, d, post->dp[cm->M][j][d]); */
      }
    }
    if(((sc - 0.) > 0.01) || ((sc - 0.) < -0.01))
      ESL_FAIL(eslFAIL, errbuf, "residue position %d has summed prob of %5.4f (2^%5.4f) in posterior cube.\n", k, (sreEXP2(sc)), sc);
    /*printf("k: %d | total: %10.2f\n", k, (sreEXP2(sc)));*/
  }  
  ESL_DPRINTF1(("CMCheckPosteriorHB() passed, all residues have summed probability of emission of 1.0.\n"));
  return eslOK;
}

/*
 * Function: CMCheckPosterior()
 * Date:     EPN 05.25.06 
 *
 * Purpose:  Given a posterior probability cube, check to make
 *           sure that for each residue k of the sequence:
 *           \sum_v p(v | k emitted from v) = 1.0
 *           To check this, we have to allow possibility that 
 *           the res at posn k was emitted from a left 
 *           emitter or a right emitter.
 *           
 * Args:     cm       - the model
 *           errbuf   - char buffer for returning error messages with ESL_FAIL
 *           i0       - first residue to check
 *           j0       - last residue to check
 *           post     - pre-filled dynamic programming cube
 *           
 * Return:   eslOK on success, eslFAIL if any residue fails check
 */
int 
CMCheckPosterior(CM_t *cm, char *errbuf, int i0, int j0, float ***post)
{
  float sc;
  int   k, v, j, d;

  /* contract check */
  if (post == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "CMCheckPosterior(), post is NULL.\n");

  for (k = i0; k <= j0; k++) {
    sc = IMPOSSIBLE;
    for (v = (cm->M - 1); v >= 0; v--) {
      {
	if((cm->sttype[v] == MP_st) ||
	   (cm->sttype[v] == ML_st) ||
	   (cm->sttype[v] == IL_st)) {
	  for (j = k; j <= j0; j++) {
	    d = j-k+1;
	    sc = FLogsum(sc, (post[v][j][d]));
	    /* printf("sc: %.4f added L v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", sc, v, (j-d+1), j, d, post[v][j][d]); */
	  }
	}
	if((cm->sttype[v] == MP_st) ||
	   (cm->sttype[v] == MR_st) ||
	   (cm->sttype[v] == IR_st)) {
	  for (d = i0; d <= k; d++) {
	    sc = FLogsum(sc, (post[v][k][d]));
	    /* printf("sc: %.4f added R v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", sc, v, (k-d+1), k, d, post[v][k][d]); */
	  }
	}
      }
    }
    /* Finally factor in possibility of a local end, i.e. that the EL state
     * may have "emitted" this residue.
     */
    if (cm->flags & CMH_LOCAL_END) {
      for (j = k; j <= j0; j++) {
	d = j-k+1;
	sc = FLogsum(sc, (post[cm->M][j][d]));
	/* printf("sc: %.4f added EL v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", sc, cm->M, (j-d+1), j, d, post[cm->M][j][d]); */
      }
    }
    if(((sc - 0.) > 0.01) || ((sc - 0.) < -0.01))
	ESL_FAIL(eslFAIL, errbuf, "residue position %d has summed prob of %5.4f (2^%5.4f) in posterior cube.\n", k, (sreEXP2(sc)), sc);
    /*printf("k: %d | total: %10.2f\n", k, (sreEXP2(sc)));*/
  }  
  ESL_DPRINTF1(("CMCheckPosterior() passed, all residues have summed probability of emission of 1.0.\n"));
  return eslOK;
}

/* Function: CMPostalCode()
 * Date:     EPN 05.25.06 based on SRE's PostalCode() 
 *           from HMMER's postprob.c
 *
 * Purpose:  Given a parse tree and a posterior
 *           probability cube, calculate two strings that
 *           represents the confidence values on each 
 *           residue in the sequence. 
 *           
 *           The code strings is 0..L-1  (L = len of target seq),
 *           so it's in the coordinate system of the sequence string;
 *           off by one from dsq; and convertible to the coordinate
 *           system of aseq using MakeAlignedString().
 *           
 *           Values are 00-99,**  
 *           for example, 93 means with >=93% posterior probabiility,
 *           residue i is aligned to the state k that it
 *           is assigned to in the given trace.
 *
 *           Because we have 2 digit precision, we need two
 *           strings, the first will be the 'tens' place of
 *           the posterior probability, '9' for the 93% example,
 *           and the second string will hold the 'ones' place,
 *           the '3' in the 93% example.
 *
 * Args:     L    - length of seq
 *           post - posterior prob cube: see CMPosterior()
 *           *tr  - parsetree to get a Postal code string for.   
 *           ret_pcode1 - 'tens' place postal code string ('9' for 93)
 *           ret_pcode2 - 'ones' place postal code string ('3' for 93)
 * Returns:  void
 *
 */
int
Fscore2postcode(float sc)
{
  int i;
  i = (int) (FScore2Prob(sc, 1.) * 100.);
  ESL_DASSERT1((i >= 0 && i <= 100)); 
  return i;
}

/* Function: FScore2Prob()
 * 
 * Purpose:  Convert a float log_2 odds score back to a probability;
 *           needs the null model probability, if any, to do the conversion.
 */
float 
FScore2Prob(float sc, float null)
{
  /*printf("in FScore2Prob: %10.2f sreEXP2: %10.2f\n", sc, (sreEXP2(sc)));*/
  if (!(NOT_IMPOSSIBLE(sc))) return 0.;
  else                       return (null * sreEXP2(sc));
}


void
CMPostalCode(CM_t *cm, int L, float ***post, Parsetree_t *tr, char **ret_pcode1, char **ret_pcode2)
{
  int status;
  int x, v, i, j, d, r;
  char *pcode1;
  char *pcode2;
  int p;

  ESL_ALLOC(pcode1, (L+1) * sizeof(char)); 
  ESL_ALLOC(pcode2, (L+1) * sizeof(char)); 

  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];
      d = j-i+1;
      /*printf("x: %2d | v: %2d | i: %2d | j: %2d | d: %2d | post[%d][%d][%d]: %f\n", x, v, i, j, d, v, j, d, post[v][j][d]);*/
      /*
       * Only P, L, R states have emissions.
       */
      if (cm->sttype[v] == MP_st) {
	p = Fscore2postcode(post[v][j][d]);
	if(p == 100) { 
	  pcode1[i-1] = pcode1[j-1] = '*';
	  pcode2[i-1] = pcode2[j-1] = '*';
	}
	else {
	  pcode1[i-1] = pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	p = Fscore2postcode(post[v][j][d]);
	if(p == 100) { 
	  pcode1[i-1] = '*';
	  pcode2[i-1] = '*';
	}
	else {
	  pcode1[i-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	p = Fscore2postcode(post[v][j][d]);
	if(p == 100) { 
	  pcode1[j-1] = '*';
	  pcode2[j-1] = '*';
	}
	else {
	  pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == EL_st) /*special case*/ {
	for(r = (i-1); r <= (j-1); r++)
	  {
	    d = j - (r+1) + 1;
	    p = Fscore2postcode(post[v][j][d]);
	    if(p == 100) { 
	      pcode1[r] = '*';
	      pcode2[r] = '*';
	    }
	    else {
	      pcode1[r] = '0' + (char) (p / 10);
	      pcode2[r] = '0' + (char) (p % 10);
	    }
	    /*printf("r: %d | post[%d][%d][%d]: %f | sc: %c\n", r, v, j, d, post[v][j][d], postcode[r]);*/
	  }
      }
    }
  pcode1[L] = '\0';
  pcode2[L] = '\0';

  *ret_pcode1 = pcode1;
  *ret_pcode2 = pcode2;
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
  return; /* never reached */
}

void
CMPostalCodeHB(CM_t *cm, int L, CM_HB_MX *post_mx, Parsetree_t *tr, char **ret_pcode1, char **ret_pcode2)
{
  int status;
  int x, v, i, j, d, r, p;
  char *pcode1;
  char *pcode2;
  int jp_v, dp_v;

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int    **hdmin = cp9b->hdmin;
  /* the DP matrix */
  float ***post  = post_mx->dp; /* pointer to the post DP matrix */

  ESL_ALLOC(pcode1, (L+1) * sizeof(char)); 
  ESL_ALLOC(pcode2, (L+1) * sizeof(char)); 

  for (x = 0; x < tr->n; x++) {
    v = tr->state[x];
    i = tr->emitl[x];
    j = tr->emitr[x];
    d = j-i+1;

    /* Only P, L, R, and EL states have emissions. */
    if(v != cm->M) { 
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    }
    else {
      jp_v = j;
      dp_v = d;
    }      

    if (cm->sttype[v] == MP_st) {
	p = Fscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[i-1] = pcode1[j-1] = '*';
	  pcode2[i-1] = pcode2[j-1] = '*';
	}
	else {
	  pcode1[i-1] = pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	p = Fscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[i-1] = '*';
	  pcode2[i-1] = '*';
	}
	else {
	  pcode1[i-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	p = Fscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[j-1] = '*';
	  pcode2[j-1] = '*';
	}
	else {
	  pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == EL_st) /*special case*/ {
	for(r = (i-1); r <= (j-1); r++) {
	  d = j - (r+1) + 1;
	  p = Fscore2postcode(post[v][j][d]);
	    if(p == 100) { 
	      pcode1[r] = '*';
	      pcode2[r] = '*';
	    }
	    else {
	      pcode1[r] = '0' + (char) (p / 10);
	      pcode2[r] = '0' + (char) (p % 10);
	    }
	    /*printf("r: %d | post[%d][%d][%d]: %f | sc: %c\n", r, v, j, d, post[v][j][d], postcode[r]);*/
	  }
      }
    }
  pcode1[L] = '\0';
  pcode2[L] = '\0';
  *ret_pcode1 = pcode1;
  *ret_pcode2 = pcode2;
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
  return; /* never reached */
}


char  **
alloc_jdbanded_vjd_yshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  int    status;
  char **a;
  int    jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;

  if(j < jmin || i > jmax)
    cm_Fail("ERROR called alloc_jdbanded_vjd_yshadow_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);

  ESL_ALLOC(a, sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  for (jp = (jlast-jfirst+1); jp <= L;     jp++) a[jp]     = NULL;

  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jp: %d | max : %d\n", jp, (jlast)); */
      ESL_DASSERT2(hdmax[jp-jmin] <= (jp+1))
      /* Based on my current understanding the above line should never be false, if it is means there's a valid d
       * in the hd band that is invalid because its > j. I think I check, or ensure, that this
       * doesn't happen when I'm constructing the d bands.
       */
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*printf("\tallocated a[%d]\n", (jp-jfirst));*/
      ESL_ALLOC(a[jp-jfirst], sizeof(char) * bw);
    }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
int** 
alloc_jdbanded_vjd_kshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  int   status;
  int **a;
  int   jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;

  if(j < jmin || i > jmax)
    cm_Fail("ERROR called alloc_jdbanded_vjd_kshadow_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);

  ESL_ALLOC(a, sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  for (jp = (jlast-jfirst+1); jp <= L;     jp++) a[jp]     = NULL;

  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      ESL_DASSERT1((hdmax[jp-jmin] <= (jp+1)));
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;
      ESL_ALLOC(a[jp-jfirst], sizeof(int) * bw);
    }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_FASTALIGN_BENCHMARK
/* gcc -g -O2 -DHAVE_CONFIG_H -I../easel  -c old_cm_dpalign.c 
 * gcc   -o benchmark-fastalign -g -O2 -I. -L. -I../easel -L../easel -DIMPL_FASTALIGN_BENCHMARK cm_dpalign.c old_cm_dpalign.o -linfernal -leasel -lm
 * mpicc -g -O2 -DHAVE_CONFIG_H -I../easel  -c old_cm_dpalign.c  
 * mpicc -o benchmark-fastalign -g -O2 -I. -L. -I../easel -L../easel -DIMPL_FASTALIGN_BENCHMARK cm_dpalign.c old_cm_dpalign.o -linfernal -leasel -lm
 * icc -g -O3 -static -DHAVE_CONFIG_H -I../easel  -c old_cm_dpalign.c 
 * icc -o benchmark-fastalign -O3 -static -I. -L. -I../easel -L../easel -DIMPL_FASTALIGN_BENCHMARK cm_dpalign.c old_cm_dpalign.o -linfernal -leasel -lm
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
#include "old_funcs.h"		/* old  function declarations (from v0.81)*/
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
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use posterior sums during HMM band calculation (widens bands)", 0 },
  { "--dlev",    eslARG_INT,    "0",   NULL, "0<=n<=3",NULL,NULL,NULL, "set verbosity of debugging print statements to <n>", 0 },
  { "--hmmcheck",eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "check that HMM posteriors are correctly calc'ed", 0 },
  { "--cmcheck", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "check that HMM posteriors are correctly calc'ed", 0 },
  { "--optacc",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute optimal accuracy HMM banded alignment alg", 0 },
  { "--post",   eslARG_NONE,    FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute fast float HMM banded Inside/Outside alignment algs", 0 },
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
  CM_HB_MX      *fout_mx;
  char           errbuf[cmERRBUFSIZE];

  int             ***oialpha;    
  float           ***ofalpha;    

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  /* determine sequence length */
  if(esl_opt_IsDefault(go, "-L")) L = cm->clen;      
  else                            L = esl_opt_GetInteger(go, "-L");

  /* configure CM for HMM banded alignment */
  cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  cm->align_opts  |= CM_ALIGN_TIME;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--sums")) cm->align_opts |= CM_ALIGN_SUMS;

  if(esl_opt_GetBoolean(go, "-l")) { 
    cm->config_opts  |= CM_CONFIG_LOCAL;
    cm->config_opts  |= CM_CONFIG_HMMLOCAL;
    cm->config_opts  |= CM_CONFIG_HMMEL;
  }
  if(esl_opt_GetBoolean(go, "--hmmcheck")) cm->align_opts |= CM_ALIGN_CHECKFB;

  ConfigCM(cm, NULL, NULL);

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
  if(esl_opt_GetBoolean(go, "--post") || esl_opt_GetBoolean(go, "--optacc")) { 
    fout_mx = cm_hb_mx_Create(cm->M);
  }
  int do_check = esl_opt_GetBoolean(go, "--cmcheck");
  for (i = 0; i < N; i++)
    {
      L = seqs_to_aln->sq[i]->n;

      esl_stopwatch_Start(w);
      if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, seqs_to_aln->sq[i]->dsq, 1, L, cm->cp9b, FALSE, 0)) != eslOK) cm_Fail(errbuf);
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "Exptl Band calc:", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      esl_stopwatch_Start(w);
      if((status = FastAlignHB(cm, errbuf, seqs_to_aln->sq[i]->dsq, L, 1, L, cm->hbmx, FALSE, NULL, NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastAlignHB() CYK:", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /* fpfast = fopen("tempfast", "w");
	 ParsetreeDump(fpfast, fasttr, cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); */

      if(esl_opt_GetBoolean(go, "-o")) {
	esl_stopwatch_Start(w);
	/* sc = CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, &slowtr, cp9b->jmin, */
	sc = CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, NULL, cm->cp9b->jmin, 
			    cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, cm->cp9b->safe_hdmin, cm->cp9b->safe_hdmax);
	printf("%4d %-30s %10.4f bits ", (i+1), "CYKInside_b_jd() SLOW:", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*fpslow = fopen("tempslow", "w");
	  ParsetreeDump(fpslow, slowtr, cm, seqs_to_aln->sq[i]->dsq, NULL, NULL);*/
      }

      if(esl_opt_GetBoolean(go, "--post")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	if((status = FastInsideAlignHB(cm, errbuf, seqs_to_aln->sq[i]->dsq, 1, L, cm->hbmx, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastInsideAlignHB():", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	if((status = FastOutsideAlignHB(cm, errbuf, seqs_to_aln->sq[i]->dsq, 1, L, fout_mx, cm->hbmx, do_check, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastOutsideAlignHB():", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      if(esl_opt_GetBoolean(go, "--optacc")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	if((status = FastAlignHB(cm, errbuf, seqs_to_aln->sq[i]->dsq, L, 1, L, cm->hbmx, TRUE,     fout_mx, NULL,   NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastAlignHB() OA:", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      /* do old int Inside/Outside if requested */
      if(esl_opt_GetBoolean(go, "--oipost")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = IInside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			     TRUE,	    /* save full alpha so we can run outside */
			     NULL, &oialpha, /* fill alpha, and return it, needed for IOutside() */
			     NULL, NULL,    /* manage your own deckpool, I don't want it */
			     esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			     cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "IInside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = IOutside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			      TRUE,	        /* save full beta */
			      NULL, NULL,	/* manage your own matrix, I don't want it */
			      NULL, NULL,	/* manage your own deckpool, I don't want it */
			      esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			      oialpha, NULL,  /* alpha matrix from IInside(), and don't save it for CMPosterior*/
			      do_check,      /* TRUE to check Outside probs agree with Inside */
			      cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "IOutside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      /* do old float Inside/Outside if requested */
      if(esl_opt_GetBoolean(go, "--ofpost")) {
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FInside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			     TRUE,	    /* save full alpha so we can run outside */
			     NULL, &ofalpha, /* fill alpha, and return it, needed for IOutside() */
			     NULL, NULL,    /* manage your own deckpool, I don't want it */
			     esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			     cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "FInside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	sc = FOutside_b_jd_me(cm, seqs_to_aln->sq[i]->dsq, 1, L,
			      TRUE,	        /* save full beta */
			      NULL, NULL,	/* manage your own matrix, I don't want it */
			      NULL, NULL,	/* manage your own deckpool, I don't want it */
			      esl_opt_GetBoolean(go, "-l"), /* TRUE to allow local begins */
			      ofalpha, NULL,  /* alpha matrix from FInside(), and don't save it */
			      do_check,      /* TRUE to check Outside probs agree with Inside */
			      cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax); /* j and d bands */
	printf("%4d %-30s %10.4f bits ", (i+1), "FOutside_b_jd_me(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }
      printf("\n");
    }

  if(esl_opt_GetBoolean(go, "--post") || esl_opt_GetBoolean(go, "--optacc")) { 
    cm_hb_mx_Destroy(fout_mx);
  }
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  FreeSeqsToAln(seqs_to_aln);

  return 0;

 ERROR:
  cm_Fail("memory allocation error");
  return 0; /* NEVERREACHED */
}
#endif /*IMPL_FASTALIGN_BENCHMARK*/

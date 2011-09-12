/* cm_dpalign_trunc.c
 * DP functions for truncated CYK and optimal accuracy alignment of a target
 * sequence to a CM.
 *
 * Many of the functions here were based on analogous ones for
 * standard (non-truncated) CYK/Inside/Outside alignment in cm_dpalign.c. 
 * Note that the 'fast_'/'Fast' prefixes were dropped.
 *
 * standard version       truncated version
 * (cm_dpalign.c)         (this file)         description
 * ----------------       -----------------   -----------
 * fast_cyk_align()                           (Tr)CYK non-banded (memory intensive)
 * fast_cyk_align_hb()    tr_cyk_align_hb()   (Tr)CYK HMM banded (memory efficient)
 * fast_alignT()                            
 * FastAlign 
 * FastAlignHB()          TrAlignHB()         (Tr)CYK/optacc alignment wrapper HMM banded
 * FastInsideAlignHB      TrInsideAlignHB()   (Tr)Inside HMM banded
 * 
 * EPN, Wed Sep  7 12:13:00 2011
 *
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "p7_config.h"
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

#define DEBUG1 0

static int tr_cyk_align_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_HB_MX *mx,
			   CM_TR_HB_SHADOW_MX *shmx, int *ret_b, float *ret_bsc, char *ret_bmode, float *ret_sc);
static int tr_cyk_align   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_MX *mx,
			   CM_TR_SHADOW_MX *shmx, int *ret_b, float *ret_bsc, char *ret_bmode, float *ret_sc);
static int tr_alignT_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, Parsetree_t *tr, 
			int i0, int j0, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, 
			int do_optacc, CM_TR_HB_MX *post_mx, float size_limit, float *ret_sc);
static int tr_alignT   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, Parsetree_t *tr, 
			int i0, int j0, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
			int do_optacc, CM_TR_MX *post_mx, float size_limit, float *ret_sc);


/* Function: tr_cyk_align_hb()
 * based on fast_cyk_align_hb) which was ...
 * baed on inside_b_me() which was ...
 * based on inside()
 *
 * Date:     EPN, Wed Sep  7 12:13:43 2011
 *           SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose: Run the inside phase of a trCYK alignment using bands in
 *           the j and d dimensions of the DP matrix. Bands were
 *           obtained from an HMM Forward-Backward parse of the target
 *           sequence. Uses float log odds scores.
 *
 *           A CM_TR_HB_MX DP matrix must be passed in. Only cells
 *           valid within the bands given in the CP9Bands_t <cm->cp9b>
 *           will be valid.
 *
 *           The optimal alignment may use a truncation begin
 *           transition, 0->b and we have to be able to trace that
 *           back.  We return a valid b (the optimal 0->b choice) and
 *           bsc (the score if 0->b is used).  If a truncated begin is
 *           part of the optimal parse tree, the optimal alignment
 *           score returned will be bsc and Jyshad[0][L][L] will be
 *           USE_TRUNC_BEGIN (actually, the equivalent in the banded
 *           matrix), and the optimal marginal alignment mode will be
 *           in bmode, telling tr_hb_alignT() to check state
 *           b's bmode matrix and start with a truncation 0->b
 *           entry transition.
 *
 *           Note, when this function was written based on
 *           fast_cyk_align_hb(), it was made less flexible: only full
 *           CM alignment 0..M-1 is allowed. I did that b/c using
 *           subgraphs was a relic from a D&C implementation of
 *           non-banded CYK, which has not been ported to the HMM
 *           banded strategy yet. If that's ever done, it might be
 *           worth trying to turn this function back into one that
 *           works on subgraphs.

 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitaized sequence [i0..j0]   
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           shmx      - the HMM banded shadow matrix to fill in, only cells within bands are valid
 *           ret_b     - best local begin state
 *           ret_bsc   - score for using ret_b
 *           ret_bmode - marginal alignment mode that gives us bsc score for using ret_b
 *           mx        - the dp matrix, only cells within bands in cm->cp9b will 
 *                       be valid. 
 *           ret_sc    - score of optimal, CYK parsetree 
 *                       
 * Returns: <ret_sc>, <ret_b>, <ret_bsc>, see 'Args'
 * 
 * Throws:  <eslOK> on success.
 *          <eslERANGE> if required CM_HB_MX size exceeds <size_limit>, in
 *                      this case, alignment has been aborted, ret_* variables are not valid
 */
int
tr_cyk_align_hb(CM_t *cm, char *errbuf,  ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_HB_MX *mx, 
		     CM_TR_HB_SHADOW_MX *shmx, int *ret_b, float *ret_bsc, char *ret_bmode, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc, Lsc, Rsc;/* temporary scores */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */
  int      bmode = TRMODE_J; /* mode TRMODE_J, TRMODE_L, TRMODE_R, or TRMODE_T truncation mode for obtaining bsc */
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
  int      dp_v, dp_y, dp_z;   /* d index for state v/y/z in alpha */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp;                 /* ESL_MAX(d-sd, 0) */
  int      dp_y_sd;            /* dp_y - sd */
  int      dp_y_sdr;           /* dp_y - sdr */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Wp;                 /* W oalso changes depending on state */
  float    tsc;                /* a transition score */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  /* variables related to truncated alignment (not in fast_cyk_align_hb() */
  float    trunc_penalty = 0.; /* penalty in bits for a truncated hit */
  int      do_J_v, do_J_y, do_J_z; /* is J matrix valid for state v, y, z? */
  int      do_L_v, do_L_y, do_L_z; /* is L matrix valid for state v, y, z? */
  int      do_R_v, do_R_y, do_R_z; /* is R matrix valid for state v, y, z? */
  int      do_T_v, do_T_y, do_T_z; /* is T matrix valid for state v, y, z? */

  /* Contract check */
  if (dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "tr_cyk_align_hb(), dsq is NULL.\n");
  if (mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "tr_cyk_align_hb(), mx is NULL.\n");
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "tr_cyk_align_hb(), cm->cp9b is NULL.\n");

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  /* the DP matrix */
  float ***Jalpha  = mx->Jdp; /* pointer to the Jalpha DP matrix */
  float ***Lalpha  = mx->Ldp; /* pointer to the Lalpha DP matrix */
  float ***Ralpha  = mx->Rdp; /* pointer to the Ralpha DP matrix */
  float ***Talpha  = mx->Tdp; /* pointer to the Talpha DP matrix */

  char  ***Jyshadow = shmx->Jyshadow; /* pointer to the Jyshadow matrix */
  char  ***Lyshadow = shmx->Lyshadow; /* pointer to the Lyshadow matrix */
  char  ***Ryshadow = shmx->Ryshadow; /* pointer to the Ryshadow matrix */
  int   ***Jkshadow = shmx->Jkshadow; /* pointer to the Jkshadow matrix */
  int   ***Lkshadow = shmx->Lkshadow; /* pointer to the Lkshadow matrix */
  int   ***Rkshadow = shmx->Rkshadow; /* pointer to the Rkshadow matrix */
  int   ***Tkshadow = shmx->Tkshadow; /* pointer to the Tkshadow matrix */
  char  ***Lkmode   = shmx->Lkmode;   /* pointer to the Lkmode matrix */
  char  ***Rkmode   = shmx->Rkmode;   /* pointer to the Rkmode matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_hb_mx_GrowTo       (cm,   mx, errbuf, cp9b, W, size_limit)) != eslOK) return status;
  if((status = cm_tr_hb_shadow_mx_GrowTo(cm, shmx, errbuf, cp9b, W, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  ESL_STOPWATCH *w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 
  if(shmx->Jy_ncells_valid > 0) for(i = 0; i < shmx->Jy_ncells_valid; i++) shmx->Jyshadow_mem[i] = USED_EL;
  if(shmx->Ly_ncells_valid > 0) for(i = 0; i < shmx->Ly_ncells_valid; i++) shmx->Lyshadow_mem[i] = USED_TRUNC_END;
  if(shmx->Ry_ncells_valid > 0) for(i = 0; i < shmx->Ry_ncells_valid; i++) shmx->Ryshadow_mem[i] = USED_TRUNC_END;
  if(shmx->Jk_ncells_valid > 0) esl_vec_ISet(shmx->Jkshadow_mem, shmx->Jk_ncells_valid, USED_EL);
  if(shmx->Lk_ncells_valid > 0) esl_vec_ISet(shmx->Lkshadow_mem, shmx->Lk_ncells_valid, USED_TRUNC_END);
  if(shmx->Rk_ncells_valid > 0) esl_vec_ISet(shmx->Rkshadow_mem, shmx->Rk_ncells_valid, USED_TRUNC_END);
  if(shmx->Tk_ncells_valid > 0) esl_vec_ISet(shmx->Tkshadow_mem, shmx->Tk_ncells_valid, USED_TRUNC_END);
  if(shmx->Lk_ncells_valid > 0) for(i = 0; i < shmx->Lk_ncells_valid; i++) shmx->Lkmode_mem[i] = TRMODE_J;
  if(shmx->Rk_ncells_valid > 0) for(i = 0; i < shmx->Rk_ncells_valid; i++) shmx->Rkmode_mem[i] = TRMODE_J;
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " Matrix init CPU time: ");

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v   = cm->oesc[v];  /* emission scores for state v */
    float const *tsc_v   = cm->tsc[v];   /* transition scores for state v */
    float const *lmesc_v = cm->lmesc[v]; /* marginal left  emission scores for state v */
    float const *rmesc_v = cm->rmesc[v]; /* marginal right emission scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];
    do_J_v = cp9b->do_J[v];
    do_L_v = cp9b->do_L[v];
    do_R_v = cp9b->do_R[v];
    do_T_v = cp9b->do_T[v];

    /* re-initialize the J deck if we can do a local end from v */
    if(do_J_v) { 
      if(NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	    dp = ESL_MAX(d-sd, 0);
	    Jalpha[v][jp_v][dp_v] = el_scA[dp] + cm->endsc[v];
	    /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	  }
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	if(do_J_v) Jalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	if(do_L_v) Lalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	if(do_R_v) Ralpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IL states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1] 
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr valid for state y? */
	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	  i    = j - d + 1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */

	  /* We need to treat R differently from and J and L here, by
	   * doing separate 'for (yoffset...' loops for J and R
	   * because we have to fully calculate Jalpha[v][jp_v][dp_v])
	   * before we can start to calculate Ralpha[v][jp_v][dp_v].
	   */
	  /* Handle J and L first */
	  if(do_J_v || do_L_v) { 
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_J_y = cp9b->do_J[y];
	      do_L_y = cp9b->do_L[y];
	      if(do_J_y || do_L_y) { 
		jp_y_sdr = j - jmin[y] - sdr;
		
		if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		  if((do_J_v && do_J_y) && 
		     ((sc = Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > Jalpha[v][jp_v][dp_v])) {
		    Jalpha[v][jp_v][dp_v]   = sc;
		    Jyshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		  }
		  if((do_L_v && do_L_y) && 
		     ((sc = Lalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > Lalpha[v][jp_v][dp_v])) {
		    Lalpha[v][jp_v][dp_v]   = sc;
		    Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_L_OFFSET;
		  }
		}
	      }
	    }
	    if(do_J_v) { 
	      Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]];
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    if(do_L_v) { 
	      if(d >= 2) { 
		Lalpha[v][jp_v][dp_v] += esc_v[dsq[i]];
	      }
	      else { 
		Lalpha[v][jp_v][dp_v]   = esc_v[dsq[i]];
		Lyshadow[v][jp_v][dp_v] = USED_TRUNC_END;
	      }
	      Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    i--;
	  }

	  if(do_R_v) { 
	    /* Handle R separately */
	    Rsc = Ralpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_R_y = cp9b->do_R[y];
	      do_J_y = cp9b->do_J[y];
	      if(do_J_y || do_R_y) { 
		jp_y_sdr = j - jmin[y] - sdr;
		
		/* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,L above) */
		if((d) >= hdmin[y][jp_y_sdr] && (d) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y = d - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));

		  if(do_J_y &&
		     ((sc = Jalpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]) > Rsc)) { 
		    Rsc = sc;
		    Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		  }
		  if(do_R_y &&
		     ((sc = Ralpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]) > Rsc)) { 
		    Rsc = sc;
		    Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_R_OFFSET;
		  }
		}
	      }
	    } /* end of for (yvalid_idx = 0... loop */
	    Ralpha[v][jp_v][dp_v] = Rsc; 
	    /* we use Rsc instead of Ralpha cell in above loop because
	     * Ralpha[v][jp_v][dp_v] may be the same cell as
	     * Ralpha[y][jp_y_sdr][dp_y] if we're an IL state
	     */
	  }
	}
      }
    }
    else if(cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IR states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1].
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */

      /* The first MR_st/IR_st 'for (j...' loop is for J and R matrices which use the same set of j values */
      if(do_J_v || do_R_v) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  j_sdr = j - sdr;
	  
	  /* determine which children y we can legally transit to for v, j */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	  
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	    /* We need to treat L differently from and J and R here, by
	     * doing separate 'for (yoffset...' loops for J because we
	     * have to fully calculate Jalpha[v][jp_v][dp_v]) before we
	     * can start to calculate Lalpha[v][jp_v][dp_v].
	     */
	    /* Handle J and R first */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_J_y = cp9b->do_J[y];
	      do_R_y = cp9b->do_R[y];
	      if(do_J_y || do_R_y) { 
		jp_y_sdr = j - jmin[y] - sdr;
		
		if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));

		  if((do_J_v && do_J_y) && 
		     ((sc = Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > Jalpha[v][jp_v][dp_v])) {
		    Jalpha[v][jp_v][dp_v]   = sc;
		    Jyshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		  }
		  if((do_R_v && do_R_y) && 
		     ((sc = Ralpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > Ralpha[v][jp_v][dp_v])) {
		    Ralpha[v][jp_v][dp_v]   = sc;
		    Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_R_OFFSET;
		  }
		}
	      }
	    }
	    if(do_J_v) { 
	      Jalpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    if(do_R_v) { 
	      if(d >= 2) { 
		Ralpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	      }
	      else { 
		Ralpha[v][jp_v][dp_v]   = esc_v[dsq[j]];
		Ryshadow[v][jp_v][dp_v] = USED_TRUNC_END;
	      }		
	      Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	  }
	}
      }

      if(do_L_v) { 
	/* The second MR_st/IR_st 'for (j...' loop is for the L matrix which use a different set of j values */
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  
	  /* determine which children y we can legally transit to for v, j */
	  /* we use 'j' and not 'j_sdr' here for the L matrix, differently from J and R matrices above */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j) >= jmin[y] && ((j) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j is valid for state y? */
	  
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	    Lsc = Lalpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_L_y = cp9b->do_L[y];
	      do_J_y = cp9b->do_J[y];
	      if(do_L_y || do_J_y) { 
		/* we use 'jp_y=j-min[y]' here, not 'jp_y_sdr=j-jmin[y]-sdr' (which we used in the corresponding loop for J,R above) */
		jp_y = j - jmin[y];
	      
		/* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,R above) */
		if((d) >= hdmin[y][jp_y] && (d) <= hdmax[y][jp_y]) { /* make sure d is valid for this v, j and y */
		  dp_y = d - hdmin[y][jp_y];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v] - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y] - hdmin[y][jp_y])));

		  if(do_J_y &&
		     (sc = Jalpha[y][jp_y][dp_y] + tsc_v[yoffset]) > Lsc) { 
		    Lsc = sc;
		    Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		  }
		  if(do_L_y &&
		     (sc = Lalpha[y][jp_y][dp_y] + tsc_v[yoffset]) > Lsc) { 
		    Lsc = sc;
		    Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_L_OFFSET;
		  }
		}
	      }
	    } /* end of for (yvalid_idx = 0... loop */
	    Lalpha[v][jp_v][dp_v] = Lsc; 
	    /* we use Lsc instead of Lalpha cell in above loop because
	     * Lalpha[v][jp_v][dp_v] may be the same cell as
	     * Lalpha[y][jp_y_sdr][dp_y] if we're an IR state
	     */
	  }
	}
      }
    }
    else if(cm->sttype[v] == MP_st) { 
      /* MP states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	do_J_y = cp9b->do_J[y];
	do_L_y = cp9b->do_L[y];
	do_R_y = cp9b->do_R[y];
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];

	/* The first MP_st 'for (jp_v...' loop is for J and R matrices which use the same set of j values */
	/* j must satisfy:
	 * j >= jmin[v]
	 * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	 * j <= jmax[v]
	 * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	 * this reduces to two ESL_MAX calls
	 */
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	/* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	
	if((do_J_v && do_J_y) || (do_R_v && (do_J_y || do_R_y))) { 
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++, jp_y++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	    
	    if(do_J_v && do_J_y) { 
	      /* J matrix: */
	      /* d must satisfy:
	       * d >= hdmin[v][jp_v]
	       * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	       * d <= hdmax[v][jp_v]
	       * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	       * this reduces to two ESL_MAX calls
	       */
	      dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	      dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	      dpn       = dn - hdmin[v][jp_v];
	      dpx       = dx - hdmin[v][jp_v];
	      dp_y_sd   = dn - hdmin[y][jp_y_sdr] - sd;
	      
	      for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
		ESL_DASSERT1((dp_v      >= 0 && dp_v       <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		ESL_DASSERT1((dp_y_sd   >= 0 && dp_y_sd    <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));

		if((sc = Jalpha[y][jp_y_sdr][dp_y_sd] + tsc) > Jalpha[v][jp_v][dp_v]) { 
		  Jalpha[v][jp_v][dp_v]   = sc;
		  Jyshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		}
	      }
	    }
	    
	    if(do_R_v && (do_R_y || do_J_y)) { 
	      /* R matrix: */
	      /* d must satisfy:
	       * d >= hdmin[v][jp_v]
	       * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	       * d <= hdmax[v][jp_v]
	       * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	       * this reduces to two ESL_MAX calls
	       */
	      dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sdr);
	      dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sdr);
	      dpn       = dn - hdmin[v][jp_v];
	      dpx       = dx - hdmin[v][jp_v];
	      dp_y_sdr  = dn - hdmin[y][jp_y_sdr] - sdr;
	      /* for {L,R}alpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	      
	      for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
		/* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
		ESL_DASSERT1((dp_y_sdr  >= 0 && dp_y_sdr   <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		if(do_J_y && 
		   ((sc = Jalpha[y][jp_y_sdr][dp_y_sdr] + tsc) > Ralpha[v][jp_v][dp_v])) { 
		  Ralpha[v][jp_v][dp_v]   = sc;
		  Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		}
		if(do_R_y && 
		   ((sc = Ralpha[y][jp_y_sdr][dp_y_sdr] + tsc) > Ralpha[v][jp_v][dp_v])) { 
		  Ralpha[v][jp_v][dp_v]   = sc;
		  Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_R_OFFSET;
		}
	      }
	    }
	  }
	}

	if(do_L_v && (do_L_y || do_J_y)) { 
	  /* The second MP_st 'for (jp_v...' loop is for L matrix, which uses a different set of j values from J and R */
	  /* j must satisfy:
	   * j >= jmin[v]
	   * j >= jmin[y] (follows from (j >= jmin[y]))
	   * j <= jmax[v]
	   * j <= jmax[y] (follows from (j <= jmax[y]))
	   * this reduces to two ESL_MAX calls
	   */
	  jn = ESL_MAX(jmin[v], jmin[y]);
	  jx = ESL_MIN(jmax[v], jmax[y]);
	  jpn = jn - jmin[v];
	  jpx = jx - jmin[v];
	  jp_y = jn - jmin[y];
	  /* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	  
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y     >= 0 && jp_y     <= (jmax[y]-jmin[y])));
	    
	    /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	   * this reduces to two ESL_MAX calls
	   */
	    dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y] + sdr);
	    dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y] + sdr);
	    dpn       = dn - hdmin[v][jp_v];
	    dpx       = dx - hdmin[v][jp_v];
	    dp_y_sdr  = dn - hdmin[y][jp_y] - sdr;
	    /* for Lalpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	    
	    for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
	      /* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
	      ESL_DASSERT1((dp_y_sdr >= 0 && dp_y_sdr  <= (hdmax[y][jp_y]     - hdmin[y][jp_y])));
	      if(do_J_y && 
		 ((sc = Jalpha[y][jp_y][dp_y_sdr] + tsc) > Lalpha[v][jp_v][dp_v])) { 
		Lalpha[v][jp_v][dp_v]  = sc;
		Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
	      }		
	      if(do_L_y && 
		 ((sc = Lalpha[y][jp_y][dp_y_sdr] + tsc) > Lalpha[v][jp_v][dp_v])) { 
		Lalpha[v][jp_v][dp_v]  = sc;
		Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_L_OFFSET;
	      }		
	    }
	  }
	}
      }
      /* add in emission score */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	i     = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v], dp_v = 0; d <= hdmax[v][jp_v]; d++, dp_v++) 
	  {
	    if(d >= 2) { 
	      if(do_J_v) Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	      if(do_L_v) Lalpha[v][jp_v][dp_v] += lmesc_v[dsq[i]];
	      if(do_R_v) Ralpha[v][jp_v][dp_v] += rmesc_v[dsq[j]];
	    }
	    else { 
	      if(do_J_v) { Jalpha[v][jp_v][dp_v] = IMPOSSIBLE; }
	      if(do_L_v) { Lalpha[v][jp_v][dp_v] = lmesc_v[dsq[i]]; Lyshadow[v][jp_v][dp_v] = USED_TRUNC_END; }
	      if(do_R_v) { Ralpha[v][jp_v][dp_v] = rmesc_v[dsq[j]]; Ryshadow[v][jp_v][dp_v] = USED_TRUNC_END; }
	    }
	    i--;
	  }
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	  if(do_J_v) Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  if(do_L_v) Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  if(do_R_v) Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is D or S */
      /* D, S states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	do_J_y = cp9b->do_J[y];
	do_L_y = cp9b->do_L[y];
	do_R_y = cp9b->do_R[y];
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	if((do_J_v && do_J_y) || (do_L_v && do_L_y) || (do_R_v && do_R_y)) { 
	  /* j must satisfy:
	   * j >= jmin[v]
	   * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	   * j <= jmax[v]
	   * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	   * this reduces to two ESL_MAX calls
	   */
	  jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	  jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	  jpn = jn - jmin[v];
	  jpx = jx - jmin[v];
	  jp_y_sdr = jn - jmin[y] - sdr;
	  
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	    
	    /* d must satisfy:
	     * d >= hdmin[v][jp_v]
	     * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	     * d <= hdmax[v][jp_v]
	     * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	     * this reduces to two ESL_MAX calls
	     */
	    dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	    dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	    dpn     = dn - hdmin[v][jp_v];
	    dpx     = dx - hdmin[v][jp_v];
	    dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	    
	    for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));

	      if((do_J_v && do_J_y) && 
		 ((sc = Jalpha[y][jp_y_sdr][dp_y_sd] + tsc) > Jalpha[v][jp_v][dp_v])) { 
		Jalpha[v][jp_v][dp_v]  = sc;
		Jyshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
	      }
	      if((do_L_v && do_L_y) && 
		 ((sc = Lalpha[y][jp_y_sdr][dp_y_sd] + tsc) > Lalpha[v][jp_v][dp_v])) { 
		Lalpha[v][jp_v][dp_v]  = sc;
		Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_L_OFFSET;
	      }
	      if((do_R_v && do_R_y) && 
		 ((sc = Ralpha[y][jp_y_sdr][dp_y_sd] + tsc) > Ralpha[v][jp_v][dp_v])) { 
		Ralpha[v][jp_v][dp_v]  = sc;
		Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_R_OFFSET;
	      }
	      /* an easy to overlook case: if d == 0, ensure L and R values are IMPOSSIBLE */
	      if(dp_v == dpn && dn == 0) { /* d is 0 */
		if(do_L_v) Lalpha[v][jp_v][dp_v] = IMPOSSIBLE;
		if(do_R_v) Ralpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }		
	    }
	  }
	}
      }
      /* no emission score to add */
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */

      do_J_y = cp9b->do_J[y];
      do_L_y = cp9b->do_L[y];
      do_R_y = cp9b->do_R[y];
      do_T_y = cp9b->do_T[y]; /* will be FALSE, y is not a B_st */

      do_J_z = cp9b->do_J[z];
      do_L_z = cp9b->do_L[z];
      do_R_z = cp9b->do_R[z];
      do_T_z = cp9b->do_T[z]; /* will be FALSE, z is not a B_st */
      
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
	i = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++, i--) {
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	  /* Find the first k value that implies a valid cell in the {J,L,R} matrix y and z decks.
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
	   *
	   * To update a cell in the T matrix with a sum of an R matrix value for y
	   * and a L matrix value for z, there are 2 additional inequalities to satisfy:
	   * (7) k != 0
	   * (8) k != d
	   * We ensure 7 and 8 in the loop below.
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
	      if((do_J_v && do_J_y && do_J_z) && 
		 ((sc = Jalpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]) > Jalpha[v][jp_v][dp_v])) { 
		Jalpha[v][jp_v][dp_v]   = sc;
		Jkshadow[v][jp_v][dp_v] = k;
	      }
	      if((do_L_v && do_J_y && do_L_z) && 
		 ((sc = Jalpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]) > Lalpha[v][jp_v][dp_v])) { 
		Lalpha[v][jp_v][dp_v]   = sc;
		Lkshadow[v][jp_v][dp_v] = k;
		Lkmode[v][jp_v][dp_v]   = TRMODE_J;
	      }
	      if((do_R_v && do_R_y && do_J_z) && 
		 ((sc = Ralpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]) > Ralpha[v][jp_v][dp_v])) { 
		Ralpha[v][jp_v][dp_v]   = sc;
		Rkshadow[v][jp_v][dp_v] = k;
		Rkmode[v][jp_v][dp_v]   = TRMODE_J;
	      }
	      if(k != 0 && k != d) {
		if((do_T_v && do_R_y && do_L_z) && 
		   ((sc = Ralpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]) > Talpha[v][jp_v][dp_v])) { 
		  Talpha[v][jp_v][dp_v]   = sc;
		  Tkshadow[v][jp_v][dp_v] = k;
		}
	      }
	    }
	  }
	}
      }

      /* two additional special cases in trCYK (these are not in standard CYK).
       * we do these in their own for(j.. { for(d.. { } } loops b/c one 
       * is independent of z, the other of y, unlike the above loop which is dependent 
       * on both.
       */
      if(do_L_v && (do_J_y || do_L_y)) { 
	jn = (jmin[v] > jmin[y]) ? jmin[v] : jmin[y];
	jx = (jmax[v] < jmax[y]) ? jmax[v] : jmax[y];
	for (j = jn; j <= jx; j++) { 
	  jp_v = j - jmin[v];
	  jp_y = j - jmin[y];
	  ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	  ESL_DASSERT1((j >= jmin[y] && j <= jmax[y]));
	  dn = (hdmin[v][jp_v] > hdmin[y][jp_y]) ? hdmin[v][jp_v] : hdmin[y][jp_y];
	  dx = (hdmax[v][jp_v] < hdmax[y][jp_y]) ? hdmax[v][jp_v] : hdmax[y][jp_y];
	  for(d = dn; d <= dx; d++) { 
	    dp_v = d - hdmin[v][jp_v];
	    dp_y = d - hdmin[y][jp_y];
	    ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	    ESL_DASSERT1((d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y]));
	    if(do_J_y &&
	       ((sc = Jalpha[y][jp_y][dp_y]) > Lalpha[v][jp_v][dp_v])) { 
	      Lalpha[v][jp_v][dp_v]   = sc;
	      Lkshadow[v][jp_v][dp_v] = 0; /* k == 0 for this case, full sequence is on left */
	      Lkmode[v][jp_v][dp_v]   = TRMODE_J;
	      /* consider making a different mode here, to let the traceback know that right child emits 0 residues,
	       * this should then effect the alignment display, no? it is a different case from the
	       * >0 residues from right child TRMODE_J case for Lalpha checked for in (for k) loop above.
	       */
	    }
	    if(do_L_y &&
	       ((sc = Lalpha[y][jp_y][dp_y]) > Lalpha[v][jp_v][dp_v])) { 
	      Lalpha[v][jp_v][dp_v]   = sc;
	      Lkshadow[v][jp_v][dp_v] = 0; /* k == 0 for this case, full sequence is on left */
	      Lkmode[v][jp_v][dp_v]   = TRMODE_L;
	      /* consider making a different mode here, to let the traceback know that right child emits 0 residues,
	       * this should then effect the alignment display, no? it is a different case from the
	       * >0 residues from right child TRMODE_L case for Lalpha checked for in (for k) loop above.
	       */
	    }
	  }
	}
      }
      if(do_R_v && (do_J_z || do_R_z)) { 
	jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
	for (j = jn; j <= jx; j++) { 
	  jp_v = j - jmin[v];
	  jp_z = j - jmin[z];
	  ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	  ESL_DASSERT1((j >= jmin[z] && j <= jmax[z]));
	  dn = (hdmin[v][jp_v] > hdmin[z][jp_z]) ? hdmin[v][jp_v] : hdmin[z][jp_z];
	  dx = (hdmax[v][jp_v] < hdmax[z][jp_z]) ? hdmax[v][jp_v] : hdmax[z][jp_z];
	  for(d = dn; d <= dx; d++) { 
	    dp_v = d - hdmin[v][jp_v];
	    dp_z = d - hdmin[z][jp_z];
	    ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	    ESL_DASSERT1((d >= hdmin[z][jp_z] && d <= hdmax[z][jp_z]));
	    if(do_J_z &&
	       ((sc = Jalpha[z][jp_z][dp_z]) > Ralpha[v][jp_v][dp_v])) { 
	      Ralpha[v][jp_v][dp_v]   = sc;
	      Rkshadow[v][jp_v][dp_v] = d; /* k == d in this case, full sequence is on right */
	      Rkmode[v][jp_v][dp_v]   = TRMODE_J;
	      /* consider making a different mode here, to let the traceback know that left child emits 0 residues,
	       * this should then effect the alignment display, no? it is a different case from the
	       * >0 residues from left child TRMODE_J case for Ralpha checked for in (for k) loop above.
	       */
	    }
	    if(do_R_z &&
	       ((sc = Ralpha[z][jp_z][dp_z]) > Ralpha[v][jp_v][dp_v])) { 
	      Ralpha[v][jp_v][dp_v]   = sc;
	      Rkshadow[v][jp_v][dp_v] = d; /* k == d in this case, full sequence is on right */
	      Rkmode[v][jp_v][dp_v]   = TRMODE_R;
	      /* consider making a different mode here, to let the traceback know that left child emits 0 residues,
	       * this should then effect the alignment display, no? it is a different case from the
	       * >0 residues from left child TRMODE_R case for Ralpha checked for in (for k) loop above.
	       */
	    }
	  }
	}
      }
    } /* finished calculating deck v. */
         
    /* The following loops originally access alpha[v][j0][W] but the index W will be
       in different positions due to the bands */
    if(j0 >= jmin[v] && j0 <= jmax[v]) { 
      jp_v = j0 - jmin[v];
      Wp = W - hdmin[v][jp_v];
      if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	/* If we get here alpha[v][jp_v][Wp] is a valid cell
	 * in the banded alpha matrix, corresponding to 
	 * alpha[v][j0][W] in the platonic matrix.
	 * 
	 * Check for truncated alignment getting us to the root.
	 * This is "off-shadow": if/when we trace back, we'll handle this
	 * case separately (and we'll know to do it because we'll immediately
	 * see a USED_TRUNC_BEGIN in the shadow matrix, 
	 * telling us to jump right to state b; see below)
	 */
#if 0
	/* Don't allow local begins, right? */
	/* check for normal local begins */
	if(do_J_v && (cm->flags & CMH_LOCAL_BEGIN) && allow_begin) {
	  if (Jalpha[v][jp_v][Wp] + cm->beginsc[v] > bsc) { 
	    b     = v;
	    bsc   = Jalpha[v][jp_v][Wp] + cm->beginsc[v];
	    bmode = TRMODE_J;
	  }
	}
#endif
	/* check for hit in J matrix (much like a normal local begin) */
	if(do_J_v && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st)) { 
	  if (Jalpha[v][jp_v][Wp] + trunc_penalty > bsc) { 
	    b     = v;
	    bsc   = Jalpha[v][jp_v][Wp] + trunc_penalty;
	    bmode = TRMODE_J;
	  }
	}
	/* check for truncated hit in L matrix */
	if(do_L_v && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st)) { 
	  if (Lalpha[v][jp_v][Wp] + trunc_penalty > bsc) { 
	    b     = v;
	    bsc   = Lalpha[v][jp_v][Wp] + trunc_penalty;
	    bmode = TRMODE_L;
	  }
	}	    
	/* check for truncated hit in R matrix */
	if(do_R_v && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st)) { 
	  if (Ralpha[v][jp_v][Wp] + trunc_penalty > bsc) { 
	    b     = v;
	    bsc   = Ralpha[v][jp_v][Wp] + trunc_penalty;
	    bmode = TRMODE_R;
	  }
	}	    
	/* check for truncated hit in T matrix */
	if(do_T_v) { 
	  if (Talpha[v][jp_v][Wp] + trunc_penalty > bsc) { 
	    b     = v;
	    bsc   = Talpha[v][jp_v][Wp] + trunc_penalty;
	    bmode = TRMODE_T;
	  }
	}	    
      }
      /* Check for whether we need to store best truncated score as
       * the optimal overall score, and if we need to put a flag in
       * the shadow matrix telling tr_alignT() to use the b we
       * return.
       */
      if (v == 0) { 
	assert(do_J_v);
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp   = W - hdmin[v][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    if (bsc > Jalpha[0][jp_v][Wp]) {
	      Jalpha[0][jp_v][Wp]   = bsc;
	      Jyshadow[v][jp_v][Wp] = USED_TRUNC_BEGIN;
	    }
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp1; fp1 = fopen("tmp.ahbmx", "w");   cm_tr_hb_mx_Dump(fp1, mx); fclose(fp1);*/
  /*FILE *fp2; fp2 = fopen("tmp.ahbshmx", "w"); cm_tr_hb_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2);*/
  
  Wp = W - hdmin[0][j0-jmin[0]];
  sc =    Jalpha[0][j0-jmin[0]][Wp]; /* this will be bsc, unless a non-truncated hit rooted at 0 is optimal */

  *ret_b     = b;    
  *ret_bsc   = bsc;  
  *ret_bmode = bmode;  
  *ret_sc    = sc;

  free(el_scA);
  free(yvalidA);

  ESL_DPRINTF1(("tr_cyk_align_hb return sc: %f\n", sc));
  printf("tr_cyk_align_hb return sc: %.4f (bsc: %.4f)\n", sc, bsc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/* Function: tr_cyk_align()
 * based on fast_cyk_align()
 *
 * Date:     EPN, Fri Sep  9 15:35:06 2011
 *
 * Note:     Very similar to inside(), but slightly more efficient.
 *           Identical to fast_cyk_align_hb() but HMM bands are NOT
 *           used.
 *
 * Purpose:  Perform trCYK alignment on a full sequence i0..j0
 *           rooted at state 0. 
 *
 *           Note, when this function was written based on
 *           fast_cyk_align(), it was made less flexible: only full CM
 *           alignment 0..M-1 is allowed. This function is mainly only
 *           useful for testing and debugging, as a non-banded version
 *           of tr_cyk_align_hb().
 *
 * Args:     cm          - the model    [0..M-1]
 *           errbuf      - char buffer for reporting errors
 *           dsq         - the digitaized sequence [i0..j0]   
 *           i0          - first position in subseq to align (1, for whole seq)
 *           j0          - last position in subseq to align (L, for whole seq)
 *           size_limit  - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx          - dp matrix 
 *           shmx        - shadow matrix
 *           ret_b       - best local begin state, or NULL if unwanted
 *           ret_bsc     - score for using ret_b, or NULL if unwanted                        
 *           ret_Jmx     - the (J)oint dp matrix, we'll allocate it here, NULL if not wanted back
 *           ret_Lmx     - the (L)oint dp matrix, we'll allocate it here, NULL if not wanted back
 *           ret_Rmx     - the (R)oint dp matrix, we'll allocate it here, NULL if not wanted back
 *           ret_Tmx     - the (T)oint dp matrix, we'll allocate it here, NULL if not wanted back
 *           ret_Jshadow - if non-NULL, caller needs shadow matrix to traceback
 *           ret_Lshadow - if non-NULL, caller needs shadow matrix to traceback
 *           ret_Rshadow - if non-NULL, caller needs shadow matrix to traceback
 *           ret_Tshadow - if non-NULL, caller needs shadow matrix to traceback
 *           ret_Lmode   - if non-NULL, caller needs shadow matrix to traceback
 *           ret_Rmode   - if non-NULL, caller needs shadow matrix to traceback
 *           ret_sc      - score of optimal, CYK parsetree 
 *                       
 * Returns: <ret_sc>, <ret_b>, <ret_bsc>, <ret_mx>, <ret_shadow>, see 'Args'.
 * 
 * Throws:  <eslOK> on success.
 *          <eslERANGE> if required DP matrix size exceeds passed in <size_limit> 
 *                      alignment has been aborted, ret_* variables are not valid
 */
int
tr_cyk_align(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
	     int *ret_b, float *ret_bsc, char *ret_bmode, float *ret_sc)
{
  int      status;          /* easel status code */
  int      v,y,z;	    /* indices for states  */
  int      j,d,i,k;	    /* indices in sequence dimensions */
  int      dp;              /* temporary d value */
  float    sc;		    /* a temporary variable holding a score */
  int      yoffset;	    /* y=base+offset -- counter in child states that v can transit to */
  int      W;	    	    /* subsequence length */
  int      b;	  	    /* best local begin state */
  float    bsc;   	    /* score for using the best local begin state */
  float   *el_scA;          /* [0..d..W-1] probability of local end emissions of length d */
  int      sd;              /* StateDelta(cm->sttype[v]) */
  int      sdl;             /* StateLeftDelta(cm->sttype[v] */
  int      sdr;             /* StateRightDelta(cm->sttype[v] */
  int      jp;              /* offset j, j = i0-1+jp */
  int      j_sdr;           /* j - sdr */
  int      d_sd;            /* d - sd */
  int      d_sdl;           /* d - sdl */
  int      d_sdr;           /* d - sdr */
  float    tsc;             /* a transition score */
  int      L = j0-i0+1;     /* length of the sequence */
  /* other variables used in truncated version, but not standard version (not in fast_cyk_align()) */
  float Lsc, Rsc;           /* temporary scores */
  float trunc_penalty = 0.; /* penalty in bits for a truncated hit */
  int   bmode = TRMODE_J;   /* mode TRMODE_J, TRMODE_L, TRMODE_R, or TRMODE_T truncation mode for obtaining bsc */

  /* Contract check */
  if(dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_cyk_inside_align(), dsq is NULL.\n");

  /* the DP matrix */
  float ***Jalpha  = mx->Jdp; /* pointer to the Jalpha DP matrix */
  float ***Lalpha  = mx->Ldp; /* pointer to the Lalpha DP matrix */
  float ***Ralpha  = mx->Rdp; /* pointer to the Ralpha DP matrix */
  float ***Talpha  = mx->Tdp; /* pointer to the Talpha DP matrix */

  char  ***Jyshadow = shmx->Jyshadow; /* pointer to the Jyshadow matrix */
  char  ***Lyshadow = shmx->Lyshadow; /* pointer to the Lyshadow matrix */
  char  ***Ryshadow = shmx->Ryshadow; /* pointer to the Ryshadow matrix */
  int   ***Jkshadow = shmx->Jkshadow; /* pointer to the Jkshadow matrix */
  int   ***Lkshadow = shmx->Lkshadow; /* pointer to the Lkshadow matrix */
  int   ***Rkshadow = shmx->Rkshadow; /* pointer to the Rkshadow matrix */
  int   ***Tkshadow = shmx->Tkshadow; /* pointer to the Tkshadow matrix */
  char  ***Lkmode   = shmx->Lkmode;   /* pointer to the Lkmode matrix */
  char  ***Rkmode   = shmx->Rkmode;   /* pointer to the Rkmode matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo       (cm, mx,   errbuf, W, size_limit)) != eslOK) return status;
  if((status = cm_tr_shadow_mx_GrowTo(cm, shmx, errbuf, W, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  ESL_STOPWATCH *w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 
  if(shmx->Jy_ncells_valid > 0) for(i = 0; i < shmx->Jy_ncells_valid; i++) shmx->Jyshadow_mem[i] = USED_EL;
  if(shmx->Ly_ncells_valid > 0) for(i = 0; i < shmx->Ly_ncells_valid; i++) shmx->Lyshadow_mem[i] = USED_TRUNC_END;
  if(shmx->Ry_ncells_valid > 0) for(i = 0; i < shmx->Ry_ncells_valid; i++) shmx->Ryshadow_mem[i] = USED_TRUNC_END;
  if(shmx->Jk_ncells_valid > 0) esl_vec_ISet(shmx->Jkshadow_mem, shmx->Jk_ncells_valid, USED_EL);
  if(shmx->Lk_ncells_valid > 0) esl_vec_ISet(shmx->Lkshadow_mem, shmx->Lk_ncells_valid, USED_TRUNC_END);
  if(shmx->Rk_ncells_valid > 0) esl_vec_ISet(shmx->Rkshadow_mem, shmx->Rk_ncells_valid, USED_TRUNC_END);
  if(shmx->Tk_ncells_valid > 0) esl_vec_ISet(shmx->Tkshadow_mem, shmx->Tk_ncells_valid, USED_TRUNC_END);
  if(shmx->Lk_ncells_valid > 0) for(i = 0; i < shmx->Lk_ncells_valid; i++) shmx->Lkmode_mem[i] = TRMODE_J;
  if(shmx->Rk_ncells_valid > 0) for(i = 0; i < shmx->Rk_ncells_valid; i++) shmx->Rkmode_mem[i] = TRMODE_J;
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " Matrix init CPU time: ");

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    float const *lmesc_v = cm->lmesc[v]; /* marginal left  emission scores for state v */
    float const *rmesc_v = cm->rmesc[v]; /* marginal right emission scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdl  = StateLeftDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);

    /* re-initialize the J deck if we can do a local end from v */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) { 
	  dp = ESL_MAX(d-sd, 0);
	  Jalpha[v][j][d] = el_scA[dp] + cm->endsc[v];
	  /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (jp = 0; jp <= W; jp++) {
	j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
	Jalpha[v][j][0] = 0.;
	Lalpha[v][j][0] = 0.;
	Ralpha[v][j][0] = 0.;
	for (d = 1; d <= jp; d++) { 
	  Jalpha[v][j][d] = IMPOSSIBLE;
	  Lalpha[v][j][d] = IMPOSSIBLE;
	  Ralpha[v][j][d] = IMPOSSIBLE;
	}
      }
    }
    else if(cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
      /* update alpha[v][j][d] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] 
       * We do ML states as well as IL states b/c they follow the same rules, 
       * and we're not worried about efficiency here.
       */
      
      /* In TrCYK: we need to treat R differently from and J and L
       * here, by doing separate 'for (yoffset...' loops for J and R
       * because we have to fully calculate Jalpha[v][j][d]) before we
       * can start to calculate Ralpha[v][j][d].
       */

      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i    = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = Jalpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > Jalpha[v][j][d]) {
	      Jalpha[v][j][d]   = sc; 
	      Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if ((sc = Lalpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > Lalpha[v][j][d]) {
	      Lalpha[v][j][d]   = sc; 
	      Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	    }
	  }
	  Jalpha[v][j][d] += esc_v[dsq[i]];
	  Jalpha[v][j][d]  = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	  if(d >= 2) { 
	    Lalpha[v][j][d] += esc_v[dsq[i]];
	  }
	  else { 
	    Lalpha[v][j][d]   = esc_v[dsq[i]];
	    Lyshadow[v][j][d] = USED_TRUNC_END;
	  }
	  Lalpha[v][j][d] = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	  i--;

	  /* handle R separately */
	  /* note we use 'd', not 'd_sd' (which we used in the corresponding loop for J,L above) */
	  Rsc = Ralpha[v][j][d]; /* this sc will be IMPOSSIBLE */
	  /* impt to use Rsc because Ralpha[v][j][d] is one of the possible child states we transit to! (if IL) */
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = Jalpha[y][j_sdr][d] + tsc_v[yoffset]) > Rsc) { 
	      Rsc = sc; 
	      Ryshadow[v][j][d]= yoffset + TRMODE_J_OFFSET;
	    }
	    if ((sc = Ralpha[y][j_sdr][d] + tsc_v[yoffset]) > Rsc) { 
	      Rsc = sc;
	      Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	    }
	  }
	  Ralpha[v][j][d] = Rsc;
	}
      }
    }
    else if(cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) { 
      /* update alpha[v][j][d] cells, for IR states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1].
       * We do MR states as well as IR states b/c they follow the same rules, 
       * and we're not worried about efficiency here.
       */

      /* In TrCYK: we need to treat L differently from and J and R
       * here, by doing separate 'for (yoffset...' loops for J and R
       * because we have to fully calculate Jalpha[v][j][d]) before we
       * can start to calculate Lalpha[v][j][d].
       */
      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = Jalpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > Jalpha[v][j][d]) {
	      Jalpha[v][j][d]   = sc; 
	      Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if ((sc = Ralpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > Ralpha[v][j][d]) {
	      Ralpha[v][j][d]   = sc; 
	      Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	    }
	  }
	  Jalpha[v][j][d] += esc_v[dsq[j]];
	  Jalpha[v][j][d]  = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	  if(d >= 2) { 
	    Ralpha[v][j][d] += esc_v[dsq[j]];
	  }
	  else { 
	    Ralpha[v][j][d]   = esc_v[dsq[j]];
	    Ryshadow[v][j][d] = USED_TRUNC_END;
	  }
	  Ralpha[v][j][d] = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);

	  /* handle L separately */
	  /* note we use 'j' and 'd', not 'j_sdr' and 'd_sd' (which we used in the corresponding loop for J,R above) */
	  Lsc = Lalpha[v][j][d]; /* this sc will be IMPOSSIBLE */
	  /* impt to use Lsc because Lalpha[v][j][d] is one of the possible child states we transit to! (if IR) */
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = Jalpha[y][j][d] + tsc_v[yoffset]) > Lsc) { 
	      Lalpha[v][j][d]   = Lsc; 
	      Lyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if ((sc = Lalpha[y][j][d] + tsc_v[yoffset]) > Lsc) { 
	      Lalpha[v][j][d]   = Lsc; 
	      Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	    }
	  }
	  Lalpha[v][j][d] = Lsc;
	}
      }
    }
    else if(cm->sttype[v] == MP_st) { 
      /* MP states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];

	for (jp = sdr; jp <= W; jp++) {
	  j = i0-1+jp;
	  j_sdr = j - sdr;

	  for (d = sd; d <= jp; d++) { /* sd == 2 for MP state */
	    d_sd = d-sd;
	    if((sc = Jalpha[y][j_sdr][d_sd] + tsc) > Jalpha[v][j][d]) {
	      Jalpha[v][j][d]   = sc;
	      Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	  }
	  /* note we use 'd_sdl' not 'd_sd' for L, plus minimum d is sdl (1) */
	  for (d = sdl; d <= jp; d++) { /* sdl == 1 for MP state */
	    d_sdl = d-sdl;
	    if((sc = Jalpha[y][j_sdr][d_sdl] + tsc) > Lalpha[v][j][d]) {
	      Lalpha[v][j][d]   = sc;
	      Lyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if((sc = Lalpha[y][j_sdr][d_sdl] + tsc) > Lalpha[v][j][d]) {
	      Lalpha[v][j][d]   = sc;
	      Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	    }
	  }
	  /* note we use 'j' and 'd_sdr' not 'j_sdr' and 'd_sd' for R, plus minimum d is sdr (1) */
	  for (d = sdr; d <= jp; d++) { /* sdr == 1 for MP state */
	    d_sdr = d - sdr;
	    if((sc = Jalpha[y][j][d_sdr] + tsc) > Ralpha[v][j][d]) {
	      Ralpha[v][j][d]   = sc;
	      Ryshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if((sc = Ralpha[y][j][d_sdr] + tsc) > Ralpha[v][j][d]) {
	      Ralpha[v][j][d]   = sc;
	      Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	    }
	  }
	}
      }
      /* add in emission score */
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	i = j;
	Jalpha[v][j][1] = IMPOSSIBLE;
	Lalpha[v][j][1] = lmesc_v[dsq[i]];
	Ralpha[v][j][1] = rmesc_v[dsq[j]];
	i--;
	for (d = 2; d <= jp; d++) {
	  Jalpha[v][j][d] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	  Lalpha[v][j][d] += lmesc_v[dsq[i]];
	  Ralpha[v][j][d] += rmesc_v[dsq[j]];
	  i--;
	}
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 1; d <= jp; d++) {
	  Jalpha[v][j][d] = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	  Lalpha[v][j][d] = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	  Ralpha[v][j][d] = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is D or S */
      /* D, S states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	for (jp = sdr; jp <= W; jp++) {
	  j = i0-1+jp;
	  j_sdr = j - sdr;
	  
	  for (d = sd; d <= jp; d++) {
	    d_sd = d-sd;
	    if((sc = Jalpha[y][j_sdr][d_sd] + tsc) > Jalpha[v][j][d]) {
	      Jalpha[v][j][d]   = sc;
	      Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if((sc = Lalpha[y][j_sdr][d_sd] + tsc) > Lalpha[v][j][d]) {
	      Lalpha[v][j][d]   = sc;
	      Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	    }
	    if((sc = Ralpha[y][j_sdr][d_sd] + tsc) > Ralpha[v][j][d]) {
	      Ralpha[v][j][d]   = sc;
	      Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	    }
	  }
	  /* an easy to overlook case: if d == 0, ensure L and R values are IMPOSSIBLE */
	  Lalpha[v][j][0] = IMPOSSIBLE;
	  Ralpha[v][j][0] = IMPOSSIBLE;
	}
      }
      /* no emission score to add */
    }
    else { /* B_st */
      assert(cm->sttype[v] == B_st);
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */

      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  for (k = 1; k < d; k++) {
	    if((sc = Jalpha[y][j-k][d-k] + Jalpha[z][j][k]) > Jalpha[v][j][d]) { 
	      Jalpha[v][j][d]   = sc;
	      Jkshadow[v][j][d] = k;
	    }
	    if((sc = Jalpha[y][j-k][d-k] + Lalpha[z][j][k]) > Lalpha[v][j][d]) { 
	      Lalpha[v][j][d]   = sc;
	      Lkshadow[v][j][d] = k;
	      Lkmode[v][j][d]   = TRMODE_J;
	    }
	    if((sc = Ralpha[y][j-k][d-k] + Jalpha[z][j][k]) > Ralpha[v][j][d]) { 
	      Ralpha[v][j][d]   = sc;
	      Rkshadow[v][j][d] = k;
	      Rkmode[v][j][d]   = TRMODE_J;
	    }
	    /*if((k != i-1) && (k != j)) {*/
	    if((sc = Ralpha[y][j-k][d-k] + Lalpha[z][j][k]) > Talpha[v][j][d]) { 
	      Talpha[v][j][d]   = sc;
	      Tkshadow[v][j][d] = k;
	      /*}*/
	    }
	  }
	  /* two additional special cases in trCYK (these are not in standard CYK) */
	  /* special case 1: k == 0 (full sequence aligns to BEGL_S left child */
	  if((sc = Jalpha[y][j][d]) > Lalpha[v][j][d]) { 
	    Lalpha[v][j][d]   = sc;
	    Lkshadow[v][j][d] = 0; /* k == 0 for this case, full sequence is on left */
	    Lkmode[v][j][d]   = TRMODE_J;
	  }
	  if((sc = Lalpha[y][j][d]) > Lalpha[v][j][d]) { 
	    Lalpha[v][j][d]   = sc;
	    Lkshadow[v][j][d] = 0; /* k == 0 for this case, full sequence is on left */
	    Lkmode[v][j][d]   = TRMODE_L;
	  }
	  /* special case 2: k == d (full sequence aligns to BEGR_S right child */
	  if((sc = Jalpha[z][j][d]) > Ralpha[v][j][d]) { 
	      Ralpha[v][j][d]   = sc;
	      Rkshadow[v][j][d] = d; /* k == d in this case, full sequence is on right */
	      Rkmode[v][j][d]   = TRMODE_J;
	  }
	  if((sc = Ralpha[z][j][d]) > Ralpha[v][j][d]) { 
	      Ralpha[v][j][d]   = sc;
	      Rkshadow[v][j][d] = d; /* k == d in this case, full sequence is on right */
	      Rkmode[v][j][d]   = TRMODE_R;
	  }
	}
      }
    } /* end of B_st recursion */

    /* Check for truncated alignment getting us to the root.
     * This is "off-shadow": if/when we trace back, we'll handle this
     * case separately (and we'll know to do it because we'll immediately
     * see a USED_TRUNC_BEGIN in the shadow matrix, 
     * telling us to jump right to state b; see below)
     */
    /* check for hit in J matrix (much like a normal local begin) */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st) { 
      if (Jalpha[v][j0][W] + trunc_penalty > bsc) { 
	b     = v;
	bsc   = Jalpha[v][j0][W] + trunc_penalty;
	bmode = TRMODE_J;
      }
    }
    /* check for truncated hit in L matrix */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st) { 
      if (Lalpha[v][j0][W] + trunc_penalty > bsc) { 
	b     = v;
	bsc   = Lalpha[v][j0][W] + trunc_penalty;
	bmode = TRMODE_L;
      }
      }
    /* check for truncated hit in R matrix */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st) { 
      if (Ralpha[v][j0][W] + trunc_penalty > bsc) { 
	b     = v;
	bsc   = Ralpha[v][j0][W] + trunc_penalty;
	bmode = TRMODE_R;
      }
    }	    
    /* check for truncated hit in T matrix */
    if(cm->sttype[v] == B_st) { 
      if (Talpha[v][j0][W] + trunc_penalty > bsc) { 
	b     = v;
	bsc   = Talpha[v][j0][W] + trunc_penalty;
	bmode = TRMODE_T;
      }	    
    }
    /* Check for whether we need to store best truncated score as
     * the optimal overall score, and if we need to put a flag in
     * the shadow matrix telling tr_alignT() to use the b we
     * return.
     */
    if (v == 0) { 
      if (bsc > Jalpha[0][j0][W]) {
	Jalpha[0][j0][W]   = bsc;
	Jyshadow[v][j0][W] = USED_TRUNC_BEGIN;
      }
    }
  } /* end loop over all v */
  /*FILE *fp1; fp1 = fopen("tmp.amx", "w");   cm_tr_mx_Dump(fp1, mx); fclose(fp1);*/
  /*FILE *fp2; fp2 = fopen("tmp.ashmx", "w"); cm_tr_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2);*/
  
  sc = Jalpha[0][j0][W];

  *ret_b     = b;    
  *ret_bsc   = bsc;  
  *ret_bmode = bmode;  
  *ret_sc    = sc;

  free(el_scA);

  ESL_DPRINTF1(("tr_cyk_align return sc: %f\n", sc));
  printf("tr_cyk_align return sc: %f\n", sc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: tr_alignT_hb()
 * Date:     EPN, Thu Sep  8 07:59:10 2011
 *           EPN 03.29.06 (fast_alignT_hb()
 *
 * Note:     based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *
 * Purpose:  Call either tr_cyk_align_hb() (if !<do_optacc>), 
 *           or tr_optimal_accuracy_align_hb()  (if  <do_optacc>),
 *           fill banded vjd shadow matrix in <shmx>; then trace back. 
 *           Append the trace to a given traceback, which already has state r at tr->n-1.
 *        
 *           If (<do_optacc>) then post_mx must != NULL.
 *
 * Returns:  <ret_sc>: if(!do_optacc): score of appended parsetree
 *                     if( do_optacc): avg posterior probability of all i0..j0 residues
 *                                     in optimally accurate alignment in tr 
 *           
 * Throws:  <eslOK>     on success
 *          <eslERANGE> if required CM_HB_MX exceeds <size_limit>, in 
 *                      this case, alignment has been aborted, tr has not been appended to
 *          <eslEINCOMPAT> if CM_HB_SHADOW_MX wasn't built for this CM 
 *
 */
int
tr_alignT_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, Parsetree_t *tr, 
		  int i0, int j0, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, 
		  int do_optacc, CM_TR_HB_MX *post_mx, float size_limit, float *ret_sc)
{
  int status;
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda_i;             /* stack that tracks bifurc parent of a right start */
  ESL_STACK *pda_c;             /* stack that tracks mode of bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, seq positions */
  int       k;			/* right subtree len for bifurcs */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  int       b;                  /* local begin state */
  float     bsc;                /* local begin score */
  int       jp_v;               /* j-jmin[v] for current j, and current v */
  int       dp_v;               /* d-hdmin[v][jp_v] for current j, current v, current d*/
  char      mode;               /* current truncation mode: TRMODE_J | TRMODE_L | TRMODE_R | TRMODE_T */
  char      prvmode, nxtmode;   /* previous, next truncation mode */
  char      bmode;              /* mode for best truncated CYK alignment */
  int       allow_S_trunc_end;  /* set to true to allow d==0 BEGL_S and BEGR_S truncated ends */

  /* contract check */
  if(dsq == NULL)                  ESL_FAIL(eslEINCOMPAT, errbuf, "tr_alignT_hb(), dsq == NULL.\n");
  if(cm->cp9b == NULL)             ESL_FAIL(eslEINCOMPAT, errbuf, "tr_alignT_hb(), cm->cp9b == NULL.\n");
  if(do_optacc && post_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "tr_alignT_hb(), do_optacc == TRUE but post_mx == NULL.\n");
			 
  /* pointers to cp9b data for convenience */
  CP9Bands_t  *cp9b = cm->cp9b;
  int         *jmin = cp9b->jmin;
  int         *jmax = cp9b->jmax;
  int       **hdmin = cp9b->hdmin;
  int       **hdmax = cp9b->hdmax;

  if(do_optacc) {
    cm_Fail("tr_optimal_accuracy_align_hb() not yet implemented");
#if 0
    status = tr_optimal_accuracy_align_hb(cm, errbuf, dsq, L, i0, j0, size_limit,
					  shmx,	     /* the banded shadow matrix, to expand and fill-in */
					  &b, &bsc,     /* if allow_begin is TRUE, gives info on optimal b */
					  mx,           /* the HMM banded mx to fill-in */
					  post_mx,      /* pre-calc'ed posterior matrix */
					  &sc);         /* avg post prob of all emissions in optimally accurate parsetree */
#endif
  }
  else {
    status = tr_cyk_align_hb(cm, errbuf, dsq, i0, j0, 
				  size_limit,       /* max size of DP matrix */
				  mx,               /* the HMM banded mx */
				  shmx,	            /* the HMM banded shadow matrix */
				  &b, &bsc, &bmode, /* optimal truncated begin state, score, and truncated mode */
				  &sc);             /* score of CYK parsetree */
  }
  if(status != eslOK) return status;

  pda_i = esl_stack_ICreate();
  pda_c = esl_stack_CCreate();
  if(pda_i == NULL) goto ERROR;
  if(pda_c == NULL) goto ERROR;
  j = j0;
  i = i0;
  d = j0-i0+1;
  v = 0;
  jp_v = j - jmin[v];
  assert(j >= jmin[v] && j <= jmax[v]);
  dp_v = d - hdmin[v][jp_v];
  assert(d >= cp9b->hdmin[v][jp_v] && d <= hdmax[v][jp_v]);

  /* Determine the root state v and marginal mode of the entire parse */
  v    = (shmx->Jyshadow[0][jp_v][dp_v] == USED_TRUNC_BEGIN) ? b     : 0;
  mode = (shmx->Jyshadow[0][jp_v][dp_v] == USED_TRUNC_BEGIN) ? bmode : shmx->Jyshadow[0][jp_v][dp_v];

  /* initialize the parsetree with the root */
  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, v, mode);

  while (1) {
    if(cm->sttype[v] != EL_st) printf("v: %4d  mode: %4d  j: %4d (%4d..%4d)  d: %4d", v, mode, j, jmin[v], jmax[v], d);
    else                       printf("v: %4d  mode: %4d  j: %4d             d: %4d EL\n", v, mode, j, d);

    if(cm->sttype[v] == S_st && d == 0 && (j < jmin[v] || j > jmax[v])) {
      /* special case: v is a BEGL_S or BEGR_S and j is outside v's j
       * bands, but we allow it if d == 0 b/c we're doing a truncated
       * end out of this state immediately.  This occurs if the parent
       * bif state emitted its full sequence via the other child
       * (BEGR_S or BEGL_S)
       */
      assert((cm->stid[v] == BEGL_S && mode == TRMODE_R) || (cm->stid[v] == BEGR_S && mode == TRMODE_L));
      ESL_DASSERT1(((cm->stid[v] == BEGL_S && mode == TRMODE_R) || (cm->stid[v] == BEGR_S && mode == TRMODE_L)));
      allow_S_trunc_end = TRUE; /* this sets yoffset to USED_TRUNC_END in the final 'else' of below code block */
    }
    else if (cm->sttype[v] != EL_st){ /* normal case, determine jp_v, dp_v, j, d offset values given bands */
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
      allow_S_trunc_end = FALSE;
      assert(j >= jmin[v]        && j <= jmax[v]);
      printf("(%4d..%4d)\n", hdmin[v][jp_v], hdmax[v][jp_v]);
      assert(d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]);
      ESL_DASSERT1((j >= jmin[v]        && j <= jmax[v]));
      ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
    }
    if (cm->sttype[v] == B_st) {
      /* get k, the len of right fragment */
      if     (mode == TRMODE_J) k = shmx->Jkshadow[v][jp_v][dp_v];
      else if(mode == TRMODE_L) k = shmx->Lkshadow[v][jp_v][dp_v];
      else if(mode == TRMODE_R) k = shmx->Rkshadow[v][jp_v][dp_v];
      else if(mode == TRMODE_T) k = shmx->Tkshadow[v][jp_v][dp_v];
      else                           ESL_FAIL(eslEINVAL, errbuf, "bogus truncation mode for B state: %d\n", mode);
      /* if k is 0, right fragment is of length 0 */
      /* determine mode of right child */
      prvmode = mode;
      if     (mode == TRMODE_J) ; /* do nothing, in J mode, right child mode remains TRMODE_J */
      else if(mode == TRMODE_L) mode = TRMODE_L; /* in TRMODE_L, right child is always Left marginal */
      else if(mode == TRMODE_R) mode = shmx->Rkmode[v][jp_v][dp_v];     
      else if(mode == TRMODE_T) mode = TRMODE_L; /* in TRMODE_T, right child is always Left marginal */

      /* Store info about the right fragment that we'll retrieve later:
       */
      if((status = esl_stack_CPush(pda_c, mode))    != eslOK) goto ERROR;  /* remember the mode of right child */
      if((status = esl_stack_IPush(pda_i, j))       != eslOK) goto ERROR;  /* remember the end j    */
      if((status = esl_stack_IPush(pda_i, k))       != eslOK) goto ERROR;  /* remember the subseq length k for right child */
      if((status = esl_stack_IPush(pda_i, tr->n-1)) != eslOK) goto ERROR;  /* remember the trace index of the parent B state */

      /* Determine mode of left start state */
      if     (prvmode == TRMODE_J) mode = TRMODE_J;
      else if(prvmode == TRMODE_L) mode = shmx->Lkmode[v][jp_v][dp_v]; 
      else if(prvmode == TRMODE_R) mode = TRMODE_R; /* for R mode, left child is always Right marginal */
      else if(prvmode == TRMODE_T) mode = TRMODE_R; /* for T mode, left child is always Right marginal */

      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
   
      y = cm->cfirst[v];
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y, mode);
      v = y;

#if DEBUG1
      printf("KACHOW added BEGL_S, dumping parsetree (prvmode: %d mode: %d:\n", prvmode, mode);
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda_i, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda_i, &d);
      esl_stack_IPop(pda_i, &j);
      esl_stack_CPop(pda_c, &mode);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNodewithMode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y, mode);
#if DEBUG1
      printf("KACHOW added E or EL, dumping parsetree:\n");
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif

      v = y;
    } 
    else {
      /* get yoffset */
      if(allow_S_trunc_end) { 
	yoffset = USED_TRUNC_END; /* nxt mode is irrelevant in this case */
      }
      else { 
	if     (mode == TRMODE_J) yoffset = shmx->Jyshadow[v][jp_v][dp_v];
	else if(mode == TRMODE_L) yoffset = shmx->Lyshadow[v][jp_v][dp_v];
	else if(mode == TRMODE_R) yoffset = shmx->Ryshadow[v][jp_v][dp_v];
	else if(mode == TRMODE_T) ESL_FAIL(eslEINVAL, errbuf, "truncation mode T for non B states");
	else                      ESL_FAIL(eslEINVAL, errbuf, "bogus truncation mode %d\n", mode);
      }
#if DEBUG1
      printf("v: %d std mode: %d yoffset: %d ", v, mode, yoffset);
#endif
      /* determine nxtmode, and correct yoffset */
      if     (yoffset == USED_TRUNC_END)   { yoffset = USED_TRUNC_END; } /* nxtmode is irrelevant in this case */
      else if(yoffset == USED_EL)          { yoffset = USED_EL;        } /* nxtmode is irrelevant in this case */
      else if(yoffset >= TRMODE_R_OFFSET)  { nxtmode = TRMODE_R; yoffset -= TRMODE_R_OFFSET; }
      else if(yoffset >= TRMODE_L_OFFSET)  { nxtmode = TRMODE_L; yoffset -= TRMODE_L_OFFSET; }
      else if(yoffset >= TRMODE_J_OFFSET)  { nxtmode = TRMODE_J; yoffset -= TRMODE_J_OFFSET; }
      else                                  ESL_FAIL(eslEINVAL, errbuf, "yoffset out of bounds: %d\n", yoffset);
#if DEBUG1
      printf("new yoffset: %d nxtmode: %d\n", yoffset, nxtmode);
      if(mode == TRMODE_J) printf("HEYA J v: %4d j: %4d d: %4d mode: %4d yoffset: %4d nxtmode: %4d\n", v, j, d, mode, yoffset, nxtmode);
      if(mode == TRMODE_L) printf("HEYA L v: %4d j: %4d d: %4d mode: %4d yoffset: %4d nxtmode: %4d\n", v, j, d, mode, yoffset, nxtmode);
      if(mode == TRMODE_R) printf("HEYA R v: %4d j: %4d d: %4d mode: %4d yoffset: %4d nxtmode: %4d\n", v, j, d, mode, yoffset, nxtmode);
#endif
      switch (cm->sttype[v]) { 
      case  D_st:
	break;
      case MP_st:
	if ( mode == TRMODE_J )          i++;
	if ( mode == TRMODE_L && d > 0 ) i++;
	if ( mode == TRMODE_J )          j--;
	if ( mode == TRMODE_R && d > 0 ) j--;
	break;
      case ML_st:
	if ( mode == TRMODE_J )          i++;
	if ( mode == TRMODE_L && d > 0 ) i++;
	break;
      case MR_st:
	if ( mode == TRMODE_J )          j--;
	if ( mode == TRMODE_R && d > 0 ) j--;
	break;
      case IL_st:
	if ( mode == TRMODE_J )          i++;
	if ( mode == TRMODE_L && d > 0 ) i++;
	break;
      case IR_st:
	if ( mode == TRMODE_J )          j--;
	if ( mode == TRMODE_R && d > 0 ) j--;
	break;
      case  S_st:
	break;
      default: ESL_FAIL(eslEINVAL, errbuf, "bogus state type %d \n", cm->sttype[v]);
      }
      d = j-i+1;

      if (yoffset == USED_EL || yoffset == USED_TRUNC_END) 
	{	/* a local alignment end  or a truncation end */
	  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M, mode);
	  v = cm->M;		/* now we're in EL (if USED_TRUNC_END, we act like we are) */
#if DEBUG1
	  printf("KACHOW added USED_EL or USED_TRUNC_END, dumping parsetree:\n");
	  ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
	}
      else if (yoffset == USED_TRUNC_BEGIN) 
	{ /* truncated begin; can only happen once, from root */
	  mode = nxtmode; /* this was set as bmode above */
	  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b, mode);
	  v    = b;
#if DEBUG1
	  printf("KACHOW added USED_TRUNC_BEGIN, dumping parsetree:\n");
	  ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; shouldn't happen in truncated alignment */
	  ESL_FAIL(status, errbuf, "local begin seen in truncated alignment, shouldn't happen.");
	}
      else 
	{
	  mode = nxtmode;
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y, mode);
#if DEBUG1
	  printf("STD yoffset: %d\n", yoffset);
	  printf("KACHOW added standard, dumping parsetree:\n");
	  ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
	  v    = y;
	}
    }
  }
  esl_stack_Destroy(pda_i);  /* it should be empty; we could check; naaah. */
  esl_stack_Destroy(pda_c);  /* it should be empty; we could check; naaah. */
  if(ret_sc != NULL) *ret_sc = sc;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVERREACHED */
}

/* Function: tr_alignT()
 * Date:     EPN, Sat Sep 10 11:25:37 2011
 *           EPN, Sun Nov 18 19:21:30 2007 [fast_alignT()]
 *
 * Note:     based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *
 * Purpose:  Call either tr_cyk_align() (if !<do_optacc>), 
 *           or tr_optimal_accuracy_align()  (if  <do_optacc>),
 *           fill vjd shadow matrix in <shadow>; then trace back. 
 *           and append to an existing but empty parsetree tr.
 *        
 *           If (<do_optacc>) then post_mx must != NULL.
 *
 * Returns:  <ret_sc>: if(!do_optacc): score of appended parsetree
 *                     if( do_optacc): avg posterior probability of all i0..j0 residues
 *                                     in optimally accurate alignment in tr 
 *           
 * Throws:  <eslOK>     on success
 *          <eslERANGE> if required DP matrix exceeds <size_limit>, in 
 *                      this case, alignment has been aborted, tr has not been appended to
 *
 */
int
tr_alignT(CM_t *cm, char *errbuf, ESL_DSQ *dsq, Parsetree_t *tr, 
	  int i0, int j0, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
	  int do_optacc, CM_TR_MX *post_mx, float size_limit, float *ret_sc)
{
  int       status;
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda_i;             /* stack that tracks bifurc parent of a right start */
  ESL_STACK *pda_c;             /* stack that tracks mode of bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, seq positions */
  int       k;			/* right subtree len for bifurcs */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  int       b;                  /* local begin state */
  float     bsc;                /* local begin score */
  /* variables specific to truncated version */
  char      mode;               /* current truncation mode: TRMODE_J | TRMODE_L | TRMODE_R | TRMODE_T */
  char      prvmode, nxtmode;   /* previous, next truncation mode */
  char      bmode;              /* mode for best truncated CYK alignment */
#if 0 
  int       allow_S_trunc_end;  /* set to true to allow d==0 BEGL_S and BEGR_S truncated ends */
#endif

  /* contract check */
  if(dsq == NULL)                  ESL_FAIL(eslEINCOMPAT, errbuf, "tr_alignT_hb(), dsq == NULL.\n");
  if(do_optacc && post_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "tr_alignT_hb(), do_optacc == TRUE but post_mx == NULL.\n");
			 
  if(do_optacc) {
    cm_Fail("tr_optimal_accuracy_align() not yet implemented");
#if 0
    status = tr_optimal_accuracy_align(cm, errbuf, dsq, L, i0, j0, size_limit,
				       shmx,	     /* the banded shadow matrix, to expand and fill-in */
				       &b, &bsc,     /* if allow_begin is TRUE, gives info on optimal b */
				       mx,           /* the HMM banded mx to fill-in */
				       post_mx,      /* pre-calc'ed posterior matrix */
				       &sc);         /* avg post prob of all emissions in optimally accurate parsetree */
#endif
  }
  else {
    status = tr_cyk_align(cm, errbuf, dsq, i0, j0, 
			  size_limit,       /* max size of DP matrix */
			  mx,               /* the HMM banded mx */
			  shmx,	            /* the HMM banded shadow matrix */
			  &b, &bsc, &bmode, /* optimal truncated begin state, score, and truncated mode */
			  &sc);             /* score of CYK parsetree */
  }
  if(status != eslOK) return status;

  pda_i = esl_stack_ICreate();
  pda_c = esl_stack_CCreate();
  if(pda_i == NULL) goto ERROR;
  if(pda_c == NULL) goto ERROR;
  j = j0;
  i = i0;
  d = j0-i0+1;
  v = 0;

  /* Determine the root state v and marginal mode of the entire parse */
  v    = (shmx->Jyshadow[0][j][d] == USED_TRUNC_BEGIN) ? b     : 0;
  mode = (shmx->Jyshadow[0][j][d] == USED_TRUNC_BEGIN) ? bmode : shmx->Jyshadow[0][j][d];

  /* initialize the parsetree with the root */
  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, v, mode);

  while (1) {
    if(cm->sttype[v] != EL_st) printf("v: %4d  mode: %4d  j: %4d d: %4d\n", v, mode, j, d);
    else                       printf("v: %4d  mode: %4d  j: %4d d: %4d EL\n", v, mode, j, d);

#if 0
    if(cm->sttype[v] == S_st && d == 0) { 
      /* special case: v is a BEGL_S or BEGR_S and j is outside v's j
       * bands, but we allow it if d == 0 b/c we're doing a truncated
       * end out of this state immediately.  This occurs if the parent
       * bif state emitted its full sequence via the other child
       * (BEGR_S or BEGL_S)
       */
      assert((cm->stid[v] == BEGL_S && mode == TRMODE_R) || (cm->stid[v] == BEGR_S && mode == TRMODE_L));
      ESL_DASSERT1(((cm->stid[v] == BEGL_S && mode == TRMODE_R) || (cm->stid[v] == BEGR_S && mode == TRMODE_L)));
      allow_S_trunc_end = TRUE; /* this sets yoffset to USED_TRUNC_END in the final 'else' of below code block */
    }
    else if (cm->sttype[v] != EL_st){ /* normal case */
      allow_S_trunc_end = FALSE;
    }
#endif
    if (cm->sttype[v] == B_st) {
      /* get k, the len of right fragment */
      if     (mode == TRMODE_J) k = shmx->Jkshadow[v][j][d];
      else if(mode == TRMODE_L) k = shmx->Lkshadow[v][j][d];
      else if(mode == TRMODE_R) k = shmx->Rkshadow[v][j][d];
      else if(mode == TRMODE_T) k = shmx->Tkshadow[v][j][d];
      else                           ESL_FAIL(eslEINVAL, errbuf, "bogus truncation mode for B state: %d\n", mode);
      /* if k is 0, right fragment is of length 0 */
      /* determine mode of right child */
      prvmode = mode;
      if     (mode == TRMODE_J) ; /* do nothing, in J mode, right child mode remains TRMODE_J */
      else if(mode == TRMODE_L) mode = TRMODE_L; /* in TRMODE_L, right child is always Left marginal */
      else if(mode == TRMODE_R) mode = shmx->Rkmode[v][j][d];     
      else if(mode == TRMODE_T) mode = TRMODE_L; /* in TRMODE_T, right child is always Left marginal */

      /* Store info about the right fragment that we'll retrieve later:
       */
      if((status = esl_stack_CPush(pda_c, mode))    != eslOK) goto ERROR;  /* remember the mode of right child */
      if((status = esl_stack_IPush(pda_i, j))       != eslOK) goto ERROR;  /* remember the end j    */
      if((status = esl_stack_IPush(pda_i, k))       != eslOK) goto ERROR;  /* remember the subseq length k for right child */
      if((status = esl_stack_IPush(pda_i, tr->n-1)) != eslOK) goto ERROR;  /* remember the trace index of the parent B state */

      /* Determine mode of left start state */
      if     (prvmode == TRMODE_J) mode = TRMODE_J;
      else if(prvmode == TRMODE_L) mode = shmx->Lkmode[v][j][d]; 
      else if(prvmode == TRMODE_R) mode = TRMODE_R; /* for R mode, left child is always Right marginal */
      else if(prvmode == TRMODE_T) mode = TRMODE_R; /* for T mode, left child is always Right marginal */

      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
   
      y = cm->cfirst[v];
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y, mode);
      v = y;

#if DEBUG1
      printf("KACHOW added BEGL_S, dumping parsetree (prvmode: %d mode: %d:\n", prvmode, mode);
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda_i, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda_i, &d);
      esl_stack_IPop(pda_i, &j);
      esl_stack_CPop(pda_c, &mode);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNodewithMode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y, mode);
#if DEBUG1
      printf("KACHOW added E or EL, dumping parsetree:\n");
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif

      v = y;
    } 
    else {
      /* get yoffset */
#if 0
      if(allow_S_trunc_end) { 
	yoffset = USED_TRUNC_END; /* nxt mode is irrelevant in this case */
      }
      else { 
#endif
	if     (mode == TRMODE_J) yoffset = shmx->Jyshadow[v][j][d];
	else if(mode == TRMODE_L) yoffset = shmx->Lyshadow[v][j][d];
	else if(mode == TRMODE_R) yoffset = shmx->Ryshadow[v][j][d];
	else if(mode == TRMODE_T) ESL_FAIL(eslEINVAL, errbuf, "truncation mode T for non B states");
	else                      ESL_FAIL(eslEINVAL, errbuf, "bogus truncation mode %d\n", mode);
#if 0 
      }
#endif
#if DEBUG1
      printf("v: %d std mode: %d yoffset: %d ", v, mode, yoffset);
#endif
      /* determine nxtmode, and correct yoffset */
      if     (yoffset == USED_TRUNC_END)   { yoffset = USED_TRUNC_END; } /* nxtmode is irrelevant in this case */
      else if(yoffset == USED_EL)          { yoffset = USED_EL;        } /* nxtmode is irrelevant in this case */
      else if(yoffset >= TRMODE_R_OFFSET)  { nxtmode = TRMODE_R; yoffset -= TRMODE_R_OFFSET; }
      else if(yoffset >= TRMODE_L_OFFSET)  { nxtmode = TRMODE_L; yoffset -= TRMODE_L_OFFSET; }
      else if(yoffset >= TRMODE_J_OFFSET)  { nxtmode = TRMODE_J; yoffset -= TRMODE_J_OFFSET; }
      else                                  ESL_FAIL(eslEINVAL, errbuf, "yoffset out of bounds: %d\n", yoffset);
#if DEBUG1
      printf("new yoffset: %d nxtmode: %d\n", yoffset, nxtmode);
      if(mode == TRMODE_J) printf("HEYA J v: %4d j: %4d d: %4d mode: %4d yoffset: %4d nxtmode: %4d\n", v, j, d, mode, yoffset, nxtmode);
      if(mode == TRMODE_L) printf("HEYA L v: %4d j: %4d d: %4d mode: %4d yoffset: %4d nxtmode: %4d\n", v, j, d, mode, yoffset, nxtmode);
      if(mode == TRMODE_R) printf("HEYA R v: %4d j: %4d d: %4d mode: %4d yoffset: %4d nxtmode: %4d\n", v, j, d, mode, yoffset, nxtmode);
#endif
      switch (cm->sttype[v]) { 
      case  D_st:
	break;
      case MP_st:
	if ( mode == TRMODE_J )          i++;
	if ( mode == TRMODE_L && d > 0 ) i++;
	if ( mode == TRMODE_J )          j--;
	if ( mode == TRMODE_R && d > 0 ) j--;
	break;
      case ML_st:
	if ( mode == TRMODE_J )          i++;
	if ( mode == TRMODE_L && d > 0 ) i++;
	break;
      case MR_st:
	if ( mode == TRMODE_J )          j--;
	if ( mode == TRMODE_R && d > 0 ) j--;
	break;
      case IL_st:
	if ( mode == TRMODE_J )          i++;
	if ( mode == TRMODE_L && d > 0 ) i++;
	break;
      case IR_st:
	if ( mode == TRMODE_J )          j--;
	if ( mode == TRMODE_R && d > 0 ) j--;
	break;
      case  S_st:
	break;
      default: ESL_FAIL(eslEINVAL, errbuf, "bogus state type %d \n", cm->sttype[v]);
      }
      d = j-i+1;

      if (yoffset == USED_EL || yoffset == USED_TRUNC_END) 
	{	/* a local alignment end  or a truncation end */
	  if(mode == TRMODE_J) { /* TRMODE_J should be only way USED_EL is possible */
	    assert(yoffset == USED_EL);
	    InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M, mode);
#if DEBUG1
	    printf("KACHOW added USED_EL or USED_TRUNC_END, dumping parsetree:\n");
	    ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
	  }
	  v = cm->M;		/* now we're in EL (if USED_TRUNC_END, we act like we are) */
	}
      else if (yoffset == USED_TRUNC_BEGIN) 
	{ /* truncated begin; can only happen once, from root */
	  mode = nxtmode; /* this was set as bmode above */
	  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b, mode);
	  v    = b;
#if DEBUG1
	  printf("KACHOW added USED_TRUNC_BEGIN, dumping parsetree:\n");
	  ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; shouldn't happen in truncated alignment */
	  ESL_FAIL(status, errbuf, "local begin seen in truncated alignment, shouldn't happen.");
	}
      else 
	{
	  mode = nxtmode;
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y, mode);
#if DEBUG1
	  printf("STD yoffset: %d\n", yoffset);
	  printf("KACHOW added standard, dumping parsetree:\n");
	  ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
#endif
	  v    = y;
	}
    }
  }
  esl_stack_Destroy(pda_i);  /* it should be empty; we could check; naaah. */
  esl_stack_Destroy(pda_c);  /* it should be empty; we could check; naaah. */
  if(ret_sc != NULL) *ret_sc = sc;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVERREACHED */
}

/* Function: TrAlignHB()
 * Incept:   EPN, Thu Sep  8 08:55:26 2011
 *           EPN, Fri Oct 26 09:31:43 2007 [FastAlignHB()]
 * 
 * Note:     based on CYKInside_b_jd() [11.04.05] which was based on CYKInside_b() 
 *           which was based on CYKInside() [SRE, Sun Jun  3 19:48:33 2001 [St. Louis]]
 *
 * Purpose: Wrapper for the tr_alignT_hb() routine - solve a full
 *           alignment problem using trCYK, truncated optimal accuracy
 *           or sampling and return the traceback and the score,
 *           without dividing & conquering, but by using bands on the
 *           j and d dimensions of the DP matrix.  Bands derived by
 *           HMM Forward/Backward runs. Optionally return a posterior
 *           code string.
 *           
 *           Identical to TrAlign() but HMM bands are used here.
 * 
 *           Input arguments allow this function to be run in 6 'modes':
 *
 *           mode      returns                 arguments
 *           ----  ----------------  ----------------------------------------
 *                 tr        pcodes  do_optacc  do_sample post_mx   ret_pcode
 *                 ----------------  ----------------------------------------
 *              1. CYK       no      FALSE      FALSE      NULL      NULL
 *              2. CYK       yes     FALSE      FALSE     !NULL     !NULL
 *              3. Opt acc   no      TRUE       FALSE     !NULL      NULL
 *              4. Opt acc   yes     TRUE       FALSE     !NULL     !NULL
 *              5. sampled   no      FALSE      TRUE       NULL      NULL
 *              6. sampled   yes     FALSE      TRUE      !NULL     !NULL
 *
 *           CYK parsetrees are most likely parsetree, 'Opt acc' parsetrees
 *           are Holmes/Durbin optimally accurate parsetrees, the parse that
 *           maximizes the posterior probability that goes through the parses cells
 *           of the DP matrix. A sampled parsetree is a parsetree sampled from
 *           an Inside matrix based on it's probability.
 *
 *           Note: if ret_tr is NULL, parsetree is not returned.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           r         - source of randomness, must be non-NULL only if do_sample==TRUE
 *           dsq       - the digitized sequence, 1..L
 *           L         - length of sequence 
 *           i0        - start of target subsequence (often 1, beginning of dsq)
 *           j0        - end of target subsequence (often L, end of dsq)
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx        - the main dp matrix, only cells within bands in cm->cp9b will be valid. 
 *           shmx      - the HMM banded shadow matrix to fill in and traceback, same cells as mx are valid.
 *           do_optacc - TRUE to not do CYK alignment, determine the Holmes/Durbin optimally 
 *                       accurate parsetree in ret_tr, requires post_mx != NULL
 *           do_sample - TRUE to sample a parsetree from the Inside matrix
 *           post_mx   - dp matrix for posterior calculation, can be NULL only if !do_optacc
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_pcode - RETURN: posterior code, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_sc    - if(!do_optacc): score of the alignment in bits.
 *                       if( do_optacc): average posterior probability of all L aligned residues 
 *                       in optimally accurate alignment
 *           ret_ins_sc- if(do_optacc || ret_pcode != NULL): inside score of sequence in bits
 *                       else: must be NULL (inside will not be run)
 * 
 * Returns: <ret_tr>, <ret_pcode>, <ret_sc>, see 'Args' section
 * 
 * Throws:  <eslOK> on success; 
 *          <eslERANGE> if required CM_HB_MX for FastInsideAlignHB(), FastOutsideAlignHB() or
 *                      cyk_align_hb() exceeds <size_limit>, in this 
 *                      case, alignment has been aborted, ret_* variables are not valid 
 */

int
TrAlignHB(CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, 
	  int do_optacc, int do_sample, CM_TR_HB_MX *post_mx, Parsetree_t **ret_tr, char **ret_pcode, float *ret_sc, float *ret_ins_sc)
{
  int          status;
  Parsetree_t *tr;
  float        sc;
  float        ins_sc; /* inside score */
  int          do_post;
  char        *pcode;
  int          have_pcodes;
  have_pcodes = (ret_pcode != NULL) ? TRUE : FALSE;
  do_post = (do_optacc || have_pcodes) ? TRUE : FALSE;

  /* Contract check */
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), mx is NULL.");
  if(post_mx == NULL && have_pcodes)   ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), post_mx == NULL but ret_pcode{1|2} != NULL.");
  if(do_optacc && post_mx == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), do_optacc is TRUE, but post_mx == NULL.");
  if((!do_post) && ret_ins_sc != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), do_post is FALSE, but ret_ins_sc != NULL.");
  if(shmx == NULL)                     ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), shadow matrix shmx == NULL.");
  if(do_optacc && do_sample)           ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), do_optacc and do_sample are both TRUE.");
  if(do_sample && r == NULL)           ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), do_sample but r is NULL.");
  if(!do_sample && r != NULL)          ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), do_sample is FALSE, but r is non-NULL.");
  if(do_sample && i0 != 1)             ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlignHB(), do_sample but i0!=1.");
  /* PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, (j0-i0+1)); */

  /* if do_post, fill Inside, Outside, Posterior matrices, in that order */
  /* if do_sample (and !do_post) fill Inside and sample from it */
  if(do_post || do_sample) { 
    cm_Fail("ERROR, do_post || do_sample not yet implemented");
#if 0 
    if((status = TrInsideAlignHB (cm, errbuf, dsq, i0, j0, size_limit, mx, &ins_sc)) != eslOK) return status;
    if(do_sample) { 
      if((status = SampleFromInsideHB(r, cm, errbuf, dsq, j0-i0+1, mx, &tr, &sc)) != eslOK) return status; 
    }
    if(do_post) { /* Inside was called above, now do Outside, then Posterior */
      if((status = TrOutsideAlignHB(cm, errbuf, dsq, i0, j0, size_limit, post_mx, mx, ((cm->align_opts & CM_ALIGN_CHECKINOUT) && (! cm->flags & CMH_LOCAL_END)), NULL)) != eslOK) return status;
      /* Note: we can only check the posteriors in TrOutsideAlignHB() if local begin/ends are off */
      if((status = CMPosteriorHB(cm, errbuf, i0, j0, size_limit, mx, post_mx, post_mx)) != eslOK) return status;   
      if(cm->align_opts & CM_ALIGN_CHECKINOUT) { 
	if((status = CMCheckPosteriorHB(cm, errbuf, i0, j0, post_mx)) != eslOK) return status;
	printf("\nHMM banded posteriors checked.\n\n");
      }
      if(ret_ins_sc != NULL) *ret_ins_sc = ins_sc; 
    }
#endif
  }

  if(!do_sample) { /* if do_sample, we already have a parsetree */
    /* Create the parse tree, we don't initialize here (as we do in
     * FastAlignHB()) b/c we need to know optimal begin state first, 
     * we initialize after we've performed CYK or optimal accuracy 
     */
    tr = CreateParsetree(100);
    /* Fill in the parsetree (either CYK or optimally accurate (if do_optacc)) */
    if((status = tr_alignT_hb(cm, errbuf, dsq, tr, i0, j0, mx, shmx, do_optacc, post_mx, size_limit, &sc)) != eslOK) return status;
  }

  if(have_pcodes) {
    cm_Fail("ERROR, have_pcodes not yet implemented");
#if 0
    if((status = CMPostCodeHB(cm, errbuf, i0, j0, post_mx, tr, TRUE, &pcode, (do_optacc ? &sc : NULL))) != eslOK) return status;
    *ret_pcode = pcode;
#endif
  }
  else if(do_optacc) { /* call CMPostCodeHB() to get the average residue posterior probability label, but not post codes */ 
    cm_Fail("ERROR, do_optacc not yet implemented");
#if 0
    if((status = CMPostCodeHB(cm, errbuf, i0, j0, post_mx, tr, TRUE, NULL, &sc)) != eslOK) return status;
#endif
  }

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  if (ret_sc != NULL) *ret_sc = sc;
  ESL_DPRINTF1(("returning from TrAlignHB() sc : %f\n", sc)); 
  return eslOK;
}


/* Function: TrAlign()
 * Incept:   EPN, Sat Sep 10 12:58:09 2011
 * 
 * Purpose: Wrapper for the tr_alignTb() routine - solve a full
 *           alignment problem using trCYK, truncated optimal accuracy
 *           or sampling and return the traceback and the score.
 *           
 *           Identical to TrAlignHB() but HMM bands are not used here.
 * 
 *           Input arguments allow this function to be run in 6 'modes':
 *
 *           mode      returns                 arguments
 *           ----  ----------------  ----------------------------------------
 *                 tr        pcodes  do_optacc  do_sample post_mx   ret_pcode
 *                 ----------------  ----------------------------------------
 *              1. CYK       no      FALSE      FALSE      NULL      NULL
 *              2. CYK       yes     FALSE      FALSE     !NULL     !NULL
 *              3. Opt acc   no      TRUE       FALSE     !NULL      NULL
 *              4. Opt acc   yes     TRUE       FALSE     !NULL     !NULL
 *              5. sampled   no      FALSE      TRUE       NULL      NULL
 *              6. sampled   yes     FALSE      TRUE      !NULL     !NULL
 *
 *           CYK parsetrees are most likely parsetree, 'Opt acc' parsetrees
 *           are Holmes/Durbin optimally accurate parsetrees, the parse that
 *           maximizes the posterior probability that goes through the parses cells
 *           of the DP matrix. A sampled parsetree is a parsetree sampled from
 *           an Inside matrix based on it's probability.
 *
 *           Note: if ret_tr is NULL, parsetree is not returned.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           r         - source of randomness, must be non-NULL only if do_sample==TRUE
 *           dsq       - the digitized sequence, 1..L
 *           L         - length of sequence 
 *           i0        - start of target subsequence (often 1, beginning of dsq)
 *           j0        - end of target subsequence (often L, end of dsq)
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx        - the main dp matrix, only cells within bands in cm->cp9b will be valid. 
 *           shmx      - the HMM banded shadow matrix to fill in and traceback, same cells as mx are valid.
 *           do_optacc - TRUE to not do CYK alignment, determine the Holmes/Durbin optimally 
 *                       accurate parsetree in ret_tr, requires post_mx != NULL
 *           do_sample - TRUE to sample a parsetree from the Inside matrix
 *           post_mx   - dp matrix for posterior calculation, can be NULL only if !do_optacc
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_pcode - RETURN: posterior code, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_sc    - if(!do_optacc): score of the alignment in bits.
 *                       if( do_optacc): average posterior probability of all L aligned residues 
 *                       in optimally accurate alignment
 *           ret_ins_sc- if(do_optacc || ret_pcode != NULL): inside score of sequence in bits
 *                       else: must be NULL (inside will not be run)
 * 
 * Returns: <ret_tr>, <ret_pcode>, <ret_sc>, see 'Args' section
 * 
 * Throws:  <eslOK> on success; 
 *          <eslERANGE> if required CM_TR_MX for TrInsideAlign(), TrOutsideAlign() or
 *                      tr_cyk_align() exceeds <size_limit>, in this 
 *                      case, alignment has been aborted, ret_* variables are not valid 
 */
int
TrAlign(CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
	int do_optacc, int do_sample, CM_TR_MX *post_mx, Parsetree_t **ret_tr, char **ret_pcode, float *ret_sc, float *ret_ins_sc)
{
  int          status;
  Parsetree_t *tr;
  float        sc;
  float        ins_sc; /* inside score */
  int          do_post;
  char        *pcode;
  int          have_pcodes;
  have_pcodes = (ret_pcode != NULL) ? TRUE : FALSE;
  do_post = (do_optacc || have_pcodes) ? TRUE : FALSE;

  /* Contract check */
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), mx is NULL.");
  if(post_mx == NULL && have_pcodes)   ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), post_mx == NULL but ret_pcode{1|2} != NULL.");
  if(do_optacc && post_mx == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), do_optacc is TRUE, but post_mx == NULL.");
  if((!do_post) && ret_ins_sc != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), do_post is FALSE, but ret_ins_sc != NULL.");
  if(shmx == NULL)                     ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), shadow matrix shmx == NULL.");
  if(do_optacc && do_sample)           ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), do_optacc and do_sample are both TRUE.");
  if(do_sample && r == NULL)           ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), do_sample but r is NULL.");
  if(!do_sample && r != NULL)          ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), do_sample is FALSE, but r is non-NULL.");
  if(do_sample && i0 != 1)             ESL_FAIL(eslEINCOMPAT, errbuf, "TrAlign(), do_sample but i0!=1.");

  /* if do_post, fill Inside, Outside, Posterior matrices, in that order */
  /* if do_sample (and !do_post) fill Inside and sample from it */
  if(do_post || do_sample) { 
    cm_Fail("ERROR, do_post || do_sample not yet implemented");
#if 0 
    if((status = TrInsideAlign(cm, errbuf, dsq, i0, j0, size_limit, mx, &ins_sc)) != eslOK) return status;
    if(do_sample) { 
      if((status = SampleFromTrInside(r, cm, errbuf, dsq, j0-i0+1, mx, &tr, &sc)) != eslOK) return status; 
    }
    if(do_post) { /* Inside was called above, now do Outside, then Posterior */
      if((status = TrOutsideAlign(cm, errbuf, dsq, i0, j0, size_limit, post_mx, mx, ((cm->align_opts & CM_ALIGN_CHECKINOUT) && (! cm->flags & CMH_LOCAL_END)), NULL)) != eslOK) return status;
      /* Note: we can only check the posteriors in FastOutsideAlignHB() if local begin/ends are off */
      if((status = CMTrPosterior(cm, errbuf, i0, j0, size_limit, mx, post_mx, post_mx)) != eslOK) return status;   
      if(cm->align_opts & CM_ALIGN_CHECKINOUT) { 
	if((status = CMTrCheckPosteriorHB(cm, errbuf, i0, j0, post_mx)) != eslOK) return status;
	printf("\nHMM banded posteriors checked.\n\n");
      }
      if(ret_ins_sc != NULL) *ret_ins_sc = ins_sc; 
    }
#endif
  }

  if(!do_sample) { /* if do_sample, we already have a parsetree */
    /* Create the parse tree, we don't initialize here (as we do in
     * FastAlign()) b/c we need to know optimal begin state first, 
     * we initialize after we've performed CYK or optimal accuracy 
     */
    tr = CreateParsetree(100);
    /* Fill in the parsetree (either CYK or optimally accurate (if do_optacc)) */
    if((status = tr_alignT(cm, errbuf, dsq, tr, i0, j0, mx, shmx, do_optacc, post_mx, size_limit, &sc)) != eslOK) return status;
  }

  if(have_pcodes) {
    cm_Fail("ERROR, have_pcodes not yet implemented");
#if 0
    if((status = CMPostCodeHB(cm, errbuf, i0, j0, post_mx, tr, TRUE, &pcode, (do_optacc ? &sc : NULL))) != eslOK) return status;
    *ret_pcode = pcode;
#endif
  }
  else if(do_optacc) { /* call CMPostCodeHB() to get the average residue posterior probability label, but not post codes */ 
    cm_Fail("ERROR, do_optacc not yet implemented");
#if 0
    if((status = CMPostCodeHB(cm, errbuf, i0, j0, post_mx, tr, TRUE, NULL, &sc)) != eslOK) return status;
#endif
  }

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  if (ret_sc != NULL) *ret_sc = sc;
  ESL_DPRINTF1(("returning from TrAlignHB() sc : %f\n", sc)); 
  return eslOK;
}


/* Function: TrInsideAlignHB()
 * Date:     EPN, Mon Sep 12 04:32:00 2011   
 *
 * Purpose: Run the truncated inside algorithm on a target sequence
 *           using bands in the j and d dimensions of the DP
 *           matrix. Bands were obtained from an HMM Forward-Backward
 *           parse of the target sequence. Uses float log odds scores.
 * 
 *           Very similar to tr_cyk_align_hb(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - we do Inside, not CYK
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *
 *           This function complements TrOutsideAlignHB().
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score
 * 
 * Returns:  <ret_sc>
 *
 * Throws:  <eslOK> on success
 *          <eslERANGE> if required CM_HB_MX for exceeds <size_limit>, 
 *                      in this case, alignment has been aborted, ret_sc is not valid
 */
int
TrInsideAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_HB_MX *mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,jp,d,i,k;	/* indices in sequence dimensions */
  float    sc, Lsc, Rsc;/* temporary scores */
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
  int      dp_v, dp_y, dp_z;   /* d index for state v/y/z in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp;                 /* ESL_MAX(d-sd, 0) */
  int      dp_y_sd;            /* dp_y - sd */
  int      dp_y_sdr;           /* dp_y - sdr */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Wp;                 /* W oalso changes depending on state */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  /* variables related to truncated alignment (not in FastInsideAlignHB()) */
  float    trunc_penalty = 0.; /* penalty in bits for a truncated hit */
  int      do_J_v, do_J_y, do_J_z; /* is J matrix valid for state v, y, z? */
  int      do_L_v, do_L_y, do_L_z; /* is L matrix valid for state v, y, z? */
  int      do_R_v, do_R_y, do_R_z; /* is R matrix valid for state v, y, z? */
  int      do_T_v, do_T_y, do_T_z; /* is T matrix valid for state v, y, z? */

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
  float ***Jalpha  = mx->Jdp; /* pointer to the Jalpha DP matrix */
  float ***Lalpha  = mx->Ldp; /* pointer to the Lalpha DP matrix */
  float ***Ralpha  = mx->Rdp; /* pointer to the Ralpha DP matrix */
  float ***Talpha  = mx->Tdp; /* pointer to the Talpha DP matrix */

  /* Allocations and initializations */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */

  /* grow the matrix based on the current sequence and bands */
  if((status = cm_tr_hb_mx_GrowTo(cm, mx, errbuf, cp9b, W, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  ESL_STOPWATCH *w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(mx->Jncells_valid > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(mx->Lncells_valid > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(mx->Rncells_valid > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(mx->Tncells_valid > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " Matrix init CPU time: ");

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) { /* all the way down to root, different from other scanners */
    float const *esc_v   = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v   = cm->tsc[v];  /* transition scores for state v */
    float const *lmesc_v = cm->lmesc[v];
    float const *rmesc_v = cm->rmesc[v];
    sd     = StateDelta(cm->sttype[v]);
    sdr    = StateRightDelta(cm->sttype[v]);
    jn     = jmin[v];
    jx     = jmax[v];
    do_J_v = cp9b->do_J[v];
    do_L_v = cp9b->do_L[v];
    do_R_v = cp9b->do_R[v];
    do_T_v = cp9b->do_T[v];

    /* re-initialize the J deck if we can do a local end from v */
    if(do_J_v) { 
      if(NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	    dp = ESL_MAX(d-sd, 0);
	    Jalpha[v][jp_v][dp_v] = el_scA[dp] + cm->endsc[v];
	    /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	  }
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	if(do_J_v) Jalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	if(do_L_v) Lalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	if(do_R_v) Ralpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IL states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1] 
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr valid for state y? */
	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	  i    = j - d + 1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */

	  /* We need to treat R differently from and J and L here, by
	   * doing separate 'for (yoffset...' loops for J and R
	   * because we have to fully calculate Jalpha[v][jp_v][dp_v])
	   * before we can start to calculate Ralpha[v][jp_v][dp_v].
	   */
	  /* Handle J and L first */
	  if(do_J_v || do_L_v) { 
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_J_y = cp9b->do_J[y];
	      do_L_y = cp9b->do_L[y];
	      if(do_J_y || do_L_y) { 
		jp_y_sdr = j - jmin[y] - sdr;
		
		if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		  if(do_J_v && do_J_y) Jalpha[v][jp_v][dp_v] = FLogsum(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
		  if(do_L_v && do_L_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
		}
	      }
	    }
	    if(do_J_v) { 
	      Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]];
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    if(do_L_v) { 
	      Lalpha[v][jp_v][dp_v] = (d >= 2) ? Lalpha[v][jp_v][dp_v] + esc_v[dsq[i]] : esc_v[dsq[i]];
	      Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    i--;
	  }

	  if(do_R_v) { 
	    /* Handle R separately */
	    Rsc = Ralpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_R_y = cp9b->do_R[y];
	      do_J_y = cp9b->do_J[y];
	      if(do_J_y || do_R_y) { 
		jp_y_sdr = j - jmin[y] - sdr;
		
		/* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,L above) */
		if((d) >= hdmin[y][jp_y_sdr] && (d) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y = d - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		  if(do_J_y) Rsc = FLogsum(Rsc, Jalpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]);
		  if(do_R_y) Rsc = FLogsum(Rsc, Ralpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]);
		}
	      }
	    } /* end of for (yvalid_idx = 0... loop */
	    Ralpha[v][jp_v][dp_v] = Rsc; 
	    /* we use Rsc instead of Ralpha cell in above loop because
	     * Ralpha[v][jp_v][dp_v] may be the same cell as
	     * Ralpha[y][jp_y_sdr][dp_y] if we're an IL state
	     */
	  }
	}
      }
    }
    else if(cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IR states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1].
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */

      /* The first MR_st/IR_st 'for (j...' loop is for J and R matrices which use the same set of j values */
      if(do_J_v || do_R_v) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  j_sdr = j - sdr;
	  
	  /* determine which children y we can legally transit to for v, j */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	  
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	    /* We need to treat L differently from and J and R here, by
	     * doing separate 'for (yoffset...' loops for J because we
	     * have to fully calculate Jalpha[v][jp_v][dp_v]) before we
	     * can start to calculate Lalpha[v][jp_v][dp_v].
	     */
	    /* Handle J and R first */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_J_y = cp9b->do_J[y];
	      do_R_y = cp9b->do_R[y];
	      if(do_J_y || do_R_y) { 
		jp_y_sdr = j - jmin[y] - sdr;
		
		if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		  dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		  if(do_J_v && do_J_y) Jalpha[v][jp_v][dp_v] = FLogsum(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
		  if(do_R_v && do_R_y) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
		}
	      }
	    }
	    if(do_J_v) { 
	      Jalpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    if(do_R_v) { 
	      Ralpha[v][jp_v][dp_v] = (d >= 2) ? Ralpha[v][jp_v][dp_v] + esc_v[dsq[j]] : esc_v[dsq[j]];
	      Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	  }
	}
      }

      if(do_L_v) { 
	/* The second MR_st/IR_st 'for (j...' loop is for the L matrix which use a different set of j values */
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  
	  /* determine which children y we can legally transit to for v, j */
	  /* we use 'j' and not 'j_sdr' here for the L matrix, differently from J and R matrices above */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j) >= jmin[y] && ((j) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j is valid for state y? */
	  
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	    Lsc = Lalpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      do_L_y = cp9b->do_L[y];
	      do_J_y = cp9b->do_J[y];
	      if(do_L_y || do_J_y) { 
		/* we use 'jp_y=j-min[y]' here, not 'jp_y_sdr=j-jmin[y]-sdr' (which we used in the corresponding loop for J,R above) */
		jp_y = j - jmin[y];
	      
		/* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,R above) */
		if((d) >= hdmin[y][jp_y] && (d) <= hdmax[y][jp_y]) { /* make sure d is valid for this v, j and y */
		  dp_y = d - hdmin[y][jp_y];
		  ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v] - hdmin[v][jp_v])));
		  ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y] - hdmin[y][jp_y])));
		  if(do_J_y) Lsc = FLogsum(Lsc, Jalpha[y][jp_y][dp_y] + tsc_v[yoffset]);
		  if(do_L_y) Lsc = FLogsum(Lsc, Lalpha[y][jp_y][dp_y] + tsc_v[yoffset]);
		}
	      }
	    } /* end of for (yvalid_idx = 0... loop */
	    Lalpha[v][jp_v][dp_v] = Lsc; 
	    /* we use Lsc instead of Lalpha cell in above loop because
	     * Lalpha[v][jp_v][dp_v] may be the same cell as
	     * Lalpha[y][jp_y_sdr][dp_y] if we're an IR state
	     */
	  }
	}
      }
    }
    else if(cm->sttype[v] == MP_st) { 
      /* MP states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	do_J_y = cp9b->do_J[y];
	do_L_y = cp9b->do_L[y];
	do_R_y = cp9b->do_R[y];
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];

	/* The first MP_st 'for (jp_v...' loop is for J and R matrices which use the same set of j values */
	/* j must satisfy:
	 * j >= jmin[v]
	 * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	 * j <= jmax[v]
	 * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	 * this reduces to two ESL_MAX calls
	 */
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	/* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	
	if((do_J_v && do_J_y) || (do_R_v && (do_J_y || do_R_y))) { 
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++, jp_y++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	    
	    if(do_J_v && do_J_y) { 
	      /* J matrix: */
	      /* d must satisfy:
	       * d >= hdmin[v][jp_v]
	       * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	       * d <= hdmax[v][jp_v]
	       * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	       * this reduces to two ESL_MAX calls
	       */
	      dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	      dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	      dpn       = dn - hdmin[v][jp_v];
	      dpx       = dx - hdmin[v][jp_v];
	      dp_y_sd   = dn - hdmin[y][jp_y_sdr] - sd;
	      
	      for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
		ESL_DASSERT1((dp_v      >= 0 && dp_v       <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		ESL_DASSERT1((dp_y_sd   >= 0 && dp_y_sd    <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		Jalpha[v][jp_v][dp_v] = FLogsum(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	      }
	    }
	    
	    if(do_R_v && (do_R_y || do_J_y)) { 
	      /* R matrix: */
	      /* d must satisfy:
	       * d >= hdmin[v][jp_v]
	       * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	       * d <= hdmax[v][jp_v]
	       * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	       * this reduces to two ESL_MAX calls
	       */
	      dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sdr);
	      dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sdr);
	      dpn       = dn - hdmin[v][jp_v];
	      dpx       = dx - hdmin[v][jp_v];
	      dp_y_sdr  = dn - hdmin[y][jp_y_sdr] - sdr;
	      /* for {L,R}alpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	      
	      for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
		/* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
		ESL_DASSERT1((dp_y_sdr  >= 0 && dp_y_sdr   <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		if(do_J_y) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sdr] + tsc); 
		if(do_R_y) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y_sdr] + tsc); 
	      }
	    }
	  }
	}

	if(do_L_v && (do_L_y || do_J_y)) { 
	  /* The second MP_st 'for (jp_v...' loop is for L matrix, which uses a different set of j values from J and R */
	  /* j must satisfy:
	   * j >= jmin[v]
	   * j >= jmin[y] (follows from (j >= jmin[y]))
	   * j <= jmax[v]
	   * j <= jmax[y] (follows from (j <= jmax[y]))
	   * this reduces to two ESL_MAX calls
	   */
	  jn = ESL_MAX(jmin[v], jmin[y]);
	  jx = ESL_MIN(jmax[v], jmax[y]);
	  jpn = jn - jmin[v];
	  jpx = jx - jmin[v];
	  jp_y = jn - jmin[y];
	  /* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	  
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y     >= 0 && jp_y     <= (jmax[y]-jmin[y])));
	    
	    /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	   * this reduces to two ESL_MAX calls
	   */
	    dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y] + sdr);
	    dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y] + sdr);
	    dpn       = dn - hdmin[v][jp_v];
	    dpx       = dx - hdmin[v][jp_v];
	    dp_y_sdr  = dn - hdmin[y][jp_y] - sdr;
	    /* for Lalpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	    
	    for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
	      /* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
	      ESL_DASSERT1((dp_y_sdr >= 0 && dp_y_sdr  <= (hdmax[y][jp_y]     - hdmin[y][jp_y])));
	      if(do_J_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Jalpha[y][jp_y][dp_y_sdr] + tsc);
	      if(do_L_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y][dp_y_sdr] + tsc);
	    }
	  }
	}
      }
      /* add in emission score */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	i     = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v], dp_v = 0; d <= hdmax[v][jp_v]; d++, dp_v++) 
	  {
	    /*if(i < i0 || j > j0) { 
	      printf("dsq[i:%d]: %d\n", i, dsq[i]);
	      printf("dsq[j:%d]: %d\n", j, dsq[j]);
	      printf("esc_v[%d]: %.5f\n", dsq[i]*cm->abc->Kp+dsq[j], esc_v[dsq[i]*cm->abc->Kp+dsq[j]]);;
	      printf("i0: %d j0: %d\n", i0, j0);
	      }*/
	    if(d >= 2) { 
	      if(do_J_v) Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	      if(do_L_v) Lalpha[v][jp_v][dp_v] += lmesc_v[dsq[i]];
	      if(do_R_v) Ralpha[v][jp_v][dp_v] += rmesc_v[dsq[j]];
	    }
	    else { 
	      if(do_J_v) Jalpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      if(do_L_v) Lalpha[v][jp_v][dp_v] = lmesc_v[dsq[i]];
	      if(do_R_v) Ralpha[v][jp_v][dp_v] = rmesc_v[dsq[j]];
	    }
	    i--;
	  }
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	  if(do_J_v) Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  if(do_L_v) Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  if(do_R_v) Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is D or S (! E && ! B && ! ML && ! IL && ! MR && ! IR) */
      /* D, S states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	do_J_y = cp9b->do_J[y];
	do_L_y = cp9b->do_L[y];
	do_R_y = cp9b->do_R[y];
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	if((do_J_v && do_J_y) || (do_L_v && do_L_y) || (do_R_v && do_R_y)) { 
	  /* j must satisfy:
	   * j >= jmin[v]
	   * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	   * j <= jmax[v]
	   * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	   * this reduces to two ESL_MAX calls
	   */
	  jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	  jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	  jpn = jn - jmin[v];
	  jpx = jx - jmin[v];
	  jp_y_sdr = jn - jmin[y] - sdr;
	  
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	    
	    /* d must satisfy:
	     * d >= hdmin[v][jp_v]
	     * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	     * d <= hdmax[v][jp_v]
	     * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	     * this reduces to two ESL_MAX calls
	     */
	    dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	    dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	    dpn     = dn - hdmin[v][jp_v];
	    dpx     = dx - hdmin[v][jp_v];
	    dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	    
	    for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	      if(do_J_v && do_J_y) Jalpha[v][jp_v][dp_v] = FLogsum(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	      if(do_L_v && do_L_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	      if(do_R_v && do_R_y) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y_sd] + tsc);

	      /* an easy to overlook case: if d == 0, set L and R values to IMPOSSIBLE */
	      if(dp_v == dpn && dn == 0) { /* d is 0 */
		if(do_L_v) Lalpha[v][jp_v][dp_v] = IMPOSSIBLE;
		if(do_R_v) Ralpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }		
	    }
	  }
	}
      }
      /* no emission score to add */
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */

      do_J_y = cp9b->do_J[y];
      do_L_y = cp9b->do_L[y];
      do_R_y = cp9b->do_R[y];
      do_T_y = cp9b->do_T[y]; /* will be FALSE, y is not a B_st */

      do_J_z = cp9b->do_J[z];
      do_L_z = cp9b->do_L[z];
      do_R_z = cp9b->do_R[z];
      do_T_z = cp9b->do_T[z]; /* will be FALSE, z is not a B_st */
      
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
	i = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++, i--) {
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	  /* Find the first k value that implies a valid cell in the {J,L,R} matrix y and z decks.
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
	   *
	   * To update a cell in the T matrix with a sum of an R matrix value for y
	   * and a L matrix value for z, there are 2 additional inequalities to satisfy:
	   * (7) k != 0
	   * (8) k != d
	   * We ensure 7 and 8 in the loop below.
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
	      if(do_J_v && do_J_y && do_J_z) Jalpha[v][jp_v][dp_v] = FLogsum(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]);
	      if(do_L_v && do_J_y && do_L_z) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Jalpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]);
	      if(do_R_v && do_R_y && do_J_z) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]);
	      if((k != 0) && (k != d)) {
		if(do_T_v && do_R_y && do_L_z) Talpha[v][jp_v][dp_v] = FLogsum(Talpha[v][jp_v][dp_v], Ralpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]);
	      }
	    }
	  }
	}
      }

      /* two additional special cases in trCYK (these are not in standard CYK).
       * we do these in their own for(j.. { for(d.. { } } loops b/c one 
       * is independent of z, the other of y, unlike the above loop which is dependent 
       * on both.
       */
      if(do_L_v && (do_J_y || do_L_y)) { 
	jn = (jmin[v] > jmin[y]) ? jmin[v] : jmin[y];
	jx = (jmax[v] < jmax[y]) ? jmax[v] : jmax[y];
	for (j = jn; j <= jx; j++) { 
	  jp_v = j - jmin[v];
	  jp_y = j - jmin[y];
	  ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	  ESL_DASSERT1((j >= jmin[y] && j <= jmax[y]));
	  dn = (hdmin[v][jp_v] > hdmin[y][jp_y]) ? hdmin[v][jp_v] : hdmin[y][jp_y];
	  dx = (hdmax[v][jp_v] < hdmax[y][jp_y]) ? hdmax[v][jp_v] : hdmax[y][jp_y];
	  for(d = dn; d <= dx; d++) { 
	    dp_v = d - hdmin[v][jp_v];
	    dp_y = d - hdmin[y][jp_y];
	    ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	    ESL_DASSERT1((d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y]));
	    if(do_J_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Jalpha[y][jp_y][dp_y]);
	    if(do_L_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y][dp_y]);
	  }
	}
      }
      if(do_R_v && (do_J_z || do_R_z)) { 
	jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
	for (j = jn; j <= jx; j++) { 
	  jp_v = j - jmin[v];
	  jp_z = j - jmin[z];
	  ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	  ESL_DASSERT1((j >= jmin[z] && j <= jmax[z]));
	  dn = (hdmin[v][jp_v] > hdmin[z][jp_z]) ? hdmin[v][jp_v] : hdmin[z][jp_z];
	  dx = (hdmax[v][jp_v] < hdmax[z][jp_z]) ? hdmax[v][jp_v] : hdmax[z][jp_z];
	  for(d = dn; d <= dx; d++) { 
	    dp_v = d - hdmin[v][jp_v];
	    dp_z = d - hdmin[z][jp_z];
	    ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	    ESL_DASSERT1((d >= hdmin[z][jp_z] && d <= hdmax[z][jp_z]));
	    if(do_J_z) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Jalpha[z][jp_z][dp_z]);
	    if(do_R_z) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Ralpha[z][jp_z][dp_z]);
	  }
	}
      }
    } /* finished calculating deck v. */
           
    /* allow for truncated begins into v that emit the full sequence i0..j0 */
    if(j0 >= jmin[v] && j0 <= jmax[v]) { 
      jp_v = j0 - jmin[v];
      Wp = W - hdmin[v][jp_v];
      if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	/* If we get here alpha[v][jp_v][Wp] is a valid cell
	 * in the banded alpha matrix, corresponding to 
	 * alpha[v][j0][W] in the platonic matrix.
	 * Check for truncated alignment getting us to the root.
	 */
#if 0
	/* Don't allow local begins, right? */
	/* check for normal local begins */
	if(do_J_v && (cm->flags & CMH_LOCAL_BEGIN) && allow_begin) {
	  bsc = FLogsum(bsc, (Jalpha[v][jp_v][Wp] + cm->beginsc[v]));
	}
#endif
	/* include hits in J matrix (much like a normal local begin) */
	if(do_J_v && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st)) { 
	  bsc = FLogsum(bsc, (Jalpha[v][jp_v][Wp] + trunc_penalty));
	}
	/* include truncated hits in L matrix */
	if(do_L_v && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st)) { 
	  bsc = FLogsum(bsc, (Lalpha[v][jp_v][Wp] + trunc_penalty));
	}	    
	/* check for truncated hit in R matrix */
	if(do_R_v && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st)) { 
	  bsc = FLogsum(bsc, (Ralpha[v][jp_v][Wp] + trunc_penalty));
	}	    
	/* check for truncated hit in T matrix */
	if(do_T_v) { 
	  bsc = FLogsum(bsc, (Talpha[v][jp_v][Wp] + trunc_penalty));
	}	    
      }
      /* include bsc into alpha[0][jp_v][Wp] (full sequence aligned to ROOT_S) */
      if (v == 0) { 
	assert(do_J_v);
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp   = W - hdmin[v][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    Jalpha[0][jp_v][Wp] = FLogsum(Jalpha[0][jp_v][Wp], bsc);
	  }
	}
      } 
    }
  } /* end of for (v = cm->M-1; v > 0; v--) */
  /*FILE *fp1; fp1 = fopen("tmp.iamx", "w");   cm_tr_hb_mx_Dump(fp1, mx); fclose(fp1);*/
  
  Wp = W - hdmin[0][j0-jmin[0]];
  sc =    Jalpha[0][j0-jmin[0]][Wp]; /* this will be bsc, unless a non-truncated hit rooted at 0 is optimal */

  free(el_scA);
  free(yvalidA);

  *ret_sc = sc;

  ESL_DPRINTF1(("TrInsideAlignHB() return sc: %f\n", sc));
  printf("TrInsideAlignHB() return sc: %.4f\n", sc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: TrInsideAlign()
 * Date:     EPN, Mon Sep 12 04:31:43 2011
 *
 * Purpose: Run the truncated inside algorithm on a target sequence.
 *          Uses float log odds scores. Very similar to TrInsideAlign()
 *          except doesn't use HMM bands.
 * 
 *           Very similar to tr_cyk_align(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - we do Inside, not CYK
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *
 *           This function complements TrOutsideAlign().
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score
 * 
 * Returns:  <ret_sc>
 *
 * Throws:  <eslOK> on success
 *          <eslERANGE> if required CM_TR_MX for exceeds <size_limit>, 
 *                      in this case, alignment has been aborted, ret_sc is not valid
 */
int
TrInsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_TR_MX *mx, float *ret_sc)
{
  int      status;          /* easel status code */
  int      v,y,z;	    /* indices for states  */
  int      j,d,i,k;	    /* indices in sequence dimensions */
  int      dp;              /* temporary d value */
  float    sc;		    /* a temporary variable holding a score */
  int      yoffset;	    /* y=base+offset -- counter in child states that v can transit to */
  int      W;	    	    /* subsequence length */
  int      b;	  	    /* best local begin state */
  float    bsc;   	    /* score for using the best local begin state */
  float   *el_scA;          /* [0..d..W-1] probability of local end emissions of length d */
  int      sd;              /* StateDelta(cm->sttype[v]) */
  int      sdl;             /* StateLeftDelta(cm->sttype[v] */
  int      sdr;             /* StateRightDelta(cm->sttype[v] */
  int      jp;              /* offset j, j = i0-1+jp */
  int      j_sdr;           /* j - sdr */
  int      d_sd;            /* d - sd */
  int      d_sdl;           /* d - sdl */
  int      d_sdr;           /* d - sdr */
  float    tsc;             /* a transition score */
  int      L = j0-i0+1;     /* length of the sequence */
  int      have_el;         /* TRUE if local ends are on */
  /* other variables used in truncated version, but not standard version (not in fast_cyk_align()) */
  float Lsc, Rsc;           /* temporary scores */
  float trunc_penalty = 0.; /* penalty in bits for a truncated hit */
  int   bmode = TRMODE_J;   /* mode TRMODE_J, TRMODE_L, TRMODE_R, or TRMODE_T truncation mode for obtaining bsc */

  /* Contract check */
  if(dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "fast_cyk_inside_align(), dsq is NULL.\n");

  /* the DP matrix */
  float ***Jalpha  = mx->Jdp; /* pointer to the Jalpha DP matrix */
  float ***Lalpha  = mx->Ldp; /* pointer to the Lalpha DP matrix */
  float ***Ralpha  = mx->Rdp; /* pointer to the Ralpha DP matrix */
  float ***Talpha  = mx->Tdp; /* pointer to the Talpha DP matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo       (cm, mx,   errbuf, W, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  ESL_STOPWATCH *w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " Matrix init CPU time: ");

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    float const *lmesc_v = cm->lmesc[v]; /* marginal left  emission scores for state v */
    float const *rmesc_v = cm->rmesc[v]; /* marginal right emission scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdl  = StateLeftDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);

    /* re-initialize the J deck if we can do a local end from v */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) { 
	  dp = ESL_MAX(d-sd, 0);
	  Jalpha[v][j][d] = el_scA[dp] + cm->endsc[v];
	  /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (jp = 0; jp <= W; jp++) {
	j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
	Jalpha[v][j][0] = 0.;
	Lalpha[v][j][0] = 0.;
	Ralpha[v][j][0] = 0.;
	for (d = 1; d <= jp; d++) { 
	  Jalpha[v][j][d] = IMPOSSIBLE;
	  Lalpha[v][j][d] = IMPOSSIBLE;
	  Ralpha[v][j][d] = IMPOSSIBLE;
	}
      }
    }
    else if(cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
      /* update alpha[v][j][d] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] 
       * We do ML states as well as IL states b/c they follow the same rules, 
       * and we're not worried about efficiency here.
       */
      
      /* In TrCYK: we need to treat R differently from and J and L
       * here, by doing separate 'for (yoffset...' loops for J and R
       * because we have to fully calculate Jalpha[v][j][d]) before we
       * can start to calculate Ralpha[v][j][d].
       */

      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i    = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	  }
	  Jalpha[v][j][d] += esc_v[dsq[i]];
	  Lalpha[v][j][d]  = (d >= 2) ? Lalpha[v][j][d] + esc_v[dsq[i]] : esc_v[dsq[i]];

	  Jalpha[v][j][d]  = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	  Lalpha[v][j][d]  = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	  i--;

	  /* handle R separately */
	  /* note we use 'd', not 'd_sd' (which we used in the corresponding loop for J,L above) */
	  Rsc = Ralpha[v][j][d]; /* this sc will be IMPOSSIBLE */
	  /* impt to use Rsc because Ralpha[v][j][d] is one of the possible child states we transit to! (if IL) */
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    Rsc = FLogsum(Rsc, Jalpha[y][j_sdr][d] + tsc_v[yoffset]);
	    Rsc = FLogsum(Rsc, Ralpha[y][j_sdr][d] + tsc_v[yoffset]);
	  }
	  Ralpha[v][j][d] = ESL_MAX(Rsc, IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) { 
      /* update alpha[v][j][d] cells, for IR states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1].
       * We do MR states as well as IR states b/c they follow the same rules, 
       * and we're not worried about efficiency here.
       */

      /* In TrCYK: we need to treat L differently from and J and R
       * here, by doing separate 'for (yoffset...' loops for J and R
       * because we have to fully calculate Jalpha[v][j][d]) before we
       * can start to calculate Lalpha[v][j][d].
       */
      for (jp = sdr; jp <= W; jp++) {
	j = i0-1+jp;
	j_sdr = j - sdr;
	for (d = sd; d <= jp; d++) {
	  d_sd = d - sd;
	  i = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	    Ralpha[v][j][d] = FLogsum(Lalpha[v][j][d], Ralpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	  }

	  Jalpha[v][j][d] += esc_v[dsq[i]];
	  Ralpha[v][j][d]  = (d >= 2) ? Ralpha[v][j][d] + esc_v[dsq[j]] : esc_v[dsq[j]];

	  Jalpha[v][j][d]  = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	  Ralpha[v][j][d]  = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);

	  /* handle L separately */
	  /* note we use 'j' and 'd', not 'j_sdr' and 'd_sd' (which we used in the corresponding loop for J,R above) */
	  Lsc = Lalpha[v][j][d]; /* this sc will be IMPOSSIBLE */
	  /* impt to use Lsc because Lalpha[v][j][d] is one of the possible child states we transit to! (if IR) */
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    Lsc = FLogsum(Lsc, Jalpha[y][j][d] + tsc_v[yoffset]);
	    Lsc = FLogsum(Lsc, Lalpha[y][j][d] + tsc_v[yoffset]);
	  }
	  Lalpha[v][j][d] = ESL_MAX(Lsc, IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == MP_st) { 
      /* MP states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];

	for (jp = sdr; jp <= W; jp++) {
	  j = i0-1+jp;
	  j_sdr = j - sdr;

	  for (d = sd; d <= jp; d++) { /* sd == 2 for MP state */
	    d_sd = d - sd;
	    Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	  }
	  /* note we use 'd-sdl' not 'd-sd' for L, plus minimum d is sdl (1) */
	  for (d = sdl; d <= jp; d++) { /* sdl == 1 for MP state */
	    d_sdl = d-sdl;
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Jalpha[y][j_sdr][d_sdl] + tsc_v[yoffset]);
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[y][j_sdr][d_sdl] + tsc_v[yoffset]);
	  }
	  /* note we use 'j' and 'd_sdr' not 'j_sdr' and 'd_sd' for R, plus minimum d is sdr (1) */
	  for (d = sdr; d <= jp; d++) { /* sdr == 1 for MP state */
	    d_sdr = d - sdr;
	    Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Jalpha[y][j][d_sdr] + tsc_v[yoffset]);
	    Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[y][j][d_sdr] + tsc_v[yoffset]);
	  }
	}
      }
      /* add in emission score */
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	i = j;
	Jalpha[v][j][1] = IMPOSSIBLE;
	Lalpha[v][j][1] = lmesc_v[dsq[i]];
	Ralpha[v][j][1] = rmesc_v[dsq[j]];
	i--;
	for (d = 2; d <= jp; d++) {
	  Jalpha[v][j][d] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	  Lalpha[v][j][d] += lmesc_v[dsq[i]];
	  Ralpha[v][j][d] += rmesc_v[dsq[j]];
	  i--;
	}
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 1; d <= jp; d++) {
	  Jalpha[v][j][d] = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	  Lalpha[v][j][d] = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	  Ralpha[v][j][d] = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is D or S */
      /* D, S states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	for (jp = sdr; jp <= W; jp++) {
	  j = i0-1+jp;
	  j_sdr = j - sdr;
	  
	  for (d = sd; d <= jp; d++) {
	    d_sd = d-sd;
	    Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j_sdr][d_sd] + tsc);
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[y][j_sdr][d_sd] + tsc);
	    Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[y][j_sdr][d_sd] + tsc);
	  }
	  /* an easy to overlook case: if d == 0, ensure L and R values are IMPOSSIBLE */
	  Lalpha[v][j][0]  = IMPOSSIBLE;
	  Ralpha[v][j][0] = IMPOSSIBLE;
	}
      }
      /* no emission score to add */
    }
    else { /* B_st */
      assert(cm->sttype[v] == B_st);
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */

      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  for (k = 1; k < d; k++) {
	    Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j-k][d-k] + Jalpha[z][j][k]);
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Jalpha[y][j-k][d-k] + Lalpha[z][j][k]);
	    Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[y][j-k][d-k] + Jalpha[z][j][k]);
	    /*if((k != i-1) && (k != j)) {*/
	    Talpha[v][j][d] = FLogsum(Talpha[v][j][d], Ralpha[y][j-k][d-k] + Lalpha[z][j][k]);
	    /*}*/
	  }
	  /* two additional special cases in trCYK (these are not in standard CYK) */
	  /* special case 1: k == 0 (full sequence aligns to BEGL_S left child */
	  Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Jalpha[v][j][d]);
	  Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[v][j][d]);
	  /* special case 2: k == d (full sequence aligns to BEGR_S right child */
	  Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Jalpha[z][j][d]);
	  Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[z][j][d]);
	}
      }
    } /* end of B_st recursion */

    /* allow for truncated begins into v that emit the full sequence i0..j0 */
#if 0
    /* Don't allow local begins, right? */
    /* check for normal local begins */
    if((cm->flags & CMH_LOCAL_BEGIN) && allow_begin) {
      bsc = FLogsum(bsc, (Jalpha[v][j0][W] + cm->beginsc[v]));
    }
#endif
    /* include hits in J matrix (much like a normal local begin) */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st) { 
      bsc = FLogsum(bsc, Jalpha[v][j0][L] + trunc_penalty);
    }
    /* include truncated hits in L matrix */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st) { 
      bsc = FLogsum(bsc, Lalpha[v][j0][L] + trunc_penalty);
    }	    
    /* check for truncated hit in R matrix */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st) { 
      bsc = FLogsum(bsc, Ralpha[v][j0][L] + trunc_penalty);
    }	    
    /* check for truncated hit in T matrix */
    if(cm->sttype[v] == B_st) { 
      bsc = FLogsum(bsc, Talpha[v][j0][L] + trunc_penalty);
    }	    
    /* include bsc into alpha[0][j0][L] (full sequence aligned to ROOT_S) */
    if (v == 0) { 
      Jalpha[0][j0][L] = FLogsum(Jalpha[0][j0][L], bsc);
    }
  } /* end of for (v = cm->M-1; v > 0; v--) */
  /*FILE *fp1; fp1 = fopen("tmp.amx", "w");   cm_tr_mx_Dump(fp1, mx); fclose(fp1);*/
  
  sc = Jalpha[0][j0][L]; 

  free(el_scA);

  *ret_sc = sc;

  ESL_DPRINTF1(("TrInsideAlign() return sc: %f\n", sc));
  printf("TrInsideAlign() return sc: %.4f\n", sc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: TrOutsideAlign()
 * Date:     EPN, Mon Sep 12 05:32:44 2011
 *
 * Purpose:  Run the truncated outside algorithm. Non-banded version.
 *           A CM_TR_MX DP matrix must be passed in. 
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align  (L, for whole seq)
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_check  - TRUE to attempt to check matrices for correctness
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ins_mx    - the dp matrix from the Inside run calculation (required)
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *                       ends are on (see *** comment towards end of function).
 *
 * Returns:  <ret_sc>
 *
 * Throws:  <eslOK> on success
 *          <eslERANGE> if required CM_TR_MX for exceeds <size_limit>, 
 *                      in this case, alignment has been aborted, ret_sc is not valid *                       
 */
int
TrOutsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, int do_check, 
	       CM_TR_MX *mx, CM_TR_MX *ins_mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float    fsc;     	       /* a temporary variable holding a float score */
  float    escore;	       /* an emission score, tmp variable */
  int      W;		       /* subsequence length */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      jp;		       /* j': relative position in the subsequence  */
  float    bsc;		       /* total score for using local begin states */
  float    freturn_sc;         /* P(S|M)/P(S|R), a float (Scorified ireturn_sc) */
  int      fail_flag = FALSE;  /* set to TRUE if do_check and we see a problem */
  /* variables used only if do_check */
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */
  /* indices used in the depths of the DP recursion */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdl;                /* StateLeftDelta(cm->sttype[y] */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */
  int      ip_L;               /* i index in L matrix */

  /* Contract check */
  if (dsq == NULL)                                     ESL_FAIL(eslEINCOMPAT, errbuf, "TrOutsideAlign(), dsq is NULL.\n");
  if (mx == NULL)                                      ESL_FAIL(eslEINCOMPAT, errbuf, "TrOutsideAlign(), mx is NULL.\n");
  if (ins_mx == NULL)                                  ESL_FAIL(eslEINCOMPAT, errbuf, "TrOutsideAlign(), ins_mx is NULL.\n");
  if (cm->flags & CMH_LOCAL_END) do_check = FALSE; /* Code for checking doesn't apply in local mode. See below. */

  /* DP matrix variables */
  float ***Jbeta   = mx->Jdp;     /* pointer to the outside Jbeta DP matrix */
  float ***Lbeta   = mx->Ldp;     /* pointer to the outside Lbeta DP matrix */
  float ***Rbeta   = mx->Rdp;     /* pointer to the outside Rbeta DP matrix */
  float ***Tbeta   = mx->Tdp;     /* pointer to the outside Tbeta DP matrix */

  float ***Jalpha  = ins_mx->Jdp; /* pointer to the precalc'ed inside Jalpha DP matrix */
  float ***Lalpha  = ins_mx->Ldp; /* pointer to the precalc'ed inside Lalpha DP matrix */
  float ***Ralpha  = ins_mx->Rdp; /* pointer to the precalc'ed inside Ralpha DP matrix */
  float ***Talpha  = ins_mx->Tdp; /* pointer to the precalc'ed inside Talpha DP matrix */

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;              /* the summed prob of all local begins */
  W   = j0-i0+1;		 /* the length of the subsequence -- used in many loops  */
				 /* if caller didn't give us a deck pool, make one */

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo(cm, mx, errbuf, W, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 

  /* now set Jbeta[0][j0][W] to 0., all J parses must end there */
  Jbeta[0][j0][W] = 0.;

  /* for any MP, ML, MR, B state, any parse can also end for any j and d, these are truncated ends */
  for(v = 1; v < cm->M; v++) { 
    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == B_st) { 
      for(jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) {
	  Lalpha[v][j][d] = 0.;
	  Ralpha[v][j][d] = 0.;
	}
      }
    }
  }

  /* No local begins in truncated mode */
  #if 0 
  /* If we can do a local begin into v, overwrite IMPOSSIBLE with the local begin score. 
   * By definition, beta[0][j0][W] == 0.
   */ 
  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (v = 1; v < cm->M; v++) {
      if(NOT_IMPOSSIBLE(cm->beginsc[v])) {
	Jbeta[v][j0][W] = cm->beginsc[v];
      }
    }
  }
#endif

  /* main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    sd  = StateDelta(cm->sttype[v]);
    sdr = StateRightDelta(cm->sttype[v]);

    /* mini-recursion for L matrix: because for any cell in L all
     * 'outside' residues must be emitted leftwise, we don't need to
     * store varying d values, i.e. j defines i and j because (for L)
     * j can never change again. Put another way, we are now behaving
     * like an HMM.  This means we only have to calculate and store L
     * matrix values for a single d instead of for all d 0..j.
     */
    if (cm->stid[v] == BEGL_S) { 
      cm_Fail("Outside BEGL_S not yet implemented");
    }
    else if(cm->stid[v] == BEGR_S) { 
      cm_Fail("Outside BEGR_S not yet implemented");
    }
    else { /* ! BEGL_S and ! BEGR_S */
      for(jp = 0; jp <= W; jp++) { 
	j = i0+jp; /* note off-by-one relative to mini-recursion for R matrix */
	if (j > i0) { 
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    sd  = StateDelta(cm->sttype[y]);
	    sdl = StateLeftDelta(cm->sttype[y]);
	    if     (cm->sttype[y] == MP_st)                           escore = cm->lmesc[y][dsq[j-1]];
	    else if(cm->sttype[y] == ML_st || cm->sttype[y] == IL_st) escore = cm->oesc[y][dsq[j-1]];
	    else                                                      escore = 0.;

	    Lbeta[v][jp][0] = FLogsum(Lbeta[v][jp][0], (Lbeta[y][jp-sdl][0] + cm->tsc[y][voffset] + escore));
	    /* sdl is 1 for MP, ML, IL, else it's 0 */
	  }
	  Lbeta[v][jp][0] = ESL_MAX(Lbeta[v][j][0], IMPOSSIBLE);
	}
      }
    }
    /* mini-recursion for R matrix: because for any cell in R all
     * 'outside' residues must be emitted rightwise, we don't need to
     * store varying d values, i.e. j defines i and j because (for R)
     * i can never change again. Put another way, we are now behaving
     * like an HMM.  This means we only have to calculate and store R
     * matrix values for a single d instead of for all d 0..j.
     */
    if (cm->stid[v] == BEGL_S) { 
      cm_Fail("Outside BEGL_S not yet implemented");
    }
    else if(cm->stid[v] == BEGR_S) { 
      cm_Fail("Outside BEGR_S not yet implemented");
    }
    else { /* ! BEGL_S and ! BEGR_S */
      for (jp = W; jp >= 0; jp--) {
	j = i0-1+jp; /* note off-by-one relative to mini-recursion for L matrix */
	if(j < j0) { 
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    sd  = StateDelta(cm->sttype[y]);
	    sdr = StateRightDelta(cm->sttype[y]);
	    if     (cm->sttype[y] == MP_st)                           escore = cm->rmesc[y][dsq[j+1]];
	    else if(cm->sttype[y] == MR_st || cm->sttype[y] == IR_st) escore = cm->oesc[y][dsq[j+1]];
	    else                                                      escore = 0.;
	    Rbeta[v][j][0] = FLogsum(Rbeta[v][j][0], (Rbeta[y][j+sdr][0] + cm->tsc[y][voffset] + escore));
	    /* sdr is 1 for MP, MR, IR, else it's 0 */
	  }
	  Rbeta[v][j][0] = ESL_MAX(Rbeta[v][j][0], IMPOSSIBLE);
	}
      }
    }

    /* main recursion, for J matrix, just like Inside/CYK, we need to do for(j ... for(d ... here */
    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      cm_Fail("Outside BEGL_S not yet implemented");
#if 0
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
#endif
    } /* end of 'if (cm->stid[v] == BEGL_S */
    else if (cm->stid[v] == BEGR_S) {
      cm_Fail("Outside BEGR_S not yet implemented");
#if 0
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
#endif
    } /* end of 'else if (cm->stid[v] == BEGR_S */
    else { /* (cm->sttype[v] != BEGL_S && cm->sttype[v] != BEGR_S */ 
      for (jp = W; jp >= 0; jp--) {
	j    = i0-1+jp;
	i    = j-jp+1;
	ip_L = i - 1; /* offset in L matrix */
	for (d = jp; d >= 0; d--, i++, ip_L++) {
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    sd  = StateDelta(cm->sttype[y]);
	    sdl = StateLeftDelta(cm->sttype[y]);
	    sdr = StateRightDelta(cm->sttype[y]);
	    switch(cm->sttype[y]) {
	    case MP_st: 
	      if (j != j0 && d != jp) { 
		escore = cm->oesc[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j-sdl][0]    + cm->tsc[y][voffset] + escore);
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Lbeta[y][j+sdr][0]    + cm->tsc[y][voffset] + escore);
	      }
	      break;
	    case ML_st:
	    case IL_st: 
	      if (d != jp) { 
		escore = cm->oesc[y][dsq[i-1]];
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j][d+sd]     + cm->tsc[y][voffset] + escore);
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Lbeta[y][ip_L-sdl][0] + cm->tsc[y][voffset] + escore);
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Rbeta[y][j][0]        + cm->tsc[y][voffset] + escore);
	      }
	      break;
	    case MR_st:
	    case IR_st:
	      if (j != j0) { 
		escore = cm->oesc[y][dsq[j+1]];
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Lbeta[y][ip_L][0]     + cm->tsc[y][voffset] + escore);
		Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Rbeta[y][j+sdr][0]    + cm->tsc[y][voffset] + escore);
	      }
	      break;
	    case S_st:
	    case E_st:
	    case D_st:
	      Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j][d] + cm->tsc[y][voffset]);
	      break;
	    } /* end of switch(cm->sttype[y] */  
	  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	  if (Jbeta[v][j][d] < IMPOSSIBLE) Jbeta[v][j][d] = IMPOSSIBLE;
	} /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
      } /* end loop over jp. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] != BEGL_S, BEGR_S */
    /* we're done calculating deck v for everything but local ends */

    /* deal with local alignment end transitions v->EL J matrix only (EL = deck at M.) */
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
	    escore = cm->oesc[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == jp) continue;	
	    escore = cm->oesc[v][dsq[i-1]];
	    Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == j0) continue;
	    escore = cm->oesc[v][dsq[j+1]];
	    Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case S_st:
	  case D_st:
	  case E_st:
	    Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v]));
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
	Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], Jbeta[cm->M][j][d+1] + cm->el_selfsc);
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
#if 0
  if(!(cm->flags & CMH_LOCAL_END)) { 
    freturn_sc = IMPOSSIBLE;
    v = cm->M-1;
    for(jp = 0; jp <= W; jp++) { 
      j = i0+jp; 
      printf("\tJalpha[%3d][%3d][%3d]: %5.2f | Jbeta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, Jalpha[(cm->M-1)][jp][0], (cm->M-1), (j), 0, Jbeta[(cm->M-1)][jp][0]);
      printf("\tLalpha[%3d][%3d][%3d]: %5.2f | Lbeta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, Lalpha[(cm->M-1)][jp][0], (cm->M-1), (j), 0, Lbeta[(cm->M-1)][jp][0]);
      printf("\tRalpha[%3d][%3d][%3d]: %5.2f | Rbeta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, Ralpha[(cm->M-1)][jp][0], (cm->M-1), (j), 0, Rbeta[(cm->M-1)][jp][0]);
      freturn_sc = FLogsum(freturn_sc, (Jbeta[v][jp][0]));
      freturn_sc = FLogsum(freturn_sc, (Lbeta[v][jp][0]));
      freturn_sc = FLogsum(freturn_sc, (Rbeta[v][jp][0]));
    }
  }
  else { /* return_sc = P(S|M) / P(S|R) from Inside() */
    freturn_sc = Jalpha[0][(j0-jmin[0])][Wp];
  }
#endif
  if(fail_flag) ESL_FAIL(eslFAIL, errbuf, "Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END)) ESL_DPRINTF1(("\tTrOutsideAlign() sc : %f\n", freturn_sc));
  else                             ESL_DPRINTF1(("\tTrOutsideAlign() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  if (ret_sc != NULL) *ret_sc = freturn_sc;
  return eslOK;
}  


/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_TRUNC_ALIGN_BENCHMARK
/* Next line is optimized (debugging on) on MacBook Pro:
 * gcc   -o benchmark-trunc-align -std=gnu99 -g -Wall -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DIMPL_TRUNC_ALIGN_BENCHMARK cm_dpalign_trunc.c -linfernal -lhmmer -leasel -lm
 * Next line is optimized (debugging not on) on wyvern:
 * gcc   -o benchmark-trunc-align -std=gnu99 -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -pthread -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DIMPL_TRUNC_ALIGN_BENCHMARK cm_dpalign_trunc.c -linfernal -lhmmer -leasel -lm 
 * ./benchmark-trunc-align <cmfile>
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_randomseq.h>
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
  { "-s",        eslARG_INT,    "181", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>, '0' for one-time arbitrary", 0 },
  { "-e",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "emit sequences from CM, don't randomly create them", 0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search in glocal mode [default: local]", 0 },
  { "-L",        eslARG_INT,     NULL, NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs, default: consensus length", 0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of target seqs",                          0 },
  { "--orig",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also do search with original trCYK",             0},
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "don't use QDBs", 0},
  { "--infile",  eslARG_INFILE,  NULL, NULL, NULL,  NULL,  NULL, "-L,-N,-e", "read sequences to search from file <s>", 2 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use posterior sums during HMM band calculation (widens bands)", 0 },
  { "--onlyhb",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only run HMM banded scanning trCYK", 2 },
  { "--tau",     eslARG_REAL,   "5e-6",NULL, "0<x<1",NULL, NULL, NULL, "set tail loss prob for HMM bands to <x>", 2 },
  { "--cp9noel", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-g",           "turn OFF local ends in cp9 HMMs", 2 },
  { "--cp9gloc", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-g,--cp9noel", "configure CP9 HMM in glocal mode", 2 },
  { "--thresh1", eslARG_REAL,  "0.01", NULL, NULL,  NULL,  NULL,  NULL, "set HMM bands thresh1 to <x>", 2 },
  { "--thresh2", eslARG_REAL,  "0.99", NULL, NULL,  NULL,  NULL,  NULL, "set HMM bands thresh2 to <x>", 2 },
  { "--sizelimit",eslARG_REAL, "128.", NULL, "x>0", NULL,  NULL,  NULL, "set maximum allowed size of HB matrices to <x> Mb", 2 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "benchmark driver for truncated alignment implementations";

int 
main(int argc, char **argv)
{
  int             status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  CM_t           *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  int64_t         L;
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq;
  int             i;
  float           sc;
  char           *cmfile = esl_opt_GetArg(go, 1);
  CM_FILE        *cmfp;	/* open input CM file stream */
  int            *dmin;
  int            *dmax;
  int             do_random;
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, either randomly created, or emitted from CM (if -e) */
  char            errbuf[cmERRBUFSIZE];
  TrScanMatrix_t *trsmx = NULL;
  ESL_SQFILE     *sqfp  = NULL;        /* open sequence input file stream */
  CMConsensus_t  *cons  = NULL;
  Parsetree_t    *tr    = NULL;
  CM_TR_HB_MX           *trhbmx= NULL;
  CM_TR_HB_SHADOW_MX    *trhbshmx= NULL;
  CM_TR_MX              *trmx= NULL;
  CM_TR_MX              *out_trmx= NULL;
  CM_TR_SHADOW_MX       *trshmx= NULL;
  float           size_limit = esl_opt_GetReal(go, "--sizelimit");
  float           save_tau, save_cp9b_thresh1, save_cp9b_thresh2;
  float           hbmx_Mb, trhbmx_Mb;
  float           bsc;
  float           parsetree_sc, parsetree_struct_sc;
  int             v;
  int             b;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  if ((status = cm_file_Open(cmfile, NULL, FALSE, &(cmfp), errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = cm_file_Read(cmfp, TRUE, &abc, &cm))                != eslOK) cm_Fail(cmfp->errbuf);
  cm_file_Close(cmfp);

  /* determine sequence length */
  if(esl_opt_IsUsed(go, "-L")) L = esl_opt_GetInteger(go, "-L");
  else                         L = cm->clen;      

  /* configure CM for HMM banded alignment */
  cm->align_opts  |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--sums")) cm->align_opts |= CM_ALIGN_SUMS;

  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      if(! esl_opt_GetBoolean(go, "--cp9noel")) cm->config_opts |= CM_CONFIG_HMMEL; 
    }
  }

  esl_stopwatch_Start(w);
  printf("%-30s", "Configuring CM...");
  fflush(stdout);
  ConfigCM(cm, errbuf, FALSE, NULL, NULL); /* FALSE says: don't calculate W */
  printf("done.  ");
  fflush(stdout);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " CPU time: ");

  if (esl_opt_IsUsed(go, "--thresh1")) { cm->cp9b->thresh1 = esl_opt_GetReal(go, "--thresh1"); }
  if (esl_opt_IsUsed(go, "--thresh2")) { cm->cp9b->thresh2 = esl_opt_GetReal(go, "--thresh2"); }

  if (esl_opt_GetBoolean(go, "--noqdb")) { 
    if(cm->dmin != NULL) { free(cm->dmin); cm->dmin = NULL; }
    if(cm->dmax != NULL) { free(cm->dmax); cm->dmax = NULL; }
  }
  dmin = cm->dmin; 
  dmax = cm->dmax; 
  cm->tau = esl_opt_GetReal(go, "--tau");

  printf("%-30s", "Creating tr hb matrix...");
  fflush(stdout);
  esl_stopwatch_Start(w);
  trhbmx   = cm_tr_hb_mx_Create(cm);
  trhbshmx = cm_tr_hb_shadow_mx_Create(cm);
  printf("done.  ");
  fflush(stdout);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " CPU time: ");
  printf("\n\n");

  if(! esl_opt_GetBoolean(go, "--onlyhb")) { 
    printf("%-30s", "Creating tr matrix...");
    fflush(stdout);
    esl_stopwatch_Start(w);
    trmx   = cm_tr_mx_Create(cm);
    trshmx = cm_tr_shadow_mx_Create(cm);
    printf("done.  ");
    fflush(stdout);
    esl_stopwatch_Stop(w);
    esl_stopwatch_Display(stdout, w, " CPU time: ");
    printf("\n\n");

    printf("%-30s", "Creating tr scan matrix...");
    esl_stopwatch_Start(w);
    trsmx = cm_CreateTrScanMatrix(cm, cm->W, dmax, cm->beta_W, cm->beta_qdb, 
				  (dmin == NULL && dmax == NULL) ? FALSE : TRUE,
				  TRUE, TRUE); /* do_float, do_int */
    if(trsmx == NULL) esl_fatal("Problem creating trsmx");
    printf("done.  ");
    fflush(stdout);
    esl_stopwatch_Stop(w);
    esl_stopwatch_Display(stdout, w, " CPU time: ");
  }
  
  /* get sequences */
  if(esl_opt_IsUsed(go, "--infile")) { 
    /* read sequences from a file */
    status = esl_sqfile_OpenDigital(cm->abc, esl_opt_GetString(go, "--infile"), eslSQFILE_UNKNOWN, NULL, &sqfp);
    if (status == eslENOTFOUND)    esl_fatal("File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--infile"));
    else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--infile"));
    else if (status == eslEINVAL)  esl_fatal("Can't autodetect stdin or .gz."); 
    else if (status != eslOK)      esl_fatal("Sequence file open failed with error %d.\n", status);

    seqs_to_aln = CreateSeqsToAln(100, FALSE);
    if((status = ReadSeqsToAln(cm->abc, sqfp, 0, seqs_to_aln, FALSE)) != eslEOF)
      esl_fatal("Error reading sqfile: %s\n", esl_opt_GetString(go, "--infile"));
    esl_sqfile_Close(sqfp);
    N = seqs_to_aln->nseq;
  }
  else if(esl_opt_IsUsed(go, "-L")) {
     double *dnull;
     ESL_DSQ *randdsq = NULL;
     ESL_ALLOC(randdsq, sizeof(ESL_DSQ)* (L+2));
     ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
     for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
     esl_vec_DNorm(dnull, cm->abc->K);
     seqs_to_aln = CreateSeqsToAln(N, FALSE);

     for (i = 0; i < N; i++) {
       if (esl_rsq_xIID(r, dnull, cm->abc->K, L, randdsq)  != eslOK) cm_Fail("Failure creating random sequence.");
       if((seqs_to_aln->sq[i] = esl_sq_CreateDigitalFrom(abc, NULL, randdsq, L, NULL, NULL, NULL)) == NULL)
         cm_Fail("Failure digitizing/copying random sequence.");
     }
  }
  else if(do_random) {
    double *dnull;
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);
    /* get gamma[0] from the QDB calc alg, which will serve as the length distro for random seqs */
    int safe_windowlen = cm->clen * 2;
    double **gamma = NULL;
    while(!(BandCalculationEngine(cm, safe_windowlen, DEFAULT_BETA, TRUE, NULL, NULL, &(gamma), NULL))) {
      safe_windowlen *= 2;
      /* This is a temporary fix becuase BCE() overwrites gamma, leaking memory
       * Probably better long-term for BCE() to check whether gamma is already allocated
       */
      FreeBandDensities(cm, gamma);
      if(safe_windowlen > (cm->clen * 1000)) cm_Fail("Error trying to get gamma[0], safe_windowlen big: %d\n", safe_windowlen);
    }
    seqs_to_aln = RandomEmitSeqsToAln(r, cm->abc, dnull, 1, N, gamma[0], safe_windowlen, FALSE);
    FreeBandDensities(cm, gamma);
    free(dnull);
  }
  else /* don't randomly generate seqs, emit them from the CM */
    seqs_to_aln = CMEmitSeqsToAln(r, cm, 1, N, FALSE, NULL, FALSE);

  CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &cons);
  
  save_tau = cm->tau;
  save_cp9b_thresh1 = cm->cp9b->thresh1;
  save_cp9b_thresh2 = cm->cp9b->thresh2;

  for (i = 0; i < N; i++)
    {
      L = seqs_to_aln->sq[i]->n;
      dsq = seqs_to_aln->sq[i]->dsq;
      cm->search_opts &= ~CM_SEARCH_INSIDE;

      cm->tau = save_tau;
      cm->cp9b->thresh1 = save_cp9b_thresh1;
      cm->cp9b->thresh2 = save_cp9b_thresh2;

      cm->align_opts  |= CM_ALIGN_HBANDED;


      /*********************Begin FastCYKScanHB****************************/
      esl_stopwatch_Start(w);
      cm->tau = save_tau;
      while(1) { 
	if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, 
				   TRUE,  /* doing search? */
				   FALSE,  /* doing trCYK/Inside/Outside? */
				   0)) != eslOK) cm_Fail(errbuf);
	if((status = cm_hb_mx_SizeNeeded(cm, errbuf, cm->cp9b, L, NULL, &hbmx_Mb)) != eslOK) return status; 
	if(hbmx_Mb < size_limit) break; /* our matrix will be small enough, break out of while(1) */
	if(cm->tau > 0.01)         cm_Fail("tau reached limit, unable to create matrix smaller than size limit of %.2f Mb\n", size_limit);
	printf("  CYK 0 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb);
	cm->tau *= 2.;
      }

      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "HMM Band calc:", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, L);
      
      esl_stopwatch_Start(w);
      if((status = FastCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, FALSE, cm->hbmx, size_limit, 0., NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScanHB(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      /*********************End FastCYKScanHB****************************/

      /*********************Begin FastAlignHB()***************************/
      esl_stopwatch_Start(w);
      while(1) { 
	if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, 
				   FALSE,  /* doing search? */
				   FALSE,  /* doing trCYK/Inside/Outside? */
				   0)) != eslOK) cm_Fail(errbuf);
	if((status = cm_hb_mx_SizeNeeded(cm, errbuf, cm->cp9b, L, NULL, &hbmx_Mb)) != eslOK) return status; 
	if(hbmx_Mb < size_limit) break; /* our matrix will be small enough, break out of while(1) */
	if(cm->tau > 0.01)         cm_Fail("tau reached limit, unable to create matrix smaller than size limit of %.2f Mb\n", size_limit);
	printf("  CYK 0 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb);
	cm->tau *= 2.;
      }
      
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "HMM Band calc:", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, L);
      
      esl_stopwatch_Start(w);
      if((status = FastAlignHB(cm, errbuf, NULL, dsq, L, 1, L, size_limit, cm->hbmx, cm->shmx, FALSE, FALSE, NULL, NULL, NULL, &sc, NULL)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastAlignHB(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End FastAlignHB()***************************/


      /*********************Begin TrCYKScanHB****************************/
      esl_stopwatch_Start(w);
      /* Calculate HMM bands. We'll tighten tau and recalculate bands until 
       * the resulting HMM banded matrix is under our size limit. 
      */
      cm->tau = save_tau;
      while(1) { 
	if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, 
				   TRUE,  /* doing search? */
				   TRUE,  /* doing trCYK/Inside/Outside? */
				   0)) != eslOK) cm_Fail(errbuf);
	if((status = cm_tr_hb_mx_SizeNeeded(cm, errbuf, cm->cp9b, L, NULL, NULL, NULL, NULL, &trhbmx_Mb)) != eslOK) return status; 
	if(trhbmx_Mb < size_limit) break; /* our matrix will be small enough, break out of while(1) */
	if(cm->tau > 0.01)         cm_Fail("tau reached limit, unable to create matrix smaller than size limit of %.2f Mb\n", size_limit);
	printf("TrCYK 0 tau: %10g  thresh1: %10g  thresh2: %10g  trhbmx_Mb: %10.2f\n", cm->tau, cm->cp9b->thresh1, cm->cp9b->thresh2, trhbmx_Mb);
	cm->tau *= 2.;
	cm->cp9b->thresh1 *= 2.; 
	cm->cp9b->thresh2 -= (1.0-cm->cp9b->thresh2); 
	cm->cp9b->thresh1 = ESL_MIN(0.25, cm->cp9b->thresh1);
	cm->cp9b->thresh2 = ESL_MAX(0.25, cm->cp9b->thresh2);
      }	  
      printf("TrCYK 1 tau: %10g  thresh1: %10g  thresh2: %10g  trhbmx_Mb: %10.2f\n", cm->tau, cm->cp9b->thresh1, cm->cp9b->thresh2, trhbmx_Mb);
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "HMM Band calc:", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, L);
      
      cm_tr_hb_mx_Destroy(trhbmx);
      trhbmx = cm_tr_hb_mx_Create(cm);
      esl_stopwatch_Start(w);
      if((status = TrCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, FALSE, trhbmx, size_limit, 0.,  NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "TrCYKScanHB(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End TrCYKScanHB****************************/

      /*********************Begin FTrInsideScanHB****************************/
      cm_tr_hb_mx_Destroy(trhbmx);
      trhbmx = cm_tr_hb_mx_Create(cm);
      esl_stopwatch_Start(w);
      if((status = FTrInsideScanHB(cm, errbuf, dsq, 1, L, 0., NULL, FALSE, trhbmx, size_limit, 0.,  NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FTrInsideScanHB(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End TrCYKScanHB****************************/

      /*********************Begin TrAlignHB()****************************/
      esl_stopwatch_Start(w);
      /* Calculate HMM bands. We'll tighten tau and recalculate bands until 
       * the resulting HMM banded matrix is under our size limit.
       */
      cm->tau = save_tau;
      while(1) { 
	if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, 
				   TRUE,  /* doing search? */
				   TRUE,  /* doing trCYK/Inside/Outside? */
				   0)) != eslOK) cm_Fail(errbuf);
	if((status = cm_tr_hb_mx_SizeNeeded(cm, errbuf, cm->cp9b, L, NULL, NULL, NULL, NULL, &trhbmx_Mb)) != eslOK) return status; 
	if(trhbmx_Mb < size_limit) break; /* our matrix will be small enough, break out of while(1) */
	if(cm->tau > 0.01)         cm_Fail("tau reached limit, unable to create matrix smaller than size limit of %.2f Mb\n", size_limit);
	printf("TrCYK 0 tau: %10g  thresh1: %10g  thresh2: %10g  trhbmx_Mb: %10.2f\n", cm->tau, cm->cp9b->thresh1, cm->cp9b->thresh2, trhbmx_Mb);
	cm->tau *= 2.;
	cm->cp9b->thresh1 *= 2.; 
	cm->cp9b->thresh2 -= (1.0-cm->cp9b->thresh2); 
	cm->cp9b->thresh1 = ESL_MIN(0.25, cm->cp9b->thresh1);
	cm->cp9b->thresh2 = ESL_MAX(0.25, cm->cp9b->thresh2);
      }	  
      printf("TrCYK 1 tau: %10g  thresh1: %10g  thresh2: %10g  trhbmx_Mb: %10.2f\n", cm->tau, cm->cp9b->thresh1, cm->cp9b->thresh2, trhbmx_Mb);
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "HMM Band calc:", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, L);
      
      cm_tr_hb_mx_Destroy(trhbmx);
      cm_tr_hb_shadow_mx_Destroy(trhbshmx);
      trhbmx   = cm_tr_hb_mx_Create(cm);
      trhbshmx = cm_tr_hb_shadow_mx_Create(cm);
      esl_stopwatch_Start(w);
      if((status = TrAlignHB(cm, errbuf, NULL, dsq, 1, L, size_limit, trhbmx, trhbshmx, FALSE, FALSE, NULL, &tr, NULL, &sc, NULL)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "TrAlignHB(): ", sc);
      esl_stopwatch_Stop(w);
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
      ParsetreeScore(cm, NULL, NULL, tr, dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
      printf("Parsetree score      : %.4f           (FULL LENGTH CYK)\n", parsetree_sc);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End TrAlignHB*****************************/

      /*********************Begin TrInsideAlignHB()****************************/
      cm_tr_hb_mx_Destroy(trhbmx);
      trhbmx = cm_tr_hb_mx_Create(cm);
      esl_stopwatch_Start(w);
      if((status = TrInsideAlignHB(cm, errbuf, dsq, 1, L, size_limit, trhbmx, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits (FULL LENGTH INSIDE)", (i+1), "TrInsideAlignHB(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End TrInsideAlignHB*****************************/

      if(! esl_opt_GetBoolean(go, "--onlyhb")) { 
	/*********************Begin TrAlign****************************/
	cm_tr_mx_Destroy(trmx);
	cm_tr_shadow_mx_Destroy(trshmx);
	trmx   = cm_tr_mx_Create(cm);
	trshmx = cm_tr_shadow_mx_Create(cm);
	esl_stopwatch_Start(w);
	if((status = TrAlign(cm, errbuf, NULL, dsq, 1, L, size_limit, trmx, trshmx, FALSE, FALSE, NULL, &tr, NULL, &sc, NULL)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits (FULL LENGTH CYK)", (i+1), "TrAlign(): ", sc);
	esl_stopwatch_Stop(w);
	ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
	ParsetreeScore(cm, NULL, NULL, tr, dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
	printf("Parsetree score      : %.4f           (FULL LENGTH CYK)\n", parsetree_sc);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*********************End TrAlign****************************/
      
	/*********************Begin TrInsideAlign()****************************/
	cm_tr_mx_Destroy(trmx);
	trmx   = cm_tr_mx_Create(cm);
	esl_stopwatch_Start(w);
	if((status = TrInsideAlign(cm, errbuf, dsq, 1, L, size_limit, trmx, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits (FULL LENGTH INSIDE)", (i+1), "TrInsideAlign(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*********************End TrInsideAlign*****************************/
      
	/*********************Begin TrOutsideAlign()****************************/
	out_trmx   = cm_tr_mx_Create(cm);
	esl_stopwatch_Start(w);
	if((status = TrOutsideAlign(cm, errbuf, dsq, 1, L, size_limit, TRUE, out_trmx, trmx, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits (FULL LENGTH INSIDE)", (i+1), "TrOutsideAlign(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*********************End TrOutsideAlign*****************************/

	/*********************Begin RefTrCYKScan****************************/
	esl_stopwatch_Start(w);
	if((status = RefTrCYKScan(cm, errbuf, trsmx, dsq, 1, L, 0., NULL, FALSE, 0., NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "RefTrCYKScan(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*********************End RefTrCYKScan****************************/

	/*********************Begin RefITrInsideScan****************************/
	cm->search_opts |= CM_SEARCH_INSIDE;
	esl_stopwatch_Start(w);
	if((status = RefITrInsideScan(cm, errbuf, trsmx, dsq, 1, L, 0., NULL, FALSE, 0., NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "RefITrInsideScan(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	cm->search_opts &= ~CM_SEARCH_INSIDE;
	/*********************End RefITrInsideScan****************************/
	
	if(esl_opt_GetBoolean(go, "--orig")) { 
	  /*********************Begin TrCYK_Inside****************************/
	  esl_stopwatch_Start(w);
	  sc = TrCYK_Inside(cm, dsq, L, 0, 1, L, FALSE, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "TrCYK_Inside():   ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	  /*********************End TrCYK_Inside****************************/
	}
      }
      printf("\n");
    }
  FreeCM(cm);
  FreeSeqsToAln(seqs_to_aln);
  cm_tr_hb_mx_Destroy(trhbmx);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  cm_Fail("memory allocation error");
  return 0; /* never reached */
}
#endif /*IMPL_TRUNC_ALIGN_BENCHMARK*/




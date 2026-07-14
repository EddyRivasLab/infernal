/* cm_dpalign.c
 * 
 * DP functions for standard (non-truncated) HMM banded and
 * non-banded, non-D&C CM alignment of a full target sequence.
 * 
 * All functions use a DP matrix and or shadow matrix, either
 * non-banded (CM_MX, CM_SHADOW_MX) or HMM banded (CM_HB_MX,
 * CM_HB_SHADOW_MX).  The HMM banded matrices only have cells within
 * bands allocated. The bands derived from a HMM Forward/Backward
 * alignment of the target sequence and are stored in a CP9Bands_t
 * object, a pointer to which must exist in the cm (CM_t object).
 * 
 * The non-banded, non-D&C alignment functions are mainly useful for
 * understanding and/or debugging the HMM banded versions.  These are
 * consistent (same logic/code organization) with their HMM banded
 * counterparts. They are memory intensive. For small memory
 * non-banded alignment functions see cm_dpsmall.c. For truncated
 * alignment functions (both non-banded and HMM banded) see
 * cm_dpalign_trunc.c.
 *
 * List of functions: 
 * non-banded version        HMM banded version
 * -----------------------   ------------------------
 * cm_alignT()               cm_alignT_hb()
 * cm_AlignSizeNeeded()      cm_AlignSizeNeededHB()
 * cm_Align()                cm_AlignHB()
 * cm_CYKInsideAlign()       cm_CYKInsideAlignHB()
 * cm_InsideAlign()          cm_InsideAlignHB()
 * cm_OptAccAlign()          cm_OptAccAlignHB()
 * cm_CYKOutsideAlign()*     cm_CYKOutsideAlignHB()*
 * cm_OutsideAlign()         cm_OutsideAlignHB()
 * cm_Posterior()            cm_PosteriorHB()  
 *
 * * cm_CYKOutsideAlign() and cm_CYKOutsideAlignHB() are for reference
 * and debugging only they're not called by any of the main Infernal
 * programs, only by test programs.
 * 
 * EPN, Wed Sep 14 05:31:02 2011 Note: post version 1.0.2, the
 * 'Fast'/'fast_' prefix was dropped from many of these functions and
 * the cm_ prefix was added. Also 'optimal_accuracy' was shortened to
 * 'optacc'. At the same time, CM_MX and CM_SHADOW_MX data structures
 * were introduced to replace the multidimensional float/void arrays
 * previously used in the non-banded functions.
 *
 * EPN, Thu Sep 29 10:01:48 2011 Note: post version 1.0.2, all
 * functions were simplified to take the target sequence length L
 * instead of start and end positions i0 and j0. Now, i0 is implicitly
 * 1 and j0 is implicitly L. To align a subsequence i..j of a larger
 * sequence the caller need only pass dsq+i as dsq and j-i+1 as L.
 * The old method of passing i0 and j0 is leftover from the D&C
 * functions in cm_dpsmall.c upon which many of the functions here
 * were based.
 * 
 * EPN, Thu Sep 29 10:44:19 2011
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

static int   cm_alignT   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, CM_MX    *mx, CM_SHADOW_MX    *shmx, CM_EMIT_MX    *emit_mx, Parsetree_t **ret_tr, float *ret_sc_or_pp);
int          cm_alignT_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, CM_HB_EMIT_MX *emit_mx, Parsetree_t **ret_tr, float *ret_sc_or_pp);


/* Function: cm_alignT()
 * Date:     EPN, Sun Nov 18 19:21:30 2007
 * 
 * Note:     Based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *           Renamed from fast_alignT() [EPN, Wed Sep 14 06:04:39 2011].
 *
 * Purpose:  Call either cm_CYKInsideAlign() (if !<do_optacc>), 
 *           or cm_OptAccAlign()  (if  <do_optacc>),
 *           get vjd shadow matrix; then trace back and
 *           append to an existing but empty parsetree tr.
 *           The full sequence 1..L will be aligned.
 *        
 *           If (<do_optacc>) then emit_mx must != NULL.
 *
 *           Very similar to cm_dpsmall.c:insideT() in case of 
 *           CYK alignment, but uses more efficient implementation
 *           of CYK alignment (cm_CYKInsideAlign()) as opposed to
 *           inside()). 
 *
 * Args:     cm           - the model 
 *           errbuf       - char buffer for reporting errors
 *           dsq          - the digitized sequence [1..L]   
 *           L            - length of the dsq to align
 *           size_limit   - max size in Mb for DP matrix
 *           do_optacc    - TRUE to align with optimal accuracy, else use CYK
 *           mx           - the DP matrix to fill in
 *           shmx         - the shadow matrix to fill in
 *           emit_mx      - the pre-filled emit matrix, must be non-NULL if do_optacc
 *           ret_tr       - RETURN: the optimal parsetree
 *           ret_sc_or_pp - RETURN: optimal score (CYK if !do_optacc, else avg PP of all 1..L residues) 
 * 
 * Returns:  <eslOK>     on success.
 * Throws:   <eslERANGE> if required DP matrix size exceeds <size_limit>, in 
 *                       this case, alignment has been aborted, ret_* variables are not valid
 *           <eslEINVAL> on traceback problem: bogus state
 */
int
cm_alignT(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, 
	  CM_MX *mx, CM_SHADOW_MX *shmx, CM_EMIT_MX *emit_mx, Parsetree_t **ret_tr, float *ret_sc_or_pp)
{
  int       status;
  Parsetree_t *tr = NULL;       /* the parsetree */
  float     sc;			/* the score of the CYK alignment */
  float     pp;			/* avg pp of all emitted residues in optacc alignment */
  ESL_STACK *pda;               /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			/* subseq len for bifurcs */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  int       b;                  /* local begin state */

  if(do_optacc) { if((status = cm_OptAccAlign   (cm, errbuf, dsq, L, size_limit, mx, shmx, emit_mx, &b, &pp)) != eslOK) return status; }
  else          { if((status = cm_CYKInsideAlign(cm, errbuf, dsq, L, size_limit, mx, shmx,          &b, &sc)) != eslOK) return status; };

  /* Create and initialize the parsetree */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  pda = esl_stack_ICreate();
  if(pda == NULL) goto ERROR;
  v = 0;
  i = 1;
  j = d = L;

  while (1) {
    if (cm->sttype[v] == B_st) {
      k = shmx->kshadow[v][j][d];   /* k = len of right fragment */

      /* Store info about the right fragment that we'll retrieve later:
       */
      /* remember the end j */
      if((status = esl_stack_IPush(pda, j))       != eslOK) goto ERROR;	/* remember the end j    */
      if((status = esl_stack_IPush(pda, k))       != eslOK) goto ERROR;	/* remember the subseq length k */
      if((status = esl_stack_IPush(pda, tr->n-1)) != eslOK) goto ERROR;	/* remember the trace index of the parent B state */

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
      yoffset = shmx->yshadow[v][j][d];

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
      default:    ESL_FAIL(eslEINVAL, errbuf, "bogus state type in cm_alignT()");
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

  if(ret_tr       != NULL) *ret_tr = tr; else FreeParsetree(tr);
  if(ret_sc_or_pp != NULL) *ret_sc_or_pp = do_optacc ? pp : sc;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVERREACHED */
}


/* Function: cm_alignT_hb()
 * Date:     EPN 03.29.06
 * 
 * Note:     Based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *           Renamed from fast_alignT_hb() [EPN, Wed Sep 14 06:00:51 2011].
 *
 * Purpose: Call either cm_CYKInsideAlignHB() (if !<do_optacc>), or
 *           cm_OptAccAlignHB() (if <do_optacc>), fill banded vjd
 *           shadow matrix in <shmx>; then trace back.  Append the
 *           trace to a given traceback, which already has state 0 at
 *           tr->n-1. 
 *        
 *           If (<do_optacc>) then emit_mx must != NULL.
 *
 * Args:     cm           - the model 
 *           errbuf       - char buffer for reporting errors
 *           dsq          - the digitized sequence [1..L]   
 *           L            - length of the dsq to align
 *           size_limit   - max size in Mb for DP matrix
 *           do_optacc    - TRUE to align with optimal accuracy, else use CYK
 *           mx           - the DP matrix to fill in
 *           shmx         - the shadow matrix to fill in
 *           emit_mx      - the pre-filled emit matrix, must be non-NULL if do_optacc
 *           ret_tr       - RETURN: the optimal parsetree
 *           ret_sc_or_pp - RETURN: optimal score (CYK if !do_optacc, else avg PP of all 1..L residues) 
 *
 *           
 * Throws:  <eslOK>     on success
 *          <eslERANGE> if required CM_HB_MX exceeds <size_limit>
 *          <eslEINVAL> on traceback problem: bogus state
 */
int
cm_alignT_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc,
	     CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, CM_HB_EMIT_MX *emit_mx, Parsetree_t **ret_tr, float *ret_sc_or_pp)
{
  int       status;
  Parsetree_t *tr = NULL;       /* the parsetree */
  float     sc;			/* the score of the CYK alignment */
  float     pp;			/* avg pp of all emitted residues in optacc alignment */
  ESL_STACK *pda;               /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			/* subseq len for bifurcs */
  /*int       z;*/              /* state index */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  int       b;                  /* local begin state */
  int       jp_v;               /* j-jmin[v] for current j, and current v */
  int       dp_v;               /* d-hd_min(cp9b, v, jp_v) for current j, current v, current d*/
  int       allow_S_local_end;  /* set to true to allow d==0 BEGL_S and BEGR_S local ends if(do_optacc) */

  /* pointers to cp9b data for convenience */
  CP9Bands_t  *cp9b = cm->cp9b;
  int         *jmin = cp9b->jmin;
  int         *jmax = cp9b->jmax;
  int       **hdmin = cp9b->hdmin;
  int       **hdmax = cp9b->hdmax;

  if(do_optacc) { if((status = cm_OptAccAlignHB   (cm, errbuf, dsq, L, size_limit, mx, shmx, emit_mx, &b, &pp)) != eslOK) return status; }
  else          { if((status = cm_CYKInsideAlignHB(cm, errbuf, dsq, L, size_limit, mx, shmx,	      &b, &sc)) != eslOK) return status; }

  /* Create and initialize the parsetree */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  pda = esl_stack_ICreate();
  if(pda == NULL) goto ERROR;
  v = 0;
  i = 1;
  j = d = L;

  while (1) {
    /* special case for HMM banded optimal accuracy, explained below, after the crazy if */
    if(do_optacc && d == 0 && (cm->stid[v] == BEGL_S || cm->stid[v] == BEGR_S) && 
       ((j < jmin[v]             || j > jmax[v]) ||              /* j is outside v's j band */
	(d < hd_min(cp9b, v, j-jmin[v]) || d > hd_max(cp9b, v, j-jmin[v])))) { /* j is within v's j band, but d is outside j's d band */
      /* special case: doing optimal accuracy and v is a BEGL_S or
       * BEGR_S and d is 0 and j is outside v's j band or j is within
       * the band but d is outside j's d band.  We allow this case
       * because although this implies a cell outside the bands, in
       * optimal accuracy only emissions add to the score and we to
       * initialize all cells to IMPOSSIBLE. This means when d==0, we
       * have no way of distinguishing those cells that have been
       * reset to IMPOSSIBLE because they correspond to a valid cell
       * (with valid cells in the B deck, BEGL_S and BEGR_S decks) and
       * those that do not correspond to a valid cell and were never
       * changed since initialization (i.e. this case). So we allow it
       * to prevent an out-of-bounds error. We have to catch it though
       * so we don't try to determine jp_v and dp_v below. We even use
       * USED_EL here if we're not in local mode. You could argue
       * either way whether we should or shouldn't allow this (e.g. we
       * already allow illegal parsetrees in optimal accuracy), but a
       * big reason I decided to allow it is that it is difficult
       * implement a way of disallowing it. Plus the goal of optimal
       * accuracy is to show the alignment that has the maximum
       * average PP on emitted residues within the bands. By allowing
       * this, we also consider a few possible alignments that violate
       * the bands, which I think is okay.
       */
      allow_S_local_end = TRUE; /* this sets yoffset to USED_LOCAL_END in the final 'else' of below code block */
    }
    else if (cm->sttype[v] != EL_st) { /* normal case, determine jp_v, dp_v, j, d offset values given bands */
      jp_v = j - jmin[v];
      dp_v = d - hd_min(cp9b, v, jp_v);
      allow_S_local_end = FALSE;
      assert(j >= jmin[v]        && j <= jmax[v]);
      assert(d >= hd_min(cp9b, v, jp_v) && d <= hd_max(cp9b, v, jp_v));
      ESL_DASSERT1((j >= jmin[v]        && j <= jmax[v]));
      ESL_DASSERT1((d >= hd_min(cp9b, v, jp_v) && d <= hd_max(cp9b, v, jp_v)));
    }

    if (cm->sttype[v] == B_st) {
      k = shmx->kshadow[v][jp_v][dp_v];   /* k = offset len of right fragment */
      /*z = cm->cnum[v];*/
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      if((status = esl_stack_IPush(pda, j)) != eslOK)       goto ERROR;	/* remember the end j    */
      if((status = esl_stack_IPush(pda, k)) != eslOK)       goto ERROR;	/* remember the subseq length k */
      if((status = esl_stack_IPush(pda, tr->n-1)) != eslOK) goto ERROR; /* remember the trace index of the parent B state */
      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    } 
    else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
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
    } 
    else { 
      /* get yoffset */
      if (allow_S_local_end) { 
	yoffset = USED_EL;
      }
      else {
	yoffset = shmx->yshadow[v][jp_v][dp_v];
      }
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    ESL_FAIL(eslEINVAL, errbuf, "Bogus state type in cm_alignT_hb()");
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
      /*ParsetreeDump(stdout, tr, cm, dsq);*/
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */

  /*ParsetreeDump(stdout, tr, cm, dsq);*/

  if(ret_tr       != NULL) *ret_tr = tr; else FreeParsetree(tr);
  if(ret_sc_or_pp != NULL) *ret_sc_or_pp = do_optacc ? pp : sc;
  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "out of memory");
  return status; /* NEVERREACHED */
}

/* Function: cm_AlignSizeNeeded()
 * Date:     EPN, Thu Jan 12 09:51:11 2012
 *
 * Purpose:  Determine size in Mb required to successfully call
 *           cm_Align() for a given model <cm>, sequence length
 *           <L> and alignment options in <do_sample> and <do_post>.
 *
 *           Return <eslERANGE> if required size exceeds size_limit.
 *           
 * Args:     cm         - the covariance model
 *           errbuf     - char buffer for reporting errors
 *           L          - length of sequence 
 *           size_limit - max size in Mb for all required matrices, return eslERANGE if exceeded
 *           do_sample  - TRUE to sample a parsetree from the Inside matrix
 *           do_post    - TRUE to do posteriors
 *           ret_mxmb   - RETURN: size in Mb of required CM_MX (we'll need 2 of these if do_post)
 *           ret_emxmb  - RETURN: size in Mb of required CM_EMIT_MX   (0. if we won't need one) 
 *           ret_shmxmb - RETURN: size in Mb of required CM_SHADOW_MX (0. if we won't need one)
 *           ret_totmb  - RETURN: size in Mb of all required matrices 
 * 
 * Returns: <eslOK> on success.
 * 
 * Throws:  <eslEINVAL> on contract violation
 *          <eslERANGE> if total size of all matrices exceeds <size_limit>
 */
int
cm_AlignSizeNeeded(CM_t *cm, char *errbuf, int L, float size_limit, int do_sample, int do_post,
		   float *ret_mxmb, float *ret_emxmb, float *ret_shmxmb, float *ret_totmb)
{
  int          status;
  float        totmb    = 0.;  /* total Mb required for all matrices (that must be simultaneously in memory) */
  float        mxmb     = 0.;  /* Mb required for CM_MX */
  float        emxmb    = 0.;  /* Mb required for CM_EMIT_MX */
  float        shmxmb   = 0.;  /* Mb required for CM_SHADOW_MX */

  /* we pass NULL values to the *_mx_SizeNeeded() functions because we don't care about cell counts */

  /* we will always need an Inside or CYK matrix */
  if((status = cm_mx_SizeNeeded(cm, errbuf, L, NULL, &mxmb)) != eslOK) return status;
  totmb = mxmb;

  /* if calc'ing posteriors, we'll also need an Outside matrix (which
   * we'll reuse as the Posterior matrix, so only count it once) and
   * an emit matrix.
   */
  if(do_post) { 
    totmb += mxmb; 
    if((status = cm_emit_mx_SizeNeeded(cm, errbuf, L, NULL, NULL, &emxmb)) != eslOK) return status;
    totmb += emxmb;
  }

  /* if we're not sampling an alignment, we'll also need a shadow
   * matrix for the traceback.
   */
  if(! do_sample) { /* if do_sample, we won't need a shadow matrix */
    if((status = cm_shadow_mx_SizeNeeded(cm, errbuf, L, NULL, NULL, &shmxmb)) != eslOK) return status;
    totmb += shmxmb;
  }

  if (ret_mxmb   != NULL) *ret_mxmb    = mxmb;
  if (ret_emxmb  != NULL) *ret_emxmb   = emxmb;
  if (ret_shmxmb != NULL) *ret_shmxmb  = shmxmb;
  if (ret_totmb  != NULL) *ret_totmb   = totmb;

  if(totmb > size_limit) {
    int recommended_mxsize = (int)(ceil(totmb / 1024.0) * 1024.0);
    ESL_FAIL(eslERANGE, errbuf, "non-banded alignment mxes need %.2f Mb > %.2f Mb limit. Use --mxsize %d, --maxtau or --tau.", totmb, (float) size_limit, recommended_mxsize);
  }

  return eslOK;
}

/* Function: cm_AlignSizeNeededHB()
 * Date:     EPN, Thu Jan 12 10:06:20 2012
 *
 * Purpose:  Determine size in Mb required to successfully call
 *           cm_AlignHB() for a given model <cm>, sequence length
 *           <L>, HMM bands <cm->cp9b> and alignment options 
 *           in <do_sample> and <do_post>.
 *
 *           Return <eslERANGE> if required size exceeds size_limit.
 *           
 * Args:     cm          - the covariance model
 *           errbuf      - char buffer for reporting errors
 *           L           - length of sequence 
 *           size_limit  - max size in Mb for all required matrices, return eslERANGE if exceeded
 *           do_sample   - TRUE to sample a parsetree from the Inside matrix
 *           do_post     - TRUE to do posteriors
 *           ret_mxmb    - RETURN: size in Mb of required CM_HB_MX (we'll need 2 of these if do_post)
 *           ret_emxmb   - RETURN: size in Mb of required CM_HB_EMIT_MX   (0. if we won't need one) 
 *           ret_shmxmb  - RETURN: size in Mb of required CM_HB_SHADOW_MX (0. if we won't need one)
 *           ret_cp9mxmb - RETURN: size in Mb of all CP9 matrices needed (fwd and possibly bck) 
 *           ret_cmtotmb - RETURN: size in Mb of all CM matrices (ret_mxmb + ret_emxmb + ret_shmxmb) 
 *           ret_totmb   - RETURN: size in Mb of all required matrices 
 * 
 * Returns: <eslOK> on success.
 * 
 * Throws:  <eslEINVAL> on contract violation
 *          <eslERANGE> if total size of all matrices exceeds <size_limit>
 */
int
cm_AlignSizeNeededHB(CM_t *cm, char *errbuf, int L, float size_limit, int do_sample, int do_post,
		     float *ret_mxmb, float *ret_emxmb, float *ret_shmxmb, float *ret_cp9mxmb, 
                     float *ret_cmtotmb, float *ret_totmb)
{
  int          status;
  float        totmb    = 0.;  /* total Mb required for all matrices including cp9 matrices (that must be simultaneously in memory) */
  float        cmtotmb  = 0.;  /* total Mb required for all CM matrices */
  float        mxmb     = 0.;  /* Mb required for CM_MX */
  float        emxmb    = 0.;  /* Mb required for CM_EMIT_MX */
  float        shmxmb   = 0.;  /* Mb required for CM_SHADOW_MX */
  float        cp9mxmb  = 0.;  /* Mb required for two CP9_MX */

  /* we pass NULL values to the *_mx_SizeNeeded() functions because we don't care about cell counts */

  /* we will always need an Inside or CYK matrix */
  if((status = cm_hb_mx_SizeNeeded(cm, errbuf, cm->cp9b, L, NULL, &mxmb)) != eslOK) return status;
  cmtotmb = mxmb;

  /* if calc'ing posteriors, we'll also need an Outside matrix (which
   * we'll reuse as the Posterior matrix, so only count it once) and
   * an emit matrix.
   */
  if(do_post) { 
    cmtotmb += mxmb; 
    if((status = cm_hb_emit_mx_SizeNeeded(cm, errbuf, cm->cp9b, L, NULL, NULL, &emxmb)) != eslOK) return status;
    cmtotmb += emxmb;
  }

  /* if we're not sampling an alignment, we'll also need a shadow
   * matrix for the traceback.
   */
  if(! do_sample) { 
    if((status = cm_hb_shadow_mx_SizeNeeded(cm, errbuf, cm->cp9b, NULL, NULL, &shmxmb)) != eslOK) return status;
    cmtotmb += shmxmb;
  }

  cp9mxmb = SizeNeededCP9Matrix(L, cm->cp9->M, NULL, NULL);
  cp9mxmb += cp9mxmb; /* add size of bck matrix */

  totmb = cmtotmb + cp9mxmb;

  if (ret_mxmb    != NULL) *ret_mxmb    = mxmb;
  if (ret_emxmb   != NULL) *ret_emxmb   = emxmb;
  if (ret_shmxmb  != NULL) *ret_shmxmb  = shmxmb;
  if (ret_cp9mxmb != NULL) *ret_cp9mxmb = cp9mxmb;
  if (ret_cmtotmb != NULL) *ret_cmtotmb = cmtotmb;
  if (ret_totmb   != NULL) *ret_totmb   = totmb;

#if eslDEBUGLEVEL >= 1  
  printf("#DEBUG: cm_AlignSizeNeededHB()\n");
  printf("#DEBUG: \t mxmb:    %.2f\n", mxmb);
  printf("#DEBUG: \t emxmb:   %.2f\n", emxmb);
  printf("#DEBUG: \t shmxmb:  %.2f\n", shmxmb);
  printf("#DEBUG: \t cp9mxmb: %.2f\n", cp9mxmb);
  printf("#DEBUG: \t cmtotmb: %.2f\n", cmtotmb);
  printf("#DEBUG: \t totmb:   %.2f\n", totmb);
  printf("#DEBUG: \t limit:   %.2f\n", size_limit);
#endif

  if(cmtotmb > size_limit) {
    int recommended_mxsize = (int)(ceil(cmtotmb / 1024.0) * 1024.0);
    ESL_FAIL(eslERANGE, errbuf, "HMM-banded DP mxes need %.1f>%.1f Mb limit (HMM mxes need extra %.1f Mb). Use --mxsize %d.", cmtotmb, (float) size_limit, cp9mxmb, recommended_mxsize);
  }

  return eslOK;
}

/* Function: cm_Align()
 * Date:     EPN, Sun Nov 18 19:26:45 2007
 *
 * Note:     Very similar to cm_dpsmall.c:CYKInside() for case
 *           of CYK alignment, but uses slightly more efficient
 *           implementation (cm_CYKInsideAlign() instead of inside()).
 *           Renamed from FastAlign() [EPN, Wed Sep 14 06:12:46 2011].
 *
 * Purpose: Wrapper for the cm_alignT() routine - solve a full
 *           alignment problem either by CYK, using optimal accuracy,
 *           or sampling, and return the traceback and the score,
 *           without dividing & conquering. Optionally return a
 *           posterior code string.
 *           
 *           Input arguments allow this function to be run in 6 'modes':
 *
 *           mode      returns                 arguments
 *           ----  ---------------  ----------------------------------------
 *                 tr        ppstr  do_optacc  do_sample post_mx   ret_ppstr
 *                 ---------------  ----------------------------------------
 *              1. CYK       no      FALSE      FALSE      NULL      NULL
 *              2. CYK       yes     FALSE      FALSE     !NULL     !NULL
 *              3. Opt acc   no      TRUE       FALSE     !NULL      NULL
 *              4. Opt acc   yes     TRUE       FALSE     !NULL     !NULL
 *              5. sampled   no      FALSE      TRUE       NULL      NULL
 *              6. sampled   yes     FALSE      TRUE      !NULL     !NULL
 *
 *           CYK parsetrees are most the likely parsetree, 'Opt acc'
 *           parsetrees are Holmes/Durbin optimally accurate
 *           parsetrees, the parse that maximizes the summed posterior
 *           probability of emitted residues. A sampled parsetree
 *           is a parsetree sampled from an Inside matrix based on
 *           it's probability.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence, 1..L
 *           L         - length of sequence 
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_optacc - TRUE: do optimal accuracy alignment, not CYK, requires post_mx != NULL
 *           do_sample - TRUE to sample a parsetree from the Inside matrix
 *           mx        - the main dp matrix, grown and filled here, must be non-NULL
 *           shmx      - the shadow matrix, grown and filled here
 *           post_mx   - dp matrix for posterior calculation, grown and filled here, can be NULL only if !do_optacc
 *           emit_mx   - emit matrix to fill
 *           r         - source of randomness, must be non-NULL only if do_sample==TRUE
 *           ret_ppstr - RETURN: posterior code 1, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_avgpp - RETURN: avg PP of emitted residues in parsetree (CYK or optacc) if ret_ppstr == NULL, set as 0.
 *           ret_sc    - RETURN: score of the alignment in bits (Inside score if do_optacc) 
 * 
 * Returns: <eslOK> on success.
 * 
 * Throws:  <eslEINVAL> on contract violation
 *          <eslERANGE> if required CM_MX for Inside/Outside/CYK/Posterior exceeds <size_limit>
 */
int
cm_Align(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, int do_sample,
	 CM_MX *mx, CM_SHADOW_MX *shmx, CM_MX *post_mx, CM_EMIT_MX *emit_mx, ESL_RANDOMNESS *r, 
	 char **ret_ppstr, Parsetree_t **ret_tr, float *ret_avgpp, float *ret_sc)
{
  int          status;
  Parsetree_t *tr = NULL;
  float        sc       = 0.;
  float        avgpp    = 0.;
  float        ins_sc   = 0.;
  int          do_post;
  char        *ppstr = NULL;
  int          have_ppstr;

  have_ppstr = (ret_ppstr != NULL)       ? TRUE : FALSE;
  do_post    = (do_optacc || have_ppstr) ? TRUE : FALSE;

  /* Contract check */
  if(do_optacc && do_sample)         ESL_FAIL(eslEINCOMPAT, errbuf, "cm_Align(), do_optacc and do_sample are both TRUE.");
  if(do_optacc && post_mx == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_Align(), do_optacc is TRUE, but post_mx == NULL.\n");
  if(do_sample && r       == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_Align(), do_sample but r is NULL.");

  /* if do_post:   fill Inside, Outside, Posterior matrices, in that order.
   * if do_sample: fill Inside and sample from it.
   */
  if(do_post || do_sample) { 
    if((status = cm_InsideAlign (cm, errbuf, dsq, L, size_limit, mx,  &ins_sc)) != eslOK) return status;
    if(do_sample) { 
      if((status = cm_StochasticParsetree(cm, errbuf, dsq, L, mx, r, &tr, &sc)) != eslOK) return status; 
    }
    if(do_post) { /* Inside was called above, now do Outside, then Posterior */
      if((status = cm_OutsideAlign(cm, errbuf, dsq, L, size_limit, ((cm->align_opts & CM_ALIGN_CHECKINOUT) && (! (cm->flags & CMH_LOCAL_END))), post_mx, mx, NULL)) != eslOK) return status;
      /* Note: we can only check the posteriors in cm_OutsideAlign() if local begin/ends are off */
      if((status = cm_Posterior       (cm, errbuf, L, size_limit, mx, post_mx, post_mx)) != eslOK) return status;   
      if((status = cm_EmitterPosterior(cm, errbuf, L, size_limit, post_mx, emit_mx, (cm->align_opts & CM_ALIGN_CHECKINOUT))) != eslOK) return status;   
    }
  }

  if(!do_sample) { /* if do_sample, we already have a parsetree */
    if((status = cm_alignT(cm, errbuf, dsq, L, size_limit, do_optacc, mx, shmx, emit_mx, &tr, (do_optacc) ? NULL : &sc)) != eslOK) return status;
  }
  
  if(have_ppstr || do_optacc) { /* call cm_PostCode to get average PP and optionally a PP string (if have_ppstr) */
    if((status = cm_PostCode(cm, errbuf, L, emit_mx, tr, (have_ppstr) ? &ppstr : NULL, &avgpp)) != eslOK) return status;
  }

  if (ret_ppstr  != NULL) *ret_ppstr  = ppstr; else free(ppstr);
  if (ret_tr     != NULL) *ret_tr     = tr;    else FreeParsetree(tr);
  if (ret_avgpp  != NULL) *ret_avgpp  = avgpp;
  if (ret_sc     != NULL) *ret_sc     = (do_optacc) ? ins_sc : sc;

  ESL_DPRINTF1(("#DEBUG: returning from cm_Align() sc : %f\n", sc)); 
  return eslOK;
}



/* Function: cm_AlignHB()
 * Incept:   EPN, Fri Oct 26 09:31:43 2007
 * 
 * Note:     Based on CYKInside_b_jd() [11.04.05] which was based on CYKInside_b() 
 *           which was based on CYKInside() [SRE, Sun Jun  3 19:48:33 2001 [St. Louis]]
 *           Renamed from cm_AlignHB() [EPN, Wed Sep 14 06:09:51 2011].
 *
 * Purpose: Wrapper for the cm_alignT() routine - solve a full
 *           alignment problem either by CYK, using optimal accuracy,
 *           or sampling, and return the traceback and the score,
 *           without dividing & conquering. Optionally return a
 *           posterior code string.
 *           
 *           Identical to cm_Align() but HMM bands are used here.
 *           See that function's 'Purpose' for more details.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence, 1..L
 *           L         - length of sequence 
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_optacc - TRUE: do optimal accuracy alignment, not CYK, requires post_mx != NULL
 *           do_sample - TRUE: sample a parsetree from the Inside matrix
 *           mx        - the main dp matrix, grown and filled here, must be non-NULL
 *           shmx      - the shadow matrix, grown and filled here
 *           post_mx   - dp matrix for posterior calculation, grown and filled here, can be NULL only if !do_optacc
 *           emit_mx   - emit matrix to fill
 *           r         - source of randomness, must be non-NULL only if do_sample==TRUE
 *           ret_ppstr - RETURN: posterior code 1, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_avgpp - RETURN: avg PP of emitted residues in parsetree (CYK or optacc) if ret_ppstr == NULL, set as 0.
 *           ret_sc    - RETURN: score of the alignment in bits (Inside score if do_optacc) 
 * 
 * Returns: <eslOK> on success
 * 
 * Throws:  <eslEINVAL> on contract violation
 *          <eslERANGE> if required CM_HB_MX for Inside/Outside/CYK/Posterior exceeds <size_limit>
 */

int
cm_AlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, int do_sample, 
	   CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, CM_HB_MX *post_mx, CM_HB_EMIT_MX *emit_mx, ESL_RANDOMNESS *r, 
	   char **ret_ppstr, Parsetree_t **ret_tr, float *ret_avgpp, float *ret_sc)
{
  int          status;
  Parsetree_t *tr = NULL;
  float        sc     = 0.;
  float        avgpp  = 0.;
  float        ins_sc = 0.;
  int          do_post;
  char        *ppstr = NULL;
  int          have_ppstr;

  have_ppstr = (ret_ppstr != NULL)       ? TRUE : FALSE;
  do_post    = (do_optacc || have_ppstr) ? TRUE : FALSE;

  /* Contract check */
  if(do_optacc && do_sample)         ESL_FAIL(eslEINCOMPAT, errbuf, "cm_AlignHB(), do_optacc and do_sample are both TRUE.");
  if(do_optacc && post_mx == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_AlignHB(), do_optacc is TRUE, but post_mx == NULL.\n");
  if(do_sample && r       == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_AlignHB(), do_sample but r is NULL.");

  /* PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, L); */

  /* Sub-stage timing instrumentation.
   * Outputs "#CM_ALIGN_HB_SUBSTAGE M=<M> L=<L> stage=<name> s=<seconds>" to stderr.
   */
  struct timespec _ta_sub, _tb_sub;
#define CM_ALIGN_SUBSTAGE_START()  clock_gettime(CLOCK_MONOTONIC, &_ta_sub)
#define CM_ALIGN_SUBSTAGE_END(tag) do { \
    clock_gettime(CLOCK_MONOTONIC, &_tb_sub); \
    double _ss = (_tb_sub.tv_sec - _ta_sub.tv_sec) + (_tb_sub.tv_nsec - _ta_sub.tv_nsec)/1e9; \
    fprintf(stderr, "#CM_ALIGN_HB_SUBSTAGE M=%d L=%d stage=%s s=%.6f\n", \
            (cm->fp7 ? cm->fp7->M : 0), L, (tag), _ss); \
    fflush(stderr); \
  } while(0)

  /* if do_post:   fill Inside, Outside, Posterior matrices, in that order.
   * if do_sample: fill Inside and sample from it.
   */
  if(do_post || do_sample) {
    CM_ALIGN_SUBSTAGE_START();
    if((status = cm_InsideAlignHB (cm, errbuf, dsq, L, size_limit, mx, &ins_sc)) != eslOK) return status;
    CM_ALIGN_SUBSTAGE_END("inside");
    if(do_sample) {
      CM_ALIGN_SUBSTAGE_START();
      if((status = cm_StochasticParsetreeHB(cm, errbuf, dsq, L, mx, r, &tr, &sc)) != eslOK) return status;
      CM_ALIGN_SUBSTAGE_END("sample_trace");
    }
    if(do_post) { /* Inside was called above, now do Outside, then Posterior */
      CM_ALIGN_SUBSTAGE_START();
      if((status = cm_OutsideAlignHB(cm, errbuf, dsq, L, size_limit, ((cm->align_opts & CM_ALIGN_CHECKINOUT) && (! (cm->flags & CMH_LOCAL_END))), post_mx, mx, NULL)) != eslOK) return status;
      CM_ALIGN_SUBSTAGE_END("outside");
      /* Note: we can only check the posteriors in cm_OutsideAlignHB() if local begin/ends are off */
      CM_ALIGN_SUBSTAGE_START();
      if((status = cm_PosteriorHB       (cm, errbuf, L, size_limit, mx, post_mx, post_mx)) != eslOK) return status;
      if((status = cm_EmitterPosteriorHB(cm, errbuf, L, size_limit, post_mx, emit_mx, (cm->align_opts & CM_ALIGN_CHECKINOUT))) != eslOK) return status;
      CM_ALIGN_SUBSTAGE_END("posterior");
    }
  }

  if(!do_sample) { /* if do_sample, we already have a parsetree */
    CM_ALIGN_SUBSTAGE_START();
    if((status = cm_alignT_hb(cm, errbuf, dsq, L, size_limit, do_optacc, mx, shmx, emit_mx, &tr, (do_optacc) ? NULL : &sc)) != eslOK) return status;
    CM_ALIGN_SUBSTAGE_END("traceback");
  }

  if(have_ppstr || do_optacc) {
    CM_ALIGN_SUBSTAGE_START();
    if((status = cm_PostCodeHB(cm, errbuf, L, emit_mx, tr, (have_ppstr) ? &ppstr : NULL, &avgpp)) != eslOK) return status;
    CM_ALIGN_SUBSTAGE_END("postcode");
  }

#undef CM_ALIGN_SUBSTAGE_START
#undef CM_ALIGN_SUBSTAGE_END

  /* Uncomment to dump emit map and parse tree */
  /* CMEmitMap_t *emap;
     emap = CreateEmitMap(cm);
     DumpEmitMap(stdout, emap, cm);
     FreeEmitMap(emap);
     ParsetreeDump(stdout, tr, cm, dsq);
  */

  if (ret_ppstr  != NULL) *ret_ppstr  = ppstr; else free(ppstr);
  if (ret_tr     != NULL) *ret_tr     = tr;    else FreeParsetree(tr);
  if (ret_avgpp  != NULL) *ret_avgpp  = avgpp;
  if (ret_sc     != NULL) *ret_sc     = (do_optacc) ? ins_sc : sc;

  ESL_DPRINTF1(("#DEBUG: returning from cm_AlignHB() sc : %f\n", sc)); 
  return eslOK;
}

/*****************************************************************
 * Checkpointed (sqrt(M)-memory) HMM-banded posterior + OptAcc
 * alignment, for NON-truncated, GLOBAL, bps=0 (pure left-emitting
 * MATL chain) CMs.
 *
 * Brief 26_0610-028: library port of the validated standalone drivers
 *   ckpt_drv.c   (022) -- checkpointed Inside + Outside + fused posterior
 *   ckptoa_drv.c (023) -- checkpointed OptAcc max-DP + checkpointed traceback
 *
 * The recurrences mirror cm_InsideAlignHB / cm_OutsideAlignHB /
 * cm_EmitterPosteriorHB / cm_OptAccAlignHB / cm_alignT_hb cell-for-cell
 * (FLogsum/max call order preserved) for the state types that occur in a
 * pure MATL chain (S, IL, IR, ML, D, E) in global mode (no local
 * begins/ends, no EL).  Only deck *storage* is windowed (sqrt(M) seed decks
 * + recompute), so the output (Z, l_pp/r_pp, parsetree, per-residue PP) is
 * byte-for-byte identical to the stock path.
 *
 * Routed by the CM_ALIGN_CHECKPT align_opt; cm_CheckptAlignHB_Qualifies()
 * gates engagement (stock path used otherwise).  These functions do NOT
 * touch cm->hb_mx / hb_omx / hb_shmx (they allocate their own sqrt(M)
 * decks); they DO grow + fill the passed emit_mx, which is the deliverable
 * read by cm_PostCodeHB().
 *****************************************************************/

/* per-call context (replaces the drivers' file-scope statics; keeps the
 * library reentrant / thread-safe -- cmalign may run multithreaded). */
typedef struct ckpt_ctx_s {
  CM_t    *cm;
  ESL_DSQ *dsq;
  int      L;
  int      M;
  int     *jmin, *jmax, *imin, *imax;
  int    **hdmin, **hdmax;
  int64_t *deck_nc;    /* [v] number of cells in deck v */
  int     *deck_njr;   /* [v] number of j-rows in deck v (jmax-jmin+1) */
  float  **my_lpp;     /* alias to emit_mx->l_pp (left-emitter posteriors)  */
  float  **my_rpp;     /* alias to emit_mx->r_pp (right-emitter posteriors) */
  int64_t  cur_bytes;  /* live CM-DP cell bytes (working-set tracking)       */
  int64_t  peak_bytes; /* high-water mark of cur_bytes                       */
  int     *kpin;       /* rung-3: per-B pinned right-frag length k* (or NULL)*/
  float ***ifull;      /* rung-3: full Inside alpha, for Outside sibling reads at B (or NULL) */
  /* R-L.2: EL (local-end) support, all NULL/0 unless CMH_LOCAL_END.  The two
   * data-dependent EL decks (Outside beta[cm->M], OptAcc alpha[cm->M]) are
   * banded on the UPPER d-edge: row j stores d=0..eldmax[j] (cells above are
   * provably IMPOSSIBLE).  The Inside/forward EL needs no deck (fixed ramp). */
  int      have_el;    /* TRUE iff cm->flags & CMH_LOCAL_END                  */
  float    el_selfsc;  /* cm->el_selfsc (EL self-loop emit score)             */
  int     *eldmax;     /* [0..L] upper d-band edge for EL rows; -1 = empty    */
  float  **elbeta;     /* [0..L][0..eldmax[j]] Outside EL deck (Step B)       */
  float  **elalpha;    /* [0..L][0..eldmax[j]] OptAcc EL prefix-sum (Step OA) */
  float   *el_esc;     /* [0..M-1] all-deletes-to-E score (OA d==0 EL route)  */
  float    el_endsc;   /* the (shared) v->EL transition score                 */
  /* R-L.2b: local BEGIN support, all 0/-1/IMPOSSIBLE unless CMH_LOCAL_BEGIN.
   * Begins have no per-cell deck (unlike EL): they're a ROOT_S reduction over
   * begin-state root-diagonal cells (v's cell (L,L)), so this is just a few
   * scalars threaded through the per-deck sweeps.  Forward (bsc_fwd) sums via
   * FLogsum over all v visited before ROOT_S (non-idempotent -> needs the
   * fwd_begin_applied guard so a later re-materialization of deck 0 doesn't
   * re-apply it).  OA (begin_bsc/begin_b) is a running max -> idempotent, no
   * guard needed; the same scalar persists correctly across Step OA and the
   * later Step TB block-recompute of deck 0. */
  int      have_local_begin;  /* TRUE iff cm->flags & CMH_LOCAL_BEGIN          */
  float    bsc_fwd;           /* Step A: running FLogsum(alpha[v][L][L]+beginsc[v]) */
  int      fwd_begin_applied; /* Step A: TRUE once applied to alpha[0][L][L]   */
  float    begin_bsc;         /* Step OA/TB: running max alpha[v][L][L] over begin states */
  int      begin_b;           /* Step OA/TB: argmax begin state (-1 if none)   */
} CKPT_CTX;

/* deck = float** of j-rows; row[0] is start of the contiguous cell block. */
static float **
ckpt_deck_alloc(CKPT_CTX *cx, int v)
{
  int     njr = cx->deck_njr[v];
  int64_t nc  = cx->deck_nc[v];
  float **row = malloc(sizeof(float*) * (njr > 0 ? njr : 1));
  float  *mem = malloc(sizeof(float)  * (nc  > 0 ? nc  : 1));
  int64_t off = 0;
  int     jp;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  if (row == NULL || mem == NULL) cm_Fail("ckpt_deck_alloc OOM v=%d", v);
  row[0] = mem; /* ensure row[0]==mem even when njr==0 (so free(row[0]) is valid) */
  for (jp = 0; jp < njr; jp++) {
    int w = hd_max(cp9b, v, jp) - hd_min(cp9b, v, jp) + 1;
    row[jp] = mem + off;
    if (w < 0) w = 0;
    off += w;
  }
  cx->cur_bytes += nc * (int64_t)sizeof(float);
  if (cx->cur_bytes > cx->peak_bytes) cx->peak_bytes = cx->cur_bytes;
  return row;
}
static void
ckpt_deck_free(CKPT_CTX *cx, int v, float **row)
{
  if (row == NULL) return;
  free(row[0]);
  free(row);
  cx->cur_bytes -= cx->deck_nc[v] * (int64_t)sizeof(float);
}
static void
ckpt_deck_init_impossible(CKPT_CTX *cx, int v, float **row)
{
  if (cx->deck_nc[v] > 0) esl_vec_FSet(row[0], (int) cx->deck_nc[v], IMPOSSIBLE);
}
/* char deck (yshadow), same row layout as a float deck */
static char **
ckpt_cdeck_alloc(CKPT_CTX *cx, int v)
{
  int     njr = cx->deck_njr[v];
  int64_t nc  = cx->deck_nc[v];
  char  **row = malloc(sizeof(char*) * (njr > 0 ? njr : 1));
  char   *mem = malloc(sizeof(char)  * (nc  > 0 ? nc  : 1));
  int64_t off = 0;
  int     jp;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  if (row == NULL || mem == NULL) cm_Fail("ckpt_cdeck_alloc OOM v=%d", v);
  row[0] = mem;
  for (jp = 0; jp < njr; jp++) {
    int w = hd_max(cp9b, v, jp) - hd_min(cp9b, v, jp) + 1;
    row[jp] = mem + off;
    if (w < 0) w = 0;
    off += w;
  }
  cx->cur_bytes += nc * (int64_t)sizeof(char);
  if (cx->cur_bytes > cx->peak_bytes) cx->peak_bytes = cx->cur_bytes;
  return row;
}
static void
ckpt_cdeck_free(CKPT_CTX *cx, int v, char **row)
{
  if (row == NULL) return;
  free(row[0]);
  free(row);
  cx->cur_bytes -= cx->deck_nc[v] * (int64_t)sizeof(char);
}

/* ---- R-L.2 EL (local-end) deck helpers --------------------------------- *
 * The two data-dependent EL decks live at the (non-banded) EL state cm->M,
 * indexed by ABSOLUTE (j,d) (j=0..L).  In stock these are full O(L^2/2)
 * lower-triangular decks; here we band the UPPER d-edge to eldmax[j], the
 * largest d any local-end-capable state references at EL-row j.  Cells with
 * d > eldmax[j] are provably IMPOSSIBLE in stock (no v->EL write, and the
 * self-transition/prefix-sum only descend from a written cell), so dropping
 * them is byte-exact (FLogsum(x,IMPOSSIBLE)==x exactly since IMPOSSIBLE=-1e36).
 */

/* Fill cx->eldmax[0..L]: for each local-end state v, EL-row r = (v's j) - sdr
 * references d-indices up to hdmax[v][jp] - sd.  This single bound covers BOTH
 * the Outside v->EL write range and the OptAcc alpha[cm->M] read range (same
 * (r,d) geometry).  -1 marks an EL-row no local end reaches. */
static void
ckpt_el_compute_dmax(CKPT_CTX *cx)
{
  CM_t *cm = cx->cm;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  int   L  = cx->L, M = cx->M;
  int   v, jp, r;
  for (r = 0; r <= L; r++) cx->eldmax[r] = -1;
  for (v = 0; v < M; v++) {
    if (! NOT_IMPOSSIBLE(cm->endsc[v])) continue;
    int sd  = StateDelta(cm->sttype[v]);
    int sdr = StateRightDelta(cm->sttype[v]);
    int njr = cx->deck_njr[v];
    for (jp = 0; jp < njr; jp++) {
      int j_band = cx->jmin[v] + jp;   /* v's own j */
      r = j_band - sdr;                /* EL-row */
      if (r < 0 || r > L) continue;
      int dmax_here = hd_max(cp9b, v, jp) - sd;
      if (dmax_here > r)        dmax_here = r;   /* d <= j on the EL diagonal */
      if (dmax_here > cx->eldmax[r]) cx->eldmax[r] = dmax_here;
    }
  }
}

/* Allocate a banded EL deck [0..L][0..eldmax[j]] (NULL rows where eldmax<0),
 * initialized to IMPOSSIBLE; account its cells in the working-set tracker. */
static float **
ckpt_el_deck_alloc(CKPT_CTX *cx)
{
  int     L = cx->L, r, d;
  float **deck = malloc(sizeof(float *) * (L+1));
  if (deck == NULL) cm_Fail("ckpt_el_deck_alloc OOM");
  for (r = 0; r <= L; r++) {
    if (cx->eldmax[r] < 0) { deck[r] = NULL; continue; }
    int n = cx->eldmax[r] + 1;
    deck[r] = malloc(sizeof(float) * n);
    if (deck[r] == NULL) cm_Fail("ckpt_el_deck_alloc OOM row %d", r);
    for (d = 0; d < n; d++) deck[r][d] = IMPOSSIBLE;
    cx->cur_bytes += (int64_t) n * (int64_t)sizeof(float);
  }
  if (cx->cur_bytes > cx->peak_bytes) cx->peak_bytes = cx->cur_bytes;
  return deck;
}
static void
ckpt_el_deck_free(CKPT_CTX *cx, float **deck)
{
  int L = cx->L, r;
  if (deck == NULL) return;
  for (r = 0; r <= L; r++) {
    if (deck[r] != NULL) { cx->cur_bytes -= (int64_t)(cx->eldmax[r]+1) * (int64_t)sizeof(float); free(deck[r]); }
  }
  free(deck);
}

/* Precompute cx->el_esc[0..M-1] (the all-deletes-to-E path score) + cx->el_endsc,
 * mirroring cm_InitializeOptAccShadowDZeroHB:7664-7724 (bps=0: no B states).
 * Used by ckpt_optacc_deck's d==0 EL routing.  Caller allocs cx->el_esc. */
static void
ckpt_el_compute_esc(CKPT_CTX *cx)
{
  CM_t *cm = cx->cm;
  int   M = cx->M, v;
  esl_vec_FSet(cx->el_esc, M, IMPOSSIBLE);
  v = 0; while (v < M && ! NOT_IMPOSSIBLE(cm->endsc[v])) v++;
  cx->el_endsc = (v < M) ? cm->endsc[v] : IMPOSSIBLE;
  for (v = M-1; v >= 0; v--) {
    if (! cm->cp9b->Jvalid[v]) continue;
    if (cm->sttype[v] == E_st) { cx->el_esc[v] = 0.; }
    else {
      int y = cm->cfirst[v];
      while (StateDelta(cm->sttype[y]) != 0) y++;
      cx->el_esc[v] = cx->el_esc[y] + cm->tsc[v][y - cm->cfirst[v]];
    }
  }
}

/* R-L.2b local-begin, Forward/Inside half: mirrors cm_InsideAlignHB's
 * "for (v...) { ...; allow local begins ... bsc = FLogsum(bsc, alpha[v][L][L]
 * + beginsc[v]); } ; alpha[0][jp_0][Lp_0] = FLogsum(alpha[0][jp_0][Lp_0], bsc);"
 * Called once per deck v, AFTER that deck's cells are fully computed (right
 * before ckpt_inside_deck's return).  cx->bsc_fwd accumulates across all v
 * visited before ROOT_S in a single Step-A-style descending sweep (v=M-1..0);
 * frozen (via fwd_begin_applied) the first time it's folded into alpha[0], so
 * a later re-materialization of deck 0 (Step B's per-block ba[] recompute,
 * whose value is otherwise unused) doesn't re-add it and corrupt the sum. */
static void
ckpt_apply_begin_fwd(CKPT_CTX *cx, int v, float **av)
{
  CM_t *cm = cx->cm;
  int   L  = cx->L;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  if (! cx->have_local_begin) return;
  if (v == 0) {
    if (! cx->fwd_begin_applied) {
      int jp_0 = L - cx->jmin[0];
      int Lp_0 = L - hd_min(cp9b, 0, jp_0);
      if (NOT_IMPOSSIBLE(cx->bsc_fwd)) av[jp_0][Lp_0] = FLogsum(av[jp_0][Lp_0], cx->bsc_fwd);
      cx->fwd_begin_applied = TRUE;
    }
    return;
  }
  if (cx->fwd_begin_applied) return; /* Step B/TB re-materialization: already folded, skip */
  if (NOT_IMPOSSIBLE(cm->beginsc[v]) && L >= cx->jmin[v] && L <= cx->jmax[v]) {
    int jp_v = L - cx->jmin[v];
    if (L >= hd_min(cp9b, v, jp_v) && L <= hd_max(cp9b, v, jp_v)) {
      int Lp = L - hd_min(cp9b, v, jp_v);
      cx->bsc_fwd = FLogsum(cx->bsc_fwd, av[jp_v][Lp] + cm->beginsc[v]);
    }
  }
}

/* R-L.2b local-begin, OptAcc half: mirrors cm_OptAccAlignHB's "allow local
 * begins ... if (alpha[v][jp_v][Lp] > bsc) { b=v; bsc=alpha[v][jp_v][Lp]; }"
 * followed by (after the full v loop) "alpha[0][jp_0][Lp_0]=bsc;
 * yshadow[0][jp_0][Lp_0]=USED_LOCAL_BEGIN;" (unconditional overwrite -- OA's
 * recursion never adds tsc, so ROOT_S's own value doesn't need comparing).
 * A running MAX is idempotent under re-comparison, so unlike the Forward
 * half this needs no freeze: Step TB's block recompute of deck 0 re-derives
 * (or re-confirms) the same (begin_b, begin_bsc) and reapplies the identical
 * override to tysh[0], which is exactly what the traceback needs to see. */
static void
ckpt_apply_begin_oa(CKPT_CTX *cx, int v, float **av, char **ysh)
{
  CM_t *cm = cx->cm;
  int   L  = cx->L;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  if (! cx->have_local_begin) return;
  if (NOT_IMPOSSIBLE(cm->beginsc[v]) && L >= cx->jmin[v] && L <= cx->jmax[v]) {
    int jp_v = L - cx->jmin[v];
    if (L >= hd_min(cp9b, v, jp_v) && L <= hd_max(cp9b, v, jp_v)) {
      int Lp = L - hd_min(cp9b, v, jp_v);
      if (av[jp_v][Lp] > cx->begin_bsc) { cx->begin_bsc = av[jp_v][Lp]; cx->begin_b = v; }
    }
  }
  if (v == 0 && NOT_IMPOSSIBLE(cx->begin_bsc)) {
    int jp_0 = L - cx->jmin[0];
    int Lp_0 = L - hd_min(cp9b, 0, jp_0);
    av[jp_0][Lp_0] = cx->begin_bsc;
    if (ysh != NULL) ysh[jp_0][Lp_0] = (char) USED_LOCAL_BEGIN;
  }
}

/* Brief 26_0610-084 (R1-L) local-begin, CYK half: mirrors ckpt_apply_begin_oa
 * but for the CYK MAX recurrence, transcribing cm_CYKInsideAlignHB's own
 * local-begin block (cm_dpalign.c:4241-4272).  Two differences vs the OA half:
 *   (a) the running max tracks alpha[v][L][L] + cm->beginsc[v] (CYK adds the
 *       0->b begin transition score; OA is posterior-mass and adds nothing), and
 *   (b) the ROOT_S override is CONDITIONAL (only if the best begin beats the
 *       root's own non-begin score), exactly as cm_CYKInsideAlignHB's
 *       "if (bsc > alpha[0][jp_0][Lp_0])" -- OA overrides unconditionally
 *       because its ROOT_S recursion never adds a transition score.
 * Like the OA half, a running MAX is idempotent under re-comparison, so Step TB's
 * block recompute of deck 0 re-derives the same (begin_b, begin_bsc) from the
 * still-persisted cx fields and reapplies the identical override to tysh[0]. */
static void
ckpt_apply_begin_cyk(CKPT_CTX *cx, int v, float **av, char **ysh)
{
  CM_t *cm = cx->cm;
  int   L  = cx->L;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  if (! cx->have_local_begin) return;
  if (NOT_IMPOSSIBLE(cm->beginsc[v]) && L >= cx->jmin[v] && L <= cx->jmax[v]) {
    int jp_v = L - cx->jmin[v];
    if (L >= hd_min(cp9b, v, jp_v) && L <= hd_max(cp9b, v, jp_v)) {
      int Lp = L - hd_min(cp9b, v, jp_v);
      float cand = av[jp_v][Lp] + cm->beginsc[v];
      if (cand > cx->begin_bsc) { cx->begin_bsc = cand; cx->begin_b = v; }
    }
  }
  if (v == 0 && NOT_IMPOSSIBLE(cx->begin_bsc)) {
    int jp_0 = L - cx->jmin[0];
    int Lp_0 = L - hd_min(cp9b, 0, jp_0);
    if (cx->begin_bsc > av[jp_0][Lp_0]) {   /* CONDITIONAL override (CYK), unlike OA's unconditional */
      av[jp_0][Lp_0] = cx->begin_bsc;
      if (ysh != NULL) ysh[jp_0][Lp_0] = (char) USED_LOCAL_BEGIN;
    }
  }
}

/* Inside deck v: mirrors cm_InsideAlignHB for S/IL/IR/ML/D/E, global mode.
 * Reads children from ba[] (in-block decks) or ck[] (checkpoint seeds). */
static void
ckpt_inside_deck(CKPT_CTX *cx, int v, float ***ba, float ***ck)
{
#define IA(vv) (ba[vv] ? ba[vv] : (ck ? ck[vv] : NULL))
  CM_t *cm = cx->cm;
  ESL_DSQ *dsq = cx->dsq;
  int   *jmin = cx->jmin, *jmax = cx->jmax;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  float **av = ba[v];
  float const *esc_v = cm->oesc[v];
  float const *tsc_v = cm->tsc[v];
  int sd  = StateDelta(cm->sttype[v]);
  int sdr = StateRightDelta(cm->sttype[v]);
  int j, d, i, y, yoffset, jp_v, dp_v, jp_y_sdr, dp_y_sd, j_sdr;
  int yvalidA[MAXCONNECT], yvalid_ct, yvalid_idx;

  ckpt_deck_init_impossible(cx, v, av);

  /* R-L.2 EL: re-init this state's deck if a local end from v is allowed.  The
   * forward EL is the FIXED ramp el_scA[d-sd]=el_selfsc*(d-sd) -- NO alpha[cm->M]
   * deck is read (mirrors cm_InsideAlignHB:3576-3583, which fills the EL deck only
   * "for completeness").  Placed before the per-state recurrence so children
   * FLogsum onto, and emissions add to, this EL base -- exactly as stock. */
  if (cx->have_el && NOT_IMPOSSIBLE(cm->endsc[v])) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      for (dp_v = 0, d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); dp_v++, d++)
        av[jp_v][dp_v] = cx->el_selfsc * (d - sd) + cm->endsc[v];
    }
  }

  if (cm->sttype[v] == E_st) {
    for (j = jmin[v]; j <= jmax[v]; j++) { jp_v = j - jmin[v]; av[jp_v][0] = 0.; }
    ckpt_apply_begin_fwd(cx, v, av);
    return;
  }
  else if (cm->sttype[v] == IL_st) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      yvalid_ct = 0; j_sdr = j - sdr;
      for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++)
        if ((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset;
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
        i = j - d + 1;
        dp_v = d - hd_min(cp9b, v, jp_v);
        for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) {
          yoffset = yvalidA[yvalid_idx];
          y = cm->cfirst[v] + yoffset;
          jp_y_sdr = j - jmin[y] - sdr;
          if ((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) {
            dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
            av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], IA(y)[jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
          }
        }
        av[jp_v][dp_v] += esc_v[dsq[i--]];
        av[jp_v][dp_v] = ESL_MAX(av[jp_v][dp_v], IMPOSSIBLE);
      }
    }
    ckpt_apply_begin_fwd(cx, v, av);
    return;
  }
  else if (cm->sttype[v] == IR_st) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      yvalid_ct = 0; j_sdr = j - sdr;
      for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++)
        if ((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset;
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
        dp_v = d - hd_min(cp9b, v, jp_v);
        for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) {
          yoffset = yvalidA[yvalid_idx];
          y = cm->cfirst[v] + yoffset;
          jp_y_sdr = j - jmin[y] - sdr;
          if ((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) {
            dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
            av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], IA(y)[jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
          }
        }
        av[jp_v][dp_v] += esc_v[dsq[j]];
        av[jp_v][dp_v] = ESL_MAX(av[jp_v][dp_v], IMPOSSIBLE);
      }
    }
    ckpt_apply_begin_fwd(cx, v, av);
    return;
  }
  else if (cm->sttype[v] == B_st) { /* rung-3: pinned bifurcation (single k*) */
    int z = cm->cnum[v];            /* right (BEGR) subtree; y=cfirst[v]=left (BEGL) */
    int kpinned = (cx->kpin != NULL) ? cx->kpin[v] : -1;
    int jp_z, jp_y, kp_z, k, kn, kx, dp_y;
    y = cm->cfirst[v];
    int jnn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
    int jxx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
    for (j = jnn; j <= jxx; j++) {
      jp_v = j - jmin[v];
      jp_y = j - jmin[y];
      jp_z = j - jmin[z];
      kn = ((j-jmax[y]) > (hd_min(cp9b, z, jp_z))) ? (j-jmax[y]) : hd_min(cp9b, z, jp_z);
      kn = ESL_MAX(kn, 0);
      kx = ( jp_y       < (hd_max(cp9b, z, jp_z))) ?  jp_y       : hd_max(cp9b, z, jp_z);
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
        dp_v = d - hd_min(cp9b, v, jp_v);
        /* R-L.4b (065): if CYK's pin walk never visited this B state (kpinned<0)
         * AND more than one k in [kn,kx] actually satisfies the per-cell band
         * constraint below for THIS (j,d), do NOT fall back to an unrestricted
         * full split-search -- that would compute a real value with no
         * shadow/traceback record of which k won, and OA's own traceback
         * could then have no pin to reconstruct from if it ever tried to step
         * into this cell.  Leave the cell at its IMPOSSIBLE-initialized
         * default instead; CYK's own reconstruction (built only from pinned B
         * states) is always a valid, fully-computable alternative, so no
         * legal traceback is ever forced through here.
         * EXCEPTION: when exactly one k satisfies the constraint (checked
         * per-d, NOT just via the j-only [kn,kx] window -- kn==kx is only a
         * necessary, not sufficient, proxy: e.g. every d==0 cell forces a
         * unique k=0 even when [kn,kx] spans many values), the split is
         * unambiguous and needs no pin -- always compute it (this is NOT an
         * untracked choice; ckpt_optacc_traceback re-derives the same unique
         * k the same way, see there). */
        int klo, khi;
        if (kpinned >= 0) {
          if (kpinned < kn || kpinned > kx) continue;
          klo = khi = kpinned;
        } else {
          int uniq_k = -1, n_valid = 0;
          for (k = kn; k <= kx; k++) {
            if ((k >= d - hd_max(cp9b, y, jp_y-k)) && k <= d - hd_min(cp9b, y, jp_y-k)) {
              n_valid++; uniq_k = k;
              if (n_valid > 1) break;
            }
          }
          if (n_valid != 1) continue;
          klo = khi = uniq_k;
        }
        for (k = klo; k <= khi; k++) {
          if ((k >= d - hd_max(cp9b, y, jp_y-k)) && k <= d - hd_min(cp9b, y, jp_y-k)) {
            kp_z = k-hd_min(cp9b, z, jp_z);
            dp_y = d-hd_min(cp9b, y, jp_y-k);
            av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], IA(y)[jp_y-k][dp_y - k] + IA(z)[jp_z][kp_z]);
          }
        }
      }
    }
    ckpt_apply_begin_fwd(cx, v, av);
    return;
  }
  else { /* ML, MP, MR, D, S (no self-transit, no B) */
    int jn, jx, jpn, jpx, dn, dx, dpn, dpx;
    float tsc;
    for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
      yoffset = y - cm->cfirst[v];
      tsc = tsc_v[yoffset];
      jn = ESL_MAX(jmin[v], jmin[y] + sdr);
      jx = ESL_MIN(jmax[v], jmax[y] + sdr);
      jpn = jn - jmin[v];
      jpx = jx - jmin[v];
      jp_y_sdr = jn - jmin[y] - sdr;
      for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
        dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y_sdr) + sd);
        dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y_sdr) + sd);
        dpn = dn - hd_min(cp9b, v, jp_v);
        dpx = dx - hd_min(cp9b, v, jp_v);
        dp_y_sd = dn - hd_min(cp9b, y, jp_y_sdr) - sd;
        for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) {
          av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], (IA(y)[jp_y_sdr][dp_y_sd] + tsc));
        }
      }
    }
    switch (cm->sttype[v]) {        /* rung-3: ML (left), MR (right), MP (pair) emission */
    case ML_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        i = j - hd_min(cp9b, v, jp_v) + 1;
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
          av[jp_v][dp_v] += esc_v[dsq[i--]];
      }
      break;
    case MR_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
          av[jp_v][dp_v] += esc_v[dsq[j]];
      }
      break;
    case MP_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        i = j - hd_min(cp9b, v, jp_v) + 1;
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
          av[jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
      }
      break;
    default: break;
    }
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
        av[jp_v][dp_v] = ESL_MAX(av[jp_v][dp_v], IMPOSSIBLE);
    }
    (void) jn; (void) jx; (void) jpn; (void) jpx; (void) dn; (void) dx; (void) dpn; (void) dpx;
    ckpt_apply_begin_fwd(cx, v, av);
    return;
  }
#undef IA
}

/* Outside deck v: mirrors cm_OutsideAlignHB for S/IL/IR/ML/D/E, global mode.
 * Reads beta from bb[] (parents, lower index).  v==0 (ROOT_S) boundary:
 * only beta[0][L][L]=0. */
static void
ckpt_outside_deck(CKPT_CTX *cx, int v, float ***bb, int jp_0, int Lp_0)
{
  CM_t *cm = cx->cm;
  ESL_DSQ *dsq = cx->dsq;
  int   L = cx->L;
  int  *jmin = cx->jmin, *jmax = cx->jmax;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  float **bv = bb[v];
  float **esc_vAA_y;
  int j, d, i, y, voffset, jp_v, jp_y, dp_v, dp_y, sd, sdr, emitmode, jn, jx, dn, dx;
  float escore;

  ckpt_deck_init_impossible(cx, v, bv);

  if (v == 0) { bv[jp_0][Lp_0] = 0.; return; }

  /* R-L.2b local-begin seed: mirrors cm_OutsideAlignHB's upfront "for (v=1;
   * v<cm->M; v++) if (beginsc[v] valid) beta[v][L][L] = beginsc[v]" pass,
   * done once before the whole v-ascending recursion.  Here it's transcribed
   * per-deck instead (each v's deck is only ever (re)computed as a complete
   * unit), which is equivalent: the seed only touches v's OWN cell (L,L) and
   * has no cross-v ordering dependency, so seeding it immediately before
   * this deck's own parent-contribution accumulation (which FLogsums onto
   * whatever's already there) reproduces the oracle byte-for-byte. */
  if (cx->have_local_begin && NOT_IMPOSSIBLE(cm->beginsc[v]) && L >= jmin[v] && L <= jmax[v]) {
    jp_v = L - jmin[v];
    if (L >= hd_min(cp9b, v, jp_v) && L <= hd_max(cp9b, v, jp_v)) {
      int Lp = L - hd_min(cp9b, v, jp_v);
      bv[jp_v][Lp] = cm->beginsc[v];
    }
  }

  /* rung-3: BEGL_S / BEGR_S children of a pinned bifurcation.  Reads parent
   * beta from bb[] and sibling Inside alpha from cx->ifull[]; pinned to k*. */
  if (cm->stid[v] == BEGL_S || cm->stid[v] == BEGR_S) {
    float ***alpha = cx->ifull;
    int z, jp_z, k, kmin, kmax, kp_z;
    int yB = cm->plast[v];                          /* parent bifurcation */
    int kpinned = (cx->kpin != NULL) ? cx->kpin[yB] : -1;
    if (cm->stid[v] == BEGL_S) {
      y = yB; z = cm->cnum[yB];                     /* z = right (BEGR) sibling */
      for (j = jmax[v]; j >= jmin[v]; j--) {
        jp_v = j - jmin[v]; jp_y = j - jmin[y]; jp_z = j - jmin[z];
        for (d = hd_max(cp9b, v, jp_v); d >= hd_min(cp9b, v, jp_v); d--) {
          dp_v = d - hd_min(cp9b, v, jp_v);
          kmin = ESL_MAX(jmin[y], jmin[z]) - j;
          kmax = ESL_MIN(jmax[y], jmax[z]) - j;
          if (kpinned >= 0) { if (kpinned < kmin || kpinned > kmax) continue; kmin = kmax = kpinned; }
          for (k = kmin; k <= kmax; k++) {
            if (k < (hd_min(cp9b, y, jp_y+k) - d) || k > (hd_max(cp9b, y, jp_y+k) - d)) continue;
            if (k < (hd_min(cp9b, z, jp_z+k))     || k > (hd_max(cp9b, z, jp_z+k)))     continue;
            kp_z = k-hd_min(cp9b, z, jp_z+k);
            dp_y = d-hd_min(cp9b, y, jp_y+k);
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y+k][dp_y+k] + alpha[z][jp_z+k][kp_z]));
          }
        }
      }
    }
    else { /* BEGR_S */
      y = yB; z = cm->cfirst[yB];                   /* z = left (BEGL) sibling */
      jn = ESL_MAX(jmin[v], jmin[y]);
      jx = ESL_MIN(jmax[v], jmax[y]);
      for (j = jx; j >= jn; j--) {
        jp_v = j - jmin[v]; jp_y = j - jmin[y]; jp_z = j - jmin[z];
        dn = ESL_MAX(hd_min(cp9b, v, jp_v), j-jmax[z]);
        dx = ESL_MIN(hd_max(cp9b, v, jp_v), jp_z);
        for (d = dx; d >= dn; d--) {
          dp_v = d - hd_min(cp9b, v, jp_v);
          if (kpinned >= 0 && d != kpinned) continue;   /* BEGR own-d == right-frag length */
          kmin = ESL_MAX((hd_min(cp9b, y, jp_y)-d), (hd_min(cp9b, z, jp_z-d)));
          kmax = ESL_MIN((hd_max(cp9b, y, jp_y)-d), (hd_max(cp9b, z, jp_z-d)));
          for (k = kmin; k <= kmax; k++) {
            kp_z = k-hd_min(cp9b, z, jp_z-d);
            dp_y = d-hd_min(cp9b, y, jp_y);
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y][dp_y+k] + alpha[z][jp_z-d][kp_z]));
          }
        }
      }
    }
    return;
  }

  if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
    for (j = jmax[v]; j >= jmin[v]; j--) {
      jp_v = j - jmin[v];
      for (d = hd_max(cp9b, v, jp_v); d >= hd_min(cp9b, v, jp_v); d--) {
        i = j - d + 1;
        dp_v = d - hd_min(cp9b, v, jp_v);
        for (y = cm->plast[v]; y > cm->plast[v] - cm->pnum[v]; y--) {
          voffset = v - cm->cfirst[y];
          switch (cm->sttype[y]) {
          case MP_st: /* rung-3: MP parent of an IL/IR child */
            if (j == L || d == j) continue;
            if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue;
            jp_y = j - jmin[y];
            if ((d+2) < hd_min(cp9b, y, (jp_y+1)) || (d+2) > hd_max(cp9b, y, (jp_y+1))) continue;
            dp_y = d - hd_min(cp9b, y, jp_y+1);
            escore = cm->oesc[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y+1][dp_y+2] + cm->tsc[y][voffset] + escore));
            break;
          case ML_st:
          case IL_st:
            if (d == j) continue;
            if (j < jmin[y] || j > jmax[y]) continue;
            jp_y = j - jmin[y];
            if ((d+1) < hd_min(cp9b, y, jp_y) || (d+1) > hd_max(cp9b, y, jp_y)) continue;
            dp_y = d - hd_min(cp9b, y, jp_y);
            escore = cm->oesc[y][dsq[i-1]];
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y][dp_y+1] + cm->tsc[y][voffset] + escore));
            break;
          case MR_st:
          case IR_st:
            if (j == L) continue;
            if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue;
            jp_y = j - jmin[y];
            if ((d+1) < hd_min(cp9b, y, (jp_y+1)) || (d+1) > hd_max(cp9b, y, (jp_y+1))) continue;
            dp_y = d - hd_min(cp9b, y, (jp_y+1));
            escore = cm->oesc[y][dsq[j+1]];
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y+1][dp_y+1] + cm->tsc[y][voffset] + escore));
            break;
          case S_st:
          case E_st:
          case D_st:
            if (j < jmin[y] || j > jmax[y]) continue;
            jp_y = j - jmin[y];
            if (d < hd_min(cp9b, y, jp_y) || d > hd_max(cp9b, y, jp_y)) continue;
            dp_y = d - hd_min(cp9b, y, jp_y);
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y][dp_y] + cm->tsc[y][voffset]));
            break;
          }
        }
        if (bv[jp_v][dp_v] < IMPOSSIBLE) bv[jp_v][dp_v] = IMPOSSIBLE;
      }
    }
  }
  else {
    esc_vAA_y = cm->oesc;
    for (y = cm->plast[v]; y > cm->plast[v] - cm->pnum[v]; y--) {
      voffset = v - cm->cfirst[y];
      sdr = StateRightDelta(cm->sttype[y]);
      sd  = StateDelta(cm->sttype[y]);
      emitmode = Emitmode(cm->sttype[y]);
      jn = ESL_MAX(jmin[v], jmin[y] - sdr);
      jx = ESL_MIN(jmax[v], jmax[y] - sdr);
      for (j = jx; j >= jn; j--) {
        jp_v = j - jmin[v];
        jp_y = j - jmin[y];
        dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y + sdr) - sd);
        dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y + sdr) - sd);
        dp_v = dx - hd_min(cp9b, v, jp_v);
        dp_y = dx - hd_min(cp9b, y, jp_y + sdr);
        i    = j - dx + 1;
        switch (emitmode) {
        case EMITPAIR: /* rung-3: MP parent (pair emission) */
          for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) {
            escore = esc_vAA_y[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y + sdr][dp_y + sd] + cm->tsc[y][voffset] + escore));
          }
          break;
        case EMITLEFT:
          for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) {
            escore = esc_vAA_y[y][dsq[i-1]];
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y + sdr][dp_y + sd] + cm->tsc[y][voffset] + escore));
          }
          break;
        case EMITRIGHT:
          escore = esc_vAA_y[y][dsq[j+1]];
          for (d = dx; d >= dn; d--, dp_v--, dp_y--) {
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y + sdr][dp_y + sd] + cm->tsc[y][voffset] + escore));
          }
          break;
        case EMITNONE:
          for (d = dx; d >= dn; d--, dp_v--, dp_y--) {
            bv[jp_v][dp_v] = FLogsum(bv[jp_v][dp_v], (bb[y][jp_y + sdr][dp_y + sd] + cm->tsc[y][voffset]));
          }
          break;
        }
      }
    }
  }

  /* R-L.2 EL tail: deal with local-end transitions v->EL (EL deck = cx->elbeta,
   * indexed by ABSOLUTE (j,d)).  Mirrors cm_OutsideAlignHB:6356-6402 verbatim;
   * the only change is the banded target row (d<=eldmax[j], guaranteed to cover
   * the [dn,dx] write range).  The EL->EL self-transition + posterior fold are
   * done once, after the whole Step B sweep, in cm_CheckptAlignHB. */
  if (cx->have_el && NOT_IMPOSSIBLE(cm->endsc[v])) {
    int sd_v  = StateDelta(cm->sttype[v]);
    int sdr_v = StateRightDelta(cm->sttype[v]);
    int emm   = Emitmode(cm->sttype[v]);
    float **esc_vAA = cm->oesc;
    int jn2 = jmin[v] - sdr_v;
    int jx2 = jmax[v] - sdr_v;
    for (j = jn2; j <= jx2; j++) {
      jp_v = j - jmin[v];
      dn   = hd_min(cp9b, v, jp_v + sdr_v) - sd_v;
      dx   = hd_max(cp9b, v, jp_v + sdr_v) - sd_v;
      i    = j - dn + 1;
      dp_v = dn - hd_min(cp9b, v, jp_v + sdr_v);
      switch (emm) {
      case EMITPAIR:
        for (d = dn; d <= dx; d++, dp_v++, i--) {
          escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
          cx->elbeta[j][d] = FLogsum(cx->elbeta[j][d], (bb[v][jp_v+sdr_v][dp_v+sd_v] + cm->endsc[v] + escore));
        }
        break;
      case EMITLEFT:
        for (d = dn; d <= dx; d++, dp_v++, i--) {
          escore = esc_vAA[v][dsq[i-1]];
          cx->elbeta[j][d] = FLogsum(cx->elbeta[j][d], (bb[v][jp_v+sdr_v][dp_v+sd_v] + cm->endsc[v] + escore));
        }
        break;
      case EMITRIGHT:
        escore = esc_vAA[v][dsq[j+1]];
        for (d = dn; d <= dx; d++, dp_v++) {
          cx->elbeta[j][d] = FLogsum(cx->elbeta[j][d], (bb[v][jp_v+sdr_v][dp_v+sd_v] + cm->endsc[v] + escore));
        }
        break;
      case EMITNONE:
        for (d = dn; d <= dx; d++, dp_v++) {
          cx->elbeta[j][d] = FLogsum(cx->elbeta[j][d], (bb[v][jp_v+sdr_v][dp_v+sd_v] + cm->endsc[v]));
        }
        break;
      }
    }
  }
}

/* OptAcc deck v: mirrors cm_OptAccAlignHB's per-cell recurrence + FLogsum/max
 * order VERBATIM for S/IL/IR/ML/D/E, global (have_el=FALSE).  Reads emit
 * posteriors from cx->my_lpp / my_rpp; children OA-alpha from oa[] (block) or
 * ck[] (seeds).  Writes OA-alpha into oa[v] and (if ysh != NULL) yshadow. */
static void
ckpt_optacc_deck(CKPT_CTX *cx, int v, float ***oa, float ***ck, char **ysh)
{
#define OA(vv) (oa[vv] ? oa[vv] : (ck ? ck[vv] : NULL))
  CM_t *cm = cx->cm;
  int  *jmin = cx->jmin, *jmax = cx->jmax, *imin = cx->imin;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  float **av = oa[v];
  int sd  = StateDelta(cm->sttype[v]);
  int sdr = StateRightDelta(cm->sttype[v]);
  int j, d, i, y, yoffset, jp_v, dp_v, ip_v, jp_y_sdr, dp_y_sd, j_sdr;
  int yvalidA[MAXCONNECT], yvalid_ct, yvalid_idx;
  float sc;

  ckpt_deck_init_impossible(cx, v, av);
  if (ysh != NULL && cx->deck_nc[v] > 0) memset(ysh[0], (int) ((char) USED_EL), (size_t) cx->deck_nc[v]);

  if (cm->sttype[v] == E_st) {
    ckpt_apply_begin_oa(cx, v, av, ysh);
    return; /* OA: E cells remain IMPOSSIBLE */
  }

  /* d==0 (zero-length subtree) shadow init, cm_InitializeOptAccShadowDZeroHB:
   * route the d==sd cell to the unique StateDelta==0 child (the delete path) --
   * OR to USED_EL if a local end beats the all-deletes path (have_el branch,
   * :7698-7708). */
  if (ysh != NULL && cm->sttype[v] != B_st && cm->cp9b->Jvalid[v]) {
    int yy = cm->cfirst[v];
    while (StateDelta(cm->sttype[yy]) != 0) yy++;
    int yoff = yy - cm->cfirst[v];
    if (cx->have_el && cx->el_endsc > cx->el_esc[v]) yoff = (int) USED_EL;
    for (j = ESL_MAX(sd, jmin[v]); j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      if (hd_min(cp9b, v, jp_v) <= hd_max(cp9b, v, jp_v)) {
        if ((j - sdr) >= jmin[yy] && (j - sdr) <= jmax[yy]) {
          int jp_y = j - sdr - jmin[yy];
          if (sd >= hd_min(cp9b, v, jp_v) && sd <= hd_max(cp9b, v, jp_v) &&
              0  >= hd_min(cp9b, yy, jp_y) && 0 <= hd_max(cp9b, yy, jp_y)) {
            dp_v = sd - hd_min(cp9b, v, jp_v);
            ysh[jp_v][dp_v] = (char) yoff;
          }
        }
      }
    }
  }

  /* R-L.2 EL: re-init this state's OA deck from the EL prefix-sum (substitute,
   * don't omit -- the 51/053 segfault was omitting this).  alpha[v][j][d] =
   * alpha[cm->M][j-sdr][d-sd] = elalpha[j-sdr][d-sd]; yshadow stays USED_EL.
   * Mirrors cm_OptAccAlignHB:4319-4329. */
  if (cx->have_el && NOT_IMPOSSIBLE(cm->endsc[v])) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
        dp_v = d - hd_min(cp9b, v, jp_v);
        av[jp_v][dp_v] = cx->elalpha[j-sdr][d-sd];  /* ysh remains USED_EL */
      }
    }
  }

  if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      yvalid_ct = 0; j_sdr = j - sdr;
      for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++)
        if ((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset;
      i = j - hd_min(cp9b, v, jp_v) + 1;
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++, i--) {
        ip_v = i - imin[v];
        dp_v = d - hd_min(cp9b, v, jp_v);
        for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) {
          yoffset = yvalidA[yvalid_idx];
          y = cm->cfirst[v] + yoffset;
          jp_y_sdr = j - jmin[y] - sdr;
          if ((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) {
            dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
            if ((sc = OA(y)[jp_y_sdr][dp_y_sd]) > av[jp_v][dp_v]) {
              av[jp_v][dp_v] = sc;
              if (ysh != NULL) ysh[jp_v][dp_v] = (char) yoffset;
            }
          }
        }
        if (cm->sttype[v] == IL_st) av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], cx->my_lpp[v][ip_v]);
        else                        av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], cx->my_rpp[v][jp_v]);
        av[jp_v][dp_v] = ESL_MAX(av[jp_v][dp_v], IMPOSSIBLE);
        if ((! cx->have_el) && ysh != NULL && ysh[jp_v][dp_v] == (char) USED_EL && d > sd)
          av[jp_v][dp_v] = IMPOSSIBLE;
      }
    }
    ckpt_apply_begin_oa(cx, v, av, ysh);
    return;
  }
  else if (cm->sttype[v] == B_st) { /* rung-3: pinned bifurcation (single k*), OA = FLogsum of subtrees */
    int z = cm->cnum[v];            /* right (BEGR) subtree; y=cfirst[v]=left (BEGL) */
    int kpinned = (cx->kpin != NULL) ? cx->kpin[v] : -1;
    int jp_z, jp_y, kp_z, k, kn, kx, dp_y, jp_y_minus_k, dp_y_minus_k;
    y = cm->cfirst[v];
    int jnn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
    int jxx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
    for (j = jnn; j <= jxx; j++) {
      jp_v = j - jmin[v];
      jp_y = j - jmin[y];
      jp_z = j - jmin[z];
      kn = ((j-jmax[y]) > (hd_min(cp9b, z, jp_z))) ? (j-jmax[y]) : hd_min(cp9b, z, jp_z);
      kn = ESL_MAX(kn, 0);
      kx = ( jp_y       < (hd_max(cp9b, z, jp_z))) ?  jp_y       : hd_max(cp9b, z, jp_z);
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
        dp_v = d - hd_min(cp9b, v, jp_v);
        /* R-L.4b (065): if CYK's pin walk never visited this B state (kpinned<0)
         * AND more than one k satisfies the per-cell band constraint below for
         * THIS (j,d), do NOT fall back to an unrestricted full split-search --
         * see the matching comment in ckpt_inside_deck's B_st combine for the
         * full rationale (this is the OA-side half of the same fix; traceback
         * here likewise never stores a kshadow for the unpinned+ambiguous
         * case, so an untracked fallback value would be unreconstructable).
         * EXCEPTION: when exactly one k satisfies the constraint (checked
         * per-d, not just via the coarser j-only [kn,kx] window -- e.g. every
         * d==0 cell forces a unique k=0), the split is unambiguous and needs
         * no pin -- always compute it; ckpt_optacc_traceback re-derives the
         * same unique k the same way. */
        int klo, khi;
        if (kpinned >= 0) {
          if (kpinned < kn || kpinned > kx) continue;
          klo = khi = kpinned;
        } else {
          int uniq_k = -1, n_valid = 0;
          for (k = kn; k <= kx; k++) {
            if ((k >= d - hd_max(cp9b, y, jp_y-k)) && k <= d - hd_min(cp9b, y, jp_y-k)) {
              n_valid++; uniq_k = k;
              if (n_valid > 1) break;
            }
          }
          if (n_valid != 1) continue;
          klo = khi = uniq_k;
        }
        for (k = klo; k <= khi; k++) {
          if ((k >= d - hd_max(cp9b, y, jp_y-k)) && k <= d - hd_min(cp9b, y, jp_y-k)) {
            kp_z = k-hd_min(cp9b, z, jp_z);
            dp_y = d-hd_min(cp9b, y, jp_y-k);
            jp_y_minus_k = jp_y-k;
            dp_y_minus_k = dp_y-k;
            /* OA: accumulate (logsum) the posterior of both subtrees; record k
             * via the i15 validity guard (left/right IMPOSSIBLE allowed only at
             * length 0).  Traceback splits at kpin[v], so no kshadow is stored. */
            if ((sc = FLogsum(OA(y)[jp_y_minus_k][dp_y_minus_k], OA(z)[jp_z][kp_z])) > av[jp_v][dp_v]) {
              if (((d == k) || (NOT_IMPOSSIBLE(OA(y)[jp_y_minus_k][dp_y_minus_k]))) &&
                  ((k == 0) || (NOT_IMPOSSIBLE(OA(z)[jp_z][kp_z])))) {
                av[jp_v][dp_v] = sc;
              }
            }
          }
        }
      }
    }
    ckpt_apply_begin_oa(cx, v, av, ysh);
    return;
  }
  else { /* ML, MP, MR, D, S (non-self, non-B); E already returned */
    int jn, jx, jpn, jpx, dn, dx, dpn, dpx;
    for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
      yoffset = y - cm->cfirst[v];
      jn = ESL_MAX(jmin[v], jmin[y] + sdr);
      jx = ESL_MIN(jmax[v], jmax[y] + sdr);
      jpn = jn - jmin[v];
      jpx = jx - jmin[v];
      jp_y_sdr = jn - jmin[y] - sdr;
      for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
        dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y_sdr) + sd);
        dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y_sdr) + sd);
        dpn = dn - hd_min(cp9b, v, jp_v);
        dpx = dx - hd_min(cp9b, v, jp_v);
        dp_y_sd = dn - hd_min(cp9b, y, jp_y_sdr) - sd;
        for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) {
          if ((sc = OA(y)[jp_y_sdr][dp_y_sd]) > av[jp_v][dp_v]) {
            av[jp_v][dp_v] = sc;
            if (ysh != NULL) ysh[jp_v][dp_v] = (char) yoffset;
          }
        }
      }
    }
    switch (cm->sttype[v]) {        /* rung-3: ML (left), MR (right), MP (pair) posterior */
    case ML_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        i = j - hd_min(cp9b, v, jp_v) + 1;
        ip_v = i - imin[v];
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++, ip_v--)
          av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], cx->my_lpp[v][ip_v]);
      }
      break;
    case MR_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
          av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], cx->my_rpp[v][jp_v]);
      }
      break;
    case MP_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        i = j - hd_min(cp9b, v, jp_v) + 1;
        ip_v = i - imin[v];
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++, ip_v--)
          av[jp_v][dp_v] = FLogsum(av[jp_v][dp_v], FLogsum(cx->my_lpp[v][ip_v], cx->my_rpp[v][jp_v]));
      }
      break;
    default: break;
    }
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
        av[jp_v][dp_v] = ESL_MAX(av[jp_v][dp_v], IMPOSSIBLE);
    }
    if ((! cx->have_el) && sd > 0 && ysh != NULL) { /* emitters only (ML here); local: EL kept */
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        d = ESL_MAX(sd+1, hd_min(cp9b, v, jp_v));
        dp_v = d - hd_min(cp9b, v, jp_v);
        for (; d <= hd_max(cp9b, v, jp_v); d++, dp_v++)
          if (ysh[jp_v][dp_v] == (char) USED_EL) av[jp_v][dp_v] = IMPOSSIBLE;
      }
    }
    (void) jn; (void) jx; (void) jpn; (void) jpx; (void) dn; (void) dx; (void) dpn; (void) dpx;
    ckpt_apply_begin_oa(cx, v, av, ysh);
    return;
  }
#undef OA
}

/* Function: cm_CheckptAlignHB_Qualifies()
 * Purpose:  Return TRUE iff <cm> is a pure left-emitting MATL chain
 *           (0 B_st, no MP/MR, only S/IL/IR/ML/D/E).  These are the
 *           conditions under which the checkpointed engines
 *           (cm_CheckptAlignHB) reproduce the stock non-truncated OptAcc
 *           path byte-for-byte.  Local begins/ends are supported (R-L.2/2b:
 *           the engine consumes CMH_LOCAL_BEGIN/CMH_LOCAL_END via its
 *           have_local_begin/have_el context), so this gate does NOT reject
 *           local CMs.  The caller (DispatchSqAlignment) falls back to the
 *           stock path when FALSE.
 */
int
cm_CheckptAlignHB_Qualifies(CM_t *cm)
{
  int v;
  for (v = 0; v < cm->M; v++) {
    int st = cm->sttype[v];
    if (st == B_st || st == MP_st || st == MR_st) return FALSE;
    if (! (st==S_st || st==IL_st || st==IR_st || st==ML_st || st==D_st || st==E_st)) return FALSE;
  }
  return TRUE;
}

/* Function: cm_CheckptOptAccAlignHB_Qualifies()
 * Purpose:  Return TRUE iff <cm> is a STRUCTURED CM (>= 1 B/MP/MR) whose
 *           state set is the supported rung-3 surface (S/IL/IR/ML/MR/MP/D/E/B).
 *           These are the conditions under which the rung-3 pipeline
 *           (cm_CheckptPostAlignHB + cm_CheckptOptAccAlignHB) reproduces the
 *           stock NON-truncated OptAcc path's alignment (modulo the brief-26_0610-032
 *           accuracy-neutral pin-B flips), with a sqrt(M) working set.  Local
 *           begins/ends are supported (R-L.4/4b: the rung-3 engine consumes
 *           CMH_LOCAL_BEGIN/CMH_LOCAL_END), so this gate does NOT reject local
 *           CMs.  bps=0 CMs return FALSE here (they use the fused
 *           cm_CheckptAlignHB instead).  Deliberately SEPARATE from
 *           cm_CheckptAlignHB_Qualifies so the truncated qualifier
 *           (cm_CheckptTrAlignHB_Qualifies, which delegates to that one) is NOT
 *           relaxed -- truncated bps>0 --ckpt must fall back to the stock engine,
 *           not mis-route to the bps=0 truncated engine. */
int
cm_CheckptOptAccAlignHB_Qualifies(CM_t *cm)
{
  int v, has_bps = FALSE;
  for (v = 0; v < cm->M; v++) {
    int st = cm->sttype[v];
    if (st == B_st || st == MP_st || st == MR_st) has_bps = TRUE;
    if (! (st==S_st || st==IL_st || st==IR_st || st==ML_st || st==MR_st ||
           st==MP_st || st==D_st || st==E_st || st==B_st)) return FALSE;
  }
  return has_bps;
}

/* Function: cm_CheckptAlignHB()
 * Incept:   Brief 26_0610-028 (library port of drivers 022/023)
 *
 * Purpose:  Checkpointed (sqrt(M)-memory) HMM-banded optimal-accuracy
 *           alignment for a NON-truncated, GLOBAL, pure-MATL-chain CM.
 *           Drop-in replacement for the do_optacc=TRUE path of cm_AlignHB(),
 *           producing byte-identical Inside score, per-residue PPs, parsetree
 *           and avg PP, but with a sqrt(M) CM-DP working set instead of two
 *           full HMM-banded cubes.
 *
 *           Pipeline (each step byte-exact vs its stock counterpart):
 *             A : checkpointed Inside  -> Z, sqrt(M) alpha seed decks
 *             B : checkpointed Outside + fused posterior -> emit_mx l_pp/r_pp
 *             OA: checkpointed OptAcc max-DP -> OA score, sqrt(M) OA seed decks
 *             TB: checkpointed OptAcc traceback (block-recompute) -> parsetree
 *           Then cm_PostCodeHB() reads emit_mx for the PP string / avg PP.
 *
 *           cm->cp9b must already hold valid bands for <dsq>,<L> (the caller,
 *           DispatchSqAlignment, derives them).  Engagement must be gated by
 *           cm_CheckptAlignHB_Qualifies(); behavior is undefined otherwise.
 *
 *           Does NOT touch cm->hb_mx / hb_omx / hb_shmx.  Grows + fills the
 *           passed <emit_mx> (the deliverable).
 *
 * Args:     cm        - the covariance model (pure MATL chain, global)
 *           errbuf    - for error messages
 *           dsq       - digitized sequence 1..L
 *           L         - length of dsq
 *           size_limit- max Mb for emit_mx (passed to cm_hb_emit_mx_GrowTo)
 *           emit_mx   - emit matrix, grown + filled here
 *           ret_ppstr - RETURN: PP code string (NULL if not wanted)
 *           ret_tr    - RETURN: parsetree (NULL if not wanted)
 *           ret_avgpp - RETURN: avg PP of emitted residues
 *           ret_sc    - RETURN: alignment score in bits (Inside score Z)
 *
 * Returns:  <eslOK> on success.
 * Throws:   <eslEINCOMPAT> on contract violation (bands invalid for L);
 *           <eslERANGE> if emit_mx exceeds size_limit.
 */
int
cm_CheckptAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,
                  CM_HB_EMIT_MX *emit_mx, char **ret_ppstr, Parsetree_t **ret_tr,
                  float *ret_avgpp, float *ret_sc)
{
  int      status;
  CKPT_CTX cx;
  int      M = cm->M;
  int      v, k;
  float    Z_ckpt = 0.;
  Parsetree_t *tr = NULL;
  char    *ppstr = NULL;
  float    avgpp = 0.;
  float ***Astore = NULL, ***ba = NULL, ***bb = NULL, ***OAstore = NULL;
  float ***tba = NULL;
  char  ***tysh = NULL;

  /* contract: optacc requires emit_mx */
  if (emit_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptAlignHB(): emit_mx is NULL");

  /* zero the whole context first: cheap defensive reset so any field this
   * (or a future) rung doesn't explicitly set below defaults to 0/NULL
   * rather than reading stack garbage (R-L.2b). */
  memset(&cx, 0, sizeof(cx));

  /* set up the per-call context from cm->cp9b bands */
  cx.cm    = cm;
  cx.dsq   = dsq;
  cx.L     = L;
  cx.M     = M;
  cx.jmin  = cm->cp9b->jmin;  cx.jmax  = cm->cp9b->jmax;
  cx.imin  = cm->cp9b->imin;  cx.imax  = cm->cp9b->imax;
  cx.hdmin = cm->cp9b->hdmin; cx.hdmax = cm->cp9b->hdmax;
  CP9Bands_t *cp9b = cm->cp9b;
  cx.cur_bytes = cx.peak_bytes = 0;
  cx.deck_nc = NULL; cx.deck_njr = NULL; /* set below */
  cx.kpin = NULL; cx.ifull = NULL;       /* bps=0 path: no bifurcations */
  cx.have_el    = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;  /* R-L.2 */
  cx.el_selfsc  = cm->el_selfsc;
  cx.eldmax = NULL; cx.elbeta = NULL; cx.elalpha = NULL;
  cx.el_esc = NULL; cx.el_endsc = IMPOSSIBLE;
  /* R-L.2b: local BEGIN support */
  cx.have_local_begin = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;
  cx.bsc_fwd = IMPOSSIBLE; cx.fwd_begin_applied = FALSE;
  cx.begin_bsc = IMPOSSIBLE; cx.begin_b = -1;

  /* ROOT_S band sanity */
  if (cx.jmin[0] > L || cx.jmax[0] < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptAlignHB(): L outside ROOT_S j band");
  int jp_0 = L - cx.jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptAlignHB(): L outside ROOT_S d band");
  int Lp_0 = L - hd_min(cp9b, 0, jp_0);

  /* deck geometry + child/parent reach */
  ESL_ALLOC(cx.deck_nc,  sizeof(int64_t) * M);
  ESL_ALLOC(cx.deck_njr, sizeof(int)     * M);
  int64_t full_cube_cells = 0;
  int Delta = 0, Delta_p = 0;
  for (v = 0; v < M; v++) {
    int njr = cx.jmax[v] - cx.jmin[v] + 1; if (njr < 0) njr = 0;
    cx.deck_njr[v] = njr;
    int64_t nc = 0;
    int jp;
    for (jp = 0; jp < njr; jp++) { int w = hd_max(cp9b, v, jp)-hd_min(cp9b, v, jp)+1; if (w>0) nc += w; }
    cx.deck_nc[v] = nc;
    full_cube_cells += nc;
    if (cm->sttype[v] != E_st) {
      int ymax = cm->cfirst[v] + cm->cnum[v] - 1;
      if (ymax - v > Delta) Delta = ymax - v;
    }
    if (cm->pnum[v] > 0) {
      int ymin_par = cm->plast[v] - cm->pnum[v] + 1;
      if (v - ymin_par > Delta_p) Delta_p = v - ymin_par;
    }
  }
  int B = (int) (sqrt((double)M) + 0.5); if (B < 1) B = 1;

  /* R-L.2: EL setup (local mode only).  Compute the per-EL-row upper d-band
   * edge now; the banded EL decks are allocated lazily at the passes that need
   * them (elbeta at Step B, elalpha at Step OA). */
  if (cx.have_el) {
    ESL_ALLOC(cx.eldmax, sizeof(int) * (L+1));
    ckpt_el_compute_dmax(&cx);
    ESL_ALLOC(cx.el_esc, sizeof(float) * M);
    ckpt_el_compute_esc(&cx);
  }

  /* ============================================================= */
  /* STEP A: checkpointed Inside -> Z_ckpt + sqrt(M) seed store     */
  /* ============================================================= */
  ESL_ALLOC(Astore, sizeof(float**) * M);
  for (v = 0; v < M; v++) Astore[v] = NULL;
  for (v = M-1; v >= 0; v--) {
    Astore[v] = ckpt_deck_alloc(&cx, v);
    ckpt_inside_deck(&cx, v, Astore, NULL);
    int y = v + Delta;
    if (y < M && Astore[y] != NULL) {
      if ((y % B) < Delta) { /* retain as checkpoint seed */ }
      else { ckpt_deck_free(&cx, y, Astore[y]); Astore[y] = NULL; }
    }
  }
  Z_ckpt = Astore[0][jp_0][Lp_0];

  /* ============================================================= */
  /* STEP B: checkpointed Outside + fused posterior -> emit_mx      */
  /* ============================================================= */
  /* grow + initialize the emit matrix (the deliverable); alias into ctx */
  if ((status = cm_hb_emit_mx_GrowTo(cm, emit_mx, errbuf, cm->cp9b, L, size_limit)) != eslOK) goto ERROR;
  esl_vec_FSet(emit_mx->l_pp_mem, emit_mx->l_ncells_valid, IMPOSSIBLE);
  esl_vec_FSet(emit_mx->r_pp_mem, emit_mx->r_ncells_valid, IMPOSSIBLE);
  cx.my_lpp = emit_mx->l_pp;
  cx.my_rpp = emit_mx->r_pp;

  /* R-L.2: allocate the banded Outside EL deck (persists across the whole Step
   * B sweep, accumulating v->EL contributions from every local-end state). */
  if (cx.have_el) cx.elbeta = ckpt_el_deck_alloc(&cx);

  ESL_ALLOC(ba, sizeof(float**) * M);
  ESL_ALLOC(bb, sizeof(float**) * M);
  for (v = 0; v < M; v++) { ba[v] = NULL; bb[v] = NULL; }

  int nblocks = (M + B - 1) / B;
  for (k = 0; k < nblocks; k++) {
    int lo = k * B;
    int hi = ESL_MIN((k+1)*B - 1, M-1);
    /* recompute alpha for [lo..hi], descending; children in-block (ba) or seeds (Astore) */
    for (v = hi; v >= lo; v--) { ba[v] = ckpt_deck_alloc(&cx, v); ckpt_inside_deck(&cx, v, ba, Astore); }
    /* beta ascending through block, fold posterior into l_pp/r_pp, free alpha as we go */
    for (v = lo; v <= hi; v++) {
      int st = cm->sttype[v];
      bb[v] = ckpt_deck_alloc(&cx, v);
      ckpt_outside_deck(&cx, v, bb, jp_0, Lp_0);
      if (st==ML_st || st==IL_st) {
        int j;
        for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) {
          int jp = j - cx.jmin[v];
          int d;
          for (d = hd_min(cp9b, v, jp); d <= hd_max(cp9b, v, jp); d++) {
            int dp = d - hd_min(cp9b, v, jp);
            int i = j - d + 1;
            float postcell = ba[v][jp][dp] + bb[v][jp][dp] - Z_ckpt;
            int ip = i - cx.imin[v];
            emit_mx->l_pp[v][ip] = FLogsum(emit_mx->l_pp[v][ip], postcell);
          }
        }
      }
      if (st==IR_st) {
        int j;
        for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) {
          int jp = j - cx.jmin[v];
          int d;
          emit_mx->r_pp[v][jp] = ba[v][jp][0] + bb[v][jp][0] - Z_ckpt;
          for (d = hd_min(cp9b, v, jp)+1; d <= hd_max(cp9b, v, jp); d++) {
            int dp = d - hd_min(cp9b, v, jp);
            float postcell = ba[v][jp][dp] + bb[v][jp][dp] - Z_ckpt;
            emit_mx->r_pp[v][jp] = FLogsum(emit_mx->r_pp[v][jp], postcell);
          }
        }
      }
      ckpt_deck_free(&cx, v, ba[v]); ba[v] = NULL;
      int old = v - Delta_p - 1;
      if (old >= 0 && bb[old] != NULL) { ckpt_deck_free(&cx, old, bb[old]); bb[old] = NULL; }
    }
  }
  for (v = 0; v < M; v++) if (bb[v]) { ckpt_deck_free(&cx, v, bb[v]); bb[v] = NULL; }

  /* R-L.2 EmitterPosterior step 1 (EL): the EL->EL self-transition over the
   * accumulated elbeta, then fold the EL posterior into l_pp[cm->M] (1-D).
   * Mirrors cm_OutsideAlignHB:6408-6412 + cm_EmitterPosteriorHB:6961-6968.
   * post[cm->M][j][d] = el_scA[d] + elbeta[j][d] - Z; cells d>eldmax[j] are
   * IMPOSSIBLE (skipped, byte-exact: FLogsum(x,IMPOSSIBLE)==x). */
  if (cx.have_el) {
    int j, d;
    /* self-transition: descend d within each EL-row (band edge: top cell d=eldmax
     * keeps its value; its absent d+1 neighbour is IMPOSSIBLE, as in stock). */
    for (j = L; j >= 1; j--) {
      if (cx.eldmax[j] < 1) continue;
      for (d = cx.eldmax[j]-1; d >= 0; d--)
        cx.elbeta[j][d] = FLogsum(cx.elbeta[j][d], (cx.elbeta[j][d+1] + cx.el_selfsc));
    }
    /* fold to l_pp[cm->M][i], i=j-d+1, d>=1 (ascending j, ascending d == stock) */
    for (j = 1; j <= L; j++) {
      int dx = (cx.eldmax[j] < j) ? cx.eldmax[j] : j;
      int i = j;
      for (d = 1; d <= dx; d++, i--) {
        float postcell = (cx.el_selfsc * d) + cx.elbeta[j][d] - Z_ckpt;
        emit_mx->l_pp[cm->M][i] = FLogsum(emit_mx->l_pp[cm->M][i], postcell);
      }
    }
    ckpt_el_deck_free(&cx, cx.elbeta); cx.elbeta = NULL;  /* done with Outside EL */
  }

  /* EmitterPosterior step 2: normalize (mirror stock order exactly).  MATP
   * combine (step 3) doesn't apply to a MATL chain; EL (step 1) handled above. */
  esl_vec_FSet(emit_mx->sum, (L+1), IMPOSSIBLE);
  for (v = 0; v < M; v++) {
    if (emit_mx->l_pp[v] != NULL) { int i; for (i = cx.imin[v]; i <= cx.imax[v]; i++) { int ip=i-cx.imin[v]; emit_mx->sum[i]=FLogsum(emit_mx->sum[i], emit_mx->l_pp[v][ip]); } }
    if (emit_mx->r_pp[v] != NULL) { int j; for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) { int jp=j-cx.jmin[v]; emit_mx->sum[j]=FLogsum(emit_mx->sum[j], emit_mx->r_pp[v][jp]); } }
  }
  /* EL contributes to sum[i] LAST (mirror cm_EmitterPosteriorHB:6995-6999) */
  if (cx.have_el && emit_mx->l_pp[cm->M] != NULL) {
    int i; for (i = 1; i <= L; i++) emit_mx->sum[i] = FLogsum(emit_mx->sum[i], emit_mx->l_pp[cm->M][i]);
  }
  for (v = 0; v < M; v++) {
    if (emit_mx->l_pp[v] != NULL) { int i; for (i = cx.imin[v]; i <= cx.imax[v]; i++) { int ip=i-cx.imin[v]; emit_mx->l_pp[v][ip] -= emit_mx->sum[i]; } }
    if (emit_mx->r_pp[v] != NULL) { int j; for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) { int jp=j-cx.jmin[v]; emit_mx->r_pp[v][jp] -= emit_mx->sum[j]; } }
  }
  /* normalize EL row too (mirror cm_EmitterPosteriorHB:6832-6843 v==M case) */
  if (cx.have_el && emit_mx->l_pp[cm->M] != NULL) {
    int i; for (i = 1; i <= L; i++) emit_mx->l_pp[cm->M][i] -= emit_mx->sum[i];
  }

  /* R-L.2: build the OptAcc EL prefix-sum deck elalpha from the (normalized)
   * 1-D l_pp[cm->M] (mirror cm_OptAccAlignHB:4298-4306).  elalpha[j][d] =
   * logsum(l_pp[M][0], l_pp[M][j..j-d+1]); banded to 0..eldmax[j].  Persists
   * through Step OA + Step TB (read by ckpt_optacc_deck's EL re-init). */
  if (cx.have_el) {
    int j, d;
    cx.elalpha = ckpt_el_deck_alloc(&cx);
    for (j = 0; j <= L; j++) {
      if (cx.eldmax[j] < 0) continue;
      int i = j;
      cx.elalpha[j][0] = emit_mx->l_pp[cm->M][0];   /* = IMPOSSIBLE (i==0 never folded) */
      for (d = 1; d <= cx.eldmax[j]; d++)
        cx.elalpha[j][d] = FLogsum(cx.elalpha[j][d-1], emit_mx->l_pp[cm->M][i--]);
    }
  }

  /* ============================================================= */
  /* STEP OA: checkpointed OptAcc max-DP -> sqrt(M) OA seed decks   */
  /* ============================================================= */
  ESL_ALLOC(OAstore, sizeof(float**) * M);
  for (v = 0; v < M; v++) OAstore[v] = NULL;
  for (v = M-1; v >= 0; v--) {
    char **ysh = ckpt_cdeck_alloc(&cx, v);
    OAstore[v] = ckpt_deck_alloc(&cx, v);
    ckpt_optacc_deck(&cx, v, OAstore, NULL, ysh); /* children all in-window (ck=NULL) */
    ckpt_cdeck_free(&cx, v, ysh);                 /* yshadow only needed inside this call */
    int y = v + Delta;
    if (y < M && OAstore[y] != NULL) {
      if ((y % B) < Delta) { /* retain seed */ }
      else { ckpt_deck_free(&cx, y, OAstore[y]); OAstore[y] = NULL; }
    }
  }

  /* ============================================================= */
  /* STEP TB: checkpointed OptAcc traceback (block-recompute)      */
  /* ============================================================= */
  ESL_ALLOC(tba,  sizeof(float**) * M);
  ESL_ALLOC(tysh, sizeof(char**)  * M);
  for (v = 0; v < M; v++) { tba[v] = NULL; tysh[v] = NULL; }
  int cur_blk = -1, blk_lo = 0, blk_hi = -1;

  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0);
  {
    int i = 1, j = L, d = L, y;
    v = 0;
    while (1) {
      if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) break;
      int blk = v / B;
      if (blk != cur_blk) {
        int w;
        for (w = blk_lo; w <= blk_hi; w++) { if (tba[w]) { ckpt_deck_free(&cx,w,tba[w]); tba[w]=NULL; } if (tysh[w]) { ckpt_cdeck_free(&cx,w,tysh[w]); tysh[w]=NULL; } }
        cur_blk = blk; blk_lo = blk * B; blk_hi = ESL_MIN((blk+1)*B - 1, M-1);
        for (w = blk_hi; w >= blk_lo; w--) {
          tba[w]  = ckpt_deck_alloc(&cx, w);
          tysh[w] = ckpt_cdeck_alloc(&cx, w);
          ckpt_optacc_deck(&cx, w, tba, OAstore, tysh[w]);
        }
      }
      int jp_v = j - cx.jmin[v], dp_v = d - hd_min(cp9b, v, jp_v);
      int yoffset = tysh[v][jp_v][dp_v];
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default: ESL_XFAIL(eslEINCOMPAT, errbuf, "cm_CheckptAlignHB(): bogus state type in traceback v=%d", v);
      }
      d = j - i + 1;
      if (yoffset == (char) USED_EL) { InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M); v = cm->M; }
      else if (yoffset == (char) USED_LOCAL_BEGIN) { /* R-L.2b: local begin, can only happen once, from ROOT_S */
        InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cx.begin_b); v = cx.begin_b;
      }
      else { y = cm->cfirst[v] + yoffset; InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y); v = y; }
    }
    {
      int w;
      for (w = blk_lo; w <= blk_hi; w++) { if (tba[w]) { ckpt_deck_free(&cx,w,tba[w]); tba[w]=NULL; } if (tysh[w]) { ckpt_cdeck_free(&cx,w,tysh[w]); tysh[w]=NULL; } }
    }
  }

  /* per-residue PP string + avg PP from the emit matrix */
  if ((status = cm_PostCodeHB(cm, errbuf, L, emit_mx, tr, (ret_ppstr != NULL) ? &ppstr : NULL, &avgpp)) != eslOK) goto ERROR;

  if (getenv("INFERNAL_CKPT_VERBOSE")) {
    double full_mb = 2.0 * full_cube_cells * 4 / (1024.0*1024.0);
    double peak_mb = cx.peak_bytes / (1024.0*1024.0);
    fprintf(stderr, "# cm_CheckptAlignHB engaged: M=%d L=%d B=%d Z=%.4f  CM-DP peak=%.2f Mb  full-cube(2x)=%.2f Mb  win~%.1fx\n",
            M, L, B, Z_ckpt, peak_mb, full_mb, (peak_mb>0.) ? full_mb/peak_mb : 0.);
    if (cx.have_el && cx.eldmax != NULL) { /* R-L.2: EL-deck banding footprint */
      int64_t elcells = 0; int j;
      for (j = 0; j <= L; j++) if (cx.eldmax[j] >= 0) elcells += (cx.eldmax[j] + 1);
      double unbanded_mb = ((double)(L+1) * (double)(L+2) / 2.0) * 4 / (1024.0*1024.0); /* full O(L^2/2) per EL deck */
      double banded_mb   = (double) elcells * 4 / (1024.0*1024.0);
      fprintf(stderr, "#   EL-deck: banded=%.3f Mb/deck (%lld cells) vs unbanded O(L^2/2)=%.3f Mb/deck  reduction~%.1fx\n",
              banded_mb, (long long) elcells, unbanded_mb, (banded_mb>0.)? unbanded_mb/banded_mb : 0.);
    }
  }

  /* free seed stores + scratch */
  for (v = 0; v < M; v++) { if (Astore[v]) ckpt_deck_free(&cx, v, Astore[v]); if (OAstore[v]) ckpt_deck_free(&cx, v, OAstore[v]); }
  free(Astore); free(ba); free(bb); free(OAstore); free(tba); free(tysh);
  if (cx.elalpha) ckpt_el_deck_free(&cx, cx.elalpha);   /* R-L.2 EL */
  if (cx.elbeta)  ckpt_el_deck_free(&cx, cx.elbeta);    /* (freed earlier on success, here for safety) */
  if (cx.eldmax)  free(cx.eldmax);
  if (cx.el_esc)  free(cx.el_esc);
  free(cx.deck_nc); free(cx.deck_njr);

  if (ret_ppstr != NULL) *ret_ppstr = ppstr; else free(ppstr);
  if (ret_tr    != NULL) *ret_tr    = tr;    else FreeParsetree(tr);
  if (ret_avgpp != NULL) *ret_avgpp = avgpp;
  if (ret_sc    != NULL) *ret_sc    = Z_ckpt;
  return eslOK;

 ERROR:
  if (Astore)  { for (v = 0; v < M; v++) if (Astore[v])  ckpt_deck_free(&cx, v, Astore[v]);  free(Astore); }
  if (OAstore) { for (v = 0; v < M; v++) if (OAstore[v]) ckpt_deck_free(&cx, v, OAstore[v]); free(OAstore); }
  if (ba)   { for (v = 0; v < M; v++) if (ba[v])   ckpt_deck_free(&cx, v, ba[v]);   free(ba); }
  if (bb)   { for (v = 0; v < M; v++) if (bb[v])   ckpt_deck_free(&cx, v, bb[v]);   free(bb); }
  if (tba)  { for (v = 0; v < M; v++) if (tba[v])  ckpt_deck_free(&cx, v, tba[v]);  free(tba); }
  if (tysh) { for (v = 0; v < M; v++) if (tysh[v]) ckpt_cdeck_free(&cx, v, tysh[v]); free(tysh); }
  if (cx.elalpha) ckpt_el_deck_free(&cx, cx.elalpha);   /* R-L.2 EL */
  if (cx.elbeta)  ckpt_el_deck_free(&cx, cx.elbeta);
  if (cx.eldmax)  free(cx.eldmax);
  if (cx.el_esc)  free(cx.el_esc);
  if (cx.deck_nc)  free(cx.deck_nc);
  if (cx.deck_njr) free(cx.deck_njr);
  if (tr)    FreeParsetree(tr);
  if (ppstr) free(ppstr);
  return status;
}

/* is v a chain-root START state?  Its Inside root deck must survive from the
 * Inside pass to the Outside pass (read by its parent B's Inside recurrence and
 * by its sibling chain's BEGL_S/BEGR_S Outside recurrence). */
static int
ckpt_is_chain_root(CM_t *cm, int v)
{
  int s = cm->stid[v];
  return (s == ROOT_S || s == BEGL_S || s == BEGR_S);
}

/* Function: cm_PinPostAlignHB()
 * Incept:   Brief 26_0610-037 (rung-3 pinned posterior, milestone 1: full storage)
 *
 * Purpose:  Pinned (bps>0, structured-RNA) HMM-banded posterior, FULL storage.
 *           Builds the pinned Inside and Outside CM-DP matrices using the
 *           extended deck recurrences (MP/MR + single-k* bifurcation), then
 *           calls the STOCK cm_PosteriorHB() + cm_EmitterPosteriorHB() to
 *           produce <emit_mx> (l_pp/r_pp) and returns the Inside score Z.
 *
 *           This is the byte-exact-correctness milestone for the rung-3
 *           recurrence port: with identical <kpin> it must reproduce the
 *           non-checkpointed pinned reference cell-for-cell.  Memory is full
 *           (no checkpointing); cm_CheckptPostAlignHB() adds the sqrt(M) win.
 *
 *           Non-truncated only.  R-L.4: supports full local (EL+begin,
 *           have_el/have_local_begin from cm->flags); GLOBAL when the CM
 *           isn't in local mode.  Production dispatch still requires
 *           GLOBAL (cm_CheckptOptAccAlignHB_Qualifies() rejects local CMs
 *           -- relaxing that gate is R-L.6's job).
 *
 * Args:     cm, errbuf, dsq, L, size_limit - usual
 *           emit_mx - the deliverable, grown + filled here
 *           kpin    - [0..M-1] per-B pinned right-frag length k* (B states
 *                     only; -1 elsewhere).  NULL => unpinned (full sum).
 *           ret_sc  - RETURN: Inside score Z (bits)
 */
int
cm_PinPostAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,
                  CM_HB_EMIT_MX *emit_mx, int *kpin, float *ret_sc)
{
  int      status;
  CKPT_CTX cx;
  int      M = cm->M;
  int      v, jp;
  CM_HB_MX *imx = NULL, *omx = NULL;

  CP9Bands_t *cp9b = cm->cp9b;
  memset(&cx, 0, sizeof(cx));
  cx.cm = cm; cx.dsq = dsq; cx.L = L; cx.M = M;
  cx.jmin = cm->cp9b->jmin; cx.jmax = cm->cp9b->jmax;
  cx.imin = cm->cp9b->imin; cx.imax = cm->cp9b->imax;
  cx.hdmin = cm->cp9b->hdmin; cx.hdmax = cm->cp9b->hdmax;
  cx.cur_bytes = cx.peak_bytes = 0;
  cx.my_lpp = cx.my_rpp = NULL;
  cx.kpin = kpin; cx.ifull = NULL;
  cx.deck_nc = NULL; cx.deck_njr = NULL;
  /* R-L.4: local (EL+begin), mirrors cm_CheckptAlignHB's setup exactly */
  cx.have_el          = (cm->flags & CMH_LOCAL_END)   ? TRUE : FALSE;
  cx.el_selfsc        = cm->el_selfsc;
  cx.eldmax = NULL; cx.elbeta = NULL; cx.elalpha = NULL;
  cx.el_esc = NULL; cx.el_endsc = IMPOSSIBLE;
  cx.have_local_begin = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;
  cx.bsc_fwd = IMPOSSIBLE; cx.fwd_begin_applied = FALSE;
  cx.begin_bsc = IMPOSSIBLE; cx.begin_b = -1;

  if (cx.jmin[0] > L || cx.jmax[0] < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_PinPostAlignHB(): L outside ROOT_S j band");
  int jp_0 = L - cx.jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_PinPostAlignHB(): L outside ROOT_S d band");
  int Lp_0 = L - hd_min(cp9b, 0, jp_0);

  /* deck geometry (used by ckpt_deck_init_impossible inside the deck helpers) */
  ESL_ALLOC(cx.deck_nc,  sizeof(int64_t) * M);
  ESL_ALLOC(cx.deck_njr, sizeof(int)     * M);
  for (v = 0; v < M; v++) {
    int njr = cx.jmax[v] - cx.jmin[v] + 1; if (njr < 0) njr = 0;
    cx.deck_njr[v] = njr;
    int64_t nc = 0;
    for (jp = 0; jp < njr; jp++) { int w = hd_max(cp9b, v, jp)-hd_min(cp9b, v, jp)+1; if (w>0) nc += w; }
    cx.deck_nc[v] = nc;
  }

  if (cx.have_el) {
    ESL_ALLOC(cx.eldmax, sizeof(int) * (L+1));
    ckpt_el_compute_dmax(&cx);
    ESL_ALLOC(cx.el_esc, sizeof(float) * M);
    ckpt_el_compute_esc(&cx);
  }

  /* full Inside, into imx->dp (same banded layout the deck helpers assume) */
  imx = cm_hb_mx_Create(M);
  if ((status = cm_hb_mx_GrowTo(cm, imx, errbuf, cm->cp9b, L, size_limit)) != eslOK) goto ERROR;
  for (v = M-1; v >= 0; v--) ckpt_inside_deck(&cx, v, imx->dp, NULL);
  float Z = imx->dp[0][jp_0][Lp_0];

  /* R-L.4: fill alpha[cm->M] (EL "for completeness" ramp), mirrors
   * cm_InsideAlignHB:2966-2975 exactly -- cm_PosteriorHB() reads this
   * unbanded deck directly, so it must be materialized for real here (the
   * forward EL substitution in ckpt_inside_deck never touches deck M). */
  if (cx.have_el) {
    int j, d;
    for (j = 0; j <= L; j++)
      for (d = 0; d <= j; d++) imx->dp[cm->M][j][d] = cx.el_selfsc * d;
  }

  /* full Outside, into omx->dp; sibling Inside reads come from cx.ifull=imx->dp */
  omx = cm_hb_mx_Create(M);
  if ((status = cm_hb_mx_GrowTo(cm, omx, errbuf, cm->cp9b, L, size_limit)) != eslOK) goto ERROR;
  cx.ifull = imx->dp;
  if (cx.have_el) cx.elbeta = ckpt_el_deck_alloc(&cx);
  for (v = 0; v < M; v++) ckpt_outside_deck(&cx, v, omx->dp, jp_0, Lp_0);

  /* R-L.4: EL self-transition (mirrors cm_CheckptAlignHB's Step-B EL tail,
   * cm_dpalign.c:1967-1975 / cm_OutsideAlignHB:6408-6412) then materialize
   * the real beta[cm->M] deck cm_PosteriorHB() reads, banded cells from
   * cx.elbeta, everything else IMPOSSIBLE (provably correct per eldmax's
   * own banding argument). */
  if (cx.have_el) {
    int j, d;
    for (j = L; j >= 1; j--) {
      if (cx.eldmax[j] < 1) continue;
      for (d = cx.eldmax[j]-1; d >= 0; d--)
        cx.elbeta[j][d] = FLogsum(cx.elbeta[j][d], (cx.elbeta[j][d+1] + cx.el_selfsc));
    }
    for (j = 0; j <= L; j++) {
      for (d = 0; d <= j; d++) omx->dp[cm->M][j][d] = IMPOSSIBLE;
      if (cx.eldmax[j] >= 0) for (d = 0; d <= cx.eldmax[j]; d++) omx->dp[cm->M][j][d] = cx.elbeta[j][d];
    }
    ckpt_el_deck_free(&cx, cx.elbeta); cx.elbeta = NULL;
  }

  /* stock posterior + emitter (byte-exact with the reference, which uses these too) */
  if ((status = cm_PosteriorHB(cm, errbuf, L, size_limit, imx, omx, omx)) != eslOK) goto ERROR;
  if ((status = cm_EmitterPosteriorHB(cm, errbuf, L, size_limit, omx, emit_mx, FALSE)) != eslOK) goto ERROR;

  if (getenv("INFERNAL_CKPT_VERBOSE")) {
    int64_t cube=0; for (v=0;v<M;v++) cube += cx.deck_nc[v];
    int64_t ecells = emit_mx->l_ncells_valid + emit_mx->r_ncells_valid;
    /* largest bifurcation-free chain: contiguous non-B state run, cell count */
    int64_t chain=0, best=0; for (v=0;v<M;v++){ if (cm->sttype[v]==B_st){ if(chain>best)best=chain; chain=0;} else chain+=cx.deck_nc[v]; } if(chain>best)best=chain;
    fprintf(stderr, "# cm_PinPostAlignHB: M=%d L=%d Z=%.4f  full-IO-cube(2x)=%.2f Mb  emit_mx=%.2f Mb  largest-chain-cube=%.2f Mb (1x)\n",
            M, L, Z, 2.0*cube*4/(1024.0*1024.0), ecells*4.0/(1024.0*1024.0), best*4.0/(1024.0*1024.0));
  }

  cm_hb_mx_Destroy(imx); cm_hb_mx_Destroy(omx);
  free(cx.deck_nc); free(cx.deck_njr);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);
  if (ret_sc != NULL) *ret_sc = Z;
  return eslOK;

 ERROR:
  if (imx) cm_hb_mx_Destroy(imx);
  if (omx) cm_hb_mx_Destroy(omx);
  if (cx.elbeta) ckpt_el_deck_free(&cx, cx.elbeta);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);
  if (cx.deck_nc)  free(cx.deck_nc);
  if (cx.deck_njr) free(cx.deck_njr);
  return status;
}

/* Function: cm_CheckptPostAlignHB()
 * Incept:   Brief 26_0610-037 (rung-3 pinned posterior, milestone 2: sqrt(M) checkpointed)
 *
 * Purpose:  sqrt(M)-memory pinned posterior (the rung-3 deliverable).  Same
 *           result as cm_PinPostAlignHB() but with a sqrt(M)-bounded CM-DP
 *           working set per bifurcation-free chain (tree-of-chains).
 *
 *           Design (brief 26_0610-038, "global two-pass" tree-of-chains):
 *             STEP A : checkpointed Inside (descending sweep, exactly the bps=0
 *                      scheme) -> Z + sqrt(M) seed decks.  Two bifurcation tweaks:
 *                      (i) the linear child-reach Delta EXCLUDES B states
 *                      (cnum[B] is the BEGR index, not a child count -- 037 §9);
 *                      (ii) chain-root (START) Inside decks are NEVER auto-freed
 *                      -- they are read later by their parent B (Inside) and by
 *                      their sibling chain (Outside).  After STEP A, Astore holds
 *                      Z, all chain roots, and the global sqrt(M) seeds.
 *             STEP B : checkpointed Outside (ascending block sweep, bps=0 scheme)
 *                      + MATP-aware fused posterior into emit_mx.  Tweaks:
 *                      (i) the parent-reach Delta_p EXCLUDES B-parent edges
 *                      (BEGL_S/BEGR_S have a far B parent);
 *                      (ii) cx.ifull = Astore, so BEGL_S/BEGR_S Outside reads the
 *                      sibling's retained Inside root;
 *                      (iii) B-state beta decks are not auto-freed by the linear
 *                      rule -- they feed both child chains' Outside;
 *                      (iv) bounded co-floor reclaim: once a BEGR_S Outside is
 *                      done, its parent B-beta and BOTH child Inside roots are
 *                      consumed, so free them (bounds the retained-root/B-beta
 *                      co-floor to O(tree depth) instead of O(#bifurcations));
 *                      (v) the per-state posterior fold is extended to MP (folds
 *                      BOTH l_pp and r_pp) and MR (r_pp), and a final MATP l_pp/
 *                      r_pp merge (cm_EmitterPosteriorHB step 3) is reproduced.
 *
 *           Byte-exact vs cm_PinPostAlignHB() by construction: identical deck
 *           recurrences (ckpt_inside_deck / ckpt_outside_deck) and identical
 *           FLogsum/normalization order; only deck STORAGE differs.
 *
 *           Non-truncated only.  R-L.4: supports full local (EL+begin,
 *           have_el/have_local_begin from cm->flags); GLOBAL when the CM
 *           isn't in local mode.  Production dispatch still requires
 *           GLOBAL (cm_CheckptOptAccAlignHB_Qualifies() rejects local CMs
 *           -- relaxing that gate is R-L.6's job).
 */
int
cm_CheckptPostAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,
                      CM_HB_EMIT_MX *emit_mx, int *kpin, float *ret_sc)
{
  int      status;
  CKPT_CTX cx;
  int      M = cm->M;
  int      v, k;
  float    Z_ckpt = 0.;
  float ***Astore = NULL, ***ba = NULL, ***bb = NULL;

  CP9Bands_t *cp9b = cm->cp9b;
  memset(&cx, 0, sizeof(cx));
  if (emit_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptPostAlignHB(): emit_mx is NULL");

  cx.cm = cm; cx.dsq = dsq; cx.L = L; cx.M = M;
  cx.jmin  = cm->cp9b->jmin;  cx.jmax  = cm->cp9b->jmax;
  cx.imin  = cm->cp9b->imin;  cx.imax  = cm->cp9b->imax;
  cx.hdmin = cm->cp9b->hdmin; cx.hdmax = cm->cp9b->hdmax;
  cx.cur_bytes = cx.peak_bytes = 0;
  cx.deck_nc = NULL; cx.deck_njr = NULL;
  cx.kpin = kpin; cx.ifull = NULL;
  cx.my_lpp = cx.my_rpp = NULL;
  /* R-L.4: local (EL+begin), mirrors cm_CheckptAlignHB's setup exactly */
  cx.have_el          = (cm->flags & CMH_LOCAL_END)   ? TRUE : FALSE;
  cx.el_selfsc        = cm->el_selfsc;
  cx.eldmax = NULL; cx.elbeta = NULL; cx.elalpha = NULL;
  cx.el_esc = NULL; cx.el_endsc = IMPOSSIBLE;
  cx.have_local_begin = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;
  cx.bsc_fwd = IMPOSSIBLE; cx.fwd_begin_applied = FALSE;
  cx.begin_bsc = IMPOSSIBLE; cx.begin_b = -1;

  if (cx.jmin[0] > L || cx.jmax[0] < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptPostAlignHB(): L outside ROOT_S j band");
  int jp_0 = L - cx.jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptPostAlignHB(): L outside ROOT_S d band");
  int Lp_0 = L - hd_min(cp9b, 0, jp_0);

  /* deck geometry + linear child/parent reach (B-aware) */
  ESL_ALLOC(cx.deck_nc,  sizeof(int64_t) * M);
  ESL_ALLOC(cx.deck_njr, sizeof(int)     * M);
  int64_t full_cube_cells = 0, root_cells = 0, bstate_cells = 0;
  int Delta = 0, Delta_p = 0;
  for (v = 0; v < M; v++) {
    int njr = cx.jmax[v] - cx.jmin[v] + 1; if (njr < 0) njr = 0;
    cx.deck_njr[v] = njr;
    int64_t nc = 0; int jp;
    for (jp = 0; jp < njr; jp++) { int w = hd_max(cp9b, v, jp)-hd_min(cp9b, v, jp)+1; if (w>0) nc += w; }
    cx.deck_nc[v] = nc;
    full_cube_cells += nc;
    if (ckpt_is_chain_root(cm, v)) root_cells   += nc;
    if (cm->sttype[v] == B_st)     bstate_cells += nc;
    /* child reach: skip B (cnum[B] = BEGR index, not a count) */
    if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) {
      int ymax = cm->cfirst[v] + cm->cnum[v] - 1;
      if (ymax - v > Delta) Delta = ymax - v;
    }
    /* parent reach: skip B-parent edges (only BEGL_S/BEGR_S have a B parent) */
    if (cm->pnum[v] > 0 && cm->stid[v] != BEGL_S && cm->stid[v] != BEGR_S) {
      int ymin_par = cm->plast[v] - cm->pnum[v] + 1;
      if (v - ymin_par > Delta_p) Delta_p = v - ymin_par;
    }
  }
  int B = (int) (sqrt((double)M) + 0.5); if (B < 1) B = 1;

  if (cx.have_el) {
    ESL_ALLOC(cx.eldmax, sizeof(int) * (L+1));
    ckpt_el_compute_dmax(&cx);
    ESL_ALLOC(cx.el_esc, sizeof(float) * M);
    ckpt_el_compute_esc(&cx);
  }

  /* ============================================================= */
  /* STEP A: checkpointed Inside -> Z_ckpt + sqrt(M) seeds + roots  */
  /* ============================================================= */
  ESL_ALLOC(Astore, sizeof(float**) * M);
  for (v = 0; v < M; v++) Astore[v] = NULL;
  for (v = M-1; v >= 0; v--) {
    Astore[v] = ckpt_deck_alloc(&cx, v);
    ckpt_inside_deck(&cx, v, Astore, NULL); /* B reads its child roots from Astore (retained) */
    int y = v + Delta;
    if (y < M && Astore[y] != NULL && ! ckpt_is_chain_root(cm, y)) {
      if ((y % B) < Delta) { /* retain checkpoint seed */ }
      else { ckpt_deck_free(&cx, y, Astore[y]); Astore[y] = NULL; }
    }
  }
  Z_ckpt = Astore[0][jp_0][Lp_0];

  /* ============================================================= */
  /* STEP B: checkpointed Outside + MATP-aware fused posterior      */
  /* ============================================================= */
  if ((status = cm_hb_emit_mx_GrowTo(cm, emit_mx, errbuf, cm->cp9b, L, size_limit)) != eslOK) goto ERROR;
  esl_vec_FSet(emit_mx->l_pp_mem, emit_mx->l_ncells_valid, IMPOSSIBLE);
  esl_vec_FSet(emit_mx->r_pp_mem, emit_mx->r_ncells_valid, IMPOSSIBLE);
  cx.my_lpp = emit_mx->l_pp;
  cx.my_rpp = emit_mx->r_pp;
  cx.ifull  = Astore;  /* BEGL_S/BEGR_S Outside reads sibling Inside root from here */
  if (cx.have_el) cx.elbeta = ckpt_el_deck_alloc(&cx);

  ESL_ALLOC(ba, sizeof(float**) * M);
  ESL_ALLOC(bb, sizeof(float**) * M);
  for (v = 0; v < M; v++) { ba[v] = NULL; bb[v] = NULL; }

  int nblocks = (M + B - 1) / B;
  for (k = 0; k < nblocks; k++) {
    int lo = k * B;
    int hi = ESL_MIN((k+1)*B - 1, M-1);
    /* recompute alpha for [lo..hi], descending; children in-block (ba) or seeds/roots (Astore) */
    for (v = hi; v >= lo; v--) { ba[v] = ckpt_deck_alloc(&cx, v); ckpt_inside_deck(&cx, v, ba, Astore); }
    /* beta ascending; fold posterior into l_pp/r_pp; free alpha + spent betas/roots as we go */
    for (v = lo; v <= hi; v++) {
      int st = cm->sttype[v];
      bb[v] = ckpt_deck_alloc(&cx, v);
      ckpt_outside_deck(&cx, v, bb, jp_0, Lp_0);

      /* --- fused posterior, step 1 (mirror cm_EmitterPosteriorHB exactly) --- */
      if (st==MP_st || st==ML_st || st==IL_st) { /* leftwise emitter -> l_pp */
        int j;
        for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) {
          int jp = j - cx.jmin[v]; int d;
          for (d = hd_min(cp9b, v, jp); d <= hd_max(cp9b, v, jp); d++) {
            int dp = d - hd_min(cp9b, v, jp);
            int i  = j - d + 1;
            int ip = i - cx.imin[v];
            float postcell = ba[v][jp][dp] + bb[v][jp][dp] - Z_ckpt;
            emit_mx->l_pp[v][ip] = FLogsum(emit_mx->l_pp[v][ip], postcell);
          }
        }
      }
      if (st==MP_st || st==MR_st || st==IR_st) { /* rightwise emitter -> r_pp */
        int j;
        for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) {
          int jp = j - cx.jmin[v]; int d;
          emit_mx->r_pp[v][jp] = ba[v][jp][0] + bb[v][jp][0] - Z_ckpt; /* peel d=hdmin */
          for (d = hd_min(cp9b, v, jp)+1; d <= hd_max(cp9b, v, jp); d++) {
            int dp = d - hd_min(cp9b, v, jp);
            float postcell = ba[v][jp][dp] + bb[v][jp][dp] - Z_ckpt;
            emit_mx->r_pp[v][jp] = FLogsum(emit_mx->r_pp[v][jp], postcell);
          }
        }
      }

      /* free this state's recomputed alpha */
      ckpt_deck_free(&cx, v, ba[v]); ba[v] = NULL;
      /* free a non-B beta now out of linear parent reach */
      int old = v - Delta_p - 1;
      if (old >= 0 && bb[old] != NULL && cm->sttype[old] != B_st) { ckpt_deck_free(&cx, old, bb[old]); bb[old] = NULL; }
      /* bounded co-floor reclaim: once BOTH child chains of a bifurcation have had
       * their root-state Outside computed, the parent B-beta and both child Inside
       * roots have had their last read.  The two children are cfirst[yB] (BEGL) and
       * cnum[yB] (BEGR); the model does NOT guarantee BEGL < BEGR, so trigger on the
       * higher-indexed sibling (processed second in this ascending sweep). */
      if (cm->stid[v] == BEGL_S || cm->stid[v] == BEGR_S) {
        int yB   = cm->plast[v];   /* parent bifurcation               */
        int begl = cm->cfirst[yB]; /* left  child Inside root          */
        int begr = cm->cnum[yB];   /* right child Inside root          */
        int second = (begl > begr) ? begl : begr; /* sibling processed last */
        if (v == second) {
          if (bb[yB]       != NULL) { ckpt_deck_free(&cx, yB,   bb[yB]);       bb[yB]       = NULL; }
          if (Astore[begl] != NULL) { ckpt_deck_free(&cx, begl, Astore[begl]); Astore[begl] = NULL; }
          if (Astore[begr] != NULL) { ckpt_deck_free(&cx, begr, Astore[begr]); Astore[begr] = NULL; }
        }
      }
    }
  }
  for (v = 0; v < M; v++) if (bb[v]) { ckpt_deck_free(&cx, v, bb[v]); bb[v] = NULL; }

  /* R-L.4 EmitterPosterior step 1 (EL): EL->EL self-transition over the
   * accumulated elbeta, then fold into l_pp[cm->M] (1-D, bifurcation-
   * independent -- mirrors cm_CheckptAlignHB:1967-1985 verbatim). */
  if (cx.have_el) {
    int j, d;
    for (j = L; j >= 1; j--) {
      if (cx.eldmax[j] < 1) continue;
      for (d = cx.eldmax[j]-1; d >= 0; d--)
        cx.elbeta[j][d] = FLogsum(cx.elbeta[j][d], (cx.elbeta[j][d+1] + cx.el_selfsc));
    }
    for (j = 1; j <= L; j++) {
      int dx = (cx.eldmax[j] < j) ? cx.eldmax[j] : j;
      int i = j;
      for (d = 1; d <= dx; d++, i--) {
        float postcell = (cx.el_selfsc * d) + cx.elbeta[j][d] - Z_ckpt;
        emit_mx->l_pp[cm->M][i] = FLogsum(emit_mx->l_pp[cm->M][i], postcell);
      }
    }
    ckpt_el_deck_free(&cx, cx.elbeta); cx.elbeta = NULL;
  }

  /* EmitterPosterior step 2: normalize (mirror stock order exactly) */
  esl_vec_FSet(emit_mx->sum, (L+1), IMPOSSIBLE);
  for (v = 0; v < M; v++) {
    if (emit_mx->l_pp[v] != NULL) { int i; for (i = cx.imin[v]; i <= cx.imax[v]; i++) { int ip=i-cx.imin[v]; emit_mx->sum[i]=FLogsum(emit_mx->sum[i], emit_mx->l_pp[v][ip]); } }
    if (emit_mx->r_pp[v] != NULL) { int j; for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) { int jp=j-cx.jmin[v]; emit_mx->sum[j]=FLogsum(emit_mx->sum[j], emit_mx->r_pp[v][jp]); } }
  }
  /* EL contributes to sum[i] LAST (mirror cm_EmitterPosteriorHB:6995-6999) */
  if (cx.have_el && emit_mx->l_pp[cm->M] != NULL) {
    int i; for (i = 1; i <= L; i++) emit_mx->sum[i] = FLogsum(emit_mx->sum[i], emit_mx->l_pp[cm->M][i]);
  }
  for (v = 0; v < M; v++) {
    if (emit_mx->l_pp[v] != NULL) { int i; for (i = cx.imin[v]; i <= cx.imax[v]; i++) { int ip=i-cx.imin[v]; emit_mx->l_pp[v][ip] -= emit_mx->sum[i]; } }
    if (emit_mx->r_pp[v] != NULL) { int j; for (j = cx.jmin[v]; j <= cx.jmax[v]; j++) { int jp=j-cx.jmin[v]; emit_mx->r_pp[v][jp] -= emit_mx->sum[j]; } }
  }
  /* normalize EL row too (mirror cm_EmitterPosteriorHB:6832-6843 v==M case) */
  if (cx.have_el && emit_mx->l_pp[cm->M] != NULL) {
    int i; for (i = 1; i <= L; i++) emit_mx->l_pp[cm->M][i] -= emit_mx->sum[i];
  }

  /* EmitterPosterior step 3: combine l_pp for MATP_MP(v)/MATP_ML(v+1) and r_pp
   * for MATP_MP(v)/MATP_MR(v+2) in the same node (mirror stock exactly). */
  for (v = 0; v < M; v++) {
    if (cm->sttype[v] == MP_st) {
      if (cx.imax[v] >= 1 && cx.imax[v+1] >= 1) {
        int in = ESL_MAX(cx.imin[v], cx.imin[v+1]);
        int ix = ESL_MIN(cx.imax[v], cx.imax[v+1]);
        int i;
        for (i = in; i <= ix; i++) {
          int ip_v = i - cx.imin[v], ip_v2 = i - cx.imin[v+1];
          emit_mx->l_pp[v][ip_v]    = FLogsum(emit_mx->l_pp[v][ip_v], emit_mx->l_pp[v+1][ip_v2]);
          emit_mx->l_pp[v+1][ip_v2] = emit_mx->l_pp[v][ip_v];
        }
      }
      if (cx.jmax[v] >= 1 && cx.jmax[v+2] >= 1) {
        int jn = ESL_MAX(cx.jmin[v], cx.jmin[v+2]);
        int jx = ESL_MIN(cx.jmax[v], cx.jmax[v+2]);
        int j;
        for (j = jn; j <= jx; j++) {
          int jp_v = j - cx.jmin[v], jp_v2 = j - cx.jmin[v+2];
          emit_mx->r_pp[v][jp_v]    = FLogsum(emit_mx->r_pp[v][jp_v], emit_mx->r_pp[v+2][jp_v2]);
          emit_mx->r_pp[v+2][jp_v2] = emit_mx->r_pp[v][jp_v];
        }
      }
    }
  }

  if (getenv("INFERNAL_CKPT_VERBOSE")) {
    double  full_mb = 2.0 * full_cube_cells * 4 / (1024.0*1024.0);
    double  peak_mb = cx.peak_bytes / (1024.0*1024.0);
    int64_t ecells  = emit_mx->l_ncells_valid + emit_mx->r_ncells_valid;
    fprintf(stderr, "# cm_CheckptPostAlignHB: M=%d L=%d B=%d Z=%.4f  CM-DP peak=%.2f Mb  full-IO-cube(2x)=%.2f Mb  win~%.1fx  emit_mx=%.2f Mb  root-cells(1x)=%.2f Mb  Bstate-cells(1x)=%.2f Mb\n",
            M, L, B, Z_ckpt, peak_mb, full_mb, (peak_mb>0.) ? full_mb/peak_mb : 0.,
            ecells*4.0/(1024.0*1024.0), root_cells*4.0/(1024.0*1024.0), bstate_cells*4.0/(1024.0*1024.0));
  }

  for (v = 0; v < M; v++) if (Astore[v]) ckpt_deck_free(&cx, v, Astore[v]);
  free(Astore); free(ba); free(bb);
  free(cx.deck_nc); free(cx.deck_njr);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);
  if (ret_sc != NULL) *ret_sc = Z_ckpt;
  return eslOK;

 ERROR:
  if (Astore) { for (v = 0; v < M; v++) if (Astore[v]) ckpt_deck_free(&cx, v, Astore[v]); free(Astore); }
  if (ba) { for (v = 0; v < M; v++) if (ba[v]) ckpt_deck_free(&cx, v, ba[v]); free(ba); }
  if (bb) { for (v = 0; v < M; v++) if (bb[v]) ckpt_deck_free(&cx, v, bb[v]); free(bb); }
  if (cx.elbeta) ckpt_el_deck_free(&cx, cx.elbeta);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);
  if (cx.deck_nc)  free(cx.deck_nc);
  if (cx.deck_njr) free(cx.deck_njr);
  return status;
}

/*****************************************************************
 * rung-3 OptAcc + pinned-tree traceback (brief 26_0610-039)
 *
 * The OptAcc max-DP runs over the SAME extended deck recurrence
 * (ckpt_optacc_deck, now with MP/MR + single-k* B), reading the emit_mx
 * posteriors that cm_PinPostAlignHB / cm_CheckptPostAlignHB produce.  The
 * pinned-tree traceback is pin-B's payoff (030 sec Q3.3): the stock
 * B-traceback's only use of the 3D shadow is reading one integer k
 * (cm_alignT_hb), and kpin[v] already supplies k*.  So the traceback is a
 * deterministic descent of the known bifurcation tree (split each B at k*) +
 * the per-chain linear yshadow descent -- no OA-outside, no 3D OA kshadow.
 *****************************************************************/

/* yshadow fetch callback: returns yshadow[v][jp_v][dp_v].  Two implementations:
 * full-storage (cm_PinOptAccAlignHB) just indexes; checkpointed
 * (cm_CheckptOptAccAlignHB) block-recomputes the OA deck + yshadow on demand. */
typedef char (*ckpt_ysh_fetch_fn)(void *ctx, int v, int jp_v, int dp_v);

/* full-storage fetch: all yshadow decks resident */
typedef struct { char ***ysh; } pin_ysh_ctx;
static char
pin_ysh_fetch(void *p, int v, int jp_v, int dp_v)
{
  pin_ysh_ctx *c = (pin_ysh_ctx *) p;
  return c->ysh[v][jp_v][dp_v];
}

/* checkpointed fetch: block-recompute the OA deck + yshadow for v's block,
 * reading children from OAstore (retained chain roots + sqrt(M) seeds).
 * Handles non-monotonic visiting (the bifurcation DFS jumps around) by
 * reloading whenever v's block differs from the loaded block. */
typedef struct {
  CKPT_CTX *cx; int M, B;
  float ***OAstore;      /* retained roots + sqrt(M) OA seeds */
  float ***tba; char ***tysh;
  int cur_blk, blk_lo, blk_hi;
} ckpt_ysh_ctx;
static char
ckpt_ysh_fetch(void *p, int v, int jp_v, int dp_v)
{
  ckpt_ysh_ctx *c = (ckpt_ysh_ctx *) p;
  int blk = v / c->B;
  if (blk != c->cur_blk) {
    int w;
    for (w = c->blk_lo; w <= c->blk_hi; w++) {
      if (c->tba[w])  { ckpt_deck_free (c->cx, w, c->tba[w]);  c->tba[w]  = NULL; }
      if (c->tysh[w]) { ckpt_cdeck_free(c->cx, w, c->tysh[w]); c->tysh[w] = NULL; }
    }
    c->cur_blk = blk; c->blk_lo = blk * c->B; c->blk_hi = ESL_MIN((blk+1)*c->B - 1, c->M-1);
    for (w = c->blk_hi; w >= c->blk_lo; w--) {
      c->tba[w]  = ckpt_deck_alloc(c->cx, w);
      c->tysh[w] = (c->cx->cm->sttype[w] == B_st) ? NULL : ckpt_cdeck_alloc(c->cx, w);
      ckpt_optacc_deck(c->cx, w, c->tba, c->OAstore, c->tysh[w]);
    }
  }
  return c->tysh[v][jp_v][dp_v];
}

/* Shared pinned-tree OptAcc traceback.  Mirrors stock cm_alignT_hb's descent +
 * bifurcation stack exactly, EXCEPT a B splits at kpin[v] (not a kshadow read)
 * and yshadow is obtained via <fetch>.  GLOBAL mode (every B is on the parse, so
 * every B is pinned by the CYK parse). */
static int
ckpt_optacc_traceback(CM_t *cm, char *errbuf, int L, int *kpin,
                      int *jmin, int *jmax, CP9Bands_t *cp9b,
                      ckpt_ysh_fetch_fn fetch, void *fctx, int begin_b, Parsetree_t **ret_tr)
{
  int status;
  Parsetree_t *tr  = NULL;
  ESL_STACK   *pda = NULL;
  int v = 0, i = 1, j = L, d = L, k, y, yoffset, bifparent;
  int jp_v = 0, dp_v = 0, allow_S_local_end;

  tr = CreateParsetree(100);
  if (tr == NULL) { status = eslEMEM; goto ERROR; }
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0);
  pda = esl_stack_ICreate();
  if (pda == NULL) { status = eslEMEM; goto ERROR; }

  while (1) {
    /* OptAcc d==0 BEGL_S/BEGR_S special case (mirror cm_alignT_hb): a zero-length
     * START subtree whose cell is outside the bands -> route to USED_EL. */
    if (d == 0 && (cm->stid[v] == BEGL_S || cm->stid[v] == BEGR_S) &&
        ((j < jmin[v] || j > jmax[v]) ||
         (d < hd_min(cp9b, v, j-jmin[v]) || d > hd_max(cp9b, v, j-jmin[v])))) {
      allow_S_local_end = TRUE;
    }
    else if (cm->sttype[v] != EL_st) {
      jp_v = j - jmin[v];
      dp_v = d - hd_min(cp9b, v, jp_v);
      allow_S_local_end = FALSE;
    }

    if (cm->sttype[v] == B_st) {
      k = (kpin != NULL) ? kpin[v] : -1;     /* pinned right-fragment length */
      if (k < 0) {
        /* R-L.4b (065): CYK's pin walk never visited this B state.  If exactly
         * one k satisfies the per-cell band constraint for THIS (j,d), this
         * isn't really "unpinned" -- it's forced by band geometry alone,
         * exactly mirroring the ckpt_*_deck B_st combine's own n_valid==1
         * exception (checked per-d, not just via the coarser j-only [kn,kx]
         * window -- e.g. every d==0 cell forces a unique k=0 even when
         * [kn,kx] itself spans many values).  Re-derive it the same way the
         * deck did rather than failing. */
        int zz = cm->cnum[v], yy = cm->cfirst[v];
        if (j >= jmin[yy] && j <= jmax[yy] && j >= jmin[zz] && j <= jmax[zz]) {
          int jp_yy = j - jmin[yy], jp_zz = j - jmin[zz];
          int kn2 = ESL_MAX(ESL_MAX(j - jmax[yy], hd_min(cp9b, zz, jp_zz)), 0);
          int kx2 = ESL_MIN(jp_yy, hd_max(cp9b, zz, jp_zz));
          int kk, uniq_k = -1, n_valid = 0;
          for (kk = kn2; kk <= kx2; kk++) {
            if ((kk >= d - hd_max(cp9b, yy, jp_yy-kk)) && (kk <= d - hd_min(cp9b, yy, jp_yy-kk))) {
              n_valid++; uniq_k = kk;
              if (n_valid > 1) break;
            }
          }
          if (n_valid == 1) k = uniq_k;
        }
      }
      if (k < 0) ESL_XFAIL(eslEINCOMPAT, errbuf, "ckpt_optacc_traceback: B state v=%d not pinned (global mode expected)", v);
      if ((status = esl_stack_IPush(pda, j))       != eslOK) goto ERROR;
      if ((status = esl_stack_IPush(pda, k))       != eslOK) goto ERROR;
      if ((status = esl_stack_IPush(pda, tr->n-1)) != eslOK) goto ERROR;
      j = j - k;
      d = d - k;
      i = j - d + 1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    }
    else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;  /* traceback complete */
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];
      y = cm->cnum[v];                       /* right START state */
      i = j - d + 1;
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
    }
    else {
      if (allow_S_local_end) yoffset = USED_EL;
      else                   yoffset = fetch(fctx, v, jp_v, dp_v);
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default: ESL_XFAIL(eslEINVAL, errbuf, "ckpt_optacc_traceback: bogus state type v=%d", v);
      }
      d = j - i + 1;
      if (yoffset == (char) USED_EL) {
        InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
        v = cm->M;
      }
      else if (yoffset == (char) USED_LOCAL_BEGIN) { /* R-L.4: local begin, can only happen once, from ROOT_S */
        InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, begin_b);
        v = begin_b;
      }
      else {
        y = cm->cfirst[v] + yoffset;
        InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
        v = y;
      }
    }
  }
  esl_stack_Destroy(pda);
  *ret_tr = tr;
  return eslOK;

 ERROR:
  if (pda) esl_stack_Destroy(pda);
  if (tr)  FreeParsetree(tr);
  return status;
}

/* Function: cm_PinOptAccAlignHB()
 * Incept:   Brief 26_0610-039 (rung-3 OptAcc, milestone 1: full storage)
 *
 * Purpose:  Full-storage pinned OptAcc alignment (the byte-exact correctness
 *           anchor for the sqrt(M) wrap).  Runs the extended OA max-DP
 *           (ckpt_optacc_deck, with MP/MR + single-k* B) into full per-state
 *           decks reading the pre-filled <emit_mx>, then a pinned-tree
 *           traceback -> parsetree + per-residue PP.  GLOBAL, non-truncated.
 */
int
cm_PinOptAccAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,
                    CM_HB_EMIT_MX *emit_mx, int *kpin,
                    char **ret_ppstr, Parsetree_t **ret_tr, float *ret_avgpp, float *ret_pp)
{
  int      status;
  CKPT_CTX cx;
  int      M = cm->M;
  int      v, jp;
  float  ***OA  = NULL;
  char   ***ysh = NULL;
  Parsetree_t *tr = NULL;
  char    *ppstr  = NULL;
  float    avgpp  = 0., pp = 0.;
  pin_ysh_ctx fctx;

  if (emit_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_PinOptAccAlignHB(): emit_mx is NULL");

  CP9Bands_t *cp9b = cm->cp9b;
  memset(&cx, 0, sizeof(cx));
  cx.cm = cm; cx.dsq = dsq; cx.L = L; cx.M = M;
  cx.jmin = cm->cp9b->jmin; cx.jmax = cm->cp9b->jmax;
  cx.imin = cm->cp9b->imin; cx.imax = cm->cp9b->imax;
  cx.hdmin = cm->cp9b->hdmin; cx.hdmax = cm->cp9b->hdmax;
  cx.cur_bytes = cx.peak_bytes = 0;
  cx.kpin = kpin; cx.ifull = NULL;        /* OA reads child OA decks, not Inside */
  cx.my_lpp = emit_mx->l_pp; cx.my_rpp = emit_mx->r_pp;
  cx.deck_nc = NULL; cx.deck_njr = NULL;
  /* R-L.4: local (EL+begin), mirrors cm_CheckptAlignHB's setup exactly */
  cx.have_el          = (cm->flags & CMH_LOCAL_END)   ? TRUE : FALSE;
  cx.el_selfsc        = cm->el_selfsc;
  cx.eldmax = NULL; cx.elbeta = NULL; cx.elalpha = NULL;
  cx.el_esc = NULL; cx.el_endsc = IMPOSSIBLE;
  cx.have_local_begin = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;
  cx.bsc_fwd = IMPOSSIBLE; cx.fwd_begin_applied = FALSE;
  cx.begin_bsc = IMPOSSIBLE; cx.begin_b = -1;

  if (cx.jmin[0] > L || cx.jmax[0] < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_PinOptAccAlignHB(): L outside ROOT_S j band");
  int jp_0 = L - cx.jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_PinOptAccAlignHB(): L outside ROOT_S d band");
  int Lp_0 = L - hd_min(cp9b, 0, jp_0);

  ESL_ALLOC(cx.deck_nc,  sizeof(int64_t) * M);
  ESL_ALLOC(cx.deck_njr, sizeof(int)     * M);
  for (v = 0; v < M; v++) {
    int njr = cx.jmax[v] - cx.jmin[v] + 1; if (njr < 0) njr = 0;
    cx.deck_njr[v] = njr;
    int64_t nc = 0;
    for (jp = 0; jp < njr; jp++) { int w = hd_max(cp9b, v, jp)-hd_min(cp9b, v, jp)+1; if (w>0) nc += w; }
    cx.deck_nc[v] = nc;
  }

  if (cx.have_el) {
    ESL_ALLOC(cx.eldmax, sizeof(int) * (L+1));
    ckpt_el_compute_dmax(&cx);
    ESL_ALLOC(cx.el_esc, sizeof(float) * M);
    ckpt_el_compute_esc(&cx);
    /* R-L.4: OptAcc EL prefix-sum deck, built from the (normalized) 1-D
     * l_pp[cm->M] the preceding posterior pass filled -- mirrors
     * cm_CheckptAlignHB:2008-2022 exactly (bifurcation-independent: this
     * is a pure per-residue prefix sum, not read/written by the B-combine). */
    int j, d;
    cx.elalpha = ckpt_el_deck_alloc(&cx);
    for (j = 0; j <= L; j++) {
      if (cx.eldmax[j] < 0) continue;
      int ii = j;
      cx.elalpha[j][0] = cx.my_lpp[cm->M][0];
      for (d = 1; d <= cx.eldmax[j]; d++)
        cx.elalpha[j][d] = FLogsum(cx.elalpha[j][d-1], cx.my_lpp[cm->M][ii--]);
    }
  }

  /* full OA max-DP (v = M-1 .. 0) + per-state yshadow (NULL for B; B uses kpin) */
  ESL_ALLOC(OA,  sizeof(float**) * M);
  ESL_ALLOC(ysh, sizeof(char**)  * M);
  for (v = 0; v < M; v++) { OA[v] = NULL; ysh[v] = NULL; }
  for (v = 0; v < M; v++) {
    OA[v]  = ckpt_deck_alloc(&cx, v);
    ysh[v] = (cm->sttype[v] == B_st) ? NULL : ckpt_cdeck_alloc(&cx, v);
  }
  for (v = M-1; v >= 0; v--) ckpt_optacc_deck(&cx, v, OA, NULL, ysh[v]);

  pp = sreEXP2(OA[0][jp_0][Lp_0]) / (float) L;

  /* pinned-tree traceback (full-storage yshadow fetch) */
  fctx.ysh = ysh;
  if ((status = ckpt_optacc_traceback(cm, errbuf, L, kpin, cx.jmin, cx.jmax, cp9b,
                                      pin_ysh_fetch, &fctx, cx.begin_b, &tr)) != eslOK) goto ERROR;

  /* per-residue PP string + avg PP from the emit matrix */
  if ((status = cm_PostCodeHB(cm, errbuf, L, emit_mx, tr, (ret_ppstr != NULL) ? &ppstr : NULL, &avgpp)) != eslOK) goto ERROR;

  if (getenv("INFERNAL_CKPT_VERBOSE")) {
    int64_t cube=0; for (v=0;v<M;v++) cube += cx.deck_nc[v];
    fprintf(stderr, "# cm_PinOptAccAlignHB: M=%d L=%d pp=%.5f  full-OA-cube(1x)=%.2f Mb (+ yshadow ~%.2f Mb)\n",
            M, L, pp, cube*4.0/(1024.0*1024.0), cube*1.0/(1024.0*1024.0));
  }

  for (v = 0; v < M; v++) { if (OA[v]) ckpt_deck_free(&cx, v, OA[v]); if (ysh[v]) ckpt_cdeck_free(&cx, v, ysh[v]); }
  free(OA); free(ysh);
  free(cx.deck_nc); free(cx.deck_njr);
  if (cx.elalpha) ckpt_el_deck_free(&cx, cx.elalpha);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);

  if (ret_ppstr != NULL) *ret_ppstr = ppstr; else free(ppstr);
  if (ret_tr    != NULL) *ret_tr    = tr;    else FreeParsetree(tr);
  if (ret_avgpp != NULL) *ret_avgpp = avgpp;
  if (ret_pp    != NULL) *ret_pp    = pp;
  return eslOK;

 ERROR:
  if (OA)  { for (v = 0; v < M; v++) if (OA[v])  ckpt_deck_free(&cx, v, OA[v]);  free(OA); }
  if (ysh) { for (v = 0; v < M; v++) if (ysh[v]) ckpt_cdeck_free(&cx, v, ysh[v]); free(ysh); }
  if (cx.elalpha) ckpt_el_deck_free(&cx, cx.elalpha);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);
  if (cx.deck_nc)  free(cx.deck_nc);
  if (cx.deck_njr) free(cx.deck_njr);
  if (tr)    FreeParsetree(tr);
  if (ppstr) free(ppstr);
  return status;
}

/* Function: cm_CheckptOptAccAlignHB()
 * Incept:   Brief 26_0610-039 (rung-3 OptAcc, milestone 2: sqrt(M) checkpointed)
 *
 * Purpose:  sqrt(M)-memory pinned OptAcc alignment (the rung-3 standalone
 *           aligner's max-DP half).  Same parsetree + PP as cm_PinOptAccAlignHB
 *           but with a sqrt(M)-bounded OA working set per bifurcation-free chain.
 *
 *           STEP OA : checkpointed OA max-DP (descending sweep, exactly the
 *                     bps=0 scheme + the 038 STEP-A bifurcation retention) ->
 *                     OAstore = {all chain roots} + {global sqrt(M) seeds}.
 *                     Chain-root OA decks are never auto-freed (a parent B reads
 *                     its child roots; the traceback block-recompute reads them).
 *           STEP TB : pinned-tree traceback (ckpt_optacc_traceback) with the
 *                     yshadow block-recomputed on demand from OAstore.  B splits
 *                     at kpin[v]; no OA-outside, no 3D OA kshadow.
 *
 *           Byte-exact vs cm_PinOptAccAlignHB() by construction: identical OA
 *           deck recurrence (ckpt_optacc_deck) and identical descent; only deck
 *           STORAGE differs.  GLOBAL, non-truncated only.
 */
int
cm_CheckptOptAccAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,
                        CM_HB_EMIT_MX *emit_mx, int *kpin,
                        char **ret_ppstr, Parsetree_t **ret_tr, float *ret_avgpp, float *ret_pp)
{
  int      status;
  CKPT_CTX cx;
  int      M = cm->M;
  int      v;
  float  ***OAstore = NULL, ***tba = NULL;
  char   ***tysh = NULL;
  Parsetree_t *tr = NULL;
  char    *ppstr  = NULL;
  float    avgpp  = 0., pp = 0.;
  ckpt_ysh_ctx fctx;

  CP9Bands_t *cp9b = cm->cp9b;
  memset(&cx, 0, sizeof(cx));
  if (emit_mx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptOptAccAlignHB(): emit_mx is NULL");

  cx.cm = cm; cx.dsq = dsq; cx.L = L; cx.M = M;
  cx.jmin = cm->cp9b->jmin; cx.jmax = cm->cp9b->jmax;
  cx.imin = cm->cp9b->imin; cx.imax = cm->cp9b->imax;
  cx.hdmin = cm->cp9b->hdmin; cx.hdmax = cm->cp9b->hdmax;
  cx.cur_bytes = cx.peak_bytes = 0;
  cx.kpin = kpin; cx.ifull = NULL;
  cx.my_lpp = emit_mx->l_pp; cx.my_rpp = emit_mx->r_pp;
  cx.deck_nc = NULL; cx.deck_njr = NULL;
  /* R-L.4: local (EL+begin), mirrors cm_CheckptAlignHB's setup exactly */
  cx.have_el          = (cm->flags & CMH_LOCAL_END)   ? TRUE : FALSE;
  cx.el_selfsc        = cm->el_selfsc;
  cx.eldmax = NULL; cx.elbeta = NULL; cx.elalpha = NULL;
  cx.el_esc = NULL; cx.el_endsc = IMPOSSIBLE;
  cx.have_local_begin = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;
  cx.bsc_fwd = IMPOSSIBLE; cx.fwd_begin_applied = FALSE;
  cx.begin_bsc = IMPOSSIBLE; cx.begin_b = -1;

  if (cx.jmin[0] > L || cx.jmax[0] < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptOptAccAlignHB(): L outside ROOT_S j band");
  int jp_0 = L - cx.jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptOptAccAlignHB(): L outside ROOT_S d band");
  int Lp_0 = L - hd_min(cp9b, 0, jp_0);

  /* deck geometry + B-aware linear child reach (037/038 Delta-fix) */
  ESL_ALLOC(cx.deck_nc,  sizeof(int64_t) * M);
  ESL_ALLOC(cx.deck_njr, sizeof(int)     * M);
  int64_t full_cube_cells = 0, root_cells = 0;
  int Delta = 0;
  for (v = 0; v < M; v++) {
    int njr = cx.jmax[v] - cx.jmin[v] + 1; if (njr < 0) njr = 0;
    cx.deck_njr[v] = njr;
    int64_t nc = 0; int jp;
    for (jp = 0; jp < njr; jp++) { int w = hd_max(cp9b, v, jp)-hd_min(cp9b, v, jp)+1; if (w>0) nc += w; }
    cx.deck_nc[v] = nc;
    full_cube_cells += nc;
    if (ckpt_is_chain_root(cm, v)) root_cells += nc;
    if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) {
      int ymax = cm->cfirst[v] + cm->cnum[v] - 1;
      if (ymax - v > Delta) Delta = ymax - v;
    }
  }
  int B = (int) (sqrt((double)M) + 0.5); if (B < 1) B = 1;

  if (cx.have_el) {
    ESL_ALLOC(cx.eldmax, sizeof(int) * (L+1));
    ckpt_el_compute_dmax(&cx);
    ESL_ALLOC(cx.el_esc, sizeof(float) * M);
    ckpt_el_compute_esc(&cx);
    /* R-L.4: OptAcc EL prefix-sum deck (mirror cm_CheckptAlignHB:2008-2022) */
    int j, d;
    cx.elalpha = ckpt_el_deck_alloc(&cx);
    for (j = 0; j <= L; j++) {
      if (cx.eldmax[j] < 0) continue;
      int ii = j;
      cx.elalpha[j][0] = cx.my_lpp[cm->M][0];
      for (d = 1; d <= cx.eldmax[j]; d++)
        cx.elalpha[j][d] = FLogsum(cx.elalpha[j][d-1], cx.my_lpp[cm->M][ii--]);
    }
  }

  /* ============================================================= */
  /* STEP OA: checkpointed OA max-DP -> roots + sqrt(M) OA seeds    */
  /* ============================================================= */
  ESL_ALLOC(OAstore, sizeof(float**) * M);
  for (v = 0; v < M; v++) OAstore[v] = NULL;
  for (v = M-1; v >= 0; v--) {
    char **ysh = (cm->sttype[v] == B_st) ? NULL : ckpt_cdeck_alloc(&cx, v);
    OAstore[v] = ckpt_deck_alloc(&cx, v);
    ckpt_optacc_deck(&cx, v, OAstore, NULL, ysh);  /* B reads child roots from OAstore (retained) */
    if (ysh) ckpt_cdeck_free(&cx, v, ysh);          /* forward-sweep yshadow not kept (recomputed in TB) */
    int y = v + Delta;
    if (y < M && OAstore[y] != NULL && ! ckpt_is_chain_root(cm, y)) {
      if ((y % B) < Delta) { /* retain checkpoint seed */ }
      else { ckpt_deck_free(&cx, y, OAstore[y]); OAstore[y] = NULL; }
    }
  }
  pp = sreEXP2(OAstore[0][jp_0][Lp_0]) / (float) L;

  /* ============================================================= */
  /* STEP TB: pinned-tree traceback (block-recompute yshadow)       */
  /* ============================================================= */
  ESL_ALLOC(tba,  sizeof(float**) * M);
  ESL_ALLOC(tysh, sizeof(char**)  * M);
  for (v = 0; v < M; v++) { tba[v] = NULL; tysh[v] = NULL; }
  fctx.cx = &cx; fctx.M = M; fctx.B = B; fctx.OAstore = OAstore;
  fctx.tba = tba; fctx.tysh = tysh;
  fctx.cur_blk = -1; fctx.blk_lo = 0; fctx.blk_hi = -1;

  if ((status = ckpt_optacc_traceback(cm, errbuf, L, kpin, cx.jmin, cx.jmax, cp9b,
                                      ckpt_ysh_fetch, &fctx, cx.begin_b, &tr)) != eslOK) goto ERROR;

  /* free the loaded traceback block */
  { int w; for (w = fctx.blk_lo; w <= fctx.blk_hi; w++) {
      if (tba[w])  { ckpt_deck_free (&cx, w, tba[w]);  tba[w]  = NULL; }
      if (tysh[w]) { ckpt_cdeck_free(&cx, w, tysh[w]); tysh[w] = NULL; } } }

  if ((status = cm_PostCodeHB(cm, errbuf, L, emit_mx, tr, (ret_ppstr != NULL) ? &ppstr : NULL, &avgpp)) != eslOK) goto ERROR;

  if (getenv("INFERNAL_CKPT_VERBOSE")) {
    double full_mb = full_cube_cells * 4 / (1024.0*1024.0);   /* OA is 1 cube (vs posterior 2) */
    double peak_mb = cx.peak_bytes / (1024.0*1024.0);
    int64_t ecells = emit_mx->l_ncells_valid + emit_mx->r_ncells_valid;
    fprintf(stderr, "# cm_CheckptOptAccAlignHB: M=%d L=%d B=%d pp=%.5f  OA-DP peak=%.2f Mb  full-OA-cube(1x)=%.2f Mb  win~%.1fx  emit_mx=%.2f Mb  root-cells(1x)=%.2f Mb\n",
            M, L, B, pp, peak_mb, full_mb, (peak_mb>0.) ? full_mb/peak_mb : 0.,
            ecells*4.0/(1024.0*1024.0), root_cells*4.0/(1024.0*1024.0));
  }

  for (v = 0; v < M; v++) if (OAstore[v]) ckpt_deck_free(&cx, v, OAstore[v]);
  free(OAstore); free(tba); free(tysh);
  free(cx.deck_nc); free(cx.deck_njr);
  if (cx.elalpha) ckpt_el_deck_free(&cx, cx.elalpha);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);

  if (ret_ppstr != NULL) *ret_ppstr = ppstr; else free(ppstr);
  if (ret_tr    != NULL) *ret_tr    = tr;    else FreeParsetree(tr);
  if (ret_avgpp != NULL) *ret_avgpp = avgpp;
  if (ret_pp    != NULL) *ret_pp    = pp;
  return eslOK;

 ERROR:
  if (OAstore) { for (v = 0; v < M; v++) if (OAstore[v]) ckpt_deck_free(&cx, v, OAstore[v]); free(OAstore); }
  if (tba)  { for (v = 0; v < M; v++) if (tba[v])  ckpt_deck_free(&cx, v, tba[v]);  free(tba); }
  if (tysh) { for (v = 0; v < M; v++) if (tysh[v]) ckpt_cdeck_free(&cx, v, tysh[v]); free(tysh); }
  if (cx.elalpha) ckpt_el_deck_free(&cx, cx.elalpha);
  if (cx.eldmax) free(cx.eldmax);
  if (cx.el_esc) free(cx.el_esc);
  if (cx.deck_nc)  free(cx.deck_nc);
  if (cx.deck_njr) free(cx.deck_njr);
  if (tr)    FreeParsetree(tr);
  if (ppstr) free(ppstr);
  return status;
}

/*****************************************************************
 * Brief 26_0610-078 R1: checkpointed, bifurcation (k*)-DISCOVERING
 * CYK max-DP, GLOBAL/non-truncated/no-EL/no-local-begin only.
 *
 * Small-scale proof of the core NEW piece brief 076 (design pass) grounded:
 * a B_st combine that SEARCHES k on already-resident checkpoint decks,
 * instead of reading a pre-supplied kpin[v] like ckpt_optacc_deck's B_st
 * branch does.  Transcribes ckpt_optacc_deck's checkpointing SHAPE (deck
 * alloc/free, chain-root retention via ckpt_is_chain_root, block-recompute
 * traceback) but swaps the OA/FLogsum recurrence for cm_CYKInsideAlignHB's
 * real max/+ recurrence (cm_dpalign.c's cm_CYKInsideAlignHB, ~3443-3835;
 * B_st k-loop ~3736-3799, the transcription source for the search bounds
 * below).  Because a B_st's children (y=BEGL_S, z=BEGR_S) are themselves
 * always chain roots (ckpt_is_chain_root), their decks are NEVER freed --
 * so unlike ckpt_optacc_deck (kpin supplied once, no kshadow stored), this
 * engine doesn't need to persist a k* per B-state cell either: the same
 * bounded search is deterministic and cheap to simply re-run at traceback
 * time against the still-resident child decks (ckpt_cyk_bsearch(), shared
 * by both the forward sweep and the traceback).
 *
 * R1 is scoped tight (076's own assessment: "arguably skippable... the
 * exact mechanism is already proven in production, just not yet used for
 * k-*discovery*") -- no EL, no local begin, GLOBAL only.  R2 (combined-mode
 * TRUNCATED, J/L/R/T, full local support) is the actual rung-4 deliverable.
 *****************************************************************/

/* B_st k-search: given already-resident child CYK decks (via CY(), the same
 * cy[]-then-ck[]-fallback pattern ckpt_optacc_deck's OA() macro uses),
 * find k* maximizing cy[y][j-k][d-k] + cy[z][j][k] for state v's (j,d) cell.
 * Bounds transcribed verbatim from cm_CYKInsideAlignHB's B_st loop
 * (cm_dpalign.c ~3744-3797).  Shared by ckpt_cyk_deck (forward fill) and
 * ckpt_cyk_traceback (B-state split) -- same search, re-run each time
 * rather than storing a kshadow, since child decks are always resident. */
static int
ckpt_cyk_bsearch(CM_t *cm, int *jmin, int *jmax, CP9Bands_t *cp9b,
                  float ***cy, float ***ck, int v, int j, int d, float *ret_sc)
{
#define CY(vv) (cy[vv] ? cy[vv] : (ck ? ck[vv] : NULL))
  int y = cm->cfirst[v], z = cm->cnum[v];
  int jp_y = j - jmin[y], jp_z = j - jmin[z];
  int kn = ESL_MAX(ESL_MAX(j-jmax[y], hd_min(cp9b, z, jp_z)), 0);
  int kx = ESL_MIN(jp_y, hd_max(cp9b, z, jp_z));
  int k, kbest = -1; float scbest = IMPOSSIBLE, sc;
  for (k = kn; k <= kx; k++) {
    if ((k >= d - hd_max(cp9b, y, jp_y-k)) && (k <= d - hd_min(cp9b, y, jp_y-k))) {
      int kp_z = k - hd_min(cp9b, z, jp_z);
      int dp_y = d - hd_min(cp9b, y, jp_y-k);
      if ((sc = CY(y)[jp_y-k][dp_y-k] + CY(z)[jp_z][kp_z]) > scbest) { scbest = sc; kbest = k; }
    }
  }
  if (ret_sc) *ret_sc = scbest;
  return kbest;
#undef CY
}

/* CYK deck v: mirrors cm_CYKInsideAlignHB's per-cell recurrence VERBATIM for
 * S/IL/IR/ML/MR/MP/D/E (real tsc_v[]+oesc_v[], not posteriors/FLogsum) --
 * global only (no EL reinit, no local-begin reduction; R1's scope).  For
 * B_st, calls ckpt_cyk_bsearch() instead of reading a pre-supplied kpin[v]
 * -- the new piece.  Writes CYK-alpha into cy[v] and (if ysh != NULL,
 * non-B states only) yshadow (yoffset, mirrors cm_CYKInsideAlignHB's
 * yshadow -- NOT a kshadow; B-state k is re-derived, never stored). */
static void
ckpt_cyk_deck(CKPT_CTX *cx, int v, float ***cy, float ***ck, char **ysh)
{
#define CY(vv) (cy[vv] ? cy[vv] : (ck ? ck[vv] : NULL))
  CM_t *cm = cx->cm;
  ESL_DSQ *dsq = cx->dsq;
  int  *jmin = cx->jmin, *jmax = cx->jmax;
  CP9Bands_t *cp9b = cx->cm->cp9b;
  float **av = cy[v];
  float const *esc_v = cm->oesc[v];
  float const *tsc_v = cm->tsc[v];
  int sd  = StateDelta(cm->sttype[v]);
  int sdr = StateRightDelta(cm->sttype[v]);
  int j, d, i, y, yoffset, jp_v, dp_v, jp_y_sdr, dp_y_sd, j_sdr;
  int yvalidA[MAXCONNECT], yvalid_ct, yvalid_idx;
  float sc;

  ckpt_deck_init_impossible(cx, v, av);

  /* R1-L (084) EL: yshadow default = USED_EL, so any cell the child recurrence
   * never overwrites traces back to the local end (mirrors cm_CYKInsideAlignHB
   * :3954).  Gated on have_el to leave the validated GLOBAL path byte-identical
   * (global never sets/reads a USED_EL shadow). */
  if (cx->have_el && ysh != NULL && cx->deck_nc[v] > 0)
    memset(ysh[0], (int) ((char) USED_EL), (size_t) cx->deck_nc[v]);

  /* R1-L (084) EL: re-init this state's CYK deck if a local end from v is
   * allowed.  Forward EL score is the FIXED ramp el_selfsc*(d-sd) plus the
   * v->EL transition cm->endsc[v] (mirrors cm_CYKInsideAlignHB:3980-4001; this
   * closed form equals that function's el_scA[d-sd]+endsc[v] for d>=sd, with no
   * alpha[cm->M] EL deck read).  Placed before the per-state recurrence so the
   * child max compares against, and emissions add onto, this EL base -- exactly
   * as stock CYK.  yshadow stays USED_EL wherever EL wins. */
  if (cx->have_el && NOT_IMPOSSIBLE(cm->endsc[v])) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      for (dp_v = 0, d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); dp_v++, d++)
        if (d >= sd) av[jp_v][dp_v] = cx->el_selfsc * (d - sd) + cm->endsc[v];
    }
  }

  if (cm->sttype[v] == E_st) {
    for (j = jmin[v]; j <= jmax[v]; j++) { jp_v = j - jmin[v]; av[jp_v][0] = 0.; }
    ckpt_apply_begin_cyk(cx, v, av, ysh);
    return;
  }

  if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      yvalid_ct = 0; j_sdr = j - sdr;
      for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++)
        if ((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset;
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
        i = j - d + 1;
        dp_v = d - hd_min(cp9b, v, jp_v);
        for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) {
          yoffset = yvalidA[yvalid_idx];
          y = cm->cfirst[v] + yoffset;
          jp_y_sdr = j - jmin[y] - sdr;
          if ((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) {
            dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
            if ((sc = CY(y)[jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > av[jp_v][dp_v]) {
              av[jp_v][dp_v] = sc;
              if (ysh != NULL) ysh[jp_v][dp_v] = (char) yoffset;
            }
          }
        }
        av[jp_v][dp_v] += (cm->sttype[v] == IL_st) ? esc_v[dsq[i]] : esc_v[dsq[j]];
        av[jp_v][dp_v] = ESL_MAX(av[jp_v][dp_v], IMPOSSIBLE);
      }
    }
    ckpt_apply_begin_cyk(cx, v, av, ysh);
    return;
  }
  else if (cm->sttype[v] == B_st) {
    int z = cm->cnum[v];
    y = cm->cfirst[v];
    int jnn = ESL_MAX(jmin[v], jmin[z]);
    int jxx = ESL_MIN(jmax[v], jmax[z]);
    for (j = jnn; j <= jxx; j++) {
      jp_v = j - jmin[v];
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
        dp_v = d - hd_min(cp9b, v, jp_v);
        sc = IMPOSSIBLE;
        if (ckpt_cyk_bsearch(cm, jmin, jmax, cp9b, cy, ck, v, j, d, &sc) >= 0)
          av[jp_v][dp_v] = sc;
      }
    }
    ckpt_apply_begin_cyk(cx, v, av, ysh);
    return;
  }
  else { /* ML, MR, MP, D, S (non-self, non-B); E already returned */
    int jn, jx, jpn, jpx, dn, dx, dpn, dpx;
    for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
      yoffset = y - cm->cfirst[v];
      float tsc = tsc_v[yoffset];
      jn = ESL_MAX(jmin[v], jmin[y]+sdr);
      jx = ESL_MIN(jmax[v], jmax[y]+sdr);
      jpn = jn - jmin[v];
      jpx = jx - jmin[v];
      jp_y_sdr = jn - jmin[y] - sdr;
      for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
        dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y_sdr) + sd);
        dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y_sdr) + sd);
        dpn = dn - hd_min(cp9b, v, jp_v);
        dpx = dx - hd_min(cp9b, v, jp_v);
        dp_y_sd = dn - hd_min(cp9b, y, jp_y_sdr) - sd;
        for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) {
          if ((sc = CY(y)[jp_y_sdr][dp_y_sd] + tsc) > av[jp_v][dp_v]) {
            av[jp_v][dp_v] = sc;
            if (ysh != NULL) ysh[jp_v][dp_v] = (char) yoffset;
          }
        }
      }
    }
    switch (cm->sttype[v]) {
    case ML_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v]; i = j - hd_min(cp9b, v, jp_v) + 1;
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++, i--)
          av[jp_v][dp_v] += esc_v[dsq[i]];
      }
      break;
    case MR_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v];
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
          av[jp_v][dp_v] += esc_v[dsq[j]];
      }
      break;
    case MP_st:
      for (j = jmin[v]; j <= jmax[v]; j++) {
        jp_v = j - jmin[v]; i = j - hd_min(cp9b, v, jp_v) + 1;
        for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++, i--)
          av[jp_v][dp_v] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
      }
      break;
    default: break;
    }
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
        av[jp_v][dp_v] = ESL_MAX(av[jp_v][dp_v], IMPOSSIBLE);
    }
    ckpt_apply_begin_cyk(cx, v, av, ysh);
    return;
  }
#undef CY
}

/* checkpointed fetch: block-recompute the CYK deck + yshadow for v's block,
 * reading children from CYstore (retained chain roots + sqrt(M) seeds).
 * Mirrors ckpt_ysh_fetch/ckpt_ysh_ctx (cm_dpalign.c ~2693-2718) verbatim,
 * swapping ckpt_optacc_deck for ckpt_cyk_deck. */
typedef struct {
  CKPT_CTX *cx; int M, B;
  float ***CYstore;
  float ***tba; char ***tysh;
  int cur_blk, blk_lo, blk_hi;
} ckpt_cyk_ysh_ctx;
static char
ckpt_cyk_ysh_fetch(void *p, int v, int jp_v, int dp_v)
{
  ckpt_cyk_ysh_ctx *c = (ckpt_cyk_ysh_ctx *) p;
  int blk = v / c->B;
  if (blk != c->cur_blk) {
    int w;
    for (w = c->blk_lo; w <= c->blk_hi; w++) {
      if (c->tba[w])  { ckpt_deck_free (c->cx, w, c->tba[w]);  c->tba[w]  = NULL; }
      if (c->tysh[w]) { ckpt_cdeck_free(c->cx, w, c->tysh[w]); c->tysh[w] = NULL; }
    }
    c->cur_blk = blk; c->blk_lo = blk * c->B; c->blk_hi = ESL_MIN((blk+1)*c->B - 1, c->M-1);
    for (w = c->blk_hi; w >= c->blk_lo; w--) {
      c->tba[w]  = ckpt_deck_alloc(c->cx, w);
      c->tysh[w] = (c->cx->cm->sttype[w] == B_st) ? NULL : ckpt_cdeck_alloc(c->cx, w);
      ckpt_cyk_deck(c->cx, w, c->tba, c->CYstore, c->tysh[w]);
    }
  }
  return c->tysh[v][jp_v][dp_v];
}

/* Checkpointed CYK traceback.  Mirrors ckpt_optacc_traceback's descent +
 * bifurcation stack (cm_dpalign.c ~2724-2839), except a B_st splits at a k
 * found via ckpt_cyk_bsearch() (live search against the retained CYstore
 * decks) instead of a pre-supplied kpin[v].  Brief 26_0610-084 (R1-L) adds
 * EL (USED_EL -> descend to cm->M) and local-begin (USED_LOCAL_BEGIN ->
 * descend to begin_b) handling, transcribed from cm_alignT_hb's own CYK
 * (do_optacc=FALSE) traceback (cm_dpalign.c:374-411).  NOTE: unlike
 * ckpt_optacc_traceback, there is NO d==0 BEGL_S/BEGR_S "allow_S_local_end"
 * special case -- that case is gated on do_optacc in cm_alignT_hb (:300), so
 * the CYK oracle never uses it and neither must this engine. */
static int
ckpt_cyk_traceback(CM_t *cm, char *errbuf, int L, float ***CYstore,
                    int *jmin, int *jmax, CP9Bands_t *cp9b,
                    ckpt_ysh_fetch_fn fetch, void *fctx, int begin_b, Parsetree_t **ret_tr)
{
  int status;
  Parsetree_t *tr  = NULL;
  ESL_STACK   *pda = NULL;
  int v = 0, i = 1, j = L, d = L, k, y, yoffset, bifparent;
  int jp_v, dp_v;

  tr = CreateParsetree(100);
  if (tr == NULL) { status = eslEMEM; goto ERROR; }
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0);
  pda = esl_stack_ICreate();
  if (pda == NULL) { status = eslEMEM; goto ERROR; }

  while (1) {
    if (cm->sttype[v] != EL_st) {
      jp_v = j - jmin[v];
      dp_v = d - hd_min(cp9b, v, jp_v);
    }

    if (cm->sttype[v] == B_st) {
      k = ckpt_cyk_bsearch(cm, jmin, jmax, cp9b, CYstore, NULL, v, j, d, NULL);
      if (k < 0) ESL_XFAIL(eslEINCOMPAT, errbuf, "ckpt_cyk_traceback: B state v=%d: no valid k (band inconsistency)", v);
      if ((status = esl_stack_IPush(pda, j))       != eslOK) goto ERROR;
      if ((status = esl_stack_IPush(pda, k))       != eslOK) goto ERROR;
      if ((status = esl_stack_IPush(pda, tr->n-1)) != eslOK) goto ERROR;
      j = j - k;
      d = d - k;
      i = j - d + 1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    }
    else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;  /* traceback complete */
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];
      y = cm->cnum[v];                       /* right START state */
      i = j - d + 1;
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
    }
    else {
      yoffset = fetch(fctx, v, jp_v, dp_v);
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default: ESL_XFAIL(eslEINVAL, errbuf, "ckpt_cyk_traceback: bogus state type v=%d", v);
      }
      d = j - i + 1;
      if (yoffset == (char) USED_EL) {              /* R1-L (084): local end */
        InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
        v = cm->M;
      }
      else if (yoffset == (char) USED_LOCAL_BEGIN) { /* R1-L (084): local begin, once, from ROOT_S */
        InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, begin_b);
        v = begin_b;
      }
      else {
        y = cm->cfirst[v] + (int) yoffset;
        InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
        v = y;
      }
    }
  }
  esl_stack_Destroy(pda);
  *ret_tr = tr;
  return eslOK;

 ERROR:
  if (pda) esl_stack_Destroy(pda);
  if (tr)  FreeParsetree(tr);
  return status;
}

/* Function: cm_CheckptCYKAlignHB()
 * Incept:   Brief 26_0610-078 R1
 *
 * Purpose:  sqrt(M)-memory checkpointed CYK max-DP alignment, discovering
 *           its OWN bifurcation k* pins (no externally-supplied kpin[], unlike
 *           the rung-3 OA engines above) via ckpt_cyk_bsearch().  Non-truncated;
 *           GLOBAL or LOCAL (EL + local begin).  Brief 26_0610-084 (R1-L) added
 *           local support (transcribed from cm_CheckptOptAccAlignHB's own
 *           have_el/have_local_begin plumbing + cm_CYKInsideAlignHB's CYK-
 *           specific EL/local-begin semantics), so cm_CheckptCYKAlignHB now
 *           serves as do_checkpt_r3's pass-1 in both configs.
 *
 *           STEP CYK: checkpointed CYK max-DP (descending sweep, identical
 *                     retention scheme to cm_CheckptOptAccAlignHB's STEP OA)
 *                     -> CYstore = {all chain roots} + {sqrt(M) seeds}.
 *           STEP TB : traceback (ckpt_cyk_traceback) with yshadow
 *                     block-recomputed on demand (ckpt_cyk_ysh_fetch) and
 *                     B-state k found by re-running ckpt_cyk_bsearch against
 *                     the still-resident CYstore chain-root decks.
 *
 *           Byte-exact vs CYKDivideAndConquerHB expected (both are exact CYK
 *           argmax under the same HMM bands) modulo float-sum-order rounding
 *           between D&C recombination and this engine's single-sweep fill --
 *           this project's own established tolerance for that comparison
 *           (hbdnc_drv.c) applies, not literal bitwise equality.
 */
int
cm_CheckptCYKAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,
                     Parsetree_t **ret_tr, float *ret_sc)
{
  int      status;
  CKPT_CTX cx;
  int      M = cm->M;
  int      v;
  float  ***CYstore = NULL, ***tba = NULL;
  char   ***tysh = NULL;
  Parsetree_t *tr = NULL;
  ckpt_cyk_ysh_ctx fctx;

  /* Brief 26_0610-084 (R1-L): the brief-078 global-only fail-fast gate is
   * removed -- local (EL + local begin) is now validated byte-exact vs
   * cm_AlignHB(do_optacc=FALSE) across 5S/RNaseP/hairpin(bifs=0) in both
   * configs, incl. genuine internal local-begin entries on fragment seqs. */

  CP9Bands_t *cp9b = cm->cp9b;
  memset(&cx, 0, sizeof(cx));
  cx.cm = cm; cx.dsq = dsq; cx.L = L; cx.M = M;
  cx.jmin = cm->cp9b->jmin; cx.jmax = cm->cp9b->jmax;
  cx.imin = cm->cp9b->imin; cx.imax = cm->cp9b->imax;
  cx.hdmin = cm->cp9b->hdmin; cx.hdmax = cm->cp9b->hdmax;
  cx.cur_bytes = cx.peak_bytes = 0;
  cx.kpin = NULL; cx.ifull = NULL;
  cx.my_lpp = NULL; cx.my_rpp = NULL;
  cx.deck_nc = NULL; cx.deck_njr = NULL;
  /* R1-L (084): CYK local support.  Unlike the OA half, CYK EL is the closed-form
   * ramp el_selfsc*(d-sd)+endsc[v] (no elalpha/eldmax/el_esc decks needed -- see
   * ckpt_cyk_deck), so only have_el/el_selfsc + have_local_begin are consumed. */
  cx.have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE; cx.el_selfsc = cm->el_selfsc;
  cx.eldmax = NULL; cx.elbeta = NULL; cx.elalpha = NULL;
  cx.el_esc = NULL; cx.el_endsc = IMPOSSIBLE;
  cx.have_local_begin = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;
  cx.bsc_fwd = IMPOSSIBLE; cx.fwd_begin_applied = FALSE;
  cx.begin_bsc = IMPOSSIBLE; cx.begin_b = -1;

  if (cx.jmin[0] > L || cx.jmax[0] < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptCYKAlignHB(): L outside ROOT_S j band");
  int jp_0 = L - cx.jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_CheckptCYKAlignHB(): L outside ROOT_S d band");
  int Lp_0 = L - hd_min(cp9b, 0, jp_0);

  ESL_ALLOC(cx.deck_nc,  sizeof(int64_t) * M);
  ESL_ALLOC(cx.deck_njr, sizeof(int)     * M);
  int Delta = 0;
  for (v = 0; v < M; v++) {
    int njr = cx.jmax[v] - cx.jmin[v] + 1; if (njr < 0) njr = 0;
    cx.deck_njr[v] = njr;
    int64_t nc = 0; int jp;
    for (jp = 0; jp < njr; jp++) { int w = hd_max(cp9b, v, jp)-hd_min(cp9b, v, jp)+1; if (w>0) nc += w; }
    cx.deck_nc[v] = nc;
    /* max forward reach (child index - v) over non-E/non-B states, exactly as
     * cm_CheckptOptAccAlignHB computes it (cm_dpalign.c ~3040-3045) -- B_st is
     * excluded because bifurcation children (BEGL_S/BEGR_S) are chain roots,
     * retained permanently regardless of index distance, not windowed by Delta. */
    if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) {
      int ymax = cm->cfirst[v] + cm->cnum[v] - 1;
      if (ymax - v > Delta) Delta = ymax - v;
    }
  }
  int B = (int) (sqrt((double)M) + 0.5); if (B < 1) B = 1;

  /* ============================================================= */
  /* STEP CYK: checkpointed CYK max-DP -> roots + sqrt(M) CYK seeds */
  /* ============================================================= */
  ESL_ALLOC(CYstore, sizeof(float**) * M);
  for (v = 0; v < M; v++) CYstore[v] = NULL;
  for (v = M-1; v >= 0; v--) {
    char **ysh = (cm->sttype[v] == B_st) ? NULL : ckpt_cdeck_alloc(&cx, v);
    CYstore[v] = ckpt_deck_alloc(&cx, v);
    ckpt_cyk_deck(&cx, v, CYstore, NULL, ysh);
    if (ysh) ckpt_cdeck_free(&cx, v, ysh);
    int y = v + Delta;
    if (y < M && CYstore[y] != NULL && ! ckpt_is_chain_root(cm, y)) {
      if ((y % B) < Delta) { /* retain checkpoint seed */ }
      else { ckpt_deck_free(&cx, y, CYstore[y]); CYstore[y] = NULL; }
    }
  }
  float sc = CYstore[0][jp_0][Lp_0];

  /* ============================================================= */
  /* STEP TB: traceback (block-recompute yshadow, live B-state k-search) */
  /* ============================================================= */
  ESL_ALLOC(tba,  sizeof(float**) * M);
  ESL_ALLOC(tysh, sizeof(char**)  * M);
  for (v = 0; v < M; v++) { tba[v] = NULL; tysh[v] = NULL; }
  fctx.cx = &cx; fctx.M = M; fctx.B = B; fctx.CYstore = CYstore;
  fctx.tba = tba; fctx.tysh = tysh;
  fctx.cur_blk = -1; fctx.blk_lo = 0; fctx.blk_hi = -1;

  if ((status = ckpt_cyk_traceback(cm, errbuf, L, CYstore, cx.jmin, cx.jmax, cp9b,
                                   ckpt_cyk_ysh_fetch, &fctx, cx.begin_b, &tr)) != eslOK) goto ERROR;

  { int w; for (w = fctx.blk_lo; w <= fctx.blk_hi; w++) {
      if (tba[w])  { ckpt_deck_free (&cx, w, tba[w]);  tba[w]  = NULL; }
      if (tysh[w]) { ckpt_cdeck_free(&cx, w, tysh[w]); tysh[w] = NULL; } } }

  if (getenv("INFERNAL_CKPT_VERBOSE")) {
    int64_t full_cube_cells = 0; for (v = 0; v < M; v++) full_cube_cells += cx.deck_nc[v];
    double full_mb = full_cube_cells * 4 / (1024.0*1024.0);
    double peak_mb = cx.peak_bytes / (1024.0*1024.0);
    fprintf(stderr, "# cm_CheckptCYKAlignHB (R1): M=%d L=%d B=%d sc=%.5f  CYK-DP peak=%.2f Mb  full-CYK-cube(1x)=%.2f Mb  win~%.1fx\n",
            M, L, B, sc, peak_mb, full_mb, (peak_mb>0.) ? full_mb/peak_mb : 0.);
  }

  for (v = 0; v < M; v++) if (CYstore[v]) ckpt_deck_free(&cx, v, CYstore[v]);
  free(CYstore); free(tba); free(tysh);
  free(cx.deck_nc); free(cx.deck_njr);

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  if (ret_sc != NULL) *ret_sc = sc;
  return eslOK;

 ERROR:
  if (CYstore) { for (v = 0; v < M; v++) if (CYstore[v]) ckpt_deck_free(&cx, v, CYstore[v]); free(CYstore); }
  if (tba)  { for (v = 0; v < M; v++) if (tba[v])  ckpt_deck_free(&cx, v, tba[v]);  free(tba); }
  if (tysh) { for (v = 0; v < M; v++) if (tysh[v]) ckpt_cdeck_free(&cx, v, tysh[v]); free(tysh); }
  if (cx.deck_nc)  free(cx.deck_nc);
  if (cx.deck_njr) free(cx.deck_njr);
  if (tr) FreeParsetree(tr);
  return status;
}

/* Function: cm_CYKInsideAlign()
 * Date:     EPN, Sun Nov 18 19:37:39 2007
 *           
 * Purpose:  Run the inside phase of a CYK alignment. Non-banded
 *           version. See cm_CYKInsideAlignHB() for HMM banded version.
 *         
 *           This function must perform a complete alignment, aligning
 *           the full sequence 1..L to the ROOT_S state 0 of the model.
 *
 *           We deal with local begins by keeping track of the optimal
 *           state that we could enter and account for the whole target 
 *           sequence: b = argmax_v  alpha_v(1,L) + log t_0(v),
 *           and bsc is the score for that. 
 *
 *           If local begins are on (cm->flags & CMH_LOCAL_BEGIN), the
 *           optimal alignment must use a local begin transition,
 *           0->b, and we have to be able to trace that back. If local
 *           begins are on, we return a valid b (the optimal 0->b
 *           choice), yshad[0][L][L] will be USE_LOCAL_BEGIN, telling
 *           cm_alignT() to check b and start with a local 0->b entry
 *           transition. 
 *
 *           Note on history of this function: It was previously
 *           fast_cyk_align() (up to Infernal 1.0.2), which was
 *           based on inside() from cm_dpsmall.c.
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence [1..L]   
 *           L         - length of the dsq to align
 *           size_limit- max size in Mb for DP matrix
 *           mx        - the DP matrix to fill in
 *           shmx      - the shadow matrix to fill in
 *           ret_b     - RETURN: local begin state if local begins are on
 *           ret_sc    - RETURN: score of optimal, CYK parsetree 
 *                       
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslERANGE> if required mx or shmx size exceeds <size_limit>
 *           In this case alignment has been aborted, <ret_*> variables are not valid
 */
int
cm_CYKInsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, 
		  CM_MX *mx, CM_SHADOW_MX *shmx, int *ret_b, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;       /* j - sdr */
  int      d_sd;        /* d - sd */
  float    tsc;         /* a transition score */
  int64_t  c;           /* 64-bit int counter */

  /* the DP matrix */
  float ***alpha   = mx->dp;        /* pointer to the alpha DP matrix */
  char  ***yshadow = shmx->yshadow; /* pointer to the yshadow matrix */
  int   ***kshadow = shmx->kshadow; /* pointer to the kshadow matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;

  /* grow the matrices based on the current sequence */
  if((status = cm_mx_GrowTo       (cm,   mx, errbuf, L, size_limit)) != eslOK) return status;
  if((status = cm_shadow_mx_GrowTo(cm, shmx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE, all cells of shadow matrix to USED_EL */
  esl_vec_FSet(mx->dp_mem, mx->ncells_valid, IMPOSSIBLE);
  for(c = 0; c < shmx->y_ncells_valid; c++) shmx->yshadow_mem[c] = USED_EL;
  esl_vec_ISet(shmx->kshadow_mem, shmx->k_ncells_valid, USED_EL);

  /* precalcuate all possible local end scores, for local end emits of 1..L residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* if local ends are on, replace the EL deck IMPOSSIBLEs with EL scores */
  if(cm->flags & CMH_LOCAL_END) { 
    for (j = 0; j <= L; j++) {
      for (d = 0;  d <= j; d++) alpha[cm->M][j][d] = el_scA[d];
    }
  }

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);

    /* re-initialize the J deck if we can do a local end from v */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = 0; j <= L; j++) { 
	for (d = sd; d <= j; d++) { 
	  alpha[v][j][d] = el_scA[d-sd] + cm->endsc[v];
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */
    
    if(cm->sttype[v] == E_st) { 
      for (j = 0; j <= L; j++) {
	alpha[v][j][0] = 0.;
	/* rest of deck remains IMPOSSIBLE */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][j][d] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (j = sdr; j <= L; j++) {
	j_sdr = j - sdr;
	for (d = sd; d <= j; d++) {
	  d_sd = d - sd;
	  i    = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = alpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc; 
	      yshadow[v][j][d]    = yoffset;
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
      for (j = sdr; j <= L; j++) {
	j_sdr = j - sdr;
	for (d = sd; d <= j; d++) {
	  d_sd = d - sd;
	  i = j - d + 1;
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = alpha[y][j_sdr][d_sd] + tsc_v[yoffset]) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc; 
	      yshadow[v][j][d]    = yoffset;
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

	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;

	  for (d = sd; d <= j; d++) {
	    if((sc = alpha[y][j_sdr][d - sd] + tsc) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc;
	      yshadow[v][j][d]    = yoffset;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = 0; j <= L; j++) {
	  i = j - 1;
	  for (d = sd; d <= j; d++) 
	    alpha[v][j][d] += esc_v[dsq[j-d+1]];
	}
	break;
      case MR_st:
	for (j = 0; j <= L; j++) {
	  for (d = sd; d <= j; d++)
	    alpha[v][j][d] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = 0; j <= L; j++) {
	  i = j - 1;
	  for (d = sd; d <= j; d++)
	    alpha[v][j][d] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = 0; j <= L; j++) {
	for (d = 0; d <= j; d++)
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  for (k = 0; k <= d; k++) {
	    if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) { 
	      alpha[v][j][d] = sc;
	      kshadow[v][j][d] = k;
	    }
	  }
	}
      }
    }
      
    /* allow local begins, if nec */
    if ((cm->flags & CMH_LOCAL_BEGIN) && 
	(NOT_IMPOSSIBLE(cm->beginsc[v])) && 
	(alpha[v][L][L] + cm->beginsc[v] > bsc)) {
      b   = v;
      bsc = alpha[v][L][L] + cm->beginsc[v];
    }
  } /* finished calculating deck v. */
  
  /* Check for whether we need to store an optimal local begin score
   * as the optimal overall score, and if we need to put a flag
   * in the shadow matrix telling cm_alignT() to use the b we return.
   */
  if (bsc > alpha[0][L][L]) {
    alpha[0][L][L] = bsc;
    yshadow[0][L][L] = USED_LOCAL_BEGIN;
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_cykmx",   "w"); cm_mx_Dump(fp1, mx); fclose(fp1); */ 
  /* FILE *fp2; fp2 = fopen("tmp.std_cykshmx", "w"); cm_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2); */
#endif
  
  sc = alpha[0][L][L];

  free(el_scA);

  if (ret_b   != NULL) *ret_b  = b;    /* b is -1 if local begins are off */
  if (ret_sc  != NULL) *ret_sc = sc;

  ESL_DPRINTF1(("#DEBUG: cm_CYKInsideAlign return sc: %f\n", sc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_CYKInsideAlignHB()
 * Date:     EPN 03.29.06 [EPN started] 
 *           SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phase of a CYK alignment using bands 
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 *           Otherwise, (meant to be) identical to cm_CYKInsideAlign()
 *           see that function for more information.
 *
 *           A CM_HB_MX DP matrix must be passed in. Only cells valid
 *           within the bands given in the CP9Bands_t <cm->cp9b> will
 *           be valid.
 *
 *           Note on history of this function: It was previously
 *           fast_cyk_align_hb() (up to Infernal 1.0.2), which was
 *           based on inside_qdb_me() which was based on inside().
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence [1..L]   
 *           L         - length of the dsq to align
 *           size_limit- max size in Mb for DP matrix
 *           mx        - the DP matrix to fill in, only cells within bands are valid
 *           shmx      - the shadow matrix to fill in, only cells within bands are valid
 *           ret_b     - RETURN: best local begin state, or NULL if unwanted
 *           ret_sc    - RETURN: score of optimal, CYK parsetree 
 *                       
 * Returns: <eslOK> on success.
 * 
 * Throws:  <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *          <eslEINVAL> if the full sequence is not within the bands for state 0
 *          In either case alignment has been aborted, ret_* variables are not valid
 * 
 */
int
cm_CYKInsideAlignHB(CM_t *cm, char *errbuf,  ESL_DSQ *dsq, int L, float size_limit, 
		    CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, int *ret_b, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..L-1] probability of local end emissions of length d */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;       /* j - sdr */
  int      c;           /* 64-bit int counter */

  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y;         /* d index for state v/y in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Lp;                 /* L index also changes depending on state */
  float    tsc;                /* a transition score */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  int      jp_0;               /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;               /* L offset in ROOT_S's (v==0) d band */

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b    = cm->cp9b;
  int        *jmin    = cp9b->jmin;  
  int        *jmax    = cp9b->jmax;
  int       **hdmin   = cp9b->hdmin;
  int       **hdmax   = cp9b->hdmax;
  float    ***alpha   = mx->dp;        /* pointer to the alpha DP matrix */
  char     ***yshadow = shmx->yshadow; /* pointer to the yshadow matrix */
  int      ***kshadow = shmx->kshadow; /* pointer to the kshadow matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  /* ensure a full alignment to ROOT_S (v==0) is allowed by the bands */
  if (cp9b->jmin[0] > L || cp9b->jmax[0] < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, cp9b->jmin[0], cp9b->jmax[0]);
  jp_0 = L - jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) 
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, hd_min(cp9b, 0, jp_0), hd_max(cp9b, 0, jp_0));
  Lp_0 = L - hd_min(cp9b, 0, jp_0);

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo       (cm,   mx, errbuf, cp9b, L, size_limit)) != eslOK) return status;
  if((status = cm_hb_shadow_mx_GrowTo(cm, shmx, errbuf, cp9b, L, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..L residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) 
   */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE);
  if(shmx->y_ncells_valid > 0) for(c = 0; c < shmx->y_ncells_valid; c++) shmx->yshadow_mem[c] = USED_EL;
  /* for B states, shadow matrix holds k, length of right fragment, this will be overwritten */
  if(shmx->k_ncells_valid > 0) esl_vec_ISet(shmx->kshadow_mem, shmx->k_ncells_valid, 0);

  /* EL deck optimization: we skip filling alpha[cm->M] here because
   * alpha[cm->M][j][d] == el_scA[d] always and we substitute el_scA[d]
   * directly at the read sites in the main recursion below.
   */

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];

    /* re-initialize if we can do a local end from v */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	if(hd_min(cp9b, v, jp_v) >= sd) { 
	  d    = hd_min(cp9b, v, jp_v);
	  dp_v = 0;
	}
	else { 
	  d    = sd;
	  dp_v = sd - hd_min(cp9b, v, jp_v);
	}
	for (; d <= hd_max(cp9b, v, jp_v); dp_v++, d++) {
	  if(d >= sd) {
	    alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
	  }
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hd_min(cp9b, v, jp_v) == 0));
	ESL_DASSERT1((hd_max(cp9b, v, jp_v) == 0));
	alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr valid for state y? */
	
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) { /* for each valid d for v, j */
	  i = j - d + 1;
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > alpha[v][jp_v][dp_v])
		{
		  alpha[v][jp_v][dp_v] = sc; 
		  yshadow[v][jp_v][dp_v]    = yoffset;
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
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) { /* for each valid d for v, j */
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > alpha[v][jp_v][dp_v])
		{
		  alpha[v][jp_v][dp_v] = sc; 
		  yshadow[v][jp_v][dp_v]    = yoffset;
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
	  ESL_DASSERT1((jp_v     >= 0 && jp_v     <= (jmax[v]-jmin[v])));
	  ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	  
	/* d must satisfy:
	 * d >= hdmin[v][jp_v]
	 * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	 * d <= hdmax[v][jp_v]
	 * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	 * this reduces to two ESL_MAX calls
	 */
	  dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y_sdr) + sd);
	  dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y_sdr) + sd);
	  dpn     = dn - hd_min(cp9b, v, jp_v);
	  dpx     = dx - hd_min(cp9b, v, jp_v);
	  dp_y_sd = dn - hd_min(cp9b, y, jp_y_sdr) - sd;
	  	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
	    if((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc) > alpha[v][jp_v][dp_v]) {
	      alpha[v][jp_v][dp_v] = sc;
	      yshadow[v][jp_v][dp_v]    = yoffset;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hd_min(cp9b, v, jp_v) + 1;
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hd_min(cp9b, v, jp_v) + 1;
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
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
	kn = ((j-jmax[y]) > (hd_min(cp9b, z, jp_z))) ? (j-jmax[y]) : hd_min(cp9b, z, jp_z);
        kn = ESL_MAX(kn, 0); /* kn must be non-negative, added with fix to bug i36 */
        /* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hd_max(cp9b, z, jp_z))) ?  jp_y       : hd_max(cp9b, z, jp_z);
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
	      
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
	    if((k >= d - hd_max(cp9b, y, jp_y-k)) && k <= d - hd_min(cp9b, y, jp_y-k)) {
	      /* for current k, all 6 inequalities have been satisified 
	       * so we know the cells corresponding to the platonic 
	       * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
	       * alpha[z][j][k] are all within the bands. These
	       * cells correspond to alpha[v][jp_v][dp_v], 
	       * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
	       * and alpha[z][jp_z][k-hdmin[jp_z]];
	       */
	      kp_z = k-hd_min(cp9b, z, jp_z);
	      dp_y = d-hd_min(cp9b, y, jp_y-k);

	      if ((sc = alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]) 
		  > alpha[v][jp_v][dp_v]) { 
		alpha[v][jp_v][dp_v] = sc;
		kshadow[v][jp_v][dp_v] = k;
	      }
	    }
	  }
	}
      }
    } /* finished calculating deck v. */
         
    /* allow local begins, if nec */
    if(cm->flags & CMH_LOCAL_BEGIN) { 
      if(L >= jmin[v] && L <= jmax[v]) { 
	jp_v = L - jmin[v];
	Lp   = L - hd_min(cp9b, v, jp_v);
	if(L >= hd_min(cp9b, v, jp_v) && L <= hd_max(cp9b, v, jp_v)) { 
	/* If we get here alpha[v][jp_v][Lp] is a valid cell
	 * in the banded alpha matrix, corresponding to 
	 * alpha[v][L][L] in the platonic matrix.
	 */
	/* Check for local begin getting us to the root.
	 * This is "off-shadow": if/when we trace back, we'll handle this
	 * case separately (and we'll know to do it because we'll immediately
	 * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	 * to jump right to state b; see below)
	 */
	  if (NOT_IMPOSSIBLE(cm->beginsc[v]) && 
	      (alpha[v][jp_v][Lp] + cm->beginsc[v] > bsc)) {
	    b   = v;
	    bsc = alpha[v][jp_v][Lp] + cm->beginsc[v];
	  }
	}
      }
    }
  } /* end loop over all v */
  /* Check for whether we need to store an optimal local begin score
   * as the optimal overall score, and if we need to put a flag
   * in the shadow matrix telling cm_alignT() to use the b we return.
   */
  if (NOT_IMPOSSIBLE(bsc) && (bsc > alpha[0][jp_0][Lp_0])) {
    alpha[0][jp_0][Lp_0] = bsc;
    yshadow[0][jp_0][Lp_0] = USED_LOCAL_BEGIN;
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_cykhbmx", "w");   cm_hb_mx_Dump(fp1, mx); fclose(fp1); */
  /* FILE *fp2; fp2 = fopen("tmp.std_cykhbshmx", "w"); cm_hb_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2); */
#endif
  
  sc = alpha[0][jp_0][Lp_0];

  free(el_scA);
  free(yvalidA);

  if (ret_b != NULL)  *ret_b   = b;    /* b is -1 if local begins are off */
  if (ret_sc != NULL) *ret_sc = sc;

  ESL_DPRINTF1(("#DEBUG: cm_CYKInsideAlignHB return sc: %f\n", sc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_InsideAlign()
 * Date:     EPN, Mon Nov 19 06:21:51 2007
 *
 * Purpose:  Run the inside algorithm on a target sequence 
 *           without using bands. The full target sequence
 *           1..L is aligned (only full alignments will
 *           contribute to the Inside score).
 *
 *           Identical to cm_InsideAlignHB() but no bands
 *           are used.
 * 
 *           Very similar to cm_CYKInsideAlign(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - we do Inside, not CYK
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *
 *           This function complements cm_OutsideAlign().
 *
 *           Note: renamed from FastInsideAlign() [EPN, Wed Sep 14 06:13:37 2011].
 *
 * Args:     cm         - the model
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitized sequence
 *           L          - target sequence length
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx         - the dp matrix, grown and filled here
 *           ret_sc     - RETURN: log P(S|M)/P(S|R), as a bit score
 *                       
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslERANGE> if required CM_MX size exceeds <size_limit>
 *           In this case alignment has been aborted, ret_sc is not valid
 */
int
cm_InsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, CM_MX *mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* the final score */
  float    tsc;         /* a temporary variable holding a transition score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  float    bsc;		/* summed score for using all local begins */
  float   *el_scA;      /* [0..d..L-1] probability of local end emissions of length d */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdl;         /* StateLeftDelta(cm->sttype[v] */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;       /* j - sdr */
  int      d_sd;        /* d - sd */

  /* the DP matrix */
  float ***alpha = mx->dp;     /* pointer to the alpha DP matrix */

  /* Allocations and initializations */
  bsc = IMPOSSIBLE;
  
  /* grow the matrix based on the current sequence */
  if((status = cm_mx_GrowTo(cm, mx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* precalcuate all possible local end scores, for local end emits of 1..L residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* if local ends are on, replace the EL deck IMPOSSIBLEs with EL scores */
  if(cm->flags & CMH_LOCAL_END) { 
    for (j = 0; j <= L; j++) {
      for (d = 0;  d <= j; d++) alpha[cm->M][j][d] = el_scA[d];
    }
  }
  
  /* Main recursion  */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; 
    float const *tsc_v = cm->tsc[v];
    sd   = StateDelta(cm->sttype[v]);
    sdl  = StateLeftDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    
    /* re-initialize the J deck if we can do a local end from v */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = 0; j <= L; j++) {
	for (d = sd; d <= j; d++) alpha[v][j][d] = el_scA[d-sd] + cm->endsc[v];
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    /* E_st: easy, no children, and d must be 0 for all valid j */
    if(cm->sttype[v] == E_st) { 
      for (j = 0; j <= L; j++) {
	alpha[v][j][0] = 0.;
	/* rest of deck remains IMPOSSIBLE */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (j = sdr; j <= L; j++) {
	j_sdr = j - sdr;
	for (d = sd; d <= j; d++) {
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
      for (j = sdr; j <= L; j++) {
	j_sdr = j - sdr;
	for (d = sd; d <= j; d++) {
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

	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++) {
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], (alpha[y][j_sdr][d-sd] + tsc));;
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = 0; j <= L; j++) {
	  i = j - sdl;
	  for (d = sd; d <= j; d++) 
	    alpha[v][j][d] += esc_v[dsq[j-d+1]];
	}
	break;
      case MR_st:
	for (j = 0; j <= L; j++) {
	  for (d = sd; d <= j; d++)
	    alpha[v][j][d] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = 0; j <= L; j++) {
	  i = j - sdl;
	  for (d = sd; d <= j; d++)
	    alpha[v][j][d] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = 0; j <= L; j++) {
	for (d = 0; d <= j; d++)
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  for (k = 0; k <= d; k++) {
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], alpha[y][j-k][d-k] + alpha[z][j][k]); 
	  }
	}
      }
    }

    /* allow local begins, if nec */
    if ((cm->flags & CMH_LOCAL_BEGIN) && 
	(NOT_IMPOSSIBLE(cm->beginsc[v]))) { 
      /* add in score for local begin getting us to the root. */
      bsc = FLogsum(bsc, alpha[v][L][L] + cm->beginsc[v]);
    }
  } /* finished calculating deck v. */

  /* include the bsc as part of alpha[0][L][L] */
  alpha[0][L][L] = FLogsum(alpha[0][L][L], bsc);

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_imx", "w");   cm_mx_Dump(fp1, mx); fclose(fp1); */
#endif

  sc =  alpha[0][L][L];

  free(el_scA);

  if(ret_sc != NULL) *ret_sc = sc;

  ESL_DPRINTF1(("#DEBUG: cm_InsideAlign() return sc: %f\n", sc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_InsideAlignHB()
 * Date:     EPN, Thu Nov  8 18:24:41 2007
 *
 * Purpose:  Run the inside algorithm on a target sequence using bands 
 *           in the j and d dimensions of the DP matrix. Bands
 *           were obtained from an HMM Forward-Backward parse
 *           of the target sequence. Uses float log odds scores.
 *           The full target sequence 1..L is aligned (only full
 *           alignments will contribute to the Inside score).
 * 
 *           Very similar to cm_CYKInsideAlignHB(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - we do Inside, not CYK
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *
 *           This function complements cm_OutsideAlignHB().
 *
 *           Note: renamed from FastInsideAlignHB() [EPN, Wed Sep 14 06:13:08 2011].
 *
 * Args:     cm         - the model
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitized sequence
 *           L          - target sequence length
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx         - the dp matrix, only cells within bands in cp9b will be valid
 *           ret_sc     - RETURN: log P(S|M)/P(S|R) (given bands), as a bit score
 * 
 * Returns:  <eslOK> on success.
 *
 * Throws:  <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *          <eslEINVAL> if the full sequence is not within the bands for state 0
 *          In either case alignment has been aborted, ret_sc is not valid
 */
int
cm_InsideAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, CM_HB_MX *mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* the final score */
  float    tsc;         /* a temporary variable holding a transition score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  float    bsc;		/* summed score for using all local begins */
  float   *el_scA;      /* [0..d..L-1] probability of local end emissions of length d */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;       /* j - sdr */

  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int     *yvalidA;            /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y;         /* d index for state v/y in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Lp;                 /* L also changes depending on state */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  int      jp_0;               /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;               /* L offset in ROOT_S's (v==0) d band */

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
  /* ensure a full alignment to ROOT_S (v==0) is allowed by the bands */
  if (cp9b->jmin[0] > L || cp9b->jmax[0] < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_InsideAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, cp9b->jmin[0], cp9b->jmax[0]);
  jp_0 = L - jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) 
    ESL_FAIL(eslEINVAL, errbuf, "cm_InsideAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, hd_min(cp9b, 0, jp_0), hd_max(cp9b, 0, jp_0));
  Lp_0 = L - hd_min(cp9b, 0, jp_0);

  /* grow the matrix based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo(cm, mx, errbuf, cp9b, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands)
   */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);


  /* EL deck optimization: skip filling alpha[cm->M]. The main recursion
   * uses el_scA[d-sd] directly at read sites (see ~line 1873). The
   * PosteriorHB function uses el_scA[d] directly instead of alpha[cm->M][j][d].
   */

  /* Main recursion  */
  for (v = cm->M-1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v]; 
    float const *tsc_v = cm->tsc[v];
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];
    
    /* re-initialize the J deck if we can do a local end from v */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v  = j - jmin[v];
	for (dp_v = 0, d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); dp_v++, d++) 
	  alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    /* E_st: easy, no children, and d must be 0 for all valid j */
    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hd_min(cp9b, v, jp_v) == 0));
	ESL_DASSERT1((hd_max(cp9b, v, jp_v) == 0));
	alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	/* rest of deck remains IMPOSSIBLE */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) { /* for each valid d for v, j */
	  i = j - d + 1;
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
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
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) { /* for each valid d for v, j */
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
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
	  
	  dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y_sdr) + sd);
	  dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y_sdr) + sd);
	  dpn     = dn - hd_min(cp9b, v, jp_v);
	  dpx     = dx - hd_min(cp9b, v, jp_v);
	  dp_y_sd = dn - hd_min(cp9b, y, jp_y_sdr) - sd;
	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y_sdr][dp_y_sd] + tsc));;
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hd_min(cp9b, v, jp_v) + 1;
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hd_min(cp9b, v, jp_v) + 1;
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	  }
      default: /* no emission */
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
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
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	kn = ((j-jmax[y]) > (hd_min(cp9b, z, jp_z))) ? (j-jmax[y]) : hd_min(cp9b, z, jp_z);
        kn = ESL_MAX(kn, 0); /* kn must be non-negative, added with fix to bug i36 */
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hd_max(cp9b, z, jp_z))) ?  jp_y       : hd_max(cp9b, z, jp_z);
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
	      
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
	    if((k >= d - hd_max(cp9b, y, jp_y-k)) && k <= d - hd_min(cp9b, y, jp_y-k)) {
	      /* for current k, all 6 inequalities have been satisified 
	       * so we know the cells corresponding to the platonic 
	       * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
	       * alpha[z][j][k] are all within the bands. These
	       * cells correspond to alpha[v][jp_v][dp_v], 
	       * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
	       * and alpha[z][jp_z][k-hdmin[jp_z]];
	       */
	      kp_z = k-hd_min(cp9b, z, jp_z);
	      dp_y = d-hd_min(cp9b, y, jp_y-k);

	      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]); 
	    }
	  }
	}
      }
    }
      
    /* allow local begins, if nec */
    if((cm->flags & CMH_LOCAL_BEGIN) && 
       NOT_IMPOSSIBLE(cm->beginsc[v])) { 
      if(L >= jmin[v] && L <= jmax[v]) { 
	jp_v = L - jmin[v];
	Lp   = L - hd_min(cp9b, v, jp_v);
	if(L >= hd_min(cp9b, v, jp_v) && L <= hd_max(cp9b, v, jp_v)) { 
	  /* If we get here alpha[v][jp_v][Lp] is a valid cell
	   * in the banded alpha matrix, corresponding to 
	   * alpha[v][L][L] in the platonic matrix.
	   */
	  /* Check for local begin getting us to the root.
	   */
	  bsc = FLogsum(bsc, (alpha[v][jp_v][Lp] + cm->beginsc[v]));
	}
      }
    }
  } /* end loop over all v */

  /* include the bsc as part of alpha[0][jp_0][Lp_0] */
  if (NOT_IMPOSSIBLE(bsc)) { 
    alpha[0][jp_0][Lp_0] = FLogsum(alpha[0][jp_0][Lp_0], bsc);
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp; fp = fopen("tmp.std_ihbmx", "w"); cm_hb_mx_Dump(fp, mx); fclose(fp); */
#endif

  sc = alpha[0][jp_0][Lp_0];

  free(el_scA);
  free(yvalidA);

  if(ret_sc != NULL) *ret_sc = sc;

  ESL_DPRINTF1(("#DEBUG: cm_InsideAlignHB() return sc: %f\n", sc));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_OptAccAlign()
 * Date:     EPN, Sun Nov 18 20:45:22 2007
 *           EPN, Sat Oct  1 05:57:49 2011 (updated to use emit matrices)
 *   
 * Purpose:  Run the Holmes/Durbin optimal accuracy algorithm 
 *           on a full target sequence 1..L, given a pre-filled
 *           posterior matrix. Uses float log odds scores.
 *           Non-banded version. See cm_OptAccAlignHB() for 
 *           HMM banded version.
 * 
 *           A CM_EMIT_MX matrix <emit_mx> must be passed in, filled by
 *           cm_EmitterPosterior(), with values:
 *
 *           l_pp[v][i]: log of the posterior probability that state v
 *           emitted residue i leftwise either at (if a match state)
 *           or *after* (if an insert state) the left consensus
 *           position modeled by state v's node.
 *
 *           r_pp[v][i]: log of the posterior probability that state v
 *           emitted residue i rightwise either at (if a match
 *           state) or *before* (if an insert state) the right
 *           consensus position modeled by state v's node.
 *
 *           l_pp[v] is NULL for states that do not emit leftwise
 *           r_pp[v] is NULL for states that do not emit rightwise
 *
 *           Additionally, a CM_MX DP matrix <mx> and CM_SHADOW_MX
 *           <shmx> must be passed in. <shmx> will be expanded and
 *           filled here with traceback pointers to allow the
 *           optimally accurate parsetree to be recovered in
 *           cm_alignT() and <mx> will be expanded and filled with the
 *           optimal accuracy scores, where:
 *
 *           mx->dp[v][j][d]: log of the sum of the posterior
 *           probabilities of emitting residues i..j in the subtree
 *           rooted at v.
 *
 *           The optimally accurate parsetree, i.e. the parsetree that
 *           maximizes the sum of the posterior probabilities of all
 *           1..L emitted residues, will be found.
 *
 *           Previously (infernal versions 1.0->1.0.2) this function
 *           (then named optimal_accuracy_align()) used the posterior
 *           matrix instead of the emit matrices used here, and thus
 *           did not determine (or at least was not guaranteed to
 *           determine) the optimally accurate parsetree as defined
 *           above. Instead it determined the parsetree that maximized
 *           the probability mass that passed through emitting states.
 *
 *           Local begins are handled the same as they are in
 *           cm_CYKInsideAlign(), see that function's purpose for specifics.
 *
 *           Note: Renamed from optimal_accuracy_align() [EPN, Wed Sep
 *	     14 06:16:38 2011].  Corrected to use emit matrices
 *	     instead of a posterior matrix [EPN, Sat Oct 1 06:04:34
 *	     2011].
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitaized sequence [i0..j0]   
 *           L         - length of the dsq
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx        - the DP matrix to fill in
 *           shmx      - the shadow matrix to fill in
 *           emit_mx   - pre-filled emit matrix
 *           ret_b     - RETURN: local begin state if local begins are on
 *           ret_pp    - RETURN: average posterior probability of aligned residues
 *                       in the optimally accurate parsetree
 *
 * Returns: <eslOK>     on success.
 * Throws:  <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *          If !eslOK: alignment has been aborted, ret_* variables are not valid
 */
int
cm_OptAccAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, CM_MX *mx, CM_SHADOW_MX *shmx, 
	       CM_EMIT_MX *emit_mx, int *ret_b, float *ret_pp)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  float    pp;		/* average posterior probability of all emitted residues */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;       /* j - sdr */
  int      d_sd;        /* d - sd */
  int      have_el;     /* TRUE if CM has local ends on, otherwise FALSE */
  int      c;           /* 64-bit counter */

  /* the DP matrices */
  float ***alpha   = mx->dp;       /* pointer to the alpha DP matrix, we'll store optimal parse in  */
  float  **l_pp    = emit_mx->l_pp; /* pointer to the prefilled posterior values for left  emitters */
  float  **r_pp    = emit_mx->r_pp; /* pointer to the prefilled posterior values for right emitters */
  char  ***yshadow = shmx->yshadow; /* pointer to the yshadow matrix */
  int   ***kshadow = shmx->kshadow; /* pointer to the kshadow matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_mx_GrowTo       (cm, mx,   errbuf, L, size_limit)) != eslOK) return status;
  if((status = cm_shadow_mx_GrowTo(cm, shmx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix */
  if(  mx->ncells_valid   > 0) esl_vec_FSet(mx->dp_mem, mx->ncells_valid, IMPOSSIBLE);
  if(shmx->y_ncells_valid > 0) for(c = 0; c < shmx->y_ncells_valid; c++) shmx->yshadow_mem[c] = USED_EL;
  /* for B states, shadow matrix holds k, length of right fragment, this will almost certainly be overwritten */
  if(shmx->k_ncells_valid > 0) esl_vec_ISet(shmx->kshadow_mem, shmx->k_ncells_valid, 0); 

  /* a special optimal accuracy specific step, initialize yshadow intelligently for d == 0 
   * (necessary b/c zero length parsetees have 0 emits and so always score IMPOSSIBLE)
   */
  if((status = cm_InitializeOptAccShadowDZero(cm, errbuf, yshadow, L)) != eslOK) return status;

  /* start with the EL state */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el && l_pp[cm->M] != NULL) { 
    for (j = 0; j <= L; j++) {
      alpha[cm->M][j][0] = l_pp[cm->M][0];
      i = j; 
      for (d = 1; d <= j; d++) { 
	alpha[cm->M][j][d] = FLogsum(alpha[cm->M][j][d-1], l_pp[cm->M][i--]);
      }
    }
  }

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);

    /* re-initialize if we can do a local end from v */
    if(have_el && NOT_IMPOSSIBLE(cm->endsc[v])) { 
      for (j = 0; j <= L; j++) {
	/* copy values from saved EL deck */
	for (d = sd; d <= j; d++) { 
	  alpha[v][j][d] = alpha[cm->M][j-sdr][d-sd];
	  /* yshadow[v][j][d] remains USED_EL */
	}
      }
    }
    /* note there's no E state update here, those cells all remain IMPOSSIBLE */

    /* we have to separate out IL_st and IR_st because IL use emit_mx->l_pp and IR use emit_mx->r_pp */
    if(cm->sttype[v] == IL_st) {
      /* update alpha[v][j][d] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (j = 1; j <= L; j++) {
	i    = j;
	d_sd = 0;
	for (d = 1; d <= j; d++, d_sd++, i--) {
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = alpha[y][j][d_sd]) > alpha[v][j][d]) {
	      alpha[v][j][d]   = sc; 
	      yshadow[v][j][d] = yoffset;
	    }
	  }
	  alpha[v][j][d] = FLogsum(alpha[v][j][d], l_pp[v][i]);
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
	  /* special case if local ends are off: explicitly disallow transitions to EL that require EL emissions 
	   * (we do allow an 'illegal' transition to EL in OptAcc but only if no EL emissions are req'd)
	   */
	  if((! have_el) && yshadow[v][j][d] == USED_EL && d > sd) { 
	    alpha[v][j][d] = IMPOSSIBLE;
	  }
	}
      }
    }
    if(cm->sttype[v] == IR_st) {
      /* IR: same loop nesting order as for IL for same reason, see IL comment above */
      j_sdr = 0;
      for (j = 1; j <= L; j++, j_sdr++) {
	d_sd  = 0;
	for (d = 1; d <= j; d++, d_sd++) {
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = cm->cfirst[v] + yoffset; 
	    if ((sc = alpha[y][j_sdr][d_sd]) > alpha[v][j][d]) {
	      alpha[v][j][d]   = sc; 
	      yshadow[v][j][d] = yoffset;
	    }
	  }
	  alpha[v][j][d] = FLogsum(alpha[v][j][d], r_pp[v][j]);
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
	  /* special case if local ends are off: explicitly disallow transitions to EL that require EL emissions 
	   * (we do allow an 'illegal' transition to EL in OptAcc but only if no EL emissions are req'd)
	   */
	  if((! have_el) && yshadow[v][j][d] == USED_EL && d > sd) { 
	    alpha[v][j][d] = IMPOSSIBLE;
	  }
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, so all cells
       * in alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }.
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	j_sdr = 0;
	for (j = sdr; j <= L; j++, j_sdr++) {
	  d_sd = 0;
	  for (d = sd; d <= j; d++, d_sd++) {
	    if((sc = alpha[y][j_sdr][d_sd]) > alpha[v][j][d]) {
	      alpha[v][j][d] = sc;
	      yshadow[v][j][d] = yoffset;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = 1; j <= L; j++) {
	  i = j;
	  for (d = sd; d <= j; d++, i--) {
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], l_pp[v][i]);
	  }
	}
	break;
      case MR_st:
	for (j = 1; j <= L; j++) {
	  for (d = sd; d <= j; d++) {
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], r_pp[v][j]);
	  }
	}
	break;
      case MP_st:
	for (j = 2; j <= L; j++) {
	  i = j-1;
	  for (d = sd; d <= j; d++, i--) {
	    alpha[v][j][d] = FLogsum(alpha[v][j][d], FLogsum(l_pp[v][i], r_pp[v][j])); 
	  }
	}
	break;
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = 0; j <= L; j++) {
	for (d = 0; d <= j; d++)
	  alpha[v][j][d] = ESL_MAX(alpha[v][j][d], IMPOSSIBLE);
      }
      /* special case if local ends are off: explicitly disallow transitions to EL that require EL emissions 
       * (we do allow an 'illegal' transition to EL in OptAcc but only if no EL emissions are req'd)
       */
      if(! have_el && sd > 0) { /* this is only necessary for emitters (MP, ML, MR in this context) */
	for (j = 0; j <= L; j++) {
	  for (d = sd+1; d <= j; d++) {
	    if(yshadow[v][j][d] == USED_EL) alpha[v][j][d] = IMPOSSIBLE;
	  }
	}
      }
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  for (k = 0; k <= d; k++) {
	    if ((sc = FLogsum(alpha[y][j-k][d-k], alpha[z][j][k])) > alpha[v][j][d])
	      {
		if(((d == k) || (NOT_IMPOSSIBLE(alpha[y][j-k][d-k]))) && /* left  subtree can only be IMPOSSIBLE if it has length 0 (in which case d==k, and d-k=0) */
		   ((k == 0) || (NOT_IMPOSSIBLE(alpha[z][j][k]))))       /* right subtree can only be IMPOSSIBLE if it has length 0 (in which case k==0) */
		  {
		    alpha[v][j][d]   = sc;
		    kshadow[v][j][d] = k;

		    /* Note: we take the logsum here, because we're keeping track of the
		     * log of the summed probability of emitting all residues up to this
		     * point, (from i..j) from left subtree (i=j-d+1..j-k) and from the 
		     * right subtree. (j-k+1..j)
		     * 
		     * EPN, Tue Nov 17 10:53:13 2009 Bug fix post infernal-1.0.2 release in "if(((sc = FLogsum..."
		     * statement above.  This is i15 in BUGTRAX, fixed as of svn revision 3056 in infernal 1.0 release
		     * branch, and revision 3057 in infernal trunk.  Bug description: See analogous section and comment
		     * in cm_OptAccAlignHB() above. In that function, in very rare cases (1 case in the 1.1 million SSU
		     * sequences in release 10_15 of RDP), this step will add two alpha values (alpha[y][j-k][d-k] for
		     * left subtree, and alpha[z][j][k] for right subtree) where one of them is IMPOSSIBLE and the
		     * corresponding subtree length ('d-k' in left subtree, or 'k' if right subtree) is non-zero, yet
		     * their FLogsum (which equals the value of the non-IMPOSSIBLE cell) is sufficiently high to be part
		     * of the optimally accurate traceback. This will probably cause a seg fault later b/c it implies a
		     * left or right subtree that is IMPOSSIBLE. It is okay if an IMPOSSIBLE scoring subtree has length
		     * 0 b/c 0 residues will contribute nothing to the summed log probability (nothing corresponds to a
		     * score of IMPOSSIBLE). We handle this case here by explicitly checking if either left or right
		     * subtree cell is IMPOSSIBLE with non-zero length before reassigning alpha[v][j][d].  I'm not sure
		     * if this is even possible in the non-banded function (this function), but I included the analogous
		     * fix here (the NOT_IMPOSSIBLE() calls) in case it was ever possible. This will slow down the
		     * implementation, but I'd rather err on the side of caution here, since we don't care so much about
		     * speed in the non-banded function, and b/c finding this bug again if the non-banded function can
		     * have the bug would be a pain in the ass.  
		     */
		  }
	      }
	  }
	}
      }
    }
    /* allow local begins, if nec */
    if((cm->flags & CMH_LOCAL_BEGIN) && (NOT_IMPOSSIBLE(cm->beginsc[v]))) {
      if (alpha[v][L][L] > bsc) {
	b   = v;
	bsc = alpha[v][L][L];
      }
    }
  } /* finished calculating deck v. */

  /* If local begins are on, the only way out of ROOT_S is via a local
   * begin, so update the optimal score and put a flag in the shadow
   * matrix telling cm_alignT() to use the b we return. 
   *
   * Note that because we're in OptAcc alpha[0][L][L] will already be
   * equal to bsc because transition scores (and thus impossible
   * transitions out of ROOT_S) have no effect on the score, so
   * whereas we can check to see if 'bsc > alpha[0][L][L]' at an
   * analogous point in CYK before setting the USED_LOCAL_BEGIN
   * flag, we can't here because it would be FALSE.
   */
  if(NOT_IMPOSSIBLE(bsc) && (cm->flags & CMH_LOCAL_BEGIN)) {
    alpha[0][L][L]   = bsc;
    yshadow[0][L][L] = USED_LOCAL_BEGIN;
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_oamx",   "w"); cm_mx_Dump(fp1, mx); fclose(fp1); */
  /* FILE *fp2; fp2 = fopen("tmp.std_oashmx", "w"); cm_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2); */ 
#endif
  
  sc = alpha[0][L][L];

  /* convert sc, a log probability, into the average posterior probability of all L aligned residues */
  pp = sreEXP2(sc) / (float) L;

  if(ret_b != NULL)  *ret_b  = b;    /* b is -1 if local ends are off */
  if(ret_pp != NULL) *ret_pp = pp;

  ESL_DPRINTF1(("#DEBUG: cm_OptAccAlign return pp: %f\n", pp));
  return eslOK;
}


/* Function: cm_OptAccAlignHB()
 * Date:     EPN, Thu Nov 15 10:48:37 2007
 *
 * Purpose:  Same as cm_OptAccAlign() but HMM bands are used.
 *           See cm_OptAccAlign()'s Purpose for more information.
 *
 *           Note: Renamed from optimal_accuracy_align_hb() [EPN, Wed Sep 14 06:16:06 2011].
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitaized sequence [i0..j0]   
 *           L         - length of the dsq
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx        - the DP matrix to fill in
 *           shmx      - the shadow matrix to fill in
 *           emit_mx   - pre-filled emit matrix
 *           ret_b     - RETURN: local begin state if local begins are on
 *           ret_pp    - RETURN: average posterior probability of aligned residues
 *                       in the optimally accurate parsetree
 *
 * Returns: <eslOK> on success.
 * 
 * Throws:  <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *          If !eslOK: alignment has been aborted, ret_* variables are not valid
 */
int
cm_OptAccAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, 
		 CM_HB_EMIT_MX *emit_mx, int *ret_b, float *ret_pp)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  float    pp;		/* average posterior probability of all emitted residues */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  int      jp_0;        /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;        /* L offset in ROOT_S's (v==0) d band */
  int      Lp;          /* L offset in any state v's d band */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      have_el;     /* TRUE if CM has local ends on, otherwise FALSE */
  int64_t  c;           /* 64-bit int counter */

  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      ip_v;               /* offset i index for state v */
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
  int      jp_y_minus_k;       /* jp_y - k, used in one loop, stored to avoid calc'ing twice */
  int      dp_y_minus_k;       /* dp_y - k, used in one loop, stored to avoid calc'ing twice */
  int      kn, kx;             /* current minimum/maximum k value */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b  = cm->cp9b;
  int        *imin  = cp9b->imin;  
  int        *imax  = cp9b->imax;
  int        *jmin  = cp9b->jmin;  
  int        *jmax  = cp9b->jmax;
  int       **hdmin = cp9b->hdmin;
  int       **hdmax = cp9b->hdmax;

  /* the DP matrices */
  float ***alpha   = mx->dp;        /* pointer to the alpha DP matrix, we'll store optimal parse in  */
  float  **l_pp    = emit_mx->l_pp; /* pointer to the prefilled posterior values for left  emitters */
  float  **r_pp    = emit_mx->r_pp; /* pointer to the prefilled posterior values for right emitters */
  char  ***yshadow = shmx->yshadow; /* pointer to the yshadow matrix */
  int   ***kshadow = shmx->kshadow; /* pointer to the kshadow matrix */

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  if (cp9b->jmin[0] > L || cp9b->jmax[0] < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_OptAccAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, cp9b->jmin[0], cp9b->jmax[0]);
  jp_0 = L - jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) 
    ESL_FAIL(eslEINVAL, errbuf, "cm_OptAccAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, hd_min(cp9b, 0, jp_0), hd_max(cp9b, 0, jp_0));
  Lp_0 = L - hd_min(cp9b, 0, jp_0);

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo       (cm, mx,   errbuf, cp9b, L, size_limit)) != eslOK) return status;
  if((status = cm_hb_shadow_mx_GrowTo(cm, shmx, errbuf, cp9b, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix */
  if(  mx->ncells_valid   > 0) esl_vec_FSet(mx->dp_mem, mx->ncells_valid, IMPOSSIBLE);
  if(shmx->y_ncells_valid > 0) for(c = 0; c < shmx->y_ncells_valid; c++) shmx->yshadow_mem[c] = USED_EL;
  /* for B states, shadow matrix holds k, length of right fragment, this will almost certainly be overwritten */
  if(shmx->k_ncells_valid > 0) esl_vec_ISet(shmx->kshadow_mem, shmx->k_ncells_valid, 0); 

  /* a special optimal accuracy specific step, initialize yshadow intelligently for d == 0 
   * (necessary b/c zero length parsetees have 0 emits and so always score IMPOSSIBLE)
   */
  if((status = cm_InitializeOptAccShadowDZeroHB(cm, cp9b, errbuf, yshadow, L)) != eslOK) return status;

  /* start with the EL state (remember, cm->M deck is non-banded) */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el && l_pp[cm->M] != NULL) { 
    for (j = 0; j <= L; j++) {
      alpha[cm->M][j][0] = l_pp[cm->M][0];
      i = j; 
      for (d = 1; d <= j; d++) { 
	alpha[cm->M][j][d] = FLogsum(alpha[cm->M][j][d-1], l_pp[cm->M][i--]);
      }
    }
  }

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) {
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);

    /* re-initialize if we can do a local end from v */
    if(have_el && NOT_IMPOSSIBLE(cm->endsc[v])) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	/* copy values from saved EL deck */
	for(d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	  dp_v = d - hd_min(cp9b, v, jp_v);
	  alpha[v][jp_v][dp_v] = alpha[cm->M][j-sdr][d-sd];
	  /* yshadow[v][jp_v][dp_v] remains USED_EL */
	}
      }
    }
    /* note there's no E state update here, those cells all remain IMPOSSIBLE */

    /* we could separate out IL_st and IR_st, but I don't it makes a significant difference in run time */
    if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
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
	
	i = j - hd_min(cp9b, v, jp_v) + 1;
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++, i--) { /* for each valid d for v, j */
	  /*printf("v: %4d j: %4d (%4d..%4d) d: %4d (%4d..%4d) i: %4d (%4d..%4d)\n", 
	    v, j, jmin[v], jmax[v], d, hd_min(cp9b, v, jp_v), hd_max(cp9b, v, jp_v), i, imin[v], imax[v]);*/
	  assert(i >= imin[v] && i <= imax[v]);
	  ESL_DASSERT1((i >= imin[v] && i <= imax[v]));
	  ip_v = i - imin[v];         /* i index for state v in emit_mx->l_pp */
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hd_min(cp9b, y, jp_y_sdr) && (d-sd) <= hd_max(cp9b, y, jp_y_sdr)) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hd_min(cp9b, y, jp_y_sdr);
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd]) > alpha[v][jp_v][dp_v])
		{
		  alpha[v][jp_v][dp_v] = sc; 
		  yshadow[v][jp_v][dp_v]    = yoffset;
		}
	    }
	  }
	  if(cm->sttype[v] == IL_st) alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], l_pp[v][ip_v]);
	  else                       alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], r_pp[v][jp_v]);
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	  /* special case if local ends are off: explicitly disallow transitions to EL that require EL emissions 
	   * (we do allow an 'illegal' transition to EL in OptAcc but only if no EL emissions are req'd)
	   */
	  if((! have_el) && yshadow[v][jp_v][dp_v] == USED_EL && d > sd) { 
	    alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	  }
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
	
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	
	for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
	  ESL_DASSERT1((jp_v     >= 0 && jp_v     <= (jmax[v]-jmin[v])));
	  ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	  
	  dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y_sdr) + sd);
	  dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y_sdr) + sd);
	  dpn     = dn - hd_min(cp9b, v, jp_v);
	  dpx     = dx - hd_min(cp9b, v, jp_v);
	  dp_y_sd = dn - hd_min(cp9b, y, jp_y_sdr) - sd;
	  	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hd_max(cp9b, v, jp_v)     - hd_min(cp9b, v, jp_v))));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hd_max(cp9b, y, jp_y_sdr) - hd_min(cp9b, y, jp_y_sdr))));
	    if((sc = alpha[y][jp_y_sdr][dp_y_sd]) > alpha[v][jp_v][dp_v]) {
	      alpha[v][jp_v][dp_v] = sc;
	      yshadow[v][jp_v][dp_v] = yoffset;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v = j - jmin[v];
	  i    = j - hd_min(cp9b, v, jp_v) + 1;
	  ip_v = i - imin[v];
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++, ip_v--) { 
	    /*printf("v: %4d j: %4d (%4d..%4d) d: %4d (%4d..%4d) i: %4d (%4d..%4d)\n", 
	      v, j, jmin[v], jmax[v], d, hd_min(cp9b, v, jp_v), hd_max(cp9b, v, jp_v), i, imin[v], imax[v]);*/
	    assert(ip_v >= 0 && ip_v <= (imax[v] - imin[v]));
	    ESL_DASSERT1((ip_v >= 0 && ip_v <= (imax[v] - imin[v])));
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], l_pp[v][ip_v]);
	  }
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++) { 
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], r_pp[v][jp_v]);
	  }
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v = j - jmin[v];
	  i    = j - hd_min(cp9b, v, jp_v) + 1;
	  ip_v = i - imin[v];
	  for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++, ip_v--) {
	    assert(ip_v >= 0 && ip_v <= (imax[v] - imin[v]));
	    ESL_DASSERT1((ip_v >= 0 && ip_v <= (imax[v] - imin[v])));
	    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], FLogsum(l_pp[v][ip_v], r_pp[v][jp_v])); 
	  }
	}
	break;
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j - jmin[v];
	for (dp_v = 0; dp_v <= (hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v)); dp_v++)
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
      }
      /* special case if local ends are off: explicitly disallow transitions to EL that require EL emissions 
       * (we do allow an 'illegal' transition to EL in OptAcc but only if no EL emissions are req'd)
       */
      if(! have_el && sd > 0) { /* this is only necessary for emitters (MP, ML, MR in this context) */
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v = j - jmin[v];
	  d = ESL_MAX(sd+1, hd_min(cp9b, v, jp_v));
	  dp_v = d - hd_min(cp9b, v, jp_v);
	  for (; d <= hd_max(cp9b, v, jp_v); d++) { 
	    if(yshadow[v][jp_v][dp_v] == USED_EL) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	    dp_v++;
	  }
	}
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
	kn = ((j-jmax[y]) > (hd_min(cp9b, z, jp_z))) ? (j-jmax[y]) : hd_min(cp9b, z, jp_z);
        kn = ESL_MAX(kn, 0); /* kn must be non-negative, added with fix to bug i36 */
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hd_max(cp9b, z, jp_z))) ?  jp_y       : hd_max(cp9b, z, jp_z);
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
	      
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
	    if((k >= d - hd_max(cp9b, y, jp_y-k)) && k <= d - hd_min(cp9b, y, jp_y-k)) {
	      /* for current k, all 6 inequalities have been satisified 
	       * so we know the cells corresponding to the platonic 
	       * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
	       * alpha[z][j][k] are all within the bands. These
	       * cells correspond to alpha[v][jp_v][dp_v], 
	       * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
	       * and alpha[z][jp_z][k-hdmin[jp_z]];
	       */
	      kp_z = k-hd_min(cp9b, z, jp_z);
	      dp_y = d-hd_min(cp9b, y, jp_y-k);
	      jp_y_minus_k = jp_y-k;
	      dp_y_minus_k = dp_y-k;

	      if((sc = FLogsum(alpha[y][jp_y_minus_k][dp_y_minus_k], alpha[z][jp_z][kp_z])) > alpha[v][jp_v][dp_v]) 
		{
		  if(((d == k) || (NOT_IMPOSSIBLE(alpha[y][jp_y_minus_k][dp_y_minus_k]))) && /* left subtree can only be IMPOSSIBLE if it has length 0 (in which case d==k, and d-k=0) */
		     ((k == 0) || (NOT_IMPOSSIBLE(alpha[z][jp_z][kp_z]))))                   /* right subtree can only be IMPOSSIBLE if it has length 0 (in which case k==0) */
		    { 
		      alpha[v][jp_v][dp_v] = sc;
		      kshadow[v][jp_v][dp_v] = k;
		      /* Note: we take the logsum here, because we're
		       * keeping track of the log of the summed probability
		       * of emitting all residues up to this point, (from
		       * i..j) from left subtree (i=j-d+1..j-k) and from the
		       * right subtree. (j-k+1..j).  
		       * 
		       * EPN, Tue Nov 17 09:57:59 2009: 
		       * Bug fix post infernal-1.0.2 release.
		       * Addition of 2-line if statement beginning
		       * "if(((d == k...)"  This is i15 in BUGTRAX,
		       * fixed as of svn revision 3056 in infernal 1.0
		       * release branch, and revision 3057 in infernal
		       * trunk.  
		       * Bug description: In very rare cases (1 case
		       * in the 1.1 million SSU sequences in release
		       * 10_15 of RDP), this step will add two alpha
		       * values (alpha[y][jp_y_minus_k][dp_y_minus_k]
		       * for left subtree, and alpha[z][jp_z][kp_z]
		       * for right subtree) where one of them is
		       * IMPOSSIBLE and the corresponding subtree
		       * length ('d-k' in left subtree, or 'k' if right
		       * subtree) is non-zero, yet their FLogsum
		       * (which equals the value of the non-IMPOSSIBLE
		       * cell) is sufficiently high to be part of the
		       * optimally accurate traceback. This will
		       * probably cause a seg fault later b/c it
		       * implies a left or right subtree that is
		       * IMPOSSIBLE. It is okay if an IMPOSSIBLE
		       * scoring subtree has length 0 b/c 0 residues
		       * will contribute nothing to the summed log
		       * probability (nothing corresponds to a score
		       * of IMPOSSIBLE). We handle this case here by
		       * explicitly checking if either left or right
		       * subtree cell is IMPOSSIBLE with non-zero
		       * length before reassigning
		       * alpha[v][jp_v][dp_v].
		       */
		    }
		}
	    }
	  }
	}
      }
    }
    /* allow local begins, if nec */
    if((cm->flags & CMH_LOCAL_BEGIN) && (NOT_IMPOSSIBLE(cm->beginsc[v]))) {
      if(L >= jmin[v] && L <= jmax[v]) { 
	jp_v = L - jmin[v];
	Lp   = L - hd_min(cp9b, v, jp_v);
	if(L >= hd_min(cp9b, v, jp_v) && L <= hd_max(cp9b, v, jp_v)) { 
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
	  if (alpha[v][jp_v][Lp] > bsc) { 
	    b   = v;
	    bsc = alpha[v][jp_v][Lp];
	  }
	}
      }
    }
  } /* end loop over all v */

  /* If local begins are on, the only way out of ROOT_S is via a local
   * begin, so update the optimal score and put a flag in the shadow
   * matrix telling cm_alignT() to use the b we return. 
   *
   * Note that because we're in OptAcc alpha[0][L][L] will already be
   * equal to bsc because transition scores (and thus impossible
   * transitions out of ROOT_S) have no effect on the score, so
   * whereas we can check to see if 'bsc > alpha[0][L][L]' at an
   * analogous point in CYK before setting the USED_LOCAL_BEGIN
   * flag, we can't here because it would be FALSE.
   */
  if(NOT_IMPOSSIBLE(bsc) && (cm->flags & CMH_LOCAL_BEGIN)) {
    alpha[0][jp_0][Lp_0]   = bsc;
    yshadow[0][jp_0][Lp_0] = USED_LOCAL_BEGIN;
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_oahbmx",   "w"); cm_hb_mx_Dump(fp1, mx); fclose(fp1); */
  /* FILE *fp2; fp2 = fopen("tmp.std_oahbshmx", "w"); cm_hb_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2); */
#endif
  
  sc = alpha[0][jp_0][Lp_0];

  /* convert sc, a log probability, into the average posterior probability of all L aligned residues */
  pp = sreEXP2(sc) / (float) L;

  free(yvalidA);

  if (ret_b  != NULL)   *ret_b  = b;   /* b is -1 if local begins are off */
  if (ret_pp != NULL)   *ret_pp = pp;  

  ESL_DPRINTF1(("#DEBUG: cm_OptAccAlignHB return pp: %f\n", pp));
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
  return status; /* never reached */
}

/* Function: cm_CYKOutsideAlign()
 * Date:     EPN, Wed Sep 14 14:01:36 2011
 *
 * Purpose:  Run the outside CYK algorithm on a target sequence.
 *           Non-banded version. See cm_CYKOutsideAlignHB() for
 *           the HMM banded version. The full target sequence
 *           1..L is aligned. 
 * 
 *           Very similar to cm_OutsideAlign() but calculates
 *           beta[v][j][d]: log probability of the most likely parse
 *           that emits 1..i-1 and j+1..L and passes through v at j,d
 *           (where i = j-d+1) instead of the log of the summed
 *           probability of all such parses. This means max operations
 *           are used instead of logsums.
 *
 *           This function complements cm_CYKInsideAlign() but is
 *           mainly useful for testing and reference. It can be used
 *           with do_check=TRUE to verify that the implementation of
 *           CYKInside and CYKOutside are consistent.  Because the
 *           structure of CYKInside and Inside, and CYKOutside and
 *           Outside are so similar and the CYK variants are easier to
 *           debug (because only the optimal parsetree is considered
 *           instead of all possible parsetrees) this function can be
 *           useful for finding bugs in Outside.  It is currently not
 *           hooked up to any of the main Infernal programs.
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           L         - length of the dsq to align
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_check  - TRUE to attempt to check 
 *           mx        - the dp matrix, grown and filled here
 *           inscyk_mx - the pre-filled dp matrix from the CYK Inside calculation 
 *                       (performed by cm_CYKInsideAlign(), required)
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score, this is from 
 *                       inscyk_mx IF local ends are on (see comments towards 
 *                       end of function).
 *
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *           <eslEMEM>   if we run out of memory
 *           <eslFAIL>   if <do_check>==TRUE and we fail a test
 *           In any of these cases, alignment has been aborted, ret_sc is not valid.
 */
int 
cm_CYKOutsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check,
		   CM_MX *mx, CM_MX *inscyk_mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float    sc;     	       /* a temporary score */
  float  **esc_vAA;            /* ptr to cm->oesc, optimized emission scores */
  float    escore;	       /* an emission score, tmp variable */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */

  /* variables used only if do_check==TRUE */
  int      fail1_flag = FALSE; /* set to TRUE if do_check and we see a problem with check 1*/
  int      fail2_flag = FALSE; /* set to TRUE if do_check and we see a problem with check 2*/
  int      fail3_flag = FALSE; /* set to TRUE if do_check and we see a problem with check 3*/
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */
  int      vmax;               /* i, offset in the matrix */
  float    tol;                /* tolerance for differences in bit scores */
  int     *optseen = NULL;     /* [1..i..W] TRUE is residue i is accounted for in optimal parse */

  /* the DP matrices */
  float ***beta  = mx->dp;        /* pointer to the Oustide DP mx */
  float ***alpha = inscyk_mx->dp; /* pointer to the CYK Inside DP mx (already calc'ed and passed in) */

  /* Allocations and initializations */
  esc_vAA = cm->oesc;            /* a ptr to the optimized emission scores */

  /* grow the matrix based on the current sequence */
  if((status = cm_mx_GrowTo(cm, mx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(beta[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* now set beta[0][L][L] to 0., all (valid) parses must end there */
  beta[0][L][L] = 0.;

  /* initialize local begin cells for emitting full seq (j==L && d == L) */
  if (cm->flags & CMH_LOCAL_BEGIN) { 
    for (v = 1; v < cm->M; v++) 
      beta[v][L][L] = cm->beginsc[v];
  }

  /* Main recursion */
  for (v = 1; v < cm->M; v++) { /* start at state 1 because we set all values for ROOT_S state 0 above */
    sd  = StateDelta(cm->sttype[v]);
    sdr = StateRightDelta(cm->sttype[v]);

    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      y = cm->plast[v];	/* the parent bifurcation    */
      z = cm->cnum[y];	/* the other (right) S state */
      for(j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  for (k = 0; k <= (L-j); k++) {
	    beta[v][j][d] = ESL_MAX(beta[v][j][d], (beta[y][j+k][d+k] + alpha[z][j+k][k]));
	  }
	}
      }
    } /* end of 'if (cm->stid[v] == BEGL_S */
    else if (cm->stid[v] == BEGR_S) {
      y = cm->plast[v];	  /* the parent bifurcation    */
      z = cm->cfirst[y];  /* the other (left) S state  */
      for(j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  for (k = 0; k <= (j-d); k++) {
 	    beta[v][j][d] = ESL_MAX(beta[v][j][d], (beta[y][j][d+k] + alpha[z][j-d][k]));
	  }
	}
      }
    } /* end of 'else if (cm->stid[v] == BEGR_S */
    else { /* (cm->sttype[v] != BEGL_S && cm->sttype[v] != BEGR_S */ 
      for (j = L; j >= 0; j--) {
	i = 1;
	for (d = j; d >= 0; d--, i++) {
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    sd  = StateDelta(cm->sttype[y]);
	    sdr = StateRightDelta(cm->sttype[y]);
	    switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == L || d == j) continue; /* boundary condition */
		escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		beta[v][j][d] = ESL_MAX(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;

	      case ML_st:
	      case IL_st: 
		if (d == j) continue;	/* boundary condition (note when j=0, d=0*/
		escore = esc_vAA[y][dsq[i-1]];
		beta[v][j][d] = ESL_MAX(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == L) continue;
		escore = esc_vAA[y][dsq[j+1]];
		beta[v][j][d] = ESL_MAX(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		beta[v][j][d] = ESL_MAX(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset]));
		break;
	    } /* end of switch(cm->sttype[y] */  
	  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	  if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;
	} /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
      } /* end loop over j. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] != BEGL_S, BEGR_S */
    /* we're done calculating deck v for everything but local begins */

    /* deal with local alignment end transitions v->EL (EL = deck at M.) */
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  i = j-d+1;
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    if (j == L || d == j) continue; /* boundary condition */
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v] 
							    + escore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == j) continue;	
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v] 
							    + escore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == L) continue;
	    escore = esc_vAA[v][dsq[j+1]];
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v]
							    + escore));
	    break;
	  case S_st:
	  case D_st:
	  case B_st:
	  case E_st:
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v]));
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
    for (j = L; j > 0; j--) { /* careful w/ boundary here */
      for (d = j-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[cm->M][j][d+1] + cm->el_selfsc));
    }
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_ocykmx", "w");   cm_mx_Dump(fp1, mx); fclose(fp1); */
#endif

  fail1_flag = FALSE;
  fail2_flag = FALSE;
  fail3_flag = FALSE;
  if(do_check) {
    /* Check for consistency between the Inside alpha matrix and the
     * Outside beta matrix. We assume the Inside CYK parse score
     * (optsc) is the optimal score, so for all v,j,d:
     * 
     * Jalpha[v][j][d] + Jbeta[v][j][d] <= optsc
     * 
     * Further, we know that each residue must be emitted by a state
     * in the optimal parse. So as we do the above check, we determine
     * when we're in a cell that may be involved in the optimal parse
     * (the sum of the Inside and Outside scores are equal to the
     * optimal parse score), if that cell corresponds to a left
     * emitter emitting position i, we know an emitted i has been
     * observed in an optimal parse and set optseen[i] to TRUE.
     * Likewise, if that cell corresponds to a right emitter emitting
     * position j, we update optseen[j] to TRUE. At the end of the
     * check optseen[i] should be TRUE for all i in the range
     * [1..L].
     *
     * Note that we don't ensure that all of our presumed optimal
     * cells make up a valid parse, so it is possible we could pass
     * this check even if the Inside and Outside matrices are
     * inconsistent (i.e. there's a bug in the implementation of one
     * and or the other) but that should be extremely unlikely.  If we
     * do this test many times for many different models and pass, we
     * should be confident we have consistent implementations.
     * 
     * This is an expensive check and should only be done while
     * debugging.
     *
     * Another test we could do but do not is to determine the CYK
     * parse by tracing back the CYK Inside matrix, then ensure that
     * for each cell in that parse alpha[v][j][d]+beta[v][j][d] ==
     * optsc.
     */
    ESL_ALLOC(optseen, sizeof(int) * (L+1));
    esl_vec_ISet(optseen, L+1, FALSE);
    vmax = (cm->flags & CMH_LOCAL_END) ? cm->M : cm->M-1;
    /* define bit score difference tolerance, somewhat arbitrarily: 
     * clen <= 200: tolerance is 0.001; then a function of clen: 
     * clen == 1000 tolerance is 0.005, 
     * clen == 2000, tolerance is 0.01.
     *
     * I did this b/c with tests with SSU_rRNA_eukarya I noticed
     * failures with bit score differences up to 0.004 or so.  This
     * could mean a bug, but I couldn't get any average sized model to
     * fail with a difference above 0.001, so I blamed it on
     * precision. I'm not entirely convinced it isn't a bug but
     * until I see a failure on a smaller model it seems precision
     * is the most likely explanation, right?  
     */ 
    tol = ESL_MAX(1e-3, (float) cm->clen / 200000.); 
    for(v = 0; v <= vmax; v++) { 
      for(j = 1; j <= L; j++) { 
	for(d = 0; d <= j; d++) { 
	  sc  = (alpha[v][j][d] + beta[v][j][d]) - alpha[0][L][L];
	  if(sc > tol) { 
	    printf("Check 1 failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, alpha[v][j][d], beta[v][j][d], alpha[v][j][d] + beta[v][j][d], alpha[0][L][L]);
	    fail1_flag = TRUE;
	  }
	  if(fabs(sc) < tol) { /* this cell is involved in a parse with the optimal score */
	    i  = j-d+1;
	    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st || (cm->sttype[v] == EL_st && d >0)) { 
	      /* i is accounted for by a parse with an optimal score */
	      optseen[i] = TRUE;
	      /*printf("\tResidue %4d possibly accounted for by Left  emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", i, Statetype(cm->sttype[v]), v, j, d);*/
	    }
	    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
	      /* j is accounted for by a parse with an optimal score */
	      optseen[j] = TRUE;
	      /*printf("\tResidue %4d possibly accounted for by Right emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", j, Statetype(cm->sttype[v]), v, j, d);*/
	    }
	  }
	}
      }
    }
    for(j = 1; j <= L; j++) { 
      if(optseen[j] == FALSE) { 
	printf("Check 2 failure: residue %d not emitted in the optimal parsetree\n", j);
	fail2_flag = TRUE;
      }	      
    }
    free(optseen);
  }
  /* Another test that we can only do if local ends are OFF */
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     *
     * Determine P(pi, S|M) / P(S|R) (probability of the sequence and most likely parse
     * tree pi given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(pi, S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Max_v [ Max_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                = P(pi, S|M) / P(S|R)
     */
    for(n = 0; n < cm->nodes; n++) {
      sc = IMPOSSIBLE;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = 0; j <= L; j++) {
	  for (d = 0; d <= j; d++) {
	    sc = ESL_MAX(sc, (alpha[v][j][d] + beta[v][j][d]));
	  }
	}
      }
      /*printf("checking node: %d | sc: %.6f\n", n, sc);*/
      diff = sc - alpha[0][L][L];
      if(diff > 0.01 || diff < -0.01) { 
	fail3_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, sc, alpha[0][L][L], diff);
      }
    }
  }
  /* Finally, calculate the optimal score, but this only works if
   * we're not in local mode:
   * 
   * If local ends are off, we know the optimal parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Max_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Max_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * *** If local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END)) { 
    sc = IMPOSSIBLE;
    v = cm->M-1;
    for (j = 0; j <= L; j++) {
      sc = ESL_MAX(sc, (beta[v][j][0]));
      /*printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
    }
  }
  else { /* return sc = P(S|M) / P(S|R) from Inside() */
    sc = alpha[0][L][L];
  }

  if(do_check) { 
    if     (fail1_flag) ESL_FAIL(eslFAIL, errbuf, "CYK Inside/Outside check1 FAILED.");
    else if(fail2_flag) ESL_FAIL(eslFAIL, errbuf, "CYK Inside/Outside check2 FAILED.");
    else if(fail3_flag) ESL_FAIL(eslFAIL, errbuf, "CYK Inside/Outside check3 FAILED.");
    ESL_DPRINTF1(("#DEBUG: SUCCESS! CYK Inside/Outside checks PASSED.\n"));
  }

  if(!(cm->flags & CMH_LOCAL_END)) ESL_DPRINTF1(("#DEBUG: \tcm_CYKOutsideAlign() sc : %f\n", sc));
  else                             ESL_DPRINTF1(("#DEBUG: \tcm_CYKOutsideAlign() sc : %f (LOCAL mode; sc is from Inside)\n", sc));

  if(ret_sc != NULL) *ret_sc = sc;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Out of memory");
  return status; /* NEVER REACHED */
}

/* Function: cm_CYKOutsideAlignHB()
 * Date:     EPN, Fri Sep 30 10:12:51 2011
 *
 * Purpose:  Run the outside CYK algorithm on a target sequence.
 *           HMM banded version. See cm_CYKOutsideAlign() for
 *           the non-banded version. The full target sequence
 *           1..L is aligned. 
 *
 *           Very similar to cm_OutsideAlignHB() but calculates
 *           beta[v][j][d]: log probability of the most likely parse
 *           that emits 1..i-1 and j+1..L and passes through v at j,d
 *           (where i = j-d+1) instead of the log of the summed
 *           probability of all such parses. This means max operations
 *           are used instead of logsums.
 *
 *           This function complements cm_CYKInsideAlignHB() but is
 *           mainly useful for testing and reference. It can be used
 *           with do_check=TRUE to verify that the implementation of
 *           CYKInsideHB and CYKOutsideHB are consistent.  Because the
 *           structure of CYKInsideHB and InsideHB, and CYKOutsideHB
 *           and OutsideHB are so similar and the CYK variants are
 *           easier to debug (because only the optimal parsetree is
 *           considered instead of all possible parsetrees) this
 *           function can be useful for finding bugs in OutsideHB.  It
 *           is currently not hooked up to any of the main Infernal
 *           programs.
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           L         - length of the dsq to align
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_check  - TRUE to attempt to check 
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ins_mx    - the dp matrix from the Inside run calculation (required)
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *                       ends are on (see *** comment towards end of function).
 *
 * Returns:  <eslOK> on success
 *
 * Throws:   <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *           <eslFAIL>   if <do_check>==TRUE and we fail a test
 *           In either of these cases, alignment has been aborted, ret_sc is not valid.
 */
int
cm_CYKOutsideAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, 
		     CM_HB_MX *mx, CM_HB_MX *ins_mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float  **esc_vAA;            /* ptr to cm->oesc, optimized emission scores */
  float    sc;		       /* a temporary score */
  float    escore;	       /* an emission score, tmp variable */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */

  /* variables used only if do_check */
  int      fail1_flag = FALSE; /* set to TRUE if do_check and we see a problem with check 1*/
  int      fail2_flag = FALSE; /* set to TRUE if do_check and we see a problem with check 2*/
  int      fail3_flag = FALSE; /* set to TRUE if do_check and we see a problem with check 3*/
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */
  int      vmax;               /* i, offset in the matrix */
  float    tol;                /* tolerance for differences in bit scores */
  int     *optseen = NULL;     /* [1..i..W] TRUE is residue i is accounted for in optimal parse */

  /* band related variables */
  int      dp_v;               /* d index for state v in alpha w/mem eff bands */
  int      dp_y;               /* d index for state y in alpha w/mem eff bands */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Lp;                 /* L index also changes depending on state */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      kmin, kmax;         /* temporary minimum/maximum allowed k */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      jp_0;               /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;               /* L offset in ROOT_S's (v==0) d band */

  /* the DP matrices */
  float ***beta  = mx->dp;     /* pointer to the Oustide DP mx */
  float ***alpha = ins_mx->dp; /* pointer to the Inside DP mx (already calc'ed and passed in) */

  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cm->cp9b->jmin;  
  int     *jmax  = cm->cp9b->jmax;
  CP9Bands_t *cp9b = cm->cp9b;  /* brief 26_0430-157: needed by hd_min()/hd_max() */
  int    **hdmin = cm->cp9b->hdmin;
  int    **hdmax = cm->cp9b->hdmax;

  /* Allocations and initializations */
  esc_vAA = cm->oesc;            /* a ptr to the optimized emission scores */

  /* grow the matrix based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo(cm, mx, errbuf, cm->cp9b, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(beta[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* ensure a full alignment to ROOT_S (v==0) is allowed by the bands */
  if (jmin[0] > L || jmax[0] < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, jmin[0], jmax[0]);
  jp_0 = L - jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L) 
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, hd_min(cp9b, 0, jp_0), hd_max(cp9b, 0, jp_0));
  Lp_0 = L - hd_min(cp9b, 0, jp_0);
  /* set the offset banded cell corresponding to beta[0][L][L] to 0., all parses must end there */
  beta[0][jp_0][Lp_0] = 0.;

  /* If we can do a local begin into v, overwrite IMPOSSIBLE with the local begin score. */
  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (v = 1; v < cm->M; v++) {
      if(NOT_IMPOSSIBLE(cm->beginsc[v])) {
	if((L >= jmin[v]) && (L <= jmax[v])) {
	  jp_v = L - jmin[v];
	  if((L >= hd_min(cp9b, v, jp_v)) && L <= hd_max(cp9b, v, jp_v)) {
	    Lp = L - hd_min(cp9b, v, jp_v);
	    beta[v][jp_v][Lp] = cm->beginsc[v];
	  }
	}
      }
    }
  }
  /* done allocation/initialization */

  /* Recursion: main loop down through the decks */
  for (v = 1; v < cm->M; v++) { /* start at state 1 because we set all values for ROOT_S state 0 above */
    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      y = cm->plast[v];	/* the parent bifurcation    */
      z = cm->cnum[y];	/* the other (right) S state */
      for (j = jmax[v]; j >= jmin[v]; j--) {
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;
	for (d = hd_max(cp9b, v, jp_v); d >= hd_min(cp9b, v, jp_v); d--) {
	  dp_v = d - hd_min(cp9b, v, jp_v);
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
	    if(k < (hd_min(cp9b, y, jp_y+k) - d) || k > (hd_max(cp9b, y, jp_y+k) - d)) continue; 
	    /* above line continues if inequality 5 or 6 is violated */
	    if(k < (hd_min(cp9b, z, jp_z+k))     || k > (hd_max(cp9b, z, jp_z+k)))     continue; 
	    /* above line continues if inequality 7 or 8 is violated */
		  
	    /* if we get here for current k, all 8 inequalities have been satisified 
	     * so we know the cells corresponding to the platonic 
	     * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
	     * alpha[z][j+k][k] are all within the bands. These
	     * cells correspond to beta[v][jp_v][dp_v], 
	     * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
	     * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
	     */
	    kp_z = k-hd_min(cp9b, z, jp_z+k);
	    dp_y = d-hd_min(cp9b, y, jp_y+k);
	    beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] + alpha[z][jp_z+k][kp_z]));
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
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;

	dn = ESL_MAX(hd_min(cp9b, v, jp_v), j-jmax[z]);
	dx = ESL_MIN(hd_max(cp9b, v, jp_v), jp_z);
	/* above makes sure that j,d are valid for state z: (jmin[z] + d) >= j >= (jmax[z] + d) */
	for (d = dx; d >= dn; d--) {
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
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
	  kmin = ESL_MAX((hd_min(cp9b, y, jp_y)-d), (hd_min(cp9b, z, jp_z-d)));
	  kmax = ESL_MIN((hd_max(cp9b, y, jp_y)-d), (hd_max(cp9b, z, jp_z-d)));
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
	    kp_z = k-hd_min(cp9b, z, jp_z-d);
	    dp_y = d-hd_min(cp9b, y, jp_y);
	    beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
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
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	for (d = hd_max(cp9b, v, jp_v); d >= hd_min(cp9b, v, jp_v); d--) {
	  i = j-d+1;
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
	  
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    
	    /* Note: this looks like it can be optimized, I tried but my 'optimization' slowed the code, so I reverted [EPN] */
	    switch(cm->sttype[y]) {
	    case MP_st: 
	      if (j == L || d == j) continue; /* boundary condition */
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+2) < hd_min(cp9b, y, (jp_y+1)) || (d+2) > hd_max(cp9b, y, (jp_y+1))) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
		   */
	      dp_y = d - hd_min(cp9b, y, jp_y+1);  /* d index for state y */
	      escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case ML_st:
	    case IL_st: 
	      if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
	      if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hd_min(cp9b, y, jp_y) || (d+1) > hd_max(cp9b, y, jp_y)) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j][d+1] in the platonic matrix.
	       */
	      dp_y = d - hd_min(cp9b, y, jp_y);  /* d index for state y */
	      escore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case MR_st:
	    case IR_st:
	      if (j == L) continue;
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hd_min(cp9b, y, (jp_y+1)) || (d+1) > hd_max(cp9b, y, (jp_y+1))) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
	       */
	      dp_y = d - hd_min(cp9b, y, (jp_y+1));  /* d index for state y */
	      escore = esc_vAA[y][dsq[j+1]];
	      /*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case S_st:
	    case E_st:
	    case D_st:
	      if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if (d < hd_min(cp9b, y, jp_y) || d > hd_max(cp9b, y, jp_y)) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
	       * corresponding to alpha[y][j][d] in the platonic matrix.
	       */
	      dp_y = d - hd_min(cp9b, y, jp_y);  /* d index for state y */
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->tsc[y][voffset])); 
	      break;
	    } /* end of switch(cm->sttype[y] */  
	  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	  if (beta[v][jp_v][dp_v] < IMPOSSIBLE) beta[v][jp_v][dp_v] = IMPOSSIBLE;
	} /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
      } /* end loop over jp. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] == IL_st || cm->sttype[v] == IR_st' */
    else { /* state v is not BEGL_S, BEGL_R IL nor IR (must be ML, MP, MR, D, S, B or E) */
      /* ML, MP, MR, D, S, B, E states cannot self transit, this means that all cells
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
	  ESL_DASSERT1((j >= 0 && j <= L));
	  jp_v = j - jmin[v];
	  jp_y = j - jmin[y];
	  ESL_DASSERT1((j+sdr >= jmin[y] && j+sdr <= jmax[y]));
	  
	  /* determine min d (dn) and max d (dx) that are valid for v and y and j */
	  dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y + sdr) - sd);
	  dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y + sdr) - sd);
	  dp_v = dx - hd_min(cp9b, v, jp_v);
	  dp_y = dx - hd_min(cp9b, y, jp_y + sdr);
	  i    = j-dx+1;
	  
	  /* for each emit mode, update beta[v][jp_v][dp_v] for all valid d = dp_v */
	  switch(emitmode) { 
	  case EMITPAIR:  /* MP_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
	      escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITLEFT:  /* ML_st, IL_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
	      escore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITRIGHT:  /* MR_st, IR_st */
	    escore = esc_vAA[y][dsq[j+1]]; /* not dependent on i */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITNONE:  /* D_st, S_st, E_st*/
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
	      beta[v][jp_v][dp_v] = ESL_MAX(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset]));
	    }
	    break;
	  } /* end of switch(emitmode) */
	} /* end of for j = jx; j >= jn; j-- */
      } /* end of for y = plast[v]... */
    } /* ends else entered for non-BEGL_S/BEGR_S/IL/IR states*/	
    /* we're done calculating deck v for everything but local begins */

    /* deal with local alignment end transitions v->EL (EL = deck at M.) */
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr      = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd       = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);        /* note emitmode is for state v */
      
      jn = jmin[v] - sdr;
      jx = jmax[v] - sdr;
      for (j = jn; j <= jx; j++) {
	jp_v =  j - jmin[v];
	dn   = hd_min(cp9b, v, jp_v + sdr) - sd;
	dx   = hd_max(cp9b, v, jp_v + sdr) - sd;
	i    = j-dn+1;                     /* we'll decrement this in for (d... loops inside switch below */
	dp_v = dn - hd_min(cp9b, v, jp_v + sdr);  /* we'll increment this in for (d... loops inside switch below */

	switch (emitmode) {
	case EMITPAIR:
	  for (d = dn; d <= dx; d++, dp_v++, i--) {
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v] 
								    + escore));
	  }
	  break;
	case EMITLEFT:
	  for (d = dn; d <= dx; d++, dp_v++, i--) {
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v] 
								    + escore));
	  }
	  break;
	  
	case EMITRIGHT:
	  escore = esc_vAA[v][dsq[j+1]];
	  for (d = dn; d <= dx; d++, dp_v++) {
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v]
								    + escore));
	  }
	  break;
	  
	case EMITNONE:
	  for (d = dn; d <= dx; d++, dp_v++) {
	    beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v]));
	  }
	  break;
	}
      }
    }
  } /* end loop over decks v. */

  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, EL->EL transitions. (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (j = L; j > 0; j--) { /* careful w/ boundary here */
      for (d = j-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = ESL_MAX(beta[cm->M][j][d], (beta[cm->M][j][d+1] + cm->el_selfsc));
    }
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.stdocykhbmx", "w");   cm_hb_mx_Dump(fp1, mx); fclose(fp1); */
#endif

  fail1_flag = FALSE;
  fail2_flag = FALSE;
  fail3_flag = FALSE;
  if(do_check) {
    /* Check for consistency between the Inside alpha matrix and the
     * Outside beta matrix. We assume the Inside CYK parse score
     * (optsc) is the optimal score, so for all v,j,d:
     * 
     * Jalpha[v][j][d] + Jbeta[v][j][d] <= optsc
     * 
     * Further, we know that each residue must be emitted by a state
     * in the optimal parse. So as we do the above check, we determine
     * when we're in a cell that may be involved in the optimal parse
     * (the sum of the Inside and Outside scores are equal to the
     * optimal parse score), if that cell corresponds to a left
     * emitter emitting position i, we know an emitted i has been
     * observed in an optimal parse and set optseen[i] to TRUE.
     * Likewise, if that cell corresponds to a right emitter emitting
     * position j, we update optseen[j] to TRUE. At the end of the
     * check optseen[i] should be TRUE for all i in the range
     * [1..L].
     *
     * Note that we don't ensure that all of our presumed optimal
     * cells make up a valid parse, so it is possible we could pass
     * this check even if the Inside and Outside matrices are
     * inconsistent (i.e. there's a bug in the implementation of one
     * and or the other) but that should be extremely unlikely.  If we
     * do this test many times for many different models and pass, we
     * should be confident we have consistent implementations.
     * 
     * This is an expensive check and should only be done while
     * debugging.
     *
     * Another test we could do but do not is to determine the CYK
     * parse by tracing back the CYK Inside matrix, then ensure that
     * for each cell in that parse alpha[v][j][d]+beta[v][j][d] ==
     * optsc.
     */
    ESL_ALLOC(optseen, sizeof(int) * (L+1));
    esl_vec_ISet(optseen, L+1, FALSE);
    vmax = (cm->flags & CMH_LOCAL_END) ? cm->M : cm->M-1;
    /* define bit score difference tolerance, somewhat arbitrarily: 
     * clen <= 200: tolerance is 0.001; then a function of clen: 
     * clen == 1000 tolerance is 0.005, 
     * clen == 2000, tolerance is 0.01.
     *
     * I did this b/c with tests with SSU_rRNA_eukarya I noticed
     * failures with bit score differences up to 0.004 or so.  This
     * could mean a bug, but I couldn't get any average sized model to
     * fail with a difference above 0.001, so I blamed it on
     * precision. I'm not entirely convinced it isn't a bug but
     * until I see a failure on a smaller model it seems precision
     * is the most likely explanation, right?  
     */ 
    tol = ESL_MAX(1e-3, (float) cm->clen / 200000.); 
    for(v = 0; v <= vmax; v++) { 
      jn = (v == cm->M) ? 1 : jmin[v];
      jx = (v == cm->M) ? L : jmax[v];
      for(j = jn; j <= jx; j++) { 
	jp_v = (v == cm->M) ? j : j - jmin[v];
	dn   = (v == cm->M) ? 0 : hd_min(cp9b, v, jp_v);
	dx   = (v == cm->M) ? j : hd_max(cp9b, v, jp_v);
	for(d = dn; d <= dx; d++) { 
	  dp_v = (v == cm->M) ? d : d - hd_min(cp9b, v, jp_v);
	  sc  = (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]) - alpha[0][jp_0][Lp_0];
	  if(sc > tol) { 
	    printf("Check 1 failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v], alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v], alpha[0][jp_0][Lp_0]);
	    fail1_flag = TRUE;
	  }
	  if(fabs(sc) < tol) { /* this cell is involved in a parse with the optimal score */
	    i  = j-d+1;
	    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st || (cm->sttype[v] == EL_st && d >0)) { 
	      /* i is accounted for by a parse with an optimal score */
	      optseen[i] = TRUE;
	      /*printf("\tResidue %4d possibly accounted for by Left  emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", i, Statetype(cm->sttype[v]), v, j, d);*/
	    }
	    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
	      /* j is accounted for by a parse with an optimal score */
	      optseen[j] = TRUE;
	      /*printf("\tResidue %4d possibly accounted for by Right emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", j, Statetype(cm->sttype[v]), v, j, d);*/
	    }
	  }
	}
      }
    }
    for(j = 1; j <= L; j++) { 
      if(optseen[j] == FALSE) { 
	printf("Check 2 failure: residue %d not emitted in the optimal parsetree\n", j);
	fail2_flag = TRUE;
      }	      
    }
    free(optseen);
  }
  /* Another test that we can only do if local ends are OFF */
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     *
     * Determine P(pi, S|M) / P(S|R) (probability of the sequence and most likely parse
     * tree pi given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(pi, S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Max_v [ Max_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                = P(pi, S|M) / P(S|R)
     */
    for(n = 0; n < cm->nodes; n++) {
      sc = IMPOSSIBLE;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	    dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
	    sc = ESL_MAX(sc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
	  }
	}
      }
      /*printf("checking node: %d | sc: %.6f\n", n, sc);*/
      diff = sc - alpha[0][jp_0][Lp_0];
      if(diff > 0.01 || diff < -0.01) { 
	fail3_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, sc, alpha[0][jp_0][Lp_0], diff);
      }
    }
  }

  /* Finally, calculate the optimal score, but this only works if
   * we're not in local mode:
   * 
   * If local ends are off, we know the optimal parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Max_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Max_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * *** If local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END)) { 
    sc = IMPOSSIBLE;
    v = cm->M-1;
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      assert(hd_min(cp9b, v, jp_v) == 0);
      sc = ESL_MAX(sc, (beta[v][jp_v][0]));
      /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
    }
  }
  else { /* return sc = P(S|M) / P(S|R) from Inside() */
    sc = alpha[0][jp_0][Lp_0];
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp; fp = fopen("tmp.std_ocykhbmx", "w"); cm_hb_mx_Dump(fp, mx); fclose(fp); */
#endif

  if(do_check) {
    if     (fail1_flag) ESL_FAIL(eslFAIL, errbuf, "CYK Inside/Outside HB check1 FAILED.");
    else if(fail2_flag) ESL_FAIL(eslFAIL, errbuf, "CYK Inside/Outside HB check2 FAILED.");
    else if(fail3_flag) ESL_FAIL(eslFAIL, errbuf, "CYK Inside/Outside HB check3 FAILED.");
    ESL_DPRINTF1(("#DEBUG: SUCCESS! CYK Inside/Outside HB checks PASSED.\n"));
  }

  if(!(cm->flags & CMH_LOCAL_END)) ESL_DPRINTF1(("#DEBUG: \tcm_CYKOutsideAlignHB() sc : %f\n", sc));
  else                             ESL_DPRINTF1(("#DEBUG: \tcm_CYKOutsideAlignHB() sc : %f (LOCAL mode; sc is from Inside)\n", sc));

  if(ret_sc != NULL) *ret_sc = sc;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Out of memory");
  return status; /* NEVER REACHED */
}  

/* Function: cm_OutsideAlign()
 * Date:     EPN, Mon Nov 19 07:00:37 2007
 *
 * Purpose:  Run the outside algorithm on a target sequence.
 *           Non-banded version. See cm_OutsideAlignHB() for
 *           the HMM banded version. The full target sequence
 *           1..L is aligned. 
 *
 *           Very similar to cm_CYKOutsideAlign() but calculates
 *           beta[v][j][d]: log of the summed probability of all
 *           parsetrees that emits 1..i-1 and j+1..L and pass through
 *           v at j,d (where i = j-d+1) instead of the log of the
 *           probability of the most likely (CYK) parse. This means
 *           logsum operations are used instead of max operations.
 *
 *           For debugging this function, the cm_CYKOutsideAlign() can
 *           be useful, because it has a very similar organization but
 *           is easier to debug because only the most likely parsetree
 *           is considered. cm_CYKOutsideAlign() also allows a more
 *           stringent test for the consistency of the CYKInside and
 *           CYKOutside matrices.
 *
 *           If <do_check> is TRUE (and the CM is not in local mode)
 *           we check that the outside matrix values are consistent
 *           with the inside matrix values (in ins_mx).  This check is
 *           described in comments towards the end of the function.
 *
 *           Note: renamed from FastOutsideAlign() [EPN, Wed Sep 14 06:13:53 2011].
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           L         - length of the dsq to align
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_check  - TRUE to attempt to check 
 *           mx        - the dp matrix, grown and filled here
 *           ins_mx    - the pre-filled dp matrix from the Inside run calculation (required)
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *                       ends are on (see *** comment towards end of function).
 *
 * Returns:  <eslOK> on success
 *
 * Throws:   <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *           <eslFAIL>   if <do_check>==TRUE and we fail a test
 *           In either of these cases, alignment has been aborted, ret_sc is not valid.

 */
int 
cm_OutsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, 
		CM_MX *mx, CM_MX *ins_mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float    sc;     	       /* a temporary score */
  float  **esc_vAA;            /* ptr to cm->oesc, optimized emission scores */
  float    escore;	       /* an emission score, tmp variable */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */

  /* variables used only if do_check==TRUE */
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */
  int      fail_flag = FALSE;  /* set to TRUE if do_check and we see a problem */

  /* the DP matrices */
  float ***beta  = mx->dp;     /* pointer to the Oustide DP mx */
  float ***alpha = ins_mx->dp; /* pointer to the Inside DP mx (already calc'ed and passed in) */

  /* Allocations and initializations */
  esc_vAA = cm->oesc;            /* a ptr to the optimized emission scores */

  /* grow the matrix based on the current sequence */
  if((status = cm_mx_GrowTo(cm, mx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(beta[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* now set beta[0][L][L] to 0., all parses must end there */
  beta[0][L][L] = 0.;

  /* initialize local begin cells for emitting full seq (j==L && d == L) */
  if (cm->flags & CMH_LOCAL_BEGIN) { 
    for (v = 1; v < cm->M; v++) 
      beta[v][L][L] = cm->beginsc[v];
  }

  /* Main recursion */
  for (v = 1; v < cm->M; v++) { /* start at state 1 because we set all values for ROOT_S state 0 above */
    sd  = StateDelta(cm->sttype[v]);
    sdr = StateRightDelta(cm->sttype[v]);

    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      y = cm->plast[v];	/* the parent bifurcation    */
      z = cm->cnum[y];	/* the other (right) S state */
      for(j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  for (k = 0; k <= (L-j); k++) {
	    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+k][d+k] + alpha[z][j+k][k]));
	  }
	}
      }
    } /* end of 'if (cm->stid[v] == BEGL_S */
    else if (cm->stid[v] == BEGR_S) {
      y = cm->plast[v];	  /* the parent bifurcation    */
      z = cm->cfirst[y];  /* the other (left) S state  */
      for(j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  for (k = 0; k <= (j-d); k++) {
 	    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j][d+k] + alpha[z][j-d][k]));
	  }
	}
      }
    } /* end of 'else if (cm->stid[v] == BEGR_S */
    else { /* (cm->sttype[v] != BEGL_S && cm->sttype[v] != BEGR_S */ 
      for (j = L; j >= 0; j--) {
	i = 1;
	for (d = j; d >= 0; d--, i++) {
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    sd  = StateDelta(cm->sttype[y]);
	    sdr = StateRightDelta(cm->sttype[y]);
	    switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == L || d == j) continue; /* boundary condition */
		escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;

	      case ML_st:
	      case IL_st: 
		if (d == j) continue;	/* boundary condition (note when j=0, d=0*/
		escore = esc_vAA[y][dsq[i-1]];
		beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore));
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == L) continue;
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
      } /* end loop over j. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] != BEGL_S, BEGR_S */
    /* we're done calculating deck v for everything but local begins */

    /* deal with local alignment end transitions v->EL (EL = deck at M.) */
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  i = j-d+1;
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    if (j == L || d == j) continue; /* boundary condition */
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v] 
							    + escore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == j) continue;	
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v] 
							    + escore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == L) continue;
	    escore = esc_vAA[v][dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+sdr][d+sd] + cm->endsc[v]
							    + escore));
	    break;
	  case S_st:
	  case D_st:
	  case B_st:
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
    for (j = L; j > 0; j--) { /* careful w/ boundary here */
      for (d = j-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1] + cm->el_selfsc));
    }
  }

  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     *
     * Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                = P(S|M) / P(S|R)
     */
    for(n = 0; n < cm->nodes; n++) {
      sc = IMPOSSIBLE;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = 0; j <= L; j++) {
	  for (d = 0; d <= j; d++) {
	    sc = FLogsum(sc, (alpha[v][j][d] + beta[v][j][d]));
	  }
	}
      }
      /*printf("checking node: %d | sc: %.6f\n", n, sc);*/
      diff = sc - alpha[0][L][L];
      if(diff > 0.01 || diff < -0.01) { 
	fail_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, sc, alpha[0][L][L], diff);
      }
    }
    if(! fail_flag) { 
      ESL_DPRINTF1(("#DEBUG: SUCCESS! all nodes passed error check (cm_OutsideAlign())\n"));
    }
  }

  /* Finally, calculate the optimal score, but this only works if
   * we're not in local mode:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * *** If local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END)) { 
    sc = IMPOSSIBLE;
    v = cm->M-1;
    for (j = 0; j <= L; j++) {
      sc = FLogsum(sc, (beta[v][j][0]));
      /*printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
    }
  }
  else { /* sc = P(S|M) / P(S|R) from Inside() */
    sc = alpha[0][L][L];
  }

  if(fail_flag) ESL_FAIL(eslFAIL, errbuf, "Not all nodes passed posterior check.");

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_omx", "w");   cm_mx_Dump(fp1, mx); fclose(fp1); */
#endif

  if(!(cm->flags & CMH_LOCAL_END)) ESL_DPRINTF1(("#DEBUG: \tcm_OutsideAlign() sc : %f\n", sc));
  else                             ESL_DPRINTF1(("#DEBUG: \tcm_OutsideAlign() sc : %f (LOCAL mode; sc is from Inside)\n", sc));

  if(ret_sc != NULL) *ret_sc = sc;

  return eslOK;
}  

/* Function: cm_OutsideAlignHB()
 * Date:     EPN, Thu Nov  8 18:40:05 2007
 *
 * Purpose:  Run the outside algorithm on a target sequence.
 *           HMM banded version. See cm_OutsideAlign() for
 *           the non-banded version. The full target sequence
 *           1..L is aligned. 
 *
 *           Very similar to cm_CYKOutsideAlignHB() but calculates
 *           beta[v][j][d]: log of the summed probability of all
 *           parsetrees that emits 1..i-1 and j+1..L and pass through
 *           v at j,d (where i = j-d+1) instead of the log of the
 *           probability of the most likely (CYK) parse. This means
 *           logsum operations are used instead of max operations.
 *
 *           For debugging this function, the cm_CYKOutsideAlign() can
 *           be useful, because it has a very similar organization but
 *           is easier to debug because only the most likely parsetree
 *           is considered. cm_CYKOutsideAlign() also allows a more
 *           stringent test for the consistency of the CYKInside and
 *           CYKOutside matrices.
  *
 *           The DP recursion has been 'optimized' for all state types
 *           except IL, IR, BEGL_S, BEGR_S. The main optimization
 *           is a change in nesting order of the for loops:
 *           optimized order:     for v { for y { for j { for d {}}}}
 *           non-optimized order: for v { for j { for d { for y {}}}}
 * 
 *           ILs and IRs are not optimized because they can self
 *           transit so mx[v][j][d] must be fully calc'ed before
 *           mx[v][j][d+1] can be calced. BEGL_S and BEGR_S are not
 *           optimized b/c they require searching for optimal d and k,
 *           which complicates the enforcement of the bands and makes
 *           this optimization strategy impossible.
 *
 *           If <do_check> is TRUE (and the CM is not in local mode) 
 *           we check that the outside matrix values are consistent
 *           with the inside matrix values (in ins_mx).  This check is
 *           described in comments towards the end of the function.
 *
 *           Note: renamed from FastOutsideAlignHB() [EPN, Wed Sep 14 06:13:53 2011].
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           L         - length of the dsq to align
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_check  - TRUE to attempt to check 
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ins_mx    - the dp matrix from the Inside run calculation (required)
 *           ret_sc    - RETURN: log P(S|M)/P(S|R), as a bit score, this is from ins_mx IF local
 *                       ends are on (see *** comment towards end of function).
 *
 * Returns:  <eslOK> on success
 *
 * Throws:   <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *           <eslFAIL>   if <do_check>==TRUE and we fail a test
 *           In either of these cases, alignment has been aborted, ret_sc is not valid.
 */
int
cm_OutsideAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, 
		  CM_HB_MX *mx, CM_HB_MX *ins_mx, float *ret_sc)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float  **esc_vAA;            /* ptr to cm->oesc, optimized emission scores */
  float    sc;		       /* a temporary score */
  float    escore;	       /* an emission score, tmp variable */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */

  /* variables used only if do_check */
  int      fail_flag = FALSE;  /* set to TRUE if do_check and we see a problem */
  int      n;                  /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states;   /* temp variable used only if do_check = TRUE */
  float    diff;               /* temp variable used only if do_check = TRUE */

  /* band related variables */
  int      dp_v;               /* d index for state v in alpha w/mem eff bands */
  int      dp_y;               /* d index for state y in alpha w/mem eff bands */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Lp;                 /* L index also changes depending on state */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      kmin, kmax;         /* temporary minimum/maximum allowed k */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      jp_0;               /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;               /* L offset in ROOT_S's (v==0) d band */
  int     *d_max_written   = NULL; /* [0..L]: max d of write RANGE written to beta[cm->M][j][d]; -1 if no write */
  int     *d_max_nonimpos  = NULL; /* [0..L]: max d where written value is actually non-IMPOSSIBLE; -1 if none */

  /* the DP matrices */
  float ***beta  = mx->dp;     /* pointer to the Oustide DP mx */
  float ***alpha = ins_mx->dp; /* pointer to the Inside DP mx (already calc'ed and passed in) */

  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cm->cp9b->jmin;
  int     *jmax  = cm->cp9b->jmax;
  CP9Bands_t *cp9b = cm->cp9b;  /* brief 26_0430-157: needed by hd_min()/hd_max() */
  int    **hdmin = cm->cp9b->hdmin;
  int    **hdmax = cm->cp9b->hdmax;

  /* Allocations and initializations */
  esc_vAA = cm->oesc;            /* a ptr to the optimized emission scores */

  /* grow the matrix based on the current sequence and bands */
  if((status = cm_hb_mx_GrowTo(cm, mx, errbuf, cm->cp9b, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(beta[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* allocate and initialize d_max_written/d_max_nonimpos for the sparse EL self-transition */
  if (cm->flags & CMH_LOCAL_END) {
    ESL_ALLOC(d_max_written,  sizeof(int) * (L+1));
    ESL_ALLOC(d_max_nonimpos, sizeof(int) * (L+1));
    esl_vec_ISet(d_max_written,  L+1, -1);
    esl_vec_ISet(d_max_nonimpos, L+1, -1);
  }

  /* ensure a full alignment to ROOT_S (v==0) is allowed by the bands */
  if (jmin[0] > L || jmax[0] < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, jmin[0], jmax[0]);
  jp_0 = L - jmin[0];
  if (hd_min(cp9b, 0, jp_0) > L || hd_max(cp9b, 0, jp_0) < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, hd_min(cp9b, 0, jp_0), hd_max(cp9b, 0, jp_0));
  Lp_0 = L - hd_min(cp9b, 0, jp_0);
  /* set the offset banded cell corresponding to beta[0][L][L] to 0., all parses must end there */
  beta[0][jp_0][Lp_0] = 0.;

  /* If we can do a local begin into v, overwrite IMPOSSIBLE with the local begin score. */
  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (v = 1; v < cm->M; v++) {
      if(NOT_IMPOSSIBLE(cm->beginsc[v])) {
	if((L >= jmin[v]) && (L <= jmax[v])) {
	  jp_v = L - jmin[v];
	  if((L >= hd_min(cp9b, v, jp_v)) && L <= hd_max(cp9b, v, jp_v)) {
	    Lp = L - hd_min(cp9b, v, jp_v);
	    beta[v][jp_v][Lp] = cm->beginsc[v];
	  }
	}
      }
    }
  }
  /* done allocation/initialization */

  /* Recursion: main loop down through the decks */
  for (v = 1; v < cm->M; v++) { /* start at state 1 because we set all values for ROOT_S state 0 above */
    if (cm->stid[v] == BEGL_S) { /* BEGL_S */
      y = cm->plast[v];	/* the parent bifurcation    */
      z = cm->cnum[y];	/* the other (right) S state */
      for (j = jmax[v]; j >= jmin[v]; j--) {
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;
	for (d = hd_max(cp9b, v, jp_v); d >= hd_min(cp9b, v, jp_v); d--) {
	  dp_v = d - hd_min(cp9b, v, jp_v);
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
	    if(k < (hd_min(cp9b, y, jp_y+k) - d) || k > (hd_max(cp9b, y, jp_y+k) - d)) continue; 
	    /* above line continues if inequality 5 or 6 is violated */
	    if(k < (hd_min(cp9b, z, jp_z+k))     || k > (hd_max(cp9b, z, jp_z+k)))     continue; 
	    /* above line continues if inequality 7 or 8 is violated */
		  
	    /* if we get here for current k, all 8 inequalities have been satisified 
	     * so we know the cells corresponding to the platonic 
	     * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
	     * alpha[z][j+k][k] are all within the bands. These
	     * cells correspond to beta[v][jp_v][dp_v], 
	     * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
	     * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
	     */
	    kp_z = k-hd_min(cp9b, z, jp_z+k);
	    dp_y = d-hd_min(cp9b, y, jp_y+k);
	    beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] + alpha[z][jp_z+k][kp_z]));
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
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	i = j-d+1;

	dn = ESL_MAX(hd_min(cp9b, v, jp_v), j-jmax[z]);
	dx = ESL_MIN(hd_max(cp9b, v, jp_v), jp_z);
	/* above makes sure that j,d are valid for state z: (jmin[z] + d) >= j >= (jmax[z] + d) */
	for (d = dx; d >= dn; d--) {
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
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
	  kmin = ESL_MAX((hd_min(cp9b, y, jp_y)-d), (hd_min(cp9b, z, jp_z-d)));
	  kmax = ESL_MIN((hd_max(cp9b, y, jp_y)-d), (hd_max(cp9b, z, jp_z-d)));
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
	    kp_z = k-hd_min(cp9b, z, jp_z-d);
	    dp_y = d-hd_min(cp9b, y, jp_y);
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
	ESL_DASSERT1((j >= 0 && j <= L));
	jp_v = j - jmin[v];
	for (d = hd_max(cp9b, v, jp_v); d >= hd_min(cp9b, v, jp_v); d--) {
	  i = j-d+1;
	  dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
	  
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
	    
	    /* Note: this looks like it can be optimized, I tried but my 'optimization' slowed the code, so I reverted [EPN] */
	    switch(cm->sttype[y]) {
	    case MP_st: 
	      if (j == L || d == j) continue; /* boundary condition */
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+2) < hd_min(cp9b, y, (jp_y+1)) || (d+2) > hd_max(cp9b, y, (jp_y+1))) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
		   */
	      dp_y = d - hd_min(cp9b, y, jp_y+1);  /* d index for state y */
	      escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case ML_st:
	    case IL_st: 
	      if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
	      if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hd_min(cp9b, y, jp_y) || (d+1) > hd_max(cp9b, y, jp_y)) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j][d+1] in the platonic matrix.
	       */
	      dp_y = d - hd_min(cp9b, y, jp_y);  /* d index for state y */
	      escore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
								  + cm->tsc[y][voffset] + escore));
	      break;
	      
	    case MR_st:
	    case IR_st:
	      if (j == L) continue;
	      if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
	      jp_y = j - jmin[y];
	      if ((d+1) < hd_min(cp9b, y, (jp_y+1)) || (d+1) > hd_max(cp9b, y, (jp_y+1))) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
	       * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
	       */
	      dp_y = d - hd_min(cp9b, y, (jp_y+1));  /* d index for state y */
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
	      if (d < hd_min(cp9b, y, jp_y) || d > hd_max(cp9b, y, jp_y)) continue; /* enforces d is valid for state y */
	      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
	       * corresponding to alpha[y][j][d] in the platonic matrix.
	       */
	      dp_y = d - hd_min(cp9b, y, jp_y);  /* d index for state y */
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->tsc[y][voffset])); 
	      break;
	    } /* end of switch(cm->sttype[y] */  
	  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	  if (beta[v][jp_v][dp_v] < IMPOSSIBLE) beta[v][jp_v][dp_v] = IMPOSSIBLE;
	} /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
      } /* end loop over jp. We know beta for this whole state */
    } /* end of 'else if cm->sttype[v] == IL_st || cm->sttype[v] == IR_st' */
    else { /* state v is not BEGL_S, BEGL_R IL nor IR (must be ML, MP, MR, D, S, B or E) */
      /* ML, MP, MR, D, S, B, E states cannot self transit, this means that all cells
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
	  ESL_DASSERT1((j >= 0 && j <= L));
	  jp_v = j - jmin[v];
	  jp_y = j - jmin[y];
	  ESL_DASSERT1((j+sdr >= jmin[y] && j+sdr <= jmax[y]));
	  
	  /* determine min d (dn) and max d (dx) that are valid for v and y and j */
	  dn = ESL_MAX(hd_min(cp9b, v, jp_v), hd_min(cp9b, y, jp_y + sdr) - sd);
	  dx = ESL_MIN(hd_max(cp9b, v, jp_v), hd_max(cp9b, y, jp_y + sdr) - sd);
	  dp_v = dx - hd_min(cp9b, v, jp_v);
	  dp_y = dx - hd_min(cp9b, y, jp_y + sdr);
	  i    = j-dx+1;
	  
	  /* for each emit mode, update beta[v][jp_v][dp_v] for all valid d = dp_v */
	  switch(emitmode) { 
	  case EMITPAIR:  /* MP_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
	      escore = esc_vAA[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITLEFT:  /* ML_st, IL_st */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--, i++) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
	      escore = esc_vAA[y][dsq[i-1]];
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITRIGHT:  /* MR_st, IR_st */
	    escore = esc_vAA[y][dsq[j+1]]; /* not dependent on i */
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
	      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y + sdr][dp_y + sd] 
								  + cm->tsc[y][voffset] + escore));
	    }
	    break;
	  case EMITNONE:  /* D_st, S_st, E_st*/
	    for (d = dx; d >= dn; d--, dp_v--, dp_y--) { 
	      ESL_DASSERT1((  d       >= hd_min(cp9b, v, jp_v)        &&   d       <= hd_max(cp9b, v, jp_v)));
	      ESL_DASSERT1((((d + sd) >= hd_min(cp9b, y, jp_y + sdr)) && ((d + sd) <= hd_max(cp9b, y, jp_y + sdr))));
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
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr      = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd       = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);        /* note emitmode is for state v */
      
      jn = jmin[v] - sdr;
      jx = jmax[v] - sdr;
      for (j = jn; j <= jx; j++) {
	jp_v =  j - jmin[v];
	dn   = hd_min(cp9b, v, jp_v + sdr) - sd;
	dx   = hd_max(cp9b, v, jp_v + sdr) - sd;
	i    = j-dn+1;                     /* we'll decrement this in for (d... loops inside switch below */
	dp_v = dn - hd_min(cp9b, v, jp_v + sdr);  /* we'll increment this in for (d... loops inside switch below */

	switch (emitmode) {
	case EMITPAIR:
	  for (d = dn; d <= dx; d++, dp_v++, i--) {
	    escore = esc_vAA[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v] 
								    + escore));
	  }
	  break;
	case EMITLEFT:
	  for (d = dn; d <= dx; d++, dp_v++, i--) {
	    escore = esc_vAA[v][dsq[i-1]];
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v] 
								    + escore));
	  }
	  break;
	  
	case EMITRIGHT:
	  escore = esc_vAA[v][dsq[j+1]];
	  for (d = dn; d <= dx; d++, dp_v++) {
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v]
								    + escore));
	  }
	  break;
	  
	case EMITNONE:
	  for (d = dn; d <= dx; d++, dp_v++) {
	    beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+sdr][dp_v+sd] + cm->endsc[v]));
	  }
	  break;
	}
	/* update d_max_written for the sparse self-transition */
	if (dx >= dn && dx > d_max_written[j]) d_max_written[j] = dx;
	/* update d_max_nonimpos: highest d with actual non-IMPOSSIBLE write for this j */
	if (dx >= dn) {
	  int d_tmp;
	  int d_lo = (dn >= 0) ? dn : 0;
	  for (d_tmp = dx; d_tmp >= d_lo; d_tmp--) {
	    if (NOT_IMPOSSIBLE(beta[cm->M][j][d_tmp])) {
	      if (d_tmp > d_max_nonimpos[j]) d_max_nonimpos[j] = d_tmp;
	      break;
	    }
	  }
	}
      }
    }
  } /* end loop over decks v. */

  /* Deal with last step needed for local alignment
   * w.r.t. ends: left-emitting, EL->EL transitions. (EL = deck at M.)
   * Sparse optimization: skip j values where d_max_nonimpos[j] < 0 (no non-IMPOSSIBLE
   * EL writes at j; all self-transitions would be IMPOSSIBLE → no-op). For j values
   * with non-IMPOSSIBLE writes, start the d-sweep at d_max_nonimpos[j]-1 instead of
   * j-1: cells at d > d_max_nonimpos[j] are IMPOSSIBLE and propagate nothing.
   */
  if (cm->flags & CMH_LOCAL_END) {
    if (getenv("EL_DMAX_DIAG")) {
      int n_written_j = 0, d_max_total = 0, d_max_max = 0;
      int n_nonimpos_j = 0, d_ni_total = 0, d_ni_max = 0;
      for (j = 0; j <= L; j++) {
        if (d_max_written[j] >= 0) {
          n_written_j++;
          d_max_total += d_max_written[j];
          if (d_max_written[j] > d_max_max) d_max_max = d_max_written[j];
        }
        if (d_max_nonimpos[j] >= 0) {
          n_nonimpos_j++;
          d_ni_total += d_max_nonimpos[j];
          if (d_max_nonimpos[j] > d_ni_max) d_ni_max = d_max_nonimpos[j];
        }
      }
      fprintf(stderr, "#EL_DMAX    M=%d L=%d n_written_j=%d avg_dmax_written=%.1f max_dmax_written=%d\n",
              cm->M, L, n_written_j, n_written_j>0 ? (double)d_max_total/n_written_j : 0.0, d_max_max);
      fprintf(stderr, "#EL_NONIMPOS M=%d L=%d n_nonimpos_j=%d avg_dmax_nonimpos=%.1f max_dmax_nonimpos=%d\n",
              cm->M, L, n_nonimpos_j, n_nonimpos_j>0 ? (double)d_ni_total/n_nonimpos_j : 0.0, d_ni_max);
      fflush(stderr);
    }
    /* DEBUG: verify d_max_nonimpos is correct before using it */
    if (getenv("EL_DMAX_VERIFY")) {
      for (j = 0; j <= L; j++) {
        int true_dmax = -1;
        for (d = j; d >= 0; d--) {
          if (NOT_IMPOSSIBLE(beta[cm->M][j][d])) { true_dmax = d; break; }
        }
        if (true_dmax != d_max_nonimpos[j]) {
          fprintf(stderr, "BUG: j=%d d_max_nonimpos=%d true_dmax=%d\n",
                  j, d_max_nonimpos[j], true_dmax);
          fflush(stderr);
        }
      }
    }
    /* run the REFERENCE sweep (d_max_written) and save beta copy for comparison */
    if (getenv("EL_SWEEP_COMPARE")) {
      /* reference sweep */
      for (j = L; j > 0; j--) {
        if (d_max_written[j] < 0) continue;
        for (d = d_max_written[j] - 1; d >= 0; d--)
          beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1] + cm->el_selfsc));
      }
    } else {
      for (j = L; j > 0; j--) {
        if (d_max_written[j] < 0) continue; /* no EL writes for this j; skip */
        if (d_max_nonimpos[j] < 0) {
          /* all writes were IMPOSSIBLE; self-transition has no effect; skip */
          continue;
        }
        /* sweep from d_max_nonimpos[j]-1 downward: the step at d_max_nonimpos[j] is always a
         * no-op (beta[j][d_max_nonimpos[j]+1] is IMPOSSIBLE), so skip it and start one below.
         * This also avoids a potential out-of-bounds access when d_max_nonimpos[j] == j. */
        for (d = d_max_nonimpos[j] - 1; d >= 0; d--)
          beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1] + cm->el_selfsc));
      }
    }
    free(d_max_written);
    d_max_written = NULL;
    free(d_max_nonimpos);
    d_max_nonimpos = NULL;
  }

  if(do_check && (!(cm->flags & CMH_LOCAL_END))) {
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     *    
     * Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
     * using both the Outside (beta) and Inside (alpha) matrices,
     * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
     * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
     *                                                    = P(S|M) / P(S|R)
     */
    
    for(n = 0; n < cm->nodes; n++) {
      sc = IMPOSSIBLE;
      num_split_states = SplitStatesInNode(cm->ndtype[n]);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	    dp_v = d - hd_min(cp9b, v, jp_v);  /* d index for state v in alpha w/mem eff bands */
	    sc = FLogsum(sc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
	    /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
	      printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	  }
	}
      }
      /*printf("checking node: %d | sc: %.6f\n", n, sc);*/
      diff = sc - alpha[0][jp_0][Lp_0];
      if(diff > 0.01 || diff < -0.01) { 
	fail_flag = TRUE;
	printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
	       n, sc, alpha[0][jp_0][Lp_0], diff);
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
    sc = IMPOSSIBLE;
    v = cm->M-1;
    for (j = jmin[v]; j <= jmax[v]; j++) {
      jp_v = j - jmin[v];
      assert(hd_min(cp9b, v, jp_v) == 0);
      sc = FLogsum(sc, (beta[v][jp_v][0]));
      /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
    }
  }
  else { /* return_sc = P(S|M) / P(S|R) from Inside() */
    sc = alpha[0][jp_0][Lp_0];
  }

  if(fail_flag) ESL_FAIL(eslFAIL, errbuf, "Not all nodes passed posterior check.");

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_ohbmx", "w");   cm_hb_mx_Dump(fp1, mx); fclose(fp1); */
#endif


  if(!(cm->flags & CMH_LOCAL_END)) ESL_DPRINTF1(("#DEBUG: \tcm_OutsideAlignHB() sc : %f\n", sc));
  else                             ESL_DPRINTF1(("#DEBUG: \tcm_OutsideAlignHB() sc : %f (LOCAL mode; sc is from Inside)\n", sc));

  if (ret_sc != NULL) *ret_sc = sc;
  return eslOK;

 ERROR:
  if (d_max_written)  free(d_max_written);
  if (d_max_nonimpos) free(d_max_nonimpos);
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_Posterior()
 * Date:     EPN, Mon Nov 19 09:02:12 2007
 * Note:     based on Ian Holmes' P7EmitterPosterior() from HMMER's 2.x postprob.c
 *           Renamed from CMPosterior() [EPN, Wed Sep 14 06:15:22 2011].
 *
 * Purpose: Combines non-banded Inside and Outside matrices into a
 *           posterior probability matrix. The value in post[v][j][d]
 *           is the log of the posterior probability of a parse
 *           subtree rooted at v emitting the subsequence i..j
 *           (i=j-d+1).  The caller must provide a <post> float
 *           matrix, but this matrix may be the same matrix as that
 *           provided as Outside <out_mx>, (overwriting it will not
 *           compromise the algorithm). Posteriors are calculated
 *           for the full sequence 1..L.
 *
 *           
 * Args:     cm         - the model
 *           errbuf     - char buffer for reporting errors
 *           L          - length of the dsq to align
 *           size_limit - max number of Mb for DP matrix
 *           ins_mx     - pre-calculated Inside matrix 
 *           out_mx     - pre-calculated Outside matrix
 *           post_mx    - pre-allocated matrix for Posteriors 
 *
 * Returns:  <eslOK>     on success.
 * Throws:   <eslERANGE> if required DP matrix size exceeds <size_limit>, in 
 *                       this case, post_mx is not filled.
 */
int
cm_Posterior(CM_t *cm, char *errbuf, int L, float size_limit, CM_MX *ins_mx, CM_MX *out_mx, CM_MX *post_mx)
{
  int   status;
  int   v, j, d; /* state, position, subseq length */
  int   vmax;    /* cm->M if local ends on, else cm->M-1 */
  float sc;      /* optimal score, from Inside matrix */
  
  /* the DP matrices */
  float ***alpha = ins_mx->dp; /* pointer to the alpha DP matrix */
  float ***beta  = out_mx->dp; /* pointer to the beta DP matrix */
  float ***post  = post_mx->dp; /* pointer to the post DP matrix */

  /* grow our post matrix, but only if isn't also our out_mx in which
   * case we know we're already big enought (also in that case we
   * don't want to call GrowTo b/c it can potentially free the DP
   * matrix memory and reallocate it, which would be bad b/c we 
   * need the out_mx!) 
   */
  if(post_mx != out_mx) { 
    if((status = cm_mx_GrowTo(cm, post_mx, errbuf, L, size_limit)) != eslOK) return status;
  }

  sc = ins_mx->dp[0][L][L];

  /* If local ends are on, start with the EL state (cm->M), otherwise
   * its not a valid deck. 
   */
  vmax = (cm->flags & CMH_LOCAL_END) ? cm->M : cm->M-1;
  for (v = vmax; v >= 0; v--) {
    for (j = 0; j <= L; j++) {
      for (d = 0; d <= j; d++) {
	post[v][j][d] = alpha[v][j][d] + beta[v][j][d] - sc;
      }
    }
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_pmx", "w");   cm_mx_Dump(fp1, post_mx); fclose(fp1); */
#endif

  return eslOK;
}

/* Function: cm_PosteriorHB()
 * Date:     EPN 05.27.06 
 * Note:     based on Ian Holmes' P7EmitterPosterior() from HMMER's 2.x postprob.c
 *           Renamed from CMPosteriorHB() [EPN, Wed Sep 14 06:14:48 2011].
 *
 * Purpose: Combines HMM banded Inside and Outside matrices into a
 *           posterior probability matrix. Any cells outside of HMM
 *           bands do not exist in memory. The value in
 *           post[v][jp_v][dp_v] is the log of the posterior
 *           probability of a parse subtree rooted at v emitting the
 *           subsequence i..j (i=j-d+1). Where j = jp_v + jmin[v], and
 *           d = dp_v + hdmin[v][jp_v]. The caller must provide a
 *           <post> CM_HB_MX matrix, but this matrix may be the same
 *           matrix as that provided as Outside <out_mx>, (overwriting
 *           it will not compromise the algorithm). Posteriors are
 *           calculated for the full sequence 1..L.
 *           
 * Args:     cm         - the model
 *           errbuf     - char buffer for reporting errors
 *           L          - length of the dsq to align
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           ins_mx     - pre-calculated Inside matrix 
 *           out_mx     - pre-calculated Outside matrix
 *           post_mx    - pre-allocated matrix for Posteriors 
 *
 * Returns:  <eslOK>     on success.
 * Throws:   <eslERANGE> if required DP matrix size exceeds <size_limit>
 *           <eslEINVAL> if the full sequence is not within the bands for state 0
 *           In either case the post_mx is not filled
 */
int
cm_PosteriorHB(CM_t *cm, char *errbuf, int L, float size_limit, CM_HB_MX *ins_mx, CM_HB_MX *out_mx, CM_HB_MX *post_mx)
{
  int      status;
  int      v, j, d; /* state, position, position, subseq length */
  float    sc;      /* total score, the log probability of the current seq  */
  int      jp_v;    /* j index for state v in alpha/beta with HMM bands */
  int      dp_v;    /* d index for state v in alpha/beta with HMM bands */
  int      jp_0;        /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;        /* L offset in ROOT_S's (v==0) d band */

  /* the DP matrices */
  float ***alpha = ins_mx->dp; /* pointer to the alpha DP matrix */
  float ***beta  = out_mx->dp; /* pointer to the beta DP matrix */
  float ***post  = post_mx->dp; /* pointer to the post DP matrix */

  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cm->cp9b->jmin;  
  int     *jmax  = cm->cp9b->jmax;
  CP9Bands_t *cp9b = cm->cp9b;  /* brief 26_0430-157: needed by hd_min()/hd_max() */
  int    **hdmin = cm->cp9b->hdmin;
  int    **hdmax = cm->cp9b->hdmax;

  /* ensure a full alignment to ROOT_S (v==0) is allowed by the bands */
  if (cm->cp9b->jmin[0] > L || cm->cp9b->jmax[0] < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, cm->cp9b->jmin[0], cm->cp9b->jmax[0]);
  jp_0 = L - jmin[0];
  if (hd_min(cm->cp9b, 0, jp_0) > L || hd_max(cm->cp9b, 0, jp_0) < L) 
    ESL_FAIL(eslEINVAL, errbuf, "cm_CYKInsideAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, hd_min(cm->cp9b, 0, jp_0), hd_max(cm->cp9b, 0, jp_0));
  Lp_0 = L - hd_min(cp9b, 0, jp_0);

  sc = alpha[0][jp_0][Lp_0];

  /* grow our post matrix, but only if isn't also our out_mx in which
   * case we know we're already big enought (also in that case we
   * don't want to call GrowTo b/c it can potentially free the DP
   * matrix memory and reallocate it, which would be bad b/c we 
   * need the out_mx!) 
   */
  if(post_mx != out_mx) { 
    if((status = cm_hb_mx_GrowTo(cm, post_mx, errbuf, cm->cp9b, L, size_limit)) != eslOK) return status; 
  }

  /* If local ends are on, fill the EL state (cm->M) posterior deck.
   * EL optimization: alpha[cm->M][j][d] == el_scA[d] always, so use
   * el_scA[d] directly instead of reading the Inside EL deck
   * (which was not filled, per Changes 1+2).
   */
  if (cm->flags & CMH_LOCAL_END) {
    float *el_scA_post;
    ESL_ALLOC(el_scA_post, sizeof(float) * (L+1));
    for (d = 0; d <= L; d++) el_scA_post[d] = cm->el_selfsc * d;
    for (j = 0; j <= L; j++) {
      if (!NOT_IMPOSSIBLE(beta[cm->M][j][0])) continue; /* j-skip: all beta IMPOSSIBLE, post stays IMPOSSIBLE */
      for (d = 0; d <= j; d++) {
	if (!NOT_IMPOSSIBLE(beta[cm->M][j][d])) break; /* contiguous strip: IMPOSSIBLE means done for this j */
	post[cm->M][j][d] = el_scA_post[d] + beta[cm->M][j][d] - sc;
      }
    }
    free(el_scA_post);
  }

  for (v = (cm->M-1); v >= 0; v--) {
    for (j = jmin[v]; j <= jmax[v]; j++) {
      ESL_DASSERT1((j >= 0 && j <= L));
      jp_v = j - jmin[v];
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	dp_v = d - hd_min(cp9b, v, jp_v);
	post[v][jp_v][dp_v] = alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v] - sc;
	/*printf("v: %3d | jp_v: %3d | dp_v: %3d | alpha: %5.2f | beta: %5.2f\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
      }
    }
  }
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_EmitterPosterior()
 * Date:     EPN, Fri Sep 30 13:53:57 2011
 *
 * Purpose: Given a posterior probability cube, where the value in
 *           post[v][j][d] is the log of the posterior probability of
 *           a parse subtree rooted at v emitting the subsequence i..j
 *           (i=j-d+1), fill a CM_EMIT_MX <emit_mx> with two 2-dimensional
 *           matrices with values:
 *
 *           emit_mx->l_pp[v][i]: log of the posterior probability that
 *           state v emitted residue i leftwise either at (if a match
 *           state) or *after* (if an insert state) the left consensus
 *           position modeled by state v's node.
 *
 *           emit_mx->r_pp[v][i]: log of the posterior probability that
 *           state v emitted residue i rightwise either at (if a match
 *           state) or *before* (if an insert state) the right
 *           consensus position modeled by state v's node.
 *
 *           l_pp[v] is NULL for states that do not emit leftwise 
 *           r_pp[v] is NULL for states that do not emit rightwise
 *
 *          This is done in 3 steps:
 *          1. Fill l_pp[v][i] and r_pp[v][i] with the posterior
 *             probability that state v emitted residue i either
 *             leftwise (l_pp) or rightwise (r_pp).
 *
 *          2. Normalize l_pp and r_pp so that probability that
 *             each residue was emitted by any state is exactly
 *             1.0.
 *
 *          3. Combine l_pp values for MATP_MP (v) and MATP_ML (y=v+1)
 *             states in the same node so they give the value defined
 *             above (i.e. l_pp[v] == l_pp[y] = the PP that either v
 *             or y emitted residue i) instead of l_pp[v] = PP that v
 *             emitted i, and l_pp[y] = PP that y emitted i.  And
 *             combine r_pp values for MATP_MP (v) and MATP_MR (y=v+2)
 *             states in an analogous way.
 *             
 *          If <do_check> we check to make sure the summed probability
 *          of any residue is > 0.98 and < 1.02 prior the step 2
 *          normalization, and throw eslFAIL if not. 
 *          
 *          Note: A failure of this test does not necessarily mean a
 *          bug in the code, because this check is known to fail for
 *          some cases with parsetrees that contain inserts of 100s of
 *          residues from the same IL or IR state (that utilize 100s
 *          of IL->IL or IR->IR self transitions). These cases were
 *          looked at in detail to determine if they were due to a bug
 *          in the DP code. This was logged in
 *          ~nawrockie/notebook/8_1016_inf-1rc3_bug_alignment/00LOG.
 *          The conclusion was that the failure of the posterior check
 *          is due completely to lack of precision in the float scores
 *          (not just in the logsum look-up table but also with using
 *          real log() and exp() calls). If this function returns an
 *          error, please check to see if the parsetree has a large
 *          insertion in it, if so you can expect probabilities up to
 *          1.03 due solely to this precision issue. See the notebook
 *          00LOG for more, included a check I performed to change the
 *          relevant IL->IL transition probability by very small
 *          values (~0.0001) and you can observe the posteriors change
 *          dramatically which demonstrates that precision of floats
 *          is the culprit.  (EPN, Sun Oct 26 14:54:31 2008
 *          (originally added to cm_Posterior() function 'Purpose'
 *          function which no longer exists, having been replaced by
 *          this function.)
 *
 * 
 * Args:     cm         - the model
 *           errbuf     - for error messages
 *           L          - length of the sequence
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           post       - pre-filled posterior cube
 *           emit_mx     - pre-allocated emit matrix, grown and filled-in here
 *           do_check   - if TRUE, return eslEFAIL if summed prob of any residue 
 *                        (before normalization) is < 0.98 or > 1.02.
 * 
 * Returns:  <eslOK>     on success.
 * Throws:   <eslERANGE> if required DP matrix size exceeds <size_limit>
 *           <eslFAIL>   if (do_check) and any residue check fails
 *           <eslEMEM>   if we run out of memory. 
 *           If !eslOK the l_pp and r_pp values are invalid.
 */
int 
cm_EmitterPosterior(CM_t *cm, char *errbuf, int L, float size_limit, CM_MX *post, CM_EMIT_MX *emit_mx, int do_check)
{
  int    status;
  int    v, j, d; /* state, position, subseq length */
  int    i;       /* sequence position */
  int    sd;      /* StateDelta(v) */
  
  /* grow the emit matrices based on the current sequence */
  if((status = cm_emit_mx_GrowTo(cm, emit_mx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the emit matrices to IMPOSSIBLE */
  esl_vec_FSet(emit_mx->l_pp_mem, emit_mx->l_ncells_valid, IMPOSSIBLE);
  esl_vec_FSet(emit_mx->r_pp_mem, emit_mx->r_ncells_valid, IMPOSSIBLE);

  /* Step 1. Fill l_pp[v][i] and r_pp[v][i] with the posterior
   *         probability that state v emitted residue i either
   *         leftwise (l_pp) or rightwise (r_pp).
   */
  for(v = 0; v < cm->M; v++) { 
    sd = StateDelta(cm->sttype[v]);
    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
      for(j = 1; j <= L; j++) { 
	i = j-sd+1;
	for(d = sd; d <= j; d++, i--) { 
	  emit_mx->l_pp[v][i] = FLogsum(emit_mx->l_pp[v][i], post->dp[v][j][d]);
	}
      }
    }
    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
      for(j = 1; j <= L; j++) { 
	sd = StateDelta(cm->sttype[v]);
	for(d = sd; d <= j; d++) { 
	  emit_mx->r_pp[v][j] = FLogsum(emit_mx->r_pp[v][j], post->dp[v][j][d]);
	}
      }
    }
  }
  /* factor in contribution of local ends, the EL state may have emitted this residue. */
  if (cm->flags & CMH_LOCAL_END) {
    for (j = 1; j <= L; j++) { 
      i = j;
      for (d = 1; d <= j; d++, i--) { /* note: d >= 1, b/c EL emits 1 residue */
	emit_mx->l_pp[cm->M][i] = FLogsum(emit_mx->l_pp[cm->M][i], post->dp[cm->M][j][d]);
      }
    }
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_unnorm_emitmx",  "w"); cm_emit_mx_Dump(fp1, cm, emit_mx); fclose(fp1); */
#endif

  /* Step 2. Normalize l_pp and r_pp so that probability that
   *         each residue was emitted by any state is exactly
   *         1.0.
   */
  esl_vec_FSet(emit_mx->sum, (L+1), IMPOSSIBLE);
  for(v = 0; v <= cm->M; v++) { 
    if(emit_mx->l_pp[v] != NULL) {
      for(i = 1; i <= L; i++) { 
	emit_mx->sum[i] = FLogsum(emit_mx->sum[i], emit_mx->l_pp[v][i]);
      }
    }
    if(emit_mx->r_pp[v] != NULL) {
      for(j = 1; j <= L; j++) { 
	emit_mx->sum[j] = FLogsum(emit_mx->sum[j], emit_mx->r_pp[v][j]);
      }
    }
  }
  /* perform the check, if nec */
  if(do_check) { 
    for(i = 1; i <= L; i++) { 
      if((sreEXP2(emit_mx->sum[i]) < 0.98) || (sreEXP2(emit_mx->sum[i]) > 1.02)) { 
	ESL_FAIL(eslFAIL, errbuf, "residue %d has summed prob of %5.4f (2^%5.4f).\nMay not be a DP coding bug, see 'Note:' on precision in cm_EmitterPosterior().\n", i, (sreEXP2(emit_mx->sum[i])), emit_mx->sum[i]);
      }
      printf("i: %d | total: %10.4f\n", i, (sreEXP2(emit_mx->sum[i])));
    }
    ESL_DPRINTF1(("#DEBUG: cm_EmitterPosterior() check passed, all residues have summed probability of emission of between 0.98 and 1.02.\n"));
  }  

  /* normalize, using the sum vector */
  for(v = 0; v <= cm->M; v++) { 
    if(emit_mx->l_pp[v] != NULL) {
      for(i = 1; i <= L; i++) { 
	emit_mx->l_pp[v][i] -= emit_mx->sum[i];
      }
    }
    if(emit_mx->r_pp[v] != NULL) {
      for(j = 1; j <= L; j++) { 
	emit_mx->r_pp[v][j] -= emit_mx->sum[j];
      }
    }
  }

  /* Step 3. Combine l_pp values for MATP_MP (v) and MATP_ML (y=v+1)
   *         states in the same node so they give the value defined
   *         above (i.e. l_pp[v] == l_pp[y] = the PP that either v or
   *         y emitted residue i) instead of l_pp[v] = PP that v
   *         emitted i, and l_pp[y] = PP that y emitted i.  And
   *         combine r_pp values for MATP_MP (v) and MATP_MR (y=v+2)
   *         states in an analogous way.
   */
  for(v = 0; v <= cm->M; v++) { 
    if(cm->sttype[v] == MP_st) { 
      for(i = 1; i <= L; i++) { 
	emit_mx->l_pp[v][i]   = FLogsum(emit_mx->l_pp[v][i], emit_mx->l_pp[v+1][i]); 
	emit_mx->l_pp[v+1][i] = emit_mx->l_pp[v][i];
      }
      for(j = 1; j <= L; j++) { 
	emit_mx->r_pp[v][j]   = FLogsum(emit_mx->r_pp[v][j], emit_mx->r_pp[v+2][j]); 
	emit_mx->r_pp[v+2][j] = emit_mx->r_pp[v][j];
      }
    }
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp2; fp2 = fopen("tmp.std_emitmx",  "w"); cm_emit_mx_Dump(fp2, cm, emit_mx); fclose(fp2); */
#endif

  return eslOK;
}


/* Function: cm_EmitterPosteriorHB()
 * Date:     EPN, Thu Oct  6 06:59:53 2011
 *
 * Purpose: Same as cm_EmitterPosterior() except HMM banded matrices
 *          are used. The main difference is that we have to be careful
 *          to stay within the bands because matrix cells outside 
 *          the bands do not exist (are not allocated). This requires
 *          keeping careful track of our offsets between the sequence
 *          position index and the corresponding indices in the matrix.
 * 
 * Args:     cm         - the model
 *           errbuf     - for error messages
 *           L          - length of the sequence
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           post       - pre-filled posterior cube
 *           emit_mx     - pre-allocated emit matrix, grown and filled-in here
 *           do_check   - if TRUE, return eslEFAIL if summed prob of any residue 
 *                        (before normalization) is < 0.98 or > 1.02.
 * 
 * Returns:  <eslOK>     on success.
 * Throws:   <eslERANGE> if required DP matrix size exceeds <size_limit>
 *           <eslFAIL>   if (do_check) and any residue check fails
 *           <eslEMEM>   if we run out of memory. 
 *           If !eslOK the l_pp and r_pp values are invalid.
 */
int 
cm_EmitterPosteriorHB(CM_t *cm, char *errbuf, int L, float size_limit, CM_HB_MX *post, CM_HB_EMIT_MX *emit_mx, int do_check)
{
  int    status;
  int    v, j, d; /* state, position, subseq length */
  int    i;       /* sequence position */
  int    ip_v;    /* offset i in banded matrix */
  int    ip_v2;   /* another offset i in banded matrix */
  int    jp_v;    /* offset j in banded matrix */
  int    jp_v2;   /* another offset j in banded matrix */
  int    dp_v;    /* offset d in banded matrix */
  int    in, ix;  /* temp min/max i */
  int    jn, jx;  /* temp min/max j */

  /* ptrs to band info, for convenience */
  int     *imin  = cm->cp9b->imin;  
  int     *imax  = cm->cp9b->imax;
  int     *jmin  = cm->cp9b->jmin;  
  int     *jmax  = cm->cp9b->jmax;
  CP9Bands_t *cp9b = cm->cp9b;  /* brief 26_0430-157: needed by hd_min()/hd_max() */
  int    **hdmin = cm->cp9b->hdmin;
  int    **hdmax = cm->cp9b->hdmax;
  
  /* grow the emit matrices based on the current sequence */
  if((status = cm_hb_emit_mx_GrowTo(cm, emit_mx, errbuf, cm->cp9b, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the emit matrices to IMPOSSIBLE */
  esl_vec_FSet(emit_mx->l_pp_mem, emit_mx->l_ncells_valid, IMPOSSIBLE);
  esl_vec_FSet(emit_mx->r_pp_mem, emit_mx->r_ncells_valid, IMPOSSIBLE);

  /* Step 1. Fill l_pp[v][i] and r_pp[v][i] with the posterior
   *         probability that state v emitted residue i either
   *         leftwise (l_pp) or rightwise (r_pp).
   */
  for(v = 0; v < cm->M; v++) { 
    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
      for(j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j - jmin[v];
	for(d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) { 
	  dp_v = d-hd_min(cp9b, v, jp_v);
	  i    = j-d+1;
	  assert(i >= imin[v] && i <= imax[v]);
	  ip_v = i - imin[v];
	  emit_mx->l_pp[v][ip_v] = FLogsum(emit_mx->l_pp[v][ip_v], post->dp[v][jp_v][dp_v]);
	}
      }
    }
    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
      for(j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	/* Peel first d: assign directly (avoids FLogsum(IMPOSSIBLE, x) == x) */
	emit_mx->r_pp[v][jp_v] = post->dp[v][jp_v][0];
	for(d = hd_min(cp9b, v, jp_v)+1; d <= hd_max(cp9b, v, jp_v); d++) {
	  dp_v = d-hd_min(cp9b, v, jp_v);
	  emit_mx->r_pp[v][jp_v] = FLogsum(emit_mx->r_pp[v][jp_v], post->dp[v][jp_v][dp_v]);
	}
      }
    }
  }
  /* factor in contribution of local ends, the EL state may have emitted this residue.
   * Note, the M deck is non-banded
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (j = 1; j <= L; j++) {
      i = j;
      for (d = 1; d <= j; d++, i--) { /* note: d >= 1, b/c EL emits 1 residue */
	emit_mx->l_pp[cm->M][i] = FLogsum(emit_mx->l_pp[cm->M][i], post->dp[cm->M][j][d]);
      }
    }
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp1; fp1 = fopen("tmp.std_unnorm_hbemitmx",  "w"); cm_hb_emit_mx_Dump(fp1, cm, emit_mx); fclose(fp1); */
#endif

  /* Step 2. Normalize l_pp and r_pp so that probability that
   *         each residue was emitted by any state is exactly
   *         1.0.
   */
  esl_vec_FSet(emit_mx->sum, (L+1), IMPOSSIBLE);
  for(v = 0; v < cm->M; v++) { /* we'll handle EL special */
    if(emit_mx->l_pp[v] != NULL) {
      for(i = imin[v]; i <= imax[v]; i++) { 
	ip_v = i - imin[v];
	emit_mx->sum[i] = FLogsum(emit_mx->sum[i], emit_mx->l_pp[v][ip_v]);
      }
    }
    if(emit_mx->r_pp[v] != NULL) {
      for(j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j - jmin[v];
	emit_mx->sum[j] = FLogsum(emit_mx->sum[j], emit_mx->r_pp[v][jp_v]);
      }
    }
  }
  /* Handle EL deck, remember it is non-banded */
  if(emit_mx->l_pp[cm->M] != NULL) { 
    for(i = 1; i <= L; i++) { 
      emit_mx->sum[i] = FLogsum(emit_mx->sum[i], emit_mx->l_pp[v][i]);
    }
  }
  
  /* perform the check, if nec */
  if(do_check) { 
    for(i = 1; i <= L; i++) { 
      if((sreEXP2(emit_mx->sum[i]) < 0.98) || (sreEXP2(emit_mx->sum[i]) > 1.02)) { 
	ESL_FAIL(eslFAIL, errbuf, "residue %d has summed prob of %5.4f (2^%5.4f).\nMay not be a DP coding bug, see 'Note:' on precision in cm_EmitterPosterior().\n", i, (sreEXP2(emit_mx->sum[i])), emit_mx->sum[i]);
      }
      printf("HB i: %d | total: %10.4f\n", i, (sreEXP2(emit_mx->sum[i])));
    }
    ESL_DPRINTF1(("#DEBUG: cm_EmitterPosteriorHB() check passed, all residues have summed probability of emission of between 0.98 and 1.02.\n"));
  }  

  /* normalize, using the sum vector */
  for(v = 0; v < cm->M; v++) { 
    if(emit_mx->l_pp[v] != NULL) {
      for(i = imin[v]; i <= imax[v]; i++) { 
	ip_v = i - imin[v];
	emit_mx->l_pp[v][ip_v] -= emit_mx->sum[i];
      }
    }
    if(emit_mx->r_pp[v] != NULL) {
      for(j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j - jmin[v];
	emit_mx->r_pp[v][jp_v] -= emit_mx->sum[j];
      }
    }
  }
  /* Handle EL deck, remember it is non-banded */
  if(emit_mx->l_pp[cm->M] != NULL) { 
    for(i = 1; i <= L; i++) { 
      emit_mx->l_pp[cm->M][i] -= emit_mx->sum[i];
    }
  }

  /* Step 3. Combine l_pp values for MATP_MP (v) and MATP_ML (y=v+1)
   *         states in the same node so they give the value defined
   *         above (i.e. l_pp[v] == l_pp[y] = the PP that either v or
   *         y emitted residue i) instead of l_pp[v] = PP that v
   *         emitted i, and l_pp[y] = PP that y emitted i.  And
   *         combine r_pp values for MATP_MP (v) and MATP_MR (y=v+2)
   *         states in an analogous way.
   */
  for(v = 0; v <= cm->M; v++) { 
    if(cm->sttype[v] == MP_st) { 
      /* we only change l_pp[v][i] and l_pp[v+1][i] if i is within both
       * state v and v+1's i band.
       */
      if(imax[v] >= 1 && imax[v+1] >= 1) { 
	in = ESL_MAX(imin[v], imin[v+1]); 
	ix = ESL_MIN(imax[v], imax[v+1]);
	for(i = in; i <= ix; i++) { 
	  ip_v  = i - imin[v];
	  ip_v2 = i - imin[v+1];
	  emit_mx->l_pp[v][ip_v]    = FLogsum(emit_mx->l_pp[v][ip_v], emit_mx->l_pp[v+1][ip_v2]); 
	  emit_mx->l_pp[v+1][ip_v2] = emit_mx->l_pp[v][ip_v];
	}
      }
      /* we only change r_pp[v][j] and r_pp[v+2][j] if j is within both
       * state v and v+2's j band.
       */
      if(jmax[v] >= 1 && jmax[v+2] >= 1) { 
	jn = ESL_MAX(jmin[v], jmin[v+2]); 
	jx = ESL_MIN(jmax[v], jmax[v+2]);
	for(j = jn; j <= jx; j++) { 
	  jp_v  = j - jmin[v];
	  jp_v2 = j - jmin[v+2];
	  emit_mx->r_pp[v][jp_v]    = FLogsum(emit_mx->r_pp[v][jp_v], emit_mx->r_pp[v+2][jp_v2]); 
	  emit_mx->r_pp[v+2][jp_v2] = emit_mx->r_pp[v][jp_v];
	}
      }
    }
  }

#if eslDEBUGLEVEL >= 3
  /* Uncomment to dump matrix to file. This could be very large, so be careful. */
  /* FILE *fp2; fp2 = fopen("tmp.std_hbemitmx",  "w"); cm_hb_emit_mx_Dump(fp2, cm, emit_mx); fclose(fp2); */
#endif

  return eslOK;
}

/* Function: cm_PostCode()
 * Date:     EPN 05.25.06 based on SRE's Postcode() 
 *           from HMMER's postprob.c
 *
 * Purpose: Given a parse tree and a filled emit matrix calculate two
 *           strings that represents the confidence values on each
 *           aligned residue in the sequence.
 *           
 *           The emit_mx values are:
 *           l_pp[v][i]: log of the posterior probability that state v emitted
 *                       residue i leftwise either at (if a match state) or
 *                       *after* (if an insert state) the left consensus
 *                       position modeled by state v's node.
 *
 *           r_pp[v][i]: log of the posterior probability that state v emitted
 *                       residue i rightwise either at (if a match state) or
 *                       *before* (if an insert state) the right consensus
 *                       position modeled by state v's node.
 *
 *           l_pp[v] is NULL for states that do not emit leftwise  (B,S,D,E,IR,MR)
 *           r_pp[v] is NULL for states that do not emit rightwise (B,S,D,E,IL,ML)
 *
 *           The PP string is 0..L-1  (L = len of target seq),
 *           so its in the coordinate system of the sequence string;
 *           off by one from dsq.
 *           
 *           Values are 0,1,2,3,4,5,6,7,8,9,*:
 *           '0' = [0.00-0.05)
 *           '1' = [0.05-0.15)
 *           '2' = [0.15-0.25)
 *           '3' = [0.25-0.35)
 *           '4' = [0.35-0.45)
 *           '5' = [0.45-0.55)
 *           '6' = [0.55-0.65)
 *           '7' = [0.65-0.75)
 *           '8' = [0.75-0.85)
 *           '9' = [0.85-0.95)
 *           '*' = [0.95-1.00)
 *
 *           cm_PostCodeHB() is nearly the same function with the
 *           difference that HMM bands were used for the alignment,
 *           so we have to deal with offset issues.
 *
 *           Renamed from CMPostCode() [EPN, Wed Sep 14 06:20:35 2011].
 *
 * Args:     cm         - the model 
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitized sequence [1..L]   
 *           L          - length of the dsq to align
 *           emit_mx    - the pre-filled emit matrix, must be non-NULL if do_optacc
 *           tr         - the parstree with the emissions we're setting PPs for
 *           ret_ppstr  - RETURN: a string of the PP code values (0..L-1) 
 *           ret_avgp   - RETURN: the average PP of all aligned residues
 *
 * Returns:  <eslOK>     on success.
 * Throws:   <eslEINVAL> if a posterior probability is > 1.01 or less than -0.01. 
 */
char
Fscore2postcode(float sc)
{
  float p = FScore2Prob(sc, 1.);
  return (p + 0.05 >= 1.0) ? '*' :  (char) ((p + 0.05) * 10.0) + '0';
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

int
cm_PostCode(CM_t *cm, char *errbuf, int L, CM_EMIT_MX *emit_mx, Parsetree_t *tr, char **ret_ppstr, float *ret_avgp)
{
  int   status;
  int   x, v, i, j, r; /* counters */
  char *ppstr;       /* the PP string, created here */
  float p;           /* a probability */
  float sum_logp;    /* log of summed probability of all residues emitted thus far */

  ESL_ALLOC(ppstr, (L+1) * sizeof(char)); 
  sum_logp = IMPOSSIBLE;

  /* go through each node of the parsetree and determine post code for emissions */
  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];

      /* Only P, L, R, and EL states have emissions. */
      if(cm->sttype[v] == EL_st) { /* EL state, we have to handle this guy special */
	for(r = i; r <= j; r++) { /* we have to annotate from residues i..j */
	  ppstr[r-1] = Fscore2postcode(emit_mx->l_pp[v][r]);
	  sum_logp   = FLogsum(sum_logp, emit_mx->l_pp[v][r]);
	  /* make sure we've got a valid probability */
	  p = FScore2Prob(emit_mx->l_pp[v][r], 1.);
	  if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for EL state v: %d residue r: %d > 1.00 (%.2f)", v, r, p);
	  if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for EL state v: %d residue r: %d < 0.00 (%.2f)", v, r, p);
	}
      }
      if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) { 
	ppstr[i-1] = Fscore2postcode(emit_mx->l_pp[v][i]);
	sum_logp   = FLogsum(sum_logp, emit_mx->l_pp[v][i]);
	/* make sure we've got a valid probability */
	p = FScore2Prob(emit_mx->l_pp[v][i], 1.);
	if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for left state v: %d residue i: %d > 1.00 (%.2f)", v, i, p);
	if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for left state v: %d residue i: %d < 0.00 (%.2f)", v, i, p);
      }
      if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
	ppstr[j-1] = Fscore2postcode(emit_mx->r_pp[v][j]);
	sum_logp   = FLogsum(sum_logp, emit_mx->r_pp[v][j]);
	/* make sure we've got a valid probability */
	p = FScore2Prob(emit_mx->r_pp[v][j], 1.);
	if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for right state v: %d residue i: %d > 1.00 (%.2f)", v, j, p);
	if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for right state v: %d residue i: %d < 0.00 (%.2f)", v, j, p);
      }
    }
  ppstr[L] = '\0';

  if(ret_ppstr != NULL) *ret_ppstr = ppstr; else free(ppstr);
  if(ret_avgp  != NULL) *ret_avgp  = sreEXP2(sum_logp) / (float) L;
  return eslOK;
  
 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "cm_Postcode(): Memory allocation error.");
  return status; /* never reached */
}


int
cm_PostCodeHB(CM_t *cm, char *errbuf, int L, CM_HB_EMIT_MX *emit_mx, Parsetree_t *tr, char **ret_ppstr, float *ret_avgp)
{
  int   status;
  int   x, v, i, j, r; /* counters */
  char *ppstr;       /* the PP string, created here */
  float p;           /* a probability */
  float sum_logp;    /* log of summed probability of all residues emitted thus far */

  /* variables used for HMM bands */
  int ip_v, jp_v; /* i, j offset within bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b = cm->cp9b;
  int     *imin  = cp9b->imin;  
  int     *imax  = cp9b->imax;  
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;  

  ESL_ALLOC(ppstr, (L+1) * sizeof(char)); 
  sum_logp = IMPOSSIBLE;

  /* go through each node of the parsetree and determine post code for emissions */
  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];

      /* Only P, L, R, and EL states have emissions. */
      if(cm->sttype[v] == EL_st) { /* EL state, we have to handle this guy special */
	for(r = i; r <= j; r++) { /* we have to annotate from residues i..j */
	  /* remember the EL deck is non-banded */
	  ppstr[r-1] = Fscore2postcode(emit_mx->l_pp[v][r]);
	  sum_logp   = FLogsum(sum_logp, emit_mx->l_pp[v][r]);
	  /* make sure we've got a valid probability */
	  p = FScore2Prob(emit_mx->l_pp[v][r], 1.);
	  if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for EL state v: %d residue r: %d > 1.00 (%.2f)", v, r, p);
	  if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for EL state v: %d residue r: %d < 0.00 (%.2f)", v, r, p);
	}
      }
      if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) { 
	ip_v = i - imin[v];
	assert(i >= imin[v] && i <= imax[v]);
	ESL_DASSERT1((i >= imin[v] && i <= imax[v]));
	ppstr[i-1] = Fscore2postcode(emit_mx->l_pp[v][ip_v]);
	sum_logp   = FLogsum(sum_logp, emit_mx->l_pp[v][ip_v]);
	/* make sure we've got a valid probability */
	p = FScore2Prob(emit_mx->l_pp[v][ip_v], 1.);
	if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for left state v: %d residue i: %d > 1.00 (%.2f)", v, i, p);
	if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for left state v: %d residue i: %d < 0.00 (%.2f)", v, i, p);
      }
      if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
	jp_v = j - jmin[v];
	assert(j >= jmin[v] && j <= jmax[v]);
	ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	ppstr[j-1] = Fscore2postcode(emit_mx->r_pp[v][jp_v]);
	sum_logp   = FLogsum(sum_logp, emit_mx->r_pp[v][jp_v]);
	/* make sure we've got a valid probability */
	p = FScore2Prob(emit_mx->r_pp[v][jp_v], 1.);
	if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for right state v: %d residue i: %d > 1.00 (%.2f)", v, j, p);
	if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_PostCode(): probability for right state v: %d residue i: %d < 0.00 (%.2f)", v, j, p);
      }
    }
  ppstr[L] = '\0';

  if(ret_ppstr != NULL) *ret_ppstr = ppstr; else free(ppstr);
  if(ret_avgp  != NULL) *ret_avgp  = sreEXP2(sum_logp) / (float) L;
  ESL_DPRINTF1(("#DEBUG: cm_PostcodeHB(): average pp %.4f\n", sreEXP2(sum_logp) / (float) L));
  /*printf("cm_PostcodeHB(): average pp %.4f\n", sreEXP2(sum_logp) / (float) L);*/

  return eslOK;
  
 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "cm_PostcodeHB(): Memory allocation error.");
  return status; /* never reached */
}

/* Function: cm_InitializeOptAccShadowDZero()
 * Date:     EPN, Fri Nov 11 13:09:14 2011
 *
 * Purpose:  Initialize a optimal accuracy shadow (traceback) matrix
 *           for d == 0, based on transition scores. Optimal accuracy
 *           traceback matrices are special when d==0 because only
 *           emissions contribute to score so the value when d==0 is
 *           always IMPOSSIBLE. So d==0 cells are never modified
 *           during the OA DP recursion.
 *
 *           In this function we determine the appropriate state to
 *           traceback to for d==0 for all states v and endpoints
 *           j. If local ends are off, this is trivial; it is simply
 *           the child state y for which StateDelta(y) == 0 (there is
 *           always exactly 1 such child state for each v). If any
 *           such state is entered for d == 0 in the optimally
 *           accurate parsetree, the parse will continue along
 *           delete->delete transitions (all with d==0) until an E_st
 *           (or E_st's if we go through a B_st) is reached.
 *      
 *           If local ends are on, it is more complex because we
 *           could do a local end instead of a string of deletes until
 *           an end is reached. We determine the score of the
 *           transitions from the current state v through y to the
 *           nearest E_st(s) and if it is less than the score for
 *           entering a EL state we set the shadow matrix to y, else
 *           we set it to USED_EL.
 *           
 *           In some cases, we initialize to USED_EL for states v for
 *           which ELs are illegal (not a MATP_MP, MATL_ML, MATR_MR,
 *           BEGL_S or BEGR_nd). This means that the eventual optimal
 *           accuracy parsetree may contain an illegal EL, but I think
 *           this is unavoidable. 
 *           
 *           Upon entrance, yshadow should be initialized to USED_EL
 *           for all values.
 *
 *           Note that if we didn't call this function, the optimally
 *           accurate parsetree would not be affected, nor its score.
 *           This function is only useful because it affects the
 *           output of the parsetree's alignment by only using a zero
 *           length EL transitions only when it is less expensive than
 *           a string of deletes.
 *        
 *           If called by a truncated optimal accuracy function
 *           (cm_TrOptAccAlign()), yshadow is really a <Jyshadow>
 *           matrix from a CM_TR_SHADOW_MX object. Otherwise it is a
 *           <yshadow> matrix from a CM_SHADOW_MX object.
 *
 * Args:     cm         - the model, used only for its alphabet and null model
 *           errbuf     - for reporting errors 
 *           yshadow    - the shadow matrix to updated, only values for which
 *                        d==0 will be modified. 
 *           L          - length of the sequence we're aligning
 * 
 *    
 * Returns:  eslOK on success
 *
 * Throws:   eslEMEM on memory error.
 */
int
cm_InitializeOptAccShadowDZero(CM_t *cm, char *errbuf, char ***yshadow, int L)
{
  int   status;
  float *esc;  /* [0..v..M-1] summed transition score for getting from v to nearest E_st(s) 
		* through only delete states */
  float endsc; /* score for transitioning to an EL state */
  int have_el; /* are local ends on? */
  int v;       /* state counter */
  int j;       /* sequence position */
  int y, z;    /* BEGL_S and BEGR_S states */
  int sd;      /* StateDelta(v) */
  int yoffset; /* child state index */

  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el) {
    ESL_ALLOC(esc, sizeof(float) * cm->M);
    esl_vec_FSet(esc, cm->M, 0.);
    /* determine score for transitioning to an EL (same for all legal states) */
    v = 0; while(! NOT_IMPOSSIBLE(cm->endsc[v])) v++;
    endsc = cm->endsc[v];
    /*printf("endsc: %.4f end %.4f\n", endsc, cm->end[v]);*/
  }
  else { 
    esc = NULL;
    endsc = IMPOSSIBLE;
  }

  for(v = cm->M-1; v >= 0; v--) { 
    sd = StateDelta(cm->sttype[v]);
    if(cm->sttype[v] == E_st) { 
      if(esc != NULL) esc[v] = 0.;
    }
    else if(cm->sttype[v] == B_st) { 
      if(esc != NULL) { 
	y = cm->cfirst[v]; /* left  subtree */
	z = cm->cnum[v];   /* right subtree */
	esc[v] = esc[y] + esc[z];
      }
    }
    else { 
      /* determine the one and only child state y for which StateDelta(y) == 0 */
      y = cm->cfirst[v];
      while(StateDelta(cm->sttype[y]) != 0) y++;
      yoffset = y-cm->cfirst[v];
      assert(cm->ndidx[v] == (cm->ndidx[y]-1));
      if(esc != NULL) { 
	esc[v] = esc[y] + cm->tsc[v][yoffset];
	if(endsc > esc[v]) yoffset = USED_EL;
	/* else yoffset is not changed */

	/*printf("EL: %10.4f  d->d->e %10.4f  ", endsc, esc[v]);
	  if(yoffset != USED_EL) printf("  path for v: %4d %4s %2s is through deletes!\n", v, Nodetype(cm->ndtype[cm->ndidx[v]]), Statetype(cm->sttype[v]));
	  else printf("\n");
	*/
      }
      for(j = sd; j <= L; j++) yshadow[v][j][sd] = yoffset;
    }
  }
  
  if(esc != NULL) free(esc);
  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Out of memory");
}


/* Function: cm_InitializeOptAccShadowDZeroHB()
 * Date:     EPN, Fri Nov 11 14:00:55 2011
 *
 * Purpose:  Same as cm_InitializeOptAccShadowDZero() but for HMM
 *           banded matrices, see that function for more information.
 *
 * Args:     cm         - the model, used only for its alphabet and null model
 *           cp9b       - CP9 Bands for current sequence
 *           errbuf     - for reporting errors 
 *           yshadow    - the shadow matrix to updated, only values for which
 *                        d==0 will be modified. 
 *           L          - length of the sequence we're aligning
 * 
 *    
 * Returns:  eslOK on success
 *
 * Throws:   eslEMEM on memory error.
 */
int
cm_InitializeOptAccShadowDZeroHB(CM_t *cm, CP9Bands_t *cp9b, char *errbuf, char ***yshadow, int L)
{
  int   status;
  float *esc;  /* [0..v..M-1] summed transition score for getting from v to nearest E_st(s) 
		* through only delete states */
  float endsc; /* score for transitioning to an EL state */
  int have_el; /* are local ends on? */
  int v;       /* state counter */
  int j;       /* sequence position */
  int y, z;    /* BEGL_S and BEGR_S states */
  int sd;      /* StateDelta(v) */
  int yoffset; /* child state index */

  /* variables needed because we've got HMM bands */
  int sdr;     /* StateRightDelta(v) */
  int jp_v;    /* j offset for state v given HMM bands */
  int jp_y;    /* j offset for state y given HMM bands */
  int dp_v;    /* d offset for state v given HMM bands */

  /* pointers to cp9b data for convenience */
  int         *jmin = cp9b->jmin;
  int         *jmax = cp9b->jmax;
  int       **hdmin = cp9b->hdmin;
  int       **hdmax = cp9b->hdmax;

  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el) {
    ESL_ALLOC(esc, sizeof(float) * cm->M);
    esl_vec_FSet(esc, cm->M, IMPOSSIBLE);
    /* determine score for transitioning to an EL (same for all legal states) */
    v = 0; while(! NOT_IMPOSSIBLE(cm->endsc[v])) v++;
    endsc = cm->endsc[v];
    /*printf("endsc: %.4f end %.4f\n", endsc, cm->end[v]);*/
  }
  else { 
    esc = NULL;
    endsc = IMPOSSIBLE;
  }

  for(v = cm->M-1; v >= 0; v--) { 
    if(cm->cp9b->Jvalid[v]) { /* only valid v values will have non-impossible esc[v] values */
      sd  = StateDelta(cm->sttype[v]);
      sdr = StateRightDelta(cm->sttype[v]);
      if(cm->sttype[v] == E_st) { 
	if(esc != NULL) esc[v] = 0.;
      }
      else if(cm->sttype[v] == B_st) { 
	if(esc != NULL) {
	  y = cm->cfirst[v]; /* left  subtree */
	  z = cm->cnum[v];   /* right subtree */
	  esc[v] = esc[y] + esc[z];
	}
      }
      else { 
	/* determine the one and only child state y for which StateDelta(y) == 0 */
	y = cm->cfirst[v];
	while(StateDelta(cm->sttype[y]) != 0) y++;
	yoffset = y-cm->cfirst[v];
	assert(cm->ndidx[v] == (cm->ndidx[y]-1));
	if(esc != NULL) { 
	esc[v] = esc[y] + cm->tsc[v][yoffset];
	if(endsc > esc[v]) yoffset = USED_EL;
	/* else yoffset is not changed */

#if 0 
	printf("EL: %10.4f  d->d->e %10.4f  ", endsc, esc[v]);
	if(yoffset != USED_EL) printf("  path for v %4d %4s %2s is through deletes!\n", v, Nodetype(cm->ndtype[cm->ndidx[v]]), Statetype(cm->sttype[v]));
	else printf("\n");
#endif

	}
	for(j = ESL_MAX(sd, jmin[v]); j <= jmax[v]; j++) { 
	  jp_v = j-jmin[v];
	  if(hd_min(cp9b, v, jp_v) <= hd_max(cp9b, v, jp_v)) { /* at least one valid d exists for this v and j */
	    if((j-sdr) >= jmin[y] && (j-sdr) <= jmax[y]) { /* j-sdr is valid for state y */
	      jp_y = j - sdr - jmin[y]; 
	      if(sd >= hd_min(cp9b, v, jp_v) && sd <= hd_max(cp9b, v, jp_v) && /* d==sd is valid for state v and end posn j */
		 0  >= hd_min(cp9b, y, jp_y) &&  0 <= hd_max(cp9b, y, jp_y)) { /* d==0  is valid for state y and end posn j-sdr */
		dp_v = sd - hd_min(cp9b, v, jp_v);
		yshadow[v][jp_v][dp_v] = yoffset;
	      }
	    }
	  }
	}
      }
    }
  }
  
  if(esc != NULL) free(esc);
  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Out of memory");
}


/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_ALIGN_BENCHMARK
/* Next line is not optimized (debugging on) on MacBook Pro:
 * gcc   -o benchmark-align -std=gnu99 -g -Wall -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DIMPL_ALIGN_BENCHMARK cm_dpalign.c -linfernal -lhmmer -leasel -lm
 * Next line is optimized (debugging not on) on wyvern:
 * gcc   -o benchmark-align -std=gnu99 -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -pthread -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DIMPL_ALIGN_BENCHMARK cm_dpalign.c -linfernal -lhmmer -leasel -lm 
 * ./benchmark-align <cmfile>
 */

#include <esl_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_sqio.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "hmmer.h"

#include "infernal.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "configure CM/HMM for local alignment", 0 },
  { "--cykout",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "run CYKOutside, to make sure it agrees with CYK (Inside)", 0 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use posterior sums during HMM band calculation (widens bands)", 0 },
  { "--dlev",    eslARG_INT,    "0",   NULL, "0<=n<=3",NULL,NULL,NULL, "set verbosity of debugging print statements to <n>", 0 },
  { "--hmmcheck",eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "check that HMM posteriors are correctly calc'ed", 0 },
  { "--cmcheck", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "check that CM posteriors are correctly calc'ed", 0 },
  { "--optacc",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute optimal accuracy HMM banded alignment alg", 0 },
  { "--tau",     eslARG_REAL,   "5e-6",NULL, "0<x<1",NULL, NULL, NULL, "set tail loss prob for HMM bands to <x>", 0 },
  { "--post",   eslARG_NONE,    FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute fast float HMM banded Inside/Outside alignment algs", 0 },
  { "--mxsize",  eslARG_REAL, "256.0", NULL, "x>0.",NULL,  NULL, NULL, "set maximum allowable DP matrix size to <x> (Mb)", 0 },
  { "--nonbanded",eslARG_NONE,  FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute non-banded alignment algorithms", 0 },
  { "--tr",       eslARG_NONE,  FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump parsetrees to stdout", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "benchmark driver for fast HMM banded CYK alignment and scanning algorithm";

int 
main(int argc, char **argv)
{
  int status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  CM_t           *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_ALPHABET   *abc     = NULL;
  int             i;
  float           sc;
  float           pp;
  char           *cmfile  = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  CM_FILE        *cmfp  = NULL;  /* open input CM file stream */
  ESL_SQFILE     *sqfp  = NULL;  /* open sequence input file stream */
  ESL_SQ         *sq    = NULL;  /* a sequence */
  int             L;             /* length of sequence */
  char            errbuf[eslERRBUFSIZE];
  float           size_limit = esl_opt_GetReal(go, "--mxsize");
  int             do_check   = esl_opt_GetBoolean(go, "--cmcheck");
  float           parsetree_sc, parsetree_struct_sc;
  Parsetree_t    *tr    = NULL;

  /* open CM file */
  if ((status = cm_file_Open(cmfile, NULL, FALSE, &(cmfp), errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = cm_file_Read(cmfp, TRUE, &abc, &cm))                != eslOK) cm_Fail(cmfp->errbuf);
  cm_file_Close(cmfp);

  /* open the sequence file */
  status = esl_sqfile_OpenDigital(cm->abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp);
  if (status == eslENOTFOUND)    esl_fatal("File %s doesn't exist or is not readable\n", seqfile);
  else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of sequence file %s\n", seqfile);
  else if (status == eslEINVAL)  esl_fatal("Can't autodetect stdin or .gz."); 
  else if (status != eslOK)      esl_fatal("Sequence file open failed with error %d.\n", status);
  
  /* configure CM */
  cm->align_opts  |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--sums")) cm->align_opts |= CM_ALIGN_SUMS;
  if(esl_opt_GetBoolean(go, "-l")) { 
    cm->config_opts  |= CM_CONFIG_LOCAL;
    cm->config_opts  |= CM_CONFIG_HMMLOCAL;
    cm->config_opts  |= CM_CONFIG_HMMEL;
  }
  if(esl_opt_GetBoolean(go, "--hmmcheck")) cm->align_opts |= CM_ALIGN_CHECKFB;
  if(esl_opt_GetBoolean(go, "--cmcheck"))  cm->align_opts |= CM_ALIGN_CHECKINOUT;
  cm->tau = esl_opt_GetReal(go, "--tau");

  if((status = cm_Configure(cm, errbuf, -1)) != eslOK) cm_Fail(errbuf);

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  i = 0;
  sq = esl_sq_CreateDigital(cm->abc);
  while((status = esl_sqio_Read(sqfp, sq)) == eslOK) { 
    i++;
    L = sq->n;

    esl_stopwatch_Start(w);
    if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK) cm_Fail(errbuf);
    esl_stopwatch_Stop(w);
    printf("%4d %-30s %17s", i, "Exptl Band calc:", "");
    esl_stopwatch_Display(stdout, w, "CPU time: ");
      
    esl_stopwatch_Start(w);
    if((status = cm_AlignHB(cm, errbuf, sq->dsq, L, size_limit, FALSE, FALSE, cm->hb_mx, cm->hb_shmx, NULL, NULL, NULL, NULL, &tr, &pp, &sc)) != eslOK) cm_Fail(errbuf);
    printf("%4d %-30s %10.4f bits ", (i), "cm_AlignHB() CYK:", sc);
    esl_stopwatch_Stop(w);
    esl_stopwatch_Display(stdout, w, " CPU time: ");

    if(esl_opt_GetBoolean(go, "--tr")) ParsetreeDump(stdout, tr, cm, sq->dsq);
    ParsetreeScore(cm, NULL, NULL, tr, sq->dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
    FreeParsetree(tr);
    printf("Parsetree score      : %.4f           (FULL LENGTH CYK)\n", parsetree_sc);

    if(esl_opt_GetBoolean(go, "--cykout")) { 
      esl_stopwatch_Start(w);
      if((status = cm_CYKOutsideAlignHB(cm, errbuf, sq->dsq, L, size_limit, TRUE, cm->hb_omx, cm->hb_mx, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i), "cm_Align() CYK:", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
    }

    if(esl_opt_GetBoolean(go, "--nonbanded")) {
      esl_stopwatch_Start(w);
      if((status = cm_Align(cm, errbuf, sq->dsq, L, size_limit, FALSE, FALSE, cm->nb_mx, cm->nb_shmx, NULL, cm->nb_emx, NULL, NULL, &tr, &pp, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i), "cm_Align() CYK:", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      if(esl_opt_GetBoolean(go, "--tr")) ParsetreeDump(stdout, tr, cm, sq->dsq);
      ParsetreeScore(cm, NULL, NULL, tr, sq->dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
      FreeParsetree(tr);
      printf("Parsetree score      : %.4f           (FULL LENGTH CYK)\n", parsetree_sc);

      if(esl_opt_GetBoolean(go, "--cykout")) { 
	esl_stopwatch_Start(w);
	if((status = cm_CYKOutsideAlign(cm, errbuf, sq->dsq, L, size_limit, TRUE, cm->nb_omx, cm->nb_mx, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i), "cm_Align() CYK:", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }
    }
    printf("\n");

    if(esl_opt_GetBoolean(go, "--post")) {
      esl_stopwatch_Start(w);
      /* need alpha matrix from Inside to do Outside */
      if((status = cm_InsideAlignHB(cm, errbuf, sq->dsq, L, size_limit, cm->hb_mx, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i), "cm_InsideAlignHB():", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      esl_stopwatch_Start(w);
      /* need alpha matrix from Inside to do Outside */
      if((status = cm_OutsideAlignHB(cm, errbuf, sq->dsq, L, size_limit, do_check, cm->hb_omx, cm->hb_mx, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i), "cm_OutsideAlignHB():", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      if(esl_opt_GetBoolean(go, "--nonbanded")) { 
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	if((status = cm_InsideAlign(cm, errbuf, sq->dsq, L, size_limit, cm->nb_mx, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i), "cm_InsideAlign():", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	  
	esl_stopwatch_Start(w);
	/* need alpha matrix from Inside to do Outside */
	if((status = cm_OutsideAlign(cm, errbuf, sq->dsq, L, size_limit, do_check, cm->nb_omx, cm->nb_mx, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i), "cm_OutsideAlign():", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }
    }

    if(esl_opt_GetBoolean(go, "--optacc")) {
      esl_stopwatch_Start(w);
      if((status = cm_AlignHB(cm, errbuf, sq->dsq, L, size_limit, TRUE, FALSE, cm->hb_mx, cm->hb_shmx, cm->hb_omx, cm->hb_emx, NULL, NULL, &tr, &pp, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f avgpp ", (i), "cm_AlignHB() OA:", pp);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      if(esl_opt_GetBoolean(go, "--tr")) ParsetreeDump(stdout, tr, cm, sq->dsq);
      ParsetreeScore(cm, NULL, NULL, tr, sq->dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
      FreeParsetree(tr);
      printf("Parsetree score      : %.4f           (FULL LENGTH OPTACC)\n", parsetree_sc);

      if(esl_opt_GetBoolean(go, "--nonbanded")) { 
	esl_stopwatch_Start(w);
	if((status = cm_Align(cm, errbuf, sq->dsq, L, size_limit, TRUE, FALSE, cm->nb_mx, cm->nb_shmx, cm->nb_omx, cm->nb_emx, NULL, NULL, &tr, &pp, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f avgpp ", (i), "cm_Align() OA:", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	if(esl_opt_GetBoolean(go, "--tr")) ParsetreeDump(stdout, tr, cm, sq->dsq);
	ParsetreeScore(cm, NULL, NULL, tr, sq->dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
	FreeParsetree(tr);
	printf("Parsetree score      : %.4f           (FULL LENGTH OPTACC)\n", parsetree_sc);
      }
    }
    printf("\n");
    esl_sq_Reuse(sq);
  }
  if(status != eslEOF) cm_Fail("ERROR reading sequence file, sequence number %d\n", i);

  FreeCM(cm);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  esl_sqfile_Close(sqfp);

  return 0;
}
#endif /*IMPL_ALIGN_BENCHMARK*/

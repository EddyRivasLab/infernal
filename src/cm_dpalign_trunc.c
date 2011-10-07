/* cm_dpalign_trunc.c
 * 
 * DP functions for truncated HMM banded and non-banded, non-D&C CM
 * alignment of a full target sequence. 
 *
 * All functions use a DP matrix and or shadow matrix, either
 * non-banded (CM_TR_MX, CM_TR_SHADOW_MX) or HMM banded (CM_TR_HB_MX,
 * CM_TR_HB_SHADOW_MX).  The HMM banded matrices only have cells
 * within bands allocated, and further for each state v, J, L, R, or T
 * matrix cells are only available if cp9b->do_{J,L,R,T}[v]. The bands
 * are derived from a HMM Forward/Backward alignment of the target
 * sequence and are stored in a CP9Bands_t object, a pointer to which
 * must exist in the cm (CM_t object).
 * 
 * The non-banded, non-D&C alignment functions are mainly useful for
 * understanding and/or debugging the HMM banded versions.  These are
 * consistent (same logic/code organization) with their HMM banded
 * counterparts. They are memory intensive. For small memory
 * non-banded alignment functions see truncyk.c.
 *
 * Many of the functions here were based on analogous ones for
 * standard (non-truncated) CYK/Inside/Outside alignment in
 * cm_dpalign.c.
 *
 * List of functions: 
 * non-banded version          HMM banded version
 * -------------------------   -----------------------
 * cm_tr_alignT()              cm_tr_alignT_hb()
 * cm_TrAlign()                cm_TrAlignHB()
 * cm_TrCYKInsideAlign()       cm_TrCYKInsideAlignHB()
 * cm_TrInsideAlign()          cm_TrInsideAlignHB()
 * cm_TrOptAccAlign()          cm_TrOptAccAlignHB()
 * cm_TrCYKOutsideAlign()*     cm_TrCYKOutsideAlignHB()*
 * cm_TrOutsideAlign()         cm_TrOutsideAlignHB()
 * cm_TrPosterior()            cm_TrPosteriorHB()  
 * cm_SampleTrParsetree()      cm_SampleTrParsetreeHB()
 *
 * * cm_TrCYKOutsideAlign() and cm_TrCYKOutsideAlignHB() are for
 * reference and debugging only they're not called by any of the main
 * Infernal programs, only by test programs.
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

#define DEBUG1  0
#define DEBUG2  0
#define DOTRACE 0

static int cm_tr_alignT   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
			   CM_TR_EMIT_MX *emit_mx, char opt_mode, Parsetree_t **ret_tr, float *ret_sc);
static int cm_tr_alignT_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, 
			   CM_TR_HB_EMIT_MX *emit_mx, char opt_mode, Parsetree_t **ret_tr, float *ret_sc);

/* Function: cm_tr_alignT()
 * Date:     EPN, Sat Sep 10 11:25:37 2011
 *           EPN, Sun Nov 18 19:21:30 2007 [cm_alignT()]
 *
 * Note:     based on insideT() [SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]]
 *
 * Purpose: Call either cm_TrCYKInsideAlign() (if !<do_optacc>), or
 *           cm_TrOptAccAlign() (if <do_optacc>), get vjd shadow
 *           matrix; then trace back and append to an existing but
 *           empty parsetree tr.  The full sequence 1..L will be
 *           aligned.
 *        
 *           If <do_optacc>==TRUE then emit_mx must != NULL and
 *           <opt_mode> must not be TRMODE_UNKNOWN, it will have been
 *           determined by caller from a cm_TrInsideAlign() call and
 *           passed in. If <do_optacc> is FALSE, then we're doing CYK
 *           alignment and we don't know the truncation mode of the
 *           alignment mode yet, so <opt_mode> will be passed in as
 *           TRMODE_UNKNOWN.
 *
 * Args:     cm         - the model 
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitized sequence [1..L]   
 *           L          - length of the dsq to align
 *           size_limit - max size in Mb for DP matrix
 *           do_optacc  - TRUE to align with optimal accuracy, else use CYK
 *           mx         - the DP matrix to fill in
 *           shmx       - the shadow matrix to fill in
 *           emit_mx    - the pre-filled emit matrix, must be non-NULL if do_optacc
 *           opt_mode   - the optimal alignment mode, if unknown TRMODE_UNKNOWN is passed
 *           ret_tr     - RETURN: the optimal parsetree
 *           ret_sc     - RETURN: optimal score (CYK if !do_optacc, else avg PP of all 1..L residues) 
 * 
 * Returns:  <eslOK>     on success.
 * Throws:   <eslERANGE> if required DP matrix size exceeds <size_limit>, in 
 *                       this case, alignment has been aborted, ret_* variables are not valid
 *           <eslEINVAL> on traceback problem: bogus state
 */
int
cm_tr_alignT(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
	     CM_TR_EMIT_MX *emit_mx, char opt_mode, Parsetree_t **ret_tr, float *ret_sc)
{
  int       status;
  Parsetree_t *tr = NULL;       /* the parsetree */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda_i;             /* stack that tracks bifurc parent of a right start */
  ESL_STACK *pda_c;             /* stack that tracks mode of bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, seq positions */
  int       k;			/* right subtree len for bifurcs */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  /* variables specific to truncated version */
  char      mode;               /* current truncation mode: TRMODE_J | TRMODE_L | TRMODE_R | TRMODE_T */
  char      prvmode, nxtmode;   /* previous, next truncation mode */
  int       b, Jb, Lb, Rb, Tb;  /* entry states for best overall, J,L,R,T alignment */

  if(do_optacc) {
    if((status = cm_TrOptAccAlign(cm, errbuf, dsq, L, 
				  size_limit,   /* max size of DP matrix */
				  opt_mode,     /* marginal mode of optimal alignment */
				  mx,	        /* the DP matrix, to expand and fill-in */
				  shmx,	        /* the shadow matrix, to expand and fill-in */
				  emit_mx,      /* pre-calc'ed emit matrix */
				  &b,           /* the entry point for optimal alignment */
				  &sc))         /* avg post prob of all emissions in optimally accurate parsetree */
       != eslOK) return status;
    /* set root state of entire parse */
    v    = b;
    mode = opt_mode;
  }
  else { 
    if((status = cm_TrCYKInsideAlign(cm, errbuf, dsq, L,
				     size_limit,         /* max size of DP matrix */
				     mx,                 /* the HMM banded mx */
				     shmx,	         /* the HMM banded shadow matrix */
				     &Jb, &Lb, &Rb, &Tb, /* entry point for optimal J, L, R, T alignment */
				     &mode, &sc))        /* mode (J,L,R or T) and score of CYK parsetree */
       != eslOK) return status; 
    if     (mode == TRMODE_J) v = Jb;
    else if(mode == TRMODE_L) v = Lb;
    else if(mode == TRMODE_R) v = Rb;
    else if(mode == TRMODE_T) v = Tb;
    else ESL_FAIL(eslEINVAL, errbuf, "cm_tr_alignT(), bogus initial alignment mode returned by TrCYK: %d\n", mode);
  }

  /* Allocations and initialization */
  pda_i = esl_stack_ICreate();
  pda_c = esl_stack_CCreate();
  if(pda_i == NULL) goto ERROR;
  if(pda_c == NULL) goto ERROR;
  j = L;
  i = 1;
  d = L;
  tr = CreateParsetree(100);
  InsertTraceNodewithMode(tr, -1, TRACE_LEFT_CHILD, i, j, v, mode);

  while (1) {
    if(cm->sttype[v] != EL_st) printf("v: %4d  mode: %4d  j: %4d d: %4d\n", v, mode, j, d);
    else                       printf("v: %4d  mode: %4d  j: %4d d: %4d EL\n", v, mode, j, d);

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
l       * traceback altogether. This is the only way to break the
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
      if     (mode == TRMODE_J) yoffset = shmx->Jyshadow[v][j][d];
      else if(mode == TRMODE_L) yoffset = shmx->Lyshadow[v][j][d];
      else if(mode == TRMODE_R) yoffset = shmx->Ryshadow[v][j][d];
      else if(mode == TRMODE_T) ESL_FAIL(eslEINVAL, errbuf, "truncation mode T for non B states");
      else                      ESL_FAIL(eslEINVAL, errbuf, "bogus truncation mode %d\n", mode);
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
#if 0
      else if (yoffset == USED_TRUNC_BEGIN) 
	{ /* truncated begin; can only happen once, from root */
	  mode = nxtmode; /* this was set as bmode above */
	  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b, mode);
	  v = b;
	}
#endif
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
	  v = y;
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

/* Function: cm_tr_alignT_hb()
 * Date:     EPN, Thu Sep  8 07:59:10 2011
 *           EPN 03.29.06 (cm_alignT_hb()
 *
 * Purpose: Call either cm_TrCYKInsideAlignHB() (if !<do_optacc>), or
 *           cm_TrOptAccAlignHB() (if <do_optacc>), get vjd shadow
 *           matrix; then trace back and append to an existing but
 *           empty parsetree tr.  The full sequence 1..L will be
 *           aligned.
 *        
 *           If <do_optacc>==TRUE then emit_mx must != NULL and
 *           <opt_mode> must not be TRMODE_UNKNOWN, it will have been
 *           determined by caller from a cm_TrInsideAlignHB() call and
 *           passed in. If <do_optacc> is FALSE, then we're doing CYK
 *           alignment and we don't know the truncation mode of the
 *           alignment mode yet, so <opt_mode> will be passed in as
 *           TRMODE_UNKNOWN.
 *
 * Args:     cm         - the model 
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitized sequence [1..L]   
 *           L          - length of the dsq to align
 *           size_limit - max size in Mb for DP matrix
 *           do_optacc  - TRUE to align with optimal accuracy, else use CYK
 *           mx         - the DP matrix to fill in
 *           shmx       - the shadow matrix to fill in
 *           emit_mx    - the pre-filled emit matrix, must be non-NULL if do_optacc
 *           opt_mode   - the optimal alignment mode, if unknown TRMODE_UNKNOWN is passed
 *           ret_tr     - RETURN: the optimal parsetree
 *           ret_sc     - RETURN: optimal score (CYK if !do_optacc, else avg PP of all 1..L residues) 
 *
 * Returns:  <eslOK>     on success.
 * Throws:   <eslERANGE> if required DP matrix size exceeds <size_limit>, in 
 *                       this case, alignment has been aborted, ret_* variables are not valid
 *           <eslEINVAL> on traceback problem: bogus state
 */
int
cm_tr_alignT_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, 
		CM_TR_HB_EMIT_MX *emit_mx, char opt_mode, Parsetree_t **ret_tr, float *ret_sc)
{
  int       status;         
  Parsetree_t *tr = NULL;       /* the parsetree */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda_i;             /* stack that tracks bifurc parent of a right start */
  ESL_STACK *pda_c;             /* stack that tracks mode of bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, seq positions */
  int       k;			/* right subtree len for bifurcs */
  int       y, yoffset;         /* child state y, it's offset */
  int       bifparent;          /* B_st parent */
  int       jp_v;               /* j-jmin[v] for current j, and current v */
  int       dp_v;               /* d-hdmin[v][jp_v] for current j, current v, current d*/
  char      mode;               /* current truncation mode: TRMODE_J | TRMODE_L | TRMODE_R | TRMODE_T */
  char      prvmode, nxtmode;   /* previous, next truncation mode */
  int       allow_S_trunc_end;  /* set to true to allow d==0 BEGL_S and BEGR_S truncated ends */
  int       Jb, Lb, Rb, Tb;     /* entry states for best J,L,R,T alignment */

  /* pointers to cp9b data for convenience */
  CP9Bands_t  *cp9b = cm->cp9b;
  int         *jmin = cp9b->jmin;
  int         *jmax = cp9b->jmax;
  int       **hdmin = cp9b->hdmin;
  int       **hdmax = cp9b->hdmax;

  if (cp9b->jmin[0]             > L || cp9b->jmax[0]             < L) ESL_FAIL(eslEINVAL, errbuf, "cm_tr_alignT_hb(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, cp9b->jmin[0], cp9b->jmax[0]);
  if (cp9b->hdmin[0][L-jmin[0]] > L || cp9b->hdmax[0][L-jmin[0]] < L) ESL_FAIL(eslEINVAL, errbuf, "cm_tr_alignT_hb(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, cp9b->hdmin[0][L-jmin[0]], cp9b->hdmax[0][L-jmin[0]]);

  if(do_optacc) {
    cm_Fail("cm_TrOptAccAlignHB() not yet implemented");
#if 0
    if((status = cm_TrOptAccAlignHB(cm, errbuf, dsq, L,
				  size_limit,   /* max size of DP matrix */
				  opt_mode,     /* marginal mode of optimal alignment */
				  mx,	        /* the DP matrix, to expand and fill-in */
				  shmx,	        /* the shadow matrix, to expand and fill-in */
				  emit_mx,      /* pre-calc'ed emit matrix */
				  &b,           /* the entry point for optimal alignment */
				  &sc))         /* avg post prob of all emissions in optimally accurate parsetree */
       != eslOK) return status;
    /* set root state of entire parse */
    v    = b;
    mode = opt_mode;
#endif
  }
  else {
    if((status = cm_TrCYKInsideAlignHB(cm, errbuf, dsq, L,
				       size_limit,         /* max size of DP matrix */
				       mx,                 /* the HMM banded mx */
				       shmx,	         /* the HMM banded shadow matrix */
				       &Jb, &Lb, &Rb, &Tb, /* entry point for optimal J, L, R, T alignment */
				       &mode, &sc))        /* mode (J,L,R or T) and score of CYK parsetree */
	!= eslOK) return status;
    if     (mode == TRMODE_J) v = Jb;
    else if(mode == TRMODE_L) v = Lb;
    else if(mode == TRMODE_R) v = Rb;
    else if(mode == TRMODE_T) v = Tb;
    else ESL_FAIL(eslEINVAL, errbuf, "cm_tr_alignT_hb(), bogus initial alignment mode returned by TrCYKHB: %d\n", mode);
  }

  /* Allocations and initialization */
  pda_i = esl_stack_ICreate();
  pda_c = esl_stack_CCreate();
  if(pda_i == NULL) goto ERROR;
  if(pda_c == NULL) goto ERROR;
  j = L;
  i = 1;
  d = L;
  tr = CreateParsetree(100);
  InsertTraceNodewithMode(tr, -1, TRACE_LEFT_CHILD, i, j, v, mode);

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
#if 0
      else if (yoffset == USED_TRUNC_BEGIN) 
	{ /* truncated begin; can only happen once, from root */
	  mode = nxtmode; /* this was set as bmode above */
	  InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b, mode);
	  v    = b;
	}
#endif
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


/* Function: cm_TrAlign()
 * Incept:   EPN, Sat Sep 10 12:58:09 2011
 * 
 * Purpose: Wrapper for the cm_tr_alignTb() routine - solve a full
 *           alignment problem using trCYK, truncated optimal accuracy
 *           or sampling and return the traceback and the score.
 *           
 *           Identical to cm_TrAlignHB() but HMM bands are not used here.
 * 
 *           Input arguments allow this function to be run in 6 'modes':
 *
 *           mode      returns                 arguments
 *           ----  ----------------  ----------------------------------------
 *                 tr        ppstrs  do_optacc  do_sample post_mx   ret_ppstr
 *                 ----------------  ----------------------------------------
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
 *           do_optacc - TRUE to not do CYK alignment, determine the Holmes/Durbin optimally 
 *                       accurate parsetree in ret_tr, requires post_mx != NULL
 *           do_sample - TRUE to sample a parsetree from the Inside matrix
 *           mx        - the main dp matrix, only cells within bands in cm->cp9b will be valid. 
 *           shmx      - the HMM banded shadow matrix to fill in and traceback, same cells as mx are valid.
 *           post_mx   - dp matrix for posterior calculation, can be NULL only if !do_optacc
 *           emit_mx   - emit matrix to fill
 *           r         - source of randomness, must be non-NULL only if do_sample==TRUE
 *           ret_ppstr - RETURN: posterior code 1, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_ins_sc- RETURN: if(do_optacc || ret_ppstr != NULL): inside score of sequence in bits
 *                               else: should be NULL (inside will not be run)
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_sc    - RETURN: if(!do_optacc): score of the alignment in bits.
 *                               if( do_optacc): avg PP of all L aligned residues in optacc parsetree
 * 
 * Returns: <eslOK> on success.
 * 
 * Throws:  <eslEINVAL> on contract violation
 *          <eslERANGE> if required CM_TR_MX for Inside/Outside/CYK/Posterior exceeds <size_limit>
 */
int
cm_TrAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, int do_sample,
	   CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, CM_TR_MX *post_mx, CM_TR_EMIT_MX *emit_mx, ESL_RANDOMNESS *r, 
	   char **ret_ppstr, float *ret_ins_sc, Parsetree_t **ret_tr, float *ret_sc)
{
  int          status;
  Parsetree_t *tr = NULL;
  float        sc;
  float        ins_sc; /* inside score */
  int          do_post;
  char        *ppstr;
  int          have_ppstr;
  char         opt_mode; /* the optimal mode */

  have_ppstr = (ret_ppstr != NULL)       ? TRUE : FALSE;
  do_post    = (do_optacc || have_ppstr) ? TRUE : FALSE;

  /* Contract check */
  if(do_optacc && do_sample)         ESL_FAIL(eslEINCOMPAT, errbuf, "cm_TrAlign(), do_optacc and do_sample are both TRUE.");
  if(do_optacc && post_mx == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_TrAlign(), do_optacc is TRUE, but post_mx == NULL.\n");
  if(do_sample && r       == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_TrAlign(), do_sample but r is NULL.");

  /* if do_post:   fill Inside, Outside, Posterior matrices, in that order.
   * if do_sample: fill Inside and sample from it.
   */
  if(do_post || do_sample) { 
    if((status = cm_TrInsideAlign(cm, errbuf, dsq, L, size_limit, mx, &opt_mode, &ins_sc)) != eslOK) return status;
    if(do_sample) { 
      cm_Fail("ERROR, do_sample not yet implemented");
      //if((status = SampleFromTrInside(r, cm, errbuf, dsq, L, mx, &tr, &sc)) != eslOK) return status; 
    }
    if(do_post) { /* Inside was called above, now do Outside, then Posterior */
      if((status = cm_TrOutsideAlign(cm, errbuf, dsq, L, size_limit, (cm->align_opts & CM_ALIGN_CHECKINOUT), opt_mode, post_mx, mx)) != eslOK) return status;
      if((status = cm_TrPosterior(cm, errbuf, L, size_limit, opt_mode, mx, post_mx, post_mx)) != eslOK) return status;   
      if((status = cm_TrEmitterPosterior(cm, errbuf, L, size_limit, post_mx, emit_mx, opt_mode, (cm->align_opts & CM_ALIGN_CHECKINOUT))) != eslOK) return status;   
    }
  }

  if(!do_sample) { /* if do_sample, we already have a parsetree */
    if((status = cm_tr_alignT(cm, errbuf, dsq, L, size_limit, do_optacc, mx, shmx, emit_mx, opt_mode, &tr, &sc)) != eslOK) return status;
    float parsetree_sc, parsetree_struct_sc;
    ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
    ParsetreeScore(cm, NULL, NULL, tr, dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
    printf("Parsetree score      : %.4f           (FULL LENGTH OPTACC)\n", parsetree_sc);
  }

  if(have_ppstr || do_optacc) { /* call cm_PostCode to get average PP and optionally a PP string (if have_ppstr) */
    if((status = cm_TrPostCode(cm, errbuf, L, emit_mx, tr, 
			       (have_ppstr ? &ppstr : NULL), 
			       (do_optacc  ? &sc    : NULL))) != eslOK) return status;
  }

  if (ret_ppstr  != NULL) *ret_ppstr  = ppstr; else free(ppstr);
  if (ret_tr     != NULL) *ret_tr     = tr;    else FreeParsetree(tr);
  if (ret_ins_sc != NULL) *ret_ins_sc = ins_sc; 
  if (ret_sc     != NULL) *ret_sc     = sc;

  ESL_DPRINTF1(("returning from cm_TrAlign() sc : %f\n", sc)); 
  return eslOK;
}

/* Function: cm_TrAlignHB()
 * Incept:   EPN, Thu Sep  8 08:55:26 2011
 *           EPN, Fri Oct 26 09:31:43 2007 [FastAlignHB()]
 * 
 * Purpose: Wrapper for the cm_tr_alignT_hb() routine - solve a full
 *           alignment problem using trCYK, truncated optimal accuracy
 *           or sampling and return the traceback and the score,
 *           without dividing & conquering, but by using bands on the
 *           j and d dimensions of the DP matrix.  Bands derived by
 *           HMM Forward/Backward runs. Optionally return a posterior
 *           code string.
 *           
 *           Identical to cm_TrAlign() but HMM bands are used here.
 *           See that function's 'Purpose' for more details.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence, 1..L
 *           L         - length of sequence 
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_optacc - TRUE to not do CYK alignment, determine the Holmes/Durbin optimally 
 *                       accurate parsetree in ret_tr, requires post_mx != NULL
 *           do_sample - TRUE to sample a parsetree from the Inside matrix
 *           mx        - the main dp matrix, only cells within bands in cm->cp9b will be valid. 
 *           shmx      - the HMM banded shadow matrix to fill in and traceback, same cells as mx are valid.
 *           post_mx   - dp matrix for posterior calculation, can be NULL only if !do_optacc
 *           emit_mx   - emit matrix to fill
 *           r         - source of randomness, must be non-NULL only if do_sample==TRUE
 *           ret_ppstr - RETURN: posterior code 1, (pass NULL if not wanted, must be NULL if post_mx == NULL)
 *           ret_ins_sc- RETURN: if(do_optacc || ret_ppstr != NULL): inside score of sequence in bits
 *                               else: should be NULL (inside will not be run)
 *           ret_tr    - RETURN: traceback (pass NULL if trace isn't wanted)
 *           ret_sc    - RETURN: if(!do_optacc): score of the alignment in bits.
 *                               if( do_optacc): avg PP of all L aligned residues in optacc parsetree
 * 
 * Returns: <ret_tr>, <ret_ppstr>, <ret_sc>, see 'Args' section
 * 
 * Returns: <eslOK> on success
 * 
 * Throws:  <eslEINVAL> on contract violation
 *          <eslERANGE> if required CM_HB_MX for Inside/Outside/CYK/Posterior exceeds <size_limit>
 */
int
cm_TrAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, int do_sample,
	     CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, CM_TR_HB_MX *post_mx, CM_TR_HB_EMIT_MX *emit_mx, ESL_RANDOMNESS *r, 
	     char **ret_ppstr, float *ret_ins_sc, Parsetree_t **ret_tr, float *ret_sc)
{
  int          status;
  Parsetree_t *tr = NULL;
  float        sc     = IMPOSSIBLE;
  float        ins_sc = IMPOSSIBLE; /* inside score */
  int          do_post;
  char        *ppstr;
  int          have_ppstr;
  char         opt_mode = TRMODE_UNKNOWN;

  have_ppstr = (ret_ppstr != NULL) ? TRUE : FALSE;
  do_post = (do_optacc || have_ppstr) ? TRUE : FALSE;

  /* Contract check */
  if(do_optacc && do_sample)         ESL_FAIL(eslEINCOMPAT, errbuf, "cm_TrAlignHB(), do_optacc and do_sample are both TRUE.");
  if(do_optacc && post_mx == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_TrAlignHB(), do_optacc is TRUE, but post_mx == NULL.\n");
  if(do_sample && r       == NULL)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_TrAlignHB(), do_sample but r is NULL.");

  /* if do_post, fill Inside, Outside, Posterior matrices, in that order */
  /* if do_sample (and !do_post) fill Inside and sample from it */
  if(do_post || do_sample) { 
    cm_Fail("ERROR, do_post || do_sample not yet implemented");
#if 0 
    if((status = cm_TrInsideAlignHB (cm, errbuf, dsq, L, size_limit, mx, &opt_mode, &ins_sc)) != eslOK) return status;
    if(do_sample) { 
      if((status = cm_SampleParsetreeHB(cm, errbuf, dsq, L, mx, r, &tr, &sc)) != eslOK) return status; 
    }
    if(do_post) { /* Inside was called above, now do Outside, then Posterior */
      if((status = cm_OutsideAlignHB(cm, errbuf, dsq, L, size_limit, (cm->align_opts & CM_ALIGN_CHECKINOUT), opt_mode, post_mx, mx, NULL)) != eslOK) return status;
      if((status = cm_TrPosteriorHB (cm, errbuf, L, size_limit, opt_mode, mx, post_mx, post_mx)) != eslOK) return status;   
      if((status = cm_TrEmitterPosteriorHB(cm, errbuf, L, size_limit, post_mx, emit_mx, opt_mode, (cm->align_opts & CM_ALIGN_CHECKINOUT))) != eslOK) return status;   
    }
#endif
  }

  if(!do_sample) { /* if do_sample, we already have a parsetree */
    if((status = cm_tr_alignT_hb(cm, errbuf, dsq, L, size_limit, do_optacc, mx, shmx, emit_mx, opt_mode, &tr, &sc)) != eslOK) return status;
    float parsetree_sc, parsetree_struct_sc;
    ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
    ParsetreeScore(cm, NULL, NULL, tr, dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
    printf("Parsetree score      : %.4f           (FULL LENGTH OPTACC)\n", parsetree_sc);
  }

  if(have_ppstr || do_optacc) {
    cm_Fail("ERROR, have_ppstr not yet implemented");
#if 0
    if((status = cm_TrPostCodeHB(cm, errbuf, L, emit_mx, tr, 
				 (have_ppstr ? &ppstr : NULL), 
				 (do_optacc  ? &sc    : NULL))) != eslOK) return status;
#endif
  }

  if (ret_ppstr  != NULL) *ret_ppstr  = ppstr; else free(ppstr);
  if (ret_tr     != NULL) *ret_tr     = tr;    else FreeParsetree(tr);
  if (ret_ins_sc != NULL) *ret_ins_sc = ins_sc; 
  if (ret_sc     != NULL) *ret_sc     = sc;

  ESL_DPRINTF1(("returning from cm_TrAlignHB() sc : %f\n", sc)); 
  return eslOK;
}

/* Function: cm_TrCYKInsideAlign()
 * based on cm_CYKInsideAlign()
 *
 * Date:     EPN, Fri Sep  9 15:35:06 2011
 *
 * Note:     Very similar to inside(), but slightly more efficient.
 *           Identical to cm_TrCYKInsideAlignHB() but HMM bands are not
 *           used.
 *
 * Purpose:  Perform trCYK alignment on a full sequence 1..L
 *           rooted at state 0. Very similar to cm_CYKInsideAlign()
 *           except we're doing truncated alignment and marginal 
 *           alignment modes are possible.
 *
 *           In truncated alignment, standard local begins are not
 *           possible.  We return the best internal entry state for
 *           each marginal mode in <ret_{J,L,R,T}b>. The optimal score
 *           for each marginal mode is in:
 *           {J,L,R,T}alpha[0][L][L]. The optimal mode (the mode that
 *           gives the max score in {J,L,R,T}alpha[0][L][L] is
 *           returned in <ret_mode>, and that max score is returned in
 *           <ret_sc>.
 *
 * Args:     cm          - the model    [0..M-1]
 *           errbuf      - char buffer for reporting errors
 *           dsq         - the digitaized sequence [1..L]   
 *           L           - length of target sequence, we align 1..L
 *           size_limit  - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx          - dp matrix 
 *           shmx        - shadow matrix
 *           ret_Jb      - best internal entry state for Joint    marginal mode
 *           ret_Lb      - best internal entry state for Left     marginal mode
 *           ret_Rb      - best internal entry state for Right    marginal mode
 *           ret_Tb      - best internal entry state for Terminal marginal mode (-1 if CM has 0 B states)
 *                         mx->{J,L,R,T}alpha[ret_{J,L,R,T}b][L][L] are the optimal {J,L,R,T} alignment scores
 *                         these are equal to the corresponding mx->{J,L,R,T}alpha[0][L][L] values.
 *           ret_mode    - mode of optimal CYK parsetree (TRMODE_J | TRMODE_L | TRMODE_R | TRMODE_T)
 *           ret_sc      - score of optimal, CYK parsetree in any mode (max of mx->{J,L,R,T}alpha[0][L][L])
 *                       
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslERANGE> if required mx or shmx size exceeds <size_limit>
 *           In this case alignment has been aborted, <ret_*> variables are not valid
 */
int
cm_TrCYKInsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
		    int *ret_Jb, int *ret_Lb, int *ret_Rb, int *ret_Tb, char *ret_mode, float *ret_sc)
{
  int      status;          /* easel status code */
  int      v,y,z;	    /* indices for states  */
  int      j,d,i,k;	    /* indices in sequence dimensions */
  float    sc;		    /* a temporary variable holding a score */
  int      yoffset;	    /* y=base+offset -- counter in child states that v can transit to */
  float   *el_scA;          /* [0..d..W-1] probability of local end emissions of length d */
  int      sd;              /* StateDelta(cm->sttype[v]) */
  int      sdl;             /* StateLeftDelta(cm->sttype[v] */
  int      sdr;             /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;           /* j - sdr */
  int      d_sd;            /* d - sd */
  int      d_sdl;           /* d - sdl */
  int      d_sdr;           /* d - sdr */
  float    tsc;             /* a transition score */

  /* other variables used in truncated version, but not standard version (not in cm_CYKAlign()) */
  float trunc_penalty = 0.; /* penalty in bits for a truncated hit */
  int   Jb, Lb, Rb, Tb;     /* state rooting {J,L,R,T} optimal parsetrees */
  int   mode = TRMODE_J;    /* truncation mode for obtaining optimal score <ret_sc> */
  int   have_el;            /* TRUE if local ends are on */
  int   Lyoffset0;          /* first yoffset to use for updating L matrix in IR/MR states, 1 if IR, 0 if MR */
  int   Ryoffset0;          /* first yoffset to use for updating R matrix in IL/ML states, 1 if IL, 0 if ML */

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
  Jb = Lb = Rb = Tb = -1;

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo       (cm, mx,   errbuf, L, size_limit)) != eslOK) return status;
  if((status = cm_tr_shadow_mx_GrowTo(cm, shmx, errbuf, L, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..L residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 
  if(shmx->Jy_ncells_valid > 0) for(i = 0; i < shmx->Jy_ncells_valid; i++) shmx->Jyshadow_mem[i] = USED_EL;
  if(shmx->Ly_ncells_valid > 0) for(i = 0; i < shmx->Ly_ncells_valid; i++) shmx->Lyshadow_mem[i] = USED_TRUNC_END;
  if(shmx->Ry_ncells_valid > 0) for(i = 0; i < shmx->Ry_ncells_valid; i++) shmx->Ryshadow_mem[i] = USED_TRUNC_END;
  /* for B states, shadow matrix holds k, length of right fragment, this will almost certainly be overwritten */
  if(shmx->Jk_ncells_valid > 0) esl_vec_ISet(shmx->Jkshadow_mem, shmx->Jk_ncells_valid, 0);
  if(shmx->Lk_ncells_valid > 0) esl_vec_ISet(shmx->Lkshadow_mem, shmx->Lk_ncells_valid, 0);
  if(shmx->Rk_ncells_valid > 0) esl_vec_ISet(shmx->Rkshadow_mem, shmx->Rk_ncells_valid, 0);
  if(shmx->Tk_ncells_valid > 0) esl_vec_ISet(shmx->Tkshadow_mem, shmx->Tk_ncells_valid, 0);
  if(shmx->Lk_ncells_valid > 0) for(i = 0; i < shmx->Lk_ncells_valid; i++) shmx->Lkmode_mem[i] = TRMODE_J;
  if(shmx->Rk_ncells_valid > 0) for(i = 0; i < shmx->Rk_ncells_valid; i++) shmx->Rkmode_mem[i] = TRMODE_J;

  /* if local ends are on, replace the EL deck IMPOSSIBLEs with EL scores */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el) { 
    for (j = 0; j <= L; j++) {
      for (d = 0;  d <= j; d++) Jalpha[cm->M][j][d] = el_scA[d];
    }
  }

  /* Main recursion */
  for (v = cm->M-1; v > 0; v--) { /* almost to ROOT, ROOT is special because all hits must use truncated begins */
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
	for (d = sd; d <= j; d++) { 
	  Jalpha[v][j][d] = el_scA[d-sd] + cm->endsc[v];
	  /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (j = 0; j <= L; j++) {
	Jalpha[v][j][0] = 0.;
	Lalpha[v][j][0] = 0.;
	Ralpha[v][j][0] = 0.;
	/* rest of deck remains IMPOSSIBLE */
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

      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
	Ryoffset0 = cm->sttype[v] == IL_st ? 1 : 0; /* don't allow IL self transits in R mode */
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++) {
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
	    for (yoffset = Ryoffset0; yoffset < cm->cnum[v]; yoffset++) { /* using Ryoffset0 instead of 0 disallows IL self transits in R mode */
	      y = cm->cfirst[v] + yoffset; 
	      if ((sc = Jalpha[y][j_sdr][d] + tsc_v[yoffset]) > Ralpha[v][j][d]) { 
		Ralpha[v][j][d] = sc; 
		Ryshadow[v][j][d]= yoffset + TRMODE_J_OFFSET;
	      }
	      if ((sc = Ralpha[y][j_sdr][d] + tsc_v[yoffset]) > Ralpha[v][j][d]) { 
		Ralpha[v][j][d] = sc;
		Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	      }
	    }
	    Ralpha[v][j][d] = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
	  }
	}
      } /* end of if(! StateIsDetached(cm,v )) */
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

      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
	Lyoffset0 = cm->sttype[v] == IR_st ? 1 : 0; /* don't allow IR self transits in L mode */
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++) {
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
	    for (yoffset = Lyoffset0; yoffset < cm->cnum[v]; yoffset++) { /* using Lyoffset0, instead of 0 disallows IR self transits in L mode */
	      y = cm->cfirst[v] + yoffset; 
	      if ((sc = Jalpha[y][j][d] + tsc_v[yoffset]) > Lalpha[v][j][d]) { 
		Lalpha[v][j][d] = sc;
		Lyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	      }
	      if ((sc = Lalpha[y][j][d] + tsc_v[yoffset]) > Lalpha[v][j][d]) { 
		Lalpha[v][j][d] = sc;
		Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	      }
	    }
	    Lalpha[v][j][d] = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	  }
	}
      } /* end of if(! StateIsDetached(cm, v)) */
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

	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;

	  for (d = sd; d <= j; d++) { /* sd == 2 for MP state */
	    d_sd = d-sd;
	    if((sc = Jalpha[y][j_sdr][d_sd] + tsc) > Jalpha[v][j][d]) {
	      Jalpha[v][j][d]   = sc;
	      Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	  }
	  /* note we use 'j' and 'd_sdl' not 'j_sdr' for 'd_sd' for L, plus minimum d is sdl (1) */
	  for (d = sdl; d <= j; d++) { /* sdl == 1 for MP state */
	    d_sdl = d-sdl;
	    if((sc = Jalpha[y][j][d_sdl] + tsc) > Lalpha[v][j][d]) {
	      Lalpha[v][j][d]   = sc;
	      Lyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if((sc = Lalpha[y][j][d_sdl] + tsc) > Lalpha[v][j][d]) {
	      Lalpha[v][j][d]   = sc;
	      Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	    }
	  }
	  /* note we use 'd_sdr' not 'd_sd' for R, plus minimum d is sdr (1) */
	  for (d = sdr; d <= j; d++) { /* sdr == 1 for MP state */
	    d_sdr = d - sdr;
	    if((sc = Jalpha[y][j_sdr][d_sdr] + tsc) > Ralpha[v][j][d]) {
	      Ralpha[v][j][d]   = sc;
	      Ryshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	    if((sc = Ralpha[y][j_sdr][d_sdr] + tsc) > Ralpha[v][j][d]) {
	      Ralpha[v][j][d]   = sc;
	      Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	    }
	  }
	}
      }
      /* add in emission score */
      for (j = 0; j <= L; j++) {
	i = j;
	Jalpha[v][j][1] = IMPOSSIBLE;
	Lalpha[v][j][1] = lmesc_v[dsq[i]];
	Lyshadow[v][j][1] = USED_TRUNC_END;
	Ralpha[v][j][1] = rmesc_v[dsq[j]];
	Ryshadow[v][j][1] = USED_TRUNC_END;
	i--;
	for (d = 2; d <= j; d++) {
	  Jalpha[v][j][d] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	  Lalpha[v][j][d] += lmesc_v[dsq[i]];
	  Ralpha[v][j][d] += rmesc_v[dsq[j]];
	  i--;
	}
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = 0; j <= L; j++) {
	for (d = 1; d <= j; d++) {
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
	
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  
	  for (d = sd; d <= j; d++) {
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

      for (j = 0; j <= L; j++) {
	for (d = 0; d <= j; d++) {
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

    /* Now handle truncated begin transitions from ROOT_S, state 0 */
    /* Standard local begins are not allowed in truncated mode (they're kind of like J alignments though) */
    /* check if we have a new optimally scoring Joint alignment in J matrix (much like a standard local begin) */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st
       || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
      sc = Jalpha[v][L][L] + trunc_penalty;
      if (sc > Jalpha[0][L][L]) { 
	Jalpha[0][L][L] = sc;
	Jb = v;
      }
      /* check if we have a new optimally scoring Left alignment in L matrix */
      if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) { 
	sc = Lalpha[v][L][L] + trunc_penalty;
	if (sc > Lalpha[0][L][L]) { 
	  Lalpha[0][L][L] = sc;
	  Lb = v;
	}
      }	    
      /* check if we have a new optimally scoring Right alignment in R matrix */
      if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
	sc = Ralpha[v][L][L] + trunc_penalty;
	if (sc > Ralpha[0][L][L]) { 
	  Ralpha[0][L][L] = sc;
	  Rb = v;
	}
      }	    
      /* check if we have a new optimally scoring Terminal alignment in T matrix */
      if(cm->sttype[v] == B_st) { 
	sc = Talpha[v][L][L] + trunc_penalty;
	if (sc > Talpha[0][L][L]) { 
	  Talpha[0][L][L] = sc;
	  Tb = v;
	}
      }	    
    }
  } /* end loop for (v = cm->M-1; v > 0; v--) */
  FILE *fp1; fp1 = fopen("tmp.trcykmx", "w");   cm_tr_mx_Dump(fp1, mx); fclose(fp1);
  FILE *fp2; fp2 = fopen("tmp.trcykshmx", "w"); cm_tr_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2);

  sc   = Jalpha[0][L][L];
  mode = TRMODE_J;
  if (Lalpha[0][L][L] > sc) { 
    sc   = Lalpha[0][L][L];
    mode = TRMODE_L;
  }
  if (Ralpha[0][L][L] > sc) { 
    sc   = Ralpha[0][L][L];
    mode = TRMODE_R;
  }
  if (Talpha[0][L][L] > sc) { 
    sc   = Talpha[0][L][L];
    mode = TRMODE_T;
  }

  if(ret_Jb   != NULL) *ret_Jb   = Jb;    
  if(ret_Lb   != NULL) *ret_Lb   = Lb;    
  if(ret_Rb   != NULL) *ret_Rb   = Rb;    
  if(ret_Tb   != NULL) *ret_Tb   = Tb;    
  if(ret_mode != NULL) *ret_mode = mode;    
  if(ret_sc   != NULL) *ret_sc   = sc;

  free(el_scA);

  ESL_DPRINTF1(("cm_TrCYKInsideAlign return sc: %f\n", sc));
  printf("cm_TrCYKInsideAlign return sc: %f\n", sc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_TrCYKInsideAlignHB()
 *
 * Date:     EPN, Wed Sep  7 12:13:43 2011
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
 *           In truncated alignment, standard local begins are not
 *           possible.  We return the best internal entry state for
 *           each marginal mode in <ret_{J,L,R,T}b>. The optimal score
 *           for each marginal mode is in:
 *           {J,L,R,T}alpha[0][L][L]. The optimal mode (the mode that
 *           gives the max score in {J,L,R,T}alpha[0][L][L] is
 *           returned in <ret_mode>, and that max score is returned in
 *           <ret_sc>.
 *
 * Args:     cm         - the model    [0..M-1]
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitaized sequence [1..L]   
 *           L          - length of target sequence, we align 1..L
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx         - the dp matrix, only cells within bands in cm->cp9b will be valid. 
 *           shmx       - the HMM banded shadow matrix to fill in, only cells within bands are valid
 *           ret_Jb     - best internal entry state for Joint    marginal mode
 *           ret_Lb     - best internal entry state for Left     marginal mode
 *           ret_Rb     - best internal entry state for Right    marginal mode
 *           ret_Tb     - best internal entry state for Terminal marginal mode (-1 if CM has 0 B states)
 *                        mx->{J,L,R,T}alpha[ret_{J,L,R,T}b][L][L] are the optimal {J,L,R,T} alignment scores
 *                        these are equal to the corresponding mx->{J,L,R,T}alpha[0][L][L] values.
 *           ret_mode   - mode of optimal CYK parsetree (TRMODE_J | TRMODE_L | TRMODE_R | TRMODE_T)
 *           ret_sc     - score of optimal, CYK parsetree in any mode (max of mx->{J,L,R,T}alpha[0][L][L])
 *                       
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslERANGE> if required mx or shmx size exceeds <size_limit>
 *           <eslEINVAL> if the full sequence is not within the bands for state 0
 *           In either case alignment has been aborted, ret_* variables are not valid
 */
int
cm_TrCYKInsideAlignHB(CM_t *cm, char *errbuf,  ESL_DSQ *dsq, int L, float size_limit, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, 
		      int *ret_Jb, int *ret_Lb, int *ret_Rb, int *ret_Tb, char *ret_mode, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;          /* temporary score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  int      sd;          /* StateDelta(cm->sttype[v]) */
  int      sdr;         /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;       /* j - sdr */

  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y, dp_z;   /* d index for state v/y/z in alpha */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dp_y_sdr;           /* dp_y - sdr */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Lp;                 /* L index also changes depending on state */
  float    tsc;                /* a transition score */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  int      jp_0;               /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;               /* L offset in ROOT_S's (v==0) d band */

  /* variables related to truncated alignment (not in cm_CYKAlignHB() */
  float    trunc_penalty = 0.;     /* penalty in bits for a truncated hit */
  int      Jb, Lb, Rb, Tb;         /* state rooting {J,L,R,T} optimal parsetrees */
  int      mode = TRMODE_J;        /* truncation mode for obtaining optimal score <ret_sc> */
  int      do_J_v, do_J_y, do_J_z; /* is J matrix valid for state v, y, z? */
  int      do_L_v, do_L_y, do_L_z; /* is L matrix valid for state v, y, z? */
  int      do_R_v, do_R_y, do_R_z; /* is R matrix valid for state v, y, z? */
  int      do_T_v, do_T_y, do_T_z; /* is T matrix valid for state v, y, z? */

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
  Jb = Lb = Rb = Tb = -1;
  if (cp9b->jmin[0] > L || cp9b->jmax[0] < L)             ESL_FAIL(eslEINVAL, errbuf, "cm_TrCYKInsideAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, cp9b->jmin[0], cp9b->jmax[0]);
  jp_0 = L - jmin[0];
  if (cp9b->hdmin[0][jp_0] > L || cp9b->hdmax[0][jp_0] < L) ESL_FAIL(eslEINVAL, errbuf, "cm_TrCYKInsideAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, cp9b->hdmin[0][jp_0], cp9b->hdmax[0][jp_0]);
  Lp_0 = L - hdmin[0][jp_0];

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_hb_mx_GrowTo       (cm,   mx, errbuf, cp9b, L, size_limit)) != eslOK) return status;
  if((status = cm_tr_hb_shadow_mx_GrowTo(cm, shmx, errbuf, cp9b, L, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..L residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to IMPOSSIBLE, all cells of shadow matrix to USED_EL */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 
  if(shmx->Jy_ncells_valid > 0) for(i = 0; i < shmx->Jy_ncells_valid; i++) shmx->Jyshadow_mem[i] = USED_EL;
  if(shmx->Ly_ncells_valid > 0) for(i = 0; i < shmx->Ly_ncells_valid; i++) shmx->Lyshadow_mem[i] = USED_TRUNC_END;
  if(shmx->Ry_ncells_valid > 0) for(i = 0; i < shmx->Ry_ncells_valid; i++) shmx->Ryshadow_mem[i] = USED_TRUNC_END;
  /* for B states, shadow matrix holds k, length of right fragment, this will be overwritten */
  if(shmx->Jk_ncells_valid > 0) esl_vec_ISet(shmx->Jkshadow_mem, shmx->Jk_ncells_valid, 0);
  if(shmx->Lk_ncells_valid > 0) esl_vec_ISet(shmx->Lkshadow_mem, shmx->Lk_ncells_valid, 0);
  if(shmx->Rk_ncells_valid > 0) esl_vec_ISet(shmx->Rkshadow_mem, shmx->Rk_ncells_valid, 0);
  if(shmx->Tk_ncells_valid > 0) esl_vec_ISet(shmx->Tkshadow_mem, shmx->Tk_ncells_valid, 0);
  if(shmx->Lk_ncells_valid > 0) for(i = 0; i < shmx->Lk_ncells_valid; i++) shmx->Lkmode_mem[i] = TRMODE_J;
  if(shmx->Rk_ncells_valid > 0) for(i = 0; i < shmx->Rk_ncells_valid; i++) shmx->Rkmode_mem[i] = TRMODE_J;

  /* Note, at this point in the non-banded version (cm_TrCYKInsideAlign()) we
   * replace EL (cm->M) deck IMPOSSIBLEs with EL scores.  But we don't
   * here. It's actually not necessary there either, but it is done
   * for completeness and so if we check that matrix in combination
   * with an Outside CYK matrix, then scores in the EL deck are valid.
   * But we won't do that here, and we're concerned with efficiency,
   * so we don't waste time with resetting the EL deck.
   */

  /* Main recursion */
  for (v = cm->M-1; v > 0; v--) { /* almost to ROOT, ROOT is special because all hits must use truncated begins */
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
	  if(hdmin[v][jp_v] >= sd) { 
	    d    = hdmin[v][jp_v];
	    dp_v = 0;
	  }
	  else { 
	    d    = sd;
	    dp_v = sd - hdmin[v][jp_v];
	  }
	  for (; d <= hdmax[v][jp_v]; dp_v++, d++) {
	    if(d >= sd) { 
	      Jalpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
	      /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	    }
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
    else if(cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IL states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1] 
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */
      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
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
	      for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
		yoffset = yvalidA[yvalid_idx];
		y = cm->cfirst[v] + yoffset;
		do_R_y = cp9b->do_R[y];
		do_J_y = cp9b->do_J[y];
		if((do_J_y || do_R_y) && (y != v)) { /* (y != v) part is to disallow IL self transits in R mode */
		  jp_y_sdr = j - jmin[y] - sdr;
		  
		  /* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,L above) */
		  if((d) >= hdmin[y][jp_y_sdr] && (d) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		    dp_y = d - hdmin[y][jp_y_sdr];
		    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		    ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		    
		    if(do_J_y &&
		       ((sc = Jalpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]) > Ralpha[v][jp_v][dp_v])) { 
		      Ralpha[v][jp_v][dp_v] = sc;
		      Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		    }
		    if(do_R_y && 
		       ((sc = Ralpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]) > Ralpha[v][jp_v][dp_v])) { 
		      Ralpha[v][jp_v][dp_v] = sc;
		      Ryshadow[v][jp_v][dp_v] = yoffset + TRMODE_R_OFFSET;
		    }
		  }
		}
	      } /* end of for (yvalid_idx = 0... loop */
	    }
	  }
	}
      } /* end of if(! StateIsDetached(cm, v)) */
    }
    else if(cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) { 
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IR states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1].
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */

      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
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
	/* Handle L separately */
	if(do_L_v) { 
	  /* The second MR_st/IR_st 'for (j...' loop is for the L matrix which use a different set of j values */
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    jp_v = j - jmin[v];
	    yvalid_ct = 0;
	  
	    /* determine which children y we can legally transit to for v, j */
	    /* we use 'j' and not 'j_sdr' here for the L matrix, differently from J and R matrices above */
	    for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	      if(y != v       && /* y == v when yoffset == 0 && v is an IR state: we don't want to allow IR self transits in L mode */
		 j >= jmin[y] && j <= jmax[y]) yvalidA[yvalid_ct++] = yoffset; /* is j is valid for state y? */
	  
	    for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	      for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
		/* Note if we're an IL state, we can't self transit in R mode, this was ensured above when we set up yvalidA[] (xref:ELN3,p5)*/
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
		       (sc = Jalpha[y][jp_y][dp_y] + tsc_v[yoffset]) > Lalpha[v][jp_v][dp_v]) { 
		      Lalpha[v][jp_v][dp_v] = sc;
		      Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_J_OFFSET;
		    }
		    if(do_L_y &&
		       (sc = Lalpha[y][jp_y][dp_y] + tsc_v[yoffset]) > Lalpha[v][jp_v][dp_v]) { 
		      Lalpha[v][jp_v][dp_v] = sc;
		      Lyshadow[v][jp_v][dp_v] = yoffset + TRMODE_L_OFFSET;
		    }
		  }
		}
	      } /* end of for (yvalid_idx = 0... loop */
	    }
	  }
	}
      } /* end of if(! StateIsDetached(cm, v) */
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
           
    /* Now handle truncated begin transitions from ROOT_S, state 0 */
    if(L >= jmin[v] && L <= jmax[v]) { 
      jp_v = L - jmin[v];
      Lp   = L - hdmin[v][jp_v];
      if(L >= hdmin[v][jp_v] && L <= hdmax[v][jp_v]) {
	/* If we get here alpha[v][jp_v][Lp] and alpha[0][jp_0][Lp_0]
	 * are valid cells in the banded alpha matrix, corresponding to 
	 * alpha[v][L][L] and alpha[0][L][L] in the platonic matrix.
	 * (Le've already made sure alpha[0][jp_0][Lp_0] was valid 
	 * at the beginning of the function.)
	 *
	 * We don't allow standard local begins in truncated mode. 
	 */

	/* check if we have a new optimally scoring Joint alignment in J matrix (much like a standard local begin) */
	if(do_J_v && cp9b->do_J[0] && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st
		      || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)) { 
	  sc = Jalpha[v][jp_v][Lp] + trunc_penalty;
	  if (sc > Jalpha[0][jp_0][Lp_0]) { 
	    Jalpha[0][jp_0][Lp_0] = sc;
	    Jb = v;
	  }
	}
	/* check if we have a new optimally scoring Left alignment in L matrix */
	if(do_L_v && cp9b->do_L[0] && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st)) { 
	  sc = Lalpha[v][jp_v][Lp] + trunc_penalty;
	  if (sc > Lalpha[0][jp_0][Lp_0]) { 
	    Lalpha[0][jp_0][Lp_0] = sc;
	    Lb = v;
	  }
	}	    
	/* check if we have a new optimally scoring Right alignment in R matrix */
	if(do_R_v && cp9b->do_R[0] && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st)) { 
	  sc = Ralpha[v][jp_v][Lp] + trunc_penalty;
	  if (sc > Ralpha[0][jp_0][Lp_0]) { 
	    Ralpha[0][jp_0][Lp_0] = sc;
	    Rb = v;
	  }
	}	    
	/* check if we have a new optimally scoring Terminal alignment in T matrix */
	if(do_T_v && cp9b->do_T[0]) { 
	  sc = Talpha[v][jp_v][Lp] + trunc_penalty;
	  if (sc > Talpha[0][jp_0][Lp_0]) { 
	    Talpha[0][jp_0][Lp_0] = sc;
	    Tb = v;
	  }
	}	    
      }
    }
  } /* end loop for (v = cm->M-1; v > 0; v--) */
  /*FILE *fp1; fp1 = fopen("tmp.ahbmx", "w");   cm_tr_hb_mx_Dump(fp1, mx); fclose(fp1);*/
  /*FILE *fp2; fp2 = fopen("tmp.ahbshmx", "w"); cm_tr_hb_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2);*/

  sc   = IMPOSSIBLE;
  if (cp9b->do_J[0] && Jalpha[0][jp_0][Lp_0] > sc) { 
    sc   = Jalpha[0][jp_0][Lp_0];
    mode = TRMODE_J;
  }
  if (cp9b->do_L[0] && Lalpha[0][jp_0][Lp_0] > sc) { 
    sc   = Lalpha[0][jp_0][Lp_0];
    mode = TRMODE_L;
  }
  if (cp9b->do_R[0] && Ralpha[0][jp_0][Lp_0] > sc) { 
    sc   = Ralpha[0][jp_0][Lp_0];
    mode = TRMODE_R;
  }
  if (cp9b->do_T[0] && Talpha[0][jp_0][Lp_0] > sc) { 
    sc   = Talpha[0][jp_0][Lp_0];
    mode = TRMODE_T;
  }

  if(ret_Jb   != NULL) *ret_Jb   = Jb;    
  if(ret_Lb   != NULL) *ret_Lb   = Lb;    
  if(ret_Rb   != NULL) *ret_Rb   = Rb;    
  if(ret_Tb   != NULL) *ret_Tb   = Tb;    
  if(ret_mode != NULL) *ret_mode = mode;    
  if(ret_sc   != NULL) *ret_sc   = sc;

  free(el_scA);
  free(yvalidA);

  ESL_DPRINTF1(("cm_TrCYKInsideAlignHB return sc: %f\n", sc));
  printf("cm_TrCYKInsideAlignHB return sc: %.4f\n", sc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Function: cm_TrInsideAlign()
 * Date:     EPN, Mon Sep 12 04:31:43 2011
 *
 * Purpose: Run the truncated inside algorithm on a target sequence
 *          without using bands. The full target sequence 1..L is 
 *          aligned (only full alignments will contribute to the 
 *          Inside score).
 *     
 *          Identical to cm_InsideAlign() but no bands are used.
 * 
 *           Very similar to cm_TrCYKInsideAlign(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - we do TrInside, not TrCYK
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return bsc, b info about local begins 
 *
 *           This function complements cm_TrOutsideAlign().
 *
 * Args:     cm         - the model
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitized sequence
 *           L          - target sequence length
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx         - the dp matrix, grown and filled here
 *           ret_mode   - RETURN: mode of optimal truncation mode, TRMODE_{J,L,R,T} if {J,L,R,T}alpha[0][L][L] is max scoring.
 *           ret_sc     - RETURN: log P(S|M)/P(S|R), as a bit score
 *                        NOTE: we don't sum over different marginal modes, we pick the highest scoring
 *                        one (J,L,R or T) and return {J,L,R,T}alpha[0][L][L] the sum of all complete 
 *                        J,L,R, or T alignments.
 *
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslERANGE> if required CM_TR_MX size exceeds <size_limit>
 *           In this case alignment has been aborted, ret_sc is not valid
 */
int
cm_TrInsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, CM_TR_MX *mx, char *ret_mode, float *ret_sc)
{
  int      status;          /* easel status code */
  int      v,y,z;	    /* indices for states  */
  int      j,d,i,k;	    /* indices in sequence dimensions */
  float    sc;		    /* a temporary variable holding a score */
  int      yoffset;	    /* y=base+offset -- counter in child states that v can transit to */
  float   *el_scA;          /* [0..d..W-1] probability of local end emissions of length d */
  int      sd;              /* StateDelta(cm->sttype[v]) */
  int      sdl;             /* StateLeftDelta(cm->sttype[v] */
  int      sdr;             /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;           /* j - sdr */
  int      d_sd;            /* d - sd */
  int      d_sdl;           /* d - sdl */
  int      d_sdr;           /* d - sdr */
  float    tsc;             /* a transition score */
  int      have_el;         /* TRUE if local ends are on */

  /* other variables used in truncated version, but not standard version (not in cm_CYKAlign()) */
  float    trunc_penalty=0.;/* penalty in bits for a truncated hit */
  int      mode = TRMODE_J; /* truncation mode for obtaining optimal score <ret_sc> */
  int      Lyoffset0;       /* first yoffset to use for updating L matrix in IR/MR states, 1 if IR, 0 if MR */
  int      Ryoffset0;       /* first yoffset to use for updating R matrix in IL/ML states, 1 if IL, 0 if ML */

  /* the DP matrix */
  float ***Jalpha  = mx->Jdp; /* pointer to the Jalpha DP matrix */
  float ***Lalpha  = mx->Ldp; /* pointer to the Lalpha DP matrix */
  float ***Ralpha  = mx->Rdp; /* pointer to the Ralpha DP matrix */
  float ***Talpha  = mx->Tdp; /* pointer to the Talpha DP matrix */

  /* Allocations and initializations  */

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo       (cm, mx,   errbuf, L, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..L residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 

  /* if local ends are on, replace the EL deck IMPOSSIBLEs with EL scores */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el) { 
    for (j = 0; j <= L; j++) {
      for (d = 0;  d <= j; d++) Jalpha[cm->M][j][d] = el_scA[d];
    }
  }

  /* Main recursion */
  for (v = cm->M-1; v > 0; v--) { /* almost to ROOT, ROOT is special because all hits must use truncated begins */
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
	for (d = sd; d <= j; d++) { 
	  Jalpha[v][j][d] = el_scA[d-sd] + cm->endsc[v];
	  /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (j = 0; j <= L; j++) {
	Jalpha[v][j][0] = 0.;
	Lalpha[v][j][0] = 0.;
	Ralpha[v][j][0] = 0.;
	/* rest of deck remains IMPOSSIBLE */
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
      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
	Ryoffset0 = cm->sttype[v] == IL_st ? 1 : 0; /* don't allow IL self transits in R mode */
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++) {
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
	    for (yoffset = Ryoffset0; yoffset < cm->cnum[v]; yoffset++) { /* using Ryoffset0 instead of 0 disallows IL self transits in R mode */
	      y = cm->cfirst[v] + yoffset; 
	      Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Jalpha[y][j_sdr][d] + tsc_v[yoffset]);
	      Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[y][j_sdr][d] + tsc_v[yoffset]);
	    }
	    Ralpha[v][j][d] = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
	  }
	}
      } /* end of if(! StateIsDetached(cm, v)) */
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
      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
	Lyoffset0 = cm->sttype[v] == IR_st ? 1 : 0; /* don't allow IR self transits in L mode */
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++) {
	    d_sd = d - sd;
	    i = j - d + 1;
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	      y = cm->cfirst[v] + yoffset; 
	      Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	      Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	    }
	    
	    Jalpha[v][j][d] += esc_v[dsq[j]];
	    Ralpha[v][j][d]  = (d >= 2) ? Ralpha[v][j][d] + esc_v[dsq[j]] : esc_v[dsq[j]];
	    
	    Jalpha[v][j][d]  = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	    Ralpha[v][j][d]  = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
	    
	    /* handle L separately */
	    /* note we use 'j' and 'd', not 'j_sdr' and 'd_sd' (which we used in the corresponding loop for J,R above) */
	    for (yoffset = Lyoffset0; yoffset < cm->cnum[v]; yoffset++) { /* using Lyoffset0, instead of 0 disallows IR self transits in L mode */
	      y = cm->cfirst[v] + yoffset; 
	      Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Jalpha[y][j][d] + tsc_v[yoffset]);
	      Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[y][j][d] + tsc_v[yoffset]);
	    }
	    Lalpha[v][j][d] = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	  }
	}
      } /* end of if(! StateIsDetached(cm, v)) */
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

	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;

	  for (d = sd; d <= j; d++) { /* sd == 2 for MP state */
	    d_sd = d - sd;
	    Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j_sdr][d_sd] + tsc_v[yoffset]);
	  }
	  /* note we use 'j' and 'd_sdl' not 'j_sdr' for 'd_sd' for L, plus minimum d is sdl (1) */
	  for (d = sdl; d <= j; d++) { /* sdl == 1 for MP state */
	    d_sdl = d-sdl;
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Jalpha[y][j][d_sdl] + tsc_v[yoffset]);
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[y][j][d_sdl] + tsc_v[yoffset]);
	  }
	  /* note we use 'd_sdr' not 'd_sd' for R, plus minimum d is sdr (1) */
	  for (d = sdr; d <= j; d++) { /* sdr == 1 for MP state */
	    d_sdr = d - sdr;
	    Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Jalpha[y][j_sdr][d_sdr] + tsc_v[yoffset]);
	    Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[y][j_sdr][d_sdr] + tsc_v[yoffset]);
	  }
	}
      }
      /* add in emission score */
      for (j = 0; j <= L; j++) {
	i = j;
	Jalpha[v][j][1] = IMPOSSIBLE;
	Lalpha[v][j][1] = lmesc_v[dsq[i]];
	Ralpha[v][j][1] = rmesc_v[dsq[j]];
	i--;
	for (d = 2; d <= j; d++) {
	  Jalpha[v][j][d] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	  Lalpha[v][j][d] += lmesc_v[dsq[i]];
	  Ralpha[v][j][d] += rmesc_v[dsq[j]];
	  i--;
	}
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = 0; j <= L; j++) {
	for (d = 1; d <= j; d++) {
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
	
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  
	  for (d = sd; d <= j; d++) {
	    d_sd = d-sd;
	    Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], Jalpha[y][j_sdr][d_sd] + tsc);
	    Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[y][j_sdr][d_sd] + tsc);
	    Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[y][j_sdr][d_sd] + tsc);
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

      for (j = 0; j <= L; j++) {
	for (d = 0; d <= j; d++) {
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
	  Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Jalpha[y][j][d]);
	  Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Lalpha[y][j][d]);
	  /* special case 2: k == d (full sequence aligns to BEGR_S right child */
	  Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Jalpha[z][j][d]);
	  Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Ralpha[z][j][d]);
	}
      }
    } /* end of B_st recursion */

    /* Now handle truncated begin transitions from ROOT_S, state 0 */
    /* Standard local begins are not allowed in truncated mode (they're kind of like J alignments though) */

    /* include full length hits in J matrix (much like a normal local begin) */
    if(cm->sttype[v] == B_st  || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || 
       cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
      Jalpha[0][L][L] = FLogsum(Jalpha[0][L][L], Jalpha[v][L][L] + trunc_penalty);
    }
    /* include full length truncated hits in L matrix */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) { 
      Lalpha[0][L][L] = FLogsum(Lalpha[0][L][L], Lalpha[v][L][L] + trunc_penalty);
    }	    
    /* include full length truncated hits in R matrix */
    if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
      Ralpha[0][L][L] = FLogsum(Ralpha[0][L][L], Ralpha[v][L][L] + trunc_penalty);
    }	    
    /* include full length truncated hits in T matrix */
    if(cm->sttype[v] == B_st) { 
      Talpha[0][L][L] = FLogsum(Talpha[0][L][L], Talpha[v][L][L] + trunc_penalty);
    }
  } /* end of for (v = cm->M-1; v > 0; v--) */
  FILE *fp1; fp1 = fopen("tmp.trimx", "w");   cm_tr_mx_Dump(fp1, mx); fclose(fp1);

  sc   = Jalpha[0][L][L];
  mode = TRMODE_J;
  if (Lalpha[0][L][L] > sc) { 
    sc   = Lalpha[0][L][L];
    mode = TRMODE_L;
  }
  if (Ralpha[0][L][L] > sc) { 
    sc   = Ralpha[0][L][L];
    mode = TRMODE_R;
  }
  if (Talpha[0][L][L] > sc) { 
    sc   = Talpha[0][L][L];
    mode = TRMODE_T;
  }

  if(ret_mode != NULL) *ret_mode = mode;    
  if(ret_sc   != NULL) *ret_sc   = sc;
  
  free(el_scA);

  ESL_DPRINTF1(("cm_TrInsideAlign() return sc: %f\n", sc));
  printf("cm_TrInsideAlign() return sc: %.4f\n", sc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/* Function: cm_TrInsideAlignHB()
 * Date:     EPN, Mon Sep 12 04:32:00 2011   
 *
 * Purpose: Run the truncated inside algorithm on a target sequence
 *           using bands in the j and d dimensions of the DP
 *           matrix. Bands were obtained from an HMM Forward-Backward
 *           parse of the target sequence. Uses float log odds scores.
 *           The full target sequence 1..L is aligned (only full
 *           alignments will contribute to the Inside score).
 *
 *           Very similar to cm_TrCYKInsideAlignHB(), see 'Purpose'
 *           of that function for more details. Only differences with
 *           that function is:
 *           - we do TrInside, not TrCYK
 *           - can't return a shadow matrix (we're not aligning)
 *           - doesn't return b, info about local begins 
 *
 *           This function complements cm_TrOutsideAlignHB().
 *
 * Args:     cm         - the model    [0..M-1]
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitized sequence
 *           L          - target sequence length
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           mx         - the dp matrix, only cells within bands in cp9b will be valid
 *           ret_mode   - RETURN: mode of optimal truncation mode, TRMODE_{J,L,R,T} if {J,L,R,T}alpha[0][L][L] is max scoring.
 *           ret_sc     - RETURN: log P(S|M)/P(S|R), as a bit score
 *                        NOTE: we don't sum over different marginal modes, we pick the highest scoring
 *                        one (J,L,R or T) and return {J,L,R,T}alpha[0][L][L] the sum of all complete 
 *                        J,L,R, or T alignments.
 * Returns:  <eslOK> on success.
 *
 * Throws:  <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *          <eslEINVAL> if the full sequence is not within the bands for state 0
 *          In either case alignment has been aborted, ret_sc is not valid
 */
int
cm_TrInsideAlignHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, CM_TR_HB_MX *mx, char *ret_mode, float *ret_sc)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;          /* temporary score */
  float    tsc;         /* a temporary variable holding a transition score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  float    bsc;		/* summed score for using all local begins */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */

  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      j_sdr;              /* j - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y, dp_z;   /* d index for state v/y/z in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dp_y_sdr;           /* dp_y - sdr */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Lp;                 /* L also changes depending on state */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  int      jp_0;               /* L offset in ROOT_S's (v==0) j band */
  int      Lp_0;               /* L offset in ROOT_S's (v==0) d band */

  /* variables related to truncated alignment (not in FastInsideAlignHB()) */
  float    trunc_penalty = 0.; /* penalty in bits for a truncated hit */
  int      mode = TRMODE_J;    /* truncation mode for obtaining optimal score <ret_sc> */
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
  /* ensure a full alignment to ROOT_S (v==0) is allowed by the bands */
  if (cp9b->jmin[0] > L || cp9b->jmax[0] < L)
    ESL_FAIL(eslEINVAL, errbuf, "cm_TrInsideAlignHB(): L (%d) is outside ROOT_S's j band (%d..%d)\n", L, cp9b->jmin[0], cp9b->jmax[0]);
  jp_0 = L - jmin[0];
  if (cp9b->hdmin[0][jp_0] > L || cp9b->hdmax[0][jp_0] < L) 
    ESL_FAIL(eslEINVAL, errbuf, "cm_TrInsideAlignHB(): L (%d) is outside ROOT_S's d band (%d..%d)\n", L, cp9b->hdmin[0][jp_0], cp9b->hdmax[0][jp_0]);
  Lp_0 = L - hdmin[0][jp_0];

  /* grow the matrix based on the current sequence and bands */
  if((status = cm_tr_hb_mx_GrowTo(cm, mx, errbuf, cp9b, L, size_limit)) != eslOK) return status;

  /* precalcuate all possible local end scores, for local end emits of 1..L residues */
  ESL_ALLOC(el_scA, sizeof(float) * (L+1));
  for(d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(mx->Jncells_valid > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(mx->Lncells_valid > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(mx->Rncells_valid > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(mx->Tncells_valid > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 

  /* Note, at this point in the non-banded version (cm_TrInsideAlign()) we
   * replace EL (cm->M) deck IMPOSSIBLEs with EL scores.  But we don't
   * here. It's actually not necessary there either, but it is done
   * for completeness and so if we check that matrix in combination
   * with an Outside matrix, then scores in the EL deck are valid.
   * But we won't do that here, and we're concerned with efficiency,
   * so we don't waste time with resetting the EL deck.
   */

  /* Main recursion */
  for (v = cm->M-1; v > 0; v--) { /* almost to ROOT, ROOT is special because all hits must use truncated begins */
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
	  if(hdmin[v][jp_v] >= sd) { 
	    d    = hdmin[v][jp_v];
	    dp_v = 0;
	  }
	  else { 
	    d    = sd;
	    dp_v = sd - hdmin[v][jp_v];
	  }
	  for (; d <= hdmax[v][jp_v]; dp_v++, d++) {
	    Jalpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
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
    else if(cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) { 
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IL states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1] 
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */
      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
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
	      for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
		yoffset = yvalidA[yvalid_idx];
		y = cm->cfirst[v] + yoffset;
		do_R_y = cp9b->do_R[y];
		do_J_y = cp9b->do_J[y];
		if((do_J_y || do_R_y) && (y != v)) { /* (y != v) part is to disallow IL self transits in R mode */
		  jp_y_sdr = j - jmin[y] - sdr;
		
		  /* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,L above) */
		  if((d) >= hdmin[y][jp_y_sdr] && (d) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		    dp_y = d - hdmin[y][jp_y_sdr];
		    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		    ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		    if(do_J_y) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]);
		    if(do_R_y) Ralpha[v][jp_v][dp_v] = FLogsum(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]);
		  }
		}
	      } /* end of for (yvalid_idx = 0... loop */
	      Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	  }
	}
      } /* end of if(! StateIsDetached(cm, v) */
    }
    else if(cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) { 
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IR states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1].
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */

      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
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
	/* Handle L separately */
	if(do_L_v) { 
	  /* The second MR_st/IR_st 'for (j...' loop is for the L matrix which use a different set of j values */
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    jp_v = j - jmin[v];
	    yvalid_ct = 0;
	  
	    /* determine which children y we can legally transit to for v, j */
	    /* we use 'j' and not 'j_sdr' here for the L matrix, differently from J and R matrices above */
	    for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	      if((y != v) && /* y == v when yoffset == 0 && v is an IR state: we don't want to allow IR self transits in L mode */
		 (j) >= jmin[y] && ((j) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j is valid for state y? */
	  
	    for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	      for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
		/* Note if we're an IL state, we can't self transit in R mode, this was ensured above when we set up yvalidA[] (xref:ELN3,p5)*/
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
		    if(do_J_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Jalpha[y][jp_y][dp_y] + tsc_v[yoffset]);
		    if(do_L_y) Lalpha[v][jp_v][dp_v] = FLogsum(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y][dp_y] + tsc_v[yoffset]);
		  }
		}
	      } /* end of for (yvalid_idx = 0... loop */
	      Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	  }
	}
      } /* end of if(! StateIsDetached(cm, v)) */
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
           
    /* Now handle truncated begin transitions from ROOT_S, state 0 */
    if(L >= jmin[v] && L <= jmax[v]) { 
      jp_v = L - jmin[v];
      Lp   = L - hdmin[v][jp_v];
      if(L >= hdmin[v][jp_v] && L <= hdmax[v][jp_v]) {
	/* If we get here alpha[v][jp_v][Lp] and alpha[0][jp_0][Lp0]
	 * are valid cells in the banded alpha matrix, corresponding to 
	 * alpha[v][L][L] and alpha[0][L][L] in the platonic matrix.
	 * (We've already made sure alpha[0][jp_0][Lp_0] was valid 
	 * at the beginning of the function.)
	 *
	 * Standard local begins are not allowed in truncated mode (they're kind of like J alignments though). 
	 */

	/* include full length hits in J matrix (much like a normal local begin) */
	if(do_J_v && (cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st
		      || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)) { 
	  Jalpha[0][jp_0][Lp_0] = FLogsum(Jalpha[0][jp_0][Lp_0], Jalpha[v][jp_v][Lp] + trunc_penalty);
	}
	/* include full length truncated hits in L matrix */
	if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) { 
	  Lalpha[0][jp_0][Lp_0] = FLogsum(Lalpha[0][jp_0][Lp_0], Lalpha[v][jp_0][Lp_0] + trunc_penalty);
	}	    
	/* include full length truncated hits in R matrix */
	if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
	  Ralpha[0][jp_0][Lp_0] = FLogsum(Ralpha[0][jp_0][Lp_0], Ralpha[v][jp_0][Lp_0] + trunc_penalty);
	}	    
	/* include full length truncated hits in T matrix */
	if(cm->sttype[v] == B_st) { 
	  Talpha[0][jp_0][Lp_0] = FLogsum(Talpha[0][jp_0][Lp_0], Talpha[v][jp_0][Lp_0] + trunc_penalty);
	}
      }
    }
  } /* end of for (v = cm->M-1; v > 0; v--) */
  /*FILE *fp1; fp1 = fopen("tmp.iamx", "w");   cm_tr_hb_mx_Dump(fp1, mx); fclose(fp1);*/
  
  sc   = Jalpha[0][jp_0][Lp_0];
  mode = TRMODE_J;
  if (Lalpha[0][jp_0][Lp_0] > sc) { 
    sc   = Lalpha[0][jp_0][Lp_0];
    mode = TRMODE_L;
  }
  if (Ralpha[0][jp_0][Lp_0] > sc) { 
    sc   = Ralpha[0][jp_0][Lp_0];
    mode = TRMODE_R;
  }
  if (Talpha[0][jp_0][Lp_0] > sc) { 
    sc   = Talpha[0][jp_0][Lp_0];
    mode = TRMODE_T;
  }

  if(ret_mode != NULL) *ret_mode = mode;    
  if(ret_sc   != NULL) *ret_sc   = sc;

  free(el_scA);
  free(yvalidA);

  if(ret_sc != NULL) *ret_sc = sc;

  ESL_DPRINTF1(("cm_TrInsideAlignHB() return sc: %f\n", sc));
  printf("cm_TrInsideAlignHB() return sc: %.4f\n", sc);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/* Function: cm_TrOptAccAlign()
 * based on cm_OptAccAlign()
 *
 * Date:     EPN, Wed Sep 28 13:16:12 2011
 *
 * Purpose: Run the truncated version of the Holmes/Durbin optimal
 *           accuracy algorithm on a full target sequence 1..L, given
 *           a pre-filled posterior matrix. Uses float log odds
 *           scores.  Non-banded version. See cm_OptAccAlignHB() for
 *           HMM banded version.
 * 
 *           A CM_TR_EMIT_MX matrix <emit_mx> must be passed in,
 *           filled by cm_TrEmitterPosterior(), with values:
 *
 *           {J,L}l_pp[v][i]: log of the posterior probability that
 *           state v emitted residue i leftwise either at (if a match
 *           state) or *after* (if an insert state) the left consensus
 *           position modeled by state v's node, either in Joint 
 *           marginal or Left marginal mode {J or L}.
 *
 *           {J,R}r_pp[v][i]: log of the posterior probability that
 *           state v emitted residue i rightwise either at (if a match
 *           state) or *before* (if an insert state) the right
 *           consensus position modeled by state v's node, either in
 *           Joint marginal or Right marginal mode {J or R}.
 *
 *           {J,L}l_pp[v] is NULL for states that do not emit leftwise
 *           {J,r}r_pp[v] is NULL for states that do not emit rightwise
 *
 *           Additionally, a CM_TR_MX DP matrix <mx> and
 *           CM_TR_SHADOW_MX <shmx> must be passed in. <shmx> will be
 *           expanded and filled here with traceback pointers to allow
 *           the optimally accurate parsetree to be recovered in
 *           cm_alignT() and <mx> will be expanded and filled with the
 *           optimal accuracy scores, where:
 *
 *           mx->{J,L,R,T}dp[v][j][d]: log of the sum of the posterior
 *           probabilities of emitting residues i..j in the subtree
 *           rooted at v given that v is in marginal mode J,L,R, or T.
 *
 *           The optimally accurate parsetree in marginal mode
 *           <opt_mode>, i.e. the parsetree that maximizes the sum of
 *           the posterior probabilities of all 1..L emitted residues,
 *           will be found.
 *
 * Args:     cm         - the model
 *           errbuf     - char buffer for reporting errors
 *           dsq        - the digitaized sequence [1..L]   
 *           L          - length of the dsq
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           opt_mode   - the optimal alignment mode, if unknown TRMODE_UNKNOWN is passed
 *           mx         - the DP matrix to fill in
 *           shmx       - the shadow matrix to fill in
 *           emit_mx    - pre-filled emit matrix
 *           ret_b      - optimal entry point (state) for the alignment
 *           ret_sc     - RETURN: average probability mass that goes through 
 *                        a cell of the optimally accurate parse
 *
 * Returns: <eslOK>     on success.
 * Throws:  <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *          If !eslOK: alignment has been aborted, ret_* variables are not valid
 */
int
cm_TrOptAccAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, CM_TR_MX *mx, CM_TR_SHADOW_MX *shmx, 
		 CM_TR_EMIT_MX *emit_mx, int *ret_b, float *ret_sc)
{
  int      status;          /* easel status code */
  int      v,y,z;	    /* indices for states  */
  int      j,d,i,k;	    /* indices in sequence dimensions */
  float    sc;		    /* temporary log odds score */
  int      yoffset;	    /* y=base+offset -- counter in child states that v can transit to */
  int      sd;              /* StateDelta(cm->sttype[v]) */
  int      sdl;             /* StateLeftDelta(cm->sttype[v] */
  int      sdr;             /* StateRightDelta(cm->sttype[v] */
  int      j_sdr;           /* j - sdr */
  int      d_sd;            /* d - sd */
  int      d_sdl;           /* d - sdl */
  int      d_sdr;           /* d - sdr */

  /* other variables used in truncated version, but not standard version (not in cm_CYKAlign()) */
  int   b;		    /* best truncated entry state */
  float bsc;		    /* score for using the best truncated entry state */
  int   have_el;            /* TRUE if local ends are on */
  int   Lyoffset0;          /* first yoffset to use for updating L matrix in IR/MR states, 1 if IR, 0 if MR */
  int   Ryoffset0;          /* first yoffset to use for updating R matrix in IL/ML states, 1 if IL, 0 if ML */
  int   fill_L, fill_R, fill_T; /* must we fill in the L, R, and T matrices? */

  /* the DP matrices */
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

  float  **Jl_pp    = emit_mx->Jl_pp; /* pointer to the prefilled posterior values for left  emitters in Joint mode */
  float  **Ll_pp    = emit_mx->Ll_pp; /* pointer to the prefilled posterior values for left  emitters in Left  mode */
  float  **Jr_pp    = emit_mx->Jr_pp; /* pointer to the prefilled posterior values for right emitters in Joint mode */
  float  **Rr_pp    = emit_mx->Rr_pp; /* pointer to the prefilled posterior values for right emitters in Right mode */

  /* Determine which matrices we need to fill in, based on <opt_mode> */
  if((status = cm_TrFillFromMode(opt_mode, &fill_L, &fill_R, &fill_T)) != eslOK) ESL_FAIL(status, errbuf, "cm_TrOptAccAlign(), bogus mode: %d", opt_mode);
  printf("in cm_TrOptAcc(), mode: %d fill_L: %d fill_R: %d fill_T: %d\n", 
	 opt_mode, fill_L, fill_R, fill_T);

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo       (cm, mx,   errbuf, L, size_limit)) != eslOK) return status;
  if((status = cm_tr_shadow_mx_GrowTo(cm, shmx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 
  if(shmx->Jy_ncells_valid > 0) for(i = 0; i < shmx->Jy_ncells_valid; i++) shmx->Jyshadow_mem[i] = USED_EL;
  if(shmx->Ly_ncells_valid > 0) for(i = 0; i < shmx->Ly_ncells_valid; i++) shmx->Lyshadow_mem[i] = USED_TRUNC_END;
  if(shmx->Ry_ncells_valid > 0) for(i = 0; i < shmx->Ry_ncells_valid; i++) shmx->Ryshadow_mem[i] = USED_TRUNC_END;
  /* for B states, shadow matrix holds k, length of right fragment, this will almost certainly be overwritten */
  if(shmx->Jk_ncells_valid > 0) esl_vec_ISet(shmx->Jkshadow_mem, shmx->Jk_ncells_valid, 0); 
  if(shmx->Lk_ncells_valid > 0) esl_vec_ISet(shmx->Lkshadow_mem, shmx->Lk_ncells_valid, 0);
  if(shmx->Rk_ncells_valid > 0) esl_vec_ISet(shmx->Rkshadow_mem, shmx->Rk_ncells_valid, 0);
  if(shmx->Tk_ncells_valid > 0) esl_vec_ISet(shmx->Tkshadow_mem, shmx->Tk_ncells_valid, 0);
  if(shmx->Lk_ncells_valid > 0) for(i = 0; i < shmx->Lk_ncells_valid; i++) shmx->Lkmode_mem[i] = TRMODE_J;
  if(shmx->Rk_ncells_valid > 0) for(i = 0; i < shmx->Rk_ncells_valid; i++) shmx->Rkmode_mem[i] = TRMODE_J;

  /* if local ends are on, replace the EL deck IMPOSSIBLEs with EL scores */
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  if(have_el) { 
    for (j = 0; j <= L; j++) {
      Jalpha[cm->M][j][0] = IMPOSSIBLE;
      for (d = 1;  d <= j; d++) { 
	Jalpha[cm->M][j][d] = FLogsum(Jalpha[cm->M][j][d-1], Jl_pp[cm->M][j]); /* optimal (and only) parse for EL is to emit all d residues */
      }
    }
  }

  /* Main recursion */
  for (v = cm->M-1; v > 0; v--) { /* almost to ROOT, ROOT is special because all hits must use truncated begins */
    sd   = StateDelta(cm->sttype[v]);
    sdl  = StateLeftDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);

    /* re-initialize if we can do a local end from v and check for a
     * special optimal-accuracy-specific initialization case
     */
    if((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = 0; j <= L; j++) {
	/* copy values from saved EL deck */
	for (d = sd; d <= j; d++) { 
	  Jalpha[v][j][d] = Jalpha[cm->M][j-sdr][d-sd];
	  /* Jyshadow[v][j][d] remains USED_EL. 
	   * L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL 
	   */
	}
      }
    }
    else if(cm->sttype[v] != B_st && cm->sttype[v] != E_st) { /* && cm->endsc[v] == IMPOSSIBLE */
      for (j = 0; j <= L; j++) {
	/* Check for special initialization case, specific to
	 * optimal_accuracy alignment, normally (with TrCYK for
	 * example) we init shadow matrix to USED_EL for all cells
	 * b/c we know that will be overwritten for the most
	 * likely transition, but with optimal accuracy, only
	 * emissions add to the score, so when d == sd, we know
	 * we'll emit sd residues from v, so the initialization
	 * will NOT be overwritten. We get around this for
	 * cells for which  d == sd and v is a state that has 
	 * a StateDelta=0 child y (DELETE or END) by initializing
	 * that transition to y is most likely.
	 */
	for (d = 0; d <= sd; d++) { 
	  y = cm->cfirst[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
	    if(StateDelta(cm->sttype[y+yoffset]) == 0) { 
	      Jyshadow[v][j][d] = yoffset;
	    }
	  }
	}
	/* for the L matrix, the same logic applies as above, but now we 
	 * have to be careful to include any state that doesn't emit left
	 */
	if(fill_L) { 
	  for (d = 0; d <= sdl; d++) { 
	    y = cm->cfirst[v];
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
	      if(StateLeftDelta(cm->sttype[y+yoffset]) == 0) { 
		Lyshadow[v][j][d] = yoffset;
	      }
	    }
	  }
	}
	/* for the R matrix, the same logic applies as above, but now we 
	 * have to be careful to include any state that doesn't emit right
	 */
	if(fill_R) { 
	  for (d = 0; d <= sdr; d++) { 
	    y = cm->cfirst[v];
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
	      if(StateRightDelta(cm->sttype[y+yoffset]) == 0) { 
		Ryshadow[v][j][d] = yoffset;
	      }
	    }
	  }
	}
      }
    }

    /* note there's no E state update here, those cells all remain IMPOSSIBLE */

    if(cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) { 
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

      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
	Ryoffset0 = cm->sttype[v] == IL_st ? 1 : 0; /* don't allow IL self transits in R mode */
	for (j = sdr; j <= L; j++) {
	  i = j;
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++, i--) {
	    d_sd = d - sd;
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	      y = cm->cfirst[v] + yoffset; 
	      if ((sc = Jalpha[y][j_sdr][d_sd]) > Jalpha[v][j][d]) {
		Jalpha[v][j][d]   = sc; 
		Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	      }
	      if (fill_L && ((sc = Lalpha[y][j_sdr][d_sd]) > Lalpha[v][j][d])) {
		Lalpha[v][j][d]   = sc; 
		Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	      }
	    }
	    Jalpha[v][j][d]  = FLogsum(Jalpha[v][j][d], Jl_pp[v][i]);
	    Jalpha[v][j][d]  = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);

	    if(fill_L) { 
	      if(d >= 2) { 
		Lalpha[v][j][d]  = FLogsum(Lalpha[v][j][d], Ll_pp[v][i]);
	      }
	      else { 
		Lalpha[v][j][d]   = Ll_pp[v][i]; /* actually I think this will give the same value as d >= 2 case above */
		Lyshadow[v][j][d] = USED_TRUNC_END;
	      }
	      Lalpha[v][j][d]  = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	    }

	    /* handle R separately */
	    if(fill_R) { 
	      /* note we use 'd', not 'd_sd' (which we used in the corresponding loop for J,L above) */
	      for (yoffset = Ryoffset0; yoffset < cm->cnum[v]; yoffset++) { /* using Ryoffset0 instead of 0 disallows IL self transits in R mode */
		y = cm->cfirst[v] + yoffset; 
		if ((sc = Jalpha[y][j_sdr][d]) > Ralpha[v][j][d]) { 
		  Ralpha[v][j][d] = sc; 
		  Ryshadow[v][j][d]= yoffset + TRMODE_J_OFFSET;
		}
		if ((sc = Ralpha[y][j_sdr][d]) > Ralpha[v][j][d]) { 
		  Ralpha[v][j][d] = sc;
		  Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
		}
	      }
	      /* no residue was emitted if we're in R mode */
	      Ralpha[v][j][d] = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
	    }
	  }
	}
      } /* end of if(! StateIsDetached(cm,v )) */
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

      if(! StateIsDetached(cm, v)) { /* if we're detached (unreachable), leave all {J,L,R}alpha values as they were initialized, as IMPOSSIBLE */
	Lyoffset0 = cm->sttype[v] == IR_st ? 1 : 0; /* don't allow IR self transits in L mode */
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++) {
	    d_sd = d - sd;
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	      y = cm->cfirst[v] + yoffset; 
	      if ((sc = Jalpha[y][j_sdr][d_sd]) > Jalpha[v][j][d]) {
		Jalpha[v][j][d]   = sc; 
		Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	      }
	      if (fill_R && ((sc = Ralpha[y][j_sdr][d_sd]) > Ralpha[v][j][d])) {
		Ralpha[v][j][d]   = sc; 
		Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	      }
	    }
	    Jalpha[v][j][d]  = FLogsum(Jalpha[v][j][d], Jr_pp[v][j]);
	    Jalpha[v][j][d]  = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);

	    if(fill_R) { 
	      if(d >= 2) { 
		Ralpha[v][j][d]  = FLogsum(Ralpha[v][j][d], Rr_pp[v][j]);
	      }
	      else { 
		Ralpha[v][j][d]   = Rr_pp[v][j]; /* actually I think this will give the same value as d >= 2 case above */
		Ryshadow[v][j][d] = USED_TRUNC_END;
	      }
	      Ralpha[v][j][d]  = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
	    }

	    /* handle L separately */
	    if(fill_L) { 
	      /* note we use 'j' and 'd', not 'j_sdr' and 'd_sd' (which we used in the corresponding loop for J,R above) */
	      for (yoffset = Lyoffset0; yoffset < cm->cnum[v]; yoffset++) { /* using Lyoffset0, instead of 0 disallows IR self transits in L mode */
		y = cm->cfirst[v] + yoffset; 
		if ((sc = Jalpha[y][j][d]) > Lalpha[v][j][d]) { 
		  Lalpha[v][j][d] = sc;
		  Lyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
		}
		if ((sc = Lalpha[y][j][d]) > Lalpha[v][j][d]) { 
		  Lalpha[v][j][d] = sc;
		  Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
		}
	      }
	      /* no residue was emitted if we're in R mode */
	      Lalpha[v][j][d] = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	    }
	  }
	}
      } /* end of if(! StateIsDetached(cm, v)) */
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

	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  for (d = sd; d <= j; d++) { /* sd == 2 for MP state */
	    d_sd = d-sd;
	    if((sc = Jalpha[y][j_sdr][d_sd]) > Jalpha[v][j][d]) {
	      Jalpha[v][j][d]   = sc;
	      Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	  }
	  if(fill_L) { 
	    /* note we use 'j' and 'd_sdl' not 'j_sdr' for 'd_sd' for L, plus minimum d is sdl (1) */
	    for (d = sdl; d <= j; d++) { /* sdl == 1 for MP state */
	      d_sdl = d-sdl;
	      if((sc = Jalpha[y][j][d_sdl]) > Lalpha[v][j][d]) {
		Lalpha[v][j][d]   = sc;
		Lyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	      }
	      if((sc = Lalpha[y][j][d_sdl]) > Lalpha[v][j][d]) {
		Lalpha[v][j][d]   = sc;
		Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	      }
	    }
	  }
	  if(fill_R) { 
	    /* note we use 'd_sdr' not 'd_sd' for R, plus minimum d is sdr (1) */
	    for (d = sdr; d <= j; d++) { /* sdr == 1 for MP state */
	      d_sdr = d - sdr;
	      if((sc = Jalpha[y][j_sdr][d_sdr]) > Ralpha[v][j][d]) {
		Ralpha[v][j][d]   = sc;
		Ryshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	      }
	      if((sc = Ralpha[y][j_sdr][d_sdr]) > Ralpha[v][j][d]) {
		Ralpha[v][j][d]   = sc;
		Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	      }
	    }
	  }
	}
      }
      /* add in emission score */
      for (j = 0; j <= L; j++) {
	i = j;
	Jalpha[v][j][1] = IMPOSSIBLE;
	if(fill_L) { 
	  Lalpha[v][j][1]   = Ll_pp[v][i];
	  Lyshadow[v][j][1] = USED_TRUNC_END;
	}
	if(fill_R) { 
	  Ralpha[v][j][1]   = Rr_pp[v][j];
	  Ryshadow[v][j][1] = USED_TRUNC_END;
	}
	for (d = 2; d <= j; d++, i--) { 
	  Jalpha[v][j][d] = FLogsum(Jalpha[v][j][d], FLogsum(Jl_pp[v][i], Jr_pp[v][j]));
	}
	if(fill_L) for (d = 2; d <= j; d++) Lalpha[v][j][d] = FLogsum(Lalpha[v][j][d], Ll_pp[v][i]);
	if(fill_R) for (d = 2; d <= j; d++) Ralpha[v][j][d] = FLogsum(Ralpha[v][j][d], Rr_pp[v][j]);
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = 0; j <= L; j++) {
	for (d = 1; d <= j; d++) Jalpha[v][j][d] = ESL_MAX(Jalpha[v][j][d], IMPOSSIBLE);
	if(fill_L) for (d = 1; d <= j; d++) Lalpha[v][j][d] = ESL_MAX(Lalpha[v][j][d], IMPOSSIBLE);
	if(fill_R) for (d = 1; d <= j; d++) Ralpha[v][j][d] = ESL_MAX(Ralpha[v][j][d], IMPOSSIBLE);
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
	
	for (j = sdr; j <= L; j++) {
	  j_sdr = j - sdr;
	  
	  for (d = sd; d <= j; d++) {
	    d_sd = d-sd;
	    if((sc = Jalpha[y][j_sdr][d_sd]) > Jalpha[v][j][d]) {
	      Jalpha[v][j][d]   = sc;
	      Jyshadow[v][j][d] = yoffset + TRMODE_J_OFFSET;
	    }
	  }
	  if(fill_L) {
	    for (d = sd; d <= j; d++) { 
	      d_sd = d-sd;
	      if((sc = Lalpha[y][j_sdr][d_sd]) > Lalpha[v][j][d]) {
		Lalpha[v][j][d]   = sc;
		Lyshadow[v][j][d] = yoffset + TRMODE_L_OFFSET;
	      }
	    }
	  }
	  if(fill_R) {
	    for (d = sd; d <= j; d++) { 
	      d_sd = d-sd;
	      if((sc = Ralpha[y][j_sdr][d_sd]) > Ralpha[v][j][d]) {
		Ralpha[v][j][d]   = sc;
		Ryshadow[v][j][d] = yoffset + TRMODE_R_OFFSET;
	      }
	    }
	  }
	  /* an easy to overlook case: if d == 0, ensure L and R values are IMPOSSIBLE */
	  if(fill_L) Lalpha[v][j][0] = IMPOSSIBLE;
	  if(fill_R) Ralpha[v][j][0] = IMPOSSIBLE;
	}
      }
      /* no emission score to add */
    }
    else { /* B_st */
      assert(cm->sttype[v] == B_st);
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */

      for (j = 0; j <= L; j++) {
	for (d = 0; d <= j; d++) {
	  for (k = 1; k < d; k++) {
	    if((NOT_IMPOSSIBLE(Jalpha[y][j-k][d-k])) && 
	       (NOT_IMPOSSIBLE(Jalpha[z][j][k])) && 
	       ((sc = Jalpha[y][j-k][d-k] + Jalpha[z][j][k]) > Jalpha[v][j][d])) { 
	      Jalpha[v][j][d]   = sc;
	      Jkshadow[v][j][d] = k;
	    }
	  }
	  if(fill_L) { 
	    for (k = 1; k < d; k++) {
	      if((NOT_IMPOSSIBLE(Jalpha[y][j-k][d-k])) && 
		 (NOT_IMPOSSIBLE(Lalpha[z][j][k])) && 
		 ((sc = Jalpha[y][j-k][d-k] + Lalpha[z][j][k]) > Lalpha[v][j][d])) { 
		Lalpha[v][j][d]   = sc;
		Lkshadow[v][j][d] = k;
		Lkmode[v][j][d]   = TRMODE_J;
	      }
	    }
	  }
	  if(fill_R) { 
	    for (k = 1; k < d; k++) {
	      if((NOT_IMPOSSIBLE(Ralpha[y][j-k][d-k])) && 
		 (NOT_IMPOSSIBLE(Jalpha[z][j][k])) && 
		 ((sc = Ralpha[y][j-k][d-k] + Jalpha[z][j][k]) > Ralpha[v][j][d])) { 
		Ralpha[v][j][d]   = sc;
		Rkshadow[v][j][d] = k;
		Rkmode[v][j][d]   = TRMODE_J;
	      }
	    }
	  }
	  if(fill_T) { 
	    for (k = 1; k < d; k++) {
	      /*if((k != i-1) && (k != j)) {*/
	      if((NOT_IMPOSSIBLE(Ralpha[y][j-k][d-k])) && 
		 (NOT_IMPOSSIBLE(Lalpha[z][j][k])) && 
		 ((sc = Ralpha[y][j-k][d-k] + Lalpha[z][j][k]) > Talpha[v][j][d])) { 
		Talpha[v][j][d]   = sc;
		Tkshadow[v][j][d] = k;
	      }
	      /*}*/
	    }
	  }
	  /* two additional special cases in trCYK (these are not in standard CYK) */
	  /* special case 1: k == 0 (full sequence aligns to BEGL_S left child */
	  if(fill_L) { 
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
	  }
	  /* special case 2: k == d (full sequence aligns to BEGR_S right child */
	  if(fill_R) { 
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
      }
    }/* end of B_st recursion */

    /* Now handle truncated begin transitions from ROOT_S, state 0 */
    /* We know the optimal mode, it was passed in as <opt_mode> */

    if(opt_mode == TRMODE_J) { 
      /* check if we have a new optimally scoring Joint alignment in J matrix (much like a standard local begin) */
      if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st
	 || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
	if(bsc < Jalpha[v][L][L]) { 
	  bsc = Jalpha[v][L][L];
	  b   = v;
	}
      }
    }
    else if(opt_mode == TRMODE_L) { 
      /* check if we have a new optimally scoring Left alignment in L matrix */
      if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) { 
	if(bsc < Lalpha[v][L][L]) { 
	  bsc = Lalpha[v][L][L];
	  b   = v;
	}
      }
    }
    else if(opt_mode == TRMODE_R) { 
      /* check if we have a new optimally scoring Right alignment in R matrix */
      if(cm->sttype[v] == B_st || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
	if(bsc < Ralpha[v][L][L]) { 
	  bsc = Ralpha[v][L][L];
	  b  = v;
	}
      }	    
    }
    else if(opt_mode == TRMODE_T) { 
      /* check if we have a new optimally scoring Terminal alignment in T matrix */
      if(cm->sttype[v] == B_st) { 
	if(bsc < Talpha[v][L][L]) { 
	  bsc = Talpha[v][L][L];
	  b  = v;
	}
      }	    
    }
  } /* end loop for (v = cm->M-1; v > 0; v--) */
  FILE *fp1; fp1 = fopen("tmp.troamx", "w");   cm_tr_mx_Dump(fp1, mx); fclose(fp1);
  FILE *fp2; fp2 = fopen("tmp.troashmx", "w"); cm_tr_shadow_mx_Dump(fp2, cm, shmx); fclose(fp2);

  if(opt_mode == TRMODE_J) Jalpha[0][L][L] = bsc;
  if(opt_mode == TRMODE_L) Lalpha[0][L][L] = bsc;
  if(opt_mode == TRMODE_R) Ralpha[0][L][L] = bsc;
  if(opt_mode == TRMODE_T) Talpha[0][L][L] = bsc;

  /* convert bsc, a log probability, into the average posterior probability of all L aligned residues */
  sc = bsc;
  sc = sreEXP2(sc) / (float) L;

  if(ret_b  != NULL) *ret_b  = b;    
  if(ret_sc != NULL) *ret_sc = sc;

  ESL_DPRINTF1(("cm_TrOptAccAlign return pp: %f\n", sc));
  printf("cm_TrOptAccAlign return pp: %f\n", sc);
  return eslOK;
}

/* Function: cm_TrCYKOutsideAlign()
 * Date:     EPN, Wed Sep 14 14:20:20 2011
 *
 * Purpose:  Run the outside TrCYK algorithm on a target sequence.
 *           Non-banded version. See cm_TrCYKOutsideAlignHB() for
 *           the HMM banded version. The full target sequence
 *           1..L is aligned. 
 * 
 *           Very similar to cm_TrOutsideAlign() but calculates
 *           beta[v][j][d]: log probability of the most likely parse
 *           that emits 1..i-1 and j+1..L and passes through v at j,d
 *           (where i = j-d+1) instead of the log of the summed
 *           probability of all such parses. This means max operations
 *           are used instead of logsums.
 *
 *           Meaning of cells:
 *
 *           Jbeta[v][j][d]: log prob of the most likely parse that
 *                           emits 1..i-1 and j+1..L and passes through 
 *                           v in Joint marginal mode at j,d.
 *           Lbeta[v][j][d]: log prob of the most likely parse that
 *                           emits 1..i-1 and j+1..L and passes through
 *                           v in Left marginal mode at j,d.
 *           Rbeta[v][j][d]: log prob of the most likely parse that
 *                           emits 1..i-1 and j+1..L and passes through 
 *                           v in Right marginal mode at j,d.
 *
 *           This function complements cm_TrCYKInsideAlign() but is
 *           mainly useful for testing and reference. It can be used
 *           with do_check=TRUE to verify that the implementation of
 *           CYKTrInside and CYKTrOutside are consistent.  Because the
 *           structure of CYKTrInside and TrInside, and CYKTrOutside
 *           and TrOutside are so similar and the CYK variants are
 *           easier to debug (because only the optimal parsetree is
 *           considered instead of all possible parsetrees) this
 *           function can be useful for finding bugs in Outside.  It
 *           is currently not hooked up to any of the main Infernal
 *           programs.
 *
 * Args:     cm        - the model
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           L         - length of the dsq to align
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_check  - TRUE to attempt to check 
 *           opt_mode  - TRMODE_J, TRMODE_L, TRMODE_R, or TRMODE_T, the optimal
 *                       alignment mode, we'll only allow alignments in this mode.
 *           mx        - the dp matrix, grown and filled here
 *           inscyk_mx - the pre-filled dp matrix from the CYK Inside calculation 
 *                       (performed by cm_CYKInsideAlign(), required)
 *
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *           <eslEMEM>   if we run out of memory
 *           <eslFAIL>   if <do_check>==TRUE and we fail a test
 */
int
cm_TrCYKOutsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, 
		     char opt_mode, CM_TR_MX *mx, CM_TR_MX *inscyk_mx)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float    Jsc,Lsc,Rsc,Tsc;    /* a temporary variable holding a float score */
  float    optsc;              /* the optimal score from the Inside matrix */
  float    escore;	       /* an emission score, tmp variable */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  /* variables used only if do_check */
  int      fail1_flag = FALSE; /* set to TRUE if do_check and we see a problem */
  int      fail2_flag = FALSE; /* set to TRUE if do_check and we see a problem */
  int      vmax;               /* i, offset in the matrix */
  int     *optseen = NULL;     /* [1..i..L] TRUE is residue i is accounted for in optimal parse */
  /* indices used in the depths of the DP recursion */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdl;                /* StateLeftDelta(cm->sttype[y] */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */

  /* DP matrix variables */
  float ***Jbeta   = mx->Jdp;     /* pointer to the outside Jbeta DP matrix */
  float ***Lbeta   = mx->Ldp;     /* pointer to the outside Lbeta DP matrix */
  float ***Rbeta   = mx->Rdp;     /* pointer to the outside Rbeta DP matrix */
  float ***Tbeta   = mx->Tdp;     /* pointer to the outside Tbeta DP matrix */

  float ***Jalpha  = inscyk_mx->Jdp; /* pointer to the precalc'ed inside Jalpha DP matrix */
  float ***Lalpha  = inscyk_mx->Ldp; /* pointer to the precalc'ed inside Lalpha DP matrix */
  float ***Ralpha  = inscyk_mx->Rdp; /* pointer to the precalc'ed inside Ralpha DP matrix */
  /* Talpha is not used */

  /* Allocations and initializations */

  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo(cm, mx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 

  /* set cells in the special ROOT_S deck corresponding to full sequence alignments to 0. */
  if     (opt_mode == TRMODE_J) Jbeta[0][L][L] = 0.; /* a full Joint    alignment is outside this cell */
  else if(opt_mode == TRMODE_L) Lbeta[0][L][L] = 0.; /* a full Left     alignment is outside this cell */
  else if(opt_mode == TRMODE_R) Rbeta[0][L][L] = 0.; /* a full Right    alignment is outside this cell */
  else if(opt_mode == TRMODE_T) Tbeta[0][L][L] = 0.; /* a full Terminal alignment is outside this cell */

  /* we also need to allow full truncated alignments, rooted 
   * at internal truncated entry/exit points.
   */ 
  for(v = 0; v < cm->M; v++) { 
    if(! StateIsDetached(cm, v)) { 
      if(opt_mode == TRMODE_J && 
	 (cm->sttype[v] == B_st  || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || 
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)) { 
	Jbeta[v][L][L] = 0.; /* a full Joint alignment is outside this cell */
      }
      if(opt_mode == TRMODE_L && 
	 (cm->sttype[v] == B_st  || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st)) { 
	if(opt_mode == TRMODE_L) Lbeta[v][L][L] = 0.; /* a full Left  alignment is outside this cell */
      }
      if(opt_mode == TRMODE_R && 
	 (cm->sttype[v] == B_st  || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st)) { 
	if(opt_mode == TRMODE_R) Rbeta[v][L][L] = 0.; /* a full Left  alignment is outside this cell */
      }
      if(opt_mode == TRMODE_T && 
	 cm->sttype[v] == B_st) { 
	if(opt_mode == TRMODE_T) Tbeta[v][L][L] = 0.; /* a full Terminal alignment is outside this cell */
      }
    }
  }

  /* main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    if(! StateIsDetached(cm, v)) { 
      sd  = StateDelta(cm->sttype[v]);
      sdr = StateRightDelta(cm->sttype[v]);

      if (cm->stid[v] == BEGL_S) { /* BEGL_S */
	y = cm->plast[v];	/* the parent bifurcation    */
	z = cm->cnum[y];	/* the other (right) S state */
	for(j = 0; j <= L; j++) { 
	  for (d = 0; d <= j; d++) {
	    for (k = 0; k <= (L-j); k++) {
	      Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Jbeta[y][j+k][d+k] + Jalpha[z][j+k][k]); /* A */
	      Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Lbeta[y][j+k][d+k] + Lalpha[z][j+k][k]); /* B */
	      Rbeta[v][j][d] = ESL_MAX(Rbeta[v][j][d], Rbeta[y][j+k][d+k] + Jalpha[z][j+k][k]); /* C */
	      if(d == j && (j+k) == L) { 
		Rbeta[v][j][d] = ESL_MAX(Rbeta[v][j][d], Tbeta[y][j+k][d+k] + Lalpha[z][j+k][k]); /* D */
		/* Note: Tbeta[y][j+k==L][d+k==L] will be 0.0 because
		 * it was initialized that way. That T cell includes the
		 * full target 1..L (any valid T alignment must
		 * because we must account for the full target) rooted
		 * at a B state, and a transition from that B state to
		 * this BEGL_S is always probability 1.0.
		 */
	      }
	    }
	    Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Lbeta[y][j][d]); /* entire sequence on left, no sequence on right, k == 0 */
	    Lbeta[v][j][d] = ESL_MAX(Lbeta[v][j][d], Lbeta[y][j][d]); /* entire sequence on left, no sequence on right, k == 0 */
	  }
	}
      } /* end of 'if (cm->stid[v] == BEGL_S */
      else if (cm->stid[v] == BEGR_S) {
	y = cm->plast[v];	  /* the parent bifurcation    */
	z = cm->cfirst[y];  /* the other (left) S state  */
	for(j = 0; j <= L; j++) { 
	  for (d = 0; d <= j; d++) {
	    i = j-d+1;
	    for (k = 0; k <= (j-d); k++) {
	      Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Jbeta[y][j][d+k] + Jalpha[z][j-d][k]); /* A */
	      Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Rbeta[y][j][d+k] + Ralpha[z][j-d][k]); /* C */
	      Lbeta[v][j][d] = ESL_MAX(Lbeta[v][j][d], Lbeta[y][j][d+k] + Jalpha[z][j-d][k]); /* B */
	      if(k == (i-1) && j == L) { 
		Lbeta[v][j][d] = ESL_MAX(Lbeta[v][j][d], Tbeta[y][j][d+k] + Ralpha[z][j-d][k]); /* D */
		/* Note: Tbeta[y][j==L][d+k==L] will be 0.0 because it
		 * was initialized that way. That T cell includes the
		 * full target 1..L (any valid T alignment must
		 * because we must account for the full target) rooted
		 * at a B state, and a transition from that B state to
		 * this BEGR_S is always probability 1.0.
		 */
	      }
	    }
	    Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Rbeta[y][j][d]); /* entire sequence on right, no sequence on left, k == 0 */
	    Rbeta[v][j][d] = ESL_MAX(Rbeta[v][j][d], Rbeta[y][j][d]); /* entire sequence on right, no sequence on left, k == 0 */
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
	      sdl = StateLeftDelta(cm->sttype[y]);
	      sdr = StateRightDelta(cm->sttype[y]);
	      switch(cm->sttype[y]) {
	      case MP_st: 
		if(j != L && d != j) { 
		  escore = cm->oesc[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Jbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		}
		if(j == L && d != j) {
		  escore = cm->lmesc[y][dsq[i-1]];
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Lbeta[y][j][d+sdl]     + cm->tsc[y][voffset] + escore);
		  Lbeta[v][j][d] = ESL_MAX(Lbeta[v][j][d], Lbeta[y][j][d+sdl]     + cm->tsc[y][voffset] + escore);
		}
		if(i == 1 && j != L) { 
		  escore = cm->rmesc[y][dsq[j+1]];
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Rbeta[y][j+sdr][d+sdr] + cm->tsc[y][voffset] + escore);
		  Rbeta[v][j][d] = ESL_MAX(Rbeta[v][j][d], Rbeta[y][j+sdr][d+sdr] + cm->tsc[y][voffset] + escore);
		}
		break;
	      case ML_st:
	      case IL_st: 
		if (d != j) { 
		  escore = cm->oesc[y][dsq[i-1]];
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Jbeta[y][j][d+sd]     + cm->tsc[y][voffset] + escore);
		  Lbeta[v][j][d] = ESL_MAX(Lbeta[v][j][d], Lbeta[y][j][d+sd]     + cm->tsc[y][voffset] + escore);
		}
		if(i == 1 && /* only allow transition from R if we're emitting first residue 1 from y  */
		   v != y) {  /* will only happen if v == IL, we don't allow silent self transitions from IL->IL */
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Rbeta[y][j][d]        + cm->tsc[y][voffset]);
		  Rbeta[v][j][d] = ESL_MAX(Rbeta[v][j][d], Rbeta[y][j][d]        + cm->tsc[y][voffset]);
		}
		break;
	      case MR_st:
	      case IR_st:
		if (j != L) { 
		  escore = cm->oesc[y][dsq[j+1]];
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Jbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		  Rbeta[v][j][d] = ESL_MAX(Rbeta[v][j][d], Rbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		}
		if(j == L && /* only allow transition from L if we're emitting final residue L from y */ 
		   v != y) {  /* will only happen if v == IR, we don't allow silent self transitions from IR->IR */
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Lbeta[y][j][d]        + cm->tsc[y][voffset]);
		  Lbeta[v][j][d] = ESL_MAX(Lbeta[v][j][d], Lbeta[y][j][d]        + cm->tsc[y][voffset]);
		}
		break;
	      case S_st:
	      case E_st:
	      case D_st:
		if(y != 0) { 
		  Jbeta[v][j][d] = ESL_MAX(Jbeta[v][j][d], Jbeta[y][j][d] + cm->tsc[y][voffset]);
		  Lbeta[v][j][d] = ESL_MAX(Lbeta[v][j][d], Lbeta[y][j][d] + cm->tsc[y][voffset]);
		  Rbeta[v][j][d] = ESL_MAX(Rbeta[v][j][d], Rbeta[y][j][d] + cm->tsc[y][voffset]);
		  /* y can't be 0 because transitions from ROOT_S are
		   * free in truncated mode, and we've already handled
		   * them by initializing all valid full alignments to 0.0.
		   */
		}
		break;
	      } /* end of switch(cm->sttype[y] */  
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (Jbeta[v][j][d] < IMPOSSIBLE) Jbeta[v][j][d] = IMPOSSIBLE;
	  } /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
	} /* end loop over j. We know beta for this whole state */
      } /* end of 'else' (if cm->sttype[v] != BEGL_S, BEGR_S) */
    } /* end of 'if(! StateIsDetached(cm, v))' */
    /* we're done calculating deck v for everything but local ends */

    /* deal with local alignment end transitions v->EL J matrix only (EL = deck at M.) */
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);   /* note emitmode is for state v */
      
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  i = j-d+1;
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    if (j == L || d == j) continue; /* boundary condition */
	    escore = cm->oesc[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    Jbeta[cm->M][j][d] = ESL_MAX(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == j) continue;	
	    escore = cm->oesc[v][dsq[i-1]];
	    Jbeta[cm->M][j][d] = ESL_MAX(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == L) continue;
	    escore = cm->oesc[v][dsq[j+1]];
	    Jbeta[cm->M][j][d] = ESL_MAX(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case S_st:
	  case D_st:
	  case E_st:
	    Jbeta[cm->M][j][d] = ESL_MAX(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v]));
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
	Jbeta[cm->M][j][d] = ESL_MAX(Jbeta[cm->M][j][d], Jbeta[cm->M][j][d+1] + cm->el_selfsc);
    }
  }

  fail1_flag = FALSE;
  fail2_flag = FALSE;
  printf("DO CHECK: %d\n", do_check);
  if(do_check) {
    /* Check for consistency between the Inside alpha matrix and the
     * Outside beta matrix. we assume the Inside CYK parse score
     * (optsc) is the optimal score, so for all v,j,d:
     * 
     * Jalpha[v][j][d] + Jbeta[v][j][d] <= optsc
     * Lalpha[v][j][d] + Lbeta[v][j][d] <= optsc
     * Ralpha[v][j][d] + Rbeta[v][j][d] <= optsc
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
     * and/or the other) but that should be extremely unlikely.  If we
     * do this test many times for many different models and pass, we
     * should be confident we have consistent implementations.
     * 
     * This is an expensive check and should only be done while
     * debugging.
     */
    ESL_ALLOC(optseen, sizeof(int) * (L+1));
    esl_vec_ISet(optseen, L+1, FALSE);
    vmax = (cm->flags & CMH_LOCAL_END) ? cm->M : cm->M-1;
    optsc = ESL_MAX(Jalpha[0][L][L], ESL_MAX(Lalpha[0][L][L], ESL_MAX(Ralpha[0][L][L], inscyk_mx->Tdp[0][L][L])));
    for(v = 0; v <= vmax; v++) { 
      for(j = 1; j <= L; j++) { 
	for(d = 0; d <= j; d++) { 
	  Jsc  = Jalpha[v][j][d] + Jbeta[v][j][d] - optsc;
	  Lsc  = (v == cm->M)            ? IMPOSSIBLE : Lalpha[v][j][d] + Lbeta[v][j][d] - optsc;
	  Rsc  = (v == cm->M)            ? IMPOSSIBLE : Ralpha[v][j][d] + Rbeta[v][j][d] - optsc;
	  Tsc  = (cm->sttype[v] != B_st) ? IMPOSSIBLE : inscyk_mx->Tdp[v][j][d] + Tbeta[v][j][d] - optsc;
	  if(Jsc > 0.001) { 
	    printf("Check 1 J failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, Jalpha[v][j][d], Jbeta[v][j][d], Jalpha[v][j][d] + Jbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	  if(Lsc > 0.001) { 
	    printf("Check 1 L failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, Lalpha[v][j][d], Lbeta[v][j][d], Lalpha[v][j][d] + Lbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	  if(Rsc > 0.001) { 
	    printf("Check 1 R failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, Ralpha[v][j][d], Rbeta[v][j][d], Ralpha[v][j][d] + Rbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	  if(cm->sttype[v] == B_st && Tsc > 0.001) { 
	    printf("Check 1 T failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, inscyk_mx->Tdp[v][j][d], Tbeta[v][j][d], inscyk_mx->Tdp[v][j][d] + Tbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	  if((fabs(Jsc) < 0.001 || fabs(Lsc) < 0.001) && 
	     (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st || (cm->sttype[v] == EL_st && d >0))) { 
	    i = j-d+1;
	    /* i is accounted for by a parse with an optimal score */
	    optseen[i] = TRUE;
	    if     (fabs(Jsc) < 0.001) printf("\tResidue %4d possibly accounted for by J matrix Left  emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", i, Statetype(cm->sttype[v]), v, j, d);
	    else if(fabs(Lsc) < 0.001) printf("\tResidue %4d possibly accounted for by L matrix Left  emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", i, Statetype(cm->sttype[v]), v, j, d);
	  }
	  if((fabs(Jsc) < 0.001 || fabs(Rsc) < 0.001) && 
	     (cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st)) { 
	    /* j is accounted for by a parse with an optimal score */
	    optseen[j] = TRUE;
	    if     (fabs(Jsc) < 0.001) printf("\tResidue %4d possibly accounted for by J matrix Right emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", i, Statetype(cm->sttype[v]), v, j, d);
	    else if(fabs(Rsc) < 0.001) printf("\tResidue %4d possibly accounted for by R matrix Right emitter %2s cell [v:%4d][j:%4d][d:%4d]\n", i, Statetype(cm->sttype[v]), v, j, d);
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
  if(fail1_flag || fail2_flag) for(j = 1; j <= L; j++) printf("dsq[%4d]: %4d\n", j, dsq[j]);

  FILE *fp1; fp1 = fopen("tmp.trocykmx", "w");   cm_tr_mx_Dump(fp1, mx); fclose(fp1);
  if     (fail1_flag) ESL_FAIL(eslFAIL, errbuf, "TrCYK Inside/Outside check1 FAILED.");
  else if(fail2_flag) ESL_FAIL(eslFAIL, errbuf, "TrCYK Inside/Outside check2 FAILED.");
  else                printf("SUCCESS! TrCYK Inside/Outside checks PASSED.\n");

  ESL_DPRINTF1(("\tcm_TrCYKOutsideAlign() sc : %f (sc is from Inside!)\n", 
		ESL_MAX(Jalpha[0][L][L], ESL_MAX(Lalpha[0][L][L], ESL_MAX(Ralpha[0][L][L], inscyk_mx->Tdp[0][L][L])))));
  
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Out of memory");
}  

/* Function: cm_TrOutsideAlign()
 * Date:     EPN, Mon Sep 19 14:48:30 2011
 *
 * Purpose: Run the truncated outside algorithm. Non-banded version.
 *           A CM_TR_MX DP matrix must be passed in.  Very similar to
 *           cm_TrCYKOutsideAlign() but calculates summed log probs of
 *           all likely parses instead of the most likely parse.
 *           i.e. uses log sum operations instead of max's.  Meaning of
 *           cells:
 *
 *           Jbeta[v][j][d]: summed log prob of all parsetrees that
 *                           emit 1..i-1 and j+1..L and pass through 
 *                           v in Joint marginal mode at j,d.
 *           Lbeta[v][j][d]: summed log prob of all parsetrees that
 *                           emit 1..i-1 and j+1..L and pass through
 *                           v in Left marginal mode at j,d.
 *           Rbeta[v][j][d]: summed log prob of all parsetrees that
 *                           emit 1..i-1 and j+1..L and pass through 
 *                           v in Right marginal mode at j,d.
 *
 * Args:     cm        - the model    [0..M-1]
 *           errbuf    - char buffer for reporting errors
 *           dsq       - the digitized sequence
 *           L         - length of the dsq to align
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           do_check  - TRUE to attempt to check matrices for correctness
 *           opt_mode  - TRMODE_J, TRMODE_L, TRMODE_R, or TRMODE_T, the optimal
 *                       alignment mode, we'll only allow alignments in this mode.
 *           mx        - the dp matrix, only cells within bands in cp9b will be valid
 *           ins_mx    - the dp matrix from the CYK Inside run calculation 
 *                       (performed by cm_TrCYKInsideAlign(), required)
 *
 * Returns:  <eslOK> on success
 *
 * Throws:   <eslERANGE> if required CM_HB_MX size exceeds <size_limit>
 *           <eslFAIL>   if <do_check>==TRUE and we fail a test
 *           In either of these cases, alignment has been aborted, ret_sc is not valid.
 */
int
cm_TrOutsideAlign(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, char opt_mode,
		  CM_TR_MX *mx, CM_TR_MX *ins_mx)
{
  int      status;
  int      v,y,z;	       /* indices for states */
  int      j,d,i,k;	       /* indices in sequence dimensions */
  float    Jsc,Lsc,Rsc,Tsc;    /* a temporary variable holding a float score */
  float    optsc;              /* the Inside score */
  float    escore;	       /* an emission score, tmp variable */
  int      voffset;	       /* index of v in t_v(y) transition scores */
  /* variables used only if do_check */
  int      fail1_flag = FALSE; /* set to TRUE if do_check and we see a problem */
  int      fail2_flag = FALSE; /* set to TRUE if do_check and we see a problem */
  int      vmax;               /* i, offset in the matrix */
  /* indices used in the depths of the DP recursion */
  int      emitmode;           /* EMITLEFT, EMITRIGHT, EMITPAIR, EMITNONE, for state y */
  int      sd;                 /* StateDelta(cm->sttype[y]) */
  int      sdl;                /* StateLeftDelta(cm->sttype[y] */
  int      sdr;                /* StateRightDelta(cm->sttype[y] */

  /* DP matrix variables */
  float ***Jbeta   = mx->Jdp;     /* pointer to the outside Jbeta DP matrix */
  float ***Lbeta   = mx->Ldp;     /* pointer to the outside Lbeta DP matrix */
  float ***Rbeta   = mx->Rdp;     /* pointer to the outside Rbeta DP matrix */
  float ***Tbeta   = mx->Tdp;     /* pointer to the outside Tbeta DP matrix */

  float ***Jalpha  = ins_mx->Jdp; /* pointer to the precalc'ed inside Jalpha DP matrix */
  float ***Lalpha  = ins_mx->Ldp; /* pointer to the precalc'ed inside Lalpha DP matrix */
  float ***Ralpha  = ins_mx->Rdp; /* pointer to the precalc'ed inside Ralpha DP matrix */
  /* Talpha is not used */

  /* Allocations and initializations */
  /* grow the matrices based on the current sequence and bands */
  if((status = cm_tr_mx_GrowTo(cm, mx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the matrix to IMPOSSIBLE */
  if(  mx->Jncells_valid   > 0) esl_vec_FSet(mx->Jdp_mem, mx->Jncells_valid, IMPOSSIBLE);
  if(  mx->Lncells_valid   > 0) esl_vec_FSet(mx->Ldp_mem, mx->Lncells_valid, IMPOSSIBLE);
  if(  mx->Rncells_valid   > 0) esl_vec_FSet(mx->Rdp_mem, mx->Rncells_valid, IMPOSSIBLE);
  if(  mx->Tncells_valid   > 0) esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid, IMPOSSIBLE); 

  /* set cells in the special ROOT_S deck corresponding to full sequence alignments to 0. */
  if     (opt_mode == TRMODE_J) Jbeta[0][L][L] = 0.; /* a full Joint    alignment is outside this cell */
  else if(opt_mode == TRMODE_L) Lbeta[0][L][L] = 0.; /* a full Left     alignment is outside this cell */
  else if(opt_mode == TRMODE_R) Rbeta[0][L][L] = 0.; /* a full Right    alignment is outside this cell */
  else if(opt_mode == TRMODE_T) Tbeta[0][L][L] = 0.; /* a full Terminal alignment is outside this cell */

  /* we also need to allow full truncated alignments, rooted 
   * at internal truncated entry/exit points.
   */ 
  for(v = 0; v < cm->M; v++) { 
    if(! StateIsDetached(cm, v)) { 
      if(opt_mode == TRMODE_J && 
	 (cm->sttype[v] == B_st  || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || 
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)) { 
	Jbeta[v][L][L] = 0.; /* a full Joint alignment is outside this cell */
      }
      if(opt_mode == TRMODE_L && 
	 (cm->sttype[v] == B_st  || cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st)) { 
	if(opt_mode == TRMODE_L) Lbeta[v][L][L] = 0.; /* a full Left  alignment is outside this cell */
      }
      if(opt_mode == TRMODE_R && 
	 (cm->sttype[v] == B_st  || cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st)) { 
	if(opt_mode == TRMODE_R) Rbeta[v][L][L] = 0.; /* a full Left  alignment is outside this cell */
      }
      if(opt_mode == TRMODE_T && 
	 cm->sttype[v] == B_st) { 
	if(opt_mode == TRMODE_T) Tbeta[v][L][L] = 0.; /* a full Terminal alignment is outside this cell */
      }
    }
  }

  /* main loop down through the decks */
  for (v = 1; v < cm->M; v++) {
    if(! StateIsDetached(cm, v)) { /* skip detached inserts, they're cells will remain IMPOSSIBLE */
      sd  = StateDelta(cm->sttype[v]);
      sdr = StateRightDelta(cm->sttype[v]);
      
      if (cm->stid[v] == BEGL_S) { /* BEGL_S */
	y = cm->plast[v];	/* the parent bifurcation    */
	z = cm->cnum[y];	/* the other (right) S state */
	for(j = 0; j <= L; j++) { 
	  for (d = 0; d <= j; d++) {
	    for (k = 0; k <= (L-j); k++) {
	      Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j+k][d+k] + Jalpha[z][j+k][k]); /* A */
	      Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Lbeta[y][j+k][d+k] + Lalpha[z][j+k][k]); /* B */
	      Rbeta[v][j][d] = FLogsum(Rbeta[v][j][d], Rbeta[y][j+k][d+k] + Jalpha[z][j+k][k]); /* C */
	      if(d == j && (j+k) == L) { 
		Rbeta[v][j][d] = FLogsum(Rbeta[v][j][d], Tbeta[y][j+k][d+k] + Lalpha[z][j+k][k]); /* D */
		/* Note: Tbeta[y][j+k==L][d+k==L] will be 0.0 because
		 * it was initialized that way. That T cell includes the
		 * full target 1..L (any valid T alignment must
		 * because we must account for the full target) rooted
		 * at a B state, and a transition from that B state to
		 * this BEGL_S is always probability 1.0.
		 */
	      }
	    }
	    Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Lbeta[y][j][d]); /* entire sequence on left, no sequence on right, k == 0 */
	    Lbeta[v][j][d] = FLogsum(Lbeta[v][j][d], Lbeta[y][j][d]); /* entire sequence on left, no sequence on right, k == 0 */
	  }
	}
      } /* end of 'if (cm->stid[v] == BEGL_S */
      else if (cm->stid[v] == BEGR_S) {
	y = cm->plast[v];	  /* the parent bifurcation    */
	z = cm->cfirst[y];  /* the other (left) S state  */
	for(j = 0; j <= L; j++) { 
	  for (d = 0; d <= j; d++) {
	    i = j-d+1;
	    for (k = 0; k <= (j-d); k++) {
	      Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j][d+k] + Jalpha[z][j-d][k]); /* A */
	      Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Rbeta[y][j][d+k] + Ralpha[z][j-d][k]); /* C */
	      Lbeta[v][j][d] = FLogsum(Lbeta[v][j][d], Lbeta[y][j][d+k] + Jalpha[z][j-d][k]); /* B */
	      if(k == (i-1) && j == L) { 
		Lbeta[v][j][d] = FLogsum(Lbeta[v][j][d], Tbeta[y][j][d+k] + Ralpha[z][j-d][k]); /* D */
		/* Note: Tbeta[y][j==L][d+k==L] will be 0.0 because it
		 * was initialized that way. That T cell includes the
		 * full target 1..L (any valid T alignment must
		 * because we must account for the full target) rooted
		 * at a B state, and a transition from that B state to
		 * this BEGR_S is always probability 1.0.
		 */
	      }
	    }
	    Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Rbeta[y][j][d]); /* entire sequence on right, no sequence on left, k == 0 */
	    Rbeta[v][j][d] = FLogsum(Rbeta[v][j][d], Rbeta[y][j][d]); /* entire sequence on right, no sequence on left, k == 0 */
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
	      sdl = StateLeftDelta(cm->sttype[y]);
	      sdr = StateRightDelta(cm->sttype[y]);
	      switch(cm->sttype[y]) {
	      case MP_st: 
		if(j != L && d != j) { 
		  escore = cm->oesc[y][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		}
		if(j == L && d != j)  { /* only allow transition from L if we haven't emitted any residues rightwise (j==L) */
		  escore = cm->lmesc[y][dsq[i-1]];
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Lbeta[y][j][d+sdl]     + cm->tsc[y][voffset] + escore);
		  Lbeta[v][j][d] = FLogsum(Lbeta[v][j][d], Lbeta[y][j][d+sdl]     + cm->tsc[y][voffset] + escore);
		}
		if(i == 1 && j != L) { /* only allow transition from R if we haven't emitted any residues leftwise (i==1) */
		  escore = cm->rmesc[y][dsq[j+1]];
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Rbeta[y][j+sdr][d+sdr] + cm->tsc[y][voffset] + escore);
		  Rbeta[v][j][d] = FLogsum(Rbeta[v][j][d], Rbeta[y][j+sdr][d+sdr] + cm->tsc[y][voffset] + escore);
		}
		break;
	      case ML_st:
	      case IL_st: 
		if (d != j) { 
		  escore = cm->oesc[y][dsq[i-1]];
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j][d+sd]     + cm->tsc[y][voffset] + escore);
		  Lbeta[v][j][d] = FLogsum(Lbeta[v][j][d], Lbeta[y][j][d+sd]     + cm->tsc[y][voffset] + escore);
		}
		if(i == 1 && /* only allow transition from R if we haven't emitted any residues leftwise (i==1) */
		   v != y ) { /* will only happen if v == IL, we don't allow silent self transitions from IL->IL */
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Rbeta[y][j][d]        + cm->tsc[y][voffset]);
		  Rbeta[v][j][d] = FLogsum(Rbeta[v][j][d], Rbeta[y][j][d]        + cm->tsc[y][voffset]);
		}
		break;
	      case MR_st:
	      case IR_st:
		if (j != L) { 
		  escore = cm->oesc[y][dsq[j+1]];
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		  Rbeta[v][j][d] = FLogsum(Rbeta[v][j][d], Rbeta[y][j+sdr][d+sd] + cm->tsc[y][voffset] + escore);
		}
		if(j == L && /* only allow transition from R if we haven't emitted any residues rightwise (j==L) */
		   v != y) {  /* will only happen if v == IR, we don't allow silent self transitions from IR->IR */
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Lbeta[y][j][d]        + cm->tsc[y][voffset]);
		  Lbeta[v][j][d] = FLogsum(Lbeta[v][j][d], Lbeta[y][j][d]        + cm->tsc[y][voffset]);
		}
		break;
	      case S_st:
	      case E_st:
	      case D_st:
		if(y != 0) {
		  Jbeta[v][j][d] = FLogsum(Jbeta[v][j][d], Jbeta[y][j][d] + cm->tsc[y][voffset]);
		  Lbeta[v][j][d] = FLogsum(Lbeta[v][j][d], Lbeta[y][j][d] + cm->tsc[y][voffset]);
		  Rbeta[v][j][d] = FLogsum(Rbeta[v][j][d], Rbeta[y][j][d] + cm->tsc[y][voffset]);
		  /* y can't be 0 because transitions from ROOT_S are
		   * free in truncated mode, and we've already handled
		   * them by initializing all valid full alignments to 0.0
		   */
		}
		break;
	      } /* end of switch(cm->sttype[y] */  
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (Jbeta[v][j][d] < IMPOSSIBLE) Jbeta[v][j][d] = IMPOSSIBLE;
	  } /* ends loop over d. We know all beta[v][j][d] in this row j and state v */
	} /* end loop over j. We know beta for this whole state */
      } /* end of 'else' (if cm->sttype[v] != BEGL_S, BEGR_S) */
    } /* end of 'if(! StateIsDetached(cm, v))' */
    /* we're done calculating deck v for everything but local ends */
      
    /* deal with local end transitions v->EL J matrix only (EL = deck at M.) */
    if ((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) {
      sdr = StateRightDelta(cm->sttype[v]); /* note sdr is for state v */
      sd  = StateDelta(cm->sttype[v]);      /* note sd  is for state v */
      emitmode = Emitmode(cm->sttype[v]);   /* note emitmode is for state v */
      
      for (j = 0; j <= L; j++) { 
	for (d = 0; d <= j; d++) {
	  i = j-d+1;
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    if (j == L || d == j) continue; /* boundary condition */
	    escore = cm->oesc[v][dsq[i-1]*cm->abc->Kp+dsq[j+1]];
	    Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case ML_st:
	  case IL_st:
	    if (d == j) continue;	
	    escore = cm->oesc[v][dsq[i-1]];
	    Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case MR_st:
	  case IR_st:
	    if (j == L) continue;
	    escore = cm->oesc[v][dsq[j+1]];
	    Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], (Jbeta[v][j+sdr][d+sd] + cm->endsc[v] + escore));
	    break;
	  case S_st:
	  case E_st:
	  case D_st:
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
    for (j = L; j > 0; j--) { /* careful w/ boundary here */
      for (d = j-1; d >= 0; d--) /* careful w/ boundary here */
	Jbeta[cm->M][j][d] = FLogsum(Jbeta[cm->M][j][d], Jbeta[cm->M][j][d+1] + cm->el_selfsc);
    }
  }

  fail1_flag = FALSE;
  fail2_flag = FALSE;
  printf("DO CHECK: %d\n", do_check);
  if(do_check) {
    /* Check for consistency between the Inside alpha matrix and the
     * Outside beta matrix. If the Inside score (optsc) really is
     * the log sum of all possible parsetrees that emit the full
     * target sequence 1..L, then for all v,j,d: 
     * 
     * Jalpha[v][j][d] + Jbeta[v][j][d] <= optsc
     * Lalpha[v][j][d] + Lbeta[v][j][d] <= optsc
     * Ralpha[v][j][d] + Rbeta[v][j][d] <= optsc
     * Talpha[v][j][d] + Tbeta[v][j][d] <= optsc
     *
     * We do a more extensive check in cm_TrCYKOutsideAlign(), but
     * it doesn't apply here, because we've summed all parsetrees
     * instead of finding only the optimal one.
     *
     * This is an expensive check and should only be done while
     * debugging.
     */
    vmax = (cm->flags & CMH_LOCAL_END) ? cm->M : cm->M-1;
    optsc = ESL_MAX(Jalpha[0][L][L], ESL_MAX(Lalpha[0][L][L], ESL_MAX(Ralpha[0][L][L], ins_mx->Tdp[0][L][L])));
    for(v = 0; v <= vmax; v++) { 
      for(j = 1; j <= L; j++) { 
	for(d = 0; d <= j; d++) { 
	  Jsc  = Jalpha[v][j][d] + Jbeta[v][j][d] - optsc;
	  Lsc  = (v == cm->M) ? IMPOSSIBLE : Lalpha[v][j][d] + Lbeta[v][j][d] - optsc;
	  Rsc  = (v == cm->M) ? IMPOSSIBLE : Ralpha[v][j][d] + Rbeta[v][j][d] - optsc;
	  Tsc  = (cm->sttype[v] != B_st) ? IMPOSSIBLE : ins_mx->Tdp[v][j][d] + Tbeta[v][j][d] - optsc;
	  if(Jsc > 0.001) { 
	    printf("Check J failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, Jalpha[v][j][d], Jbeta[v][j][d], Jalpha[v][j][d] + Jbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	  if(Lsc > 0.001) { 
	    printf("Check L failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, Lalpha[v][j][d], Lbeta[v][j][d], Lalpha[v][j][d] + Lbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	  if(Rsc > 0.001) { 
	    printf("Check R failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, Ralpha[v][j][d], Rbeta[v][j][d], Ralpha[v][j][d] + Rbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	  if(cm->sttype[v] == B_st && Tsc > 0.001) { 
	    printf("Check 1 T failure: v: %4d j: %4d d: %4d (%.4f + %.4f) %.4f > %.4f\n", 
		   v, j, d, ins_mx->Tdp[v][j][d], Tbeta[v][j][d], ins_mx->Tdp[v][j][d] + Tbeta[v][j][d], optsc);
	    fail1_flag = TRUE;
	  }
	}
      }
    }
    if(fail1_flag) for(j = 1; j <= L; j++) printf("dsq[%4d]: %4d\n", j, dsq[j]);
  }
  FILE *fp1; fp1 = fopen("tmp.tromx", "w");   cm_tr_mx_Dump(fp1, mx); fclose(fp1);

  if     (fail1_flag) ESL_FAIL(eslFAIL, errbuf, "Tr Inside/Outside check1 FAILED.");
  else if(fail2_flag) ESL_FAIL(eslFAIL, errbuf, "Tr Inside/Outside check2 FAILED.");
  else                printf("SUCCESS! Tr Inside/Outside checks PASSED.\n");

  ESL_DPRINTF1(("\tTrOutsideAlign() sc : %f (sc is from Inside!)\n", 
		ESL_MAX(Jalpha[0][L][L], ESL_MAX(Lalpha[0][L][L], ESL_MAX(Ralpha[0][L][L], ins_mx->Tdp[0][L][L])))));

  return eslOK;
}

/* Function: cm_TrPosterior() 
 * Date:     EPN, Tue Sep 13 16:18:25 2011
 * Note:     based on Ian Holmes' P7EmitterPosterior() from HMMER's 2.x postprob.c
 *
 * Purpose: Combines non-banded cm_TrInside and cm_TrOutside matrices
 *           into a posterior probability matrix. The value in
 *           post->{J,L,R}[v][j][d] is the log of the posterior
 *           probability of a parse subtree rooted at v emitting the
 *           subsequence i..j (i=j-d+1) and being in J, L, or R mode
 *           at at state v.  The caller must provide a <post> float
 *           matrix, but this matrix may be the same matrix as that
 *           provided as Outside <out_mx>, (overwriting it will not
 *           compromise the algorithm).
 *
 * Args:     cm         - the model
 *           errbuf     - char buffer for reporting errors
 *           L          - length of the target sequence we're aligning
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           opt_mode   - mode of optimal alignment: TRMODE_J, TRMODE_L, TRMODE_R, or TRMODE_T
 *           ins_mx     - pre-calculated Inside matrix 
 *           out_mx     - pre-calculated Outside matrix
 *           post_mx    - pre-allocated matrix for Posteriors 
 *
 * Return:   eslOK on success, eslEINCOMPAT on contract violation
 */
int
cm_TrPosterior(CM_t *cm, char *errbuf, int L, float size_limit, char opt_mode, CM_TR_MX *ins_mx, CM_TR_MX *out_mx, CM_TR_MX *post_mx)
{
  int   status;
  int   v, j, d, i;
  float sc;
  int   fill_L, fill_R, fill_T; /* should we fill-in values for L, R, T? (we always fill in J) */

  /* Determine which matrices we need to fill-in, and the optimal score */
  if((status = cm_TrFillFromMode(opt_mode, &fill_L, &fill_R, &fill_T)) != eslOK) 
    ESL_FAIL(status, errbuf, "cm_TrPosterior, bogus opt_mode: %d", opt_mode);
  if(opt_mode == TRMODE_J) sc = ins_mx->Jdp[0][L][L];
  if(opt_mode == TRMODE_L) sc = ins_mx->Ldp[0][L][L];
  if(opt_mode == TRMODE_R) sc = ins_mx->Rdp[0][L][L];
  if(opt_mode == TRMODE_T) sc = ins_mx->Tdp[0][L][L];

  printf("HEYOO in cm_TrPosterior() mode: %d sc: %.4f fill{L,R,T}: %d%d%d\n", opt_mode, sc, fill_L, fill_R, fill_T);

  /* If local ends are on, start with the EL state (cm->M), otherwise
   * it's not a valid deck. */
  if(cm->flags & CMH_LOCAL_END) { 
    for (j = L; j >= 0; j--) {
      for (d = 0; d <= j; d++, i--) {
	post_mx->Jdp[cm->M][j][d] = ins_mx->Jdp[cm->M][j][d] + out_mx->Jdp[cm->M][j][d] - sc;
      }
    }
  }
  /* Fill in the rest of the matrices */
  for (v = cm->M-1; v >= 0; v--) { 
    for (j = L; j >= 0; j--) {
      i = j;
      for (d = 0; d <= j; d++, i--) {
	post_mx->Jdp[v][j][d] = ins_mx->Jdp[v][j][d] + out_mx->Jdp[v][j][d] - sc;
      }
    }
  }
  if (fill_L) { 
    for (v = cm->M-1; v >= 0; v--) { 
      for (j = L; j >= 0; j--) {
	i = j;
	for (d = 0; d <= j; d++, i--) {
	  post_mx->Ldp[v][j][d] = ins_mx->Ldp[v][j][d] + out_mx->Ldp[v][j][d] - sc;
	}
      }
    }
  }
  if (fill_R) { 
    for (v = cm->M-1; v >= 0; v--) { 
      for (j = L; j >= 0; j--) {
	i = j;
	for (d = 0; d <= j; d++, i--) {
	  post_mx->Rdp[v][j][d] = ins_mx->Rdp[v][j][d] + out_mx->Rdp[v][j][d] - sc;
	}
      }
    }
  }
  if (fill_T) { 
    for (v = cm->M-1; v >= 0; v--) { 
      if (v == 0 || cm->sttype[v] == B_st) { 
	for (j = L; j >= 0; j--) {
	  i = j;
	  for (d = 0; d <= j; d++, i--) {
	    post_mx->Tdp[v][j][d] = ins_mx->Tdp[v][j][d] + out_mx->Tdp[v][j][d] - sc;
	  }
	}
      }
    }
  }
  FILE *fp1; fp1 = fopen("tmp.trpmx", "w");   cm_tr_mx_Dump(fp1, post_mx); fclose(fp1);
  return eslOK;
}

/* Function: cm_TrModeFromAlphas()
 * Date:     EPN, Wed Sep 28 04:56:29 2011
 *
 * Purpose: Given a filled Inside truncated alpha matrix <mx> (a
 *           CM_TR_MX), determine the optimal mode as the mode (J,L,R
 *           or T) that has the highest mx->{J,L,R,T}dp[0][W][W] value
 *           and return it in <*ret_mode> and return the optimal score
 *           in <*ret_sc>. 
 *              
 *           If all of mx->{J,L,R,T}dp[0][L][L] have IMPOSSIBLE scores
 *           return TRMODE_UNKNOWN in <ret_mode>.
 *
 *           Important: this function only works on alpha matrices
 *           filled by 'alignment' TrCYK/TrInside functions where the
 *           full target 1..L must be accounted for in every valid
 *           parsetree. It is not designed for use with matrices
 *           filled by 'scanning' TrCYK/TrInside functions that allow
 *           valid parsetrees that only emit a subsequence of the full
 *           target.
 *
 * Args:     mx         - the Inside matrix, pre-filled.
 *           L          - length of the sequence mx was filled for
 *           ret_mode   - RETURN: optimal mode
 *           ret_sc     - RETURN: optimal score
 */
void
cm_TrModeFromAlphas(CM_TR_MX *mx, int L, char *ret_mode, float *ret_sc)
{
  char mode; 
  float sc = IMPOSSIBLE;

  sc     = mx->Jdp[0][L][L];
  mode   = TRMODE_J;

  if(mx->Ldp[0][L][L] > sc) { 
    sc     = mx->Ldp[0][L][L];
    mode   = TRMODE_L;
  }
  if(mx->Rdp[0][L][L] > sc) { 
    sc     = mx->Rdp[0][L][L];
    mode   = TRMODE_R;
  }
  if(mx->Tdp[0][L][L] > sc) { 
    sc     = mx->Tdp[0][L][L];
    mode   = TRMODE_T;
  }

  /* check for special case that all (allowed) scores were IMPOSSIBLE */
  if(! NOT_IMPOSSIBLE(sc)) { 
    mode   = TRMODE_UNKNOWN;
  }

  if(ret_sc     != NULL) *ret_sc     = sc;
  if(ret_mode   != NULL) *ret_mode   = mode;

  return;
}

/* Function: cm_TrFillFromMode()
 * Date:     EPN, Wed Sep 28 05:29:19 2011
 *
 * Purpose: Given an optimal marginal alignment mode from a 
 *          truncated matrix filled in the Inside direction (either
 *          CYK or Inside), determine which of the marginal matrices 
 *          in a truncated Outside or Posterior matrix we need to 
 *          fill in to accompany that Inside matrix.
 *
 *          If mode == TRMODE_J: fill J matrix only
 *          If mode == TRMODE_L: fill J and L matrices only
 *          If mode == TRMODE_R: fill J and R  matrices only
 *          If mode == TRMODE_T: fill J, L, R, and T matrices
 *          If mode == TRMODE_UNKNOWN: fill J, L, R, and T matrices
 *
 *          Return TRUE/FALSE values in <ret_fill_{L,R,T}>.
 *          Note that we always must fill in J matrices so a fill_J 
 *          value is unnecessary, it's implicitly true.
 *
 * Args:     mode       - optimal mode
 *           ret_fill_L - RETURN: should we fill in L based on <ret_mode>?
 *           ret_fill_R - RETURN: should we fill in R based on <ret_mode>?
 *           ret_fill_T - RETURN: should we fill in T based on <ret_mode>?
 *
 * Throws:   eslEINVAL if mode is not TRMODE_J, TRMODE_L, TRMODE_R, TRMODE_T nor TRMODE_UNKNOWN.
 */
int
cm_TrFillFromMode(char mode, int *ret_fill_L, int *ret_fill_R, int *ret_fill_T)
{
  int fill_L, fill_R, fill_T;
  int invalid_mode = FALSE;

  fill_L = fill_R = fill_T = FALSE;
  switch(mode) {
  case TRMODE_J:
    break;
  case TRMODE_L:
    fill_L = TRUE;
    break;
  case TRMODE_R:
    fill_R = TRUE;
    break;
  case TRMODE_T:
  case TRMODE_UNKNOWN:
    fill_L = fill_R = fill_T = TRUE;
    break;
  default: 
    invalid_mode = TRUE;
    break;
  }

  if(ret_fill_L != NULL) *ret_fill_L = fill_L;
  if(ret_fill_R != NULL) *ret_fill_R = fill_R;
  if(ret_fill_T != NULL) *ret_fill_T = fill_T;

  if(invalid_mode) return eslEINVAL;
  return eslOK;
}

/* Function: cm_TrEmitterPosterior()
 * Date:     EPN, Fri Oct  7 05:30:31 2011
 *
 * Purpose: Given a posterior probability cube, where the value in
 *           post[v][j][d] is the log of the posterior probability of
 *           a parse subtree rooted at v emitting the subsequence i..j
 *           (i=j-d+1), fill a CM_EMIT_MX <emit_mx> with 
 *           matrices with values:
 *
 *           emit_mx->*l_pp[v][i]: log of the posterior probability
 *           that state v emitted residue i leftwise while in * (J or
 *           L, Joint of Left) marginal mode either at (if a match
 *           state) or *after* (if an insert state) the left consensus
 *           position modeled by state v's node. 
 *
 *           emit_mx->*r_pp[v][i]: log of the posterior probability
 *           that state v emitted residue i rightwise while in * (J
 *           or R, Joint or Right) marginal mode either at (if a match
 *           state) or *before* (if an insert state) the right
 *           consensus position modeled by state v's node.
 *
 *           *l_pp[v] is NULL for states that do not emit leftwise 
 *           *r_pp[v] is NULL for states that do not emit rightwise
 *
 *           We only need to fill a subset of the *l_pp and *r_pp
 *           matrices, depending on the <mode> of the alignment
 *           while is known and passed in:
 *           <mode> == TRMODE_J, fill Jl_pp, Jr_pp
 *           <mode> == TRMODE_L, fill Jl_pp, Jr_pp and Ll_pp
 *           <mode> == TRMODE_R, fill Jl_pp, Jr_pp and Rr_pp
 *           <mode> == TRMODE_T, fill Jl_pp, Jr_pp, Ll_pp, and Rr_pp
 * 
 *          This is done in 3 steps:
 *          1. Fill *l_pp[v][i] and *r_pp[v][i] with the posterior
 *             probability that state v emitted residue i either
 *             leftwise (l_pp) or rightwise (r_pp).
 *
 *          2. Normalize *l_pp and *r_pp so that probability that
 *             each residue was emitted by any state is exactly
 *             1.0.
 *
 *          3. Combine *l_pp values for MATP_MP (v) and MATP_ML (y=v+1)
 *             states in the same node so they give the value defined
 *             above (i.e. *l_pp[v] == *l_pp[y] = the PP that either v
 *             or y emitted residue i) instead of *l_pp[v] = PP that v
 *             emitted i, and *l_pp[y] = PP that y emitted i.  And
 *             combine *r_pp values for MATP_MP (v) and MATP_MR (y=v+2)
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
 * Args:     cm         - the model
 *           errbuf     - for error messages
 *           L          - length of the sequence
 *           size_limit - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           post       - pre-filled posterior cube
 *           emit_mx    - pre-allocated emit matrix, grown and filled-in here
 *           opt_mode   - known optimal mode of the alignment 
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
cm_TrEmitterPosterior(CM_t *cm, char *errbuf, int L, float size_limit, CM_TR_MX *post, CM_TR_EMIT_MX *emit_mx, char opt_mode, int do_check)
{
  int    status;
  int    v, j, d; /* state, position, subseq length */
  int    i;       /* sequence position */
  int    sd;      /* StateDelta(v) */
  int    fill_L, fill_R; /* do we need to fill Ll_pp/Rr_pp matrices? */
  
  /* grow the emit matrices based on the current sequence */
  if((status = cm_tr_emit_mx_GrowTo(cm, emit_mx, errbuf, L, size_limit)) != eslOK) return status;

  /* initialize all cells of the emit matrices to IMPOSSIBLE */
  esl_vec_FSet(emit_mx->Jl_pp_mem, emit_mx->l_ncells_valid, IMPOSSIBLE);
  esl_vec_FSet(emit_mx->Ll_pp_mem, emit_mx->l_ncells_valid, IMPOSSIBLE);
  esl_vec_FSet(emit_mx->Jr_pp_mem, emit_mx->r_ncells_valid, IMPOSSIBLE);
  esl_vec_FSet(emit_mx->Rr_pp_mem, emit_mx->r_ncells_valid, IMPOSSIBLE);

  if((status = cm_TrFillFromMode(opt_mode, &fill_L, &fill_R, NULL)) != eslOK) ESL_FAIL(status, errbuf, "cm_TrCheckFromPosterior, bogus mode: %d", opt_mode);
  printf("in cm_TrEmitterPosterior(), mode: %d fill_L: %d fill_R: %d\n", 
	 opt_mode, fill_L, fill_R);
     
  /* Step 1. Fill *l_pp[v][i] and *r_pp[v][i] with the posterior
   *         probability that state v emitted residue i either
   *         leftwise (*l_pp) or rightwise (*r_pp).
   */
  for(v = 0; v < cm->M; v++) { 
    sd = StateDelta(cm->sttype[v]);
    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
      for(j = 1; j <= L; j++) { 
	i = j-sd+1;
	for(d = sd; d <= j; d++, i--) { 
	  emit_mx->Jl_pp[v][i] = FLogsum(emit_mx->Jl_pp[v][i], post->Jdp[v][j][d]);
	}
	if(fill_L) { 
	  for(d = sd; d <= j; d++, i--) { 
	    emit_mx->Ll_pp[v][i] = FLogsum(emit_mx->Ll_pp[v][i], post->Ldp[v][j][d]);
	  }
	}
      }
    }
    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
      for(j = 1; j <= L; j++) { 
	sd = StateDelta(cm->sttype[v]);
	for(d = sd; d <= j; d++) { 
	  emit_mx->Jr_pp[v][j] = FLogsum(emit_mx->Jr_pp[v][j], post->Jdp[v][j][d]);
	}
	if(fill_R) { 
	  for(d = sd; d <= j; d++) { 
	    emit_mx->Rr_pp[v][j] = FLogsum(emit_mx->Rr_pp[v][j], post->Rdp[v][j][d]);
	  }
	}
      }
    }
  }
  /* factor in contribution of local ends, the EL state may have emitted this residue. */
  if (cm->flags & CMH_LOCAL_END) {
    for (j = 1; j <= L; j++) { 
      i = j;
      for (d = 1; d <= j; d++, i--) { /* note: d >= 1, b/c EL emits 1 residue */
	emit_mx->Jl_pp[cm->M][i] = FLogsum(emit_mx->Jl_pp[cm->M][i], post->Jdp[cm->M][j][d]);
      }
    }
  }
  FILE *fp1; fp1 = fopen("tmp.trunc_unnorm_emitmx",  "w"); cm_tr_emit_mx_Dump(fp1, cm, emit_mx); fclose(fp1);

  /* Step 2. Normalize *l_pp and *r_pp so that probability that
   *         each residue was emitted by any state is exactly
   *         1.0.
   */
  esl_vec_FSet(emit_mx->sum, (L+1), IMPOSSIBLE);
  for(v = 0; v <= cm->M; v++) { 
    if(emit_mx->Jl_pp[v] != NULL) {
      for(i = 1; i <= L; i++) { 
	emit_mx->sum[i] = FLogsum(emit_mx->sum[i], emit_mx->Jl_pp[v][i]);
      }
    }
    if(emit_mx->Ll_pp[v] != NULL && fill_L) { 
      for(i = 1; i <= L; i++) { 
	emit_mx->sum[i] = FLogsum(emit_mx->sum[i], emit_mx->Ll_pp[v][i]);
      }
    }
    if(emit_mx->Jr_pp[v] != NULL) {
      for(j = 1; j <= L; j++) { 
	emit_mx->sum[j] = FLogsum(emit_mx->sum[j], emit_mx->Jr_pp[v][j]);
      }
    }
    if(emit_mx->Rr_pp[v] != NULL && fill_R) { 
      for(j = 1; j <= L; j++) { 
	emit_mx->sum[j] = FLogsum(emit_mx->sum[j], emit_mx->Rr_pp[v][j]);
      }
    }
  }
  /* perform the check, if nec */
  if(do_check) { 
    for(i = 1; i <= L; i++) { 
      if((sreEXP2(emit_mx->sum[i]) < 0.98) || (sreEXP2(emit_mx->sum[i]) > 1.02)) { 
	ESL_FAIL(eslFAIL, errbuf, "residue %d has summed prob of %5.4f (2^%5.4f).\nMay not be a DP coding bug, see 'Note:' on precision in cm_TrEmitterPosterior().\n", i, (sreEXP2(emit_mx->sum[i])), emit_mx->sum[i]);
      }
      printf("i: %d | total: %10.4f\n", i, (sreEXP2(emit_mx->sum[i])));
    }
    ESL_DPRINTF1(("cm_TrEmitterPosterior() check passed, all residues have summed probability of emission of between 0.98 and 1.02.\n"));
  }  

  /* normalize, using the sum vector */
  for(v = 0; v <= cm->M; v++) { 
    if(emit_mx->Jl_pp[v] != NULL) {
      for(i = 1; i <= L; i++) emit_mx->Jl_pp[v][i] -= emit_mx->sum[i];
    }
    if(emit_mx->Ll_pp[v] != NULL && fill_L) { 
      for(i = 1; i <= L; i++) emit_mx->Ll_pp[v][i] -= emit_mx->sum[i];
    }
    if(emit_mx->Jr_pp[v] != NULL) {
      for(j = 1; j <= L; j++) emit_mx->Jr_pp[v][j] -= emit_mx->sum[j];
    }
    if(emit_mx->Rr_pp[v] != NULL && fill_R) { 
      for(j = 1; j <= L; j++) emit_mx->Rr_pp[v][j] -= emit_mx->sum[j];
    }
  }

  /* Step 3. Combine *l_pp values for MATP_MP (v) and MATP_ML (y=v+1)
   *         states in the same node so they give the value defined
   *         above (i.e. *l_pp[v] == *l_pp[y] = the PP that either v or
   *         y emitted residue i) instead of *l_pp[v] = PP that v
   *         emitted i, and *l_pp[y] = PP that y emitted i.  And
   *         combine *r_pp values for MATP_MP (v) and MATP_MR (y=v+2)
   *         states in an analogous way.
   */
  for(v = 0; v <= cm->M; v++) { 
    if(cm->sttype[v] == MP_st) { 
      for(i = 1; i <= L; i++) { 
	emit_mx->Jl_pp[v][i]   = FLogsum(emit_mx->Jl_pp[v][i], emit_mx->Jl_pp[v+1][i]); 
	emit_mx->Jl_pp[v+1][i] = emit_mx->Jl_pp[v][i];
      }
      if(fill_L) { 
	for(i = 1; i <= L; i++) { 
	  emit_mx->Ll_pp[v][i]   = FLogsum(emit_mx->Ll_pp[v][i], emit_mx->Ll_pp[v+1][i]); 
	  emit_mx->Ll_pp[v+1][i] = emit_mx->Ll_pp[v][i];
	}
      }
      for(j = 1; j <= L; j++) { 
	emit_mx->Jr_pp[v][j]   = FLogsum(emit_mx->Jr_pp[v][j], emit_mx->Jr_pp[v+2][j]); 
	emit_mx->Jr_pp[v+2][j] = emit_mx->Jr_pp[v][j];
      }
      if(fill_R){ 
	for(j = 1; j <= L; j++) { 
	  emit_mx->Rr_pp[v][j]   = FLogsum(emit_mx->Rr_pp[v][j], emit_mx->Rr_pp[v+2][j]); 
	  emit_mx->Rr_pp[v+2][j] = emit_mx->Rr_pp[v][j];
	}
      }
    }
  }
  FILE *fp2; fp2 = fopen("tmp.trunc_emitmx",  "w"); cm_tr_emit_mx_Dump(fp2, cm, emit_mx); fclose(fp2);

  return eslOK;
}

/* Function: cm_TrPostCode()
 * Date:     EPN, Fri Oct  7 14:30:32 2011
 *
 * Purpose: Given a parse tree and a filled emit matrix calculate two
 *           strings that represents the confidence values on each
 *           aligned residue in the sequence.
 *           
 *           The emit_mx values are:
 *           {J,L}l_pp[v][i]: log of the posterior probability that state v emitted
 *                            residue i leftwise either at (if a match state) or
 *                            *after* (if an insert state) the left consensus
 *                            position modeled by state v's node in Joint marginal
 *                            mode (for Jl_pp) or Left marginal mode (for Ll_pp).
 *
 *           {J,R}r_pp[v][i]: log of the posterior probability that state v emitted
 *                            residue i rightwise either at (if a match state) or
 *                            *before* (if an insert state) the right consensus
 *                            position modeled by state v's node in Joint marginal
 *                            mode (for Jr_pp) or Right marginal mode (for Rr_pp).
 *
 *           {J,L}l_pp[v] is NULL for states that do not emit leftwise  (B,S,D,E,IR,MR)
 *           {J,R}r_pp[v] is NULL for states that do not emit rightwise (B,S,D,E,IL,ML)
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
 *           cm_TrPostCodeHB() is nearly the same function with the
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
 *                       or if we get a marginal mode in the parsetree that doesn't
 *                       make sense.
 */
int
cm_TrPostCode(CM_t *cm, char *errbuf, int L, CM_TR_EMIT_MX *emit_mx, Parsetree_t *tr, char **ret_ppstr, float *ret_avgp)
{
  int   status;
  int   x, v, i, j, d, r; /* counters */
  char *ppstr;       /* the PP string, created here */
  float p;           /* a probability */
  float sum_logp;    /* log of summed probability of all residues emitted thus far */
  float cur_log_pp;  /* current log probability of emitting a residue */
  char  mode;        /* marginal mode: TRMODE_J, TRMODE_L or TRMODE_R */

  ESL_ALLOC(ppstr, (L+1) * sizeof(char)); 
  sum_logp = IMPOSSIBLE;

  /* go through each node of the parsetree and determine post code for emissions */
  for (x = 0; x < tr->n; x++) { 
    v    = tr->state[x];
    i    = tr->emitl[x];
    j    = tr->emitr[x];
    d    = j-i+1;
    mode = tr->mode[x];

    /* Only P, L, R, and EL states have emissions. */
    if(cm->sttype[v] == EL_st) { /* EL state, we have to handle this guy special */
      if(mode == TRMODE_J) {
	for(r = i; r <= j; r++) { /* we have to annotate from residues i..j */
	  cur_log_pp = emit_mx->Jl_pp[v][r];
	  ppstr[r-1] = Fscore2postcode(cur_log_pp);
	  sum_logp   = FLogsum(sum_logp, cur_log_pp);
	  /* make sure we've got a valid probability */
	  p = FScore2Prob(cur_log_pp, 1.);
	  if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): probability for EL state v: %d residue r: %d > 1.00 (%.2f)", v, r, p);
	  if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): probability for EL state v: %d residue r: %d < 0.00 (%.2f)", v, r, p);
	}
      }
      else ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): invalid mode for EL state in the parsetree: %d\n", mode);
    }
    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) { 
      if(mode == TRMODE_J || mode == TRMODE_L) {  
	cur_log_pp = (mode == TRMODE_J) ? emit_mx->Jl_pp[v][i] : emit_mx->Ll_pp[v][i];
	ppstr[i-1] = Fscore2postcode(cur_log_pp);
	sum_logp   = FLogsum(sum_logp, cur_log_pp);
	/* make sure we've got a valid probability */
	p = FScore2Prob(cur_log_pp, 1.);
	if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): probability for left state v: %d residue i: %d > 1.00 (%.2f)", v, i, p);
	if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): probability for left state v: %d residue i: %d < 0.00 (%.2f)", v, i, p);
      }
      else if(mode != TRMODE_R) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): non-sensical mode for MP, ML, IL state: %d", mode);
    }
    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
      if(mode == TRMODE_J || mode == TRMODE_R) {  
	cur_log_pp = (mode == TRMODE_J) ? emit_mx->Jr_pp[v][j] : emit_mx->Rr_pp[v][j];
	ppstr[j-1] = Fscore2postcode(cur_log_pp);
	sum_logp   = FLogsum(sum_logp, cur_log_pp);
	/* make sure we've got a valid probability */
	p = FScore2Prob(cur_log_pp, 1.);
	if(p >  1.01) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): probability for right state v: %d residue i: %d > 1.00 (%.2f)", v, j, p);
	if(p < -0.01) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): probability for right state v: %d residue i: %d < 0.00 (%.2f)", v, j, p);
      }
      else if(mode != TRMODE_L) ESL_FAIL(eslEINVAL, errbuf, "cm_TrPostCode(): non-sensical mode for MP, MR, IR state: %d", mode);
    }
  }
  ppstr[L] = '\0';

  if(ret_ppstr != NULL) *ret_ppstr = ppstr; else free(ppstr);
  if(ret_avgp  != NULL) *ret_avgp  = sreEXP2(sum_logp) / (float) L;
  return eslOK;
  
 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "cm_TrPostcode(): Memory allocation error.");
  return status; /* never reached */
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
  { "--cykout",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "run TrCYKOutside, to make sure it agrees with TrCYK (Inside)", 0 },
  { "--std",     eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also do standard (non-truncated) alignments",    2},
  { "--orig",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,"--search", NULL, "also do search with original trCYK",         2},
  { "--hb",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also do HMM banded alignments",                   2},
  { "--failok",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "allow failures of Inside vs Outside checks",      2},
  { "--search",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also run search algorithms",                   2},
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "don't use QDBs", 2},
  { "--noinside",eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "don't run Inside", 2},
  { "--infile",  eslARG_INFILE,  NULL, NULL, NULL,  NULL,  NULL, "-L,-N,-e", "read sequences to search from file <s>", 2 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use posterior sums during HMM band calculation (widens bands)", 2 },
  { "--onlyhb",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only run HMM banded scanning trCYK", 2 },
  { "--oa",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "run optimal accuracy alignment instead of CYK", 2 },
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
  CM_TR_HB_MX           *out_trhbmx= NULL;
  CM_TR_HB_EMIT_MX      *trhbemit_mx = NULL;
  CM_TR_HB_SHADOW_MX    *trhbshmx= NULL;
  CM_TR_MX              *trmx= NULL;
  CM_TR_MX              *out_trmx= NULL;
  CM_TR_EMIT_MX         *tremit_mx= NULL;
  CM_TR_SHADOW_MX       *trshmx= NULL;
  float           size_limit = esl_opt_GetReal(go, "--sizelimit");
  float           save_tau, save_cp9b_thresh1, save_cp9b_thresh2;
  float           hbmx_Mb, trhbmx_Mb;
  float           parsetree_sc, parsetree_struct_sc;
  char            mode;   
  /* variables related to non-banded cyk/inside/outside */
  CM_MX             *mx   = NULL;       /* alpha DP matrix for non-banded CYK/Inside() */
  CM_MX             *out_mx = NULL;     /* outside matrix for HMM banded Outside() */
  CM_SHADOW_MX      *shmx = NULL;       /* shadow matrix for non-banded tracebacks */

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  if ((status = cm_file_Open(cmfile, NULL, FALSE, &(cmfp), errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = cm_file_Read(cmfp, TRUE, &abc, &cm))                != eslOK) cm_Fail(cmfp->errbuf);
  cm_file_Close(cmfp);

  if(esl_opt_GetBoolean(go, "--std")) { 
    mx     = cm_mx_Create(cm);
    out_mx = cm_mx_Create(cm);
    shmx   = cm_shadow_mx_Create(cm);
  }

  if((esl_opt_GetBoolean(go, "--search") && esl_opt_GetBoolean(go, "--std"))) { 
    trsmx = cm_CreateTrScanMatrix(cm, cm->W, dmax, cm->beta_W, cm->beta_qdb, 
				  (dmin == NULL && dmax == NULL) ? FALSE : TRUE,
				  TRUE, TRUE); /* do_float, do_int */
  }
    
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
  cm->align_opts |= CM_ALIGN_CHECKINOUT;

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
    trmx       = cm_tr_mx_Create(cm);
    out_trmx   = cm_tr_mx_Create(cm);
    trshmx     = cm_tr_shadow_mx_Create(cm);
    tremit_mx  = cm_tr_emit_mx_Create(cm);
    printf("done.  ");
    fflush(stdout);
    esl_stopwatch_Stop(w);
    esl_stopwatch_Display(stdout, w, " CPU time: ");
    printf("\n\n");
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

  for (i = 0; i < N; i++) { 
    if(L == 0) continue;
    L = seqs_to_aln->sq[i]->n;
    dsq = seqs_to_aln->sq[i]->dsq;
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    /* int x; for(x = 1; x <= L; x++) printf("dsq[%4d]: %4d\n", x, dsq[x]); */

    cm->tau = save_tau;
    cm->cp9b->thresh1 = save_cp9b_thresh1;
    cm->cp9b->thresh2 = save_cp9b_thresh2;

    cm->align_opts  |= CM_ALIGN_HBANDED;

    /* 1. non-banded truncated alignment, unless --onlyhb
     * 2. non-banded standard  alignment, if requested
     * 3. HMM banded truncated alignment, if requested
     * 4. HMM banded standard  alignment, if requested
     * 5. non-banded truncated search,    if requested
     * 6. non-banded standard  search,    if requested
     * 7. HMM banded truncated search,    if requested
     * 8. HMM banded standard  search,    if requested 
     */

    float cyksc = IMPOSSIBLE;
    /* 1. non-banded truncated alignment, unless --onlyhb */
    if(! esl_opt_GetBoolean(go, "--onlyhb")) { 
      /*********************Begin cm_TrAlign****************************/
      cm_tr_mx_Destroy(trmx);
      cm_tr_shadow_mx_Destroy(trshmx);
      trmx   = cm_tr_mx_Create(cm);
      trshmx = cm_tr_shadow_mx_Create(cm);
      esl_stopwatch_Start(w);
      if(esl_opt_GetBoolean(go, "--oa")) { 
	if((status = cm_TrAlign(cm, errbuf, dsq, L, size_limit, TRUE, FALSE, trmx, trshmx, out_trmx, tremit_mx, NULL, NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f PP   (FULL LENGTH OPTACC)", (i+1), "cm_TrAlign(): ", sc);
      }
      else { 
	if((status = cm_TrAlign(cm, errbuf, dsq, L, size_limit, FALSE, FALSE, trmx, trshmx, out_trmx, tremit_mx, NULL, NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits (FULL LENGTH CYK)", (i+1), "cm_TrAlign(): ", sc);
      }
      esl_stopwatch_Stop(w);
      cyksc = sc;
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
      ParsetreeScore(cm, NULL, NULL, tr, dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
      printf("Parsetree score      : %.4f           (FULL LENGTH CYK)\n", parsetree_sc);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End cm_TrAlign****************************/

      if(esl_opt_GetBoolean(go, "--cykout")) { 
	/*********************Begin cm_TrCYKOutsideAlign****************************/
	esl_stopwatch_Start(w);
	/* determine optimal mode from the trmx */
	sc   = trmx->Jdp[0][L][L];
	mode = TRMODE_J;
	if(trmx->Ldp[0][L][L] > sc) { 
	  sc   = trmx->Ldp[0][L][L];
	  mode = TRMODE_L;
	}
	if(trmx->Rdp[0][L][L] > sc) { 
	  sc   = trmx->Rdp[0][L][L];
	  mode = TRMODE_R;
	}
	if(trmx->Tdp[0][L][L] > sc) { 
	  sc   = trmx->Tdp[0][L][L];
	  mode = TRMODE_T;
	}
	status = cm_TrCYKOutsideAlign(cm, errbuf, seqs_to_aln->sq[i]->dsq,  L, size_limit, TRUE, mode, out_trmx, trmx);
	if     (status != eslOK && esl_opt_GetBoolean(go, "--failok")) printf("%s\nError detected, but continuing thanks to --failok\n", errbuf);
	else if(status != eslOK)                                       cm_Fail(errbuf);
	printf("%4d %-30s %10s bits ", (i+1), "cm_TrCYKOutsideAlign() CYK:", "?");
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*********************End cm_TrCYKOutsideAlign****************************/
      }

      if(! esl_opt_GetBoolean(go, "--noinside")) { 
	/*********************Begin cm_TrInsideAlign()****************************/
	cm_tr_mx_Destroy(trmx);
	trmx   = cm_tr_mx_Create(cm);
	esl_stopwatch_Start(w);
	if((status = cm_TrInsideAlign(cm, errbuf, dsq, L, size_limit, trmx, &mode, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits (FULL LENGTH INSIDE)", (i+1), "cm_TrInsideAlign(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	if(sc < cyksc) { 
	  int j;
	  for(j = 1; j <= L; j++) printf("dsq[%4d]: %4d\n", j, dsq[j]);
	  cm_Fail("cm_TrInsideAlign() score: %.4f < CYK score : %.4f\n", sc, cyksc);
	}
	/*********************End cm_TrInsideAlign*****************************/
	
	/*********************Begin cm_TrOutsideAlign()****************************/
	esl_stopwatch_Start(w);
	if((status = cm_TrOutsideAlign(cm, errbuf, dsq, L, size_limit, TRUE, mode, out_trmx, trmx)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10s bits (FULL LENGTH OUTSIDE)", (i+1), "cm_TrOutsideAlign(): ", "?");
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*********************End cm_TrOutsideAlign*****************************/
      }
    }

    /* 2. non-banded standard (non-truncated) alignment, if requested */
    if(esl_opt_GetBoolean(go, "--std") && (! esl_opt_GetBoolean(go, "--onlyhb"))) { 
      /*********************Begin cm_Align()****************************/
      esl_stopwatch_Start(w);
      if((status = cm_Align  (cm, errbuf, dsq, L, size_limit, FALSE, FALSE, mx, shmx, NULL, NULL, NULL, NULL, NULL, &tr, &sc)) != eslOK) return status;
      printf("%4d %-30s %10.4f bits (FULL LENGTH CYK)", (i+1), "cm_Align(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End cm_Align*****************************/

      if(esl_opt_GetBoolean(go, "--cykout")) { 
	/*********************Begin cm_CYKOutsideAlign****************************/
	esl_stopwatch_Start(w);
	status = cm_CYKOutsideAlign(cm, errbuf, seqs_to_aln->sq[i]->dsq, L, size_limit, TRUE, out_mx, mx, &sc);
	if     (status != eslOK && esl_opt_GetBoolean(go, "--failok")) printf("%s\nError detected, but continuing thanks to --failok\n", errbuf);
	else if(status != eslOK)                                       cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "cm_CYKOutsideAlign() CYK:", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	/*********************End cm_CYKOutsideAlign****************************/
      }


      /*********************Begin cm_InsideAlign()****************************/
      esl_stopwatch_Start(w);
      if((status = cm_InsideAlign (cm, errbuf, dsq, L, size_limit, mx, &sc)) != eslOK) return status;
      printf("%4d %-30s %10.4f bits (FULL LENGTH INSIDE)", (i+1), "cm_InsideAlign(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End cm_InsideAlign*****************************/
	
      /*********************Begin cm_OutsideAlign*****************************/
      esl_stopwatch_Start(w);
      if((status = cm_OutsideAlign(cm, errbuf, dsq, L, size_limit, TRUE, out_mx, mx, &sc)) != eslOK) return status;
      printf("%4d %-30s %10.4f bits (FULL LENGTH OUTSIDE)", (i+1), "cm_OutsideAlign(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End cm_OutsideAlign*****************************/
      if((status = cm_Posterior(cm, errbuf, L, size_limit, mx, out_mx, out_mx)) != eslOK) return status;   
    }
      
    /* 3. HMM banded truncated alignment, if requested */
    if(esl_opt_GetBoolean(go, "--hb") || esl_opt_GetBoolean(go, "--onlyhb")) { 
      /*********************Begin cm_TrAlignHB()****************************/
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
      cm_tr_hb_mx_Destroy(out_trhbmx);
      cm_tr_hb_emit_mx_Destroy(trhbemit_mx);
      cm_tr_hb_shadow_mx_Destroy(trhbshmx);
      trhbmx      = cm_tr_hb_mx_Create(cm);
      out_trhbmx  = cm_tr_hb_mx_Create(cm);
      trhbshmx    = cm_tr_hb_shadow_mx_Create(cm);
      trhbemit_mx = cm_tr_hb_emit_mx_Create(cm);
      esl_stopwatch_Start(w);
      if((status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, FALSE, FALSE, trhbmx, trhbshmx, out_trhbmx, trhbemit_mx, NULL, NULL, NULL, &tr, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits (FULL LENGTH CYK)", (i+1), "cm_TrAlign(): ", sc);

      printf("%4d %-30s %10.4f bits ", (i+1), "cm_TrAlignHB(): ", sc);
      esl_stopwatch_Stop(w);
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
      ParsetreeScore(cm, NULL, NULL, tr, dsq, FALSE, &parsetree_sc, &parsetree_struct_sc, NULL, NULL, NULL);
      printf("Parsetree score      : %.4f           (FULL LENGTH CYK)\n", parsetree_sc);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End cm_TrAlignHB*****************************/
	
      /*********************Begin cm_TrInsideAlignHB()****************************/
      cm_tr_hb_mx_Destroy(trhbmx);
      trhbmx = cm_tr_hb_mx_Create(cm);
      esl_stopwatch_Start(w);
      if((status = cm_TrInsideAlignHB(cm, errbuf, dsq, L, size_limit, trhbmx, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits (FULL LENGTH INSIDE)", (i+1), "cm_TrInsideAlignHB(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End cm_TrInsideAlignHB*****************************/
    }
      
    /* 4. HMM banded standard alignment, if requested */
    if(esl_opt_GetBoolean(go, "--std") && (esl_opt_GetBoolean(go, "--hb") || esl_opt_GetBoolean(go, "--onlyhb"))) { 
      /*********************Begin cm_AlignHB()***************************/
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
      if((status = cm_AlignHB(cm, errbuf, dsq, L, size_limit, FALSE, FALSE, cm->hbmx, cm->shhbmx, NULL, NULL, NULL, NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "cm_AlignHB(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End cm_AlignHB()***************************/
    }

    if(esl_opt_GetBoolean(go, "--search")) { 
      /* 5. non-banded truncated search, if requested */
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

      /* 6. non-banded standard search, if requested */
      if(cm->smx == NULL) cm_CreateScanMatrixForCM(cm, TRUE, TRUE); /* impt to do this after QDBs set up in ConfigCM() */
      /*********************Begin FastCYKScan****************************/
      esl_stopwatch_Start(w);
      if((status = FastCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, 0., NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      /*********************End FastCYKScan****************************/
	  
      /*********************Begin FastIInsideScan****************************/
      cm->search_opts |= CM_SEARCH_INSIDE;
      esl_stopwatch_Start(w);
      if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastIInsideScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      cm->search_opts &= ~CM_SEARCH_INSIDE;
      /*********************End RefITrInsideScan****************************/

      /* 7. HMM banded truncated search, if requested */
      if(esl_opt_GetBoolean(go, "--hb") || esl_opt_GetBoolean(go, "--onlyhb")) { 
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
	/*********************End FTrInsideScanHB***********************/

	/* 8. HMM banded standard search, if requested */
	if(esl_opt_GetBoolean(go, "--std") && (esl_opt_GetBoolean(go, "--hb") || esl_opt_GetBoolean(go, "--onlyhb"))) { 
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
	}
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
  if(mx != NULL)     cm_mx_Destroy(mx);
  if(out_mx != NULL) cm_mx_Destroy(out_mx);
  if(shmx != NULL)   cm_shadow_mx_Destroy(shmx);
  if(tremit_mx != NULL) cm_tr_emit_mx_Destroy(tremit_mx);

  return 0;

 ERROR:
  cm_Fail("memory allocation error");
  return 0; /* never reached */
}
#endif /*IMPL_TRUNC_ALIGN_BENCHMARK*/




/* cm_parsetree.c
 * cove 1.0: Mon May 17 09:38:14 1993
 * moved to cove 2.0, Mon Sep  6 13:34:55 1993
 * cove4: SRE 29 Feb 2000 [Seattle]
 * infernal: SRE, Fri Jul 28 08:55:47 2000 [StL]
 * SVN $Id$
 * 
 * Unlike a traceback of a normal HMM alignment, which is linear,
 * the traceback of a CM is a tree structure. Here
 * we provide support for the traceback data structure.
 * 
 * Non-BIFURC states have a NULL right branch. 
 * 
 * The pushdown stack structure has a dummy begin node, and the
 * end is signified by a final NULL ptr.
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/* Function: CreateParsetree()
 * Incept:   SRE 29 Feb 2000 [Seattle] from cove2.0 code.
 * 
 * Purpose:  Creates a parse tree structure.
 *           The first operation on a newly created tree is
 *           generally to add the root:
 *           InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 0, L-1, 0);
 * 
 * Return:   ptr to the new tree.
 */          
Parsetree_t * 
CreateParsetree(int size)
{
  int status;
  Parsetree_t *new;

  ESL_ALLOC(new, sizeof(Parsetree_t));
  new->memblock = 25;		/* allocation block size can be optimized here if you want. */
  new->nalloc   = size;
  ESL_ALLOC(new->emitl, sizeof(int) * new->nalloc);
  ESL_ALLOC(new->emitr, sizeof(int) * new->nalloc);
  ESL_ALLOC(new->state, sizeof(int) * new->nalloc);
  ESL_ALLOC(new->mode,  sizeof(char) * new->nalloc);
  ESL_ALLOC(new->nxtl,  sizeof(int) * new->nalloc);
  ESL_ALLOC(new->nxtr,  sizeof(int) * new->nalloc);
  ESL_ALLOC(new->prv,   sizeof(int) * new->nalloc);
  new->n = 0;
  new->is_std = TRUE;
  return new;
 ERROR:
  cm_Fail("ERROR allocated parsetree.\n");
  return NULL; /* never reached */
}

/* Function: GrowParsetree()
 * Incept:   SRE 1 March 2000 [Seattle]
 * 
 * Purpose:  Increase the number of available nodes in a parse tree.
 */
void
GrowParsetree(Parsetree_t *tr)
{
  int   status;
  void *tmp;
  tr->nalloc += tr->memblock;
  ESL_RALLOC(tr->emitl, tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->emitr, tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->state, tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->mode,  tmp, sizeof(char) * tr->nalloc);
  ESL_RALLOC(tr->nxtl,  tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->nxtr,  tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->prv,   tmp, sizeof(int) * tr->nalloc);
  return;
  
 ERROR:
  cm_Fail("ERROR growing parsetree.\n");
}

/* Function: FreeParsetree()
 * Incept:   SRE 1 March 2000 [Seattle]
 *
 * Purpose:  Destroy a parse tree.
 */
void
FreeParsetree(Parsetree_t *tr)
{
  free(tr->emitl);
  free(tr->emitr);
  free(tr->state);
  free(tr->mode);
  free(tr->nxtl);
  free(tr->nxtr);
  free(tr->prv);
  free(tr);
}

/* Function: SizeofParsetree()
 * Incept:   EPN, Wed Dec  2 17:26:45 2009
 * 
 * Purpose:  Return the allocated size of a parsetree in Mb.
 */
float
SizeofParsetree(Parsetree_t *tr)
{
  float Mb;
  
  Mb = 0.;
  Mb += 3 * sizeof(int); /* n, nalloc, memblock */
  Mb += tr->nalloc * (sizeof(int)) * 6; 
  Mb += tr->nalloc * (sizeof(char)) * 1; 
  Mb /= 1000000.;
  /* 6 = emitl,emitr,state,nxtl,nxtr,prv */
  /* 1 = mode */
  return Mb;
}

/* Function: InsertTraceNodewithMode()
 * Incept:   SRE 1 March 2000 [Seattle]
 * 
 * Purpose:  Insert a new node in a trace tree, attached to node y,
 *           either TRACE_LEFT_CHILD or TRACE_RIGHT_CHILD.
 *   
 *           Before:                             After:
 *                 y                                  y
 *               /   \                              /   \
 *              a     b                            n     b
 *                                                / \
 *                                               a   -
 *           The new node has index tr->n.
 *           GrowTrace() if necessary.
 *           The new node n gets connectivity:
 *                  l = a
 *                  r = 1 (a dummy state, e.g. nothing)
 *                prv = y
 *           The old node y gets connectivity :
 *             l or r = n   
 *           The downstream node a gets a new parent:
 *                if (a != 1) a's prv = n   
 *
 *           Usually we're attaching a node, so a and b are the
 *           terminal dummy state 1, which does not remember its
 *           parents.
 *           
 *           For the special case of initializing the root node, use y==-1
 *           and whichway==TRACE_LEFT_CHILD. 
 *
 *           <mode> is the alignment mode, either TRMODE_J, TRMODE_L, TRMODE_R, or TRMODE_T.
 *           
 * Returns:  index of new node.
 */          
int
InsertTraceNodewithMode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state, char mode)
{
  int a;
  int n;

  n = tr->n;
	/* a==-1 unless we're inserting a node into an existing tree, which is rare */
  if (y >= 0)
    a = (whichway == TRACE_LEFT_CHILD ? tr->nxtl[y] : tr->nxtr[y]);
  else 
    a = -1;			/* special case of initializing the root. */

  if (tr->n == tr->nalloc) GrowParsetree(tr);
				/* information in new node */
  tr->emitl[n] = emitl;
  tr->emitr[n] = emitr;
  tr->state[n] = state;
  tr->mode[n]  = mode;
				/* connectivity of new node */
  tr->nxtl[n]  = a;
  tr->nxtr[n]  = -1;
  tr->prv[n]   = y;
				/* connectivity of parent   */
  if (y >= 0) {
    if (whichway == TRACE_LEFT_CHILD)  tr->nxtl[y] = n;
    else                               tr->nxtr[y] = n;
  }
				/* connectivity of child, 
				   if we're inserting instead of just adding  */
  if (a != -1)  tr->prv[a] = n;
				/* bump counter, return index of new node */
  tr->n++;
  return n;
}

/* Function: InsertTraceNode()
 *
 * Purpose:  Standard, non-mode-aware API
 *           Calls InsertTraceNodewithMode()
 *           with default mode value
 *
 * Returns:  index of new node
 */
int
InsertTraceNode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state)
{
   int n;

   n = InsertTraceNodewithMode(tr, y, whichway, emitl, emitr, state, TRMODE_J);

   return n;
}

/* Function: ParsetreeCount()
 * Date:     SRE, Mon Jul 31 19:19:08 2000 [St. Louis]
 *
 * Purpose:  Count a parsetree into a counts-based CM structure,
 *           in the course of estimating new CM probability parameters.
 *
 * Args:     cm   - CM to collect counts in
 *           tr   - the parse tree to collect from.
 *           dsq  - digitized sequence that we're counting symbols from
 *           wgt  - weight on this sequence (often just 1.0)
 *
 * Returns:  (void)
 */
void
ParsetreeCount(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, float wgt)
{
  int tidx;			/* counter through positions in the parsetree        */
  int v,z;			/* parent, child state index in CM                   */

		/* trivial preorder traverse, since we're already numbered that way */
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    if (v != cm->M && cm->sttype[v] != E_st && cm->sttype[v] != B_st) 
      {
	z = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */

	if (z == cm->M)                
	  cm->end[v] += wgt;
	else if (v == 0 && z - cm->cfirst[v] >= cm->cnum[v])
	  cm->begin[z] += wgt;
	else
	  cm->t[v][z - cm->cfirst[v]] += wgt; 

	if (cm->sttype[v] == MP_st) 
	  PairCount(cm->abc, cm->e[v], dsq[tr->emitl[tidx]], dsq[tr->emitr[tidx]], wgt);
	else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	  esl_abc_FCount(cm->abc, cm->e[v], dsq[tr->emitl[tidx]], wgt);
	else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	  esl_abc_FCount(cm->abc, cm->e[v], dsq[tr->emitr[tidx]], wgt);
      }
  }
}    
    
/* Function: ParsetreeScore()
 * Date:     SRE, Wed Aug  2 13:54:07 2000 [St. Louis]
 *
 * Purpose:  Calculate the score of a given parse tree for a sequence,
 *           given a CM that's prepared in log-odds form. Also calculate
 *           the contribution of structure to that score, by summing the 
 *           difference in MP emissions and the marginalized left and right
 *           scores. Also calculate the primary sequence score, the sum of
 *           all singlet emissions, plus marginalized base pair emissions.
 *           Also determine the first (ret_spos) and last (ret_epos) consensus
 *           columns occupied in the parsetree.
 *
 * Returns:  eslOK on success
 *           eslEINCOMPAT on contract violation.
 */
int 
ParsetreeScore(CM_t *cm, CMEmitMap_t *emap, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int do_null2, float *ret_sc, float *ret_struct_sc, float *ret_primary_sc, int *ret_spos, int *ret_epos)
{
  int status;                   /* Easel status code */
  int tidx;			/* counter through positions in the parsetree        */
  int v,y;			/* parent, child state index in CM                   */
  ESL_DSQ symi, symj;		/* symbol indices for emissions, 0..cm->abc->Kp-1    */
  float sc;			/* the log-odds score of the parse tree */
  char mode;
  float struct_sc;              /* contribution of the structure to the score */
  float primary_sc;             /* contribution of primary sequence emissions to the score */
  float lsc, rsc;
  int   nd;

  /* contract check */
  if(dsq    == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "ParsetreeScore(): dsq == NULL.");
  if(ret_sc == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "ParsetreeScore(): ret_sc == NULL.");
  if(emap   == NULL && (ret_spos != NULL || ret_epos != NULL)) ESL_FAIL(eslEINCOMPAT, errbuf, "ParsetreeScore(): emap == NULL but ret_spos and ret_epos != NULL.");

  int spos = cm->clen + 1;
  int epos = 0;

		/* trivial preorder traverse, since we're already numbered that way */
  sc = struct_sc = primary_sc = 0.;
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    mode = tr->mode[tidx];
    if (v == cm->M) continue;      	/* special case: v is EL, local alignment end */
    nd = cm->ndidx[v];
    if (cm->sttype[v] != E_st && cm->sttype[v] != B_st && tr->nxtl[tidx] != -1) /* no scores in B,E or if nxtl is -1 */
      {
	y = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */
	if(v == 0) { 
	  if(cm->flags & CMH_LOCAL_BEGIN)
	    sc += (tr->is_std) ? cm->beginsc[y] : cm->trbeginsc[y];
	  else if(tr->mode[tidx] == TRMODE_T && cm->sttype[y] == B_st)
	    sc += 0.; /* special case in *glocal* truncated alignment, 'local' begins into B states are free to allow for T alignment */                          
	  else
	    sc += cm->tsc[v][y - cm->cfirst[v]]; /* non-local transition out of ROOT_S */
	}
	else if (y == cm->M) /* CMH_LOCAL_END is presumably set, else this wouldn't happen */
	  sc += cm->endsc[v] + (cm->el_selfsc * (tr->emitr[tidx] - tr->emitl[tidx] + 1 - StateDelta(cm->sttype[v])));
	else 		/* y - cm->first[v] gives us the offset in the transition vector */
	  sc += cm->tsc[v][y - cm->cfirst[v]];
	
	if (cm->sttype[v] == MP_st) 
	  {
	    symi = dsq[tr->emitl[tidx]];
	    symj = dsq[tr->emitr[tidx]];
            if (mode == TRMODE_J)
              {
  	        if (symi < cm->abc->K && symj < cm->abc->K) { 
	          sc += cm->esc[v][(int) (symi*cm->abc->K+symj)];
		  struct_sc += cm->esc[v][(int) (symi*cm->abc->K+symj)];
		}
	        else { 
	          sc += DegeneratePairScore(cm->abc, cm->esc[v], symi, symj);
		  struct_sc += cm->esc[v][(int) (symi*cm->abc->K+symj)];
		}
		lsc = cm->lmesc[v][symi];
		rsc = cm->rmesc[v][symj];
		struct_sc -= lsc;  /* subtract left  marginalized score */
		struct_sc -= rsc; /* subtract right marginalized score */
		primary_sc += lsc;
		primary_sc += rsc;
	      }
            else if (mode == TRMODE_L)
              sc += cm->lmesc[v][symi];
            else if (mode == TRMODE_R)
              sc += cm->rmesc[v][symj];
	    if(emap != NULL) { 
	      spos = ESL_MIN(spos, emap->lpos[nd]);
	      epos = ESL_MAX(epos, emap->rpos[nd]);
	    }
	  } 
	else if ( (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) && (mode == TRMODE_J || mode == TRMODE_L) )
	  {
	    symi = dsq[tr->emitl[tidx]];
	    if (symi < cm->abc->K) lsc = cm->esc[v][(int) symi];
	    else                   lsc = esl_abc_FAvgScore(cm->abc, symi, cm->esc[v]);
	    sc += lsc;
	    primary_sc += lsc;
	    if(emap != NULL && cm->stid[v] == MATL_ML) { 
	      spos = ESL_MIN(spos, emap->lpos[nd]);
	      epos = ESL_MAX(epos, emap->lpos[nd]);
	    }
	  } 
	else if ( (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) && (mode == TRMODE_J || mode == TRMODE_R) )
	  {
	    symj = dsq[tr->emitr[tidx]];
	    if (symj < cm->abc->K) rsc = cm->esc[v][(int) symj];
	    else                   rsc = esl_abc_FAvgScore(cm->abc, symj, cm->esc[v]);
	    sc += rsc;
	    primary_sc += rsc;
	    if(emap != NULL && cm->stid[v] == MATR_MR) { 
	      spos = ESL_MIN(spos, emap->rpos[nd]);
	      epos = ESL_MAX(epos, emap->rpos[nd]);
	    }
	  }
      }
  }

  if(do_null2) { 
    float corr_sc;
    if((status = ParsetreeScoreCorrectionNull2(cm, errbuf, tr, dsq, 0, cm->null2_omega, &corr_sc)) != eslOK) return status;
    sc -= corr_sc;
    primary_sc -= corr_sc;
    /* don't subtract corr_sc from struct_sc, b/c we would have subtracted it from 
     * both the marginalized and non-marginalized MP scores, thus it cancels out for struct_sc 
     */
  }
  if(ret_sc != NULL)         *ret_sc        = sc;
  if(ret_struct_sc != NULL)  *ret_struct_sc = struct_sc;
  if(ret_primary_sc != NULL) *ret_primary_sc = primary_sc;

  if(spos == cm->clen+1) spos = -1;
  if(epos == 0)          epos = -1;

  if(ret_spos != NULL)  *ret_spos = spos;
  if(ret_epos != NULL)  *ret_epos = epos;
  return eslOK;
}




/* Function: PrintParsetree()
 * Date:     SRE, Fri Jul 28 12:47:06 2000 [St. Louis]
 *
 * Purpose:  Debugging: show a tabular representation of a
 *           parsetree structure.
 *           
 *           This just shows information in the
 *           parsetree structure itself. ParsetreeDump() 
 *           is more detailed, showing sequence information
 *           aligned to the tree. PrintParsetree() is
 *           called by cmbuild.c to print a master guide
 *           tree, which doesn't have an individual 
 *           sequence aligned to it.
 *
 * Args:     fp  - output stream (stdout?)
 *           tr  - the tree to show
 *
 * Returns:  void
 */
void
PrintParsetree(FILE *fp, Parsetree_t *tr)
{
  int x;

  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s %5s\n",
	  " idx ","emitl", "emitr", "state", " nxtl", " nxtr", " prv ", " mode");
  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s %5s\n",
	  "-----", "-----", "-----", "-----", "-----","-----", "-----", "-----");
  for (x = 0; x < tr->n; x++)
    fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5s\n",
	    x, tr->emitl[x], tr->emitr[x], tr->state[x], 
	    tr->nxtl[x], tr->nxtr[x], tr->prv[x], MarginalMode(tr->mode[x]));
  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s %5s\n",
	  "-----", "-----", "-----", "-----","-----","-----", "-----", "-----");

  fprintf(fp, "n      = %d\n", tr->n);
  fprintf(fp, "nalloc = %d\n", tr->nalloc);
  fprintf(fp, "block  = %d\n", tr->memblock);
}

/* Function: ParsetreeDump()
 * Date:     SRE, Fri Aug  4 10:43:20 2000 [St. Louis]
 *
 * Purpose:  Generate a detailed picture of a parsetree data structure,
 *           annotated with relevant information from the sequence
 *           and the model.
 *
 * Args:    fp    - FILE to write output to.
 *          tr    - parsetree to examine.
 *          cm    - model that was aligned to dsq to generate the parsetree
 *          dsq   - digitized sequence that was aligned to cm to generate the parsetree
 *          dmin  - minimum subseq length for each state; [0..v..M-1] NULL for non-banded output
 *          dmax  - maximum subseq length for each state; [0..v..M-1] NULL for non-banded output
 *
 * Returns:  (void)
 */
void
ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax)
{
  int   x;
  char  syml, symr;
  float tsc;
  float esc;
  int   v,y;
  char  mode;
  int   do_banded;
  int   L;

  /* Contract check */
  if(dmin == NULL && dmax != NULL)
    cm_Fail("In ParsetreeDump(), dmin is NULL, dmax is not.\n");
  if(dmin != NULL && dmax == NULL)
    cm_Fail("In ParsetreeDump(), dmax is NULL, dmin is not.\n");
  if(dsq == NULL)
    cm_Fail("In ParsetreeDump(), dsq is NULL");

  if(dmin != NULL && dmax != NULL) do_banded = TRUE;
  else                             do_banded = FALSE;

  if(do_banded)
    {
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	      " idx ", "emitl", "emitr", "state", " mode", " nxtl", " nxtr", " prv ", " tsc ", " esc ", 
	      " L   ", " dmin", " dmax");
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	      "-----", "------", "------", "-------", "-----", "-----","-----", "-----","-----", "-----",
	      "-----", "-----", "-----");
    }
  else
    {
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s\n",
	      " idx ","emitl", "emitr", "state", " mode", " nxtl", " nxtr", " prv ", " tsc ", " esc ");
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s\n",
	      "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----", "-----");
    }
  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      mode = tr->mode[x];

      /* Set syml, symr: one char representation of what we emit, or ' '.
       * Set esc:        emission score, or 0.
       * Only P, L, R states have emissions.
       */
      syml = symr = ' ';
      esc = 0.;
      if (cm->sttype[v] == MP_st) {
	if (mode == TRMODE_J || mode == TRMODE_L) syml = cm->abc->sym[dsq[tr->emitl[x]]]; 
	if (mode == TRMODE_J || mode == TRMODE_R) symr = cm->abc->sym[dsq[tr->emitr[x]]];
	if      (mode == TRMODE_J) esc = DegeneratePairScore(cm->abc, cm->esc[v], dsq[tr->emitl[x]], dsq[tr->emitr[x]]);
        else if (mode == TRMODE_L) esc = cm->lmesc[v][dsq[tr->emitl[x]]];
        else if (mode == TRMODE_R) esc = cm->rmesc[v][dsq[tr->emitr[x]]];
      } else if ( (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) && (mode == TRMODE_J || mode == TRMODE_L) ) {
	syml = cm->abc->sym[dsq[tr->emitl[x]]];
	esc  = esl_abc_FAvgScore(cm->abc, dsq[tr->emitl[x]], cm->esc[v]);
      } else if ( (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) && (mode == TRMODE_J || mode == TRMODE_R) ) {
	symr = cm->abc->sym[dsq[tr->emitr[x]]];
	esc  = esl_abc_FAvgScore(cm->abc, dsq[tr->emitr[x]], cm->esc[v]);
      }

      /* Set tsc: transition score, or 0.
       * B, E, and the special EL state (M, local end) have no transitions.
       */
      tsc = 0.;
      if (v != cm->M && cm->sttype[v] != B_st && cm->sttype[v] != E_st && tr->nxtl[x] != -1) {
	y = tr->state[tr->nxtl[x]];
	if(v == 0) { 
	  if(cm->flags & CMH_LOCAL_BEGIN)
	    tsc = (tr->is_std) ? cm->beginsc[y] : cm->trbeginsc[y];
	  else if(tr->mode[x] == TRMODE_T && cm->sttype[y] == B_st)
	    tsc = 0.; /* special case in *glocal* truncated alignment, 'local' begins into B states are free to allow for T alignment */                          
	  else
	    tsc = cm->tsc[v][y - cm->cfirst[v]]; /* non-local transition out of ROOT_S */
	}
	else if (y == cm->M) /* CMH_LOCAL_END is presumably set, else this wouldn't happen */
	  tsc = cm->endsc[v] + (cm->el_selfsc * (tr->emitr[x] - tr->emitl[x] + 1 - StateDelta(cm->sttype[v])));
	else 		/* y - cm->first[v] gives us the offset in the transition vector */
	  tsc = cm->tsc[v][y - cm->cfirst[v]];
      }

      /* Print the info line for this state
       */
      if(do_banded)
	{
	  L = tr->emitr[x]-tr->emitl[x]+1;
	  fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5s %5d %5d %5d %5.2f %5.2f %5d %5d %5d %2s\n",
		  x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], 
		  Statetype(cm->sttype[v]), MarginalMode(tr->mode[x]), 
		  tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc,
		  L, dmin[v], dmax[v],
		  (L >= dmin[v] && L <= dmax[v]) ? "" : "!!");
	}
      else
	{
	  fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5s %5d %5d %5d %5.2f %5.2f\n",
		  x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], 
		  Statetype(cm->sttype[v]), MarginalMode(tr->mode[x]),
		  tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc);
	}
    }
  if(do_banded)
    fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	    "-----", "------", "------", "-------", "-----","-----","-----", "-----","-----", "-----",
	    "-----", "-----", "-----");
  else
    fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s\n",
	    "-----", "------", "------", "-------", "-----","-----","-----", "-----","-----", "-----");
  fflush(fp);
} 


/* Function: ParsetreeCompare()
 * Date:     SRE, Sat Aug 12 22:05:38 2000 [Titusville]
 *
 * Purpose:  Compare two parse trees to each other, for debugging
 *           purposes. If they are not exactly alike, return 0.
 *           Else return 1.
 */
int
ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2)
{
  int x;

  if (t1->n != t2->n) return 0;
  for (x = 0; x < t1->n; x++) 
    {
      if (t1->emitl[x] != t2->emitl[x]) return 0;
      if (t1->emitr[x] != t2->emitr[x]) return 0;
      if (t1->state[x] != t2->state[x]) return 0;
      if (t1->mode[x]  != t2->mode[x])  return 0;
      if (t1->nxtl[x]  != t2->nxtl[x])  return 0;
      if (t1->nxtr[x]  != t2->nxtr[x])  return 0;
    }
  return 1;
}


/* Function: SummarizeMasterTrace()
 * Date:     SRE, Fri Jul 28 13:42:30 2000 [St. Louis]
 *
 * Purpose:  Debugging: count the nodes used in a master trace.
 *           Note that it takes advantage of the overloading of
 *           tr->state; in a master trace, this is a node type
 *           (e.g. MATP_nd), not a state index.
 *
 * Args:     fp - output file (e.g. stdout)
 *           tr - master trace to summarize
 *
 * Returns:  void
 */
void
SummarizeMasterTrace(FILE *fp, Parsetree_t *tr)
{
  int x;
  int count[NODETYPES];
  
  for (x = 0; x < NODETYPES; x++) count[x] = 0;
  for (x = 0; x < tr->n; x++)     count[tr->state[x]]++;

  fprintf(fp, "Summary report for the master trace:\n");
  fprintf(fp, "------------------------------------\n");
  fprintf(fp, "Total nodes:  %d\n", tr->n);
  fprintf(fp, "Bifurcations: %d\n", count[0]);
  fprintf(fp, "MATP:         %d\n", count[1]);
  fprintf(fp, "MATL:         %d\n", count[2]);
  fprintf(fp, "MATR:         %d\n", count[3]);
  fprintf(fp, "BEGL:         %d\n", count[4]);
  fprintf(fp, "BEGR:         %d\n", count[5]);
  fprintf(fp, "ROOT:         %d\n", count[6]);
  fprintf(fp, "END:          %d\n", count[7]);
}

/* Function: MasterTraceDisplay()
 * Date:     SRE, Mon Aug  7 10:05:16 2000 [St. Louis]
 *
 * Purpose:  prettified display of a master trace, for debugging
 *           and planning purposes. works by recursively calling
 *           mtd_visit_node().
 */
static void
mtd_visit_node(FILE *fp, Parsetree_t *mtr, CM_t *cm, int v, int depth)
{
  int y;
				/* find next start states in "binary tree" */
  for (y = v+1; y < mtr->n; y++)
    if (mtr->state[y] == END_nd || mtr->state[y] == BIF_nd) break;
				/* visit right */
  if (mtr->state[y] == BIF_nd)
    mtd_visit_node(fp, mtr, cm, mtr->nxtr[y], depth+1);
				/* deal with root */
  fprintf(fp, "%*s%d: %d[%d]: %d..%d, %d nt\n", depth*6, "", depth, v, cm->nodemap[v], mtr->emitl[v], mtr->emitr[v], mtr->emitr[v] - mtr->emitl[v] +1);
				/* visit left */
  if (mtr->state[y] == BIF_nd)
    mtd_visit_node(fp, mtr, cm, mtr->nxtl[y], depth+1);
}
void
MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm)
{
  mtd_visit_node(fp, mtr, cm, 0, 0);
}


/* Function : Parsetrees2Alignment()
 *
 * Purpose:   Creates a MSA from a set of parsetrees and a CM. If (do_matchonly)
 *            the MSA will only include consensus columns. If (do_resonly) only
 *            columns including at least 1 residue are included (by default, all
 *            consensus columns are included even if 0 seqs have a residue in them).
 * 
 * Args:     cm         - the CM the CP9 was built from, needed to get emitmap,
 *                        so we know where to put EL transitions
 *           errbuf     - for error messages.
 *           abc        - alphabet to use to create the return MSA
 *           sq         - sequences, must be digitized (we check for it)
 *           wgt        - weights for seqs, NULL for none
 *           tr         - array of tracebacks
 *           postcode   - posterior code string (NULL for none)
 *           nseq       - number of sequences
 *           insertfp   - file to print per-seq insert information to (NULL if none)
 *           elfp       - file to print per-seq EL insert information to (NULL if none)
 *           do_full    - TRUE to always include all match columns in alignment
 *           do_matchonly - TRUE to ONLY include match columns
 *           be_efficient - TRUE to free a parsetree as soon as we've created an aligned seq
 *           ret_msa    - MSA, alloc'ed/created here
 *
 * Returns:   eslOK on success, eslEMEM on memory error, eslEINVALID on contract violation.
 *            Also ret_msa is filled with a new MSA.
 *
 */
int
Parsetrees2Alignment(CM_t *cm, char *errbuf, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
		     Parsetree_t **tr, char **postcode, int nseq,
		     FILE *insertfp, FILE *elfp, int do_full, int do_matchonly, int be_efficient, ESL_MSA **ret_msa)
{
  int          status;          /* easel status flag */
  ESL_MSA     *msa   = NULL;    /* multiple sequence alignment */
  CMEmitMap_t *emap  = NULL;    /* consensus emit map for the CM */
  int          i;               /* counter over traces */
  int          v, nd;           /* state, node indices */
  int          cpos;            /* counter over consensus positions (0)1..clen */
  int         *matuse= NULL;    /* TRUE if we need a cpos in mult alignment */
  int         *iluse = NULL;    /* # of IL insertions after a cpos for 1 trace */
  int         *eluse = NULL;    /* # of EL insertions after a cpos for 1 trace */
  int         *iruse = NULL;    /* # of IR insertions after a cpos for 1 trace */
  int         *maxil = NULL;    /* max # of IL insertions after a cpos */
  int         *maxel = NULL;    /* max # of EL insertions after a cpos */
  int         *maxir = NULL;    /* max # of IR insertions after a cpos */
  int	      *matmap= NULL;    /* apos corresponding to a cpos */
  int         *ilmap = NULL;    /* first apos for an IL following a cpos */
  int         *elmap = NULL;    /* first apos for an EL following a cpos */
  int         *irmap = NULL;    /* first apos for an IR following a cpos */
  int          alen;	        /* length of msa in columns */
  int          apos;	        /* position in an aligned sequence in MSA */
  int          rpos;	        /* position in an unaligned sequence in dsq */
  int         *ifirst = NULL;   /* first uapos (unaligned position) for an insert following a cpos in cur seq */
  int         *elfirst = NULL;  /* first uapos (unaligned position) for an EL following a cpos in cur seq */
  int          tpos;            /* position in a parsetree */
  int          el_len;	        /* length of an EL insertion in residues */
  CMConsensus_t *con = NULL;    /* consensus information for the CM */
  int          prvnd;	        /* keeps track of previous node for EL */
  int          nins;            /* insert counter used for splitting inserts */
  int          do_post;         /* TRUE if postcode != NULL (we should write at least 1 sequence's posteriors) */
  int          do_cur_post;     /* TRUE if we're writing posteriors for the current sequence inside a loop */
  char        *tmp_aseq = NULL; /* will hold aligned sequence for current sequence */
  char        *tmp_apc  = NULL; /* will hold aligned postcode characters for current sequence */
  /* s_cposA and e_cposA are only alloc'ed and filled if insertfp || elfp != NULL, b/c they're only purpose is to be written to those files */
  int         *s_cposA = NULL;  /* [0..nseq-1] the first consensus position filled by a nongap for seq i */
  int         *e_cposA = NULL;  /* [0..nseq-1] the final consensus position filled by a nongap for seq i */
  int          s_cpos;          /* first consensus position filled by a nongap for current seq */
  int          e_cpos;          /* final consensus position filled by a nongap for current seq */

  /* Contract check. We allow the caller to specify the alphabet they want the 
   * resulting MSA in, but it has to make sense (see next few lines). */
  if(cm->abc->type == eslRNA) {  
      if(abc->type != eslRNA && abc->type != eslDNA)
	ESL_XFAIL(eslEINVAL, errbuf, "Error in Parsetrees2Alignment(), cm alphabet is RNA, but requested output alphabet is neither DNA nor RNA.");
  }
  else if(cm->abc->K != abc->K) {
    ESL_XFAIL(eslEINVAL, errbuf, "Error in Parsetrees2Alignment(), cm alphabet size is %d, but requested output alphabet size is %d.", cm->abc->K, abc->K);
  }

  /* Determine if we are doing posteriors. We should just be able to 
   * check if postcode array is NULL, but sometimes it is non-NULL
   * yet all postcode[] els are NULL (e.g. for cmalign --mpi --no-prob).
   */
  do_post = FALSE;
  if(postcode != NULL) { 
    for(i = 0; i < nseq; i++) {
      if(postcode[i] != NULL) { 
	do_post = TRUE;
	break;
      }
    }
  }

  emap = CreateEmitMap(cm);

  ESL_ALLOC(matuse,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(iluse,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(eluse,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(iruse,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxil,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxel,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxir,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(matmap,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(ilmap,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(elmap,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(irmap,   sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(ifirst,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(elfirst, sizeof(int)*(emap->clen+1));   

  if(insertfp != NULL || elfp != NULL) { 
    ESL_ALLOC(s_cposA,  sizeof(int)*nseq);
    ESL_ALLOC(e_cposA,  sizeof(int)*nseq);
    esl_vec_ISet(s_cposA, nseq, 0); /* all nseq values will be overwritten, actually this initialization is unnec */
    esl_vec_ISet(e_cposA, nseq, 0); /* all nseq values will be overwritten, actually this initialization is unnec */
  }

  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      if(!do_full || cpos == 0)
	matuse[cpos] = 0;
      else
	matuse[cpos] = 1;
      maxil[cpos] = maxel[cpos] = maxir[cpos] = 0;
      ilmap[cpos] = elmap[cpos] = irmap[cpos] = 0;
    }

  /* Look at all the traces; find maximum length of
   * insert needed at each of the clen+1 possible gap
   * points. (There are three types of insert, IL/EL/IR.)
   * Also find whether we don't need some of the match
   * (consensus) columns.
   */
  for (i = 0; i < nseq; i++) 
    {
      for (cpos = 0; cpos <= emap->clen; cpos++) 
	iluse[cpos] = eluse[cpos] = iruse[cpos] = 0;

      for (tpos = 0; tpos < tr[i]->n; tpos++)
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];
	  
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    matuse[emap->lpos[nd]] = 1;
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case ML_st:
	    matuse[emap->lpos[nd]] = 1;
	    break;
	  case MR_st:
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case IL_st:
	    iluse[emap->lpos[nd]]++;
	    break;
	  case IR_st:		
            /* remember, convention on rpos is that IR precedes this
             * cpos. Make it after the previous cpos, hence the -1. 
	     */
	    iruse[emap->rpos[nd]-1]++;
	    break;
	  case EL_st:
	    el_len = tr[i]->emitr[tpos] - tr[i]->emitl[tpos] + 1;
	    eluse[emap->epos[nd]] = el_len;
            /* not possible to have >1 EL in same place; could assert this */
	    break;
	  }

	  prvnd = nd;
	} /* end looking at trace i */

      for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
	  if (iluse[cpos] > maxil[cpos]) maxil[cpos] = iluse[cpos];
	  if (eluse[cpos] > maxel[cpos]) maxel[cpos] = eluse[cpos];
	  if (iruse[cpos] > maxir[cpos]) maxir[cpos] = iruse[cpos];
	}
    } /* end calculating lengths used by all traces */

  /* Now we can calculate the total length of the multiple alignment, alen;
   * and the maps ilmap, elmap, and irmap that turn a cpos into an apos
   * in the multiple alignment: e.g. for an IL that follows consensus position
   * cpos, put it at or after apos = ilmap[cpos] in aseq[][].
   * IR's are filled in backwards (3'->5') and rightflushed.
   */
  alen = 0;
  for (cpos = 0; cpos <= emap->clen; cpos++)
    {
      if (matuse[cpos]) {
	matmap[cpos] = alen; 
	alen++;
      } else 
	matmap[cpos] = -1;

      ilmap[cpos] = alen; 
      if(! do_matchonly) alen += maxil[cpos];
      elmap[cpos] = alen; 
      if(! do_matchonly) alen += maxel[cpos];
      if(! do_matchonly) alen += maxir[cpos]; 
      irmap[cpos] = alen-1; 
      /* note: if do_matchonly, no inserts are printed, ilmap, elmap, irmap are irrelevant */
    }

  /* We're getting closer.
   * Now we know the size of the MSA, but we only allocate a chunk of seqs
   * at a time, so that if <be_efficient>, we can free each parsetree as 
   * we create an aligned seq, to save memory.
   */
  msa = esl_msa_Create(nseq, -1);
  if(msa == NULL) goto ERROR;
  msa->nseq = nseq;
  msa->alen = alen;
  msa->abc  = NULL;
  if(do_post) { 
    ESL_ALLOC(msa->pp, sizeof(char *) * msa->nseq);
  }

  /* we reuse aseq, copying it to the msa after it's completed for each seq */
  ESL_ALLOC(tmp_aseq, sizeof(char) * (msa->alen+1));
  if(do_post) { /* these will be reused for each sequence (the aligned postcode arrays are all the same length) */
    ESL_ALLOC(tmp_apc, sizeof(char) * (msa->alen+1));
  }

  for (i = 0; i < nseq; i++)
    {
      s_cpos = emap->clen+1; /* an ESL_MIN with any cpos <= clen will replace this */
      e_cpos = 0;            /* an ESL_MAX with any cpos  > 0 will replace this */

      /* Contract check */
      if(sq[i]->dsq == NULL) ESL_XFAIL(eslEINVAL, errbuf, "Error in Parsetrees2Alignment(), sq %d is not digitized.\n", i);
      do_cur_post = FALSE;
      if(do_post) { 
	do_cur_post = (postcode[i] != NULL) ? TRUE : FALSE;
      }

      /* Initialize the aseq with all pads '.' (in insert cols) 
       * and deletes '-' (in match cols).
       */
      if(! do_matchonly) { /* if do_matchonly, we'll have no insert columns */
	for (apos = 0; apos < alen; apos++) { 
	  tmp_aseq[apos] = '.';
	}
	if(do_cur_post) { 
	  for (apos = 0; apos < alen; apos++) { 
	    tmp_apc[apos] = '.';
	  }
	}
      }
      for (cpos = 0; cpos <= emap->clen; cpos++) {
	if (matmap[cpos] != -1) { 
	  tmp_aseq[matmap[cpos]] = '-';
	}
      }
      if(do_cur_post) { 
	for (cpos = 0; cpos <= emap->clen; cpos++) {
	  if (matmap[cpos] != -1) { 
	    tmp_apc[matmap[cpos]] = '.';
	  }
	}
      }
      tmp_aseq[alen] = '\0';
      if(do_cur_post) tmp_apc[alen] = '\0';

      /* Traverse this guy's trace, and place all his
       * emitted residues and posteriors.
       */
      for (cpos = 0; cpos <= emap->clen; cpos++) { 
	iluse[cpos] = iruse[cpos] = eluse[cpos] = 0;
	ifirst[cpos] = elfirst[cpos] = -1;
      }

      for (tpos = 0; tpos < tr[i]->n; tpos++) 
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];

	  switch (cm->sttype[v]) {
	  case MP_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    tmp_aseq[apos] = abc->sym[sq[i]->dsq[rpos]];
	    if(do_cur_post) { 
	      tmp_apc[apos] = postcode[i][rpos-1];
	    }
	    s_cpos = ESL_MIN(s_cpos, cpos);
	    e_cpos = ESL_MAX(e_cpos, cpos);

	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    tmp_aseq[apos] = abc->sym[sq[i]->dsq[rpos]];
	    if(do_cur_post) { 
	      tmp_apc[apos] = postcode[i][rpos-1];
	    }
	    s_cpos = ESL_MIN(s_cpos, cpos);
	    e_cpos = ESL_MAX(e_cpos, cpos);
	    break;
	    
	  case ML_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    tmp_aseq[apos] = abc->sym[sq[i]->dsq[rpos]];
	    if(do_cur_post) { 
	      tmp_apc[apos] = postcode[i][rpos-1];
	    }
	    s_cpos = ESL_MIN(s_cpos, cpos);
	    e_cpos = ESL_MAX(e_cpos, cpos);
	    break;

	  case MR_st:
	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    tmp_aseq[apos] = abc->sym[sq[i]->dsq[rpos]];
	    if(do_cur_post) { 
	      tmp_apc[apos] = postcode[i][rpos-1];
	    }
	    s_cpos = ESL_MIN(s_cpos, cpos);
	    e_cpos = ESL_MAX(e_cpos, cpos);
	    break;

	  case IL_st:
	    cpos = emap->lpos[nd];
	    apos = ilmap[cpos] + iluse[cpos];
	    rpos = tr[i]->emitl[tpos];
	    if(iluse[cpos] == 0) ifirst[cpos] = rpos; /* only update ifirst if this is the first insert for this IL */
	    iluse[cpos]++;
	    if(do_matchonly) break; /* we don't break until this point in case we're writing to insertfp, in which case we need ifirst[] and iluse[] */
	    tmp_aseq[apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
	    if(do_cur_post) { 
	      tmp_apc[apos] = postcode[i][rpos-1];
	    }
	    break;

	  case EL_st: 
            /* we can assert eluse[cpos] always == 0 when we enter,
	     * because we can only have one EL insertion event per 
             * cpos. If we ever decide to regularize (split) insertions,
             * though, we'll want to calculate eluse in the rpos loop.
             */
	    cpos = emap->epos[nd]; 
	    apos = elmap[cpos]; 
	    eluse[cpos] = tr[i]->emitr[tpos] - tr[i]->emitl[tpos] + 1;
	    elfirst[cpos] = tr[i]->emitl[tpos];
	    if(do_matchonly) break; /* we don't break until this point in case we're writing to elfp, in which case we need elfirst[] and eluse[] */
	    for (rpos = tr[i]->emitl[tpos]; rpos <= tr[i]->emitr[tpos]; rpos++)
	      {
		tmp_aseq[apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
		if(do_cur_post) { 
		  tmp_apc[apos] = postcode[i][rpos-1];
		}
		apos++;
	      }
	    break;

	  case IR_st: 
	    cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
	    apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
	    rpos = tr[i]->emitr[tpos];
	    ifirst[cpos] = rpos; /* by overwriting each time we will end up with min rpos used by this IR */
	    iruse[cpos]++;
	    if(do_matchonly) break; /* we don't break until this point in case we're writing to elfp, in which case we need ifirst[] and iruse[] */
	    tmp_aseq[apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
	    if(do_cur_post) { 
	      tmp_apc[apos] = postcode[i][rpos-1];
	    }
	    break;

	  case D_st:
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATL_D) 
	      {
		cpos = emap->lpos[nd];
		if (matuse[cpos]) tmp_aseq[matmap[cpos]] = '-';
	      }
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATR_D) 
	      {
		cpos = emap->rpos[nd];
		if (matuse[cpos]) tmp_aseq[matmap[cpos]] = '-';
	      }
	    break;

	  } /* end of the switch statement */


	  prvnd = nd;
	  /* we only set s_cposA[i] and e_cposA[i] if we'll output them to an insertfp or elfp */
	  if(insertfp != NULL || elfp != NULL) { 
	    s_cposA[i] = (s_cpos == (emap->clen+1)) ? -1 : s_cpos;
	    e_cposA[i] = (e_cpos == 0)              ? -1 : e_cpos;
	  }
	} /* end traversal over trace i. */

      /* copy tmp_aseq to the msa */
      if((status = esl_strdup(tmp_aseq, msa->alen, &(msa->aseq[i]))) != eslOK) goto ERROR;
      /* if be_efficient, free tr[i], we're done with it */
      if(be_efficient) { 
	FreeParsetree(tr[i]);
	tr[i] = NULL; 
      }

      /* add tmp_apc posterior probabilities to msa, if nec */
      if(do_cur_post) { 
	if((status = esl_strdup(tmp_apc, msa->alen, &(msa->pp[i]))) != eslOK) goto ERROR;
      }
      else if (do_post) { msa->pp[i] = NULL; }

      /* rejustify inserts and posteriors (if nec) */
      if(! do_matchonly) { 
	/* IL/EL Insertions are currently flush-left and IR insertions are currently flush-right.
	 * This is pre-1.0 Infernal behavior. If(cm->align_opts & CM_ALIGN_FLUSHINSERTS) we leave them all alone,
	 * otherwise we regularize (split) the internal inserts, we flush the 5' inserts right and the 3'
	 * inserts left (note: pre 1.0 behavior does the opposite, flushes 5' left (assuming they're ROOT_ILs)
	 * and flushes 3' right (assuming they're ROOT_IRs).
	 *
	 * We have to be careful about EL's. We don't want to group IL/IR's and EL's together and then split them
	 * because we need to annotate IL/IR's as '.'s in the consensus structure and EL's as '~'. So we split
	 * each group separately. There should only be either IL or IR's at any position (b/c assuming we've
	 * detached the CM grammar ambiguity (which is default in cmbuild)). But we don't count on it here.
	 */
	if(! (cm->align_opts & CM_ALIGN_FLUSHINSERTS)) /* default behavior, split insert in half */
	  {
	    /* Deal with inserts before first consensus position, ILs, then ELs, then IRs
	     * IL's are flush left, we want flush right */
	    rightjustify(abc, msa->aseq[i], maxil[0]);
	    if(do_cur_post) rightjustify(abc, msa->pp[i], maxil[0]);

	    /* EL's are flush left, we want flush right I think these are impossible, but just in case... */
	    rightjustify(abc, msa->aseq[i]+maxil[0], maxel[0]);
	    if(do_cur_post) rightjustify(abc, msa->pp[i]+maxil[0], maxel[0]);
	    /* IR's are flush right, we want flush right, do nothing */
	    
	    /* split all internal insertions */
	    for (cpos = 1; cpos < emap->clen; cpos++) 
	      {
		if(maxil[cpos] > 1) /* we're flush LEFT, want to split */
		  {
		    apos = matmap[cpos]+1;
		    for (nins = 0; islower((int) (msa->aseq[i][apos])); apos++)
		      nins++;
		    nins /= 2;		/* split the insertion in half */
		    rightjustify(abc, msa->aseq[i]+matmap[cpos]+1+nins, maxil[cpos]-nins);
		    if(do_cur_post) rightjustify(abc, msa->pp[i]+matmap[cpos]+1+nins, maxil[cpos]-nins);
		  }
		if(maxel[cpos] > 1) /* we're flush LEFT, want to split */
		  {
		    apos = matmap[cpos]+1 + maxil[cpos];
		    for (nins = 0; islower((int) (msa->aseq[i][apos])); apos++)
		      nins++;
		    nins /= 2;		/* split the insertion in half */
		    rightjustify(abc, msa->aseq[i]+matmap[cpos]+1+maxil[cpos]+nins, maxel[cpos]-nins);
		    if(do_cur_post) rightjustify(abc, msa->pp[i]+matmap[cpos]+1+maxil[cpos]+nins, maxel[cpos]-nins);
		  }
		if(maxir[cpos] > 1) /* we're flush RIGHT, want to split */
		  {
		    apos = matmap[cpos+1]-1;
		    for (nins = 0; islower((int) (msa->aseq[i][apos])); apos--)
		      nins++;
		    nins ++; nins /= 2;		/* split the insertion in half (++ makes it same behavior as IL/EL */
		    leftjustify(abc, msa->aseq[i]+matmap[cpos]+1 + maxil[cpos] + maxel[cpos], maxir[cpos]-nins);
		    if(do_cur_post) leftjustify(abc, msa->pp[i]+matmap[cpos]+1 + maxil[cpos] + maxel[cpos], maxir[cpos]-nins);
		  }
	      }
	    /* Deal with inserts after final consensus position, IL's then EL's, then IR's
	     * IL's are flush left, we want flush left, do nothing 
	     * EL's are flush left, we want flush left, do nothing 
	     * IR's are flush right, we want flush left */
	    leftjustify(abc, msa->aseq[i]+matmap[emap->clen]+1 + maxil[emap->clen] + maxel[emap->clen], maxir[emap->clen]);
	    if(do_cur_post) leftjustify(abc, msa->pp[i]+matmap[emap->clen]+1 + maxil[emap->clen] + maxel[emap->clen], maxir[emap->clen]);
	  }
      }
    
      /* output insert and/or EL info to the insertfp and elfp output files, if nec */
      if(insertfp != NULL || elfp != NULL) { 
	if(insertfp != NULL) { fprintf(insertfp, "%s %" PRId64 " %d %d", sq[i]->name, sq[i]->n, s_cposA[i], e_cposA[i]); }
	if(elfp != NULL)     { fprintf(elfp,     "%s %" PRId64 " %d %d", sq[i]->name, sq[i]->n, s_cposA[i], e_cposA[i]); }
	for (cpos = 0; cpos <= emap->clen; cpos++) 
	  {
	    if((insertfp != NULL) && ((iluse[cpos] + iruse[cpos]) > 0)) { 
	      fprintf(insertfp, "  %d %d %d", cpos, ifirst[cpos], (iluse[cpos] + iruse[cpos])); /* note cpos+1 puts cpos from 1..clen, ifirst[] is already 1..sq->n */
	      /* Note: only 1 of iluse[cpos] or iruse[cpos] should be != 0 (I think) */
	    }
	    if((elfp != NULL) && (eluse[cpos] > 0)) { 
	      fprintf(elfp, "  %d %d %d", cpos, elfirst[cpos], eluse[cpos]); /* note cpos+1 puts cpos from 1..clen, ifirst[] is already 1..sq->n */
	    }
	  }
	if(insertfp != NULL) { fprintf(insertfp, "\n"); }
	if(elfp != NULL)     { fprintf(elfp,     "\n"); }
      }
    } /* end loop over all parsetrees */
  if(be_efficient) { 
    free(tr);
    *(&tr) = NULL;
  }

  /* Gee, wasn't that easy?
   * Add the rest of the ("optional") information to the MSA.
   */
  CreateCMConsensus(cm, abc, 3.0, 1.0, &con);

  /* "author" info */
  ESL_ALLOC(msa->au, sizeof(char) * (strlen(INFERNAL_VERSION)+10));
  sprintf(msa->au, "Infernal %s", INFERNAL_VERSION);

  /* per-seq info */
  for (i = 0; i < nseq; i++)
    {
      esl_msa_SetSeqName(msa, i, sq[i]->name, -1);
      if (sq[i]->acc[0]  != '\0') esl_msa_SetSeqAccession  (msa, i, sq[i]->acc,  -1);
      if (sq[i]->desc[0] != '\0') esl_msa_SetSeqDescription(msa, i, sq[i]->desc, -1);
      if (msa->sqlen != NULL) msa->sqlen[i] = sq[i]->n;
      if (wgt == NULL) msa->wgt[i] = 1.0;
      else             msa->wgt[i] = wgt[i];
      /* TODO: individual SS annotations */
    }

  /* Construct the secondary structure consensus line, msa->ss_cons:
   *       IL, IR are annotated as .
   *       EL is annotated as ~
   *       and match columns use the structure code.
   * Also the primary sequence consensus/reference coordinate system line,
   * msa->rf.
   */
  ESL_ALLOC(msa->ss_cons, (sizeof(char) * (alen+1)));
  ESL_ALLOC(msa->rf,      (sizeof(char) * (alen+1)));
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      if (matuse[cpos]) 
	{ /* CMConsensus is off-by-one right now, 0..clen-1 relative to cpos's 1..clen */

	  /* bug i1, xref STL7 p.12. Before annotating something as a base pair,
	   * make sure the paired column is also present.
	   */
	  if (con->ct[cpos-1] != -1 && matuse[con->ct[cpos-1]+1] == 0) {
	    msa->ss_cons[matmap[cpos]] = '.';
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  } else {
	    msa->ss_cons[matmap[cpos]] = con->cstr[cpos-1];	
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  }
	}
      if ((maxil[cpos] > 0) && (! do_matchonly)) 
	for (apos = ilmap[cpos]; apos < ilmap[cpos] + maxil[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
      if ((maxel[cpos] > 0) && (! do_matchonly))
	for (apos = elmap[cpos]; apos < elmap[cpos] + maxel[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '~';
	    msa->rf[apos] = '~';
	  }
      if ((maxir[cpos] > 0) && (! do_matchonly)) /* remember to write backwards */
	for (apos = irmap[cpos]; apos > irmap[cpos] - maxir[cpos]; apos--)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
    }
  msa->ss_cons[alen] = '\0';
  msa->rf[alen] = '\0';

  if(tmp_aseq != NULL) free(tmp_aseq);
  if(tmp_apc  != NULL) free(tmp_apc);
  if(s_cposA  != NULL) free(s_cposA);
  if(e_cposA  != NULL) free(e_cposA);
  FreeCMConsensus(con);
  FreeEmitMap(emap);
  free(matuse);
  free(iluse);
  free(eluse);
  free(iruse);
  free(maxil);
  free(maxel);
  free(maxir);
  free(matmap);
  free(ilmap);
  free(elmap);
  free(irmap);
  free(ifirst);
  free(elfirst);
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if(tmp_apc != NULL) free(tmp_apc);
  if(con   != NULL)  FreeCMConsensus(con);
  if(emap  != NULL)  FreeEmitMap(emap);
  if(matuse!= NULL)  free(matuse);
  if(iluse != NULL)  free(iluse);
  if(eluse != NULL)  free(eluse);
  if(iruse != NULL)  free(iruse);
  if(maxil != NULL)  free(maxil);
  if(maxel != NULL)  free(maxel);
  if(maxir != NULL)  free(maxir);
  if(matmap!= NULL)  free(matmap);
  if(ilmap != NULL)  free(ilmap);
  if(elmap != NULL)  free(elmap);
  if(irmap != NULL)  free(irmap);
  if(msa   != NULL)  esl_msa_Destroy(msa);
  return status;
}

/* Function: ParsetreeScore_Global2Local()
 * Date:     EPN, Wed May 23 09:57:38 2007
 *
 * Purpose:  Given a parsetree of dsq that corresponds to a globally
 *           configured CM, return the highest scoring local parsetree 
 *           of dsq or a subsequence of dsq (due to local begins)
 *           that is consistent with it. The hope is that this
 *           score will *be close* to the optimal local parse of dsq so we
 *           can calculate CP9 filter thresholds without the need to 
 *           search for the optimal local parse. 
 * 
 *           All residues 1..L must exist in the local parse emitted 
 *           from the same states they were emitted in the global
 *           parse unless (1) the local parse contains a local begin into
 *           state v at parstree node t, where tr->emitl[t] > 1 and/or
 *           tr->emitr[t] < L, (2) residues were emitted from an EL
 *           state because it was higher scoring than the subtree
 *           of the global parse.
 *           
 */
float
ParsetreeScore_Global2Local(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, int print_flag)
{
  int   status;
  int tidx;			/* counter through positions in the parsetree        */
  int v,y;			/* parent, child state index in CM                   */
  ESL_DSQ symi, symj;		/* symbol indices for emissions, 0..Kp-1             */
  char mode;
  int    tp;                    /* trace index offset, for v's with > tidx (IL or IR)*/
  float *tr_esc;                /* [0..tr->n-1] score of emissions from each trace node */
  float *tr_tsc;                /* [0..tr->n-1] score of transitions from each trace node */
  int   *v2n_map;               /* [0..cm->M-1], the trace node each state v corresponds to 
				 * -1 if none */
  int   *v2n_ct;                /* [0..cm->M-1], # of trace nodes state v corresponds to */
  float *lsc;                   /* [0..tr->n-1], the score of the best local parse *
				 * rooted at v = v2n_map[tidx] for trace node tidx *
				 * -1 if none */
  float max_local_sc;           /* the best local parse score consistent with tr */
  float below_me_sc;            /* score of tr-consistent best local parse under v */
  float tmp_endsc;              /* score of jumping out of v to EL */
  /* Contract check, CM must be LOCALLY configured, (could config to global, but
   * we assume we'll be calling this function serially for many parses and don't
   * want to need to switch CM back and forth from local/global) */
  if((!(cm->flags & CMH_LOCAL_BEGIN)) || (!(cm->flags & CMH_LOCAL_END)))
    cm_Fail("ERROR in ParsetreeScore_Global2Local() CM is not in local mode.\n");
  if(dsq == NULL)
    cm_Fail("ERROR in ParsetreeScore_Global2Local(), dsq is NULL.\n");

  /* Allocate and initialize */
  ESL_ALLOC(v2n_map, sizeof(int)   * cm->M); 
  ESL_ALLOC(v2n_ct,  sizeof(int)   * cm->M); 
  ESL_ALLOC(lsc,     sizeof(float) * tr->n);
  ESL_ALLOC(tr_esc,  sizeof(float) * tr->n); 
  ESL_ALLOC(tr_tsc,  sizeof(float) * tr->n); 
  esl_vec_ISet(v2n_map, cm->M, -1);
  esl_vec_ISet(v2n_ct,  cm->M, 0);
  esl_vec_FSet(lsc, tr->n, 0.);
  esl_vec_FSet(tr_tsc, tr->n, 0.);
  esl_vec_FSet(tr_esc, tr->n, 0.);

  /* Determine the score that each trace node contributes to the overall parsetree score */

  for (tidx = 0; tidx < tr->n; tidx++) 
    {
      v = tr->state[tidx];        	/* index of parent state in CM */
      v2n_map[v] = tidx;
      v2n_ct[v]++; /* insert states could be visited > once */
      mode = tr->mode[tidx];
      if (v == cm->M) 
	cm_Fail("ERROR in ParsetreeScore_Global2Local(), EL in parse, but it should be global!\n");
      if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) /* no scores in B,E */
	{
	  y = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */

	  if (y == cm->M) 
	    cm_Fail("ERROR in ParsetreeScore_Global2Local(), EL in parse, but it should be global!\n");
	  if (v == 0 && y > cm->cnum[0])
	    cm_Fail("ERROR in ParsetreeScore_Global2Local(), we did a local begin in the parse, but it should be global!\n");
	  /* for v == 0, we don't care that transition score has changed from global CM that
	   * was used to generate the parsetree, because the transition from root is not
	   * considered when we look for best local parse below. */

	  /* y - cm->first[v] gives us the offset in the transition vector */
	  tr_tsc[tidx] = cm->tsc[v][y - cm->cfirst[v]];
	
	  if (cm->sttype[v] == MP_st) 
	    {
	      symi = dsq[tr->emitl[tidx]];
	      symj = dsq[tr->emitr[tidx]];
	      if (mode == TRMODE_J)
		{
		  if (symi < cm->abc->K && symj < cm->abc->K)
		    tr_esc[tidx] = cm->esc[v][(int) (symi*cm->abc->K+symj)];
		  else
		    tr_esc[tidx] = DegeneratePairScore(cm->abc, cm->esc[v], symi, symj);
		}
	      else if (mode == TRMODE_L)
		tr_esc[tidx] = cm->lmesc[v][symi];
	      else if (mode == TRMODE_R)
		tr_esc[tidx] = cm->rmesc[v][symj];
	    } 
	  else if ( (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) && (mode == TRMODE_J || mode == TRMODE_L) )
	    {
	      symi = dsq[tr->emitl[tidx]];
	      if (symi < cm->abc->K) tr_esc[tidx] = cm->esc[v][(int) symi];
	      else                   tr_esc[tidx] = esl_abc_FAvgScore(cm->abc, symi, cm->esc[v]);
	    } 
	  else if ( (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) && (mode == TRMODE_J || mode == TRMODE_L) )
	    {
	      symj = dsq[tr->emitr[tidx]];
	      if (symj < cm->abc->K) tr_esc[tidx] = cm->esc[v][(int) symj];
	      else                   tr_esc[tidx] = esl_abc_FAvgScore(cm->abc, symj, cm->esc[v]);
	    }
	}
    }

  /* Now traverse CM from inside-out, for each v in the parse, 
   * keep track of the best local CM score of the parse rooted 
   * at v, with the possibility of local ends. Keep track 
   * of maximum score considering all possible local begins */
  max_local_sc = IMPOSSIBLE;
  for(v = cm->M-1; v > 0; v--)
    {
      for(tp = 0; tp < v2n_ct[v]; tp++) 
	{
	  tidx = v2n_map[v] - tp; 
	  if(print_flag) 
	    printf("loop start tidx: %d esc: %f tsc: %f\n", tidx, tr_esc[tidx], tr_tsc[tidx]);

	  if(cm->sttype[v] == B_st)
	    {
	      below_me_sc = lsc[tr->nxtl[tidx]] + lsc[tr->nxtr[tidx]];
	      if(print_flag) 
		{
		  printf("B state L %d: %f R %d: %f\n", tr->nxtl[tidx], lsc[tr->nxtl[tidx]], tr->nxtr[tidx], lsc[tr->nxtr[tidx]]); 
		}
	    }
	  else if(cm->sttype[v] == E_st)
	    below_me_sc = 0.;
	  else
	    {
	      below_me_sc = lsc[tr->nxtl[tidx]];
	      if(print_flag) printf("non B: below_me_sc %d: %f\n", tr->nxtl[tidx], lsc[tr->nxtl[tidx]]);
	    }
	  /* Check if we could've jumped to an EL instead of traversing the 
	   * subparse rooted here at v, would it have been worth it? */
	  tmp_endsc = cm->endsc[v] + /* score of transition to EL */
	    cm->el_selfsc *  /* score of emitting 1 residue */
	    ((tr->emitr[tidx] - StateRightDelta(cm->sttype[v])) - 
	     (tr->emitl[tidx] + StateLeftDelta(cm->sttype[v])) + 1); /* number of residues EL must emit */
	  if(below_me_sc < (tmp_endsc - tr_tsc[tidx])) /* careful to consider sc of transition out of v */
	    {
	      below_me_sc = tmp_endsc - tr_tsc[tidx]; /* we'll add tr_tsc[idx] back in next */
	      if(print_flag) 
		{
		  printf("\nTOOK LOCAL END!\n");
		  printf("tmp_endsc: %f cm->endsc: %f + %d emits\n", tmp_endsc, cm->endsc[v], ((tr->emitr[tidx] - StateRightDelta(cm->sttype[v])) - (tr->emitl[tidx] + StateLeftDelta(cm->sttype[v])) + 1)); 
		}
	    }
	  lsc[tidx] = tr_esc[tidx] + tr_tsc[tidx] + below_me_sc; /* note we add in tsc even if local end taken */
	  
	  if(print_flag) 
	    {
	      printf("tidx: %d\nv: %d\nlsc[tidx]: %f\nbegin_sc: %f\nmax_local_sc: %f\ntmp_endsc: %f\n\n", tidx, v, lsc[tidx], cm->beginsc[v], max_local_sc, tmp_endsc);
	      printf("tr_esc[tidx]: %f\ntr_tsc[tidx]: %f\nbelow_me_sc: %f\n", tr_esc[tidx], tr_tsc[tidx], below_me_sc);
	    }
	  /* Could we have jumped into this state from ROOT_S? Would it have
	   * been worth it (based on what I've seen so far) */
	  if(print_flag) 
	    printf("cur max_local_sc: %f\n\n", max_local_sc);
	  if(max_local_sc < (cm->beginsc[v] + lsc[tidx]))
	    {
	      max_local_sc = cm->beginsc[v] + lsc[tidx];
	      if(print_flag) 
		printf("\nNEW max_local_sc: %f\n\n", max_local_sc);
	    }	  
	}
    }
  free(v2n_map);
  free(v2n_ct);
  free(lsc);
  free(tr_esc);
  free(tr_tsc);
  /*printf("in ParsetreeScore_Global2Local() returning sc: %f\n", max_local_sc);*/
  return max_local_sc;

 ERROR: 
  cm_Fail("ERROR in ParsetreeScore_Global2Local()\n");
  return -1.;
}

/*
 * Function: Parsetree2CP9trace()
 * Incept:   EPN, Wed May 30 09:33:01 2007
 *
 * Purpose: Convert a CM parsetree into it's implicit CP9 trace.
 * Returns: eslOK on success
 *
 * Args:    
 * CM_t  *cm                - the CM, with valid cp9loc, and cp9glb if CM has local begins off
 * Parsetree_t *cm_tr       - valid parsetree to convert
 * cp9trace_s *ret_cp9_tr   - the CP9 trace to return, alloc'ed here
 */
int
Parsetree2CP9trace(CM_t *cm, Parsetree_t *tr, CP9trace_t **ret_cp9_tr)
{
  /* Check the contract */
  if(cm->cp9loc == NULL || (!(cm->flags & CMH_CP9LOC)))
    cm_Fail("In Parsetree2CP9trace, cm->cp9loc is not valid.\n");
  if(cm->cp9map == NULL)
    cm_Fail("In Parsetree2CP9trace, cm->cp9map is NULL.\n");

  int status;                    /* Easel status                            */
  CP9trace_t *cp9_tr;            /* the CP9 trace we're creating            */
  int **ks_ct = NULL;            /* [0..2][0..cp9->M] number of times each state was used
				  * 1st D: 0 = MATCH, 1 = INSERT, 2 = DELETE */
  int  tidx;                     /* counter over parsetree nodes */
  int  v;                        /* CM state index */
  int  k, ks;                    /* HMM nodes and state indices */
  int  i;                        /* generic counter */
  int  cp9_tr_size;              /* number of nodes we'll need for cp9_tr */
  int  lmost_k;                  /* left most HMM node visited in parse (often 1) */
  int  rmost_k;                  /* right most HMM node visited in parse (often M) */
  int  ip;
  int  ins_ct = 0;               /* total number of inserts */
  CP9_t *cp9 = NULL;

  cp9 = (cm->flags & CMH_LOCAL_BEGIN) ? cm->cp9loc : cm->cp9glb;
  if(cp9 == NULL) cm_Fail("In Parsetree2CP9trace, relevant cp9 does not exist\n");

  lmost_k = cp9->M + 1; 
  rmost_k = 0;              

  ESL_ALLOC(ks_ct,           sizeof(int *) * 3);
  for(ks = 0; ks < 3; ks++)
    {
      ESL_ALLOC(ks_ct[ks], sizeof(int) * (cp9->M+1));
      esl_vec_ISet(ks_ct[ks], cp9->M+1, 0);
    }

  /* Traverse parsetree, keeping track of implied HMM states used by each HMM node. */
  v = tr->state[0]; 
  if(v != 0) cm_Fail("ERROR in Parsetree2CP9Trace(), first Parsetree node not root.\n");
  /* we leave ks_ct[HMMMATCH][0] as 0 for convenience later, we know it was used. */

  for (tidx = 1; tidx < tr->n; tidx++) 
    {
      v  = tr->state[tidx];        	/* index of parent state in CM */
      for(i = 0; i < 2; i++) /* each CM state maps to 0, 1 or 2 HMM states */
	{ 
	  k  = cm->cp9map->cs2hn[v][i];
	  ks = cm->cp9map->cs2hs[v][i];
	  if(k == -1) continue; /* when HMM EL's are implemented, we'll have to have a special
				 * case for them, but for now we visit deletes in between. */
	  ks_ct[ks][k]++;
	  if(ks == HMMINSERT) ins_ct++;
	}
    }
  /* Determine the first (leftmost) node used and last (rightmost) node used, 
   * anything else was skipped by a smith-waterman local begin or end. */
  lmost_k = 1; 
  rmost_k = cp9->M; 
  if(cp9->flags & CPLAN9_LOCAL_BEGIN)
    {
      while((ks_ct[HMMMATCH][lmost_k] + ks_ct[HMMINSERT][lmost_k] + ks_ct[HMMDELETE][lmost_k]) == 0)
	lmost_k++;
    }
  if(cp9->flags & CPLAN9_LOCAL_END)
    {
      while((ks_ct[HMMMATCH][rmost_k] + ks_ct[HMMINSERT][rmost_k] + ks_ct[HMMDELETE][rmost_k]) == 0)
	rmost_k--;
    }
  /* Now build the CP9 trace */
  cp9_tr_size = (rmost_k - lmost_k + 1) + ins_ct + 2; /* number of match/deletes we'll visit plus
						       * number of inserts + begin/end */
  CP9AllocTrace(cp9_tr_size, &cp9_tr);  /* allow room for B & E */
  /* start at node 0 with the begin */
  cp9_tr->statetype[0] = CSTB;
  cp9_tr->nodeidx[0]   = 0;
  cp9_tr->pos[0]       = 0;
  tidx = 1;
  i    = 1;
  /* are there inserts from node 0? */
  for(ip = 0; ip < ks_ct[HMMINSERT][0]; ip++)
    {
      cp9_tr->statetype[tidx] = CSTI;
      cp9_tr->nodeidx[tidx]   = 0;
      cp9_tr->pos[tidx]       = i++;
      tidx++;
    }
  /* now go through nodes 1..M */
  for(k = lmost_k; k <= rmost_k; k++)
    {
      if(ks_ct[HMMMATCH][k])
	{
	  cp9_tr->statetype[tidx] = CSTM;
	  cp9_tr->nodeidx[tidx]   = k;
	  cp9_tr->pos[tidx]       = i++;
	  tidx++;
	}
      else if(ks_ct[HMMDELETE][k]) 
	{
	  cp9_tr->statetype[tidx] = CSTD;
	  cp9_tr->nodeidx[tidx]   = k;
	  cp9_tr->pos[tidx]       = 0;
	  tidx++;
	}
      else /* skipped due to local end, treat as delete for now */
	{
	  cp9_tr->statetype[tidx] = CSTD;
	  cp9_tr->nodeidx[tidx]   = k;
	  cp9_tr->pos[tidx]       = 0;
	  tidx++;
	}
      for(ip = 0; ip < ks_ct[HMMINSERT][k]; ip++)
	{
	  cp9_tr->statetype[tidx] = CSTI;
	  cp9_tr->nodeidx[tidx]   = k;
	  cp9_tr->pos[tidx]       = i++;
	  tidx++;
	}
    }
  /* all traces end with E state */
  cp9_tr->statetype[tidx]  = CSTE;
  cp9_tr->nodeidx[tidx]    = 0;
  cp9_tr->pos[tidx]        = 0;
  tidx++;
  cp9_tr->tlen = tidx;

  *ret_cp9_tr = cp9_tr;

  for(ks = 0; ks < 3; ks++)
    if(ks_ct[ks] != NULL) free(ks_ct[ks]);
  if(ks_ct != NULL) free(ks_ct);
  return eslOK;

 ERROR:
  for(ks = 0; ks < 3; ks++)
    if(ks_ct[ks] != NULL) free(ks_ct[ks]);
  if(ks_ct != NULL) free(ks_ct);
  return eslFAIL;
}

/* Function: rightjustify()
 * 
 * Purpose:  Given a gap-containing string of length n,
 *           pull all the non-gap characters as far as
 *           possible to the right, leaving gaps on the
 *           left side. Used to rearrange the positions
 *           of insertions in CM generated alignments.
 */
void
rightjustify(const ESL_ALPHABET *abc, char *s, int n)
{
  int npos;
  int opos;

  npos = n-1;
  opos = n-1;
  while (opos >= 0) {
    if (esl_abc_CIsGap(abc, s[opos]))
      opos--;
    else
      s[npos--]=s[opos--];  
  }
  while (npos >= 0) 
    s[npos--] = '.';
}

/* Function: leftjustify()
 * 
 * Purpose:  Given a gap-containing string of length n,
 *           pull all the non-gap characters as far as
 *           possible to the left, leaving gaps on the
 *           right side. Used to rearrange the positions
 *           of insertions in CM generated alignments.
 */
void
leftjustify(const ESL_ALPHABET *abc, char *s, int n)
{
  int npos;
  int opos;

  npos = 0;
  opos = 0;
  while (opos < n) {
    if (esl_abc_CIsGap(abc, s[opos]))
      opos++;
    else
      s[npos++]=s[opos++];  
  }
  while (npos < n) 
    s[npos++] = '.';
}



/* Function:  EmitParsetree()
 * Incept:    SRE, Mon Oct 13 22:35:46 2003 [Rams whupping Falcons, Monday Night Football]
 *            Easel'ed: EPN, Fri Aug  3 08:15:12 2007
 *
 * Purpose:   Sample a parsetree and sequence from the joint distribution
 *            Prob(sequence, parsetree | CM).
 *            
 *            Be careful screwing with the logic in here. You've got
 *            two tree traversals going simultaneously: a traversal of
 *            the CM, and a traversal of the growing parsetree. It
 *            wasn't obvious how to get it all to work in step
 *            together. Remember, one of your constraints is that the
 *            parsetree is numbered in preorder traversal - so you
 *            must push and defer the right child of a bifurcation,
 *            rather than attaching it immediately.  Another
 *            constraint is that you must set emitr in the parsetree
 *            even for nonemitting states, so you must always push a right
 *            marker along with a parsetree node index tpos, for deferred
 *            assignment of tr->emitr[tpos]. And since you don't
 *            know tpos until you've attached the state, you have to
 *            push the right marker after your deferred attachment of v - not 
 *            when v was produced - which is why you have a double
 *            deferral of the right emission or marker: you produce
 *            a V b, push that info onto the pda, pop it back off,
 *            attach V, store a, push b back onto the pda (now storing
 *            the trace position tpos for V), then produce from V.
 *            Yeesh.
 *
 *            Added capacity for local begins/ends. [EPN, Wed May  2 05:59:19 2007]
 *
 * Args:      cm      - covariance model to generate from
 *            errbuf  - for error messages
 *            r       - source of randomness
 *            name    - name for the sequence (ESL_SQ name field is mandatory)
 *            do_digital - TRUE to digitize sq before returning, FALSE not to
 *            ret_tr  - RETURN: generated parse tree. Pass NULL if unwanted.
 *            ret_sq  - RETURN: generated sequence
 *            ret_N   - RETURN: length of generated sequence.
 *
 * Returns:   eslOK on success; eslEMEM on memory error;
 *            eslEINCONCEIVABLE if something inconceivable happens.
 *            tr, sq are allocated here; whichever ones the caller
 *            requests (with non-NULL ret_ pointers) the caller is responsible
 *            for free'ing:
 *               FreeParsetree(tr); esl_sq_Destroy(sq);
 */
int
EmitParsetree(CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, char *name, int do_digital, Parsetree_t **ret_tr, ESL_SQ **ret_sq, int *ret_N)
{
  int status;
  Parsetree_t *tr = NULL;       /* parse tree under construction */
  ESL_STACK *pda = NULL;        /* pushdown automaton for traversing parse tree */              
  ESL_STACK *gsq = NULL;        /* growing sequence under construction */
  ESL_SQ    *sq  = NULL;        /* finished sequence, initially normal alphabet form */
  char      *seq;               /* alphabetic sequence to build sq with */
  int N;			/* current emitted sequence length */
  int tparent;			/* parent node index, last attached to parse tree */
  int tpos;			/* child node index, just attached to parse tree */
  int v;			/* index of current state */
  int y,z;			/* indices for next state(s)    */
  int type;			/* PDA_RESIDUE or PDA_STATE */
  int lchar, rchar;		/* index of emitted chars in cm->abc->sym[], or -1 for nothing */
  int whichway;			/* how to attach: TRACE_LEFT_CHILD or TRACE_RIGHT_CHILD */
  int x;			/* tmp variable for sampling MP emission */
  int lpos;                     /* tmp variable for inserting EL trace node */
  float *tmp_tvec = NULL;       /* tmp transition vector to choose from, 
				 * for dealing with local end transitions */
  /* Contract check */
  if(cm->flags & CMH_LOCAL_END && (fabs(sreEXP2(cm->el_selfsc) - 1.0) < 0.01))
    ESL_FAIL(eslEINVAL, errbuf, "EL self transition probability %f is too high, would emit long (too long) EL insertions.", sreEXP2(cm->el_selfsc));
  if(cm->abc == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "CM does not have a valid alphabet.");
  if(ret_sq != NULL && name == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "EmitParsetree requires a sequence name for the sequence it's creating.");

  tr  = CreateParsetree(100);
  if((pda = esl_stack_ICreate()) == NULL) goto ERROR;
  if((gsq = esl_stack_CCreate()) == NULL) goto ERROR;
  N   = 0;			
  ESL_ALLOC(tmp_tvec, sizeof(float) * (MAXCONNECT+1)); /* enough room for max possible transitions, plus
							* a local end transition */
  /* Init by pushing root state's info onto pda
   */
  if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* doesn't emit an rchar */
  if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* doesn't emit an lchar either */
  if((status = esl_stack_IPush(pda, TRACE_LEFT_CHILD)) != eslOK) goto ERROR;
  if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* attach this state to parsetree node -1 (init) */  
  if((status = esl_stack_IPush(pda, 0)) != eslOK) goto ERROR;		/* it's the root state, v=0 */
  if((status = esl_stack_IPush(pda, PDA_STATE)) != eslOK) goto ERROR;

  /* Iterate until the pda is empty...
   */
  while (esl_stack_IPop(pda, &type) != eslEOD) 
    {
      if (type == PDA_RESIDUE)
	{
	  esl_stack_IPop(pda, &tpos);
	  esl_stack_IPop(pda, &rchar);

	  if (rchar != -1) {
	    if((status = esl_stack_CPush(gsq, cm->abc->sym[rchar])) != eslOK) goto ERROR;
	    N++;
	  }
	  tr->emitr[tpos] = N;
	}
      else if (type == PDA_STATE) 
	{
	  esl_stack_IPop(pda, &v);
	  esl_stack_IPop(pda, &tparent);
	  esl_stack_IPop(pda, &whichway);
	  esl_stack_IPop(pda, &lchar);
	  esl_stack_IPop(pda, &rchar);

	  /* Attach state v to the parent parsetree node that generated it,
	   * which is tparent. Set emitl now; emitr gets deferred and set later.
           * The insertion function returns tpos, the index of the node in the
           * parse tree that we just created.
	   */
	  tpos = InsertTraceNode(tr, tparent, whichway, N+1, -1, v);	    

	  /* If v emitted left: add that symbol to the growing seq.
	   */
	  if (lchar != -1)
	    {
	      if((status = esl_stack_CPush(gsq, cm->abc->sym[lchar])) != eslOK) goto ERROR;
	      N++;
	    }

	  /* Push right emission info for state v onto the pda, now
           * that we know tpos for where v is in the parsetree. We have
           * to do this even if rchar is -1, to be sure that we will set the emitr
           * bound properly even for nonemitting states in the parsetree.
	   */
	  if((status = esl_stack_IPush(pda, rchar)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, tpos)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, PDA_RESIDUE)) != eslOK) goto ERROR;

	  /* Decide what state we're going to next.
           * B is special case of a bifurcation to two S states. 
	   */
	  if (cm->sttype[v] == B_st)
	    {
	      y = cm->cfirst[v];	/* left child  */
	      z = cm->cnum[v];	        /* right child */
	  
	      /* Push the right start state's info
	       */
	      if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* doesn't emit right */
	      if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* doesn't emit left */
	      if((status = esl_stack_IPush(pda, TRACE_RIGHT_CHILD)) != eslOK) goto ERROR; /* attach as right child of the B */
	      if((status = esl_stack_IPush(pda, tpos)) != eslOK) goto ERROR;		/* attach it to B, which is tpos in parsetree*/
	      if((status = esl_stack_IPush(pda, z)) != eslOK) goto ERROR;		/* state z */
	      if((status = esl_stack_IPush(pda, PDA_STATE)) != eslOK) goto ERROR;

	      /* Push the left start state's info
	       */
	      if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* doesn't emit right */
	      if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* doesn't emit left */
	      if((status = esl_stack_IPush(pda, TRACE_LEFT_CHILD)) != eslOK) goto ERROR; /* attach as left child of the B */
	      if((status = esl_stack_IPush(pda, tpos)) != eslOK) goto ERROR;		/* attach it to B, which is tpos in parsetree*/
	      if((status = esl_stack_IPush(pda, y)) != eslOK) goto ERROR;		/* state z */
	      if((status = esl_stack_IPush(pda, PDA_STATE)) != eslOK) goto ERROR;
	    }
	  else
	    {
	      if(v == 0 && cm->flags & CMH_LOCAL_BEGIN)	{ /* ROOT_S with local begins, special */
		if(cm->flags & CM_EMIT_NO_LOCAL_BEGINS) { /* even though local begins are on, we don't allow them during emission */
		  if(cm->root_trans == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "EmitParsetree(), cm->flags CM_EMIT_NO_LOCAL_BEGINS and CM_EMIT_GLOBAL flags raised, but cm->root_trans is NULL.");
		  y = cm->cfirst[v] + esl_rnd_FChoose(r, cm->root_trans, cm->cnum[0]); /* choose next state, y, from 0's children using initial transitions (those from global model, read from CM file) */
		}		  
		else { 
		  y = esl_rnd_FChoose(r, cm->begin, cm->M); /* choose next state, y */
		}
	      }
	      else if(cm->flags & CMH_LOCAL_END) /* special case, we could transit to EL, if CM_EMIT_NO_LOCAL_ENDS flag is down */
		{
		  if(cm->flags & CM_EMIT_NO_LOCAL_ENDS) { /* even though local ends are on, we dont allow them during emission */
		    /* create temporary vector for choosing transition */
		    esl_vec_FSet(tmp_tvec, (MAXCONNECT+1), 0.);
		    esl_vec_FCopy(cm->t[v], cm->cnum[v], tmp_tvec);
		    esl_vec_FNorm(tmp_tvec, cm->cnum[v]);
		    y = cm->cfirst[v] + esl_rnd_FChoose(r, tmp_tvec, cm->cnum[v]); /* choose next state, y, but don't include a local end as a possibility */
		  }		  
		  else { /* we may choose a child of v, or a local end */
		    esl_vec_FSet(tmp_tvec, (MAXCONNECT+1), 0.);
		    esl_vec_FCopy(cm->t[v], cm->cnum[v], tmp_tvec);
		    tmp_tvec[cm->cnum[v]] = cm->end[v];
		    y = esl_rnd_FChoose(r, tmp_tvec, (cm->cnum[v]+1)); /* choose next state, y's offset */
		    if(y == cm->cnum[v]) y = cm->M; /* local end */
		    else y += cm->cfirst[v];        
		  }
		}		  
	      else
		y = cm->cfirst[v] + esl_rnd_FChoose(r, cm->t[v], cm->cnum[v]); /* choose next state, y */

	      switch (cm->sttype[y]) {
	      case MP_st: 
		x     = esl_rnd_FChoose(r, cm->e[y], cm->abc->K*cm->abc->K);
		lchar = x / cm->abc->K;
		rchar = x % cm->abc->K;
		break;
	      case ML_st:
	      case IL_st:
		lchar = esl_rnd_FChoose(r, cm->e[y], cm->abc->K);
		rchar = -1;
		break;
	      case MR_st:
	      case IR_st:
		lchar = -1;
		rchar = esl_rnd_FChoose(r, cm->e[y], cm->abc->K);
		break;
	      case EL_st: /* EL emits on transition, here we don't emit */
		lchar = -1;
		rchar = -1;
		break;
	      default:
		lchar = -1;
		rchar = -1;
	      }
	      if (cm->sttype[y] == E_st)
		{
		  /*InsertTraceNode(tr, tpos, TRACE_LEFT_CHILD, -1, -1, y);*/
		  InsertTraceNode(tr, tpos, TRACE_LEFT_CHILD, N+1, N, y);
		} 
	      else if(cm->sttype[y] == EL_st) /* y == cm->M */
		{
		  lpos = N+1; /* remember lpos, we need it after we emit from EL */
		  /* Now choose number of residues emitted from EL, could be 0.
		   * We do this here b/c convention for EL is to have a single trace node,
		   * even if multiple residues are emitted. */
		  esl_vec_FSet(tmp_tvec, (MAXCONNECT+1), 0.);
		  tmp_tvec[0] = sreEXP2(cm->el_selfsc); /* EL self probability */
		  tmp_tvec[1] = 1. - tmp_tvec[0];       /* probability of going to implicit END */
		  y = esl_rnd_FChoose(r, tmp_tvec, 2); /* choose next state, either EL or implicit END */
		  while(y == 0) /* we've self-transitioned, emit 1 res from NULL distro */
		    {
		      lchar = esl_rnd_FChoose(r, cm->null, cm->abc->K);
		      if((status = esl_stack_CPush(gsq, cm->abc->sym[lchar])) != eslOK) goto ERROR;
		      N++;
		      y = esl_rnd_FChoose(r, tmp_tvec, 2); /* choose next state, either EL or implicit END */
		    }
		  InsertTraceNode(tr, tpos, TRACE_LEFT_CHILD, lpos, N, cm->M); /* careful to reset y to cm->M */
		}
	      else 
		{
		  if((status = esl_stack_IPush(pda, rchar)) != eslOK) goto ERROR;		/* does it emit right? */
		  if((status = esl_stack_IPush(pda, lchar)) != eslOK) goto ERROR;		/* does it emit left? */
		  if((status = esl_stack_IPush(pda, TRACE_LEFT_CHILD)) != eslOK) goto ERROR; /* non-B's: attach as left child by conv */
		  if((status = esl_stack_IPush(pda, tpos)) != eslOK) goto ERROR;		/* attach it to v, which is tpos in parsetree*/
		  if((status = esl_stack_IPush(pda, y)) != eslOK) goto ERROR;		/* next state we're going to */
		  if((status = esl_stack_IPush(pda, PDA_STATE)) != eslOK) goto ERROR;
		}
	    } /* end of PDA_STATE logic */  
	} /* end of else (which we enter if v not a B state) */
    } /* end of main "while esl_stack_IPop()" loop */

  if((seq = esl_stack_Convert2String(gsq)) == NULL) goto ERROR; /* this destroys gsq char stack */
  if(name != NULL) sq  = esl_sq_CreateFrom(name, seq, NULL, NULL, NULL);
  else             sq  = esl_sq_CreateFrom("seq", seq, NULL, NULL, NULL); 
  if(sq == NULL) goto ERROR;
  free(seq); /* we made a copy of this when creating sq */
  /* name can only be NULL if ret_sq == NULL, so we're throwing it away anyway */

  /* digitize if nec */
  if(do_digital) 
    if((status = esl_sq_Digitize(cm->abc, sq)) != eslOK) goto ERROR;
  /*ParsetreeDump(stdout, tr, cm, dsq);*/ 

  esl_stack_Destroy(pda);

  free(tmp_tvec);
  if (ret_tr  != NULL) *ret_tr  = tr;  else FreeParsetree(tr);
  if (ret_sq  != NULL) *ret_sq  = sq;  else esl_sq_Destroy(sq);
  if (ret_N   != NULL) *ret_N   = N; 
  return eslOK;
  
 ERROR:
  if(tr  != NULL) FreeParsetree(tr);
  if(gsq != NULL) esl_stack_Destroy(gsq);
  if(pda != NULL) esl_stack_Destroy(pda);
  if(sq  != NULL) esl_sq_Destroy(sq);
  if(tmp_tvec != NULL) free(tmp_tvec);
  return status;
}
  
/* Function: ParsetreeScoreCorrectionNull2()
 * based on     TraceScoreCorrection() from HMMER:
 * EPN 08.24.06 Janelia
 * 
 * Purpose:  Calculate a correction (in integer log_2 odds) to be
 *           applied to a sequence, using a second null model, 
 *           based on a traceback. All emissions are corrected;
 *           The null model is constructed /post hoc/ as the
 *           average over all the emission distributions used by the trace.
 *           
 *           <omega> is the prior probability of the null3 model. 
 *           The log of this value is the penalty for the null3 
 *           correction. If <omega> is 1./65536., which it is 
 *           by default (1/2^16), we apply a 16 bit penalty.
 *
 * Return:   ret_sc: the log_2-odds score correction.          
 *           eslEINCOMPAT on contract violation
 */
int 
ParsetreeScoreCorrectionNull2(CM_t *cm, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int start, float omega, float *ret_sc)
{
  int status;
  float *p;		/* null2 model distribution */
  float *sc;	        /* null2 model scores       */
  int   a,b;            /* residue index counters */
  int   v;              /* state index counter */
  int   i, j;           /* seq posn counter */
  int   tidx;
  float score;
  float struct_score;   /* structure contribution to the score */

  if(ret_sc == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "ParsetreeScoreCorrectionNull2() ret_sc is NULL.");

  /* Rarely, the alignment was totally impossible, and tr is NULL.
   */
  if (tr == NULL) return 0.0;
  
  /* Set up model: average over the emission distributions of
   * all M, I states that appear in the trace. Ad hoc? Sure, you betcha. 
   */
  /* trivial preorder traverse, since we're already numbered that way */
  ESL_ALLOC(p, sizeof(float) * cm->abc->K);
  esl_vec_FSet(p, cm->abc->K, 0.0);
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    if(cm->sttype[v] == MP_st) { 
      /* we treat this as two match states. */
      for(a = 0; a < cm->abc->K; a++) { 
	/* first add contribution to null2 for left half. */
	for(b = (a * cm->abc->K); b < ((a+1) * cm->abc->K); b++) p[a] += cm->e[v][b]; 
	/* now add contribution for right half. */
	for(b = a; b < (cm->abc->K * cm->abc->K); b += cm->abc->K) p[a] += cm->e[v][b]; 
      }
    }
    else if(cm->sttype[v] == ML_st || cm->sttype[v] == IL_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
      esl_vec_FAdd(p, cm->e[v], cm->abc->K);
    }
  }
  esl_vec_FNorm(p, cm->abc->K);

  ESL_ALLOC(sc,  sizeof(float) * (cm->abc->Kp));
  /* calculate null2 scores of each possible emission, first the base alphabet */
  for (a = 0; a < cm->abc->K; a++) sc[a] = sreLOG2(p[a] / cm->null[a]);
  /* the ambiguities */
  for (a = cm->abc->K+1; a < cm->abc->Kp-1; a++) sc[a] = esl_abc_FAvgScore(cm->abc, a, sc);  

  /* Score all the state emissions that appear in the trace.
   */
  score = struct_score = 0;
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    i = tr->emitl[tidx];
    j = tr->emitr[tidx];
    if (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) score += sc[dsq[i+start-1]];
    if (cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) score += sc[dsq[j+start-1]];
  }
   /* Apply an ad hoc log(omega) fudge factor penalty;
    * interpreted as a prior, saying that the third null model is 
    * <omega> as likely as the standard null model. 
    * <omega> is by default 1/(2^16), so this is by default a
    * 16 bit penalty.
    */
  score += sreLOG2(omega);
  
  /* Return the correction to the bit score. */
  ESL_DPRINTF1(("ParsetreeScoreCorrectionNull2 return sc: %f\n", LogSum2(0., score)));
  free(sc);
  free(p);
  score = LogSum2(0., score);
  *ret_sc = score;
  return eslOK;
  
 ERROR:
  ESL_FAIL(status, errbuf, "ParsetreeScoreCorrectionNull2(): memory allocation error.");
  return status; /* NEVERREACHED*/
}

  
/* Function: ParsetreeScoreCorrectionNull3()
 * Incept:   EPN, Sat May  3 15:38:24 2008
 * 
 * Purpose:  Calculate a correction (in integer log_2 odds) to be
 *           applied to a sequence, using a third null model, the
 *           composition of the target sequence. 
 *           All emissions are corrected;
 *           The null model is constructed /post hoc/ as the
 *           distribution of the target sequence; if the target
 *           sequence is 40% A, 5% C, 5% G, 40% U, then the null 
 *           model is (0.4, 0.05, 0.05, 0.4).
 *           
 *           NOTE: (start) is offset in dsq such that tr->emitl[0] corresponds
 *           to the residue in dsq[1];
 *
 *           <omega> is the prior probability of the null3 model. 
 *           The log of this value is the penalty for the null3 
 *           correction. If <omega> is 1./65536., which it is 
 *           by default (1/2^16), we apply a 16 bit penalty.
 * 
 * Return:   ret_sc: the log_2-odds score correction.          
 *           eslEINCOMPAT on contract violation
 */
int 
ParsetreeScoreCorrectionNull3(CM_t *cm, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int start, float omega, float *ret_sc)
{
  int status;
  float *p;		/* null3 model distribution */
  float *sc;	        /* null3 model scores       */
  int   a;              /* residue index counters */
  int   v;              /* state index counter */
  int   i, j;           /* seq posn counter */
  int   tidx;
  float score;
  float struct_score;   /* structure contribution to the score */
  
  if(ret_sc == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "ParsetreeScoreCorrectionNull3() ret_sc is NULL.");
  /* Rarely, the alignment was totally impossible, and tr is NULL.
   */
  if (tr == NULL) return 0.0;
  
  /* get composition of full subseq in parse from tr->emitl[0]..tr->emitr[0], 
   * starting at dsq+start-1 b/c coords in tr->emit* are offset relative to dsq by start
   * such that tr->emitl[0] always equals 1. 
   * Note: we're INcluding any EL emissions here in ACGU composition calculation but then
   * EXcluding them when we correct the score below (not sure if this is right),
   * in a way we're always assuming scores of ELs are 0.0, meaning there's no difference
   * between the probability they're emitted by the model and any possible NULL model.
   */
  get_alphabet_comp(cm->abc, dsq+start-1, tr->emitl[0], tr->emitr[0], &p);
  ESL_ALLOC(sc,  sizeof(float) * (cm->abc->Kp));
  /* calculate null3 scores of each possible emission, first the base alphabet */
  for (a = 0; a < cm->abc->K; a++) { 
    sc[a] = sreLOG2(p[a] / cm->null[a]);
    /*printf("p[%d]: %.3f sc %.3f\n", a, p[a], sc[a]);*/
  }
  /* the ambiguities */
  for (a = cm->abc->K+1; a < cm->abc->Kp-1; a++) sc[a] = esl_abc_FAvgScore(cm->abc, a, sc);  

  /* Score all the state emissions that appear in the trace.
   */
  score = struct_score = 0.;
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    i = tr->emitl[tidx];
    j = tr->emitr[tidx];
    if (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) score += sc[dsq[i+start-1]];
    if (cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) score += sc[dsq[j+start-1]];
  }
   /* Apply an ad hoc log(omega) fudge factor penalty;
    * interpreted as a prior, saying that the third null model is 
    * <omega> as likely as the standard null model. 
    * <omega> is by default 1/(2^16), so this is by default a
    * 16 bit penalty.
    */
  score += sreLOG2(omega);

  /* Return the correction to the bit score. */
  /*printf("ParsetreeScoreCorrectionNull3 return sc: %f\n", LogSum2(0., score));*/
  ESL_DPRINTF1(("ParsetreeScoreCorrectionNull3 return sc: %f\n", LogSum2(0., score)));
  free(sc);
  free(p);
  score = LogSum2(0., score);
  *ret_sc = score;
  return eslOK;
  
 ERROR:
   ESL_FAIL(status, errbuf, "ParsetreeScoreCorrectionNull3(): memory allocation error.");
   return status; /* NEVERREACHED*/
}
  
/* Function: ScoreCorrectionNull3()
 * Incept:   EPN, Sat May 10 17:58:03 2008
 * 
 * Purpose:  Calculate a correction (in integer log_2 odds) to be
 *           applied to a sequence, using a third null model, the
 *           composition of the target sequence. 
 *           All emissions are corrected;
 *           The null model is constructed /post hoc/ as the
 *           distribution of the target sequence; if the target
 *           sequence is 40% A, 5% C, 5% G, 40% U, then the null 
 *           model is (0.4, 0.05, 0.05, 0.4).
 * 
 *           Note: no trace or parsetree is needed. The bit score correction
 *           can be derived solely by the nucleotide composition of the hit and
 *           it's length. 
 *           
 *           <omega> is the prior probability of the null3 model. 
 *           The log of this value is the penalty for the null3 
 *           correction. If <omega> is 1./65536., which it is 
 *           by default (1/2^16), we apply a 16 bit penalty.
 * 
 * Args:     abc  - alphabet for hit (only used to get alphabet size, which is size of <comp>)
 *           null0- the first null model used when building the CM, usually cm->null or cm->cp9->null
 *           comp - [0..a..abc->K-1] frequency of residue a in the hit we're correcting the score for, 
 *                  passed in b/c we can efficiently compute this during scanning DP funcs instead of
 *                  calcing it each time per hit which is wasteful for many, possibly overlapping hits
 *                  which is the case during model calibration with cmcalibrate.
 *           len  - length of the hit
 *           ret_sc- RETURN: the correction to the score, caller subtracts this from hit score to get 
 *                   corrected score.
 *
 * Return:   void, ret_sc: the log_2-odds score correction.          
 */
void
ScoreCorrectionNull3(const ESL_ALPHABET *abc, float *null0, float *comp, int len, float omega, float *ret_sc)
{
  int   a;              /* residue index counters */
  float score = 0.;

  /*printf("\n");
    esl_vec_FDump(stdout, comp, abc->K, NULL);*/
  
  for (a = 0; a < abc->K; a++) score += sreLOG2(comp[a] / null0[a]) * comp[a] * len;

   /* Apply an ad hoc log(omega) fudge factor penalty;
    * interpreted as a prior, saying that the third null model is 
    * <omega> as likely as the standard null model. 
    * <omega> is by default 1/(2^16), so this is by default a
    * 16 bit penalty.
    */
  score += sreLOG2(omega);

  /* Return the correction to the bit score. */

  /* Return the correction to the bit score. */
  /*printf("ScoreCorrectionNull3 return sc: %.3f\n", LogSum2(0., score));*/
  ESL_DPRINTF3(("ScoreCorrectionNull3 return sc: %f\n", LogSum2(0., score)));
  score = LogSum2(0., score);
  *ret_sc = score;
  return;
}

  
/* Function: ScoreCorrectionNull3CompUnknown()
 * Incept:   EPN, Thu May 22 13:16:04 2008
 * 
 * Purpose:  Calculate a correction (in integer log_2 odds) to be
 *           applied to a sequence, using a third null model, the
 *           composition of the target sequence. 
 *           All emissions are corrected;
 *           The null model is constructed /post hoc/ as the
 *           distribution of the target sequence; if the target
 *           sequence is 40% A, 5% C, 5% G, 40% U, then the null 
 *           model is (0.4, 0.05, 0.05, 0.4).
 * 
 *           Same as ScoreCorrectionNull3() except that no <comp> vector is needed,
 *           the composition is determined within this function. 
 *           
 *           <omega> is the prior probability of the null3 model. 
 *           The log of this value is the penalty for the null3 
 *           correction. If <omega> is 1./65536., which it is 
 *           by default (1/2^16), we apply a 16 bit penalty.
 * 
 * Args:     abc  - alphabet for hit (only used to get alphabet size, which is size of <comp>)
 *           null0- the first null model used when building the CM, usually cm->null or cm->cp9->null
 *           dsq  - the sequence the hit resides in
 *           start- start position of hit in dsq
 *           end  - end   position of hit in dsq
 *           ret_sc- RETURN: the correction to the score, caller subtracts this from hit score to get 
 *                   corrected score.
 * Return:   void, ret_sc: the log_2-odds score correction.          
 */
void
ScoreCorrectionNull3CompUnknown(const ESL_ALPHABET *abc, float *null0, ESL_DSQ *dsq, int start, int stop, float omega, float *ret_sc)
{
  float score = 0.;
  float *comp;		/* null3 model distribution */

  get_alphabet_comp(abc, dsq, start, stop, &comp);
  ScoreCorrectionNull3(abc, null0, comp, (stop-start+1), omega, &score);
  free(comp);
  *ret_sc = score;
  return;
}


    
/* Function: ParsetreeCountMPEmissions()
 * Date:     EPN, Thu May 22 14:11:28 2008
 *
 * Purpose:  Given a parsetree, return the number of residues emitted by MP states.
 *
 * Returns:  Number of residues emitted by MP states in <tr>.
 */
int 
ParsetreeCountMPEmissions(CM_t *cm, Parsetree_t *tr)
{
  int tidx;
  int nres_by_mp = 0;

  for (tidx = 0; tidx < tr->n; tidx++) {  
    if(cm->sttype[tr->state[tidx]] == MP_st) nres_by_mp += 2;
  }
  return nres_by_mp;
}

/* Function: Alignment2Parsetrees()
 * EPN, Fri Jul 11 09:49:50 2008
 *
 * Purpose:  Given a MSA <msa>, a CM <cm> and a guidetree <mtr> for <cm>,
 *           Determine the implicit parsetrees of the sequences in the
 *           MSA to the CM. Return the parsetrees in <ret_tr> if non-NULL, 
 *           sequence objects in <ret_sq> if non-null. 
 *
 *           Dealign the MSA seqs in <ret_sq> and convert from aligned to
 *           unaligned coordinates in <ret_tr>.
 *
 * Args:     msa          - MSA we want to infer parsetrees from
 *           cm           - CM we're aligning to 
 *           mtr          - master parsetree, guide tree for CM 
 *           errbuf       - easel error message
 *           ret_sq       - Return: dealigned msa seqs in digital form
 *           ret_tr       - Return: parsetree for seqs in dealigned coords
 * 
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation, eslEMEM on memory error
 *           <ret_tr>, <ret_sq>, see 'Purpose'.
 */
int 
Alignment2Parsetrees(ESL_MSA *msa, CM_t *cm, Parsetree_t *mtr, char *errbuf, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr)
{
  int           status;
  int           i;	      /* counter over aseqs       */
  int           apos;         /*   aligned position index */
  int           uapos;        /* unaligned position index */
  int           x;            /* counter of parsetree nodes */
  int          *map   = NULL; /* for current seq, [0..msa->alen] map from aligned posns to unaligned (non-gap) posns */
  char         *uaseq = NULL; /* current seq, dealigned from the MSA */
  char         *aseq  = NULL; /* current seq, aligned text */
  Parsetree_t **tr    = NULL; /* [0..msa->nseq-1] new parsetrees, one per seq in msa */
  ESL_SQ      **sq    = NULL; /* [0..msa->nseq-1] new ESL_SQ objects, one per seq in msa */

  /* Contract check */
  if(msa == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "Alignment2Parsetrees() msa is NULL.\n");
  if(! (msa->flags & eslMSA_DIGITAL))  ESL_FAIL(eslEINCOMPAT, errbuf, "Alignment2Parsetrees() msa is not digitized.\n");
  if(ret_tr == NULL && ret_sq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "Alignment2Parsetrees() ret_sq and ret_tr both NULL.");

  if(ret_tr != NULL) ESL_ALLOC(tr, (sizeof(Parsetree_t *) * msa->nseq));
  if(ret_sq != NULL) ESL_ALLOC(sq, (sizeof(ESL_SQ *)      * msa->nseq));
  ESL_ALLOC(aseq,  sizeof(char) * (msa->alen+1));
  ESL_ALLOC(map,   sizeof(int)  * (msa->alen+1));
  map[0] = -1; /* invalid */

  for (i = 0; i < msa->nseq; i++) { 
    uapos = 1;
    /* map aligned to dealigned coords (digitized coords, 1..alen) for this seq
     * map is needed b/c we want the parsetree in dealigned coords so we can
     * call Parsetrees2Alignment with it, but mtr is in aligned coords, and
     * Transmogrify works in aligned coords, so after calling Transmogrify
     * we have to convert tr->emitl and tr->emitr to dealigned coords using map.
     */
    for(apos = 1; apos <= msa->alen; apos++) 
      map[apos] = esl_abc_XIsGap(msa->abc, msa->ax[i][apos]) ? -1 : uapos++;
    /* get text seq, we need digitized AND text seqs for Transmogrify */
    esl_abc_Textize(msa->abc, msa->ax[i], msa->alen, aseq);
    esl_strdup(aseq, -1, &uaseq);
    /* dealign seq */
    esl_strdealign(uaseq, uaseq, "-_.~", NULL);
    /* Transmogrify the aligned seq to get a parsetree */
    if(ret_tr != NULL) { 
      tr[i] = Transmogrify(cm, mtr, msa->ax[i], aseq, msa->alen);
      /*ParsetreeDump(stdout, tr[i], cm, msa->ax[i], NULL, NULL);*/
      /* tr[i] is in alignment coords, convert it to unaligned coords, */
      for(x = 0; x < tr[i]->n; x++) { 
	/*printf("i: %d x: %d emitl %d emitr %d\n", i, x, tr[i]->emitl[x], tr[i]->emitr[x]);*/
	if(tr[i]->emitl[x] != -1) { 
	  /*printf("\tmapl: %d\n", map[i][tr[i]->emitl[x]]);*/
	  tr[i]->emitl[x] = map[tr[i]->emitl[x]];
	}
	if(tr[i]->emitr[x] != -1) { 
	  /*printf("\tmapr: %d\n", map[i][tr[i]->emitr[x]]);*/
	  tr[i]->emitr[x] = map[tr[i]->emitr[x]];
	}
      }
    }
    if(ret_sq != NULL) { 
      sq[i] = esl_sq_CreateFrom(msa->sqname[i], uaseq, NULL, NULL, NULL);
      esl_sq_Digitize(cm->abc, sq[i]);
    }
    free(uaseq); /* this gets reallocated and filled per seq in esl_strdup() call above */
  }
  free(aseq);
  free(map);

  /* tr and sq are only allocated if ret_tr and ret_sq were non-null */
  if(ret_tr != NULL) *ret_tr = tr;
  if(ret_sq != NULL) *ret_sq = sq;

  return eslOK;

 ERROR:
  if(map != NULL )  free(map);
  if(uaseq != NULL) free(uaseq);
  if(aseq != NULL)  free(aseq);
  return status;
}


/* Function: ParsetreeMode()
 * Incept:   EPN, Thu Nov 10 11:31:04 2011
 * 
 * Purpose:  Return the alignment mode of a parsetree.
 */
char 
ParsetreeMode(Parsetree_t *tr)
{
  return tr->mode[0];
}

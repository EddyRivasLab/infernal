/* parsetree.c
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
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_sqio.h"

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
  new->memblock = 100;		/* allocation block size can be optimized here if you want. */
  new->nalloc   = size;
  ESL_ALLOC(new->emitl, sizeof(int) * new->nalloc);
  ESL_ALLOC(new->emitr, sizeof(int) * new->nalloc);
  ESL_ALLOC(new->state, sizeof(int) * new->nalloc);
  ESL_ALLOC(new->mode,  sizeof(int) * new->nalloc);
  ESL_ALLOC(new->nxtl,  sizeof(int) * new->nalloc);
  ESL_ALLOC(new->nxtr,  sizeof(int) * new->nalloc);
  ESL_ALLOC(new->prv,   sizeof(int) * new->nalloc);
  new->n = 0;
  return new;
 ERROR:
  esl_fatal("ERROR allocated parsetree.\n");
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
  ESL_RALLOC(tr->mode,  tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->nxtl,  tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->nxtr,  tmp, sizeof(int) * tr->nalloc);
  ESL_RALLOC(tr->prv,   tmp, sizeof(int) * tr->nalloc);
  return;
  
 ERROR:
  esl_fatal("ERROR growing parsetree.\n");
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
 * Returns:  index of new node.
 */          
int
InsertTraceNodewithMode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state, int mode)
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

   n = InsertTraceNodewithMode(tr, y, whichway, emitl, emitr, state, 3);

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
 *           given a CM that's prepared in log-odds form.
 */
float
ParsetreeScore(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, int do_null2)
{
  /* contract check */
  if(dsq == NULL)
    esl_fatal("ERROR in ParsetreeScore(), dsq is NULL.\n");

  int tidx;			/* counter through positions in the parsetree        */
  int v,y;			/* parent, child state index in CM                   */
  ESL_DSQ symi, symj;		/* symbol indices for emissions, 0..cm->abc->Kp-1    */
  float sc;			/* the log-odds score of the parse tree */
  int mode;

		/* trivial preorder traverse, since we're already numbered that way */
  sc = 0.;
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    mode = tr->mode[tidx];
    if (v == cm->M) continue;      	/* special case: v is EL, local alignment end */
    if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) /* no scores in B,E */
      {
	y = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */

        if (tr->nxtl[tidx] == -1)
          ;
	else if (v == 0 && (cm->flags & CM_LOCAL_BEGIN))
	  sc += cm->beginsc[y];
	else if (y == cm->M) /* CM_LOCAL_END is presumably set, else this wouldn't happen */
	  sc += cm->endsc[v] + (cm->el_selfsc * (tr->emitr[tidx] - tr->emitl[tidx] + 1 - StateDelta(cm->sttype[v])));
	else 		/* y - cm->first[v] gives us the offset in the transition vector */
	  sc += cm->tsc[v][y - cm->cfirst[v]];
	
	if (cm->sttype[v] == MP_st) 
	  {
	    symi = dsq[tr->emitl[tidx]];
	    symj = dsq[tr->emitr[tidx]];
            if (mode == 3)
              {
  	        if (symi < cm->abc->K && symj < cm->abc->K)
	          sc += cm->esc[v][(int) (symi*cm->abc->K+symj)];
	        else
	          sc += DegeneratePairScore(cm->abc, cm->esc[v], symi, symj);
              }
            else if (mode == 2)
              sc += LeftMarginalScore(cm->abc, cm->esc[v], symi);
            else if (mode == 1)
              sc += RightMarginalScore(cm->abc, cm->esc[v], symj);
	  } 
	else if ( (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) && (mode == 3 || mode == 2) )
	  {
	    symi = dsq[tr->emitl[tidx]];
	    if (symi < cm->abc->K) sc += cm->esc[v][(int) symi];
	    else                   sc += esl_abc_FAvgScore(cm->abc, symi, cm->esc[v]);
	  } 
	else if ( (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) && (mode == 3 || mode == 1) )
	  {
	    symj = dsq[tr->emitr[tidx]];
	    if (symj < cm->abc->K) sc += cm->esc[v][(int) symj];
	    else                   sc += esl_abc_FAvgScore(cm->abc, symj, cm->esc[v]);
	  }
      }
  }

  if(do_null2)
    sc -= CM_TraceScoreCorrection(cm, tr, dsq);

  return sc;
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

  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s\n",
	  " idx ","emitl", "emitr", "state", " nxtl", " nxtr", " prv ");
  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s\n",
	 "-----", "-----", "-----", "-----", "-----","-----", "-----");
  for (x = 0; x < tr->n; x++)
    fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d\n",
	   x, tr->emitl[x], tr->emitr[x], tr->state[x], 
	   tr->nxtl[x], tr->nxtr[x], tr->prv[x]);
  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s\n",
	 "-----", "-----", "-----", "-----","-----","-----", "-----");

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
 *          gamma - cumulative subsequence length probability distributions
 *                  used to generate the bands; from BandDistribution(); [0..v..M-1][0..W]
 *          W     - maximum window length W (gamma distributions range up to this)        
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
  int   mode;
  int   do_banded;
  int   L;

  /* Contract check */
  if(dmin == NULL && dmax != NULL)
    esl_fatal("In ParsetreeDump(), dmin is NULL, dmax is not.\n");
  if(dmin != NULL && dmax == NULL)
    esl_fatal("In ParsetreeDump(), dmax is NULL, dmin is not.\n");
  if(dsq == NULL)
    esl_fatal("In ParsetreeDump(), dsq is NULL");

  if(dmin != NULL && dmax != NULL) do_banded = TRUE;
  else                             do_banded = FALSE;

  if(do_banded)
    {
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	      " idx ", "emitl", "emitr", "state", " nxtl", " nxtr", " prv ", " tsc ", " esc ", 
	      " L   ", " dmin", " dmax");
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	      "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----",
	      "-----", "-----", "-----");
    }
  else
    {
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s\n",
	      " idx ","emitl", "emitr", "state", " nxtl", " nxtr", " prv ", " tsc ", " esc ");
      fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s\n",
	      "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----");
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
	if (mode == 3 || mode == 2) syml = cm->abc->sym[dsq[tr->emitl[x]]]; 
	if (mode == 3 || mode == 1) symr = cm->abc->sym[dsq[tr->emitr[x]]];
	if      (mode == 3) esc = DegeneratePairScore(cm->abc, cm->esc[v], dsq[tr->emitl[x]], dsq[tr->emitr[x]]);
        else if (mode == 2) esc =   LeftMarginalScore(cm->abc, cm->esc[v], dsq[tr->emitl[x]]);
        else if (mode == 1) esc =  RightMarginalScore(cm->abc, cm->esc[v],                        dsq[tr->emitr[x]]);
      } else if ( (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) && (mode == 3 || mode == 2) ) {
	syml = cm->abc->sym[dsq[tr->emitl[x]]];
	esc  = esl_abc_FAvgScore(cm->abc, dsq[tr->emitl[x]], cm->esc[v]);
      } else if ( (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) && (mode == 3 || mode == 1) ) {
	symr = cm->abc->sym[dsq[tr->emitr[x]]];
	esc  = esl_abc_FAvgScore(cm->abc, dsq[tr->emitr[x]], cm->esc[v]);
      }

      /* Set tsc: transition score, or 0.
       * B, E, and the special EL state (M, local end) have no transitions.
       */
      tsc = 0.;
      if (v != cm->M && cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
	y = tr->state[tr->nxtl[x]];

        if (tr->nxtl[x] == -1)
          ;
	else if (v == 0 && (cm->flags & CM_LOCAL_BEGIN))
	  tsc = cm->beginsc[y];
	else if (y == cm->M) /* CM_LOCAL_END is presumably set, else this wouldn't happen */
	  tsc = cm->endsc[v] + (cm->el_selfsc * (tr->emitr[x] - tr->emitl[x] + 1 - StateDelta(cm->sttype[v])));
	else 		/* y - cm->first[v] gives us the offset in the transition vector */
	  tsc = cm->tsc[v][y - cm->cfirst[v]];
      }

      /* Print the info line for this state
       */
      if(do_banded)
	{
	  L = tr->emitr[x]-tr->emitl[x]+1;
	  fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5d %5d %5d %5.2f %5.2f %5d %5d %5d %2s\n",
		  x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], 
		  Statetype(cm->sttype[v]), tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc,
		  L, dmin[v], dmax[v],
		  (L >= dmin[v] && L <= dmax[v]) ? "" : "!!");
	}
      else
	{
	  fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5d %5d %5d %5.2f %5.2f\n",
		  x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], 
		  Statetype(cm->sttype[v]), tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc);
	}
    }
  if(do_banded)
    fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	    "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----",
	    "-----", "-----", "-----");
  else
    fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s\n",
	    "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----");
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


/*
 * Function : Parsetrees2Alignment()
 *
 * Purpose:   Creates a MSA from a set of parsetrees and a CM.
 *
 * 
 * Args:     cm         - the CM the CP9 was built from, needed to get emitmap,
 *                        so we know where to put EL transitions
 *           abc        - alphabet to use to create the return MSA
 *           sq         - sequences, must be digitized (we check for it)
 *           wgt        - weights for seqs, NULL for none
 *           nseq       - number of sequences
 *           tr         - array of tracebacks
 *           do_full    - TRUE to always include all match columns in alignment
 *           do_matchonly - TRUE to ONLY include match columns
 *           ret_msa    - MSA, alloc'ed/created here
 *
 * Return:   eslOK on succes, eslEMEM on memory error.
 *           MSA structure in ret_msa, caller responsible for freeing.
 *
 * Returns:   eslOK on success, eslEMEM on memory error, 
 *            Also ret_msa is filled with a new MSA.
 *
 */
int
Parsetrees2Alignment(CM_t *cm, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
		     Parsetree_t **tr, int nseq, int do_full, int do_matchonly, 
		     ESL_MSA **ret_msa)
{
  /* Contract check. We allow the caller to specify the alphabet they want the 
   * resulting MSA in, but it has to make sense (see next few lines). */
  if(cm->abc->type == eslRNA)
    { 
      if(abc->type != eslRNA && abc->type != eslDNA)
	esl_fatal("ERROR in Parsetrees2Alignment(), cm alphabet is RNA, but requested output alphabet is neither DNA nor RNA.");
    }
  else if(cm->abc->K != abc->K)
    esl_fatal("ERROR in Parsetrees2Alignment(), cm alphabet size is %d, but requested output alphabet size is %d.", cm->abc->K, abc->K);

  int          status;       /* easel status flag */
  ESL_MSA     *msa   = NULL; /* multiple sequence alignment */
  CMEmitMap_t *emap  = NULL; /* consensus emit map for the CM */
  int          i;            /* counter over traces */
  int          v, nd;        /* state, node indices */
  int          cpos;         /* counter over consensus positions (0)1..clen */
  int         *matuse= NULL; /* TRUE if we need a cpos in mult alignment */
  int         *iluse = NULL; /* # of IL insertions after a cpos for 1 trace */
  int         *eluse = NULL; /* # of EL insertions after a cpos for 1 trace */
  int         *iruse = NULL; /* # of IR insertions after a cpos for 1 trace */
  int         *maxil = NULL; /* max # of IL insertions after a cpos */
  int         *maxel = NULL; /* max # of EL insertions after a cpos */
  int         *maxir = NULL; /* max # of IR insertions after a cpos */
  int	      *matmap= NULL; /* apos corresponding to a cpos */
  int         *ilmap = NULL; /* first apos for an IL following a cpos */
  int         *elmap = NULL; /* first apos for an EL following a cpos */
  int         *irmap = NULL; /* first apos for an IR following a cpos */
  int          alen;	     /* length of msa in columns */
  int          apos;	     /* position in an aligned sequence in MSA */
  int          rpos;	     /* position in an unaligned sequence in dsq */
  int          tpos;         /* position in a parsetree */
  int          el_len;	     /* length of an EL insertion in residues */
  CMConsensus_t *con = NULL; /* consensus information for the CM */
  int          prvnd;	     /* keeps track of previous node for EL */
  int          nins;          /* insert counter used for splitting inserts */

  emap = CreateEmitMap(cm);

  ESL_ALLOC(matuse, sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(iluse,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(eluse,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(iruse,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxil,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxel,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxir,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(matmap, sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(ilmap,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(elmap,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(irmap,  sizeof(int)*(emap->clen+1));   
  
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

      ilmap[cpos] = alen; alen += maxil[cpos];
      elmap[cpos] = alen; alen += maxel[cpos];
      alen += maxir[cpos]; irmap[cpos] = alen-1; 
    }

  /* We're getting closer.
   * Now we can allocate for the MSA.
   */
  msa = esl_msa_Create(nseq, alen);
  msa->nseq = nseq;
  msa->alen = alen;
  msa->abc  = (ESL_ALPHABET *) abc;

  for (i = 0; i < nseq; i++)
    {
      /* Contract check */
      if(! (sq[i]->flags & eslSQ_DIGITAL))
	esl_fatal("ERROR in Parsetrees2Alignment(), sq %d is not digitized.\n", i);

      /* Initialize the aseq with all pads '.' (in insert cols) 
       * and deletes '-' (in match cols).
       */
      for (apos = 0; apos < alen; apos++)
	msa->aseq[i][apos] = '.';
      for (cpos = 0; cpos <= emap->clen; cpos++)
	if (matmap[cpos] != -1) msa->aseq[i][matmap[cpos]] = '-';
      msa->aseq[i][alen] = '\0';

      /* Traverse this guy's trace, and place all his
       * emitted residues.
       */
      for (cpos = 0; cpos <= emap->clen; cpos++)
	iluse[cpos] = iruse[cpos] = 0;

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
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];

	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];
	    break;
	    
	  case ML_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];
	    break;

	  case MR_st:
	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];
	    break;

	  case IL_st:
	    cpos = emap->lpos[nd];
	    apos = ilmap[cpos] + iluse[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
	    iluse[cpos]++;
	    break;

	  case EL_st: 
            /* we can assert eluse[cpos] always == 0 when we enter,
	     * because we can only have one EL insertion event per 
             * cpos. If we ever decide to regularize (split) insertions,
             * though, we'll want to calculate eluse in the rpos loop.
             */
	    cpos = emap->epos[nd]; 
	    apos = elmap[cpos]; 
	    for (rpos = tr[i]->emitl[tpos]; rpos <= tr[i]->emitr[tpos]; rpos++)
	      {
		msa->aseq[i][apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
		apos++;
	      }
	    break;

	  case IR_st: 
	    cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
	    apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
	    iruse[cpos]++;
	    break;

	  case D_st:
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATL_D) 
	      {
		cpos = emap->lpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATR_D) 
	      {
		cpos = emap->rpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    break;

	  } /* end of the switch statement */


	  prvnd = nd;
	} /* end traversal over trace i. */

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
	  rightjustify(msa->abc, msa->aseq[i], maxil[0]);
	  /* EL's are flush left, we want flush right I think these are impossible, but just in case... */
	  rightjustify(msa->abc, msa->aseq[i]+maxil[0], maxel[0]);
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
		  rightjustify(msa->abc, msa->aseq[i]+matmap[cpos]+1+nins, maxil[cpos]-nins);
		}
	      if(maxel[cpos] > 1) /* we're flush LEFT, want to split */
		{
		  apos = matmap[cpos]+1 + maxil[cpos];
		  for (nins = 0; islower((int) (msa->aseq[i][apos])); apos++)
		    nins++;
		  nins /= 2;		/* split the insertion in half */
		  rightjustify(msa->abc, msa->aseq[i]+matmap[cpos]+1+maxil[cpos]+nins, maxel[cpos]-nins);
		}
	      if(maxir[cpos] > 1) /* we're flush RIGHT, want to split */
		{
		  apos = matmap[cpos+1]-1;
		  for (nins = 0; islower((int) (msa->aseq[i][apos])); apos--)
		    nins++;
		  nins ++; nins /= 2;		/* split the insertion in half (++ makes it same behavior as IL/EL */
		  leftjustify(msa->abc, msa->aseq[i]+matmap[cpos]+1 + maxil[cpos] + maxel[cpos], maxir[cpos]-nins);
		}
	    }
	  /* Deal with inserts after final consensus position, IL's then EL's, then IR's
	   * IL's are flush left, we want flush left, do nothing 
	   * EL's are flush left, we want flush left, do nothing 
	   * IR's are flush right, we want flush left */
	  leftjustify(msa->abc, msa->aseq[i]+matmap[emap->clen]+1 + maxil[emap->clen] + maxel[emap->clen], maxir[emap->clen]);
	}
    } /* end loop over all parsetrees */


  /* Gee, wasn't that easy?
   * Add the rest of the ("optional") information to the MSA.
   */
  CreateCMConsensus(cm, abc, 3.0, 1.0, &con);

  /* "author" info */
  ESL_ALLOC(msa->au, sizeof(char) * (strlen(PACKAGE_VERSION)+10));
  sprintf(msa->au, "Infernal %s", PACKAGE_VERSION);

  for (i = 0; i < nseq; i++)
    {
      esl_strdup(sq[i]->name, -1, &(msa->sqname[i]));
      /* TODO: individual SS annotations
       */
      if (wgt == NULL) msa->wgt[i] = 1.0;
      else             msa->wgt[i] = wgt[i];
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
      if (maxil[cpos] > 0) 
	for (apos = ilmap[cpos]; apos < ilmap[cpos] + maxil[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
      if (maxel[cpos] > 0)
	for (apos = elmap[cpos]; apos < elmap[cpos] + maxel[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '~';
	    msa->rf[apos] = '~';
	  }
      if (maxir[cpos] > 0)	/* remember to write backwards */
	for (apos = irmap[cpos]; apos > irmap[cpos] - maxir[cpos]; apos--)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
    }
  msa->ss_cons[alen] = '\0';
  msa->rf[alen] = '\0';

  /* If we only want the match columns, shorten the alignment
   * by getting rid of the inserts. (Alternatively we could probably
   * simplify the building of the alignment, but all that pretty code
   * above already existed, so we do this post-msa-building shortening).
   */
  if(do_matchonly)
    {
      int *useme;
      ESL_ALLOC(useme, sizeof(int) * (msa->alen));
      esl_vec_ISet(useme, msa->alen, FALSE);
      for(cpos = 0; cpos <= emap->clen; cpos++)
	if(matmap[cpos] != -1) useme[matmap[cpos]] = TRUE;
      esl_msa_ColumnSubset(msa, useme);
      free(useme);
    }

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
  *ret_msa = msa;
  return eslOK;

 ERROR:
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
  ESL_DSQ symi, symj;		/* symbol indices for emissions, 0..Alphabet_iupac-1 */
  int mode;
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
  /* Contract check, CM must be LOCALLY configured, (could demand global, but
   * we assume we'll be calling this function serially for many parses and don't
   * want to need to switch CM back and forth from local/global */
  if((!(cm->flags & CM_LOCAL_BEGIN)) || (!(cm->flags & CM_LOCAL_END)))
    esl_fatal("ERROR in ParsetreeScore_Global2Local() CM is not in local mode.\n");
  if(dsq == NULL)
    esl_fatal("ERROR in ParsetreeScore_Global2Local(), dsq is NULL.\n");

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
	esl_fatal("ERROR in ParsetreeScore_Global2Local(), EL in parse, but it should be global!\n");
      if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) /* no scores in B,E */
	{
	  y = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */

	  if (y == cm->M) 
	    esl_fatal("ERROR in ParsetreeScore_Global2Local(), EL in parse, but it should be global!\n");
	  if (v == 0 && y > cm->cnum[0])
	    esl_fatal("ERROR in ParsetreeScore_Global2Local(), we did a local begin in the parse, but it should be global!\n");
	  /* for v == 0, we don't care that transition score has changed from global CM that
	   * was used to generate the parsetree, because the transition from root is not
	   * considered when we look for best local parse below. */

	  /* y - cm->first[v] gives us the offset in the transition vector */
	  tr_tsc[tidx] = cm->tsc[v][y - cm->cfirst[v]];
	
	  if (cm->sttype[v] == MP_st) 
	    {
	      symi = dsq[tr->emitl[tidx]];
	      symj = dsq[tr->emitr[tidx]];
	      if (mode == 3)
		{
		  if (symi < cm->abc->K && symj < cm->abc->K)
		    tr_esc[tidx] = cm->esc[v][(int) (symi*cm->abc->K+symj)];
		  else
		    tr_esc[tidx] = DegeneratePairScore(cm->abc, cm->esc[v], symi, symj);
		}
	      else if (mode == 2)
		tr_esc[tidx] = LeftMarginalScore(cm->abc, cm->esc[v], symi);
	      else if (mode == 1)
		tr_esc[tidx] = RightMarginalScore(cm->abc, cm->esc[v], symj);
	    } 
	  else if ( (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) && (mode == 3 || mode == 2) )
	    {
	      symi = dsq[tr->emitl[tidx]];
	      if (symi < cm->abc->K) tr_esc[tidx] = cm->esc[v][(int) symi];
	      else                   tr_esc[tidx] = esl_abc_FAvgScore(cm->abc, symi, cm->esc[v]);
	    } 
	  else if ( (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) && (mode == 3 || mode == 2) )
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
  esl_fatal("ERROR in ParsetreeScore_Global2Local()\n");
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
 * CM_t  *cm                - the CM, must have a valid cm->cp9
 * Parsetree_t *cm_tr       - valid parsetree to convert
 * cp9trace_s *ret_cp9_tr   - the CP9 trace to return, alloc'ed here
 */
int
Parsetree2CP9trace(CM_t *cm, Parsetree_t *tr, CP9trace_t **ret_cp9_tr)
{
  /* Check the contract */
  if(cm->cp9 == NULL || (!(cm->flags & CM_CP9)))
    esl_fatal("In Parsetree2CP9trace, cm->cp9 is not valid.\n");
  if(cm->cp9map == NULL)
    esl_fatal("In Parsetree2CP9trace, cm->cp9map is NULL.\n");

  int status;                /* Easel status                            */
  CP9trace_t *cp9_tr; /* the CP9 trace we're creating            */
  int **ks_ct = NULL; /* [0..2][0..cp9->M] number of times each state was used
		       * 1st D: 0 = MATCH, 1 = INSERT, 2 = DELETE */
  int  tidx;                 /* counter over parsetree nodes */
  int  v;                    /* CM state index */
  int  k, ks;                /* HMM nodes and state indices */
  int  i;                    /* generic counter */
  int  cp9_tr_size;          /* number of nodes we'll need for cp9_tr */
  int  lmost_k = cm->cp9->M + 1; /* left most HMM node visited in parse (often 1) */
  int  rmost_k = 0;          /* right most HMM node visited in parse (often M) */
  int  ip;
  int  ins_ct = 0;               /* total number of inserts */
  ESL_ALLOC(ks_ct,           sizeof(int *) * 3);
  for(ks = 0; ks < 3; ks++)
    {
      ESL_ALLOC(ks_ct[ks], sizeof(int) * (cm->cp9->M+1));
      esl_vec_ISet(ks_ct[ks], cm->cp9->M+1, 0);
    }

  /* Traverse parsetree, keeping track of implied HMM states used by each HMM node. */
  v = tr->state[0]; 
  if(v != 0) esl_fatal("ERROR in Parsetree2CP9Trace(), first Parsetree node not root.\n");
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
  rmost_k = cm->cp9->M; 
  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN)
    {
      while((ks_ct[HMMMATCH][lmost_k] + ks_ct[HMMINSERT][lmost_k] + ks_ct[HMMDELETE][lmost_k]) == 0)
	lmost_k++;
    }
  if(cm->cp9->flags & CPLAN9_LOCAL_END)
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


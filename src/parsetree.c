/* parsetree.c
 * cove 1.0: Mon May 17 09:38:14 1993
 * moved to cove 2.0, Mon Sep  6 13:34:55 1993
 * cove4: SRE 29 Feb 2000 [Seattle]
 * infernal: SRE, Fri Jul 28 08:55:47 2000 [StL]
 * CVS $Id$
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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "structs.h"
#include "funcs.h"
#include "squid.h"

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
CreateParsetree(void)
{
  Parsetree_t *new;

  new           = MallocOrDie (sizeof(Parsetree_t));
  new->memblock = 100;		/* allocation block size can be optimized here if you want. */
  new->nalloc   = new->memblock;
  new->emitl    = MallocOrDie(sizeof(int) * new->nalloc);
  new->emitr    = MallocOrDie(sizeof(int) * new->nalloc);
  new->state    = MallocOrDie(sizeof(int) * new->nalloc);
  new->nxtl     = MallocOrDie(sizeof(int) * new->nalloc);
  new->nxtr     = MallocOrDie(sizeof(int) * new->nalloc);
  new->prv      = MallocOrDie(sizeof(int) * new->nalloc);
  new->n = 0;
  return new;
}

/* Function: GrowParsetree()
 * Incept:   SRE 1 March 2000 [Seattle]
 * 
 * Purpose:  Increase the number of available nodes in a parse tree.
 */
void
GrowParsetree(Parsetree_t *tr)
{
  tr->nalloc += tr->memblock;
  tr->emitl = ReallocOrDie(tr->emitl, sizeof(int) * tr->nalloc);
  tr->emitr = ReallocOrDie(tr->emitr, sizeof(int) * tr->nalloc);
  tr->state = ReallocOrDie(tr->state, sizeof(int) * tr->nalloc);
  tr->nxtl  = ReallocOrDie(tr->nxtl,  sizeof(int) * tr->nalloc);
  tr->nxtr  = ReallocOrDie(tr->nxtr,  sizeof(int) * tr->nalloc);
  tr->prv   = ReallocOrDie(tr->prv,   sizeof(int) * tr->nalloc);
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
  free(tr->nxtl);
  free(tr->nxtr);
  free(tr->prv);
  free(tr);
}

/* Function: InsertTraceNode()
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
InsertTraceNode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state)
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
ParsetreeCount(CM_t *cm, Parsetree_t *tr, char *dsq, float wgt)
{
  int tidx;			/* counter through positions in the parsetree        */
  int v,z;			/* parent, child state index in CM                   */

		/* trivial preorder traverse, since we're already numbered that way */
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) /* no parameters estimated from B,E */
      {
	z = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */
			/* z - cm->first[v] gives us the offset in the transition vector */
	cm->t[v][z - cm->cfirst[v]] += wgt;

	if (cm->sttype[v] == MP_st) 
	  PairCount(cm->e[v], dsq[tr->emitl[tidx]], dsq[tr->emitr[tidx]], wgt);
	else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	  SingletCount(cm->e[v], dsq[tr->emitl[tidx]], wgt);
	else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	  SingletCount(cm->e[v], dsq[tr->emitr[tidx]], wgt);
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
ParsetreeScore(CM_t *cm, Parsetree_t *tr, char *dsq)
{
  int tidx;			/* counter through positions in the parsetree        */
  int v,y;			/* parent, child state index in CM                   */
  char symi, symj;		/* symbol indices for emissions, 0..Alphabet_iupac-1 */
  float sc;			/* the log-odds score of the parse tree */

		/* trivial preorder traverse, since we're already numbered that way */
  sc = 0.;
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    if (cm->sttype[v] != E_st && cm->sttype[v] != B_st) /* no scores in B,E */
      {
	y = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */
			/* y - cm->first[v] gives us the offset in the transition vector */
	sc += cm->tsc[v][y - cm->cfirst[v]];
	
	if (cm->sttype[v] == MP_st) 
	  {
	    symi = dsq[tr->emitl[tidx]];
	    symj = dsq[tr->emitr[tidx]];
	    if (symi < Alphabet_size && symj < Alphabet_size)
	      sc += cm->esc[v][(int) (symi*Alphabet_size+symj)];
	    else
	      sc += DegeneratePairScore(cm->esc[v], symi, symj);
	  } 
	else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	  {
	    symi = dsq[tr->emitl[tidx]];
	    if (symi < Alphabet_size) sc += cm->esc[v][(int) symi];
	    else                      sc += DegenerateSingletScore(cm->esc[v], symi);
	  } 
	else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	  {
	    symj = dsq[tr->emitr[tidx]];
	    if (symi < Alphabet_size) sc += cm->esc[v][(int) symj];
	    else                      sc += DegenerateSingletScore(cm->esc[v], symj);
	  }
      }
  }
  return sc;
}




/* Function: PrintParsetree()
 * Date:     SRE, Fri Jul 28 12:47:06 2000 [St. Louis]
 *
 * Purpose:  Debugging: show a tabular representation of a
 *           parsetree structure.
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
 * Args:    
 *
 * Returns:  
 *
 * Example:  
 */
void
ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, char *dsq)
{
  int   x;
  char  syml, symr;
  float tsc;
  float esc;
  int   v;

  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s\n",
	  " idx ","emitl", "emitr", "state", " nxtl", " nxtr", " prv ", " tsc ", " esc ");
  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s\n",
	 "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----");
  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      syml = symr = ' ';
      esc = tsc = 0.;
      if (cm->sttype[v] == MP_st) {
	syml = Alphabet[(int)dsq[tr->emitl[x]]]; 
	symr = Alphabet[(int)dsq[tr->emitr[x]]];
	esc  = DegeneratePairScore(cm->esc[v], dsq[tr->emitl[x]], dsq[tr->emitr[x]]);
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	syml = Alphabet[(int)dsq[tr->emitl[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitl[x]]);
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	symr = Alphabet[(int)dsq[tr->emitr[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitr[x]]);
      }
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	tsc = cm->tsc[v][tr->state[tr->nxtl[x]] - cm->cfirst[v]];

      fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5d %5d %5d %5.2f %5.2f\n",
	      x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], Statetype(cm->sttype[v]),
	      tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc);

    }
  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	  "-----", "-----", "-----", "-----", "-----","-----", "-----","-----", "-----");
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

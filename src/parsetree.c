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
 * The trace tree structure has a dummy node at its beginning,
 * and dummy end nodes at the termination of each branch. Non-BIFURC
 * states have a NULL right branch. 
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
 * Purpose:  Creates a parse tree tree structure.
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
  new->type     = MallocOrDie(sizeof(int) * new->nalloc);
  new->nxtl     = MallocOrDie(sizeof(int) * new->nalloc);
  new->nxtr     = MallocOrDie(sizeof(int) * new->nalloc);
  new->prv      = MallocOrDie(sizeof(int) * new->nalloc);
  
				/* state 0 is always the dummy beginning. */
  new->state[0] = 0;
  new->type[0]  = DUMMY;
  new->emitl[0] = -1;
  new->emitr[0] = -1;
  new->nxtl[0]  = 1;		/* connects to dummy end */
  new->nxtr[0]  = 1;		/* connects to dummy end */
  new->prv[0]   = -1;		/* no parent */

				/* state 1 is always the dummy end. */
  new->state[1] = 1;
  new->type[1]  = DUMMY;
  new->emitl[1] = -1;
  new->emitr[1] = -1;
  new->nxtl[1]  = -1;
  new->nxtr[1]  = -1;
  new->prv[1]   = -1;		/* end state is shared: has multiple parents, and forgets them. */

  new->n = 2;

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
  tr->type  = ReallocOrDie(tr->type,  sizeof(int) * tr->nalloc);
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
  free(tr->type);
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
 * Returns:  index of new node.
 */          
int
InsertTraceNode(Parsetree_t *tr, int y, int whichway,
		int emitl, int emitr, int state, int type)
{
  int a;
  int n;

  n = tr->n;
  a = (whichway == TRACE_LEFT_CHILD ? tr->nxtl[y] : tr->nxtr[y]);

  if (tr->n == tr->nalloc) GrowParsetree(tr);
				/* information in new node */
  tr->emitl[n] = emitl;
  tr->emitr[n] = emitr;
  tr->state[n] = state;
  tr->type[n]  = type;
				/* connectivity of new node */
  tr->nxtl[n]  = a;
  tr->nxtr[n]  = 1;
  tr->prv[n]   = y;
				/* connectivity of parent   */
  if (whichway == TRACE_LEFT_CHILD)  tr->nxtl[y] = n;
  else                               tr->nxtr[y] = n;
				/* connectivity of child    */
  if (a != 1)  tr->prv[a] = n;
				/* bump counter, return index of new node */
  tr->n++;
  return n;
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

  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s %5s\n",
	  " idx ","emitl", "emitr", "state", " type", " nxtl", " nxtr", " prv ");
  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s %5s\n",
	 "-----", "-----", "-----", "-----","-----","-----","-----", "-----");
  for (x = 0; x < tr->n; x++)
    fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %5d\n",
	   x, tr->emitl[x], tr->emitr[x], tr->state[x], tr->type[x],
	   tr->nxtl[x], tr->nxtr[x], tr->prv[x]);
  fprintf(fp, "%5s %5s %5s %5s %5s %5s %5s %5s\n",
	 "-----", "-----", "-----", "-----","-----","-----","-----", "-----");

  fprintf(fp, "n      = %d\n", tr->n);
  fprintf(fp, "nalloc = %d\n", tr->nalloc);
  fprintf(fp, "block  = %d\n", tr->memblock);
}

/* Function: SummarizeMasterTrace()
 * Date:     SRE, Fri Jul 28 13:42:30 2000 [St. Louis]
 *
 * Purpose:  Debugging: count the nodes used in a master trace.
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
				/* start at 2 to skip the dummies */
  for (x = 2; x < tr->n; x++)
    count[tr->type[x]]++;

  fprintf(fp, "Bifurcations: %d\n", count[0]);
  fprintf(fp, "MATP:         %d\n", count[1]);
  fprintf(fp, "MATL:         %d\n", count[2]);
  fprintf(fp, "MATR:         %d\n", count[3]);
  fprintf(fp, "BEGL:         %d\n", count[4]);
  fprintf(fp, "BEGR:         %d\n", count[5]);
  fprintf(fp, "ROOT:         %d\n", count[6]);
  fprintf(fp, "END:          %d\n", count[7]);
}

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


#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "squid.h"

#include "structs.h"
#include "funcs.h"
#include "cplan9.h"

#include "easel.h"
#include "esl_vectorops.h"

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
  new->mode     = MallocOrDie(sizeof(int) * new->nalloc);
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
  tr->mode  = ReallocOrDie(tr->mode,  sizeof(int) * tr->nalloc);
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
ParsetreeCount(CM_t *cm, Parsetree_t *tr, char *dsq, float wgt)
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
ParsetreeScore(CM_t *cm, Parsetree_t *tr, char *dsq, int do_null2)
{
  int tidx;			/* counter through positions in the parsetree        */
  int v,y;			/* parent, child state index in CM                   */
  char symi, symj;		/* symbol indices for emissions, 0..Alphabet_iupac-1 */
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

	if (v == 0 && (cm->flags & CM_LOCAL_BEGIN))
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
  	        if (symi < Alphabet_size && symj < Alphabet_size)
	          sc += cm->esc[v][(int) (symi*Alphabet_size+symj)];
	        else
	          sc += DegeneratePairScore(cm->esc[v], symi, symj);
              }
            else if (mode == 2)
              sc += LeftMarginalScore(cm->esc[v], symi);
            else if (mode == 1)
              sc += RightMarginalScore(cm->esc[v], symj);
	  } 
	else if ( (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) && (mode == 3 || mode == 2) )
	  {
	    symi = dsq[tr->emitl[tidx]];
	    if (symi < Alphabet_size) sc += cm->esc[v][(int) symi];
	    else                      sc += DegenerateSingletScore(cm->esc[v], symi);
	  } 
	else if ( (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) && (mode == 3 || mode == 2) )
	  {
	    symj = dsq[tr->emitr[tidx]];
	    if (symj < Alphabet_size) sc += cm->esc[v][(int) symj];
	    else                      sc += DegenerateSingletScore(cm->esc[v], symj);
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
  int   v,y;
  int   mode;

  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s\n",
	  " idx ","emitl", "emitr", "state", " nxtl", " nxtr", " prv ", " tsc ", " esc ");
  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s\n",
	 "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----");
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
	if (mode == 3 || mode == 2) syml = Alphabet[(int)dsq[tr->emitl[x]]]; 
	if (mode == 3 || mode == 1) symr = Alphabet[(int)dsq[tr->emitr[x]]];
	if      (mode == 3) esc = DegeneratePairScore(cm->esc[v], dsq[tr->emitl[x]], dsq[tr->emitr[x]]);
        else if (mode == 2) esc =   LeftMarginalScore(cm->esc[v], dsq[tr->emitl[x]]);
        else if (mode == 1) esc =  RightMarginalScore(cm->esc[v],                    dsq[tr->emitr[x]]);
      } else if ( (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) && (mode == 3 || mode == 2) ) {
	syml = Alphabet[(int)dsq[tr->emitl[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitl[x]]);
      } else if ( (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) && (mode == 3 || mode == 1) ) {
	symr = Alphabet[(int)dsq[tr->emitr[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitr[x]]);
      }

      /* Set tsc: transition score, or 0.
       * B, E, and the special EL state (M, local end) have no transitions.
       */
      tsc = 0.;
      if (v != cm->M && cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
	y = tr->state[tr->nxtl[x]];

	if (v == 0 && (cm->flags & CM_LOCAL_BEGIN))
	  tsc = cm->beginsc[y];
	else if (y == cm->M) /* CM_LOCAL_END is presumably set, else this wouldn't happen */
	  tsc = cm->endsc[v] + (cm->el_selfsc * (tr->emitr[x] - tr->emitl[x] + 1 - StateDelta(cm->sttype[v])));
	else 		/* y - cm->first[v] gives us the offset in the transition vector */
	  tsc = cm->tsc[v][y - cm->cfirst[v]];
      }

      /* Print the info line for this state
       */
      fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5d %5d %5d %5.2f %5.2f\n",
	      x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], 
	      Statetype(cm->sttype[v]), tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc);
    }
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


/***********************************************************************************
 * Function : Parsetrees2Alignment()
 *            
 *            EPN Modified to optionally print out all consensus columns (even those
 *            that have all gaps aligned to them) if the --full option was used.
 ***********************************************************************************/

MSA *
Parsetrees2Alignment(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
		     Parsetree_t **tr, int nseq, int do_full)
{
  MSA         *msa;          /* multiple sequence alignment */
  CMEmitMap_t *emap;         /* consensus emit map for the CM */
  int          i;            /* counter over traces */
  int          v, nd;        /* state, node indices */
  int          cpos;         /* counter over consensus positions (0)1..clen */
  int         *matuse;       /* TRUE if we need a cpos in mult alignment */
  int         *iluse;        /* # of IL insertions after a cpos for 1 trace */
  int         *eluse;        /* # of EL insertions after a cpos for 1 trace */
  int         *iruse;        /* # of IR insertions after a cpos for 1 trace */
  int         *maxil;        /* max # of IL insertions after a cpos */
  int         *maxel;        /* max # of EL insertions after a cpos */
  int         *maxir;        /* max # of IR insertions after a cpos */
  int	      *matmap;       /* apos corresponding to a cpos */
  int         *ilmap;        /* first apos for an IL following a cpos */
  int         *elmap;        /* first apos for an EL following a cpos */
  int         *irmap;        /* first apos for an IR following a cpos */
  int          alen;	     /* length of msa in columns */
  int          apos;	     /* position in an aligned sequence in MSA */
  int          rpos;	     /* position in an unaligned sequence in dsq */
  int          tpos;         /* position in a parsetree */
  int          el_len;	     /* length of an EL insertion in residues */
  CMConsensus_t *con;        /* consensus information for the CM */
  int          prvnd;	     /* keeps track of previous node for EL */

  emap = CreateEmitMap(cm);

  matuse = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  eluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iruse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxil  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxel  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxir  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  matmap = MallocOrDie(sizeof(int)*(emap->clen+1));   
  ilmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  elmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  irmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      /* EPN modified only following 4 lines to allow --full option */
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
  msa = MSAAlloc(nseq, alen);
  msa->nseq = nseq;
  msa->alen = alen;

  for (i = 0; i < nseq; i++)
    {
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
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];

	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;
	    
	  case ML_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;

	  case MR_st:
	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;

	  case IL_st:
	    cpos = emap->lpos[nd];
	    apos = ilmap[cpos] + iluse[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
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
		msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
		apos++;
	      }
	    break;

	  case IR_st: 
	    cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
	    apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
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

      /* Here is where we could put some insert-regularization code
       * a la HMMER: reach into each insert, find a random split point,
       * and shove part of it flush-right. But, for now, don't bother.
       */

    } /* end loop over all parsetrees */


  /* Gee, wasn't that easy?
   * Add the rest of the ("optional") information to the MSA.
   */
  con = CreateCMConsensus(cm, 3.0, 1.0);

  /* "author" info */
  msa->au   = MallocOrDie(sizeof(char) * (strlen(PACKAGE_VERSION)+10));
  sprintf(msa->au, "Infernal %s", PACKAGE_VERSION);

  for (i = 0; i < nseq; i++)
    {
      msa->sqname[i] = sre_strdup(sqinfo[i].name, -1);
      msa->sqlen[i]  = sqinfo[i].len;
      if (sqinfo[i].flags & SQINFO_ACC)
        MSASetSeqAccession(msa, i, sqinfo[i].acc);
      if (sqinfo[i].flags & SQINFO_DESC)
        MSASetSeqDescription(msa, i, sqinfo[i].desc);

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
  msa->ss_cons = MallocOrDie(sizeof(char) * (alen+1));
  msa->rf = MallocOrDie(sizeof(char) * (alen+1));
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
  return msa;
}

/***********************************************************************************
 * Function : ESL_Parsetrees2Alignment()
 *            
 *            EPN Modified to optionally print out all consensus columns (even those
 *            that have all gaps aligned to them) if the --full option was used.
 ***********************************************************************************/

MSA *
ESL_Parsetrees2Alignment(CM_t *cm, ESL_SQ **sq, float *wgt, 
			 Parsetree_t **tr, int nseq, int do_full)
{
  MSA         *msa;          /* multiple sequence alignment */
  CMEmitMap_t *emap;         /* consensus emit map for the CM */
  int          i;            /* counter over traces */
  int          v, nd;        /* state, node indices */
  int          cpos;         /* counter over consensus positions (0)1..clen */
  int         *matuse;       /* TRUE if we need a cpos in mult alignment */
  int         *iluse;        /* # of IL insertions after a cpos for 1 trace */
  int         *eluse;        /* # of EL insertions after a cpos for 1 trace */
  int         *iruse;        /* # of IR insertions after a cpos for 1 trace */
  int         *maxil;        /* max # of IL insertions after a cpos */
  int         *maxel;        /* max # of EL insertions after a cpos */
  int         *maxir;        /* max # of IR insertions after a cpos */
  int	      *matmap;       /* apos corresponding to a cpos */
  int         *ilmap;        /* first apos for an IL following a cpos */
  int         *elmap;        /* first apos for an EL following a cpos */
  int         *irmap;        /* first apos for an IR following a cpos */
  int          alen;	     /* length of msa in columns */
  int          apos;	     /* position in an aligned sequence in MSA */
  int          rpos;	     /* position in an unaligned sequence in dsq */
  int          tpos;         /* position in a parsetree */
  int          el_len;	     /* length of an EL insertion in residues */
  CMConsensus_t *con;        /* consensus information for the CM */
  int          prvnd;	     /* keeps track of previous node for EL */

  emap = CreateEmitMap(cm);

  matuse = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  eluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iruse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxil  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxel  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxir  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  matmap = MallocOrDie(sizeof(int)*(emap->clen+1));   
  ilmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  elmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  irmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      /* EPN modified only following 4 lines to allow --full option */
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
  msa = MSAAlloc(nseq, alen);
  msa->nseq = nseq;
  msa->alen = alen;

  for (i = 0; i < nseq; i++)
    {
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
	    msa->aseq[i][apos] = Alphabet[(int) sq[i]->dsq[rpos]];

	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) sq[i]->dsq[rpos]];
	    break;
	    
	  case ML_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) sq[i]->dsq[rpos]];
	    break;

	  case MR_st:
	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) sq[i]->dsq[rpos]];
	    break;

	  case IL_st:
	    cpos = emap->lpos[nd];
	    apos = ilmap[cpos] + iluse[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) sq[i]->dsq[rpos]]);
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
		msa->aseq[i][apos] = tolower((int) Alphabet[(int) sq[i]->dsq[rpos]]);
		apos++;
	      }
	    break;

	  case IR_st: 
	    cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
	    apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) sq[i]->dsq[rpos]]);
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

      /* Here is where we could put some insert-regularization code
       * a la HMMER: reach into each insert, find a random split point,
       * and shove part of it flush-right. But, for now, don't bother.
       */

    } /* end loop over all parsetrees */


  /* Gee, wasn't that easy?
   * Add the rest of the ("optional") information to the MSA.
   */
  con = CreateCMConsensus(cm, 3.0, 1.0);

  /* "author" info */
  msa->au   = MallocOrDie(sizeof(char) * (strlen(PACKAGE_VERSION)+10));
  sprintf(msa->au, "Infernal %s", PACKAGE_VERSION);

  for (i = 0; i < nseq; i++)
    {
      msa->sqname[i] = sre_strdup(sq[i]->name, -1);
      msa->sqlen[i]  = sq[i]->n;
      /*if (sqinfo[i].flags & SQINFO_ACC)
        MSASetSeqAccession(msa, i, sqinfo[i].acc);
      if (sqinfo[i].flags & SQINFO_DESC)
        MSASetSeqDescription(msa, i, sqinfo[i].desc);
      */

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
  msa->ss_cons = MallocOrDie(sizeof(char) * (alen+1));
  msa->rf = MallocOrDie(sizeof(char) * (alen+1));
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
  return msa;
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
ParsetreeScore_Global2Local(CM_t *cm, Parsetree_t *tr, char *dsq)
{
  int   status;
  int tidx;			/* counter through positions in the parsetree        */
  int v,y;			/* parent, child state index in CM                   */
  char symi, symj;		/* symbol indices for emissions, 0..Alphabet_iupac-1 */
  int mode;
  float *tr_esc;                /* [0..tr->n-1] score of emissions from each trace node */
  float *tr_tsc;                /* [0..tr->n-1] score of transitions from each trace node */
  float *v2n_map;               /* [0..cm->M-1], the trace node each state v corresponds to 
				 * -1 if none */
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

  /* Allocate and initialize */
  ESL_ALLOC(v2n_map, sizeof(int)   * cm->M); 
  ESL_ALLOC(lsc,     sizeof(float) * tr->n);
  ESL_ALLOC(tr_esc,  sizeof(float) * tr->n); 
  ESL_ALLOC(tr_tsc,  sizeof(float) * tr->n); 
  esl_vec_FSet(v2n_map, cm->M, -1);
  esl_vec_FSet(lsc, tr->n, 0.);
  esl_vec_FSet(tr_tsc, tr->n, 0.);
  esl_vec_FSet(tr_esc, tr->n, 0.);

  /* Determine the score that each trace node contributes to the overall parsetree score */

  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    v2n_map[v] = tidx;
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

	/* y - cm->first[v] gives us the offset in the transition vector */
	tr_tsc[tidx] = cm->tsc[v][y - cm->cfirst[v]];
	
	if (cm->sttype[v] == MP_st) 
	  {
	    symi = dsq[tr->emitl[tidx]];
	    symj = dsq[tr->emitr[tidx]];
            if (mode == 3)
              {
  	        if (symi < Alphabet_size && symj < Alphabet_size)
	          tr_esc[tidx] = cm->esc[v][(int) (symi*Alphabet_size+symj)];
	        else
	          tr_esc[tidx] = DegeneratePairScore(cm->esc[v], symi, symj);
              }
            else if (mode == 2)
              tr_esc[tidx] = LeftMarginalScore(cm->esc[v], symi);
            else if (mode == 1)
              tr_esc[tidx] = RightMarginalScore(cm->esc[v], symj);
	  } 
	else if ( (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) && (mode == 3 || mode == 2) )
	  {
	    symi = dsq[tr->emitl[tidx]];
	    if (symi < Alphabet_size) tr_esc[tidx] = cm->esc[v][(int) symi];
	    else                      tr_esc[tidx] = DegenerateSingletScore(cm->esc[v], symi);
	  } 
	else if ( (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) && (mode == 3 || mode == 2) )
	  {
	    symj = dsq[tr->emitr[tidx]];
	    if (symj < Alphabet_size) tr_esc[tidx] = cm->esc[v][(int) symj];
	    else                      tr_esc[tidx] = DegenerateSingletScore(cm->esc[v], symj);
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
      tidx = v2n_map[v];
      if(tidx == -1) continue; /* v's not in the parse */

      if(cm->sttype[v] == B_st)
	below_me_sc = lsc[tr->nxtl[tidx]] + lsc[tr->nxtr[tidx]];
      else if(cm->sttype[v] == E_st)
	below_me_sc = 0.;
      else
	below_me_sc = lsc[tr->nxtl[tidx]];

      /* Check if we could've jumped to an EL instead of traversing the 
       * subparse rooted here at v, would it have been worth it? */
      tmp_endsc = cm->endsc[v] + 
	cm->el_selfsc *  /* score of emitting 1 residue */
	((tr->emitr[tidx] - StateRightDelta(cm->sttype[v])) - 
	 (tr->emitl[tidx] + StateLeftDelta(cm->sttype[v])) + 1); /* number of residues EL must emit */
      if(below_me_sc < tmp_endsc -tr_tsc[tidx]) /* careful to consider sc of transition out of v */
	below_me_sc = tmp_endsc - tr_tsc[tidx]; /* we'll add tr_tsc[idx] back in next */
	
      lsc[tidx] = tr_esc[tidx] + tr_tsc[tidx] + below_me_sc; /* note we add in tsc even if local end taken */

      /*printf("tidx: %d\nv: %d\nlsc[tidx]: %f\nbegin_sc: %f\nmax_local_sc: %f\ntmp_endsc: %f\n\n", tidx, v, lsc[tidx], cm->beginsc[v], max_local_sc, tmp_endsc);
	printf("tr_esc[tidx]: %f\ntr_tsc[tidx]: %f\nbelow_me_sc: %f\n", tr_esc[tidx], tr_tsc[tidx], below_me_sc);*/

      /* Could we have jumped into this state from ROOT_S? Would it have
       * been worth it (based on what I've seen so far) */
      if((cm->beginsc[v] + lsc[tidx]) > max_local_sc)
	max_local_sc = cm->beginsc[v] + lsc[tidx];
    }

  free(v2n_map);
  free(lsc);
  free(tr_esc);
  free(tr_tsc);
  printf("in ParsetreeScore_Global2Local() returning sc: %f\n", max_local_sc);
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
Parsetree2CP9trace(CM_t *cm, Parsetree_t *tr, struct cp9trace_s **ret_cp9_tr)
{
  /* Check the contract */
  if(cm->cp9 == NULL || (!(cm->flags & CM_CP9)))
    esl_fatal("In Parsetree2CP9trace, cm->cp9 is not valid.\n");
  if(cm->cp9map == NULL)
    esl_fatal("In Parsetree2CP9trace, cm->cp9map is NULL.\n");

  int status;                /* Easel status                            */
  struct cp9trace_s *cp9_tr; /* the CP9 trace we're creating            */
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
  return -1;
}

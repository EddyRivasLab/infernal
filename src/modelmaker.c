/* modelmaker.c
 * SRE, 28 Feb 2000
 * SVN $Id$
 *
 * Construct a model from an alignment. 
 *
 * Outline of the process:
 *    1. construct a "guide" tree (gtr) for the alignment, 
 *       specifying which columns are match vs. insert and
 *       how the model tree branches. 
 *    2. The guide tree is converted to a CM by cm_from_guide().
 *    3. Individual tracebacks are constructed from individual 
 *       aligned sequences by transmogrify().
 *    4. The individual tracebacks are counted into a new model 
 *       with ParsetreeCount().
 *
 * The CM containing counts is returned. The caller has to assign a 
 * prior to it, and convert it to probabilities; then assign a null 
 * model to it and convert to log-odds scores.  
 *
 * The "guide tree" is a special use of a Parsetree_t structure. 
 * - tr->state contains a node type (e.g. MATP_nd), not a state index.
 * - The numbering of the guide tree is a preorder traverse, identical to
 *   the numbering in the final CM.
 * - emitl and emitr are relative to the alignment columns, not individual
 *   sequence positions.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */


#include <stdlib.h>

#include "squid.h"		
#include "msa.h"		/* multiple sequence alignments */
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"




static void         cm_from_guide(CM_t *cm, Parsetree_t *gtr);

/* Function: HandModelmaker()
 * Incept:   SRE 29 Feb 2000 [Seattle]; from COVE 2.0 code
 * 
 * Purpose:  The customer always knows best.
 * 
 *           Construct a model given a stated structure. The structure
 *           is provided via a "ss_cons" (consensus structure) line, as would
 *           occur in an annotated SELEX or Stockholm file. Only > and < characters
 *           in this line are interpreted (as base pairs). Pseudoknots, 
 *           if annotated, are ignored.
 *           
 *           Match vs. insert can be determined one of two ways. By default,
 *           the assignment is made by "gapthresh"; for columns with
 *           fractional occurence of gaps greater than this, the column
 *           is assigned to insert. If "use_rf" is TRUE, the rf (reference)
 *           line is interpreted as the assignment -- columns with non-space
 *           characters in the rf line are assigned to MATCH.
 *           
 *           Both rf and cs are provided in the msa structure.
 *           
 * Args:     msa       - multiple alignment to build model from
 *           dsq       - digitized aligned sequences            
 *           use_rf    - TRUE to use RF annotation to determine match/insert
 *           gapthresh - fraction of gaps to allow in a match column (if use_rf=FALSE)
 *           ret_cm    - RETURN: new model                      (maybe NULL)
 *           ret_gtr   - RETURN: guide tree for alignment (maybe NULL)
 *           
 * Return:   void
 *           cm is allocated here. FreeCM(*ret_cm).
 *           gtr is allocated here. FreeTrace().
 */
void
HandModelmaker(MSA *msa, char **dsq, int use_rf, float gapthresh, 
	       CM_t **ret_cm, Parsetree_t **ret_gtr)
{
  CM_t           *cm;		/* new covariance model                       */
  Parsetree_t    *gtr;		/* guide tree for alignment                   */
  Nstack_t       *pda;		/* pushdown stack used in building gtr        */
  int            *matassign;	/* 0..alen-1 array; 0=insert col, 1=match col */
  int            *ct;		/* 0..alen-1 base pair partners array         */
  int             apos;		/* counter over columns of alignment          */
  int             idx;		/* counter over sequences in the alignment    */
  int             v;		/* index of current node                      */
  int             i,j,k;	/* subsequence indices                        */
  int  type;			/* type of node we're working on              */
  int  diff, bestdiff, bestk;   /* used while finding optimal split points    */   
  int  nnodes;			/* number of nodes in CM                      */
  int  nstates;			/* number of states in CM                     */

  if (msa->ss_cons == NULL)
    Die("No consensus structure annotation available for that alignment.");
  if (use_rf && msa->rf == NULL) 
    Die("No reference annotation available for that alignment.");

  /* 1. Determine match/insert assignments
   *    matassign is 1..alen. Values are 1 if a match column, 0 if insert column.
   */
  matassign = MallocOrDie(sizeof(int) * (msa->alen+1));
  if (use_rf)
    {
      for (apos = 1; apos <= msa->alen; apos++)
	matassign[apos] = (isgap(msa->rf[apos-1]) ? 0 : 1);
    }
  else
    {
      int gaps;
      for (apos = 1; apos <= msa->alen; apos++)
	{
	  for (gaps = 0, idx = 0; idx < msa->nseq; idx++)
	    if (dsq[idx][apos] == DIGITAL_GAP) gaps++;
	  matassign[apos] = ((double) gaps / (double) msa->nseq > gapthresh) ? 0 : 1;
	}
    }

  /* 2. Determine a "ct" array, base-pairing partners for each position.
   *    Disallow/ignore pseudoknots. (That's what the FALSE flag does.)
   *    ct[] values give the index of a base pairing partner, or 0 for unpaired positions.
   *    Even though msa->ss_cons is in the 0..alen-1 coord system of msa, ct[]
   *    comes back in the 1..alen coord system of dsq.
   */
  if (! WUSS2ct(msa->ss_cons, msa->alen, FALSE, &ct))  
    Die("Consensus structure string is inconsistent"); 

  /* 3. Make sure the consensus structure "ct" is consistent with the match assignments.
   *    Wipe out all structure in insert columns; including the base-paired 
   *    partner of insert-assigned columns.
   */
  for (apos = 1; apos <= msa->alen; apos++)
    if (! matassign[apos])
      { 
	if (ct[apos] != 0)  ct[ct[apos]] = 0;
	ct[apos] = 0;
      }

  /* 4. Construct a guide tree.
   *    This code is borrowed from yarn's KHS2Trace().
   *    
   *    We also keep track of how many states we'll need in the final CM,
   *    so we'll know how much to allocate -- and the number of nodes,
   *    for informational purposes.
   */
  nstates = nnodes = 0;
  gtr = CreateParsetree();	/* the parse tree we'll grow        */
  pda = CreateNstack();		/* a pushdown stack for our indices */

  /* Construction strategy has to make sure we number the nodes in
   * preorder traversal: for bifurcations, we can't attach the right 
   * child until we've fully traversed the left side. Therefore, we have
   * to push what we intend to attach, and pop it later. And since we
   * don't know an index for the node until we attach it, we have no
   * place to put the node's data except the stack -- so we have to
   * push several numbers onto the stack: what type of node, what
   * subseq it's responsible for (emitl...emitr), and what node
   * index it attaches to.
   * 
   * Note that we have to deal with the fact that ct is off-by-one
   * in both indices and values: e.g. the base pairing partner 
   * j of residue i is ct[i-1]-1. 
   */
  PushNstack(pda, -1);		/* what node it's attached to */
  PushNstack(pda, 1);		/* emitl */
  PushNstack(pda, msa->alen);	/* emitr */
  PushNstack(pda, ROOT_nd);	/* "state" (e.g. node type) */

  while (PopNstack(pda, &type))	/* pop a node type to attach */
    {
      PopNstack(pda, &j);
      PopNstack(pda, &i);	/* i..j == subseq we're responsible for */
      PopNstack(pda, &v);	/* v = index of parent node in gtr */

      /* This node accounts for i..j, but we usually don't know how yet.
       * Six possibilities:
       *    i > j; this is an END state; do nothing.
       *    this is already assigned as a BEGIN; push i,j
       *    i is unpaired; this is a MATL state; push i+1, j
       *    j is unpaired; this is a MATR state; push i,j-1
       *    i,j pair to each other; this is a MATP state; push i+1,j-1
       *    i,j pair but not to each other; this is a BIFURC state;
       *        pick mid ip <= mid < jp; push BEGIN i,mid and working i,mid,
       *        and push BEGIN mid+1,j and working mid+1,j
       */
      if (i > j) {
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, END_nd);
	nstates += 1;		/* END_nd -> E_st */
	nnodes++;
      }

      else if (type == ROOT_nd) { /* try to push i,j; but deal with INSL and INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, ROOT_nd);
	for (; i <= j; i++) if (matassign[i]) break;
	for (; j >= i; j--) if (matassign[j]) break;
	PushNstack(pda, v);	/* here v==0 always. */
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* ROOT_nd -> S_st, IL_st, IR_st */
	nnodes++;
      }

      else if (type == BEGL_nd) {    /* no inserts */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BEGL_nd);
	PushNstack(pda, v);	
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 1;		/* BEGL_nd -> S_st */
	nnodes++;
      }

      else if (type == BEGR_nd)  { /* look for INSL */
	v = InsertTraceNode(gtr, v, TRACE_RIGHT_CHILD, i, j, BEGR_nd);
	for (; i <= j; i++) if (matassign[i]) break; 
	PushNstack(pda, v);	
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 2;		/* BEGR_nd -> S_st IL_st */
	nnodes++;
      }

      else if (ct[i] == 0) {
	 	/* i unpaired. This is a MATL node; allow INSL */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
	for (i = i+1; i <= j; i++)  if (matassign[i]) break;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* MATL_nd -> ML_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[j] == 0) { 	/* j unpaired. MATR node. Deal with INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATR_nd);
	for (j = j-1; j >= i; j--) if (matassign[j]) break;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* MATR_nd -> MR_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[i] == j) { /* i,j paired to each other. MATP. deal with INSL, INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATP_nd);
	for (i = i+1; i <= j; i++) if (matassign[i]) break;
	for (j = j-1; j >= i; j--) if (matassign[j]) break;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 6;		/* MATP_nd -> MP_st, ML_st, MR_st, D_st, IL_st, IR_st */
	nnodes++;
      }

      else /* i,j paired but not to each other. BIFURC. no INS. */
	{
	  /* Here's the first of two places where we can optimize the topology 
           * of a CM. (The other comes from choosing a state traversal order when 
           * building the CM.) Imagine a multifurcation of four domains:
           *   [1]..[2]..[3]..[4]
           * The "default leftwise" rule means that BEGL/INSL generates the 
           * intervening sequences, so we must model the four stems as:
           *   [1],    ..[2],    ..[3],   ..[4]
           * but we have a choice of how we bifurcate:
           *   (1,2)(3,4)    1,(2,(3,4))   ((1,2),3),4 
           * Our choice affects the time and memory requirements of a divide
           * conquer alignment algorithm. (1,(2,(3,4)) is most efficient for
           * memory; (1,2)(3,4) is most efficient for time.
           * 
           * So we may want to choose carefully from several possible split 
           * points k (3, in the above example). A priori we only know one 
           * possible midpoint precisely: ct[i]+1, the next base after closing 
           * domain 1. We can find the others by scanning for them, and we can 
           * be reasonably efficient about scanning by using ct[] to instantly 
           * skip subdomains.
	   */
	  /* One possible rule: optimize by finding most balanced split.
           * Each stop of the following loop gives a possible midpoint k, which is
           * then evaluated, keeping track of the best split so far.
           */
	  v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BIF_nd);

	  bestk    = ct[i]+1;
	  bestdiff = msa->alen;
	  for (k = ct[i] + 1; k < ct[j]; k = ct[k] + 1) 
	    {
	      diff = abs(i+j-2*k); /* = len2-len1-1, where len2 = j-k+1, len1= k-i */
	      if (diff < bestdiff) {
		bestdiff = diff; 
		bestk    = k;
	      }
	      while (ct[k] == 0) k++;
	    }
				/* push the right BEGIN node first */
	  PushNstack(pda, v);	
	  PushNstack(pda, bestk);
	  PushNstack(pda, j);
	  PushNstack(pda, BEGR_nd);
				/* then push the left BEGIN node */
	  PushNstack(pda, v);	
	  PushNstack(pda, i);
	  PushNstack(pda, bestk-1);
	  PushNstack(pda, BEGL_nd);
	  nstates += 1;		/* BIF_nd -> B_st */
	  nnodes++;
	}

    }	/* while something's on the stack */
  FreeNstack(pda);
  free(ct);

  /* OK, we've converted ct into gtr -- gtr is a tree structure telling us the
   * arrangement of consensus nodes. Now do the drill for constructing a full model 
   * using this guide tree.
   */
  cm = CreateCM(nnodes, nstates);
  cm_from_guide(cm, gtr);
  CMZero(cm);

  free(matassign);
  if (ret_cm  != NULL) *ret_cm  = cm;  else FreeCM(cm);
  if (ret_gtr != NULL) *ret_gtr = gtr; else FreeParsetree(gtr);
}



/* Function: cm_from_guide()
 * Date:     SRE, Sat Jul 29 09:25:49 2000 [St. Louis]
 *
 * Purpose:  given a guide tree and an allocated CM, 
 *           fill in all the structural information of the CM.
 *           
 * Args:     cm  - allocated cm to construct
 *           gtr - guide tree
 *
 * Returns:  (void)
 */
static void
cm_from_guide(CM_t *cm, Parsetree_t *gtr)
{
  Nstack_t   *pda;              /* pushdown stack used for traversing gtr */
  int         v;		/* what node we're working on (in gtr index system)*/
  int         node;		/* what node (preorder traversal numbering of CM) */
  int         state;		/* what state (preorder traversal numbering of CM) */
  int  nxtnodetype;		/* type of a child node (e.g. MATP_nd) */
  int  prvnodetype;		/* type of a parent node (e.g. MATP_nd) */

  /* Some CM structural configuration info:
   * child_count[] gives how many states are connectable in a child node. 
   * parent_count[] gives how many states are connectable in a parent node.
   */
 				/* BIF, MATP, MATL, MATR, BEGL, BEGR, ROOT, END */  
  int child_count[] =             {  1,    4,    2,    2,    1,    1,    0,   1};
  int parent_count[] =            {  1,    6,    3,    3,    1,    2,    3,   0};

  node = state = 0;
  pda = CreateNstack();
  PushNstack(pda, 0);		/* push ROOT_nd onto the stack */
  while (PopNstack(pda, &v))
    {
      if      (gtr->state[v] == BIF_nd) {
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = BIF_nd;

	cm->sttype[state] = B_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = BIF_B;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = -1; /* we fill this in later, when we see the BEGR... */
	PushNstack(pda, state);	/* ... the trick we use to remember the connection */
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;
	
	node++;
	PushNstack(pda, gtr->nxtr[v]);
	PushNstack(pda, gtr->nxtl[v]);
      }

      else if (gtr->state[v] == MATP_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = MATP_nd;

	cm->sttype[state] = MP_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATP_MP;
	cm->cfirst[state] = state+4;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = ML_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATP_ML;
	cm->cfirst[state] = state+3;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-2;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = MR_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATP_MR;
	cm->cfirst[state] = state+2;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-3;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = D_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATP_D;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-4;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATP_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 5;
	state++;

	cm->sttype[state] = IR_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATP_IR;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 6;
	state++;

	node++;
	PushNstack(pda, gtr->nxtl[v]);
      }

      else if (gtr->state[v] == MATL_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = MATL_nd;

	cm->sttype[state] = ML_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATL_ML;
	cm->cfirst[state] = state+2;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = D_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATL_D;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-2;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATL_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 3;
	state++;

	node++;
	PushNstack(pda, gtr->nxtl[v]);
      }
      
      else if (gtr->state[v] == MATR_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = MATR_nd;

	cm->sttype[state] = MR_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATR_MR;
	cm->cfirst[state] = state+2;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = D_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATR_D;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-2;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = IR_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = MATR_IR;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 3;
	state++;

	node++;
	PushNstack(pda, gtr->nxtl[v]);
      }

      else if (gtr->state[v] == BEGL_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node]  = BEGL_nd;

	cm->sttype[state] = S_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = BEGL_S;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = 1;
	state++;

	node++;
	PushNstack(pda, gtr->nxtl[v]);
      }

      else if (gtr->state[v] == BEGR_nd) {
	int bifparent;

	nxtnodetype = gtr->state[gtr->nxtl[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node]  = BEGR_nd;

	/* A trick: we need to attach this start state to the previous
	 * bifurcation. We stored the bif state index by pushing it onto
	 * the pda -- retrieve it now.
	 */
	PopNstack(pda, &bifparent);
	cm->cnum[bifparent] = state; /* remember, cnum overloaded for bif: idx of right child */

	cm->sttype[state] = S_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = BEGR_S;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = bifparent;
	cm->pnum[state]   = 1;
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = BEGR_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 2;
	state++;

	node++;
	PushNstack(pda, gtr->nxtl[v]);
      }

      else if (gtr->state[v] == ROOT_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node]  = ROOT_nd;

	cm->sttype[state] = S_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = ROOT_S;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 2 + child_count[nxtnodetype]; 
	cm->plast[state] = -1;
	cm->pnum[state]   = 0;
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = ROOT_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 2 + child_count[nxtnodetype]; 
	cm->plast[state] = state;
	cm->pnum[state]   = 2;
	state++;

	cm->sttype[state] = IR_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = ROOT_IR;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype]; 
	cm->plast[state] = state;
	cm->pnum[state]   = 3;
	state++;

	node++;
	PushNstack(pda, gtr->nxtl[v]);
      }

      else if (gtr->state[v] == END_nd) {
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node]  = END_nd;

	cm->sttype[state] = E_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = END_E;
	cm->cfirst[state] = -1;
	cm->cnum[state]   = 0;
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	node++;
      }
    }
  FreeNstack(pda);
  cm->M     = state;
  cm->nodes = node;
}




/* Function: Transmogrify()
 * Date:     SRE, Thu May 30 15:10:22 2002 [a coffee shop in Madison]
 *
 * Purpose:  Construct a "fake" tree for a given aligned sequence (aseq)
 *           and its digitized form dsq, given a new CM structure (cm) and
 *           a model guide tree (gtr). Same as transmogrify(), above,
 *           except this version can deal w/ local alignments.
 *           
 *           Keep in mind that aseq is 0..alen-1, dsq is 1..alen, and 
 *           gtr is working in dsq's coordinates - hence the -1's in aseq
 *           indexing.
 *           
 *           We need aseq because we encode local alignment there: all
 *           non-insert columns marked '~' are local deletions. (Marking
 *           of gaps in insert columns is ignored.)
 *           
 *           Assumes that local begins do not go to insert states, and local
 *           ends do not come from insert states. To assert otherwise
 *           is an invalid input.
 *           
 *           It expects that local alignment transmogrification is only
 *           done in two situations: debugging, and training. Therefore,
 *           users don't provide local alignments to this function, and
 *           it's polite to Die() on any kind of input error.
 * 
 * Args:     cm    - the newly built covariance model, corresponding to gtr
 *           gtr   - guide tree
 *           dsq   - a digitized aligned sequence [1..alen]
 *           aseq  - aligned sequence itself [0..alen-1]
 *           alen  - length of alignment
 *
 * Returns:  the individual parse tree. 
 *           Caller is responsible for free'ing this w/ FreeParsetree().
 */
Parsetree_t *
Transmogrify(CM_t *cm, Parsetree_t *gtr, char *dsq, char *aseq, int alen)
{
  Parsetree_t *tr;
  int          node;		/* index of node in *gtr* we're working on */
  int          state;		/* index of a state in the *CM*            */
  int          type;		/* a unique statetype                      */
  Nstack_t    *pda;             /* pushdown automaton for positions in tr  */
  int          tidx;		/* index *in parsetree tr* of state        */
  int          i,j;		/* coords in aseq                          */
  int          started;		/* TRUE if we've transited out of ROOT     */
  int          ended;		/* TRUE if we've transited to EL and ended */
  int          nstarts;         /* # of local transits out of ROOT: <= 1   */
  int         *localrun;        /* local alignment gap run lengths         */
  int          need_leftside;
  int          need_rightside;

  tr  = CreateParsetree();
  pda = CreateNstack();
  
  started = FALSE;
  ended   = FALSE;
  nstarts = 0;

  /* We preprocess the aseq to help with local alignment.
   */
  localrun = MallocOrDie(sizeof(int) * (alen+1));
  localrun[0] = 0;
  for (i = 0; i <= alen; i++)
    if (i > 0 && aseq[i-1] == '~') localrun[i] = localrun[i-1]+1;
    else                           localrun[i] = 0;

  /* Because the gtr is already indexed in a preorder traversal,
   * we can preorder traverse it easily w/ a for loop...
   */
  tidx = -1;	   /* first state to attach to; -1=special case for attaching root */
  for (node = 0; node < cm->nodes; node++)
    {
      /* A generic sanity check: we can't end if we haven't started.
       */
      if (ended && ! started) goto FAILURE;

      /* A (big) switch on node type.
       */
      switch (gtr->state[node]) { 

	/* The root node.
	 * Assume ROOT_S=0, ROOT_IL=1, ROOT_IR=2.
	 * started/ended are always FALSE when we get here.
	 * we can init a local start (in IL or IR).
	 */
      case ROOT_nd:
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
			       gtr->emitl[node], gtr->emitr[node], 0);
	for (i = gtr->emitl[node]; i < gtr->emitl[gtr->nxtl[node]]; i++)
	  if (dsq[i] != DIGITAL_GAP) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   i, gtr->emitr[node], 1);
	    if (! started) { started = TRUE; nstarts++; }
	  }
	for (j = gtr->emitr[node]; j > gtr->emitr[gtr->nxtl[node]]; j--)
	  if (dsq[j] != DIGITAL_GAP) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, j, 2);	
	    if (! started) { started = TRUE; nstarts++; }
	  }
	break;

	/* A bifurcation node.
	 * Assume that we'll process the BEGL node next; push info
	 * for BEGR onto the PDA.
	 * If we ended above here, the B doesn't go into the parsetree.
	 * If we didn't start yet, the B doesn't go into the parsetree.
	 */
      case BIF_nd:
	if (ended) {
	  if (aseq[gtr->emitl[node]-1] == '~' && aseq[gtr->emitr[node]-1] == '~') 
	    break;
	  else 
	    goto FAILURE;
	}

	i = gtr->emitl[gtr->nxtl[node]];
	j = gtr->emitr[gtr->nxtl[node]];
	need_leftside = (localrun[j] - localrun[i-1] != j - i + 1) ? TRUE : FALSE; 

	i = gtr->emitl[gtr->nxtr[node]];
	j = gtr->emitr[gtr->nxtr[node]];
	need_rightside = (localrun[j] - localrun[i-1] != j - i + 1) ? TRUE : FALSE; 
	  
	if (need_leftside && need_rightside) {
	  state = CalculateStateIndex(cm, node, BIF_B);
	  tidx  = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				  gtr->emitl[node], gtr->emitr[node], state);
	  if (! started) { started = TRUE; nstarts++; }
	} 
	PushNstack(pda, ended);   /* remember our ending status */
	PushNstack(pda, started); /* remember our start status */
	PushNstack(pda, tidx);    /* remember index in tr; we pop in BEGR */
	break;

	/* A MATP node.
	 * If we see *,* in the seq, this is a local deletion.
         *    If we haven't started yet, just skip the node.
	 *    If we have ended already, just skip the node; 
         *    If we haven't ended yet, end on an EL.
         * (* in only one position is invalid input.)
         * Else, this is a real state: emission or deletion.
         *    If we thought we ended, that's invalid input.
         *    Else, attach this guy. If it's a new start, it gets attached
         *      to ROOT, and we bump nstarts; we should only do this once
         *      on valid input. 
	 */
      case MATP_nd:
	if (dsq[gtr->emitl[node]] == DIGITAL_GAP) {
	  if (dsq[gtr->emitr[node]] == DIGITAL_GAP) type = MATP_D;
	  else                                      type = MATP_MR;
	} else {
	  if (dsq[gtr->emitr[node]] == DIGITAL_GAP) type = MATP_ML;
	  else                                      type = MATP_MP;
	}

	if (type == MATP_D 
	    && aseq[gtr->emitl[node]-1] == '~' 
	    && aseq[gtr->emitr[node]-1] == '~')
	  {
	    if (! started || ended)  break;
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   gtr->emitl[node], gtr->emitr[node], cm->M);
	    ended = TRUE;
	    break;
	  }
	if (aseq[gtr->emitl[node]-1] == '~' || aseq[gtr->emitr[node]-1] == '~')
	  goto FAILURE;

	if (ended) goto FAILURE;
	state = CalculateStateIndex(cm, node, type);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
			       gtr->emitl[node], gtr->emitr[node], state);	      
	if (! started) { started = TRUE; nstarts++; }

	state = CalculateStateIndex(cm, node, MATP_IL);
	for (i = gtr->emitl[node]+1; i < gtr->emitl[gtr->nxtl[node]]; i++)
	  if (dsq[i] != DIGITAL_GAP) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   i, gtr->emitr[node]-1, state);
	    if (! started) goto FAILURE;
	  }

	state = CalculateStateIndex(cm, node, MATP_IR);
	for (j = gtr->emitr[node]-1; j > gtr->emitr[gtr->nxtl[node]]; j--)
	  if (dsq[j] != DIGITAL_GAP) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, j, state);	
	    if (! started) goto FAILURE;
	  }
	break;

	/* A MATL node.
	 * If we see * in the seq, this is a local deletion.
	 *   If we haven't started yet, skip the node.
	 *   If we have ended already, skip the node.
	 *   If we haven't ended yet, end on an EL.
	 * Else, this is a real state (emission or deletion).
	 *   If we thought we ended, this is invalid input.
	 *   Else, attach this guy. If it's a new start, it is
	 *   attached to root (tidx == -1), and we bump nstarts.
	 */
      case MATL_nd:
	if (dsq[gtr->emitl[node]] == DIGITAL_GAP) type = MATL_D;
	else                                      type = MATL_ML;

	if (type == MATL_D && aseq[gtr->emitl[node]-1] == '~')
	  {
	    if (! started || ended) break;
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   gtr->emitl[node], gtr->emitr[node], cm->M);
	    ended = TRUE;
	    break;
	  }

	if (ended) goto FAILURE;
	state = CalculateStateIndex(cm, node, type);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
			       gtr->emitl[node], gtr->emitr[node], state);
	if (! started) { started = TRUE; nstarts++; }

	state = CalculateStateIndex(cm, node, MATL_IL);
	for (i = gtr->emitl[node]+1; i < gtr->emitl[gtr->nxtl[node]]; i++)
	  if (dsq[i] != DIGITAL_GAP) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   i, gtr->emitr[node], state);
	    if (! started) goto FAILURE;
	  }
	break;

	/* MATR node. 
	 * Similar logic as MATL above.
	 */
      case MATR_nd:
	if (dsq[gtr->emitr[node]] == DIGITAL_GAP) type = MATR_D;
	else                                      type = MATR_MR;

	if (type == MATR_D && aseq[gtr->emitl[node]-1] == '~')
	  {
	    if (! started || ended)  break;
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   gtr->emitl[node], gtr->emitr[node], cm->M);
	    ended = TRUE;
	    break;
	  }

	if (ended) goto FAILURE;	
	state = CalculateStateIndex(cm, node, type);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
			       gtr->emitl[node], gtr->emitr[node], state);
	if (! started) { started = TRUE; nstarts++; }

	state = CalculateStateIndex(cm, node, MATR_IR);
	for (j = gtr->emitr[node]-1; j > gtr->emitr[gtr->nxtl[node]]; j--)
	  if (dsq[j] != DIGITAL_GAP) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   gtr->emitl[node], j, state);
	    if (! started) goto FAILURE;
	  }
	break;

	/* BEGL_nd. 
	 * If not started, or ended, skip node. Else, attach it.
	 */
      case BEGL_nd:
	if (! started || ended) break;
	state = CalculateStateIndex(cm, node, BEGL_S);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
			       gtr->emitl[node], gtr->emitr[node], state);
	break;

	/* BEGR_nd.
	 * Pop off info on whether we started or ended above this
	 *    node in the CM.
	 * Logic different than BEGL above, because BEGR is dealing
	 * with an insert left state:
	 * If we've started, and not ended, attach the node.
	 * In dealing with inserts, if we think we've ended already,
	 * that's an invalid input.
	 */
      case BEGR_nd:
	PopNstack(pda, &tidx);	  /* recover parent bifurcation's index in trace */
	PopNstack(pda, &started); /* did we start above here? */
	PopNstack(pda, &ended);   /* did we end above here? */

	if (started && !ended) 
	  {
	    state = CalculateStateIndex(cm, node, BEGR_S);
	    tidx = InsertTraceNode(tr, tidx, TRACE_RIGHT_CHILD, 
				   gtr->emitl[node], gtr->emitr[node], state);
	  }
	state = CalculateStateIndex(cm, node, BEGR_IL);
	for (i = gtr->emitl[node]; i < gtr->emitl[gtr->nxtl[node]]; i++)
	  if (dsq[i] != DIGITAL_GAP) {
	    if (ended) goto FAILURE;
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, 
				   gtr->emitr[node], state);
	    if (! started) goto FAILURE;
	  }
	break;

	/* An END node.
	 * If we've already ended (on EL), skip. 
	 */
      case END_nd:
	if (started && ! ended) {
	  state = CalculateStateIndex(cm, node, END_E);
	  tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, -1, -1, state);
	}
	break;

      default: 
	Die("bogus node type %d in transmogrify()", gtr->state[node]);
      }
    }
  if (nstarts > 1) goto FAILURE;
  free(localrun);
  FreeNstack(pda);
  return tr;

 FAILURE:
  free(localrun);
  FreeNstack(pda);
  FreeParsetree(tr);
  Die("transmogrification failed: bad input sequence.");
  return NULL;			/* not reached */
}


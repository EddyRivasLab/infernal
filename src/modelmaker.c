/* modelmaker.c
 * SRE, 28 Feb 2000
 * CVS $Id$
 *
 * Construct a model from an alignment. 
 *
 * Outline of the process:
 *    1. construct a "master" traceback (mtr) for the alignment, 
 *       specifying which columns are match vs. insert and
 *       how the model tree branches. 
 *    2. This traceback is assigned a numbering system by NumberMasterTrace(),
 *       which returns the number of nodes; the caller then allocates a new CM. 
 *    3. This new model is numbered (assigned a branching structure) by TopofyNewCM(). 
 *    4. Individual tracebacks are constructed from individual aligned sequences 
 *       by Transmogrify().
 *    5. The individual tracebacks are counted into a new model with TraceCount().
 *    6. The counts converted to probabilities with ProbifyCM().
 *
 * The "master trace" is a special use of a Parsetree_t structure. 
 * - tr->state contains a node type (e.g. MATP_nd), not a state index.
 * - The numbering of the master trace is a postorder traverse, identical to
 *   the numbering in the final CM.
 * - emitl and emitr are relative to the alignment columns, not individual
 *   sequence positions.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */


#include <stdlib.h>
#include "structs.h"
#include "funcs.h"
#include "nstack.h"
#include "squid.h"		
#include "msa.h"		/* multiple sequence alignments */


static void         cm_from_master(CM_t *cm, Parsetree_t *mtr);
static Parsetree_t *transmogrify(CM_t *cm, Parsetree_t *mtr, char *aseq);

/* Function: HandModelmaker()
 * Incept:   SRE 29 Feb 2000 [Seattle]; from COVE 2.0 code
 * 
 * Purpose:  The customer always knows best.
 * 
 *           Construct a model given a stated structure. The structure
 *           is provided via a "ss_cons" (consensus structure) line, as would
 *           occur in an annotated SELEX or Stockoholm file. Only > and < characters
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
 *           use_rf    - TRUE to use RF annotation to determine match/insert
 *           gapthresh - fraction of gaps to allow in a match column (if use_rf is FALSE)
 *           ret_cm    - RETURN: new model                      (maybe NULL)
 *           ret_mtr   - RETURN: master traceback for alignment (maybe NULL)
 *           
 * Return:   void
 *           cm is allocated here. FreeCM(*ret_cm).
 *           tr is allocated here. FreeTrace() on each one, then free(*ret_tr).
 */
void
HandModelmaker(MSA *msa, int use_rf, float gapthresh, CM_t **ret_cm, Parsetree_t **ret_mtr)
{
  CM_t           *cm;		/* new covariance model                       */
  Parsetree_t    *mtr;		/* master traceback tree for alignment        */
  Parsetree_t    *tr;		/* an individual parse tree                   */
  Nstack_t       *pda;		/* pushdown stack used in building mtr        */
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

  if (msa->ss_cons == NULL)      Die("No consensus structure annotation available for that alignment.");
  if (use_rf && msa->rf == NULL) Die("No reference annotation available for that alignment.");

  /* 1. Determine match/insert assignments
   *    matassign is 0..alen-1. Values are 1 if a match column, 0 if insert column.
   */
  matassign = MallocOrDie(sizeof(int) * msa->alen);
  if (use_rf)
    {
      for (apos = 0; apos < msa->alen; apos++)
	matassign[apos] = (isgap(msa->rf[apos]) ? 0 : 1);
    }
  else
    {
      int gaps;
      for (apos = 0; apos < msa->alen; apos++)
	{
	  for (gaps = 0, idx = 0; idx < msa->nseq; idx++)
	    if (isgap(msa->aseq[idx][apos])) gaps++;
	  matassign[apos] = ((double) gaps / (double) msa->nseq > gapthresh) ? 0 : 1;
	}
    }

  /* 2. Determine a "ct" array, base-pairing partners for each position.
   *    Disallow/ignore pseudoknots. (That's what the FALSE flag does.)
   *    ct[] values give the index of a base pairing partner, or -1 for unpaired positions.
   */
  if (! KHS2ct(msa->ss_cons, msa->alen, FALSE, &ct))  
    Die("Consensus structure string is inconsistent"); 

  /* 3. Make sure the consensus structure "ct" is consistent with the match assignments.
   *    Wipe out all structure in insert columns; including the base-paired 
   *    partner of insert-assigned columns.
   */
  for (apos = 0; apos < msa->alen; apos++)
    if (! matassign[apos])
      { 
	if (ct[apos] != -1)  ct[ct[apos]] = -1;
	ct[apos] = -1;
      }

  /* 4. Construct a master traceback tree.
   *    This code is borrowed from yarn's KHS2Trace().
   *    
   *    We also keep track of how many states we'll need in the final CM,
   *    so we'll know how much to allocate -- and the number of nodes,
   *    for informational purposes.
   */
  nstates = nnodes = 0;
  mtr = CreateParsetree();	/* the parse tree we'll grow        */
  pda = CreateNstack();		/* a pushdown stack for our indices */

  /* Construction strategy has to make sure we number the nodes in
   * postorder traversal: for bifurcations, we can't attach the right 
   * child until we've fully traversed the left side. Therefore, we have
   * to push what we intend to attach, and pop it later. And since we
   * don't know an index for the node until we attach it, we have no
   * place to put the node's data except the stack -- so we have to
   * push several numbers onto the stack: what type of node, what
   * subseq it's responsible for (emitl...emitr), and what node
   * index it attaches to.
   */
  PushNstack(pda, -1);		/* what state it's attached to */
  PushNstack(pda, 0);		/* emitl */
  PushNstack(pda, msa->alen-1);	/* emitr */
  PushNstack(pda, ROOT_nd);	/* "state" (e.g. node type) */

  while (PopNstack(pda, &type))	/* pop a node type to attach */
    {
      PopNstack(pda, &j);
      PopNstack(pda, &i);	/* i..j == subseq we're responsible for */
      PopNstack(pda, &v);	/* v = index of parent node in mtr */

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
	v = InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, END_nd);
	nstates += 1;		/* END_nd -> E_st */
	nnodes++;
      }

      else if (type == ROOT_nd) { /* try to push i,j; but deal with INSL and INSR */
	v = InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, ROOT_nd);
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
	v = InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, BEGL_nd);
	PushNstack(pda, v);	
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 1;		/* BEGL_nd -> S_st */
	nnodes++;
      }

      else if (type == BEGR_nd)  { /* look for INSL */
	v = InsertTraceNode(mtr, v, TRACE_RIGHT_CHILD, i, j, BEGR_nd);
	for (; i <= j; i++) if (matassign[i]) break; 
	PushNstack(pda, v);	
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 2;		/* BEGR_nd -> S_st IL_st */
	nnodes++;
      }

      else if (ct[i] == -1) {
	 	/* i unpaired. This is a MATL node; allow INSL */
	v = InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
	for (i = i+1; i <= j; i++)  if (matassign[i]) break;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* MATL_nd -> ML_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[j] == -1) { 	/* j unpaired. MATR node. Deal with INSR */
	v = InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, MATR_nd);
	for (j = j-1; j >= i; j--) if (matassign[j]) break;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* MATR_nd -> MR_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[i] == j) { 	/* i,j paired to each other. MATP. deal with INSL, INSR */
	v = InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, MATP_nd);
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
	  /* Here's the first of two places where we can optimize the topology of a CM.
           * (The other comes from choosing a state traversal order when building the CM.)
           * Imagine a multifurcation of four domains:
           *   [1]..[2]..[3]..[4]
           * The "default leftwise" rule means that BEGL/INSL generates the intervening
           * sequences, so we must model the four stems as:
           *   [1],    ..[2],    ..[3],   ..[4]
           * but we have a choice of how we bifurcate:
           *   (1,2)(3,4)    1,(2,(3,4))   ((1,2),3),4 
           * Our choice may affect the time and memory requirements of a divide
           * conquer alignment algorithm.
           * 
           * So we may have to choose carefully from several possible split points k (3,
           * in the above example). A priori we only know one possible midpoint precisely:
           * ct[i]+1, the next base after closing domain 1. We can find the others
           * by scanning for them, and we can be reasonably efficient about scanning
           * by using ct[] to instantly skip subdomains.
	   */
	  /* One possible rule: optimize by finding most balanced split.
           * Each stop of the following loop gives a possible midpoint k, which is
           * then evaluated, keeping track of the best split so far.
           */
	  v = InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, BIF_nd);
	  bestdiff = msa->alen;
	  bestk    = ct[i]+1;
	  for (k = ct[i] + 1; k < ct[j]; k = ct[k] + 1) 
	    {
	      diff = abs(i+j-2*k); /* = len2-len1-1, where len2 = j-k+1, len1= k-i */
	      if (diff < bestdiff) {
		bestdiff = diff; 
		bestk    = k;
	      }
	      while (ct[k] == -1) k++;
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

  /* OK, we've converted ct into mtr -- mtr is a tree structure telling us the
   * arrangement of consensus nodes. Now do the drill for constructing a full model 
   * using this master trace.
   */
  cm = CreateCM(nnodes, nstates);
  cm_from_master(cm, mtr);

  for (idx = 0; idx < msa->nseq; idx++)
    {
      tr = transmogrify(cm, mtr, msa->aseq[idx]);
      /* printf("### parse tree for sequence %d\n", idx);
	 PrintParsetree(stdout, tr); 
      */
      ParsetreeCount(cm, tr, msa->aseq[idx], msa->wgt[idx]);
      FreeParsetree(tr);
    }
  CMSimpleProbify(cm);

  free(matassign);
  if (ret_cm  != NULL) *ret_cm  = cm;  else FreeCM(cm);
  if (ret_mtr != NULL) *ret_mtr = mtr; else FreeParsetree(mtr);
}



/* Function: cm_from_master()
 * Date:     SRE, Sat Jul 29 09:25:49 2000 [St. Louis]
 *
 * Purpose:  given a master traceback and an allocated CM, 
 *           fill in all the structural information of the CM.
 *           
 * Args:     cm  - allocated cm to construct
 *           mtr - master trace
 *
 * Returns:  (void)
 */
static void
cm_from_master(CM_t *cm, Parsetree_t *mtr)
{
  Nstack_t   *pda;              /* pushdown stack used for traversing mtr */
  int         v;		/* what node we're working on (in mtr index system)*/
  int         node;		/* what node this is (preorder traversal numbering of CM) */
  int         state;		/* what state this is (preorder traversal numbering of CM) */
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
      if      (mtr->state[v] == BIF_nd) {
	prvnodetype = mtr->state[mtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = BIF_nd;

	cm->sttype[state] = B_st;
	cm->ndidx[state]  = node;
	cm->stid[state]   = BIF_B;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = -1; /* we have to fill this in later, when we see the BEGR... */
	PushNstack(pda, state);	/* ...and this is the trick we use to remember the connection */
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;
	
	node++;
	PushNstack(pda, mtr->nxtr[v]);
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->state[v] == MATP_nd) {
	nxtnodetype = mtr->state[mtr->nxtl[v]];
	prvnodetype = mtr->state[mtr->prv[v]];

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
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->state[v] == MATL_nd) {
	nxtnodetype = mtr->state[mtr->nxtl[v]];
	prvnodetype = mtr->state[mtr->prv[v]];

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
	PushNstack(pda, mtr->nxtl[v]);
      }
      
      else if (mtr->state[v] == MATR_nd) {
	nxtnodetype = mtr->state[mtr->nxtl[v]];
	prvnodetype = mtr->state[mtr->prv[v]];

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
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->state[v] == BEGL_nd) {
	nxtnodetype = mtr->state[mtr->nxtl[v]];

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
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->state[v] == BEGR_nd) {
	int bifparent;

	nxtnodetype = mtr->state[mtr->nxtl[v]];

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
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->state[v] == ROOT_nd) {
	nxtnodetype = mtr->state[mtr->nxtl[v]];

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
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->state[v] == END_nd) {
	prvnodetype = mtr->state[mtr->prv[v]];

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



/* Function: transmogrify()
 * Date:     SRE, Mon Jul 31 14:30:58 2000 [St. Louis]
 *
 * Purpose:  Construct a "fake" parsetree for a given aligned sequence (aseq),
 *           given a new CM structure (cm) and a master structure tree (mtr).
 *
 * Args:     cm    - the new covariance model
 *           mtr   - master parsetree
 *           aseq  - an aligned sequence
 *
 * Returns:  the individual parse tree. 
 *           Caller is responsible for free'ing this.
 */
static Parsetree_t *
transmogrify(CM_t *cm, Parsetree_t *mtr, char *aseq)
{
  Parsetree_t *tr;
  int          node;		/* index of the node in *mtr* we're currently working on */
  int          state;		/* index of a state in the *CM*                          */
  int          type;		/* a unique statetype                                    */
  Nstack_t    *pda;             /* pushdown automaton for remembering positions in tr    */
  int          tidx;		/* index *in the parsetree tr* of the state we're supposed to attach to next */
  int          i,j;		/* coords in aseq */

  tr  = CreateParsetree();
  pda = CreateNstack();

  /* Because the mtr is already indexed in a postorder traversal,
   * we can postorder traverse it easily w/ a for loop...
   */
  tidx = -1;			/* first state to attach to; -1 is special case for attaching root */
  for (node = 0; node < cm->nodes; node++)
    {
      switch (mtr->state[node]) { /* e.g. switch on node type: */
      case BIF_nd:
	state = CalculateStateIndex(cm, node, BIF_B);
	tidx  = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, mtr->emitl[node], mtr->emitr[node], state);
	PushNstack(pda, tidx); /* remember this index in tr, we'll use it for BEGL and BEGR */
	PushNstack(pda, tidx);
	break;

      case MATP_nd:
	if (isgap(mtr->emitl[node])) {
	  if (isgap(aseq[mtr->emitr[node]])) type = MATP_D;
	  else                               type = MATP_MR;
	} else {
	  if (isgap(aseq[mtr->emitr[node]])) type = MATP_ML;
	  else                               type = MATP_MP;
	}
	state = CalculateStateIndex(cm, node, type);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, mtr->emitl[node], mtr->emitr[node], state);

	state = CalculateStateIndex(cm, node, MATP_IL);
	for (i = mtr->emitl[node]+1; i < mtr->emitl[mtr->nxtl[node]]; i++)
	  if (! isgap(aseq[i]))
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, mtr->emitr[node]-1, state);

	state = CalculateStateIndex(cm, node, MATP_IR);
	for (j = mtr->emitr[node]-1; j > mtr->emitr[mtr->nxtr[node]]; j--)
	  if (! isgap(aseq[j]))
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, j, state);	
	break;

      case MATL_nd:
	if (isgap(aseq[mtr->emitl[node]])) type = MATL_D;
	else                               type = MATL_ML;

	state = CalculateStateIndex(cm, node, type);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, mtr->emitl[node], mtr->emitr[node], state);

	state = CalculateStateIndex(cm, node, MATL_IL);
	for (i = mtr->emitl[node]+1; i < mtr->emitl[mtr->nxtl[node]]; i++)
	  if (! isgap(aseq[i]))
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, mtr->emitr[node], state);
	break;

      case MATR_nd:
	if (isgap(aseq[mtr->emitr[node]])) type = MATR_D;
	else                               type = MATR_MR;

	state = CalculateStateIndex(cm, node, type);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, mtr->emitl[node], mtr->emitr[node], state);

	state = CalculateStateIndex(cm, node, MATR_IR);
	for (j = mtr->emitr[node]-1; j > mtr->emitr[mtr->nxtl[node]]; j--)
	  if (! isgap(aseq[j]))
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, mtr->emitl[node], j, state);
	break;

      case BEGL_nd:
	PopNstack(pda, &tidx);	/* recover parent bifurcation's index in trace */
	
	state = CalculateStateIndex(cm, node, BEGL_S);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, mtr->emitl[node], mtr->emitr[node], state);
	break;

      case BEGR_nd:
	PopNstack(pda, &tidx);	/* recover parent bifurcation's index in trace */
	
	state = CalculateStateIndex(cm, node, BEGR_S);
	tidx = InsertTraceNode(tr, tidx, TRACE_RIGHT_CHILD, mtr->emitl[node], mtr->emitr[node], state);

	state = CalculateStateIndex(cm, node, BEGR_IL);
	for (i = mtr->emitl[node]; i < mtr->emitl[mtr->nxtl[node]]; i++)
	  if (! isgap(aseq[i]))
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, mtr->emitr[node], state);
	break;

      case ROOT_nd:
				/* we assume ROOT_S == 0, ROOT_IL == 1, ROOT_IR == 2 in the CM */
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, mtr->emitl[node], mtr->emitr[node], 0);
	for (i = mtr->emitl[node]; i < mtr->emitl[mtr->nxtl[node]]; i++)
	  if (! isgap(aseq[i]))
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, mtr->emitr[node], 1);
	for (j = mtr->emitr[node]; j > mtr->emitr[mtr->nxtr[node]]; j--)
	  if (! isgap(aseq[j]))
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, j, 2);	
	break;

      case END_nd:
	state = CalculateStateIndex(cm, node, END_E);
	tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, -1, -1, state);
	break;

      default: 
	Die("bogus node type %d in transmogrify()", mtr->state[node]);
      }
    }
  FreeNstack(pda);
  return tr;
}

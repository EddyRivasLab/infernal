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


static void cm_from_master(CM_t *cm, Parsetree_t *mtr);

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
  Nstack_t       *pda;		/* pushdown stack used in building mtr        */
  int            *matassign;	/* 0..alen-1 array; 0=insert col, 1=match col */
  int            *ct;		/* 0..alen-1 base pair partners array         */
  int             apos;		/* counter over columns of alignment          */
  int             idx;		/* counter over sequences in the alignment    */
  int             v;		/* index of current node                      */
  int             i,j,k;	/* subsequence indices                        */
  int  diff, bestdiff, bestk;   /* used while finding optimal split points    */   
  int  nxti, nxtj;		/* used when skipping over insertions         */
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
   *    mtr's emitl, emitr, and type are properly set by this section.
   *    We also keep track of how many states we'll need in the final CM,
   *    so we'll know how much to allocate -- and the number of nodes,
   *    for informational purposes.
   */
  nstates = nnodes = 0;
  mtr = CreateParsetree();	/* the parse tree we'll grow        */
  pda = CreateNstack();		/* a pushdown stack for our indices */
				/* attach the root, push onto stack */
  v = InsertTraceNode(mtr, 0, TRACE_LEFT_CHILD, 0, msa->alen-1, -1, ROOT_nd);
  PushNstack(pda, v);
  while (PopNstack(pda, &v))
    {
      i = mtr->emitl[v];
      j = mtr->emitr[v];
      
      /* This node accounts for i..j, but we don't know how yet.
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
	mtr->type[v] = END_nd;
	nstates += 1;		/* END_nd -> E_st */
	nnodes++;
      }

      else if (mtr->type[v] == ROOT_nd) { /* try to push i,j; but deal with INSL and INSR */
	for (; i <= j; i++) if (matassign[i]) break;
	for (; j >= i; j--) if (matassign[j]) break;
	PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));
	nstates += 3;		/* ROOT_nd -> S_st, IL_st, IR_st */
	nnodes++;
      }

      else if (mtr->type[v] == BEGL_nd) {    /* no inserts */
	PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));
	nstates += 1;		/* BEGL_nd -> S_st */
	nnodes++;
      }


      else if (mtr->type[v] == BEGR_nd)  { /* look for INSL */
	for (; i <= j; i++) if (matassign[i]) break; 
	PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));
	nstates += 2;		/* BEGR_nd -> S_st IL_st */
	nnodes++;
      }

      else if (ct[i] == -1) {
	 	/* i unpaired. This is a MATL node; allow INSL */
	mtr->type[v] = MATL_nd;
	for (i = i+1; i <= j; i++)  if (matassign[i]) break;
	PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));
	nstates += 3;		/* MATL_nd -> ML_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[j] == -1) { 	/* j unpaired. MATR node. Deal with INSR */
	mtr->type[v] = MATR_nd;
	for (nxtj = j-1; nxtj >= i; nxtj--) if (matassign[nxtj]) break;
	PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, nxtj, -1, DUMMY_nd));
	nstates += 3;		/* MATR_nd -> MR_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[i] == j) { 	/* i,j paired to each other. MATP. deal with INSL, INSR */
	mtr->type[v] = MATP_nd;
	for (nxti = i+1; nxti <= j; nxti++)    if (matassign[nxti]) break;
	for (nxtj = j-1; nxtj >= nxti; nxtj--) if (matassign[nxtj]) break;
	PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, nxti, nxtj, -1,DUMMY_nd));
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

	  mtr->type[v] = BIF_nd;
	  PushNstack(pda, InsertTraceNode(mtr, v, TRACE_RIGHT_CHILD, bestk, j,   -1,BEGR_nd));
	  PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD,  i, bestk-1, -1,BEGL_nd));
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
  cm = CreateCM(nstates);
  cm_from_master(cm, mtr);

#if 0
  for (idx = 0; idx < nseq; idx++)
    {
      Transmogrify(mtr, aseq[idx], &tr, &pool);
      if (! TraceCount(cm, aseq[idx], 
		       (ainfo->sqinfo[idx].flags & SQINFO_WGT) ? (double) ainfo->sqinfo[idx].weight : 1.0,
		       tr))
	Die("TraceCount() failed");
      FreeTrace(tr, pool);
    }
  ProbifyCM(cm, prior);
  
  free(matassign);
#endif /*0*/
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
  PushNstack(pda, 2);		/* push the ROOT_nd onto the stack. (0,1 are dummies) */
  while (PopNstack(pda, &v))
    {
      if      (mtr->type[v] == BIF_nd) {
	prvnodetype = mtr->type[mtr->prv[v]];

	cm->sttype[state] = B_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = BIF_nd;
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

      else if (mtr->type[v] == MATP_nd) {
	nxtnodetype = mtr->type[mtr->nxtl[v]];
	prvnodetype = mtr->type[mtr->prv[v]];

	cm->sttype[state] = MP_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATP_nd;
	cm->stid[state]   = MATP_MP;
	cm->cfirst[state] = state+4;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = ML_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATP_nd;
	cm->stid[state]   = MATP_ML;
	cm->cfirst[state] = state+3;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-2;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = MR_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATP_nd;
	cm->stid[state]   = MATP_MR;
	cm->cfirst[state] = state+2;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-3;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = D_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATP_nd;
	cm->stid[state]   = MATP_D;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state-4;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATP_nd;
	cm->stid[state]   = MATP_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 2 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 5;
	state++;

	cm->sttype[state] = IR_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATP_nd;
	cm->stid[state]   = MATP_IR;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 6;
	state++;

	node++;
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->type[v] == MATL_nd) {
	nxtnodetype = mtr->type[mtr->nxtl[v]];
	prvnodetype = mtr->type[mtr->prv[v]];

	cm->sttype[state] = ML_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATL_nd;
	cm->stid[state]   = MATL_ML;
	cm->cfirst[state] = state+2;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = D_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATL_nd;
	cm->stid[state]   = MATL_D;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-2;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATL_nd;
	cm->stid[state]   = MATL_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 3;
	state++;

	node++;
	PushNstack(pda, mtr->nxtl[v]);
      }
      
      else if (mtr->type[v] == MATR_nd) {
	nxtnodetype = mtr->type[mtr->nxtl[v]];
	prvnodetype = mtr->type[mtr->prv[v]];

	cm->sttype[state] = MR_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATR_nd;
	cm->stid[state]   = MATR_MR;
	cm->cfirst[state] = state+2;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = D_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATR_nd;
	cm->stid[state]   = MATR_D;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state-2;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;

	cm->sttype[state] = IR_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = MATR_nd;
	cm->stid[state]   = MATR_IR;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 3;
	state++;

	node++;
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->type[v] == BEGL_nd) {
	nxtnodetype = mtr->type[mtr->nxtl[v]];

	cm->sttype[state] = S_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = BEGL_nd;
	cm->stid[state]   = BEGL_S;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = child_count[nxtnodetype];
	cm->plast[state] = state-1;
	cm->pnum[state]   = 1;
	state++;

	node++;
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->type[v] == BEGR_nd) {
	int bifparent;

	nxtnodetype = mtr->type[mtr->nxtl[v]];

	/* A trick: we need to attach this start state to the previous
	 * bifurcation. We stored the bif state index by pushing it onto
	 * the pda -- retrieve it now.
	 */
	PopNstack(pda, &bifparent);
	cm->cnum[bifparent] = state; /* remember, cnum overloaded for bif: idx of right child */

	cm->sttype[state] = S_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = BEGR_nd;
	cm->stid[state]   = BEGR_S;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = bifparent;
	cm->pnum[state]   = 1;
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = BEGR_nd;
	cm->stid[state]   = BEGR_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype];
	cm->plast[state] = state;
	cm->pnum[state]   = 2;
	state++;

	node++;
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->type[v] == ROOT_nd) {
	nxtnodetype = mtr->type[mtr->nxtl[v]];

	cm->sttype[state] = S_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = ROOT_nd;
	cm->stid[state]   = ROOT_S;
	cm->cfirst[state] = state+1;
	cm->cnum[state]   = 2 + child_count[nxtnodetype]; 
	cm->plast[state] = -1;
	cm->pnum[state]   = 0;
	state++;

	cm->sttype[state] = IL_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = ROOT_nd;
	cm->stid[state]   = ROOT_IL;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 2 + child_count[nxtnodetype]; 
	cm->plast[state] = state;
	cm->pnum[state]   = 2;
	state++;

	cm->sttype[state] = IR_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = ROOT_nd;
	cm->stid[state]   = ROOT_IR;
	cm->cfirst[state] = state;
	cm->cnum[state]   = 1 + child_count[nxtnodetype]; 
	cm->plast[state] = state;
	cm->pnum[state]   = 3;
	state++;

	node++;
	PushNstack(pda, mtr->nxtl[v]);
      }

      else if (mtr->type[v] == END_nd) {
	prvnodetype = mtr->type[mtr->prv[v]];

	cm->sttype[state] = E_st;
	cm->ndidx[state]  = node;
	cm->ndtype[state] = END_nd;
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








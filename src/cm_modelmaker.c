/* cm_modelmaker.c
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
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "easel.h"		
#include "esl_msa.h"		
#include "esl_stack.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "funcs.h"
#include "structs.h"

static int check_for_pknots(char *cs, int alen);

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
 *           errbuf    - for error messages
 *           use_rf    - TRUE to use RF annotation to determine match/insert
 *           gapthresh - fraction of gaps to allow in a match column (if use_rf=FALSE)
 *           ret_cm    - RETURN: new model                      (maybe NULL)
 *           ret_gtr   - RETURN: guide tree for alignment (maybe NULL)
 *           
 * Return:   eslOK on success;
 *           eslEINCOMPAT on contract violation
 */
int
HandModelmaker(ESL_MSA *msa, char *errbuf, int use_rf, float gapthresh, CM_t **ret_cm, Parsetree_t **ret_gtr)
{
  int             status;
  CM_t           *cm;		/* new covariance model                       */
  Parsetree_t    *gtr;		/* guide tree for alignment                   */
  ESL_STACK      *pda;		/* pushdown stack used in building gtr        */
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
  int  clen;                    /* consensus length of the model              */

  if (msa->ss_cons == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "HandModelMaker(): No consensus structure annotation available for that alignment.");
  if (! (msa->flags & eslMSA_DIGITAL)) ESL_FAIL(eslEINCOMPAT, errbuf, "HandModelMaker(): MSA is not digitized.");
  if (use_rf && msa->rf == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "HandModelMaker(): No reference annotation available for the alignment.");

  /* 1. Determine match/insert assignments
   *    matassign is 1..alen. Values are 1 if a match column, 0 if insert column.
   */
  ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));

  /* Watch for off-by-one. rf is [0..alen-1]; matassign is [1..alen] */
  if (use_rf)
    {
      for (apos = 1; apos <= msa->alen; apos++)
	matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1]) ? FALSE : TRUE);
    }
  else
    {
      int gaps;
      for (apos = 1; apos <= msa->alen; apos++)
	{
	  for (gaps = 0, idx = 0; idx < msa->nseq; idx++)
	    if (esl_abc_XIsGap(msa->abc, msa->ax[idx][apos])) gaps++;
	  matassign[apos] = ((double) gaps / (double) msa->nseq > gapthresh) ? 0 : 1;
	}
    }

  /* 2. Determine a "ct" array, base-pairing partners for each position.
   *    Disallow/ignore pseudoknots by removing them prior to making the ct array.
   *    ct[] values give the index of a base pairing partner, or 0 for unpaired positions.
   *    Even though msa->ss_cons is in the 0..alen-1 coord system of msa, ct[]
   *    comes back in the 1..alen coord system of dsq.
   */
  esl_wuss_nopseudo(msa->ss_cons, msa->ss_cons); /* remove pknots in place */
  ESL_ALLOC(ct, (msa->alen+1) * sizeof(int));
  if (esl_wuss2ct(msa->ss_cons, msa->alen, ct) == eslESYNTAX)  
    cm_Fail("Consensus structure string is inconsistent"); 
  else if (esl_wuss2ct(msa->ss_cons, msa->alen, ct) != eslOK)  goto ERROR;

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
  gtr = CreateParsetree(100);	/* the parse tree we'll grow        */
  pda = esl_stack_ICreate();    /* a pushdown stack for our indices */
  if(pda == NULL) goto ERROR;
  clen = 0;

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
  if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* what node it's attached to */
  if((status = esl_stack_IPush(pda, 1))  != eslOK) goto ERROR;		/* emitl */
  if((status = esl_stack_IPush(pda, msa->alen)) != eslOK) goto ERROR;	/* emitr */
  if((status = esl_stack_IPush(pda, ROOT_nd)) != eslOK)   goto ERROR;	/* "state" (e.g. node type) */

  while (esl_stack_IPop(pda, &type) != eslEOD) /* pop a node type to attach */
    {
      esl_stack_IPop(pda, &j);
      esl_stack_IPop(pda, &i); /* i..j == subseq we're responsible for */
      esl_stack_IPop(pda, &v); /* v = index of parent node in gtr */

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
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	/* here v==0 always. */
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 3;		/* ROOT_nd -> S_st, IL_st, IR_st */
	nnodes++;
      }

      else if (type == BEGL_nd) {    /* no inserts */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BEGL_nd);
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 1;		/* BEGL_nd -> S_st */
	nnodes++;
      }

      else if (type == BEGR_nd)  { /* look for INSL */
	v = InsertTraceNode(gtr, v, TRACE_RIGHT_CHILD, i, j, BEGR_nd);
	for (; i <= j; i++) if (matassign[i]) break; 
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 2;		/* BEGR_nd -> S_st IL_st */
	nnodes++;
      }

      else if (ct[i] == 0) {
	 	/* i unpaired. This is a MATL node; allow INSL */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
	for (i = i+1; i <= j; i++)  if (matassign[i]) break;
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 3;		/* MATL_nd -> ML_st, D_st, IL_st */
	nnodes++;
	clen += 1;
      }

      else if (ct[j] == 0) { 	/* j unpaired. MATR node. Deal with INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATR_nd);
	for (j = j-1; j >= i; j--) if (matassign[j]) break;
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 3;		/* MATR_nd -> MR_st, D_st, IL_st */
	nnodes++;
	clen += 1;
      }

      else if (ct[i] == j) { /* i,j paired to each other. MATP. deal with INSL, INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATP_nd);
	for (i = i+1; i <= j; i++) if (matassign[i]) break;
	for (j = j-1; j >= i; j--) if (matassign[j]) break;
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 6;		/* MATP_nd -> MP_st, ML_st, MR_st, D_st, IL_st, IR_st */
	nnodes++;
	clen += 2;
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
	  if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	  if((status = esl_stack_IPush(pda, bestk)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, BEGR_nd)) != eslOK) goto ERROR;
				/* then push the left BEGIN node */
	  if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	  if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, bestk-1)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, BEGL_nd)) != eslOK) goto ERROR;
	  nstates += 1;		/* BIF_nd -> B_st */
	  nnodes++;
	}

    }	/* while something's on the stack */
  esl_stack_Destroy(pda);
  free(ct);

  /* OK, we've converted ct into gtr -- gtr is a tree structure telling us the
   * arrangement of consensus nodes. Now do the drill for constructing a full model 
   * using this guide tree.
   */
  cm = CreateCM(nnodes, nstates, msa->abc);
  if((status = cm_from_guide(cm, errbuf, gtr, FALSE)) != eslOK) return status; /* FALSE says, we're not building a sub CM that will never be localized */
  CMZero(cm);
  cm->clen = clen;

  free(matassign);
  if (ret_cm  != NULL) *ret_cm  = cm;  else FreeCM(cm);
  if (ret_gtr != NULL) *ret_gtr = gtr; else FreeParsetree(gtr);
  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "HandModelMaker(): memory allocation error.");
  return eslEMEM; /* never reached */
}


/* Function: cm_from_guide()
 * Date:     SRE, Sat Jul 29 09:25:49 2000 [St. Louis]
 *
 * Purpose:  given a guide tree and an allocated CM, 
 *           fill in all the structural information of the CM.
 *           
 * Args:     cm  - allocated cm to construct
 *           errbuf - for error messages
 *           gtr - guide tree
 *           will_never_localize- TRUE if we're building a sub CM that we will never localize.
 *                                This is only relevant b/c we can allow 'invalid' CMs in this case.
 *                                An invalid CM is one that, if localized, could not generate all 
 *                                possible sequences (see comments in code below). 
 *                                We allow sub CMs to be invalid b/c we don't want cmalign to die
 *                                when a target seq results in a sub CM that is invalid in the middle
 *                                of a run. This is a pure hack and relies UNSAFELY on the assumption
 *                                that the sub CM will never be localized (though in the current
 *                                implementation sub CMs are only used by cmalign in global mode, thus
 *                                they never are localized). Still, we don't raise a flag in the 
 *                                CM to prevent downstream localization, which is dangerous - if the
 *                                implementation changes to allow sub CMs to become localized. Even then
 *                                though the risk is small b/c an invalid CM only is a problem if we
 *                                try to align a single residue sequence to it (again, see comments below
 *                                for more explanation).
 *
 * Returns:  eslOK on success;
 */
int
cm_from_guide(CM_t *cm, char *errbuf, Parsetree_t *gtr, int will_never_localize)
{
  int         status;
  ESL_STACK  *pda;              /* pushdown stack used for traversing gtr */
  int         v;		/* what node we're working on (in gtr index system)*/
  int         node;		/* what node (preorder traversal numbering of CM) */
  int         state;		/* what state (preorder traversal numbering of CM) */
  int         clen;		/* current count of consensus length   */
  int  nxtnodetype;		/* type of a child node (e.g. MATP_nd) */
  int  prvnodetype;		/* type of a parent node (e.g. MATP_nd) */

  /* Some CM structural configuration info:
   * child_count[] gives how many states are connectable in a child node. 
   * parent_count[] gives how many states are connectable in a parent node.
   */
 				/* BIF, MATP, MATL, MATR, BEGL, BEGR, ROOT, END */  
  int child_count[] =             {  1,    4,    2,    2,    1,    1,    0,   1};
  int parent_count[] =            {  1,    6,    3,    3,    1,    2,    3,   0};

  node = state = clen = 0;
  pda = esl_stack_ICreate();
  if((status = esl_stack_IPush(pda, 0)) != eslOK) goto ERROR;		/* push ROOT_nd onto the stack */
  while (esl_stack_IPop(pda, &v) != eslEOD)
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
	if((status = esl_stack_IPush(pda, state)) != eslOK) goto ERROR;	/* ... the trick we use to remember the connection */
	cm->plast[state] = state-1;
	cm->pnum[state]   = parent_count[prvnodetype];
	state++;
	
	node++;
	if((status = esl_stack_IPush(pda, gtr->nxtr[v])) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, gtr->nxtl[v])) != eslOK) goto ERROR;
      }

      else if (gtr->state[v] == MATP_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = MATP_nd;
	clen             += 2;

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
	if((status = esl_stack_IPush(pda, gtr->nxtl[v])) != eslOK) goto ERROR;
      }

      else if (gtr->state[v] == MATL_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = MATL_nd;
	clen             += 1;

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
	if((status = esl_stack_IPush(pda, gtr->nxtl[v])) != eslOK) goto ERROR;
      }
      
      else if (gtr->state[v] == MATR_nd) {
	nxtnodetype = gtr->state[gtr->nxtl[v]];
	prvnodetype = gtr->state[gtr->prv[v]];

	cm->nodemap[node] = state;
	cm->ndtype[node ] = MATR_nd;
	clen             += 1;

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
	if((status = esl_stack_IPush(pda, gtr->nxtl[v])) != eslOK) goto ERROR;
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
	if((status = esl_stack_IPush(pda, gtr->nxtl[v])) != eslOK) goto ERROR;
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
	esl_stack_IPop(pda, &bifparent);
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
	if((status = esl_stack_IPush(pda, gtr->nxtl[v])) != eslOK) goto ERROR;
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
	if((status = esl_stack_IPush(pda, gtr->nxtl[v])) != eslOK) goto ERROR;
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
  esl_stack_Destroy(pda);
  cm->M     = state;
  cm->nodes = node;
  cm->clen  = clen;

  /* A couple of checks to make sure our CM is valid for local alignment/search.
   * The following is invalid:
   * 1. CMs with exactly 3 nodes. This must be either {ROOT, MATL, END} or
   *    {ROOT, MATP, END}. Either way a local end is impossible b/c local ends
   *    from nodes adjacent to end states are impossible. This is bad. Even
   *    worse is a {ROOT, MATL, END} model can't emit/align more than a single
   *    residue in local mode (ROOT_IL, ROOT_IR are unreachable, and so is MATL_IL,
   *    b/c it was detached to remove an ambiguity with ROOT_IR).
   * 2. CMs with 0 MATL, MATR and BIF nodes. The reason is because such a CM only has
   *    a ROOT, a bunch of MATPs and an END, and it is impossible to align a single
   *    residue to such a model when it's in local mode. 
   */
  if(cm->nodes == 3) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_from_guide(), it's illegal to construct a CM of only 3 nodes."); 
  if((CMCountNodetype(cm, MATL_nd) == 0) && (CMCountNodetype(cm, MATR_nd) == 0) && (CMCountNodetype(cm,BIF_nd) == 0)) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_from_guide(), it's illegal to construct a CM with 0 MATL, MATR and BIF nodes."); 
  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "cm_from_guide(): memory allocation error.");
  return eslEMEM; /* NEVERREACHED */
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
 *           it's polite to cm_Fail() on any kind of input error.
 * 
 * Args:     cm    - the newly built covariance model, corresponding to gtr
 *           gtr   - guide tree
 *           ax    - a digitized aligned sequence [1..alen]
 *           aseq  - aligned sequence itself [0..alen-1]
 *           alen  - length of alignment
 *
 * Returns:  the individual parse tree. 
 *           Caller is responsible for free'ing this w/ FreeParsetree().
 */
Parsetree_t *
Transmogrify(CM_t *cm, Parsetree_t *gtr, ESL_DSQ *ax, char *aseq, int alen)
{
  int          status;
  Parsetree_t *tr;
  int          node;		/* index of node in *gtr* we're working on */
  int          state;		/* index of a state in the *CM*            */
  int          type;		/* a unique statetype                      */
  ESL_STACK   *pda;             /* pushdown automaton for positions in tr  */
  int          tidx;		/* index *in parsetree tr* of state        */
  int          i,j;		/* coords in aseq                          */
  int          started;		/* TRUE if we've transited out of ROOT     */
  int          ended;		/* TRUE if we've transited to EL and ended */
  int          nstarts;         /* # of local transits out of ROOT: <= 1   */
  int         *localrun;        /* local alignment gap run lengths         */
  int          need_leftside;
  int          need_rightside;

  tr  = CreateParsetree(100);
  pda = esl_stack_ICreate();
  
  started = FALSE;
  ended   = FALSE;
  nstarts = 0;

  /* We preprocess the aseq to help with local alignment.
   */
  ESL_ALLOC(localrun, sizeof(int) * (alen+1));
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
	  if (!esl_abc_XIsGap(cm->abc, ax[i])) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   i, gtr->emitr[node], 1);
	    if (! started) { started = TRUE; nstarts++; }
	  }
	for (j = gtr->emitr[node]; j > gtr->emitr[gtr->nxtl[node]]; j--)
	  if (!esl_abc_XIsGap(cm->abc, ax[j])) {
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
	if((status = esl_stack_IPush(pda, ended)) != eslOK) goto ERROR;   /* remember our ending status */
	if((status = esl_stack_IPush(pda, started)) != eslOK) goto ERROR; /* remember our start status */
	if((status = esl_stack_IPush(pda, tidx)) != eslOK) goto ERROR;    /* remember index in tr; we pop in BEGR */
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
	if (esl_abc_XIsGap(cm->abc, ax[gtr->emitl[node]])) {
	  if (esl_abc_XIsGap(cm->abc, ax[gtr->emitr[node]])) type = MATP_D;
	  else                                               type = MATP_MR;
	} else {
	  if (esl_abc_XIsGap(cm->abc, ax[gtr->emitr[node]])) type = MATP_ML;
	  else                                               type = MATP_MP;
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
	  if (!esl_abc_XIsGap(cm->abc, ax[i])) {
	    tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				   i, gtr->emitr[node]-1, state);
	    if (! started) goto FAILURE;
	  }

	state = CalculateStateIndex(cm, node, MATP_IR);
	if (cm->ndtype[gtr->nxtl[node]] != END_nd) { /* MATP_IR's before an END are detached, and transitions are impossible, this is a hack to deal with ambiguity in CM parsetree */
	  for (j = gtr->emitr[node]-1; j > gtr->emitr[gtr->nxtl[node]]; j--)
	    if (!esl_abc_XIsGap(cm->abc, ax[j])) {
	      tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, i, j, state);	
	      if (! started) goto FAILURE;
	    }
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
	if (esl_abc_XIsGap(cm->abc, ax[gtr->emitl[node]])) type = MATL_D;
	else                                               type = MATL_ML;

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
	if (cm->ndtype[gtr->nxtl[node]] != END_nd) { /* MATL_IL's before an END are detached, and transitions are impossible, this is a hack to deal with ambiguity in CM parsetree */
	  for (i = gtr->emitl[node]+1; i < gtr->emitl[gtr->nxtl[node]]; i++)
	    if (!esl_abc_XIsGap(cm->abc, ax[i])) {
	      tidx = InsertTraceNode(tr, tidx, TRACE_LEFT_CHILD, 
				     i, gtr->emitr[node], state);
	      if (! started) goto FAILURE;
	    }
	}
	break;

	/* MATR node. 
	 * Similar logic as MATL above.
	 */
      case MATR_nd:
	if (esl_abc_XIsGap(cm->abc, ax[gtr->emitr[node]])) type = MATR_D;
	else                                               type = MATR_MR;

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
	  if (!esl_abc_XIsGap(cm->abc, ax[j])) {
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
	esl_stack_IPop(pda, &tidx);	  /* recover parent bifurcation's index in trace */
	esl_stack_IPop(pda, &started); /* did we start above here? */
	esl_stack_IPop(pda, &ended);   /* did we end above here? */

	if (started && !ended) 
	  {
	    state = CalculateStateIndex(cm, node, BEGR_S);
	    tidx = InsertTraceNode(tr, tidx, TRACE_RIGHT_CHILD, 
				   gtr->emitl[node], gtr->emitr[node], state);
	  }
	state = CalculateStateIndex(cm, node, BEGR_IL);
	for (i = gtr->emitl[node]; i < gtr->emitl[gtr->nxtl[node]]; i++)
	  if (!esl_abc_XIsGap(cm->abc, ax[i])) {
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
	cm_Fail("bogus node type %d in transmogrify()", gtr->state[node]);
      }
    }
  if (nstarts > 1) goto FAILURE;
  free(localrun);
  esl_stack_Destroy(pda);
  return tr;

 FAILURE:
  free(localrun);
  esl_stack_Destroy(pda);
  FreeParsetree(tr);
  cm_Fail("transmogrification failed: bad input sequence.");
  return NULL;			/* not reached */

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL;			/* not reached */

}

/* Function: ConsensusModelmaker()
 * EPN 08.29.06 based closely on HandModelMaker:
 *              SRE 29 Feb 2000 [Seattle]; from COVE 2.0 code
 * 
 * Purpose:  Construct a model given a stated structure. The structure
 *           is provided via a "ss_cons" (consensus structure) line, as would
 *           occur in an annotated SELEX or Stockholm file. Only > and < characters
 *           in this line are interpreted (as base pairs). Pseudoknots, 
 *           if annotated, are ignored.
 *           
 *           All positions/columns of the given structure are considered 
 *           consensus, and this is the difference b/t this function and
 *           HandModelmaker. Also, this function does not take in a MSA 
 *           data structure. It was originally written for building a 
 *           new CM (a sub CM) that models a contiguous subset of columns
 *           of it's template (mother) CM.
 *           
 * Args:     abc       - the alphabet
 *           errbuf    - for error messages
 *           ss_cons   - input consensus structure string 
 *           clen      - length of ss_cons, number of consensus columns
 *           building_sub_model - TRUE if building a sub CM (usually TRUE)
 *           ret_cm    - RETURN: new model                      (maybe NULL)
 *           ret_gtr   - RETURN: guide tree for alignment (maybe NULL)
 *           
 * Return:   eslOK on success;
 *           eslEINCOMPAT on contract violation
 */
int
ConsensusModelmaker(const ESL_ALPHABET *abc, char *errbuf, char *ss_cons, int clen, int building_sub_model, CM_t **ret_cm, Parsetree_t **ret_gtr)
{
  int             status;
  CM_t           *cm;		/* new covariance model                       */
  Parsetree_t    *gtr;		/* guide tree for alignment                   */
  ESL_STACK      *pda;		/* pushdown stack used in building gtr        */
  int            *ct;		/* 0..alen-1 base pair partners array         */
  int             v;		/* index of current node                      */
  int             i,j,k;	/* subsequence indices                        */
  int  type;			/* type of node we're working on              */
  int  diff, bestdiff, bestk;   /* used while finding optimal split points    */   
  int  nnodes;			/* number of nodes in CM                      */
  int  nstates;			/* number of states in CM                     */
  int  obs_clen;                /* observed (MATL+MATR+2*MATP) consensus len  */

  if (ss_cons == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "No consensus structure annotation available in ConsensusModelmaker().");

  /* 1. Determine a "ct" array, base-pairing partners for each position.
   *    Disallow/ignore pseudoknots by removing them prior to making the ct array.
   *    ct[] values give the index of a base pairing partner, or 0 for unpaired positions.
   *    Even though ss_cons is in the 0..clen-1 coord system of msa, ct[]
   *    comes back in the 1..alen coord system of the sequence.
   */
  esl_wuss_nopseudo(ss_cons, ss_cons); /* remove pknots in place */
  ESL_ALLOC(ct, (clen+1) * sizeof(int));
  if ((status = esl_wuss2ct(ss_cons, clen, ct)) != eslOK) ESL_FAIL(status, errbuf, "Consensus string is inconsisent in ConsensusModelMaker().", status);

  /* 2. Construct a guide tree. 
   *    This codes is borrowed from HandModelmaker(), where it
   *    was originally borrowed from yarn's KHS2Trace().
   *    
   *    We also keep track of how many states we'll need in the final CM,
   *    so we'll know how much to allocate -- and the number of nodes,
   *    for informational purposes.
   */
  nstates = nnodes = 0;
  gtr = CreateParsetree(100);	/* the parse tree we'll grow        */
  pda = esl_stack_ICreate();    /* a pushdown stack for our indices */
  obs_clen = 0;

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
  if((status = esl_stack_IPush(pda, -1)) != eslOK) goto ERROR;		/* what node it's attached to */
  if((status = esl_stack_IPush(pda, 1)) != eslOK) goto ERROR;		/* emitl */
  if((status = esl_stack_IPush(pda, clen)) != eslOK) goto ERROR;	/* emitr */
  if((status = esl_stack_IPush(pda, ROOT_nd)) != eslOK) goto ERROR;	/* "state" (e.g. node type) */

  while (esl_stack_IPop(pda, &type) != eslEOD) /* pop a node type to attach */
    {
      esl_stack_IPop(pda, &j);
      esl_stack_IPop(pda, &i); /* i..j == subseq we're responsible for */
      esl_stack_IPop(pda, &v); /* v = index of parent node in gtr */

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
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	/* here v==0 always. */
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 3;		/* ROOT_nd -> S_st, IL_st, IR_st */
	nnodes++;
      }

      else if (type == BEGL_nd) {    /* no inserts */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BEGL_nd);
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 1;		/* BEGL_nd -> S_st */
	nnodes++;
      }

      else if (type == BEGR_nd)  { /* look for INSL */
	v = InsertTraceNode(gtr, v, TRACE_RIGHT_CHILD, i, j, BEGR_nd);
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 2;		/* BEGR_nd -> S_st IL_st */
	nnodes++;
      }

      else if (ct[i] == 0) {
	 	/* i unpaired. This is a MATL node; allow INSL */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
	i++;
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 3;		/* MATL_nd -> ML_st, D_st, IL_st */
	nnodes++;
	obs_clen++;
      }

      else if (ct[j] == 0) { 	/* j unpaired. MATR node. Deal with INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATR_nd);
	j--;
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 3;		/* MATR_nd -> MR_st, D_st, IL_st */
	nnodes++;
	obs_clen++;
      }

      else if (ct[i] == j) { /* i,j paired to each other. MATP. deal with INSL, INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATP_nd);
	i++;
	j--;
	if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, DUMMY_nd)) != eslOK) goto ERROR; /* we don't know yet what the next node will be */
	nstates += 6;		/* MATP_nd -> MP_st, ML_st, MR_st, D_st, IL_st, IR_st */
	nnodes++;
	obs_clen += 2;
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
	  bestdiff = clen;
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
	  if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	  if((status = esl_stack_IPush(pda, bestk)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, BEGR_nd)) != eslOK) goto ERROR;
				/* then push the left BEGIN node */
	  if((status = esl_stack_IPush(pda, v)) != eslOK) goto ERROR;	
	  if((status = esl_stack_IPush(pda, i)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, bestk-1)) != eslOK) goto ERROR;
	  if((status = esl_stack_IPush(pda, BEGL_nd)) != eslOK) goto ERROR;
	  nstates += 1;		/* BIF_nd -> B_st */
	  nnodes++;
	}
    }	/* while something's on the stack */
  if(obs_clen != clen) cm_Fail("ConsensusModelMaker(): obs_clen: %d != passed in clen: %d\n", obs_clen, clen);
  esl_stack_Destroy(pda);
  free(ct);

  /* OK, we've converted ct into gtr -- gtr is a tree structure telling us the
   * arrangement of consensus nodes. Now do the drill for constructing a full model 
   * using this guide tree.
   */
  cm = CreateCM(nnodes, nstates, abc);
  if((status = cm_from_guide(cm, errbuf, gtr, building_sub_model)) != eslOK) return status;
  CMZero(cm);
  cm->clen = clen;

  if (ret_cm  != NULL) *ret_cm  = cm;  else FreeCM(cm);
  if (ret_gtr != NULL) *ret_gtr = gtr; else FreeParsetree(gtr);
  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "ConsensusModelMaker(): memory allocation error.");
  return eslEMEM; /* never reached */

}

/**************************************************************************
 * EPN 09.25.06 [Parkway Hotel St. Louis, MO]
 * cm_find_and_detach_dual_inserts()
 *
 * Given a CM (potentially in counts form), find cases where two 
 * insert states insert at the same position (due to an ambiguity in the 
 * CM architecture). We know from the way CMs are constructed in 
 * HandModelMaker() that one of these states must be an IL or IR state 
 * immediately prior to an END_E state, and by the way counts are 
 * collected in ParseTreeCount() that this END_E-1 state will 
 * not be filled with any counts from the input seed sequences. 
 * However, to be safe, there's an option to this function
 * to check to make sure both of these guarantees hold. 
 * 
 * Usually if this option is enabled with do_check=TRUE, the 
 * CM is in counts form so its possible to check to make sure
 * the END_E-1 has 0 counts. Also, in this case the other option,
 * do_detach is set to FALSE to tell the function not to detach
 * the insert quite yet, we want to wait until the model has
 * been priorified - and once it has we revisit this function
 * with the do_check option as FALSE and do_detach as TRUE.
 * With do_detach == TRUE, we 'detach' the END_E-1 state by
 * setting all transitions into it as 0.0, making it impossible
 * to reach.
 *
 * There should be exactly 1 dual insert for every END_E state.
 *
 * Args:    
 * CM_t  cm,
 * int   do_check;           TRUE to check the 
 * Returns: TRUE on success if all dual inserts are found and detached.
 */
int
cm_find_and_detach_dual_inserts(CM_t *cm, int do_check, int do_detach)
{

  int          status;
  CMEmitMap_t *emap;         /* consensus emit map for the cm */
  int *cc2lins_map;
  int *cc2rins_map;
  int cc;
  int nd;
  int end_e_ct;
  int detach_ct;
  int v;

  end_e_ct = 0;
  detach_ct = 0;

  /* Determine the number of END_E states in the model, this 
   * will be the number of inserts we want to detach.
   */
  for(v = 0; v <= cm->M; v++)
    if(cm->sttype[v] == E_st)
      end_e_ct++;

  emap = CreateEmitMap(cm);
  /*DumpEmitMap(stdout, emap, cm);*/

  /* Based on the emitmap, make map of which nodes have an insert state that
   * inserts AFTER (in case of *lmap) and BEFORE (in case of *rmap)
   * each consensus node. 
   * cc2lins_map[cc] = v, where cm state v is an IL_st that
   *                           emits after consensus column cc.
   * cc2rins_map[cc] = v, where cm state v is an IR_st that
   *                           emits after consensus column cc.
   * if no such state exists the value will be -1.
   */

  /* Allocate and initialize */
  ESL_ALLOC(cc2lins_map, sizeof(int) * (emap->clen + 1));
  ESL_ALLOC(cc2rins_map, sizeof(int) * (emap->clen + 1));
  for(cc = 0; cc <= emap->clen; cc++)
    {
      cc2lins_map[cc] = -1;
      cc2rins_map[cc] = -1;
    }
  /* fill in the map */
  /* ROOT is special */
  cc2lins_map[0] = 1; /* ROOT_IL */
  cc2rins_map[emap->clen] = 2; /* ROOT_IR */
  for(nd = 0; nd < cm->nodes; nd++)
    {
      switch (cm->ndtype[nd]) {
      case MATP_nd:
	cc2lins_map[emap->lpos[nd]] = cm->nodemap[nd] + 4; /* MATP_IL */
	cc2rins_map[emap->rpos[nd] - 1] = cm->nodemap[nd] + 5; /* MATP_IR */
	break;
	
      case MATL_nd:
	cc2lins_map[emap->lpos[nd]] = cm->nodemap[nd] + 2; /* MATL_IL */
	break;
	
      case MATR_nd:
	cc2rins_map[emap->rpos[nd] - 1] = cm->nodemap[nd] + 2; /* MATR_IR */
	break;
	
      case BEGR_nd:
	cc2lins_map[emap->lpos[nd]] = cm->nodemap[nd] + 1; /* BEGR_IL */
	break;
	
      default: {} /*do nothing*/
      }
    }

  for(cc = 0; cc <= emap->clen; cc++)
    {
      if(cc2lins_map[cc] != -1 && cc2rins_map[cc] != -1)
	{
	  detach_ct++;
	  /* Found a dual insert. */
	  if(do_check)
	    {
	      if(!(cm_check_before_detaching(cm, cc2lins_map[cc], cc2rins_map[cc])))
		cm_Fail("ERROR cm_check_before_detaching() returned false\n");		 
	    }
	  if(do_detach)
	    if(!(cm_detach_state(cm, cc2lins_map[cc], cc2rins_map[cc])))
	      cm_Fail("ERROR cm_detach_state() returned false\n");		 
	}
    }

  FreeEmitMap(emap);
  free(cc2lins_map);
  free(cc2rins_map);

  if(detach_ct != end_e_ct)
    return FALSE;
  else
    return TRUE;

 ERROR:
  cm_Fail("Memory allocation error.");
  return FALSE; /* never reached */
}

/**************************************************************************
 * EPN 09.18.06
 * cm_detach_state()
 *
 * Given two insert states that map to the same state in the original,
 * template CM, detach one of the them from the rest of the model by 
 * setting all transitions into it to 0.0. We choose the state to detach
 * as the state that is immediately prior to an END_E in the sub_cm. 
 * This will always be the case because the sole source of alignment 
 * ambiguity in the CM architecture always involves 1 insert state 
 * immediately before an END_E.
 * 
 * Args:    
 * CM_t  cm,
 * int   insert1; 
 * int   insert2; 
 * Returns: TRUE on success if insert1 or insert2 is a state immediately 
 *          before an END_E, and we've detached it.
 *          FALSE otherwise
 */
int
cm_detach_state(CM_t *cm, int insert1, int insert2)
{
  /*printf("\t**in cm_detach_state: insert1: %d | insert2: %d\n", insert1, insert2);*/

  int ret_val;
  int x, y;
  int to_detach;
  int x_offset;

  ret_val = FALSE;

  if(insert1 == insert2)
    cm_Fail("ERROR in cm_detach_state: insert1==insert2:%d\n", insert1);

  if(cm->sttype[insert1+1] == E_st)
    {
      ret_val = TRUE;
      to_detach = insert1;
    }
  else
    {
      if(cm->sttype[insert2+1] != E_st)
	cm_Fail("ERROR: in cm_detach_state insert1: %d and insert2: %d neither map to END_E-1 states.\n", insert1, insert2);
      if(ret_val)
	cm_Fail("ERROR: in cm_detach_state insert1: %d and insert2: %d both map to END_E-1 states.\n", insert1, insert2);
      ret_val = TRUE;
      to_detach = insert2;
    }
  if(ret_val)
    {
      /* Determine if we're detaching an IL_st, or the rare case of a MATP_IR st */
      if(cm->sttype[to_detach] == IL_st)
	x_offset = 0;
      else
	{
	  if(cm->stid[to_detach] != MATP_IR)
	    cm_Fail("ERROR: in cm_detach_state trying to detach a non-IL, non-MATP_IR state!\n");
	  x_offset = 1; /* MATP_* -> MATP_IR is second possible transition for MATP_*,
			 * unless * == MATP_IR, but we don't get there in for loop below. */
	}
      for (y = cm->pnum[to_detach]-1; y >= 1; y--)  
	/* y >= 1 means we never get to 
	 * to_detach->to_detach prob, which is irrelevant. */
	{
	  x = cm->plast[to_detach] - y;
	  cm->t[x][x_offset] = 0.0; /* x is a split set state in same node
			      * as insert1, we're setting transition
			      * from x -> to_detach as impossible.
			      */
	  /* Renormalize transitions out of x */
	  esl_vec_FNorm(cm->t[x], cm->cnum[x]);
	  /*printf("****setting transition probabilitity of x: %d to to_detach: %d cm->t[x][%d] as 0.0\n", x, to_detach, x_offset);*/
	}
    }
  return ret_val;
}

/**************************************************************************
 * EPN 09.25.06
 * cm_check_before_detaching()
 *
 * Given two insert states that map to the same state in a CM, and
 * given the CM in counts form (after being filled with counts from 
 * the parses implicit in the seed alignment), check the following
 * two guarantees are met:
 * (a) exactly one of the two inserts is immediately prior to an END_E
 *     state (END_E - 1).
 * (b) the END_E - 1 state has been parameterized with 0 counts from
 *     the input alignment.
 *
 * NOTE: A special case is when insert1 is the MATP_IL and insert2 the 
 *       MATP_IR of the same MATP node. This is the only case where
 *       the insert state to be detached is not an IL state, but rather
 *       the MATP_IR (but the guarantees still hold, the MATP_IR is 
 *       END_E - 1, and always gets 0 counts)
 *
 * Args:    
 * CM_t  cm,
 * int   insert1; 
 * int   insert2; 
 * Returns: TRUE on success if insert1 or insert2 follow the guarantees
 *          FALSE otherwise
 */
int
cm_check_before_detaching(CM_t *cm, int insert1, int insert2)
{
  int ret_val;
  int i, yoffset;
  int to_detach;
  int to_keep;
  float diff;
  
  ret_val = FALSE;

  if(insert1 == insert2)
    cm_Fail("ERROR in cm_check_before_detaching(), insert1==insert2 (%d)\n", insert1);

  if(cm->sttype[insert1+1] == E_st)
    {
      ret_val = TRUE;
      to_detach = insert1;
      to_keep   = insert2;
    }
  if(cm->sttype[insert2+1] == E_st)
    {
      if(ret_val)
	cm_Fail("ERROR: in cm_check_before_detaching() insert1: %d and insert2: %d both map to END_E-1 states.\n", insert1, insert2);
      ret_val = TRUE;
      to_detach = insert2;
      to_keep   = insert1;
    }

  /* check to make sure we have 0.0 counts in to_detach */
  if(ret_val)
    {
      for(i = 0; i < MAXABET; i++)
	{
	  if(cm->e[to_detach][i] >= 0.)
	    diff = cm->e[to_detach][i] - 0.;
	  else
	    diff = 0. - cm->e[to_detach][i];
	  if(diff > 0.000001)
	    cm_Fail("ERROR, to_detach state: %d e->[%d] is non-zero but rather %f\n", to_detach, i, cm->e[to_detach][i]);
	}
      for(yoffset = 0; yoffset < cm->cnum[to_detach]; yoffset++)
	{
	  /*printf("to_detach t[%d] cts: %f\n", yoffset, cm->t[to_detach][yoffset]);*/
	  if(cm->t[to_detach][yoffset] >= 0.)
	    diff = cm->t[to_detach][yoffset] - 0.;
	  else
	    diff = 0. - cm->t[to_detach][yoffset];
	  if(diff > 0.000001)
	    cm_Fail("ERROR, to_detach state: %d t->[%d] is non-zero but rather %f\n", to_detach, yoffset, cm->t[to_detach][yoffset]);
	}
    }
  return ret_val;
}

/* Functions: clean_cs()
 * Date:      SRE, Fri May 17 14:52:42 2002 [St. Louis]
 *
 * Purpose:   Verify and (if needed) clean the consensus structure annotation.
 */
int
clean_cs(char *cs, int alen, int be_quiet)
{
  int   status;
  int   i;
  int  *ct;
  int   nright = 0;
  int   nleft = 0;
  int   nbad = 0;
  char  example;
  int   first;
  int   has_pseudoknots = FALSE;

  /* 1. Check if we have a good CS line with >= 0 pseudoknotted
   *    base pairs.
   */
  ESL_ALLOC(ct, (alen+1) * sizeof(int));
  if (esl_wuss2ct(cs, alen, ct) != eslOK)  
    cm_Fail("Consensus structure string is inconsistent"); 
  free(ct);

  /* 2. CS line is good, check for and remove pseudoknots 
   *    if necessary. */
  if (check_for_pknots(cs, alen)) {
    has_pseudoknots = TRUE; 
    if(!be_quiet) printf("    [Consensus structure has annotated pseudoknots that will be ignored.]\n");
    fflush(stdout);
  }
  else return TRUE; /* we're good, no need to clean it, there's no 
		     * pseudoknots */

  /* 3. Delete everything we don't recognize.
   */
  for (i = 0; i < alen; i++)
    {
      if      (strchr("{[(<", cs[i]) != NULL) nleft++;  
      else if (strchr(">)]}", cs[i]) != NULL) nright++; 
      else if (strchr(":_-,.~", cs[i]) != NULL) ;
      else if (has_pseudoknots && isalpha((int) cs[i])) cs[i] = '.';
      else {	/* count bad chars; remember first one; replace w/gap */
	if (nbad == 0) { example = cs[i]; first = i; }
	nbad++;
	cs[i] = '.';
      }
    }
  if (nbad > 0) {
    if(!be_quiet) printf("    [Removed %d bad chars from consensus line. Example: a %c at position %d.]\n",
	   nbad, example, first);
    fflush(stdout);
  }

  /* Check it again.
   */
  ESL_ALLOC(ct, (alen+1) * sizeof(int));
  status = esl_wuss2ct(cs, alen, ct);  
  free(ct);
  if(status == eslOK) 
    return TRUE;
  printf("    [Failed to parse the consensus structure line.]\n"); 
  return FALSE;

 ERROR:
  cm_Fail("Memory allocation error.");
  return FALSE; /* never reached */
}

/* Functions: check_for_pknots()
 * Date:      EPN, Mon Aug  6 14:46:24 2007
 *
 * Purpose:   Simple check for pseudoknots in a consensus structure annotation.
 *            ASSUMES: CS has already been checked for consistency.
 */
static int
check_for_pknots(char *cs, int alen)
{
  int i;
  for (i = 0; i < alen; i++)
    {
      if (isalpha((int) cs[i]))
	return TRUE; /* assumes we know the CS is consistent */
    }
  return FALSE;
}



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
   */
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
      if (i > j) mtr->type[v] = END_nd;

      else if (mtr->type[v] == ROOT_nd)
	{ /* try to push i,j; but deal with INSL and INSR */
	  for (; i <= j; i++) if (matassign[i]) break;
	  for (; j >= i; j--) if (matassign[j]) break;
	  PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));
	}

      else if (mtr->type[v] == BEGL_nd) /* no inserts */
	PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));

      else if (mtr->type[v] == BEGR_nd) /* look for INSL */
	{
	  for (; i <= j; i++) if (matassign[i]) break;
	  PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));
	}

      else if (ct[i] == -1) 	/* i unpaired. This is a MATL node; allow INSL */
	{
	  mtr->type[v] = MATL_nd;
	  for (i = i+1; i <= j; i++)  if (matassign[i]) break;
	  PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, j, -1, DUMMY_nd));
	}

      else if (ct[j] == -1) 	/* j unpaired. MATR node. Deal with INSR */
	{
	  mtr->type[v] = MATR_nd;
	  for (nxtj = j-1; nxtj >= i; nxtj--) if (matassign[nxtj]) break;
	  PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, i, nxtj, -1, DUMMY_nd));
	}

      else if (ct[i] == j) 	/* i,j paired to each other. MATP. deal with INSL, INSR */
	{
	  mtr->type[v] = MATP_nd;
	  for (nxti = i+1; nxti <= j; nxti++)    if (matassign[nxti]) break;
	  for (nxtj = j-1; nxtj >= nxti; nxtj--) if (matassign[nxtj]) break;
	  PushNstack(pda, InsertTraceNode(mtr, v, TRACE_LEFT_CHILD, nxti, nxtj, -1,DUMMY_nd));
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
	}

    }	/* while something's on the stack */
  FreeNstack(pda);
  free(ct);

#if 0
  /* OK, we've converted ct into mtr -- mtr is a tree structure telling us the
   * arrangement of consensus nodes. Now do the drill for constructing a full model 
   * using this master trace. First find out how many states we need:
   */
  NumberMasterTrace(mtr, &nodes);
  if ((cm = AllocCM(nodes)) == NULL)
    Die("failed to allocate for new model of %d nodes\n", nodes);
  TopofyNewCM(cm, mtr);

  
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
  if (ret_cm  != NULL) *ret_cm  = cm;  else FreeCM(cm);
#endif /*0*/
  if (ret_mtr != NULL) *ret_mtr = mtr; else FreeParsetree(mtr);
  PrintParsetree(stdout, mtr);
  PrintParsetree(stdout, *ret_mtr);
}

/* emit.c
 * SRE, Mon Oct 13 12:36:35 2003 [St. Louis]
 * 
 * Sampling sequences or parse trees from a CM. 
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 * CVS $Id$
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

/* Function:  EmitParsetree()
 * Incept:    SRE, Mon Oct 13 22:35:46 2003 [Rams whupping Falcons, Monday Night Football]
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
 * Args:      cm      - covariance model to generate from
 *            ret_tr  - RETURN: generated parse tree. Pass NULL if unwanted.
 *            ret_seq - RETURN: generated sequence (alphabetic). Pass NULL if unwanted.
 *            ret_dsq - RETURN: generated sequence (digitized). Pass NULL if unwanted.
 *            ret_N   - RETURN: length of generated sequence.
 *
 * Returns:   (void)
 *            tr, seq, dsq are allocated here; whichever ones the caller
 *            requests (with non-NULL ret_ pointers) the caller is responsible
 *            for free'ing:
 *               FreeParsetree(tr); free(seq); free(dsq);
 */
void
EmitParsetree(CM_t *cm, Parsetree_t **ret_tr, char **ret_seq, char **ret_dsq, int *ret_N)
{
  Parsetree_t *tr;              /* parse tree under construction */
  Nstack_t *pda;                /* pushdown automaton for traversing parse tree */              
  Cstack_t *gsq;                /* growing sequence under construction */
  char     *seq;                /* finished sequence, normal alphabet form, 0..N-1 */
  char     *dsq;                /* digitized sequence, 1..N */
  int N;			/* current emitted sequence length */
  int tparent;			/* parent node index, last attached to parse tree */
  int tpos;			/* child node index, just attached to parse tree */
  int v;			/* index of current state */
  int y,z;			/* indices for next state(s)    */
  int type;			/* PDA_RESIDUE or PDA_STATE */
  int lchar, rchar;		/* index of emitted chars in Alphabet[], or -1 for nothing */
  int whichway;			/* how to attach: TRACE_LEFT_CHILD or TRACE_RIGHT_CHILD */
  int x;			/* tmp variable for sampling MP emission */

  tr  = CreateParsetree();
  pda = CreateNstack();
  gsq = CreateCstack();
  N   = 0;			

  /* Init by pushing root state's info onto pda
   */
  PushNstack(pda, -1);		/* doesn't emit an rchar */
  PushNstack(pda, -1);		/* doesn't emit an lchar either */
  PushNstack(pda, TRACE_LEFT_CHILD);
  PushNstack(pda, -1);		/* attach this state to parsetree node -1 (init) */  
  PushNstack(pda, 0);		/* it's the root state, v=0 */
  PushNstack(pda, PDA_STATE);

  /* Iterate until the pda is empty...
   */
  while (PopNstack(pda, &type)) 
    {
      if (type == PDA_RESIDUE)
	{
	  PopNstack(pda, &tpos);
	  PopNstack(pda, &rchar);

	  if (rchar != -1) {
	    PushCstack(gsq, Alphabet[rchar]);
	    N++;
	  }
	  tr->emitr[tpos] = N;
	}
      else if (type == PDA_STATE) 
	{
	  PopNstack(pda, &v);
	  PopNstack(pda, &tparent);
	  PopNstack(pda, &whichway);
	  PopNstack(pda, &lchar);
	  PopNstack(pda, &rchar);

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
	      PushCstack(gsq, Alphabet[lchar]);
	      N++;
	    }

	  /* Push right emission info for state v onto the pda, now
           * that we know tpos for where v is in the parsetree. We have
           * to do this even if rchar is -1, to be sure that we will set the emitr
           * bound properly even for nonemitting states in the parsetree.
	   */
	  PushNstack(pda, rchar);
	  PushNstack(pda, tpos);
	  PushNstack(pda, PDA_RESIDUE);

	  /* Decide what state we're going to next.
           * B is special case of a bifurcation to two S states. 
	   */
	  if (cm->sttype[v] == B_st)
	    {
	      y = cm->cfirst[v];	/* left child  */
	      z = cm->cnum[v];	        /* right child */
	  
	      /* Push the right start state's info
	       */
	      PushNstack(pda, -1);		/* doesn't emit right */
	      PushNstack(pda, -1);		/* doesn't emit left */
	      PushNstack(pda, TRACE_RIGHT_CHILD); /* attach as right child of the B */
	      PushNstack(pda, tpos);		/* attach it to B, which is tpos in parsetree*/
	      PushNstack(pda, z);		/* state z */
	      PushNstack(pda, PDA_STATE);

	      /* Push the left start state's info
	       */
	      PushNstack(pda, -1);		/* doesn't emit right */
	      PushNstack(pda, -1);		/* doesn't emit left */
	      PushNstack(pda, TRACE_LEFT_CHILD); /* attach as left child of the B */
	      PushNstack(pda, tpos);		/* attach it to B, which is tpos in parsetree*/
	      PushNstack(pda, y);		/* state z */
	      PushNstack(pda, PDA_STATE);
	    }
	  else			/* all other parent states v: */
	    {
	      y = cm->cfirst[v] + FChoose(cm->t[v], cm->cnum[v]); /* choose next state, y */

	      switch (cm->sttype[y]) {
	      case MP_st: 
		x     = FChoose(cm->e[y], Alphabet_size*Alphabet_size);
		lchar = x / Alphabet_size;
		rchar = x % Alphabet_size;
		break;
	      case ML_st:
	      case IL_st:
		lchar = FChoose(cm->e[y], Alphabet_size);
		rchar = -1;
		break;
	      case MR_st:
	      case IR_st:
		lchar = -1;
		rchar = FChoose(cm->e[y], Alphabet_size);
		break;
	      default:
		lchar = -1;
		rchar = -1;
	      }

	      if (cm->sttype[y] == E_st) {
		InsertTraceNode(tr, tpos, TRACE_LEFT_CHILD, -1, -1, y);
	      } else {
		PushNstack(pda, rchar);		/* does it emit right? */
		PushNstack(pda, lchar);		/* does it emit left? */
		PushNstack(pda, TRACE_LEFT_CHILD); /* non-B's: attach as left child by conv */
		PushNstack(pda, tpos);		/* attach it to v, which is tpos in parsetree*/
		PushNstack(pda, y);		/* next state we're going to */
		PushNstack(pda, PDA_STATE);
	      }
	    }
	} /* end of PDA_STATE logic */  
    } /* end of main "while PopNstack()" loop */

  seq = CstackString(gsq);
  dsq = DigitizeSequence(seq, N);
  FreeNstack(pda);

  if (ret_tr  != NULL) *ret_tr  = tr;  else FreeParsetree(tr);
  if (ret_seq != NULL) *ret_seq = seq; else free(seq);
  if (ret_dsq != NULL) *ret_dsq = dsq; else free(dsq);
  if (ret_N   != NULL) *ret_N   = N; 
}
  


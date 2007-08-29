/* emit.c
 * SRE, Mon Oct 13 12:36:35 2003 [St. Louis]
 * 
 * Sampling sequences or parse trees from a CM. 
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 * SVN $Id$
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_vectorops.h"
#include "esl_random.h"
#include "esl_stack.h"
#include "esl_sqio.h"

#include "structs.h"
#include "funcs.h"

/* Function:  EmitParsetree()
 * Incept:    SRE, Mon Oct 13 22:35:46 2003 [Rams whupping Falcons, Monday Night Football]
 *            Easel'ed: EPN, Fri Aug  3 08:15:12 2007
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
 *            Added capacity for local begins/ends. [EPN, Wed May  2 05:59:19 2007]
 *
 * Args:      cm      - covariance model to generate from
 *            r       - source of randomness
 *            name    - name for the sequence (ESL_SQ name field is mandatory)
 *            do_digital - TRUE to digitize sq before returning, FALSE not to
 *            ret_tr  - RETURN: generated parse tree. Pass NULL if unwanted.
 *            ret_sq  - RETURN: generated sequence
 *            ret_N   - RETURN: length of generated sequence.
 *
 * Returns:   eslOK on success; eslEMEM on memory error;
 *            tr, sq are allocated here; whichever ones the caller
 *            requests (with non-NULL ret_ pointers) the caller is responsible
 *            for free'ing:
 *               FreeParsetree(tr); esl_sq_Destroy(sq);
 */
int
EmitParsetree(CM_t *cm, ESL_RANDOMNESS *r, char *name, int do_digital, Parsetree_t **ret_tr, ESL_SQ **ret_sq, int *ret_N)
{
  int status;
  Parsetree_t *tr = NULL;       /* parse tree under construction */
  ESL_STACK *pda = NULL;        /* pushdown automaton for traversing parse tree */              
  ESL_STACK *gsq = NULL;        /* growing sequence under construction */
  ESL_SQ    *sq  = NULL;        /* finished sequence, initially normal alphabet form */
  char      *seq;               /* alphabetic sequence to build sq with */
  int N;			/* current emitted sequence length */
  int tparent;			/* parent node index, last attached to parse tree */
  int tpos;			/* child node index, just attached to parse tree */
  int v;			/* index of current state */
  int y,z;			/* indices for next state(s)    */
  int type;			/* PDA_RESIDUE or PDA_STATE */
  int lchar, rchar;		/* index of emitted chars in cm->abc->sym[], or -1 for nothing */
  int whichway;			/* how to attach: TRACE_LEFT_CHILD or TRACE_RIGHT_CHILD */
  int x;			/* tmp variable for sampling MP emission */
  int lpos;                     /* tmp variable for inserting EL trace node */
  float *tmp_tvec = NULL;       /* tmp transition vector to choose from, 
				 * for dealing with local end transitions */
  /* Contract check */
  if(cm->flags & CM_LOCAL_END && (fabs(sreEXP2(cm->el_selfsc) - 1.0) < 0.01))
    ESL_XEXCEPTION(eslEINVAL, "EL self transition probability %f is too high, would emit long (too long) EL insertions.", sreEXP2(cm->el_selfsc));
  if(cm->abc == NULL)
    ESL_XEXCEPTION(eslEINVAL, "CM does not have a valid alphabet.");
  if(ret_sq != NULL && name == NULL)
    ESL_XEXCEPTION(eslEINVAL, "EmitParsetree requires a sequence name for the sequence it's creating.");

  tr  = CreateParsetree(100);
  pda = esl_stack_ICreate();
  gsq = esl_stack_CCreate();
  N   = 0;			
  ESL_ALLOC(tmp_tvec, sizeof(float) * (MAXCONNECT+1)); /* enough room for max possible transitions, plus
							* a local end transition */
  /* Init by pushing root state's info onto pda
   */
  esl_stack_IPush(pda, -1);		/* doesn't emit an rchar */
  esl_stack_IPush(pda, -1);		/* doesn't emit an lchar either */
  esl_stack_IPush(pda, TRACE_LEFT_CHILD);
  esl_stack_IPush(pda, -1);		/* attach this state to parsetree node -1 (init) */  
  esl_stack_IPush(pda, 0);		/* it's the root state, v=0 */
  esl_stack_IPush(pda, PDA_STATE);

  /* Iterate until the pda is empty...
   */
  while (esl_stack_IPop(pda, &type) != eslEOD) 
    {
      if (type == PDA_RESIDUE)
	{
	  esl_stack_IPop(pda, &tpos);
	  esl_stack_IPop(pda, &rchar);

	  if (rchar != -1) {
	    esl_stack_CPush(gsq, cm->abc->sym[rchar]);
	    N++;
	  }
	  tr->emitr[tpos] = N;
	}
      else if (type == PDA_STATE) 
	{
	  esl_stack_IPop(pda, &v);
	  esl_stack_IPop(pda, &tparent);
	  esl_stack_IPop(pda, &whichway);
	  esl_stack_IPop(pda, &lchar);
	  esl_stack_IPop(pda, &rchar);

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
	      esl_stack_CPush(gsq, cm->abc->sym[lchar]);
	      N++;
	    }

	  /* Push right emission info for state v onto the pda, now
           * that we know tpos for where v is in the parsetree. We have
           * to do this even if rchar is -1, to be sure that we will set the emitr
           * bound properly even for nonemitting states in the parsetree.
	   */
	  esl_stack_IPush(pda, rchar);
	  esl_stack_IPush(pda, tpos);
	  esl_stack_IPush(pda, PDA_RESIDUE);

	  /* Decide what state we're going to next.
           * B is special case of a bifurcation to two S states. 
	   */
	  if (cm->sttype[v] == B_st)
	    {
	      y = cm->cfirst[v];	/* left child  */
	      z = cm->cnum[v];	        /* right child */
	  
	      /* Push the right start state's info
	       */
	      esl_stack_IPush(pda, -1);		/* doesn't emit right */
	      esl_stack_IPush(pda, -1);		/* doesn't emit left */
	      esl_stack_IPush(pda, TRACE_RIGHT_CHILD); /* attach as right child of the B */
	      esl_stack_IPush(pda, tpos);		/* attach it to B, which is tpos in parsetree*/
	      esl_stack_IPush(pda, z);		/* state z */
	      esl_stack_IPush(pda, PDA_STATE);

	      /* Push the left start state's info
	       */
	      esl_stack_IPush(pda, -1);		/* doesn't emit right */
	      esl_stack_IPush(pda, -1);		/* doesn't emit left */
	      esl_stack_IPush(pda, TRACE_LEFT_CHILD); /* attach as left child of the B */
	      esl_stack_IPush(pda, tpos);		/* attach it to B, which is tpos in parsetree*/
	      esl_stack_IPush(pda, y);		/* state z */
	      esl_stack_IPush(pda, PDA_STATE);
	    }
	  else
	    {
	      if(v == 0 && cm->flags & CM_LOCAL_BEGIN)		/* ROOT_S with local begins, special */
		y = esl_rnd_FChoose(r, cm->begin, cm->M); /* choose next state, y */
	      else if(cm->flags & CM_LOCAL_END) /* special case, we could transit to EL */
		{
		  esl_vec_FSet(tmp_tvec, (MAXCONNECT+1), 0.);
		  esl_vec_FCopy(cm->t[v], cm->cnum[v], tmp_tvec);
		  tmp_tvec[cm->cnum[v]] = cm->end[v];
		  y = esl_rnd_FChoose(r, tmp_tvec, (cm->cnum[v]+1)); /* choose next state, y's offset */
		  if(y == cm->cnum[v]) y = cm->M; /* local end */
		  else y += cm->cfirst[v];        
		}		  
	      else
		y = cm->cfirst[v] + esl_rnd_FChoose(r, cm->t[v], cm->cnum[v]); /* choose next state, y */

	      switch (cm->sttype[y]) {
	      case MP_st: 
		x     = esl_rnd_FChoose(r, cm->e[y], cm->abc->K*cm->abc->K);
		lchar = x / cm->abc->K;
		rchar = x % cm->abc->K;
		break;
	      case ML_st:
	      case IL_st:
		lchar = esl_rnd_FChoose(r, cm->e[y], cm->abc->K);
		rchar = -1;
		break;
	      case MR_st:
	      case IR_st:
		lchar = -1;
		rchar = esl_rnd_FChoose(r, cm->e[y], cm->abc->K);
		break;
	      case EL_st: /* EL emits on transition, here we don't emit */
		lchar = -1;
		rchar = -1;
		break;
	      default:
		lchar = -1;
		rchar = -1;
	      }
	      if (cm->sttype[y] == E_st)
		{
		  /*InsertTraceNode(tr, tpos, TRACE_LEFT_CHILD, -1, -1, y);*/
		  InsertTraceNode(tr, tpos, TRACE_LEFT_CHILD, N+1, N, y);
		} 
	      else if(cm->sttype[y] == EL_st) /* y == cm->M */
		{
		  lpos = N+1; /* remember lpos, we need it after we emit from EL */
		  /* Now choose number of residues emitted from EL, could be 0.
		   * We do this here b/c convention for EL is to have a single trace node,
		   * even if multiple residues are emitted. */
		  esl_vec_FSet(tmp_tvec, (MAXCONNECT+1), 0.);
		  tmp_tvec[0] = sreEXP2(cm->el_selfsc); /* EL self probability */
		  tmp_tvec[1] = 1. - tmp_tvec[0];       /* probability of going to implicit END */
		  y = esl_rnd_FChoose(r, tmp_tvec, 2); /* choose next state, either EL or implicit END */
		  while(y == 0) /* we've self-transitioned, emit 1 res from NULL distro */
		    {
		      lchar = esl_rnd_FChoose(r, cm->null, cm->abc->K);
		      esl_stack_CPush(gsq, cm->abc->sym[lchar]);
		      N++;
		      y = esl_rnd_FChoose(r, tmp_tvec, 2); /* choose next state, either EL or implicit END */
		    }
		  InsertTraceNode(tr, tpos, TRACE_LEFT_CHILD, lpos, N, cm->M); /* careful to reset y to cm->M */
		}
	      else 
		{
		  esl_stack_IPush(pda, rchar);		/* does it emit right? */
		  esl_stack_IPush(pda, lchar);		/* does it emit left? */
		  esl_stack_IPush(pda, TRACE_LEFT_CHILD); /* non-B's: attach as left child by conv */
		  esl_stack_IPush(pda, tpos);		/* attach it to v, which is tpos in parsetree*/
		  esl_stack_IPush(pda, y);		/* next state we're going to */
		  esl_stack_IPush(pda, PDA_STATE);
		}
	    } /* end of PDA_STATE logic */  
	} /* end of else (which we enter if v not a B state) */
    } /* end of main "while esl_stack_IPop()" loop */

  seq = esl_stack_Convert2String(gsq); /* this destroys gsq char stack */
  if(name != NULL) sq  = esl_sq_CreateFrom(name, seq, NULL, NULL, NULL);
  else             sq  = esl_sq_CreateFrom("seq", seq, NULL, NULL, NULL); 
  free(seq); /* we made a copy of this when creating sq */
  /* name can only be NULL if ret_sq == NULL, so we're throwing it away anyway */

  /* digitize if nec */
  if(do_digital) 
    if((status = esl_sq_Digitize(cm->abc, sq)) != eslOK) goto ERROR;
  /*ParsetreeDump(stdout, tr, cm, dsq);*/ 

  esl_stack_Destroy(pda);

  free(tmp_tvec);
  if (ret_tr  != NULL) *ret_tr  = tr;  else FreeParsetree(tr);
  if (ret_sq  != NULL) *ret_sq  = sq;  else esl_sq_Destroy(sq);
  if (ret_N   != NULL) *ret_N   = N; 
  return eslOK;
  
 ERROR:
  if(tr  != NULL) FreeParsetree(tr);
  if(gsq != NULL) esl_stack_Destroy(gsq);
  if(pda != NULL) esl_stack_Destroy(pda);
  if(sq  != NULL) esl_sq_Destroy(sq);
  if(tmp_tvec != NULL) free(tmp_tvec);
  return status;
}
  


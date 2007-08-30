/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_masks.c
 * EPN, 08.24.06 Janelia 
 * based on HMMER 2.3.2's masks.c, ported to CMs for the TraceScoreCorrection()
 * function.
 * Most original SRE comments untouched.
 * HMMER's masks.c: SRE, Tue Nov 18 10:12:28 1997
 * 
 * Corrections for biased composition
 * target sequences. 
 * 
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/* Function: CM_TraceScoreCorrection()
 * based on     TraceScoreCorrection() from HMMER:
 * EPN 08.24.06 Janelia
 * 
 * Purpose:  Calculate a correction (in integer log_2 odds) to be
 *           applied to a sequence, using a second null model, 
 *           based on a traceback. All emissions are corrected;
 *           The null model is constructed /post hoc/ as the
 *           average over all the emission distributions used by the trace.
 *           
 * Return:   the log_2-odds score correction.          
 */
float
CM_TraceScoreCorrection(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq)
{
  float p[MAXABET];		/* null2 model distribution */
  float sc[MAXDEGEN];		/* null model scores       */
  int   a,b;
  int   v;                      /* state index counter */
  int   i, j;                   /* seq posn counter */
  int   tidx;
  float score;

  /* Rarely, the alignment was totally impossible, and tr is NULL.
   */
  if (tr == NULL) return 0.0;
  
  /* Set up model: average over the emission distributions of
   * all M, I states that appear in the trace. Ad hoc? Sure, you betcha. 
   */
		/* trivial preorder traverse, since we're already numbered that way */
  esl_vec_FSet(p, MAXABET, 0.0);
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    if(cm->sttype[v] == MP_st)
      {
	/* we treat this as two match states. */
	for(a = 0; a < MAXABET; a++)
	  {
	    /* first add contribution to null2 for left half. 
	       There's probably a simpler (vectorops.c func call) way to do this. */
	    for(b = (a * MAXABET); b < ((a+1) * MAXABET); b++)
	      p[a] += cm->e[v][b]; 
	    /* now add contribution for right half. */
	    for(b = a; b < (MAXABET * MAXABET); b += MAXABET)
	      p[a] += cm->e[v][b]; 
	  }
      }
    else if(cm->sttype[v] == ML_st || cm->sttype[v] == IL_st ||
	    cm->sttype[v] == MR_st || cm->sttype[v] == IR_st)
      esl_vec_FAdd(p, cm->e[v], MAXABET);
  }

  esl_vec_FNorm(p, MAXABET);

  for (a = 0; a < MAXABET; a++)
    sc[a] = sreLOG2(p[a] / cm->null[a]);
				/* could avoid this chunk if we knew
				   we didn't need any degenerate char scores */
  for (a = MAXABET; a < MAXDEGEN; a++)
    sc[a] = esl_abc_FAvgScore(cm->abc, a, p);  /* not completely sure about this line EPN */

  /* Score all the state emissions that appear in the trace.
   */
  score = 0;
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    i = tr->emitl[tidx];
    j = tr->emitr[tidx];
    if (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || 
	cm->sttype[v] == IL_st)
      score += sc[dsq[i]];
    if (cm->sttype[v] == MP_st || cm->sttype[v] == MR_st || 
	cm->sttype[v] == IR_st)
      score += sc[dsq[j]];
  }
   /* Apply an ad hoc 8 bit fudge factor penalty;
    * interpreted as a prior, saying that the second null model is 
    * 1/2^8 (1/256) as likely as the standard null model
    */
   score -= 8.;	

   /* Return the correction to the bit score.
    */
   /*printf("returning from CM_TraceScoreCorrection: %f\n", LogSum2(0,score));*/
   return LogSum2(0,score);
}

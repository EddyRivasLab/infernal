/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_eweight.c [EPN 11.07.05]
 * based on: HMMER 2.4devl's lsj_eweight.c
 * Most original comments from lsj_eweight.c untouched.
 * 
 * LSJ, Wed Feb  4 15:03:58 CST 2004
 * 
 * entropy targeting:
 * Code for setting effective sequence number (in hmmbuild) by
 * achieving a certain target entropy loss, relative to background
 * null distribution.
 *
 * CVS $Id: lsj_eweight.c 950 2004-06-16 14:46:04Z eddy $
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <string.h>

#include "prior.h"

#include "structs.h"
#include "funcs.h"
/*#include "squid.h"*/
/*#include "vectorops.h"*/
#include <esl_vectorops.h>
#include "cm_eweight.h"

/* Function: CM_Eweight [EPN]
 * based on:
 * Eweight() LSJ 2/6/04
 * 
 * Purpose:  Main entropy-based weighting function. 
 *           
 * Args:  
 *              cm       - the model
 *           **pri       - Model priors.
 *       numb_seqs       - Number of sequences in alignment.
 *       targetent       - Target mean match state entropy. 
 *           
 * Return: eff_no        - New effective sequence number.                         
 */
double
CM_Eweight(CM_t *cm, Prior_t *pri, float numb_seqs, 
	float targetent)
{
  int i;
  int j;
  float eff_no;                  /* New effective sequence number */
  double current;                /* Current mean match state entropy */
  double prevent;                 /* Previous mean match state entropy */
  float scale;                   /* Current model counts scaling factor */
  float leftscale;               /* Bracket scaling value used in binary search. Lowest mean entropy value. */
  float rightscale;              /* Bracket scaling value used in binary search. Highest mean entropy value. */

  double *ent;                    /* Match state entropy values */
  int count;                     /* Counter for binary search */
  int flag;                      /* Used to detect entropy adjustment failure */

  int nmatch_cols;               /* num MATL_nd + MATR_nd + 2 * MATP_nd in CM */
  
  /* analags of parameters from Infernal's prior.c()'s PriorifyCM().*/
  double *counts;                 /* Temp array of match state counts */
  double *probs;                  /* Temp array of match state probs */
  double *mixq;                   /* posterior probs of mixture components, P(q | c) */


  /**************
   * Allocations
   **************/
  ent    = MallocOrDie((cm->nodes) * sizeof(double));
  counts = MallocOrDie(sizeof(double) * pri->maxnalpha);
  probs  = MallocOrDie(sizeof(double) * pri->maxnalpha);
  mixq   = MallocOrDie(sizeof(double) * pri->maxnq);
	  	  
  /*****************
   * Initializations 
   *****************/
  current  = 0.;
  scale    = 1.;
  count    = 0;
  flag     = 0;
  nmatch_cols = 0;

  for(i = 0; i < cm->nodes; i++)
    ent[i] = 0.;

  /***************************************
   * Calculate the starting model entropy 
   ***************************************/

  /* Copy model match state probabilities into our temporary counts[]
   * (Current implementation only considers MATP_MP as a match state,
   *  for MATP nodes, not MATP_ML or MATP_MR).
   * For nodes i with no match state (BEGL, BEGR, ROOT, BIF and END)
   * ent[i] is left as its initialized value; 0.0. This effectively
   * eliminates any contribution to 'current' from such nodes.
   * Remember our CM is still in counts form, so cm->e[][] is a count
   * not a probability. 
   */
  for(i = 0; i < cm->nodes; i++)
    { 
      if(cm->ndtype[i] == MATP_nd)
	{
	  nmatch_cols += 2; /* two match columns */
	  for(j = 0; j < (MAXABET * MAXABET); j++)
	    counts[j] = cm->e[cm->nodemap[i]][j];
	  /* cm->nodemap[i] = first state, node i (here, MP state) */

	  /* Add priors to the current match state prob dist. (easel/esl_dirichlet.c) */
	  esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
				     pri->mbp,
				     mixq, probs);
	  /* ent[] is assigned the current MP_st state emission entropy. */
	  ent[i] = esl_vec_DEntropy(probs, (MAXABET * MAXABET));
	}
      else if ((cm->ndtype[i] == MATL_nd) ||
	       (cm->ndtype[i] == MATR_nd))
	{
	  nmatch_cols++;
	  for(j = 0; j < MAXABET; j++)
	    counts[j] = cm->e[cm->nodemap[i]][j];
	  /* cm->nodemap[i] = first state, node i (here, ML or MR state) */

	  /* Add priors to the current match state prob dist. (easel/esl_dirichlet.c) */
	  esl_mixdchlet_MPParameters(counts, MAXABET,
				     pri->mnt,
				     mixq, probs);
	  /* ent[] is assigned the current consensus singlet emission entropy. */
	  ent[i] = esl_vec_DEntropy(probs, MAXABET);
	}
      /* other nodes are skipped, ent[i] for these nodes remains 0.0 */
    }
  /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
  current = esl_vec_DSum(ent, cm->nodes)/nmatch_cols;
  /*printf("target ent: %f\n", targetent);*/
  /*printf("0 current: %f\n", current);*/

  /****************************************
   * Initialize binary search bracket values
   *****************************************/

  /* The reason the values seem backwards is because I'm trying to
     bracket my target mean entropy with model count scaling
     factors. A higher scaling factor generally produces a lower
     Entropy and a lower scaling factor produces a higher
     entropy. Thus, the leftscale produces the lowest mean entropy
     bracket and rightscale produces the highest mean entropy
     bracket */
  if(current < targetent){
    leftscale  = 1; 
    rightscale = 0; 
  } 
  else{
    /* Current model has a higher entropy than our target.
       Calculated effective seq numb <= Number of seqs. Design decision.
    */
    printf("[scale=%.2f] [e=%.2f >= %.2f] ...", scale, current, targetent);
    free(mixq);
    free(counts);
    free(probs);
    free(ent);
    return(numb_seqs);
  }
  /***************************************
   * Binary search for target mean entropy
   ***************************************/
  /* Check to see if the current model mean entropy is within 0.01 bits of our target */
  while((current < targetent - 0.01) || (current > targetent + 0.01))
    {
      count++;
      nmatch_cols = 0;
    
    /* Emergency brake in case there is a bug in our binary search */
    if(count > 50){
      printf("\nBUG: Problem with adjusting the model entropy. Please report.\n");
      break;
    }

    /* Calculate current scaling factor based on bracket values */
    scale = (leftscale + rightscale)/2;
    
    prevent = current;

    /*******************************************
     * Scale the counts and re-calc the entropy
     *******************************************/
    /* Re-copy match state probabilities into counts[] */
    for(i = 0; i < cm->nodes; i++)
      { 
	if(cm->ndtype[i] == MATP_nd)
	  {
	    nmatch_cols += 2; /* two match columns */
	    for(j = 0; j < (MAXABET * MAXABET); j++)
	      counts[j] = cm->e[cm->nodemap[i]][j];
	    /* cm->nodemap[i] = first state, node i (here, MP state) */

	    /* Re-scale the current counts by the previously determined amount. 
	     * (easel/esl_vectorops.c) 
	     */
	    DScale(counts, (MAXABET*MAXABET), scale);

	    /* Re-add priors to these scaled counts. (easel/esl_dirichlet.c) */
	    esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
				       pri->mbp,
				       mixq, probs);
	    /* Again, ent[] is assigned the current match emission entropy */
	    ent[i] = esl_vec_DEntropy(probs, (MAXABET * MAXABET));
	  }
	else if ((cm->ndtype[i] == MATL_nd) ||
		 (cm->ndtype[i] == MATR_nd))
	  {
	    nmatch_cols++;
	    for(j = 0; j < MAXABET; j++)
	      counts[j] = cm->e[cm->nodemap[i]][j];
	    /* cm->nodemap[i] = first state, node i (here, ML or MR state) */

	    /* Re-scale the current counts by the previously determined amount. 
	     * (easel/esl_vectorops.c) 
	     */
	    DScale(counts, MAXABET, scale);

	    /* Re-add the priors to these scaled counts. (easel/esl_dirichlet.c) */
	    esl_mixdchlet_MPParameters(counts, MAXABET,
				       pri->mnt,
				       mixq, probs);

	    /* Again, ent[] is assigned the current match emission entropy */
	    ent[i] = esl_vec_DEntropy(probs, MAXABET);
	  }
	/* other nodes are skipped, ent[i] for these nodes remains 0.0 */
      }
    /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
    current = esl_vec_DSum(ent, cm->nodes)/nmatch_cols;
    /*    printf("current : %f\n", current);*/

    /* Adjust the brackets according to the new mean entropy value */
    if(current < targetent){
      leftscale = scale;
    }
    else{
      /* We overshot the target. Replace right bracket with the current scale */
      rightscale = scale;
    }
  }
    free(mixq);
    free(counts);
    free(probs);
    free(ent);
  /**********************************************************************************************
   * End of binary search
   *********************************************************************************************/
  eff_no = numb_seqs * scale;
  printf("[scale=%.2f] returning eff_no : %f\n", scale, eff_no);
  return(eff_no);
}



/************************************************/
/* Functions just used in debugging/calibrating */ 
/************************************************/

/* Function: ModelContent() LSJ 10/14/03
 * 
 * Purpose:  This is a highly mutable grab-bag function I use  
 *           in benchmarking/debugging to examine model guts.
 *           
 * Args:     
 *           *ent1       - Column entropies for count data.
 *           *ent2       - Column entropies for count+prior data.
 *           M           - number of states in model
 *           
 * Return:   (void)                         
 */
void ModelContent(float *ent1, float *ent2, int M)
{
  int i;
  float sum1, sum2, sum3;
  float mean1, mean2, mean3;

  sum1  = 0;
  sum2  = 0;
  sum3  = 0;
  mean1 = 0;
  mean2 = 0;
  mean3 = 0;

  for(i = 1; i < M+1; i++){
    sum1 += ent1[i];
    sum2 += ent2[i];
    /*    sum3 += relent[i];
     */
    printf("%d\t%2.4f %2.4f %2.4f\n", i, ent1[i], ent2[i], (ent2[i] - ent1[i]));
  }
  mean1 = sum1/M;
  mean2 = sum2/M;
  /*  mean3 = sum3/M;
  fprintf(fp, "Mean Relative Entropy/Column: %2.4f\n", mean3);
  */
  printf("Counts Mean Entropy/Column: %2.4f\n", mean1);
  printf("Counts+Priors Mean Entropy/Column: %2.4f\n", mean2);
  printf("Diff: %2.4f\n", (mean2-mean1));
}

/* Function:  CMRescale() 
 *            (should this go in cm.c?)
 * Incept:    EPN 11.07.05
 * based on:  HMMER's plan7.c's Plan7Rescale() (Steve Johnson)
 *
 * Purpose:   Scale a counts-based CM by some factor, for
 *            adjusting counts to a new effective sequence number.
 *
 * Args:      cm         - counts based CM.
 *            scale      - scaling factor (e.g. eff_nseq/nseq); 1.0= no scaling.
 *
 * Returns:   (void)
 */
void 
CMRescale(CM_t *cm, float scale)
{
  int v;

  for (v = 0; v < cm->M; v++)
    {
      /* Scale transition counts vector if not a BIF or E state */
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	{
	  /* Number of transitions is cm->cnum[v] */
	  esl_vec_FScale(cm->t[v], cm->cnum[v], scale);
	}
      /* Scale emission counts vectors */
      if (cm->sttype[v] == MP_st)
	{       /* Consensus base pairs */
	  esl_vec_FScale(cm->e[v], (MAXABET*MAXABET), scale);
	}
      else if ((cm->sttype[v] == ML_st) ||
	       (cm->sttype[v] == MR_st) ||
	       (cm->sttype[v] == IL_st) ||
	       (cm->sttype[v] == IR_st))
	{      /* singlets (some consensus, some not)*/
	  esl_vec_FScale(cm->e[v], MAXABET, scale);
	}
    }/* end loop over states v */

  /* begin, end transitions; only valid [0..M-1] */
  FScale(cm->begin, cm->M, scale);
  FScale(cm->end,   cm->M, scale);
  
  return;
}

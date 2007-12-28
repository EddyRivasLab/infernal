/************************************************************
 * @LICENSE@
 ************************************************************/

/* eweight.c [EPN 11.07.05]
 * based on: HMMER 2.4devl's lsj_eweight.c
 * Most original comments from lsj_eweight.c untouched.
 * 
 * LSJ, Wed Feb  4 15:03:58 CST 2004
 * 
 * entropy targeting:
 * Code for setting effective sequence number (in cmbuild) by
 * achieving a certain target entropy loss, relative to background
 * null distribution.
 *
 * SVN $Id$
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

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
CM_Eweight(CM_t *cm, const Prior_t *pri, float numb_seqs, 
	float targetent)
{
  int status;
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
  ESL_ALLOC(ent,      sizeof(double) * cm->nodes);
  ESL_ALLOC(counts,   sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(probs,    sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(mixq,     sizeof(double) * pri->maxnq);
	  	  
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
   *  for MATP nodes, not MATP_ML or MATP_MR (MATL_ML and MATR_MR are
   *  also considered match states)).
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
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						  pri->mbp,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
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
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						  pri->mnt,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
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
    
    /* Emergency brake in case there is a bug in our binary search.
     * Its more likely that the target entropy is unattainable. */
      if(count > 50){
	printf("\nThe requested target entropy of %f is unattainable. [scale=%.2f] \n", targetent, scale);
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
	      esl_vec_DScale(counts, (MAXABET*MAXABET), scale);
	      
	      /* Re-add priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						      pri->mbp,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
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
	      esl_vec_DScale(counts, MAXABET, scale);
	      
	      /* Re-add the priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						      pri->mnt,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	      
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
  /*printf("[scale=%.2f] ", scale);*/
  return(eff_no);

 ERROR: 
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
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
  esl_vec_FScale(cm->begin, cm->M, scale);
  esl_vec_FScale(cm->end,   cm->M, scale);
  
  return;
}


/* Function:  CPlan9Rescale() 
 *            EPN based on Steve Johnsons plan 7 version
 *
 * Purpose:   Scale a counts-based HMM by some factor, for
 *            adjusting counts to a new effective sequence number.
 *
 * Args:      hmm        - counts based HMM.
 *            scale      - scaling factor (e.g. eff_nseq/nseq); 1.0= no scaling.
 *
 * Returns:   (void)
 */
void 
CPlan9Rescale(CP9_t *hmm, float scale)
{
  int k;

  /* emissions and transitions in the main model.
   * Note that match states are 1..M, insert states are 0..M,
   * and deletes are 0..M-1
   */
  for(k = 1; k <= hmm->M; k++) 
    esl_vec_FScale(hmm->mat[k], hmm->abc->K, scale);
  for(k = 0; k <=  hmm->M; k++) 
    esl_vec_FScale(hmm->ins[k], hmm->abc->K, scale);
  for(k = 0; k <  hmm->M; k++) 
    esl_vec_FScale(hmm->t[k],   cp9_NTRANS,             scale);

  /* begin, end transitions; only valid [1..M] */
  esl_vec_FScale(hmm->begin+1, hmm->M, scale);
  esl_vec_FScale(hmm->end+1,   hmm->M, scale);
  
  return;
}

/* Function: CM_Eweight_RE [EPN]
 * based on:
 * Eweight() LSJ 2/6/04
 * 
 * Purpose:  Main entropy-based weighting function. Calculates
 *           relative entropy (RE) instead of entropy. Requires background
 *           distribution. 
 *           
 * Args:  
 *              cm       - the model
 *           **pri       - Model priors.
 *       numb_seqs       - Number of sequences in alignment.
 *     target_relent     - Target mean match state relative entropy. 
 * randomseq[MAXABET]    - null sequence model
 * 
 * Return: eff_no        - New effective sequence number.                         
 */
double
CM_Eweight_RE(CM_t *cm, const Prior_t *pri, float numb_seqs, 
	      float target_relent, float *randomseq)
{
  int status;
  int i;
  int j;
  float eff_no;                  /* New effective sequence number */
  double current;                /* Current mean match state entropy */
  double prevent;                 /* Previous mean match state entropy */
  float scale;                   /* Current model counts scaling factor */
  float leftscale;               /* Bracket scaling value used in binary search. Lowest mean entropy value. */
  float rightscale;              /* Bracket scaling value used in binary search. Highest mean entropy value. */

  double *rel_ent;                    /* Match state relative entropy values */
  int count;                     /* Counter for binary search */
  int flag;                      /* Used to detect entropy adjustment failure */

  int nmatch_cols;               /* num MATL_nd + MATR_nd + 2 * MATP_nd in CM */
  
  /* analags of parameters from Infernal's prior.c()'s PriorifyCM().*/
  double *counts;                 /* Temp array of match state counts */
  double *probs;                  /* Temp array of match state probs */
  double *mixq;                   /* posterior probs of mixture components, P(q | c) */
  double Drandomseq[MAXABET];    /* the randomseq background prob dist, in doubles*/
  double Drandomseq_bp[MAXABET*MAXABET]; /* the randomseq BP background 
					    prob dist, in doubles*/

  /**************
   * Allocations
   **************/
  ESL_ALLOC(rel_ent, sizeof(double) * (cm->nodes));
  ESL_ALLOC(counts,  sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(probs,   sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(mixq,    sizeof(double) * pri->maxnq);
	  	  
  /*****************
   * Initializations 
   *****************/
  current  = 0.;
  scale    = 1.;
  count    = 0;
  flag     = 0;
  nmatch_cols = 0;

  for(i = 0; i < cm->nodes; i++)
    rel_ent[i] = 0.;

  for(i = 0; i < MAXABET; i++)
    Drandomseq[i] = (double) randomseq[i];
  
  for(i = 0; i < MAXABET; i++)
    for(j = 0; j < MAXABET; j++)
      Drandomseq_bp[i*MAXABET+j] = Drandomseq[i] * Drandomseq[j];

  /***************************************
   * Calculate the starting model entropy 
   ***************************************/

  /* Copy model match state probabilities into our temporary counts[]
   * (Current implementation only considers MATP_MP as a match state,
   *  for MATP nodes, not MATP_ML or MATP_MR (MATL_ML and MATR_MR are
   *  also considered match states)).
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
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						  pri->mbp,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	  /* rel_ent[] is assigned the current MP_st state emission relative
	   * entropy. */
	  rel_ent[i] = DRelEntropy(probs, Drandomseq_bp, (MAXABET * MAXABET));
	}
      else if ((cm->ndtype[i] == MATL_nd) ||
	       (cm->ndtype[i] == MATR_nd))
	{
	  nmatch_cols++;
	  for(j = 0; j < MAXABET; j++)
	    counts[j] = cm->e[cm->nodemap[i]][j];
	  /* cm->nodemap[i] = first state, node i (here, ML or MR state) */

	  /* Add priors to the current match state prob dist. (easel/esl_dirichlet.c) */
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						  pri->mnt,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	  /* rel_ent[] is assigned the current consensus singlet emission 
	     relative entropy. */
	  rel_ent[i] = DRelEntropy(probs, Drandomseq, MAXABET);
	  /*printf("rel_ent[%d] : %f\n", i, rel_ent[i]);*/
	}
      /* other nodes are skipped, ent[i] for these nodes remains 0.0 */
    }
  /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
  current = esl_vec_DSum(rel_ent, cm->nodes)/nmatch_cols;
  /*printf("target rel ent: %f\n", target_relent);*/
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
  if(current > target_relent){
    leftscale  = 1; 
    rightscale = 0; 
  } 
  else{
    /* Current model has a lower relative entropy than our target.
       Calculated effective seq numb <= Number of seqs. Design decision.
    */
    /*printf("[scale=%.2f] [re=%.2f >= %.2f] ...", scale, current, target_relent);*/
    free(mixq);
    free(counts);
    free(probs);
    free(rel_ent);
    return(numb_seqs);
  }
  /***************************************
   * Binary search for target mean entropy
   ***************************************/
  /* Check to see if the current model mean entropy is within 0.01 bits of our target */
  while((current < target_relent - 0.01) || (current > target_relent + 0.01))
    {
      count++;
      nmatch_cols = 0;
    
    /* Emergency brake in case there is a bug in our binary search.
     * Its more likely that the target entropy is unattainable. */
      if(count > 50){
	/*printf("\nThe requested target relative entropy of %f is unattainable. [scale=%.2f] \n", target_relent, scale);*/
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
	      esl_vec_DScale(counts, (MAXABET*MAXABET), scale);
	      
	      /* Re-add priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						      pri->mbp,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	      /* Again, rel_ent[] is assigned the current match emission 
		 relative entropy */
	      rel_ent[i] = DRelEntropy(probs, Drandomseq_bp, (MAXABET * MAXABET));
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
	      esl_vec_DScale(counts, MAXABET, scale);
	      
	      /* Re-add the priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						      pri->mnt,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	      
	      /* rel_ent[] is assigned the current consensus singlet emission 
		 relative entropy. */
	      rel_ent[i] = DRelEntropy(probs, Drandomseq, MAXABET);
	    }
	  /* other nodes are skipped, ent[i] for these nodes remains 0.0 */
	}
      /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
      current = esl_vec_DSum(rel_ent, cm->nodes)/nmatch_cols;
      /*printf("current : %f\n", current);*/
      
      /* Adjust the brackets according to the new mean entropy value */
      if(current > target_relent){
	leftscale = scale;
      }
      else{
	/* Replace right bracket with the current scale */
	rightscale = scale;
      }
    }
  free(mixq);
  free(counts);
  free(probs);
  free(rel_ent);
  /**********************************************************************************************
   * End of binary search
   *********************************************************************************************/
  eff_no = numb_seqs * scale;
  /*printf("[scale=%.2f] ", scale);*/
  return(eff_no);

 ERROR: 
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* Function:  DRelEntropy()
 *
 * Purpose:   Returns the relative entropy (KL distance) 
 *            between probability vector <p> and <f>
 *            in bits ($\log_2$).
 *
 */
double DRelEntropy(double *p, double *f, int n)
{
  int    i;
  double rel_entropy;
  double eps;
  double temp;

  eps = 0.0000000001;

  rel_entropy = 0.;
  for(i = 0; i < n; i++)
    {
      if (f[i] > 0.) temp = f[i]; else temp = f[i] * -1;
      if (temp < (0. + eps)) 
	{ 
	  printf("error in DRelEntropy(), f[%d] is %f\nuh not sure what to do if f[x] is 0! Abort!\n", i, f[i]); 
	  exit(1); 
	}
      if (p[i] > 0.) rel_entropy += p[i] * log(p[i] / f[i]);
    }
  return(1.44269504 * rel_entropy); /* converts to bits */
}

/***********************************************************************
 * Function: CMAverageMatchEntropy
 * Incept:   EPN, Tue May  1 14:06:37 2007
 * 
 * Purpose:  Return the average match state entropy for a CM. 
 *           Calculated as summed entropy of all match states (including
 *           MATP_MPs) divided by (2*|MATP_MP| + |MATL_ML| + |MATR_MR|).
 *
 * Args:    cm       - the covariance model
 *
 * Returns: (void) 
 */
float
CMAverageMatchEntropy(CM_t *cm)
{
  float summed_entropy, denom; 
  int v;
  summed_entropy = denom = 0.;
  for(v = 0; v < cm->M; v++)
    { 
      if(cm->stid[v] == MATP_MP)
	{
	  denom          += 2.; /* two match columns */
	  summed_entropy += esl_vec_FEntropy(cm->e[v], (MAXABET * MAXABET));
	}
      else if(cm->stid[v] == MATL_ML || 
	      cm->stid[v] == MATR_MR)
	{
	  denom          += 1.; /* two match columns */
	  summed_entropy += esl_vec_FEntropy(cm->e[v], MAXABET);
	}
    }
  /*printf("\nCM  average entropy: %.5f / %.2f = %.5f\n", summed_entropy, denom, (summed_entropy/denom));*/
  return (summed_entropy / denom);
}

/***********************************************************************
 * Function: CP9AverageMatchEntropy
 * Incept:   EPN, Tue May  1 14:13:39 2007
 * 
 * Purpose:  Return the average match state entropy for a CP9 HMM. 
 *           Calculated as summed entropy of all match states divided
 *           by number of match states. One match state per node.
 * 
 * Args:    cp9       - the CM Plan 9 HMM
 *
 * Returns: (void) 
 */
float
CP9AverageMatchEntropy(CP9_t *cp9)
{
  float summed_entropy;
  int k;
  summed_entropy = 0.;
  for(k = 1; k <= cp9->M; k++)
    summed_entropy += esl_vec_FEntropy(cp9->mat[k], MAXABET);

  /*printf("\nCP9 average entropy: %.5f / %.2f = %.5f\n", summed_entropy, (float) cp9->M, (summed_entropy/cp9->M));*/
  return (summed_entropy / ((float) cp9->M));
}


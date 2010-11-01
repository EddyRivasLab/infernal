/* p7_modelmaker.c
 * EPN, Tue Aug  5 15:32:34 2008
 * SVN $Id: cm_modelmaker.c 2327 2008-02-13 22:09:06Z nawrockie $
 *
 * Construct a p7 model from CM and it's CP9 HMM.
 *
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "easel.h"		
#include "esl_exponential.h"		
#include "esl_msa.h"		
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stack.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "hmmer.h"

#include "funcs.h"
#include "structs.h"

/* EPN, Mon Aug 25 09:00:10 2008
 * Had difficulty compiling infernal with hmmer as a subdir and with
 * impl_sse. Because I currently don't need OPROFILEs I #if 0ed out the code
 * here and made a new version of this func which doesn't return an omx (BELOW).
 */
#if 0
/* Function: BuildP7HMM_MatchEmitsOnly()
 * Incept:   EPN, Tue Aug  5 15:33:00 2008
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and it's CP9 HMM.
 *           Copy only the match emissions of the CP9 HMM, the rest of 
 *           the p7 model parameters are irrelevant. 
 * 
 * Args:     cm        - the cm
 *           ret_p7    - RETURN: new p7 model 
 *           ret_gm    - RETURN: new p7 generic profile
 *           ret_om    - RETURN: new p7 optimized profile
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7, P7_PROFILE **ret_gm, P7_OPROFILE **ret_om)
{
  int        status;
  P7_HMM     *hmm = NULL;        /* RETURN: new hmm */
  P7_PROFILE  *gm = NULL;        /* RETURN: new generic profile */
  P7_OPROFILE *om = NULL;        /* RETURN: new optimized profile */
  int        k;

  if(cm->cp9 == NULL)         return eslEINCOMPAT; 
  if(! (cm->flags & CMH_CP9)) return eslEINCOMPAT; 

  if ((hmm    = p7_hmm_Create(cm->clen, cm->abc)) == NULL)  return eslEMEM;
  if ((status = p7_hmm_Zero(hmm))                 != eslOK) return status;

  /* copy only match emissions */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->mat[k], cm->abc->K, hmm->mat[k]);

  /* parameterize, hacked from hmmer/p7_prior.c::p7_ParameterEstimation() */
  /* match transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k],   3); 

  /* insert transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k]+3, 2); 

  /* delete transitions */
  for (k = 1; k < hmm->M; k++) esl_vec_FNorm(hmm->t[k]+5, 2); 
  /* For k=0, which is unused; convention sets TMM=1.0, TMD=0.0
   * For k=M, TMM = 1.0 (to the E state) and TMD=0.0 (no next D; must go to E).
   */
  hmm->t[0][p7H_DM] = hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[0][p7H_DD] = hmm->t[hmm->M][p7H_DD] = 0.0;

  /* insert emissions */
  for (k = 0; k <= hmm->M; k++) esl_vec_FNorm(hmm->ins[k], hmm->abc->K); /* normalize inserts (0.25 each) */

  p7_hmm_SetName(hmm, cm->name);
  p7_hmm_SetAccession(hmm, cm->acc);
  p7_hmm_SetDescription(hmm, cm->desc);
  p7_hmm_SetCtime(hmm);
  if(cm->comlog != NULL && cm->comlog->bcom != NULL) { 
    ESL_ALLOC(hmm->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
    strcat(hmm->comlog, cm->comlog->bcom);
  }
  else hmm->comlog = NULL;

  hmm->eff_nseq = cm->eff_nseq;
  hmm->nseq     = cm->nseq;
  hmm->checksum = 0;

  /* make the profiles */
  gm = p7_profile_Create (hmm->M, hmm->abc);
  om = p7_oprofile_Create(hmm->M, hmm->abc);

  *ret_p7 = hmm;
  *ret_gm = gm;
  *ret_om = om;

  return eslOK;

 ERROR: 
  if(hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}
#endif

/* Function: BuildP7HMM_MatchEmitsOnly()
 * Incept:   EPN, Tue Aug  5 15:33:00 2008
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and it's CP9 HMM.
 *           Copy only the match emissions of the CP9 HMM, the rest of 
 *           the p7 model parameters are irrelevant. 
 * 
 * Args:     cm        - the cm
 *           ret_p7    - RETURN: new p7 model 
 *           ret_gm    - RETURN: new p7 generic profile
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7, P7_PROFILE **ret_gm)
{
  int        status;
  P7_HMM     *hmm = NULL;        /* RETURN: new hmm */
  P7_PROFILE  *gm = NULL;        /* RETURN: new generic profile */
  int        k;

  if(cm->cp9 == NULL)         return eslEINCOMPAT; 
  if(! (cm->flags & CMH_CP9)) return eslEINCOMPAT; 

  if ((hmm    = p7_hmm_Create(cm->clen, cm->abc)) == NULL)  return eslEMEM;
  if ((status = p7_hmm_Zero(hmm))                 != eslOK) return status;

  /* copy only match emissions */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->mat[k], cm->abc->K, hmm->mat[k]);

  /* parameterize, hacked from hmmer/p7_prior.c::p7_ParameterEstimation() */
  /* match transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k],   3); 

  /* insert transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k]+3, 2); 

  /* delete transitions */
  for (k = 1; k < hmm->M; k++) esl_vec_FNorm(hmm->t[k]+5, 2); 
  /* For k=0, which is unused; convention sets TMM=1.0, TMD=0.0
   * For k=M, TMM = 1.0 (to the E state) and TMD=0.0 (no next D; must go to E).
   */
  hmm->t[0][p7H_DM] = hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[0][p7H_DD] = hmm->t[hmm->M][p7H_DD] = 0.0;

  /* insert emissions */
  for (k = 0; k <= hmm->M; k++) esl_vec_FNorm(hmm->ins[k], hmm->abc->K); /* normalize inserts (0.25 each) */

  p7_hmm_SetName(hmm, cm->name);
  p7_hmm_SetAccession(hmm, cm->acc);
  p7_hmm_SetDescription(hmm, cm->desc);
  p7_hmm_SetCtime(hmm);
  if(cm->comlog != NULL && cm->comlog->bcom != NULL) { 
    ESL_ALLOC(hmm->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
    strcat(hmm->comlog, cm->comlog->bcom);
  }
  else hmm->comlog = NULL;

  hmm->eff_nseq = cm->eff_nseq;
  hmm->nseq     = cm->nseq;
  hmm->checksum = 0;

  /* make the profiles */
  gm = p7_profile_Create (hmm->M, hmm->abc);

  *ret_p7 = hmm;
  *ret_gm = gm;


  return eslOK;

 ERROR: 
  if(hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}

/* Function: CP9_to_P7()
 * Incept:   EPN, Fri Sep 24 13:46:37 2010
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and it's CP9 HMM.
 * 
 * Args:     cm        - the cm, must have a cp9 model in it.
 *           ret_p7    - RETURN: new p7 model 
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
CP9_to_P7(CM_t *cm, P7_HMM **ret_p7)
{
  int        status;
  P7_HMM     *hmm = NULL;        /* RETURN: new hmm */
  int        k;

  if(cm->cp9 == NULL)         return eslEINCOMPAT; 
  if(! (cm->flags & CMH_CP9)) return eslEINCOMPAT; 

  if ((hmm    = p7_hmm_Create(cm->clen, cm->abc)) == NULL)  return eslEMEM;
  if ((status = p7_hmm_Zero(hmm))                 != eslOK) return status;

  /* copy transitions */
  for (k = 0; k <= hmm->M; k++) { 
    hmm->t[k][p7H_MM] = cm->cp9->t[k][CTMM];
    hmm->t[k][p7H_MI] = cm->cp9->t[k][CTMI];
    hmm->t[k][p7H_MD] = cm->cp9->t[k][CTMD];
    hmm->t[k][p7H_IM] = cm->cp9->t[k][CTIM];
    hmm->t[k][p7H_II] = cm->cp9->t[k][CTII];
    hmm->t[k][p7H_DM] = cm->cp9->t[k][CTDM];
    hmm->t[k][p7H_DD] = cm->cp9->t[k][CTDD];
    /* note: the cp9 CTDI and CTID transitions do not exist the p7 model */
  }
  /* normalize match transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k],  3); 
  /* normalize insert transitions */
  for (k = 0; k < hmm->M; k++) esl_vec_FNorm(hmm->t[k]+3, 2); 
  /* normalize delete transitions */
  for (k = 1; k < hmm->M; k++) esl_vec_FNorm(hmm->t[k]+5, 2); 

  /* enforce HMMER conventions */
  hmm->t[hmm->M][p7H_MD] = 0.0;
  esl_vec_FNorm(hmm->t[hmm->M], 3);
  hmm->t[0][p7H_DM] = hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[0][p7H_DD] = hmm->t[hmm->M][p7H_DD] = 0.0;

  /* enforce INFERNAL CP9 convention, the 0'th node's MM transition is really begin[0] */
  hmm->t[0][p7H_MM] = cm->cp9->begin[1];
  esl_vec_FNorm(hmm->t[0], 3);

  /* match emissions: copy, then normalize (should be unnec actually) */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->mat[k], cm->abc->K, hmm->mat[k]);
  for (k = 1; k <= cm->clen; k++) esl_vec_FNorm(cm->cp9->mat[k], cm->abc->K);
  /* special case */
  esl_vec_FSet(hmm->mat[0], hmm->abc->K, 0.);
  hmm->mat[0][0] = 1.0;

  /* insert emissions: copy, then normalize (should be unnec actually) */
  for (k = 0; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->ins[k], cm->abc->K, hmm->ins[k]);
  for (k = 0; k <= cm->clen; k++) esl_vec_FNorm(cm->cp9->ins[k], cm->abc->K);

  /* copy the all important max length parameter */
  hmm->max_length = cm->W;

  p7_hmm_SetName       (hmm, cm->name);
  p7_hmm_SetAccession  (hmm, cm->acc);
  p7_hmm_SetDescription(hmm, cm->desc);
  p7_hmm_SetCtime      (hmm);
  if(cm->comlog != NULL && cm->comlog->bcom != NULL) { 
    ESL_ALLOC(hmm->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
    strcat(hmm->comlog, cm->comlog->bcom);
  }
  else hmm->comlog = NULL;

  hmm->eff_nseq = cm->eff_nseq;
  hmm->nseq     = cm->nseq;
  hmm->checksum = 0;

  /* set the model composition */
  if ((status = p7_hmm_SetComposition(hmm)) != eslOK) goto ERROR;

  /* finally, calibrate the model */
  P7_OPROFILE    *om = NULL;
  P7_BG          *bg = NULL;
  P7_PROFILE     *gm = NULL;
  ESL_RANDOMNESS *r  = NULL;

  if((status = p7_Calibrate(hmm, NULL, &r, &bg, &gm, &om)) != eslOK) goto ERROR; 

  /* experimental step: try to determine a forward tau for glocal mode */
  double mu;
  double lambda;
  /*if ((r = esl_randomness_CreateFast(42)) == NULL) { status = eslEMEM; goto ERROR; }
    if ((bg     = p7_bg_Create(hmm->abc)) == NULL)   { status = eslEMEM; goto ERROR; }
    if ((om     = p7_oprofile_Create(hmm->M, hmm->abc)) == NULL) { status = eslEMEM; goto ERROR; }
    if ((status = p7_oprofile_Convert(gm, om))         != eslOK) goto ERROR; */

  if ((status = p7_ProfileConfig(hmm, bg, gm, cm->W*2, p7_GLOCAL)) != eslOK) goto ERROR; 
  if ((status = p7_GlocalLambdaMu(r, gm, bg, cm->W*2, 10000, 0.002, &lambda, &mu)) != eslOK) goto ERROR; 

  printf("\n\n\n p7 glocal lambda: %g  mu: %g\n", lambda, mu);
  printf("cp9  local lambda: %g  mu: %g\n", cm->stats->expAA[EXP_CP9_GF][0]->lambda, cm->stats->expAA[EXP_CP9_GF][0]->mu_extrap);

  double sc = -10.0; 
  for(sc = -10.0; sc < 40.0; sc += 1.0) { 
    ;/*printf("%5.2f  %15g  %15g  %15g\n", sc, 
	   esl_exp_surv(sc,  hmm->evparam[p7_FTAU],  hmm->evparam[p7_FLAMBDA]),
	   esl_exp_surv(sc,  mu, lambda),
	   esl_exp_surv(sc,  cm->stats->expAA[EXP_CP9_GF][0]->mu_extrap, cm->stats->expAA[EXP_CP9_GF][0]->lambda));*/
  }
  cm->p7_glambda = lambda;
  cm->p7_gmu     = mu;

  esl_randomness_Destroy(r); 
  p7_bg_Destroy(bg);         
  p7_oprofile_Destroy(om);   
  p7_profile_Destroy(gm);   
  
  *ret_p7 = hmm;

  return eslOK;

 ERROR: 
  if(hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}


/* Function: dump_p7()
 * Incept:   EPN, Fri Sep 24 14:22:49 2010
 * 
 * Purpose:  Dump parameters of a p7 HMM to a file.
 * 
 * Args:     hmm       - the p7 HMM
 *           fp        - the file to print to
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
dump_p7(P7_HMM *hmm, FILE *fp)
{
  return p7_hmmfile_WriteASCII(fp, -1, hmm);
}


/* Function:  p7_GlocalLambdaMu()
 * Synopsis:  Determine Forward lambda and tau for glocal mode by simulation.
 * Incept:    EPN, Wed Oct 27 06:49:17 2010
 *            SRE, Thu Aug  9 15:08:39 2007 [p7_Tau()]
 *
 * Purpose:   This function is a modified version of p7_Tau(), with
 *            changes only made as necessary to use scores from a
 *            glocally configured HMM, instead of a locally configured
 *            on. One important difference is that lambda is estimated
 *            not passed in, another is that a generic profile <gm>
 *            is used here instead of the optimized one used by
 *            p7_Tau().  This is necessary b/c glocal scores cannot be
 *            calculated using the optimized routines because they
 *            make assumptions that are violated by a glocal
 *            configuration.  All notes below are from p7_Tau():
 *
 *            Determine the <tau> parameter for an exponential tail
 *            fit to the Forward score distribution for model <om>, on
 *            random sequences with the composition of the background
 *            model <bg>. This <tau> parameter is for an exponential
 *            distribution anchored from $P=1.0$, so it's not really a
 *            tail per se; but it's only an accurate fit in the tail
 *            of the Forward score distribution, from about $P=0.001$
 *            or so.
 *            
 *            The determination of <tau> is done by a brief simulation
 *            in which we fit a Gumbel distribution to a small number
 *            of Forward scores of random sequences, and use that to
 *            predict the location of the tail at probability <tailp>.
 *            
 *            The Gumbel is of course inaccurate, but we can use it
 *            here solely as an empirical distribution to determine
 *            the location of a reasonable <tau> more accurately on a
 *            smaller number of samples than we could do with raw
 *            order statistics. 
 *            
 *            Typical choices are L=100, N=200, tailp=0.04, which
 *            typically yield estimates $\hat{\mu}$ with a precision
 *            (standard deviation) of $\pm$ 0.2 bits, corresponding to
 *            a $\pm$ 15\% error in E-values. See [J1/135].
 *            
 *            The use of Gumbel fitting to a small number of $N$
 *            samples and the extrapolation of $\hat{\mu}$ from the
 *            estimated location of the 0.04 tail mass are both
 *            empirical and carefully optimized against several
 *            tradeoffs. Most importantly, around this choice of tail
 *            probability, a systematic error introduced by the use of
 *            the Gumbel fit is being cancelled by systematic error
 *            introduced by the use of a higher tail probability than
 *            the regime in which the exponential tail is a valid
 *            approximation. See [J1/135] for discussion.
 *            
 *            This function changes the length configuration of both
 *            <om> and <bg>. The caller must remember to reconfigure
 *            both of their length models appropriately for any
 *            subsequent alignments.
 *            
 * Args:      r      : source of randomness
 *            gm     : configured profile to score sequences with
 *            bg     : null model (for background residue frequencies)
 *            L      : mean length model for seq emission from profile
 *            N      : number of sequences to generate
 *            tailp  : tail mass from which we will extrapolate tau
 *            ret_lambda: RETURN: estimate for the Forward lambda
 *            ret_mu:     RETURN: estimate for the Forward mu
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_tau> is 0.
 */
int
p7_GlocalLambdaMu(ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_BG *bg, int L, int N, double tailp, double *ret_lambda, double *ret_mu)
{
  P7_GMX  *gx      = p7_gmx_Create(gm->M, L); /* DP matrix: for ForwardParser,  L rows */
  ESL_DSQ *dsq     = NULL;
  double  *xv      = NULL;
  float    fsc, nullsc;		                  
  double   gmu, glam;
  int      status;
  int      i;
  int      n;
  ESL_HISTOGRAM *h = NULL;

  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if ((h = esl_histogram_CreateFull(-50., 50., 0.2)) == NULL) { status = eslEMEM; goto ERROR; }
  if (gx == NULL) { status = eslEMEM; goto ERROR; }

  p7_ReconfigLength(gm, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rsq_xfIID(r, bg->f, gm->abc->K, L, dsq)) != eslOK) goto ERROR;
      if ((status = p7_GForward(dsq, L, gm, gx, &fsc))           != eslOK) goto ERROR;
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;   
      esl_histogram_Add(h, ((fsc - nullsc) / eslCONST_LOG2));
    }

  /*esl_histogram_Print(stdout, h);*/

  esl_histogram_GetTailByMass(h, tailp, &xv, &n, NULL);
  if ((status = esl_exp_FitComplete(xv, n, &gmu, &glam)) != eslOK) goto ERROR;

  /* Explanation of the eqn below: first find the x at which the Gumbel tail
   * mass is predicted to be equal to tailp. Then back up from that x
   * by log(tailp)/lambda to set the origin of the exponential tail to 1.0
   * instead of tailp.
   */
  *ret_mu = gmu - log(1./tailp) / glam;

  *ret_lambda =  glam;
  
  free(dsq);
  p7_gmx_Destroy(gx);
  esl_histogram_Destroy(h);
  return eslOK;

 ERROR:
  *ret_mu     = 0.;
  *ret_lambda = 0.;
  if (xv  != NULL) free(xv);
  if (dsq != NULL) free(dsq);
  if (gx  != NULL) p7_gmx_Destroy(gx);
  if (h   != NULL) esl_histogram_Destroy(h);
  return status;
}


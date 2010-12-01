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

/* Function: cm_cp9_to_p7()
 * Incept:   EPN, Fri Sep 24 13:46:37 2010
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and it's CP9 HMM.
 * 
 * Args:     cm        - the cm, must have a cp9 model in it.
 *
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
cm_cp9_to_p7(CM_t *cm)
{
  int        status;
  int        k;

  if(cm->cp9 == NULL)         return eslEINCOMPAT; 
  if(! (cm->flags & CMH_CP9)) return eslEINCOMPAT; 
  if(cm->mlp7 != NULL)          return eslEINCOMPAT;

  if ((cm->mlp7 = p7_hmm_Create(cm->clen, cm->abc)) == NULL)  return eslEMEM;
  if ((status = p7_hmm_Zero(cm->mlp7))                 != eslOK) return status;

  /* copy transitions */
  for (k = 0; k <= cm->mlp7->M; k++) { 
    cm->mlp7->t[k][p7H_MM] = cm->cp9->t[k][CTMM];
    cm->mlp7->t[k][p7H_MI] = cm->cp9->t[k][CTMI];
    cm->mlp7->t[k][p7H_MD] = cm->cp9->t[k][CTMD];
    cm->mlp7->t[k][p7H_IM] = cm->cp9->t[k][CTIM];
    cm->mlp7->t[k][p7H_II] = cm->cp9->t[k][CTII];
    cm->mlp7->t[k][p7H_DM] = cm->cp9->t[k][CTDM];
    cm->mlp7->t[k][p7H_DD] = cm->cp9->t[k][CTDD];
    /* note: the cp9 CTDI and CTID transitions do not exist the p7 model */
  }
  /* normalize match transitions */
  for (k = 1; k <= cm->mlp7->M; k++) esl_vec_FNorm(cm->mlp7->t[k],  3); 
  /* normalize insert transitions */
  for (k = 0; k < cm->mlp7->M; k++) esl_vec_FNorm(cm->mlp7->t[k]+3, 2); 
  /* normalize delete transitions */
  for (k = 1; k < cm->mlp7->M; k++) esl_vec_FNorm(cm->mlp7->t[k]+5, 2); 

  /* enforce HMMER conventions */
  cm->mlp7->t[cm->mlp7->M][p7H_MD] = 0.0;
  esl_vec_FNorm(cm->mlp7->t[cm->mlp7->M], 3);
  cm->mlp7->t[0][p7H_DM] = cm->mlp7->t[cm->mlp7->M][p7H_DM] = 1.0;
  cm->mlp7->t[0][p7H_DD] = cm->mlp7->t[cm->mlp7->M][p7H_DD] = 0.0;

  /* enforce INFERNAL CP9 convention, the 0'th node's MM transition is really begin[0] */
  cm->mlp7->t[0][p7H_MM] = cm->cp9->begin[1];
  esl_vec_FNorm(cm->mlp7->t[0], 3);

  /* match emissions: copy, then normalize (should be unnec actually) */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->mat[k], cm->abc->K, cm->mlp7->mat[k]);
  for (k = 1; k <= cm->clen; k++) esl_vec_FNorm(cm->cp9->mat[k], cm->abc->K);
  /* special case */
  esl_vec_FSet(cm->mlp7->mat[0], cm->mlp7->abc->K, 0.);
  cm->mlp7->mat[0][0] = 1.0;

  /* insert emissions: copy, then normalize (should be unnec actually) */
  for (k = 0; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->ins[k], cm->abc->K, cm->mlp7->ins[k]);
  for (k = 0; k <= cm->clen; k++) esl_vec_FNorm(cm->cp9->ins[k], cm->abc->K);

  /* copy the max length parameter */
  cm->mlp7->max_length = cm->W;

  p7_hmm_SetName       (cm->mlp7, cm->name);
  p7_hmm_SetAccession  (cm->mlp7, cm->acc);
  p7_hmm_SetDescription(cm->mlp7, cm->desc);
  p7_hmm_SetCtime      (cm->mlp7);
  if(cm->comlog != NULL && cm->comlog->bcom != NULL) { 
    ESL_ALLOC(cm->mlp7->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
    *(cm->mlp7->comlog) = '\0'; /* need this to make strcat work */
    strcat(cm->mlp7->comlog, cm->comlog->bcom);
  }
  else cm->mlp7->comlog = NULL;

  cm->mlp7->eff_nseq = cm->eff_nseq;
  cm->mlp7->nseq     = cm->nseq;
  cm->mlp7->checksum = 0;

  /* set the model composition */
  if ((status = p7_hmm_SetComposition(cm->mlp7)) != eslOK) goto ERROR;

  cm->flags |= CMH_MLP7; /* raise the P7 flag */

  return eslOK;

 ERROR:
  if(cm->mlp7 != NULL) { p7_hmm_Destroy(cm->mlp7); cm->mlp7 = NULL; }
  return status;
}

/* Function: cm_p7_Calibrate()
 * Incept:   EPN, Tue Nov  9 06:16:57 2010
 * 
 * Purpose:  Calibrate a CM's p7 HMM for local MSV, Viterbi, Forward and 
 *           also glocal Forward. 
 * 
 * Args:     cm        - the cm, cm->mlp7 must be non-NULL
 *           errbuf    - for error messages
 *           ElmL      - length of sequences to sample for local MSV
 *           ElvL      - length of sequences to sample for local Vit
 *           ElfL      - length of sequences to sample for local Fwd
 *           EgfL      - length of sequences to sample for glocal Fwd
 *           ElmN      - number of sequences to sample for local MSV
 *           ElvN      - number of sequences to sample for local Vit
 *           ElfN      - number of sequences to sample for local Fwd
 *           EgfN      - number of sequences to sample for glocal Fwd
 *           ElfT      - fraction of tail mass to fit for  local Fwd (usually (HMMER3 is) 0.04)
 *           EgfT      - fraction of tail mass to fit for glocal Fwd 
 *           do_real   - TRUE to sample realistic genomic seqs, not IID
 *           do_null3  - TRUE to use null3 correction on scores before tail fit
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
cm_p7_Calibrate(CM_t *cm, char *errbuf, 
		int ElmL, int ElvL, int ElfL, int EgfL, 
		int ElmN, int ElvN, int ElfN, int EgfN, 
		double ElfT, double EgfT, int do_real, int do_null3, float n3omega)
{
  int        status;
  P7_OPROFILE    *om = NULL;
  P7_BG          *bg = NULL;
  P7_PROFILE     *gm = NULL;
  ESL_RANDOMNESS *r  = NULL;
  double lmmu, lvmu, lftau, gfmu;
  double lambda, gflambda;

  printf("cm_p7_Calibrate:\n\tElmL: %d\n\tElvL: %d\n\tElfL: %d\n\tEgfL: %d\n\tElmN: %d\n\tElvN: %d\n\tElfN: %d\n\tEgfN: %d\n\tElfT: %f\n\tEgfT: %f\n\tdo_real: %d\n\tdo_null3: %d\n", ElmL, ElvL, ElfL, EgfL, ElmN, ElvN, ElfN, EgfN, ElfT, EgfT, do_real, do_null3);

  if(cm->mlp7 == NULL)       { status = eslEINCOMPAT; goto ERROR; }
  if(! cm->flags & CMH_MLP7) { status = eslEINCOMPAT; goto ERROR; }

  /* most of this code stolen from hmmer's evalues.c::p7_Calibrate() */
  if ((r      = esl_randomness_CreateFast(42)) == NULL)                    ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to create RNG");
  if ((bg     = p7_bg_Create(cm->mlp7->abc)) == NULL)                        ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to allocate background");
  if ((gm     = p7_profile_Create(cm->mlp7->M, cm->mlp7->abc))  == NULL)       ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to allocate profile");
  if ((status = p7_ProfileConfig(cm->mlp7, bg, gm, ElmL, p7_LOCAL)) != eslOK) ESL_XFAIL(status,  errbuf, "cm_p7_Calibrate(): failed to configure profile");
  if ((om     = p7_oprofile_Create(cm->mlp7->M, cm->mlp7->abc)) == NULL)       ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to create optimized profile");
  if ((status = p7_oprofile_Convert(gm, om)) != eslOK)                     ESL_XFAIL(status,  errbuf, "cm_p7_Calibrate(): failed to convert to optimized profile");


  /* The calibration steps themselves */
  if ((status = p7_Lambda      (cm->mlp7, bg, &lambda))                         != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine lambda");
  if ((status = cm_p7_MSVMu    (r, errbuf, om, bg, ElmL, ElmN, lambda, do_real, do_null3, n3omega, &lmmu))        != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine msv mu");
  if (ElvL != ElmL) p7_oprofile_ReconfigLength(om, ElvL);
  if ((status = cm_p7_ViterbiMu(r, errbuf, om, bg, ElvL, ElvN, lambda, do_real, do_null3, n3omega, &lvmu))        != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine vit mu");
  if (ElfL != ElvL) p7_oprofile_ReconfigLength(om, ElfL);
  if ((status = cm_p7_Tau      (r, errbuf, om, bg, ElfL, ElfN, lambda, ElfT, do_real, do_null3, n3omega, &lftau)) != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine fwd tau");

  /*if ((status = p7_Calibrate(cm->mlp7, NULL, &r, &bg, &gm, &om)) != eslOK) goto ERROR;*/

  /* finally, determine Glocal Forward stats */
  if ((status = p7_ProfileConfig(cm->mlp7, bg, gm, EgfL, p7_GLOCAL)) != eslOK) goto ERROR; 
  if ((status = p7_GlocalLambdaMu(cm, r, gm, bg, do_real, do_null3, n3omega, EgfL, EgfN, EgfT, errbuf, &gflambda, &gfmu)) != eslOK) goto ERROR; 

  /* copy the p7's evparam[], which were set in p7_Calibrate() to the cm's p7_evparam,
   * which will additionally store the Glocal Lambda and Mu for Forward */
  cm->mlp7_evparam[CM_p7_LMMU]     = cm->mlp7->evparam[p7_MMU]     = lmmu;
  cm->mlp7_evparam[CM_p7_LMLAMBDA] = cm->mlp7->evparam[p7_MLAMBDA] = lambda;
  cm->mlp7_evparam[CM_p7_LVMU]     = cm->mlp7->evparam[p7_VMU]     = lvmu;  
  cm->mlp7_evparam[CM_p7_LVLAMBDA] = cm->mlp7->evparam[p7_VLAMBDA] = lambda;
  cm->mlp7_evparam[CM_p7_LFTAU]    = cm->mlp7->evparam[p7_FTAU]    = lftau; 
  cm->mlp7_evparam[CM_p7_LFLAMBDA] = cm->mlp7->evparam[p7_FLAMBDA] = lambda;  

  cm->mlp7_evparam[CM_p7_GFMU]     = gfmu;
  cm->mlp7_evparam[CM_p7_GFLAMBDA] = gflambda;

  cm->flags |= CMH_MLP7_STATS;

  printf("p7 glocal lambda: %g  mu: %g\n", gflambda, gfmu);

  esl_randomness_Destroy(r); 
  p7_bg_Destroy(bg);         
  p7_oprofile_Destroy(om);   
  p7_profile_Destroy(gm);   
  
  return eslOK;

 ERROR: 
  if(r != NULL)  esl_randomness_Destroy(r); 
  if(bg != NULL) p7_bg_Destroy(bg);         
  if(om != NULL) p7_oprofile_Destroy(om);   
  if(gm != NULL) p7_profile_Destroy(gm);   
  return status;
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
 * Args:      cm     : the model
 *            r      : source of randomness
 *            gm     : configured profile to score sequences with
 *            bg     : null model (for background residue frequencies)
 *            do_real: sample realistic genomic sequences, don't use iid
 *            do_null3: TRUE to use null3 correction on scores, FALSE not to
 *            L      : mean length model for seq emission from profile
 *            N      : number of sequences to generate
 *            tailp  : tail mass from which we will extrapolate tau
 *            errbuf : for error messages
 *            ret_lambda: RETURN: estimate for the Forward lambda
 *            ret_mu:     RETURN: estimate for the Forward mu
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_tau> is 0.
 */
int
p7_GlocalLambdaMu(CM_t *cm, ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_BG *bg, int do_real, int do_null3, float n3omega, int L, int N, double tailp, char *errbuf, double *ret_lambda, double *ret_mu)
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
  float    null3sc = 0.;
  float    sc;

  /* the HMM that generates sequences */
  int     ghmm_nstates = 0;       /* number of states in the HMM */
  double  *ghmm_sA  = NULL;       /* start probabilities [0..ghmm_nstates-1] */
  double **ghmm_tAA = NULL;       /* transition probabilities [0..nstates-1][0..nstates-1] */
  double **ghmm_eAA = NULL;       /* emission probabilities   [0..nstates-1][0..abc->K-1] */

  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if ((h = esl_histogram_CreateFull(-50., 50., 0.2)) == NULL) { status = eslEMEM; goto ERROR; }
  if (gx == NULL) { status = eslEMEM; goto ERROR; }

  if(do_real) { 
    if((status = CreateGenomicHMM(cm->abc, errbuf, &ghmm_sA, &ghmm_tAA, &ghmm_eAA, &ghmm_nstates)) != eslOK) goto ERROR;
  }

  p7_ReconfigLength(gm, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if(do_real) { if((status = SampleGenomicSequenceFromHMM(r, cm->abc, errbuf, ghmm_sA, ghmm_tAA, ghmm_eAA, ghmm_nstates, L, &dsq) != eslOK)) goto ERROR; }
      else        { if((status = esl_rsq_xfIID(r, bg->f, gm->abc->K, L, dsq)) != eslOK) goto ERROR; }
      if ((status = p7_GForward(dsq, L, gm, gx, &fsc))           != eslOK) goto ERROR;
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;   
      sc = ((fsc-nullsc) / eslCONST_LOG2);
      if(do_null3) { 
	ScoreCorrectionNull3CompUnknown(cm->abc, cm->null, dsq, 1, L, n3omega, &null3sc);
	null3sc *= (float) cm->clen / (float) L; /* assume hit would be of length clen, not full window len */
	sc -= null3sc;
      }
      esl_histogram_Add(h, sc);
      if(do_real) { free(dsq); dsq = NULL; }
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

  /* free HMM if nec */
  if(do_real) { 
    for(i = 0; i < ghmm_nstates; i++) { 
      free(ghmm_eAA[i]); 
      free(ghmm_tAA[i]); 
    }
    free(ghmm_eAA);
    free(ghmm_tAA);
    free(ghmm_sA);
  }
  
  if(dsq != NULL) free(dsq);
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
  for(i = 0; i < ghmm_nstates; i++) { 
    if(ghmm_eAA != NULL) free(ghmm_eAA[i]); 
    if(ghmm_tAA != NULL) free(ghmm_tAA[i]); 
  }
  if(ghmm_eAA != NULL) free(ghmm_eAA);
  if(ghmm_tAA != NULL) free(ghmm_tAA);
  if(ghmm_sA  != NULL) free(ghmm_sA);
  return status;
}


/* Function:  cm_p7_MSVMu()
 * Synopsis:  Determines the local MSV Gumbel mu parameter for a model.
 * Incept:    SRE, Mon Aug  6 13:00:57 2007 [Janelia] (p7_MSVMu())
 *
 * Purpose:   Identical to p7_MSVMu() except that options <do_real> and
 *            <do_null3> allow target sequences to be generated by 
 *            a 5-state HMM that generates genome-like background sequence
 *            (if <do_real> is TRUE) and applies a null3 penalty
 *            (if <do_null3> is TRUE). See hmmer/evalues.c::cm_p7_MSVMu
 *            for additional information.
 *            
 * Args:      r       :  source of random numbers
 *            om      :  score profile (length config is changed upon return!)
 *            bg      :  null model    (length config is changed upon return!)
 *            L       :  length of sequences to simulate
 *            N	      :  number of sequences to simulate		
 *            lambda  :  known Gumbel lambda parameter
 *            do_real :  TRUE to generate target seqs from genomic HMM, FALSE to do iid
 *            do_null3:  TRUE to apply null3 penalty with <n3omega> omega.
 *            n3omega :  omega for null3, irrelevant if do_null3.
 *            ret_mmu :  RETURN: ML estimate of location param mu
 *
 * Returns:   <eslOK> on success, and <ret_mu> contains the ML estimate
 *            of $\mu$.
 *
 * Throws:    (no abnormal error conditions)
 * 
 * Note:      The FitCompleteLoc() function is simple, and it's tempting
 *            to inline it here and save the <xv> working memory. However,
 *            the FitCompleteLoc() function is vulnerable
 *            to under/overflow error, and we'll probably fix it
 *            eventually - need to be sure that fix applies here too.
 */
int
cm_p7_MSVMu(ESL_RANDOMNESS *r, char *errbuf, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda, int do_real, int do_null3, float n3omega, double *ret_mmu)
{
  P7_OMX  *ox      = p7_omx_Create(om->M, 0, 0); /* DP matrix: 1 row version */
  ESL_DSQ *dsq     = NULL;
  double  *xv      = NULL;
  int      i;
  float    sc, nullsc, null3sc;
#ifndef p7_IMPL_DUMMY
  float    maxsc   = (255 - om->base_b) / om->scale_b; /* if score overflows, use this */
#endif
  int      status;

  /* the HMM that generates sequences if do_real==TRUE */
  int     ghmm_nstates = 0;       /* number of states in the HMM */
  double  *ghmm_sA  = NULL;       /* start probabilities [0..ghmm_nstates-1] */
  double **ghmm_tAA = NULL;       /* transition probabilities [0..nstates-1][0..nstates-1] */
  double **ghmm_eAA = NULL;       /* emission probabilities   [0..nstates-1][0..abc->K-1] */

  if (ox == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  if(do_real) { 
    if((status = CreateGenomicHMM(bg->abc, errbuf, &ghmm_sA, &ghmm_tAA, &ghmm_eAA, &ghmm_nstates)) != eslOK) goto ERROR;
  }

  p7_oprofile_ReconfigLength(om, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if(do_real) { if((status = SampleGenomicSequenceFromHMM(r, bg->abc, errbuf, ghmm_sA, ghmm_tAA, ghmm_eAA, ghmm_nstates, L, &dsq) != eslOK)) goto ERROR; }
      else        { if((status = esl_rsq_xfIID(r, bg->f, bg->abc->K, L, dsq)) != eslOK) goto ERROR; }
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;   

      status = p7_MSVFilter(dsq, L, om, ox, &sc); 
#ifndef p7_IMPL_DUMMY
      if (status == eslERANGE) { sc = maxsc; status = eslOK; }
#endif
      if (status != eslOK)     goto ERROR;
      sc = (sc - nullsc) / eslCONST_LOG2;
      if(do_null3) { 
	ScoreCorrectionNull3CompUnknown(bg->abc, bg->f, dsq, 1, L, n3omega, &null3sc);
	null3sc *= (float) om->M / (float) L; /* assume hit would be of length clen, not full window len */
	sc -= null3sc;
      }
      xv[i] = sc;
    }

  if ((status = esl_gumbel_FitCompleteLoc(xv, N, lambda, ret_mmu))  != eslOK) goto ERROR;
  p7_omx_Destroy(ox);
  free(xv);
  free(dsq);
  /* free HMM if nec */
  if(do_real) { 
    for(i = 0; i < ghmm_nstates; i++) { 
      free(ghmm_eAA[i]); 
      free(ghmm_tAA[i]); 
    }
    free(ghmm_eAA);
    free(ghmm_tAA);
    free(ghmm_sA);
  }
  
  return eslOK;

 ERROR:
  *ret_mmu = 0.0;
  if (ox  != NULL) p7_omx_Destroy(ox);
  if (xv  != NULL) free(xv);
  if (dsq != NULL) free(dsq);
  return status;
}

/* Function:  cm_p7_ViterbiMu()
 * Synopsis:  Determines the local Viterbi Gumbel mu parameter for a model.
 * Incept:    SRE, Tue May 19 10:26:19 2009 [Janelia] (p7_ViterbiMu())
 *
 * Purpose:   Identical to cm_p7_ViterbiMu() except that options
 *            <do_real> and <do_null3> allow target sequences to be
 *            generated by a 5-state HMM that generates genome-like
 *            background sequence (if <do_real> is TRUE) and applies a
 *            null3 penalty (if <do_null3> is TRUE). See
 *            hmmer/evalues.c::cm_p7_ViterbiMu for additional information.
 *
 * Args:      r       :  source of random numbers
 *            om      :  score profile (length config is changed upon return!)
 *            bg      :  null model    (length config is changed upon return!)
 *            L       :  length of sequences to simulate
 *            N	      :  number of sequences to simulate		
 *            lambda  :  known Gumbel lambda parameter
 *            do_real :  TRUE to generate target seqs from genomic HMM, FALSE to do iid
 *            do_null3:  TRUE to apply null3 penalty with <n3omega> omega.
 *            n3omega :  omega for null3, irrelevant if do_null3.
 *            ret_vmu :  RETURN: ML estimate of location param mu
 *
 * Returns:   <eslOK> on success, and <ret_mu> contains the ML estimate
 *            of $\mu$.
 *
 * Throws:    (no abnormal error conditions)
 */
int
cm_p7_ViterbiMu(ESL_RANDOMNESS *r, char *errbuf,P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda, int do_real, int do_null3, float n3omega, double *ret_vmu)
{
  P7_OMX  *ox      = p7_omx_Create(om->M, 0, 0); /* DP matrix: 1 row version */
  ESL_DSQ *dsq     = NULL;
  double  *xv      = NULL;
  int      i;
  float    sc, nullsc, null3sc;
#ifndef p7_IMPL_DUMMY
  float    maxsc   = (32767.0 - om->base_w) / om->scale_w; /* if score overflows, use this [J4/139] */
#endif
  int      status;

  /* the HMM that generates sequences if do_real==TRUE */
  int     ghmm_nstates = 0;       /* number of states in the HMM */
  double  *ghmm_sA  = NULL;       /* start probabilities [0..ghmm_nstates-1] */
  double **ghmm_tAA = NULL;       /* transition probabilities [0..nstates-1][0..nstates-1] */
  double **ghmm_eAA = NULL;       /* emission probabilities   [0..nstates-1][0..abc->K-1] */

  if (ox == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  if(do_real) { 
    if((status = CreateGenomicHMM(bg->abc, errbuf, &ghmm_sA, &ghmm_tAA, &ghmm_eAA, &ghmm_nstates)) != eslOK) goto ERROR;
  }

  p7_oprofile_ReconfigLength(om, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if(do_real) { if((status = SampleGenomicSequenceFromHMM(r, bg->abc, errbuf, ghmm_sA, ghmm_tAA, ghmm_eAA, ghmm_nstates, L, &dsq) != eslOK)) goto ERROR; }
      else        { if((status = esl_rsq_xfIID(r, bg->f, bg->abc->K, L, dsq)) != eslOK) goto ERROR; }
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;   

      status = p7_ViterbiFilter(dsq, L, om, ox, &sc); 
#ifndef p7_IMPL_DUMMY
      if (status == eslERANGE) { sc = maxsc; status = eslOK; }
#endif
      if (status != eslOK)     goto ERROR;
      sc = (sc - nullsc) / eslCONST_LOG2;
      if(do_null3) { 
	ScoreCorrectionNull3CompUnknown(bg->abc, bg->f, dsq, 1, L, n3omega, &null3sc);
	null3sc *= (float) om->M / (float) L; /* assume hit would be of length clen, not full window len */
	sc -= null3sc;
      }
      xv[i] = sc;
    }

  if ((status = esl_gumbel_FitCompleteLoc(xv, N, lambda, ret_vmu))  != eslOK) goto ERROR;
  p7_omx_Destroy(ox);
  free(xv);
  free(dsq);
  /* free HMM if nec */
  if(do_real) { 
    for(i = 0; i < ghmm_nstates; i++) { 
      free(ghmm_eAA[i]); 
      free(ghmm_tAA[i]); 
    }
    free(ghmm_eAA);
    free(ghmm_tAA);
    free(ghmm_sA);
  }
  return eslOK;

 ERROR:
  *ret_vmu = 0.0;
  if (ox  != NULL) p7_omx_Destroy(ox);
  if (xv  != NULL) free(xv);
  if (dsq != NULL) free(dsq);
  return status;

}


/* Function:  cm_p7_Tau()
 * Synopsis:  Determine Forward tau by brief simulation.
 * Incept:    SRE, Thu Aug  9 15:08:39 2007 [Janelia] (p7_Tau())
 *
 * Purpose:   Identical to p7_Tau() except that options <do_real> and
 *            <do_null3> allow target sequences to be generated by 
 *            a 5-state HMM that generates genome-like background sequence
 *            (if <do_real> is TRUE) and applies a null3 penalty
 *            (if <do_null3> is TRUE). See hmmer/evalues.c::cm_p7_Tau
 *            for additional information.
 *            
 * Args:      r      : source of randomness
 *            om     : configured profile to sample sequences from
 *            bg     : null model (for background residue frequencies)
 *            L      : mean length model for seq emission from profile
 *            N      : number of sequences to generate
 *            lambda : expected slope of the exponential tail (from p7_Lambda())
 *            tailp  : tail mass from which we will extrapolate mu
 *            do_real :  TRUE to generate target seqs from genomic HMM, FALSE to do iid
 *            do_null3:  TRUE to apply null3 penalty with <n3omega> omega.
 *            n3omega :  omega for null3, irrelevant if do_null3.
 *            ret_mu : RETURN: estimate for the Forward mu (base of exponential tail)
 *
 * Returns:   <eslOK> on success, and <*ret_fv> is the score difference
 *            in bits.
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_fv> is 0.
 */
int
cm_p7_Tau(ESL_RANDOMNESS *r, char *errbuf, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda, double tailp, int do_real, int do_null3, float n3omega, double *ret_tau)
{
  P7_OMX  *ox      = p7_omx_Create(om->M, 0, L);     /* DP matrix: for ForwardParser,  L rows */
  ESL_DSQ *dsq     = NULL;
  double  *xv      = NULL;
  float    sc, fsc, nullsc, null3sc;		                  
  double   gmu, glam;
  int      status;
  int      i;

  /* the HMM that generates sequences if do_real==TRUE */
  int     ghmm_nstates = 0;       /* number of states in the HMM */
  double  *ghmm_sA  = NULL;       /* start probabilities [0..ghmm_nstates-1] */
  double **ghmm_tAA = NULL;       /* transition probabilities [0..nstates-1][0..nstates-1] */
  double **ghmm_eAA = NULL;       /* emission probabilities   [0..nstates-1][0..abc->K-1] */

  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if (ox == NULL) { status = eslEMEM; goto ERROR; }

  if(do_real) { 
    if((status = CreateGenomicHMM(bg->abc, errbuf, &ghmm_sA, &ghmm_tAA, &ghmm_eAA, &ghmm_nstates)) != eslOK) goto ERROR;
  }

  p7_oprofile_ReconfigLength(om, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if(do_real) { if((status = SampleGenomicSequenceFromHMM(r, bg->abc, errbuf, ghmm_sA, ghmm_tAA, ghmm_eAA, ghmm_nstates, L, &dsq) != eslOK)) goto ERROR; }
      else        { if((status = esl_rsq_xfIID(r, bg->f, bg->abc->K, L, dsq)) != eslOK) goto ERROR; }
      if ((status = p7_ForwardParser(dsq, L, om, ox, &fsc))      != eslOK) goto ERROR;
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;   
      sc = (fsc - nullsc) / eslCONST_LOG2;
      if(do_null3) { 
	ScoreCorrectionNull3CompUnknown(bg->abc, bg->f, dsq, 1, L, n3omega, &null3sc);
	null3sc *= (float) om->M / (float) L; /* assume hit would be of length clen, not full window len */
	sc -= null3sc;
      }
      xv[i] = sc;
    }
  if ((status = esl_gumbel_FitComplete(xv, N, &gmu, &glam)) != eslOK) goto ERROR;

  /* Explanation of the eqn below: first find the x at which the Gumbel tail
   * mass is predicted to be equal to tailp. Then back up from that x
   * by log(tailp)/lambda to set the origin of the exponential tail to 1.0
   * instead of tailp.
   */
  *ret_tau =  esl_gumbel_invcdf(1.0-tailp, gmu, glam) + (log(tailp) / lambda);
  
  free(xv);
  free(dsq);
  p7_omx_Destroy(ox);
  /* free HMM if nec */
  if(do_real) { 
    for(i = 0; i < ghmm_nstates; i++) { 
      free(ghmm_eAA[i]); 
      free(ghmm_tAA[i]); 
    }
    free(ghmm_eAA);
    free(ghmm_tAA);
    free(ghmm_sA);
  }
  return eslOK;

 ERROR:
  *ret_tau = 0.;
  if (xv  != NULL) free(xv);
  if (dsq != NULL) free(dsq);
  if (ox  != NULL) p7_omx_Destroy(ox);
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


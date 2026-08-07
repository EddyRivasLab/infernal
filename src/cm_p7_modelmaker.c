/* cm_p7_modelmaker.c
 * EPN, Tue Aug  5 15:32:34 2008
 *
 * Construct a p7 model from CM and its CP9 HMM.
 */
#include <esl_config.h>
#include <p7_config.h>
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
#ifdef HMMER_THREADS
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "hmmer.h"

#include "infernal.h"

/* Function: BuildP7HMM_MatchEmitsOnly()
 * Incept:   EPN, Tue Aug  5 15:33:00 2008
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and its CP9 HMM.
 *           Copy only the match emissions of the CP9 HMM, the rest of 
 *           the p7 model parameters are irrelevant. 
 * 
 * Args:     cm        - the cm
 *           cp9       - the cp9 HMM to build the p7 HMM from (usually cm->cp9loc)
 *           ret_p7    - RETURN: new p7 model 
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
BuildP7HMM_MatchEmitsOnly(CM_t *cm, CP9_t *cp9, P7_HMM **ret_p7)
{
  int        status;
  P7_HMM     *hmm = NULL;        /* RETURN: new hmm */
  int        k;

  if(cp9 == NULL)         return eslEINCOMPAT; 
  if(cp9->M != cm->clen)  return eslEINCOMPAT;

  if ((hmm    = p7_hmm_Create(cm->clen, cm->abc)) == NULL)  return eslEMEM;
  if ((status = p7_hmm_Zero(hmm))                 != eslOK) return status;

  /* copy only match emissions */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cp9->mat[k], cm->abc->K, hmm->mat[k]);

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
  if((status = p7_hmm_SetConsensus(hmm, NULL)) != eslOK) goto ERROR;
  if(cm->comlog != NULL) { 
    if((status = esl_strdup(cm->comlog, -1, &(hmm->comlog))) != eslOK) goto ERROR;
  }
  else { 
    hmm->comlog = NULL;
  }

  hmm->eff_nseq = cm->eff_nseq;
  hmm->nseq     = cm->nseq;
  hmm->checksum = 0;

  *ret_p7 = hmm;

  return eslOK;

 ERROR: 
  if(hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}

/* Function: cm_cp9_to_p7()
 * Incept:   EPN, Fri Sep 24 13:46:37 2010
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and a CP9 HMM.
 * 
 * Args:     cm     - the cm, must have a cp9 model in it.
 *           cp9    - the CP9 HMM, usually cm->cp9loc
 *           errbuf - for error messages
 *
 * Return:   eslOK on success
 * Throws:   eslEINCOMPAT on contract violation, errbuf filled.
 *           eslEMEM on memory error, errbuf filled.
 */
int
cm_cp9_to_p7(CM_t *cm, CP9_t *cp9, char *errbuf)
{
  int        status;
  int        k;

  /* contract check */
  if(cp9 == NULL)            ESL_XFAIL(eslEINCOMPAT, errbuf, "trying to create a p7 from cp9 HMM, but cp9 is NULL");
  if(cm->mlp7 != NULL)       ESL_XFAIL(eslEINCOMPAT, errbuf, "trying to create ml p7, but it already exists");
  if(cm->W == 0)             ESL_XFAIL(eslEINCOMPAT, errbuf, "trying to create ml p7, cm->W is 0");
  if(cp9->M != cm->clen)     ESL_XFAIL(eslEINCOMPAT, errbuf, "trying to create ml p7, cm->clen != cp9->M");
  if(cm->cmcons == NULL)     ESL_XFAIL(eslEINCOMPAT, errbuf, "trying to create ml p7, cm->cmcons is NULL, we need it's structure string");

  if ((cm->mlp7 = p7_hmm_Create(cm->clen, cm->abc)) == NULL) ESL_XFAIL(eslEMEM, errbuf, "out of memory");
  p7_hmm_Zero(cm->mlp7);

  /* copy transitions */
  for (k = 0; k <= cm->mlp7->M; k++) { 
    cm->mlp7->t[k][p7H_MM] = cp9->t[k][CTMM];
    cm->mlp7->t[k][p7H_MI] = cp9->t[k][CTMI];
    cm->mlp7->t[k][p7H_MD] = cp9->t[k][CTMD];
    cm->mlp7->t[k][p7H_IM] = cp9->t[k][CTIM];
    cm->mlp7->t[k][p7H_II] = cp9->t[k][CTII];
    cm->mlp7->t[k][p7H_DM] = cp9->t[k][CTDM];
    cm->mlp7->t[k][p7H_DD] = cp9->t[k][CTDD];
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
  cm->mlp7->t[0][p7H_MM] = cp9->begin[1];
  esl_vec_FNorm(cm->mlp7->t[0], 3);

  /* match emissions: copy, then normalize (should be unnec actually) */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cp9->mat[k], cm->abc->K, cm->mlp7->mat[k]);
  for (k = 1; k <= cm->clen; k++) esl_vec_FNorm(cm->mlp7->mat[k], cm->abc->K);
  /* special case */
  esl_vec_FSet(cm->mlp7->mat[0], cm->mlp7->abc->K, 0.);
  cm->mlp7->mat[0][0] = 1.0;

  /* insert emissions: copy, then normalize (should be unnec actually) */
  for (k = 0; k <= cm->clen; k++) esl_vec_FCopy(cp9->ins[k], cm->abc->K, cm->mlp7->ins[k]);
  for (k = 0; k <= cm->clen; k++) esl_vec_FNorm(cm->mlp7->ins[k], cm->abc->K);

  /* copy cm->W as max_length */
  cm->mlp7->max_length = cm->W;

  p7_hmm_SetName       (cm->mlp7, cm->name);
  p7_hmm_SetAccession  (cm->mlp7, cm->acc);
  p7_hmm_SetDescription(cm->mlp7, cm->desc);
  p7_hmm_SetCtime      (cm->mlp7);
  if((status = p7_hmm_SetConsensus(cm->mlp7, NULL)) != eslOK) ESL_XFAIL(status, errbuf, "out of memory");
  if(cm->comlog != NULL) { 
    if((status = esl_strdup(cm->comlog, -1, &(cm->mlp7->comlog))) != eslOK) goto ERROR;
  }
  else { 
    cm->mlp7->comlog = NULL;
  }

  /* copy CM's RF annotation to mlp7 */
  if(cm->flags & CMH_RF && cm->rf != NULL) { 
    ESL_ALLOC(cm->mlp7->rf, sizeof(char) * (cm->clen+2));
    strcpy(cm->mlp7->rf, cm->rf);
    cm->mlp7->flags |= p7H_RF;
  }

  /* copy CM's consensus structure annotation to mlp7 */
  if(cm->cmcons != NULL) { 
    ESL_ALLOC(cm->mlp7->cs, sizeof(char) * (cm->clen+2));
    cm->mlp7->cs[0] = ' ';
    for(k = 1; k <= cm->clen; k++) { 
      cm->mlp7->cs[k] = cm->cmcons->cstr[k-1]; /* cmcons->cstr is 0..cm->clen-1, mlp7->cs is 1..cm->clen */
    }
    cm->mlp7->cs[cm->clen+1] = '\0';
    cm->mlp7->flags |= p7H_CS;
  }    

  cm->mlp7->eff_nseq = cm->eff_nseq;
  cm->mlp7->nseq     = cm->nseq;
  if(cm->flags & CMH_CHKSUM) { 
    cm->mlp7->checksum = cm->checksum;
    cm->mlp7->flags |= p7H_CHKSUM;
  }    
  else { 
    cm->mlp7->checksum = 0;
  }

  /* set the model composition */
  if ((status = p7_hmm_SetComposition(cm->mlp7)) != eslOK) ESL_XFAIL(status, errbuf, "out of memory");

  cm->flags |= CMH_MLP7; /* raise the P7 flag */

  return eslOK;

 ERROR:
  if(cm->mlp7 != NULL) { p7_hmm_Destroy(cm->mlp7); cm->mlp7 = NULL; }
  return status;
}

/* Cheap glocal-Forward lambda predictor (briefs 26_0719-053/26_0719-055,
 * deployed by brief 26_0719-054; supersedes the 2-feature 046/048/049 form).
 *
 * The glocal Forward E-value slope (GFLAMBDA) is predicted in closed form
 * from four model features instead of being reused from the local Forward
 * lambda:
 *     z0 = log(clen)                    z1 = mean_H
 *     z2 = mean_H^2                     z3 = log(min(eff_nseq, 20))
 * where mean_H is the mean per-column relative entropy (bits) of the match
 * emissions vs a uniform background (mean_relentropy_bits()) and eff_nseq is
 * the *CM's* effective sequence count (see the note in cm_p7_Calibrate()).
 * The coefficients are a z-scored OLS linear fit in log space (natural
 * log/exp), refit by brief 26_0719-055 on the combined pool of brief
 * 26_0719-053's 1053 held-out-validated multi-sequence families plus 64 real
 * nseq=1 models. Source of truth:
 * brief055_run/lambda_refit_nseq1.json -> results.H_mandatory.spec.
 *
 * Why this form: lambda error is depth-amplified (tau error is not), so
 * lambda accuracy dominates the deep tail where the real user-facing
 * `cmsearch --trmF5` E-values live. Adding mean_H^2 and the eff_nseq term
 * takes the held-out deep-tail (P=1e-8) median error to 0.354 log10 units
 * (2.3x) at only N=4 samples -- better than the old 2-feature form at
 * --EgfN 50. The 055 refit additionally brings the single-sequence (nseq=1)
 * class to population parity (5.4x -> 2.3x) at no cost to the multi-seq
 * population. mean_H^2 is quadratic in a *bounded* feature (not in clen), so
 * it does not blow up on extrapolation at large clen: brief 053 gate D
 * measured the large-clen filter-safety envelope IMPROVING, max|dS*|
 * 11.83 -> 10.78 bits.
 *
 * Kept as a clean swappable static const block: mu/sd are the training-set
 * feature means/sds, coef[0] is the intercept and coef[1..4] multiply
 * z0..z3 in order.
 */
static const double gfcalib_feat_mu[4] = {
  4.77095051518072,     /* log(clen)                */
  0.5804811877658898,   /* mean_H                   */
  0.3882958223057928,   /* mean_H^2                 */
  0.90508110760581      /* log(min(eff_nseq, 20))   */
};
static const double gfcalib_feat_sd[4] = {
  0.7059490720195424,
  0.22657760912255817,
  0.335344075259694,
  0.7676273070778664
};
static const double gfcalib_coef[5]    = {
  -0.725985532232747,    /* intercept */
  -0.2118010089398733,   /* z0 */
  -0.39601282842140223,  /* z1 */
   0.2423273564376488,   /* z2 */
  -0.06810158052837548   /* z3 */
};

/* Learned tau shrinkage correction (briefs 26_0719-053/26_0719-055, deployed
 * by brief 26_0719-054). Applied on top of the all-order-statistic raw tau
 * returned by cm_p7_Tau(); see gfcalib_shrink_tau() below for the formula and
 * the exact parametrization. Source of truth:
 * brief055_run/shrinkage_spec_combined.json -> spec_all, fit by
 * scripts/fit_task055_shrinkage_combined.py on the combined
 * {1097 multi-seq + 64 nseq=1} pool.
 */
static const double gfcalib_shrink_mu[4] = {
   5.109825312361632,    /* log(clen)                */
   0.5437359378271508,   /* mean_H                   */
   1.1712644136150876,   /* log(min(eff_nseq, 20))   */
 -10.058235004116202     /* tau_raw * lambda         */
};
static const double gfcalib_shrink_sd[4] = {
   1.0783021543068307,
   0.21815714303833725,
   0.9871147041905497,
  10.076412813149693
};
static const double gfcalib_shrink_coef[5] = {
   1.092634838132801,    /* intercept */
  -0.3867448779584632,   /* z0 */
  -0.2187994840669399,   /* z1 */
  -0.19904086095964266,  /* z2 */
  -0.5429849757828361    /* z3 */
};

/* mean_relentropy_bits()
 * Mean over the M match columns of the relative entropy (bits) of the match
 * emission distribution vs a uniform 1/K background. Matches the training
 * feature exactly (scripts/extract_p7_features.py:42,86-91): uniform 1/K
 * background (0.25 for RNA), per-column renormalize for fp roundoff, skip
 * p<=0 terms. brief 26_0719-046.
 */
static double
mean_relentropy_bits(const P7_HMM *hmm)
{
  int    k, a;
  int    K  = hmm->abc->K;
  double bg = 1.0 / (double) K;   /* uniform background, matches training (not bg->f) */
  double sum_H = 0.;

  for (k = 1; k <= hmm->M; k++) {
    double s = 0.;
    double h = 0.;
    for (a = 0; a < K; a++) s += hmm->mat[k][a];
    if (s <= 0.) continue;
    for (a = 0; a < K; a++) {
      double p = hmm->mat[k][a] / s;    /* renormalize (fp roundoff), as in training */
      if (p > 0.) h += p * (log(p / bg) / eslCONST_LOG2);  /* log2 */
    }
    sum_H += h;
  }
  return sum_H / (double) hmm->M;
}

/* gfcalib_effn_feature()
 * The eff_nseq feature shared by the lambda predictor and the tau shrinkage:
 * log(min(eff_nseq, 20)) -- exactly as fit
 * (scripts/fit_task055_shrinkage_combined.py::corr_features/lam_features_H).
 *
 * The 20.0 cap is part of the fitted form: real SEEDs saturate around
 * eff_nseq 3-17, the training pool has little support above 20 (one model sits
 * at 10000), so the cap keeps a deep alignment from extrapolating off the fit.
 *
 * There is deliberately NO floor at 1.0. Entropy weighting routinely produces
 * eff_nseq < 1 -- 51/1097 of the brief-053 multi-seq training pool (min 0.30)
 * and 20/64 of the brief-055 nseq=1 pool (min 0.72) -- and the fit consumed
 * those raw values. Flooring at 1 would silently make the deployed code
 * disagree with the fitted formula on a real minority of models. The only
 * guard is a tiny epsilon against a nonpositive eff_nseq (log domain error);
 * it is unreachable for a real CM, and gate 0 of brief 26_0719-054 checks
 * eff_nseq > 0 at the call site. brief 26_0719-054.
 */
static double
gfcalib_effn_feature(double eff_nseq)
{
  double e = eff_nseq;
  if (e < 1e-3) e = 1e-3;    /* log-domain guard only; unreachable for a real CM */
  if (e > 20.0) e = 20.0;
  return log(e);
}

/* predict_glocal_lambda()
 * Closed-form glocal Forward lambda predictor, 5-feature form
 * (briefs 26_0719-053/26_0719-055, deployed by brief 26_0719-054):
 *     lambda = exp(c0 + c1*z0 + c2*z1 + c3*z2 + c4*z3)
 * clen = hmm->M; mean_H from mean_relentropy_bits(); eff_nseq is the CM's
 * effective sequence count.
 */
static double
predict_glocal_lambda(int clen, double mean_H, double eff_nseq)
{
  double x0 = log((double) clen);           /* natural log            */
  double x1 = mean_H;
  double x2 = mean_H * mean_H;
  double x3 = gfcalib_effn_feature(eff_nseq);
  double z0 = (x0 - gfcalib_feat_mu[0]) / gfcalib_feat_sd[0];
  double z1 = (x1 - gfcalib_feat_mu[1]) / gfcalib_feat_sd[1];
  double z2 = (x2 - gfcalib_feat_mu[2]) / gfcalib_feat_sd[2];
  double z3 = (x3 - gfcalib_feat_mu[3]) / gfcalib_feat_sd[3];
  return exp(gfcalib_coef[0] + gfcalib_coef[1] * z0 + gfcalib_coef[2] * z1
                             + gfcalib_coef[3] * z2 + gfcalib_coef[4] * z3);
}

/* gfcalib_shrink_tau()
 * Learned shrinkage correction applied to the raw all-order-statistic tau
 * returned by cm_p7_Tau() (briefs 26_0719-053/26_0719-055, deployed by brief
 * 26_0719-054).
 *
 * The raw all-order average is the *worst* unshrunk estimator of the family
 * tried in brief 26_0719-041, but the best once shrunk: averaging all N
 * anchors minimizes variance, and the learned correction removes the bias
 * that averaging the low-rank anchors introduces. The last feature
 * (tau_raw*lambda) is what does the shrinking -- it partly replaces the noisy
 * N=4 sample estimate with a deterministic feature-based one.
 *
 * Parametrization matches scripts/fit_task055_shrinkage_combined.py exactly:
 * the OLS target there is d = lambda*(tau_gt - tau_raw) (a dimensionless
 * quantity), so the predicted d_hat must be divided by lambda to get a
 * correction in bits before adding it to tau_raw (fit script line 178,
 * `te_raw + apply_ols(spec, Xs) / lam`). Getting that factor of lambda wrong
 * is the easy mistake here; gate 2 of brief 26_0719-054 checks it.
 *
 * Deterministic in (features, sorted sample), so GFMU stays byte-identical
 * across --cpu.
 */
static double
gfcalib_shrink_tau(double tau_raw, double lambda, int clen, double mean_H, double eff_nseq)
{
  double x0 = log((double) clen);
  double x1 = mean_H;
  double x2 = gfcalib_effn_feature(eff_nseq);
  double x3 = tau_raw * lambda;
  double z0 = (x0 - gfcalib_shrink_mu[0]) / gfcalib_shrink_sd[0];
  double z1 = (x1 - gfcalib_shrink_mu[1]) / gfcalib_shrink_sd[1];
  double z2 = (x2 - gfcalib_shrink_mu[2]) / gfcalib_shrink_sd[2];
  double z3 = (x3 - gfcalib_shrink_mu[3]) / gfcalib_shrink_sd[3];
  double d_hat = gfcalib_shrink_coef[0] + gfcalib_shrink_coef[1] * z0
                                        + gfcalib_shrink_coef[2] * z1
                                        + gfcalib_shrink_coef[3] * z2
                                        + gfcalib_shrink_coef[4] * z3;
  return tau_raw + d_hat / lambda;    /* d_hat is in units of lambda*bits */
}

/* Function: cm_p7_Calibrate()
 * Incept:   EPN, Tue Nov  9 06:16:57 2010
 *
 * Purpose:  Calibrate a p7 HMM for local MSV, Viterbi, Forward and
 *           also glocal Forward.
 * 
 * Args:     hmm       - the hmm
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
 *           seed      - RNG seed for calibration (0=one-time arbitrary)
 *           ncpus     - number of CPUs for threaded glocal Fwd calibration (0=serial)
 *           eff_nseq  - the CM's effective sequence count (cm->eff_nseq); a
 *                       feature of both glocal Fwd (tau,lambda) predictors.
 *                       Passed in rather than read off <hmm>: see the note at
 *                       the glocal block below.
 *           ret_gfmu  - RETURN: mu for glocal forward
 *           ret_gflambda - RETURN: lambda for glocal forward
 *
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
cm_p7_Calibrate(P7_HMM *hmm, char *errbuf,
		int ElmL, int ElvL, int ElfL, int EgfL,
		int ElmN, int ElvN, int ElfN, int EgfN,
		double ElfT, double EgfT,
		int seed, int ncpus, double eff_nseq,
		double *ret_gfmu, double *ret_gflambda)
{
  int        status;
  P7_OPROFILE    *om = NULL;
  P7_BG          *bg = NULL;
  P7_PROFILE     *gm = NULL;
  ESL_RANDOMNESS *r  = NULL;
  double lmmu, lvmu, lftau, gfmu;
  double lmlam, lvlam, lflam, gflambda, lambda;

  /*printf("cm_p7_Calibrate:\n\tElmL: %d\n\tElvL: %d\n\tElfL: %d\n\tEgfL: %d\n\tElmN: %d\n\tElvN: %d\n\tElfN: %d\n\tEgfN: %d\n\tElfT: %f\n\tEgfT: %f\n\n", ElmL, ElvL, ElfL, EgfL, ElmN, ElvN, ElfN, EgfN, ElfT, EgfT, do_real, do_null3, do_fitlam, do_bias);*/

  /* most of this code stolen from hmmer's evalues.c::p7_Calibrate() */
  if (seed > 0) { if ((r = esl_randomness_CreateFast(seed)) == NULL) ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to create RNG"); }
  else          { if ((r = esl_randomness_Create(0))       == NULL) ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to create RNG"); }
  if ((bg     = p7_bg_Create(hmm->abc)) == NULL)                          ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to allocate background");
  if ((gm     = p7_profile_Create(hmm->M, hmm->abc))  == NULL)            ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to allocate profile");
  if ((status = p7_ProfileConfig(hmm, bg, gm, ElmL, p7_LOCAL)) != eslOK)  ESL_XFAIL(status,  errbuf, "cm_p7_Calibrate(): failed to configure profile");
  if ((om     = p7_oprofile_Create(hmm->M, hmm->abc)) == NULL)            ESL_XFAIL(eslEMEM, errbuf, "cm_p7_Calibrate(): failed to create optimized profile");
  if ((status = p7_oprofile_Convert(gm, om)) != eslOK)                    ESL_XFAIL(status,  errbuf, "cm_p7_Calibrate(): failed to convert to optimized profile");

  /* The calibration steps themselves */
  lambda = lmlam = lvlam = lflam = 0.;
  if ((status = p7_Lambda      (hmm, bg, &lambda))                         != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine lambda");
  if ((status = p7_MSVMu    (r, om, bg, ElmL, ElmN, lambda, &lmmu))        != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine msv mu");
  if ((status = p7_ViterbiMu(r, om, bg, ElvL, ElvN, lambda, &lvmu))        != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine vit mu");
  if ((status = p7_Tau      (r, om, bg, ElfL, ElfN, lambda, ElfT, &lftau)) != eslOK)   ESL_XFAIL(status,  errbuf, "failed to determine fwd tau");

  /* set the p7's evparam[] */
  hmm->evparam[p7_MMU]     = lmmu;
  hmm->evparam[p7_MLAMBDA] = lambda;
  hmm->evparam[p7_VMU]     = lvmu;  
  hmm->evparam[p7_VLAMBDA] = lambda;
  hmm->evparam[p7_FTAU]    = lftau; 
  hmm->evparam[p7_FLAMBDA] = lambda;
  hmm->flags              |= p7H_STATS;

  /* finally, determine Glocal Forward stats (briefs 26_0719-053/26_0719-055,
   * deployed by brief 26_0719-054).
   *
   * GFLAMBDA is predicted in closed form from (clen, mean_H, mean_H^2,
   * eff_nseq) -- NOT reused from the local Forward lambda. GFMU (tau) is then
   * the all-order-statistic average returned by cm_p7_Tau() at the predicted
   * lambda, plus a learned shrinkage correction applied here. EgfT (tailp) is
   * unused on this path; the tailp choice (0.015) is baked into the trained
   * lambda predictor.
   *
   * The shrinkage is applied here rather than inside cm_p7_Tau() so that all
   * the feature machinery lives in one place and cm_p7_Tau() keeps its
   * signature: everything the correction needs (clen, mean_H, eff_nseq, the
   * predicted lambda) is already in hand at this point.
   *
   * NOTE on <eff_nseq>: this is the *CM's* eff_nseq, passed in by the caller,
   * NOT hmm->eff_nseq. Those are genuinely different numbers -- on cmbuild's
   * default path the filter HMM's eff_nseq comes from a temporary CM built to
   * match the filter HMM's relative-entropy target (cmbuild.c
   * ::build_and_calibrate_p7_filter()), and it can differ from the CM's own
   * eff_nseq by more than 2x (e.g. Rfam 6C: CM 2.398 vs filter HMM 5.552).
   * Both predictors were trained against the EFFN on the CM header line
   * (scripts/build_task053_pool_effnseq.py, scripts/build_task055_singleseq_panel.py),
   * so reading hmm->eff_nseq here would silently feed them the wrong feature.
   */
  {
    double mean_H  = mean_relentropy_bits(hmm);
    double tau_raw;

    gflambda = predict_glocal_lambda(hmm->M, mean_H, eff_nseq);
    if ((status = p7_ProfileConfig(hmm, bg, gm, EgfL, p7_GLOCAL)) != eslOK) goto ERROR;
    if ((status = cm_p7_Tau(r, errbuf, NULL, gm, bg, EgfL, EgfN, gflambda, EgfT, ncpus, &tau_raw)) != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine fwd tau");
    gfmu = gfcalib_shrink_tau(tau_raw, gflambda, hmm->M, mean_H, eff_nseq);
  }

  esl_randomness_Destroy(r); 
  p7_bg_Destroy(bg);         
  p7_oprofile_Destroy(om);   
  p7_profile_Destroy(gm);   
  
  if(ret_gfmu != NULL)     *ret_gfmu = gfmu;
  if(ret_gflambda != NULL) *ret_gflambda = gflambda;

  return eslOK;

 ERROR: 
  if(r != NULL)  esl_randomness_Destroy(r); 
  if(bg != NULL) p7_bg_Destroy(bg);         
  if(om != NULL) p7_oprofile_Destroy(om);   
  if(gm != NULL) p7_profile_Destroy(gm);   
  if(ret_gfmu != NULL) *ret_gfmu = 0.;
  if(ret_gflambda != NULL) *ret_gflambda = 0.;
  return status;
}

/* Function:  cm_p7_GForwardScoreOnly()
 * Synopsis:  Two-row generic Forward, returning score only.
 * Incept:    EPN*, Thu Mar 19 2026
 *
 * Purpose:   Compute the Forward score for digital sequence <dsq> of
 *            length <L> against profile <gm>, using only two rows of
 *            DP memory instead of the full L x M matrix used by
 *            p7_GForward(). Only the final Forward score is returned;
 *            no DP matrix is retained.
 *
 *            Adapted from p7_GForward() (hmmer/src/generic_fwdback.c)
 *            and forward_row() (hmmer/src/generic_fwdback_chk.c).
 *
 *            This is intended for use in calibration (cm_p7_Tau()),
 *            where we need Forward scores for many random sequences
 *            but never need the full DP matrix. For large models
 *            (M=35000, L=70000), this reduces memory from ~29 GB
 *            to ~840 KB.
 *
 * Args:      dsq    - digital sequence, 1..L
 *            L      - length of dsq
 *            gm     - profile (configured for length L)
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Returns:   <eslOK> on success, <*opt_sc> is the Forward score in nats.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
cm_p7_GForwardScoreOnly(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, float *opt_sc)
{
  int          status;
  float const *tsc  = gm->tsc;
  int          M    = gm->M;
  float        esc  = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  int          rowsize = (M+1) * p7G_NSCELLS + p7G_NXCELLS;  /* MID states + specials per row */
  float       *mem  = NULL;    /* allocated memory for two rows */
  float       *prev = NULL;    /* pointer to previous row */
  float       *cur  = NULL;    /* pointer to current row */
  float       *tmp;
  int          i, k;

  /* Macros for accessing states in a flat row.
   * MID states are at row[k * p7G_NSCELLS + state].
   * Specials are at row[(M+1) * p7G_NSCELLS + special].
   */
#define ROWMX(row,k,s) ((row)[(k) * p7G_NSCELLS + (s)])
#define ROWXM(row,s)   ((row)[(M+1) * p7G_NSCELLS + (s)])

  /* NOTE (brief 26_0719-046): the p7_FLogsum() lookup table must already be
   * initialized by the caller (cm_p7_Tau() does this once in the main thread).
   * We do NOT call p7_FLogsumInit() here: it rewrites a global static table,
   * and doing so per-call would race with concurrent p7_FLogsum() reads in the
   * threaded (--cpu>0) calibration path, making scores (and thus GFMU)
   * nondeterministic across thread counts and run-to-run.
   */

  ESL_ALLOC(mem, sizeof(float) * 2 * rowsize);
  prev = mem;
  cur  = mem + rowsize;

  /* Initialization of row 0 */
  for (k = 0; k <= M; k++)
    ROWMX(prev, k, p7G_M) = ROWMX(prev, k, p7G_I) = ROWMX(prev, k, p7G_D) = -eslINFINITY;
  ROWXM(prev, p7G_N) = 0;
  ROWXM(prev, p7G_B) = gm->xsc[p7P_N][p7P_MOVE];
  ROWXM(prev, p7G_E) = ROWXM(prev, p7G_C) = ROWXM(prev, p7G_J) = -eslINFINITY;

  /* Recursion */
  for (i = 1; i <= L; i++)
    {
      float const *rsc = gm->rsc[dsq[i]];
      float sc;

      ROWMX(cur, 0, p7G_M) = ROWMX(cur, 0, p7G_I) = ROWMX(cur, 0, p7G_D) = -eslINFINITY;
      ROWXM(cur, p7G_E) = -eslINFINITY;

      for (k = 1; k < M; k++)
	{
	  /* match state */
	  sc = p7_FLogsum(p7_FLogsum(ROWMX(prev,k-1,p7G_M) + TSC(p7P_MM,k-1),
				     ROWMX(prev,k-1,p7G_I) + TSC(p7P_IM,k-1)),
			  p7_FLogsum(ROWMX(prev,k-1,p7G_D) + TSC(p7P_DM,k-1),
				     ROWXM(prev,p7G_B)      + TSC(p7P_BM,k-1)));
	  ROWMX(cur, k, p7G_M) = sc + MSC(k);

	  /* insert state */
	  sc = p7_FLogsum(ROWMX(prev,k,p7G_M) + TSC(p7P_MI,k),
			  ROWMX(prev,k,p7G_I) + TSC(p7P_II,k));
	  ROWMX(cur, k, p7G_I) = sc + ISC(k);

	  /* delete state */
	  ROWMX(cur, k, p7G_D) = p7_FLogsum(ROWMX(cur,k-1,p7G_M) + TSC(p7P_MD,k-1),
					     ROWMX(cur,k-1,p7G_D) + TSC(p7P_DD,k-1));

	  /* E state update */
	  ROWXM(cur, p7G_E) = p7_FLogsum(p7_FLogsum(ROWMX(cur,k,p7G_M) + esc,
						     ROWMX(cur,k,p7G_D) + esc),
					  ROWXM(cur, p7G_E));
	}

      /* unrolled match state M_M */
      sc = p7_FLogsum(p7_FLogsum(ROWMX(prev,M-1,p7G_M) + TSC(p7P_MM,M-1),
				 ROWMX(prev,M-1,p7G_I) + TSC(p7P_IM,M-1)),
		      p7_FLogsum(ROWMX(prev,M-1,p7G_D) + TSC(p7P_DM,M-1),
				 ROWXM(prev,p7G_B)      + TSC(p7P_BM,M-1)));
      ROWMX(cur, M, p7G_M) = sc + MSC(M);
      ROWMX(cur, M, p7G_I) = -eslINFINITY;

      /* unrolled delete state D_M */
      ROWMX(cur, M, p7G_D) = p7_FLogsum(ROWMX(cur,M-1,p7G_M) + TSC(p7P_MD,M-1),
					 ROWMX(cur,M-1,p7G_D) + TSC(p7P_DD,M-1));

      /* unrolled E state update */
      ROWXM(cur, p7G_E) = p7_FLogsum(p7_FLogsum(ROWMX(cur,M,p7G_M),
						 ROWMX(cur,M,p7G_D)),
				      ROWXM(cur, p7G_E));

      /* J state */
      ROWXM(cur, p7G_J) = p7_FLogsum(ROWXM(prev, p7G_J) + gm->xsc[p7P_J][p7P_LOOP],
				      ROWXM(cur,  p7G_E) + gm->xsc[p7P_E][p7P_LOOP]);
      /* C state */
      ROWXM(cur, p7G_C) = p7_FLogsum(ROWXM(prev, p7G_C) + gm->xsc[p7P_C][p7P_LOOP],
				      ROWXM(cur,  p7G_E) + gm->xsc[p7P_E][p7P_MOVE]);
      /* N state */
      ROWXM(cur, p7G_N) = ROWXM(prev, p7G_N) + gm->xsc[p7P_N][p7P_LOOP];

      /* B state */
      ROWXM(cur, p7G_B) = p7_FLogsum(ROWXM(cur, p7G_N) + gm->xsc[p7P_N][p7P_MOVE],
				      ROWXM(cur, p7G_J) + gm->xsc[p7P_J][p7P_MOVE]);

      /* swap rows */
      tmp = prev; prev = cur; cur = tmp;
    }

  /* after the swap, prev holds the final row L */
  if (opt_sc != NULL) *opt_sc = ROWXM(prev, p7G_C) + gm->xsc[p7P_C][p7P_MOVE];

  free(mem);

#undef ROWMX
#undef ROWXM
  return eslOK;

 ERROR:
  if (mem != NULL) free(mem);
  if (opt_sc != NULL) *opt_sc = 0.;
  return status;
}


/* Structure for passing work units through the work queue
 * in the threaded glocal Forward calibration path.
 * A small pool of these structs cycles between reader and workers.
 */
typedef struct {
  ESL_DSQ *dsq;     /* digital sequence buffer, 1..L (owned by this struct) */
  int      L;       /* sequence length; 0 = sentinel (stop signal) */
  double   sc;      /* RETURN: bit score (fwd - null) / log2 */
  int      idx;     /* sequence index in xv[] array; -1 for sentinel */
} CM_P7_TAU_WORK;

#ifdef HMMER_THREADS
/* Per-worker data for threaded glocal Forward calibration.
 * Each worker accumulates scores in its own local scA[] array,
 * which the main thread merges after all workers finish.
 * This follows the cmcalibrate pattern.
 */
typedef struct {
  P7_PROFILE      *gm;
  P7_BG           *bg;
  ESL_WORK_QUEUE  *queue;
  double          *scA;     /* worker-local score array, pre-allocated by main thread */
  int              nsc;     /* number of scores stored in scA */
} CM_P7_TAU_WINFO;

/* cm_p7_tau_thread_worker()
 * Worker function for threaded glocal Forward calibration.
 * Each worker pulls sequences from the queue, runs
 * cm_p7_GForwardScoreOnly(), and accumulates bit scores
 * in its own local scA[] array.
 */
static void
cm_p7_tau_thread_worker(void *arg)
{
  ESL_THREADS      *obj = (ESL_THREADS *) arg;
  int               workeridx;
  CM_P7_TAU_WINFO  *winfo;
  CM_P7_TAU_WORK   *work = NULL;
  void             *newwork;

  esl_threads_Started(obj, &workeridx);
  winfo = (CM_P7_TAU_WINFO *) esl_threads_GetData(obj, workeridx);

  winfo->nsc = 0;

  esl_workqueue_WorkerUpdate(winfo->queue, NULL, &newwork);
  work = (CM_P7_TAU_WORK *) newwork;

  while (work->L > 0)   /* sentinel: L==0 means stop */
    {
      float fsc, nullsc;

      cm_p7_GForwardScoreOnly(work->dsq, work->L, winfo->gm, &fsc);
      p7_bg_NullOne(winfo->bg, work->dsq, work->L, &nullsc);

      winfo->scA[winfo->nsc] = (double)((fsc - nullsc) / eslCONST_LOG2);
      winfo->nsc++;

      esl_workqueue_WorkerUpdate(winfo->queue, work, &newwork);
      work = (CM_P7_TAU_WORK *) newwork;
    }
  esl_workqueue_WorkerUpdate(winfo->queue, work, NULL);
  esl_threads_Finished(obj, workeridx);
}
#endif /* HMMER_THREADS */


/* Function:  cm_p7_Tau()
 * Synopsis:  Determine Forward tau by brief simulation.
 * Incept:    SRE, Thu Aug  9 15:08:39 2007 [Janelia] (p7_Tau())
 *
 * Purpose:   Identical to p7_Tau() except that it can handle
 *            either an optimized profile or a generic profile,
 *            the latter of which is used for glocal Forward.
 *            See hmmer/evalues.c::cm_p7_Tau for additional information.
 *
 *            When <ncpus> > 0 and <gm> is non-NULL (generic/glocal
 *            path), the N Forward evaluations are parallelized
 *            across <ncpus> worker threads. Sequences are
 *            pre-generated by the main thread for reproducibility.
 *
 * Args:      r      : source of randomness
 *            errbuf : for error messages
 *            om     : configured profile (optimized), if non-NULL, <gm> must be NULL
 *            gm     : configured profile (generic),   if non-NULL, <om> must be NULL
 *            bg     : null model (for background residue frequencies)
 *            L      : mean length model for seq emission from profile
 *            N      : number of sequences to generate
 *            lambda : expected slope of the exponential tail (from p7_Lambda())
 *            tailp  : tail mass from which we will extrapolate mu
 *            ncpus  : number of CPUs for threaded glocal Fwd (0=serial)
 *            ret_tau : RETURN: estimate for the Forward tau (base of
 *                      exponential tail). On the glocal path this is the RAW
 *                      all-order-statistic tau -- the caller is expected to
 *                      apply gfcalib_shrink_tau() to it to get the final
 *                      GFMU. See brief 26_0719-054.
 *
 * Returns:   <eslOK> on success, and <*ret_tau> is the tau estimate.
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_tau> is 0.
 */
int
cm_p7_Tau(ESL_RANDOMNESS *r, char *errbuf, P7_OPROFILE *om, P7_PROFILE *gm, P7_BG *bg, int L, int N, double lambda, double tailp, int ncpus, double *ret_tau)
{
  P7_OMX  *ox = NULL;

  ESL_DSQ *dsq     = NULL;
  double  *xv      = NULL;
  float    fsc, nullsc;
  int      status;
  int      i;
  int do_generic;

  if(om == NULL && gm == NULL) { status = eslEINVAL; goto ERROR; }
  if(om != NULL && gm != NULL) { status = eslEINVAL; goto ERROR; }
  do_generic = (gm != NULL) ? TRUE : FALSE;

  ESL_ALLOC(xv,  sizeof(double)  * N);

  /* Initialize the global p7_FLogsum() lookup table ONCE, here in the main
   * thread, before any worker scores a sequence (brief 26_0719-046). The
   * per-call init previously inside cm_p7_GForwardScoreOnly() raced with
   * concurrent reads under --cpu>0 and made GFMU nondeterministic. It writes
   * the same deterministic values every time, so a single up-front init makes
   * the threaded and serial scoring paths bit-identical. */
  p7_FLogsumInit();

  if(do_generic) p7_ReconfigLength(gm, L);
  else           p7_oprofile_ReconfigLength(om, L);
  p7_bg_SetLength(bg, L);

#ifdef HMMER_THREADS
  /* Threaded path for generic (glocal) Forward calibration.
   * Pattern: main thread pre-generates all N sequences, then
   * uses a small pool of recycling work items to feed workers.
   */
  if(do_generic && ncpus > 0)
    {
      ESL_THREADS      *threadObj = NULL;
      ESL_WORK_QUEUE   *queue     = NULL;
      CM_P7_TAU_WINFO  *winfo    = NULL;
      CM_P7_TAU_WORK   *wpool    = NULL;     /* small recycling pool */
      ESL_DSQ         **dsqpool   = NULL;     /* all N pre-generated sequences */
      int               npool;                /* size of recycling pool */
      int               next_seq;             /* next sequence to assign */
      int               sentinels_sent;
      int               j;
      void             *newptr;
      CM_P7_TAU_WORK   *work;

      /* 1. Pre-generate all N sequences deterministically */
      ESL_ALLOC(dsqpool, sizeof(ESL_DSQ *) * N);
      for (i = 0; i < N; i++) {
	ESL_ALLOC(dsqpool[i], sizeof(ESL_DSQ) * (L+2));
	if((status = esl_rsq_xfIID(r, bg->f, bg->abc->K, L, dsqpool[i])) != eslOK) goto ERROR;
      }

      /* 2. Create recycling pool of work items */
      npool = ncpus * 2;
      ESL_ALLOC(wpool, sizeof(CM_P7_TAU_WORK) * npool);
      for (j = 0; j < npool; j++) {
	wpool[j].dsq = NULL;
	wpool[j].L   = 0;
	wpool[j].sc  = 0.;
	wpool[j].idx = -1;
      }

      /* 3. Set up threads and work queue */
      threadObj = esl_threads_Create(&cm_p7_tau_thread_worker);
      queue     = esl_workqueue_Create(npool);

      ESL_ALLOC(winfo, sizeof(CM_P7_TAU_WINFO) * ncpus);
      for (j = 0; j < ncpus; j++) {
	winfo[j].gm    = gm;
	winfo[j].bg    = bg;
	winfo[j].queue = queue;
	winfo[j].nsc   = 0;
	ESL_ALLOC(winfo[j].scA, sizeof(double) * N);  /* pre-alloc to max possible */
	esl_threads_AddThread(threadObj, &winfo[j]);
      }

      /* 4. Initialize queue with empty work items */
      for (j = 0; j < npool; j++)
	esl_workqueue_Init(queue, &wpool[j]);

      /* 5. Reader loop: fill work items with pre-generated sequences
       * and push to workers. Workers accumulate scores in their
       * own local scA[] arrays.
       */
      esl_workqueue_Reset(queue);
      esl_threads_WaitForStart(threadObj);

      next_seq       = 0;
      sentinels_sent = 0;

      /* get first empty work item */
      status = esl_workqueue_ReaderUpdate(queue, NULL, &newptr);
      if (status != eslOK) goto ERROR;
      work = (CM_P7_TAU_WORK *) newptr;

      while (sentinels_sent < ncpus)
	{
	  /* fill work item with next sequence, or make it a sentinel */
	  if (next_seq < N) {
	    work->dsq = dsqpool[next_seq];
	    work->L   = L;
	    work->idx = next_seq;
	    next_seq++;
	  } else {
	    work->dsq = NULL;
	    work->L   = 0;    /* sentinel */
	    work->idx = -1;
	    sentinels_sent++;
	  }

	  /* send filled/sentinel item, get back a recycled one */
	  status = esl_workqueue_ReaderUpdate(queue, work, &newptr);
	  if (status != eslOK) goto ERROR;
	  work = (CM_P7_TAU_WORK *) newptr;
	}

      esl_threads_WaitForFinish(threadObj);
      esl_workqueue_Complete(queue);

      /* 6. Merge per-worker scores into xv[], following cmcalibrate pattern */
      {
	int n = 0;
	for (j = 0; j < ncpus; j++) {
	  for (i = 0; i < winfo[j].nsc; i++)
	    xv[n++] = winfo[j].scA[i];
	}
      }

      /* 7. Cleanup thread resources */
      for (i = 0; i < N; i++)
	free(dsqpool[i]);
      free(dsqpool);
      free(wpool);
      for (j = 0; j < ncpus; j++)
	free(winfo[j].scA);
      free(winfo);
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
  else
#endif /* HMMER_THREADS */
    {
      /* Serial path (original behavior) */
      if(! do_generic) {
	ox = p7_omx_Create(om->M, 0, L);
	if (ox == NULL) { status = eslEMEM; goto ERROR; }
      }

      ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

      for (i = 0; i < N; i++)
	{
	  if((status = esl_rsq_xfIID(r, bg->f, bg->abc->K, L, dsq)) != eslOK) goto ERROR;
	  if(do_generic) {
	    if ((status = cm_p7_GForwardScoreOnly(dsq, L, gm, &fsc))   != eslOK) goto ERROR;
	  }
	  else {
	    if ((status = p7_ForwardParser(dsq, L, om, ox, &fsc))      != eslOK) goto ERROR;
	  }
	  if((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;
	  /* keep full double precision (match the threaded worker exactly, no
	   * intermediate float rounding) so serial and threaded xv[] -- and thus
	   * GFMU -- are bit-identical across --cpu (brief 26_0719-046). */
	  xv[i] = (double)((fsc - nullsc) / eslCONST_LOG2);
	}

      free(dsq); dsq = NULL;
      if (ox != NULL) { p7_omx_Destroy(ox); ox = NULL; }
    }

  /* known-lambda all-order-statistic tau estimator (briefs 26_0719-041/
   * 26_0719-053, deployed by brief 26_0719-054; replaces the top-half
   * estimator of brief 26_0719-046). Given the N sampled bit scores and the
   * *predicted* glocal lambda passed in <lambda>, sort ascending and use the
   * distribution-free order-statistic identity E[S(X_(k))] = (N-k+1)/(N+1)
   * under an exponential tail S(x)=exp(-lambda*(x-tau)):
   *     tau_hat_(k) = X_(k) + log((N-k+1)/(N+1)) / lambda
   * then average the anchors over ALL k = 1..N (brief 26_0719-041's
   * "known_allavg"; matches Python tau_k.mean()).
   *
   * On its own this is the *worst* raw estimator of the family -- averaging
   * in the low-rank anchors is biased -- but it has the lowest variance, and
   * it is the best of the family once the caller applies the learned
   * shrinkage correction that removes that bias (gfcalib_shrink_tau()). So
   * what this function returns is a raw tau, not the final GFMU.
   *
   * The sort makes the result independent of the threaded worker merge order,
   * so (tau,lambda) is byte-identical across --cpu. <tailp> (EgfT) is
   * intentionally UNUSED here -- the tailp choice is baked into the predicted
   * lambda.
   */
  esl_vec_DSortIncreasing(xv, N);
  {
    double tau_sum = 0.;
    int    k;
    for (k = 0; k < N; k++) {            /* k is 0-based rank */
      double p_k = (double) (N - k) / (double) (N + 1);   /* (N-(k+1)+1)/(N+1) */
      tau_sum += xv[k] + log(p_k) / lambda;
    }
    *ret_tau = tau_sum / (double) N;
  }

  free(xv);
  return eslOK;

 ERROR:
  *ret_tau = 0.;
  if (xv  != NULL) free(xv);
  if (dsq != NULL) free(dsq);
  if (ox  != NULL) p7_omx_Destroy(ox);
  return status;
}

/* Function: cm_SetFilterHMM()
 * Incept:   EPN, Mon Dec 27 07:59:47 2010
 * 
 * Purpose:  Assign a p7 HMM as a CM's filter hmm (cm->fp7)
 * 
 * Args:     cm       - the CM
 *           hmm      - the HMM to add
 *           gfmu     - glocal forward mu parameter
 *           gflambda - glocal forward lambda parameter
 *           errbuf   - for error messages 
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
cm_SetFilterHMM(CM_t *cm, P7_HMM *hmm, double gfmu, double gflambda)
{
  if(cm->fp7 != NULL) { 
    p7_hmm_Destroy(cm->fp7);
  }
  cm->fp7 = hmm;

  if(hmm->flags & p7H_STATS) {
    cm->fp7_evparam[CM_p7_LMMU]     = hmm->evparam[p7_MMU];
    cm->fp7_evparam[CM_p7_LMLAMBDA] = hmm->evparam[p7_MLAMBDA];
    cm->fp7_evparam[CM_p7_LVMU]     = hmm->evparam[p7_VMU];
    cm->fp7_evparam[CM_p7_LVLAMBDA] = hmm->evparam[p7_VLAMBDA];
    cm->fp7_evparam[CM_p7_LFTAU]    = hmm->evparam[p7_FTAU];
    cm->fp7_evparam[CM_p7_LFLAMBDA] = hmm->evparam[p7_FLAMBDA];
    cm->fp7_evparam[CM_p7_GFMU]     = gfmu;
    cm->fp7_evparam[CM_p7_GFLAMBDA] = gflambda;
  }    
  else { /* this should never happen */
    cm->fp7_evparam[CM_p7_LMMU]     = 0.;
    cm->fp7_evparam[CM_p7_LMLAMBDA] = 0.;
    cm->fp7_evparam[CM_p7_LVMU]     = 0.;
    cm->fp7_evparam[CM_p7_LVLAMBDA] = 0.;
    cm->fp7_evparam[CM_p7_LFTAU]    = 0.;
    cm->fp7_evparam[CM_p7_LFLAMBDA] = 0.;
    cm->fp7_evparam[CM_p7_GFMU]     = 0.;
    cm->fp7_evparam[CM_p7_GFLAMBDA] = 0.;
  }
  cm->flags |= CMH_FP7; /* raise the FP7 flag */

  return eslOK;
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


/* Function: cm_p7_hmm_Sizeof()
 * Incept:   EPN, Wed Jan 18 10:10:10 2012
 * 
 * Purpose:  Calculate and return size of a P7_HMM
 *           in Mb.
 * 
 * Args:     hmm       - the p7 HMM
 *           
 * Return:   size of hmm in Mb
 */
float
cm_p7_hmm_Sizeof(P7_HMM *hmm)
{
  float bytes = 0.;

  bytes = sizeof(P7_HMM);

  if(hmm->M > 0 && hmm->abc != NULL) 
  /* following from p7_hmm_CreateBody() */
  bytes += sizeof(float *) * (hmm->M+1); /* t */
  bytes += sizeof(float *) * (hmm->M+1); /* mat */
  bytes += sizeof(float *) * (hmm->M+1); /* ins */

  bytes += sizeof(float) * p7H_NTRANSITIONS*(hmm->M+1); /* t */
  bytes += sizeof(float *) * hmm->abc->K * (hmm->M+1); /* mat */
  bytes += sizeof(float *) * hmm->abc->K * (hmm->M+1); /* ins */

  if(hmm->rf        != NULL) bytes += sizeof(char) * (hmm->M+2); 
  if(hmm->consensus != NULL) bytes += sizeof(char) * (hmm->M+2); 
  if(hmm->cs        != NULL) bytes += sizeof(char) * (hmm->M+2); 
  if(hmm->ca        != NULL) bytes += sizeof(char) * (hmm->M+2); 
  if(hmm->map       != NULL) bytes += sizeof(int)  * (hmm->M+1); 

  return bytes / 1000000.;
}

/* Function:  cm_p7_hmm_SetConsensus()
 * Incept:    EPN, Wed May  9 14:13:37 2012 
 * Synopsis:  Set the consensus residue line of the HMM.
 *
 * Purpose:   Sets the consensus annotation line of the model <hmm>.
 *            
 *            Based on p7_hmm_SetConsensus() which is flexible to
 *            setting the consensus as a single sequence or a
 *            consensus from a multiple sequence alignment.  Here,
 *            only the latter case is handled, i.e. but in the future,
 *            we should relax this to allow for single sequence
 *            models.
 *
 *            This function only exists because p7_hmm_SetConsensus()
 *            uses a threshold probability of 0.9 for setting
 *            a consensus residue as uppercase, while we want 
 *            to be able to use 0.5 since that's what we use with
 *            single stranded CM positions. (Actually we use 1.0 
 *            bits, which equates to a 0.5 probability for a default
 *            null1 model (so using 0.5 will be wrong for non-standard
 *            null models...)).
 *
 *            The most likely (highest emission probability) residue
 *            is the consensus at each position.  If the emission
 *            probability is $\geq$ certain threshold (0.5), the
 *            residue is upper cased.
 *            
 * Args:      hmm - model with valid probability parameters mat[1..M][x]
 *           
 * Returns:   <eslOK> on success. The <p7H_CONS> flag on the <hmm> is raised
 *            if it wasn't already. The <hmm->consensus> line is set.
 *
 * Throws:    <eslEMEM> on allocation error. The <p7H_CONS> is dropped, even
 *            if it was up to begin with, and the <hmm->consensus> is <NULL>,
 *            even if we had one to begin with.
 *
 */
int
cm_p7_hmm_SetConsensus(P7_HMM *hmm)
{
  int   k, x;
  float mthresh = 0.5;
  int   status;
  
  /* allocation, if needed */
  if (! hmm->consensus) ESL_ALLOC(hmm->consensus, sizeof(char) * (hmm->M+2));

  /* set our arbitrary threshold for upper/lower casing */

  hmm->consensus[0] = ' ';
  for (k = 1; k <= hmm->M; k++) 
    {
      x = esl_vec_FArgMax(hmm->mat[k], hmm->abc->K);
      hmm->consensus[k] = ((hmm->mat[k][x] >= mthresh) ? toupper(hmm->abc->sym[x]) : tolower(hmm->abc->sym[x]));
    }
  hmm->consensus[hmm->M+1] = '\0';
  hmm->flags  |= p7H_CONS;	
  return eslOK;

 ERROR:
  if (hmm->consensus) free(hmm->consensus);
  hmm->consensus = NULL;
  hmm->flags    &= (~p7H_CONS);	
  return status;
}

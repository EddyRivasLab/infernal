/* cm_filtercutoff.c
 *
 * Per-CM F1/F2/F3 P-value cutoff calibration.
 *
 * cm_CalibrateFilterPvalCutoffs() emits N sequences from a CM, runs the
 * F1 (MSV), F2 (Viterbi), and F3 (local Forward) HMM filter stages on each
 * emission, and picks the M-th-rank P-value at each stage as the per-CM
 * cutoff (M = ceil(0.99*N), giving 99% retention of CM-emitted sequences).
 * Floor (loosest = global default) and ceiling (tightest = default/10) are
 * applied; monotonicity F1 >= F2 >= F3 is enforced. Results are stored on
 * the CM struct (cm->F{1,2,3}_pcutoff) with CMH_FILTER_PVAL_CUTOFFS raised.
 *
 * EPN, 2026-04-28
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

#define CM_FILTER_PVAL_RETENTION    0.99   /* 99% of CM-emitted seqs survive */

/* Logistic CLEN-scaled ceiling factor (used at SEARCH TIME by
 * cm_pli_NewModel(), not at calibration time). factor(clen) =
 *   1 + 29 / (1 + exp(-5 * (log10(clen) - 2.666)))
 * Passes through 50→1.2, 100→2, 200→5, 463→15, 3000→29.5; asymptotes at 30.
 *
 * This stays here (rather than cm_pipeline.c) so the curve definition
 * and the cmbuild-time emission code are in the same translation unit
 * if/when calibration ever wants to know it.
 */
double
cm_filter_ceiling_factor_clen(int clen)
{
  double lc = log10((double) clen);
  double f  = 1.0 + 29.0 / (1.0 + exp(-5.0 * (lc - 2.666)));
  if (f < 1.0) f = 1.0;
  return f;
}

static int
cmp_double_asc(const void *a, const void *b)
{
  double da = *(const double *) a;
  double db = *(const double *) b;
  if (da < db) return -1;
  if (da > db) return  1;
  return 0;
}

/* Function: cm_CalibrateFilterPvalCutoffs()
 * Synopsis: Compute per-CM F1/F2/F3 P-value cutoffs from N CM emissions.
 *
 * Purpose:  Emit <N> sequences from <cm>, run the local-mode F1 (MSV),
 *           F2 (Viterbi), and F3 (local Forward) HMM filter stages on each,
 *           then pick the M-th-rank P-value at each stage (M = ceil(0.99*N))
 *           as the per-CM cutoff. Store the **raw** quantile in
 *           cm->F{1,2,3}_pcutoff and raise CMH_FILTER_PVAL_CUTOFFS.
 *
 *           Floor (= pli->F*_orig) and ceiling (= pli->F*_orig / factor(clen))
 *           are NOT applied here — they're applied at SEARCH TIME in
 *           cm_pli_NewModel() where pli->F*_orig is known (depends on the
 *           Z-tier the user chose with --FZ). Storing the raw quantile keeps
 *           the CM file invariant under future search-time policy tweaks
 *           (different ceiling factor curves, different Z tiers, etc.).
 *           Likewise F1≥F2≥F3 monotonicity is re-enforced at search time
 *           after clamping; raw quantiles can violate monotonicity.
 *
 *           Caller must ensure <cm> has a valid cm->fp7 (CMH_FP7) with
 *           local-mode evparam fields (CM_p7_LMMU/LMLAMBDA/LVMU/LVLAMBDA/
 *           LFTAU/LFLAMBDA) set.
 *
 * Args:     cm     - CM with fp7 attached
 *           r      - random number generator (drives EmitParsetree)
 *           N      - number of sequences to emit (e.g. 1000)
 *           errbuf - for error messages
 *
 * Returns:  <eslOK> on success.
 *           <eslEINVAL> if cm->fp7 is missing (with errbuf message).
 *           <eslEMEM> on allocation failure.
 */
int
cm_CalibrateFilterPvalCutoffs(CM_t *cm, ESL_RANDOMNESS *r, int N,
                              char *errbuf)
{
  int           status;
  P7_BG        *bg   = NULL;
  P7_PROFILE   *gm   = NULL;
  P7_OPROFILE  *om   = NULL;
  P7_OMX       *oxf  = NULL;
  ESL_SQ       *sq   = NULL;
  Parsetree_t  *tr   = NULL;
  double       *F1ps = NULL;
  double       *F2ps = NULL;
  double       *F3ps = NULL;
  int           i;
  int           M_idx;
  double        F1cut, F2cut, F3cut;

  if (cm->fp7 == NULL || ! (cm->flags & CMH_FP7))
    ESL_XFAIL(eslEINVAL, errbuf, "cm_CalibrateFilterPvalCutoffs: CM has no valid fp7 filter HMM");
  if (N < 10)
    ESL_XFAIL(eslEINVAL, errbuf, "cm_CalibrateFilterPvalCutoffs: N must be >= 10 (got %d)", N);

  ESL_ALLOC(F1ps, sizeof(double) * N);
  ESL_ALLOC(F2ps, sizeof(double) * N);
  ESL_ALLOC(F3ps, sizeof(double) * N);

  bg = p7_bg_Create(cm->abc);
  gm = p7_profile_Create(cm->fp7->M, cm->abc);
  om = p7_oprofile_Create(cm->fp7->M, cm->abc);
  /* oxf must be created with M and a target length; reuse one matrix
   * across all emissions; resize per-sequence via p7_omx_GrowTo if needed.
   */
  oxf = p7_omx_Create(cm->fp7->M, 0, 0);
  if (bg == NULL || gm == NULL || om == NULL || oxf == NULL)
    ESL_XFAIL(eslEMEM, errbuf, "cm_CalibrateFilterPvalCutoffs: alloc failure");

  /* Local p7 profile, matching what F1/F2/F3 use in the search pipeline. */
  if ((status = p7_ProfileConfig(cm->fp7, bg, gm, 100, p7_LOCAL))     != eslOK)
    ESL_XFAIL(status, errbuf, "p7_ProfileConfig failed");
  if ((status = p7_oprofile_Convert(gm, om))                          != eslOK)
    ESL_XFAIL(status, errbuf, "p7_oprofile_Convert failed");

  for (i = 0; i < N; i++)
    {
      float  mfsc, vfsc, fwdsc, nullsc;
      char   name[32];
      int    L;

      snprintf(name, sizeof(name), "emit%d", i);
      if ((status = EmitParsetree(cm, errbuf, r, name, TRUE, &tr, &sq, NULL)) != eslOK) goto ERROR;
      L = sq->n;

      p7_oprofile_ReconfigLength(om, L);
      p7_bg_SetLength(bg, L);
      p7_bg_NullOne(bg, sq->dsq, L, &nullsc);

      if ((status = p7_omx_GrowTo(oxf, om->M, 0, L)) != eslOK)
        ESL_XFAIL(status, errbuf, "p7_omx_GrowTo failed at emission %d", i);

      /* F1: MSV. eslERANGE means score overflowed — P-value essentially 0. */
      status = p7_MSVFilter(sq->dsq, L, om, oxf, &mfsc);
      if (status == eslERANGE)      F1ps[i] = 0.0;
      else if (status != eslOK)     ESL_XFAIL(status, errbuf, "p7_MSVFilter failed at emission %d (status %d)", i, status);
      else                          F1ps[i] = esl_gumbel_surv((mfsc - nullsc) / eslCONST_LOG2,
                                                              cm->fp7_evparam[CM_p7_LMMU],
                                                              cm->fp7_evparam[CM_p7_LMLAMBDA]);

      /* F2: Viterbi. eslERANGE handled the same way. */
      status = p7_ViterbiFilter(sq->dsq, L, om, oxf, &vfsc);
      if (status == eslERANGE)      F2ps[i] = 0.0;
      else if (status != eslOK)     ESL_XFAIL(status, errbuf, "p7_ViterbiFilter failed at emission %d (status %d)", i, status);
      else                          F2ps[i] = esl_gumbel_surv((vfsc - nullsc) / eslCONST_LOG2,
                                                              cm->fp7_evparam[CM_p7_LVMU],
                                                              cm->fp7_evparam[CM_p7_LVLAMBDA]);

      /* F3: local Forward (parser). Forward uses double-precision, so no
       * range overflow expected, but be defensive.
       */
      status = p7_ForwardParser(sq->dsq, L, om, oxf, &fwdsc);
      if (status == eslERANGE)      F3ps[i] = 0.0;
      else if (status != eslOK)     ESL_XFAIL(status, errbuf, "p7_ForwardParser failed at emission %d (status %d)", i, status);
      else                          F3ps[i] = esl_exp_surv((fwdsc - nullsc) / eslCONST_LOG2,
                                                           cm->fp7_evparam[CM_p7_LFTAU],
                                                           cm->fp7_evparam[CM_p7_LFLAMBDA]);
      status = eslOK;

      FreeParsetree(tr); tr = NULL;
      esl_sq_Destroy(sq); sq = NULL;
    }

  /* M-th-rank: 99%-quantile of P-values per stage.
   * Sort ascending; index = ceil(0.99*N) - 1 picks the value such that
   * 99% of P-values are <= that value (i.e., 99% retention).
   */
  qsort(F1ps, N, sizeof(double), cmp_double_asc);
  qsort(F2ps, N, sizeof(double), cmp_double_asc);
  qsort(F3ps, N, sizeof(double), cmp_double_asc);

  M_idx = (int) ceil(CM_FILTER_PVAL_RETENTION * (double) N) - 1;
  if (M_idx < 0)  M_idx = 0;
  if (M_idx >= N) M_idx = N - 1;

  /* Store the RAW M-th-rank P-value at each stage. Floor, ceiling, and
   * F1>=F2>=F3 monotonicity are applied at SEARCH TIME (cm_pli_NewModel)
   * using pli->F*_orig (Z-tier-dependent) so policy tracks --FZ.
   */
  F1cut = F1ps[M_idx];
  F2cut = F2ps[M_idx];
  F3cut = F3ps[M_idx];

  cm->F1_pcutoff = (float) F1cut;
  cm->F2_pcutoff = (float) F2cut;
  cm->F3_pcutoff = (float) F3cut;
  cm->flags |= CMH_FILTER_PVAL_CUTOFFS;

  free(F1ps); free(F2ps); free(F3ps);
  if (oxf != NULL) p7_omx_Destroy(oxf);
  if (om  != NULL) p7_oprofile_Destroy(om);
  if (gm  != NULL) p7_profile_Destroy(gm);
  if (bg  != NULL) p7_bg_Destroy(bg);
  return eslOK;

 ERROR:
  if (F1ps != NULL) free(F1ps);
  if (F2ps != NULL) free(F2ps);
  if (F3ps != NULL) free(F3ps);
  if (tr   != NULL) FreeParsetree(tr);
  if (sq   != NULL) esl_sq_Destroy(sq);
  if (oxf  != NULL) p7_omx_Destroy(oxf);
  if (om   != NULL) p7_oprofile_Destroy(om);
  if (gm   != NULL) p7_profile_Destroy(gm);
  if (bg   != NULL) p7_bg_Destroy(bg);
  return status;
}

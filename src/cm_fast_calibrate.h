/* cm_fast_calibrate.h
 * Fast CM calibration via ridge regression on compiled-in JSON models.
 *
 * Phase 1 stub: cm_FastCalibrate() returns eslFAIL with no side effects.
 * JSON parsing and model loading are implemented in Phase 1; real prediction
 * lands in Phase 4.
 *
 * Phase 2: feature extraction (cm_FastCalibrate_ExtractFeatures and friends).
 *
 * See CPORT_SPEC.md for the full specification.
 */
#ifndef cm_FAST_CALIBRATE_INCLUDED
#define cm_FAST_CALIBRATE_INCLUDED
#include "infernal.h"  /* CM_t */

/* cm_FastCalibrate()
 *   Compute lambda, mu_extrap, mu_orig for the 4 ECM modes via ridge prediction.
 *   Populates cm->expA[ECM_LC|ECM_LI|ECM_GC|ECM_GI] in place.
 *   Returns eslOK on success, eslFAIL on prediction failure.
 *
 *   Phase 1 stub: returns eslFAIL with no side effects. Real
 *   implementation lands in Phase 4.
 */
extern int  cm_FastCalibrate(CM_t *cm);

/* cm_FastCalibrateCleanup()
 *   Free model memory. Idempotent.
 */
extern void cm_FastCalibrateCleanup(void);

/* cm_FastCalibrate_PrintModels()
 *   Phase 1 debug helper: dump loaded ridge structures to fp in a
 *   format that can be diff'd against the source JSONs.
 *   Triggers lazy model load on first call.
 *   Returns eslOK after successful dump, eslFAIL on load error.
 */
extern int  cm_FastCalibrate_PrintModels(FILE *fp);

/* Phase 2 feature extraction.
 * fills feats[] with all features computable without topo_fraglen_v2.
 * Caller passes a buffer >= FAST_CAL_NFEAT_PHASE2 doubles; on return, feats[]
 * contains values indexed by the FAST_CAL_FEAT_* enum.
 *
 * Returns eslOK on success, eslFAIL if features can't be computed
 * (e.g., CM has clen <= 1).
 */
extern int cm_FastCalibrate_ExtractFeatures(CM_t *cm, double *feats);

/* Feature index enumeration. Phase 2 covers indices 0..19 (clen + 4
 * noss-fraglen + 7 STR + 8 C2). Phase 3 will extend to indices 20..26
 * (topo-noend-basic + C1_OLD).
 */
enum {
    FAST_CAL_FEAT_clen = 0,
    FAST_CAL_FEAT_mean_L_noss,
    FAST_CAL_FEAT_var_L_noss,
    FAST_CAL_FEAT_KL_noss_to_unif,
    FAST_CAL_FEAT_p_full_length,
    FAST_CAL_FEAT_n_matp,
    FAST_CAL_FEAT_pct_matp,
    FAST_CAL_FEAT_bp_density,
    FAST_CAL_FEAT_mean_matp_relent,
    FAST_CAL_FEAT_max_matp_relent,
    FAST_CAL_FEAT_sum_matp_relent,
    FAST_CAL_FEAT_mean_ml_relent,
    FAST_CAL_FEAT_mean_node_mean_g,
    FAST_CAL_FEAT_mean_node_var_g,
    FAST_CAL_FEAT_ES_full_g,
    FAST_CAL_FEAT_VarS_full_g,
    FAST_CAL_FEAT_mean_ES_g,
    FAST_CAL_FEAT_var_ES_g,
    FAST_CAL_FEAT_mean_VarS_g,
    FAST_CAL_FEAT_cov_L_ES_g,
    FAST_CAL_NFEAT_PHASE2  /* sentinel = 20 */
};

/* Feature name lookup — used by ridge dispatch to map JSON feature
 * names to indices.
 */
extern const char *cm_FastCalibrate_FeatureName(int idx);
extern int         cm_FastCalibrate_FeatureIndex(const char *name);

#endif /* cm_FAST_CALIBRATE_INCLUDED */

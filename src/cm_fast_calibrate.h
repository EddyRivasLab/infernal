/* cm_fast_calibrate.h
 * Fast CM calibration via ridge regression on compiled-in JSON models.
 *
 * Phase 1 stub: cm_FastCalibrate() returns eslFAIL with no side effects.
 * JSON parsing and model loading are implemented in Phase 1; real prediction
 * lands in Phase 4.
 *
 * Phase 2: feature extraction (cm_FastCalibrate_ExtractFeatures and friends).
 * Phase 3: topo_fraglen_v2 (noend basic) and C1_OLD legacy features.
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

/* cm_LocalMu()
 *   Fixed-lambda mini-simulation for local mu_extrap estimation.
 *   Overwrites cm->expA[EXP_CM_LC/LI]->mu_extrap and ->mu_orig.
 *   Glocal modes (EXP_CM_GC/GI) are unchanged.
 *   Called internally by cm_FastCalibrate() when g_localmu_on == 1.
 */
extern int  cm_LocalMu(CM_t *cm, ESL_RANDOMNESS *rng, int N, int use_wcap, char *errbuf);

/* Globals controlling cm_LocalMu() — set by cmbuild option parsing.
 * Defaults: N=200, seed=42, wcap=1 (on), on=1.
 * cmbuild flags: --localmu-N, --localmu-seed, --localmu-nowcap, --no-localmu
 */
extern int    g_localmu_N;
extern int    g_localmu_seed;
extern int    g_localmu_wcap;
extern int    g_localmu_on;
extern double g_localmu_lambda_lc;   /* if >0, override regression lambda for EXP_CM_LC */
extern double g_localmu_lambda_li;   /* if >0, override regression lambda for EXP_CM_LI */
extern int    g_localmu_L;           /* if >0, override per-seq L (default 2*W_eff) */
extern double g_localmu_beta;        /* QDB beta for cm_LocalMu's clone */
extern char  *g_localmu_score_dump;  /* if non-NULL, dump all hit scores to this TSV */
extern int    g_localmu_K_from_sim;  /* if 1, replace nrandhits with sim-derived K (v14 expt A) */

/* Brief 23 small-CM lambda controls.
 * cmbuild flags: --smallcm-lambda, --smallcm-clenmax, --localmu-smallonly, --localmu-fitlambda
 */
extern double g_smallcm_lambda;      /* if >0, override local-mode ridge lambda for clen < g_smallcm_clen_max */
extern int    g_smallcm_clen_max;    /* clen threshold for small-CM lambda override + --localmu-smallonly */
extern int    g_localmu_smallonly;   /* if 1, run cm_LocalMu only for clen < g_smallcm_clen_max */
extern int    g_localmu_fitlambda;   /* if 1, refit lambda in cm_LocalMu (default); 0 = hold ridge lambda */

/* cm_FastCalibrateCleanup()
 *   Free model memory. Idempotent.
 */
extern void cm_FastCalibrateCleanup(void);

/* cm_fastcal_gc_emit()
 *   brief 67: corrected consensus-emission GC fraction of a (built) CM —
 *   one consensus emission per node (MATP->MP marginalized, MATL->ML,
 *   MATR->MR), excluding MATP-internal fallback states. Used as the
 *   Option-E AT-rich routing gate; also intended for reuse by the planned
 *   null3 store-both work. If opt_ncons != NULL it receives the consensus-
 *   position count, which must equal cm->clen. Returns gc in [0,1], or
 *   -1.0 if the CM has no consensus emission mass.
 */
extern double cm_fastcal_gc_emit(const CM_t *cm, int *opt_ncons);

/* cm_FastCalibrate_PrintModels()
 *   Phase 1 debug helper: dump loaded ridge structures to fp in a
 *   format that can be diff'd against the source JSONs.
 *   Triggers lazy model load on first call.
 *   Returns eslOK after successful dump, eslFAIL on load error.
 */
extern int  cm_FastCalibrate_PrintModels(FILE *fp);

/* cm_FastCalibrate_ExtractFeatures()
 * Fills feats[0..FAST_CAL_NFEAT-1] from cm (all 27 features, Phases 2+3).
 * Caller passes a buffer >= FAST_CAL_NFEAT doubles.
 *
 * Returns eslOK on success, eslFAIL/eslEMEM on error.
 */
extern int cm_FastCalibrate_ExtractFeatures(CM_t *cm, double *feats);

/* Feature index enumeration.
 * Phase 2: indices 0..19 (clen + 4 noss-fraglen + 7 STR + 8 C2).
 * Phase 3: indices 20..26 (topo-noend-basic + C1_OLD legacy).
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
    FAST_CAL_FEAT_cov_L_ES_g,            /* index 19 */

    /* Phase 3 — topo_fraglen_v2 (noend basic) */
    FAST_CAL_FEAT_noend_mean_L,          /* 20 */
    FAST_CAL_FEAT_noend_var_L,
    FAST_CAL_FEAT_noend_KL_to_unif,
    FAST_CAL_FEAT_noend_p_full_length,   /* 23 */

    /* Phase 3 — C1_OLD legacy */
    FAST_CAL_FEAT_mean_L_str,            /* 24 */
    FAST_CAL_FEAT_var_L_str,
    FAST_CAL_FEAT_KL_str_to_unif,        /* 26 */

    /* Phase 4 — K-ridge derived feature.
     * log_clen is a simple derived value (log of clen) used as the sole
     * feature in the v4.x-converted K-ridge JSONs (clen power law for
     * glocal modes). Added as a distinct slot so ridge_predict() can
     * look it up by name like any other feature.
     */
    FAST_CAL_FEAT_log_clen,              /* 27 */

    /* Phase 5 — v5.5 feature set widening.
     * Family A: bulk IC features (8 features)
     */
    FAST_CAL_FEAT_ic_mean,              /* 28 */
    FAST_CAL_FEAT_ic_var,               /* 29 */
    FAST_CAL_FEAT_ic_p10,               /* 30 */
    FAST_CAL_FEAT_ic_p50,               /* 31 */
    FAST_CAL_FEAT_ic_p90,               /* 32 */
    FAST_CAL_FEAT_ic_skew,              /* 33 */
    FAST_CAL_FEAT_ic_mean_singlet,      /* 34 */
    FAST_CAL_FEAT_ic_mean_pair,         /* 35 */

    /* Family E: state-type ratio features (2 features) */
    FAST_CAL_FEAT_n_pair_frac,          /* 36 */
    FAST_CAL_FEAT_n_singlet_frac,       /* 37 */

    /* Family B: spatial IC features (5 features) */
    FAST_CAL_FEAT_ic_spatial_entropy,      /* 38 */
    FAST_CAL_FEAT_ic_runs_above_median,    /* 39 */
    FAST_CAL_FEAT_ic_autocorr_lag1,        /* 40 */
    FAST_CAL_FEAT_ic_autocorr_lag5,        /* 41 */
    FAST_CAL_FEAT_max_consecutive_low_ic,  /* 42 */

    /* Family C: fragment-score-weighted features (5 features) */
    FAST_CAL_FEAT_frag_score_mean,         /* 43 */
    FAST_CAL_FEAT_frag_score_var,          /* 44 */
    FAST_CAL_FEAT_frag_score_per_pos_mean, /* 45 */
    FAST_CAL_FEAT_frag_score_p90,          /* 46 */
    FAST_CAL_FEAT_cov_S_L,                 /* 47 */

    /* Family D: withend_rich topology features (11 features) */
    FAST_CAL_FEAT_withend_rich_mean_L,        /* 48 */
    FAST_CAL_FEAT_withend_rich_var_L,         /* 49 */
    FAST_CAL_FEAT_withend_rich_KL_to_unif,    /* 50 */
    FAST_CAL_FEAT_withend_rich_p_full_length, /* 51 */
    FAST_CAL_FEAT_withend_rich_skew_L,        /* 52 */
    FAST_CAL_FEAT_withend_rich_kurt_L,        /* 53 */
    FAST_CAL_FEAT_withend_rich_P10_L,         /* 54 */
    FAST_CAL_FEAT_withend_rich_P25_L,         /* 55 */
    FAST_CAL_FEAT_withend_rich_P50_L,         /* 56 */
    FAST_CAL_FEAT_withend_rich_P75_L,         /* 57 */
    FAST_CAL_FEAT_withend_rich_P90_L,         /* 58 */

    /* Group D: composition-aware features (10 features; tiny bucket only in v5.5) */
    FAST_CAL_FEAT_ic_real_mean,             /* 59 */
    FAST_CAL_FEAT_ic_real_var,              /* 60 */
    FAST_CAL_FEAT_ic_real_p10,             /* 61 */
    FAST_CAL_FEAT_ic_real_p50,             /* 62 */
    FAST_CAL_FEAT_ic_real_p90,             /* 63 */
    FAST_CAL_FEAT_ic_real_skew,            /* 64 */
    FAST_CAL_FEAT_KL_cm_uniform,           /* 65 */
    FAST_CAL_FEAT_KL_cm_genomic,           /* 66 */
    FAST_CAL_FEAT_expected_null3_lw,       /* 67 */
    FAST_CAL_FEAT_expected_null3_frag_var, /* 68 */

    /* Brief 46: NOSS hybrid predictor needs effective-sequence-number.
     * Used by base.mu_orig.ECMGC, base.K.ECMGI, and BOTH
     * largehuge_override.K cells. Appended to keep all existing indices
     * stable. Sourced from cm->eff_nseq (float32 EFFN header value). */
    FAST_CAL_FEAT_effn,                    /* 69 */

    FAST_CAL_NFEAT = 70,                 /* sentinel; total feature count */
    FAST_CAL_NFEAT_PHASE2 = 20          /* Phase 2 backward-compat alias */
};

/* Feature name lookup — used by ridge dispatch to map JSON feature
 * names to indices.
 */
extern const char *cm_FastCalibrate_FeatureName(int idx);
extern int         cm_FastCalibrate_FeatureIndex(const char *name);

#endif /* cm_FAST_CALIBRATE_INCLUDED */

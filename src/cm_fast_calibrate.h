/* cm_fast_calibrate.h
 * Fast CM calibration via ridge regression on compiled-in JSON models.
 *
 * Phase 1 stub: cm_FastCalibrate() returns eslFAIL with no side effects.
 * JSON parsing and model loading are implemented in Phase 1; real prediction
 * lands in Phase 4.
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

#endif /* cm_FAST_CALIBRATE_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmband.h
 * Eric Nawrocki 12.16.05 
 * (many functions older than 12.16.05)
 * Functions to support either CM Plan 9 (CP9) or HMMER 2.x plan 7
 * (P7) HMMs for band calculation. Includes functions for 
 * calc'ing posteriors. For Plan 7 posterior functions, see
 * HMMER's postprob.c.
 * 
 * Naming conventions:
 * CP9_* : CM Plan 9 functions
 * P7_*  : Plan 7 functions
 * 
 * Note: many functions have both CP9 and P7 versions and differ
 * only very slightly to deal with diffs b/t the two architecture
 * types (ex. P7 HMMs have no D_1 or D_M states, while CP9's do).
 *
 * Any function without "CP9" or "P7" at the beginning of the
 * function name can be used for either HMM type.
 */

#include "config.h"
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */
#include "structs.h"		/* data structures, macros, #define's   */
#include "hmmer_funcs.h"
#include "hmmer_structs.h"
#include "sre_stack.h"
#include "cplan9.h"

extern double dbl_Score2Prob(int sc, float null);

/* Functions for getting posterior probabilities from CP9 HMMs 
 * based on Ian Holmes' hmmer/src/postprob.c functions 
 * P7Forward() is in HMMER's core_algorithms.c 
 * and P7Backward() is in HMMER's postprob.c*/
extern float CP9Forward(unsigned char *dsq, int i0, int j0, struct cplan9_s *hmm, 
			struct cp9_dpmatrix_s **ret_mx);
extern float CP9Viterbi(unsigned char *dsq, int i0, int j0, struct cplan9_s *hmm, struct cp9_dpmatrix_s *mx);
extern float CP9Backward(unsigned char *dsq, int i0, int j0, struct cplan9_s *hmm, struct cp9_dpmatrix_s **ret_mx);
extern void  CP9FullPosterior(unsigned char *dsq, int i0, int j0,
			      struct cplan9_s *hmm,
			      struct cp9_dpmatrix_s *fmx,
			      struct cp9_dpmatrix_s *bmx,
			      struct cp9_dpmatrix_s *mx);
extern void CP9_ifill_post_sums(struct cp9_dpmatrix_s *post, int i0, int j0, int M,
				int *isum_pn_m, int *isum_pn_i, int *isum_pn_d);

/* Functions to determine HMM bands */
extern void CP9_hmm_band_bounds(int **post, int i0, int j0, int M, int *isum_pn, int *pn_min, int *pn_max, double p_thresh, 
				int state_type, int debug_level);

/* Functions to go from HMM bands to i and j bands on a CM */
extern void hmm2ij_bands(CM_t *cm, CP9Map_t *cp9map, int i0, int j0, int *pn_min_m, 
			 int *pn_max_m, int *pn_min_i, int *pn_max_i, int *pn_min_d, 
			 int *pn_max_d, int *imin, int *imax, int *jmin, int *jmax, 
			 int debug_level);

/* Helper functions for *_hmm2ij_bands() */
extern void hmm2ij_prestate_step0_initialize(int n, int *nss_max_imin, int *nss_min_jmax, int i0, int j0);
extern void hmm2ij_prestate_step1_set_node_inserts(int n, int *nis_imin, int *nis_imax, 
						   int *nis_jmin, int *nis_jmax,
						   int *nss_imin, int *nss_imax, 
						   int *nss_jmin, int *nss_jmax,
						   int *pn_min_i, int *pn_max_i, 
						   CP9Map_t *cp9map);
extern void hmm2ij_prestate_step2_determine_safe(int n, 	
						 int nss_max_imin_np1, int nss_min_jmax_np1,
						 int nis_imin_n, 
						 int nis_jmax_n,
						 int *safe_imax, int *safe_jmin);
extern void hmm2ij_prestate_step3_preset_node_splits(int n, int *nis_imin, int *nis_imax, 
						     int *nis_jmin, int *nis_jmax,
						     int *nss_imin, int *nss_imax, 
						     int *nss_jmin, int *nss_jmax,
						     int *pn_min_m, int *pn_max_m, 
						     int *pn_min_d, int *pn_max_d, 
						     CP9Map_t *cp9map);
extern void hmm2ij_split_state_step1_set_state_bands(int v, int n, 
						     int tmp_imin, int tmp_imax, 
						     int tmp_jmin, int tmp_jmax,
						     int *imin, int *imax, int *jmin, int *jmax,
						     int *nss_imin, int *nss_imax,
						     int *nss_jmin, int *nss_jmax);
extern void hmm2ij_insert_state_step1_set_state_bands(int v, 
						      int tmp_imin, int tmp_imax, 
						      int tmp_jmin, int tmp_jmax,
						      int *imin, int *imax, int *jmin, int *jmax);
extern void hmm2ij_state_step2_enforce_safe_trans(CM_t *cm, int v, int n, int *imax, int *jmin,
						  int *nss_imax, int *nss_jmin, 
						  int safe_imax, int safe_jmin);
extern void hmm2ij_state_step3_enforce_state_delta(CM_t *cm, int v, int *jmin, int *jmax);
extern void hmm2ij_state_step4_update_safe_holders(int v, int n, int imin_v, int jmax_v, int *nss_max_imin, 
						   int *nss_min_jmax);
extern void hmm2ij_state_step5_non_emitter_d0_hack(int v, int imax_v, int *jmin);

/* Debugging print functions */
extern void print_hmm_bands(FILE *ofp, int L, int M, int *pn_min_m, int *pn_max_m,
			    int *pn_min_i, int *pn_max_i, int *pn_min_d,
			    int *pn_max_d, double hmm_bandp, int debug_level);
extern void ij_banded_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
				      int *jmin, int *jmax, int debug_level);
extern void ijd_banded_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
				      int *jmin, int *jmax, int **hdmin, int **hdmax, 
				      int debug_level);
extern void debug_check_CP9_FB(struct cp9_dpmatrix_s *fmx, 
			       struct cp9_dpmatrix_s *bmx, 
			       struct cplan9_s *hmm, float sc, int i0, int j0,
			       unsigned char *dsq);

/* Other misc. functions */
extern void relax_root_bands(int *imin, int *imax, int *jmin, int *jmax);

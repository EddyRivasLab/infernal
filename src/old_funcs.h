#include "esl_config.h"
#include "config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sqio.h"
#include "esl_msa.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dirichlet.h"
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_fileparser.h"
#include "esl_gamma.h"
#include "esl_getopts.h"
#include "esl_gev.h"
#include "esl_gumbel.h"
#include "esl_histogram.h"
#include "esl_hyperexp.h"
#include "esl_keyhash.h"
#include "esl_minimizer.h"
#include "esl_mixgev.h"
/*#include "esl_mpi.h"*/
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msaweight.h"
#include "esl_normal.h"
#include "esl_paml.h"
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_regexp.h"
#include "esl_rootfinder.h"
#include "esl_scorematrix.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_stopwatch.h"
#include "esl_stretchexp.h"
/*#include "esl_swat.h"*/
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_weibull.h"
#include "esl_wuss.h"

#include "funcs.h"
#include "structs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

/* from old_cm_dpsearch.c */
extern float CYKScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, 
		      float cutoff, search_results_t *results);
extern float CYKBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
			   float cutoff, search_results_t *results);

extern float CYKScanRequires(CM_t *cm, int L, int W);

extern float InsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, 
			 float cutoff, search_results_t *results);
extern float InsideBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
			       float cutoff, search_results_t *results);
extern void  InsideBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax,
				 int i0, int j0, int W, 
				 int *ret_nhits, int **ret_hitr, 
				 int **ret_hiti, int **ret_hitj, float **ret_hitsc,
				 float min_thresh);
extern float iInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, 
			 float cutoff, search_results_t *results);
extern float iInsideBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
			       float cutoff, search_results_t *results);

/* from old_cm_dpalign.c */
extern int OldActuallyAlignTargets(CM_t *cm, seqs_to_aln_t *seqs_to_aln, ESL_DSQ *dsq, search_results_t *results, 
				   int first_result, int bdump_level, int debug_level, int silent_mode, ESL_RANDOMNESS *r);
extern float FInside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		     float ***alpha, float ****ret_alpha, 
		     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		     int allow_begin);
extern float IInside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		     int ***alpha, int ****ret_alpha, 
		     struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
		     int allow_begin);
extern float FOutside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		      float ***beta, float ****ret_beta, 
		      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		      int allow_begin, float ***alpha, float ****ret_alpha, int do_check);
extern float IOutside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		      int ***beta, int ****ret_beta, 
		      struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
		      int allow_begin, int ***alpha, int ****ret_alpha, int do_check);
extern void  FCMPosterior(int L, CM_t *cm, float ***alpha, float ****ret_alpha, float ***beta, 
			  float ****ret_beta, float ***post, float ****ret_post);
extern void  ICMPosterior(int L, CM_t *cm, int ***alpha, int ****ret_alpha, int ***beta, 
			 int ****ret_beta, int ***post, int ****ret_post);
extern void ICMPostalCode(CM_t *cm, int L, int ***post, Parsetree_t *tr, char **ret_pcode1, char **ret_pcode2);
extern int  Iscore2postcode(int sc);
extern void ICMCheckPosterior(int L, CM_t *cm, int ***post);
extern float FInside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
			     float ***alpha, float ****ret_alpha, 
			     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			     int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern float IInside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
			     int ***alpha, int ****ret_alpha, 
			     struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
			     int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern float FOutside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
			      float ***beta, float ****ret_beta, 
			      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			      int allow_begin, float ***alpha, float ****ret_alpha, 
			      int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern float IOutside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
			      int ***beta, int ****ret_beta, 
			      struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
			      int allow_begin, int ***alpha, int ****ret_alpha, 
			      int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void  CMPosterior_b_jd_me(int L, CM_t *cm, float ***alpha, float ****ret_alpha, 
				 float ***beta, float ****ret_beta, float ***post, float ****ret_post,
				 int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void ICMPosterior_b_jd_me(int L, CM_t *cm, int ***alpha, int ****ret_alpha, 
				 int ***beta, int ****ret_beta, int ***post, int ****ret_post,
				 int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void CMPostalCode_b_jd_me(CM_t *cm, int L, float ***post, Parsetree_t *tr,
				   int *jmin, int *jmax, int **hdmin, int **hdmax, char **ret_pcode1, char **ret_pcode2);
extern void ICMPostalCode_b_jd_me(CM_t *cm, int L, int ***post, Parsetree_t *tr,
				  int *jmin, int *jmax, int **hdmin, int **hdmax, char **ret_pcode1, char **ret_pcode2);
extern float ParsetreeSampleFromIInside(ESL_RANDOMNESS *r, CM_t *cm, ESL_DSQ *dsq, int L, int ***alpha, Parsetree_t **ret_tr,
					int ****ret_alpha);
extern float ParsetreeSampleFromIInside_b_jd_me(ESL_RANDOMNESS *r, CM_t *cm, ESL_DSQ *dsq, int L, int ***alpha, CP9Bands_t *cp9b, 
						Parsetree_t **ret_tr, int ****ret_alpha);
extern float CYKInside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, 
			    Parsetree_t **ret_tr, int *jmin, int *jmax, 
			    int **hdmin, int **hdmax, int *dmin, int *dmax);
extern void debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
					int **hdmin, int **hdmax);
extern float ** alloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
extern float CYKBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
			      int j0, int W, float cutoff, search_results_t *results);
extern float iInsideBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
				  int j0, int W, float cutoff, search_results_t *results);
extern int   ** alloc_jdbanded_vjd_kshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
extern char  ** alloc_jdbanded_vjd_yshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
/* memory management routines analogous to those in smallcyk.c for
 * handling scaled int log odds scores instead of floats. */
extern Ideckpool_t *Ideckpool_create(void);
extern void    Ideckpool_push(struct Ideckpool_s *dpool, int **deck);
extern int     Ideckpool_pop(struct Ideckpool_s *d, int ***ret_deck);
extern void    Ideckpool_free(struct Ideckpool_s *d);
extern int   **Ialloc_vjd_deck(int L, int i, int j);
extern int     Isize_vjd_deck(int L, int i, int j);
extern void    Ifree_vjd_deck(int **a, int i, int j);
extern void    Ifree_vjd_matrix(int ***a, int M, int i, int j);
extern int ** Ialloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);


/* from old_cp9_dp.c */
extern float CP9Viterbi(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
			int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
			int be_efficient, CP9_MX **ret_mx, CP9trace_t **ret_tr);
extern float CP9Forward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_isc, 
			int *ret_maxres, search_results_t *results, int do_scan, int doing_align, 
			int be_efficient, CP9_MX **ret_mx);
extern float CP9Backward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_isc, 
			 int *ret_maxres, search_results_t *results, int do_scan, int doing_align, 
			 int be_efficient, CP9_MX **ret_mx);
extern float CP9ForwardScanDemands(CP9_t *cp9, int L);
extern float CP9ViterbiAlign (ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx, CP9trace_t **ret_tr);
extern float CP9ForwardAlign (ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx);
extern float CP9BackwardAlign(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx);

/* Functions below from old_miscfuncs.c are not compiled, they're not used for 
 * any currently working code, bu are kept solely for reference */
#if 0
/* from old_miscfuncs.c */
extern float CP9Scan_dispatch(CM_t *cm,ESL_DSQ *dsq, int i0, int j0, int W, float cm_cutoff, 
			      float cp9_cutoff, search_results_t *results, int doing_cp9_stats, int *ret_flen);
extern float RescanFilterSurvivors(CM_t *cm, ESL_DSQ *dsq, search_results_t *hmm_results, int i0, 
				   int j0, int W, int padmode, int ipad, int jpad, int do_collapse,
				   float cm_cutoff, float cp9_cutoff, search_results_t *results, 
				   int *ret_flen);
extern void  CP9ScanPosterior(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *mx)
extern float FindCP9FilterThreshold(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
				    float Fmin, float Smin, float Starget, float Spad, int N, 
				    int use_cm_cutoff, float cm_ecutoff, int db_size, 
				    int emit_mode, int fthr_mode, int hmm_gum_mode, 
				    int do_fastfil, int do_Fstep, int my_rank, int nproc, 
				    int do_mpi, char *histfile, FILE *Rpts_fp, float *ret_F);
extern int   ConfigForGumbelMode(CM_t *cm, int statmode);
#endif

/* Old functions from cplan9.c for setting CP9 attributes.
 * These attributes are never used anymore, as of 01.04.07
 */
#if 0
extern void   CPlan9SetName(CP9_t *hmm, char *name);
extern void   CPlan9SetAccession(CP9_t *hmm, char *acc);
extern void   CPlan9SetDescription(CP9_t *hmm, char *desc);
extern void   CPlan9ComlogAppend(CP9_t *hmm, int argc, char **argv);
extern void   CPlan9SetCtime(CP9_t *hmm);
#endif

/* Old MPI functions from mpisupport.c */
#if 0
extern void mpi_worker_search_target(CM_t *cm, int my_rank);
extern void mpi_worker_cm_and_cp9_search(CM_t *cm, int do_fast, int my_rank);
extern void mpi_worker_cm_and_cp9_search_maxsc(CM_t *cm, int do_fast, int do_minmax, int my_rank);
extern int dsq_MPISend(ESL_DSQ *dsq, int L, int dest);
extern int dsq_MPIRecv(ESL_DSQ **ret_dsq, int *ret_L);
extern int dsq_maxsc_MPISend(char *dsq, int L, float maxsc, int dest);
extern int dsq_maxsc_MPIRecv(char **ret_dsq, int *ret_L, float *ret_maxsc);
#endif

#if 0
/* old hmmband.c funcs */
extern void simple_cp9_HMM2ijBands(CM_t *cm, int ncc, int *node_cc_left, int *node_cc_right, 
				   int *pn_min_m, int *pn_max_m, int *pn_min_i, int *pn_max_i, 
				   int *pn_min_d, int *pn_max_d, int *imin, int *imax, 
				   int *jmin, int *jmax, int **cs2hn_map, int **cs2hs_map, 
				   int debug_level);
/* old hmmband.c Plan 7 HMM banding functions */
extern void P7_hmm_band_bounds(int **post, int L, int M, int *isum_pn, int *pn_min, int *pn_max, 
			       double p_thresh, int state_type, int debug_level);
extern void P7_map_cm2hmm_and_hmm2cm(CM_t *cm, struct plan7_s *hmm, int *node_cc_left, int *node_cc_right, int *cc_node_map, int ***ret_cs2hn_map, int ***ret_cs2hs_map, int ****ret_hns2cs_map, int debug_level);
extern void P7_last_hmm_insert_state_hack(int M,  int *pn_min_m, int *pn_max_m, int *pn_min_i, int *pn_max_i);
extern void P7_last_and_first_hmm_delete_state_hack(int M,  int *pn_min_m, int *pn_max_m, int *pn_min_d, int *pn_max_d, int L);
extern void P7FullPosterior(int L, struct plan7_s *hmm, struct dpmatrix_s *forward, struct dpmatrix_s *backward, struct dpmatrix_s *mx);
extern void P7_ifill_post_sums(struct dpmatrix_s *post, int L, int M, int *isum_pn_m, int *isum_pn_i, int *isum_pn_d);
extern void P7_debug_print_post_decode(int L, int M, struct dpmatrix_s *posterior);
extern void P7_debug_print_dp_matrix(int L, int M, struct dpmatrix_s *mx);
#endif

/* old sub_cm.c functions */
#if 0 
extern int  check_sub_cm_by_sampling2(CM_t *orig_cm, CM_t *sub_cm, ESL_RANDOMNESS *r, int spos, int epos, int nseq);
#endif

#if 0
/* old cmcalibrate.c functiions */
extern seqs_to_aln_t *cm_emit_seqs_to_aln_above_cutoff(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int nseq);
extern float  search_target_cm_calibration(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float **ret_vsc);
extern int process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			    int emit_from_cm, float ***ret_vscAA, float **ret_cp9scA, float **ret_other_cp9scA);
#endif

#if 0
/* EPN, Fri Jan  4 14:27:42 2008, realized we shouldn't allow zeroing of inserts outside cmbuild, b/c this would invalidate E-values */
extern void  CMHackInsertScores(CM_t *cm);
extern void  CP9HackInsertScores(CP9_t *cp9);
#endif

#if 0
/* EPN, Mon Jan  7 14:10:43 2008, changed entropy weighting functions to HMMER3 entropy weighting functions */
extern double CM_Eweight(CM_t *cm,  const Prior_t *pri, 
			 float numb_seqs, float targetent);
extern void   ModelContent(float *ent1, float *ent2, int M);
extern double CM_Eweight_RE(CM_t *cm, const Prior_t *pri, float numb_seqs, 
			    float target_relent, float *randomseq);
extern double DRelEntropy(double *p, double *f, int n);
#endif
 

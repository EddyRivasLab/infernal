/* hmmband.c
 * EPN 12.16.05 (many functions older than this)
 * 
 * Functions to support CM Plan 9 (CP9) HMMs for band 
 * calculation. Includes functions for calc'ing 
 * posteriors. 
 *
 * Note: At the end of the function are old versions
 *       of plan 7 HMM functions related to the CP9 
 *       functions. Plan 7 HMMs for banding are no
 *       longer supported as of 10.26.06; these
 *       functions are kept here for reference.
 * 
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "easel.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* helper functions for cp9_HMM2ijBands() */
static void hmm2ij_prestate_step0_initialize(int n, int *nss_max_imin, int *nss_min_jmax, int i0, int j0);
static void hmm2ij_prestate_step1_set_node_inserts(int n, int *nis_imin, int *nis_imax, 
						   int *nis_jmin, int *nis_jmax,
						   int *nss_imin, int *nss_imax, 
						   int *nss_jmin, int *nss_jmax,
						   int *pn_min_i, int *pn_max_i, 
						   CP9Map_t *cp9map);
static void hmm2ij_prestate_step2_determine_safe(int n, 	
						 int nss_max_imin_np1, int nss_min_jmax_np1,
						 int nis_imin_n, 
						 int nis_jmax_n,
						 int *safe_imax, int *safe_jmin);
static void hmm2ij_prestate_step3_preset_node_splits(int n, int *nis_imin, int *nis_imax, 
						     int *nis_jmin, int *nis_jmax,
						     int *nss_imin, int *nss_imax, 
						     int *nss_jmin, int *nss_jmax,
						     int *pn_min_m, int *pn_max_m, 
						     int *pn_min_d, int *pn_max_d, 
						     CP9Map_t *cp9map);
static void hmm2ij_split_state_step1_set_state_bands(int v, int n, 
						     int tmp_imin, int tmp_imax, 
						     int tmp_jmin, int tmp_jmax,
						     int *imin, int *imax, int *jmin, int *jmax,
						     int *nss_imin, int *nss_imax,
						     int *nss_jmin, int *nss_jmax);
static void hmm2ij_insert_state_step1_set_state_bands(int v, 
						      int tmp_imin, int tmp_imax, 
						      int tmp_jmin, int tmp_jmax,
						      int *imin, int *imax, int *jmin, int *jmax);
static void hmm2ij_state_step2_enforce_safe_trans(CM_t *cm, int v, int n, int *imax, int *jmin,
						  int *nss_imax, int *nss_jmin, 
						  int safe_imax, int safe_jmin);
static void hmm2ij_state_step3_enforce_state_delta(CM_t *cm, int v, int *jmin, int *jmax);
static void hmm2ij_state_step4_update_safe_holders(int v, int n, int imin_v, int jmax_v, int *nss_max_imin, 
						   int *nss_min_jmax);
static void hmm2ij_state_step5_non_emitter_d0_hack(int v, int imax_v, int *jmin);


/**************************************************************************
 * EPN 10.28.06
 * Function: AllocCP9Bands()
 * 
 * Purpose:  Allocate the arrays needed for creating i and j
 *           bands on a CM based on a CP9 parse. See structs.h
 *           for description of this structure.
 *
 * Args:    
 * CM_t *cm            - the CM
 * cplan9_s *hmm       - the CP9 HMM for the CM
 * Returns: (void) 
 *
 */

CP9Bands_t *
AllocCP9Bands(CM_t *cm, struct cplan9_s *hmm)
{
  int status;
  CP9Bands_t  *cp9bands;

  ESL_ALLOC(cp9bands, sizeof(CP9Bands_t));

  cp9bands->cm_M  = cm->M;
  cp9bands->hmm_M = hmm->M;
  
  ESL_ALLOC(cp9bands->pn_min_m, sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->pn_max_m, sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->pn_min_i, sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->pn_max_i, sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->pn_min_d, sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->pn_max_d, sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->isum_pn_m,sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->isum_pn_i,sizeof(int) * (cp9bands->hmm_M+1));
  ESL_ALLOC(cp9bands->isum_pn_d, sizeof(int) * (cp9bands->hmm_M+1));

  ESL_ALLOC(cp9bands->imin,       sizeof(int)   * cp9bands->cm_M);
  ESL_ALLOC(cp9bands->imax,       sizeof(int)   * cp9bands->cm_M);
  ESL_ALLOC(cp9bands->jmin,       sizeof(int)   * cp9bands->cm_M);
  ESL_ALLOC(cp9bands->jmax,       sizeof(int)   * cp9bands->cm_M);
  ESL_ALLOC(cp9bands->safe_hdmin, sizeof(int)   * cp9bands->cm_M);
  ESL_ALLOC(cp9bands->safe_hdmax, sizeof(int)   * cp9bands->cm_M);
  ESL_ALLOC(cp9bands->hdmin,      sizeof(int *) * cp9bands->cm_M);
  ESL_ALLOC(cp9bands->hdmax,      sizeof(int *) * cp9bands->cm_M);
  cp9bands->hdmin_mem = NULL;
  cp9bands->hdmax_mem = NULL;
  /* NOTE: cp9bands->hdmin and hdmax are 2D arrays, the ptrs are 
   * alloc'ed here, but the actually memory is alloc'ed by
   * hmmband.c:cp9_Seq2Bands() with a call to hmmband.c:cp9_GrowHDBands(). 
   * We set hdmin[0] = hdmax[0] = NULL so we know not to free them 
   * if they were never alloc'ed.
   */
  cp9bands->hd_needed  = 0;
  cp9bands->hd_alloced = 0;
  return cp9bands;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return NULL; /* never reached */
}

/* Function: FreeCP9Bands() 
 * Returns: (void) 
 */
void 
FreeCP9Bands(CP9Bands_t *cp9bands)
{
  free(cp9bands->imin);
  free(cp9bands->imax);
  free(cp9bands->jmin);
  free(cp9bands->jmax);
  free(cp9bands->safe_hdmin);
  free(cp9bands->safe_hdmax);
  if(cp9bands->hdmin_mem != NULL)
    free(cp9bands->hdmin_mem); /* all v were malloc'ed as a block */
  if(cp9bands->hdmax_mem != NULL) 
    free(cp9bands->hdmax_mem); /* all v were malloc'ed as a block */
  free(cp9bands->hdmin);
  free(cp9bands->hdmax);

  free(cp9bands->pn_min_m);
  free(cp9bands->pn_max_m);
  free(cp9bands->pn_min_i);
  free(cp9bands->pn_max_i);
  free(cp9bands->pn_min_d);
  free(cp9bands->pn_max_d);
  free(cp9bands->isum_pn_m);
  free(cp9bands->isum_pn_i);
  free(cp9bands->isum_pn_d);
  free(cp9bands);
}

/* Function: DScore2Prob()
 * 
 * Purpose:  Convert an integer log_2 odds score back to a probability;
 *           needs the null model probability, if any, to do the conversion.
 */
double 
DScore2Prob(int sc, float null)
{
  if (sc == -INFTY) return 0.;
  else              return (null * sreEXP2((double) sc / INTSCALE));
}


/* Function: cp9_Seq2Bands
 * Date    : EPN, Mon Jan  8 07:23:34 2007
 *           EPN, Wed Oct 17 04:53:58 2007  [updated/optimized]
 *
 * Purpose:  Given a CM with precalc'ed CP9 HMM and CP9Map, a sequence and 
 *           a CP9Bands_t structure, calculate the HMM bands and store them
 *           in the CP9Bands_t structure.
 *           
 * Args:     cm          - the covariance model
 *           errbuf      - char buffer for reporting errors
 *           dsq         - sequence in digitized form
 *           i0          - start of target subsequence (often 1, beginning of sq)
 *           j0          - end of target subsequence (often L, end of sq)
 *           cp9b        - PRE-ALLOCATED, the HMM bands for this sequence, filled here.
 *           doing_search- TRUE if we're going to use these HMM bands for search, not alignment
 *           debug_level - verbosity level for debugging printf()s
 * Return:  eslOK on success;
 * 
 */
int
cp9_Seq2Bands(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, int doing_search, int debug_level)
{
  int             status;
  int             use_sums; /* TRUE to fill and use posterior sums during HMM band calc  *
			     * leads to wider bands                                      */
  int             sc;
  float           fsc;
  CP9_dpmatrix_t *cp9_fwd;   /* growable DP matrix for forward                       */
  CP9_dpmatrix_t *cp9_bck;   /* growable DP matrix for backward                      */
  int do_scan2bands;         /* TRUE to use scanning Forward/Backward to get posteriors */
  int be_safe = TRUE;        /* TEMPORARY, pass this in after calcing it once in actually_align_targets() */

  /* Contract checks */
  if(cm->cp9 == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_seq2bands, but cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_seq2bands, but cm->cp9map is NULL.\n");
  if(dsq == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_seq2bands, dsq is NULL.");
  if((cm->align_opts & CM_ALIGN_HBANDED) && (cm->search_opts & CM_SEARCH_HBANDED))           ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_seq2bands, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both up, exactly 1 must be up.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED)))        ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_seq2bands, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down, exactly 1 must be up.\n");
  if((cm->search_opts & CM_SEARCH_HMMSCANBANDS) && (!(cm->search_opts & CM_SEARCH_HBANDED))) ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_seq2bands, CM_SEARCH_HMMSCANBANDS flag raised, but not CM_SEARCH_HBANDED flag, this doesn't make sense\n");
  
  use_sums = ((cm->align_opts & CM_ALIGN_SUMS) || (cm->search_opts & CM_SEARCH_SUMS)) ? TRUE : FALSE;
    
  /* Step 1: Get HMM Forward/Backward DP matrices.
   * Step 2: F/B       -> HMM bands.
   * Step 3: HMM bands -> CM bands.
   */

  /* Step 1: Get HMM Forward/Backward DP matrices. */

  do_scan2bands = (cm->search_opts & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;
  /* TO DO have these guys return status codes */
  fsc = Xcp9_FastForward(cm, dsq, i0, j0, j0-i0+1, 0., NULL, NULL, NULL,
			 do_scan2bands, /* are we using scanning Forward/Backward */
			 TRUE,      /* we are going to use posteriors to align */
			 FALSE,     /* we're not rescanning */
			 FALSE,     /* don't be memory efficient */
			 be_safe,   /* can we accelerate w/ no -INFTY logsum funcs? */
			 &cp9_fwd); /* give the DP matrix back */

  fsc = Xcp9_FastBackward(cm, dsq, i0, j0, (j0-i0+1), 0, NULL, NULL, NULL,
			  do_scan2bands, /* are we using scanning Forward/Backward */
			  TRUE,  /* we are going to use posteriors to align */
			  FALSE, /* we're not rescanning */
			  FALSE, /* don't be memory efficient */
			  &cp9_bck); /* give the DP matrix back */

  if(cm->align_opts & CM_ALIGN_CHECKFB) cp9_DebugCheckFB(cp9_fwd, cp9_bck, cm->cp9, sc, i0, j0, dsq);

  /* Step 2: F/B -> HMM bands. */
  if(use_sums){
    if((status = cp9_FB2HMMBandsWithSums(cm->cp9, errbuf, dsq, cp9_fwd, cp9_bck, cp9b, i0, j0, cp9b->hmm_M,
					 (1.-cm->tau), (cm->search_opts & CM_SEARCH_HMMSCANBANDS), debug_level)) != eslOK) return status;
  }
  else {
    if((status = cp9_FB2HMMBands(cm->cp9, errbuf, dsq, cp9_fwd, cp9_bck, cp9b, i0, j0, cp9b->hmm_M,
				 (1.-cm->tau), (cm->search_opts & CM_SEARCH_HMMSCANBANDS), debug_level)) != eslOK) return status;
  }
  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, j0, cp9b, cm->tau, 1);

  /* Step 3: HMM bands  ->  CM bands. */
  cp9_HMM2ijBands(cm, errbuf, cm->cp9map, i0, j0, cp9b->pn_min_m, cp9b->pn_max_m, 
		  cp9b->pn_min_i, cp9b->pn_max_i, cp9b->pn_min_d, cp9b->pn_max_d, 
		  cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, debug_level);
  
  /* Use the CM bands on i and j to get bands on d, specific to j. */
  if(doing_search) cp9_RelaxRootBandsForSearch(cm, cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax);
  /* cp9_GrowHDBands() must be called before ij2d_bands() so hdmin, hdmax are adjusted for new seq */
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) return status; 
  ij2d_bands(cm, (j0-i0+1), cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, debug_level);
  
#if eslDEBUGLEVEL >= 1
  if((status = cp9_ValidateBands(cm, errbuf, cp9b, i0, j0)) != eslOK) return status;
#endif

  if(debug_level > 0) PrintDPCellsSaved_jd(cm, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, (j0-i0+1));

  /* free matrices, currently no option of returning them */
  FreeCPlan9Matrix(cp9_fwd);
  FreeCPlan9Matrix(cp9_bck);

  return eslOK;
}


/*
 * Function: cp9_Seq2Posteriors
 * Date    : EPN, Mon Jan  8 07:27:21 2007
 *
 * Purpose:  Given a CM with precalc'ed CP9 HMM and CP9Map, and a sequence,
 *           run HMM Forward and Backward algorithms, and return a CP9 posterior
 *           matrix.
 *           
 * Args:     cm           - the covariance model
 *           errbuf       - char buffer for error messages
 *           dsq          - sequence in digitized form
 *           i0           - start of target subsequence (often 1, beginning of dsq)
 *           j0           - end of target subsequence (often L, end of dsq)
 *           ret_cp9_post - RETURN: the HMM posterior matrix for this sequence
 *           debug_level  - verbosity level for debugging printf()s
 *           
 * Return:  eslOK on success
 */
int
cp9_Seq2Posteriors(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, CP9_dpmatrix_t **ret_cp9_post,
		   int debug_level)
{
  /*CP9_dpmatrix_t *cp9_mx;*/    /* growable DP matrix for viterbi                       */
  CP9_dpmatrix_t *cp9_fwd;       /* growable DP matrix for forward                       */
  CP9_dpmatrix_t *cp9_bck;       /* growable DP matrix for backward                      */
  CP9_dpmatrix_t *cp9_post;      /* growable DP matrix for posterior decode              */
  float sc;
  int do_scan2bands;             /* TRUE to use scanning Forward/Backward to get posteriors
				  * that we'll use for a CM scan */
  int be_safe = TRUE;        /* TEMPORARY, pass this in after calcing it once in actually_align_targets() */

  /* Contract checks */
  if(dsq == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors(), dsq is NULL.");
  if(cm->cp9 == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, but cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, but cm->cp9map is NULL.\n");
  if((cm->align_opts & CM_ALIGN_HBANDED) && (cm->search_opts & CM_SEARCH_HBANDED)) 
    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both up, exactly 1 must be up.\n");
  if((cm->align_opts & CM_SEARCH_HMMSCANBANDS) && (! (cm->search_opts & CM_SEARCH_HBANDED))) 
    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, CM_SEARCH_HMMSCANBANDS flag raised, but not CM_SEARCH_HBANDED flag, this doesn't make sense\n");

  do_scan2bands = (cm->search_opts & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;

  /* Step 1: Get HMM posteriors.*/
  sc = Xcp9_FastForward(cm, dsq, i0, j0, j0-i0+1, 0., NULL, NULL, NULL,
			do_scan2bands, /* are we using scanning Forward/Backward */
			TRUE,      /* we are going to use posteriors to align */
			FALSE,     /* we're not rescanning */
			FALSE,     /* don't be memory efficient */
			be_safe,   /* can we accelerate w/ no -INFTY logsum funcs? */
			&cp9_fwd); /* give the DP matrix back */

  if(debug_level > 0) printf("CP9 Forward  score : %.4f\n", sc);
  sc = Xcp9_FastBackward(cm, dsq, i0, j0, (j0-i0+1), 0, NULL, NULL, NULL,
			 do_scan2bands, /* are we using scanning Forward/Backward */
			 TRUE,  /* we are going to use posteriors to align */
			 FALSE, /* we're not rescanning */
			 FALSE, /* don't be memory efficient */
			 &cp9_bck); /* give the DP matrix back */
  if(debug_level > 0) printf("CP9 Backward score : %.4f\n", sc);

  if(cm->align_opts & CM_ALIGN_CHECKFB) 
    cp9_DebugCheckFB(cp9_fwd, cp9_bck, cm->cp9, sc, i0, j0, dsq);

  /* Get posteriors */
  cp9_post = cp9_bck;
  cp9_Posterior(dsq, i0, j0, cm->cp9, cp9_fwd, cp9_bck, cp9_post, do_scan2bands);

  FreeCPlan9Matrix(cp9_fwd);
  if(ret_cp9_post != NULL)
    *ret_cp9_post = cp9_post;
  else
    FreeCPlan9Matrix(cp9_post);
  return eslOK;
}


/*
 * Function: cp9_FB2HMMBands()
 * Date:     EPN, 04.03.06
 *           EPN, Mon Oct 15 18:20:42 2007 [updated/optimized]
 *
 * Purpose: Determine the band on all HMM states given a Forward and
 *          Backward matrix. Do this by calculating and summing log posterior
 *          probabilities that each state emitted/was visited at each posn,
 *          starting at the sequence ends, and creeping in, until the half the
 *          maximum allowable probability excluded is reached on each side.
 *
 * Args:
 *
 * CP9_t hmm        the HMM
 * errbuf           char buffer for error messages
 * dsq              the digitized sequence
 * CP9_dpmatrix_t fmx: forward DP matrix, already calc'ed
 * CP9_dpmatrix_t bmx: backward DP matrix, already calc'ed
 * CP9Bands_t cp9b  CP9 bands data structure
 * int i0           start of target subsequence (often 1, beginning of dsq)
 * int j0           end of target subsequence (often L, end of dsq)
 * int   M          number of nodes in HMM (num columns of post matrix)
 * double p_thresh  the probability mass we're requiring is within each band
 * int did_scan     TRUE if Forward/Backward were run in 'scan mode'
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 * 
 * Returns: eslOK on success;
 */
int
cp9_FB2HMMBands(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_dpmatrix_t *fmx, CP9_dpmatrix_t *bmx, CP9Bands_t *cp9b, 
		 int i0, int j0, int M, double p_thresh, int did_scan, int debug_level)
{
  int status;
  int k;                                  /* counter over nodes of the model */
  int L = j0-i0+1;                        /* length of sequence */
  int thresh = Prob2Score(((1. - p_thresh)/2.), 1.); /* allowable prob mass excluded on each side */
  int max;                                /* temporary max value */
  int pnmax;                              /* position that gives max */

  /* *_m = match, *_i = insert, *_d = delete */
  int *nset_m, *nset_i, *nset_d;          /* [0..k..hmm->M], has minimum been set for this state? */
  int *xset_m, *xset_i, *xset_d;          /* [0..k..hmm->M], has maximum been set for this state? */
  int *mass_m, *mass_i, *mass_d;          /* [0..k..hmm->M], summed log prob of post->mx[i][k] from 0..k or k..L */
  int i, ip;                              /* actual position and relative position in sequence, ip = i-i0+1 */
  CP9_dpmatrix_t *post;                   /* posterior dp matrix for CP9 HMM, alloc'ed, filled and free'd here */
  int sc;                                 /* summed score of all parses (derived from backward matrix) 
					   * if(cm->search_opts & CM_SEARCH_HMMSCANBANDS) Forward and Backward
					   * were run in 'scan mode' where each residue can be begin/end of a parse,
					   * so we have to sum up parses that end at each posn, 
					   * if ! (cm->search_opts & CM_SEARCH_HMMSCANBANDS) we know we have 
					   * to start at residue i0 and end at residue j0, so sc is simply bmx->mmx[0][0]
					   */
  /* printf("in cp9_FB2HMMBands()\n"); */
  post = AllocCPlan9Matrix(L+1, M, NULL, NULL, NULL, NULL, NULL);
  ESL_ALLOC(nset_m, sizeof(int) * (M+1) * 9);
  ESL_ALLOC(nset_i, sizeof(int) * (M+1));
  ESL_ALLOC(nset_d, sizeof(int) * (M+1));
  ESL_ALLOC(xset_m, sizeof(int) * (M+1));
  ESL_ALLOC(xset_i, sizeof(int) * (M+1));
  ESL_ALLOC(xset_d, sizeof(int) * (M+1));
  ESL_ALLOC(mass_m, sizeof(int) * (M+1));
  ESL_ALLOC(mass_i, sizeof(int) * (M+1));
  ESL_ALLOC(mass_d, sizeof(int) * (M+1));  

  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  esl_vec_ISet(nset_m, M+1, FALSE);
  esl_vec_ISet(nset_i, M+1, FALSE);
  esl_vec_ISet(nset_d, M+1, FALSE);
  esl_vec_ISet(xset_m, M+1, FALSE);
  esl_vec_ISet(xset_i, M+1, FALSE);
  esl_vec_ISet(xset_d, M+1, FALSE);

  if(did_scan) { /* Forward/Backward run in 'scan mode' which allow parses to start/stop anywhere */
    sc = -INFTY;
    for (ip = 0; ip <= L; ip++) {
      /*printf("bmx->mmx[i:%d][0]: %d\n", i, bmx->mmx[ip][0]); */
      sc = ILogsum(sc, (bmx->mmx[ip][0])); 
    }
  }
  else sc = bmx->mmx[0][0]; /* Forward/Backward run in 'align mode' parses must start at i0, end at j0 */
  /* sc is summed log prob of all possible parses of seq i0..j0 */

  /* comment *: off-by-one issue with non-emitters (includes all D states and M_0): 
   * pn_min_d[k] = i, means posn i was last residue emitted
   * prior to entering node k's delete state. However, for a CM,
   * if a delete states sub-parsetree is bounded by i' and j', then
   * positions i' and j' HAVE YET TO BE EMITTED.
   * For M_0, so we don't have to check each node to see if k == 0, we
   * do the off-by-one correction at the end of the function.
   */
  
  /* note boundary conditions, ip = 0, i = i0-1 */
  post->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc; /* fmx->mmx[0][0] is 0, bmx->mmx[1][0] is overall score */
  post->imx[0][0] = -INFTY; /*need seq to get here*/
  post->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  if((mass_m[0] = post->mmx[0][0]) > thresh) { 
    cp9b->pn_min_m[0] = 0; /* "off-by-one" (from comment * above) is dealt with special at end of function for M_0 */
    nset_m[0] = TRUE; 
  }
  mass_i[0] = -INFTY; /* b/c post->imx[0][0] is -INFTY, set above */
  mass_d[0] = -INFTY; /* b/c post->dmx[0][0] is -INFTY, set above */

  for (k = 1; k <= M; k++) {
    post->mmx[0][k] = -INFTY; /*need seq to get here*/
    post->imx[0][k] = -INFTY; /*need seq to get here*/
    post->dmx[0][k] = fmx->dmx[0][k] + bmx->dmx[0][k] - sc;
    /* mass_m[k] doesn't change b/c post->mmx[0][k] is -INFTY */
    /* mass_i[k] doesn't change b/c post->mmx[0][k] is -INFTY */
    if((mass_d[k] = post->dmx[0][k]) > thresh) { 
      cp9b->pn_min_d[k] = 0+1; /* see "off-by-one" comment * above */
      nset_d[k] = TRUE; 
    }
  }
  /* Find minimum position in band for each state (M,I,D) of each node (0..M) */
  for (ip = 1; ip <= L; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      post->mmx[ip][0] = -INFTY; /* M_0 doesn't emit */
      post->imx[ip][0] = ESL_MAX(fmx->imx[ip][0] + bmx->imx[ip][0] - hmm->isc[dsq[i]][0] - sc, -INFTY);
      /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
      post->dmx[ip][0] = -INFTY; /* D_0 doesn't exist */

      for(k = 1; k <= M; k++)
	{
	  post->mmx[ip][k] = ESL_MAX(fmx->mmx[ip][k] + bmx->mmx[ip][k] - hmm->msc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  post->imx[ip][k] = ESL_MAX(fmx->imx[ip][k] + bmx->imx[ip][k] - hmm->isc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  post->dmx[ip][k] = ESL_MAX(fmx->dmx[ip][k] + bmx->dmx[ip][k] - sc, -INFTY);

	  if(! nset_m[k]) { 
	    if((mass_m[k] = ILogsum(mass_m[k], post->mmx[ip][k])) > thresh) { 
	      cp9b->pn_min_m[k] = i;
	      nset_m[k] = TRUE; 
	    }
	  }
	  if(! nset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], post->imx[ip][k])) > thresh) { 
	      cp9b->pn_min_i[k] = i;
	      nset_i[k] = TRUE; 
	    }
	  }
	  if(! nset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], post->dmx[ip][k])) > thresh) { 
	      cp9b->pn_min_d[k] = i+1; /* see "off-by-one" comment * above */
	      nset_d[k] = TRUE; 
	    }
	  }
	}
    }	  

  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  /* Find maximum position in band for each state (M,I,D) of each node (0..M)
   * by moving from L down to 1 */
  for (ip = L; ip >= 1; ip--) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      for(k = 0; k <= M; k++)
	{
	  if(! xset_m[k]) { 
	    if((mass_m[k] = ILogsum(mass_m[k], post->mmx[ip][k])) > thresh) { 
	      cp9b->pn_max_m[k] = i;
	      xset_m[k] = TRUE; 
	    }
	  }
	  if(! xset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], post->imx[ip][k])) > thresh) { 
	      cp9b->pn_max_i[k] = i;
	      xset_i[k] = TRUE; 
	    }
	  }
	  if(! xset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], post->dmx[ip][k])) > thresh) { 
	      cp9b->pn_max_d[k] = i+1; /* see "off-by-one" comment * above */
	      xset_d[k] = TRUE; 
	    }
	  }
	}
    }	  
  /* note boundary conditions, ip = 0, i = i0-1 */
  if(! xset_m[0]) { 
    if((mass_m[0] = ILogsum(mass_m[0], post->mmx[0][0])) > thresh) { 
      cp9b->pn_max_m[0] = 0; /* "off-by-one" (from comment * above) is dealt with special at end of function for M_0 */
      xset_m[0] = TRUE; 
    }
  }
  /* mass_i[0] is unchaged because b/c post->imx[0][0] is -INFTY, set above */
  /* mass_d[0] is unchaged because b/c post->dmx[0][0] is -INFTY, set above */
  for (k = 1; k <= M; k++) {
    /* mass_m[k] doesn't change b/c post->mmx[0][k] is -INFTY */
    /* mass_i[k] doesn't change b/c post->mmx[0][k] is -INFTY */
    if(!xset_d[k]) { 
      if((mass_d[k] = ILogsum(mass_d[k], post->dmx[0][k])) > thresh) { 
	cp9b->pn_max_d[k] = 0+1; /* see "off-by-one" comment * above */
	xset_d[k] = TRUE; 
      }
    }
  }	 
  /* Some states may not have had their min/max set. This occurs if the entire
   * state is outside the band (i.e. the summed probablity the state is entered for ANY i
   * is less than our threshold. Current strategy in this situation is to set the
   * band to width 1 of the most likely position for that state, but to do that we
   * need to find what the most likely posn is, we could do this in the loop above,
   * but this is a rare situation, and so that turns out to be wasteful.
   */
  for(k = 0; k <= M; k++)
    {
      /* theoretically either nset_*[k] and xset_*[k] should be either both TRUE or both
       * FALSE, but I'm slightly worried about rare precision issues, so we check if one 
       * or the other is unset, and if so, we set both to argmax position */
      if((! nset_m[k]) || (! xset_m[k])) { 
	max = post->mmx[0][k];
	for(ip = 1; ip <= L; ip++)
	  if(post->mmx[ip][k] > max) { pnmax = i; max = post->mmx[ip][k]; }
	cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = pnmax;
      }
      if((! nset_i[k]) || (! xset_i[k])) { 
	max = post->imx[0][k];
	for(ip = 1; ip <= L; ip++)
	  if(post->imx[ip][k] > max) { pnmax = i; max = post->imx[ip][k]; }
	cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = pnmax;
      }
      if((! nset_d[k]) || (! xset_d[k])) { 
	max = post->dmx[0][k];
	for(ip = 1; ip <= L; ip++)
	  if(post->dmx[ip][k] > max) { pnmax = i; max = post->dmx[ip][k]; }
	cp9b->pn_min_d[k] = cp9b->pn_max_d[k] = pnmax + 1; /* see "off-by-one" comment * above */
      }
    }
  /* correct for M_0 off-by-one explained in comment * above */
  cp9b->pn_min_m[0]++;
  cp9b->pn_max_m[0]++;
  cp9b->pn_min_d[0] = -1; /* D_0 doesn't exist */
  cp9b->pn_max_d[0] = -1; /* D_0 doesn't exist */

  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, j0, cp9b, (1.-p_thresh), 1);
  free(mass_m);
  free(mass_i);
  free(mass_d);
  free(nset_m);
  free(nset_i);
  free(nset_d);
  free(xset_m);
  free(xset_i);
  free(xset_d);
  FreeCPlan9Matrix(post);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/*****************************************************************************
 * Function: cp9_FB2HMMBandsWithSums()
 * Date:     EPN, Wed Oct 17 10:22:44 2007
 *
 * Purpose: Determine the band on all HMM states given a Forward and
 *          Backward matrix. Do this by calculating and summing log posterior
 *          probabilities that each state emitted/was visited at each posn,
 *          starting at the sequence ends, and creeping in, until the half the
 *          maximum allowable probability excluded is reached on each side.
 *
 * CP9_t hmm        the HMM
 * errbuf           char buffer for error messages
 * dsq              the digitized sequence
 * CP9_dpmatrix_t fmx: forward DP matrix, already calc'ed
 * CP9_dpmatrix_t bmx: backward DP matrix, already calc'ed
 * CP9Bands_t cp9b  CP9 bands data structure
 * int i0           start of target subsequence (often 1, beginning of dsq)
 * int j0           end of target subsequence (often L, end of dsq)
 * int   M          number of nodes in HMM (num columns of post matrix)
 * double p_thresh  the probability mass we're requiring is within each band
 * int did_scan     TRUE if Forward/Backward were run in 'scan mode'
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *
 * Returns: eslOK on success;
 *
 *****************************************************************************/
int
cp9_FB2HMMBandsWithSums(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_dpmatrix_t *fmx, CP9_dpmatrix_t *bmx, CP9Bands_t *cp9b, 
			int i0, int j0, int M, double p_thresh, int did_scan, int debug_level)
{
  int status;
  int k;                                  /* counter over nodes of the model */
  int L = j0-i0+1;                        /* length of sequence */
  int thresh = Prob2Score(((1. - p_thresh)/2.), 1.); /* allowable prob mass excluded on each side */

  /* *_m = match, *_i = insert, *_d = delete */
  int *kthresh_m, *kthresh_i, *kthresh_d; /* [0..k..hmm->M], individual thresholds for each state */
  int *nset_m, *nset_i, *nset_d;          /* [0..k..hmm->M], has minimum been set for this state? */
  int *xset_m, *xset_i, *xset_d;          /* [0..k..hmm->M], has maximum been set for this state? */
  int *mass_m, *mass_i, *mass_d;          /* [0..k..hmm->M], summed log prob of post->mx[i][k] from 0..k or k..L */
  int i, ip;                              /* actual position and relative position in sequence, ip = i-i0+1 */
  CP9_dpmatrix_t *post;                   /* posterior dp matrix for CP9 HMM */
  
  /* printf("in cp9_FB2HMMBandsWithSums()\n"); */
  /* allocations and initializations */
  ESL_ALLOC(nset_m, sizeof(int) * (M+1));
  ESL_ALLOC(nset_i, sizeof(int) * (M+1));
  ESL_ALLOC(nset_d, sizeof(int) * (M+1));
  ESL_ALLOC(xset_m, sizeof(int) * (M+1));
  ESL_ALLOC(xset_i, sizeof(int) * (M+1));
  ESL_ALLOC(xset_d, sizeof(int) * (M+1));
  ESL_ALLOC(mass_m, sizeof(int) * (M+1));
  ESL_ALLOC(mass_i, sizeof(int) * (M+1));
  ESL_ALLOC(mass_d, sizeof(int) * (M+1));  
  ESL_ALLOC(kthresh_m, sizeof(int) * (M+1));
  ESL_ALLOC(kthresh_i, sizeof(int) * (M+1));
  ESL_ALLOC(kthresh_d, sizeof(int) * (M+1));  

  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  esl_vec_ISet(nset_m, M+1, FALSE);
  esl_vec_ISet(nset_i, M+1, FALSE);
  esl_vec_ISet(nset_d, M+1, FALSE);
  esl_vec_ISet(xset_m, M+1, FALSE);
  esl_vec_ISet(xset_i, M+1, FALSE);
  esl_vec_ISet(xset_d, M+1, FALSE);

  /* get the posterior matrix first, we need it b/c each state will have a different log prob threshold */
  post = bmx;
  cp9_Posterior(dsq, i0, j0, hmm, fmx, bmx, post, did_scan);

  /* fill ipost_sums in cp9bands data structure */
  cp9_IFillPostSums(post, cp9b, i0, j0);

  /* set state dependent cutoff thresholds for log prob mass we need on each side (this is unique to 
   * WithSums() function */
  for(k = 0; k <= M; k++) {
    kthresh_m[k] = thresh + cp9b->isum_pn_m[k];
    kthresh_i[k] = thresh + cp9b->isum_pn_i[k];
    kthresh_d[k] = thresh + cp9b->isum_pn_d[k];
  }
  /* comment *: off-by-one issue with non-emitters (includes all D states and M_0): 
   * pn_min_d[k] = i, means posn i was last residue emitted
   * prior to entering node k's delete state. However, for a CM,
   * if a delete states sub-parsetree is bounded by i' and j', then
   * positions i' and j' HAVE NOT YET BEEN EMITTED.
   * For M_0, so we don't have to check each node to see if k == 0, we
   * do the off-by-one correction at the end of the function.
   */

  /* Find minimum position in band for each state (M,I,D) of each node (0..M) */
  for (ip = 0; ip <= L; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      for(k = 0; k <= M; k++)
	{
	  if(! nset_m[k]) { 
	    if((mass_m[k] = ILogsum(mass_m[k], post->mmx[ip][k])) > kthresh_m[k]) { 
	      cp9b->pn_min_m[k] = i;
	      nset_m[k] = TRUE; 
	    }
	  }
	  if(! nset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], post->imx[ip][k])) > kthresh_i[k]) { 
	      cp9b->pn_min_i[k] = i;
	      nset_i[k] = TRUE; 
	    }
	  }
	  if(! nset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], post->dmx[ip][k])) > kthresh_d[k]) { 
	      cp9b->pn_min_d[k] = i+1; /* see "off-by-one" comment * above */
	      nset_d[k] = TRUE; 
	    }
	  }
	}
    }	  
  /* Find maximum position in band for each state (M,I,D) of each node (0..M)
   * by moving from L down to 0 */
  /* reset mass_* arrays */
  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  for (ip = L; ip >= 0; ip--) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      for(k = 0; k <= M; k++)
	{
	  if(! xset_m[k]) { 
	    if((mass_m[k] = ILogsum(mass_m[k], post->mmx[ip][k])) > kthresh_m[k]) { 
	      cp9b->pn_max_m[k] = i;
	      xset_m[k] = TRUE; 
	    }
	  }
	  if(! xset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], post->imx[ip][k])) > kthresh_i[k]) { 
	      cp9b->pn_max_i[k] = i;
	      xset_i[k] = TRUE; 
	    }
	  }
	  if(! xset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], post->dmx[ip][k])) > kthresh_d[k]) { 
	      cp9b->pn_max_d[k] = i+1; /* see "off-by-one" comment * above */
	      xset_d[k] = TRUE; 
	    }
	  }
	}
    }	  

  /* all states should have their min/max set because we've normalized the probability
   * of entering each state to 1.0, so we assert this to be true */
  /* assert(nset_m[0]);
     assert(nset_i[0]);
     assert(xset_m[0]);
     assert(xset_i[0]); */
  ESL_DASSERT1((nset_m[0]));
  ESL_DASSERT1((nset_i[0]));
  ESL_DASSERT1((xset_m[0]));
  ESL_DASSERT1((xset_i[0]));
  /* D_0 state does not exist */
  for(k = 1; k <= M; k++)
    {
      /*assert(nset_m[k]);
	assert(nset_i[k]);
	assert(nset_d[k]);
	assert(xset_m[k]);
	assert(xset_i[k]);
	assert(xset_d[k]);*/
      ESL_DASSERT1((nset_m[k]));
      ESL_DASSERT1((nset_i[k]));
      ESL_DASSERT1((nset_d[k]));
      ESL_DASSERT1((xset_m[k]));
      ESL_DASSERT1((xset_i[k]));
      ESL_DASSERT1((xset_d[k]));
    }
  /* correct for M_0 off-by-one explained in comment * above */
  cp9b->pn_min_m[0]++;
  cp9b->pn_max_m[0]++;
  cp9b->pn_min_d[0] = -1; /* D_0 doesn't exist */
  cp9b->pn_max_d[0] = -1; /* D_0 doesn't exist */

  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, j0, cp9b, (1.-p_thresh), 1);
  free(mass_m);
  free(mass_i);
  free(mass_d);
  free(nset_m);
  free(nset_i);
  free(nset_d);
  free(xset_m);
  free(xset_i);
  free(xset_d);
  free(kthresh_m);
  free(kthresh_i);
  free(kthresh_d);
  /* careful not to free post, it is bmx, but overwritten, and the caller
   * is responsible for free'ing it. */
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/* Functions for getting posterior probabilities from the HMMs 
 * based on Ian Holmes' hmmer/src/postprob.c functions 
 * CP9ForwardAlign()
 * CP9ViterbiAlign()
 * CP9BackwardAlign()
 * CP9Posterior()
 * CP9_ifill_post_sums()
 */

/***********************************************************************
 * Function: CP9ForwardAlign
 * based on  P7Forward() <-- this function's comments below  
 *           from HMMER 2.4devl core_algorithms.c
 *
 * Purpose:  The Forward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Note:     Differs from CP9_scan.c:CP9Forward() in that this function
 *           does not work in 'scan' mode. CP9_scan.c:CP9Forward() can
 *           do scanning OR alignment, but is somewhat more complex than
 *           this function.
 *
 * Args:     dsq    - sequence in digitized form
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
CP9ForwardAlign(ESL_DSQ *dsq, int i0, int j0, struct cplan9_s *hmm, 
	      struct cp9_dpmatrix_s **ret_mx)
{
  if(dsq == NULL) 
    cm_Fail("In CP9ForwardAlign() dsq is NULL.");
  struct cp9_dpmatrix_s *mx;
  int **mmx;
  int **imx;
  int **dmx;
  int **elmx;
  int  *erow;
  int   i,k;
  int   sc;
  int   L;		/* subsequence length */
  int   ip;		/* i': relative position in the subsequence  */

  printf("in CP9Forward()\n");
  L  = j0-i0+1;		/* the length of the subsequence */

  /* Allocate a DP matrix with 0..L rows, 0..M-1 
   */ 
  mx = AllocCPlan9Matrix(L+1, hmm->M, &mmx, &imx, &dmx, &elmx, &erow);

  /* Initialization of the zero row.
   */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= hmm->M; k++)
    dmx[0][k] = ILogsum(ILogsum(mmx[0][k-1] + hmm->tsc[CTMD][k-1],
				imx[0][k-1] + hmm->tsc[CTID][k-1]),
			dmx[0][k-1] + hmm->tsc[CTDD][k-1]);
  
  erow[0] = dmx[0][hmm->M] + hmm->tsc[CTDM][hmm->M]; 
  /*printf("i: %d score(erow[ip]): %f erow[ip]: %d\n", 0, Scorify(erow[0]), erow[0]);*/
  /* Recursion. Done as a pull.
   */
  for (ip = 1; ip <= L; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */

      mmx[ip][0] = dmx[ip][0] = -INFTY;  /*M_0 (B) and D_0 (non-existent)
					 don't emit.
				       */
      imx[ip][0]  = ILogsum(ILogsum(mmx[ip-1][0] + hmm->tsc[CTMI][0],
				    imx[ip-1][0] + hmm->tsc[CTII][0]),
			    dmx[ip-1][0] + hmm->tsc[CTDI][0]);
      imx[ip][0] += hmm->isc[dsq[i]][0];
      for (k = 1; k <= hmm->M; k++)
	{
	  mmx[ip][k]  = ILogsum(ILogsum(mmx[ip-1][k-1] + hmm->tsc[CTMM][k-1],
				       imx[ip-1][k-1] + hmm->tsc[CTIM][k-1]),
			       ILogsum(mmx[ip-1][0] + hmm->bsc[k],
				       dmx[ip-1][k-1] + hmm->tsc[CTDM][k-1]));
	  mmx[ip][k] += hmm->msc[dsq[i]][k];

	  dmx[ip][k]  = ILogsum(ILogsum(mmx[ip][k-1] + hmm->tsc[CTMD][k-1],
				       imx[ip][k-1] + hmm->tsc[CTID][k-1]),
			       dmx[ip][k-1] + hmm->tsc[CTDD][k-1]);

	  imx[ip][k]  = ILogsum(ILogsum(mmx[ip-1][k] + hmm->tsc[CTMI][k],
				       imx[ip-1][k] + hmm->tsc[CTII][k]),
			       dmx[ip-1][k] + hmm->tsc[CTDI][k]);
	  imx[ip][k] += hmm->isc[dsq[i]][k];
	}
      erow[ip] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	erow[ip] = ILogsum(erow[ip], mmx[ip][k] + hmm->esc[k]);
      erow[ip] = ILogsum(erow[ip], dmx[ip][hmm->M] + hmm->tsc[CTDM][hmm->M]); 
      erow[ip] = ILogsum(erow[ip], imx[ip][hmm->M] + hmm->tsc[CTIM][hmm->M]); 
		       /* transition from D_M -> end */
      /*printf("i: %d score(erow[ip]): %f erow[ip]: %d\n", i, Scorify(erow[ip]), erow[ip]);*/
    }		
  sc = erow[L];
  /*printf("F erow[%d]: %d\n", i, erow[L]);*/
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);

  return Scorify(sc);		/* the total Forward score. */
}

/* Function: CP9ViterbiAlign()
 * based on  P7Viterbi() <-- this function's comments below  
 *           from HMMER 2.4devl core_algorithms.c
 *
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's. 
 *           
 * Args:     dsq    - sequence in digitized form
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           hmm    - the model
 *           mx     - reused DP matrix
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
CP9ViterbiAlign(ESL_DSQ *dsq, int i0, int j0, struct cplan9_s *hmm, struct cp9_dpmatrix_s *mx,
	      CP9trace_t **ret_tr)
{
  if(dsq == NULL)
    cm_Fail("in CP9ViterbiAlign(), dsq is NULL.");

  CP9trace_t  *tr;
  int **mmx;
  int **imx;
  int **dmx;
  int **elmx;
  int  *erow;
  int   i,k,c;
  int   sc;
  int   W;		/* subsequence length */
  int   ip;		/* i': relative position in the subsequence  */

  W  = j0-i0+1;		/* the length of the subsequence */

  /* Allocate a DP matrix with 0..W rows, 0..M-1 columns.
   */ 
  ResizeCPlan9Matrix(mx, W, hmm->M, &mmx, &imx, &dmx, &elmx, &erow);
  /* Initialization of the zero row.
   */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  elmx[0][0]= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY; 

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= hmm->M; k++)
    {
      dmx[0][k]  = -INFTY;
      if((sc = mmx[0][k-1] + hmm->tsc[CTMD][k-1]) > dmx[0][k])
	dmx[0][k] = sc;
      if((sc = dmx[0][k-1] + hmm->tsc[CTDD][k-1]) > dmx[0][k])
	dmx[0][k] = sc;
    }

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   */
  for (ip = 1; ip <= W; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      mmx[ip][0] = dmx[ip][0] = -INFTY;  /* M_0 (B) and D_0 (non-existent)
					  * don't emit. */
      elmx[ip][0] = -INFTY; /* there is no EL state for node 0 */
      imx[ip][0] = -INFTY;
      if((sc = mmx[ip-1][0] + hmm->tsc[CTMI][0]) > imx[ip][0])
	imx[ip][0] = sc;
      if((sc = imx[ip-1][0] + hmm->tsc[CTII][0]) > imx[ip][0])
	imx[ip][0] = sc;
      if((sc = dmx[ip-1][0] + hmm->tsc[CTDI][0]) > imx[ip][0])
	imx[ip][0] = sc;
      if(imx[ip][0] != -INFTY)
	imx[ip][0] += hmm->isc[dsq[i]][0];
      else 
	imx[ip][0] = -INFTY;

      for (k = 1; k <= hmm->M; k++)
	{
	  /*match state*/
	  mmx[ip][k] = -INFTY;
	  if((sc = mmx[ip-1][k-1] + hmm->tsc[CTMM][k-1]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  if((sc = imx[ip-1][k-1] + hmm->tsc[CTIM][k-1]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  if((sc = mmx[ip-1][0] + hmm->bsc[k]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  if((sc = dmx[ip-1][k-1] + hmm->tsc[CTDM][k-1]) > mmx[ip][k])
	    mmx[ip][k] = sc;
	  /* Check if we came from an EL state */
	  if(hmm->flags & CPLAN9_EL) /* no need to waste time */
	    {
	      for(c = 0; c < hmm->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
		{
		  /* transition penalty to EL incurred when EL was entered */
		  if((sc = elmx[ip-1][hmm->el_from_idx[k][c]]) > mmx[ip][k])
		    mmx[ip][k] = sc;
		}
	    }
	  if(mmx[ip][k] != -INFTY)
	    mmx[ip][k] += hmm->msc[dsq[i]][k];
	  else 
	    mmx[ip][k] = -INFTY;

	  /*insert state*/
	  imx[ip][k] = -INFTY;
	  if((sc = mmx[ip-1][k] + hmm->tsc[CTMI][k]) > imx[ip][k])
	    imx[ip][k] = sc;
	  if((sc = imx[ip-1][k] + hmm->tsc[CTII][k]) > imx[ip][k])
	    imx[ip][k] = sc;
	  if((sc = dmx[ip-1][k] + hmm->tsc[CTDI][k]) > imx[ip][k])
	    imx[ip][k] = sc;
	  if(imx[ip][k] != -INFTY)
	    imx[ip][k] += hmm->isc[dsq[i]][k];
	  else 
	    imx[ip][k] = -INFTY;

	  /*delete state*/
	  dmx[ip][k] = -INFTY;
	  if((sc = mmx[ip][k-1] + hmm->tsc[CTMD][k-1]) > dmx[ip][k])
	    dmx[ip][k] = sc;
	  if((sc = imx[ip][k-1] + hmm->tsc[CTID][k-1]) > dmx[ip][k])
	    dmx[ip][k] = sc;
	  if((sc = dmx[ip][k-1] + hmm->tsc[CTDD][k-1]) > dmx[ip][k])
	    dmx[ip][k] = sc;

	  /*EL (end-local) state*/
	  elmx[ip][k] = -INFTY;
	  if((hmm->flags & CPLAN9_EL) && hmm->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
							    HMM nodes that map to right half of a MATP_MP) */
	    {
	      if((sc = mmx[ip][k] + hmm->tsc[CTMEL][k]) > elmx[ip][k])
		elmx[ip][k] = sc; /* transitioned from cur node's match state */
	      if((sc = elmx[ip-1][k] + hmm->el_selfsc) > elmx[ip][k])
		elmx[ip][k] = sc; /* transitioned from cur node's EL state emitted ip on transition */
	    }
	}
      erow[ip] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	if ((sc = mmx[ip][k] + hmm->esc[k]) > erow[ip])
	  erow[ip] = sc;
      if ((sc =  dmx[ip][hmm->M] + hmm->tsc[CTDM][hmm->M]) > erow[ip]) /* transition from D_M -> end */
	erow[ip] = sc;
      /* check if we came from an EL */
      if(hmm->flags & CPLAN9_EL) /* no need to waste time */
	{
	  for(c = 0; c < hmm->el_from_ct[hmm->M+1]; c++) /* el_from_ct[hmm->M+1] holds # ELs that can go to END */
	    {
	      /* transition penalty to EL incurred when EL was entered */
	      if((sc = elmx[ip][hmm->el_from_idx[hmm->M+1][c]]) > erow[ip])
		erow[ip] = sc;
	    }
	}
    } 
  sc = erow[W];
  /* printf("returning sc: %d from CP9ViterbiAlign()\n", sc); */
  
  if (ret_tr != NULL) {
    CP9ViterbiTrace(hmm, dsq, i0, j0, mx, &tr);
    /* CP9PrintTrace(stdout, tr, hmm, dsq); */
    *ret_tr = tr;
  }
  return Scorify(sc);		/* the total Viterbi score. */
}

/*********************************************************************
 * Function: CP9BackwardAlign
 * based on  P7Backward() <-- this function's comments below
 *           from HMMER 2.4devl core_algorithms.c
 * 
 * Note:     Differs from CP9_scan.c:CP9Backward() in that this function
 *           does not work in 'scan' mode. CP9_scan.c:CP9Backward() can
 *           do scanning OR alignment, but is somewhat more complex than
 *           this function.
 * 
 * Purpose:  The Backward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
CP9BackwardAlign(ESL_DSQ *dsq, int i0, int j0, struct cplan9_s *hmm, struct cp9_dpmatrix_s **ret_mx)
{
  if(dsq == NULL)
    cm_Fail("in CP9BackwardAlign(), dsq is NULL.");

  struct cp9_dpmatrix_s *mx;
  int **mmx;
  int **imx;
  int **dmx;
  int **elmx;
  int  *erow;
  int   i,k;
  int   sc;

  int   W;		/* subsequence length */
  int   ip;		/* i': relative position in the subsequence  */

  W  = j0-i0+1;		/* the length of the subsequence */

  /* Allocate a DP matrix with 0..W rows, 0..M-1 columns.
   */ 
  mx = AllocCPlan9Matrix(W+1, hmm->M, &mmx, &imx, &dmx, &elmx, &erow);

  /* Initialization of the W row.
   */
  i = j0;

  erow[W] = 0; /*have to end in E*/

  mmx[W][hmm->M] = erow[W] + hmm->esc[hmm->M]; /* M<-E ...                   */
  mmx[W][hmm->M] += hmm->msc[dsq[i]][hmm->M]; /* ... + emitted match symbol */
  /* can't come from I_M b/c we've emitted a single residue, L from M_M */
  imx[W][hmm->M] = erow[W] + hmm->tsc[CTIM][hmm->M];   /* I_M(C)<-E ... */
  imx[W][hmm->M] += hmm->isc[dsq[i]][hmm->M];           /* ... + emitted match symbol */
  dmx[W][hmm->M] = erow[W] + hmm->tsc[CTDM][hmm->M];    /* D_M<-E */
  for (k = hmm->M-1; k >= 1; k--)
    {
      mmx[W][k]  = hmm->esc[k] + erow[W];
      mmx[W][k]  = ILogsum(mmx[W][k], dmx[W][k+1] + hmm->tsc[CTMD][k]);
      mmx[W][k] += hmm->msc[dsq[i]][k];

      imx[W][k] = dmx[W][k+1] + hmm->tsc[CTID][k];
      imx[W][k] += hmm->isc[dsq[i]][k];

      dmx[W][k] = dmx[W][k+1] + hmm->tsc[CTDD][k];
    }
  
  mmx[W][0] = -INFTY; /*M_0 doesn't emit*/
  imx[W][0] = dmx[W][1] + hmm->tsc[CTID][0];
  imx[W][0] += hmm->isc[dsq[i]][hmm->M];    
  dmx[W][0] = -INFTY; /*D_0 doesn't exist*/

  /* Recursion. Done as a pull.
   */

  for (ip = W-1; ip >= 1; ip--) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from j0 down to i0 */
      erow[ip] = -INFTY;
      
      /* Now the main states. Note the boundary conditions at M.
       */
      mmx[ip][hmm->M] = imx[ip+1][hmm->M] + hmm->tsc[CTMI][hmm->M];
      mmx[ip][hmm->M] += hmm->msc[dsq[i]][hmm->M];
      imx[ip][hmm->M] = imx[ip+1][hmm->M] + hmm->tsc[CTII][hmm->M];
      imx[ip][hmm->M] += hmm->isc[dsq[i]][hmm->M];
      dmx[ip][hmm->M] = imx[ip+1][hmm->M] + hmm->tsc[CTDI][hmm->M];  /* * */
      for (k = hmm->M-1; k >= 1; k--)
	{
	  mmx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTMM][k],
				       imx[ip+1][k] + hmm->tsc[CTMI][k]),
				dmx[ip][k+1] + hmm->tsc[CTMD][k]);
	  mmx[ip][k] += hmm->msc[dsq[i]][k];
	  
	  imx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTIM][k],
				       imx[ip+1][k] + hmm->tsc[CTII][k]),
			       dmx[ip][k+1] + hmm->tsc[CTID][k]);
	  imx[ip][k] += hmm->isc[dsq[i]][k];
	  
	  dmx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTDM][k],
				       imx[ip+1][k] + hmm->tsc[CTDI][k]),
			       dmx[ip][k+1] + hmm->tsc[CTDD][k]);
	}

      imx[ip][0]  = ILogsum(ILogsum(mmx[ip+1][1] + hmm->tsc[CTIM][0],
				    imx[ip+1][0] + hmm->tsc[CTII][0]),
			    dmx[ip][1] + hmm->tsc[CTID][0]);
      imx[ip][0] += hmm->isc[dsq[i]][0];
      mmx[ip][0] = -INFTY;
      /*for (k = hmm->M-1; k >= 1; k--)*/ /*M_0 is the B state, it doesn't emit*/
      for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	mmx[ip][0] = ILogsum(mmx[ip][0], mmx[ip+1][k] + hmm->bsc[k]);
      /* 04.17.07 NOTE: ERROR! mmx[ip][k==0] should takes mmx[ip+0(not 1)][k] b/c
       * M_0 does not emit, and ip was accounted for in mmx[ip][k] */
      mmx[ip][0] = ILogsum(mmx[ip][0], imx[ip+1][k] + hmm->tsc[CTMI][k]);
      mmx[ip][0] = ILogsum(mmx[ip][0], dmx[ip][k+1] + hmm->tsc[CTMD][k]);

      dmx[ip][0] = -INFTY; /* D_0 does not exist */
    }
  /* case when ip = 0 */
  erow[0] = -INFTY;
  mmx[0][hmm->M] = -INFTY; /* need seq to get here */
  imx[0][hmm->M] = -INFTY; /* need seq to get here */
  dmx[0][hmm->M] = imx[1][hmm->M] + hmm->tsc[CTDI][hmm->M];  /* * */
  for (k = hmm->M-1; k >= 1; k--)
    {
      mmx[0][k] = -INFTY; /* need seq to get here */
      imx[0][k] = -INFTY; /* need seq to get here */
      dmx[0][k]  = ILogsum(ILogsum(mmx[1][k+1] + hmm->tsc[CTDM][k],
				   imx[1][k] + hmm->tsc[CTDI][k]),
			   dmx[0][k+1] + hmm->tsc[CTDD][k]);
    }
  imx[0][0] = -INFTY; /*need seq to get here*/
  dmx[0][0] = -INFTY; /*D_0 doesn't exist*/

  mmx[0][0] = -INFTY;
  for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
    {
      mmx[0][0] = ILogsum(mmx[0][0], mmx[1][k] + hmm->bsc[k]);
    }
  mmx[0][0] = ILogsum(mmx[0][0], imx[1][0] + hmm->tsc[CTMI][0]);
  mmx[0][0] = ILogsum(mmx[0][0], dmx[0][1] + hmm->tsc[CTMD][0]);

  sc = mmx[0][0];
  /*printf("B final score: %f (i: %d)\n", Scorify(sc), sc);*/

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);

  return Scorify(sc);		/* the total Backward score. */
}

/* Function: cp9_Posterior()
 * based on Ian Holmes' hmmer/src/postprob.c::P7EmitterPosterior()
 *
 * Purpose:  Combines Forward and Backward matrices into a posterior
 *           probability matrix. For emitters (match and inserts) the 
 *           entries in row i of this matrix are the logs of the posterior 
 *           probabilities of each state emitting symbol i of the sequence. 
 *           For non-emitters the entries in row i of this matrix are the 
 *           logs of the posterior probabilities of each state being 'visited' 
 *           when the last emitted residue in the parse was symbol i of the
 *           sequence.
 *           The last point distinguishes this function from P7EmitterPosterior() 
 *           which set all posterior values for for non-emitting states to -INFTY.
 *           The caller must allocate space for the matrix, although the
 *           backward matrix can be used instead (overwriting it will not
 *           compromise the algorithm).
 *
 *           if(did_scan == TRUE) forward/backward run in scan mode, which allow
 *           parses to start/stop at any position of sequence, this changes how
 *           we calculate summed prob of all parses (calculation of 'sc', see code).
 *           
 * Args:     dsq      - sequence in digitized form
 *           i0       - start of target subsequence (often 1, beginning of dsq)
 *           j0       - end of target subsequence (often L, end of dsq)
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           mx       - pre-allocated dynamic programming matrix
 *           did_scan - TRUE if Forward/Backward were run in 'scan' mode, which means
 *                      parses can start and end at any position of the sequence
 *           
 * Return:   void
 */
void
cp9_Posterior(ESL_DSQ *dsq, int i0, int j0,
	      struct cplan9_s *hmm,
	      struct cp9_dpmatrix_s *fmx,
	      struct cp9_dpmatrix_s *bmx,
	      struct cp9_dpmatrix_s *mx,
	      int did_scan)
{
  if(dsq == NULL)
    cm_Fail("in cp9_posterior(), dsq is NULL.");

  int i;
  int k;
  int sc;
  int L;		/* subsequence length */
  int ip;		/* i': relative position in the subsequence  */
  /*float temp_sc;*/

  L  = j0-i0+1;		/* the length of the subsequence */

  if(did_scan) { /* parses could start/stop anywhere */
    sc = -INFTY;
    for (ip = 0; ip <= L; ip++) {
      /*printf("bmx->mmx[i:%d][0]: %d\n", i, bmx->mmx[ip][0]); */
      sc = ILogsum(sc, (bmx->mmx[ip][0])); 
    }
  } /* parses must start/stop at (i = i0)/(j = j0) */
  else sc = bmx->mmx[0][0];

  /* note boundary conditions, case by case by case... */
  mx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc; /* fmx->mmx[0][0] is 0, bmx->mmx[1][0] is overall score */
  mx->imx[0][0] = -INFTY; /*need seq to get here*/
  mx->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  for (k = 1; k <= hmm->M; k++) {
      mx->mmx[0][k] = -INFTY; /*need seq to get here*/
      mx->imx[0][k] = -INFTY; /*need seq to get here*/
      mx->dmx[0][k] = fmx->dmx[0][k] + bmx->dmx[0][k] - sc;
  }
      
  for (ip = 1; ip <= L; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      mx->mmx[ip][0] = -INFTY; /*M_0 does not emit*/
      mx->imx[ip][0] = fmx->imx[ip][0] + bmx->imx[ip][0] - hmm->isc[dsq[i]][0] - sc;
      /*hmm->isc[dsq[i]][0] will have been counted in both fmx->imx and bmx->imx*/
      mx->dmx[ip][0] = -INFTY; /*D_0 does not exist*/

      /*printf("fmx->mmx[ip:%d][0]: %d\n bmx->mmx[ip:%d][0]: %d\n", ip, fmx->mmx[ip][0], ip, bmx->mmx[ip][0]);
	printf("fmx->imx[ip:%d][0]: %d\n bmx->imx[ip:%d][0]: %d\n", ip, fmx->imx[ip][0], ip, bmx->imx[ip][0]);
	printf("fmx->dmx[ip:%d][0]: %d\n bmx->dmx[ip:%d][0]: %d\n", ip, fmx->dmx[ip][0], ip, bmx->dmx[ip][0]);*/
      for (k = 1; k <= hmm->M; k++) 
	{
	  mx->mmx[ip][k] = ESL_MAX(fmx->mmx[ip][k] + bmx->mmx[ip][k] - hmm->msc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  mx->imx[ip][k] = ESL_MAX(fmx->imx[ip][k] + bmx->imx[ip][k] - hmm->isc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
	  mx->dmx[ip][k] = ESL_MAX(fmx->dmx[ip][k] + bmx->dmx[ip][k] - sc, -INFTY);
	  /*printf("fmx->mmx[ip:%d][%d]: %d\n bmx->mmx[ip:%d][%d]: %d\n", ip, k, fmx->mmx[ip][k], ip, k, bmx->mmx[ip][k]);
	  printf("fmx->imx[ip:%d][%d]: %d\n bmx->imx[ip:%d][%d]: %d\n", ip, k, fmx->imx[ip][k], ip, k, bmx->imx[ip][k]);
	  printf("fmx->dmx[ip:%d][%d]: %d\n bmx->dmx[ip:%d][%d]: %d\n\n", ip, k, fmx->dmx[ip][k], ip, k, bmx->dmx[ip][k]);*/
	}	  
    }

  /*
    float temp_sc;
    for(i = 0; i <= L; i++)
    {
    for(k = 0; k <= hmm->M; k++)
    {
    temp_sc = Score2Prob(mx->mmx[i][k], 1.);
    if(temp_sc > .0001)
    printf("mx->mmx[%3d][%3d]: %9d | %8f\n", i, k, mx->mmx[i][k], temp_sc);
    temp_sc = Score2Prob(mx->imx[i][k], 1.);
    if(temp_sc > .0001)
    printf("mx->imx[%3d][%3d]: %9d | %8f\n", i, k, mx->imx[i][k], temp_sc);
    temp_sc = Score2Prob(mx->dmx[i][k], 1.);
    if(temp_sc > .0001)
    printf("mx->dmx[%3d][%3d]: %9d | %8f\n", i, k, mx->dmx[i][k], temp_sc);
    }
    }*/
}

/*****************************************************************************
 * EPN 03.23.06
 * Function: cp9_IFillPostSums()
 * based on: ifill_post_sums_del() (deprecated) 11.23.05
 *
 * Purpose:  Given a posterior matrix post, where post->mmx[i][k]
 *           is the log odds score of the probability that
 *           match state k emitted position i of the sequence,
 *           sum the log probabilities that each state emitted
 *           each position. Do this for inserts, matches, and
 *           and deletes.
 * 
 * arguments:
 * cp9_dpmatrix_s *post  dpmatrix_s posterior matrix, xmx, mmx, imx, dmx 
 *                       2D int arrays. [0.1..N][0.1..M]
 * CP9Bands_t *cp9b - the cp9 bands data structure
 * int  i0          start of target subsequence (often 1, beginning of dsq)
 * int  j0          end of target subsequence (often L, end of dsq)
 *****************************************************************************/
void
cp9_IFillPostSums(struct cp9_dpmatrix_s *post, CP9Bands_t *cp9b, int i0, int j0)
{
  int i;            /* counter over positions of the sequence */
  int k;            /* counter over nodes of the model */
  int L;	    /* subsequence length */
  int M;            /* consensus length of cp9 */
  M = cp9b->hmm_M;
  L  = j0-i0+1;		/* the length of the subsequence */
  
  /* step through each node, fill the post sum structures */
  for(k = 0; k <= M; k++)
    {
      cp9b->isum_pn_m[k] = -INFTY;
      cp9b->isum_pn_i[k] = -INFTY;
      cp9b->isum_pn_d[k] = -INFTY;
      for(i = 0; i <= L; i++) {
	cp9b->isum_pn_m[k] = ILogsum(cp9b->isum_pn_m[k], post->mmx[i][k]);
	cp9b->isum_pn_i[k] = ILogsum(cp9b->isum_pn_i[k], post->imx[i][k]);
	cp9b->isum_pn_d[k] = ILogsum(cp9b->isum_pn_d[k], post->dmx[i][k]);
      }
    }
}

/*****************************************************************************
 * Functions to go from HMM bands to i and j bands on a CM 
 * cp9_HMM2ijBands()
 */
/*
 * Function: cp9_HMM2ijBands()
 *           EPN 12.21.05
 * 
 * Purpose:  Determine the band for each cm state v on i (the band on the 
 *           starting index in the subsequence emitted from the subtree rooted
 *           at state v), and on j (the band on the ending index in the
 *           subsequence emitted from the subtree rooted at state v). 
 * 
 *           Some i and d bands are calculated from HMM bands on match and insert 
 *           and delete states from each node of the HMM that maps to a left emitting
 *           node of the CM (including MATP nodes). The HMM bands were
 *           calculated previously from the posterior matrices for mmx,
 *           imx and dmx from a CP9 HMM.
 * 
 *           Some j bands are calculated from HMM bands on match and insert and
 *           delete states from each node of the HMM that maps to a right emitting
 *           node of the CM (including MATP nodes). 
 * 
 *           i and j bands that cannot be directly determined from the
 *           HMM bands are inferred based on the constraints imposed
 *           on them by the i and j bands that CAN be determined from
 *           the HMM bands.
 *             
 *           Specifically, this function implements strategy 3 (EPN)
 *           of the hmm bands to i and j bands problem. The main idea
 *           of strategy 3 is to set i and j bands for each state v
 *           such that at least one state y (y \in C_v (y is reachable
 *           from v)) can be reached from v while staying within the i
 *           and j bands for v and y.  This constraint is enforced by
 *           determining the min and max i and j bands across all
 *           states y (into safe* data structures) for a given v, and
 *           then enforcing that at least one cell in the i and j
 *           bands of v can transit to at least one cell in a band for
 *           a y state after accounting for the direction specific
 *           StateDelta() values for v.
 *           
 *           This function needs to be called only once, it determines
 *           bands for ALL states. Its unclear the best way to handle
 *           any states that don't have an explicit mapping to an HMM
 *           state that we have a band on (i.e. all delete states, and
 *           ROOT_IR, ROOT_IL, BEGR_IL, BIF_B, and start states).
 *           (11.02.05) I take a simple approach, and set the bands on i
 *           for such states to the same as those for states in a close
 *           proximity. (see code for exact definitions)
 * 
 *           This function uses HMM derived bands on delete states I'm
 *           not sure if the method used to get these delete states is
 *           100% sound (as it required me writing a new function
 *           P7FullPosterior() to derive them from the HMMER forwards
 *           and backwards parses).  
 *
 * arguments:
 *
 * CM_t *cm         the CM 
 * errbuf           char buffer for error messages
 * CP9Map_t *cp9map map from CM to CP9 HMM and vice versa
 * int i0           start of target subsequence (often 1, beginning of dsq)
 * int j0           end of target subsequence (often L, end of dsq)
 * int *pn_min_m    pn_min_m[k] = first position in HMM band for match state of HMM node k
 * int *pn_max_m    pn_max_m[k] = last position in HMM band for match state of HMM node k
 * int *pn_min_i    pn_min_i[k] = first position in HMM band for insert state of HMM node k
 * int *pn_max_i    pn_max_i[k] = last position in HMM band for insert state of HMM node k
 * int *pn_min_d    pn_min_d[k] = first position in HMM band for delete state of HMM node k
 * int *pn_max_d    pn_max_d[k] = last position in HMM band for delete state of HMM node k
 * int *imin        imin[v] = first position in band on i for state v
 *                  to be filled in this function. [1..M]
 * int *imax        imax[v] = last position in band on i for state v
 *                  to be filled in this function. [1..M]
 * int *jmin        jmin[v] = first position in band on j for state v
 *                  to be filled in this function. [1..M]
 * int *jmax        jmax[v] = last position in band on j for state v
 *                  to be filled in this function. [1..M]
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *
 * Returns: eslOK on success;
 */
int
cp9_HMM2ijBands(CM_t *cm, char *errbuf, CP9Map_t *cp9map, int i0, int j0, int *pn_min_m, 
		int *pn_max_m, int *pn_min_i, int *pn_max_i, int *pn_min_d, 
		int *pn_max_d, int *imin, int *imax, int *jmin, int *jmax, 
		int debug_level)
{
  int status;
  int v;            /* counter over states of the CM */

  int *nss_imin;      /* nss_imin[n] = imin of each split set state in node n*/
  int *nss_imax;      /* nss_imax[n] = imax of each split set state in node n*/
  int *nss_jmin;      /* nss_jmin[n] = jmin of each split set state in node n*/
  int *nss_jmax;      /* nss_jmax[n] = jmax of each split set state in node n*/

  int *nis_imin;      /* nss_imin[n] = imin of each insert set state in node n*/
  int *nis_imax;      /* nss_imax[n] = imax of each insert set state in node n*/
  int *nis_jmin;      /* nss_jmin[n] = jmin of each insert set state in node n*/
  int *nis_jmax;      /* nss_jmax[n] = jmax of each insert set state in node n*/

  int *nss_max_imin;  /* nss_max_imin[n] = max imin over split set states in node n*/
  int *nss_min_jmax;  /* nss_min_jmax[n] = min jmax over split set states in node n*/

  int safe_imax; 
  int safe_jmin; 

  int tmp_imin;
  int tmp_imax;
  int tmp_jmin;
  int tmp_jmax;

  int n;            /* counter over CM nodes. */
  int doing_search; /* TRUE if bands will be used to accelerate search, FALSE for alignment */
  /* Contract checks */
  if((cm->align_opts & CM_ALIGN_HBANDED) && (cm->search_opts & CM_SEARCH_HBANDED))    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_HMM2ijBands(), CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both up, exactly 1 must be up.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED))) ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_HMM2ijBands(), CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down, exactly 1 must be up.\n");

  doing_search = (cm->search_opts & CM_SEARCH_HBANDED) ? TRUE : FALSE;

  ESL_ALLOC(nss_imin, sizeof(int) * cm->nodes);
  ESL_ALLOC(nss_imax, sizeof(int) * cm->nodes);
  ESL_ALLOC(nss_jmin, sizeof(int) * cm->nodes);
  ESL_ALLOC(nss_jmax, sizeof(int) * cm->nodes);

  ESL_ALLOC(nis_imin, sizeof(int) * cm->nodes);
  ESL_ALLOC(nis_imax, sizeof(int) * cm->nodes);
  ESL_ALLOC(nis_jmin, sizeof(int) * cm->nodes);
  ESL_ALLOC(nis_jmax, sizeof(int) * cm->nodes);

  ESL_ALLOC(nss_max_imin, sizeof(int) * cm->nodes);
  ESL_ALLOC(nss_min_jmax, sizeof(int) * cm->nodes);

  /* Initialize all bands to -1. */
  esl_vec_ISet(imin, cm->M, -1);
  esl_vec_ISet(imax, cm->M, -1);
  esl_vec_ISet(jmin, cm->M, -1);
  esl_vec_ISet(jmax, cm->M, -1);

  /* We go node by node, bottom up, and fill in the bands on each
   * state for each node. Keeping track of the node split set min and max i's 
   * and j's, as well as the node insert set's
   * also because they influence all nodes above (until a BEGL or BEGR at least).
   */

  /* For match nodes (MATP, MATL, MATR):
   * First calc the split set node mins and maxes, then impose these
   * on each state v in the split set of the node, requiring that any valid
   * d resulting from the i and j bands on state v
   * is least dv = StateDelta(v).
   * This is done by ensuring that jmin[v] >= dv & jmax[v] >= dv.
   * (We don't have to worry about i as we check again when we create
   *  the d bands from the i and j bands in ij2d_bands()).
   * We really only have to enforce the StateDelta issue here so we 
   * don't run into d band on j that is 0 cells in ij2d_bands(). 
   * Alternatively, we could ignore the StateDelta() issue here, and
   * allow ij2d_bands() to modify j bands when it enforces the StateDelta()
   * issue.
   */

   for(n = (cm->nodes-1); n >= 0; n--) {
     switch (cm->ndtype[n]) { 
     case END_nd:
       /* Special case, we need to know the bands on the states
	* in the node ABOVE this one. Node above MUST be MATP, MATL
	* or MATR. For END states, the band on i = the band on j,
	* this is because d must be 0, so i must be (j+1), so its pointless
	* to allow an i value that (j+1) is not allowed to be or vice versa.
	* If the node above is MATL, we use the HMM band that maps
	* to the ML state - these correspond to bands on i. If its a MATR, 
	* we use the HMM band that maps to the MR state - these correspond
	* to bands on j. If its a MATP, we get fancy (see below).
	*/
       v = cm->nodemap[n];
       if(cm->ndtype[n-1] == MATL_nd) {
	 /* tricky. we keep the n_*m** structures ignorant of the fact that we're in
	  * an end state, i.e. we don't force a d=0 (j-i+1=0). This way when
	  * the node immediately above the end (the MATL) looks at it when its determining
	  * the correct bands on i, it doesn't get screwed up (as it would if j < i).
	  */
	 
	 /*minimum of delete and match states of node above*/
	 nss_imin[n] = (pn_min_m[cp9map->nd2lpos[n-1]] <= (pn_min_d[cp9map->nd2lpos[n-1]])) ? 
	   pn_min_m[cp9map->nd2lpos[n-1]] : (pn_min_d[cp9map->nd2lpos[n-1]]);
	 /*for the max, we must allow possibility of inserts and deletes.*/
	 nss_imax[n] = (pn_max_m[cp9map->nd2lpos[n-1]] >= pn_max_i[cp9map->nd2lpos[n-1]]) ? 
	   pn_max_m[cp9map->nd2lpos[n-1]] : pn_max_i[cp9map->nd2lpos[n-1]];
	 /* deletes max bands may always be less than match max bands...(not sure)*/
	 if(nss_imax[n] < (pn_max_d[cp9map->nd2lpos[n-1]]))
	   nss_imax[n] = (pn_max_d[cp9map->nd2lpos[n-1]]);
	 
	 nss_jmin[n] = nss_imin[n];
	 nss_jmax[n] = nss_imax[n];
	 
	 imin[v] = nss_imin[n];
	 imax[v] = nss_imax[n] + 1; /* we add 1 because we have to figure in the emission
				     * of the MATL_ML (or final MATL_IL), which would increase
				     * i by 1 potentially relative to the imax of that state.
				     */
	 jmin[v] = imin[v] - 1; /* d must be 0 for end states. */
	 jmax[v] = imax[v] - 1; /* d must be 0 for end states. */
	 
	 nss_max_imin[n] = imin[v];
	 nss_min_jmax[n] = jmax[v];
       }
       else if(cm->ndtype[n-1] == MATR_nd) {
	 /* tricky. we keep the nss_*m** structures ignorant of the fact that we're in
	  * an end state, i.e. we don't force a d=0 (j-i+1=0). This way when
	  * the node immediately above the end (the MATR) looks at it when its determining
	  * the correct bands on i, it doesn't get screwed up (as it would if j < i).
	  */
	 
	 /*minimum of delete and match states of node above */
	 nss_jmin[n] = (pn_min_m[cp9map->nd2rpos[n-1]] <= pn_min_d[cp9map->nd2rpos[n-1]]) ? 
	   pn_min_m[cp9map->nd2rpos[n-1]] : pn_min_d[cp9map->nd2rpos[n-1]];
	 /*for the max, we must allow possibility of inserts.*/
	 nss_jmax[n] = (pn_max_m[cp9map->nd2rpos[n-1]] >= pn_max_i[cp9map->nd2rpos[n-1]]) ? 
	   pn_max_m[cp9map->nd2rpos[n-1]] : pn_max_i[cp9map->nd2rpos[n-1]];
	 /* deletes max bands may always be less than match max bands...(not sure)*/
	 if(nss_jmax[n] < pn_max_d[cp9map->nd2rpos[n-1]])
	   nss_jmax[n] = pn_max_d[cp9map->nd2rpos[n-1]];
	 nss_imin[n] = nss_jmin[n];
	 nss_imax[n] = nss_jmax[n];
	 
	 jmin[v] = nss_jmin[v] - 1; /* we subtract 1 because of we have to figure 
				     * in the emission of the MATR_MR (or final MATR_IR), which would 
				     * decrease j by 1 potentially relative to jmin of that state.
				     */
	 jmax[v] = nss_jmax[n];
	 imin[v] = jmin[v] + 1; /*d (j-i+1) must be 0 for end states*/
	 imax[v] = jmax[v] + 1; /*d (j-i+1) must be 0 for end states*/
	 
	 nss_max_imin[n] = imin[v];
	 nss_min_jmax[n] = jmax[v];
       }
       else if(cm->ndtype[n-1] == MATP_nd) { 
	 /* Very rare case, only if the last bp in a stem is the last left consensus
	  * column (respecting gap_thresh) in that alignment. Does happen though, 
	  * (at least in RFAM 6.1) because the training counts for transition priors
	  * had counts for MATP_* state -> END_nd transition sets.
	  */
	 
	 /* tricky. we keep the nss_*m** structures ignorant of the fact that we're in
	  * an end state, i.e. we don't force a d=0 (j-i+1=0). This way when
	  * the node immediately above the end (the MATP) looks at it when its determining
	  * the correct bands on j, it doesn't get screwed up (as it would if j < i).
	  */
	 /*minimum of delete and match states of node above*/
	 nss_imin[n] = (pn_min_m[cp9map->nd2lpos[n-1]] <= (pn_min_d[cp9map->nd2lpos[n-1]])) ? 
	   pn_min_m[cp9map->nd2lpos[n-1]] : (pn_min_d[cp9map->nd2lpos[n-1]]);
	 /*for the max, we must allow possibility of inserts and deletes.*/
	 nss_imax[n] = (pn_max_m[cp9map->nd2lpos[n-1]] >= pn_max_i[cp9map->nd2lpos[n-1]]) ? 
	   pn_max_m[cp9map->nd2lpos[n-1]] : pn_max_i[cp9map->nd2lpos[n-1]];
	 /* deletes max bands may always be less than match max bands...(not sure)*/
	 if(nss_imax[n] < (pn_max_d[cp9map->nd2lpos[n-1]]))
	   nss_imax[n] = (pn_max_d[cp9map->nd2lpos[n-1]]);
	 
	 /*minimum of delete and match states of node above*/
	 nss_jmin[n] = (pn_min_m[cp9map->nd2rpos[n-1]] <= pn_min_d[cp9map->nd2rpos[n-1]]) ? 
	   pn_min_m[cp9map->nd2rpos[n-1]] : pn_min_d[cp9map->nd2rpos[n-1]];
	 /*for the max, we must allow possibility of inserts.*/
	 nss_jmax[n] = (pn_max_m[cp9map->nd2rpos[n-1]] >= pn_max_i[cp9map->nd2rpos[n-1]]) ? 
	   pn_max_m[cp9map->nd2rpos[n-1]] : pn_max_i[cp9map->nd2rpos[n-1]];
	 /* deletes max bands may always be less than match max bands...(not sure)*/
	 if(nss_jmax[n] < pn_max_d[cp9map->nd2rpos[n-1]])
	   nss_jmax[n] = pn_max_d[cp9map->nd2rpos[n-1]];
	 
	 /* unique situation. end's d must be 0, so we are constrained on what 
	  * i can be relative to j, and j can be relative to i, but what we want
	  * are the constraints on what i can be, and j can be. 
	  * because d=0 => j-i+1 = 0. then imin should equal = jmin + 1 and imax = jmax + 1.
	  * so we really just want to know a min over i and j, and a max over i and j.
	  * below we take min of imin and jmin (should always be imin i think) as the min, 
	  * and max of imax and jmax (should always be jmax i think) after accounting for 
	  * the possibility that a single base was just emitted left and/or right.
	  */
	 imax[v] = ((nss_imax[n] + 1) > nss_jmax[n]) ? 
	   (nss_imax[n] + 1) : nss_jmax[n];
	 imin[v] = ((nss_imin[n]) < (nss_jmin[n] - 1)) ? 
	   (nss_imin[n]) : (nss_jmin[n] - 1);
	 /* we can't have an i < i0 */
	 imin[v] = ESL_MAX(imin[v], i0);
	 imax[v] = ESL_MAX(imax[v], i0);
	 jmin[v] = imin[v] - 1; /* d must be 0 for end states. */
	 jmax[v] = imax[v] - 1; /* d must be 0 for end states. */

	 nss_max_imin[n] = imin[v];
	 nss_min_jmax[n] = jmax[v];
       }
       break;

	case MATP_nd:
	  hmm2ij_prestate_step0_initialize(n, nss_max_imin, nss_min_jmax, i0, j0);
	  hmm2ij_prestate_step1_set_node_inserts(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						 nss_imin, nss_imax, nss_jmin, nss_jmax,
						 pn_min_i, pn_max_i, cp9map);

	  hmm2ij_prestate_step2_determine_safe(n, nss_max_imin[n+1], nss_min_jmax[n+1],
					       nis_imin[n], nis_jmax[n],
					       &safe_imax, &safe_jmin);
	  hmm2ij_prestate_step3_preset_node_splits(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax,
						   pn_min_m, pn_max_m, pn_min_d, pn_max_d,
						   cp9map);
	  /* 6 states MATP_MP, MATP_ML, MATP_MR, MATP_D, MATP_IL, MATP_IR */
	  v = cm->nodemap[n]; /* MATP_MP */
	  /* Determine implied v bands using hmm for mapped 'direction(s)' and 
	   * next node's bands for non-mapped direction(s).
	   */
	  tmp_imin = pn_min_m[cp9map->nd2lpos[n]]; 
	  tmp_imax = pn_max_m[cp9map->nd2lpos[n]];
	  tmp_jmin = pn_min_m[cp9map->nd2rpos[n]];
	  tmp_jmax = pn_max_m[cp9map->nd2rpos[n]];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);

	  v++; /*MATP_ML*/
	  /* Determine implied v bands using hmm for mapped 'direction(s)' and 
	   * next node's bands for non-mapped direction(s).
	   */
	  tmp_imin = pn_min_m[cp9map->nd2lpos[n]]; 
	  tmp_imax = pn_max_m[cp9map->nd2lpos[n]];
	  /* 12.19.05 - trying to deal with the right delete off-by-one
	   * inverted relative to left delete issue.
	   */
	  tmp_jmin = (pn_min_d[cp9map->nd2rpos[n]] < nss_jmin[n+1]) ?
	    pn_min_d[cp9map->nd2rpos[n]] : nss_jmin[n+1];
	  tmp_jmax = (pn_max_d[cp9map->nd2rpos[n]] > nss_jmax[n+1]) ? 
	    pn_max_d[cp9map->nd2rpos[n]] : nss_jmax[n+1]; 

	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  
	  v++; /*MATP_MR*/
	  /* this D-left state gets the delete band from the HMM node
	   * that maps to the left side.
	   */
	  tmp_imin = pn_min_d[cp9map->nd2lpos[n]]; 
	  tmp_imax = pn_max_d[cp9map->nd2lpos[n]];
	  tmp_jmin = pn_min_m[cp9map->nd2rpos[n]];
	  tmp_jmax = pn_max_m[cp9map->nd2rpos[n]];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  
	  v++; /*MATP_D*/
	  tmp_imin = pn_min_d[cp9map->nd2lpos[n]]; 
	  tmp_imax = pn_max_d[cp9map->nd2lpos[n]];
	  /* 12.19.05 - trying to deal with the right delete off-by-one
	   * inverted relative to left delete issue.
	   */
	  tmp_jmin = (pn_min_d[cp9map->nd2rpos[n]] < nss_jmin[n+1]) ?
	    pn_min_d[cp9map->nd2rpos[n]] : nss_jmin[n+1];
	  tmp_jmax = (pn_max_d[cp9map->nd2rpos[n]] > nss_jmax[n+1]) ? 
	    pn_max_d[cp9map->nd2rpos[n]] : nss_jmax[n+1]; 
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  hmm2ij_state_step5_non_emitter_d0_hack(v, imax[v], jmin);
	  
	  v++; /*MATP_IL*/
	  /* This state maps to the insert state of HMM node cp9map->cs2hn[v][0]*/
	  tmp_imin = pn_min_i[cp9map->cs2hn[v][0]]; /* insert states can only map to 1 HMM node */
	  tmp_imax = pn_max_i[cp9map->cs2hn[v][0]]; /* insert states can only map to 1 HMM node */
	  tmp_jmin = nss_jmin[n];
	  tmp_jmax = nss_jmax[n];
	  hmm2ij_insert_state_step1_set_state_bands(v, tmp_imin, tmp_imax, tmp_jmin,
						    tmp_jmax, imin, imax, jmin, jmax);
	  /* Enforce safe transitions, this makes sure that at least one state
	   * y \in C_v is reachable from v. And further (special case for inserts)
	   * make sure that we don't consider v as a possible y.  IF we did, we might
	   * be faced with a situation where v could only transit to itself, and then
	   * we'd be in the same situation, but we may no longer be able to transit ANYWHERE
	   * including to itself.
	   */
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, (nss_max_imin[n+1]), nss_min_jmax[n+1]);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  
	  v++; /*MATP_IR*/
	  /* skip detached inserts */
	  if(cp9map->cs2hn[v][0] == -1)
	    continue;
	  /* Special case, one of only two situations (other is ROOT_IR)
	   * we could have come where v is an insert, and a possible
	   * state x that we came from is an insert, but x != y (x can be the MATP_IL).
	   * So we have to determine imin and imax carefully.
	   */
	  tmp_imin = (nss_imin[n] < imin[v-1]) ? 
	    nss_imin[n] : imin[v-1];
	  tmp_imax = (nss_imax[n] > imax[v-1]) ? 
	    nss_imax[n] : imax[v-1];
	  /* This state maps to the insert state of HMM node cp9map->cs2hn[v][0]*/
	  tmp_jmin = pn_min_i[cp9map->cs2hn[v][0]]; 
	  tmp_jmax = pn_max_i[cp9map->cs2hn[v][0]];
	  hmm2ij_insert_state_step1_set_state_bands(v, tmp_imin, tmp_imax, tmp_jmin,
						    tmp_jmax, imin, imax, jmin, jmax);
	  /* Enforce safe transitions, this makes sure that at least one state
	   * y \in C_v is reachable from v. And further (special case for inserts)
	   * make sure that we don't consider v as a possible y.  IF we did, we might
	   * be faced with a situation where v could only transit to itself, and then
	   * we'd be in the same situation, but we may no longer be able to transit ANYWHERE
	   * including to itself.
	   */
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, (nss_max_imin[n+1]), nss_min_jmax[n+1]);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  break;

	case MATL_nd:
	  hmm2ij_prestate_step0_initialize(n, nss_max_imin, nss_min_jmax, i0, j0);
	  hmm2ij_prestate_step1_set_node_inserts(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						 nss_imin, nss_imax, nss_jmin, nss_jmax,
						 pn_min_i, pn_max_i, cp9map);

	  hmm2ij_prestate_step2_determine_safe(n, nss_max_imin[n+1], nss_min_jmax[n+1],
					       nis_imin[n], nis_jmax[n],
					       &safe_imax, &safe_jmin);
	  hmm2ij_prestate_step3_preset_node_splits(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax,
						   pn_min_m, pn_max_m, pn_min_d, pn_max_d,
						   cp9map);

	  /* 3 states MATL_ML, MATL_D, MATL_IL */
	  v = cm->nodemap[n]; /* MATL_ML */
	  tmp_imin = pn_min_m[cp9map->nd2lpos[n]]; 
	  tmp_imax = pn_max_m[cp9map->nd2lpos[n]];
	  tmp_jmin = nss_jmin[n+1];
	  tmp_jmax = nss_jmax[n+1];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);

	  v++; /*MATL_D*/
	  /* this D-left state gets the delete band from the HMM node
	   * that maps to the left side.
	   */
	  tmp_imin = pn_min_d[cp9map->nd2lpos[n]]; 
	  tmp_imax = pn_max_d[cp9map->nd2lpos[n]];
	  tmp_jmin = nss_jmin[n+1];
	  tmp_jmax = nss_jmax[n+1];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  hmm2ij_state_step5_non_emitter_d0_hack(v, imax[v], jmin);

	  v++; /*MATL_IL*/
	  /* skip detached inserts */
	  if(cp9map->cs2hn[v][0] == -1)
	    continue;
	  /* This state maps to the insert state of HMM node cp9map->cs2hn[v][0]*/
	  tmp_imin = pn_min_i[cp9map->cs2hn[v][0]];
	  tmp_imax = pn_max_i[cp9map->cs2hn[v][0]];
	  tmp_jmin = nss_jmin[n];
	  tmp_jmax = nss_jmax[n];

	  hmm2ij_insert_state_step1_set_state_bands(v, tmp_imin, tmp_imax, tmp_jmin,
						    tmp_jmax, imin, imax, jmin, jmax);
	  /* Enforce safe transitions, this makes sure that at least one state
	   * y \in C_v is reachable from v. And further (special case for inserts)
	   * make sure that we don't consider v as a possible y.  IF we did, we might
	   * be faced with a situation where v could only transit to itself, and then
	   * we'd be in the same situation, but we may no longer be able to transit ANYWHERE
	   * including to itself.
	   */
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, (nss_max_imin[n+1]), nss_min_jmax[n+1]);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  break;

	case MATR_nd:
	  hmm2ij_prestate_step0_initialize(n, nss_max_imin, nss_min_jmax, i0, j0);
	  hmm2ij_prestate_step1_set_node_inserts(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						 nss_imin, nss_imax, nss_jmin, nss_jmax,
						 pn_min_i, pn_max_i, cp9map);

	  hmm2ij_prestate_step2_determine_safe(n, nss_max_imin[n+1], nss_min_jmax[n+1],
					       nis_imin[n], nis_jmax[n],
					       &safe_imax, &safe_jmin);
	  hmm2ij_prestate_step3_preset_node_splits(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax,
						   pn_min_m, pn_max_m, pn_min_d, pn_max_d,
						   cp9map);

	  /* 3 states MATR_MR, MATR_D, MATR_IR */
	  v = cm->nodemap[n]; /* MATR_MR */
	  tmp_imin = nss_imin[n+1];
	  tmp_imax = nss_imax[n+1];
	  tmp_jmin = pn_min_m[cp9map->nd2rpos[n]]; 
	  tmp_jmax = pn_max_m[cp9map->nd2rpos[n]];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  
	  v++; /*MATR_D*/
	  /* this D-left state gets the delete band from the HMM node
	   * that maps to the left side.
	   */
	  tmp_imin = nss_imin[n+1];
	  tmp_imax = nss_imax[n+1];
	  /* 12.19.05 - trying to deal with the right delete off-by-one
	   * inverted relative to left delete issue.
	   */
	  tmp_jmin = (pn_min_d[cp9map->nd2rpos[n]] < nss_jmin[n+1]) ?
	    pn_min_d[cp9map->nd2rpos[n]] : nss_jmin[n+1];
	  tmp_jmax = (pn_max_d[cp9map->nd2rpos[n]] > nss_jmax[n+1]) ? 
	    pn_max_d[cp9map->nd2rpos[n]] : nss_jmax[n+1]; 
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  hmm2ij_state_step5_non_emitter_d0_hack(v, imax[v], jmin);

	  v++; /*MATR_IR*/
	  /* skip detached inserts */
	  if(cp9map->cs2hn[v][0] == -1)
	    continue;
	  tmp_imin = nss_imin[n];
	  tmp_imax = nss_imax[n];
	  /* This state maps to the insert state of HMM node cshn_map[v]*/
	  tmp_jmin = pn_min_i[cp9map->cs2hn[v][0]];
	  tmp_jmax = pn_max_i[cp9map->cs2hn[v][0]];
	  hmm2ij_insert_state_step1_set_state_bands(v, tmp_imin, tmp_imax, tmp_jmin,
						    tmp_jmax, imin, imax, jmin, jmax);
	  /* Enforce safe transitions, this makes sure that at least one state
	   * y \in C_v is reachable from v. And further (special case for inserts)
	   * make sure that we don't consider v as a possible y.  IF we did, we might
	   * be faced with a situation where v could only transit to itself, and then
	   * we'd be in the same situation, but we may no longer be able to transit ANYWHERE
	   * including to itself.
	   */
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, (nss_max_imin[n+1]), nss_min_jmax[n+1]);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  break;

	case ROOT_nd:
	  hmm2ij_prestate_step0_initialize(n, nss_max_imin, nss_min_jmax, i0, j0);
	  hmm2ij_prestate_step1_set_node_inserts(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						 nss_imin, nss_imax, nss_jmin, nss_jmax,
						 pn_min_i, pn_max_i, cp9map);

	  hmm2ij_prestate_step2_determine_safe(n, nss_max_imin[n+1], nss_min_jmax[n+1],
					       nis_imin[n], nis_jmax[n],
					       &safe_imax, &safe_jmin);
	  hmm2ij_prestate_step3_preset_node_splits(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax,
						   pn_min_m, pn_max_m, pn_min_d, pn_max_d,
						   cp9map);
	  /* 3 states, ROOT_S, ROOT_IL, and ROOT_IR*/
	  v = cm->nodemap[n]; /* ROOT_S SPECIAL CASE */
	  if(doing_search) { /* we're doing search, ROOT_S doesn't necessarily emit full sequence */
	    tmp_imin = nss_imin[n+1];
	    tmp_imax = nss_imax[n+1];
	    tmp_jmin = nss_jmin[n+1];
	    tmp_jmax = nss_jmax[n+1];
	  }
	  else { /* we're doing alignment, enforce ROOT_S emits full sequence */
	    tmp_imin = i0;
	    tmp_imax = i0;
	    tmp_jmin = j0;
	    tmp_jmax = j0;
	  }
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  
	  v++; /*ROOT_IL SPECIAL CASE*/
	  /* This state maps to the insert state of HMM node cp9map->cs2hn[v][0], which is HMM node 0*/
	  if(doing_search)
	    tmp_imin =  pn_min_i[cp9map->cs2hn[v][0]]; /* should this be imin[0]? */
	  else
	    tmp_imin =  i0; /* Have to be able to transit here from ROOT_S */
	  tmp_imax = nss_imax[n+1];
	  if(doing_search) {
	    tmp_jmin = nss_jmin[n+1];
	    tmp_jmax = nss_jmax[n+1];
	  }
	  else {
	    tmp_jmin = j0; /* we never emit to the right in this state */
	    tmp_jmax = j0; /* we never emit to the right in this state */
	  }
	  hmm2ij_insert_state_step1_set_state_bands(v, tmp_imin, tmp_imax, tmp_jmin,
						    tmp_jmax, imin, imax, jmin, jmax);
	  /* Enforce safe transitions, this makes sure that at least one state
	   * y \in C_v is reachable from v. And further (special case for inserts)
	   * make sure that we don't consider v as a possible y.  IF we did, we might
	   * be faced with a situation where v could only transit to itself, and then
	   * we'd be in the same situation, but we may no longer be able to transit ANYWHERE
	   * including to itself.
	   */
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, (nss_max_imin[n+1]), nss_min_jmax[n+1]);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);

	  v++; /*ROOT_IR SPECIAL CASE analagous to ROOT_IL*/
	  if(doing_search)
	    tmp_imin = nss_imin[n+1]; /* same tmp_imin as ROOT_S */
	  else
	    tmp_imin = i0; /* we never emit to the left in this state */
	  tmp_imax = nss_imax[n+1]; 
	  tmp_jmin = nss_jmin[n+1];
	  tmp_jmax = j0; /* Have to be able to transit here from ROOT_S */
	  hmm2ij_insert_state_step1_set_state_bands(v, tmp_imin, tmp_imax, tmp_jmin,
						    tmp_jmax, imin, imax, jmin, jmax);
	  /* Enforce safe transitions, this makes sure that at least one state
	   * y \in C_v is reachable from v. And further (special case for inserts)
	   * make sure that we don't consider v as a possible y.  IF we did, we might
	   * be faced with a situation where v could only transit to itself, and then
	   * we'd be in the same situation, but we may no longer be able to transit ANYWHERE
	   * including to itself.
	   */
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, (nss_max_imin[n+1]), nss_min_jmax[n+1]);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  break;
	  
	case BEGL_nd:
	  hmm2ij_prestate_step0_initialize(n, nss_max_imin, nss_min_jmax, i0, j0);
	  hmm2ij_prestate_step1_set_node_inserts(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						 nss_imin, nss_imax, nss_jmin, nss_jmax,
						 pn_min_i, pn_max_i, cp9map);

	  hmm2ij_prestate_step2_determine_safe(n, nss_max_imin[n+1], nss_min_jmax[n+1],
					       nis_imin[n], nis_jmax[n],
					       &safe_imax, &safe_jmin);
	  hmm2ij_prestate_step3_preset_node_splits(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax,
						   pn_min_m, pn_max_m, pn_min_d, pn_max_d,
						   cp9map);
	  /* 1 state BEGL_S */
	  v = cm->nodemap[n];
	  /* The next node MUST be a match node (MATP 
	   * specifically due to model building
	   * algorithm) or a BIF node. We derive imin, imax,
	   * jmin and jmax from that node.
	   */
	  /* Use the next nodes split set band, which
	   * will be wider of match and delete states bands
	   * for split set states in next node. 
	   */
	  tmp_imin = nss_imin[n+1];
	  tmp_imax = nss_imax[n+1];
	  tmp_jmin = nss_jmin[n+1];
	  tmp_jmax = nss_jmax[n+1];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  hmm2ij_state_step5_non_emitter_d0_hack(v, imax[v], jmin);
	  break;

	case BEGR_nd:
	  hmm2ij_prestate_step0_initialize(n, nss_max_imin, nss_min_jmax, i0, j0);
	  hmm2ij_prestate_step1_set_node_inserts(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						 nss_imin, nss_imax, nss_jmin, nss_jmax,
						 pn_min_i, pn_max_i, cp9map);

	  hmm2ij_prestate_step2_determine_safe(n, nss_max_imin[n+1], nss_min_jmax[n+1],
					       nis_imin[n], nis_jmax[n],
					       &safe_imax, &safe_jmin);
	  hmm2ij_prestate_step3_preset_node_splits(n, nis_imin, nis_imax, nis_jmin, nis_jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax,
						   pn_min_m, pn_max_m, pn_min_d, pn_max_d,
						   cp9map);
	  /* 2 states BEGR_S and BEGR_IL */
	  v = cm->nodemap[n]; /*BEGR_S*/	  
	  /* Use either the next nodes split set band, which
	   * will be wider of match and delete states bands
	   * for split set states in next node OR
	   * the band on the insert state that maps to the
	   * BEGR_IL, erring on the safe side (wider band).
	   */
	  tmp_imin = nss_imin[n+1];
	  tmp_imax = nss_imax[n+1];
	  if(pn_min_i[cp9map->cs2hn[v+1][0]] < tmp_imin)
	    tmp_imin = pn_min_i[cp9map->cs2hn[v+1][0]];
	  if(pn_max_i[cp9map->cs2hn[v+1][0]] > tmp_imax)
	    tmp_imax = pn_max_i[cp9map->cs2hn[v+1][0]];
	  tmp_jmin = nss_jmin[n+1];
	  tmp_jmax = nss_jmax[n+1];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  hmm2ij_state_step5_non_emitter_d0_hack(v, imax[v], jmin);

	  v++; /*BEGR_IL*/
	  /* This state maps to the insert state of HMM node cp9map->cs2hn[v][0]*/
	  tmp_imin = pn_min_i[cp9map->cs2hn[v][0]];
	  tmp_imax = pn_max_i[cp9map->cs2hn[v][0]];
	  tmp_jmin = nss_jmin[n+1];
	  tmp_jmax = nss_jmax[n+1];
	  hmm2ij_insert_state_step1_set_state_bands(v, tmp_imin, tmp_imax, tmp_jmin,
						    tmp_jmax, imin, imax, jmin, jmax);
	  /* Enforce safe transitions, this makes sure that at least one state
	   * y \in C_v is reachable from v. And further (special case for inserts)
	   * make sure that we don't consider v as a possible y.  IF we did, we might
	   * be faced with a situation where v could only transit to itself, and then
	   * we'd be in the same situation, but we may no longer be able to transit ANYWHERE
	   * including to itself.
	   */
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, (nss_max_imin[n+1]), nss_min_jmax[n+1]);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  break;

	case BIF_nd:
	  hmm2ij_prestate_step0_initialize(n, nss_max_imin, nss_min_jmax, i0, j0);

	  /* 1 state BIF_B */
	  v = cm->nodemap[n]; /*BIF_B*/
	  /* The only two connected states are BEGL_S and BEGR_S.
	   * We can derive our imin, imax, jmin, and jmax from 
	   * those two states.
	   * cm->cfirst[v] is the state index of the left child.
	   * cm->cnum[v] is the state index of the right child.
	   */
	  nis_imin[n] = imin[cm->cfirst[v]];
	  nis_imax[n] = imax[cm->cfirst[v]];
	  nis_jmin[n] = jmin[cm->cnum[v]];
	  nis_jmax[n] = jmax[cm->cnum[v]];

	  nss_imin[n] = imin[cm->cfirst[v]];
	  nss_imax[n] = imax[cm->cfirst[v]];
	  nss_jmin[n] = jmin[cm->cnum[v]];
	  nss_jmax[n] = jmax[cm->cnum[v]];

	  hmm2ij_prestate_step2_determine_safe(n, nss_max_imin[n+1], nss_min_jmax[n+1],
					       nis_imin[n], nis_jmax[n],
					       &safe_imax, &safe_jmin);
	  tmp_imin = imin[cm->cfirst[v]];
	  tmp_imax = imax[cm->cfirst[v]];
	  tmp_jmin = jmin[cm->cnum[v]];
	  tmp_jmax = jmax[cm->cnum[v]];
	  hmm2ij_split_state_step1_set_state_bands(v, n, tmp_imin, tmp_imax, tmp_jmin,
						   tmp_jmax, imin, imax, jmin, jmax,
						   nss_imin, nss_imax, nss_jmin, nss_jmax);
	  hmm2ij_state_step2_enforce_safe_trans(cm, v, n, imax, jmin, nss_imax, nss_jmin, safe_imax, 
	  					safe_jmin);
	  hmm2ij_state_step3_enforce_state_delta(cm, v, jmin, jmax);
	  hmm2ij_state_step4_update_safe_holders(v, n, imin[v], jmax[v], nss_max_imin, nss_min_jmax);
	  hmm2ij_state_step5_non_emitter_d0_hack(v, imax[v], jmin);
	  break;
	}
    }

   /* Tie up some loose ends: 
    * 1. Set detached inserts states to imin=imax=jmin=jmax=i0 to avoid 
    *    problems in downstream functions. These states WILL NEVER BE ENTERED 
    * 2. Do a quick check to make sure we've assigned the bands
    *    on i and j for all states to positive values (none were
    *    left as -1 EXCEPT for end states which should have i bands left as -1).
    * 3. Ensure that all *max[v] and *min[v] values are <= L, values greater
    *    than this don't make sense.
    */

  for(v = 0; v < cm->M; v++) {
    /* set bands for detached inserts */
    if(cm->sttype[v+1] == E_st) imin[v] = imax[v] = jmin[v] = jmax[v] = i0;

    /* Ensure: for all i imin[v]..i..imax[v]
     *             i0 <= i <= j0+1
     *         for all j jmin[v]..j..jmax[v]
     *             i0 <= j <= j0
     * Note: i can be j0+1 to allow delete states to be entered with 
     * d = 0, after the entire seq has been emitted.
     */
    imin[v] = ESL_MAX(imin[v], i0);
    imin[v] = ESL_MIN(imin[v], j0+1);
    imax[v] = ESL_MAX(imax[v], i0);
    imax[v] = ESL_MIN(imax[v], j0+1);
    jmin[v] = ESL_MAX(jmin[v], i0);
    jmin[v] = ESL_MIN(jmin[v], j0);
    jmax[v] = ESL_MAX(jmax[v], i0);
    jmax[v] = ESL_MIN(jmax[v], j0);

    /* Ensure: jmax[v] - jmin[v] + 1 >= 0 
     *         imax[v] - imin[v] + 1 >= 0 
     * jmax[v] - jmin[v] + 1 == 0 means there are no valid j's for state v,
     * so state v is not allowed to be in the parse, we allow this (maybe we shouldn't)
     */
    imin[v] = ESL_MIN(imin[v], imax[v]+1);
    jmin[v] = ESL_MIN(jmin[v], jmax[v]+1);

    ESL_DASSERT1((! ((cm->sttype[v] != E_st) && (imin[v] == -1))));
    ESL_DASSERT1((! ((cm->sttype[v] != E_st) && (imax[v] == -1))));
    ESL_DASSERT1((! ((cm->sttype[v] != E_st) && (jmin[v] == -1))));
    ESL_DASSERT1((! ((cm->sttype[v] != E_st) && (jmax[v] == -1))));
  }

  if(debug_level > 0) { 
    printf("bands on i\n");
    debug_print_bands(stdout, cm, imin, imax);
    
    printf("bands on j\n");
    debug_print_bands(stdout, cm, jmin, jmax);
  }
  free(nss_imin);
  free(nss_imax);
  free(nss_jmin);
  free(nss_jmax);
  free(nis_imin);
  free(nis_imax);
  free(nis_jmin);
  free(nis_jmax); 
  free(nss_max_imin);
  free(nss_min_jmax);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}

/**************************************************************************
 * Helper functions for *_cp9_HMM2ijBands() 
 *  hmm2ij_prestate_step0_initialize()
 *  hmm2ij_prestate_step1_set_node_inserts()
 *  hmm2ij_prestate_step2_determine_safe()
 *  hmm2ij_prestate_step3_preset_node_splits()
 *  hmm2ij_split_state_step1_set_state_bands()
 *  hmm2ij_insert_state_step1_set_state_bands()
 *  hmm2ij_state_step2_enforce_safe_trans()
 *  hmm2ij_state_step5_non_emitter_d0_hack()
 */

/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_prestate_step0_initialize
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void 
hmm2ij_prestate_step0_initialize(int n, int *nss_max_imin, int *nss_min_jmax, int i0, int j0)
{
  nss_max_imin[n] = i0-1;
  nss_min_jmax[n] = j0;
}

/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_prestate_step1_set_node_inserts
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void 
hmm2ij_prestate_step1_set_node_inserts(int n, int *nis_imin, int *nis_imax, 
				       int *nis_jmin, int *nis_jmax,
				       int *nss_imin, int *nss_imax, 
				       int *nss_jmin, int *nss_jmax,
				       int *pn_min_i, int *pn_max_i, 
				       CP9Map_t *cp9map)

{
  if(cp9map->nd2lpos[n] != -1)
    {  
      nis_imin[n] = pn_min_i[cp9map->nd2lpos[n]];
      nis_imax[n] = pn_max_i[cp9map->nd2lpos[n]];
    }
  else
    {
      nis_imin[n] = nss_imin[n+1];
      nis_imax[n] = nss_imax[n+1];
    }
  if(cp9map->nd2rpos[n] != -1)
    {  
      nis_jmin[n] = pn_min_i[cp9map->nd2rpos[n]];
      nis_jmax[n] = pn_max_i[cp9map->nd2rpos[n]];
    }
  else
    {
      nis_jmin[n] = nss_jmin[n+1];
      nis_jmax[n] = nss_jmax[n+1];
    }
}
/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_prestate_step1_set_node_inserts
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void 
hmm2ij_prestate_step2_determine_safe(int n, 	
				     int nss_max_imin_np1, int nss_min_jmax_np1,
				     int nis_imin_n, 
				     int nis_jmax_n,
				     int *safe_imax, int *safe_jmin)
{
  *safe_imax = (nss_max_imin_np1 < nis_imin_n) ? 
    nss_max_imin_np1 : nis_imin_n;
  *safe_jmin = (nss_min_jmax_np1 > nis_jmax_n) ? 
    nss_min_jmax_np1 : nis_jmax_n;
}
/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_prestate_step1_set_node_inserts
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void 
hmm2ij_prestate_step3_preset_node_splits(int n, int *nis_imin, int *nis_imax, 
					 int *nis_jmin, int *nis_jmax,
					 int *nss_imin, int *nss_imax, 
					 int *nss_jmin, int *nss_jmax,
					 int *pn_min_m, int *pn_max_m, 
					 int *pn_min_d, int *pn_max_d, 
					 CP9Map_t *cp9map)
{
  if(cp9map->nd2lpos[n] != -1)
    {  
      nss_imin[n] = (pn_min_m[cp9map->nd2lpos[n]] < (pn_min_d[cp9map->nd2lpos[n]])) ?
	pn_min_m[cp9map->nd2lpos[n]] : (pn_min_d[cp9map->nd2lpos[n]]);
      nss_imax[n] = (pn_max_m[cp9map->nd2lpos[n]] > (pn_max_d[cp9map->nd2lpos[n]])) ?
	pn_max_m[cp9map->nd2lpos[n]] : (pn_max_d[cp9map->nd2lpos[n]]);
    }
  else
    {
      nss_imin[n] = nss_imin[n+1];
      nss_imax[n] = nss_imax[n+1];
    }
  if(cp9map->nd2rpos[n] != -1)
    {  
      nss_jmin[n] = (pn_min_m[cp9map->nd2rpos[n]] < pn_min_d[cp9map->nd2rpos[n]]) ?
	pn_min_m[cp9map->nd2rpos[n]] : pn_min_d[cp9map->nd2rpos[n]];
      nss_jmax[n] = (pn_max_m[cp9map->nd2rpos[n]] > pn_max_d[cp9map->nd2rpos[n]]) ?
	pn_max_m[cp9map->nd2rpos[n]] : pn_max_d[cp9map->nd2rpos[n]];
    }
  else
    {
      nss_jmin[n] = nss_jmin[n+1];
      nss_jmax[n] = nss_jmax[n+1];
    }
}

/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_split_state_step1_set_state_bands
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void 
hmm2ij_split_state_step1_set_state_bands(int v, int n, 
					 int tmp_imin, int tmp_imax, 
					 int tmp_jmin, int tmp_jmax,
					 int *imin, int *imax, int *jmin, int *jmax,
					 int *nss_imin, int *nss_imax,
					 int *nss_jmin, int *nss_jmax)
{
  imin[v] = tmp_imin;
  imax[v] = tmp_imax;
  jmin[v] = tmp_jmin;
  jmax[v] = tmp_jmax;
  if(imin[v] < nss_imin[n])
    nss_imin[n] = imin[v];
  if(imax[v] > nss_imax[n])
    nss_imax[n] = imax[v];
  if(jmin[v] < nss_jmin[n])
    nss_jmin[n] = jmin[v];
  if(jmax[v] > nss_jmax[n])
    nss_jmax[n] = jmax[v];

}
/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_prestate_step1_set_node_inserts
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void hmm2ij_insert_state_step1_set_state_bands(int v,
					       int tmp_imin, int tmp_imax, 
					       int tmp_jmin, int tmp_jmax,
					       int *imin, int *imax, int *jmin, int *jmax)
{
  imin[v] = tmp_imin;
  imax[v] = tmp_imax;
  jmin[v] = tmp_jmin;
  jmax[v] = tmp_jmax;
}
/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_state_step2_enforce_safe_trans
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void 
hmm2ij_state_step2_enforce_safe_trans(CM_t *cm, int v, int n, int *imax, int *jmin,
				      int *nss_imax, int *nss_jmin, 
				      int safe_imax, int safe_jmin)
{
  int dv_l;
  int dv_r;
  if((cm->sttype[v] == ML_st) ||
     (cm->sttype[v] == IL_st) ||
     (cm->sttype[v] == MP_st))     
    dv_l = 1;
  else
    dv_l = 0;
  if((cm->sttype[v] == MR_st) ||
     (cm->sttype[v] == IR_st) ||
     (cm->sttype[v] == MP_st))     
    dv_r = 1;
  else
    dv_r = 0;
  if(imax[v] < safe_imax - dv_l)
    {
      imax[v] = safe_imax - dv_l;
      if(imax[v] > nss_imax[n]) 
	nss_imax[n] = imax[v];
    }
  if(jmin[v] > safe_jmin + dv_r)
    {
      jmin[v] = safe_jmin + dv_r;
      if(jmin[v] < nss_jmin[n])
	nss_jmin[n] = jmin[v];
    }
}

/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_state_step3_enforce_state_delta
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void 
hmm2ij_state_step3_enforce_state_delta(CM_t *cm, int v, int *jmin, int *jmax)
{
  int dv_l;
  int dv_r;
  if((cm->sttype[v] == ML_st) ||
     (cm->sttype[v] == IL_st) ||
     (cm->sttype[v] == MP_st))     
    dv_l = 1;
  else
    dv_l = 0;
  if((cm->sttype[v] == MR_st) ||
     (cm->sttype[v] == IR_st) ||
     (cm->sttype[v] == MP_st))     
    dv_r = 1;
  else
    dv_r = 0;
  if(jmin[v] < (dv_l + dv_r))
     jmin[v] = dv_l + dv_r;
  if(jmax[v] < (dv_l + dv_r))
    jmax[v] = dv_l + dv_r;
}
/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_state_step4_update_safe_holders
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void
hmm2ij_state_step4_update_safe_holders(int v, int n, int imin_v, int jmax_v, int *nss_max_imin, 
				       int *nss_min_jmax)
{
  if(imin_v > nss_max_imin[n])
    nss_max_imin[n] = imin_v;
  if(jmax_v < nss_min_jmax[n])
    nss_min_jmax[n] = jmax_v;
}

/*****************************************************************************
 * EPN 12.21.05
 * Function: hmm2ij_state_step5_non_emitter_d0_hack
 *
 * Purpose:  cp9_HMM2ijBands*() function helper function. 
 * 
 *****************************************************************************/
void
hmm2ij_state_step5_non_emitter_d0_hack(int v, int imax_v, int *jmin)
{
  /* allow for possibility that d=0 for delete states*/
  if(jmin[v] <= imax_v && jmin[v] > 0)
    jmin[v]--;
  /* if imax = L, allow possibility for 
  if(imax[v] == Limax_v && jmin[v] > 0)
  jmin[v]--;*/
}

/****************************************************************************
 * Debugging print functions
 * cp9_DebugPrintHMMBands()
 * ijBandedTraceInfoDump()
 * ijdBandedTraceInfoDump()
 * cp9_DebugCheckFB()
 */
/**************************************************************
 * EPN 12.18.05
 * cp9_DebugPrintHMMBands()
 * based loosely on: cmbuild.c's
 * Function: model_trace_info_dump
 *
 * Purpose:  Print out the bands derived from the posteriors for the
 *           insert and match states of each HMM node.
 * 
 * Args:    
 * FILE *ofp      - filehandle to print to (can by STDOUT)
 * int L          - length of sequence
 * CP9Bands_t     - the CP9 bands data structure
 * double hmm_bandp - fraction of probability mass allowed outside each band.
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 * Returns: (void) 
 */

void
cp9_DebugPrintHMMBands(FILE *ofp, int L, CP9Bands_t *cp9b, double hmm_bandp, int debug_level)
{
  int M;
  int k;
  int cells_in_bands_m; /* number of cells within all the bands for match states*/
  int cells_in_bands_i; /* number of cells within all the bands for insert states*/
  int cells_in_bands_d; /* number of cells within all the bands for delete states*/
  int cells_in_bands_all; /* number of cells within all the bands for match and insert states*/
  /*double fraction_excluded_m;*/ /* fraction of the dp matrix excluded from bands for all match states*/
  /*double fraction_excluded_i;*/ /* fraction of the dp matrix excluded from bands for all insert states*/
  /*double fraction_excluded_d;*/ /* fraction of the dp matrix excluded from bands for all delete states*/
  /*double fraction_excluded_all;*/ /* fraction of the dp matrix excluded from bands for all match and insert states*/

  M = cp9b->hmm_M;
  cells_in_bands_m = cells_in_bands_i = cells_in_bands_d = cells_in_bands_all = 0;

  /* first print the bands on the match states */
  fprintf(ofp, "***********************************************************\n");
  if(debug_level > 0)
    fprintf(ofp, "printing hmm bands\n");
  fprintf(ofp, "hmm_bandp: %f\n", hmm_bandp);
  if(debug_level > 0)
    {    
      fprintf(ofp, "\n");
      fprintf(ofp, "match states\n");
    }
  for(k = 0; k <= cp9b->hmm_M; k++)
    {
      if(debug_level > 0 || debug_level == -1)
	fprintf(ofp, "M node: %3d | min %3d | max %3d | w %3d \n", k, cp9b->pn_min_m[k], cp9b->pn_max_m[k], (cp9b->pn_max_m[k] - cp9b->pn_min_m[k]+1));
      cells_in_bands_m += cp9b->pn_max_m[k] - cp9b->pn_min_m[k] + 1;
    }
  if(debug_level > 0)
    fprintf(ofp, "\n");
  if(debug_level > 0)
    fprintf(ofp, "insert states\n");
  for(k = 0; k <= cp9b->hmm_M; k++)
    {
      if(debug_level > 0 || debug_level == -1)
	fprintf(ofp, "I node: %3d | min %3d | max %3d | w %3d\n", k, cp9b->pn_min_i[k], cp9b->pn_max_i[k], (cp9b->pn_max_i[k] - cp9b->pn_min_i[k]+1));
      cells_in_bands_i += cp9b->pn_max_i[k] - cp9b->pn_min_i[k] + 1;
    }
  if(debug_level > 0)
    fprintf(ofp, "\n");
  if(debug_level > 0)
    fprintf(ofp, "delete states\n");
  for(k = 1; k <= cp9b->hmm_M; k++)
    {
      if(debug_level > 0 || debug_level == -1)
	fprintf(ofp, "D node: %3d | min %3d | max %3d | w %3d\n", k, cp9b->pn_min_d[k], cp9b->pn_max_d[k], (cp9b->pn_max_d[k] - cp9b->pn_min_d[k]+1));
      cells_in_bands_d += cp9b->pn_max_d[k] - cp9b->pn_min_d[k] + 1;
    }
  if(debug_level > 0)
    {
      fprintf(ofp, "\n");
      printf("cells_in_bands_m : %d\n", cells_in_bands_m);
      printf("cells_in_bands_i : %d\n", cells_in_bands_i);
      printf("cells_in_bands_d : %d\n", cells_in_bands_d);
    }

  cells_in_bands_all = cells_in_bands_m + cells_in_bands_i + cells_in_bands_d;
  printf("fraction match excluded  : %f\n", (1 - ((float) cells_in_bands_m / (M * L))));
  printf("fraction insert excluded : %f\n", (1 - ((float) cells_in_bands_i / ((M-1) * L))));
  printf("fraction delete excluded : %f\n", (1 - ((float) cells_in_bands_d / ((M-1) * L))));
  printf("fraction total excluded  : %f\n", (1 - ((float) (cells_in_bands_all) / (((M-1) * L) + ((M-1) * L) + (M *L)))));
  fprintf(ofp, "***********************************************************\n");
	 
}

/**************************************************************
 * cp9_compare_bands()
 * based loosely on: cmbuild.c's
 * Function: model_trace_info_dump
 *
 * Purpose:  Compare 2 CP9Bands_t objects, die if they're at all different.
 * 
 * Returns: (void) 
 */
void
cp9_compare_bands(CP9Bands_t *cp9b1, CP9Bands_t *cp9b2)
{
  int k, v, d;

  if(cp9b1->hmm_M != cp9b2->hmm_M) cm_Fail("cp9_compare_bands(): cp9b1->hmm_M: %d != cp9b2->hmm_M: %d\n", cp9b1->hmm_M, cp9b2->hmm_M);
  if(cp9b1->cm_M  != cp9b2->cm_M)  cm_Fail("cp9_compare_bands(): cp9b1->cm_M: %d != cp9b2->cm_M: %d\n", cp9b1->cm_M, cp9b2->cm_M);
  for(k = 0; k <= cp9b1->hmm_M; k++)
    {
      if(cp9b1->pn_min_m[k] != cp9b2->pn_min_m[k]) cm_Fail("cp9_compare_bands(): cp9b1->pn_min_m[%d]: %d != cp9b2->pn_min_m[%d]: %d\n", k, cp9b1->pn_min_m[k], k, cp9b2->pn_min_m[k]);
      if(cp9b1->pn_min_i[k] != cp9b2->pn_min_i[k]) cm_Fail("cp9_compare_bands(): cp9b1->pn_min_i[%d]: %d != cp9b2->pn_min_i[%d]: %d\n", k, cp9b1->pn_min_i[k], k, cp9b2->pn_min_i[k]);
      if(cp9b1->pn_min_d[k] != cp9b2->pn_min_d[k]) cm_Fail("cp9_compare_bands(): cp9b1->pn_min_d[%d]: %d != cp9b2->pn_min_d[%d]: %d\n", k, cp9b1->pn_min_d[k], k, cp9b2->pn_min_d[k]);

      if(cp9b1->pn_max_m[k] != cp9b2->pn_max_m[k]) cm_Fail("cp9_compare_bands(): cp9b1->pn_max_m[%d]: %d != cp9b2->pn_max_m[%d]: %d\n", k, cp9b1->pn_max_m[k], k, cp9b2->pn_max_m[k]);
      if(cp9b1->pn_max_i[k] != cp9b2->pn_max_i[k]) cm_Fail("cp9_compare_bands(): cp9b1->pn_max_i[%d]: %d != cp9b2->pn_max_i[%d]: %d\n", k, cp9b1->pn_max_i[k], k, cp9b2->pn_max_i[k]);
      if(cp9b1->pn_max_d[k] != cp9b2->pn_max_d[k]) cm_Fail("cp9_compare_bands(): cp9b1->pn_max_d[%d]: %d != cp9b2->pn_max_d[%d]: %d\n", k, cp9b1->pn_max_d[k], k, cp9b2->pn_max_d[k]);
    }
  for(v = 0; v < cp9b1->cm_M; v++)
    {
      if(cp9b1->imin[v] != cp9b2->imin[v]) cm_Fail("cp9_compare_bands(): cp9b1->imin[%d]: %d != cp9b2->imin[%d]: %d\n", v, cp9b1->imin[v], v, cp9b2->imin[v]);
      if(cp9b1->imax[v] != cp9b2->imax[v]) cm_Fail("cp9_compare_bands(): cp9b1->imax[%d]: %d != cp9b2->imax[%d]: %d\n", v, cp9b1->imax[v], v, cp9b2->imax[v]);

      if(cp9b1->jmin[v] != cp9b2->jmin[v]) cm_Fail("cp9_compare_bands(): cp9b1->jmin[%d]: %d != cp9b2->jmin[%d]: %d\n", v, cp9b1->jmin[v], v, cp9b2->jmin[v]);
      if(cp9b1->jmax[v] != cp9b2->jmax[v]) cm_Fail("cp9_compare_bands(): cp9b1->jmax[%d]: %d != cp9b2->jmax[%d]: %d\n", v, cp9b1->jmax[v], v, cp9b2->jmax[v]);

      for(d = 0; d <= (cp9b1->jmax[v] - cp9b1->jmin[v]); d++) {
	if(cp9b1->hdmin[v][d] != cp9b2->hdmin[v][d]) cm_Fail("cp9_compare_bands(): cp9b1->hdmin[%d][%d]: %d != cp9b2->hdmin[%d][%d]: %d\n", v, d, cp9b1->hdmin[v][d], v, d, cp9b2->hdmin[v][d]);
	if(cp9b1->hdmax[v][d] != cp9b2->hdmax[v][d]) cm_Fail("cp9_compare_bands(): cp9b1->hdmax[%d][%d]: %d != cp9b2->hdmax[%d][%d]: %d\n", v, d, cp9b1->hdmax[v][d], v, d, cp9b2->hdmax[v][d]);
      }
      /* don't compare safe_hdmin, safe_hdmax, b/c we want to be able to compare cp9bands objects
       * before safe_hdmin, safe_hdmax are calced
       */
    }
}      

/* EPN 11.03.05
 * Function: ijBandedTraceInfoDump()
 *
 * Purpose:  Experimental HMMERNAL function used in development.
 *           This function determines how close the
 *           trace was to the bands for i and j at each state in the trace,
 *           and prints out that information in differing levels
 *           of verbosity depending on an input parameter 
 *           (debug_level).
 * 
 * Args:    cm       - the CM (useful for determining which states are E states)
 *          tr       - the parsetree (trace)
 *          imin     - minimum i bound for each state v; [0..v..M-1]
 *          imax     - maximum i bound for each state v; [0..v..M-1]
 *          jmin     - minimum j bound for each state v; [0..v..M-1]
 *          jmax     - maximum j bound for each state v; [0..v..M-1]
 *          debug_level - level of verbosity
 * Returns: (void) 
 */

void
ijBandedTraceInfoDump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
		      int *jmin, int *jmax, int debug_level)
{
  int status;
  int v, i, j, d, tpos;
  int imindiff;            /* i - imin[v] */
  int imaxdiff;            /* imax[v] - i */
  int jmindiff;            /* j - jmin[v] */
  int jmaxdiff;            /* jmax[v] - j */
  int imin_out;
  int imax_out;
  int jmin_out;
  int jmax_out;
  char **sttypes;
  char **nodetypes;

  ESL_ALLOC(sttypes, sizeof(char *) * 10);
  sttypes[0] = "D";
  sttypes[1] = "MP";
  sttypes[2] = "ML";
  sttypes[3] = "MR";
  sttypes[4] = "IL";
  sttypes[5] = "IR";
  sttypes[6] = "S";
  sttypes[7] = "E";
  sttypes[8] = "B";
  sttypes[9] = "EL";

  ESL_ALLOC(nodetypes, sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  imin_out = 0;
  imax_out = 0;
  jmin_out = 0;
  jmax_out = 0;

  debug_level = 2;

  for (tpos = 0; tpos < tr->n; tpos++)
    {
      v  = tr->state[tpos];
      i = tr->emitl[tpos];
      j = tr->emitr[tpos];
      d = j-i+1;
      imindiff = i-imin[v];
      imaxdiff = imax[v]-i;
      jmindiff = j-jmin[v];
      jmaxdiff = jmax[v]-j;
      if(cm->sttype[v] != E_st)
	{
	  if(imindiff < 0)
	    imin_out++;
	  if(imaxdiff < 0)
	    imax_out++;
	  if(jmindiff < 0)
	    jmin_out++;
	  if(jmaxdiff < 0)
	    jmax_out++;
	  
	  if(debug_level > 1 || ((imindiff < 0) || (imaxdiff < 0) || (jmindiff < 0) || (jmaxdiff < 0)))
	    {
	      printf("v: %4d %-4s %-2s | d: %4d | i: %4d | in: %4d | ix: %4d | %3d | %3d |\n", v, nodetypes[(int) cm->ndtype[(int) cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], d, i, imin[v], imax[v], imindiff, imaxdiff);
	      printf("                          | j: %4d | jn: %4d | jx: %4d | %3d | %3d |\n", j, jmin[v], jmax[v], jmindiff, jmaxdiff);

	    }
	}
      else if(cm->sttype[v] == E_st)
	{
	  if(debug_level > 1)
	    {
	      printf("v: %4d %-4s %-2s | d: %4d | i: %4d | in: %4d | ix: %4d | %3d | %3d |\n", v, nodetypes[(int) cm->ndtype[(int) cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], d, i, imin[v], imax[v], imindiff, imaxdiff);
	      printf("                          | j: %4d | jn: %4d | jx: %4d | %3d | %3d |\n", j, jmin[v], jmax[v], jmindiff, jmaxdiff);
	    }
	}
    }
  printf("\nimin out: %d\n", imin_out);
  printf("imax out: %d\n", imax_out);
  printf("jmin out: %d\n", jmin_out);
  printf("jmax out: %d\n", jmax_out);
  
  if((imin_out + imax_out + jmin_out + jmax_out) > 0)
    {
      printf("ERROR, some of the i and j bands are going to prevent optimal alignment. Sorry.\n");
    }

  free(sttypes);
  free(nodetypes);
  return;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
}


/* EPN 11.03.05
 * Function: ijdBandedTraceInfoDump()
 *
 * Purpose:  Experimental HMMERNAL function used in development.
 *           This function determines how close the
 *           trace was to the bands for i and j and d at each state in the trace,
 *           and prints out that information in differing levels
 *           of verbosity depending on an input parameter 
 *           (debug_level).
 * 
 * Args:    cm       - the CM (useful for determining which states are E states)
 *          tr       - the parsetree (trace)
 *          imin     - minimum i bound for each state v; [0..v..M-1]
 *          imax     - maximum i bound for each state v; [0..v..M-1]
 *          jmin     - minimum j bound for each state v; [0..v..M-1]
 *          jmax     - maximum j bound for each state v; [0..v..M-1]
 *          hdmin    - minimum d bound for each state v and offset j;
 *                     [0..v..M-1][0..(jmax[v]-jmin[v])]
 *          hdmax    - maximum d bound for each state v and offset j; 
 *                     [0..v..M-1][0..(jmax[v]-jmin[v])]
 *          debug_level - level of verbosity
 * Returns: (void) 
 */

void
ijdBandedTraceInfoDump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
		       int *jmin, int *jmax, int **hdmin, int **hdmax, int debug_level)
{
  int status;
  int v, i, j, d, tpos;
  int imindiff;            /* i - imin[v] */
  int imaxdiff;            /* imax[v] - i */
  int jmindiff;            /* j - jmin[v] */
  int jmaxdiff;            /* jmax[v] - j */
  int hdmindiff;           /* d - hdmin[v][j] */
  int hdmaxdiff;           /* hdmax[v][j] - d */

  int imin_out;
  int imax_out;
  int jmin_out;
  int jmax_out;
  int hdmin_out;
  int hdmax_out;
  int local_used;

  char **sttypes;
  char **nodetypes;

  ESL_ALLOC(sttypes, sizeof(char *) * 10);
  sttypes[0] = "D";
  sttypes[1] = "MP";
  sttypes[2] = "ML";
  sttypes[3] = "MR";
  sttypes[4] = "IL";
  sttypes[5] = "IR";
  sttypes[6] = "S";
  sttypes[7] = "E";
  sttypes[8] = "B";
  sttypes[9] = "EL";

  ESL_ALLOC(nodetypes, sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  imin_out = 0;
  imax_out = 0;
  jmin_out = 0;
  jmax_out = 0;
  hdmin_out = 0;
  hdmax_out = 0;
  local_used = 0;

  debug_level = 2;

  for (tpos = 0; tpos < tr->n; tpos++)
    {
      v  = tr->state[tpos];
      i = tr->emitl[tpos];
      j = tr->emitr[tpos];
      d = j-i+1;
      if(cm->sttype[v] == EL_st) /*END LOCAL state*/
	{
	  if(debug_level > 1)
	    {
	      printf("v: %4d NA   %-2s (  NA) | d: %4d | i: %4d | in: NA    | ix: NA   | NA  | NA  |\n", v, sttypes[(int) cm->sttype[v]], d, i);
	  printf("                                 | j: %4d | jn: NA   | jx: NA  | NA  | NA  |\n", j);
	  printf("                                 | d: %4d | dn: NA   | dx: NA   | NA  | NA  |\n", d);
	  
	  local_used++;
	    }
	}
      else
	{
	  imindiff = i-imin[v];
	  imaxdiff = imax[v]-i;
	  jmindiff = j-jmin[v];
	  jmaxdiff = jmax[v]-j;
	  if(j >= jmin[v] && j <= jmax[v])
	    {
	      hdmindiff = d - hdmin[v][j-jmin[v]];
	      hdmaxdiff = hdmax[v][j-jmin[v]] - d;
	    }  
	  else
	    {
	      hdmindiff = -1000;
	      hdmaxdiff = -1000;
	    }
	  if(imindiff < 0)
	    imin_out++;
	  if(imaxdiff < 0)
	    imax_out++;
	  if(jmindiff < 0)
	    jmin_out++;
	  if(jmaxdiff < 0)
	    jmax_out++;
	  if(hdmindiff < 0)
	    hdmin_out++;
	  if(hdmaxdiff < 0)
	    hdmax_out++;
	  
	  if(debug_level > 1 || ((imindiff < 0) || (imaxdiff < 0) || (jmindiff < 0) || (jmaxdiff < 0) || 
				 (hdmindiff < 0) || (hdmaxdiff < 0)))
	    {
	      printf("v: %4d %-4s %-2s (%4d) | d: %4d | i: %4d | in: %4d | ix: %4d | %3d | %3d |\n", v, nodetypes[(int) cm->ndtype[(int) cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], cm->ndidx[v], d, i, imin[v], imax[v], imindiff, imaxdiff);
	      printf("                                 | j: %4d | jn: %4d | jx: %4d | %3d | %3d |\n", j, jmin[v], jmax[v], jmindiff, jmaxdiff);
	      if(j >= jmin[v] && j <= jmax[v])
		{
		  printf("                                 | d: %4d | dn: %4d | dx: %4d | %3d | %3d |\n", d, hdmin[v][j-jmin[v]], hdmax[v][j-jmin[v]], hdmindiff, hdmaxdiff);
		}	  
	      else
		{
		  printf("                                 | d: %4d | dn: jout | dx: jout | %3d | %3d |\n", d, hdmindiff, hdmaxdiff);
		}
	    }
	}
    }
  printf("\nimin out  : %d\n", imin_out);
  printf("imax out  : %d\n", imax_out);
  printf("jmin out  : %d\n", jmin_out);
  printf("jmax out  : %d\n", jmax_out);
  printf("hdmin out : %d\n", hdmin_out);
  printf("hdmax out : %d\n", hdmax_out);
  printf("local used: %d\n", local_used);
  
  if((imin_out + imax_out + jmin_out + jmax_out) > 0)
    {
      printf("ERROR, some of the i and j bands are going to prevent optimal alignment. Sorry.\n");
    }

  free(sttypes);
  free(nodetypes);
  return;

 ERROR:
  cm_Fail("Memory allocation error.\n");
}


/*********************************************************************
 * Function: cp9_DebugCheckFB()
 * 
 * Purpose:  Debugging function to make sure CP9Forward() and 
 *           CP9Backward are working by checking:
 *           For all positions i, and states k:
 *             sum_k f[i][k] * b[i][k] = P(x|hmm)
 *           
 * Args:     fmx    - forward dp matrix, already filled
 *           bmx    - backward dp matrix, already filled
 *           hmm    - the model
 *           sc     - P(x|hmm) the probability of the entire
 *                    seq given the model
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           dsq    - the digitized sequence
 *           
 * Note about sequence position indexing: although this function
 * works on a subsequence from i0 to j0, fmx and bmx have offset indices,
 * from 1 to W, with W = j0-i0+1.
 * Return:   (void) Exits if any errors are found.
 */
void
cp9_DebugCheckFB(struct cp9_dpmatrix_s *fmx, struct cp9_dpmatrix_s *bmx, 
		 struct cplan9_s *hmm, float sc, int i0, int j0, ESL_DSQ *dsq)
{
  if(dsq == NULL)
    cm_Fail("in cp9_DebugCheckFB(), dsq is NULL.");

  int k, i;
  float max_diff;  /* maximum allowed difference between sc and 
		    * sum_k f[i][k] * b[i][k] for any i */
  float diff;
  int fb_sum;
  float fb_sc;
  int   W;		/* subsequence length */
  int   ip;		/* i': relative position in the subsequence  */
  int to_add;

  W  = j0-i0+1;		/* the length of the subsequence */

  max_diff = 0.01;
  /*printf("sc: %f\n", sc);*/

  /* In all possible paths through the model, each residue of the sequence must have 
   * been emitted by exactly 1 insert or match state. */
  for (ip = 1; ip <= W; ip++)
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      fb_sum = -INFTY;
      //fb_sum = ILogsum(fb_sum, (fmx->mmx[ip][0] + bmx->mmx[ip][0]));
      //fb_sum = ILogsum(fb_sum, (fmx->imx[ip][0] + bmx->imx[ip][0]));
      for (k = 0; k <= hmm->M; k++) 
	{
	  if(fmx->mmx[ip][k] == -INFTY) to_add = -INFTY;
	  else if(bmx->mmx[ip][k] == -INFTY) to_add = -INFTY;
	  else 
	    {
	      to_add = fmx->mmx[ip][k] + bmx->mmx[ip][k];
	      if(k > 0) to_add -= hmm->msc[dsq[i]][k];
	    }
	  /* hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx
	   * unless, we're talking about M_0, the B state, it doesn't emit */
	  fb_sum = ILogsum(fb_sum, to_add);
	  printf("fmx->mmx[ip:%d][k:%d]: %d\n", ip, k, fmx->mmx[ip][k]);
	  printf("bmx->mmx[ip:%d][k:%d]: %d sum: %d\n", ip, k, (bmx->mmx[ip][k]-hmm->msc[dsq[i]][k]), fb_sum);
	  if(fmx->imx[ip][k] == -INFTY) to_add = -INFTY;
	  else if(bmx->imx[ip][k] == -INFTY) to_add = -INFTY;
	  else 
	    {
	      to_add = fmx->imx[ip][k] + bmx->imx[ip][k]; 
	      to_add -= hmm->isc[dsq[i]][k];
	    }
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  fb_sum = ILogsum(fb_sum, to_add);
	  printf("fmx->imx[ip:%d][k:%d]: %d\n", ip, k, fmx->imx[ip][k]);
	  printf("bmx->imx[ip:%d][k:%d]: %d sum: %d\n", ip, k, (bmx->imx[ip][k]-hmm->isc[dsq[i]][k]), fb_sum);
	  /*fb_sum = ILogsum(fb_sum, fmx->dmx[ip][k] + bmx->dmx[ip][k]);*/
	  printf("fmx->dmx[ip:%d][k:%d]: %d\n",   ip, k, fmx->dmx[ip][k]);
	  printf("bmx->dmx[ip:%d][k:%d]: %d\n\n", ip, k, bmx->dmx[ip][k]);
	}
      fb_sc  = Scorify(fb_sum);
      diff = sc - fb_sc;
      printf("ip: %d sc: %f diff: %f\n", ip, fb_sc, diff);
      if(diff < 0.) diff *= -1.;
      if(diff > max_diff)
	{
	  printf("ERROR, fb_sc[%d]: %f too different from P(x|hmm): %f\n", i, fb_sc, sc);
	  exit(1);
	}
    }
}


/*********************************************************************
 * Function: cp9_RelaxRootBandsForSearch()
 * 
 * Purpose:  In cp9_HMM2ijBands(), ROOT_S (state 0) sets imin[0]=imax[0]=i0,
 *           and jmin[0]=jmax[0]=j0, which is important for alignment,
 *           but during search enforces that the optimal alignment start
 *           at i0 and end at j0, but when searching we want to relax this
 *           requirement in case a higher scoring parse has different endpoints.
 *           This function sets imax[0] as maximum i = imax[1] (ROOT_IL) and
 *           jmin[0] = jmin[2] (ROOT_IR).
 *           
 * Args:
 * cm               the cm
 * int *imin        imin[v] = first position in band on i for state v
 * int *imax        imax[v] = last position in band on i for state v
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 */
void
cp9_RelaxRootBandsForSearch(CM_t *cm, int *imin, int *imax, int *jmin, int *jmax)
{
  int y, yoffset;
  /* look at all children y of ROOT_S (v == 0) and set:
   * imin[0] = min_y imin[y];
   * imax[0] = max_y imax[y];
   * jmin[0] = min_y jmin[y];
   * jmax[0] = max_y jmax[y];
   */
  /* First look at children of 0 (these probs will be 0. if local begins on, but it doesn't matter for our purposes here) */
  for (yoffset = 0; yoffset < cm->cnum[0]; yoffset++) {
    y = cm->cnum[0] + yoffset;
    imin[0] = ESL_MIN(imin[0], imin[y]);
    imax[0] = ESL_MAX(imax[0], imax[y]);
    jmin[0] = ESL_MIN(jmin[0], jmin[y]);
    jmax[0] = ESL_MAX(jmax[0], jmax[y]);
  }
  /* now for possible local begins */
  if(cm->flags & CMH_LOCAL_BEGIN) {
    for (y = 1; y < cm->M; y++) {
      if(NOT_IMPOSSIBLE(cm->beginsc[y])) { 
	imin[0] = ESL_MIN(imin[0], imin[y]);
	imax[0] = ESL_MAX(imax[0], imax[y]);
	jmin[0] = ESL_MIN(jmin[0], jmin[y]);
	jmax[0] = ESL_MAX(jmax[0], jmax[y]);
      }
    }
  }
}


/*
 * Function: cp9_ValidateBands()
 * Incept:   EPN, Wed Nov 14 15:49:08 2007
 * Purpose:  Validate the info in CP9Bands_t data structure is internally
 *           consistent.
 *           
 * Args:     cm     the cm
 *           errbuf char buffer for error message
 *           cp9b   the CP9 bands object 
 *           i0     first residue we can possibly allow as valid j
 *           j0     final residue we can possibly allow as valid j
 *
 * Returns: eslOK, or, if error, other status code and filled errbuf
 */
int
cp9_ValidateBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int i0, int j0)
{
  int v;            /* counter over states of the CM */
  int jp;           /* counter over valid j's, but offset. jp+jmin[v] = actual j */
  int sd;           /* minimum d allowed for a state, ex: MP_st = 2, ML_st = 1. etc. */
  int hd_needed;
  int j;

  if(cm->M    != cp9b->cm_M)  ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cm->M != cp9b->cm_M\n");
  if(cm->clen != cp9b->hmm_M) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cm->clen != cp9b->hmm_M\n");
  
  hd_needed = 0; 
  for(v = 0; v < cp9b->cm_M; v++) 
    hd_needed += cp9b->jmax[v] - cp9b->jmin[v] + 1;
  if(hd_needed != cp9b->hd_needed) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cp9b->hd_needed inconsistent.");

  for(v = 0; v < cm->M; v++) {
    if(cm->sttype[v] == E_st) {
      for(jp = 0; jp <= (cp9b->jmax[v]-cp9b->jmin[v]); jp++) {
	if(cp9b->hdmin[v][jp] != 0) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cp9b->hdmin for E state is inconsistent.");
	if(cp9b->hdmax[v][jp] != 0) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cp9b->hdmin for E state is inconsistent.");
      }
    }
    else {
      sd = StateDelta(cm->sttype[v]);
      for(jp = 0; jp <= (cp9b->jmax[v]-cp9b->jmin[v]); jp++) {
	j = jp+cp9b->jmin[v];
	if(cp9b->hdmin[v][jp] != ESL_MAX((j - cp9b->imax[v] + 1), sd)) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cp9b->hdmin %d (sd: %d) for state %d, j: %d imax[v]: %d is inconsistent.", cp9b->hdmin[v][jp], sd, v, j, cp9b->imax[v]);
	if(cp9b->hdmax[v][jp] != ESL_MAX((j - cp9b->imin[v] + 1), sd)) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cp9b->hdmax %d (sd: %d) for state %d, j: %d imin[v]: %d is inconsistent.", cp9b->hdmax[v][jp], sd, v, j, cp9b->imin[v]);
      }
    }

    for(j = cp9b->jmin[v]; j <= cp9b->jmax[v]; j++) {
      if(j < i0) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), j: %d outside i0:%d..j0:%d is within v's j band: jmin[%d]: %d jmax[%d]: %d\n", j, i0, j0, v, cp9b->jmin[v], v, cp9b->jmax[v]);
      if(j > j0) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), j: %d outside i0:%d..j0:%d is within v's j band: jmin[%d]: %d jmax[%d]: %d\n", j, i0, j0, v, cp9b->jmin[v], v, cp9b->jmax[v]);
      if(cp9b->hdmin[v][(j-cp9b->jmin[v])] < StateDelta(cm->sttype[v])) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), v: %d j: %d hdmin[v][jp_v:%d] : %d less than StateDelta for v: %d\n", v, j, (j-cp9b->jmin[v]), cp9b->hdmin[v][(j-cp9b->jmin[v])], StateDelta(cm->sttype[v]));
     if(cp9b->hdmax[v][(j-cp9b->jmin[v])] < StateDelta(cm->sttype[v])) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), v: %d j: %d hdmax[v][jp_v:%d] : %d less than StateDelta for v: %d\n", v, j, (j-cp9b->jmin[v]), cp9b->hdmax[v][(j-cp9b->jmin[v])], StateDelta(cm->sttype[v]));
    }

    if(cp9b->imin[v] < cp9b->imin[0]) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cp9b->imin[v:%d]: %d less than cp9b->imin[0]: %d.", v, cp9b->imin[v], cp9b->imin[0]);
    if(cp9b->jmax[v] > cp9b->jmax[0]) ESL_FAIL(eslEINVAL, errbuf, "cp9_ValidateBands(), cp9b->jmax[v:%d]: %d greater than cp9b->jmax[0]: %d.", v, cp9b->jmax[v], cp9b->jmax[0]);
  }
  return eslOK;
}

/*
 * Function: cp9_GrowHDBands()
 * 
 * Incept:   EPN, Thu Oct 25 13:24:29 2007
 * Purpose:  Rearrange CP9 hdmin and hdmax pointers for a new sequence
 *           based on j bands (jmin and jmax). If the currently allocated
 *           size for hdmin, hdmax is not big enough, reallocate them.
 *
 * Args:
 * CP9Bands_t cp9b    the CP9 Bands object.
 * errbuf   char buffer for error messages
 *           
 * Returns: eslOK on success, eslEMEM if memory allocation error
 */
int 
cp9_GrowHDBands(CP9Bands_t *cp9b, char *errbuf)
{
  int status;
  int v;
  int cur_size = 0;
  int jbw;

  /* count size we need for hdmin/hdmax given current jmin, jmax */
  cp9b->hd_needed = 0; /* we'll rewrite this */
  for(v = 0; v < cp9b->cm_M; v++) 
    cp9b->hd_needed += cp9b->jmax[v] - cp9b->jmin[v] + 1;

  if(cp9b->hd_alloced < cp9b->hd_needed) {
    void *tmp;
    if(cp9b->hdmin_mem == NULL) ESL_ALLOC(cp9b->hdmin_mem, sizeof(int) * cp9b->hd_needed);
    else                        ESL_RALLOC(cp9b->hdmin_mem, tmp, sizeof(int) * cp9b->hd_needed);
    if(cp9b->hdmax_mem == NULL) ESL_ALLOC(cp9b->hdmax_mem, sizeof(int) * cp9b->hd_needed);
    else                        ESL_RALLOC(cp9b->hdmax_mem, tmp, sizeof(int) * cp9b->hd_needed);
  }
 
  /* set pointers */
  cur_size = 0;
  for(v = 0; v < cp9b->cm_M; v++) { 
    cp9b->hdmin[v] = cp9b->hdmin_mem + cur_size;
    cp9b->hdmax[v] = cp9b->hdmax_mem + cur_size;
    jbw = cp9b->jmax[v] - cp9b->jmin[v] + 1;
    ESL_DASSERT1((jbw >= 0));
    cur_size += jbw;
  }
  cp9b->hd_alloced = cur_size;
  ESL_DASSERT1((cp9b->hd_alloced == cp9b->hd_needed));
  return eslOK;
  
 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.");
}

#if 0
/* Here are the P7 versions of the functions, for reference */

/*****************************************************************************
 * EPN 04.03.06
 * Function: P7_hmm_band_bounds()
 *
 * Purpose:  Determine the band on all HMM states given the posterior
 *           matrices. Do this by summing log probabilities, starting
 *           at the sequence ends, and creeping in, until the half the
 *           maximum allowable probability excluded is reached on each
 *           side respectively.
 * 
 * below * = 'i', 'm' or 'd', for either (i)nsert, (m)atch or (d)elete states
 * arguments:
 *
 * int post         posterior matrix for *mx (matches, inserts or deletes)
 *                  2D int array [0.1..N][0.1..M] M = num nodes in HMM
 * int   L          length of sequence (num rows of post matrix)
 * int   M          number of nodes in HMM (num columns of post matrix)
 * int  *isum_pn    [1..M] sum_pn[k] = sum over i of log probabilities
 *                  from post->*mx[i][k]
 *                  if NULL: don't use sums, just use raw log probs
 * int pn_min       pn_min[k] = first position in band for * state of node k
 *                  to be filled in this function.
 * int pn_max       pn_max[k] = last position in band for * state of node k
 *                  to be filled in this function.
 * double p_thresh  the probability mass we're requiring is within each band
 * int state_type   HMMMATCH, HMMINSERT, or HMMDELETE, for deletes we have to deal
 *                  with the CM->HMM delete off-by-one issue (see code below).
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *****************************************************************************/
void
P7_hmm_band_bounds(int **post, int L, int M, int *isum_pn, int *pn_min, int *pn_max, 
		   double p_thresh, int state_type, int debug_level)
{
  int k;         /* counter over nodes of the model */
  int lmass_exc; /* the log of the probability mass currently excluded on the left*/
  int rmass_exc; /* the log of the probability mass currently excluded on the right*/
  int log_p_side;/* the log probability we're allowed to exclude on each side */
  int curr_log_p_side; /* the log probability we're allowed to exclude on each side for the current state */
  int argmax_pn; /* for curr state, the state with the highest log p, 
	          * IFF we determine the entire sequence is outside the
		  * band for a state, we set the band to a single position,
		  * the most likely one. Therefore its value is only 
		  * relevant (and valid!) if pmin[k] == L. 
		  * (otherwise we'd have some positions within the band).*/
  int max_post;  /* post[argmax_pn][k] for current k */
  /* NOTE: all *log_p* structures, and other structures that hold log probs
   * don't actually hold log probs. but scores, which are scaled up 1000X (INTSCALE)
   */

  log_p_side = Prob2Score(((1. - p_thresh)/2.), 1.); /* allowable prob mass excluded on each side */

  /* step through each node */
  for(k = 1; k <= M; k++)
    {
      curr_log_p_side = log_p_side; 
      if(isum_pn != NULL) 
	curr_log_p_side += isum_pn[k]; /* if we use sums strategy, normalize
					* so total prob of entering k = 1. */
      argmax_pn = 1;
      max_post = post[1][k];
      pn_min[k] = 2;
      pn_max[k] = L-1;
      lmass_exc = post[(pn_min[k]-1)][k];
      rmass_exc = post[(pn_max[k]+1)][k];
      /*creep in on the left, until we exceed our allowable prob mass to exclude.*/
      while(pn_min[k] <= L && lmass_exc <= (curr_log_p_side))
	{
	  if(post[pn_min[k]][k] > max_post) /* save info on most likely posn 
					     * in case whole seq is outside band */
	    {
	      max_post = post[pn_min[k]][k];
	      argmax_pn = pn_min[k];
	    }
	  lmass_exc = ILogsum(lmass_exc, post[pn_min[k]][k]);
	  pn_min[k]++;
	}
      /* we went one posn too far, back up*/
      pn_min[k]--;
      
      /*creep in on the right, until we exceed our allowable prob mass to exclude.*/
      while(pn_max[k] >= 1 && rmass_exc <= (curr_log_p_side))
	{
	  rmass_exc = ILogsum(rmass_exc, post[pn_max[k]][k]);
	  pn_max[k]--;
	}
      /* we went one posn too far, back up*/
      pn_max[k]++;
      
      if(pn_min[k] > pn_max[k])
	{
	  /* The sum of the posteriors for all posns for this state
	   * is less than tau. Current strategy, set band to a single
	   * cell, the most likely posn found when creeping in from left.
	   */
	  pn_min[k] = argmax_pn;
	  pn_max[k] = argmax_pn;
	}
      if(state_type == HMMDELETE)
	{
	  /* We have to deal with off by ones in the delete states 
	   * e.g. pn_min_d[k] = i, means posn i was last residue emitted
	   * prior to entering node k's delete state. However, for a CM,
	   * if a delete states sub-parsetree is bounded by i' and j', then
	   * positions i' and j' HAVE YET TO BE EMITTED.
	   */
	  pn_min[k]++;
	  pn_max[k]++;
	  /* In plan 7 HMMs, a delete state can only be entered after
	   * visiting at least one match state (M_1). But in a CM we 
	   * can start in deletes, so we explicitly check and fix. 
	   */
	  if(pn_min[k] == 2) pn_min[k] = 1; 
	}
    }
}


/**************************************************************************
 * P7_* functions no longer supported as of 10.26.06, 
 *      They remain here for reference.
 *      This code was written before the CMCP9Map_t data structure
 *      was introduced.
 * 
 * simple_cp9_HMM2ijBands() is an attempt I made to simplify the horribly
 *   convoluted cp9_HMM2ijBands() function, but it wasn't nearly as effective,
 *   and often obscured the optimal alignment, so it was abandoned.
 */
/**************************************************************************
 * EPN 03.26.06
 * P7_map_cm2hmm_and_hmm2cm()
 *
 * Purpose:  Determine maps between a CM and an HMM by filling 3 multi-dimensional
 *           arrays. All arrays must be pre-allocated and freed by caller.
 * Args:    
 * CM_t *cm          - the CM
 * cplan9_s *hmm     - the HMM
 * int *node_cc_left - consensus column each node's left emission maps to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *node_cc_right- consensus column each node's right emission corresponds to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *cc_node_map  - node that each consensus column maps to (is modelled by)
 *                     [1..hmm_ncc]
 * int **cs2hn_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
 *                     that contains state that maps to CM state, -1 if none.
 * int **cs2hs_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 2 elements for up to 2 matching HMM states, 
 *                     value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state, -1 if none.
 * 
 *                     For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
 *                                  maps to CM state v.
 * 
 * int ***hns2cs_map  - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                      HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                      2 matching CM states, value: CM states that map, -1 if none.
 *              
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * Returns: (void) 
 */
void
P7_map_cm2hmm_and_hmm2cm(CM_t *cm, struct plan7_s *hmm, int *node_cc_left, int *node_cc_right, int *cc_node_map, int ***ret_cs2hn_map, int ***ret_cs2hs_map, int ****ret_hns2cs_map, int debug_level)
{

  int status;
  int k; /* HMM node counter */
  int ks; /* HMM state counter (0(Match) 1(insert) or 2(delete)*/
  int n; /* CM node that maps to HMM node k */
  int nn; /* CM node index */
  int n_begr; /* CM node index */
  int is_left; /* TRUE if HMM node k maps to left half of CM node n */
  int is_right; /* TRUE if HMM node k maps to right half of CM node n */
  int v; /* state index in CM */
  int v1, v2;
  int **cs2hn_map;
  int **cs2hs_map;
  int ***hns2cs_map;

  /* Allocate the maps */
  ESL_ALLOC(cs2hn_map, sizeof(int *) * (cm->M+1));
  ESL_ALLOC(cs2hn_map[0], sizeof(int) * 2 * (cm->M+1));
  for(v = 0; v <= cm->M; v++) cs2hn_map[v]     = cs2hn_map[0] + v * 2; 
  
  ESL_ALLOC(cs2hs_map, sizeof(int *) * (cm->M+1));
  ESL_ALLOC(cs2hs_map[0], sizeof(int) * 2 * (cm->M+1));
  for(v = 0; v <= cm->M; v++) cs2hs_map[v]     = cs2hs_map[0] + v * 2; 

  ESL_ALLOC(hns2cs_map, sizeof(int *) * (hmm->M+1));
  ESL_ALLOC(hns2cs_map[0], sizeof(int) * 3 * 2 * (cm->M+1));
  for(k = 0; k <= hmm->M; k++) 
    {
      hns2cs_map[k] = hns2cs_map[0] + k * 3 * 2; 
      for(ks = 0; ks < 3; ks++) 
	hns2cs_map[k][ks]    = hns2cs_map[k] + ks; 
    }	
      
  /* Initialize the maps */
  for(v = 0; v <= cm->M; v++)
    {
      cs2hn_map[v][0] = -1;
      cs2hn_map[v][1] = -1;
      cs2hs_map[v][0] = -1;
      cs2hs_map[v][1] = -1;
    }
  for(k = 0; k <= hmm->M; k++)
    {
      for(ks = 0; ks < 3; ks++)
	{
	  hns2cs_map[k][ks][0] = -1;
	  hns2cs_map[k][ks][1] = -1;
	}
    }

  /* One of the few differences b/t this version (p7) of the function
   * and the CP9 version, P7 HMMs have no node 0. We say nothing
   * maps to the ROOT node states. (even though B maps to ROOT_S,
   * N 'sort of' maps to ROOT_IL. Another difference is that Plan7 
   * HMMs don't have  a D_1, I_M and D_M state, so they are not 
   * mapped inside the following for loop.
   */
  
  /* Step through HMM nodes, filling in maps as we go */
  for(k = 1; k <= hmm->M; k++)
    {
      n = cc_node_map[k];
      if(node_cc_left[n] == k)
	{
	  is_left = TRUE;
	  is_right = FALSE;
	}
      else if(node_cc_right[n] == k)
	{
	  is_left = FALSE;
	  is_right = TRUE;
	}
      switch(cm->ndtype[n])
	{
	case ROOT_nd:
	case BIF_nd:
	case BEGL_nd:
	case BEGR_nd:
	case END_nd:
	  printf("ERROR: HMM node k doesn't map to MATP, MATR or MATL\n");
	  exit(1);
	  break;
	  
	case MATP_nd:
	  if(is_left)
	    {
	      ks = 0; /*match*/
	      v = cm->nodemap[n]; /*MATP_MP*/
	      map_helper(cp9map, k, ks, v);
	      v = cm->nodemap[n] + 1; /*MATP_ML*/
	      map_helper(cp9map, k, ks, v);

	      ks = 1; /*insert*/
	      if(k != hmm->M)
		{
		  v = cm->nodemap[n] + 4; /*MATP_IL*/
		  map_helper(cp9map, k, ks, v);
		}
	      ks = 2; /*delete*/
	      if(k != hmm->M && k != 1)
		{
		  v = cm->nodemap[n] + 2; /*MATP_MR*/
		  map_helper(cp9map, k, ks, v);
		  v = cm->nodemap[n] + 3; /*MATP_D*/
		  map_helper(cp9map, k, ks, v);
		}
	    }
	  else if(is_right)
	    {
	      ks = 0; /*match*/
	      v = cm->nodemap[n]; /*MATP_MP*/
	      map_helper(cp9map, k, ks, v);
	      v = cm->nodemap[n] + 2; /*MATP_MR*/
	      map_helper(cp9map, k, ks, v);

	      ks = 1; /*insert*/
	      /* whoa... careful, we want the CM state that will insert to the RIGHT
	       * of column k (the right consensus column modelled by MATP node n),
	       * but MATP_IR inserts to the LEFT of column k.
	       * What we need to determine is the CM node nn that models column k+1,
	       * and further which half (left or right) of nn models k+1, then
	       * we can map the HMM state to the correct CM state (see code).
	       */
	      if(k != hmm->M) /* There is o P7 I_M state */
		{
		  nn = cc_node_map[k+1];
		  if(node_cc_left[nn] == (k+1))
		    {
		      /* find the closest BEGR node above node nn */
		      n_begr = nn;
		      while(n_begr >= 0 && (cm->ndtype[n_begr] != BEGR_nd))
			n_begr--;
		      if(n_begr == -1)
			{
			  printf("ERROR: can't find BEGR node above node %d\n", nn);
			  printf("k is %d\n", k);
			  exit(1);
			}
		      v = cm->nodemap[n_begr] + 1; /*BEGR_IL*/
		      map_helper(cp9map, k, ks, v);
		    }
		  else if(node_cc_right[nn] == (k+1))
		    {
		      /*simple*/
		      if(cm->ndtype[nn] == MATP_nd)
			{
			  v = cm->nodemap[nn] + 5; /*MATP_IR*/
			  map_helper(cp9map, k, ks, v);
			}
		      else if(cm->ndtype[nn] == MATR_nd)
			{
			  v = cm->nodemap[nn] + 2; /*MATR_IR*/
			  map_helper(cp9map, k, ks, v);
			}
		    }
		} /* end of if (k != hmm->M) */
	      /* NOT DONE YET, the MATP_IR has to map to an HMM state,
	       * if the previous column (k-1) is modelled by a CM MATR or 
	       * MATP node, then the above block will take care of this situation
	       * (in the previous iteration of this loop when k = k-1), 
	       * HOWEVER, if (k-1) is modelled by a MATL, then this 
	       * MATP_IR's contribution to the HMM will be ignored, 
	       * unless we do something about it. 
	       */ 
	      if(node_cc_left[cc_node_map[k-1]] == (k-1)) /*k-1 modelled by MATL*/
		{
		  if(cm->ndtype[cc_node_map[k-1]] != MATL_nd)
		    {
		      if(cm->ndtype[cc_node_map[k-1]] == MATP_nd)
			{
			  /* A rare, but possible case. Previous column
			   * k-1 column is modelled by left half of the MATP
			   * node whose right half models column k.
			   * Proceed below. 
			   */
			}
		      else
			{
			  printf("ERROR, full understanding of the CM architecture remains elusive 0)...\n");
			  exit(1);
			}
		    }
		  v = cm->nodemap[n] + 5; /*MATP_IR*/
		  if(k != hmm->M) /* There is no plan 7 I_M state */
		    map_helper(cp9map, (k-1), ks, v);
		}
	      
	      ks = 2; /*delete*/
	      v = cm->nodemap[n] + 1; /*MATP_ML*/
	      if(k != 1 && k != hmm->M) /* There is no plan 7 D_1 or D_M state */
		map_helper(cp9map, k, ks, v);
	      
	      v = cm->nodemap[n] + 3; /*MATP_D*/
	      if(k != 1 && k != hmm->M) /* There is no plan 7 D_1 or D_M state */
		map_helper(cp9map, k, ks, v);
	    }
	  break;

	case MATL_nd:
	  ks = 0; /*match*/
	  v = cm->nodemap[n]; /*MATL_ML*/
	  map_helper(cp9map, k, ks, v);

	  ks = 1; /*insert*/
	  v = cm->nodemap[n] + 2; /*MATL_IL*/
	  if(k != hmm->M) /* There is no plan 7 I_M state */
	    map_helper(cp9map, k, ks, v);

	  ks = 2; /*delete*/
	  v = cm->nodemap[n] + 1; /*MATL_D*/
	  if(k != 1 && k != hmm->M) /* There is no plan 7 D_1 or D_M state */
	    map_helper(cp9map, k, ks, v);

	  break;

	case MATR_nd:
	  ks = 0; /*match*/
	  v = cm->nodemap[n]; /*MATR_MR*/
	  map_helper(cp9map, k, ks, v);

	  ks = 1; /*insert*/
	  /* whoa... careful, we want the CM state that will insert to the RIGHT
	   * of column k (the consensus column modelled by MATR node n),
	   * but MATR_IR inserts to the LEFT of column k.
	   * What we need to determine is the CM node nn that models column k+1,
	   * and further which half (left or right) of nn models k+1, then
	   * we can map the HMM state to the correct CM state (see code).
	   */
	  if(k != hmm->M) /* There is no plan 7 I_M state */
	    {
	      nn = cc_node_map[k+1];
	      if(node_cc_left[nn] == (k+1))
		{
		  /* find the closest BEGR node above node nn */
		  n_begr = nn;
		  while((cm->ndtype[n_begr] != BEGR_nd) && n_begr >= 0)
		    n_begr--;
		  if(n_begr == -1)
		    {
		      printf("ERROR: can't find BEGR node above node %d\n", nn);
		      exit(1);
		    }
		  v = cm->nodemap[n_begr] + 1; /*BEGR_IL*/
		  map_helper(cp9map, k, ks, v);
		}
	      else if(node_cc_right[nn] == (k+1))
		{
		  /*simple*/
		  if(cm->ndtype[nn] == MATP_nd)
		    {
		      v = cm->nodemap[nn] + 5;
		      map_helper(cp9map, k, ks, v);
		    }
		  else if(cm->ndtype[nn] == MATR_nd)
		    {
		      v = cm->nodemap[nn] + 2; /*MATP_IR*/
		      map_helper(cp9map, k, ks, v);
		    }
		}
	    } /* end of if (k != hmm->M) */
	  if(node_cc_left[cc_node_map[k-1]] == (k-1)) /*k-1 modelled by MATL*/
	    {
	      //	      printf("ERROR, full understanding of the CM architecture remains elusive (1)...\n");
	      exit(1);
	    }
	  
	  ks = 2; /*delete*/
	  v = cm->nodemap[n] + 1; /*MATR_D*/
	  if(k != 1 && k != hmm->M) /* There is no plan 7 D_1 or D_M state */
	    map_helper(cp9map, k, ks, v);
	  break;
	}
    }

  /* Check to make sure that insert states map to only 1 HMM node state. */
  for(v = 0; v <= cm->M; v++)
    {
      if((cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) && cs2hn_map[v][1] != -1)
	{
	  printf("ERROR during cs2hn_map construction\ncs2hn_map[%d][0]: %d | cs2hn_map[%d][1]: %d AND v is an INSERT state\n", v, cs2hn_map[v][0], v, cs2hn_map[v][1]);
	  exit(1);
	}
      /* each CM state should map to only 1 HMM state. */
    }


  /* print hns2cs_map, checking consistency with cs2hn_map and cs2hs_map along
     the way.  */
  for(k = 1; k <= hmm->M; k++)
    {
      for(ks = 0; ks < 3; ks++)
	{
	  v1 = hns2cs_map[k][ks][0];
	  v2 = hns2cs_map[k][ks][1];
	  if(debug_level > 1)
	    printf("hns2cs[%3d][%3d][0]: %3d | hns2cs[%3d][%3d[1]: %3d\n", k, ks, v1, k, ks, v2);
	  if(v1 != -1 && (cs2hn_map[v1][0] == k && cs2hs_map[v1][0] == ks))
	    {
	      /* okay */
	    }
	  else if(v1 != -1 && (cs2hn_map[v1][1] == k && cs2hs_map[v1][1] == ks))
	    {
	      /* okay */
	    }
	  else if(v2 != -1 && (cs2hn_map[v2][0] == k && cs2hs_map[v2][0] == ks))
	    {
	      /* okay */
	    }
	  else if(v2 != -1 && (cs2hn_map[v2][1] == k && cs2hs_map[v2][1] == ks))
	    {
	      /* okay */
	    }
	  else if(v1 == -1 && v2 == -1 && (k == 1 && ks == 2)) 
	    {
	      /*okay - D_0 maps to nothing */
	    }
	  else if(v1 == -1 && v2 == -1 && (k == hmm->M && ks == 1)) 
	    {
	      /*okay - D_0 maps to nothing */
	    }
	  else if(v1 == -1 && v2 == -1 && (k == hmm->M && ks == 2)) 
	    {
	      /*okay - D_0 maps to nothing */
	    }
	  else if(v1 == -1 && v2 == -1)
	    /* only D_1, D_M and I_M should map to nothing*/
	    {
	      /* not okay */
	      printf("maps inconsistent case 1, HMM node state (non D_0) maps to no CM state, v1: %d | v2: %d k: %d | ks: %d\n", v1, v2, k, ks);
	      exit(1);
	    }	      
	  else
	    {
	      /* not okay */
	      printf("maps inconsistent case 2 v1: %d | v2: %d k: %d | ks: %d\n", v1, v2, k, ks);
	      exit(1);
	    }
	}
    }
  *ret_cs2hn_map = cs2hn_map;
  *ret_cs2hs_map = cs2hs_map;
  *ret_hns2cs_map = hns2cs_map;
  return;
 ERROR:
  cm_Fail("Memory allocation error.\n");
}

/*****************************************************************************
 * EPN 11.03.05
 * Function: P7_last_hmm_insert_state_hack
 *
 * Purpose:  HMMER plan 7 doesn't have an insert node in
 *           its final node, so we have to use a hack
 *           to fill pn_min_i[hmm->M] and pn_max_i[hmm->M].
 *           We simply use the match node's band.
 * 
 * arguments:
 * int   M          number of nodes in HMM (num columns of post matrix)
 * int *pn_min_m    pn_min[k] = first position in band for match state of node k
 * int *pn_max_m    pn_max[k] = last position in band for match state of node k
 * int *pn_min_i    pn_min[k] = first position in band for insert state of node k
 *                  array element M is the only one filled in this function .
 * int *pn_max_i    pn_max[k] = last position in band for insert state of node k
 *                  array element M is the only one filled in this function .
 *****************************************************************************/
void
P7_last_hmm_insert_state_hack(int M,  int *pn_min_m, int *pn_max_m, int *pn_min_i, int *pn_max_i)
{
  pn_min_i[M] = pn_min_m[M];
  pn_max_i[M] = pn_max_m[M];
}
/*****************************************************************************
 * EPN 12.18.05
 * Function: P7_last_and_first_delete_state_hack
 *
 * Purpose:  P7 HMMs don't have delete states in their 
 *           first or final node, so we have to use a hack
 *           to fill pn_min_d[1], pn_max_d[1], pn_min_d[hmm->M] and 
 *           pn_max_d[hmm->M]. We use the match state bands.
 * 
 * arguments:
 * int   M          number of nodes in HMM (num columns of post matrix)
 * int *pn_min_m    pn_min[k] = first position in band for match state of node k
 * int *pn_max_m    pn_max[k] = last position in band for match state of node k
 * int *pn_min_d    pn_min[k] = first position in band for delete state of node k
 *                  array element M is the only one filled in this function .
 * int *pn_max_d    pn_max[k] = last position in band for deletea state of node k
 *                  array element M is the only one filled in this function .
 * int L            length of target database sequence.
 *****************************************************************************/
void
P7_last_and_first_hmm_delete_state_hack(int M,  int *pn_min_m, int *pn_max_m, int *pn_min_d, int *pn_max_d, int L)
{
  /* Maybe I should be basing the delete bands for the first state on the
   * N state posteriors in the HMM, and for the last state on the 
   * C state posteriors in the HMM.
   */
  pn_min_d[1] = pn_min_m[1];
  pn_max_d[1] = pn_max_m[1]; 
  pn_min_d[M] = pn_max_m[M];
  pn_max_d[M] = pn_max_m[M];
}

/* Function: P7FullPosterior()
 * based on Ian Holmes' hmmer/src/postprob.c::P7EmitterPosterior()
 *
 * Purpose:  Combines Forward and Backward matrices into a posterior
 *           probability matrix.
 *           For emitters (match and inserts) the entries in row i of this 
 *           matrix are the logs of the posterior probabilities of each state 
 *           emitting symbol i of the sequence. For non-emitters 
 *           (main node deletes only (not XME, XMB)
 *           the entries in row i of this matrix are the logs of the posterior
 *           probabilities of each state being 'visited' when the last emitted
 *           residue in the parse was symbol i of the sequence (I think this
 *           is valid, but not sure (EPN)). The last point distinguishes this
 *           function from P7EmitterPosterior() which set all posterior values
 *           for for non-emitting states to -INFTY.
 *           The caller must allocate space for the matrix, although the
 *           backward matrix can be used instead (overwriting it will not
 *           compromise the algorithm).
 *           
 * Args:     L        - length of sequence
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           mx       - pre-allocated dynamic programming matrix
 *           
 * Return:   void
 */
void
P7FullPosterior(int L,
		 struct plan7_s *hmm,
		 struct dpmatrix_s *forward,
		 struct dpmatrix_s *backward,
		 struct dpmatrix_s *mx)
{
  int i;
  int k;
  int sc;

  sc = backward->xmx[0][XMN];
  for (i = L; i >= 1; i--)
    {
      mx->xmx[i][XMC] = forward->xmx[i-1][XMC] + hmm->xsc[XTC][LOOP] + backward->xmx[i][XMC] - sc;
      
      mx->xmx[i][XMJ] = forward->xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP] + backward->xmx[i][XMJ] - sc;
 
      mx->xmx[i][XMN] = forward->xmx[i-1][XMN] + hmm->xsc[XTN][LOOP] + backward->xmx[i][XMN] - sc;

      mx->xmx[i][XMB] = mx->xmx[i][XME] = -INFTY;
      
      for (k = 1; k < hmm->M; k++) {
	mx->mmx[i][k]  = backward->mmx[i][k];
	/* we don't want to account for possibility we came from a delete, first of all
	 * they don't emit so that means we'd be interested in the i index of the previous
	 * delete state (forward->dmx[i][k-1], however, that value has already accounted
	 * for the emission of position i; so its not of interest here.
	 */
	mx->mmx[i][k] += ILogsum(ILogsum(forward->mmx[i-1][k-1] + hmm->tsc[TMM][k-1],
					 forward->imx[i-1][k-1] + hmm->tsc[TIM][k-1]),
				 ILogsum(forward->xmx[i-1][XMB] + hmm->bsc[k],
					 forward->dmx[i-1][k-1] + hmm->tsc[TDM][k-1]));
	mx->mmx[i][k] -= sc;
	
	mx->imx[i][k]  = backward->imx[i][k];
	mx->imx[i][k] += ILogsum(forward->mmx[i-1][k] + hmm->tsc[TMI][k],
				 forward->imx[i-1][k] + hmm->tsc[TII][k]);
	mx->imx[i][k] -= sc;
	
	/*mx->dmx[i][k] = -INFTY;*/
	/* I think this is right?? */
	mx->dmx[i][k]  = backward->dmx[i][k];
	mx->dmx[i][k] += forward->dmx[i][k];
	mx->dmx[i][k] -= sc;
      }
      mx->mmx[i][hmm->M]  = backward->mmx[i][hmm->M];
      mx->mmx[i][hmm->M] += ILogsum(ILogsum(forward->mmx[i-1][hmm->M-1] + hmm->tsc[TMM][hmm->M-1],
					    forward->imx[i-1][hmm->M-1] + hmm->tsc[TIM][hmm->M-1]),
				    ILogsum(forward->xmx[i-1][XMB] + hmm->bsc[hmm->M],
					    forward->dmx[i-1][hmm->M-1] + hmm->tsc[TDM][hmm->M-1]));
      mx->mmx[i][hmm->M] -= sc;

      mx->imx[i][hmm->M] = mx->dmx[i][hmm->M] = mx->dmx[i][0] = -INFTY;
      
    }
  /*  for(i = 1; i <= L; i++)
    {
      for(k = 1; k < hmm->M; k++)
	{
	  temp_sc = Score2Prob(mx->mmx[i][k], 1.);
	  if(temp_sc > .0001)
	    printf("p7 mx->mmx[%3d][%3d]: %9d | %8f\n", i, k, mx->mmx[i][k], temp_sc);
	  temp_sc = Score2Prob(mx->imx[i][k], 1.);
	  if(temp_sc > .0001)
	    printf("p7 mx->imx[%3d][%3d]: %9d | %8f\n", i, k, mx->imx[i][k], temp_sc);
	  temp_sc = Score2Prob(mx->dmx[i][k], 1.);
	  if(temp_sc > .0001)
	    printf("p7 mx->dmx[%3d][%3d]: %9d | %8f\n", i, k, mx->dmx[i][k], temp_sc);
	}
    }
  */
}




/*****************************************************************************
 * EPN 12.16.05
 * Function: P7_ifill_post_sums()
 * based on: 
 * Function: ifill_post_sums()
 * EPN 11.23.05
 *
 * Purpose:  Given a posterior matrix post, where post->mmx[i][k]
 *           is the log odds score of the probability that
 *           match state k emitted position i of the sequence,
 *           sum the log probabilities that each state emitted
 *           each position. Do this for inserts and matches, and
 *           and deletes.
 * 
 * arguments:
 * dpmatrix_s *post dpmatrix_s posterior matrix, xmx, mmx, imx, dmx 
 *                  2D int arrays. [0.1..N][0.1..M]
 * int   L          length of sequence (num rows of matrix)
 * int   M          number of nodes in HMM (num columns of matrix)
 * int  *isum_pn_m  [1..M] sum_pn_m[k] = sum over i of log probabilities
 *                  from post->mmx[i][k]
 *                  filled in this function, must be freed by caller.
 * int  *isum_pn_i  [1..M] sum_pn_m[k] = sum over i of log probabilities
 *                  from post->imx[i][k]
 *                  filled in this function, must be freed by caller.
 * int  *isum_pn_d  [1..M] sum_pn_m[k] = sum over i of log probabilities
 *                  from post->dmx[i][k]
 *                  filled in this function, must be freed by caller.
 *****************************************************************************/
void
P7_ifill_post_sums(struct dpmatrix_s *post, int L, int M,
		  int *isum_pn_m, int *isum_pn_i, int *isum_pn_d)
{
  int i;            /* counter over positions of the sequence */
  int k;            /* counter over nodes of the model */
  
  /* step through each node, fill the post sum structures */
  for(k = 1; k <= M; k++)
    {
      isum_pn_m[k] = post->mmx[1][k];
      isum_pn_i[k] = post->imx[1][k];
      isum_pn_d[k] = post->dmx[1][k];
      for(i = 2; i <= L; i++)
	{
	  isum_pn_m[k] = ILogsum(isum_pn_m[k], post->mmx[i][k]);
	  isum_pn_i[k] = ILogsum(isum_pn_i[k], post->imx[i][k]);
	  isum_pn_d[k] = ILogsum(isum_pn_d[k], post->dmx[i][k]);
	} 
    }
}

/*************************************************************
 * EPN 10.10.05
 * Function: P7_debug_print_post_decode()
 *
 * Purpose:  Print the post decode matrix.
 * 
 * arguments:
 * L          length of sequence (num rows of matrix)
 * M          number of nodes in HMM (num cols of matrix)
 */

void
P7_debug_print_post_decode(int L, int M, struct dpmatrix_s *posterior)
{
  int i, k;
  printf("\nPrinting post decode matrix :\n");
  printf("************************************\n");
  for(i = 1; i <= L; i++)
    {
      printf("====================================\n");
      printf("score_pd:xmx[%3d][%3d(XMB)] : %.15f\n", i, XMB, DScore2Prob(posterior->xmx[i][XMB], 1.));
      printf("score_pd:xmx[%3d][%3d(XME)] : %.15f\n", i, XME, DScore2Prob(posterior->xmx[i][XME], 1.));
      printf("score_pd:xmx[%3d][%3d(XMC)] : %.15f\n", i, XMC, DScore2Prob(posterior->xmx[i][XMC], 1.));
      printf("score_pd:xmx[%3d][%3d(XMJ)] : %.15f\n", i, XMJ, DScore2Prob(posterior->xmx[i][XMJ], 1.));
      printf("score_pd:xmx[%3d][%3d(XMN)] : %.15f\n", i, XMN, DScore2Prob(posterior->xmx[i][XMN], 1.));
      for(k = 1; k <= M; k++)
	{
	  printf("------------------------------------\n");
	  printf("score_pd:mmx[%3d][%3d] : %.15f\n", i, k, DScore2Prob(posterior->mmx[i][k], 1.));
	  printf("score_pd:imx[%3d][%3d] : %.15f\n", i, k, DScore2Prob(posterior->imx[i][k], 1.));
	  printf("score_pd:dmx[%3d][%3d] : %.15f\n", i, k, DScore2Prob(posterior->dmx[i][k], 1.));
	  printf("pd:mmx[%3d][%3d] : %d\n", i, k, posterior->mmx[i][k]);
	  printf("pd:imx[%3d][%3d] : %d\n", i, k, posterior->imx[i][k]);
	  printf("pd:dmx[%3d][%3d] : %d\n", i, k, posterior->dmx[i][k]);
	}
    }
    printf("****************\n\n");
}

/* EPN 10.10.05
 * Function: P7_debug_print_dp_matrix()
 *
 * Purpose:  Print a dp matrix.
 * 
 * arguments:
 * L          length of sequence (num rows of matrix)
 * M          number of nodes in HMM (num cols of matrix)
 */

void
P7_debug_print_dp_matrix(int L, int M, struct dpmatrix_s *mx)
{
  int i, k;
  printf("\nPrinting dp matrix :\n");
  printf("************************************\n");
  for(i = 1; i <= L; i++)
    {
      printf("====================================\n");
      printf("score_dp:xmx[%3d][%3d (XMB)] : %.15f\n", i, XMB, DScore2Prob(mx->xmx[i][XMB], 1.));
      printf("score_dp:xmx[%3d][%3d (XME)] : %.15f\n", i, XME, DScore2Prob(mx->xmx[i][XME], 1.));
      printf("score_dp:xmx[%3d][%3d (XMC)] : %.15f\n", i, XMC, DScore2Prob(mx->xmx[i][XMC], 1.));
      printf("score_dp:xmx[%3d][%3d (XMJ)] : %.15f\n", i, XMJ, DScore2Prob(mx->xmx[i][XMJ], 1.));
      printf("score_dp:xmx[%3d][%3d (XMN)] : %.15f\n", i, XMN, DScore2Prob(mx->xmx[i][XMN], 1.));
      for(k = 1; k <= M; k++)
	{
	  printf("------------------------------------\n");
	  printf("score_dp:mmx[%3d][%3d] : %.15f\n", i, k, DScore2Prob(mx->mmx[i][k], 1.));
	  printf("score_dp:imx[%3d][%3d] : %.15f\n", i, k, DScore2Prob(mx->imx[i][k], 1.));
	  printf("score_dp:dmx[%3d][%3d] : %.15f\n", i, k, DScore2Prob(mx->dmx[i][k], 1.));
	  printf("dp:mmx[%3d][%3d] : %d\n", i, k, mx->mmx[i][k]);
	  printf("dp:imx[%3d][%3d] : %d\n", i, k, mx->imx[i][k]);
	  printf("dp:dmx[%3d][%3d] : %d\n", i, k, mx->dmx[i][k]);
	}
    }
    printf("****************\n\n");
}

/*****************************************************************************
 * ABANDONED 04.20.06 - in practice doesn't work as well as cp9_HMM2ijBands()
 *                      where precautions taken to ensure "safe" paths 
 *                      make it more robust to finding the optimal alignment.
 *                      This function, which doesn't take any precautions,
 *                      often messes up alignment. Though messy, cp9_HMM2ijBands()
 *                      is not worth overhauling right now. 
 * 
 * Function: simple_cp9_HMM2ijBands()
 * 
 * Purpose:  Determine the band for each cm state v on i (the band on the 
 *           starting index in the subsequence emitted from the subtree rooted
 *           at state v), and on j (the band on the ending index in the
 *           subsequence emitted from the subtree rooted at state v). 
 * 
 *           Some i and d bands are calculated from HMM bands on match and insert 
 *           and delete states from each node of the HMM that maps to a left emitting
 *           node of the CM (including MATP nodes). The HMM bands were
 *           calculated previously from the posterior matrices for mmx,
 *           imx and dmx from HMMER.
 * 
 *           Some j bands are calculated from HMM bands on match and insert and
 *           delete states from each node of the HMM that maps to a right emitting
 *           node of the CM (including MATP nodes). 
 * 
 *           i and j bands that cannot be directly determined from the
 *           HMM bands are inferred based on the constraints imposed
 *           on them by the i and j bands that CAN be determined from
 *           the HMM bands.
 *             
 * arguments:
 *
 * CM_t *cm         the CM 
 * int  ncc         number of consensus columns in HMM and CM seed alignment
 * int *node_cc_left consensus column each node's left emission maps to
 *                   [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *node_cc_right consensus column each node's right emission corresponds to
 *                    [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *cc_node_map  node that each consensus column maps to (is modelled by)
 *                     [1..ncc]
 * int  L           length of sequence we're aligning
 * int *pn_min_m    pn_min_m[k] = first position in HMM band for match state of HMM node k
 * int *pn_max_m    pn_max_m[k] = last position in HMM band for match state of HMM node k
 * int *pn_min_i    pn_min_i[k] = first position in HMM band for insert state of HMM node k
 * int *pn_max_i    pn_max_i[k] = last position in HMM band for insert state of HMM node k
 * int *pn_min_d    pn_min_d[k] = first position in HMM band for delete state of HMM node k
 * int *pn_max_d    pn_max_d[k] = last position in HMM band for delete state of HMM node k
 * int *imin        imin[v] = first position in band on i for state v
 *                  to be filled in this function. [1..M]
 * int *imax        imax[v] = last position in band on i for state v
 *                  to be filled in this function. [1..M]
 * int *jmin        jmin[v] = first position in band on j for state v
 *                  to be filled in this function. [1..M]
 * int *jmax        jmax[v] = last position in band on j for state v
 *                  to be filled in this function. [1..M]
 * int **cs2hn_map  2D CM state to HMM node map, 1st D - CM state index
 *                  2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
 *                  that contains state that maps to CM state, -1 if none.
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *****************************************************************************/
void
simple_cp9_HMM2ijBands(CM_t *cm, int ncc, int *node_cc_left, int *node_cc_right, 
		    int *pn_min_m, int *pn_max_m, int *pn_min_i, int *pn_max_i, 
		    int *pn_min_d, int *pn_max_d, int *imin, int *imax, 
		    int *jmin, int *jmax, int **cs2hn_map, int **cs2hs_map, 
		    int debug_level)
{

  int v;
  int n;
  int prev_imin, prev_imax;
  int prev_jmin, prev_jmax;
  int is_left, is_right;

  int k_left, k_right;
  int ks_left, ks_right;
  int careful;

  int left_unknown;
  int right_unknown; 
  int last_E;
  int y;
  int at_new_node;
  int prev_n;
  int set_left_known;
  int set_right_known;
  careful = TRUE;
  set_left_known = FALSE;
  set_right_known = FALSE;

  prev_n = cm->ndidx[cm->M-1];
  for (v = (cm->M-1); v >= 0; v--)
    {
      is_left  = FALSE;
      is_right = FALSE;
      n        = cm->ndidx[v];    /* CM node that contains state v */

      if(prev_n != n)
	{
	  at_new_node = TRUE;
	  if(set_left_known == TRUE)
	    {
	      left_unknown = FALSE;
	      for(y = v+1; y <= last_E; y++)
		{
		  imin[y] = prev_imin;
		  imax[y] = prev_imax;
		}
	      set_left_known = FALSE;
	    }
	  if(set_right_known == TRUE)
	    {
	      right_unknown = FALSE;
	      for(y = v+1; y <= last_E; y++)
		{
		  jmin[y] = prev_jmin;
		  jmax[y] = prev_jmax;
		}
	      set_right_known = FALSE;
	    }
	}
      else
	{
	  at_new_node = FALSE;
	}
      prev_n = n;

      if(cm->sttype[v] == E_st)
	{
	  last_E = v;
	  left_unknown  = TRUE;
	  right_unknown = TRUE;
	  prev_imin     = 0;
	  prev_imax     = 0;
	  prev_jmin     = 0;
	  prev_jmax     = 0;
	}
      else if(cm->sttype[v] == B_st)
	{
	  /* special case, use i band from left child and j band from right child */
	  imin[v] = imin[cm->cfirst[v]];
	  imax[v] = imax[cm->cfirst[v]];
	  jmin[v] = jmin[cm->cnum[v]];
	  jmax[v] = jmax[cm->cnum[v]];
	}
      else
	{
	  if((cm->sttype[v] == IR_st) ||
	     (cm->sttype[v] == MR_st && cm->ndtype[n] != MATP_nd) ||
	     (cm->sttype[v] == D_st && cm->ndtype[n] == MATR_nd))
	    {
	      /* only for these cases will the CM state v map to only 1 HMM state that
	       * is only giving information for bands on j, and not on i. 
	       */
	      k_right  = cs2hn_map[v][0]; /* HMM state 1 (0(match), 1(insert), or 2(delete) that maps to v */
	      ks_right = cs2hs_map[v][0]; /* HMM state 2 (0(match), 1(insert), or 2(delete) that maps to v 
					   * (-1 if none) */
	      k_left = -1;
	      ks_left = -1;
	    }
	  else
	    {
	      /* for these cases, the CM state v either maps to 1 or 2 HMM states, and
	       * either gives us info on bands on i, j or both. 
	       */
	      k_left   = cs2hn_map[v][0]; /* HMM node 1 that maps to v */
	      ks_left  = cs2hs_map[v][0]; /* HMM state 1 (0(match), 1(insert), or 2(delete) that maps to v */
	      k_right  = cs2hn_map[v][1]; /* HMM node 2 that maps to v (-1 if none)*/
	      ks_right = cs2hs_map[v][1]; /* HMM state 2 (0(match), 1(insert), or 2(delete) that maps to v 
					   * (-1 if none) */
	      /* the way cs2hs_map is constructed, by moving left to right in the HMM, guarantees
	       * that FOR THESE CASES: cs2hn_map[v][0] < cs2hn_map[v][1]; 
	       * so cs2hn_map[v][0] = the HMM node mapping to the left half of CM state v (or -1 if none)
	       *  & cs2hn_map[v][1] = the HMM node mapping to the right half of CM state v (or -1 if none)
	       */
	      //printf("normal case v: %d | k_left: %d | k_right: %d\n", v, k_left, k_right);
	    }

	  if(k_left != -1)
	    is_left = TRUE;
	  if(k_right != -1)
	    is_right = TRUE;

	  if(is_left)
	    {
	      if(ks_left == HMMMATCH) /* match */
		{
		  imin[v] = pn_min_m[k_left];
		  imax[v] = pn_max_m[k_left];
		}		
	      if(ks_left == HMMINSERT) /* insert */
		{
		  imin[v] = pn_min_i[k_left];
		  imax[v] = pn_max_i[k_left];
		}		
	      if(ks_left == HMMDELETE) /* delete */
		{
		  imin[v] = pn_min_d[k_left];
		  imax[v] = pn_max_d[k_left];
		}		
	      if(left_unknown) 
		{
		  /* we're going to  need to fill in imin and imax values 
		   * up to the nearest E state, but not yet, wait til we
		   * reach a new node, so we can be sure of safe imin and imax values.
		   */
		  set_left_known = TRUE;
		  prev_imin = imin[v];
		  prev_imax = imax[v];
		}
	      else if(at_new_node)
		{
		  prev_imin = imin[v];
		  prev_imax = imax[v];
		}
	      else
		{	      
		  if(imin[v] < prev_imin)
		    prev_imin = imin[v];
		  if(imax[v] > prev_imax)
		    prev_imax = imax[v];
		}
	    }	  
	  else
	    {
	      //printf("setting imin[%d] to prev_imin%d\n", imin[v], prev_imin);
	      //printf("setting imax[%d] to prev_imax%d\n", imax[v], prev_imax);
	      imin[v] = prev_imin;
	      imax[v] = prev_imax;
	    }

	  if(is_right)
	    {
	      if(ks_right == HMMMATCH) /* match */
		{
		  jmin[v] = pn_min_m[k_right];
		  jmax[v] = pn_max_m[k_right];
		}		
	      if(ks_right == HMMINSERT) /* insert */
		{
		  jmin[v] = pn_min_i[k_right];
		  jmax[v] = pn_max_i[k_right];
		}		
	      if(ks_right == HMMDELETE) /* delete */
		{
		  jmin[v] = pn_min_d[k_right];
		  jmax[v] = pn_max_d[k_right];
		  //printf("set v: %d | right delete | jmin[v]: %d | jmax[v]: %d\n", v, jmin[v], jmax[v]);
		}		
	      if(right_unknown) 
		{
		  /* we're going to  need to fill in jmin and jmax values 
		   * up to the nearest E state, but not yet, wait til we
		   * reach a new node, so we can be sure of safe imin and imax values.
		   */
		  set_right_known = TRUE;
		  prev_jmin = jmin[v];
		  prev_jmax = jmax[v];
		}
	      else if(at_new_node) 
		{
		  prev_jmin = jmin[v];
		  prev_jmax = jmax[v];
		}
	      else
		{
		  if(jmin[v] < prev_jmin)
		    prev_jmin = jmin[v];
		  if(jmax[v] > prev_jmax)
		    prev_jmax = jmax[v];
		}
	    }	  
	  else
	    {
	      jmin[v] = prev_jmin;
	      jmax[v] = prev_jmax;
	    }

	  if(!(is_left) && !(is_right)) /* v is a B, S or E state, which doesn't map to anything */
	    {
	      if(cm->sttype[v] != B_st && cm->sttype[v] != E_st && cm->sttype[v] != S_st)
		{
		  printf("ERROR: v: %d not B, S or E, but doesn't map to HMM. Exiting...\n", v);
		  exit(1);
		}
	      imin[v] = prev_imin;
	      imax[v] = prev_imax;
	      jmin[v] = prev_jmin;
	      jmax[v] = prev_jmax;
	      //printf("no map for v: %d | imin: %d | imax: %d | jmin: %d | jmax: %d\n", v, imin[v], imax[v], jmin[v], jmax[v]);
	    }
	}
    }      
}

#endif


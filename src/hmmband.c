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
AllocCP9Bands(CM_t *cm, CP9_t *hmm)
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
 *           fmx          - CP9 dp matrix for Forward()
 *           bmx          - CP9 dp matrix for Backward()
 *           pmx          - CP9 dp matrix to fill with posteriors, can == bmx
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
cp9_Seq2Bands(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, int doing_search, int debug_level)
{
  int             status;
  int             use_sums; /* TRUE to fill and use posterior sums during HMM band calc, yields wider bands  */
  float           sc;
  int do_scan2bands;         /* TRUE to use scanning Forward/Backward to get posteriors */

  /* Contract checks */
  if(cm->cp9 == NULL)    ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Bands, but cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Bands, but cm->cp9map is NULL.\n");
  if(dsq == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Bands, dsq is NULL.");
  if(i0 > j0)            ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Bands, i0: %d > j0: %d\n", i0, j0);
  if((cm->align_opts & CM_ALIGN_HBANDED) && (cm->search_opts & CM_SEARCH_HBANDED))           ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Bands, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both up, exactly 1 must be up.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED)))        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Bands, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down, exactly 1 must be up.\n");
  if((cm->search_opts & CM_SEARCH_HMMSCANBANDS) && (!(cm->search_opts & CM_SEARCH_HBANDED))) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Bands, CM_SEARCH_HMMSCANBANDS flag raised, but not CM_SEARCH_HBANDED flag, this doesn't make sense\n");
  
  use_sums = ((cm->align_opts & CM_ALIGN_SUMS) || (cm->search_opts & CM_SEARCH_SUMS)) ? TRUE : FALSE;
    
  /* Step 1: Get HMM Forward/Backward DP matrices.
   * Step 2: F/B       -> HMM bands.
   * Step 3: HMM bands -> CM bands.
   */

  /* Step 1: Get HMM Forward/Backward DP matrices. */


  do_scan2bands = (cm->search_opts & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;
  if((status = cp9_Forward(cm, errbuf, fmx, dsq, i0, j0, j0-i0+1, 0., NULL,
			   do_scan2bands, /* are we using scanning Forward/Backward */
			   TRUE,      /* we are going to use posteriors to align */
			   FALSE,     /* don't be memory efficient */
			   NULL, NULL,
			   &sc)) != eslOK) return status;
  if((status = cp9_Backward(cm, errbuf, bmx, dsq, i0, j0, (j0-i0+1), 0, NULL, 
			    do_scan2bands, /* are we using scanning Forward/Backward */
			    TRUE,  /* we are going to use posteriors to align */
			    FALSE, /* don't be memory efficient */
			    NULL, NULL,
			    &sc)) != eslOK) return status;

  if(cm->align_opts & CM_ALIGN_CHECKFB) { 
    if((status = cp9_CheckFB(fmx, bmx, cm->cp9, errbuf, sc, i0, j0, dsq)) != eslOK) return status;
    printf("Forward/Backward matrices checked.\n");
  }


  /* Step 2: F/B -> HMM bands. */
  if(use_sums){
    if((status = cp9_FB2HMMBandsWithSums(cm->cp9, errbuf, dsq, fmx, bmx, pmx, cp9b, i0, j0, cp9b->hmm_M,
					 (1.-cm->tau), (cm->search_opts & CM_SEARCH_HMMSCANBANDS), debug_level)) != eslOK) return status;
  }
  else {
    if((status = cp9_FB2HMMBands(cm->cp9, errbuf, dsq, fmx, bmx, pmx, cp9b, i0, j0, cp9b->hmm_M,
				 (1.-cm->tau), (cm->search_opts & CM_SEARCH_HMMSCANBANDS), debug_level)) != eslOK) return status;
  }
  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, j0, cp9b, cm->tau, 1);

  /* Step 3: HMM bands  ->  CM bands. */
  if((status = cp9_HMM2ijBands(cm, errbuf, cm->cp9b, cm->cp9map, i0, j0, debug_level)) != eslOK) return status;
  
  /* Use the CM bands on i and j to get bands on d, specific to j. */
  if(doing_search) cp9_RelaxRootBandsForSearch(cm, cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax);
  /* cp9_GrowHDBands() must be called before ij2d_bands() so hdmin, hdmax are adjusted for new seq */
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) return status; 
  ij2d_bands(cm, (j0-i0+1), cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, debug_level);
  
#if eslDEBUGLEVEL >= 1
  if((status = cp9_ValidateBands(cm, errbuf, cp9b, i0, j0)) != eslOK) return status;
#endif

  if(debug_level > 0) PrintDPCellsSaved_jd(cm, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, (j0-i0+1));

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
 *           fmx          - CP9 dp matrix for Forward()
 *           bmx          - CP9 dp matrix for Backward()
 *           pmx          - CP9 dp matrix to fill with posteriors, can == bmx
 *           dsq          - sequence in digitized form
 *           i0           - start of target subsequence (often 1, beginning of dsq)
 *           j0           - end of target subsequence (often L, end of dsq)
 *           debug_level  - verbosity level for debugging printf()s
 *           
 * Return:  eslOK on success
 */
int
cp9_Seq2Posteriors(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, int debug_level)
{
  /*CP9_MX *cp9_mx;*/    /* growable DP matrix for viterbi                       */
  int status;
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
  if((cm->search_opts & CM_SEARCH_HMMSCANBANDS) && (! (cm->search_opts & CM_SEARCH_HBANDED))) 
    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, CM_SEARCH_HMMSCANBANDS flag raised, but not CM_SEARCH_HBANDED flag, this doesn't make sense\n");

  do_scan2bands = (cm->search_opts & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;

  /* Step 1: Get HMM posteriors.*/
  if((status = cp9_FastForward(cm, errbuf, fmx, dsq, i0, j0, j0-i0+1, 0., NULL,
			       do_scan2bands, /* are we using scanning Forward/Backward */
				TRUE,      /* we are going to use posteriors to align */
			       FALSE,     /* don't be memory efficient */
			       be_safe,   /* can we accelerate w/ no -INFTY logsum funcs? */
			       NULL, NULL,
			       &sc)) != eslOK) return status;
  if(debug_level > 0) printf("CP9 Forward  score : %.4f\n", sc);
  if((status = cp9_Backward(cm, errbuf, bmx, dsq, i0, j0, (j0-i0+1), 0, NULL, 
			    do_scan2bands, /* are we using scanning Forward/Backward */
			    TRUE,  /* we are going to use posteriors to align */
			    FALSE, /* don't be memory efficient */
			    NULL, NULL,
			    &sc)) != eslOK) return status;
  if(debug_level > 0) printf("CP9 Backward  score : %.4f\n", sc);

  if(cm->align_opts & CM_ALIGN_CHECKFB) {
    if((status = cp9_CheckFB(fmx, bmx, cm->cp9, errbuf, sc, i0, j0, dsq)) != eslOK) return status;
    printf("Forward/Backward matrices checked.\n");
  }

  /* Get posteriors */
  cp9_Posterior(dsq, i0, j0, cm->cp9, fmx, bmx, pmx, do_scan2bands);

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
 * CP9_MX fmx:      forward DP matrix, already calc'ed
 * CP9_MX bmx:      backward DP matrix, already calc'ed
 * CP9_MX pmx:      DP matrix for posteriors, filled here, can == bmx
 * dsq              the digitized sequence
 * CP9Bands_t cp9b  CP9 bands data structure
 * int i0           start of target subsequence (often 1, beginning of dsq)
 * int j0           end of target subsequence (often L, end of dsq)
 * int   M          number of nodes in HMM (num columns of pmx matrix)
 * double p_thresh  the probability mass we're requiring is within each band
 * int did_scan     TRUE if Forward/Backward were run in 'scan mode'
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 * 
 * Returns: eslOK on success;
 */
int
cp9_FB2HMMBands(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b, 
		 int i0, int j0, int M, double p_thresh, int did_scan, int debug_level)
{
  int status;
  int k;                                  /* counter over nodes of the model */
  int L = j0-i0+1;                        /* length of sequence */
  int thresh = Prob2Score(((1. - p_thresh)/2.), 1.); /* allowable prob mass excluded on each side */
  int max;                                /* temporary max value */
  int pnmax;                              /* position that gives max */

  /* *_m = match, *_i = insert, *_d = delete */
  int *kthresh_m, *kthresh_i, *kthresh_d; /* [0..k..hmm->M], individual thresholds for each state */
  int *nset_m, *nset_i, *nset_d;          /* [0..k..hmm->M], has minimum been set for this state? */
  int *xset_m, *xset_i, *xset_d;          /* [0..k..hmm->M], has maximum been set for this state? */
  int *mass_m, *mass_i, *mass_d;          /* [0..k..hmm->M], summed log prob of pmx->mx[i][k] from 0..k or k..L */
  int i, ip;                              /* actual position and relative position in sequence, ip = i-i0+1 */
  int sc;                                 /* summed score of all parses (derived from backward matrix) 
					   * if(cm->search_opts & CM_SEARCH_HMMSCANBANDS) Forward and Backward
					   * were run in 'scan mode' where each residue can be begin/end of a parse,
					   * so we have to sum up parses that end at each posn, 
					   * if ! (cm->search_opts & CM_SEARCH_HMMSCANBANDS) we know we have 
					   * to start at residue i0 and end at residue j0, so sc is simply bmx->mmx[0][0]
					   */
  if(bmx != pmx) GrowCP9Matrix(pmx, errbuf, L, M, NULL, NULL, NULL, NULL, NULL);

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
  pmx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc; /* fmx->mmx[0][0] is 0, bmx->mmx[1][0] is overall score */
  pmx->imx[0][0] = -INFTY; /*need seq to get here*/
  pmx->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  if((mass_m[0] = pmx->mmx[0][0]) > thresh) { 
    cp9b->pn_min_m[0] = 0; /* "off-by-one" (from comment * above) is dealt with special at end of function for M_0 */
    nset_m[0] = TRUE; 
  }
  mass_i[0] = -INFTY; /* b/c pmx->imx[0][0] is -INFTY, set above */
  mass_d[0] = -INFTY; /* b/c pmx->dmx[0][0] is -INFTY, set above */

  for (k = 1; k <= M; k++) {
    pmx->mmx[0][k] = -INFTY; /*need seq to get here*/
    pmx->imx[0][k] = -INFTY; /*need seq to get here*/
    pmx->dmx[0][k] = fmx->dmx[0][k] + bmx->dmx[0][k] - sc;
    /* mass_m[k] doesn't change b/c pmx->mmx[0][k] is -INFTY */
    /* mass_i[k] doesn't change b/c pmx->mmx[0][k] is -INFTY */
    if((mass_d[k] = pmx->dmx[0][k]) > thresh) { 
      cp9b->pn_min_d[k] = 0+1; /* see "off-by-one" comment * above */
      nset_d[k] = TRUE; 
    }
  }
  /* Find minimum position in band for each state (M,I,D) of each node (0..M) */
  for (ip = 1; ip <= L; ip++) /* ip is the relative position in the seq */
    {
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      pmx->mmx[ip][0] = -INFTY; /* M_0 doesn't emit */
      pmx->imx[ip][0] = ESL_MAX(fmx->imx[ip][0] + bmx->imx[ip][0] - hmm->isc[dsq[i]][0] - sc, -INFTY);
      /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
      pmx->dmx[ip][0] = -INFTY; /* D_0 doesn't exist */

      for(k = 1; k <= M; k++)
	{
	  pmx->mmx[ip][k] = ESL_MAX(fmx->mmx[ip][k] + bmx->mmx[ip][k] - hmm->msc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  pmx->imx[ip][k] = ESL_MAX(fmx->imx[ip][k] + bmx->imx[ip][k] - hmm->isc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  pmx->dmx[ip][k] = ESL_MAX(fmx->dmx[ip][k] + bmx->dmx[ip][k] - sc, -INFTY);

	  if(! nset_m[k]) { 
	    if((mass_m[k] = ILogsum(mass_m[k], pmx->mmx[ip][k])) > thresh) { 
	      cp9b->pn_min_m[k] = i;
	      nset_m[k] = TRUE; 
	    }
	  }
	  if(! nset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], pmx->imx[ip][k])) > thresh) { 
	      cp9b->pn_min_i[k] = i;
	      nset_i[k] = TRUE; 
	    }
	  }
	  if(! nset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[ip][k])) > thresh) { 
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
	    if((mass_m[k] = ILogsum(mass_m[k], pmx->mmx[ip][k])) > thresh) { 
	      cp9b->pn_max_m[k] = i;
	      xset_m[k] = TRUE; 
	    }
	  }
	  if(! xset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], pmx->imx[ip][k])) > thresh) { 
	      cp9b->pn_max_i[k] = i;
	      xset_i[k] = TRUE; 
	    }
	  }
	  if(! xset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[ip][k])) > thresh) { 
	      cp9b->pn_max_d[k] = i+1; /* see "off-by-one" comment * above */
	      xset_d[k] = TRUE; 
	    }
	  }
	}
    }	  
  /* note boundary conditions, ip = 0, i = i0-1 */
  if(! xset_m[0]) { 
    if((mass_m[0] = ILogsum(mass_m[0], pmx->mmx[0][0])) > thresh) { 
      cp9b->pn_max_m[0] = 0; /* "off-by-one" (from comment * above) is dealt with special at end of function for M_0 */
      xset_m[0] = TRUE; 
    }
  }
  /* mass_i[0] is unchaged because b/c pmx->imx[0][0] is -INFTY, set above */
  /* mass_d[0] is unchaged because b/c pmx->dmx[0][0] is -INFTY, set above */
  for (k = 1; k <= M; k++) {
    /* mass_m[k] doesn't change b/c pmx->mmx[0][k] is -INFTY */
    /* mass_i[k] doesn't change b/c pmx->mmx[0][k] is -INFTY */
    if(!xset_d[k]) { 
      if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[0][k])) > thresh) { 
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
	max = pmx->mmx[0][k];
	for(ip = 1; ip <= L; ip++)
	  if(pmx->mmx[ip][k] > max) { pnmax = i; max = pmx->mmx[ip][k]; }
	cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = pnmax;
      }
      if((! nset_i[k]) || (! xset_i[k])) { 
	max = pmx->imx[0][k];
	for(ip = 1; ip <= L; ip++)
	  if(pmx->imx[ip][k] > max) { pnmax = i; max = pmx->imx[ip][k]; }
	cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = pnmax;
      }
      if((! nset_d[k]) || (! xset_d[k])) { 
	max = pmx->dmx[0][k];
	for(ip = 1; ip <= L; ip++)
	  if(pmx->dmx[ip][k] > max) { pnmax = i; max = pmx->dmx[ip][k]; }
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
  free(kthresh_m);
  free(kthresh_i);
  free(kthresh_d);

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
 * CP9_MX fmx:      forward DP matrix, already calc'ed
 * CP9_MX bmx:      backward DP matrix, already calc'ed
 * CP9_MX pmx:      DP matrix for posteriors, filled here, can == bmx
 * dsq              the digitized sequence
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
cp9_FB2HMMBandsWithSums(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b, 
			int i0, int j0, int M, double p_thresh, int did_scan, int debug_level)
{
  int status;
  int k;                                  /* counter over nodes of the model */
  int L = j0-i0+1;                        /* length of sequence */
  int thresh = Prob2Score(((1. - p_thresh)/2.), 1.); /* allowable prob mass excluded on each side */

  /* *_m = match, *_i = insert, *_d = delete */
  int i, ip;                              /* actual position and relative position in sequence, ip = i-i0+1 */
  int *kthresh_m, *kthresh_i, *kthresh_d; /* [0..k..hmm->M], individual thresholds for each state */
  int *nset_m, *nset_i, *nset_d;          /* [0..k..hmm->M], has minimum been set for this state? */
  int *xset_m, *xset_i, *xset_d;          /* [0..k..hmm->M], has maximum been set for this state? */
  int *mass_m, *mass_i, *mass_d;          /* [0..k..hmm->M], summed log prob of pmx->mx[i][k] from 0..k or k..L */
  
  if(bmx != pmx) GrowCP9Matrix(pmx, errbuf, L, M, NULL, NULL, NULL, NULL, NULL);

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
  cp9_Posterior(dsq, i0, j0, hmm, fmx, bmx, pmx, did_scan);

  /* fill ipost_sums in cp9bands data structure */
  cp9_IFillPostSums(pmx, cp9b, i0, j0);

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
	    if((mass_m[k] = ILogsum(mass_m[k], pmx->mmx[ip][k])) > kthresh_m[k]) { 
	      cp9b->pn_min_m[k] = i;
	      nset_m[k] = TRUE; 
	    }
	  }
	  if(! nset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], pmx->imx[ip][k])) > kthresh_i[k]) { 
	      cp9b->pn_min_i[k] = i;
	      nset_i[k] = TRUE; 
	    }
	  }
	  if(! nset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[ip][k])) > kthresh_d[k]) { 
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
	    if((mass_m[k] = ILogsum(mass_m[k], pmx->mmx[ip][k])) > kthresh_m[k]) { 
	      cp9b->pn_max_m[k] = i;
	      xset_m[k] = TRUE; 
	    }
	  }
	  if(! xset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], pmx->imx[ip][k])) > kthresh_i[k]) { 
	      cp9b->pn_max_i[k] = i;
	      xset_i[k] = TRUE; 
	    }
	  }
	  if(! xset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[ip][k])) > kthresh_d[k]) { 
	      cp9b->pn_max_d[k] = i+1; /* see "off-by-one" comment * above */
	      xset_d[k] = TRUE; 
	    }
	  }
	}
    }	  

#if eslDEBUGLEVEL >= 1
  /* all states should have their min/max set because we've normalized the probability
   * of entering each state to 1.0, so we assert this to be true */
  ESL_DASSERT1((nset_m[0]));
  ESL_DASSERT1((nset_i[0]));
  ESL_DASSERT1((xset_m[0]));
  ESL_DASSERT1((xset_i[0]));
  /* D_0 state does not exist */
  for(k = 1; k <= M; k++)
    {
      ESL_DASSERT1((nset_m[k]));
      ESL_DASSERT1((nset_i[k]));
      ESL_DASSERT1((nset_d[k]));
      ESL_DASSERT1((xset_m[k]));
      ESL_DASSERT1((xset_i[k]));
      ESL_DASSERT1((xset_d[k]));
    }
#endif

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

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
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
	      CP9_t *hmm,
	      CP9_MX *fmx,
	      CP9_MX *bmx,
	      CP9_MX *mx,
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
cp9_IFillPostSums(CP9_MX *post, CP9Bands_t *cp9b, int i0, int j0)
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
 * CM_t *cm         the CM, must have valid cp9b (CP9 bands object)
 * errbuf           char buffer for error messages
 * CP9Bands_t *cp9b the CP9 bands object, usually cm->cp9b
 * CP9Map_t *cp9map map from CM to CP9 HMM and vice versa
 * int i0           start of target subsequence (often 1, beginning of dsq)
 * int j0           end of target subsequence (often L, end of dsq)
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *
 * Returns: eslOK on success;
 */
int
cp9_HMM2ijBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int debug_level)
{
  int v;              /* counter over states of the CM */

  int status;
  int safe_imax; 
  int safe_jmin; 

  int tmp_imin;
  int tmp_imax;
  int tmp_jmin;
  int tmp_jmax;

  /* ptrs to cp9b data, for convenience */
  int *pn_min_m;      /* pn_min_m[k] = first position in HMM band for match state of HMM node k */
  int *pn_max_m;      /* pn_max_m[k] = last position in HMM band for match state of HMM node k */
  int *pn_min_i;      /* pn_min_i[k] = first position in HMM band for insert state of HMM node k */
  int *pn_max_i;      /* pn_max_i[k] = last position in HMM band for insert state of HMM node k */
  int *pn_min_d;      /* pn_min_d[k] = first position in HMM band for delete state of HMM node k */
  int *pn_max_d;      /* pn_max_d[k] = last position in HMM band for delete state of HMM node k */
  int *imin;          /* imin[v] = first position in band on i for state v to be filled in this function. [1..M] */
  int *imax;          /* imax[v] = last position in band on i for state v to be filled in this function. [1..M] */
  int *jmin;          /* jmin[v] = first position in band on j for state v to be filled in this function. [1..M] */
  int *jmax;          /* jmax[v] = last position in band on j for state v to be filled in this function. [1..M] */
  
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

  int n;            /* counter over CM nodes. */
  int doing_search; /* TRUE if bands will be used to accelerate search, FALSE for alignment */

  /* Contract checks */
  if (cp9b == NULL)                                                                   ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), cp9b is NULL.\n");
  if((cm->align_opts & CM_ALIGN_HBANDED) && (cm->search_opts & CM_SEARCH_HBANDED))    ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both up, exactly 1 must be up.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED))) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down, exactly 1 must be up.\n");

  doing_search = (cm->search_opts & CM_SEARCH_HBANDED) ? TRUE : FALSE;

  /* set pointers to cp9b data
   * note: these arrays used to be allocated here, but that was wasteful, now it's allocated
   * once per model (instead of once per sequence) in AllocCP9Bands()  
   */

  pn_min_m = cp9b->pn_min_m;
  pn_max_m = cp9b->pn_max_m;
  pn_min_i = cp9b->pn_min_i;
  pn_max_i = cp9b->pn_max_i;
  pn_min_d = cp9b->pn_min_d;
  pn_max_d = cp9b->pn_max_d;
  imin     = cp9b->imin;
  imax     = cp9b->imax;
  jmin     = cp9b->jmin;
  jmax     = cp9b->jmax;

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

  esl_vec_ISet(nss_imin, cm->nodes, -1);
  esl_vec_ISet(nss_imax, cm->nodes, -1);
  esl_vec_ISet(nss_jmin, cm->nodes, -1);
  esl_vec_ISet(nss_jmax, cm->nodes, -1);

  esl_vec_ISet(nis_imin, cm->nodes, -1);
  esl_vec_ISet(nis_imax, cm->nodes, -1);
  esl_vec_ISet(nis_jmin, cm->nodes, -1);
  esl_vec_ISet(nis_jmax, cm->nodes, -1);

  esl_vec_ISet(nss_max_imin, cm->nodes, -1);
  esl_vec_ISet(nss_min_jmax, cm->nodes, -1);

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

    /* Ensure: for all v imin[v] >= imin[0],
     *                   jmax[v] <= jmax[0].
     */
    imin[v] = ESL_MAX(imin[v], imin[0]);
    imax[v] = ESL_MAX(imax[v], imin[0]);
    jmax[v] = ESL_MIN(jmax[v], jmax[0]);
    jmin[v] = ESL_MIN(jmin[v], jmax[0]);

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


/*****************************************************************************
 * EPN 11.03.05
 * Function: ij2d_bands()
 *
 * Purpose:  Determine the band for each cm state v on d (the band on the 
 *           length of the subsequence emitted from the subtree rooted
 *           at state v). These are easily calculated given the bands on i
 *           and j.
 * 
 * arguments:
 *
 * CM_t *cm         the CM 
 * int  W           length of sequence we're aligning
 * int *imin        imin[v] = first position in band on i for state v
 * int *imax        imax[v] = last position in band on i for state v
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 * int **hdmin      hdmin[v][j0] = first position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Filled in this function.
 * int **hdmax      hdmax[v][j0] = last position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Filled in this function.
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *****************************************************************************/
void
ij2d_bands(CM_t *cm, int W, int *imin, int *imax, int *jmin, int *jmax,
	   int **hdmin, int **hdmax, int debug_level)
{
  int v;            /* counter over states of the CM */
  int j0;           /* counter over valid j's, but offset. j0+jmin[v] = actual j */
  int j;            /* actual j */
  int sd;           /* minimum d allowed for a state, ex: MP_st = 2, ML_st = 1. etc. */
  for(v = 0; v < cm->M; v++) {
    if(cm->sttype[v] == E_st) {
      for(j0 = 0; j0 <= (jmax[v]-jmin[v]); j0++) {
	hdmin[v][j0] = 0;
	hdmax[v][j0] = 0;
      }
    }
    else {
      sd = StateDelta(cm->sttype[v]);
      for(j0 = 0; j0 <= (jmax[v]-jmin[v]); j0++) {
	j = j0+jmin[v];
	hdmin[v][j0] = ESL_MAX((j - imax[v] + 1), sd);
	hdmax[v][j0] = ESL_MAX((j - imin[v] + 1), sd);
	/*printf("hd[%d][j=%d]: min: %d | max: %d\n", v, (j0+jmin[v]), hdmin[v][j0], hdmax[v][j0]);*/
      }
    }
  }
}

/*****************************************************************************
 * EPN, Thu Apr 26 13:27:16 2007
 * Function: combine_qdb_hmm_d_bands()
 *
 * Purpose:  Given hdmin and hdmax bands, and query dependent bands (QDBs)
 *           in cm->dmin and cm->dmax, combine them by redefining the 
 *           hdmin and hdmax bands where necessary:
 *           hdmin[v][j] = max(hdmin[v][j], dmin[v])
 *           hdmax[v][j] = min(hdmin[v][j], dmin[v])
 * 
 * arguments:
 *
 * CM_t *cm         the CM 
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 * int **hdmin      hdmin[v][j0] = first position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Redefined in this function.
 * int **hdmax      hdmax[v][j0] = last position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Redefined in this function.
 *****************************************************************************/
void
combine_qdb_hmm_d_bands(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int v;            /* counter over states of the CM */
  int jp;           /* counter over valid j's, but offset. jp+jmin[v] = actual j */

  /* Contract check */
  if(!(cm->flags & CMH_QDB))
    esl_fatal("ERROR, in combine_qdb_hmm_d_bands(), CM QDBs invalid.\n");
  if(cm->dmin == NULL || cm->dmax == NULL)
    esl_fatal("ERROR, in combine_qdb_hmm_d_bands() but cm->dmin and/or cm->dmax is NULL.\n");

  for(v = 0; v < cm->M; v++)
    {
      for(jp = 0; jp <= (jmax[v]-jmin[v]); jp++)
	{
	  hdmin[v][jp] = hdmin[v][jp] > cm->dmin[v] ? hdmin[v][jp] : cm->dmin[v];
	  hdmax[v][jp] = hdmax[v][jp] < cm->dmax[v] ? hdmax[v][jp] : cm->dmax[v];
	}
    }
}


/*****************************************************************************
 * EPN 11.04.05
 * Function: hd2safe_hd_bands
 *
 * Purpose:  HMMERNAL milestone 4 function. Given 
 *           hdmin and hdmax 2D arrays, simply
 *           fill safe_hdmin and safe_hdmax (1D arrays):
 *           safe_hdmin[v] = min_d (hdmin[v][j0])
 *           safe_hdmax[v] = max_d (hdmax[v][j0])
 * 
 * arguments:
 * int M            num states in the CM.
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 * int **hdmin      hdmin[v][j0] = first position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 * int **hdmax      hdmax[v][j0] = last position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 * int *safe_hdmin  safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0)
 *                  filled in this function.
 * int *safe_hdmax  safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0)
 *                  filled in this function.
 *****************************************************************************/
void
hd2safe_hd_bands(int M, int *jmin, int *jmax, int **hdmin, int **hdmax,
		 int *safe_hdmin, int *safe_hdmax)

{
  int v;            /* counter over states of the CM */
  int j0;           /* counter over valid j's, but offset. j0+jmin[v] = actual j */

  for(v = 0; v < M; v++)
    {
      safe_hdmin[v] = hdmin[v][0];
      safe_hdmax[v] = hdmax[v][0];
      /*printf("j0: %2d | j: %2d | v: %3d | smin %d | smax : %d\n", 0, (jmin[v]), v, safe_hdmin[v], safe_hdmax[v]);*/
      for(j0 = 1; j0 <= (jmax[v]-jmin[v]); j0++)
	{
	  if(hdmin[v][j0] < safe_hdmin[v])
	    safe_hdmin[v] = hdmin[v][j0];
	  if(hdmax[v][j0] > safe_hdmax[v])
	    safe_hdmax[v] = hdmax[v][j0];
	  /*printf("j0: %2d | j: %2d | v: %3d | smin %d | smax : %d\n", j0, (j0+jmin[v]), v, safe_hdmin[v], safe_hdmax[v]);*/
	}
    }

  for(v = 0; v < M; v++)
    {
      for(j0 = 1; j0 <= (jmax[v]-jmin[v]); j0++)
	{
	  /* check to make sure that hdmax[v][j0] doesn't exceed j */
	  if(hdmax[v][j0] > (j0+jmin[v]))
	    {
	      printf("ERROR: hd2safe_bands.1\n");
	      /*exit(1);*/
	    }
	}
    }
}

/*****************************************************************************
 * EPN 01.18.06
 * Function: debug_print_hd_bands
 *
 * Purpose:  Print out the v and j dependent hd bands.
 * 
 *****************************************************************************/
void
debug_print_hd_bands(CM_t *cm, int **hdmin, int **hdmax, int *jmin, int *jmax)
{
  int status;
  int v, j;
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

  printf("\nPrinting hd bands :\n");
  printf("****************\n");
  for(v = 0; v < cm->M; v++)
   {
     for(j = jmin[v]; j <= jmax[v]; j++) 
       {
	 printf("band v:%d j:%d n:%d %-4s %-2s min:%d max:%d\n", v, j, cm->ndidx[v], nodetypes[(int) cm->ndtype[(int) cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], hdmin[v][j-jmin[v]], hdmax[v][j-jmin[v]]);
       }
     printf("\n");
   }
  printf("****************\n\n");

  free(sttypes);
  free(nodetypes);
  return;

 ERROR:
  esl_fatal("Memory allocation error.");
}


/* Function: PrintDPCellsSaved_jd()
 * Prints out an estimate of the speed up due to j and d bands */
void
PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax,
		     int W)
{
  int v;
  int j;
  int max;
  double after, before;

  printf("Printing DP cells saved using j and d bands:\n");
  before = after = 0;
  for (v = 0; v < cm->M; v++) 
    {
      for(j = 0; j <= W; j++)
	if (cm->sttype[v] != E_st) 
	  before += j + 1;
      for(j = jmin[v]; j <= jmax[v]; j++)
	if (cm->sttype[v] != E_st) 
	  {
	    max = (j < hdmax[v][j-jmin[v]]) ? j : hdmax[v][j-jmin[v]];
	    after += max - hdmin[v][j-jmin[v]] + 1;
	  }
    }
  printf("Before:  something like %.0f\n", before);
  printf("After:   something like %.0f\n", after);
  printf("Speedup: maybe %.2f fold\n\n", (float) before / (float) after);
}

/* old_miscfuncs.c
 * EPN, Wed Dec  5 10:35:55 2007
 *
 * Options from Infernal codebase as of revision 2243 that
 * are not used in any circumstance. These are typically not
 * old DP functions which can be found in old_cm_dpsearch.c, 
 * old_cm_dpalign.c and old_cp9_dp.c. 
 * 
 * This code is kept solely for reference, and is not
 * compiled/compilable.
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_histogram.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_stats.h"

#include "funcs.h"
#include "structs.h"
  

/***********************************************************************
 * Function: CP9Scan_dispatch()
 * Incept:   EPN, Tue Jan  9 06:28:49 2007
 * 
 * Purpose:  Scan a sequence with a CP9, potentially rescan CP9 hits with CM.
 *
 *           3 possible modes:
 *
 *           Mode 1: Filter mode with default pad 
 *                   (IF cm->search_opts & CM_SEARCH_HMMFILTER and 
 *                      !cm->search_opts & CM_SEARCH_HMMPAD)
 *                   Scan with CP9Forward() to get likely endpoints (j) of 
 *                   hits, for each j do a CP9Backward() from j-W+1..j to get
 *                   the most likely start point i for this j. 
 *                   Set i' = j-W+1 and
 *                       j' = i+W-1. 
 *                   Each i'..j' subsequence is rescanned with the CM.
 * 
 *           Mode 2: Filter mode with user-defined pad 
 *                   (IF cm->search_opts & CM_SEARCH_HMMFILTER and 
 *                       cm->search_opts & CM_SEARCH_HMMPAD)
 *                   Same as mode 1, but i' and j' defined differently:
 *                   Set i' = i - (cm->hmmpad) and
 *                       j' = j + (cm->hmmpad) 
 *                   Each i'..j' subsequence is rescanned with the CM.
 * 
 *           Mode 3: HMM only mode (IF cm->search_opts & CM_SEARCH_HMMONLY)
 *                   Hit boundaries i and j are defined the same as in mode 1, but
 *                   no rescanning with CM takes place. i..j hits are reported 
 *                   (note i' and j' are never calculated).
 * 
 * Args:     
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report (or keep if filtering)
 *           results    - search_results_t to add to, only passed to 
 *                        OldActuallySearchTarget()
 *           doing_cp9_stats- TRUE if we're calc'ing stats for the CP9, in this 
 *                            case we never rescan with CM
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc, score of maximally scoring end position j 
 */
float
CP9Scan_dispatch(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cm_cutoff, 
		 float cp9_cutoff, search_results_t *results, int doing_cp9_stats,
		 int *ret_flen)
{
  int h;
  int i;
  int min_i;
  float best_hmm_sc;
  float best_hmm_fsc;
  float cur_best_hmm_bsc;
  float best_cm_sc;
  int   flen;
  float ffrac;
  int do_collapse;
  int ipad;
  int jpad;
  int padmode;
  search_results_t *fwd_results;
  search_results_t *bwd_results;

  /* check contract */
  if(cm->cp9 == NULL)
    cm_Fail("ERROR in CP9Scan_dispatch(), cm->cp9 is NULL\n");
  if((cm->search_opts & CM_SEARCH_HMMPAD) &&
     (!(cm->search_opts & CM_SEARCH_HMMFILTER)))
     cm_Fail("ERROR in CP9Scan_dispatch(), CM_SEARCH_HMMPAD flag up, but CM_SEARCH_HMMFILTER flag down.\n");
  if(!doing_cp9_stats && (!((cm->search_opts & CM_SEARCH_HMMFILTER) || 
			    (cm->search_opts & CM_SEARCH_HMMONLY))))
    cm_Fail("ERROR in CP9Scan_dispatch(), not doing CP9 stats and neither CM_SEARCH_HMMFILTER nor CM_SEARCH_HMMONLY flag is up.\n");
  if(dsq == NULL)
    cm_Fail("ERROR in CP9Scan_dispatch, dsq is NULL.");

  /*printf("in CP9Scan_dispatch(), i0: %d j0: %d\n", i0, j0);
    printf("cp9_cutoff: %f\n", cp9_cutoff);*/

  best_cm_sc = best_hmm_sc = IMPOSSIBLE;
  /* set up options for RescanFilterSurvivors() if we're filtering */
  if(cm->search_opts & CM_SEARCH_HMMFILTER)
    {
      if(cm->search_opts & CM_SEARCH_HMMPAD) /* mode 2 */
	{
	  padmode = PAD_SUBI_ADDJ;
	  ipad = jpad = cm->hmmpad; /* subtract hmmpad from i, add hmmpad to j */
	}
      else /* mode 1 */
	{
	  padmode = PAD_ADDI_SUBJ;
	  ipad = jpad = W-1; /* subtract W-1 from j, add W-1 to i */
	}
      if(cm->search_opts && CM_SEARCH_HBANDED)
	do_collapse = FALSE;
      else
	do_collapse = TRUE;
    }
  
  /* Scan the (sub)seq w/Forward, getting j end points of hits above cutoff */
  fwd_results = CreateResults(INIT_RESULTS);
  best_hmm_fsc = CP9Forward(cm, dsq, i0, j0, W, cp9_cutoff, NULL, NULL, fwd_results,
			    TRUE,   /* we're scanning */
			    FALSE,  /* we're not ultimately aligning */
			    FALSE,  /* we're not rescanning */
			    TRUE,   /* be memory efficient */
			    NULL);  /* don't want the DP matrix back */
  best_hmm_sc = best_hmm_fsc;

  /* Remove overlapping hits, if we're being greedy */
  if(cm->search_opts & CM_SEARCH_HMMGREEDY) /* resolve overlaps by being greedy */
    {
      assert(i0 == 1); 
      remove_overlapping_hits (fwd_results, i0, j0);
    }

  /* Determine start points (i) of the hits based on Backward scan starting at j,
   * report hits IFF CM_SEARCH_HMMONLY */
  bwd_results = CreateResults(INIT_RESULTS);
  for(h = 0; h < fwd_results->num_results; h++) 
    {
      min_i = (fwd_results->data[h].stop - W + 1) >= 1 ? (fwd_results->data[h].stop - W + 1) : 1;
      cur_best_hmm_bsc = CP9Backward(cm, dsq, min_i, fwd_results->data[h].stop, W, cp9_cutoff, 
				     NULL, /* don't care about score of each posn */
				     &i,   /* set i as the best scoring start point from j-W..j */
				     ((cm->search_opts & CM_SEARCH_HMMONLY) ? results : bwd_results),  
				     TRUE,  /* we're scanning */
				     /*FALSE,*/  /* we're not scanning */
				     FALSE, /* we're not ultimately aligning */
				     FALSE, /* don't rescan */
				     TRUE,  /* be memory efficient */
				     NULL); /* don't want the DP matrix back */
      //FALSE,  /* don't be memory efficient */
      //&bmx); /* give the DP matrix back */
      /* this only works if we've saved the matrices, and didn't do scan mode
       * for both Forward and Backward:
       * debug_check_CP9_FB(fmx, bmx, cm->cp9, cur_best_hmm_bsc, i0, j0, dsq); */
      
      if(cur_best_hmm_bsc > best_hmm_sc) best_hmm_sc = cur_best_hmm_bsc;
      /*printf("cur_best_hmm_bsc: %f\n", cur_best_hmm_bsc);*/
    }	  
  /* Rescan with CM if we're filtering and not doing cp9 stats */
  if(!doing_cp9_stats && (cm->search_opts & CM_SEARCH_HMMFILTER))
    {
      /* Remove overlapping hits, if we're being greedy */
      if(cm->search_opts & CM_SEARCH_HMMGREEDY) 
	{
	  assert(i0 == 1); 
	  remove_overlapping_hits (bwd_results, i0, j0);
	}
      best_cm_sc = RescanFilterSurvivors(cm, dsq, bwd_results, i0, j0, W, 
					 padmode, ipad, jpad, 
					 do_collapse, cm_cutoff, cp9_cutoff, 
					 results, &flen);
      if(flen == 0) ffrac = 100.;
      else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
      /*if(!(cm->search_opts & CM_SEARCH_HMMONLY))
	printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);*/
    }
  FreeResults (fwd_results);
  FreeResults (bwd_results);

  /*printf("in CP9Scan_dispatch, returning best_hmm_sc: %f\n", best_hmm_sc);*/
  if(doing_cp9_stats || cm->search_opts & CM_SEARCH_HMMONLY)
    return best_hmm_sc;
  else
    return best_cm_sc;
}

/* Function: CP9ScanPosterior()
 * based on Ian Holmes' hmmer/src/postprob.c::P7EmitterPosterior()
 *
 * Purpose:  Combines Forward and Backward scanning matrices into a posterior
 *           probability matrix. 
 *
 *           The main difference between this function and CP9Posterior()
 *           in hmmband.c is that this func takes matrices from CP9ForwardBackwardScan()
 *           in which parses are allowed to start and end in any residue.
 *           In CP9Posterior(), the matrices are calc'ed in CP9Forward()
 *           and CP9Backward() which force all parses considered to start at posn
 *           1 and end at L. This means here we have to calculate probability
 *           that each residue from 1 to L is contained in any parse prior
 *           to determining the posterior probability it was emitted from
 *           each state.
 * 
 *           For emitters (match and inserts) the 
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
 * Args:     dsq      - sequence in digitized form
 *           i0       - start of target subsequence
 *           j0       - end of target subsequence 
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           mx       - pre-allocated dynamic programming matrix
 *           
 * Return:   void
 */
void
CP9ScanPosterior(ESL_DSQ *dsq, int i0, int j0,
		     CP9_t *hmm,
		     CP9_MX *fmx,
		     CP9_MX *bmx,
		     CP9_MX *mx)
{
  int i;
  int ip;
  int k;
  int fb_sum; /* tmp value, the probability that the current residue (i) was
	       * visited in any parse */
  /* contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in CP9ScanPosterior(), dsq is NULL.");

  /*printf("\n\nin CP9ScanPosterior() i0: %d, j0: %d\n", i0, j0);*/
  fb_sum = -INFTY;
  for (i = i0-1; i <= j0; i++) 
    {
      ip = i-i0+1; /* ip is relative position in seq, 0..L */
      /*printf("bmx->mmx[i:%d][0]: %d\n", i, bmx->mmx[ip][0]); */
      fb_sum = ILogsum(fb_sum, (bmx->mmx[ip][0])); 
    }
  /*printf("fb_sc: %f\n", Scorify(fb_sum));*/
    /* fb_sum is the probability of all parses */

    /*for(k = 1; k <= hmm->M; k++)*/
  /*{*/
      /*fbsum_ = ILogsum(fmx->mmx[0][k] + bmx->mmx[0][k]))*/; /* residue 0 can't be emitted
								    * but we can start in BEGIN,
								    * before any residues */
      /*fb_sum = ILogsum(fb_sum, (fmx->imx[0][k] + bmx->imx[0][k]))*/; /* these will be all -INFTY */
  /*}*/      

  /* note boundary conditions, case by case by case... */
  mx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - fb_sum;
  mx->imx[0][0] = -INFTY; /*need seq to get here*/
  mx->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  for (k = 1; k <= hmm->M; k++) 
    {
      mx->mmx[0][k] = -INFTY; /*need seq to get here*/
      mx->imx[0][k] = -INFTY; /*need seq to get here*/
      mx->dmx[0][k] = fmx->dmx[0][k] + bmx->dmx[0][k] - fb_sum;
    }
      
  for (i = i0; i <= j0; i++)
    {
      ip = i-i0+1; /* ip is relative position in seq, 0..L */
      /*fb_sum = -INFTY;*/ /* this will be probability of seeing residue i in any parse */
      /*for (k = 0; k <= hmm->M; k++) 
	{
	fb_sum = ILogsum(fb_sum, (fmx->mmx[i][k] + bmx->mmx[i][k] - hmm->msc[dsq[i]][k]));*/
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
      /*fb_sum = ILogsum(fb_sum, (fmx->imx[i][k] + bmx->imx[i][k] - hmm->isc[dsq[i]][k]));*/
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
      /*}*/
      mx->mmx[ip][0] = -INFTY; /*M_0 does not emit*/
      mx->imx[ip][0] = fmx->imx[ip][0] + bmx->imx[ip][0] - hmm->isc[dsq[i]][0] - fb_sum;
      /*hmm->isc[dsq[i]][0] will have been counted in both fmx->imx and bmx->imx*/
      mx->dmx[ip][0] = -INFTY; /*D_0 does not exist*/
      for (k = 1; k <= hmm->M; k++) 
	{
	  mx->mmx[ip][k] = fmx->mmx[ip][k] + bmx->mmx[ip][k] - hmm->msc[dsq[i]][k] - fb_sum;
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  mx->imx[ip][k] = fmx->imx[ip][k] + bmx->imx[ip][k] - hmm->isc[dsq[i]][k] - fb_sum;
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
	  mx->dmx[ip][k] = fmx->dmx[ip][k] + bmx->dmx[ip][k] - fb_sum;
	}	  
    }
  /*
  float temp_sc;
  for (i = i0; i <= j0; i++)
    {
      ip = i-i0+1; 
      for(k = 0; k <= hmm->M; k++)
	{
	  temp_sc = Score2Prob(mx->mmx[ip][k], 1.);
	  if(temp_sc > .0001)
	    printf("mx->mmx[%3d][%3d]: %9d | %8f\n", i, k, mx->mmx[ip][k], temp_sc);
	  temp_sc = Score2Prob(mx->imx[ip][k], 1.);
	  if(temp_sc > .0001)
	    printf("mx->imx[%3d][%3d]: %9d | %8f\n", i, k, mx->imx[ip][k], temp_sc);
	  temp_sc = Score2Prob(mx->dmx[ip][k], 1.);
	  if(temp_sc > .0001)
	    printf("mx->dmx[%3d][%3d]: %9d | %8f\n", i, k, mx->dmx[ip][k], temp_sc);
	}
    }
  */
}

/***********************************************************************
 * Function: CP9Scan_dispatch()
 * Incept:   EPN, Tue Jan  9 06:28:49 2007
 * 
 * Purpose:  Scan a sequence with a CP9, potentially rescan CP9 hits with CM.
 *
 *           3 possible modes:
 *
 *           Mode 1: Filter mode with default pad 
 *                   (IF cm->search_opts & CM_SEARCH_HMMFILTER and 
 *                      !cm->search_opts & CM_SEARCH_HMMPAD)
 *                   Scan with CP9Forward() to get likely endpoints (j) of 
 *                   hits, for each j do a CP9Backward() from j-W+1..j to get
 *                   the most likely start point i for this j. 
 *                   Set i' = j-W+1 and
 *                       j' = i+W-1. 
 *                   Each i'..j' subsequence is rescanned with the CM.
 * 
 *           Mode 2: Filter mode with user-defined pad 
 *                   (IF cm->search_opts & CM_SEARCH_HMMFILTER and 
 *                       cm->search_opts & CM_SEARCH_HMMPAD)
 *                   Same as mode 1, but i' and j' defined differently:
 *                   Set i' = i - (cm->hmmpad) and
 *                       j' = j + (cm->hmmpad) 
 *                   Each i'..j' subsequence is rescanned with the CM.
 * 
 *           Mode 3: HMM only mode (IF cm->search_opts & CM_SEARCH_HMMONLY)
 *                   Hit boundaries i and j are defined the same as in mode 1, but
 *                   no rescanning with CM takes place. i..j hits are reported 
 *                   (note i' and j' are never calculated).
 * 
 * Args:     
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report (or keep if filtering)
 *           results    - search_results_t to add to, only passed to 
 *                        OldActuallySearchTarget()
 *           doing_cp9_stats- TRUE if we're calc'ing stats for the CP9, in this 
 *                            case we never rescan with CM
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc, score of maximally scoring end position j 
 */
float
CP9Scan_dispatch(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cm_cutoff, 
		 float cp9_cutoff, search_results_t *results, int doing_cp9_stats,
		 int *ret_flen)
{
  int h;
  int i;
  int min_i;
  float best_hmm_sc;
  float best_hmm_fsc;
  float cur_best_hmm_bsc;
  float best_cm_sc;
  int   flen;
  float ffrac;
  int do_collapse;
  int ipad;
  int jpad;
  int padmode;
  search_results_t *fwd_results;
  search_results_t *bwd_results;

  /* check contract */
  if(cm->cp9 == NULL)
    cm_Fail("ERROR in CP9Scan_dispatch(), cm->cp9 is NULL\n");
  if((cm->search_opts & CM_SEARCH_HMMPAD) &&
     (!(cm->search_opts & CM_SEARCH_HMMFILTER)))
     cm_Fail("ERROR in CP9Scan_dispatch(), CM_SEARCH_HMMPAD flag up, but CM_SEARCH_HMMFILTER flag down.\n");
  if(!doing_cp9_stats && (!((cm->search_opts & CM_SEARCH_HMMFILTER) || 
			    (cm->search_opts & CM_SEARCH_HMMONLY))))
    cm_Fail("ERROR in CP9Scan_dispatch(), not doing CP9 stats and neither CM_SEARCH_HMMFILTER nor CM_SEARCH_HMMONLY flag is up.\n");
  if(dsq == NULL)
    cm_Fail("ERROR in CP9Scan_dispatch, dsq is NULL.");

  /*printf("in CP9Scan_dispatch(), i0: %d j0: %d\n", i0, j0);
    printf("cp9_cutoff: %f\n", cp9_cutoff);*/

  best_cm_sc = best_hmm_sc = IMPOSSIBLE;
  /* set up options for RescanFilterSurvivors() if we're filtering */
  if(cm->search_opts & CM_SEARCH_HMMFILTER)
    {
      if(cm->search_opts & CM_SEARCH_HMMPAD) /* mode 2 */
	{
	  padmode = PAD_SUBI_ADDJ;
	  ipad = jpad = cm->hmmpad; /* subtract hmmpad from i, add hmmpad to j */
	}
      else /* mode 1 */
	{
	  padmode = PAD_ADDI_SUBJ;
	  ipad = jpad = W-1; /* subtract W-1 from j, add W-1 to i */
	}
      if(cm->search_opts && CM_SEARCH_HBANDED)
	do_collapse = FALSE;
      else
	do_collapse = TRUE;
    }
  
  /* Scan the (sub)seq w/Forward, getting j end points of hits above cutoff */
  fwd_results = CreateResults(INIT_RESULTS);
  best_hmm_fsc = CP9Forward(cm, dsq, i0, j0, W, cp9_cutoff, NULL, NULL, fwd_results,
			    TRUE,   /* we're scanning */
			    FALSE,  /* we're not ultimately aligning */
			    FALSE,  /* we're not rescanning */
			    TRUE,   /* be memory efficient */
			    NULL);  /* don't want the DP matrix back */
  best_hmm_sc = best_hmm_fsc;

  /* Remove overlapping hits, if we're being greedy */
  if(cm->search_opts & CM_SEARCH_HMMGREEDY) /* resolve overlaps by being greedy */
    {
      assert(i0 == 1); 
      remove_overlapping_hits (fwd_results, i0, j0);
    }

  /* Determine start points (i) of the hits based on Backward scan starting at j,
   * report hits IFF CM_SEARCH_HMMONLY */
  bwd_results = CreateResults(INIT_RESULTS);
  for(h = 0; h < fwd_results->num_results; h++) 
    {
      min_i = (fwd_results->data[h].stop - W + 1) >= 1 ? (fwd_results->data[h].stop - W + 1) : 1;
      cur_best_hmm_bsc = CP9Backward(cm, dsq, min_i, fwd_results->data[h].stop, W, cp9_cutoff, 
				     NULL, /* don't care about score of each posn */
				     &i,   /* set i as the best scoring start point from j-W..j */
				     ((cm->search_opts & CM_SEARCH_HMMONLY) ? results : bwd_results),  
				     TRUE,  /* we're scanning */
				     /*FALSE,*/  /* we're not scanning */
				     FALSE, /* we're not ultimately aligning */
				     FALSE, /* don't rescan */
				     TRUE,  /* be memory efficient */
				     NULL); /* don't want the DP matrix back */
      //FALSE,  /* don't be memory efficient */
      //&bmx); /* give the DP matrix back */
      /* this only works if we've saved the matrices, and didn't do scan mode
       * for both Forward and Backward:
       * debug_check_CP9_FB(fmx, bmx, cm->cp9, cur_best_hmm_bsc, i0, j0, dsq); */
      
      if(cur_best_hmm_bsc > best_hmm_sc) best_hmm_sc = cur_best_hmm_bsc;
      /*printf("cur_best_hmm_bsc: %f\n", cur_best_hmm_bsc);*/
    }	  
  /* Rescan with CM if we're filtering and not doing cp9 stats */
  if(!doing_cp9_stats && (cm->search_opts & CM_SEARCH_HMMFILTER))
    {
      /* Remove overlapping hits, if we're being greedy */
      if(cm->search_opts & CM_SEARCH_HMMGREEDY) 
	{
	  assert(i0 == 1); 
	  remove_overlapping_hits (bwd_results, i0, j0);
	}
      best_cm_sc = RescanFilterSurvivors(cm, dsq, bwd_results, i0, j0, W, 
					 padmode, ipad, jpad, 
					 do_collapse, cm_cutoff, cp9_cutoff, 
					 results, &flen);
      if(flen == 0) ffrac = 100.;
      else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
      /*if(!(cm->search_opts & CM_SEARCH_HMMONLY))
	printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);*/
    }
  FreeResults (fwd_results);
  FreeResults (bwd_results);

  /*printf("in CP9Scan_dispatch, returning best_hmm_sc: %f\n", best_hmm_sc);*/
  if(doing_cp9_stats || cm->search_opts & CM_SEARCH_HMMONLY)
    return best_hmm_sc;
  else
    return best_cm_sc;
}

/***********************************************************************
 * Function: RescanFilterSurvivors()
 * Incept:   EPN, Wed Apr 11 05:51:55 2007
 * 
 * Purpose:  Given start and end points of hits that have survived
 *           a CP9 filter, pad a given amount of residues on 
 *           on each side, and rescan with a CM. Optionally,
 *           collapse overlapping subseqs into one larger subseq before
 *           rescanning (we don't always do this b/c we may want to
 *           derive HMM bands for a subseq from a Forward/Backward scan).
 * 
 *           Can be run in 2 modes, depending on input variable padmode:
 *           Mode 1: padmode = PAD_SUBI_ADDJ
 *                   For each i,j pair in hiti, hitj: 
 *                   set i' = i - ipad; and j' = j + jpad, 
 *                   Rescan subseq from i' to j'.
 *           Mode 2: padmode = PAD_ADDI_SUBJ
 *                   For each i,j pair in hiti, hitj: 
 *                   set j' = i + ipad; and i' = j - jpad, 
 *                   ensure j' >= j and i' <= i. 
 *                   Rescan subseq from i' to j'.
 *
 * Args:     
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           hmm_results- info on HMM hits that survived filter
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           padmode    - PAD_SUBI_ADDJ or PAD_ADDI_SUBJ (see above)
 *           ipad       - number of residues to subtract/add from each i 
 *           jpad       - number of residues to add/subtract from each j 
 *           do_collapse- TRUE: collapse overlapping hits (after padding) into 1
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report 
 *           results    - search_results_t to add to, only passed to 
 *                        OldActuallySearchTarget()
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc found when rescanning with CM 
 */
float
RescanFilterSurvivors(CM_t *cm, ESL_DSQ *dsq, search_results_t *hmm_results, int i0, int j0,
		      int W, int padmode, int ipad, int jpad, int do_collapse,
		      float cm_cutoff, float cp9_cutoff, search_results_t *results, int *ret_flen)
{
  int h;
  int i, j;
  float best_cm_sc;
  float cm_sc;
  int   flen;
  int   prev_j = j0;
  int   next_j;
  int   nhits;

  /* check contract */
  if(padmode != PAD_SUBI_ADDJ && padmode != PAD_ADDI_SUBJ)
    ESL_EXCEPTION(eslEINCOMPAT, "can't determine mode.");
  if(dsq == NULL)
    cm_Fail("ERROR in RescanFilterSurvivors(), dsq is NULL.\n");

  best_cm_sc = IMPOSSIBLE;
  flen = 0;

  /*if(padmode == PAD_SUBI_ADDJ)
    printf("in RescanFilterSurvivors(), mode: PAD_SUBI_ADDJ\n");
    else
    printf("in RescanFilterSurvivors(), mode: PAD_ADDI_SUBJ\n");
    printf("\tipad: %d, jpad: %d collapse: %d\n", ipad, jpad, do_collapse);*/

  /* For each hit, add pad according to mode and rescan by calling OldActuallySearchTarget(). 
   * If do_collapse, collapse multiple overlapping hits into 1 before rescanning */
  /* hits should always be sorted by decreasing j, if this is violated - die. */
  nhits = hmm_results->num_results;
  for(h = 0; h < nhits; h++) 
    {
      if(hmm_results->data[h].stop > prev_j) 
	ESL_EXCEPTION(eslEINCOMPAT, "j's not in descending order");

      prev_j = hmm_results->data[h].stop;

      /* add pad */
      if(padmode == PAD_SUBI_ADDJ)
	{
	  i = ((hmm_results->data[h].start - ipad) >= 1)    ? (hmm_results->data[h].start - ipad) : 1;
	  j = ((hmm_results->data[h].stop  + jpad) <= j0)   ? (hmm_results->data[h].stop  + jpad) : j0;
	  if((h+1) < nhits)
	    next_j = ((hmm_results->data[h+1].stop + jpad) <= j0)   ? (hmm_results->data[h+1].stop + jpad) : j0;
	  else
	    next_j = -1;
	}
      else if(padmode == PAD_ADDI_SUBJ)
	{
	  i = ((hmm_results->data[h].stop  - jpad) >= 1)    ? (hmm_results->data[h].stop  - jpad) : 1;
	  j = ((hmm_results->data[h].start + ipad) <= j0)   ? (hmm_results->data[h].start + ipad) : j0;
	  if((h+1) < nhits)
	    next_j = ((hmm_results->data[h+1].start + ipad) <= j0)   ? (hmm_results->data[h+1].start + ipad) : j0;
	  else
	    next_j = -1;
	}
      /*printf("subseq: hit %d i: %d (%d) j: %d (%d)\n", h, i, hmm_results->data[h].start[h], j, hmm_results->data[h].stop[h]);*/

      if(do_collapse) /* collapse multiple hits that overlap after padding on both sides into a single hit */
	{
	  while(((h+1) < nhits) && (next_j >= i))
	    {
	      /* suck in hit */
	      h++;
	      if(padmode == PAD_SUBI_ADDJ)
		{
		  i = ((hmm_results->data[h].start - ipad) >= 1)    ? (hmm_results->data[h].start - ipad) : 1;
		  if((h+1) < nhits)
		    next_j = ((hmm_results->data[h+1].stop + jpad) <= j0)   ? (hmm_results->data[h+1].stop + jpad) : j0;
		  else
		    next_j = -1;
		}
	      else if(padmode == PAD_ADDI_SUBJ)
		{
		  i = ((hmm_results->data[h].stop - jpad) >= 1)    ? (hmm_results->data[h].stop - jpad) : 1;
		  if((h+1) < nhits)
		    next_j = ((hmm_results->data[h+1].start + ipad) <= j0)   ? (hmm_results->data[h+1].start + ipad) : j0;
		  else
		    next_j = -1;
		}
	      /*printf("\tsucked in subseq: hit %d new_i: %d j (still): %d\n", h, i, j);*/
	    }
	}
      /*printf("in RescanFilterSurvivors(): calling OldActuallySearchTarget: %d %d h: %d nhits: %d\n", i, j, h, nhits);*/
      cm_sc =
	OldActuallySearchTarget(cm, dsq, i, j, cm_cutoff, cp9_cutoff,
				results, /* keep results                                 */
				FALSE,   /* don't filter, we already have                */
				FALSE,   /* we're not building a histogram for CM stats  */
				FALSE,   /* we're not building a histogram for CP9 stats */
				NULL,    /* filter fraction N/A                          */
				FALSE);  /* DO NOT align the hits in this recursive call */
      flen += (j - i + 1);
      if(cm_sc > best_cm_sc) best_cm_sc = cm_sc;
    }

  //if(flen == 0) ffrac = 100.;
  //else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
  if(ret_flen != NULL) *ret_flen = flen;
  return best_cm_sc;
}


/*
 * Function: FindCP9FilterThreshold()
 * Incept:   EPN, Wed May  2 10:00:45 2007
 *
 * Purpose:  Sample sequences from a CM and determine the CP9 HMM E-value
 *           threshold necessary to recognize a specified fraction of those
 *           hits. Sequences are sampled from the CM until N with a E-value
 *           better than cm_ecutoff are sampled (those with worse E-values
 *           are rejected). CP9 scans are carried out in either local or
 *           glocal mode depending on hmm_gum_mode. CM is configured in 
 *           local/glocal and sampled seqs are scored in CYK/inside depending
 *           on fthr_mode (4 possibilities). E-values are determined using
 *           lambda from cm->stats, and a recalc'ed mu using database size
 *           of 'db_size'.
 *
 *           If do_fastfil and fthr_mode is local or glocal CYK, parsetree 
 *           scores of emitted sequences are assumed to an optimal CYK scores 
 *           (not nec true). This saves a lot of time b/c no need to
 *           scan emitted seqs, but it's statistically wrong. 
 *           If do_fastfil and fthr_mode is local or glocal inside, 
 *           contract is violated and we Die. Current strategy *outside*
 *           of this function is to copy HMM filtering thresholds from
 *           CYK for Inside cases.
 *
 *           If !do_fastfil this function takes much longer b/c emitted 
 *           parsetree score is not necessarily (a) optimal nor (b) highest 
 *           score returned from a CYK scan (a subseq of the full seq could
 *           score higher even if parsetree was optimal). For Inside, Inside
 *           scan will always be higher than parstree score. This means we have
 *           to scan each emitted seq with CYK or Inside.
 *
 * Args:
 *           cm           - the CM
 *           cmstats      - CM stats object we'll get Gumbel stats from
 *           r            - source of randomness for EmitParsetree()
 *           Fmin         - minimum target fraction of CM hits to detect with CP9
 *           Smin         - minimum filter survival fraction, for lower fractions, we'll
 *                          increase F. 
 *           Starget      - target filter survival fraction, for lower fractions that 
 *                          satisfy Fmin, we'll increase F until Starget is reached
 *           Spad         - fraction of (sc(S) - sc(max(Starget, Smin))) to subtract from sc(S)
 *                          in case 2. 0.0 = fast, 1.0 = safe.
 *           N            - number of sequences to sample from CM better than cm_minsc
 *           use_cm_cutoff- TRUE to only accept CM parses w/E-values above cm_ecutoff
 *           cm_ecutoff   - minimum CM E-value to accept 
 *           db_size      - DB size (nt) to use w/cm_ecutoff to calc CM bit score cutoff 
 *           emit_mode    - CM_LC or CM_GC, CM mode to emit with
 *           fthr_mode    - gives CM search strategy to use, and Gumbel to use
 *           hmm_gum_mode - CP9_L to search with  local HMM (we're filling a fthr->l_eval)
 *                          CP9_G to search with glocal HMM (we're filling a fthr->g_eval)
 *           do_fastfil   - TRUE to use fast method: assume parsetree score
 *                          is optimal CYK score
 *           do_Fstep     - TRUE to step from F towards 1.0 while S < Starget in case 2.
 *           my_rank      - MPI rank, 0 if serial
 *           nproc        - number of processors in MPI rank, 1 if serial
 *           do_mpi       - TRUE if we're doing MPI, FALSE if not
 *           histfile     - root string for histogram files, we'll 4 of them
 *           Rpts_fp      - open file ptr for optimal HMM/CM score pts
 *           ret_F        - the fraction of observed CM hits we've scored with the HMM better
 *                          than return value
 * 
 * Returns: HMM E-value cutoff above which the HMM scores (ret_F * N) CM 
 *          hits with CM E-values better than cm_ecutoff 
 * 
 */
float FindCP9FilterThreshold(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
			     float Fmin, float Smin, float Starget, float Spad, int N, 
			     int use_cm_cutoff, float cm_ecutoff, int db_size, 
			     int emit_mode, int fthr_mode, int hmm_gum_mode, 
			     int do_fastfil, int do_Fstep, int my_rank, int nproc, int do_mpi, 
			     char *histfile, FILE *Rpts_fp, float *ret_F)
{

  /* Contract checks */
  if (!(cm->flags & CMH_CP9) || cm->cp9 == NULL) 
    cm_Fail("ERROR in FindCP9FilterThreshold() CP9 does not exist\n");
  if (Fmin < 0. || Fmin > 1.)  
    cm_Fail("ERROR in FindCP9FilterThreshold() Fmin is %f, should be [0.0..1.0]\n", Fmin);
  if (Smin < 0. || Smin > 1.)  
    cm_Fail("ERROR in FindCP9FilterThreshold() Smin is %f, should be [0.0..1.0]\n", Smin);
  if (Starget < 0. || Starget > 1.)  
    cm_Fail("ERROR in FindCP9FilterThreshold() Starget is %f, should be [0.0..1.0]\n", Starget);
  if((fthr_mode != CM_LI) && (fthr_mode != CM_GI) && (fthr_mode != CM_LC) && (fthr_mode != CM_GC))
    cm_Fail("ERROR in FindCP9FilterThreshold() fthr_mode not CM_LI, CM_GI, CM_LC, or CM_GC\n");
  if(hmm_gum_mode != CP9_L && hmm_gum_mode != CP9_G)
    cm_Fail("ERROR in FindCP9FilterThreshold() hmm_gum_mode not CP9_L or CP9_G\n");
  if(do_fastfil && (fthr_mode == CM_LI || fthr_mode == CM_GI))
    cm_Fail("ERROR in FindCP9FilterThreshold() do_fastfil TRUE, but fthr_mode CM_GI or CM_LI\n");
  if(my_rank > 0 && !do_mpi)
    cm_Fail("ERROR in FindCP9FilterThreshold() my_rank is not 0, but do_mpi is FALSE\n");
  if(emit_mode != CM_GC && emit_mode != CM_LC)
    cm_Fail("ERROR in FindCP9FilterThreshold() emit_mode not CM_LC or CM_GC\n");
  if(emit_mode == CM_LC && (fthr_mode == CM_GC || fthr_mode == CM_GI))
    cm_Fail("ERROR in FindCP9FilterThreshold() emit_mode CM_LC but score mode CM_GC or CM_GI.\n");
  if(Spad < 0 || Spad > 1.0)
    cm_Fail("ERROR in FindCP9FilterThreshold() Spad %f not between 0.0 and 1.0\n");

#if defined(USE_MPI)  && defined(NOTDEFINED)
  /* If a worker in MPI mode, we go to worker function mpi_worker_cm_and_cp9_search */
  if(my_rank > 0) 
    {
      /* Configure the CM based on the fthr mode COULD BE DIFFERENT FROM MASTER! */
      ConfigForGumbelMode(cm, fthr_mode);
      /* Configure the HMM based on the hmm_gum_mode */
      if(hmm_gum_mode == CP9_L)
	{
	  CPlan9SWConfig(cm->cp9, cm->pbegin, cm->pbegin);
	  CPlan9ELConfig(cm);
	}
      else /* hmm_gum_mode == CP9_G (it's in the contract) */
	CPlan9GlobalConfig(cm->cp9);
      CP9Logoddsify(cm->cp9);

      //mpi_worker_cm_and_cp9_search(cm, do_fastfil, my_rank);
      mpi_worker_cm_and_cp9_search_maxsc(cm, do_fastfil, do_minmax, my_rank);

      *ret_F = 0.0; /* this return value is irrelevant */
      return 0.0;   /* this return value is irrelevant */
    }
#endif 
  int            status;         /* Easel status */
  CM_t          *cm_for_scoring; /* used to score parsetrees, nec b/c  *
				  * emitting mode may != scoring mode   */
  Parsetree_t   *tr = NULL;      /* parsetree                           */
  ESL_SQ        *sq = NULL;      /* digitized sequence                  */
  float         *tr_sc;          /* scores of all parsetrees sampled    */
  float         *cm_sc_p;        /* CM scores of samples above thresh   */
  float         *hmm_sc_p;       /* HMM scores of samples above thresh  */
  float         *hmm_eval_p;     /* HMM E-values of samples above thresh*/  
  int            i  = 0;         /* counter over samples                */
  int            ip = 0;         /* counter over samples above thresh   */
  int            imax = 500 * N; /* max num samples                     */
  int            p;              /* counter over partitions             */
  int            L;              /* length of a sample                  */
  float          F;              /* fraction of hits found by HMM >= Fmin*/
  float          E;              /* HMM CP9 E-value cutoff to return    */
  float          Emin;           /* E-value that corresponds to Smin */
  float          Etarget;        /* E-value that corresponds to Starget */
  double        *cm_mu;          /* mu for each partition's CM Gumbel   */
  double        *hmm_mu;         /* mu for each partition's HMM Gumbel  */
  float         *cm_minbitsc = NULL; /* minimum CM bit score to allow to pass for each partition */
  float         *cm_maxbitsc = NULL; /* maximum CM bit score to allow to pass for each partition */
  double         tmp_K;          /* used for recalc'ing Gumbel stats for DB size */
  int            was_hmmonly;    /* flag for if HMM only search was set */
  int            nalloc;         /* for cm_sc, hmm_sc, hmm_eval, num alloc'ed */
  int            chunksize;      /* allocation chunk size               */
  float         *scores;         /* CM and HMM score returned from worker*/
  void          *tmp;            /* temp void pointer for ESL_RALLOC() */
  int            clen;           /* consensus length of CM             */
  float          avg_hit_len;    /* crude estimate of average hit len  */
  int            Fidx;           /* index within hmm_eval              */
  float          S;              /* predicted survival fraction        */
  int            init_flag;      /* used for finding F                 */ 
  int            passed_flag = FALSE;
  float          cm_sc;
  float          orig_tau;
  float          hb_sc= 0.; 
  float          S_sc = 0.;
  float          Starget_sc = 0.;

#if defined(USE_MPI)  && defined(NOTDEFINED)
  int            have_work;      /* MPI: do we still have work to send out?*/
  int            nproc_working;  /* MPI: number procs currently doing work */
  int            wi;             /* MPI: worker index                      */
  ESL_SQ       **sqlist = NULL; /* MPI: queue of digitized seqs being worked on, 1..nproc-1 */
  int           *plist = NULL;   /* MPI: queue of partition indices of seqs being worked on, 1..nproc-1 */
  Parsetree_t  **trlist = NULL;  /* MPI: queue of traces of seqs being worked on, 1..nproc-1 */
  MPI_Status     mstatus;        /* useful info from MPI Gods         */
#endif


  /* TEMPORARY! */
  do_mpi = FALSE;
  printf("in FindCP9FilterThreshold fthr_mode: %d hmm_gum_mode: %d emit_mode: %d\n", fthr_mode, 
	 hmm_gum_mode, emit_mode);

  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  cm->search_opts &= ~CM_SEARCH_HMMONLY;

  chunksize = 5 * N;
  nalloc    = chunksize;
  ESL_ALLOC(tr_sc,       sizeof(float) * nalloc);
  ESL_ALLOC(hmm_eval_p,  sizeof(float) * N);
  ESL_ALLOC(hmm_sc_p,    sizeof(float) * N);
  ESL_ALLOC(cm_sc_p,     sizeof(float) * N);
  ESL_ALLOC(cm_minbitsc, sizeof(float) * cmstats->np);
  ESL_ALLOC(cm_mu,       sizeof(double)* cmstats->np);
  ESL_ALLOC(hmm_mu,      sizeof(double)* cmstats->np);
  ESL_ALLOC(scores,      sizeof(float) * 2);
  
#if defined(USE_MPI) && defined(NOTDEFINED)
  if(do_mpi)
    {
      ESL_ALLOC(sqlist,      sizeof(ESL_SQ *)      * nproc);
      ESL_ALLOC(plist,       sizeof(int)           * nproc);
      ESL_ALLOC(trlist,      sizeof(Parsetree_t *) * nproc);
    }
#endif

  if(use_cm_cutoff) printf("CM E cutoff: %f\n", cm_ecutoff);
  else              printf("Not using CM cutoff\n");

  /* Configure the CM based on the emit mode COULD DIFFERENT FROM WORKERS! */
  ConfigForGumbelMode(cm, emit_mode);
  /* Copy the CM into cm_for_scoring, and reconfigure it if nec.,
   * We do this, so we change emission modes of the original CM */
  cm_for_scoring = DuplicateCM(cm); 
  /*if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))*/
  ConfigForGumbelMode(cm_for_scoring, fthr_mode);

  cm_CreateScanMatrixForCM(cm_for_scoring, TRUE, TRUE);

  /* Configure the HMM based on the hmm_gum_mode */
  if(hmm_gum_mode == CP9_L)
    {
      CPlan9SWConfig(cm_for_scoring->cp9, cm_for_scoring->pbegin, cm_for_scoring->pbegin);
      if(! (cm_for_scoring->flags & CMH_LOCAL_END))
	ConfigLocal(cm_for_scoring, cm_for_scoring->pbegin, cm_for_scoring->pend); 	/* need CM in local mode to calculate HMM EL probs, sloppy */
      CPlan9ELConfig(cm_for_scoring);
      if(! (cm_for_scoring->flags & CMH_LOCAL_END))
	ConfigGlobal(cm_for_scoring); 	/* return CM back to global mode, sloppy */
    }
  else /* hmm_gum_mode == CP9_G (it's in the contract) */
    CPlan9GlobalConfig(cm_for_scoring->cp9);
  CP9Logoddsify(cm_for_scoring->cp9);
  if(cm_for_scoring->config_opts & CM_CONFIG_ZEROINSERTS)
    CP9HackInsertScores(cm_for_scoring->cp9);

  /* Determine bit cutoff for each partition, calc'ed from cm_ecutoff */
  for (p = 0; p < cmstats->np; p++)
    {
      /* First determine mu based on db_size */
      tmp_K      = exp(cmstats->gumAA[hmm_gum_mode][p]->mu * cmstats->gumAA[hmm_gum_mode][p]->lambda) / 
	cmstats->gumAA[hmm_gum_mode][p]->L;
      hmm_mu[p]  = log(tmp_K * ((double) db_size)) / cmstats->gumAA[hmm_gum_mode][p]->lambda;
      tmp_K      = exp(cmstats->gumAA[fthr_mode][p]->mu * cmstats->gumAA[fthr_mode][p]->lambda) / 
	cmstats->gumAA[fthr_mode][p]->L;
      cm_mu[p]   = log(tmp_K  * ((double) db_size)) / cmstats->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cm_minbitsc[p] = cm_mu[p] - (log(cm_ecutoff) / cmstats->gumAA[fthr_mode][p]->lambda);
      if(use_cm_cutoff)
	printf("E: %f p: %d %d--%d bit score: %f\n", cm_ecutoff, p, 
	       cmstats->ps[p], cmstats->pe[p], cm_minbitsc[p]);
    }
  
  /* Do the work, emit parsetrees and collect the scores 
   * either in serial or MPI depending on do_mpi flag.
   */
  /*********************SERIAL BLOCK*************************************/
  int nleft = 0; /* number of seqs with scores < min CM score */
  int tr_np, tr_na, s1_np, s1_na, s2_np, s2_na, s3_np, s3_na;
  int do_slow = FALSE;
  char *name;
  if(Rpts_fp != NULL) do_slow = TRUE; /* we'll always find optimal CM parse to use as point for R 2D plot */

  tr_np = tr_na = s1_np = s1_na = s2_np = s2_na = s3_np = s3_na = 0;
  printf("06.11.07 Min np: %5d Smin: %12f Starget: %f Spad: %.3f do_Fstep: %d\n", N, Smin, Starget, Spad, do_Fstep);
  if(!(do_mpi))
    {
      
      orig_tau = cm_for_scoring->tau;

      /* Serial strategy: 
       * Emit sequences one at a time and search them with the CM and HMM,
       * keeping track of scores. If do_fastfil we don't have to search 
       * with the CM we just keep track of the parsetree score. 
       *
       * If do_fastfil && emit_mode is different than scoring mode we 
       * score with 'cm_for_scoring', instead of with 'cm'.
       */

      while(ip < N) /* while number seqs passed CM score threshold (ip) < N */
	{
	  ESL_ALLOC(name, sizeof(char) * 50);
	  sprintf(name, "seq%d", ip+1);
	  EmitParsetree(cm, r, name, TRUE, &tr, &sq, &L); /* TRUE: digitize the seq */
	  while(L == 0) { FreeParsetree(tr); esl_sq_Destroy(sq); EmitParsetree(cm, r, name, TRUE, &tr, &sq, &L); }
	  tr_na++;
	  passed_flag = FALSE;

	  p = cmstats->gc2p[(get_gc_comp(sq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */
	  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
	    tr_sc[i] = ParsetreeScore_Global2Local(cm_for_scoring, tr, sq->dsq, FALSE);
	  else
	    tr_sc[i] = ParsetreeScore(cm_for_scoring, tr, sq->dsq, FALSE); 
	  /*
	  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, seq);
	  ParsetreeDump(stdout, tr, cm, sq);
	  printf("%d Parsetree Score: %f\n\n", i, tr_sc[i]);
	  */

	  /* If do_minmax, check if the parsetree score less than maximum allowed */
	  if(tr_sc[i] > cm_minbitsc[p] && !do_slow) /* we know we've passed */
	    {
	      tr_np++;
	      passed_flag = TRUE;
	      cm_sc_p[ip] = tr_sc[i];
	      //printf("TR P (P: %5d L: %5d)\n", ip, nleft);
	    }
	  else
	    {
	      /* we're not sure if our optimal score exceeds cm_minbitsc */
	      /* STAGE 1 */
	      
	      /* For speed first see if a strict (high tau) HMM banded search finds a 
	       * conditional optimal parse with score > min score */
	      s1_na++;
	      cm_for_scoring->search_opts |= CM_SEARCH_HBANDED;
	      //cm_for_scoring->tau = 0.1;
	      cm_for_scoring->tau = 0.01;
	      //cm_for_scoring->tau = 0.001;
	      
	      hb_sc = OldActuallySearchTarget(cm_for_scoring, sq->dsq, 1, L,
					      0.,    /* cutoff is 0 bits (actually we'll find highest
						     * negative score if it's < 0.0) */
					      0.,    /* CP9 cutoff is 0 bits */
					      NULL,  /* don't keep results */
					      FALSE, /* don't filter with a CP9 HMM */
					      FALSE, /* we're not calcing CM  stats */
					      FALSE, /* we're not calcing CP9 stats */
					      NULL,  /* filter fraction N/A */
					      FALSE);/* do NOT align the hits */
	      //if(!do_fastfil) printf("%4d %5d %d T: %10.4f BC: %10.4f ", ip, i, passed_flag, tr_sc[i], hb_sc);
	      if(hb_sc > cm_minbitsc[p] && !do_slow)
		{
		  s1_np++;
		  passed_flag = TRUE;
		  cm_sc_p[ip] = hb_sc;
		  //printf("S1 P (P: %5d L: %5d)\n", ip, nleft);
		}
	      else /* Stage 2, search with another, less strict (lower tau)  HMM banded parse */
		{
		  s2_na++;
		  cm_for_scoring->search_opts |= CM_SEARCH_HMMSCANBANDS;
		  cm_for_scoring->tau = 1e-15;
		  cm_sc = OldActuallySearchTarget(cm_for_scoring, sq->dsq, 1, L,
						 0.,    /* cutoff is 0 bits (actually we'll find highest
							 * negative score if it's < 0.0) */
						 0.,    /* CP9 cutoff is 0 bits */
						 NULL,  /* don't keep results */
						 FALSE, /* don't filter with a CP9 HMM */
						 FALSE, /* we're not calcing CM  stats */
						 FALSE, /* we're not calcing CP9 stats */
						 NULL,  /* filter fraction N/A */
						 FALSE);/* do NOT align the hits */
		  if(cm_sc > cm_minbitsc[p] && !do_slow)
		    {
		      s2_np++;
		      passed_flag = TRUE;
		      cm_sc_p[ip] = cm_sc;
		      //printf("S2 P (P: %5d L: %5d)\n", ip, nleft);
		    }
		  else /* search for the optimal parse */
		    {
		      s3_na++;
		      /* Stage 3 do QDB CYK */
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HBANDED;
		      cm_for_scoring->search_opts &= ~CM_SEARCH_HMMSCANBANDS;
		      cm_sc = OldActuallySearchTarget(cm_for_scoring, sq->dsq, 1, L,
						     0.,    /* cutoff is 0 bits (actually we'll find highest
							     * negative score if it's < 0.0) */
						     0.,    /* CP9 cutoff is 0 bits */
						     NULL,  /* don't keep results */
						     FALSE, /* don't filter with a CP9 HMM */
						     FALSE, /* we're not calcing CM  stats */
						     FALSE, /* we're not calcing CP9 stats */
						     NULL,  /* filter fraction N/A */
						     FALSE);/* do NOT align the hits */
		      if(cm_sc > cm_minbitsc[p])
			{
			  s3_np++;
			  passed_flag = TRUE;
			  cm_sc_p[ip] = cm_sc;
			  //printf("S3 P (P: %5d L: %5d)\n", ip, nleft);
			}
		    }
		}
	    }
	  if(!passed_flag) 
	    {
	      nleft++;
	      //printf("LEFT (P: %5d L: %5d)\n", ip, nleft);
	    }
	  else if(passed_flag)
	    {
	      /* Scan seq with HMM */
	      /* DO NOT CALL OldActuallySearchTarget b/c that will run Forward then 
	       * Backward to get score of best hit, but we'll be detecting by a
	       * Forward scan (then running Backward only on hits above our threshold).
	       */
	      hmm_sc_p[ip] = CP9Forward(cm_for_scoring, sq->dsq, 1, L, cm_for_scoring->W, 0., 
					NULL,   /* don't return scores of hits */
					NULL,   /* don't return posns of hits */
					NULL,   /* don't keep track of hits */
					TRUE,   /* we're scanning */
					FALSE,  /* we're not ultimately aligning */
					FALSE,  /* we're not rescanning */
					TRUE,   /* be memory efficient */
					NULL);  /* don't want the DP matrix back */
	      hmm_eval_p[ip] = RJK_ExtremeValueE(hmm_sc_p[ip], hmm_mu[p], cmstats->gumAA[hmm_gum_mode][p]->lambda);
	      if(Rpts_fp != NULL)
		fprintf(Rpts_fp, "%.15f %.15f\n", hmm_sc_p[ip], cm_sc);
	      ip++; /* increase counter of seqs passing threshold */
	    }
	  esl_sq_Destroy(sq);
	  
	  /* Check if we need to reallocate */
	  i++;
	  if(i > imax) cm_Fail("ERROR number of attempts exceeded 500 times number of seqs.\n");
	  if (i == nalloc) 
	    {
	      nalloc += chunksize;
	      ESL_RALLOC(tr_sc,    tmp, nalloc * sizeof(float));
	    }
	  if(ip % N == 0)
	    {
	      ESL_RALLOC(hmm_eval_p, tmp, (ip + N) * sizeof(float));
	      ESL_RALLOC(hmm_sc_p,   tmp, (ip + N) * sizeof(float));
	      ESL_RALLOC(cm_sc_p,    tmp, (ip + N) * sizeof(float));
	    }
	}
    }
  N = ip; /* update N based on number of seqs sampled */
  printf("06.11.07 np: %5d nleft: %5d\n", ip, nleft);
  printf("06.11.07 tr_a: %5d tr_p: %5d\n06.11.07 s1_a: %5d s1_p: %5d\n06.11.07 s2_a: %5d s2_p: %5d\n06.11.07 s3_a: %5d s3_p: %5d\n", tr_na, tr_np, s1_na, s1_np, s2_na, s2_np, s3_na, s3_np);
  
  /*************************END OF SERIAL BLOCK****************************/
#if defined(USE_MPI)  && defined(NOTDEFINED)
  /*************************MPI BLOCK****************************/
  if(do_mpi)
    {
      /* MPI Strategy: 
       * Emit seqs one at a time from the CM, if the CM parse tree score is greater
       * than our max E-value (if do_minmax, otherwise there is no maximum),
       * send it to a worker. The worker then tries to quickly determine if
       * the sequence is within the acceptable E-value range by doing HMM
       * banded searches. If the sequence is within the E-value range,
       * it is searched with an HMM. Both the optimal CM parse scores and HMM
       * scores are passed back to the master, whose keeping track of how
       * many seqs have been sampled within the E-value range.
       *
       * If do_fastfil, the worker skips the CM search and returns 0. bits
       * as CM score, which master replaces with parsetree score. (This is 
       * an old strategy I'm probably about to deprecate (06.07.07))
       *
       * Sean's design pattern for data parallelization in a master/worker model:
       * three phases: 
       *  1. load workers;
       *  2. recv result/send work loop;
       *  3. collect remaining results
       * but implemented in a single while loop to avoid redundancy.
       */
      have_work     = TRUE;
      nproc_working = 0;
      wi            = 1;
      i             = 0;
      while (have_work || nproc_working)
	{
	  /* Get next work unit. */
	  if(ip < N) /* if number seqs passed CM score threshold (ip) < N */
	    {
	      tr_passed_flag = FALSE;
	      while(!tr_passed_flag)
		{
		  sprintf(name, "seq%d", ip+1);
		  EmitParsetree(cm, r, name, TRUE, &tr, &sq, &L); /* TRUE: digitize the seq */
		  while(L == 0) { FreeParsetree(tr); esl_sq_Destroy(sq); EmitParsetree(cm, r, name, TRUE, &tr, &sq, &L); }

		  p = cmstats->gc2p[(get_gc_comp(sq, 1, L))]; /* in get_gc_comp() should be i and j of best hit */
		  free(seq);
		  /*ParsetreeDump(stdout, tr, cm, sq);
		    printf("%d Parsetree Score: %f\n\n", (nattempts), ParsetreeScore(cm, tr, dsq, FALSE)); */
		  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
		    tr_sc[i] = ParsetreeScore_Global2Local(cm_for_scoring, trlist[wi], dsqlist[wi], FALSE);
		  else
		    tr_sc[i] = ParsetreeScore(cm_for_scoring, tr, dsq, FALSE); 
		  FreeParsetree(tr);
		  /* If do_minmax, check if the parsetree score less than maximum allowed */
		  if((!do_minmax || (do_minmax && (tr_sc[i] <= cm_maxbitsc[p]))))
		    tr_passed_flag = TRUE;
		}		    
	      if(!do_fastfil)
		FreeParsetree(tr);
	    }
	  else have_work = FALSE;
	  /* If we have work but no free workers, or we have no work but workers
	   * Are still working, then wait for a result to return from any worker.
	   */
	  if ( (have_work && nproc_working == nproc-1) || (! have_work && nproc_working > 0))
	    {
	      MPI_Recv(scores,  2, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mstatus);
	      cm_sc    = scores[0];
	      hmm_sc   = scores[1];
	      wi = mstatus.MPI_SOURCE;

	      if(do_fastfil)
		{
		  if(emit_mode == CM_GC && (fthr_mode == CM_LC || fthr_mode == CM_LI))
		    cm_sc = ParsetreeScore_Global2Local(cm_for_scoring, trlist[wi], dsqlist[wi], FALSE);
		  else
		    cm_sc = ParsetreeScore(cm_for_scoring, trlist[wi], dsqlist[wi], FALSE); 
		  FreeParsetree(trlist[wi]);
		  trlist[wi] = NULL;
		}
	      if(ip < N && 
		 ( do_minmax && (cm_sc >= cm_minbitsc[plist[wi]] && cm_sc <= cm_maxbitsc[plist[wi]])) ||
		 (!do_minmax && (cm_sc >= cm_minbitsc[plist[wi]])))
		{
		  cm_sc_p[ip]    = cm_sc;
		  hmm_sc_p[ip]   = hmm_sc;
		  hmm_eval_p[ip] = RJK_ExtremeValueE(hmm_sc, hmm_mu[plist[wi]], cmstats->gumAA[hmm_gum_mode][plist[wi]]->lambda);
		  ip++; /* increase counter of seqs passing threshold */
		}
	      nproc_working--;
	      free(dsqlist[wi]);
	      dsqlist[wi] = NULL;
	      
	      /* Check if we need to reallocate */
	      i++;
	      if(i > imax) cm_Fail("ERROR number of attempts exceeded 500 times number of seqs.\n");
	      if (i == nalloc) 
		{
		  nalloc += chunksize;
		  ESL_RALLOC(tr_sc,    tmp, nalloc * sizeof(float));
		}
	    }
	  /* If we have work, assign it to a free worker;
	   * else, terminate the free worker.
	   */
	  if (have_work) 
	    {
	      dsq_maxsc_MPISend(dsq, L, cm_maxbitsc[p], wi);
	      dsqlist[wi] = dsq;
	      plist[wi]   = p;
	      if(do_fastfil)
		trlist[wi] = tr;
	      wi++;
	      nproc_working++;
	    }
	  else 
	    dsq_maxsc_MPISend(NULL, -1, -1, wi);	
	}
    }
  /*********************END OF MPI BLOCK*****************************/
#endif /* if defined(USE_MPI) */
  /* Sort the HMM E-values with quicksort */
  esl_vec_FSortIncreasing(hmm_eval_p, N);
  esl_vec_FSortDecreasing(hmm_sc_p, N); /* TEMPORARY FOR PRINTING */

  /* Determine E-value to return based on predicted filter survival fraction */
  clen = (2*CMCountStatetype(cm, MATP_MP) + 
	  CMCountStatetype(cm, MATL_ML) + 
	  CMCountStatetype(cm, MATR_MR));
  avg_hit_len = clen; /* currently, we don't correct for local hits probably being shorter */

  /* Strategy 1, enforce minimum F, but look for highest that gives survival fraction of at least
   * Starget */
  init_flag  = TRUE;
  S          = IMPOSSIBLE;
  Emin    = ((double) db_size * Smin)    / ((2. * cm_for_scoring->W) - avg_hit_len);
  Etarget = ((double) db_size * Starget) / ((2. * cm_for_scoring->W) - avg_hit_len);
  printf("Smin: %12f\nStarget: %f\n", Smin, Starget);
  printf("Emin: %12f\nEtarget: %f\n", Emin, Etarget);

  Fidx  = (int) (Fmin * (float) N) - 1; /* off by one, ex: 95% cutoff in 100 len array is index 94 */
  E     = hmm_eval_p[Fidx];
  F     = Fmin;
  /* Are we case 1 or case 2? 
   * Case 1: our E is greater than our Etarget, so we cannot satisfy Starget AND capture F fraction of
   *         the true CM hits, it's more impt we capture F fraction, so return S > Starget.
   * Case 2: our E is less than our Etarget so we can satisfy both criteria, we have to choose an S
   *         now. One option is If (do_Fstep) increase F toward 1.0 while E < Etarget. Then add
   *         Spad fraction of (bitscore(S) - bitscore(max(Starget, Smin))) to bitscore(S), and calculate
   *         the new S. If Spad == 0., speed trumps filter safety, if S == 1.0, we want the filter
   *         as safe as possible, we're returning S = Starget. 
   */
  if(E > Etarget) 
    {
      F = ((float) (Fidx + 1.)) / (float) N;
      S = (E * ((2. * cm_for_scoring->W) - avg_hit_len)) / ((double) db_size); 
      if(E < Emin)  /* no need for E < Emin */
	  printf("Case 1A rare case: init Emin %12f > E %f > Etarget %f F: %f S: %.12f Spad: %.3f\n", Emin, E, Etarget, F, S, Spad);
      else
	  printf("Case 1B bad  case: init E %f > Etarget %f && E > Emin %12f F: %f S: %.12f Spad: %.3f\n", E, Etarget, Emin, F, S, Spad);
      if(E < Emin) E = Emin; /* No reason to have an E below Emin */
    }
  else
    {
      /* If(do_Fstep): increase F until we're about to go over Etarget (which gives our target 
	 survival fraction Starget). */
      if(do_Fstep)
	while(((Fidx+1) < N) && hmm_eval_p[(Fidx+1)] < Etarget)
	  Fidx++;
      F = ((float) (Fidx + 1.)) / (float) N;

      /* Subtract Spad * (score(S) - score(max(Starget, Smin))) from score(S) and recalculate E&S */
      /* Use partition that covers GC = 50 */
      p = cm->stats->gc2p[50];
      E = hmm_eval_p[Fidx];
      printf("Before Spad E: %f\n", E);
      S_sc =    (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
		 (log(E) / cm->stats->gumAA[hmm_gum_mode][p]->lambda)); 
      Starget_sc = (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
		    (log(Etarget) / cm->stats->gumAA[hmm_gum_mode][p]->lambda));

      printf("S_sc 0: %.3f - %.3f * %.3f = ", S_sc, Spad, S_sc - Starget_sc);
      S_sc -= Spad * (S_sc - Starget_sc); /* Spad may be 0. */

      printf(" S1: %.3f\n", S_sc);
      /* now recalculate what E and S should be based on S_sc */
      E = RJK_ExtremeValueE(S_sc, cm->stats->gumAA[hmm_gum_mode][p]->mu, 
			    cmstats->gumAA[hmm_gum_mode][p]->lambda);
      S = (E * ((2. * cm_for_scoring->W) - avg_hit_len)) / ((double) db_size); 

      if(fabs(E - Emin) < 0.01 || E < Emin)  /* no need for E < Emin */
	printf("Case 2A: best case Emin %12f > E %f < Etarget %f F: %f S: %.12f Spad: %.3f\n", Emin, E, Etarget, F, S, Spad);
       else
	printf("Case 2B: good case Emin %12f < E %f < Etarget %f F: %f S: %.12f\n Spad: %.3f", Emin, E, Etarget, F, S, Spad);
      if(E < Emin) E = Emin; /* No reason to have an E below Emin */
    }
  /* Print cutoff info to Rpts file for 2D plot if nec */
  if(Rpts_fp != NULL)
    {
      fprintf(Rpts_fp, "TSC %.15f\n", (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
				       (log(Etarget) / cm->stats->gumAA[hmm_gum_mode][p]->lambda)));
      fprintf(Rpts_fp, "MSC %.15f\n", (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
				       (log(Emin) / cm->stats->gumAA[hmm_gum_mode][p]->lambda)));
      fprintf(Rpts_fp, "FSC %.15f\n", (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
					(log(E) / cm->stats->gumAA[hmm_gum_mode][p]->lambda)));
      fprintf(Rpts_fp, "FMINSC %.15f\n", (cm->stats->gumAA[hmm_gum_mode][p]->mu - 
					(log(hmm_eval_p[(int) (Fmin * (float) N)-1]) / cm->stats->gumAA[hmm_gum_mode][p]->lambda)));
      fprintf(Rpts_fp, "F   %.15f\n", F);
      fprintf(Rpts_fp, "CMGUM %.15f %15f\n", cm->stats->gumAA[emit_mode][p]->lambda, cm->stats->gumAA[emit_mode][p]->mu);
      fprintf(Rpts_fp, "HMMGUM %.15f %15f\n", cm->stats->gumAA[hmm_gum_mode][p]->lambda, cm->stats->gumAA[hmm_gum_mode][p]->mu);
      fprintf(Rpts_fp, "STARGET %.15f\n", Starget);
      fprintf(Rpts_fp, "SMIN %.15f\n", Smin);
      fprintf(Rpts_fp, "SPAD %.15f\n", Spad);
      fprintf(Rpts_fp, "FMIN %.15f\n", Fmin);
      fprintf(Rpts_fp, "FSTEP %d\n", do_Fstep);
      fprintf(Rpts_fp, "ECUTOFF %.15f\n", cm_ecutoff);
      for (p = 0; p < cmstats->np; p++)
	fprintf(Rpts_fp, "SCCUTOFF %d %.15f\n", p, cm_minbitsc[p]);
      fclose(Rpts_fp);
    }	  

  /* Make sure our E is less than the DB size and greater than Emin */
  if(E > ((float) db_size)) /* E-val > db_size is useless */
    {
      printf("Case 3 : worst case E (%f) > db_size (%f)\n", E, (double) db_size);
      E = (float) db_size;
    }  
  
  /* Informative, temporary print statements */
  for (i = ((int) (Fmin  * (float) N) -1); i < N; i++)
    printf("%d i: %4d hmm sc: %10.4f hmm E: %10.4f\n", hmm_gum_mode, i, hmm_sc_p[i], hmm_eval_p[i]);
  printf("\nSummary: %d %d %d %d %f %f\n", fthr_mode, hmm_gum_mode, do_fastfil, emit_mode,
	 (hmm_sc_p[Fidx]), (hmm_eval_p[Fidx]));
  printf("05.21.07 %d %d %f %f\n", fthr_mode, hmm_gum_mode, hmm_sc_p[Fidx], E);

  /* Reset CM_SEARCH_HMMONLY search option as it was when function was entered */
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;

  /* Clean up and exit */
  free(tr_sc);
  free(hmm_eval_p);
  free(hmm_sc_p);
  free(cm_sc_p);
  free(hmm_mu);
  free(cm_mu);
  free(cm_minbitsc);
  free(cm_maxbitsc);
  free(scores);
  cm_for_scoring->tau = orig_tau;
  FreeCM(cm_for_scoring);
  /* Return threshold */
  *ret_F = F;
  printf("F: %f\n", *ret_F);
  return E;

 ERROR:
  cm_Fail("Reached ERROR in FindCP9FilterThreshold()\n");
  return 0.;
}



/* Following functions for CPlan9 HMMs were deprecated 01.04.07,
 * we never use these aspects of a CP9 HMM.
 */
#if 0
/* Function: CPlan9SetName()
 * 
 * Purpose:  Change the name of a CPlan9 HMM. Convenience function.
 *      
 * Note:     Trailing whitespace and \n's are chopped.     
 */
void
CPlan9SetName(CP9_t *hmm, char *name)
{
  if (hmm->name != NULL) free(hmm->name);
  hmm->name = Strdup(name);
  StringChop(hmm->name);
}
/* Function: Cplan9SetAccession()
 * 
 * Purpose:  Change the accession number of a Cplan9 HMM. Convenience function.
 *      
 * Note:     Trailing whitespace and \n's are chopped.     
 */
void
CPlan9SetAccession(CP9_t *hmm, char *acc)
{
  if (hmm->acc != NULL) free(hmm->acc);
  hmm->acc = Strdup(acc);
  StringChop(hmm->acc);
  hmm->flags |= CPLAN9_ACC;
}

/* Function: CPlan9SetDescription()
 * 
 * Purpose:  Change the description line of a Cplan9 HMM. Convenience function.
 * 
 * Note:     Trailing whitespace and \n's are chopped.
 */
void
CPlan9SetDescription(CP9_t *hmm, char *desc)
{
  if (hmm->desc != NULL) free(hmm->desc);
  hmm->desc = Strdup(desc);
  StringChop(hmm->desc); 
  hmm->flags |= CPLAN9_DESC;
}

/* Function: CPlan9ComlogAppend()
 * Date:     SRE, Wed Oct 29 09:57:30 1997 [TWA 721 over Greenland] 
 * 
 * Purpose:  Concatenate command line options and append to the
 *           command line log.
 */
void
CPlan9ComlogAppend(CP9_t *hmm, int argc, char **argv)
{
  int len;
  int i;

  /* figure out length of command line, w/ spaces and \n */
  len = argc;
  for (i = 0; i < argc; i++)
    len += strlen(argv[i]);

  /* allocate */
  if (hmm->comlog != NULL)
    {
      len += strlen(hmm->comlog);
      ESL_RALLOC(hmm->comlog, tmp, sizeof(char)* (len+1));
    }
  else
    {
      ESL_ALLOC(hmm->comlog, sizeof(char)* (len+1));
      *(hmm->comlog) = '\0'; /* need this to make strcat work */
    }

  /* append */
  strcat(hmm->comlog, "\n");
  for (i = 0; i < argc; i++)
    {
      strcat(hmm->comlog, argv[i]);
      if (i < argc-1) strcat(hmm->comlog, " ");
    }
}

/* Function: CPlan9SetCtime()
 * Date:     SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 * 
 * Purpose:  Set the ctime field in a new HMM to the current time.
 */
void
CPlan9SetCtime(CP9_t *hmm)
{
  time_t date = time(NULL);
  if (hmm->ctime != NULL) free(hmm->ctime);
  hmm->ctime = Strdup(ctime(&date));
  StringChop(hmm->ctime);
}
#endif

/* OLD MPI functions */
#if 0
/**************************************************************************************/
/* EPN, Thu May 10 10:11:18 2007 New functions roughly following Easel/H3 conventions */
/* Function: mpi_worker_search_target()
 * Incept:   EPN, Wed May  9 17:07:48 2007
 * Purpose:  The main control for an MPI worker process for searching sequences. 
 *           Worker receives CM, then loops over receipt of sequences, returning
 *           best score and results data structure for each.
 *           Never do revcomp, we'll call this function twice once with 
 *           plus once with minus strand.
 */
void
mpi_worker_search_target(CM_t *cm, int my_rank)
{
  ESL_DSQ *dsq = NULL;
  int   L;
  float best_sc;

  int doing_cm_stats  = FALSE;
  int doing_cp9_stats = FALSE;

  if(cm->search_opts & CM_SEARCH_HMMONLY) doing_cp9_stats = TRUE;
  else doing_cm_stats = TRUE;
  /* Main loop */
  while (dsq_MPIRecv(&dsq, &L) == eslOK)
    {
      best_sc = actually_search_target(cm, dsq, 1, L, 
				       0.,    /* minimum CM bit cutoff, irrelevant (?) */
				       0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
				       NULL,  /* do not keep results */
				       FALSE, /* do not filter with a CP9 HMM */
				       doing_cm_stats, doing_cp9_stats,
				       NULL,  /* filter fraction, nobody cares */
				       FALSE);/* don't align hits */
      
      MPI_Send(&(best_sc), 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
      free(dsq);
    }
  return;

}

/* Function:  dsq_MPISend()
 * Incept:    EPN, Wed May  9 17:30:14 2007
 *
 * Purpose:   Send sequence <dsq> to processor <dest>.
 *            
 * Returns:   eslOK on success; eslEINVAL if <dsq> is NULL
 *            and eslESYS if there is an MPI error.
 */
int
dsq_MPISend(ESL_DSQ *dsq, int L, int dest)
{
  int status;

  if(dsq == NULL) { status = eslESYS; goto ERROR; }

  if(MPI_Send(&L, 1, MPI_INT, dest, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed.");
  /* receiver will now allocate storage, before reading on...*/
  if(MPI_Send(dsq, (L+2), MPI_BYTE, dest, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed.");
  return eslOK;

 ERROR: 
  return status;
}

/* Function: mpi_worker_cm_and_cp9_search()
 * Incept:   EPN, Thu May 10 10:04:02 2007
 * Purpose:  The main control for an MPI worker process for searching sequences
 *           twice, once with a CM and once with a CP9, both scores are returned.
 *           Called in mpi_FindCP9FilterThreshold9).
 * Args:
 *           cm       - the covariance model
 *           do_fast  - don't search with CM, only do CP9 search
 *           my_rank  - my MPI rank
 */
void
mpi_worker_cm_and_cp9_search(CM_t *cm, int do_fast, int my_rank)
{
  int status;
  ESL_DSQ *dsq = NULL;
  int   L;
  float *scores = NULL;
  ESL_ALLOC(scores, sizeof(float) * 2);
  int was_hmmonly;
  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  /* Main loop */
  while (dsq_MPIRecv(&dsq, &L) == eslOK)
    {
      /* Do the CM search first */
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      if(do_fast)
	scores[0] = 0.;
      else
	scores[0] = actually_search_target(cm, dsq, 1, L, 
					   0.,    /* minimum CM bit cutoff, irrelevant (?) */
					   0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
					   NULL,  /* do not keep results */
					   FALSE, /* do not filter with a CP9 HMM */
					   FALSE, FALSE, /* not doing CM or CP9 Gumbel calcs */
					   NULL,  /* filter fraction, nobody cares */
					   FALSE);/* don't align hits */
      /* DO NOT CALL actually_search_target b/c that will run Forward then 
       * Backward to get score of best hit, but we'll be detecting by a
       * Forward scan (then running Backward only on hits above our threshold),
       * since we're calc'ing the threshold here it's impt we only do Forward.
       */
      scores[1] =  CP9Forward(cm, dsq, 1, L, cm->W, 0., 
			      NULL,   /* don't return scores of hits */
			      NULL,   /* don't return posns of hits */
			      NULL,   /* don't keep track of hits */
			      TRUE,   /* we're scanning */
			      FALSE,  /* we're not ultimately aligning */
			      FALSE,  /* we're not rescanning */
			      TRUE,   /* be memory efficient */
			      NULL);  /* don't want the DP matrix back */
      
      MPI_Send(scores, 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); /* send together so results don't interleave */
      free(dsq);
    }
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;
  if(scores != NULL) free(scores);
  /*printf("\trank: %d RETURNING!\n", my_rank);*/
  return;
  
 ERROR:
  if (dsq != NULL) free(dsq);
  if (scores != NULL) free(scores);
  return;
}

/* Function: mpi_worker_cm_and_cp9_search_maxsc()
 * Incept:   EPN, Thu Jun  7 15:02:54 2007   
 * Purpose:  The main control for an MPI worker process for searching sequences
 *           with decreasingly fast techniques, quitting if any technique 
 *           returns a score greater than some specified maximum bit score. 
 *           The goal is to determine if the optimal parse is within a 
 *           given range during a empirical HMM filter threshold calculation. 
 *           Called in mpi_FindCP9FilterThreshold().
 * Args:
 *           cm       - the covariance model
 *           do_fast  - don't search with CM, only do CP9 search
 *           my_rank  - my MPI rank
 */
void
mpi_worker_cm_and_cp9_search_maxsc(CM_t *cm, int do_fast, int do_minmax, int my_rank)
{
  int status;
  char *dsq = NULL;
  int   L;
  float maxsc;
  float *scores = NULL;
  ESL_ALLOC(scores, sizeof(float) * 2);
  int was_hmmonly;
  int was_hbanded;
  float orig_tau;
  orig_tau = cm->tau;

  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  if(cm->search_opts & CM_SEARCH_HBANDED) was_hbanded = TRUE;
  else was_hbanded = FALSE;
  /* Main loop */
  while (dsq_maxsc_MPIRecv(&dsq, &L, &maxsc) == eslOK)
    {
      /* Do the CM search first */
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      if(do_fast)
	scores[0] = 0.;
      else if(do_minmax)
	{
	  cm->search_opts |= CM_SEARCH_HBANDED;
	  cm->tau = 0.1;
	  scores[0] = actually_search_target(cm, dsq, 1, L,
					     0.,    /* cutoff is 0 bits (actually we'll find highest
						     * negative score if it's < 0.0) */
					     0.,    /* CP9 cutoff is 0 bits */
					     NULL,  /* don't keep results */
					     FALSE, /* don't filter with a CP9 HMM */
					     FALSE, /* we're not calcing CM  stats */
					     FALSE, /* we're not calcing CP9 stats */
					     NULL,  /* filter fraction N/A */
					     FALSE);/* don't align hits */
	  
	  if(scores[0] < maxsc) /* search with another, less strict (lower tau)  HMM banded parse */
	    {
	      cm->tau = 1e-10;
	      scores[0] = actually_search_target(cm, dsq, 1, L,
						 0.,    /* cutoff is 0 bits (actually we'll find highest
							 * negative score if it's < 0.0) */
						 0.,    /* CP9 cutoff is 0 bits */
						 NULL,  /* don't keep results */
						 FALSE, /* don't filter with a CP9 HMM */
						 FALSE, /* we're not calcing CM  stats */
						 FALSE, /* we're not calcing CP9 stats */
						 NULL,  /* filter fraction N/A */
						 FALSE);/* don't align hits */
	    }
	}
      else
	scores[0] = actually_search_target(cm, dsq, 1, L, 
					   0.,    /* minimum CM bit cutoff, irrelevant (?) */
					   0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
					   NULL,  /* do not keep results */
					   FALSE, /* do not filter with a CP9 HMM */
					   FALSE, FALSE, /* not doing CM or CP9 Gumbel calcs */
					   NULL,  /* filter fraction, nobody cares */
					   FALSE);/* don't align hits */
      
      /* Now do HMM search, but if do_minmax, only do HMM search 
       * if our CM score hasn't exceeded the max */
      if(do_minmax && scores[0] > maxsc)
	scores[1] = 0.;
      else
	/* DO NOT CALL actually_search_target b/c that will run Forward then 
	 * Backward to get score of best hit, but we'll be detecting by a
	 * Forward scan (then running Backward only on hits above our threshold),
	 * since we're calc'ing the threshold here it's impt we only do Forward.
	 */
	scores[1] =  CP9Forward(cm, dsq, 1, L, cm->W, 0., 
				NULL,   /* don't return scores of hits */
				NULL,   /* don't return posns of hits */
				NULL,   /* don't keep track of hits */
				TRUE,   /* we're scanning */
				FALSE,  /* we're not ultimately aligning */
				FALSE,  /* we're not rescanning */
				TRUE,   /* be memory efficient */
				NULL);  /* don't want the DP matrix back */
      MPI_Send(scores, 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); /* send together so results don't interleave */
      free(dsq);
    }
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;
  if(was_hbanded) cm->search_opts |= CM_SEARCH_HBANDED;
  else cm->search_opts &= ~CM_SEARCH_HBANDED;
  
  if(scores != NULL) free(scores);
  /*printf("\trank: %d RETURNING!\n", my_rank);*/
  return;
  
 ERROR:
  if (dsq != NULL) free(dsq);
  if (scores != NULL) free(scores);
  return;
}

/* Function:  dsq_MPIRecv()
 * Incept:    EPN, Wed May  9 17:34:43 2007
 *
 * Purpose:   Receive a sequence sent from the master MPI process (src=0)
 *            on a worker MPI process. 
 *            
 *            If it receives an end-of-data signal, returns <eslEOD>.
 */
int
dsq_MPIRecv(ESL_DSQ **ret_dsq, int *ret_L)
{
  int status;
  char *dsq = NULL;
  MPI_Status mpistatus;
  int L;
  
  if(MPI_Recv(&L, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_EXCEPTION(eslESYS, "mpi receive failed.");
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if(MPI_Recv(dsq, (L+2), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0); ESL_EXCEPTION(eslESYS, "mpi receive failed.");
  dsq[0] = dsq[(L+1)] = eslDSQ_SENTINEL;
  *ret_L   = L;
  *ret_dsq = dsq;
  
  return eslOK;

 ERROR:
  return status;
}

/* Function:  dsq_maxsc_MPIRecv()
 * Incept:    EPN, Thu Jun  7 15:00:29 2007    
 *
 * Purpose:   Receive a sequence and maximum score 
 *            sent from the master MPI process (src=0)
 *            on a worker MPI process. 
 *            
 *            If it receives an end-of-data signal, returns <eslEOD>.
 */
int
dsq_maxsc_MPIRecv(char **ret_dsq, int *ret_L, float *ret_maxsc)
{
  int status;
  char *dsq = NULL;
  MPI_Status mpistatus;
  int L;
  float maxsc;
  
  MPI_Recv(&L, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  if (L == -1) return eslEOD;
  ESL_ALLOC(dsq, sizeof(char) * (L+2));
  MPI_Recv(dsq, (L+2), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(&maxsc, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  *ret_L   = L;
  *ret_dsq = dsq;
  *ret_maxsc = maxsc;
  return eslOK;

 ERROR:
  return status;
}

#endif


#if 0
/* Here are the P7 versions of the HMM banding related
 * functions, for reference */

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

/* EPN, Mon Dec 10 09:20:02 2007, simplified use of gumbel modes in cmcalibrate,
 * making ConfigForGumbelMode unnecessary */
#if 0
/*
 * Function: ConfigForGumbelMode
 * Date:     EPN, Thu May  3 09:05:48 2007
 * Purpose:  Configure a CM and it's CP9 for determining statistics for 
 *           a specific 'gum_mode'.
 *
 *           0. CM_LC : !cm->search_opts & CM_SEARCH_INSIDE  w/  local CM
 *           1. CM_GC : !cm->search_opts & CM_SEARCH_INSIDE  w/ glocal CM
 *           2. CM_LI :  cm->search_opts & CM_SEARCH_INSIDE  w/  local CM
 *           3. CM_GI :  cm->search_opts & CM_SEARCH_INSIDE  w/ glocal CM
 *           4. CP9_LV:  cm->search_opts & CM_SEARCH_HMMVITERBI
 *                      !cm->search_opts & CM_SEARCH_HMMFORWARD w/  local CP9 HMM
 *           5. CP9_LV:  cm->search_opts & CM_SEARCH_HMMVITERBI
 *                      !cm->search_opts & CM_SEARCH_HMMFORWARD w/ glocal CP9 HMM
 *           6. CP9_LF: !cm->search_opts & CM_SEARCH_HMMVITERBI
 *                       cm->search_opts & CM_SEARCH_HMMFORWARD w/  local CP9 HMM
 *           7. CP9_LF: !cm->search_opts & CM_SEARCH_HMMVITERBI
 *                       cm->search_opts & CM_SEARCH_HMMFORWARD w/ glocal CP9 HMM
 * 
 * Args:
 *           CM           - the covariance model
 *           gum_mode     - the mode 0..7
 */
int
ConfigForGumbelMode(CM_t *cm, int gum_mode)
{
  int do_cm_local  = FALSE;
  int do_cp9_local = FALSE;

  /*printf("in ConfigForGumbelMode\n");*/
  /* First set search opts and flags based on gum_mode */
  switch (gum_mode) {
  case CM_LC: /* local CYK */
    /*printf("CM_LC\n");*/
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    do_cm_local  = TRUE;
    break;
  case CM_GC: /* glocal CYK */
    /*printf("CM_GC\n");*/
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    do_cm_local  = FALSE;
    break;
  case CM_LI: /* local inside */
    /*printf("CM_LI\n");*/
    cm->search_opts |= CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    do_cm_local  = TRUE;
    break;
  case CM_GI: /* glocal inside */
    /*printf("CM_GI\n");*/
    cm->search_opts |= CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    do_cm_local  = FALSE;
    break;
  case CP9_LV: /* local CP9 Viterbi */
    /*printf("CP9_LV\n");*/
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts |= CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    do_cm_local   = TRUE; /* need CM local ends to make CP9 local ends */
    do_cp9_local  = TRUE;
    break;
  case CP9_GV: /* glocal CP9 Viterbi */
    /*printf("CP9_GV\n");*/
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts |=  CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    do_cp9_local  = FALSE;
    break;
  case CP9_LF: /* local CP9 Forward */
    /*printf("CP9_LF\n");*/
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts |= CM_SEARCH_HMMFORWARD;
    do_cm_local   = TRUE; /* need CM local ends to make CP9 local ends */
    do_cp9_local  = TRUE;
    break;
  case CP9_GF: /* glocal CP9 Forward */
    /*printf("CP9_GF\n");*/
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts |= CM_SEARCH_HMMFORWARD;
    do_cp9_local  = FALSE;
    break;
  default: 
    cm_Fail("ERROR unrecognized gum_mode: %d in ConfigForGumbelMode");
  }
  /* configure CM and, if needed, CP9 */
  if(do_cm_local || do_cp9_local) 
    {
      /* If we're in local, wastefully convert to global, 
       * then back to local, so we follow our rule that ConfigLocal()
       * cannot be called with a model already locally configured.
       * That rule was put in place to force caller to understand what
       * it's doing. */
      if(cm->flags & CMH_LOCAL_BEGIN || cm->flags & CMH_LOCAL_END) 
	ConfigGlobal(cm);
      ConfigLocal(cm, cm->pbegin, cm->pend);
    }
  else if(cm->flags & CMH_LOCAL_BEGIN || cm->flags & CMH_LOCAL_END) /* these *should* both either be up or down */
    ConfigGlobal(cm);
  CMLogoddsify(cm);
  if(cm->config_opts & CM_CONFIG_ZEROINSERTS)
    CMHackInsertScores(cm);	    /* insert emissions are all equiprobable,
				     * makes all CP9 (if non-null) inserts equiprobable*/
  if((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD))
    {
      if(!(cm->flags & CMH_CP9) || cm->cp9 == NULL) /* error, we should have one */
	cm_Fail("CP9 must already be built in ConfigForGumbelMode()\n");
      if(do_cp9_local)
	{
	  /* To do: Make the CP9 local to match the CM, as close as we can */
	  /*CPlan9CMLocalBeginConfig(cm); <-- Finish this function */
	  CPlan9SWConfig(cm->cp9, cm->pbegin, cm->pbegin); 
	  CPlan9ELConfig(cm);

	  /* CPlan9SWConfig(cm->cp9, ((cm->cp9->M)-1.)/cm->cp9->M,*/  /* all start pts equiprobable, including 1 */
	  /*                  ((cm->cp9->M)-1.)/cm->cp9->M);*/          /* all end pts equiprobable, including M */
	}
      else
	CPlan9GlobalConfig(cm->cp9);
      CP9Logoddsify(cm->cp9);
      if(cm->config_opts & CM_CONFIG_ZEROINSERTS)
	CP9HackInsertScores(cm->cp9);
    }
  return eslOK;
}


/**************************************************************************
 * EPN 09.06.06
 * Function: check_sub_cm_by_sampling2()
 *
 * Purpose:  Given a CM and a sub CM that is supposed to mirror 
 *           the CM as closely as possible between two given consensus
 *           columns (spos and epos), check that the sub_cm was correctly 
 *           constructed. 
 *           
 *           The approach is to sample from the CM and the sub_cm 
 *           and use those samples to build two CP9 HMMs, then 
 *           compare those two CP9 HMMs.
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t  *sub_cm     - the sub CM built from the orig_cm
 * ESL_RANDOMNESS *r - source of randomness
 * int spos          - first consensus column in cm that hmm models (often 1)
 * int epos          -  last consensus column in cm that hmm models 
 * int nseq          - number of sequences to sample to build the new HMMs.
 *
 * Returns: TRUE: if CM and sub CM are "close enough" (see code)
 *          FALSE: otherwise
 */
int 
check_sub_cm_by_sampling2(CM_t *orig_cm, CM_t *sub_cm, ESL_RANDOMNESS *r, int spos, int epos, int nseq)
{
  int status;
  CP9_t       *orig_hmm;/* constructed CP9 HMM from the sub_cm */
  CP9_t       *sub_hmm; /* constructed CP9 HMM from the sub_cm */
  int ret_val;                    /* return value */
  Parsetree_t **tr;               /* Parsetrees of emitted aligned sequences */
  ESL_SQ  **sq;                   /* sequences */
  ESL_MSA  *msa;                  /* alignment */
  float    *wgt;
  char     *name;                 /* name for emitted seqs */
  int i;
  int L;
  int apos;
  int *matassign;
  int *useme;
  CP9trace_t **cp9_tr;          /* fake tracebacks for each seq            */
  int msa_nseq;                 /* this is the number of sequences per MSA,
				 * current strategy is to sample (nseq/nseq_per_msa)
				 * alignments from the CM, and add counts from
				 * each to the shmm in counts form (to limit memory)
				 */
  int nsampled;                 /* number of sequences sampled thus far */
  int debug_level;
  int cc;
  char         *tmp_name;           /* name for seqs */
  char         *tmp_text_sq;        /* text seqs */
  char errbuf[cmERRBUFSIZE];
  
  debug_level = 0;
  ret_val = TRUE;
  msa_nseq = 1000;
  msa = NULL;
  
  /* Build two CP9 HMMs */
  /* the orig_hmm only models consensus positions spos to epos of the orig_cm */
  orig_hmm = AllocCPlan9((epos-spos+1), orig_cm->abc);
  ZeroCPlan9(orig_hmm);
  CPlan9SetNullModel(orig_hmm, orig_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */
  
  sub_hmm = AllocCPlan9((epos-spos+1), orig_cm->abc);
  ZeroCPlan9(sub_hmm);
  CPlan9SetNullModel(sub_hmm, sub_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */
  
  /* First sample from the orig_cm and use the samples to fill in orig_hmm
   * sample MSA(s) from the CM 
   */
  nsampled = 0;
  ESL_ALLOC(sq, sizeof(ESL_SQ *)     * msa_nseq);
  ESL_ALLOC(tr, (sizeof(Parsetree_t) * msa_nseq));
  ESL_ALLOC(wgt,(sizeof(float)       * msa_nseq));
  esl_vec_FSet(wgt, msa_nseq, 1.0);
  int namelen = 3 + IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */

  while(nsampled < nseq)
    {
      /*printf("nsampled: %d\n", nsampled);*/
      if(nsampled != 0)
	{
	  /* clean up from previous MSA */
	  esl_msa_Destroy(msa);
	  free(matassign);
	  free(useme);
	  for (i = 0; i < msa_nseq; i++)
	    {
	      CP9FreeTrace(cp9_tr[i]);
	      FreeParsetree(tr[i]);
	      esl_sq_Reuse(sq[i]);
	    }
	  free(cp9_tr);
	}
      /* Emit msa_nseq parsetrees from the CM */
      if(nsampled + msa_nseq > nseq)
	msa_nseq = nseq - nsampled;
      for (i = 0; i < msa_nseq; i++)
	{
	  ESL_ALLOC(name, sizeof(char) * namelen);
	  sprintf(name, "seq%d", i+1);
	  if((status = EmitParsetree(orig_cm, errbuf, r, name, FALSE, &(tr[i]), &(sq[i]), &L)) != eslOK) cm_Fail(errbuf);
	  free(name);
	}
      /* Build a new MSA from these parsetrees */
      Parsetrees2Alignment(orig_cm, orig_cm->abc, sq, NULL, tr, msa_nseq, TRUE, FALSE, &msa);
      /* MSA should be in text mode, not digitized */
      if(msa->flags & eslMSA_DIGITAL)
	cm_Fail("ERROR in sub_cm_check_by_sampling(), sampled MSA should NOT be digitized.\n");
      
      /* Truncate the alignment prior to consensus column spos and after 
	 consensus column epos */
      ESL_ALLOC(useme, sizeof(int) * (msa->alen+1));
      for (apos = 0, cc = 0; apos < msa->alen; apos++)
	{
	  /* Careful here, placement of cc++ increment is impt, 
	   * we want all inserts between cc=spos-1 and cc=spos,
	   * and between cc=epos and cc=epos+1.
	   */
	  if(cc < (spos-1) || cc > epos)
	    useme[apos] = 0;
	  else
	    useme[apos] = 1;
	  if (! esl_abc_CIsGap(msa->abc, msa->rf[apos]))
	    { 
	      cc++; 
	      if(cc == (epos+1))
		useme[apos] = 0; 
	      /* we misassigned this guy, overwrite */ 
	    }
	}
      esl_msa_ColumnSubset(msa, useme);
      
      /* Shorten the dsq's */
      for (i = 0; i < msa_nseq; i++)
	{
	  MakeDealignedString(msa->abc, msa->aseq[i], msa->alen, msa->aseq[i], &(tmp_text_sq)); 
	  ESL_ALLOC(tmp_name, sizeof(char) * namelen);
	  sprintf(tmp_name, "seq%d", i+1);
	  esl_sq_CreateFrom(tmp_name, tmp_text_sq, NULL, NULL, NULL);
	  free(tmp_text_sq);
	  if(esl_sq_Digitize(msa->abc, sq[i]) != eslOK)
	    cm_Fail("ERROR digitizing sequence in CP9_check_by_sampling().\n");
	}
      
      /* Determine match assignment from RF annotation
       */
      ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));
      matassign[0] = 0;
      for (apos = 0; apos < msa->alen; apos++)
	{
	  matassign[apos+1] = 0;
	  if (!esl_abc_CIsGap(msa->abc, msa->rf[apos])) 
	    matassign[apos+1] = 1;
	}
      /* make fake tracebacks for each seq */
      CP9_fake_tracebacks(msa, matassign, &cp9_tr);
      
      /* build model from tracebacks (code from HMMER's modelmakers.c::matassign2hmm() */
      for (i = 0; i < msa->nseq; i++) {
	CP9TraceCount(orig_hmm, sq[i]->dsq, msa->wgt[i], cp9_tr[i]);
      }
      nsampled += msa_nseq;
    }
  /*Next, renormalize the orig_hmm and logoddisfy it */
  CPlan9Renormalize(orig_hmm);
  CP9Logoddsify(orig_hmm);
  
  /* clean up from previous MSA */
  esl_msa_Destroy(msa);
  free(matassign);
  free(useme);
  for (i = 0; i < msa_nseq; i++)
    {
      CP9FreeTrace(cp9_tr[i]);
      FreeParsetree(tr[i]);
      esl_sq_Destroy(sq[i]);
    }
  free(cp9_tr);
  
  /* Now for the sub_hmm. Sample from the sub_cm and use the 
   * samples to fill in sub_hmm sample MSA(s) from the CM */
  nsampled = 0;
  ESL_ALLOC(sq, sizeof(ESL_SQ *)     * msa_nseq);
  ESL_ALLOC(tr, (sizeof(Parsetree_t) * msa_nseq));
  ESL_ALLOC(wgt,(sizeof(float)       * msa_nseq));
  esl_vec_FSet(wgt, msa_nseq, 1.0);
  
  while(nsampled < nseq)
    {
      /*printf("nsampled: %d\n", nsampled);*/
      if(nsampled != 0)
	{
	  /* clean up from previous MSA */
	  esl_msa_Destroy(msa);
	  free(matassign);
	  for (i = 0; i < msa_nseq; i++)
	    {
	      CP9FreeTrace(cp9_tr[i]);
	      FreeParsetree(tr[i]);
	      esl_sq_Reuse(sq[i]);
	    }
	  free(cp9_tr);
	}
      /* Emit msa_nseq parsetrees from the CM */
      if(nsampled + msa_nseq > nseq)
	msa_nseq = nseq - nsampled;
      for (i = 0; i < msa_nseq; i++)
	{
	  ESL_ALLOC(name, sizeof(char) * namelen);
	  sprintf(name, "seq%d", i+1);
	  if((status = EmitParsetree(sub_cm, errbuf, r, name, FALSE, &(tr[i]), &(sq[i]), &L)) != eslOK) cm_Fail(errbuf);
	  free(name);
	}
      /* Build a new MSA from these parsetrees */
      Parsetrees2Alignment(sub_cm, sub_cm->abc, sq, NULL, tr, msa_nseq, TRUE, FALSE, &msa);
      /* MSA should be in text mode, not digitized */
      if(msa->flags & eslMSA_DIGITAL)
	cm_Fail("ERROR in sub_cm_check_by_sampling(), sampled MSA should NOT be digitized.\n");
      
      /* Determine match assignment from RF annotation
       */
      ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));
      matassign[0] = 0;
      for (apos = 0; apos < msa->alen; apos++)
	{
	  matassign[apos+1] = 0;
	  if (!esl_abc_CIsGap(msa->abc, msa->rf[apos])) 
	    matassign[apos+1] = 1;
	}
      /* make fake tracebacks for each seq */
      CP9_fake_tracebacks(msa, matassign, &cp9_tr);
      
      /* build model from tracebacks (code from HMMER's modelmakers.c::matassign2hmm() */
      for (i = 0; i < msa->nseq; i++) {
	CP9TraceCount(orig_hmm, sq[i]->dsq, msa->wgt[i], cp9_tr[i]);
      }
      nsampled += msa_nseq;
    }
  
  /*Next, renormalize the sub_hmm and logoddisfy it */
  CPlan9Renormalize(sub_hmm);
  CP9Logoddsify(sub_hmm);

  /* clean up from previous MSA */
  esl_msa_Destroy(msa);
  free(matassign);
  free(useme);
  for (i = 0; i < msa_nseq; i++)
    {
      CP9FreeTrace(cp9_tr[i]);
      FreeParsetree(tr[i]);
      esl_sq_Destroy(sq[i]);
    }
  free(cp9_tr);
  /**************************************************/
  
  printf("PRINTING SAMPLED ORIG HMM PARAMS:\n");
  debug_print_cp9_params(stdout, orig_hmm, TRUE);
  printf("DONE PRINTING SAMPLED ORIG HMM PARAMS:\n");
  

  printf("PRINTING SAMPLED SUB HMM PARAMS:\n");
  debug_print_cp9_params(stdout, sub_hmm, TRUE);
  printf("DONE PRINTING SAMPLED SUB HMM PARAMS:\n");
  
  FreeCPlan9(orig_hmm);
  FreeCPlan9(sub_hmm);
  return TRUE;

 ERROR:
  cm_Fail("Memory allocation error.");
  return FALSE; /* never reached */
}

#if 0
/*
 * Function: cm_emit_seqs_to_aln_above_cutoff()
 * Date:     EPN, Mon Sep 10 17:31:36 2007
 *
 * Purpose:  Create a seqs_to_aln object by generating sequences
 *           from a CM. Only accept sequences that have a CM hit
 *           within them above a bit score cutoff.
 *
 * Args:     go              - getopts
 *           cfg             - cmcalibrate's configuration
 *           cm              - CM to emit from
 *           nseq            - number of seqs to emit
 *           cutoff          - bit score cutoff 
 *
 * Returns:  Ptr to a newly allocated seqs_to_aln object with nseq sequences to align.
 *           Dies immediately on failure with informative error message.
 */
seqs_to_aln_t *cm_emit_seqs_to_aln_above_cutoff(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int nseq)
{
  int status;
  seqs_to_aln_t *seqs_to_aln = NULL;
  char *name = NULL;
  int namelen;
  int L;
  int i;
  int do_cyk = TRUE;
  Parsetree_t *tr = NULL;

  if(cm->dmin == NULL || cm->dmax == NULL) cm_Fail("cm_emit_seqs_to_aln_above_cutoff(), dmin, dmax are NULL.");
  if(cm->search_opts & CM_SEARCH_NOQDB)    cm_Fail("cm_emit_seqs_to_aln_above_cutoff(), search opt NOQDB enabled.");

  seqs_to_aln = CreateSeqsToAln(nseq, FALSE);

  namelen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
  if(cm->name != NULL) namelen += strlen(cm->name) + 1;
  ESL_ALLOC(name, sizeof(char) * namelen);

  while(i < nseq)
    {
      if(cm->name != NULL) sprintf(name, "%s-%d", cm->name, i+1);
      else                 sprintf(name, "%d-%d", cfg->ncm-1, i+1);
      L = 0; 
      EmitParsetree(cm, cfg->r, name, TRUE, &tr, &(seqs_to_aln->sq[i]), &L);
      while(L == 0) { FreeParsetree(tr); esl_sq_Destroy(seqs_to_aln->sq[i]); EmitParsetree(cm, cfg->r, name, TRUE, &tr, &(seqs_to_aln->sq[i]), &L); }
      p = cfg->cmstatsA[(ncm-1)]->gc2p[(get_gc_comp(seqs_to_aln->sq[i], 1, L))]; /* in get_gc_comp() should be i and j of best hit */

      sc = ParsetreeScore(cm, tr, seqs_to_aln->sq[i]->dsq, FALSE); 
      FreeParsetree(tr);
      if(sc > cfg->cutoffA[p]) { i++; continue; }

      /* If we get here, parse score is not above cfg->cutoffA[p], we want to determine if 
       * this sequence has a hit in it above the cfg->cutoffA[p] as quickly as possible. 
       *
       * Stage 1: HMM banded search tau = 1e-2
       * Stage 2: HMM banded search with scanning bands, tau = 1e-10
       * Stage 3: QDB search (CYK or inside), beta = cm->beta (should be default beta)
       *
       * If we find a hit > cfg->cutoffA[p] at any stage, we accept the seq, increment i and move on.
       * We don't do a full non-banded parse to ensure that we don't exceed the cfg->cutoffA[p], 
       * because QDB is on in cmsearch by default.
       */

      /* stage 1 */
      cm->search_opts |= CM_SEARCH_HBANDED;
      cm->tau = 0.01;
      if((sc = actually_search_target(cm, seqs_to_aln->sq[i]->dsq, 1, L, 0., 0., NULL, FALSE, FALSE, FALSE, NULL, FALSE)) > cfg->cutoffA[p]) 
	{ i++; break; }
      s1_np++;
      /* stage 2 */
      cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
      cm->tau = 1e-10;
      if((sc = actually_search_target(cm, seqs_to_aln->sq[i]->dsq, 1, L, 0., 0., NULL, FALSE, FALSE, FALSE, NULL, FALSE)) > cfg->cutoffA[p]) 
	{ i++; break; }
      s2_np++;
      /* stage 3 */
      cm->search_opts &= ~CM_SEARCH_HBANDED;
      cm->search_opts &= ~CM_SEARCH_HBANDED;
      if((sc = search_target_cm_calibration(cm, seqs_to_aln->sq[i]->dsq, cm->dmin, cm->dmax, 1, seqs_to_aln->sq[i]->n, cm->W, NULL)) > cfg->cutoffA[p])
	{ i++; break; }
      s3_np++;
      if(s3_np > (1000 * nseq)) cm_Fail("cm_emit_seqs_to_aln_above_cutoff(), wanted %d seqs above cutoff: %d bits, reached limit of %d seqs\n", nseq, cfg->cutoffA[p], (1000 * nseq));

      /* didn't pass */
      esl_sq_Destroy(seqs_to_aln->sq[i]);
    }

  seqs_to_aln->nseq = nseq;

  free(name);
  return seqs_to_aln;

 ERROR:
  cm_Fail("memory allocation error");
  return NULL;
}


/* Function: search_target_cm_calibration() based on bandcyk.c:CYKBandedScan()
 * Date:     EPN, Sun Sep  9 19:05:07 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm. If bands are NULL, reverts to non-banded
 *           (scancyk.c:CYKScan()). 
 *
 *           Special local cmcalibrate function that only cares about
 *           collecting the best score at each state for the sequence.
 *
 * Args:     cm        - the covariance model
 *           dsq       - the digitized sequence
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           ret_vsc  - RETURN: [0..v..M-1] best score at each state v
 *
 * Returns:  score of best overall hit (vsc[0]).
 *           dies immediately if some error occurs.
 */
float 
search_target_cm_calibration(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float **ret_vsc)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int    ***ialpha;             /* Inside DP score matrix, [v][j][d] */
  float    *vsc;           /* best score for each state (float) */
  float    *ivsc;          /* best score for each state (int, only used if do_inside) */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_w;  	        /* offset j for state w */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       dn;                 /* temporary value for min d in for loops */
  int       dx;                 /* temporary value for max d in for loops */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */

  int       do_cyk    = FALSE;  /* TRUE: do cyk; FALSE: do_inside */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  float     ret_val;

  /* determine if we're doing cyk/inside banded/non-banded */
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) do_cyk    = TRUE;
  if(dmin != NULL && dmax != NULL)           do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 

  ESL_ALLOC(vsc, sizeof(float) * cm->M);
  esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);

  if(do_cyk) { 

    /*****************************************************************
     * alpha allocations.
     * The scanning matrix is indexed [v][j][d]. 
     *    v ranges from 0..M-1 over states in the model.
     *    j takes values 0 or 1: only the previous (prv) or current (cur) row
     *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
     *      deck in memory, and j ranges from 0..W, and yes it must be square
     *      because we'll use a rolling pointer trick thru it
     *    d ranges from 0..W over subsequence lengths.
     * Note that E memory is shared: all E decks point at M-1 deck.
     *****************************************************************/
    ESL_ALLOC(alpha, (sizeof(float **) * cm->M));
    for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
      if (cm->stid[v] == BEGL_S)
	{
	  ESL_ALLOC(alpha[v], (sizeof(float *) * (W+1)));
	  for (j = 0; j <= W; j++)
	    ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
	}
      else if (cm->sttype[v] == E_st && v < cm->M-1) 
	alpha[v] = alpha[cm->M-1];
      else 
	{
	  ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	  for (j = 0; j < 2; j++) 
	    ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
	}
    }

    /*****************************************************************
     * alpha initializations.
     * We initialize on d=0, subsequences of length 0; these are
     * j-independent. Any generating state (P,L,R) is impossible on d=0.
     * E=0 for d=0. B,S,D must be calculated. 
     * Also, for MP, d=1 is impossible.
     * Also, for E, all d>0 are impossible.
     *
     * and, for banding: any cell outside our bands is impossible.
     * These inits are never changed in the recursion, so even with the
     * rolling, matrix face reuse strategy, this works.
     *****************************************************************/ 
    for (v = cm->M-1; v >= 0; v--)
      {
	alpha[v][0][0] = IMPOSSIBLE;

	if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
	else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = IMPOSSIBLE;
	else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	  {
	    y = cm->cfirst[v];
	    alpha[v][0][0] = cm->endsc[v];
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	      alpha[v][0][0] = ESL_MAX(alpha[v][0][0], (alpha[y+yoffset][0][0] + cm->tsc[v][yoffset]));
	    /* ...we don't bother to look at local alignment starts here... */
	    alpha[v][0][0] = ESL_MAX(alpha[v][0][0], IMPOSSIBLE);
	  }
	else if (cm->sttype[v] == B_st) 
	  {
	    w = cm->cfirst[v];
	    y = cm->cnum[v];
	    alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	  }

	alpha[v][1][0] = alpha[v][0][0];
	if (cm->stid[v] == BEGL_S) 
	  for (j = 2; j <= W; j++) 
	    alpha[v][j][0] = alpha[v][0][0];
      }
    /* Impose the bands.
     *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
     *    the first loop will be skipped, the second initializes the E states.)
     */
    if(do_banded) { 
      for (v = 0; v < cm->M; v++) {
	jmax = (cm->stid[v] == BEGL_S) ? W : 1;
	for (d = 0; d < dmin[v] && d <=W; d++) 
	  for(j = 0; j <= jmax; j++)
	    alpha[v][j][d] = IMPOSSIBLE;
	for (d = dmax[v]+1; d <= W;      d++) 
	  for(j = 0; j <= jmax; j++)
	    alpha[v][j][d] = IMPOSSIBLE;
      }
    }

    /* The main loop: scan the sequence from position i0 to j0.
     */
    for (j = i0; j <= j0; j++) 
      {
	cur = j%2;
	prv = (j-1)%2;
	for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	  {
	    /* determine min/max d we're allowing for this state v and this position j */
	    if(do_banded) { 
	      dn = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	      dx = ESL_MIN((j-i0+1), dmax[v]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { 
	      dn = (cm->sttype[v] == MP_st) ? 2 : 1;
	      dx = ESL_MIN((j-i0+1), W); 
	    }

	    jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	    jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	    sd   = StateDelta(cm->sttype[v]);

	    if(cm->sttype[v] == B_st) {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      for (d = dn; d <= dx; d++) {
		/* k is the length of the right fragment */
		/* Careful, make sure k is consistent with bands in state w and state y. */
		if(do_banded) {
		  kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		  kmin = ESL_MAX(kmin, 0);
		  kmax = ESL_MIN(dmax[y], (d-dmin[w]));
		}
		else { kmin = 0; kmax = d; }

		alpha[v][jp_v][d] = ESL_MAX(IMPOSSIBLE, cm->endsc[v] + (cm->el_selfsc * (d - sd)));
		for (k = kmin; k <= kmax; k++) { 
		  jp_w = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      alpha[v][jp_v][d] = ESL_MAX(alpha[v][jp_v][d], (alpha[w][jp_w][d-k] + alpha[y][jp_y][k]));
		}
		vsc[v] = ESL_MAX(vsc[v], alpha[v][jp_v][d]);
	      }
	    }
	    else { /* if cm->sttype[v] != B_st */
	      for (d = dn; d <= dx; d++) {
		alpha[v][jp_v][d] = ESL_MAX (IMPOSSIBLE, (cm->endsc[v] + (cm->el_selfsc * (d - sd))));
		y = cm->cfirst[v];
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  alpha[v][jp_v][d] = ESL_MAX (alpha[v][jp_v][d], (alpha[y+yoffset][jp_y][d - sd] + cm->tsc[v][yoffset]));
		
		/* add in emission score, if any */
		i = j-d+1;
		switch (cm->sttype[v]) {
		case MP_st: 
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][jp_v][d] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  break;
		case ML_st:
		case IL_st:
		  alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  break;
		case MR_st:
		case IR_st:
		  alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  break;
		} /* end of switch */
		vsc[v] = ESL_MAX(vsc[v], alpha[v][jp_v][d]);
	      } /* end of d = dn; d <= dx; d++ */
	    } /* end of else (v != B_st) */
	  } /*loop over decks v>0 */

	/* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
	 * 
	 * If local begins are off, the hit must be rooted at v=0.
	 * With local begins on, the hit is rooted at the second state in
	 * the traceback (e.g. after 0), the internal entry point. Divide & conquer
	 * can only handle this if it's a non-insert state; this is guaranteed
	 * by the way local alignment is parameterized (other transitions are
	 * -INFTY), which is probably a little too fragile of a method. 
	 */

	/* determine min/max d we're allowing for the root state and this position j */
	if(do_banded) { 
	  dn = ESL_MAX(dmin[0], 1); 
	  dx = ESL_MIN((j-i0+1), dmax[0]); 
	  dx = ESL_MIN(dx, W);
	}
	else { 
	  dn = 1; 
	  dx = ESL_MIN((j-i0+1), W); 
	}
	jp_v = cur;
	jp_y = cur;
	for (d = dn; d <= dx; d++) {
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = ESL_MAX(IMPOSSIBLE, alpha[y][cur][d] + cm->tsc[0][0]);
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	    alpha[0][cur][d] = ESL_MAX (alpha[0][cur][d], (alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]));
	  vsc[0] = ESL_MAX(vsc[0], alpha[0][cur][d]);
	}
	
	if (cm->flags & CM_LOCAL_BEGIN) {
	  for (y = 1; y < cm->M; y++) {
	    if(do_banded) {
	      dn = (cm->sttype[y] == MP_st) ? ESL_MAX(dmin[y], 2) : ESL_MAX(dmin[y], 1); 
	      dn = ESL_MAX(dn, dmin[0]);
	      dx = ESL_MIN((j-i0+1), dmax[y]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { 
	      dn = 1; 
	      dx = ESL_MIN((j-i0+1), W); 
	    }
	    jp_y = (cm->stid[y] == BEGL_S) ? (j % (W+1)) : cur;
	    for (d = dn; d <= dx; d++) {
	      alpha[0][cur][d] = ESL_MAX(alpha[0][cur][d], alpha[y][jp_y][d] + cm->beginsc[y]);
	      vsc[0] = ESL_MAX(vsc[0], alpha[0][cur][d]);
	    }
	  }
	}
      } /* end loop over end positions j */
    /* free alpha, we only care about vsc 
     */
    for (v = 0; v < cm->M; v++) 
      {
	if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	  for (j = 0; j <= W; j++) free(alpha[v][j]);
	  free(alpha[v]);
	} else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	  continue;
	} else {
	  free(alpha[v][0]);
	  free(alpha[v][1]);
	  free(alpha[v]);
	}
      }
    free(alpha);
  }
  /*********************
   * end of if(do_cyk) *
   *********************/
  else { /* ! do_cyk, do_inside, with scaled int log odds scores instead of floats */

    ESL_ALLOC(ivsc, sizeof(int) * cm->M);
    esl_vec_FSet(ivsc, cm->M, -INFTY);
    
    /* ialpha allocations. (see comments for do_cyk section */ 
    ESL_ALLOC(ialpha, sizeof(int **) * cm->M);
    for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(ialpha[v], sizeof(int *) * (W+1));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(ialpha[v][j], sizeof(int) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      ialpha[v] = ialpha[cm->M-1];
    else 
      {
	ESL_ALLOC(ialpha[v], sizeof(int *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(ialpha[v][j], sizeof(int) * (W+1));
      }
    }
    /* ialpha initializations. (see comments for do_cyk section */
    for (v = cm->M-1; v >= 0; v--)  {
	ialpha[v][0][0] = -INFTY;
	if      (cm->sttype[v] == E_st)  ialpha[v][0][0] = 0;
	else if (cm->sttype[v] == MP_st) ialpha[v][0][1] = ialpha[v][1][1] = -INFTY;
	else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	  {
	    y = cm->cfirst[v];
	    ialpha[v][0][0] = cm->iendsc[v];
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	      ialpha[v][0][0] = ILogsum(ialpha[v][0][0], (ialpha[y+yoffset][0][0] 
							+ cm->itsc[v][yoffset]));
	    /* ...we don't bother to look at local alignment starts here... */
	    /* ! */
	    if (ialpha[v][0][0] < -INFTY) ialpha[v][0][0] = -INFTY;	
	  }
	else if (cm->sttype[v] == B_st)  {
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  ialpha[v][0][0] = ialpha[w][0][0] + ialpha[y][0][0]; 
	}
      ialpha[v][1][0] = ialpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  ialpha[v][j][0] = ialpha[v][0][0];
    }
    /* Impose the bands.
     *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
     *    the first loop will be skipped, the second initializes the E states.)
     */
    if(do_banded) {
      for (v = 0; v < cm->M; v++) {
	if(cm->stid[v] == BEGL_S) jmax = W; 
	else jmax = 1;
	
	dx = ESL_MIN(dmin[v], W);
	for (d = 0; d < dx; d++) 
	  for(j = 0; j <= jmax; j++)
	    ialpha[v][j][d] = -INFTY;
	
	for (d = dmax[v]+1; d <= W;      d++) 
	  for(j = 0; j <= jmax; j++)
	    ialpha[v][j][d] = -INFTY;
      }
    }

    /* The main loop: scan the sequence from position i0 to j0.
     */
    for (j = i0; j <= j0; j++) 
      {
	cur = j%2;
	prv = (j-1)%2;
	for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	  {
	    /* determine min/max d we're allowing for this state v and this position j */
	    if(do_banded) { 
	      dn = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	      dx = ESL_MIN((j-i0+1), dmax[v]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { 
	      dn = (cm->sttype[v] == MP_st) ? 2 : 1;
	      dx = ESL_MIN((j-i0+1), W); 
	    }

	    jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	    jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	    sd   = StateDelta(cm->sttype[v]);

	    if(cm->sttype[v] == B_st) {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      for (d = dn; d <= dx; d++) {
		/* k is the length of the right fragment */
		/* Careful, make sure k is consistent with bands in state w and state y. */
		if(do_banded) {
		  kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		  kmin = ESL_MAX(kmin, 0);
		  kmax = ESL_MIN(dmax[y], (d-dmin[w]));
		}
		else { kmin = 0; kmax = d; }

		ialpha[v][jp_v][d] = ESL_MAX(-INFTY, cm->iendsc[v] + (cm->iel_selfsc * (d - sd)));
		for (k = kmin; k <= kmax; k++) { 
		  jp_w = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      ialpha[v][jp_v][d] = ESL_MAX(ialpha[v][jp_v][d], (ialpha[w][jp_w][d-k] + ialpha[y][jp_y][k]));
		}
		ivsc[v] = ESL_MAX(ivsc[v], ialpha[v][jp_v][d]);
	      }
	    }
	    else { /* if cm->sttype[v] != B_st */
	      for (d = dn; d <= dx; d++) {
		ialpha[v][jp_v][d] = ESL_MAX (-INFTY, (cm->iendsc[v] + (cm->iel_selfsc * (d - sd))));
		y = cm->cfirst[v];
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  ialpha[v][jp_v][d] = ESL_MAX (ialpha[v][jp_v][d], (ialpha[y+yoffset][jp_y][d - sd] + cm->itsc[v][yoffset]));
		
		/* add in emission score, if any */
		i = j-d+1;
		switch (cm->sttype[v]) {
		case MP_st: 
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    ialpha[v][jp_v][d] += cm->iesc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		  else
		    ialpha[v][cur][d] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);
		  break;
		case ML_st:
		case IL_st:
		  ialpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);
		  break;
		case MR_st:
		case IR_st:
		  ialpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);
		  break;
		} /* end of switch */
		ivsc[v] = ESL_MAX(ivsc[v], ialpha[v][jp_v][d]);
	      } /* end of d = dn; d <= dx; d++ */
	    } /* end of else (v != B_st) */
	  } /*loop over decks v>0 */
	/* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
	 * 
	 * If local begins are off, the hit must be rooted at v=0.
	 * With local begins on, the hit is rooted at the second state in
	 * the traceback (e.g. after 0), the internal entry point. Divide & conquer
	 * can only handle this if it's a non-insert state; this is guaranteed
	 * by the way local alignment is parameterized (other transitions are
	 * -INFTY), which is probably a little too fragile of a method. 
	 */

	/* determine min/max d we're allowing for the root state and this position j */
	if(do_banded) { 
	  dn = ESL_MAX(dmin[0], 1); 
	  dx = ESL_MIN((j-i0+1), dmax[0]); 
	  dx = ESL_MIN(dx, W);
	}
	else { 
	  dn = 1; 
	  dx = ESL_MIN((j-i0+1), W); 
	}
	jp_v = cur;
	jp_y = cur;
	for (d = dn; d <= dx; d++) {
	  y = cm->cfirst[0];
	  ialpha[0][cur][d] = ESL_MAX(IMPOSSIBLE, ialpha[y][cur][d] + cm->itsc[0][0]);
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	    ialpha[0][cur][d] = ESL_MAX (ialpha[0][cur][d], (ialpha[y+yoffset][cur][d] + cm->itsc[0][yoffset]));
	  ivsc[0] = ESL_MAX(ivsc[0], ialpha[0][cur][d]);
	}
	
	if (cm->flags & CM_LOCAL_BEGIN) {
	  for (y = 1; y < cm->M; y++) {
	    if(do_banded) {
	      dn = ESL_MAX(1,  dmin[y]);
	      dn = ESL_MAX(dn, dmin[0]);
	      dx = ESL_MIN((j-i0+1), dmax[y]); 
	      dx = ESL_MIN(dx, W);
	    }
	    else { dn = 1; dx = W; }
	    jp_y = (cm->stid[y] == BEGL_S) ? (j % (W+1)) : cur;
	    for (d = dn; d <= dx; d++) {
	      ialpha[0][cur][d] = ESL_MAX(ialpha[0][cur][d], ialpha[y][jp_y][d] + cm->ibeginsc[y]);
	      ivsc[0] = ESL_MAX(ivsc[0], ialpha[0][cur][d]);
	    }
	  }
	}
      } /* end loop over end positions j */
    /* free ialpha, we only care about ivsc 
     */
    for (v = 0; v < cm->M; v++) 
      {
	if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	  for (j = 0; j <= W; j++) free(ialpha[v][j]);
	  free(ialpha[v]);
	} else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	  continue;
	} else {
	  free(ialpha[v][0]);
	  free(ialpha[v][1]);
	  free(ialpha[v]);
	}
      }
    free(ialpha);
    /* convert ivsc to floats in vsc */
    ESL_ALLOC(vsc, sizeof(float) * cm->M);
    for(v = 0; v < cm->M; v++)
      vsc[v] = Scorify(ivsc[v]);
    free(ivsc);
  }
  /**************************
   * end of else (do_inside)
   **************************/

  ret_val = vsc[0];
  if (ret_vsc != NULL) *ret_vsc = vsc;
  else free(vsc);
  
  return ret_val;

  ERROR:
    cm_Fail("Memory allocation error.\n");
    return 0.; /* NEVERREACHED */
}
#endif

/* EPN, Mon Dec 10 13:15:32 2007, old process_workunit() function */
#if 0

/* Function: process_workunit()
 * Date:     EPN, Mon Sep 10 16:55:09 2007
 *
 * Purpose:  A work unit consists of a CM, a int specifying a number of 
 *           sequences <nseq>, and a flag indicated how to generate those
 *           sequences. The job is to generate <nseq> sequences and search
 *           them with a CM and/or CP9, saving scores, which are returned.
 *
 *           This function can be run in 1 of 3 modes, determined by the
 *           status of the input variables:
 *         
 *           Mode 1. Gumbel calculation for CM. 
 *           <emit_from_cm> is FALSE, <ret_vscAA> != NULL, <ret_cp9scA> == NULL
 *           Emit randomly and search only with the CM. <ret_vscAA> is filled
 *           with the best CM score at each state for each sequence.
 *
 *           Mode 2. Gumbel calculation for CP9.
 *           <emit_from_cm> is FALSE, <ret_vscAA> == NULL, <ret_cp9scA> != NULL
 *           Emit randomly and search only with the CP9. <ret_cp9scA> is filled
 *           with the best CP9 score for each sequence.
 *
 *           Mode 3. Scores will eventually be used to determine filter thresholds.
 *           <emit_from_cm> is TRUE, <ret_vscAA> != NULL, <ret_cp9scA> != NULL, 
 *           <ret_other_cp9scA> != NULL.
 *           Emit from the CM (which is already configured how we want it). Search
 *           with the CM first, then with the CP9 twice, first w/Viterbi then w/Forward
 *           <ret_vscAA> filled with the best CM score at each state for each sequence,
 *           <ret_cp9scA> filled with the best CP9 Viterbi score for each sequence,
 *           <ret_other_cp9scA> filled with the best CP9 Forward score for each sequence,
 *           Importantly, in this mode, each sequence must have a NON-BANDED CM scan 
 *           (either CYK or Inside) hit above a given cutoff. That cutoff is given
 *           as a bit score in cfg->cutoffA[p], where p is the partition for the
 *           sequence (p is determined in get_cmemit_dsq() called from this function). 
 *           Sequences that have no hit better than cutoff are not accepted, (they're
 *           rejected and not searched, and another seq is emitted). The cutoff[p] 
 *           value is assumed to be already set before this function is entered.
 *
 *           The ability to run 3 different modes complicates the code a bit,
 *           but I prefered it to making 2 separate functions b/c a significant
 *           part of those 2 functions would have identical code. Also it makes
 *           the MPI implementation a bit easier because the workers can always
 *           call this function, whether they're calcing Gumbels or filter thresholds.
 *
 * Args:     go           - getopts
 *           cfg          - cmcalibrate's configuration
 *           errbuf       - for writing out error messages
 *           cm           - the CM (already configured as we want it)
 *           nseq         - number of seqs to generate
 *           emit_from_cm - TRUE to emit from CM; FALSE emit random 
 *           ret_vscAA    - RETURN: [0..v..cm->M-1][0..nseq-1] best 
 *                                  score at each state v for each seq
 *           ret_cp9scA   - RETURN: [0..nseq-1] best CP9 score for each seq
 *                                  if (emit_from_cm) these will be Viterbi scores, else 
 *                                  could be Viterbi or Forward
 *           ret_other_cp9scA - RETURN: [0..nseq-1] best CP9 score for each seq
 *                                      if (emit_from_cm) these will be Forward scores, else 
 *                                      it will == NULL
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
		 int emit_from_cm, float ***ret_vscAA, float **ret_cp9scA, float **ret_other_cp9scA)
{
  int            status;
  int            mode; /* 1, 2, or 3, determined by status of input args, as explained in 'Purpose' above. */
  float        **vscAA        = NULL;  /* [0..v..cm->M-1][0..i..nseq-1] best CM score for each state, each seq */
  float         *cur_vscA     = NULL;  /* [0..v..cm->M-1]               best CM score for each state cur seq */
  float         *cp9scA       = NULL;  /*                [0..i..nseq-1] best CP9 score for each seq, 
					*                               if (emit_from_cm) these will be Viterbi scores,
                                        *                               else they could be Viterbi or Forward */
  float         *other_cp9scA = NULL;  /*                [0..i..nseq-1] best CP9 Forward score for each seq 
					*                               only if (emit_from_cm), else stays NULL */
  double        *dnull        = NULL; /* double version of cm->null, for generating random seqs */
  int            p;                   /* what partition we're in, not used unless emit_from_cm = TRUE */
  int            i;
  int            v;
  int            L;
  int            nfailed = 0;
  Parsetree_t   *tr;
  ESL_DSQ       *dsq;
  float          sc;
  float         *fwd_sc_ptr;

  /* determine mode, and enforce mode-specific contract */
  if     (ret_vscAA != NULL && ret_cp9scA == NULL)     mode = 1; /* calcing CM  gumbel stats */
  else if(ret_vscAA == NULL && ret_cp9scA != NULL)     mode = 2; /* calcing CP9 gumbel stats */
  else if(ret_vscAA != NULL && ret_cp9scA != NULL && 
	  ret_other_cp9scA != NULL && emit_from_cm) mode = 3; /* collecting filter threshold stats */
  else ESL_FAIL(eslEINCOMPAT, errbuf, "can't determine mode in process_workunit.");
  if(emit_from_cm && mode != 3) ESL_FAIL(eslEINCOMPAT, errbuf, "emit_from_cm is TRUE, but mode is: %d (should be 3)\n", mode);

  ESL_DPRINTF1(("in process_workunit nseq: %d mode: %d\n", nseq, mode));

  int do_cyk     = FALSE;
  int do_inside  = FALSE;
  int do_viterbi = FALSE;
  int do_forward = FALSE;
  /* determine algs we'll use and allocate the score arrays we'll pass back */
  if(mode == 1 || mode == 3) {
    if(cm->search_opts & CM_SEARCH_INSIDE) do_inside = TRUE;
    else                                   do_cyk    = TRUE;
    ESL_ALLOC(vscAA, sizeof(float *) * cm->M);
    for(v = 0; v < cm->M; v++) ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    ESL_ALLOC(cur_vscA, sizeof(float) * cm->M);
  }
  if(mode == 2) {
    if(cm->search_opts & CM_SEARCH_HMMVITERBI) do_viterbi = TRUE;
    if(cm->search_opts & CM_SEARCH_HMMFORWARD) do_forward = TRUE;
    if((do_viterbi + do_forward) > 1) ESL_FAIL(eslEINVAL, errbuf, "process_workunit, mode 2, and cm->search_opts CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags both raised.");
    ESL_ALLOC(cp9scA, sizeof(float) * nseq); /* will hold Viterbi or Forward scores */
  }
  if(mode == 3) {
    do_viterbi = do_forward = TRUE;
    ESL_ALLOC(cp9scA,       sizeof(float) * nseq); /* will hold Viterbi scores */
    ESL_ALLOC(other_cp9scA, sizeof(float) * nseq); /* will hold Forward scores */
  }
  ESL_DPRINTF1(("do_cyk:     %d\ndo_inside:  %d\ndo_viterbi: %d\ndo_forward: %d\n", do_cyk, do_inside, do_viterbi, do_forward)); 
  
  /* fill dnull, a double version of cm->null, but only if we're going to need it to generate random seqs */
  if(!emit_from_cm && cfg->pgc_freq == NULL) {
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);    
  }
  
  /* generate dsqs one at a time and collect best CM scores at each state and/or best overall CP9 score */
  for(i = 0; i < nseq; i++) {
    if(emit_from_cm) { /* if emit_from_cm == TRUE, use_cm == TRUE */
      if(nfailed > 1000 * nseq) { cm_Fail("Max number of failures (%d) reached while trying to emit %d seqs.\n", nfailed, nseq); }
      dsq = get_cmemit_dsq(cfg, cm, &L, &p, &tr);
      /* we only want to use emitted seqs with a sc > cutoff, cm_find_hit_above_cutoff returns false if no such hit exists in dsq */
      if((status = cm_find_hit_above_cutoff(go, cfg, errbuf, cm, dsq, tr, L, cfg->cutoffA[p], &sc)) != eslOK) return status;
      while(sc < cfg->cutoffA[p]) { 
	free(dsq); 	
	/* parsetree tr is freed in cm_find_hit_above_cutoff() */
	dsq = get_cmemit_dsq(cfg, cm, &L, &p, &tr);
	nfailed++;
	if((status = cm_find_hit_above_cutoff(go, cfg, errbuf, cm, dsq, tr, L, cfg->cutoffA[p], &sc)) != eslOK) return status;
      }
      ESL_DPRINTF1(("i: %d nfailed: %d\n", i, nfailed));
    }
    else { 
      dsq = get_random_dsq(cfg, cm, dnull, cm->W*2); 
      L = cm->W*2; 
    }
    /* if nec, search with CM */
    if (do_cyk)    if((status = FastCYKScan    (cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc)) != eslOK) return status;
    if (do_inside) if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc)) != eslOK) return status;
    /* if nec, search with CP9 */
    if (do_viterbi) 
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       NULL,   /* don't want traces back */
			       &(cp9scA[i]))) != eslOK) return status;
    if (do_forward) {
      if      (mode == 2) fwd_sc_ptr = &(cp9scA[i]);       /* fill cp9scA[i] */
      else if (mode == 3) fwd_sc_ptr = &(other_cp9scA[i]); /* fill other_cp9scA[i] */
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       fwd_sc_ptr)) != eslOK) return status;
    }
    free(dsq);
    if (cur_vscA != NULL) /* will be NULL if do_cyk == do_inside == FALSE (mode 2) */
      for(v = 0; v < cm->M; v++) vscAA[v][i] = cur_vscA[v];
    free(cur_vscA);
  }

  if(dnull != NULL) free(dnull);
  if(ret_vscAA  != NULL)       *ret_vscAA  = vscAA;
  if(ret_cp9scA != NULL)       *ret_cp9scA = cp9scA;
  if(ret_other_cp9scA != NULL) *ret_other_cp9scA = other_cp9scA;
  return eslOK;

 ERROR:
  return status;
}


/* Function: calc_best_filter()
 * Date:     EPN, Thu Nov  1 15:05:03 2007
 *
 * Purpose:  Given a CM and scores for a CP9 and CM scan of target seqs
 *           determine the best filter we could use, either an HMM only
 *           or a hybrid scan with >= 1 sub CM roots.
 *            
 * Returns:  eslOK on success;
 *           Dies immediately on an error.
 */
int
calc_best_filter(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, float **fil_vscAA, float *fil_vit_cp9scA, float *fil_fwd_cp9scA)
{
  int    status;
  int    v;
  float  *sorted_fil_vit_cp9scA;
  float  *sorted_fil_fwd_cp9scA;
  float **sorted_fil_vscAA;
  float **sorted_fil_EAA;
  int    filN  = esl_opt_GetInteger(go, "--filN");
  int    Fidx;
  float  vit_sc, fwd_sc, sc;
  float  E;
  float  fil_calcs;
  float  surv_calcs;
  float  fil_plus_surv_calcs;
  float  nonfil_calcs;
  float  spdup;
  int    i;
  int    cmi = cfg->ncm-1;
  int    cp9_vit_mode, cp9_fwd_mode;
  
  float F = esl_opt_GetReal(go, "--F");
  Fidx  = (int) ((1. - F) * (float) filN);

  if(cfg->cmstatsA[cfg->ncm-1]->np != 1) cm_Fail("calc_sub_filter_sets(), not yet implemented for multiple partitions.\nYou'll need to keep track of partition of each sequence OR\nstore E-values not scores inside process_workunit.");

  /* Determine the predicted CP9 filter speedup */
  ESL_ALLOC(sorted_fil_vit_cp9scA, sizeof(float) * filN);
  esl_vec_FCopy(fil_vit_cp9scA, filN, sorted_fil_vit_cp9scA); 
  esl_vec_FSortIncreasing(sorted_fil_vit_cp9scA, filN);
  vit_sc = sorted_fil_vit_cp9scA[Fidx];

  ESL_ALLOC(sorted_fil_fwd_cp9scA, sizeof(float) * filN);
  esl_vec_FCopy(fil_fwd_cp9scA, filN, sorted_fil_fwd_cp9scA); 
  esl_vec_FSortIncreasing(sorted_fil_fwd_cp9scA, filN);
  fwd_sc = sorted_fil_fwd_cp9scA[Fidx];

  printf("\n\n***********************************\n\n");
  for(i = 0; i < filN; i++)
    printf("HMM i: %4d vit sc: %10.4f fwd sc: %10.4f\n", i, sorted_fil_vit_cp9scA[i], sorted_fil_fwd_cp9scA[i]);
  printf("***********************************\n\n");

  cp9_vit_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? GUM_CP9_LV : GUM_CP9_GV;
  cp9_fwd_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? GUM_CP9_LF : GUM_CP9_GF;

  /* print out predicted speed up with Viterbi filter */
  /* E is expected number of hits for db of length 2 * cm->W */
  /* EPN, Sun Dec  9 16:40:39 2007
   * idea: calculate E for each partition, then take weighted average, assuming each GC segment is equally likely (or some other weighting) */
  E  = RJK_ExtremeValueE(sc, cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->mu, cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->lambda);
  fil_calcs  = cfg->hsi->full_cp9_ncalcs;
  surv_calcs = (E / cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->L) * /* calcs are in units of millions of dp calcs per residue */
    (cfg->avglen[0]) * cfg->hsi->full_cm_ncalcs; /* cfg->avglen[0] is average length of subseq for subtree rooted at v==0, for current gumbel mode configuration */
  fil_plus_surv_calcs = fil_calcs + surv_calcs;
  nonfil_calcs = cfg->hsi->full_cm_ncalcs;
  spdup = nonfil_calcs / fil_plus_surv_calcs; 
  printf("HMM(vit) sc: %10.4f E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", vit_sc, E, fil_calcs, surv_calcs, fil_plus_surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup);

  /* print out predicted speed up with Forward filter */
  /* EPN, Sun Dec  9 16:40:39 2007
   * idea: calculate E for each partition, then take weighted average, assuming each GC segment is equally likely (or some other weighting) */
  E  = RJK_ExtremeValueE(sc, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->mu, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->lambda);
  fil_calcs  = cfg->hsi->full_cp9_ncalcs;
  surv_calcs = (E / cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->L) * /* calcs are in units of millions of dp calcs per residue */
    (cfg->avglen[0]) * cfg->hsi->full_cm_ncalcs; /* cfg->avglen[0] is average length of subseq for subtree rooted at v==0, for current gumbel mode configuration */
  fil_plus_surv_calcs = fil_calcs + surv_calcs;
  nonfil_calcs = cfg->hsi->full_cm_ncalcs;
  spdup = nonfil_calcs / (fil_plus_surv_calcs * 2.); /* the '* 2.' is to correct for fact that Forward is about 2X slower than Viterbi, due to the logsum() instead of ESL_MAX() calculations */
  printf("HMM(fwd) sc: %10.4f E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", fwd_sc, E, fil_calcs, surv_calcs, fil_plus_surv_calcs, 2.*fil_plus_surv_calcs, nonfil_calcs, spdup);

  exit(1);
  /*****************TEMPORARY PRINTF BEGIN***************************/
  ESL_ALLOC(sorted_fil_vscAA, sizeof(float *) * cm->M);
  ESL_ALLOC(sorted_fil_EAA, sizeof(float *) * cm->M);

  /*printf("\n\n***********************************\nvscAA[0] scores:\n");
  for(i = 0; i < filN; i++)
    printf("i: %4d sc: %10.4f\n", i, fil_vscAA[0][i]);
    printf("***********************************\n\n");*/

  for(v = 0; v < cm->M; v++)
    {
      ESL_ALLOC(sorted_fil_vscAA[v], sizeof(float) * filN);
      esl_vec_FCopy(fil_vscAA[v], filN, sorted_fil_vscAA[v]); 
      esl_vec_FSortIncreasing(sorted_fil_vscAA[v], filN);
    }
  for(v = 0; v < cm->M; v++) {
    if(cfg->hsi->iscandA[v]) {
      sc = sorted_fil_vscAA[v][Fidx];
      /* E is expected number of hits for db of length 2 * cm->W */
      E  = RJK_ExtremeValueE(sc, cfg->vmuAA[0][v], cfg->vlambdaAA[0][v]);
      /* note partition = 0, this is bogus if more than 1 partition, that's why we die if there are more (see above). */
      fil_calcs  = cfg->hsi->cm_vcalcs[v];
      surv_calcs = E * (cm->W * 2) * cfg->hsi->full_cm_ncalcs;
      printf("SUB %3d sg: %2d sc: %10.4f E: %10.4f filt: %10.4f ", v, cfg->hsi->startA[v], sc, E, fil_calcs);
      fil_calcs += surv_calcs;
      nonfil_calcs = cfg->hsi->full_cm_ncalcs;
      spdup = nonfil_calcs / fil_calcs;
      printf("surv: %10.4f sum : %10.4f full: %10.4f spdup %10.4f\n", surv_calcs, fil_calcs, nonfil_calcs, spdup);
    }  
  }
  
  for(v = 0; v < cm->M; v++) 
    {
      ESL_ALLOC(sorted_fil_EAA[v], sizeof(float *) * filN);
      /* E is expected number of hits for db of length 2 * cm->W */
      /* assumes only 1 partition */
      /*sorted_fil_EAA[v] = RJK_ExtremeValueE(sorted_fil_vscAA[v], cfg->vmuAA[0][v], cfg->vlambdaAA[0][v]);*/
    }      

  /*****************TEMPORARY PRINTF END***************************/

  return eslOK;

 ERROR:
  cm_Fail("calc_best_filter(), memory allocation error.");
  return status; /* NEVERREACHED */
}

#endif 


#if 0
/* Function:  CMHackInsertScores()
 * Incept:    SRE, Wed Jul 24 09:48:22 2002 [St. Louis]
 *
 * Purpose:   Temporary (I hope): make all insert scores 0.
 *            If you let inserts train on the data, you can get
 *            positive insert emission scores. Local alignments,
 *            in particular, can then consist of just a couple of
 *            consensus states and a long string of insert 
 *            states, hitting base-composition-biased sequence
 *            with very high score. This is a Bad Thing.
 *            
 *            The long term solution for this problem will
 *            go in with mixture Dirichlet priors, but for now
 *            (with only Laplace coded), this'll appease the
 *            pitchfork and torches mob at Cambridge.
 *
 * Args:      cm - the model 
 *
 * Returns:   (void)
 *
 * Xref:      STL6 p.93.
 */
void
CMHackInsertScores(CM_t *cm)
{
  int v, x;
  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	for (x = 0; x < cm->abc->K; x++)
	  {
	    cm->esc[v][x]  = 0.;
	    cm->iesc[v][x] = 0;
	  }
    }
  if(cm->ioesc != NULL) { 
    for (v = 0; v < cm->M; v++)
      {
	if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
	  for (x = 0; x < cm->abc->K; x++)
	    cm->ioesc[v][x]  = 0;
	  for(x = cm->abc->K+1; x < cm->abc->Kp-1; x++) /* note boundary conditions, gap, missing data symbols stay IMPOSSIBLE */
	    cm->ioesc[v][x]  = 0;
	}
      }
  }
  if(cm->oesc != NULL) { 
    for (v = 0; v < cm->M; v++)
      {
	if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
	  for (x = 0; x < cm->abc->K; x++)
	    cm->oesc[v][x]  = 0.;
	  for(x = cm->abc->K+1; x < cm->abc->Kp-1; x++) /* note boundary conditions, gap, missing data symbols stay IMPOSSIBLE */
	    cm->oesc[v][x]  = 0.;
	}
      }
  }

  if(cm->cp9 != NULL)
    CP9HackInsertScores(cm->cp9);
}

/* Function:  CP9HackInsertScores()
 * Incept:    EPN, Fri Feb  9 10:59:12 2007
 *
 * Purpose:   Make all inserts 0. Usually called from CMHackInsertScores()
 *            to make the HMM inserts match the CM inserts.
 *
 * Args:      cp9 - the CP9 HMM 
 *
 * Returns:   (void)
 */
void
CP9HackInsertScores(CP9_t *cp9)
{
  int k, x;
  for (k = 0; k <= cp9->M; k++)
    /* CP9 HMMs have insert states in nodes 0 and M */
    for (x = 0; x < MAXDEGEN; x++)
      cp9->isc[x][k] = 0;
}

#endif 



#if 0
/*
 * Function: DuplicateCM
 * Date:     EPN, Thu May 24 09:57:12 2007
 * Purpose:  Given a template CM 'cm', copy it's params into a
 *           new CM which is allocated here and must be
 *           freed by the caller. 
 * 
 * Args:
 *           src          - the template covariance model
 */
CM_t *
DuplicateCM(CM_t *cm)
{
  cm_Fail("Duplicate CM is deprecated, you can undeprecate it, but then you have to figure out how to deal with cm->si SearchInfo_t\n");
       
  int       status;
  int       v;	          /* counter over states */
  int       x;		  /* counter over transitions, residues, nodes */
  CM_t     *new;
  ESL_ALPHABET *abc;
  abc = esl_alphabet_Create(cm->abc->type);
  if(abc == NULL) goto ERROR;

  /* Create the new model and copy everything over except the cp9, stats and ScanMatrix */
  new = CreateCM(cm->nodes, cm->M, cm->abc);
  esl_strdup(cm->name, -1, &(new->name));
  esl_strdup(cm->acc,  -1, &(new->acc));
  esl_strdup(cm->desc, -1, &(new->desc));
  new->flags       = cm->flags;
  new->search_opts = cm->search_opts;
  new->align_opts  = cm->align_opts;
  new->config_opts = cm->config_opts;

  new->nodes      = cm->nodes;
  for(x = 0; x < cm->nodes; x++)
    {
      new->nodemap[x]   = cm->nodemap[x];
      new->ndtype[x]    = cm->ndtype[x];
    }
  for(v = 0; v < cm->M; v++)
    {
      new->sttype[v]  = cm->sttype[v];
      new->ndidx[v]   = cm->ndidx[v];
      new->stid[v]    = cm->stid[v];

      new->cfirst[v]  = cm->cfirst[v];
      new->cnum[v]    = cm->cnum[v];

      new->pnum[v]    = cm->pnum[v];
      new->plast[v]   = cm->plast[v];

      new->begin[v]   = cm->begin[v];
      new->beginsc[v] = cm->beginsc[v];
      new->ibeginsc[v]= cm->ibeginsc[v];
      new->end[v]     = cm->end[v];
      new->endsc[v]   = cm->endsc[v];
      new->iendsc[v]  = cm->iendsc[v];

      /* copy transitions and emissions*/
      for (x = 0; x < MAXCONNECT; x++)
	{
	  new->t[v][x]     = cm->t[v][x];
	  new->tsc[v][x]   = cm->tsc[v][x];
	  new->itsc[v][x]  = cm->itsc[v][x];
	}
      for (x = 0; x < cm->abc->K * cm->abc->K; x++)
	{
	  new->e[v][x]     = cm->e[v][x];
	  new->esc[v][x]   = cm->esc[v][x];
	  new->iesc[v][x]  = cm->iesc[v][x];
	}
    }      
  if(cm->dmin != NULL && cm->dmax != NULL)
    {
      ESL_ALLOC(new->dmin, sizeof(int) * cm->M);
      ESL_ALLOC(new->dmax, sizeof(int) * cm->M);
      for(v = 0; v < cm->M; v++)
	{
	  new->dmin[v] = cm->dmin[v];
	  new->dmax[v] = cm->dmax[v];
	}
    }
  else 
    {
      new->dmin = new->dmax = NULL;
      new->flags &= ~CMH_QDB;
    }
  new->W      = cm->W;
  new->el_selfsc  = cm->el_selfsc;
  new->iel_selfsc = cm->iel_selfsc;
  new->beta  = cm->beta;
  new->tau   = cm->tau;
  new->enf_start = cm->enf_start;
  if(cm->enf_seq != NULL)
    if((status = esl_strdup(cm->enf_seq, -1, &(new->enf_seq))) != eslOK) goto ERROR;
  else new->enf_seq = NULL;
  new->enf_scdiff = cm->enf_scdiff;
  new->ffract     = cm->ffract;
  if(cm->root_trans == NULL)
    new->root_trans = NULL;
  else
    {
      ESL_ALLOC(new->root_trans, sizeof(float) * cm->cnum[0]);
      for (v = 0; v < cm->cnum[0]; v++)
	new->root_trans[v] = cm->root_trans[v];
    }

  new->cp9  = NULL;

  /* calculate the ScanMatrix, if it exists and is valid */
  if(cm->flags & CMH_SCANMATRIX) {
    cm_CreateScanMatrixForCM(new, (cm->smx->flags & cmSMX_HAS_FLOAT), (cm->smx->flags & cmSMX_HAS_INT));
  }

  /* create HMM banded matrix */
  new->hbmx = cm_hb_mx_Create(cm->M);

  /* create CP9 matrix, if it exists in cm */
  if(cm->cp9_mx  != NULL) new->cp9_mx = CreateCP9Matrix(1, cm->clen);
  if(cm->cp9_bmx != NULL) new->cp9_bmx = CreateCP9Matrix(1, cm->clen);

  /* we can't copy the CM stats if they exist */
  if(cm->flags & CMH_GUMBEL_STATS)
    {
      new->stats = AllocCMStats(cm->stats->np);
      CopyCMStats(cm->stats, new->stats);
    }

  /* Copy the CP9 if it exists */
  if(cm->flags & CMH_CP9)
    {
      DuplicateCP9(cm, new);
      new->flags |= CMH_CP9; /* raise the CP9 flag */
    }

  return new;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return NULL; /* never reached */
}
#endif

#if 0
/* Function: cm_AppendComlog()
 * Incept:   SRE, Mon Jan  1 18:23:42 2007 [Casa de Gatos]
 * 
 * Purpose:  Concatenate command line options and append as a new line in the
 *           command line log. Command line log is multiline, with each line
 *           ending in newline char, except for last line.
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
cm_AppendComlog(CM_t *cm, int argc, char **argv)
{
  int   status;
  void *tmp;
  int   n;
  int   i;

  /* figure out length of added command line, and (re)allocate comlog */
  n = argc-1;	/* account for 1 space per arg, except last one */
  for (i = 0; i < argc; i++)
    n += strlen(argv[i]);

  if (cm->comlog != NULL) {
    n += strlen(cm->comlog) + 1; /* +1 for the \n we're going to add to the old comlog */
    ESL_RALLOC(cm->comlog, tmp, sizeof(char)* (n+1));
    strcat(cm->comlog, "\n");
  } else {
    ESL_ALLOC(cm->comlog, sizeof(char)* (n+1));
    *(cm->comlog) = '\0'; /* need this to make strcat work */
  }

  for (i = 0; i < argc-1; i++)
    {
      strcat(cm->comlog, argv[i]);
      strcat(cm->comlog, " ");
    }
  strcat(cm->comlog, argv[argc-1]);
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_SetCdate()
 * Date:     SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 * 
 * Purpose:  Set the <cdate> field in a new CM to the current time.
 *
 *           This function is not reentrant and not threadsafe, because
 *           it calls the nonreentrant ANSI C cdate() function.
 * 
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure. <eslESYS> if the <time()>
 *           system call fails to obtain the calendar time.
 */
int
cm_SetCdate(CM_t *cm)
{
  int    status;
  char  *s = NULL;
  time_t date;

  if ((date   = time(NULL))                       == -1) { status = eslESYS; goto ERROR; }
  if ((status = esl_strdup(cdate(&date), -1, &s)) != eslOK) goto ERROR;
  if ((status = esl_strchop(s, -1))               != eslOK) goto ERROR;
  
  if (cm->cdate != NULL) free(cm->cdate);
  cm->cdate = s;
  return eslOK;

 ERROR:
  if (s != NULL) free(s);
  return status;
}
#endif

/* EPN, Mon Jan  7 14:07:53 2008
 * prior to commit 2291, replaced main entropy weighting functions
 * with HMMER3 entropy weighting functions which are Easelfied.
 */
#if 0

/* Function: CM_Eweight [EPN]
 * based on:
 * Eweight() LSJ 2/6/04
 * 
 * Purpose:  Main entropy-based weighting function. 
 *           
 * Args:  
 *              cm       - the model
 *           **pri       - Model priors.
 *       numb_seqs       - Number of sequences in alignment.
 *       targetent       - Target mean match state entropy. 
 *           
 * Return: eff_no        - New effective sequence number.                         
 */
double
CM_Eweight(CM_t *cm, const Prior_t *pri, float numb_seqs, 
	float targetent)
{
  int status;
  int i;
  int j;
  float eff_no;                  /* New effective sequence number */
  double current;                /* Current mean match state entropy */
  double prevent;                 /* Previous mean match state entropy */
  float scale;                   /* Current model counts scaling factor */
  float leftscale;               /* Bracket scaling value used in binary search. Lowest mean entropy value. */
  float rightscale;              /* Bracket scaling value used in binary search. Highest mean entropy value. */

  double *ent;                    /* Match state entropy values */
  int count;                     /* Counter for binary search */
  int flag;                      /* Used to detect entropy adjustment failure */

  int nmatch_cols;               /* num MATL_nd + MATR_nd + 2 * MATP_nd in CM */
  
  /* analags of parameters from Infernal's prior.c()'s PriorifyCM().*/
  double *counts;                 /* Temp array of match state counts */
  double *probs;                  /* Temp array of match state probs */
  double *mixq;                   /* posterior probs of mixture components, P(q | c) */


  /**************
   * Allocations
   **************/
  ESL_ALLOC(ent,      sizeof(double) * cm->nodes);
  ESL_ALLOC(counts,   sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(probs,    sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(mixq,     sizeof(double) * pri->maxnq);
	  	  
  /*****************
   * Initializations 
   *****************/
  current  = 0.;
  scale    = 1.;
  count    = 0;
  flag     = 0;
  nmatch_cols = 0;

  for(i = 0; i < cm->nodes; i++)
    ent[i] = 0.;

  /***************************************
   * Calculate the starting model entropy 
   ***************************************/

  /* Copy model match state probabilities into our temporary counts[]
   * (Current implementation only considers MATP_MP as a match state,
   *  for MATP nodes, not MATP_ML or MATP_MR (MATL_ML and MATR_MR are
   *  also considered match states)).
   * For nodes i with no match state (BEGL, BEGR, ROOT, BIF and END)
   * ent[i] is left as its initialized value; 0.0. This effectively
   * eliminates any contribution to 'current' from such nodes.
   * Remember our CM is still in counts form, so cm->e[][] is a count
   * not a probability. 
   */
  for(i = 0; i < cm->nodes; i++)
    { 
      if(cm->ndtype[i] == MATP_nd)
	{
	  nmatch_cols += 2; /* two match columns */
	  for(j = 0; j < (MAXABET * MAXABET); j++)
	    counts[j] = cm->e[cm->nodemap[i]][j];
	  /* cm->nodemap[i] = first state, node i (here, MP state) */

	  /* Add priors to the current match state prob dist. (easel/esl_dirichlet.c) */
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						  pri->mbp,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	  /* ent[] is assigned the current MP_st state emission entropy. */
	  ent[i] = esl_vec_DEntropy(probs, (MAXABET * MAXABET));
	}
      else if ((cm->ndtype[i] == MATL_nd) ||
	       (cm->ndtype[i] == MATR_nd))
	{
	  nmatch_cols++;
	  for(j = 0; j < MAXABET; j++)
	    counts[j] = cm->e[cm->nodemap[i]][j];
	  /* cm->nodemap[i] = first state, node i (here, ML or MR state) */

	  /* Add priors to the current match state prob dist. (easel/esl_dirichlet.c) */
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						  pri->mnt,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	  /* ent[] is assigned the current consensus singlet emission entropy. */
	  ent[i] = esl_vec_DEntropy(probs, MAXABET);
	}
      /* other nodes are skipped, ent[i] for these nodes remains 0.0 */
    }
  /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
  current = esl_vec_DSum(ent, cm->nodes)/nmatch_cols;
  /*printf("target ent: %f\n", targetent);*/
  /*printf("0 current: %f\n", current);*/

  /****************************************
   * Initialize binary search bracket values
   *****************************************/

  /* The reason the values seem backwards is because I'm trying to
     bracket my target mean entropy with model count scaling
     factors. A higher scaling factor generally produces a lower
     Entropy and a lower scaling factor produces a higher
     entropy. Thus, the leftscale produces the lowest mean entropy
     bracket and rightscale produces the highest mean entropy
     bracket */
  if(current < targetent){
    leftscale  = 1; 
    rightscale = 0; 
  } 
  else{
    /* Current model has a higher entropy than our target.
       Calculated effective seq numb <= Number of seqs. Design decision.
    */
    printf("[scale=%.2f] [e=%.2f >= %.2f] ...", scale, current, targetent);
    free(mixq);
    free(counts);
    free(probs);
    free(ent);
    return(numb_seqs);
  }
  /***************************************
   * Binary search for target mean entropy
   ***************************************/
  /* Check to see if the current model mean entropy is within 0.01 bits of our target */
  while((current < targetent - 0.01) || (current > targetent + 0.01))
    {
      count++;
      nmatch_cols = 0;
    
    /* Emergency brake in case there is a bug in our binary search.
     * Its more likely that the target entropy is unattainable. */
      if(count > 50){
	printf("\nThe requested target entropy of %f is unattainable. [scale=%.2f] \n", targetent, scale);
	break;
      }
      
      /* Calculate current scaling factor based on bracket values */
      scale = (leftscale + rightscale)/2;
      
      prevent = current;
      
      /*******************************************
       * Scale the counts and re-calc the entropy
       *******************************************/
      /* Re-copy match state probabilities into counts[] */
      for(i = 0; i < cm->nodes; i++)
	{ 
	  if(cm->ndtype[i] == MATP_nd)
	    {
	      nmatch_cols += 2; /* two match columns */
	      for(j = 0; j < (MAXABET * MAXABET); j++)
		counts[j] = cm->e[cm->nodemap[i]][j];
	      /* cm->nodemap[i] = first state, node i (here, MP state) */
	      
	      /* Re-scale the current counts by the previously determined amount. 
	       * (easel/esl_vectorops.c) 
	       */
	      esl_vec_DScale(counts, (MAXABET*MAXABET), scale);
	      
	      /* Re-add priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						      pri->mbp,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	      /* Again, ent[] is assigned the current match emission entropy */
	      ent[i] = esl_vec_DEntropy(probs, (MAXABET * MAXABET));
	    }
	  else if ((cm->ndtype[i] == MATL_nd) ||
		   (cm->ndtype[i] == MATR_nd))
	    {
	      nmatch_cols++;
	      for(j = 0; j < MAXABET; j++)
		counts[j] = cm->e[cm->nodemap[i]][j];
	      /* cm->nodemap[i] = first state, node i (here, ML or MR state) */
	      
	      /* Re-scale the current counts by the previously determined amount. 
	       * (easel/esl_vectorops.c) 
	       */
	      esl_vec_DScale(counts, MAXABET, scale);
	      
	      /* Re-add the priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						      pri->mnt,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	      
	      /* Again, ent[] is assigned the current match emission entropy */
	      ent[i] = esl_vec_DEntropy(probs, MAXABET);
	    }
	  /* other nodes are skipped, ent[i] for these nodes remains 0.0 */
	}
      /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
      current = esl_vec_DSum(ent, cm->nodes)/nmatch_cols;
      /*    printf("current : %f\n", current);*/
      
      /* Adjust the brackets according to the new mean entropy value */
      if(current < targetent){
	leftscale = scale;
      }
      else{
	/* We overshot the target. Replace right bracket with the current scale */
	rightscale = scale;
      }
    }
  free(mixq);
  free(counts);
  free(probs);
  free(ent);
  /**********************************************************************************************
   * End of binary search
   *********************************************************************************************/
  eff_no = numb_seqs * scale;
  /*printf("[scale=%.2f] ", scale);*/
  return(eff_no);

 ERROR: 
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}



/************************************************/
/* Functions just used in debugging/calibrating */ 
/************************************************/

/* Function: ModelContent() LSJ 10/14/03
 * 
 * Purpose:  This is a highly mutable grab-bag function I use  
 *           in benchmarking/debugging to examine model guts.
 *           
 * Args:     
 *           *ent1       - Column entropies for count data.
 *           *ent2       - Column entropies for count+prior data.
 *           M           - number of states in model
 *           
 * Return:   (void)                         
 */
void ModelContent(float *ent1, float *ent2, int M)
{
  int i;
  float sum1, sum2, sum3;
  float mean1, mean2, mean3;

  sum1  = 0;
  sum2  = 0;
  sum3  = 0;
  mean1 = 0;
  mean2 = 0;
  mean3 = 0;

  for(i = 1; i < M+1; i++){
    sum1 += ent1[i];
    sum2 += ent2[i];
    /*    sum3 += relent[i];
     */
    printf("%d\t%2.4f %2.4f %2.4f\n", i, ent1[i], ent2[i], (ent2[i] - ent1[i]));
  }
  mean1 = sum1/M;
  mean2 = sum2/M;
  /*  mean3 = sum3/M;
  fprintf(fp, "Mean Relative Entropy/Column: %2.4f\n", mean3);
  */
  printf("Counts Mean Entropy/Column: %2.4f\n", mean1);
  printf("Counts+Priors Mean Entropy/Column: %2.4f\n", mean2);
  printf("Diff: %2.4f\n", (mean2-mean1));
}


/* Function: CM_Eweight_RE [EPN]
 * based on:
 * Eweight() LSJ 2/6/04
 * 
 * Purpose:  Main entropy-based weighting function. Calculates
 *           relative entropy (RE) instead of entropy. Requires background
 *           distribution. 
 *           
 * Args:  
 *              cm       - the model
 *           **pri       - Model priors.
 *       numb_seqs       - Number of sequences in alignment.
 *     target_relent     - Target mean match state relative entropy. 
 * randomseq[MAXABET]    - null sequence model
 * 
 * Return: eff_no        - New effective sequence number.                         
 */
double
CM_Eweight_RE(CM_t *cm, const Prior_t *pri, float numb_seqs, 
	      float target_relent, float *randomseq)
{
  int status;
  int i;
  int j;
  float eff_no;                  /* New effective sequence number */
  double current;                /* Current mean match state entropy */
  double prevent;                 /* Previous mean match state entropy */
  float scale;                   /* Current model counts scaling factor */
  float leftscale;               /* Bracket scaling value used in binary search. Lowest mean entropy value. */
  float rightscale;              /* Bracket scaling value used in binary search. Highest mean entropy value. */

  double *rel_ent;                    /* Match state relative entropy values */
  int count;                     /* Counter for binary search */
  int flag;                      /* Used to detect entropy adjustment failure */

  int nmatch_cols;               /* num MATL_nd + MATR_nd + 2 * MATP_nd in CM */
  
  /* analags of parameters from Infernal's prior.c()'s PriorifyCM().*/
  double *counts;                 /* Temp array of match state counts */
  double *probs;                  /* Temp array of match state probs */
  double *mixq;                   /* posterior probs of mixture components, P(q | c) */
  double Drandomseq[MAXABET];    /* the randomseq background prob dist, in doubles*/
  double Drandomseq_bp[MAXABET*MAXABET]; /* the randomseq BP background 
					    prob dist, in doubles*/

  /**************
   * Allocations
   **************/
  ESL_ALLOC(rel_ent, sizeof(double) * (cm->nodes));
  ESL_ALLOC(counts,  sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(probs,   sizeof(double) * pri->maxnalpha);
  ESL_ALLOC(mixq,    sizeof(double) * pri->maxnq);
	  	  
  /*****************
   * Initializations 
   *****************/
  current  = 0.;
  scale    = 1.;
  count    = 0;
  flag     = 0;
  nmatch_cols = 0;

  for(i = 0; i < cm->nodes; i++)
    rel_ent[i] = 0.;

  for(i = 0; i < MAXABET; i++)
    Drandomseq[i] = (double) randomseq[i];
  
  for(i = 0; i < MAXABET; i++)
    for(j = 0; j < MAXABET; j++)
      Drandomseq_bp[i*MAXABET+j] = Drandomseq[i] * Drandomseq[j];

  /***************************************
   * Calculate the starting model entropy 
   ***************************************/

  /* Copy model match state probabilities into our temporary counts[]
   * (Current implementation only considers MATP_MP as a match state,
   *  for MATP nodes, not MATP_ML or MATP_MR (MATL_ML and MATR_MR are
   *  also considered match states)).
   * For nodes i with no match state (BEGL, BEGR, ROOT, BIF and END)
   * ent[i] is left as its initialized value; 0.0. This effectively
   * eliminates any contribution to 'current' from such nodes.
   * Remember our CM is still in counts form, so cm->e[][] is a count
   * not a probability. 
   */
  for(i = 0; i < cm->nodes; i++)
    { 
      if(cm->ndtype[i] == MATP_nd)
	{
	  nmatch_cols += 2; /* two match columns */
	  for(j = 0; j < (MAXABET * MAXABET); j++)
	    counts[j] = cm->e[cm->nodemap[i]][j];
	  /* cm->nodemap[i] = first state, node i (here, MP state) */

	  /* Add priors to the current match state prob dist. (easel/esl_dirichlet.c) */
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						  pri->mbp,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	  /* rel_ent[] is assigned the current MP_st state emission relative
	   * entropy. */
	  rel_ent[i] = DRelEntropy(probs, Drandomseq_bp, (MAXABET * MAXABET));
	}
      else if ((cm->ndtype[i] == MATL_nd) ||
	       (cm->ndtype[i] == MATR_nd))
	{
	  nmatch_cols++;
	  for(j = 0; j < MAXABET; j++)
	    counts[j] = cm->e[cm->nodemap[i]][j];
	  /* cm->nodemap[i] = first state, node i (here, ML or MR state) */

	  /* Add priors to the current match state prob dist. (easel/esl_dirichlet.c) */
	  if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						  pri->mnt,
						  mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	  /* rel_ent[] is assigned the current consensus singlet emission 
	     relative entropy. */
	  rel_ent[i] = DRelEntropy(probs, Drandomseq, MAXABET);
	  /*printf("rel_ent[%d] : %f\n", i, rel_ent[i]);*/
	}
      /* other nodes are skipped, ent[i] for these nodes remains 0.0 */
    }
  /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
  current = esl_vec_DSum(rel_ent, cm->nodes)/nmatch_cols;
  printf("target rel ent: %f\n", target_relent);
  printf("0 current: %f\n", current);

  /****************************************
   * Initialize binary search bracket values
   *****************************************/

  /* The reason the values seem backwards is because I'm trying to
     bracket my target mean entropy with model count scaling
     factors. A higher scaling factor generally produces a lower
     Entropy and a lower scaling factor produces a higher
     entropy. Thus, the leftscale produces the lowest mean entropy
     bracket and rightscale produces the highest mean entropy
     bracket */
  if(current > target_relent){
    leftscale  = 1; 
    rightscale = 0; 
  } 
  else{
    /* Current model has a lower relative entropy than our target.
       Calculated effective seq numb <= Number of seqs. Design decision.
    */
    printf("[scale=%.2f] [re=%.2f >= %.2f] ...", scale, current, target_relent);
    free(mixq);
    free(counts);
    free(probs);
    free(rel_ent);
    return(numb_seqs);
  }
  /***************************************
   * Binary search for target mean entropy
   ***************************************/
  /* Check to see if the current model mean entropy is within 0.01 bits of our target */
  ///while((current < target_relent - 0.01) || (current > target_relent + 0.01))
  while((current < target_relent - 0.001) || (current > target_relent + 0.001))
    {
      count++;
      nmatch_cols = 0;
    
    /* Emergency brake in case there is a bug in our binary search.
     * Its more likely that the target entropy is unattainable. */
      if(count > 50){
	printf("\nThe requested target relative entropy of %f is unattainable. [scale=%.2f] \n", target_relent, scale);
	break;
      }
      
      /* Calculate current scaling factor based on bracket values */
      scale = (leftscale + rightscale)/2;
      
      prevent = current;
      
      /*******************************************
       * Scale the counts and re-calc the entropy
       *******************************************/
      /* Re-copy match state probabilities into counts[] */
      for(i = 0; i < cm->nodes; i++)
	{ 
	  if(cm->ndtype[i] == MATP_nd)
	    {
	      nmatch_cols += 2; /* two match columns */
	      for(j = 0; j < (MAXABET * MAXABET); j++)
		counts[j] = cm->e[cm->nodemap[i]][j];
	      /* cm->nodemap[i] = first state, node i (here, MP state) */
	      
	      /* Re-scale the current counts by the previously determined amount. 
	       * (easel/esl_vectorops.c) 
	       */
	      esl_vec_DScale(counts, (MAXABET*MAXABET), scale);
	      
	      /* Re-add priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
						      pri->mbp,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	      /* Again, rel_ent[] is assigned the current match emission 
		 relative entropy */
	      rel_ent[i] = DRelEntropy(probs, Drandomseq_bp, (MAXABET * MAXABET));
	    }
	  else if ((cm->ndtype[i] == MATL_nd) ||
		   (cm->ndtype[i] == MATR_nd))
	    {
	      nmatch_cols++;
	      for(j = 0; j < MAXABET; j++)
		counts[j] = cm->e[cm->nodemap[i]][j];
	      /* cm->nodemap[i] = first state, node i (here, ML or MR state) */
	      
	      /* Re-scale the current counts by the previously determined amount. 
	       * (easel/esl_vectorops.c) 
	       */
	      esl_vec_DScale(counts, MAXABET, scale);
	      
	      /* Re-add the priors to these scaled counts. (easel/esl_dirichlet.c) */
	      if((status = esl_mixdchlet_MPParameters(counts, MAXABET,
						      pri->mnt,
						      mixq, probs)) != eslOK) cm_Fail("esl_mixdchlet_MPParameters() call failed with status: %d\n", status);
	      
	      /* rel_ent[] is assigned the current consensus singlet emission 
		 relative entropy. */
	      rel_ent[i] = DRelEntropy(probs, Drandomseq, MAXABET);
	    }
	  /* other nodes are skipped, ent[i] for these nodes remains 0.0 */
	}
      /* Calculate the mean match state entropy. (easel/esl_vectorops.c::DSum) */
      current = esl_vec_DSum(rel_ent, cm->nodes)/nmatch_cols;
      printf("current : %f\n", current);
      
      /* Adjust the brackets according to the new mean entropy value */
      if(current > target_relent){
	leftscale = scale;
      }
      else{
	/* Replace right bracket with the current scale */
	rightscale = scale;
      }
    }
  free(mixq);
  free(counts);
  free(probs);
  free(rel_ent);
  /**********************************************************************************************
   * End of binary search
   *********************************************************************************************/
  eff_no = numb_seqs * scale;
  printf("[scale=%.2f] ", scale);
  return(eff_no);

 ERROR: 
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* Function:  DRelEntropy()
 *
 * Purpose:   Returns the relative entropy (KL distance) 
 *            between probability vector <p> and <f>
 *            in bits ($\log_2$).
 *
 */
double DRelEntropy(double *p, double *f, int n)
{
  int    i;
  double rel_entropy;
  double eps;
  double temp;

  eps = 0.0000001;

  rel_entropy = 0.;
  for(i = 0; i < n; i++)
    {
      if (f[i] > 0.) temp = f[i]; else temp = f[i] * -1;
      if (temp < (0. + eps)) 
	{ 
	  printf("error in DRelEntropy(), f[%d] is %f\nuh not sure what to do if f[x] is 0! Abort!\n", i, f[i]); 
	  exit(1); 
	}
      if (p[i] > 0.) rel_entropy += p[i] * log(p[i] / f[i]);
    }
  return(1.44269504 * rel_entropy); /* converts to bits */
}
#endif


#if 0
/* Function: CopyCMStats()
 * Incept:   EPN, Tue May 29 06:00:41 2007
 * 
 * Purpose:  Copy the Gumbel and possibly best filter info
 *           in a source CMStats_t object into
 *           a pre-alloc'ed destination CMStats_t object.
 */
int CopyCMStats(CMStats_t *src, CMStats_t *dest)
{
  /* Check contract */
  if(src->np != dest->np)
    cm_Fail("ERROR in CopyCMStats() src->np: %d not equal to alloc'ed dest->np: %d\n", src->np, dest->np);
  
  CopyCMStatsGumbel(src, dest);
  CopyFThrInfo(src->fthrA[FTHR_CM_LC], dest->fthrA[FTHR_CM_LC]);
  CopyFThrInfo(src->fthrA[FTHR_CM_GC], dest->fthrA[FTHR_CM_GC]);
  CopyFThrInfo(src->fthrA[FTHR_CM_LI], dest->fthrA[FTHR_CM_LI]);
  CopyFThrInfo(src->fthrA[FTHR_CM_GI], dest->fthrA[FTHR_CM_GI]);
  return eslOK;
}


/* Function: CopyFThrInfo()
 * Incept:   EPN, Fri May  4 15:54:51 2007
 */
int CopyFThrInfo(CP9FilterThr_t *src, CP9FilterThr_t *dest)
{
  dest->N           = src->N;
  dest->cm_eval     = src->cm_eval;
  dest->l_eval      = src->l_eval;
  dest->l_F         = src->l_F;
  dest->g_eval      = src->g_eval;
  dest->g_F         = src->g_F;
  dest->db_size     = src->db_size;
  dest->was_fast    = src->was_fast;
  return eslOK;
}

/* Function: DuplicateCMStatsGumbel()
 * Incept:   EPN, Mon May  7 06:04:58 2007
 * 
 * Purpose:  Copy the Gumbel stats in a source CMStats_t object into
 *           a pre-alloc'ed destination CMStats_t object.
 */
int DuplicateCMStatsGumbel(CMStats_t *src, CMStats_t *dest)
{
  int i, p;

  /* Check contract */
  if(src->np != dest->np)
    cm_Fail("ERROR in CopyCMStatsGumbel() src->np: %d not equal to alloc'ed dest->np: %d\n", src->np, dest->np);

  for(p = 0; p < src->np; p++) {
    dest->ps[p] = src->ps[p];
    dest->pe[p] = src->pe[p];
  }
  for(i = 0; i < GC_SEGMENTS; i++)
    dest->gc2p[i] = src->gc2p[i]; 

  for(i = 0; i < GUM_NMODES; i++) {
    for(p = 0; p < src->np; p++) {
      dest->gumAA[i][p]->N      = src->gumAA[i][p]->N;
      dest->gumAA[i][p]->L      = src->gumAA[i][p]->L;
      dest->gumAA[i][p]->mu     = src->gumAA[i][p]->mu;
      dest->gumAA[i][p]->lambda = src->gumAA[i][p]->lambda;
    }
  }
  return eslOK;
}

/* Function: debug_print_filterthrinfo
 */
int debug_print_filterthrinfo(CMStats_t *cmstats, CP9FilterThr_t *fthr)
{
  if(! fthr->isvalid) { printf("invalid (not yet set)\n"); return eslOK; }

  double l_x;
  double g_x;
  double tmp_K, tmp_mu;
  tmp_K = exp(cmstats->gumAA[GUM_CP9_GF][0]->mu * cmstats->gumAA[GUM_CP9_GF][0]->lambda) / 
    cmstats->gumAA[GUM_CP9_GF][0]->L;
  tmp_mu = log(tmp_K * ((double) fthr->db_size)) / cmstats->gumAA[GUM_CP9_GF][0]->lambda;
  g_x = tmp_mu - (log(fthr->g_eval) / cmstats->gumAA[GUM_CP9_GF][0]->lambda);

  tmp_K = exp(cmstats->gumAA[GUM_CP9_LF][0]->mu * cmstats->gumAA[GUM_CP9_LF][0]->lambda) / 
    cmstats->gumAA[GUM_CP9_LF][0]->L;
  tmp_mu = log(tmp_K * ((double) fthr->db_size)) / cmstats->gumAA[GUM_CP9_LF][0]->lambda;
  l_x = tmp_mu - (log(fthr->l_eval) / cmstats->gumAA[GUM_CP9_LF][0]->lambda);
  printf("\tN: %d gsc: %.5f F: %.5f (%.5f bits) lsc: %.5f F: %.5f (%.5f bits)\n\tcmsc: %.5f db_size: %d was_fast: %d\n",
	 fthr->N, fthr->g_eval, fthr->g_F, g_x, fthr->l_eval, fthr->l_F, l_x, fthr->cm_eval, fthr->db_size, fthr->was_fast);
  return eslOK;
}

#endif

/* EPN, Tue Jan 15 09:02:07 2008
 * Between revisions 2298-2299:
 * Major overhaul of how cmcalibrate does filter thresholds. 
 * --hybrid option abandoned (temporarily?) 
 * and HMM threshold for ALL relevant CM cutoffs implemented, instead of only doing a single CM cutoff.
 * 
 * The revision 2298 functions are copied below: 
 * These functions are (possibly) also in cmcalibrate-hybrid.c. 

static int process_gumbel_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq, int L,
				   float ***ret_vscAA, float **ret_cp9scA, float **ret_hybscA);
static int process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
				   float ***ret_vscAA, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, float **ret_hyb_scA, int **ret_partA);
static int cm_find_hit_above_cutoff(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq, Parsetree_t *tr, int L, float cutoff, float *ret_sc);
static int predict_cp9_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_vit_cp9scA, float *fil_fwd_cp9scA, int *fil_partA, BestFilterInfo_t *bf);
static int predict_best_sub_cm_roots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **fil_vscAA, int **ret_best_sub_roots);
static int predict_hybrid_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_hybscA, int *fil_partA, GumbelInfo_t **gum_hybA, BestFilterInfo_t *bf, int *ret_getting_faster);
static int  update_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int fthr_mode);
static int cm_fit_histograms(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **vscA, int nscores, int p);

And old mpisupport.c functions used with cmcalibrate --hybrid that were unnecessarily complex.
extern int cmcalibrate_cm_gumbel_results_MPIPackSize(float **vscAA, int nseq, int M, MPI_Comm comm, int *ret_n);
extern int cmcalibrate_cm_gumbel_results_MPIPack(float **vscAA, int nseq, int M, char *buf, int n, int *position, MPI_Comm comm);
extern int cmcalibrate_cm_gumbel_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, int M, float ***ret_vscAA, int *ret_nseq);
extern int cmcalibrate_cp9_filter_results_hyb_MPIPackSize(int nseq, int M, MPI_Comm comm, int *ret_n);;
extern int cmcalibrate_cp9_filter_results_hyb_MPIPack(float **vscAA, float *vit_cp9scA, float *fwd_cp9scA, int *partA, int nseq, int M, char *buf, int n, int *position, MPI_Comm comm);
extern int cmcalibrate_cp9_filter_results_hyb_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, int M, float ***ret_vscAA, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, int **ret_partA, int *ret_nseq);
*/
#if 0


/* Function: process_gumbel_workunit()
 * Date:     EPN, Mon Dec 10 06:09:09 2007
 *
 * Purpose:  A gumbel work unit consists of a CM, and an int specifying a 
 *           number of sequences <nseq>. The job is to randomly generate <nseq> 
 *           sequences using the cm->null background distribution, and 
 *           search them with either (a) the CM, (b) the CM's CP9 HMM, or
 *           (c) a hybrid CM/CP9 CYK/Viterbi scanning algorithm, with hybrid
 *           scanning info in cfg->hsi.
 *
 *           Thus, this function can be run in 1 of 3 modes, determined by the
 *           status of the input variables:
 *         
 *           Mode 1. Gumbel calculation for CM. 
 *           <ret_vscAA> != NULL, <ret_cp9scA> == NULL, <ret_hybscA> == NULL.
 *           Search random sequences with only the CM, either CYK or Inside
 *           (as specified by cm->search_opts>. <ret_vscAA> is filled
 *           with the best CM score at each state for each sequence.
 *
 *           Mode 2. Gumbel calculation for the CP9. 
 *           <ret_vscAA> == NULL, <ret_cp9scA> != NULL, <ret_hybscA> == NULL.
 *           Search random sequences with only the CP9, either Viterbi or Forward
 *           (as specified by cm->search_opts). <ret_cp9scA> is filled
 *           with the best CP9 score for each sequence.
 *
 *           Mode 3. Gumbel calculation for hybrid scanner.
 *           <ret_vscAA> == NULL, <ret_cp9scA> == NULL, <ret_hybscA> != NULL.
 *           Search random sequences with only a hybrid CM/CP9 scanner, 
 *           using hybrid info in cfg->hsi. <ret_hybscA> is filled
 *           with the best hybrid score for each sequence.
 *
 * Args:     go           - getopts
 *           cfg          - cmcalibrate's configuration
 *           errbuf       - for writing out error messages
 *           cm           - the CM (already configured as we want it)
 *           nseq         - number of seqs to generate
 *           L            - length of sequences to search, L==cm->W*2 unless --gumL enabled, in which case
 *                          L = ESL_MAX(cm->W*2, esl_opt_GetInteger(go, "--gumL")
 *           ret_vscAA    - RETURN: [0..v..cm->M-1][0..nseq-1] best score at each state v for each seq
 *           ret_cp9scA   - RETURN: [0..nseq-1] best CP9 score for each seq
 *           ret_hybscA   - RETURN: [0..nseq-1] best hybrid score for each seq
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_gumbel_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			int L, float ***ret_vscAA, float **ret_cp9scA, float **ret_hybscA)
{
  int            status;
  int            mode; /* 1, 2, or 3, determined by status of input args, as explained in 'Purpose' above. */
  float        **vscAA        = NULL;  /* [0..v..cm->M-1][0..i..nseq-1] best CM score for each state, each seq */
  float         *cur_vscA     = NULL;  /* [0..v..cm->M-1]               best CM score for each state cur seq */
  float         *cp9scA       = NULL;  /*                [0..i..nseq-1] best CP9 score for each seq, */
  float         *hybscA       = NULL;  /*                [0..i..nseq-1] best hybrid score for each seq */
  double        *dnull        = NULL; /* double version of cm->null, for generating random seqs */
  int            i;
  int            v;
  ESL_DSQ       *dsq;
  float          sc;
  float          update_i = nseq / 20.;

  /* determine mode, and enforce mode-specific contract */
  if     (ret_vscAA != NULL && ret_cp9scA == NULL && ret_hybscA == NULL) mode = 1; /* calcing CM     gumbel stats */
  else if(ret_vscAA == NULL && ret_cp9scA != NULL && ret_hybscA == NULL) mode = 2; /* calcing CP9    gumbel stats */
  else if(ret_vscAA == NULL && ret_cp9scA == NULL && ret_hybscA != NULL) mode = 3; /* calcing hybrid gumbel stats */
  else ESL_FAIL(eslEINCOMPAT, errbuf, "can't determine mode in process_gumbel_workunit.");

  ESL_DPRINTF1(("in process_gumbel_workunit nseq: %d L: %d mode: %d\n", nseq, L, mode));

  int do_cyk     = FALSE;
  int do_inside  = FALSE;
  int do_viterbi = FALSE;
  int do_forward = FALSE;
  int do_hybrid  = FALSE;
  /* determine algs we'll use and allocate the score arrays we'll pass back */
  if(mode == 1) {
    if(cm->search_opts & CM_SEARCH_INSIDE) do_inside = TRUE;
    else                                   do_cyk    = TRUE;
    ESL_ALLOC(vscAA, sizeof(float *) * cm->M);
    for(v = 0; v < cm->M; v++) ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    ESL_ALLOC(cur_vscA, sizeof(float) * cm->M);
  }
  else if(mode == 2) {
    if(cm->search_opts & CM_SEARCH_HMMVITERBI) do_viterbi = TRUE;
    if(cm->search_opts & CM_SEARCH_HMMFORWARD) do_forward = TRUE;
    if((do_viterbi + do_forward) > 1) ESL_FAIL(eslEINVAL, errbuf, "process_gumbel_workunit, mode 2, and cm->search_opts CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags both raised.");
    ESL_ALLOC(cp9scA, sizeof(float) * nseq); /* will hold Viterbi or Forward scores */
  }
#ifdef HAVE_DEVOPTS
  else if(mode == 3) {
    do_hybrid = TRUE;
    ESL_ALLOC(hybscA,       sizeof(float) * nseq); /* will hold hybrid scores */
  }
#endif
  else if(mode == 3) { /* never entered if HAVE_DEVOPTS is defined */
    ESL_FAIL(eslEINCOMPAT, errbuf, "process_gumbel_workunit(), mode 3 unavailable (HAVE_DEVOPTS is undefined)");
  }

  ESL_DPRINTF1(("do_cyk:     %d\ndo_inside:  %d\ndo_viterbi: %d\ndo_forward: %d\ndo_hybrid: %d", do_cyk, do_inside, do_viterbi, do_forward, do_hybrid)); 
  
  /* fill dnull, a double version of cm->null, but only if we're going to need it to generate random seqs */
  if(cfg->pgc_freq == NULL) {
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);    
  }
  
  /* generate dsqs one at a time and collect best CM scores at each state and/or best overall CP9 score */
  for(i = 0; i < nseq; i++) {
    if(cfg->my_rank == 0 && i > update_i) { /* print status update to stdout */
      printf("=");
      fflush(stdout); 
      update_i += nseq / 20.; 
    }
    if((status = get_random_dsq(cfg, errbuf, cm, dnull, L, &dsq)) != eslOK) return status; 

    /* if nec, search with CM */
    if (do_cyk)    if((status = FastCYKScan    (cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc)) != eslOK) return status;
    if (do_inside) if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, &(cur_vscA), &sc)) != eslOK) return status;
    /* if nec, search with CP9 */
    if (do_viterbi) 
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       NULL,   /* don't want traces back */
			       &(cp9scA[i]))) != eslOK) return status;
    if (do_forward) {
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       &(cp9scA[i]))) != eslOK) return status;
    }
#ifdef HAVE_DEVOPTS
    if (do_hybrid) {
      if((status = cm_cp9_HybridScan(cm, errbuf, cm->cp9_mx, dsq, cfg->hsi, 1, L, cfg->hsi->W, 0., 
				     NULL, /* don't report results */
				     NULL, /* don't want best score at each posn back */
				     NULL, /* don't want the max scoring posn back */
				     &(hybscA[i]))) != eslOK) return status;
    }
#endif

    /*to print seqs to stdout uncomment this block 
    ESL_SQ *tmp;
    tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
    esl_sq_Textize(tmp);
    printf(">seq%d\n%s\n", i, tmp->seq);
    esl_sq_Destroy(tmp);
    */

    free(dsq);
    if (cur_vscA != NULL) /* will be NULL if do_cyk == do_inside == FALSE (mode 2) */
      for(v = 0; v < cm->M; v++) vscAA[v][i] = cur_vscA[v];
    free(cur_vscA);
  }
  if(cfg->my_rank == 0) { printf("=]"); }

  if(dnull != NULL) free(dnull);
  if(ret_vscAA  != NULL) *ret_vscAA  = vscAA;
  if(ret_cp9scA != NULL) *ret_cp9scA = cp9scA;
  if(ret_hybscA != NULL) *ret_hybscA = hybscA;
  return eslOK;

 ERROR:
  return status;
  }

/* Function: process_filter_workunit()
 * Date:     EPN, Mon Dec 10 05:48:35 2007
 *
 * Purpose:  A filter work unit consists of a CM, an int specifying a 
 *           number of sequences <nseq>, and a flag indicating how to search
 *           the sequences. The job is to generate <nseq> sequences from the
 *           CM and search them, the way they're searched is mode dependent
 *           (see below).  with either (a) the CM using bands from
 *           hybrid scanning info in cfg->hsi, then the CP9 HMM with Viterbi and 
 *           Forward or (b) using the hybrid CM/CP9 CYK/Viterbi algorithm
 *           with the hybrid scanning info in cfg->hsi.
 *
 *           Thus, this function can be run in 1 of 3 modes, determined by the
 *           status of the input variables. Note modes 2 and 3 are only possible
 *           if the --hybrid option is enabled, which is only even available if
 *           HAVE_DEVOPTS is defined.
 *         
 *           Mode 1. Scores will be used for calc'ing filter threshold of CP9 HMM.
 *           <ret_vscAA> == NULL, <ret_vit_cp9scA> != NULL, <ret_fwd_cp9scA> != NULL, <ret_hyb_cmscA> == NULL
 *           Emit from CM and score with CP9 Viterbi and Forward, <ret_vit_cp9scA> 
 *           and <ret_fwd_cp9scA> are filled with the best CP9 Viterbi/HMM score 
 *           for each sequence.
 *
 *           Mode 2. Scores will be used for calc'ing filter threshold of CP9 HMM
 *           and CM scores will be used to predict which sub CM roots will be good 
 *           at filtering.
 *           <ret_vscAA> != NULL, <ret_vit_cp9scA> != NULL, <ret_fwd_cp9scA> != NULL, <ret_hyb_cmscA> == NULL
 *           Emit from CM and search first with CM using QDBs from hybrid scanning
 *           info in cfg->hsi, best score from each state of the CM for each 
 *           seq is stored in >ret_vscAA>. Then search (same seq) with CM CP9 Viterbi and Forward, 
 *           <ret_vit_cp9scA> and <ret_fwd_cp9scA> are filled with the best CP9 
 *           Viterbi/HMM score for each sequence.
 *
 *           Mode 3. Scores will be used for calc'ing filter threshold of hybrid scanner.
 *           <ret_vscAA> == NULL, <ret_vit_cp9scA> == NULL, <ret_fwd_cp9scA> == NULL, <ret_hybscA> != NULL
 *           Emit from CM and score with hybrid CM/CP9 CYK/Viterbi scanner, 
 *           <ret_hybscA> are filled with the best hybrid scanner scores 
 *           for each sequence.
 *
 * Args:     go             - getopts
 *           cfg            - cmcalibrate's configuration
 *           errbuf         - for writing out error messages
 *           cm             - the CM (already configured as we want it)
 *           nseq           - number of seqs to generate
 *           ret_vscAA      - RETURN: [0..v..cm->M-1][0..nseq-1] best score at each state v for each seq
 *           ret_vit_cp9scA - RETURN: [0..nseq-1] best Viterbi CP9 score for each seq
 *           ret_fwd_cp9scA - RETURN: [0..nseq-1] best Forward CP9 score for each seq
 *           ret_hybscA     - RETURN: [0..nseq-1] best Hybrid CM/CP9 score for each seq
 *           ret_partA      - RETURN: [0..nseq-1] partition of each seq 
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			float ***ret_vscAA, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, float **ret_hybscA, int **ret_partA)
{
  int            status;
  int            mode; /* 1 or 2 determined by status of input args, as explained in 'Purpose' above. */
  float        **vscAA      = NULL;  /* [0..v..cm->M-1][0..i..nseq-1] best CM score for each state, each seq */
  float         *cur_vscA   = NULL;  /* [0..v..cm->M-1]               best CM score for each state cur seq */
  float         *vit_cp9scA = NULL;  /* [0..i..nseq-1] best CP9 Viterbi score for each seq */
  float         *fwd_cp9scA = NULL;  /* [0..i..nseq-1] best CP9 Forward score for each seq */
  float         *hybscA     = NULL;  /* [0..i..nseq-1] best hybrid CM/CP9 scanner score for each seq */
  int           *partA      = NULL;  /* [0..i..nseq-1] partitions of each seq */
  int            p;                  /* what partition we're in, not used unless emit_from_cm = TRUE */
  int            i, v;
  int            L;
  int            nfailed = 0;
  Parsetree_t   *tr;
  ESL_DSQ       *dsq;
  float          sc;
  int            inside_flag_raised = FALSE;
  float          update_i = nseq / 20.;


  /* determine mode, and enforce mode-specific contract */
  if     (ret_vscAA == NULL && ret_vit_cp9scA != NULL && ret_fwd_cp9scA != NULL && ret_hybscA == NULL) mode = 1; /* running CP9 Viterbi and Forward */
#if HAVE_DEVOPTS
  else if(ret_vscAA != NULL && ret_vit_cp9scA != NULL && ret_fwd_cp9scA != NULL && ret_hybscA == NULL) mode = 2; /* running CM CYK and CP9 Viterbi and Forward */
  else if(ret_vscAA == NULL && ret_vit_cp9scA == NULL && ret_fwd_cp9scA == NULL && ret_hybscA != NULL) mode = 3; /* running hybrid CM/CP9 scanner */
#endif
  else ESL_FAIL(eslEINCOMPAT, errbuf, "can't determine mode in process_filter_workunit.");
  if (ret_partA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "ret_partA is NULL in process_filter_workunit.");

  ESL_DPRINTF1(("in process_filter_workunit nseq: %d mode: %d\n", nseq, mode));
  /* if we get this far, if HAVE_DEVOPTS is undefined, mode MUST be 1 */

#ifndef HAVE_DEVOPTS
  if(mode != 1) ESL_FAIL(eslEINCOMPAT, errbuf, "HAVE_DEVOPTS is undefined, but mode is not 1 in process_filter_workunit(). This shoudln't happen.");
#endif

  /* determine algs we'll use and allocate the score arrays we'll pass back */
  ESL_ALLOC(partA, sizeof(int) * nseq); /* will hold partitions */

  if(mode == 1 || mode == 2) {
    ESL_ALLOC(vit_cp9scA, sizeof(float) * nseq); /* will hold Viterbi scores */
    ESL_ALLOC(fwd_cp9scA, sizeof(float) * nseq); /* will hold Forward scores */
    if(mode == 2) { 
      ESL_ALLOC(vscAA, sizeof(float *) * cm->M);
      for(v = 0; v < cm->M; v++) ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
      ESL_ALLOC(cur_vscA, sizeof(float) * cm->M);
    }
    if(cm->search_opts & CM_SEARCH_INSIDE) { inside_flag_raised = TRUE; cm->search_opts &= ~CM_SEARCH_INSIDE; }
    else inside_flag_raised = FALSE;
  }
  else  /* mode == 3 */
    ESL_ALLOC(hybscA, sizeof(float) * nseq); /* will hold hybrid scores */

  /* generate dsqs one at a time and collect best CM scores at each state and/or best overall CP9 score */
  for(i = 0; i < nseq; i++) {
    if(cfg->my_rank == 0 && i > update_i) { /* print status update to stdout */
      printf("=");
      fflush(stdout); 
      update_i += nseq / 20.; 
    }
    if((status = get_cmemit_dsq(cfg, errbuf, cm, &L, &p, &tr, &dsq)) != eslOK) return status;
    /* we only want to use emitted seqs with a sc > cutoff */
    if((status = cm_find_hit_above_cutoff(go, cfg, errbuf, cm, dsq, tr, L, cfg->cutoffA[p], &sc)) != eslOK) return status;
    while(sc < cfg->cutoffA[p]) { 
      free(dsq); 	
      /* parsetree tr is freed in cm_find_hit_above_cutoff() */
      if((status = get_cmemit_dsq(cfg, errbuf, cm, &L, &p, &tr, &dsq)) != eslOK) return status;
      nfailed++;
      if(nfailed > 1000 * nseq) ESL_FAIL(eslERANGE, errbuf, "process_filter_workunit(), max number of failures (%d) reached while trying to emit %d seqs.\n", nfailed, nseq);
      if((status = cm_find_hit_above_cutoff(go, cfg, errbuf, cm, dsq, tr, L, cfg->cutoffA[p], &sc)) != eslOK) return status;
    }

    /*to print seqs to stdout uncomment this block  
    ESL_SQ *tmp;
    tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
    esl_sq_Textize(tmp);
    printf(">seq%d\n%s\n", i, tmp->seq);
    esl_sq_Destroy(tmp);
    */

    partA[i] = p;
    assert(partA[i] < cfg->np);
    ESL_DPRINTF1(("i: %d nfailed: %d cutoff: %.3f p: %d\n", i, nfailed, cfg->cutoffA[p], p));

    /* search dsq with mode-specific search algs */
    if(mode == 1 || mode == 2) {
      /* Note: in mode 2, with FastCYKScan, we use cfg->hsi->smx scan matrix, which may have qdbs calc'ed differently than cm->smx */
      if(mode == 2) { if((status = FastCYKScan(cm, errbuf, cfg->hsi->smx, dsq, 1, L, 0., NULL, &(cur_vscA), NULL)) != eslOK) return status; }
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       NULL,   /* don't want traces back */
			       &(vit_cp9scA[i]))) != eslOK) return status;
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			       TRUE,   /* yes, we are scanning */
			       FALSE,  /* no, we are not aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL,   /* don't want best score at each posn back */
			       NULL,   /* don't want the max scoring posn back */
			       &(fwd_cp9scA[i]))) != eslOK) return status;
    }
    else { /* mode == 3 */
      if((status = cm_cp9_HybridScan(cm, errbuf, cm->cp9_mx, dsq, cfg->hsi, 1, L, cfg->hsi->W, 0., 
				     NULL, /* don't report results */
				     NULL, /* don't want best score at each posn back */
				     NULL, /* don't want the max scoring posn back */
				     &(hybscA[i]))) != eslOK) return status;
    }
    free(dsq);
    if (cur_vscA != NULL) /* will be NULL if do_cyk == do_inside == FALSE (mode 3) */
      for(v = 0; v < cm->M; v++) vscAA[v][i] = cur_vscA[v];
    free(cur_vscA);
  }
  if(cfg->my_rank == 0) { printf("=]"); }
  *ret_partA = partA;
  if(ret_vscAA      != NULL)  *ret_vscAA      = vscAA;
  if(ret_vit_cp9scA != NULL)  *ret_vit_cp9scA = vit_cp9scA;
  if(ret_fwd_cp9scA != NULL)  *ret_fwd_cp9scA = fwd_cp9scA;
  if(ret_hybscA != NULL)      *ret_hybscA     = hybscA;

  if(inside_flag_raised) cm->search_opts |= CM_SEARCH_INSIDE;

  return eslOK;

 ERROR:
  return status;
}


/*
 * Function: cm_find_hit_above_cutoff()
 * Date:     EPN, Wed Sep 12 04:59:08 2007
 *
 * Purpose:  Given a CM, a sequence, and a cutoff, try to 
 *           *quickly* answer the question: Does this sequence 
 *           contain a hit to the CM above the cutoff?
 *           To do this we first check the parsetree score, and
 *           then do do up to 3 iterations of search.
 *           The first 2 are performend with j and d bands 
 *           (of decreasing tightness), then default 
 *           search (with QDB unless --noqdb enabled) is done.
 *           We return TRUE if any search finds a hit above
 *           cutoff, and FALSE otherwise.
 *
 * Args:     go              - getopts
 *           cfg             - cmcalibrate's configuration
 *           errbuf          - char buffer for error message
 *           cm              - CM to emit from
 *           dsq             - the digitized sequence to search
 *           tr              - parsetree for dsq
 *           L               - length of sequence
 *           cutoff          - bit score cutoff 
 *           ret_sc          - score of a hit within dsq, if < cutoff,
 *                             this is score of best hit within dsq, which
 *                             means no hit with sc > cutoff exists. If > cutoff,
 *                             not necessarily score of best hit within dsq.
 *
 * Returns:  eslOK on success. other status code upon failure, errbuf filled with error message.
 */
int 
cm_find_hit_above_cutoff(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq,
			 Parsetree_t *tr, int L, float cutoff, float *ret_sc)
{
  int status;
  int turn_qdb_back_on = FALSE;
  int turn_hbanded_back_off = FALSE;
  int turn_hmmscanbands_back_off = FALSE;
  double orig_tau = cm->tau;
  float sc;
  float size_limit = esl_opt_GetReal(go, "--mxsize");

#if eslDEBUGLEVEL >= 1
  int init_flags       = cm->flags;
  int init_search_opts = cm->search_opts;
#endif

  if(ret_sc == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_find_hit_above_cutoff(), ret_sc == NULL.\n");

  /* Determine if this sequence has a hit in it above the cutoff as quickly as possible. 
   * Stage 0: Check parsetree score
   * Stage 1: HMM banded search tau = 1e-2
   * Stage 2: HMM banded search with scanning bands, tau = 1e-10
   * Stage 3: QDB search (CYK or inside), beta = --beta, (THIS IS MOST LENIENT SEARCH WE'LL DO)
   *
   * The earliest stage at which we find a hit > cutoff at any stage, we return cm->flags, cm->search_opts
   * to how they were when we entered, and return TRUE.
   *
   * NOTE: We don't do a full non-banded parse to be 100% sure we don't exceed the cutoff, 
   * unless --noqdb was enabled (ScanMatrix_t *smx stores dn/dx (min/max d) for each state), 
   * because we assume the --beta value used in *this* cmcalibrate 
   * run will also be used for any cmsearch runs.
   */

  sc = ParsetreeScore(cm, tr, dsq, FALSE); 
  FreeParsetree(tr);
  if(sc > cutoff || L == 0) { /* parse score exceeds cutoff, or zero length sequence (only 1 path is possible, must be parse score) */
    ESL_DASSERT1((cm->flags       == init_flags));
    ESL_DASSERT1((cm->search_opts == init_search_opts));
    /* printf("0 sc: %10.4f\n", sc); */
    *ret_sc = sc;
    return eslOK;
  } 

  if(!(cm->search_opts & CM_SEARCH_NOQDB))        turn_qdb_back_on = TRUE;
  if(!(cm->search_opts & CM_SEARCH_HBANDED))      turn_hbanded_back_off = TRUE;
  if(!(cm->search_opts & CM_SEARCH_HMMSCANBANDS)) turn_hmmscanbands_back_off = TRUE;

  cm->search_opts |= CM_SEARCH_NOQDB;

  /* stage 1 */
  cm->search_opts |= CM_SEARCH_HBANDED;
  cm->tau = 0.01;
  if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, TRUE, 0)) != eslOK) return status;
  status = FastCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, cm->hbmx, size_limit, &sc);
  if(status == eslOK) { /* FastCYKScanHB() successfully finished */
    if(sc > cutoff) { /* score exceeds cutoff, we're done, reset search_opts and return */
      if(turn_qdb_back_on)        cm->search_opts &= ~CM_SEARCH_NOQDB; 
      if(turn_hbanded_back_off) { cm->search_opts &= ~CM_SEARCH_HBANDED; cm->tau = orig_tau; }
      ESL_DASSERT1((cm->flags       == init_flags));
      ESL_DASSERT1((cm->search_opts == init_search_opts));
      *ret_sc = sc;
      return eslOK;
    }
  }
  else if (status != eslERANGE) return status; /* else if status == eslERANGE, FastCYKScanHB() couldn't grow its DP matrix big enough, move onto next stage */

  /* stage 2 */
  cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
  cm->tau = 1e-10;
  if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, TRUE, 0)) != eslOK) return status;
  status = FastCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, cm->hbmx, size_limit, &sc);
  if(status == eslOK) { /* FastCYKScanHB() successfully finished */
    if(sc > cutoff) { /* score exceeds cutoff, we're done, reset search_opts and return */
      if(turn_qdb_back_on)             cm->search_opts &= ~CM_SEARCH_NOQDB; 
      if(turn_hbanded_back_off)      { cm->search_opts &= ~CM_SEARCH_HBANDED;      cm->tau = orig_tau; }
      if(turn_hmmscanbands_back_off) { cm->search_opts &= ~CM_SEARCH_HMMSCANBANDS; cm->tau = orig_tau; }
      ESL_DASSERT1((cm->flags       == init_flags));
      ESL_DASSERT1((cm->search_opts == init_search_opts));
      *ret_sc = sc;
      return eslOK;
    }
  }
  else if (status != eslERANGE) return status; /* else if status == eslERANGE, FastCYKScanHB() couldn't grow its DP matrix big enough, move onto next stage */

  /* stage 3, use 'default' dmin, dmax (which could be NULL) CYK or Inside */
  cm->search_opts &= ~CM_SEARCH_HBANDED;
  cm->search_opts &= ~CM_SEARCH_HMMSCANBANDS;
  if(turn_qdb_back_on) cm->search_opts &= ~CM_SEARCH_NOQDB; 

  if(cm->search_opts & CM_SEARCH_INSIDE) {
    if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &sc)) != eslOK) return status;
  }
  else { 
    if((status = FastCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, NULL, &sc)) != eslOK) return status;
  }
  if(!turn_hbanded_back_off)      { cm->search_opts |= CM_SEARCH_HBANDED;      cm->tau = orig_tau; }
  if(!turn_hmmscanbands_back_off) { cm->search_opts |= CM_SEARCH_HMMSCANBANDS; cm->tau = orig_tau; }
  ESL_DASSERT1((cm->flags       == init_flags));
  ESL_DASSERT1((cm->search_opts == init_search_opts));

  /*if(sc > cutoff) { printf("3 sc: %10.4f\n", sc); }*/
  *ret_sc = sc;
  return eslOK;
}


/* Function: predict_cp9_filter_speedup()
 * Date:     EPN, Mon Dec 10 11:55:24 2007
 *
 * Purpose:  Given a CM and scores for a CP9 Viterbi and Forward scan
 *           of target seqs predict the speedup with an HMM filter, Forward and
 *           Viterbi, then update a BestFilterInfo_t object to 
 *           hold info on the faster of the two.
 *            
 * Args:     go  - command line options
 *           cfg - cmcalibrate's cfg object, mucho data (probably too much)
 *           errbuf - for printing error messages
 *           cm - the model
 *           fil_vit_cp9scA - [0..i..filN-1] best Viterbi score in sequence i 
 *           fil_fwd_cp9scA - [0..i..filN-1] best Foward score in sequence i 
 *           fil_partA      - [0..i..filN-1] partition of sequence i 
 *           bf             - BestFilterInfo_t object, we'll update this to hold info on a Viterbi or Forward filter strategy 
 *
 * Returns:  Updates BestFilterInfo_t object <bf> to hold info on fastest HMM filter, Viterbi or Forward
 *           eslOK on success;
 *           Other easel status code on an error with errbuf filled with error message.
 */
int
predict_cp9_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_vit_cp9scA, float *fil_fwd_cp9scA, int *fil_partA, BestFilterInfo_t *bf)
{
  int    status;
  float *sorted_fil_vit_EA;       /* sorted Viterbi E-values, so we can easily choose a threshold */
  float *sorted_fil_fwd_EA;       /* sorted Forward E-values, so we can easily choose a threshold */
  float  vit_E, fwd_E;            /* a Viterbi and Forward E value */
  int    cp9_vit_mode, cp9_fwd_mode; /* a Viterbi, Forward Gumbel mode, respectively  */
  int    evalue_L;                /* database length used for calcing E-values in CP9 gumbels from cfg->cmstats */
  float  fil_calcs;               /* number of million dp calcs predicted for the HMM filter scan */
  float  vit_surv_calcs;          /* number of million dp calcs predicted for the CM scan of Viterbi filter survivors */
  float  fwd_surv_calcs;          /* number of million dp calcs predicted for the CM scan of Forward filter survivors */
  float  vit_fil_plus_surv_calcs; /* Viterbi filter calcs plus survival calcs */ 
  float  fwd_fil_plus_surv_calcs; /* Foward filter calcs plus survival calcs */ 
  float  nonfil_calcs;            /* number of million dp calcs predicted for a full CM scan */
  float  vit_spdup, fwd_spdup;    /* predicted speedups for Viterbi and Forward */
  int    i, p;                    /* counters */
  int    cmi = cfg->ncm-1;        /* CM index we're on */
  int    Fidx;                    /* index in sorted scores that threshold will be set at (1-F) * N */
  float  filE = esl_opt_GetReal(go, "--filE"); /* E-value cutoff for accepting CM hits for filter test */
  float  F = esl_opt_GetReal(go, "--F"); /* fraction of CM seqs we require filter to let pass */
  int    filN  = esl_opt_GetInteger(go, "--filN"); /* number of sequences we emitted from CM for filter test */
  float  Starg   = esl_opt_GetReal(go, "--starg"); /* target survival fraction */
  float  E_Starg;                 /* E-value threshold that exactly satisifies Starg survival fraction */
  float  vit_E_F1;                /* viterbi E value if F were == 1.0 */
  float  fwd_E_F1;                /* forward E value if F were == 1.0 */
  float  surv_res_per_hit;        /* expected number of residues to survive filter from DB for each hit 2*W-avglen[0], twice W minus the average lenght of a hit from QDB calc */
  float  Smin;                    /* minimally useful survival fraction, any less than this and our filter would be (predicted to be) doing more than 10X the work of the CM */
  float  E_min;                   /* E-value threshold that exactly satisifies Smin survival fraction */
  float  vit_surv_fract;          /* predicted survival fraction for Viterbi filter */
  float  fwd_surv_fract;          /* predicted survival fraction for Forward filter */
  
  cp9_vit_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? GUM_CP9_LV : GUM_CP9_GV;
  cp9_fwd_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? GUM_CP9_LF : GUM_CP9_GF;

  /* contract checks */
  if(! (cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), gumbel stats for CP9 viterbi mode: %d are not valid.\n", cp9_vit_mode);
  if(! (cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), gumbel stats for CP9 forward mode: %d are not valid.\n", cp9_fwd_mode);
  if(cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->L != cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), db length for gumbel stats for CP9 viterbi (%d) and forward (%d) differ.\n", cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->L, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][0]->L);

  evalue_L = cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][0]->L;

  /* contract checks specific to case when there is more than 1 partition */
  if(cfg->cmstatsA[cfg->ncm-1]->np != 1) { 
    for(p = 1; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) {
      if(evalue_L != cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), partition %d db length (%d) for Viterbi gumbel stats differ than from partition 1 Viterbi db length (%d).\n", p, cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->L, evalue_L);
      if(evalue_L != cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_cp9_filter_speedup(), partition %d db length (%d) for Forward gumbel stats differ than from partition 1 Viterbi db length (%d).\n", p, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->L, evalue_L);
    }
  }

  Fidx  = (int) ((1. - F) * (float) filN);

  /* convert bit scores to E-values and sort them */
  ESL_ALLOC(sorted_fil_vit_EA, sizeof(float) * filN);
  ESL_ALLOC(sorted_fil_fwd_EA, sizeof(float) * filN);
  for(i = 0; i < filN; i++) { 
    p = fil_partA[i];
    sorted_fil_vit_EA[i] = RJK_ExtremeValueE(fil_vit_cp9scA[i], cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->mu, cfg->cmstatsA[cmi]->gumAA[cp9_vit_mode][p]->lambda); 
    sorted_fil_fwd_EA[i] = RJK_ExtremeValueE(fil_fwd_cp9scA[i], cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->mu, cfg->cmstatsA[cmi]->gumAA[cp9_fwd_mode][p]->lambda); 
  }
  esl_vec_FSortDecreasing(sorted_fil_vit_EA, filN);
  vit_E = sorted_fil_vit_EA[Fidx];
  esl_vec_FSortDecreasing(sorted_fil_fwd_EA, filN);
  fwd_E = sorted_fil_fwd_EA[Fidx];

  /* now vit_E and fwd_E are expected number of CP9 Viterbi/Forward hits with score above threshold 
   * (Fidx'th best score) in a sequence DB of length evalue_L, convert that DB size to cfg->dbsize */
  vit_E *= cfg->dbsize / evalue_L;
  fwd_E *= cfg->dbsize / evalue_L;

  fil_calcs  = cfg->full_cp9_ncalcs;  /* fil_calcs is millions of DP calcs for CP9 scan of 1 residue */
  fil_calcs *= cfg->dbsize;           /* fil_calcs is millions of DP calcs for CP9 scan of length cfg->dbsize */
  nonfil_calcs = cfg->full_vcalcs[0]; /* total number of millions of DP calculations for full CM scan of 1 residue */
  nonfil_calcs *= cfg->dbsize;        /* now nonfil-calcs corresponds to cfg->dbsize */

  /* determine our thresholds:
   * 1. if vit_E and/or fwd_E yield predicted survival fractions that are less than Starg:
   *    rewrite vit_E or fwd_E as the maximum E-value threshold that sees ALL hits (when F==1.0,
   *    this is {vit,fwd}_E_F1 below), and the E-value threshold that exactly satisfies Starg (E_Starg below).
   * 2. if after 1, vit_E and/or fwd_E still yield predicted survival fractions less than Smin,
   *    the survival fraction at which the number of DP calcs for the survivors of the filter is only
   *    10% the number of dp calcs for the filter, then we set vit_E or fwd_E to the E value that
   *    exactly satisfies Smin.
   */
  surv_res_per_hit = (2. * cm->W - (cfg->avglen[0])); /* avg length of surviving fraction of db from a single hit (cfg->avglen[0] is avg subseq len in subtree rooted at v==0, from QDB calculation) */
  vit_E_F1= (sorted_fil_vit_EA[0] * (cfg->dbsize / evalue_L));
  fwd_E_F1= (sorted_fil_fwd_EA[0] * (cfg->dbsize / evalue_L));
  E_Starg = (Starg * (float) cfg->dbsize) / surv_res_per_hit;
  Smin    = fil_calcs / (10. * nonfil_calcs);
  E_min   = (Smin * (float) cfg->dbsize) / surv_res_per_hit;
  if(Starg < Smin) { /* we never go less than Smin */
    Starg = Smin;
    E_Starg = E_min;
  }
  
  vit_surv_fract = (vit_E * surv_res_per_hit) / (float) cfg->dbsize;
  fwd_surv_fract = (fwd_E * surv_res_per_hit) / (float) cfg->dbsize;
  
  ESL_DPRINTF1(("vit_E:    %.5f\n", vit_E));
  ESL_DPRINTF1(("vit_surv: %.10f\n", vit_surv_fract));
  ESL_DPRINTF1(("fwd_E:    %.5f\n", fwd_E));
  ESL_DPRINTF1(("fwd_surv: %.10f\n", fwd_surv_fract));

  ESL_DPRINTF1(("vit_E_F1: %.5f\n", vit_E_F1));
  ESL_DPRINTF1(("fwd_E_F1: %.5f\n", fwd_E_F1));
  ESL_DPRINTF1(("Starg:    %.5f\n", Starg));
  ESL_DPRINTF1(("E_Starg:  %.5f\n", E_Starg));
  ESL_DPRINTF1(("Smin:     %.5f\n", Smin));
  ESL_DPRINTF1(("E_min:    %.5f\n", E_min));

  if(vit_surv_fract < Starg) { 
    if(vit_E_F1 < E_Starg) { 
      vit_E = vit_E_F1;
      if(cfg->be_verbose) printf("set vit_E as vit_E_F1: %.5f\n", vit_E);
      else { ESL_DPRINTF1(("set vit_E as vit_E_F1: %.5f\n", vit_E)); }
    }
    else { 
      vit_E = E_Starg;
      if(cfg->be_verbose) printf("set vit_E as E_Starg: %.5f\n", vit_E);
      else { ESL_DPRINTF1(("set vit_E as E_Starg: %.5f\n", vit_E)); }
    }
    if(vit_E < E_min) {
      vit_E = E_min;
      if(cfg->be_verbose) printf("set vit_E as E_min: %.5f\n", vit_E);
      else { ESL_DPRINTF1(("set vit_E as E_min: %.5f\n", vit_E)); }
    }      
  }

  if(fwd_surv_fract < Starg) { 
    if(fwd_E_F1 < E_Starg) { 
      fwd_E = fwd_E_F1;
      if(cfg->be_verbose) printf("set fwd_E as fwd_E_F1: %.5f\n", fwd_E);
      else { ESL_DPRINTF1(("set fwd_E as fwd_E_F1: %.5f\n", fwd_E)); }
    }
    else { 
      fwd_E = E_Starg;
      if(cfg->be_verbose) printf("set fwd_E as E_Starg: %.5f\n", fwd_E);
      else { ESL_DPRINTF1(("set fwd_E as E_Starg: %.5f\n", fwd_E)); }
    }
    if(fwd_E < E_min) {
      fwd_E = E_min;
      if(cfg->be_verbose) printf("set fwd_E as E_min: %.5f\n", fwd_E);
      else { ESL_DPRINTF1(("set fwd_E as E_min: %.5f\n", fwd_E)); }
    }      
  }

  for(i = 0; i < filN; i++) ESL_DPRINTF1(("HMM i: %4d vit E: %15.10f fwd E: %15.10f\n", i, sorted_fil_vit_EA[i], sorted_fil_fwd_EA[i]));

  /* calculate speedup for Viterbi */
  vit_surv_calcs = vit_E *     /* number of hits expected to survive filter */
    surv_res_per_hit *         /* avg length of surviving fraction of db from a single hit */
    cfg->full_vcalcs[0];       /* number of calculations for full CM scan of 1 residue */

  fwd_surv_calcs = fwd_E *     /* number of hits expected to survive filter */
    surv_res_per_hit *         /* avg length of surviving fraction of db from a single hit */
    cfg->full_vcalcs[0];       /* number of calculations for full CM scan of 1 residue */

  vit_fil_plus_surv_calcs = fil_calcs + vit_surv_calcs; /* total number of millions of DP calculations expected using the CP9 viterbi filter for scan of cfg->dbsize */
  vit_spdup = nonfil_calcs / vit_fil_plus_surv_calcs;
  fwd_fil_plus_surv_calcs = (fil_calcs * 2.) + fwd_surv_calcs; /* total number of millions of DP calculations expected using the CP9 forward filter for scan of cfg->dbsize (logsum corrected, Forward calcs *= 2.) */
  fwd_spdup = nonfil_calcs / fwd_fil_plus_surv_calcs;
  /* We multiply number of forward calculations by 2.0 to correct for the fact that Forward takes about 2X as long as Viterbi, b/c it requires logsum operations instead of ESL_MAX's,
   * so we factor this in when calc'ing the predicted speedup. */

  if(cfg->be_verbose) { 
    printf("\nHMM(vit) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", vit_E, fil_calcs, vit_surv_calcs, vit_fil_plus_surv_calcs, nonfil_calcs, vit_spdup);
    printf("HMM(fwd) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", fwd_E, fil_calcs, fwd_surv_calcs, fwd_fil_plus_surv_calcs, nonfil_calcs, fwd_spdup);
  }
  else {
    ESL_DPRINTF1(("\nHMM(vit) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", vit_E, fil_calcs, vit_surv_calcs, vit_fil_plus_surv_calcs, nonfil_calcs, vit_spdup));
    ESL_DPRINTF1(("HMM(fwd) E: %15.10f filt: %10.4f surv: %10.4f logsum corrected sum: %10.4f full CM: %10.4f spdup %10.4f\n", fwd_E, fil_calcs, fwd_surv_calcs, fwd_fil_plus_surv_calcs, nonfil_calcs, fwd_spdup));
  }

  if(esl_opt_GetBoolean(go, "--fviterbi")) { /* user specified Viterbi */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_VITERBI, vit_E, fil_calcs, vit_fil_plus_surv_calcs)) != eslOK) return status;
  }
  else if(esl_opt_GetBoolean(go, "--fforward")) { /* user specified Forward */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_FORWARD, fwd_E, fil_calcs, fwd_fil_plus_surv_calcs)) != eslOK) return status;
  }
  else if (vit_spdup > fwd_spdup) { /* Viterbi is winner */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_VITERBI, vit_E, fil_calcs, vit_fil_plus_surv_calcs)) != eslOK) return status;
  }
  else { /* Forward is winner */
    if((status = SetBestFilterInfoHMM(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, FILTER_WITH_HMM_FORWARD, fwd_E, fil_calcs, fwd_fil_plus_surv_calcs)) != eslOK) return status;
  }
  return eslOK;

 ERROR:
  return status; 
}


/* Function: predict_best_sub_cm_roots()
 * Date:     EPN, Mon Dec 10 15:56:00 2007
 *
 * Purpose:  Given a CM and scores for a CM scan of target seqs
 *           predict the best sub CM roots we could use to 
 *           filter with.
 *            
 * Returns:  eslOK on success;
 *           Other status code on error, with error message in errbuf.
 */
int 
predict_best_sub_cm_roots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float **fil_vscAA, int **ret_sorted_best_roots_v)
{
  int    status;
  float **sorted_fil_vscAA;       /* [0..v..cm->M-1][0..filN-1] best score for each state v, each target seq */
  float  fil_calcs;               /* number of million dp calcs predicted for the HMM filter scan */
  float  surv_calcs;              /* number of million dp calcs predicted for the CM scan of filter survivors */
  float  fil_plus_surv_calcs;     /* filter calcs plus survival calcs */ 
  float  nonfil_calcs;            /* number of million dp calcs predicted for a full CM scan */
  float  spdup;                   /* predicted speedup a sub CM filter */
  float  E, tmp_E;                /* E value */
  float  sc;                      /* bit score */
  int    i, p, s, v;              /* counters */
  int    Fidx;                    /* index in sorted scores that threshold will be set at (1-F) * N */
  float  F = esl_opt_GetReal(go, "--F"); /* fraction of CM seqs we require filter to let pass */
  int    filN  = esl_opt_GetInteger(go, "--filN"); /* number of sequences we emitted from CM for filter test */
  int    nstarts;                  /* # start states (and start groups) in the CM, from cfg->hsi */                                 
  int   *best_per_start_v;         /* sub CM filter state v that gives best speedup per start group */
  float *best_per_start_spdup;     /* best sub CM filter state speedup per start group */

  if(ret_sorted_best_roots_v == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_best_sub_cm_roots, ret_sorted_best_roots_v == NULL.\n");

  Fidx  = (int) ((1. - F) * (float) filN);

  ESL_ALLOC(sorted_fil_vscAA, sizeof(float *) * cm->M);
  /*ESL_ALLOC(sorted_fil_EAA, sizeof(float *) * cm->M);*/

  nstarts = cfg->hsi->nstarts;
  ESL_ALLOC(best_per_start_v,     sizeof(int)   * nstarts);
  ESL_ALLOC(best_per_start_spdup, sizeof(float) * nstarts);
  for(s = 0; s < nstarts; s++) {
    best_per_start_v[s] = -1;
    best_per_start_spdup[s] = -eslINFINITY;
  }

  for(v = 0; v < cm->M; v++) {
    ESL_ALLOC(sorted_fil_vscAA[v], sizeof(float) * filN);
    esl_vec_FCopy(fil_vscAA[v], filN, sorted_fil_vscAA[v]); 
    esl_vec_FSortIncreasing(sorted_fil_vscAA[v], filN);
  }
  ESL_DPRINTF1(("vscAA[0] scores:\n"));
  for(i = 0; i < filN; i++) ESL_DPRINTF1(("i: %4d sc: %10.4f\n", i, sorted_fil_vscAA[0][i]));

  for(v = 0; v < cm->M; v++) {
    if(cfg->hsi->iscandA[v]) {
      sc = sorted_fil_vscAA[v][Fidx];
      /* set E as E-value for sc from partition that gives sc the lowest E-value (conservative, sc will be at least as significant as E across all partitions) */
      E  = RJK_ExtremeValueE(sc, cfg->vmuAA[0][v], cfg->vlambdaAA[0][v]); 
      for(p = 1; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) {
	tmp_E = RJK_ExtremeValueE(sc, cfg->vmuAA[p][v], cfg->vlambdaAA[p][v]); 
	if(tmp_E < E) E = tmp_E;
      }
      /* E is now expected number of hits for db of cfg->length 2 * cm->W */
      E *= cfg->dbsize / (cm->W * 2.);
      /* E is now expected number of hits for db of cfg->dbsize */

      fil_calcs   = cfg->hsi->cm_vcalcs[v]; /* fil_calcs is millions of DP calcs for sub CM (root = v) scan of 1 residue */
      fil_calcs  *= cfg->dbsize;            /* fil_calcs is millions of DP calcs for sub CM (root = v) scan of length cfg->dbsize */
      surv_calcs = E *     /* number of hits expected to survive filter */
	(2. * cm->W - (cfg->avglen[v])) * /* average length of surviving fraction of db from a single hit (cfg->avglen[v] is avg subseq len in  subtree rooted at v */
	cfg->hsi->full_cm_ncalcs; /* number of calculations for full CM scan of 1 residue */
      fil_plus_surv_calcs = fil_calcs + surv_calcs;
      nonfil_calcs = cfg->hsi->full_cm_ncalcs;      /* total number of millions of DP calculations for full CM scan of 1 residue */
      nonfil_calcs *= cfg->dbsize;                  /* now nonfil-calcs corresponds to cfg->dbsize */
      spdup = nonfil_calcs / fil_plus_surv_calcs;
      if(cfg->be_verbose) { printf("SUB %3d sg: %2d sc: %10.4f E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", v, cfg->hsi->startA[v], sc, E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup); } 
      else                { ESL_DPRINTF1(("SUB %3d sg: %2d sc: %10.4f E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", v, cfg->hsi->startA[v], sc, E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup)); } 
      s = cfg->hsi->startA[v];
      if(spdup > best_per_start_spdup[s] && cm->ndidx[v] != 0) { /* can't filter with a state in node 0 */
	best_per_start_v[s]     = v;
	best_per_start_spdup[s] = spdup;
      }	
    }  
  }
  for(s = 0; s < nstarts; s++) { 
    if(cfg->be_verbose) printf("START %d v: %d spdup: %10.4f\n", s, best_per_start_v[s], best_per_start_spdup[s]);
    else {       ESL_DPRINTF1(("START %d v: %d spdup: %10.4f\n", s, best_per_start_v[s], best_per_start_spdup[s])); }
    ESL_DASSERT1((best_per_start_v[s] != -1));
  }

  /* sort the best sub CM roots (1 per start group) by their speedup,
   * this is an embarassing N^2 sorting, but biggest RNAs have ~ 100 starts, so this is okay I guess (LSU has ~140 starts) 
   */
  int *sorted_best_roots_v; 
  int *sorted_best_roots_start; 
  float *sorted_best_roots_spdup;
  int *already_chosen;
  int s1, s2;
  int best_cur_v;
  int best_cur_start;
  float best_cur_spdup;

  ESL_ALLOC(sorted_best_roots_v,     sizeof(int) * nstarts);
  ESL_ALLOC(sorted_best_roots_start, sizeof(int) * nstarts);
  ESL_ALLOC(sorted_best_roots_spdup, sizeof(int) * nstarts);
  ESL_ALLOC(already_chosen,          sizeof(int) * nstarts);
  esl_vec_ISet(already_chosen, nstarts, FALSE);
  for(s1 = 0; s1 < nstarts; s1++) {
    best_cur_v = -1;
    best_cur_start = -1;
    best_cur_spdup = -eslINFINITY;
    for(s2 = 0; s2 < nstarts; s2++) { 
      if(! already_chosen[s2]) {
	if(best_per_start_spdup[s2] > best_cur_spdup) { 
	  best_cur_v = best_per_start_v[s2];
	  best_cur_start = s2;
	  best_cur_spdup = best_per_start_spdup[s2];
	}
      }
    }
    sorted_best_roots_v[s1] = best_cur_v;
    sorted_best_roots_start[s1] = best_cur_start;
    sorted_best_roots_spdup[s1] = best_cur_spdup;
    already_chosen[best_cur_start] = TRUE;
  }
  for(s1 = 0; s1 < nstarts; s1++) {
    if(cfg->be_verbose) printf("SORTED rank: %d v: %d spdup: %.5f start: %d\n", s1, sorted_best_roots_v[s1], sorted_best_roots_spdup[s1], sorted_best_roots_start[s1]);
    else { ESL_DPRINTF1(("SORTED rank: %d v: %d spdup: %.5f start: %d\n", s1, sorted_best_roots_v[s1], sorted_best_roots_spdup[s1], sorted_best_roots_start[s1])); } 
    ESL_DASSERT1((sorted_best_roots_v[s1] != -1));
  }
  *ret_sorted_best_roots_v = sorted_best_roots_v;

  free(sorted_best_roots_start);
  free(sorted_best_roots_spdup);
  free(already_chosen);
  free(best_per_start_v);
  free(best_per_start_spdup);
  for(v = 0; v < cm->M; v++) free(sorted_fil_vscAA[v]);
  free(sorted_fil_vscAA);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "predict_best_sub_cm_roots(), memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: predict_hybrid_filter_speedup()
 * Date:     EPN, Tue Dec 11 04:56:39 2007
 *
 * Purpose:  Given a CM and scores for a hybrid CYK/Viterbi scan
 *           of target seqs predict the speedup with a hybrid filter,
 *           then if it's faster than the existing best filter in
 *           BestFilterInfo_t object <bf>, update <bf> to hold info 
 *           on the hybrid filter.
 *            
 * Args:     go  - command line options
 *           cfg - cmcalibrate's cfg object, mucho data (probably too much)
 *           errbuf - for printing error messages
 *           cm - the model
 *           fil_hybscA     - [0..i..filN-1] best Foward score in sequence i 
 *           gum_hybA       - [0..cfg->np]   hybrid gumbels for each partition
 *           fil_partA      - [0..i..filN-1] partition of sequence i 
 *           bf             - BestFilterInfo_t object, we'll update this to hold info on a Viterbi or Forward filter strategy 
 * 
 * Returns:  possibly updates BestFilterInfo_t object <bf> to hold info on hybrid filter
 *           eslOK on success;
 *           Other easel status code on an error with errbuf filled with error message.
 *           <ret_getting_faster> set to TRUE if hybrid scanner replaced previous best filter,
 *           FALSE if not.
 */
int
predict_hybrid_filter_speedup(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *fil_hybscA, int *fil_partA, GumbelInfo_t **gum_hybA, BestFilterInfo_t *bf, int *ret_getting_faster)
{
  int    status;
  float  *sorted_fil_hybEA;       /* sorted hybrid E-values, so we can easily choose a threshold */
  float  E;                       /* E-value */
  int    evalue_L;                /* length used for calc'ing E values */
  float  fil_calcs;               /* number of million dp calcs predicted for the hybrid scan */
  float  surv_calcs;              /* number of million dp calcs predicted for the CM scan of filter survivors */
  float  fil_plus_surv_calcs;     /* filter calcs plus survival calcs */ 
  float  nonfil_calcs;            /* number of million dp calcs predicted for a full CM scan */
  float  spdup;                   /* predicted speedups for Viterbi and Forward, and a temporary one */
  int    i, p;                    /* counters */
  int    Fidx;                     /* index in sorted scores that threshold will be set at (1-F) * N */
  float  F = esl_opt_GetReal(go, "--F"); /* fraction of CM seqs we require filter to let pass */
  float  filE = esl_opt_GetReal(go, "--filE"); /* E-value of cutoff for accepting CM seqs */
  int    filN  = esl_opt_GetInteger(go, "--filN"); /* number of sequences we emitted from CM for filter test */

  /* contract checks */
  if(! (gum_hybA[0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_hybrid_filter_speedup(), gumbel stats for hybrid scanner are not valid.\n");
  evalue_L = gum_hybA[0]->L;
  /* contract checks specific to case when there is more than 1 partition */
  if(cfg->cmstatsA[cfg->ncm-1]->np != 1) { 
    for(p = 1; p < cfg->cmstatsA[cfg->ncm-1]->np; p++) {
      if(evalue_L != gum_hybA[p]->L) ESL_FAIL(eslEINCOMPAT, errbuf, "predict_hybrid_filter_speedup(), partition %d db length (%d) for hybrid gumbel stats differ than from partition 1 hybrid db length (%d).\n", p, gum_hybA[p]->L, evalue_L);
    }
  }

  Fidx  = (int) ((1. - F) * (float) filN);

  /* convert bit scores to E-values and sort them */
  ESL_ALLOC(sorted_fil_hybEA, sizeof(float) * filN);
  for(i = 0; i < filN; i++) { 
    p = fil_partA[i];
    sorted_fil_hybEA[i] = RJK_ExtremeValueE(fil_hybscA[i], gum_hybA[p]->mu, gum_hybA[p]->lambda); 
  }
  esl_vec_FSortDecreasing(sorted_fil_hybEA, filN);
  E = sorted_fil_hybEA[Fidx];

  for(i = 0; i < filN; i++) ESL_DPRINTF1(("HYBRID i: %4d E: %10.4f\n", i, sorted_fil_hybEA[i]));
  
  /* calculate speedup */
  /* E is expected number of hybrid hits with score above threshold (Fidx'th best score) at least sc in sequence DB of length evalue_L */
  E *= cfg->dbsize / evalue_L;
  /* E is now expected number of CP9 Viterbi or Forward hits with score above threshold in sequence DB of length cfg->dbsize */
  fil_calcs  = cfg->hsi->hybrid_ncalcs; /* fil_calcs is millions of DP calcs for hybrid scan of 1 residue */
  fil_calcs *= cfg->dbsize;             /* fil_calcs is millions of DP calcs for hybrid scan of length cfg->dbsize */
  surv_calcs = E *     /* number of hits expected to survive filter */
    (2. * cm->W - (cfg->avglen[0])) * /* average length of surviving fraction of db from a single hit (cfg->avglen[0] is avg subseq len in  subtree rooted at v==0, from QDB calculation, so slightly inappropriate b/c we're concerned with hybrid hits here) */
    cfg->hsi->full_cm_ncalcs; /* number of calculations for full CM scan of 1 residue */
  fil_plus_surv_calcs = fil_calcs  + surv_calcs; /* total number of millions of DP calculations expected using the hybrid filter for scan of cfg->dbsize (logsum corrected, Forward calcs *= 2.) */
  nonfil_calcs = cfg->hsi->full_cm_ncalcs;      /* total number of millions of DP calculations for full CM scan of 1 residue */
  nonfil_calcs *= cfg->dbsize;                  /* now nonfil-calcs corresponds to cfg->dbsize */
  spdup = nonfil_calcs / fil_plus_surv_calcs;
  
  if(cfg->be_verbose) printf("HYBRID E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup);
  else { ESL_DPRINTF1(("HYBRID E: %10.4f filt: %10.4f surv: %10.4f sum: %10.4f full CM: %10.4f spdup %10.4f\n", E, fil_calcs, surv_calcs, fil_plus_surv_calcs, nonfil_calcs, spdup)); }

  if(spdup > (bf->full_cm_ncalcs / bf->fil_plus_surv_ncalcs)) { /* hybrid is best filter strategy so far */
    if((status = SetBestFilterInfoHybrid(bf, errbuf, cm->M, filE, F, filN, cfg->dbsize, nonfil_calcs, E, fil_calcs, fil_plus_surv_calcs, cfg->hsi, cfg->cmstatsA[cfg->ncm-1]->np, gum_hybA)) != eslOK) return status;
    *ret_getting_faster = TRUE;
  }
  else *ret_getting_faster = FALSE;

  return eslOK;

 ERROR:
  return status; 
}

/* update_cutoffs()
 * Update the cfg->cutoffA array to have the bit score cutoff for each partition
 * for the 'current' cm (number ncm-1).
 */
static int
update_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int fthr_mode)
{
  double         tmp_K;          /* used for recalc'ing Gumbel stats for DB size */
  double         e_cutoff;       /* E-value cutoff for each partition */
  int            sc_cutoff;      /* bit score cutoff for each partition */
  int            p;              /* partition index */
  double         mu;             /* mu for a requested db size (which is 1Mb unless --db enabled) */
  int            revert_to_default_filE; /* if --ga, --nc, or --tc enabled but CM does not have GA, NC or TC cutoff in the CM file,
					  * and we've already calibrated at least 1 CM in this CM file, pretend like --ga, --nc, or --tc,
					  * was not enabled by reverting to the default --fil-E value of 0.1
					  */
  revert_to_default_filE = FALSE;
  /* if --ga, --nc, or --tc: 
   * cfg->filE:                set as E-value (in 1Mb db) that the GA, NC, or TC bit score corresponds to for each partition in this fthr_mode
   * cfg->cutoffA[0..p..np-1]: set as GA, NC, or TC bit score for all p for this fthr_mode
   */
  if ((esl_opt_GetBoolean(go, "--ga")) || (esl_opt_GetBoolean(go, "--nc")) || (esl_opt_GetBoolean(go, "--tc"))) {
    if(esl_opt_GetBoolean(go, "--ga")) { 
      if(! (cm->flags & CMH_GA)) {
	if(cfg->ncm > 1) revert_to_default_filE = TRUE;
	else             cm_Fail("--ga enabled but first CM in CM file does not have a Rfam GA cutoff.");
      }
      else sc_cutoff = cm->ga;
    }
    if(esl_opt_GetBoolean(go, "--nc")) { 
      if(! (cm->flags & CMH_NC)) {
	if(cfg->ncm > 1) revert_to_default_filE = TRUE;
	else             cm_Fail("--nc enabled but first CM in CM file does not have a Rfam NC cutoff.");
      }
      else sc_cutoff = cm->nc;
    }
    if(esl_opt_GetBoolean(go, "--tc")) { 
      if(! (cm->flags & CMH_TC)) {
	if(cfg->ncm > 1) revert_to_default_filE = TRUE;
	else             cm_Fail("--tc enabled but first CM in CM file does not have a Rfam TC cutoff.");
      }
      else sc_cutoff = cm->tc;
    }
    if(! revert_to_default_filE) { /* we've set sc_cutoff above, now determine e_cutoff for each partition */
      for (p = 0; p < cfg->np; p++) 
	cfg->cutoffA[p] = sc_cutoff; /* either cm->ga, cm->nc, or cm->tc as set above */
      return eslOK; /* we're done */
    }
  }

  if(esl_opt_GetBoolean(go, "--all")) {
    for (p = 0; p < cfg->np; p++)
      cfg->cutoffA[p] = -eslINFINITY;
  }
  else { /* either none of: --fil-E, --ga, --nc, --tc, --all were enabled, or --fil-E was enabled, 
	  * or --ga, --nc, --tc were enabled, but CM does not have cm->ga, cm->nc, or cm->tc and CM is not first in file */
    e_cutoff = esl_opt_GetReal(go, "--fil-E"); 
    for (p = 0; p < cfg->np; p++) {
      /* first determine mu based on db_size */
      tmp_K = exp(cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->mu * cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda) / 
	cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->L;
      mu = log(tmp_K  * ((double) cfg->dbsize)) / cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda;
      /* Now determine bit score */
      cfg->cutoffA[p] = mu - (log(e_cutoff) / cfg->cmstatsA[cfg->ncm-1]->gumAA[fthr_mode][p]->lambda);
    }
  }
  return eslOK;
}  

/* cm_fit_histograms()
 * We want gumbels for each cm state we can do a legal local begin into.
 * Call fit_histogram() for each such state.
 */
static int
cm_fit_histograms(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, 
		  float **vscA, int nscores, int p)
{
  int status;
  
  if(cfg->vmuAA[p]     != NULL) free(cfg->vmuAA[p]);
  if(cfg->vlambdaAA[p] != NULL) free(cfg->vlambdaAA[p]);
  
  ESL_ALLOC(cfg->vmuAA[p],     sizeof(double) * 1);
  ESL_ALLOC(cfg->vlambdaAA[p], sizeof(double) * 1);
  if((status = fit_histogram(go, cfg, errbuf, vscA[0], nscores, &(cfg->vmuAA[p][0]), &(cfg->vlambdaAA[p][0]))) != eslOK) return status;

  return eslOK;

 ERROR:
  return status;
}

#endif
#if 0

/* Function:  cmcalibrate_cm_gumbel_results_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            results for CM scan for cmcalibrate.
 * Incept:    EPN, Thu Dec  6 16:44:17 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cmcalibrate_cm_gumbel_results_MPIPack()> will need 
 *            to pack it's results in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cmcalibrate_cm_gumbel_results_MPIPack()>.
 */
int
cmcalibrate_cm_gumbel_results_MPIPackSize(float **vscAA, int nseq, int M, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);   n += sz;      if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size(M, MPI_FLOAT, comm, &sz); n += nseq*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmcalibrate_cm_gumbel_results_MPIPack()
 * Synopsis:  Packs CM vscAA scores into MPI buffer.
 * Incept:    EPN, Thu Dec  6 16:47:58 2007
 *
 * Purpose:   Packs <vscAA> into an MPI packed message 
 *            buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 *            Note: <vscAA> is a 2D array, vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan (CYK/Inside) of sequence i.
 *            But we send it as a 1D array, vscA, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <tr>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <rnode> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cmcalibrate_cm_gumbel_results_MPIPack(float **vscAA, int nseq, int M, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int i,v,idx;
  float *vscA = NULL;

  ESL_DPRINTF2(("cmcalibrate_cm_gumbel_results_MPIPack(): ready.\n"));

  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  idx = 0;
  for(v = 0; v < M; v++) 
    for(i = 0; i < nseq; i++)
      vscA[idx++] = vscAA[v][i];

  status = MPI_Pack((int *) &(nseq), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(vscA,  (M*nseq), MPI_FLOAT,  buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cmcalibrate_cm_gumbel_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;

 ERROR:
  if(vscA  != NULL) free(vscA);
  return status;
}

/* Function:  cmcalibrate_cm_gumbel_results_MPIUnpack()
 * Synopsis:  Unpacks <vscAA> from an MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:10:20 2007
 *
 * Purpose:   Unpack a newly allocated set of scores <vscAA> from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 *            Note: We return <ret_vscAA> as a 2D array, 
 *            ret_vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan (CYK/Inside) of sequence i.
 *            But vscA is sent as a 1D array, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_tr>
 *            contains a newly allocated parsetree, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_vscAA> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cmcalibrate_cm_gumbel_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, int M, float ***ret_vscAA, int *ret_nseq)
{
  int status;
  float  *vscA  = NULL;
  float **vscAA = NULL;
  int nseq = 0;
  int i, v, idx;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  status = MPI_Unpack (buf, n, pos, vscA, (M*nseq), MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(vscAA, sizeof(float *) * (M));
  idx = 0;
  for(v = 0; v < M; v++) { 
    ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    for(i = 0; i < nseq; i++)
      vscAA[v][i] = vscA[idx++];
  }
  ESL_DASSERT1((idx == (M*nseq)));

  free(vscA);
  *ret_vscAA = vscAA;
  *ret_nseq = nseq;
  return eslOK;

 ERROR:
  if(vscA  != NULL) free(vscA);
  if(vscAA != NULL) { 
    for(i = 0; i < nseq; i++)
      free(vscAA[i]);
    free(vscAA);
  }
  *ret_vscAA = NULL;
  *ret_nseq = 0;
  return status;
}


/* Function:  cmcalibrate_cp9_filter_results_hyb_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack 
 *            CP9 filter results for cmcalibrate with --hybrid.
 *            enabled. 
 *           
 *            Differs from cmcalibrate_cp9_filter_results_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            is packed. these scores are eventually used to calculate
 *            Gumbels for each possible sub CM state.
 *
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPackSize() funcs above.
 *            
 * Incept:    EPN, Wed Dec 12 16:30:20 2007
 *           
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cmcalibrate_cp9_filter_results_hyb_MPIPack()>.
 */
int
cmcalibrate_cp9_filter_results_hyb_MPIPackSize(int nseq, int M, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);        n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for nseq */
  status = MPI_Pack_size(M*nseq, MPI_FLOAT, comm, &sz); n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for vscAA (we'll send it as a 1D array of M * nseq floats */
  status = MPI_Pack_size(nseq, MPI_FLOAT, comm, &sz);   n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for vit_cp9scA, fwd_cp9scA */
  status = MPI_Pack_size(nseq, MPI_INT,   comm, &sz);   n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for partA */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmcalibrate_cp9_filter_results_hyb_MPIPack()
 * Synopsis:  Packs cmcalibrate CP9 filter results into MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPack() funcs above.
 * 
 * Incept:    EPN, Wed Dec 12 16:36:02 2007
 *
 *            Differs from cmcalibrate_cp9_filter_results_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            is packed. these scores are eventually used to calculate
 *            Gumbels for each possible sub CM state.
 *
 *            Note: <vscAA> is a 2D array, vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan of sequence i.
 *            But we send it as a 1D array, vscA, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
 *
 */
int
cmcalibrate_cp9_filter_results_hyb_MPIPack(float **vscAA, float *vit_cp9scA, float *fwd_cp9scA, int *partA, int nseq, int M, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int i,v,idx;
  float *vscA = NULL;

  ESL_DPRINTF2(("cmcalibrate_cp9_filter_results_hyb_MPIPack(): ready.\n"));

  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  idx = 0;
  for(v = 0; v < M; v++) 
    for(i = 0; i < nseq; i++)
      vscA[idx++] = vscAA[v][i];

  status = MPI_Pack((int *) &(nseq), 1,        MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(vscA,            (M*nseq), MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(vit_cp9scA,      nseq,     MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(fwd_cp9scA,      nseq,     MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(partA,           nseq,     MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cmcalibrate_cp9_filter_results_hyb_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;

 ERROR:
  if(vscA  != NULL) free(vscA);
  return status;
}

/* Function:  cmcalibrate_cp9_filter_results_hyb_MPIUnpack()
 * Synopsis:  Unpacks cmcalibrate cp9 filter results from an MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIUnpack() funcs above.
 * Incept:    EPN, Wed Dec 12 16:38:15 2007
 *
 *            Differs from cmcalibrate_cp9_filter_results_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            is packed. these scores are eventually used to calculate
 *            Gumbels for each possible sub CM state.
 *
 *            Note: We return <ret_vscAA> as a 2D array, 
 *            ret_vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan of sequence i.
 *            But vscA is sent as a 1D array, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
 *
 */
int
cmcalibrate_cp9_filter_results_hyb_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, int M, float ***ret_vscAA, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, int **ret_partA, int *ret_nseq)
{
  int status;
  float  *vit_cp9scA  = NULL;
  float  *fwd_cp9scA  = NULL;
  int    *partA       = NULL;
  float  *vscA  = NULL;
  float **vscAA = NULL;
  int nseq = 0;
  int i, v, idx;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  status = MPI_Unpack (buf, n, pos, vscA, (M*nseq), MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(vscAA, sizeof(float *) * (M));
  idx = 0;
  for(v = 0; v < M; v++) { 
    ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    for(i = 0; i < nseq; i++)
      vscAA[v][i] = vscA[idx++];
  }
  ESL_DASSERT1((idx == (M*nseq)));

  ESL_ALLOC(vit_cp9scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, vit_cp9scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(fwd_cp9scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, fwd_cp9scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(partA,      sizeof(int) * nseq);
  status = MPI_Unpack (buf, n, pos, partA, nseq, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_vscAA = vscAA;
  *ret_vit_cp9scA = vit_cp9scA;
  *ret_fwd_cp9scA = fwd_cp9scA;
  *ret_partA      = partA;
  *ret_nseq = nseq;
  return eslOK;

 ERROR:
  if(vit_cp9scA != NULL) free(vit_cp9scA);
  if(vit_cp9scA != NULL) free(fwd_cp9scA);
  if(partA      != NULL) free(partA);

  if(vscA  != NULL) free(vscA);
  if(vscAA != NULL) { 
    for(i = 0; i < nseq; i++)
      free(vscAA[i]);
    free(vscAA);
  }
  *ret_vscAA = NULL;
  *ret_nseq = 0;
  return status;
}

#endif

#if 0

/* Function:  summarize_alignment()
 * Incept:    
 *
 * Purpose:   Summarize alignment statistics to varying extents
 *            based on command-line options.
 */
int
summarize_alignment(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w) 
{
  /* HERE: do HMM banded alignment stats
   * sample N=100 seqs, and calculate posteriors, determine new
   * number of CYK DP calcs AND CP9 F/B calcs to get bands. */
  int status;
  float dpc;  /* # DP calcs for non-banded alignment of consensus */
  CMConsensus_t *con = NULL;            /* consensus info for the CM */
  ESL_SQ *csq = NULL;
  float t_dc; /* user seconds time for D&C alignment */
  float t_hb; /* user seconds time for HMM banded alignment */
  float mc_s; /* million calcs/second */
  float size_limit = esl_opt_GetReal(go, "--mxsize");

  /* Create and align consensus sequence for D&C stats */
  CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &con);
  if((csq = esl_sq_CreateFrom("consensus", con->cseq, NULL, NULL, NULL)) == NULL)
    { status = eslEMEM; goto ERROR; }
  esl_sq_Digitize(cm->abc, csq);
  dpc = count_align_dp_calcs(cm, csq->n) / 1000000.;

  /* cyk inside (score only) */
  esl_stopwatch_Start(w);
  /*CYKDivideAndConquer(cm, csq->dsq, csq->n, 0, 1, csq->n, NULL, NULL, NULL);*/
  CYKInsideScore(cm, csq->dsq, csq->n, 0, 1, csq->n, 
		 NULL, NULL); /* don't do QDB mode */
  esl_stopwatch_Stop(w);
  t_dc = w->user;
  mc_s = dpc / t_dc;

  /* HMM banded */
  /* Emit N seqs, and align them, to get total time up to reasonable level,
   * and to average out tightness of bands */
  int N = esl_opt_GetInteger(go, "-N");
  seqs_to_aln_t *seqs_to_aln = NULL;
  ESL_SQ **sq = NULL;
  ESL_ALLOC(sq, sizeof(ESL_SQ *) * N);
  int L; 
  float L_avg = 0.; 
  int i;
  for(i = 0; i < N; i++)
    {
      if((status = EmitParsetree(cm, errbuf, r, "seq", TRUE, NULL, &(sq[i]), &L)) != eslOK) goto ERROR;
      /*esl_sqio_Write(stdout, sq[i], eslSQFILE_FASTA);*/
      L_avg += L;
    }
  L_avg /= (float) N;
  cm->align_opts |= CM_ALIGN_HBANDED;
  esl_stopwatch_Start(w);
  seqs_to_aln = CreateSeqsToAlnFromSq(sq, N, FALSE);
  if((status = DispatchAlignments(cm, errbuf, seqs_to_aln, NULL, NULL, 0, 0, 0, TRUE, NULL, size_limit)) != eslOK) goto ERROR;
  esl_stopwatch_Stop(w);
  t_hb = w->user / (float) N;
  FreeSeqsToAln(seqs_to_aln);

  fprintf(ofp, "#\n");
  fprintf(ofp, "#\t\t\t Alignment statistics:\n");
  fprintf(ofp, "#\t\t\t %7s %6s %6s %8s %8s %8s\n",             "alg",     "Mc",     "L",     "Mc/s",    "s/seq",        "accel");
  fprintf(ofp, "#\t\t\t %7s %6s %6s %8s %8s %8s\n",             "-------", "------", "------","--------", "--------", "--------");
  /*mc_s = dpc / t_dc; */
  fprintf(ofp, " \t\t\t %7s %6.1f %6d %8.1f %8.3f %8s\n",       "cyk",      dpc,      csq->n,  mc_s,      t_dc,       "-");
  fprintf(ofp, " \t\t\t %7s %6s %6.0f %8s %8.3f %8.2f\n",       "hb cyk",   "?",      L_avg,   "?",       t_hb,       (t_dc/t_hb));

  esl_sq_Destroy(csq);
  FreeCMConsensus(con);
  return eslOK;
 ERROR:
  cm_Fail("ERROR code %d in summarize_stats().", status);
  return status; /* NOTREACHED */
}

/* Function: count_align_dp_calcs()
 * Date:     EPN, Wed Aug 22 09:08:03 2007
 *
 * Purpose:  Count all non-d&c inside DP calcs for a CM 
 *           alignment of a seq of length L. Similar to cm_dpsmall.c's
 *           CYKDemands() but takes into account number of
 *           transitions from each state, and is concerned
 *           with a scanning dp matrix, not an alignment matrix.
 *
 * Args:     cm     - the model
 *           L      - length of sequence
 *
 * Returns: (float) the total number of DP calculations.
 */
float count_align_dp_calcs(CM_t *cm, int L)
{
  int v, j;
  float dpcalcs = 0.;
  float dpcalcs_bif = 0.;
  
  float  dpcells     = 0.;
  float  dpcells_bif = 0.;

  dpcells = (L+2) * (L+1) * 0.5; /* fillable dp cells per state (deck) */
  for (j = 0; j <= L; j++)
    dpcells_bif += (j+2) * (j+1) * .5;
  dpcalcs_bif = CMCountStatetype(cm, B_st) * dpcells_bif; /* no choice of transitions */
  for(v = 0; v < cm->M; v++)
    if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
      dpcalcs += dpcells * cm->cnum[v]; /* cnum choices of transitions */

  return dpcalcs + dpcalcs_bif;
}
#endif

/* EPN, Sun Jan 20 11:14:02 2008 
 * from cmcalibrate.c: revision 2305
 */
#if 0

/* Function: estimate_workunit_time()
 * Date:     EPN, Thu Nov  1 17:57:20 2007
 * 
 * Purpose:  Estimate time req'd for a cmcalibrate workunit
 *
 * Returns:  eslOK on success;
 */
void
estimate_workunit_time(const ESL_GETOPTS *go, const struct cfg_s *cfg, int nseq, int L, int gum_mode)
{
  /* these are ballparks for a 3 GHz machine with optimized code */
  float cyk_megacalcs_per_sec = 275.;
  float ins_megacalcs_per_sec =  75.;
  float fwd_megacalcs_per_sec = 175.;
  float vit_megacalcs_per_sec = 380.;
  
  float seconds = 0.;

  if(! esl_opt_IsDefault(go, "--gum-L")) L = ESL_MAX(L, esl_opt_GetInteger(go, "--gum-L")); /* minimum L we allow is 2 * cm->W (L is sent into this func as 2 * cm->W), this is enforced silently (!) */

  switch(gum_mode) { 
  case GUM_CM_LC: 
  case GUM_CM_GC: 
    seconds = cfg->gum_cm_ncalcs * (float) L * (float) nseq / cyk_megacalcs_per_sec;
    break;
  case GUM_CM_LI:
  case GUM_CM_GI:
    seconds = cfg->gum_cm_ncalcs * (float) L * (float) nseq / ins_megacalcs_per_sec;
    break;
  case GUM_CP9_LV: 
  case GUM_CP9_GV: 
    seconds = cfg->cp9_ncalcs * (float) L * (float) nseq / vit_megacalcs_per_sec;
    break;
  case GUM_CP9_LF: 
  case GUM_CP9_GF: 
    seconds = cfg->cp9_ncalcs * (float) L * (float) nseq / fwd_megacalcs_per_sec;
    break;
  }
  printf("Estimated time for this workunit: %10.2f seconds\n", seconds);

  return;
}
#endif

/* EPN, Sun Jan 20 14:19:42 2008
 * Decided to remove the --enforce options. 
 * Here's the related functions.
 * These are revision 2305.
 */
#if 0
/*extern void  ConfigCMEnforce(CM_t *cm);*/
/*extern void  ConfigLocalEnforce(CM_t *cm, float p_internal_start, float p_internal_exit);*/
/*extern int   EnforceSubsequence(CM_t *cm);*/
/*extern float EnforceScore(CM_t *cm);*/
/*extern int   EnforceFindEnfStart(CM_t *cm, int enf_cc_start);*/
/*extern void  CPlan9SWConfigEnforce(CP9_t *hmm, float pentry, float pexit, int enf_start_pos, int enf_end_pos);*/
/*extern void  CP9EnforceHackMatchScores(CP9_t *cp9, int enf_start_pos, int enf_end_pos);*/

/*
 * Function: ConfigCMEnforce
 * Date:     EPN, Wed Feb 14 12:57:21 2007
 * Purpose:  Configure a CM for enforcing a subsequence for search or 
 *           alignment. 
 * 
 * Args:     CM           - the covariance model
 */
void
ConfigCMEnforce(CM_t *cm)
{
  int do_build_cp9  = FALSE;
  int enf_start_pos;            /* consensus left position node enf_start emits to   */
  int enf_end_pos;              /* consensus left position node enf_end   emits to   */
  int enf_end;                  /* last node we're enforcing                         */
  CMEmitMap_t *emap;            /* consensus emit map for the CM, used iff enforcing */
  float nonenf_sc;              /* score of cm->enfseq we're about to enforce before *
				 * we reparameterize the CM                          */
  float enf_sc;                 /* score of cm->enfseq subseq after CM is            *
				 * reparameterized to enforce it                     */

  /* Contract checks */
  if(!(cm->config_opts & CM_CONFIG_ENFORCE))
    cm_Fail("ERROR in ConfigCMEnforce() trying to enforce a subsequence but CM_CONFIG_ENFORCE flag is down.");
  if(cm->flags & CM_ENFORCED)
    cm_Fail("ERROR in ConfigCMEnforce() trying to enforce a subsequence but CM_IS_ENFORCED flag is up.");
  /* Can't enforce in RSEARCH mode yet */  
  if(cm->flags & CM_IS_RSEARCH)
    cm_Fail("ERROR in ConfigCMEnforce() trying to enforce a subsequence in RSEARCH mode, not yet implemented.");
  /* Can't enforce in sub mode */  
  if(cm->align_opts & CM_ALIGN_SUB)
    cm_Fail("ERROR in ConfigCMEnforce() can't enforce a subsequence in sub alignment mode.");

  /* First, get the score of the enforced subseq for the non-enforced model */
  nonenf_sc = EnforceScore(cm);

  /* IMPORTANT: if CM has local begins, make it global, we'll relocalize 
   * it later based on cm->config_opts, cm->search_opts, and/or cm->align_opts,
   * we need to do this so we can build a CP9 (which can't be done with local CMs yet)*/
  if(cm->flags & CMH_LOCAL_BEGIN)
    ConfigGlobal(cm);

  /* Enforce the sequence */
  EnforceSubsequence(cm);

  /* if we have a CP9, free it, and build a new one, (this one will automatically
   * have the subseq enforced b/c it's built from the reparam'ized CM) */
  if(cm->flags & CMH_CP9)
    {
      FreeCPlan9(cm->cp9);
      cm->flags &= ~CMH_CP9; /* drop the CP9 flag */
      do_build_cp9 = TRUE;
    }
  else if (cm->config_opts & CM_CONFIG_ENFORCEHMM)
    {
      /* we didn't have a CP9 before, but we need one now */
      do_build_cp9 = TRUE;
    }
  if(do_build_cp9)
    {
      if(!(build_cp9_hmm(cm, &(cm->cp9), &(cm->cp9map), 
			TRUE, /* b/c we're enforcing, check CP9 mirrors CM */
			0.0001, 0)))
	cm_Fail("Couldn't build a CP9 HMM from the CM\n");
      cm->flags |= CMH_CP9; /* raise the CP9 flag */
    }

  /* Configure the CM for local alignment . */
  if (cm->config_opts & CM_CONFIG_LOCAL)
    { 
      ConfigLocalEnforce(cm, cm->pbegin, cm->pend); /* even in local we require each parse 
						     * go through the enforced subseq */
      CMLogoddsify(cm);
    }
  /* Possibly configure the CP9 for local alignment
   * Note: CP9 local/glocal config does not necessarily match CM config 
   *       in fact cmsearch default is local CM, glocal CP9 */
  if((cm->flags & CMH_CP9) && (cm->config_opts & CM_CONFIG_HMMLOCAL))
    {
      /* Set up the CP9 locality to enforce a subseq */
      emap = CreateEmitMap(cm); 
      enf_end = cm->enf_start + strlen(cm->enf_seq) - 1;
      enf_start_pos = emap->lpos[cm->enf_start];
      enf_end_pos   = emap->lpos[enf_end];
      FreeEmitMap(emap);
      CPlan9SWConfigEnforce(cm->cp9, cm->pbegin, cm->pbegin, enf_start_pos, enf_end_pos);
      CP9Logoddsify(cm->cp9);
    }
     
  if(cm->config_opts & CM_CONFIG_ENFORCEHMM)
    {
      if(!(cm->flags & CMH_CP9))
	cm_Fail("ERROR trying to configure the HMM for naive enforcement, but the cm's CMH_CP9 flag is down.\n");
      /* We make the HMM ignorant of any sequence conservation besides
       * the enforced subseq. This way ALL subseqs with the enforced
       * subseq will be recognized as high scoring by the HMM and 
       * be passed to the CM (if filtering (which is default in this mode)).
       * To achieve this, make all emissions (match and insert) score 0,
       * except for the few match emissions that model the enforced subseq */
      emap = CreateEmitMap(cm); 
      enf_end = cm->enf_start + strlen(cm->enf_seq) - 1;
      enf_start_pos = emap->lpos[cm->enf_start];
      enf_end_pos   = emap->lpos[enf_end];
      FreeEmitMap(emap);
      CP9EnforceHackMatchScores(cm->cp9, enf_start_pos, enf_end_pos);
    }	

  /* Determine the score of the enforced subseq for the enforced model */
  enf_sc = EnforceScore(cm);
  
  cm->enf_scdiff = enf_sc - nonenf_sc;
  cm->flags |= CM_ENFORCED; /* raise the enforced flag */

  if(cm->flags & CMH_SCANMATRIX)
    cm->flags &= ~CMH_SCANMATRIX; /* enforcement invalidates ScanMatrix */

  return; 
}

/* EPN, Thu Jan  4 10:10:07 2007
 * 
 * Function: ConfigLocalEnforce
 * 
 * Purpose:  Given a CM with valid cm->enf_start and cm->enf_seq variables,
 *           modify local entries and exits so that the nodes starting
 *           at cm->enf_start and going to cm->enf_start + strlen(cm->enf_seq)
 *           must be entered, i.e. disallow any local pass that omits them.
 * 
 * Args:     CM           - the covariance model
 *           p_internal_start - total prob of a local begin to spread 
 *           p_internal_exit  - total prob of a local end to spread
 */
void
ConfigLocalEnforce(CM_t *cm, float p_internal_start, float p_internal_exit)
{
  int v;			/* counter over states */
  int nd;			/* counter over nodes */
  int nstarts;			/* number of possible internal starts */
  int enf_start_pos;            /* consensus left position node enf_start emits to */
  int enf_end_pos;              /* consensus left position node enf_end   emits to */
  int nexits;			/* number of possible internal ends */
  float denom;
  CMEmitMap_t *emap;            /* consensus emit map for the CM */
  int enf_end;

  /* Contract checks */
  if(cm->enf_seq == NULL || cm->enf_start == 0)
    cm_Fail("ERROR, in ConfigLocalEnforce, but no subseq to enforce.\n");
  if(cm->flags & CMH_LOCAL_BEGIN)
    cm_Fail("ERROR in ConfigLocalEnforce() CMH_LOCAL_BEGIN flag already up.\n");
  if(cm->flags & CMH_LOCAL_END)
    cm_Fail("ERROR in ConfigLocalEnforce() CMH_LOCAL_END flag already up.\n");

  enf_end = cm->enf_start + strlen(cm->enf_seq) - 1;
  /* We want every parse to go through the MATL stretch from enf_start
   * to enf_end. To enforce this we disallow local begin and ends that
   * would allow parses to miss these nodes. */
  for(nd = cm->enf_start; nd <= enf_end; nd++)
    {
      if(cm->ndtype[nd] != MATL_nd)
	cm_Fail("ERROR, trying to enforce a non-MATL stretch (node: %d not MATL).\n", nd);
    }
  emap = CreateEmitMap(cm); /* diff from ConfigLocalEnds() */
  enf_start_pos = emap->lpos[cm->enf_start];
  enf_end_pos   = emap->lpos[enf_end];

  /* The following code is copied from ConfigLocal() and ConfigLocalEnds()
   * with modification to disallow local begins before enf_start and 
   * disallow exits from between closest_start and enf_end. This 
   * implementation sets local entry to the first node as 1-p_internal_start,
   * and end from final node as 1-p_internal_exit. */

  /*****************************************************************
   * Internal entry.
   *****************************************************************/
  /* Count "internal" nodes: MATP, MATL, MATR, and BIF nodes.
   * Ignore all start nodes, and also node 1 (which is always the
   * "first" node and gets an entry prob of 1-p_internal_start).
   */
  nstarts = 0;
  for (nd = 2; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
      if(emap->lpos[nd] <= enf_start_pos &&
	 emap->rpos[nd] >= enf_end_pos) /* diff from ConfigLocalEnds() */
	nstarts++;
  }

  /* Zero everything.
   */
  for (v = 0; v < cm->M; v++)  cm->begin[v] = 0.;

  /* Erase the previous transition p's from node 0. The only
   * way out of node 0 is going to be local begin transitions
   * from the root v=0 directly to MATP_MP, MATR_MR, MATL_ML,
   * and BIF_B states.
   */
  for (v = 0; v < cm->cnum[0]; v++)  cm->t[0][v] = 0.;

  /* Node 1 gets prob 1-p_internal_start.
   */
  cm->begin[cm->nodemap[1]] = 1.-p_internal_start;

  /* Remaining nodes share p_internal_start.
   */
  for (nd = 2; nd < cm->nodes; nd++) 
    {
      if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	  cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
	{
	  if(emap->lpos[nd] <= enf_start_pos &&
	     emap->rpos[nd] >= enf_end_pos) /* diff from ConfigLocalEnds() */
	    {
	      /*printf("enabling local begin into nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	      cm->begin[cm->nodemap[nd]] = p_internal_start/(float)nstarts;
	    }
	  else
	    ;/*printf("NOT enabling local begin into nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	}
    }
  cm->flags |= CMH_LOCAL_BEGIN;
  
  /*****************************************************************
   * Internal exit.
   *****************************************************************/
  /* Count internal nodes MATP, MATL, MATR, BEGL, BEGR that aren't
   * adjacent to END nodes.
   */
  nexits = 0;
  for (nd = 1; nd < cm->nodes; nd++) {
    if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	 cm->ndtype[nd] == BEGR_nd) && 
	cm->ndtype[nd+1] != END_nd)
      if(emap->lpos[nd] >= enf_end_pos || 
	 emap->rpos[nd] <  enf_start_pos) /* diff from ConfigLocalEnds() */
	nexits++;
  } 
  /* Spread the exit probability across internal nodes.
   * Currently does not compensate for the decreasing probability
   * of reaching a node, the way HMMER does: therefore the probability
   * of exiting at later nodes is actually lower than the probability 
   * of exiting at earlier nodes. This should be a small effect.
   */
  for (v = 0; v < cm->M; v++) cm->end[v] = 0.;
  for (nd = 1; nd < cm->nodes; nd++) {
    if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	 cm->ndtype[nd] == BEGR_nd) && 
	cm->ndtype[nd+1] != END_nd)
      {
	v = cm->nodemap[nd];
	if(emap->lpos[nd] >= enf_end_pos || 
	   emap->rpos[nd] <  enf_start_pos) /* diff from ConfigLocalEnds() */
	  {
	    /*printf("enabling local end from nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	    cm->end[v] = p_internal_exit / (float) nexits;
	  }
	else
	  {
	    ;/*printf("NOT enabling local end from nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	  }
	/* renormalize the main model transition distribution,
	 * it's important to do this for all states that
	 * may have had a local end possible prior to this function call*/
	denom = esl_vec_FSum(cm->t[v], cm->cnum[v]);
	denom += cm->end[v];
	esl_vec_FScale(cm->t[v], cm->cnum[v], 1./denom);
      }
  }
  cm->flags |= CMH_LOCAL_END;
  FreeEmitMap(emap);

  return;
}

/*******************************************************************************
 * Function: EnforceSubsequence()
 * Date:     EPN, Thu Jan  4 10:13:08 2007
 * Purpose:  Modify CM probabilities so that if a particular subsequence (cm->enf_subseq)
 *           is not emitted, a big bit score penalty is incurred. Specifically designed
 *           for enforcing the telomerase RNA template sequence.
 */
int  
EnforceSubsequence(CM_t *cm)
{
  int status;
  int nd;
  float small_chance = 1e-15; /* any parse not including the enforced path includes
			       * an emission or transition with a -45 bit score */
  int   enf_end;
  int v;
  int a;
  float *nt;
  ESL_SQ *enf_sq = NULL;     /* We'll fill this with enf_seq and digitize it */

  ESL_ALLOC(nt, sizeof(float) * cm->abc->K);

  enf_end = cm->enf_start + strlen(cm->enf_seq) - 1;
  /*printf("in EnforceSubsequence, start posn: %d cm->enf_seq: %s\n", cm->enf_start, cm->enf_seq);*/
  for(nd = (cm->enf_start-1); nd <= enf_end; nd++)
    {
      if(cm->ndtype[nd] != MATL_nd)
	cm_Fail("ERROR, trying to enforce a non-MATL stretch (node: %d not MATL).\n", nd);
    }

  /* Go through each node and enforce the template by changing the
   * emission and transition probabilities as appropriate. */

  /* First deal with node before cm->enf_start, we want to ensure that cm->enf_start is
   * entered. We know cm->enf_start - 1 and cm->enf_start are both MATL nodes */
  nd = cm->enf_start - 1;
  v  = cm->nodemap[nd];       /* MATL_ML*/
  cm->t[v][2] = small_chance; /* ML->D  */
  v++;                        /* MATL_D */
  cm->t[v][2] = small_chance; /*  D->D  */
  v++;                        /* MATL_IL*/
  cm->t[v][2] = small_chance; /*  IL->D */

  /* Now move on to the MATL nodes we're enforcing emits the cm->enf_seq */
  enf_sq  = esl_sq_CreateFrom("enforced", cm->enf_seq, NULL, NULL, NULL);
  if(enf_sq == NULL) goto ERROR;
  if((status = esl_sq_Digitize(cm->abc, enf_sq)) != eslOK) goto ERROR;

  for(nd = cm->enf_start; nd <= enf_end; nd++) 
    {
      /*printf("enforcing subseq for node: %d\n", nd);*/
      /* Enforce the transitions, unless we're the last node of the stretch */
      v  = cm->nodemap[nd];       /* MATL_ML*/
      if(nd < enf_end)
	{
	  cm->t[v][0] = small_chance; /* ML->IL */
	  cm->t[v][2] = small_chance; /* ML->D  */
	}
      /* Enforce the emission. Taking into account ambiguities. */
      esl_vec_FSet(nt, cm->abc->K, 0.);
      /*printf("enf_dsq[%d]: %d\n", (nd-cm->enf_start+1), (int) (enf_dsq[(nd-cm->enf_start+1)]));*/
      esl_abc_FCount(cm->abc, nt, enf_sq->dsq[(nd - cm->enf_start + 1)], 1.);
      /* nt is now a count vector norm'ed to 1.0 with relative contributions 
       * of each (A,C,G,U) nucleotides towards the (potentially ambiguous)
       * residue in enf_dsq[(nd-cm->enf_start+1)]) 
       */

      for(a = 0; a < cm->abc->K; a++)
	{
	  /* start out by setting each residue to 'small_chance' */
	  cm->e[v][a] =  small_chance;
	  cm->e[v][a] += nt[a];
	}
    }
  free(nt);
  esl_sq_Destroy(enf_sq);

  CMRenormalize(cm);
  /* new probs invalidate log odds scores */
  cm->flags &= ~CMH_BITS;
  /* Recalc QDBs if they exist */
  if(cm->flags & CMH_QDB) {
    free(cm->dmin);
    free(cm->dmax);
    cm->dmin = NULL;
    cm->dmax = NULL;
    cm->flags &= ~CMH_QDB;
    ConfigQDB(cm);
  }      
  /* free and rebuild scan matrix to correspond to new QDBs, if it exists */
  if(cm->flags & CMH_SCANMATRIX) {
    int do_float = cm->smx->flags & cmSMX_HAS_FLOAT;
    int do_int   = cm->smx->flags & cmSMX_HAS_INT;
    cm_FreeScanMatrixForCM(cm);
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
  }

  CMLogoddsify(cm); /* QDB calculation invalidates log odds scores */

  /*for(nd = cm->enf_start; nd <= enf_end; nd++) 
    {
      v  = cm->nodemap[nd];      
      for(a = 0; a < cm->abc->K; a++)
	printf("cm->e[v:%d][a:%d]: %f sc: %f\n", v, a, cm->e[v][a], cm->esc[v][a]);
    }
  printf("\n");*/

  return eslOK;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return status; /* never reached */
}

/*******************************************************************************
 * Function: EnforceScore()
 * Date:     EPN, Wed Feb 14 16:19:22 2007
 * Purpose:  Determine the subparse score of aligning cm->enfseq to the MATL_ML 
 *           states of consecutive MATL nodes starting at cm->enfstart. This
 *           function can be called before and after enforcing the subseq 
 *           via reparameterization of the relevant nodes, to determine the
 *           score difference of cm->enfseq b/t the non-enforced and enforced
 *           CMs.
 */
float
EnforceScore(CM_t *cm)
{
  int     status;
  ESL_SQ *enf_sq;/* a digitized version of cm->enf_seq */
  int   enf_end; /* last node to be enforced */
  int   nd;      /* node index  */
  int   v;       /* state index */
  int   i;       /* sequence position index */
  float score;   /* score of subparse that starts in first MATL_ML of enforced stretch,
		  * goes through each MATL_ML and emits the enforced residues (which
		  * can be ambiguous). */

  /* Contract check. */
  if(!(cm->config_opts & CM_CONFIG_ENFORCE))
    cm_Fail("ERROR in EnforceScore(), cm->config_opt CM_CONFIG_ENFORCE not raised.\n");

  enf_end = cm->enf_start + strlen(cm->enf_seq) - 1;
  /*printf("in EnforceScore(), start posn: %d cm->enf_seq: %s\n", cm->enf_start, cm->enf_seq);*/
  for(nd = (cm->enf_start-1); nd <= enf_end; nd++)
    {
      if(cm->ndtype[nd] != MATL_nd)
	cm_Fail("ERROR, trying to enforce a non-MATL stretch (node: %d not MATL).\n", nd);
    }

  /* Go through each node and determine the score of the subparse that
   * goes through the nodes that are/will be enforced.  To start, we
   * have to transit to MATL_ML of cm->enf_start from either MATL_ML
   * of MATL_ML or MATL_IL of nd=cm->enf_start-1, but we don't know
   * which.  Can't think of robust way of handling this, current
   * strategy is to take the average of the two transition scores.
   * (this is hacky, but should have small effect on cm->enf_scdiff).
   */
  nd = cm->enf_start - 1;
  v  = cm->nodemap[nd];       /* MATL_ML*/
  score =  (cm->tsc[v][1] + cm->tsc[v+2][1]) / 2;
  /*printf("init v: %d ML->ML: %f IL->ML: %f avg: %f\n", v, cm->tsc[v][1], cm->tsc[(v+2)][1], score);*/

  /* Now move on to the MATL nodes we're enforcing emits the cm->enf_seq */
  enf_sq  = esl_sq_CreateFrom("enforced", cm->enf_seq, NULL, NULL, NULL);
  if(enf_sq == NULL) goto ERROR;
  if((status = esl_sq_Digitize(cm->abc, enf_sq)) != eslOK) goto ERROR;

  for(nd = cm->enf_start; nd <= enf_end; nd++) 
    {
      i = nd - cm->enf_start+1; /* enf_dsq goes 1..(strlen(cm->enf_seq)) 
				 * bordered by sentinels */
      /* Add score for the MATL_ML->MATL_ML transition, 
       * unless we're the last node of the stretch */
      v  = cm->nodemap[nd];       /* MATL_ML*/
      if(nd < enf_end)
	score += cm->tsc[v][1]; /* ML->ML */

      /* Add score for the emission. Taking into account ambiguities. */
      if (enf_sq->dsq[i] < cm->abc->K)
	score += cm->esc[v][enf_sq->dsq[i]];
      else
	score += esl_abc_FAvgScore(cm->abc, enf_sq->dsq[i], cm->esc[v]);

    }
  /*printf("in EnforceScore() returning sc: %f\n", score);*/
  esl_sq_Destroy(enf_sq);

  return score;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return status; /* never reached */
}

/*******************************************************************************
 * Function: EnforceFindEnfStart()
 * Date:     EPN, Fri Feb  9 10:32:44 2007
 * Purpose:  Determine the node cm->enf_start given the consensus column it 
 *           models, and check that it's a MATL node (this requirement could 
 *           be relaxed in the future).
 * Returns:  (int) the CM MATL node index that emits to consensus column 
 *           enf_cc_start. Dies if there's no such node.
 */
int  
EnforceFindEnfStart(CM_t *cm, int enf_cc_start)
{
  CMEmitMap_t *emap;            /* consensus emit map for the CM */
  int enf_start;                /* CM MATL node that emits to enf_cc_start */
  int nd;                       /* counter over nodes */
  
  emap      = CreateEmitMap(cm); 
  enf_start = -1;
  if(enf_cc_start > emap->clen)
    cm_Fail("ERROR --enfstart <n>, there's only %d columns, you chose column %d\n", 
	enf_cc_start, emap->clen);
  for(nd = 0; nd < cm->nodes; nd++)
    {
      if(emap->lpos[nd] == enf_cc_start) 
	{
	  if(cm->ndtype[nd] == MATL_nd)	      
	    {
	      enf_start = nd;
	      break;
	    }
	  else if(cm->ndtype[nd] == MATP_nd)	      
	    cm_Fail("ERROR --enfstart <n>, <n> must correspond to MATL modelled column\nbut %d is modelled by a MATP node.\n", enf_cc_start);
	}
      else if(emap->rpos[nd] == enf_cc_start)
	{
	  if(cm->ndtype[nd] == MATR_nd)	      
	    cm_Fail("ERROR --enfstart <n>, <n> must correspond to MATL modelled column\nbut %d is modelled by a MATR node.\n", enf_cc_start);
	  if(cm->ndtype[nd] == MATP_nd)	      
	    cm_Fail("ERROR --enfstart <n>, <n> must correspond to MATL modelled column\nbut %d is modelled by the right half of a MATP node.\n", enf_cc_start);
	}	      
    }
  if(enf_start == -1)
    cm_Fail("ERROR trying to determine the start node for the enforced subsequence.\n");
  FreeEmitMap(emap);
  return(enf_start);
}


/* Function:  CP9EnforceHackMatchScores()
 * Incept:    EPN, Fri Feb  9 11:06:31 2007
 *
 * Purpose:   Make all match emissions 0, except those enforce
 *            a specified subsequence (it's assumed the CP9 
 *            is already set up for this enforcement). 
 *
 * Args:      cp9           - the CP9 HMM 
 *            enf_start_pos - first posn of enforced subseq
 *            enf_end_pos   - last  posn of enforced subseq
 * Returns:   (void)
 */
void
CP9EnforceHackMatchScores(CP9_t *cp9, int enf_start_pos, int enf_end_pos)
{
  int k, x;
  for (k = 1; k < enf_start_pos; k++) /* M_0 is the begin state, it's silent */
    for (x = 0; x < MAXDEGEN; x++)
      cp9->msc[x][k] = 0.;
  for (k = enf_end_pos+1; k <= cp9->M; k++)
    for (x = 0; x < MAXDEGEN; x++)
      cp9->msc[x][k] = 0.;
}

/* Function: CPlan9SWConfigEnforce()
 * EPN, Fri Feb  9 05:47:37 2007
 * based on SRE's Plan7SWConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to hmmsw (Smith/Waterman) configuration.
 *           Same as CPlan9SWConfig but enforces a contiguous subset of
 *           nodes start at x, ending at y must be entered by forbidding
 *           local entries after x and local exits before y.
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over (enf_start) 
 *                    match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over (M-enf_end) match 
 *                    states.
 *           enf_start_pos - HMM node where enforced node subset begins         
 *           enf_end_pos   - HMM node where enforced node subset ends
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9SWConfigEnforce(CP9_t *hmm, float pentry, float pexit,
		      int enf_start_pos, int enf_end_pos)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */

  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* Configure entry.
   * To match CM, we enforce the only way out of the B state (M_0)
   * is through a local begin into a match state 
   */
  hmm->t[0][CTMI] = hmm->t[0][CTMD] = hmm->t[0][CTMEL] = 0.;
  hmm->begin[1] = 1. - pentry;
  for (k = 2; k <= enf_start_pos; k++)
    hmm->begin[k] = pentry / (float)(enf_start_pos-1);

  /* OLD WAY (more smith-waterman-like, less CM-like) EPN, Thu Jun 21 15:30:46 2007
     hmm->begin[1] = (1. - pentry) * (1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD])); 
  for (k = 2; k <= enf_start_pos; k++)
    hmm->begin[k] = (pentry * (1.- (hmm->t[0][CTMI] + hmm->t[0][CTMD]))) / (float)(enf_start_pos-1);
  for (k = (enf_start_pos+1); k <= hmm->M; k++)
    hmm->begin[k] = 0.;
  */
    
  /* Configure exit.
   * Don't touch hmm->end[hmm->M]
   */
  if(enf_end_pos == hmm->M) /* no local exit possible */
    basep = 0.0;
  else
    basep = pexit / (float) (hmm->M-enf_end_pos);

  for (k = 0; k < enf_end_pos; k++)
    hmm->end[k] = 0.;
  for (k = enf_end_pos; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) ((k-enf_end_pos)-1));
  CPlan9RenormalizeExits(hmm, 1);
  hmm->flags       &= ~CPLAN9_HASBITS;     /* reconfig invalidates log-odds scores */
  hmm->flags       |= CPLAN9_LOCAL_BEGIN; /* local begins now on */
  hmm->flags       |= CPLAN9_LOCAL_END;   /* local ends now on */
}
#endif

/* EPN, Mon Jan 21 10:20:50 2008
 * reorg of cmsearch, deemed some functions unnec: 
 * 
 * static int set_window(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
 * static int read_qdb_file(FILE *fp, CM_t *cm, int *dmin, int *dmax);
 * 
 * and others were rewritten:
 * static int set_searchinfo_OLD(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
 */

#if 0 

/* A CP9 filter work unit consists of a CM and an int (nseq).
 * The job is to emit nseq sequences with a score better than cutoff (rejecting
 * those that are worse), and then search those seqs with a CP9, returning the scores of the
 * best CP9 hit within each sequence.
 */
static int
process_cp9filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq)
{
  /*int status;*/
  cm_Fail("WRITE process_cp9filter_workunit()");
  return eslOK;
  
  /* ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_cp9filter_workunit\n", cfg->my_rank));
  return status;*/
}
#endif
#if 0
/* set_searchinfo_OLD()
 * Determine how many rounds of searching we will do (all rounds but last
 * round are filters), and set the relevant info in the SearchInfo_t <cm->si>
 * object, including cutoffs.
 */
static int
set_searchinfo_OLD(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;
  int n;
  int stype;
  int add_cyk_filter     = FALSE;
  int add_inside_filter  = FALSE;
  int add_viterbi_filter = FALSE;
  int add_forward_filter = FALSE;
  int search_opts;
  int use_hmmonly;
  int fthr_mode;
  int cutoff_type;
  float sc_cutoff = -1.;
  float e_cutoff = -1.;
  int  *dmin, *dmax; /* these become QDBs if we add a CM_FILTER */
  ScanMatrix_t *fsmx; 
  int safe_windowlen;
  float surv_fract;
  int cm_mode, cp9_mode;
  int cut_point;

  if(cm->si != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo(), cm->si is not NULL, shouldn't happen.\n");

  /* Create SearchInfo, specifying no filtering, we change the threshold below */
  CreateSearchInfo(cm, SCORE_CUTOFF, 0., -1.);
  if(cm->si == NULL) cm_Fail("set_searchinfo(), CreateSearchInfo() call failed.");
  SearchInfo_t *si = cm->si; 
  if(si->nrounds > 0) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_search_info(), si->nrounds (%d) > 0\n", si->nrounds);

  /*************************************************************************************
   * Filter related options:
   *
   * User can specify 0,1 or 2 rounds of filtering and cutoffs on the command line, 
   * 0 or 1 rounds can be CM  filters, with --fcyk or --finside, --fcmT,  and --fcmE (req's Gumbels)
   *                                        --finside specifies use inside, not CYK
   * 0 or 1 rounds can by HMM filters, with --fhmmviterbi or --fhmmforward, --fhmmT, and --fhmmE (req's Gumbels)
   *                                        --fhmmforward specifies use forward, not viterbi
   * Or user can specify that an HMM filter as described in the the CM file 
   * be used with option --fgiven. --fgiven is incompatible with --fhmmviterbi and --fhmmforward
   * but not with --fcyk and --finside
   *
   *************************************************************************************
   * Final round related options (after all filtering is complete):
   *
   * --cyk:        search with CM CYK (TRUE by default)
   * --inside:     search with CM inside 
   * -T:           CM bit score threshold
   * -E:           CM E-value threshold (requires Gumbel info in CM file)
   * --ga:         use Rfam gathering threshold (bit sc) from CM file
   * --tc:         use Rfam trusted cutoff      (bit sc) from CM file
   * --nc:         use Rfam noise cutoff        (bit sc) from CM file
   *
   * --viterbi: search with HMM viterbi
   * --forward: search with HMM forward
   * --hmmT:       bit score threshold for --viterbi or --forward
   * --hmmE:       E-value threshold (requires Gumbel info in CM file)
   *
   *************************************************************************************
   */
  
  /* First, set up cutoff for final round, this will be round 0, unless filter info was read from the CM file */
  n           = si->nrounds;
  stype       = si->stype[n];
  search_opts = si->search_opts[n];

  /* determine configuration of CM and CP9 based on cm->flags & cm->search_opts */
  CM2Gumbel_mode(cm, search_opts, &cm_mode, &cp9_mode); 

  use_hmmonly = ((search_opts & CM_SEARCH_HMMVITERBI) || (search_opts & CM_SEARCH_HMMFORWARD));
  if(! use_hmmonly) {
    if(stype != SEARCH_WITH_CM) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), search_opts for final round of search does not have HMMVITERBI or HMMFORWARD flags raised, but is not of type SEARCH_WITH_CM.");
    /* set up CM cutoff, either 0 or 1 of 6 options is enabled. 
     * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
     * We will NOT use this if --viterbi
     */
    if(esl_opt_IsDefault(go, "-E") && 
       esl_opt_IsDefault(go, "-T") && 
       esl_opt_IsDefault(go, "--ga") && 
       esl_opt_IsDefault(go, "--tc") && 
       esl_opt_IsDefault(go, "--nc")) { 
      /* Choose from, in order of priority:
       * 1. default CM E value if CM file has Gumbel stats
       * 3. default CM bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CM E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "-E");
	if((status = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "-T");
	e_cutoff    = -1.; /* invalid, we'll never use it */  
      }
    }
    else if(! esl_opt_IsDefault(go, "-E")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "-E requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "-E");
      if((status = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
    }
    else if(! esl_opt_IsDefault(go, "-T")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "-T");
      e_cutoff    = -1.; /* invalid, we'll never use it */  
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--greedy"))) ESL_FAIL(eslEINVAL, errbuf, "with -T <x> option, <x> can only be less than 0. if --greedy also enabled.");
    }
    else if(! esl_opt_IsDefault(go, "--ga")) {
      if(! (cm->flags & CMH_GA))
	ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = cm->ga;
      e_cutoff    = -1.; /* we'll never use it */
    }
    else if(! esl_opt_IsDefault(go, "--tc")) {
      if(! (cm->flags & CMH_TC))
	ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = cm->tc;
      e_cutoff    = -1.; /* we'll never use it */
    }
    else if(! esl_opt_IsDefault(go, "--nc")) {
      if(! (cm->flags & CMH_NC))
	ESL_FAIL(eslEINVAL, errbuf, "No NC noise cutoff in CM file, can't use --nc.");
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = cm->nc;
      e_cutoff    = -1.; /* we'll never use it */
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM cutoff selected. This shouldn't happen.");
  } /* end of if(! use_hmmonly) */
  else { 
    if(stype != SEARCH_WITH_HMM) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "search_opts for final round of search has HMMVITERBI or HMMFORWARD flags raised, but is not of type SEARCH_WITH_HMM.");
    /* Set up CP9 HMM cutoff, either 0 or 1 of 2 options is enabled 
     * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
     */
    if(esl_opt_IsDefault(go, "--hmmE") && 
       esl_opt_IsDefault(go, "--hmmT")) {
      /* Choose from, in order of priority:
       * 1. default CP9 E value if CM file has Gumbel stats
       * 2. default CP9 bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CP9 E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "--hmmE"); 
	if((status = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
	e_cutoff    = -1; /* we'll never use it */
      }
    }
    else if(! esl_opt_IsDefault(go, "--hmmE")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--hmmE requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "--hmmE");
      if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
    }
    else if(! esl_opt_IsDefault(go, "--hmmT")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
      e_cutoff    = -1.; /* we'll never use this */
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--hmmgreedy"))) ESL_FAIL(eslEINVAL, errbuf, "with --hmmT <x> option, <x> can only be less than 0. if --hmmgreedy also enabled.");
    }
  }
  /* update the search info, which holds the thresholds */
  UpdateSearchInfoCutoff(cm, cm->si->nrounds, cutoff_type, sc_cutoff, e_cutoff);   
  ValidateSearchInfo(cm, cm->si);
  /* DumpSearchInfo(cm->si); */
  /* done with threshold for final round */

  /* Set up the filters and their thresholds 
   * 1. add a CM  filter, if necessary
   * 2. add a HMM filter, if necessary
   */

  /* CM filter */
  add_cyk_filter    = esl_opt_GetBoolean(go, "--fcyk");
  add_inside_filter = esl_opt_GetBoolean(go, "--finside");
  ESL_DASSERT1((!(add_cyk_filter && add_inside_filter))); /* should be enforced by getopts */
  if(add_cyk_filter || add_inside_filter) { /* determine thresholds for filters */
    /* set up CM filter cutoff, either 0 or 1 of 2 options is enabled. 
     * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
     */
    if(esl_opt_IsDefault(go, "--fE") && 
       esl_opt_IsDefault(go, "--fT")) {
      /* Choose from, in order of priority:
       * 1. default CM filter E value if CM file has Gumbel stats
       * 3. default CM filter bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CM E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "--fE");
	if((status  = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "--fT");
	e_cutoff    = -1.; /* invalid, we'll never use it */  
      }
    }
    else if(! esl_opt_IsDefault(go, "--fE")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--fE requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "--fE");
      if((status  = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
    }
    else if(! esl_opt_IsDefault(go, "--fT")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--fT");
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--fgreedy"))) ESL_FAIL(eslEINVAL, errbuf, "with --fT <x> option, <x> can only be less than 0. if --fgreedy also enabled.");
      e_cutoff    = -1.; /* we'll never use it */
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM filter cutoff selected. This shouldn't happen.");
    
    /* build the ScanMatrix_t for this round, requires calcing dmin, dmax */
    safe_windowlen = cm->W * 2;
    while(!(BandCalculationEngine(cm, safe_windowlen, esl_opt_GetReal(go, "--fbeta"), FALSE, &dmin, &dmax, NULL, NULL))) {
      free(dmin);
      free(dmax);
      dmin = NULL;
      dmax = NULL;
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), band calculation safe_windowlen big: %d\n", safe_windowlen);
    }
    fsmx = cm_CreateScanMatrix(cm, dmax[0], dmin, dmax, esl_opt_GetReal(go, "--fbeta"), TRUE, add_cyk_filter, add_inside_filter);
    /* add the filter */
    AddFilterToSearchInfo(cm, add_cyk_filter, add_inside_filter, FALSE, FALSE, FALSE, fsmx, NULL, cutoff_type, sc_cutoff, e_cutoff);
    ValidateSearchInfo(cm, cm->si);
    /* DumpSearchInfo(cm->si); */
  }
  else if (! esl_opt_IsDefault(go, "--fbeta")) ESL_FAIL(eslEINCOMPAT, errbuf, "--fbeta has an effect with --fcyk or --finside");

  /* HMM filter */
  /* if --fgiven was enabled, --fhmmviterbi, --fhmmforward could NOT have been selected, 
   * so we won't enter any of the loops below.
   */
  add_viterbi_filter = esl_opt_GetBoolean(go, "--fhmmviterbi");
  add_forward_filter = esl_opt_GetBoolean(go, "--fhmmforward");
  ESL_DASSERT1((!(add_viterbi_filter && add_forward_filter))); /* should be enforced by getopts */

  if(esl_opt_GetBoolean(go, "--fgiven")) {
    /* determine filter threshold mode, the mode of final stage of searching, either FTHR_CM_LC,
     * FTHR_CM_LI, FTHR_CM_GC, FTHR_CM_GI, (can't be an HMM mode b/c getopts enforces --fgiven incompatible with
     * --viterbi and --forward). 
     */

    if((status = CM2FthrMode(cm, errbuf, cm->search_opts, &fthr_mode)) != eslOK) return status;
    HMMFilterInfo_t *hfi_ptr = cm->stats->hfiA[fthr_mode]; /* for convenience */
    if(!(cm->flags & CMH_FILTER_STATS))   ESL_FAIL(eslEINCOMPAT, errbuf,      "set_searchinfo(), --fgiven enabled, but cm's CMH_FILTER_STATS flag is down.");
    if(hfi_ptr->is_valid == FALSE)        ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), --fgiven enabled, cm's CMH_FILTER_STATS is raised, but best filter info for fthr_mode %d is invalid.", fthr_mode);
    ESL_DASSERT1((!add_viterbi_filter)); /* --fhmmviterbi is incompatible with --fgiven, enforced by getopts */
    ESL_DASSERT1((!add_forward_filter)); /* --fhmmforward is incompatible with --fgiven, enforced by getopts */

    cutoff_type = E_CUTOFF;
    
    /* determine the appropriate filter cut point <cut_point> to use, for details
     * see comments in function: searchinfo.c:GetHMMFilterFwdECutGivenCME() for details */
    if(cm->si->cutoff_type[cm->si->nrounds] == SCORE_CUTOFF) { 
      if((status = GetHMMFilterFwdECutGivenCMBitScore(hfi_ptr, errbuf, cm->si->sc_cutoff[cm->si->nrounds], cfg->dbsize, &cut_point, cm, cm_mode)) != eslOK) return status; 
      /* cm->si->sc_cutoff[cm->si->nrounds] is bit score cutoff used for the final round of searching (when we'll be done filtering) */
    }
    else { 
      if((status = GetHMMFilterFwdECutGivenCME(hfi_ptr, errbuf, cm->si->e_cutoff[cm->si->nrounds], cfg->dbsize, &cut_point)) != eslOK) return status; 
      /* cm->si->e_cutoff[cm->si->nrounds] is E cutoff used for the final round of searching (when we'll be done filtering) */
    }
    if(cut_point != -1) { /* it's worth it to filter */
      e_cutoff = hfi_ptr->fwd_E_cut[cut_point] * ((double) cfg->dbsize / (double) hfi_ptr->dbsize); 
      /* check if --hmmEmax applies */
      if(! (esl_opt_IsDefault(go, "--hmmEmax"))) {
	e_cutoff = ESL_MIN(e_cutoff, esl_opt_GetReal(go, "--hmmEmax"));
      }
      if((status  = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; /* note: use cm_mode, not fthr_mode */
      add_forward_filter = TRUE;
      /* TEMPORARY! */
      /* Predict survival fraction from filter based on E-value, assume average hit length is cfg->avg_hit_len (from QD band calc) */
      surv_fract = (e_cutoff * ((2. * cm->W) - cfg->avg_hit_len)) / ((double) cfg->dbsize); 
      printf("\n\nsurv_fract: %f\n\n\n", surv_fract);
    }
    else { /* it's not worth it to filter, our HMM filter cutoff would be so low, 
	    * letting so much of the db survive, the filter is a waste of time */
      add_forward_filter = FALSE;
      printf("cut_point -1, always_better FALSE\n");
    }
  }
  else if(add_viterbi_filter || add_forward_filter) { /* determine thresholds for filters */
    /* Set up HMM cutoff, either 0 or 1 of 3 options is enabled */
    ESL_DASSERT1((! use_hmmonly)); /* should be enforced by getopts */
    if(esl_opt_IsDefault(go, "--hmmcalcthr") && 
       esl_opt_IsDefault(go, "--hmmE") && 
       esl_opt_IsDefault(go, "--hmmT")) {
      /* Choose from, in order of priority:
       * 1. default CP9 E value if CM file has Gumbel stats
       * 2. default CP9 bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CM E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "--hmmE");
	if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
	e_cutoff    = -1.; /* we'll never use this */
      }
    }
    else if(! esl_opt_IsDefault(go, "--hmmcalcthr")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--hmmcalcthr requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      /* this gets overwritten later after threshold is calculated */
      e_cutoff = esl_opt_GetReal(go, "--hmmE");
      if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
    }
    else if(! esl_opt_IsDefault(go, "--hmmE")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--hmmE requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "--hmmE");
      if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
    }
    else if(! esl_opt_IsDefault(go, "--hmmT")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
      e_cutoff    = -1.; /* we'll never use this */
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--hmmgreedy"))) ESL_FAIL(eslEINVAL, errbuf, "with --hmmT <x> option, <x> can only be less than 0. if --hmmgreedy also enabled.");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No HMM filter cutoff selected. This shouldn't happen.");
  }
  if(add_viterbi_filter || add_forward_filter) {
    /* add the filter */
    AddFilterToSearchInfo(cm, FALSE, FALSE, add_viterbi_filter, add_forward_filter, FALSE, NULL, NULL, cutoff_type, sc_cutoff, e_cutoff);
    ValidateSearchInfo(cm, cm->si);
    /*DumpSearchInfo(cm->si); */
  }

  return eslOK;
}
#endif
#if 0
/* set_window()
 * Set cm->W, the window size for scanning.
 *
 * 1. cm->W is set to dmax[0] if --noqdb, --viterbi, or --forward enabled after calc'ing QDBs 
 *    with beta == cm->beta of (esl_opt_GetReal(go, "--beta")) if that was enabled. 
 * 3. else cm->W was set to cm->dmax[0] in ConfigCM()'s call to ConfigQDB(), 
 *    which is what it should be.
 */
static int
set_window(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int use_hmmonly;
  int do_float;
  int do_int;

  use_hmmonly = (esl_opt_GetBoolean(go, "--viterbi") || esl_opt_GetBoolean(go, "--forward")) ? TRUE : FALSE;

  if(esl_opt_GetBoolean(go, "--noqdb") || use_hmmonly) {
    if(cm->dmin != NULL || cm->dmax != NULL) 
      ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--viterbi, --forward or --noqdb enabled, but cm->dmin and cm->dmax non-null. This shouldn't happen.");
    int *dmin;
    int *dmax;
    int safe_windowlen = cm->clen * 2;
    while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta, 0, &(dmin), &(dmax), NULL, NULL)))
      {
	free(dmin);
	free(dmax);
	safe_windowlen *= 2;
	if(safe_windowlen > (cm->clen * 1000))
	  ESL_FAIL(eslEINVAL, errbuf, "ERROR in set_window, safe_windowlen big: %d\n", safe_windowlen);
      }
    cm->W = dmax[0];
    free(dmin);
    free(dmax);
    CMLogoddsify(cm); /* QDB calculation invalidates log odds scores */
  }

  /* Setup ScanMatrix for CYK/Inside scanning functions, we can't 
   * do it in initialize_cm(), b/c it's W dependent; W was just set.
   * We don't need it if we're only using an HMM though.
   */
  if(use_hmmonly) cm->smx = NULL;
  else { 
    do_float = TRUE;
    do_int   = FALSE;
    if(cm->search_opts & CM_SEARCH_INSIDE) { do_float = FALSE; do_int = TRUE; }
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
    if(cm->smx == NULL) cm_Fail("set_window(), use_hmmonly is FALSE, CreateScanMatrixForCM() call failed, mx is NULL.");
  }

  return eslOK;
}
#endif
#if 0

/* read_qdb_file()
 * Read QDBs from a file outputted from cmbuild. Only useful for testing/debugging,
 */
static int  
read_qdb_file(FILE *fp, CM_t *cm, int *dmin, int *dmax)
{
  int     status;
  char   *buf;
  int     n;			/* length of buf */
  char   *s;
  int     M;			/* number of states in model */
  int     v;		        /* counter for states */
  char   *tok;
  int     toklen;
  int     read_v;

  /* format of QDB file: 
   * line  1        :<cm->M>
   * lines 2 -> M+1 :<v> <dmin> <dmax> */

  buf = NULL;
  n   = 0;
  if (feof(fp) || (status = esl_fgets(&buf, &n, fp)) != eslOK) goto ERROR;

  s   = buf;
  if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;
  if (! is_integer(tok))                                    goto ERROR;
  M = atoi(tok);
  if(M != cm->M) goto ERROR;

  v = 0;
  while ((status = esl_fgets(&buf, &n, fp)) == eslOK) 
    {
      if (strncmp(buf, "//", 2) == 0) 
	break;
      s   = buf;
      if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;      
      if (! is_integer(tok)) { status = eslEINVAL;                    goto ERROR; }
      read_v = atoi(tok);
      if(v != read_v) goto ERROR;

      if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;      
      if (! is_integer(tok)) { status = eslEINVAL;                    goto ERROR; }
      dmin[v] = atoi(tok);

      if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;      
      if (! is_integer(tok)) {                                        goto ERROR; }
      dmax[v] = atoi(tok);

      v++;
    }
  if(v != M) { status = eslEINVAL; goto ERROR; }
  if(status != eslOK) goto ERROR;

  if (buf != NULL) free(buf);
  return eslOK;

 ERROR:
  if (cm != NULL)  FreeCM(cm);
  if (buf != NULL) free(buf);
  return status;
}
#endif
#if
//extern float      MinScCutoff (CM_t *cm, SearchInfo_t *si, int n);
/*
 * Function: MinScCutoff
 * Date:     EPN, Mon May  7 17:36:56 2007
 *
 * Purpose:  Return the minimum bit score cutoff for CM
 *           for round n in SearchInfo_t si.
 *           Trivial if si->cutoff_type[n] == SCORE_CUTOFF,
 *           but if E_CUTOFF return minimal bit score across 
 *           all partitions for the E cutoff in the 
 *           appropriate search algorithm.
 *
 */
float MinScCutoff (CM_t *cm, SearchInfo_t *si, int n)
{
  float E, low_sc, sc;
  int cm_mode, cp9_mode, gum_mode;
  int p; 

  /* contract check */
  if(si == NULL)      cm_Fail("MinCMScCutoff(), si == NULL.\n");
  if(n > si->nrounds) cm_Fail("MinCMScCutoff(), n (%d) > si->nrounds\n", n, si->nrounds);

  if(si->cutoff_type[n] == SCORE_CUTOFF) return si->sc_cutoff[n];
  
  /* if we get here, cutoff_type is E_CUTOFF we better have stats */
  ESL_DASSERT1((si->cutoff_type[n] == E_CUTOFF));
  if(!(cm->flags & CMH_GUMBEL_STATS)) cm_Fail("ERROR in MinScCutoff, cutoff type E value, but no stats.\n");

  /* Determine appropriate Gumbel mode */
  CM2Gumbel_mode(cm, si->search_opts[n], &cm_mode, &cp9_mode);
  E = si->e_cutoff[n];

  if(si->stype[n] == SEARCH_WITH_HMM) {
    ESL_DASSERT1(((si->search_opts[n] & CM_SEARCH_HMMVITERBI) || (si->search_opts[n] & CM_SEARCH_HMMFORWARD)));
    gum_mode = cp9_mode; 
  }
  else if (si->stype[n] == SEARCH_WITH_CM) {
    ESL_DASSERT1((! ((si->search_opts[n] & CM_SEARCH_HMMVITERBI) || (si->search_opts[n] & CM_SEARCH_HMMFORWARD))));
    gum_mode = cm_mode; 
  }
  else cm_Fail("MinScCutoff(), asking for E-value cutoff for SEARCH_WITH_HYBRID search round.\n");

  low_sc = cm->stats->gumAA[cm_mode][0]->mu - 
    (log(E) / cm->stats->gumAA[gum_mode][0]->lambda);
  for (p = 1; p < cm->stats->np; p++) {
    sc = cm->stats->gumAA[gum_mode][p]->mu - 
      (log(E) / cm->stats->gumAA[gum_mode][p]->lambda);
    if (sc < low_sc) low_sc = sc;
  }
  return (low_sc);
}
#endif


/* EPN, Mon Feb  4 15:55:43 2008
 * Working in ~/notebook/7_0127_inf_hmm2ij_simple/
 * cleaning up old versions of functions that were developed
 * during that notebook directory's duration of work.
 */
/*********************************************************************
 * Function: hmm2ij_newer()
 * 
 *
 * Args:
 * cm               the cm
 * i0               first position of seq
 * j0               last position of seq
 * int *imin        imin[v] = first position in band on i for state v
 * int *imax        imax[v] = last position in band on i for state v
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 */
int
hmm2ij_newer(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int doing_search, int debug_level)
{

  int status;
  int v;

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
  
  int nd;           /* counter over CM nodes. */
  int y;   /* counters over children states */
  /* Contract checks */

  if (cp9b == NULL)                                                                   ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), cp9b is NULL.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED))) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down, exactly 1 must be up.\n");
  if(i0 < 1) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), i0 < 1: %d\n", i0);
  if(j0 < 1) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), j0 < 1: %d\n", j0);
  if(j0 < i0) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_HMM2ijBands(), i0 (%d) < j0 (%d)\n", i0, j0);

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

#if 0 
  /* Initialize all *min bands to INT_MAX and *max bands to -INT_MAX */
  esl_vec_ISet(imin, cm->M,  INT_MAX);
  esl_vec_ISet(imax, cm->M, -INT_MAX);
  esl_vec_ISet(jmin, cm->M,  INT_MAX);
  esl_vec_ISet(jmax, cm->M, -INT_MAX);
#endif
  /* Initialize all *min bands to INT_MAX and *max bands to -INT_MAX */
  esl_vec_ISet(imin, cm->M, -1);
  esl_vec_ISet(imax, cm->M, -1);
  esl_vec_ISet(jmin, cm->M, -1);
  esl_vec_ISet(jmax, cm->M, -1);

  /* TEMPORARY to deal with incorrect handling of delete off-by-one in cp9_FB2HMMBandsWithSums() */
  int k;
  for(k = 0; k <= cm->cp9->M; k++) { 
    if(pn_min_d[k] != -1) pn_min_d[k]--;
    if(pn_min_d[k] != -1) pn_max_d[k]--;
  }

  ////cp9_DebugPrintHMMBands(stdout, j0, cp9b, cm->tau, 1);

  /* Step 2: Traverse nodes from cpos = 1..clen using a stack and filling
   *         in i and j bands (imin, imax, jmin, jmax) for all states that
   *         don't map to HMM states.
   */
  ESL_STACK   *pda;
  ESL_STACK   *pda2;
  int          on_right;
  int          abv_ss_imin, abv_ss_imax, abv_ss_jmin, abv_ss_jmax;
  int          blw_ss_imin, blw_ss_imax, blw_ss_jmin, blw_ss_jmax;
  int          ss_imin, ss_imax, ss_jmin, ss_jmax;
  int          sd;
  int          w;
  int          just_saw_begl_on_right = FALSE;
#if 0
  abv_ss_imin = blw_ss_imin = i0;
  abv_ss_imax = blw_ss_imax = j0+1;
  abv_ss_jmin = blw_ss_jmin = i0-1; 
  abv_ss_jmax = blw_ss_jmax = j0;
#endif
  abv_ss_imin = abv_ss_jmin = blw_ss_imin = blw_ss_jmin =  INT_MAX;
  abv_ss_imax = abv_ss_jmax = blw_ss_imax = blw_ss_jmax = -INT_MAX;
  ///if(! doing_search) { abv_ss_imin = i0; abv_ss_jmax = j0; } /* this will make ROOT_S necessarily have imin[0] = i0, jmax[0] = j0, which
  ///* is mandatory if we're aligning, b/c we have to align the full seq */
  
  nd   = 0;
  if ((pda  = esl_stack_ICreate()) == NULL) goto ERROR;
  if ((pda2 = esl_stack_ICreate()) == NULL) goto ERROR;
  if ((status = esl_stack_IPush(pda, 0)) != eslOK) goto ERROR;		/* 0 = left side. 1 would = right side. */
  if ((status = esl_stack_IPush(pda, nd)) != eslOK) goto ERROR;
  while (esl_stack_IPop(pda, &nd) != eslEOD)
    {
      esl_stack_IPop(pda, &on_right);

      /* begin temporary block */
      if (just_saw_begl_on_right) { if((cm->ndtype[nd] != BEGR_nd) || on_right)     cm_Fail("you don't get it (1)."); }
      else                        { if((cm->ndtype[nd] == BEGR_nd) && (! on_right)) cm_Fail("you don't get it (2)."); }
      just_saw_begl_on_right = FALSE;
      /* end temporary block */

      if (on_right) 
	{
	  ////printf("on right nd: %d\n", nd);
	  if(cm->ndtype[nd] == BIF_nd) { /* special case, set i bands based on left child, j bands based on right child */
	    esl_stack_IPop(pda2, &blw_ss_imax);
	    esl_stack_IPop(pda2, &blw_ss_imin);

	    v = cm->nodemap[nd];
	    w = cm->cfirst[v]; /* BEGL_S */
	    y = cm->cnum[v];   /* BEGR_S */
	    imin[v] = imin[w];
	    imax[v] = imax[w];
	    jmin[v] = jmin[y];
	    jmax[v] = jmax[y];
	  }
	  else { 
	    HMM2ijBandsForNode(cm, cp9b, cp9map, nd, i0, &ss_imin, &ss_imax, &ss_jmin, &ss_jmax);
	    
	    if(cm->ndtype[nd] != BEGR_nd) {  /* HACK! justify this, reason is that BEGR only has an insert, and we can set
					      * that insert in HMM2ijBnadsForNode if it has an explicit band from the HMM,
					      * but that shouldn't impact our implicit band on BEGR...
					      */
	      if(ss_imin != -1) blw_ss_imin = ss_imin;
	      if(ss_imax != -1) blw_ss_imax = ss_imax;
	      if(ss_jmin != -1) blw_ss_jmin = ss_jmin;
	      if(ss_jmax != -1) blw_ss_jmax = ss_jmax;
	    }
	    /* pop off {s,i}_{i,j}{min,max} if they're not -1, update cur_{s,i}_{i,j}{min, max} */
	    esl_stack_IPop(pda, &abv_ss_jmax);
	    esl_stack_IPop(pda, &abv_ss_jmin);
	    esl_stack_IPop(pda, &abv_ss_imax);
	    esl_stack_IPop(pda, &abv_ss_imin);
	    ss_imin = ESL_MIN(abv_ss_imin, blw_ss_imin);
	    ss_imax = ESL_MAX(abv_ss_imax, blw_ss_imax);
	    ss_jmin = ESL_MIN(abv_ss_jmin, blw_ss_jmin);
	    ss_jmax = ESL_MAX(abv_ss_jmax, blw_ss_jmax);
	    if(ss_imin ==  INT_MAX) ss_imin = i0;
	    if(ss_imax == -INT_MAX) ss_imax = j0+1;
	    if(ss_jmin ==  INT_MAX) ss_jmin = i0-1;
	    if(ss_jmax == -INT_MAX) ss_jmax = j0;

	    ////printf("nd: %4d abv  (%11d %11d  %11d %11d)\n", nd, abv_ss_imin, abv_ss_imax, abv_ss_jmin, abv_ss_jmax);
	    ////printf("nd: %4d blw  (%11d %11d  %11d %11d)\n", nd, blw_ss_imin, blw_ss_imax, blw_ss_jmin, blw_ss_jmax);
	    ////printf("nd: %4d bst  (%11d %11d  %11d %11d)\n", nd, ss_imin, ss_imax, ss_jmin, ss_jmax);

	    if(cm->ndtype[nd] == BEGL_nd) { /* reset blw_ss_{i,j}{min,max} to their initial values for next stem */
	      //// printf("HEYA resetting blw_ss_imin: %d blw_ss_imax: %d blw_ss_jmin: %d blw_ss_jmax: %d\n", blw_ss_imin, blw_ss_imax, blw_ss_jmin, blw_ss_jmax);
	      if ((status = esl_stack_IPush(pda2, blw_ss_imin)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda2, blw_ss_imax)) != eslOK) goto ERROR;
	      abv_ss_imin = abv_ss_jmin = blw_ss_imin = blw_ss_jmin =  INT_MAX;
	      abv_ss_imax = abv_ss_jmax = blw_ss_imax = blw_ss_jmax = -INT_MAX;
	      ///blw_ss_imin = blw_ss_jmin =  INT_MAX;
	      ///blw_ss_imax = blw_ss_jmax = -INT_MAX;
	      just_saw_begl_on_right = TRUE;
	    }
	  }	    
	  /* set bands for states in node nd */
	  /* the split set states */
	  for(v = cm->nodemap[nd]; v < (cm->nodemap[nd] + TotalStatesInNode(cm->ndtype[nd])); v++) { 
	    sd = StateDelta(cm->sttype[v]);
#if 0 
	    if(imin[v] == -1) { 
	      printf("!i v: %d setting imin[v] as MIN(ss_imin: %d (jmin[v]-sd+1): %d)\n", v, ss_imin, jmin[v]-sd+1);
	      imin[v] = (jmin[v] == -1) ? ss_imin : ESL_MIN(ss_imin, (jmin[v] - sd + 1)); /* if jmin != INT_MAX, we set jmin explicitly based on an HMM band and we enforce that jmin[v] - imin[v] +1 >= sd */
	    }
	    if(imax[v] == -1) imax[v] = ss_imax;
	    if(jmin[v] == -1) { 
	      printf("!j v: %d setting jmin[v] as MIN(ss_jmin: %d (imin[v]+sd-1): %d)\n", v, ss_jmin, imin[v]+sd-1);
	      jmin[v] = ESL_MAX(ss_jmin, (imin[v] + sd - 1)); /* enforces that jmin[v] - imin[v] +1 >= sd, if jmin[v] - imin[v] + 1 < sd, there is no valid i for which j == jmin[v] is valid (d will be < sd) */
	    }
	    if(jmax[v] == -1) jmax[v] = ss_jmax;
#endif
#if  1
	    if(imin[v] == -1) imin[v] = ss_imin;
	    if(imax[v] == -1) imax[v] = ss_imax;
	    if(jmin[v] == -1) jmin[v] = ESL_MAX(ss_jmin, i0+sd-1); /* j must be at least i0+sd-1, that is i0 for singlet emitters, i0+1 for MP states */
	    if(jmax[v] == -1) jmax[v] = ESL_MAX(ss_jmax, i0+sd-1); /* j must be at least i0+sd-1, that is i0 for singlet emitters, i0+1 for MP states */
#endif
	  }
	}
      else
	{
	  if (cm->ndtype[nd] == BIF_nd)
	    {
	      /* BIF is special, we infer bands from left and right child, no pushing nec */
  			    /* push the BIF back on for its right side  */
	      if ((status = esl_stack_IPush(pda, 1)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda, nd)) != eslOK) goto ERROR;
                            /* push node index for right child */
	      if ((status = esl_stack_IPush(pda, 0)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda, cm->ndidx[cm->cnum[cm->nodemap[nd]]])) != eslOK) goto ERROR;   
                            /* push node index for left child */
	      if ((status = esl_stack_IPush(pda, 0)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda, cm->ndidx[cm->cfirst[cm->nodemap[nd]]])) != eslOK) goto ERROR; 
	    }
	  else
	    {
	      HMM2ijBandsForNode(cm, cp9b, cp9map, nd, i0, &ss_imin, &ss_imax, &ss_jmin, &ss_jmax);
	      
	      if(ss_imin != -1) abv_ss_imin = ss_imin;
	      if(ss_imax != -1) abv_ss_imax = ss_imax;
	      if(ss_jmin != -1) abv_ss_jmin = ss_jmin;
	      if(ss_jmax != -1) abv_ss_jmax = ss_jmax;
	      /* push bands to stack */
	      if ((status = esl_stack_IPush(pda, abv_ss_imin)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda, abv_ss_imax)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda, abv_ss_jmin)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda, abv_ss_jmax)) != eslOK) goto ERROR;
				/* push the node back on for right side */
	      if ((status = esl_stack_IPush(pda, 1)) != eslOK) goto ERROR;
	      if ((status = esl_stack_IPush(pda, nd)) != eslOK) goto ERROR;
                 	      /* push child node on */
	      if (cm->ndtype[nd] != END_nd) {
		if ((status = esl_stack_IPush(pda, 0)) != eslOK) goto ERROR;
		if ((status = esl_stack_IPush(pda, nd+1)) != eslOK) goto ERROR;
	      }
	    }
	}
    }
  if(! doing_search) { /* if we're aligning the full seq must be aligned at the root state */
    imin[0] = imin[1] = i0; /* first residue must be in subtree of ROOT_S and IFF ROOT_IL is used, of ROOT_IL also */
    jmax[0] = jmax[1] = jmax[2] = j0; /* final residue must be in subtree of ROOT_S and IFF ROOT_IL is used and/or ROOT_IR is used, of them also */
  }

  for (v = 0; v < cm->M; v++) { 
    /* TEMPORARY */ if(StateIsDetached(cm, v)) imin[v] = imax[v] = jmin[v] = jmax[v] = i0; 
  }
  ////debug_print_ij_bands(cm); 


  /* AT THIS POINT ALL BANDS SHOULD BE WITHIN VALID RANGE! */
  /* EVENTUALLY WHEN --simple IS DEFAULT, PUT THIS IN cp9_validateBands() */
  //ESL_DPRINTF1(("i0: %d j0: %d\n", i0, j0));
  ////printf("i0: %d j0: %d\n", i0, j0);
  for (v = 0; v < cm->M; v++) { 
    sd = StateDelta(cm->sttype[v]);

    //ESL_DPRINTF1(("v: %4d in: %4d ix: %4d  jn: %4d jx: %4d\n", v, imin[v], imax[v], jmin[v], jmax[v]));
    ////printf("v: %4d in: %4d ix: %4d  jn: %4d jx: %4d\n", v, imin[v], imax[v], jmin[v], jmax[v]);
    ESL_DASSERT1((imin[v] >= (i0)));    /* min allowed i is i0 */
    ESL_DASSERT1((imax[v] >= (i0)));    /* min allowed i is i0 */
    ESL_DASSERT1((imin[v] <= (j0+1)));  /* max allowed i is j0+1 */
    ESL_DASSERT1((imax[v] <= (j0+1)));  /* max allowed i is j0+1 */
    ESL_DASSERT1((jmin[v] >= (i0-1)));  /* min allowed j is i0-1 */
    ESL_DASSERT1((jmax[v] >= (i0-1)));  /* min allowed j is i0-1 */
    ESL_DASSERT1((jmin[v] <= (j0)));    /* max allowed j is j0 */
    ESL_DASSERT1((jmax[v] <= (j0)));    /* max allowed j is j0 */

    ESL_DASSERT1((imax[v] >= imin[v])); /* imax[v] >= imin[v] */
    ESL_DASSERT1((jmax[v] >= jmin[v])); /* jmax[v] >= jmin[v] */

    ESL_DASSERT1(((jmax[v]-imin[v]+1) <= (j0-i0+1))); /* largest  possible d is max j - max i + 1 = j0-i0+1 */

    assert(imin[v] >= (i0));    /* min allowed i is i0 */
    assert(imax[v] >= (i0));    /* min allowed i is i0 */
    assert(imin[v] <= (j0+1));  /* max allowed i is j0+1 */
    assert(imax[v] <= (j0+1));  /* max allowed i is j0+1 */
    assert(jmin[v] >= (i0-1));  /* min allowed j is i0-1 */
    assert(jmax[v] >= (i0-1));  /* min allowed j is i0-1 */
    assert(jmin[v] <= (j0));    /* max allowed j is j0 */
    assert(jmax[v] <= (j0));    /* max allowed j is j0 */

    assert(imax[v] >= imin[v]); /* imax[v] >= imin[v] */
    assert(jmax[v] >= jmin[v]); /* jmax[v] >= jmin[v] */

    assert((jmax[v]-imin[v]+1) <= (j0-i0+1)); /* largest  possible d is max j - min i + 1 = j0-i0+1 */
  }

  /* Step 3: Enforce safe transitions.
   *         Goal is to enforce there's at least 1 valid path through the model.
   */
  int child_imin, child_imax, child_jmin, child_jmax;
  int *v_is_r,  *r_imin, *r_imax, *r_jmin, *r_jmax; /* r = reachable */
  int *nd_is_r, *nd_r_imin, *nd_r_imax, *nd_r_jmin, *nd_r_jmax; /* r = reachable */
  int sdl, sdr;
  int y_nd, w_nd;
  ESL_ALLOC(v_is_r,   sizeof(int) * cm->M);
  ESL_ALLOC(r_imin, sizeof(int) * cm->M);
  ESL_ALLOC(r_imax, sizeof(int) * cm->M);
  ESL_ALLOC(r_jmin, sizeof(int) * cm->M);
  ESL_ALLOC(r_jmax, sizeof(int) * cm->M);
  ESL_ALLOC(nd_is_r,   sizeof(int) * cm->nodes);
  ESL_ALLOC(nd_r_imin, sizeof(int) * cm->nodes);
  ESL_ALLOC(nd_r_imax, sizeof(int) * cm->nodes);
  ESL_ALLOC(nd_r_jmin, sizeof(int) * cm->nodes);
  ESL_ALLOC(nd_r_jmax, sizeof(int) * cm->nodes);

  esl_vec_ISet(v_is_r, cm->M, FALSE);
  esl_vec_ISet(nd_is_r, cm->nodes, FALSE);

#if 0
  for (v = 0; v < cm->M; v++) { 
    r_imin[v] = imin[v];
    r_imax[v] = imax[v];
    r_jmin[v] = jmin[v];
    r_jmax[v] = jmax[v];
  }

  for (nd = 0; nd < cm->nodes; nd++) { 
    v0 = cm->nodemap[nd];
    nd_r_imin[nd] = imin[v0];
    nd_r_imax[nd] = imax[v0];
    nd_r_jmin[nd] = jmin[v0];
    nd_r_jmax[nd] = jmax[v0];
  }
#endif
  for (v = 0; v < cm->M; v++) { 
    r_imin[v] =  INT_MAX;
    r_imax[v] = -INT_MAX;
    r_jmin[v] =  INT_MAX;
    r_jmax[v] = -INT_MAX;
  }

  for (nd = 0; nd < cm->nodes; nd++) { 
    nd_r_imin[nd] =  INT_MAX;
    nd_r_imax[nd] = -INT_MAX;
    nd_r_jmin[nd] =  INT_MAX;
    nd_r_jmax[nd] = -INT_MAX;
  }

  ////printf("safe parse check\n");
  nd_is_r[0] = TRUE;
  v_is_r[0] = TRUE;
  r_imin[0] = nd_r_imin[0] = imin[0];
  r_imax[0] = nd_r_imax[0] = imax[0];
  r_jmin[0] = nd_r_jmin[0] = jmin[0];
  r_jmax[0] = nd_r_jmax[0] = jmax[0];
  
  for (nd = 0; nd < cm->nodes; nd++) { 
    for (v = cm->nodemap[nd]; v < (cm->nodemap[nd] + TotalStatesInNode(cm->ndtype[nd])); v++) { 
      if(StateIsDetached(cm, v)) { ; } 
      else if(cm->sttype[v] == E_st) { 
	if((r_imin[v] <= r_imax[v] && r_jmin[v] <= r_jmax[v]) && ((r_imax[v] - r_jmin[v] - 1) >= 0)) { 
	  /* END state v is reachable for some i, j such that j-i+1 = d = 0 (which is required for E states) */
	  v_is_r[v] = TRUE;
	  nd_is_r[nd] = TRUE;
	}
	else { /* no possible i, j such that j-i+1 = d = 0, this means E is unreachable (d must be 0 for end states */	  
	  v_is_r[v] = FALSE;
	  nd_is_r[nd] = FALSE;
	  //printf("END_E nd: %d v: %d unreachable!\n", nd, v);
	  //return eslOK;
	}
      }	      
      else if (cm->sttype[v] == B_st) { 
	/* same loop as if v != B_st, (the else case below) but we know sdl = sdr = 0, and we have two children BEGL_S and BEGR_S */
	if((r_imin[v] <= r_imax[v] && r_jmin[v] <= r_jmax[v]) && ((r_jmax[v] - r_imin[v] + 1) >= sd)) { 
	  /* v is reachable for some i, j */
	  v_is_r[v] = TRUE; 
	  nd_is_r[nd] = TRUE;
	  w = cm->cfirst[v]; /* BEGL_S */
	  y = cm->cnum[v];   /* BEGR_S */

	  r_imin[w] = ESL_MAX(imin[w], imin[v]);
	  r_imax[w] = ESL_MIN(imax[w], imax[v]);
	  r_jmin[w] = jmin[w];
	  r_jmax[w] = jmax[w];
	  w_nd = cm->ndidx[w];
	  nd_r_imin[w_nd] = r_imin[w];
	  nd_r_imax[w_nd] = r_imax[w];
	  nd_r_jmin[w_nd] = r_jmin[w];
	  nd_r_jmax[w_nd] = r_jmax[w];

	  r_imin[y] = imin[y];
	  r_imax[y] = imax[y];
	  r_jmin[y] = ESL_MAX(jmin[y], jmin[v]);
	  r_jmax[y] = ESL_MIN(jmax[y], jmax[v]);
	  y_nd = cm->ndidx[y];
	  nd_r_imin[y_nd] = r_imin[y];
	  nd_r_imax[y_nd] = r_imax[y];
	  nd_r_jmin[y_nd] = r_jmin[y];
	  nd_r_jmax[y_nd] = r_jmax[y];

	  if(r_jmax[w] < (r_imin[y]-1)) { 
	    cm_Fail("BEGL_S state w: %d nd: %d and BEGR_S state y: %d nd: %d  bands fail to touch, residues %d to %d cannot be emitted!\n", w, w_nd, y, y_nd, r_jmax[w]+1, r_imin[y]-1);
	    ////printf("BEGL_S state w: %d nd: %d and BEGR_S state y: %d nd: %d  bands fail to touch, residues %d to %d cannot be emitted! (unreachable cm) \n", w, w_nd, y, y_nd, r_jmax[w]+1, r_imin[y]-1);
	    return eslOK;
	  }	     
	}
	else { /* v is not reachable for any i, j */
	  v_is_r[v] = FALSE; 
	}
      }
      else { 
	sdl = StateLeftDelta(cm->sttype[v]);
	sdr = StateRightDelta(cm->sttype[v]);
	sd  = sdl + sdr;
	if((r_imin[v] <= r_imax[v] && r_jmin[v] <= r_jmax[v]) && ((r_jmax[v] - r_imin[v] + 1) >= sd)) { 
	  /* v is reachable for some i, j */
	  v_is_r[v] = TRUE; 
	  nd_is_r[nd] = TRUE;
	  child_imin = r_imin[v] + sdl;
	  child_imax = r_imax[v] + sdl;
	  child_jmin = r_jmin[v] - sdr;
	  child_jmax = r_jmax[v] - sdr;
	  if(cm->sttype[v] == IL_st) child_imax = ESL_MAX(child_imax, (imax[v]+1));
	  if(cm->sttype[v] == IR_st) child_jmin = ESL_MIN(child_jmin, (jmin[v]-1));
	  for(y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) { 
	    r_imin[y] = ESL_MIN(r_imin[y], ESL_MAX(imin[y], child_imin));
	    r_imax[y] = ESL_MAX(r_imax[y], ESL_MIN(imax[y], child_imax));
	    r_jmin[y] = ESL_MIN(r_jmin[y], ESL_MAX(jmin[y], child_jmin));
	    r_jmax[y] = ESL_MAX(r_jmax[y], ESL_MIN(jmax[y], child_jmax));
	    ////printf("\tv: %d y: %d (%11d %11d  %11d %11d)\n", v, y, r_imin[y], r_imax[y], r_jmin[y], r_jmax[y]);
	    y_nd = cm->ndidx[y];
	    nd_r_imin[y_nd] = ESL_MIN(nd_r_imin[y_nd], r_imin[y]);
	    nd_r_imax[y_nd] = ESL_MAX(nd_r_imax[y_nd], r_imax[y]);
	    nd_r_jmin[y_nd] = ESL_MIN(nd_r_jmin[y_nd], r_jmin[y]);
	    nd_r_jmax[y_nd] = ESL_MAX(nd_r_jmax[y_nd], r_jmax[y]);
	    if(cm->sttype[v] == IL_st && cm->sttype[y] == IR_st) { imax[y] = ESL_MAX(imax[y], imax[v]+1); } /* talk about a hack... ILs and IRs in the same state cause some weird situations, 
													     * here I'm enforcing that we can reach an IR from an IL in the same state no matter
													     * what, this is to avoid problems downstream */
	  }
	}
	else { /* v is not reachable for any i, j */
	  v_is_r[v] = FALSE; 
	}
      } /* end of if v != B_st */
      
      ////printf("ck v  %4s %2s %4d R %d (%11d %11d  %11d %11d)\n", Nodetype(cm->ndtype[nd]), Statetype(cm->sttype[v]), v, v_is_r[v], r_imin[v], r_imax[v], r_jmin[v], r_jmax[v]);
    } /* end of for (v) loop */
    ////printf("ck nd %4s    %4d R %d (%11d %11d  %11d %11d)\n", Nodetype(cm->ndtype[nd]), nd, nd_is_r[nd], nd_r_imin[nd], nd_r_imax[nd], nd_r_jmin[nd], nd_r_jmax[nd]);
    ////printf("\n");
    if(nd_is_r[nd] == FALSE) { 
      cm_Fail("node %d is unreachable\n", nd);
      ////printf("! node %d is unreachable cm\n", nd);
      return eslOK;
    }
  }
  ////printf("all nodes reachable\n");

  /* a final check to make sure we didn't mess up */
  /* make sure each residue can be emitted by at least one emitter */
  for(v = 0; v < cm->M; v++) { 
	sdl = StateLeftDelta(cm->sttype[v]);
	sdr = StateRightDelta(cm->sttype[v]);
	/* HERE HERE HERE */
  }


  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/* Function: HMM2ijBandsForNode()
 * 
 * Purpose:  Given HMM bands, determine the min and max bands on 
 *           i for a CM node. 
 */
int
HMM2ijBandsForNode(CM_t *cm, CP9Bands_t *cp9b, CP9Map_t *cp9map, int nd, int i0,
		   int *ret_ss_imin, int *ret_ss_imax, int *ret_ss_jmin, int *ret_ss_jmax)
{
  int ss_imin, ss_imax, ss_jmin, ss_jmax;
  int v;
  
  ss_imin = ss_jmin =  INT_MAX;
  ss_imax = ss_jmax = -INT_MAX;

  for(v = cm->nodemap[nd]; v < cm->nodemap[nd] + TotalStatesInNode(cm->ndtype[nd]); v++) { 
    if(cp9map->cs2hn[v][0] != -1) { 
      /* MATCH state */
      if(StateMapsMatch(cm->stid[v])) {
	if(StateMapsLeft(cm->stid[v])) { 
	  cp9b->imin[v] = cp9b->pn_min_m[cp9map->cs2hn[v][0]];
	  cp9b->imax[v] = cp9b->pn_max_m[cp9map->cs2hn[v][0]];
	  if(cp9b->imin[v] != -1) ss_imin = ESL_MIN(ss_imin, cp9b->imin[v]);
	  if(cp9b->imax[v] != -1) ss_imax = ESL_MAX(ss_imax, cp9b->imax[v]+1); 
	  /* +1 above allows left emission, i will increment by 1, note: we don't 
	   * have to ensure cp9b->imax[v]+1 <= (j0+1) (i think) b/c pn_max_m[k] <= j0 */
	}
	else { /* state must map right */
	  ESL_DASSERT1((StateMapsRight(cm->stid[v])));
	  cp9b->jmin[v] = cp9b->pn_min_m[cp9map->cs2hn[v][0]];
	  cp9b->jmax[v] = cp9b->pn_max_m[cp9map->cs2hn[v][0]];
	  if(cp9b->jmin[v] != -1) ss_jmin = ESL_MIN(ss_jmin, cp9b->jmin[v]-1);
	  /* -1 above allows right emission, j will decrement by 1, note: we don't 
	   * have to ensure cp9b->jmin[v]-1 >= (i0-1) (i think) b/c pn_min_m[k] >= i0 */
	  if(cp9b->jmax[v] != -1) ss_jmax = ESL_MAX(ss_jmax, cp9b->jmax[v]);
	}
	if(cm->stid[v] == MATP_MP) { /* special case, only match state type that can map to two HMM states, 
				      * also special [*1*] b/c jmin[v] >= i0+1 and jmax[v] >= i0+1 b/c 
				      * MATP_MP's emit left and right, so i emitted from MP >= i0, thus j emitted from 
				      * MP >= i0+1. 
				      */
	  ESL_DASSERT1((cp9map->cs2hn[v][1] != -1));
	  if(cp9b->pn_max_m[cp9map->cs2hn[v][1]] == i0) { /* special case (see [*1*] above), HMM tells us right half of MP state must emit first residue 
							   * in sequence (i0, this is probably very rare); in this case we ignore the HMM */
	    ESL_DASSERT1((pn_min_m[cp9map->cs2hn[v][1]] == i0));
	    cp9b->jmin[v] = -1; /* ignore hmm */
	    cp9b->jmax[v] = -1; /* ignore hmm */
	  }
	  else if (cp9b->pn_min_m[cp9map->cs2hn[v][1]] == i0) { /* special case (see [*1*] above), HMM tells us right half of MP state could possibly
						 		 * emit first residue (i0), we ignore it and say the leftmost residue it could emit is i0+1 */
	    cp9b->jmin[v] = i0 + 1; /* pad 1 onto what the hmm thought */
	    cp9b->jmax[v] = cp9b->pn_max_m[cp9map->cs2hn[v][1]]; 
	  }
	  else { 
	    cp9b->jmin[v] = cp9b->pn_min_m[cp9map->cs2hn[v][1]];
	    cp9b->jmax[v] = cp9b->pn_max_m[cp9map->cs2hn[v][1]];
	    if(cp9b->jmin[v] != -1) ss_jmin = ESL_MIN(ss_jmin, cp9b->jmin[v]-1);
	    /* -1 above allows right emission, j will decrement by 1, note: we don't 
	     * have to ensure cp9b->jmin[v]-1 >= (i0-1) (i think) b/c pn_min_m[k] >= i0 */
	    if(cp9b->jmax[v] != -1) ss_jmax = ESL_MAX(ss_jmax, cp9b->jmax[v]);
	  }
	}
      }
      /* INSERT state */
      else if (StateMapsInsert(cm->stid[v])) { 
	if(StateMapsLeft(cm->stid[v])) { 
	  cp9b->imin[v] = cp9b->pn_min_i[cp9map->cs2hn[v][0]];
	  cp9b->imax[v] = cp9b->pn_max_i[cp9map->cs2hn[v][0]];
	  if(cp9b->imin[v] != -1) ss_imin = ESL_MIN(ss_imin, cp9b->imin[v]);
	  if(cp9b->imax[v] != -1) ss_imax = ESL_MAX(ss_imax, cp9b->imax[v]+1); 
	  /* +1 above allows left emission, i will increment by 1, note: we don't 
	   * have to ensure cp9b->imax[v]+1 <= (j0+1) (i think) b/c pn_max_i[k] <= j0 */
	}
	else { /* state must map right */
	  ESL_DASSERT1((StateMapsRight(cm->stid[v])));
	  cp9b->jmin[v] = cp9b->pn_min_i[cp9map->cs2hn[v][0]];
	  cp9b->jmax[v] = cp9b->pn_max_i[cp9map->cs2hn[v][0]];
	  if(cp9b->jmin[v] != -1) ss_jmin = ESL_MIN(ss_jmin, cp9b->jmin[v]-1);
	  /* -1 above allows right emission, j will decrement by 1, note: we don't 
	   * have to ensure cp9b->jmin[v]-1 >= (i0-1) (i think) b/c pn_min_i[k] >= i0 */
	  if(cp9b->jmax[v] != -1) ss_jmax = ESL_MAX(ss_jmax, cp9b->jmax[v]);
	}
      }
      /* DELETE state */
      else if(StateMapsDelete(cm->stid[v])) { 
	if(StateMapsLeft(cm->stid[v])) { 
	  cp9b->imin[v] = (cp9b->pn_min_d[cp9map->cs2hn[v][0]] == -1) ? -1 : cp9b->pn_min_d[cp9map->cs2hn[v][0]] + 1;
	  cp9b->imax[v] = (cp9b->pn_max_d[cp9map->cs2hn[v][0]] == -1) ? -1 : cp9b->pn_max_d[cp9map->cs2hn[v][0]] + 1;
	  if(cp9b->imin[v] != -1) ss_imin = ESL_MIN(ss_imin, cp9b->imin[v]);
	  if(cp9b->imax[v] != -1) ss_imax = ESL_MAX(ss_imax, cp9b->imax[v]);
	}
	else { /* state must map right */
	  ESL_DASSERT1((StateMapsRight(cm->stid[v])));
	  cp9b->jmin[v] = cp9b->pn_min_d[cp9map->cs2hn[v][0]];
	  cp9b->jmax[v] = cp9b->pn_max_d[cp9map->cs2hn[v][0]];
	  if(cp9b->jmin[v] != -1) ss_jmin = ESL_MIN(ss_jmin, cp9b->jmin[v]);
	  if(cp9b->jmax[v] != -1) ss_jmax = ESL_MAX(ss_jmax, cp9b->jmax[v]);
	}

	if(cm->ndtype[cm->ndidx[v]] == MATP_nd) { /* special case 3 delete state types (MATP_ML, MATP_MR, MATP_D) that can map to two HMM states */
	  ESL_DASSERT1((cp9map->cs2hn[v][1] != -1));
	  if(cm->stid[v] == MATP_ML || cm->stid[v] == MATP_MR) { 
	    /* really special case [*2*] because v emits (left in case of MATP_ML, right in case of MATP_MR), 
	     * so jmin[v] >= i0 and jmax[v] >= i0 (can't be i0-1). j == i0-1 is non-sensical because we have 
	     * to emit at least 1  residue, so d >= 1 thus i <= j, so if j == i0-1, i would be i0-1 which doesn't make sense. 
	     * This is similar to the case for MATP_MP above (see [*1*]) 
	     */
 	    if(cp9b->pn_max_d[cp9map->cs2hn[v][1]] == (i0-1)) { /* special case (see [*2*] above), HMM tells us this ML or MR state can't emit! 
								 * (this is probably very rare); in this case we ignore the HMM */
	      ESL_DASSERT1((pn_min_d[cp9map->cs2hn[v][1]] == i0-1));
	      cp9b->jmin[v] = -1; /* ignore hmm */
	      cp9b->jmax[v] = -1; /* ignore hmm */
	    }
	    else if (cp9b->pn_min_d[cp9map->cs2hn[v][1]] == i0-1) { /* special case (see [*2*] above), HMM tells us this ML or MR state could possibly
								     * emit 0 residues, we ignore it and say the leftmost residue it could emit is i0 */
	      cp9b->jmin[v] = i0; /* pad 1 onto what the hmm thought */
	      cp9b->jmax[v] = cp9b->pn_max_d[cp9map->cs2hn[v][1]];
	    }
	    else { 
	      cp9b->jmin[v] = cp9b->pn_min_d[cp9map->cs2hn[v][1]];
	      cp9b->jmax[v] = cp9b->pn_max_d[cp9map->cs2hn[v][1]];
	      if(cp9b->jmin[v] != -1) ss_jmin = ESL_MIN(ss_jmin, cp9b->jmin[v]);
	      if(cp9b->jmax[v] != -1) ss_jmax = ESL_MAX(ss_jmax, cp9b->jmax[v]);
	    }
	  } /* end of if cm->stid[v] == MATP_ML || cm->stid[v] == MATP_MR */
	  else { /* node type is MATP_nd, but we're not a MATP_ML or MATP_MR, we must be a MATP_D, don't have same issues as MATP_ML or MATP_MR b/c we don't emit */
	    cp9b->jmin[v] = cp9b->pn_min_d[cp9map->cs2hn[v][1]];
	    cp9b->jmax[v] = cp9b->pn_max_d[cp9map->cs2hn[v][1]];
	    if(cp9b->jmin[v] != -1) ss_jmin = ESL_MIN(ss_jmin, cp9b->jmin[v]);
	    if(cp9b->jmax[v] != -1) ss_jmax = ESL_MAX(ss_jmax, cp9b->jmax[v]);
	  }
	}
      }
    }
    printf("OLD v:  %4d h2ij (%11d %11d  %11d %11d)\n", v, cp9b->imin[v], cp9b->imax[v], cp9b->jmin[v], cp9b->jmax[v]);
  }

  printf("OLD nd: %4d h2ij (%11d %11d  %11d %11d)\n", nd, ss_imin, ss_imax, ss_jmin, ss_jmax);

  if(ret_ss_imin != NULL) *ret_ss_imin = (ss_imin !=  INT_MAX) ? ss_imin : -1; 
  if(ret_ss_imax != NULL) *ret_ss_imax = (ss_imax != -INT_MAX) ? ss_imax : -1; 
  if(ret_ss_jmin != NULL) *ret_ss_jmin = (ss_jmin !=  INT_MAX) ? ss_jmin : -1; 
  if(ret_ss_jmax != NULL) *ret_ss_jmax = (ss_jmax != -INT_MAX) ? ss_jmax : -1; 
  return eslOK;
}

  

/* Function: HMM2ijBandsForNode_new()
 * Incept:   EPN, Fri Feb  1 14:31:10 2008
 * 
 * Purpose:  Set minimum and maximum i and j bands for all states in 
 *           node <nd>. This is done using the HMM bands in <cp9b>.
 *           We consider CM states as having two 'sides' a left and
 *           a right. Some CM states map to two HMM states, one on the
 *           left side and one on the right (MATP_MP, MATP_ML, MATP_MP
 *           and MATP_D). Some CM states map to one HMM state on the
 *           left and zero on the right (MATL_ML, MATL_D, and all IL 
 *           states). Some states map to one HMM state on the right
 *           and zero on the left (MATR_MR, MATR_D, and all IR states).
 *           And some CM states map to zero states on both the left and
 *           the right (S, B, and E states), for these states no 
 *           action is taken. 
 *
 *           First, i bands (imin/imax) are set for CM states v which
 *           map to a HMM state on the left. Then j bands (jmin/jmax) are 
 *           set for v which map to an HMM state on the right.
 */
int
HMM2ijBandsForNode_new(CM_t *cm, CP9Bands_t *cp9b, CP9Map_t *cp9map, int nd, int i0, int *ret_nd_imin, int *ret_nd_imax, int *ret_nd_jmin, int *ret_nd_jmax)
{
  int nd_imin, nd_imax, nd_jmin, nd_jmax; /* min/max i/j values in this node <nd>, 
					     set based on min/max i/j values for states v within <nd> */
  int v;   /* state index, ranges across all states in node <nd> */
  int hn1; /* first  HMM node that v maps to, will be -1 (flag for invalid) if v 
	    * is not an emitter or delete (-1 for S, B and E states) */
  int hn2; /* second HMM node that v maps to, will be -1 (flag for invalid) if v 
	      is a non-insert state in a MATP node, only MATP_MP, MATP_ML, MATP_MR, MATP_D map 
	      to 2 different HMM nodes, and these sometimes require special care (see special cases below) */
  int sd;  /* number of residues emitted from an HMM state that v maps to, either 1 or 0 */

  nd_imin = nd_jmin =  INT_MAX; 
  nd_imax = nd_jmax = -INT_MAX;

  for(v = cm->nodemap[nd]; v < cm->nodemap[nd] + TotalStatesInNode(cm->ndtype[nd]); v++) { 
    if(!StateIsDetached(cm, v)) { /* detached inserts are unreachable, skip them */
      hn1 = cp9map->cs2hn[v][0]; /* -1 only if v is a S, B, or E state */
      hn2 = cp9map->cs2hn[v][1]; /* -1 unless v is in non-insert state in a MATP node (either MATP_MP, MATP_ML, MATP_MR, MATP_MP */
      /* set i bands for state v, if left side of v maps to an HMM state */
      switch(cm->stid[v]) { 
      case MATP_MP:
      case MATP_ML:
      case MATL_ML:
	cp9b->imin[v] = cp9b->pn_min_m[hn1];
	cp9b->imax[v] = cp9b->pn_max_m[hn1];
	sd = 1;
	break;
	
      case MATP_IL:
      case MATL_IL:
      case ROOT_IL: 
      case BEGR_IL: 
	cp9b->imin[v] = cp9b->pn_min_i[hn1];
	cp9b->imax[v] = cp9b->pn_max_i[hn1];
	sd = 1;
	break;
	
      case MATP_D:
      case MATP_MR:
      case MATL_D:
	cp9b->imin[v] = (cp9b->pn_min_d[hn1] == -1) ? -1 : cp9b->pn_min_d[hn1] + 1;
	cp9b->imax[v] = (cp9b->pn_max_d[hn1] == -1) ? -1 : cp9b->pn_max_d[hn1] + 1;
	/* the plus 1 is to deal with an off-by-one between the HMM and CM, an HMM parse that has a delete state D_k
	 * entered at residue i means that residue i has already been emitted by a match or insert state to the left
	 * of D_k. In a CM, if a delete state v has a subtree from i..j, this means that residue i has yet to be 
	 * emitted, so we need to correct for this with the + 1 above (but only for imin, imax, not for jmin, jmax below). 
	 */
	sd = 0;
	break;
	
	if(cp9b->imin[v] != -1) nd_imin = ESL_MIN(nd_imin, cp9b->imin[v]);
	if(cp9b->imax[v] != -1) nd_imax = ESL_MAX(nd_imax, cp9b->imax[v]+sd); 
      }

      /* set j bands for state v, if right side of v maps to an HMM state */
      switch(cm->stid[v]) { 
      case MATR_MR:
	cp9b->jmin[v] = cp9b->pn_min_m[hn1];
	cp9b->jmax[v] = cp9b->pn_max_m[hn1];
	sd = 1;
	break;
	
      case MATP_IR:
      case MATR_IR: 
      case ROOT_IR: 
	cp9b->jmin[v] = cp9b->pn_min_i[hn1];
	cp9b->jmax[v] = cp9b->pn_max_i[hn1];
	sd = 1;
	break;
	
      case MATR_D:
	cp9b->jmin[v] = cp9b->pn_min_d[hn1];
	cp9b->jmax[v] = cp9b->pn_max_d[hn1];
	sd = 0;
	break;
	
      case MATP_MR: /* maps to 2 HMM states, right side (this side) maps to cp9map->cs2hn[v][*1*] */
	cp9b->jmin[v] = cp9b->pn_min_m[hn2];
	cp9b->jmax[v] = cp9b->pn_max_m[hn2];
	sd = 1;
	break;

      case MATP_D: /* maps to 2 HMM states, right side (this side) maps to cp9map->cs2hn[v][*1*] */
	cp9b->jmin[v] = cp9b->pn_min_d[hn2];
	cp9b->jmax[v] = cp9b->pn_max_d[hn2];
	sd = 0;
      break;

      /* the two special cases, MATP_MP and MATP_ML, which also map to HMM states that emit on the
       * left, so we have to be careful about boundary conditions (see comments below).
       */
      case MATP_MP:
	/* special case [*1*]: v emits left and right, so jmin[v] >= i0+1 and jmax[v] >= i0+1
	 * b/c i emitted from MATP_MP >= i0, thus j emitted from MATP_MP >= i0+1. 
	 */
	if(cp9b->pn_max_m[hn2] == i0) { /* HMM tells us right half of MP state must emit first residue in the sequence,
					 * but we know it can't because the left half of this MATP_MP state must emit 1 residue, 
					 * which can't be before the first one. In this case we ignore the HMM and set
					 * the j band to dummy values which means it's unset, i.e. it doesn't yet exist. */
	  ESL_DASSERT1((pn_min_m[hn2] == i0));
	  cp9b->jmin[v] = -1; /* ignore hmm */
	  cp9b->jmax[v] = -1; /* ignore hmm */
	}
	else if (cp9b->pn_min_m[hn2] == i0) { /* HMM tells us right half of MP state could possibly
					       * emit first residue (i0), but we know it can't (see comment above).
					       * we ignore it and say the leftmost residue it could emit is i0+1 */
	  cp9b->jmin[v] = i0 + 1; /* pad 1 onto what the hmm thought */
	  cp9b->jmax[v] = cp9b->pn_max_m[hn2]; /* leave this, we konw it's not == i0, we checked for that case above */
	}
	else { 
	  cp9b->jmin[v] = cp9b->pn_min_m[hn2]; 
	  cp9b->jmax[v] = cp9b->pn_max_m[hn2];
	}
	sd = 1;
	break;
	
      case MATP_ML:
	/* special case [*2*]: v emits left, so jmin[v] >= i0 and jmax[v] >= i0 
	 * b/c i emitted from MATP_ML >= i0, thus j must be >= i0
	 * This is similar to the case for MATP_MP above (see [*1*]) 
	 */
	if(cp9b->pn_max_d[hn2] == (i0-1)) { /* HMM tells us we can only enter this state having emitted 0 residues,
					     * we know better, in this case we ignore the HMM, and set the j band
					     * to dummy values, which means it's unset and doesn't yet exist
					     */
	  ESL_DASSERT1((pn_min_d[hn2] == i0-1));
	  cp9b->jmin[v] = -1; /* ignore hmm */
	  cp9b->jmax[v] = -1; /* ignore hmm */
	}
	else if (cp9b->pn_min_d[hn2] == (i0-1)) { /* HMM tells us right half of MATP_ML state could be entered
						   * having emitted 0 residues, but we know it can't (see comment above).
						   * we ignore it and say we must have at least emitted residue i0.
						   */
	  cp9b->jmin[v] = i0; /* pad 1 onto what the hmm thought */
	  cp9b->jmax[v] = cp9b->pn_max_d[hn2]; /* leave this, we know it's not == i0-1, we checked for that case above */
	}
	else { 
	  cp9b->jmin[v] = cp9b->pn_min_d[hn2];
	  cp9b->jmax[v] = cp9b->pn_max_d[hn2];
	}
        sd = 0; 
	break;
      }
      
      if(cp9b->jmin[v] != -1) nd_jmin = ESL_MIN(nd_jmin, cp9b->jmin[v]-sd);
      if(cp9b->jmax[v] != -1) nd_jmax = ESL_MAX(nd_jmax, cp9b->jmax[v]);
      /* '-sd' in the jmin setting above allows right emission for emitters, j will decrement by 1, 
       * note: we don't have to ensure cp9b->jmin[v]-sd >= (i0-1) (i think) b/c pn_min_{m,i}[k] >= i0 for all k >= 1 */
      ////printf("NEW v:  %4d h2ij (%11d %11d  %11d %11d)\n", v, cp9b->imin[v], cp9b->imax[v], cp9b->jmin[v], cp9b->jmax[v]);
    }
  }
  
  ////printf("NEW nd: %4d h2ij (%11d %11d  %11d %11d)\n", nd, nd_imin, nd_imax, nd_jmin, nd_jmax);
  ////HMM2ijBandsForNode(cm, cp9b, cp9map, nd, i0, NULL, NULL, NULL, NULL);
  ////printf("\n");

  if(ret_nd_imin != NULL) *ret_nd_imin = (nd_imin !=  INT_MAX) ? nd_imin : -1; 
  if(ret_nd_imax != NULL) *ret_nd_imax = (nd_imax != -INT_MAX) ? nd_imax : -1; 
  if(ret_nd_jmin != NULL) *ret_nd_jmin = (nd_jmin !=  INT_MAX) ? nd_jmin : -1; 
  if(ret_nd_jmax != NULL) *ret_nd_jmax = (nd_jmax != -INT_MAX) ? nd_jmax : -1; 
  return eslOK;
}

/* Function: HMM2ijBandsForCpos()
 * Incept:   EPN, Mon Feb  4 05:31:43 2008
 * 
 */
int
HMM2ijBandsForCpos(CM_t *cm, CP9Bands_t *cp9b, CP9Map_t *cp9map, int k, int i0, int j0, int r_mn, int r_mx, int r_in, int r_ix, int r_dn, int r_dx, int *r_nxt_nnA, int *r_nxt_nxA)
{
  int imin, imax, jmin, jmax;
  int v1, v2;   
  int r_nxt_nn, r_nxt_nx;
  int sd;  /* number of residues emitted from an HMM state that v maps to, either 1 or 0 */
  int v, nd;

  r_nxt_nn =  INT_MAX;
  r_nxt_nx = -INT_MAX;

  /* match state */
  if(r_mn <= r_mx) { /* M_k is reachable for i = r_mn[k]..r_mx[k] */
    imin = jmin = r_mn;
    imax = jmax = r_mx;
  }
  else { 
#if 0
    imin = i0+1; 
    imax = i0; 
    jmin = j0+1;
    jmax = j0;
#endif 
    imin = imax = jmin = jmax = -2;
  }
  
  v1 = cp9map->hns2cs[k][HMMMATCH][0];
  v2 = cp9map->hns2cs[k][HMMMATCH][1];
  printf("k: %4d MATCH  v1: %4d %4s %2s v2: %4d %4s %2s\n", k, v1, ((v1 == -1) ? "" : Nodetype(cm->ndtype[cm->ndidx[v1]])), ((v1 == -1) ? "" : Statetype(cm->sttype[v1])), v2, ((v2 == -1) ? "" : Nodetype(cm->ndtype[cm->ndidx[v2]])), ((v2 == -1) ? "" : Statetype(cm->sttype[v2])));
  
  switch(cm->sttype[v1]) { 
  case S_st: /* special case the ROOT_S CM state maps to the M_0 HMM state */
    break;
  case MP_st: 
    if(cp9map->nd2lpos[cm->ndidx[v1]] == k) { 
      cp9b->imin[v1] = imin;
      cp9b->imax[v1] = imax;
      r_nxt_nn = ESL_MIN(r_nxt_nn, imin + 1);
      r_nxt_nx = ESL_MAX(r_nxt_nx, imax + 1);
    }
    else { 
      assert(cp9map->nd2rpos[cm->ndidx[v1]] == k);
      cp9b->jmin[v1] = jmin;
      cp9b->jmax[v1] = jmax;
      r_nxt_nn = ESL_MIN(r_nxt_nn, jmin - 1);
      r_nxt_nx = ESL_MAX(r_nxt_nx, jmax - 1);
    }
    break;
  case ML_st: 
    cp9b->imin[v1] = imin;
    cp9b->imax[v1] = imax;
    r_nxt_nn = ESL_MIN(r_nxt_nn, imin + 1);
    r_nxt_nx = ESL_MAX(r_nxt_nx, imax + 1);
    break;
  case MR_st: 
    cp9b->jmin[v1] = jmin;
    cp9b->jmax[v1] = jmax;
    r_nxt_nn = ESL_MIN(r_nxt_nn, jmin - 1);
    r_nxt_nx = ESL_MAX(r_nxt_nx, jmax - 1);
    assert(cm->ndtype[cm->ndidx[v1]] == MATR_nd);
    cp9b->imin[v] = r_nxt_nnA[cm->ndidx[v1]];
    cp9b->imax[v] = r_nxt_nxA[cm->ndidx[v1]];
    break;
  default:
    cm_Fail("0 bogus state type (%d) for state %d\n", cm->sttype[v1], v1);
  }
  
  if(v2 != -1) { 
    switch(cm->sttype[v2]) { 
    case ML_st: 
      cp9b->imin[v2] = imin;
      cp9b->imax[v2] = imax;
      r_nxt_nn = ESL_MIN(r_nxt_nn, imin + 1);
      r_nxt_nx = ESL_MAX(r_nxt_nx, imax + 1);
      break;
    case MR_st: 
      cp9b->jmin[v2] = jmin;
      cp9b->jmax[v2] = jmax;
      r_nxt_nn = ESL_MIN(r_nxt_nn, jmin - 1);
      r_nxt_nx = ESL_MAX(r_nxt_nx, jmax - 1);
      break;
    default:
      cm_Fail("1 bogus state type (%d) for state %d\n", cm->sttype[v1], v1);
    }
  }

  /* insert state */
  v1 = cp9map->hns2cs[k][HMMINSERT][0];
  assert(cp9map->hns2cs[k][HMMINSERT][1] == -1);
  printf("k: %4d INSERT v1: %4d %4s %2s v2: %4d %4s %2s\n", k, v1, ((v1 == -1) ? "" : Nodetype(cm->ndtype[cm->ndidx[v1]])), ((v1 == -1) ? "" : Statetype(cm->sttype[v1])), v2, ((v2 == -1) ? "" : Nodetype(cm->ndtype[cm->ndidx[v2]])), ((v2 == -1) ? "" : Statetype(cm->sttype[v2])));

  if(r_in <= r_ix) { /* I_k is reachable for i = r_in..r_ix */
    imin = jmin = r_in;
    imax = jmax = r_ix;
  }
  else { 
#if 0
    imin = i0+1; 
    imax = i0; 
    jmin = j0+1;
    jmax = j0;
#endif 
    imin = imax = jmin = jmax = -2;
  }

  switch(cm->sttype[v1]) { 
  case IL_st: 
    cp9b->imin[v1] = imin;
    cp9b->imax[v1] = imax;
    r_nxt_nn = ESL_MIN(r_nxt_nn, imin + 1);
    r_nxt_nx = ESL_MAX(r_nxt_nx, imax + 1);
    break;
  case IR_st: 
    cp9b->jmin[v1] = jmin;
    cp9b->jmax[v1] = jmax;
    r_nxt_nn = ESL_MIN(r_nxt_nn, jmin - 1);
    r_nxt_nx = ESL_MAX(r_nxt_nx, jmax - 1);
    break;
  default:
    cm_Fail("2 bogus state type (%d) for state %d\n", cm->sttype[v1], v1);
  }

  if(k > 0) { /* no D_0 state */
    /* delete state */
    v1 = cp9map->hns2cs[k][HMMDELETE][0];
    v2 = cp9map->hns2cs[k][HMMDELETE][1];
    printf("k: %4d DELETE v1: %4d %4s %2s v2: %4d %4s %2s\n", k, v1, ((v1 == -1) ? "" : Nodetype(cm->ndtype[cm->ndidx[v1]])), ((v1 == -1) ? "" : Statetype(cm->sttype[v1])), v2, ((v2 == -1) ? "" : Nodetype(cm->ndtype[cm->ndidx[v2]])), ((v2 == -1) ? "" : Statetype(cm->sttype[v2])));

    if(r_dn <= r_dx) { /* M_k is reachable for i = r_dn[k]..r_dx[k] */
      imin = jmin = r_dn;
      imax = jmax = r_dx;
    }
    else { 
#if 0
      imin = i0+1; 
      imax = i0; 
      jmin = j0+1;
      jmax = j0;
#endif 
      imin = imax = jmin = jmax = -2;
    }
  
    switch(cm->sttype[v1]) { 
    case D_st: 
      if(cp9map->nd2lpos[cm->ndidx[v1]] == k) { 
	cp9b->imin[v1] = imin;
	cp9b->imax[v1] = imax;
	r_nxt_nn = ESL_MIN(r_nxt_nn, imin);
	r_nxt_nx = ESL_MAX(r_nxt_nx, imax);
      }
      else { 
	assert(cp9map->nd2rpos[cm->ndidx[v1]] == k);
	cp9b->jmin[v1] = jmin;
	cp9b->jmax[v1] = jmax;
	r_nxt_nn = ESL_MIN(r_nxt_nn, jmin);
	r_nxt_nx = ESL_MAX(r_nxt_nx, jmax);
      }
      break;
    case ML_st: 
      cp9b->jmin[v1] = imin;
      cp9b->jmax[v1] = imax;
      r_nxt_nn = ESL_MIN(r_nxt_nn, jmin);
      r_nxt_nx = ESL_MAX(r_nxt_nx, jmax);
      break;
    case MR_st: 
      cp9b->imin[v1] = imin;
      cp9b->imax[v1] = imax;
      r_nxt_nn = ESL_MIN(r_nxt_nn, imin);
      r_nxt_nx = ESL_MAX(r_nxt_nx, imax);
      break;
    default:
      cm_Fail("3 bogus state type (%d) for state %d\n", cm->sttype[v1], v1);
    }
  
    if(v2 != -1) { 
      switch(cm->sttype[v2]) { 
      case D_st: 
	if(cp9map->nd2lpos[cm->ndidx[v2]] == k) { 
	  cp9b->imin[v2] = imin;
	  cp9b->imax[v2] = imax;
	  r_nxt_nn = ESL_MIN(r_nxt_nn, imin);
	  r_nxt_nx = ESL_MAX(r_nxt_nx, imax);
	}
	else { 
	  assert(cp9map->nd2rpos[cm->ndidx[v2]] == k);
	  cp9b->jmin[v2] = jmin;
	  cp9b->jmax[v2] = jmax;
	  r_nxt_nn = ESL_MIN(r_nxt_nn, jmin);
	  r_nxt_nx = ESL_MAX(r_nxt_nx, jmax);
	}
      break;
      case ML_st: 
	cp9b->imin[v2] = imin;
	cp9b->imax[v2] = imax;
	r_nxt_nn = ESL_MIN(r_nxt_nn, imin);
	r_nxt_nx = ESL_MAX(r_nxt_nx, imax);
	break;
      case MR_st: 
	cp9b->jmin[v2] = jmin;
	cp9b->jmax[v2] = jmax;
	r_nxt_nn = ESL_MIN(r_nxt_nn, jmin);
	r_nxt_nx = ESL_MAX(r_nxt_nx, jmax);
	break;
      default:
	cm_Fail("4 bogus state type (%d) for state %d\n", cm->sttype[v2], v2);
      }
    }
  }
  /* now update the r_nxt_nnA and r_nxt_nxA arrays, which keep track of which residues are reachable in the next HMM node,
   * and fill in j bands for any CM MATL nodes in the CM between this */
  int nxt_nd;
  if(k < cp9b->hmm_M) nxt_nd = cp9map->pos2nd[k+1];
  else                nxt_nd = cm->nodes;
  printf("! k: %4d (nd: %4d, nxt: %4d)\n", k, cp9map->pos2nd[k], nxt_nd);
  ///if(cm->ndtype[nd] == MATP_nd) { /* special case, if we're on the right side of a MATP we can now set the MATP_IL's jmin/jmax and MATP_IR's imin/imax */ 
  nd = cp9map->pos2nd[k] + 1; 
  while(nd != nxt_nd) { 
    r_nxt_nnA[nd] = r_nxt_nn;
    r_nxt_nxA[nd] = r_nxt_nx;
    if(cm->ndtype[nd] == MATL_nd) { /* set jmin, jmax */
      printf("\tsetting jmin/jmax for MATL node: %d\n", nd);
      for(v = cm->nodemap[nd]; v < cm->nodemap[nd] + TotalStatesInNode(cm->ndtype[nd]); v++) { 
	if(cp9b->imax[v] - cp9b->imin[v] >= 0) { 
	  cp9b->jmin[v] = r_nxt_nn;
	  cp9b->jmax[v] = r_nxt_nx;
	}
      }
    }
    if(nxt_nd < nd) nd--;
    else nd++;
  }
  return eslOK;
}


#endif
/* EPN, Thu Feb  7 17:56:18 2008
   Got rid of: 
   extern void  CPlan9GlobalConfig(CP9_t *hmm);
   
   The reason I wanted to get rid of this CPlan9GlobalConfig() call is b/c I've changed
    how CP9's are locally configured, and the M_0->I_0, M_0->D_1, M_M->I_M and D_M->I_M transitions
    are all set to IMPOSSIBLE (to make a local CP9 more like a local CM), and thus it makes it
    hard to change a locally configured CP9 back to global mode b/c the initial values of those
    transitions is lost (the solution would be to add vectors to the cp9 data structure that
    remember these transition probs, but I don't have to do that if I NEVER need to globalize 
    a locally configured CP9.

*/
#if 0
/* Function: CPlan9GlobalConfig()
 * EPN 09.24.06
 * based on SRE's Plan7GlobalConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to (Needleman/Wunsch) configuration.
 *           Make all transitions to EL states impossible.
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9GlobalConfig(CP9_t *hmm)
{
  cm_Fail("you can't use CPlan9GlobalConfig(), unless you've stored what transitions M_0->I_0, M_0->D_1, M_M->I_M and D_M->I_M were before you localized this HMM!\n");
  int k;
  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* Configure entry.
   * Exactly 3 ways to start, B->M_1 (hmm->begin[1]), B->I_0 (hmm->t[0][CTMI]),
   *                      and B->D_1 (hmm->t[0][CTMD])
   */
  hmm->begin[1] = 1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTMEL]); 
  /* this is okay, hmm->t[0] is never changed, even during local
   * configuration */
  esl_vec_FSet(hmm->begin+2, hmm->M-1, 0.);
  
  hmm->end[hmm->M] = 1. - hmm->t[hmm->M][CTMI];
  esl_vec_FSet(hmm->end+1, hmm->M-1, 0.);
  CPlan9RenormalizeExits(hmm, 1);

  /* Make all transitions to EL impossible, node 0, M special and should 
   * always have CTMEL transition as impossible. */
  for(k = 1; k < hmm->M; k++)
    {
      hmm->t[k][CTMEL] = 0.;
      esl_vec_FNorm(hmm->t[k], 4); /* renormalize transitions out of node k */
    }
  hmm->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
  hmm->flags       &= ~CPLAN9_LOCAL_BEGIN; /* local begins now off */
  hmm->flags       &= ~CPLAN9_LOCAL_END;   /* local ends now off */
  hmm->flags       &= ~CPLAN9_EL;          /* EL end locals now off */

  CP9Logoddsify(hmm);
}
#endif
/* EPN, Fri Feb 15 12:37:01 2008
 * Changed fitting gumbels in cmcalibrate to fitting exponential tails. 
 * Some Gumbel functions are now unused.
 */
#if 0


/*
 * Function: RJK_ExtremeValueE
 * Date:     RJK, Mon Sep 30, 2002 [St. Louis]
 * Purpose:  Given a score (x), mu, and lambda, calculates 
 *           E=exp(-1*lambda(x-mu)) using first part of code from Sean's
 *           ExtremeValueP
 */
double RJK_ExtremeValueE (float x, double mu, double lambda) {
                        /* avoid underflow fp exceptions near P=0.0*/
  if ((lambda * (x - mu)) >= 2.3 * (double) DBL_MAX_10_EXP) 
    return 0.0;
  else 
    return(exp(-1. * lambda * (x - mu)));
}

/* fit_histogram()
 * Create, fill and fit a histogram to a gumbel. Data to fill the histogram
 * is given as <data>.
 */
static int
fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores,
	      double *ret_mu, double *ret_lambda)
{
  int status;
  double mu;
  double lambda;
  int i;
  double *xv;         /* raw data from histogram */
  int     n,z;  
  float tailfit;
  double mufix;

  ESL_HISTOGRAM *h = NULL;       /* histogram of scores */


  /* Initialize histogram; these numbers are guesses */
  if((h = esl_histogram_CreateFull(-100., 100., .1)) == NULL) return eslEMEM;    

  /* fill histogram */
  for(i = 0; i < nscores; i++)
    if((status = esl_histogram_Add(h, scores[i])) != eslOK) ESL_FAIL(status, errbuf, "fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);

  /* fit scores to a exp tail */
  tailfit = esl_opt_GetReal(go, "--exp-tail");
  esl_histogram_GetTailByMass(h, tailfit, &xv, &n, &z); /* fit to right 'tailfit' fraction, 0.5 by default */
  esl_exptail_FitCensored(xv, n, z, xv[0], &mu, &lambda);
  esl_exptail_FitCensoredLoc(xv, n, z, xv[0], 0.693147, &mufix);

  /* print to output files if nec */
  if(cfg->gumhfp != NULL)
    esl_histogram_Plot(cfg->gumhfp, h);
  if(cfg->gumqfp != NULL) {
      double  params[2];  
      params[0] = mu;
      params[1] = lambda;
      esl_histogram_PlotQQ(cfg->gumqfp, h, &esl_exp_generic_invcdf, params);
  }

  if (cfg->gumsfp != NULL) {
    esl_histogram_PlotSurvival(cfg->gumsfp, h);
    esl_exp_Plot(cfg->gumsfp, mu,    lambda,   esl_exp_surv, h->xmin - 5., h->xmax + 5., 0.1);
    esl_exp_Plot(cfg->gumsfp, mufix, 0.693147, esl_exp_surv, h->xmin - 5., h->xmax + 5., 0.1);
  }

  esl_histogram_Destroy(h);

  *ret_mu     = mu;
  *ret_lambda = lambda;
  return eslOK;
}

#endif

#if 0

/* Structure CP9FThresh_t: CP9 HMM filter thresholds, determined empirically
 * by sampling from the CM
 */
typedef struct cp9filterthr_s {
  int   N;             /* number of CM hits used to get threshold ((N*F) passed)*/
  float cm_eval;       /* CM E-value threshold, we rejected worse than   */
  float l_eval;        /*  local CP9 scanning E-value threshold    */
  float g_eval;        /* glocal CP9 scanning E-value threshold    */
  float l_F;           /* fraction of empirical CM hits survive filter for l_eval cutoff */
  float g_F;           /* fraction of empirical CM hits survive filter for g_eval cutoff */
  int   db_size;       /* db size used to calculate exp tail mu for *_eval calculations */
  int   was_fast;      /* TRUE if hacky fast method for calcing thresholds was used */
  int   isvalid;       /* TRUE if values have been set, FALSE if not */
} CP9FilterThr_t;


/* Structure SubFilterInfo_t: Information on possible sub CM filters for a CM.                           
 * States of a CM are grouped into 'start groups'. There is one start group                              
 * for each start state of the CM. A 'start group' begins with a start state and ends                    
 * with a E or B state, and includes all states in between.                                              
 */                                                                                                      
typedef struct subfilterinfo_s {
  int    M;            /* # states in the CM */                                                          
  int    nstarts;      /* # start states (and start groups) in the CM */                                 
  int    ncands;       /* number of candidate states, these *could* be sub CM roots */                   
  double beta;         /* beta used for calculating avglenA */                                           
  float  minlen;       /* minimum average length (avglen) a candidate state must have */                 
  int   *iscandA;      /* [0..v..cm->M-1] TRUE if state v is a candidate sub CM root, FALSE otherwise */   
  float *avglenA;      /* [0..v..cm->M-1] average length of a hit rooted at v (from QDB) */                
  int   *startA;       /* [0..i..cm->M-1] start group this state belongs to */                               
  int   *firstA;       /* [0..i..nstarts-1], first state in start state i's group */                     
  int   *lastA;        /* [0..i..nstarts-1], last state in start state i's group */                      
  int  **withinAA;     /* [0..i..nstarts-1][0..j..nstarts-1] = TRUE if start state j's group             
                        * is within start state i's group.                                               
                        *  emap->startA[cm->nodemap[i]]->lpos < emap->startA[cm->nodemap[j]]->lpos  &&   
                        *  emap->endA  [cm->nodemap[i]]->rpos > emap->endA  [cm->nodemap[j]]->rpos       
                        */			
} SubFilterInfo_t;

#endif

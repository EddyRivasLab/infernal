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

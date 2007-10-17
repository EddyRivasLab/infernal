/* cm_fastalign.c
 * EPN, Wed Oct 10 07:20:48 2007
 * 
 * Fast versions of CYK and Inside search functions.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

#define TSC(s,k) (tsc[(v) * MAXCONNECT + (s)])
#define AMX(j,v,d) (alphap[(j * cm->M * (W+1)) + ((v) * (W+1) + d)])



/*
 * Function: Xcp9_seq2bands
 * Date    : EPN, Mon Jan  8 07:23:34 2007
 *           EPN, Wed Oct 17 04:53:58 2007  [updated/optimized]
 *
 * Purpose:  Given a CM with precalc'ed CP9 HMM and CP9Map, a sequence and 
 *           a CP9Bands_t structure, calculate the HMM bands and store them
 *           in the CP9Bands_t structure.
 *           
 * Args:     cm          - the covariance model
 *           dsq         - sequence in digitized form
 *           i0          - start of target subsequence (often 1, beginning of sq)
 *           j0          - end of target subsequence (often L, end of sq)
 *           cp9b        - PRE-ALLOCATED, the HMM bands for this sequence, filled here.
 *           ret_cp9_post - RETURN: the HMM posterior matrix (NULL if not wanted)
 *           debug_level - verbosity level for debugging printf()s
 * Return:  void
 */
void 
Xcp9_seq2bands(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, 
	       CP9_dpmatrix_t **ret_cp9_post, int debug_level)
{
  int             status;
  ESL_STOPWATCH  *watch;    /* for timings if cm->align_opts & CM_ALIGN_TIME             */
  int             use_sums; /* TRUE to fill and use posterior sums during HMM band calc  *
			     * leads to wider bands                                      */
  int             v;        /* state index                                               */
  int             sc;

  /* Contract checks */
  if(cm->cp9 == NULL)
    cm_Fail("in Xcp9_seq2bands, but cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL)
    cm_Fail("in Xcp9_seq2bands, but cm->cp9map is NULL.\n");
  if((cm->align_opts & CM_ALIGN_HBANDED) && (cm->search_opts & CM_SEARCH_HBANDED)) 
    cm_Fail("in Xcp9_seq2bands, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both up, exactly 1 must be up.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED))) 
    cm_Fail("in Xcp9_seq2bands, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down, exactly 1 must be up.\n");
  if((cm->search_opts & CM_SEARCH_HMMSCANBANDS) && (!(cm->search_opts & CM_SEARCH_HBANDED))) 
    cm_Fail("in Xcp9_seq2bands, CM_SEARCH_HMMSCANBANDS flag raised, but not CM_SEARCH_HBANDED flag, this doesn't make sense\n");
  if(dsq == NULL) 
    cm_Fail("in Xcp9_seq2bands, dsq is NULL.");
  
  use_sums = FALSE;
  if((cm->align_opts & CM_ALIGN_SUMS) || (cm->search_opts & CM_SEARCH_SUMS))
    use_sums = TRUE;
    
  /* Step 1: Get HMM Forward/Backward DP matrices.
   * Step 2: F/B       -> HMM bands.
   * Step 3: HMM bands -> CM bands.
   */

  /* Step 1: Get HMM Forward/Backward DP matrices. */
  CP9_dpmatrix_t *cp9_fwd;       /* growable DP matrix for forward                       */
  CP9_dpmatrix_t *cp9_bck;       /* growable DP matrix for backward                      */
  int do_scan2bands;             /* TRUE to use scanning Forward/Backward to get posteriors */
  if(cm->align_opts & CM_ALIGN_TIME) {
  watch = esl_stopwatch_Create();
        esl_stopwatch_Start(watch);
  }
  do_scan2bands = (cm->search_opts & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;
  int be_safe = FALSE; /* TEMPORARY, pass this in after calcing it once in actually_align_targets() */
  //int be_safe = TRUE;
  sc = Xcp9_FastForward(cm, dsq, i0, j0, j0-i0+1, 0., NULL, NULL, NULL,
			do_scan2bands, /* are we using scanning Forward/Backward */
			TRUE,      /* we are going to use posteriors to align */
			FALSE,     /* we're not rescanning */
			FALSE,     /* don't be memory efficient */
			be_safe,   /* can we accelerate w/ no -INFTY logsum funcs? */
			&cp9_fwd); /* give the DP matrix back */
  
  if(debug_level >= 0) printf("CP9 Forward  score : %.4f\n", sc);
  if(cm->align_opts & CM_ALIGN_TIME) {
      esl_stopwatch_Stop(watch);
      esl_stopwatch_Display(stdout, watch, "CP9 bands Forward CPU time: ");
      esl_stopwatch_Start(watch);
  }
  sc = Xcp9_FastBackward(cm, dsq, i0, j0, (j0-i0+1), 
			 0,    /* cp9_cutoff score, irrelevant */
			 NULL,  /* don't care about score of each posn */
			 NULL,  /* don't care about best scoring start point */
			 NULL,  /* don't report hits to results data structure */
			 do_scan2bands, /* are we using scanning Forward/Backward */
			 TRUE,  /* we are going to use posteriors to align */
			 FALSE, /* we're not rescanning */
			 FALSE, /* don't be memory efficient */
			 &cp9_bck); /* give the DP matrix back */
  if(debug_level >= 0) printf("CP9 Backward score : %.4f\n", sc);
  if(cm->align_opts & CM_ALIGN_TIME) {
      esl_stopwatch_Stop(watch);
      esl_stopwatch_Display(stdout, watch, "CP9 bands Backward CPU time: ");
      esl_stopwatch_Start(watch);
  }

  /* Step 2: F/B -> HMM bands. */
  if(use_sums)
    Xcp9_FB2HMMBandsWithSums(cm->cp9, dsq, cp9_fwd, cp9_bck, cp9b, i0, j0, cp9b->hmm_M,
			     (1.-cm->tau), (cm->search_opts & CM_SEARCH_HMMSCANBANDS), debug_level);
  else
    Xcp9_FB2HMMBands(cm->cp9, dsq, cp9_fwd, cp9_bck, cp9b, i0, j0, cp9b->hmm_M,
		     (1.-cm->tau), (cm->search_opts & CM_SEARCH_HMMSCANBANDS), debug_level);
  
  if(cm->align_opts & CM_ALIGN_TIME) {
      esl_stopwatch_Stop(watch);
      esl_stopwatch_Display(stdout, watch, "CP9 band bounds calculation CPU time: ");
      esl_stopwatch_Destroy(watch);
  }

  if(debug_level > 0) debug_print_hmm_bands(stdout, j0, cp9b, cm->tau, 1);

  if(cm->align_opts & CM_ALIGN_TIME) {
  watch = esl_stopwatch_Create();
        esl_stopwatch_Start(watch);
  }
  /* Step 3: HMM bands  ->  CM bands. */
  hmm2ij_bands(cm, cm->cp9map, i0, j0, cp9b->pn_min_m, cp9b->pn_max_m, 
	       cp9b->pn_min_i, cp9b->pn_max_i, cp9b->pn_min_d, cp9b->pn_max_d, 
	       cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, debug_level);
  
  if(cm->align_opts & CM_ALIGN_TIME) {
      esl_stopwatch_Stop(watch);
      esl_stopwatch_Display(stdout, watch, "hmm2ij bands calculation CPU time: ");
      esl_stopwatch_Destroy(watch);
  }
  /* Use the CM bands on i and j to get bands on d, specific to j. */
  for(v = 0; v < cm->M; v++) {
    ESL_ALLOC(cp9b->hdmin[v], sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
    ESL_ALLOC(cp9b->hdmax[v], sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
  }
  ij2d_bands(cm, (j0-i0+1), cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax,
	     cp9b->hdmin, cp9b->hdmax, debug_level);
  
  if(debug_level > 0)
    PrintDPCellsSaved_jd(cm, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
			 (j0-i0+1));

  return;

 ERROR:
  cm_Fail("Memory allocation error.\n");
}

/*****************************************************************************
 * Function: Xcp9_FB2HMMBands()
 * Date:     EPN, 04.03.06
 *           EPN, Mon Oct 15 18:20:42 2007 [updated/optimized]
 *
 * Purpose: Determine the band on all HMM states given a Forward and
 *          Backward matrix. Do this by calculating and summing log posterior
 *          probabilities that each state emitted/was visited at each posn,
 *          starting at the sequence ends, and creeping in, until the half the
 *          maximum allowable probability excluded is reached on each side.
 *
 * CP9_t hmm        the HMM
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
 *****************************************************************************/
void
Xcp9_FB2HMMBands(CP9_t *hmm, ESL_DSQ *dsq, CP9_dpmatrix_t *fmx, CP9_dpmatrix_t *bmx, CP9Bands_t *cp9b, 
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

  printf("in Xcp9_FB2HMMBands()\n");
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
    for (ip = 0; i <= L; i++) {
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
      for(k = 0; k <= M; k++)
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
	     
	  /*mass_m[k] = ILogsum(mass_m[k], post->mmx[ip][k]);
	    mass_i[k] = ILogsum(mass_i[k], post->imx[ip][k]);
	    mass_d[k] = ILogsum(mass_d[k], post->dmx[ip][k]);

	    if((! nset_m[k]) && (mass_m[k] > thresh)) { cp9b->pn_min_m[k] = i;   nset_m[k] = TRUE; }
	    if((! nset_i[k]) && (mass_i[k] > thresh)) { cp9b->pn_min_i[k] = i;   nset_i[k] = TRUE; }
	    if((! nset_d[k]) && (mass_d[k] > thresh)) { cp9b->pn_min_d[k] = i-1; nset_d[k] = TRUE; } */
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
	  ///mass_m[k] = ILogsum(mass_m[k], post->mmx[ip][k]);
	  ///mass_i[k] = ILogsum(mass_i[k], post->imx[ip][k]);
	  ///mass_d[k] = ILogsum(mass_d[k], post->dmx[ip][k]);

	  ///if((! xset_m[k]) && (mass_m[k] > thresh)) { cp9b->pn_min_m[k] = i;   xset_m[k] = TRUE; }
	  ///if((! xset_i[k]) && (mass_i[k] > thresh)) { cp9b->pn_min_i[k] = i;   xset_i[k] = TRUE; }
	  ///if((! xset_d[k]) && (mass_d[k] > thresh)) { cp9b->pn_min_d[k] = i-1; xset_d[k] = TRUE; }
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
      /* TEMPORARY this is how I used to do it with hmm_band_bounds() is this right?
       * maybe it should set max = min if min was set but not max and vice versa? 
       */
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

  if(debug_level > 0) debug_print_hmm_bands(stdout, j0, cp9b, (1.-p_thresh), 1);
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
  return;

 ERROR:
  cm_Fail("memory allocation error.");
}


/*****************************************************************************
 * Function: Xcp9_FB2HMMBandsWithSums()
 * Date:     EPN, Wed Oct 17 10:22:44 2007
 *
 * Purpose: Determine the band on all HMM states given a Forward and
 *          Backward matrix. Do this by calculating and summing log posterior
 *          probabilities that each state emitted/was visited at each posn,
 *          starting at the sequence ends, and creeping in, until the half the
 *          maximum allowable probability excluded is reached on each side.
 *
 * CP9_t hmm        the HMM
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
 *****************************************************************************/
void
Xcp9_FB2HMMBandsWithSums(CP9_t *hmm, ESL_DSQ *dsq, CP9_dpmatrix_t *fmx, CP9_dpmatrix_t *bmx, CP9Bands_t *cp9b, 
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
  
  printf("in Xcp9_FB2HMMBandsWithSums()\n");
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
  CP9Posterior(dsq, i0, j0, hmm, fmx, bmx, post, did_scan);

  /* fill ipost_sums in cp9bands data structure */
  CP9_ifill_post_sums(post, cp9b, i0, j0);

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
  assert(nset_m[0]);
  assert(nset_i[0]);
  assert(xset_m[0]);
  assert(xset_i[0]);
  /* D_0 state does not exist */
  for(k = 1; k <= M; k++)
    {
      assert(nset_m[k]);
      assert(nset_i[k]);
      assert(nset_d[k]);
      assert(xset_m[k]);
      assert(xset_i[k]);
      assert(xset_d[k]);
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

  if(debug_level > 0) debug_print_hmm_bands(stdout, j0, cp9b, (1.-p_thresh), 1);
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
  return;

 ERROR:
  cm_Fail("memory allocation error.");
}

/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_FASTALIGN_BENCHMARK
/* gcc -o benchmark-fastalign -g -O2 -I. -L. -I../easel -L../easel -DIMPL_FASTALIGN_BENCHMARK cm_fastalign.c -linfernal -leasel -lm
 * ./benchmark-fastalign <cmfile> <seqfile>
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_sqio.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "33", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-e",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "emit sequences from CM, don't randomly create them", 0 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "configure CM/HMM for local alignment", 0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of target seqs",                          0 },
  { "-L",        eslARG_INT,     NULL, NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs, default: consensus length", 0 },
  { "-o",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute original CYK HMM banded alignment implementation", 0 },
  { "--scan",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "run in scan mode, not alignment mode", 0 },
  { "--compare", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "compare HMM bands calc'ed with new optimized imp vs 0.8 imp", 0 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use posterior sums during HMM band calculation (widens bands)", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "benchmark driver for fast HMM banded CYK alignment and scanning algorithm";

int 
main(int argc, char **argv)
{
  int status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  CM_t           *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  int             i;
  float           sc, bsc;
  int             b;
  char           *cmfile = esl_opt_GetArg(go, 1);
  CMFILE         *cmfp;	    /* open input CM file stream */
  ESL_DSQ        *dsq;
  int             L;        /* length of sequence */
  int             do_random;
  int             N = esl_opt_GetInteger(go, "-N");
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, either randomly created, or emitted from CM (if -e) */
  int             do_align;
  
  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  do_align = TRUE;
  if(esl_opt_GetBoolean(go, "--scan")) do_align = FALSE; 

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  /* determine sequence length */
  if(esl_opt_IsDefault(go, "-L")) L = cm->clen;      
  else                            L = esl_opt_GetInteger(go, "-L");

  /* configure CM for HMM banded alignment */
  //cm->config_opts |= CM_CONFIG_QDB;
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  
  cm->align_opts  |= CM_ALIGN_TIME;
  if(do_align) { 
    cm->align_opts  |= CM_ALIGN_NOSMALL;
    cm->align_opts  |= CM_ALIGN_HBANDED;
    if(esl_opt_GetBoolean(go, "--sums")) cm->align_opts |= CM_ALIGN_SUMS;
  }
  else /* don't align, scan */ {
    cm->search_opts  |= CM_SEARCH_HBANDED;
    cm->search_opts  |= CM_SEARCH_HMMSCANBANDS;
    if(esl_opt_GetBoolean(go, "--sums")) cm->search_opts |= CM_SEARCH_SUMS;
  }
  if(esl_opt_GetBoolean(go, "-l")) { 
    cm->config_opts  |= CM_CONFIG_LOCAL;
    cm->config_opts  |= CM_CONFIG_HMMLOCAL;
    cm->config_opts  |= CM_CONFIG_HMMEL;
  }

  ConfigCM(cm, NULL, NULL);
  CP9Bands_t *cp9b;
  cp9b = AllocCP9Bands(cm, cm->cp9);

  /* get sequences */
  if(do_random) {
    double *dnull;
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);
    seqs_to_aln = RandomEmitSeqsToAln(r, cm->abc, dnull, 1, N, L, FALSE);
    free(dnull);
  }
  else /* don't randomly generate seqs, emit them from the CM */
    seqs_to_aln = CMEmitSeqsToAln(r, cm, 1, N, FALSE);

  for (i = 0; i < N; i++)
    {
      esl_stopwatch_Start(w);
      Xcp9_seq2bands(cm, seqs_to_aln->sq[i]->dsq, 1, seqs_to_aln->sq[i]->n, cp9b,  
		     NULL, /* we don't want posterior matrix back */
		     0);
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "Exptl Band calc", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");
      
      if(esl_opt_GetBoolean(go, "--compare")) {
	printf("\n\n");
	CP9Bands_t *cp9b_old;
	cp9b_old = AllocCP9Bands(cm, cm->cp9);
	CP9_seq2bands(cm, seqs_to_aln->sq[i]->dsq, 1, seqs_to_aln->sq[i]->n, cp9b_old,  
		      NULL, /* we don't want posterior matrix back */
		      0);
	esl_stopwatch_Stop(w);
	printf("%4d %-30s %17s", i+1, "Old Band calc", "");
	esl_stopwatch_Display(stdout, w, "CPU time: ");
	printf("\n\n");
	printf("Comparing bands...");
	cp9_compare_bands(cp9b, cp9b_old);
	printf("done. [passed]\n");
	FreeCP9Bands(cp9b_old);
      }
      esl_stopwatch_Start(w);
      /*esl_stopwatch_Start(w);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastInside_b_jd_me(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");*/
      /*sc = CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, seqs_to_aln->sq[i]->n, 0, 1, L, NULL, cp9b->jmin, 
			    cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);
      
      printf("%4d %-30s %10.4f bits ", (i+1), "CYKInside_b_jd(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");*/
      //      }      
      printf("\n");
    }

  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);

  printf("success\n");
  return 0;

 ERROR:
  cm_Fail("memory allocation error");
}
#endif /*IMPL_FASTALIGN_BENCHMARK*/

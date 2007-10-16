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
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

#define TSC(s,k) (tsc[(v) * MAXCONNECT + (s)])
#define AMX(j,v,d) (alphap[(j * cm->M * (W+1)) + ((v) * (W+1) + d)])


/*****************************************************************************/
/* Functions to determine HMM bands 
 * CP9_hmm_band_bounds()
 */
/*****************************************************************************
 * EPN, Mon Oct 15 18:20:42 2007
 * originally started: EPN 04.03.06
 * Function: cp9_EXPTL_hmm_band_bounds()
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
 * int i0           start of target subsequence (often 1, beginning of dsq)
 * int j0           end of target subsequence (often L, end of dsq)
 * int   M          number of nodes in HMM (num columns of post matrix)
 * CP9Bands_t cp9b  CP9 bands data structure
 * double p_thresh  the probability mass we're requiring is within each band
 * int use_sums     TRUE to use posterior sums instead of raw posteriors 
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *****************************************************************************/
void
cp9_EXPTL_S_hmm_band_bounds(int **post, int i0, int j0, int M, CP9Bands_t *cp9b, 
			    double p_thresh, int use_sums, int debug_level)
{
  int k;         /* counter over nodes of the model */
  int lmass_exc; /* the log of the probability mass currently excluded on the left*/
  int rmass_exc; /* the log of the probability mass currently excluded on the right*/
  int log_p_side; /* the log probability we're allowed to exclude on each side */
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

  printf("in cp9_EXPTL_hmm_band_bounds()\n");

  log_p_side = Prob2Score(((1. - p_thresh)/2.), 1.); /* allowable prob mass excluded on each side */

  /* step through each node */
  for(k = 0; k <= M; k++)
    {
      curr_log_p_side_m = curr_log_p_side_i = curr_log_p_side_d = log_p_side; 
      if(use_sums) {
	curr_log_p_side_m += cp9b->isum_pn_m[k]; /* if we use sums strategy, normalize
						  * so total prob of entering k = 1. */
	curr_log_p_side_i += cp9b->isum_pn_i[k];
	curr_log_p_side_d += cp9b->isum_pn_d[k];

      argmax_pn_m = i0;
      argmax_pn_i = i0;
      argmax_pn_d = i0;
      max_post_m = post->mmx[1][k]; /* post[1][k] corresponds to i0 for node k */
      max_post_i = post->imx[1][k]; /* post[1][k] corresponds to i0 for node k */
      max_post_d = post->dmx[1][k]; /* post[1][k] corresponds to i0 for node k */
      cp9b->pn_min_m[k] = cp9b->pn_min_i[k] = cp9b->pn_min_d[k] = i0; 
      cp9b->pn_max_m[k] = cp9b->pn_max_i[k] = cp9b->pn_max_d[k] = j0-1;
      /* i' = i - i0 + 1; i' is offset index for first dimension of post[][] */
      lmass_exc_m = post->mmx[(cp9b->pn_min_m[k]-1)-i0+1][k];
      rmass_exc_m = post->mmx[(cp9b->pn_max_m[k]+1)-i0+1][k];
      /*creep in on the left, until we exceed our allowable prob mass to exclude.*/ 
      while(pn_min[k] <= j0 && lmass_exc <= (curr_log_p_side))
	{
	  //if(post[(pn_min[k])-i0+1][k] > max_post) /* save info on most likely posn 
	  //* in case whole seq is outside band */
	  //{
	  //max_post = post[(pn_min[k])-i0+1][k];
	  //argmax_pn = pn_min[k];
	  //}
	  lmass_exc = ILogsum(lmass_exc, post[(pn_min[k])-i0+1][k]);
	  pn_min[k]++;
	}
      /* we went one posn too far, back up*/
      pn_min[k]--;
      
      /*creep in on the right, until we exceed our allowable prob mass to exclude.*/
      while(pn_max[k] >= (i0-1) && rmass_exc <= (curr_log_p_side))
	{
	  rmass_exc = ILogsum(rmass_exc, post[(pn_max[k]-i0+1)][k]);
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
      if(state_type == HMMDELETE || (k == 0 && state_type == HMMMATCH))
	{
	  /* We have to deal with off by ones in the delete states (which includes
	   * the HMM state M_0 which is silent - mapping to ROOT_S of the CM)
	   * e.g. pn_min_d[k] = i, means posn i was last residue emitted
	   * prior to entering node k's delete state. However, for a CM,
	   * if a delete states sub-parsetree is bounded by i' and j', then
	   * positions i' and j' HAVE YET TO BE EMITTED.
	   */
	  pn_min[k]++;
	  pn_max[k]++;
	}
    }
}


/* Function: CP9EXPTLPosterior()
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
 * Args:     dsq      - sequence in digitized form
 *           i0       - start of target subsequence (often 1, beginning of dsq)
 *           j0       - end of target subsequence (often L, end of dsq)
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           mx       - pre-allocated dynamic programming matrix
 *           
 * Return:   void
 */
void
CP9EXPTL_S_Posterior(ESL_DSQ *dsq, int i0, int j0,
		     struct cplan9_s *hmm,
		     struct cp9_dpmatrix_s *fmx,
		     struct cp9_dpmatrix_s *bmx,
		     struct cp9_dpmatrix_s *mx,
		     double p_thresh,
		     CP9Bands_t *cp9b)
{
  int i;
  int status;
  int k;
  int sc;
  int W;		/* subsequence length */
  int ip;		/* i': relative position in the subsequence  */
  /*float temp_sc;*/
  int tmin; /* the log probability we're allowed to exclude on each side */
  int tmax;
  int *mxset, *ixset, *dxset, *mnset, *inset, *dnset;

  if(dsq == NULL)
    cm_Fail("in CP9EXPTLPosterior(), dsq is NULL.");
  printf("in CP9EXPTLPosterior()\n");

  W  = j0-i0+1;		/* the length of the subsequence */

  sc = bmx->mmx[0][0];

  tmin = Prob2Score((1. - p_thresh)/2., 1.); /*      allowable prob mass excluded on each side, for mins  */
  tmax = Prob2Score(      p_thresh /2., 1.); /* 1. - allowable prob mass excluded on each side, for maxes */
  printf("tmin: %d\ntmax: %d\n", tmin, tmax);
  printf("p_thresh: %f\n", p_thresh);
  ESL_ALLOC(mxset, sizeof(int) * (hmm->M+1));
  ESL_ALLOC(ixset, sizeof(int) * (hmm->M+1));
  ESL_ALLOC(dxset, sizeof(int) * (hmm->M+1));
  ESL_ALLOC(mnset, sizeof(int) * (hmm->M+1));
  ESL_ALLOC(inset, sizeof(int) * (hmm->M+1));
  ESL_ALLOC(dnset, sizeof(int) * (hmm->M+1));
  esl_vec_ISet(mxset, hmm->M+1, FALSE);
  esl_vec_ISet(ixset, hmm->M+1, FALSE);
  esl_vec_ISet(dxset, hmm->M+1, FALSE);
  esl_vec_ISet(mnset, hmm->M+1, FALSE);
  esl_vec_ISet(inset, hmm->M+1, FALSE);
  esl_vec_ISet(dnset, hmm->M+1, FALSE);
  
  /* note boundary conditions, case by case by case... */
  mx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc; /* fmx->mmx[0][0] is 0, bmx->mmx[1][0] is overall score */
  mx->imx[0][0] = -INFTY; /*need seq to get here*/
  mx->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  for (k = 1; k <= hmm->M; k++) 
    {
      mx->mmx[0][k] = -INFTY; /*need seq to get here*/
      mx->imx[0][k] = -INFTY; /*need seq to get here*/
      mx->dmx[0][k] = mx->dmx[0][k-1], fmx->dmx[0][k] + bmx->dmx[0][k] - sc;
      if(! mnset[k] && mx->mmx[0][k] > tmin) { cp9b->pn_min_m[k] = 0; mnset[k] = TRUE; }
      if(! inset[k] && mx->imx[0][k] > tmin) { cp9b->pn_min_i[k] = 0; inset[k] = TRUE; }
      if(! dnset[k] && mx->dmx[0][k] > tmin) { cp9b->pn_min_d[k] = 0; dnset[k] = TRUE; }
      if(! mxset[k] && mx->mmx[0][k] > tmax) { cp9b->pn_max_m[k] = 0; mxset[k] = TRUE; }
      if(! ixset[k] && mx->imx[0][k] > tmax) { cp9b->pn_max_i[k] = 0; ixset[k] = TRUE; }
      if(! dxset[k] && mx->dmx[0][k] > tmax) { cp9b->pn_max_d[k] = 0; dxset[k] = TRUE; }
    }
      
  for (ip = 1; ip <= W; ip++) /* ip is the relative position in the seq */
    {
      int temp_sc;
      i = i0+ip-1;		/* e.g. i is actual index in dsq, runs from i0 to j0 */
      mx->mmx[ip][0] = -INFTY; /*M_0 does not emit*/
      mx->imx[ip][0] = ILogsum(mx->imx[ip-1][0], fmx->imx[ip][0] + bmx->imx[ip][0] - hmm->isc[dsq[i]][0] - sc);
      /*hmm->isc[dsq[i]][0] will have been counted in both fmx->imx and bmx->imx*/
      mx->dmx[ip][0] = -INFTY; /*D_0 does not exist*/

      /*printf("fmx->mmx[ip:%d][0]: %d\n bmx->mmx[ip:%d][0]: %d\n", ip, fmx->mmx[ip][0], ip, bmx->mmx[ip][0]);
	printf("fmx->imx[ip:%d][0]: %d\n bmx->imx[ip:%d][0]: %d\n", ip, fmx->imx[ip][0], ip, bmx->imx[ip][0]);
	printf("fmx->dmx[ip:%d][0]: %d\n bmx->dmx[ip:%d][0]: %d\n", ip, fmx->dmx[ip][0], ip, bmx->dmx[ip][0]);*/
      for (k = 1; k <= hmm->M; k++) 
	{
	  temp_sc         = fmx->mmx[ip][k] + bmx->mmx[ip][k] - hmm->msc[dsq[i]][k] - sc;
	  mx->mmx[ip][k]  = ILogsum(mx->mmx[ip-1][k], temp_sc);
	  fmx->mmx[ip][k] = temp_sc;
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  temp_sc         = fmx->imx[ip][k] + bmx->imx[ip][k] - hmm->isc[dsq[i]][k] - sc;
	  mx->imx[ip][k]  = ILogsum(mx->imx[ip-1][k], temp_sc);
	  fmx->imx[ip][k] = temp_sc;
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
	  temp_sc         = fmx->dmx[ip][k] + bmx->dmx[ip][k] - sc;
	  mx->dmx[ip][k]  = ILogsum(mx->dmx[ip-1][k], temp_sc);
	  fmx->dmx[ip][k] = temp_sc;

	  if(! mnset[k] && mx->mmx[ip][k] > tmin) { cp9b->pn_min_m[k] = ip; mnset[k] = TRUE; }
	  if(! inset[k] && mx->imx[ip][k] > tmin) { cp9b->pn_min_i[k] = ip; inset[k] = TRUE; }
	  if(! dnset[k] && mx->dmx[ip][k] > tmin) { cp9b->pn_min_d[k] = ip-1; dnset[k] = TRUE; }
	  if(! mxset[k] && mx->mmx[ip][k] > tmax) { cp9b->pn_max_m[k] = ip; mxset[k] = TRUE; }
	  if(! ixset[k] && mx->imx[ip][k] > tmax) { cp9b->pn_max_i[k] = ip; ixset[k] = TRUE; }
	  if(! dxset[k] && mx->dmx[ip][k] > tmax) { cp9b->pn_max_d[k] = ip-1; dxset[k] = TRUE; }
	  /*printf("fmx->mmx[ip:%d][%d]: %d\n bmx->mmx[ip:%d][%d]: %d\n", ip, k, fmx->mmx[ip][k], ip, k, bmx->mmx[ip][k]);
	  printf("fmx->imx[ip:%d][%d]: %d\n bmx->imx[ip:%d][%d]: %d\n", ip, k, fmx->imx[ip][k], ip, k, bmx->imx[ip][k]);
	  printf("fmx->dmx[ip:%d][%d]: %d\n bmx->dmx[ip:%d][%d]: %d\n\n", ip, k, fmx->dmx[ip][k], ip, k, bmx->dmx[ip][k]);*/
	}	  
    }
  /* Some states may still have min/max unfilled, this is true if whole sequence is outside band for
   * that state (meaning that state is very unlikely to be in optimal parse). Current strategy is
   * to set such states to band of width 1, the most likely position that state accounts for.
   */
  for(k = 0; k <= hmm->M; k++) {
    int best_pos;
    int best_sc;
    if(mnset[k]) assert(mxset[k]);
    if(mxset[k]) assert(mnset[k]);
    if(inset[k]) assert(ixset[k]);
    if(ixset[k]) assert(inset[k]);
    if(dnset[k]) assert(dxset[k]);
    if(dxset[k]) assert(dnset[k]);

    if(! mnset[k]) { 
      assert(! mxset[k]); 
      best_pos = 0;
      best_sc = fmx->mmx[0][k];
      for(ip = 1; ip <= W; ip++) {
	if(fmx->mmx[ip][k] > best_sc) {
	  best_pos = ip;
	  best_sc  = fmx->mmx[ip][k];
	}
      }
      cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = best_pos;
    }

    if(! inset[k]) { 
      assert(! ixset[k]); 
      best_pos = 0;
      best_sc = fmx->imx[0][k];
      for(ip = 1; ip <= W; ip++) {
	if(fmx->imx[ip][k] > best_sc) {
	  best_pos = ip;
	  best_sc  = fmx->imx[ip][k];
	}
      }
      cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = best_pos;
    }

    if(! dnset[k]) { 
      assert(! dxset[k]); 
      best_pos = 0;
      best_sc = fmx->dmx[1][k];
      for(ip = 2; ip <= W; ip++) { /* note we start at 2, b/c we have to subtract 1 for deletes */
	if(fmx->dmx[ip][k] > best_sc) {
	  best_pos = ip - 1;
	  best_sc  = fmx->dmx[ip-1][k];
	}
      }
      cp9b->pn_min_d[k] = cp9b->pn_max_d[k] = best_pos;
    }
  }
  /* We have to deal with off by ones in the delete states (which includes
   * the HMM state M_0 which is silent - mapping to ROOT_S of the CM)
   * e.g. pn_min_d[k] = i, means posn i was last residue emitted
   * prior to entering node k's delete state. However, for a CM,
   * if a delete states sub-parsetree is bounded by i' and j', then
   * positions i' and j' HAVE YET TO BE EMITTED.
   */
  //cp9b->pn_min_m[0] = 1;
  //cp9b->pn_max_m[0] = 1;
  /*  float temp_sc;
      for(i = 0; i <= W; i++)
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
  free(mnset);
  free(inset);
  free(dnset);
  free(mxset);
  free(ixset);
  free(dxset);
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
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
  { "-o",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute original, non-hybrid jd banded CYK alignment implementation", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "benchmark driver for fast HMM banded CYK alignment algorithm";

int 
main(int argc, char **argv)
{
  int status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  CM_t            *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_ALPHABET   *abc     = NULL;
  int             i;
  float           sc, bsc;
  int             b;
  char            *cmfile = esl_opt_GetArg(go, 1);
  char            *sqfile = esl_opt_GetArg(go, 2);
  ESL_SQFILE      *sqfp;     /* open sequence input file stream */
  CMFILE          *cmfp;     /* open input CM file stream */
  seqs_to_aln_t   *seqs_to_aln;  /* sequences to align, holds seqs, parsetrees, CP9 traces, postcodes */
  int              L;        /* length of sequence */

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  /* open input sequence file */
  status = esl_sqfile_Open(sqfile, eslSQFILE_UNKNOWN, NULL, &(sqfp));
  if (status == eslENOTFOUND)    cm_Fail("File %s doesn't exist or is not readable\n", sqfile);
  else if (status == eslEFORMAT) cm_Fail("Couldn't determine format of sequence file %s\n", sqfile);
  else if (status == eslEINVAL)  cm_Fail("Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      cm_Fail("Sequence file open failed with error %d\n", status);

  /* configure CM for HMM banded alignment */
  //cm->config_opts |= CM_CONFIG_QDB;
  cm->align_opts  |= CM_ALIGN_NOSMALL;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  cm->align_opts  |= CM_ALIGN_TIME;
  //cm->config_opts  |= CM_CONFIG_LOCAL;
  //cm->config_opts  |= CM_CONFIG_HMMLOCAL;
  ConfigCM(cm, NULL, NULL);

  /* read in all sequences, this is wasteful, but Parsetrees2Alignment() requires all seqs in memory */
  seqs_to_aln = CreateSeqsToAln(100, FALSE);
  if((status = ReadSeqsToAln(cm->abc, sqfp, 0, TRUE, seqs_to_aln, FALSE)) != eslEOF) cm_Fail("Error reading sqfile: %s\n", sqfile);
  
  CP9Bands_t *cp9b;
  cp9b = AllocCP9Bands(cm, cm->cp9);
  for (i = 0; i < seqs_to_aln->nseq; i++)
    {
      esl_stopwatch_Start(w);
      L = seqs_to_aln->sq[i]->n;
      CP9_seq2bands(cm, seqs_to_aln->sq[i]->dsq, 1, L, cp9b,  
		    NULL, /* we don't want posterior matrix back */
		    0);
      esl_stopwatch_Stop(w);
      printf("%4d %-30s %17s", i+1, "Band calc", "");
      esl_stopwatch_Display(stdout, w, "CPU time: ");

      /*esl_stopwatch_Start(w);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastInside_b_jd_me(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");*/
      
      //if(esl_opt_GetBoolean(go, "-o")) {
      esl_stopwatch_Start(w);

	sc = CYKInside_b_jd(cm, seqs_to_aln->sq[i]->dsq, L, 0, 1, L, NULL, cp9b->jmin, 
			    cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);

	printf("%4d %-30s %10.4f bits ", (i+1), "CYKInside_b_jd(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	//      }      
	printf("\n");
    }

  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IMPL_FASTALIGN_BENCHMARK*/

/* CP9scan.c 
 * EPN
 * 
 * Scanning algorithms for CM Plan 9 HMMs.  These algorithms align
 * subsequences of the target sequence to the model (e.g. glocal or
 * local alignment) Global alignment algorithms are in hmmband.c.
 *
 * These functions are still under development and are very
 * experimental. No guarantees here.
 *
 *################################################################
 * CP9ForwardScan()      - Scan input sequence for high scoring
 *                         forward hits to the model.
 *################################################################
 * 
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

#include "stopwatch.h"          /* squid's process timing module        */
#include "hmmband.h"
#include "cm_wrappers.h"
  
/***********************************************************************
 * Function: CP9ForwardScan()
 * 
 * Purpose:  A Forward dynamic programming algorithm that scans an input sequence.
 *           Only allocates 2 rows of the Forward matrix.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 * 
 * Args:     
 *           cm        - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq       - sequence in digitized form
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - the maximum size of a hit (often cm->W)
 *           cutoff    - minimum score to report
 *           ret_isc   - RETURN: int log odds Forward score for each end point [0..(j0-i0+1)]
 *           ret_nhits - RETURN: number of hits above cutoff saved.
 *           ret_hitj  - RETURN: end positions of hits, 0..ret_nhits-1
 *           results   - scan_results_t to add to; if NULL, don't keep results
 *
 * Returns:  best_sc, score of maximally scoring end position j 
 */
float
CP9ForwardScan(CM_t *cm, char *dsq, int i0, int j0, int W, float cutoff, int **ret_isc, 
	       int **ret_hitj, int *ret_nhits, scan_results_t *results)
{
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /*     j-W: position in the subsequence                         */
  int          ip;          /* i': relative position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  int          L;           /* j0-i0+1: subsequence length                                  */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int         *isc;         /* prob (seq from j0..jp | HMM) [0..jp..cm->cp9->M]             */
  float        fsc;         /* float log odds score                                         */
  float        curr_sc;     /* temporary score used for filling in gamma                    */
  float       *gamma;       /* SHMM DP matrix for optimum nonoverlap resolution [0..j0-i0+1]*/
  int         *gback;       /* traceback pointers for SHMM                      [0..j0-i0+1]*/ 
  float       *savesc;      /* saves score of hit added to best parse at j      [0..j0-i0+1]*/ 
  int          nhits;       /* number of hits found w/bit score >= cutoff                   */
  int         *hitj;        /* end positions (j) of hits [0..nhits-1]                       */
  int          alloc_nhits; /* used to grow the hitj array                                  */
  float        best_sc;     /* Best overall score from semi-HMM to return                   */
  float        best_negsc;  /* Best score overall score to return, used if all scores < 0.  */

  if(cm->cp9 == NULL)
    Die("ERROR in CP9ForwardScan, but cm->cp9 is NULL.\n");

  /*printf("CP9 Forward memory  :   %8.2f MB\n", CP9ForwardScanRequires(cm->cp9, (j0-i0+1), W));*/

  best_sc     = IMPOSSIBLE;
  best_negsc  = IMPOSSIBLE;
  L = j0-i0+1;

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+1));
  gamma[0] = 0.;
  gback    = MallocOrDie(sizeof(int)   * (L+1));
  gback[0] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+1));

  /* Allocate a DP matrix with 2 rows, 0..M 
   */ 
  mmx   = MallocOrDie(sizeof(int *) * 2);
  imx   = MallocOrDie(sizeof(int *) * 2);
  dmx   = MallocOrDie(sizeof(int *) * 2);
  for(j = 0; j < 2; j++)
    {
      mmx[j]   = MallocOrDie(sizeof(int) * (cm->cp9->M+1));
      imx[j]   = MallocOrDie(sizeof(int) * (cm->cp9->M+1));
      dmx[j]   = MallocOrDie(sizeof(int) * (cm->cp9->M+1));
    }
  /* isc will hold P(seq up to j | Model) in int log odds form */
  isc   = MallocOrDie(sizeof(int) * (j0-i0+2));
			
  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= cm->cp9->M; k++)
    mmx[0][k] = imx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= cm->cp9->M; k++)
    dmx[0][k] = ILogsum(ILogsum(mmx[0][k-1] + cm->cp9->tsc[CTMD][k-1],
				imx[0][k-1] + cm->cp9->tsc[CTID][k-1]),
			dmx[0][k-1] + cm->cp9->tsc[CTDD][k-1]);

  /* We can do a full parse through all delete states. */
  isc[0] = dmx[0][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]; 

  /*****************************************************************
   * The main loop: scan the sequence from position i0 to j0.
   *****************************************************************/
  /* Recursion. */
  for (j = i0; j <= j0; j++)
    {
      jp = j-i0+1;     /* jp is relative position in the sequence */
      cur = (j-i0+1)%2;
      prv = (j-i0)%2;

      mmx[cur][0] = 0;  /* This is the 1 difference between a Forward scanner and the 
			 * regular Forward (CP9Forward), in CP9Forward, this cell is 
			 * set to -INFTY b/c we have to start at posn 0, 
			 * but here, we can start at any position in the seq.
			 */
      dmx[cur][0] = -INFTY;  /*D_0 is non-existent*/
      imx[cur][0]  = ILogsum(ILogsum(mmx[prv][0] + cm->cp9->tsc[CTMI][0],
				     imx[prv][0] + cm->cp9->tsc[CTII][0]),
			     dmx[prv][0] + cm->cp9->tsc[CTDI][0]);
      imx[cur][0] += cm->cp9->isc[(int) dsq[j]][0];
      
      for (k = 1; k <= cm->cp9->M; k++)
	{
	  mmx[cur][k]  = ILogsum(ILogsum(mmx[prv][k-1] + cm->cp9->tsc[CTMM][k-1],
				       imx[prv][k-1] + cm->cp9->tsc[CTIM][k-1]),
			       ILogsum(mmx[prv][0] + cm->cp9->bsc[k],
				       dmx[prv][k-1] + cm->cp9->tsc[CTDM][k-1]));
	  mmx[cur][k] += cm->cp9->msc[(int) dsq[j]][k];
	  
	  dmx[cur][k]  = ILogsum(ILogsum(mmx[cur][k-1] + cm->cp9->tsc[CTMD][k-1],
					imx[cur][k-1] + cm->cp9->tsc[CTID][k-1]),
				dmx[cur][k-1] + cm->cp9->tsc[CTDD][k-1]);
	  
	  imx[cur][k]  = ILogsum(ILogsum(mmx[prv][k] + cm->cp9->tsc[CTMI][k],
				       imx[prv][k] + cm->cp9->tsc[CTII][k]),
			       dmx[prv][k] + cm->cp9->tsc[CTDI][k]);
	  imx[cur][k] += cm->cp9->isc[(int) dsq[j]][k];
	  /*printf("mmx[jp:%d][%d]: %d\n", jp, k, mmx[cur][k]);
	    printf("imx[jp:%d][%d]: %d\n", jp, k, imx[cur][k]);
	    printf("dmx[jp:%d][%d]: %d\n", jp, k, dmx[cur][k]);*/
	}
      /* determine isc, the int score of all possible parses ending at the current
       * position (j) of the target sequence. */
      isc[jp] = -INFTY;
      for (k = 1; k <= cm->cp9->M; k++)
	isc[jp] = ILogsum(isc[jp], mmx[cur][k] + cm->cp9->esc[k]);
      isc[jp] = ILogsum(isc[jp], dmx[cur][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]); 
      isc[jp] = ILogsum(isc[jp], imx[cur][cm->cp9->M] + cm->cp9->tsc[CTIM][cm->cp9->M]); 
      /* transition from D_M -> end */
      fsc = Scorify(isc[jp]);

      if(fsc > best_negsc) 
	best_negsc = fsc;

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[jp]  = gamma[jp-1] + 0; /* extend without adding a new hit */
      gback[jp]  = -1;
      savesc[jp] = IMPOSSIBLE;
      i = ((j-W+1)> i0) ? (j-W+1) : i0;
      ip = i-i0+1;
      curr_sc = gamma[ip-1] + fsc + cm->cp9_sc_boost;
      /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
       * value is 0.0 if technique not used. */
      if (curr_sc > gamma[jp])
	{
	  gamma[jp]  = curr_sc;
	  gback[jp]  = i;
	  savesc[jp] = fsc;
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one and report them.
   *****************************************************************/ 
  alloc_nhits = 10;
  hitj        = MallocOrDie(sizeof(int)   * alloc_nhits);
  nhits       = 0;
  j           = j0;
  while (j > i0) {
    jp = j-i0+1;
    if (gback[jp] == -1) /* no hit */
      j--; 
    else                /* a hit, a palpable hit */
      {
	if(savesc[jp] > best_sc) 
	  best_sc = savesc[jp];
	if(savesc[jp] >= cutoff)
	{
	  if(results != NULL) /* report the hit */
	    {
	      report_hit(gback[jp], j, -1, savesc[jp], results); 
	      /* -1 is for saver, which is irrelevant for HMM hits */
	    }
	  hitj[nhits] = j;
	  nhits++;
	  if (nhits == alloc_nhits) 
	    {
	      alloc_nhits += 10;
	      hitj = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits));
	    }
	}
	j = gback[jp]-1;
      }
  }

  /* clean up and exit */
  free(gback);
  free(gamma);
  free(savesc);

  /*printf("returning from CP9ForwardScan()\n");*/
  if(ret_isc != NULL) *ret_isc = isc;
  else free(isc);
  if(ret_hitj != NULL) *ret_hitj = hitj;
  else free(hitj);
  if(ret_nhits != NULL) *ret_nhits = nhits;
  if(best_sc <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_sc = best_negsc;
  return best_sc;
}

/* Function: CP9ForwardScanRequires()  EPN 04.21.06
 * based on  CYKScanRequires() by:
 * Date:     SRE, Mon May  6 18:48:18 2002 [St. Louis]
 *
 * Purpose:  Return the memory required by CP9ForwardScan(), in megabytes.
 */
float
CP9ForwardScanRequires(CP9_t *hmm, int L, int W)
{
  float ram;

  ram =  (float) (2 * 3 * (sizeof(int) * (hmm->M+1))); /* the DP matrix */
  ram += (float) (6 * (sizeof (int *)));  /* DP matrix pointers */ 
  ram += (int)   (sizeof(int) * (L+1)); /* the isc scores */
  ram += (float) (sizeof(float) * (L+1)); /* gamma allocation */
  ram += (float) (sizeof(int)   * (L+1)); /* gback allocation */
  ram += (float) (sizeof(float) * (L+1)); /* savesc allocation */
  return (ram / 1000000.);
}

/***********************************************************************
 * Function: CP9FilteredScan()
 * Incept:   EPN, Tue Jan  9 06:28:49 2007
 * 
 * Purpose:  Scan a sequence with a CP9, filtering for promising subseqs
 *           that could be hits to the CM the CP9 was derived from.
 *           Pass filtered subseqs to actually_search_target() to be.
 *           scanned with a CM scan. 
 * 
 * Args:     
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report (or keep if filtering)
 *           results    - scan_results_t to add to, only passed to 
 *                        actually_search_target()
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc, score of maximally scoring end position j 
 */
float
CP9FilteredScan(CM_t *cm, char *dsq, int i0, int j0, int W, float cm_cutoff, 
		float cp9_cutoff, scan_results_t *results, int *ret_flen)
{
  printf("in CP9FilteredScan()\n");
  int *hitj;
  int  nhits;
  int h;
  int curr_i0, curr_j0, next_i0, next_j0;
  float best_hmm_sc;
  float best_cm_sc;
  float cm_sc;
  int   flen;
  float ffrac;

  /* Scan the (sub)seq w/Forward, getting j end points of hits above cutoff */
  best_hmm_sc = CP9ForwardScan(cm, dsq, i0, j0, W, cp9_cutoff, NULL, &hitj, &nhits, NULL);
  /* Send promising subseqs to actually_search_target(): send subseq [j-W..j..j+W] 
   * for each hit endpoint j returned from CP9ForwardScan, or if there's overlap, send 
   * minimal subseq that encompasses all overlapping hits with W residue 'pad' on 
   * both sides.
   */
  best_cm_sc = IMPOSSIBLE;
  flen = 0;
  for(h = 0; h <= nhits-1; h++) 
    {
      curr_j0 = ((hitj[h] + W) <= j0)    ? (hitj[h] + W)     : j0;
      curr_i0 = ((curr_j0 - (2*W)) >= 1) ? (curr_j0 - (2*W)) : 1;
      if((h+1) != nhits)
	{	
	  next_j0 = ((hitj[h] + W) <= j0)    ? (hitj[h] + W)     : j0;
	  next_i0 = ((next_j0 - (2*W)) >= 1) ? (next_j0 - (2*W)) : 1;
	}
      else next_i0 = next_j0 = -1;
      printf("hit: %d j: %d j-W: %d j+W: %d\n", h, hitj[h], (hitj[h]-W), (hitj[h]+W));
      while(curr_i0 <= next_j0)
	{
	  curr_i0 = next_i0;
	  h++;
	  printf("\tsucked in hit: %d i0: %d j0: %d\n", h, curr_i0, curr_j0);
	  if((h+1) != nhits)
	    {	
	      next_j0 = ((hitj[h] + W) <= j0)    ? (hitj[h] + W)     : j0;
	      next_i0 = ((next_j0 - (2*W)) >= 1) ? (next_j0 - (2*W)) : 1;
	    }
	  else next_i0 = next_j0 = -1;
	}
      printf("calling actually_search_target: %d %d\n", curr_i0, curr_j0);
      cm_sc =
	actually_search_target(cm, dsq, curr_i0, curr_j0, cm_cutoff, cp9_cutoff,
			       results, /* keep results                                 */
			       FALSE,   /* don't filter, we already have                */
			       FALSE,   /* we're not building a histogram for CM stats  */
			       FALSE,   /* we're not building a histogram for CP9 stats */
			       NULL);   /* filter fraction N/A                          */
      flen += (curr_j0 - curr_i0 + 1);
      if(cm_sc > best_cm_sc) best_cm_sc = cm_sc;
    }

  free(hitj);
  if(flen == 0) ffrac = 100.;
  else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
  printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);
  if(ret_flen != NULL) *ret_flen = flen;
  return best_cm_sc;
}
/***********************************************************************
 * Function: CP9ForwardBackwardScan()
 * 
 * Purpose:  Runs the Forward dynamic programming algorithm that scans an 
 *           input subsequence (i0-j0). And then the Backward DP algorithm on the
 *           same sequence. Forward finds likely endpoints of hits, 
 *           and Backward finds likely startpoints of hits.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 *           This function is messy in that it does not check that
 *           the log-odds scores are -INFTY (-987654321) before
 *           adding them. This means that -INFTY is interpreted
 *           as a possible transition, with a very low probability
 *           (2^-987654321).
 * 
 * Args:     dsq    - sequence in digitized form
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - max d: max size of a hit
 *           hmm    - the model
 *           ret_fmx - the Forward scanning matrix
 *           ret_bmx - the Backward scanning matrix
 *           ret_nhits - RETURN: number of hits
 *           ret_hitr  - RETURN: start states of hits, 0..nhits-1
 *           ret_hiti  - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj  - RETURN: end positions of hits, 0..nhits-1
 *           ret_hitsc - RETURN: scores of hits, 0..nhits-1            
 *           min_thresh- minimum score to report (EPN via Alex Coventry 03.11.06)
 *           pad       - number of nucleotides to add on to start and end points
 *
 * Returns:  
 *           hiti, hitj, hitsc are allocated in a helper function; caller free's w/ free().
 */
float
CP9ForwardBackwardScan(char *dsq, int i0, int j0, int W, CP9_t *hmm, CP9_dpmatrix_t **ret_fmx,
		       CP9_dpmatrix_t **ret_bmx, int *ret_nhits, int **ret_hitr, int **ret_hiti, 
		       int **ret_hitj, float **ret_hitsc, float min_thresh, int pad)
{
  CP9_dpmatrix_t *fmx;
  CP9_dpmatrix_t *bmx;
  int **fmmx;
  int **fimx;
  int **fdmx;
  int **femx;
  int **bmmx;
  int **bimx;
  int **bdmx;
  int **bemx;
  int   i,j,k;
  float     sc;     
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */

  int       alloc_nhits;	/* used to grow the hit arrays */

  int fwd_sum;                  /* calc'ed from the Forward matrix, the sum over all
				 * positions 1..L of all parses that start and end 
				 * anywhere from 1..L */
  int bck_sum;                  /* same as above, but calc'ed from the Backward matrix*/
  float fwd_sc;                 /* fwd_sum as a score */
  float bck_sc;                 /* bck_sum as a score */
  
  int fwd_nhits;                /* number of hits from Forward */
  int *fwd_hitr;                /* start states of hits from Forward, 0..nhits-1 (ALWAYS 0 as implemented)*/
  int *fwd_hiti;                /* start positions of hits from Forward, 0..nhits-1, (calc'ed as end points - W) */
  int *fwd_hitj;                /* end positions of hits from Forward, 0..nhits-1 */
  float *fwd_hitsc;             /* scores of hits, 0..nhits-1 */

  int bck_nhits;                /* number of hits from Backward */
  int *bck_hitr;                /* start states of hits from Backward, 0..nhits-1 (ALWAYS 0 as implemented)*/
  int *bck_hiti;                /* start positions of hits from Backward, 0..nhits-1 */
  int *bck_hitj;                /* end positions of hits from Backward, 0..nhits-1 (calc'ed as start points + W) */
  float *bck_hitsc;             /* scores of hits, 0..nhits-1 */
  int      L;                   /* j0-i0+1: subsequence length */
  int      jp;		        /* j': relative position in the subsequence  */
  int      ip;		        /* i': relative position in the subsequence  */


  fwd_sum = -INFTY;
  bck_sum = -INFTY;

  L = j0-i0+1;
  if (W > L) W = L; 

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. We use it first for
   * the Forward scan, and reuse it for the Backward scan.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+2));
  gamma[0] = 0.;
  gamma[L+1] = 0.;
  gback    = MallocOrDie(sizeof(int)   * (L+2));
  gback[0] = -1;
  gback[L+1] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+2));
  saver    = MallocOrDie(sizeof(int)   * (L+2));

  /* Allocate a DP matrix with 0..L rows, 0..M-1 
   */ 
  fmx = AllocCPlan9Matrix(L+1, hmm->M, &fmmx, &fimx, &fdmx, &femx);
  bmx = AllocCPlan9Matrix(L+1, hmm->M, &bmmx, &bimx, &bdmx, &bemx);

  /*****************************************************************
   * The main Forward loop: scan the sequence from position 1 to L.
   *****************************************************************/
  /* Initialization of the zero row.
   */
  fmmx[0][0] = 0;      /* M_0 is state B, we can start here */
  fimx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  fdmx[0][0] = -INFTY; /* D_0 doesn't exist. */

  /* Because there's a D state for every node 1..M, 
     fdmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= hmm->M; k++)
    fmmx[0][k] = fimx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= hmm->M; k++)
    fdmx[0][k] = ILogsum(ILogsum(fmmx[0][k-1] + hmm->tsc[CTMD][k-1],
				fimx[0][k-1] + hmm->tsc[CTID][k-1]),
			 fdmx[0][k-1] + hmm->tsc[CTDD][k-1]);
  femx[0][0] = fdmx[0][hmm->M] + hmm->tsc[CTDM][hmm->M]; 
  
  fwd_sum = ILogsum(fwd_sum, femx[0][0]);
  /* Recursion. Done as a pull.
   */
  for (j = i0; j <= j0; j++) 
    {
      jp = j-i0+1;     /* e.g. jp is relative position in j */
      fmmx[jp][0] = 0; /* This is the 1 difference between a Forward scanner and the 
		       * regular Forward (CP9Forward), 
		       * in CP9Forward, this cell is set to -INFTY.
		       * Here, we can start at any position in the seq,
		       * this isn't true in global alignment. 
		       */
      fdmx[jp][0] = -INFTY;  /*D_0 is non-existent*/
      fimx[jp][0]  = ILogsum(fmmx[jp-1][0] + hmm->tsc[CTMI][0],
			    fimx[jp-1][0] + hmm->tsc[CTII][0]);
      fimx[jp][0] += hmm->isc[(int) dsq[j]][0];
      
      for (k = 1; k <= hmm->M; k++)
	{
	  fmmx[jp][k]  = ILogsum(ILogsum(fmmx[jp-1][k-1] + hmm->tsc[CTMM][k-1],
				       fimx[jp-1][k-1] + hmm->tsc[CTIM][k-1]),
			       ILogsum(fmmx[jp-1][0] + hmm->bsc[k],
				       fdmx[jp-1][k-1] + hmm->tsc[CTDM][k-1]));
	  fmmx[jp][k] += hmm->msc[(int) dsq[j]][k];
	  
	  fdmx[jp][k]  = ILogsum(ILogsum(fmmx[jp][k-1] + hmm->tsc[CTMD][k-1],
				       fimx[jp][k-1] + hmm->tsc[CTID][k-1]),
			       fdmx[jp][k-1] + hmm->tsc[CTDD][k-1]);
	  
	  fimx[jp][k]  = ILogsum(ILogsum(fmmx[jp-1][k] + hmm->tsc[CTMI][k],
				       fimx[jp-1][k] + hmm->tsc[CTII][k]),
			       fdmx[jp-1][k] + hmm->tsc[CTDI][k]);
	  fimx[jp][k] += hmm->isc[(int) dsq[j]][k];
	  //printf("fmmx[%d][%d]: %d\n", i, k, fmmx[jp][k]);
	  //printf("fimx[%d][%d]: %d\n", i, k, fimx[jp][k]);
	  //printf("fdmx[%d][%d]: %d\n", i, k, fdmx[jp][k]);
	}
      
      femx[0][jp] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	{
	  femx[0][jp] = ILogsum(femx[0][jp], fmmx[jp][k] + hmm->esc[k]);
	  /*if(j == L) printf("k: %d | femx[0][L]: %d | fmmx[L][k]: %d | hmm->esc[k]: %d\n", k, femx[0][L], fmmx[jp][k], hmm->esc[k]);*/
	}
      femx[0][jp] = ILogsum(femx[0][jp], fdmx[jp][hmm->M] + hmm->tsc[CTDM][hmm->M]); 
      /*if(j == L) printf("added D trans | femx[0][L]: %d\n", femx[0][L]);*/
      /* transition from D_M -> end */
      femx[0][jp] = ILogsum(femx[0][jp], fimx[jp][hmm->M] + hmm->tsc[CTIM][hmm->M]); 
      /*if(j == L) printf("added I trans | femx[0][L]: %d\n", femx[0][L]);*/
      /* transition from I_M -> end */

      fwd_sc = Scorify(femx[0][jp]);
      /*printf("femx[0][%3d]: %10.5f\n", j, fwd_sc);*/

      fwd_sum = ILogsum(fwd_sum, femx[0][jp]);

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[jp]  = gamma[jp-1] + 0; /* extend without adding a new hit */
      gback[jp]  = -1;
      savesc[jp] = IMPOSSIBLE;
      saver[jp]  = -1;
      i = ((j-W+1)> i0) ? (j-W+1) : i0;
      sc = gamma[i-1] + Scorify(femx[0][jp]) - min_thresh;
      /*printf("j: %d | gamma[jp]: %f | sc: %f\n", j, gamma[jp], sc);*/
      if (sc > gamma[jp])
	{
	  gamma[jp]  = sc;
	  gback[jp]  = i;
	  savesc[jp] = Scorify(femx[0][jp]);
	  saver[jp]  = 0;
	  //saver[jp]  = bestr[d];
	}
    } /* end loop over end positions j */

  /*fwd_sc = Scorify(femx[0][L]);*/
  /*printf("Total Forward  score (summed over all posns): %f\n", Scorify(fwd_sum));*/
  fwd_sc = Scorify(fwd_sum); /* the total Forward score. */

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  alloc_nhits = 10;
  fwd_hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
  fwd_hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
  fwd_hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
  fwd_hitsc = MallocOrDie(sizeof(float) * alloc_nhits);

  j     = j0;
  fwd_nhits = 0;
  while (j >= i0) {
    jp = j-i0+1;
    if (gback[jp] == -1) /* no hit */
      j--; 
    else                /* a hit, a palpable hit */
      {
	fwd_hitr[fwd_nhits]   = saver[jp];
	fwd_hitj[fwd_nhits]   = j;
	fwd_hiti[fwd_nhits]   = gback[jp];
	fwd_hitsc[fwd_nhits]  = savesc[jp];
	fwd_nhits++;
	/*printf("fwd_nhits: %d\n", fwd_nhits);*/
	j = gback[jp]-1;
	
	if (fwd_nhits == alloc_nhits) {
	  fwd_hitr  = ReallocOrDie(fwd_hitr,  sizeof(int)   * (alloc_nhits + 10));
	  fwd_hitj  = ReallocOrDie(fwd_hitj,  sizeof(int)   * (alloc_nhits + 10));
	  fwd_hiti  = ReallocOrDie(fwd_hiti,  sizeof(int)   * (alloc_nhits + 10));
	  fwd_hitsc = ReallocOrDie(fwd_hitsc, sizeof(float) * (alloc_nhits + 10));
	  alloc_nhits += 10;
	}
      }
  }

  /*****************************************************************
   * The main Backward loop: scan the sequence from position L to 1.
   *****************************************************************/

  /* Initialization of the L row.
   */
  bemx[0][L] = 0; /*have to end in E*/

  bmmx[L][hmm->M] = bemx[0][L] + hmm->esc[hmm->M]; /* M<-E ...                   */
  bmmx[L][hmm->M] += hmm->msc[(int) dsq[j0]][hmm->M]; /* ... + emitted match symbol */
  bimx[L][hmm->M] = bemx[0][L] + hmm->tsc[CTIM][hmm->M];   /* I_M(C)<-E ... */
  bimx[L][hmm->M] += hmm->isc[(int) dsq[j0]][hmm->M];           /* ... + emitted match symbol */
  bdmx[L][hmm->M] = bemx[0][L] + hmm->tsc[CTDM][hmm->M];    /* D_M<-E */
  for (k = hmm->M-1; k >= 1; k--)
    {
      bmmx[L][k]  = bemx[0][L] + hmm->esc[k];
      bmmx[L][k]  = ILogsum(bmmx[L][k], bdmx[L][k+1] + hmm->tsc[CTMD][k]);
      bmmx[L][k] += hmm->msc[(int) dsq[j0]][k];

      bimx[L][k] = bdmx[L][k+1] + hmm->tsc[CTID][k];
      bimx[L][k] += hmm->isc[(int) dsq[j0]][k];

      bdmx[L][k] = bdmx[L][k+1] + hmm->tsc[CTDD][k];
    }
  
  bmmx[L][0] = -INFTY; /* no esc[0] b/c its impossible */
  bimx[L][0] = bdmx[L][1] + hmm->tsc[CTID][0];
  bimx[L][0] += hmm->isc[(int) dsq[j0]][hmm->M];    
  bdmx[L][0] = -INFTY; /*D_0 doesn't exist*/

  bck_sum = ILogsum(bck_sum, bmmx[L][0]);

  bck_sc = Scorify(bmmx[L][0]);
  /*printf("bmmx[L][  0]: %10.5f\n", bck_sc);*/
  /*printf("tmp_bck_sum L: %d\n", tmp_bck_sum);*/

  /* Recursion. Done as a pull.
   */
  for (i = j0-1; i >= i0; i--)
    {
      ip = i-i0 + 1;		/* ip is relative index in dsq (1 to L) */
      if(i > 0)
	{
	  bemx[0][ip] = 0; /* This is 1 of 3 differences between a Backward scanner and the 
			   * regular Backward (CP9Backward), 
			   * in CP9Backward, this cell is set to -INFTY.
			   * Here, we can end at the END state at any position,
			   * this isn't true in global alignment. 
			   * Second and third differences are related, 
			   * see the four comments 
			   * labelled "2nd diff" or "3rd diff" below
			   */
	  
	  /* Now the main states. Note the boundary conditions at M.
	   */
	  bmmx[ip][hmm->M] = ILogsum(bemx[0][i+1] + hmm->esc[hmm->M], /* M<-E ... 2nd diff w/non-scanner*/
				    bimx[ip+1][hmm->M] + hmm->tsc[CTMI][hmm->M]);
	  bmmx[ip][hmm->M] += hmm->msc[(int) dsq[i]][hmm->M];
	  bimx[ip][hmm->M] = ILogsum(bemx[0][ip] + hmm->tsc[CTIM][hmm->M],    /* I_M(C)<-E ... */
				    bimx[ip+1][hmm->M] + hmm->tsc[CTII][hmm->M]);
	  bimx[ip][hmm->M] += hmm->isc[(int) dsq[i]][hmm->M];
	  bdmx[ip][hmm->M] = ILogsum(bemx[0][ip] + hmm->tsc[CTDM][hmm->M], /* D_M<-E */
				    bimx[ip+1][hmm->M] + hmm->tsc[CTDI][hmm->M]);  
	  for (k = hmm->M-1; k >= 1; k--)
	    {
	      bmmx[ip][k]  = ILogsum(ILogsum(bemx[0][i+1] + hmm->esc[k], /* M<-E ... 2nd diff w/non-scanner*/
					    bmmx[ip+1][k+1] + hmm->tsc[CTMM][k]),
				    ILogsum(bimx[ip+1][k] + hmm->tsc[CTMI][k],
					    bdmx[ip][k+1] + hmm->tsc[CTMD][k]));	  
	      /*if(k == 1 && i == 1) printf("\tk: %d | bemx[0][i+1]: %d | hmm->esc[k]: %d\n\t\tbmmx[ip+1][k+1]: %d | hmm->tsc[CTMM][k]: %d\n\t\tbimx[ip+1][k]: %d | hmm->tsc[CTMI][k]: %d\n\t\tbdmx[ip][k+1]: %d | hmm->tsc[CTMD][k]: %d\n", k, bemx[0][(i+1)], hmm->esc[k], bmmx[(i+1)][(k+1)], hmm->tsc[CTMM][k], bimx[(i+1)][k], hmm->tsc[CTMI][k], bdmx[ip][(k+1)], hmm->tsc[CTMD][k]);*/
	      
	      bmmx[ip][k] += hmm->msc[(int) dsq[i]][k];
	      
	      bimx[ip][k]  = ILogsum(ILogsum(bmmx[ip+1][k+1] + hmm->tsc[CTIM][k],
					    bimx[ip+1][k] + hmm->tsc[CTII][k]),
				    bdmx[ip][k+1] + hmm->tsc[CTID][k]);
	      bimx[ip][k] += hmm->isc[(int) dsq[i]][k];
	      
	      bdmx[ip][k]  = ILogsum(ILogsum(bmmx[ip+1][k+1] + hmm->tsc[CTDM][k],
					    bimx[ip+1][k] + hmm->tsc[CTDI][k]),
				    bdmx[ip][k+1] + hmm->tsc[CTDD][k]);
	    }
	  
	  bimx[ip][0]  = ILogsum(ILogsum(bmmx[ip+1][1] + hmm->tsc[CTIM][0],
					bimx[ip+1][0] + hmm->tsc[CTII][0]),
				bdmx[ip][1] + hmm->tsc[CTID][0]);
	  bimx[ip][0] += hmm->isc[(int) dsq[i]][0];
	  bmmx[ip][0] = -INFTY;
	  /*for (k = hmm->M-1; k >= 1; k--)*/ /*M_0 is the B state, it doesn't emit*/
	  for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	    {
	      bmmx[ip][0] = ILogsum(bmmx[ip][0], bmmx[ip+1][k] + hmm->bsc[k]);
	      /*printf("hmm->bsc[%d]: %d | bmmx[ip+1][k]: %d\n", k, hmm->bsc[k], bmmx[(i+1)][k]);
		printf("bmmx[%3d][%3d]: %d | k: %3d\n", i, 0, bmmx[ip][0], k);*/
	    }
	  bmmx[ip][0] = ILogsum(bmmx[ip][0], bimx[ip+1][k] + hmm->tsc[CTMI][k]);
	  bmmx[ip][0] = ILogsum(bmmx[ip][0], bdmx[ip][k+1] + hmm->tsc[CTMD][k]);
	  
	  bdmx[ip][0] = -INFTY; /* D_0 does not exist */
	  
	}
      else if(i == 0)
	{
	  /* case when i = 0 */
	  bemx[0][0] = 0.; /* 3rd diff with non-scanner */
	  bmmx[0][hmm->M] = -INFTY; /* need seq to get here */
	  bimx[0][hmm->M] = -INFTY; /* need seq to get here */
	  bdmx[0][hmm->M] = ILogsum(bemx[0][0] + hmm->tsc[CTDM][hmm->M],  /*D_M<-E*/ /* 3rd diff with non-scanner */
				    bimx[1][hmm->M] + hmm->tsc[CTDI][hmm->M]);  
	  for (k = hmm->M-1; k >= 1; k--)
	    {
	      bmmx[0][k] = -INFTY; /* need seq to get here */
	      bimx[0][k] = -INFTY; /* need seq to get here */
	      bdmx[0][k]  = ILogsum(ILogsum(bmmx[1][k+1] + hmm->tsc[CTDM][k],
					    bimx[1][k] + hmm->tsc[CTDI][k]),
				    bdmx[0][k+1] + hmm->tsc[CTDD][k]);
	    }
	  bimx[0][0] = -INFTY; /*need seq to get here*/
	  bdmx[0][0] = -INFTY; /*D_0 doesn't exist*/
	  
	  bmmx[0][0] = -INFTY;
	  for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	    {
	      bmmx[0][0] = ILogsum(bmmx[0][0], bmmx[1][k] + hmm->bsc[k]);
	      /*printf("k: %d | bmmx[0][0]: %d | bmmx[1][k]: %d | hmm->bsc[k]: %d\n", k, bmmx[0][0], bmmx[1][k], hmm->bsc[k]);*/
	    }
	  bmmx[0][0] = ILogsum(bmmx[0][0], bdmx[0][1] + hmm->tsc[CTMD][0]);
	  /*printf("added D trans | bmmx[0][0]: %d \n", bmmx[0][0]);*/
	  bmmx[0][0] = ILogsum(bmmx[0][0], bimx[1][0] + hmm->tsc[CTMI][0]);
	  /*printf("added I trans | bmmx[0][0]: %d \n", bmmx[0][0]);*/
	  
	  bck_sc = Scorify(bmmx[0][0]);
	  /*printf("bmmx[0][  0]: %10.5f\n", bck_sc);*/
	  
	  bck_sum = ILogsum(bck_sum, bmmx[0][0]);
	  /*printf("Total Backward score (summed over all posns): %f\n", Scorify(bck_sum));		*/
	  bck_sc  = Scorify(bck_sum);  /* the total Backward score. */
	}
      /* The little semi-Markov model that deals with multihit parsing:
       */
      /* Here's a vicious off-by-one with HMMs and CMs:
       * an HMM parse can begin at position 0, i.e. backward[0][0] is the probability
       * of the entire sequence being emitted starting at the Begin state at position 0.
       * But a CM always has parses that start at position 1. So we have to translate
       * from HMM to CM in this way, and that's why there's a i+1 in the following loop.
       */
      gamma[ip+1]  = gamma[ip+2] + 0; /* extend without adding a new hit */
      /*printf("i: %d | gamma[i]: %f | gamma[i+1]: %f\n", i, gamma[i], gamma[i+1]);*/
      gback[ip+1]  = -1;
      savesc[ip+1] = IMPOSSIBLE;
      saver[ip+1]  = -1;
      j = ((i+1+W-1) < j0) ? (i+1+W-1) : j0;
      jp = j-i0+1;
      sc = gamma[jp+1] + Scorify(bmmx[ip][0]) - min_thresh;

      /*printf("bmmx[%3d][0]: %10.5f\n", i, Scorify(bmmx[ip][0]));*/
      bck_sum = ILogsum(bck_sum, bmmx[ip][0]);

      bck_sc = Scorify(bmmx[ip][0]);
      /*printf("bmmx[%3d][0]: %10.5f\n", i, bck_sc);*/
      /*printf("i: %d | bmmx[ip][0]: %f | j: %d | gamma[i]: %f | sc: %f\n", i, Scorify(bmmx[ip][0]), j, gamma[i], sc);*/
      if (sc > gamma[ip+1])
	{
	  gamma[ip+1]  = sc;
	  gback[ip+1]  = j;
	  savesc[ip+1] = Scorify(bmmx[ip][0]);
	  saver[ip+1]  = 0;
	  //saver[ip+1]  = bestr[d];
	}
    }
  /* End of Backward() code */

  /*****************************************************************
   * Traceback stage for Backward.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 

  alloc_nhits = 10;
  bck_hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
  bck_hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
  bck_hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
  bck_hitsc = MallocOrDie(sizeof(float) * alloc_nhits);

  i     = i0;
  bck_nhits = 0;
  ip = i-i0+1;
  while (i <= L) {
    if (gback[ip] == -1) /* no hit */
      { i++; ip++; }
    else                /* a hit, a palpable hit */
      {
	bck_hitr[bck_nhits]   = saver[ip];
	bck_hitj[bck_nhits]   = gback[ip];
	bck_hiti[bck_nhits]   = i;
	bck_hitsc[bck_nhits]  = savesc[ip];
	bck_nhits++;
	i = gback[ip]+1;
	ip = i-i0+1;
	if (bck_nhits == alloc_nhits) {
	  bck_hitr  = ReallocOrDie(bck_hitr,  sizeof(int)   * (alloc_nhits + 10));
	  bck_hitj  = ReallocOrDie(bck_hitj,  sizeof(int)   * (alloc_nhits + 10));
	  bck_hiti  = ReallocOrDie(bck_hiti,  sizeof(int)   * (alloc_nhits + 10));
	  bck_hitsc = ReallocOrDie(bck_hitsc, sizeof(float) * (alloc_nhits + 10));
	  alloc_nhits += 10;
	}
      }
  }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  /*debug_check_CP9_FBscan(fmx, bmx, hmm, Scorify(femx[0][L]), L, dsq);*/

  CP9_combine_FBscan_hits(i0, j0, W, fwd_nhits, fwd_hitr, fwd_hiti, fwd_hitj, fwd_hitsc, 
			        bck_nhits, bck_hitr, bck_hiti, bck_hitj, bck_hitsc,
			        ret_nhits, ret_hitr, ret_hiti, ret_hitj, ret_hitsc, pad);

  free(fwd_hitr);
  free(fwd_hitj);
  free(fwd_hiti);
  free(fwd_hitsc);

  free(bck_hitr);
  free(bck_hitj);
  free(bck_hiti);
  free(bck_hitsc);

  if (ret_fmx != NULL) *ret_fmx = fmx;
  else                FreeCPlan9Matrix(fmx);

  if (ret_bmx != NULL) *ret_bmx = bmx;
  else                FreeCPlan9Matrix(bmx);

  return fwd_sc;		/* the total Forward score. */
}


/* Function: CP9ScanFullPosterior()
 * based on Ian Holmes' hmmer/src/postprob.c::P7EmitterPosterior()
 *
 * Purpose:  Combines Forward and Backward scanning matrices into a posterior
 *           probability matrix. 
 *
 *           The main difference between this function and CP9FullPosterior()
 *           in hmmband.c is that this func takes matrices from CP9ForwardBackwardScan()
 *           in which parses are allowed to start and end in any residue.
 *           In CP9FullPosterior(), the matrices are calc'ed in CP9Forward()
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
 *           sequence (I think this is valid, but not sure (EPN)). 
 *           The last point distinguishes this function from P7EmitterPosterior() 
 *           which set all posterior values for for non-emitting states to -INFTY.
 *           The caller must allocate space for the matrix, although the
 *           backward matrix can be used instead (overwriting it will not
 *           compromise the algorithm).
 *           
 * Args:     dsq      - sequence in digitized form
 *           L        - length of sequence
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           mx       - pre-allocated dynamic programming matrix
 *           
 * Return:   void
 */
void
CP9ScanFullPosterior(char *dsq, int L,
		     CP9_t *hmm,
		     CP9_dpmatrix_t *fmx,
		     CP9_dpmatrix_t *bmx,
		     CP9_dpmatrix_t *mx)
{
  int i;
  int k;
  int fb_sum; /* tmp value, the probability that the current residue (i) was
	       * visited in any parse */
  fb_sum = -INFTY;
  for (i = 0; i <= L; i++) 
    {
      fb_sum = ILogsum(fb_sum, (fmx->emx[0][i]));
    }
  printf("fb_sc: %f\n", Scorify(fb_sum));
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
      
  for (i = 1; i <= L; i++)
    {
      /*fb_sum = -INFTY;*/ /* this will be probability of seeing residue i in any parse */
      /*for (k = 0; k <= hmm->M; k++) 
	{
	fb_sum = ILogsum(fb_sum, (fmx->mmx[i][k] + bmx->mmx[i][k] - hmm->msc[dsq[i]][k]));*/
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
      /*fb_sum = ILogsum(fb_sum, (fmx->imx[i][k] + bmx->imx[i][k] - hmm->isc[dsq[i]][k]));*/
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
      /*}*/
      mx->mmx[i][0] = -INFTY; /*M_0 does not emit*/
      mx->imx[i][0] = fmx->imx[i][0] + bmx->imx[i][0] - hmm->isc[(int) dsq[i]][0] - fb_sum;
      /*hmm->isc[dsq[i]][0] will have been counted in both fmx->imx and bmx->imx*/
      mx->dmx[i][0] = -INFTY; /*D_0 does not exist*/
      for (k = 1; k <= hmm->M; k++) 
	{
	  mx->mmx[i][k] = fmx->mmx[i][k] + bmx->mmx[i][k] - hmm->msc[(int) dsq[i]][k] - fb_sum;
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  mx->imx[i][k] = fmx->imx[i][k] + bmx->imx[i][k] - hmm->isc[(int) dsq[i]][k] - fb_sum;
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
	  mx->dmx[i][k] = fmx->dmx[i][k] + bmx->dmx[i][k] - fb_sum;
	}	  
    }

  /*  for(i = 0; i <= L; i++)
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
    }
  */
}

/***********************************************************************
 * Function: CP9_combine_FBscan_hit()
 * 
 * Purpose:  Given separate hits ((i, j, sc) triples) from a Forward scan
 *           and Backward scan of the same sequence, use the fact that
 *           the Forward scan probably got close to the correct END point (j),
 *           and the Backward scan probably got close to the correct START
 *           point (i) for each sequence to combine the hits from both 
 *           scans into (i, j, sc) triples we can then search with a CM.
 * 
 * 
 * Args:     i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - window length (maximum size of a hit considered)
 *           fwd_nhits - Forward: number of hits
 *           fwd_hitr  - Forward: start states of hits, 0..nhits-1
 *           fwd_hiti  - Forward: start positions of hits, 0..nhits-1
 *           fwd_hitj  - Forward: end positions of hits, 0..nhits-1
 *           fwd_hitsc - Forward: scores of hits, 0..nhits-1            
 *           bck_nhits - Backward: number of hits
 *           bck_hitr  - Backward: start states of hits, 0..nhits-1
 *           bck_hiti  - Backward: start positions of hits, 0..nhits-1
 *           bck_hitj  - Backward: end positions of hits, 0..nhits-1
 *           bck_hitsc - Backward: scores of hits, 0..nhits-1            
 *           ret_nhits - RETURN: number of hits
 *           ret_hitr  - RETURN: start states of hits, 0..nhits-1
 *           ret_hiti  - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj  - RETURN: end positions of hits, 0..nhits-1
 *           ret_hitsc - RETURN: scores of hits, 0..nhits-1            
 *           pad       - number of nucleotides to add on to start and end points
 * Returns:  
 *           hiti, hitj, hitsc are allocated here; caller free's w/ free().
 */
void
CP9_combine_FBscan_hits(int i0, int j0, int W, int fwd_nhits, int *fwd_hitr, int *fwd_hiti, 
			int *fwd_hitj, float *fwd_hitsc, int bck_nhits, 
			int *bck_hitr, int *bck_hiti, int *bck_hitj, float *bck_hitsc, 
			int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, 
			float **ret_hitsc, int pad)
{
  int i,j;
  int       nhits;		/* # of hits in optimal parse */
  int      *hitr;		/* initial state indices of hits in optimal parse */
  int      *hiti;		/* initial state indices of hits in optimal parse */
  int      *hitj;               /* end positions of hits in optimal parse */
  float    *hitsc;              /* scores of hits in optimal parse */
  int       alloc_nhits;	/* used to grow the hit arrays */
  int       fwd_ctr;
  int       bck_ctr;
  
  alloc_nhits = 10;
  hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitsc = MallocOrDie(sizeof(float) * alloc_nhits);

  /*if(fwd_nhits != bck_nhits)
    {
      printf("ERROR: fwd_nhits: %d | bck_nhits: %d\n", fwd_nhits, bck_nhits);
    }
  */

  /*
  for(i = 0; i < fwd_nhits; i++)
    {
      printf("FWD hit %3d | i: %3d | j: %3d | sc: %f\n", i, fwd_hiti[i], fwd_hitj[i], fwd_hitsc[i]);
    }
  printf("\n\n");
  for(i = 0; i < bck_nhits; i++)
    {
      printf("BCK hit %3d | i: %3d | j: %3d | sc: %f\n", i, bck_hiti[i], bck_hitj[i], bck_hitsc[i]);
    }
  printf("\n\n");
  */

  bck_ctr = 0;
  fwd_ctr = fwd_nhits - 1;
  nhits = 0;
  while(bck_ctr < bck_nhits || fwd_ctr >= 0)
    {
      /*printf("fwd_ctr: %d | bck_ctr: %d\n", fwd_ctr, bck_ctr);*/
      if(bck_ctr > (bck_nhits+1))
	{
	  printf("ERROR: bck_ctr: %d | bck_nhits: %d\n", bck_ctr, bck_nhits);
	  exit(1);
	}
      if(fwd_ctr < -1)
	{
	  printf("ERROR: fwd_ctr: %d\n", fwd_ctr);
	  exit(1);
	}
      if(fwd_ctr == -1 && bck_ctr == bck_nhits)
	{
	  printf("ERROR: bck_ctr: %d | fwd_ctr: %d\n", bck_ctr, fwd_ctr);
	  exit(1);
	}
      if(fwd_ctr == -1) 
	{
	  i = (bck_hiti[bck_ctr] - pad >= i0) ? (bck_hiti[bck_ctr] - pad) : i0;
	  j = (bck_hitj[bck_ctr] + pad <= j0) ? (bck_hitj[bck_ctr] + pad) : j0;
	  while ((j-i+1) > W)
	    { i++; if((j-i+1) > W) j--; }
	  hiti[nhits]   = i;
	  hitj[nhits]   = j;
	  hitsc[nhits]  = bck_hitsc[bck_ctr];
	  hitr[nhits]   = bck_hitr[bck_ctr];
	  bck_ctr++;
	}
      else if(bck_ctr == bck_nhits) 
	{
	  i = (fwd_hiti[fwd_ctr] - pad >= i0) ? (fwd_hiti[fwd_ctr] - pad) : i0;
	  j = (fwd_hitj[fwd_ctr] + pad <= j0) ? (fwd_hitj[fwd_ctr] + pad) : j0;
	  while ((j-i+1) > W)
	    { i++; if((j-i+1) > W) j--; }
	  hiti[nhits]   = i;
	  hitj[nhits]   = j;
	  hitsc[nhits]  = fwd_hitsc[fwd_ctr];
	  hitr[nhits]   = fwd_hitr[fwd_ctr];
	  fwd_ctr--;
	}
      else
	{
	  i = (bck_hiti[bck_ctr] - pad >= i0) ? (bck_hiti[bck_ctr] - pad) : i0;
	  j = (fwd_hitj[fwd_ctr] + pad <= j0) ? (fwd_hitj[fwd_ctr] + pad) : j0;
	  
	  if((j-i+1 > 0) && (j-i+1 <= (W + 2 * pad)))
	    {
	      /* creep in on i and j until we're at the window size */
	      while ((j-i+1) > W)
		{ i++; if((j-i+1) > W) j--; }
	      hiti[nhits]   = i;
	      hitj[nhits]   = j;
	      hitsc[nhits]  = 0.5 * (fwd_hitsc[fwd_ctr] + bck_hitsc[bck_ctr]);
	      if(fwd_hitr[fwd_ctr] != bck_hitr[bck_ctr]) { printf("ERROR: hitr's don't match!\n"); exit(1); }
	      hitr[nhits] = fwd_hitr[fwd_ctr];
	      fwd_ctr--;
	      bck_ctr++;
	    }
	  else if(j > i) /* hit from Backward() unmatched by a Forward hit */
	    {
	      hiti[nhits]   = i;
	      hitj[nhits]   = ((i+W-1) < j0) ? (i+W-1) : j0;
	      /*hitj[nhits]   = bck_hitj[bck_ctr]; *//* this will be (i+W) */
	      hitsc[nhits]  = bck_hitsc[bck_ctr];
	      hitr[nhits]   = bck_hitr[bck_ctr];
	      bck_ctr++;
	    }
	  else if(i > j) /* hit from Backward() unmatched by a Forward hit */
	    {
	      hitj[nhits]   = j;
	      hiti[nhits]   = ((j-W+1) > i0) ? (j-W+1) : i0;
	      /*hiti[nhits]   = fwd_hiti[fwd_ctr]; *//* this will be (j-W) */
	      hitsc[nhits]  = fwd_hitsc[fwd_ctr];
	      hitr[nhits]   = fwd_hitr[fwd_ctr];
	      fwd_ctr--;
	    }
	  else
	    {
	      printf("Uh... this shouldn't happen.\n");
	      exit(1);
	    }
	}

      nhits++;
      if (nhits == alloc_nhits) 
	{
	  hitr  = ReallocOrDie(hitr,  sizeof(int)   * (alloc_nhits + 10));
	  hitj  = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits + 10));
	  hiti  = ReallocOrDie(hiti,  sizeof(int)   * (alloc_nhits + 10));
	  hitsc = ReallocOrDie(hitsc, sizeof(float) * (alloc_nhits + 10));
	  alloc_nhits += 10;
	}
    }
  
  /*printf("\n\n");
  for(i = 0; i < nhits; i++)
    {
      printf("COMBINED hit %3d | i: %3d | j: %3d | sc: %f\n", i, hiti[i], hitj[i], hitsc[i]);
    }
  */

  *ret_nhits = nhits;
  *ret_hitr  = hitr;
  *ret_hiti  = hiti;
  *ret_hitj  = hitj;
  *ret_hitsc = hitsc;
}



/***********************************************************************
 * Function: CP9BackwardScan()
 * 
 * Purpose:  Runs the Backward dynamic programming algorithm that scans an 
 *           input subsequence (i0-j0). Complements CP9ForwardScan().
 *           Very similar to CP9Backward with the only difference being
 *           that parses can end at any position j, not only j0 (the endpoint).
 *           
 *           This function is messy in that it does not check that
 *           the log-odds scores are -INFTY (-987654321) before
 *           adding them. This means that -INFTY is interpreted
 *           as a possible transition, with a very low probability
 *           (2^-987654321).
 * 
 * Args:     dsq    - sequence in digitized form
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - max d: max size of a hit
 *           hmm    - the model
 *           ret_bmx - the Backward scanning matrix
 *           ret_nhits - RETURN: number of hits
 *           ret_hitr  - RETURN: start states of hits, 0..nhits-1
 *           ret_hiti  - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj  - RETURN: end positions of hits, 0..nhits-1
 *           ret_hitsc - RETURN: scores of hits, 0..nhits-1            
 *           min_thresh- minimum score to report (EPN via Alex Coventry 03.11.06)
 *
 * Returns:  
 *           hiti, hitj, hitsc are allocated in a helper function; caller free's w/ free().
 */
float
CP9BackwardScan(char *dsq, int i0, int j0, int W, CP9_t *hmm,
		       CP9_dpmatrix_t **ret_mx, int *ret_nhits, int **ret_hitr, int **ret_hiti, 
		       int **ret_hitj, float **ret_hitsc, float min_thresh)
{
  CP9_dpmatrix_t *mx;
  int **mmx;
  int **imx;
  int **dmx;
  int **emx;
  int   i,j,k;
  float     sc;     
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */

  int       nhits;		/* # of hits in optimal parse */
  int      *hitr;		/* initial state indices of hits in optimal parse */
  int      *hiti;		/* initial state indices of hits in optimal parse */
  int      *hitj;               /* end positions of hits in optimal parse */
  float    *hitsc;              /* scores of hits in optimal parse */
  int       alloc_nhits;	/* used to grow the hit arrays */

  int      L;                   /* j0-i0+1: subsequence length */
  int      jp;		        /* j': relative position in the subsequence  */
  int      ip;		        /* i': relative position in the subsequence  */
  /*int bck_sum;*/                  /* calc'ed from the Backward matrix, the sum over all
				 * positions i0..j0 of all parses that start and end 
				 * anywhere from i0..j0 */

  /*bck_sum = -INFTY;*/

  L = j0-i0+1;
  if (W > L) W = L; 

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. We use it first for
   * the Forward scan, and reuse it for the Backward scan.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+2));
  gamma[0] = 0.;
  gamma[L+1] = 0.;
  gback    = MallocOrDie(sizeof(int)   * (L+2));
  gback[0] = -1;
  gback[L+1] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+2));
  saver    = MallocOrDie(sizeof(int)   * (L+2));

  /* Allocate a DP matrix with 0..L rows, 0..M-1 
   */ 
  mx = AllocCPlan9Matrix(L+1, hmm->M, &mmx, &imx, &dmx, &emx);

  /*****************************************************************
   * The main loop: scan the sequence from position j0 to i0.
   *****************************************************************/

  /* Initialization of the L row.
   */
  emx[0][L] = 0; /*have to end in E*/

  mmx[L][hmm->M] = emx[0][L] + hmm->esc[hmm->M]; /* M<-E ...                   */
  mmx[L][hmm->M] += hmm->msc[(int) dsq[j0]][hmm->M]; /* ... + emitted match symbol */
  imx[L][hmm->M] = emx[0][L] + hmm->tsc[CTIM][hmm->M];   /* I_M(C)<-E ... */
  imx[L][hmm->M] += hmm->isc[(int) dsq[j0]][hmm->M];           /* ... + emitted match symbol */
  dmx[L][hmm->M] = emx[0][L] + hmm->tsc[CTDM][hmm->M];    /* D_M<-E */
  for (k = hmm->M-1; k >= 1; k--)
    {
      mmx[L][k]  = emx[0][L] + hmm->esc[k];
      mmx[L][k]  = ILogsum(mmx[L][k], dmx[L][k+1] + hmm->tsc[CTMD][k]);
      mmx[L][k] += hmm->msc[(int) dsq[j0]][k];

      imx[L][k] = dmx[L][k+1] + hmm->tsc[CTID][k];
      imx[L][k] += hmm->isc[(int) dsq[j0]][k];

      dmx[L][k] = dmx[L][k+1] + hmm->tsc[CTDD][k];
    }
  
  mmx[L][0] = -INFTY; /* no esc[0] b/c its impossible */
  imx[L][0] = dmx[L][1] + hmm->tsc[CTID][0];
  imx[L][0] += hmm->isc[(int) dsq[j0]][hmm->M];    
  dmx[L][0] = -INFTY; /*D_0 doesn't exist*/

  /*bck_sum = ILogsum(bck_sum, mmx[L][0]);*/

  /*printf("mmx[L][  0]: %10.5f\n", sc);*/
  /*printf("tmp_bck_sum L: %d\n", tmp_bck_sum);*/

  /* Recursion. Done as a pull.
   */
  for (i = j0-1; i >= i0-1; i--)
    {
      ip = i-i0 + 1;		/* ip is relative index in dsq (1 to L) */
      if(ip > 0)
	{
	  emx[0][ip] = 0; /* This is 1 of 3 differences between a Backward scanner and the 
			   * regular Backward (CP9Backward), 
			   * in CP9Backward, this cell is set to -INFTY.
			   * Here, we can end at the END state at any position,
			   * this isn't true in global alignment. 
			   * Second and third differences are related, 
			   * see the four comments 
			   * labelled "2nd diff" or "3rd diff" below
			   */
	  
	  /* Now the main states. Note the boundary conditions at M.
	   */
	  mmx[ip][hmm->M] = ILogsum(emx[0][ip+1] + hmm->esc[hmm->M], /* M<-E ... 2nd diff w/non-scanner*/
				    imx[ip+1][hmm->M] + hmm->tsc[CTMI][hmm->M]);
	  mmx[ip][hmm->M] += hmm->msc[(int) dsq[i]][hmm->M];
	  imx[ip][hmm->M] = ILogsum(emx[0][ip] + hmm->tsc[CTIM][hmm->M],    /* I_M(C)<-E ... */
				    imx[ip+1][hmm->M] + hmm->tsc[CTII][hmm->M]);
	  imx[ip][hmm->M] += hmm->isc[(int) dsq[i]][hmm->M];
	  dmx[ip][hmm->M] = ILogsum(emx[0][ip] + hmm->tsc[CTDM][hmm->M], /* D_M<-E */
				    imx[ip+1][hmm->M] + hmm->tsc[CTDI][hmm->M]);  
	  for (k = hmm->M-1; k >= 1; k--)
	    {
	      mmx[ip][k]  = ILogsum(ILogsum(emx[0][ip+1] + hmm->esc[k], /* M<-E ... 2nd diff w/non-scanner*/
					    mmx[ip+1][k+1] + hmm->tsc[CTMM][k]),
				    ILogsum(imx[ip+1][k] + hmm->tsc[CTMI][k],
					    dmx[ip][k+1] + hmm->tsc[CTMD][k]));	  
	      
	      mmx[ip][k] += hmm->msc[(int) dsq[i]][k];
	      
	      imx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTIM][k],
					    imx[ip+1][k] + hmm->tsc[CTII][k]),
				    dmx[ip][k+1] + hmm->tsc[CTID][k]);
	      imx[ip][k] += hmm->isc[(int) dsq[i]][k];
	      
	      dmx[ip][k]  = ILogsum(ILogsum(mmx[ip+1][k+1] + hmm->tsc[CTDM][k],
					    imx[ip+1][k] + hmm->tsc[CTDI][k]),
				    dmx[ip][k+1] + hmm->tsc[CTDD][k]);
	    }
	  imx[ip][0]  = ILogsum(ILogsum(mmx[ip+1][1] + hmm->tsc[CTIM][0],
					imx[ip+1][0] + hmm->tsc[CTII][0]),
				dmx[ip][1] + hmm->tsc[CTID][0]);
	  imx[ip][0] += hmm->isc[(int) dsq[i]][0];
	  mmx[ip][0] = -INFTY;
	  /*for (k = hmm->M-1; k >= 1; k--)*/ /*M_0 is the B state, it doesn't emit*/
	  for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	    {
	      mmx[ip][0] = ILogsum(mmx[ip][0], mmx[ip+1][k] + hmm->bsc[k]);
	      /*printf("hmm->bsc[%d]: %d | mmx[ip+1][k]: %d\n", k, hmm->bsc[k], mmx[(i+1)][k]);
		printf("mmx[%3d][%3d]: %d | k: %3d\n", i, 0, mmx[ip][0], k);*/
	    }
	  mmx[ip][0] = ILogsum(mmx[ip][0], imx[ip+1][k] + hmm->tsc[CTMI][k]);
	  mmx[ip][0] = ILogsum(mmx[ip][0], dmx[ip][k+1] + hmm->tsc[CTMD][k]);
	  
	  dmx[ip][0] = -INFTY; /* D_0 does not exist */
	  
	}
      if(ip == 0)
	{
	  /* case when ip = 0 */
	  emx[0][0] = 0.; /* 3rd diff with non-scanner */
	  mmx[0][hmm->M] = -INFTY; /* need seq to get here */
	  imx[0][hmm->M] = -INFTY; /* need seq to get here */
	  dmx[0][hmm->M] = ILogsum(emx[0][0] + hmm->tsc[CTDM][hmm->M],  /*D_M<-E*/ /* 3rd diff with non-scanner */
				   imx[1][hmm->M] + hmm->tsc[CTDI][hmm->M]);  
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
	      /*printf("k: %d | mmx[0][0]: %d | mmx[1][k]: %d | hmm->bsc[k]: %d\n", k, mmx[0][0], mmx[1][k], hmm->bsc[k]);*/
	    }
	  mmx[0][0] = ILogsum(mmx[0][0], dmx[0][1] + hmm->tsc[CTMD][0]);
	  /*printf("added D trans | mmx[0][0]: %d \n", mmx[0][0]);*/
	  mmx[0][0] = ILogsum(mmx[0][0], imx[1][0] + hmm->tsc[CTMI][0]);
	  /*printf("added I trans | mmx[0][0]: %d \n", mmx[0][0]);*/
	  
	  sc = Scorify(mmx[0][0]);   /* the total Backward score. */
	  /*printf("mmx[0][  0]: %10.5f\n", sc);*/
	  
	  /*bck_sum = ILogsum(bck_sum, mmx[0][0]);*/
	  /*printf("Total Backward score (summed over all posns): %f\n", Scorify(bck_sum));		*/
	  /*sc  = Scorify(bck_sum);*/
	}
      /* The little semi-Markov model that deals with multihit parsing:
       */
      /* Here's a vicious off-by-one with HMMs and CMs:
       * an HMM parse can begin at position 0, i.e. backward[0][0] is the probability
       * of the entire sequence being emitted starting at the Begin state at position 0.
       * But a CM always has parses that start at position 1. So we have to translate
       * from HMM to CM in this way, and that's why there's a i+1 in the following loop.
       */
      gamma[ip+1]  = gamma[ip+2] + 0; /* extend without adding a new hit */
      /*printf("i: %d | gamma[i]: %f | gamma[i+1]: %f\n", i, gamma[i], gamma[i+1]);*/
      gback[ip+1]  = -1;
      savesc[ip+1] = IMPOSSIBLE;
      saver[ip+1]  = -1;
      j = ((i+1+W-1) < j0) ? (i+1+W-1) : j0;
      jp = j-i0+1;
      sc = gamma[jp+1] + Scorify(mmx[ip][0]) - min_thresh;
      
      /*printf("mmx[%3d][0]: %10.5f\n", i, Scorify(mmx[ip][0]));*/
      /*bck_sum = ILogsum(bck_sum, mmx[ip][0]);*/
      
      sc = Scorify(mmx[ip][0]);
      /*printf("mmx[%3d][0]: %10.5f\n", i, sc);*/
      /*printf("i: %d | mmx[ip][0]: %f | j: %d | gamma[i]: %f | sc: %f\n", i, Scorify(mmx[ip][0]), j, gamma[i], sc);*/
      if (sc > gamma[ip+1])
	{
	  gamma[ip+1]  = sc;
	  gback[ip+1]  = j;
	  savesc[ip+1] = Scorify(mmx[ip][0]);
	  saver[ip+1]  = 0;
	  //saver[ip+1]  = bestr[d];
	}
    }
  /* End of Backward() code */

  /*****************************************************************
   * Traceback stage for Backward.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 

  alloc_nhits = 10;
  hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitsc = MallocOrDie(sizeof(float) * alloc_nhits);

  i     = i0;
  nhits = 0;
  ip = i-i0+1;
  while (i <= L) {
    if (gback[ip] == -1) /* no hit */
      { i++; ip++; }
    else                /* a hit, a palpable hit */
      {
	hitr[nhits]   = saver[ip];
	hitj[nhits]   = gback[ip];
	hiti[nhits]   = i;
	hitsc[nhits]  = savesc[ip];
	nhits++;
	i = gback[ip]+1;
	ip = i-i0+1;
	if (nhits == alloc_nhits) {
	  hitr  = ReallocOrDie(hitr,  sizeof(int)   * (alloc_nhits + 10));
	  hitj  = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits + 10));
	  hiti  = ReallocOrDie(hiti,  sizeof(int)   * (alloc_nhits + 10));
	  hitsc = ReallocOrDie(hitsc, sizeof(float) * (alloc_nhits + 10));
	  alloc_nhits += 10;
	}
      }
  }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  *ret_nhits = nhits;
  *ret_hitr  = hitr;
  *ret_hiti  = hiti;
  *ret_hitj  = hitj;
  *ret_hitsc = hitsc;

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);
  return sc;		/* the total Backward score. */
}

#if 0 
/* Below are CP9 functions for scanning a sequence that are memory inefficient, 
 * in that they keep the entire matrix in memory. I originally wrote these
 * functions b/c so I wouldn't need to 'rescan' the hits that survive the 
 * CP9 filter to get HMM bands. The currently implemented versions of
 * these functions are memory efficient in that they only allocate
 * and fill 2 rows of a DP matrix. I think these are more useful right now,
 * but the old versions are kept here for reference 
 * (EPN, Mon Jan  8 15:57:13 2007)
 */
/***********************************************************************
 * Function: CP9ForwardScan()
 * 
 * Purpose:  A Forward dynamic programming algorithm that scans an input sequence.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 * 
 * Args:     
 *           cm      - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq     - sequence in digitized form
 *           i0      - start of target subsequence (1 for beginning of dsq)
 *           j0      - end of target subsequence (L for end of dsq)
 *           W       - the maximum size of a hit (often cm->W)
 *           cutoff  - minimum score to report
 *           ret_mx  - RETURN: the Forward matrix alloc'ed here (NULL if not wanted)
 *           results - scan_results_t to add to; if NULL, don't keep results
 *
 * Returns:  best_score, score of maximally scoring end position j 
 */
float
CP9ForwardScan(CM_t *cm, char *dsq, int i0, int j0, int W, 
	       float cutoff, CP9_dpmatrix_t **ret_mx, 
	       scan_results_t *results)
{
  CP9_dpmatrix_t *mx;
  int **mmx;
  int **imx;
  int **dmx;
  int **emx;
  int   i,j,k;
  float sc;
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */

  int      L;                   /* j0-i0+1: subsequence length */
  int      jp;		        /* j': relative position in the subsequence  */
  int      ip;		        /* i': relative position in the subsequence  */

  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0. */

  if(cm->cp9 == NULL)
    Die("ERROR in CP9ForwardScan, but cm->cp9 is NULL.\n");

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+1));
  gamma[0] = 0.;
  gback    = MallocOrDie(sizeof(int)   * (L+1));
  gback[0] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+1));
  saver    = MallocOrDie(sizeof(int)   * (L+1));

  /* Allocate a DP matrix with 0..L rows, 0..M-1 
   */ 
  mx = AllocCPlan9Matrix(L+1, cm->cp9->M, &mmx, &imx, &dmx, &emx);

  /* Initialization of the zero row.
   */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= cm->cp9->M; k++)
    mmx[0][k] = imx[0][k] = -INFTY;      /* need seq to get here */
  for (k = 1; k <= cm->cp9->M; k++)
    dmx[0][k] = ILogsum(ILogsum(mmx[0][k-1] + cm->cp9->tsc[CTMD][k-1],
				imx[0][k-1] + cm->cp9->tsc[CTID][k-1]),
			dmx[0][k-1] + cm->cp9->tsc[CTDD][k-1]);
  
  emx[0][0] = dmx[0][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]; 

  /*****************************************************************
   * The main loop: scan the sequence from position i0 to j0.
   *****************************************************************/
  
  /* Recursion. Done as a pull.
   */
  for (j = i0; j <= j0; j++)
    {
      jp = j-i0+1;     /* jp is relative position in the sequence */
      mmx[jp][0] = 0;  /* This is the 1 difference between a Forward scanner and the 
			* regular Forward (CP9Forward), 
			* in CP9Forward, this cell is set to -INFTY.
			* Here, we can start at any position in the seq,
			* this isn't true in global alignment. 
			*/
      dmx[jp][0] = -INFTY;  /*D_0 is non-existent*/
      imx[jp][0]  = ILogsum(ILogsum(mmx[jp-1][0] + cm->cp9->tsc[CTMI][0],
				    imx[jp-1][0] + cm->cp9->tsc[CTII][0]),
			    dmx[jp-1][0] + cm->cp9->tsc[CTDI][0]);
      imx[jp][0] += cm->cp9->isc[(int) dsq[j]][0];
      
      for (k = 1; k <= cm->cp9->M; k++)
	{
	  mmx[jp][k]  = ILogsum(ILogsum(mmx[jp-1][k-1] + cm->cp9->tsc[CTMM][k-1],
				       imx[jp-1][k-1] + cm->cp9->tsc[CTIM][k-1]),
			       ILogsum(mmx[jp-1][0] + cm->cp9->bsc[k],
				       dmx[jp-1][k-1] + cm->cp9->tsc[CTDM][k-1]));
	  mmx[jp][k] += cm->cp9->msc[(int) dsq[j]][k];
	  
	  dmx[jp][k]  = ILogsum(ILogsum(mmx[jp][k-1] + cm->cp9->tsc[CTMD][k-1],
					imx[jp][k-1] + cm->cp9->tsc[CTID][k-1]),
				dmx[jp][k-1] + cm->cp9->tsc[CTDD][k-1]);
	  
	  imx[jp][k]  = ILogsum(ILogsum(mmx[jp-1][k] + cm->cp9->tsc[CTMI][k],
				       imx[jp-1][k] + cm->cp9->tsc[CTII][k]),
			       dmx[jp-1][k] + cm->cp9->tsc[CTDI][k]);
	  imx[jp][k] += cm->cp9->isc[(int) dsq[j]][k];
	  //printf("mmx[%d][%d]: %d\n", i, k, mmx[jp][k]);
	  //printf("imx[%d][%d]: %d\n", i, k, imx[jp][k]);
	  //printf("dmx[%d][%d]: %d\n", i, k, dmx[jp][k]);

	}
      
      emx[0][jp] = -INFTY;
      for (k = 1; k <= cm->cp9->M; k++)
	emx[0][jp] = ILogsum(emx[0][jp], mmx[jp][k] + cm->cp9->esc[k]);
      emx[0][jp] = ILogsum(emx[0][jp], dmx[jp][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]); 
      emx[0][jp] = ILogsum(emx[0][jp], imx[jp][cm->cp9->M] + cm->cp9->tsc[CTIM][cm->cp9->M]); 
      /* transition from D_M -> end */
      sc = Scorify(emx[0][jp]);

      if(sc > best_neg_score) 
	best_neg_score = sc;

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[jp]  = gamma[jp-1] + 0; /* extend without adding a new hit */
      gback[jp]  = -1;
      savesc[jp] = IMPOSSIBLE;
      saver[jp]  = -1;
      i = ((j-W+1)> i0) ? (j-W+1) : i0;
      ip = i-i0+1;
      sc = gamma[ip-1] + Scorify(emx[0][jp]) + cm->cp9_sc_boost;
      /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
       * value is 0.0 if technique not used. */
      if (sc > gamma[jp])
	{
	  gamma[jp]  = sc;
	  gback[jp]  = i;
	  savesc[jp] = Scorify(emx[0][jp]);
	  saver[jp]  = 0;
	  //saver[jp]  = bestr[d];
	}
    } /* end loop over end positions j */
  sc = Scorify(emx[0][L]);
  //printf("END %d %10.2f\n", j, sc);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j > i0) {
    jp = j-i0+1;
    if (gback[jp] == -1) /* no hit */
      j--; 
    else                /* a hit, a palpable hit */
      {
	if(savesc[jp] > best_score) 
	  best_score = savesc[jp];
	if(savesc[jp] >= cutoff && results != NULL) /* report the hit */
	  report_hit(gback[jp], j, saver[jp], savesc[jp], results);
	j = gback[jp]-1;
      }
  }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);
  
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeCPlan9Matrix(mx);
  
  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;
  return best_score;
}

/* This function is never called, but could be useful for debugging
 * in the future */

static void debug_check_CP9_FBscan(CP9_dpmatrix_t *fmx, CP9_dpmatrix_t *bmx, 
				   CP9_t *hmm, float sc, int L, char *dsq);

/*********************************************************************
 * Function: debug_CP9_check_FBscan()
 * 
 * Purpose:  Debugging function for CP9ForwardBackwardScan(),
 *           print probability each position of the target
 *           sequence is emitted by some state of the model.
 *           Notably different from debug_CP9_check_FB() 
 *           because parses are allowed to start at any
 *           position and end at any position (during alignment
 *           its assumed the entire sequence is involved in 
 *           each parse).
 *           
 * Args:     fmx    - forward dp matrix, already filled
 *           bmx    - backward dp matrix, already filled
 *           hmm    - the model
 *           sc     - ~P(x|hmm) the summed probability of 
 *                    subsequences of the target sequence
 *                    given the model
 *           L      - length of the sequence
 *           
 * Return:   (void) Exits if any errors are found.
 */
static void
debug_check_CP9_FBscan(CP9_dpmatrix_t *fmx, CP9_dpmatrix_t *bmx, 
		       CP9_t *hmm, float sc, int L, char *dsq)
{
  int k, i;
  float max_diff;  /* maximum allowed difference between sc and 
		    * sum_k f[i][k] * b[i][k] for any i */
  int fb_sum;
  float fb_sc;

  max_diff = 0.01;
  /*printf("sc: %f\n", sc);*/

  /* In all possible paths through the model, each residue of the sequence must have 
   * been emitted by exactly 1 insert or match state. */
  for (i = 0; i <= L; i++)
    {
      fb_sum = -INFTY;
      for (k = 0; k <= hmm->M; k++) 
	{
	  if(i == 0)
	    {
	      printf("fmx->mmx[%d][%d]: %d\n", i, k, fmx->mmx[i][k]);
	      printf("bmx->mmx[%d][%d]: %d\n", i, k, bmx->mmx[i][k]);
	      printf("hmm->msc[dsq[i]][%d]: %d\n", k, hmm->msc[(int) dsq[i]][k]);
	    }
	  fb_sum = ILogsum(fb_sum, (fmx->mmx[i][k] + bmx->mmx[i][k] - hmm->msc[(int) dsq[i]][k]));
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  fb_sum = ILogsum(fb_sum, (fmx->imx[i][k] + bmx->imx[i][k] - hmm->isc[(int) dsq[i]][k]));
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->imx and bmx->imx*/
	  /*fb_sum = ILogsum(fb_sum, fmx->dmx[i][k] + bmx->dmx[i][k]);*/
	}
      fb_sc  = Scorify(fb_sum);
      printf("position %5d | sc: %f\n", i, fb_sc);
      /*      diff = sc - fb_sc;
      if(diff < 0.) diff *= -1.;
      if(diff > max_diff)
	{
	  printf("ERROR, fb_sc[%d]: %f too different from P(x|hmm): %f\n", i, fb_sc, sc);
	  exit(1);
	}
      */
    }
}
#endif

/* CP9_scan.c 
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
 * CP9BackwardScan()     - Scan input sequence for high scoring
 *                         backward hits to the model.
 * CP9FilteredScan()     - use CP9ForwardScan() and potentially CP9BackwardScan()
 *                         (if CM_SEARCH_HMMFB) to filter DB, 
 *                         and pass promising subseqs to CM search,
 *                         can be run in different modes.
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
#include "cm_dispatch.h"
  
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
 *           ret_maxres- RETURN: start position that gives maximum score max argmax_i isc[i]
 *           results   - scan_results_t to add to; if NULL, don't keep results
 *           doing_rescan - TRUE if we've called this function recursively, and we don't
 *                          want to again
 *
 * Returns:  best_sc, score of maximally scoring end position j 
 */
float
CP9ForwardScan(CM_t *cm, char *dsq, int i0, int j0, int W, float cutoff, int **ret_isc, 
	       int **ret_hitj, int *ret_nhits, int *ret_bestpos, scan_results_t *results, 
	       int doing_rescan)
{
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /* j-W: position in the subsequence                             */
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
  float        best_pos;    /* residue giving best score overall                            */
  float        best_negsc;  /* Best score overall score to return, used if all scores < 0.  */
  float        best_negpos; /* residue giving best score overall score used if all scores < 0.*/
  int          accept;      /* Flag used if cm->search_opts & CM_SEARCH_HMMRESCAN           */
  float        temp_sc;     /* temporary score                                              */

  /* Contract checks */
  if(cm->cp9 == NULL)
    Die("ERROR in CP9ForwardScan, cm->cp9 is NULL.\n");

  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_negsc  = IMPOSSIBLE;
  best_negpos = -1;
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
  /*printf("mmx[jp:%d][%d]: %d\n", 0, 0, mmx[0][0]);
    printf("imx[jp:%d][%d]: %d\n", 0, 0, imx[0][0]);
    printf("dmx[jp:%d][%d]: %d\n", 0, 0, dmx[0][0]);*/

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= cm->cp9->M; k++)
    {
      mmx[0][k] = imx[0][k] = -INFTY;      /* need seq to get here */
      dmx[0][k] = ILogsum(ILogsum(mmx[0][k-1] + cm->cp9->tsc[CTMD][k-1],
				  imx[0][k-1] + cm->cp9->tsc[CTID][k-1]),
			  dmx[0][k-1] + cm->cp9->tsc[CTDD][k-1]);
      /*printf("mmx[jp:%d][%d]: %d\n", 0, k, mmx[0][k]);
	printf("imx[jp:%d][%d]: %d\n", 0, k, imx[0][k]);
	printf("dmx[jp:%d][%d]: %d\n", 0, k, dmx[0][k]);*/
    }
  /* We can do a full parse through all delete states. */
  isc[0] = dmx[0][cm->cp9->M] + cm->cp9->tsc[CTDM][cm->cp9->M]; 
  /*printf("jp: %d fsc: %d\n\n", 0, isc[0]);*/

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
      /*printf("mmx[jp:%d][%d]: %d\n", jp, 0, mmx[cur][0]);
	printf("imx[jp:%d][%d]: %d\n", jp, 0, imx[cur][0]);
	printf("dmx[jp:%d][%d]: %d\n", jp, 0, dmx[cur][0]);*/
      
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
      /*printf("jp: %d fsc: %f\n\n", jp, fsc);*/

      if(fsc > best_negsc) 
	{
	  best_negsc = fsc;
	  best_negpos= j;
	}
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
	  {
	    best_sc = savesc[jp];
	    best_pos= j;
	  }
	if(savesc[jp] >= cutoff)
	{
	  accept = TRUE;
	  /* Potentially rescan just the subseq that is the hit we're about to report.
	   * Implemented to deal with fact that --enfseq option was enforcing the subseq
	   * to have hit pass filter b/c this Forward scanning function is 'infinite' length
	   * (Weinberg). Sometimes the subseq we're about to report has a really
	   * crappy score, even though the cumulative Forward score (starting at i0) is good. */
	  if(cm->search_opts & CM_SEARCH_HMMRESCAN && doing_rescan == FALSE)
	    {
	      /*printf("rechecking hit from %d to %d\n", gback[jp], j);*/
	      temp_sc = CP9ForwardScan(cm, dsq, gback[jp], j, cm->W, cutoff, 
				       NULL,  /* don't care about scores of each pos */
				       NULL,  /* don't care about locations of hits */
				       NULL,  /* don't care about num hits found */
				       NULL,  /* don't care about best scoring position */
				       NULL,  /* don't report hits */
				       TRUE); /* set the doing_rescan arg to TRUE, 
						 so we don't potentially infinitely recurse */
	      /*printf("new score: %f old score %f\n", temp_sc, savesc[jp]);*/
	      if(temp_sc >= cutoff) 
		{ 
		  accept = TRUE; 
		  /*printf("rechecked hit from %d to %d\n", gback[jp], j);
		    printf("new score: %f old score %f\n", temp_sc, savesc[jp]);*/
		}
	      else accept = FALSE;
	      savesc[jp] = temp_sc;
	    }
	  if(accept)
	    {
	      if(results != NULL) /* report the hit */
		{
		  report_hit(gback[jp], j, 0, savesc[jp], results); 
		  /* 0 is for saver, which is irrelevant for HMM hits */
		}
	      hitj[nhits] = j;
	      nhits++;
	      if (nhits == alloc_nhits) 
		{
		  alloc_nhits += 10;
		  hitj = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits));
		}
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
    {
      best_sc = best_negsc;
      if(ret_bestpos != NULL) *ret_bestpos = best_negpos;
    }
  else if(ret_bestpos != NULL) *ret_bestpos = best_pos;
  return best_sc;
}

/***********************************************************************
 * Function: CP9BackwardScan()
 * 
 * Purpose:  Runs the Backward dynamic programming algorithm that scans an 
 *           input subsequence (i0-j0). Complements CP9ForwardScan().
 *           Only allocates 2 rows of the Backward matrix.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           Very similar to CP9Backward with the only difference being
 *           that parses can end at any position j, not only j0 (the endpoint).
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
 *           ret_hiti  - RETURN: start positions of hits, 0..ret_nhits-1
 *           ret_bestpos- RETURN: start position that gives maximum score max argmax_i isc[i]
 *           results   - scan_results_t to add to; if NULL, don't keep results
 *           doing_rescan - TRUE if we've called this function recursively, and we don't
 *                          want to again
 *
 * Returns:  best_sc, score of maximally scoring start position i
 */
float
CP9BackwardScan(CM_t *cm, char *dsq, int i0, int j0, int W, float cutoff, int **ret_isc, 
		int **ret_hiti, int *ret_nhits, int *ret_bestpos, scan_results_t *results, 
		int doing_rescan)
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
  int         *hiti;        /* start positions (i) of hits [0..nhits-1]                     */
  int          alloc_nhits; /* used to grow the hitj array                                  */
  float        best_sc;     /* Best overall score from semi-HMM to return                   */
  float        best_pos;    /* residue giving best score overall                            */
  float        best_negsc;  /* Best score overall score to return, used if all scores < 0.  */
  float        best_negpos; /* residue giving best score overall score used if all scores < 0.*/
  int          accept;      /* Flag used if cm->search_opts & CM_SEARCH_HMMRESCAN           */
  float        temp_sc;     /* temporary score                                              */

  /* Contract checks */
  if(cm->cp9 == NULL)
    Die("ERROR in CP9BackwardScan, cm->cp9 is NULL.\n");

  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_negsc  = IMPOSSIBLE;
  best_negpos = -1;
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
  /* isc will hold P(seq from i..j0 | Model) for each i in int log odds form */
  isc   = MallocOrDie(sizeof(int) * (j0-i0+2));

  /* Initialization of the L (i = j0, cur = i-j0%2 = 0) row. */
  cur = 0;
  mmx[cur][cm->cp9->M]  = cm->cp9->esc[cm->cp9->M];                 /* M<-E ...                   */
  mmx[cur][cm->cp9->M] += cm->cp9->msc[(int) dsq[j0]][cm->cp9->M];  /* ... + emitted match symbol */
  imx[cur][cm->cp9->M]  = cm->cp9->tsc[CTIM][cm->cp9->M];           /* I_M(C)<-E ...              */
  imx[cur][cm->cp9->M] += cm->cp9->isc[(int) dsq[j0]][cm->cp9->M];  /* ... + emitted match symbol */
  dmx[cur][cm->cp9->M]  = cm->cp9->tsc[CTDM][cm->cp9->M];           /* D_M<-E                     */
  /*printf("mmx[jp:%d][%d]: %d\n", 0, cm->cp9->M, mmx[cur][cm->cp9->M]);
    printf("imx[jp:%d][%d]: %d\n", 0, cm->cp9->M, imx[cur][cm->cp9->M]);
    printf("dmx[jp:%d][%d]: %d\n", 0, cm->cp9->M, dmx[cur][cm->cp9->M]);*/
  for (k = cm->cp9->M-1; k >= 1; k--)
    {
      mmx[cur][k]  = cm->cp9->esc[k];
      mmx[cur][k]  = ILogsum(mmx[cur][k], dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]);
      mmx[cur][k] += cm->cp9->msc[(int) dsq[j0]][k];

      imx[cur][k] = dmx[cur][k+1] + cm->cp9->tsc[CTID][k];
      imx[cur][k] += cm->cp9->isc[(int) dsq[j0]][k];

      dmx[cur][k] = dmx[cur][k+1] + cm->cp9->tsc[CTDD][k];
      /*printf("mmx[ip:%d][%d]: %d\n", 0, k, mmx[cur][k]);
	printf("imx[ip:%d][%d]: %d\n", 0, k, imx[cur][k]);
	printf("dmx[ip:%d][%d]: %d\n", 0, k, dmx[cur][k]);*/
    }
  
  mmx[cur][0]  = -INFTY; /* no esc[0] b/c its impossible */
  imx[cur][0]  = dmx[cur][1] + cm->cp9->tsc[CTID][0];
  imx[cur][0] += cm->cp9->isc[(int) dsq[j0]][cm->cp9->M];    
  dmx[cur][0]  = -INFTY; /*D_0 doesn't exist*/

  /* We can do a full parse through all delete states, prob of getting
   * to D_1 time prob of transiting to D_1 from M_0 (which is the B state) */
  isc[L] = dmx[cur][1] + cm->cp9->tsc[CTMD][0]; 
  /*printf("ip: %d fsc: %d\n\n", L, isc[L]);*/

  /*****************************************************************
   * The main loop: scan the sequence from position j0 to i0.
   *****************************************************************/
  /* Recursion */
  for (i = j0-1; i >= (i0-1); i--)
    {
      ip = i-(i0-1);		/* ip is relative index in dsq (0 to L) */
      cur = (j0-i)  % 2;
      prv = (j0-i+1)% 2;
      /* A main difference between a Backward scanner and the 
       * regular Backward (CP9Backward()) is that a scanner
       * can end at the END state at any position,
       * this isn't true in global alignment. Relevant code
       * that is unique to a scanner relative to global is marked with
       * 'unique to scanner' below.
       */
      if(ip > 0)
	{
	  /* Now the main states. Note the boundary conditions at M. */
	  mmx[cur][cm->cp9->M] =  ILogsum(cm->cp9->esc[cm->cp9->M], /* M_M<-E ... unique to scanner */ 
				      imx[prv][cm->cp9->M] + cm->cp9->tsc[CTMI][cm->cp9->M]);
	  mmx[cur][cm->cp9->M] += cm->cp9->msc[(int) dsq[i]][cm->cp9->M];
	  imx[cur][cm->cp9->M] =  ILogsum(cm->cp9->tsc[CTIM][cm->cp9->M],    /* I_M(C)<-E ... unique to scanner */
				      imx[prv][cm->cp9->M] + cm->cp9->tsc[CTII][cm->cp9->M]);
	  imx[cur][cm->cp9->M] += cm->cp9->isc[(int) dsq[i]][cm->cp9->M];
	  dmx[cur][cm->cp9->M] =  ILogsum(cm->cp9->tsc[CTDM][cm->cp9->M], /* D_M<-E unique to scanner */
				      imx[prv][cm->cp9->M] + cm->cp9->tsc[CTDI][cm->cp9->M]);  
	  /*printf("mmx[ip:%d][%d]: %d\n", ip, cm->cp9->M, mmx[cur][cm->cp9->M]);
	    printf("imx[ip:%d][%d]: %d\n", ip, cm->cp9->M, imx[cur][cm->cp9->M]);
	    printf("dmx[ip:%d][%d]: %d\n", ip, cm->cp9->M, dmx[cur][cm->cp9->M]);*/
      
	  for (k = cm->cp9->M-1; k >= 1; k--)
	    {
	      mmx[cur][k]  = ILogsum(ILogsum(cm->cp9->esc[k], /* M<-E ... unique to scanner */
					     mmx[prv][k+1] + cm->cp9->tsc[CTMM][k]),
				     ILogsum(imx[prv][k] + cm->cp9->tsc[CTMI][k],
					     dmx[cur][k+1] + cm->cp9->tsc[CTMD][k]));	  
	  
	      mmx[cur][k] += cm->cp9->msc[(int) dsq[i]][k];
	  
	      imx[cur][k]  = ILogsum(ILogsum(mmx[prv][k+1] + cm->cp9->tsc[CTIM][k],
					     imx[prv][k] + cm->cp9->tsc[CTII][k]),
				     dmx[cur][k+1] + cm->cp9->tsc[CTID][k]);
	      imx[cur][k] += cm->cp9->isc[(int) dsq[i]][k];
	  
	      dmx[cur][k]  = ILogsum(ILogsum(mmx[prv][k+1] + cm->cp9->tsc[CTDM][k],
					     imx[prv][k] + cm->cp9->tsc[CTDI][k]),
				     dmx[cur][k+1] + cm->cp9->tsc[CTDD][k]);
	      /*printf("mmx[ip:%d][%d]: %d\n", ip, k, mmx[cur][k]);
		printf("imx[ip:%d][%d]: %d\n", ip, k, imx[cur][k]);
		printf("dmx[ip:%d][%d]: %d\n", ip, k, dmx[cur][k]);*/
	    }
	  imx[cur][0]  = ILogsum(ILogsum(mmx[prv][1] + cm->cp9->tsc[CTIM][0],
					 imx[prv][0] + cm->cp9->tsc[CTII][0]),
				 dmx[cur][1] + cm->cp9->tsc[CTID][0]);
	  imx[cur][0] += cm->cp9->isc[(int) dsq[i]][0];
	}

      if(ip == 0)
	{
	  /* special case, no sequence has been emitted */
	  mmx[cur][cm->cp9->M] = -INFTY; /* need seq to get here */
	  imx[cur][cm->cp9->M] = -INFTY; /* need seq to get here */
	  dmx[cur][cm->cp9->M] = ILogsum(cm->cp9->tsc[CTDM][cm->cp9->M],  /*D_M<-E* ... unique to scanner */
				     imx[prv][cm->cp9->M] + cm->cp9->tsc[CTDI][cm->cp9->M]);  
	  for (k = cm->cp9->M-1; k >= 1; k--)
	    {
	      mmx[cur][k] = -INFTY; /* need seq to get here */
	      imx[cur][k] = -INFTY; /* need seq to get here */
	      dmx[cur][k]  = ILogsum(ILogsum(mmx[prv][k+1] + cm->cp9->tsc[CTDM][k],
					     imx[prv][k] + cm->cp9->tsc[CTDI][k]),
				     dmx[cur][k+1] + cm->cp9->tsc[CTDD][k]);
	      /*printf("mmx[ip:%d][%d]: %d\n", ip, k, mmx[cur][k]);
		printf("imx[ip:%d][%d]: %d\n", ip, k, imx[cur][k]);
		printf("dmx[ip:%d][%d]: %d\n", ip, k, dmx[cur][k]);*/
	    }
	  imx[cur][0] = -INFTY; /*need seq to get here*/
	}

      /* Following block executed for all ip >= 0*/
      dmx[cur][0] = -INFTY; /* D_0 does not exist */
      mmx[cur][0] = -INFTY;
      /*for (k = cm->cp9->M-1; k >= 1; k--)*/ /*M_0 is the B state, it doesn't emit*/
      for (k = cm->cp9->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	{
	  mmx[cur][0] = ILogsum(mmx[cur][0], mmx[prv][k] + cm->cp9->bsc[k]);
	  /*printf("cm->cp9->bsc[%d]: %d | mmx[ip+1][k]: %d\n", k, cm->cp9->bsc[k], mmx[(i+1)][k]);
	    printf("mmx[%3d][%3d]: %d | k: %3d\n", i, 0, mmx[ip][0], k);*/
	}
      mmx[cur][0] = ILogsum(mmx[cur][0], imx[prv][0] + cm->cp9->tsc[CTMI][0]);
      mmx[cur][0] = ILogsum(mmx[cur][0], dmx[cur][1] + cm->cp9->tsc[CTMD][0]);

      /*printf("mmx[ip:%d][%d]: %d\n", ip, 0, mmx[cur][0]);
	printf("imx[ip:%d][%d]: %d\n", ip, 0, imx[cur][0]);
	printf("dmx[ip:%d][%d]: %d\n", ip, 0, dmx[cur][0]);*/

      /* determine isc, the int score of all possible parses starting at the current
       * position (i) of the target sequence. */
      isc[ip] = -INFTY;

      for (k = cm->cp9->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	isc[ip] = ILogsum(isc[ip], mmx[cur][k] + cm->cp9->bsc[k]); /* B->M_k */
      isc[ip] = ILogsum(isc[ip], dmx[cur][1] + cm->cp9->tsc[CTMD][0]); /* B->D_1 */
      isc[ip] = ILogsum(isc[ip], imx[cur][0] + cm->cp9->tsc[CTMI][0]); /* B->I_1 */
      fsc = Scorify(isc[ip]);
      /*printf("ip: %d fsc: %f\n\n", ip, fsc);*/

      if(fsc > best_negsc) 
	{
	  best_negsc = fsc;
	  best_negpos= i;
	}

      /* The little semi-Markov model that deals with multihit parsing:
       */
      /* Here's a vicious off-by-one with HMMs and CMs:
       * an HMM parse can begin at position 0, i.e. backward[0][0] is the probability
       * of the entire sequence being emitted starting at the Begin state at position 0.
       * But a CM always has parses that start at position 1. So we have to translate
       * from HMM to CM in this way, and that's why there's a ip+1 in the following loop.
       */
      gamma[ip+1]  = gamma[ip+2] + 0; /* extend without adding a new hit */
      /*printf("i: %d | gamma[i]: %f | gamma[i+1]: %f\n", i, gamma[i], gamma[i+1]);*/
      gback[ip+1]  = -1;
      savesc[ip+1] = IMPOSSIBLE;
      j = ((i+1+W-1) < j0) ? (i+1+W-1) : j0;
      jp = j-i0+1;
      curr_sc = gamma[jp+1] + fsc + cm->cp9_sc_boost;
      /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
       * value is 0.0 if technique not used. */
      if (curr_sc > gamma[ip+1])
	{
	  gamma[ip+1]  = curr_sc;
	  gback[ip+1]  = j;
	  savesc[ip+1] = fsc;
	}
    }
  /* End of Backward() code */

  /*****************************************************************
   * Traceback stage for Backward.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  alloc_nhits = 10;
  hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
  nhits = 0;
  i     = i0;
  while (i < j0) 
    {
      ip    = i-i0+1;
      if (gback[ip] == -1) /* no hit */
	i++; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[ip] > best_sc) 
	    {
	      best_sc = savesc[ip];
	      best_pos= i;
	    }
	  if(savesc[ip] >= cutoff)
	    {
	      accept = TRUE;
	      /* Potentially rescan just the subseq that is the hit we're about to report.
	       * Implemented to deal with fact that --enfseq option was enforcing the subseq
	       * to have hit pass filter b/c this Forward scanning function is 'infinite' length
	       * (Weinberg). Sometimes the subseq we're about to report has a really
	       * crappy score, even though the cumulative Forward score (starting at i0) is good. */
	      if(cm->search_opts & CM_SEARCH_HMMRESCAN && doing_rescan == FALSE)
		{
		  /*printf("rechecking hit from %d to %d\n", i, gback[ip]);*/
	      temp_sc = CP9ForwardScan(cm, dsq, i, gback[ip], cm->W, cutoff, 
				       NULL,  /* don't care about scores of each pos */
				       NULL,  /* don't care about locations of hits */
				       NULL,  /* don't care about num hits found */
				       NULL,  /* don't care about best scoring position */
				       NULL,  /* don't report hits */
				       TRUE); /* set the doing_rescan arg to TRUE, 
						 so we don't potentially infinitely recurse */
		  /*printf("new score: %f old score %f\n", temp_sc, savesc[ip]);*/
		  if(temp_sc >= cutoff) 
		    { 
		      accept = TRUE; 
		      /*printf("rechecked hit from %d to %d\n", i, gback[ip]);
			printf("new score: %f old score %f\n", temp_sc, savesc[ip]);*/
		    }
		  else accept = FALSE;
		  savesc[ip] = temp_sc;
		}
	      if(accept)
		{
		  if(results != NULL) /* report the hit */
		    {
		      report_hit(i, gback[ip], 0, savesc[ip], results); 
		      /* 0 is for saver, which is irrelevant for HMM hits */
		    }
		  hiti[nhits] = i;
		  nhits++;
		  if (nhits == alloc_nhits) 
		    {
		      alloc_nhits += 10;
		      hiti = ReallocOrDie(hiti,  sizeof(int)   * (alloc_nhits));
		    }
		}
	    }
	  i = gback[ip]+1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);

  /*printf("returning from CP9BackwardScan()\n");*/
  if(ret_isc != NULL) *ret_isc = isc;
  else free(isc);
  if(ret_hiti != NULL) *ret_hiti = hiti;
  else free(hiti);
  if(ret_nhits != NULL) *ret_nhits = nhits;
  if(best_sc <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    {
      best_sc = best_negsc;
      if(ret_bestpos != NULL) *ret_bestpos = best_negpos;
    }
  else if(ret_bestpos != NULL) *ret_bestpos = best_pos;
  return best_sc;
}

/***********************************************************************
 * Function: CP9Scan_dispatch()
 * Incept:   EPN, Tue Jan  9 06:28:49 2007
 * 
 * Purpose:  Scan a sequence with a CP9, potentially rescan CP9 hits with CM.
 *
 *           3 possible modes:
 *
 *           Mode 1: Filter Weinberg style (IF cm->search_opts & CM_SEARCH_HMMWEINBERG)
 *                   Scan with CP9ForwardScan() to get likely endpoints (j) of 
 *                   hits, set i for each j as j-W+1. Then, set i' = i-W+1 j'= j+W-1.
 *                   Each i'..j' CP9 subequence is then rescanned with the CM. 
 *
 *           Mode 2: Filter FB style (IF cm->search_opts & CM_SEARCH_HMMFB)
 *                   Scan with CP9ForwardScan() to get likely endpoints (j) of 
 *                   hits, for each j do a CP9BackwardScan() from j-W+1..j to get
 *                   the most likely start point i for this j. Set i' = j-W+1 and
 *                   j' = i+W-1. Each i'..j' subsequence is rescanned with the CM.
 * 
 *           Mode 3: HMM only mode (IF cm->search_opts & CM_SEARCH_HMMONLY)
 *                   Hit boundaries i and j are defined the same as in mode 2, but
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
 *           results    - scan_results_t to add to, only passed to 
 *                        actually_search_target()
 *           doing_cp9_stats- TRUE if we're calc'ing stats for the CP9, in this 
 *                            case we never rescan with CM
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc, score of maximally scoring end position j 
 */
float
CP9Scan_dispatch(CM_t *cm, char *dsq, int i0, int j0, int W, float cm_cutoff, 
		 float cp9_cutoff, scan_results_t *results, int doing_cp9_stats,
		 int *ret_flen)
{
  int  f_nhits; /* number of hits above threshold from Forward scan */
  int *f_hitj;  /* end points of hits from Forward scan */
  int  nhits;   /* number of hits */
  int *hiti;    /* start points of hits from Backward scan */
  int *hitj;    /* end points j of hits from Forward scan for which 
		 * Backward found a hit j-W+1..j with sc >= cp9_cutoff */
  float *hitsc; /* scores of hits */
  int alloc_nhits; /* for growing hit arrays */
  int h;
  int i, j;
  int min_i;
  float best_hmm_sc;
  float best_hmm_fsc;
  float best_hmm_bsc;
  float cur_best_hmm_bsc;
  float best_cm_sc;
  float cm_sc;
  int   flen;
  float ffrac;
  int do_collapse;
  int i_lpad;
  int i_rpad;
  int j_lpad;
  int j_rpad;

  /*printf("in CP9Scan_dispatch(), i0: %d j0: %d\n", i0, j0);
    printf("cp9_cutoff: %f\n", cp9_cutoff);*/

  best_cm_sc = best_hmm_sc = IMPOSSIBLE;
  /* set up options for RescanFilterSurvivors() if we're filtering */
  if(cm->search_opts & CM_SEARCH_HMMWEINBERG)
    {
      i_lpad = j_rpad = W;
      i_rpad = j_lpad = 0;
      do_collapse = TRUE;
    }
  else if(cm->search_opts & CM_SEARCH_HMMFB)
    {
      i_lpad = j_rpad = 0;
      i_rpad = j_lpad = W-1;
      do_collapse = TRUE;
    }

  /* Scan the (sub)seq w/Forward, getting j end points of hits above cutoff */
  best_hmm_fsc = CP9ForwardScan(cm, dsq, i0, j0, W, cp9_cutoff, NULL, &f_hitj, &f_nhits, NULL, NULL, FALSE);

  /* Determine start points (i) of the hits */
  if(cm->search_opts & CM_SEARCH_HMMWEINBERG)
    {
      /* i is j - W + 1 */
      hiti = MallocOrDie(sizeof(int) * f_nhits);
      hitj = MallocOrDie(sizeof(int) * f_nhits);
      hitsc = MallocOrDie(sizeof(float) * f_nhits); /* only used in mode 3 HMMONLY */
      nhits = f_nhits;
      for(h = 0; h < f_nhits; h++) 
	{
	  hitj[h] = f_hitj[h];
	  hiti[h] = (hitj[h] - W + 1) >= 1 ? (hitj[h] - W + 1) : 1;
	}
      best_hmm_sc = best_hmm_fsc;
    }
  else if((cm->search_opts & CM_SEARCH_HMMFB) || (cm->search_opts & CM_SEARCH_HMMONLY))
    {
      /* determine i based on Backward scan starting at j */

      alloc_nhits = 10 < f_nhits ? 10 : f_nhits;
      hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
      hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
      hitsc = MallocOrDie(sizeof(float) * alloc_nhits);

      nhits = 0;
      for(h = 0; h < f_nhits; h++) 
	{
	  min_i = (f_hitj[h] - W + 1) >= 1 ? (f_hitj[h] - W + 1) : 1;
	  cur_best_hmm_bsc = CP9BackwardScan(cm, dsq, min_i, f_hitj[h], W, cp9_cutoff, 
					     NULL, /* don't care about score of each posn */
					     NULL, /* don't care about locations of hits */
					     NULL, /* don't care about number of hits */
					     &i,   /* set i as the best scoring start point from j-W..j */
					     NULL,  /* don't report hits */
					     FALSE); /* don't rescan */
	  if(cur_best_hmm_bsc > best_hmm_sc) best_hmm_sc = cur_best_hmm_bsc;
	  if(cur_best_hmm_bsc >= cp9_cutoff)
	    {
	      /* this is a real hit, Backward found hit from j-W+1..j >= cutoff */
	      hiti[nhits]  = i;
	      hitj[nhits]  = f_hitj[h];
	      hitsc[nhits] = cur_best_hmm_bsc; /* only used if CM_SEARCH_HMMONLY */
	      nhits++;
	      if (nhits == alloc_nhits) 
		{
		  alloc_nhits += 10;
		  hiti  = ReallocOrDie(hiti,  sizeof(int)   * alloc_nhits);
		  hitj  = ReallocOrDie(hitj,  sizeof(int)   * alloc_nhits);
		  hitsc = ReallocOrDie(hitsc, sizeof(float) * alloc_nhits);
		}
	    }
	}	  
    }
  else Die("ERROR neither filter mode selected and CM_SEARCH_HMMONLY is false.\n");

  /* Rescan with CM if we're filtering and not doing cp9 stats */
  if(!doing_cp9_stats && ((cm->search_opts & CM_SEARCH_HMMWEINBERG) ||
			  (cm->search_opts & CM_SEARCH_HMMFB)))
    {
      best_cm_sc = RescanFilterSurvivors(cm, dsq, hiti, hitj, nhits, i0, j0, W, 
					 i_lpad, i_rpad, j_lpad, j_rpad,
					 do_collapse, cm_cutoff, cp9_cutoff, 
					 results, ret_flen);
    }
  else if(cm->search_opts & CM_SEARCH_HMMONLY)
    {
      /* report hits as i,j pairs i from CP9ForwardScan() and j from CP9BackwardScan() */
      for(h = 0; h <= nhits-1; h++) 
	{
	  report_hit(hiti[h], hitj[h], 0, hitsc[h], results); 
	  /* 0 is for saver, which is irrelevant for HMM hits */
	}
    }
  free(hiti);
  free(hitj);
  free(hitsc);
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
 *           Can be run in 2 modes, depending on input variables: i_lpad,
 *           i_rpad, j_lpad, j_rpad. 
 *           Mode 1: i_lpad == j_rpad == 0
 *                   For each i,j pair in hiti, hitj: 
 *                   set j' = i + i_rpad; and i' = j - j_lpad, 
 *                   ensure j' >= j and i' <= i. 
 *                   Rescan subseq from i' to j'.
 *           Mode 2: i_rpad == j_lpad == 0
 *                   For each i,j pair in hiti, hitj: 
 *                   set i' = i - i_lpad; and j' = j + j_rpad, 
 *                   Rescan subseq from i' to j'.
 *
 * Args:     
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           hiti       - start of subseqs that survived filter [0..nhits-1]
 *           hitj       - end of subseqs that survived filter   [0..nhits-1]
 *           nhits      - number of hits that survived filter
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           i_lpad     - number of residues to subtract from each i 
 *           i_rpad     - number of residues to add      to   each i to get j
 *           j_lpad     - number of residues to subtract from each j to get i
 *           j_rpad     - number of residues to add      to   each j
 *           do_collapse- TRUE: collapse overlapping hits (after padding) into 1
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report 
 *           results    - scan_results_t to add to, only passed to 
 *                        actually_search_target()
 *           ret_flen   - RETURN: subseq len that survived filter
 * Returns:  best_sc found when rescanning with CM 
 */
float
RescanFilterSurvivors(CM_t *cm, char *dsq, int *hiti, int *hitj, int nhits, int i0, int j0,
		      int W, int i_lpad, int i_rpad, int j_lpad, int j_rpad, int do_collapse,
		      float cm_cutoff, float cp9_cutoff, scan_results_t *results, int *ret_flen)
{
  int h;
  int i, j;
  float best_cm_sc;
  float cm_sc;
  int   flen;
  float ffrac;
  int   mode;
  int   prev_j;
  int   next_j;

  best_cm_sc = IMPOSSIBLE;
  flen = 0;

  /* check contract */
  if(i_lpad == 0 && j_rpad == 0)      mode = 1; /* it's possible all 4 (or 3) are zero in mode 1 */
  else if(i_rpad == 0 && j_lpad == 0) mode = 2;
  else ESL_EXCEPTION(eslEINCOMPAT, "can't determine mode, neither is true: i_lpad == j_rpad == 0; i_rpad == j_lpad == 0.");

  /*printf("in RescanFilterSurvivors(), mode: %d i_lpad: %d, i_rpad: %d, j_lpad: %d, j_rpad: %d collapse: %d\n", mode, i_lpad, i_rpad, j_lpad, j_rpad, do_collapse);*/

  /* For each hit, add pad according to mode and rescan by calling actually_search_target(). 
   * If do_collapse, collapse multiple overlapping hits into 1 before rescanning */
  /* hits should always be sorted by decreasing j, if this is violated - die. */
  for(h = 0; h <= nhits-1; h++) 
    {
      if(h != 0 && hitj[h] > prev_j) 
	ESL_EXCEPTION(eslEINCOMPAT, "j's not in descending order");
      prev_j = hitj[h];

      /* add pad */
      if(mode == 1)
	{
	  i = ((hitj[h] - j_lpad) >= 1)    ? (hitj[h] - j_lpad) : 1;
	  j = ((hiti[h] + i_rpad) <= j0)   ? (hiti[h] + i_rpad) : j0;
	  if((h+1) < nhits)
	    next_j = ((hiti[h+1] + i_rpad) <= j0)   ? (hiti[h+1] + i_rpad) : j0;
	  else
	    next_j = -1;
	}
      else /* mode == 2 */
	{
	  i = ((hiti[h] - i_lpad) >= 1)    ? (hiti[h] - i_lpad) : 1;
	  j = ((hitj[h] + j_rpad) <= j0)   ? (hitj[h] + j_rpad) : j0;
	  if((h+1) < nhits)
	    next_j = ((hitj[h+1] + j_rpad) <= j0)   ? (hitj[h+1] + j_rpad) : j0;
	  else
	    next_j = -1;
	}
      /*printf("subseq: hit %d i: %d (%d) j: %d (%d)\n", h, i, hiti[h], j, hitj[h]);*/
      while(((h+1) < nhits) && (next_j >= i))
      {
	/* suck in hit */
	h++;
	if(mode == 1) 
	  {
	    i = ((hitj[h] - j_lpad) >= 1)    ? (hitj[h] - j_lpad) : 1;
	    if((h+1) < nhits)
	      next_j = ((hiti[h+1] + i_rpad) <= j0)   ? (hiti[h+1] + i_rpad) : j0;
	    else
	      next_j = -1;
	  }
	if(mode == 2)
	  {
	    i = ((hiti[h] - i_lpad) >= 1)    ? (hiti[h] - i_lpad) : 1;
	    if((h+1) < nhits)
	      next_j = ((hitj[h+1] + j_rpad) <= j0)   ? (hitj[h+1] + j_rpad) : j0;
	    else
	      next_j = -1;
	  }
	/*printf("\tsucked in subseq: hit %d new_i: %d j (still): %d\n", h, i, j);*/
      }
      /*printf("calling actually_search_target: %d %d h: %d nhits: %d\n", i, j, h, nhits);*/
      cm_sc =
	actually_search_target(cm, dsq, i, j, cm_cutoff, cp9_cutoff,
			       results, /* keep results                                 */
			       FALSE,   /* don't filter, we already have                */
			       FALSE,   /* we're not building a histogram for CM stats  */
			       FALSE,   /* we're not building a histogram for CP9 stats */
			       NULL);   /* filter fraction N/A                          */
      flen += (j - i + 1);
      if(cm_sc > best_cm_sc) best_cm_sc = cm_sc;
    }

  if(flen == 0) ffrac = 100.;
  else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
  /*printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);*/
  if(ret_flen != NULL) *ret_flen = flen;
  return best_cm_sc;
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
  /*printf("fb_sc: %f\n", Scorify(fb_sum));*/
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
 * Purpose:  Given likely end points of hits from a CP9ForwardScan() and
 *           likely end points from a CP9BackwardScan() combine them in 
 *           a greedy manner to get subseqs to research with CM.
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
 * Function: CP9FilteredWeinbergScan()
 * Incept:   EPN, Tue Jan  9 06:28:49 2007
 * 
 * Purpose:  Scan a sequence with a CP9, filtering for promising subseqs
 *           that could be hits to the CM the CP9 was derived from.
 *           Pass filtered subseqs to actually_search_target() to be
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
CP9FilteredWeinbergScan(CM_t *cm, char *dsq, int i0, int j0, int W, float cm_cutoff, 
			float cp9_cutoff, scan_results_t *results, int *ret_flen)
{
  int *hitj;
  int  nhits;
  int h;
  int i, j;
  float best_hmm_sc;
  float best_cm_sc;
  float cm_sc;
  int   flen;
  float ffrac;
  /*printf("in CP9FilteredWeinbergScan(), i0: %d j0: %d\n", i0, j0);*/
  /*printf("cp9_cutoff: %f\n", cp9_cutoff);*/
  /* Scan the (sub)seq w/Forward, getting j end points of hits above cutoff */
  best_hmm_sc = CP9ForwardScan(cm, dsq, i0, j0, W, cp9_cutoff, NULL, &hitj, &nhits, NULL, NULL, FALSE);
  /* Send promising subseqs to actually_search_target(): send subseq [j-W..j..j+W] 
   * for each hit endpoint j returned from CP9ForwardScan, or if there's overlap, send 
   * minimal subseq that encompasses all overlapping hits with W residue 'pad' on 
   * both sides.
   */
  best_cm_sc = IMPOSSIBLE;
  flen = 0;
  /* hits are always sorted by decreasing j */
  for(h = 0; h <= nhits-1; h++) 
    {
      j = ((hitj[h] + W) <= j0) ? (hitj[h] + W) : j0;
      i = ((j - (2*W)) >= 1)    ? (j - (2*W))   : 1;
      /*printf("subseq: hit %d j: %d i: %d j: %d\n", h, hitj[h], i, j);*/
      while(((h+1) < nhits) && ((hitj[(h+1)]+W) >= i))
      {
	/* suck in hit */
	h++;
	i = ((hitj[h]-W) >= 1) ? (hitj[h]-W) : 1;
	/*printf("\tsucked in subseq: hit %d new_i: %d j (still): %d\n", h, i, j);*/
      }
      /*printf("calling actually_search_target: %d %d h: %d nhits: %d\n", i, j, h, nhits);*/
      cm_sc =
	actually_search_target(cm, dsq, i, j, cm_cutoff, cp9_cutoff,
			       results, /* keep results                                 */
			       FALSE,   /* don't filter, we already have                */
			       FALSE,   /* we're not building a histogram for CM stats  */
			       FALSE,   /* we're not building a histogram for CP9 stats */
			       NULL);   /* filter fraction N/A                          */
      flen += (j - i + 1);
      if(cm_sc > best_cm_sc) best_cm_sc = cm_sc;
    }

  free(hitj);
  if(flen == 0) ffrac = 100.;
  else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
  /*printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);*/
  if(ret_flen != NULL) *ret_flen = flen;
  return best_cm_sc;
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
 * Args:     
 *           fmx    - forward dp matrix, already filled
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
  printf("sc: %f\n", sc);

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
      diff = sc - fb_sc;
      if(diff < 0.) diff *= -1.;
      if(diff > max_diff)
	{
	  printf("ERROR, fb_sc[%d]: %f too different from P(x|hmm): %f\n", i, fb_sc, sc);
	  exit(1);
	}
    }
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
 * Args: 
 *           cm         - the covariance model, includes cm->cp9: a CP9 HMM
 *           dsq        - sequence in digitized form
 *           i0         - start of target subsequence (1 for beginning of dsq)
 *           j0         - end of target subsequence (L for end of dsq)
 *           W          - the maximum size of a hit (often cm->W)
 *           cp9_cutoff - minimum CP9 score to report (or keep if filtering)
 *           ret_hitsc  - RETURN: scores of hits, 0..nhits-1            
 *           ret_nhits  - RETURN: number of hits
 *           ret_hiti   - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj   - RETURN: end positions of hits, 0..nhits-1
 *           results    - scan_results_t to add to, only passed to 
 *                        actually_search_target()
 *           doing_rescan - TRUE if we've called this function recursively, and we don't
 *                          want to again
 *           pad        - number of nucleotides to add on to start and end points
 *
 * Returns:  
 *           hiti, hitj, hitsc are allocated in a helper function; caller free's w/ free().
 */
float
CP9ForwardBackwardScan(CM_t *cm, char *dsq, int i0, int j0, int W, 
		       int **ret_hitsc, int *ret_nhits, int **ret_hiti, int **ret_hitj, 
		       scan_results_t *results, int doing_rescan, int pad)
{
  int      L;                   /* j0-i0+1: subsequence length */
  int      jp;		        /* j': relative position in the subsequence  */
  int      ip;		        /* i': relative position in the subsequence  */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int         *isc;         /* prob (seq from j0..jp | HMM) [0..jp..cm->cp9->M]             */
  int   i,j,k;
  float     fsc;     
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int       alloc_nhits;	/* used to grow the hit arrays */

  int fwd_sum;                  /* calc'ed from the Forward matrix, the sum over all
				 * positions 1..L of all parses that start and end 
				 * anywhere from 1..L */
  int bck_sum;                  /* same as above, but calc'ed from the Backward matrix*/
  float fwd_sc;                 /* fwd_sum as a score */
  float bck_sc;                 /* bck_sum as a score */
  
  int fwd_nhits;                /* number of hits from Forward */
  int *fwd_hiti;                /* start positions of hits from Forward, 0..nhits-1, (calc'ed as end points - W) */
  int *fwd_hitj;                /* end positions of hits from Forward, 0..nhits-1 */
  float *fwd_hitsc;             /* scores of hits, 0..nhits-1 */

  int bck_nhits;                /* number of hits from Backward */
  int *bck_hiti;                /* start positions of hits from Backward, 0..nhits-1 */
  int *bck_hitj;                /* end positions of hits from Backward, 0..nhits-1 (calc'ed as start points + W) */
  float *bck_hitsc;             /* scores of hits, 0..nhits-1 */

  /* Contract checks */
  if(cm->cp9 == NULL)
    Die("ERROR in CP9ForwardBackwardScan, cm->cp9 is NULL.\n");

  best_sc     = IMPOSSIBLE;
  best_negsc  = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  fwd_sum = -INFTY;
  bck_sum = -INFTY;

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
			

  /*****************************************************************
   * The main Forward loop: scan the sequence from position 1 to L.
   *****************************************************************/
  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */

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


/***********************************************************************
 * Function: CP9FilteredFBScan()
 * Incept:   EPN, Tue Apr 10 06:25:09 2007
 * 
 * Purpose:  Scan a sequence with a CP9, using Forward() to get end points
 *           and Backward() to get startpoints of CP9 HMM hits that could 
 *           hits to the CM the CP9 was derived from.
 *           Pass filtered subseqs to actually_search_target() to be
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
CP9FilteredFBScan(CM_t *cm, char *dsq, int i0, int j0, int W, float cm_cutoff, 
		  float cp9_cutoff, scan_results_t *results, int *ret_flen)
{
  int *hiti;
  int *hitj;
  int  nhits;
  int h;
  int i, j;
  float best_hmm_fsc;
  float best_hmm_bsc;
  float best_cm_sc;
  float cm_sc;
  int   flen;
  float ffrac;

  printf("in CP9FilteredFBScan(), i0: %d j0: %d\n", i0, j0);
  printf("cp9_cutoff: %f\n", cp9_cutoff);
  /* Scan the (sub)seq w/Forward getting likely j end points */
  best_hmm_fsc = CP9ForwardScan(cm, dsq, i0, j0, W, cp9_cutoff, 
				NULL,  /* don't care about score of each posn */
				&hitj, /* keep track of locations of hits */
				&nhits,/* keep track of number of hits */
				NULL,  /* don't care about best scoring start point */
				NULL,  /* don't report hits */
				FALSE);/* don't rescan */
  
  /* OLD STRATEGY AS OF 04.10.07 */
  /* Combine likely start and end points with greedy strategy */
  /*  CP9_combine_FBscan_hits(i0, j0, W, fwd_nhits, fwd_hitr, fwd_hiti, fwd_hitj, fwd_hitsc, 
      bck_nhits, bck_hitr, bck_hiti, bck_hitj, bck_hitsc,
      ret_nhits, ret_hitr, ret_hiti, ret_hitj, ret_hitsc, pad);*/

  /* NEW STRATEGY AS OF 04.10.07 */
  /* For each likely end point j, find the most likely start point with CP9BackwardScan() */
  for(h = 0; h <= nhits-1; h++) 
    {
      best_hmm_bsc = CP9BackwardScan(cm, dsq, (hitj[h]-W+1), hitj[h], W, cp9_cutoff, 
				     NULL, /* don't care about score of each posn */
				     NULL, /* don't care about locations of hits */
				     NULL, /* don't care about number of hits */
				     &i,   /* set i as the best scoring start point from j-W..j */
				     NULL,  /* don't report hits */
				     FALSE); /* don't rescan */
      if(best_hmm_bsc > cp9_cutoff)
	{
	  /* Rescan this hit with the CM. 
	   * OPTIONS FOR IMPLEMENTATION.
	   * 1. pad 'pad' residues on left of i and right of j
	   * 2. if do_hbanded_scan, first get full scanning matrices 
	   *    (this will require new code can I modify CP9ForwardScan 
	   *     and CP9BackwardScan() to keep matrix?), 
	   *    then call (after modifying) CP9ScanFullPosterior to get
	   *    bands, then make actually_search_target optionally use
	   *    HMM bands.
	   */
	  cm_sc =
	    actually_search_target(cm, dsq, i, j, cm_cutoff, cp9_cutoff,
				   results, /* keep results                                 */
				   FALSE,   /* don't filter, we already have                */
				   FALSE,   /* we're not building a histogram for CM stats  */
				   FALSE,   /* we're not building a histogram for CP9 stats */
				   NULL);   /* filter fraction N/A                          */
	  flen += (j - i + 1);
	  if(cm_sc > best_cm_sc) best_cm_sc = cm_sc;
	}
    }
  free(hitj);
  if(flen == 0) ffrac = 100.;
  else ffrac = 1. - (((float) flen) / (((float) (j0-i0+1))));
  /*printf("orig_len: %d flen: %d fraction %6.2f\n", (j0-i0+1), (flen), ffrac);*/
  if(ret_flen != NULL) *ret_flen = flen;
  return best_cm_sc;
}

#endif

/* CP9scan.c 
 * 
 * Scanning algorithms for CM Plan 9 HMMs.
 * These algorithms align to the subsequences of the target
 * sequence to the model (e.g. glocal or local alignment)
 * Global alignment algorithms are in hmmband.c.
 *
 *################################################################
 * CP9ForwardScan()      - Scan input sequence for high scoring
 *                          forward hits to the model .
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
#include "hmmer_funcs.h"
#include "hmmer_structs.h"
#include "hmmband.h"

static void debug_check_CP9_FBscan(struct cp9_dpmatrix_s *fmx, struct cp9_dpmatrix_s *bmx, 
				   struct cplan9_s *hmm, float sc, int L, unsigned char *dsq);

/***********************************************************************
 * Function: CP9ForwardScan()
 * 
 * Purpose:  A Forward dynamic programming algorithm that scans an input sequence.
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
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_fmx - the Forward scanning matrix
 *           ret_nhits - RETURN: number of hits
 *           ret_hitr  - RETURN: start states of hits, 0..nhits-1
 *           ret_hiti  - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj  - RETURN: end positions of hits, 0..nhits-1
 *           ret_hitsc - RETURN: scores of hits, 0..nhits-1            
 *           min_thresh- minimum score to report (EPN via Alex Coventry 03.11.06)
 *
 * Returns:  
 *           hiti, hitj, hitsc are allocated here; caller free's w/ free().
 */
float
CP9ForwardScan(unsigned char *dsq, int L, int W, struct cplan9_s *hmm, struct cp9_dpmatrix_s **ret_mx,
	       int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc, 
	       float min_thresh)
{
  struct cp9_dpmatrix_s *mx;
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

  int       nhits;		/* # of hits in optimal parse */
  int      *hitr;		/* initial state indices of hits in optimal parse */
  int      *hiti;		/* initial state indices of hits in optimal parse */
  int      *hitj;               /* end positions of hits in optimal parse */
  float    *hitsc;              /* scores of hits in optimal parse */
  int       alloc_nhits;	/* used to grow the hit arrays */

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
  mx = AllocCPlan9Matrix(L+1, hmm->M, &mmx, &imx, &dmx, &emx);

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
  
  emx[0][0] = dmx[0][hmm->M] + hmm->tsc[CTDM][hmm->M]; 


  
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
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
  
  emx[0][0] = dmx[0][hmm->M] + hmm->tsc[CTDM][hmm->M]; 
  
  /* Recursion. Done as a pull.
   */
  for (j = 1; j <= L; j++)
    {
      //cur = j%2;
      //prv = (j-1)%2;
      mmx[j][0] = 0;  /* This is the 1 difference between a Forward scanner and the 
		       * regular Forward (CP9Forward), 
		       * in CP9Forward, this cell is set to -INFTY.
		       * Here, we can start at any position in the seq,
		       * this isn't true in global alignment. 
		       */
      dmx[j][0] = -INFTY;  /*D_0 is non-existent*/
      imx[j][0]  = ILogsum(ILogsum(mmx[j-1][0] + hmm->tsc[CTMI][0],
				   imx[j-1][0] + hmm->tsc[CTII][0]),
			   dmx[j-1][0] + hmm->tsc[CTDI][0]);
      imx[j][0] += hmm->isc[dsq[j]][0];
      
      for (k = 1; k <= hmm->M; k++)
	{
	  mmx[j][k]  = ILogsum(ILogsum(mmx[j-1][k-1] + hmm->tsc[CTMM][k-1],
				       imx[j-1][k-1] + hmm->tsc[CTIM][k-1]),
			       ILogsum(mmx[j-1][0] + hmm->bsc[k],
				       dmx[j-1][k-1] + hmm->tsc[CTDM][k-1]));
	  mmx[j][k] += hmm->msc[dsq[j]][k];
	  
	  dmx[j][k]  = ILogsum(ILogsum(mmx[j][k-1] + hmm->tsc[CTMD][k-1],
				       imx[j][k-1] + hmm->tsc[CTID][k-1]),
			       dmx[j][k-1] + hmm->tsc[CTDD][k-1]);
	  
	  imx[j][k]  = ILogsum(ILogsum(mmx[j-1][k] + hmm->tsc[CTMI][k],
				       imx[j-1][k] + hmm->tsc[CTII][k]),
			       dmx[j-1][k] + hmm->tsc[CTDI][k]);
	  imx[j][k] += hmm->isc[dsq[j]][k];
	  //printf("mmx[%d][%d]: %d\n", i, k, mmx[j][k]);
	  //printf("imx[%d][%d]: %d\n", i, k, imx[j][k]);
	  //printf("dmx[%d][%d]: %d\n", i, k, dmx[j][k]);

	}
      
      emx[0][j] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	emx[0][j] = ILogsum(emx[0][j], mmx[j][k] + hmm->esc[k]);
      emx[0][j] = ILogsum(emx[0][j], dmx[j][hmm->M] + hmm->tsc[CTDM][hmm->M]); 
      emx[0][j] = ILogsum(emx[0][j], imx[j][hmm->M] + hmm->tsc[CTIM][hmm->M]); 
      /* transition from D_M -> end */
      sc = Scorify(emx[0][j]);
      if(sc > 0)
	{
	  //printf("%d %10.2f\n", j, sc);
	}

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[j]  = gamma[j-1] + 0; /* extend without adding a new hit */
      gback[j]  = -1;
      savesc[j] = IMPOSSIBLE;
      saver[j]  = -1;
      i = ((j-W+1)> 0) ? (j-W+1) : 1;
      sc = gamma[i-1] + Scorify(emx[0][j]) - min_thresh;
      /*printf("j: %d | gamma[j]: %f | sc: %f\n", j, gamma[j], sc);*/
      if (sc > gamma[j])
	{
	  gamma[j]  = sc;
	  gback[j]  = i;
	  savesc[j] = Scorify(emx[0][j]);
	  saver[j]  = 0;
	  //saver[j]  = bestr[d];
	}
    } /* end loop over end positions j */
  sc = Scorify(emx[0][L]);
  //printf("END %d %10.2f\n", j, sc);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  alloc_nhits = 10;
  hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitsc = MallocOrDie(sizeof(float) * alloc_nhits);

  j     = L;
  nhits = 0;
  while (j > 0) {
    if (gback[j] == -1) /* no hit */
      j--; 
    else                /* a hit, a palpable hit */
      {
	hitr[nhits]   = saver[j];
	hitj[nhits]   = j;
	hiti[nhits]   = gback[j];
	hitsc[nhits]  = savesc[j];
	nhits++;
	j = gback[j]-1;
	
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
  return sc;		/* the total Forward score. */
}
  
/* Function: CP9ForwardScanRequires()  EPN 04.21.06
 * based on  CYKScanRequires() by:
 * Date:     SRE, Mon May  6 18:48:18 2002 [St. Louis]
 *
 * Purpose:  Return the memory required by CP9ForwardScan(), in megabytes.
 */
float
CP9ForwardScanRequires(struct cplan9_s *hmm, int L, int W)
{
  float ram;
				/* alpha allocations */
  ram = 1000000. * SizeCPlan9Matrix(L, hmm->M);
  ram += (float) (sizeof(float) * (L+1)); /* gamma allocation */
  ram += (float) (sizeof(int)   * (L+1)); /* gback allocation */
  ram += (float) (sizeof(float) * (L+1)); /* savesc allocation */
  return (ram / 1000000.);
}

/***********************************************************************
 * Function: CP9ForwardBackwardScan()
 * 
 * Purpose:  Runs the Forward dynamic programming algorithm that scans an 
 *           input sequence. And then the Backward DP algorithm on the
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
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_fmx - the Forward scanning matrix
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
CP9ForwardBackwardScan(unsigned char *dsq, int L, int W, struct cplan9_s *hmm, struct cp9_dpmatrix_s **ret_fmx,
		       struct cp9_dpmatrix_s **ret_bmx, int *ret_nhits, int **ret_hitr, int **ret_hiti, 
		       int **ret_hitj, float **ret_hitsc, float min_thresh)
{
  struct cp9_dpmatrix_s *fmx;
  struct cp9_dpmatrix_s *bmx;
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
  for (j = 1; j <= L; j++)
    {
      //cur = j%2;
      //prv = (j-1)%2;
      fmmx[j][0] = 0; /* This is the 1 difference between a Forward scanner and the 
		       * regular Forward (CP9Forward), 
		       * in CP9Forward, this cell is set to -INFTY.
		       * Here, we can start at any position in the seq,
		       * this isn't true in global alignment. 
		       */
      fdmx[j][0] = -INFTY;  /*D_0 is non-existent*/
      fimx[j][0]  = ILogsum(fmmx[j-1][0] + hmm->tsc[CTMI][0],
			    fimx[j-1][0] + hmm->tsc[CTII][0]);
      fimx[j][0] += hmm->isc[dsq[j]][0];
      
      for (k = 1; k <= hmm->M; k++)
	{
	  fmmx[j][k]  = ILogsum(ILogsum(fmmx[j-1][k-1] + hmm->tsc[CTMM][k-1],
				       fimx[j-1][k-1] + hmm->tsc[CTIM][k-1]),
			       ILogsum(fmmx[j-1][0] + hmm->bsc[k],
				       fdmx[j-1][k-1] + hmm->tsc[CTDM][k-1]));
	  /*if(k == 3 && j == L) printf("\tk: %d | fmmx[j-1][k-1]: %d | hmm->tsc[CTMM][k-1]: %d\n\t\tfimx[j-1][k-1]: %d | hmm->tsc[CTIM][k-1]: %d\n\t\tfmmx[j-1][0]: %d | bsc: %d\n\t\tfdmx[j-1][k-1]: %d | hmm->tsc[CTDM][k-1]: %d\n", k, fmmx[(j-1)][(k-1)], hmm->tsc[CTMM][(k-1)], fimx[(j-1)][(k-1)], hmm->tsc[CTIM][(k-1)], fmmx[(j-1)][0], hmm->bsc[k], fdmx[(j-1)][(k-1)], hmm->tsc[CTDM][(k-1)]);*/
	  fmmx[j][k] += hmm->msc[dsq[j]][k];
	  
	  fdmx[j][k]  = ILogsum(ILogsum(fmmx[j][k-1] + hmm->tsc[CTMD][k-1],
				       fimx[j][k-1] + hmm->tsc[CTID][k-1]),
			       fdmx[j][k-1] + hmm->tsc[CTDD][k-1]);
	  
	  fimx[j][k]  = ILogsum(ILogsum(fmmx[j-1][k] + hmm->tsc[CTMI][k],
				       fimx[j-1][k] + hmm->tsc[CTII][k]),
			       fdmx[j-1][k] + hmm->tsc[CTDI][k]);
	  fimx[j][k] += hmm->isc[dsq[j]][k];
	  //printf("fmmx[%d][%d]: %d\n", i, k, fmmx[j][k]);
	  //printf("fimx[%d][%d]: %d\n", i, k, fimx[j][k]);
	  //printf("fdmx[%d][%d]: %d\n", i, k, fdmx[j][k]);
	}
      
      femx[0][j] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	{
	  femx[0][j] = ILogsum(femx[0][j], fmmx[j][k] + hmm->esc[k]);
	  /*if(j == L) printf("k: %d | femx[0][L]: %d | fmmx[L][k]: %d | hmm->esc[k]: %d\n", k, femx[0][L], fmmx[j][k], hmm->esc[k]);*/
	}
      femx[0][j] = ILogsum(femx[0][j], fdmx[j][hmm->M] + hmm->tsc[CTDM][hmm->M]); 
      /*if(j == L) printf("added D trans | femx[0][L]: %d\n", femx[0][L]);*/
      /* transition from D_M -> end */
      femx[0][j] = ILogsum(femx[0][j], fimx[j][hmm->M] + hmm->tsc[CTIM][hmm->M]); 
      /*if(j == L) printf("added I trans | femx[0][L]: %d\n", femx[0][L]);*/
      /* transition from I_M -> end */

      fwd_sc = Scorify(femx[0][j]);
      /*printf("femx[0][%3d]: %10.5f\n", j, fwd_sc);*/

      fwd_sum = ILogsum(fwd_sum, femx[0][j]);

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[j]  = gamma[j-1] + 0; /* extend without adding a new hit */
      gback[j]  = -1;
      savesc[j] = IMPOSSIBLE;
      saver[j]  = -1;
      i = ((j-W+1)> 0) ? (j-W+1) : 1;
      sc = gamma[i-1] + Scorify(femx[0][j]) - min_thresh;
      /*printf("j: %d | gamma[j]: %f | sc: %f\n", j, gamma[j], sc);*/
      if (sc > gamma[j])
	{
	  gamma[j]  = sc;
	  gback[j]  = i;
	  savesc[j] = Scorify(femx[0][j]);
	  saver[j]  = 0;
	  //saver[j]  = bestr[d];
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

  j     = L;
  fwd_nhits = 0;
  while (j > 0) {
    if (gback[j] == -1) /* no hit */
      j--; 
    else                /* a hit, a palpable hit */
      {
	fwd_hitr[fwd_nhits]   = saver[j];
	fwd_hitj[fwd_nhits]   = j;
	fwd_hiti[fwd_nhits]   = gback[j];
	fwd_hitsc[fwd_nhits]  = savesc[j];
	fwd_nhits++;
	/*printf("fwd_nhits: %d\n", fwd_nhits);*/
	j = gback[j]-1;
	
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
  bmmx[L][hmm->M] += hmm->msc[dsq[L]][hmm->M]; /* ... + emitted match symbol */
  bimx[L][hmm->M] = bemx[0][L] + hmm->tsc[CTIM][hmm->M];   /* I_M(C)<-E ... */
  bimx[L][hmm->M] += hmm->isc[dsq[L]][hmm->M];           /* ... + emitted match symbol */
  bdmx[L][hmm->M] = bemx[0][L] + hmm->tsc[CTDM][hmm->M];    /* D_M<-E */
  for (k = hmm->M-1; k >= 1; k--)
    {
      bmmx[L][k]  = bemx[0][L] + hmm->esc[k];
      bmmx[L][k]  = ILogsum(bmmx[L][k], bdmx[L][k+1] + hmm->tsc[CTMD][k]);
      bmmx[L][k] += hmm->msc[dsq[L]][k];

      bimx[L][k] = bdmx[L][k+1] + hmm->tsc[CTID][k];
      bimx[L][k] += hmm->isc[dsq[L]][k];

      bdmx[L][k] = bdmx[L][k+1] + hmm->tsc[CTDD][k];
    }
  
  bmmx[L][0] = -INFTY; /* no esc[0] b/c its impossible */
  bimx[L][0] = bdmx[L][1] + hmm->tsc[CTID][0];
  bimx[L][0] += hmm->isc[dsq[L]][hmm->M];    
  bdmx[L][0] = -INFTY; /*D_0 doesn't exist*/

  bck_sum = ILogsum(bck_sum, bmmx[L][0]);

  bck_sc = Scorify(bmmx[L][0]);
  /*printf("bmmx[L][  0]: %10.5f\n", bck_sc);*/
  /*printf("tmp_bck_sum L: %d\n", tmp_bck_sum);*/

  /* Recursion. Done as a pull.
   */
  for (i = L-1; i >= 0; i--)
    {
      if(i > 0)
	{
	  bemx[0][i] = 0; /* This is 1 of 3 differences between a Backward scanner and the 
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
	  bmmx[i][hmm->M] = ILogsum(bemx[0][i+1] + hmm->esc[hmm->M], /* M<-E ... 2nd diff w/non-scanner*/
				    bimx[i+1][hmm->M] + hmm->tsc[CTMI][hmm->M]);
	  bmmx[i][hmm->M] += hmm->msc[dsq[i]][hmm->M];
	  bimx[i][hmm->M] = ILogsum(bemx[0][i] + hmm->tsc[CTIM][hmm->M],    /* I_M(C)<-E ... */
				    bimx[i+1][hmm->M] + hmm->tsc[CTII][hmm->M]);
	  bimx[i][hmm->M] += hmm->isc[dsq[i]][hmm->M];
	  bdmx[i][hmm->M] = ILogsum(bemx[0][i] + hmm->tsc[CTDM][hmm->M], /* D_M<-E */
				    bimx[i+1][hmm->M] + hmm->tsc[CTDI][hmm->M]);  
	  for (k = hmm->M-1; k >= 1; k--)
	    {
	      bmmx[i][k]  = ILogsum(ILogsum(bemx[0][i+1] + hmm->esc[k], /* M<-E ... 2nd diff w/non-scanner*/
					    bmmx[i+1][k+1] + hmm->tsc[CTMM][k]),
				    ILogsum(bimx[i+1][k] + hmm->tsc[CTMI][k],
					    bdmx[i][k+1] + hmm->tsc[CTMD][k]));	  
	      /*if(k == 1 && i == 1) printf("\tk: %d | bemx[0][i+1]: %d | hmm->esc[k]: %d\n\t\tbmmx[i+1][k+1]: %d | hmm->tsc[CTMM][k]: %d\n\t\tbimx[i+1][k]: %d | hmm->tsc[CTMI][k]: %d\n\t\tbdmx[i][k+1]: %d | hmm->tsc[CTMD][k]: %d\n", k, bemx[0][(i+1)], hmm->esc[k], bmmx[(i+1)][(k+1)], hmm->tsc[CTMM][k], bimx[(i+1)][k], hmm->tsc[CTMI][k], bdmx[i][(k+1)], hmm->tsc[CTMD][k]);*/
	      
	      bmmx[i][k] += hmm->msc[dsq[i]][k];
	      
	      bimx[i][k]  = ILogsum(ILogsum(bmmx[i+1][k+1] + hmm->tsc[CTIM][k],
					    bimx[i+1][k] + hmm->tsc[CTII][k]),
				    bdmx[i][k+1] + hmm->tsc[CTID][k]);
	      bimx[i][k] += hmm->isc[dsq[i]][k];
	      
	      bdmx[i][k]  = ILogsum(ILogsum(bmmx[i+1][k+1] + hmm->tsc[CTDM][k],
					    bimx[i+1][k] + hmm->tsc[CTDI][k]),
				    bdmx[i][k+1] + hmm->tsc[CTDD][k]);
	    }
	  
	  bimx[i][0]  = ILogsum(ILogsum(bmmx[i+1][1] + hmm->tsc[CTIM][0],
					bimx[i+1][0] + hmm->tsc[CTII][0]),
				bdmx[i][1] + hmm->tsc[CTID][0]);
	  bimx[i][0] += hmm->isc[dsq[i]][0];
	  bmmx[i][0] = -INFTY;
	  /*for (k = hmm->M-1; k >= 1; k--)*/ /*M_0 is the B state, it doesn't emit*/
	  for (k = hmm->M; k >= 1; k--) /*M_0 is the B state, it doesn't emit*/
	    {
	      bmmx[i][0] = ILogsum(bmmx[i][0], bmmx[i+1][k] + hmm->bsc[k]);
	      /*printf("hmm->bsc[%d]: %d | bmmx[i+1][k]: %d\n", k, hmm->bsc[k], bmmx[(i+1)][k]);
		printf("bmmx[%3d][%3d]: %d | k: %3d\n", i, 0, bmmx[i][0], k);*/
	    }
	  bmmx[i][0] = ILogsum(bmmx[i][0], bimx[i+1][k] + hmm->tsc[CTMI][k]);
	  bmmx[i][0] = ILogsum(bmmx[i][0], bdmx[i][k+1] + hmm->tsc[CTMD][k]);
	  
	  bdmx[i][0] = -INFTY; /* D_0 does not exist */
	  
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
      gamma[i+1]  = gamma[i+2] + 0; /* extend without adding a new hit */
      /*printf("i: %d | gamma[i]: %f | gamma[i+1]: %f\n", i, gamma[i], gamma[i+1]);*/
      gback[i+1]  = -1;
      savesc[i+1] = IMPOSSIBLE;
      saver[i+1]  = -1;
      j = ((i+1+W-1) < L) ? (i+1+W-1) : L;
      sc = gamma[j+1] + Scorify(bmmx[i][0]) - min_thresh;

      /*printf("bmmx[%3d][0]: %10.5f\n", i, Scorify(bmmx[i][0]));*/
      bck_sum = ILogsum(bck_sum, bmmx[i][0]);

      bck_sc = Scorify(bmmx[i][0]);
      /*printf("bmmx[%3d][0]: %10.5f\n", i, bck_sc);*/
      /*printf("i: %d | bmmx[i][0]: %f | j: %d | gamma[i]: %f | sc: %f\n", i, Scorify(bmmx[i][0]), j, gamma[i], sc);*/
      if (sc > gamma[i+1])
	{
	  gamma[i+1]  = sc;
	  gback[i+1]  = j;
	  savesc[i+1] = Scorify(bmmx[i][0]);
	  saver[i+1]  = 0;
	  //saver[i+1]  = bestr[d];
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

  i     = 1;
  bck_nhits = 0;
  while (i <= L) {
    if (gback[i] == -1) /* no hit */
      i++; 
    else                /* a hit, a palpable hit */
      {
	bck_hitr[bck_nhits]   = saver[i];
	bck_hitj[bck_nhits]   = gback[i];
	bck_hiti[bck_nhits]   = i;
	bck_hitsc[bck_nhits]  = savesc[i];
	bck_nhits++;
	i = gback[i]+1;
	
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

  CP9_combine_FBscan_hits(L, W, fwd_nhits, fwd_hitr, fwd_hiti, fwd_hitj, fwd_hitsc, 
			        bck_nhits, bck_hitr, bck_hiti, bck_hitj, bck_hitsc,
			        ret_nhits, ret_hitr, ret_hiti, ret_hitj, ret_hitsc, 0);

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
debug_check_CP9_FBscan(struct cp9_dpmatrix_s *fmx, struct cp9_dpmatrix_s *bmx, 
		       struct cplan9_s *hmm, float sc, int L, unsigned char *dsq)
{
  int k, i;
  float max_diff;  /* maximum allowed difference between sc and 
		    * sum_k f[i][k] * b[i][k] for any i
		    */
  float diff;
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
	      printf("hmm->msc[dsq[i]][%d]: %d\n", k, hmm->msc[dsq[i]][k]);
	    }
	  fb_sum = ILogsum(fb_sum, (fmx->mmx[i][k] + bmx->mmx[i][k] - hmm->msc[dsq[i]][k]));
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  fb_sum = ILogsum(fb_sum, (fmx->imx[i][k] + bmx->imx[i][k] - hmm->isc[dsq[i]][k]));
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
CP9ScanFullPosterior(unsigned char *dsq, int L,
		     struct cplan9_s *hmm,
		     struct cp9_dpmatrix_s *fmx,
		     struct cp9_dpmatrix_s *bmx,
		     struct cp9_dpmatrix_s *mx)
{
  int i;
  int k;
  int fb_sum; /* tmp value, the probability that the current residue (i) was
	       * visited in any parse */
  float temp_sc;

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
      mx->imx[i][0] = fmx->imx[i][0] + bmx->imx[i][0] - hmm->isc[dsq[i]][0] - fb_sum;
      /*hmm->isc[dsq[i]][0] will have been counted in both fmx->imx and bmx->imx*/
      mx->dmx[i][0] = -INFTY; /*D_0 does not exist*/
      for (k = 1; k <= hmm->M; k++) 
	{
	  mx->mmx[i][k] = fmx->mmx[i][k] + bmx->mmx[i][k] - hmm->msc[dsq[i]][k] - fb_sum;
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  mx->imx[i][k] = fmx->imx[i][k] + bmx->imx[i][k] - hmm->isc[dsq[i]][k] - fb_sum;
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
 * Args:     L      - length of dsq
 *           W      - window length (maximum size of a hit considered)
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
CP9_combine_FBscan_hits(int L, int W, int fwd_nhits, int *fwd_hitr, int *fwd_hiti, 
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

  /*  for(i = 0; i < fwd_nhits; i++)
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
	  i = bck_hiti[bck_ctr];
	  j = bck_hitj[bck_ctr];
	  hiti[nhits]   = i;
	  hitj[nhits]   = j;
	  hitsc[nhits]  = bck_hitsc[bck_ctr];
	  hitr[nhits]   = bck_hitr[bck_ctr];
	  bck_ctr++;
	}
      else if(bck_ctr == bck_nhits) 
	{
	  i = fwd_hiti[fwd_ctr];
	  j = fwd_hitj[fwd_ctr];
	  hiti[nhits]   = i;
	  hitj[nhits]   = j;
	  hitsc[nhits]  = fwd_hitsc[fwd_ctr];
	  hitr[nhits]   = fwd_hitr[fwd_ctr];
	  fwd_ctr--;
	}
      else
	{
	  i = bck_hiti[bck_ctr];
	  j = fwd_hitj[fwd_ctr];
	  if((j-i+1 > 0) && (j-i+1 <= W))
	    {
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
	      hitj[nhits]   = bck_hitj[bck_ctr]; /* this will be i+W */
	      hitsc[nhits]  = bck_hitsc[bck_ctr];
	      hitr[nhits]   = bck_hitr[bck_ctr];
	      bck_ctr++;
	    }
	  else if(i > j) /* hit from Backward() unmatched by a Forward hit */
	    {
	      hitj[nhits]   = j;
	      hiti[nhits]   = fwd_hiti[fwd_ctr]; /* this will be i+W */
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



/* bandcyk.c
 * SRE, Wed Nov 20 07:46:56 2002 [flight home from Airlie mtg]
 * CVS $Id$
 * 
 * Banded CYK implementation.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "squid.h"
#include "vectorops.h"
#include "structs.h"
#include "funcs.h"


/* Function: BandDistribution()
 * Date:     SRE, Wed Nov 20 07:15:31 2002 [flight home from Washington DC]
 *
 * Purpose:  Given a model, and a DP maximum window size W; calculate
 *           cumulative distribution g_v(n) = P(len <= n | subgraph rooted at v).
 *           Return as matrix g[v][n].
 *
 * Args:     cm   - CM to calculate length distribution for.
 *           W    - maximum DP window size.       
 *
 * Returns:  gamma[v][n] (0..M-1, 0..W).
 *           Caller free's w/ DMX2Free(gamma).
 */
double **
BandDistribution(CM_t *cm, int W)
{
  double **gamma;            /* gamma[v][n] = log P(length n | state v); [0..W][0..M-1] */
  int      n,x;
  int      v,y;

  n     = MAX(MAXCONNECT, W+1);
  gamma = DMX2Alloc(cm->M, W+1);
  
  for (n = 0; n <= W; n++)
    for (v = cm->M-1; v >= 0; v--)
      {
	gamma[v][n] = 0.;

	switch (cm->sttype[v]) {
	case S_st:
	case D_st:
	  for (y = 0; y < cm->cnum[v]; y++)
	    gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n];
	  break;

	case ML_st:
	case MR_st:
	case IL_st:
	case IR_st:
	  if (n >= 1) {
	    for (y = 0; y < cm->cnum[v]; y++)
	      gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-1];
	  }
	  break;

	case MP_st:
	  if (n >= 2) {
	    for (y = 0; y < cm->cnum[v]; y++)
	      gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-2];
	  } 
	  break;
  
	case B_st:
	  for (x = 0; x <= n; x++) /* x = length of left child */
	    gamma[v][n] += gamma[cm->cfirst[v]][x] * gamma[cm->cnum[v]][n-x];
	  break;

	case E_st:
	  if (n == 0) gamma[v][n] = 1.;
	  break;
	
	default: 
	  Die("gamma on fire");
	}
      }

  /* Reduce numerical imprecision issues: renormalize the distributions.
   * Don't do this if you're debugging; it makes all other problems go away too.
   */
  for (v = 0; v < cm->M; v++)
    DNorm(gamma[v], W+1);

  /* Convert to cumulative distribution, P(len <= n) 
   * Squash any numerical imprecision that overflows past 1.0.
   */
  for (v = 0; v < cm->M; v++)
    for (n = 1; n <= W; n++)
      {
	gamma[v][n] = gamma[v][n] + gamma[v][n-1];
	if (gamma[v][n] > 1.) gamma[v][n] = 1.;
      }
				
  return gamma;
}

/* Function: BandBounds()
 * Date:     SRE, Fri Sep 26 10:01:41 2003 [AA 886 from Salk]
 *
 * Purpose: Given a cumulative probability distribution gamma[n] =
 *          P(length <= n) for each state v and a tail
 *          cutoff probability p.
 * 
 *          Find dmin and dmax for each v, such that the probability
 *          of missing a hit is <= p on either the low side or high side:
 *             gamma[dmin-1]   <= p
 *             1 - gamma[dmax] <= p
 * 
 *          The total probability mass these bounds encompass is:
 *             gamma[dmax] - gamma[dmin-1] >= 1-2p   
 *
 * Args:    gamma    - cumulative probability distribution P(length <= n) for state v;
 *                     [0..v..M-1][0..W] 
 *          M        - # of states in CM
 *          W        - maximum subsequence length in DP         
 *          p        - tail cutoff probability
 *          ret_dmin - RETURN: minimum d bound for each state v; [0..v..M-1]
 *          ret_dmax - RETURN: maximum d bound for each state v; [0..v..M-1]
 *
 * Returns: (void) 
 *          dmin, dmax are allocated here. Caller must free.
 */
void
BandBounds(double **gamma, int M, int W, double p, int **ret_dmin, int **ret_dmax)
{
  int     *dmin, *dmax;
  int      v;

#if (SRE_CONLEVEL >= 1)
  for (v = 0; v < M; v++) 
    {
      assert(gamma[v][W] <= 1.0 && gamma[v][0] >= 0.0);
      assert(gamma[v][W] >= gamma[v][0]);
      assert(p >= 0. && p <= 0.5);
    }
#endif

  dmin = MallocOrDie(sizeof(int) * M);
  dmax = MallocOrDie(sizeof(int) * M);

  for (v = 0; v < M; v++)
    {
      dmin[v] = 0; while (dmin[v] <= W && gamma[v][dmin[v]]   <= p)    dmin[v]++;
      dmax[v] = W; while (dmax[v] > 0  && gamma[v][dmax[v]-1] >= 1.-p) dmax[v]--;
    }

#if (SRE_CONLEVEL >= 1)
  for (v = 0; v < M; v++)
    {
      if (dmin[v] > 0) 
	assert(gamma[v][dmax[v]] - gamma[v][dmin[v]-1] >= 1.-2*p);
      else
	assert(gamma[v][dmax[v]] >= 1.-2*p);
      assert(dmin[v] >= 0 && dmin[v] <= W);
      assert(dmax[v] >= 0 && dmax[v] <= W);
      assert(dmax[v] >= dmin[v]);
    }
#endif  
  
  *ret_dmin = dmin;
  *ret_dmax = dmax;
}
      

/* A couple of quick hacks. ...
 * Print an XMGRACE xy file for a specified v, showing the
 * cumulative distribution. Needed this for the R01 renewal.
 * SRE, Wed Feb 19 08:35:32 2003
 */
void
PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int W)
{
  int n;

  for (n = 0; n <= W; n++)
    fprintf(fp, "%d %.6f\n", n, gamma[v][n]);
  fprintf(fp, "&\n");
  fprintf(fp, "%d  0\n",   min[v]);
  fprintf(fp, "%d  1.0\n", min[v]);
  fprintf(fp, "&\n");
  fprintf(fp, "%d  0\n",   max[v]);
  fprintf(fp, "%d  1.0\n", max[v]);
  fprintf(fp, "&\n");
}
/* ... and, estimate the total savings in DP cells filled.
 */
void
PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W)
{
  int v;
  int after, before;

  before = after = 0;
  for (v = 0; v <= cm->M; v++) 
    {
      if (cm->sttype[v] != E_st) {
	after  += max[v] - min[v] + 1;
	before += W;
      }
    }
  printf("Before:  something like %d\n", before);
  printf("After:   something like %d\n", after);
  printf("Speedup: maybe %.2f fold\n", (float) before / (float) after);
}


/* Function: CYKBandedScan()
 * Date:     SRE, Fri Sep 26 09:53:42 2003 [AA 886, returning from Salk Institute]
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm.
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *
 *           dmin,dmax set the bounds. Both are arrays, 0..v..cm->M. 
 *           The band for v is dmin[v]..dmax[v], inclusive; that is,
 *           dmin[v] is the minimum allowed d for state v (inclusive);
 *           dmax[v] is the maximum; that is, v can only account for
 *           subsequences of length >= dmin[v] and <= dmax[v].
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; 1..L
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           L         - length of sequence
 *           W         - max d: max size of a hit
 *           ret_nhits - RETURN: number of hits
 *           ret_hitr  - RETURN: start states of hits, 0..nhits-1
 *           ret_hiti  - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj  - RETURN: end positions of hits, 0..nhits-1
 *           ret_hitsc - RETURN: scores of hits, 0..nhits-1            
 *
 * Returns:  
 *           hiti, hitj, hitsc are allocated here; caller free's w/ free().
 */
void
CYKBandedScan(CM_t *cm, char *dsq, int *dmin, int *dmax, int L, int W, 
	      int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc)
{
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       nhits;		/* # of hits in optimal parse */
  int      *hitr;		/* initial state indices of hits in optimal parse */
  int      *hiti;               /* start positions of hits in optimal parse */
  int      *hitj;               /* end positions of hits in optimal parse */
  float    *hitsc;              /* scores of hits in optimal parse */
  int       alloc_nhits;	/* used to grow the hit arrays */

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
  alpha = MallocOrDie (sizeof(float **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	alpha[v] = MallocOrDie(sizeof(float *) * (W+1));
	for (j = 0; j <= W; j++)
	  alpha[v][j] = MallocOrDie(sizeof(float) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	alpha[v] = MallocOrDie(sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  alpha[v][j] = MallocOrDie(sizeof(float) * (W+1));
      }
  }
  bestr = MallocOrDie(sizeof(int) * (W+1));

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
	    if ((sc = alpha[y+yoffset][0][0] + cm->tsc[v][yoffset]) > alpha[v][0][0]) 
	      alpha[v][0][0] = sc;
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  if (alpha[v][0][0] < IMPOSSIBLE) alpha[v][0][0] = IMPOSSIBLE;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j < W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  /* Impose the bands.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  for (v = 0; v < cm->M; v++)
    {
      for (d = 0;         d < dmin[v]; d++) alpha[v][0][d] = alpha[v][1][d] = IMPOSSIBLE;
      for (d = dmax[v]+1; d <= W;      d++) alpha[v][0][d] = alpha[v][1][d] = IMPOSSIBLE;
    }

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+1));
  gamma[0] = 0;
  gback    = MallocOrDie(sizeof(int)   * (L+1));
  gback[0] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+1));
  saver    = MallocOrDie(sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = 1; j <= L; j++) 
    {
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = dmin[v]; d <= dmax[v] && d <= j; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->endsc[v]; 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[v][yoffset]) > alpha[v][jp][d]) 
		      alpha[v][jp][d] = sc;
		  if (alpha[v][jp][d] < IMPROBABLE) alpha[v][jp][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = dmin[v]; d <= dmax[v] && d <= j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v];
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = dmin[v]; d <= dmax[v] && d <= j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v];
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = dmin[v]; d <= dmax[v] && d <= j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v];
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  if (dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = dmin[v]; d <= dmax[v] && d <= j; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v];
		  for (k = 0; k <= d; k++) /* k is length of right fragment */
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      if ((sc = alpha[w][jp][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			alpha[v][cur][d] = sc;
		    }
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      for (d = dmin[0]; d <= dmax[0] && d <=j; d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
	      alpha[0][cur][d] = sc;

	  if (cm->flags & CM_LOCAL_BEGIN) {
	    for (y = 1; y < cm->M; y++) {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	    }
	  }
	  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	}

      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[j]  = gamma[j-1] + 0; /* extend without adding a new hit */
      gback[j]  = -1;
      savesc[j] = IMPOSSIBLE;
      saver[j]  = -1;
      for (d = dmin[0]; d <= dmax[0] && d <= j; d++) 
	{
	  i = j-d+1;
	  if (i == 0) sc = alpha[0][cur][d];
	  else        sc = gamma[i-1] + alpha[0][cur][d];
	  if (sc > gamma[j]) 
	    {
	      gamma[j]  = sc;
	      gback[j]  = i;
	      savesc[j] = alpha[0][cur][d]; 
	      saver[j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
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
  free(bestr);

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
  return;
}


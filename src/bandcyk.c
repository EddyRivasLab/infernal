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
#include <float.h>

#include "squid.h"
#include "vectorops.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

static int band_calculation(CM_t *cm, int W, double p_thresh, 
			    int **ret_dmin, int **ret_dmax);
static int  ok_truncation_error(double *density, int b, int W);
static void print_band_distribution(double *density, int W);

void
BandExperiment(CM_t *cm)
{
  int  W;
  int *dmin, *dmax;

  W = 1000;
  while (! band_calculation(cm, W, 0.00001, &dmin, &dmax))
    {
      W += 1000;
      SQD_DPRINTF1(("increasing W to %d, redoing band calculation...\n", W));
    }
}

/* Function:  
 * Incept:    SRE, Sat Oct 11 14:17:40 2003 [St. Louis]
 *
 * Purpose:   
 * 
 *            Let L_v(n) be the cumulative probability distribution,
 *            P(length <= n), for state v:
 *                L_v(n) = \sum_{i=0}^{n}  \gamma_v(i)
 *                
 *            For each state v, find a dmin such that the probability
 *            of missing a hit is <= p on the low side:
 *                L_v(dmin-1) <= p
 *                
 *            On the high side, let H_v(n) be 1-L_v(n): e.g.
 *            P(length > n) for state v. But it is important not 
 *            to calculate this as 1-L_v(n); because of numerical
 *            roundoff issues, it must be done as:
 *                H_v(n) = \sum_{i=n+1}{\infty} \gamma_v(i)
 *                
 *            Then for each state v, find a dmax such that the 
 *            probability of missing a hit is >= p:
 *               G_v(dmax) >= p
 *               
 *            But we can't calculate the sum to \infty; we have
 *            to truncate somewhere. If we know that our width W
 *            is wide enough that we can see a probability mass of
 *            p in the right tail, it is sufficient to be sure
 *            that \gamma[v][W] < DBL_EPSILON * p. (because
 *            1+x = 1 for x < DBL_EPSILON, p + xp = p.)

 *
 * Args:      
 *
 * Returns:   
 *
 * Xref: */
static int
band_calculation(CM_t *cm, int W, double p_thresh, int **ret_dmin, int **ret_dmax)
{
  double **gamma;               /* P(length = n) for each state v                 */
  double  *tmp;
  int     *dmin;                /* lower bound for band. */
  int     *dmax;                /* upper bound for band. */
  int      v;			/* counter over states, 0..M-1                    */
  int      y;			/* counter over connected states                  */
  int      n;			/* counter over lengths, 0..W */
  int      dv;			/* Delta for state v */
  int      leftn;		/* length of left subsequence under a bifurc      */
  double   pdf; 		/* P(<=n) or P(>=n) for this state v */
  int     *touch;               /* touch[y] = # of higher states that depend on y */
  Mstack_t *beamstack;          /* pool of beams we can reuse  */
  int      status;		/* return status. */

  /* gamma[v][n] is Prob(state v generates subseq of length n)
   */
  gamma = MallocOrDie(sizeof(double *) * cm->M);        
  for (v = 0; v < cm->M; v++) gamma[v] = NULL;

  /* dmin[v] and dmax[v] are the determined bounds that we return.
   */
  dmin = MallocOrDie(sizeof(int) * cm->M);
  dmax = MallocOrDie(sizeof(int) * cm->M);  
  *ret_dmin = NULL;
  *ret_dmax = NULL;

  /* beamstack is a trick for reusing memory: a pushdown stack of 
   * "beams" (gamma[v] rows) we can reuse.
   */
  beamstack = CreateMstack();

  /* The second component of memory saving is the "touch" array.
   * touch[y] is the number of states above state [y] that will
   * depend on y but haven't been calculated yet. When we're done
   * calculating a new state v, we decrement touch[y] for all
   * y \in C_v. Any time touch[y] reaches 0, we put that beam
   * back into the pool for reuse.
   */
  touch  = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++) touch[v] = cm->pnum[v];

  /* Allocate and initialize the shared end beam.
   */
  gamma[cm->M-1] = MallocOrDie(sizeof(double) * (W+1));
  DSet(gamma[cm->M-1], W+1, 0.);
  gamma[cm->M-1][0] = 1.0;

  for (v = cm->M-1; v >= 0; v--)
    {
      /* Get a beam of memory from somewhere.
       *   1. If we're an E_st, we're sharing the end beam, and
       *      it's already initialized for us; don't do anything
       *      else to it.
       *   2. If there's a beam in the pool we can reuse, take it
       *      and set it back to 0's.
       *   3. Else, allocate and initialize to 0's.
       */
      if (cm->sttype[v] == E_st) {
	gamma[v] = gamma[cm->M-1];
	continue;
      }

      if ((gamma[v] = PopMstack(beamstack)) == NULL)
	gamma[v] = MallocOrDie(sizeof(double) * (W+1));		    
      DSet(gamma[v], W+1, 0.);


      /* Recursively calculate the probability density P(length=n) for this state v.
       * (The heart of the algorithm is right here.)
       */
      if (cm->sttype[v] == B_st) 
	{
	  pdf = 0.;
	  for (n = 0; n <= W; n++)
	    {
	      for (leftn = 0; leftn <= n; leftn++) 
		gamma[v][n] += gamma[cm->cfirst[v]][leftn]*gamma[cm->cnum[v]][n-leftn];
	      pdf += gamma[v][n];
	    }
	}
      else 
	{
	  pdf = 0.;
	  dv = StateDelta(cm->sttype[v]);
	  for (n = dv; n <= W; n++)
	    {
	      for (y = 0; y < cm->cnum[v]; y++)
		gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-dv];
	      pdf += gamma[v][n];
	    }
	}

      /* Make sure we've captured "enough" of the distribution (e.g.,
       * we have captured the right tail; our truncation error is
       * negligible).
       *   Of our 3 criteria, we apply two to every state:
       *     1. we're on the right side of the density (pdf is > 0.5
       *        would be enough, but we use .999)
       *     2. gamma_v(W) < p * DBL_EPSILON
       *        This must be true if \sum_{i=W+1...\infty} g(i) < p*DBL_EPSILON,
       *        which is really what we're trying to prove
       */
      if (pdf < 0.999 || gamma[v][W-1000] > p_thresh * DBL_EPSILON)
	{
	  status = 0;  	/* fail; truncation error unacceptable; caller must increase W and rerun. */
	  goto CLEANUP;
	}
      
      /* Renormalize this beam. (Should we really be doing this?)
       */
      if (pdf > 1.0) DNorm(gamma[v], W+1);

      /* Determine our left bound, dmin.
       */
      pdf = 0.;
      for (n = 0; n <= W; n++)
	{
	  pdf += gamma[v][n];
	  if (pdf > p_thresh) { dmin[v] = n; break; }
	}

      /* And our right bound, dmax.
       */
      pdf = 0.;
      for (n = W; n >= 0; n--)
	{
	  pdf += gamma[v][n];
	  if (pdf > p_thresh) { dmax[v] = n; break; }
	}

      /* Reuse memory where possible, using the "touch" trick:
       *   look at all children y \in C_v.
       *   decrement touch[y]
       *   if touch[y] reaches 0, no higher state v depends on this
       *     state's numbers; release the memory.
       *   we're reusing the end state for every E, so don't free it
       *   'til we're done.
       */
      if (cm->sttype[v] == B_st)
	{  /* the connected children of a B st are handled specially, remember */
	  y = cm->cfirst[v];  PushMstack(beamstack, gamma[y]); gamma[y] = NULL;
	  y = cm->cnum[v];    PushMstack(beamstack, gamma[y]); gamma[y] = NULL;
	}
      else
	{
	  for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	    {
	      touch[y]--;
	      if (touch[y] == 0 && cm->sttype[y] != E_st) {
		PushMstack(beamstack, gamma[y]); 
		gamma[y] = NULL;
	      }
	    }
	}
    }


  print_band_distribution(gamma[0],W);
  printf ("%d  %d\n", dmin[0], -800);
  printf ("%d  %d\n", dmin[0], 0);
  printf ("&\n");
  printf ("%d  %d\n", dmax[0], -800);
  printf ("%d  %d\n", dmax[0], 0);
  printf ("&\n");
  printf ("%d  %d\n", W-1000, -800);
  printf ("%d  %d\n", W-1000, 0);
  printf ("&\n");

  if (! ok_truncation_error(gamma[0], dmax[0], W-1000)) { status = 0; goto CLEANUP; }

  *ret_dmin = dmin;
  *ret_dmax = dmax;
  status = 1;

 CLEANUP:
  free(touch);

  /* Free whatever's left in the gamma beams.
   */
  for (v = 0; v < cm->M; v++) 
    if (cm->sttype[v] != E_st && gamma[v] != NULL) 
      { free(gamma[v]); gamma[v] = NULL; }
  free(gamma[cm->M-1]);		/* free the end state */
  free(gamma);
  
  /* Free the reused stack of beams.
   */
  while ((tmp = PopMstack(beamstack)) != NULL) free(tmp);
  FreeMstack(beamstack);

  return status;
}


static int
ok_truncation_error(double *density, int b, int W)
{
  double logbeta;
  double beta;			/* geometric decay parameter */
  double C;			/* area under density from b+1..W inclusive */
  double D;			/* area under unseen density from W+1..\infty */
  int    i;
  
  C = 0.;
  for (i = b+1; i <= W; i++) C += density[i];

  /* if density is falling off as a geometric, log(beta) is 
   * the slope of log(p). Estimate slope by a simple 2-point
   * fit at our boundaries.
   */
  logbeta = (log(density[W]) - log(density[b+1])) / (W - b - 1);
  beta = exp(logbeta);
	     
  D = (beta / (1.-beta)) * density[W];

  printf("%d  %g\n", b+1, log(density[b+1]));
  printf("%d  %g\n", W, log(density[b+1]) + (W-b-1) * logbeta);
  printf("%d  %g\n", W+1000, log(density[b+1]) + (W-b-1+1000) * logbeta);
  printf("&\n");

  if (D < C * DBL_EPSILON) return 1;
  else                     return 0;
}  
  


static void
print_band_distribution(double *density, int W)
{
  int    i;
  double pdf;

  /* Set 0. The density.
   */
  for (i = 0; i <= W; i++)
    printf("%d\t%g\n", i, density[i]);
  printf("&\n");

  /* Set 1. The cumulative pdf.
   */
  pdf = 0.;
  for (i = 0; i <= W; i++) {
    pdf += density[i];
    printf("%d\t%g\n", i, pdf);
  }
  printf("&\n");

  /* Set 2. log(density). For fitting & determining geometric 
   *        decay parameter.
   */
  for (i = 0; i <= W; i++) 
    printf("%d\t%g\n", i, log(density[i]));
  printf("&\n");
}


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


/* Function: BandedParsetreeDump()
 * Date:     SRE, Sat Sep 27 15:44:22 2003 [St. Louis]
 * Xref:     STL7 pg. 119
 *
 * Purpose:  Generate a detailed picture of a parsetree data structure,
 *           annotated with relevant information from the sequence
 *           and the model -- and show agreement (or disagreement)
 *           with calculated bands.
 *           
 *           Modified from ParsetreeDump().
 *
 * Args:    fp    - FILE to write output to.
 *          tr    - parsetree to examine.
 *          cm    - model that was aligned to dsq to generate the parsetree
 *          dsq   - sequence that was aligned to cm to generate the parsetree
 *          gamma - cumulative subsequence length probability distributions
 *                  used to generate the bands; from BandDistribution(); [0..v..M-1][0..W]
 *          W     - maximum window length W (gamma distributions range up to this)        
 *          dmin  - minimum subseq length for each state; from BandBounds(); [0..v..M-1]
 *          dmax  - maximum ""
 *
 * Returns:  (void)
 */
void
BandedParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, char *dsq, 
		    double **gamma, int W, int *dmin, int *dmax)
{
  int   x;
  char  syml, symr;
  float tsc;
  float esc;
  int   v,y;
  int   L;

  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	  " idx ", "emitl", "emitr", "state", " nxtl", " nxtr", " prv ", " tsc ", " esc ", 
	  " L   ", " dmin", " dmax");
  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	  "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----",
	  "-----", "-----", "-----");
  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];

      /* Set syml, symr: one char representation of what we emit, or ' '.
       * Set esc:        emission score, or 0.
       * Only P, L, R states have emissions.
       */
      syml = symr = ' ';
      esc = 0.;
      if (cm->sttype[v] == MP_st) {
	syml = Alphabet[(int)dsq[tr->emitl[x]]]; 
	symr = Alphabet[(int)dsq[tr->emitr[x]]];
	esc  = DegeneratePairScore(cm->esc[v], dsq[tr->emitl[x]], dsq[tr->emitr[x]]);
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	syml = Alphabet[(int)dsq[tr->emitl[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitl[x]]);
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	symr = Alphabet[(int)dsq[tr->emitr[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitr[x]]);
      }

      /* Set tsc: transition score, or 0.
       * B, E, and the special EL state (M, local end) have no transitions.
       */
      tsc = 0.;
      if (v != cm->M && cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
	y = tr->state[tr->nxtl[x]];

	if (v == 0 && (cm->flags & CM_LOCAL_BEGIN))
	  tsc = cm->beginsc[y];
	else if (y == cm->M) /* CM_LOCAL_END is presumably set, else this wouldn't happen */
	  tsc = cm->endsc[v];
	else 		/* y - cm->first[v] gives us the offset in the transition vector */
	  tsc = cm->tsc[v][y - cm->cfirst[v]];
      }

      /* Print the info line for this state
       */
      L = tr->emitr[x]-tr->emitl[x]+1;
      fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5d %5d %5d %5.2f %5.2f %5d %5d %5d %2s\n",
	      x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], 
	      Statetype(cm->sttype[v]), tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc,
	      L, dmin[v], dmax[v],
	      (L >= dmin[v] && L <= dmax[v]) ? "" : "!!");
    }
  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	  "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----",
	  "-----", "-----", "-----");
  fflush(fp);
} 

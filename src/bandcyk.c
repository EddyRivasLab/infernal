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

#include <math.h>
#include "squid.h"
#include "vectorops.h"
#include "structs.h"

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
 *           Caller free's w/ DMX2Free(mx).
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
   */
  for (v = 0; v < cm->M; v++)
    for (n = 1; n <= W; n++)
      gamma[v][n] = gamma[v][n] + gamma[v][n-1];
				
  return gamma;
}


void
BandBounds(double **gamma, int M, int W, double p, int **ret_min, int **ret_max)
{
  int     *min, *max;
  double   mincut, maxcut;
  int      v;

  min = MallocOrDie(sizeof(int) * (M+1));
  max = MallocOrDie(sizeof(int) * (M+1));

  mincut = p;
  maxcut = 1-p;

  for (v = 0; v < M; v++)
    {
      min[v] = 0; while (gamma[v][min[v]]   <= mincut) min[v]++;
      max[v] = W; while (gamma[v][max[v]-1] >= maxcut) max[v]--;
    }
  *ret_min = min;
  *ret_max = max;
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

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
  
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == E_st) gamma[v][0] = 1.0;
    else gamma[v][0] = 0.0;

  for (n = 1; n <= W; n++)
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
	  for (y = 0; y < cm->cnum[v]; y++)
	    gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-1];
	  break;

	case MP_st:
	  if (n > 1) {
	    for (y = 0; y < cm->cnum[v]; y++)
	      gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-2];
	  } 
	  break;
  
	case B_st:
	  for (x = 0; x <= n; x++) /* x = length of left child */
	    gamma[v][n] += gamma[cm->cfirst[v]][x] * gamma[cm->cnum[v]][n-x];
	  break;

	case E_st:
	  break;
	
	default: 
	  Die("gamma on fire");
	}
      }

  /* Convert to cumulative distribution, P(len <= n) 
   */
  for (v = 0; v < cm->M; v++)
    for (n = 1; n <= W; n++)
      gamma[v][n] = gamma[v][n] + gamma[v][n-1];
				
  return gamma;
}


void
BandBounds(double **gamma, int M, int W, double p)
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
      min[v] = 0; while (gamma[v][min[v]]   <= mincut)   min[v]++;
      max[v] = W; while (gamma[v][max[v]-1] >= maxcut) max[v]--;
    }

  free(min);
  free(max);
}
      

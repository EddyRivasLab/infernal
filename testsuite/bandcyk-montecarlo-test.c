/* bandcyk-montecarlo-test.c
 * SRE, Tue Oct 21 17:00:44 2003
 * From notebook/1015-infernal-bands; xref STL7 p.130
 * Easelfied: EPN, Fri Nov 30 09:48:32 2007
 *
 * Produce densities both by Monte Carlo and by calculation engine.
 * Compare each density distribution by chi-square test.
 *
 * Small numbers (close to zero) indicate that the two densities
 * are significantly different.
 * 
 * If any chi-square is less than a threshold (0.01), fail.
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 * SVN $Id$
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stats.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                             docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",          0 },
  { "-n",        eslARG_INT,  "10000", NULL, NULL,  NULL,  NULL, NULL, "number of monte carlo samples to do",           0 },
  { "-s",        eslARG_INT,     NULL, NULL, "n>0", NULL,  NULL, NULL, "set random number seed for Monte Carlo to <n>", 0 },
  { "-t",        eslARG_REAL,   "0.01",NULL, "x>0.",NULL,  NULL, NULL, "threshold for rejecting hypothesis that distros are identical ", 0 },
  { "-W",        eslARG_INT,   "1000", NULL, "n>0", NULL,  NULL, NULL, "set maximum W (subseq length) to <n>",          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "Monte Carlo sampling test program for query-dependent band calculation";

static double DChiSquareFit(double *d1, double *d2, int N);

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go    = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char    *cmfile = esl_opt_GetArg(go, 1);
  CMFILE  *cmfp;		/* open CM file for reading */
  CM_t    *cm;			/* a covariance model       */
  int      v;			/* counter over states */
  
  double **mc_gamma;		/* densities from monte carlo   */
  double **gamma;		/* densities from calculation engine  */
  double   p;			/* p from chi squared test */
  double   threshold;		/* probability threshold for rejecting */

  int      maxW;		/* maximum length that densities are calc'ed for */
  int      mc_nsample;		/* # of monte carlo samples to do */
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;

  /*********************************************** 
   * Parse command line
   ***********************************************/

  maxW           = esl_opt_GetInteger(go, "-W");
  mc_nsample     = esl_opt_GetInteger(go, "-n");
  threshold      = esl_opt_GetReal   (go, "-t");

  /* create RNG */
  if (! esl_opt_IsDefault(go, "-s")) 
    r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
  else r = esl_randomness_CreateTimeseeded();

  /*********************************************** 
   * Preliminaries: get our CM
   ***********************************************/

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  /*****************************************************************
   * Do the band calculations
   ****************************************************************/

  /* BandMonteCarlo() collects "density" as unnormalized counts
   */
  if (! BandMonteCarlo(cm, mc_nsample, maxW, &mc_gamma))
    cm_Fail("Your maxW (%d) must be too small, sorry...\n", maxW);

  /* BandCalculationEngine() calculates a real density for each state v
   */
  if (! BandCalculationEngine(cm, maxW, 0.001, TRUE, NULL, NULL, &gamma, NULL))
    cm_Fail("Your maxW (%d) must be too small, sorry...\n", maxW);

  for (v = 0; v < cm->M; v++)
    {
      esl_vec_DScale(gamma[v],    maxW+1, esl_vec_DSum(mc_gamma[v], maxW+1)); /* convert to #'s */
      p = DChiSquareFit(gamma[v], mc_gamma[v], maxW+1);	      /* compare #'s    */

      if (cm->sttype[v] != E_st 
	  && cm->ndtype[cm->ndidx[v]+1] != END_nd /* skip nodes with unreachable inserts */
	  && p < threshold)
	cm_Fail("Rejected band distribution for state %d: chi-squared p = %f\n", v, p);
    }
  FreeBandDensities(cm, mc_gamma);
  FreeBandDensities(cm, gamma);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);

  exit(0);
}

static double 
DChiSquareFit(double *d1, double *d2, int N)
{
  int    i;
  double diff;
  double chisq = 0.0;
  int    n;
  double qax;
  
  n = 0;
  for (i = 0; i < N; i++)
    {
      if (d1[i] == 0. && d2[i] == 0.) continue;
      diff = d1[i] - d2[i];
      chisq += diff * diff / (d1[i]+d2[i]);
      n++;
    }

  if (n > 1) {
    if(esl_stats_IncompleteGamma(((double) n-1.)/2., chisq/2., NULL, &qax) != eslOK)
      cm_Fail("DChiSquareFit() call to esl_stats_IncompleteGamma() failed.");
    return qax;
  }
  else return -1.;
}

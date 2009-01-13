/* bandcyk-truncation-test.c
 * SRE, Tue Oct 21 17:07:56 2003 [St. Louis]
 * Adapted from check-bandtruncation.c; xref notebook/1015-infernal-bands, STL7 p. 130.
 * Easelfied: EPN, Fri Nov 30 10:04:16 2007
 * 
 * Implements two tests of the truncation error:
 *   1. Calculate bands for two different sizes of W; verify that bands
 *      are identical.
 *      
 *   2. Verifies that the geometric decay assumption used to estimate
 *      the truncation error was indeed an upper bound, by comparing
 *      the predicted density from a shorter length W1 to the observed
 *      density for a longer length W2, in the interval n=W1+1..W2.
 * 
 * Band calculation is supposed to be guaranteed to be independent
 * of input W -- the BandCalculationEngine() returns an error whenever
 * the input W is too short to guarantee this. If truncation error calculation 
 * in BandCalculationEngine() is wrong, bands may differ.
 * 
 * Example tests;
 * using models built w/ --rf from alignments in intro/ subdirectory:
 * ./bandcyk-truncation-test trna.cm 160 1000
 * ./bandcyk-truncation-test rp.cm   525 2000
 * ./bandcyk-truncation-test ssu.cm 1687 3000
 * 
 * xref STL7 p.130.
 ******************************************************************
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

#include "easel.h"
#include "esl_getopts.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name        type         default  env  range toggles reqs incomp  help                                            docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",                0 },
  { "-p",        eslARG_REAL,  "1E-5", NULL, "x>0.",NULL,  NULL, NULL, "set tail probability threshold to <p>",               0 },
  { "--verbose", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show output (default: silently return 0 on success)", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <W1> <W2> where width W1 < W2.";
static char banner[] = "Check of truncation error calculation in BandCalculationEngine()";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 3, argc, argv, banner, usage);

  char    *cmfile = esl_opt_GetArg(go, 1);
  int      W1     = atoi(esl_opt_GetArg(go, 2));
  int      W2     = atoi(esl_opt_GetArg(go, 3));
  CMFILE  *cmfp;		/* open CM file for reading */
  CM_t    *cm;			/* a covariance model       */
  int      v;			/* counter over states */
  int      n;			/* counter over lengths 0..W */

  double   p_thresh;            /* tail probability threshold for banding */
  int     *dmin1, *dmax1;       /* calculated bands from pass 1 */
  int     *dmin2, *dmax2;       /* calculated bands from pass 2 */
  double **gamma1, **gamma2;    /* probability densities from pass 1, 2 */
  double   beta1;		/* geometric decay constant, truncation error calc */
  double   g;			/* log of an estimated gamma[n] */
  int      be_verbose;
  ESL_ALPHABET   *abc     = NULL;

  /*********************************************** 
   * Parse command line
   ***********************************************/

  p_thresh   = esl_opt_GetReal   (go, "-p");
  be_verbose = esl_opt_GetBoolean(go, "--verbose");

  if (W1 >= W2) cm_Fail("Please set a width W1 < width W2, else the check won't work right.");

  /* Get our CM
   */

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if ((CMFileRead(cmfp, NULL, &abc, &cm)) != eslOK) cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  /* Do two band calculations with the different W's.
   * Save the gamma_0 densities for root state 0.
   */
  if (! BandCalculationEngine(cm, W1, p_thresh, FALSE, &dmin1, &dmax1, &gamma1, NULL))
    cm_Fail("Your W1 (%d) must be too small, sorry.\n", W1);
  if (! BandCalculationEngine(cm, W2, p_thresh, FALSE, &dmin2, &dmax2, &gamma2, NULL))
    cm_Fail("Your W2 (%d) must be too small, sorry.\n", W2);

  /* Verify that the bands are all identical, regardless of choice of W.
   */
  for (v = 0; v < cm->M; v++)
    if (dmin1[v] != dmin2[v] || dmax1[v] != dmax2[v])
      cm_Fail("failed at v=%d: Band for W1=%d is %d..%d; for W2=%d is %d..%d\n",
	  v, W1, dmin1[v], dmax1[v], W2, dmin2[v], dmax2[v]);
      
  /* Verify that the geometric decay predicted for the shorter width W1
   * was indeed an upper bound for what we saw in the longer width W2,
   * in the interval W1+1..W2.
   */
  if (! BandTruncationNegligible(gamma1[0], dmax1[0], W1, &beta1))
    cm_Fail("shouldn't happen, because we already checked this in the Engine().");
  for (n = W1+1; n <= W2; n++)
    {
      g = log(gamma1[0][W1]) + (n-W1) * log(beta1);   /* g = log(estimated gamma[n]); should be upper bound */
      if (g < log(gamma2[0][n])) 
	cm_Fail("truncation error test failed: geometric is not an upper bound on tail of gamma[n]");
    }
  
  /* Clean up and exit.
   */
  FreeBandDensities(cm, gamma1);
  FreeBandDensities(cm, gamma2);
  free(dmin1);
  free(dmax1);
  free(dmin2);
  free(dmax2);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  exit(0);
}

/* bandcyk-montecarlo-test.c
 * SRE, Tue Oct 21 17:00:44 2003
 * From notebook/1015-infernal-bands; xref STL7 p.130
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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "squid.h"
#include "sqfuncs.h"
#include "dirichlet.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

static char usage[] = "\
Usage: bandcyk-montecarlo-test [-options] <cmfile>\n\
  where options are:\n\
  -n <n> : number of monte carlo samples to do [default 10000]\n\
  -s <n> : set random number seed for Monte Carlo\n\
  -t <p> : threshold for rejecting hypothesis that distros are identical [default: 0.01]\n\
  -W <n> : set maximum W (subseq length) to <n> [default=1000]\n\
";

static struct opt_s OPTIONS[] = { 
  { "-n", TRUE, sqdARG_INT },
  { "-s", TRUE, sqdARG_INT },
  { "-t", TRUE, sqdARG_FLOAT },
  { "-W", TRUE, sqdARG_INT },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static double DChiSquareFit(double *d1, double *d2, int N);

int
main(int argc, char **argv)
{
  char    *cmfile;		/* file to read CM from */	
  CMFILE  *cmfp;		/* open CM file for reading */
  CM_t    *cm;			/* a covariance model       */
  int      v;			/* counter over states */
  
  double **mc_gamma;		/* densities from monte carlo   */
  double **gamma;		/* densities from calculation engine  */
  double   p;			/* p from chi squared test */
  double   threshold;		/* probability threshold for rejecting */

  int      maxW;		/* maximum length that densities are calc'ed for */
  int      mc_nsample;		/* # of monte carlo samples to do */
  int      seed;		/* random number seed for MC */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  maxW           = 1000;
  mc_nsample     = 10000;
  seed           = (int) time ((time_t *) NULL);
  threshold      = 0.01;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n") == 0) mc_nsample     = atoi(optarg);
    else if (strcmp(optname, "-s") == 0) seed           = atoi(optarg);
    else if (strcmp(optname, "-t") == 0) threshold      = atof(optarg);
    else if (strcmp(optname, "-W") == 0) maxW           = atoi(optarg);
  }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
 

  /*********************************************** 
   * Preliminaries: get our CM
   ***********************************************/

  sre_srandom(seed);

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);

  /*****************************************************************
   * Do the band calculations
   ****************************************************************/

  /* BandMonteCarlo() collects "density" as unnormalized counts
   */
  if (! BandMonteCarlo(cm, mc_nsample, maxW, &mc_gamma))
    Die("Your maxW (%d) must be too small, sorry...\n", maxW);

  /* BandCalculationEngine() calculates a real density for each state v
   */
  if (! BandCalculationEngine(cm, maxW, 0.001, TRUE, NULL, NULL, &gamma, FALSE))
    Die("Your maxW (%d) must be too small, sorry...\n", maxW);

  int n;
  for (v = 0; v < cm->M; v++)
    {
      DScale(gamma[v],    maxW+1, DSum(mc_gamma[v], maxW+1)); /* convert to #'s */
      p = DChiSquareFit(gamma[v], mc_gamma[v], maxW+1);	      /* compare #'s    */

      if (cm->sttype[v] != E_st 
	  && cm->ndtype[cm->ndidx[v]+1] != END_nd /* skip nodes with unreachable inserts */
	  && p < threshold)
	Die("Rejected band distribution for state %d: chi-squared p = %f\n", v, p);
    }

  FreeBandDensities(cm, mc_gamma);
  FreeBandDensities(cm, gamma);
  FreeCM(cm);
  SqdClean();
  exit(0);
}



static double 
DChiSquareFit(double *d1, double *d2, int N)
{
  int    i;
  double diff;
  double chisq = 0.0;
  int    n;
  
  n = 0;
  for (i = 0; i < N; i++)
    {
      if (d1[i] == 0. && d2[i] == 0.) continue;
      diff = d1[i] - d2[i];
      chisq += diff * diff / (d1[i]+d2[i]);
      n++;
    }

  if (n > 1) 
    return (IncompleteGamma(((double) n-1.)/2., chisq/2.));
  else 
    return -1.;
}

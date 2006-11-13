/* cp9-test.c
 * EPN, Mon Nov 13 17:49:07 2006
 * 
 * Test the CM -> CP9 HMM construction procedure.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "squid.h"
#include "sqfuncs.h"
#include "sre_stack.h"
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"
#include "funcs.h"

static char banner[] = "cp9-test - test CP9 HMM construction procedure";

static char usage[] = "\
Usage: cp9-test [-options] <cmfile>\n\
  where options are:\n\
  -h     : help; print brief help on version and usage\n\
  -s <n> : set random number seed\n\
  -t <p> : probability threshold for reporting violations [default: 1E-5]\n\
";

static char experts[] = "\
  --psionly   : only check that psi and phi values match\n\
  --nseq <n>  : use <n> samples to build CP9 HMM for [df: 50000]\n\
  --chi <f>   : fail sampling check if any chi-square test < <f> [df: 0.01]\n\
  --dlev <n>  : set verbosity of debugging print statements to <n> (1..3)\n\
";

static struct opt_s OPTIONS[] = { 
  { "-h", TRUE, sqdARG_NONE }, 
  { "-s", TRUE, sqdARG_INT },
  { "-t", TRUE, sqdARG_FLOAT },
  { "--psionly",   FALSE, sqdARG_NONE },
  { "--nseq",      FALSE, sqdARG_INT },
  { "--chi",       FALSE, sqdARG_FLOAT },
  { "--dlev",     FALSE, sqdARG_INT },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char              *cmfile;	  /* file to read CM from */	
  CMFILE            *cmfp;        /* open CM file for reading */
  CM_t              *cm;          /* a covariance model       */
  struct cplan9_s   *hmm;         /* constructed CP9 HMM; written to hmmfile              */
  CP9Map_t          *cp9map;      /* maps the hmm to the cm and vice versa */
  Stopwatch_t       *watch;       /* for timings */
  double             pthresh;     /* psi threshold for calling violations */
  int                seed;	  /* random number seed for MC */
  char              *optname;     /* name of option found by Getopt()        */
  char              *optarg;      /* argument found by Getopt()              */
  int                optind;      /* index in argv[]                         */
  int                nsamples;    /* Number of samples to build sampled HMM check. */
  float              chi_thresh;  /* if any chi-square test during the sampling check 
				   * is below this threshold, fail. */
  int                do_psionly;  /* don't do a sampling check only compare the expected
				   * number of times each HMM and CM state is entered */
  int                debug_level; /* verbosity of debugging printf statements */
  /*********************************************** 
   * Parse command line
   ***********************************************/
  pthresh        = 0.00001;
  do_psionly     = FALSE;
  nsamples       = 50000;
  chi_thresh     = 0.01;
  debug_level    = 0;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-s")          == 0) seed       = atoi(optarg);
    else if (strcmp(optname, "-t")          == 0) pthresh    = atof(optarg);
    else if (strcmp(optname, "--psionly")   == 0) do_psionly = TRUE;
    else if (strcmp(optname, "--nseq")      == 0) nsamples   = atoi(optarg);
    else if (strcmp(optname, "--chi")       == 0) chi_thresh = atof(optarg);
    else if (strcmp(optname, "--dlev")      == 0) debug_level= atoi(optarg);
    else if (strcmp(optname, "-h")          == 0) {
      MainBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
 
  if(nsamples < 10000)
    Die("Minimum number of samples allowed is 10,000.\n");

  /********************************************`*** 
   * Preliminaries: get our CM
   ***********************************************/

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);

  watch = StopwatchCreate(); 
  StopwatchZero(watch);
  StopwatchStart(watch);
  if(!build_cp9_hmm(cm, &hmm, &cp9map, debug_level))
    Die("CM Plan 9 HMM fails the psi/phi comparison test.\n");
  StopwatchStop(watch);
  StopwatchDisplay(stdout, "CP9 construction CPU time: ", watch);

  StopwatchZero(watch);
  StopwatchStart(watch);
  sre_srandom(seed);
  if(!(CP9_check_cp9_by_sampling(cm, hmm, 
				 NULL,     /* Don't keep track of failures (sub_cm feature) */
				 1, hmm->M, chi_thresh, nsamples, debug_level)))
    Die("CM Plan 9 fails sampling check!\n");
  else
    printf("CM Plan 9 passed sampling check.\n");

  StopwatchStop(watch);
  StopwatchDisplay(stdout, "CP9 sampling check CPU time: ", watch);

  /* clean up and exit */
  StopwatchFree(watch);
  FreeCP9Map(cp9map);
  FreeCPlan9(hmm);
  FreeCM(cm);
  return 0;
}

/* sub_cm-test.c
 * EPN, 09.18.06
 * 
 * Build many submodels from a template CM by choosing
 * random model start and end positions.
 * Compare the submodels to the corresponding stretch of
 * main model by determining the expected number of
 * times each state is entered. 
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
#include "dirichlet.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

static char banner[] = "sub_cm-test - test sub CM construction procedure";

static char usage[] = "\
Usage: sub_cm-test [-options] <cmfile>\n\
  where options are:\n\
  -h     : help; print brief help on version and usage\n\
  -n <n> : number of sub CMs to build and test [default 100]\n\
  -s <n> : set random number seed\n\
  -b <n> : set sub CM begin consensus (match) column as <n>\n\
  -e <n> : set sub CM end consensus (match) column as <n>\n\
  -t <p> : probability threshold for reporting violations [default: 1E-5]\n\
";

static char experts[] = "\
  --psionly   : only check that psi values match (don't build HMMs)\n\
  --sample    : build and check two CP9 HMMs (one an ML HMM via sampling)\n\
  --nseq <n>  : use <n> samples to build ML HMM for --samplecp9  [50000]\n\
  --chi <f>   : fail sampling check if any chi-square test < <f> [0.01]\n\
  --exhaust   : build and check every possible sub CM (all (N^2+N)/2)\n\
  --full      : build sub CM(s) with ONLY structure removed\n\
  --debug     : turn debugging print statements ON\n\
";

static struct opt_s OPTIONS[] = { 
  { "-h", TRUE, sqdARG_NONE }, 
  { "-n", TRUE, sqdARG_INT },
  { "-s", TRUE, sqdARG_INT },
  { "-t", TRUE, sqdARG_FLOAT },
  { "-b", TRUE, sqdARG_INT },
  { "-e", TRUE, sqdARG_INT },
  { "--psionly",   FALSE, sqdARG_NONE },
  { "--samplecp9", FALSE, sqdARG_NONE },
  { "--nseq",      FALSE, sqdARG_INT },
  { "--chi",       FALSE, sqdARG_FLOAT },
  { "--sample",    FALSE, sqdARG_NONE },
  { "--exhaust",   FALSE, sqdARG_NONE },
  { "--full",      FALSE, sqdARG_NONE },
  { "--debug",     FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *cmfile;		/* file to read CM from */	
  CMFILE  *cmfp;		/* open CM file for reading */
  CM_t    *cm;			/* a covariance model       */
  CM_t    *sub_cm;              /* sub covariance model     */
  int      v;			/* counter over states */
  int      nmodels;             /* number of sub CMs to build */
  int      spos;                /* start position for sub CM */
  int      epos;                /* end position for sub CM */
  int      ncols;               /* number of consensus (match) columns in CM */
  int      i;                   /* counter over sub CMs */
  int    **orig2sub_smap;       /* 2D state map from orig_cm (template) to sub_cm.
	                         * 1st dimension - state index in orig_cm 
				 * 2nd D - 2 elements for up to 2 matching sub_cm states */
  int    **sub2orig_smap;       /* 2D state map from orig_cm (template) to sub_cm.
				 * 1st dimension - state index in sub_cm (0..sub_cm->M-1)
				 * 2nd D - 2 elements for up to 2 matching orig_cm states */
  int      temp;

  double   threshold;		/* psi threshold for calling violations */
  int      seed;		/* random number seed for MC */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */
  int   begin_set;              /* TRUE if -b entered at command line */
  int   end_set;                /* TRUE if -e entered at command line */
  int  *imp_cc;                 /* imp_cc[k] = 1 if CP9 node k is an impossible case to get 
				 * the right transition distros for the sub_cm. */
  int do_atest;                 /* TRUE to build 2 ML HMMs, one from the CM and one from
				 * the sub_cm, analytically, and check to make sure
				 * the corresponding parameters of these two HMMS
				 * are within 'threshold' of each other.
				 */
  int do_stest;                 /* TRUE to build an ML HMM from a truncated MSA emitted from the
				 * original CM and test it via chi-squared tests against 
				 * the CP9 analytically built from the sub_cm.
				 */
  int do_exhaust;               /* TRUE to build every possible sub_cm */
  int do_fullsub;               /* TRUE to build sub CM(s) that model same number of columns
				 * as the template CM, with structure outside spos..epos
				 * removed.                          */
  int ndone;                    /* number of models built so far */
  int print_flag;               /* TRUE to print debug statements */
  int *awrong_predict_ct;       /* For 1 'atest's: the  number of times we predict we'll 
				 * fail the test for an HMM node and we get it right, for 
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  int *swrong_predict_ct;       /* For 1 sampling test: the  number of times we predict we'll 
				 * fail the test for an HMM node and we get it right, for 
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  int *awrong_predict_total_ct; /* For ALL 'atest's: the  number of times we predict we'll 
				 * fail the test for an HMM node and we get it right, for 
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  int *swrong_predict_total_ct; /* For ALL sampling tests: the  number of times we predict we'll 
				 * fail the test for an HMM node and we get it right, for 
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  int *apredict_ct;             /* For 1 'atest's: the  number of times we predict we'll 
				 * fail the test for an HMM node for
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  int *spredict_ct;             /* For 1 sampling test: the  number of times we predict we'll 
				 * fail the test for an HMM node for
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  int *apredict_total_ct;       /* For ALL 'atest's: the  number of times we predict we'll 
				 * fail the test for an HMM node for
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  int *spredict_total_ct;       /* For ALL sampling tests: the  number of times we predict we'll 
				 * fail the test for an HMM node for
				 * each of 5 cases - 5 different reasons we predict we'll fail.
				 */
  float chi_threshold;          /* if any chi-square test (which we haven't deemed 'impossible' 
				 * during the sampling check is below this threshold, fail. 
				 */
  int nsamples;                 /* Number of samples to build the ML HMM with during a sampling
				 * check.
				 */
  int npredict_cases;           /* Number of different cases for predicting a node's transitions
				 * will be impossible to match b/t the two HMMs. 
				 */
  
  /*********************************************** 
   * Parse command line
   ***********************************************/

  nmodels        = 100;
  seed           = (int) time ((time_t *) NULL);
  threshold      = 0.00001;
  begin_set      = FALSE;
  end_set        = FALSE;
  do_atest       = TRUE;
  do_stest       = FALSE;
  do_exhaust     = FALSE;
  do_fullsub     = FALSE;
  nsamples       = 50000;
  chi_threshold  = 0.01;
  npredict_cases = 6;
  print_flag     = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n") == 0) nmodels        = atoi(optarg);
    else if (strcmp(optname, "-s") == 0) seed           = atoi(optarg);
    else if (strcmp(optname, "-t") == 0) threshold      = atof(optarg);
    else if (strcmp(optname, "-b") == 0) { begin_set = TRUE; spos = atoi(optarg); nmodels = 1; }
    else if (strcmp(optname, "-e") == 0) { end_set   = TRUE; epos = atoi(optarg); }
    else if (strcmp(optname, "--psionly")   == 0) do_atest   = FALSE;
    else if (strcmp(optname, "--sample")    == 0) do_stest   = TRUE;
    else if (strcmp(optname, "--nseq")      == 0) nsamples = atoi(optarg);
    else if (strcmp(optname, "--chi")       == 0) chi_threshold = atof(optarg);
    else if (strcmp(optname, "--sample")    == 0) do_stest   = TRUE;
    else if (strcmp(optname, "--exhaust")   == 0) do_exhaust = TRUE;
    else if (strcmp(optname, "--full")      == 0) do_fullsub = TRUE;
    else if (strcmp(optname, "--debug")     == 0) print_flag = TRUE;
    else if (strcmp(optname, "-h") == 0) {
      MainBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
 
  if(begin_set && !end_set || !begin_set && end_set)
    Die("Must use both -b and -e or neither.\n");
  if(do_exhaust && begin_set)
    Die("--exhaust doesn't make sense with -b and -e.\n");
  if(do_exhaust && do_stest)
    Warn("--exhaust and --sample might take a long time...\n");
  if(begin_set && nmodels != 1)
    Die("-n does not make sense with -b and -e.\n");
  if(begin_set && spos > epos)
    Die("For -b <x> and -e <y> y must be >= x.\n");
  if(begin_set && spos > epos)
    Die("For -b <x> and -e <y> y must be >= x.\n");
  if(nsamples < 10000)
    Die("Minimum number of samples allowed is 10,000.\n");

  
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

  /* Allocate and initialize our *wrong_predict_total_ct arrays */
  apredict_total_ct       = MallocOrDie(sizeof(int) * (npredict_cases+1));
  spredict_total_ct       = MallocOrDie(sizeof(int) * (npredict_cases+1));
  awrong_predict_total_ct = MallocOrDie(sizeof(int) * (npredict_cases+1));
  swrong_predict_total_ct = MallocOrDie(sizeof(int) * (npredict_cases+1));
  for(i = 0; i <= npredict_cases; i++)
    {
      apredict_total_ct[i] = 0;
      spredict_total_ct[i] = 0;
      awrong_predict_total_ct[i] = 0;
      swrong_predict_total_ct[i] = 0;
    }

  /***********************************************************
   * Stategy: 
   * If do_atest (default), we do:
   *  1. Build a CP9 HMM (cp9_1) from the sub_cm (this is done).
   *  2. Build a CP9 HMM (cp9_2) from the full cm.
   *  3. Reconfig cp9_2 so start node is spos and end node is epos.
   *  4. Check corresponding parameters of cp9_1 and cp9_2 to make
   *     sure they're within the threshold.
   *
   * If do_stest, we also do:
   *  1. Build a CP9 HMM (cp9_1) from the sub_cm (this is done.)
   *  2. Sample a deep MSA from the CM. 
   *  3. Truncate the MSA between spos and epos.
   *  4. Perform chi-squared tests to see if sample from 
   *     (3) could have come from cp9_1. We do this by first 
   *     building a ML CP9 HMM (cp9_2) from the counts in 
   *     the truncated MSA.
   * 
   * No matter what, we always check to make sure that the
   * corresponding psi values of the orig_cm and sub_cm
   * are within 'threshold' of each other.
   *
   * This is all done within the build_sub_cm() function. 
   *********************************************************/

  /* Determine number of consensus columns modelled by CM */
  ncols = 0;
  for(v = 0; v <= cm->M; v++)
    {
      if(cm->stid[v] ==  MATP_MP)
	ncols += 2;
      else if(cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
	ncols++;
    }
  ndone = 0;
  if(do_exhaust) /* Build every possible sub CM. */
    {
      nmodels = (ncols * ncols + ncols) / 2;
      printf("Building and checking all possible sub CM (%5d different start/end positions):\n", ncols);

      for(spos = 1; spos <= ncols; spos++)
	{
	  printf("\tBuilding models with start pos: %5d (%5d / %5d completed)\n", spos, ndone, nmodels);
	  for(epos = spos; epos <= ncols; epos++)
	    {
	      if(!(build_sub_cm(cm, &sub_cm, spos, epos, NULL, NULL, &imp_cc, &apredict_ct, 
				&awrong_predict_ct, &spredict_ct, &swrong_predict_ct, threshold, do_fullsub, 
				do_atest, do_stest, chi_threshold, nsamples, print_flag)))
	      {
		printf("\nSub CM construction for spos: %4d epos: %4d failed one of the following tests:\n", spos, epos);
		printf("\tpsi test;            but this should never happen.\n");
		if(do_atest)
		  printf("\tanalytical HMM test; but this should never happen.\n");
		if(do_stest)
		  printf("\tsampling   HMM test; but this should never happen.\n");
		Die("\tLooks like there's a bug...\n");
	      }
	      /* keep track of number of each case of wrong prediction */
	      for(i = 1; i <= npredict_cases; i++)
		{
		  apredict_total_ct[i] += apredict_ct[i];
		  spredict_total_ct[i] += spredict_ct[i];
		  awrong_predict_total_ct[i] += awrong_predict_ct[i];
		  swrong_predict_total_ct[i] += swrong_predict_ct[i];
		}		  
	      FreeCM(sub_cm);
	      free(imp_cc);
	      free(apredict_ct);
	      free(spredict_ct);
	      free(awrong_predict_ct);
	      free(swrong_predict_ct);
	      ndone++;
	    }	      
	}
      printf("\nDone. %5d sub CM(s) were constructed and passed the following tests:\n", ndone);
      printf("\tpsi test\n");
      if(do_atest)
	printf("\tanalytical HMM test\n");
      if(do_stest)
	printf("\tsampling   HMM test\n");
    }	 
  else /* Build models with either preset begin point (spos) and end points (epos) 
	* or randomly chosen ones*/
    {
      if(begin_set && end_set)
	{
	  if(spos < 1) spos = 1;
	  if(epos > ncols) epos = ncols;
	  printf("Building a single sub CM with spos: %4d and epos: %4d ... ", spos, epos);
	}
      else
	printf("\tBuilding models with random start and end positions:\n");
      for(i = 0; i < nmodels; i++) /* if begin_set && end_set, nmodels is 1 */
	{
	  if(!(begin_set && end_set))
	    {
	      /* Randomly pick a start and end between 1 and ncols, inclusive */
	      spos = ((int) (sre_random() * ncols)) + 1;
	      epos = ((int) (sre_random() * ncols)) + 1;
	      if(spos > epos)
		{
		  temp = spos;
		  spos = epos;
		  epos = temp;
		}	      
	    }
	  if(!(build_sub_cm(cm, &sub_cm, spos, epos, NULL, NULL, &imp_cc, &apredict_ct, 
			    &awrong_predict_ct, &spredict_ct, &swrong_predict_ct, threshold, do_fullsub, 
			    do_atest, do_stest, chi_threshold, nsamples, print_flag)))
	    {
	      printf("\nDone. Sub CM construction for spos: %4d epos: %4d failed one of the following tests:\n", spos, epos);
	      printf("\tpsi test; but this should never happen.\n");
	      if(do_atest)
		printf("\tanalytical HMM test; but this should never happen.\n");
	      if(do_stest)
		printf("\tsampling   HMM test; but this should never happen.\n");
	      Die("\tLooks like there's a bug...\n");
	    }
	  /* keep track of number of each case of wrong prediction */
	  for(i = 1; i <= npredict_cases; i++)
	    {
	      apredict_total_ct[i] += apredict_ct[i];
	      spredict_total_ct[i] += spredict_ct[i];
	      awrong_predict_total_ct[i] += awrong_predict_ct[i];
	      swrong_predict_total_ct[i] += swrong_predict_ct[i];
	    }		  

	  FreeCM(sub_cm);
	  free(imp_cc);
	  free(apredict_ct);
	  free(spredict_ct);
	  free(awrong_predict_ct);
	  free(swrong_predict_ct);
	  ndone++;
	}
      printf("done.\n%5d sub CMs were constructed and passed the following tests:\n", ndone);
      printf("\tpsi test\n");
      if(do_atest)
	printf("\tanalytical HMM test\n");
      if(do_stest)
	printf("\tsampling   HMM test\n");
    }
  if(do_atest)
    {
      printf("\nPrinting summary of HMM nodes predicted to fail the analytical test:\n");
      for(i = 1; i <= npredict_cases; i++)
	printf("\tcase %d: %6d (%6d passed)\n", i, apredict_total_ct[i], awrong_predict_total_ct[i]);
    }
  if(do_stest)
    {
      printf("\nPrinting summary of HMM nodes predicted to fail the sampling test:\n");
      for(i = 0; i <= 6; i++)
	printf("\tcase %d: %6d (%6d passed)\n", i, spredict_total_ct[i], swrong_predict_total_ct[i]);
    }
  printf("\n");
  free(apredict_total_ct);
  free(spredict_total_ct);
  free(awrong_predict_total_ct);
  free(swrong_predict_total_ct);

  FreeCM(cm);
  return 0;
}

/* sub_cm-psi-test.c
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

static char usage[] = "\
Usage: sub_cm-psi-test [-options] <cmfile>\n\
  where options are:\n\
  -n <n> : number of sub CMs to build and test [default 100]\n\
  -s <n> : set random number seed\n\
  -t <p> : threshold for reporting violations in psi [default: 0.001]\n\
  -b <n> : set sub CM begin consensus (match) column as <n>\n\
  -e <n> : set sub CM end consensus (match) column as <n>\n\
";

static struct opt_s OPTIONS[] = { 
  { "-n", TRUE, sqdARG_INT },
  { "-s", TRUE, sqdARG_INT },
  { "-t", TRUE, sqdARG_FLOAT },
  { "-b", TRUE, sqdARG_INT },
  { "-e", TRUE, sqdARG_INT },
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

  /*********************************************** 
   * Parse command line
   ***********************************************/

  nmodels        = 100;
  seed           = (int) time ((time_t *) NULL);
  threshold      = 0.0001;
  begin_set      = FALSE;
  end_set        = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n") == 0) nmodels        = atoi(optarg);
    else if (strcmp(optname, "-s") == 0) seed           = atoi(optarg);
    else if (strcmp(optname, "-t") == 0) threshold      = atof(optarg);
    else if (strcmp(optname, "-b") == 0) { begin_set = TRUE; spos = atoi(optarg); }
    else if (strcmp(optname, "-e") == 0) { end_set   = TRUE; epos = atoi(optarg); }
  }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
 
  if(begin_set && !end_set || !begin_set && end_set)
    Die("Must use both -b and -e or neither.\n");

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

  /* Determine number of consensus columns modelled by CM */
  ncols = 0;
  for(v = 0; v <= cm->M; v++)
    {
      if(cm->stid[v] ==  MATP_MP)
	ncols += 2;
      else if(cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
	ncols++;
    }

  if(begin_set && end_set)
    {
      /* build one model, as specified on command-line */
      if(spos < 1) spos = 1;
      if(epos > ncols) epos = ncols;
      build_sub_cm(cm, &sub_cm, spos, epos, spos, epos, orig2sub_smap, sub2orig_smap);

      /* check_sub_cm_by_sampling() call builds a CP9 HMM from the sub_cm and checks to make 
       * sure this CP9 HMM is correct. This check is done by sampling a deep MSA from the CM, 
       * truncating it before hmm_start_node and after hmm_end_node and then doing chi-squared
       * tests to see if the samples came from the CP9 HMM distribution.
       */
      check_sub_cm_by_sampling(cm, sub_cm, spos, epos, 0.01, 100000);    
    }
  else
    {
      for(i = 0; i < nmodels; i++)
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
	  printf("$$$$$$$$$$$$$$$ MODELO %d $$$$$$$$$$$$$$$$$$$$\n", i);
	  printf("$$$$$$$$$$$$$$$ S: %d E: %d **#$$$$$$$$$$$$$$$\n", spos, epos);
	  /* Build a sub CM between spos and epos, inclusive */
	  build_sub_cm(cm, &sub_cm, spos, epos, spos, epos, orig2sub_smap, sub2orig_smap);

	  /* check_sub_cm_by_sampling() call builds a CP9 HMM from the sub_cm and checks to make 
	   * sure this CP9 HMM is correct. This check is done by sampling a deep MSA from the CM, 
	   * truncating it before hmm_start_node and after hmm_end_node and then doing chi-squared
	   * tests to see if the samples came from the CP9 HMM distribution.
	   */
	  check_sub_cm_by_sampling(cm, sub_cm, spos, epos, 0.01, 100000);
	}
    }
  FreeCM(cm);
  exit(0);
}

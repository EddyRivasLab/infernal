/* cp9-test.c
 * EPN, Mon Nov 13 17:49:07 2006
 * Easelification: EPN, Fri Nov 30 10:14:17 2007
 * 
 * Test the CM -> CP9 HMM construction procedure.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
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
#include "esl_random.h"
#include "esl_stopwatch.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name        type         default  env  range toggles reqs incomp  help                                            docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     NULL, NULL, "n>0", NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-t",        eslARG_REAL,   "1E-4",NULL, "x>0.",NULL,  NULL, NULL, "probability threshold for reporting violations", 0 },
  { "--psionly", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only check that psi and phi values match",       1 },
  { "--nseq",    eslARG_INT,  "500000",NULL, "n>=1000", NULL,  NULL, NULL, "use <n> samples to build CP9 HMM from",          1 },
  { "--chi",     eslARG_REAL,   ".01", NULL, "x>0.",NULL,  NULL, NULL, "fail sampling check if any chi-square test < <f>", 1},
  { "--dlev",    eslARG_INT,     NULL, NULL, "0<n<4",NULL, NULL, NULL, "set verbosity of debugging print statements to <n>", 1},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "test CP9 HMM construction procedure";

int
main(int argc, char **argv)
{
  ESL_GETOPTS       *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char              *cmfile = esl_opt_GetArg(go, 1);
  CMFILE            *cmfp;        /* open CM file for reading */
  CM_t              *cm;          /* a covariance model       */
  CP9_t             *hmm;         /* constructed CP9 HMM; written to hmmfile              */
  CP9Map_t          *cp9map;      /* maps the hmm to the cm and vice versa */
  double             pthresh;     /* psi threshold for calling violations */
  int                nsamples;    /* Number of samples to build sampled HMM check. */
  float              chi_thresh;  /* if any chi-square test during the sampling check 
				   * is below this threshold, fail. */
  int                do_psionly;  /* don't do a sampling check only compare the expected
				   * number of times each HMM and CM state is entered */
  int                debug_level; /* verbosity of debugging printf statements */
  ESL_STOPWATCH     *w    = NULL; /* for timings */
  ESL_RANDOMNESS    *r    = NULL; /* source of randomness */
  ESL_ALPHABET      *abc  = NULL; /* alphabet, for the CM */

  /*********************************************** 
   * Parse command line
   ***********************************************/
  pthresh        = esl_opt_GetReal   (go, "-t");
  do_psionly     = esl_opt_GetBoolean(go, "--psionly");
  nsamples       = esl_opt_GetInteger(go, "--nseq");
  chi_thresh     = esl_opt_GetReal   (go, "--chi");
  if(esl_opt_IsDefault(go, "--dlev")) debug_level = 0;
  else                                debug_level = esl_opt_GetInteger(go, "--dlev");
 
  /********************************************`*** 
   * Preliminaries: get our CM
   ***********************************************/

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if ((CMFileRead(cmfp, NULL, &abc, &cm)) != eslOK) cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  w  = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  if(!build_cp9_hmm(cm, &hmm, &cp9map, TRUE, pthresh, debug_level))
    cm_Fail("CM Plan 9 HMM fails the psi/phi comparison test.\n");
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "CP9 construction CPU time: ");
  
  if(!do_psionly) {
    esl_stopwatch_Start(w);
    if (! esl_opt_IsDefault(go, "-s")) 
      r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
    else r = esl_randomness_CreateTimeseeded();
    
    if(!(CP9_check_by_sampling(cm, hmm, r,
			       NULL,     /* Don't keep track of failures (sub_cm feature) */
			       1, hmm->M, chi_thresh, nsamples, debug_level)))
	cm_Fail("CP9 HMM fails sampling check!\n");
      else
	printf("CP9 HMM passed sampling check.\n");
      
    esl_stopwatch_Stop(w);
    esl_stopwatch_Display(stdout, w, "CP9 sampling check CPU time: ");
  }
  /* clean up and exit */
  FreeCP9Map(cp9map);
  FreeCPlan9(hmm);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}

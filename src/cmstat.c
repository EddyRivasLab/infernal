/************************************************************
 * @LICENSE@
 ************************************************************/

/* cmstat.c
 * EPN, Tue Aug 21 12:50:34 2007
 *
 * Display summary statistics for a CM or CM database 
 * (such as Rfam). 
 *
 * Based on SRE's hmmstat.c from HMMER3.
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sqio.h"
#include "esl_stats.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

#define CMCUTOPTS    "-E,-T,--ga,--tc,--nc"                       /* exclusive choice for CM cutoff */

#define DEFAULT_DBSIZE 10000000 /* 10 Mb */

static ESL_OPTIONS options[] = {
  /* name           type      default      env  range     toggles      reqs     incomp    help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "display summary statistics for CMs";

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing   */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory error, stopwatch not created.\n");
  esl_stopwatch_Start(w);

  ESL_ALPHABET    *abc = NULL;  /* alphabet                  */
  char            *cmfile;	/* name of input CM file     */ 
  CM_FILE         *cmfp;	/* open input CM file stream */
  CM_t            *cm;          /* CM most recently read     */
  int              ncm;         /* CM index                  */
  char             errbuf[cmERRBUFSIZE]; /* for error messages */
  int              status;      /* easel status */

  /* Process the command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nOptions:");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=docgroup, 2 = indentation; 80=textwidth*/
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if ((cmfile = esl_opt_GetArg(go, 1)) == NULL) 
    {
      puts("Failed to read <cmfile> argument from command line.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  cm_banner(stdout, argv[0], banner);

  /* Initializations: open the CM file
   */
  status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf);
  if      (status == eslENOTFOUND) cm_Fail("File existence/permissions problem in trying to open CM file %s.\n%s\n", cmfile, errbuf);
  else if (status == eslEFORMAT)   cm_Fail("File format problem in trying to open CM file %s.\n%s\n",                cmfile, errbuf);
  else if (status != eslOK)        cm_Fail("Unexpected error %d in opening CM file %s.\n%s\n",               status, cmfile, errbuf);  

  /* Main body: read CMs one at a time, print stats 
   */
  fprintf(stdout, "# %-4s  %-20s  %-12s  %8s  %8s  %4s  %4s  %4s  %4s  %5s  %12s\n",    "",      "",                     "",             "",         "",         "",     "",      "",      "", "",    "rel entropy");
  fprintf(stdout, "# %-4s  %-20s  %-12s  %8s  %8s  %4s  %4s  %4s  %4s  %5s  %12s\n",    "",      "",                     "",             "",         "",         "",     "",      "",      "", "",    "------------");
  fprintf(stdout, "# %-4s  %-20s  %-12s  %8s  %8s  %4s  %4s  %4s  %4s  %5s  %5s  %5s\n", "idx",  "name",                 "accession",    "nseq",     "eff_nseq", "clen", "bps",   "bifs",  "W",     "M",     "CM",     "HMM");
  fprintf(stdout, "# %-4s  %-20s  %-12s  %8s  %8s  %4s  %4s  %4s  %4s  %5s  %5s  %5s\n", "----", "--------------------", "------------", "--------", "--------", "----", "----", "----",   "----", "-----", "-----", "-----");

  ncm = 0;
  while ((status = cm_file_Read(cmfp, TRUE, &abc, &cm)) != eslEOF) 
    {
      if      (status == eslEOD)  cm_Fail("read failed, CM file %s may be truncated?", cmfile);
      else if (status != eslOK)   cm_Fail(cmfp->errbuf);
      ncm++;

      /* build the cp9 HMM */
      if(!(build_cp9_hmm(cm, &(cm->cp9), &(cm->cp9map), FALSE, 0.0001, 0))) cm_Fail("Couldn't build a CP9 HMM from the CM\n");
      
      fprintf(stdout, "%6d  %-20s  %-12s  %8d  %8.2f  %4d  %4d  %4d  %4d  %5d  %5.3f  %5.3f\n",
	      ncm,
	      cm->name,
	      cm->acc == NULL ? "-" : cm->acc,
	      cm->nseq,
	      cm->eff_nseq,
	      cm->clen,
	      CMCountStatetype(cm, MP_st),
	      CMCountStatetype(cm, B_st),
	      cm->W,
	      cm->M,
	      cm_MeanMatchRelativeEntropy(cm),
	      cp9_MeanMatchRelativeEntropy(cm));
      
      FreeCM(cm);
    }
  
  fprintf(stdout, "#\n");
  esl_alphabet_Destroy(abc);
  cm_file_Close(cmfp);
  esl_getopts_Destroy(go);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  esl_stopwatch_Destroy(w);
  exit(0);
}

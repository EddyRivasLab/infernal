#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_randomseq.h>
#include <esl_sqio.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "--cutoff",  eslARG_REAL,   "0.0", NULL, NULL,  NULL,  NULL, NULL, "set bitscore cutoff to <x>",                     0 },
  { "--S_Sa",    eslARG_REAL,  "0.50", NULL, NULL,  NULL,  NULL, NULL, "S emit single residue probability <x>",          0 },
  { "--S_SM",    eslARG_REAL,  "0.20", NULL, NULL,  NULL,  NULL, NULL, "S add model segment probability <x>",            0 },
  { "--S_e",     eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "S end probability (unsettable, 1-S_Sa-S_SM)",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "test driver for MSCYK implementation";

int 
main(int argc, char **argv)
{
  int             status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  CM_t            *cm;
  ESL_ALPHABET   *abc     = NULL;
  int             i;
  float           sc;
  char           *cmfile = esl_opt_GetArg(go, 1);
  char           *seqfile= esl_opt_GetArg(go, 2);
  CMFILE         *cmfp;	/* open input CM file stream */
  ESL_SQFILE     *sqfp;
  ESL_SQ         *seq;
  char            errbuf[cmERRBUFSIZE];
  CM_CONSENSUS   *ccm;
  uint8_t         cutoff;
  float           f_cutoff;
  float           f_S_Sa, f_S_SM, f_S_e, f_M_S;
  int             format = eslSQFILE_UNKNOWN;
  int             max, imax, jmax;
  search_results_t *results = NULL;

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if ((status = CMFileRead(cmfp, errbuf, &abc, &cm) != eslOK))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  if (esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK) cm_Fail("Failed to open sequence database file\n");

  cm->config_opts |= CM_CONFIG_LOCAL;
  cm->search_opts |= CM_SEARCH_NOQDB;
  ConfigCM(cm, errbuf, TRUE, NULL, NULL); /* TRUE says: calculate W */
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE); /* impt to do this after QDBs set up in ConfigCM() */
  ccm = cm_consensus_Convert(cm);

if(ccm == NULL) cm_Fail("ccm NULL!\n");
if(ccm->oesc == NULL) cm_Fail("oesc NULL!\n");

  /* Set meta-model parameters */
  f_S_Sa = esl_opt_GetReal(go, "--S_Sa");
  if (f_S_Sa <= 0 || f_S_Sa >= 1) cm_Fail("S->Sa must be between zero and 1\n");
  f_S_SM = esl_opt_GetReal(go, "--S_SM");
f_S_SM = (cm->clen - f_S_Sa*(cm->clen + 1))/(2*cm->clen + ccm->e_fraglen);
  if (f_S_SM <= 0 || f_S_SM >= 1) cm_Fail("S->SM must be between zero and 1\n");
  f_S_e = 1. - f_S_Sa - f_S_SM;
  if (f_S_e <= 0) cm_Fail("Error: S->e out of range\n");
  ccm->tsb_S_Sa = unbiased_byteify(ccm,sreLOG2(f_S_Sa));
  ccm->tsb_S_SM = unbiased_byteify(ccm,sreLOG2(f_S_SM));
  ccm->tsb_S_e  = unbiased_byteify(ccm,sreLOG2(f_S_e ));
  ccm->tsb_M_S  = unbiased_byteify(ccm,ccm->sc_frag);

  f_cutoff = esl_opt_GetReal(go, "--cutoff");
  /* Need to scale and offset, but not change sign -> switch sign ahead of time */
  cutoff = ccm->base_b + unbiased_byteify(ccm,-f_cutoff);

  seq = esl_sq_Create();
  while ( esl_sqio_Read(sqfp, seq) == eslOK)
  {
    if (seq->n == 0) continue;
    if (seq->dsq == NULL) esl_sq_Digitize(abc, seq);
    fprintf(stdout,"%s\t",seq->name);
  
    results = CreateResults(INIT_RESULTS);
    if((status = SSE_MSCYK(ccm, errbuf, cm->smx->W, seq->dsq, 1, seq->n, cutoff, results, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);

    /* Rudimentary output of results */
    max = 0;
    for (i = 0; i < results->num_results; i++) {
      if ((int)results->data[i].score > max) {
        max = results->data[i].score;
        imax= results->data[i].start;
        jmax= results->data[i].stop;
      }
    }
    fprintf(stdout,"%4d %4d %4d %6f\n",imax,jmax,max,sc);
    fflush(stdout);

    FreeResults(results);

    esl_sq_Destroy(seq);
    seq = esl_sq_Create();
  }
  esl_sq_Destroy(seq);


  if (ccm != NULL) cm_consensus_Free(ccm);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);
  esl_getopts_Destroy(go);

  return 0;

 ERROR:
  cm_Fail("memory allocation error");
  return 0;
}

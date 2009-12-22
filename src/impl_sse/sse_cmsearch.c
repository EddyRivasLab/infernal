#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sqio.h"
#include "esl_stats.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  /* Basic options */
  { "-h",       eslARG_NONE,  NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           1 },
  /* Cutoff level options */
  { "--s1-T",eslARG_REAL,  NULL, NULL,   NULL,  NULL,  NULL, "--s1-E", "set stage 1 cutoff to bitscore <x>",             2 },
  { "--s1-E",eslARG_REAL, "0.1", NULL,"0<x<1",  NULL,  NULL, "--s1-T", "set stage 1 cutoff to probability <x> per kb",   2 },
  { "--s2-T",eslARG_REAL,  NULL, NULL,   NULL,  NULL,  NULL, "--s2-E", "set stage 2 cutoff to bitscore <x>",             2 },
  { "--s2-E",eslARG_REAL,"1e-3", NULL,"0<x<1",  NULL,  NULL, "--s2-T", "set stage 2 cutoff to probability <x>",          2 },
  { "--s3-T",eslARG_REAL,  NULL, NULL,   NULL,  NULL,  NULL, "--s3-E", "set stage 3 cutoff to bitscore <x>",             2 },
  { "--s3-E",eslARG_REAL,"1e-5", NULL,"0<x<1",  NULL,  NULL, "--s3-T", "set stage 3 cutoff to probability <x>",          2 },
  /* Meta-model parameters */
  { "--mmp_auto",eslARG_NONE,"default",NULL, NULL, "--mmp_manual",NULL,NULL,"automatic determination of MSCYK model parameters",3},
  { "--mmp_manual",eslARG_NONE,  FALSE,NULL, NULL, "--mmp_auto",  NULL,NULL,"manual determination of MSCYK model parameters",3},
  { "--S_Sa",    eslARG_REAL,  "0.50", NULL, NULL,  NULL,"--mmp_manual","--mmp_auto", "S emit single residue probability <x>",          3 },
  { "--S_SM",    eslARG_REAL,    NULL, NULL, NULL,  NULL,"--mmp_manual","--mmp_auto", "S add model segment probability <x>",            3 },
  { "--S_e",     eslARG_NONE,    NULL, NULL, NULL,  NULL,"--mmp_manual","--mmp_auto", "S end probability (unsettable, 1-S_Sa-S_SM)",    3 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "multistage SSE pipeline for sequence database search with an RNA CM";

int 
main(int argc, char **argv)
{
  int             status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  CM_t            *cm;
  ESL_ALPHABET   *abc     = NULL;
  int             i;
  float           sc, sc3;
  int             sc2, br2, bsc2;
  char           *cmfile = esl_opt_GetArg(go, 1);
  char           *seqfile= esl_opt_GetArg(go, 2);
  CMFILE         *cmfp;	/* open input CM file stream */
  ESL_SQFILE     *sqfp;
  ESL_SQ         *seq;
  char            errbuf[cmERRBUFSIZE];
  CM_CONSENSUS   *ccm = NULL;
  CM_OPTIMIZED   *ocm = NULL;
  uint8_t         s1_cutoff;
  int16_t         s2_cutoff;
  float           s3_cutoff;
  float           e_cutoff, f_cutoff;
  float           f_S_Sa, f_S_SM, f_S_e;
  int             format = eslSQFILE_UNKNOWN;
  int             max, imax, jmax;

  /* vars for MSCYK score distribution fitting */
  int L = 1000;	/* P-values will be per-kb, so simulate 1 kb at a time */
  int samples = 10000;
  float *x;
  double *xv;
  int n;
  double ccm_mu, ccm_lambda;
  double score, eval;
  double param[2];
  ESL_RANDOMNESS *r;
  double background[4] = {0.25, 0.25, 0.25, 0.25};
  ESL_HISTOGRAM *hist;
    
  search_results_t *results = NULL; /* First stage results from MSCYK */
  search_results_t *windows = NULL; /* Expanded and merged hit candidate windows after first stage */

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if ((status = CMFileRead(cmfp, errbuf, &abc, &cm) != eslOK))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  if (esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK) cm_Fail("Failed to open sequence database file\n");

  cm->config_opts |= CM_CONFIG_LOCAL;
  cm->search_opts |= CM_SEARCH_NOQDB;
  ConfigCM(cm, errbuf, TRUE, NULL, NULL); /* TRUE says: calculate W */
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE); /* impt to do this after QDBs set up in ConfigCM() */
  ccm = cm_consensus_Convert(cm);
  ocm = cm_optimized_Convert(cm);

  if(ccm == NULL) cm_Fail("ccm NULL!\n");
  if(ccm->oesc == NULL) cm_Fail("oesc NULL!\n");

  /* Set meta-model parameters */
  if (esl_opt_IsOn(go, "--mmp_auto")) {
    f_S_Sa = esl_opt_GetReal(go, "--S_Sa");
    f_S_SM = (cm->clen - f_S_Sa*(cm->clen + 1))/(2*cm->clen + ccm->e_fraglen);
    ccm->r = ((float) cm->clen)/(cm->clen+1.);
    f_S_e = 1. - f_S_Sa - f_S_SM;
  }
  else {
    float nullL;
    f_S_Sa = esl_opt_GetReal(go, "--S_Sa");
    if (esl_opt_IsOn(go, "--S_SM")) f_S_SM = esl_opt_GetReal(go, "--S_SM");
    else                            f_S_SM = 0.24;
    f_S_e = 1. - f_S_Sa - f_S_SM;
    nullL = MSCYK_explen(ccm->e_fraglen,f_S_Sa,f_S_SM,f_S_e);
    ccm->r = nullL/(nullL+1.);
  }
  if (f_S_Sa <= 0 || f_S_Sa >= 1) cm_Fail("S->Sa must be between zero and 1\n");
  if (f_S_SM <= 0 || f_S_SM >= 1) cm_Fail("S->SM must be between zero and 1\n");
  if (f_S_e <= 0)                 cm_Fail("Error: S->e out of range\n");
  ccm->tsb_S_Sa = unbiased_byteify(ccm,sreLOG2(f_S_Sa));
  ccm->tsb_S_SM = unbiased_byteify(ccm,sreLOG2(f_S_SM));
  ccm->tsb_S_e  = unbiased_byteify(ccm,sreLOG2(f_S_e ));
  ccm->tsb_M_S  = unbiased_byteify(ccm,ccm->sc_frag);

  /* Set filtering cutoff parameters */
  if (esl_opt_IsOn(go, "--s1-T")) f_cutoff = esl_opt_GetReal(go, "--s1-T");
  else {
    /* Simulation to determine score distribution tail for MSCYK/CM combined model */
    hist = esl_histogram_CreateFull(-50., 50., 0.5);
    seq = esl_sq_Create();
    r = esl_randomness_CreateTimeseeded();
    ESL_ALLOC(x, sizeof(float) * samples);
    ESL_ALLOC(seq->dsq, sizeof(ESL_DSQ) * (L+2));
    for (int k = 0; k < samples; k++) {
      esl_rsq_xIID(r, background, 4, L, seq->dsq);
      if((status = SSE_MSCYK(ccm, errbuf, cm->smx->W, seq->dsq, 1, L, -45, NULL, FALSE, NULL, &(x[k]))) != eslOK) cm_Fail(errbuf);
      esl_histogram_Add(hist, x[k]);
    }
    esl_histogram_GetTailByMass(hist, 0.1, &xv, &n, NULL);
    esl_exp_FitComplete(xv, n, &ccm_mu, &ccm_lambda);
    param[0] = ccm_mu; param[1] = ccm_lambda;
    esl_histogram_SetExpectedTail(hist, ccm_mu, 0.1, &esl_exp_generic_cdf, &param);
    esl_histogram_PlotSurvival(stdout, hist);
    printf("ccm_mu = %7.3f\t ccm_lambda = %7.3f\n",ccm_mu,ccm_lambda);
/*
esl_histogram_GetRank(hist, 10, &score);
eval  = esl_exp_surv(score, ccm_mu, ccm_lambda);
printf("E-val at rank 10 = %1.3e\t",eval);
*/

    e_cutoff = esl_opt_GetReal(go,"--s1-E");
    f_cutoff = esl_exp_invcdf(e_cutoff,ccm_mu,ccm_lambda);
    fprintf(stderr,"Stage 1: P-value cutoff %.2e per kb -> score cutoff %.2f\n",e_cutoff,f_cutoff);
    free(x);
    esl_randomness_Destroy(r);
    esl_sq_Destroy(seq);
    esl_histogram_Destroy(hist);
  }
  /* Need to scale and offset, but not change sign -> switch sign ahead of time */
  s1_cutoff = ccm->base_b + unbiased_byteify(ccm,-f_cutoff);

  if (esl_opt_IsOn(go, "--s2-T")) f_cutoff = esl_opt_GetReal(go, "--s2-T");
  else {
    f_cutoff = 0.0;
    fprintf(stderr,"WARNING: P-value cutoffs not implemented, setting stage 2 bitscore cutoff to %.1f\n",f_cutoff);
  }
  /* Need to scale */
  s2_cutoff = wordify(ocm->scale_w, f_cutoff);

  if (esl_opt_IsOn(go, "--s3-T")) f_cutoff = esl_opt_GetReal(go, "--s3-T");
  else {
    e_cutoff = esl_opt_GetReal(go,"--s3-E");
//FIXME!!  Always picking the first partition, here, which probably isn't what we want.
//FIXME!!  We don't necessarily know how %GC space was divided for calibration.
//FIXME!!  Also, need to actually check that the cm has been calibrated before we start pulling data from cm->stats
    //f_cutoff = esl_gumbel_invcdf(e_cutoff,cm->stats->expAA[EXP_CM_LC][0]->mu_extrap,cm->stats->expAA[EXP_CM_LC][0]->lambda);
    E2MinScore(cm,errbuf,EXP_CM_LC,e_cutoff,&f_cutoff);
    //fprintf(stderr,"WARNING: P-value cutoffs not implemented, setting stage 3 bitscore cutoff to %.1f\n",f_cutoff);
    fprintf(stderr,"Stage 3: P-value cutoff %.2e -> score cutoff %.2f\n",e_cutoff,f_cutoff);
  }
  s3_cutoff = f_cutoff;

  seq = esl_sq_Create();
  while ( esl_sqio_Read(sqfp, seq) == eslOK)
  {
    if (seq->n == 0) continue;
    if (seq->dsq == NULL) esl_sq_Digitize(abc, seq);
    fprintf(stdout,"%s\t",seq->name);
  
    /* Stage 1: MSCYK */
    results = CreateResults(INIT_RESULTS);
    if((status = SSE_MSCYK(ccm, errbuf, cm->smx->W, seq->dsq, 1, seq->n, s1_cutoff, results, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);

    /* Rudimentary output of results */
    max = 0; imax = -1; jmax = -1;
    for (i = 0; i < results->num_results; i++) {
      if ((int)results->data[i].score > max) {
        max = results->data[i].score;
        imax= results->data[i].start;
        jmax= results->data[i].stop;
      }
    }
    fprintf(stdout,"%4d %4d %4d %6f\n",imax,jmax,max,sc);
    fflush(stdout);

    /* Convert hits to windows for next stage */
    windows = ResolveMSCYK(results, 1, seq->n, cm->smx->W);
    FreeResults(results);

    for (i = 0; i < windows->num_results; i++) {
      /* Stage 2: medium-precision CYK */
      sc2 = SSE_CYKFilter_epi16(ocm, seq->dsq, seq->n, 0, ocm->M, windows->data[i].start, windows->data[i].stop, TRUE, &br2, &bsc2);

      if (sc2 > s2_cutoff) {
        /* Stage 3: full-precision CYK */
        sc3 = SSE_CYKInsideScore(cm, seq->dsq, seq->n, 0, windows->data[i].start, windows->data[i].stop);
      }
    }

    esl_sq_Destroy(seq);
    seq = esl_sq_Create();
  }
  esl_sq_Destroy(seq);


  if (ccm != NULL) cm_consensus_Free(ccm);
  if (ocm != NULL) cm_optimized_Free(ocm); free(ocm);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);
  esl_getopts_Destroy(go);

  return 0;

 ERROR:
  cm_Fail("memory allocation error");
  return 0;
}

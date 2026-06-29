/* test_wviterbi.c -- brief 169 windowed-Viterbi band harness.
 *
 * The brief-169 windowed-Viterbi band is, by construction:
 *     band = (Viterbi MAP trace i2k)  +/- (per-node pad nodepad[k])
 * The MAP trace i2k is exactly the per-row argmax-k pin the IBV deriver
 * already produces (ibv_through_scan).  The per-node pad has two sources:
 *   --calib : the brief-169 F+B-halfwidth p95 calibration (cm_ComputeP7WVNodePad),
 *             which reproduces the prototype's pn_p95 band (the GO config); OR
 *   default : cm->p7_cm_nodepad -- a DIFFERENT (much narrower) calibration
 *             (cm_ComputeP7CMNodePad: Viterbi-pin deficit p99), which REGRESSES
 *             accuracy here.  Use --calib for the real WV band.
 *
 * This tool:
 *   (1) derives i2k via p7_Seq2BandsIBV[_dnc]  (the exact MAP-trace oracle),
 *   (2) builds the WV band via p7_pins2bands_nodepad(i2k, nodepad+padplus),
 *   (3) optionally dumps each seq's band in the PB_LOAD_BAND TSV format so it
 *       can be fed through the brief-168 cmalign pinbridge accuracy harness,
 *   (4) reports band stats and (with -c) WV-vs-D&C i2k agreement.
 *
 *   gcc -O3 -msse2 -I. -I../easel -I../hmmer/src ... -o test_wviterbi \
 *       test_wviterbi.c libinfernal.a ../hmmer/src/libhmmer.a ../easel/libeasel.a -lm
 */
#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "infernal.h"

static ESL_OPTIONS options[] = {
  { "-h", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show help",                              0 },
  { "-o", eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL, "dump per-seq band TSVs to dir <s>",      0 },
  { "-p", eslARG_INT,     "0", NULL, "n>=0",NULL, NULL, NULL, "padplus: add <n> to every node pad",     0 },
  { "--flat", eslARG_NONE,FALSE,NULL, NULL, NULL, NULL, NULL, "use flat IBV (not D&C) for i2k",         0 },
  { "--wv",   eslARG_NONE,FALSE,NULL, NULL, NULL, NULL, NULL, "use windowed-Viterbi kernel for i2k",    0 },
  { "-c", eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "cross-check WV i2k vs D&C i2k",          0 },
  { "--trunc", eslARG_NONE,FALSE,NULL, NULL, NULL, NULL, NULL, "brief171: Tgm begin/end-anywhere bands", 0 },
  { "--delta", eslARG_INT,"20000",NULL,"n>=0",NULL,NULL,NULL,"IBV delta milli-bits (for delta band)",  0 },
  { "--calib", eslARG_NONE,FALSE,NULL, NULL, NULL, NULL, NULL, "use F+B-halfwidth WV pad (not cm nodepad)",0 },
  { "--nsamp", eslARG_INT, "40", NULL, "n>0", NULL, NULL, NULL, "calib: # CM-emitted samples",          0 },
  { "--q",     eslARG_REAL,"0.99",NULL,"0<x<=1",NULL,NULL,NULL,"calib: pad quantile",                   0 },
  { "--floor", eslARG_INT,  "2", NULL, "n>=0",NULL, NULL, NULL, "calib: floor pad",                     0 },
  { "--seed",  eslARG_INT,"181", NULL, "n>=0",NULL, NULL, NULL, "calib: RNG seed",                      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "brief169 windowed-Viterbi band harness";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char         *cmfile  = esl_opt_GetArg(go, 1);
  char         *seqfile = esl_opt_GetArg(go, 2);
  char         *outdir  = esl_opt_IsOn(go, "-o") ? esl_opt_GetString(go, "-o") : NULL;
  int           padplus = esl_opt_GetInteger(go, "-p");
  int           use_flat= esl_opt_GetBoolean(go, "--flat");
  int           use_wv  = esl_opt_GetBoolean(go, "--wv");
  int           do_cc   = esl_opt_GetBoolean(go, "-c");
  int           do_trunc= esl_opt_GetBoolean(go, "--trunc");
  int           delta   = esl_opt_GetInteger(go, "--delta");
  char          errbuf[eslERRBUFSIZE];

  ESL_ALPHABET *abc  = NULL;
  CM_FILE      *cmfp = NULL;
  CM_t         *cm   = NULL;
  P7_HMM       *hmm  = NULL;
  ESL_SQFILE   *sqfp = NULL;
  ESL_SQ       *sq   = NULL;
  int           status;
  int          *nodepad = NULL;
  int           M, k;

  if (cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf) != eslOK) p7_Fail("cm_file_Open: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm)              != eslOK) p7_Fail("cm_file_Read failed");
  cm_file_Close(cmfp);
  if (cm_Configure(cm, errbuf, -1) != eslOK) p7_Fail("cm_Configure: %s", errbuf);
  if (cm->fp7 == NULL) p7_Fail("CM has no fp7 filter HMM");
  hmm = cm->fp7;
  M   = hmm->M;

  if (esl_opt_GetBoolean(go, "--calib")) {
    ESL_RANDOMNESS *rng = esl_randomness_Create((uint32_t) esl_opt_GetInteger(go, "--seed"));
    int *wvpad = NULL;
    if (cm_ComputeP7WVNodePad(cm, errbuf, rng, esl_opt_GetInteger(go, "--nsamp"),
                              esl_opt_GetReal(go, "--q"), 20000,
                              esl_opt_GetInteger(go, "--floor"), &wvpad) != eslOK)
      p7_Fail("cm_ComputeP7WVNodePad: %s", errbuf);
    esl_randomness_Destroy(rng);
    ESL_ALLOC(nodepad, sizeof(int) * (M + 1));
    for (k = 0; k <= M; k++) nodepad[k] = wvpad[k] + padplus;
    free(wvpad);
  } else {
    if (! (cm->flags & CMH_P7NODEPAD))
      p7_Fail("CM has no P7NODEPAD; use --calib to compute the WV pad instead");
    ESL_ALLOC(nodepad, sizeof(int) * (M + 1));
    for (k = 0; k <= M; k++) nodepad[k] = cm->p7_cm_nodepad[k] + padplus;
  }

  if (esl_sqfile_Open(seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) p7_Fail("open seqfile failed");
  sq = esl_sq_CreateDigital(abc);

  printf("# CM=%s  M=%d  padplus=%d  i2k_src=%s\n", cmfile, M,
         padplus, use_wv ? "wv" : (use_flat ? "flat" : "dnc"));
  printf("# %-26s %7s %6s %9s %9s %7s %7s\n",
         "seq", "L", "npins", "wv_ncells", "dlt_ncells", "wv_abw", "i2kdiff");

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int   L = sq->n;
    int  *i2k=NULL,*kmin_d=NULL,*kmax_d=NULL,nc_d=0;
    int  *i2k_wv=NULL,*kmin=NULL,*kmax=NULL,nc=0;
    int  *i2k_dnc=NULL,*kd2=NULL,*kx2=NULL,ncd2=0;
    int   npins=0, i, i2kdiff=-1;

    /* (1) derive i2k (the MAP trace) */
    if (use_wv) {
      if ((status = p7_Seq2BandsWV(cm, errbuf, sq->dsq, L, nodepad, do_trunc,
                                   &i2k, &kmin, &kmax, &nc)) != eslOK)
        p7_Fail("p7_Seq2BandsWV failed on %s: %s", sq->name, errbuf);
      /* WV returns the nodepad band directly; also keep i2k for cross-check. */
    } else {
      if (use_flat)
        status = p7_Seq2BandsIBV(cm, errbuf, sq->dsq, L, delta, do_trunc,
                                 P7IBV_MODE_DELTA, 0, &i2k, &kmin_d, &kmax_d, &nc_d);
      else
        status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, L, delta, 0, FALSE, do_trunc,
                                     P7IBV_MODE_DELTA, 0, &i2k, &kmin_d, &kmax_d, &nc_d);
      if (status != eslOK) p7_Fail("IBV deriver failed on %s: %s", sq->name, errbuf);

      /* (2) build the WV band from i2k + nodepad (i2k is pruned in place; copy) */
      ESL_ALLOC(i2k_wv, sizeof(int) * (L + 1));
      memcpy(i2k_wv, i2k, sizeof(int) * (L + 1));
      if ((status = p7_pins2bands_nodepad(i2k_wv, errbuf, L, M, nodepad, 0,
                                          &kmin, &kmax, &nc)) != eslOK)
        p7_Fail("p7_pins2bands_nodepad failed on %s: %s", sq->name, errbuf);
    }

    for (i = 1; i <= L; i++) if (i2k[i] != -1) npins++;

    /* optional cross-check WV i2k vs D&C i2k */
    if (do_cc) {
      if ((status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, L, delta, 0, FALSE, do_trunc,
                                        P7IBV_MODE_DELTA, 0, &i2k_dnc, &kd2, &kx2, &ncd2)) != eslOK)
        p7_Fail("D&C oracle failed on %s: %s", sq->name, errbuf);
      i2kdiff = 0;
      for (i = 1; i <= L; i++) if (i2k[i] != i2k_dnc[i]) i2kdiff++;
    }

    /* (3) dump band TSV in PB_LOAD_BAND format */
    if (outdir != NULL) {
      char path[2048];
      char safe[1024];
      int  z;
      for (z = 0; sq->name[z] != '\0' && z < (int)sizeof(safe)-1; z++)
        safe[z] = (sq->name[z] == '/') ? '_' : sq->name[z];
      safe[z] = '\0';
      snprintf(path, sizeof(path), "%s/%s.band.tsv", outdir, safe);
      FILE *fp = fopen(path, "w");
      if (fp == NULL) p7_Fail("cannot open %s for writing", path);
      fprintf(fp, "# M=%d L=%d\n# i\tkmin\tkmax\n", M, L);
      for (i = 1; i <= L; i++) fprintf(fp, "%d\t%d\t%d\n", i, kmin[i], kmax[i]);
      fclose(fp);
    }

    double wv_abw = (double) nc / (double) (L > 0 ? L : 1);
    printf("  %-26s %7d %6d %9d %9d %7.1f %7d\n",
           sq->name, L, npins, nc, nc_d, wv_abw, i2kdiff);

    free(i2k); if (i2k_wv) free(i2k_wv);
    if (kmin_d) free(kmin_d); if (kmax_d) free(kmax_d);
    free(kmin); free(kmax);
    if (i2k_dnc) free(i2k_dnc); if (kd2) free(kd2); if (kx2) free(kx2);
    esl_sq_Reuse(sq);
  }

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  free(nodepad);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  p7_Fail("allocation failure");
  return status;
}

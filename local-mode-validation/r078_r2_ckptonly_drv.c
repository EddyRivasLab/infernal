/* r078_r2_ckptonly_drv.c : brief 078 R2 -- ckpt-only timing probe (no D&C
 * oracle), for cases where TrCYKDivideAndConquerHB itself is too slow to
 * finish in reasonable time (e.g. wide-band divergent genome-scale seqs) --
 * see brief 076's own finding that D&C's cost grows badly with band width.
 * Reports cm_CheckptTrCYKAlignHB's own wall-time and score/mode only. */
#include <esl_config.h>
#include <p7_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "infernal.h"

static ESL_OPTIONS options[] = {
  { "-h",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "help", 0 },
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

static double now_sec(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return (double)ts.tv_sec + (double)ts.tv_nsec/1.0e9; }

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK) { printf("bad opts\n"); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2) { printf("usage: r078_r2_ckptonly_drv [--tau x] <cmfile> <seqfile>\n"); exit(0); }
  float tau = esl_opt_GetReal(go, "--tau");
  float size_limit = esl_opt_GetReal(go, "--mxsize");
  char *cmfile = esl_opt_GetArg(go, 1);
  char *seqfile = esl_opt_GetArg(go, 2);
  int pass_idx = PLI_PASS_5P_AND_3P_ANY;

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->config_opts |= CM_CONFIG_TRUNC;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  printf("# CM=%s M=%d clen=%d tau=%g\n", cm->name, cm->M, cm->clen, tau);

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n; ESL_DSQ *dsq = sq->dsq;
    double t0 = now_sec();
    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);
    double t_bands = now_sec() - t0;

    char mode = TRMODE_UNKNOWN; float sc = 0.; Parsetree_t *tr = NULL;
    t0 = now_sec();
    if ((status = cm_CheckptTrCYKAlignHB(cm, errbuf, dsq, L, size_limit, pass_idx, &tr, &mode, &sc)) != eslOK)
      cm_Fail("cm_CheckptTrCYKAlignHB: %s", errbuf);
    double t_ckpt = now_sec() - t0;

    printf("# CKPTONLY seq=%s L=%d mode=%s sc=%.6f nodes=%d t_bands=%.4fs t_ckpt=%.4fs\n",
           sq->name, L, MarginalMode(mode), sc, tr->n, t_bands, t_ckpt);
    FreeParsetree(tr);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

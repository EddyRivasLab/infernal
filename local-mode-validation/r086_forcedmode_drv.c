/* r086_forcedmode_drv.c : brief 26_0610-086 validation.
 *
 * Directly exercises each marginal plane of TrCYKDivideAndConquerHB (D&C) by
 * FORCING each root-valid mode (J/L/R/T) and comparing byte-exact against the
 * trusted monolithic oracle cm_TrAlignHB() forced to the SAME preset_mode.
 * This isolates the J-plane fix (086) rather than letting the argmax hide it
 * (the argmax winner may be a mode that was already correct).
 *
 * Gate per (seq,mode): |sc_dnc - sc_oracle| <= 0.01 AND ParsetreeCompare == 0.
 *
 * Usage: r086_forcedmode_drv [-g] [--tau x] <cmfile> <seqfile>
 */
#include <esl_config.h>
#include <p7_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "hmmer.h"
#include "infernal.h"

static ESL_OPTIONS options[] = {
  { "-h",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "help", 0 },
  { "-g",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "global config (default local)", 0 },
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { "-v",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "verbose per-seq", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: r086_forcedmode_drv [-g] [--tau x] [-v] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  int    verbose   = esl_opt_GetBoolean(go, "-v");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);
  int    pass_idx  = PLI_PASS_5P_AND_3P_ANY;

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->config_opts |= CM_CONFIG_TRUNC;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();

  printf("# CM=%s M=%d mode=%s tau=%g LOCAL_BEGIN=%d LOCAL_END=%d\n", cm->name, cm->M,
         do_global?"global":"local", tau, (cm->flags&CMH_LOCAL_BEGIN)?1:0, (cm->flags&CMH_LOCAL_END)?1:0);

  char modes[4] = { TRMODE_J, TRMODE_L, TRMODE_R, TRMODE_T };
  int n_seq = 0, n_checks = 0, n_fail = 0, n_maxdev_seq = 0;
  float maxdev = 0.;

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;
    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);
    n_seq++;

    int m;
    for (m = 0; m < 4; m++) {
      char mode = modes[m];
      /* only root-valid modes */
      if (mode==TRMODE_J && !cm->cp9b->Jvalid[0]) continue;
      if (mode==TRMODE_L && !cm->cp9b->Lvalid[0]) continue;
      if (mode==TRMODE_R && !cm->cp9b->Rvalid[0]) continue;
      if (mode==TRMODE_T && !cm->cp9b->Tvalid[0]) continue;

      /* D&C forced mode */
      char dnc_mode = TRMODE_UNKNOWN; Parsetree_t *tr_dnc = NULL;
      float sc_dnc = TrCYKDivideAndConquerHB(cm, dsq, L, 0, 1, L, pass_idx, mode, &dnc_mode, &tr_dnc, cm->cp9b);

      /* oracle forced mode */
      CM_TR_HB_MX *mx = cm_tr_hb_mx_Create(cm);
      CM_TR_HB_SHADOW_MX *shmx = cm_tr_hb_shadow_mx_Create(cm);
      Parsetree_t *tr_or = NULL; char or_mode = TRMODE_UNKNOWN; float sc_or = IMPOSSIBLE;
      status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, mode, pass_idx, FALSE, FALSE,
                            mx, shmx, NULL, NULL, NULL, NULL, &tr_or, &or_mode, NULL, &sc_or);
      if (status != eslOK) cm_Fail("oracle forced mode %c failed seq %s: %s", "JLRT"[m], sq->name, errbuf);

      float dev = fabsf(sc_dnc - sc_or);
      int peq = ParsetreeCompare(tr_or, tr_dnc);
      int ok = (dev <= 0.01) && peq;
      n_checks++;
      if (dev > maxdev) maxdev = dev;
      if (!ok) {
        n_fail++;
        printf("# FAIL seq=%s mode=%c sc_dnc=%.6f sc_oracle=%.6f dev=%.3e ParsetreeCompare=%s nodes_dnc=%d nodes_or=%d\n",
               sq->name, "JLRT"[m], sc_dnc, sc_or, dev, peq?"IDENTICAL":"DIFFER", tr_dnc?tr_dnc->n:-1, tr_or?tr_or->n:-1);
      } else if (verbose) {
        printf("# ok   seq=%s mode=%c sc=%.6f nodes=%d\n", sq->name, "JLRT"[m], sc_dnc, tr_dnc?tr_dnc->n:-1);
      }
      if (tr_dnc) FreeParsetree(tr_dnc);
      if (tr_or)  FreeParsetree(tr_or);
      cm_tr_hb_mx_Destroy(mx); cm_tr_hb_shadow_mx_Destroy(shmx);
    }
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  printf("# SUMMARY: n_seq=%d n_checks=%d n_fail=%d maxdev=%.3e VERDICT=%s\n",
         n_seq, n_checks, n_fail, maxdev, (n_fail==0)?"PASS":"FAIL");
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

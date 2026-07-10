/* r078_step0_moderesolve_drv.c : brief 26_0610-078 Step 0 insurance check.
 *
 * cm_alndata.c's do_trckpt_r4 branch (~502-536) resolves the marginal mode
 * for a STRUCTURED (bps>0) truncated --ckpt alignment by running
 * TrCYKDivideAndConquerHB() once per root-valid candidate mode (J/L/R/T) and
 * taking the argmax over the returned scores.  The in-tree comment at
 * cm_alndata.c:508-510 claims "each TrCYKDivideAndConquerHB() score == the
 * oracle cm_TrCYKInsideAlignHB()'s {J,L,R,T}alpha[0][L][L], so the argmax ==
 * stock's CYK mode resolution" -- but per brief 078 this claim was never
 * independently re-verified.  This driver checks it empirically: it runs
 * (a) today's ncand-loop (mirrors rl5_ckpttrpinalign_drv.c's
 * resolve_trpins()) and (b) cm_TrCYKInsideAlignHB(preset_mode=TRMODE_UNKNOWN)
 * -- ONE combined-fill call -- on the same sequence, and compares resolved
 * mode + score.
 *
 * Usage: r078_step0_moderesolve_drv [-g] [--tau x] <cmfile> <seqfile>
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
  { 0,0,0,0,0,0,0,0,0,0 },
};

/* today's production ncand-loop, mirrors cm_alndata.c's do_trckpt_r4 (~502-536)
 * and rl5_ckpttrpinalign_drv.c's resolve_trpins() (duplicated here, file-static there too) */
static int
loop_resolve(CM_t *cm, ESL_DSQ *dsq, int L, int pass_idx, char *ret_mode, float *ret_sc, int *ret_ncand)
{
  Parsetree_t *tr_best = NULL;
  char r4_mode = TRMODE_UNKNOWN;
  float r4_cyk = IMPOSSIBLE;
  char cand[4]; int ncand = 0, m;
  if (cm->cp9b->Jvalid[0]) cand[ncand++] = TRMODE_J;
  if (cm->cp9b->Lvalid[0]) cand[ncand++] = TRMODE_L;
  if (cm->cp9b->Rvalid[0]) cand[ncand++] = TRMODE_R;
  if (cm->cp9b->Tvalid[0]) cand[ncand++] = TRMODE_T;
  for (m = 0; m < ncand; m++) {
    Parsetree_t *tr_m = NULL; char mm = TRMODE_UNKNOWN;
    float sc_m = TrCYKDivideAndConquerHB(cm, dsq, L, 0, 1, L, pass_idx, cand[m], &mm, &tr_m, cm->cp9b);
    if (sc_m > r4_cyk) { r4_cyk = sc_m; r4_mode = cand[m]; if (tr_best) FreeParsetree(tr_best); tr_best = tr_m; }
    else if (tr_m) FreeParsetree(tr_m);
  }
  if (tr_best == NULL) return eslFAIL;
  FreeParsetree(tr_best);
  *ret_mode = r4_mode; *ret_sc = r4_cyk; *ret_ncand = ncand;
  return eslOK;
}

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: r078_step0_moderesolve_drv [-g] [--tau x] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);
  int    pass_idx  = PLI_PASS_5P_AND_3P_ANY;

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->config_opts |= CM_CONFIG_TRUNC;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  cm->align_opts  |= CM_ALIGN_OPTACC;
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  int M = cm->M;

  int nb = 0, v;
  for (v = 0; v < M; v++) if (cm->sttype[v] == B_st) nb++;
  printf("# CM=%s M=%d clen=%d bifs=%d mode=%s tau=%g\n",
         cm->name, cm->M, cm->clen, nb, do_global ? "global" : "local(full)", tau);

  int n_total = 0, n_agree_mode = 0, n_agree_sc = 0;

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    char loop_mode = TRMODE_UNKNOWN; float loop_sc = IMPOSSIBLE; int ncand = 0;
    int lstatus = loop_resolve(cm, dsq, L, pass_idx, &loop_mode, &loop_sc, &ncand);

    char oracle_mode = TRMODE_UNKNOWN; float oracle_sc = IMPOSSIBLE; int b = -1;
    int ostatus = cm_TrCYKInsideAlignHB(cm, errbuf, dsq, L, size_limit, TRMODE_UNKNOWN, pass_idx,
                                         cm->trhb_mx, cm->trhb_shmx, &b, &oracle_mode, &oracle_sc);

    n_total++;
    if (lstatus == eslOK && ostatus == eslOK) {
      float dsc = fabsf(loop_sc - oracle_sc);
      int mode_agree = (loop_mode == oracle_mode);
      /* tolerance matches this project's own established D&C-vs-full-matrix-fill
       * convention (hbdnc_drv.c:192, fabs(scB-scA) > 0.01 is the mismatch gate) --
       * TrCYKDivideAndConquerHB (pass-1) and cm_TrCYKInsideAlignHB (full combined
       * fill, the oracle here) sum terms in different orders, so exact float
       * equality is not expected even when both are correct. */
      int sc_agree   = (dsc <= 0.01);
      if (mode_agree) n_agree_mode++;
      if (sc_agree)   n_agree_sc++;
      printf("# STEP0 seq=%s L=%d ncand=%d loop_mode=%s loop_sc=%.6f oracle_mode=%s oracle_sc=%.6f |dsc|=%.3e mode_agree=%s sc_agree=%s\n",
             sq->name, L, ncand, MarginalMode(loop_mode), loop_sc, MarginalMode(oracle_mode), oracle_sc, dsc,
             mode_agree?"YES":"NO", sc_agree?"YES":"NO");
    } else {
      printf("# STEP0 seq=%s L=%d ERROR lstatus=%d ostatus=%d (%s)\n", sq->name, L, lstatus, ostatus, errbuf);
    }
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  printf("# STEP0 SUMMARY: n=%d mode_agree=%d sc_agree=%d VERDICT=%s\n",
         n_total, n_agree_mode, n_agree_sc, (n_agree_mode==n_total && n_agree_sc==n_total) ? "PASS" : "FAIL");
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

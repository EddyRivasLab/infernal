/* r078_r2_ckpttrcyk_drv.c : brief 26_0610-078 R2 validation.
 *
 * Validates cm_CheckptTrCYKAlignHB() -- the new checkpointed, COMBINED-MODE
 * (J/L/R/T) bifurcation (k*)-DISCOVERING truncated CYK max-DP engine
 * (cm_dpalign_trunc.c, brief 078 R2) -- against TODAY'S PRODUCTION pass-1:
 * cm_alndata.c's do_trckpt_r4 ncand-loop, which runs TrCYKDivideAndConquerHB
 * once per root-valid marginal mode candidate and takes the argmax (mirrored
 * here as loop_resolve(), duplicated from Step 0's driver / rl5's
 * resolve_trpins(), since it's file-static in cm_alndata.c).
 *
 * GLOBAL, structured (bps>0) CMs only (R2 scope).
 *
 * Gate: resolved mode MUST match; score agreement within this project's own
 * established D&C-vs-full-fill tolerance (0.01, hbdnc_drv.c's convention --
 * NOT literal bitwise equality); ParsetreeCompare()==0 (byte-identical
 * parse, including marginal mode per node AND every discovered k*, per
 * rl4/rl5's convention -- ParsetreeCompare checks tr->mode[] too).
 *
 * Also reports wall-time (via a simple in-process clock loop over N reps)
 * for both engines, to gauge R2's improvement over the D&C-loop cost per
 * brief 076's 6.6x-21.6x baseline measurement.
 *
 * Usage: r078_r2_ckpttrcyk_drv [--tau x] [--reps N] <cmfile> <seqfile>
 */
#include <esl_config.h>
#include <p7_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
  { "--reps", eslARG_INT,  "1",   NULL, NULL, NULL, NULL, NULL, "wall-time reps per seq", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

/* today's production ncand-loop (cm_alndata.c ~502-536 / rl5's resolve_trpins) */
static int
loop_resolve(CM_t *cm, ESL_DSQ *dsq, int L, int pass_idx, char *ret_mode, float *ret_sc, Parsetree_t **ret_tr)
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
  *ret_mode = r4_mode; *ret_sc = r4_cyk; *ret_tr = tr_best;
  return eslOK;
}

static double
now_sec(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double) ts.tv_sec + (double) ts.tv_nsec / 1.0e9;
}

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: r078_r2_ckpttrcyk_drv [--tau x] [--reps N] <cmfile> <seqfile>\n"); exit(0); }
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  int    reps      = esl_opt_GetInteger(go, "--reps");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);
  int    pass_idx  = PLI_PASS_5P_AND_3P_ANY;

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->config_opts |= CM_CONFIG_TRUNC;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  /* R2 scope: GLOBAL only -- do NOT set CM_CONFIG_LOCAL/HMMLOCAL/HMMEL */
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  int M = cm->M;

  int nb = 0, v;
  for (v = 0; v < M; v++) if (cm->sttype[v] == B_st) nb++;
  if (nb < 1) cm_Fail("CM has no bifurcations (bps=0): not an R2 target");
  printf("# CM=%s M=%d clen=%d bifs=%d mode=global tau=%g reps=%d\n", cm->name, cm->M, cm->clen, nb, tau, reps);

  int n_total = 0, n_pass = 0, n_mode_agree = 0;
  double t_loop_tot = 0., t_ckpt_tot = 0.;

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- oracle: today's production ncand-loop ---- */
    char loop_mode = TRMODE_UNKNOWN; float loop_sc = IMPOSSIBLE; Parsetree_t *tr_loop = NULL;
    double t0 = now_sec();
    int r;
    for (r = 0; r < reps; r++) {
      if (tr_loop) FreeParsetree(tr_loop);
      if (loop_resolve(cm, dsq, L, pass_idx, &loop_mode, &loop_sc, &tr_loop) != eslOK)
        cm_Fail("loop_resolve failed for seq %s", sq->name);
    }
    double t_loop = (now_sec() - t0) / reps;

    /* ---- R2 engine: cm_CheckptTrCYKAlignHB ---- */
    char ckpt_mode = TRMODE_UNKNOWN; float ckpt_sc = 0.; Parsetree_t *tr_ckpt = NULL;
    t0 = now_sec();
    for (r = 0; r < reps; r++) {
      if (tr_ckpt) FreeParsetree(tr_ckpt);
      if ((status = cm_CheckptTrCYKAlignHB(cm, errbuf, dsq, L, size_limit, pass_idx, &tr_ckpt, &ckpt_mode, &ckpt_sc)) != eslOK)
        cm_Fail("cm_CheckptTrCYKAlignHB: %s", errbuf);
    }
    double t_ckpt = (now_sec() - t0) / reps;

    n_total++;
    t_loop_tot += t_loop; t_ckpt_tot += t_ckpt;
    float dsc = fabsf(loop_sc - ckpt_sc);
    int mode_agree = (loop_mode == ckpt_mode);
    int peq = ParsetreeCompare(tr_loop, tr_ckpt);
    int pass = mode_agree && (dsc <= 0.01) && peq;
    if (mode_agree) n_mode_agree++;
    if (pass) n_pass++;
    printf("# R2 seq=%s L=%d loop_mode=%s ckpt_mode=%s mode_agree=%s sc_loop=%.6f sc_ckpt=%.6f |dsc|=%.3e "
           "nodes_loop=%d nodes_ckpt=%d ParsetreeCompare=%s t_loop=%.4fs t_ckpt=%.4fs speedup=%.2fx VERDICT=%s\n",
           sq->name, L, MarginalMode(loop_mode), MarginalMode(ckpt_mode), mode_agree?"YES":"NO",
           loop_sc, ckpt_sc, dsc, tr_loop->n, tr_ckpt->n, peq?"IDENTICAL":"DIFFER",
           t_loop, t_ckpt, (t_ckpt>0.)?t_loop/t_ckpt:0., pass?"PASS":"FAIL");

    FreeParsetree(tr_loop); FreeParsetree(tr_ckpt);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  printf("# R2 SUMMARY: n=%d mode_agree=%d pass=%d VERDICT=%s  t_loop_tot=%.4fs t_ckpt_tot=%.4fs aggregate_speedup=%.2fx\n",
         n_total, n_mode_agree, n_pass, (n_pass==n_total && n_total>0) ? "PASS" : "FAIL",
         t_loop_tot, t_ckpt_tot, (t_ckpt_tot>0.)?t_loop_tot/t_ckpt_tot:0.);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

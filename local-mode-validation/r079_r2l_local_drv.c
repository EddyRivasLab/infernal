/* r079_r2l_local_drv.c : brief 26_0610-079 (R2-L) validation.
 *
 * Validates cm_CheckptTrCYKAlignHB()'s new EL (local end) + local-begin
 * support against TODAY'S PRODUCTION pass-1: cm_alndata.c's do_trckpt_r4
 * ncand-loop (TrCYKDivideAndConquerHB run once per root-valid mode
 * candidate, argmax) -- run in LOCAL config (the common/default case R2 was
 * previously unable to serve).  Directly descended from
 * r078_r2_ckpttrcyk_drv.c (GLOBAL-only oracle comparison + wall-time), with
 * local-config wiring mirrored from rl5_ckpttrpinalign_drv.c (-g toggle,
 * CM_CONFIG_LOCAL|HMMLOCAL|HMMEL, default local) and EL-node counting
 * (count_el) mirrored from rl1's hbdnc_drv.c convention.
 *
 * Gate: resolved mode MUST match; score agreement within this project's own
 * established D&C-vs-full-fill tolerance (0.01); ParsetreeCompare()==0
 * (byte-identical parse incl. marginal mode per node + every discovered k*
 * + every USED_EL node, since ParsetreeCompare checks tr->state[] too).
 *
 * Requires the engine be compiled with -DCKPT_TRCYK_R2L_ALLOW_LOCAL (the
 * brief-079 fail-fast gate is only lifted under this macro until validation
 * passes and the gate is removed for good).
 *
 * Usage: r079_r2l_local_drv [-g] [--tau x] [--reps N] <cmfile> <seqfile>
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
  { "-g",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "global config (default local)", 0 },
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { "--reps", eslARG_INT,  "1",   NULL, NULL, NULL, NULL, NULL, "wall-time reps per seq", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

/* count USED_EL nodes (parsetree nodes whose state index == cm->M) */
static int
count_el(Parsetree_t *tr, int M)
{
  int n = 0, x;
  for (x = 0; x < tr->n; x++) if (tr->state[x] == M) n++;
  return n;
}

/* report the root-child (local-begin) entry state: tr->state[] of the first
 * node inserted as ROOT_S's LEFT_CHILD (node index 1, mirrors rl5's
 * entry_state() -- node 0 is the [-1,ROOT] sentinel, node 1 is v==0's own
 * traceback attach... actually the FIRST real entry is tr->state[1] since
 * InsertTraceNodewithMode(tr,-1,...) creates node 0 = ROOT_S itself, and the
 * v==0 case attaches child b at node 1). */
static int
entry_state(Parsetree_t *tr)
{
  return (tr->n > 1) ? tr->state[1] : -1;
}

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

/* trusted monolithic (non-checkpointed, non-D&C) banded truncated CYK oracle:
 * cm_TrAlignHB(do_optacc=FALSE, do_sample=FALSE) internally calls the file-static
 * cm_tr_alignT_hb() -> cm_TrCYKInsideAlignHB() (the SAME well-established engine
 * production uses for preset_mode==UNKNOWN CYK) and does its own traceback -> a
 * real Parsetree_t, directly ParsetreeCompare-able against tr_ckpt (unlike the
 * score-only cross-check this was promoted from -- see brief 079 summary: this
 * became the PRIMARY local-mode gate after loop_resolve/TrCYKDivideAndConquerHB
 * was caught reporting an unreachable, too-high score in local config). */
static int
mono_resolve(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int pass_idx,
             char *ret_mode, float *ret_sc, Parsetree_t **ret_tr)
{
  int status;
  CM_TR_HB_MX *mx = cm_tr_hb_mx_Create(cm);
  CM_TR_HB_SHADOW_MX *shmx = cm_tr_hb_shadow_mx_Create(cm);
  Parsetree_t *tr = NULL; char mode = TRMODE_UNKNOWN; float sc = IMPOSSIBLE;
  status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, TRMODE_UNKNOWN, pass_idx,
                         FALSE /*do_optacc*/, FALSE /*do_sample*/, mx, shmx, NULL /*post_mx*/,
                         NULL /*emit_mx*/, NULL /*r*/, NULL /*ret_ppstr*/, &tr, &mode, NULL /*ret_avgpp*/, &sc);
  cm_tr_hb_mx_Destroy(mx); cm_tr_hb_shadow_mx_Destroy(shmx);
  if (status != eslOK) return status;
  *ret_mode = mode; *ret_sc = sc; *ret_tr = tr;
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
    { printf("usage: r079_r2l_local_drv [-g] [--tau x] [--reps N] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
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
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  int M = cm->M;

  int nb = 0, v;
  for (v = 0; v < M; v++) if (cm->sttype[v] == B_st) nb++;
  if (nb < 1) cm_Fail("CM has no bifurcations (bps=0): not an R2-L target");
  printf("# CM=%s M=%d clen=%d bifs=%d mode=%s tau=%g reps=%d\n", cm->name, cm->M, cm->clen, nb,
         do_global ? "global" : "local(full:begins+ends)", tau, reps);
  printf("# flags: LOCAL_BEGIN=%d LOCAL_END=%d\n",
         (cm->flags & CMH_LOCAL_BEGIN)?1:0, (cm->flags & CMH_LOCAL_END)?1:0);

  int n_total = 0, n_pass = 0, n_mode_agree = 0;
  int n_el_loop = 0, n_el_ckpt = 0, n_el_mismatch = 0, n_dnc_mismatch = 0;
  int n_nonroot_entry = 0;
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

    /* ---- trusted oracle: cm_tr_alignT_hb / cm_TrCYKInsideAlignHB (mono) ---- */
    char mono_mode = TRMODE_UNKNOWN; float mono_sc = IMPOSSIBLE; Parsetree_t *tr_mono = NULL;
    if (mono_resolve(cm, errbuf, dsq, L, size_limit, pass_idx, &mono_mode, &mono_sc, &tr_mono) != eslOK)
      cm_Fail("mono_resolve (cm_tr_alignT_hb) failed for seq %s: %s", sq->name, errbuf);

    /* ---- R2-L engine: cm_CheckptTrCYKAlignHB (now local-capable) ---- */
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
    /* PRIMARY gate: ckpt vs trusted mono oracle (NOT vs the D&C ncand-loop --
     * see mono_resolve()'s comment: the D&C loop is informational-only here). */
    float dsc = fabsf(mono_sc - ckpt_sc);
    int mode_agree = (mono_mode == ckpt_mode);
    int peq = ParsetreeCompare(tr_mono, tr_ckpt);
    int pass = mode_agree && (dsc <= 0.01) && peq;
    if (mode_agree) n_mode_agree++;
    if (pass) n_pass++;
    /* informational: does this case ALSO agree with the D&C ncand-loop? */
    int loop_agree = (loop_mode == ckpt_mode) && (fabsf(loop_sc - ckpt_sc) <= 0.01) && ParsetreeCompare(tr_loop, tr_ckpt);
    if (! loop_agree) n_dnc_mismatch++;

    int elA = count_el(tr_mono, M), elB = count_el(tr_ckpt, M);
    n_el_loop += elA; n_el_ckpt += elB;
    if (elA != elB) n_el_mismatch++;
    int entryA = entry_state(tr_mono);
    int entryB = entry_state(tr_ckpt);
    if (entryB != 0) n_nonroot_entry++;

    printf("# R2L seq=%s L=%d mono_mode=%s ckpt_mode=%s mode_agree=%s sc_mono=%.6f sc_ckpt=%.6f |dsc|=%.3e "
           "nodes_mono=%d nodes_ckpt=%d ParsetreeCompare=%s elA=%d elB=%d entryA=%d entryB=%d "
           "t_loop=%.4fs t_ckpt=%.4fs speedup=%.2fx VERDICT=%s dnc_loop_agree=%s sc_dnc=%.6f\n",
           sq->name, L, MarginalMode(mono_mode), MarginalMode(ckpt_mode), mode_agree?"YES":"NO",
           mono_sc, ckpt_sc, dsc, tr_mono->n, tr_ckpt->n, peq?"IDENTICAL":"DIFFER",
           elA, elB, entryA, entryB,
           t_loop, t_ckpt, (t_ckpt>0.)?t_loop/t_ckpt:0., pass?"PASS":"FAIL",
           loop_agree?"YES":"NO", loop_sc);

    FreeParsetree(tr_loop); FreeParsetree(tr_ckpt); FreeParsetree(tr_mono);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  printf("# R2L SUMMARY: n=%d mode_agree=%d pass=%d VERDICT=%s  t_loop_tot=%.4fs t_ckpt_tot=%.4fs aggregate_speedup=%.2fx "
         "el_loop=%d el_ckpt=%d el_mismatch=%d nonroot_entry=%d dnc_loop_mismatch=%d\n",
         n_total, n_mode_agree, n_pass, (n_pass==n_total && n_total>0) ? "PASS" : "FAIL",
         t_loop_tot, t_ckpt_tot, (t_ckpt_tot>0.)?t_loop_tot/t_ckpt_tot:0.,
         n_el_loop, n_el_ckpt, n_el_mismatch, n_nonroot_entry, n_dnc_mismatch);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

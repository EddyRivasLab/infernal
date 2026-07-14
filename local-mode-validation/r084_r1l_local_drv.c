/* r084_r1l_local_drv.c : brief 26_0610-084 (R1-L) validation.
 *
 * Validates cm_CheckptCYKAlignHB()'s NEW EL (local end) + local-begin support
 * (brief 084) against the trusted monolithic banded-CYK oracle:
 * cm_AlignHB(do_optacc=FALSE) -> cm_alignT_hb() -> cm_CYKInsideAlignHB() (the
 * SAME full-matrix banded CYK engine + traceback production uses), run in the
 * SAME config.  A monolithic oracle (not the D&C CYKDivideAndConquerHB) is used
 * deliberately: brief 079's R2-L work caught the truncated D&C oracle reporting
 * an unreachable, too-high score in local config, so this project now prefers
 * the monolithic full-matrix engine for local-mode parse comparison.
 *
 * Non-truncated, structured (bps>0) CMs.  Default LOCAL (begins+ends); -g for
 * GLOBAL regression.  Also usable on bps>0/bifs=0 CMs (single-hairpin, MP/MR
 * but no B) to exercise the do_checkpt_r3 qualifier-gate landmine (brief 084 §3).
 *
 * Gate: score agreement within this project's own established D&C-vs-full-fill
 * tolerance (0.01) AND ParsetreeCompare()==0 (byte-identical parse incl. every
 * discovered k*, every USED_EL node -- ParsetreeCompare checks tr->state[]).
 *
 * Usage: r084_r1l_local_drv [-g] [--tau x] [--reps N] <cmfile> <seqfile>
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

/* first real entry state (node 1); != 0 iff a local begin was used */
static int
entry_state(Parsetree_t *tr)
{
  return (tr->n > 1) ? tr->state[1] : -1;
}

/* trusted monolithic (non-checkpointed, non-D&C) banded CYK oracle:
 * cm_AlignHB(do_optacc=FALSE) -> cm_alignT_hb() -> cm_CYKInsideAlignHB(). */
static int
mono_resolve(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,
             float *ret_sc, Parsetree_t **ret_tr)
{
  int status;
  CM_HB_MX          *mx   = cm_hb_mx_Create(cm->M);
  CM_HB_SHADOW_MX   *shmx = cm_hb_shadow_mx_Create(cm);
  Parsetree_t *tr = NULL; float sc = IMPOSSIBLE;
  status = cm_AlignHB(cm, errbuf, dsq, L, size_limit, FALSE /*do_optacc*/, FALSE /*do_sample*/,
                       mx, shmx, NULL /*post_mx*/, NULL /*emit_mx*/, NULL /*r*/,
                       NULL /*ret_ppstr*/, &tr, NULL /*ret_avgpp*/, &sc);
  cm_hb_mx_Destroy(mx); cm_hb_shadow_mx_Destroy(shmx);
  if (status != eslOK) return status;
  *ret_sc = sc; *ret_tr = tr;
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
    { printf("usage: r084_r1l_local_drv [-g] [--tau x] [--reps N] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  int    reps      = esl_opt_GetInteger(go, "--reps");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->align_opts  |= CM_ALIGN_HBANDED;
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  int M = cm->M;

  int nb = 0, nbps = 0, v;
  for (v = 0; v < M; v++) { if (cm->sttype[v] == B_st) nb++;
                            if (cm->sttype[v]==B_st||cm->sttype[v]==MP_st||cm->sttype[v]==MR_st) nbps++; }
  printf("# CM=%s M=%d clen=%d bifs=%d bps_states=%d mode=%s tau=%g reps=%d\n", cm->name, cm->M, cm->clen, nb, nbps,
         do_global ? "global" : "local(full:begins+ends)", tau, reps);
  printf("# flags: LOCAL_BEGIN=%d LOCAL_END=%d\n",
         (cm->flags & CMH_LOCAL_BEGIN)?1:0, (cm->flags & CMH_LOCAL_END)?1:0);

  int n_total = 0, n_pass = 0, n_el_o = 0, n_el_c = 0, n_el_mismatch = 0, n_nonroot_entry = 0;
  double t_o_tot = 0., t_c_tot = 0.;

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- oracle: cm_AlignHB(do_optacc=FALSE) monolithic banded CYK ---- */
    float sc_o = IMPOSSIBLE; Parsetree_t *tr_o = NULL;
    double t0 = now_sec(); int r;
    for (r = 0; r < reps; r++) {
      if (tr_o) FreeParsetree(tr_o);
      if (mono_resolve(cm, errbuf, dsq, L, size_limit, &sc_o, &tr_o) != eslOK)
        cm_Fail("mono_resolve (cm_AlignHB) failed for seq %s: %s", sq->name, errbuf);
    }
    double t_o = (now_sec() - t0) / reps;

    /* ---- R1-L engine: cm_CheckptCYKAlignHB (now local-capable) ---- */
    Parsetree_t *tr_c = NULL; float sc_c = 0.;
    t0 = now_sec();
    for (r = 0; r < reps; r++) {
      if (tr_c) FreeParsetree(tr_c);
      if ((status = cm_CheckptCYKAlignHB(cm, errbuf, dsq, L, size_limit, &tr_c, &sc_c)) != eslOK)
        cm_Fail("cm_CheckptCYKAlignHB: %s", errbuf);
    }
    double t_c = (now_sec() - t0) / reps;

    n_total++;
    t_o_tot += t_o; t_c_tot += t_c;
    float dsc = fabsf(sc_o - sc_c);
    int peq = ParsetreeCompare(tr_o, tr_c);
    int pass = (dsc <= 0.01) && peq;
    if (pass) n_pass++;
    int elA = count_el(tr_o, M), elB = count_el(tr_c, M);
    n_el_o += elA; n_el_c += elB;
    if (elA != elB) n_el_mismatch++;
    int entryB = entry_state(tr_c);
    if (entryB != 0) n_nonroot_entry++;

    printf("# R1L seq=%s L=%d sc_oracle=%.6f sc_ckpt=%.6f |dsc|=%.3e nodes_o=%d nodes_c=%d "
           "ParsetreeCompare=%s elA=%d elB=%d entryA=%d entryB=%d t_o=%.4fs t_c=%.4fs speedup=%.2fx VERDICT=%s\n",
           sq->name, L, sc_o, sc_c, dsc, tr_o->n, tr_c->n, peq?"IDENTICAL":"DIFFER",
           elA, elB, entry_state(tr_o), entryB, t_o, t_c, (t_c>0.)?t_o/t_c:0., pass?"PASS":"FAIL");

    FreeParsetree(tr_o); FreeParsetree(tr_c);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  printf("# R1L SUMMARY: n=%d pass=%d VERDICT=%s  t_oracle_tot=%.4fs t_ckpt_tot=%.4fs aggregate_speedup=%.2fx "
         "el_oracle=%d el_ckpt=%d el_mismatch=%d nonroot_entry=%d\n",
         n_total, n_pass, (n_pass==n_total && n_total>0) ? "PASS" : "FAIL",
         t_o_tot, t_c_tot, (t_c_tot>0.)?t_o_tot/t_c_tot:0.,
         n_el_o, n_el_c, n_el_mismatch, n_nonroot_entry);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

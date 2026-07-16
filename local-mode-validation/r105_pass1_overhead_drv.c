/* r105_pass1_overhead_drv.c : brief 26_0610-105 Task B.
 *
 * Measures rung-3's pass-1 CYK cost as a fraction of end-to-end time on
 * bp>0/bif==0 CMs (a simple hairpin/stem-loop, no multifurcating junctions),
 * where 26_0430's code-reading finding says pass-1's full parse is discarded
 * unused (rung3_kpin_from_cyk only extracts pins from B_st nodes, which don't
 * exist here -- kpin[] stays all -1, and pass-2 runs as an ordinary unpinned
 * sweep regardless of what pass-1 computed).
 *
 * For each seq, times:
 *   t_cyk   : pass-1 alone, cm_CheckptCYKAlignHB (R1, the checkpointed CYK
 *             that brief 084 wired in as the current default pass-1 engine)
 *   t_dnc   : pass-1 alone, the D&C alternative CYKDivideAndConquerHB
 *             (INFERNAL_CKPT_FORCE_DNC path)
 *   t_pass2 : pass-2 alone (cm_CheckptPostAlignHB + cm_CheckptOptAccAlignHB),
 *             fed an all-(-1) kpin[] -- exactly what pass-2 receives today
 *             for a bif==0 CM regardless of which pass-1 engine ran.  This
 *             also serves as the ckpt_old-equivalent cost proxy (Task A found
 *             ckpt_old cannot safely handle MP_st/MR_st internally, so it
 *             can't be measured directly; pass-2 alone is structurally the
 *             same "single unpinned checkpointed sweep" shape ckpt_old uses).
 *   t_total : t_cyk + t_pass2 (what today's dispatch actually pays, R1 pass-1)
 *
 * Usage: r105_pass1_overhead_drv [-g] [--tau x] [--reps N] <cmfile> <seqfile>
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
#include "esl_vectorops.h"

#include "hmmer.h"
#include "infernal.h"

static ESL_OPTIONS options[] = {
  { "-h",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "help", 0 },
  { "-g",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "global config (default local)", 0 },
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { "--reps", eslARG_INT,  "3",   NULL, NULL, NULL, NULL, NULL, "wall-time reps per seq", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

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
    { printf("usage: r105_pass1_overhead_drv [-g] [--tau x] [--reps N] <cmfile> <seqfile>\n"); exit(0); }
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
  if (nb != 0) { printf("# WARNING: bifs=%d != 0, this driver assumes bif==0 (kpin all -1)\n", nb); }

  int *kpin = malloc(sizeof(int) * M);
  for (v = 0; v < M; v++) kpin[v] = -1;

  CM_HB_EMIT_MX *emit_mx = cm_hb_emit_mx_Create(cm);

  double t_cyk_tot = 0., t_dnc_tot = 0., t_pass2_tot = 0.;
  int n_total = 0;

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- pass-1, R1 (cm_CheckptCYKAlignHB, today's default) ---- */
    Parsetree_t *tr_cyk = NULL; float sc_cyk = 0.;
    double t0 = now_sec(); int r;
    for (r = 0; r < reps; r++) {
      if (tr_cyk) FreeParsetree(tr_cyk);
      if ((status = cm_CheckptCYKAlignHB(cm, errbuf, dsq, L, size_limit, &tr_cyk, &sc_cyk)) != eslOK)
        cm_Fail("cm_CheckptCYKAlignHB: %s", errbuf);
    }
    double t_cyk = (now_sec() - t0) / reps;
    if (tr_cyk) { FreeParsetree(tr_cyk); tr_cyk = NULL; }

    /* ---- pass-1, D&C alternative (CYKDivideAndConquerHB) ---- */
    double t0d = now_sec();
    Parsetree_t *tr_dnc = NULL;
    for (r = 0; r < reps; r++) {
      if (tr_dnc) FreeParsetree(tr_dnc);
      tr_dnc = NULL;
      (void) CYKDivideAndConquerHB(cm, dsq, L, 0, 1, L, &tr_dnc, cm->cp9b);
    }
    double t_dnc = (now_sec() - t0d) / reps;
    if (tr_dnc) { FreeParsetree(tr_dnc); tr_dnc = NULL; }

    /* ---- pass-2 alone, all-(-1) kpin (what pass-2 gets today for bif==0;
     *      also the ckpt_old-equivalent proxy per Task A) ---- */
    double t0p = now_sec();
    for (r = 0; r < reps; r++) {
      float r3_Z = 0.;
      esl_vec_FSet(emit_mx->l_pp_mem, emit_mx->l_ncells_valid, IMPOSSIBLE);
      esl_vec_FSet(emit_mx->r_pp_mem, emit_mx->r_ncells_valid, IMPOSSIBLE);
      if ((status = cm_CheckptPostAlignHB(cm, errbuf, dsq, L, size_limit, emit_mx, kpin, &r3_Z)) != eslOK)
        cm_Fail("cm_CheckptPostAlignHB: %s", errbuf);
      Parsetree_t *tr2 = NULL; float pp2 = 0.;
      if ((status = cm_CheckptOptAccAlignHB(cm, errbuf, dsq, L, size_limit, emit_mx, kpin, NULL, &tr2, NULL, &pp2)) != eslOK)
        cm_Fail("cm_CheckptOptAccAlignHB: %s", errbuf);
      if (tr2) FreeParsetree(tr2);
    }
    double t_pass2 = (now_sec() - t0p) / reps;

    n_total++;
    t_cyk_tot += t_cyk; t_dnc_tot += t_dnc; t_pass2_tot += t_pass2;
    double t_total_r1  = t_cyk + t_pass2;
    double t_total_dnc = t_dnc + t_pass2;
    printf("# R105 seq=%s L=%d t_cyk(R1)=%.5fs t_dnc=%.5fs t_pass2=%.5fs t_total(R1+p2)=%.5fs "
           "cyk_frac=%.4f dnc_frac=%.4f\n",
           sq->name, L, t_cyk, t_dnc, t_pass2, t_total_r1,
           (t_total_r1>0.)?t_cyk/t_total_r1:0., (t_total_dnc>0.)?t_dnc/t_total_dnc:0.);

    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  double t_total_tot = t_cyk_tot + t_pass2_tot;
  printf("# R105 SUMMARY: n=%d t_cyk_tot=%.5fs t_dnc_tot=%.5fs t_pass2_tot=%.5fs t_total_tot(R1+p2)=%.5fs "
         "cyk_frac=%.4f dnc_frac=%.4f\n",
         n_total, t_cyk_tot, t_dnc_tot, t_pass2_tot, t_total_tot,
         (t_total_tot>0.)?t_cyk_tot/t_total_tot:0.,
         ((t_dnc_tot+t_pass2_tot)>0.)?t_dnc_tot/(t_dnc_tot+t_pass2_tot):0.);
  free(kpin);
  cm_hb_emit_mx_Destroy(emit_mx);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

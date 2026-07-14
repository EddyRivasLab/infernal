/* r090_forked_drv.c : brief 26_0610-090 crash-tolerant byte-exact panel driver.
 *
 * Same gate as r086_argmax_drv (mode match AND |dev|<=0.01 AND ParsetreeCompare==0),
 * but FORKS a child per sequence so the known pre-existing stock segfault (067/068,
 * ~18% on very short truncated-local seqs) on one seq does not abort the whole panel.
 * Parent prints one line per seq: PASS / FAIL / CRASH. This lets a pre-fix vs post-fix
 * fail-set diff be computed reliably despite the crashes.
 *
 * The CM is configured once in the parent; each child inherits it via fork (COW) and
 * does the per-seq band derivation + DP + oracle + compare, then _exit(0/1).
 *
 * Usage: r090_forked_drv [-g] [--tau x] <cmfile> <seqfile>
 */
#include <esl_config.h>
#include <p7_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
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

/* child: process one seq, print result line, _exit(0). On crash the parent sees the signal. */
static void
child_process(CM_t *cm, ESL_ALPHABET *abc, ESL_SQ *sq, float size_limit, int pass_idx)
{
  char errbuf[eslERRBUFSIZE];
  int L = sq->n;
  ESL_DSQ *dsq = sq->dsq;
  if (cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0) != eslOK)
    { printf("BANDFAIL seq=%s L=%d\n", sq->name, L); fflush(stdout); _exit(0); }

  char loop_mode = TRMODE_UNKNOWN; float loop_sc = IMPOSSIBLE; Parsetree_t *tr_loop = NULL;
  if (loop_resolve(cm, dsq, L, pass_idx, &loop_mode, &loop_sc, &tr_loop) != eslOK)
    { printf("LOOPFAIL seq=%s L=%d\n", sq->name, L); fflush(stdout); _exit(0); }

  CM_TR_HB_MX *mx = cm_tr_hb_mx_Create(cm);
  CM_TR_HB_SHADOW_MX *shmx = cm_tr_hb_shadow_mx_Create(cm);
  Parsetree_t *tr_or = NULL; char or_mode = TRMODE_UNKNOWN; float sc_or = IMPOSSIBLE;
  int status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, TRMODE_UNKNOWN, pass_idx, FALSE, FALSE,
                            mx, shmx, NULL, NULL, NULL, NULL, &tr_or, &or_mode, NULL, &sc_or);
  if (status != eslOK) { printf("ORACLEFAIL seq=%s L=%d: %s\n", sq->name, L, errbuf); fflush(stdout); _exit(0); }

  float dev = fabsf(loop_sc - sc_or);
  int peq = ParsetreeCompare(tr_or, tr_loop);
  int ok = (loop_mode==or_mode) && (dev <= 0.01) && peq;
  printf("%s seq=%s L=%d loop_mode=%s or_mode=%s sc_loop=%.6f sc_oracle=%.6f dev=%.3e ParsetreeCompare=%s nodes_loop=%d nodes_or=%d\n",
         ok?"PASS":"FAIL", sq->name, L, MarginalMode(loop_mode), MarginalMode(or_mode), loop_sc, sc_or, dev,
         peq?"IDENTICAL":"DIFFER", tr_loop?tr_loop->n:-1, tr_or?tr_or->n:-1);
  fflush(stdout);
  _exit(0);
}

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: r090_forked_drv [-g] [--tau x] <cmfile> <seqfile>\n"); exit(0); }
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
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();

  int nb = 0, v;
  for (v = 0; v < cm->M; v++) if (cm->sttype[v] == B_st) nb++;
  printf("# CM=%s M=%d bifs=%d mode=%s tau=%g LOCAL_END=%d\n", cm->name, cm->M, nb, do_global?"global":"local", tau, (cm->flags&CMH_LOCAL_END)?1:0);
  fflush(stdout);

  int n=0, npass=0, nfail=0, ncrash=0, nother=0;
  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    n++;
    fflush(stdout);
    pid_t pid = fork();
    if (pid == 0) { child_process(cm, abc, sq, size_limit, pass_idx); /* never returns */ }
    else if (pid < 0) cm_Fail("fork failed");
    int wstatus;
    waitpid(pid, &wstatus, 0);
    if (WIFSIGNALED(wstatus)) { ncrash++; printf("CRASH seq=%s L=%d signal=%d\n", sq->name, (int)sq->n, WTERMSIG(wstatus)); fflush(stdout); }
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  printf("# DONE n=%d (per-seq PASS/FAIL/CRASH lines above; crashes=%d)\n", n, ncrash);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

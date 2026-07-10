/* r078_r1_ckptcyk_drv.c : brief 26_0610-078 R1 validation.
 *
 * Validates cm_CheckptCYKAlignHB() -- the new checkpointed, bifurcation
 * (k*)-DISCOVERING CYK max-DP engine (cm_dpalign.c, brief 078 R1) -- against
 * CYKDivideAndConquerHB() (the existing D&C CYK oracle, cm_dpsmall.c:497-561,
 * used elsewhere in this project as a valid reference regardless of whether
 * it's in cmalign's current production dispatch -- see brief 078 citations).
 *
 * GLOBAL, non-truncated, structured (bps>0) CMs only (R1 scope).
 *
 * Gate: score agreement within this project's own established D&C-vs-full-
 * fill tolerance (0.01, hbdnc_drv.c's convention -- NOT literal bitwise
 * equality, since D&C recombination and a single-sweep fill sum terms in
 * different orders) AND ParsetreeCompare()==0 (byte-identical parse,
 * including every discovered k* -- B nodes' TRACE_RIGHT_CHILD i/j and tree
 * structure are part of the node-for-node comparison, per rl4/rl5's
 * convention).
 *
 * Usage: r078_r1_ckptcyk_drv [--tau x] <cmfile> <seqfile>
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
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
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
    { printf("usage: r078_r1_ckptcyk_drv [--tau x] <cmfile> <seqfile>\n"); exit(0); }
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->align_opts  |= CM_ALIGN_HBANDED;
  /* R1 scope: GLOBAL only -- do NOT set CM_CONFIG_LOCAL/HMMLOCAL/HMMEL */
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  int M = cm->M;

  int nb = 0, v;
  for (v = 0; v < M; v++) if (cm->sttype[v] == B_st) nb++;
  if (nb < 1) cm_Fail("CM has no bifurcations (bps=0): not an R1 target");
  printf("# CM=%s M=%d clen=%d bifs=%d mode=global tau=%g\n", cm->name, cm->M, cm->clen, nb, tau);

  int n_total = 0, n_pass = 0;

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- oracle: CYKDivideAndConquerHB ---- */
    Parsetree_t *tr_o = NULL;
    float sc_o = CYKDivideAndConquerHB(cm, dsq, L, 0, 1, L, &tr_o, cm->cp9b);

    /* ---- R1 engine: cm_CheckptCYKAlignHB ---- */
    Parsetree_t *tr_c = NULL; float sc_c = 0.;
    if ((status = cm_CheckptCYKAlignHB(cm, errbuf, dsq, L, size_limit, &tr_c, &sc_c)) != eslOK)
      cm_Fail("cm_CheckptCYKAlignHB: %s", errbuf);

    n_total++;
    float dsc = fabsf(sc_o - sc_c);
    int peq = ParsetreeCompare(tr_o, tr_c);
    int pass = (dsc <= 0.01) && peq;
    if (pass) n_pass++;
    printf("# R1 seq=%s L=%d sc_oracle=%.6f sc_ckpt=%.6f |dsc|=%.3e nodes_oracle=%d nodes_ckpt=%d ParsetreeCompare=%s VERDICT=%s\n",
           sq->name, L, sc_o, sc_c, dsc, tr_o->n, tr_c->n, peq?"IDENTICAL":"DIFFER", pass?"PASS":"FAIL");

    FreeParsetree(tr_o); FreeParsetree(tr_c);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  printf("# R1 SUMMARY: n=%d pass=%d VERDICT=%s\n", n_total, n_pass, (n_pass==n_total && n_total>0) ? "PASS" : "FAIL");
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

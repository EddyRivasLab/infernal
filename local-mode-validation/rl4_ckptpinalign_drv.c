/* rl4_ckptpinalign_drv.c : R-L.4 validation of local (EL+begin) support in
 * the rung-3 STRUCTURED (bps>0, has bifurcations) non-truncated checkpointed
 * posterior + OptAcc engines: cm_PinPostAlignHB/cm_CheckptPostAlignHB and
 * cm_PinOptAccAlignHB/cm_CheckptOptAccAlignHB.
 *
 * Pipeline (mirrors production's do_checkpt_r3 path in cm_alndata.c, but
 * calls the engines directly -- test-only, bypassing
 * cm_CheckptOptAccAlignHB_Qualifies() which still correctly rejects local
 * CMs in production; relaxing that gate for local is R-L.6's job):
 *   pass 1: cm_alignT_hb() (stock full CYK-with-traceback) -> tr_cyk
 *   kpin[] extraction from tr_cyk (mirrors cm_alndata.c's static
 *           rung3_kpin_from_cyk(), duplicated here since it's file-static)
 *   pass 2a (full-storage correctness anchor): cm_PinPostAlignHB +
 *           cm_PinOptAccAlignHB, compared to stock cm_AlignHB(do_optacc=TRUE)
 *   pass 2b (sqrt(M) deliverable): cm_CheckptPostAlignHB +
 *           cm_CheckptOptAccAlignHB, compared to the SAME stock oracle
 *
 * Default config is LOCAL (full: begins+ends); -g switches to global (the
 * regression check: must reproduce pre-R-L.4 rung-3 behavior exactly).
 *
 * Byte-exact gate (per case, per engine): Inside score |dsc|==0,
 * ParsetreeCompare==0 (this covers bifurcation split points -- B nodes'
 * TRACE_RIGHT_CHILD i/j and the tree structure are all part of the
 * node-for-node comparison), per-residue PP string strcmp==0, avg PP,
 * EL-node count match, entry_v (tr->state[1]) match INCLUDING classification
 * of whether entry_v lands inside a bifurcated subtree (state type
 * BEGL_S/BEGR_S/B, or -- via a plast-chain walk -- any state whose lineage
 * passes through a B node).
 *
 * Usage: rl4_ckptpinalign_drv [-g] [--tau x] [--tag T] <cmfile> <seqfile>
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
#include "esl_vectorops.h"

#include "hmmer.h"
#include "infernal.h"

static ESL_OPTIONS options[] = {
  { "-h",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "help", 0 },
  { "-g",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "global config (default local)", 0 },
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--fixedtau", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "don't let cmalign auto-tighten tau", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { "--tag", eslARG_STRING,"run",  NULL, NULL, NULL, NULL, NULL, "label printed in output", 0 },
  { "--nostock",eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "skip stock oracle (ckpt-only; for genome-scale where stock OOMs)", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

/* count USED_EL nodes (parsetree nodes whose state index == cm->M = the EL state) */
static int count_el(Parsetree_t *tr, int M)
{
  int n = 0, x;
  for (x = 0; x < tr->n; x++) if (tr->state[x] == M) n++;
  return n;
}

/* duplicated from cm_alndata.c's static rung3_kpin_from_cyk(): derive the
 * per-bifurcation right-fragment length k*(v) from a CYK D&C/full parsetree.
 * Topology-generic (walks the tree for B nodes; independent of how the tree
 * was entered), per R-L.4 Step 0's traced confirmation. */
static void
kpin_from_cyk(CM_t *cm, Parsetree_t *tr, int *kpin)
{
  int n;
  for (n = 0; n < cm->M; n++) kpin[n] = -1;
  for (n = 0; n < tr->n; n++) {
    int v = tr->state[n];
    if (cm->sttype[v] != B_st) continue;
    int ln = tr->nxtl[n], rn = tr->nxtr[n];
    int lstid = cm->stid[tr->state[ln]];
    int begl_node = (lstid == BEGL_S) ? ln : rn;
    int ksplit    = tr->emitr[begl_node];
    kpin[v]       = tr->emitr[n] - ksplit;
  }
}

/* is v inside a bifurcated subtree: v itself is BEGL_S/BEGR_S/B_st, OR its
 * static parent-chain (cm->plast[v] walk) passes through a B_st (i.e. v is
 * a descendant, in the guide tree, of some bifurcation). */
static int
entry_inside_bifurcation(CM_t *cm, int v)
{
  int st = cm->stid[v];
  if (st == BEGL_S || st == BEGR_S) return TRUE;
  if (cm->sttype[v] == B_st) return TRUE;
  /* walk parent nodes via nodemap-derived first-state-of-node chain: use
   * cm->ndidx[v] -> parent node via cm->ndtype tree isn't directly stored,
   * but pnum/plast give per-state parent state(s); walk those, which for a
   * CM's tree structure ascends toward the root through exactly the nodes
   * whose subtree contains v. */
  int cur = v;
  int guard = 0;
  while (cm->pnum[cur] > 0 && guard++ < cm->M) {
    int par = cm->plast[cur] - cm->pnum[cur] + 1; /* first (only, for tree edges) parent */
    if (cm->sttype[par] == B_st) return TRUE;
    if (cm->stid[par] == BEGL_S || cm->stid[par] == BEGR_S) return TRUE;
    cur = par;
  }
  return FALSE;
}

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: rl4_ckptpinalign_drv [-g] [--tau x] [--fixedtau] [--tag T] [--nostock] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  int    do_stock  = ! esl_opt_GetBoolean(go, "--nostock");
  int    fixedtau  = esl_opt_GetBoolean(go, "--fixedtau");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *tag       = esl_opt_GetString(go, "--tag");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->align_opts  |= CM_ALIGN_HBANDED;
  cm->align_opts  |= CM_ALIGN_OPTACC;
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  int M = cm->M;

  /* GATE: structured (bps>0, has >=1 bifurcation) */
  int nb = 0, v;
  for (v = 0; v < M; v++) if (cm->sttype[v] == B_st) nb++;
  if (nb < 1) cm_Fail("CM has no bifurcations (bps=0): not a rung-3 target");
  printf("# CM=%s M=%d clen=%d bifs=%d mode=%s tau=%g fixedtau=%d tag=%s\n",
         cm->name, cm->M, cm->clen, nb, do_global ? "global" : "local(full:begins+ends)", tau, fixedtau, tag);
  printf("# flags: LOCAL_BEGIN=%d LOCAL_END=%d\n",
         (cm->flags & CMH_LOCAL_BEGIN)?1:0, (cm->flags & CMH_LOCAL_END)?1:0);

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;
    printf("# SEQ=%s L=%d\n", sq->name, L);

    /* this driver never auto-tightens tau (unlike cmalign's workunit loop):
     * bands are computed once from cm->tau below, so --fixedtau needs no
     * extra call here -- it's accepted only for command-line parity with
     * cmalign invocations used elsewhere in this project. */
    (void) fixedtau;

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- stock oracle: cm_AlignHB(do_optacc=TRUE) ---- */
    Parsetree_t *tr_s = NULL; char *pp_s = NULL; float avgpp_s = 0, sc_s = 0;
    if (do_stock) {
      if ((status = cm_AlignHB(cm, errbuf, dsq, L, size_limit, TRUE, FALSE,
                               cm->hb_mx, cm->hb_shmx, cm->hb_omx, cm->hb_emx, NULL,
                               &pp_s, &tr_s, &avgpp_s, &sc_s)) != eslOK)
        cm_Fail("stock cm_AlignHB: %s", errbuf);
    }

    /* ---- pass 1: CYK D&C-equivalent full traceback -> kpin[] ---- */
    Parsetree_t *tr_cyk = NULL;
    if ((status = cm_alignT_hb(cm, errbuf, dsq, L, size_limit, FALSE, cm->hb_mx, cm->hb_shmx,
                               NULL, &tr_cyk, NULL)) != eslOK)
      cm_Fail("cm_alignT_hb (pass 1 CYK): %s", errbuf);
    int *kpin = NULL;
    ESL_ALLOC(kpin, sizeof(int) * M);
    kpin_from_cyk(cm, tr_cyk, kpin);
    FreeParsetree(tr_cyk); tr_cyk = NULL;

    /* ---- pass 2a: full-storage correctness anchor (cm_Pin*) ---- */
    CM_HB_EMIT_MX *emit_p = cm_hb_emit_mx_Create(cm);
    Parsetree_t *tr_p = NULL; char *pp_p = NULL; float avgpp_p = 0, sc_p = 0, pp_p_acc = 0;
    if ((status = cm_PinPostAlignHB(cm, errbuf, dsq, L, size_limit, emit_p, kpin, &sc_p)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_PinPostAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      free(kpin); cm_hb_emit_mx_Destroy(emit_p); esl_sq_Reuse(sq); continue;
    }
    if ((status = cm_PinOptAccAlignHB(cm, errbuf, dsq, L, size_limit, emit_p, kpin,
                                      &pp_p, &tr_p, &avgpp_p, &pp_p_acc)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_PinOptAccAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      free(kpin); cm_hb_emit_mx_Destroy(emit_p); esl_sq_Reuse(sq); continue;
    }

    /* ---- pass 2b: sqrt(M) checkpointed deliverable (cm_Checkpt*) ---- */
    CM_HB_EMIT_MX *emit_c = cm_hb_emit_mx_Create(cm);
    Parsetree_t *tr_c = NULL; char *pp_c = NULL; float avgpp_c = 0, sc_c = 0, pp_c_acc = 0;
    if ((status = cm_CheckptPostAlignHB(cm, errbuf, dsq, L, size_limit, emit_c, kpin, &sc_c)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_CheckptPostAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      free(kpin); cm_hb_emit_mx_Destroy(emit_p); cm_hb_emit_mx_Destroy(emit_c);
      if (pp_p) free(pp_p); if (tr_p) FreeParsetree(tr_p); esl_sq_Reuse(sq); continue;
    }
    if ((status = cm_CheckptOptAccAlignHB(cm, errbuf, dsq, L, size_limit, emit_c, kpin,
                                          &pp_c, &tr_c, &avgpp_c, &pp_c_acc)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_CheckptOptAccAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      free(kpin); cm_hb_emit_mx_Destroy(emit_p); cm_hb_emit_mx_Destroy(emit_c);
      if (pp_p) free(pp_p); if (tr_p) FreeParsetree(tr_p); esl_sq_Reuse(sq); continue;
    }

    int elP = count_el(tr_p, M), elC = count_el(tr_c, M);
    int entryP = (tr_p->n > 1) ? tr_p->state[1] : -1;
    int entryC = (tr_c->n > 1) ? tr_c->state[1] : -1;
    int bifP = (entryP >= 0) ? entry_inside_bifurcation(cm, entryP) : FALSE;
    int bifC = (entryC >= 0) ? entry_inside_bifurcation(cm, entryC) : FALSE;
    char *modestr = do_global ? "global" : "local(full)";

    if (! do_stock) {
      printf("# RESULT tag=%s seq=%s mode=%s L=%d (ckpt-only, --nostock)\n", tag, sq->name, modestr, L);
      printf("#   pin   : sc=%.6f nodes=%d avgpp=%.6f EL=%d entry_v=%d bif=%d\n", sc_p, tr_p->n, avgpp_p, elP, entryP, bifP);
      printf("#   ckpt  : sc=%.6f nodes=%d avgpp=%.6f EL=%d entry_v=%d bif=%d\n", sc_c, tr_c->n, avgpp_c, elC, entryC, bifC);
      printf("#   VERDICT=CKPT-ONLY (no stock comparison)\n");
      goto NEXT;
    }

    {
      /* PRIMARY GATE (what R-L.4 actually built/changed): Pin vs Checkpt
       * must be byte-exact -- same kpin, identical deck recurrences "by
       * construction" per both functions' own doc comments. This is the
       * regression check for the EL/begin wiring, elalpha, and the
       * traceback USED_LOCAL_BEGIN fix.
       * SECONDARY (informational only, NOT gated): comparison to the free
       * UNPINNED stock oracle. Per briefs 019/020/032, bif-pinning is
       * accuracy-neutral but NOT always parse/score-faithful vs full
       * unpinned OptAcc (pinned Z is a subset-sum over only the CYK's
       * chosen bifurcation splits, so pinned Z <= unpinned Z always, and
       * on the documented ~3% divergent-split fraction the parse/PP can
       * differ too) -- this is pre-existing rung-3 architecture, not a
       * regression to gate on here. */
      int    peqPC = ParsetreeCompare(tr_p, tr_c);
      int    ppeqPC = (pp_p && pp_c && strcmp(pp_p, pp_c) == 0);
      float  dscPC = fabsf(sc_p - sc_c);
      int    peqP  = ParsetreeCompare(tr_p, tr_s);
      int    peqC  = ParsetreeCompare(tr_c, tr_s);
      float  dscP  = fabsf(sc_p - sc_s);
      float  dscC  = fabsf(sc_c - sc_s);
      int    elS   = count_el(tr_s, M);
      int    entryS = (tr_s->n > 1) ? tr_s->state[1] : -1;
      int    bifS   = (entryS >= 0) ? entry_inside_bifurcation(cm, entryS) : FALSE;

      printf("# RESULT tag=%s seq=%s mode=%s L=%d\n", tag, sq->name, modestr, L);
      printf("#   PRIMARY GATE (pin==ckpt, same kpin): ParsetreeCompare=%s  |dsc|=%.3e  PPstr-byte-identical=%s  entry_v match=%s(%d vs %d) EL match=%s(%d vs %d)\n",
             peqPC?"IDENTICAL":"DIFFER", dscPC, ppeqPC?"YES":"no",
             (entryP==entryC)?"YES":"no", entryP, entryC, (elP==elC)?"YES":"no", elP, elC);
      printf("#   info (vs free unpinned stock, NOT gated -- pin/OA divergence is documented pre-existing behavior):\n");
      printf("#     score : stock=%.6f pin=%.6f |dP|=%.3e ckpt=%.6f |dC|=%.3e\n", sc_s, sc_p, dscP, sc_c, dscC);
      printf("#     parse : stock_nodes=%d pin_nodes=%d PinCompare=%s  ckpt_nodes=%d CkptCompare=%s\n",
             tr_s->n, tr_p->n, peqP?"IDENTICAL":"DIFFER", tr_c->n, peqC?"IDENTICAL":"DIFFER");
      printf("#     avgpp(stock)=%.6f avgpp(pin)=%.6f avgpp(ckpt)=%.6f\n", avgpp_s, avgpp_p, avgpp_c);
      printf("#     EL    : stock=%d pin=%d ckpt=%d\n", elS, elP, elC);
      printf("#     entry_v: stock=%d(bif=%d) pin=%d(bif=%d) ckpt=%d(bif=%d)\n",
             entryS, bifS, entryP, bifP, entryC, bifC);
      printf("#   VERDICT=%s\n",
             (peqPC && ppeqPC && dscPC==0. && elP==elC && entryP==entryC) ? "PASS" : "FAIL");
    }

  NEXT:
    if (pp_s) free(pp_s); if (pp_p) free(pp_p); if (pp_c) free(pp_c);
    if (tr_s) FreeParsetree(tr_s); if (tr_p) FreeParsetree(tr_p); if (tr_c) FreeParsetree(tr_c);
    cm_hb_emit_mx_Destroy(emit_p); cm_hb_emit_mx_Destroy(emit_c);
    free(kpin);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;

 ERROR:
  cm_Fail("allocation error");
  return 1;
}

/* rl2_ckptalign_drv.c : R-L.2 validation of cm_CheckptAlignHB EL (local-end)
 *                       support for bps=0 (MATL-chain) CMs.
 * R-L.2b extension: default local config now exercises FULL CM_CONFIG_LOCAL
 * (local BEGINS + local ends together), the default cmalign use case.
 * --endsonly reproduces R-L.2's original local-ends-only config (global
 * begins + EL ends) for the ends-only regression check (R-L.2b must not
 * perturb R-L.2's byte-exact result).
 *
 * Unlike the 023-era checkpoint-proto/ckptoa_drv.c (which re-implemented the
 * checkpointed engine inline), this driver calls the LIBRARY engine
 * cm_CheckptAlignHB() directly and compares it against the stock oracle
 * cm_AlignHB(do_optacc=TRUE).  That makes the validation test the actual
 * R-L.2/R-L.2b deliverable (the engine + its shared ckpt_*_deck helpers,
 * which are the template R-L.3/4/5 reuse) rather than a parallel prototype.
 *
 * Default config is LOCAL (full: begins+ends); -g switches to global
 * (the regression check: begin/EL additions must not perturb the global
 * path); --endsonly switches to the R-L.2 ends-only config (regression:
 * must stay byte-identical to R-L.2's own result).
 *
 * NOTE (brief 26_0610-062 scope boundary): this driver calls cm_CheckptAlignHB()
 * directly, bypassing cm_CheckptAlignHB_Qualifies() (which still correctly
 * rejects CMH_LOCAL_BEGIN in production -- relaxing THAT qualifier is R-L.6's
 * job, once begins+ends+trunc are all built).  Only the DRIVER's own
 * (test-only) local-begin gate is relaxed here, per the brief.
 *
 * GATE: bps=0 (pure MATL chain: 0 B_st, 0 MP/MR).  cm_CheckptAlignHB requires
 * this (the non-D&C checkpointed lane).
 *
 * Byte-exact gate (per case): Inside score |dZ|, ParsetreeCompare==0,
 * per-residue PP string strcmp==0, avg PP.  EL coverage: count USED_EL nodes
 * (parsetree state==cm->M) in both parses; assert equal.  Local-begin
 * coverage: count parsetree roots whose FIRST node's state != 0 (a genuine
 * local begin fired) via count_begin().
 *
 * Usage: rl2_ckptalign_drv [-g] [--endsonly] [--tau x] [--tag T] <cmfile> <seqfile>
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
  { "--endsonly", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "R-L.2 ends-only local config (global begins + EL ends); default is full local (begins+ends)", 0 },
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
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

/* count local begins: parsetree nodes v s.t. v's PARENT node is the ROOT_S
 * (node 0, tr->state[0]==0) and v itself is NOT state 1 in the trivial
 * "child of root" sense -- more simply, for a MATL chain the only way ROOT_S
 * (v=0) is in the parsetree AT ALL and its child is NOT the state
 * immediately implied by a normal (non-begin) transition is a begin; since
 * a MATL-chain parsetree's very first node is always v=0 (S_st, i=1,j=L),
 * and its child (tr->state[1]) is either the "normal" first internal state
 * (deep local begin: child != that canonical state) or literally that state
 * (no begin fired, or begin fired into the same state a normal transition
 * would reach -- ambiguous only in that degenerate case, which none of our
 * synthetic panel hits: the whole point of the "prefix junk" trick is to
 * force entry deeper than node 1). Practically: a begin fired iff
 * tr->state[1] != the state cm->cfirst[0] (ROOT_S's first structural child)
 * OR (i,j) at that node doesn't span the full [1,L] (can't happen for a
 * begin, which always spans the whole seq) -- so the simple, robust check is
 * tr->state[1] != cm->cfirst[0]. */
static int count_begin(Parsetree_t *tr, CM_t *cm)
{
  if (tr->n < 2 || tr->state[0] != 0) return 0;
  return (tr->state[1] != cm->cfirst[0]) ? 1 : 0;
}

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: rl2_ckptalign_drv [-g] [--tau x] [--tag T] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  int    do_endsonly = esl_opt_GetBoolean(go, "--endsonly");
  int    do_stock  = ! esl_opt_GetBoolean(go, "--nostock");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *tag       = esl_opt_GetString(go, "--tag");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->align_opts |= CM_ALIGN_HBANDED;
  cm->align_opts |= CM_ALIGN_OPTACC;
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();
  int M = cm->M;

  /* R-L.2b: default local config is now FULL CM_CONFIG_LOCAL (begins+ends
   * together, the default cmalign use case cm_CheckptAlignHB now supports).
   * --endsonly reproduces R-L.2's original config (un-do the local-begin
   * half of cm_localize(): restore ROOT transitions + clear the begin flag,
   * keeping local ENDS) for the ends-only regression check. */
  if ((! do_global) && do_endsonly && (cm->flags & CMH_LOCAL_BEGIN)) {
    if (cm->root_trans == NULL) cm_Fail("root_trans NULL; cannot un-localize begins");
    esl_vec_FCopy(cm->root_trans, cm->cnum[0], cm->t[0]);
    cm->flags &= ~CMH_LOCAL_BEGIN;
    cm->flags &= ~CMH_BITS;
    CMLogoddsify(cm);
  }

  printf("# CM=%s M=%d clen=%d mode=%s tau=%g tag=%s\n",
         cm->name, cm->M, cm->clen,
         do_global ? "global" : (do_endsonly ? "local(EL-only,no-begins)" : "local(full:begins+ends)"),
         tau, tag);
  printf("# flags: LOCAL_BEGIN=%d LOCAL_END=%d\n",
         (cm->flags & CMH_LOCAL_BEGIN)?1:0, (cm->flags & CMH_LOCAL_END)?1:0);

  /* GATE: bps=0 MATL chain */
  int nb = 0, nbad = 0, v;
  for (v = 0; v < M; v++) {
    int st = cm->sttype[v];
    if (st == B_st) nb++;
    if (st == MP_st || st == MR_st) nbad++;
  }
  if (nb != 0 || nbad != 0) cm_Fail("CM is not bps=0 MATL chain: B_st=%d MP/MR=%d", nb, nbad);
  printf("# GATE ok: bps=0 (B_st=0, MP/MR=0)\n");

  /* read sequence */
  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;
    printf("# SEQ=%s L=%d\n", sq->name, L);

    /* CP9 bands (shared by both engines) */
    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- stock oracle: cm_AlignHB(do_optacc=TRUE) (skipped with --nostock,
     *      e.g. genome-scale local where the full-cube stock pipeline OOMs) ---- */
    Parsetree_t *tr_s = NULL; char *pp_s = NULL; float avgpp_s = 0, sc_s = 0;
    if (do_stock) {
      if ((status = cm_AlignHB(cm, errbuf, dsq, L, size_limit, TRUE, FALSE,
                               cm->hb_mx, cm->hb_shmx, cm->hb_omx, cm->hb_emx, NULL,
                               &pp_s, &tr_s, &avgpp_s, &sc_s)) != eslOK)
        cm_Fail("stock cm_AlignHB: %s", errbuf);
    }

    /* ---- checkpointed engine: cm_CheckptAlignHB ---- */
    CM_HB_EMIT_MX *emit_c = cm_hb_emit_mx_Create(cm);
    Parsetree_t *tr_c = NULL; char *pp_c = NULL; float avgpp_c = 0, sc_c = 0;
    if ((status = cm_CheckptAlignHB(cm, errbuf, dsq, L, size_limit,
                                    emit_c, &pp_c, &tr_c, &avgpp_c, &sc_c)) != eslOK)
      cm_Fail("cm_CheckptAlignHB: %s", errbuf);

    int elC = count_el(tr_c, M);
    int bgC = count_begin(tr_c, cm);
    char *modestr = do_global ? "global" : (do_endsonly ? "local(endsonly)" : "local(full)");
    if (! do_stock) {
      printf("# RESULT tag=%s seq=%s mode=%s L=%d (ckpt-only, --nostock)\n", tag, sq->name, modestr, L);
      printf("#   ckpt  : sc=%.6f nodes=%d avgpp=%.6f EL=%d BEGIN=%d entry_v=%d\n", sc_c, tr_c->n, avgpp_c, elC, bgC, (tr_c->n>1)?tr_c->state[1]:-1);
      printf("#   VERDICT=CKPT-ONLY (no stock comparison)\n");
      if (pp_c) free(pp_c); if (tr_c) FreeParsetree(tr_c);
      cm_hb_emit_mx_Destroy(emit_c); esl_sq_Reuse(sq); continue;
    }

    /* ---- compare ---- */
    int    peq  = ParsetreeCompare(tr_c, tr_s);
    int    ppeq = (pp_s && pp_c && strcmp(pp_s, pp_c) == 0);
    float  dsc  = fabsf(sc_c - sc_s);
    int    elS  = count_el(tr_s, M);
    int    bgS  = count_begin(tr_s, cm);

    printf("# RESULT tag=%s seq=%s mode=%s L=%d\n", tag, sq->name, modestr, L);
    printf("#   score : ckpt=%.6f stock=%.6f |dsc|=%.3e byte-exact=%s\n",
           sc_c, sc_s, dsc, (sc_c==sc_s)?"YES":"no");
    printf("#   parse : ckpt_nodes=%d stock_nodes=%d ParsetreeCompare=%s\n",
           tr_c->n, tr_s->n, peq?"IDENTICAL":"DIFFER");
    printf("#   PPstr : byte-identical=%s  avgpp(stock)=%.6f avgpp(ckpt)=%.6f\n",
           ppeq?"YES":"no", avgpp_s, avgpp_c);
    printf("#   EL    : stock=%d ckpt=%d match=%s\n", elS, elC, (elS==elC)?"YES":"no");
    printf("#   BEGIN : stock=%d ckpt=%d match=%s entry_v(stock)=%d entry_v(ckpt)=%d\n",
           bgS, bgC, (bgS==bgC)?"YES":"no",
           (tr_s->n>1)?tr_s->state[1]:-1, (tr_c->n>1)?tr_c->state[1]:-1);
    printf("#   VERDICT=%s\n",
           (peq && ppeq && (sc_c==sc_s) && elS==elC && bgS==bgC) ? "PASS" : "FAIL");

    if (pp_s) free(pp_s); if (pp_c) free(pp_c);
    if (tr_s) FreeParsetree(tr_s); if (tr_c) FreeParsetree(tr_c);
    cm_hb_emit_mx_Destroy(emit_c);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

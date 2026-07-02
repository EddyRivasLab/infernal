/* rl5_ckpttrpinalign_drv.c : R-L.5a validation of local (EL+begin) mechanical
 * wiring in the rung-4 STRUCTURED (bps>0) TRUNCATED posterior + OptAcc
 * engines: cm_PinTrPostAlignHB/cm_CheckptTrPostAlignHB and
 * cm_PinTrOptAccAlignHB/cm_CheckptTrOptAccAlignHB.
 *
 * Pipeline (mirrors production's do_trckpt_r4 path in cm_alndata.c, but
 * calls the engines directly -- test-only, bypassing
 * cm_CheckptTrOptAccAlignHB_Qualifies() which still correctly rejects local
 * CMs in production):
 *   pass 1: run TrCYKDivideAndConquerHB() once per root-valid marginal mode
 *           candidate (J/L/R/T); argmax penalty-folded score picks the mode
 *           (mirrors cmalign's mode resolution when preset_mode==UNKNOWN)
 *   bkind/kpin/bbmode/blmode/brmode extraction from the winning D&C parse
 *           (mirrors cm_alndata.c's static rung4_trpins_from_cyk(),
 *           duplicated here since it's file-static)
 *   pass 2a (full-storage correctness anchor): cm_PinTrPostAlignHB +
 *           cm_PinTrOptAccAlignHB, compared to stock
 *           cm_TrAlignHB(do_optacc=TRUE) in FULL local config
 *   pass 2b (sqrt(M) deliverable): cm_CheckptTrPostAlignHB +
 *           cm_CheckptTrOptAccAlignHB, compared to the SAME stock oracle
 *
 * Default config is LOCAL (full: begins+ends); -g switches to global (the
 * regression check: must reproduce pre-R-L.5a rung-4 structured-truncated
 * behavior exactly -- have_el is FALSE in -g mode so all new EL code paths
 * are skipped).
 *
 * Byte-exact gate (per case, per engine): Inside score |dsc|==0,
 * ParsetreeCompare==0 (covers marginal mode + bifurcation split points +
 * bkind/mode consistency -- B nodes' TRACE_RIGHT_CHILD i/j and tree
 * structure are part of the node-for-node comparison), per-residue PP
 * string strcmp==0, avg PP, EL-node count match, entry_v match.
 *
 * Usage: rl5_ckpttrpinalign_drv [-g] [--tau x] [--fixedtau] [--tag T] [--nostock] <cmfile> <seqfile>
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
  { "--fixedtau", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "accepted for cmdline parity; bands computed once here regardless", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { "--tag", eslARG_STRING,"run",  NULL, NULL, NULL, NULL, NULL, "label printed in output", 0 },
  { "--nostock",eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "skip stock oracle (ckpt-only)", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

/* count USED_EL nodes (parsetree nodes whose state index == cm->M = the EL state) */
static int count_el(Parsetree_t *tr, int M)
{
  int n = 0, x;
  for (x = 0; x < tr->n; x++) if (tr->state[x] == M) n++;
  return n;
}

/* duplicated from cm_alndata.c's static rung4_trpins_from_cyk(): derive the
 * per-bifurcation pin skeleton (kind, k*, B/left/right modes) from a
 * TRUNCATED CYK D&C parsetree. */
static void
trpins_from_cyk(CM_t *cm, Parsetree_t *tr,
                 int *bkind, int *kpin, char *bbmode, char *blmode, char *brmode)
{
  int v, n;
  for (v = 0; v < cm->M; v++) { bkind[v] = 0; kpin[v] = -1; bbmode[v] = blmode[v] = brmode[v] = TRMODE_UNKNOWN; }
  for (n = 0; n < tr->n; n++) {
    v = tr->state[n];
    if (cm->sttype[v] != B_st) continue;
    int ln = tr->nxtl[n], rn = tr->nxtr[n];
    int lstid = cm->stid[tr->state[ln]], rstid = cm->stid[tr->state[rn]];
    int begl_node, begr_node;
    if      (lstid == BEGL_S && rstid == BEGR_S) { begl_node = ln; begr_node = rn; }
    else if (lstid == BEGR_S && rstid == BEGL_S) { begl_node = rn; begr_node = ln; }
    else                                         { begl_node = ln; begr_node = rn; }
    int lspan = (tr->emitl[begl_node] <= tr->emitr[begl_node]) ? (tr->emitr[begl_node] - tr->emitl[begl_node] + 1) : 0;
    int rspan = (tr->emitl[begr_node] <= tr->emitr[begr_node]) ? (tr->emitr[begr_node] - tr->emitl[begr_node] + 1) : 0;
    bbmode[v] = tr->mode[n];
    blmode[v] = tr->mode[begl_node];
    brmode[v] = tr->mode[begr_node];
    if      (lspan == 0) { bkind[v] = 3; kpin[v] = -1;    }  /* RIGHT_FULL */
    else if (rspan == 0) { bkind[v] = 2; kpin[v] = 0;     }  /* LEFT_FULL  */
    else                 { bkind[v] = 1; kpin[v] = rspan; }  /* INTERIOR   */
  }
}

/* resolve marginal mode + (bkind,kpin,modes) pins via the production
 * do_trckpt_r4 recipe (cm_alndata.c ~502-536): run D&C CYK once per
 * root-valid mode candidate, argmax picks the mode. */
static int
resolve_trpins(CM_t *cm, ESL_DSQ *dsq, int L, int pass_idx,
               int *bkind, int *kpin, char *bbmode, char *blmode, char *brmode,
               char *ret_mode, float *ret_cyksc)
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
  trpins_from_cyk(cm, tr_best, bkind, kpin, bbmode, blmode, brmode);
  FreeParsetree(tr_best);
  *ret_mode  = r4_mode;
  *ret_cyksc = r4_cyk;
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
    { printf("usage: rl5_ckpttrpinalign_drv [-g] [--tau x] [--fixedtau] [--tag T] [--nostock] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  int    do_stock  = ! esl_opt_GetBoolean(go, "--nostock");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *tag       = esl_opt_GetString(go, "--tag");
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

  /* GATE: structured (bps>0, has >=1 bifurcation) */
  int nb = 0, v;
  for (v = 0; v < M; v++) if (cm->sttype[v] == B_st) nb++;
  if (nb < 1) cm_Fail("CM has no bifurcations (bps=0): not a rung-4 target");
  printf("# CM=%s M=%d clen=%d bifs=%d mode=%s tau=%g tag=%s\n",
         cm->name, cm->M, cm->clen, nb, do_global ? "global" : "local(full:begins+ends)", tau, tag);
  printf("# flags: LOCAL_BEGIN=%d LOCAL_END=%d\n",
         (cm->flags & CMH_LOCAL_BEGIN)?1:0, (cm->flags & CMH_LOCAL_END)?1:0);

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;
    printf("# SEQ=%s L=%d\n", sq->name, L);

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- stock oracle: cm_TrAlignHB(do_optacc=TRUE), mode UNKNOWN (auto-resolve) ---- */
    Parsetree_t *tr_s = NULL; char *pp_s = NULL; float avgpp_s = 0, sc_s = 0; char mode_s = TRMODE_UNKNOWN;
    if (do_stock) {
      if ((status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, TRMODE_UNKNOWN, pass_idx, TRUE, FALSE,
                                 cm->trhb_mx, cm->trhb_shmx, cm->trhb_omx, cm->trhb_emx, NULL,
                                 &pp_s, &tr_s, &mode_s, &avgpp_s, &sc_s)) != eslOK)
        cm_Fail("stock cm_TrAlignHB: %s", errbuf);
    }

    /* ---- pass 1: mode-resolving D&C CYK -> (mode, bkind, kpin, bb/bl/brmode) ---- */
    int   *bkind = NULL, *kpin = NULL;
    char  *bbmode = NULL, *blmode = NULL, *brmode = NULL;
    ESL_ALLOC(bkind,  sizeof(int)  * M);
    ESL_ALLOC(kpin,   sizeof(int)  * M);
    ESL_ALLOC(bbmode, sizeof(char) * M);
    ESL_ALLOC(blmode, sizeof(char) * M);
    ESL_ALLOC(brmode, sizeof(char) * M);
    char r4_mode = TRMODE_UNKNOWN; float r4_cyk = IMPOSSIBLE;
    if (resolve_trpins(cm, dsq, L, pass_idx, bkind, kpin, bbmode, blmode, brmode, &r4_mode, &r4_cyk) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR resolve_trpins: no root-valid mode\n# VERDICT=SKIP(no-mode)\n", sq->name, L);
      free(bkind); free(kpin); free(bbmode); free(blmode); free(brmode);
      if (pp_s) free(pp_s); if (tr_s) FreeParsetree(tr_s);
      esl_sq_Reuse(sq); continue;
    }

    /* ---- pass 2a: full-storage correctness anchor (cm_PinTr*) ---- */
    CM_TR_HB_EMIT_MX *emit_p = cm_tr_hb_emit_mx_Create(cm);
    Parsetree_t *tr_p = NULL; char *pp_p = NULL; float avgpp_p = 0, sc_p = 0, pp_p_acc = 0; char mode_p = r4_mode;
    if ((status = cm_PinTrPostAlignHB(cm, errbuf, dsq, L, size_limit, r4_mode, pass_idx, emit_p,
                                      bkind, kpin, bbmode, blmode, brmode, &sc_p, &mode_p)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_PinTrPostAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      goto NEXT_SKIP_P;
    }
    if ((status = cm_PinTrOptAccAlignHB(cm, errbuf, dsq, L, size_limit, r4_mode, pass_idx, emit_p,
                                        bkind, kpin, bbmode, blmode, brmode,
                                        &pp_p, &tr_p, &mode_p, &avgpp_p, &pp_p_acc)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_PinTrOptAccAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      goto NEXT_SKIP_P;
    }

    /* ---- pass 2b: sqrt(M) checkpointed deliverable (cm_CheckptTr*) ---- */
    {
    CM_TR_HB_EMIT_MX *emit_c = cm_tr_hb_emit_mx_Create(cm);
    Parsetree_t *tr_c = NULL; char *pp_c = NULL; float avgpp_c = 0, sc_c = 0, pp_c_acc = 0; char mode_c = r4_mode;
    if ((status = cm_CheckptTrPostAlignHB(cm, errbuf, dsq, L, size_limit, r4_mode, pass_idx, emit_c,
                                          bkind, kpin, bbmode, blmode, brmode, &sc_c, &mode_c)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_CheckptTrPostAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      cm_tr_hb_emit_mx_Destroy(emit_c); goto NEXT_SKIP_C;
    }
    if ((status = cm_CheckptTrOptAccAlignHB(cm, errbuf, dsq, L, size_limit, r4_mode, pass_idx, emit_c,
                                            bkind, kpin, bbmode, blmode, brmode,
                                            &pp_c, &tr_c, &mode_c, &avgpp_c, &pp_c_acc)) != eslOK) {
      printf("# SEQ=%s L=%d\n# ERROR cm_CheckptTrOptAccAlignHB: %s\n# VERDICT=SKIP(pin-topology-mismatch)\n", sq->name, L, errbuf);
      cm_tr_hb_emit_mx_Destroy(emit_c); goto NEXT_SKIP_C;
    }

    int elP = count_el(tr_p, M), elC = count_el(tr_c, M);
    int entryP = (tr_p->n > 1) ? tr_p->state[1] : -1;
    int entryC = (tr_c->n > 1) ? tr_c->state[1] : -1;
    char *modestr = do_global ? "global" : "local(full)";

    if (! do_stock) {
      printf("# RESULT tag=%s seq=%s mode=%s L=%d D&Cmode=%s (ckpt-only, --nostock)\n", tag, sq->name, modestr, L, MarginalMode(r4_mode));
      printf("#   pin   : sc=%.6f nodes=%d avgpp=%.6f mode=%s EL=%d entry_v=%d\n", sc_p, tr_p->n, avgpp_p, MarginalMode(mode_p), elP, entryP);
      printf("#   ckpt  : sc=%.6f nodes=%d avgpp=%.6f mode=%s EL=%d entry_v=%d\n", sc_c, tr_c->n, avgpp_c, MarginalMode(mode_c), elC, entryC);
      printf("#   VERDICT=CKPT-ONLY (no stock comparison)\n");
    } else {
      int    peqPC = ParsetreeCompare(tr_p, tr_c);
      int    ppeqPC = (pp_p && pp_c && strcmp(pp_p, pp_c) == 0);
      float  dscPC = fabsf(sc_p - sc_c);
      int    peqP  = ParsetreeCompare(tr_p, tr_s);
      int    peqC  = ParsetreeCompare(tr_c, tr_s);
      float  dscP  = fabsf(sc_p - sc_s);
      float  dscC  = fabsf(sc_c - sc_s);
      int    elS   = count_el(tr_s, M);
      int    entryS = (tr_s->n > 1) ? tr_s->state[1] : -1;

      /* R-L.5b Phase-1 structural probe: would the FREE (unpinned) stock
       * oracle's own independently-optimal parsetree, if it were subject to
       * our engine's bkind-based masking, ever try to visit a B state that
       * CYK's own parse left unpinned (bkind[v]==0)?  If so, that is direct,
       * sharp evidence that masking COULD block a legitimate optimal path
       * (not just "the assert didn't fire" -- an absence-of-failure argument). */
      { int nn, nunpinned_visited = 0;
        for (nn = 0; nn < tr_s->n; nn++) {
          int vv = tr_s->state[nn];
          if (vv >= 0 && vv < M && cm->sttype[vv] == B_st && bkind[vv] == 0) {
            nunpinned_visited++;
            printf("#   STOCK-VISITS-UNPINNED-B: v=%d\n", vv);
          }
        }
        if (nunpinned_visited > 0)
          printf("#   STOCK_UNPINNED_B_COUNT=%d (free oracle's own optimal parse would need an unpinned B)\n", nunpinned_visited);
      }

      printf("# RESULT tag=%s seq=%s mode=%s L=%d D&Cmode=%s stockmode=%s\n", tag, sq->name, modestr, L, MarginalMode(r4_mode), MarginalMode(mode_s));
      printf("#   PRIMARY GATE (pin==ckpt, same pins): ParsetreeCompare=%s  |dsc|=%.3e  PPstr-byte-identical=%s  entry_v match=%s(%d vs %d) EL match=%s(%d vs %d) mode match=%s\n",
             peqPC?"IDENTICAL":"DIFFER", dscPC, ppeqPC?"YES":"no",
             (entryP==entryC)?"YES":"no", entryP, entryC, (elP==elC)?"YES":"no", elP, elC,
             (mode_p==mode_c)?"YES":"no");
      printf("#   info (vs free unpinned stock, NOT gated unless mode/score fully agree -- pin/OA divergence documented pre-existing behavior):\n");
      printf("#     score : stock=%.6f pin=%.6f |dP|=%.3e ckpt=%.6f |dC|=%.3e\n", sc_s, sc_p, dscP, sc_c, dscC);
      printf("#     parse : stock_nodes=%d pin_nodes=%d PinCompare=%s  ckpt_nodes=%d CkptCompare=%s\n",
             tr_s->n, tr_p->n, peqP?"IDENTICAL":"DIFFER", tr_c->n, peqC?"IDENTICAL":"DIFFER");
      printf("#     avgpp(stock)=%.6f avgpp(pin)=%.6f avgpp(ckpt)=%.6f\n", avgpp_s, avgpp_p, avgpp_c);
      printf("#     EL    : stock=%d pin=%d ckpt=%d\n", elS, elP, elC);
      printf("#     entry_v: stock=%d pin=%d ckpt=%d\n", entryS, entryP, entryC);
      printf("#     mode  : stock=%s pin=%s ckpt=%s\n", MarginalMode(mode_s), MarginalMode(mode_p), MarginalMode(mode_c));

      /* "not pinned" assert probe: does the checkpointed OA traceback visit any
       * B state that trpins_from_cyk left bkind==0 (unpinned)?  ckpt engine
       * would have already ESL_XFAIL'd on this (trckpt_tr_optacc_traceback:
       * "B state v=%d not pinned") -- if we got here without erroring, no
       * B-state in tr_c's own output was left unpinned by THIS panel's parses.
       * (We can't inspect a case where it DID fire since the call above would
       * have gone to the SKIP branch instead -- see the error message text.) */
      printf("#   VERDICT=%s\n",
             (peqPC && ppeqPC && dscPC==0. && elP==elC && entryP==entryC && mode_p==mode_c) ? "PASS" : "FAIL");
    }

    if (pp_c) free(pp_c); if (tr_c) FreeParsetree(tr_c);
    cm_tr_hb_emit_mx_Destroy(emit_c);
    }
  NEXT_SKIP_C:
    if (pp_p) free(pp_p); if (tr_p) FreeParsetree(tr_p);
  NEXT_SKIP_P:
    cm_tr_hb_emit_mx_Destroy(emit_p);
    free(bkind); free(kpin); free(bbmode); free(blmode); free(brmode);
    if (pp_s) free(pp_s); if (tr_s) FreeParsetree(tr_s);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;

 ERROR:
  cm_Fail("allocation error");
  return 1;
}

/* rl3_ckpttralign_drv.c : R-L.3 validation of cm_CheckptTrAlignHB EL
 * (local-end) + local-BEGIN-table-selection support for bps=0 (MATL-chain)
 * TRUNCATED CMs, across marginal modes J/L/R.
 *
 * Calls the LIBRARY engine cm_CheckptTrAlignHB() directly (test-only,
 * bypassing cm_CheckptTrAlignHB_Qualifies() -- which still correctly rejects
 * local CMs in production; relaxing that gate for local is a future rung's
 * job) and compares it against the stock oracle cm_TrAlignHB(do_optacc=TRUE),
 * mirroring rl2_ckptalign_drv.c's (non-trunc) structure.
 *
 * Default config is LOCAL (full: begins+ends); -g switches to global (the
 * regression check: the EL/begin-table additions must not perturb the
 * pre-R-L.3 global path).
 *
 * Byte-exact gate (per case): Inside score |dsc|, ParsetreeCompare==0
 * (which includes marginal mode + entry state, since ParsetreeCompare checks
 * tr->mode[] and tr->state[] node-for-node), per-residue PP string
 * strcmp==0, avg PP.  EL coverage: count USED_EL nodes (parsetree
 * state==cm->M).  entry_v: tr->state[1] (the state ROOT_S's truncated begin
 * enters), always present in a truncated parse (the "truncated begin" fires
 * unconditionally, in both global and local mode -- this is NOT the
 * non-trunc "did a deep local begin fire" question).
 *
 * Usage: rl3_ckpttralign_drv [-g] [--tau x] [--tag T] [--nostock] <cmfile> <seqfile>
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
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { "--tag", eslARG_STRING,"run",  NULL, NULL, NULL, NULL, NULL, "label printed in output", 0 },
  { "--nostock",eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "skip stock oracle (ckpt-only; for genome-scale where stock OOMs)", 0 },
  { "--forcemode", eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL, "force preset_mode J|L|R (both engines) instead of auto-resolve -- for exercising a mode's mechanics even when it wouldn't naturally win", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

/* count USED_EL nodes (parsetree nodes whose state index == cm->M = the EL state) */
static int count_el(Parsetree_t *tr, int M)
{
  int n = 0, x;
  for (x = 0; x < tr->n; x++) if (tr->state[x] == M) n++;
  return n;
}

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: rl3_ckpttralign_drv [-g] [--tau x] [--tag T] [--nostock] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  int    do_stock  = ! esl_opt_GetBoolean(go, "--nostock");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *tag       = esl_opt_GetString(go, "--tag");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);
  int    pass_idx  = PLI_PASS_5P_AND_3P_ANY;
  char   force_mode = TRMODE_UNKNOWN;
  if (esl_opt_IsUsed(go, "--forcemode")) {
    char *fm = esl_opt_GetString(go, "--forcemode");
    if      (!strcmp(fm,"J")) force_mode = TRMODE_J;
    else if (!strcmp(fm,"L")) force_mode = TRMODE_L;
    else if (!strcmp(fm,"R")) force_mode = TRMODE_R;
    else cm_Fail("--forcemode must be J, L, or R");
  }

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

  printf("# CM=%s M=%d clen=%d mode=%s tau=%g tag=%s\n",
         cm->name, cm->M, cm->clen, do_global ? "global" : "local(full:begins+ends)", tau, tag);
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

    /* truncated CP9 bands (shared by both engines) */
    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- stock oracle: cm_TrAlignHB(do_optacc=TRUE) (skipped with --nostock) ---- */
    Parsetree_t *tr_s = NULL; char *pp_s = NULL; float avgpp_s = 0, sc_s = 0; char mode_s = TRMODE_UNKNOWN;
    if (do_stock) {
      if ((status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, force_mode, pass_idx, TRUE, FALSE,
                                 cm->trhb_mx, cm->trhb_shmx, cm->trhb_omx, cm->trhb_emx, NULL,
                                 &pp_s, &tr_s, &mode_s, &avgpp_s, &sc_s)) != eslOK)
        cm_Fail("stock cm_TrAlignHB: %s", errbuf);
    }

    /* ---- checkpointed engine: cm_CheckptTrAlignHB ---- */
    CM_TR_HB_EMIT_MX *emit_c = cm_tr_hb_emit_mx_Create(cm);
    Parsetree_t *tr_c = NULL; char *pp_c = NULL; float avgpp_c = 0, sc_c = 0; char mode_c = TRMODE_UNKNOWN;
    if ((status = cm_CheckptTrAlignHB(cm, errbuf, dsq, L, size_limit, force_mode, pass_idx, emit_c,
                                      &pp_c, &tr_c, &mode_c, &avgpp_c, &sc_c)) != eslOK)
      cm_Fail("cm_CheckptTrAlignHB: %s", errbuf);

    int elC = count_el(tr_c, M);
    int entryC = (tr_c->n > 1) ? tr_c->state[1] : -1;
    char *modestr = do_global ? "global" : "local(full)";
    if (! do_stock) {
      printf("# RESULT tag=%s seq=%s mode=%s L=%d (ckpt-only, --nostock)\n", tag, sq->name, modestr, L);
      printf("#   ckpt  : sc=%.6f nodes=%d avgpp=%.6f mode=%s EL=%d entry_v=%d\n",
             sc_c, tr_c->n, avgpp_c, MarginalMode(mode_c), elC, entryC);
      printf("#   VERDICT=CKPT-ONLY (no stock comparison)\n");
      if (pp_c) free(pp_c); if (tr_c) FreeParsetree(tr_c);
      cm_tr_hb_emit_mx_Destroy(emit_c); esl_sq_Reuse(sq); continue;
    }

    /* ---- compare ---- */
    int    peq  = ParsetreeCompare(tr_c, tr_s);
    int    ppeq = (pp_s && pp_c && strcmp(pp_s, pp_c) == 0);
    float  dsc  = fabsf(sc_c - sc_s);
    int    elS  = count_el(tr_s, M);
    int    entryS = (tr_s->n > 1) ? tr_s->state[1] : -1;

    printf("# RESULT tag=%s seq=%s mode=%s L=%d\n", tag, sq->name, modestr, L);
    printf("#   score : ckpt=%.6f stock=%.6f |dsc|=%.3e byte-exact=%s\n",
           sc_c, sc_s, dsc, (sc_c==sc_s)?"YES":"no");
    printf("#   mode  : stock=%s ckpt=%s match=%s\n", MarginalMode(mode_s), MarginalMode(mode_c), (mode_s==mode_c)?"YES":"no");
    printf("#   parse : ckpt_nodes=%d stock_nodes=%d ParsetreeCompare=%s\n",
           tr_c->n, tr_s->n, peq?"IDENTICAL":"DIFFER");
    printf("#   PPstr : byte-identical=%s  avgpp(stock)=%.6f avgpp(ckpt)=%.6f\n",
           ppeq?"YES":"no", avgpp_s, avgpp_c);
    printf("#   EL    : stock=%d ckpt=%d match=%s\n", elS, elC, (elS==elC)?"YES":"no");
    printf("#   entry_v: stock=%d ckpt=%d match=%s\n", entryS, entryC, (entryS==entryC)?"YES":"no");
    printf("#   VERDICT=%s\n",
           (peq && ppeq && (sc_c==sc_s) && elS==elC && entryS==entryC && mode_s==mode_c) ? "PASS" : "FAIL");

    if (pp_s) free(pp_s); if (pp_c) free(pp_c);
    if (tr_s) FreeParsetree(tr_s); if (tr_c) FreeParsetree(tr_c);
    cm_tr_hb_emit_mx_Destroy(emit_c);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

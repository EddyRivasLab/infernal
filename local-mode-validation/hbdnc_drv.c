/* hbdnc_drv.c : isolated correctness + memory driver for HMM-banded (CP9)
 * divide-and-conquer CYK (brief 26_0610-007, Stage 1a.1).
 *
 * For each sequence in <seqfile>, derives CP9 bands and aligns it:
 *   route A: cm_AlignHB() CYK (full banded matrix CM_HB_MX)   [oracle: parsetree A]
 *   route B: CYKDivideAndConquerHB() (banded D&C, under test) [parsetree B]
 *   route X: CYKInside() with NULL bands (exact non-banded full CYK) [clip detection]
 *
 * Pass criterion (Stage 1a.1): score(B)==score(A) ALWAYS; and on NON-CLIP
 * sequences (banded score == exact score) parsetree B must be byte-identical
 * to parsetree A. For clipping seqs B may differ from A (class-2 unbanded);
 * those are recorded, not failed.
 *
 * Also reports, per sequence, the D&C-HB peak live-deck working set
 * (CYKDeckTrackMaxMb) vs the full CM_HB_MX size (cm->hb_mx->size_Mb) and the
 * x-factor (capacity-win evidence).
 *
 * Read-only wrt infernal source; links libinfernal.
 * Usage: hbdnc_drv [-g] [--tau <x>] <cmfile> <seqfile>
 *   default: local config (HMM local + EL), matching cmsearch; -g: global.
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
  { "-h",    eslARG_NONE,   FALSE,   NULL, NULL, NULL, NULL, NULL, "help", 0 },
  { "-g",        eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "global config (default local)", 0 },
  { "--tau",     eslARG_REAL, "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--noexact", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "skip exact non-banded CYK (route X) [for large CMs]", 0 },
  { "--nofull",  eslARG_NONE, FALSE, NULL, NULL, NULL,"--noexact",NULL, "skip full banded CYK (route A) too; run only D&C (route B)", 0 },
  { "--mxsize",  eslARG_REAL, "4096",NULL, NULL, NULL, NULL, NULL, "size limit (Mb) for banded matrices", 0 },
  { "--sizeonly",eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "only derive bands + allocate CM_HB_MX, print size_Mb, skip alignment", 0 },
  { "--dumpB",   eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL, "dump D&C parsetree B (state/emitl/emitr per node) to <f> for byte-exact regression", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

static int trees_equal(Parsetree_t *a, Parsetree_t *b)
{
  int i;
  if (a == NULL || b == NULL) return 0;
  if (a->n != b->n) return 0;
  for (i = 0; i < a->n; i++) {
    if (a->state[i] != b->state[i]) return 0;
    if (a->emitl[i] != b->emitl[i]) return 0;
    if (a->emitr[i] != b->emitr[i]) return 0;
  }
  return 1;
}

/* count_el: number of EL local-end nodes in a parsetree (state == cm->M,
 * the EL state recorded by InsertTraceNode at cm_dpsmall.c:7705/8497).
 * This is the unambiguous "EL local-end taken" signal (USED_EL transition). */
static int count_el(Parsetree_t *tr, int M)
{
  int i, n = 0;
  if (tr == NULL) return 0;
  for (i = 0; i < tr->n; i++) if (tr->state[i] == M) n++;
  return n;
}

/* local-begin entry state: in a parsetree, node 0 is ROOT_S (state 0); node 1
 * is its child. A global parse enters via ROOT_IL/IR/the first split state
 * (state <= 2 or the canonical next state); a local begin (0->r) enters at an
 * internal subtree state r. We report that entry state so the caller can flag
 * whether a local-begin was actually exercised. Returns -1 if no child. */
static int entry_state(Parsetree_t *tr)
{
  if (tr == NULL || tr->n < 2) return -1;
  return tr->state[1];
}

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go) != eslOK) { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h")) { printf("usage: hbdnc_drv [-g] [--tau <x>] <cmfile> <seqfile>\n"); exit(0); }
  int   do_global = esl_opt_GetBoolean(go, "-g");
  float tau       = esl_opt_GetReal(go, "--tau");
  int   do_exact  = ! esl_opt_GetBoolean(go, "--noexact");
  int   do_full   = ! esl_opt_GetBoolean(go, "--nofull");
  if (esl_opt_ArgNumber(go) != 2) { printf("usage: hbdnc_drv [-g] [--tau <x>] <cmfile> <seqfile>\n"); exit(1); }
  char *cmfile  = esl_opt_GetArg(go, 1);
  char *seqfile = esl_opt_GetArg(go, 2);

  /* read CM */
  CM_FILE *cmfp = NULL;
  ESL_ALPHABET *abc = NULL;
  CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);

  /* configure for HMM banded alignment; local unless -g */
  cm->align_opts |= CM_ALIGN_HBANDED;
  if (! do_global) {
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);

  init_ilogsum();
  FLogsumInit();

  float size_limit = esl_opt_GetReal(go, "--mxsize");

  printf("# CM: %s  M=%d clen=%d  mode=%s  tau=%g\n",
         cm->name, cm->M, cm->clen, do_global?"global":"local", tau);
  printf("# %-22s %10s %10s %10s  %-26s %10s %10s %8s\n",
         "seqname/L", "hbAlign(A)", "hbDnC(B)", "exact(X)", "verdict",
         "DnCpeakMb", "HBmxMb", "xfactor");

  FILE *dumpfp = NULL;
  char *dumpfile = esl_opt_GetString(go, "--dumpB");
  if (dumpfile != NULL) {
    if ((dumpfp = fopen(dumpfile, "w")) == NULL) cm_Fail("can't open --dumpB file %s", dumpfile);
  }

  ESL_SQFILE *sqfp = NULL;
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK)
    cm_Fail("open seqfile %s", seqfile);
  ESL_SQ *sq = esl_sq_CreateDigital(abc);

  int n_nonclip = 0, n_nonclip_treeok = 0, n_score_mismatch = 0, n_clip = 0, n_bug = 0;
  int n_el_mismatch = 0;          /* cases where elA != elB (EL structure diverged) */
  int max_el_A = 0, max_el_B = 0; /* panel-wide max EL nodes seen (coverage) */
  int n_localbegin = 0;          /* cases whose D&C parse entered via a local begin (entry state > 2) */

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    Parsetree_t *trA = NULL, *trB = NULL, *trX = NULL;
    float scA = -9999., scB, scX = -9999., ppA;
    float hbmx_Mb = 0.0;

    /* derive CP9 bands for this seq -> cm->cp9b */
    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                sq->dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    if (esl_opt_GetBoolean(go, "--sizeonly")) {
      if ((status = cm_hb_mx_GrowTo(cm, cm->hb_mx, errbuf, cm->cp9b, L, esl_opt_GetReal(go,"--mxsize"))) != eslOK)
        cm_Fail("cm_hb_mx_GrowTo: %s", errbuf);
      char nm[40]; snprintf(nm, sizeof(nm), "%.16s/%d", sq->name, L);
      printf("  %-22s  CM_HB_MX size_Mb = %.4f\n", nm, cm->hb_mx->size_Mb);
      esl_sq_Reuse(sq);
      continue;
    }

    /* route A: full banded CYK (oracle). Grows cm->hb_mx to the banded size. */
    if (do_full) {
      if ((status = cm_AlignHB(cm, errbuf, sq->dsq, L, size_limit, FALSE, FALSE,
                               cm->hb_mx, cm->hb_shmx, NULL, NULL, NULL, NULL,
                               &trA, &ppA, &scA)) != eslOK)
        cm_Fail("cm_AlignHB: %s", errbuf);
      hbmx_Mb = cm->hb_mx->size_Mb;
    }

    /* route B: banded D&C (under test), instrumented for peak live-deck bytes */
    CYKDeckTrackReset();
    scB = CYKDivideAndConquerHB(cm, sq->dsq, L, 0, 1, L, &trB, cm->cp9b);
    double dnc_peak_Mb = CYKDeckTrackMaxMb();
    double dnc_vjd_Mb  = CYKDeckTrackVjdAtPeakMb(); /* class-1 banded decks at peak */
    double dnc_vji_Mb  = CYKDeckTrackVjiAtPeakMb(); /* class-2 V-problem decks at peak */

    /* route X: exact non-banded full CYK (clip detection) */
    if (do_exact)
      scX = CYKInside(cm, sq->dsq, L, 0, 1, L, &trX, NULL, NULL);

    /* classify */
    char verdict[64];
    int is_clip = (do_exact && scB < scX - 0.01);  /* banded lost the optimum */
    if (do_exact && scB > scX + 0.01) {
      snprintf(verdict, sizeof(verdict), "BUG:B>X(%.3f)", scB-scX);
      n_bug++;
    } else if (do_full && fabs(scB - scA) > 0.01) {
      snprintf(verdict, sizeof(verdict), "BUG:B!=A(%.3f)", scB-scA);
      n_score_mismatch++;
    } else if (is_clip) {
      snprintf(verdict, sizeof(verdict), "clip(B-X=%.3f)%s", scB-scX,
               do_full ? (trees_equal(trB, trA) ? " tree==" : " tree!=") : "");
      n_clip++;
    } else if (! do_full) {
      snprintf(verdict, sizeof(verdict), "B-only(score%.3f)", scB);
    } else {
      /* non-clip: require byte-identical parsetree B == A */
      if (trees_equal(trB, trA)) { snprintf(verdict, sizeof(verdict), "OK tree=="); n_nonclip_treeok++; }
      else                       { snprintf(verdict, sizeof(verdict), "FAIL tree!="); }
      n_nonclip++;
    }

    /* EL local-end coverage: count USED_EL nodes (state==cm->M) in each parse.
     * For non-clip cases the parses are byte-identical so elA==elB; we still
     * track it for the per-case "EL taken?" column and panel-wide coverage. */
    int elA = do_full ? count_el(trA, cm->M) : -1;
    int elB = count_el(trB, cm->M);
    int esB = entry_state(trB);
    if (do_full && elA != elB && !is_clip) n_el_mismatch++;
    if (elA > max_el_A) max_el_A = elA;
    if (elB > max_el_B) max_el_B = elB;
    if (esB > 2) n_localbegin++;   /* entered at an internal state => local begin */

    double xfac = (dnc_peak_Mb > 0.0) ? (hbmx_Mb / dnc_peak_Mb) : 0.0;
    char nm[40]; snprintf(nm, sizeof(nm), "%.16s/%d", sq->name, L);
    printf("  %-22s %10.3f %10.3f %10.3f  %-26s %10.4f %10.4f %7.1fx  vjd=%.4f vji=%.4f  elA=%d elB=%d entryB=%d\n",
           nm, scA, scB, scX, verdict, dnc_peak_Mb, hbmx_Mb, xfac, dnc_vjd_Mb, dnc_vji_Mb, elA, elB, esB);

    if (dumpfp != NULL) {
      fprintf(dumpfp, ">%s/%d n=%d scB=%.3f\n", sq->name, L, trB ? trB->n : -1, scB);
      if (trB != NULL) {
        int ti;
        for (ti = 0; ti < trB->n; ti++)
          fprintf(dumpfp, "%d %d %d %d\n", ti, trB->state[ti], trB->emitl[ti], trB->emitr[ti]);
      }
    }

    if (trA) FreeParsetree(trA);
    if (trB) FreeParsetree(trB);
    if (trX) FreeParsetree(trX);
    esl_sq_Reuse(sq);
  }

  printf("# summary: non-clip=%d (tree==: %d), clip=%d, score-mismatch(B!=A)=%d, bug(B>X)=%d\n",
         n_nonclip, n_nonclip_treeok, n_clip, n_score_mismatch, n_bug);
  printf("# EL coverage: max EL-nodes/parse  A=%d  B=%d ; EL-structure mismatches(non-clip)=%d ; local-begin cases=%d\n",
         max_el_A, max_el_B, n_el_mismatch, n_localbegin);
  if (n_nonclip == n_nonclip_treeok && n_score_mismatch == 0 && n_bug == 0 && n_el_mismatch == 0)
    printf("# RESULT: PASS (byte-exact on all non-clip; B==A score on all; EL structure matches)\n");
  else
    printf("# RESULT: FAIL\n");

  if (dumpfp != NULL) fclose(dumpfp);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

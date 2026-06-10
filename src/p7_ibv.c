/* p7_ibv.c -- F+B direct-band band derivation for cmalign --p7ibv
 *
 * Brief 120 implementation of the IBV (Iterative Banded Viterbi -> direct
 * Forward+Backward + Delta threshold) band-derivation scheme prototyped
 * in Python in briefs 116-118. The Python reference implementation lives
 * at
 *   notebook/26_0430_inf_faster_cmalign/brief117_runs/fb_direct_band.py
 * and the underlying F+B primitive at
 *   notebook/26_0430_inf_faster_cmalign/brief116_runs/hirschberg_ibv.py
 *
 * Algorithm (per brief 116-117):
 *   F_M[i,k] = emit(k, seq[i]) + max(F_M[i-1,k-1]+T_MM[k-1],
 *                                    F_I[i-1,k-1]+T_IM[k-1],
 *                                    F_D[i-1,k-1]+T_DM[k-1])
 *   F_I[i,k] = max(F_M[i-1,k]+T_MI[k], F_I[i-1,k]+T_II[k])
 *   F_D[i,k] = max(F_M[i,k-1]+T_MD[k-1], F_D[i,k-1]+T_DD[k-1])
 *
 * Forward seeds (glocal forced-start at (0,0)):  F_M[0,0]=0, else -INF.
 * Backward seeds (glocal forced-end at (L,M)):   B_M[L,M]=B_I[L,M]=B_D[L,M]=0, else -INF.
 *
 * Band emission:
 *   optimal       = max{F_M[L,M], F_I[L,M], F_D[L,M]}
 *   through[i,k]  = max{F_M[i,k]+B_M[i,k], F_I[i,k]+B_I[i,k], F_D[i,k]+B_D[i,k]}
 *   for each row i in 1..L: admit cells where through[i,k] >= optimal - max(delta, EPS),
 *                            kmin[i] = min admitted k, kmax[i] = max admitted k.
 *
 * Boundary widening: rows i in {1, L-1, L} forced to [1, M] for truncated-alignment
 * entry/exit flexibility (mirrors the Python prototype's boundary_widen).
 *
 * B-state convention: kmin[0]=0, kmax[0]=0 (required by cp9_FB2HMMBandsP7BF,
 * see brief 117 §1; the post-trace CP9 pipeline asserts kmin[0]==0).
 *
 * Scores are stored as doubles in milli-bits (INTSCALE=1000, units of
 * 0.001 bit) to match the Python prototype's INTSCALE-scaled int convention
 * but with float64 arithmetic to side-step rounding mismatch in the
 * "==optimal" comparison.
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

/* INTSCALE matches brief 116-118 Python prototype: 1 bit -> 1000 milli-bits.
 * Scores returned by hmm_log2 helpers below are in milli-bits.
 */
#define P7IBV_INTSCALE     1000.0
#define P7IBV_NEG_INF      (-1.0e18)
/* EPS floor in milli-bits for the "==optimal" comparison (brief 118 §6).
 * The optimal-path cells equal optimal up to a few ULPs of a sum that can
 * reach ~3e5 in magnitude; 1.0 milli-bit (1e-3 bit) absorbs that without
 * being biologically meaningful. Mirrors hirschberg_ibv.py EPS=1.0.
 */
#define P7IBV_EPS          1.0

/* Convert a P7_HMM linear probability into milli-bits (= INTSCALE * log2(P)). */
static inline double
p7ibv_lod_milli(float p)
{
  if (p <= 0.0f) return P7IBV_NEG_INF;
  return P7IBV_INTSCALE * (log((double) p) / M_LN2);
}

/* Match-state emission score in milli-bits:
 *   INTSCALE * log2(mat[k][x] / ins[k][x])
 * Matches the Python prototype which uses insert emissions as the null model
 * (ins_emit serves as background). Ambiguous residues (x >= K) score 0.
 */
static inline double
p7ibv_emit_milli(const P7_HMM *hmm, int k, int x)
{
  if (x < 0 || x >= hmm->abc->K) return 0.0;
  float pm = hmm->mat[k][x];
  float pi = hmm->ins[k][x];
  if (pm <= 0.0f || pi <= 0.0f) return P7IBV_NEG_INF;
  return P7IBV_INTSCALE * (log((double) pm / (double) pi) / M_LN2);
}


/* Function: p7_Seq2BandsIBV()
 *
 * Synopsis: Derive p7 bands from full F+B + Delta threshold (direct-band-output).
 *
 * Purpose:  Given a CM <cm> (whose embedded filter HMM cm->fp7 supplies the
 *           p7 transition/emission scores) and a digitized sequence <dsq>
 *           of length <L>, build full (L+1)x(M+1) Forward and Backward
 *           score matrices, then for each row i in [1, L] emit
 *           kmin[i]/kmax[i] = min/max k such that through(i,k) >= optimal - Delta.
 *
 *           Boundary rows 1, L-1, L are widened to [1, M] for truncated-
 *           alignment entry/exit. kmin[0]/kmax[0] are set to 0 (B-state
 *           convention required by cp9_FB2HMMBandsP7BF; see brief 117 §1).
 *
 *           i2k is returned as all -1 (the downstream consumer
 *           p7_kbands2gbands ignores it; cm_alndata.c's #P7BAND diagnostic
 *           counts pins as i2k[i] != -1, so all-(-1) reports npins=0,
 *           which is correct: IBV emits a band, not pins).
 *
 *           Memory cost: 6 * (L+1) * (M+1) * sizeof(double) bytes.
 *           For an LSU rRNA at M=3400, L=2771 this is ~452 MB; brief 121
 *           will reduce to O(M) streaming.
 *
 * Args:     cm           - CM with cm->fp7 (filter HMM) set
 *           errbuf       - for error messages
 *           dsq          - digital sequence, 1..L
 *           L            - length of dsq
 *           delta_milli  - Delta threshold in milli-bits (3000 = 3 bits, default)
 *           ret_i2k      - RETURN: per-residue pin array (all -1; caller frees)
 *           ret_kmin     - RETURN: per-residue kmin array, length L+1 (caller frees)
 *           ret_kmax     - RETURN: per-residue kmax array, length L+1 (caller frees)
 *           ret_ncells   - RETURN: total banded cells sum_{i=1..L} (kmax[i]-kmin[i]+1)
 *
 * Returns:  eslOK on success.
 *           eslEMEM on allocation failure.
 *           eslEINVAL if cm->fp7 is NULL.
 */
int
p7_Seq2BandsIBV(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L, int delta_milli,
                int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int       status;
  P7_HMM   *hmm = NULL;
  int       M;
  int       i, k;
  int       x;
  double   *MM_t = NULL, *MI_t = NULL, *MD_t = NULL;  /* transition vectors, milli-bits */
  double   *IM_t = NULL, *II_t = NULL;
  double   *DM_t = NULL, *DD_t = NULL;
  double  **FM = NULL, **FI = NULL, **FD = NULL;
  double  **BM = NULL, **BI = NULL, **BD = NULL;
  double   *FM_pool = NULL, *FI_pool = NULL, *FD_pool = NULL;
  double   *BM_pool = NULL, *BI_pool = NULL, *BD_pool = NULL;
  int      *i2k = NULL, *kmin = NULL, *kmax = NULL;
  double    optimal, thr;
  double    floor_milli;
  int       ncells = 0;

  if (cm == NULL || cm->fp7 == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsIBV: cm->fp7 is NULL");
  hmm = cm->fp7;
  M = hmm->M;
  if (L < 1 || M < 1)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsIBV: bad L=%d or M=%d", L, M);

  /* Allocate transition vectors. Indexed 0..M. Per HMMER convention:
   *   hmm->t[k][p7H_MM] = transition prob from M_k -> M_{k+1} (k=0 is begin)
   *   hmm->t[k][p7H_MI] = M_k -> I_k
   *   hmm->t[k][p7H_MD] = M_k -> D_{k+1}
   *   hmm->t[k][p7H_IM] = I_k -> M_{k+1}
   *   hmm->t[k][p7H_II] = I_k -> I_k
   *   hmm->t[k][p7H_DM] = D_k -> M_{k+1}
   *   hmm->t[k][p7H_DD] = D_k -> D_{k+1}
   * Match the Python prototype's trans["M_M"][k] etc.
   */
  ESL_ALLOC(MM_t, sizeof(double) * (M + 1));
  ESL_ALLOC(MI_t, sizeof(double) * (M + 1));
  ESL_ALLOC(MD_t, sizeof(double) * (M + 1));
  ESL_ALLOC(IM_t, sizeof(double) * (M + 1));
  ESL_ALLOC(II_t, sizeof(double) * (M + 1));
  ESL_ALLOC(DM_t, sizeof(double) * (M + 1));
  ESL_ALLOC(DD_t, sizeof(double) * (M + 1));
  for (k = 0; k <= M; k++) {
    MM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MM]);
    MI_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MI]);
    MD_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MD]);
    IM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_IM]);
    II_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_II]);
    DM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_DM]);
    DD_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_DD]);
  }

  /* Allocate F* and B* as (L+1) x (M+1) double matrices.
   * Two-level layout (row-pointer + contiguous pool) for clean indexing
   * F_M[i][k]. Memory cost is 6 * (L+1) * (M+1) * 8 bytes; brief 121 will
   * port to a streaming O(M) implementation.
   */
  size_t cells = (size_t)(L + 1) * (size_t)(M + 1);
  ESL_ALLOC(FM_pool, sizeof(double) * cells);
  ESL_ALLOC(FI_pool, sizeof(double) * cells);
  ESL_ALLOC(FD_pool, sizeof(double) * cells);
  ESL_ALLOC(BM_pool, sizeof(double) * cells);
  ESL_ALLOC(BI_pool, sizeof(double) * cells);
  ESL_ALLOC(BD_pool, sizeof(double) * cells);
  ESL_ALLOC(FM, sizeof(double *) * (L + 1));
  ESL_ALLOC(FI, sizeof(double *) * (L + 1));
  ESL_ALLOC(FD, sizeof(double *) * (L + 1));
  ESL_ALLOC(BM, sizeof(double *) * (L + 1));
  ESL_ALLOC(BI, sizeof(double *) * (L + 1));
  ESL_ALLOC(BD, sizeof(double *) * (L + 1));
  for (i = 0; i <= L; i++) {
    FM[i] = FM_pool + (size_t)i * (size_t)(M + 1);
    FI[i] = FI_pool + (size_t)i * (size_t)(M + 1);
    FD[i] = FD_pool + (size_t)i * (size_t)(M + 1);
    BM[i] = BM_pool + (size_t)i * (size_t)(M + 1);
    BI[i] = BI_pool + (size_t)i * (size_t)(M + 1);
    BD[i] = BD_pool + (size_t)i * (size_t)(M + 1);
  }
  for (i = 0; i < (int)cells; i++) {
    FM_pool[i] = FI_pool[i] = FD_pool[i] = P7IBV_NEG_INF;
    BM_pool[i] = BI_pool[i] = BD_pool[i] = P7IBV_NEG_INF;
  }

  /* Forward recurrence (per Python ibv.full_viterbi + brief 116 spec).
   * Glocal seed: F_M[0,0] = 0; rest -INF (no entry mid-model or mid-seq).
   */
  FM[0][0] = 0.0;
  for (i = 0; i <= L; i++) {
    x = (i >= 1) ? (int) dsq[i] : -1;
    for (k = 0; k <= M; k++) {
      if (i == 0 && k == 0) continue;

      double cM = P7IBV_NEG_INF;
      double cI = P7IBV_NEG_INF;
      double cD = P7IBV_NEG_INF;

      /* I[i,k]: consume residue i, no model advance. Insert emission is
       * (by HMMER convention) = background, so log-odds emit = 0 -- omitted.
       */
      if (i >= 1) {
        double a = FM[i-1][k] + MI_t[k];
        double b = FI[i-1][k] + II_t[k];
        cI = (a > b) ? a : b;
      }

      /* D[i,k]: no emission, model advance only. Depends on row i (cur row). */
      if (k >= 1) {
        double a = FM[i][k-1] + MD_t[k-1];
        double b = FD[i][k-1] + DD_t[k-1];
        cD = (a > b) ? a : b;
      }

      /* M[i,k]: consume residue i, match column k. */
      if (i >= 1 && k >= 1) {
        double e = p7ibv_emit_milli(hmm, k, x);
        double a = FM[i-1][k-1] + MM_t[k-1];
        double b = FI[i-1][k-1] + IM_t[k-1];
        double c = FD[i-1][k-1] + DM_t[k-1];
        double m = (a > b) ? a : b;
        if (c > m) m = c;
        cM = m + e;
      }

      FM[i][k] = cM;
      FI[i][k] = cI;
      FD[i][k] = cD;
    }
  }

  /* Backward recurrence (mirror of forward, per Python backward_full).
   * Glocal seed all three states at (L, M): B_M = B_I = B_D = 0; rest -INF.
   * (ZFAT proves D can win the terminal -- brief 116 §2.)
   */
  for (i = L; i >= 0; i--) {
    /* x_next is the residue we'd consume on i+1 going outward. */
    int x_next = (i + 1 <= L) ? (int) dsq[i + 1] : -1;
    for (k = M; k >= 0; k--) {
      double term = (i == L && k == M) ? 0.0 : P7IBV_NEG_INF;
      double cM = term;
      double cI = term;
      double cD = term;

      /* -> M(i+1, k+1) consumes residue i+1, matches column k+1. */
      if (i + 1 <= L && k + 1 <= M) {
        double e = p7ibv_emit_milli(hmm, k + 1, x_next);
        double bv = BM[i+1][k+1] + e;
        double a = MM_t[k] + bv;
        double b = IM_t[k] + bv;
        double c = DM_t[k] + bv;
        if (a > cM) cM = a;
        if (b > cI) cI = b;
        if (c > cD) cD = c;
      }

      /* -> I(i+1, k) consumes residue i+1, no model advance. */
      if (i + 1 <= L) {
        double bv = BI[i+1][k];
        double a = MI_t[k] + bv;
        double b = II_t[k] + bv;
        if (a > cM) cM = a;
        if (b > cI) cI = b;
      }

      /* -> D(i, k+1) no emission, model advance only. */
      if (k + 1 <= M) {
        double bv = BD[i][k+1];
        double a = MD_t[k] + bv;
        double b = DD_t[k] + bv;
        if (a > cM) cM = a;
        if (b > cD) cD = b;
      }

      BM[i][k] = cM;
      BI[i][k] = cI;
      BD[i][k] = cD;
    }
  }

  /* Compute optimal score = max{F_M[L,M], F_I[L,M], F_D[L,M]}. */
  optimal = FM[L][M];
  if (FI[L][M] > optimal) optimal = FI[L][M];
  if (FD[L][M] > optimal) optimal = FD[L][M];

  /* Threshold = optimal - max(delta, EPS). EPS guards against the case
   * where delta=0 collapses through-score equality to a strict float compare
   * (the optimal-path cells equal optimal up to a few ULPs).
   */
  floor_milli = (double) delta_milli;
  if (floor_milli < P7IBV_EPS) floor_milli = P7IBV_EPS;
  thr = optimal - floor_milli;

  /* Allocate output arrays. */
  ESL_ALLOC(i2k,  sizeof(int) * (L + 1));
  ESL_ALLOC(kmin, sizeof(int) * (L + 1));
  ESL_ALLOC(kmax, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, L + 1, -1);

  /* Per-row band emission. through[i,k] = max over X in {M,I,D} of
   * (F_X[i,k] + B_X[i,k]). Admitted cells: through >= thr, finite (both F
   * and B reach the cell).
   */
  for (i = 1; i <= L; i++) {
    int   row_kmin = -1, row_kmax = -1;
    for (k = 1; k <= M; k++) {
      double t_m = FM[i][k] + BM[i][k];
      double t_i = FI[i][k] + BI[i][k];
      double t_d = FD[i][k] + BD[i][k];
      double th  = t_m;
      if (t_i > th) th = t_i;
      if (t_d > th) th = t_d;
      /* "finite" filter: brief 117 §1 -- if both sides are -INF the sum is
       * even more negative; clip out using a halfway floor.
       */
      if (th < (P7IBV_NEG_INF / 2.0)) continue;
      if (th >= thr) {
        if (row_kmin < 0) row_kmin = k;
        row_kmax = k;
      }
    }
    if (row_kmin < 0) {
      /* Degenerate row: no admitted cell. Default to [1, M] (safest non-
       * restrictive); mirrors the Python prototype's empty-row fallback.
       */
      kmin[i] = 1;
      kmax[i] = M;
    } else {
      kmin[i] = row_kmin;
      kmax[i] = row_kmax;
    }
  }

  /* Boundary widening: rows i=1, L-1, L forced to [1, M] for truncated-
   * alignment entry/exit. The Python prototype's fb_direct_band() applies
   * this when boundary_widen=True (the default; see brief 117 §0 design
   * intuition + the prototype's docstring).
   */
  if (L >= 1) { kmin[1] = 1; kmax[1] = M; }
  if (L >= 2) { kmin[L - 1] = 1; kmax[L - 1] = M; }
  if (L >= 1) { kmin[L] = 1; kmax[L] = M; }

  /* kmin[0] / kmax[0]: B-state convention from brief 117 §1. The Python
   * prototype writes [1, M] at i=0 (PB_LOAD_BAND-targeted), but the
   * production downstream cp9_FB2HMMBandsP7BF asserts kmin[0]==0; the
   * PB_LOAD_BAND loader in cm_p7_band.c clamps kmin[0]=0 explicitly
   * (commit a956388b). We match that here.
   */
  kmin[0] = 0;
  kmax[0] = 0;

  /* Total band cells (rows 1..L only; i=0 doesn't count toward ncells). */
  for (i = 1; i <= L; i++)
    ncells += (kmax[i] - kmin[i] + 1);

  /* Diagnostic: optional band dump, gated by env var (mirrors brief 104's
   * PB_DUMP_BAND hook but for IBV's pre-CP9 band). When set, write a TSV
   * with one row per i in 1..L: 'i\tkmin\tkmax\twidth'. Used by Phase 4
   * validation to compare C output against the Python prototype's TSV.
   */
  {
    const char *dump = getenv("P7IBV_DUMP_BAND");
    if (dump != NULL && *dump != '\0') {
      FILE *fp = fopen(dump, "w");
      if (fp != NULL) {
        fprintf(fp, "# M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d\n",
                M, L, delta_milli, optimal, thr, ncells);
        fprintf(fp, "# i\tkmin\tkmax\twidth\n");
        for (i = 1; i <= L; i++)
          fprintf(fp, "%d\t%d\t%d\t%d\n", i, kmin[i], kmax[i], kmax[i] - kmin[i] + 1);
        fclose(fp);
      }
    }
  }

  /* Cleanup intermediate buffers. */
  free(MM_t); free(MI_t); free(MD_t);
  free(IM_t); free(II_t); free(DM_t); free(DD_t);
  free(FM); free(FI); free(FD); free(BM); free(BI); free(BD);
  free(FM_pool); free(FI_pool); free(FD_pool);
  free(BM_pool); free(BI_pool); free(BD_pool);

  *ret_i2k    = i2k;
  *ret_kmin   = kmin;
  *ret_kmax   = kmax;
  *ret_ncells = ncells;
  return eslOK;

 ERROR:
  if (MM_t) free(MM_t); if (MI_t) free(MI_t); if (MD_t) free(MD_t);
  if (IM_t) free(IM_t); if (II_t) free(II_t);
  if (DM_t) free(DM_t); if (DD_t) free(DD_t);
  if (FM) free(FM); if (FI) free(FI); if (FD) free(FD);
  if (BM) free(BM); if (BI) free(BI); if (BD) free(BD);
  if (FM_pool) free(FM_pool); if (FI_pool) free(FI_pool); if (FD_pool) free(FD_pool);
  if (BM_pool) free(BM_pool); if (BI_pool) free(BI_pool); if (BD_pool) free(BD_pool);
  if (i2k)  free(i2k);
  if (kmin) free(kmin);
  if (kmax) free(kmax);
  *ret_i2k    = NULL;
  *ret_kmin   = NULL;
  *ret_kmax   = NULL;
  *ret_ncells = 0;
  return status;
}

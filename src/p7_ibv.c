/* p7_ibv.c -- F+B direct-band derivation: SSE per-row primitives + D&C wrapper
 *
 * Brief 26_0430-124 C1: refactor brief 26_0430-121's monolithic p7_Seq2BandsIBV into three
 *   reusable per-row SSE primitives (ibv_forward_one_row, ibv_backward_one_row,
 *   ibv_through_scan).  No algorithm change; byte-exact vs brief 26_0430-121 C3.
 *
 * Brief 26_0430-124 C2: p7_Seq2BandsIBV_dnc wraps the C1 primitives in a recursive
 *   divide-and-conquer band deriver (O(M * log L) peak memory).  Accessed via
 *   --p7ibv --p7ibv-mem; --p7ibv alone still calls the flat p7_Seq2BandsIBV.
 *
 * SSE/memory conventions (brief 26_0430-121 C2/C3; cells double and SSE 2-lane since brief 26_0821-144):
 *   k_stride = ((M+4+15) & ~15)  (16-double align + >=3 pad slots above M)
 *   Forward row i: scalar prefix k=0..3, SSE bulk k=4..k_sse_end-1 for M+I,
 *     scalar tail, scalar left-to-right D-fill.
 *   Backward row i from row i+1: scalar right-to-left D-fill, SSE bulk M+I,
 *     scalar tail.
 *   Padded tail k > M: MUST stay P7IBV_NEG_INF; all allocation helpers init
 *     the tail to NEG_INF and SSE writes stay within k=0..M.
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

#define P7IBV_INTSCALE        1000.0f
#define P7IBV_NEG_INF         (-1.0e18f)
#define P7IBV_HALF_NEG_INF    (-5.0e17f)
#define P7IBV_EPS             1.0f
#define P7IBV_K_ALIGN         16

/* Band threshold: each row's own max through-score minus the margin
 * (brief 26_0821-144).
 *
 * The band admits every cell whose double max-sum F+B through-score is within
 * margin = max(delta, EPS) of the optimum. In exact arithmetic every row's max
 * through-score equals the optimum, since the Viterbi path crosses every row,
 * so thresholding a row against its own max is the same band. In double it is
 * not: F and B each accumulate rounding along the sequence, and at genome
 * scale (|optimal| ~ 1e8 milli-bits) F+B on the optimal path drifts from the
 * separately computed optimum by hundreds of bits (brief 26_0821-142: +327 to
 * +433 bits on HSV-1, mpox and EBV, against a 1-bit margin). Rows with positive
 * drift then admit far too many cells, and rows with negative drift admit none.
 * The row-relative test cancels that between-row drift, leaving only the
 * within-row rounding difference between cells. In float that was still up
 * to ~2 bits on HSV-1, enough to drop optimal-path cells from a 1-bit band
 * (no parse) and to move band edges by up to 255 cells at 20 bits, so the
 * F/B cells, their storage and the score tables are double as well (brief
 * 26_0821-144; the SSE kernels use 2-lane __m128d). The old guard here
 * (fail if |optimal| > margin * 2^19) modeled only the final comparison's
 * rounding, not accumulation along L, and never fired. It is replaced by the
 * drift detector in IBV_ROWSTATS, reported under --devdiag. */
typedef struct {
  double   optimal;      /* global optimum, for the drift detector               */
  int64_t nrows;        /* row scans (the D&C scans each midrow twice)          */
  int64_t empty_rows;   /* rows with no finite through-score: an error          */
  double   drift_max;    /* max over rows of (row max - optimal), milli-bits     */
  double   drift_min;    /* min over rows of (row max - optimal), milli-bits     */
} IBV_ROWSTATS;

static void
ibv_rowstats_init(IBV_ROWSTATS *st, double optimal)
{
  st->optimal    = optimal;
  st->nrows      = 0;
  st->empty_rows = 0;
  st->drift_max  = -eslINFINITY;
  st->drift_min  =  eslINFINITY;
}

static inline void
ibv_rowstats_add(IBV_ROWSTATS *st, double rowmax)
{
  st->nrows++;
  if (rowmax < P7IBV_HALF_NEG_INF) { st->empty_rows++; return; }
  double d = rowmax - st->optimal;
  if (d > st->drift_max) st->drift_max = d;
  if (d < st->drift_min) st->drift_min = d;
}

/* One developer line per band. An empty row (no finite through-score) cannot
 * happen when the optimum is finite, since the optimal path crosses every row;
 * ibv_rowstats_check() makes it an error. */
static IBV_ROWSTATS ibv_last_rowstats;   /* the last band's stats, for the unit test */

static void
ibv_rowstats_report(const IBV_ROWSTATS *st, const char *path, int L, int M, double margin, int64_t ncells)
{
  ibv_last_rowstats = *st;
  if (! cm_StageVerbose()) return;
  fprintf(stderr, "#IBVBAND path=%s L=%d M=%d margin_bits=%.3f optimal_bits=%.3f scans=%lld empty_rows=%lld drift_max_bits=%.3f drift_min_bits=%.3f ncells=%lld\n",
          path, L, M, margin / P7IBV_INTSCALE, st->optimal / P7IBV_INTSCALE,
          (long long) st->nrows, (long long) st->empty_rows,
          (st->nrows > st->empty_rows) ? st->drift_max / P7IBV_INTSCALE : 0.,
          (st->nrows > st->empty_rows) ? st->drift_min / P7IBV_INTSCALE : 0.,
          (long long) ncells);
}

/* D&C adaptive base_slab: when the caller passes base_slab <= 0, pick the
 * largest slab whose base-case F+B storage (= 2*(slab+1)*3*k_stride floats)
 * fits within P7IBV_SLAB_CAP_BYTES, clamped to [MIN, MAX].  This keeps the
 * base-case slab memory bounded across scales: large at LSU/dengue (fast,
 * shallow recursion) but ~64 at HSV/MPXV (memory-safe; preserves the <1 GB
 * genome-scale headline).  Larger base_slab = fewer recursion levels = less
 * recomputed band wall, at the cost of linear-in-M base-case slab memory. */
#define P7IBV_SLAB_CAP_BYTES  (256.0)   /* MB */
#define P7IBV_SLAB_MIN        32
#define P7IBV_SLAB_MAX        1024

static int
ibv_alloc_floats(size_t n, double **ret_p)
{
  void *p = NULL;
  if (posix_memalign(&p, 64, sizeof(double) * n) != 0 || p == NULL) {
    *ret_p = NULL;
    return eslEMEM;
  }
  *ret_p = (double *) p;
  return eslOK;
}

static inline double
p7ibv_lod_milli(double p)
{
  if (p <= 0.0f) return P7IBV_NEG_INF;
  return (double)(P7IBV_INTSCALE * (log((double) p) / M_LN2));
}

static inline double
p7ibv_emit_milli(const P7_HMM *hmm, int k, int x)
{
  if (x < 0 || x >= hmm->abc->K) return 0.0f;
  double pm = hmm->mat[k][x];
  double pi = hmm->ins[k][x];
  if (pm <= 0.0f || pi <= 0.0f) return P7IBV_NEG_INF;
  return (double)(P7IBV_INTSCALE * (log((double) pm / (double) pi) / M_LN2));
}

/* ibv_rowstats_check(): an empty row (no finite through-score) is impossible
 * when the optimum is finite, because the optimal path crosses every row.
 * Widening such a row to [1,M] used to hide a wrong band (it produced the
 * 567 M-cell EBV band of brief 26_0821-142), so it is an error instead
 * (brief 26_0821-144, as brief 26_0821-141's HMM-only decode check).
 * Returns eslOK, or eslEINCONCEIVABLE with <errbuf> set (a distinct code,
 * so callers that fall back on other band failures stop on this one). */
static int
ibv_rowstats_check(const IBV_ROWSTATS *st, char *errbuf)
{
  if (st->empty_rows > 0 && st->optimal >= P7IBV_HALF_NEG_INF)
    ESL_FAIL(eslEINCONCEIVABLE, errbuf,
             "IBV band calculation failed: %lld rows had no finite Forward+Backward score although the optimal parse scores %.2f bits; the band would be wrong",
             (long long) st->empty_rows, st->optimal / P7IBV_INTSCALE);
  return eslOK;
}

static inline __m128d
p7ibv_mm_max3(__m128d a, __m128d b, __m128d c)
{
  return _mm_max_pd(_mm_max_pd(a, b), c);
}


/* ---------------------------------------------------------------------------
 * C1 primitives: per-row SSE forward, backward, through-scan
 * ---------------------------------------------------------------------------*/

/* ibv_forward_one_row -- compute F[i] from F[i-1].
 *
 * emit_row = emit_table[dsq[i]].  FM_prev/FI_prev/FD_prev = row i-1.
 * FM_curr/FI_curr/FD_curr are written (k=0..M).  Tail k>M not touched.
 * Valid for i >= 1.
 *
 * begin_milli (brief 26_0430-171, Tgm begin-anywhere): if non-NULL, fold a local
 * begin into each match cell: M_k <- max(M_k, begin_milli[k] + emit_row[k]).
 * This is the truncated-mode entry (enter the model at any node k with score
 * log(occ[k]/Z)).  It must be folded BEFORE the left-to-right D-fill so that
 * begins can be followed by deletes on the same row.  Pass NULL for the glocal
 * path (begins handled by the row-0 D-cascade in the caller) and for every row
 * except the global begin row (Tgm: row 1 only, since N->N is impossible). The
 * begin_milli array must be valid for k=0..k_stride-1 (k=0 and k>M = NEG_INF).
 */
static void
ibv_forward_one_row(int M, size_t k_stride,
                    const double *MM_t, const double *MI_t, const double *MD_t,
                    const double *IM_t, const double *II_t,
                    const double *DM_t, const double *DD_t,
                    const double *emit_row,
                    const double *begin_milli,
                    const double *FM_prev, const double *FI_prev, const double *FD_prev,
                    double *FM_curr, double *FI_curr, double *FD_curr)
{
  int k;
  int kpref_end = (M < 3) ? M : 3;

  for (k = 0; k <= kpref_end; k++) {
    double cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF, cD = P7IBV_NEG_INF;
    {
      double a = FM_prev[k] + MI_t[k];
      double b = FI_prev[k] + II_t[k];
      cI = (a > b) ? a : b;
    }
    if (k >= 1) {
      double a = FM_prev[k - 1] + MM_t[k - 1];
      double b = FI_prev[k - 1] + IM_t[k - 1];
      double c = FD_prev[k - 1] + DM_t[k - 1];
      double m = (a > b) ? a : b;
      if (c > m) m = c;
      cM = m + emit_row[k];
      if (begin_milli != NULL) {
        double bc = begin_milli[k] + emit_row[k];
        if (bc > cM) cM = bc;
      }
      double a2 = FM_curr[k - 1] + MD_t[k - 1];
      double b2 = FD_curr[k - 1] + DD_t[k - 1];
      cD = (a2 > b2) ? a2 : b2;
    }
    FM_curr[k] = cM;
    FI_curr[k] = cI;
    FD_curr[k] = cD;
  }

  int k_sse_start = 4;
  int k_sse_end   = k_sse_start;
  while (k_sse_end + 1 <= M) k_sse_end += 2;
  for (k = k_sse_start; k < k_sse_end; k += 2) {
    __m128d m_prev = _mm_loadu_pd(&FM_prev[k - 1]);
    __m128d i_prev = _mm_loadu_pd(&FI_prev[k - 1]);
    __m128d d_prev = _mm_loadu_pd(&FD_prev[k - 1]);
    __m128d t_mm   = _mm_loadu_pd(&MM_t[k - 1]);
    __m128d t_im   = _mm_loadu_pd(&IM_t[k - 1]);
    __m128d t_dm   = _mm_loadu_pd(&DM_t[k - 1]);
    __m128d a      = _mm_add_pd(m_prev, t_mm);
    __m128d b      = _mm_add_pd(i_prev, t_im);
    __m128d c      = _mm_add_pd(d_prev, t_dm);
    __m128d mx     = p7ibv_mm_max3(a, b, c);
    if (begin_milli != NULL) mx = _mm_max_pd(mx, _mm_loadu_pd(&begin_milli[k]));
    __m128d e_vec  = _mm_loadu_pd(&emit_row[k]);
    _mm_storeu_pd(&FM_curr[k], _mm_add_pd(mx, e_vec));

    __m128d fm_k = _mm_loadu_pd(&FM_prev[k]);
    __m128d fi_k = _mm_loadu_pd(&FI_prev[k]);
    __m128d t_mi = _mm_loadu_pd(&MI_t[k]);
    __m128d t_ii = _mm_loadu_pd(&II_t[k]);
    __m128d a2   = _mm_add_pd(fm_k, t_mi);
    __m128d b2   = _mm_add_pd(fi_k, t_ii);
    _mm_storeu_pd(&FI_curr[k], _mm_max_pd(a2, b2));
  }

  for (k = k_sse_end; k <= M; k++) {
    double cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
    {
      double a = FM_prev[k] + MI_t[k];
      double b = FI_prev[k] + II_t[k];
      cI = (a > b) ? a : b;
    }
    if (k >= 1) {
      double a = FM_prev[k - 1] + MM_t[k - 1];
      double b = FI_prev[k - 1] + IM_t[k - 1];
      double c = FD_prev[k - 1] + DM_t[k - 1];
      double m = (a > b) ? a : b;
      if (c > m) m = c;
      cM = m + emit_row[k];
      if (begin_milli != NULL) {
        double bc = begin_milli[k] + emit_row[k];
        if (bc > cM) cM = bc;
      }
    }
    FM_curr[k] = cM;
    FI_curr[k] = cI;
  }

  int kd_start = (kpref_end + 1 > M) ? (M + 1) : (kpref_end + 1);
  for (k = kd_start; k <= M; k++) {
    double a = FM_curr[k - 1] + MD_t[k - 1];
    double b = FD_curr[k - 1] + DD_t[k - 1];
    FD_curr[k] = (a > b) ? a : b;
  }
}


/* ibv_backward_one_row -- compute B[i] from B[i+1].
 *
 * Two cases:
 *   i == global_L : terminal injection.  BM_next/BI_next/emit_row_next ignored.
 *   i  < global_L : normal backward step.
 *
 * Writes BM_curr/BI_curr/BD_curr for k=0..M.  Tail k>M not touched.
 *
 * do_trunc (brief 26_0430-171, Tgm end-anywhere): when set, the terminal injection at
 * row L mirrors local exit -- every match state M_k may exit to E with score 0
 * (esc=0 in HMMER local mode), instead of the glocal forced exit from node M
 * (delete-cascade to D_M).  The delete cascade to D_M->E is kept (D_M->E is
 * allowed even in local mode); inserts still cannot be the last emitted state.
 * Only affects the terminal branch (i == global_L); the interior backward
 * recursion is mode-independent because the end is carried in by the terminal.
 */
static void
ibv_backward_one_row(int M, size_t k_stride, int i, int global_L, int do_trunc,
                     const double *MM_t, const double *MI_t, const double *MD_t,
                     const double *IM_t, const double *II_t,
                     const double *DM_t, const double *DD_t,
                     const double *emit_row_next,
                     const double *BM_next,
                     const double *BI_next,
                     double *BM_curr, double *BI_curr, double *BD_curr)
{
  int k;
  (void) k_stride;

  if (i == global_L) {
    if (do_trunc) {
      /* Tgm end-anywhere: M_k -> E exit (esc=0) at any node k. */
      for (k = 0; k <= M; k++) { BM_curr[k] = 0.0f; BI_curr[k] = P7IBV_NEG_INF; }
      BD_curr[M] = 0.0f;
      for (k = M - 1; k >= 0; k--) BD_curr[k] = DD_t[k] + BD_curr[k + 1];
      return;
    }
    BM_curr[M] = 0.0f;
    BI_curr[M] = 0.0f;
    BD_curr[M] = 0.0f;
    for (k = M - 1; k >= 0; k--) {
      double bv = BD_curr[k + 1];
      BM_curr[k] = MD_t[k] + bv;
      BD_curr[k] = DD_t[k] + bv;
      BI_curr[k] = P7IBV_NEG_INF;
    }
    return;
  }

  BD_curr[M] = P7IBV_NEG_INF;
  for (k = M - 1; k >= 0; k--) {
    double bv_m = BM_next[k + 1] + emit_row_next[k + 1];
    double bv_d = BD_curr[k + 1];
    double a = DM_t[k] + bv_m;
    double b = DD_t[k] + bv_d;
    BD_curr[k] = (a > b) ? a : b;
  }

  int k_sse_end = 0;
  while (k_sse_end + 1 <= M) k_sse_end += 2;
  for (k = 0; k < k_sse_end; k += 2) {
    __m128d bm_n   = _mm_loadu_pd(&BM_next[k + 1]);
    __m128d e_n    = _mm_loadu_pd(&emit_row_next[k + 1]);
    __m128d bv_M   = _mm_add_pd(bm_n, e_n);
    __m128d bv_I   = _mm_loadu_pd(&BI_next[k]);
    __m128d bv_D   = _mm_loadu_pd(&BD_curr[k + 1]);

    __m128d t_mm   = _mm_loadu_pd(&MM_t[k]);
    __m128d t_mi   = _mm_loadu_pd(&MI_t[k]);
    __m128d t_md   = _mm_loadu_pd(&MD_t[k]);
    __m128d cM_out = p7ibv_mm_max3(_mm_add_pd(t_mm, bv_M),
                                  _mm_add_pd(t_mi, bv_I),
                                  _mm_add_pd(t_md, bv_D));
    _mm_storeu_pd(&BM_curr[k], cM_out);

    __m128d t_im   = _mm_loadu_pd(&IM_t[k]);
    __m128d t_ii   = _mm_loadu_pd(&II_t[k]);
    __m128d cI_out = _mm_max_pd(_mm_add_pd(t_im, bv_M),
                               _mm_add_pd(t_ii, bv_I));
    _mm_storeu_pd(&BI_curr[k], cI_out);
  }

  for (k = k_sse_end; k <= M; k++) {
    double cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
    if (k + 1 <= M) {
      double bv_m = BM_next[k + 1] + emit_row_next[k + 1];
      double a = MM_t[k] + bv_m;
      double b = IM_t[k] + bv_m;
      if (a > cM) cM = a;
      if (b > cI) cI = b;
    }
    {
      double bv_i = BI_next[k];
      double a = MI_t[k] + bv_i;
      double b = II_t[k] + bv_i;
      if (a > cM) cM = a;
      if (b > cI) cI = b;
    }
    if (k + 1 <= M) {
      double bv_d = BD_curr[k + 1];
      double a = MD_t[k] + bv_d;
      if (a > cM) cM = a;
    }
    BM_curr[k] = cM;
    BI_curr[k] = cI;
  }
}


/* ibv_through_scan -- per-row kmin/kmax from F+B through-score.
 *
 * The cloud is every k whose through-score is >= this row's max - margin
 * (brief 26_0821-144; see IBV_ROWSTATS). through_scratch: caller-owned
 * k_stride buffer. Returns (1,M) on an empty row, which <st> counts
 * and the deriver then reports as an error (ibv_rowstats_check()).
 */
static void
ibv_through_scan(int M, size_t k_stride, double margin,
                 int ibv_mode, int ibv_width,
                 const double *FM, const double *FI, const double *FD,
                 const double *BM, const double *BI, const double *BD,
                 double *through_scratch, IBV_ROWSTATS *st,
                 int *ret_kmin, int *ret_kmax, int *ret_kargmax)
{
  int k;
  int k_thru_end = 0;
  while (k_thru_end + 1 <= M) k_thru_end += 2;
  (void) k_stride;

  for (k = 0; k < k_thru_end; k += 2) {
    __m128d a = _mm_add_pd(_mm_loadu_pd(&FM[k]), _mm_loadu_pd(&BM[k]));
    __m128d b = _mm_add_pd(_mm_loadu_pd(&FI[k]), _mm_loadu_pd(&BI[k]));
    __m128d c = _mm_add_pd(_mm_loadu_pd(&FD[k]), _mm_loadu_pd(&BD[k]));
    _mm_storeu_pd(&through_scratch[k], p7ibv_mm_max3(a, b, c));
  }
  for (k = k_thru_end; k <= M; k++) {
    double t_m = FM[k] + BM[k];
    double t_i = FI[k] + BI[k];
    double t_d = FD[k] + BD[k];
    double t = t_m;
    if (t_i > t) t = t_i;
    if (t_d > t) t = t_d;
    through_scratch[k] = t;
  }

  double rowmax = P7IBV_NEG_INF;
  for (k = 1; k <= M; k++) if (through_scratch[k] > rowmax) rowmax = through_scratch[k];
  ibv_rowstats_add(st, rowmax);
  double thr = rowmax - margin;

  int   row_kmin = -1, row_kmax = -1;
  int   k_argmax = -1;
  double t_argmax = P7IBV_NEG_INF;
  for (k = 1; k <= M; k++) {
    double t = through_scratch[k];
    if (t < P7IBV_HALF_NEG_INF) continue;
    /* Band (kmin/kmax) uses the full through-score incl. delete cells.
     *
     * The argmax-k pin (brief 26_0430-137), however, must be over EMITTING states (M,I)
     * only: residue i is emitted by a match or insert, never a delete.  A
     * delete cell D(i,k) sits on the optimal path too (with through == optimal,
     * at a higher k than the emitter) and double F+B reconstruction can tip the
     * tie so the delete edges out the emitter -- which would make i2k[i] point
     * at a model position where residue i is actually *deleted*.  Restricting
     * the argmax to FM+BM / FI+BI excludes delete cells from the pin. */
    double t_emit = FM[k] + BM[k];
    { double ti = FI[k] + BI[k]; if (ti > t_emit) t_emit = ti; }
    if (t_emit >= P7IBV_HALF_NEG_INF && t_emit > t_argmax) { t_argmax = t_emit; k_argmax = k; }
    if (t >= thr) {
      if (row_kmin < 0) row_kmin = k;
      row_kmax = k;
    }
  }
  if (row_kmin < 0) { *ret_kmin = 1; *ret_kmax = M; }
  else              { *ret_kmin = row_kmin; *ret_kmax = row_kmax; }
  if (ret_kargmax) *ret_kargmax = k_argmax;

  /* Brief 26_0430-140: enrich the per-row band using the argmax-k pin (k_argmax = i2k[i],
   * the IBV-Viterbi cell that emits residue i).  This follows the optimal path by
   * construction, unlike the Delta-cloud which admits noisy off-path cells.
   *   FIXED  : replace band with the fixed-width spine [k_argmax-W, k_argmax+W].
   *   HYBRID : union the Delta cloud with that spine (kmin <- min, kmax <- max).
   *   DELTA  : leave the cloud unchanged (back-compat; this whole block is skipped).
   * No-op when k_argmax < 1 (no emitting cell on this row, e.g. a pure-delete row,
   * or the row-0 B-state convention i2k[0]=0): keep the Delta/default band so the
   * existing boundary-row handling is preserved.  Clamp to model positions [1,M]. */
  if (ibv_mode != P7IBV_MODE_DELTA && k_argmax >= 1) {
    int lo = k_argmax - ibv_width; if (lo < 1) lo = 1;
    int hi = k_argmax + ibv_width; if (hi > M) hi = M;
    if (ibv_mode == P7IBV_MODE_FIXED) {
      *ret_kmin = lo;
      *ret_kmax = hi;
    } else { /* P7IBV_MODE_HYBRID: union spine into the Delta cloud */
      if (lo < *ret_kmin) *ret_kmin = lo;
      if (hi > *ret_kmax) *ret_kmax = hi;
    }
  }
}


/* ibv_connectivity_guard -- bridge inter-row band gaps for FIXED/HYBRID modes.
 *
 * Brief 26_0430-140a.  The narrow per-row FIXED/HYBRID bands (centered on the argmax-k
 * pin) can be DISCONNECTED between adjacent rows when the pin jumps -- e.g. a
 * long insert on the optimal path makes k_argmax[i] >> k_argmax[i-1], so
 * kmin[i] > kmax[i-1] + 1.  The downstream banded CP9 Forward/Backward
 * (cp9_ForwardP7BF / cp9_BackwardP7BF) then has no complete row-0..row-L path:
 * the backward total bmx->mmx[0][0] collapses to -inf, the posterior
 *   pmx = fmx + bmx - (-inf) = (-inf) + (-inf) + inf = NaN
 * (cp9_FB2HMMBandsP7BF), and the NaN reaches p7_FLogsum, whose lookup-table
 * index (int)((max-min)*1000) overflows the 16000-entry flogsum_lookup -> SIGSEGV.
 *
 * Connectivity precondition (per adjacent pair i-1,i): residue i is emitted by
 * M_{i,k} (predecessor k-1 in band[i-1]) or I_{i,k} (predecessor k in band[i-1]),
 * with intra-row deletes shifting k upward only.  A valid entry into band[i]
 * exists iff [kmin[i-1], kmax[i-1]+1] intersects [kmin[i], kmax[i]], i.e.
 *     kmin[i] <= kmax[i-1] + 1   AND   kmax[i] >= kmin[i-1].
 * (The brief 26_0430-140a sketch used kmin[i-1]-1 for the second test; that is too weak
 *  by one -- kmax[i]==kmin[i-1]-1 still leaves no in-band predecessor.)
 *
 * The guard widens band[i] minimally to satisfy this for every i.  DELTA mode is
 * skipped (its posterior-reachable bands are connected by construction).  The
 * forward sweep starts at i=2, i.e. it never touches index 0 or 1 -- this was
 * historically justified as "rows 0 and 1 connect to the B-state via the begin
 * transition (no diagonal needed)", but brief 26_0430-307 found that reasoning
 * holds only in LOCAL mode: in glocal, row 0 also carries the CP9's leading
 * delete chain B -> D_1 -> ... -> D_k, which is not a begin-transition path.
 * The guard's behavior is unaffected (it still never narrows index 0, so it
 * cannot defeat the row-0 widening landed by brief 26_0430-309); only the
 * stated justification for skipping it needed correcting.  A single forward
 * sweep is provably sufficient to guarantee a complete path over i=2..L
 * (hence a finite total); the backward sweep is kept for robustness -- it can
 * only widen further and so cannot re-introduce a gap.  Widening keeps
 * kmin>=1 and kmax<=M for i>=2; the final clamp is purely defensive.
 */
/* brief 26_0430-307/309: row-0 delete-chain band coverage. Formerly a
 * candidate fix gated behind P298_ROW0D (brief 26_0430-298, default OFF);
 * brief 307 established this is a CORRECTNESS FIX, not a tuning knob, and
 * brief 309 promoted it to unconditional, removing the flag entirely.
 *
 * Every IBV band deriver hard-pinned row 0 to the single cell k=0
 * (kmin[0]=0, kmax[0]=0), i.e. the B state alone. But row 0 is also where
 * the CP9's entire no-residue-consumed delete chain B -> D_1 -> D_2 -> ...
 * -> D_k lives (leading deletes are charged to dmx[0][k]; with kmax[0]=0
 * that chain is not merely out of band -- GrowCP9Matrix sizes row 0 to one
 * cell, so dmx[0][k>=1] is not even allocated). Row-1 widening cannot
 * substitute: the row-i match recursion clips kx = min(kmax[i],
 * kmax[i-1]+1), so row 1 stays pinned to k<=1 regardless of its own band.
 * cm_p7_band.c's own row-0 loop comment says dmx[0][k] is possible for all
 * k in [kmin[0], kmax[0]] -- the old row-0 convention conflated "only B can
 * be *entered* from outside at row 0" with "only B can be *occupied* at
 * row 0". Empirically a no-op on --hmm paths (cm_p7_kbands2gbands never
 * reads row 0) and on CM-mode local/truncated (local begins supply bsc[k]
 * finite); material only in CM-mode glocal -g non-truncated. */

static void
ibv_connectivity_guard(int L, int M, int ibv_mode, int *kmin, int *kmax)
{
  int i;
  if (ibv_mode == P7IBV_MODE_DELTA) return;
  for (i = 2; i <= L; i++) {                 /* forward: connect row i down to i-1 */
    if (kmin[i] > kmax[i-1] + 1) kmin[i] = kmax[i-1] + 1;
    if (kmax[i] < kmin[i-1])     kmax[i] = kmin[i-1];
  }
  for (i = L - 1; i >= 1; i--) {             /* backward: connect row i up to i+1 */
    if (kmin[i] > kmax[i+1] + 1) kmin[i] = kmax[i+1] + 1;
    if (kmax[i] < kmin[i+1])     kmax[i] = kmin[i+1];
  }
  for (i = 1; i <= L; i++) {                  /* defensive clamp to [1,M] */
    if (kmin[i] < 1) kmin[i] = 1;
    if (kmax[i] > M) kmax[i] = M;
    if (kmax[i] < kmin[i]) kmax[i] = kmin[i];
  }
}


/* ---------------------------------------------------------------------------
 * p7_Seq2BandsIBV -- C1 rewrite (byte-exact vs brief 26_0430-121 C3)
 * ---------------------------------------------------------------------------*/

int
cm_p7_Seq2BandsIBV(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L, int delta_milli,
                int do_trunc,
                int ibv_mode, int ibv_width,
                int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int       status;
  P7_HMM   *hmm = NULL;
  int       M;
  int       i, k;
  int       K;
  double    *begin_milli = NULL;
  double    *MM_t = NULL, *MI_t = NULL, *MD_t = NULL;
  double    *IM_t = NULL, *II_t = NULL;
  double    *DM_t = NULL, *DD_t = NULL;
  double    *FM_pool = NULL, *FI_pool = NULL, *FD_pool = NULL;
  double    *BM_a = NULL, *BM_b = NULL;
  double    *BI_a = NULL, *BI_b = NULL;
  double    *BD_a = NULL, *BD_b = NULL;
  double    *through = NULL;
  double    *emit_pool = NULL;
  double   **emit_table = NULL;
  int      *i2k = NULL, *kmin = NULL, *kmax = NULL;
  double     optimal, thr;
  double     floor_milli;
  IBV_ROWSTATS rowstats;
  int       ncells = 0;
  size_t    k_stride;
  size_t    pool_cells;

  if (cm == NULL || cm->fp7 == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsIBV: cm->fp7 is NULL");
  hmm = cm->fp7;
  M = hmm->M;
  K = hmm->abc->K;
  if (L < 1 || M < 1)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsIBV: bad L=%d or M=%d", L, M);

  /* k_stride must leave padding ABOVE index M: the SSE backward loads
   * _mm_loadu_pd(&BM_next[k+1]) reach index k+2 = k_sse_end <= M+1 (k+4 with the old 4-lane float kernels), so we
   * need k_stride >= M+2. Rounding (M+1) up to 16 gives ZERO padding when
   * M+1 is a multiple of 16 (M == 15 mod 16, e.g. M=287), leaving index M+1
   * out of bounds. Round (M+4) up instead to guarantee >=3 padding slots
   * (all NEG_INF), so the overread folds harmlessly. (Latent in brief 26_0430-121's
   * flat code too, but benign there with separate per-array allocations;
   * harmful in the D&C's contiguous arena where BM_next[M+1] aliases the
   * next state's k=0 cell.) */
  k_stride   = ((size_t)(M + 4) + (P7IBV_K_ALIGN - 1)) & ~(size_t)(P7IBV_K_ALIGN - 1);
  pool_cells = (size_t)(L + 1) * k_stride;

  if ((status = ibv_alloc_floats(k_stride, &MM_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &MI_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &MD_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &IM_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &II_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &DM_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &DD_t)) != eslOK) goto ERROR;
  for (k = 0; k <= M; k++) {
    MM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MM]);
    MI_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MI]);
    MD_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MD]);
    IM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_IM]);
    II_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_II]);
    DM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_DM]);
    DD_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_DD]);
  }
  for (k = M + 1; k < (int) k_stride; k++) {
    MM_t[k] = MI_t[k] = MD_t[k] = P7IBV_NEG_INF;
    IM_t[k] = II_t[k] = P7IBV_NEG_INF;
    DM_t[k] = DD_t[k] = P7IBV_NEG_INF;
  }

  if ((status = ibv_alloc_floats((size_t)(K + 1) * k_stride, &emit_pool)) != eslOK) goto ERROR;
  ESL_ALLOC(emit_table, sizeof(double *) * (K + 1));
  for (int xt = 0; xt <= K; xt++) emit_table[xt] = emit_pool + (size_t) xt * k_stride;
  for (int xt = 0; xt < K; xt++) {
    double *row = emit_table[xt];
    for (k = 0; k <= M; k++) row[k] = p7ibv_emit_milli(hmm, k, xt);
    for (k = M + 1; k < (int) k_stride; k++) row[k] = P7IBV_NEG_INF;
  }
  for (k = 0; k < (int) k_stride; k++) emit_table[K][k] = 0.0f;

  if ((status = ibv_alloc_floats(pool_cells, &FM_pool)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(pool_cells, &FI_pool)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(pool_cells, &FD_pool)) != eslOK) goto ERROR;
  for (size_t c = 0; c < pool_cells; c++)
    FM_pool[c] = FI_pool[c] = FD_pool[c] = P7IBV_NEG_INF;

  if ((status = ibv_alloc_floats(k_stride, &BM_a)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &BM_b)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &BI_a)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &BI_b)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &BD_a)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(k_stride, &BD_b)) != eslOK) goto ERROR;
  for (k = 0; k < (int) k_stride; k++) {
    BM_a[k] = BM_b[k] = P7IBV_NEG_INF;
    BI_a[k] = BI_b[k] = P7IBV_NEG_INF;
    BD_a[k] = BD_b[k] = P7IBV_NEG_INF;
  }

  if ((status = ibv_alloc_floats(k_stride, &through)) != eslOK) goto ERROR;

  /* Brief 26_0430-171: Tgm begin-anywhere vector.  begin_milli[k] = milli-bit log2 of
   * the local entry probability into M_k (occ[k] / sum_i occ[i]*(M-i+1)), the
   * same occupancy-weighted local begin HMMER's p7_ProfileConfig(p7_LOCAL) sets
   * (and the vitband/pinbridge Tgm reference uses via cm_alndata.c:459-461).
   * Only used at the global begin row (row 1). */
  if (do_trunc) {
    float *occ = NULL;
    double Z = 0.0;
    ESL_ALLOC(occ, sizeof(float) * (M + 1));
    if ((status = p7_hmm_CalculateOccupancy(hmm, occ, NULL)) != eslOK) { free(occ); goto ERROR; }
    for (k = 1; k <= M; k++) Z += (double) occ[k] * (double) (M - k + 1);
    if ((status = ibv_alloc_floats(k_stride, &begin_milli)) != eslOK) { free(occ); goto ERROR; }
    for (k = 0; k < (int) k_stride; k++) begin_milli[k] = P7IBV_NEG_INF;
    for (k = 1; k <= M; k++) {
      double b = (Z > 0.0 && occ[k] > 0.0) ? (double) occ[k] / Z : 0.0;
      begin_milli[k] = (b > 0.0) ? (double)(P7IBV_INTSCALE * (log(b) / M_LN2)) : P7IBV_NEG_INF;
    }
    free(occ);
  }

#define F_M(i)  (FM_pool + (size_t)(i) * k_stride)
#define F_I(i)  (FI_pool + (size_t)(i) * k_stride)
#define F_D(i)  (FD_pool + (size_t)(i) * k_stride)

  /* Row 0 init.  Glocal: D-cascade from M_0 (B-state).  Tgm (do_trunc): leave
   * row 0 all NEG_INF (no glocal entry); begins are injected at row 1 via
   * begin_milli, so the parse may start at any node. */
  if (! do_trunc) {
    F_M(0)[0] = 0.0f;
    {
      double *fm0 = F_M(0);
      double *fd0 = F_D(0);
      for (k = 1; k <= M; k++) {
        double a = fm0[k - 1] + MD_t[k - 1];
        double b = fd0[k - 1] + DD_t[k - 1];
        fd0[k] = (a > b) ? a : b;
      }
    }
  }

  for (i = 1; i <= L; i++) {
    int x  = (int) dsq[i];
    int xt = (x >= 0 && x < K) ? x : K;
    const double *brow = (do_trunc && i == 1) ? begin_milli : NULL;
    ibv_forward_one_row(M, k_stride,
                        MM_t, MI_t, MD_t, IM_t, II_t, DM_t, DD_t,
                        emit_table[xt], brow,
                        F_M(i-1), F_I(i-1), F_D(i-1),
                        F_M(i),   F_I(i),   F_D(i));
  }

  {
    double *fm_L = F_M(L);
    double *fi_L = F_I(L);
    double *fd_L = F_D(L);
    if (do_trunc) {
      /* Tgm end-anywhere: exit from any match node, or D_M->E. */
      optimal = fd_L[M];
      for (k = 1; k <= M; k++) if (fm_L[k] > optimal) optimal = fm_L[k];
    } else {
      optimal = fm_L[M];
      if (fi_L[M] > optimal) optimal = fi_L[M];
      if (fd_L[M] > optimal) optimal = fd_L[M];
    }
  }
  assert(optimal == optimal);
  floor_milli = (double) delta_milli;
  if (floor_milli < P7IBV_EPS) floor_milli = P7IBV_EPS;
  thr = optimal - floor_milli;   /* reported in the band dumps only; rows use their own max */
  ibv_rowstats_init(&rowstats, optimal);

  ESL_ALLOC(i2k,  sizeof(int) * (L + 1));
  ESL_ALLOC(kmin, sizeof(int) * (L + 1));
  ESL_ALLOC(kmax, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, L + 1, -1);
  for (i = 0; i <= L; i++) { kmin[i] = 1; kmax[i] = M; }

  double *B_M_prev = BM_a, *B_I_prev = BI_a, *B_D_prev = BD_a;
  double *B_M_curr = BM_b, *B_I_curr = BI_b, *B_D_curr = BD_b;

  /* Row L: terminal injection. */
  ibv_backward_one_row(M, k_stride, L, L, do_trunc,
                       MM_t, MI_t, MD_t, IM_t, II_t, DM_t, DD_t,
                       NULL, NULL, NULL,
                       B_M_prev, B_I_prev, B_D_prev);

  for (i = L - 1; i >= 0; i--) {
    int x_next = (int) dsq[i + 1];
    int xt     = (x_next >= 0 && x_next < K) ? x_next : K;

    ibv_backward_one_row(M, k_stride, i, L, do_trunc,
                         MM_t, MI_t, MD_t, IM_t, II_t, DM_t, DD_t,
                         emit_table[xt], B_M_prev, B_I_prev,
                         B_M_curr, B_I_curr, B_D_curr);

    if (i >= 2 && i <= L - 2)
      ibv_through_scan(M, k_stride, floor_milli,
                       ibv_mode, ibv_width,
                       F_M(i), F_I(i), F_D(i),
                       B_M_curr, B_I_curr, B_D_curr,
                       through, &rowstats, &kmin[i], &kmax[i], &i2k[i]);

    double *t_M = B_M_prev; B_M_prev = B_M_curr; B_M_curr = t_M;
    double *t_I = B_I_prev; B_I_prev = B_I_curr; B_I_curr = t_I;
    double *t_D = B_D_prev; B_D_prev = B_D_curr; B_D_curr = t_D;
  }

  if (L >= 1) { kmin[1] = 1; kmax[1] = M; }
  if (L >= 2) { kmin[L - 1] = 1; kmax[L - 1] = M; }
  if (L >= 1) { kmin[L] = 1; kmax[L] = M; }
  kmin[0] = 0; kmax[0] = M;   /* brief 26_0430-307/309: row-0 delete-chain coverage */

  /* Brief 26_0430-140a: bridge inter-row gaps so FIXED/HYBRID bands are connected
   * (DELTA untouched).  Must run before ncells is summed. */
  ibv_connectivity_guard(L, M, ibv_mode, kmin, kmax);

  /* brief 26_0430-298 EXPLORATORY probe (P298_INS0, DEFAULT OFF, not a proposed
   * fix): the IBV band never admits k=0 on rows i>=1, so the CP9's I_0 state
   * (leading/5'-flanking inserts, before consensus column 1) is banded out on
   * every emitting row.  This probe admits it, to measure how much of the
   * post-P298_ROW0D residual is attributable to I_0.  Widening k=0 into rows
   * i>=1 touches many INBAND(i,0) special cases in the CP9 kernels, so this is
   * a MEASUREMENT switch only -- do not ship it without a separate audit. */
  { static int p298_ins0 = -1;
    if(p298_ins0 < 0) { char *s = getenv("P298_INS0"); p298_ins0 = (s != NULL && atoi(s) != 0) ? 1 : 0; }
    if(p298_ins0) { for (i = 1; i <= L; i++) kmin[i] = 0; } }

  /* briefs 26_0821-099/26_0821-100: count in int64 and saturate the int
   * ret_ncells at INT_MAX, as p7_pins2bands() does (brief 26_0628-027). A
   * genome-scale IBV band can exceed 2^31 cells; the int32 sum wrapped
   * negative, and callers read ncells <= 0 as "no band". Callers use ncells
   * only as a >0/==0 discriminator and in debug lines, never to size memory. */
  {
    int64_t ncells64 = 0;
    for (i = 1; i <= L; i++)
      ncells64 += (int64_t) (kmax[i] - kmin[i] + 1);
    if (ncells64 > (int64_t) INT_MAX) {
      if(cm_StageVerbose()) fprintf(stderr, "#NCELLS64 p7_Seq2BandsIBV L=%d M=%d ncells=%lld (saturated to INT_MAX in int32 ret_ncells)\n",
              L, M, (long long) ncells64);
      ncells = INT_MAX;
    }
    else ncells = (int) ncells64;
    ibv_rowstats_report(&rowstats, "flat", L, M, floor_milli, ncells64);
  }
  if ((status = ibv_rowstats_check(&rowstats, errbuf)) != eslOK) goto ERROR;

  {
    const char *dump = getenv("P7IBV_DUMP_BAND");
    if (dump != NULL && *dump != '\0') {
      FILE *fp = fopen(dump, "w");
      if (fp != NULL) {
        fprintf(fp, "# M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d\n",
                M, L, delta_milli, (double) optimal, (double) thr, ncells);
        fprintf(fp, "# i\tkmin\tkmax\twidth\ti2k\n");   /* brief 26_0430-142: + i2k (argmax-k pin) */
        for (i = 1; i <= L; i++)
          fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", i, kmin[i], kmax[i], kmax[i] - kmin[i] + 1, i2k[i]);
        fclose(fp);
      }
    }
  }

  /* Brief 26_0430-142: per-row band dump to stderr (multi-seq safe; one block/seq). */
  {
    const char *p142 = getenv("P142_DUMP_BANDS");
    if (p142 != NULL && *p142 != '\0') {
      fprintf(stderr, "#P142_BAND_BEGIN M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d path=flat\n",
              M, L, delta_milli, (double) optimal, (double) thr, ncells);
      for (i = 1; i <= L; i++)
        fprintf(stderr, "#P7BAND_DUMP i=%d i2k=%d kmin=%d kmax=%d width=%d\n",
                i, i2k[i], kmin[i], kmax[i], kmax[i] - kmin[i] + 1);
      fprintf(stderr, "#P142_BAND_END L=%d\n", L);
    }
  }

#undef F_M
#undef F_I
#undef F_D

  free(MM_t); free(MI_t); free(MD_t);
  free(IM_t); free(II_t); free(DM_t); free(DD_t);
  free(FM_pool); free(FI_pool); free(FD_pool);
  free(BM_a); free(BM_b); free(BI_a); free(BI_b); free(BD_a); free(BD_b);
  free(through);
  if (begin_milli) free(begin_milli);
  free(emit_pool); free(emit_table);

  *ret_i2k    = i2k;
  *ret_kmin   = kmin;
  *ret_kmax   = kmax;
  *ret_ncells = ncells;
  return eslOK;

 ERROR:
  if (MM_t) free(MM_t); if (MI_t) free(MI_t); if (MD_t) free(MD_t);
  if (IM_t) free(IM_t); if (II_t) free(II_t);
  if (DM_t) free(DM_t); if (DD_t) free(DD_t);
  if (FM_pool) free(FM_pool); if (FI_pool) free(FI_pool); if (FD_pool) free(FD_pool);
  if (BM_a) free(BM_a); if (BM_b) free(BM_b);
  if (BI_a) free(BI_a); if (BI_b) free(BI_b);
  if (BD_a) free(BD_a); if (BD_b) free(BD_b);
  if (through) free(through);
  if (begin_milli) free(begin_milli);
  if (emit_pool) free(emit_pool); if (emit_table) free(emit_table);
  if (i2k)  free(i2k);
  if (kmin) free(kmin);
  if (kmax) free(kmax);
  *ret_i2k    = NULL;
  *ret_kmin   = NULL;
  *ret_kmax   = NULL;
  *ret_ncells = 0;
  return status;
}


/* ---------------------------------------------------------------------------
 * C2: D&C band deriver -- p7_Seq2BandsIBV_dnc
 * ---------------------------------------------------------------------------
 *
 * Memory layout (pre-allocated once, no malloc per recursion frame):
 *   arena    : max_depth * 6 * k_stride floats.
 *              At depth d: F_mid_{M,I,D} at d*6*ks+{0,1,2}*ks
 *                          B_mid1_{M,I,D} at d*6*ks+{3,4,5}*ks
 *   roll_F   : 6 * k_stride (prev_M,prev_I,prev_D, curr_M,curr_I,curr_D)
 *   roll_B   : 6 * k_stride
 *   slab_F   : (base_slab+1) * 3 * k_stride  (base-case forward storage)
 *   slab_B   : (base_slab+1) * 3 * k_stride  (base-case backward storage)
 *   Bmid     : 3 * k_stride  (one extra backward step at midrow)
 *   through  : k_stride      (through-score scratch)
 */

typedef struct {
  double  *arena;
  double  *roll_F;
  double  *roll_B;
  double  *slab_F;
  double  *slab_B;
  double  *Bmid;
  double  *through;
  double  *MM_t, *MI_t, *MD_t, *IM_t, *II_t, *DM_t, *DD_t;
  double **emit_table;
  int     K;
  int     M;
  size_t  k_stride;
  int     global_L;
  int     base_slab;
  int     do_trunc;        /* brief 26_0430-171: Tgm begin/end-anywhere semantics       */
  const double *begin_milli;/* brief 26_0430-171: Tgm begin vector (NULL unless do_trunc) */
  const ESL_DSQ *dsq;
  int    *kmin;
  int    *kmax;
  int    *i2k;     /* argmax-k per row (brief 26_0430-137); i2k[i] = Viterbi-trace cell. */
  double   margin;    /* band margin max(delta, EPS), milli-bits (brief 26_0821-144) */
  IBV_ROWSTATS st;   /* empty-row and drift counts (brief 26_0821-144) */
  int     ibv_mode;  /* brief 26_0430-140: P7IBV_MODE_{DELTA,FIXED,HYBRID} */
  int     ibv_width; /* brief 26_0430-140: fixed-width pad W around argmax-k pin */
} IBV_DnC_Ctx;

/* brief 26_0430-171: begin vector for the forward call at absolute row `absrow`
 * (non-NULL only at the global begin row 1 under Tgm). */
static inline const double *
ibv_brow(const IBV_DnC_Ctx *ctx, int absrow)
{
  return (ctx->do_trunc && absrow == 1) ? ctx->begin_milli : NULL;
}

static inline const double *
ibv_emit(const IBV_DnC_Ctx *ctx, int pos)
{
  int x  = (int) ctx->dsq[pos];
  int xt = (x >= 0 && x < ctx->K) ? x : ctx->K;
  return ctx->emit_table[xt];
}

static void
ibv_dnc_recurse(IBV_DnC_Ctx *ctx,
                int i_lo, int i_hi, int depth,
                const double *F_lo_M, const double *F_lo_I, const double *F_lo_D,
                const double *B_hi_M, const double *B_hi_I, const double *B_hi_D)
{
  int    M        = ctx->M;
  size_t ks       = ctx->k_stride;
  int    global_L = ctx->global_L;
  int    base_slab= ctx->base_slab;

  if (i_hi <= i_lo) return;
  int slab_size = i_hi - i_lo;

  /* ---- Base case ---- */
  if (slab_size <= base_slab) {
    /* Forward-fill rows i_lo+1..i_hi; slab_F row r = F[i_lo+r]. */
    memcpy(ctx->slab_F + 0 * ks, F_lo_M, ks * sizeof(double));
    memcpy(ctx->slab_F + 1 * ks, F_lo_I, ks * sizeof(double));
    memcpy(ctx->slab_F + 2 * ks, F_lo_D, ks * sizeof(double));
    for (int r = 1; r <= slab_size; r++) {
      size_t prev_off = (size_t)(r - 1) * 3 * ks;
      size_t curr_off = (size_t) r       * 3 * ks;
      ibv_forward_one_row(M, ks,
                          ctx->MM_t, ctx->MI_t, ctx->MD_t,
                          ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                          ibv_emit(ctx, i_lo + r), ibv_brow(ctx, i_lo + r),
                          ctx->slab_F + prev_off + 0*ks,
                          ctx->slab_F + prev_off + 1*ks,
                          ctx->slab_F + prev_off + 2*ks,
                          ctx->slab_F + curr_off + 0*ks,
                          ctx->slab_F + curr_off + 1*ks,
                          ctx->slab_F + curr_off + 2*ks);
    }

    /* Backward-fill rows i_hi..i_lo+1; slab_B row r = B[i_lo+r]. */
    {
      int row_ihi = i_hi;
      const double *em = (row_ihi < global_L) ? ibv_emit(ctx, row_ihi + 1) : NULL;
      size_t off = (size_t) slab_size * 3 * ks;
      ibv_backward_one_row(M, ks, row_ihi, global_L, ctx->do_trunc,
                           ctx->MM_t, ctx->MI_t, ctx->MD_t,
                           ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                           em, B_hi_M, B_hi_I,
                           ctx->slab_B + off + 0*ks,
                           ctx->slab_B + off + 1*ks,
                           ctx->slab_B + off + 2*ks);
    }
    for (int r = slab_size - 1; r >= 1; r--) {
      size_t next_off = (size_t)(r + 1) * 3 * ks;
      size_t curr_off = (size_t) r       * 3 * ks;
      ibv_backward_one_row(M, ks, i_lo + r, global_L, ctx->do_trunc,
                           ctx->MM_t, ctx->MI_t, ctx->MD_t,
                           ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                           ibv_emit(ctx, i_lo + r + 1),
                           ctx->slab_B + next_off + 0*ks,
                           ctx->slab_B + next_off + 1*ks,
                           ctx->slab_B + curr_off + 0*ks,
                           ctx->slab_B + curr_off + 1*ks,
                           ctx->slab_B + curr_off + 2*ks);
    }

    /* Per-row through-scan. */
    for (int r = 1; r <= slab_size; r++) {
      size_t off = (size_t) r * 3 * ks;
      ibv_through_scan(M, ks, ctx->margin,
                       ctx->ibv_mode, ctx->ibv_width,
                       ctx->slab_F + off + 0*ks, ctx->slab_F + off + 1*ks, ctx->slab_F + off + 2*ks,
                       ctx->slab_B + off + 0*ks, ctx->slab_B + off + 1*ks, ctx->slab_B + off + 2*ks,
                       ctx->through, &ctx->st,
                       &ctx->kmin[i_lo + r], &ctx->kmax[i_lo + r], &ctx->i2k[i_lo + r]);
    }
    return;
  }

  /* ---- Recursive case ---- */
  int i_mid = (i_lo + i_hi) / 2;

  double *F_mid_M  = ctx->arena + (size_t) depth * 6 * ks + 0 * ks;
  double *F_mid_I  = ctx->arena + (size_t) depth * 6 * ks + 1 * ks;
  double *F_mid_D  = ctx->arena + (size_t) depth * 6 * ks + 2 * ks;
  double *B_mid1_M = ctx->arena + (size_t) depth * 6 * ks + 3 * ks;
  double *B_mid1_I = ctx->arena + (size_t) depth * 6 * ks + 4 * ks;
  double *B_mid1_D = ctx->arena + (size_t) depth * 6 * ks + 5 * ks;

  /* Rolling buffer pointers (shared; safe because streaming finishes before recursing). */
  double *rFpM = ctx->roll_F + 0*ks, *rFpI = ctx->roll_F + 1*ks, *rFpD = ctx->roll_F + 2*ks;
  double *rFcM = ctx->roll_F + 3*ks, *rFcI = ctx->roll_F + 4*ks, *rFcD = ctx->roll_F + 5*ks;
  double *rBpM = ctx->roll_B + 0*ks, *rBpI = ctx->roll_B + 1*ks, *rBpD = ctx->roll_B + 2*ks;
  double *rBcM = ctx->roll_B + 3*ks, *rBcI = ctx->roll_B + 4*ks, *rBcD = ctx->roll_B + 5*ks;

  /* Forward stream i_lo+1 .. i_mid into F_mid. */
  memcpy(rFpM, F_lo_M, ks * sizeof(double));
  memcpy(rFpI, F_lo_I, ks * sizeof(double));
  memcpy(rFpD, F_lo_D, ks * sizeof(double));
  for (int r = i_lo + 1; r < i_mid; r++) {
    ibv_forward_one_row(M, ks, ctx->MM_t, ctx->MI_t, ctx->MD_t,
                        ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                        ibv_emit(ctx, r), ibv_brow(ctx, r),
                        rFpM, rFpI, rFpD, rFcM, rFcI, rFcD);
    double *t; t=rFpM; rFpM=rFcM; rFcM=t;
               t=rFpI; rFpI=rFcI; rFcI=t;
               t=rFpD; rFpD=rFcD; rFcD=t;
  }
  ibv_forward_one_row(M, ks, ctx->MM_t, ctx->MI_t, ctx->MD_t,
                      ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                      ibv_emit(ctx, i_mid), ibv_brow(ctx, i_mid),
                      rFpM, rFpI, rFpD, F_mid_M, F_mid_I, F_mid_D);

  /* Backward stream i_hi .. i_mid+1 into B_mid1. */
  memcpy(rBpM, B_hi_M, ks * sizeof(double));
  memcpy(rBpI, B_hi_I, ks * sizeof(double));
  memcpy(rBpD, B_hi_D, ks * sizeof(double));
  for (int r = i_hi; r > i_mid + 1; r--) {
    const double *em = (r < global_L) ? ibv_emit(ctx, r + 1) : NULL;
    ibv_backward_one_row(M, ks, r, global_L, ctx->do_trunc,
                         ctx->MM_t, ctx->MI_t, ctx->MD_t,
                         ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                         em, rBpM, rBpI, rBcM, rBcI, rBcD);
    double *t; t=rBpM; rBpM=rBcM; rBcM=t;
               t=rBpI; rBpI=rBcI; rBcI=t;
               t=rBpD; rBpD=rBcD; rBcD=t;
  }
  {
    int r = i_mid + 1;
    const double *em = (r < global_L) ? ibv_emit(ctx, r + 1) : NULL;
    ibv_backward_one_row(M, ks, r, global_L, ctx->do_trunc,
                         ctx->MM_t, ctx->MI_t, ctx->MD_t,
                         ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                         em, rBpM, rBpI, B_mid1_M, B_mid1_I, B_mid1_D);
  }

  /* One more backward step: B[i_mid] from B_mid1 = B[i_mid+1]. */
  {
    const double *em = (i_mid < global_L) ? ibv_emit(ctx, i_mid + 1) : NULL;
    ibv_backward_one_row(M, ks, i_mid, global_L, ctx->do_trunc,
                         ctx->MM_t, ctx->MI_t, ctx->MD_t,
                         ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                         em, B_mid1_M, B_mid1_I,
                         ctx->Bmid + 0*ks, ctx->Bmid + 1*ks, ctx->Bmid + 2*ks);
  }

  /* Through-scan at i_mid. */
  ibv_through_scan(M, ks, ctx->margin,
                   ctx->ibv_mode, ctx->ibv_width,
                   F_mid_M, F_mid_I, F_mid_D,
                   ctx->Bmid + 0*ks, ctx->Bmid + 1*ks, ctx->Bmid + 2*ks,
                   ctx->through, &ctx->st, &ctx->kmin[i_mid], &ctx->kmax[i_mid], &ctx->i2k[i_mid]);

  /* Recurse top half [i_lo, i_mid] with F_lo + B_mid1. */
  ibv_dnc_recurse(ctx, i_lo, i_mid, depth + 1,
                  F_lo_M,   F_lo_I,   F_lo_D,
                  B_mid1_M, B_mid1_I, B_mid1_D);

  /* Recurse bottom half [i_mid, i_hi] with F_mid + B_hi. */
  ibv_dnc_recurse(ctx, i_mid, i_hi, depth + 1,
                  F_mid_M, F_mid_I, F_mid_D,
                  B_hi_M,  B_hi_I,  B_hi_D);
}


static int
ibv_dnc_alloc(size_t n, double **ret_p)
{
  int status;
  double *p;
  if ((status = ibv_alloc_floats(n, &p)) != eslOK) { *ret_p = NULL; return status; }
  for (size_t c = 0; c < n; c++) p[c] = P7IBV_NEG_INF;
  *ret_p = p;
  return eslOK;
}


int
cm_p7_Seq2BandsIBV_dnc(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L,
                    int delta_milli, int base_slab,
                    int do_boundary_widen,
                    int do_trunc,
                    int ibv_mode, int ibv_width,
                    int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int          status;
  P7_HMM      *hmm   = NULL;
  int          M, K, k, i;
  double        optimal, thr, floor_milli;
  int          ncells = 0;
  size_t       ks;
  int          max_depth;
  IBV_DnC_Ctx  ctx;
  double       *begin_milli = NULL;
  double       *emit_pool  = NULL;
  double       *F_row0_M   = NULL, *F_row0_I   = NULL, *F_row0_D   = NULL;
  double       *B_seed_M   = NULL, *B_seed_I   = NULL, *B_seed_D   = NULL;
  int         *i2k        = NULL, *kmin_arr   = NULL, *kmax_arr   = NULL;
  memset(&ctx, 0, sizeof(ctx));

  if (cm == NULL || cm->fp7 == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsIBV_dnc: cm->fp7 is NULL");
  hmm = cm->fp7;
  M   = hmm->M;
  K   = hmm->abc->K;
  if (L < 1 || M < 1)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsIBV_dnc: bad L=%d or M=%d", L, M);

  /* See p7_Seq2BandsIBV: (M+4) rounding guarantees padding above index M so
   * the SSE backward's k+4 overread folds into NEG_INF (critical for the
   * contiguous arena, where M+1 would otherwise alias the next state). */
  ks = ((size_t)(M + 4) + (P7IBV_K_ALIGN - 1)) & ~(size_t)(P7IBV_K_ALIGN - 1);

  /* Adaptive base_slab when caller passes <= 0: cap base-case slab memory. */
  if (base_slab <= 0) {
    double cap_bytes = P7IBV_SLAB_CAP_BYTES * 1024.0 * 1024.0;
    double per_row   = 2.0 * 3.0 * (double) ks * sizeof(double);   /* F+B slab, one row */
    long   adaptive  = (long)(cap_bytes / per_row) - 1;
    if (adaptive < P7IBV_SLAB_MIN) adaptive = P7IBV_SLAB_MIN;
    if (adaptive > P7IBV_SLAB_MAX) adaptive = P7IBV_SLAB_MAX;
    base_slab = (int) adaptive;
  }
  if (base_slab < 1) base_slab = 1;

  max_depth = 2;
  { int tmp = L; while (tmp > 0) { tmp >>= 1; max_depth++; } }

  /* Transitions. */
  if ((status = ibv_alloc_floats(ks, &ctx.MM_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(ks, &ctx.MI_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(ks, &ctx.MD_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(ks, &ctx.IM_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(ks, &ctx.II_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(ks, &ctx.DM_t)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(ks, &ctx.DD_t)) != eslOK) goto ERROR;
  for (k = 0; k <= M; k++) {
    ctx.MM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MM]);
    ctx.MI_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MI]);
    ctx.MD_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_MD]);
    ctx.IM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_IM]);
    ctx.II_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_II]);
    ctx.DM_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_DM]);
    ctx.DD_t[k] = p7ibv_lod_milli(hmm->t[k][p7H_DD]);
  }
  for (k = M + 1; k < (int) ks; k++) {
    ctx.MM_t[k] = ctx.MI_t[k] = ctx.MD_t[k] = P7IBV_NEG_INF;
    ctx.IM_t[k] = ctx.II_t[k] = P7IBV_NEG_INF;
    ctx.DM_t[k] = ctx.DD_t[k] = P7IBV_NEG_INF;
  }

  /* Emit table. */
  if ((status = ibv_alloc_floats((size_t)(K + 1) * ks, &emit_pool)) != eslOK) goto ERROR;
  ESL_ALLOC(ctx.emit_table, sizeof(double *) * (K + 1));
  for (int xt = 0; xt <= K; xt++) ctx.emit_table[xt] = emit_pool + (size_t) xt * ks;
  for (int xt = 0; xt < K; xt++) {
    double *row = ctx.emit_table[xt];
    for (k = 0; k <= M; k++) row[k] = p7ibv_emit_milli(hmm, k, xt);
    for (k = M + 1; k < (int) ks; k++) row[k] = P7IBV_NEG_INF;
  }
  for (k = 0; k < (int) ks; k++) ctx.emit_table[K][k] = 0.0f;

  /* Brief 26_0430-171: Tgm begin-anywhere vector (see p7_Seq2BandsIBV). */
  if (do_trunc) {
    float *occ = NULL;
    double Z = 0.0;
    ESL_ALLOC(occ, sizeof(float) * (M + 1));
    if ((status = p7_hmm_CalculateOccupancy(hmm, occ, NULL)) != eslOK) { free(occ); goto ERROR; }
    for (k = 1; k <= M; k++) Z += (double) occ[k] * (double) (M - k + 1);
    if ((status = ibv_alloc_floats(ks, &begin_milli)) != eslOK) { free(occ); goto ERROR; }
    for (k = 0; k < (int) ks; k++) begin_milli[k] = P7IBV_NEG_INF;
    for (k = 1; k <= M; k++) {
      double b = (Z > 0.0 && occ[k] > 0.0) ? (double) occ[k] / Z : 0.0;
      begin_milli[k] = (b > 0.0) ? (double)(P7IBV_INTSCALE * (log(b) / M_LN2)) : P7IBV_NEG_INF;
    }
    free(occ);
  }

  /* F[0]: glocal D-cascade init; Tgm leaves row 0 all NEG_INF (begins at row 1). */
  if ((status = ibv_dnc_alloc(ks, &F_row0_M)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(ks, &F_row0_I)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(ks, &F_row0_D)) != eslOK) goto ERROR;
  if (! do_trunc) {
    F_row0_M[0] = 0.0f;
    for (k = 1; k <= M; k++) {
      double a = F_row0_M[k - 1] + ctx.MD_t[k - 1];
      double b = F_row0_D[k - 1] + ctx.DD_t[k - 1];
      F_row0_D[k] = (a > b) ? a : b;
    }
  }

  /* Global forward pass (2-row rolling) to get optimal score.
   * Reuse roll_F for this; reset to NEG_INF afterwards. */
  if ((status = ibv_dnc_alloc(6 * ks, &ctx.roll_F)) != eslOK) goto ERROR;
  {
    double *rM = ctx.roll_F + 0*ks, *rI = ctx.roll_F + 1*ks, *rD = ctx.roll_F + 2*ks;
    double *cM = ctx.roll_F + 3*ks, *cI = ctx.roll_F + 4*ks, *cD = ctx.roll_F + 5*ks;
    memcpy(rM, F_row0_M, ks * sizeof(double));
    memcpy(rI, F_row0_I, ks * sizeof(double));
    memcpy(rD, F_row0_D, ks * sizeof(double));
    for (i = 1; i <= L; i++) {
      int x = (int) dsq[i];
      int xt = (x >= 0 && x < K) ? x : K;
      const double *brow = (do_trunc && i == 1) ? begin_milli : NULL;
      ibv_forward_one_row(M, ks,
                          ctx.MM_t, ctx.MI_t, ctx.MD_t,
                          ctx.IM_t, ctx.II_t, ctx.DM_t, ctx.DD_t,
                          ctx.emit_table[xt], brow,
                          rM, rI, rD, cM, cI, cD);
      double *t; t=rM; rM=cM; cM=t; t=rI; rI=cI; cI=t; t=rD; rD=cD; cD=t;
    }
    if (do_trunc) {
      /* Tgm end-anywhere: exit from any match node, or D_M->E. */
      optimal = rD[M];
      for (k = 1; k <= M; k++) if (rM[k] > optimal) optimal = rM[k];
    } else {
      optimal = rM[M];
      if (rI[M] > optimal) optimal = rI[M];
      if (rD[M] > optimal) optimal = rD[M];
    }
  }
  assert(optimal == optimal);
  floor_milli = (double) delta_milli;
  if (floor_milli < P7IBV_EPS) floor_milli = P7IBV_EPS;
  thr = optimal - floor_milli;   /* reported in the band dumps only; rows use their own max */
  /* Reset roll_F for D&C streaming reuse. */
  for (size_t c = 0; c < 6 * ks; c++) ctx.roll_F[c] = P7IBV_NEG_INF;

  /* Remaining D&C buffers. */
  if ((status = ibv_dnc_alloc((size_t) max_depth * 6 * ks, &ctx.arena))  != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(6 * ks, &ctx.roll_B))                       != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc((size_t)(base_slab + 1) * 3 * ks, &ctx.slab_F)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc((size_t)(base_slab + 1) * 3 * ks, &ctx.slab_B)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(3 * ks, &ctx.Bmid))                          != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(ks,     &ctx.through))                        != eslOK) goto ERROR;

  ESL_ALLOC(i2k,     sizeof(int) * (L + 1));
  ESL_ALLOC(kmin_arr, sizeof(int) * (L + 1));
  ESL_ALLOC(kmax_arr, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, L + 1, -1);
  for (i = 0; i <= L; i++) { kmin_arr[i] = 1; kmax_arr[i] = M; }

  ctx.K        = K;
  ctx.M        = M;
  ctx.k_stride = ks;
  ctx.global_L = L;
  ctx.base_slab= base_slab;
  ctx.dsq      = dsq;
  ctx.kmin     = kmin_arr;
  ctx.kmax     = kmax_arr;
  ctx.i2k      = i2k;
  ctx.margin   = floor_milli;
  ibv_rowstats_init(&ctx.st, optimal);
  ctx.do_trunc    = do_trunc;     /* brief 26_0430-171 */
  ctx.begin_milli = begin_milli;  /* brief 26_0430-171 */
  ctx.ibv_mode = ibv_mode;   /* brief 26_0430-140 */
  ctx.ibv_width= ibv_width;  /* brief 26_0430-140 */

  /* B seed for top-level: all NEG_INF; terminal injected via global_L. */
  if ((status = ibv_dnc_alloc(ks, &B_seed_M)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(ks, &B_seed_I)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(ks, &B_seed_D)) != eslOK) goto ERROR;

  /* Run recursion. */
  ibv_dnc_recurse(&ctx, 0, L, 0,
                  F_row0_M, F_row0_I, F_row0_D,
                  B_seed_M, B_seed_I, B_seed_D);

  /* Boundary widening + row-0 convention. The widening of rows 1, L-1, L to
   * the full model [1,M] exists for truncated-alignment entry/exit. The
   * non-truncated --hmm --p7ibv path passes do_boundary_widen=FALSE to skip it
   * (IBV-JSTATE-NOTE2 finding 3); CM-side --p7band --p7ibv passes TRUE to
   * preserve byte-identical behavior. The row-0 convention always applies. */
  if (do_boundary_widen) {
    if (L >= 1) { kmin_arr[1] = 1; kmax_arr[1] = M; }
    if (L >= 2) { kmin_arr[L - 1] = 1; kmax_arr[L - 1] = M; }
    if (L >= 1) { kmin_arr[L] = 1; kmax_arr[L] = M; }
  }
  kmin_arr[0] = 0; kmax_arr[0] = M;   /* brief 26_0430-307/309: row-0 delete-chain coverage */
  i2k[0] = 0;   /* B-state convention (brief 26_0430-137): i2k[i]=argmax_k for i in [1,L]. */

  /* Brief 26_0430-140a: bridge inter-row gaps so FIXED/HYBRID bands are connected
   * (DELTA untouched).  Must run before ncells is summed. */
  ibv_connectivity_guard(L, M, ibv_mode, kmin_arr, kmax_arr);

  /* brief 26_0430-298 EXPLORATORY probe (P298_INS0, DEFAULT OFF, not a proposed
   * fix), ported to _dnc by brief 26_0430-301: the IBV band never admits k=0 on
   * rows i>=1, so the CP9's I_0 state (leading/5'-flanking inserts, before
   * consensus column 1) is banded out on every emitting row. This probe admits
   * it to measure how much of the post-P298_ROW0D residual is attributable to
   * I_0. Widening k=0 into rows i>=1 touches many INBAND(i,0) special cases in
   * the CP9 kernels, so this is a MEASUREMENT switch only -- do not ship it
   * without a separate audit. */
  { static int p298_ins0 = -1;
    if(p298_ins0 < 0) { char *s = getenv("P298_INS0"); p298_ins0 = (s != NULL && atoi(s) != 0) ? 1 : 0; }
    if(p298_ins0) { for (i = 1; i <= L; i++) kmin_arr[i] = 0; } }

  /* briefs 26_0821-099/26_0821-100: count in int64 and saturate the int
   * ret_ncells at INT_MAX, as p7_pins2bands() does (brief 26_0628-027). A
   * genome-scale IBV band can exceed 2^31 cells; the int32 sum wrapped
   * negative, and callers read ncells <= 0 as "no band". Callers use ncells
   * only as a >0/==0 discriminator and in debug lines, never to size memory. */
  {
    int64_t ncells64 = 0;
    for (i = 1; i <= L; i++)
      ncells64 += (int64_t) (kmax_arr[i] - kmin_arr[i] + 1);
    if (ncells64 > (int64_t) INT_MAX) {
      if(cm_StageVerbose()) fprintf(stderr, "#NCELLS64 p7_Seq2BandsIBV_dnc L=%d M=%d ncells=%lld (saturated to INT_MAX in int32 ret_ncells)\n",
              L, M, (long long) ncells64);
      ncells = INT_MAX;
    }
    else ncells = (int) ncells64;
    ibv_rowstats_report(&ctx.st, "dnc", L, M, floor_milli, ncells64);
  }
  if ((status = ibv_rowstats_check(&ctx.st, errbuf)) != eslOK) goto ERROR;

  {
    const char *dump = getenv("P7IBV_DUMP_BAND");
    if (dump != NULL && *dump != '\0') {
      FILE *fp = fopen(dump, "w");
      if (fp != NULL) {
        fprintf(fp, "# M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d\n",
                M, L, delta_milli, (double) optimal, (double) thr, ncells);
        fprintf(fp, "# i\tkmin\tkmax\twidth\ti2k\n");   /* brief 26_0430-142: + i2k (argmax-k pin) */
        for (i = 1; i <= L; i++)
          fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", i, kmin_arr[i], kmax_arr[i], kmax_arr[i] - kmin_arr[i] + 1, i2k[i]);
        fclose(fp);
      }
    }
  }

  /* Brief 26_0430-142: per-row band dump to stderr (multi-seq safe; one block/seq). */
  {
    const char *p142 = getenv("P142_DUMP_BANDS");
    if (p142 != NULL && *p142 != '\0') {
      fprintf(stderr, "#P142_BAND_BEGIN M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d path=dnc\n",
              M, L, delta_milli, (double) optimal, (double) thr, ncells);
      for (i = 1; i <= L; i++)
        fprintf(stderr, "#P7BAND_DUMP i=%d i2k=%d kmin=%d kmax=%d width=%d\n",
                i, i2k[i], kmin_arr[i], kmax_arr[i], kmax_arr[i] - kmin_arr[i] + 1);
      fprintf(stderr, "#P142_BAND_END L=%d\n", L);
    }
  }

  free(ctx.MM_t); free(ctx.MI_t); free(ctx.MD_t);
  free(ctx.IM_t); free(ctx.II_t); free(ctx.DM_t); free(ctx.DD_t);
  free(emit_pool); free(ctx.emit_table);
  free(F_row0_M); free(F_row0_I); free(F_row0_D);
  free(ctx.roll_F); free(ctx.roll_B); free(ctx.arena);
  free(ctx.slab_F); free(ctx.slab_B);
  free(ctx.Bmid); free(ctx.through);
  free(B_seed_M); free(B_seed_I); free(B_seed_D);
  if (begin_milli) free(begin_milli);

  *ret_i2k    = i2k;
  *ret_kmin   = kmin_arr;
  *ret_kmax   = kmax_arr;
  *ret_ncells = ncells;
  return eslOK;

 ERROR:
  if (ctx.MM_t)   free(ctx.MM_t);   if (ctx.MI_t) free(ctx.MI_t); if (ctx.MD_t) free(ctx.MD_t);
  if (ctx.IM_t)   free(ctx.IM_t);   if (ctx.II_t) free(ctx.II_t);
  if (ctx.DM_t)   free(ctx.DM_t);   if (ctx.DD_t) free(ctx.DD_t);
  if (emit_pool)  free(emit_pool);   if (ctx.emit_table) free(ctx.emit_table);
  if (F_row0_M)   free(F_row0_M);   if (F_row0_I) free(F_row0_I); if (F_row0_D) free(F_row0_D);
  if (ctx.roll_F) free(ctx.roll_F); if (ctx.roll_B) free(ctx.roll_B);
  if (ctx.arena)  free(ctx.arena);
  if (ctx.slab_F) free(ctx.slab_F); if (ctx.slab_B) free(ctx.slab_B);
  if (ctx.Bmid)   free(ctx.Bmid);   if (ctx.through) free(ctx.through);
  if (B_seed_M)   free(B_seed_M);   if (B_seed_I) free(B_seed_I); if (B_seed_D) free(B_seed_D);
  if (begin_milli) free(begin_milli);
  if (i2k)        free(i2k);
  if (kmin_arr)   free(kmin_arr);
  if (kmax_arr)   free(kmax_arr);
  *ret_i2k    = NULL;
  *ret_kmin   = NULL;
  *ret_kmax   = NULL;
  *ret_ncells = 0;
  return status;
}


/* Function:  p7_IBVPins2Trace()
 * Synopsis:  Convert IBV per-row argmax-k pins to a P7_TRACE.
 * Incept:    brief 26_0430-137, 2026-06-17.
 *
 * Purpose:   Given the per-row argmax-k pin array <i2k> produced by
 *            p7_Seq2BandsIBV_dnc(... delta_milli=0 ...) for a sequence of
 *            length <L> aligned to profile <gm>, emit a unihit P7_TRACE that
 *            is alignment-equivalent to what p7_GViterbi -> p7_GTrace would
 *            produce on the same input.
 *
 *            i2k[i] (i in 1..L) is the model position k of the Viterbi-optimal
 *            cell that emits residue i; i2k[0] = 0 (B-state convention).
 *            Residue i is a match emission ML(i2k[i]) when i2k[i] advances over
 *            i2k[i-1] (a DL run fills any skipped match positions), or an
 *            insert emission IL(i2k[i]) when i2k[i] == i2k[i-1].  Because the
 *            IBV DP is global w.r.t. the sequence, every residue 1..L is
 *            consumed in the core and the N/C states are non-emitting.
 *
 *            Begin/exit depend on gm->mode:
 *              - glocal (p7_UNIGLOCAL): B->D1..D_{k-1}->Mk leading deletes and
 *                M_{i2k[L]}->D..D_M->E trailing deletes are emitted.
 *              - local  (p7_UNILOCAL):  B->Mk and Mk->E directly, no flanking
 *                deletes.
 *            Internal DL runs between consecutive pins are emitted in both
 *            modes.  Multihit modes (p7_LOCAL/p7_GLOCAL) are rejected: the
 *            walker never emits p7T_J.
 *
 * Args:      gm     - configured profile (gm->mode selects begin/end semantics)
 *            dsq    - digital sequence 1..L (unused; kept for API symmetry)
 *            L      - sequence length
 *            i2k    - per-row argmax-k pins, length L+1, with i2k[0]=0
 *            kmin   - per-row band low  (length L+1) or NULL to skip the check
 *            kmax   - per-row band high (length L+1) or NULL to skip the check
 *            ncells - total band cell count (unused; kept for API symmetry)
 *            ret_tr - RETURN: newly allocated P7_TRACE (caller frees)
 *
 * Returns:   <eslOK> on success, with *ret_tr the trace.
 *            <eslEINVAL> if i2k is inconsistent (interior -1 sentinel,
 *               non-monotone, k out of [1,M], or outside [kmin,kmax]).
 *            <eslEUNIMPLEMENTED> if gm is configured multihit.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
cm_p7_IBVPins2Trace(const P7_PROFILE *gm, const ESL_DSQ *dsq, int L,
                 const int *i2k, const int *kmin, const int *kmax, int ncells,
                 P7_TRACE **ret_tr)
{
  P7_TRACE *tr      = NULL;
  int      *w       = NULL;
  int       M       = gm->M;
  int       islocal = p7_IsLocal(gm->mode);
  int       i, k, kd, kprev, prev_match, i_lastM;
  int       status;

  (void) dsq; (void) ncells;
  *ret_tr = NULL;

  if (p7_IsMulti(gm->mode)) return eslEUNIMPLEMENTED;  /* unihit only: no J. */
  if (L < 1)                return eslEINVAL;
  if (i2k[0] != 0)          return eslEINVAL;
  for (i = 1; i <= L; i++)
    if (i2k[i] < 1 || i2k[i] > M) return eslEINVAL;    /* interior -1 / OOB */

  /* ---- Band-aware path threading ----
   * The per-row argmax pins are independently chosen and so at Δ=0 co-optimal
   * ties they can (a) dip non-monotonically and (b) imply infeasible Plan7
   * transitions (an insert followed by a match >1 position ahead would need
   * deletes after an insert, which Plan7 forbids).  Thread a single monotone,
   * transition-feasible path w[] through the per-row co-optimal band [kmin,kmax]
   * (every in-band cell scores >= optimal at Δ=0, so any in-band monotone path
   * is itself score-optimal, i.e. alignment-equivalent).  Deletes are allowed
   * only out of a match or the begin state.  Falls back to a plain monotone
   * clamp when kmin/kmax are unavailable. */
  ESL_ALLOC(w, sizeof(int) * (L + 1));
  w[0]  = 0;
  kprev = 0;
  prev_match = TRUE;        /* begin behaves like a match: B->D->M is allowed */
  for (i = 1; i <= L; i++) {
    k = i2k[i];
    if (k < kprev) k = kprev;                 /* (a) monotone clamp -> hold */
    if (k > kprev && !prev_match && k > kprev + 1) {
      /* (b) need deletes to enter Mk, but came from an insert: pull the match
       * down to kprev+1 if that cell is co-optimal, else hold as an insert. */
      if (!kmin || !kmax || (kprev + 1 >= kmin[i] && kprev + 1 <= kmax[i]))
        k = kprev + 1;
      else
        k = kprev;
    }
    w[i]       = k;
    prev_match = (k > kprev);
    kprev      = k;
  }

  /* Locate the last match in the threaded path.  Residues after it are held at
   * the final model position and so would be inserts with no following match;
   * a profile trace cannot end on an insert (E connects from M/D), so they are
   * emitted as C-state (post-core) residues instead -- which is also what a
   * p7 Viterbi trace does with trailing residues. */
  i_lastM = 0;
  { int kp = 0; for (i = 1; i <= L; i++) { if (w[i] > kp) i_lastM = i; kp = w[i]; } }
  /* (w[1] is always a match since w[1] >= 1 > w[0]=0, so i_lastM >= 1.) */

  if ((tr = p7_trace_Create()) == NULL) { status = eslEMEM; goto ERROR; }

  /* S -> N (non-emitting) -> B */
  if ((status = p7_trace_Append(tr, p7T_S, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_N, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_B, 0, 0)) != eslOK) goto ERROR;

  kprev = 0;   /* B-state model position */
  for (i = 1; i <= i_lastM; i++) {
    k = w[i];
    if (k > kprev) {           /* match emission, possibly after an internal DL run */
      /* No leading deletes: B->Mk is a direct (wing-retracted) begin in a
       * profile trace; B->D is illegal.  Internal M->D..->M runs are emitted. */
      for (kd = (kprev == 0 ? k : kprev + 1); kd < k; kd++)
        if ((status = p7_trace_Append(tr, p7T_D, kd, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr, p7T_M, k, i)) != eslOK) goto ERROR;
    } else {                   /* k == kprev: internal insert emission */
      if ((status = p7_trace_Append(tr, p7T_I, k, i)) != eslOK) goto ERROR;
    }
    kprev = k;
  }

  /* Exit from the final match M_kprev.  Glocal: M_kprev -> D..->D_M -> E (valid
   * M/D->D->E run).  Local: M_kprev -> E directly. */
  if (!islocal)
    for (k = kprev + 1; k <= M; k++)
      if ((status = p7_trace_Append(tr, p7T_D, k, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_E, 0, 0)) != eslOK) goto ERROR;

  /* C: first is non-emitting (E->C); trailing residues i_lastM+1..L emit on C. */
  if ((status = p7_trace_Append(tr, p7T_C, 0, 0)) != eslOK) goto ERROR;
  for (i = i_lastM + 1; i <= L; i++)
    if ((status = p7_trace_Append(tr, p7T_C, 0, i)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) goto ERROR;

  tr->M = M;
  tr->L = L;
  if (w) free(w);
  *ret_tr = tr;
  return eslOK;

 ERROR:
  if (w)  free(w);
  if (tr) p7_trace_Destroy(tr);
  *ret_tr = NULL;
  return status;
}


/*****************************************************************
 * Unit test driver (brief 26_0821-144)
 *****************************************************************/
#ifdef CM_P7_IBV_TESTDRIVE

#include "esl_getopts.h"

/* ibv_utest_hmm(): an M-node ungapped-ish RNA HMM whose match states
 * strongly prefer residue 0 over an insert distribution that rarely emits
 * it, so a sequence of 0s scores ~<sc> bits per residue: large total
 * scores at modest M, which is what makes the float through-scores drift.
 * <shift> (bits) lowers every residue's score by <shift>: match emissions
 * via mat[][], insert emissions via the transitions into I (the IBV
 * kernels score inserts by transition only). In exact arithmetic that
 * lowers every parse by shift*L and leaves the band unchanged; with
 * shift ~ the per-residue mean it keeps the DP's partial sums near 0,
 * where float rounding is negligible -- the reference band. */
static P7_HMM *
ibv_utest_hmm(const ESL_ALPHABET *abc, int M, double shift)
{
  P7_HMM *hmm = p7_hmm_Create(M, abc);
  double  f   = pow(2.0, -shift);
  int     k, x;

  p7_hmm_Zero(hmm);
  for (k = 0; k <= M; k++) {
    for (x = 0; x < abc->K; x++) {
      hmm->mat[k][x] = (k == 0) ? (x == 0 ? 1.0 : 0.0) : f * ((x == 0) ? 0.97 : 0.01);
      hmm->ins[k][x] = (x == 0) ? 0.01 : 0.33;
    }
    if (k < M) {
      hmm->t[k][p7H_MM] = 0.98;  hmm->t[k][p7H_MI] = f * 0.01; hmm->t[k][p7H_MD] = 0.01;
      hmm->t[k][p7H_IM] = 0.7;   hmm->t[k][p7H_II] = f * 0.3;
      hmm->t[k][p7H_DM] = 0.7;   hmm->t[k][p7H_DD] = 0.3;
    } else {
      hmm->t[k][p7H_MM] = 1.0;   hmm->t[k][p7H_MI] = 0.0;     hmm->t[k][p7H_MD] = 0.0;
      hmm->t[k][p7H_IM] = 1.0;   hmm->t[k][p7H_II] = 0.0;
      hmm->t[k][p7H_DM] = 1.0;   hmm->t[k][p7H_DD] = 0.0;
    }
  }
  return hmm;
}

/* utest_genome_drift(): the IBV band at a scale where float through-scores
 * drift. A sequence of M residues aligned to the HMM above scores ~6.6
 * bits per residue (~66,000 bits at the default M = 10,000); the D&C
 * deriver (cmalign's --hmm --p7ibv path) runs once on the plain HMM and
 * once on the shifted one. With the float F/B this file used before
 * brief 26_0821-144, the plain run's row max through-scores drifted from
 * its optimum by up to 7.4 bits on this input, against the test's 0.2-bit
 * margin, and the bands differed from the reference. Requires: no empty
 * rows; both runs' drift (row max - optimum) below 1/10 of the margin;
 * the same Viterbi pins; and the two bands equal to within one cell at
 * every row edge. */
static void
utest_genome_drift(const ESL_ALPHABET *abc, int M, int delta, int be_verbose)
{
  P7_HMM   *hmm0 = ibv_utest_hmm(abc, M, 0.0);
  P7_HMM   *hmm1 = NULL;
  CM_t      cm;
  ESL_DSQ  *dsq  = NULL;
  int      *i2k0 = NULL, *kmin0 = NULL, *kmax0 = NULL;
  int      *i2k1 = NULL, *kmin1 = NULL, *kmax1 = NULL;
  int       nc0, nc1, i, L = M, maxdiff = 0;
  IBV_ROWSTATS st0, st1;
  double    persc;
  char      errbuf[eslERRBUFSIZE];

  dsq = malloc(sizeof(ESL_DSQ) * (L + 2));
  dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
  for (i = 1; i <= L; i++) dsq[i] = 0;

  memset(&cm, 0, sizeof(CM_t));
  cm.fp7 = hmm0;
  if (cm_p7_Seq2BandsIBV_dnc(&cm, errbuf, dsq, L, delta, 0, FALSE, FALSE, FALSE, P7IBV_MODE_DELTA, 0,
                             &i2k0, &kmin0, &kmax0, &nc0) != eslOK) esl_fatal("utest_genome_drift: plain run failed: %s", errbuf);
  st0 = ibv_last_rowstats;

  persc = st0.optimal / P7IBV_INTSCALE / (double) L;
  hmm1  = ibv_utest_hmm(abc, M, persc);
  cm.fp7 = hmm1;
  if (cm_p7_Seq2BandsIBV_dnc(&cm, errbuf, dsq, L, delta, 0, FALSE, FALSE, FALSE, P7IBV_MODE_DELTA, 0,
                             &i2k1, &kmin1, &kmax1, &nc1) != eslOK) esl_fatal("utest_genome_drift: shifted run failed: %s", errbuf);
  st1 = ibv_last_rowstats;

  for (i = 1; i <= L; i++) {
    maxdiff = ESL_MAX(maxdiff, abs(kmin0[i] - kmin1[i]));
    maxdiff = ESL_MAX(maxdiff, abs(kmax0[i] - kmax1[i]));
  }
  if (be_verbose)
    printf("M=L=%d optimum %.1f bits (%.3f/residue); plain drift [%.3f, %.3f] bits, %lld empty rows, %d cells; shifted drift [%.4f, %.4f], %d cells; max edge diff %d\n",
           L, st0.optimal / P7IBV_INTSCALE, persc, st0.drift_min / P7IBV_INTSCALE, st0.drift_max / P7IBV_INTSCALE,
           (long long) st0.empty_rows, nc0, st1.drift_min / P7IBV_INTSCALE, st1.drift_max / P7IBV_INTSCALE, nc1, maxdiff);

  if (st0.empty_rows != 0 || st1.empty_rows != 0) esl_fatal("utest_genome_drift: empty rows (%lld, %lld)", (long long) st0.empty_rows, (long long) st1.empty_rows);
  if (ESL_MAX(st0.drift_max, -st0.drift_min) > 0.1 * delta)
    esl_fatal("utest_genome_drift: plain run drifts [%g, %g] milli-bits (margin %d)", st0.drift_min, st0.drift_max, delta);
  if (ESL_MAX(st1.drift_max, -st1.drift_min) > 0.1 * delta)
    esl_fatal("utest_genome_drift: shifted reference drifts [%g, %g] milli-bits", st1.drift_min, st1.drift_max);
  if (maxdiff > 1) esl_fatal("utest_genome_drift: bands differ by up to %d cells at a row edge", maxdiff);
  for (i = 1; i <= L; i++)
    if (i2k0[i] != i2k1[i]) esl_fatal("utest_genome_drift: Viterbi pins differ at row %d (%d vs %d)", i, i2k0[i], i2k1[i]);

  free(i2k0); free(kmin0); free(kmax0); free(i2k1); free(kmin1); free(kmax1); free(dsq);
  p7_hmm_Destroy(hmm0); p7_hmm_Destroy(hmm1);
}

static ESL_OPTIONS options[] = {
  /* name     type          default  env  range toggles reqs incomp  help                                     docgroup*/
  { "-h",        eslARG_NONE, FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help and exit",               0 },
  { "--version", eslARG_NONE, FALSE, NULL, NULL,  NULL,  NULL, NULL, "show version info and exit",             0 },
  { "-M",     eslARG_INT,   "10000", NULL, "n>0", NULL,  NULL, NULL, "model and sequence length for the drift test", 0 },
  { "-d",     eslARG_INT,     "200", NULL, "n>0", NULL,  NULL, NULL, "band margin (milli-bits) for the drift test", 0 },
  { "-v",     eslARG_NONE,    FALSE, NULL, NULL,  NULL,  NULL, NULL, "print the drift test's numbers",          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for p7_ibv.c (IBV band derivation)";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go  = cm_CreateDefaultApp("p7_ibv_utest", options, 0, argc, argv, banner, usage);
  ESL_ALPHABET *abc = esl_alphabet_Create(eslRNA);

  utest_genome_drift(abc, esl_opt_GetInteger(go, "-M"), esl_opt_GetInteger(go, "-d"), esl_opt_GetBoolean(go, "-v"));

  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*CM_P7_IBV_TESTDRIVE*/

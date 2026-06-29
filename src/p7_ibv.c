/* p7_ibv.c -- F+B direct-band derivation: SSE per-row primitives + D&C wrapper
 *
 * Brief 124 C1: refactor brief 121's monolithic p7_Seq2BandsIBV into three
 *   reusable per-row SSE primitives (ibv_forward_one_row, ibv_backward_one_row,
 *   ibv_through_scan).  No algorithm change; byte-exact vs brief 121 C3.
 *
 * Brief 124 C2: p7_Seq2BandsIBV_dnc wraps the C1 primitives in a recursive
 *   divide-and-conquer band deriver (O(M * log L) peak memory).  Accessed via
 *   --p7ibv --p7ibv-mem; --p7ibv alone still calls the flat p7_Seq2BandsIBV.
 *
 * SSE/memory conventions (unchanged from brief 121 C2/C3):
 *   k_stride = ((M+4+15) & ~15)  (16-float align + >=3 pad slots above M)
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
#define P7IBV_OPTIMAL_SANITY  1.0e7f

/* Float-precision guard on the optimal score.  The band threshold test
 * (through >= optimal - margin, margin = max(delta, EPS)) stays reliable
 * while the float ULP at the score magnitude is << margin.  ULP(through)
 * ~= 2*|optimal|*2^-23, so the test is safe while
 *   |optimal| <= margin * 2^23 / (2 * SAFETY).
 * With SAFETY=8 the bound is margin * 2^19.  At the default delta=3000 that
 * is ~1.6e9 milli-bits, which admits genome-scale models (HSV optimal is
 * O(1e7-1e8)).  The old fixed 1e7 guard (tuned to ULP < EPS=1) was far too
 * conservative -- it tripped on HSV/MPXV even though the real band margin is
 * delta, not EPS.  Floor the guard at OPTIMAL_SANITY so tiny deltas still
 * permit reasonable scores. */
#define P7IBV_ULP_SAFETY_SHIFT 19   /* margin << 2^19 ULPs */

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
ibv_alloc_floats(size_t n, float **ret_p)
{
  void *p = NULL;
  if (posix_memalign(&p, 64, sizeof(float) * n) != 0 || p == NULL) {
    *ret_p = NULL;
    return eslEMEM;
  }
  *ret_p = (float *) p;
  return eslOK;
}

static inline float
p7ibv_lod_milli(float p)
{
  if (p <= 0.0f) return P7IBV_NEG_INF;
  return (float)(P7IBV_INTSCALE * (log((double) p) / M_LN2));
}

static inline float
p7ibv_emit_milli(const P7_HMM *hmm, int k, int x)
{
  if (x < 0 || x >= hmm->abc->K) return 0.0f;
  float pm = hmm->mat[k][x];
  float pi = hmm->ins[k][x];
  if (pm <= 0.0f || pi <= 0.0f) return P7IBV_NEG_INF;
  return (float)(P7IBV_INTSCALE * (log((double) pm / (double) pi) / M_LN2));
}

static inline __m128
p7ibv_mm_max3(__m128 a, __m128 b, __m128 c)
{
  return _mm_max_ps(_mm_max_ps(a, b), c);
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
 * begin_milli (brief 171, Tgm begin-anywhere): if non-NULL, fold a local
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
                    const float *MM_t, const float *MI_t, const float *MD_t,
                    const float *IM_t, const float *II_t,
                    const float *DM_t, const float *DD_t,
                    const float *emit_row,
                    const float *begin_milli,
                    const float *FM_prev, const float *FI_prev, const float *FD_prev,
                    float *FM_curr, float *FI_curr, float *FD_curr)
{
  int k;
  int kpref_end = (M < 3) ? M : 3;

  for (k = 0; k <= kpref_end; k++) {
    float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF, cD = P7IBV_NEG_INF;
    {
      float a = FM_prev[k] + MI_t[k];
      float b = FI_prev[k] + II_t[k];
      cI = (a > b) ? a : b;
    }
    if (k >= 1) {
      float a = FM_prev[k - 1] + MM_t[k - 1];
      float b = FI_prev[k - 1] + IM_t[k - 1];
      float c = FD_prev[k - 1] + DM_t[k - 1];
      float m = (a > b) ? a : b;
      if (c > m) m = c;
      cM = m + emit_row[k];
      if (begin_milli != NULL) {
        float bc = begin_milli[k] + emit_row[k];
        if (bc > cM) cM = bc;
      }
      float a2 = FM_curr[k - 1] + MD_t[k - 1];
      float b2 = FD_curr[k - 1] + DD_t[k - 1];
      cD = (a2 > b2) ? a2 : b2;
    }
    FM_curr[k] = cM;
    FI_curr[k] = cI;
    FD_curr[k] = cD;
  }

  int k_sse_start = 4;
  int k_sse_end   = k_sse_start;
  while (k_sse_end + 3 <= M) k_sse_end += 4;
  for (k = k_sse_start; k < k_sse_end; k += 4) {
    __m128 m_prev = _mm_loadu_ps(&FM_prev[k - 1]);
    __m128 i_prev = _mm_loadu_ps(&FI_prev[k - 1]);
    __m128 d_prev = _mm_loadu_ps(&FD_prev[k - 1]);
    __m128 t_mm   = _mm_loadu_ps(&MM_t[k - 1]);
    __m128 t_im   = _mm_loadu_ps(&IM_t[k - 1]);
    __m128 t_dm   = _mm_loadu_ps(&DM_t[k - 1]);
    __m128 a      = _mm_add_ps(m_prev, t_mm);
    __m128 b      = _mm_add_ps(i_prev, t_im);
    __m128 c      = _mm_add_ps(d_prev, t_dm);
    __m128 mx     = p7ibv_mm_max3(a, b, c);
    if (begin_milli != NULL) mx = _mm_max_ps(mx, _mm_loadu_ps(&begin_milli[k]));
    __m128 e_vec  = _mm_loadu_ps(&emit_row[k]);
    _mm_storeu_ps(&FM_curr[k], _mm_add_ps(mx, e_vec));

    __m128 fm_k = _mm_loadu_ps(&FM_prev[k]);
    __m128 fi_k = _mm_loadu_ps(&FI_prev[k]);
    __m128 t_mi = _mm_loadu_ps(&MI_t[k]);
    __m128 t_ii = _mm_loadu_ps(&II_t[k]);
    __m128 a2   = _mm_add_ps(fm_k, t_mi);
    __m128 b2   = _mm_add_ps(fi_k, t_ii);
    _mm_storeu_ps(&FI_curr[k], _mm_max_ps(a2, b2));
  }

  for (k = k_sse_end; k <= M; k++) {
    float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
    {
      float a = FM_prev[k] + MI_t[k];
      float b = FI_prev[k] + II_t[k];
      cI = (a > b) ? a : b;
    }
    if (k >= 1) {
      float a = FM_prev[k - 1] + MM_t[k - 1];
      float b = FI_prev[k - 1] + IM_t[k - 1];
      float c = FD_prev[k - 1] + DM_t[k - 1];
      float m = (a > b) ? a : b;
      if (c > m) m = c;
      cM = m + emit_row[k];
      if (begin_milli != NULL) {
        float bc = begin_milli[k] + emit_row[k];
        if (bc > cM) cM = bc;
      }
    }
    FM_curr[k] = cM;
    FI_curr[k] = cI;
  }

  int kd_start = (kpref_end + 1 > M) ? (M + 1) : (kpref_end + 1);
  for (k = kd_start; k <= M; k++) {
    float a = FM_curr[k - 1] + MD_t[k - 1];
    float b = FD_curr[k - 1] + DD_t[k - 1];
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
 * do_trunc (brief 171, Tgm end-anywhere): when set, the terminal injection at
 * row L mirrors local exit -- every match state M_k may exit to E with score 0
 * (esc=0 in HMMER local mode), instead of the glocal forced exit from node M
 * (delete-cascade to D_M).  The delete cascade to D_M->E is kept (D_M->E is
 * allowed even in local mode); inserts still cannot be the last emitted state.
 * Only affects the terminal branch (i == global_L); the interior backward
 * recursion is mode-independent because the end is carried in by the terminal.
 */
static void
ibv_backward_one_row(int M, size_t k_stride, int i, int global_L, int do_trunc,
                     const float *MM_t, const float *MI_t, const float *MD_t,
                     const float *IM_t, const float *II_t,
                     const float *DM_t, const float *DD_t,
                     const float *emit_row_next,
                     const float *BM_next,
                     const float *BI_next,
                     float *BM_curr, float *BI_curr, float *BD_curr)
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
      float bv = BD_curr[k + 1];
      BM_curr[k] = MD_t[k] + bv;
      BD_curr[k] = DD_t[k] + bv;
      BI_curr[k] = P7IBV_NEG_INF;
    }
    return;
  }

  BD_curr[M] = P7IBV_NEG_INF;
  for (k = M - 1; k >= 0; k--) {
    float bv_m = BM_next[k + 1] + emit_row_next[k + 1];
    float bv_d = BD_curr[k + 1];
    float a = DM_t[k] + bv_m;
    float b = DD_t[k] + bv_d;
    BD_curr[k] = (a > b) ? a : b;
  }

  int k_sse_end = 0;
  while (k_sse_end + 3 <= M) k_sse_end += 4;
  for (k = 0; k < k_sse_end; k += 4) {
    __m128 bm_n   = _mm_loadu_ps(&BM_next[k + 1]);
    __m128 e_n    = _mm_loadu_ps(&emit_row_next[k + 1]);
    __m128 bv_M   = _mm_add_ps(bm_n, e_n);
    __m128 bv_I   = _mm_loadu_ps(&BI_next[k]);
    __m128 bv_D   = _mm_loadu_ps(&BD_curr[k + 1]);

    __m128 t_mm   = _mm_loadu_ps(&MM_t[k]);
    __m128 t_mi   = _mm_loadu_ps(&MI_t[k]);
    __m128 t_md   = _mm_loadu_ps(&MD_t[k]);
    __m128 cM_out = p7ibv_mm_max3(_mm_add_ps(t_mm, bv_M),
                                  _mm_add_ps(t_mi, bv_I),
                                  _mm_add_ps(t_md, bv_D));
    _mm_storeu_ps(&BM_curr[k], cM_out);

    __m128 t_im   = _mm_loadu_ps(&IM_t[k]);
    __m128 t_ii   = _mm_loadu_ps(&II_t[k]);
    __m128 cI_out = _mm_max_ps(_mm_add_ps(t_im, bv_M),
                               _mm_add_ps(t_ii, bv_I));
    _mm_storeu_ps(&BI_curr[k], cI_out);
  }

  for (k = k_sse_end; k <= M; k++) {
    float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
    if (k + 1 <= M) {
      float bv_m = BM_next[k + 1] + emit_row_next[k + 1];
      float a = MM_t[k] + bv_m;
      float b = IM_t[k] + bv_m;
      if (a > cM) cM = a;
      if (b > cI) cI = b;
    }
    {
      float bv_i = BI_next[k];
      float a = MI_t[k] + bv_i;
      float b = II_t[k] + bv_i;
      if (a > cM) cM = a;
      if (b > cI) cI = b;
    }
    if (k + 1 <= M) {
      float bv_d = BD_curr[k + 1];
      float a = MD_t[k] + bv_d;
      if (a > cM) cM = a;
    }
    BM_curr[k] = cM;
    BI_curr[k] = cI;
  }
}


/* ibv_through_scan -- per-row kmin/kmax from F+B through-score.
 *
 * through_scratch: caller-owned k_stride buffer.  Returns (1,M) on empty row.
 */
static void
ibv_through_scan(int M, size_t k_stride, float thr,
                 int ibv_mode, int ibv_width,
                 const float *FM, const float *FI, const float *FD,
                 const float *BM, const float *BI, const float *BD,
                 float *through_scratch,
                 int *ret_kmin, int *ret_kmax, int *ret_kargmax)
{
  int k;
  int k_thru_end = 0;
  while (k_thru_end + 3 <= M) k_thru_end += 4;
  (void) k_stride;

  for (k = 0; k < k_thru_end; k += 4) {
    __m128 a = _mm_add_ps(_mm_loadu_ps(&FM[k]), _mm_loadu_ps(&BM[k]));
    __m128 b = _mm_add_ps(_mm_loadu_ps(&FI[k]), _mm_loadu_ps(&BI[k]));
    __m128 c = _mm_add_ps(_mm_loadu_ps(&FD[k]), _mm_loadu_ps(&BD[k]));
    _mm_storeu_ps(&through_scratch[k], p7ibv_mm_max3(a, b, c));
  }
  for (k = k_thru_end; k <= M; k++) {
    float t_m = FM[k] + BM[k];
    float t_i = FI[k] + BI[k];
    float t_d = FD[k] + BD[k];
    float t = t_m;
    if (t_i > t) t = t_i;
    if (t_d > t) t = t_d;
    through_scratch[k] = t;
  }

  int   row_kmin = -1, row_kmax = -1;
  int   k_argmax = -1;
  float t_argmax = P7IBV_NEG_INF;
  for (k = 1; k <= M; k++) {
    float t = through_scratch[k];
    if (t < P7IBV_HALF_NEG_INF) continue;
    /* Band (kmin/kmax) uses the full through-score incl. delete cells.
     *
     * The argmax-k pin (brief 137), however, must be over EMITTING states (M,I)
     * only: residue i is emitted by a match or insert, never a delete.  A
     * delete cell D(i,k) sits on the optimal path too (with through == optimal,
     * at a higher k than the emitter) and float F+B reconstruction can tip the
     * tie so the delete edges out the emitter -- which would make i2k[i] point
     * at a model position where residue i is actually *deleted*.  Restricting
     * the argmax to FM+BM / FI+BI excludes delete cells from the pin. */
    float t_emit = FM[k] + BM[k];
    { float ti = FI[k] + BI[k]; if (ti > t_emit) t_emit = ti; }
    if (t_emit >= P7IBV_HALF_NEG_INF && t_emit > t_argmax) { t_argmax = t_emit; k_argmax = k; }
    if (t >= thr) {
      if (row_kmin < 0) row_kmin = k;
      row_kmax = k;
    }
  }
  if (row_kmin < 0) { *ret_kmin = 1; *ret_kmax = M; }
  else              { *ret_kmin = row_kmin; *ret_kmax = row_kmax; }
  if (ret_kargmax) *ret_kargmax = k_argmax;

  /* Brief 140: enrich the per-row band using the argmax-k pin (k_argmax = i2k[i],
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
 * Brief 140a.  The narrow per-row FIXED/HYBRID bands (centered on the argmax-k
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
 * (The brief 140a sketch used kmin[i-1]-1 for the second test; that is too weak
 *  by one -- kmax[i]==kmin[i-1]-1 still leaves no in-band predecessor.)
 *
 * The guard widens band[i] minimally to satisfy this for every i.  DELTA mode is
 * skipped (its posterior-reachable bands are connected by construction).  Rows 0
 * and 1 connect to the B-state via the begin transition (no diagonal needed), so
 * the forward sweep starts at i=2.  A single forward sweep is provably sufficient
 * to guarantee a complete path (hence a finite total); the backward sweep is kept
 * for robustness -- it can only widen further and so cannot re-introduce a gap.
 * Widening keeps kmin>=1 and kmax<=M; the final clamp is purely defensive.
 */
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
 * p7_Seq2BandsIBV -- C1 rewrite (byte-exact vs brief 121 C3)
 * ---------------------------------------------------------------------------*/

int
p7_Seq2BandsIBV(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L, int delta_milli,
                int do_trunc,
                int ibv_mode, int ibv_width,
                int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int       status;
  P7_HMM   *hmm = NULL;
  int       M;
  int       i, k;
  int       K;
  float    *begin_milli = NULL;
  float    *MM_t = NULL, *MI_t = NULL, *MD_t = NULL;
  float    *IM_t = NULL, *II_t = NULL;
  float    *DM_t = NULL, *DD_t = NULL;
  float    *FM_pool = NULL, *FI_pool = NULL, *FD_pool = NULL;
  float    *BM_a = NULL, *BM_b = NULL;
  float    *BI_a = NULL, *BI_b = NULL;
  float    *BD_a = NULL, *BD_b = NULL;
  float    *through = NULL;
  float    *emit_pool = NULL;
  float   **emit_table = NULL;
  int      *i2k = NULL, *kmin = NULL, *kmax = NULL;
  float     optimal, thr;
  float     floor_milli;
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
   * _mm_loadu_ps(&BM_next[k+1]) reach index k+4 = k_sse_end <= M+1, so we
   * need k_stride >= M+2. Rounding (M+1) up to 16 gives ZERO padding when
   * M+1 is a multiple of 16 (M == 15 mod 16, e.g. M=287), leaving index M+1
   * out of bounds. Round (M+4) up instead to guarantee >=3 padding slots
   * (all NEG_INF), so the overread folds harmlessly. (Latent in brief 121's
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
  ESL_ALLOC(emit_table, sizeof(float *) * (K + 1));
  for (int xt = 0; xt <= K; xt++) emit_table[xt] = emit_pool + (size_t) xt * k_stride;
  for (int xt = 0; xt < K; xt++) {
    float *row = emit_table[xt];
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

  /* Brief 171: Tgm begin-anywhere vector.  begin_milli[k] = milli-bit log2 of
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
      begin_milli[k] = (b > 0.0) ? (float)(P7IBV_INTSCALE * (log(b) / M_LN2)) : P7IBV_NEG_INF;
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
      float *fm0 = F_M(0);
      float *fd0 = F_D(0);
      for (k = 1; k <= M; k++) {
        float a = fm0[k - 1] + MD_t[k - 1];
        float b = fd0[k - 1] + DD_t[k - 1];
        fd0[k] = (a > b) ? a : b;
      }
    }
  }

  for (i = 1; i <= L; i++) {
    int x  = (int) dsq[i];
    int xt = (x >= 0 && x < K) ? x : K;
    const float *brow = (do_trunc && i == 1) ? begin_milli : NULL;
    ibv_forward_one_row(M, k_stride,
                        MM_t, MI_t, MD_t, IM_t, II_t, DM_t, DD_t,
                        emit_table[xt], brow,
                        F_M(i-1), F_I(i-1), F_D(i-1),
                        F_M(i),   F_I(i),   F_D(i));
  }

  {
    float *fm_L = F_M(L);
    float *fi_L = F_I(L);
    float *fd_L = F_D(L);
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
  floor_milli = (float) delta_milli;
  if (floor_milli < P7IBV_EPS) floor_milli = P7IBV_EPS;
  {
    float guard = floor_milli * (float)(1u << P7IBV_ULP_SAFETY_SHIFT);
    if (guard < P7IBV_OPTIMAL_SANITY) guard = P7IBV_OPTIMAL_SANITY;
    if (fabsf(optimal) > guard)
      ESL_FAIL(eslEINVAL, errbuf,
               "p7_Seq2BandsIBV: |optimal|=%g milli-bits exceeds guard %g "
               "(float ULP approaches delta=%d margin; need double precision)",
               (double) optimal, (double) guard, delta_milli);
  }
  thr = optimal - floor_milli;

  ESL_ALLOC(i2k,  sizeof(int) * (L + 1));
  ESL_ALLOC(kmin, sizeof(int) * (L + 1));
  ESL_ALLOC(kmax, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, L + 1, -1);
  for (i = 0; i <= L; i++) { kmin[i] = 1; kmax[i] = M; }

  float *B_M_prev = BM_a, *B_I_prev = BI_a, *B_D_prev = BD_a;
  float *B_M_curr = BM_b, *B_I_curr = BI_b, *B_D_curr = BD_b;

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
      ibv_through_scan(M, k_stride, thr,
                       ibv_mode, ibv_width,
                       F_M(i), F_I(i), F_D(i),
                       B_M_curr, B_I_curr, B_D_curr,
                       through, &kmin[i], &kmax[i], &i2k[i]);

    float *t_M = B_M_prev; B_M_prev = B_M_curr; B_M_curr = t_M;
    float *t_I = B_I_prev; B_I_prev = B_I_curr; B_I_curr = t_I;
    float *t_D = B_D_prev; B_D_prev = B_D_curr; B_D_curr = t_D;
  }

  if (L >= 1) { kmin[1] = 1; kmax[1] = M; }
  if (L >= 2) { kmin[L - 1] = 1; kmax[L - 1] = M; }
  if (L >= 1) { kmin[L] = 1; kmax[L] = M; }
  kmin[0] = 0; kmax[0] = 0;

  /* Brief 140a: bridge inter-row gaps so FIXED/HYBRID bands are connected
   * (DELTA untouched).  Must run before ncells is summed. */
  ibv_connectivity_guard(L, M, ibv_mode, kmin, kmax);

  for (i = 1; i <= L; i++)
    ncells += (kmax[i] - kmin[i] + 1);

  {
    const char *dump = getenv("P7IBV_DUMP_BAND");
    if (dump != NULL && *dump != '\0') {
      FILE *fp = fopen(dump, "w");
      if (fp != NULL) {
        fprintf(fp, "# M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d\n",
                M, L, delta_milli, (double) optimal, (double) thr, ncells);
        fprintf(fp, "# i\tkmin\tkmax\twidth\ti2k\n");   /* brief 142: + i2k (argmax-k pin) */
        for (i = 1; i <= L; i++)
          fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", i, kmin[i], kmax[i], kmax[i] - kmin[i] + 1, i2k[i]);
        fclose(fp);
      }
    }
  }

  /* Brief 142: per-row band dump to stderr (multi-seq safe; one block/seq). */
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
  float  *arena;
  float  *roll_F;
  float  *roll_B;
  float  *slab_F;
  float  *slab_B;
  float  *Bmid;
  float  *through;
  float  *MM_t, *MI_t, *MD_t, *IM_t, *II_t, *DM_t, *DD_t;
  float **emit_table;
  int     K;
  int     M;
  size_t  k_stride;
  int     global_L;
  int     base_slab;
  int     do_trunc;        /* brief 171: Tgm begin/end-anywhere semantics       */
  const float *begin_milli;/* brief 171: Tgm begin vector (NULL unless do_trunc) */
  const ESL_DSQ *dsq;
  int    *kmin;
  int    *kmax;
  int    *i2k;     /* argmax-k per row (brief 137); i2k[i] = Viterbi-trace cell. */
  float   thr;
  int     ibv_mode;  /* brief 140: P7IBV_MODE_{DELTA,FIXED,HYBRID} */
  int     ibv_width; /* brief 140: fixed-width pad W around argmax-k pin */
  int     kband_pad; /* brief 172: k-band child-narrowing pad (do_kband path) */
  int     wide_thresh; /* brief 173: route nodes with band width >= this to the
                        * SSE full-M primitives; narrower nodes use scalar _b.   */
} IBV_DnC_Ctx;

/* brief 171: begin vector for the forward call at absolute row `absrow`
 * (non-NULL only at the global begin row 1 under Tgm). */
static inline const float *
ibv_brow(const IBV_DnC_Ctx *ctx, int absrow)
{
  return (ctx->do_trunc && absrow == 1) ? ctx->begin_milli : NULL;
}

static inline const float *
ibv_emit(const IBV_DnC_Ctx *ctx, int pos)
{
  int x  = (int) ctx->dsq[pos];
  int xt = (x >= 0 && x < ctx->K) ? x : ctx->K;
  return ctx->emit_table[xt];
}

static void
ibv_dnc_recurse(IBV_DnC_Ctx *ctx,
                int i_lo, int i_hi, int depth,
                const float *F_lo_M, const float *F_lo_I, const float *F_lo_D,
                const float *B_hi_M, const float *B_hi_I, const float *B_hi_D)
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
    memcpy(ctx->slab_F + 0 * ks, F_lo_M, ks * sizeof(float));
    memcpy(ctx->slab_F + 1 * ks, F_lo_I, ks * sizeof(float));
    memcpy(ctx->slab_F + 2 * ks, F_lo_D, ks * sizeof(float));
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
      const float *em = (row_ihi < global_L) ? ibv_emit(ctx, row_ihi + 1) : NULL;
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
      ibv_through_scan(M, ks, ctx->thr,
                       ctx->ibv_mode, ctx->ibv_width,
                       ctx->slab_F + off + 0*ks, ctx->slab_F + off + 1*ks, ctx->slab_F + off + 2*ks,
                       ctx->slab_B + off + 0*ks, ctx->slab_B + off + 1*ks, ctx->slab_B + off + 2*ks,
                       ctx->through,
                       &ctx->kmin[i_lo + r], &ctx->kmax[i_lo + r], &ctx->i2k[i_lo + r]);
    }
    return;
  }

  /* ---- Recursive case ---- */
  int i_mid = (i_lo + i_hi) / 2;

  float *F_mid_M  = ctx->arena + (size_t) depth * 6 * ks + 0 * ks;
  float *F_mid_I  = ctx->arena + (size_t) depth * 6 * ks + 1 * ks;
  float *F_mid_D  = ctx->arena + (size_t) depth * 6 * ks + 2 * ks;
  float *B_mid1_M = ctx->arena + (size_t) depth * 6 * ks + 3 * ks;
  float *B_mid1_I = ctx->arena + (size_t) depth * 6 * ks + 4 * ks;
  float *B_mid1_D = ctx->arena + (size_t) depth * 6 * ks + 5 * ks;

  /* Rolling buffer pointers (shared; safe because streaming finishes before recursing). */
  float *rFpM = ctx->roll_F + 0*ks, *rFpI = ctx->roll_F + 1*ks, *rFpD = ctx->roll_F + 2*ks;
  float *rFcM = ctx->roll_F + 3*ks, *rFcI = ctx->roll_F + 4*ks, *rFcD = ctx->roll_F + 5*ks;
  float *rBpM = ctx->roll_B + 0*ks, *rBpI = ctx->roll_B + 1*ks, *rBpD = ctx->roll_B + 2*ks;
  float *rBcM = ctx->roll_B + 3*ks, *rBcI = ctx->roll_B + 4*ks, *rBcD = ctx->roll_B + 5*ks;

  /* Forward stream i_lo+1 .. i_mid into F_mid. */
  memcpy(rFpM, F_lo_M, ks * sizeof(float));
  memcpy(rFpI, F_lo_I, ks * sizeof(float));
  memcpy(rFpD, F_lo_D, ks * sizeof(float));
  for (int r = i_lo + 1; r < i_mid; r++) {
    ibv_forward_one_row(M, ks, ctx->MM_t, ctx->MI_t, ctx->MD_t,
                        ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                        ibv_emit(ctx, r), ibv_brow(ctx, r),
                        rFpM, rFpI, rFpD, rFcM, rFcI, rFcD);
    float *t; t=rFpM; rFpM=rFcM; rFcM=t;
               t=rFpI; rFpI=rFcI; rFcI=t;
               t=rFpD; rFpD=rFcD; rFcD=t;
  }
  ibv_forward_one_row(M, ks, ctx->MM_t, ctx->MI_t, ctx->MD_t,
                      ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                      ibv_emit(ctx, i_mid), ibv_brow(ctx, i_mid),
                      rFpM, rFpI, rFpD, F_mid_M, F_mid_I, F_mid_D);

  /* Backward stream i_hi .. i_mid+1 into B_mid1. */
  memcpy(rBpM, B_hi_M, ks * sizeof(float));
  memcpy(rBpI, B_hi_I, ks * sizeof(float));
  memcpy(rBpD, B_hi_D, ks * sizeof(float));
  for (int r = i_hi; r > i_mid + 1; r--) {
    const float *em = (r < global_L) ? ibv_emit(ctx, r + 1) : NULL;
    ibv_backward_one_row(M, ks, r, global_L, ctx->do_trunc,
                         ctx->MM_t, ctx->MI_t, ctx->MD_t,
                         ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                         em, rBpM, rBpI, rBcM, rBcI, rBcD);
    float *t; t=rBpM; rBpM=rBcM; rBcM=t;
               t=rBpI; rBpI=rBcI; rBcI=t;
               t=rBpD; rBpD=rBcD; rBcD=t;
  }
  {
    int r = i_mid + 1;
    const float *em = (r < global_L) ? ibv_emit(ctx, r + 1) : NULL;
    ibv_backward_one_row(M, ks, r, global_L, ctx->do_trunc,
                         ctx->MM_t, ctx->MI_t, ctx->MD_t,
                         ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                         em, rBpM, rBpI, B_mid1_M, B_mid1_I, B_mid1_D);
  }

  /* One more backward step: B[i_mid] from B_mid1 = B[i_mid+1]. */
  {
    const float *em = (i_mid < global_L) ? ibv_emit(ctx, i_mid + 1) : NULL;
    ibv_backward_one_row(M, ks, i_mid, global_L, ctx->do_trunc,
                         ctx->MM_t, ctx->MI_t, ctx->MD_t,
                         ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                         em, B_mid1_M, B_mid1_I,
                         ctx->Bmid + 0*ks, ctx->Bmid + 1*ks, ctx->Bmid + 2*ks);
  }

  /* Through-scan at i_mid. */
  ibv_through_scan(M, ks, ctx->thr,
                   ctx->ibv_mode, ctx->ibv_width,
                   F_mid_M, F_mid_I, F_mid_D,
                   ctx->Bmid + 0*ks, ctx->Bmid + 1*ks, ctx->Bmid + 2*ks,
                   ctx->through, &ctx->kmin[i_mid], &ctx->kmax[i_mid], &ctx->i2k[i_mid]);

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
ibv_dnc_alloc(size_t n, float **ret_p)
{
  int status;
  float *p;
  if ((status = ibv_alloc_floats(n, &p)) != eslOK) { *ret_p = NULL; return status; }
  for (size_t c = 0; c < n; c++) p[c] = P7IBV_NEG_INF;
  *ret_p = p;
  return eslOK;
}


/* ---------------------------------------------------------------------------
 * Brief 172: k-banded scalar primitives + k-banded D&C  (the windowed-Viterbi
 *            i2k speed lever)
 * ---------------------------------------------------------------------------
 *
 * The optimal (Viterbi) path's model column k is monotone non-decreasing in a
 * left-right profile HMM (brief 168, verified: 0 backward steps).  So between
 * two EXACT through-scan pins at rows i_lo and i_hi (k_lo = i2k[i_lo] <= k_hi =
 * i2k[i_hi]) the path's k stays in [k_lo, k_hi].  Banding each D&C node's
 * forward/backward streams to its boundary-pin k-range is therefore EXACT for
 * the per-row argmax-k pin i2k: the path cell is interior to the band, and all
 * of a path cell's in-band predecessors are themselves interior (a path cell
 * at the band floor k_lo can only be reached by an I-transition that keeps k =
 * k_lo, never by an M/D from k_lo-1, which would mean k decreased below the
 * segment floor).  Out-of-band predecessors are read as NEG_INF via the
 * (k-1>=k_lo)/(k+1<=k_hi) guards, so the primitives are self-contained and do
 * not depend on neighbor-cell hygiene.
 *
 * The top level is the full model [0,M] (begin/end anywhere); children narrow
 * to the monotone tube around the just-computed exact midline pin.  Deep nodes
 * are thus cheap, and the per-level full-M factor of the unbanded D&C (~log L
 * full sweeps) collapses to ~1-2 sweeps of work.  These primitives are scalar
 * (not SSE): correctness-first, and the geometric majority of nodes are deep
 * narrow-band ones.  KPAD widens each child band by a few nodes as cheap
 * insurance against float ties at band edges (the exactness proof needs none).
 * In truncated (Tgm) mode the per-row argmax oracle wobbles between near-optimal
 * registers over a wider k-spread, so KPAD=64 (not 8) is needed for exact i2k
 * match there: empirically i2kdiff drops 92->0 on a truncated genome fragment
 * (M=152k) going 8->64, at negligible speed cost (the full-M top levels dominate).
 *
 * Only i2k is consumed by the windowed-Viterbi band (p7_Seq2BandsWV rebuilds
 * the band from i2k +/- nodepad), so the banded delta-cloud kmin/kmax (which a
 * narrow band would shrink relative to the full-M cloud) is irrelevant here.
 */
#define P7IBV_KPAD 64

static void
ibv_fwd_row_b(int M, int k_lo, int k_hi,
              const float *MM_t, const float *MI_t, const float *MD_t,
              const float *IM_t, const float *II_t,
              const float *DM_t, const float *DD_t,
              const float *emit_row, const float *begin_milli,
              const float *FM_prev, const float *FI_prev, const float *FD_prev,
              float *FM_curr, float *FI_curr, float *FD_curr)
{
  int k;
  if (k_lo < 0) k_lo = 0;
  if (k_hi > M) k_hi = M;
  for (k = k_lo; k <= k_hi; k++) {
    float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
    {
      float a = FM_prev[k] + MI_t[k];
      float b = FI_prev[k] + II_t[k];
      cI = (a > b) ? a : b;
    }
    if (k >= 1 && (k - 1) >= k_lo) {
      float a = FM_prev[k - 1] + MM_t[k - 1];
      float b = FI_prev[k - 1] + IM_t[k - 1];
      float c = FD_prev[k - 1] + DM_t[k - 1];
      float m = (a > b) ? a : b; if (c > m) m = c;
      cM = m + emit_row[k];
    }
    if (begin_milli != NULL && k >= 1) {
      float bc = begin_milli[k] + emit_row[k];
      if (bc > cM) cM = bc;
    }
    FM_curr[k] = cM;
    FI_curr[k] = cI;
  }
  for (k = k_lo; k <= k_hi; k++) {
    if (k >= 1 && (k - 1) >= k_lo) {
      float a = FM_curr[k - 1] + MD_t[k - 1];
      float b = FD_curr[k - 1] + DD_t[k - 1];
      FD_curr[k] = (a > b) ? a : b;
    } else {
      FD_curr[k] = P7IBV_NEG_INF;
    }
  }
}

static void
ibv_bwd_row_b(int M, int k_lo, int k_hi, int i, int global_L, int do_trunc,
              const float *MM_t, const float *MI_t, const float *MD_t,
              const float *IM_t, const float *II_t,
              const float *DM_t, const float *DD_t,
              const float *emit_row_next,
              const float *BM_next, const float *BI_next,
              float *BM_curr, float *BI_curr, float *BD_curr)
{
  int k;
  if (k_lo < 0) k_lo = 0;
  if (k_hi > M) k_hi = M;

  if (i == global_L) {
    if (do_trunc) {
      for (k = k_lo; k <= k_hi; k++) { BM_curr[k] = 0.0f; BI_curr[k] = P7IBV_NEG_INF; }
      for (k = k_hi; k >= k_lo; k--) {
        if      (k == M)            BD_curr[k] = 0.0f;
        else if ((k + 1) <= k_hi)   BD_curr[k] = DD_t[k] + BD_curr[k + 1];
        else                        BD_curr[k] = P7IBV_NEG_INF;
      }
      return;
    }
    for (k = k_hi; k >= k_lo; k--) {
      if (k == M) { BM_curr[k] = 0.0f; BI_curr[k] = 0.0f; BD_curr[k] = 0.0f; }
      else {
        float bv = ((k + 1) <= k_hi) ? BD_curr[k + 1] : P7IBV_NEG_INF;
        BM_curr[k] = MD_t[k] + bv;
        BD_curr[k] = DD_t[k] + bv;
        BI_curr[k] = P7IBV_NEG_INF;
      }
    }
    return;
  }

  for (k = k_hi; k >= k_lo; k--) {
    float bv_m = ((k + 1) <= k_hi) ? (BM_next[k + 1] + emit_row_next[k + 1]) : P7IBV_NEG_INF;
    float bv_d = ((k + 1) <= k_hi) ? BD_curr[k + 1] : P7IBV_NEG_INF;
    float a = DM_t[k] + bv_m;
    float b = DD_t[k] + bv_d;
    BD_curr[k] = (a > b) ? a : b;
  }
  for (k = k_lo; k <= k_hi; k++) {
    float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
    if ((k + 1) <= k_hi) {
      float bv_m = BM_next[k + 1] + emit_row_next[k + 1];
      float a = MM_t[k] + bv_m; if (a > cM) cM = a;
      float b = IM_t[k] + bv_m; if (b > cI) cI = b;
    }
    {
      float bv_i = BI_next[k];
      float a = MI_t[k] + bv_i; if (a > cM) cM = a;
      float b = II_t[k] + bv_i; if (b > cI) cI = b;
    }
    if ((k + 1) <= k_hi) {
      float bv_d = BD_curr[k + 1];
      float a = MD_t[k] + bv_d; if (a > cM) cM = a;
    }
    BM_curr[k] = cM;
    BI_curr[k] = cI;
  }
}

static void
ibv_through_b(int M, int k_lo, int k_hi, float thr,
              const float *FM, const float *FI, const float *FD,
              const float *BM, const float *BI, const float *BD,
              int *ret_kmin, int *ret_kmax, int *ret_kargmax)
{
  int   k, row_kmin = -1, row_kmax = -1, k_argmax = -1;
  float t_argmax = P7IBV_NEG_INF;
  if (k_lo < 1) k_lo = 1;
  if (k_hi > M) k_hi = M;
  for (k = k_lo; k <= k_hi; k++) {
    float t_m = FM[k] + BM[k];
    float t_i = FI[k] + BI[k];
    float t_d = FD[k] + BD[k];
    float t = t_m; if (t_i > t) t = t_i; if (t_d > t) t = t_d;
    if (t < P7IBV_HALF_NEG_INF) continue;
    float t_emit = (t_m > t_i) ? t_m : t_i;
    if (t_emit >= P7IBV_HALF_NEG_INF && t_emit > t_argmax) { t_argmax = t_emit; k_argmax = k; }
    if (t >= thr) { if (row_kmin < 0) row_kmin = k; row_kmax = k; }
  }
  if (row_kmin < 0) { *ret_kmin = 1; *ret_kmax = M; }
  else              { *ret_kmin = row_kmin; *ret_kmax = row_kmax; }
  if (ret_kargmax) *ret_kargmax = k_argmax;
}

/* ---------------------------------------------------------------------------
 * Brief 173 Part B: SSE-route the wide k-banded D&C levels.
 * ---------------------------------------------------------------------------
 *
 * The scalar _b primitives compute over [k_lo,k_hi] one cell at a time; the
 * full-M SSE primitives (ibv_forward_one_row / ibv_backward_one_row /
 * ibv_through_scan) compute over [0,M] in __m128 chunks (~4x/op faster).  Per
 * the 172 re-profile the WIDE top D&C levels (band ~= [1,M]) dominate the cost,
 * so we route them through the SSE full-M primitives and keep the scalar _b
 * primitives only for the narrow deep nodes.
 *
 * Exactness: the full-M SSE primitive computes a SUPERSET of the band
 * ([0,M] >= [k_lo,k_hi]).  For i2k -- the per-row EMITTING argmax -- the
 * monotone-k property (brief 172) guarantees the optimal path's cell at every
 * row is interior to [k_lo,k_hi], and full-M Viterbi computes the true optimal
 * value at that cell, so argmax over [1,M] == argmax over [k_lo,k_hi] == the
 * path cell.  (A non-path cell's full-M through-score is <= the global optimum,
 * so it cannot outscore the in-band path cell.)  The through-scan range
 * difference ([1,M] vs [k_lo,k_hi]) is therefore i2k-invariant.
 *
 * Safety invariant: band width is non-increasing down the recursion (each child
 * band is a subset of its parent's), so {nodes routed to SSE} = {nodes with
 * width >= wide_thresh} form a connected TOP-PREFIX of the tree.  Every
 * SSE-routed node thus receives full-M-valid boundary rows (F_lo/B_hi) from an
 * equally-or-wider SSE ancestor (or the root seeds), so the SSE forward/backward
 * never read stale out-of-band cells.  Narrow (scalar) descendants only read
 * their own band (guards read NEG_INF beyond k_lo-1/k_hi+1), so they are
 * unaffected by the extra full-M values a wide parent leaves in the buffers.
 *
 * We reuse the EXISTING SSE primitives unchanged (do NOT hand-roll a banded SSE
 * primitive), so we inherit the brief-124/125 k_stride overread fix for
 * M == 15 (mod 16) for free.
 */
static inline void
ibv_fwd_row_dispatch(IBV_DnC_Ctx *ctx, int wide, int k_lo, int k_hi, int absrow,
                     const float *FM_prev, const float *FI_prev, const float *FD_prev,
                     float *FM_curr, float *FI_curr, float *FD_curr)
{
  if (wide)
    ibv_forward_one_row(ctx->M, ctx->k_stride,
                        ctx->MM_t, ctx->MI_t, ctx->MD_t, ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                        ibv_emit(ctx, absrow), ibv_brow(ctx, absrow),
                        FM_prev, FI_prev, FD_prev, FM_curr, FI_curr, FD_curr);
  else
    ibv_fwd_row_b(ctx->M, k_lo, k_hi,
                  ctx->MM_t, ctx->MI_t, ctx->MD_t, ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                  ibv_emit(ctx, absrow), ibv_brow(ctx, absrow),
                  FM_prev, FI_prev, FD_prev, FM_curr, FI_curr, FD_curr);
}

static inline void
ibv_bwd_row_dispatch(IBV_DnC_Ctx *ctx, int wide, int k_lo, int k_hi, int i,
                     const float *emit_row_next,
                     const float *BM_next, const float *BI_next,
                     float *BM_curr, float *BI_curr, float *BD_curr)
{
  if (wide)
    ibv_backward_one_row(ctx->M, ctx->k_stride, i, ctx->global_L, ctx->do_trunc,
                         ctx->MM_t, ctx->MI_t, ctx->MD_t, ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                         emit_row_next, BM_next, BI_next, BM_curr, BI_curr, BD_curr);
  else
    ibv_bwd_row_b(ctx->M, k_lo, k_hi, i, ctx->global_L, ctx->do_trunc,
                  ctx->MM_t, ctx->MI_t, ctx->MD_t, ctx->IM_t, ctx->II_t, ctx->DM_t, ctx->DD_t,
                  emit_row_next, BM_next, BI_next, BM_curr, BI_curr, BD_curr);
}

static inline void
ibv_through_dispatch(IBV_DnC_Ctx *ctx, int wide, int k_lo, int k_hi,
                     const float *FM, const float *FI, const float *FD,
                     const float *BM, const float *BI, const float *BD,
                     int *ret_kmin, int *ret_kmax, int *ret_kargmax)
{
  if (wide)
    ibv_through_scan(ctx->M, ctx->k_stride, ctx->thr, ctx->ibv_mode, ctx->ibv_width,
                     FM, FI, FD, BM, BI, BD, ctx->through, ret_kmin, ret_kmax, ret_kargmax);
  else
    ibv_through_b(ctx->M, k_lo, k_hi, ctx->thr,
                  FM, FI, FD, BM, BI, BD, ret_kmin, ret_kmax, ret_kargmax);
}

/* k-banded mirror of ibv_dnc_recurse.  k_lo/k_hi bound the optimal path's
 * model column over rows [i_lo,i_hi] (monotone-k tube).  Children narrow the
 * band around the exact midline pin i2k[i_mid].  Brief 173: wide nodes route to
 * the SSE full-M primitives, narrow nodes to the scalar _b primitives. */
static void
ibv_dnc_recurse_banded(IBV_DnC_Ctx *ctx, int i_lo, int i_hi, int depth,
                       int k_lo, int k_hi,
                       const float *F_lo_M, const float *F_lo_I, const float *F_lo_D,
                       const float *B_hi_M, const float *B_hi_I, const float *B_hi_D)
{
  int    M        = ctx->M;
  size_t ks       = ctx->k_stride;
  int    global_L = ctx->global_L;
  int    base_slab= ctx->base_slab;

  if (i_hi <= i_lo) return;
  int slab_size = i_hi - i_lo;

  /* Brief 173: this node is "wide" if its band spans >= wide_thresh model
   * columns; wide nodes route to the SSE full-M primitives, narrow nodes to the
   * scalar _b primitives.  One decision per node (the band [k_lo,k_hi] is fixed
   * for all of this node's forward/backward/through streams). */
  int wide = (k_hi - k_lo + 1) >= ctx->wide_thresh;

  if (slab_size <= base_slab) {
    memcpy(ctx->slab_F + 0 * ks, F_lo_M, ks * sizeof(float));
    memcpy(ctx->slab_F + 1 * ks, F_lo_I, ks * sizeof(float));
    memcpy(ctx->slab_F + 2 * ks, F_lo_D, ks * sizeof(float));
    for (int r = 1; r <= slab_size; r++) {
      size_t po = (size_t)(r - 1) * 3 * ks, co = (size_t) r * 3 * ks;
      ibv_fwd_row_dispatch(ctx, wide, k_lo, k_hi, i_lo + r,
                    ctx->slab_F + po + 0*ks, ctx->slab_F + po + 1*ks, ctx->slab_F + po + 2*ks,
                    ctx->slab_F + co + 0*ks, ctx->slab_F + co + 1*ks, ctx->slab_F + co + 2*ks);
    }
    {
      int rh = i_hi;
      const float *em = (rh < global_L) ? ibv_emit(ctx, rh + 1) : NULL;
      size_t off = (size_t) slab_size * 3 * ks;
      ibv_bwd_row_dispatch(ctx, wide, k_lo, k_hi, rh,
                    em, B_hi_M, B_hi_I,
                    ctx->slab_B + off + 0*ks, ctx->slab_B + off + 1*ks, ctx->slab_B + off + 2*ks);
    }
    for (int r = slab_size - 1; r >= 1; r--) {
      size_t no = (size_t)(r + 1) * 3 * ks, co = (size_t) r * 3 * ks;
      ibv_bwd_row_dispatch(ctx, wide, k_lo, k_hi, i_lo + r,
                    ibv_emit(ctx, i_lo + r + 1),
                    ctx->slab_B + no + 0*ks, ctx->slab_B + no + 1*ks,
                    ctx->slab_B + co + 0*ks, ctx->slab_B + co + 1*ks, ctx->slab_B + co + 2*ks);
    }
    for (int r = 1; r <= slab_size; r++) {
      size_t off = (size_t) r * 3 * ks;
      ibv_through_dispatch(ctx, wide, k_lo, k_hi,
                    ctx->slab_F + off + 0*ks, ctx->slab_F + off + 1*ks, ctx->slab_F + off + 2*ks,
                    ctx->slab_B + off + 0*ks, ctx->slab_B + off + 1*ks, ctx->slab_B + off + 2*ks,
                    &ctx->kmin[i_lo + r], &ctx->kmax[i_lo + r], &ctx->i2k[i_lo + r]);
    }
    return;
  }

  int i_mid = (i_lo + i_hi) / 2;
  float *F_mid_M  = ctx->arena + (size_t) depth * 6 * ks + 0 * ks;
  float *F_mid_I  = ctx->arena + (size_t) depth * 6 * ks + 1 * ks;
  float *F_mid_D  = ctx->arena + (size_t) depth * 6 * ks + 2 * ks;
  float *B_mid1_M = ctx->arena + (size_t) depth * 6 * ks + 3 * ks;
  float *B_mid1_I = ctx->arena + (size_t) depth * 6 * ks + 4 * ks;
  float *B_mid1_D = ctx->arena + (size_t) depth * 6 * ks + 5 * ks;

  float *rFpM = ctx->roll_F + 0*ks, *rFpI = ctx->roll_F + 1*ks, *rFpD = ctx->roll_F + 2*ks;
  float *rFcM = ctx->roll_F + 3*ks, *rFcI = ctx->roll_F + 4*ks, *rFcD = ctx->roll_F + 5*ks;
  float *rBpM = ctx->roll_B + 0*ks, *rBpI = ctx->roll_B + 1*ks, *rBpD = ctx->roll_B + 2*ks;
  float *rBcM = ctx->roll_B + 3*ks, *rBcI = ctx->roll_B + 4*ks, *rBcD = ctx->roll_B + 5*ks;

  memcpy(rFpM, F_lo_M, ks * sizeof(float));
  memcpy(rFpI, F_lo_I, ks * sizeof(float));
  memcpy(rFpD, F_lo_D, ks * sizeof(float));
  for (int r = i_lo + 1; r < i_mid; r++) {
    ibv_fwd_row_dispatch(ctx, wide, k_lo, k_hi, r,
                  rFpM, rFpI, rFpD, rFcM, rFcI, rFcD);
    float *t; t=rFpM; rFpM=rFcM; rFcM=t; t=rFpI; rFpI=rFcI; rFcI=t; t=rFpD; rFpD=rFcD; rFcD=t;
  }
  ibv_fwd_row_dispatch(ctx, wide, k_lo, k_hi, i_mid,
                rFpM, rFpI, rFpD, F_mid_M, F_mid_I, F_mid_D);

  memcpy(rBpM, B_hi_M, ks * sizeof(float));
  memcpy(rBpI, B_hi_I, ks * sizeof(float));
  memcpy(rBpD, B_hi_D, ks * sizeof(float));
  for (int r = i_hi; r > i_mid + 1; r--) {
    const float *em = (r < global_L) ? ibv_emit(ctx, r + 1) : NULL;
    ibv_bwd_row_dispatch(ctx, wide, k_lo, k_hi, r,
                  em, rBpM, rBpI, rBcM, rBcI, rBcD);
    float *t; t=rBpM; rBpM=rBcM; rBcM=t; t=rBpI; rBpI=rBcI; rBcI=t; t=rBpD; rBpD=rBcD; rBcD=t;
  }
  {
    int r = i_mid + 1;
    const float *em = (r < global_L) ? ibv_emit(ctx, r + 1) : NULL;
    ibv_bwd_row_dispatch(ctx, wide, k_lo, k_hi, r,
                  em, rBpM, rBpI, B_mid1_M, B_mid1_I, B_mid1_D);
  }
  {
    const float *em = (i_mid < global_L) ? ibv_emit(ctx, i_mid + 1) : NULL;
    ibv_bwd_row_dispatch(ctx, wide, k_lo, k_hi, i_mid,
                  em, B_mid1_M, B_mid1_I, ctx->Bmid + 0*ks, ctx->Bmid + 1*ks, ctx->Bmid + 2*ks);
  }

  ibv_through_dispatch(ctx, wide, k_lo, k_hi,
                F_mid_M, F_mid_I, F_mid_D,
                ctx->Bmid + 0*ks, ctx->Bmid + 1*ks, ctx->Bmid + 2*ks,
                &ctx->kmin[i_mid], &ctx->kmax[i_mid], &ctx->i2k[i_mid]);

  int kmid = ctx->i2k[i_mid];
  int top_hi, bot_lo;
  if (kmid >= 1) {
    top_hi = kmid + ctx->kband_pad; if (top_hi > k_hi) top_hi = k_hi;
    bot_lo = kmid - ctx->kband_pad; if (bot_lo < k_lo) bot_lo = k_lo;
  } else {
    top_hi = k_hi; bot_lo = k_lo;   /* no emitting pin: keep parent band */
  }
  ibv_dnc_recurse_banded(ctx, i_lo, i_mid, depth + 1, k_lo, top_hi,
                         F_lo_M, F_lo_I, F_lo_D, B_mid1_M, B_mid1_I, B_mid1_D);
  ibv_dnc_recurse_banded(ctx, i_mid, i_hi, depth + 1, bot_lo, k_hi,
                         F_mid_M, F_mid_I, F_mid_D, B_hi_M, B_hi_I, B_hi_D);
}


int
p7_Seq2BandsIBV_dnc(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L,
                    int delta_milli, int base_slab,
                    int do_boundary_widen,
                    int do_kband,
                    int do_trunc,
                    int ibv_mode, int ibv_width,
                    int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int          status;
  P7_HMM      *hmm   = NULL;
  int          M, K, k, i;
  float        optimal, thr, floor_milli;
  int          ncells = 0;
  size_t       ks;
  int          max_depth;
  IBV_DnC_Ctx  ctx;
  float       *begin_milli = NULL;
  float       *emit_pool  = NULL;
  float       *F_row0_M   = NULL, *F_row0_I   = NULL, *F_row0_D   = NULL;
  float       *B_seed_M   = NULL, *B_seed_I   = NULL, *B_seed_D   = NULL;
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
    double per_row   = 2.0 * 3.0 * (double) ks * 4.0;   /* F+B slab, one row */
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
  ESL_ALLOC(ctx.emit_table, sizeof(float *) * (K + 1));
  for (int xt = 0; xt <= K; xt++) ctx.emit_table[xt] = emit_pool + (size_t) xt * ks;
  for (int xt = 0; xt < K; xt++) {
    float *row = ctx.emit_table[xt];
    for (k = 0; k <= M; k++) row[k] = p7ibv_emit_milli(hmm, k, xt);
    for (k = M + 1; k < (int) ks; k++) row[k] = P7IBV_NEG_INF;
  }
  for (k = 0; k < (int) ks; k++) ctx.emit_table[K][k] = 0.0f;

  /* Brief 171: Tgm begin-anywhere vector (see p7_Seq2BandsIBV). */
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
      begin_milli[k] = (b > 0.0) ? (float)(P7IBV_INTSCALE * (log(b) / M_LN2)) : P7IBV_NEG_INF;
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
      float a = F_row0_M[k - 1] + ctx.MD_t[k - 1];
      float b = F_row0_D[k - 1] + ctx.DD_t[k - 1];
      F_row0_D[k] = (a > b) ? a : b;
    }
  }

  /* Global forward pass (2-row rolling) to get optimal score.
   * Reuse roll_F for this; reset to NEG_INF afterwards. */
  if ((status = ibv_dnc_alloc(6 * ks, &ctx.roll_F)) != eslOK) goto ERROR;
  {
    float *rM = ctx.roll_F + 0*ks, *rI = ctx.roll_F + 1*ks, *rD = ctx.roll_F + 2*ks;
    float *cM = ctx.roll_F + 3*ks, *cI = ctx.roll_F + 4*ks, *cD = ctx.roll_F + 5*ks;
    memcpy(rM, F_row0_M, ks * sizeof(float));
    memcpy(rI, F_row0_I, ks * sizeof(float));
    memcpy(rD, F_row0_D, ks * sizeof(float));
    for (i = 1; i <= L; i++) {
      int x = (int) dsq[i];
      int xt = (x >= 0 && x < K) ? x : K;
      const float *brow = (do_trunc && i == 1) ? begin_milli : NULL;
      ibv_forward_one_row(M, ks,
                          ctx.MM_t, ctx.MI_t, ctx.MD_t,
                          ctx.IM_t, ctx.II_t, ctx.DM_t, ctx.DD_t,
                          ctx.emit_table[xt], brow,
                          rM, rI, rD, cM, cI, cD);
      float *t; t=rM; rM=cM; cM=t; t=rI; rI=cI; cI=t; t=rD; rD=cD; cD=t;
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
  floor_milli = (float) delta_milli;
  if (floor_milli < P7IBV_EPS) floor_milli = P7IBV_EPS;
  {
    float guard = floor_milli * (float)(1u << P7IBV_ULP_SAFETY_SHIFT);
    if (guard < P7IBV_OPTIMAL_SANITY) guard = P7IBV_OPTIMAL_SANITY;
    if (fabsf(optimal) > guard)
      ESL_FAIL(eslEINVAL, errbuf,
               "p7_Seq2BandsIBV_dnc: |optimal|=%g milli-bits exceeds guard %g "
               "(float ULP approaches delta=%d margin; need double precision)",
               (double) optimal, (double) guard, delta_milli);
  }
  thr = optimal - floor_milli;
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
  ctx.thr      = thr;
  ctx.do_trunc    = do_trunc;     /* brief 171 */
  ctx.begin_milli = begin_milli;  /* brief 171 */
  ctx.ibv_mode = ibv_mode;   /* brief 140 */
  ctx.ibv_width= ibv_width;  /* brief 140 */
  ctx.kband_pad = P7IBV_KPAD;  /* brief 172 */
  { const char *kp = getenv("P7IBV_KBAND_PAD");
    if (kp && *kp) { int v = atoi(kp); if (v >= 0) ctx.kband_pad = v; } }
  /* Brief 173 Part B: route nodes whose band width >= wide_thresh to the SSE
   * full-M primitives.  Crossover ~ M/4 (SSE is ~4x/cell, so full-M SSE beats
   * scalar-over-band once the band exceeds ~M/4 columns); default M/2 is a
   * safe, slightly-conservative start.  P7WV_WIDE_THRESH overrides (absolute
   * # columns) for tuning; a huge value (> M) disables SSE routing entirely. */
  ctx.wide_thresh = (M / 2 > 1) ? (M / 2) : 1;
  { const char *wt = getenv("P7WV_WIDE_THRESH");
    if (wt && *wt) { int v = atoi(wt); if (v >= 1) ctx.wide_thresh = v; } }

  /* B seed for top-level: all NEG_INF; terminal injected via global_L. */
  if ((status = ibv_dnc_alloc(ks, &B_seed_M)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(ks, &B_seed_I)) != eslOK) goto ERROR;
  if ((status = ibv_dnc_alloc(ks, &B_seed_D)) != eslOK) goto ERROR;

  /* Run recursion.  Brief 172: with do_kband, band each node's streams to the
   * monotone-k tube (exact for i2k; the windowed-Viterbi genome speed lever).
   * Top band is the full model [0,M] (begin/end anywhere). */
  if (do_kband)
    ibv_dnc_recurse_banded(&ctx, 0, L, 0, 0, M,
                           F_row0_M, F_row0_I, F_row0_D,
                           B_seed_M, B_seed_I, B_seed_D);
  else
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
  kmin_arr[0] = 0; kmax_arr[0] = 0;
  i2k[0] = 0;   /* B-state convention (brief 137): i2k[i]=argmax_k for i in [1,L]. */

  /* Brief 140a: bridge inter-row gaps so FIXED/HYBRID bands are connected
   * (DELTA untouched).  Must run before ncells is summed. */
  ibv_connectivity_guard(L, M, ibv_mode, kmin_arr, kmax_arr);

  for (i = 1; i <= L; i++)
    ncells += (kmax_arr[i] - kmin_arr[i] + 1);

  {
    const char *dump = getenv("P7IBV_DUMP_BAND");
    if (dump != NULL && *dump != '\0') {
      FILE *fp = fopen(dump, "w");
      if (fp != NULL) {
        fprintf(fp, "# M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d\n",
                M, L, delta_milli, (double) optimal, (double) thr, ncells);
        fprintf(fp, "# i\tkmin\tkmax\twidth\ti2k\n");   /* brief 142: + i2k (argmax-k pin) */
        for (i = 1; i <= L; i++)
          fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", i, kmin_arr[i], kmax_arr[i], kmax_arr[i] - kmin_arr[i] + 1, i2k[i]);
        fclose(fp);
      }
    }
  }

  /* Brief 142: per-row band dump to stderr (multi-seq safe; one block/seq). */
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
 * Incept:    brief 137, 2026-06-17.
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
p7_IBVPins2Trace(const P7_PROFILE *gm, const ESL_DSQ *dsq, int L,
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


/* ---------------------------------------------------------------------------
 * cm_ComputeP7WVNodePad -- brief 169 per-node pad calibration (F+B-halfwidth p95)
 * ---------------------------------------------------------------------------
 *
 * The windowed-Viterbi band uses a per-node pad calibrated as the <quantile>
 * (default p95) of the F+B Delta-band half-width observed at the Viterbi MAP
 * trace cell, over a Monte-Carlo sample of CM-emitted sequences.  This is the
 * brief-168 prototype's calibrate_pernode_hw / pad_from_hw, ported to C:
 *
 *   for s in 1..nsamples:
 *     emit a sequence from the CM (EmitParsetree)
 *     run the flat F+B IBV deriver at <delta_milli> -> (i2k, kmin, kmax)
 *     for each residue i with c = i2k[i] >= 1:
 *       hw[c].append( max(c-kmin[i], kmax[i]-c) )       # the F+B half-width
 *   pad[k] = max(floorpad, ceil( quantile(hw[k]) ))     # linear interp, numpy-style
 *
 * IMPORTANT: this is a DIFFERENT calibration from cm->p7_cm_nodepad
 * (cm_ComputeP7CMNodePad), which is the p99 of a *Viterbi-pin deficit* against
 * the embedded true alignment -- a much narrower quantity tuned for the
 * p7_Seq2BandsVit pin band.  Using cm->p7_cm_nodepad here regresses accuracy
 * (brief 169 spot-check: MISL -0.19, Bp1 -0.06 vs F+B); the wider F+B-halfwidth
 * pad reproduces F+B accuracy (the prototype's pn_p95 result).
 *
 * Caller owns the returned pad[0..M].  Align-time use: compute once per CM and
 * cache; works on existing CMs with no rebuild (only needs cm->fp7).
 */
static int
p7wv_cmp_int(const void *a, const void *b)
{
  int x = *(const int *) a, y = *(const int *) b;
  return (x > y) - (x < y);
}

int
cm_ComputeP7WVNodePad(CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, int nsamples,
                      double quantile, int delta_milli, int floorpad,
                      int **ret_nodepad)
{
  int       status;
  P7_HMM   *hmm = NULL;
  int       M, k, s, i;
  int     **hw  = NULL;   /* hw[k] = collected half-widths at node k */
  int      *hwn = NULL;   /* hw[k] count */
  int      *hwa = NULL;   /* hw[k] alloc */
  int      *pad = NULL;
  /* brief 171: calibrate with the same begin/end semantics the align-time band
   * will use, so the F+B half-widths match the deployed deriver. */
  int       do_trunc = (cm->align_opts & CM_ALIGN_TRUNC) ? TRUE : FALSE;

  if (cm == NULL || cm->fp7 == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "cm_ComputeP7WVNodePad: cm->fp7 is NULL");
  if (r == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "cm_ComputeP7WVNodePad: RNG is NULL");
  hmm = cm->fp7;
  M   = hmm->M;
  if (floorpad < 0) floorpad = 0;

  ESL_ALLOC(hw,  sizeof(int *) * (M + 1));
  ESL_ALLOC(hwn, sizeof(int)   * (M + 1));
  ESL_ALLOC(hwa, sizeof(int)   * (M + 1));
  for (k = 0; k <= M; k++) { hw[k] = NULL; hwn[k] = 0; hwa[k] = 0; }

  for (s = 0; s < nsamples; s++) {
    Parsetree_t *tr  = NULL;
    ESL_SQ      *esq = NULL;
    int          L   = 0;
    char         name[32];
    int         *i2k = NULL, *kmin = NULL, *kmax = NULL, nc = 0;

    snprintf(name, sizeof(name), "wv%d", s);
    if ((status = EmitParsetree(cm, errbuf, r, name, TRUE, &tr, &esq, &L)) != eslOK) goto ERROR;

    /* Brief 172: the pad calibration measures the FULL-M delta-band half-width,
     * so it must NOT use the k-banded deriver (whose delta cloud is tube-clipped).
     * Use the flat deriver when its O(L*M) pool fits ~2 GB (small/viral emits),
     * else the memory-bounded UNbanded D&C (do_kband=FALSE) -- correct delta
     * cloud, genome-safe.  (Flat would need ~466 GB at genome-length emits.) */
    int wv_derr = eslFAIL;
    if (L >= 3) {
      double pool = 12.0 * (double)(L + 1) * (double)(M + 4);
      /* P7WV_FAST_CALIB: use the k-banded D&C for the calibration deriver at
       * genome scale.  For typical (glocal/non-truncated) emits the KPAD=64
       * monotone tube is wider than the Delta cloud, so the banded delta band
       * equals the full one and the pad is unchanged -- but the per-emit deriver
       * is ~order(s) faster, making genome-scale calibration tractable.  Off by
       * default (exact full-cloud pad). */
      int calib_kband = FALSE;
      { const char *fc = getenv("P7WV_FAST_CALIB"); if (fc && *fc && *fc != '0') calib_kband = TRUE; }
      if (pool <= 2.0e9)
        wv_derr = p7_Seq2BandsIBV(cm, errbuf, esq->dsq, L, delta_milli, do_trunc,
                                  P7IBV_MODE_DELTA, 0, &i2k, &kmin, &kmax, &nc);
      else
        wv_derr = p7_Seq2BandsIBV_dnc(cm, errbuf, esq->dsq, L, delta_milli, 0,
                                      FALSE, calib_kband, do_trunc, P7IBV_MODE_DELTA, 0,
                                      &i2k, &kmin, &kmax, &nc);
    }
    if (wv_derr == eslOK) {
      for (i = 1; i <= L; i++) {
        int c = i2k[i];
        int d1, d2, h;
        if (c < 1 || c > M) continue;
        d1 = c - kmin[i];
        d2 = kmax[i] - c;
        h  = (d1 > d2) ? d1 : d2;
        if (h < 0) h = 0;
        if (hwn[c] >= hwa[c]) {
          int ns = hwa[c] ? hwa[c] * 2 : 8;
          ESL_REALLOC(hw[c], sizeof(int) * ns);
          hwa[c] = ns;
        }
        hw[c][hwn[c]++] = h;
      }
      free(i2k); free(kmin); free(kmax);
    }
    FreeParsetree(tr);
    esl_sq_Destroy(esq);
  }

  ESL_ALLOC(pad, sizeof(int) * (M + 1));
  pad[0] = 0;
  for (k = 1; k <= M; k++) {
    if (hwn[k] > 0) {
      double idx, frac;
      int    lo, hi, v;
      qsort(hw[k], hwn[k], sizeof(int), p7wv_cmp_int);
      idx  = quantile * (double)(hwn[k] - 1);   /* numpy 'linear' percentile */
      lo   = (int) floor(idx);
      hi   = lo + 1;
      frac = idx - (double) lo;
      if (hi >= hwn[k]) v = hw[k][hwn[k] - 1];
      else              v = (int) ceil((double) hw[k][lo] + frac * (double)(hw[k][hi] - hw[k][lo]));
      pad[k] = (v > floorpad) ? v : floorpad;
    } else {
      pad[k] = floorpad;
    }
  }

  for (k = 0; k <= M; k++) if (hw[k]) free(hw[k]);
  free(hw); free(hwn); free(hwa);
  *ret_nodepad = pad;
  return eslOK;

 ERROR:
  if (hw)  { for (k = 0; k <= M; k++) if (hw[k]) free(hw[k]); free(hw); }
  if (hwn) free(hwn);
  if (hwa) free(hwa);
  if (pad) free(pad);
  *ret_nodepad = NULL;
  return status;
}


/* ---------------------------------------------------------------------------
 * p7_Seq2BandsWV -- brief 169 windowed-Viterbi band deriver
 * ---------------------------------------------------------------------------
 *
 * The windowed-Viterbi band (brief 168 prototype, GO verdict) is, by
 * construction:
 *
 *     band[i] = [ i2k[i] - nodepad[i2k[i]] ,  i2k[i] + nodepad[i2k[i]] ]
 *
 * where i2k[] is the Viterbi MAP trace (the model column the optimal path
 * occupies at residue i) and nodepad[k] is the per-node pad.  The prototype
 * proved this reproduces the F+B Delta-band's alignment accuracy on 183/189
 * sequences (96.8%; exact on all rmark + dossier) -- see brief 168 summary.
 *
 * Two facts make this a thin composition of existing, validated machinery:
 *   (a) The MAP trace i2k is exactly the per-row argmax-k pin the IBV deriver
 *       already returns (ibv_through_scan restricts the argmax to emitting
 *       M/I cells, brief 137).
 *   (b) The per-node pad is exactly cm->p7_cm_nodepad -- cm_ComputeP7CMNodePad
 *       calibrates it by the same emit-from-CM Monte-Carlo deficit-quantile
 *       method the prototype reinvented (calibrate_pernode_hw / pad_from_hw).
 *       It is stored on the CM file (CMH_P7NODEPAD), so it works on existing
 *       Rfam/VADR CMs with NO rebuild.
 *   The band itself is then the same i2k+nodepad construction the
 *   --p7pinbridge / vitband paths already use via p7_pins2bands_nodepad
 *   (D-state bridge / connectivity guard built in).
 *
 * So the ONLY genuinely new lever is *how fast* and *at what memory* we obtain
 * i2k.  This v1 obtains i2k from the existing memory-bounded D&C deriver
 * (p7_Seq2BandsIBV_dnc) -- correct, genome-capable, accuracy-validatable, but
 * NOT yet faster than --p7ibv-mem (same ~35-sweep i2k cost).  The brief-169
 * speed win (single windowed forward-Viterbi + traceback for i2k, ~3 sweeps or
 * less) is a drop-in replacement for the i2k source below; it must reproduce
 * this i2k byte-for-byte (the monotone-k trace is unique up to float ties).
 *
 * <nodepad> is the caller-owned [0..M] per-node pad array (typically
 * cm->p7_cm_nodepad[k] + cm->p7bpad).  Returns (i2k, kmin, kmax, ncells) with
 * the same conventions as p7_Seq2BandsIBV.
 */
int
p7_Seq2BandsWV(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L, int *nodepad,
               int do_trunc,
               int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int   status;
  int   M;
  int  *i2k = NULL, *kmin_tmp = NULL, *kmax_tmp = NULL;
  int  *i2k_band = NULL, *kmin = NULL, *kmax = NULL;
  int   nc_tmp = 0, ncells = 0;

  if (cm == NULL || cm->fp7 == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsWV: cm->fp7 is NULL");
  if (nodepad == NULL)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsWV: nodepad is NULL (CM lacks P7NODEPAD?)");
  M = cm->fp7->M;

  /* (1) Exact MAP trace i2k.  We keep only i2k (the per-row argmax pin); the
   *     Delta band is discarded.  i2k is threshold-independent so the delta
   *     value does not matter.
   *
   *     Speed/memory: the FLAT deriver is ~2 sweeps (vs the D&C's ~35) but
   *     needs an O(L*M) forward pool (~12*L*M bytes).  Use it when that pool
   *     fits a ~2 GB budget (covers small/viral), else fall back to the
   *     memory-bounded D&C (genome).  This delivers the de-recursion speedup
   *     wherever memory allows; the genome-scale windowed forward-Viterbi
   *     (single bounded sweep) is the remaining brief-169 speed lever and is a
   *     drop-in replacement for this block (must reproduce this i2k). */
  {
    double flat_pool_bytes = 12.0 * (double)(L + 1) * (double)(M + 4);
    double FLAT_BUDGET = 2.0e9;
    /* Brief 172: P7WV_FORCE_KBAND forces the k-banded D&C i2k path even when the
     * flat pool would fit, so the banded kernel can be validated (exact i2k vs
     * the unbanded D&C oracle) on small/viral seqs via test_wviterbi --wv -c. */
    { const char *fk = getenv("P7WV_FORCE_KBAND"); if (fk && *fk && *fk != '0') FLAT_BUDGET = 0.0; }
    if (flat_pool_bytes <= FLAT_BUDGET)
      status = p7_Seq2BandsIBV(cm, errbuf, dsq, L, cm->p7_ibv_delta, do_trunc,
                               P7IBV_MODE_DELTA, 0,
                               &i2k, &kmin_tmp, &kmax_tmp, &nc_tmp);
    else
      /* Brief 172: genome scale -- k-banded D&C (do_kband=TRUE) restricts each
       * node's streams to the monotone-k tube, collapsing the unbanded D&C's
       * ~log L full-M sweeps to ~1-2 sweeps of work.  Exact for i2k. */
      status = p7_Seq2BandsIBV_dnc(cm, errbuf, dsq, L,
                                   cm->p7_ibv_delta, cm->p7_ibv_base_slab,
                                   FALSE, TRUE, do_trunc, P7IBV_MODE_DELTA, 0,
                                   &i2k, &kmin_tmp, &kmax_tmp, &nc_tmp);
  }
  if (status != eslOK)
    return status;
  free(kmin_tmp); kmin_tmp = NULL;
  free(kmax_tmp); kmax_tmp = NULL;

  /* (2) Windowed-Viterbi band = i2k +/- nodepad.  p7_pins2bands_nodepad prunes
   *     i2k in place (non-monotone pins), so build the band from a copy and
   *     return the unpruned i2k to the caller. */
  ESL_ALLOC(i2k_band, sizeof(int) * (L + 1));
  memcpy(i2k_band, i2k, sizeof(int) * (L + 1));
  if ((status = p7_pins2bands_nodepad(i2k_band, errbuf, L, M, nodepad, 0,
                                      &kmin, &kmax, &ncells)) != eslOK)
    goto ERROR;
  free(i2k_band); i2k_band = NULL;

  *ret_i2k    = i2k;
  *ret_kmin   = kmin;
  *ret_kmax   = kmax;
  *ret_ncells = ncells;
  return eslOK;

 ERROR:
  if (i2k)      free(i2k);
  if (i2k_band) free(i2k_band);
  if (kmin)     free(kmin);
  if (kmax)     free(kmax);
  if (kmin_tmp) free(kmin_tmp);
  if (kmax_tmp) free(kmax_tmp);
  *ret_i2k = NULL; *ret_kmin = NULL; *ret_kmax = NULL; *ret_ncells = 0;
  return status;
}

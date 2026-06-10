/* p7_ibv.c -- F+B direct-band band derivation for cmalign --p7ibv
 *
 * Brief 121 C3: F-stored, B-streamed memory pattern.
 *   - F_M / F_I / F_D : stored for all rows 0..L (needed by inline
 *     through-score scan when B[i] is computed during backward sweep).
 *   - B_M / B_I / B_D : 2-row rolling buffer per state (prev = i+1, curr = i).
 *     6 single-row buffers of k_stride floats each.
 *   - Through-score scan + per-row kmin/kmax extraction integrated into
 *     the backward sweep, right after B[i] is computed. Allows B[i+1]
 *     to be discarded as soon as B[i] is produced.
 *
 * Memory budget (LSU M=3400 L=2771, k_stride=3408, float):
 *   - F_M + F_I + F_D       : 3 * 2772 * 3408 * 4 ~= 113 MB  (vs C2 = 226 MB)
 *   - 6 B rolling rows      : 6 * 3408 * 4         ~= 82 KB
 *   - emit_table (K+1 rows) : 5 * 3408 * 4         ~= 68 KB
 *   - through scratch       : 3408 * 4             ~= 14 KB
 *   - kmin/kmax/i2k arrays  : 3 * 2772 * 4         ~= 33 KB
 *   - total                 : ~113 MB              (vs C2 ~226 MB, 2x reduction)
 *
 * SSE M+I along k, scalar D-fill (unchanged from C2):
 *   - Forward: scalar prefix k=0..3 (M, I, D); SSE bulk k=4..k_sse_end-1
 *     for M and I; scalar tail; scalar left-to-right D-fill k=4..M.
 *   - Backward: scalar right-to-left D-fill k=M..0; SSE bulk M and I
 *     from k=0..k_sse_end-1; scalar tail.
 *
 * Padded tail cells (k > M) MUST stay -INF so subsequent SSE loads at
 * row edges fold harmlessly. Backward rolling buffers initialize their
 * tails to -INF once at allocation, then are reused without re-init
 * because SSE writes stay within [0, M] and the scalar D-fill writes
 * k=0..M-1 and explicit -INF for k=M.
 *
 * ULP analysis (float, INTSCALE=1000 milli-bits):
 *   Cumulative path scores in practice are O(1e5) milli-bits.
 *   Float ULP at 3e5 is ~0.04 milli-bit; EPS=1 absorbs 25x margin.
 *   Worst case (cumulative ~1e7): ULP ~0.6 milli-bit; EPS=1 still 1.5x
 *   above; tight but safe. Asserted optimal < 1e7 as a sanity guard.
 *
 * Algorithm (per brief 116-117), unchanged from brief 120:
 *   F_M[i,k] = emit(k, seq[i]) + max(F_M[i-1,k-1]+T_MM[k-1],
 *                                    F_I[i-1,k-1]+T_IM[k-1],
 *                                    F_D[i-1,k-1]+T_DM[k-1])
 *   F_I[i,k] = max(F_M[i-1,k]+T_MI[k], F_I[i-1,k]+T_II[k])
 *   F_D[i,k] = max(F_M[i,k-1]+T_MD[k-1], F_D[i,k-1]+T_DD[k-1])
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
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

#define P7IBV_INTSCALE        1000.0f
#define P7IBV_NEG_INF         (-1.0e18f)
#define P7IBV_HALF_NEG_INF    (-5.0e17f)
#define P7IBV_EPS             1.0f
#define P7IBV_K_ALIGN         16
#define P7IBV_OPTIMAL_SANITY  1.0e7f

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


int
p7_Seq2BandsIBV(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L, int delta_milli,
                int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int       status;
  P7_HMM   *hmm = NULL;
  int       M;
  int       i, k;
  int       K;
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

  k_stride   = ((size_t)(M + 1) + (P7IBV_K_ALIGN - 1)) & ~(size_t)(P7IBV_K_ALIGN - 1);
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

  /* F: stored, full (L+1) rows per state. */
  if ((status = ibv_alloc_floats(pool_cells, &FM_pool)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(pool_cells, &FI_pool)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(pool_cells, &FD_pool)) != eslOK) goto ERROR;
  for (size_t c = 0; c < pool_cells; c++)
    FM_pool[c] = FI_pool[c] = FD_pool[c] = P7IBV_NEG_INF;

  /* B: 2-row rolling buffers per state. _a and _b swap roles each backward iter. */
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

  /* Through-score scratch (1 row). */
  if ((status = ibv_alloc_floats(k_stride, &through)) != eslOK) goto ERROR;

#define F_M(i)  (FM_pool + (size_t)(i) * k_stride)
#define F_I(i)  (FI_pool + (size_t)(i) * k_stride)
#define F_D(i)  (FD_pool + (size_t)(i) * k_stride)

  /* ---------- Forward (unchanged from C2) ---------- */
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

  for (i = 1; i <= L; i++) {
    int   x      = (int) dsq[i];
    int   xt     = (x >= 0 && x < K) ? x : K;
    float *emit_row = emit_table[xt];
    float *fm_i     = F_M(i);
    float *fi_i     = F_I(i);
    float *fd_i     = F_D(i);
    float *fm_im1   = F_M(i - 1);
    float *fi_im1   = F_I(i - 1);
    float *fd_im1   = F_D(i - 1);

    int kpref_end = (M < 3) ? M : 3;
    for (k = 0; k <= kpref_end; k++) {
      float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF, cD = P7IBV_NEG_INF;
      {
        float a = fm_im1[k] + MI_t[k];
        float b = fi_im1[k] + II_t[k];
        cI = (a > b) ? a : b;
      }
      if (k >= 1) {
        float a = fm_im1[k - 1] + MM_t[k - 1];
        float b = fi_im1[k - 1] + IM_t[k - 1];
        float c = fd_im1[k - 1] + DM_t[k - 1];
        float m = (a > b) ? a : b;
        if (c > m) m = c;
        cM = m + emit_row[k];
        float a2 = fm_i[k - 1] + MD_t[k - 1];
        float b2 = fd_i[k - 1] + DD_t[k - 1];
        cD = (a2 > b2) ? a2 : b2;
      }
      fm_i[k] = cM;
      fi_i[k] = cI;
      fd_i[k] = cD;
    }

    int k_sse_start = 4;
    int k_sse_end   = k_sse_start;
    while (k_sse_end + 3 <= M) k_sse_end += 4;
    for (k = k_sse_start; k < k_sse_end; k += 4) {
      __m128 m_prev = _mm_loadu_ps(&fm_im1[k - 1]);
      __m128 i_prev = _mm_loadu_ps(&fi_im1[k - 1]);
      __m128 d_prev = _mm_loadu_ps(&fd_im1[k - 1]);
      __m128 t_mm   = _mm_loadu_ps(&MM_t[k - 1]);
      __m128 t_im   = _mm_loadu_ps(&IM_t[k - 1]);
      __m128 t_dm   = _mm_loadu_ps(&DM_t[k - 1]);
      __m128 a      = _mm_add_ps(m_prev, t_mm);
      __m128 b      = _mm_add_ps(i_prev, t_im);
      __m128 c      = _mm_add_ps(d_prev, t_dm);
      __m128 mx     = p7ibv_mm_max3(a, b, c);
      __m128 e_vec  = _mm_loadu_ps(&emit_row[k]);
      _mm_storeu_ps(&fm_i[k], _mm_add_ps(mx, e_vec));

      __m128 fm_k = _mm_loadu_ps(&fm_im1[k]);
      __m128 fi_k = _mm_loadu_ps(&fi_im1[k]);
      __m128 t_mi = _mm_loadu_ps(&MI_t[k]);
      __m128 t_ii = _mm_loadu_ps(&II_t[k]);
      __m128 a2   = _mm_add_ps(fm_k, t_mi);
      __m128 b2   = _mm_add_ps(fi_k, t_ii);
      _mm_storeu_ps(&fi_i[k], _mm_max_ps(a2, b2));
    }

    for (k = k_sse_end; k <= M; k++) {
      float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
      {
        float a = fm_im1[k] + MI_t[k];
        float b = fi_im1[k] + II_t[k];
        cI = (a > b) ? a : b;
      }
      if (k >= 1) {
        float a = fm_im1[k - 1] + MM_t[k - 1];
        float b = fi_im1[k - 1] + IM_t[k - 1];
        float c = fd_im1[k - 1] + DM_t[k - 1];
        float m = (a > b) ? a : b;
        if (c > m) m = c;
        cM = m + emit_row[k];
      }
      fm_i[k] = cM;
      fi_i[k] = cI;
    }

    for (k = (kpref_end + 1 > M) ? (M + 1) : (kpref_end + 1); k <= M; k++) {
      float a = fm_i[k - 1] + MD_t[k - 1];
      float b = fd_i[k - 1] + DD_t[k - 1];
      fd_i[k] = (a > b) ? a : b;
    }
  }

  /* ---------- optimal + threshold (between forward and backward) ---------- */
  {
    float *fm_L = F_M(L);
    float *fi_L = F_I(L);
    float *fd_L = F_D(L);
    optimal = fm_L[M];
    if (fi_L[M] > optimal) optimal = fi_L[M];
    if (fd_L[M] > optimal) optimal = fd_L[M];
  }
  assert(optimal == optimal);
  if (fabsf(optimal) > P7IBV_OPTIMAL_SANITY)
    ESL_FAIL(eslEINVAL, errbuf,
             "p7_Seq2BandsIBV: |optimal|=%g milli-bits exceeds %g (float ULP > EPS); "
             "consider double-precision fallback",
             (double) optimal, (double) P7IBV_OPTIMAL_SANITY);
  floor_milli = (float) delta_milli;
  if (floor_milli < P7IBV_EPS) floor_milli = P7IBV_EPS;
  thr = optimal - floor_milli;

  /* Output arrays. */
  ESL_ALLOC(i2k,  sizeof(int) * (L + 1));
  ESL_ALLOC(kmin, sizeof(int) * (L + 1));
  ESL_ALLOC(kmax, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, L + 1, -1);
  /* Default per-row: no admitted cells -> [1, M] fallback (matches scalar). */
  for (i = 0; i <= L; i++) { kmin[i] = 1; kmax[i] = M; }

  /* ---------- Backward (streamed) + inline through-scan ---------- */
  /* Seed row i=L. */
  float *B_M_prev = BM_a;
  float *B_I_prev = BI_a;
  float *B_D_prev = BD_a;
  float *B_M_curr = BM_b;
  float *B_I_curr = BI_b;
  float *B_D_curr = BD_b;

  /* B[L]: seed at (L,M) = 0, deletion cascade for M and D, I = -INF. */
  B_M_prev[M] = 0.0f;
  B_I_prev[M] = 0.0f;
  B_D_prev[M] = 0.0f;
  for (k = M - 1; k >= 0; k--) {
    float bv = B_D_prev[k + 1];
    B_M_prev[k] = MD_t[k] + bv;
    B_D_prev[k] = DD_t[k] + bv;
    B_I_prev[k] = P7IBV_NEG_INF;
  }
  /* Row i=L band: boundary-widened to [1, M] (already set above); skip scan. */

  /* Iterate i = L-1 downto 0. */
  for (i = L - 1; i >= 0; i--) {
    int   x_next  = (int) dsq[i + 1];
    int   xt      = (x_next >= 0 && x_next < K) ? x_next : K;
    float *emit_row_next = emit_table[xt];
    float *bm_ip1 = B_M_prev;
    float *bi_ip1 = B_I_prev;
    float *bd_ip1 = B_D_prev;  /* not used directly; in-row D comes from B_D_curr */
    float *bm_i   = B_M_curr;
    float *bi_i   = B_I_curr;
    float *bd_i   = B_D_curr;
    (void) bd_ip1;

    /* D scalar fill, right-to-left, k=M..0. */
    bd_i[M] = P7IBV_NEG_INF;
    for (k = M - 1; k >= 0; k--) {
      float bv_m = bm_ip1[k + 1] + emit_row_next[k + 1];
      float bv_d = bd_i[k + 1];
      float a = DM_t[k] + bv_m;
      float b = DD_t[k] + bv_d;
      bd_i[k] = (a > b) ? a : b;
    }

    /* SSE bulk M+I, k=0..k_sse_end-1. */
    int k_sse_end = 0;
    while (k_sse_end + 3 <= M) k_sse_end += 4;
    for (k = 0; k < k_sse_end; k += 4) {
      __m128 bm_n   = _mm_loadu_ps(&bm_ip1[k + 1]);
      __m128 e_n    = _mm_loadu_ps(&emit_row_next[k + 1]);
      __m128 bv_M   = _mm_add_ps(bm_n, e_n);
      __m128 bv_I   = _mm_loadu_ps(&bi_ip1[k]);
      __m128 bv_D   = _mm_loadu_ps(&bd_i[k + 1]);

      __m128 t_mm   = _mm_loadu_ps(&MM_t[k]);
      __m128 t_mi   = _mm_loadu_ps(&MI_t[k]);
      __m128 t_md   = _mm_loadu_ps(&MD_t[k]);
      __m128 cM_out = p7ibv_mm_max3(_mm_add_ps(t_mm, bv_M),
                                    _mm_add_ps(t_mi, bv_I),
                                    _mm_add_ps(t_md, bv_D));
      _mm_storeu_ps(&bm_i[k], cM_out);

      __m128 t_im   = _mm_loadu_ps(&IM_t[k]);
      __m128 t_ii   = _mm_loadu_ps(&II_t[k]);
      __m128 cI_out = _mm_max_ps(_mm_add_ps(t_im, bv_M),
                                 _mm_add_ps(t_ii, bv_I));
      _mm_storeu_ps(&bi_i[k], cI_out);
    }
    /* Scalar tail M+I, k=k_sse_end..M. */
    for (k = k_sse_end; k <= M; k++) {
      float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
      if (k + 1 <= M) {
        float bv_m = bm_ip1[k + 1] + emit_row_next[k + 1];
        float a = MM_t[k] + bv_m;
        float b = IM_t[k] + bv_m;
        if (a > cM) cM = a;
        if (b > cI) cI = b;
      }
      {
        float bv_i = bi_ip1[k];
        float a = MI_t[k] + bv_i;
        float b = II_t[k] + bv_i;
        if (a > cM) cM = a;
        if (b > cI) cI = b;
      }
      if (k + 1 <= M) {
        float bv_d = bd_i[k + 1];
        float a = MD_t[k] + bv_d;
        if (a > cM) cM = a;
      }
      bm_i[k] = cM;
      bi_i[k] = cI;
    }

    /* Through-score scan for row i. Skip rows that are boundary-widened
     * (i in {1, L-1}) or below the band-emission range (i == 0). Row i=L
     * was handled above. The kmin/kmax for skipped rows stay at their
     * default [1, M] (or [0, 0] for i=0; reset below).
     */
    if (i >= 2 && i <= L - 2) {
      float *fm = F_M(i);
      float *fi = F_I(i);
      float *fd = F_D(i);
      /* SSE through[k] = max(fm+bm, fi+bi, fd+bd) along k. */
      int k_thru_end = 0;
      while (k_thru_end + 3 <= M) k_thru_end += 4;
      for (k = 0; k < k_thru_end; k += 4) {
        __m128 a = _mm_add_ps(_mm_loadu_ps(&fm[k]), _mm_loadu_ps(&bm_i[k]));
        __m128 b = _mm_add_ps(_mm_loadu_ps(&fi[k]), _mm_loadu_ps(&bi_i[k]));
        __m128 c = _mm_add_ps(_mm_loadu_ps(&fd[k]), _mm_loadu_ps(&bd_i[k]));
        _mm_storeu_ps(&through[k], p7ibv_mm_max3(a, b, c));
      }
      for (k = k_thru_end; k <= M; k++) {
        float t_m = fm[k] + bm_i[k];
        float t_i = fi[k] + bi_i[k];
        float t_d = fd[k] + bd_i[k];
        float t = t_m;
        if (t_i > t) t = t_i;
        if (t_d > t) t = t_d;
        through[k] = t;
      }
      /* Scalar scan for kmin/kmax in [1, M]. */
      int row_kmin = -1, row_kmax = -1;
      for (k = 1; k <= M; k++) {
        float t = through[k];
        if (t < P7IBV_HALF_NEG_INF) continue;
        if (t >= thr) {
          if (row_kmin < 0) row_kmin = k;
          row_kmax = k;
        }
      }
      if (row_kmin < 0) {
        kmin[i] = 1; kmax[i] = M;
      } else {
        kmin[i] = row_kmin;
        kmax[i] = row_kmax;
      }
    }

    /* Swap prev <-> curr for next iteration. */
    float *t_M = B_M_prev; B_M_prev = B_M_curr; B_M_curr = t_M;
    float *t_I = B_I_prev; B_I_prev = B_I_curr; B_I_curr = t_I;
    float *t_D = B_D_prev; B_D_prev = B_D_curr; B_D_curr = t_D;
  }

  /* Boundary rows: i in {1, L-1, L} forced to [1, M]; i=0 forced to [0, 0]. */
  if (L >= 1) { kmin[1] = 1; kmax[1] = M; }
  if (L >= 2) { kmin[L - 1] = 1; kmax[L - 1] = M; }
  if (L >= 1) { kmin[L] = 1; kmax[L] = M; }
  kmin[0] = 0;
  kmax[0] = 0;

  for (i = 1; i <= L; i++)
    ncells += (kmax[i] - kmin[i] + 1);

  {
    const char *dump = getenv("P7IBV_DUMP_BAND");
    if (dump != NULL && *dump != '\0') {
      FILE *fp = fopen(dump, "w");
      if (fp != NULL) {
        fprintf(fp, "# M=%d L=%d delta=%d optimal_milli=%.6f thr_milli=%.6f ncells=%d\n",
                M, L, delta_milli, (double) optimal, (double) thr, ncells);
        fprintf(fp, "# i\tkmin\tkmax\twidth\n");
        for (i = 1; i <= L; i++)
          fprintf(fp, "%d\t%d\t%d\t%d\n", i, kmin[i], kmax[i], kmax[i] - kmin[i] + 1);
        fclose(fp);
      }
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

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
 *   k_stride = ((M+1+15) & ~15)  (16-float = 64-byte alignment)
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


/* ---------------------------------------------------------------------------
 * C1 primitives: per-row SSE forward, backward, through-scan
 * ---------------------------------------------------------------------------*/

/* ibv_forward_one_row -- compute F[i] from F[i-1].
 *
 * emit_row = emit_table[dsq[i]].  FM_prev/FI_prev/FD_prev = row i-1.
 * FM_curr/FI_curr/FD_curr are written (k=0..M).  Tail k>M not touched.
 * Valid for i >= 1.
 */
static void
ibv_forward_one_row(int M, size_t k_stride,
                    const float *MM_t, const float *MI_t, const float *MD_t,
                    const float *IM_t, const float *II_t,
                    const float *DM_t, const float *DD_t,
                    const float *emit_row,
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
 */
static void
ibv_backward_one_row(int M, size_t k_stride, int i, int global_L,
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
                 const float *FM, const float *FI, const float *FD,
                 const float *BM, const float *BI, const float *BD,
                 float *through_scratch,
                 int *ret_kmin, int *ret_kmax)
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

  int row_kmin = -1, row_kmax = -1;
  for (k = 1; k <= M; k++) {
    float t = through_scratch[k];
    if (t < P7IBV_HALF_NEG_INF) continue;
    if (t >= thr) {
      if (row_kmin < 0) row_kmin = k;
      row_kmax = k;
    }
  }
  if (row_kmin < 0) { *ret_kmin = 1; *ret_kmax = M; }
  else              { *ret_kmin = row_kmin; *ret_kmax = row_kmax; }
}


/* ---------------------------------------------------------------------------
 * p7_Seq2BandsIBV -- C1 rewrite (byte-exact vs brief 121 C3)
 * ---------------------------------------------------------------------------*/

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

#define F_M(i)  (FM_pool + (size_t)(i) * k_stride)
#define F_I(i)  (FI_pool + (size_t)(i) * k_stride)
#define F_D(i)  (FD_pool + (size_t)(i) * k_stride)

  /* Row 0: D-cascade init. */
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
    int x  = (int) dsq[i];
    int xt = (x >= 0 && x < K) ? x : K;
    ibv_forward_one_row(M, k_stride,
                        MM_t, MI_t, MD_t, IM_t, II_t, DM_t, DD_t,
                        emit_table[xt],
                        F_M(i-1), F_I(i-1), F_D(i-1),
                        F_M(i),   F_I(i),   F_D(i));
  }

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

  ESL_ALLOC(i2k,  sizeof(int) * (L + 1));
  ESL_ALLOC(kmin, sizeof(int) * (L + 1));
  ESL_ALLOC(kmax, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, L + 1, -1);
  for (i = 0; i <= L; i++) { kmin[i] = 1; kmax[i] = M; }

  float *B_M_prev = BM_a, *B_I_prev = BI_a, *B_D_prev = BD_a;
  float *B_M_curr = BM_b, *B_I_curr = BI_b, *B_D_curr = BD_b;

  /* Row L: terminal injection. */
  ibv_backward_one_row(M, k_stride, L, L,
                       MM_t, MI_t, MD_t, IM_t, II_t, DM_t, DD_t,
                       NULL, NULL, NULL,
                       B_M_prev, B_I_prev, B_D_prev);

  for (i = L - 1; i >= 0; i--) {
    int x_next = (int) dsq[i + 1];
    int xt     = (x_next >= 0 && x_next < K) ? x_next : K;

    ibv_backward_one_row(M, k_stride, i, L,
                         MM_t, MI_t, MD_t, IM_t, II_t, DM_t, DD_t,
                         emit_table[xt], B_M_prev, B_I_prev,
                         B_M_curr, B_I_curr, B_D_curr);

    if (i >= 2 && i <= L - 2)
      ibv_through_scan(M, k_stride, thr,
                       F_M(i), F_I(i), F_D(i),
                       B_M_curr, B_I_curr, B_D_curr,
                       through, &kmin[i], &kmax[i]);

    float *t_M = B_M_prev; B_M_prev = B_M_curr; B_M_curr = t_M;
    float *t_I = B_I_prev; B_I_prev = B_I_curr; B_I_curr = t_I;
    float *t_D = B_D_prev; B_D_prev = B_D_curr; B_D_curr = t_D;
  }

  if (L >= 1) { kmin[1] = 1; kmax[1] = M; }
  if (L >= 2) { kmin[L - 1] = 1; kmax[L - 1] = M; }
  if (L >= 1) { kmin[L] = 1; kmax[L] = M; }
  kmin[0] = 0; kmax[0] = 0;

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

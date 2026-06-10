/* p7_ibv.c -- F+B direct-band band derivation for cmalign --p7ibv
 *
 * Brief 121 C2: SSE M+I vectorized along k, scalar D-fill (in-row dep).
 *   - State-major padded-k pools from C1 unchanged.
 *   - Forward: per row i>=1, scalar prefix k=0..3 for M/I/D; SSE bulk
 *     for M and I from k=4 while k+3<=M; scalar tail for M and I;
 *     scalar left-to-right D-fill for k=4..M (k=0..3 already done in
 *     prefix). Row i=0 (no residue) stays scalar (deletion cascade only).
 *   - Backward: per row i<L, scalar right-to-left D-fill k=M..0; SSE
 *     bulk for M and I from k=0 while k+3<=M; scalar tail for M and I.
 *     Row i=L stays scalar (seed at (L,M) + deletion cascade on the
 *     terminal row).
 *   - Through-score scan stays scalar (cheap, 1 pass per row); C3 will
 *     integrate it into the backward sweep when streaming B.
 *
 * D-state in-row dep: D[i,k] depends on D[i,k-1] (forward) or
 * D[i,k+1] (backward) on the SAME row. NOT SIMD-able with _mm_max_ps
 * along k. Per brief 120 §8.3 and brief 121 §4.2, accept scalar D-fill;
 * D is a small fraction of per-row work.
 *
 * Padded tails of every row (k > M) MUST stay -INF, so future row's
 * SSE loads of k+1..k+3 at the right edge pick up -INF instead of stale
 * data. SSE writes are kept inside [0, M].
 *
 * ULP analysis (float, INTSCALE=1000 milli-bits):
 *   Cumulative path scores in practice are O(1e5) milli-bits for
 *   biologically meaningful matches; brief 120 estimated ~3e5 in
 *   magnitude. Float ULP at magnitude 3e5 is 3e5 / 2^23 ~= 0.04
 *   milli-bit. EPS = 1 milli-bit is 25x above that, so the
 *   "==optimal" through-score equality test absorbs the rounding cleanly.
 *
 *   Worst case (very long M+L with cumulative magnitude approaching 1e7):
 *   float ULP rises to ~0.6 milli-bit. EPS still 1.5x above; tight but
 *   safe. We assert optimal < 1e7 milli-bits as a sanity check.
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
#define P7IBV_K_ALIGN         16              /* k-stride alignment in floats (= 64 bytes) */
#define P7IBV_OPTIMAL_SANITY  1.0e7f

/* Aligned float-buffer allocator. 64-byte alignment = 1 cache line = 4 SSE registers. */
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

/* INTSCALE * log2(P) in milli-bits. */
static inline float
p7ibv_lod_milli(float p)
{
  if (p <= 0.0f) return P7IBV_NEG_INF;
  return (float)(P7IBV_INTSCALE * (log((double) p) / M_LN2));
}

/* Match-state emission log-odds in milli-bits, using insert emit as background. */
static inline float
p7ibv_emit_milli(const P7_HMM *hmm, int k, int x)
{
  if (x < 0 || x >= hmm->abc->K) return 0.0f;
  float pm = hmm->mat[k][x];
  float pi = hmm->ins[k][x];
  if (pm <= 0.0f || pi <= 0.0f) return P7IBV_NEG_INF;
  return (float)(P7IBV_INTSCALE * (log((double) pm / (double) pi) / M_LN2));
}

/* 4-way SSE max. */
static inline __m128
p7ibv_mm_max3(__m128 a, __m128 b, __m128 c)
{
  return _mm_max_ps(_mm_max_ps(a, b), c);
}


/* Function: p7_Seq2BandsIBV()
 *
 * Synopsis: Derive p7 bands from full F+B + Delta threshold.
 *
 * Args:     cm           - CM with cm->fp7 set
 *           errbuf       - for error messages
 *           dsq          - digital sequence, 1..L
 *           L            - length of dsq
 *           delta_milli  - Delta threshold in milli-bits (3000 = 3 bits, default)
 *           ret_i2k      - RETURN: per-residue pin array (all -1; caller frees)
 *           ret_kmin     - RETURN: per-residue kmin array, length L+1 (caller frees)
 *           ret_kmax     - RETURN: per-residue kmax array, length L+1 (caller frees)
 *           ret_ncells   - RETURN: total banded cells sum_{i=1..L} (kmax[i]-kmin[i]+1)
 *
 * Returns:  eslOK on success; eslEMEM/eslEINVAL on error.
 */
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
  float    *BM_pool = NULL, *BI_pool = NULL, *BD_pool = NULL;
  float    *emit_pool = NULL;
  float   **emit_table = NULL;          /* [K+1][k_stride], last row all-zeros (ambig) */
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

  /* Precompute emit log-odds per (residue, k). Last row (index K) is
   * all-zeros for ambiguous/unknown residues. Pads beyond M+1 to -INF
   * so SSE loads in the M-output at k=M-3..M (which look at emit[k+1])
   * fold harmlessly into the boundary -INF.
   */
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
  if ((status = ibv_alloc_floats(pool_cells, &BM_pool)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(pool_cells, &BI_pool)) != eslOK) goto ERROR;
  if ((status = ibv_alloc_floats(pool_cells, &BD_pool)) != eslOK) goto ERROR;
  for (size_t c = 0; c < pool_cells; c++) {
    FM_pool[c] = FI_pool[c] = FD_pool[c] = P7IBV_NEG_INF;
    BM_pool[c] = BI_pool[c] = BD_pool[c] = P7IBV_NEG_INF;
  }

#define F_M(i)  (FM_pool + (size_t)(i) * k_stride)
#define F_I(i)  (FI_pool + (size_t)(i) * k_stride)
#define F_D(i)  (FD_pool + (size_t)(i) * k_stride)
#define B_M(i)  (BM_pool + (size_t)(i) * k_stride)
#define B_I(i)  (BI_pool + (size_t)(i) * k_stride)
#define B_D(i)  (BD_pool + (size_t)(i) * k_stride)

  /* ---------- Forward ---------- */
  F_M(0)[0] = 0.0f;
  /* Row i=0: no residue. F_M[0,0]=0 set; F_M[0,k>=1] = -INF; F_I[0,*] = -INF;
   * F_D[0, k>=1] = max(F_M[0, k-1] + MD_t[k-1], F_D[0, k-1] + DD_t[k-1])
   * (deletion cascade through the model). Scalar fill; only this row needs it.
   */
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

    /* Scalar prefix: k=0..3 (handles k=0 boundary, supplies fm_i[3] etc. for
     * the SSE block's k=4 load of fm_im1[3] -- not needed but D-fill at k=4
     * needs fd_i[3] which is set here).
     */
    int kpref_end = (M < 3) ? M : 3;
    for (k = 0; k <= kpref_end; k++) {
      float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF, cD = P7IBV_NEG_INF;
      /* I: i>=1 always holds in this loop. */
      {
        float a = fm_im1[k] + MI_t[k];
        float b = fi_im1[k] + II_t[k];
        cI = (a > b) ? a : b;
      }
      if (k >= 1) {
        /* M */
        float a = fm_im1[k - 1] + MM_t[k - 1];
        float b = fi_im1[k - 1] + IM_t[k - 1];
        float c = fd_im1[k - 1] + DM_t[k - 1];
        float m = (a > b) ? a : b;
        if (c > m) m = c;
        cM = m + emit_row[k];
        /* D (in-row) */
        float a2 = fm_i[k - 1] + MD_t[k - 1];
        float b2 = fd_i[k - 1] + DD_t[k - 1];
        cD = (a2 > b2) ? a2 : b2;
      }
      fm_i[k] = cM;
      fi_i[k] = cI;
      fd_i[k] = cD;
    }

    /* SSE bulk for M and I: k=4..k_sse_end-1 in chunks of 4 where k+3 <= M. */
    int k_sse = (M < 3) ? 4 : 4;
    int k_sse_end = k_sse;
    while (k_sse_end + 3 <= M) k_sse_end += 4;
    for (k = k_sse; k < k_sse_end; k += 4) {
      /* M_out: max(F_M[i-1,k-1..k+2]+MM[k-1..k+2], +IM, +DM) + emit[k..k+3] */
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

      /* I_out: max(F_M[i-1,k..k+3]+MI[k..k+3], F_I[i-1,k..k+3]+II[k..k+3]) */
      __m128 fm_k = _mm_loadu_ps(&fm_im1[k]);
      __m128 fi_k = _mm_loadu_ps(&fi_im1[k]);
      __m128 t_mi = _mm_loadu_ps(&MI_t[k]);
      __m128 t_ii = _mm_loadu_ps(&II_t[k]);
      __m128 a2   = _mm_add_ps(fm_k, t_mi);
      __m128 b2   = _mm_add_ps(fi_k, t_ii);
      _mm_storeu_ps(&fi_i[k], _mm_max_ps(a2, b2));
    }

    /* Scalar tail for M and I from k=k_sse_end to k=M. */
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

    /* D scalar fill, left-to-right, k=4..M (k=0..3 done in prefix). */
    for (k = (kpref_end + 1 > M) ? (M + 1) : (kpref_end + 1); k <= M; k++) {
      float a = fm_i[k - 1] + MD_t[k - 1];
      float b = fd_i[k - 1] + DD_t[k - 1];
      fd_i[k] = (a > b) ? a : b;
    }
  }

  /* ---------- Backward ---------- */
  /* Row i=L (terminal): seed (L,M)=0, plus in-row D-deletion cascade for
   * B_D[L, k<M] = DD_t[k]+B_D[L,k+1] (starting from B_D[L,M]=0). And in-row
   * to B_M via MD: B_M[L, k<M] from MD_t[k]+B_D[L,k+1]. B_I[L,*]=-INF for k<M.
   * Scalar; one row only.
   */
  {
    float *bm_L = B_M(L);
    float *bi_L = B_I(L);
    float *bd_L = B_D(L);
    bm_L[M] = 0.0f;
    bi_L[M] = 0.0f;
    bd_L[M] = 0.0f;
    for (k = M - 1; k >= 0; k--) {
      float bv = bd_L[k + 1];
      float a = MD_t[k] + bv;
      float b = DD_t[k] + bv;
      bm_L[k] = a;
      bd_L[k] = b;
      bi_L[k] = P7IBV_NEG_INF;
    }
  }

  for (i = L - 1; i >= 0; i--) {
    int   x_next  = (int) dsq[i + 1];   /* i+1<=L always here */
    int   xt      = (x_next >= 0 && x_next < K) ? x_next : K;
    float *emit_row_next = emit_table[xt];
    float *bm_i   = B_M(i);
    float *bi_i   = B_I(i);
    float *bd_i   = B_D(i);
    float *bm_ip1 = B_M(i + 1);
    float *bi_ip1 = B_I(i + 1);

    /* D scalar fill, right-to-left, k=M..0.
     * cD = max(DM_t[k] + bm_ip1[k+1] + emit_row_next[k+1],   [if k+1<=M]
     *          DD_t[k] + bd_i[k+1])                          [if k+1<=M]
     * else -INF.
     */
    bd_i[M] = P7IBV_NEG_INF;  /* k=M has neither contribution */
    for (k = M - 1; k >= 0; k--) {
      float bv_m = bm_ip1[k + 1] + emit_row_next[k + 1];
      float bv_d = bd_i[k + 1];
      float a = DM_t[k] + bv_m;
      float b = DD_t[k] + bv_d;
      bd_i[k] = (a > b) ? a : b;
    }

    /* SSE bulk M+I: k=0..k_sse_end-1 in chunks of 4 where k+3 <= M.
     * For each k-block, we need:
     *   bv_M_vec = bm_ip1[k+1..k+4] + emit_row_next[k+1..k+4]
     *   bv_I_vec = bi_ip1[k..k+3]
     *   bv_D_vec = bd_i[k+1..k+4]   (just filled by D scalar pass)
     */
    int k_sse_end = 0;
    while (k_sse_end + 3 <= M) k_sse_end += 4;
    for (k = 0; k < k_sse_end; k += 4) {
      __m128 bm_n     = _mm_loadu_ps(&bm_ip1[k + 1]);
      __m128 e_n      = _mm_loadu_ps(&emit_row_next[k + 1]);
      __m128 bv_M     = _mm_add_ps(bm_n, e_n);
      __m128 bv_I     = _mm_loadu_ps(&bi_ip1[k]);
      __m128 bv_D     = _mm_loadu_ps(&bd_i[k + 1]);

      __m128 t_mm     = _mm_loadu_ps(&MM_t[k]);
      __m128 t_mi     = _mm_loadu_ps(&MI_t[k]);
      __m128 t_md     = _mm_loadu_ps(&MD_t[k]);
      __m128 a_M      = _mm_add_ps(t_mm, bv_M);
      __m128 b_M      = _mm_add_ps(t_mi, bv_I);
      __m128 c_M      = _mm_add_ps(t_md, bv_D);
      __m128 cM_out   = p7ibv_mm_max3(a_M, b_M, c_M);
      _mm_storeu_ps(&bm_i[k], cM_out);

      __m128 t_im     = _mm_loadu_ps(&IM_t[k]);
      __m128 t_ii     = _mm_loadu_ps(&II_t[k]);
      __m128 a_I      = _mm_add_ps(t_im, bv_M);
      __m128 b_I      = _mm_add_ps(t_ii, bv_I);
      __m128 cI_out   = _mm_max_ps(a_I, b_I);
      _mm_storeu_ps(&bi_i[k], cI_out);
    }

    /* Scalar tail for M and I from k=k_sse_end..M. */
    for (k = k_sse_end; k <= M; k++) {
      float cM = P7IBV_NEG_INF, cI = P7IBV_NEG_INF;
      if (k + 1 <= M) {
        float bv_m = bm_ip1[k + 1] + emit_row_next[k + 1];
        float a = MM_t[k] + bv_m;
        float b = IM_t[k] + bv_m;
        if (a > cM) cM = a;
        if (b > cI) cI = b;
      }
      /* I-block always contributes (i+1<=L, k always in range for bi_ip1[k]) */
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
  }

  /* ---------- optimal score + threshold ---------- */
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

  /* ---------- per-row band emission ---------- */
  for (i = 1; i <= L; i++) {
    int   row_kmin = -1, row_kmax = -1;
    float *fm = F_M(i), *fi = F_I(i), *fd = F_D(i);
    float *bm = B_M(i), *bi = B_I(i), *bd = B_D(i);
    for (k = 1; k <= M; k++) {
      float t_m = fm[k] + bm[k];
      float t_i = fi[k] + bi[k];
      float t_d = fd[k] + bd[k];
      float th  = t_m;
      if (t_i > th) th = t_i;
      if (t_d > th) th = t_d;
      if (th < P7IBV_HALF_NEG_INF) continue;
      if (th >= thr) {
        if (row_kmin < 0) row_kmin = k;
        row_kmax = k;
      }
    }
    if (row_kmin < 0) {
      kmin[i] = 1;
      kmax[i] = M;
    } else {
      kmin[i] = row_kmin;
      kmax[i] = row_kmax;
    }
  }

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
#undef B_M
#undef B_I
#undef B_D

  free(MM_t); free(MI_t); free(MD_t);
  free(IM_t); free(II_t); free(DM_t); free(DD_t);
  free(FM_pool); free(FI_pool); free(FD_pool);
  free(BM_pool); free(BI_pool); free(BD_pool);
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
  if (BM_pool) free(BM_pool); if (BI_pool) free(BI_pool); if (BD_pool) free(BD_pool);
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

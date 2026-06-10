/* p7_ibv.c -- F+B direct-band band derivation for cmalign --p7ibv
 *
 * Brief 121 layout: float (not double) cells, state-major padded-k
 * pools sized to a 16-float (64-byte) stride. Preserves brief 120's
 * algorithm exactly; only the score type and per-row stride change.
 *
 *   - F_M[i][k], F_I[i][k], F_D[i][k] : (L+1) rows x k_stride floats
 *   - B_M[i][k], B_I[i][k], B_D[i][k] : same
 *   - k_stride = ((M+1 + 15) & ~15)   (16-float = 64-byte = 1 cache line)
 *   - Pools allocated 64-byte aligned via posix_memalign.
 *   - Padded tail cells (k in [M+1, k_stride-1]) initialized to -INF, so
 *     future SSE max-reductions (C2) won't pick them up.
 *
 * ULP analysis (float, INTSCALE=1000 milli-bits):
 *   Cumulative path scores in practice are O(1e5) milli-bits for
 *   biologically meaningful matches; the brief-120 comment estimated
 *   ~3e5 in magnitude. Float ULP at magnitude 3e5 is 3e5 / 2^23 ~= 0.04
 *   milli-bit. EPS = 1 milli-bit is 25x above that, so the
 *   "==optimal" through-score equality test absorbs the rounding cleanly.
 *
 *   Worst case (very long M+L with cumulative magnitude approaching 1e7):
 *   float ULP rises to ~0.6 milli-bit. EPS still 1.5x above; tight but
 *   safe. We assert optimal < 1e7 milli-bits as a sanity check; LSU
 *   M=3400 L=2771 gives ~1e5-3e5, well below the assert.
 *
 * Algorithm (per brief 116-117), unchanged from brief 120:
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
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

#define P7IBV_INTSCALE     1000.0f
#define P7IBV_NEG_INF      (-1.0e18f)
#define P7IBV_HALF_NEG_INF (-5.0e17f)
#define P7IBV_EPS          1.0f
#define P7IBV_K_ALIGN      16              /* k-stride alignment in floats (= 64 bytes) */
#define P7IBV_OPTIMAL_SANITY  1.0e7f

/* Aligned float-buffer allocator. Returns eslOK and sets *ret_p, or eslEMEM.
 * 64-byte alignment = 1 cache line = 4 SSE registers.
 */
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
 *           Memory cost (brief 121 C1): 6 * (L+1) * k_stride * 4 bytes
 *           (state-major padded-k float pools). For LSU M=3400 L=2771,
 *           k_stride=3408 -> ~226 MB. Brief 121 C3 will halve B to a
 *           rolling 2-row buffer (-> ~113 MB total).
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
  float    *MM_t = NULL, *MI_t = NULL, *MD_t = NULL;
  float    *IM_t = NULL, *II_t = NULL;
  float    *DM_t = NULL, *DD_t = NULL;
  float    *FM_pool = NULL, *FI_pool = NULL, *FD_pool = NULL;
  float    *BM_pool = NULL, *BI_pool = NULL, *BD_pool = NULL;
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
  if (L < 1 || M < 1)
    ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsIBV: bad L=%d or M=%d", L, M);

  /* State-major padded-k layout. k_stride padded up to multiple of
   * P7IBV_K_ALIGN floats so each row's start is 64-byte aligned and
   * the row's SSE-padded tail (k in [M+1, k_stride-1]) is a full set
   * of 4-float SSE registers that can be filled with -INF and ignored
   * by future SSE max-reductions.
   */
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
  /* Pad transition tails to -INF so SSE loads off the end can't pick up
   * a "good" transition that wraps the row boundary.
   */
  for (k = M + 1; k < (int) k_stride; k++) {
    MM_t[k] = MI_t[k] = MD_t[k] = P7IBV_NEG_INF;
    IM_t[k] = II_t[k] = P7IBV_NEG_INF;
    DM_t[k] = DD_t[k] = P7IBV_NEG_INF;
  }

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

  /* Forward recurrence. Glocal seed: F_M[0,0] = 0; rest -INF. */
  F_M(0)[0] = 0.0f;
  for (i = 0; i <= L; i++) {
    x = (i >= 1) ? (int) dsq[i] : -1;
    float *fm_i   = F_M(i);
    float *fi_i   = F_I(i);
    float *fd_i   = F_D(i);
    float *fm_im1 = (i >= 1) ? F_M(i - 1) : NULL;
    float *fi_im1 = (i >= 1) ? F_I(i - 1) : NULL;
    float *fd_im1 = (i >= 1) ? F_D(i - 1) : NULL;
    for (k = 0; k <= M; k++) {
      if (i == 0 && k == 0) continue;

      float cM = P7IBV_NEG_INF;
      float cI = P7IBV_NEG_INF;
      float cD = P7IBV_NEG_INF;

      if (i >= 1) {
        float a = fm_im1[k] + MI_t[k];
        float b = fi_im1[k] + II_t[k];
        cI = (a > b) ? a : b;
      }

      if (k >= 1) {
        float a = fm_i[k - 1] + MD_t[k - 1];
        float b = fd_i[k - 1] + DD_t[k - 1];
        cD = (a > b) ? a : b;
      }

      if (i >= 1 && k >= 1) {
        float e = p7ibv_emit_milli(hmm, k, x);
        float a = fm_im1[k - 1] + MM_t[k - 1];
        float b = fi_im1[k - 1] + IM_t[k - 1];
        float c = fd_im1[k - 1] + DM_t[k - 1];
        float m = (a > b) ? a : b;
        if (c > m) m = c;
        cM = m + e;
      }

      fm_i[k] = cM;
      fi_i[k] = cI;
      fd_i[k] = cD;
    }
  }

  /* Backward recurrence. Glocal seed at (L, M): B_M=B_I=B_D=0; rest -INF. */
  for (i = L; i >= 0; i--) {
    int x_next = (i + 1 <= L) ? (int) dsq[i + 1] : -1;
    float *bm_i   = B_M(i);
    float *bi_i   = B_I(i);
    float *bd_i   = B_D(i);
    float *bm_ip1 = (i + 1 <= L) ? B_M(i + 1) : NULL;
    float *bi_ip1 = (i + 1 <= L) ? B_I(i + 1) : NULL;
    for (k = M; k >= 0; k--) {
      float term = (i == L && k == M) ? 0.0f : P7IBV_NEG_INF;
      float cM = term;
      float cI = term;
      float cD = term;

      if (i + 1 <= L && k + 1 <= M) {
        float e = p7ibv_emit_milli(hmm, k + 1, x_next);
        float bv = bm_ip1[k + 1] + e;
        float a = MM_t[k] + bv;
        float b = IM_t[k] + bv;
        float c = DM_t[k] + bv;
        if (a > cM) cM = a;
        if (b > cI) cI = b;
        if (c > cD) cD = c;
      }

      if (i + 1 <= L) {
        float bv = bi_ip1[k];
        float a = MI_t[k] + bv;
        float b = II_t[k] + bv;
        if (a > cM) cM = a;
        if (b > cI) cI = b;
      }

      if (k + 1 <= M) {
        float bv = bd_i[k + 1];
        float a = MD_t[k] + bv;
        float b = DD_t[k] + bv;
        if (a > cM) cM = a;
        if (b > cD) cD = b;
      }

      bm_i[k] = cM;
      bi_i[k] = cI;
      bd_i[k] = cD;
    }
  }

  /* optimal = max{F_M[L,M], F_I[L,M], F_D[L,M]}. */
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
  if (i2k)  free(i2k);
  if (kmin) free(kmin);
  if (kmax) free(kmax);
  *ret_i2k    = NULL;
  *ret_kmin   = NULL;
  *ret_kmax   = NULL;
  *ret_ncells = 0;
  return status;
}

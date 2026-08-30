/* cm_p7_band_chk.c
 *
 * Checkpointed banded CP9 P7B Forward/Backward + band reduction, double
 * precision (brief 26_0430-150/153/154). Produces cp9b bands with
 * O(sqrt(L)*avg_bw) working memory instead of full ncells-sized banded CP9
 * P7B matrices.
 *
 * History: this file originally also held a parallel INTEGER-precision
 * checkpointed F/B family (brief 26_0430-146, the original 144-B path). That
 * int family's roundoff caused genome-scale posterior collapse (brief
 * 26_0430-154), it was superseded end-to-end by the double-precision family
 * below, it had acquired zero callers, and it was deleted in brief
 * 26_0430-306. To retrieve it: `git log -S cp9_chk_bwd_row -- src/cm_p7_band_chk.c`
 * locates the deletion commit; read that commit message in full before
 * reviving anything from it — the int family never received the M_k <- EL_k
 * banded-Backward EL fold that every live kernel in this file has, so
 * reviving it as-is reintroduces an inflated-posterior bug.
 *
 * Design (see project memory project_ckpt_dual_sweep):
 *   Bands are derived via TWO opposite-direction cumulative dlogsum sweeps
 *   (MIN ascending i; MAX descending i). A single HMMER-style descending
 *   decode does not map. So:
 *     1. Forward-checkpointed fill: store F planes at block boundaries.
 *     2. Backward-checkpointed fill: store B planes at boundaries; also yields
 *        sc = bmx->mmx[0][0] (the posterior denominator).
 *     3. MIN sweep: segments left->right; per segment recompute fwd & bck from
 *        the boundary checkpoints, build posterior rows, accumulate min.
 *     4. MAX sweep: segments right->left; accumulate max.
 *   The float/truncated family below (brief 26_0430-150, promoted to double
 *   precision by brief 154) is a verbatim transcription of the
 *   non-checkpointed float originals cp9_ForwardP7BF / cp9_BackwardP7BF /
 *   cp9_PredictStartAndEndPositionsP7BF (cm_p7_band.c), so checkpointed bands
 *   match that non-checkpointed path modulo the double promotion. The
 *   originals are kept untouched as the structural reference.
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

/* macros local to cm_p7_band.c; re-defined here. INBAND needs locals named
 * kmin/kmax; CP9TSC needs a local named tsc. */
#define CP9TSC(s,k) (tsc[(k) * cp9O_NTRANS + (s)])
#define INBAND(i,k) ((k >= kmin[i]) && (k <= kmax[i]))

/*****************************************************************
 * FLOAT/TRUNCATED PATH (brief 26_0430-150, Phase 2).
 *
 * A parallel set of float-typed checkpointed kernels mirroring the int
 * (non-trunc) machinery above, but transcribed VERBATIM from the float
 * originals cp9_ForwardP7BF / cp9_BackwardP7BF / cp9_PosteriorP7BF (inline in
 * cp9_FB2HMMBandsP7BF) / cp9_PredictStartAndEndPositionsP7BF (cm_p7_band.c),
 * so that the checkpointed float bands are BYTE-IDENTICAL to the
 * non-checkpointed float truncated path. Two float-specific structural
 * features that the int path lacks and MUST be reproduced:
 *   (a) cp9_ForwardP7BF has the brief-134 BM-coverage supplementary pass
 *       (begin + EL-from-into-M for the full row band minus the match range);
 *       it writes mmx/elmx cells the posterior reads.
 *   (b) the erow EL-from-M+1 accumulation is NOT gated on INBAND(i,M).
 *       (erow is dead for band derivation, but reproduced for fidelity.)
 *
 * Function/struct suffix F. Substitutions vs the originals:
 *   mmx[i]->mc, mmx[i-1]->mp, mmx[i+1]->mn (and i/d/el likewise).
 *
 * BRIEF 154: this FLOAT-suffix family is now DOUBLE precision (cells +
 * accumulation + checkpoint storage end-to-end). The non-checkpointed
 * cp9_ForwardP7BF/cp9_BackwardP7BF (cm_p7_band.c) stay float; the structural
 * 1:1 correspondence with those float originals is preserved so they can be
 * diffed (only the cell type and the logsum changed). Brief 26_0430-153 proved the
 * genome-scale -g truncated band collapse is float32 accumulation roundoff in
 * this CP9 F/B (~-35 nats at HSV, ~-45 at MPXV); double drives the F/B gap to
 * ~0. p7_FLogsum is a float LUT returning float, so it cannot carry double
 * precision; cp9_chk_dlogsum (exact log1p/exp, double-valued) replaces it.
 *****************************************************************/

/* Brief 26_0430-193: double-precision LUT logsum, the production replacement for
 * the exact log1p/exp transcendental in the checkpointed-truncated CP9 F/B.
 * Mirrors p7_FLogsum's float LUT (hmmer/src/logsum.c) but double-valued with a
 * widened range: double eps (2.2e-16) pushes the "just return max" cutoff from
 * float's 15.7 nats to CP9_DLOGSUM_CUTOFF (36.7 nats). Rests on the brief
 * 26_0430-153/154 finding that logsum *discretization* error was never the source
 * of the genome-scale float32 collapse (float32 *storage* accumulation was); the
 * matrix storage stays double (brief 154), only the per-call transcendental is
 * replaced by an ~5-10 cy table lookup. Validated in brief 26_0430-193 Task 1:
 * byte-identical alignments to the exact path at genome scale (HSV L=146678, MPXV
 * L=197226) and dengue/SARS kmerchain; the only divergence is within-noise low-PP
 * near-tie tipping (norovirus kmerchain, ~4/30 seqs, alignment structure preserved).
 * Task 2: 2.3-2.65x on stage c (cp9_IterateSeq2BandsP7B), 2.2-2.3x end-to-end for
 * compute-bound kmerchain dengue/SARS. Chosen config = brief's candidate A
 * (SCALE=1000 → 0.001-nat buckets, CUTOFF=36.7 → ~296 KB L2-resident table; the
 * finer 10x table was measurably slower with no accuracy benefit at genome scale).
 *
 * CP9_DLOGSUM_EXACT=1 restores the exact log1p/exp transcendental (debug/regression
 * escape hatch). Table built once via a load-time constructor (single-threaded,
 * before cmalign spawns workers) so the hot path reads a read-only table. */
#define CP9_DLOGSUM_SCALE  1000.0
/* BASE 2. Cutoff is where 2^-x drops under double epsilon: 2^-53 = 1.1e-16.
 * (The former base-e cutoff 36.7 was the same threshold expressed in nats,
 * e^-36.7 = 1.1e-16 -- so this is the identical precision bound, restated in
 * the base these scores are actually on.) */
#define CP9_DLOGSUM_CUTOFF 53.0
#define CP9_DLOGSUM_TBLN   53002    /* (int)(CUTOFF*SCALE)+2 = 53000+2 */
static double  cp9_dlogsum_tbl[CP9_DLOGSUM_TBLN];
static int     cp9_dlogsum_exact = 0;   /* set by CP9_DLOGSUM_EXACT=1 */

__attribute__((constructor)) static void
cp9_chk_dlogsum_lut_init(void)
{
  char *s;
  int i;
  /* BASE 2, matching ILogsum's sreLOG2(1.+sreEXP2(-i/INTSCALE)) exactly. */
  for(i = 0; i < CP9_DLOGSUM_TBLN; i++) cp9_dlogsum_tbl[i] = log2(1. + exp2((double) -i / CP9_DLOGSUM_SCALE));
  if((s = getenv("CP9_DLOGSUM_EXACT")) != NULL && atoi(s) != 0) cp9_dlogsum_exact = 1;
}

/* Double-precision log-sum for the checkpointed double-trunc CP9 F/B kernels
 * (brief 26_0430-153/154/193). -inf-guarded. LUT by default; exact transcendental
 * under CP9_DLOGSUM_EXACT=1.
 *
 * BASE 2 (briefs 26_0821-024/025/026). Every value this combines is a
 * Scorify()d score, i.e. BITS -- Infernal has used bits, not nats, since 2007
 * (see the EPN note in logsum.c). This was transcribed from the int kernels
 * with ILogsum -> p7_FLogsum, which silently changed the base while the values
 * stayed in bits; it under-added every sum by up to 1-ln2 = 0.307 bits,
 * one-directionally, compounding row by row. */
static inline double
cp9_chk_dlogsum(double a, double b)
{
  if(a == -eslINFINITY) return b;
  if(b == -eslINFINITY) return a;
  if(cp9_dlogsum_exact) {
    if(a > b) return a + log2(1. + exp2(b - a));
    else      return b + log2(1. + exp2(a - b));
  }
  const double max = (a > b) ? a : b;
  const double min = (a > b) ? b : a;
  return ((max - min) >= CP9_DLOGSUM_CUTOFF) ? max
         : max + cp9_dlogsum_tbl[(int)((max - min) * CP9_DLOGSUM_SCALE)];
}

/* CP9_DMX: file-local double mirror of CP9_FMX, used ONLY by the env-gated
 * CP9_CKPTF_FBCMP debug block below (the chk kernels write double planes, so the
 * scratch matrix that block hands them must be double). Full (unbanded) layout,
 * matching CreateCP9FMatrix(L,M)'s rows-at-i*(M+1) scratch use. */
typedef struct {
  double **mmx, **imx, **dmx, **elmx;
  double  *mmx_mem, *imx_mem, *dmx_mem, *elmx_mem;
  int      M, rows;
} CP9_DMX;

static CP9_DMX *
CreateCP9DMatrix(int N, int M)
{
  int status;
  CP9_DMX *mx;
  int i;
  ESL_ALLOC(mx,      sizeof(CP9_DMX));
  ESL_ALLOC(mx->mmx, sizeof(double *) * (N+1));
  ESL_ALLOC(mx->imx, sizeof(double *) * (N+1));
  ESL_ALLOC(mx->dmx, sizeof(double *) * (N+1));
  ESL_ALLOC(mx->elmx,sizeof(double *) * (N+1));
  ESL_ALLOC(mx->mmx_mem, sizeof(double) * ((N+1)*(M+1)));
  ESL_ALLOC(mx->imx_mem, sizeof(double) * ((N+1)*(M+1)));
  ESL_ALLOC(mx->dmx_mem, sizeof(double) * ((N+1)*(M+1)));
  ESL_ALLOC(mx->elmx_mem,sizeof(double) * ((N+1)*(M+1)));
  memset(mx->mmx_mem,  0, sizeof(double) * ((N+1)*(M+1)));
  memset(mx->imx_mem,  0, sizeof(double) * ((N+1)*(M+1)));
  memset(mx->dmx_mem,  0, sizeof(double) * ((N+1)*(M+1)));
  memset(mx->elmx_mem, 0, sizeof(double) * ((N+1)*(M+1)));
  mx->mmx[0] = mx->mmx_mem; mx->imx[0] = mx->imx_mem;
  mx->dmx[0] = mx->dmx_mem; mx->elmx[0]= mx->elmx_mem;
  for (i = 1; i <= N; i++) {
    mx->mmx[i] = mx->mmx[0] + (i*(M+1));
    mx->imx[i] = mx->imx[0] + (i*(M+1));
    mx->dmx[i] = mx->dmx[0] + (i*(M+1));
    mx->elmx[i]= mx->elmx[0]+ (i*(M+1));
  }
  mx->M = M; mx->rows = N;
  return mx;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL;
}

static void
FreeCP9DMatrix(CP9_DMX *mx)
{
  if(mx == NULL) return;
  free(mx->mmx_mem); free(mx->imx_mem); free(mx->dmx_mem); free(mx->elmx_mem);
  free(mx->mmx); free(mx->imx); free(mx->dmx); free(mx->elmx);
  free(mx);
}

/* Forward row 0 (init). Verbatim from cp9_ForwardP7BF i=0 block.
 * Planes mc/ic/dc/ec are row-0 (offset 0 = kmin[0]==0). */
static void
cp9_chk_fwd_row0F(CP9_t *cp9, int *kmin, int *kmax, int M,
                  double *mc, double *ic, double *dc, double *ec, double *ret_erow0)
{
  int const *tsc = cp9->otsc;
  int k, kn, kp;
  double sc, erow0;

  mc[0] = 0.;       /* M_0 is state B */
  ic[0] = -eslINFINITY;
  dc[0] = -eslINFINITY;
  ec[0] = -eslINFINITY;

  kn = ESL_MAX(1, kmin[0]);
  kp = kn - kmin[0];
  for (k = kn; k <= kmax[0]; k++, kp++) {
    mc[kp] = ic[kp] = ec[kp] = -eslINFINITY;
    sc = -eslINFINITY;
    if(kp > 0) {
      sc = cp9_chk_dlogsum(cp9_chk_dlogsum(mc[kp-1] + Scorify(CP9TSC(cp9O_MD,k-1)),
                                 ic[kp-1] + Scorify(CP9TSC(cp9O_ID,k-1))),
                      dc[kp-1] + Scorify(CP9TSC(cp9O_DD,k-1)));
    }
    dc[kp] = sc;
  }
  erow0 = -eslINFINITY;
  if(INBAND(0, M)) { erow0 = dc[M] + Scorify(CP9TSC(cp9O_DM,M)); } /* kmin[0]==0 so abs M == rel M */
  if(ret_erow0) *ret_erow0 = erow0;
}

/* Forward row i (1..L). Verbatim from cp9_ForwardP7BF main loop body, INCLUDING
 * the brief-134 BM-coverage supplementary pass. Cur planes mc/ic/dc/ec (offset
 * 0 = kmin[i]); prev planes mp/ip/dp/ep (offset 0 = kmin[i-1]). */
static void
cp9_chk_fwd_rowF(CP9_t *cp9, ESL_DSQ *dsq, int i, int *kmin, int *kmax, int M,
                 double *mp, double *ip, double *dp, double *ep,
                 double *mc, double *ic, double *dc, double *ec, double *ret_erow)
{
  int const *tsc = cp9->otsc;
  int const *isc_i = cp9->isc[dsq[i]];
  int const *msc_i = cp9->msc[dsq[i]];
  double endsc = -eslINFINITY;
  double sc;
  int k, kn, kx, kpcur, kpprv;

  if(kmin[i] == 0) {
    mc[0]  = -eslINFINITY;
    dc[0]  = -eslINFINITY;
    ec[0]  = -eslINFINITY;
    sc = cp9_chk_dlogsum(cp9_chk_dlogsum(mp[0] + Scorify(CP9TSC(cp9O_MI,0)),
                               ip[0] + Scorify(CP9TSC(cp9O_II,0))),
                    dp[0] + Scorify(CP9TSC(cp9O_DI,0)));
    ic[0] = sc + Scorify(isc_i[0]);
    kn = 1;
  }
  else {
    kn = kmin[i];
  }

  /* match */
  kn = ESL_MAX(kn, (kmin[i-1]+1));
  kx = ESL_MIN(kmax[i], kmax[i-1]+1);

  for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mc[kpcur] = -eslINFINITY;
  for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mc[kpcur] = -eslINFINITY;
  for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) ec[kpcur] = -eslINFINITY;
  for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) ec[kpcur] = -eslINFINITY;

  kpcur = kn - kmin[i];
  kpprv = kn - kmin[i-1];
  for (k = kn; k <= kx; k++, kpcur++, kpprv++) {
    sc = cp9_chk_dlogsum(cp9_chk_dlogsum(mp[kpprv-1] + Scorify(CP9TSC(cp9O_MM,k-1)),
                               ip[kpprv-1] + Scorify(CP9TSC(cp9O_IM,k-1))),
                    dp[kpprv-1] + Scorify(CP9TSC(cp9O_DM,k-1)));
    if(INBAND(i-1, 0)) {
      assert(kmin[(i-1)] == 0);
      if(mp[0] != -eslINFINITY)
        sc = cp9_chk_dlogsum(sc, mp[0] + Scorify(CP9TSC(cp9O_BM,k)));
    }
    if (cp9->flags & CPLAN9_EL) {
      int c_el, kpprv_el;
      for (c_el = 0; c_el < cp9->el_from_ct[k]; c_el++) {
        if (INBAND(i-1, cp9->el_from_idx[k][c_el])) {
          kpprv_el = cp9->el_from_idx[k][c_el] - kmin[i-1];
          sc = cp9_chk_dlogsum(sc, ep[kpprv_el]);
        }
      }
    }
    if(sc != -eslINFINITY) {
      mc[kpcur] = sc + Scorify(msc_i[k]);
      endsc = cp9_chk_dlogsum(endsc, mc[kpcur] + Scorify(CP9TSC(cp9O_ME,k)));
    }
    else {
      mc[kpcur] = -eslINFINITY;
    }
    {
      double el_sc = -eslINFINITY;
      if ((cp9->flags & CPLAN9_EL) && cp9->has_el[k]) {
        el_sc = mc[kpcur] + Scorify(CP9TSC(cp9O_MEL, k)); /* M_k -> EL_k */
        if (INBAND(i-1, k)) {                              /* EL self-loop */
          int kpprv_el = k - kmin[i-1];
          el_sc = cp9_chk_dlogsum(el_sc, ep[kpprv_el] + Scorify(cp9->el_selfsc));
        }
      }
      ec[kpcur] = el_sc;
    }
  }

  /* brief 26_0430-134 BM-coverage supplementary pass (double-only; absent in int kernel).
   * Full row band [max(1,kmin[i]),kmax[i]] MINUS the [kn,kx] match range, adding
   * begin + EL-from-into-M (no in-band diagonal predecessor here). */
  if(INBAND(i-1, 0) && mp[0] != -eslINFINITY) {
    int k_lo = ESL_MAX(1, kmin[i]);
    int k_hi = kmax[i];
    for (k = k_lo; k <= k_hi; k++) {
      if(k >= kn && k <= kx) continue;   /* already filled by the match loop */
      kpcur = k - kmin[i];

      sc = mp[0] + Scorify(CP9TSC(cp9O_BM,k));  /* begin; diagonal predecessor out of band here */

      if (cp9->flags & CPLAN9_EL) {
        int c_el, kpprv_el;
        for (c_el = 0; c_el < cp9->el_from_ct[k]; c_el++) {
          if (INBAND(i-1, cp9->el_from_idx[k][c_el])) {
            kpprv_el = cp9->el_from_idx[k][c_el] - kmin[i-1];
            sc = cp9_chk_dlogsum(sc, ep[kpprv_el]);
          }
        }
      }

      if(sc != -eslINFINITY) {
        mc[kpcur] = sc + Scorify(msc_i[k]);
        endsc = cp9_chk_dlogsum(endsc, mc[kpcur] + Scorify(CP9TSC(cp9O_ME,k)));
      }
      else {
        mc[kpcur] = -eslINFINITY;
      }

      {
        double el_sc = -eslINFINITY;
        if ((cp9->flags & CPLAN9_EL) && cp9->has_el[k]) {
          el_sc = mc[kpcur] + Scorify(CP9TSC(cp9O_MEL, k)); /* M_k -> EL_k */
          if (INBAND(i-1, k)) {                              /* EL self-loop */
            int kpprv_el = k - kmin[i-1];
            el_sc = cp9_chk_dlogsum(el_sc, ep[kpprv_el] + Scorify(cp9->el_selfsc));
          }
        }
        ec[kpcur] = el_sc;
      }
    }
  }

  /* insert */
  kn = ESL_MAX(kmin[i], kmin[i-1]);
  kx = ESL_MIN(kmax[i], kmax[i-1]);
  for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) ic[kpcur] = -eslINFINITY;
  for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) ic[kpcur] = -eslINFINITY;
  kpcur = kn - kmin[i];
  kpprv = kn - kmin[i-1];
  for (k = kn; k <= kx; k++, kpcur++, kpprv++) {
    sc = cp9_chk_dlogsum(cp9_chk_dlogsum(mp[kpprv] + Scorify(CP9TSC(cp9O_MI,k)),
                               ip[kpprv] + Scorify(CP9TSC(cp9O_II,k))),
                    dp[kpprv] + Scorify(CP9TSC(cp9O_DI,k)));
    if(sc != -eslINFINITY) ic[kpcur] = sc + Scorify(isc_i[k]);
    else                   ic[kpcur] = -eslINFINITY;
  }

  /* delete */
  kn = kmin[i]+1;
  for (kpcur = 0; kpcur < (kn - kmin[i]); kpcur++) dc[kpcur] = -eslINFINITY;
  kpcur = kn - kmin[i];
  for (k = kn; k <= kmax[i]; k++, kpcur++) {
    sc = cp9_chk_dlogsum(cp9_chk_dlogsum(mc[kpcur-1] + Scorify(CP9TSC(cp9O_MD,k-1)),
                               ic[kpcur-1] + Scorify(CP9TSC(cp9O_ID,k-1))),
                    dc[kpcur-1] + Scorify(CP9TSC(cp9O_DD,k-1)));
    dc[kpcur] = sc;
  }

  if(INBAND(i, M)) {
    endsc = cp9_chk_dlogsum(cp9_chk_dlogsum(endsc, dc[M-kmin[i]] + Scorify(CP9TSC(cp9O_DM,M))),
                       ic[M-kmin[i]] + Scorify(CP9TSC(cp9O_IM,M)));
  }
  /* erow[i] = endsc; EL-from-M+1 accumulation is UNGATED in the double original
   * (does not require INBAND(i,M)). erow is dead for band derivation. */
  if (cp9->flags & CPLAN9_EL) {
    int c_el;
    for (c_el = 0; c_el < cp9->el_from_ct[M+1]; c_el++) {
      if (INBAND(i, cp9->el_from_idx[M+1][c_el])) {
        int kpel = cp9->el_from_idx[M+1][c_el] - kmin[i];
        endsc = cp9_chk_dlogsum(endsc, ec[kpel]);
      }
    }
  }
  if(ret_erow) *ret_erow = endsc;
}

/* Backward row L (init). Verbatim from cp9_BackwardP7BF i=L block. */
static void
cp9_chk_bwd_rowLF(CP9_t *cp9, ESL_DSQ *dsq, int L, int *kmin, int *kmax, int M,
                  double *mc, double *ic, double *dc, double *ec)
{
  int const *tsc = cp9->otsc;
  int i = L;
  int k, c, kpcur, kpcur_el, kx;

  kpcur = 0;
  for (k = kmin[i]; k <= kmax[i]; k++, kpcur++) ec[kpcur] = -eslINFINITY;
  if(cp9->flags & CPLAN9_EL) {
    for(c = 0; c < cp9->el_from_ct[cp9->M+1]; c++)
      if(INBAND(i, cp9->el_from_idx[M+1][c])) {
        kpcur_el = cp9->el_from_idx[M+1][c] - kmin[i];
        ec[kpcur_el] = 0.;
      }
  }

  if(INBAND(i, M)) {
    assert(M == kmax[i]);
    kpcur = M-kmin[i];
    mc[kpcur]  = 0. +
      cp9_chk_dlogsum(ec[kpcur] + Scorify(CP9TSC(cp9O_MEL, M)),
                 Scorify(CP9TSC(cp9O_ME,M)));
    mc[kpcur] += Scorify(cp9->msc[dsq[i]][M]);
    ic[kpcur]  = 0. + Scorify(CP9TSC(cp9O_IM,M));
    ic[kpcur] += Scorify(cp9->isc[dsq[i]][M]);
    dc[kpcur]  = Scorify(CP9TSC(cp9O_DM,M));
    kx = M-1;
    kpcur--;
  }
  else { kx = kmax[i]; kpcur = kmax[i]-kmin[i]; }

  for (k = kx; k >= kmin[i]; k--, kpcur--)
    {
      mc[kpcur]  = 0 + Scorify(CP9TSC(cp9O_ME,k));
      if(INBAND(i, k+1)) {
        mc[kpcur]  = cp9_chk_dlogsum(mc[kpcur], dc[kpcur+1] + Scorify(CP9TSC(cp9O_MD,k)));
      }
      if(cp9->flags & CPLAN9_EL)
        mc[kpcur]  = cp9_chk_dlogsum(mc[kpcur], ec[kpcur] + Scorify(CP9TSC(cp9O_MEL,k)));
      mc[kpcur] += Scorify(cp9->msc[dsq[i]][k]);

      if(INBAND(i, k+1)) {
        ic[kpcur]  = dc[kpcur+1] + Scorify(CP9TSC(cp9O_ID,k));
        ic[kpcur] += Scorify(cp9->isc[dsq[i]][k]);
        dc[kpcur]  = dc[kpcur+1] + Scorify(CP9TSC(cp9O_DD,k));
      }
      else {
        ic[kpcur] = -eslINFINITY;
        dc[kpcur] = -eslINFINITY;
      }
    }

  if(INBAND(i, 0)) {
    mc[0]  = dc[1] + Scorify(CP9TSC(cp9O_MD,0));
    ic[0]  = dc[1] + Scorify(CP9TSC(cp9O_ID,0));
    ic[0] += Scorify(cp9->isc[dsq[i]][0]);
    dc[0]   = -eslINFINITY;
    ec[0]  = -eslINFINITY;
  }
}

/* Backward row i (1..L-1). Verbatim from cp9_BackwardP7BF main loop body.
 * Cur planes mc/ic/dc/ec (row i); next planes mn/in/dn/en (row i+1). */
static void
cp9_chk_bwd_rowF(CP9_t *cp9, ESL_DSQ *dsq, int i, int *kmin, int *kmax, int M,
                 double *mc, double *ic, double *dc, double *ec,
                 double *mn, double *in, double *dn, double *en)
{
  int const *tsc = cp9->otsc;
  int k, c, kpcur, kpprv, kpcur_el, kn, kx, kprv, kprvn, kprvx;

  kpcur = 0;
  for (k = kmin[i]; k <= kmax[i]; k++, kpcur++) ec[kpcur] = -eslINFINITY;

  if(INBAND(i, M)) {
    kpcur = M-kmin[i];
    if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
      ec[kpcur] = ec[kpcur] + Scorify(cp9->el_selfsc);

    if(INBAND(i+1, M)) {
      kpprv = M-kmin[i+1];
      mc[kpcur]  = in[kpprv] + Scorify(CP9TSC(cp9O_MI,M));
      mc[kpcur] += Scorify(cp9->msc[dsq[i]][M]);
      ic[kpcur]  = in[kpprv] + Scorify(CP9TSC(cp9O_II,M));
      ic[kpcur] += Scorify(cp9->isc[dsq[i]][M]);
      dc[kpcur]  = in[kpprv] + Scorify(CP9TSC(cp9O_DI,M));
    }
    else {
      mc[kpcur] = ic[kpcur] = dc[kpcur] = -eslINFINITY;
    }

    if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
      mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], ec[kpcur] + Scorify(CP9TSC(cp9O_MEL,M)));

    if(INBAND(i+1, M)) {
      if(cp9->flags & CPLAN9_EL) {
        for(c = 0; c < cp9->el_from_ct[M]; c++)
          if(INBAND(i, cp9->el_from_idx[M][c])) {
            kpcur_el = cp9->el_from_idx[M][c] - kmin[i];
            ec[kpcur_el] = cp9_chk_dlogsum(ec[kpcur_el], mn[kpprv]);
          }
      }
    }
  }

  /* MATCH: *_k <- M_k+1 */
  kn = ESL_MAX(kmin[i], kmin[i+1]-1);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]-1);

  for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mc[kpcur] = ic[kpcur] = dc[kpcur] = ec[kpcur] = -eslINFINITY;
  for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mc[kpcur] = ic[kpcur] = dc[kpcur] = ec[kpcur] = -eslINFINITY;

  kpcur = kx - kmin[i];
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    {
      if(cp9->flags & CPLAN9_EL) {
        for(c = 0; c < cp9->el_from_ct[k]; c++) {
          if(INBAND(i, cp9->el_from_idx[k][c])) {
            kpcur_el = cp9->el_from_idx[k][c] - kmin[i];
            ec[kpcur_el] = cp9_chk_dlogsum(ec[kpcur_el], mn[kpprv]);
          }
        }
      }
      if(INBAND(i+1, k)) {
        if((cp9->flags & CPLAN9_EL) && (cp9->has_el[k]))
          ec[kpcur] = cp9_chk_dlogsum(ec[kpcur], en[kpprv] + Scorify(cp9->el_selfsc));
      }
      mc[kpcur] = mn[kpprv+1] + Scorify(CP9TSC(cp9O_MM,k));
      ic[kpcur] = mn[kpprv+1] + Scorify(CP9TSC(cp9O_IM,k));
      dc[kpcur] = mn[kpprv+1] + Scorify(CP9TSC(cp9O_DM,k));
    }

  /* INSERTIONS: *_k <- I_k+1 */
  kn = ESL_MAX(kmin[i], kmin[i+1]);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]);
  kpcur = kx - kmin[i];
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    {
      mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], in[kpprv] + Scorify(CP9TSC(cp9O_MI,k)));
      ic[kpcur] = cp9_chk_dlogsum(ic[kpcur], in[kpprv] + Scorify(CP9TSC(cp9O_II,k)));
      dc[kpcur] = cp9_chk_dlogsum(dc[kpcur], in[kpprv] + Scorify(CP9TSC(cp9O_DI,k)));
    }

  /* DELETIONS: *_k <- D_k+1 */
  kn = ESL_MAX(kmin[i], kmin[i]-1);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i]-1);
  kpcur = kx - kmin[i];
  for (k = kx; k >= kn; k--, kpcur--)
    {
      mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], dc[kpcur+1] + Scorify(CP9TSC(cp9O_MD,k)));
      ic[kpcur] = cp9_chk_dlogsum(ic[kpcur], dc[kpcur+1] + Scorify(CP9TSC(cp9O_ID,k)));
      dc[kpcur] = cp9_chk_dlogsum(dc[kpcur], dc[kpcur+1] + Scorify(CP9TSC(cp9O_DD,k)));
      if((cp9->flags & CPLAN9_EL) && cp9->has_el[k])
        mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], ec[kpcur] + Scorify(CP9TSC(cp9O_MEL,k)));
      mc[kpcur] += Scorify(cp9->msc[dsq[i]][k]);
      ic[kpcur] += Scorify(cp9->isc[dsq[i]][k]);
    }
  for(k = kx+1; k <= kmax[i]; k++) {
    kpcur = k - kmin[i];
    if((cp9->flags & CPLAN9_EL) && cp9->has_el[k])
      mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], ec[kpcur] + Scorify(CP9TSC(cp9O_MEL,k)));
    mc[kpcur] += Scorify(cp9->msc[dsq[i]][k]);
    ic[kpcur] += Scorify(cp9->isc[dsq[i]][k]);
  }

  /* special case k == 0 */
  kpcur = 0;
  kpprv = 0 - kmin[i+1];
  if(INBAND(i, 0)) {
    assert(kmin[i] == 0);
    dc[kpcur]  = -eslINFINITY;
    ec[kpcur] = -eslINFINITY;

    ic[kpcur] = -eslINFINITY;
    if(INBAND(i+1, 1)) {
      if(mn[kpprv+1] != -eslINFINITY)
        ic[kpcur] = cp9_chk_dlogsum(ic[kpcur], mn[kpprv+1] + Scorify(CP9TSC(cp9O_IM,0)));
    }
    if(INBAND(i+1, 0)) {
      if(in[kpprv] != -eslINFINITY)
        ic[kpcur] = cp9_chk_dlogsum(ic[kpcur], in[kpprv] + Scorify(CP9TSC(cp9O_II,0)));
    }
    if(INBAND(i, 1)) {
      if(dc[kpcur+1] != -eslINFINITY)
        ic[kpcur] = cp9_chk_dlogsum(ic[kpcur], dc[kpcur+1] + Scorify(CP9TSC(cp9O_ID,0)));
    }

    kprvn = ESL_MAX(1, kmin[i+1]);
    kprvx = kmax[i+1];
    kpprv = kprvx - kmin[i+1];
    mc[kpcur] = -eslINFINITY;
    for(kprv = kprvx; kprv >= kprvn; kprv--, kpprv--) {
      if(mn[kpprv] != -eslINFINITY)
        mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], (mn[kpprv] + Scorify(CP9TSC(cp9O_BM,kprv))));
    }
    k = 0;
    if(INBAND(i+1, 0)) {
      kpprv = k - kmin[i+1];
      if(in[kpprv] != -eslINFINITY) {
        mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], (in[kpprv] + Scorify(CP9TSC(cp9O_MI,0))));
      }
    }
    if(INBAND(i, 1)) {
      if(dc[kpcur+1] != -eslINFINITY) {
        mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], (dc[kpcur+1] + Scorify(CP9TSC(cp9O_MD,0))));
      }
    }
  }
}

/* Backward row 0. Verbatim from cp9_BackwardP7BF i==0 block.
 * Cur planes mc/ic/dc/ec (row 0); next planes mn/in/dn/en (row 1). */
static void
cp9_chk_bwd_row0F(CP9_t *cp9, ESL_DSQ *dsq, int *kmin, int *kmax, int M,
                  double *mc, double *ic, double *dc, double *ec,
                  double *mn, double *in, double *dn, double *en)
{
  int const *tsc = cp9->otsc;
  int i = 0;
  int k, kpcur, kpprv, kn, kx, kprv, kprvn, kprvx;

  for (kpcur = 0; kpcur <= kmax[0] - kmin[0]; kpcur++) mc[kpcur] = ic[kpcur] = dc[kpcur] = ec[kpcur] = -eslINFINITY;

  /* D_M(i==0) <- I_M(i==1) */
  if(INBAND(i, M)) {
    kpcur = M - kmin[i];
    kpprv = M - kmin[i+1];
    if(INBAND(i+1, M)) {
      dc[kpcur]  = in[kpprv] + Scorify(CP9TSC(cp9O_DI,M));
    }
  }

  /* D_k(i==0) <- M_k+1(i==1) */
  kn = ESL_MAX(kmin[i], kmin[i+1]-1);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]-1);
  kpcur = kx - kmin[i];
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    dc[kpcur]  = mn[kpprv+1] + Scorify(CP9TSC(cp9O_DM,k));

  /* D_k(i==0) <- I_k(i==1) */
  kn = ESL_MAX(kmin[i], kmin[i+1]);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]);
  kpcur = kx - kmin[i];
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    dc[kpcur] = cp9_chk_dlogsum(dc[kpcur], in[kpprv] + Scorify(CP9TSC(cp9O_DI,k)));

  /* D_k(i==0) <- D_k+1(i==0) */
  kn = ESL_MAX(kmin[i], kmin[i]-1);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i]-1);
  kpcur = kx - kmin[i];
  for (k = kx; k >= kn; k--, kpcur--)
    dc[kpcur] = cp9_chk_dlogsum(dc[kpcur], dc[kpcur+1] + Scorify(CP9TSC(cp9O_DD,k)));

  /* k == 0 */
  k = 0;
  if(INBAND(i, 0)) {
    assert(kmin[i]  == 0);
    ic[0] = -eslINFINITY;
    dc[0]   = -eslINFINITY;
    ec[0]  = -eslINFINITY;
    mc[0] = -eslINFINITY;

    kprvn = ESL_MAX(1, kmin[i+1]);
    kprvx = kmax[i+1];
    kpcur = 0;
    kpprv = kprvx - kmin[i+1];
    mc[kpcur] = -eslINFINITY;
    for(kprv = kprvx; kprv >= kprvn; kprv--, kpprv--) {
      if(mn[kpprv] != -eslINFINITY)
        mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], (mn[kpprv] + Scorify(CP9TSC(cp9O_BM,kprv))));
    }
    k = 0;
    if(INBAND(i+1, 0)) {
      kpprv = k - kmin[i+1];
      if(in[kpprv] != -eslINFINITY) {
        mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], (in[kpprv] + Scorify(CP9TSC(cp9O_MI,0))));
      }
    }
    if(INBAND(i, 1)) {
      if(dc[kpcur+1] != -eslINFINITY) {
        mc[kpcur] = cp9_chk_dlogsum(mc[kpcur], (dc[kpcur+1] + Scorify(CP9TSC(cp9O_MD,0))));
      }
    }
  }
}

/*****************************************************************
 * Checkpoint store + fills.
 *****************************************************************/

typedef struct {
  int   L;
  int   M;
  int   blk;       /* block size B (~sqrt(L)) */
  int   nbnd;      /* # boundaries */
  int  *bnd;       /* [0..nbnd-1] boundary row indices */
  int64_t *off;    /* [0..nbnd-1] cell offset of each boundary row's planes */
  int64_t ncells;  /* total stored cells */
  double *fmmx, *fimx, *fdmx, *felmx;   /* Forward checkpoint planes */
  double *bmmx, *bimx, *bdmx, *belmx;   /* Backward checkpoint planes */
  double fsc;      /* brief 26_0430-154 diag: forward total (erow at i==L), for P154_FBDUMP */
} cp9chkF_t;

static void
cp9chkF_Destroy(cp9chkF_t *s)
{
  if(s == NULL) return;
  if(s->bnd)  free(s->bnd);
  if(s->off)  free(s->off);
  if(s->fmmx) free(s->fmmx);
  if(s->fimx) free(s->fimx);
  if(s->fdmx) free(s->fdmx);
  if(s->felmx)free(s->felmx);
  if(s->bmmx) free(s->bmmx);
  if(s->bimx) free(s->bimx);
  if(s->bdmx) free(s->bdmx);
  if(s->belmx)free(s->belmx);
  free(s);
}

static cp9chkF_t *
cp9chkF_Create(int L, int M, int *kmin, int *kmax, char *errbuf)
{
  int status;
  cp9chkF_t *s = NULL;
  int blk, cap, j, r;
  int64_t cum;

  ESL_ALLOC(s, sizeof(cp9chkF_t));
  s->bnd = NULL; s->off = NULL;
  s->fmmx = s->fimx = s->fdmx = s->felmx = NULL;
  s->bmmx = s->bimx = s->bdmx = s->belmx = NULL;
  s->L = L; s->M = M;

  blk = (int) sqrt((double) (L > 0 ? L : 1));
  if(blk < 1) blk = 1;
  { const char *bs = getenv("CP9_CKPT_BLK"); if(bs != NULL) { blk = atoi(bs); if(blk < 1) blk = 1; } }
  s->blk = blk;

  cap = (L / blk) + 3;
  ESL_ALLOC(s->bnd, sizeof(int)     * cap);
  ESL_ALLOC(s->off, sizeof(int64_t) * cap);

  j = 0; r = 0;
  while(r < L) { s->bnd[j++] = r; r += blk; }
  s->bnd[j++] = L;
  s->nbnd = j;

  cum = 0;
  for(j = 0; j < s->nbnd; j++) {
    r = s->bnd[j];
    s->off[j] = cum;
    cum += (kmax[r] - kmin[r] + 1);
  }
  s->ncells = cum;

  ESL_ALLOC(s->fmmx,  sizeof(double) * s->ncells);
  ESL_ALLOC(s->fimx,  sizeof(double) * s->ncells);
  ESL_ALLOC(s->fdmx,  sizeof(double) * s->ncells);
  ESL_ALLOC(s->felmx, sizeof(double) * s->ncells);
  ESL_ALLOC(s->bmmx,  sizeof(double) * s->ncells);
  ESL_ALLOC(s->bimx,  sizeof(double) * s->ncells);
  ESL_ALLOC(s->bdmx,  sizeof(double) * s->ncells);
  ESL_ALLOC(s->belmx, sizeof(double) * s->ncells);

  return s;

 ERROR:
  cp9chkF_Destroy(s);
  if(errbuf) sprintf(errbuf, "cp9chkF_Create: OOM");
  return NULL;
}

/* Forward-checkpointed fill (double): roll 2 rows over 1..L, store F planes at
 * each boundary. Buffers memset-0 (all-bytes-zero = double 0.0f) before each
 * kernel call, replicating GrowCP9FMatrix's memset-0 (recurrences read some
 * cells before assigning them). */
static int
cp9chkF_FwdFill(cp9chkF_t *s, CP9_t *cp9, ESL_DSQ *dsq, int *kmin, int *kmax, char *errbuf)
{
  int status;
  int M = s->M, L = s->L;
  double *m0,*i0,*d0,*e0, *m1,*i1,*d1,*e1;
  double *mp,*ip,*dp,*ep, *mc,*ic,*dc,*ec;
  int i, j, w;

  ESL_ALLOC(m0, sizeof(double)*(M+1)); ESL_ALLOC(i0, sizeof(double)*(M+1));
  ESL_ALLOC(d0, sizeof(double)*(M+1)); ESL_ALLOC(e0, sizeof(double)*(M+1));
  ESL_ALLOC(m1, sizeof(double)*(M+1)); ESL_ALLOC(i1, sizeof(double)*(M+1));
  ESL_ALLOC(d1, sizeof(double)*(M+1)); ESL_ALLOC(e1, sizeof(double)*(M+1));

  mc=m0; ic=i0; dc=d0; ec=e0;
  memset(mc,0,sizeof(double)*(M+1)); memset(ic,0,sizeof(double)*(M+1)); memset(dc,0,sizeof(double)*(M+1)); memset(ec,0,sizeof(double)*(M+1));
  cp9_chk_fwd_row0F(cp9, kmin, kmax, M, mc, ic, dc, ec, NULL);
  j = 0;
  if(s->bnd[j] == 0) {
    w = kmax[0]-kmin[0]+1;
    memcpy(s->fmmx + s->off[j], mc, sizeof(double)*w);
    memcpy(s->fimx + s->off[j], ic, sizeof(double)*w);
    memcpy(s->fdmx + s->off[j], dc, sizeof(double)*w);
    memcpy(s->felmx+ s->off[j], ec, sizeof(double)*w);
    j++;
  }
  mp=m0; ip=i0; dp=d0; ep=e0;

  for(i = 1; i <= L; i++) {
    if((i & 1) == 1) { mc=m1; ic=i1; dc=d1; ec=e1; }
    else             { mc=m0; ic=i0; dc=d0; ec=e0; }
    memset(mc,0,sizeof(double)*(M+1)); memset(ic,0,sizeof(double)*(M+1)); memset(dc,0,sizeof(double)*(M+1)); memset(ec,0,sizeof(double)*(M+1));
    cp9_chk_fwd_rowF(cp9, dsq, i, kmin, kmax, M, mp, ip, dp, ep, mc, ic, dc, ec, (i==L ? &s->fsc : NULL));
    if(j < s->nbnd && s->bnd[j] == i) {
      w = kmax[i]-kmin[i]+1;
      memcpy(s->fmmx + s->off[j], mc, sizeof(double)*w);
      memcpy(s->fimx + s->off[j], ic, sizeof(double)*w);
      memcpy(s->fdmx + s->off[j], dc, sizeof(double)*w);
      memcpy(s->felmx+ s->off[j], ec, sizeof(double)*w);
      j++;
    }
    mp=mc; ip=ic; dp=dc; ep=ec;
  }

  free(m0);free(i0);free(d0);free(e0);free(m1);free(i1);free(d1);free(e1);
  return eslOK;
 ERROR:
  if(errbuf) sprintf(errbuf, "cp9chkF_FwdFill: OOM");
  return status;
}

/* Backward-checkpointed fill (double): roll 2 rows over L..0, store B planes at
 * each boundary. Returns sc = bmx->mmx[0][0]. */
static int
cp9chkF_BwdFill(cp9chkF_t *s, CP9_t *cp9, ESL_DSQ *dsq, int *kmin, int *kmax, double *ret_sc, char *errbuf)
{
  int status;
  int M = s->M, L = s->L;
  double *m0,*i0,*d0,*e0, *m1,*i1,*d1,*e1;
  double *mn,*in,*dn,*en, *mc,*ic,*dc,*ec;
  int i, j, w;

  ESL_ALLOC(m0, sizeof(double)*(M+1)); ESL_ALLOC(i0, sizeof(double)*(M+1));
  ESL_ALLOC(d0, sizeof(double)*(M+1)); ESL_ALLOC(e0, sizeof(double)*(M+1));
  ESL_ALLOC(m1, sizeof(double)*(M+1)); ESL_ALLOC(i1, sizeof(double)*(M+1));
  ESL_ALLOC(d1, sizeof(double)*(M+1)); ESL_ALLOC(e1, sizeof(double)*(M+1));

  j = s->nbnd - 1;

  /* row L. Buffer = i%2. */
  if(L & 1) { mc=m1; ic=i1; dc=d1; ec=e1; } else { mc=m0; ic=i0; dc=d0; ec=e0; }
  memset(mc,0,sizeof(double)*(M+1)); memset(ic,0,sizeof(double)*(M+1)); memset(dc,0,sizeof(double)*(M+1)); memset(ec,0,sizeof(double)*(M+1));
  cp9_chk_bwd_rowLF(cp9, dsq, L, kmin, kmax, M, mc, ic, dc, ec);
  if(s->bnd[j] == L) {
    w = kmax[L]-kmin[L]+1;
    memcpy(s->bmmx + s->off[j], mc, sizeof(double)*w);
    memcpy(s->bimx + s->off[j], ic, sizeof(double)*w);
    memcpy(s->bdmx + s->off[j], dc, sizeof(double)*w);
    memcpy(s->belmx+ s->off[j], ec, sizeof(double)*w);
    j--;
  }
  mn=mc; in=ic; dn=dc; en=ec; /* next row (L) = whatever buffer row L used (Phase-1 fix: NOT hardcoded m0) */

  for(i = L-1; i >= 1; i--) {
    if(i & 1) { mc=m1; ic=i1; dc=d1; ec=e1; }
    else      { mc=m0; ic=i0; dc=d0; ec=e0; }
    memset(mc,0,sizeof(double)*(M+1)); memset(ic,0,sizeof(double)*(M+1)); memset(dc,0,sizeof(double)*(M+1)); memset(ec,0,sizeof(double)*(M+1));
    cp9_chk_bwd_rowF(cp9, dsq, i, kmin, kmax, M, mc, ic, dc, ec, mn, in, dn, en);
    if(j >= 0 && s->bnd[j] == i) {
      w = kmax[i]-kmin[i]+1;
      memcpy(s->bmmx + s->off[j], mc, sizeof(double)*w);
      memcpy(s->bimx + s->off[j], ic, sizeof(double)*w);
      memcpy(s->bdmx + s->off[j], dc, sizeof(double)*w);
      memcpy(s->belmx+ s->off[j], ec, sizeof(double)*w);
      j--;
    }
    mn=mc; in=ic; dn=dc; en=ec;
  }

  /* row 0 (buffer 0 = m0; next=row1 used buffer 1 -> differ) */
  mc=m0; ic=i0; dc=d0; ec=e0;
  memset(mc,0,sizeof(double)*(M+1)); memset(ic,0,sizeof(double)*(M+1)); memset(dc,0,sizeof(double)*(M+1)); memset(ec,0,sizeof(double)*(M+1));
  cp9_chk_bwd_row0F(cp9, dsq, kmin, kmax, M, mc, ic, dc, ec, mn, in, dn, en);
  *ret_sc = mc[0]; /* bmx->mmx[0][0] */
  if(j >= 0 && s->bnd[j] == 0) {
    w = kmax[0]-kmin[0]+1;
    memcpy(s->bmmx + s->off[j], mc, sizeof(double)*w);
    memcpy(s->bimx + s->off[j], ic, sizeof(double)*w);
    memcpy(s->bdmx + s->off[j], dc, sizeof(double)*w);
    memcpy(s->belmx+ s->off[j], ec, sizeof(double)*w);
    j--;
  }

  free(m0);free(i0);free(d0);free(e0);free(m1);free(i1);free(d1);free(e1);
  return eslOK;
 ERROR:
  if(errbuf) sprintf(errbuf, "cp9chkF_BwdFill: OOM");
  return status;
}

/*****************************************************************
 * FLOAT: Segment materialize + posterior + dual-sweep reduction.
 *****************************************************************/

typedef struct {
  int     maxrows;
  int64_t maxcells;
  double  *fm,*fi,*fd,*fe;
  double  *bm,*bi,*bd,*be;
  double **fmr,**fir,**fdr,**fer;
  double **bmr,**bir,**bdr,**ber;
  double  *pm,*pi,*pd;
} cp9segF_t;

static void
cp9segF_Destroy(cp9segF_t *g)
{
  if(g == NULL) return;
  if(g->fm)free(g->fm); if(g->fi)free(g->fi); if(g->fd)free(g->fd); if(g->fe)free(g->fe);
  if(g->bm)free(g->bm); if(g->bi)free(g->bi); if(g->bd)free(g->bd); if(g->be)free(g->be);
  if(g->fmr)free(g->fmr); if(g->fir)free(g->fir); if(g->fdr)free(g->fdr); if(g->fer)free(g->fer);
  if(g->bmr)free(g->bmr); if(g->bir)free(g->bir); if(g->bdr)free(g->bdr); if(g->ber)free(g->ber);
  if(g->pm)free(g->pm); if(g->pi)free(g->pi); if(g->pd)free(g->pd);
  free(g);
}

static cp9segF_t *
cp9segF_Create(cp9chkF_t *s, int *kmin, int *kmax, char *errbuf)
{
  int status;
  cp9segF_t *g = NULL;
  int j, r, maxw = 0;
  int64_t maxcells = 0;

  ESL_ALLOC(g, sizeof(cp9segF_t));
  g->fm=g->fi=g->fd=g->fe=NULL; g->bm=g->bi=g->bd=g->be=NULL;
  g->fmr=g->fir=g->fdr=g->fer=NULL; g->bmr=g->bir=g->bdr=g->ber=NULL;
  g->pm=g->pi=g->pd=NULL;

  for(j = 1; j < s->nbnd; j++) {
    int a = s->bnd[j-1], b = s->bnd[j];
    int64_t cells = 0;
    for(r = a; r <= b; r++) cells += (kmax[r]-kmin[r]+1);
    if(cells > maxcells) maxcells = cells;
  }
  for(r = 0; r <= s->L; r++) { int w = kmax[r]-kmin[r]+1; if(w > maxw) maxw = w; }
  g->maxrows  = s->blk + 1;
  g->maxcells = maxcells;

  ESL_ALLOC(g->fm, sizeof(double)*maxcells); ESL_ALLOC(g->fi, sizeof(double)*maxcells);
  ESL_ALLOC(g->fd, sizeof(double)*maxcells); ESL_ALLOC(g->fe, sizeof(double)*maxcells);
  ESL_ALLOC(g->bm, sizeof(double)*maxcells); ESL_ALLOC(g->bi, sizeof(double)*maxcells);
  ESL_ALLOC(g->bd, sizeof(double)*maxcells); ESL_ALLOC(g->be, sizeof(double)*maxcells);
  ESL_ALLOC(g->fmr, sizeof(double*)*(g->maxrows+1)); ESL_ALLOC(g->fir, sizeof(double*)*(g->maxrows+1));
  ESL_ALLOC(g->fdr, sizeof(double*)*(g->maxrows+1)); ESL_ALLOC(g->fer, sizeof(double*)*(g->maxrows+1));
  ESL_ALLOC(g->bmr, sizeof(double*)*(g->maxrows+1)); ESL_ALLOC(g->bir, sizeof(double*)*(g->maxrows+1));
  ESL_ALLOC(g->bdr, sizeof(double*)*(g->maxrows+1)); ESL_ALLOC(g->ber, sizeof(double*)*(g->maxrows+1));
  ESL_ALLOC(g->pm, sizeof(double)*maxw); ESL_ALLOC(g->pi, sizeof(double)*maxw); ESL_ALLOC(g->pd, sizeof(double)*maxw);

  return g;
 ERROR:
  cp9segF_Destroy(g);
  if(errbuf) sprintf(errbuf, "cp9segF_Create: OOM");
  return NULL;
}

/* Materialize double segment j = rows [a..b]: fwd from F-ckpt at a, bck from
 * B-ckpt at b. Buffers memset-0 to replicate GrowCP9FMatrix's memset-0. */
static void
cp9segF_Fill(cp9segF_t *g, cp9chkF_t *s, int j, CP9_t *cp9, ESL_DSQ *dsq, int *kmin, int *kmax)
{
  int M = s->M, L = s->L;
  int a = s->bnd[j-1], b = s->bnd[j];
  int r, w;
  int64_t fo, bo, segcells;

  segcells = 0;
  for(r = a; r <= b; r++) segcells += (kmax[r]-kmin[r]+1);
  memset(g->fm,0,sizeof(double)*segcells); memset(g->fi,0,sizeof(double)*segcells);
  memset(g->fd,0,sizeof(double)*segcells); memset(g->fe,0,sizeof(double)*segcells);
  memset(g->bm,0,sizeof(double)*segcells); memset(g->bi,0,sizeof(double)*segcells);
  memset(g->bd,0,sizeof(double)*segcells); memset(g->be,0,sizeof(double)*segcells);

  fo = 0; bo = 0;
  for(r = a; r <= b; r++) {
    w = kmax[r]-kmin[r]+1;
    g->fmr[r-a] = g->fm + fo; g->fir[r-a] = g->fi + fo; g->fdr[r-a] = g->fd + fo; g->fer[r-a] = g->fe + fo;
    g->bmr[r-a] = g->bm + bo; g->bir[r-a] = g->bi + bo; g->bdr[r-a] = g->bd + bo; g->ber[r-a] = g->be + bo;
    fo += w; bo += w;
  }

  /* forward: row a = stored F-ckpt[j-1], recompute a+1..b */
  w = kmax[a]-kmin[a]+1;
  memcpy(g->fmr[0], s->fmmx + s->off[j-1], sizeof(double)*w);
  memcpy(g->fir[0], s->fimx + s->off[j-1], sizeof(double)*w);
  memcpy(g->fdr[0], s->fdmx + s->off[j-1], sizeof(double)*w);
  memcpy(g->fer[0], s->felmx+ s->off[j-1], sizeof(double)*w);
  for(r = a+1; r <= b; r++) {
    cp9_chk_fwd_rowF(cp9, dsq, r, kmin, kmax, M,
                     g->fmr[r-1-a], g->fir[r-1-a], g->fdr[r-1-a], g->fer[r-1-a],
                     g->fmr[r-a],   g->fir[r-a],   g->fdr[r-a],   g->fer[r-a], NULL);
  }

  /* backward: row b = stored B-ckpt[j], recompute b-1..a */
  w = kmax[b]-kmin[b]+1;
  memcpy(g->bmr[b-a], s->bmmx + s->off[j], sizeof(double)*w);
  memcpy(g->bir[b-a], s->bimx + s->off[j], sizeof(double)*w);
  memcpy(g->bdr[b-a], s->bdmx + s->off[j], sizeof(double)*w);
  memcpy(g->ber[b-a], s->belmx+ s->off[j], sizeof(double)*w);
  for(r = b-1; r >= a; r--) {
    if(r == 0)
      cp9_chk_bwd_row0F(cp9, dsq, kmin, kmax, M,
                        g->bmr[0], g->bir[0], g->bdr[0], g->ber[0],
                        g->bmr[1], g->bir[1], g->bdr[1], g->ber[1]);
    else if(r == L)
      cp9_chk_bwd_rowLF(cp9, dsq, L, kmin, kmax, M,
                        g->bmr[r-a], g->bir[r-a], g->bdr[r-a], g->ber[r-a]);
    else
      cp9_chk_bwd_rowF(cp9, dsq, r, kmin, kmax, M,
                       g->bmr[r-a],   g->bir[r-a],   g->bdr[r-a],   g->ber[r-a],
                       g->bmr[r+1-a], g->bir[r+1-a], g->bdr[r+1-a], g->ber[r+1-a]);
  }
}

/* Compute double posterior row i into g->pm/pi/pd. Verbatim from the inline
 * posterior in cp9_FB2HMMBandsP7BF / cp9_PosteriorP7BF. */
static void
cp9segF_PostRow(cp9segF_t *g, CP9_t *hmm, ESL_DSQ *dsq, int i, int *kmin, int *kmax, double sc,
                double *fmr, double *fir, double *fdr, double *bmr, double *bir, double *bdr)
{
  int k, kp, kn, kx;
  double *pm=g->pm, *pi=g->pi, *pd=g->pd;

  if(i == 0) {
    pm[0] = fmr[0] + bmr[0] - sc;
    pi[0] = -eslINFINITY;
    pd[0] = -eslINFINITY;
    kn = ESL_MAX(kmin[0], 1);
    kx = kmax[0];
    kp = kn - kmin[0];
    for(k = kn; k <= kx; k++, kp++) {
      pm[kp] = -eslINFINITY;
      pi[kp] = -eslINFINITY;
      pd[kp] = fdr[kp] + bdr[kp] - sc;
    }
    return;
  }

  k = 0;
  if(INBAND(i,0)) {
    kp = 0;
    pm[kp] = ESL_MAX(fmr[kp] + bmr[kp] - sc, -eslINFINITY);
    pi[kp] = ESL_MAX(fir[kp] + bir[kp] - Scorify(hmm->isc[dsq[i]][0]) - sc, -eslINFINITY);
    pd[kp] = -eslINFINITY;
  }
  kn = ESL_MAX(kmin[i], 1);
  kx = kmax[i];
  kp = kn - kmin[i];
  for(k = kn; k <= kx; k++, kp++) {
    pm[kp] = ESL_MAX(fmr[kp] + bmr[kp] - Scorify(hmm->msc[dsq[i]][k]) - sc, -eslINFINITY);
    pi[kp] = ESL_MAX(fir[kp] + bir[kp] - Scorify(hmm->isc[dsq[i]][k]) - sc, -eslINFINITY);
    pd[kp] = ESL_MAX(fdr[kp] + bdr[kp] - sc, -eslINFINITY);
  }
}

/* brief 26_0430-310 (M3): detect a no-parse band instead of letting it poison
 * every posterior.
 *
 * A band deriver can emit a p7 band (kmin/kmax) that contains NO complete CP9
 * parse of the sequence.  When that happens the checkpointed CP9 F/B returns a
 * total score of -inf, every posterior below becomes -inf - (-inf) = NaN, no
 * band edge ever crosses threshold, and cp9b is left as all -1 sentinels.
 * Nothing on this path notices; the first thing that does is an M0 sanity check
 * in a different file, hundreds of lines away, whose message says nothing about
 * bands.
 *
 * THIS IS A DIAGNOSTIC, NOT A FIX.  It does not stop a deriver emitting a
 * no-parse band and it does not recover the alignment; it only makes the
 * resulting failure name its own cause at the point of detection.  The
 * underlying deriver defect is untouched.
 *
 * Why FAIL here, when the sibling p7 banded-decoding path CLAMPS instead:
 *   - p7_GDecodingBanded()/p7b_decode_row() clamp a non-finite INDIVIDUAL CELL
 *     posterior.  One DP cell is degenerate; the parse is still fine, so the
 *     right response is to clamp and continue.
 *   - here, sc/fsc are the WHOLE-SEQUENCE NORMALIZERS -- subtracted uniformly
 *     from every cell in cp9segF_PostRow(), with no per-cell role.  Non-finite
 *     means there is no parse under this band AT ALL, so every downstream
 *     number is meaningless and clamping would only manufacture a plausible
 *     wrong answer.
 * The two different responses to "a non-finite value appeared in a banded
 * posterior computation" are deliberate, not an inconsistency to harmonise.
 *
 * Unconditional failure is safe because this _chk family is not reachable from
 * the search pipeline (cm_pipeline.c makes no _chk calls; search uses the
 * non-checkpointed entry points).  In search, "no complete parse under this
 * band" is an ordinary non-match outcome and a hard failure would be wrong; in
 * alignment every sequence must produce a parse, so it is unambiguously an
 * error.
 */
static int
cp9_chk_noparse_check(char *errbuf, const char *where, int L, int M,
                      int *kmin, int *kmax, double fsc, double bsc)
{
  double ncells = 0., full;
  int    i, nempty = 0, minw = M+1, maxw = 0;

  if(isfinite(fsc) && isfinite(bsc)) return eslOK;

  for(i = 0; i <= L; i++) {
    int w = kmax[i] - kmin[i] + 1;
    if(w <= 0) { nempty++; w = 0; }
    ncells += (double) w;
    if(w < minw) minw = w;
    if(w > maxw) maxw = w;
  }
  full = ((double) L + 1.) * ((double) M + 1.);

  /* The detail goes to stderr, not errbuf: eslERRBUFSIZE is 128 bytes, far too
   * small to carry it, and silently truncating the explanation would recreate
   * the illegibility this guard exists to remove. errbuf gets a short form so
   * whatever finally prints it still names the cause. */
  fprintf(stderr,
          "\nERROR: no-parse p7 band in checkpointed CP9 band derivation (%s).\n"
          "       This p7 band contains no complete parse of the sequence, so the\n"
          "       checkpointed CP9 Forward/Backward totals are not finite (fwd=%g bwd=%g);\n"
          "       every posterior would be NaN and no CP9 band would be set.\n"
          "       band: L=%d M=%d cells=%.0f of %.0f (cover=%.6f)\n"
          "       band width: min=%d max=%d mean=%.1f, %d empty row(s)\n",
          where, fsc, bsc, L, M, ncells, full,
          (full > 0. ? ncells / full : 0.), minw, maxw,
          ncells / ((double) L + 1.), nempty);

  ESL_FAIL(eslENORESULT, errbuf,
           "no-parse p7 band: no complete parse under this band (CP9 F/B total not finite)");
}

/* brief 26_0430-310: which deriver produced the band.  Reported separately from
 * cp9_chk_noparse_check() because that function is a leaf and has no CM_t; the
 * whole failure mode this brief addresses is a symptom surfacing far from its
 * cause, so a message that names the no-parse band but not the deriver that
 * built it would reproduce half of the original problem.
 */
static const char *
cp9_chk_band_deriver_name(CM_t *cm)
{
  if(cm->p7_use_kmerchain) return "--p7kmerchain";
  if(cm->p7_use_ibv)       return "--p7ibv";
  if(cm->p7_use_pinbridge) return "--p7pinbridge";
  return "default p7 Viterbi-trace";
}

/* brief 26_0430-310: emit the no-parse diagnostic on stderr AND fold the deriver
 * name into errbuf, so the cause is visible whether the caller prints errbuf or
 * discards it.
 */
static void
cp9_chk_noparse_report(CM_t *cm, char *errbuf)
{
  const char *deriver = cp9_chk_band_deriver_name(cm);

  fprintf(stderr, "       band deriver: %s\n", deriver);
  snprintf(errbuf, eslERRBUFSIZE,
           "no-parse p7 band from %s: no complete parse under this band "
           "(CP9 F/B total not finite)", deriver);
}

/* The double checkpointed band reduction: produces cp9b pn_min/pn_max bands and
 * the per-node pocc_arr (match+delete occupancy, streamed in the MIN sweep)
 * exactly as cp9_FB2HMMBandsP7BF + cp9_PredictStartAndEndPositionsP7BF's pocc
 * loop, but with O(sqrt(L)*avg_bw) memory.
 *
 * pocc_arr: caller-allocated [0..M]. Filled here: pocc_arr[k] = sum over
 * i=0..L of (exp(pmx->mmx[i][kp]) + exp(pmx->dmx[i][kp])) for nodes k with a
 * band set, else -1.0 (matches the original's skip+sentinel behavior).
 */
int
cp9_FB2HMMBandsP7BF_chk(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9Bands_t *cp9b,
                        int L, int M, double p_thresh, int *kmin, int *kmax,
                        int debug_level, int do_pnmono, int do_pnmono_print,
                        double *pocc_arr)
{
  int status;
  double thresh = log((1. - p_thresh) / 2.);
  int *nset_m=NULL,*nset_i=NULL,*nset_d=NULL;
  int *xset_m=NULL,*xset_i=NULL,*xset_d=NULL;
  double *mass_m=NULL,*mass_i=NULL,*mass_d=NULL;
  int i, k, kp, kn, kx, j;
  double sc;
  int hmm_is_localized;
  cp9chkF_t *s = NULL;
  cp9segF_t *g = NULL;

  hmm_is_localized = ((hmm->flags & CPLAN9_LOCAL_BEGIN) || (hmm->flags & CPLAN9_LOCAL_END) || (hmm->flags & CPLAN9_EL)) ? TRUE : FALSE;

  ESL_ALLOC(nset_m, sizeof(int)*(M+1)); ESL_ALLOC(nset_i, sizeof(int)*(M+1)); ESL_ALLOC(nset_d, sizeof(int)*(M+1));
  ESL_ALLOC(xset_m, sizeof(int)*(M+1)); ESL_ALLOC(xset_i, sizeof(int)*(M+1)); ESL_ALLOC(xset_d, sizeof(int)*(M+1));
  ESL_ALLOC(mass_m, sizeof(double)*(M+1)); ESL_ALLOC(mass_i, sizeof(double)*(M+1)); ESL_ALLOC(mass_d, sizeof(double)*(M+1));
  esl_vec_DSet(mass_m, M+1, -eslINFINITY); esl_vec_DSet(mass_i, M+1, -eslINFINITY); esl_vec_DSet(mass_d, M+1, -eslINFINITY);
  esl_vec_ISet(nset_m, M+1, FALSE); esl_vec_ISet(nset_i, M+1, FALSE); esl_vec_ISet(nset_d, M+1, FALSE);
  esl_vec_ISet(xset_m, M+1, FALSE); esl_vec_ISet(xset_i, M+1, FALSE); esl_vec_ISet(xset_d, M+1, FALSE);

  /* pocc accumulators (k=0..M). k=0 is the B state, never summed. */
  for(k = 0; k <= M; k++) pocc_arr[k] = 0.0;

  if((s = cp9chkF_Create(L, M, kmin, kmax, errbuf)) == NULL) { status = eslEMEM; goto ERROR; }

  /* DEBUG: double chk kernels vs the non-checkpointed FLOAT path (cp9_*P7BF).
   * Post-154 the chk side is double and the ref is float, so this is a
   * precision cross-check with tolerance, not an exact-equality test. */
  if(getenv("CP9_CKPTF_FBCMP") != NULL) {
    CP9_FMX *fr=CreateCP9FMatrix(L,M), *br=CreateCP9FMatrix(L,M);
    CP9_DMX *fc=CreateCP9DMatrix(L,M), *bc=CreateCP9DMatrix(L,M);
    float rsc; int ii, kk, kpp, fmis=0, bmis=0;
    cp9_ForwardP7BF (hmm, errbuf, fr, dsq, L, kmin, kmax, &rsc);
    cp9_BackwardP7BF(hmm, errbuf, br, dsq, L, kmin, kmax, NULL);
    cp9_chk_fwd_row0F(hmm, kmin, kmax, M, fc->mmx[0], fc->imx[0], fc->dmx[0], fc->elmx[0], NULL);
    for(ii=1; ii<=L; ii++) cp9_chk_fwd_rowF(hmm, dsq, ii, kmin, kmax, M, fc->mmx[ii-1],fc->imx[ii-1],fc->dmx[ii-1],fc->elmx[ii-1], fc->mmx[ii],fc->imx[ii],fc->dmx[ii],fc->elmx[ii], NULL);
    cp9_chk_bwd_rowLF(hmm, dsq, L, kmin, kmax, M, bc->mmx[L],bc->imx[L],bc->dmx[L],bc->elmx[L]);
    for(ii=L-1; ii>=1; ii--) cp9_chk_bwd_rowF(hmm, dsq, ii, kmin, kmax, M, bc->mmx[ii],bc->imx[ii],bc->dmx[ii],bc->elmx[ii], bc->mmx[ii+1],bc->imx[ii+1],bc->dmx[ii+1],bc->elmx[ii+1]);
    cp9_chk_bwd_row0F(hmm, dsq, kmin, kmax, M, bc->mmx[0],bc->imx[0],bc->dmx[0],bc->elmx[0], bc->mmx[1],bc->imx[1],bc->dmx[1],bc->elmx[1]);
    for(ii=0; ii<=L; ii++){ kpp=0; for(kk=kmin[ii]; kk<=kmax[ii]; kk++,kpp++){
#define FBCMP_NE(a,b) (fabs((double)(a)-(double)(b)) > 1e-2)
      if(FBCMP_NE(fr->mmx[ii][kpp],fc->mmx[ii][kpp])||FBCMP_NE(fr->imx[ii][kpp],fc->imx[ii][kpp])||FBCMP_NE(fr->dmx[ii][kpp],fc->dmx[ii][kpp])||FBCMP_NE(fr->elmx[ii][kpp],fc->elmx[ii][kpp])){ if(fmis<12) fprintf(stderr,"#FBCMPF F i=%d k=%d ref[M%g I%g D%g E%g] chk[M%g I%g D%g E%g]\n",ii,kk,fr->mmx[ii][kpp],fr->imx[ii][kpp],fr->dmx[ii][kpp],fr->elmx[ii][kpp],fc->mmx[ii][kpp],fc->imx[ii][kpp],fc->dmx[ii][kpp],fc->elmx[ii][kpp]); fmis++; }
      if(FBCMP_NE(br->mmx[ii][kpp],bc->mmx[ii][kpp])||FBCMP_NE(br->imx[ii][kpp],bc->imx[ii][kpp])||FBCMP_NE(br->dmx[ii][kpp],bc->dmx[ii][kpp])||FBCMP_NE(br->elmx[ii][kpp],bc->elmx[ii][kpp])){ if(bmis<12) fprintf(stderr,"#FBCMPF B i=%d k=%d ref[M%g I%g D%g E%g] chk[M%g I%g D%g E%g]\n",ii,kk,br->mmx[ii][kpp],br->imx[ii][kpp],br->dmx[ii][kpp],br->elmx[ii][kpp],bc->mmx[ii][kpp],bc->imx[ii][kpp],bc->dmx[ii][kpp],bc->elmx[ii][kpp]); bmis++; }
#undef FBCMP_NE
    }}
    fprintf(stderr,"#FBCMPF forward mismatches=%d backward mismatches=%d (L=%d M=%d)\n", fmis, bmis, L, M);
    FreeCP9FMatrix(fr);FreeCP9FMatrix(br);FreeCP9DMatrix(fc);FreeCP9DMatrix(bc);
  }

  if(getenv("CP9_CKPT_VERBOSE") != NULL) {
    int64_t full_ncells = 0; int r;
    for(r = 0; r <= L; r++) full_ncells += (kmax[r]-kmin[r]+1);
    fprintf(stderr, "#CP9_CKPTF L=%d M=%d blk=%d nbnd=%d | non-ckpt 3xCP9_FMX=%.1f MB | ckpt stores=%.1f MB (%.0fx smaller)\n",
            L, M, s->blk, s->nbnd,
            3.0 * 32.0 * (double) full_ncells / 1.0e6,
            2.0 * 32.0 * (double) s->ncells   / 1.0e6,
            (3.0 * (double) full_ncells) / (2.0 * (double) (s->ncells > 0 ? s->ncells : 1)));
  }

  if((status = cp9chkF_FwdFill(s, hmm, dsq, kmin, kmax, errbuf)) != eslOK) goto ERROR;
  if((status = cp9chkF_BwdFill(s, hmm, dsq, kmin, kmax, &sc, errbuf)) != eslOK) goto ERROR;
  /* brief 26_0430-310 (M3): no-parse band guard, BEFORE sc is used as the
   * posterior normalizer in cp9segF_PostRow() below. See cp9_chk_noparse_check(). */
  if((status = cp9_chk_noparse_check(errbuf, "cp9_FB2HMMBandsP7BF_chk", L, M,
                                     kmin, kmax, s->fsc, sc)) != eslOK) goto ERROR;
  /* brief 26_0430-154 diag: double checkpointed CP9 F/B totals (cf. 153 ref fwd/bwd≈7993.144,
   * gap≈+0.00003; float ckpt path gave gap≈-1.687 at norovirus, ≈-35/-45 at HSV/MPXV). */
  if(getenv("P154_FBDUMP") != NULL)
    fprintf(stderr, "P154 FBDUMP L=%d M=%d fwd_d=%.6f bwd_d=%.6f gap_d=%.6f\n",
            L, M, s->fsc, sc, s->fsc - sc);
  if((g = cp9segF_Create(s, kmin, kmax, errbuf)) == NULL) { status = eslEMEM; goto ERROR; }

  /* === MIN sweep: ascending i (0..L). Streams pocc_arr. === */
  for(j = 1; j < s->nbnd; j++) {
    int a = s->bnd[j-1], b = s->bnd[j];
    int istart = (j == 1) ? 0 : a+1;
    cp9segF_Fill(g, s, j, hmm, dsq, kmin, kmax);
    for(i = istart; i <= b; i++) {
      int ri = i - a;
      int kk, kkp;
      cp9segF_PostRow(g, hmm, dsq, i, kmin, kmax, sc,
                      g->fmr[ri], g->fir[ri], g->fdr[ri], g->bmr[ri], g->bir[ri], g->bdr[ri]);
      if(i == 0) {
        if((mass_m[0] = g->pm[0]) > thresh) { cp9b->pn_min_m[0] = 0; nset_m[0] = TRUE; }
        mass_i[0] = -eslINFINITY;
        mass_d[0] = -eslINFINITY;
        kn = ESL_MAX(kmin[0], 1); kx = kmax[0]; kp = kn - kmin[0];
        for(k = kn; k <= kx; k++, kp++) {
          if((mass_d[k] = g->pd[kp]) > thresh) { cp9b->pn_min_d[k] = 0; nset_d[k] = TRUE; }
        }
      }
      else {
        k = 0;
        if(INBAND(i,0)) {
          kp = 0;
          if(! nset_m[0]) { if((mass_m[0] = cp9_chk_dlogsum(mass_m[0], g->pm[kp])) > thresh) { cp9b->pn_min_m[0] = i; nset_m[0] = TRUE; } }
          if(! nset_i[0]) { if((mass_i[0] = cp9_chk_dlogsum(mass_i[0], g->pi[kp])) > thresh) { cp9b->pn_min_i[0] = i; nset_i[0] = TRUE; } }
        }
        kn = ESL_MAX(kmin[i], 1); kx = kmax[i]; kp = kn - kmin[i];
        for(k = kn; k <= kx; k++, kp++) {
          if(! nset_m[k]) { if((mass_m[k] = cp9_chk_dlogsum(mass_m[k], g->pm[kp])) > thresh) { cp9b->pn_min_m[k] = i; nset_m[k] = TRUE; } }
          if(! nset_i[k]) { if((mass_i[k] = cp9_chk_dlogsum(mass_i[k], g->pi[kp])) > thresh) { cp9b->pn_min_i[k] = i; nset_i[k] = TRUE; } }
          if(! nset_d[k]) { if((mass_d[k] = cp9_chk_dlogsum(mass_d[k], g->pd[kp])) > thresh) { cp9b->pn_min_d[k] = i; nset_d[k] = TRUE; } }
        }
      }
      /* pocc streaming: sum exp(pm)+exp(pd) over k>=1 in band, in TWO separate
       * adds (matches cp9_PredictStartAndEndPositionsP7BF's double add order). */
      kkp = ESL_MAX(1, kmin[i]) - kmin[i];
      for(kk = ESL_MAX(1, kmin[i]); kk <= kmax[i]; kk++, kkp++) {
        pocc_arr[kk] += exp(g->pm[kkp]);
        pocc_arr[kk] += exp(g->pd[kkp]);
      }
    }
  }

  /* === MAX sweep: descending i (L..1), then row 0 boundary === */
  esl_vec_DSet(mass_m, M+1, -eslINFINITY); esl_vec_DSet(mass_i, M+1, -eslINFINITY); esl_vec_DSet(mass_d, M+1, -eslINFINITY);
  for(j = s->nbnd - 1; j >= 1; j--) {
    int a = s->bnd[j-1], b = s->bnd[j];
    int istart = b;
    int iend   = (j == 1) ? 1 : a+1;
    cp9segF_Fill(g, s, j, hmm, dsq, kmin, kmax);
    for(i = istart; i >= iend; i--) {
      int ri = i - a;
      cp9segF_PostRow(g, hmm, dsq, i, kmin, kmax, sc,
                      g->fmr[ri], g->fir[ri], g->fdr[ri], g->bmr[ri], g->bir[ri], g->bdr[ri]);
      kp = 0;
      for(k = kmin[i]; k <= kmax[i]; k++, kp++) {
        if(! xset_m[k]) { if((mass_m[k] = cp9_chk_dlogsum(mass_m[k], g->pm[kp])) > thresh) { cp9b->pn_max_m[k] = i; xset_m[k] = TRUE; } }
        if(! xset_i[k]) { if((mass_i[k] = cp9_chk_dlogsum(mass_i[k], g->pi[kp])) > thresh) { cp9b->pn_max_i[k] = i; xset_i[k] = TRUE; } }
        if(! xset_d[k]) { if((mass_d[k] = cp9_chk_dlogsum(mass_d[k], g->pd[kp])) > thresh) { cp9b->pn_max_d[k] = i; xset_d[k] = TRUE; } }
      }
      if(j == 1 && i == 1) {
        cp9segF_PostRow(g, hmm, dsq, 0, kmin, kmax, sc,
                        g->fmr[0], g->fir[0], g->fdr[0], g->bmr[0], g->bir[0], g->bdr[0]);
        if(INBAND(0,0)) {
          if(! xset_m[0]) { if((mass_m[0] = cp9_chk_dlogsum(mass_m[0], g->pm[0])) > thresh) { cp9b->pn_max_m[0] = 0; xset_m[0] = TRUE; } }
        }
        kn = ESL_MAX(kmin[0], 1); kx = kmax[0]; kp = kn - kmin[0];
        for(k = kn; k <= kx; k++, kp++) {
          if(!xset_d[k]) { if((mass_d[k] = cp9_chk_dlogsum(mass_d[k], g->pd[kp])) > thresh) { cp9b->pn_max_d[k] = 0; xset_d[k] = TRUE; } }
        }
      }
    }
  }

  /* finalize (verbatim from cp9_FB2HMMBandsP7BF) */
  {
    int mset, dset;
    for(k = 0; k <= M; k++) {
      mset = dset = TRUE;
      if(((! nset_m[k])) || (! xset_m[k]) || (cp9b->pn_max_m[k] < cp9b->pn_min_m[k])) { cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = -1; mset = FALSE; }
      if(((! nset_i[k])) || (! xset_i[k]) || (cp9b->pn_max_i[k] < cp9b->pn_min_i[k])) { cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = -1; }
      if(((! nset_d[k])) || (! xset_d[k]) || (cp9b->pn_max_d[k] < cp9b->pn_min_d[k])) { cp9b->pn_min_d[k] = cp9b->pn_max_d[k] = -1; dset = FALSE; }
      if((!hmm_is_localized) && (mset == FALSE && dset == FALSE)) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "node: %d match nor delete HMM state bands were set in non-localized, non-scanning HMM, lower tau (should be << 0.5).\n", k);
    }
    cp9b->pn_min_d[0] = -1;
    cp9b->pn_max_d[0] = -1;
  }

  if(do_pnmono) pn_match_bands_enforce_monotone(cp9b->pn_min_m, cp9b->pn_max_m, M, L, do_pnmono_print, "fb2hmm_p7bf_chk");

  /* Mask pocc_arr: nodes with no band set -> -1.0 (sentinel), matching the
   * original's skip behavior; k=0 (B state) -> -1.0. */
  pocc_arr[0] = -1.0;
  for(k = 1; k <= M; k++) {
    if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) pocc_arr[k] = -1.0;
  }

  /* DEBUG: compare chk pn arrays vs the non-checkpointed double reduction. */
  if(getenv("CP9_CKPTF_DEBUG") != NULL) {
    int *tmn,*tmx_m,*tin,*tix,*tdn,*tdx;
    CP9_FMX *fmx=NULL,*bmx=NULL,*pmx=NULL; float rsc; int kk, nmis=0;
    /* brief 26_0430-186: optional full CSV dump (every mismatch, not just the first
     * 20) to a file, so a follow-on can characterize the divergence shape
     * (uniform vs localized vs boundary-correlated) instead of eyeballing a
     * truncated stderr sample. Env-gated, inert unless CP9_CKPTFDBG_CSV is
     * set to an output path. */
    FILE *csv = NULL;
    char *csv_path = getenv("CP9_CKPTFDBG_CSV");
    if(csv_path != NULL) csv = fopen(csv_path, "w");
    if(csv != NULL) fprintf(csv, "k,type,ref_min,ref_max,chk_min,chk_max\n");
    ESL_ALLOC(tmn, sizeof(int)*(M+1)); ESL_ALLOC(tmx_m, sizeof(int)*(M+1));
    ESL_ALLOC(tin, sizeof(int)*(M+1)); ESL_ALLOC(tix, sizeof(int)*(M+1));
    ESL_ALLOC(tdn, sizeof(int)*(M+1)); ESL_ALLOC(tdx, sizeof(int)*(M+1));
    for(kk=0;kk<=M;kk++){ tmn[kk]=cp9b->pn_min_m[kk]; tmx_m[kk]=cp9b->pn_max_m[kk]; tin[kk]=cp9b->pn_min_i[kk]; tix[kk]=cp9b->pn_max_i[kk]; tdn[kk]=cp9b->pn_min_d[kk]; tdx[kk]=cp9b->pn_max_d[kk]; }
    fmx=CreateCP9FMatrix(1,M); bmx=CreateCP9FMatrix(1,M); pmx=CreateCP9FMatrix(1,M);
    cp9_ForwardP7BF (hmm, errbuf, fmx, dsq, L, kmin, kmax, &rsc);
    cp9_BackwardP7BF(hmm, errbuf, bmx, dsq, L, kmin, kmax, NULL);
    /* brief 26_0430-186: compare float-path total F/B scores (rsc = fwd total, bmx->mmx[0][0]
     * = bwd total, both accumulated via p7_FLogsum over L residues) against the
     * checkpointed double-path totals already in scope (s->fsc, sc) -- a direct
     * check of whether the float path's cumulative FLogsum drift over this L is
     * large enough to explain the pn-band collapse below. */
    fprintf(stderr, "#CKPTFDBG_TOTALS float_fwd=%.6f float_bwd=%.6f dbl_fwd=%.6f dbl_bwd=%.6f fwd_gap=%.6f bwd_gap=%.6f (L=%d)\n",
            rsc, bmx->mmx[0][0], s->fsc, sc, (double)rsc - s->fsc, (double)bmx->mmx[0][0] - sc, L);
    cp9_FB2HMMBandsP7BF(hmm, errbuf, dsq, fmx, bmx, pmx, cp9b, L, M, p_thresh, 0, kmin, kmax, 0, do_pnmono, do_pnmono_print);
    for(kk=0;kk<=M;kk++){
      if(tmn[kk]!=cp9b->pn_min_m[kk]||tmx_m[kk]!=cp9b->pn_max_m[kk]) { if(nmis<20) fprintf(stderr,"#CKPTFDBG k=%d M ref[%d,%d] chk[%d,%d]\n",kk,cp9b->pn_min_m[kk],cp9b->pn_max_m[kk],tmn[kk],tmx_m[kk]); if(csv!=NULL) fprintf(csv,"%d,M,%d,%d,%d,%d\n",kk,cp9b->pn_min_m[kk],cp9b->pn_max_m[kk],tmn[kk],tmx_m[kk]); nmis++; }
      if(tin[kk]!=cp9b->pn_min_i[kk]||tix[kk]!=cp9b->pn_max_i[kk]) { if(nmis<20) fprintf(stderr,"#CKPTFDBG k=%d I ref[%d,%d] chk[%d,%d]\n",kk,cp9b->pn_min_i[kk],cp9b->pn_max_i[kk],tin[kk],tix[kk]); if(csv!=NULL) fprintf(csv,"%d,I,%d,%d,%d,%d\n",kk,cp9b->pn_min_i[kk],cp9b->pn_max_i[kk],tin[kk],tix[kk]); nmis++; }
      if(tdn[kk]!=cp9b->pn_min_d[kk]||tdx[kk]!=cp9b->pn_max_d[kk]) { if(nmis<20) fprintf(stderr,"#CKPTFDBG k=%d D ref[%d,%d] chk[%d,%d]\n",kk,cp9b->pn_min_d[kk],cp9b->pn_max_d[kk],tdn[kk],tdx[kk]); if(csv!=NULL) fprintf(csv,"%d,D,%d,%d,%d,%d\n",kk,cp9b->pn_min_d[kk],cp9b->pn_max_d[kk],tdn[kk],tdx[kk]); nmis++; }
    }
    fprintf(stderr,"#CKPTFDBG total pn mismatches: %d (M=%d L=%d)\n", nmis, M, L);
    if(csv != NULL) fclose(csv);
    for(kk=0;kk<=M;kk++){ cp9b->pn_min_m[kk]=tmn[kk]; cp9b->pn_max_m[kk]=tmx_m[kk]; cp9b->pn_min_i[kk]=tin[kk]; cp9b->pn_max_i[kk]=tix[kk]; cp9b->pn_min_d[kk]=tdn[kk]; cp9b->pn_max_d[kk]=tdx[kk]; }
    FreeCP9FMatrix(fmx); FreeCP9FMatrix(bmx); FreeCP9FMatrix(pmx);
    free(tmn);free(tmx_m);free(tin);free(tix);free(tdn);free(tdx);
  }

  cp9segF_Destroy(g);
  cp9chkF_Destroy(s);
  free(nset_m);free(nset_i);free(nset_d);
  free(xset_m);free(xset_i);free(xset_d);
  free(mass_m);free(mass_i);free(mass_d);
  return eslOK;

 ERROR:
  if(g) cp9segF_Destroy(g);
  if(s) cp9chkF_Destroy(s);
  if(nset_m)free(nset_m); if(nset_i)free(nset_i); if(nset_d)free(nset_d);
  if(xset_m)free(xset_m); if(xset_i)free(xset_i); if(xset_d)free(xset_d);
  if(mass_m)free(mass_m); if(mass_i)free(mass_i); if(mass_d)free(mass_d);
  return status;
}

/* Function: cp9_FB2HMMBandsP7BF_chk_multi()
 *
 * Brief 26_0430-167 (tau-ratchet single-pass): multi-threshold sibling of
 * cp9_FB2HMMBandsP7BF_chk. The checkpointed double CP9 F/B (FwdFill + BwdFill)
 * and the posterior re-materialization (cp9segF_PostRow) are tau/thresh-
 * independent (brief 26_0430-166 Q3), so we run them ONCE and thread NS distinct
 * tau-derived thresholds through the MIN/MAX sweeps simultaneously. For each
 * materialized posterior row (computed once) we apply the existing "> thresh[t]"
 * band-edge update for every step t. The per-(t,k) "if(!nset[t][k])"
 * short-circuit is preserved exactly, so each step t's mass accumulation is
 * bit-identical to its standalone single-threshold run (the load-bearing
 * determinism property; see CP9_CKPTF_DEBUG-style cross-check in the driver).
 *
 * Inputs:  p_thresh[t] = (1. - tau_t) for step t, t=0..NS-1.
 * Outputs: pn_{min,max}_{m,i,d}_out[t][0..M] = finalized, monotone-enforced
 *          per-step band edges (1..L coords, pre-i0-shift); pocc_out[t][0..M] =
 *          raw streamed pocc masked by step t's bands (matches the single
 *          kernel's masking at lines ~2340). Caller owns all output arrays.
 */
int
cp9_FB2HMMBandsP7BF_chk_multi(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9Bands_t *cp9b,
                              int L, int M, double *p_thresh, int NS, int *kmin, int *kmax,
                              int debug_level, int do_pnmono, int do_pnmono_print,
                              int **pn_min_m_out, int **pn_max_m_out,
                              int **pn_min_i_out, int **pn_max_i_out,
                              int **pn_min_d_out, int **pn_max_d_out,
                              double **pocc_out)
{
  int status;
  double *thresh = NULL;                       /* thresh[t] = log((1-p_thresh[t])/2) */
  int    **nset_m=NULL,**nset_i=NULL,**nset_d=NULL;   /* [t][k] */
  int    **xset_m=NULL,**xset_i=NULL,**xset_d=NULL;
  double **mass_m=NULL,**mass_i=NULL,**mass_d=NULL;
  double  *pocc_raw=NULL;                       /* tau-independent raw streamed pocc */
  int i, k, kp, kn, kx, j, t;
  double sc;
  int hmm_is_localized;
  cp9chkF_t *s = NULL;
  cp9segF_t *g = NULL;

  hmm_is_localized = ((hmm->flags & CPLAN9_LOCAL_BEGIN) || (hmm->flags & CPLAN9_LOCAL_END) || (hmm->flags & CPLAN9_EL)) ? TRUE : FALSE;

  ESL_ALLOC(thresh,   sizeof(double)*NS);
  for(t = 0; t < NS; t++) thresh[t] = log((1. - p_thresh[t]) / 2.);

  /* Per-step accumulators (arrays of NS pointers, each (M+1) long). */
  ESL_ALLOC(nset_m, sizeof(int*)*NS); ESL_ALLOC(nset_i, sizeof(int*)*NS); ESL_ALLOC(nset_d, sizeof(int*)*NS);
  ESL_ALLOC(xset_m, sizeof(int*)*NS); ESL_ALLOC(xset_i, sizeof(int*)*NS); ESL_ALLOC(xset_d, sizeof(int*)*NS);
  ESL_ALLOC(mass_m, sizeof(double*)*NS); ESL_ALLOC(mass_i, sizeof(double*)*NS); ESL_ALLOC(mass_d, sizeof(double*)*NS);
  for(t = 0; t < NS; t++) { nset_m[t]=nset_i[t]=nset_d[t]=NULL; xset_m[t]=xset_i[t]=xset_d[t]=NULL; mass_m[t]=mass_i[t]=mass_d[t]=NULL; }
  for(t = 0; t < NS; t++) {
    ESL_ALLOC(nset_m[t], sizeof(int)*(M+1)); ESL_ALLOC(nset_i[t], sizeof(int)*(M+1)); ESL_ALLOC(nset_d[t], sizeof(int)*(M+1));
    ESL_ALLOC(xset_m[t], sizeof(int)*(M+1)); ESL_ALLOC(xset_i[t], sizeof(int)*(M+1)); ESL_ALLOC(xset_d[t], sizeof(int)*(M+1));
    ESL_ALLOC(mass_m[t], sizeof(double)*(M+1)); ESL_ALLOC(mass_i[t], sizeof(double)*(M+1)); ESL_ALLOC(mass_d[t], sizeof(double)*(M+1));
    esl_vec_DSet(mass_m[t], M+1, -eslINFINITY); esl_vec_DSet(mass_i[t], M+1, -eslINFINITY); esl_vec_DSet(mass_d[t], M+1, -eslINFINITY);
    esl_vec_ISet(nset_m[t], M+1, FALSE); esl_vec_ISet(nset_i[t], M+1, FALSE); esl_vec_ISet(nset_d[t], M+1, FALSE);
    esl_vec_ISet(xset_m[t], M+1, FALSE); esl_vec_ISet(xset_i[t], M+1, FALSE); esl_vec_ISet(xset_d[t], M+1, FALSE);
  }
  ESL_ALLOC(pocc_raw, sizeof(double)*(M+1));
  for(k = 0; k <= M; k++) pocc_raw[k] = 0.0;

  if((s = cp9chkF_Create(L, M, kmin, kmax, errbuf)) == NULL) { status = eslEMEM; goto ERROR; }

  if((status = cp9chkF_FwdFill(s, hmm, dsq, kmin, kmax, errbuf)) != eslOK) goto ERROR;
  if((status = cp9chkF_BwdFill(s, hmm, dsq, kmin, kmax, &sc, errbuf)) != eslOK) goto ERROR;
  /* brief 26_0430-310 (M3): no-parse band guard, BEFORE sc is used as the
   * posterior normalizer in cp9segF_PostRow() below. See cp9_chk_noparse_check(). */
  if((status = cp9_chk_noparse_check(errbuf, "cp9_FB2HMMBandsP7BF_chk_multi", L, M,
                                     kmin, kmax, s->fsc, sc)) != eslOK) goto ERROR;
  /* brief 26_0430-193: same F/B total dump as P154_FBDUMP (single-threshold path,
   * line ~2243), replicated here in the MULTI (tau-ratchet) band-extraction path so
   * genome-scale runs that take this path still report the F/B gap for the LUT-vs-exact
   * precision comparison. Same env var, so only one dump fires per sequence. */
  if(getenv("P154_FBDUMP") != NULL)
    fprintf(stderr, "P154 FBDUMP L=%d M=%d fwd_d=%.6f bwd_d=%.6f gap_d=%.6f\n",
            L, M, s->fsc, sc, s->fsc - sc);
  if((g = cp9segF_Create(s, kmin, kmax, errbuf)) == NULL) { status = eslEMEM; goto ERROR; }

  /* === MIN sweep: ascending i (0..L). Streams pocc_raw once. === */
  for(j = 1; j < s->nbnd; j++) {
    int a = s->bnd[j-1], b = s->bnd[j];
    int istart = (j == 1) ? 0 : a+1;
    cp9segF_Fill(g, s, j, hmm, dsq, kmin, kmax);
    for(i = istart; i <= b; i++) {
      int ri = i - a;
      int kk, kkp;
      cp9segF_PostRow(g, hmm, dsq, i, kmin, kmax, sc,
                      g->fmr[ri], g->fir[ri], g->fdr[ri], g->bmr[ri], g->bir[ri], g->bdr[ri]);
      if(i == 0) {
        for(t = 0; t < NS; t++) {
          if((mass_m[t][0] = g->pm[0]) > thresh[t]) { pn_min_m_out[t][0] = 0; nset_m[t][0] = TRUE; }
          mass_i[t][0] = -eslINFINITY;
          mass_d[t][0] = -eslINFINITY;
        }
        kn = ESL_MAX(kmin[0], 1); kx = kmax[0]; kp = kn - kmin[0];
        for(k = kn; k <= kx; k++, kp++) {
          for(t = 0; t < NS; t++) {
            if((mass_d[t][k] = g->pd[kp]) > thresh[t]) { pn_min_d_out[t][k] = 0; nset_d[t][k] = TRUE; }
          }
        }
      }
      else {
        if(INBAND(i,0)) {
          kp = 0;
          for(t = 0; t < NS; t++) {
            if(! nset_m[t][0]) { if((mass_m[t][0] = cp9_chk_dlogsum(mass_m[t][0], g->pm[kp])) > thresh[t]) { pn_min_m_out[t][0] = i; nset_m[t][0] = TRUE; } }
            if(! nset_i[t][0]) { if((mass_i[t][0] = cp9_chk_dlogsum(mass_i[t][0], g->pi[kp])) > thresh[t]) { pn_min_i_out[t][0] = i; nset_i[t][0] = TRUE; } }
          }
        }
        kn = ESL_MAX(kmin[i], 1); kx = kmax[i]; kp = kn - kmin[i];
        for(k = kn; k <= kx; k++, kp++) {
          for(t = 0; t < NS; t++) {
            if(! nset_m[t][k]) { if((mass_m[t][k] = cp9_chk_dlogsum(mass_m[t][k], g->pm[kp])) > thresh[t]) { pn_min_m_out[t][k] = i; nset_m[t][k] = TRUE; } }
            if(! nset_i[t][k]) { if((mass_i[t][k] = cp9_chk_dlogsum(mass_i[t][k], g->pi[kp])) > thresh[t]) { pn_min_i_out[t][k] = i; nset_i[t][k] = TRUE; } }
            if(! nset_d[t][k]) { if((mass_d[t][k] = cp9_chk_dlogsum(mass_d[t][k], g->pd[kp])) > thresh[t]) { pn_min_d_out[t][k] = i; nset_d[t][k] = TRUE; } }
          }
        }
      }
      /* pocc streaming: tau-independent, accumulate once (two separate adds to
       * match the single kernel's double add order at lines ~2284-2286). */
      kkp = ESL_MAX(1, kmin[i]) - kmin[i];
      for(kk = ESL_MAX(1, kmin[i]); kk <= kmax[i]; kk++, kkp++) {
        pocc_raw[kk] += exp(g->pm[kkp]);
        pocc_raw[kk] += exp(g->pd[kkp]);
      }
    }
  }

  /* === MAX sweep: descending i (L..1), then row 0 boundary === */
  for(t = 0; t < NS; t++) { esl_vec_DSet(mass_m[t], M+1, -eslINFINITY); esl_vec_DSet(mass_i[t], M+1, -eslINFINITY); esl_vec_DSet(mass_d[t], M+1, -eslINFINITY); }
  for(j = s->nbnd - 1; j >= 1; j--) {
    int a = s->bnd[j-1], b = s->bnd[j];
    int istart = b;
    int iend   = (j == 1) ? 1 : a+1;
    cp9segF_Fill(g, s, j, hmm, dsq, kmin, kmax);
    for(i = istart; i >= iend; i--) {
      int ri = i - a;
      cp9segF_PostRow(g, hmm, dsq, i, kmin, kmax, sc,
                      g->fmr[ri], g->fir[ri], g->fdr[ri], g->bmr[ri], g->bir[ri], g->bdr[ri]);
      kp = 0;
      for(k = kmin[i]; k <= kmax[i]; k++, kp++) {
        for(t = 0; t < NS; t++) {
          if(! xset_m[t][k]) { if((mass_m[t][k] = cp9_chk_dlogsum(mass_m[t][k], g->pm[kp])) > thresh[t]) { pn_max_m_out[t][k] = i; xset_m[t][k] = TRUE; } }
          if(! xset_i[t][k]) { if((mass_i[t][k] = cp9_chk_dlogsum(mass_i[t][k], g->pi[kp])) > thresh[t]) { pn_max_i_out[t][k] = i; xset_i[t][k] = TRUE; } }
          if(! xset_d[t][k]) { if((mass_d[t][k] = cp9_chk_dlogsum(mass_d[t][k], g->pd[kp])) > thresh[t]) { pn_max_d_out[t][k] = i; xset_d[t][k] = TRUE; } }
        }
      }
      if(j == 1 && i == 1) {
        cp9segF_PostRow(g, hmm, dsq, 0, kmin, kmax, sc,
                        g->fmr[0], g->fir[0], g->fdr[0], g->bmr[0], g->bir[0], g->bdr[0]);
        if(INBAND(0,0)) {
          for(t = 0; t < NS; t++) {
            if(! xset_m[t][0]) { if((mass_m[t][0] = cp9_chk_dlogsum(mass_m[t][0], g->pm[0])) > thresh[t]) { pn_max_m_out[t][0] = 0; xset_m[t][0] = TRUE; } }
          }
        }
        kn = ESL_MAX(kmin[0], 1); kx = kmax[0]; kp = kn - kmin[0];
        for(k = kn; k <= kx; k++, kp++) {
          for(t = 0; t < NS; t++) {
            if(!xset_d[t][k]) { if((mass_d[t][k] = cp9_chk_dlogsum(mass_d[t][k], g->pd[kp])) > thresh[t]) { pn_max_d_out[t][k] = 0; xset_d[t][k] = TRUE; } }
          }
        }
      }
    }
  }

  /* Per-step finalize + monotone + mask (verbatim per-t from the single kernel). */
  for(t = 0; t < NS; t++) {
    int mset, dset;
    for(k = 0; k <= M; k++) {
      mset = dset = TRUE;
      if(((! nset_m[t][k])) || (! xset_m[t][k]) || (pn_max_m_out[t][k] < pn_min_m_out[t][k])) { pn_min_m_out[t][k] = pn_max_m_out[t][k] = -1; mset = FALSE; }
      if(((! nset_i[t][k])) || (! xset_i[t][k]) || (pn_max_i_out[t][k] < pn_min_i_out[t][k])) { pn_min_i_out[t][k] = pn_max_i_out[t][k] = -1; }
      if(((! nset_d[t][k])) || (! xset_d[t][k]) || (pn_max_d_out[t][k] < pn_min_d_out[t][k])) { pn_min_d_out[t][k] = pn_max_d_out[t][k] = -1; dset = FALSE; }
      if((!hmm_is_localized) && (mset == FALSE && dset == FALSE)) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "node: %d match nor delete HMM state bands were set in non-localized, non-scanning HMM, lower tau (should be << 0.5).\n", k);
    }
    pn_min_d_out[t][0] = -1;
    pn_max_d_out[t][0] = -1;

    if(do_pnmono) pn_match_bands_enforce_monotone(pn_min_m_out[t], pn_max_m_out[t], M, L, do_pnmono_print, "fb2hmm_p7bf_chk_multi");

    /* mask: copy raw pocc, then -1 sentinel for k=0 and nodes with no band set. */
    pocc_out[t][0] = -1.0;
    for(k = 1; k <= M; k++) {
      pocc_out[t][k] = pocc_raw[k];
      if(pn_min_m_out[t][k] == -1 && pn_min_i_out[t][k] == -1 && pn_min_d_out[t][k] == -1) pocc_out[t][k] = -1.0;
    }
  }

  cp9segF_Destroy(g);
  cp9chkF_Destroy(s);
  free(thresh); free(pocc_raw);
  for(t = 0; t < NS; t++) {
    free(nset_m[t]); free(nset_i[t]); free(nset_d[t]);
    free(xset_m[t]); free(xset_i[t]); free(xset_d[t]);
    free(mass_m[t]); free(mass_i[t]); free(mass_d[t]);
  }
  free(nset_m); free(nset_i); free(nset_d);
  free(xset_m); free(xset_i); free(xset_d);
  free(mass_m); free(mass_i); free(mass_d);
  return eslOK;

 ERROR:
  if(g) cp9segF_Destroy(g);
  if(s) cp9chkF_Destroy(s);
  if(thresh) free(thresh);
  if(pocc_raw) free(pocc_raw);
  if(nset_m) { for(t=0;t<NS;t++) if(nset_m[t]) free(nset_m[t]); free(nset_m); }
  if(nset_i) { for(t=0;t<NS;t++) if(nset_i[t]) free(nset_i[t]); free(nset_i); }
  if(nset_d) { for(t=0;t<NS;t++) if(nset_d[t]) free(nset_d[t]); free(nset_d); }
  if(xset_m) { for(t=0;t<NS;t++) if(xset_m[t]) free(xset_m[t]); free(xset_m); }
  if(xset_i) { for(t=0;t<NS;t++) if(xset_i[t]) free(xset_i[t]); free(xset_i); }
  if(xset_d) { for(t=0;t<NS;t++) if(xset_d[t]) free(xset_d[t]); free(xset_d); }
  if(mass_m) { for(t=0;t<NS;t++) if(mass_m[t]) free(mass_m[t]); free(mass_m); }
  if(mass_i) { for(t=0;t<NS;t++) if(mass_i[t]) free(mass_i[t]); free(mass_i); }
  if(mass_d) { for(t=0;t<NS;t++) if(mass_d[t]) free(mass_d[t]); free(mass_d); }
  return status;
}

/*****************************************************************
 * FLOAT: truncated start/end prediction from a precomputed pocc_arr +
 *        the checkpointed truncated band-derivation wrapper.
 *****************************************************************/

/* Parts 1-4 of cp9_PredictStartAndEndPositionsP7BF (cm_p7_band.c), taking the
 * already-streamed pocc_arr instead of recomputing it from a full pmx. Verbatim
 * from the original starting at "Part 1: sp1/sp2". */
static void
cp9_PredictStartAndEndFromPoccF(double *pocc_arr, CP9Bands_t *cp9b, int i0, int j0)
{
  int   k;
  double pocc;

  /* Part 1: sp1/sp2 — leftmost nodes with significant occupancy. */
  k = 1;
  cp9b->sp1 = cp9b->sp2 = -1;
  while(k <= cp9b->hmm_M && (cp9b->sp1 == -1 || cp9b->sp2 == -1)) {
    if(pocc_arr[k] < 0.0) { k++; }
    else {
      pocc = pocc_arr[k];
      if((cp9b->sp1 == -1) && (pocc > cp9b->thresh1)) cp9b->sp1 = k;
      if((cp9b->sp2 == -1) && (pocc > cp9b->thresh2)) cp9b->sp2 = k;
      k++;
    }
  }
  if(k == cp9b->hmm_M+1) {
    if(cp9b->sp1 == -1) { cp9b->sp1 = cp9b->hmm_M+1; }
    if(cp9b->sp2 == -1) { cp9b->sp2 = cp9b->hmm_M+1; }
  }

  /* Part 2: ep1/ep2 — rightmost nodes with significant occupancy. */
  if((cp9b->sp1 == cp9b->hmm_M+1) &&
     (cp9b->sp2 == cp9b->hmm_M+1)) {
    cp9b->ep1 = 0;
    cp9b->ep2 = 0;
  }
  else {
    cp9b->ep1 = cp9b->ep2 = -1;
    k = cp9b->hmm_M;
    while(k >= 1 && (cp9b->ep1 == -1 || cp9b->ep2 == -1)) {
      if(pocc_arr[k] < 0.0) { k--; }
      else {
        pocc = pocc_arr[k];
        if((cp9b->ep1 == -1) && (pocc > cp9b->thresh1)) cp9b->ep1 = k;
        if((cp9b->ep2 == -1) && (pocc > cp9b->thresh2)) cp9b->ep2 = k;
        k--;
      }
    }
    if(k == 0) {
      if(cp9b->ep1 == -1) { cp9b->ep1 = 0; }
      if(cp9b->ep2 == -1) { cp9b->ep2 = 0; }
    }
  }

  /* Parts 3-4: Rmarg/Lmarg derivation — identical to the double original. */

  /* Rmarg_imin */
  if(cp9b->sp1 == cp9b->hmm_M+1) { cp9b->Rmarg_imin = i0; }
  else {
    cp9b->Rmarg_imin = INT_MAX;
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_m[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_m[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_i[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_i[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_d[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_d[cp9b->sp1]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_m[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_m[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_i[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_i[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_d[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_d[cp9b->sp2]);
    if(cp9b->Rmarg_imin == INT_MAX || cp9b->sp1 == (cp9b->hmm_M+1) || cp9b->sp2 == (cp9b->hmm_M+1)) cp9b->Rmarg_imin = i0;
    cp9b->Rmarg_imin = ESL_MAX(i0,   cp9b->Rmarg_imin);
    cp9b->Rmarg_imin = ESL_MIN(j0+1, cp9b->Rmarg_imin);
  }
  /* Rmarg_imax */
  if(cp9b->sp1 == cp9b->hmm_M+1) { cp9b->Rmarg_imax = j0; }
  else {
    cp9b->Rmarg_imax = INT_MIN;
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_m[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_m[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_i[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_i[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_d[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_d[cp9b->sp1]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_m[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_m[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_i[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_i[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_d[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_d[cp9b->sp2]);
    if(cp9b->Rmarg_imax == INT_MIN || cp9b->sp1 == (cp9b->hmm_M+1) || cp9b->sp2 == (cp9b->hmm_M+1)) cp9b->Rmarg_imax = j0+1;
    cp9b->Rmarg_imax = ESL_MAX(i0,   cp9b->Rmarg_imax);
    cp9b->Rmarg_imax = ESL_MIN(j0+1, cp9b->Rmarg_imax);
  }
  /* Lmarg_jmin */
  if(cp9b->ep1 == 0) { cp9b->Lmarg_jmin = i0-1; }
  else {
    cp9b->Lmarg_jmin = INT_MAX;
    if(cp9b->ep1 != 0 && cp9b->pn_min_m[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_m[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_min_i[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_i[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_min_d[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_d[cp9b->ep1]-1);
    if(cp9b->ep2 != 0 && cp9b->pn_min_m[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_m[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_min_i[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_i[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_min_d[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_d[cp9b->ep2]-1);
    if(cp9b->Lmarg_jmin == INT_MAX || cp9b->ep1 == 0 || cp9b->ep2 == 0) cp9b->Lmarg_jmin = i0-1;
    cp9b->Lmarg_jmin = ESL_MAX(i0-1, cp9b->Lmarg_jmin);
    cp9b->Lmarg_jmin = ESL_MIN(j0,   cp9b->Lmarg_jmin);
  }
  /* Lmarg_jmax */
  if(cp9b->ep1 == 0) { cp9b->Lmarg_jmax = j0; }
  else {
    cp9b->Lmarg_jmax = INT_MIN;
    if(cp9b->ep1 != 0 && cp9b->pn_max_m[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_m[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_max_i[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_i[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_max_d[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_d[cp9b->ep1]-1);
    if(cp9b->ep2 != 0 && cp9b->pn_max_m[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_m[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_max_i[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_i[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_max_d[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_d[cp9b->ep2]-1);
    if(cp9b->Lmarg_jmax == INT_MIN || cp9b->ep1 == 0 || cp9b->ep2 == 0) cp9b->Lmarg_jmax = j0;
    cp9b->Lmarg_jmax = ESL_MAX(i0-1, cp9b->Lmarg_jmax);
    cp9b->Lmarg_jmax = ESL_MIN(j0,   cp9b->Lmarg_jmax);
  }
}

/* Function: cp9_FinishBandsFromPnPoccF_chk()
 *
 * Brief 26_0430-167: the shared "band-finishing tail" extracted verbatim from
 * cp9_FBMatrices2BandsP7BF_chk so the single-call path AND the tau-ratchet
 * single-pass driver (cp9_IterateSeq2BandsP7BF_chk_multi) run a BYTE-IDENTICAL
 * tail. Assumes cp9b->pn_{min,max}_{m,i,d} (1..L coords) and <pocc_arr> are
 * already populated for the desired ratchet step, and that the caller has
 * already set cm->tau, cp9b->thresh1/thresh2, cp9b->tau. Shifts bands to
 * i0..j0, predicts sp/ep (+ brief-149 glocal floor) & marginal candidates
 * (trunc) or sets non-trunc valid arrays, then HMM2ij -> GrowHD -> ij2d.
 */
static int
cp9_FinishBandsFromPnPoccF_chk(CM_t *cm, char *errbuf, CP9_t *cp9, CP9Bands_t *cp9b,
                               double *pocc_arr, int *kmin, int *kmax, int L, int i0, int j0,
                               int pass_idx, int debug_level)
{
  int status;
  int do_old_hmm2ij = ((cm->align_opts & CM_ALIGN_HMM2IJOLD) || (cm->search_opts & CM_SEARCH_HMM2IJOLD)) ? TRUE : FALSE;
  int do_trunc      = cm_pli_PassAllowsTruncation(pass_idx);

  /* Step 2b: shift HMM bands from 1..L to i0..j0 coords. */
  if(i0 != 1) {
    int offset = i0 - 1;
    int k;
    for(k = 0; k <= cp9b->hmm_M; k++) {
      if(cp9b->pn_min_m[k] != -1) { cp9b->pn_min_m[k] += offset; cp9b->pn_max_m[k] += offset; }
      if(cp9b->pn_min_i[k] != -1) { cp9b->pn_min_i[k] += offset; cp9b->pn_max_i[k] += offset; }
      if(cp9b->pn_min_d[k] != -1) { cp9b->pn_min_d[k] += offset; cp9b->pn_max_d[k] += offset; }
    }
  }

  /* Step 2c: marginal candidates (trunc) or non-trunc valid arrays. */
  if(do_trunc) {
    cp9_PredictStartAndEndFromPoccF(pocc_arr, cp9b, i0, j0);
    /* brief 26_0430-149 (restored for the ckpt path in brief 26_0430-162; mirrors
     * cm_p7_band.c:5937-5940): in glocal alignment the full (J-mode) parse must
     * always be geometrically available. The thresh1 escalation can retreat ep1
     * below clen (and push sp1 above 1) on models with a decaying posterior-
     * occupancy tail -- e.g. pure-MATL VADR genome models such as NC_001959,
     * where once ep1 < clen every state gets Jvalid[v] = FALSE in
     * cp9_MarginalCandidatesFromStartEndPositions(), excluding the full parse so
     * cm_TrInsideAlignHB() returns "no valid parsetree" in -g mode. Floor sp1 <= 1
     * and ep1 >= clen so the whole model stays J-valid. Scoped to glocal so local
     * mode (which already has a valid root + EL tail escape) stays byte-identical.
     * The floor changes only sp1/ep1, not the Rmarg/Lmarg fields that
     * cp9_PredictStartAndEndFromPoccF already derived from the pre-floor sp1/ep1 --
     * identical to the non-ckpt twin, where the floor likewise follows the predictor. */
    if(! (cm->flags & CMH_LOCAL_BEGIN)) {
      if(cp9b->sp1 > 1)        cp9b->sp1 = 1;
      if(cp9b->ep1 < cm->clen) cp9b->ep1 = cm->clen;
    }
    if((status = cp9_MarginalCandidatesFromStartEndPositions(cm, cp9b, pass_idx, errbuf)) != eslOK) return status;
  }
  else {
    esl_vec_ISet(cp9b->Jvalid, cm->M+1, TRUE);
    esl_vec_ISet(cp9b->Lvalid, cm->M+1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, cm->M+1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, cm->M+1, FALSE);
    /* brief 26_0430-201: do_trunc has an analogous 5'/3'-coverage-gap escape
     * hatch built from cp9_PredictStartAndEndFromPoccF's Lmarg/Rmarg marginal
     * candidates (immediately above); !do_trunc has none. When kmerchain (or
     * any p7-banded deriver) leaves a leading/trailing/internal run of HMM
     * nodes with literally zero posterior evidence (pn_min_{m,i,d}[k] == -1
     * at every substate -- e.g. an unpinned, unanchored model region), the
     * downstream cp9_HMM2ijBands() "brutal hack" (hmmband.c, gated on
     * hmm_is_localized && cm_is_fully_localized) only ever widens CM node 1
     * (or its BIF descendants) to guarantee >=1 valid parse; it does not, and
     * structurally cannot, bridge an interior run of dead nodes elsewhere in
     * the tree. A dead run collapses the model's normal left-to-right MATL/
     * MATR/MATP chain (each node's (i,j) band is a mandatory gateway for
     * every downstream node's reachability), leaving only the degenerate
     * local-begin+EL escape at node 1 as the sole valid parse -- silently
     * discarding a real, correctly-banded parse elsewhere in the tree
     * (confirmed by direct trace on the brief-201 dengue LC436672.1/
     * OR029744.1 reproducers: consensus columns 1-18 were entirely
     * unreachable while columns 19+ carried a normal, tight, correct band).
     * Fix at the source: for a locally-configured CM (matching the brutal
     * hack's own gating condition), widen any such zero-evidence node
     * instead of leaving it as an unreachable sentinel, so
     * cp9_HMM2ijBands()'s ordinary (non-hack) traversal keeps the whole
     * model chain navigable and CYK/Inside remains free to find the real,
     * higher-scoring parse. Bound each dead node by its NEAREST VALID
     * NEIGHBORS on either side (not the full i0..j0 span): an unbounded
     * i0..j0 open range for every dead node also over-widens any trailing
     * dead run (e.g. beyond the last kmerchain-pinned column), which let
     * CYK extend the alignment past the model's real, evidence-backed
     * endpoint (empirically: cm-to jumped from the correct ~10539 all the
     * way to cm_M=10723, the literal last column, when tested with a full
     * i0..j0 widening). Bounding by nearest-neighbor evidence keeps the
     * widened region's admissible (i,j) span consistent with where the
     * model actually has data, matching the legacy (pre-double-ckpt) int
     * F/B kernel's band shape much more closely. do_trunc is untouched
     * (separate branch above); this changes ONLY inputs to the
     * non-truncated (!do_trunc) path, so cp9_HMM2ijBands() and every other
     * caller of the shared pn_* arrays outside this function are
     * unaffected. */
    if(cm->flags & CMH_LOCAL_BEGIN) {
      int k;
      int hmm_M = cp9b->hmm_M;
      int *lb, *ub;
      lb = malloc(sizeof(int) * (hmm_M+1));
      ub = malloc(sizeof(int) * (hmm_M+1));
      if(lb == NULL || ub == NULL) ESL_FAIL(eslEMEM, errbuf, "cp9_FinishBandsFromPnPoccF_chk: OOM allocating dead-node bound arrays");
      /* left-to-right sweep: lb[k] = nearest valid pn_min_m to the left of (or at) k, else i0 */
      {
        int cur = i0;
        for(k = 0; k <= hmm_M; k++) {
          if(cp9b->pn_min_m[k] != -1) cur = cp9b->pn_min_m[k];
          lb[k] = cur;
        }
      }
      /* right-to-left sweep: ub[k] = nearest valid pn_max_m to the right of (or at) k, else j0 */
      {
        int cur = j0;
        for(k = hmm_M; k >= 0; k--) {
          if(cp9b->pn_max_m[k] != -1) cur = cp9b->pn_max_m[k];
          ub[k] = cur;
        }
      }
      for(k = 0; k <= hmm_M; k++) {
        if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) {
          cp9b->pn_min_m[k] = cp9b->pn_min_i[k] = cp9b->pn_min_d[k] = lb[k];
          cp9b->pn_max_m[k] = cp9b->pn_max_i[k] = cp9b->pn_max_d[k] = ub[k];
        }
      }
      free(lb);
      free(ub);
    }
  }

  /* Step 3: HMM bands -> CM bands. */
  if(do_old_hmm2ij) {
    /* brief 26_0430-162: doing_search=FALSE (alignment-mode tight j-bands). */
    if((status = cp9_HMM2ijBands_OLD(cm, errbuf, cm->cp9b, cm->cp9map, i0, j0, FALSE, debug_level)) != eslOK) return status;
  }
  else {
    if((status = cp9_HMM2ijBands(cm, errbuf, cp9, cm->cp9b, cm->cp9map, i0, j0, FALSE, do_trunc, debug_level)) != eslOK) return status;
  }
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) return status;
  ij2d_bands(cm, cp9b, do_trunc, debug_level);

  if(do_trunc && (! (cm->flags & CMH_LOCAL_BEGIN))) {
    /* brief 26_0430-185 (mirrors cm_p7_band.c's cp9_FBMatrices2BandsF, non-ckpt twin):
     * the brief-149 sp1/ep1 floor above forces Jvalid[v] = TRUE for essentially
     * every state without widening the real per-state (j,d) bands computed by
     * ij2d_bands() just above, from the un-floored 1-tau threshold signal. On a
     * sequence with a genuinely, biologically missing region this leaves
     * "phantom valid" states whose real band is empty, which a truncated-
     * alignment traceback can walk into and die on (cm_TrInsideAlignHB() "no
     * valid parsetree found"). Veto Jvalid[v] back to FALSE for any state whose
     * real band is empty at every j in its jband, via the hd_min()/hd_max()
     * recompute-on-demand accessors (brief 26_0430-157) -- never reintroduce flat
     * hdmin[v][]/hdmax[v][] reads here, they're gone. This is the production
     * --p7ibv-ckpt path (brief 26_0430-162), so this fix must land here too, not just
     * in the non-ckpt twin. */
    int v, jp, njp, found;
    for(v = 0; v < cp9b->cm_M; v++) {
      if(! cp9b->Jvalid[v]) continue;
      njp = cp9b->jmax[v] - cp9b->jmin[v] + 1;
      found = FALSE;
      for(jp = 0; jp < njp; jp++) {
        if(hd_min(cp9b, v, jp) <= hd_max(cp9b, v, jp)) { found = TRUE; break; }
      }
      if(! found) cp9b->Jvalid[v] = FALSE;
    }
  }

  return eslOK;
}

/* Function: cp9_FBMatrices2BandsP7BF_chk()
 *
 * Checkpointed double drop-in for cp9_FBMatrices2BandsF (truncated path).
 * Produces the same cp9b bands + sp/ep prediction via cp9_FB2HMMBandsP7BF_chk
 * (which streams pocc_arr) + the identical downstream finishing tail
 * (cp9_FinishBandsFromPnPoccF_chk). No full CP9_FMX matrices.
 */
int
cp9_FBMatrices2BandsP7BF_chk(CM_t *cm, char *errbuf, CP9_t *cp9, ESL_DSQ *dsq, CP9Bands_t *cp9b,
                             int *kmin, int *kmax, int L, int i0, int j0, int pass_idx,
                             int debug_level, int do_pnmono, int do_pnmono_print)
{
  int status;
  int use_sums      = ((cm->align_opts & CM_ALIGN_SUMS) || (cm->search_opts & CM_SEARCH_SUMS)) ? TRUE : FALSE;
  double *pocc_arr   = NULL;

  if(use_sums) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_FBMatrices2BandsP7BF_chk: use_sums not supported.");

  ESL_ALLOC(pocc_arr, sizeof(double) * (cp9b->hmm_M + 1));

  /* Step 1+2: checkpointed double F/B -> HMM bands + streamed pocc_arr. */
  if((status = cp9_FB2HMMBandsP7BF_chk(cp9, errbuf, dsq, cp9b, L, cp9b->hmm_M,
                                       (1.-cm->tau), kmin, kmax, debug_level,
                                       do_pnmono, do_pnmono_print, pocc_arr)) != eslOK) {
    /* brief 26_0430-310: name the deriver; the leaf has no CM_t. */
    if(status == eslENORESULT) cp9_chk_noparse_report(cm, errbuf);
    goto ERROR;
  }
  cp9b->tau = cm->tau;

  if((status = cp9_FinishBandsFromPnPoccF_chk(cm, errbuf, cp9, cp9b, pocc_arr, kmin, kmax,
                                              L, i0, j0, pass_idx, debug_level)) != eslOK) goto ERROR;

  free(pocc_arr);
  return eslOK;

 ERROR:
  if(pocc_arr) free(pocc_arr);
  return status;
}

/* Function: cp9_IterateSeq2BandsP7BF_chk_multi()
 *
 * Brief 26_0430-167: single-pass replacement for the ckpt-truncated tau-ratchet loop
 * in cp9_IterateSeq2BandsP7B (which recomputed the WHOLE checkpointed float F/B
 * on every step, up to ~26 steps). Two phases:
 *
 *   Phase 1 -- evaluate step 0 (the current cm->tau/thresh1/thresh2) via the
 *     single-call cp9_FBMatrices2BandsP7BF_chk. For the common 0-bump case this
 *     returns here, byte-identical to the old loop's first iteration, with ZERO
 *     extra memory/compute (critical: keeps the brief-165 genome capstone path
 *     and its 2-3 GB memory recipe unchanged).
 *
 *   Phase 2 -- if step 0 doesn't fit, enumerate the remaining ratchet grid
 *     (steps 1..NS, mirroring the bump logic exactly) and run ONE checkpointed
 *     F/B + MIN + MAX sweep (cp9_FB2HMMBandsP7BF_chk_multi) that evaluates all
 *     NS thresholds at once, then scan steps in order and break at the first
 *     that fits size_limit (NO binary search: size-vs-step is non-monotone via
 *     deck-validity flips, brief 26_0430-166 Q4). Breaking at first fit leaves cp9b in
 *     the selected step's state automatically (each step's tail fully rederives
 *     cp9b), so no separate "restore" is needed.
 *
 * Returns eslOK with bands fitting size_limit, or eslERANGE if even the
 * all-capped final step still exceeds it (byte-identical terminal behavior to
 * the old loop). *ret_nbump = selected step index (0 = no bump).
 */
int
cp9_IterateSeq2BandsP7BF_chk_multi(CM_t *cm, char *errbuf, CP9_t *cp9, ESL_DSQ *dsq, int L,
                                   int *kmin, int *kmax, int i0, int j0, int pass_idx,
                                   float size_limit, int doing_search, int do_sample, int do_post,
                                   double maxtau, int do_pnmono, int do_pnmono_print,
                                   int *ret_nbump, float *ret_Mb)
{
  int status;
  CP9Bands_t *cp9b = cm->cp9b;
  int M = cp9b->hmm_M;
  int do_trunc = cm_pli_PassAllowsTruncation(pass_idx);
  int debug_level = 0;
  float cp9mx_Mb = 0., hbmx_Mb = 0., tot_Mb;
  int   s, k, t, NS = 0;
  double *tau_grid = NULL, *t1_grid = NULL, *t2_grid = NULL, *p_thresh = NULL;
  int **pnmm = NULL, **pnxm = NULL, **pnmi = NULL, **pnxi = NULL, **pnmd = NULL, **pnxd = NULL;
  double **pocc = NULL;

  /* ---- Phase 1: step 0 (current tau/thresh1/thresh2). ---- */
  if((status = cp9_FBMatrices2BandsP7BF_chk(cm, errbuf, cp9, dsq, cp9b, kmin, kmax, L, i0, j0,
                                            pass_idx, debug_level, do_pnmono, do_pnmono_print)) != eslOK) return status;
  /* brief 26_0821-014: branch on do_trunc, as cp9_IterateSeq2BandsP7B() does. This
   * driver is documented (above) as running a BYTE-IDENTICAL ratchet to that loop,
   * but sized unconditionally with the TRUNCATED estimator, so a non-truncated run
   * ratcheted against a matrix with marginal planes it will never fill. hbmx_Mb is
   * not diagnostic here -- it decides eslOK vs eslERANGE, and so which engine the
   * caller escalates to. The error was conservative (over-estimate => more
   * escalation), which is why it was invisible. */
  if(doing_search) {
    if(do_trunc) { if((status = cm_tr_hb_mx_SizeNeeded(cm, errbuf, cp9b, j0-i0+1, NULL, NULL, NULL, NULL, &hbmx_Mb)) != eslOK) return status; }
    else         { if((status = cm_hb_mx_SizeNeeded   (cm, errbuf, cp9b, j0-i0+1, NULL, &hbmx_Mb)) != eslOK) return status; }
  }
  else {
    if(do_trunc) status = cm_TrAlignSizeNeededHB(cm, errbuf, j0-i0+1, size_limit, do_sample, do_post, NULL, NULL, NULL, &cp9mx_Mb, &hbmx_Mb, &tot_Mb);
    else         status = cm_AlignSizeNeededHB  (cm, errbuf, j0-i0+1, size_limit, do_sample, do_post, NULL, NULL, NULL, &cp9mx_Mb, &hbmx_Mb, &tot_Mb);
    if(status != eslOK && status != eslERANGE) return status;
  }
  if(ret_nbump != NULL) *ret_nbump = 0;
  if(hbmx_Mb < size_limit) { if(ret_Mb != NULL) *ret_Mb = hbmx_Mb; return eslOK; }

  /* brief 26_0430-271 item 1: --mxesc-fixedtau (CM_ALIGN_MXESC_FIXEDTAU) skips the
   * ratchet grid below entirely -- step 0's (current tau/thresh) bands are already
   * populated in cp9b, so the caller (cp9_IterateSeq2BandsP7B, via DispatchSqAlignment)
   * keeps valid-but-wide bands and lets mxesc per-seq engine escalation (tier b/c)
   * carry the memory instead of band-tightening. Byte-identical to the pre-271
   * ratchet when the flag is unset (the default). */
  if(cm->align_opts & CM_ALIGN_MXESC_FIXEDTAU) {
    if(ret_Mb != NULL) *ret_Mb = hbmx_Mb;
    return eslERANGE;
  }

  /* ---- Build the remaining ratchet grid (steps 1..NS), mirroring the bump
   * logic in cp9_IterateSeq2BandsP7B exactly (tau*=2 cap maxtau; thresh1 +=
   * DELTA cap MAX; thresh2 -= DELTA floor MIN; for do_trunc all three move). ---- */
  {
    double tau = cm->tau, th1 = cp9b->thresh1, th2 = cp9b->thresh2;
    int tau_lim = FALSE;
    int th1_lim = (do_trunc) ? FALSE : TRUE;
    int th2_lim = (do_trunc) ? FALSE : TRUE;
    int cap = 64; /* safety bound; real worst case ~25 (brief 26_0430-166 Q1) */
    ESL_ALLOC(tau_grid, sizeof(double)*cap);
    ESL_ALLOC(t1_grid,  sizeof(double)*cap);
    ESL_ALLOC(t2_grid,  sizeof(double)*cap);
    while(! (tau_lim && th1_lim && th2_lim)) {
      if(! tau_lim) { tau *= TAU_MULTIPLIER; if(tau >= maxtau) { tau = maxtau; tau_lim = TRUE; } }
      if(! th1_lim) { th1 += DELTA_CP9BANDS_THRESH1; if(th1 >= MAX_CP9BANDS_THRESH1) { th1 = MAX_CP9BANDS_THRESH1; th1_lim = TRUE; } }
      if(! th2_lim) { th2 -= DELTA_CP9BANDS_THRESH2; if(th2 <= MIN_CP9BANDS_THRESH2) { th2 = MIN_CP9BANDS_THRESH2; th2_lim = TRUE; } }
      if(NS >= cap) ESL_XFAIL(eslEINCONCEIVABLE, errbuf, "cp9_IterateSeq2BandsP7BF_chk_multi: ratchet grid overflow (NS=%d).", NS);
      tau_grid[NS] = tau; t1_grid[NS] = th1; t2_grid[NS] = th2;
      NS++;
    }
  }
  if(NS == 0) { /* step 0 was already all-capped: old loop would break -> eslERANGE */
    if(ret_Mb != NULL) *ret_Mb = hbmx_Mb;
    status = eslERANGE; goto DONE;
  }

  /* ---- Allocate multi-threshold outputs (NS slots). Each pointer array is
   * NULL-initialized immediately after allocation so any OOM mid-allocation
   * leaves a state the unified cleanup can free safely. ---- */
  ESL_ALLOC(p_thresh, sizeof(double)*NS);
  for(t = 0; t < NS; t++) p_thresh[t] = 1. - tau_grid[t];
  ESL_ALLOC(pnmm, sizeof(int*)*NS);    for(t=0;t<NS;t++) pnmm[t]=NULL;
  ESL_ALLOC(pnxm, sizeof(int*)*NS);    for(t=0;t<NS;t++) pnxm[t]=NULL;
  ESL_ALLOC(pnmi, sizeof(int*)*NS);    for(t=0;t<NS;t++) pnmi[t]=NULL;
  ESL_ALLOC(pnxi, sizeof(int*)*NS);    for(t=0;t<NS;t++) pnxi[t]=NULL;
  ESL_ALLOC(pnmd, sizeof(int*)*NS);    for(t=0;t<NS;t++) pnmd[t]=NULL;
  ESL_ALLOC(pnxd, sizeof(int*)*NS);    for(t=0;t<NS;t++) pnxd[t]=NULL;
  ESL_ALLOC(pocc, sizeof(double*)*NS); for(t=0;t<NS;t++) pocc[t]=NULL;
  for(t = 0; t < NS; t++) {
    ESL_ALLOC(pnmm[t], sizeof(int)*(M+1)); ESL_ALLOC(pnxm[t], sizeof(int)*(M+1));
    ESL_ALLOC(pnmi[t], sizeof(int)*(M+1)); ESL_ALLOC(pnxi[t], sizeof(int)*(M+1));
    ESL_ALLOC(pnmd[t], sizeof(int)*(M+1)); ESL_ALLOC(pnxd[t], sizeof(int)*(M+1));
    ESL_ALLOC(pocc[t], sizeof(double)*(M+1));
  }

  /* ---- Phase 2: ONE checkpointed F/B + MIN + MAX sweep over all NS steps. ---- */
  if((status = cp9_FB2HMMBandsP7BF_chk_multi(cp9, errbuf, dsq, cp9b, L, M, p_thresh, NS, kmin, kmax,
                                             debug_level, do_pnmono, do_pnmono_print,
                                             pnmm, pnxm, pnmi, pnxi, pnmd, pnxd, pocc)) != eslOK) {
    /* brief 26_0430-310: name the deriver; the leaf has no CM_t. */
    if(status == eslENORESULT) cp9_chk_noparse_report(cm, errbuf);
    goto DONE;
  }

  /* ---- G3 determinism harness (brief 26_0430-167): for every grid step, recompute the
   * pn arrays + masked pocc the OLD single-call way (cp9_FB2HMMBandsP7BF_chk,
   * which re-runs its own F/B) and assert equality vs the multi-threshold
   * sweep's slot. Catches any accumulator-order divergence. cp9b->pn_* is used
   * as scratch here; the per-step tail below re-sets it from the slots, so no
   * save/restore is needed. Env-gated; off in production. ---- */
  if(getenv("CP9_TAURATCHET_DBG") != NULL) {
    double *dpocc = NULL;
    int nmis = 0, ss, kk;
    ESL_ALLOC(dpocc, sizeof(double)*(M+1));
    for(ss = 0; ss < NS; ss++) {
      if((status = cp9_FB2HMMBandsP7BF_chk(cp9, errbuf, dsq, cp9b, L, M, p_thresh[ss], kmin, kmax,
                                           debug_level, do_pnmono, do_pnmono_print, dpocc)) != eslOK) { free(dpocc); goto DONE; }
      for(kk = 0; kk <= M; kk++) {
        if(cp9b->pn_min_m[kk]!=pnmm[ss][kk] || cp9b->pn_max_m[kk]!=pnxm[ss][kk] ||
           cp9b->pn_min_i[kk]!=pnmi[ss][kk] || cp9b->pn_max_i[kk]!=pnxi[ss][kk] ||
           cp9b->pn_min_d[kk]!=pnmd[ss][kk] || cp9b->pn_max_d[kk]!=pnxd[ss][kk] ||
           dpocc[kk]!=pocc[ss][kk]) {
          if(nmis < 20) fprintf(stderr, "#TAURATCHET_DBG MISMATCH s=%d(step%d) k=%d ref[m %d,%d|i %d,%d|d %d,%d|pocc %g] multi[m %d,%d|i %d,%d|d %d,%d|pocc %g]\n",
              ss, ss+1, kk, cp9b->pn_min_m[kk],cp9b->pn_max_m[kk],cp9b->pn_min_i[kk],cp9b->pn_max_i[kk],cp9b->pn_min_d[kk],cp9b->pn_max_d[kk],dpocc[kk],
              pnmm[ss][kk],pnxm[ss][kk],pnmi[ss][kk],pnxi[ss][kk],pnmd[ss][kk],pnxd[ss][kk],pocc[ss][kk]);
          nmis++;
        }
      }
    }
    free(dpocc);
    fprintf(stderr, "#TAURATCHET_DBG total pn/pocc mismatches across %d grid steps: %d (M=%d L=%d)\n", NS, nmis, M, L);
  }

  /* ---- Per-step tail: scan s=0..NS-1, break at first fit. ---- */
  for(s = 0; s < NS; s++) {
    cm->tau       = tau_grid[s];
    cp9b->thresh1 = t1_grid[s];
    cp9b->thresh2 = t2_grid[s];
    cp9b->tau     = tau_grid[s];
    for(k = 0; k <= M; k++) {
      cp9b->pn_min_m[k] = pnmm[s][k]; cp9b->pn_max_m[k] = pnxm[s][k];
      cp9b->pn_min_i[k] = pnmi[s][k]; cp9b->pn_max_i[k] = pnxi[s][k];
      cp9b->pn_min_d[k] = pnmd[s][k]; cp9b->pn_max_d[k] = pnxd[s][k];
    }
    if((status = cp9_FinishBandsFromPnPoccF_chk(cm, errbuf, cp9, cp9b, pocc[s], kmin, kmax,
                                                L, i0, j0, pass_idx, debug_level)) != eslOK) goto DONE;
    if(doing_search) { /* brief 26_0821-014: branch on do_trunc; see the step-0 site above */
      if(do_trunc) { if((status = cm_tr_hb_mx_SizeNeeded(cm, errbuf, cp9b, j0-i0+1, NULL, NULL, NULL, NULL, &hbmx_Mb)) != eslOK) goto DONE; }
      else         { if((status = cm_hb_mx_SizeNeeded   (cm, errbuf, cp9b, j0-i0+1, NULL, &hbmx_Mb)) != eslOK) goto DONE; }
    }
    else {
      if(do_trunc) status = cm_TrAlignSizeNeededHB(cm, errbuf, j0-i0+1, size_limit, do_sample, do_post, NULL, NULL, NULL, &cp9mx_Mb, &hbmx_Mb, &tot_Mb);
      else         status = cm_AlignSizeNeededHB  (cm, errbuf, j0-i0+1, size_limit, do_sample, do_post, NULL, NULL, NULL, &cp9mx_Mb, &hbmx_Mb, &tot_Mb);
      if(status != eslOK && status != eslERANGE) goto DONE;
    }
    if(ret_nbump != NULL) *ret_nbump = s + 1; /* grid slot s == ratchet step s+1 */
    if(hbmx_Mb < size_limit) break; /* first fit; cp9b now holds this step */
  }
  /* If no step fit, cp9b holds the final all-capped step (s == NS-1). */

  if(ret_Mb != NULL) *ret_Mb = hbmx_Mb;
  status = (hbmx_Mb > size_limit) ? eslERANGE : eslOK;

 ERROR:  /* ESL_ALLOC failures land here; fall through to the same guarded cleanup. */
 DONE:
  if(tau_grid) free(tau_grid);
  if(t1_grid)  free(t1_grid);
  if(t2_grid)  free(t2_grid);
  if(p_thresh) free(p_thresh);
  if(pnmm) { for(t=0;t<NS;t++) if(pnmm[t]) free(pnmm[t]); free(pnmm); }
  if(pnxm) { for(t=0;t<NS;t++) if(pnxm[t]) free(pnxm[t]); free(pnxm); }
  if(pnmi) { for(t=0;t<NS;t++) if(pnmi[t]) free(pnmi[t]); free(pnmi); }
  if(pnxi) { for(t=0;t<NS;t++) if(pnxi[t]) free(pnxi[t]); free(pnxi); }
  if(pnmd) { for(t=0;t<NS;t++) if(pnmd[t]) free(pnmd[t]); free(pnmd); }
  if(pnxd) { for(t=0;t<NS;t++) if(pnxd[t]) free(pnxd[t]); free(pnxd); }
  if(pocc) { for(t=0;t<NS;t++) if(pocc[t]) free(pocc[t]); free(pocc); }
  return status;
}

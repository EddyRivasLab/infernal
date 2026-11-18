/* cm_p7_band_chk.c
 *
 * Checkpointed banded CP9 P7B Forward/Backward + band reduction, double
 * precision (brief 26_0430-150/26_0430-153/26_0430-154). Produces cp9b bands with
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
 * Design (briefs 26_0430-144/26_0430-146):
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
 *   precision by brief 26_0430-154) is a verbatim transcription of the
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
 *   (a) cp9_ForwardP7BF has the brief 26_0430-134 BM-coverage supplementary pass
 *       (begin + EL-from-into-M for the full row band minus the match range);
 *       it writes mmx/elmx cells the posterior reads.
 *   (b) the erow EL-from-M+1 accumulation is NOT gated on INBAND(i,M).
 *       (erow is dead for band derivation, but reproduced for fidelity.)
 *
 * BRIEF 26_0821-207: the row kernels no longer follow the originals for the
 * local-end (EL) transitions. The originals' Forward and Backward used
 * different EL transition sets under a narrow band (see cp9_chk_fwd_rowF()'s
 * supplementary pass and the EL block at the top of cp9_chk_bwd_rowF()), so
 * their totals disagreed; both directions now use every EL transition whose
 * two cells are in band. The non-checkpointed float originals in cm_p7_band.c
 * have no live caller and are left as they were.
 *
 * Function/struct suffix F. Substitutions vs the originals:
 *   mmx[i]->mc, mmx[i-1]->mp, mmx[i+1]->mn (and i/d/el likewise).
 *
 * BRIEF 26_0430-154: this FLOAT-suffix family is now DOUBLE precision (cells +
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
 * matrix storage stays double (brief 26_0430-154), only the per-call transcendental is
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

__attribute__((constructor)) static void
cp9_chk_dlogsum_lut_init(void)
{
  int i;
  /* BASE 2, matching ILogsum's sreLOG2(1.+sreEXP2(-i/INTSCALE)) exactly. */
  for(i = 0; i < CP9_DLOGSUM_TBLN; i++) cp9_dlogsum_tbl[i] = log2(1. + exp2((double) -i / CP9_DLOGSUM_SCALE));
}

/* Double-precision log-sum for the checkpointed double-trunc CP9 F/B kernels
 * (brief 26_0430-153/26_0430-154/26_0430-193). -inf-guarded. LUT.
 *
 * BASE 2 (briefs 26_0821-024/26_0821-025/26_0821-026). Every value this combines is a
 * Scorify()d score, i.e. BITS -- Infernal has used bits, not nats, since 2007
 * (see the EPN note in logsum.c). This was transcribed from the int kernels
 * with ILogsum -> p7_FLogsum, which silently changed the base while the values
 * stayed in bits; it under-added every sum by up to 1-ln2 = 0.307 bits,
 * one-directionally, compounding row by row. */
#ifdef CM_P7_BAND_CHK_TESTDRIVE
/* unit tests only: TRUE makes cp9_chk_dlogsum() exact, so the kernels can be
 * compared cell by cell with an exact reference (brief 26_0821-207). */
static int cp9_chk_dlogsum_exact = FALSE;
#endif
static inline double
cp9_chk_dlogsum(double a, double b)
{
  if(a == -eslINFINITY) return b;
  if(b == -eslINFINITY) return a;
  const double max = (a > b) ? a : b;
  const double min = (a > b) ? b : a;
#ifdef CM_P7_BAND_CHK_TESTDRIVE
  if(cp9_chk_dlogsum_exact) return max + log2(1. + exp2(min - max));
#endif
  return ((max - min) >= CP9_DLOGSUM_CUTOFF) ? max
         : max + cp9_dlogsum_tbl[(int)((max - min) * CP9_DLOGSUM_SCALE)];
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
 * the brief 26_0430-134 BM-coverage supplementary pass. Cur planes mc/ic/dc/ec (offset
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

  /* Supplementary pass over the rest of the row band: [max(1,kmin[i]),kmax[i]]
   * MINUS the [kn,kx] match range. These cells have no in-band diagonal
   * predecessor, but M_k(i) can still be reached by a begin (B -> M_k, only
   * when row i-1 holds column 0) or by an EL return (EL_c(i-1) -> M_k), and
   * EL_k(i) by M_k(i) -> EL_k(i) or by the self-loop EL_k(i-1) -> EL_k(i).
   *
   * brief 26_0430-134 added this pass for the begin term and ran it only when
   * row i-1 holds column 0, so on every other row the EL transitions into
   * these cells were dropped while the Backward kept them. The Forward total
   * then fell short of the Backward total, posteriors (normalized by the
   * Backward total) lost occupancy after a local end, and the band predictor
   * opened the truncation window and blew up the CM bands (brief
   * 26_0821-206). The pass now runs on every row; only its begin term is
   * gated (brief 26_0821-207). */
  {
    int k_lo   = ESL_MAX(1, kmin[i]);
    int k_hi   = kmax[i];
    int do_beg = (INBAND(i-1, 0) && mp[0] != -eslINFINITY) ? TRUE : FALSE;
    for (k = k_lo; k <= k_hi; k++) {
      if(k >= kn && k <= kx) continue;   /* already filled by the match loop */
      kpcur = k - kmin[i];

      sc = do_beg ? mp[0] + Scorify(CP9TSC(cp9O_BM,k)) : -eslINFINITY;  /* begin; diagonal predecessor out of band here */

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
   * (does not require INBAND(i,M)). erow is dead for band derivation: the sum
   * is only returned (ret_erow != NULL) on row L, so only accumulate it there.
   * Walking all el_from_ct[M+1] sources on every row costs O(L*M) per pass on
   * models where nearly every node can local-end (e.g. pure-MATL genome CMs),
   * independent of the band width (brief 26_0821-107). */
  if ((cp9->flags & CPLAN9_EL) && ret_erow != NULL) {
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

  /* Stop at k=1, as the dense cp9ck_bck_rowL() does: M_0 has no emission, so
   * msc[x][0] is never assigned, and the k=0 cells are all set by the
   * INBAND(i, 0) block below (brief 26_0821-134). */
  for (k = kx; k >= ESL_MAX(kmin[i], 1); k--, kpcur--)
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

  /* EL_c(i), complete, before any cell of row i reads it. EL cells depend
   * only on row i+1, so they are filled here in one pass over every EL
   * transition whose two cells are in band, the mirror of the Forward's
   * transition set:
   *   EL_k(i) -> EL_k(i+1)  self-loop, k in the bands of rows i and i+1;
   *   EL_c(i) -> M_k(i+1)   return (no transition score), every k in row
   *                         i+1's band whose source c is in row i's band.
   * They used to be accumulated inside the match loop below, over its window
   * [max(kmin[i],kmin[i+1]-1), min(kmax[i],kmax[i+1]-1)], plus a separate
   * k == M block; that dropped the self-loop at in-band k outside the window
   * and the returns into M_k(i+1) for k outside it, a different set from the
   * Forward's, so the Forward and Backward totals disagreed (brief
   * 26_0821-206). Row i's M_k(i) reads EL_k(i) in the DELETIONS loops below
   * (brief 26_0821-207). */
  kpcur = 0;
  for (k = kmin[i]; k <= kmax[i]; k++, kpcur++) ec[kpcur] = -eslINFINITY;
  if(cp9->flags & CPLAN9_EL) {
    kn = ESL_MAX(1, ESL_MAX(kmin[i], kmin[i+1]));
    kx = ESL_MIN(kmax[i], kmax[i+1]);
    for (k = kn; k <= kx; k++)
      if(cp9->has_el[k]) ec[k-kmin[i]] = en[k-kmin[i+1]] + Scorify(cp9->el_selfsc);
    for (k = ESL_MAX(1, kmin[i+1]); k <= kmax[i+1]; k++) {
      kpprv = k - kmin[i+1];
      if(mn[kpprv] == -eslINFINITY) continue;
      for(c = 0; c < cp9->el_from_ct[k]; c++)
        if(INBAND(i, cp9->el_from_idx[k][c])) {
          kpcur_el = cp9->el_from_idx[k][c] - kmin[i];
          ec[kpcur_el] = cp9_chk_dlogsum(ec[kpcur_el], mn[kpprv]);
        }
    }
  }

  if(INBAND(i, M)) {
    kpcur = M-kmin[i];
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
  }

  /* MATCH: *_k <- M_k+1 */
  kn = ESL_MAX(kmin[i], kmin[i+1]-1);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]-1);

  for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mc[kpcur] = ic[kpcur] = dc[kpcur] = -eslINFINITY;
  for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mc[kpcur] = ic[kpcur] = dc[kpcur] = -eslINFINITY;

  kpcur = kx - kmin[i];
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    {
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

  /* D_k(i==0) <- I_k(i==1), k < M: D_M was set from I_M above, and adding
   * it again here counted that transition twice (+1 bit on D_M(0); brief
   * 26_0821-207). */
  kn = ESL_MAX(kmin[i], kmin[i+1]);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]);
  kx = ESL_MIN(kx, M-1);
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
  int i, j, w, zw;

  ESL_ALLOC(m0, sizeof(double)*(M+1)); ESL_ALLOC(i0, sizeof(double)*(M+1));
  ESL_ALLOC(d0, sizeof(double)*(M+1)); ESL_ALLOC(e0, sizeof(double)*(M+1));
  ESL_ALLOC(m1, sizeof(double)*(M+1)); ESL_ALLOC(i1, sizeof(double)*(M+1));
  ESL_ALLOC(d1, sizeof(double)*(M+1)); ESL_ALLOC(e1, sizeof(double)*(M+1));

  mc=m0; ic=i0; dc=d0; ec=e0;
  /* zero only the band + the single +1 edge cell each row kernel reads
   * (see brief 26_0821-084's audit); was M+1 (full model width). */
  zw = ESL_MIN(M+1, kmax[0]-kmin[0]+2);
  memset(mc,0,sizeof(double)*zw); memset(ic,0,sizeof(double)*zw); memset(dc,0,sizeof(double)*zw); memset(ec,0,sizeof(double)*zw);
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
    zw = ESL_MIN(M+1, kmax[i]-kmin[i]+2);
    memset(mc,0,sizeof(double)*zw); memset(ic,0,sizeof(double)*zw); memset(dc,0,sizeof(double)*zw); memset(ec,0,sizeof(double)*zw);
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
  int i, j, w, zw;

  ESL_ALLOC(m0, sizeof(double)*(M+1)); ESL_ALLOC(i0, sizeof(double)*(M+1));
  ESL_ALLOC(d0, sizeof(double)*(M+1)); ESL_ALLOC(e0, sizeof(double)*(M+1));
  ESL_ALLOC(m1, sizeof(double)*(M+1)); ESL_ALLOC(i1, sizeof(double)*(M+1));
  ESL_ALLOC(d1, sizeof(double)*(M+1)); ESL_ALLOC(e1, sizeof(double)*(M+1));

  j = s->nbnd - 1;

  /* row L. Buffer = i%2. */
  if(L & 1) { mc=m1; ic=i1; dc=d1; ec=e1; } else { mc=m0; ic=i0; dc=d0; ec=e0; }
  /* zero only the band + the single +1 edge cell each row kernel reads
   * (see brief 26_0821-084's audit); was M+1 (full model width). */
  zw = ESL_MIN(M+1, kmax[L]-kmin[L]+2);
  memset(mc,0,sizeof(double)*zw); memset(ic,0,sizeof(double)*zw); memset(dc,0,sizeof(double)*zw); memset(ec,0,sizeof(double)*zw);
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
    zw = ESL_MIN(M+1, kmax[i]-kmin[i]+2);
    memset(mc,0,sizeof(double)*zw); memset(ic,0,sizeof(double)*zw); memset(dc,0,sizeof(double)*zw); memset(ec,0,sizeof(double)*zw);
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
  zw = ESL_MIN(M+1, kmax[0]-kmin[0]+2);
  memset(mc,0,sizeof(double)*zw); memset(ic,0,sizeof(double)*zw); memset(dc,0,sizeof(double)*zw); memset(ec,0,sizeof(double)*zw);
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
 * total score of -inf or, more often, about IMPOSSIBLE (-1e36: see the test
 * below), every posterior below becomes meaningless, no
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
 * alignment every sequence must produce a parse, so this returns eslENORESULT.
 *
 * It is an error only if no caller-level fallback recovers.  Since brief
 * 26_0821-096, cmalign usually does: a genome-scale no-parse band is retried
 * once with --p7ibv-mem bands, and below genome scale the non-p7-banded CP9
 * re-derivation runs.  This leaf cannot tell which will happen, so it
 * prints nothing; when no fallback recovers, the caller fails with an
 * Error: message built from errbuf (brief 26_0821-097).
 */
static int
cp9_chk_noparse_check(char *errbuf, const char *where, int L, int M,
                      int *kmin, int *kmax, double fsc, double bsc)
{
  /* An impossible total is usually NOT non-finite here. The CP9 transition
   * scores reach these kernels through Scorify(), which maps -INFTY to the
   * finite IMPOSSIBLE (-1e36, infernal.h), so a band with no complete parse
   * typically yields a total of about IMPOSSIBLE rather than -inf (briefs
   * 26_0821-099, 26_0821-100). infernal.h guarantees that any sum involving
   * IMPOSSIBLE is below IMPROBABLE, so test against that as well. */
  if(isfinite(fsc) && isfinite(bsc) && NOT_IMPROBABLE(fsc) && NOT_IMPROBABLE(bsc)) return eslOK;

  ESL_FAIL(eslENORESULT, errbuf,
           "no-parse p7 band: no complete parse under this band (CP9 F/B total impossible)");
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
  return "default p7 Viterbi-trace";
}

/* brief 26_0430-310: fold the deriver name into errbuf, so the error the caller
 * prints names it.
 */
static void
cp9_chk_noparse_report(CM_t *cm, char *errbuf)
{
  const char *deriver = cp9_chk_band_deriver_name(cm);

  snprintf(errbuf, eslERRBUFSIZE,
           "no-parse p7 band from %s: no complete parse under this band "
           "(CP9 F/B total impossible)", deriver);
}

/* The double checkpointed band reduction: produces cp9b pn_min/pn_max bands and
 * the per-node pocc_arr (match+delete occupancy, streamed in the MIN sweep)
 * exactly as cp9_FB2HMMBandsP7BF + cp9_PredictStartAndEndPositionsP7BF's pocc
 * loop, but with O(sqrt(L)*avg_bw) memory.
 *
 * pocc_arr: caller-allocated [0..M]. Filled here: pocc_arr[k] = sum over
 * i=0..L of (exp2(pmx->mmx[i][kp]) + exp2(pmx->dmx[i][kp])) for nodes k with a
 * band set, else -1.0 (matches the original's skip+sentinel behavior).
 */
int
cp9_FB2HMMBandsP7BF_chk(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9Bands_t *cp9b,
                        int L, int M, double p_thresh, int *kmin, int *kmax,
                        int debug_level, int do_pnmono, int do_pnmono_print,
                        double *pocc_arr)
{
  int status;
  /* BASE 2 (26_0821-026): compared against mass_* values accumulated by
   * cp9_chk_dlogsum() from Scorify()d BIT scores. The int arm builds the same
   * threshold as Prob2Score(((1.-p_thresh)/2.), 1.) = INTSCALE*sreLOG2(p),
   * i.e. base 2 (hmmband.c:673, and identically in 1.1.5 at :579/:906). */
  double thresh = log2((1. - p_thresh) / 2.);
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


  if((status = cp9chkF_FwdFill(s, hmm, dsq, kmin, kmax, errbuf)) != eslOK) goto ERROR;
  if((status = cp9chkF_BwdFill(s, hmm, dsq, kmin, kmax, &sc, errbuf)) != eslOK) goto ERROR;
  /* brief 26_0430-310 (M3): no-parse band guard, BEFORE sc is used as the
   * posterior normalizer in cp9segF_PostRow() below. See cp9_chk_noparse_check(). */
  if((status = cp9_chk_noparse_check(errbuf, "cp9_FB2HMMBandsP7BF_chk", L, M,
                                     kmin, kmax, s->fsc, sc)) != eslOK) goto ERROR;
  cp9b->hmm_sc = sc; /* for centering CM alignment DP scores (DispatchSqAlignment()) */
  /* brief 26_0430-154 diag: double checkpointed CP9 F/B totals (cf. 153 ref fwd/bwd≈7993.144,
   * gap≈+0.00003; float ckpt path gave gap≈-1.687 at norovirus, ≈-35/-45 at HSV/MPXV). */
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
      /* pocc streaming: sum exp2(pm)+exp2(pd) over k>=1 in band, in TWO separate
       * adds. BASE 2 (26_0821-029): g->pm/pd are Scorify()d posteriors, i.e.
       * BITS, so the probability is 2^x, not e^x. See the commit message: this
       * correction is known to DEGRADE some glocal cells, which is evidence of
       * a second, non-independent defect -- not a reason to keep this one.
       * adds (matches cp9_PredictStartAndEndPositionsP7BF's double add order). */
      kkp = ESL_MAX(1, kmin[i]) - kmin[i];
      for(kk = ESL_MAX(1, kmin[i]); kk <= kmax[i]; kk++, kkp++) {
        pocc_arr[kk] += exp2(g->pm[kkp]);
        pocc_arr[kk] += exp2(g->pd[kkp]);
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
  double *thresh = NULL;                       /* thresh[t] = log2((1-p_thresh[t])/2) */
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
  /* BASE 2 -- see the note in cp9_FB2HMMBandsP7BF_chk (26_0821-026). */
  for(t = 0; t < NS; t++) thresh[t] = log2((1. - p_thresh[t]) / 2.);

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
  cp9b->hmm_sc = sc; /* for centering CM alignment DP scores (DispatchSqAlignment()) */
  /* brief 26_0430-193: same F/B total dump as P154_FBDUMP (single-threshold path,
   * line ~2243), replicated here in the MULTI (tau-ratchet) band-extraction path so
   * genome-scale runs that take this path still report the F/B gap for the LUT-vs-exact
   * precision comparison. Same env var, so only one dump fires per sequence. */
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
        pocc_raw[kk] += exp2(g->pm[kkp]);
        pocc_raw[kk] += exp2(g->pd[kkp]);
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
 * i0..j0, predicts sp/ep (+ brief 26_0430-149 glocal floor) & marginal candidates
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
     * (confirmed by direct trace on the brief 26_0430-201 dengue LC436672.1/
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
    cp9b->im_widened = FALSE;
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
          /* the last node's insert maps to ROOT_IR; record that it got a band
           * with no posterior evidence so the --mxsize fallback can withdraw
           * it (cp9_EndHullRetryBands(), brief 26_0821-103) */
          if(k == hmm_M) cp9b->im_widened = TRUE;
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
     * the brief 26_0430-149 sp1/ep1 floor above forces Jvalid[v] = TRUE for essentially
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
 *     extra memory/compute (critical: keeps the brief 26_0430-165 genome capstone path
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


/*****************************************************************
 * Unit tests (brief 26_0821-207)
 *   make cm_p7_band_chk_utest; ./cm_p7_band_chk_utest <cmfile>
 *
 * The CP9 is built from a real CM (with local ends), so the tests
 * need a CM file; the testsuite passes one.
 *****************************************************************/
#ifdef CM_P7_BAND_CHK_TESTDRIVE

#include "esl_getopts.h"
#include "esl_random.h"

/* Exact reference: the unbanded CP9 Forward and Backward recursions of
 * cp9_Forward()/cp9_Backward() (cp9_dp.c; alignment mode, no scanning) in
 * double with an exact log-sum, with every cell outside kmin[i]..kmax[i] set
 * to -inf as soon as it is computed. Its in-band cells therefore sum exactly
 * the parses that stay inside the band, which is what the banded kernels must
 * compute. With a full band it is the unbanded F/B. One deliberate departure
 * from cp9_Backward(): EL_M keeps its self-loop there too (cp9_Backward()
 * adds el_selfsc to a cell it has just set to -INFTY), so that Forward and
 * Backward use the same transitions. */
typedef struct {
  int     L, M;
  double *fm, *fi, *fd, *fe;
  double *bm, *bi, *bd, *be;
  double  fsc, bsc;
} refFB_t;

#define RX(a,i,k)   ((a)[(int64_t)(i) * (M+1) + (k)])
#define T(s,k)      ((double) Scorify(CP9TSC((s),(k))))
#define REF_MSC(i,k)    ((double) Scorify(cp9->msc[dsq[(i)]][(k)]))
#define REF_ISC(i,k)    ((double) Scorify(cp9->isc[dsq[(i)]][(k)]))
#define MASK(i,k,v) (INBAND((i),(k)) ? (v) : -eslINFINITY)

static double
ref_ls(double a, double b)
{
  if(a == -eslINFINITY) return b;
  if(b == -eslINFINITY) return a;
  return (a > b) ? a + log2(1. + exp2(b - a)) : b + log2(1. + exp2(a - b));
}
static double ref_ls3(double a, double b, double c) { return ref_ls(ref_ls(a, b), c); }

static void
ref_Destroy(refFB_t *r)
{
  if(r == NULL) return;
  free(r->fm); free(r->fi); free(r->fd); free(r->fe);
  free(r->bm); free(r->bi); free(r->bd); free(r->be);
  free(r);
}

static refFB_t *
ref_Create(CP9_t *cp9, ESL_DSQ *dsq, int L, int *kmin, int *kmax)
{
  int const *tsc = cp9->otsc;
  int        M   = cp9->M;
  int        el  = (cp9->flags & CPLAN9_EL) ? TRUE : FALSE;
  int64_t    n   = (int64_t) (L+1) * (M+1);
  refFB_t   *r   = malloc(sizeof(refFB_t));
  double     sc;
  int        i, k, c;

  r->L = L; r->M = M;
  r->fm = malloc(sizeof(double) * n); r->fi = malloc(sizeof(double) * n);
  r->fd = malloc(sizeof(double) * n); r->fe = malloc(sizeof(double) * n);
  r->bm = malloc(sizeof(double) * n); r->bi = malloc(sizeof(double) * n);
  r->bd = malloc(sizeof(double) * n); r->be = malloc(sizeof(double) * n);

  /* Forward */
  RX(r->fm,0,0) = MASK(0, 0, 0.);
  RX(r->fi,0,0) = RX(r->fd,0,0) = RX(r->fe,0,0) = -eslINFINITY;
  for(k = 1; k <= M; k++) {
    RX(r->fm,0,k) = RX(r->fi,0,k) = RX(r->fe,0,k) = -eslINFINITY;
    RX(r->fd,0,k) = MASK(0, k, ref_ls3(RX(r->fm,0,k-1) + T(cp9O_MD,k-1), RX(r->fi,0,k-1) + T(cp9O_ID,k-1), RX(r->fd,0,k-1) + T(cp9O_DD,k-1)));
  }
  for(i = 1; i <= L; i++) {
    RX(r->fm,i,0) = RX(r->fd,i,0) = RX(r->fe,i,0) = -eslINFINITY;
    RX(r->fi,i,0) = MASK(i, 0, ref_ls3(RX(r->fm,i-1,0) + T(cp9O_MI,0), RX(r->fi,i-1,0) + T(cp9O_II,0), RX(r->fd,i-1,0) + T(cp9O_DI,0)) + REF_ISC(i,0));
    for(k = 1; k <= M; k++) {
      sc = ref_ls(ref_ls3(RX(r->fm,i-1,k-1) + T(cp9O_MM,k-1), RX(r->fi,i-1,k-1) + T(cp9O_IM,k-1), RX(r->fd,i-1,k-1) + T(cp9O_DM,k-1)),
                  RX(r->fm,i-1,0) + T(cp9O_BM,k));
      if(el) for(c = 0; c < cp9->el_from_ct[k]; c++) sc = ref_ls(sc, RX(r->fe,i-1,cp9->el_from_idx[k][c]));
      RX(r->fm,i,k) = MASK(i, k, sc + REF_MSC(i,k));
      RX(r->fi,i,k) = MASK(i, k, ref_ls3(RX(r->fm,i-1,k) + T(cp9O_MI,k), RX(r->fi,i-1,k) + T(cp9O_II,k), RX(r->fd,i-1,k) + T(cp9O_DI,k)) + REF_ISC(i,k));
      RX(r->fd,i,k) = MASK(i, k, ref_ls3(RX(r->fm,i,k-1) + T(cp9O_MD,k-1), RX(r->fi,i,k-1) + T(cp9O_ID,k-1), RX(r->fd,i,k-1) + T(cp9O_DD,k-1)));
      RX(r->fe,i,k) = (el && cp9->has_el[k]) ? MASK(i, k, ref_ls(RX(r->fm,i,k) + T(cp9O_MEL,k), RX(r->fe,i-1,k) + Scorify(cp9->el_selfsc))) : -eslINFINITY;
    }
  }
  r->fsc = ref_ls(RX(r->fd,L,M) + T(cp9O_DM,M), RX(r->fi,L,M) + T(cp9O_IM,M));
  for(k = 1; k <= M; k++) r->fsc = ref_ls(r->fsc, RX(r->fm,L,k) + T(cp9O_ME,k));
  if(el) for(c = 0; c < cp9->el_from_ct[M+1]; c++) r->fsc = ref_ls(r->fsc, RX(r->fe,L,cp9->el_from_idx[M+1][c]));

  /* Backward, row L */
  for(k = 0; k <= M; k++) RX(r->be,L,k) = -eslINFINITY;
  if(el) for(c = 0; c < cp9->el_from_ct[M+1]; c++) { k = cp9->el_from_idx[M+1][c]; RX(r->be,L,k) = MASK(L, k, 0.); }
  RX(r->bm,L,M) = MASK(L, M, ref_ls(RX(r->be,L,M) + T(cp9O_MEL,M), T(cp9O_ME,M)) + REF_MSC(L,M));
  RX(r->bi,L,M) = MASK(L, M, T(cp9O_IM,M) + REF_ISC(L,M));
  RX(r->bd,L,M) = MASK(L, M, T(cp9O_DM,M));
  for(k = M-1; k >= 1; k--) {
    sc = ref_ls(T(cp9O_ME,k), RX(r->bd,L,k+1) + T(cp9O_MD,k));
    if(el) sc = ref_ls(sc, RX(r->be,L,k) + T(cp9O_MEL,k));
    RX(r->bm,L,k) = MASK(L, k, sc + REF_MSC(L,k));
    RX(r->bi,L,k) = MASK(L, k, RX(r->bd,L,k+1) + T(cp9O_ID,k) + REF_ISC(L,k));
    RX(r->bd,L,k) = MASK(L, k, RX(r->bd,L,k+1) + T(cp9O_DD,k));
  }
  RX(r->bm,L,0) = MASK(L, 0, RX(r->bd,L,1) + T(cp9O_MD,0));
  RX(r->bi,L,0) = MASK(L, 0, RX(r->bd,L,1) + T(cp9O_ID,0) + REF_ISC(L,0));
  RX(r->bd,L,0) = RX(r->be,L,0) = -eslINFINITY;

  /* Backward, rows L-1..1 */
  for(i = L-1; i >= 1; i--) {
    for(k = 0; k <= M; k++) RX(r->be,i,k) = -eslINFINITY;
    if(el) {
      for(k = 1; k <= M; k++) if(cp9->has_el[k]) RX(r->be,i,k) = RX(r->be,i+1,k) + Scorify(cp9->el_selfsc);
      for(k = 1; k <= M; k++)
        for(c = 0; c < cp9->el_from_ct[k]; c++) RX(r->be,i,cp9->el_from_idx[k][c]) = ref_ls(RX(r->be,i,cp9->el_from_idx[k][c]), RX(r->bm,i+1,k));
      for(k = 1; k <= M; k++) RX(r->be,i,k) = MASK(i, k, RX(r->be,i,k));
    }
    sc = RX(r->bi,i+1,M) + T(cp9O_MI,M);
    if(el && cp9->has_el[M]) sc = ref_ls(sc, RX(r->be,i,M) + T(cp9O_MEL,M));
    RX(r->bm,i,M) = MASK(i, M, sc + REF_MSC(i,M));
    RX(r->bi,i,M) = MASK(i, M, RX(r->bi,i+1,M) + T(cp9O_II,M) + REF_ISC(i,M));
    RX(r->bd,i,M) = MASK(i, M, RX(r->bi,i+1,M) + T(cp9O_DI,M));
    for(k = M-1; k >= 1; k--) {
      sc = ref_ls3(RX(r->bm,i+1,k+1) + T(cp9O_MM,k), RX(r->bi,i+1,k) + T(cp9O_MI,k), RX(r->bd,i,k+1) + T(cp9O_MD,k));
      if(el && cp9->has_el[k]) sc = ref_ls(sc, RX(r->be,i,k) + T(cp9O_MEL,k));
      RX(r->bm,i,k) = MASK(i, k, sc + REF_MSC(i,k));
      RX(r->bi,i,k) = MASK(i, k, ref_ls3(RX(r->bm,i+1,k+1) + T(cp9O_IM,k), RX(r->bi,i+1,k) + T(cp9O_II,k), RX(r->bd,i,k+1) + T(cp9O_ID,k)) + REF_ISC(i,k));
      RX(r->bd,i,k) = MASK(i, k, ref_ls3(RX(r->bm,i+1,k+1) + T(cp9O_DM,k), RX(r->bi,i+1,k) + T(cp9O_DI,k), RX(r->bd,i,k+1) + T(cp9O_DD,k)));
    }
    RX(r->bi,i,0) = MASK(i, 0, ref_ls3(RX(r->bm,i+1,1) + T(cp9O_IM,0), RX(r->bi,i+1,0) + T(cp9O_II,0), RX(r->bd,i,1) + T(cp9O_ID,0)) + REF_ISC(i,0));
    RX(r->bd,i,0) = RX(r->be,i,0) = -eslINFINITY;
    sc = ref_ls(RX(r->bi,i+1,0) + T(cp9O_MI,0), RX(r->bd,i,1) + T(cp9O_MD,0));
    for(k = M; k >= 1; k--) sc = ref_ls(sc, RX(r->bm,i+1,k) + T(cp9O_BM,k));
    RX(r->bm,i,0) = MASK(i, 0, sc);
  }

  /* Backward, row 0 */
  for(k = 0; k <= M; k++) RX(r->bm,0,k) = RX(r->bi,0,k) = RX(r->be,0,k) = -eslINFINITY;
  RX(r->bd,0,M) = MASK(0, M, RX(r->bi,1,M) + T(cp9O_DI,M));
  for(k = M-1; k >= 1; k--)
    RX(r->bd,0,k) = MASK(0, k, ref_ls3(RX(r->bm,1,k+1) + T(cp9O_DM,k), RX(r->bi,1,k) + T(cp9O_DI,k), RX(r->bd,0,k+1) + T(cp9O_DD,k)));
  RX(r->bd,0,0) = -eslINFINITY;
  sc = ref_ls(RX(r->bi,1,0) + T(cp9O_MI,0), RX(r->bd,0,1) + T(cp9O_MD,0));
  for(k = M; k >= 1; k--) sc = ref_ls(sc, RX(r->bm,1,k) + T(cp9O_BM,k));
  RX(r->bm,0,0) = MASK(0, 0, sc);
  r->bsc = RX(r->bm,0,0);
  return r;
}

/* Two cell values agree: both impossible (-inf, or below -1e5 bits, i.e.
 * reached only through an -INFTY score), or within tol. */
static int
cells_agree(double a, double b, double tol)
{
  if(a < -1e5 && b < -1e5) return TRUE;
  return (fabs(a - b) <= tol) ? TRUE : FALSE;
}

/* Run the checkpointed banded F/B on <dsq>, <kmin>/<kmax> and compare it with
 * the reference <r>: totals within <tol_tot>; with <do_cells>, every in-band
 * cell of all eight planes (rebuilt segment by segment from the checkpoints,
 * as the band derivation does) within <tol_cell>. Returns the number of
 * failed checks; *ret_maxd is the largest cell or total difference. */
static int
chk_vs_ref(CP9_t *cp9, ESL_DSQ *dsq, int L, int *kmin, int *kmax, refFB_t *r,
           int do_cells, double tol_tot, double tol_cell,
           double *ret_fsc, double *ret_bsc, double *ret_maxd, char *what, int be_verbose)
{
  char        errbuf[eslERRBUFSIZE];
  int         M = cp9->M;
  cp9chkF_t  *s = cp9chkF_Create(L, M, kmin, kmax, errbuf);
  cp9segF_t  *g = NULL;
  double      bsc, d, maxd = 0.;
  double     *kv[8], *rv[8];
  int         nfail = 0, j, rr, kp, k, p, a;

  if(s == NULL) esl_fatal("%s: cp9chkF_Create failed", what);
  if(cp9chkF_FwdFill(s, cp9, dsq, kmin, kmax, errbuf)       != eslOK) esl_fatal("%s: FwdFill failed", what);
  if(cp9chkF_BwdFill(s, cp9, dsq, kmin, kmax, &bsc, errbuf) != eslOK) esl_fatal("%s: BwdFill failed", what);
  if(! (isfinite(s->fsc) && NOT_IMPROBABLE(s->fsc))) nfail++;
  if(! (isfinite(bsc)    && NOT_IMPROBABLE(bsc)))    nfail++;
  d = fabs(s->fsc - r->fsc); if(! (d <= tol_tot)) nfail++; if(isfinite(d) && d > maxd) maxd = d; else if(! isfinite(d)) maxd = eslINFINITY;
  d = fabs(bsc    - r->bsc); if(! (d <= tol_tot)) nfail++; if(isfinite(d) && d > maxd) maxd = d; else if(! isfinite(d)) maxd = eslINFINITY;
  if(be_verbose && nfail) printf("  %s: F %.6f B %.6f ref F %.6f B %.6f\n", what, s->fsc, bsc, r->fsc, r->bsc);

  if(do_cells) {
    if((g = cp9segF_Create(s, kmin, kmax, errbuf)) == NULL) esl_fatal("%s: cp9segF_Create failed", what);
    for(j = 1; j < s->nbnd; j++) {
      cp9segF_Fill(g, s, j, cp9, dsq, kmin, kmax);
      a = s->bnd[j-1];
      for(rr = a; rr <= s->bnd[j]; rr++) {
        kv[0] = g->fmr[rr-a]; kv[1] = g->fir[rr-a]; kv[2] = g->fdr[rr-a]; kv[3] = g->fer[rr-a];
        kv[4] = g->bmr[rr-a]; kv[5] = g->bir[rr-a]; kv[6] = g->bdr[rr-a]; kv[7] = g->ber[rr-a];
        rv[0] = r->fm; rv[1] = r->fi; rv[2] = r->fd; rv[3] = r->fe;
        rv[4] = r->bm; rv[5] = r->bi; rv[6] = r->bd; rv[7] = r->be;
        for(kp = 0, k = kmin[rr]; k <= kmax[rr]; k++, kp++)
          for(p = 0; p < 8; p++)
            if(! cells_agree(kv[p][kp], RX(rv[p],rr,k), tol_cell)) {
              if(be_verbose && nfail < 5) printf("  %s: plane %d row %d col %d: kernel %.6f ref %.6f\n", what, p, rr, k, kv[p][kp], RX(rv[p],rr,k));
              nfail++;
              d = fabs(kv[p][kp] - RX(rv[p],rr,k)); maxd = isfinite(d) ? ESL_MAX(maxd, d) : eslINFINITY;
            }
      }
    }
    cp9segF_Destroy(g);
  }
  if(ret_fsc)  *ret_fsc  = s->fsc;
  if(ret_bsc)  *ret_bsc  = bsc;
  if(ret_maxd) *ret_maxd = maxd;
  cp9chkF_Destroy(s);
  return nfail;
}

/* A random parse with a local end, and a band around it.
 * The parse: begin into M_s, matches M_s..M_a, a local end EL_a for <r>
 * residues, then either EL_a -> E, or EL_a -> M_b and matches M_b..M_e
 * followed by a local exit. (a,b) is drawn from the CP9's EL table, so the
 * parse exists in the model. The band is the parse's column on each row,
 * widened by a random 0..W cells on each side (row 0 always starts at
 * column 0). W = 0 is the parse itself: the EL rows then hold one column,
 * the case where the band kernels' diagonal match window misses the EL
 * self-loop on every row. */
static int
sample_el_band(ESL_RANDOMNESS *rng, CP9_t *cp9, int W, int abcK, ESL_DSQ **ret_dsq, int **ret_kmin, int **ret_kmax)
{
  int      M = cp9->M;
  int      npair = 0, pick, a = -1, b = -1, s, e, r, L, i, k, c, nrow;
  int     *col;
  int     *kmin, *kmax;
  ESL_DSQ *dsq;

  /* (a,b) pairs with an EL state at a >= 1 */
  for(k = 1; k <= M+1; k++)
    for(c = 0; c < cp9->el_from_ct[k]; c++)
      if(cp9->el_from_idx[k][c] >= 1 && cp9->has_el[cp9->el_from_idx[k][c]]) npair++;
  if(npair == 0) return -1;
  pick = esl_rnd_Roll(rng, npair);
  for(k = 1; k <= M+1 && a < 0; k++)
    for(c = 0; c < cp9->el_from_ct[k]; c++)
      if(cp9->el_from_idx[k][c] >= 1 && cp9->has_el[cp9->el_from_idx[k][c]] && pick-- == 0) { a = cp9->el_from_idx[k][c]; b = k; break; }

  s = 1 + esl_rnd_Roll(rng, a);                   /* begin node, 1..a */
  r = 1 + esl_rnd_Roll(rng, 30);                  /* EL residues      */
  e = (b > M) ? -1 : (esl_rnd_Roll(rng, 2) ? M : b + esl_rnd_Roll(rng, M-b+1));
  L = (a-s+1) + r + ((e < 0) ? 0 : (e-b+1));

  col  = malloc(sizeof(int) * (L+1));
  kmin = malloc(sizeof(int) * (L+1));
  kmax = malloc(sizeof(int) * (L+1));
  dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  nrow = 0;
  col[nrow++] = 0;
  for(k = s; k <= a; k++)        col[nrow++] = k;
  for(i = 0; i < r; i++)         col[nrow++] = a;
  if(e > 0) for(k = b; k <= e; k++) col[nrow++] = k;
  ESL_DASSERT1((nrow == L+1));

  for(i = 0; i <= L; i++) {
    int wl = esl_rnd_Roll(rng, W+1);   /* drawn once: ESL_MIN/ESL_MAX evaluate an argument twice */
    int wr = esl_rnd_Roll(rng, W+1);
    kmin[i] = (i == 0) ? 0 : ESL_MAX(0, col[i] - wl);
    kmax[i] = ESL_MIN(M, col[i] + wr);
  }
  dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
  for(i = 1; i <= L; i++) dsq[i] = esl_rnd_Roll(rng, abcK);

  free(col);
  *ret_dsq = dsq; *ret_kmin = kmin; *ret_kmax = kmax;
  return L;
}

/* utest_narrow_band(): on random narrow bands around a parse with a local
 * end (sample_el_band()), the checkpointed banded F/B must
 *  (1) with an exact log-sum: match the exact banded reference cell by cell
 *      and in both totals, to 1e-6 bits (so Forward total = Backward total);
 *  (2) with the production log-sum table: give Forward and Backward totals
 *      that agree to 0.01 bits, and lie within 0.01 bits of the exact ones.
 * Every sampled parse is a valid parse of the model inside the band, so a
 * total of -inf (no parse found) is a failure. */
static void
utest_narrow_band(ESL_RANDOMNESS *rng, CP9_t *cp9, char *cp9name, int abcK, int N, int be_verbose)
{
  int      Wv[5] = { 0, 1, 2, 3, 6 };
  char     what[128];
  ESL_DSQ *dsq;
  int     *kmin, *kmax;
  refFB_t *r;
  double   fsc, bsc, maxd, worst_cell = 0., worst_fb = 0.;
  int      n, L, nfail = 0, nskip = 0, nrun = 0;

  for(n = 0; n < N; n++) {
    int W = Wv[n % 5];
    if((L = sample_el_band(rng, cp9, W, abcK, &dsq, &kmin, &kmax)) <= 0) esl_fatal("utest_narrow_band: %s has no local ends", cp9name);
    r = ref_Create(cp9, dsq, L, kmin, kmax);
    snprintf(what, 128, "utest_narrow_band %s trial %d (W %d, L %d)", cp9name, n, W, L);
    if(! (isfinite(r->fsc) && NOT_IMPROBABLE(r->fsc))) { nskip++; goto NEXT; }   /* parse not allowed by this configuration */
    if(fabs(r->fsc - r->bsc) > 1e-6) esl_fatal("%s: reference F %.6f != B %.6f", what, r->fsc, r->bsc);
    nrun++;

    cp9_chk_dlogsum_exact = TRUE;
    if(chk_vs_ref(cp9, dsq, L, kmin, kmax, r, TRUE, 1e-6, 1e-6, NULL, NULL, &maxd, what, be_verbose)) nfail++;
    worst_cell = ESL_MAX(worst_cell, maxd);
    cp9_chk_dlogsum_exact = FALSE;
    chk_vs_ref(cp9, dsq, L, kmin, kmax, r, FALSE, 0.01, 0., &fsc, &bsc, &maxd, what, FALSE);
    if(! (fabs(fsc - bsc) <= 0.01 && fabs(fsc - r->fsc) <= 0.01)) {
      if(be_verbose) printf("  %s: table log-sum F %.6f B %.6f exact %.6f\n", what, fsc, bsc, r->fsc);
      nfail++;
    }
    worst_fb = isfinite(fsc - bsc) ? ESL_MAX(worst_fb, fabs(fsc - bsc)) : eslINFINITY;
  NEXT:
    ref_Destroy(r); free(dsq); free(kmin); free(kmax);
  }
  if(be_verbose) printf("utest_narrow_band: %s: %d trials run, %d skipped, %d failed; max exact cell/total diff %.3g bits, max table-logsum |F-B| %.3g bits\n",
                        cp9name, nrun, nskip, nfail, worst_cell, worst_fb);
  if(nrun < N/2) esl_fatal("utest_narrow_band: %s: only %d of %d sampled parses are allowed by the model", cp9name, nrun, N);
  if(nfail)      esl_fatal("utest_narrow_band: %s: %d of %d trials failed (run with -v)", cp9name, nfail, nrun);
}

/* utest_full_band(): with every cell in band, the checkpointed banded F/B
 * must equal the unbanded F/B: cell by cell and in both totals to 1e-6 bits
 * against the exact reference (exact log-sum), and in both totals to 0.01
 * bits with the production log-sum table. The reference itself is checked
 * against the production unbanded int cp9_Forward()/cp9_Backward() totals,
 * to 2e-4 bits per row and column (those accumulate integer log-sum error,
 * ~5e-5 bits per row and column observed). */
static void
utest_full_band(ESL_RANDOMNESS *rng, CP9_t *cp9, char *cp9name, int abcK, int N, int be_verbose)
{
  char     errbuf[eslERRBUFSIZE], what[128];
  int      M = cp9->M;
  CP9_MX  *mx = CreateCP9Matrix(1, M);
  ESL_DSQ *dsq;
  int     *kmin, *kmax;
  refFB_t *r;
  float    ifsc, ibsc;
  double   fsc, bsc, maxd, worst = 0., worst_int = 0.;
  int      n, i, L, nfail = 0;

  for(n = 0; n < N; n++) {
    L    = M/2 + esl_rnd_Roll(rng, 2*M);
    dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
    kmin = malloc(sizeof(int) * (L+1));
    kmax = malloc(sizeof(int) * (L+1));
    dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
    for(i = 1; i <= L; i++) dsq[i] = esl_rnd_Roll(rng, abcK);
    for(i = 0; i <= L; i++) { kmin[i] = 0; kmax[i] = M; }
    r = ref_Create(cp9, dsq, L, kmin, kmax);
    snprintf(what, 128, "utest_full_band %s trial %d (L %d)", cp9name, n, L);
    if(! (isfinite(r->fsc) && NOT_IMPROBABLE(r->fsc))) esl_fatal("%s: reference has no parse", what);
    if(fabs(r->fsc - r->bsc) > 1e-6)                    esl_fatal("%s: reference F %.6f != B %.6f", what, r->fsc, r->bsc);

    if(cp9_Forward (cp9, errbuf, mx, dsq, 1, L, FALSE, TRUE, FALSE, NULL, NULL, &ifsc) != eslOK) esl_fatal("%s: cp9_Forward failed", what);
    if(cp9_Backward(cp9, errbuf, mx, dsq, 1, L, FALSE, TRUE, FALSE, NULL, NULL, &ibsc) != eslOK) esl_fatal("%s: cp9_Backward failed", what);
    worst_int = ESL_MAX(worst_int, ESL_MAX(fabs(ifsc - r->fsc), fabs(ibsc - r->bsc)));
    if(fabs(ifsc - r->fsc) > 2e-4 * (L+M)) esl_fatal("%s: int cp9_Forward %.4f vs reference %.4f", what, ifsc, r->fsc);
    if(! cp9->has_el[M] && fabs(ibsc - r->bsc) > 2e-4 * (L+M)) esl_fatal("%s: int cp9_Backward %.4f vs reference %.4f", what, ibsc, r->bsc);

    cp9_chk_dlogsum_exact = TRUE;
    if(chk_vs_ref(cp9, dsq, L, kmin, kmax, r, TRUE, 1e-6, 1e-6, NULL, NULL, &maxd, what, be_verbose)) nfail++;
    worst = ESL_MAX(worst, maxd);
    cp9_chk_dlogsum_exact = FALSE;
    chk_vs_ref(cp9, dsq, L, kmin, kmax, r, FALSE, 0.01, 0., &fsc, &bsc, &maxd, what, FALSE);
    if(! (fabs(fsc - bsc) <= 0.01 && fabs(fsc - r->fsc) <= 0.01)) {
      if(be_verbose) printf("  %s: table log-sum F %.6f B %.6f exact %.6f\n", what, fsc, bsc, r->fsc);
      nfail++;
    }
    ref_Destroy(r); free(dsq); free(kmin); free(kmax);
  }
  FreeCP9Matrix(mx);
  if(be_verbose) printf("utest_full_band: %s: %d trials, %d failed; max exact cell/total diff %.3g bits; reference vs int unbanded totals %.3g bits\n",
                        cp9name, N, nfail, worst, worst_int);
  if(nfail) esl_fatal("utest_full_band: %s: %d of %d trials failed (run with -v)", cp9name, nfail, N);
}

static ESL_OPTIONS options[] = {
  /* name     type          default  env  range toggles reqs incomp  help                                 docgroup*/
  { "-h",        eslARG_NONE, FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help and exit",           0 },
  { "--version", eslARG_NONE, FALSE, NULL, NULL,  NULL,  NULL, NULL, "show version info and exit",         0 },
  { "-s",        eslARG_INT,    "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",      0 },
  { "-N",        eslARG_INT,  "100", NULL, "n>0", NULL,  NULL, NULL, "number of random narrow bands",      0 },
  { "-v",        eslARG_NONE, FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "unit test driver for cm_p7_band_chk.c (checkpointed banded CP9 Forward/Backward)";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = cm_CreateDefaultApp("cm_p7_band_chk_utest", options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  int             N   = esl_opt_GetInteger(go, "-N");
  int             v   = esl_opt_GetBoolean(go, "-v");
  ESL_ALPHABET   *abc = NULL;
  CM_FILE        *cmfp = NULL;
  CM_t           *cm   = NULL;
  char            errbuf[eslERRBUFSIZE];

  init_ilogsum();   /* the int unbanded cp9_Forward()/cp9_Backward() need it */
  if(cm_file_Open(esl_opt_GetArg(go, 1), NULL, FALSE, &cmfp, errbuf) != eslOK) esl_fatal("failed to open CM file: %s", errbuf);
  if(cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) esl_fatal("failed to read CM: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  /* cmalign's default configuration: local, with CP9 local ends, truncated */
  cm->config_opts |= CM_CONFIG_LOCAL | CM_CONFIG_HMMLOCAL | CM_CONFIG_HMMEL | CM_CONFIG_TRUNC;
  if(cm_Configure(cm, errbuf, -1) != eslOK) esl_fatal("cm_Configure failed: %s", errbuf);
  if(v) printf("%s: M %d, EL flag %s\n", cm->name, cm->cp9->M, (cm->cp9->flags & CPLAN9_EL) ? "on" : "off");

  /* cm->cp9 (cmalign --notrunc) and cm->Tcp9 (default) */
  utest_narrow_band(rng, cm->cp9,  "cp9",  abc->K, N, v);
  utest_narrow_band(rng, cm->Tcp9, "Tcp9", abc->K, N, v);
  utest_full_band  (rng, cm->cp9,  "cp9",  abc->K, 5, v);
  utest_full_band  (rng, cm->Tcp9, "Tcp9", abc->K, 5, v);

  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#undef RX
#undef T
#undef REF_MSC
#undef REF_ISC
#undef MASK
#endif /*CM_P7_BAND_CHK_TESTDRIVE*/

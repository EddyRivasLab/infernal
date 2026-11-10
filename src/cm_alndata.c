/* CM_ALNDATA: data structure containing information relevant to the
 * alignment of a sequence and its output. Currently, mainly used by
 * cmalign, but also used in cmbuild with the --refine option.
 * 
 * Contents:
 *    1. The CM_ALNDATA object.
 *    2. Alignment workunit processing functions, which create and 
 *       fill CM_ALNDATA objects.
 *
 * EPN, Fri Jan  6 09:00:31 2012
 */
#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <inttypes.h>
#include <math.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

/* The largest |c| (bits per residue) center_shift() uses to center CM
 * alignment DP scores; a larger CP9 score per residue is not a usable
 * estimate of the CM's (log-odds scores per residue are a few bits). */
#define CM_CENTER_MAXSHIFT 20.

static int  plain_cp9_engine(CM_t *cm, char *errbuf, ESL_SQ *sq, float mxsize, int nmx, int do_iterate, int *ret_use_ckpt);
static void shrink_cp9_matrices(CM_t *cm);

/* rung-3: derive the per-bifurcation right-fragment length k* (kpin[v]) from a
 * CYK parsetree, exactly as the rung3_drv extract_bifs harness does.  kpin[v]
 * stays -1 for B states not on the parse (none, in global mode -- every B is
 * visited).  Used to pin cm_CheckptPostAlignHB / cm_CheckptOptAccAlignHB. */
static void
rung3_kpin_from_cyk(CM_t *cm, Parsetree_t *tr, int *kpin)
{
  int n;
  for (n = 0; n < cm->M; n++) kpin[n] = -1;
  for (n = 0; n < tr->n; n++) {
    int v = tr->state[n];
    if (cm->sttype[v] != B_st) continue;
    int ln = tr->nxtl[n], rn = tr->nxtr[n];
    int lstid = cm->stid[tr->state[ln]];
    int begl_node = (lstid == BEGL_S) ? ln : rn;      /* left  (BEGL) subtree node */
    int ksplit    = tr->emitr[begl_node];             /* last residue of left fragment */
    kpin[v]       = tr->emitr[n] - ksplit;            /* right-fragment length k* */
  }
}

/* rung-4: derive the per-bifurcation pin skeleton (kind, k*, B/left/right modes)
 * from a TRUNCATED CYK D&C parsetree, exactly as the rung-4 troa_drv extract_tr_pins
 * harness does.  Truncated analogue of rung3_kpin_from_cyk that additionally
 * records the bifurcation kind (1=INTERIOR, 2=LEFT_FULL, 3=RIGHT_FULL) and the
 * per-B marginal modes the engine traceback needs.  All five arrays are sized
 * cm->M; B states not on the parse stay {kind=0, k*=-1, modes=UNKNOWN} (none, in
 * global mode -- every B is visited).  Used to pin cm_CheckptTrPostAlignHB /
 * cm_CheckptTrOptAccAlignHB. */
static void
rung4_trpins_from_cyk(CM_t *cm, Parsetree_t *tr,
                      int *bkind, int *kpin, char *bbmode, char *blmode, char *brmode)
{
  int v, n;
  for (v = 0; v < cm->M; v++) { bkind[v] = 0; kpin[v] = -1; bbmode[v] = blmode[v] = brmode[v] = TRMODE_UNKNOWN; }
  for (n = 0; n < tr->n; n++) {
    v = tr->state[n];
    if (cm->sttype[v] != B_st) continue;
    int ln = tr->nxtl[n], rn = tr->nxtr[n];
    int lstid = cm->stid[tr->state[ln]], rstid = cm->stid[tr->state[rn]];
    int begl_node, begr_node;
    if      (lstid == BEGL_S && rstid == BEGR_S) { begl_node = ln; begr_node = rn; }
    else if (lstid == BEGR_S && rstid == BEGL_S) { begl_node = rn; begr_node = ln; }
    else                                         { begl_node = ln; begr_node = rn; }
    int lspan = (tr->emitl[begl_node] <= tr->emitr[begl_node]) ? (tr->emitr[begl_node] - tr->emitl[begl_node] + 1) : 0;
    int rspan = (tr->emitl[begr_node] <= tr->emitr[begr_node]) ? (tr->emitr[begr_node] - tr->emitl[begr_node] + 1) : 0;
    bbmode[v] = tr->mode[n];
    blmode[v] = tr->mode[begl_node];
    brmode[v] = tr->mode[begr_node];
    if      (lspan == 0) { bkind[v] = 3; kpin[v] = -1;    }  /* RIGHT_FULL: left empty (k*=d at traceback) */
    else if (rspan == 0) { bkind[v] = 2; kpin[v] = 0;     }  /* LEFT_FULL:  right empty (k*=0)             */
    else                 { bkind[v] = 1; kpin[v] = rspan; }  /* INTERIOR:   k*=right span                  */
  }
}

/* brief 26_0430-234: total d-band cell count across all CM states -- a proxy for
 * pass-2's HB DP work.  Mirrors the inline band-area loop the --cykbands pre-pass
 * used, factored out; cykbands_tighten_bounded() uses it for its blowup guard. */
static double
ckpt_cykbands_cellcount(CM_t *cm)
{
  CP9Bands_t *cp9b = cm->cp9b;
  double cells = 0.;
  int v, jp;
  for(v = 0; v < cm->M; v++)
    for(jp = 0; jp <= cp9b->jmax[v] - cp9b->jmin[v]; jp++)
      if(hd_min(cp9b, v, jp) <= hd_max(cp9b, v, jp))
        cells += hd_max(cp9b, v, jp) - hd_min(cp9b, v, jp) + 1;
  return cells;
}

/* brief 26_0430-243: env-agnostic, diagnostic-neutral CORE of the CYK-band
 * tightening + never-loosen bounding, factored out of ckpt_cykbands_tighten()
 * (below) so the standalone --cykbands path can share the *exact same* bounded
 * mechanism.  Before this, --cykbands called cm_BandsFromCYKParsetree() raw --
 * no baseline snapshot, no never-loosen intersection, no bounded hd-recompute --
 * so the unvisited-state full-envelope fallback (brief 26_0430-065)
 * could inflate the pass-2 matrix without bound and outright fail at genome scale
 * (brief 26_0430-242: band_area_ratio up to 145924x, 343 GB matrix, 4 dengue
 * cells ESL_XFAIL'd).  The --ckpt-cykbands path never blew up because it already
 * had the never-loosen bounding (brief 26_0430-238); this helper hands that same
 * bounding to --cykbands.
 *
 * <pad> and <maxratio> are plain arguments -- this core reads NO env vars.  The
 * CKPT_CYKBANDS_PAD / CKPT_CYKBANDS_MAXRATIO overrides are applied ONLY by the
 * ckpt_cykbands_tighten() wrapper; the --cykbands caller passes cm->p7_cykbands_pad
 * (from --cykpad) and a plain constant backstop, so --cykbands is NOT overridable
 * by the CKPT_* env vars (brief 26_0430-243 constraint 1).  Diagnostic counts are
 * returned via out-params (this core emits NO fprintf); each caller prints its
 * own labeled message.  All band-mutation semantics are otherwise identical to
 * the pre-refactor ckpt_cykbands_tighten() -- see that wrapper's doc comment for
 * the never-loosen / preserve_valid / blowup-guard rationale.
 *
 * <ret_n_sentinel>/<ret_n_real_empty>: never-loosen empty-intersection counts
 *   (CP9 baseline already empty vs real accuracy-risk empty; brief 26_0430-238).
 * <ret_revert_ratio>: 0.0 if the blowup-guard did not fire; else the
 *   (pre-revert tight)/orig ratio that triggered the revert-to-baseline.
 *
 * Returns eslOK (bands tightened, or reverted-to-baseline if the guard fired) or
 * a failure status (validity restored, snapshots freed). */
static int
cykbands_tighten_bounded(CM_t *cm, char *errbuf, Parsetree_t *tr, int L, int pass_idx,
                         int preserve_valid, int pad, double maxratio,
                         double *ret_orig, double *ret_tight,
                         int64_t *ret_n_sentinel, int64_t *ret_n_real_empty,
                         double *ret_revert_ratio)
{
  int    status;
  int    M   = cm->M;
  int    nv  = cm->M + 1;
  int   *Jv = NULL, *Lv = NULL, *Rv = NULL, *Tv = NULL;
  /* brief 26_0430-237: snapshot of the pre-tighten (untightened baseline) spatial
   * bands, so a low-coverage sequence whose CYK-parsetree tightening would BLOW
   * UP the band area -- most of a large CM's states unvisited, each falling back
   * to the full [i0..j0] envelope, contributing O(L^2) d-cells apiece -- can be
   * safely REVERTED to the untightened bands rather than committing pass-2 to an
   * unbounded (hundreds-of-billions-of-cell) DP matrix (hang / OOM). */
  int   *s_imin=NULL, *s_imax=NULL, *s_jmin=NULL, *s_jmax=NULL, *s_hddn=NULL;
  int64_t s_hd_needed = 0, s_hd_alloced = 0;
  int64_t n_cp9_sentinel = 0, n_real_empty = 0;
  double  orig, tight, revert_ratio = 0.;

  orig = ckpt_cykbands_cellcount(cm);

  /* snapshot untightened spatial bands (+ derived hd_dn / hd_needed) for revert */
  ESL_ALLOC(s_imin, sizeof(int)*M); ESL_ALLOC(s_imax, sizeof(int)*M);
  ESL_ALLOC(s_jmin, sizeof(int)*M); ESL_ALLOC(s_jmax, sizeof(int)*M);
  ESL_ALLOC(s_hddn, sizeof(int)*M);
  memcpy(s_imin, cm->cp9b->imin,  sizeof(int)*M);
  memcpy(s_imax, cm->cp9b->imax,  sizeof(int)*M);
  memcpy(s_jmin, cm->cp9b->jmin,  sizeof(int)*M);
  memcpy(s_jmax, cm->cp9b->jmax,  sizeof(int)*M);
  memcpy(s_hddn, cm->cp9b->hd_dn, sizeof(int)*M);
  s_hd_needed  = cm->cp9b->hd_needed;
  s_hd_alloced = cm->cp9b->hd_alloced;

  if(preserve_valid) {
    ESL_ALLOC(Jv, sizeof(int)*nv); ESL_ALLOC(Lv, sizeof(int)*nv);
    ESL_ALLOC(Rv, sizeof(int)*nv); ESL_ALLOC(Tv, sizeof(int)*nv);
    memcpy(Jv, cm->cp9b->Jvalid, sizeof(int)*nv);
    memcpy(Lv, cm->cp9b->Lvalid, sizeof(int)*nv);
    memcpy(Rv, cm->cp9b->Rvalid, sizeof(int)*nv);
    memcpy(Tv, cm->cp9b->Tvalid, sizeof(int)*nv);
  }

  status = cm_BandsFromCYKParsetree(cm, errbuf, tr, 1, L, pad, FALSE,
                                    cm->cp9b, pass_idx, 0);

  if(preserve_valid) {
    memcpy(cm->cp9b->Jvalid, Jv, sizeof(int)*nv);
    memcpy(cm->cp9b->Lvalid, Lv, sizeof(int)*nv);
    memcpy(cm->cp9b->Rvalid, Rv, sizeof(int)*nv);
    memcpy(cm->cp9b->Tvalid, Tv, sizeof(int)*nv);
    free(Jv); free(Lv); free(Rv); free(Tv);
    Jv = Lv = Rv = Tv = NULL;
  }

  if(status != eslOK) goto ERROR;

  /* brief 26_0430-238: NEVER-LOOSEN -- intersect the just-computed CYK-tightened
   * spatial bands with the snapshotted untightened (CP9-posterior) bands, per
   * state. The CYK tightening is supposed to only ever SHRINK a state's band;
   * for a state unvisited by the pass-1 CYK parse (common in local mode for a
   * low-coverage/fragment sequence -- most of a large CM's states) the
   * full-envelope fallback above instead EXPANDS it past the CP9 baseline
   * (brief 26_0430-065 / brief 26_0430-237's measured 94K->190M-cell
   * blowup on the LSU fragment repro). Intersecting bounds the result to
   * [never worse than the untightened baseline]: unvisited states'
   * full-envelope CYK band always contains the CP9 band, so the intersection
   * collapses back to exactly the CP9 band (safe, unchanged cost, and -- per
   * brief 26_0430-235 -- still covers every state OptAcc/MEA might diverge into);
   * visited states' CYK band already sits inside (or pad-widens slightly past)
   * the CP9 band, so the intersection preserves the brief 26_0430-234/26_0430-236 tightening
   * win. This makes brief 26_0430-237's revert-on-blowup ratio cap a pure backstop
   * that should no longer fire on real data -- the mechanism is now
   * intrinsically bounded rather than merely capped. */
  {
    int     v;
    for(v = 0; v < M; v++) {
      int cp9_is_sentinel = (s_imin[v] > s_imax[v] || s_jmin[v] > s_jmax[v]);
      int ni_min = ESL_MAX(cm->cp9b->imin[v], s_imin[v]);
      int ni_max = ESL_MIN(cm->cp9b->imax[v], s_imax[v]);
      int nj_min = ESL_MAX(cm->cp9b->jmin[v], s_jmin[v]);
      int nj_max = ESL_MIN(cm->cp9b->jmax[v], s_jmax[v]);
      if(ni_min > ni_max || nj_min > nj_max) {
        /* Empty intersection. Never leave a state silently unreachable here
         * (brief 26_0430-238's safety requirement) -- fall back to the CP9
         * (pre-tighten) band for this state alone, whatever it is
         * (sentinel or real). */
        ni_min = s_imin[v]; ni_max = s_imax[v];
        nj_min = s_jmin[v]; nj_max = s_jmax[v];
        if(cp9_is_sentinel) n_cp9_sentinel++;   /* CP9 baseline itself already unreachable (expected, harmless) */
        else                n_real_empty++;     /* real accuracy-risk empty (brief 26_0430-238; should be rare/never) */
      }
      cm->cp9b->imin[v] = ni_min; cm->cp9b->imax[v] = ni_max;
      cm->cp9b->jmin[v] = nj_min; cm->cp9b->jmax[v] = nj_max;
    }
    /* hd_dn[v] is a pure function of state type + do_trunc (ij2d_bands), not of
     * band width, so it is unaffected by the intersection above and does not
     * strictly need recomputing -- but call both to stay in lockstep with
     * cm_BandsFromCYKParsetree's own post-band-fill sequence and keep
     * hd_needed (diagnostic) consistent with the now-narrower jmin/jmax. */
    if((status = cp9_GrowHDBands(cm->cp9b, errbuf)) != eslOK) goto ERROR;
    ij2d_bands(cm, cm->cp9b, cm_pli_PassAllowsTruncation(pass_idx), 0);
  }

  tight = ckpt_cykbands_cellcount(cm);

  /* brief 26_0430-237: blowup guard -- if tightening inflated the band area past
   * the ratio cap, revert to the untightened baseline bands. Safe: this one
   * sequence simply forgoes the CYK-tightening speedup and aligns with exactly
   * the bands it would have used with tightening off (no accuracy change vs off).
   * With the never-loosen intersection above this is a pure backstop that should
   * no longer fire on real data. */
  if(maxratio > 0. && orig > 0. && tight > maxratio * orig) {
    memcpy(cm->cp9b->imin,  s_imin, sizeof(int)*M);
    memcpy(cm->cp9b->imax,  s_imax, sizeof(int)*M);
    memcpy(cm->cp9b->jmin,  s_jmin, sizeof(int)*M);
    memcpy(cm->cp9b->jmax,  s_jmax, sizeof(int)*M);
    memcpy(cm->cp9b->hd_dn, s_hddn, sizeof(int)*M);
    cm->cp9b->hd_needed  = s_hd_needed;
    cm->cp9b->hd_alloced = s_hd_alloced;
    revert_ratio = tight / orig;
    tight = orig;   /* effective (post-revert) cell count */
  }

  free(s_imin); free(s_imax); free(s_jmin); free(s_jmax); free(s_hddn);
  if(ret_orig)         *ret_orig         = orig;
  if(ret_tight)        *ret_tight        = tight;
  if(ret_n_sentinel)   *ret_n_sentinel   = n_cp9_sentinel;
  if(ret_n_real_empty) *ret_n_real_empty = n_real_empty;
  if(ret_revert_ratio) *ret_revert_ratio = revert_ratio;
  return eslOK;

 ERROR:
  if(Jv) free(Jv);
  if(Lv) free(Lv);
  if(Rv) free(Rv);
  if(Tv) free(Tv);
  if(s_imin) free(s_imin);
  if(s_imax) free(s_imax);
  if(s_jmin) free(s_jmin);
  if(s_jmax) free(s_jmax);
  if(s_hddn) free(s_hddn);
  return status;
}

/* brief 26_0430-234: piggyback spatial band tightening on --ckpt's OWN pass-1
 * CYK/D&C parsetree (<tr>), before it is freed, to shrink pass-2's Inside/Outside/
 * Posterior/OptAcc band area at ZERO extra CYK cost -- the parse was going to be
 * computed and discarded anyway (only its bifurcation k* pins are used today).
 * This is NOT the separate --cykbands mechanism (which runs its own extra, full-
 * memory, non-checkpointed CYK pre-pass, defeating --ckpt's whole point); this
 * reuses the checkpointed pass-1 tree already in hand.  Off by default; on with
 * --ckpt-cykbands.
 *
 * <preserve_valid> TRUE for the truncated rung-4 path.  cm_BandsFromCYKParsetree's
 * do_trunc branch declares truncated bands "not supported" and clobbers
 * cp9b->{J,L,R,T}valid to J-only, which would break rung-4 pass-2's per-mode
 * validity gating (cm_CheckptTr{Post,OptAcc}AlignHB skip cells on
 * mode==TRMODE_L && !Lvalid[v], etc).  BUT ij2d_bands recomputes the on-demand
 * hd_dn[] d-band floor from do_trunc + state type ONLY -- never from the valid
 * flags -- so the tightened spatial i/j/hd bands ARE correct for truncated mode;
 * only the mode-validity is wrong.  We therefore save cp9b's validity across the
 * call and restore it: spatial bands tighten, the resolved r4_mode stays valid.
 * For the non-truncated rung-3 path (preserve_valid FALSE) the non-trunc branch's
 * Jvalid-all-TRUE / L,R,Tvalid-FALSE already matches rung-3's expected validity,
 * so no save/restore is needed.
 *
 * Correctness of pins vs bands: the bifurcation k* pins were just extracted from
 * the SAME <tr>.  Every pinned B state is a VISITED state, so its tightened
 * band is [visited emitl/emitr] +/- pad and necessarily contains the pinned
 * cell.  Unvisited B states carry no pin (kpin=-1), so the phantom-wide-band
 * inheritance (briefs 26_0430-065/26_0430-068) can only widen an unpinned
 * state -- it can never exclude a pin.  No pin/band conflict is possible.
 *
 * brief 26_0430-237: for a low-coverage sequence (most of a large CM's states
 * unvisited by the pass-1 CYK parse) the unvisited-state full-envelope fallback
 * in cm_BandsFromCYKParsetree can INFLATE the band area by 30x-40000x (O(L^2)
 * d-cells per fallback state), producing a hundreds-of-billions-of-cell pass-2
 * matrix that hangs or OOM-kills the process. Guard: after tightening, if the
 * cell count exceeds <maxratio> x the untightened baseline (3.0),
 * REVERT to the snapshotted untightened --ckpt bands.
 *
 * Returns eslOK (bands tightened, or reverted-to-untightened if the guard fired;
 * *ret_orig/*ret_tight set to pre/post d-cell counts, *ret_tight == *ret_orig on
 * revert) or a failure status (validity restored, bands left as cm_Bands... left
 * them). */
static int
ckpt_cykbands_tighten(CM_t *cm, char *errbuf, Parsetree_t *tr, int L, int pass_idx,
                      int preserve_valid, double *ret_orig, double *ret_tight)
{
  int         status;
  int         pad      = cm->p7_cykbands_pad;   /* default 5 (cm.c); shared with --cykbands */
  double      maxratio = 3.0;                   /* brief 26_0430-237: revert if band-cells > maxratio x
                                                 * untightened baseline. Legit widening tops out ~2.24x
                                                 * (5S_rRNA pad=8, brief 26_0430-236); low-coverage blowup is
                                                 * 30x-40000x (brief 26_0430-237); 3.0 cleanly separates them.
                                                 * With never-loosen this is a pure backstop. */

  status = cykbands_tighten_bounded(cm, errbuf, tr, L, pass_idx, preserve_valid,
                                    pad, maxratio, ret_orig, ret_tight,
                                    NULL, NULL, NULL);
  return status;
}

/*****************************************************************
 * 1. The CM_ALNDATA object
 *****************************************************************/

/* Function:  cm_alndata_Create()
 * Synopsis:  Allocate a CM_ALNDATA object.
 * Incept:    EPN, Fri Jan  6 09:01:33 2012
 *            
 * Purpose:   Allocates a new <CM_ALNDATA> and returns a pointer
 *            to it.
 *
 * Throws:    <NULL> on allocation failure.
 */
/* plain_cp9_engine()
 * Incept: EPN, Mon Oct 12 2026 (w/Claude)
 *
 * The --mxsize gate for the CP9 HMM band calculation of the default
 * (not p7-banded) path, which is independent of the CM DP's: the two
 * stages do not hold their matrices at the same time (see
 * shrink_cp9_matrices()). The full-matrix calculation holds <nmx> CP9
 * matrices of (L+1) x (M+1) cells (Forward, Backward and, if bands are
 * tightened, a posterior matrix). If they exceed <mxsize>, use the
 * checkpointed calculation, which gives the same bands in O(sqrt(L) * M)
 * memory. If that exceeds <mxsize> too, fail with eslERANGE and a
 * message giving the --mxsize that would fit: --mxsize is a limit the
 * user chose, and this stage can go no lower.
 *
 * Returns eslOK, or eslERANGE (errbuf filled) if no CP9 band
 * calculation fits <mxsize>.
 */
static int
plain_cp9_engine(CM_t *cm, char *errbuf, ESL_SQ *sq, float mxsize, int nmx, int do_iterate, int *ret_use_ckpt)
{
  float       dense_Mb = (float) nmx * SizeNeededCP9Matrix(sq->L, cm->cp9->M, NULL, NULL);
  float       ckpt_Mb  = cp9_CheckptSeq2BandsMbNeeded(cm, sq->L, do_iterate, cm->maxtau);
  int         usable   = cp9_CheckptBandsUsable(cm, FALSE);
  int         use_ckpt = (usable && dense_Mb > mxsize) ? TRUE : FALSE;
  float       need_Mb  = use_ckpt ? ckpt_Mb : dense_Mb;

  if(cm_StageVerbose()) fprintf(stderr, "#CP9ENGINE seq=%s L=%d M=%d engine=%s dense_Mb=%.1f ckpt_Mb=%.1f mxsize=%.1f\n",
                                sq->name, (int) sq->L, cm->cp9->M, use_ckpt ? "ckpt" : "dense", dense_Mb, ckpt_Mb, mxsize);
  *ret_use_ckpt = use_ckpt;
  if(need_Mb > mxsize)
    ESL_FAIL(eslERANGE, errbuf, "HMM band calculation needs %.0f Mb > %.0f Mb limit.\nUse --mxsize %d.",
             need_Mb, mxsize, (int) (floor(need_Mb / 1024.) + 1.) * 1024);
  return eslOK;
}

/* shrink_cp9_matrices()
 * Incept: EPN, Mon Oct 12 2026 (w/Claude)
 *
 * Shrink cm->cp9_mx and cm->cp9_bmx back to the two rows they are
 * created with, so that a sequence's CP9 band calculation matrices are
 * not still allocated during its CM DP alignment (or a later sequence's).
 * They grow again when next needed.
 */
static void
shrink_cp9_matrices(CM_t *cm)
{
  if(cm->cp9 == NULL) return;
  if(cm->cp9_mx  != NULL && (cm->cp9_mx->rows  > 1 || cm->cp9_mx->ncells_allocated  > 2 * (int64_t) (cm->cp9->M+1))) {
    FreeCP9Matrix(cm->cp9_mx);
    cm->cp9_mx  = CreateCP9Matrix(1, cm->cp9->M);
  }
  if(cm->cp9_bmx != NULL && (cm->cp9_bmx->rows > 1 || cm->cp9_bmx->ncells_allocated > 2 * (int64_t) (cm->cp9->M+1))) {
    FreeCP9Matrix(cm->cp9_bmx);
    cm->cp9_bmx = CreateCP9Matrix(1, cm->cp9->M);
  }
}

CM_ALNDATA *
cm_alndata_Create(void)
{
  int status;
  CM_ALNDATA *data = NULL;

  ESL_ALLOC(data, sizeof(CM_ALNDATA));
  data->sq         = NULL;
  data->idx        = -1;
  data->tr         = NULL;
  data->sc         = 0.;
  data->pp         = 0.;
  data->ppstr      = NULL;
  data->spos       = -1;
  data->epos       = -1;
  data->secs_bands = 0.;
  data->secs_aln   = 0.;
  data->mb_tot     = 0.;
  data->tau        = -1.;
  
  return data;

 ERROR: 
  return NULL;
}

/* Function:  cm_alndata_Destroy()
 * Synopsis:  Free a CM_ALNDATA object.
 * Incept:    EPN, Fri Jan  6 09:10:55 2012
 *            
 * Purpose:  Frees a <CM_ALNDATA> object, but only frees the
 *           ESL_SQ <sq> if <free_sq> is TRUE. Often this is
 *           only a pointer to a sequence in another data 
 *           structure that will be free'd with that structure.
 *
 * Returns:  void.
 */
void
cm_alndata_Destroy(CM_ALNDATA *data, int free_sq)
{ 
  if(data == NULL) return;

  if(free_sq && data->sq != NULL) esl_sq_Destroy(data->sq);
  if(data->tr    != NULL)         FreeParsetree(data->tr);
  if(data->ppstr != NULL)         free(data->ppstr);
  free(data);

  return;
}

/*****************************************************************
 * 2. Alignment workunit processing functions
 *****************************************************************/

/* Function: DispatchSqBlockAlignment()
 * Date:     EPN, Fri Dec 30 14:59:43 2011
 *
 * Purpose:  Given a CM and a block of sequences, align the
 *           sequence(s) using the appropriate alignment function and
 *           return relevant data for eventual output in <ret_dataA>.
 *           This function simply calls DispatchSqAlignment() serially
 *           for each sequence in the block, and creates an array
 *           of the <ret_data> DispatchSqAlignment() returns.
 *
 *           Currently <mode>, <cp9b_valid> and <pass_idx> values sent
 *           to DispatchSqAlignment() are hard-coded to
 *           TRMODE_UNKNOWN, FALSE, and PLI_PASS_5P_AND_3P_FORCE (if
 *           cm->align_opts & CM_ALIGN_TRUNC) or PLI_PASS_STD_ANY (if
 *           (! cm->align_opts & CM_ALIGN_TRUNC)). This is because
 *           this function is only used by the alignment pipeline, in
 *           which these values are correct. If this changes, we may
 *           want caller to pass in an array of modes, cp9b_valids and
 *           pass_idx values, one per sq.
 *
 *           If (cm->flags & CM_ALIGN_XTAU) we'll potentially tighten
 *           HMM bands until the required DP matrices are below out
 *           limit (<mxsize>). cm->maxtau is the max allowed tau value
 *           during this iterative band tightening, and cm->xtau is
 *           the factor by which we multiply cm->tau at each iteration
 *           during band tightening.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           sq_block  - block of sequences to align
 *           mxsize    - max size in Mb of allowable DP mx
 *           w         - stopwatch for timing individual stages
 *           w_tot     - stopwatch for timing total time per seq
 *           r         - RNG, req'd if CM_ALIGN_SAMPLE, can be NULL otherwise
 *           ret_dataA - RETURN: newly created array of CM_ALNDATA objects
 *
 * Returns:  eslOK on success;
 *           eslEINCOMPAT on contract violation, errbuf is filled;
 *           eslEMEM if we run out of memory;
 *           <ret_dataA> is alloc'ed and filled with sq_block->count CM_ALNDATA objects.
 */
int
DispatchSqBlockAlignment(CM_t *cm, char *errbuf, ESL_SQ_BLOCK *sq_block, float mxsize, ESL_STOPWATCH *w, 
			 ESL_STOPWATCH *w_tot, ESL_RANDOMNESS *r, CM_ALNDATA ***ret_dataA)
{
  int           status;          /* easel status */
  int           j;               /* counter over parsetrees */
  CM_ALNDATA  **dataA = NULL;    /* CM_ALNDATA array we'll create and return */
  ESL_SQ       *sqp;             /* ptr to a ESL_SQ */
  int           pass_idx;        /* pass_idx passed to DispatchSqAlignment() */
  char          mode;            /* mode passed to DispatchSqAlignment() */
  int           cp9b_valid;      /* passed to DispatchSqAlignment() */

  ESL_ALLOC(dataA, sizeof(CM_ALNDATA *) * ESL_MAX(1, sq_block->count)); // avoid 0 malloc
  for(j = 0; j < sq_block->count; j++) dataA[j] = NULL;

  /* DispatchSqAligment() needs a mode, pipeline pass index, and
   * knowledge of whether cm->cp9b are valid for sequence to align
   * (see note in 'Purpose' above). Currently the relevant values 
   * for these are as follows:
   */
  mode       = TRMODE_UNKNOWN;
  pass_idx   = (cm->align_opts & CM_ALIGN_TRUNC) ? PLI_PASS_5P_AND_3P_FORCE : PLI_PASS_STD_ANY; 
  cp9b_valid = FALSE;

  /* main loop: for each sequence, call DispatchSqAlignment() to do the work */
  for(j = 0; j < sq_block->count; j++) {
    sqp = sq_block->list + j;
    if((status = DispatchSqAlignment(cm, errbuf, sqp, sq_block->first_seqidx + j, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, &(dataA[j]))) != eslOK) goto ERROR;
  }
  *ret_dataA = dataA;

  return eslOK;

 ERROR: 
  if(dataA != NULL) { 
    for(j = 0; j < sq_block->count; j++) { 
      if(dataA[j] != NULL) cm_alndata_Destroy(dataA[j], FALSE);
    }
    free(dataA);
  }
  *ret_dataA = NULL;
  fprintf(stderr, "Problem during alignment of sequence %s\n", sqp->name);
  if(status == eslEMEM) ESL_FAIL(status, errbuf, "DispatchSqBlockAlignment(), out of memory");
  else return status; /* errbuf was filled by DispatchSqAlignment() */
}

/* Local-mode EL detect-and-retry for --p7kmerchain bands (briefs
 * 26_0821-127/26_0821-128/26_0821-132): the per-call band state that
 * DispatchSqAlignment() passes to dispatch_sq_alignment_core().
 *
 * A first pass (mode ELR_MODE_NONE, capture TRUE) returns a copy of its
 * kmerchain band and the chain's pins before junction slack. A retry either
 * reuses that band with kmin lowered on some rows (ELR_MODE_EDIT; the chain is
 * not rebuilt, so the band is the first pass's band plus the edit) or aligns
 * without a p7 band, i.e. with native CP9 bands (ELR_MODE_NATIVE).
 */
#define ELR_MODE_NONE   0
#define ELR_MODE_EDIT   1
#define ELR_MODE_NATIVE 2
typedef struct {
  int  mode;      /* ELR_MODE_NONE | ELR_MODE_EDIT | ELR_MODE_NATIVE                     */
  int  capture;   /* TRUE: return the kmerchain band below (first pass only)             */
  int  L;         /* length of the arrays below, [0..L]                                  */
  int *raw;       /* chain pins before junction slack (-1 = unpinned)                    */
  int *i2k;       /* pins after junction slack, as the band was built from               */
  int *kmin;      /* the band; ELR_MODE_EDIT: the band to reuse (not owned)              */
  int *kmax;
  int *klo;       /* ELR_MODE_EDIT: [1..L] lower kmin to this where 0 < klo[i] < kmin[i] */
} ELR_BAND;

static void
elr_band_release(ELR_BAND *eb)
{
  if(eb->raw)  free(eb->raw);
  if(eb->i2k)  free(eb->i2k);
  if(eb->kmin) free(eb->kmin);
  if(eb->kmax) free(eb->kmax);
  eb->raw = eb->i2k = eb->kmin = eb->kmax = NULL;
  eb->L = 0;
}

/* ELR_MODE_EDIT: copies of <eb>'s band with the kmin edit applied, and the new cell count */
static int
elr_band_edited(ELR_BAND *eb, int L, int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int     status;
  int    *i2k = NULL, *kmin = NULL, *kmax = NULL;
  int64_t n   = 0;
  int     i;

  ESL_ALLOC(i2k,  sizeof(int) * (L+1));
  ESL_ALLOC(kmin, sizeof(int) * (L+1));
  ESL_ALLOC(kmax, sizeof(int) * (L+1));
  memcpy(i2k,  eb->i2k,  sizeof(int) * (L+1));
  memcpy(kmin, eb->kmin, sizeof(int) * (L+1));
  memcpy(kmax, eb->kmax, sizeof(int) * (L+1));
  for(i = 1; i <= L; i++) {
    if(eb->klo[i] > 0 && eb->klo[i] < kmin[i]) kmin[i] = eb->klo[i];
    n += (int64_t) (kmax[i] - kmin[i] + 1);
  }
  *ret_i2k = i2k; *ret_kmin = kmin; *ret_kmax = kmax;
  *ret_ncells = (n > (int64_t) INT_MAX) ? INT_MAX : (int) n;
  return eslOK;

 ERROR:
  if(i2k)  free(i2k);
  if(kmin) free(kmin);
  if(kmax) free(kmax);
  return eslEMEM;
}

/* First pass: keep a copy of the kmerchain band (<raw> is handed over) */
static int
elr_band_keep(ELR_BAND *eb, int L, int *raw, const int *i2k, const int *kmin, const int *kmax)
{
  int status;

  elr_band_release(eb);
  eb->raw = raw;
  ESL_ALLOC(eb->i2k,  sizeof(int) * (L+1));
  ESL_ALLOC(eb->kmin, sizeof(int) * (L+1));
  ESL_ALLOC(eb->kmax, sizeof(int) * (L+1));
  memcpy(eb->i2k,  i2k,  sizeof(int) * (L+1));
  memcpy(eb->kmin, kmin, sizeof(int) * (L+1));
  memcpy(eb->kmax, kmax, sizeof(int) * (L+1));
  eb->L = L;
  return eslOK;

 ERROR:
  elr_band_release(eb);
  return eslEMEM;
}

/* Function:  center_shift()
 * Incept:    EPN, Fri Oct 30 2026 (w/Claude)
 *
 * Purpose:   Return the per-residue shift <c> (bits) for centering the
 *            CM alignment DP scores of a sequence of length <L> (see
 *            DispatchSqAlignment()): the CP9 Forward score of the HMM
 *            band calculation that just ran, divided by <L>. This is
 *            close to the mean score per residue of the CM alignment.
 *            Return 0. (no centering) if that score is not available
 *            (HMM bands not calculated, or calculated in scan mode)
 *            or is not a usable value.
 */
static float
center_shift(CM_t *cm, int L)
{
  double c;

  if(cm->cp9b == NULL || L <= 0)         return 0.;
  if(! isfinite(cm->cp9b->hmm_sc))       return 0.;
  c = cm->cp9b->hmm_sc / (double) L;
  if(fabs(c) > CM_CENTER_MAXSHIFT)       return 0.; /* e.g. a CP9 total near -INFTY: no parse in the HMM bands */
  return (float) c;
}

/* Function:  center_emission_scores()
 * Incept:    EPN, Fri Oct 30 2026 (w/Claude)
 *
 * Purpose:   Lower the score of every residue <cm> can emit by <c>
 *            bits: the emission scores (cm->oesc, cm->esc) of ML, MR,
 *            IL and IR states and the marginal emission scores
 *            (cm->lmesc, cm->rmesc) of MP states by <c>, the pair
 *            emission scores of MP states by 2c, and the EL
 *            self-transition score, which is the score of each residue
 *            an EL state emits, by <c>. These are all the scores the
 *            CM alignment DP functions add per residue. The lmesc and
 *            rmesc of single-emitter states are copies the DP functions
 *            don't read, and are left alone.
 *
 *            The original scores are saved in <*ret_save> for
 *            restore_emission_scores(), which puts them back exactly
 *            (subtracting and adding <c> back in float would not).
 *
 * Returns:   eslOK on success, eslEMEM if out of memory.
 */
static int
center_emission_scores(CM_t *cm, float c, float **ret_save)
{
  int     status;
  int     v, a;
  int     K     = cm->abc->K;
  int     Kp    = cm->abc->Kp;
  int64_t n;
  int64_t nsave = 1;
  float  *save  = NULL;

  for(v = 0; v < cm->M; v++) {
    if     (cm->sttype[v] == MP_st)                                    nsave += Kp*Kp + K*K + 2*Kp;
    else if(cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
            cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)          nsave += Kp + K;
  }
  ESL_ALLOC(save, sizeof(float) * nsave);

  n = 0;
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case ML_st: case MR_st: case IL_st: case IR_st:
      for(a = 0; a < Kp; a++) { save[n++] = cm->oesc[v][a]; cm->oesc[v][a] -= c; }
      for(a = 0; a < K;  a++) { save[n++] = cm->esc[v][a];  cm->esc[v][a]  -= c; }
      break;
    case MP_st:
      for(a = 0; a < Kp*Kp; a++) { save[n++] = cm->oesc[v][a];  cm->oesc[v][a]  -= 2.*c; }
      for(a = 0; a < K*K;   a++) { save[n++] = cm->esc[v][a];   cm->esc[v][a]   -= 2.*c; }
      for(a = 0; a < Kp;    a++) { save[n++] = cm->lmesc[v][a]; cm->lmesc[v][a] -= c;    }
      for(a = 0; a < Kp;    a++) { save[n++] = cm->rmesc[v][a]; cm->rmesc[v][a] -= c;    }
      break;
    }
  }
  save[n++] = cm->el_selfsc;
  cm->el_selfsc -= c;

  *ret_save = save;
  return eslOK;

 ERROR:
  *ret_save = NULL;
  return status;
}

/* Function:  restore_emission_scores()
 * Incept:    EPN, Fri Oct 30 2026 (w/Claude)
 *
 * Purpose:   Put back the scores center_emission_scores() lowered,
 *            from <save>, and free <save>.
 */
static void
restore_emission_scores(CM_t *cm, float *save)
{
  int     v, a;
  int     K  = cm->abc->K;
  int     Kp = cm->abc->Kp;
  int64_t n  = 0;

  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case ML_st: case MR_st: case IL_st: case IR_st:
      for(a = 0; a < Kp; a++) cm->oesc[v][a] = save[n++];
      for(a = 0; a < K;  a++) cm->esc[v][a]  = save[n++];
      break;
    case MP_st:
      for(a = 0; a < Kp*Kp; a++) cm->oesc[v][a]  = save[n++];
      for(a = 0; a < K*K;   a++) cm->esc[v][a]   = save[n++];
      for(a = 0; a < Kp;    a++) cm->lmesc[v][a] = save[n++];
      for(a = 0; a < Kp;    a++) cm->rmesc[v][a] = save[n++];
      break;
    }
  }
  cm->el_selfsc = save[n++];
  free(save);
}

/* Function: DispatchSqAlignment()
 * Date:     EPN, Thu Jan 12 14:47:26 2012
 *
 * Purpose:  Given a CM and a sequence, align the sequence(s) using
 *           the appropriate alignment function and return relevant
 *           data for eventual output in <ret_data>. 
 *
 *           This function can be called from either an alignment
 *           pipeline (i.e. cmalign) or a search/scan pipeline
 *           (i.e. cmsearch or cmscan). <idx> is the (overloaded) flag
 *           for determining which, if -1, we're a search/scan
 *           pipeline. This is only relevant because in a search/scan
 *           pipeline we don't care about determining spos/epos so we
 *           don't call ParsetreeToCMBounds().
 *                        
 *           If (cm->flags & CM_ALIGN_XTAU) we'll potentially tighten
 *           HMM bands until the required DP matrices are below out
 *           limit (<mxsize>). cm->maxtau is the max allowed tau value
 *           during this iterative band tightening, and cm->xtau is
 *           the factor by which we multiply cm->tau at each iteration
 *           during band tightening.
 *
 *           Score centering: in cmalign (<idx> != -1), with HMM
 *           bands calculated here, every residue emission score is
 *           lowered by c = (CP9 Forward score) / L for the CM DP, and
 *           c*L is added back to the score it returns. Every
 *           alignment of the sequence emits all L residues (in
 *           truncated alignment too, and residues emitted by EL
 *           states get the shift), so this lowers every parse's score
 *           by exactly c*L: the optimal parses and posteriors are the
 *           same in exact arithmetic. But the float partial sums of
 *           the DP stay near 0 instead of growing to the size of the
 *           total score, which on long sequences (tens of kb and up)
 *           otherwise makes the reported score drift by bits to
 *           thousands of bits and moves some residues (brief
 *           26_0821-164). The scores in the CM are restored before
 *           returning. Not done in a search pipeline (<idx> == -1,
 *           whose hits are bounded by W), or when no HMM bands are
 *           calculated here (--nonbanded, --small, --qdb, or
 *           <cp9b_valid>).
 *
 * Args:     cm         - the covariance model
 *           errbuf     - char buffer for reporting errors
 *           sq         - sequence to align
 *           idx        - index of sequence (may be used to reorder data later)
 *           mxsize     - max size in Mb of allowable DP mx 
 *           mode       - preset mode of alignment (TRMODE_UNKNOWN if unknown)
 *           pass_idx   - pipeline pass index, determines trunc penalty
 *           cp9b_valid - TRUE if cm->cp9b are valid, don't compute HMM bands
 *           w          - stopwatch for timing individual stages, can be NULL
 *           w_tot      - stopwatch for timing total time per seq, can be NULL
 *           r          - RNG, req'd if CM_ALIGN_SAMPLE, can be NULL otherwise
 *           elr        - (dispatch_sq_alignment_core() only) local-mode EL
 *                        retry state, see ELR_BAND; NULL for a plain alignment.
 *           ret_data   - RETURN: newly created CM_ALNDATA object
 *
 *           The work is done by dispatch_sq_alignment_core(), which has this
 *           interface plus <elr>. DispatchSqAlignment(), at the end of this
 *           file, wraps it (briefs 26_0821-128/26_0821-132).
 *
 * Returns:  eslOK on success;
 *           eslEINCOMPAT on contract violation, errbuf is filled;
 *           eslEMEM if we run out of memory;
 *           <ret_data> is alloc'ed and filled.
 */
static int
dispatch_sq_alignment_core(CM_t *cm, char *errbuf, ESL_SQ *sq, int64_t idx, float mxsize, char mode, int pass_idx,
			   int cp9b_valid, ESL_STOPWATCH *w, ESL_STOPWATCH *w_tot, ESL_RANDOMNESS *r,
			   ELR_BAND *elr, CM_ALNDATA **ret_data)
{
  int           status;            /* easel status */
  CM_ALNDATA   *data         = NULL; /* CM_ALNDATA we'll create and fill */
  float         sc           = 0.;   /* score from alignment function */
  float         pp           = 0.;   /* average PP from alignment function */
  Parsetree_t  *tr           = NULL; /* ptr to a parsetree */
  char         *ppstr        = NULL; /* ptr to a PP string */
  float         secs_bands   = 0.;   /* seconds elapsed for band calculation */
  float         secs_aln     = 0.;   /* seconds elapsed for alignment calculation */
  float         mb_tot       = 0.;   /* size of all DP matrices used for alignment */
  double        tau          = -1.;  /* tau used for calculating bands */
  float         thresh1      = -1.;  /* cp9b->thresh1 used for calculating bands */
  float         thresh2      = -1.;  /* cp9b->thresh2 used for calculating bands */
  int           spos         = -1;   /* start posn: first non-gap CM consensus position */
  int           epos         = -1;   /* end   posn: final non-gap CM consensus position */
  double        save_tau     = cm->tau; /* cm->tau upon entrance, we restore before leaving */
  float         save_thresh1 = (cm->cp9b == NULL) ? -1. : cm->cp9b->thresh1;
  float         save_thresh2 = (cm->cp9b == NULL) ? -1. : cm->cp9b->thresh2;
  float         center_c     = 0.;   /* per-residue score shift for the CM DP (see Purpose) */
  float        *center_save  = NULL; /* saved emission scores while centered, else NULL */
  CM_t         *center_cm    = NULL; /* the CM whose scores are centered */

  /* alignment options */
  int do_nonbanded = (cm->align_opts & CM_ALIGN_NONBANDED) ? TRUE  : FALSE;
  int do_qdb       = (cm->align_opts & CM_ALIGN_QDB)       ? TRUE  : FALSE;
  int do_hbanded   = (do_nonbanded || do_qdb)              ? FALSE : TRUE;
  int do_optacc    = (cm->align_opts & CM_ALIGN_OPTACC)    ? TRUE  : FALSE;
  int do_sample    = (cm->align_opts & CM_ALIGN_SAMPLE)    ? TRUE  : FALSE;
  int do_post      = (cm->align_opts & CM_ALIGN_POST)      ? TRUE  : FALSE;
  int do_small     = (cm->align_opts & CM_ALIGN_SMALL)     ? TRUE  : FALSE;
  int do_trunc     = (cm->align_opts & CM_ALIGN_TRUNC)     ? TRUE  : FALSE;
  int do_xtau      = (cm->align_opts & CM_ALIGN_XTAU)      ? TRUE  : FALSE;
  int do_p7band    = (cm->align_opts & CM_ALIGN_P7BANDED)  ? TRUE  : FALSE;
  /* brief 26_0821-100: TRUE when the non-truncated CP9 that derives bands from a
   * p7 band can end a parse only at node M (cmalign -g --notrunc). Every p7 band
   * is then widened at row L to reach M before use; see cm_p7_bands_ReachModelEnd(). */
  int p7b_end_at_M = ((! do_trunc) && cm->cp9 != NULL &&
                      (! (cm->cp9->flags & (CPLAN9_LOCAL_END | CPLAN9_EL)))) ? TRUE : FALSE;
  /* brief 26_0430-269: --mxsize auto-escalation. When enabled (default, CM_ALIGN_MXESC
   * set, cleared by --no-mxesc) AND the user did NOT force an engine, pick the cheapest
   * engine whose estimated CM-DP peak fits --mxsize: tier (a) standard free-OptAcc,
   * else (b) checkpointed sqrt(M) pinned-OptAcc, else (c) the CYK floor (D&C).  Applies
   * only to the HB free-OptAcc path (--ckpt/--small/--nonbanded/--sample each keep
   * their own engine).  See the tier-selection block just before the HB align dispatch. */
  int do_mxesc     = ((cm->align_opts & CM_ALIGN_MXESC)    &&
                      (! (cm->align_opts & CM_ALIGN_CHECKPT)) &&
                      do_optacc && (! do_sample) && (! do_small) &&
                      (! do_nonbanded) && (! do_qdb)) ? TRUE : FALSE;
  int mxesc_tier   = 0;    /* 0=none/not-decided, 'a'/'b'/'c' once decided (per seq) */
  int eff_checkpt  = (cm->align_opts & CM_ALIGN_CHECKPT) ? TRUE : FALSE; /* effective ckpt engine choice; do_mxesc may raise it per seq */
  /* brief 26_0821-014: WILL the checkpointed engine actually run for this sequence?
   * The four Qualifies() predicates read only cm->sttype[], never cm->cp9b, so this
   * is answerable BEFORE band derivation -- which is the point: two places downstream
   * size against the full HMM-banded cube that a checkpointed run never allocates.
   *   (i)  cp9_IterateSeq2Bands()'s tau ratchet, which otherwise tightens bands (and
   *        so CHANGES THE ALIGNMENT) to fit a matrix that is never built; and
   *   (ii) the pre-dispatch cm_[Tr]AlignSizeNeededHB() eslERANGE gate below, which
   *        otherwise admits the checkpointed engine only when the full free-OptAcc
   *        cube would have fit -- the very gate checkpointing exists to avoid.
   * Mirrors the mxesc tier selector's ckpt_avail OR-of-qualifiers and the engine
   * gates at the two HB dispatch sites, so it is TRUE exactly when one of
   * do_checkpt / do_checkpt_r3 / do_trckpt / do_trckpt_r4 will be.  Note eff_checkpt
   * is still just the --ckpt flag here (do_mxesc raises it only AFTER bands exist),
   * so no non---ckpt run can reach either change. */
  char ckpt_est_mode = (mode == TRMODE_J || mode == TRMODE_L ||
                        mode == TRMODE_R || mode == TRMODE_T) ? mode : TRMODE_T; /* max-plane => safe over-estimate */
  int  ckpt_will_run = FALSE;
  int p7b_iterate_ran = FALSE; /* TRUE once cp9_IterateSeq2BandsP7B() ran for this seq (bands valid even on eslERANGE) */
  int endhull_tried   = FALSE; /* TRUE once the end-band --mxsize fallback ran for this seq (brief 26_0821-103) */
  int doing_search = FALSE;
  /* Brief 26_0430-120: IBV HMM-divergence fallback. Set when cm_TrAlignHB / cm_AlignHB
   * fails on IBV-derived bands and we've already rebuilt with vitband for
   * this sequence; prevents infinite retry.
   */
  int ibv_fallback_used = FALSE;
  /* local-mode EL retry with native CP9 bands (briefs 26_0821-128/26_0821-132) */
  if(elr != NULL && elr->mode == ELR_MODE_NATIVE) do_p7band = FALSE;

#if eslDEBUGLEVEL >= 1
  printf("#DEBUG: in DispatchSqAlignment() %s\n", sq->name);
  printf("#DEBUG: \tdo_nonbanded: %d\n", do_nonbanded);
  printf("#DEBUG: \tdo_optacc:    %d\n", do_optacc);
  printf("#DEBUG: \tdo_sample:    %d\n", do_sample);
  printf("#DEBUG: \tdo_post:      %d\n", do_post);
  printf("#DEBUG: \tdo_small:     %d\n", do_small);
  printf("#DEBUG: \tdo_trunc:     %d\n", do_trunc);
  printf("#DEBUG: \tdo_qdb:       %d\n", do_qdb);
  printf("#DEBUG: \tdoing_search: %d\n", doing_search);
#endif
  
  /* contract check */
  if(do_small  && do_hbanded)       ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do small and HMM banded alignment");
  if(do_small  && do_optacc)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do small and opt acc alignment");
  if(do_post   && do_small)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do PP and small alignment");
  if(do_optacc && do_sample)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to sample and do optacc alignment");
  if(do_sample && r == NULL)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to sample but RNG r == NULL");
  if(do_xtau   && ! do_hbanded)     ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to multiply tau without HMM banded alignment");
  if(do_xtau   && cp9b_valid)       ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to multiply tau but HMM bands already valid");
  if(do_qdb    && do_nonbanded)     ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do qdb and nonbanded alignment");
  if(do_qdb    && do_trunc)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs and truncated alignment");
  /* qdb + trunc combo disallowed only b/c no function exists for it yet */
  if(do_qdb    && (! do_small))     ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs but not divide and conquer");
  /* qdb + small combo disallowed b/c only non-HMM banded non-small alignment functions are not set up to use QDBs */
  if(do_qdb && cm->qdbinfo == NULL) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs but cm->qdbinfo is NULL");
  }
  if(do_qdb && (cm->qdbinfo->dmin2 == NULL || cm->qdbinfo->dmax2 == NULL)) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs but cm->qdbinfo is NULL");
  }
  if(do_trunc && (! cm_pli_PassAllowsTruncation(pass_idx))) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do truncated alignment, but pass_idx doesn't allow truncation (PLI_PASS_STD_ANY)");
  }
  if(pass_idx == PLI_PASS_STD_ANY && (mode == TRMODE_L || mode == TRMODE_R || mode == TRMODE_T)) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() mode is L, R, or T, but pass_idx is PLI_PASS_STD_ANY");
  }

  /* brief 26_0821-014: resolve ckpt_will_run now that every engine-selecting option
   * boolean above is known (see its declaration for why this matters). */
  ckpt_will_run = (eff_checkpt && do_optacc && (! do_sample) && do_hbanded &&
                   (! do_small) &&
                   (do_trunc ? (cm_CheckptTrAlignHB_Qualifies(cm) || cm_CheckptTrOptAccAlignHB_Qualifies(cm))
                             : (cm_CheckptAlignHB_Qualifies(cm)   || cm_CheckptOptAccAlignHB_Qualifies(cm)))) ? TRUE : FALSE;

  if(w_tot != NULL) esl_stopwatch_Start(w_tot);

  if(w != NULL) esl_stopwatch_Start(w);
  /* do small D&C alignment, if nec */
  if(do_small) { 
    if(do_trunc) { 
      sc = TrCYK_DnC(cm, sq->dsq, sq->L, 0, 1, sq->L, pass_idx, FALSE, &tr); /* FALSE: don't reproduce 1.0 behavior */
      mb_tot = 4. * CYKNonQDBSmallMbNeeded(cm, sq->L); /* not sure how accurate this is */
    }
    else { 
      /* with QDB, always use dmin2/dmax2, the looser of the two sets of QDBs in cm->qdbinfo */
      sc = CYKDivideAndConquer(cm, sq->dsq, sq->L, 0, 1, sq->L, &tr, 
			       (do_qdb) ? cm->qdbinfo->dmin2 : NULL, 
			       (do_qdb) ? cm->qdbinfo->dmax2 : NULL);
      mb_tot = CYKNonQDBSmallMbNeeded(cm, sq->L);
    }
  }
  else { /* do_small is FALSE */
    if(do_nonbanded || do_qdb) { /* do not use HMM bands */
      if(do_trunc) { 
	if((status = cm_TrAlignSizeNeeded(cm, errbuf, sq->L, mxsize, do_sample, do_post, 
					  NULL, NULL, NULL, &mb_tot)) != eslOK) goto ERROR;
	if((status = cm_TrAlign(cm, errbuf, sq->dsq, sq->L, mxsize, mode, pass_idx, 
				do_optacc, do_sample, cm->trnb_mx, cm->trnb_shmx, cm->trnb_omx, 
				cm->trnb_emx, r, do_post ? &ppstr : NULL, &tr, NULL, &pp, &sc)) != eslOK) goto ERROR;
      }
      else {
	if((status = cm_AlignSizeNeeded(cm, errbuf, sq->L, mxsize, do_sample, do_post, 
					NULL, NULL, NULL, &mb_tot)) != eslOK) goto ERROR;
	if((status = cm_Align(cm, errbuf, sq->dsq, sq->L, mxsize, do_optacc, do_sample, cm->nb_mx, cm->nb_shmx, 
			      cm->nb_omx, cm->nb_emx, r, do_post ? &ppstr : NULL, &tr, &pp, &sc)) != eslOK) goto ERROR;
      }
    }
    else { /* use HMM bands */
      /* _p7b_kind names the p7 band deriver that ran (NULL unless the
       * do_p7band branch below runs); used by #P7BAND_TIME and by the
       * no-parse fallback. */
      const char *_p7b_kind = NULL;
      if(! cp9b_valid) {
	if(cm->cp9b != NULL) cm->cp9b->hmm_sc = -eslINFINITY; /* set again by the band calculation, for centering */
	if(do_p7band && cm->fp7 != NULL) {
	  /* p7-derived bands: Viterbi trace -> kmin/kmax -> banded CP9 F/B -> CM bands */
	  P7_PROFILE *gm_p7b  = NULL;
	  P7_GMX     *gx_p7b  = NULL;
	  P7_BG      *bg_p7b  = NULL;
	  P7_TRACE   *tr_p7b  = NULL;
	  int        *p7_kmin  = NULL;
	  int        *p7_kmax  = NULL;
	  int        *p7_i2k   = NULL;
	  int         p7_ncells = 0;

	  /* Create P7 objects */
	  bg_p7b = p7_bg_Create(cm->abc);
	  gm_p7b = p7_profile_Create(cm->fp7->M, cm->abc);
	  /* Configure profile for band derivation:
	   * - Standard (non-truncated) alignment: GLOCAL (must align full model)
	   * - Truncated alignment: T-profile (LOCAL + forced full-sequence parse).
	   *   The T-profile allows entry/exit at any model node but forces the parse
	   *   to include the first and last residue (N->N and C->C set to -inf).
	   */
	  if(do_trunc) {
	    p7_ProfileConfig(cm->fp7, bg_p7b, gm_p7b, sq->L, p7_LOCAL);
	    p7_ProfileConfig5PrimeAnd3PrimeTrunc(gm_p7b, sq->L);
	  } else {
	    p7_ProfileConfig(cm->fp7, bg_p7b, gm_p7b, sq->L, p7_GLOCAL);
	  }
	  /* gx_p7b (full O(M*L) p7 matrix) is allocated lazily only where the
	   * unbanded p7_Seq2BandsVit path actually needs it (brief 26_0430-094),
	   * to avoid an eager full-matrix alloc that overflows int32 in
	   * p7_gmx_Create at M=L ~ 1.5e5 (e.g. HSV). */
	  tr_p7b = p7_trace_Create();

	  /* Build local nodepad copy with p7bpad (p7padplus) added */
	  int *local_nodepad = NULL;
	  if(cm->flags & CMH_P7NODEPAD) {
	    int k;
	    ESL_ALLOC(local_nodepad, sizeof(int) * (cm->fp7->M + 1));
	    for(k = 0; k <= cm->fp7->M; k++) local_nodepad[k] = cm->p7_cm_nodepad[k] + cm->p7bpad;
	  }

	  /* Derive p7 bands: IBV (--p7ibv), k-mer chain (--p7kmerchain), or
	   * full unbanded p7_GViterbi trace. */
	  struct timespec _ta_p7b, _tb_p7b;
	  clock_gettime(CLOCK_MONOTONIC, &_ta_p7b);
	  if (cm->p7_use_ibv) {
	    _p7b_kind = "p7ibv";
	    /* F+B direct-band band derivation (brief 26_0430-120). Does NOT take gm/gx
	     * because it extracts transitions directly from cm->fp7. We still
	     * built gm/gx above for the vitband fallback path.
	     *
	     * brief 26_0430-124: with --p7ibv-mem, dispatch to the divide-and-conquer
	     * O(M*logL) band deriver, byte-identical to the flat path but with
	     * dramatically lower peak memory at large M/L.
	     */
	    if (cm->p7_ibv_mem) {
	      _p7b_kind = "p7ibv-dnc";
	      status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, sq->L,
					   cm->p7_ibv_delta, cm->p7_ibv_base_slab,
					   TRUE, /* do_boundary_widen: CM-side preserves current behavior */
					   do_trunc, /* brief 26_0430-171 */
					   cm->p7_ibv_mode, cm->p7_ibv_width, /* brief 26_0430-140 */
					   &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    } else {
	      status = p7_Seq2BandsIBV(cm, errbuf, sq->dsq, sq->L,
				       cm->p7_ibv_delta, do_trunc, /* brief 26_0430-171 */
				       cm->p7_ibv_mode, cm->p7_ibv_width, /* brief 26_0430-140 */
				       &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    }
	    /* No internal ncells==0 fallback here: IBV always produces a band.
	     * An empty row is an error inside the kernel (brief 26_0821-144);
	     * stop on it here, before the fallback below treats any failure as
	     * "bands too wide" and replaces the cause.
	     */
	    if(status == eslEINCONCEIVABLE) goto ERROR;
	  } else if (cm->p7_use_kmerchain) {
	    /* Brief 26_0628-027: genome-scale k-mer seed-and-chain. Collect all seeds
	     * genome-wide, chain by global colinearity, emit multi-segment pins
	     * into the unmodified p7_pins2bands_nodepad. Opt-in. */
	    _p7b_kind = "kmerchain";
	    if (elr != NULL && elr->mode == ELR_MODE_EDIT) {
	      /* local-mode EL retry (briefs 26_0821-128/26_0821-132): the first
	       * pass's band with kmin lowered on some rows; the chain is not rebuilt */
	      status = elr_band_edited(elr, sq->L, &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    } else {
	      int *p7_raw = NULL;
	      status = p7_Seq2BandsKmerChain(cm, errbuf, sq->dsq, sq->L, local_nodepad,
	                                     do_trunc, /* brief 26_0628-033 */
	                                     &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells,
	                                     NULL, NULL, NULL,
	                                     (elr != NULL && elr->capture) ? &p7_raw : NULL);
	      /* first pass of a local alignment: keep the band and the raw pins for
	       * DispatchSqAlignment()'s EL check (briefs 26_0821-128/26_0821-132) */
	      if (status == eslOK && p7_raw != NULL) status = elr_band_keep(elr, sq->L, p7_raw, p7_i2k, p7_kmin, p7_kmax);
	    }
	    /* ncells==0 => M-gate/N-gate fired, or no usable chain. Brief 26_0628-047:
	     * default fallback is now --p7ibv's D&C deriver (this file's own
	     * --p7ibv-mem branch above, same defaults/params), instead of a
	     * Vit-trace band; --p7kmerchain-fbvit reverts to
	     * the old p7_Seq2BandsVit fallback. */
	    if (status == eslOK && p7_ncells == 0) {
	      /* brief 26_0430-260: chain=NONE (zero k-mer anchors) fallback is a
	       * mutually-exclusive 3-way selector: native CP9 banding (DEFAULT) |
	       * --p7kmerchain-fbibv (--p7ibv D&C deriver) | --p7kmerchain-fbvit
	       * (Vit-trace band). Native became the default here in brief
	       * 26_0430-256 (originally env-gated via P7KMERCHAIN_NATIVE_FALLBACK,
	       * now CLI-only -- the env gate is gone, per this project's
	       * documented `env VAR=` empty-vs-unset trap). chain=NONE means the
	       * p7 model found zero signal, so falling back to *another* p7
	       * deriver (p7ibv at the narrow struct-default delta) compounds the
	       * same band-coverage blind spot (brief 26_0430-255: a p7-IBV
	       * coverage failure -- the CM prefers the correct register in 60/62,
	       * the band just excludes it). Native CP9 HMM banding (the identical
	       * mechanism plain, non-p7band cmalign / R_native uses) was measured
	       * to Pareto-dominate the p7ibv fallback over the whole 3281-seq
	       * zero-anchor population (brief 26_0430-255 addendum: 62/62
	       * collapse-seqs recover, 0 regress). */
	      if (cm->p7_kmerchain_fallback_ibv) {
		_p7b_kind = "kmerchain->p7ibv";
		status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, sq->L,
					     cm->p7_ibv_delta, cm->p7_ibv_base_slab,
					     TRUE, do_trunc, cm->p7_ibv_mode, cm->p7_ibv_width,
					     &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
		if (status == eslEINCONCEIVABLE) goto ERROR; /* brief 26_0821-144: empty IBV row */
		if (status == eslOK && p7_ncells == 0) {
		  _p7b_kind = "kmerchain->p7ibv->vitband";
		  if (gx_p7b == NULL) gx_p7b = p7_gmx_Create(cm->fp7->M, sq->L);
		  status = p7_Seq2BandsVit(errbuf, gm_p7b, gx_p7b, bg_p7b, tr_p7b,
					   sq->dsq, sq->L, cm->p7bpad, local_nodepad,
					   0, 0, &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
		}
	      } else if (cm->p7_kmerchain_fallback_vit) {
		_p7b_kind = "kmerchain->vitband";
		if (gx_p7b == NULL) gx_p7b = p7_gmx_Create(cm->fp7->M, sq->L);
		status = p7_Seq2BandsVit(errbuf, gm_p7b, gx_p7b, bg_p7b, tr_p7b,
					 sq->dsq, sq->L, cm->p7bpad, local_nodepad,
					 0, 0, &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	      } else {
		_p7b_kind = "kmerchain->native";
		/* status stays eslOK, p7_ncells stays 0 -> falls through to the
		 * native cp9_Seq2Bands fallback below. Do nothing else here. */
	      }
	    }
	  } else {
	    _p7b_kind = "vitband";
	    if (gx_p7b == NULL) gx_p7b = p7_gmx_Create(cm->fp7->M, sq->L);
	    status = p7_Seq2BandsVit(errbuf, gm_p7b, gx_p7b, bg_p7b, tr_p7b,
				     sq->dsq, sq->L, cm->p7bpad,
				     local_nodepad,
				     0, 0, /* hopback=0, vitend=0 */
				     &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	  }
	  clock_gettime(CLOCK_MONOTONIC, &_tb_p7b);
	  {
	    double _p7b_s = (_tb_p7b.tv_sec - _ta_p7b.tv_sec) +
	                    (_tb_p7b.tv_nsec - _ta_p7b.tv_nsec) / 1e9;
	    if(cm_StageVerbose()) fprintf(stderr, "#P7BAND_TIME %s kind=%s L=%d M=%d t=%.6f ncells=%d ibvmode=%d ibvwidth=%d ibvdelta=%d\n",
	            sq->name, _p7b_kind, (int)sq->L, cm->fp7->M, _p7b_s,
	            p7_ncells, cm->p7_ibv_mode, cm->p7_ibv_width, cm->p7_ibv_delta); /* brief 26_0430-140: ncells = band-size discriminator */
	  }

	  /* --devdiag: report the band, or that the deriver failed */
	  if(status == eslOK && p7_ncells > 0) {
	    int dbg_npin = 0, dbg_minpin = sq->L+1, dbg_maxpin = 0;
	    int dbg_i;
	    for(dbg_i = 1; dbg_i <= sq->L; dbg_i++) {
	      if(p7_i2k[dbg_i] != -1) {
		dbg_npin++;
		if(dbg_i < dbg_minpin) dbg_minpin = dbg_i;
		if(dbg_i > dbg_maxpin) dbg_maxpin = dbg_i;
	      }
	    }
	    float dbg_avgbw = (float)p7_ncells / (float)sq->L;
	    if(cm_StageVerbose()) fprintf(stderr, "#P7BAND %s L=%d M=%d npins=%d pin_range=[%d,%d] ncells=%d avg_bw=%.1f\n",
		    sq->name, (int)sq->L, cm->fp7->M, dbg_npin, dbg_minpin, dbg_maxpin, p7_ncells, dbg_avgbw);
	  }
	  else {
	    if(cm_StageVerbose()) fprintf(stderr, "#P7BAND %s L=%d M=%d FAILED (status=%d ncells=%d)\n",
		    sq->name, (int)sq->L, cm->fp7->M, status, p7_ncells);
	  }

	  if(status == eslOK && p7_ncells > 0) {
	    /* Use p7 bands to derive CM bands via p7-banded CP9 F/B with tau-ratcheting */
	    if(p7b_end_at_M) cm_p7_bands_ReachModelEnd(p7_kmax, sq->L, cm->cp9->M); /* brief 26_0821-100 */
	    status = cp9_IterateSeq2BandsP7B(cm, errbuf, sq->dsq, sq->L, p7_kmin, p7_kmax,
					     1, sq->L, pass_idx, mxsize,
					     doing_search, do_sample, do_post,
					     cm->maxtau, 0, 0, NULL);
	    /* brief 26_0430-312: whitelist, not blacklist. eslOK and the maxtau-capped
	     * eslERANGE are the only statuses this call can return that leave cp9b
	     * validly (if widely) populated -- brief 26_0430-311 verified this holds
	     * for all three eslERANGE origins reachable here. Any other status (a
	     * contract-check eslEINCOMPAT, eslEMEM, eslENORESULT, or one invented three
	     * frames down by a future call graph change) means cp9b was NOT left in a
	     * trustworthy state, whether or not we've seen that status fire here
	     * before -- a blacklist would need re-auditing every time the call graph
	     * grows; this doesn't. */
	    p7b_iterate_ran = (status == eslOK || status == eslERANGE) ? TRUE : FALSE;
	  }
	  else {
	    /* Viterbi found no path or error; fall back to standard cp9 bands */
	    status = eslERANGE;
	  }

	  /* Free p7 objects */
	  if(local_nodepad) free(local_nodepad);
	  if(p7_i2k)  free(p7_i2k);
	  if(p7_kmin) free(p7_kmin);
	  if(p7_kmax) free(p7_kmax);
	  p7_trace_Destroy(tr_p7b);
	  p7_gmx_Destroy(gx_p7b);
	  p7_profile_Destroy(gm_p7b);
	  p7_bg_Destroy(bg_p7b);

	  /* brief 26_0430-310 (M3) narrowed by brief 26_0430-312: eslENORESULT means the
	   * p7 band contained no complete parse, so cp9b was never populated (all -1
	   * sentinels) -- not the "bands too wide" condition the block below exists
	   * for, and it does not satisfy p7b_iterate_ran's premise (already excluded
	   * from the whitelist above). At genome scale (cp9fb_Mb > mxsize, below) that
	   * now falls to the preserve-the-cause goto ERROR in the else arm rather than
	   * being laundered into eslOK or having its errbuf overwritten. At sub-genome
	   * scale it now falls through to the real fallback re-derivation below, same
	   * as it would have pre-310 -- M3's original unconditional goto ERROR here
	   * fired at every scale, silently blocking a sub-genome fallback that could
	   * have succeeded. */

	  if(status != eslOK) {
	    /* P7B bands too wide even at maxtau; fall back to standard cp9 band derivation */
	    /* Brief 26_0430-159: the standard fallback below builds a *non-banded* full CP9
	     * F/B matrix, (L+1)*(M+1) cells. At genome scale (M~L~2e5) that alone is
	     * hundreds of GB, far over --mxsize, and unlike the CM DP matrix there is
	     * no banding/tau lever to shrink the CP9 F/B. The p7-banded path already
	     * returned !eslOK (matrix too big), so attempting the even-larger
	     * non-banded CP9 F/B would only OOM (it did: brief 26_0430-157 SIGSEGV via int32
	     * overflow, now int64-correct but still 155 GB -> SIGKILL). If the
	     * non-banded CP9 F/B (fwd+bck) would itself exceed --mxsize, refuse with
	     * eslERANGE up front -- the eslERANGE-before-alloc convention used for
	     * every CM DP matrix (cm_mx.c). Only fires in this already-over-budget
	     * p7band fallback, so plain non-p7band cmalign is unaffected (its CP9 F/B
	     * is intentionally not mxsize-gated, matching cm_*AlignSizeNeededHB). */
	    float cp9fb_Mb = 2.0 * (float) SizeNeededCP9Matrix(sq->L, cm->cp9->M, NULL, NULL);
	    /* brief 26_0821-013: under fixed-tau (CM_ALIGN_MXESC_FIXEDTAU, the default
	     * inside the mxesc path since brief 26_0430-271) eslERANGE from
	     * cp9_IterateSeq2BandsP7B is not a failure -- it is the DESIGNED signal that
	     * step 0's bands are valid but wider than --mxsize, and that mxesc engine
	     * escalation, not band-tightening, is meant to carry the memory.  The
	     * cp9fb_Mb > mxsize test below predates fixed-tau: it gates that
	     * keep-the-p7-bands rescue on whether the standard fallback is AFFORDABLE,
	     * not on whether it is BETTER.  Because cp9fb_Mb depends only on L and M,
	     * raising --mxsize across it silently flips band derivation from the p7
	     * bands to a non-ratcheted cp9_Seq2Bands re-derivation, which can be far
	     * wider (measured: est_std 1720 -> 12608 Mb on a 3.4 Kb-consensus rRNA CM,
	     * tier b -> tier c, alignment -> refusal at --mxsize 398 that succeeded at
	     * 397).  Hoist the rescue out of the affordability test so more --mxsize
	     * never buys a worse deriver.  Scoped to fixed-tau so the ratchet-restored
	     * (--no-mxesc-fixedtau) path is byte-identical to before. */
	    /* brief 26_0821-096: no-parse p7 band at genome scale. eslENORESULT here
	     * means the p7 band held no complete parse (see the M3 comment above).
	     * Below genome scale the standard CP9 re-derivation further down
	     * recovers; at genome scale (cp9fb_Mb > mxsize) that path cannot run,
	     * and this used to be a hard error. Instead, unless the band already
	     * came from IBV, re-derive the p7 band with the
	     * memory-bounded IBV D&C deriver (the same call and settings as
	     * kmerchain's --p7kmerchain-fbibv no-anchor fallback above, still
	     * truncation-aware if do_trunc), and rerun the p7-banded CP9 band
	     * derivation once. A full-matrix p7 band builder is deliberately not
	     * used here. A second no-parse, or a no-parse IBV band to begin with,
	     * is an error naming both attempts. Any other status from the retry
	     * (eslOK, eslERANGE) flows into the same handling below as a first
	     * attempt's would. */
	    if(status == eslENORESULT && cp9fb_Mb > mxsize) {
	      const char *np_src = (_p7b_kind == NULL)                    ? "p7"                       :
	                           (strcmp(_p7b_kind, "kmerchain") == 0)  ? "--p7kmerchain"            :
	                           "default p7 Viterbi-trace";
	      if(cm->p7_use_ibv || (_p7b_kind != NULL && strcmp(_p7b_kind, "kmerchain->p7ibv") == 0)) {
	        ESL_XFAIL(eslENORESULT, errbuf, "no-parse p7 band from --p7ibv (CP9 F/B total impossible); no further fallback");
	      }
	      p7_i2k = p7_kmin = p7_kmax = NULL; /* freed above */
	      p7_ncells = 0;
	      if((status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, sq->L,
					       cm->p7_ibv_delta, cm->p7_ibv_base_slab,
					       TRUE, do_trunc, cm->p7_ibv_mode, cm->p7_ibv_width,
					       &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells)) != eslOK) goto ERROR;
	      if(p7_ncells == 0) ESL_XFAIL(eslENORESULT, errbuf, "--p7ibv-mem no-parse retry produced no p7 band");
	      if(p7b_end_at_M) cm_p7_bands_ReachModelEnd(p7_kmax, sq->L, cm->cp9->M); /* brief 26_0821-100 */
	      status = cp9_IterateSeq2BandsP7B(cm, errbuf, sq->dsq, sq->L, p7_kmin, p7_kmax,
					       1, sq->L, pass_idx, mxsize,
					       doing_search, do_sample, do_post,
					       cm->maxtau, 0, 0, NULL);
	      p7b_iterate_ran = (status == eslOK || status == eslERANGE) ? TRUE : FALSE; /* same whitelist as the first attempt */
	      free(p7_i2k);  p7_i2k  = NULL;
	      free(p7_kmin); p7_kmin = NULL;
	      free(p7_kmax); p7_kmax = NULL;
	      if(status == eslENORESULT) {
	        snprintf(errbuf, eslERRBUFSIZE, "no-parse p7 band from %s and from its --p7ibv-mem retry (CP9 F/B total impossible)", np_src);
	        goto ERROR;
	      }
	    }
	    if(status == eslOK) {
	      ; /* brief 26_0821-096: the --p7ibv-mem no-parse retry above succeeded; keep its bands */
	    }
	    else if(do_mxesc && p7b_iterate_ran && status == eslERANGE &&
	       (cm->align_opts & CM_ALIGN_MXESC_FIXEDTAU)) {
	      errbuf[0] = '\0';
	      status = eslOK; /* keep the valid p7-banded bands; mxesc escalates instead */
	    }
	    else if(cp9fb_Mb > mxsize) {
	      /* brief 26_0430-269: this is the genome-scale abort the framework replaces --
	       * the p7-banded matrix didn't fit --mxsize even at maxtau, and the standard
	       * non-banded CP9 F/B fallback below is itself too big. Pre-269 this ESL_XFAILs.
	       * Under do_mxesc, if the p7-banded iterate actually ran it left fully-valid
	       * (wide, maxtau) bands in cm->cp9b (cp9_IterateSeq2BandsP7BF_chk_multi leaves the
	       * final all-capped step populated); KEEP them and let the per-seq engine
	       * escalation (tier b/c, below) carry the memory instead of aborting. When the
	       * fallback WOULD fit (cp9fb_Mb <= mxsize, sub-genome), we skip this and re-derive
	       * tighter standard bands exactly as pre-269 -- so the framework is a strict
	       * superset of the old behavior, differing only where the old path aborted. */
	      if(do_mxesc && p7b_iterate_ran) {
	        errbuf[0] = '\0';
	        status = eslOK; /* keep the wide p7-banded bands; skip the standard fallback */
	      }
	      else if(status == eslERANGE)
	        ESL_XFAIL(eslERANGE, errbuf,
			"non-banded CP9 F/B band derivation needs %.1f > %.1f Mb limit.\nUse --mxsize, --maxtau or --tau (this seq needs a p7-banded/--ckpt path).",
			cp9fb_Mb, (float) mxsize);
	      else
	        /* brief 26_0430-312: status is a real error (eslEINCOMPAT, eslEMEM,
		 * eslENORESULT, ...), not the "bands too wide" eslERANGE this branch
		 * exists for. cp9_IterateSeq2BandsP7B already left its cause in errbuf;
		 * the ESL_XFAIL above would overwrite it with an unrelated matrix-size
		 * message. Preserve the real cause instead of laundering it. */
	        goto ERROR;
	    }
	    if(status != eslOK) { /* fallback re-derivation (only when cp9fb_Mb <= mxsize, or !do_mxesc kept status=eslERANGE) */
	      if(do_xtau) {
	        if((status = cp9_IterateSeq2Bands(cm, errbuf, sq->dsq, 1, sq->L, pass_idx, mxsize, doing_search, do_sample, do_post, 1,
						ckpt_will_run, ckpt_est_mode, cm->maxtau, NULL)) != eslOK) goto ERROR;
	      }
	      else {
	        if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq,
					 1, sq->L, cm->cp9b, doing_search, pass_idx, 0)) != eslOK) goto ERROR;
	      }
	    }
	  }
	}
	else if(do_xtau) { /* multiply tau (if nec) until required mx is below Mb limit (mxsize) */
	  /* brief 26_0821-114: size the EL deck as it is allocated (banded), not as the
	   * full L x L triangle that refused every local alignment above ~9.2 Kb. */
	  /* brief 26_0821-116: the full-matrix CP9 F/B holds three (L+1) x (M+1)
	   * matrices (forward, backward, posterior); over --mxsize, calculate the
	   * same bands checkpointed. */
	  int use_ckpt_cp9 = FALSE;
	  if((status = plain_cp9_engine(cm, errbuf, sq, mxsize, 3, TRUE, &use_ckpt_cp9)) != eslOK) goto ERROR;
	  if(use_ckpt_cp9) {
	    if((status = cp9_CheckptIterateSeq2BandsELBanded(cm, errbuf, sq->dsq, 1, sq->L, pass_idx, mxsize, doing_search, do_sample, do_post, 1 /*do_iterate*/,
							     ckpt_will_run, ckpt_est_mode, cm->maxtau, NULL)) != eslOK) goto ERROR;
	  }
	  else {
	    if((status = cp9_IterateSeq2BandsELBanded(cm, errbuf, sq->dsq, 1, sq->L, pass_idx, mxsize, doing_search, do_sample, do_post, 1 /*do_iterate*/,
						      ckpt_will_run, ckpt_est_mode, cm->maxtau, NULL)) != eslOK) goto ERROR;
	  }
	}
	else {
	  /* brief 26_0821-116: two (L+1) x (M+1) matrices here (posteriors overwrite backward) */
	  int use_ckpt_cp9 = FALSE;
	  if((status = plain_cp9_engine(cm, errbuf, sq, mxsize, 2, FALSE, &use_ckpt_cp9)) != eslOK) goto ERROR;
	  if(use_ckpt_cp9) {
	    if((status = cp9_CheckptSeq2Bands(cm, errbuf, sq->dsq, 1, sq->L, cm->cp9b, doing_search, pass_idx, 0)) != eslOK) goto ERROR;
	  }
	  else {
	    if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq, 
				       1, sq->L, cm->cp9b, doing_search, pass_idx, 0)) != eslOK) goto ERROR;
	  }
	}
	if(w != NULL) esl_stopwatch_Stop(w);
	secs_bands = (w == NULL) ? 0. : w->elapsed;
	tau     = cm->tau;
	thresh1 = cm->cp9b->thresh1;
	thresh2 = cm->cp9b->thresh2;
	/* note: we don't set these three if cp9b_valid is TRUE */

	/* brief 26_0821-116: band calculation is done; free its CP9 matrices before
	 * the CM DP, so peak memory is the larger of the two stages, not their sum. */
	shrink_cp9_matrices(cm);
      }

      /* Center the CM DP scores (see Purpose), from here through every CM
       * DP engine and fallback below; restored at MXESC_ALN_DONE or ERROR. */
      if(idx != -1 && ! cp9b_valid && (center_c = center_shift(cm, sq->L)) != 0.) {
	if((status = center_emission_scores(cm, center_c, &center_save)) != eslOK) goto ERROR;
	center_cm = cm;
      }

      /* CYK pre-pass: run CYK on CP9 bands, derive tighter per-state bands before Inside/Outside */
      if(cm->p7_use_cykbands) {
	Parsetree_t *_cyk_tr = NULL;
	float        _cyk_sc  = 0.;
	int          _cyk_ok  = FALSE;

	/* brief 26_0430-273: size-conditional D&C-CYK fallback for the standalone --cykbands
	 * pre-pass. Legacy behavior (still exactly what runs under --no-cykbands-dnc): always
	 * use the HMM-banded FULL-MATRIX CYK engine (cm_TrAlignHB/cm_alignT_hb); if it doesn't
	 * fit mxsize those functions return eslERANGE (not cm_Fail()) so _cyk_ok simply stays
	 * FALSE and tightening is silently skipped -- no crash, but also no tightening, at
	 * whatever genome scale the banded full matrix stops fitting.
	 * New default: estimate the full-matrix cost first (cm_[Tr]AlignSizeNeededHB, the
	 * same estimator + mxsize semantics briefs 26_0430-225/26_0430-226/26_0430-269 already use for the main
	 * alignment engine); if it doesn't fit, fall back to HMM-banded D&C-CYK
	 * (CYKDivideAndConquerHB/TrCYKDivideAndConquerHB), gated on the D&C estimate
	 * (cm_[Tr]DnCAlignSizeNeededHB, brief 26_0430-228's banded-EL version, already on this base).
	 * If even the D&C estimate doesn't fit, this degrades exactly like the legacy path:
	 * _cyk_ok stays FALSE, original bands kept, no crash. Budget = full mxsize for both
	 * tiers: the pre-pass matrix is transient (allocated, consumed for the parsetree,
	 * freed before the main Inside/Outside/OA/Posterior alignment), so it is never
	 * concurrent with the main alignment's own matrices -- sequential, not additive. */
	float _cykpp_est_std = 0., _cykpp_est_dnc = 0.;
	int   _cykpp_use_dnc = FALSE; /* FALSE: full-matrix (tier a) */

	if(! cm->p7_cykbands_no_dnc) { /* --no-cykbands-dnc: always attempt the full-matrix engine below */
	  int _cykpp_status;
	  if(do_trunc) _cykpp_status = cm_TrAlignSizeNeededHB(cm, errbuf, sq->L, mxsize, FALSE/*do_sample*/, FALSE/*do_post*/,
	                                                      NULL, NULL, NULL, NULL, &_cykpp_est_std, NULL);
	  else         _cykpp_status = cm_AlignSizeNeededHB  (cm, errbuf, sq->L, mxsize, FALSE/*do_sample*/, FALSE/*do_post*/,
	                                                      NULL, NULL, NULL, NULL, &_cykpp_est_std, NULL);
	  if(_cykpp_status != eslOK) {
	    errbuf[0] = '\0'; /* clear the eslERANGE message; tier (a) simply doesn't fit */
	    float _cykpp_dnc_vjd = 0., _cykpp_dnc_sh = 0.;
	    int _cykpp_dnc_status;
	    if(do_trunc) {
	      /* estimator is ignorant of the per-mode ncand loop below; TRMODE_T requests all
	       * three (J+L+R) planes, the safe upper bound over whichever single mode the
	       * ncand loop actually runs (same conservative choice mxesc tier (c) makes). */
	      _cykpp_dnc_status = cm_TrDnCAlignSizeNeededHB(cm, errbuf, sq->L, TRMODE_T, &_cykpp_dnc_vjd, &_cykpp_dnc_sh, &_cykpp_est_dnc);
	    }
	    else {
	      _cykpp_dnc_status = cm_DnCAlignSizeNeededHB(cm, errbuf, sq->L, &_cykpp_dnc_vjd, &_cykpp_dnc_sh, &_cykpp_est_dnc);
	    }
	    if(_cykpp_dnc_status != eslOK) goto ERROR;
	    if(_cykpp_est_dnc <= mxsize) _cykpp_use_dnc = TRUE; /* D&C (tier b); else neither fits */
	  }
	}

	if(_cykpp_use_dnc) {
	  /* HMM-banded D&C-CYK floor: CYKDivideAndConquerHB()/TrCYKDivideAndConquerHB()
	   * handle every CM type (bifurcated or bps=0), same engines/pattern as the main
	   * alignment's mxesc tier (c) (brief 26_0430-269) just above in this same function.
	   * CYKDivideAndConquerHB() has no status/errbuf return -- it cm_Fail()s internally
	   * on error, so a returned tr is always valid; _cyk_ok is unconditionally TRUE. */
	  if(do_trunc) {
	    char cand[4]; int ncand = 0, _m273;
	    Parsetree_t *tr_best = NULL; float c_cyk = IMPOSSIBLE;
	    if(cm->cp9b->Jvalid[0]) cand[ncand++] = TRMODE_J;
	    if(cm->cp9b->Lvalid[0]) cand[ncand++] = TRMODE_L;
	    if(cm->cp9b->Rvalid[0]) cand[ncand++] = TRMODE_R;
	    if(cm->cp9b->Tvalid[0]) cand[ncand++] = TRMODE_T;
	    for(_m273 = 0; _m273 < ncand; _m273++) {
	      Parsetree_t *tr_m = NULL; char rm = TRMODE_UNKNOWN;
	      float sc_m = TrCYKDivideAndConquerHB(cm, sq->dsq, sq->L, 0, 1, sq->L, pass_idx, cand[_m273], &rm, &tr_m, cm->cp9b);
	      if(sc_m > c_cyk) { c_cyk = sc_m; if(tr_best) FreeParsetree(tr_best); tr_best = tr_m; }
	      else if(tr_m) FreeParsetree(tr_m);
	    }
	    if(tr_best != NULL) { _cyk_tr = tr_best; _cyk_sc = c_cyk; _cyk_ok = TRUE; }
	    /* else: no root-valid truncation mode -- degrade gracefully, _cyk_ok stays FALSE */
	  }
	  else {
	    _cyk_sc = CYKDivideAndConquerHB(cm, sq->dsq, sq->L, 0, 1, sq->L, &_cyk_tr, cm->cp9b);
	    _cyk_ok = TRUE;
	  }
	}
	else if(do_trunc) {
	  char _cyk_mode = TRMODE_UNKNOWN;
	  float _cyk_avgpp = 0.;
	  if(cm_TrAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize,
			  TRMODE_UNKNOWN, pass_idx, FALSE/*do_optacc*/, FALSE/*do_sample*/,
			  cm->trhb_mx, cm->trhb_shmx, NULL/*post_mx*/, NULL/*emit_mx*/,
			  NULL/*r*/, NULL/*ret_ppstr*/,
			  &_cyk_tr, &_cyk_mode, &_cyk_avgpp, &_cyk_sc) == eslOK && _cyk_tr != NULL) {
	    _cyk_ok = TRUE;
	  }
	} else {
	  if(cm_alignT_hb(cm, errbuf, sq->dsq, sq->L, mxsize, FALSE/*do_optacc*/,
			  cm->hb_mx, cm->hb_shmx, NULL/*emit_mx*/,
			  &_cyk_tr, &_cyk_sc) == eslOK && _cyk_tr != NULL) {
	    _cyk_ok = TRUE;
	  }
	}

	if(_cyk_ok) {
	  /* brief 26_0430-243: tighten via the SHARED bounded helper
	   * (cykbands_tighten_bounded: baseline snapshot + never-loosen intersection
	   * + bounded hd-recompute + backstop) -- the EXACT mechanism --ckpt-cykbands
	   * uses (brief 26_0430-238).  This replaces the old raw cm_BandsFromCYKParsetree()
	   * call, whose unvisited-state full-envelope fallback could inflate the pass-2
	   * band area WITHOUT BOUND and outright fail at genome scale (brief 26_0430-242:
	   * band_area_ratio up to 145924x -> 343 GB matrix -> 4 dengue cells ESL_XFAIL'd).
	   * The intersection bounds every state's i/j band to the pre-tighten CP9/kmerchain
	   * baseline that cm->cp9b holds right here, so the result is bounded-by-construction.
	   *   pad      = cm->p7_cykbands_pad (from --cykpad).
	   *   maxratio = 3.0 plain-constant backstop.
	   *              With never-loosen this backstop should never fire.
	   *   preserve_valid = do_trunc: cm_BandsFromCYKParsetree clobbers marginal
	   *              validity to J-only in truncated mode, which pass-2's cm_TrAlignHB
	   *              needs intact (constraint 3; this also fixes a latent validity bug
	   *              in the old raw-call trunc path).
	   * Per-state pad was archived 2026-05-19 (see cm_CYKPerstatePadCompute doc
	   * comment for failure analysis). Production uses uniform pad. */
	  /* on failure the CP9 bands are kept */
	  (void) cykbands_tighten_bounded(cm, errbuf, _cyk_tr, (int) sq->L, pass_idx,
					  do_trunc/*preserve_valid*/, cm->p7_cykbands_pad, 3.0/*maxratio backstop*/,
					  NULL, NULL, NULL, NULL, NULL);
	  FreeParsetree(_cyk_tr);
	}
      }

      /* --dump-bands: write per-(v,j) band TSV before cm_AlignHB/cm_TrAlignHB */
      if(cm->p7_dump_bands_file != NULL && cm->cp9b != NULL) {
        FILE *_dbfp = fopen(cm->p7_dump_bands_file, "w");
        if(_dbfp != NULL) {
          CP9Bands_t *_dbands = cm->cp9b;
          int _dv, _dj, _djp;
          fprintf(_dbfp, "v\tstate_type\tj\thdmin\thdmax\tdwidth\n");
          for(_dv = 0; _dv < cm->M; _dv++) {
            for(_djp = 0; _djp <= _dbands->jmax[_dv] - _dbands->jmin[_dv]; _djp++) {
              _dj = _dbands->jmin[_dv] + _djp;
              int _dmin = hd_min(_dbands, _dv, _djp);
              int _dmax = hd_max(_dbands, _dv, _djp);
              int _dwidth = (_dmax >= _dmin) ? (_dmax - _dmin + 1) : 0;
              fprintf(_dbfp, "%d\t%s\t%d\t%d\t%d\t%d\n",
                      _dv, Statetype(cm->sttype[_dv]), _dj, _dmin, _dmax, _dwidth);
            }
          }
          fclose(_dbfp);
        }
      }

      if(w != NULL) esl_stopwatch_Start(w);
	CM_ALIGN_HB_RETRY:
	  /* brief 26_0430-269: --mxsize engine auto-escalation. Bands (cm->cp9b) are now
	   * derived (possibly maxtau-wide). Pick the cheapest engine whose estimated CM-DP
	   * peak fits --mxsize, per-sequence: (a) standard free-OptAcc, (b) checkpointed
	   * sqrt(M) pinned-OptAcc, (c) HMM-banded D&C-CYK floor. The comparison quantity is the CM DP
	   * matrices only (excludes the CP9 F/B), matching stock --mxsize semantics
	   * (cm_AlignSizeNeededHB:587 compares cmtotmb, not totmb). Estimators are
	   * safe-overestimates (briefs 26_0430-225/26_0430-226), so a "fits" verdict never under-provisions. */
	  if(do_mxesc) {
	    float est_std_cm = 0., est_std_tot = 0.;
	    char  est_trmode = (mode == TRMODE_J || mode == TRMODE_L ||
	                        mode == TRMODE_R || mode == TRMODE_T) ? mode : TRMODE_T; /* max-plane => safe over-estimate for the trunc ckpt/D&C estimators */
	    int   ckpt_avail = do_trunc
	                       ? (cm_CheckptTrAlignHB_Qualifies(cm) || cm_CheckptTrOptAccAlignHB_Qualifies(cm))
	                       : (cm_CheckptAlignHB_Qualifies(cm)   || cm_CheckptOptAccAlignHB_Qualifies(cm));
	    float ck_dp = 0., ck_em = 0., ck_cp9 = 0., ck_tot = 0., ck_cmmb = 0.;
	    float dnc_vjd = 0., dnc_sh = 0., dnc_tot = 0.;

	    /* tier (a): standard free-OptAcc. The estimator's own cmtotmb-vs-mxsize test
	     * (returns eslERANGE if over) is exactly the gate cm_AlignHB()/cm_TrAlignHB()
	     * apply internally, so eslOK here guarantees the tier-(a) engine won't abort. */
	    /* sparse EL: size the EL decks as cm_[tr_]hb_mx_GrowTo() lays them out (banded),
	     * not as the full O(L^2) triangle, which kept every local alignment above
	     * ~7 kb out of tier (a). */
	    if(do_trunc) status = cm_TrAlignSizeNeededHBELBanded(cm, errbuf, sq->L, mxsize, do_sample, do_post,
	                                                 NULL, NULL, NULL, NULL, &est_std_cm, &est_std_tot);
	    else         status = cm_AlignSizeNeededHBELBanded  (cm, errbuf, sq->L, mxsize, do_sample, do_post,
	                                                 NULL, NULL, NULL, NULL, &est_std_cm, &est_std_tot);
	    if(status != eslOK && status != eslERANGE) goto ERROR;
	    if(status == eslOK) { mxesc_tier = 'a'; eff_checkpt = FALSE; mb_tot = est_std_tot; }
	    else {
	      errbuf[0] = '\0'; /* clear the eslERANGE message; tier (a) simply doesn't fit */
	      /* tier (b): checkpointed sqrt(M) pinned-OptAcc. CM-DP budget = ckptdpmb+emxmb
	       * (excludes cp9mxmb; the ckpt engine's only mxsize gate is its emit_mx, which
	       * this includes). Requires a ckpt engine to exist for this (cm,mode). */
	      if(do_trunc) status = cm_CheckptTrAlignSizeNeededHB(cm, errbuf, sq->L, est_trmode, NULL, NULL,
	                                                          &ck_dp, &ck_em, &ck_cp9, &ck_tot);
	      else         status = cm_CheckptAlignSizeNeededHB  (cm, errbuf, sq->L, NULL, NULL,
	                                                          &ck_dp, &ck_em, &ck_cp9, &ck_tot);
	      if(status != eslOK) goto ERROR;
	      ck_cmmb = ck_dp + ck_em;
	      if(ckpt_avail && ck_cmmb <= mxsize) { mxesc_tier = 'b'; eff_checkpt = TRUE; mb_tot = ck_cmmb; }
	      else {
	        /* tier (c): D&C-CYK floor (CYK parse, no PP). Real HMM-banded D&C peak is small
	         * since brief 26_0430-227 (banded EL deck), and brief 26_0430-228's banded-EL estimator (cherry-
	         * picked onto this base) tracks it tightly at genome scale -- where tier (c)
	         * actually fires (tier b handles everything smaller). cm_[Tr]DnCAlignSizeNeededHB
	         * is a safe over-estimate; the truncated estimate uses est_trmode (max planes) as
	         * a conservative upper bound over the per-mode D&C runs the engine will do.
	         * CYKDivideAndConquerHB / TrCYKDivideAndConquerHB handle every CM type (bifurcated
	         * or bps=0), so no bifurcation gate is needed. */
	        if(do_trunc) status = cm_TrDnCAlignSizeNeededHB(cm, errbuf, sq->L, est_trmode, &dnc_vjd, &dnc_sh, &dnc_tot);
	        else         status = cm_DnCAlignSizeNeededHB  (cm, errbuf, sq->L, &dnc_vjd, &dnc_sh, &dnc_tot);
	        if(status != eslOK) goto ERROR;
	        mxesc_tier = 'c'; eff_checkpt = FALSE; mb_tot = dnc_tot;
	      }
	    }
	    if(cm_StageVerbose()) fprintf(stderr, "#MXESC seq=%s L=%d M=%d trunc=%d tier=%c est_std=%.1f est_ckpt=%.1f est_dnc=%.1f mxsize=%.1f ckpt_avail=%d\n",
	            sq->name, (int)sq->L, (cm->fp7 ? cm->fp7->M : cm->clen), do_trunc, mxesc_tier,
	            est_std_cm, ck_cmmb, dnc_tot, (float)mxsize, ckpt_avail);

	    if(mxesc_tier == 'c') {
	      /* brief 26_0821-103: before giving up, re-derive the CM bands once
	       * without the ROOT_IR tail route that band widenings (not the CP9
	       * posterior) added -- the end-band hull of a target that runs far
	       * past the model's 3' end. Runs only where this sequence would
	       * otherwise abort, so every alignment that fits is untouched; never
	       * on search-derived bands. */
	      if(mb_tot > mxsize && ! endhull_tried && ! cp9b_valid && cm->cp9b != NULL) {
	        int eh_changed = FALSE;
	        endhull_tried = TRUE;
	        if((status = cp9_EndHullRetryBands(cm, errbuf, cm->cp9b, sq->L, pass_idx, &eh_changed)) != eslOK) goto ERROR;
	        if(eh_changed) {
	          goto CM_ALIGN_HB_RETRY;
	        }
	      }
	      /* If even the D&C-CYK floor exceeds --mxsize, no engine fits: fail cleanly. */
	      if(mb_tot > mxsize) {
	        ESL_XFAIL(eslERANGE, errbuf,
	                  "no alignment engine fits --mxsize %.0f Mb for %s (L=%d): free-OptAcc %.0f Mb, ckpt-OptAcc %.0f Mb, D&C-CYK floor %.0f Mb. Raise --mxsize (or use --small).",
	                  (float)mxsize, sq->name, (int)sq->L, est_std_cm, ck_cmmb, dnc_tot);
	      }
	      /* Run the HMM-banded D&C-CYK floor -> CYK parsetree; ppstr stays NULL
	       * (Parsetrees2Alignment() tolerates a per-seq NULL PP under global do_post,
	       * cm_parsetree.c:1055). Truncated: pick the argmax over root-valid marginal modes
	       * (the rung-4 D&C pattern); each TrCYKDivideAndConquerHB run is a single-mode D&C
	       * whose peak is bounded by the est_trmode (all-plane) estimate above. */
	      if(do_trunc) {
	        char cand[4]; int ncand = 0, mm269;
	        Parsetree_t *tr_best = NULL; char c_mode = TRMODE_UNKNOWN; float c_cyk = IMPOSSIBLE;
	        if(cm->cp9b->Jvalid[0]) cand[ncand++] = TRMODE_J;
	        if(cm->cp9b->Lvalid[0]) cand[ncand++] = TRMODE_L;
	        if(cm->cp9b->Rvalid[0]) cand[ncand++] = TRMODE_R;
	        if(cm->cp9b->Tvalid[0]) cand[ncand++] = TRMODE_T;
	        for(mm269 = 0; mm269 < ncand; mm269++) {
	          Parsetree_t *tr_m = NULL; char rm = TRMODE_UNKNOWN;
	          float sc_m = TrCYKDivideAndConquerHB(cm, sq->dsq, sq->L, 0, 1, sq->L, pass_idx, cand[mm269], &rm, &tr_m, cm->cp9b);
	          if(sc_m > c_cyk) { c_cyk = sc_m; c_mode = cand[mm269]; if(tr_best) FreeParsetree(tr_best); tr_best = tr_m; }
	          else if(tr_m) FreeParsetree(tr_m);
	        }
	        if(tr_best == NULL) ESL_XFAIL(eslEINCOMPAT, errbuf, "mxesc tier (c): no root-valid truncation mode for %s", sq->name);
	        tr = tr_best; sc = c_cyk; (void) c_mode;
	      }
	      else {
	        /* CYKDivideAndConquerHB() cm_Fail()s internally on error; a returned tr is valid. */
	        sc = CYKDivideAndConquerHB(cm, sq->dsq, sq->L, 0, 1, sq->L, &tr, cm->cp9b);
	      }
	      goto MXESC_ALN_DONE;
	    }
	  }
	  if(do_trunc) {
		/* brief 26_0430-126 merge: route
		 * SizeNeededHB failure to CM_ALIGN_HB_CHECK_FB (IBV vitband fallback)
		 * instead of directly to ERROR, so the brief 26_0430-120 IBV fallback stays live
		 * in the trunc path. For non-IBV runs CHECK_FB falls through to ERROR.
		 * brief 26_0430-269: under do_mxesc the tier selector above already validated
		 * the chosen engine's fit and set mb_tot, so skip this eslERANGE size-gate
		 * (for tier b it would spuriously fail on the full-matrix estimate). */
		if(! do_mxesc && ! ckpt_will_run) {
		  /* brief 26_0821-114: banded EL deck, as cm_tr_hb_mx_GrowTo() allocates it */
		  status = cm_TrAlignSizeNeededHBELBanded(cm, errbuf, sq->L, mxsize, do_sample, do_post,
					    NULL, NULL, NULL, NULL, NULL, &mb_tot);
		  if(status != eslOK) goto CM_ALIGN_HB_CHECK_FB;
		}
		/* checkpointed sqrt(M)-memory TRUNCATED OptAcc path: engaged by --ckpt
		 * (CM_ALIGN_CHECKPT) for the pure-MATL-chain (bps=0) OptAcc case it
		 * supports (marginal modes J/L/R, T absent), in either local (default) or
		 * global (-g) config; stock cm_TrAlignHB() otherwise (byte-identical
		 * output, but full-cube memory). On failure fall through to
		 * CM_ALIGN_HB_CHECK_FB (IBV vitband fallback), not ERROR. */
		int do_trckpt = (eff_checkpt && do_optacc && (! do_sample) &&
				 cm_CheckptTrAlignHB_Qualifies(cm)) ? TRUE : FALSE;
		/* rung-4 checkpointed STRUCTURED (bps>0) TRUNCATED OptAcc pipeline: engaged
		 * by --ckpt for structured CMs in truncated mode (local or global).  Truncated
		 * analogue of the non-truncated rung-3 path below.  SEPARATE qualifier so the
		 * bps=0 truncated gate (cm_CheckptTrAlignHB_Qualifies, which delegates to the
		 * bps=0 cm_CheckptAlignHB_Qualifies) is NOT relaxed: bps=0 truncated stays on
		 * cm_CheckptTrAlignHB; structured truncated --ckpt routes here; structured
		 * truncated WITHOUT --ckpt falls back to stock cm_TrAlignHB. */
		int do_trckpt_r4 = (eff_checkpt && do_optacc && (! do_sample) &&
				    (! do_trckpt) && cm_CheckptTrOptAccAlignHB_Qualifies(cm)) ? TRUE : FALSE;
		if(do_trckpt) {
		  status = cm_CheckptTrAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, mode, pass_idx,
					       cm->trhb_emx, do_post ? &ppstr : NULL, &tr, NULL, &pp, &sc);
		}
		else if(do_trckpt_r4) {
		  /* brief 26_0610-094: rung-4 pass 1 = the sqrt(M) combined-mode (J/L/R/T)
		   * checkpointed CYK sweep, cm_CheckptTrCYKAlignHB() (brief 26_0610-080, R3): the
		   * only sqrt(M) way to resolve the mode for a STRUCTURED CM (an unpinned
		   * sqrt(M) bifurcation Inside hits the brief 26_0610-017 2D-coupling wall). It
		   * resolves the (mode, bkind/kpin/bbmode/blmode/brmode) pins for the
		   * pass-2 checkpointed pinned posterior + OptAcc + traceback below. */
		  int   *bkind = NULL, *kpin = NULL;
		  char  *bbmode = NULL, *blmode = NULL, *brmode = NULL;
		  Parsetree_t *tr_best = NULL;
		  char   r4_mode = TRMODE_UNKNOWN;
		  float  r4_cyk  = IMPOSSIBLE, r4_Z = 0.;
		  ESL_ALLOC(bkind,  sizeof(int)  * cm->M);
		  ESL_ALLOC(kpin,   sizeof(int)  * cm->M);
		  ESL_ALLOC(bbmode, sizeof(char) * cm->M);
		  ESL_ALLOC(blmode, sizeof(char) * cm->M);
		  ESL_ALLOC(brmode, sizeof(char) * cm->M);
		  status = cm_CheckptTrCYKAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, pass_idx, &tr_best, &r4_mode, &r4_cyk);
		  if(status != eslOK) { free(bkind); free(kpin); free(bbmode); free(blmode); free(brmode); goto CM_ALIGN_HB_CHECK_FB; }
		  rung4_trpins_from_cyk(cm, tr_best, bkind, kpin, bbmode, blmode, brmode);
		  /* brief 26_0430-234: opt-in --ckpt-cykbands -- tighten pass-2 bands from
		   * this same (truncated) pass-1 CYK tree before freeing it.  preserve_valid
		   * TRUE: cm_BandsFromCYKParsetree clobbers marginal validity to J-only in
		   * do_trunc mode; save/restore keeps the resolved r4_mode valid (see helper
		   * doc). Only spatial i/j/hd bands tighten. */
		  /* brief 26_0430-271 item 2: --ckpt-cykbands. On failure the bands are kept. */
		  if(cm->align_opts & CM_ALIGN_CKPT_CYKBANDS)
		    (void) ckpt_cykbands_tighten(cm, errbuf, tr_best, (int) sq->L, pass_idx,
						 TRUE/*preserve_valid*/, NULL, NULL);
		  FreeParsetree(tr_best); tr_best = NULL;
		  /* pass 2: checkpointed pinned truncated posterior -> emit_mx, then checkpointed
		   * pinned truncated OptAcc + pinned-tree traceback -> parsetree + PP.  sc = the
		   * truncated Inside score Z (matches stock cm_TrAlignHB's optacc ret_sc=ins_sc).
		   * Identical for both engine choices above -- only pass-1 differs. */
		  status = cm_CheckptTrPostAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, r4_mode, pass_idx, cm->trhb_emx,
						    bkind, kpin, bbmode, blmode, brmode, &r4_Z, NULL);
		  if(status == eslOK)
		    status = cm_CheckptTrOptAccAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, r4_mode, pass_idx, cm->trhb_emx,
							bkind, kpin, bbmode, blmode, brmode, do_post ? &ppstr : NULL, &tr, NULL, &pp, NULL);
		  sc = r4_Z;
		  free(bkind); free(kpin); free(bbmode); free(blmode); free(brmode);
		  if(status != eslOK) goto CM_ALIGN_HB_CHECK_FB;
		}
		else {
	      	  status = cm_TrAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, mode, pass_idx,
				    do_optacc, do_sample, cm->trhb_mx, cm->trhb_shmx, cm->trhb_omx,
				    cm->trhb_emx, r, do_post ? &ppstr : NULL, &tr, NULL, &pp, &sc);
		}
      }
      else {
	/* brief 26_0430-269: skip this eslERANGE size-gate under do_mxesc (tier selector
	 * already validated the chosen engine's fit and set mb_tot). */
	/* brief 26_0821-014: ... and skip it when the CHECKPOINTED engine will run, for
	 * the same reason: this gate measures the free-OptAcc cube that engine never
	 * allocates, so it refused checkpointed runs that fit comfortably. The
	 * checkpointed engines apply their own mxsize gate internally. */
	if(! do_mxesc && ! ckpt_will_run) {
	  /* brief 26_0821-114: banded EL deck, as cm_hb_mx_GrowTo() allocates it */
	  if((status = cm_AlignSizeNeededHBELBanded(cm, errbuf, sq->L, mxsize, do_sample, do_post,
					  NULL, NULL, NULL, NULL, NULL, &mb_tot)) != eslOK) goto CM_ALIGN_HB_CHECK_FB;
	}
	/* checkpointed sqrt(M)-memory OptAcc path: engaged by --ckpt (CM_ALIGN_CHECKPT)
	 * only for the non-truncated, pure-MATL-chain OptAcc case it supports (local
	 * or global config); stock cm_AlignHB() otherwise (byte-identical output, but
	 * full-cube memory). On failure fall through to CM_ALIGN_HB_CHECK_FB (IBV
	 * vitband fallback), not ERROR. */
	int do_checkpt = (eff_checkpt && do_optacc && (! do_sample) &&
			  cm_CheckptAlignHB_Qualifies(cm)) ? TRUE : FALSE;
	/* rung-3 checkpointed STRUCTURED (bps>0) OptAcc pipeline: engaged by --ckpt
	 * for non-truncated structured CMs (local or global).  Pass 1 checkpointed
	 * CYK (brief 26_0610-084, R1-L) supplies the bifurcation k* pins; pass 2 = checkpointed
	 * pinned posterior -> emit_mx -> checkpointed pinned OptAcc + pinned-tree
	 * traceback.  Separate qualifier so the truncated gate is NOT relaxed
	 * (truncated bps>0 falls back to the rung-4 branch above).  No engine
	 * selector here (brief 26_0610-094 scope: rung-4 only) -- carried through
	 * unchanged from ckpt-trcyk-impl. */
	int do_checkpt_r3 = (eff_checkpt && do_optacc && (! do_sample) &&
			     (! do_checkpt) && cm_CheckptOptAccAlignHB_Qualifies(cm)) ? TRUE : FALSE;
	if(do_checkpt) {
	  status = cm_CheckptAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, cm->hb_emx,
				     do_post ? &ppstr : NULL, &tr, &pp, &sc);
	}
	else if(do_checkpt_r3) {
	  int    *kpin   = NULL;
	  Parsetree_t *tr_cyk = NULL;
	  float   r3_Z   = 0.;
	  ESL_ALLOC(kpin, sizeof(int) * cm->M);
	  /* brief 26_0610-105: a bp>0/bif==0 CM (single stem-loop, no multifurcating
	   * junctions -- cm_CheckptOptAccAlignHB_Qualifies() routes it here on its
	   * MP_st/MR_st alone) has no B_st for pass-1's CYK parse to pin.
	   * rung3_kpin_from_cyk() only extracts pins from B_st nodes, so kpin[]
	   * comes back all -1 regardless of what pass-1 computes, and no pass-2
	   * consumer reads kpin[] outside a B_st context either -- pass-1's parse
	   * is discarded unused.  Measured: pass-1's CYK costs ~27-30% of
	   * end-to-end rung-3 time on such CMs, a ~constant fraction from M=40 to
	   * M=1654, local and global (real, wasted compute, not a rounding
	   * artifact).  Skip it: go straight to the all-(-1) pin set pass-2 would
	   * derive from it anyway. */
	  int has_bif = (CMCountStatetype(cm, B_st) > 0);
	  /* pass 1: HMM-banded CYK parse -> k* pins, with the checkpointed
	   * sqrt(M)-memory CYK engine (cm_CheckptCYKAlignHB, brief 26_0610-084 R1-L),
	   * which closed the pass-1 memory floor (briefs 26_0610-040/26_0610-041) so the
	   * whole rung-3 pipeline is checkpointed. The pins feed the pass-2
	   * checkpointed pinned posterior + OptAcc + traceback below via
	   * rung3_kpin_from_cyk. */
	  if(! has_bif) {
	    int kv; for (kv = 0; kv < cm->M; kv++) kpin[kv] = -1;
	    status = eslOK;
	  }
	  else {
	    status = cm_CheckptCYKAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, &tr_cyk, NULL);
	    if(status == eslOK) {
	      rung3_kpin_from_cyk(cm, tr_cyk, kpin);
	      /* brief 26_0430-234: opt-in --ckpt-cykbands -- tighten pass-2 bands from
	       * this same (non-truncated) pass-1 CYK tree before freeing it. No
	       * valid-flag save/restore needed (rung-3 is non-truncated). */
	      /* brief 26_0430-271 item 2: --ckpt-cykbands. On failure the bands are kept. */
	      if(cm->align_opts & CM_ALIGN_CKPT_CYKBANDS)
		(void) ckpt_cykbands_tighten(cm, errbuf, tr_cyk, (int) sq->L, pass_idx,
					     FALSE/*preserve_valid*/, NULL, NULL);
	      FreeParsetree(tr_cyk); tr_cyk = NULL;
	    }
	  }
	  if(status == eslOK) {
	    /* pass 2: checkpointed pinned posterior + checkpointed pinned OptAcc + traceback */
	    status = cm_CheckptPostAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, cm->hb_emx, kpin, &r3_Z);
	    if(status == eslOK)
	      status = cm_CheckptOptAccAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, cm->hb_emx, kpin,
						do_post ? &ppstr : NULL, &tr, &pp, NULL);
	  }
	  sc = r3_Z;
	  free(kpin);
	  if(status != eslOK) goto CM_ALIGN_HB_CHECK_FB;
	}
	else {
	  status = cm_AlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, do_optacc, do_sample, cm->hb_mx, cm->hb_shmx,
			      cm->hb_omx, cm->hb_emx, r, do_post ? &ppstr : NULL, &tr, &pp, &sc);
	}
      }
    CM_ALIGN_HB_CHECK_FB:
      /* Brief 26_0430-120 IBV HMM-divergence fallback: cm_TrInsideAlignHB returns
       * eslEAMBIGUOUS "no valid parsetree found" on the 3/14 brief 26_0430-117 seqs
       * where HMM-Viterbi disagrees with the CM's preferred parse (brief 26_0430-117
       * §5). Re-derive bands using the unbanded p7_Seq2BandsVit path (PAD=20
       * around the p7-Viterbi trace) and retry the CM alignment once.
       */
      if (status != eslOK && cm->p7_use_ibv && !ibv_fallback_used) {
	ibv_fallback_used = TRUE;
	errbuf[0] = '\0';
	status    = eslOK;

	/* Self-contained vitband redo: build temporary p7 profile/gx/bg/trace,
	 * call p7_Seq2BandsVit, push through cp9_IterateSeq2BandsP7B, free.
	 */
	{
	  P7_PROFILE *fb_gm = p7_profile_Create(cm->fp7->M, cm->abc);
	  P7_BG      *fb_bg = p7_bg_Create(cm->abc);
	  P7_GMX     *fb_gx = p7_gmx_Create(cm->fp7->M, sq->L);
	  P7_TRACE   *fb_tr = p7_trace_Create();
	  int        *fb_i2k = NULL, *fb_kmin = NULL, *fb_kmax = NULL;
	  int         fb_ncells = 0;
	  int        *fb_nodepad = NULL;
	  int         fbk;

	  if (do_trunc) {
	    p7_ProfileConfig(cm->fp7, fb_bg, fb_gm, sq->L, p7_LOCAL);
	    p7_ProfileConfig5PrimeAnd3PrimeTrunc(fb_gm, sq->L);
	  } else {
	    p7_ProfileConfig(cm->fp7, fb_bg, fb_gm, sq->L, p7_GLOCAL);
	  }
	  if (cm->flags & CMH_P7NODEPAD) {
	    fb_nodepad = (int *) malloc(sizeof(int) * (cm->fp7->M + 1));
	    if (fb_nodepad) {
	      for (fbk = 0; fbk <= cm->fp7->M; fbk++)
		fb_nodepad[fbk] = cm->p7_cm_nodepad[fbk] + cm->p7bpad;
	    }
	  }
	  status = p7_Seq2BandsVit(errbuf, fb_gm, fb_gx, fb_bg, fb_tr,
				   sq->dsq, sq->L, cm->p7bpad, fb_nodepad,
				   0, 0,
				   &fb_i2k, &fb_kmin, &fb_kmax, &fb_ncells);
	  if (status == eslOK && fb_ncells > 0) {
	    if (p7b_end_at_M) cm_p7_bands_ReachModelEnd(fb_kmax, sq->L, cm->cp9->M); /* brief 26_0821-100 */
	    status = cp9_IterateSeq2BandsP7B(cm, errbuf, sq->dsq, sq->L,
					     fb_kmin, fb_kmax,
					     1, sq->L, pass_idx, mxsize,
					     doing_search, do_sample, do_post,
					     cm->maxtau, 0, 0, NULL);
	  }
	  if (fb_i2k) free(fb_i2k);
	  if (fb_kmin) free(fb_kmin);
	  if (fb_kmax) free(fb_kmax);
	  if (fb_nodepad) free(fb_nodepad);
	  p7_trace_Destroy(fb_tr);
	  p7_gmx_Destroy(fb_gx);
	  p7_profile_Destroy(fb_gm);
	  p7_bg_Destroy(fb_bg);
	}
	/* brief 26_0821-013 follow-up: mirror the :1299 whitelist here. Under fixed-tau
	 * (CM_ALIGN_MXESC_FIXEDTAU, the default inside the mxesc path) eslERANGE from
	 * cp9_IterateSeq2BandsP7B is the DESIGNED signal that step 0's bands are valid
	 * but wider than --mxsize, not a failure -- and the whitelist above
	 * (brief 26_0430-312, verified by 26_0430-311) establishes that eslERANGE is
	 * one of the two statuses that leave cp9b validly, if widely, populated.
	 * Treating it as a hard error here aborts a sequence that mxesc engine
	 * escalation would have aligned. Clear ONLY eslERANGE; every other status
	 * still errors. */
	if (status == eslERANGE && (cm->align_opts & CM_ALIGN_MXESC_FIXEDTAU)) {
	  errbuf[0] = '\0';
	  status = eslOK; /* keep the valid p7-banded bands; mxesc escalates instead */
	}
	if (status != eslOK) goto ERROR;
	goto CM_ALIGN_HB_RETRY;
      }
      if (status != eslOK) goto ERROR;
    }
  }

 MXESC_ALN_DONE: /* brief 26_0430-269: tier-(c) D&C-CYK floor lands here, skipping the free/ckpt OptAcc engine dispatch */

  if(center_save != NULL) {
    restore_emission_scores(center_cm, center_save);
    center_save = NULL;
    sc = (float) ((double) sc + (double) center_c * (double) sq->L);
  }

  if(w != NULL) esl_stopwatch_Stop(w);
  secs_aln = (w == NULL) ? 0. : w->elapsed;

  /* determine start and end points of the parsetree, 
   * but only if we're not in a search/scan pipeline 
   */
  if(idx != -1) { /* we're not in a search/scan pipeline */
    if((status = ParsetreeToCMBounds(cm, tr, TRUE, TRUE, errbuf, NULL, NULL, NULL, NULL, &spos, &epos)) != eslOK) goto ERROR;
  }
  
  /* create and fill data */
  ESL_ALLOC(data, sizeof(CM_ALNDATA));
  data->sq         = sq;
  data->idx        = idx;
  data->tr         = tr;
  data->sc         = sc;
  data->pp         = (do_post)      ? pp     : 0.;
  data->ppstr      = (do_post)      ? ppstr  : NULL;
  data->spos       = spos;
  data->epos       = epos;
  data->secs_bands = (do_nonbanded) ? 0.     : secs_bands;
  data->secs_aln   = secs_aln;
  data->mb_tot     = mb_tot;
  data->tau        = tau;
  data->thresh1    = thresh1;
  data->thresh2    = thresh2;
  if(w_tot != NULL) esl_stopwatch_Stop(w_tot);
  data->secs_tot   = (w_tot == NULL) ? 0. : w_tot->elapsed;

  *ret_data = data;

  cm->tau = save_tau;
  if(cm->cp9b != NULL) { 
    cm->cp9b->thresh1 = save_thresh1;
    cm->cp9b->thresh2 = save_thresh2;
  }
  return eslOK;

 ERROR: 
  if(center_save != NULL) restore_emission_scores(center_cm, center_save);
  cm->tau = save_tau;
  if(cm->cp9b != NULL) { 
    cm->cp9b->thresh1 = save_thresh1;
    cm->cp9b->thresh2 = save_thresh2;
  }
  if(data != NULL) cm_alndata_Destroy(data, FALSE);
  *ret_data = NULL;

  if(status == eslEMEM) ESL_FAIL(status, errbuf, "DispatchSqAlignment(), out of memory");

  return status; 
}

/*****************************************************************
 * Local-mode EL detect-and-retry for --p7kmerchain bands
 * (briefs 26_0821-127/26_0821-128/26_0821-132)
 *
 * The kmerchain p7 band is built without local ends (ELs). In a local
 * alignment (no -g) on a CM with branches it can exclude the optimal EL
 * parse, e.g. when a large insertion follows a branch end, or when the
 * genome is tandem-doubled. Then the parse puts residues that the chain
 * pinned into an EL instead. DispatchSqAlignment() checks for that after
 * a local kmerchain alignment and, when it finds it, realigns once with a
 * band that admits the alternative EL.
 *
 * None of the rules below has a tuned constant: the trigger is a complete
 * run of chain pins (the chain's own unit of evidence, an exact match of
 * at least its seed length), the retry band is read off the CM topology,
 * a retry is kept only if its parse scores higher, and the pre-retry
 * estimate is a sum of model scores.
 *****************************************************************/

/* EL-capable node: the same rule cm_localize() uses for cm->end[]. */
static int
elr_node_has_el(CM_t *cm, int nd)
{
  return ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
           cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
           cm->ndtype[nd] == BEGR_nd) && nd+1 < cm->nodes && cm->ndtype[nd+1] != END_nd) ? TRUE : FALSE;
}

/* The retry column for an EL entered at node <nd>: the left bound (lpos) of
 * the branch that holds <nd>, i.e. the nearest BEGL/BEGR at or before <nd> in
 * node order (a branch is a contiguous node run from its BEG to its END/BIF).
 * A band that reaches only lpos[nd] itself is not enough: on a 32 kb insert
 * test case, kmin = lpos[nd] collapses the p7-banded CP9 posterior onto an
 * EL-to-the-end path, a kmin 3 columns lower inflates the CM bands past every
 * engine, and every kmin 7 or more columns lower recovers the optimum.
 * Covering every EL entry in the branch removes that dependence. On the ROOT
 * spine (no BEG), use lpos[nd]. */
static int
elr_retry_col(CM_t *cm, CMEmitMap_t *emap, int nd)
{
  int b;
  for(b = nd; b > 0; b--)
    if(cm->ndtype[b] == BEGL_nd || cm->ndtype[b] == BEGR_nd) return ESL_MAX(1, emap->lpos[b]);
  return ESL_MAX(1, emap->lpos[nd]);
}

/* One EL of a parse and the chain pins it swallows. */
typedef struct {
  int x, nd, e1, e2;     /* parse-tree index, entry node, EL residue span      */
  int npin, nrun;        /* raw pins / complete raw-pin runs inside [e1,e2]    */
  int ra, ca, rb, cb;    /* first/last pinned row inside the EL, and its column */
} ELR_EL;

/* Collect the ELs of <tr> and the raw pins each swallows. A run is a maximal
 * stretch of rows whose pins step +1 column per row; it counts only if it
 * starts and ends inside the EL (a run that only touches the EL boundary is
 * the parse leaving or rejoining the chain, not a swallowed one). Returns the
 * number of ELs that swallow a complete run (-1 on allocation failure); *ret_el is NULL if <tr> has no EL. */
static int
elr_detect(CM_t *cm, Parsetree_t *tr, const int *raw, int L, ELR_EL **ret_el, int *ret_nel)
{
  int     status;
  ELR_EL *el = NULL;
  int     x, i, j, nel = 0, ntrig = 0;

  *ret_el = NULL; *ret_nel = 0;
  for(x = 0; x < tr->n; x++) if(tr->state[x] == cm->M) nel++;
  if(nel == 0) return 0;
  ESL_ALLOC(el, sizeof(ELR_EL) * nel);
  nel = 0;
  for(x = 0; x < tr->n; x++) {
    if(tr->state[x] != cm->M) continue;
    ELR_EL *e = &el[nel++];
    e->x  = x;
    e->nd = cm->ndidx[tr->state[tr->prv[x]]];
    e->e1 = tr->emitl[x];
    e->e2 = tr->emitr[x];
    e->npin = e->nrun = 0; e->ra = e->ca = e->rb = e->cb = -1;
    for(i = ESL_MAX(1, e->e1); i <= ESL_MIN(L, e->e2); i++) {
      if(raw[i] == -1) continue;
      e->npin++;
      if(e->ra == -1) { e->ra = i; e->ca = raw[i]; }
      e->rb = i; e->cb = raw[i];
      if(raw[i-1] != -1 && raw[i-1] == raw[i]-1) continue;   /* not a run start (raw[0] is -1) */
      for(j = i; j+1 <= L && raw[j+1] != -1 && raw[j+1] == raw[j]+1; j++) ;
      if(j <= e->e2) e->nrun++;                                /* the run ends inside the EL */
    }
    if(e->nrun > 0) ntrig++;
  }
  *ret_el = el; *ret_nel = nel;
  return ntrig;

 ERROR:
  return -1;
}

/* Lower kmin to <K> on rows [lo,hi] of klo[]. */
static void
elr_klo_set(int *klo, int L, int lo, int hi, int K)
{
  int i;
  for(i = ESL_MAX(1, lo); i <= ESL_MIN(L, hi); i++)
    if(klo[i] == 0 || K < klo[i]) klo[i] = K;
}

/* R1, the retry band for one triggered EL: lower kmin on some rows (kmax is
 * unchanged). In the CP9, an EL from CM node u is a state at HMM node lpos[u]
 * that the path enters from M_lpos[u] and leaves to node rpos[u]
 * (CPlan9InitEL()), so on every row the EL covers, the p7 band must reach
 * column lpos[u]; elr_retry_col() gives the column actually used.
 * (a) EL before the swallowed pins: take the EL-capable u with the latest
 *     resume point epos[u] < ca (ties: latest lpos). Its EL ends u's subtree
 *     before column ca, so the parse can resume and follow the chain from ca.
 *     Lower kmin from the first row whose band reaches the retry column
 *     through ra-1 (the EL has to end before the first swallowed pin).
 * (b) EL after the swallowed pins: take the latest (max lpos) EL-capable u
 *     with lpos[u] <= cb that the band does not already admit on the EL's rows
 *     after rb (lpos[u] < their min kmin). Lower kmin from the first row whose
 *     band reaches the retry column through e2.
 * <do_a>/<do_b> say whether to edit klo[] for (a)/(b); either way *ret_ua and
 * *ret_ub return the chosen nodes (-1 if none). Returns the number of rows
 * whose kmin the edit lowers. */
static int
elr_r1_edit(CM_t *cm, CMEmitMap_t *emap, const ELR_EL *e, const int *kmin, const int *kmax, int L, int *klo,
            int do_a, int do_b, int *ret_ua, int *ret_ub, int *ret_Ka, int *ret_Kb)
{
  int nd, i, n = 0;
  int ua = -1, ub = -1, kmin_tail = INT_MAX;

  *ret_ua = *ret_ub = *ret_Ka = *ret_Kb = -1;
  if(e->ra == -1) return 0;
  for(nd = 1; nd < cm->nodes; nd++)
    if(elr_node_has_el(cm, nd) && emap->epos[nd] < e->ca &&
       (ua == -1 || emap->epos[nd] > emap->epos[ua] ||
        (emap->epos[nd] == emap->epos[ua] && emap->lpos[nd] > emap->lpos[ua]))) ua = nd;
  *ret_ua = ua;
  if(ua != -1 && do_a) {
    int lo = 1, Ka = elr_retry_col(cm, emap, ua), na = 0;
    while(lo <= L && kmax[lo] < Ka) lo++;
    for(i = lo; i < e->ra; i++) if(kmin[i] > Ka) na++;
    if(na > 0) { elr_klo_set(klo, L, lo, e->ra - 1, Ka); *ret_Ka = Ka; n += na; }
  }
  for(i = e->rb+1; i <= ESL_MIN(L, e->e2); i++) kmin_tail = ESL_MIN(kmin_tail, kmin[i]);
  if(kmin_tail != INT_MAX) {
    for(nd = 1; nd < cm->nodes; nd++)
      if(elr_node_has_el(cm, nd) && emap->lpos[nd] <= e->cb && emap->lpos[nd] < kmin_tail &&
         (ub == -1 || emap->lpos[nd] > emap->lpos[ub])) ub = nd;
    *ret_ub = ub;
    if(ub != -1 && do_b) {
      int lo = 1, Kb = elr_retry_col(cm, emap, ub), nb = 0;
      while(lo <= L && kmax[lo] < Kb) lo++;
      for(i = lo; i <= e->e2; i++) if(kmin[i] > Kb) nb++;
      if(nb > 0) { elr_klo_set(klo, L, lo, e->e2, Kb); *ret_Kb = Kb; n += nb; }
    }
  }
  return n;
}

/* The pre-retry estimate: before realigning, estimate whether R1's
 * alternative EL can beat the first parse, so a retry that cannot win is not
 * paid for. Every term is a model score:
 *   gain   = the swallowed residues aligned ungapped along the chain's
 *            diagonals (each takes the nearest preceding raw pin's diagonal;
 *            emission plus the p7 M->M transition), from the first pin after
 *            the largest pin-free stretch inside the EL (the insert or the
 *            unanchored copy; chance anchors inside it do not count), and
 *            stopping before the first column the first parse aligns after its
 *            EL; (b): only up to the last such pin at or before lpos[u];
 *   lost   = (a) only: the first parse's own score (emissions and transitions,
 *            not EL cost) on the residues the alternative moves into its EL,
 *            i.e. after the last residue aligned at or before lpos[u] and
 *            before the first parse's EL;
 *   detour = the p7 model's delete-path score over the columns the
 *            alternative must cross with no residue ((a): epos[u]+1 .. the
 *            first recovered pin's column - 1; (b): epos[u]+1 .. the first
 *            EL's resume column);
 *   dEL    = el_selfsc x (alternative EL length - first EL length).
 * delta = gain - lost + detour + dEL. */
static float
elr_col_emit(CM_t *cm, const int *colL, const int *colR, int c, ESL_DSQ sym)
{
  int v;
  if(sym >= cm->abc->K) return 0.;
  if((v = colL[c]) != -1) return (cm->sttype[v] == MP_st) ? cm->lmesc[v][sym] : cm->esc[v][sym];
  if((v = colR[c]) != -1) return (cm->sttype[v] == MP_st) ? cm->rmesc[v][sym] : cm->esc[v][sym];
  return 0.;
}

/* score contribution of parse node x (its emission and its transition out,
 * as ParsetreeScore() sums them), without the cost of entering an EL */
static float
elr_node_sc(CM_t *cm, Parsetree_t *tr, int x, const ESL_DSQ *dsq)
{
  int     v    = tr->state[x], y;
  char    mode = tr->mode[x];
  float   sc   = 0.;
  ESL_DSQ a, b;

  if(v == cm->M || cm->sttype[v] == E_st || cm->sttype[v] == B_st) return 0.;
  if(tr->nxtl[x] != -1) {
    y = tr->state[tr->nxtl[x]];
    if(v == 0)            sc += (! tr->is_std) ? tr->trpenalty : ((cm->flags & CMH_LOCAL_BEGIN) ? cm->beginsc[y] : cm->tsc[v][y - cm->cfirst[v]]);
    else if(y != cm->M)   sc += cm->tsc[v][y - cm->cfirst[v]];   /* an EL's cost is the dEL term */
  }
  if(cm->sttype[v] == MP_st) {
    a = dsq[tr->emitl[x]]; b = dsq[tr->emitr[x]];
    if(mode == TRMODE_J)      sc += (a < cm->abc->K && b < cm->abc->K) ? cm->esc[v][a*cm->abc->K+b] : DegeneratePairScore(cm->abc, cm->esc[v], a, b);
    else if(mode == TRMODE_L) sc += (a < cm->abc->K) ? cm->lmesc[v][a] : 0.;
    else if(mode == TRMODE_R) sc += (b < cm->abc->K) ? cm->rmesc[v][b] : 0.;
  }
  else if((cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) && ModeEmitsLeft(mode)) {
    a = dsq[tr->emitl[x]]; sc += (a < cm->abc->K) ? cm->esc[v][a] : esl_abc_FAvgScore(cm->abc, a, cm->esc[v]);
  }
  else if((cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) && ModeEmitsRight(mode)) {
    b = dsq[tr->emitr[x]]; sc += (b < cm->abc->K) ? cm->esc[v][b] : esl_abc_FAvgScore(cm->abc, b, cm->esc[v]);
  }
  return sc;
}

/* Score of aligning residues [r1,r2] ungapped along the chain's diagonals:
 * each residue takes the diagonal of the nearest raw pin at or before it and
 * is scored as a match (emission + the p7 model's M->M transition) at that
 * column, if the column is increasing and <= cmax. */
static float
elr_ungapped_gain(CM_t *cm, const int *colL, const int *colR, const ESL_DSQ *dsq, const int *raw, int r1, int r2, int cmax)
{
  int   i, d = 0, have = FALSE, c, last = 0;
  float sc = 0.;

  for(i = r1; i <= r2; i++) {
    if(raw[i] != -1) { d = raw[i] - i; have = TRUE; }
    if(! have) continue;
    c = i + d;
    if(c <= last || c < 1 || c > cmax) continue;
    sc += elr_col_emit(cm, colL, colR, c, dsq[i]);
    if(cm->fp7 != NULL && last > 0 && c == last + 1) sc += log2f(ESL_MAX(cm->fp7->t[last][p7H_MM], 1e-30));
    last = c;
  }
  return sc;
}

/* delete-path score of the p7 model over columns [c1,c2] (0 if empty) */
static float
elr_delete_path(CM_t *cm, int c1, int c2)
{
  P7_HMM *h  = cm->fp7;
  float   sc = 0.;
  int     k;

  if(h == NULL || c1 > c2) return 0.;
  c1 = ESL_MAX(1, c1); c2 = ESL_MIN(h->M, c2);
  if(c1 > c2) return 0.;
  sc += log2f(ESL_MAX(h->t[c1-1][p7H_MD], 1e-30));
  for(k = c1; k < c2; k++) sc += log2f(ESL_MAX(h->t[k][p7H_DD], 1e-30));
  if(c2 < h->M) sc += log2f(ESL_MAX(h->t[c2][p7H_DM], 1e-30));
  return sc;
}

/* The estimates for (a) (u = ua) and (b) (u = ub) of one EL; a missing u,
 * or no pin to recover, gives -eslINFINITY. */
static int
elr_estimate(CM_t *cm, CMEmitMap_t *emap, Parsetree_t *tr, const ESL_DSQ *dsq, int L, const int *raw,
             const ELR_EL *e, int ua, int ub, float *ret_da, float *ret_db)
{
  int    status;
  int   *colL = NULL, *colR = NULL, *res2col = NULL;
  int    nd, v, x, i, c, ru, rb2 = -1, ra2 = -1, cmax;
  float  ga, la = 0., wa, ea, gb, wb, eb;
  int    lfirst = e->e2 - e->e1 + 1;

  *ret_da = *ret_db = -eslINFINITY;
  ESL_ALLOC(colL,    sizeof(int) * (cm->clen+2));
  ESL_ALLOC(colR,    sizeof(int) * (cm->clen+2));
  ESL_ALLOC(res2col, sizeof(int) * (L+2));
  for(c = 0; c <= cm->clen+1; c++) colL[c] = colR[c] = -1;
  for(nd = 0; nd < cm->nodes; nd++) {
    v = cm->nodemap[nd];
    if(cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATP_nd) colL[emap->lpos[nd]] = v;
    if(cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == MATP_nd) colR[emap->rpos[nd]] = v;
  }
  for(i = 0; i <= L+1; i++) res2col[i] = -1;
  for(x = 0; x < tr->n; x++) {   /* residue -> the column the first parse aligns it to */
    v = tr->state[x];
    if(v == cm->M) continue;
    nd = cm->ndidx[v];
    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st) res2col[tr->emitl[x]] = emap->lpos[nd];
    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st) res2col[tr->emitr[x]] = emap->rpos[nd];
  }
  /* the swallowed block starts at the first pin after the largest pin-free stretch inside the EL */
  {
    int prev = e->e1 - 1, gap = -1;
    for(i = e->e1; i <= e->e2; i++) if(raw[i] != -1) { if(i - prev - 1 > gap) { gap = i - prev - 1; ra2 = i; } prev = i; }
  }
  /* recovered columns stop before the first column the first parse aligns after its EL */
  cmax = cm->clen;
  for(i = e->e2 + 1; i <= L; i++) if(res2col[i] > 0) cmax = ESL_MIN(cmax, res2col[i] - 1);
  if(ua != -1 && ra2 != -1) {
    for(i = e->e1 - 1; i >= 1; i--) if(res2col[i] > 0 && res2col[i] <= emap->lpos[ua]) break;
    ru = i + 1;
    ga = elr_ungapped_gain(cm, colL, colR, dsq, raw, ra2, e->e2, cmax);
    for(x = 0; x < tr->n; x++) {
      v = tr->state[x];
      if(v == cm->M) continue;
      i = (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) ? tr->emitr[x] : tr->emitl[x];
      if(i >= ru && i < e->e1) la += elr_node_sc(cm, tr, x, dsq);
    }
    wa = elr_delete_path(cm, emap->epos[ua] + 1, raw[ra2] - 1);
    ea = cm->el_selfsc * (float) (((ra2 - 1) - ru + 1) - lfirst);
    *ret_da = ga - la + wa + ea;
  }
  if(ub != -1 && ra2 != -1) {
    for(i = e->e2; i >= ra2; i--) if(raw[i] != -1 && raw[i] <= emap->lpos[ub]) { rb2 = i; break; }
    if(rb2 != -1) {
      gb = elr_ungapped_gain(cm, colL, colR, dsq, raw, ra2, rb2, emap->lpos[ub]);
      wb = elr_delete_path(cm, emap->epos[ub] + 1, emap->epos[e->nd]);
      eb = cm->el_selfsc * (float) (-(rb2 - ra2 + 1));
      *ret_db = gb + wb + eb;
    }
  }
  free(colL); free(colR); free(res2col);
  return eslOK;

 ERROR:
  if(colL)    free(colL);
  if(colR)    free(colR);
  if(res2col) free(res2col);
  return eslEMEM;
}

/* Corroboration of the first parse by the chain: the number of complete
 * raw-pin runs that lie outside every EL of <tr> and that the parse aligns
 * exactly at their pinned columns. The estimate assumes the first parse
 * follows the chain outside its ELs; with no corroborated run it does not
 * apply (the chain and the parse disagree wholesale, e.g. a one-anchor chain
 * on a small RNA, and the band itself is the suspect). */
static int
elr_corroborated_runs(CM_t *cm, CMEmitMap_t *emap, Parsetree_t *tr, const int *raw, int L, int *ret_ncorr)
{
  int  status;
  int *res2col = NULL, *inel = NULL;
  int  x, v, nd, i, j, n = 0, ok;

  ESL_ALLOC(res2col, sizeof(int) * (L+2));
  ESL_ALLOC(inel,    sizeof(int) * (L+2));
  esl_vec_ISet(res2col, L+2, 0);
  esl_vec_ISet(inel,    L+2, 0);
  for(x = 0; x < tr->n; x++) {
    v = tr->state[x];
    if(v == cm->M) { for(i = ESL_MAX(1, tr->emitl[x]); i <= ESL_MIN(L, tr->emitr[x]); i++) inel[i] = 1; continue; }
    nd = cm->ndidx[v];
    if(cm->sttype[v] == MP_st || cm->sttype[v] == ML_st) res2col[tr->emitl[x]] = emap->lpos[nd];
    if(cm->sttype[v] == MP_st || cm->sttype[v] == MR_st) res2col[tr->emitr[x]] = emap->rpos[nd];
  }
  for(i = 1; i <= L; i++) {
    if(raw[i] == -1 || (i > 1 && raw[i-1] != -1 && raw[i-1] == raw[i]-1)) continue;   /* run start */
    for(j = i, ok = TRUE; j <= L && raw[j] != -1 && (j == i || raw[j] == raw[j-1]+1); j++)
      if(inel[j] || res2col[j] != raw[j]) ok = FALSE;
    if(ok) n++;
  }
  free(res2col); free(inel);
  *ret_ncorr = n;
  return eslOK;

 ERROR:
  if(res2col) free(res2col);
  if(inel)    free(inel);
  return eslEMEM;
}

/* T, the tail rule, applied with no parse (on a refusal): the rows after the
 * chain's last raw pin can be absorbed by an EL only if the band holds that
 * EL's HMM node. If the band already holds an EL entry that can absorb the
 * tail (an EL-capable node with kmin_tail <= lpos <= the last pin's column,
 * e.g. near the end of a single-branch CM), the tail is not what the band is
 * missing: no edit. Otherwise take the latest (max lpos) EL-capable u with
 * lpos[u] <= the last pin's column and lpos[u] below the tail rows' kmin, and
 * lower kmin to elr_retry_col(u) from the first row whose band reaches it to
 * L. Returns the number of rows whose kmin the edit lowers. */
static int
elr_tail_edit(CM_t *cm, CMEmitMap_t *emap, const int *raw, const int *kmin, const int *kmax, int L, int *klo, int *ret_K)
{
  int i, nd, rb = -1, cb = -1, ub = -1, kmin_tail = INT_MAX, n = 0;

  *ret_K = -1;
  for(i = L; i >= 1; i--) if(raw[i] != -1) { rb = i; cb = raw[i]; break; }
  if(rb == -1 || rb == L) return 0;
  for(i = rb+1; i <= L; i++) kmin_tail = ESL_MIN(kmin_tail, kmin[i]);
  for(nd = 1; nd < cm->nodes; nd++)
    if(elr_node_has_el(cm, nd) && emap->lpos[nd] >= kmin_tail && emap->lpos[nd] <= cb) return 0;
  for(nd = 1; nd < cm->nodes; nd++)
    if(elr_node_has_el(cm, nd) && emap->lpos[nd] <= cb && emap->lpos[nd] < kmin_tail &&
       (ub == -1 || emap->lpos[nd] > emap->lpos[ub])) ub = nd;
  if(ub == -1) return 0;
  {
    int lo = 1, K = elr_retry_col(cm, emap, ub);
    while(lo <= L && kmax[lo] < K) lo++;
    for(i = lo; i <= L; i++) if(kmin[i] > K) n++;
    if(n > 0) { elr_klo_set(klo, L, lo, L, K); *ret_K = K; }
  }
  return n;
}

/* elr_consider(): a retry has run with status <status2> and result <*data2>.
 * Keep it iff its parse scores higher than the best so far. The comparison
 * is on the reported parse's own score (ParsetreeScore()), never on
 * data->sc: under optimal accuracy data->sc is the Inside score summed over
 * the band, which a different band changes even when the parse is the same.
 * The retry's time is added to <sum> (bands, alignment, total) either way.
 * Returns TRUE if kept. */
static int
elr_consider(CM_t *cm, ESL_SQ *sq, int status2, CM_ALNDATA **data2, CM_ALNDATA **best, float *best_psc, double *sum, float *ret_psc)
{
  char  errbuf2[eslERRBUFSIZE];
  float psc2 = -eslINFINITY;
  int   kept = FALSE;

  if(status2 == eslOK) {
    sum[0] += (*data2)->secs_bands;
    sum[1] += (*data2)->secs_aln;
    sum[2] += (*data2)->secs_tot;
    if(ParsetreeScore(cm, NULL, errbuf2, (*data2)->tr, sq->dsq, FALSE, &psc2, NULL, NULL, NULL, NULL) != eslOK) psc2 = -eslINFINITY;
    if(psc2 > *best_psc) kept = TRUE;
  }
  if(kept) { cm_alndata_Destroy(*best, FALSE); *best = *data2; *best_psc = psc2; }
  else if(*data2 != NULL) cm_alndata_Destroy(*data2, FALSE);
  *data2 = NULL;
  if(ret_psc != NULL) *ret_psc = psc2;
  return kept;
}

/* Function: DispatchSqAlignment()
 *
 * Purpose:  See the header of dispatch_sq_alignment_core() above, which does
 *           the alignment.
 *
 *           Local-mode EL detect-and-retry (briefs 26_0821-127/26_0821-128/
 *           26_0821-132): for a local alignment (CMH_LOCAL_END) in cmalign
 *           (idx != -1), not sampled (--sample), whose p7
 *           band came from --p7kmerchain itself, check
 *           each EL of the parse for a complete run of the chain's raw pins
 *           inside it. If one does, the band probably excluded the better EL
 *           placement, and the sequence is realigned:
 *             R1: the first pass's band with kmin lowered on the rows the
 *                 alternative EL needs (elr_r1_edit()), for each part whose
 *                 pre-retry estimate (elr_estimate()) is positive; the
 *                 estimate is applied only if the chain corroborates the first
 *                 parse outside its ELs (elr_corroborated_runs()), else R1
 *                 edits every part;
 *             R2: if R1 was not kept, native CP9 bands, unless the topology
 *                 offered an alternative and every estimate rejected it.
 *           A retry is kept only if its parse scores higher (elr_consider()).
 *           If the first pass is refused (eslERANGE), the tail rule T
 *           (elr_tail_edit()) gets one retry.
 *           The returned data describe the kept alignment; its times include
 *           every pass. Nothing is printed.
 */
int
DispatchSqAlignment(CM_t *cm, char *errbuf, ESL_SQ *sq, int64_t idx, float mxsize, char mode, int pass_idx,
		    int cp9b_valid, ESL_STOPWATCH *w, ESL_STOPWATCH *w_tot, ESL_RANDOMNESS *r,
		    CM_ALNDATA **ret_data)
{
  int          status, status2;
  int          L     = sq->L;
  ELR_BAND     eb;              /* the first pass's kmerchain band */
  ELR_BAND     rq;              /* a retry's request */
  CMEmitMap_t *emap  = NULL;
  ELR_EL      *el    = NULL;
  int         *klo   = NULL;
  float       *da    = NULL, *db = NULL;
  CM_ALNDATA  *best  = NULL, *data2 = NULL;
  char         errbuf2[eslERRBUFSIZE];
  float        best_psc = -eslINFINITY, psc2;
  double       sum[3]   = { 0., 0., 0. };   /* bands, alignment, total seconds over every pass */
  int          nel = 0, ntrig, k, ua, ub, a, b, nrow;
  int          ncand = 0, npos = 0, ncorr = 0, gate, nkept = 0, nretry = 0;

  /* Not for a sampled alignment (keeping the higher-scoring of two samples
   * would bias the sample). */
  if(idx == -1 || ! (cm->flags & CMH_LOCAL_END) || (cm->align_opts & CM_ALIGN_SAMPLE))
    return dispatch_sq_alignment_core(cm, errbuf, sq, idx, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, NULL, ret_data);

  memset(&eb, 0, sizeof(ELR_BAND));
  memset(&rq, 0, sizeof(ELR_BAND));
  eb.mode    = ELR_MODE_NONE;
  eb.capture = TRUE;
  status = dispatch_sq_alignment_core(cm, errbuf, sq, idx, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, &eb, ret_data);
  if(eb.raw == NULL || eb.L != L) { elr_band_release(&eb); return status; }   /* the band did not come from kmerchain */

  /* a refusal (eslERANGE): one retry with the tail rule T. It is not applied
   * before alignment, since it would change many bands that align fine. */
  if(status == eslERANGE) {
    double t_first = 0.;
    int    K;
    if(w_tot != NULL) { esl_stopwatch_Stop(w_tot); t_first = w_tot->elapsed; }   /* the refused pass stops nothing */
    if((emap = CreateEmitMap(cm)) == NULL) { status = eslEMEM; goto ERROR; }
    ESL_ALLOC(klo, sizeof(int) * (L+1));
    esl_vec_ISet(klo, L+1, 0);
    nrow = elr_tail_edit(cm, emap, eb.raw, eb.kmin, eb.kmax, L, klo, &K);
    if(nrow > 0) {
      rq.mode = ELR_MODE_EDIT; rq.L = L; rq.i2k = eb.i2k; rq.kmin = eb.kmin; rq.kmax = eb.kmax; rq.klo = klo;
      status2 = dispatch_sq_alignment_core(cm, errbuf2, sq, idx, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, &rq, &data2);
      if(status2 == eslOK) {
        data2->secs_tot += t_first;
        best = data2; data2 = NULL;
        status = eslOK; errbuf[0] = '\0';
      }
    }
    goto DONE;
  }
  if(status != eslOK) goto DONE;

  best = *ret_data;
  *ret_data = NULL;
  sum[0] = best->secs_bands;
  sum[1] = best->secs_aln;
  sum[2] = best->secs_tot;
  if((ntrig = elr_detect(cm, best->tr, eb.raw, L, &el, &nel)) < 0) { status = eslEMEM; goto ERROR; }
  if(ntrig == 0) goto DONE;

  if(ParsetreeScore(cm, NULL, errbuf2, best->tr, sq->dsq, FALSE, &best_psc, NULL, NULL, NULL, NULL) != eslOK) best_psc = -eslINFINITY;
  if((emap = CreateEmitMap(cm)) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(klo, sizeof(int) * (L+1));
  esl_vec_ISet(klo, L+1, 0);
  ESL_ALLOC(da, sizeof(float) * nel);
  ESL_ALLOC(db, sizeof(float) * nel);

  /* the pre-retry estimate of each triggered EL's alternatives */
  for(k = 0; k < nel; k++) {
    da[k] = db[k] = -eslINFINITY;
    if(el[k].nrun == 0) continue;
    elr_r1_edit(cm, emap, &el[k], eb.kmin, eb.kmax, L, klo, FALSE, FALSE, &ua, &ub, &a, &b);
    if((status = elr_estimate(cm, emap, best->tr, sq->dsq, L, eb.raw, &el[k], ua, ub, &da[k], &db[k])) != eslOK) goto ERROR;
    if(da[k] > -eslINFINITY || db[k] > -eslINFINITY) ncand++;
    if(da[k] > 0.           || db[k] > 0.)           npos++;
  }
  if((status = elr_corroborated_runs(cm, emap, best->tr, eb.raw, L, &ncorr)) != eslOK) goto ERROR;
  gate = (ncorr > 0) ? TRUE : FALSE;   /* apply the estimate only where its model holds */

  /* R1: the topology band edit */
  esl_vec_ISet(klo, L+1, 0);
  nrow = 0;
  for(k = 0; k < nel; k++) {
    if(el[k].nrun == 0) continue;
    nrow += elr_r1_edit(cm, emap, &el[k], eb.kmin, eb.kmax, L, klo, (! gate) || da[k] > 0., (! gate) || db[k] > 0., &ua, &ub, &a, &b);
  }
  if(nrow > 0) {
    rq.mode = ELR_MODE_EDIT; rq.L = L; rq.i2k = eb.i2k; rq.kmin = eb.kmin; rq.kmax = eb.kmax; rq.klo = klo;
    status2 = dispatch_sq_alignment_core(cm, errbuf2, sq, idx, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, &rq, &data2);
    nretry++;
    if(elr_consider(cm, sq, status2, &data2, &best, &best_psc, sum, &psc2)) nkept++;
  }
  /* R2: native CP9 bands, if R1 was not kept. Skipped when the estimate
   * applies, the topology offered an alternative EL and every estimate
   * rejected it: then the first parse's EL placement stands. */
  if(nkept == 0 && ! (gate && ncand > 0 && npos == 0)) {
    rq.mode = ELR_MODE_NATIVE; rq.L = 0; rq.i2k = rq.kmin = rq.kmax = rq.klo = NULL;
    status2 = dispatch_sq_alignment_core(cm, errbuf2, sq, idx, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, &rq, &data2);
    nretry++;
    if(elr_consider(cm, sq, status2, &data2, &best, &best_psc, sum, &psc2)) nkept++;
  }

 DONE:
  if(best != NULL) {
    if(nretry > 0) { best->secs_bands = sum[0]; best->secs_aln = sum[1]; best->secs_tot = sum[2]; }
    *ret_data = best;
  }
  elr_band_release(&eb);
  if(emap) FreeEmitMap(emap);
  if(el)   free(el);
  if(klo)  free(klo);
  if(da)   free(da);
  if(db)   free(db);
  return status;

 ERROR:
  elr_band_release(&eb);
  if(emap)  FreeEmitMap(emap);
  if(el)    free(el);
  if(klo)   free(klo);
  if(da)    free(da);
  if(db)    free(db);
  if(best)  cm_alndata_Destroy(best, FALSE);
  if(*ret_data != NULL) { cm_alndata_Destroy(*ret_data, FALSE); *ret_data = NULL; }
  ESL_FAIL(status, errbuf, "DispatchSqAlignment(), out of memory");
}

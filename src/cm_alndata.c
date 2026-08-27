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

#include "easel.h"

#include "hmmer.h"

#include "infernal.h"

static int sub_alignment_prep(CM_t *orig_cm, char *errbuf, ESL_SQ *sq, CMSubMap_t **ret_submap, CM_t **ret_sub_cm);

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
 * uses (cm_alndata.c ~L940), factored out so the CKPT_CYKBANDS piggyback can
 * report the pre/post-tighten reduction with identical methodology to brief 233. */
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
 * so the unvisited-state full-envelope fallback (project_cykbands_h5_cache_verdict)
 * could inflate the pass-2 matrix without bound and outright fail at genome scale
 * (brief 26_0430-242: band_area_ratio up to 145924x, 343 GB matrix, 4 dengue
 * cells ESL_XFAIL'd).  The CKPT_CYKBANDS path never blew up because it already
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

  status = cm_BandsFromCYKParsetree(cm, errbuf, tr, 1, L, pad, NULL, FALSE,
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
   * (project_cykbands_h5_cache_verdict / brief 237's measured 94K->190M-cell
   * blowup on the LSU fragment repro). Intersecting bounds the result to
   * [never worse than the untightened baseline]: unvisited states'
   * full-envelope CYK band always contains the CP9 band, so the intersection
   * collapses back to exactly the CP9 band (safe, unchanged cost, and -- per
   * brief 235 -- still covers every state OptAcc/MEA might diverge into);
   * visited states' CYK band already sits inside (or pad-widens slightly past)
   * the CP9 band, so the intersection preserves the brief 234/236 tightening
   * win. This makes brief 237's revert-on-blowup ratio cap a pure backstop
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
         * (brief 238's safety requirement) -- fall back to the CP9
         * (pre-tighten) band for this state alone, whatever it is
         * (sentinel or real). */
        ni_min = s_imin[v]; ni_max = s_imax[v];
        nj_min = s_jmin[v]; nj_max = s_jmax[v];
        if(cp9_is_sentinel) n_cp9_sentinel++;   /* CP9 baseline itself already unreachable (expected, harmless) */
        else                n_real_empty++;     /* real accuracy-risk empty (brief 238; should be rare/never) */
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

  /* brief 26_0430-243 VALIDATION instrumentation (env-gated, diagnostic-only,
   * zero effect on bands/output -- Gate 0 byte-identity is preserved).
   * CYKBANDS_SUBSET_CHECK: assert the never-loosen bounded-by-construction
   *   invariant -- every state's final i/j band must be a SUBSET of the
   *   pre-tighten baseline snapshot (imin>=base, imax<=base, jmin>=base,
   *   jmax<=base). Counts + prints any violation (should be exactly 0).
   * CYKBANDS_DUMP_TIGHT=<file>: dump per-state final tightened band AND the
   *   baseline snapshot band, so the fixed --cykbands path and the CKPT_CYKBANDS
   *   path (both call this helper) can be diffed for cross-path band equivalence. */
  if(getenv("CYKBANDS_SUBSET_CHECK") != NULL) {
    int v, nviol = 0;
    for(v = 0; v < M; v++) {
      if(cm->cp9b->imin[v] < s_imin[v] || cm->cp9b->imax[v] > s_imax[v] ||
         cm->cp9b->jmin[v] < s_jmin[v] || cm->cp9b->jmax[v] > s_jmax[v]) {
        nviol++;
        if(nviol <= 10)
          fprintf(stderr, "#CYKBANDS_SUBSET VIOLATION v=%d tight[%d,%d][%d,%d] base[%d,%d][%d,%d]\n",
                  v, cm->cp9b->imin[v], cm->cp9b->imax[v], cm->cp9b->jmin[v], cm->cp9b->jmax[v],
                  s_imin[v], s_imax[v], s_jmin[v], s_jmax[v]);
      }
    }
    fprintf(stderr, "#CYKBANDS_SUBSET_CHECK M=%d L=%d violations=%d (0 == bounded-by-construction OK)\n",
            cm->M, L, nviol);
  }
  {
    const char *df = getenv("CYKBANDS_DUMP_TIGHT");
    if(df != NULL) {
      FILE *fp = fopen(df, "w");
      if(fp != NULL) {
        int v;
        fprintf(fp, "v\ttimin\ttimax\ttjmin\ttjmax\tbimin\tbimax\tbjmin\tbjmax\n");
        for(v = 0; v < M; v++)
          fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", v,
                  cm->cp9b->imin[v], cm->cp9b->imax[v], cm->cp9b->jmin[v], cm->cp9b->jmax[v],
                  s_imin[v], s_imax[v], s_jmin[v], s_jmax[v]);
        fclose(fp);
      }
    }
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
 * reuses the checkpointed pass-1 tree already in hand.  Off by default; opt-in
 * diagnostic env var CKPT_CYKBANDS (pad override via CKPT_CYKBANDS_PAD).
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
 * inheritance (project_cykbands_h5_cache_verdict) can only widen an unpinned
 * state -- it can never exclude a pin.  No pin/band conflict is possible.
 *
 * brief 26_0430-237: for a low-coverage sequence (most of a large CM's states
 * unvisited by the pass-1 CYK parse) the unvisited-state full-envelope fallback
 * in cm_BandsFromCYKParsetree can INFLATE the band area by 30x-40000x (O(L^2)
 * d-cells per fallback state), producing a hundreds-of-billions-of-cell pass-2
 * matrix that hangs or OOM-kills the process. Guard: after tightening, if the
 * cell count exceeds <maxratio> x the untightened baseline (default 3.0, env
 * CKPT_CYKBANDS_MAXRATIO), REVERT to the snapshotted untightened --ckpt bands.
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
                                                 * (5S_rRNA pad=8, brief 236); low-coverage blowup is
                                                 * 30x-40000x (brief 237); 3.0 cleanly separates them.
                                                 * With never-loosen this is a pure backstop. */
  const char *pad_env  = getenv("CKPT_CYKBANDS_PAD");
  const char *mr_env   = getenv("CKPT_CYKBANDS_MAXRATIO"); /* <=0 disables the guard */
  int64_t     n_cp9_sentinel = 0, n_real_empty = 0;
  double      revert_ratio   = 0.;

  /* CKPT_CYKBANDS-only env overrides (the shared core is env-agnostic; brief 26_0430-243). */
  if(pad_env != NULL) pad      = atoi(pad_env);
  if(mr_env  != NULL) maxratio = atof(mr_env);

  status = cykbands_tighten_bounded(cm, errbuf, tr, L, pass_idx, preserve_valid,
                                    pad, maxratio, ret_orig, ret_tight,
                                    &n_cp9_sentinel, &n_real_empty, &revert_ratio);
  if(status != eslOK) return status;

  if((n_cp9_sentinel > 0 || n_real_empty > 0) && (getenv("INFERNAL_CKPT_VERBOSE") || getenv("CKPT_CYKBANDS_VERBOSE")))
    fprintf(stderr, "#CKPT_CYKBANDS never-loosen: %" PRId64 " state(s) CP9-sentinel (expected, zero-cost) + "
            "%" PRId64 " state(s) real-empty-intersection (accuracy-risk fallback) (M=%d L=%d)\n",
            n_cp9_sentinel, n_real_empty, cm->M, L);
  if(revert_ratio > 0. && (getenv("INFERNAL_CKPT_VERBOSE") || getenv("CKPT_CYKBANDS_VERBOSE")))
    fprintf(stderr, "#CKPT_CYKBANDS blowup-guard REVERTED: tight/orig=%.1f > %.1f, kept untightened bands (M=%d L=%d)\n",
            revert_ratio, maxratio, cm->M, L);
  return eslOK;
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

/* Function: sub_alignment_prep()
 * Date:     EPN, Mon Jan  9 05:25:26 2012
 *
 * Purpose:  Prepare for an alignment workunit in sub-mode.
 *
 * Args:     orig_cm    - the covariance model
 *           errbuf     - char buffer for reporting errors
 *           sq         - the sequence we're creating the sub CM for
 *           ret_submap - RETURN: the sub CM to original CM map, created here
 *           ret_sub_cm - RETURN: the sub CM, created here
 *
 * Returns:  eslOK on success;
 *           eslEMEM if we run out of memory;
 *           eslEINVAL on other error, errbuf is filled;
 *           <ret_dataA> is alloc'ed and filled with sq_block->count CM_ALNDATA objects.
 */
int
sub_alignment_prep(CM_t *orig_cm, char *errbuf, ESL_SQ *sq, CMSubMap_t **ret_submap, CM_t **ret_sub_cm)
{
  int          status;            /* easel status */
  CM_t        *sub_cm  = NULL;    /* the sub CM */
  CMSubMap_t  *submap  = NULL;    /* map from mother CM to sub CM, and vice versa */
  int          spos;              /* HMM node most likely to have emitted posn 1 of target seq */
  int          spos_state;        /* HMM state type for curr spos 0=match or 1=insert */
  int          epos;              /* HMM node most likely to have emitted posn L of target seq */
  int          epos_state;        /* HMM state type for curr epos 0=match or 1=insert */

  /* step 1. predict start and end positions (HMM nodes) from posterior matrix */
  if((status = cp9_Seq2Posteriors(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, sq->dsq, 1, sq->L, 0)) != eslOK) return status; 
  CP9NodeForPosn(orig_cm->cp9, 1, sq->L,     1, orig_cm->cp9_bmx, &spos, &spos_state, 0., TRUE,  0);
  CP9NodeForPosn(orig_cm->cp9, 1, sq->L, sq->L, orig_cm->cp9_bmx, &epos, &epos_state, 0., FALSE, 0);
  /* Deal with special cases for sub-CM alignment: If the most
   * likely state to have emitted the first or last residue is the
   * insert state in node 0, it only makes sense to start modelling
   * at consensus column 1. */
  if(spos == 0 && spos_state == 1) spos = 1;
  if(epos == 0 && epos_state == 1) epos = 1;
  /* If most-likely HMM node to emit final position comes BEFORE or
   * EQUALS the most-likely HMM node to emit first position, our HMM
   * alignment is crap, default to using the full CM. (note: If
   * EQUALS we could be right, but we can't build a CM from a single
   * consensus column (see notes in cm_modelmaker.c::cm_from_guide),
   * and I would argue we don't really care about getting single
   * residue alignments correct anyway. 
   */
  if(epos <= spos) { spos = 1; epos = orig_cm->cp9->M; } 
  
  /* step 2. build the sub_cm from the original CM. */
  if((status = build_sub_cm(orig_cm, errbuf, &sub_cm, 
			    spos, epos,                /* first and last col of structure kept in the sub_cm  */
			    &submap,                   /* this maps from the sub_cm to cm and vice versa      */
			    0)) != eslOK)              /* don't print debugging info */
    return status;

  /* step 3. configure the sub_cm */
  if((status = cm_ConfigureSub(sub_cm, errbuf, -1, orig_cm, submap)) != eslOK) return status; 

  *ret_sub_cm = sub_cm;
  *ret_submap = submap;

  return eslOK;
}

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
  CM_P7_OM_HOLDER om_holder;     /* reusable LOCAL p7 profile/OPROFILE for the
				  * --p7pinbridge SW scan, built once and reused
				  * across this block (brief 26_0430-090). Safe because
				  * this call owns the whole block. */

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
  cm_p7_om_holder_Init(&om_holder);
  for(j = 0; j < sq_block->count; j++) {
    sqp = sq_block->list + j;
    if((status = DispatchSqAlignment(cm, errbuf, sqp, sq_block->first_seqidx + j, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, &om_holder, &(dataA[j]))) != eslOK) { cm_p7_om_holder_Reset(&om_holder); goto ERROR; }
  }
  cm_p7_om_holder_Reset(&om_holder);
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
 *           om_holder  - reusable LOCAL p7 profile/OPROFILE holder for the
 *                        --p7pinbridge SW scan (brief 26_0430-090); built once per
 *                        worker/block and reused across sequences. Can be NULL
 *                        (then the pinbridge wrapper builds/frees its own per
 *                        call). MUST be per-thread (not shared across threads).
 *           ret_data   - RETURN: newly created CM_ALNDATA object
 *
 * Returns:  eslOK on success;
 *           eslEINCOMPAT on contract violation, errbuf is filled;
 *           eslEMEM if we run out of memory;
 *           <ret_data> is alloc'ed and filled.
 */
int
DispatchSqAlignment(CM_t *cm, char *errbuf, ESL_SQ *sq, int64_t idx, float mxsize, char mode, int pass_idx,
		    int cp9b_valid, ESL_STOPWATCH *w, ESL_STOPWATCH *w_tot, ESL_RANDOMNESS *r,
		    CM_P7_OM_HOLDER *om_holder, CM_ALNDATA **ret_data)
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

  /* alignment options */
  int do_nonbanded = (cm->align_opts & CM_ALIGN_NONBANDED) ? TRUE  : FALSE;
  int do_qdb       = (cm->align_opts & CM_ALIGN_QDB)       ? TRUE  : FALSE;
  int do_hbanded   = (do_nonbanded || do_qdb)              ? FALSE : TRUE;
  int do_optacc    = (cm->align_opts & CM_ALIGN_OPTACC)    ? TRUE  : FALSE;
  int do_sample    = (cm->align_opts & CM_ALIGN_SAMPLE)    ? TRUE  : FALSE;
  int do_post      = (cm->align_opts & CM_ALIGN_POST)      ? TRUE  : FALSE;
  int do_sub       = (cm->align_opts & CM_ALIGN_SUB)       ? TRUE  : FALSE;
  int do_small     = (cm->align_opts & CM_ALIGN_SMALL)     ? TRUE  : FALSE;
  int do_trunc     = (cm->align_opts & CM_ALIGN_TRUNC)     ? TRUE  : FALSE;
  int do_xtau      = (cm->align_opts & CM_ALIGN_XTAU)      ? TRUE  : FALSE;
  int do_p7band    = (cm->align_opts & CM_ALIGN_P7BANDED)  ? TRUE  : FALSE;
  /* brief 26_0430-269: --mxsize auto-escalation. When enabled (default, CM_ALIGN_MXESC
   * set, cleared by --no-mxesc) AND the user did NOT force an engine, pick the cheapest
   * engine whose estimated CM-DP peak fits --mxsize: tier (a) standard free-OptAcc,
   * else (b) checkpointed sqrt(M) pinned-OptAcc, else (c) the CYK floor (D&C).  Applies
   * only to the HB free-OptAcc path (--ckpt/--small/--nonbanded/--sample/--sub each keep
   * their own engine).  See the tier-selection block just before the HB align dispatch. */
  int do_mxesc     = ((cm->align_opts & CM_ALIGN_MXESC)    &&
                      (! (cm->align_opts & CM_ALIGN_CHECKPT)) &&
                      do_optacc && (! do_sample) && (! do_small) &&
                      (! do_nonbanded) && (! do_qdb) && (! do_sub)) ? TRUE : FALSE;
  int mxesc_tier   = 0;    /* 0=none/not-decided, 'a'/'b'/'c' once decided (per seq) */
  int eff_checkpt  = (cm->align_opts & CM_ALIGN_CHECKPT) ? TRUE : FALSE; /* effective ckpt engine choice; do_mxesc may raise it per seq */
  int p7b_iterate_ran = FALSE; /* TRUE once cp9_IterateSeq2BandsP7B() ran for this seq (bands valid even on eslERANGE) */
  int doing_search = FALSE;
  /* Brief 26_0430-120: IBV HMM-divergence fallback. Set when cm_TrAlignHB / cm_AlignHB
   * fails on IBV-derived bands and we've already rebuilt with vitband for
   * this sequence; prevents infinite retry.
   */
  int ibv_fallback_used = FALSE;

#if eslDEBUGLEVEL >= 1
  printf("#DEBUG: in DispatchSqAlignment() %s\n", sq->name);
  printf("#DEBUG: \tdo_nonbanded: %d\n", do_nonbanded);
  printf("#DEBUG: \tdo_optacc:    %d\n", do_optacc);
  printf("#DEBUG: \tdo_sample:    %d\n", do_sample);
  printf("#DEBUG: \tdo_post:      %d\n", do_post);
  printf("#DEBUG: \tdo_sub:       %d\n", do_sub);
  printf("#DEBUG: \tdo_small:     %d\n", do_small);
  printf("#DEBUG: \tdo_trunc:     %d\n", do_trunc);
  printf("#DEBUG: \tdo_qdb:       %d\n", do_qdb);
  printf("#DEBUG: \tdoing_search: %d\n", doing_search);
#endif
  
  /* sub-mode specific variables (wouldn't be needed if sub mode were not supported) */
  CM_t        *orig_cm = cm;      /* pointer to the original CM */
  CM_t        *sub_cm  = NULL;    /* the sub CM */
  CMSubMap_t  *submap  = NULL;    /* map from mother CM to sub CM, and vice versa */
  Parsetree_t *full_tr = NULL;    /* converted parsetree to full CM */

  /* contract check */
  if(do_small  && do_hbanded)       ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do small and HMM banded alignment");
  if(do_small  && do_optacc)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do small and opt acc alignment");
  if(do_post   && do_small)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do PP and small alignment");
  if(do_optacc && do_sample)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to sample and do optacc alignment");
  if(do_sub    && do_small)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do sub and small alignment");
  if(do_sub    && do_trunc)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do sub and truncated alignment");
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

  if(w_tot != NULL) esl_stopwatch_Start(w_tot);

  /* do sub-mode specific pre-alignment steps, if nec */
  if(do_sub) { 
    if((status = sub_alignment_prep(cm, errbuf, sq, &submap, &sub_cm)) != eslOK) goto ERROR;
    cm = sub_cm;
  }

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
      /* brief 26_0628-059: declared here (this scope encloses BOTH the
       * if(!cp9b_valid){...do_p7band derivation...} block below AND the later
       * CP9-iterate/CYK-prepass/CM-alignment-DP code that follows it -- the
       * latter runs even when cp9b_valid is TRUE, i.e. bands already computed
       * by a prior call, so declaring inside if(!cp9b_valid) put these out of
       * scope for the alignment-DP timing). _p7b_kind stays NULL unless the
       * do_p7band branch actually runs, so the final #STAGETIME print (gated
       * on _p7b_kind != NULL) only fires for the do_p7band derivers this brief
       * targets (--p7ibv/--p7kmerchain). */
      const char *_p7b_kind = NULL;
      int    _st059_on = (getenv("BRIEF059_STAGETIME") != NULL);
      double _st059_a_s = 0., _st059_b_s = 0., _st059_bd_s = 0., _st059_c_s = 0., _st059_d_s = 0., _st059_ab_total_s = 0.;
      int    _st059_ab_split = FALSE;
      if(! cp9b_valid) {
	/* TODO #9 mitigation: --p7band produces too-narrow k-envelopes for small-M
	 * models (M < ~200), causing accuracy regression on rmark4e (Lacto-usp, atoC,
	 * snoZ152, ar45, SNORA47). Fall back to unbanded CP9 F/B for small CMs; the
	 * absolute wall savings from --p7band on tiny CMs is negligible. See brief 26_0430-059. */
#define P7BAND_MIN_M 0   /* was 200; gate was introduced in brief 26_0430-062 (session 19) for
                            the pinbridge-era band-derivation path. The current F+B IBV
                            (briefs 26_0430-116-126) is a different algorithm; gate removed per
                            brief 26_0430-141. */
	if(do_p7band && cm->fp7 != NULL && cm->fp7->M < P7BAND_MIN_M) {
	  fprintf(stderr, "#P7BAND_SKIP M=%d reason=small_M_acc_gap (threshold=%d)\n",
		  cm->fp7->M, P7BAND_MIN_M);
	  do_p7band = FALSE;
	}
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
	   * unbanded p7_Seq2BandsVit path actually needs it (brief 26_0430-094). The
	   * --p7pinbridge success path never touches it, so we avoid the eager
	   * full-matrix alloc that (a) defeats pinbridge's large-M memory win and
	   * (b) overflows int32 in p7_gmx_Create at M=L ~ 1.5e5 (e.g. HSV). */
	  tr_p7b = p7_trace_Create();

	  /* Build local nodepad copy with p7bpad (p7padplus) added */
	  int *local_nodepad = NULL;
	  if(cm->flags & CMH_P7NODEPAD) {
	    int k;
	    ESL_ALLOC(local_nodepad, sizeof(int) * (cm->fp7->M + 1));
	    for(k = 0; k <= cm->fp7->M; k++) local_nodepad[k] = cm->p7_cm_nodepad[k] + cm->p7bpad;
	  }

	  /* Derive p7 bands: either SW-pinbridge prefilter + banded Viterbi (--p7pinbridge)
	   * or full unbanded p7_GViterbi (default). Pinbridge wrapper falls back to
	   * full p7_Seq2BandsVit if the banded trace fails inside the prefilter band
	   * (signaled by ncells=0). Pinbridge is correct for both truncated and
	   * non-truncated alignment.
	   *
	   * M-threshold gate: pinbridge has per-sequence overhead (build OPROFILE,
	   * SSE scan setup, LSIS allocation, banded traceback) that on small models
	   * exceeds the cost of full p7_GViterbi at O(LM). For M < 200 the full
	   * Viterbi wins; gate pinbridge on M >= 200 to capture the big-M speedup
	   * without the tiny-M tail regressions. */
	  struct timespec _ta_p7b, _tb_p7b;
	  /* _p7b_kind and the _st059_* stage-timing vars are declared up at the
	   * enclosing if(! cp9b_valid) scope (brief 26_0628-059) so they stay in
	   * scope through the shared CP9-iterate/CYK-prepass/CM-align-DP code below. */
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
	    if (cm->p7_ibv_wv) {
	      /* Brief 26_0430-169: windowed-Viterbi band = MAP-trace i2k +/- F+B-halfwidth
	       * pad (cm->p7_wv_nodepad, calibrated once at align-time setup). */
	      _p7b_kind = "p7ibv-wv";
	      int *wv_nodepad = NULL;
	      int  wk;
	      if (cm->p7_wv_nodepad == NULL)
		ESL_FAIL(eslEINVAL, errbuf, "--p7ibv-wv: cm->p7_wv_nodepad not calibrated");
	      ESL_ALLOC(wv_nodepad, sizeof(int) * (cm->fp7->M + 1));
	      for (wk = 0; wk <= cm->fp7->M; wk++) wv_nodepad[wk] = cm->p7_wv_nodepad[wk] + cm->p7bpad;
	      status = p7_Seq2BandsWV(cm, errbuf, sq->dsq, sq->L, wv_nodepad,
				      do_trunc, /* brief 26_0430-171 */
				      &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	      free(wv_nodepad);
	    } else if (cm->p7_ibv_mem) {
	      _p7b_kind = "p7ibv-dnc";
	      status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, sq->L,
					   cm->p7_ibv_delta, cm->p7_ibv_base_slab,
					   TRUE, /* do_boundary_widen: CM-side preserves current behavior */
					   FALSE, /* brief 26_0430-172: do_kband (unbanded; --p7ibv-mem keeps exact delta band) */
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
	     * Empty rows default to [1, M] inside the kernel.
	     */
	  } else if (cm->p7_use_kmerchain) {
	    /* Brief 26_0628-027: genome-scale k-mer seed-and-chain. Collect all seeds
	     * genome-wide, chain by global colinearity, emit multi-segment pins
	     * into the unmodified p7_pins2bands_nodepad. Opt-in. */
	    _p7b_kind = "kmerchain";
	    status = p7_Seq2BandsKmerChain(cm, errbuf, sq->dsq, sq->L, local_nodepad,
	                                   do_trunc, /* brief 26_0628-033 */
	                                   &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells,
	                                   _st059_on ? &_st059_a_s : NULL, _st059_on ? &_st059_b_s : NULL,
	                                   _st059_on ? &_st059_bd_s : NULL);
	    if (_st059_on) _st059_ab_split = TRUE; /* provisional; cleared below if a fallback fires */
	    /* ncells==0 => M-gate/N-gate fired, or no usable chain. Brief 26_0628-047:
	     * default fallback is now --p7ibv's D&C deriver (this file's own
	     * --p7ibv-mem branch above, same defaults/params), instead of a
	     * Vit-trace band; --p7kmerchain-fbvit reverts to
	     * the old p7_Seq2BandsVit fallback. */
	    if (status == eslOK && p7_ncells == 0) {
	      _st059_ab_split = FALSE; /* brief 26_0628-059: a_s/b_s only cover the failed kmerchain attempt, not the fallback -- report combined ab_s instead */
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
					     TRUE, FALSE, do_trunc, cm->p7_ibv_mode, cm->p7_ibv_width,
					     &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
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
		/* one-line firing diagnostic so tests can count native-fallback
		 * invocations directly (the authoritative counter; the cosmetic
		 * "#P7BAND ... FAILED" line below just reflects p7_ncells==0). */
		fprintf(stderr, "#NATIVE_FALLBACK seq=%s L=%d M=%d\n",
			sq->name, (int)sq->L, cm->fp7 ? cm->fp7->M : 0);
		/* status stays eslOK, p7_ncells stays 0 -> falls through to the
		 * native cp9_Seq2Bands fallback below. Do nothing else here. */
	      }
	    }
	  } else if (cm->p7_use_pinbridge) {
	    _p7b_kind = "pinbridge";
	    status = p7_Seq2BandsPinBridgeWrap(cm, errbuf, gm_p7b, bg_p7b, tr_p7b,
					       sq->dsq, sq->L, cm->p7bpad,
					       local_nodepad,
					       0, 0, /* hopback=0, vitend=0 */
					       om_holder,
					       &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    /* If pinbridge couldn't produce a trace (rare; band missed the trace
	     * entirely), fall back to full unbanded Viterbi for this sequence. */
	    if (status == eslOK && p7_ncells == 0) {
	      if (gx_p7b == NULL) gx_p7b = p7_gmx_Create(cm->fp7->M, sq->L);
	      status = p7_Seq2BandsVit(errbuf, gm_p7b, gx_p7b, bg_p7b, tr_p7b,
				       sq->dsq, sq->L, cm->p7bpad,
				       local_nodepad,
				       0, 0,
				       &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
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
	    _st059_ab_total_s = _p7b_s; /* brief 26_0628-059: authoritative a+b total, robust across fallbacks */
	    fprintf(stderr, "#P7BAND_TIME %s kind=%s L=%d M=%d t=%.6f ncells=%d ibvmode=%d ibvwidth=%d ibvdelta=%d\n",
	            sq->name, _p7b_kind, (int)sq->L, cm->fp7->M, _p7b_s,
	            p7_ncells, cm->p7_ibv_mode, cm->p7_ibv_width, cm->p7_ibv_delta); /* brief 26_0430-140: ncells = band-size discriminator */
	  }

	  /* Debug: report Viterbi band stats */
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
	    fprintf(stderr, "#P7BAND %s L=%d M=%d npins=%d pin_range=[%d,%d] ncells=%d avg_bw=%.1f\n",
		    sq->name, (int)sq->L, cm->fp7->M, dbg_npin, dbg_minpin, dbg_maxpin, p7_ncells, dbg_avgbw);
	  }
	  else {
	    fprintf(stderr, "#P7BAND %s L=%d M=%d FAILED (status=%d ncells=%d)\n",
		    sq->name, (int)sq->L, cm->fp7->M, status, p7_ncells);
	  }

	  if(status == eslOK && p7_ncells > 0) {
	    /* Brief 26_0430-215 Phase A: OPTIONAL kmerchain-band TIGHTENING via a
	     * Viterbi MAP trace, wired into the REAL path (feeds bands_2 into
	     * cp9_IterateSeq2BandsP7B so it propagates to the final alignment).
	     * Env-gated: P215_TIGHTEN_N=<n> (constant half-width) or
	     * P215_TIGHTEN_PERNODE=1 (the CM's stored per-node pad). When NEITHER
	     * is set, t_kmin/t_kmax alias p7_kmin/p7_kmax => byte-identical to
	     * production. Phase A reuses the existing (unbounded, full-model)
	     * p7_Seq2BandsWV purely to obtain the Viterbi MAP trace i2k; each i2k[i]
	     * is CLAMPED into kmerchain's [p7_kmin,p7_kmax] (exactly mimicking a
	     * band-bounded Viterbi on the rows where they'd differ), so bands_1 =
	     * i2k +/- N always overlaps bands_0 and the row-wise min (bands_2) is
	     * never empty -- no disjoint-row fallback needed (see brief 214 addendum).
	     * Cost of this WV call is IRRELEVANT here (Phase B builds the real
	     * bounded kernel); Phase A is the accuracy gate only. */
	    int   *t_kmin = p7_kmin, *t_kmax = p7_kmax;   /* default: untightened alias */
	    int   *tight_kmin = NULL, *tight_kmax = NULL;
	    {
	      /* brief 26_0430-262: --p7vittighten/--p7vitcloud take precedence over the
	       * P215/P216/P248 getenv() family when used (cm->p215_mode != P215_MODE_OFF);
	       * else fall back to reading the env vars exactly as before (byte-identical
	       * default behavior -- gate G3). The CLI flags cover only the two
	       * genome-validated combos (constant-N pin, and cloud); PERNODE/PADPLUS/
	       * FLAT/VALIDATE stay env-only diagnostics regardless of flag usage. */
	      int         p215_cli_pin   = (cm->p215_mode == P215_MODE_PIN);
	      int         p215_cli_cloud = (cm->p215_mode == P215_MODE_CLOUD);
	      int         p215_cli_used  = (p215_cli_pin || p215_cli_cloud);
	      int         p215_pernode = 0, p215_N = -1;
	      if(p215_cli_used) {
		p215_N = p215_cli_pin ? cm->p215_tighten_n : 20; /* cloud: any N>=0 enters this scope; unused once p248_cloud fires below */
	      } else {
		const char *e_pn = getenv("P215_TIGHTEN_PERNODE");
		const char *e_n  = getenv("P215_TIGHTEN_N");
		if(e_pn != NULL && atoi(e_pn) != 0) p215_pernode = 1;
		if(e_n  != NULL)                     p215_N       = atoi(e_n);
	      }
	      int have_pernode = (cm->flags & CMH_P7NODEPAD) && cm->p7_cm_nodepad != NULL;
	      if(cm->p7_use_kmerchain && (p215_pernode || p215_N >= 0)) {
		if(p215_pernode && !have_pernode) {
		  fprintf(stderr, "#T215 seq=%s WARN pernode requested but CM lacks P7NODEPAD; NO tightening applied\n", sq->name);
		} else {
		  int  *wv_i2k = NULL, *wv_kmin = NULL, *wv_kmax = NULL, *zero_pad = NULL;
		  int  *i2k_c  = NULL, *nodepad215 = NULL, *b1_kmin = NULL, *b1_kmax = NULL;
		  int   wv_ncells = 0, b1_ncells = 0, k215, i215, status215, M215 = cm->fp7->M;
		  /* brief 26_0430-248: P248_CLOUD=1 replaces the pin+/-N band with the
		   * extband kernel's delta-CLOUD band (robust to the truncation-mode-flip
		   * collapse; see brief 247/248).  Cloud delta (milli-bits) via
		   * P248_CLOUD_DELTA, default = cm->p7_ibv_delta (--p7ibv-delta, 20 bits). */
		  int   p248_cloud = p215_cli_used ? p215_cli_cloud : (getenv("P248_CLOUD") != NULL);
		  int  *wv_cloud_kmin = NULL, *wv_cloud_kmax = NULL;
		  int   p248_delta = p215_cli_used
		    ? (p215_cli_cloud ? cm->p215_cloud_delta : cm->p7_ibv_delta)
		    : ((getenv("P248_CLOUD_DELTA") != NULL) ? atoi(getenv("P248_CLOUD_DELTA")) : cm->p7_ibv_delta);
		  struct timespec _twv0, _twv1;
		  ESL_ALLOC(zero_pad, sizeof(int) * (M215 + 1));
		  for(k215 = 0; k215 <= M215; k215++) zero_pad[k215] = 0;
		  /* (2) exact Viterbi MAP trace i2k (unbounded; Phase B replaces this
		   *     with a band-bounded kernel -- cost irrelevant to the accuracy gate). */
		  int p215_bounded = p215_cli_used ? 1 : (getenv("P215_BOUNDED") != NULL);  /* brief 215 Phase B: Viterbi bounded to kmerchain band; both CLI combos require it */
		  /* brief 26_0430-216 Phase 3: the compact O(L*bandwidth)-storage kernel is now
		   * the default under P215_BOUNDED (validated byte-identical i2k vs the flat
		   * oracle on tRNA/5S/RNaseP/SSU/LSU/norovirus/dengue, both do_trunc branches,
		   * 273 sequences, 0 mismatches -- see brief 216 summary).  P216_FLAT=1 is an
		   * escape hatch back to the O(L*M) flat oracle (diagnostic/rollback only; it
		   * OOMs or is impractically slow at genome scale, e.g. sarscov2 M~30kb).
		   * P216_VALIDATE=1 runs BOTH kernels and diffs i2k (bring-up / spot-check). */
		  int p216_flat     = (getenv("P216_FLAT")     != NULL);
		  int p216_validate = (getenv("P216_VALIDATE") != NULL);
		  clock_gettime(CLOCK_MONOTONIC, &_twv0);
		  if(p215_bounded && p216_flat)
		    status215 = p7_Seq2BandsIBV_extband(cm, errbuf, sq->dsq, sq->L, do_trunc,
							p7_kmin, p7_kmax, p248_delta, &wv_i2k,
							p248_cloud ? &wv_cloud_kmin : NULL,
							p248_cloud ? &wv_cloud_kmax : NULL);  /* O(L*M) flat oracle (rollback) */
		  else if(p215_bounded)
		    status215 = p7_Seq2BandsIBV_extband_compact(cm, errbuf, sq->dsq, sq->L, do_trunc,
							p7_kmin, p7_kmax, p248_delta, &wv_i2k,
							p248_cloud ? &wv_cloud_kmin : NULL,
							p248_cloud ? &wv_cloud_kmax : NULL);  /* O(L*bandwidth), compact storage */
		  else
		    status215 = p7_Seq2BandsWV(cm, errbuf, sq->dsq, sq->L, zero_pad, do_trunc,
					       &wv_i2k, &wv_kmin, &wv_kmax, &wv_ncells);  /* unbounded O(L*M) */
		  clock_gettime(CLOCK_MONOTONIC, &_twv1);
		  fprintf(stderr, "#T215_WVTIME seq=%s M=%d L=%d bounded=%d wv_s=%.6f\n", sq->name, M215, (int)sq->L, p215_bounded,
			  (_twv1.tv_sec - _twv0.tv_sec) + (_twv1.tv_nsec - _twv0.tv_nsec)/1e9);
		  if(status215 == eslOK && p215_bounded && p216_validate) {
		    int  *oracle_i2k = NULL;
		    int   ostat = p7_Seq2BandsIBV_extband(cm, errbuf, sq->dsq, sq->L, do_trunc, p7_kmin, p7_kmax, p248_delta, &oracle_i2k, NULL, NULL);
		    if(ostat == eslOK) {
		      int ndiff = 0, i216;
		      for(i216 = 0; i216 <= sq->L; i216++) if(oracle_i2k[i216] != wv_i2k[i216]) ndiff++;
		      fprintf(stderr, "#P216_VALIDATE seq=%s L=%d ndiff=%d %s\n", sq->name, (int)sq->L, ndiff,
			      ndiff == 0 ? "MATCH" : "MISMATCH");
		      free(oracle_i2k);
		    } else {
		      fprintf(stderr, "#P216_VALIDATE seq=%s ORACLE_FAILED status=%d\n", sq->name, ostat);
		    }
		  }
		  if(status215 == eslOK && p248_cloud && wv_cloud_kmin != NULL && wv_cloud_kmax != NULL) {
		    /* brief 26_0430-248: bands_2 = the extband kernel's delta-CLOUD band
		     * directly (already bounded to bands_0 by the kernel, so ⊆ bands_0 --
		     * a genuine tightening).  Replaces the fragile pin+/-N band that
		     * collapses under the M>L truncation-mode flip (brief 247). */
		    long km_totw = 0, cloud_totw = 0;
		    ESL_ALLOC(tight_kmin, sizeof(int) * (sq->L + 1));
		    ESL_ALLOC(tight_kmax, sizeof(int) * (sq->L + 1));
		    for(i215 = 0; i215 <= sq->L; i215++) { tight_kmin[i215] = wv_cloud_kmin[i215]; tight_kmax[i215] = wv_cloud_kmax[i215]; }
		    for(i215 = 1; i215 <= sq->L; i215++) { km_totw += (p7_kmax[i215]-p7_kmin[i215]+1); cloud_totw += (tight_kmax[i215]-tight_kmin[i215]+1); }
		    t_kmin = tight_kmin; t_kmax = tight_kmax;
		    fprintf(stderr, "#T248_CLOUD seq=%s M=%d L=%d delta=%d km_totw=%ld cloud_totw=%ld ratio=%.4f\n",
			    sq->name, M215, (int)sq->L, p248_delta, km_totw, cloud_totw,
			    km_totw > 0 ? (double)cloud_totw/(double)km_totw : 1.0);
		  }
		  if(status215 == eslOK && !p248_cloud) {
		    /* (2b) CLAMP each pinned i2k[i] into kmerchain's [kmin,kmax] (mimics a
		     *      band-bounded Viterbi -- the bounded MAP pin lies in bands_0). */
		    ESL_ALLOC(i2k_c, sizeof(int) * (sq->L + 1));
		    for(i215 = 0; i215 <= sq->L; i215++) {
		      int kk = wv_i2k[i215];
		      if(i215 == 0 || kk == -1) { i2k_c[i215] = kk; continue; }
		      if(kk < p7_kmin[i215]) kk = p7_kmin[i215];
		      else if(kk > p7_kmax[i215]) kk = p7_kmax[i215];
		      i2k_c[i215] = kk;
		    }
		    /* (3) bands_1 = i2k +/- N via the PRODUCTION pin->band converter (yields a
		     *     monotone, connected, DP-valid band -- do NOT hand-roll this). */
		    /* brief 26_0430-219: P215_TIGHTEN_PADPLUS=<n> overrides the +cm->p7bpad
		     * (--p7padplus) offset added to the calibrated per-node pad; default is
		     * cm->p7bpad (== what brief 218 tested). Set =0 for "raw pernode" (bare
		     * calibrated pad, no padplus) to disentangle the per-node SHAPE from the
		     * uniform +p7bpad offset -- the brief 218 pernode-vs-constant-N confound. */
		    {
		      const char *e_pp215 = getenv("P215_TIGHTEN_PADPLUS");
		      int p215_padplus = (e_pp215 != NULL) ? atoi(e_pp215) : cm->p7bpad;
		      ESL_ALLOC(nodepad215, sizeof(int) * (M215 + 1));
		      for(k215 = 0; k215 <= M215; k215++)
			nodepad215[k215] = p215_pernode ? (cm->p7_cm_nodepad[k215] + p215_padplus) : p215_N;
		    }
		    status215 = p7_pins2bands_nodepad(i2k_c, errbuf, sq->L, M215, nodepad215,
						      0, cm->p7_kmerchain_ramp_alpha,
						      &b1_kmin, &b1_kmax, &b1_ncells);
		  }
		  if(status215 == eslOK && b1_kmin != NULL) {
		    /* (4) bands_2 = per-row min(bands_0, bands_1). max()/min() of two monotone
		     *     bands stays monotone; a disjoint row (rare) falls back to kmerchain. */
		    long km_totw = 0, tight_totw = 0; int n_empty = 0;
		    int p215_pure = p215_cli_used ? p215_cli_pin : (getenv("P215_PURE") != NULL);  /* brief 215: pure Viterbi band (no kmerchain intersection); --p7vittighten implies it, --p7vitcloud doesn't reach this branch */
		    ESL_ALLOC(tight_kmin, sizeof(int) * (sq->L + 1));
		    ESL_ALLOC(tight_kmax, sizeof(int) * (sq->L + 1));
		    tight_kmin[0] = p215_pure ? b1_kmin[0] : p7_kmin[0];
		    tight_kmax[0] = p215_pure ? b1_kmax[0] : p7_kmax[0];
		    for(i215 = 1; i215 <= sq->L; i215++) {
		      int a = p7_kmin[i215], b = p7_kmax[i215];
		      int lo, hi;
		      if(p215_pure) {                               /* Eric's proposal: pure Viterbi+/-N band -> CP9 F/B */
			lo = b1_kmin[i215]; hi = b1_kmax[i215];
		      } else {                                      /* bands_2 = min(bands_0, bands_1) */
			lo = ESL_MAX(a, b1_kmin[i215]);
			hi = ESL_MIN(b, b1_kmax[i215]);
			if(hi < lo) { lo = a; hi = b; n_empty++; }  /* disjoint -> kmerchain fallback */
		      }
		      tight_kmin[i215] = lo; tight_kmax[i215] = hi;
		      km_totw += (b - a + 1); tight_totw += (hi - lo + 1);
		    }
		    t_kmin = tight_kmin; t_kmax = tight_kmax;
		    fprintf(stderr, "#T215 seq=%s M=%d L=%d mode=%s N=%d km_totw=%ld tight_totw=%ld ratio=%.4f n_empty=%d\n",
			    sq->name, M215, (int)sq->L, p215_pernode ? "pernode" : "const", p215_N,
			    km_totw, tight_totw, km_totw > 0 ? (double)tight_totw/(double)km_totw : 1.0, n_empty);
		  } else {
		    fprintf(stderr, "#T215 seq=%s WV/pins2bands_FAILED status=%d; NO tightening applied\n", sq->name, status215);
		  }
		  if(zero_pad)    free(zero_pad);
		  if(wv_i2k)      free(wv_i2k);
		  if(wv_kmin)     free(wv_kmin);
		  if(wv_kmax)     free(wv_kmax);
		  if(i2k_c)       free(i2k_c);
		  if(nodepad215)  free(nodepad215);
		  if(b1_kmin)     free(b1_kmin);
		  if(b1_kmax)     free(b1_kmax);
		  if(wv_cloud_kmin) free(wv_cloud_kmin);
		  if(wv_cloud_kmax) free(wv_cloud_kmax);
		}
	      }
	    }
	    /* Use p7 bands to derive CM bands via p7-banded CP9 F/B with tau-ratcheting */
	    struct timespec _ta_cp9, _tb_cp9;
	    clock_gettime(CLOCK_MONOTONIC, &_ta_cp9);
	    status = cp9_IterateSeq2BandsP7B(cm, errbuf, sq->dsq, sq->L, t_kmin, t_kmax,
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
	    clock_gettime(CLOCK_MONOTONIC, &_tb_cp9);
	    double _cp9_s = (_tb_cp9.tv_sec - _ta_cp9.tv_sec) + (_tb_cp9.tv_nsec - _ta_cp9.tv_nsec)/1e9;
	    _st059_c_s = _cp9_s; /* brief 26_0628-059: stage (c) band construction = HMM-band -> CM v/j-band conversion */
	    fprintf(stderr, "#P7PB_POST M=%d L=%d cp9_iterate=%.4f\n", cm->fp7->M, (int)sq->L, _cp9_s);
	    if(tight_kmin) free(tight_kmin);
	    if(tight_kmax) free(tight_kmax);
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
	     * non-banded CP9 F/B would only OOM (it did: brief-157 SIGSEGV via int32
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
	    if(do_mxesc && p7b_iterate_ran && status == eslERANGE &&
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
						cm->maxtau, NULL)) != eslOK) goto ERROR;
	      }
	      else {
	        if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq,
					 1, sq->L, cm->cp9b, doing_search, pass_idx, 0)) != eslOK) goto ERROR;
	      }
	    }
	  }
	}
	else if(do_xtau) { /* multiply tau (if nec) until required mx is below Mb limit (mxsize) */
	  if((status = cp9_IterateSeq2Bands(cm, errbuf, sq->dsq, 1, sq->L, pass_idx, mxsize, doing_search, do_sample, do_post, 1 /*do_iterate*/,
					    cm->maxtau, NULL)) != eslOK) goto ERROR;
	}
	else {
	  if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq, 
				     1, sq->L, cm->cp9b, doing_search, pass_idx, 0)) != eslOK) goto ERROR;
	}
	if(w != NULL) esl_stopwatch_Stop(w);
	secs_bands = (w == NULL) ? 0. : w->elapsed;
	tau     = cm->tau;
	thresh1 = cm->cp9b->thresh1;
	thresh2 = cm->cp9b->thresh2;
	/* note: we don't set these three if cp9b_valid is TRUE */
      }

      /* DIAGNOSTIC: always print CP9 banded matrix cell count (whether
       * cykbands is on or off) so we can compare CP9-only vs cykbands runs. */
      {
        double _cp9_cells = 0.;
        CP9Bands_t *_cp9b = cm->cp9b;
        int _v, _jp;
        for(_v = 0; _v < cm->M; _v++)
          for(_jp = 0; _jp <= _cp9b->jmax[_v] - _cp9b->jmin[_v]; _jp++)
            if(hd_min(_cp9b, _v, _jp) <= hd_max(_cp9b, _v, _jp))
              _cp9_cells += hd_max(_cp9b, _v, _jp) - hd_min(_cp9b, _v, _jp) + 1;
        fprintf(stderr, "#P7PB_POST M=%d L=%d cp9_band_cells=%.0f\n",
                (cm->fp7 ? cm->fp7->M : 0), (int)sq->L, _cp9_cells);
      }

      /* CYK pre-pass: run CYK on CP9 bands, derive tighter per-state bands before Inside/Outside */
      if(cm->p7_use_cykbands) {
	struct timespec _ta_cyk, _tb_cyk;
	clock_gettime(CLOCK_MONOTONIC, &_ta_cyk);

	/* measure original band area (total d-band cells across all states) */
	double _orig_cells = 0.;
	{
	  CP9Bands_t *_cp9b = cm->cp9b;
	  int _v, _jp;
	  for(_v = 0; _v < cm->M; _v++)
	    for(_jp = 0; _jp <= _cp9b->jmax[_v] - _cp9b->jmin[_v]; _jp++)
	      if(hd_min(_cp9b, _v, _jp) <= hd_max(_cp9b, _v, _jp))
		_orig_cells += hd_max(_cp9b, _v, _jp) - hd_min(_cp9b, _v, _jp) + 1;
	}

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
	 * same estimator + mxsize semantics briefs 225/226/269 already use for the main
	 * alignment engine); if it doesn't fit, fall back to HMM-banded D&C-CYK
	 * (CYKDivideAndConquerHB/TrCYKDivideAndConquerHB), gated on the D&C estimate
	 * (cm_[Tr]DnCAlignSizeNeededHB, brief 228's banded-EL version, already on this base).
	 * If even the D&C estimate doesn't fit, this degrades exactly like the legacy path:
	 * _cyk_ok stays FALSE, original bands kept, no crash. Budget = full mxsize for both
	 * tiers: the pre-pass matrix is transient (allocated, consumed for the parsetree,
	 * freed before the main Inside/Outside/OA/Posterior alignment), so it is never
	 * concurrent with the main alignment's own matrices -- sequential, not additive. */
	char  _cykpp_engine  = 'f'; /* 'f' = full-matrix (tier a), 'd' = D&C (tier b), 'x' = neither fits */
	float _cykpp_est_std = 0., _cykpp_est_dnc = 0.;
	int   _cykpp_use_dnc = FALSE;

	if(cm->p7_cykbands_no_dnc) {
	  _cykpp_engine = 'f'; /* legacy: always attempt the full-matrix engine below */
	}
	else if(getenv("CYKPP_FORCE_DNC") != NULL) {
	  /* brief 26_0430-273 gate (a): force the D&C engine regardless of the estimate, so
	   * the choice can be isolated from --mxsize (which also drives the SEPARATE mxesc
	   * tier selector for the main alignment, brief 26_0430-269, a few hundred lines
	   * above this block in the same function) -- debug/gate-only, not for production use. */
	  _cykpp_engine = 'd'; _cykpp_use_dnc = TRUE;
	}
	else {
	  int _cykpp_status;
	  if(do_trunc) _cykpp_status = cm_TrAlignSizeNeededHB(cm, errbuf, sq->L, mxsize, FALSE/*do_sample*/, FALSE/*do_post*/,
	                                                      NULL, NULL, NULL, NULL, &_cykpp_est_std, NULL);
	  else         _cykpp_status = cm_AlignSizeNeededHB  (cm, errbuf, sq->L, mxsize, FALSE/*do_sample*/, FALSE/*do_post*/,
	                                                      NULL, NULL, NULL, NULL, &_cykpp_est_std, NULL);
	  if(_cykpp_status == eslOK) {
	    _cykpp_engine = 'f';
	  }
	  else {
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
	    if(_cykpp_est_dnc <= mxsize) { _cykpp_engine = 'd'; _cykpp_use_dnc = TRUE; }
	    else                         { _cykpp_engine = 'x'; }
	  }
	}

	fprintf(stderr, "#CYKPP_ENGINE seq=%s L=%d M=%d trunc=%d engine=%c est_std=%.1f est_dnc=%.1f mxsize=%.1f\n",
	        sq->name, (int)sq->L, (cm->fp7 ? cm->fp7->M : cm->clen), do_trunc, _cykpp_engine,
	        _cykpp_est_std, _cykpp_est_dnc, (float)mxsize);

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

	double _tight_cells = _orig_cells;
	if(_cyk_ok) {
	  /* brief 26_0430-243: tighten via the SHARED bounded helper
	   * (cykbands_tighten_bounded: baseline snapshot + never-loosen intersection
	   * + bounded hd-recompute + backstop) -- the EXACT mechanism CKPT_CYKBANDS
	   * uses (brief 26_0430-238).  This replaces the old raw cm_BandsFromCYKParsetree()
	   * call, whose unvisited-state full-envelope fallback could inflate the pass-2
	   * band area WITHOUT BOUND and outright fail at genome scale (brief 26_0430-242:
	   * band_area_ratio up to 145924x -> 343 GB matrix -> 4 dengue cells ESL_XFAIL'd).
	   * The intersection bounds every state's i/j band to the pre-tighten CP9/kmerchain
	   * baseline that cm->cp9b holds right here, so the result is bounded-by-construction.
	   *   pad      = cm->p7_cykbands_pad (from --cykpad); NOT CKPT_CYKBANDS_PAD.
	   *   maxratio = 3.0 plain-constant backstop; NOT CKPT_CYKBANDS_MAXRATIO -- the
	   *              --cykbands path must not read the CKPT_* env vars (constraint 1).
	   *              With never-loosen this backstop should never fire.
	   *   preserve_valid = do_trunc: cm_BandsFromCYKParsetree clobbers marginal
	   *              validity to J-only in truncated mode, which pass-2's cm_TrAlignHB
	   *              needs intact (constraint 3; this also fixes a latent validity bug
	   *              in the old raw-call trunc path).
	   * Per-state pad was archived 2026-05-19 (see cm_CYKPerstatePadCompute doc
	   * comment for failure analysis). Production uses uniform pad. */
	  double  _oc = 0., _tc = 0., _rr = 0.;
	  int64_t _ns = 0, _nre = 0;
	  if(cykbands_tighten_bounded(cm, errbuf, _cyk_tr, (int) sq->L, pass_idx,
				      do_trunc/*preserve_valid*/, cm->p7_cykbands_pad, 3.0/*maxratio backstop*/,
				      &_oc, &_tc, &_ns, &_nre, &_rr) == eslOK) {
	    _tight_cells = _tc;
	    if((_ns > 0 || _nre > 0 || _rr > 0.) && getenv("CYKBANDS_VERBOSE"))
	      fprintf(stderr, "#CYKBANDS never-loosen: %" PRId64 " state(s) CP9-sentinel + %" PRId64
		      " state(s) real-empty%s (M=%d L=%d)\n",
		      _ns, _nre, (_rr > 0. ? " + blowup-guard REVERTED" : ""), cm->M, (int) sq->L);
	  }
	  FreeParsetree(_cyk_tr);
	}

	clock_gettime(CLOCK_MONOTONIC, &_tb_cyk);
	double _cyk_s   = (_tb_cyk.tv_sec  - _ta_cyk.tv_sec)  + (_tb_cyk.tv_nsec  - _ta_cyk.tv_nsec) /1e9;
	double _b_ratio = (_orig_cells > 0.) ? _tight_cells / _orig_cells : 1.0;
	fprintf(stderr, "#P7PB_POST M=%d L=%d cyk_prepass=%.4f band_area_ratio=%.4f orig_cells=%.0f tight_cells=%.0f\n",
		(cm->fp7 ? cm->fp7->M : 0), (int)sq->L, _cyk_s, _b_ratio, _orig_cells, _tight_cells);
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
	  struct timespec _ta_cm, _tb_cm;
	  clock_gettime(CLOCK_MONOTONIC, &_ta_cm);
	CM_ALIGN_HB_RETRY:
	  /* brief 26_0430-269: --mxsize engine auto-escalation. Bands (cm->cp9b) are now
	   * derived (possibly maxtau-wide). Pick the cheapest engine whose estimated CM-DP
	   * peak fits --mxsize, per-sequence: (a) standard free-OptAcc, (b) checkpointed
	   * sqrt(M) pinned-OptAcc, (c) HMM-banded D&C-CYK floor. The comparison quantity is the CM DP
	   * matrices only (excludes the CP9 F/B), matching stock --mxsize semantics
	   * (cm_AlignSizeNeededHB:587 compares cmtotmb, not totmb). Estimators are
	   * safe-overestimates (briefs 225/226), so a "fits" verdict never under-provisions. */
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
	    if(do_trunc) status = cm_TrAlignSizeNeededHB(cm, errbuf, sq->L, mxsize, do_sample, do_post,
	                                                 NULL, NULL, NULL, NULL, &est_std_cm, &est_std_tot);
	    else         status = cm_AlignSizeNeededHB  (cm, errbuf, sq->L, mxsize, do_sample, do_post,
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
	         * post-brief-227 (banded EL deck), and brief 228's banded-EL estimator (cherry-
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
	    fprintf(stderr, "#MXESC seq=%s L=%d M=%d trunc=%d tier=%c est_std=%.1f est_ckpt=%.1f est_dnc=%.1f mxsize=%.1f ckpt_avail=%d\n",
	            sq->name, (int)sq->L, (cm->fp7 ? cm->fp7->M : cm->clen), do_trunc, mxesc_tier,
	            est_std_cm, ck_cmmb, dnc_tot, (float)mxsize, ckpt_avail);

	    if(mxesc_tier == 'c') {
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
	      if(getenv("INFERNAL_CKPT_VERBOSE"))
	        fprintf(stderr, "# mxesc tier (c) D&C-CYK floor engaged: M=%d L=%d trunc=%d (%s)\n",
	                (cm->fp7 ? cm->fp7->M : cm->clen), (int)sq->L, do_trunc,
	                (cm->flags & (CMH_LOCAL_BEGIN|CMH_LOCAL_END)) ? "local" : "global");
	      goto MXESC_ALN_DONE;
	    }
	  }
	  if(do_trunc) {
		/* brief 26_0430-126 merge: keep cd577024's #DBG-009 instrumentation, but route
		 * SizeNeededHB failure to CM_ALIGN_HB_CHECK_FB (IBV vitband fallback)
		 * instead of directly to ERROR, so the brief-120 IBV fallback stays live
		 * in the trunc path. For non-IBV runs CHECK_FB falls through to ERROR.
		 * brief 26_0430-269: under do_mxesc the tier selector above already validated
		 * the chosen engine's fit and set mb_tot, so skip this eslERANGE size-gate
		 * (for tier b it would spuriously fail on the full-matrix estimate). */
		if(! do_mxesc) {
		  status = cm_TrAlignSizeNeededHB(cm, errbuf, sq->L, mxsize, do_sample, do_post,
					    NULL, NULL, NULL, NULL, NULL, &mb_tot);
		  fprintf(stderr, "#DBG-009 trunc SizeNeededHB status=%d mb_tot=%.2f mxsize=%.2f do_post=%d errbuf=[%s]\n",
			status, mb_tot, (float) mxsize, do_post, errbuf);
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
		  /* brief 26_0610-094: rung-4 runtime engine selector.  Default (unset/0)
		   * = the sqrt(M) combined-mode (J/L/R/T) checkpointed CYK sweep,
		   * cm_CheckptTrCYKAlignHB() (brief 080, R3): the only sqrt(M) way to
		   * resolve the mode for a STRUCTURED CM (an unpinned sqrt(M) bifurcation
		   * Inside hits the brief-017 2D-coupling wall).  Set
		   * INFERNAL_CKPT_FORCE_DNC to force the alternative, now-fully-bug-fixed
		   * (briefs 081-093) per-candidate-mode TrCYKDivideAndConquerHB ncand-loop
		   * instead -- a smaller-footprint fallback/comparison engine, re-routed
		   * here unchanged rather than reimplemented.  Both paths resolve the same
		   * (mode, bkind/kpin/bbmode/blmode/brmode) pins and feed the identical
		   * pass-2 checkpointed pinned posterior + OptAcc + traceback below;
		   * only pass-1 mode/pin resolution differs. */
		  int   *bkind = NULL, *kpin = NULL;
		  char  *bbmode = NULL, *blmode = NULL, *brmode = NULL;
		  Parsetree_t *tr_best = NULL;
		  char   r4_mode = TRMODE_UNKNOWN;
		  float  r4_cyk  = IMPOSSIBLE, r4_Z = 0.;
		  int    use_dnc_r4 = (getenv("INFERNAL_CKPT_FORCE_DNC") != NULL) ? TRUE : FALSE;
		  ESL_ALLOC(bkind,  sizeof(int)  * cm->M);
		  ESL_ALLOC(kpin,   sizeof(int)  * cm->M);
		  ESL_ALLOC(bbmode, sizeof(char) * cm->M);
		  ESL_ALLOC(blmode, sizeof(char) * cm->M);
		  ESL_ALLOC(brmode, sizeof(char) * cm->M);
		  if(use_dnc_r4) {
		    /* D&C ncand-loop path (from trcyk-dnc-localbug-impl, briefs 081-093):
		     * run TrCYKDivideAndConquerHB once per root-valid marginal mode; the
		     * argmax penalty-folded root score picks the mode (each score ==
		     * the oracle cm_TrCYKInsideAlignHB()'s {J,L,R,T}alpha[0][L][L]),
		     * and that mode's parse supplies the (kind,k*,modes) pins. */
		    char cand[4]; int ncand = 0, m;
		    if(cm->cp9b->Jvalid[0]) cand[ncand++] = TRMODE_J;
		    if(cm->cp9b->Lvalid[0]) cand[ncand++] = TRMODE_L;
		    if(cm->cp9b->Rvalid[0]) cand[ncand++] = TRMODE_R;
		    if(cm->cp9b->Tvalid[0]) cand[ncand++] = TRMODE_T;
		    for(m = 0; m < ncand; m++) {
		      Parsetree_t *tr_m = NULL; char mm = TRMODE_UNKNOWN;
		      float sc_m = TrCYKDivideAndConquerHB(cm, sq->dsq, sq->L, 0, 1, sq->L, pass_idx, cand[m], &mm, &tr_m, cm->cp9b);
		      if(sc_m > r4_cyk) { r4_cyk = sc_m; r4_mode = cand[m]; if(tr_best) FreeParsetree(tr_best); tr_best = tr_m; }
		      else if(tr_m) FreeParsetree(tr_m);
		    }
		    if(tr_best == NULL) {
		      free(bkind); free(kpin); free(bbmode); free(blmode); free(brmode);
		      ESL_XFAIL(eslEINCOMPAT, errbuf, "rung-4 --ckpt (INFERNAL_CKPT_FORCE_DNC): no root-valid truncation mode for %s", sq->name);
		    }
		    status = eslOK;
		  }
		  else {
		    status = cm_CheckptTrCYKAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, pass_idx, &tr_best, &r4_mode, &r4_cyk);
		  }
		  if(status != eslOK) { free(bkind); free(kpin); free(bbmode); free(blmode); free(brmode); goto CM_ALIGN_HB_CHECK_FB; }
		  rung4_trpins_from_cyk(cm, tr_best, bkind, kpin, bbmode, blmode, brmode);
		  /* brief 26_0430-234: opt-in CKPT_CYKBANDS -- tighten pass-2 bands from
		   * this same (truncated) pass-1 CYK tree before freeing it.  preserve_valid
		   * TRUE: cm_BandsFromCYKParsetree clobbers marginal validity to J-only in
		   * do_trunc mode; save/restore keeps the resolved r4_mode valid (see helper
		   * doc). Only spatial i/j/hd bands tighten. */
		  /* brief 26_0430-271 item 2: --ckpt-cykbands promotes this env-only diagnostic
		   * to a real CLI flag; the env var is still honored as a sweep override. */
		  if((getenv("CKPT_CYKBANDS") != NULL) || (cm->align_opts & CM_ALIGN_CKPT_CYKBANDS)) {
		    double _oc = 0., _tc = 0.;
		    if(ckpt_cykbands_tighten(cm, errbuf, tr_best, (int) sq->L, pass_idx,
					     TRUE/*preserve_valid*/, &_oc, &_tc) == eslOK &&
		       (getenv("INFERNAL_CKPT_VERBOSE") || getenv("CKPT_CYKBANDS_VERBOSE")))
		      fprintf(stderr, "#CKPT_CYKBANDS rung=4 M=%d L=%d mode=%c band_area_ratio=%.4f orig_cells=%.0f tight_cells=%.0f\n",
			      cm->M, (int) sq->L,
			      (r4_mode==TRMODE_J)?'J':(r4_mode==TRMODE_L)?'L':(r4_mode==TRMODE_R)?'R':'T',
			      (_oc > 0. ? _tc/_oc : 1.0), _oc, _tc);
		  }
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
		  if(getenv("INFERNAL_CKPT_VERBOSE"))
		    fprintf(stderr, "# rung-4 checkpointed structured truncated OptAcc engaged: M=%d L=%d mode=%c engine=%s (%s, truncated, bps>0)\n",
			    cm->M, (int) sq->L, (r4_mode==TRMODE_J)?'J':(r4_mode==TRMODE_L)?'L':(r4_mode==TRMODE_R)?'R':'T',
			    use_dnc_r4 ? "dnc" : "ckpt",
			    (cm->flags & (CMH_LOCAL_BEGIN|CMH_LOCAL_END)) ? "local" : "global");
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
	if(! do_mxesc) {
	  if((status = cm_AlignSizeNeededHB(cm, errbuf, sq->L, mxsize, do_sample, do_post,
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
	 * CYK (brief 084, R1-L) supplies the bifurcation k* pins; pass 2 = checkpointed
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
	  /* pass 1: HMM-banded CYK parse -> k* pins.  brief 26_0610-098: rung-3
	   * runtime engine selector, UNIFIED with rung-4's (same INFERNAL_CKPT_FORCE_DNC
	   * env var -- one meaning, "use D&C", at every rung it applies to).
	   * Default (unset/0) = the checkpointed sqrt(M)-memory CYK engine
	   * (cm_CheckptCYKAlignHB, brief 26_0610-084 R1-L), which closed the pass-1
	   * memory floor (briefs 26_0610-040/041) so the whole rung-3 pipeline is
	   * checkpointed.  Set INFERNAL_CKPT_FORCE_DNC to force the D&C alternative
	   * CYKDivideAndConquerHB (brief 26_0610-041/098) -- a smaller-footprint
	   * fallback that carries the CP9 bands down a divide-and-conquer recursion.
	   * Both engines produce the SAME exact-CYK parse (single-optimum max-DP), so
	   * both feed the identical pass-2 checkpointed pinned posterior + OptAcc +
	   * traceback below via the same generic Parsetree_t -> rung3_kpin_from_cyk;
	   * only pass-1 k*-pin resolution differs. */
	  int    use_dnc_r3 = (getenv("INFERNAL_CKPT_FORCE_DNC") != NULL) ? TRUE : FALSE;
	  if(! has_bif) {
	    int kv; for (kv = 0; kv < cm->M; kv++) kpin[kv] = -1;
	    status = eslOK;
	  }
	  else {
	    if(use_dnc_r3) {
	      /* CYKDivideAndConquerHB() has no status/errbuf return -- it cm_Fail()s
	       * internally on its trust-but-verify root-state checks -- so a returned
	       * tr is always valid; treat as eslOK. */
	      (void) CYKDivideAndConquerHB(cm, sq->dsq, sq->L, 0, 1, sq->L, &tr_cyk, cm->cp9b);
	      status = eslOK;
	    }
	    else {
	      status = cm_CheckptCYKAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, &tr_cyk, NULL);
	    }
	    if(status == eslOK) {
	      rung3_kpin_from_cyk(cm, tr_cyk, kpin);
	      /* brief 26_0430-234: opt-in CKPT_CYKBANDS -- tighten pass-2 bands from
	       * this same (non-truncated) pass-1 CYK tree before freeing it. No
	       * valid-flag save/restore needed (rung-3 is non-truncated). */
	      /* brief 26_0430-271 item 2: --ckpt-cykbands promotes this env-only diagnostic
	       * to a real CLI flag; the env var is still honored as a sweep override. */
	      if((getenv("CKPT_CYKBANDS") != NULL) || (cm->align_opts & CM_ALIGN_CKPT_CYKBANDS)) {
		double _oc = 0., _tc = 0.;
		if(ckpt_cykbands_tighten(cm, errbuf, tr_cyk, (int) sq->L, pass_idx,
					 FALSE/*preserve_valid*/, &_oc, &_tc) == eslOK &&
		   (getenv("INFERNAL_CKPT_VERBOSE") || getenv("CKPT_CYKBANDS_VERBOSE")))
		  fprintf(stderr, "#CKPT_CYKBANDS rung=3 M=%d L=%d band_area_ratio=%.4f orig_cells=%.0f tight_cells=%.0f\n",
			  cm->M, (int) sq->L, (_oc > 0. ? _tc/_oc : 1.0), _oc, _tc);
	      }
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
	  if(getenv("INFERNAL_CKPT_VERBOSE"))
	    fprintf(stderr, "# rung-3 checkpointed structured OptAcc engaged: M=%d L=%d engine=%s (%s, non-truncated, bps>0)\n", cm->M, (int) sq->L,
		    (! has_bif) ? "skip-p1(bif0)" : (use_dnc_r3 ? "dnc" : "ckpt"),
		    (cm->flags & (CMH_LOCAL_BEGIN|CMH_LOCAL_END)) ? "local" : "global");
	}
	else {
	  status = cm_AlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, do_optacc, do_sample, cm->hb_mx, cm->hb_shmx,
			      cm->hb_omx, cm->hb_emx, r, do_post ? &ppstr : NULL, &tr, &pp, &sc);
	}
      }
    CM_ALIGN_HB_CHECK_FB:
      /* Brief 26_0430-120 IBV HMM-divergence fallback: cm_TrInsideAlignHB returns
       * eslEAMBIGUOUS "no valid parsetree found" on the 3/14 brief-117 seqs
       * where HMM-Viterbi disagrees with the CM's preferred parse (brief 26_0430-117
       * §5). Re-derive bands using the unbanded p7_Seq2BandsVit path (PAD=20
       * around the p7-Viterbi trace) and retry the CM alignment once.
       */
      if (status != eslOK && cm->p7_use_ibv && !ibv_fallback_used) {
	fprintf(stderr, "#P7IBV_FALLBACK %s L=%d M=%d status=%d errbuf='%s'\n",
		sq->name, (int)sq->L, (cm->fp7 ? cm->fp7->M : 0), status, errbuf);
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
	if (status != eslOK) goto ERROR;
	goto CM_ALIGN_HB_RETRY;
      }
      if (status != eslOK) goto ERROR;
      clock_gettime(CLOCK_MONOTONIC, &_tb_cm);
      double _cm_s = (_tb_cm.tv_sec - _ta_cm.tv_sec) + (_tb_cm.tv_nsec - _ta_cm.tv_nsec)/1e9;
      _st059_d_s = _cm_s; /* brief 26_0628-059: stage (d) CM alignment DP itself */
      fprintf(stderr, "#P7PB_POST M=%d L=%d cm_align_hb=%.4f%s\n",
              (cm->fp7 ? cm->fp7->M : 0), (int)sq->L, _cm_s,
              ibv_fallback_used ? " ibv_fallback=1" : "");
      /* brief 26_0628-059: single per-sequence 4-stage timing line, gated by
       * BRIEF059_STAGETIME (silent no-op by default). See the declaration
       * comment above (near _ta_p7b) for the a_s/b_s vs ab_s convention.
       * _p7b_kind != NULL guards against the do_xtau/plain-CP9 branches
       * (this shared cm-align-DP code also runs for those, but they never
       * went through a do_p7band deriver, so there's no method to report). */
      if (_st059_on && _p7b_kind != NULL) {
        if (_st059_ab_split)
          /* 2026-07-11 (brief 190 follow-up): bd_s = converting the winning
           * kmerchain chain's pins into HMM bands (p7_pins2bands_nodepad),
           * split out of what was previously folded into b_s. bd_s is 0 (and
           * meaningless) for non-kmerchain methods (p7ibv/p7ibv-wv/etc.),
           * which don't route through p7_Seq2BandsKmerChain at all. */
          fprintf(stderr, "#STAGETIME seq=%s L=%d M=%d method=%s a_s=%.6f b_s=%.6f bd_s=%.6f c_s=%.6f d_s=%.6f tot_s=%.6f\n",
                  sq->name, (int)sq->L, cm->fp7->M, _p7b_kind,
                  _st059_a_s, _st059_b_s, _st059_bd_s, _st059_c_s, _st059_d_s,
                  _st059_a_s + _st059_b_s + _st059_bd_s + _st059_c_s + _st059_d_s);
        else
          fprintf(stderr, "#STAGETIME seq=%s L=%d M=%d method=%s ab_s=%.6f c_s=%.6f d_s=%.6f tot_s=%.6f\n",
                  sq->name, (int)sq->L, cm->fp7->M, _p7b_kind,
                  _st059_ab_total_s, _st059_c_s, _st059_d_s,
                  _st059_ab_total_s + _st059_c_s + _st059_d_s);
      }
    }
  }

 MXESC_ALN_DONE: /* brief 26_0430-269: tier-(c) D&C-CYK floor lands here, skipping the free/ckpt OptAcc engine dispatch */

  if(do_sub) { /* add size of original CM's CP9 matrices used for calculating start/end position */
    mb_tot += orig_cm->cp9_mx->size_Mb;
    mb_tot += orig_cm->cp9_bmx->size_Mb;
  }

  if(w != NULL) esl_stopwatch_Stop(w);
  secs_aln = (w == NULL) ? 0. : w->elapsed;

  if(do_sub) { 
    /* convert sub cm parsetree to a full CM parsetree */
    if((status = sub_cm2cm_parsetree(orig_cm, cm, &full_tr, tr, submap, 0)) != eslOK) ESL_XFAIL(status, errbuf, "out of memory, converting sub parsetree to full parsetree");
    /* free sub data structures, we're done with them */
    FreeParsetree(tr);   tr     = full_tr;
    FreeCM(cm);          cm     = orig_cm;
    FreeSubMap(submap);  submap = NULL;
  }
  
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

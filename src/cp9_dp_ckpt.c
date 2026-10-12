/* cp9_dp_ckpt.c
 *
 * Checkpointed (memory-bounded), unbanded CP9 HMM Forward/Backward and
 * posterior HMM band calculation for CM alignment.
 *
 * This is a reference implementation of the band calculation done by
 * cp9_Forward(), cp9_Backward() and cp9_FB2HMMBands() (and the node
 * occupancies used by cp9_PredictStartAndEndPositions()), with the same
 * scaled-integer arithmetic, done in O(sqrt(L) * M) memory instead of the
 * three (L+1) x (M+1) matrices the full-matrix path holds. The HMM bands
 * it produces are identical to the full-matrix ones: every DP cell is
 * computed by the same expression, in the same order, with the same
 * ILogsum() calls, and every posterior mass is accumulated in the same
 * direction and order.
 *
 * It handles alignment only: Forward and Backward are run in 'align
 * mode' (the parse must start at i0 and end at j0; do_scan = FALSE in
 * the full-matrix functions). Search does not use it.
 *
 * Scheme: the band on each HMM state is set by two cumulative posterior
 * mass sweeps, one from the left (minimum position, ascending i) and one
 * from the right (maximum position, descending i). Integer ILogsum() is
 * not associative, so each sweep must see the posterior rows in its own
 * order. The DP rows 0..L are split into segments of R rows. Then:
 *
 *   1. Backward, rows L..0, keeping two rows; save the last row of each
 *      segment, and the total score sc = B[0].M_0.
 *   2. Left sweep: for each segment left to right, recompute its Backward
 *      rows from its saved row; continue Forward (two rows) through the
 *      segment, saving its first row; combine each Forward row with the
 *      Backward row for the same position into posteriors.
 *   3. Right sweep: for each segment right to left, recompute its Forward
 *      rows from its saved row; continue Backward (two rows) through the
 *      segment; combine.
 *
 * Five DP fills in all (Backward, Backward segments, Forward, Forward
 * segments, Backward) and two posterior passes. Memory is the saved rows
 * of both kinds plus one segment: (2 * nseg + R + 4) rows of 4 (M+1)
 * ints, smallest at R ~ sqrt(2(L+1)).
 *
 * The posterior masses do not depend on tau, only the thresholds they are
 * compared with do. So the sweeps compare each state's mass with several
 * thresholds at once, one for each tau the band tightening (tau
 * multiplication) might try, and one calculation serves every iteration
 * of the tightening loop.
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

#define CP9TSC(s,k) (tsc[(k) * cp9O_NTRANS + (s)])

/* one row of a CP9 DP matrix: match, insert, delete and EL scores [0..M] */
typedef struct {
  int *m;
  int *i;
  int *d;
  int *el;
} cp9ck_row_t;

static void
cp9ck_copy_row(cp9ck_row_t *dst, cp9ck_row_t *src, int M)
{
  memcpy(dst->m,  src->m,  sizeof(int) * (M+1));
  memcpy(dst->i,  src->i,  sizeof(int) * (M+1));
  memcpy(dst->d,  src->d,  sizeof(int) * (M+1));
  memcpy(dst->el, src->el, sizeof(int) * (M+1));
}

/* Forward row 0 (position i0-1). Same as the initialization of the zero
 * row in cp9_Forward(). */
static void
cp9ck_fwd_row0(CP9_t *cp9, cp9ck_row_t *cur)
{
  int const *tsc = cp9->otsc;
  int M = cp9->M;
  int k, sc;

  cur->m[0]  = 0;
  cur->i[0]  = -INFTY;
  cur->d[0]  = -INFTY;
  cur->el[0] = -INFTY;
  for (k = 1; k <= M; k++) {
    cur->m[k] = cur->i[k] = cur->el[k] = -INFTY;
    sc = ILogsum(ILogsum(cur->m[k-1] + CP9TSC(cp9O_MD,k-1),
                         cur->i[k-1] + CP9TSC(cp9O_ID,k-1)),
                 cur->d[k-1] + CP9TSC(cp9O_DD,k-1));
    cur->d[k] = sc;
  }
}

/* Forward row for residue j, from the row for j-1. Same as one iteration
 * of the main loop of cp9_Forward() with do_scan FALSE (the end state
 * score, which band calculation does not use, is not computed). */
static void
cp9ck_fwd_row(CP9_t *cp9, ESL_DSQ *dsq, int j, cp9ck_row_t *prv, cp9ck_row_t *cur)
{
  int const *tsc = cp9->otsc;
  int const *isc = cp9->isc[dsq[j]];
  int const *msc = cp9->msc[dsq[j]];
  int el_selfsc  = cp9->el_selfsc;
  int M = cp9->M;
  int k, c, sc;

  cur->m[0]  = -INFTY;
  cur->d[0]  = -INFTY;
  cur->el[0] = -INFTY;

  sc = ILogsum(ILogsum(prv->m[0] + CP9TSC(cp9O_MI,0),
                       prv->i[0] + CP9TSC(cp9O_II,0)),
               prv->d[0] + CP9TSC(cp9O_DI,0));
  cur->i[0] = sc + isc[0];

  for (k = 1; k <= M; k++) {
    /* match state */
    sc = ILogsum(ILogsum(prv->m[k-1] + CP9TSC(cp9O_MM,k-1),
                         prv->i[k-1] + CP9TSC(cp9O_IM,k-1)),
                 ILogsum(prv->d[k-1] + CP9TSC(cp9O_DM,k-1),
                         prv->m[0]   + CP9TSC(cp9O_BM,k  )));
    for(c = 0; c < cp9->el_from_ct[k]; c++)
      sc = ILogsum(sc, prv->el[cp9->el_from_idx[k][c]]);
    cur->m[k] = sc + msc[k];

    /* insert state */
    sc = ILogsum(ILogsum(prv->m[k] + CP9TSC(cp9O_MI,k),
                         prv->i[k] + CP9TSC(cp9O_II,k)),
                 prv->d[k] + CP9TSC(cp9O_DI,k));
    cur->i[k] = sc + isc[k];

    /* delete state */
    sc = ILogsum(ILogsum(cur->m[k-1] + CP9TSC(cp9O_MD,k-1),
                         cur->i[k-1] + CP9TSC(cp9O_ID,k-1)),
                 cur->d[k-1] + CP9TSC(cp9O_DD,k-1));
    cur->d[k] = sc;

    /* EL state */
    sc = -INFTY;
    if((cp9->flags & CPLAN9_EL) && cp9->has_el[k])
      sc = ILogsum(cur->m[k]  + CP9TSC(cp9O_MEL,k),
                   prv->el[k] + el_selfsc);
    cur->el[k] = sc;
  }
}

/* Backward row L (residue i = j0). Same as the initialization of the L
 * row in cp9_Backward(). */
static void
cp9ck_bck_rowL(CP9_t *cp9, ESL_DSQ *dsq, int i, cp9ck_row_t *cur)
{
  int const *tsc = cp9->otsc;
  int M = cp9->M;
  int k, c;

  for (k = 1; k <= M; k++)
    cur->el[k] = -INFTY;
  if(cp9->flags & CPLAN9_EL) {
    for(c = 0; c < cp9->el_from_ct[M+1]; c++)
      cur->el[cp9->el_from_idx[M+1][c]] = 0.;
  }

  cur->m[M]  = 0. +
    ILogsum(cur->el[M] + CP9TSC(cp9O_MEL,M),
            CP9TSC(cp9O_ME,M));
  cur->m[M] += cp9->msc[dsq[i]][M];
  cur->i[M]  = 0. + CP9TSC(cp9O_IM,M);
  cur->i[M] += cp9->isc[dsq[i]][M];
  cur->d[M]  = CP9TSC(cp9O_DM,M);

  for (k = M-1; k >= 1; k--) {
    cur->m[k]  = 0 + CP9TSC(cp9O_ME,k);
    cur->m[k]  = ILogsum(cur->m[k], cur->d[k+1] + CP9TSC(cp9O_MD,k));
    if(cp9->flags & CPLAN9_EL)
      cur->m[k]  = ILogsum(cur->m[k], cur->el[k] + CP9TSC(cp9O_MEL,k));
    cur->m[k] += cp9->msc[dsq[i]][k];

    cur->i[k]  = cur->d[k+1] + CP9TSC(cp9O_ID,k);
    cur->i[k] += cp9->isc[dsq[i]][k];

    cur->d[k]  = cur->d[k+1] + CP9TSC(cp9O_DD,k);
  }

  cur->m[0]  = cur->d[1] + CP9TSC(cp9O_MD,0);
  cur->i[0]  = cur->d[1] + CP9TSC(cp9O_ID,0);
  cur->i[0] += cp9->isc[dsq[i]][0];
  cur->d[0]  = -INFTY;
  cur->el[0] = -INFTY;
}

/* Backward row for residue i (i0 <= i < j0), from the row for i+1. Same
 * as one iteration of the main loop of cp9_Backward() with do_scan
 * FALSE. */
static void
cp9ck_bck_row(CP9_t *cp9, ESL_DSQ *dsq, int i, cp9ck_row_t *prv, cp9ck_row_t *cur)
{
  int const *tsc = cp9->otsc;
  int M = cp9->M;
  int k, c;

  for (k = 0; k <= M; k++) cur->el[k] = -INFTY;

  if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
    cur->el[M] = cur->el[M] + cp9->el_selfsc;

  cur->m[M]  = prv->i[M] + CP9TSC(cp9O_MI,M);
  cur->m[M] += cp9->msc[dsq[i]][M];

  if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
    cur->m[M] = ILogsum(cur->m[M], cur->el[M] + CP9TSC(cp9O_MEL,M));

  cur->i[M]  = prv->i[M] + CP9TSC(cp9O_II,M);
  cur->i[M] += cp9->isc[dsq[i]][M];

  cur->d[M]  = prv->i[M] + CP9TSC(cp9O_DI,M);

  if(cp9->flags & CPLAN9_EL) {
    for(c = 0; c < cp9->el_from_ct[M]; c++)
      cur->el[cp9->el_from_idx[M][c]] = ILogsum(cur->el[cp9->el_from_idx[M][c]], prv->m[M]);
  }

  for (k = M-1; k >= 1; k--) {
    if(cp9->flags & CPLAN9_EL) {
      for(c = 0; c < cp9->el_from_ct[k]; c++)
        cur->el[cp9->el_from_idx[k][c]] = ILogsum(cur->el[cp9->el_from_idx[k][c]], prv->m[k]);
    }
    if((cp9->flags & CPLAN9_EL) && (cp9->has_el[k]))
      cur->el[k] = ILogsum(cur->el[k], prv->el[k] + cp9->el_selfsc);

    cur->m[k]  = ILogsum(ILogsum((prv->m[k+1] + CP9TSC(cp9O_MM,k)),
                                 (prv->i[k]   + CP9TSC(cp9O_MI,k))),
                         (cur->d[k+1] + CP9TSC(cp9O_MD,k)));
    if((cp9->flags & CPLAN9_EL) && (cp9->has_el[k]))
      cur->m[k] = ILogsum(cur->m[k], cur->el[k] + CP9TSC(cp9O_MEL,k));
    cur->m[k] += cp9->msc[dsq[i]][k];

    cur->i[k]  = ILogsum(ILogsum((prv->m[k+1] + CP9TSC(cp9O_IM,k)),
                                 (prv->i[k]   + CP9TSC(cp9O_II,k))),
                         (cur->d[k+1] + CP9TSC(cp9O_ID,k)));
    cur->i[k] += cp9->isc[dsq[i]][k];

    cur->d[k]  = ILogsum(ILogsum((prv->m[k+1] + CP9TSC(cp9O_DM,k)),
                                 (prv->i[k]   + CP9TSC(cp9O_DI,k))),
                         (cur->d[k+1] + CP9TSC(cp9O_DD,k)));
  }
  /* k == 0 */
  cur->i[0]  = ILogsum(ILogsum((prv->m[1] + CP9TSC(cp9O_IM,0)),
                               (prv->i[0] + CP9TSC(cp9O_II,0))),
                       (cur->d[1] + CP9TSC(cp9O_ID,0)));
  cur->i[0] += cp9->isc[dsq[i]][0];
  cur->d[0]  = -INFTY;
  cur->el[0] = -INFTY;

  cur->m[0] = -INFTY;
  for (k = M; k >= 1; k--)
    cur->m[0] = ILogsum(cur->m[0], (prv->m[k] + CP9TSC(cp9O_BM,k)));
  cur->m[0] = ILogsum(cur->m[0], (prv->i[0] + CP9TSC(cp9O_MI,0)));
  cur->m[0] = ILogsum(cur->m[0], (cur->d[1] + CP9TSC(cp9O_MD,0)));
}

/* Backward row 0 (position i0-1), from row 1. Same as the special case
 * ip == 0 of cp9_Backward() with do_scan FALSE. */
static void
cp9ck_bck_row0(CP9_t *cp9, cp9ck_row_t *prv, cp9ck_row_t *cur)
{
  int const *tsc = cp9->otsc;
  int M = cp9->M;
  int k;

  for (k = 1; k <= M; k++) cur->el[k] = -INFTY;

  cur->m[M]  = -INFTY;
  cur->i[M]  = -INFTY;
  cur->el[M] = -INFTY;
  cur->d[M]  = prv->i[M] + CP9TSC(cp9O_DI,M);
  for (k = M-1; k >= 1; k--) {
    cur->m[k]  = -INFTY;
    cur->i[k]  = -INFTY;
    cur->el[k] = -INFTY;
    cur->d[k]  = ILogsum(ILogsum((prv->m[k+1] + CP9TSC(cp9O_DM,k)),
                                 (prv->i[k]   + CP9TSC(cp9O_DI,k))),
                         (cur->d[k+1] + CP9TSC(cp9O_DD,k)));
  }
  cur->i[0]  = -INFTY;
  cur->d[0]  = -INFTY;
  cur->el[0] = -INFTY;

  cur->m[0] = -INFTY;
  for (k = M; k >= 1; k--)
    cur->m[0] = ILogsum(cur->m[0], (prv->m[k] + CP9TSC(cp9O_BM,k)));
  cur->m[0] = ILogsum(cur->m[0], (prv->i[0] + CP9TSC(cp9O_MI,0)));
  cur->m[0] = ILogsum(cur->m[0], (cur->d[1] + CP9TSC(cp9O_MD,0)));
}

/* Backward row r (r = 0..L, residue i0+r-1) from row r+1, choosing the
 * right one of the three Backward row functions. */
static void
cp9ck_bck_any(CP9_t *cp9, ESL_DSQ *dsq, int i0, int L, int r, cp9ck_row_t *prv, cp9ck_row_t *cur)
{
  if     (r == L) cp9ck_bck_rowL(cp9, dsq, i0+r-1, cur);
  else if(r == 0) cp9ck_bck_row0(cp9, prv, cur);
  else            cp9ck_bck_row (cp9, dsq, i0+r-1, prv, cur);
}

/* Forward row r (r = 0..L, residue i0+r-1) from row r-1. */
static void
cp9ck_fwd_any(CP9_t *cp9, ESL_DSQ *dsq, int i0, int r, cp9ck_row_t *prv, cp9ck_row_t *cur)
{
  if(r == 0) cp9ck_fwd_row0(cp9, cur);
  else       cp9ck_fwd_row (cp9, dsq, i0+r-1, prv, cur);
}

/* Fold posterior <p> into cumulative mass <*mass> of one state and record,
 * at position <pos>, every threshold that the mass now exceeds for the
 * first time. <thr> is sorted ascending; <*ncross> is the number already
 * exceeded, and so the index of the smallest one not yet exceeded. Once
 * all are exceeded the mass is no longer needed. Since the thresholds are
 * ascending, a threshold above the smallest one not yet exceeded cannot
 * have been exceeded earlier, so each threshold is recorded at the first
 * position its own accumulated mass exceeds it, as in cp9_FB2HMMBands(). */
static inline void
cp9ck_fold(int p, int pos, int *mass, int *ncross, const int *thr, int nthr, int *pn, int stride)
{
  if(*ncross >= nthr) return;
  *mass = ILogsum(*mass, p);
  while(*ncross < nthr && *mass > thr[*ncross]) { pn[(*ncross) * stride] = pos; (*ncross)++; }
}

/* the same, for the boundary row where cp9_FB2HMMBands() assigns the
 * mass instead of summing into it */
static inline void
cp9ck_set(int p, int pos, int *mass, int *ncross, const int *thr, int nthr, int *pn, int stride)
{
  *mass = p;
  while(*ncross < nthr && *mass > thr[*ncross]) { pn[(*ncross) * stride] = pos; (*ncross)++; }
}

/* Function:  cp9_CheckptBandsCreate()
 * Incept:    EPN, Mon Oct 12 2026 (w/Claude)
 *
 * Purpose:   Create a CP9_CKPTBANDS object for a CP9 HMM of <M> nodes and
 *            the <nthr> integer posterior mass thresholds in <thr> (any
 *            order, duplicates allowed). It will hold, for each distinct
 *            threshold, the minimum and maximum position of every HMM
 *            state's band.
 *
 * Returns:   the new object, or NULL on allocation failure.
 */
CP9_CKPTBANDS *
cp9_CheckptBandsCreate(int M, const int *thr, int nthr)
{
  int status;
  CP9_CKPTBANDS *cb = NULL;
  int a, b, t, n;
  int64_t ncells;

  ESL_ALLOC(cb, sizeof(CP9_CKPTBANDS));
  cb->M        = M;
  cb->nthr     = 0;
  cb->thr      = NULL;
  cb->pn_min_m = cb->pn_max_m = NULL;
  cb->pn_min_i = cb->pn_max_i = NULL;
  cb->pn_min_d = cb->pn_max_d = NULL;
  cb->occ      = NULL;
  cb->sc       = 0;

  /* sort ascending, remove duplicates */
  ESL_ALLOC(cb->thr, sizeof(int) * ESL_MAX(nthr, 1));
  for(a = 0; a < nthr; a++) cb->thr[a] = thr[a];
  for(a = 1; a < nthr; a++) {
    t = cb->thr[a];
    for(b = a-1; b >= 0 && cb->thr[b] > t; b--) cb->thr[b+1] = cb->thr[b];
    cb->thr[b+1] = t;
  }
  n = 0;
  for(a = 0; a < nthr; a++) if(n == 0 || cb->thr[a] != cb->thr[n-1]) cb->thr[n++] = cb->thr[a];
  cb->nthr = n;

  ncells = (int64_t) cb->nthr * (M+1);
  ESL_ALLOC(cb->pn_min_m, sizeof(int) * ESL_MAX(ncells, 1));
  ESL_ALLOC(cb->pn_max_m, sizeof(int) * ESL_MAX(ncells, 1));
  ESL_ALLOC(cb->pn_min_i, sizeof(int) * ESL_MAX(ncells, 1));
  ESL_ALLOC(cb->pn_max_i, sizeof(int) * ESL_MAX(ncells, 1));
  ESL_ALLOC(cb->pn_min_d, sizeof(int) * ESL_MAX(ncells, 1));
  ESL_ALLOC(cb->pn_max_d, sizeof(int) * ESL_MAX(ncells, 1));
  ESL_ALLOC(cb->occ,      sizeof(int) * (M+1));
  return cb;

 ERROR:
  cp9_CheckptBandsDestroy(cb);
  return NULL;
}

/* Function:  cp9_CheckptBandsDestroy()
 * Incept:    EPN, Mon Oct 12 2026 (w/Claude)
 */
void
cp9_CheckptBandsDestroy(CP9_CKPTBANDS *cb)
{
  if(cb == NULL) return;
  if(cb->thr)      free(cb->thr);
  if(cb->pn_min_m) free(cb->pn_min_m);
  if(cb->pn_max_m) free(cb->pn_max_m);
  if(cb->pn_min_i) free(cb->pn_min_i);
  if(cb->pn_max_i) free(cb->pn_max_i);
  if(cb->pn_min_d) free(cb->pn_min_d);
  if(cb->pn_max_d) free(cb->pn_max_d);
  if(cb->occ)      free(cb->occ);
  free(cb);
}

/* Function:  cp9_CheckptBandsThreshIdx()
 * Incept:    EPN, Mon Oct 12 2026 (w/Claude)
 *
 * Purpose:   Return the index in <cb> of threshold <thresh>, or -1 if
 *            <cb> was not created with it.
 */
int
cp9_CheckptBandsThreshIdx(CP9_CKPTBANDS *cb, int thresh)
{
  int t;
  for(t = 0; t < cb->nthr; t++) if(cb->thr[t] == thresh) return t;
  return -1;
}

/* segment length, in DP rows, for a sequence of length L */
static int
cp9ck_seglen(int L)
{
  int R = (int) ceil(sqrt(2. * (double) (L+1)));
  if(R < 1)   R = 1;
  if(R > L+1) R = L+1;
  return R;
}

/* Function:  cp9_CheckptBandsMbNeeded()
 * Incept:    EPN, Mon Oct 12 2026 (w/Claude)
 *
 * Purpose:   Return the number of Mb cp9_CheckptFB2HMMBands() needs for
 *            a target sequence of length <L>, a CP9 HMM of <M> nodes and
 *            <nthr> thresholds (the CP9_CKPTBANDS object included).
 */
float
cp9_CheckptBandsMbNeeded(int L, int M, int nthr)
{
  int    R    = cp9ck_seglen(L);
  int    nseg = (L+1 + R-1) / R;
  double nrow = 2. * nseg + R + 4.;
  double b    = 0.;

  b += nrow * 4. * (M+1) * sizeof(int);                  /* saved, segment and rolling DP rows */
  b += (6. * (double) ESL_MAX(nthr, 1) + 1.) * (M+1) * sizeof(int); /* CP9_CKPTBANDS bands, occupancies */
  b += 6. * (M+1) * sizeof(int);                         /* masses and counts */
  return (float) (b / 1000000.);
}

/* Function:  cp9_CheckptFB2HMMBands()
 * Incept:    EPN, Mon Oct 12 2026 (w/Claude)
 *
 * Purpose:   Calculate the HMM bands that cp9_Forward(), cp9_Backward()
 *            and cp9_FB2HMMBands() would calculate for an alignment of
 *            target subsequence <i0>..<j0> of <dsq> to <cp9>, for every
 *            threshold in <cb>, and the node occupancies that
 *            cp9_PredictStartAndEndPositions() would calculate, in
 *            O(sqrt(L) * M) memory. See the comment at the top of this
 *            file.
 *
 *            For threshold index t and node k, on return
 *            cb->pn_min_m[t * (M+1) + k] is the minimum position of the
 *            match state band of node k, or -1 if the state's
 *            cumulative posterior mass from the left never exceeds
 *            cb->thr[t] (cp9_FB2HMMBands()'s nset_m[k] FALSE); likewise
 *            for pn_max_* (from the right) and the insert (*_i) and
 *            delete (*_d) states. cb->occ[k] is the summed posterior of
 *            node k's match and delete states over positions i0-1..j0,
 *            and cb->sc is the Backward score.
 *
 *            Unlike cp9_FB2HMMBands(), this does not reconcile the
 *            minimum and maximum of each band; the caller does that
 *            for the threshold it uses.
 *
 * Args:      cp9    - the CP9 HMM
 *            errbuf - for error messages
 *            dsq    - the digitized sequence
 *            i0, j0 - target subsequence, i0 <= j0
 *            cb     - created by cp9_CheckptBandsCreate() for cp9->M nodes; filled here
 *
 * Returns:   eslOK on success.
 *            eslEINCOMPAT on contract violation; eslEMEM on allocation
 *            failure; errbuf is filled.
 */
int
cp9_CheckptFB2HMMBands(CP9_t *cp9, char *errbuf, ESL_DSQ *dsq, int i0, int j0, CP9_CKPTBANDS *cb)
{
  int     status;
  int     M = cp9->M;
  int     L = j0-i0+1;
  int     R, nseg, nrows;
  int     s, a, b, r, k, i, pos0;
  int     sc;
  int     p_m, p_i, p_d;
  int     stride = M+1;
  int     nthr;
  const int *thr;
  int    *mem   = NULL;
  cp9ck_row_t *rows = NULL;
  cp9ck_row_t *Fck, *Bck, *seg, *froll, *broll; /* views into rows[] */
  cp9ck_row_t *F, *B;
  int    *mass_m = NULL, *mass_i = NULL, *mass_d = NULL;
  int    *nc_m   = NULL, *nc_i   = NULL, *nc_d   = NULL;
  int64_t n;

  if(cp9 == NULL)       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckptFB2HMMBands(), cp9 is NULL");
  if(dsq == NULL)       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckptFB2HMMBands(), dsq is NULL");
  if(cb == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckptFB2HMMBands(), cb is NULL");
  if(cb->M != M)        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckptFB2HMMBands(), cb->M %d != cp9->M %d", cb->M, M);
  if(i0 < 1 || L < 1)   ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckptFB2HMMBands(), i0: %d j0: %d invalid", i0, j0);
  nthr = cb->nthr;
  thr  = cb->thr;

  R     = cp9ck_seglen(L);
  nseg  = (L+1 + R-1) / R;
  nrows = 2*nseg + R + 4;
  ESL_ALLOC(rows, sizeof(cp9ck_row_t) * nrows);
  ESL_ALLOC(mem,  sizeof(int) * (int64_t) nrows * 4 * (M+1));
  for(r = 0; r < nrows; r++) {
    n = (int64_t) r * 4 * (M+1);
    rows[r].m  = mem + n;
    rows[r].i  = mem + n +   (M+1);
    rows[r].d  = mem + n + 2*(M+1);
    rows[r].el = mem + n + 3*(M+1);
  }
  Fck   = rows;             /* [0..nseg-1] Forward row a_s, the first row of segment s */
  Bck   = rows + nseg;      /* [0..nseg-1] Backward row b_s, the last row of segment s */
  seg   = rows + 2*nseg;    /* [0..R-1]    rows a_s..b_s of the segment being processed */
  froll = rows + 2*nseg + R;     /* [0..1] Forward rows, indexed by row parity */
  broll = rows + 2*nseg + R + 2; /* [0..1] Backward rows, indexed by row parity */

  ESL_ALLOC(mass_m, sizeof(int) * (M+1));
  ESL_ALLOC(mass_i, sizeof(int) * (M+1));
  ESL_ALLOC(mass_d, sizeof(int) * (M+1));
  ESL_ALLOC(nc_m,   sizeof(int) * (M+1));
  ESL_ALLOC(nc_i,   sizeof(int) * (M+1));
  ESL_ALLOC(nc_d,   sizeof(int) * (M+1));

  n = (int64_t) nthr * (M+1);
  for(r = 0; r < n; r++) {
    cb->pn_min_m[r] = cb->pn_max_m[r] = -1;
    cb->pn_min_i[r] = cb->pn_max_i[r] = -1;
    cb->pn_min_d[r] = cb->pn_max_d[r] = -1;
  }
  pos0 = ESL_MAX(i0-1, 0); /* position recorded for row 0 */

  /* 1. Backward, rows L..0, saving the last row of each segment. Row r
   *    is the last row of its segment if r == L or r % R == R-1. */
  for(r = L; r >= 0; r--) {
    cp9ck_bck_any(cp9, dsq, i0, L, r, &(broll[(r+1)%2]), &(broll[r%2]));
    if(r == L || (r % R) == R-1) cp9ck_copy_row(&(Bck[r / R]), &(broll[r%2]), M);
  }
  sc = broll[0].m[0]; /* summed score of all parses, as in cp9_FB2HMMBands() with did_fwd_scan FALSE */
  cb->sc = sc;

  /* 2. Left sweep: minimum position of each state's band, and node occupancies. */
  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  esl_vec_ISet(nc_m, M+1, 0);
  esl_vec_ISet(nc_i, M+1, 0);
  esl_vec_ISet(nc_d, M+1, 0);
  for(s = 0; s < nseg; s++) {
    a = s * R;
    b = ESL_MIN(L, a + R - 1);
    /* recompute this segment's Backward rows */
    cp9ck_copy_row(&(seg[b-a]), &(Bck[s]), M);
    for(r = b-1; r >= a; r--)
      cp9ck_bck_any(cp9, dsq, i0, L, r, &(seg[r+1-a]), &(seg[r-a]));
    for(r = a; r <= b; r++) {
      cp9ck_fwd_any(cp9, dsq, i0, r, &(froll[(r+1)%2]), &(froll[r%2]));
      if(r == a) cp9ck_copy_row(&(Fck[s]), &(froll[r%2]), M);
      F = &(froll[r%2]);
      B = &(seg[r-a]);
      if(r == 0) {
        /* boundary row, i = i0-1: only M_0 and the D states can be occupied */
        p_m = F->m[0] + B->m[0] - sc;
        cp9ck_set(p_m, pos0, &(mass_m[0]), &(nc_m[0]), thr, nthr, &(cb->pn_min_m[0]), stride);
        for(k = 1; k <= M; k++) {
          p_d = F->d[k] + B->d[k] - sc;
          cp9ck_set(p_d, pos0, &(mass_d[k]), &(nc_d[k]), thr, nthr, &(cb->pn_min_d[k]), stride);
          cb->occ[k] = ILogsum(-INFTY, ILogsum(-INFTY, p_d));
        }
      }
      else {
        i = i0 + r - 1;
        p_m = ESL_MAX(F->m[0] + B->m[0] - sc, -INFTY);
        cp9ck_fold(p_m, i, &(mass_m[0]), &(nc_m[0]), thr, nthr, &(cb->pn_min_m[0]), stride);
        p_i = ESL_MAX(F->i[0] + B->i[0] - cp9->isc[dsq[i]][0] - sc, -INFTY);
        cp9ck_fold(p_i, i, &(mass_i[0]), &(nc_i[0]), thr, nthr, &(cb->pn_min_i[0]), stride);
        for(k = 1; k <= M; k++) {
          p_m = ESL_MAX(F->m[k] + B->m[k] - cp9->msc[dsq[i]][k] - sc, -INFTY);
          p_i = ESL_MAX(F->i[k] + B->i[k] - cp9->isc[dsq[i]][k] - sc, -INFTY);
          p_d = ESL_MAX(F->d[k] + B->d[k] - sc, -INFTY);
          cp9ck_fold(p_m, i, &(mass_m[k]), &(nc_m[k]), thr, nthr, &(cb->pn_min_m[k]), stride);
          cp9ck_fold(p_i, i, &(mass_i[k]), &(nc_i[k]), thr, nthr, &(cb->pn_min_i[k]), stride);
          cp9ck_fold(p_d, i, &(mass_d[k]), &(nc_d[k]), thr, nthr, &(cb->pn_min_d[k]), stride);
          cb->occ[k] = ILogsum(cb->occ[k], ILogsum(p_m, p_d));
        }
      }
    }
  }
  cb->occ[0] = -INFTY; /* unused: node 0 has no occupancy */

  /* 3. Right sweep: maximum position of each state's band. */
  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  esl_vec_ISet(nc_m, M+1, 0);
  esl_vec_ISet(nc_i, M+1, 0);
  esl_vec_ISet(nc_d, M+1, 0);
  for(s = nseg-1; s >= 0; s--) {
    a = s * R;
    b = ESL_MIN(L, a + R - 1);
    /* recompute this segment's Forward rows */
    cp9ck_copy_row(&(seg[0]), &(Fck[s]), M);
    for(r = a+1; r <= b; r++)
      cp9ck_fwd_any(cp9, dsq, i0, r, &(seg[r-1-a]), &(seg[r-a]));
    for(r = b; r >= a; r--) {
      cp9ck_bck_any(cp9, dsq, i0, L, r, &(broll[(r+1)%2]), &(broll[r%2]));
      F = &(seg[r-a]);
      B = &(broll[r%2]);
      if(r > 0) {
        i = i0 + r - 1;
        p_m = ESL_MAX(F->m[0] + B->m[0] - sc, -INFTY);
        p_i = ESL_MAX(F->i[0] + B->i[0] - cp9->isc[dsq[i]][0] - sc, -INFTY);
        cp9ck_fold(p_m,    i, &(mass_m[0]), &(nc_m[0]), thr, nthr, &(cb->pn_max_m[0]), stride);
        cp9ck_fold(p_i,    i, &(mass_i[0]), &(nc_i[0]), thr, nthr, &(cb->pn_max_i[0]), stride);
        cp9ck_fold(-INFTY, i, &(mass_d[0]), &(nc_d[0]), thr, nthr, &(cb->pn_max_d[0]), stride); /* D_0 doesn't exist */
        for(k = 1; k <= M; k++) {
          p_m = ESL_MAX(F->m[k] + B->m[k] - cp9->msc[dsq[i]][k] - sc, -INFTY);
          p_i = ESL_MAX(F->i[k] + B->i[k] - cp9->isc[dsq[i]][k] - sc, -INFTY);
          p_d = ESL_MAX(F->d[k] + B->d[k] - sc, -INFTY);
          cp9ck_fold(p_m, i, &(mass_m[k]), &(nc_m[k]), thr, nthr, &(cb->pn_max_m[k]), stride);
          cp9ck_fold(p_i, i, &(mass_i[k]), &(nc_i[k]), thr, nthr, &(cb->pn_max_i[k]), stride);
          cp9ck_fold(p_d, i, &(mass_d[k]), &(nc_d[k]), thr, nthr, &(cb->pn_max_d[k]), stride);
        }
      }
      else {
        /* boundary row, i = i0-1 */
        p_m = F->m[0] + B->m[0] - sc;
        cp9ck_fold(p_m, pos0, &(mass_m[0]), &(nc_m[0]), thr, nthr, &(cb->pn_max_m[0]), stride);
        for(k = 1; k <= M; k++) {
          p_d = F->d[k] + B->d[k] - sc;
          cp9ck_fold(p_d, pos0, &(mass_d[k]), &(nc_d[k]), thr, nthr, &(cb->pn_max_d[k]), stride);
        }
      }
    }
  }

  free(mass_m); free(mass_i); free(mass_d);
  free(nc_m);   free(nc_i);   free(nc_d);
  free(rows);
  free(mem);
  return eslOK;

 ERROR:
  if(mass_m) free(mass_m);
  if(mass_i) free(mass_i);
  if(mass_d) free(mass_d);
  if(nc_m)   free(nc_m);
  if(nc_i)   free(nc_i);
  if(nc_d)   free(nc_d);
  if(rows)   free(rows);
  if(mem)    free(mem);
  ESL_FAIL(eslEMEM, errbuf, "cp9_CheckptFB2HMMBands(): out of memory");
}

/* Beginning modifications for a vectorized implementation ... */

/* sse_cm_dpsearch.c
 *
 * DP functions for CYK CM similarity search,
 * 4x single-precision float SSE implementation (non-striped)
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_vectorops.h"
#include "esl_sse.h"

#include "funcs.h"
#include "structs.h"
#include "impl_sse.h"

/* Function: SSE_CYKScan()
 * Date:     DLK, 17 June 2009
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           a reference CYK scanning algorithm.  
 *
 * Args:     cm              - the covariance model
 *           errbuf          - char buffer for reporting errors
 *           smx             - ScanMatrix_t for this search w/this model (incl. DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           do_null3        - TRUE to do NULL3 score correction, FALSE not to
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *           ret_sc          - RETURN: score of best overall hit (vsc[0])
 *
 * Note:     This function is heavily synchronized with RefIInsideScan() and RefCYKScan()
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  eslOK on succes;
 *           <ret_sc> is score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
int
SSE_CYKScan(CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, 
	   search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc)
{
  int       status;
  GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
  float    *vsc;                /* best score for each state (float) */
  float     vsc_root;           /* best overall score (score at ROOT_S) */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k, kp;		/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       W;                  /* max d; max size of a hit, this is min(L, smx->W) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  float   **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */
  double  **act;                /* [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1) */

  /* Contract check */
  if(! cm->flags & CMH_BITS)             ESL_FAIL(eslEINCOMPAT, errbuf, "SSE_CYKScan, CMH_BITS flag is not raised.\n");
  if(j0 < i0)                            ESL_FAIL(eslEINCOMPAT, errbuf, "SSE_CYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "SSE_CYKScan, dsq is NULL\n");
  if(smx == NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "SSE_CYKScan, smx == NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) ESL_FAIL(eslEINCOMPAT, errbuf, "SSE_CYKScan, CM_SEARCH_INSIDE flag raised");
  if(! (cm->smx->flags & cmSMX_HAS_FLOAT)) ESL_FAIL(eslEINCOMPAT, errbuf, "SSE_CYKScan, ScanMatrix's cmSMX_HAS_FLOAT flag is not raised");
  if(smx == cm->smx && (! cm->flags & CMH_SCANMATRIX)) ESL_FAIL(eslEINCOMPAT, errbuf, "SSE_CYKScan, smx == cm->smx, and cm->flags & CMH_SCANMATRIX is down, matrix is invalid.");

  /* make pointers to the ScanMatrix/CM data for convenience */
  float ***alpha      = smx->falpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  float ***alpha_begl = smx->falpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = smx->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = smx->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *bestr       = smx->bestr;       /* [0..d..W] best root state (for local begins or 0) for this d */

  /* SIMD vectors */
  int sW, delta, x;
  int *esc_stale;
  __m128   *mem_init_scAA;
  __m128  **vec_init_scAA;
  __m128   *mem_alpha;
  __m128 ***vec_alpha;
  __m128   *mem_alpha_begl;
  __m128 ***vec_alpha_begl;
  __m128   *mem_esc;
  __m128 ***vec_esc;
  __m128    zerov;
  __m128    neginfv;
  __m128    vec_tsc;
  __m128   *mem_bestr = NULL;
  __m128   *vec_bestr;
  __m128    vec_beginsc;
  __m128    mask;
  __m128    tmpv;
  __m128   *tmpary;
  __m128   *mem_tmpary;
  float     tmp_esc;
  const int vecwidth = 4;

  L = j0-i0+1;
  W = smx->W;
  if (W > L) W = L; 

  /* set vsc array */
  vsc = NULL;
  if(ret_vsc != NULL) { 
    ESL_ALLOC(vsc, sizeof(float) * cm->M);
    esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);
  }
  vsc_root    = IMPOSSIBLE;

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. */
  if(results != NULL) gamma = CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff, FALSE);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* precalculate the initial scores for all cells */
  init_scAA = FCalcInitDPScores(cm);

  /* if do_null3: allocate and initialize act vector */
  if(do_null3) { 
    ESL_ALLOC(act, sizeof(double *) * (W+1));
    for(i = 0; i <= W; i++) { 
      ESL_ALLOC(act[i], sizeof(double) * cm->abc->K);
      esl_vec_DSet(act[i], cm->abc->K, 0.);
    }
  }
  else act = NULL;

  /* Vectorize! */
  sW = W/4 + 1;
  zerov = _mm_setzero_ps();
  neginfv = _mm_set1_ps(-eslINFINITY);
  /* Memory allocation */
  ESL_ALLOC(vec_init_scAA,      sizeof(__m128  *) * cm->M);
  ESL_ALLOC(mem_init_scAA,      sizeof(__m128   ) * cm->M * (sW+1) + 15);
  ESL_ALLOC(vec_alpha,          sizeof(__m128 **) * 2);
  ESL_ALLOC(vec_alpha[0],       sizeof(__m128  *) * 2 * cm->M);
  ESL_ALLOC(mem_alpha,          sizeof(__m128   ) * 2 * cm->M * (sW+1) + 15);
  ESL_ALLOC(vec_alpha_begl,     sizeof(__m128 **) * (W+1));
  ESL_ALLOC(vec_alpha_begl[0],  sizeof(__m128  *) * (W+1) * cm->M);
  ESL_ALLOC(mem_alpha_begl,     sizeof(__m128   ) * (W+1) * cm->M * (sW+1) + 15);
  ESL_ALLOC(vec_esc,            sizeof(__m128 **) * cm->abc->Kp);
  ESL_ALLOC(vec_esc[0],         sizeof(__m128  *) * cm->abc->Kp * cm->M);
  ESL_ALLOC(mem_esc,            sizeof(__m128   ) * cm->abc->Kp * cm->M * (sW+1) + 15);
  ESL_ALLOC(esc_stale,          sizeof(int     *) * cm->abc->Kp);
  ESL_ALLOC(mem_tmpary,         sizeof(__m128   ) * (W+1) + 15);
  ESL_ALLOC(mem_bestr,          sizeof(__m128   ) * (sW+1) + 15);

  /* Set pointers for things that aren't directly allocated */
  vec_bestr            = (__m128 *) (((unsigned long int) mem_bestr + 15) & (~0xf));
  tmpary               = (__m128 *) (((unsigned long int) mem_tmpary + 15) & (~0xf));
  vec_init_scAA[0]     = (__m128 *) (((unsigned long int) mem_init_scAA + 15) & (~0xf));
  vec_alpha[1]         = vec_alpha[0] + cm->M;
  vec_alpha[0][0]      = (__m128 *) (((unsigned long int) mem_alpha + 15) & (~0xf));
  vec_alpha[1][0]      = vec_alpha[0][0] + cm->M * (sW+1);
  vec_alpha_begl[0][0] = (__m128 *) (((unsigned long int) mem_alpha_begl + 15) & (~0xf));
  vec_esc[0][0]        = (__m128 *) (((unsigned long int) mem_esc + 15) & (~0xf));
  esc_stale[0]         = 0;
  for (j = 1; j <= W; j++)
    {
      vec_alpha_begl[j]    = vec_alpha_begl[0]    + j*(cm->M);
      vec_alpha_begl[j][0] = vec_alpha_begl[0][0] + j * cm->M * (sW+1);
    }
  for (j = 1; j < cm->abc->Kp; j++)
    {
      vec_esc[j]    = vec_esc[0]    + j*(cm->M);
      vec_esc[j][0] = vec_esc[0][0] + j * cm->M * (sW+1);
      esc_stale[j]  = 0;
    }
  for (v = 1; v < cm->M; v++)
    {
      vec_init_scAA[v] = vec_init_scAA[0] + v*(sW+1);
      vec_alpha[0][v]  = vec_alpha[0][0]  + v*(sW+1);
      vec_alpha[1][v]  = vec_alpha[1][0]  + v*(sW+1);
      for (j = 0; j <=W; j++)
        vec_alpha_begl[j][v] = vec_alpha_begl[j][0] + v*(sW+1);
      for (j = 0; j < cm->abc->Kp; j++)
        vec_esc[j][v] = vec_esc[j][0] + v*(sW+1);
    }

  /* Initialize vectors */
  for (v = 0; v < cm->M; v++) {
    for (d = 0; d < sW; d++)
      {
        vec_init_scAA[v][d] = _mm_setr_ps(init_scAA[v][d*vecwidth  ], init_scAA[v][d*vecwidth+1],
                                          init_scAA[v][d*vecwidth+2], init_scAA[v][d*vecwidth+3]);
        if (cm->stid[v] == BEGL_S) {
          for (j = 0; j <= W; j++) {
              vec_alpha_begl[j][v][d] = _mm_setr_ps(alpha_begl[j][v][d*vecwidth  ],
                                                    alpha_begl[j][v][d*vecwidth+1],
                                                    alpha_begl[j][v][d*vecwidth+2],
                                                    alpha_begl[j][v][d*vecwidth+3]);
            }
          }
        else {
          vec_alpha[0][v][d] = _mm_setr_ps(alpha[0][v][d*vecwidth  ], alpha[0][v][d*vecwidth+1],
                                           alpha[0][v][d*vecwidth+2], alpha[0][v][d*vecwidth+3]);
          vec_alpha[1][v][d] = _mm_setr_ps(alpha[1][v][d*vecwidth  ], alpha[1][v][d*vecwidth+1],
                                           alpha[1][v][d*vecwidth+2], alpha[1][v][d*vecwidth+3]);
          }
      }
    /* last d, might be only partially filled up to W */
    d = sW;
    vec_init_scAA[v][d] = _mm_setr_ps((d*vecwidth     <= W) ? init_scAA[v][d*vecwidth  ] : -eslINFINITY,
                                      (d*vecwidth + 1 <= W) ? init_scAA[v][d*vecwidth+1] : -eslINFINITY,
                                      (d*vecwidth + 2 <= W) ? init_scAA[v][d*vecwidth+2] : -eslINFINITY,
                                      (d*vecwidth + 3 <= W) ? init_scAA[v][d*vecwidth+3] : -eslINFINITY);
    if (cm->stid[v] == BEGL_S) {
      for (j = 0; j <= W; j++)
        vec_alpha_begl[j][v][d] = _mm_setr_ps((d*vecwidth     <= W) ? alpha_begl[j][v][d*vecwidth  ] : -eslINFINITY,
                                              (d*vecwidth + 1 <= W) ? alpha_begl[j][v][d*vecwidth+1] : -eslINFINITY,
                                              (d*vecwidth + 1 <= W) ? alpha_begl[j][v][d*vecwidth+1] : -eslINFINITY,
                                              (d*vecwidth + 1 <= W) ? alpha_begl[j][v][d*vecwidth+1] : -eslINFINITY);
      }
    else {
        vec_alpha[0][v][d] = _mm_setr_ps((d*vecwidth      <= W) ? alpha[0][v][d*vecwidth  ] : -eslINFINITY,
                                         (d*vecwidth  + 1 <= W) ? alpha[0][v][d*vecwidth+1] : -eslINFINITY,
                                         (d*vecwidth  + 2 <= W) ? alpha[0][v][d*vecwidth+2] : -eslINFINITY,
                                         (d*vecwidth  + 3 <= W) ? alpha[0][v][d*vecwidth+3] : -eslINFINITY);
        vec_alpha[1][v][d] = _mm_setr_ps((d*vecwidth      <= W) ? alpha[1][v][d*vecwidth  ] : -eslINFINITY,
                                         (d*vecwidth  + 1 <= W) ? alpha[1][v][d*vecwidth+1] : -eslINFINITY,
                                         (d*vecwidth  + 2 <= W) ? alpha[1][v][d*vecwidth+2] : -eslINFINITY,
                                         (d*vecwidth  + 3 <= W) ? alpha[1][v][d*vecwidth+3] : -eslINFINITY);
      }
    }

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      jp_g = j-i0+1; /* j is actual index in dsq, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      /* if do_null3 (act != NULL), update act */
      if(act != NULL) { 
	esl_vec_DCopy(act[(jp_g-1)%(W+1)], cm->abc->K, act[jp_g%(W+1)]);
	esl_abc_DCount(cm->abc, act[jp_g%(W+1)], dsq[j], 1.);
	/*printf("j: %3d jp_g: %3d jp_g/W: %3d act[0]: %.3f act[1]: %.3f act[2]: %.3f act[3]: %.3f\n", j, jp_g, jp_g%(W+1), act[jp_g%(W+1)][0], act[jp_g%(W+1)][1], act[jp_g%(W+1)][2], act[jp_g%(W+1)][3]);*/
      }

      /* Build/increment pre-calc'd 3D matrix for esc */
      /* 1st dim: dsq[j] (1..abc->Kp)                 */
      /* 2nd dim: v      (0..cm->M-1)                 */
      /* 3rd dim: d      (0..min(j,W))                */
      delta = j>0 ? j - esc_stale[dsq[j]] : 0;
      if (esc_stale[dsq[j]] == 0) {
          /* Build from scratch */
          for (v = cm->M-1; v > 0; v--) {
            if (cm->sttype[v] == B_st || cm->sttype[v] == S_st || cm->sttype[v] == D_st || cm->sttype[v] == E_st)
              for (d = 0; d <= sW; d++)
                vec_esc[dsq[j]][v][d] = zerov;
            else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st)
              for (d = 0; d <= sW; d++)
                vec_esc[dsq[j]][v][d] = _mm_set1_ps(cm->oesc[v][dsq[j]]);
            else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st)
              for (d = 0; d <= sW; d++)
                vec_esc[dsq[j]][v][d] = _mm_setr_ps((d*vecwidth  <= j && ((j^j0)|d)) ? cm->oesc[v][dsq[j-(d*vecwidth  )+1]] : -eslINFINITY,
                                                    (d*vecwidth+1<= j              ) ? cm->oesc[v][dsq[j-(d*vecwidth+1)+1]] : -eslINFINITY,
                                                    (d*vecwidth+2<= j              ) ? cm->oesc[v][dsq[j-(d*vecwidth+2)+1]] : -eslINFINITY,
                                                    (d*vecwidth+3<= j              ) ? cm->oesc[v][dsq[j-(d*vecwidth+3)+1]] : -eslINFINITY);
            else if (cm->sttype[v] == MP_st)
              for (d = 0; d <= sW; d++)
                vec_esc[dsq[j]][v][d] = _mm_setr_ps((d*vecwidth  <= j && ((j^j0)|d)) ? cm->oesc[v][dsq[j-(d*vecwidth  )+1]*cm->abc->Kp+dsq[j]] : -eslINFINITY,
                                                    (d*vecwidth+1<= j              ) ? cm->oesc[v][dsq[j-(d*vecwidth+1)+1]*cm->abc->Kp+dsq[j]] : -eslINFINITY,
                                                    (d*vecwidth+2<= j              ) ? cm->oesc[v][dsq[j-(d*vecwidth+2)+1]*cm->abc->Kp+dsq[j]] : -eslINFINITY,
                                                    (d*vecwidth+3<= j              ) ? cm->oesc[v][dsq[j-(d*vecwidth+3)+1]*cm->abc->Kp+dsq[j]] : -eslINFINITY);
            else
              cm_Fail("SSE_CYKScan, building vec_esc, missed a case?");
          }
      }
      else if (delta == 1) {
        for (v = cm->M-1; v > 0; v--) {
          if (cm->sttype[v] == B_st || cm->sttype[v] == S_st  || cm->sttype[v] == D_st ||
              cm->sttype[v] == E_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) ;
          else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
            for (d = sW; d > 0; d--) 
              vec_esc[dsq[j]][v][d] = alt_rightshift_ps(vec_esc[dsq[j]][v][d], vec_esc[dsq[j]][v][d-1]);
            vec_esc[dsq[j]][v][0] = alt_rightshift_ps(vec_esc[dsq[j]][v][0], (j<j0) ? _mm_set1_ps(cm->oesc[v][dsq[j+1]]) : neginfv);
          }
          else if (cm->sttype[v] == MP_st) {
            for (d = sW; d > 0; d--) 
              vec_esc[dsq[j]][v][d] = alt_rightshift_ps(vec_esc[dsq[j]][v][d], vec_esc[dsq[j]][v][d-1]);
            vec_esc[dsq[j]][v][0] = alt_rightshift_ps(vec_esc[dsq[j]][v][0], (j<j0) ? _mm_set1_ps(cm->oesc[v][dsq[j+1]*cm->abc->Kp+dsq[j]]) : neginfv);
          }
        }
      }
      else if (delta == 2) {
        for (v = cm->M-1; v > 0; v--) {
          if (cm->sttype[v] == B_st || cm->sttype[v] == S_st  || cm->sttype[v] == D_st ||
              cm->sttype[v] == E_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) ;
          else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
            for (d = sW; d > 0; d--) {
              tmpv = _mm_movelh_ps(neginfv,vec_esc[dsq[j]][v][d]);
              vec_esc[dsq[j]][v][d] = _mm_movehl_ps(tmpv, vec_esc[dsq[j]][v][d-1]);
            }
            tmpv = _mm_setr_ps(j<j0 ? cm->oesc[v][dsq[j+1]] : -eslINFINITY, cm->oesc[v][dsq[j]], -eslINFINITY, -eslINFINITY);
            vec_esc[dsq[j]][v][0] = _mm_movelh_ps(tmpv, vec_esc[dsq[j]][v][0]);
          }
          else if (cm->sttype[v] == MP_st) {
            for (d = sW; d > 0; d--) {
              tmpv = _mm_movelh_ps(neginfv,vec_esc[dsq[j]][v][d]);
              vec_esc[dsq[j]][v][d] = _mm_movehl_ps(tmpv, vec_esc[dsq[j]][v][d-1]);
            }
            tmpv = _mm_setr_ps(j<j0 ? cm->oesc[v][dsq[j+1]*cm->abc->Kp+dsq[j]] : -eslINFINITY, cm->oesc[v][dsq[j]*cm->abc->Kp+dsq[j]], -eslINFINITY, -eslINFINITY);
            vec_esc[dsq[j]][v][0] = _mm_movelh_ps(tmpv, vec_esc[dsq[j]][v][0]);
          }
        }
      }
      else if (delta == 3) {
        for (v = cm->M-1; v > 0; v--) {
          if (cm->sttype[v] == B_st || cm->sttype[v] == S_st  || cm->sttype[v] == D_st ||
              cm->sttype[v] == E_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) ;
          else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
            for (d = sW; d > 0; d--) 
              vec_esc[dsq[j]][v][d] = esl_sse_leftshift_ps(vec_esc[dsq[j]][v][d-1], vec_esc[dsq[j]][v][d]);
            tmpv = _mm_setr_ps(-eslINFINITY, j<j0 ? cm->oesc[v][dsq[j+1]] : -eslINFINITY,
                                                    cm->oesc[v][dsq[j  ]],
                                                    cm->oesc[v][dsq[j-1]]);
            vec_esc[dsq[j]][v][0] = esl_sse_leftshift_ps(tmpv, vec_esc[dsq[j]][v][0]);
          }
          else if (cm->sttype[v] == MP_st) {
            for (d = sW; d > 0; d--) 
              vec_esc[dsq[j]][v][d] = esl_sse_leftshift_ps(vec_esc[dsq[j]][v][d-1], vec_esc[dsq[j]][v][d]);
            tmpv = _mm_setr_ps(-eslINFINITY, j<j0 ? cm->oesc[v][dsq[j+1]*cm->abc->Kp+dsq[j]] : -eslINFINITY,
                                                    cm->oesc[v][dsq[j  ]*cm->abc->Kp+dsq[j]],
                                                    cm->oesc[v][dsq[j-1]*cm->abc->Kp+dsq[j]]);
            vec_esc[dsq[j]][v][0] = esl_sse_leftshift_ps(tmpv, vec_esc[dsq[j]][v][0]);
          }
        }
      }
      else if (delta == 4) {
        for (v = cm->M-1; v > 0; v--) {
          if (cm->sttype[v] == B_st || cm->sttype[v] == S_st  || cm->sttype[v] == D_st ||
              cm->sttype[v] == E_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) ;
          else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
            for (d = sW; d > 0; d--) 
              vec_esc[dsq[j]][v][d] = vec_esc[dsq[j]][v][d-1];
            vec_esc[dsq[j]][v][0] = _mm_setr_ps(j<j0 ? cm->oesc[v][dsq[j+1]] : -eslINFINITY,
                                                       cm->oesc[v][dsq[j  ]],
                                                       cm->oesc[v][dsq[j-1]],
                                                       cm->oesc[v][dsq[j-2]]);
          }
          else if (cm->sttype[v] == MP_st) {
            for (d = sW; d > 0; d--) 
              vec_esc[dsq[j]][v][d] = vec_esc[dsq[j]][v][d-1];
            vec_esc[dsq[j]][v][0] = _mm_setr_ps(j<j0 ? cm->oesc[v][dsq[j+1]*cm->abc->Kp+dsq[j]] : -eslINFINITY,
                                                       cm->oesc[v][dsq[j  ]*cm->abc->Kp+dsq[j]],
                                                       cm->oesc[v][dsq[j-1]*cm->abc->Kp+dsq[j]],
                                                       cm->oesc[v][dsq[j-2]*cm->abc->Kp+dsq[j]]);
          }
        }
      }
      else cm_Fail("SSE_CYKScan, delta == %d",delta);
      esc_stale[dsq[j]] = j;

      /* Main loop */
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  float const *tsc_v = cm->tsc[v];
          __m128 vec_tmp_begl;
          __m128 vec_tmp_begr;

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  sd   = StateDelta(cm->sttype[v]);
	  cnum = cm->cnum[v];
	  dn   = dnA[v];
	  dx   = dxA[v];

	  if(cm->sttype[v] == B_st) {
	    w = cm->cfirst[v]; /* BEGL_S */
	    y = cm->cnum[v];   /* BEGR_S */

	    for (d = 0; d <= sW; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for B_st */
              vec_alpha[jp_v][v][d] = vec_init_scAA[v][d];
            }

/*
	    for (d = 0; d <= sW; d++) {
	      for (k = 0; k <= d; k++) 
		sc = ESL_MAX(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	    }
*/
            /* Case: k = 0 */
            k = 0;
            vec_tmp_begr = _mm_set1_ps(*((float *) &(vec_alpha[jp_y][y][0])));
            for (d = 0; d <= sW; d++) {
              vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_alpha_begl[jp_wA[k]][w][d],vec_tmp_begr));
            }
            /* Case: k = x*vecwidth */
            for (k = vecwidth; k <= W && k <= j; k+=vecwidth) {
              kp = k/vecwidth;
              vec_tmp_begr = _mm_set1_ps(*((float *) &(vec_alpha[jp_y][y][kp])));
              for (d = kp; d <= sW; d++) {
                vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_alpha_begl[jp_wA[k]][w][d-kp],vec_tmp_begr));
              }
            }
            /* Case: k = x*vecwidth+1 */
            for (k = 1; k <= W && k <= j; k+=vecwidth) {
              kp = k/vecwidth;
              vec_tmp_begr = _mm_set1_ps(*((float *) &(vec_alpha[jp_y][y][kp]) + 1));

              d = kp;
              vec_tmp_begl = esl_sse_rightshift_ps(vec_alpha_begl[jp_wA[k]][w][0], neginfv);
              vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_tmp_begl, vec_tmp_begr));
              
              for(d = kp+1; d <=sW; d++) {
                vec_tmp_begl = alt_rightshift_ps(vec_alpha_begl[jp_wA[k]][w][d-kp], vec_alpha_begl[jp_wA[k]][w][d-kp-1]);
                vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_tmp_begl, vec_tmp_begr));
              }
            }
            /* Case: k = x*vecwidth+2 */
            for (k = 2; k <= W && k <= j; k+=vecwidth) {
              kp = k/vecwidth;
              vec_tmp_begr = _mm_set1_ps(*((float *) &(vec_alpha[jp_y][y][kp]) + 2));

              d = kp;
              vec_tmp_begl = _mm_movelh_ps(neginfv, vec_alpha_begl[jp_wA[k]][w][0]);
              vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_tmp_begl, vec_tmp_begr));
              
              for(d = kp+1; d <=sW; d++) {
                vec_tmp_begl = _mm_movelh_ps(neginfv, vec_alpha_begl[jp_wA[k]][w][d-kp]);
                vec_tmp_begl = _mm_movehl_ps(vec_tmp_begl, vec_alpha_begl[jp_wA[k]][w][d-kp-1]);
                vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_tmp_begl, vec_tmp_begr));
              }
            }
            /* Case: k = x*vecwidth+3 */
            for (k = 3; k <= W && k <= j; k+=vecwidth) {
              kp = k/vecwidth;
              vec_tmp_begr = _mm_set1_ps(*((float *) &(vec_alpha[jp_y][y][kp]) + 3));

              d = kp;
              vec_tmp_begl = esl_sse_leftshift_ps(neginfv, vec_alpha_begl[jp_wA[k]][w][0]);
              vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_tmp_begl, vec_tmp_begr));
              
              for(d = kp+1; d <=sW; d++) {
                vec_tmp_begl = esl_sse_leftshift_ps(vec_alpha_begl[jp_wA[k]][w][d-kp-1], vec_alpha_begl[jp_wA[k]][w][d-kp]);
                vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_tmp_begl, vec_tmp_begr));
              }
            }
//fprintf(stderr,"j%2d v%2d ",j,v);
//for (d = 0; d <= sW && d <= j/4; d++) { vecprint_ps(vec_alpha[jp_v][v][d]); }
//fprintf(stderr,"\n");
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = 0; d <= sW; d++) {
	      vec_alpha_begl[jp_v][v][d] = vec_init_scAA[v][d]; /* state delta (sd) is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
                tmpv = _mm_add_ps(vec_alpha[jp_y][y+yoffset][d], _mm_set1_ps(tsc_v[yoffset]));
		vec_alpha_begl[jp_v][v][d] = _mm_max_ps(vec_alpha_begl[jp_v][v][d], tmpv);
              }
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
//fprintf(stderr,"j%2d v%2d ",j,v);
//for (d = 0; d <= sW && d <= j/4; d++) { vecprint_ps(vec_alpha_begl[jp_v][v][d]); }
//fprintf(stderr,"\n");
	  }
          else if (cm->sttype[v] == IL_st) { /* Special self-transition case */
            y = cm->cfirst[v];

            /* Initialize with vec_init_scAA */
            vec_alpha[jp_v][v][0] = alt_rightshift_ps(vec_init_scAA[v][0], neginfv);
            for (d = 1; d <= sW; d++) 
              vec_alpha[jp_v][v][d] = alt_rightshift_ps(vec_init_scAA[v][d], vec_init_scAA[v][d-1]);

            /* All children except the self-transition */
            for (yoffset = cm->cnum[v]-1; yoffset > 0; yoffset--) {
              vec_tsc = _mm_set1_ps(tsc_v[yoffset]);
              d = 0;
              tmpv = alt_rightshift_ps(vec_alpha[jp_y][y+yoffset][d], neginfv);
              vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(tmpv, vec_tsc));
              for (d = 1; d <= sW; d++) {
                tmpv = alt_rightshift_ps(vec_alpha[jp_y][y+yoffset][d], vec_alpha[jp_y][y+yoffset][d-1]);
                vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(tmpv, vec_tsc));
              }
            }
            for (d = 0; d <= sW; d++)
              vec_alpha[jp_v][v][d] = _mm_add_ps(vec_alpha[jp_v][v][d], vec_esc[dsq[j]][v][d]);

            /* Self-loop transition (yoffset = 0) */
            // FIXME: this is fully serialized - there's probably a lazy way to evaluate
            // FIXME: only if the values are actually updating
            vec_tsc = _mm_set1_ps(tsc_v[0]);
            d = 0;
            for (k = 2; k < vecwidth; k++) {
              tmpv = alt_rightshift_ps(vec_alpha[jp_y][y][d], neginfv);
              tmpv = _mm_add_ps(vec_esc[dsq[j]][v][d], _mm_add_ps(tmpv, vec_tsc));
              vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], tmpv);
            }
            for (d = 1; d <= sW; d++) {
              for (k = 0; k < vecwidth; k++) {
                tmpv = alt_rightshift_ps(vec_alpha[jp_y][y][d], vec_alpha[jp_y][y][d-1]);
                tmpv = _mm_add_ps(vec_esc[dsq[j]][v][d], _mm_add_ps(tmpv, vec_tsc));
                vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], tmpv);
              }
            }
//fprintf(stderr,"j%2d v%2d ",j,v);
//for (d = 0; d <= sW && d <= j/4; d++) { vecprint_ps(vec_alpha[jp_v][v][d]); }
//fprintf(stderr,"\n");
          }
	  else { /* ! B_st, ! BEGL_S st, ! IL_st */
	    y = cm->cfirst[v]; 
	    //i = j - dnA[v] + 1;

            switch (sd) {
              case 0:
                for (d = 0; d <= sW; d++)
                  vec_alpha[jp_v][v][d] = vec_init_scAA[v][d];
                for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
                  vec_tsc = _mm_set1_ps(tsc_v[yoffset]);
                  for (d = 0; d <= sW; d++) {
                    vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(vec_alpha[jp_y][y+yoffset][d], vec_tsc));
                  }
                }
                break;

              case 1:
                d = 0;
                vec_alpha[jp_v][v][d] = alt_rightshift_ps(vec_init_scAA[v][d], neginfv);
                for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
                  tmpv = alt_rightshift_ps(vec_alpha[jp_y][y+yoffset][d], neginfv);
                  vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(tmpv, _mm_set1_ps(tsc_v[yoffset])));
                }
                vec_alpha[jp_v][v][d] = _mm_add_ps(vec_alpha[jp_v][v][d], vec_esc[dsq[j]][v][d]);

                for (d = 1; d <= sW; d++) {
                  vec_alpha[jp_v][v][d] = alt_rightshift_ps(vec_init_scAA[v][d], vec_init_scAA[v][d-1]);
                  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
                    tmpv = alt_rightshift_ps(vec_alpha[jp_y][y+yoffset][d], vec_alpha[jp_y][y+yoffset][d-1]);
                    vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(tmpv, _mm_set1_ps(tsc_v[yoffset])));
                  }
                  vec_alpha[jp_v][v][d] = _mm_add_ps(vec_alpha[jp_v][v][d], vec_esc[dsq[j]][v][d]);
                }
                break;

              case 2:
                d = 0;
                vec_alpha[jp_v][v][d] = _mm_movelh_ps(neginfv, vec_init_scAA[v][d]);
                for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
                  tmpv = _mm_movelh_ps(neginfv, vec_alpha[jp_y][y+yoffset][d]);
                  vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(tmpv, _mm_set1_ps(tsc_v[yoffset])));
                }
                vec_alpha[jp_v][v][d] = _mm_add_ps(vec_alpha[jp_v][v][d], vec_esc[dsq[j]][v][d]);

                for (d = 1; d <= sW; d++) {
                  vec_alpha[jp_v][v][d] = _mm_movelh_ps(neginfv, vec_init_scAA[v][d]);
                  vec_alpha[jp_v][v][d] = _mm_movehl_ps(vec_alpha[jp_v][v][d], vec_init_scAA[v][d-1]);
                  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
                    tmpv = _mm_movelh_ps(neginfv, vec_alpha[jp_y][y+yoffset][d]);
                    tmpv = _mm_movehl_ps(tmpv, vec_alpha[jp_y][y+yoffset][d-1]);
                    vec_alpha[jp_v][v][d] = _mm_max_ps(vec_alpha[jp_v][v][d], _mm_add_ps(tmpv, _mm_set1_ps(tsc_v[yoffset])));
                  }
                  vec_alpha[jp_v][v][d] = _mm_add_ps(vec_alpha[jp_v][v][d], vec_esc[dsq[j]][v][d]);
                }
                break;

              default:
                cm_Fail("SSE_CYKScan; state delta out of range (sd = %d)",sd);
            }
//fprintf(stderr,"j%2d v%2d ",j,v);
//for (d = 0; d <= sW && d <= j/4; d++) { vecprint_ps(vec_alpha[jp_v][v][d]); }
//fprintf(stderr,"\n");
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) {
            tmpv = neginfv;
	    if(cm->stid[v] != BEGL_S) for (d = 0; d <= sW; d++) tmpv = _mm_max_ps(tmpv, vec_alpha[jp_v][v][d]);
	    else                      for (d = 0; d <= sW; d++) tmpv = _mm_max_ps(tmpv, vec_alpha_begl[jp_v][v][d]);
            esl_sse_hmax_ps(tmpv, &vsc[v]);
	  }

	} /*loop over decks v>0 */

      /* move esc vec's if their 'staleness' hits vecwidth */
      for (x = 0; x < cm->abc->Kp; x++) {
        if (esc_stale[x] == 0) { /* no values have been set at all; set d = 0 before sliding */
                                 /* and set the const-per-row emissions */
          for (v = cm->M-1; v > 0; v--) {
            if (cm->sttype[v] == B_st || cm->sttype[v] == S_st || cm->sttype[v] == D_st || cm->sttype[v] == E_st)
              for (d = 0; d <= sW; d++)
                vec_esc[x][v][d] = zerov;
            else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st)
              for (d = 0; d <= sW; d++)
                vec_esc[x][v][d] = _mm_set1_ps(cm->oesc[v][x]);
            else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st)
              vec_esc[x][v][0] = alt_rightshift_ps(neginfv, _mm_set1_ps(cm->oesc[v][dsq[1]]));
            else if (cm->sttype[v] == MP_st)
              vec_esc[x][v][0] = alt_rightshift_ps(neginfv, _mm_set1_ps(cm->oesc[v][dsq[1]*cm->abc->Kp+x]));
          }
        }

        delta = j - esc_stale[x];
        if (delta == 4) {
          for (v = cm->M-1; v > 0; v--) {
            if (cm->sttype[v] == B_st || cm->sttype[v] == S_st  || cm->sttype[v] == D_st ||
              cm->sttype[v] == E_st || cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) ;
            else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
              for (d = sW; d > 0; d--) 
                vec_esc[x][v][d] = vec_esc[x][v][d-1];
              vec_esc[x][v][0] = _mm_setr_ps(j<j0 ? cm->oesc[v][dsq[j+1]] : -eslINFINITY,
                                                    cm->oesc[v][dsq[j  ]],
                                                    cm->oesc[v][dsq[j-1]],
                                                    cm->oesc[v][dsq[j-2]]);
            }
            else if (cm->sttype[v] == MP_st) {
              for (d = sW; d > 0; d--) 
                vec_esc[x][v][d] = vec_esc[x][v][d-1];
              vec_esc[x][v][0] = _mm_setr_ps(j<j0 ? cm->oesc[v][dsq[j+1]*cm->abc->Kp+x] : -eslINFINITY,
                                                    cm->oesc[v][dsq[j  ]*cm->abc->Kp+x],
                                                    cm->oesc[v][dsq[j-1]*cm->abc->Kp+x],
                                                    cm->oesc[v][dsq[j-2]*cm->abc->Kp+x]);
            }
          }
          esc_stale[x] = j;
        }
      }
      
      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      float const *tsc_v = cm->tsc[0];
      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = 0; d <= sW; d++) {
	vec_bestr[d] = _mm_setzero_ps();	/* root of the traceback = root state 0 */
	y = cm->cfirst[0];
        vec_tsc = _mm_set1_ps(tsc_v[0]);
	vec_alpha[jp_v][0][d] = _mm_max_ps(neginfv, _mm_add_ps(vec_alpha[cur][y][d], vec_tsc));
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) {
          vec_tsc = _mm_set1_ps(tsc_v[yoffset]);
	  vec_alpha[jp_v][0][d] = _mm_max_ps(vec_alpha[jp_v][0][d], _mm_add_ps(vec_alpha[cur][y+yoffset][d], vec_tsc));
        }
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(NOT_IMPOSSIBLE(cm->beginsc[y])) {
            vec_beginsc = _mm_set1_ps(cm->beginsc[y]);
	    if(cm->stid[y] == BEGL_S)
	      {
		jp_y = j % (W+1);
		for (d = 0; d <= sW; d++) {
                  tmpv = _mm_add_ps(vec_alpha_begl[jp_y][y][d], vec_beginsc);
                  mask  = _mm_cmpgt_ps(tmpv, vec_alpha[jp_v][0][d]);
                  vec_alpha[jp_v][0][d] = _mm_max_ps(tmpv, vec_alpha[jp_v][0][d]);
                  vec_bestr[d] = esl_sse_select_ps(_mm_set1_ps((float) y), vec_bestr[d], mask);
		}
	      }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = 0; d <= sW; d++) {
		{
                  tmpv = _mm_add_ps(vec_alpha[jp_y][y][d], vec_beginsc);
                  mask  = _mm_cmpgt_ps(tmpv, vec_alpha[jp_v][0][d]);
                  vec_alpha[jp_v][0][d] = _mm_max_ps(tmpv, vec_alpha[jp_v][0][d]);
                  vec_bestr[d] = esl_sse_select_ps(_mm_set1_ps((float) y), vec_bestr[d], mask);
		}
	      }
	    }
	  }
	}
      }
      /* find the best score */
      tmpv = neginfv;
      for (d = 0; d < sW; d++)
        tmpv = _mm_max_ps(tmpv, vec_alpha[jp_v][0][d]);
      esl_sse_hmax_ps(tmpv, &tmp_esc);          /* overloaded tmp_esc, just need a temporary float vec */
      vsc_root = ESL_MAX(vsc_root, tmp_esc);
      /* for UpdateGammaHitMxCM to work, these data need to be un-vectorized: alpha[jp_v][0], bestr */
      if (results != NULL) {
        union vec_union { __m128 v; float x[4]; } utmp;
        for (d = 0; d < sW; d++) {
          utmp.v = vec_alpha[jp_v][0][d];
          alpha[jp_v][0][d*vecwidth  ] = utmp.x[0];
          alpha[jp_v][0][d*vecwidth+1] = utmp.x[1];
          alpha[jp_v][0][d*vecwidth+2] = utmp.x[2];
          alpha[jp_v][0][d*vecwidth+3] = utmp.x[3];
          utmp.v = vec_bestr[d];
          bestr[d*vecwidth  ] = utmp.x[0];
          bestr[d*vecwidth+1] = utmp.x[1];
          bestr[d*vecwidth+2] = utmp.x[2];
          bestr[d*vecwidth+3] = utmp.x[3];
        }
        utmp.v = vec_alpha[jp_v][0][sW];
        if ( sW*vecwidth   <= W ) alpha[jp_v][0][sW*vecwidth  ] = utmp.x[0];
        if ( sW*vecwidth+1 <= W ) alpha[jp_v][0][sW*vecwidth+1] = utmp.x[1];
        if ( sW*vecwidth+2 <= W ) alpha[jp_v][0][sW*vecwidth+2] = utmp.x[2];
        if ( sW*vecwidth+3 <= W ) alpha[jp_v][0][sW*vecwidth+3] = utmp.x[3];
        utmp.v = vec_bestr[sW];
        if ( sW*vecwidth   <= W ) bestr[sW*vecwidth  ] = utmp.x[0];
        if ( sW*vecwidth+1 <= W ) bestr[sW*vecwidth+1] = utmp.x[1];
        if ( sW*vecwidth+2 <= W ) bestr[sW*vecwidth+2] = utmp.x[2];
        if ( sW*vecwidth+3 <= W ) bestr[sW*vecwidth+3] = utmp.x[3];
      }
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) if((status = UpdateGammaHitMxCM(cm, errbuf, gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], FALSE, smx->bestr, results, W, act)) != eslOK) return status;
      /* cm_DumpScanMatrixAlpha(cm, si, j, i0, TRUE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in UpdateGammaHitMxCM() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    TBackGammaHitMxForward(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) FreeGammaHitMx(gamma);
  if (act != NULL) { 
    for(i = 0; i <= W; i++) free(act[i]); 
    free(act);
  }
  free(jp_wA);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc         = vsc;
  else free(vsc);
  if (ret_sc != NULL) *ret_sc = vsc_root;
  free(vec_init_scAA);     free(mem_init_scAA);
  free(vec_alpha[0]);      free(vec_alpha);      free(mem_alpha);
  free(vec_alpha_begl[0]); free(vec_alpha_begl); free(mem_alpha_begl);
  free(vec_esc[0]);        free(vec_esc);        free(mem_esc);
  free(esc_stale);
  free(mem_bestr);
  free(mem_tmpary);

  ESL_DPRINTF1(("SSE_CYKScan() return score: %10.4f\n", vsc_root)); 
  return eslOK;
  
 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_SEARCH_BENCHMARK
/* gcc -o sse-bmark -g -O2 -std=gnu99 -msse2 -mpentiumpro -I../ -L../ -I../../easel -L../../easel -I../../hmmer/src -L../../hmmer/src -DIMPL_SEARCH_BENCHMARK cm_dpsearch.c -linfernal -lhmmer -leasel -lm
 * icc -o sse-bmark -g -O3 -static -I../ -L../ -I../../easel -L../../easel -I../../hmmer/src -L../../hmmer/src -DIMPL_SEARCH_BENCHMARK sse_cm_dpsearch.c -linfernal -lhmmer -leasel -lm

 * Not updated for this file ...
 * mpicc -g -O2 -DHAVE_CONFIG_H -I../easel  -c old_cm_dpsearch.c 
 * mpicc -o benchmark-search -g -O2 -I. -L. -I../easel -L../easel -DIMPL_SEARCH_BENCHMARK cm_dpsearch.c old_cm_dpsearch.o -linfernal -leasel -lm
 * ./benchmark-search <cmfile>
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_randomseq.h>
#include <esl_sqio.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "funcs.h"		/* function declarations                */
/*     #include "old_funcs.h"	*/	/* function declarations for 0.81 versions */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "33", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-e",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "emit sequences from CM, don't randomly create them", 0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search in glocal mode [default: local]", 0 },
  { "-L",        eslARG_INT,  "10000", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-w",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute new reference CYK scan implementation", 0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute experimental CYK scan implementation", 0 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute non-banded optimized CYK scan implementation", 0 },
  { "--rsearch", eslARG_NONE,   FALSE, NULL, NULL,  NULL,"--noqdb", NULL, "also execute ported RSEARCH's CYK scan implementation", 0 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute optimized int inside scan implementation", 0 },
  { "--riins",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute reference int inside scan implementation", 0 },
  { "--fins",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute optimized float inside scan implementation", 0 },
  { "--rfins",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute reference float inside scan implementation", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "benchmark driver for an optimized scanning CYK implementation";

int 
main(int argc, char **argv)
{
  int             status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  CM_t            *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq;
  int             i;
  float           sc;
  char            *cmfile = esl_opt_GetArg(go, 1);
  CMFILE          *cmfp;	/* open input CM file stream */
  int            *dmin;
  int            *dmax;
  int             do_random;
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, either randomly created, or emitted from CM (if -e) */
  char            errbuf[cmERRBUFSIZE];

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if ((status = CMFileRead(cmfp, errbuf, &abc, &cm) != eslOK))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  if(! esl_opt_GetBoolean(go, "-g")) cm->config_opts  |= CM_CONFIG_LOCAL;
  cm->search_opts |= CM_SEARCH_NOQDB;
  ConfigCM(cm, errbuf, TRUE, NULL, NULL); /* TRUE says: calculate W */

  dmin = NULL; dmax = NULL;

  cm_CreateScanMatrixForCM(cm, TRUE, TRUE); /* impt to do this after QDBs set up in ConfigCM() */

  /* get sequences */
  if(!esl_opt_IsDefault(go, "-L")) {
     double *dnull;
     ESL_DSQ *randdsq = NULL;
     ESL_ALLOC(randdsq, sizeof(ESL_DSQ)* (L+2));
     ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
     for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
     esl_vec_DNorm(dnull, cm->abc->K);
     seqs_to_aln = CreateSeqsToAln(N, FALSE);

     for (i = 0; i < N; i++) {
       if (esl_rsq_xIID(r, dnull, cm->abc->K, L, randdsq)  != eslOK) cm_Fail("Failure creating random sequence.");
       if((seqs_to_aln->sq[i] = esl_sq_CreateDigitalFrom(abc, NULL, randdsq, L, NULL, NULL, NULL)) == NULL)
         cm_Fail("Failure digitizing/copying random sequence.");
     }
     seqs_to_aln->nseq = N;
     free(dnull);
     free(randdsq);
  }
  else if(do_random) {
    double *dnull;
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);
    /* get gamma[0] from the QDB calc alg, which will serve as the length distro for random seqs */
    int safe_windowlen = cm->clen * 2;
    double **gamma = NULL;
    while(!(BandCalculationEngine(cm, safe_windowlen, DEFAULT_HS_BETA, TRUE, NULL, NULL, &(gamma), NULL))) {
      safe_windowlen *= 2;
      FreeBandDensities(cm, gamma);
      if(safe_windowlen > (cm->clen * 1000)) cm_Fail("Error trying to get gamma[0], safe_windowlen big: %d\n", safe_windowlen);
     }
    seqs_to_aln = RandomEmitSeqsToAln(r, cm->abc, dnull, 1, N, gamma[0], safe_windowlen, FALSE);
    FreeBandDensities(cm, gamma);
    free(dnull);
  }
  else /* don't randomly generate seqs, emit them from the CM */
    seqs_to_aln = CMEmitSeqsToAln(r, cm, 1, N, FALSE, NULL, FALSE);

  for (i = 0; i < N; i++)
    {
      L = seqs_to_aln->sq[i]->n;
      dsq = seqs_to_aln->sq[i]->dsq;
      cm->search_opts  &= ~CM_SEARCH_INSIDE;

      esl_stopwatch_Start(w);
      if((status = FastCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

/*
      esl_stopwatch_Start(w);
      if((status = RefCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "RefCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
*/

      if (esl_opt_GetBoolean(go, "-w")) 
	{ 
	  esl_stopwatch_Start(w);
	  if((status = SSE_CYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	  printf("%4d %-30s %10.4f bits ", (i+1), "SSE_CYKScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}
  
      if (esl_opt_GetBoolean(go, "--rsearch")) 
	{ 
	  esl_stopwatch_Start(w);
	  if((status = rsearch_CYKScan (cm, errbuf, dsq, L, 0., cm->W, NULL, &sc)) != eslOK) cm_Fail(errbuf); 
	  printf("%4d %-30s %10.4f bits ", (i+1), "rsearch_CYKScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      /* integer inside implementations */
      if (esl_opt_GetBoolean(go, "--iins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	  printf("%4d %-30s %10.4f bits ", (i+1), "FastIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");

	  esl_stopwatch_Start(w);
	  if((status = XFastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	  printf("%4d %-30s %10.4f bits ", (i+1), "XFastIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");

	  esl_stopwatch_Start(w);
	  if((status = X2FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);;
	  printf("%4d %-30s %10.4f bits ", (i+1), "X2FastIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--riins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  if((status = RefIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	  printf("%4d %-30s %10.4f bits ", (i+1), "RefIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");

	  esl_stopwatch_Start(w);
	  if((status = XRefIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	  printf("%4d %-30s %10.4f bits ", (i+1), "XRefIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}
  
      /* float inside implementations */
      if (esl_opt_GetBoolean(go, "--fins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  if((status = FastFInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	  printf("%4d %-30s %10.4f bits ", (i+1), "FastFInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--rfins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  if((status = RefFInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	  printf("%4d %-30s %10.4f bits ", (i+1), "RefFInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}
  
      printf("\n");
    }
  FreeCM(cm);
  FreeSeqsToAln(seqs_to_aln);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  cm_Fail("memory allocation error");
  return 0; /* never reached */
}
#endif /*IMPL_SEARCH_BENCHMARK*/


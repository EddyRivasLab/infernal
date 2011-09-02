/* cm_dpsearch_trunc.c
 *
 * DP functions for truncated CYK and Inside CM similarity search.
 * 
 * RefTrCYKScan():  reference implementation of a scanning version
 *                  trCYK [Kolbe, Eddy 2009]. No FastTrCYKScan()
 *                  exists. I wrote one based on FastCYKScan() but
 *                  it was only about 5% faster and three times
 *                  as many lines of code as RefTrCYKScan(), so 
 *                  I scrapped it. Its in the subversion repository
 *                  though: r3663.
 *                  ref: ~nawrockie/notebook/11_0816_inf_banded_trcyk/00LOG
 *
 * EPN, Tue Aug 16 04:15:32 2011
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

#define PRINTFALPHA 0
#define PRINTIALPHA 0

/* Function: RefTrCYKScan()
 * Date:     EPN, Tue Aug 16 04:16:03 2011
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           a reference trCYK scanning algorithm. Query-dependent 
 *           bands are used or not used as specified in ScanMatrix_t <si>.
 *
 * Args:     cm              - the covariance model
 *           errbuf          - char buffer for reporting errors
 *           smx             - ScanMatrix_t for this search w/this model (incl. DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           cutoff          - minimum score to report
 *           th              - CM_TOPHITS to add to; if NULL, don't add to it
 *           do_null3        - TRUE to do NULL3 score correction, FALSE not to
 *           env_cutoff      - ret_envi..ret_envj will include all hits that exceed this bit sc
 *           ret_envi        - min position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted
 *           ret_envj        - max position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted 
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *           ret_sc          - RETURN: score of best overall hit (vsc[0])
 *
 * Note:     This function is heavily synchronized with RefITrInsideScan()
 *           any change to this function should be mirrored in that function.. 
 *
 * Returns:  eslOK on succes;
 *           <ret_sc> is score of best overall hit (vsc[0]). Information on hits added to <hitlist>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
int
RefTrCYKScan(CM_t *cm, char *errbuf, TrScanMatrix_t *trsmx, ESL_DSQ *dsq, int i0, int j0, float cutoff, CM_TOPHITS *hitlist,
	     int do_null3, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float **ret_vsc, float *ret_sc)
{
  int       status;
  GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
  float    *vsc;                /* best score for each state (float) */
  float     vsc_root;           /* best overall score (score at ROOT_S) */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jq_y;  	        /* offset j for state y plus 1 (if jp_y is prv, jq_y is cur, and vice versa) */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       W;                  /* max d; max size of a hit, this is min(L, smx->W) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       kn, kx;             /* minimum/maximum valid k for current d in B_st recursion */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  float   **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */
  double  **act;                /* [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1) */
  int       do_env_defn;        /* TRUE to calculate envi, envj, FALSE not to (TRUE if ret_envi != NULL or ret_envj != NULL */
  int64_t   envi, envj;         /* min/max positions that exist in any hit with sc >= env_cutoff */

  /* Contract check */
  if(! cm->flags & CMH_BITS)               ESL_FAIL(eslEINCOMPAT, errbuf, "RefTrCYKScan, CMH_BITS flag is not raised.\n");
  if(j0 < i0)                              ESL_FAIL(eslEINCOMPAT, errbuf, "RefTrCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                          ESL_FAIL(eslEINCOMPAT, errbuf, "RefTrCYKScan, dsq is NULL\n");
  if(trsmx == NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "RefTrCYKScan, trsmx == NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE)   ESL_FAIL(eslEINCOMPAT, errbuf, "RefTrCYKScan, CM_SEARCH_INSIDE flag raised");
  if(! (trsmx->flags & cmTRSMX_HAS_FLOAT)) ESL_FAIL(eslEINCOMPAT, errbuf, "RefTrCYKScan, ScanMatrix's cmTRSMX_HAS_FLOAT flag is not raised");

  /* make pointers to the ScanMatrix/CM data for convenience */
  float ***Jalpha      = trsmx->fJalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Jalpha DP matrix, NULL for v == BEGL_S */
  float ***Jalpha_begl = trsmx->fJalpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] Jalpha DP matrix, NULL for v != BEGL_S */
  float ***Lalpha      = trsmx->fLalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Lalpha DP matrix, NULL for v == BEGL_S */
  float ***Lalpha_begl = trsmx->fLalpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] Lalpha DP matrix, NULL for v != BEGL_S */
  float ***Ralpha      = trsmx->fRalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Ralpha DP matrix, NULL for v == BEGL_S */
  float ***Ralpha_begl = trsmx->fRalpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] Ralpha DP matrix, NULL for v != BEGL_S */
  float ***Talpha      = trsmx->fTalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Talpha DP matrix, NULL for v != BIF_B  */
  int    **dnAA        = trsmx->dnAA;         /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int    **dxAA        = trsmx->dxAA;         /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int     *bestr       = trsmx->bestr;        /* [0..d..W] best root state (for local begins or 0) for this d */
  int     *bestmode    = trsmx->bestmode;     /* [0..d..W] mode of best parsetree for this d */
  int     *dmax        = trsmx->dmax;         /* [0..v..cm->M-1] maximum d allowed for this state */
  float  **esc_vAA     = cm->oesc;            /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
 					       * and all possible emissions a (including ambiguities) */
  float  **lmesc_vAA   = cm->lmesc;           /* [0..v..cm->M-1][0..a..(cm->abc->Kp-1)] left  marginal emission scores for v */
  float  **rmesc_vAA   = cm->rmesc;           /* [0..v..cm->M-1][0..a..(cm->abc->Kp-1)] right marginal emission scores for v */

  /* determine if we're doing banded/non-banded */
  if(trsmx->dmax != NULL) do_banded = TRUE;
  
  L = j0-i0+1;
  W = trsmx->W;
  if (W > L) W = L; 

  /* set vsc array */
  vsc = NULL;
  ESL_ALLOC(vsc, sizeof(float) * cm->M);
  esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);
  vsc_root = IMPOSSIBLE;

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. */
  if(hitlist != NULL) gamma = CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff, FALSE);
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

  /* initialize envelope boundary variables */
  do_env_defn = (ret_envi != NULL || ret_envj != NULL) ? TRUE : FALSE;
  envi = j0+1;
  envj = i0-1;

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float Jsc, Lsc, Rsc, Tsc;
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

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  float const *esc_v   = esc_vAA[v]; 
	  float const *tsc_v   = cm->tsc[v];
	  float const *lmesc_v = lmesc_vAA[v]; 
	  float const *rmesc_v = rmesc_vAA[v]; 
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v  = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y  = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  jq_y = (StateRightDelta(cm->sttype[v]) > 0) ? cur : prv;
	  sd    = StateDelta(cm->sttype[v]);
	  cnum  = cm->cnum[v];
	  /* if we emit right, precalc score of emitting res j from state v */
	  float   esc_j = IMPOSSIBLE;
	  float rmesc_j = IMPOSSIBLE;
	  if(cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) { 
	    esc_j   =   esc_v[dsq[j]];
	    rmesc_j = rmesc_v[dsq[j]];
	  }
	  if(cm->sttype[v] == MP_st) { 
	    rmesc_j = rmesc_v[dsq[j]];
	  }

	  if(cm->sttype[v] == B_st) {
	    w = cm->cfirst[v]; /* BEGL_S */
	    y = cm->cnum[v];   /* BEGR_S */
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(0, (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], d);
	      }
	      else { kmin = 0; kmax = d; }

	      Jsc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for B_st */
	      Lsc = IMPOSSIBLE;
	      Rsc = IMPOSSIBLE;
	      Tsc = IMPOSSIBLE;

	      /* Careful with Tsc, it isn't updated for k == 0 or  k == d, 
	       * but Jsc, Lsc, Rsc, are all updated for k == 0 and k == d */
	      for (k = kmin; k <= kmax; k++) {
		Jsc = ESL_MAX(Jsc, (Jalpha_begl[jp_wA[k]][w][d-k] + Jalpha[jp_y][y][k]));
		Lsc = ESL_MAX(Lsc, (Jalpha_begl[jp_wA[k]][w][d-k] + Lalpha[jp_y][y][k]));
		Rsc = ESL_MAX(Rsc, (Ralpha_begl[jp_wA[k]][w][d-k] + Jalpha[jp_y][y][k]));
	      }		
	      kn = ESL_MAX(1,   kmin);
	      kx = ESL_MIN(d-1, kmax);
	      for (k = kn; k <= kx; k++) {
		Tsc = ESL_MAX(Tsc, (Ralpha_begl[jp_wA[k]][w][d-k] + Lalpha[jp_y][y][k]));
	      }

	      Jalpha[jp_v][v][d] = Jsc;
	      Talpha[jp_v][v][d] = Tsc;
	      if(kmin == 0) Lalpha[jp_v][v][d] = ESL_MAX(Lsc, ESL_MAX(Jalpha_begl[jp_wA[0]][w][d], Lalpha_begl[jp_wA[0]][w][d])); 
	      else          Lalpha[jp_v][v][d] = Lsc;

	      if(kmax == d) Ralpha[jp_v][v][d] = ESL_MAX(Rsc, ESL_MAX(Jalpha[jp_y][y][d], Ralpha[jp_y][y][d]));
	      else          Ralpha[jp_v][v][d] = Rsc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      Lsc = IMPOSSIBLE;
	      Rsc = IMPOSSIBLE;
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Jsc = ESL_MAX(Jsc, Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Lsc = ESL_MAX(Lsc, Lalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Rsc = ESL_MAX(Rsc, Ralpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha_begl[jp_v][v][d] = Jsc;
	      Lalpha_begl[jp_v][v][d] = Lsc;
	      Ralpha_begl[jp_v][v][d] = Rsc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else if (emitmode == EMITLEFT) {
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    assert(dnA[v] == 1);
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = IMPOSSIBLE;
	      Rsc = IMPOSSIBLE;
	      Ralpha[jp_v][v][d] = Rsc; /* this is important b/c if we're an IL, we'll access this cell in the recursion below for Ralpha */

	      /* We need to do separate 'for (yoffset...' loops for J
	       * and R matrices, because jp_v == jp_y for all states
	       * here, and for IL states, v can equal y+yoffset (when
	       * yoffset==0).  This means we have to fully calculate
	       * the Jalpha[jp_v][y+yoffset][d] cell (which is
	       * Jalpha[jp_v][v][d]) before we can start to calculate
	       * Ralpha[jp_v][v][d]. 
	       */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		Jsc = ESL_MAX(Jsc,         Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Lsc = ESL_MAX(Lsc,         Lalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha[jp_v][v][d] = Jsc + esc_v[dsq[i]];
	      Lalpha[jp_v][v][d] = (d >= 2) ? Lsc + esc_v[dsq[i]] : esc_v[dsq[i]];
	      
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		Rsc = ESL_MAX(Rsc, ESL_MAX(Jalpha[jp_y][y+yoffset][d]      + tsc_v[yoffset],
					   Ralpha[jp_y][y+yoffset][d]      + tsc_v[yoffset]));
	      }
	      Ralpha[jp_v][v][d] = Rsc;
	      i--;
	    }
	  }
	  else if (emitmode == EMITRIGHT) { 
	    y = cm->cfirst[v]; 
	    assert(dnA[v] == 1);
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = IMPOSSIBLE;
	      Rsc = IMPOSSIBLE;
	      Lalpha[jp_v][v][d] = Lsc; /* this is important b/c if we're an IR, we'll access this cell in the recursion below for Lalpha */
	      
	      /* We need to do separate 'for (yoffset...' loops for J
	       * and L matrices, because jp_v == jq_y for all states
	       * here, and for IR states, v can equal y+yoffset (when
	       * yoffset==0).  This means we have to fully calculate
	       * the Jalpha[jq_y][y+yoffset][d] cell (which is
	       * Jalpha[jp_v][v][d]) before we can start to calculate
	       * Lalpha[jp_v][v][d]. 
	       */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Jsc = ESL_MAX(Jsc,         Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Rsc = ESL_MAX(Rsc,         Ralpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha[jp_v][v][d] = Jsc + esc_j;
	      Ralpha[jp_v][v][d] = (d >= 2) ? Rsc + esc_j : esc_j;

	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Lsc = ESL_MAX(Lsc, ESL_MAX(Jalpha[jq_y][y+yoffset][d]     + tsc_v[yoffset],
					   Lalpha[jq_y][y+yoffset][d]     + tsc_v[yoffset]));
	      }
	      Lalpha[jp_v][v][d] = Lsc;
	    }
	  }
	  else if (emitmode == EMITPAIR) { 
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    assert(dnA[v] == 1);
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = IMPOSSIBLE;
	      Rsc = IMPOSSIBLE;
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Jsc = ESL_MAX(Jsc,         Jalpha[jp_y][y+yoffset][d - 2] + tsc_v[yoffset]);
		Lsc = ESL_MAX(Lsc, ESL_MAX(Jalpha[jq_y][y+yoffset][d - 1] + tsc_v[yoffset],
					   Lalpha[jq_y][y+yoffset][d - 1] + tsc_v[yoffset])),
		Rsc = ESL_MAX(Rsc, ESL_MAX(Jalpha[jp_y][y+yoffset][d - 1] + tsc_v[yoffset],
					   Ralpha[jp_y][y+yoffset][d - 1] + tsc_v[yoffset]));
	      }
	      Jalpha[jp_v][v][d] = (d >= 2) ? Jsc + esc_v[dsq[i]*cm->abc->Kp+dsq[j]] : IMPOSSIBLE;
	      Lalpha[jp_v][v][d] = (d >= 2) ? Lsc + lmesc_v[dsq[i]]                  : lmesc_v[dsq[i]];
	      Ralpha[jp_v][v][d] = (d >= 2) ? Rsc + rmesc_j                          : rmesc_j;
	      i--;
	    }
	  }
	  else { /* ! B_st && ! BEGL_S st && ! L_st && ! R_st && ! P_st (emitmode == EMITNONE) */
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = IMPOSSIBLE;
	      Rsc = IMPOSSIBLE;
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Jsc = ESL_MAX(Jsc, Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Lsc = ESL_MAX(Lsc, Lalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Rsc = ESL_MAX(Rsc, Ralpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha[jp_v][v][d] = Jsc;
	      Lalpha[jp_v][v][d] = Lsc;
	      Ralpha[jp_v][v][d] = Rsc;
	    }
	  }
	  if(vsc != NULL) {
	    if(cm->stid[v] == BIF_B) { 
	      for (d = dnA[v]; d <= dxA[v]; d++) { 
		vsc[v] = ESL_MAX(vsc[v], 
				 ESL_MAX(Jalpha[jp_v][v][d], 
					 ESL_MAX(Lalpha[jp_v][v][d], 
						 ESL_MAX(Ralpha[jp_v][v][d], Talpha[jp_v][v][d]))));
	      }
	    }
	    else if(cm->stid[v] == BEGL_S) { 
	      for (d = dnA[v]; d <= dxA[v]; d++) { 
		vsc[v] = ESL_MAX(vsc[v], 
				 ESL_MAX(Jalpha_begl[jp_v][v][d], 
					 ESL_MAX(Lalpha_begl[jp_v][v][d], Ralpha_begl[jp_v][v][d])));
	      }
	    }
	    else {
	      for (d = dnA[v]; d <= dxA[v]; d++) { 
		vsc[v] = ESL_MAX(vsc[v], 
				 ESL_MAX(Jalpha[jp_v][v][d], 
					 ESL_MAX(Lalpha[jp_v][v][d], Ralpha[jp_v][v][d])));
	      }
	    }
	  }
#if PRINTFALPHA
	  if(cm->stid[v] == BIF_B) { 
	    for(d = dnA[v]; d <= dxA[v]; d++) { 
	      printf("R j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		     j, v, d, 
		     NOT_IMPOSSIBLE(Jalpha[jp_v][v][d]) ? Jalpha[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Lalpha[jp_v][v][d]) ? Lalpha[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Ralpha[jp_v][v][d]) ? Ralpha[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Talpha[jp_v][v][d]) ? Talpha[jp_v][v][d] : -9999.9);
	    }
	  }
	  else if(cm->stid[v] == BEGL_S) { 
	    for(d = dnA[v]; d <= dxA[v]; d++) { 
	      printf("R j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		     j, v, d, 
		     NOT_IMPOSSIBLE(Jalpha_begl[jp_v][v][d]) ? Jalpha_begl[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Lalpha_begl[jp_v][v][d]) ? Lalpha_begl[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Ralpha_begl[jp_v][v][d]) ? Ralpha_begl[jp_v][v][d] : -9999.9, 
		     -9999.9);
	    }
	  }
	  else if(cm->stid[v] == BEGR_S) {
	    for(d = dnA[v]; d <= dxA[v]; d++) { 
	      printf("R j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		     j, v, d, 
		     NOT_IMPOSSIBLE(Jalpha[jp_v][v][d]) ? Jalpha[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Lalpha[jp_v][v][d]) ? Lalpha[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Ralpha[jp_v][v][d]) ? Ralpha[jp_v][v][d] : -9999.9,
		     -9999.9);
	    }
	  }
#endif 
#if 0 
	  else {
	    for(d = dnA[v]; d <= dxA[v]; d++) { 
	      printf("R j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		     j, v, d, 
		     NOT_IMPOSSIBLE(Jalpha[jp_v][v][d]) ? Jalpha[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Lalpha[jp_v][v][d]) ? Lalpha[jp_v][v][d] : -9999.9,
		     NOT_IMPOSSIBLE(Ralpha[jp_v][v][d]) ? Ralpha[jp_v][v][d] : -9999.9,
		     -9999.9);
	    }
	  }
	  printf("\n");
#endif
	}  /*loop over decks v>0 */
      
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
      for (d = dnA[0]; d <= dxA[0]; d++) {
	bestr[d]    = 0;	     /* root of the traceback = root state 0 */
	bestmode[d] = CM_HIT_MODE_J; /* mode of the parsetree, joint (J) */
	y = cm->cfirst[0];
	Jalpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, Jalpha[cur][y][d] + tsc_v[0]);
	Lalpha[jp_v][0][d] = IMPOSSIBLE;
	Ralpha[jp_v][0][d] = IMPOSSIBLE;
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) {
	  Jalpha[jp_v][0][d] = ESL_MAX (Jalpha[jp_v][0][d], (Jalpha[cur][y+yoffset][d] + tsc_v[yoffset]));
	  Lalpha[jp_v][0][d] = ESL_MAX (Lalpha[jp_v][0][d], (Lalpha[cur][y+yoffset][d] + tsc_v[yoffset]));
	  Ralpha[jp_v][0][d] = ESL_MAX (Ralpha[jp_v][0][d], (Ralpha[cur][y+yoffset][d] + tsc_v[yoffset]));
	}
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(NOT_IMPOSSIBLE(cm->beginsc[y])) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = j % (W+1);
	      /* the trCYK paper doesn't allow best_sc to be rooted at an S state, why not? */
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		if(Jalpha[jp_v][0][d] < (Jalpha_begl[jp_y][y][d] + cm->beginsc[y])) {
		  Jalpha[jp_v][0][d] = Jalpha_begl[jp_y][y][d] + cm->beginsc[y];
		  bestr[d] = y;
		}
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		if(Jalpha[jp_v][0][d] < (Jalpha[jp_y][y][d] + cm->beginsc[y])) {
		  Jalpha[jp_v][0][d] = Jalpha[jp_y][y][d] + cm->beginsc[y];
		  bestr[d] = y;
		}
	      }
	    }
	  }
	}
      }
      /* Finally, allow for truncated hits */
      for (y = 1; y < cm->M; y++) {
	jp_y = cur;
	/*if(NOT_IMPOSSIBLE(cm->beginsc[y])) {*/
	/*trunc_penalty = cm->beginsc[y];*/
	float trunc_penalty = 0.;
	if(cm->sttype[y] == B_st || cm->sttype[y] == MP_st || cm->sttype[y] == ML_st || cm->sttype[y] == MR_st) { 
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Jalpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Jalpha[jp_y][y][d] + trunc_penalty;
	      bestr[d] = y;
	      /* bestmode[d] remains CM_HIT_MODE_J */
	    }
	  }
	}
	if(cm->sttype[y] == B_st || cm->sttype[y] == MP_st || cm->sttype[y] == ML_st) { 
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Lalpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Lalpha[jp_y][y][d] + trunc_penalty;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_L;
	    }
	  }
	}
	if(cm->sttype[y] == B_st || cm->sttype[y] == MP_st || cm->sttype[y] == MR_st) { 
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Ralpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Ralpha[jp_y][y][d] + trunc_penalty;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_R;
	    }
	  }
	}
	if(cm->sttype[y] == B_st) { 
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Talpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Talpha[jp_y][y][d] + trunc_penalty;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_T;
	    }
	  }
	}
      }
      /* find the best score in J */
      for (d = dnA[0]; d <= dxA[0]; d++) {
	vsc_root = ESL_MAX(vsc_root, Jalpha[jp_v][0][d]);
      }
      /* update envi, envj, if nec */
      if(do_env_defn) { 
	for (d = dnA[0]; d <= dxA[0]; d++) {
	  if(Jalpha[jp_v][0][d] >= env_cutoff) { 
	    envi = ESL_MIN(envi, j-d+1);
	    envj = ESL_MAX(envj, j);
	  }
	}
      }

      /* update gamma, but only if we're reporting hits to th */
      if(hitlist != NULL) if((status = UpdateGammaHitMx(cm, errbuf, gamma, jp_g, Jalpha[jp_v][0], dnA[0], dxA[0], FALSE, trsmx->bestr, trsmx->bestmode, hitlist, W, act)) != eslOK) return status;

      /* cm_DumpScanMatrixAlpha(cm, si, j, i0, TRUE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* find the best score in any matrix at any state */
  float best_tr_sc = IMPOSSIBLE;
  for(v = 0; v < cm->M; v++) { 
    best_tr_sc = ESL_MAX(best_tr_sc, vsc[v]);
  }
  printf("Best truncated score: %.4f (%.4f)\n",
	 best_tr_sc, 
	 best_tr_sc + sreLOG2(2./(cm->clen * (cm->clen+1))));

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in UpdateGammaHitMx() for each position j */
  if(hitlist != NULL && gamma->iamgreedy == FALSE) 
    TBackGammaHitMx(gamma, hitlist, i0, j0);

  /* set envelope return variables if nec */
  if(ret_envi != NULL) { *ret_envi = (envi == j0+1) ? -1 : envi; }
  if(ret_envj != NULL) { *ret_envj = (envj == i0-1) ? -1 : envj; }

  /* clean up and return */
  if(gamma != NULL) FreeGammaHitMx(gamma);
  if (act != NULL) { 
    for(i = 0; i <= W; i++) free(act[i]); 
    free(act);
  }
  free(jp_wA);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  else free(vsc);
  if (ret_sc != NULL) *ret_sc = vsc_root;

  ESL_DPRINTF1(("RefTrCYKScan() return score: %10.4f\n", vsc_root)); 
  return eslOK;
  
 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

/* Function: RefITrInsideScan()
 * Date:     EPN, Wed Aug 24 15:23:37 2011
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           a reference trInside scanning algorithm. Query-dependent 
 *           bands are used or not used as specified in ScanMatrix_t <si>.
 *
 * Args:     cm              - the covariance model
 *           errbuf          - char buffer for reporting errors
 *           smx             - ScanMatrix_t for this search w/this model (incl. DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           cutoff          - minimum score to report
 *           th              - CM_TOPHITS to add to; if NULL, don't add to it
 *           do_null3        - TRUE to do NULL3 score correction, FALSE not to
 *           env_cutoff      - ret_envi..ret_envj will include all hits that exceed this bit sc
 *           ret_envi        - min position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted
 *           ret_envj        - max position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted 
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *           ret_sc          - RETURN: score of best overall hit (vsc[0])
 *
 * Note:     This function is heavily synchronized with RefTrCYKScan()
 *           any change to this function should be mirrored in that functions. 
 *
 * Returns:  eslOK on succes;
 *           <ret_sc> is score of best overall hit (vsc[0]). Information on hits added to <hitlist>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
int
RefITrInsideScan(CM_t *cm, char *errbuf, TrScanMatrix_t *trsmx, ESL_DSQ *dsq, int i0, int j0, float cutoff, CM_TOPHITS *hitlist,
		 int do_null3, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float **ret_vsc, float *ret_sc)
{
  int       status;
  GammaHitMx_t *gamma;          /* semi-HMM for hit resoultion */
  int      ivsc;                /* integer score */
  float    *vsc;                /* best score for each state (float) */
  float     vsc_root;           /* best overall score (score at ROOT_S) */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jq_y;  	        /* offset j for state y plus 1 (if jp_y is prv, jq_y is cur, and vice versa) */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       W;                  /* max d; max size of a hit, this is min(L, smx->W) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       kn, kx;             /* minimum/maximum valid k for current d in B_st recursion */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  int     **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */
  float    *gamma_row;          /* holds floatized scores for updating gamma matrix, only really used if hitlist != NULL */
  double  **act;                /* [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1) */
  int       do_env_defn;        /* TRUE to calculate envi, envj, FALSE not to (TRUE if ret_envi != NULL or ret_envj != NULL */
  int64_t   envi, envj;         /* min/max positions that exist in any hit with sc >= env_cutoff */

  /* Contract check */
  if(! cm->flags & CMH_BITS)               ESL_FAIL(eslEINCOMPAT, errbuf, "RefITrInsideScan, CMH_BITS flag is not raised.\n");
  if(j0 < i0)                              ESL_FAIL(eslEINCOMPAT, errbuf, "RefITrInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                          ESL_FAIL(eslEINCOMPAT, errbuf, "RefITrInsideScan, dsq is NULL\n");
  if(trsmx == NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "RefITrInsideScan, trsmx == NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE)   ESL_FAIL(eslEINCOMPAT, errbuf, "RefITrInsideScan, CM_SEARCH_INSIDE flag raised");
  if(! (trsmx->flags & cmTRSMX_HAS_INT))   ESL_FAIL(eslEINCOMPAT, errbuf, "RefITrInsideScan, ScanMatrix's cmTRSMX_HAS_FLOAT flag is not raised");

  /* make pointers to the ScanMatrix/CM data for convenience */
  int  ***Jalpha       = trsmx->iJalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Jalpha DP matrix, NULL for v == BEGL_S */
  int  ***Jalpha_begl  = trsmx->iJalpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] Jalpha DP matrix, NULL for v != BEGL_S */
  int  ***Lalpha       = trsmx->iLalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Lalpha DP matrix, NULL for v == BEGL_S */
  int  ***Lalpha_begl  = trsmx->iLalpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] Lalpha DP matrix, NULL for v != BEGL_S */
  int  ***Ralpha       = trsmx->iRalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Ralpha DP matrix, NULL for v == BEGL_S */
  int  ***Ralpha_begl  = trsmx->iRalpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] Ralpha DP matrix, NULL for v != BEGL_S */
  int  ***Talpha       = trsmx->iTalpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] Talpha DP matrix, NULL for v != BIF_B  */
  int   **dnAA         = trsmx->dnAA;         /* [0..j..W][0..v..cm->M-1] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA         = trsmx->dxAA;         /* [0..j..W][0..v..cm->M-1] maximum d for v, j (for j > W use [v][W]) */
  int    *bestr        = trsmx->bestr;        /* [0..d..W] best root state (for local begins or 0) for this d */
  int    *bestmode     = trsmx->bestmode;     /* [0..d..W] mode of best parsetree for this d */
  int    *dmax         = trsmx->dmax;         /* [0..v..cm->M-1] maximum d allowed for this state */
  int   **esc_vAA      = cm->ioesc;           /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
 					       * and all possible emissions a (including ambiguities) */
  int   **lmesc_vAA    = cm->ilmesc;          /* [0..v..cm->M-1][0..a..(cm->abc->Kp-1)] left  marginal emission scores for v */
  int   **rmesc_vAA    = cm->irmesc;          /* [0..v..cm->M-1][0..a..(cm->abc->Kp-1)] right marginal emission scores for v */

  /* determine if we're doing banded/non-banded */
  if(trsmx->dmax != NULL) do_banded = TRUE;
  
  L = j0-i0+1;
  W = trsmx->W;
  if (W > L) W = L; 

  /* set vsc array */
  vsc = NULL;
  ESL_ALLOC(vsc, sizeof(float) * cm->M);
  esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);
  vsc_root = IMPOSSIBLE;

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. */
  if(hitlist != NULL) gamma = CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff, FALSE);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* precalculate the initial scores for all cells */
  init_scAA = ICalcInitDPScores(cm);

  /* allocate/initialize gamma_row, only used for updating gamma if hitlist != NULL */
  ESL_ALLOC(gamma_row, sizeof(float) * (W+1));
  esl_vec_FSet(gamma_row, (W+1), IMPOSSIBLE);

  /* if do_null3: allocate and initialize act vector */
  if(do_null3) { 
    ESL_ALLOC(act, sizeof(double *) * (W+1));
    for(i = 0; i <= W; i++) { 
      ESL_ALLOC(act[i], sizeof(double) * cm->abc->K);
      esl_vec_DSet(act[i], cm->abc->K, 0.);
    }
  }
  else act = NULL;

  /* initialize envelope boundary variables */
  do_env_defn = (ret_envi != NULL || ret_envj != NULL) ? TRUE : FALSE;
  envi = j0+1;
  envj = i0-1;

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float Jsc, Lsc, Rsc, Tsc;
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

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  int const *esc_v   = esc_vAA[v]; 
	  int const *tsc_v   = cm->itsc[v];
	  int const *lmesc_v = lmesc_vAA[v]; 
	  int const *rmesc_v = rmesc_vAA[v]; 
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v  = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y  = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  jq_y = (StateRightDelta(cm->sttype[v]) > 0) ? cur : prv;
	  sd    = StateDelta(cm->sttype[v]);
	  cnum  = cm->cnum[v];
	  /* if we emit right, precalc score of emitting res j from state v */
	  int   esc_j = -INFTY;
	  int rmesc_j = -INFTY;
	  if(cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) { 
	    esc_j   =   esc_v[dsq[j]];
	    rmesc_j = rmesc_v[dsq[j]];
	  }
	  if(cm->sttype[v] == MP_st) { 
	    rmesc_j = rmesc_v[dsq[j]];
	  }

	  if(cm->sttype[v] == B_st) {
	    w = cm->cfirst[v]; /* BEGL_S */
	    y = cm->cnum[v];   /* BEGR_S */
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(0, (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], d);
	      }
	      else { kmin = 0; kmax = d; }

	      Jsc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for B_st */
	      Lsc = -INFTY;
	      Rsc = -INFTY;
	      Tsc = -INFTY;

	      /* Careful with Tsc, it isn't updated for k == 0 or  k == d, 
	       * but Jsc, Lsc, Rsc, are all updated for k == 0 and k == d */
	      for (k = kmin; k <= kmax; k++) {
		//Jsc = ESL_MAX(Jsc, (Jalpha_begl[jp_wA[k]][w][d-k] + Jalpha[jp_y][y][k]));
		//Lsc = ESL_MAX(Lsc, (Jalpha_begl[jp_wA[k]][w][d-k] + Lalpha[jp_y][y][k]));
		//Rsc = ESL_MAX(Rsc, (Ralpha_begl[jp_wA[k]][w][d-k] + Jalpha[jp_y][y][k]));
		Jsc = ILogsum(Jsc, (Jalpha_begl[jp_wA[k]][w][d-k] + Jalpha[jp_y][y][k]));
		Lsc = ILogsum(Lsc, (Jalpha_begl[jp_wA[k]][w][d-k] + Lalpha[jp_y][y][k]));
		Rsc = ILogsum(Rsc, (Ralpha_begl[jp_wA[k]][w][d-k] + Jalpha[jp_y][y][k]));
	      }		
	      kn = ESL_MAX(1,   kmin);
	      kx = ESL_MIN(d-1, kmax);
	      for (k = kn; k <= kx; k++) {
		Tsc = ILogsum(Tsc, (Ralpha_begl[jp_wA[k]][w][d-k] + Lalpha[jp_y][y][k]));
	      }

	      Jalpha[jp_v][v][d] = Jsc;
	      Talpha[jp_v][v][d] = Tsc;
	      if(kmin == 0) Lalpha[jp_v][v][d] = ILogsum(Lsc, ESL_MAX(Jalpha_begl[jp_wA[0]][w][d], Lalpha_begl[jp_wA[0]][w][d])); 
	      else          Lalpha[jp_v][v][d] = Lsc;

	      if(kmax == d) Ralpha[jp_v][v][d] = ILogsum(Rsc, ESL_MAX(Jalpha[jp_y][y][d], Ralpha[jp_y][y][d]));
	      else          Ralpha[jp_v][v][d] = Rsc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      Lsc = -INFTY;
	      Rsc = -INFTY;
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		//Jsc = ESL_MAX(Jsc, Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		//Lsc = ESL_MAX(Lsc, Lalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		//Rsc = ESL_MAX(Rsc, Ralpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Jsc = ILogsum(Jsc, Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Lsc = ILogsum(Lsc, Lalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Rsc = ILogsum(Rsc, Ralpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha_begl[jp_v][v][d] = Jsc;
	      Lalpha_begl[jp_v][v][d] = Lsc;
	      Ralpha_begl[jp_v][v][d] = Rsc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else if (emitmode == EMITLEFT) {
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    assert(dnA[v] == 1);
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = -INFTY;
	      Rsc = -INFTY;
	      Ralpha[jp_v][v][d] = Rsc; /* this is important b/c if we're an IL, we'll access this cell in the recursion below for Ralpha */

	      /* We need to do separate 'for (yoffset...' loops for J
	       * and R matrices, because jp_v == jp_y for all states
	       * here, and for IL states, v can equal y+yoffset (when
	       * yoffset==0).  This means we have to fully calculate
	       * the Jalpha[jp_v][y+yoffset][d] cell (which is
	       * Jalpha[jp_v][v][d]) before we can start to calculate
	       * Ralpha[jp_v][v][d]. 
	       */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		Jsc = ILogsum(Jsc,         Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Lsc = ILogsum(Lsc,         Lalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha[jp_v][v][d] = Jsc + esc_v[dsq[i]];
	      Lalpha[jp_v][v][d] = (d >= 2) ? Lsc + esc_v[dsq[i]] : esc_v[dsq[i]];
	      
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		Rsc = ILogsum(Rsc, ILogsum(Jalpha[jp_y][y+yoffset][d]      + tsc_v[yoffset],
					   Ralpha[jp_y][y+yoffset][d]      + tsc_v[yoffset]));
	      }
	      Ralpha[jp_v][v][d] = Rsc;
	      i--;
	    }
	  }
	  else if (emitmode == EMITRIGHT) { 
	    y = cm->cfirst[v]; 
	    assert(dnA[v] == 1);
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = -INFTY;
	      Rsc = -INFTY;
	      Lalpha[jp_v][v][d] = Lsc; /* this is important b/c if we're an IR, we'll access this cell in the recursion below for Lalpha */
	      
	      /* We need to do separate 'for (yoffset...' loops for J
	       * and L matrices, because jp_v == jq_y for all states
	       * here, and for IR states, v can equal y+yoffset (when
	       * yoffset==0).  This means we have to fully calculate
	       * the Jalpha[jq_y][y+yoffset][d] cell (which is
	       * Jalpha[jp_v][v][d]) before we can start to calculate
	       * Lalpha[jp_v][v][d]. 
	       */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Jsc = ILogsum(Jsc,         Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Rsc = ILogsum(Rsc,         Ralpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha[jp_v][v][d] = Jsc + esc_j;
	      Ralpha[jp_v][v][d] = (d >= 2) ? Rsc + esc_j : esc_j;

	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Lsc = ILogsum(Lsc, ILogsum(Jalpha[jq_y][y+yoffset][d]     + tsc_v[yoffset],
					   Lalpha[jq_y][y+yoffset][d]     + tsc_v[yoffset]));
	      }
	      Lalpha[jp_v][v][d] = Lsc;
	    }
	  }
	  else if (emitmode == EMITPAIR) { 
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    assert(dnA[v] == 1);
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = -INFTY;
	      Rsc = -INFTY;
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Jsc = ILogsum(Jsc,         Jalpha[jp_y][y+yoffset][d - 2] + tsc_v[yoffset]);
		Lsc = ILogsum(Lsc, ESL_MAX(Jalpha[jq_y][y+yoffset][d - 1] + tsc_v[yoffset],
					   Lalpha[jq_y][y+yoffset][d - 1] + tsc_v[yoffset])),
		Rsc = ILogsum(Rsc, ESL_MAX(Jalpha[jp_y][y+yoffset][d - 1] + tsc_v[yoffset],
					   Ralpha[jp_y][y+yoffset][d - 1] + tsc_v[yoffset]));
	      }
	      Jalpha[jp_v][v][d] = (d >= 2) ? Jsc + esc_v[dsq[i]*cm->abc->Kp+dsq[j]] : -INFTY;
	      Lalpha[jp_v][v][d] = (d >= 2) ? Lsc + lmesc_v[dsq[i]]                  : lmesc_v[dsq[i]];
	      Ralpha[jp_v][v][d] = (d >= 2) ? Rsc + rmesc_j                          : rmesc_j;
	      i--;
	    }
	  }
	  else { /* ! B_st && ! BEGL_S st && ! L_st && ! R_st && ! P_st (emitmode == EMITNONE) */
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      Jsc = init_scAA[v][d-sd]; 
	      Lsc = -INFTY;
	      Rsc = -INFTY;
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) { 
		Jsc = ILogsum(Jsc, Jalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Lsc = ILogsum(Lsc, Lalpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		Rsc = ILogsum(Rsc, Ralpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      }
	      Jalpha[jp_v][v][d] = Jsc;
	      Lalpha[jp_v][v][d] = Lsc;
	      Ralpha[jp_v][v][d] = Rsc;
	    }
	  }
	  if(vsc != NULL) {
	    ivsc = -INFTY;
	    if(cm->stid[v] == BIF_B) { 
	      for (d = dnA[v]; d <= dxA[v]; d++) { 
		ivsc = ESL_MAX(ivsc,
			       ESL_MAX(Jalpha[jp_v][v][d], 
				       ESL_MAX(Lalpha[jp_v][v][d],
					       ESL_MAX(Ralpha[jp_v][v][d], Talpha[jp_v][v][d]))));
	      }
	    }
	    else if(cm->stid[v] == BEGL_S) { 
	      for (d = dnA[v]; d <= dxA[v]; d++) { 
		ivsc = ESL_MAX(ivsc, 
			       ESL_MAX(Jalpha_begl[jp_v][v][d], 
				       ESL_MAX(Lalpha_begl[jp_v][v][d], Ralpha_begl[jp_v][v][d])));
	      }
	    }
	    else {
	      for (d = dnA[v]; d <= dxA[v]; d++) { 
		ivsc = ESL_MAX(ivsc, 
			       ESL_MAX(Jalpha[jp_v][v][d], 
				       ESL_MAX(Lalpha[jp_v][v][d], Ralpha[jp_v][v][d])));
	      }
	    }
	    vsc[v] = Scorify(ivsc);
	  }
#if PRINTIALPHA
	  if(cm->stid[v] == BIF_B) { 
	    for(d = dnA[v]; d <= dxA[v]; d++) { 
	      printf("R j: %3d  v: %3d  d: %3d   J: %10d  L: %10d  R: %10d  T: %10d\n", 
		     j, v, d, 
		     Jalpha[jp_v][v][d], Lalpha[jp_v][v][d], Ralpha[jp_v][v][d],
		     Talpha[jp_v][v][d]);
	    }
	  }
	  else if(cm->stid[v] == BEGL_S) { 
	    for(d = dnA[v]; d <= dxA[v]; d++) { 
	      printf("R j: %3d  v: %3d  d: %3d   J: %10d  L: %10d  R: %10d  T: %10d\n", 
		     j, v, d, 
		     Jalpha_begl[jp_v][v][d], Lalpha_begl[jp_v][v][d], Ralpha_begl[jp_v][v][d], -INFTY);
	    }
	  }
	  else {
	    for(d = dnA[v]; d <= dxA[v]; d++) { 
	      printf("R j: %3d  v: %3d  d: %3d   J: %10d  L: %10d  R: %10d\n",
		     j, v, d, 
		     Jalpha[jp_v][v][d], Lalpha[jp_v][v][d], Ralpha[jp_v][v][d]);
	    }
	  }
	  printf("\n");
#endif
	}  /*loop over decks v>0 */
      
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
      for (d = dnA[0]; d <= dxA[0]; d++) {
	bestr[d]    = 0;	     /* root of the traceback = root state 0 */
	bestmode[d] = CM_HIT_MODE_J; /* mode of the parsetree, joint (J) */
	y = cm->cfirst[0];
	Jalpha[jp_v][0][d] = ESL_MAX(-INFTY, Jalpha[cur][y][d] + tsc_v[0]);
	Lalpha[jp_v][0][d] = -INFTY;
	Ralpha[jp_v][0][d] = -INFTY;
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) {
	  Jalpha[jp_v][0][d] = ILogsum (Jalpha[jp_v][0][d], (Jalpha[cur][y+yoffset][d] + tsc_v[yoffset]));
	  Lalpha[jp_v][0][d] = ILogsum (Lalpha[jp_v][0][d], (Lalpha[cur][y+yoffset][d] + tsc_v[yoffset]));
	  Ralpha[jp_v][0][d] = ILogsum (Ralpha[jp_v][0][d], (Ralpha[cur][y+yoffset][d] + tsc_v[yoffset]));
	}
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = j % (W+1);
	      /* the trCYK paper doesn't allow best_sc to be rooted at an S state, why not? */
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		if(Jalpha[jp_v][0][d] < (Jalpha_begl[jp_y][y][d] + cm->ibeginsc[y])) {
		  Jalpha[jp_v][0][d] = Jalpha_begl[jp_y][y][d] + cm->ibeginsc[y];
		  bestr[d] = y;
		}
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		if(Jalpha[jp_v][0][d] < (Jalpha[jp_y][y][d] + cm->ibeginsc[y])) {
		  Jalpha[jp_v][0][d] = Jalpha[jp_y][y][d] + cm->ibeginsc[y];
		  bestr[d] = y;
		}
	      }
	    }
	  }
	}
      }
      /* Finally, allow for truncated hits */
      for (y = 1; y < cm->M; y++) {
	jp_y = cur;
	/*if(cm->ibeginsc[y] != -INFTY) {*/
	/*trunc_penalty = cm->ibeginsc[y];*/
	int trunc_penalty = 0;
	if(cm->sttype[y] == B_st || cm->sttype[y] == MP_st || cm->sttype[y] == ML_st || cm->sttype[y] == MR_st) { 
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Jalpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Jalpha[jp_y][y][d] + trunc_penalty;
	      bestr[d] = y;
	      /* bestmode[d] remains CM_HIT_MODE_J */
	    }
	  }
	}
	if(cm->sttype[y] == B_st || cm->sttype[y] == MP_st || cm->sttype[y] == ML_st) { 
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Lalpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Lalpha[jp_y][y][d] + trunc_penalty;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_L;
	    }
	  }
	}
	if(cm->sttype[y] == B_st || cm->sttype[y] == MP_st || cm->sttype[y] == MR_st) { 
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Ralpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Ralpha[jp_y][y][d] + trunc_penalty;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_R;
	    }
	  }
	}
	if(cm->sttype[y] == B_st) { 
	  /*if(y == 54) { printf("Talpha[%d][%d][%d] %.4f\n", jp_y, y, d, Talpha[jp_y][y][d]); }*/
	  for (d = dnA[y]; d <= dxA[y]; d++) {
	    if(Jalpha[jp_v][0][d] < (Talpha[jp_y][y][d] + trunc_penalty)) {
	      Jalpha[jp_v][0][d] = Talpha[jp_y][y][d] + trunc_penalty;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_T;
	    }
	  }
	}
      }
      /* find the best score in J */
      for (d = dnA[0]; d <= dxA[0]; d++) {
	vsc_root = ESL_MAX(vsc_root, Scorify(Jalpha[jp_v][0][d]));
      }
      /* update envi, envj, if nec */
      if(do_env_defn) { 
	for (d = dnA[0]; d <= dxA[0]; d++) {
	  if(Scorify(Jalpha[jp_v][0][d]) >= env_cutoff) { 
	    envi = ESL_MIN(envi, j-d+1);
	    envj = ESL_MAX(envj, j);
	  }
	}
      }

      /* update gamma, but only if we're reporting hits to th */
      if(hitlist != NULL) { 
	for(d = dnA[0]; d <= dxA[0]; d++) { gamma_row[d] = Scorify(Jalpha[jp_v][0][d]); }
	if((status = UpdateGammaHitMx(cm, errbuf, gamma, jp_g, gamma_row, dnA[0], dxA[0], FALSE, trsmx->bestr, trsmx->bestmode, hitlist, W, act)) != eslOK) return status;
      }

      /* cm_DumpScanMatrixAlpha(cm, si, j, i0, TRUE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* find the best score in any matrix at any state */
  float best_tr_sc = IMPOSSIBLE;
  for(v = 0; v < cm->M; v++) { 
    best_tr_sc = ESL_MAX(best_tr_sc, vsc[v]);
  }
  printf("Best truncated score: %.4f (%.4f)\n",
	 best_tr_sc, 
	 best_tr_sc + sreLOG2(2./(cm->clen * (cm->clen+1))));

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in UpdateGammaHitMx() for each position j */
  if(hitlist != NULL && gamma->iamgreedy == FALSE) 
    TBackGammaHitMx(gamma, hitlist, i0, j0);

  /* set envelope return variables if nec */
  if(ret_envi != NULL) { *ret_envi = (envi == j0+1) ? -1 : envi; }
  if(ret_envj != NULL) { *ret_envj = (envj == i0-1) ? -1 : envj; }

  /* clean up and return */
  if(gamma != NULL) FreeGammaHitMx(gamma);
  if (act != NULL) { 
    for(i = 0; i <= W; i++) free(act[i]); 
    free(act);
  }
  free(jp_wA);
  free(init_scAA[0]);
  free(init_scAA);
  free(gamma_row);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  else free(vsc);
  if (ret_sc != NULL) *ret_sc = vsc_root;

  ESL_DPRINTF1(("RefTrCYKScan() return score: %10.4f\n", vsc_root)); 
  return eslOK;
  
 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/* Function: FastTrCYKScanHB()
 * Incept:   EPN, Thu Aug 25 15:19:28 2011
 *
 * Purpose:  An HMM banded version of a scanning TrCYK algorithm. Takes
 *           a CM_TR_HB_MX data structure which is indexed [v][j][d] with
 *           only cells within the bands allocated.
 *           (different than other (non-HB) scanning function's convention 
 *            of [j][v][d]).
 *
 * Args:     cm        - the model    [0..M-1]
 *           sq        - the sequence [1..L]   
 *                     - length of the dsq
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           cutoff    - minimum score to report
 *           hitlist   - CM_TOPHITS hitlist to add to; if NULL, don't add to it
 *           do_null3  - TRUE to do NULL3 score correction, FALSE not to
 *           env_cutoff- ret_envi..ret_envj will include all hits that exceed this bit sc
 *           ret_envi  - min position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted
 *           ret_envj  - max position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted 
 *           mx        - the dp matrix, only cells within bands in cm->cp9b will 
 *                       be valid. 
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           ret_sc    - RETURN: score of best overall hit (vsc[0])
 *                       
 * Returns: eslOK on success; eslERANGE if required matrix size exceeds <size_limit>
 *          <ret_sc>: score of the best hit.
 */
int
FastTrCYKScanHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, 
		CM_TR_HB_MX *mx, float size_limit, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float *ret_sc)
{
  int      status;
  GammaHitMx_t *gamma;  /* semi-HMM for hit resoultion */
  int     *bestr;       /* best root state for d at current j */
  int     *bestmode;    /* best mode for parsetree for d at current j */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc, Jsc, Lsc, Rsc, Tsc; /* temporary scores */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      j_sdr;              /* j - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y, dp_z;   /* offset d index for states v, y, z */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp;                 /* ESL_MAX(d-sd, 0) */
  int      dp_y_sd;            /* dp_y - sd */
  int      dp_y_sdr;           /* dp_y - sdr, often for jp_y_sdr */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  float    tsc;                /* a transition score */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  float    vsc_root;           /* score of best hit */
  int      W;                  /* max d over all hdmax[v][j] for all valid v, j */
  double **act;                /* [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1) */
  int      jp;                 /* j index in act */
  int      do_env_defn;        /* TRUE to calculate envi, envj, FALSE not to (TRUE if ret_envi != NULL or ret_envj != NULL */
  int64_t  envi, envj;         /* min/max positions that exist in any hit with sc >= env_cutoff */
  float    trunc_penalty = 0.; /* penalty in bits for a truncated hit */

  /* Contract check */
  if(dsq == NULL)       ESL_FAIL(eslEINCOMPAT, errbuf, "FastTrCYKScanHB(), dsq is NULL.\n");
  if (mx == NULL)       ESL_FAIL(eslEINCOMPAT, errbuf, "FastTrCYKScanHB(), mx is NULL.\n");
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "FastTrCYKScanHB(), mx is NULL.\n");

  ESL_DPRINTF1(("cm->search_opts & CM_SEARCH_HMMALNBANDS: %d\n", cm->search_opts & CM_SEARCH_HMMALNBANDS));

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b   = cm->cp9b; 
  int        *jmin   = cp9b->jmin;  
  int        *jmax   = cp9b->jmax;
  int       **hdmin  = cp9b->hdmin;
  int       **hdmax  = cp9b->hdmax;
  /* the DP matrix */
  float ***Jalpha  = mx->Jdp; /* pointer to the Jalpha DP matrix */
  float ***Lalpha  = mx->Ldp; /* pointer to the Lalpha DP matrix */
  float ***Ralpha  = mx->Rdp; /* pointer to the Ralpha DP matrix */
  float ***Talpha  = mx->Tdp; /* pointer to the Talpha DP matrix */

  /* Allocations and initializations  */
  /* grow the matrix based on the current sequence and bands */
  if((status = cm_tr_hb_mx_GrowTo(cm, mx, errbuf, cp9b, (j0-i0+1), size_limit)) != eslOK) return status;

  /* set W as j0-i0+1 (this may exceed max size of a hit our bands will allow, 
   * but that's okay b/c W is only used for sizing of act and bestr vectors */
  W = j0-i0+1;
  /* make sure our bands won't allow a hit bigger than W (this could be modified to only execute in debugging mode) */
  for(j = jmin[0]; j <= jmax[0]; j++) {
    if(W < (hdmax[0][(j-jmin[0])])) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "FastCYKScanHB(), band allows a hit (j:%d hdmax[0][j]:%d) greater than j0-i0+1 (%d)", j, hdmax[0][(j-jmin[0])], j0-i0+1);
  }

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  ESL_STOPWATCH *w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(Jalpha[0][0], mx->JLRncells_valid, IMPOSSIBLE);
  esl_vec_FSet(Lalpha[0][0], mx->JLRncells_valid, IMPOSSIBLE);
  esl_vec_FSet(Ralpha[0][0], mx->JLRncells_valid, IMPOSSIBLE);
  if(mx->Tncells_valid > 0)   esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid,   IMPOSSIBLE); /* careful use Tdp_mem not Talpha[0][0] because Talpha[0] is NULL */
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " Matrix init CPU time: ");

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse of multiple nonoverlapping hits. */
  if(hitlist != NULL) gamma = CreateGammaHitMx(j0-i0+1, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff, FALSE);
  else                gamma = NULL;

  /* if do_null3: allocate and initialize act vector */
  if(do_null3) { 
    ESL_ALLOC(act, sizeof(double *) * (W+1));
    for(i = 0; i <= W; i++) { 
      ESL_ALLOC(act[i], sizeof(double) * cm->abc->K);
      esl_vec_DSet(act[i], cm->abc->K, 0.);
    }
    /* pre-fill act, different than non-HMM banded scanner b/c our main loop doesn't step j through residues */
    for(j = i0; j <= j0; j++) { 
      jp = j-i0+1; /* j is actual index in dsq, jp_g is offset j relative to start i0 (j index for act) */
      esl_vec_DCopy(act[(jp-1)%(W+1)], cm->abc->K, act[jp%(W+1)]);
      esl_abc_DCount(cm->abc, act[jp%(W+1)], dsq[j], 1.);
    }
  }
  else act = NULL;

  /* initialize envelope boundary variables */
  do_env_defn = (ret_envi != NULL || ret_envj != NULL) ? TRUE : FALSE;
  envi = j0+1;
  envj = i0-1;

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) { /* all the way down to root, different from other scanners */
    float const *esc_v   = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v   = cm->tsc[v];  /* transition scores for state v */
    float const *lmesc_v = cm->lmesc[v];
    float const *rmesc_v = cm->rmesc[v];
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];

    /* re-initialize the deck if we can do a local end from v */
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	  dp = ESL_MAX(d-sd, 0);
	  Jalpha[v][jp_v][dp_v] = el_scA[dp] + cm->endsc[v];
	  /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	Jalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	Lalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	Ralpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IL states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1] 
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr valid for state y? */
	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	  i = j - d + 1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */

	  /* We need to treat R differently from and J and L here, by
	   * doing separate 'for (yoffset...' loops for J because we
	   * have to fully calculate Jalpha[v][jp_v][dp_v]) before we
	   * can start to calculate Ralpha[v][jp_v][dp_v].
	   */
	  /* Handle J and L first */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	      Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]];
	  Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  Lalpha[v][jp_v][dp_v] = (d >= 2) ? Lalpha[v][jp_v][dp_v] + esc_v[dsq[i]]: esc_v[dsq[i]];
	  Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  i--;

	  /* Handle R separately */
	  Rsc = Ralpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    /* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,L above) */
	    if((d) >= hdmin[y][jp_y_sdr] && (d) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
	      dp_y = d - hdmin[y][jp_y_sdr];
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	      Rsc = ESL_MAX(Rsc, ESL_MAX(Jalpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset],
					 Ralpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]));
	    }
	  } /* end of for (yvalid_idx = 0... loop */
	  Ralpha[v][jp_v][dp_v] = Rsc; 
	  /* we use Rsc instead of Ralpha cell in above loop because
	   * Ralpha[v][jp_v][dp_v] may be the same cell as
	   * Ralpha[y][jp_y_sdr][dp_y] if we're an IL state
	   */
	}
      }
    }
    else if(cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IR states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1].
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */

      /* The first MR_st/IR_st 'for (j...' loop is for J and R matrices which use the same set of j values */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */

	  /* We need to treat L differently from and J and R here, by
	   * doing separate 'for (yoffset...' loops for J because we
	   * have to fully calculate Jalpha[v][jp_v][dp_v]) before we
	   * can start to calculate Lalpha[v][jp_v][dp_v].
	   */
	  /* Handle J and R first */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	      Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	    }
	  }
	  Jalpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	  Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  Ralpha[v][jp_v][dp_v] = (d >= 2) ? Ralpha[v][jp_v][dp_v] + esc_v[dsq[j]] : esc_v[dsq[j]];
	  Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
      
      /* The second MR_st/IR_st 'for (j...' loop is for the L matrix which use a different set of j values */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	
	/* determine which children y we can legally transit to for v, j */
	/* we use 'j' and not 'j_sdr' here for the L matrix, differently from J and R matrices above */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j) >= jmin[y] && ((j) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j is valid for state y? */
	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	  
	  Lsc = Lalpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    /* we use 'jp_y=j-min[y]' here, not 'jp_y_sdr=j-jmin[y]-sdr' (which we used in the corresponding loop for J,R above) */
	    jp_y = j - jmin[y];
	    
	    /* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,R above) */
	    if((d) >= hdmin[y][jp_y] && (d) <= hdmax[y][jp_y]) { /* make sure d is valid for this v, j and y */
	      dp_y = d - hdmin[y][jp_y];
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v] - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y] - hdmin[y][jp_y])));
	      Lsc = ESL_MAX(Lsc, ESL_MAX(Jalpha[y][jp_y][dp_y] + tsc_v[yoffset], 
					 Lalpha[y][jp_y][dp_y] + tsc_v[yoffset]));
	    }
	  } /* end of for (yvalid_idx = 0... loop */
	  Lalpha[v][jp_v][dp_v] = Lsc; 
	  /* we use Lsc instead of Lalpha cell in above loop because
	   * Lalpha[v][jp_v][dp_v] may be the same cell as
	   * Lalpha[y][jp_y_sdr][dp_y] if we're an IR state
	   */
	}
      }
    }
    else if(cm->sttype[v] == MP_st) { 
      /* MP states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	/* The first MP_st 'for (jp_v...' loop is for J and R matrices which use the same set of j values */
	/* j must satisfy:
	 * j >= jmin[v]
	 * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	 * j <= jmax[v]
	 * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	 * this reduces to two ESL_MAX calls
	 */
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	/* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	
	for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++, jp_y++) {
	  ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	  ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));

	  /* J matrix: */
	  /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	   * this reduces to two ESL_MAX calls
	   */
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	  dpn       = dn - hdmin[v][jp_v];
	  dpx       = dx - hdmin[v][jp_v];
	  dp_y_sd   = dn - hdmin[y][jp_y_sdr] - sd;

	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v      >= 0 && dp_v       <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	    ESL_DASSERT1((dp_y_sd   >= 0 && dp_y_sd    <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	  }
	  
	  /* R matrix: */
	  /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	   * this reduces to two ESL_MAX calls
	   */
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sdr);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sdr);
	  dpn       = dn - hdmin[v][jp_v];
	  dpx       = dx - hdmin[v][jp_v];
	  dp_y_sdr  = dn - hdmin[y][jp_y_sdr] - sdr;
	  /* for {L,R}alpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
	    /* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
	    ESL_DASSERT1((dp_y_sdr  >= 0 && dp_y_sdr   <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], 
					    ESL_MAX(Jalpha[y][jp_y_sdr][dp_y_sdr] + tsc, 
						    Ralpha[y][jp_y_sdr][dp_y_sdr] + tsc)); 
	  }
	}

	/* The second MP_st 'for (jp_v...' loop is for L matrix, which uses a different set of j values from J and R */
	/* j must satisfy:
	 * j >= jmin[v]
	 * j >= jmin[y] (follows from (j >= jmin[y]))
	 * j <= jmax[v]
	 * j <= jmax[y] (follows from (j <= jmax[y]))
	 * this reduces to two ESL_MAX calls
	 */
	jn = ESL_MAX(jmin[v], jmin[y]);
	jx = ESL_MIN(jmax[v], jmax[y]);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y = jn - jmin[y];
	/* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	
	for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y++) {
	  ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	  ESL_DASSERT1((jp_y     >= 0 && jp_y     <= (jmax[y]-jmin[y])));
	  
	  /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	   * this reduces to two ESL_MAX calls
	   */
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y] + sdr);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y] + sdr);
	  dpn       = dn - hdmin[v][jp_v];
	  dpx       = dx - hdmin[v][jp_v];
	  dp_y_sdr  = dn - hdmin[y][jp_y] - sdr;
	  /* for Lalpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
	    /* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
	    ESL_DASSERT1((dp_y_sdr >= 0 && dp_y_sdr  <= (hdmax[y][jp_y]     - hdmin[y][jp_y])));
	    Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], 
					    ESL_MAX(Jalpha[y][jp_y][dp_y_sdr] + tsc, 
						    Lalpha[y][jp_y][dp_y_sdr] + tsc)); 
	  }
	}
      }
      /* add in emission score */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	i     = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v], dp_v = 0; d <= hdmax[v][jp_v]; d++, dp_v++) 
	  {
	    /*if(i < i0 || j > j0) { 
	      printf("dsq[i:%d]: %d\n", i, dsq[i]);
	      printf("dsq[j:%d]: %d\n", j, dsq[j]);
	      printf("esc_v[%d]: %.5f\n", dsq[i]*cm->abc->Kp+dsq[j], esc_v[dsq[i]*cm->abc->Kp+dsq[j]]);;
	      printf("i0: %d j0: %d\n", i0, j0);
	      }*/
	    if(d >= 2) { 
	      Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	      Lalpha[v][jp_v][dp_v] += lmesc_v[dsq[i]];
	      Ralpha[v][jp_v][dp_v] += rmesc_v[dsq[j]];;
	    }
	    else { 
	      Jalpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      Lalpha[v][jp_v][dp_v] = lmesc_v[dsq[i]];
	      Ralpha[v][jp_v][dp_v] = rmesc_v[dsq[j]];
	    }
	    i--;
	  }
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	  Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is D or S (! E && ! B && ! ML && ! IL && ! MR && ! IR) */
      /* D, S states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	/* j must satisfy:
	 * j >= jmin[v]
	 * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	 * j <= jmax[v]
	 * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	 * this reduces to two ESL_MAX calls
	 */
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	
	for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
	  ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	  ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	  
	  /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	   * this reduces to two ESL_MAX calls
	   */
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	  dpn     = dn - hdmin[v][jp_v];
	  dpx     = dx - hdmin[v][jp_v];
	  dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	  	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	    Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	    Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y_sd] + tsc);
	  }
	}
      }
      /* no emission score to add */
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      /* Any valid j must be within both state v and state z's j band 
       * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
       * constructed, but we'll check anyway. 
       */
      jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
      jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
      /* the main j loop */
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	kn = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	i = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++, i--) {
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	  /* Find the first k value that implies a valid cell in the {J,L,R} matrix y and z decks.
	   * This k must satisfy the following 6 inequalities (some may be redundant):
	   * (1) k >= j-jmax[y];
	   * (2) k <= j-jmin[y]; 
	   *     1 and 2 guarantee (j-k) is within state y's j band
	   *
	   * (3) k >= hdmin[z][j-jmin[z]];
	   * (4) k <= hdmax[z][j-jmin[z]]; 
	   *     3 and 4 guarantee k is within z's j=(j), d band
	   *
	   * (5) k >= d-hdmax[y][j-jmin[y]-k];
	   * (6) k <= d-hdmin[y][j-jmin[y]-k]; 
	   *     5 and 6 guarantee (d-k) is within state y's j=(j-k) d band
	   *
	   * kn and kx were set above (outside (for (dp_v...) loop) that
	   * satisfy 1-4 (b/c 1-4 are d-independent and k-independent)
	   * RHS of inequalities 5 and 6 are dependent on k, so we check
	   * for these within the next for loop.
	   *
	   * To update a cell in the T matrix with a sum of an R matrix value for y
	   * and a L matrix value for z, there are 2 additional inequalities to satisfy:
	   * (7) k != i-1  (where i = j-d+1)
	   * (8) k != j
	   * These are checked for in the loop below as well. 
	   */
	  for(k = kn; k <= kx; k++) { 
	    if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k]) {
	      /* for current k, all 6 inequalities have been satisified 
	       * so we know the cells corresponding to the platonic 
	       * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
	       * alpha[z][j][k] are all within the bands. These
	       * cells correspond to alpha[v][jp_v][dp_v], 
	       * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
	       * and alpha[z][jp_z][k-hdmin[jp_z]];
	       */
	      kp_z = k-hdmin[z][jp_z];
	      dp_y = d-hdmin[y][jp_y-k];
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]);
	      Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], Jalpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]);
	      Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]);
	      /*if((k != i-1) && (k != j)) {*/
	      Talpha[v][jp_v][dp_v] = ESL_MAX(Talpha[v][jp_v][dp_v], Ralpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]);
		/*}*/
	    }
	  }
	}
      }

      /* two additional special cases in trCYK (these are not in standard CYK).
       * we do these in their own for(j.. { for(d.. { } } loops b/c one 
       * is independent of z, the other of y, unlike the above loop which is dependent 
       * on both.
       */
      jn = (jmin[v] > jmin[y]) ? jmin[v] : jmin[y];
      jx = (jmax[v] < jmax[y]) ? jmax[v] : jmax[y];
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	ESL_DASSERT1((j >= jmin[y] && j <= jmax[y]));
	dn = (hdmin[v][jp_v] > hdmin[y][jp_y]) ? hdmin[v][jp_v] : hdmin[y][jp_y];
	dx = (hdmax[v][jp_v] < hdmax[y][jp_y]) ? hdmax[v][jp_v] : hdmax[y][jp_y];
	for(d = dn; d <= dx; d++) { 
	  dp_v = d - hdmin[v][jp_v];
	  dp_y = d - hdmin[y][jp_y];
	  ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	  ESL_DASSERT1((d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y]));
	  Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], 
					  ESL_MAX(Jalpha[y][jp_y][dp_y],
						  Lalpha[y][jp_y][dp_y]));
	}
      }
      jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
      jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_z = j - jmin[z];
	ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	ESL_DASSERT1((j >= jmin[z] && j <= jmax[z]));
	dn = (hdmin[v][jp_v] > hdmin[z][jp_z]) ? hdmin[v][jp_v] : hdmin[z][jp_z];
	dx = (hdmax[v][jp_v] < hdmax[z][jp_z]) ? hdmax[v][jp_v] : hdmax[z][jp_z];
	for(d = dn; d <= dx; d++) { 
	  dp_v = d - hdmin[v][jp_v];
	  dp_z = d - hdmin[z][jp_z];
	  ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	  ESL_DASSERT1((d >= hdmin[z][jp_z] && d <= hdmax[z][jp_z]));
	  Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], 
					  ESL_MAX(Jalpha[z][jp_z][dp_z],
						  Ralpha[z][jp_z][dp_z]));
	}
      }
    } /* finished calculating deck v. */
#if PRINTFALPHA
	  if(cm->stid[v] == BIF_B) { 
	    /* the main j loop */
	    for (j = jmin[v]; j <= jmax[v]; j++) { 
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		printf("H j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		       j, v, d, 
		       NOT_IMPOSSIBLE(Jalpha[v][jp_v][dp_v]) ? Jalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Lalpha[v][jp_v][dp_v]) ? Lalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Ralpha[v][jp_v][dp_v]) ? Ralpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Talpha[v][jp_v][dp_v]) ? Talpha[v][jp_v][dp_v] : -9999.9);
	      }
	    }
	  }
	  if((cm->stid[v] == BEGL_S) || (cm->stid[v] == BEGR_S)) { 
	    /* the main j loop */
	    for (j = jmin[v]; j <= jmax[v]; j++) { 
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		printf("H j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		       j, v, d, 
		       NOT_IMPOSSIBLE(Jalpha[v][jp_v][dp_v]) ? Jalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Lalpha[v][jp_v][dp_v]) ? Lalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Ralpha[v][jp_v][dp_v]) ? Ralpha[v][jp_v][dp_v] : -9999.9);
	      }
	    }
	  }
#endif
#if 0 
	  else { 
	    for (j = jmin[v]; j <= jmax[v]; j++) { 
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		printf("H j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		       j, v, d, 
		       NOT_IMPOSSIBLE(Jalpha[v][jp_v][dp_v]) ? Jalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Lalpha[v][jp_v][dp_v]) ? Lalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Ralpha[v][jp_v][dp_v]) ? Ralpha[v][jp_v][dp_v] : -9999.9, 
		       -9999.9);
	      }
	    }
	  }
	  printf("\n");
#endif
  } /* end of for (v = cm->M-1; v > 0; v--) */
        
  /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
   * 
   * If local begins are off, all hits must be rooted at v=0.
   * With local begins on, the hit is rooted at the second state in
   * the traceback (e.g. after 0), the internal entry point. 
   * 
   * Hits rooted at 0 that not involved with local begins are 
   * already calc'ed from the v loop with v == 0 
   */

  /* Report all possible hits, but only after looking at local begins (if they're on) */
  v = 0;
  sd = sdr = 0;
  jpn = 0;
  jpx = jmax[v] - jmin[v];
  j   = jmin[v];
  
  ESL_ALLOC(bestr, sizeof(int) * (W+1));
  ESL_ALLOC(bestmode, sizeof(int) * (W+1));
  esl_vec_ISet(bestr,    W+1, 0);
  esl_vec_ISet(bestmode, W+1, CM_HIT_MODE_J);

  /* first report all hits with j < jmin[0] are impossible, only if we're reporting hits to hitlist */
  if(hitlist != NULL) { 
    for(j = i0; j < jmin[v]; j++) {
      if((status = UpdateGammaHitMx(cm, errbuf, gamma, j-i0+1, 
				    NULL, 0, 0,   /* alpha_row is NULL, we're telling UpdateGammaHitMx, this j can't be a hit end pt */
				    TRUE, bestr, bestmode, hitlist, W, act)) != eslOK) return status;
    }
  }

  /* consider local hits */
  v = 0;
  for (j = jmin[v]; j <= jmax[v]; j++) {
    jp_v = j - jmin[v];
    esl_vec_ISet(bestr, (W+1), 0); /* init bestr to 0, all hits are rooted at 0 unless we find a better local begin below */
    if (cm->flags & CMH_LOCAL_BEGIN) {
      for (y = 1; y < cm->M; y++) {
	if(NOT_IMPOSSIBLE(cm->beginsc[y]) && /* local begin is possible into y */
	   (j >= jmin[y] && j <= jmax[y])) { /* j is within state y's j band */
	  jp_y = j - jmin[y];
	  dn   = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y]);
	  dx   = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y]);
	  dp_v = dn - hdmin[v][jp_v];
	  dp_y = dn - hdmin[y][jp_y];
	  for(d = dn; d <= dx; d++, dp_v++, dp_y++) { 
	    sc = Jalpha[y][jp_y][dp_y] + cm->beginsc[y];
	    if(sc > Jalpha[0][jp_v][dp_v]) {
	      Jalpha[0][jp_v][dp_v] = sc;
	      bestr[d] = y;
	    }
	  }
	}
      } /* end of for(y = 1; y < cm->M; y++) */
    } /* end of if(cm->flags & CMH_LOCAL_BEGIN */
  }    

  /* Finally, allow for truncated hits */
  v = 0;
  for (j = jmin[v]; j <= jmax[v]; j++) {
    jp_v = j - jmin[v];
    esl_vec_ISet(bestmode, (W+1), CM_HIT_MODE_J); /* init bestmode to CM_HIT_MODE_J */
    for (y = 1; y < cm->M; y++) {
      if(j >= jmin[y] && j <= jmax[y]) { /* j is within state y's band */
	jp_y = j - jmin[y];
	if(cm->sttype[y] == MP_st || cm->sttype[y] == ML_st || cm->sttype[y] == MR_st) { 
	  /*if(NOT_IMPOSSIBLE(cm->beginsc[y])) {*/
	  /*trunc_penalty = cm->beginsc[y];*/
	  dn   = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y]);
	  dx   = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y]);
	  dp_v = dn - hdmin[v][jp_v];
	  dp_y = dn - hdmin[y][jp_y];
	  for(d = dn; d <= dx; d++, dp_v++, dp_y++) { 
	    Jsc = Jalpha[y][jp_y][dp_y] + trunc_penalty; 
	    Lsc = Lalpha[y][jp_y][dp_y] + trunc_penalty; 
	    Rsc = Ralpha[y][jp_y][dp_y] + trunc_penalty; 
	    if(Jsc > Jalpha[0][jp_v][dp_v]) { 
	      Jalpha[0][jp_v][dp_v] = Jsc;
	      bestr[d] = y;
	      /* bestmode[d] remains CM_HIT_MODE_J */
	    }
	    if(Lsc > Jalpha[0][jp_v][dp_v] && (cm->sttype[y] != MR_st)) { 
	      Jalpha[0][jp_v][dp_v] = Lsc;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_L;
	    }
	    if(Rsc > Jalpha[0][jp_v][dp_v] && (cm->sttype[y] != ML_st)) { 
	      Jalpha[0][jp_v][dp_v] = Rsc;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_R;;
	    }
	  }
	}
	else if(cm->sttype[y] == B_st) { 
	  dn   = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y]);
	  dx   = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y]);
	  dp_v = dn - hdmin[v][jp_v];
	  dp_y = dn - hdmin[y][jp_y];
	  for(d = dn; d <= dx; d++, dp_v++, dp_y++) { 
	    Tsc = Talpha[y][jp_y][dp_y] + trunc_penalty; 
	    if(Tsc > Jalpha[0][jp_v][dp_v]) {
	      Jalpha[0][jp_v][dp_v] = Tsc;
	      bestr[d] = y;
	      bestmode[d] = CM_HIT_MODE_T; 
	    }
	  }
	}
      }
    }
  }
    
  /* report all hits with valid d for this j, only if hitlist != NULL */
  if(hitlist != NULL) { 
    if((status = UpdateGammaHitMx(cm, errbuf, gamma, j-i0+1, Jalpha[0][jp_v], hdmin[0][jp_v], hdmax[0][jp_v], TRUE, bestr, bestmode, hitlist, W, act)) != eslOK) return status;
  }
  /* finally report all hits with j > jmax[0] are impossible, only if we're reporting hits to hitlist */
  if(hitlist != NULL) { 
    for(j = jmax[v]+1; j <= j0; j++) {
      if((status = UpdateGammaHitMx(cm, errbuf, gamma, j-i0+1, 
				    NULL, 0, 0,   /* alpha_row is NULL, we're telling UpdateGammaHitMx, this j is can't be a hit end pt */
				    TRUE, bestr, bestmode, hitlist, W, act)) != eslOK) return status;
    }
  }
  /* find the best scoring hit, and update envelope boundaries if nec */
  vsc_root = IMPOSSIBLE;
  v = 0;
  jpn = 0;
  jpx = jmax[v] - jmin[v];
  for(jp_v = jpn; jp_v <= jpx; jp_v++) {
    dpn     = 0;
    dpx     = hdmax[v][jp_v] - hdmin[v][jp_v];
    for(dp_v = dpn; dp_v <= dpx; dp_v++) {
      vsc_root = ESL_MAX(vsc_root, Jalpha[0][jp_v][dp_v]);
    }
    /* update envelope boundaries, if nec */
    if(do_env_defn) { 
      j = jp_v + jmin[v];
      for(dp_v = dpn; dp_v <= dpx; dp_v++) {
	if(Jalpha[0][jp_v][dp_v] >= env_cutoff) { 
	  i = j - (dp_v + hdmin[v][jp_v]) + 1;
	  envi = ESL_MIN(envi, i);
	  envj = ESL_MAX(envj, j);
	}
      }
    }
  }
  /* find the best score in any matrix at any state */
  float best_tr_sc = IMPOSSIBLE;
  best_tr_sc = vsc_root;
  printf("Best truncated score: %.4f (%.4f)\n",
	 best_tr_sc - trunc_penalty, 
	 best_tr_sc + sreLOG2(2./(cm->clen * (cm->clen+1))));

  free(el_scA);
  free(yvalidA);
  free(bestr);
  free(bestmode);
  if (act != NULL) { 
    for(i = 0; i <= W; i++) free(act[i]); 
    free(act);
  }

  if(hitlist != NULL && gamma->iamgreedy == FALSE) TBackGammaHitMx(gamma, hitlist, i0, j0);

  /* set envelope return variables if nec */
  if(ret_envi != NULL) { *ret_envi = (envi == j0+1) ? -1 : envi; }
  if(ret_envj != NULL) { *ret_envj = (envj == i0-1) ? -1 : envj; }

  if(gamma != NULL) FreeGammaHitMx(gamma);

  if (ret_sc != NULL) *ret_sc = vsc_root;
  ESL_DPRINTF1(("FastCYKScanHB() return sc: %f\n", vsc_root));
  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return 0.; /* never reached */
}


/* Function: FastTrCYKScanHBSkipSlow()
 * Incept:   EPN, Thu Aug 25 15:19:28 2011
 *
 * Purpose:  An HMM banded version of a scanning TrCYK algorithm. Takes
 *           a CM_TR_HB_MX data structure which is indexed [v][j][d] with
 *           only cells within the bands allocated.
 *           (different than other (non-HB) scanning function's convention 
 *            of [j][v][d]).
 *
 * Args:     cm        - the model    [0..M-1]
 *           sq        - the sequence [1..L]   
 *                     - length of the dsq
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           cutoff    - minimum score to report
 *           hitlist   - CM_TOPHITS hitlist to add to; if NULL, don't add to it
 *           do_null3  - TRUE to do NULL3 score correction, FALSE not to
 *           env_cutoff- ret_envi..ret_envj will include all hits that exceed this bit sc
 *           ret_envi  - min position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted
 *           ret_envj  - max position in any hit w/sc >= env_cutoff, set to -1 if no such hits exist, NULL if not wanted 
 *           mx        - the dp matrix, only cells within bands in cm->cp9b will 
 *                       be valid. 
 *           size_limit- max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           ret_sc    - RETURN: score of best overall hit (vsc[0])
 *                       
 * Returns: eslOK on success; eslERANGE if required matrix size exceeds <size_limit>
 *          <ret_sc>: score of the best hit.
 */
int
FastTrCYKScanHBSkipSlow(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, 
			CM_TR_HB_MX *mx, float size_limit, float env_cutoff, CMEmitMap_t *emap, int64_t *ret_envi, int64_t *ret_envj, float *ret_sc)
{
  int      status;
  GammaHitMx_t *gamma;  /* semi-HMM for hit resoultion */
  int     *bestr;       /* best root state for d at current j */
  int     *bestmode;    /* best mode for parsetree for d at current j */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc, Jsc, Lsc, Rsc, Tsc; /* temporary scores */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      j_sdr;              /* j - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y, dp_z;   /* offset d index for states v, y, z */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp;                 /* ESL_MAX(d-sd, 0) */
  int      dp_y_sd;            /* dp_y - sd */
  int      dp_y_sdr;           /* dp_y - sdr, often for jp_y_sdr */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  float    tsc;                /* a transition score */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */
  float    vsc_root;           /* score of best hit */
  int      W;                  /* max d over all hdmax[v][j] for all valid v, j */
  double **act;                /* [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1) */
  int      jp;                 /* j index in act */
  int      do_env_defn;        /* TRUE to calculate envi, envj, FALSE not to (TRUE if ret_envi != NULL or ret_envj != NULL */
  int64_t  envi, envj;         /* min/max positions that exist in any hit with sc >= env_cutoff */
  float    trunc_penalty = 0.; /* penalty in bits for a truncated hit */
  int      nd, lpos, rpos;     /* node index, left and right consensus position boundaries of subtree rooted at v */
  int      do_J, do_L, do_R, do_T; /* TRUE to fill in values for J matrix L matrix, R matrix, T matrix for current v */

  /* Contract check */
  if(dsq == NULL)       ESL_FAIL(eslEINCOMPAT, errbuf, "FastTrCYKScanHB(), dsq is NULL.\n");
  if (mx == NULL)       ESL_FAIL(eslEINCOMPAT, errbuf, "FastTrCYKScanHB(), mx is NULL.\n");
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "FastTrCYKScanHB(), mx is NULL.\n");

  ESL_DPRINTF1(("cm->search_opts & CM_SEARCH_HMMALNBANDS: %d\n", cm->search_opts & CM_SEARCH_HMMALNBANDS));

  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  CP9Bands_t *cp9b   = cm->cp9b; 
  int        *jmin   = cp9b->jmin;  
  int        *jmax   = cp9b->jmax;
  int       **hdmin  = cp9b->hdmin;
  int       **hdmax  = cp9b->hdmax;
  /* the DP matrix */
  float ***Jalpha  = mx->Jdp; /* pointer to the Jalpha DP matrix */
  float ***Lalpha  = mx->Ldp; /* pointer to the Lalpha DP matrix */
  float ***Ralpha  = mx->Rdp; /* pointer to the Ralpha DP matrix */
  float ***Talpha  = mx->Tdp; /* pointer to the Talpha DP matrix */

  /* Allocations and initializations  */
  /* grow the matrix based on the current sequence and bands */
  if((status = cm_tr_hb_mx_GrowTo(cm, mx, errbuf, cp9b, (j0-i0+1), size_limit)) != eslOK) return status;

  /* set W as j0-i0+1 (this may exceed max size of a hit our bands will allow, 
   * but that's okay b/c W is only used for sizing of act and bestr vectors */
  W = j0-i0+1;
  /* make sure our bands won't allow a hit bigger than W (this could be modified to only execute in debugging mode) */
  for(j = jmin[0]; j <= jmax[0]; j++) {
    if(W < (hdmax[0][(j-jmin[0])])) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "FastCYKScanHB(), band allows a hit (j:%d hdmax[0][j]:%d) greater than j0-i0+1 (%d)", j, hdmax[0][(j-jmin[0])], j0-i0+1);
  }

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  ESL_STOPWATCH *w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);
  /* initialize all cells of the matrix to IMPOSSIBLE */
  esl_vec_FSet(Jalpha[0][0], mx->JLRncells_valid, IMPOSSIBLE);
  esl_vec_FSet(Lalpha[0][0], mx->JLRncells_valid, IMPOSSIBLE);
  esl_vec_FSet(Ralpha[0][0], mx->JLRncells_valid, IMPOSSIBLE);
  if(mx->Tncells_valid > 0)   esl_vec_FSet(mx->Tdp_mem, mx->Tncells_valid,   IMPOSSIBLE); /* careful use Tdp_mem not Talpha[0][0] because Talpha[0] is NULL */
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " Matrix init CPU time: ");

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse of multiple nonoverlapping hits. */
  if(hitlist != NULL) gamma = CreateGammaHitMx(j0-i0+1, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff, FALSE);
  else                gamma = NULL;

  /* if do_null3: allocate and initialize act vector */
  if(do_null3) { 
    ESL_ALLOC(act, sizeof(double *) * (W+1));
    for(i = 0; i <= W; i++) { 
      ESL_ALLOC(act[i], sizeof(double) * cm->abc->K);
      esl_vec_DSet(act[i], cm->abc->K, 0.);
    }
    /* pre-fill act, different than non-HMM banded scanner b/c our main loop doesn't step j through residues */
    for(j = i0; j <= j0; j++) { 
      jp = j-i0+1; /* j is actual index in dsq, jp_g is offset j relative to start i0 (j index for act) */
      esl_vec_DCopy(act[(jp-1)%(W+1)], cm->abc->K, act[jp%(W+1)]);
      esl_abc_DCount(cm->abc, act[jp%(W+1)], dsq[j], 1.);
    }
  }
  else act = NULL;

  /* initialize envelope boundary variables */
  do_env_defn = (ret_envi != NULL || ret_envj != NULL) ? TRUE : FALSE;
  envi = j0+1;
  envj = i0-1;

  /* Main recursion */
  for (v = cm->M-1; v >= 0; v--) { /* all the way down to root, different from other scanners */
    float const *esc_v   = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v   = cm->tsc[v];  /* transition scores for state v */
    float const *lmesc_v = cm->lmesc[v];
    float const *rmesc_v = cm->rmesc[v];
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];
    nd   = cm->ndidx[v];
    /* Careful, emitmap is off-by-one for our purposes for lpos if v is not MATP_MP or MATL_ML, and rpos if v is not MATP_MP or MATR_MR */
    lpos = (cm->stid[v] == MATP_MP || cm->stid[v] == MATL_ML) ? emap->lpos[nd] : emap->lpos[nd]+1;
    rpos = (cm->stid[v] == MATP_MP || cm->stid[v] == MATR_MR) ? emap->rpos[nd] : emap->rpos[nd]-1;

    /* lpos almost certainly not used OR  rpos almost certainly not used => don't do_J, else do */
    do_J = ((cp9b->occ[lpos] < (1.0-cp9b->occthresh)) || (cp9b->occ[rpos] < (1.0-cp9b->occthresh))) ? FALSE : TRUE; 

    /* lpos almost certainly not used OR  rpos almost certainly     used => don't do_L, else do */
    do_L = ((cp9b->occ[lpos] < (1.0-cp9b->occthresh)) || (cp9b->occ[rpos] > cp9b->occthresh))       ? FALSE : TRUE;                                                     

    /* rpos almost certainly not used OR  lpos almost certainly     used => don't do_R, else do */
    do_R = ((cp9b->occ[rpos] < (1.0-cp9b->occthresh)) || (cp9b->occ[lpos] > cp9b->occthresh))       ? FALSE : TRUE;

    /* lpos almost certainly     used OR  rpos almost certainly     used => don't do_T, else do */
    /* COULD ALSO CHECK TO SEE IF ALL POSITIONS lpos..rpos are almost certain not used, if so, do_T becomes FALSE */
    do_T = ((cp9b->occ[lpos] > cp9b->occthresh)       || (cp9b->occ[rpos] > cp9b->occthresh))       ? FALSE : TRUE;

    printf("v: %4d %4s %2s [%4d..%4d] (%6.4f..%6.4f) J: %d L: %d R: %d T: %d\n", v, Nodetype(cm->ndtype[nd]), Statetype(cm->sttype[v]), lpos, rpos, cp9b->occ[lpos], cp9b->occ[rpos], do_J, do_L, do_R, do_T);

    /* re-initialize the deck if we can do a local end from v */
    if(do_J) { 
      if(NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	  dp = ESL_MAX(d-sd, 0);
	  Jalpha[v][jp_v][dp_v] = el_scA[dp] + cm->endsc[v];
	  /* L,Ralpha[v] remain IMPOSSIBLE, they can't go to EL */
	  }
	}
      }
    }
    /* otherwise this state's deck has already been initialized to IMPOSSIBLE */

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	if(do_J) Jalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	if(do_L) Lalpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
	if(do_R) Ralpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IL states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1] 
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr valid for state y? */
	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	  i = j - d + 1;
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */

	  /* We need to treat R differently from and J and L here, by
	   * doing separate 'for (yoffset...' loops for J because we
	   * have to fully calculate Jalpha[v][jp_v][dp_v]) before we
	   * can start to calculate Ralpha[v][jp_v][dp_v].
	   */
	  /* Handle J and L first */
	  if(do_J || do_L) { 
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      jp_y_sdr = j - jmin[y] - sdr;
	      
	      if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		if(do_J) Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
		if(do_L) Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	      }
	    }
	    if(do_J) { 
	      Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]];
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    if(do_L) { 
	      Lalpha[v][jp_v][dp_v] = (d >= 2) ? Lalpha[v][jp_v][dp_v] + esc_v[dsq[i]]: esc_v[dsq[i]];
	      Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    i--;
	  }

	  if(do_R) { 
	    /* Handle R separately */
	    Rsc = Ralpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      jp_y_sdr = j - jmin[y] - sdr;
	      
	      /* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,L above) */
	      if((d) >= hdmin[y][jp_y_sdr] && (d) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		dp_y = d - hdmin[y][jp_y_sdr];
		ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		Rsc = ESL_MAX(Rsc, ESL_MAX(Jalpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset],
					   Ralpha[y][jp_y_sdr][dp_y] + tsc_v[yoffset]));
	      }
	    } /* end of for (yvalid_idx = 0... loop */
	    Ralpha[v][jp_v][dp_v] = Rsc; 
	    /* we use Rsc instead of Ralpha cell in above loop because
	     * Ralpha[v][jp_v][dp_v] may be the same cell as
	     * Ralpha[y][jp_y_sdr][dp_y] if we're an IL state
	     */
	  }
	}
      }
    }
    else if(cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) { 
      /* update {J,L,R}alpha[v][jp_v][dp_v] cells, for IR states, loop
       * nesting order is: for j { for d { for y { } } } because they
       * can self transit, and a {J,L,R}alpha[v][j][d] cell must be
       * complete (that is we must have looked at all children y)
       * before can start calc'ing for {J,L,R}alpha[v][j][d+1].
       * We could be slightly more efficient if we separated out 
       * MR from IR b/c self-transits in MRs are impossible, but 
       * we don't do that here. */

      /* The first MR_st/IR_st 'for (j...' loop is for J and R matrices which use the same set of j values */
      if(do_J || do_R) { 
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  j_sdr = j - sdr;
	  
	  /* determine which children y we can legally transit to for v, j */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	  
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	    /* We need to treat L differently from and J and R here, by
	     * doing separate 'for (yoffset...' loops for J because we
	     * have to fully calculate Jalpha[v][jp_v][dp_v]) before we
	     * can start to calculate Lalpha[v][jp_v][dp_v].
	     */
	    /* Handle J and R first */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      jp_y_sdr = j - jmin[y] - sdr;
	      
	      if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
		dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
		ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		if(do_J) Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
		if(do_R) Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]);
	      }
	    }
	    if(do_J) { 
	      Jalpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	      Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	    if(do_R) { 
	      Ralpha[v][jp_v][dp_v] = (d >= 2) ? Ralpha[v][jp_v][dp_v] + esc_v[dsq[j]] : esc_v[dsq[j]];
	      Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	    }
	  }
	}
      }

      if(do_L) { 
	/* The second MR_st/IR_st 'for (j...' loop is for the L matrix which use a different set of j values */
	for (j = jmin[v]; j <= jmax[v]; j++) {
	  jp_v = j - jmin[v];
	  yvalid_ct = 0;
	  
	  /* determine which children y we can legally transit to for v, j */
	  /* we use 'j' and not 'j_sdr' here for the L matrix, differently from J and R matrices above */
	  for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	    if((j) >= jmin[y] && ((j) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j is valid for state y? */
	  
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	    dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	    
	    Lsc = Lalpha[v][jp_v][dp_v]; /* this sc will be IMPOSSIBLE */
	    for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	      yoffset = yvalidA[yvalid_idx];
	      y = cm->cfirst[v] + yoffset;
	      /* we use 'jp_y=j-min[y]' here, not 'jp_y_sdr=j-jmin[y]-sdr' (which we used in the corresponding loop for J,R above) */
	      jp_y = j - jmin[y];
	      
	      /* we use 'd' and 'dp_y' here, not 'd-sd' and 'dp_y_sd' (which we used in the corresponding loop for J,R above) */
	      if((d) >= hdmin[y][jp_y] && (d) <= hdmax[y][jp_y]) { /* make sure d is valid for this v, j and y */
		dp_y = d - hdmin[y][jp_y];
		ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v] - hdmin[v][jp_v])));
		ESL_DASSERT1((dp_y    >= 0 && dp_y     <= (hdmax[y][jp_y] - hdmin[y][jp_y])));
		Lsc = ESL_MAX(Lsc, ESL_MAX(Jalpha[y][jp_y][dp_y] + tsc_v[yoffset], 
					   Lalpha[y][jp_y][dp_y] + tsc_v[yoffset]));
	      }
	    } /* end of for (yvalid_idx = 0... loop */
	    Lalpha[v][jp_v][dp_v] = Lsc; 
	    /* we use Lsc instead of Lalpha cell in above loop because
	     * Lalpha[v][jp_v][dp_v] may be the same cell as
	     * Lalpha[y][jp_y_sdr][dp_y] if we're an IR state
	   */
	  }
	}
      }
    }
    else if(cm->sttype[v] == MP_st) { 
      /* MP states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	/* The first MP_st 'for (jp_v...' loop is for J and R matrices which use the same set of j values */
	/* j must satisfy:
	 * j >= jmin[v]
	 * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	 * j <= jmax[v]
	 * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	 * this reduces to two ESL_MAX calls
	 */
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	/* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	
	if(do_J || do_R) { 
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++, jp_y++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	    
	    if(do_J) { 
	      /* J matrix: */
	      /* d must satisfy:
	       * d >= hdmin[v][jp_v]
	       * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	       * d <= hdmax[v][jp_v]
	       * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	       * this reduces to two ESL_MAX calls
	       */
	      dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	      dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	      dpn       = dn - hdmin[v][jp_v];
	      dpx       = dx - hdmin[v][jp_v];
	      dp_y_sd   = dn - hdmin[y][jp_y_sdr] - sd;
	      
	      for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
		ESL_DASSERT1((dp_v      >= 0 && dp_v       <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
		ESL_DASSERT1((dp_y_sd   >= 0 && dp_y_sd    <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	      }
	    }
	    
	    if(do_R) { 
	      /* R matrix: */
	      /* d must satisfy:
	       * d >= hdmin[v][jp_v]
	       * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	       * d <= hdmax[v][jp_v]
	       * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	       * this reduces to two ESL_MAX calls
	       */
	      dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sdr);
	      dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sdr);
	      dpn       = dn - hdmin[v][jp_v];
	      dpx       = dx - hdmin[v][jp_v];
	      dp_y_sdr  = dn - hdmin[y][jp_y_sdr] - sdr;
	      /* for {L,R}alpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	      
	      for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
		/* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
		ESL_DASSERT1((dp_y_sdr  >= 0 && dp_y_sdr   <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
		Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], 
						ESL_MAX(Jalpha[y][jp_y_sdr][dp_y_sdr] + tsc, 
						      Ralpha[y][jp_y_sdr][dp_y_sdr] + tsc)); 
	      }
	    }
	  }
	}

	if(do_L) { 
	  /* The second MP_st 'for (jp_v...' loop is for L matrix, which uses a different set of j values from J and R */
	  /* j must satisfy:
	   * j >= jmin[v]
	   * j >= jmin[y] (follows from (j >= jmin[y]))
	   * j <= jmax[v]
	   * j <= jmax[y] (follows from (j <= jmax[y]))
	   * this reduces to two ESL_MAX calls
	   */
	  jn = ESL_MAX(jmin[v], jmin[y]);
	  jx = ESL_MIN(jmax[v], jmax[y]);
	  jpn = jn - jmin[v];
	  jpx = jx - jmin[v];
	  jp_y = jn - jmin[y];
	  /* for Lalpha, we use 'jp_y=j-min[y]' instead of 'jp_y_sdr=j-jmin[y]-sdr' */
	  
	  for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y++) {
	    ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	    ESL_DASSERT1((jp_y     >= 0 && jp_y     <= (jmax[y]-jmin[y])));
	    
	    /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sd]+sd (follows from (d-sd >= hdmin[y][jp_y_sd]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sd]+sd (follows from (d-sd <= hdmax[y][jp_y_sd]))
	   * this reduces to two ESL_MAX calls
	   */
	    dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y] + sdr);
	    dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y] + sdr);
	    dpn       = dn - hdmin[v][jp_v];
	    dpx       = dx - hdmin[v][jp_v];
	    dp_y_sdr  = dn - hdmin[y][jp_y] - sdr;
	    /* for Lalpha, we use 'dp_y_sdr' instead of 'dy_y_sd' */
	    
	    for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sdr++) { 
	      /* we use 'dp_y_sdr' here, not 'dp_y_sd' (which we used in the corresponding loop for J above) */
	      ESL_DASSERT1((dp_y_sdr >= 0 && dp_y_sdr  <= (hdmax[y][jp_y]     - hdmin[y][jp_y])));
	      Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], 
					      ESL_MAX(Jalpha[y][jp_y][dp_y_sdr] + tsc, 
						      Lalpha[y][jp_y][dp_y_sdr] + tsc)); 
	    }
	  }
	}
      }
      /* add in emission score */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	i     = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v], dp_v = 0; d <= hdmax[v][jp_v]; d++, dp_v++) 
	  {
	    /*if(i < i0 || j > j0) { 
	      printf("dsq[i:%d]: %d\n", i, dsq[i]);
	      printf("dsq[j:%d]: %d\n", j, dsq[j]);
	      printf("esc_v[%d]: %.5f\n", dsq[i]*cm->abc->Kp+dsq[j], esc_v[dsq[i]*cm->abc->Kp+dsq[j]]);;
	      printf("i0: %d j0: %d\n", i0, j0);
	      }*/
	    if(d >= 2) { 
	      if(do_J) Jalpha[v][jp_v][dp_v] += esc_v[dsq[i]*cm->abc->Kp+dsq[j]];
	      if(do_L) Lalpha[v][jp_v][dp_v] += lmesc_v[dsq[i]];
	      if(do_R) Ralpha[v][jp_v][dp_v] += rmesc_v[dsq[j]];;
	    }
	    else { 
	      if(do_J) Jalpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      if(do_L) Lalpha[v][jp_v][dp_v] = lmesc_v[dsq[i]];
	      if(do_R) Ralpha[v][jp_v][dp_v] = rmesc_v[dsq[j]];
	    }
	    i--;
	  }
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	  if(do_J) Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  if(do_L) Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], IMPOSSIBLE);
	  if(do_R) Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is D or S (! E && ! B && ! ML && ! IL && ! MR && ! IR) */
      /* D, S states cannot self transit, this means that all cells in
       * alpha[v] are independent of each other, only depending on
       * alpha[y] for previously calc'ed y.  We can do the for loops
       * in any nesting order, this implementation does what I think
       * is most efficient: for y { for j { for d { } } }
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	/* j must satisfy:
	 * j >= jmin[v]
	 * j >= jmin[y]+sdr (follows from (j-sdr >= jmin[y]))
	 * j <= jmax[v]
	 * j <= jmax[y]+sdr (follows from (j-sdr <= jmax[y]))
	 * this reduces to two ESL_MAX calls
	 */
	jn = ESL_MAX(jmin[v], jmin[y]+sdr);
	jx = ESL_MIN(jmax[v], jmax[y]+sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	
	for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
	  ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	  ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	  
	  /* d must satisfy:
	   * d >= hdmin[v][jp_v]
	   * d >= hdmin[y][jp_y_sdr]+sd (follows from (d-sd >= hdmin[y][jp_y_sdr]))
	   * d <= hdmax[v][jp_v]
	   * d <= hdmax[y][jp_y_sdr]+sd (follows from (d-sd <= hdmax[y][jp_y_sdr]))
	   * this reduces to two ESL_MAX calls
	   */
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	  dpn     = dn - hdmin[v][jp_v];
	  dpx     = dx - hdmin[v][jp_v];
	  dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	  	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    if(do_J) Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	    if(do_L) Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], Lalpha[y][jp_y_sdr][dp_y_sd] + tsc);
	    if(do_R) Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y_sdr][dp_y_sd] + tsc);
	  }
	}
      }
      /* no emission score to add */
    }
    else { /* B_st */ 
      y = cm->cfirst[v]; /* left  subtree */
      z = cm->cnum[v];   /* right subtree */
      
      /* Any valid j must be within both state v and state z's j band 
       * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
       * constructed, but we'll check anyway. 
       */
      jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
      jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
      /* the main j loop */
      for (j = jn; j <= jx; j++) { 
	jp_v = j - jmin[v];
	jp_y = j - jmin[y];
	jp_z = j - jmin[z];
	kn = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
	/* kn satisfies inequalities (1) and (3) (listed below)*/	
	kx = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
	/* kn satisfies inequalities (2) and (4) (listed below)*/	
	i = j - hdmin[v][jp_v] + 1;
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++, i--) {
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	  /* Find the first k value that implies a valid cell in the {J,L,R} matrix y and z decks.
	   * This k must satisfy the following 6 inequalities (some may be redundant):
	   * (1) k >= j-jmax[y];
	   * (2) k <= j-jmin[y]; 
	   *     1 and 2 guarantee (j-k) is within state y's j band
	   *
	   * (3) k >= hdmin[z][j-jmin[z]];
	   * (4) k <= hdmax[z][j-jmin[z]]; 
	   *     3 and 4 guarantee k is within z's j=(j), d band
	   *
	   * (5) k >= d-hdmax[y][j-jmin[y]-k];
	   * (6) k <= d-hdmin[y][j-jmin[y]-k]; 
	   *     5 and 6 guarantee (d-k) is within state y's j=(j-k) d band
	   *
	   * kn and kx were set above (outside (for (dp_v...) loop) that
	   * satisfy 1-4 (b/c 1-4 are d-independent and k-independent)
	   * RHS of inequalities 5 and 6 are dependent on k, so we check
	   * for these within the next for loop.
	   *
	   * To update a cell in the T matrix with a sum of an R matrix value for y
	   * and a L matrix value for z, there are 2 additional inequalities to satisfy:
	   * (7) k != i-1  (where i = j-d+1)
	   * (8) k != j
	   * These are checked for in the loop below as well. 
	   */
	  for(k = kn; k <= kx; k++) { 
	    if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k]) {
	      /* for current k, all 6 inequalities have been satisified 
	       * so we know the cells corresponding to the platonic 
	       * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
	       * alpha[z][j][k] are all within the bands. These
	       * cells correspond to alpha[v][jp_v][dp_v], 
	       * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
	       * and alpha[z][jp_z][k-hdmin[jp_z]];
	       */
	      kp_z = k-hdmin[z][jp_z];
	      dp_y = d-hdmin[y][jp_y-k];
	      if(do_J)Jalpha[v][jp_v][dp_v] = ESL_MAX(Jalpha[v][jp_v][dp_v], Jalpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]);
	      if(do_L) Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], Jalpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]);
	      if(do_R) Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], Ralpha[y][jp_y-k][dp_y - k] + Jalpha[z][jp_z][kp_z]);
	      /*if((k != i-1) && (k != j)) {*/
	      if(do_T) Talpha[v][jp_v][dp_v] = ESL_MAX(Talpha[v][jp_v][dp_v], Ralpha[y][jp_y-k][dp_y - k] + Lalpha[z][jp_z][kp_z]);
		/*}*/
	    }
	  }
	}
      }

      /* two additional special cases in trCYK (these are not in standard CYK).
       * we do these in their own for(j.. { for(d.. { } } loops b/c one 
       * is independent of z, the other of y, unlike the above loop which is dependent 
       * on both.
       */
      if(do_L) { 
	jn = (jmin[v] > jmin[y]) ? jmin[v] : jmin[y];
	jx = (jmax[v] < jmax[y]) ? jmax[v] : jmax[y];
	for (j = jn; j <= jx; j++) { 
	  jp_v = j - jmin[v];
	  jp_y = j - jmin[y];
	  ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	  ESL_DASSERT1((j >= jmin[y] && j <= jmax[y]));
	  dn = (hdmin[v][jp_v] > hdmin[y][jp_y]) ? hdmin[v][jp_v] : hdmin[y][jp_y];
	  dx = (hdmax[v][jp_v] < hdmax[y][jp_y]) ? hdmax[v][jp_v] : hdmax[y][jp_y];
	  for(d = dn; d <= dx; d++) { 
	    dp_v = d - hdmin[v][jp_v];
	    dp_y = d - hdmin[y][jp_y];
	    ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	    ESL_DASSERT1((d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y]));
	    Lalpha[v][jp_v][dp_v] = ESL_MAX(Lalpha[v][jp_v][dp_v], 
					    ESL_MAX(Jalpha[y][jp_y][dp_y],
						    Lalpha[y][jp_y][dp_y]));
	  }
	}
      }
      if(do_R) { 
	jn = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	jx = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];
	for (j = jn; j <= jx; j++) { 
	  jp_v = j - jmin[v];
	  jp_z = j - jmin[z];
	  ESL_DASSERT1((j >= jmin[v] && j <= jmax[v]));
	  ESL_DASSERT1((j >= jmin[z] && j <= jmax[z]));
	  dn = (hdmin[v][jp_v] > hdmin[z][jp_z]) ? hdmin[v][jp_v] : hdmin[z][jp_z];
	  dx = (hdmax[v][jp_v] < hdmax[z][jp_z]) ? hdmax[v][jp_v] : hdmax[z][jp_z];
	  for(d = dn; d <= dx; d++) { 
	    dp_v = d - hdmin[v][jp_v];
	    dp_z = d - hdmin[z][jp_z];
	    ESL_DASSERT1((d >= hdmin[v][jp_v] && d <= hdmax[v][jp_v]));
	    ESL_DASSERT1((d >= hdmin[z][jp_z] && d <= hdmax[z][jp_z]));
	    Ralpha[v][jp_v][dp_v] = ESL_MAX(Ralpha[v][jp_v][dp_v], 
					    ESL_MAX(Jalpha[z][jp_z][dp_z],
						    Ralpha[z][jp_z][dp_z]));
	  }
	}
      }
    } /* finished calculating deck v. */
#if PRINTFALPHA
	  if(cm->stid[v] == BIF_B) { 
	    /* the main j loop */
	    for (j = jmin[v]; j <= jmax[v]; j++) { 
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		printf("H j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		       j, v, d, 
		       NOT_IMPOSSIBLE(Jalpha[v][jp_v][dp_v]) ? Jalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Lalpha[v][jp_v][dp_v]) ? Lalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Ralpha[v][jp_v][dp_v]) ? Ralpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Talpha[v][jp_v][dp_v]) ? Talpha[v][jp_v][dp_v] : -9999.9);
	      }
	    }
	  }
	  if((cm->stid[v] == BEGL_S) || (cm->stid[v] == BEGR_S)) { 
	    /* the main j loop */
	    for (j = jmin[v]; j <= jmax[v]; j++) { 
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		printf("H j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		       j, v, d, 
		       NOT_IMPOSSIBLE(Jalpha[v][jp_v][dp_v]) ? Jalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Lalpha[v][jp_v][dp_v]) ? Lalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Ralpha[v][jp_v][dp_v]) ? Ralpha[v][jp_v][dp_v] : -9999.9);
	      }
	    }
	  }
#endif
#if 0 
	  else { 
	    for (j = jmin[v]; j <= jmax[v]; j++) { 
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		printf("H j: %3d  v: %3d  d: %3d   J: %10.4f  L: %10.4f  R: %10.4f  T: %10.4f\n", 
		       j, v, d, 
		       NOT_IMPOSSIBLE(Jalpha[v][jp_v][dp_v]) ? Jalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Lalpha[v][jp_v][dp_v]) ? Lalpha[v][jp_v][dp_v] : -9999.9,
		       NOT_IMPOSSIBLE(Ralpha[v][jp_v][dp_v]) ? Ralpha[v][jp_v][dp_v] : -9999.9, 
		       -9999.9);
	      }
	    }
	  }
	  printf("\n");
#endif
  } /* end of for (v = cm->M-1; v > 0; v--) */
        
  /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
   * 
   * If local begins are off, all hits must be rooted at v=0.
   * With local begins on, the hit is rooted at the second state in
   * the traceback (e.g. after 0), the internal entry point. 
   * 
   * Hits rooted at 0 that not involved with local begins are 
   * already calc'ed from the v loop with v == 0 
   */

  /* Report all possible hits, but only after looking at local begins (if they're on) */
  v = 0;
  sd = sdr = 0;
  jpn = 0;
  jpx = jmax[v] - jmin[v];
  j   = jmin[v];
  
  ESL_ALLOC(bestr, sizeof(int) * (W+1));
  ESL_ALLOC(bestmode, sizeof(int) * (W+1));
  esl_vec_ISet(bestr,    W+1, 0);
  esl_vec_ISet(bestmode, W+1, CM_HIT_MODE_J);

  /* first report all hits with j < jmin[0] are impossible, only if we're reporting hits to hitlist */
  if(hitlist != NULL) { 
    for(j = i0; j < jmin[v]; j++) {
      if((status = UpdateGammaHitMx(cm, errbuf, gamma, j-i0+1, 
				    NULL, 0, 0,   /* alpha_row is NULL, we're telling UpdateGammaHitMx, this j can't be a hit end pt */
				    TRUE, bestr, bestmode, hitlist, W, act)) != eslOK) return status;
    }
  }

  /* consider local hits */
  v = 0;
  for (j = jmin[v]; j <= jmax[v]; j++) {
    jp_v = j - jmin[v];
    esl_vec_ISet(bestr, (W+1), 0); /* init bestr to 0, all hits are rooted at 0 unless we find a better local begin below */
    if (cm->flags & CMH_LOCAL_BEGIN) {
      for (y = 1; y < cm->M; y++) {
	if(NOT_IMPOSSIBLE(cm->beginsc[y]) && /* local begin is possible into y */
	   (j >= jmin[y] && j <= jmax[y])) { /* j is within state y's j band */
	  jp_y = j - jmin[y];
	  dn   = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y]);
	  dx   = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y]);
	  dp_v = dn - hdmin[v][jp_v];
	  dp_y = dn - hdmin[y][jp_y];
	  for(d = dn; d <= dx; d++, dp_v++, dp_y++) { 
	    sc = Jalpha[y][jp_y][dp_y] + cm->beginsc[y];
	    if(sc > Jalpha[0][jp_v][dp_v]) {
	      Jalpha[0][jp_v][dp_v] = sc;
	      bestr[d] = y;
	    }
	  }
	}
      } /* end of for(y = 1; y < cm->M; y++) */
    } /* end of if(cm->flags & CMH_LOCAL_BEGIN */
  }    

  esl_stopwatch_Start(w);
  /* Finally, allow for truncated hits */
  v = 0;
  for (j = jmin[v]; j <= jmax[v]; j++) {
    jp_v = j - jmin[v];
    esl_vec_ISet(bestmode, (W+1), CM_HIT_MODE_J); /* init bestmode to CM_HIT_MODE_J */
    for (y = 1; y < cm->M; y++) {
      if(j >= jmin[y] && j <= jmax[y]) { /* j is within state y's band */
	jp_y = j - jmin[y];
	if(cm->sttype[y] == MP_st || cm->sttype[y] == ML_st || cm->sttype[y] == MR_st) { 
	  /*if(NOT_IMPOSSIBLE(cm->beginsc[y])) {*/
	  /*trunc_penalty = cm->beginsc[y];*/
	  dn   = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y]);
	  dx   = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y]);
	  dp_v = dn - hdmin[v][jp_v];
	  dp_y = dn - hdmin[y][jp_y];
	  for(d = dn; d <= dx; d++, dp_v++, dp_y++) { 
	    Jsc = Jalpha[y][jp_y][dp_y] + trunc_penalty; 
	    Lsc = Lalpha[y][jp_y][dp_y] + trunc_penalty; 
	    Rsc = Ralpha[y][jp_y][dp_y] + trunc_penalty; 
	    if(Jsc > Jalpha[0][jp_v][dp_v]) { 
	      Jalpha[0][jp_v][dp_v] = Jsc;
	      bestr[d] = y;
	      /* bestmode[d] remains CM_HIT_MODE_J */
	    }
	    if(Lsc > Jalpha[0][jp_v][dp_v] && (cm->sttype[y] != MR_st)) { 
	      Jalpha[0][jp_v][dp_v] = Lsc;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_L;
	    }
	    if(Rsc > Jalpha[0][jp_v][dp_v] && (cm->sttype[y] != ML_st)) { 
	      Jalpha[0][jp_v][dp_v] = Rsc;
	      bestr[d]    = y;
	      bestmode[d] = CM_HIT_MODE_R;;
	    }
	  }
	}
	else if(cm->sttype[y] == B_st) { 
	  dn   = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y]);
	  dx   = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y]);
	  dp_v = dn - hdmin[v][jp_v];
	  dp_y = dn - hdmin[y][jp_y];
	  for(d = dn; d <= dx; d++, dp_v++, dp_y++) { 
	    Tsc = Talpha[y][jp_y][dp_y] + trunc_penalty; 
	    if(Tsc > Jalpha[0][jp_v][dp_v]) {
	      Jalpha[0][jp_v][dp_v] = Tsc;
	      bestr[d] = y;
	      bestmode[d] = CM_HIT_MODE_T; 
	    }
	  }
	}
      }
    }
  }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " Considering truncated hits time: ");

    
  /* report all hits with valid d for this j, only if hitlist != NULL */
  if(hitlist != NULL) { 
    if((status = UpdateGammaHitMx(cm, errbuf, gamma, j-i0+1, Jalpha[0][jp_v], hdmin[0][jp_v], hdmax[0][jp_v], TRUE, bestr, bestmode, hitlist, W, act)) != eslOK) return status;
  }
  /* finally report all hits with j > jmax[0] are impossible, only if we're reporting hits to hitlist */
  if(hitlist != NULL) { 
    for(j = jmax[v]+1; j <= j0; j++) {
      if((status = UpdateGammaHitMx(cm, errbuf, gamma, j-i0+1, 
				    NULL, 0, 0,   /* alpha_row is NULL, we're telling UpdateGammaHitMx, this j is can't be a hit end pt */
				    TRUE, bestr, bestmode, hitlist, W, act)) != eslOK) return status;
    }
  }
  /* find the best scoring hit, and update envelope boundaries if nec */
  vsc_root = IMPOSSIBLE;
  v = 0;
  jpn = 0;
  jpx = jmax[v] - jmin[v];
  for(jp_v = jpn; jp_v <= jpx; jp_v++) {
    dpn     = 0;
    dpx     = hdmax[v][jp_v] - hdmin[v][jp_v];
    for(dp_v = dpn; dp_v <= dpx; dp_v++) {
      vsc_root = ESL_MAX(vsc_root, Jalpha[0][jp_v][dp_v]);
    }
    /* update envelope boundaries, if nec */
    if(do_env_defn) { 
      j = jp_v + jmin[v];
      for(dp_v = dpn; dp_v <= dpx; dp_v++) {
	if(Jalpha[0][jp_v][dp_v] >= env_cutoff) { 
	  i = j - (dp_v + hdmin[v][jp_v]) + 1;
	  envi = ESL_MIN(envi, i);
	  envj = ESL_MAX(envj, j);
	}
      }
    }
  }
  /* find the best score in any matrix at any state */
  float best_tr_sc = IMPOSSIBLE;
  best_tr_sc = vsc_root;
  printf("Best truncated score: %.4f (%.4f)\n",
	 best_tr_sc - trunc_penalty, 
	 best_tr_sc + sreLOG2(2./(cm->clen * (cm->clen+1))));

  free(el_scA);
  free(yvalidA);
  free(bestr);
  free(bestmode);
  if (act != NULL) { 
    for(i = 0; i <= W; i++) free(act[i]); 
    free(act);
  }

  if(hitlist != NULL && gamma->iamgreedy == FALSE) TBackGammaHitMx(gamma, hitlist, i0, j0);

  /* set envelope return variables if nec */
  if(ret_envi != NULL) { *ret_envi = (envi == j0+1) ? -1 : envi; }
  if(ret_envj != NULL) { *ret_envj = (envj == i0-1) ? -1 : envj; }

  if(gamma != NULL) FreeGammaHitMx(gamma);

  if (ret_sc != NULL) *ret_sc = vsc_root;
  ESL_DPRINTF1(("FastCYKScanHB() return sc: %f\n", vsc_root));
  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return 0.; /* never reached */
}


/*****************************************************************
 *   1. TrScanMatrix_t data structure functions,
 *      auxiliary info and matrix of float and/or int scores for 
 *      truncated query dependent banded or non-banded CM DP search 
 *      functions.
 *****************************************************************/

/* Function: cm_CreateTrScanMatrix()
 * Date:     EPN, Tue Aug 16 04:23:41 2011
 *
 * Purpose:  Given relevant info, allocate and initialize
 *           TrScanMatrix_t object.  Note that unlike a ScanMatrix_t,
 *           dmin is not used to set minimum values, even if we're
 *           going to use QDBs, because minimum subtree lengths are
 *           illogical with the truncated version of CYK/Inside, but
 *           maximum lengths are not, so <dmax> is considered here.
 *            
 * Returns:  eslOK on success, dies immediately on some error
 */
TrScanMatrix_t *
cm_CreateTrScanMatrix(CM_t *cm, int W, int *dmax, double beta_W, double beta_qdb, int do_banded, int do_float, int do_int)
{ 
  int status;
  TrScanMatrix_t *trsmx;
  int v,j;

  if((!do_float) && (!do_int)) cm_Fail("cm_CreateScanMatrix(), do_float and do_int both FALSE.\n");
  if(do_banded && dmax == NULL) cm_Fail("cm_CreateScanMatrix(), do_banded is TRUE, but dmax is NULL.\n");

  ESL_ALLOC(trsmx, sizeof(TrScanMatrix_t));

  trsmx->flags    = 0;
  trsmx->cm_M     = cm->M;
  trsmx->W        = W;
  trsmx->dmax     = dmax; /* could be NULL */
  trsmx->beta_W   = beta_W; 
  trsmx->beta_qdb = beta_qdb; 

  /* precalculate minimum and maximum d for each state and each sequence index (1..j..W). 
   * this is not always just 0, dmax, (for ex. if j < W). */
  ESL_ALLOC(trsmx->dnAA, sizeof(int *) * (trsmx->W+1));
  ESL_ALLOC(trsmx->dxAA, sizeof(int *) * (trsmx->W+1));
  trsmx->dnAA[0] = trsmx->dxAA[0] = NULL; /* corresponds to j == 0, which is out of bounds */
  for(j = 1; j <= trsmx->W; j++) {
    ESL_ALLOC(trsmx->dnAA[j], sizeof(int) * cm->M);
    ESL_ALLOC(trsmx->dxAA[j], sizeof(int) * cm->M);
    for(v = 0; v < cm->M; v++) {
      /* dnAA[j][v] is 1 for all states, even MATP, b/c d == 1 is valid for MATP in L,R matrices */
      trsmx->dnAA[j][v] = 1;
      if(do_banded) { 
	trsmx->dxAA[j][v] = ESL_MIN(j, trsmx->dmax[v]); 
	trsmx->dxAA[j][v] = ESL_MIN(trsmx->dxAA[j][v], trsmx->W);
      }
      else { 
	trsmx->dxAA[j][v] = ESL_MIN(j, trsmx->W); 
      }
    }
  }
  /* allocate bestr, which holds best root state at alpha[0][cur][d] */
  ESL_ALLOC(trsmx->bestr,    (sizeof(int) * (trsmx->W+1)));
  ESL_ALLOC(trsmx->bestmode, (sizeof(int) * (trsmx->W+1)));
  /* initialize bestmode to J for all positions */
  esl_vec_ISet(trsmx->bestmode, (trsmx->W+1), CM_HIT_MODE_J);

  /* Some info about the falpha/ialpha matrix
   * The alpha matrix holds data for all states EXCEPT BEGL_S states
   * The alpha scanning matrix is indexed [j][v][d]. 
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *    v ranges from 0..M-1 over states in the model.
   *    d ranges from 0..W over subsequence lengths.
   * Note if v is a BEGL_S alpha[j][v] == NULL
   * Note that old convention of sharing E memory is no longer,
   * each E state has it's own deck.
   *
   * alpha_begl matrix holds data for ONLY BEGL_S states
   *    j takes value of 0..W
   *    v ranges from 0..M-1 over states in the model
   *    d ranges from 0..W over subsequence lengths.
   * Note if v is NOT a BEGL_S alpha_begl[j][v] == NULL
   *
   * alpha and alpha_begl are allocated in contiguous blocks
   * of memory in {f,i}alpha_mem and {f,i}alpha_begl_mem
   */

  /* Some info on alpha initialization 
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   *
   * The way we initialize is just to set the entire matrix
   * to -INFTY or IMPOSSIBLE (for ints and floats, respectively),
   * and then reset those cells that should not be -INFTY or
   * IMPOSSIBLE as listed above. This way we don't have to
   * step through the bands, setting cells outside them to IMPOSSIBLE
   * or -INFY;
   */

  trsmx->fJalpha          = trsmx->fLalpha          = trsmx->fRalpha          = trsmx->fTalpha          = NULL;
  trsmx->fJalpha_begl     = trsmx->fLalpha_begl     = trsmx->fRalpha_begl     = NULL;
  trsmx->fJalpha_mem      = trsmx->fLalpha_mem      = trsmx->fRalpha_mem      = trsmx->fTalpha_mem      = NULL;
  trsmx->fJalpha_begl_mem = trsmx->fLalpha_begl_mem = trsmx->fRalpha_begl_mem = NULL;

  trsmx->iJalpha          = trsmx->iLalpha          = trsmx->iRalpha          = trsmx->iTalpha          = NULL;
  trsmx->iJalpha_begl     = trsmx->iLalpha_begl     = trsmx->iRalpha_begl     = NULL;
  trsmx->iJalpha_mem      = trsmx->iLalpha_mem      = trsmx->iRalpha_mem      = trsmx->iTalpha_mem      = NULL;
  trsmx->iJalpha_begl_mem = trsmx->iLalpha_begl_mem = trsmx->iRalpha_begl_mem = NULL;

  trsmx->ncells_alpha      = 0;
  trsmx->ncells_alpha_begl = 0;
  trsmx->ncells_Talpha     = 0;

  if(do_float) /* allocate float mx and scores */
    cm_FloatizeTrScanMatrix(cm, trsmx);
  if(do_int)   /* allocate int mx and scores */
    cm_IntizeTrScanMatrix(cm, trsmx);
  return trsmx;

 ERROR:
  cm_Fail("memory allocation error in cm_CreateTrScanMatrix().\n");
  return NULL; /* NEVERREACHED */
}

/* Function: cm_FloatizeTrScanMatrix()
 * Date:     EPN, Tue Aug 16 04:36:29 2011
 *
 * Purpose: Allocate and initialize float data structures in a
 *           TrScanMatrix_t object for <cm>.  This initializes a
 *           scanning float DP matrix for trCYK/trInside, for details
 *           on that matrix see the notes by the
 *           cm_TrFloatizeScanMatrix() function call in
 *           cm_CreateTrScanMatrix().
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FloatizeTrScanMatrix(CM_t *cm, TrScanMatrix_t *trsmx)
{
  int status;
  int j, v;
  int y, yoffset, w;
  int n_begl, n_bif;
  int n_non_begl;
  int cur_cell;

  /* contract check */
  if(trsmx->flags & cmTRSMX_HAS_FLOAT) cm_Fail("cm_FloatizeScanMatrix(), trsmx's cmTRSMX_HAS_FLOAT flag is already up.");
  if(trsmx->fJalpha != NULL)       cm_Fail("cm_FloatizeScanMatrix(), trsmx->fJalpha is not NULL.");
  if(trsmx->fJalpha_begl != NULL)  cm_Fail("cm_FloatizeScanMatrix(), trsmx->fJalpha_begl is not NULL.");
  if(trsmx->fLalpha != NULL)       cm_Fail("cm_FloatizeScanMatrix(), trsmx->fLalpha is not NULL.");
  if(trsmx->fLalpha_begl != NULL)  cm_Fail("cm_FloatizeScanMatrix(), trsmx->fLalpha_begl is not NULL.");
  if(trsmx->fRalpha != NULL)       cm_Fail("cm_FloatizeScanMatrix(), trsmx->fRalpha is not NULL.");
  if(trsmx->fRalpha_begl != NULL)  cm_Fail("cm_FloatizeScanMatrix(), trsmx->fRalpha_begl is not NULL.");
  if(trsmx->fTalpha != NULL)       cm_Fail("cm_FloatizeScanMatrix(), trsmx->fTalpha is not NULL.");
  
  /* allocate alpha 
   * we allocate only as many cells as necessary,
   * for f{J,L,R,T}alpha,      we only allocate for non-BEGL_S states,
   * for f{J,L,R,T}alpha_begl, we only allocate for     BEGL_S states
   * for fTalpha,              we only allocate for     BIF_B  states
   *
   * note: deck for the EL state, cm->M is never used for scanners
   */
  n_begl = 0;
  n_bif  = 0;
  for (v = 0; v < cm->M; v++) if (cm->stid[v] == BEGL_S) n_begl++;
  for (v = 0; v < cm->M; v++) if (cm->stid[v] == BIF_B)  n_bif++;
  n_non_begl = cm->M - n_begl;

  /* allocate f{J,L,R,T}alpha */
  /* j == 0 v == 0 cells, followed by j == 1 v == 0, then j == 0 v == 1 etc.. */
  ESL_ALLOC(trsmx->fJalpha,        sizeof(float **) * 2);
  ESL_ALLOC(trsmx->fJalpha[0],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fJalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->fJalpha[1],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fJalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->fJalpha_mem,    sizeof(float) * 2 * n_non_begl * (trsmx->W+1));

  ESL_ALLOC(trsmx->fLalpha,        sizeof(float **) * 2);
  ESL_ALLOC(trsmx->fLalpha[0],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fLalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->fLalpha[1],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fLalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->fLalpha_mem,    sizeof(float) * 2 * n_non_begl * (trsmx->W+1));

  ESL_ALLOC(trsmx->fRalpha,        sizeof(float **) * 2);
  ESL_ALLOC(trsmx->fRalpha[0],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fRalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->fRalpha[1],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fRalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->fRalpha_mem,    sizeof(float) * 2 * n_non_begl * (trsmx->W+1));

  ESL_ALLOC(trsmx->fTalpha,        sizeof(float **) * 2);
  ESL_ALLOC(trsmx->fTalpha[0],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BIF_B, fTalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->fTalpha[1],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BIF_B, fTalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->fTalpha_mem,    sizeof(float) * 2 * n_non_begl * (trsmx->W+1));

  if((trsmx->flags & cmTRSMX_HAS_INT) && ((2 * n_non_begl * (trsmx->W+1)) != trsmx->ncells_alpha)) 
    cm_Fail("cm_FloatizeScanMatrix(), cmTRSMX_HAS_INT flag raised, but trsmx->ncells_alpha %d != %d (predicted num float cells size)\n", trsmx->ncells_alpha, (2 * n_non_begl * (trsmx->W+1)));
  trsmx->ncells_alpha = 2 * n_non_begl * (trsmx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    if (cm->stid[v] != BEGL_S) {
      trsmx->fJalpha[0][v] = trsmx->fJalpha_mem + cur_cell;
      trsmx->fLalpha[0][v] = trsmx->fLalpha_mem + cur_cell;
      trsmx->fRalpha[0][v] = trsmx->fRalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
      trsmx->fJalpha[1][v] = trsmx->fJalpha_mem + cur_cell;
      trsmx->fLalpha[1][v] = trsmx->fLalpha_mem + cur_cell;
      trsmx->fRalpha[1][v] = trsmx->fRalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
    }
    else { 
      trsmx->fJalpha[0][v] = NULL;
      trsmx->fJalpha[1][v] = NULL;
      trsmx->fLalpha[0][v] = NULL;
      trsmx->fLalpha[1][v] = NULL;
      trsmx->fRalpha[0][v] = NULL;
      trsmx->fRalpha[1][v] = NULL;
    }
  }
  if(cur_cell != trsmx->ncells_alpha) cm_Fail("cm_FloatizeScanMatrix(), error allocating falpha, cell cts differ %d != %d\n", cur_cell, trsmx->ncells_alpha);

  if((trsmx->flags & cmTRSMX_HAS_INT) && ((2 * n_bif * (trsmx->W+1)) != trsmx->ncells_Talpha)) 
    cm_Fail("cm_FloatizeScanMatrix(), cmTRSMX_HAS_INT flag raised, but trsmx->ncells_Talpha %d != %d (predicted num float cells size in Talpha)\n", trsmx->ncells_Talpha, (2 * n_bif * (trsmx->W+1)));
  trsmx->ncells_Talpha = 2 * n_bif * (trsmx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    if (cm->stid[v] == BIF_B) { 
      trsmx->fTalpha[0][v] = trsmx->fTalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
      trsmx->fTalpha[1][v] = trsmx->fTalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
    }
    else { 
      trsmx->fTalpha[0][v] = NULL;
      trsmx->fTalpha[1][v] = NULL;
    }
  }
  if(cur_cell != trsmx->ncells_Talpha) cm_Fail("cm_FloatizeScanMatrix(), error allocating fTalpha, cell cts differ %d != %d\n", cur_cell, trsmx->ncells_Talpha);
  

  /* allocate falpha_begl */
  /* j == d, v == 0 cells, followed by j == d+1, v == 0, etc. */
  ESL_ALLOC(trsmx->fJalpha_begl, sizeof(float **) * (trsmx->W+1));
  ESL_ALLOC(trsmx->fLalpha_begl, sizeof(float **) * (trsmx->W+1));
  ESL_ALLOC(trsmx->fRalpha_begl, sizeof(float **) * (trsmx->W+1));
  for (j = 0; j <= trsmx->W; j++) {
    ESL_ALLOC(trsmx->fJalpha_begl[j],  sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, fJalpha_begl[0][v] will be NULL */
    ESL_ALLOC(trsmx->fLalpha_begl[j],  sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, fLalpha_begl[0][v] will be NULL */
    ESL_ALLOC(trsmx->fRalpha_begl[j],  sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, fRalpha_begl[0][v] will be NULL */
  }
  ESL_ALLOC(trsmx->fJalpha_begl_mem,   sizeof(float) * (trsmx->W+1) * n_begl * (trsmx->W+1));
  ESL_ALLOC(trsmx->fLalpha_begl_mem,   sizeof(float) * (trsmx->W+1) * n_begl * (trsmx->W+1));
  ESL_ALLOC(trsmx->fRalpha_begl_mem,   sizeof(float) * (trsmx->W+1) * n_begl * (trsmx->W+1));
  if((trsmx->flags & cmTRSMX_HAS_INT) && (((trsmx->W+1) * n_begl * (trsmx->W+1)) != trsmx->ncells_alpha_begl)) 
    cm_Fail("cm_IntizeScanMatrix(), cmTRSMX_HAS_FLOAT flag raised, but trsmx->ncells_alpha_begl %d != %d (predicted num float cells size)\n", trsmx->ncells_alpha_begl, ((trsmx->W+1) * n_begl * (trsmx->W+1)));
  trsmx->ncells_alpha_begl = (trsmx->W+1) * n_begl * (trsmx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    for (j = 0; j <= trsmx->W; j++) { 
      if (cm->stid[v] == BEGL_S) {
	trsmx->fJalpha_begl[j][v] = trsmx->fJalpha_begl_mem + cur_cell;
	trsmx->fLalpha_begl[j][v] = trsmx->fLalpha_begl_mem + cur_cell;
	trsmx->fRalpha_begl[j][v] = trsmx->fRalpha_begl_mem + cur_cell;
	cur_cell += trsmx->W+1;
      }
      else { 
	trsmx->fJalpha_begl[j][v] = NULL;
	trsmx->fLalpha_begl[j][v] = NULL;
	trsmx->fRalpha_begl[j][v] = NULL;
      }
    }
  }
  if(cur_cell != trsmx->ncells_alpha_begl) cm_Fail("cm_FloatizeScanMatrix(), error allocating falpha_begl, cell cts differ %d != %d\n", cur_cell, trsmx->ncells_alpha_begl);

  /* Initialize matrix */
  /* First, init entire matrix to IMPOSSIBLE */
  esl_vec_FSet(trsmx->fJalpha_mem,      trsmx->ncells_alpha,      IMPOSSIBLE);
  esl_vec_FSet(trsmx->fLalpha_mem,      trsmx->ncells_alpha,      IMPOSSIBLE);
  esl_vec_FSet(trsmx->fRalpha_mem,      trsmx->ncells_alpha,      IMPOSSIBLE);
  esl_vec_FSet(trsmx->fTalpha_mem,      trsmx->ncells_Talpha,     IMPOSSIBLE);
  esl_vec_FSet(trsmx->fJalpha_begl_mem, trsmx->ncells_alpha_begl, IMPOSSIBLE);
  esl_vec_FSet(trsmx->fLalpha_begl_mem, trsmx->ncells_alpha_begl, IMPOSSIBLE);
  esl_vec_FSet(trsmx->fRalpha_begl_mem, trsmx->ncells_alpha_begl, IMPOSSIBLE);
  /* Now, initialize cells that should not be IMPOSSIBLE in f{J,L,RT}alpha and f{J,L,R}alpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      if (cm->sttype[v] == E_st) { 
	trsmx->fJalpha[0][v][0] = trsmx->fJalpha[1][v][0] = 0.;
	trsmx->fLalpha[0][v][0] = trsmx->fLalpha[1][v][0] = 0.;
	trsmx->fRalpha[0][v][0] = trsmx->fRalpha[1][v][0] = 0.;
	/* rest of E deck is IMPOSSIBLE, it's already set */
      }
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	trsmx->fJalpha[0][v][0] = cm->endsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  trsmx->fJalpha[0][v][0] = ESL_MAX(trsmx->fJalpha[0][v][0], (trsmx->fJalpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	trsmx->fJalpha[0][v][0] = ESL_MAX(trsmx->fJalpha[0][v][0], IMPOSSIBLE);
	/* {L,R}alpha[0][v][0] remain IMPOSSIBLE */
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	trsmx->fJalpha[0][v][0] = trsmx->fJalpha_begl[0][w][0] + trsmx->fJalpha[0][y][0]; 
      }
      trsmx->fJalpha[1][v][0] = trsmx->fJalpha[0][v][0];
      /* {L,R,T}alpha[{0,1}][v][0] remain IMPOSSIBLE */
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      trsmx->fJalpha_begl[0][v][0] = cm->endsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	trsmx->fJalpha_begl[0][v][0] = ESL_MAX(trsmx->fJalpha_begl[0][v][0], (trsmx->fJalpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in trsmx->fJalpha */
      trsmx->fJalpha_begl[0][v][0] = ESL_MAX(trsmx->fJalpha_begl[0][v][0], IMPOSSIBLE);
      for (j = 1; j <= trsmx->W; j++) 
	trsmx->fJalpha_begl[j][v][0] = trsmx->fJalpha_begl[0][v][0];
      /* {L,R}alpha_begl[j][v][0] remain IMPOSSIBLE for all j */
    }
  }
  /* set the flag that tells us we've got valid floats */
  trsmx->flags |= cmTRSMX_HAS_FLOAT;
  return eslOK;

 ERROR: 
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: cm_IntizeTrScanMatrix()
 * Date:     EPN, Wed Aug 24 15:00:32 2011
 *
 * Purpose: Allocate and initialize integer data structures in a
 *           TrScanMatrix_t object for <cm>.  This initializes a
 *           scanning int DP matrix for trCYK/trInside, for details
 *           on that matrix see the notes by the
 *           cm_TrIntizeScanMatrix() function call in
 *           cm_CreateTrScanMatrix().
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_IntizeTrScanMatrix(CM_t *cm, TrScanMatrix_t *trsmx)
{
  int status;
  int j, v;
  int y, yoffset, w;
  int n_begl, n_bif;
  int n_non_begl;
  int cur_cell;

  /* contract check */
  if(trsmx->flags & cmTRSMX_HAS_INT) cm_Fail("cm_IntizeScanMatrix(), trsmx's cmTRSMX_HAS_INT flag is already up.");
  if(trsmx->iJalpha != NULL)         cm_Fail("cm_IntizeScanMatrix(), trsmx->iJalpha is not NULL.");
  if(trsmx->iJalpha_begl != NULL)    cm_Fail("cm_IntizeScanMatrix(), trsmx->iJalpha_begl is not NULL.");
  if(trsmx->iLalpha != NULL)         cm_Fail("cm_IntizeScanMatrix(), trsmx->iLalpha is not NULL.");
  if(trsmx->iLalpha_begl != NULL)    cm_Fail("cm_IntizeScanMatrix(), trsmx->iLalpha_begl is not NULL.");
  if(trsmx->iRalpha != NULL)         cm_Fail("cm_IntizeScanMatrix(), trsmx->iRalpha is not NULL.");
  if(trsmx->iRalpha_begl != NULL)    cm_Fail("cm_IntizeScanMatrix(), trsmx->iRalpha_begl is not NULL.");
  if(trsmx->iTalpha != NULL)         cm_Fail("cm_IntizeScanMatrix(), trsmx->iTalpha is not NULL.");
  
  /* allocate alpha 
   * we allocate only as many cells as necessary,
   * for i{J,L,R}alpha,      we only allocate for non-BEGL_S states,
   * for i{J,L,R}alpha_begl, we only allocate for     BEGL_S states
   * for iTalpha,            we only allocate for     BIF_B  states
   *
   * note: deck for the EL state, cm->M is never used for scanners
   */
  n_begl = 0;
  n_bif  = 0;
  for (v = 0; v < cm->M; v++) if (cm->stid[v] == BEGL_S) n_begl++;
  for (v = 0; v < cm->M; v++) if (cm->stid[v] == BIF_B)  n_bif++;
  n_non_begl = cm->M - n_begl;

  /* allocate f{J,L,R,T}alpha */
  /* j == 0 v == 0 cells, followed by j == 1 v == 0, then j == 0 v == 1 etc.. */
  ESL_ALLOC(trsmx->iJalpha,        sizeof(int **) * 2);
  ESL_ALLOC(trsmx->iJalpha[0],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fJalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->iJalpha[1],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fJalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->iJalpha_mem,    sizeof(int) * 2 * n_non_begl * (trsmx->W+1));

  ESL_ALLOC(trsmx->iLalpha,        sizeof(int **) * 2);
  ESL_ALLOC(trsmx->iLalpha[0],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fLalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->iLalpha[1],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fLalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->iLalpha_mem,    sizeof(int) * 2 * n_non_begl * (trsmx->W+1));

  ESL_ALLOC(trsmx->iRalpha,        sizeof(int **) * 2);
  ESL_ALLOC(trsmx->iRalpha[0],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fRalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->iRalpha[1],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, fRalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->iRalpha_mem,    sizeof(int) * 2 * n_non_begl * (trsmx->W+1));

  ESL_ALLOC(trsmx->iTalpha,        sizeof(int **) * 2);
  ESL_ALLOC(trsmx->iTalpha[0],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BIF_B, fTalpha[0][v] will be NULL */
  ESL_ALLOC(trsmx->iTalpha[1],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BIF_B, fTalpha[1][v] will be NULL */
  ESL_ALLOC(trsmx->iTalpha_mem,    sizeof(int) * 2 * n_non_begl * (trsmx->W+1));

  if((trsmx->flags & cmTRSMX_HAS_INT) && ((2 * n_non_begl * (trsmx->W+1)) != trsmx->ncells_alpha)) 
    cm_Fail("cm_IntizeScanMatrix(), cmTRSMX_HAS_INT flag raised, but trsmx->ncells_alpha %d != %d (predicted num int cells size)\n", trsmx->ncells_alpha, (2 * n_non_begl * (trsmx->W+1)));
  trsmx->ncells_alpha = 2 * n_non_begl * (trsmx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    if (cm->stid[v] != BEGL_S) {
      trsmx->iJalpha[0][v] = trsmx->iJalpha_mem + cur_cell;
      trsmx->iLalpha[0][v] = trsmx->iLalpha_mem + cur_cell;
      trsmx->iRalpha[0][v] = trsmx->iRalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
      trsmx->iJalpha[1][v] = trsmx->iJalpha_mem + cur_cell;
      trsmx->iLalpha[1][v] = trsmx->iLalpha_mem + cur_cell;
      trsmx->iRalpha[1][v] = trsmx->iRalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
    }
    else { 
      trsmx->iJalpha[0][v] = NULL;
      trsmx->iJalpha[1][v] = NULL;
      trsmx->iLalpha[0][v] = NULL;
      trsmx->iLalpha[1][v] = NULL;
      trsmx->iRalpha[0][v] = NULL;
      trsmx->iRalpha[1][v] = NULL;
    }
  }
  if(cur_cell != trsmx->ncells_alpha) cm_Fail("cm_IntizeScanMatrix(), error allocating falpha, cell cts differ %d != %d\n", cur_cell, trsmx->ncells_alpha);

  if((trsmx->flags & cmTRSMX_HAS_INT) && ((2 * n_bif * (trsmx->W+1)) != trsmx->ncells_Talpha)) 
    cm_Fail("cm_IntizeScanMatrix(), cmTRSMX_HAS_INT flag raised, but trsmx->ncells_Talpha %d != %d (predicted num int cells size in Talpha)\n", trsmx->ncells_Talpha, (2 * n_bif * (trsmx->W+1)));
  trsmx->ncells_Talpha = 2 * n_bif * (trsmx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    if (cm->stid[v] == BIF_B) { 
      trsmx->iTalpha[0][v] = trsmx->iTalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
      trsmx->iTalpha[1][v] = trsmx->iTalpha_mem + cur_cell;
      cur_cell += trsmx->W+1;
    }
    else { 
      trsmx->iTalpha[0][v] = NULL;
      trsmx->iTalpha[1][v] = NULL;
    }
  }
  if(cur_cell != trsmx->ncells_Talpha) cm_Fail("cm_IntizeScanMatrix(), error allocating iTalpha, cell cts differ %d != %d\n", cur_cell, trsmx->ncells_Talpha);
  

  /* allocate falpha_begl */
  /* j == d, v == 0 cells, followed by j == d+1, v == 0, etc. */
  ESL_ALLOC(trsmx->iJalpha_begl, sizeof(int **) * (trsmx->W+1));
  ESL_ALLOC(trsmx->iLalpha_begl, sizeof(int **) * (trsmx->W+1));
  ESL_ALLOC(trsmx->iRalpha_begl, sizeof(int **) * (trsmx->W+1));
  for (j = 0; j <= trsmx->W; j++) {
    ESL_ALLOC(trsmx->iJalpha_begl[j],  sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, fJalpha_begl[0][v] will be NULL */
    ESL_ALLOC(trsmx->iLalpha_begl[j],  sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, fLalpha_begl[0][v] will be NULL */
    ESL_ALLOC(trsmx->iRalpha_begl[j],  sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, fRalpha_begl[0][v] will be NULL */
  }
  ESL_ALLOC(trsmx->iJalpha_begl_mem,   sizeof(int) * (trsmx->W+1) * n_begl * (trsmx->W+1));
  ESL_ALLOC(trsmx->iLalpha_begl_mem,   sizeof(int) * (trsmx->W+1) * n_begl * (trsmx->W+1));
  ESL_ALLOC(trsmx->iRalpha_begl_mem,   sizeof(int) * (trsmx->W+1) * n_begl * (trsmx->W+1));
  if((trsmx->flags & cmTRSMX_HAS_INT) && (((trsmx->W+1) * n_begl * (trsmx->W+1)) != trsmx->ncells_alpha_begl)) 
    cm_Fail("cm_IntizeScanMatrix(), cmTRSMX_HAS_INT flag raised, but trsmx->ncells_alpha_begl %d != %d (predicted num int cells size)\n", trsmx->ncells_alpha_begl, ((trsmx->W+1) * n_begl * (trsmx->W+1)));
  trsmx->ncells_alpha_begl = (trsmx->W+1) * n_begl * (trsmx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    for (j = 0; j <= trsmx->W; j++) { 
      if (cm->stid[v] == BEGL_S) {
	trsmx->iJalpha_begl[j][v] = trsmx->iJalpha_begl_mem + cur_cell;
	trsmx->iLalpha_begl[j][v] = trsmx->iLalpha_begl_mem + cur_cell;
	trsmx->iRalpha_begl[j][v] = trsmx->iRalpha_begl_mem + cur_cell;
	cur_cell += trsmx->W+1;
      }
      else { 
	trsmx->iJalpha_begl[j][v] = NULL;
	trsmx->iLalpha_begl[j][v] = NULL;
	trsmx->iRalpha_begl[j][v] = NULL;
      }
    }
  }
  if(cur_cell != trsmx->ncells_alpha_begl) cm_Fail("cm_IntizeScanMatrix(), error allocating ialpha_begl, cell cts differ %d != %d\n", cur_cell, trsmx->ncells_alpha_begl);

 /* Initialize matrix */
  /* First, init entire matrix to -INFTY */
  esl_vec_ISet(trsmx->iJalpha_mem,      trsmx->ncells_alpha,      -INFTY);
  esl_vec_ISet(trsmx->iLalpha_mem,      trsmx->ncells_alpha,      -INFTY);
  esl_vec_ISet(trsmx->iRalpha_mem,      trsmx->ncells_alpha,      -INFTY);
  esl_vec_ISet(trsmx->iTalpha_mem,      trsmx->ncells_Talpha,     -INFTY);
  esl_vec_ISet(trsmx->iJalpha_begl_mem, trsmx->ncells_alpha_begl, -INFTY);
  esl_vec_ISet(trsmx->iLalpha_begl_mem, trsmx->ncells_alpha_begl, -INFTY);
  esl_vec_ISet(trsmx->iRalpha_begl_mem, trsmx->ncells_alpha_begl, -INFTY);
  /* Now, initialize cells that should not be -INFTY in f{J,L,RT}alpha and f{J,L,R}alpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      if (cm->sttype[v] == E_st) { 
	trsmx->iJalpha[0][v][0] = trsmx->iJalpha[1][v][0] = 0.;
	trsmx->iLalpha[0][v][0] = trsmx->iLalpha[1][v][0] = 0.;
	trsmx->iRalpha[0][v][0] = trsmx->iRalpha[1][v][0] = 0.;
	/* rest of E deck is -INFTY, it's already set */
      }
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	trsmx->iJalpha[0][v][0] = cm->endsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  trsmx->iJalpha[0][v][0] = ESL_MAX(trsmx->iJalpha[0][v][0], (trsmx->iJalpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	trsmx->iJalpha[0][v][0] = ESL_MAX(trsmx->iJalpha[0][v][0], -INFTY);
	/* {L,R}alpha[0][v][0] remain -INFTY */
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	trsmx->iJalpha[0][v][0] = trsmx->iJalpha_begl[0][w][0] + trsmx->iJalpha[0][y][0]; 
      }
      trsmx->iJalpha[1][v][0] = trsmx->iJalpha[0][v][0];
      /* {L,R,T}alpha[{0,1}][v][0] remain -INFTY */
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      trsmx->iJalpha_begl[0][v][0] = cm->endsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	trsmx->iJalpha_begl[0][v][0] = ESL_MAX(trsmx->iJalpha_begl[0][v][0], (trsmx->iJalpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in trsmx->iJalpha */
      trsmx->iJalpha_begl[0][v][0] = ESL_MAX(trsmx->iJalpha_begl[0][v][0], -INFTY);
      for (j = 1; j <= trsmx->W; j++) 
	trsmx->iJalpha_begl[j][v][0] = trsmx->iJalpha_begl[0][v][0];
      /* {L,R}alpha_begl[j][v][0] remain -INFTY for all j */
    }
  }
  /* set the flag that tells us we've got valid ints */
  trsmx->flags |= cmTRSMX_HAS_INT;
  return eslOK;

 ERROR: 
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: cm_FreeTrScanMatrix()
 * Date:     EPN, Wed Aug 17 14:22:45 2011
 *
 * Purpose:  Free a TrScanMatrix_t object corresponding
 *           to CM <cm>.
 *            
 * Returns:  void
 */
void
cm_FreeTrScanMatrix(CM_t *cm, TrScanMatrix_t *trsmx)
{
  int j;
  //if(! ((cm->flags & CMH_TRSCANMATRIX) && (trsmx == cm->trsmx))) { /* don't free the cm->trsmx's dmax */
  //if(trsmx->dmax != cm->dmax && trsmx->dmax != NULL) { free(trsmx->dmax); trsmx->dmax = NULL; }
  //}

  for(j = 1; j <= trsmx->W; j++) {
    free(trsmx->dnAA[j]);
    free(trsmx->dxAA[j]);
  }
  free(trsmx->dnAA);
  free(trsmx->dxAA);
  free(trsmx->bestr);
  free(trsmx->bestmode);
  
  if(trsmx->flags & cmTRSMX_HAS_FLOAT) cm_FreeFloatsFromTrScanMatrix(cm, trsmx);
  if(trsmx->flags & cmTRSMX_HAS_INT)   cm_FreeIntsFromTrScanMatrix(cm, trsmx);
  free(trsmx);
  return;
}

/* Function: cm_FreeFloatsFromTrScanMatrix()
 * Date:     EPN, Wed Aug 17 14:19:21 2011
 *
 * Purpose:  Free float data structures in a TrScanMatrix_t object 
 *           corresponding to <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeFloatsFromTrScanMatrix(CM_t *cm, TrScanMatrix_t *trsmx)
{
  int j;

  /* contract check */
  if(! trsmx->flags & cmTRSMX_HAS_FLOAT)  cm_Fail("cm_FreeFloatsFromScanMatrix(), si's cmTRSMX_HAS_FLOAT flag is down.");
  if(trsmx->fJalpha == NULL)              cm_Fail("cm_FreeFloatsFromScanMatrix(), trsmx->fJalpha is already NULL.");
  if(trsmx->fJalpha_begl == NULL)         cm_Fail("cm_FreeFloatsFromScanMatrix(), trsmx->fJalpha_begl is already NULL.");
  if(trsmx->fLalpha == NULL)              cm_Fail("cm_FreeFloatsFromScanMatrix(), trsmx->fLalpha is already NULL.");
  if(trsmx->fLalpha_begl == NULL)         cm_Fail("cm_FreeFloatsFromScanMatrix(), trsmx->fLalpha_begl is already NULL.");
  if(trsmx->fRalpha == NULL)              cm_Fail("cm_FreeFloatsFromScanMatrix(), trsmx->fRalpha is already NULL.");
  if(trsmx->fRalpha_begl == NULL)         cm_Fail("cm_FreeFloatsFromScanMatrix(), trsmx->fRalpha_begl is already NULL.");
  if(trsmx->fTalpha == NULL)              cm_Fail("cm_FreeFloatsFromScanMatrix(), trsmx->fTalpha is already NULL.");

  free(trsmx->fJalpha_mem);
  free(trsmx->fJalpha[1]);
  free(trsmx->fJalpha[0]);
  free(trsmx->fJalpha);
  trsmx->fJalpha = NULL;

  free(trsmx->fLalpha_mem);
  free(trsmx->fLalpha[1]);
  free(trsmx->fLalpha[0]);
  free(trsmx->fLalpha);
  trsmx->fLalpha = NULL;

  free(trsmx->fRalpha_mem);
  free(trsmx->fRalpha[1]);
  free(trsmx->fRalpha[0]);
  free(trsmx->fRalpha);
  trsmx->fRalpha = NULL;

  free(trsmx->fTalpha_mem);
  free(trsmx->fTalpha[1]);
  free(trsmx->fTalpha[0]);
  free(trsmx->fTalpha);
  trsmx->fTalpha = NULL;

  free(trsmx->fJalpha_begl_mem);
  for (j = 0; j <= trsmx->W; j++) free(trsmx->fJalpha_begl[j]);
  free(trsmx->fJalpha_begl);
  trsmx->fJalpha_begl = NULL;

  free(trsmx->fLalpha_begl_mem);
  for (j = 0; j <= trsmx->W; j++) free(trsmx->fLalpha_begl[j]);
  free(trsmx->fLalpha_begl);
  trsmx->fLalpha_begl = NULL;

  free(trsmx->fRalpha_begl_mem);
  for (j = 0; j <= trsmx->W; j++) free(trsmx->fRalpha_begl[j]);
  free(trsmx->fRalpha_begl);
  trsmx->fRalpha_begl = NULL;

  trsmx->flags &= ~cmTRSMX_HAS_FLOAT;
  return eslOK;
}


/* Function: cm_FreeIntsFromTrScanMatrix()
 * Date:     EPN, Wed Aug 24 14:56:18 2011
 *
 * Purpose:  Free float data structures in a TrScanMatrix_t object 
 *           corresponding to <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeIntsFromTrScanMatrix(CM_t *cm, TrScanMatrix_t *trsmx)
{
  int j;

  /* contract check */
  if(! trsmx->flags & cmTRSMX_HAS_INT)  cm_Fail("cm_FreeIntsFromScanMatrix(), si's cmTRSMX_HAS_FLOAT flag is down.");
  if(trsmx->iJalpha == NULL)            cm_Fail("cm_FreeIntsFromScanMatrix(), trsmx->iJalpha is already NULL.");
  if(trsmx->iJalpha_begl == NULL)       cm_Fail("cm_FreeIntsFromScanMatrix(), trsmx->iJalpha_begl is already NULL.");
  if(trsmx->iLalpha == NULL)            cm_Fail("cm_FreeIntsFromScanMatrix(), trsmx->iLalpha is already NULL.");
  if(trsmx->iLalpha_begl == NULL)       cm_Fail("cm_FreeIntsFromScanMatrix(), trsmx->iLalpha_begl is already NULL.");
  if(trsmx->iRalpha == NULL)            cm_Fail("cm_FreeIntsFromScanMatrix(), trsmx->iRalpha is already NULL.");
  if(trsmx->iRalpha_begl == NULL)       cm_Fail("cm_FreeIntsFromScanMatrix(), trsmx->iRalpha_begl is already NULL.");
  if(trsmx->iTalpha == NULL)            cm_Fail("cm_FreeIntsFromScanMatrix(), trsmx->iTalpha is already NULL.");

  free(trsmx->iJalpha_mem);
  free(trsmx->iJalpha[1]);
  free(trsmx->iJalpha[0]);
  free(trsmx->iJalpha);
  trsmx->iJalpha = NULL;

  free(trsmx->iLalpha_mem);
  free(trsmx->iLalpha[1]);
  free(trsmx->iLalpha[0]);
  free(trsmx->iLalpha);
  trsmx->iLalpha = NULL;

  free(trsmx->iRalpha_mem);
  free(trsmx->iRalpha[1]);
  free(trsmx->iRalpha[0]);
  free(trsmx->iRalpha);
  trsmx->iRalpha = NULL;

  free(trsmx->iTalpha_mem);
  free(trsmx->iTalpha[1]);
  free(trsmx->iTalpha[0]);
  free(trsmx->iTalpha);
  trsmx->iTalpha = NULL;

  free(trsmx->iJalpha_begl_mem);
  for (j = 0; j <= trsmx->W; j++) free(trsmx->iJalpha_begl[j]);
  free(trsmx->iJalpha_begl);
  trsmx->iJalpha_begl = NULL;

  free(trsmx->iLalpha_begl_mem);
  for (j = 0; j <= trsmx->W; j++) free(trsmx->iLalpha_begl[j]);
  free(trsmx->iLalpha_begl);
  trsmx->iLalpha_begl = NULL;

  free(trsmx->iRalpha_begl_mem);
  for (j = 0; j <= trsmx->W; j++) free(trsmx->iRalpha_begl[j]);
  free(trsmx->iRalpha_begl);
  trsmx->iRalpha_begl = NULL;

  trsmx->flags &= ~cmTRSMX_HAS_INT;
  return eslOK;
}


/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_TRUNC_SEARCH_BENCHMARK
/* Next line is optimized (debugging not on) on wyvern:
 * gcc   -o benchmark-trunc-search -std=gnu99 -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -pthread -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DIMPL_TRUNC_SEARCH_BENCHMARK cm_dpsearch_trunc.c -linfernal -lhmmer -leasel -lm
 * gcc   -o benchmark-trunc-search -std=gnu99 -g -Wall -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DIMPL_TRUNC_SEARCH_BENCHMARK cm_dpsearch_trunc.c -linfernal -lhmmer -leasel -lm
 * ./benchmark-trunc-search <cmfile>
 */

#include "esl_config.h"
#include "p7_config.h"
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
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,    "181", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>, '0' for one-time arbitrary", 0 },
  { "-e",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "emit sequences from CM, don't randomly create them", 0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search in glocal mode [default: local]", 0 },
  { "-L",        eslARG_INT,  "10000", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-T",        eslARG_REAL,    "5.", NULL, NULL,  NULL,  NULL, NULL, "set bit score reporting threshold as <x>",       0 },
  { "--orig",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also search with original trCYK",                0},
  { "--dc",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also search with D&C trCYK",                     0},
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "don't use QDBs", 0},
  { "--infile",  eslARG_INFILE,  NULL, NULL, NULL,  NULL,  NULL, "-L,-N,-e", "read sequences to search from file <s>", 2 },
  { "--i27",     eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "reproduce Kolbe, Eddy 2009 marginal score calculation", 2 },
  { "--hb",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also run HMM banded scanning trCYK", 2 },
  { "--ins",     eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also run trInside", 2 },
  { "--tau",     eslARG_REAL,   "5e-6",NULL, "0<x<1",NULL,"--hb",NULL, "set tail loss prob for --hb to <x>", 2 },
  { "--cp9noel", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-g",           "turn OFF local ends in cp9 HMMs", 2 },
  { "--cp9gloc", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-g,--cp9noel", "configure CP9 HMM in glocal mode", 2 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "benchmark driver for scanning trCYK implementations";

int 
main(int argc, char **argv)
{
  int             status;
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  CM_t           *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq;
  int             i;
  float           sc;
  char           *cmfile = esl_opt_GetArg(go, 1);
  CM_FILE        *cmfp;	/* open input CM file stream */
  int            *dmin;
  int            *dmax;
  int             do_random;
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, either randomly created, or emitted from CM (if -e) */
  char            errbuf[cmERRBUFSIZE];
  TrScanMatrix_t *trsmx = NULL;
  ESL_SQFILE     *sqfp  = NULL;        /* open sequence input file stream */
  CMConsensus_t  *cons  = NULL;
  Parsetree_t    *tr    = NULL;
  CM_TR_HB_MX    *trhbmx= NULL;
  CMEmitMap_t    *emap  = NULL;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((status = cm_file_Open(cmfile, NULL, FALSE, &(cmfp), errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = cm_file_Read(cmfp, TRUE, &abc, &cm))                != eslOK) cm_Fail(cmfp->errbuf);
  cm_file_Close(cmfp);

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      if(! esl_opt_GetBoolean(go, "--cp9noel")) cm->config_opts |= CM_CONFIG_HMMEL; 
    }
  }
  if( esl_opt_GetBoolean(go, "--noqdb")) cm->search_opts |= CM_SEARCH_NOQDB;
  cm->search_opts |= CM_SEARCH_CMGREEDY; /* greedy hit resolution is default */

  esl_stopwatch_Start(w);
  printf("%-30s", "Configuring CM...");
  fflush(stdout);
  ConfigCM(cm, errbuf, FALSE, NULL, NULL); /* FALSE says: don't calculate W */
  printf("done.  ");
  fflush(stdout);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " CPU time: ");

  if (esl_opt_GetBoolean(go, "--noqdb")) { 
    if(cm->dmin != NULL) { free(cm->dmin); cm->dmin = NULL; }
    if(cm->dmax != NULL) { free(cm->dmax); cm->dmax = NULL; }
  }
  dmin = cm->dmin; 
  dmax = cm->dmax; 
  cm->tau = esl_opt_GetReal(go, "--tau");

#if 0
  printf("%-30s", "Creating scan matrix...");
  fflush(stdout);
  esl_stopwatch_Start(w);
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE); /* impt to do this after QDBs set up in ConfigCM() */
  printf("done.  ");
  fflush(stdout);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " CPU time: ");

  printf("%-30s", "Creating tr scan matrix...");
  esl_stopwatch_Start(w);
  trsmx = cm_CreateTrScanMatrix(cm, cm->W, dmax, cm->beta_W, cm->beta_qdb, 
				(dmin == NULL && dmax == NULL) ? FALSE : TRUE,
				TRUE, TRUE); /* do_float, do_int */
  if(trsmx == NULL) esl_fatal("Problem creating trsmx");
  printf("done.  ");
  fflush(stdout);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, " CPU time: ");
#endif

  if(esl_opt_GetBoolean(go, "--hb")) { 
    printf("%-30s", "Creating tr hb matrix...");
    fflush(stdout);
    esl_stopwatch_Start(w);
    trhbmx = cm_tr_hb_mx_Create(cm);
    printf("done.  ");
    fflush(stdout);
    esl_stopwatch_Stop(w);
    esl_stopwatch_Display(stdout, w, " CPU time: ");
  }
  printf("\n\n");

  /* get sequences */
  if(esl_opt_IsUsed(go, "--infile")) { 
    /* read sequences from a file */
    status = esl_sqfile_OpenDigital(cm->abc, esl_opt_GetString(go, "--infile"), eslSQFILE_UNKNOWN, NULL, &sqfp);
    if (status == eslENOTFOUND)    esl_fatal("File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--infile"));
    else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--infile"));
    else if (status == eslEINVAL)  esl_fatal("Can't autodetect stdin or .gz."); 
    else if (status != eslOK)      esl_fatal("Sequence file open failed with error %d.\n", status);

    seqs_to_aln = CreateSeqsToAln(100, FALSE);
    if((status = ReadSeqsToAln(cm->abc, sqfp, 0, seqs_to_aln, FALSE)) != eslEOF)
      esl_fatal("Error reading sqfile: %s\n", esl_opt_GetString(go, "--infile"));
    esl_sqfile_Close(sqfp);
    N = seqs_to_aln->nseq;
  }
  else if(esl_opt_IsUsed(go, "-L")) {
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
  }
  else if(do_random) {
    double *dnull;
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);
    /* get gamma[0] from the QDB calc alg, which will serve as the length distro for random seqs */
    int safe_windowlen = cm->clen * 2;
    double **gamma = NULL;
    while(!(BandCalculationEngine(cm, safe_windowlen, DEFAULT_BETA, TRUE, NULL, NULL, &(gamma), NULL))) {
      safe_windowlen *= 2;
      /* This is a temporary fix becuase BCE() overwrites gamma, leaking memory
       * Probably better long-term for BCE() to check whether gamma is already allocated
       */
      FreeBandDensities(cm, gamma);
      if(safe_windowlen > (cm->clen * 1000)) cm_Fail("Error trying to get gamma[0], safe_windowlen big: %d\n", safe_windowlen);
    }
    seqs_to_aln = RandomEmitSeqsToAln(r, cm->abc, dnull, 1, N, gamma[0], safe_windowlen, FALSE);
    FreeBandDensities(cm, gamma);
    free(dnull);
  }
  else /* don't randomly generate seqs, emit them from the CM */
    seqs_to_aln = CMEmitSeqsToAln(r, cm, 1, N, FALSE, NULL, FALSE);

  if(esl_opt_GetBoolean(go, "--i27")) { 
    SetMarginalScores_reproduce_bug_i27(cm);
  }
  CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &cons);
  
  emap = CreateEmitMap(cm);

  for (i = 0; i < N; i++)
    {
      L = seqs_to_aln->sq[i]->n;
      dsq = seqs_to_aln->sq[i]->dsq;
      cm->search_opts &= ~CM_SEARCH_INSIDE;
      
#if 0
      esl_stopwatch_Start(w);
      if((status = FastCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, 0., NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      esl_stopwatch_Start(w);
      if((status = RefCYKScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, 0., NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "RefCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      esl_stopwatch_Start(w);
      if((status = RefTrCYKScan(cm, errbuf, trsmx, dsq, 1, L, 0., NULL, FALSE, 0., NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "RefTrCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      if(esl_opt_GetBoolean(go, "--orig")) { 
	esl_stopwatch_Start(w);
	sc = TrCYK_Inside(cm, dsq, L, 0, 1, L, FALSE, &tr);
	printf("%4d %-30s %10.4f bits ", (i+1), "TrCYK_Inside():   ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      if(esl_opt_GetBoolean(go, "--dc")) { 
	esl_stopwatch_Start(w);
	sc = TrCYK_DnC(cm, dsq, L, 0, 1, L, FALSE);
	printf("%4d %-30s %10.4f bits ", (i+1), "TrCYK_DnC():      ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }
#endif
      if(esl_opt_GetBoolean(go, "--hb")) { 
	cm->align_opts  |= CM_ALIGN_HBANDED;

	esl_stopwatch_Start(w);
	if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, 
				   TRUE,  /* doing search? */
				   FALSE,  /* doing trCYK/Inside/Outside? */
				   0)) != eslOK) cm_Fail(errbuf);
	esl_stopwatch_Stop(w);
	printf("%4d %-30s %17s", i+1, "HMM Band calc:", "");
	esl_stopwatch_Display(stdout, w, "CPU time: ");

	PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, L);

	esl_stopwatch_Start(w);
	if((status = FastCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, FALSE, cm->hbmx, 1028., 0., NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScanHB(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	esl_stopwatch_Start(w);
	if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, 
				   TRUE,  /* doing search? */
				   TRUE,  /* doing trCYK/Inside/Outside? */
				   0)) != eslOK) cm_Fail(errbuf);
	esl_stopwatch_Stop(w);
	printf("%4d %-30s %17s", i+1, "HMM Band calc:", "");
	esl_stopwatch_Display(stdout, w, "CPU time: ");

	PrintDPCellsSaved_jd(cm, cm->cp9b->jmin, cm->cp9b->jmax, cm->cp9b->hdmin, cm->cp9b->hdmax, L);

	cm_tr_hb_mx_Destroy(trhbmx);
	trhbmx = cm_tr_hb_mx_Create(cm);
	esl_stopwatch_Start(w);
	if((status = FastTrCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, FALSE, trhbmx, 1028., 0., NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastTrCYKScanHB(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");

	cm_tr_hb_mx_Destroy(trhbmx);
	trhbmx = cm_tr_hb_mx_Create(cm);
	esl_stopwatch_Start(w);
	if((status = FastTrCYKScanHBSkipSlow(cm, errbuf, dsq, 1, L, 0., NULL, FALSE, trhbmx, 1028., 0., emap, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastTrCYKScanHBSkipSlow(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }

      if(esl_opt_GetBoolean(go, "--ins")) { 
	cm->search_opts  |= CM_SEARCH_INSIDE;

	esl_stopwatch_Start(w);
	if((status = RefITrInsideScan(cm, errbuf, trsmx, dsq, 1, L, 0., NULL, FALSE, 0., NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "RefITrInsideScan(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	
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
	
	esl_stopwatch_Start(w);
	if((status = FastFInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "FastFInsideScan(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
	
	esl_stopwatch_Start(w);
	if((status = RefFInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, FALSE, NULL, &sc)) != eslOK) cm_Fail(errbuf);
	printf("%4d %-30s %10.4f bits ", (i+1), "RefFInsideScan(): ", sc);
	esl_stopwatch_Stop(w);
	esl_stopwatch_Display(stdout, w, " CPU time: ");
      }	
      /*
      fali = CreateFancyAli(cm->abc, tr, cm, cons, dsq, FALSE, NULL);
      ParsetreeDump(stdout, tr, cm, dsq, NULL, NULL);
      PrintFancyAli(stdout, fali, 0, FALSE, FALSE, 60);
      FreeParsetree(tr);
      FreeFancyAli(fali);
      */

      printf("\n");
    }
  FreeCM(cm);
  FreeSeqsToAln(seqs_to_aln);
  FreeEmitMap(emap);
  cm_tr_hb_mx_Destroy(trhbmx);
#if 0
  cm_FreeTrScanMatrix(cm, trsmx);
#endif
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  cm_Fail("memory allocation error");
  return 0; /* never reached */
}
#endif /*IMPL_TRUNC_SEARCH_BENCHMARK*/



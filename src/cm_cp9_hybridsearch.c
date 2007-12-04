/* cm_cp9_hybridsearch.c
 * EPN, Mon Oct 29 06:24:24 2007
 * 
 * Implementation of a hybrid CYK/Viterbi scanning algorithm.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

#define CP9TSC(s,k) (cp9tsc[(k) * cp9O_NTRANS + (s)])


/*
 * Function: cm_cp9_HybridScan()
 * 
 * Purpose:  Scan an digitized sequence with a hybrid HMM Viterbi/CM CYK
 *           DP algorithm. 
 * 
 * Args:     cm          - the CM
 *           errbuf    - char buffer for error messages
 *           CP9_MX     - the CP9 matrix
 *           dsq        - the digitized sequence
 *           hsi        - hybrid scan info, describing which parts of model are modelled by HMM/CM
 *           i0         - start of target subsequence (1 for full seq)
 *           j0         - end of target subsequence (L for full seq)
 *           W          - max d: max size of a hit
 *           dmin       - minimum d for each CM state 
 *           dmax       - maximum d for each CM state 
 *           cutoff     - minimum score to report
 *           results   - search_results_t to add to; if NULL, don't keep results
 *           ret_psc    - RETURN: int log odds Forward score for each end point [0..(j0-i0+1)]
 *           ret_maxres - RETURN: start position that gives maximum score max argmax_i sc[i]
 *           ret_sc     - RETURN: score of best hit
 *                        max log P(S|M)/P(S|R), for argmax subseq S of input seq i0..j0,
 *
 * Returns:  eslOK on success;
 *           eslEINCOMPAT on contract violation;
 */
int
cm_cp9_HybridScan(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, HybridScanInfo_t *hsi, int i0, int j0, int W, 
		  float cutoff, search_results_t *results, int **ret_psc, int *ret_maxres, float *ret_sc)
{
  /* HMM related variables */
  int          status;
  GammaHitMx_t *gamma;      /* semi-HMM for hit resoultion                                  */
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /* j-W: position in the subsequence                             */
  int          ip;          /* i': relative position in the subsequence                     */
  int          cm_cur, cm_prv; /* rows in alpha DP matrix 0 or 1                            */
  int          k;           /* CP9 HMM node position                                        */
  int          L;           /* j0-i0+1: subsequence length                                  */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cm->cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int         *scA;         /* prob (seq from j0..jp | HMM) [0..jp..cm->cp9->M]             */
  float        fsc;         /* float log odds score                                         */
  float        best_sc;     /* score of best hit overall                                    */
  float        best_pos;    /* residue (j) giving best_sc, where best hit ends              */
  int          c;           /* counter for EL states                                        */
  int          M;           /* cm->cp9->M, query length, number of consensus nodes of model */

  /* CM related variables */
  int       yoffset;		/* offset to a child state */
  int       d;			/* a subsequence length, 0..W */
  int       cp9_cur, cp9_prv;   /* rows in HMM DP matrix [0..M+1]*/
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  int      *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  int     **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* HMM related ontract checks */
  if(cm->cp9 == NULL)                  ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), cm->cp9 is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), mx is NULL.\n");
  if(mx->M != cm->clen)                ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), mx->M != cm->clen.\n");
  if(cm->clen != cm->cp9->M)           ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), cm->clen != cm->cp9->M.\n");
 
  /* CM related contract checks */
  if(! cm->flags & CMH_BITS)               ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANMATRIX)         ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), cm->smx is not valid.\n");
  if(j0 < i0)                              ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                          ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE)   ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), CM_SEARCH_INSIDE flag raised");
  if(hsi->smx == NULL)                     ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), hsi has no ScanMatrix");
  if(! (hsi->smx->flags & cmSMX_HAS_INT))  ESL_FAIL(eslEINCOMPAT, errbuf, "cm_cp9_HybridScan(), hsi's ScanMatrix's cmSMX_HAS_INT flag is not raised");

  /* make pointers to the ScanMatrix/CM data for convenience */
  ScanMatrix_t *smx   = hsi->smx;         /* we use hsi's ScanMatrix NOT the CM's */
  int   ***alpha      = smx->ialpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  int   ***alpha_begl = smx->ialpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = smx->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = smx->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = smx->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = smx->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  int   **esc_vAA     = cm->ioesc;        /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					   * and all possible emissions a (including ambiguities) */
  /* determine if we're doing banded/non-banded */
  if(smx->dmin != NULL && smx->dmax != NULL) do_banded = TRUE;

  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  L = j0-i0+1;
  M = cm->cp9->M;
  if (W > L) W = L; 
  if (W > smx->W) ESL_FAIL(eslEINCOMPAT, errbuf, "FastCYKScan, W: %d greater than smx->W: %d\n", W, smx->W);
  int const *cp9tsc = cm->cp9->otsc; /* ptr to efficiently ordered CP9 transition scores           */

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. */
  if(results != NULL) gamma = CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff, FALSE);
  else                gamma = NULL;

  /* Grow DP matrix, if nec, to W+1 rows, it stays M+1 columns */
  GrowCP9Matrix(mx, NULL, W, M, &mmx, &imx, &dmx, &elmx, &erow); 
  ESL_DPRINTF1(("cm_cp9_HybridScan(): CP9 matrix size: %.8f Mb rows: %d.\n", mx->size_Mb, mx->rows));

  /* scA will hold P(seq up to j | Model) in int log odds form */
  ESL_ALLOC(scA, sizeof(int) * (j0-i0+2));
			
  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  elmx[0][0]= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY;   

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(int) * (W+1)));
  esl_vec_ISet(sc_v, (W+1), -INFTY);

  /* precalculate the initial scores for all cells of the CM matrix */
  init_scAA = ICalcInitDPScores(cm);

  /* HMM: allow a full HMM parse through all delete states, with 0 emissions,
   * this isn't really a hybrid hit, but it would be nonsense if reported anyway,
   * and I can't imagine it ever being reported.
   */
  int tmp_sc;
  for (k = 1; k <= M; k++) {
    mmx[0][k] = imx[0][k] = elmx[0][k] = -INFTY;      /* need seq to get here */
    tmp_sc = ESL_MAX(mmx[0][k-1] + CP9TSC(cp9O_MD,k-1),
		 imx[0][k-1]     + CP9TSC(cp9O_ID,k-1));
    tmp_sc = ESL_MAX(tmp_sc, dmx[0][k-1] + CP9TSC(cp9O_DD,k-1));
    dmx[0][k] = tmp_sc;
  }
  /* We can do a full parse through all delete states. Hard to believe this will ever be optimal */
  erow[0]  = dmx[0][M] + CP9TSC(cp9O_DM,M);
  scA[0]   = erow[0];
  fsc      = Scorify(scA[0]);
  if(fsc > best_sc) { best_sc = fsc; best_pos= i0-1; }

  printf("dmin[hsi->v_last: %4d]: %d\n", hsi->v_last, dnAA[W][hsi->v_last]);
  printf("dmax[hsi->v_last: %4d]: %d\n", hsi->v_last, dxAA[W][hsi->v_last]);

  /*****************************************************************
   * The main loop: scan the sequence from position i0 to j0.
   *****************************************************************/
  /* Recursion. */
  for (j = i0; j <= j0; j++)
    {
      /*********************************************************************/
      /* CM section */
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cm_cur  = j%2;
      cm_prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else {          dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);
      
      for (v = hsi->v_last; v >= hsi->v_first; v = hsi->v_prv[v]) {
	int sc;
	ESL_DASSERT1((hsi->v_mb[v] == MB_CM));
	if(cm->sttype[v] == E_st) continue;
	int const *esc_v = esc_vAA[v]; 
	int const *tsc_v = cm->itsc[v];
	int emitmode = Emitmode(cm->sttype[v]);
	
	jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cm_cur;
	jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? cm_prv : cm_cur;
	sd   = StateDelta(cm->sttype[v]);
	cnum = cm->cnum[v];
	dn   = dnA[v];
	dx   = dxA[v];
	/* if we emit right, precalc score of emitting res j from state v */
	int esc_j = -INFTY;
	if(cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	  esc_j = esc_v[dsq[j]];

	if(cm->sttype[v] == B_st) {
	  w = cm->cfirst[v]; /* BEGL_S */
	  y = cm->cnum[v];   /* BEGR_S */
	  for (d = dnA[v]; d <= dxA[v]; d++) {
	    /* k is the length of the right fragment */
	    /* Careful, make sure k is consistent with bands in state w and state y. */
	    if(do_banded) {
	      kmin = ESL_MAX(dmin[y], (d-dmax[w]));
	      kmin = ESL_MAX(kmin, 0);
	      kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	    }
	    else { kmin = 0; kmax = d; }

	    sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for B_st */
	    for (k = kmin; k <= kmax; k++) 
	      sc = ESL_MAX(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	    alpha[jp_v][v][d] = sc;
	    /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	  }
	}
	else if (cm->stid[v] == BEGL_S) {
	  y = cm->cfirst[v]; 
	  for (d = dnA[v]; d <= dxA[v]; d++) {
	    sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S */
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	      sc = ESL_MAX (sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);
	    alpha_begl[jp_v][v][d] = sc;
	    /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	  }
	}
	else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
	  y    = cm->cfirst[v];
	  dp_y = dn - sd; /* initial dp_y, we increment it at end of 'for(d = ...' loop */
	  i    = j-dn+1;  /* initial i,    we decrement it when we access it, inside each possible case of the switch (cnum) below */
	  
	  int const *arow0;
	  int const *arow1;
	  int const *arow2;
	  int const *arow3;
	  int const *arow4;
	  int const *arow5;
	
	  /* Note: order of cnum cases in switch and cases in each
	   * nested emitmode switch is based on empirical
	   * frequency in large test set, more frequent guys come
	   * earlier, so average num calcs in each switch is
	   * minimized.
	   */
	
	  switch (cnum) {
	  case 3: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc = ESL_MAX(arow2[dp_y] + tsc_v[2],
			   arow1[dp_y] + tsc_v[1]);		
	      sc = ESL_MAX(sc, init_scAA[v][dp_y]);
	      sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
	    
	      switch (emitmode) {
	      case EMITLEFT:
		sc += esc_v[dsq[i--]];
		break;
	      case EMITRIGHT:
		sc += esc_j;
		break;		
	      } /* end of switch (cm->sttype[v]) */
	      alpha[jp_v][v][d] = sc;
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;
	  
	  case 6: /* necessarily 2 inserts */
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    arow3 = (int * const) alpha[jp_y][y+3];
	    arow4 = (int * const) alpha[jp_y][y+4];
	    arow5 = (int * const) alpha[jp_y][y+5];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc = ESL_MAX(arow5[dp_y] + tsc_v[5],
			   init_scAA[v][dp_y]);
	      sc = ESL_MAX(sc, arow4[dp_y] + tsc_v[4]);		
	      sc = ESL_MAX(sc, arow3[dp_y] + tsc_v[3]);		
	      sc = ESL_MAX(sc, arow2[dp_y] + tsc_v[2]);		
	      sc = ESL_MAX(sc, arow1[dp_y] + tsc_v[1]);		
	      sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
	    
	      switch (emitmode) {
	      case EMITLEFT:
		sc += esc_v[dsq[i--]];
		break;
	      case EMITRIGHT:
		sc += esc_j;
		break;		
	      } /* end of switch (cm->sttype[v]) */
	      alpha[jp_v][v][d] = sc;
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;

	  case 4: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    arow3 = (int * const) alpha[jp_y][y+3];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc = ESL_MAX(arow3[dp_y] + tsc_v[3],
			   arow2[dp_y] + tsc_v[2]);		
	      sc = ESL_MAX(sc, init_scAA[v][dp_y]);
	      sc = ESL_MAX(sc, arow1[dp_y] + tsc_v[1]);		
	      sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
		
	      switch (emitmode) {
	      case EMITLEFT:
		sc += esc_v[dsq[i--]];
		break;
	      case EMITRIGHT:
		sc += esc_j;
		break;		
	      } /* end of switch (cm->sttype[v]) */
	      alpha[jp_v][v][d] = sc;
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;

	  case 5: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    arow3 = (int * const) alpha[jp_y][y+3];
	    arow4 = (int * const) alpha[jp_y][y+4];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc = ESL_MAX(arow4[dp_y] + tsc_v[4],
			   arow3[dp_y] + tsc_v[3]);		
	      sc = ESL_MAX(sc, init_scAA[v][dp_y]);
	      sc = ESL_MAX(sc, arow1[dp_y] + tsc_v[1]);		
	      sc = ESL_MAX(sc, arow2[dp_y] + tsc_v[2]);		
	      sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		

	      switch (emitmode) {
	      case EMITRIGHT:
		sc += esc_j;
		break;		
	      case EMITLEFT:
		sc += esc_v[dsq[i--]];
		break;
		/* MP states can't have 5 children */
	      } /* end of switch (cm->sttype[v]) */
	      alpha[jp_v][v][d] = sc;
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;

	  case 2: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc = ESL_MAX(arow1[dp_y] + tsc_v[1],
			   init_scAA[v][dp_y]);
	      sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
	      switch (emitmode) {
	      case EMITLEFT:
		sc += esc_v[dsq[i--]];
		break;
	      case EMITRIGHT:
		sc += esc_j;
		break;		
	      } /* end of switch (cm->sttype[v]) */
	      alpha[jp_v][v][d] = sc;
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;
	  } /* end of switch(cnum) */
	  /* for (d = dn; d <= dx; d++) 
	     printf("alpha[j:%d][v:%d][d:%d]: %10.4f\n", j, v, d, alpha[jp_v][v][d]); */
	} /* end of else if (v == IL_st || v == IR_st) */
	else { /* this else is entered if cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st */
	  y    = cm->cfirst[v];
	  dp_y = dn - sd; /* initial dp_y, we increment it at end of 'for(d = ...' loop */
	  i    = j-dn+1;  /* initial i,    we decrement it when we access it, inside each possible case of the switch (cnum) below */

	  int const *arow0;
	  int const *arow1;
	  int const *arow2;
	  int const *arow3;
	  int const *arow4;
	  int const *arow5;

	  /* Note: order of cnum cases in switch and cases in each
	   * nested emitmode switch is based on empirical
	   * frequency in large test set, more frequent guys come
	   * earlier, so average num calcs in each switch is
	   * minimized.
	   */

	  switch (cnum) {
	  case 3: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc_v[d] = ESL_MAX(arow2[dp_y] + tsc_v[2],
				arow1[dp_y] + tsc_v[1]);		
	      sc_v[d] = ESL_MAX(sc_v[d], init_scAA[v][dp_y]);
	      sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;

	  case 6: /* necessarily 2 inserts */
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    arow3 = (int * const) alpha[jp_y][y+3];
	    arow4 = (int * const) alpha[jp_y][y+4];
	    arow5 = (int * const) alpha[jp_y][y+5];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc_v[d] = ESL_MAX(arow5[dp_y] + tsc_v[5],
				init_scAA[v][dp_y]);
	      sc_v[d] = ESL_MAX(sc_v[d], arow4[dp_y] + tsc_v[4]);		
	      sc_v[d] = ESL_MAX(sc_v[d], arow3[dp_y] + tsc_v[3]);		
	      sc_v[d] = ESL_MAX(sc_v[d], arow2[dp_y] + tsc_v[2]);		
	      sc_v[d] = ESL_MAX(sc_v[d], arow1[dp_y] + tsc_v[1]);		
	      sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;

	  case 4: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    arow3 = (int * const) alpha[jp_y][y+3];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc_v[d] = ESL_MAX(arow3[dp_y] + tsc_v[3],
				arow2[dp_y] + tsc_v[2]);		
	      sc_v[d] = ESL_MAX(sc_v[d], init_scAA[v][dp_y]);
	      sc_v[d] = ESL_MAX(sc_v[d], arow1[dp_y] + tsc_v[1]);		
	      sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	    } /* end of for(d = dn; d <= dx; d++) */
	    break;

	  case 5: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    arow2 = (int * const) alpha[jp_y][y+2];
	    arow3 = (int * const) alpha[jp_y][y+3];
	    arow4 = (int * const) alpha[jp_y][y+4];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc_v[d] = ESL_MAX(arow4[dp_y] + tsc_v[4],
				arow3[dp_y] + tsc_v[3]);		
	      sc_v[d] = ESL_MAX(sc_v[d], init_scAA[v][dp_y]);
	      sc_v[d] = ESL_MAX(sc_v[d], arow1[dp_y] + tsc_v[1]);		
	      sc_v[d] = ESL_MAX(sc_v[d], arow2[dp_y] + tsc_v[2]);		
	      sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	    } /* end of for (d = dn; d <= dx; d++, dp_y++) */
	    break;

	  case 2: 
	    arow0 = (int * const) alpha[jp_y][y];
	    arow1 = (int * const) alpha[jp_y][y+1];
	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc_v[d] = ESL_MAX(arow1[dp_y] + tsc_v[1],
				init_scAA[v][dp_y]);
	      sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	    }
	    break; 
	  } /* end of switch(cnum) */
	  /* add in emission score (if any), and set alpha[jp_v][v][d] cell */
	  switch (emitmode) {
	  case EMITLEFT:
	    for (d = dn; d <= dx; d++) {
	      alpha[jp_v][v][d] = sc_v[d] + esc_v[dsq[i--]];
	    }
	    break;
	  case EMITNONE:
	    for (d = dn; d <= dx; d++)
	      alpha[jp_v][v][d] = sc_v[d];
	    break;
	  case EMITRIGHT:
	    for (d = dn; d <= dx; d++) {
	      alpha[jp_v][v][d] = sc_v[d] + esc_j;
	    }
	    break;		
	  case EMITPAIR:
	    for (d = dn; d <= dx; d++) {
	      alpha[jp_v][v][d] = sc_v[d] + esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	    }
	    break;
	  } /* end of switch (emitmode) */
	} /* end of else (cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st) */
	/* if(cm->stid[v] != BEGL_S)
	   for (d = dn; d <= dx; d++) { printf("alpha[j:%4d][v:%4d][d:%4d]: %.5f\n", j, v, d, alpha[jp_v][v][d]); }*/
      } /* end of for (v = hsi->v_last; v >= hsi->v_first; v = hsi->v_prv[v]) */
      /*********************************************************************/
      /* HMM section */
      int const *isc = cm->cp9->isc[dsq[j]];
      int const *msc = cm->cp9->msc[dsq[j]];

      int endsc     = -INFTY;
      int el_selfsc = cm->cp9->el_selfsc;
      //int sc;

      jp = j-i0+1;     /* jp is relative position in the sequence 1..L */
      cp9_cur = (j-i0+1) % (W+1);
      cp9_prv = (j-i0)   % (W+1);

      /* The 1 difference between a Viterbi scanner and the 
       * regular Viterbi. In non-scanner parse must begin in B at
       * position 0 (i0-1), in scanner we can start at any position 
       * in the seq. */
      mmx[cp9_cur][0]  = 0;
      dmx[cp9_cur][0]  = -INFTY;  /*D_0 is non-existent*/
      elmx[cp9_cur][0] = -INFTY;  /*no EL state for node 0 */

      tmp_sc = ESL_MAX(mmx[cp9_prv][0] + CP9TSC(cp9O_MI,0),
		   imx[cp9_prv][0] + CP9TSC(cp9O_II,0));
      tmp_sc = ESL_MAX(tmp_sc, dmx[cp9_prv][0] + CP9TSC(cp9O_DI,0));
      imx[cp9_cur][0] = ESL_MAX(tmp_sc + isc[0], -INFTY);
      
      int kp;
      int dx;
      for (k = ESL_MAX(1, hsi->k_firstr); k >= 0; k = hsi->k_nxtr[k]) {
	int sc;
	/* printf("\tk: %d\n", k); */
	/*********************************************************************/
	if(hsi->k_mb[(k-1)] == MB_CP9 && hsi->k_mb[k] == MB_CP9) { 
	  /* normal Viterbi recursion, this node (k) and previous node (k-1) both modelled by the CP9 HMM */
	  /*match state*/
	  sc = ESL_MAX(    mmx[cp9_prv][k-1] + CP9TSC(cp9O_MM,k-1),
			   imx[cp9_prv][k-1] + CP9TSC(cp9O_IM,k-1));
	  sc = ESL_MAX(sc, dmx[cp9_prv][k-1] + CP9TSC(cp9O_DM,k-1));
	  sc = ESL_MAX(sc, mmx[cp9_prv][0]   + CP9TSC(cp9O_BM,k));
	  /* check possibility we came from an EL, if they're valid */
	  for(c = 0; c < cm->cp9->el_from_ct[k]; c++) { /* el_from_ct[k] is >= 0 */
	    kp = cm->cp9->el_from_idx[k][c];
	    if(hsi->k_mb[kp] == MB_CP9) /* only can come from kp's EL if kp is modelled by the CP9 */
	      sc = ESL_MAX(sc, elmx[cp9_prv][cm->cp9->el_from_idx[k][c]]);
	  }
	  /* transition penalty to EL incurred when EL was entered */
	  mmx[cp9_cur][k] = ESL_MAX(sc + msc[k], -INFTY);
	  
	  /* E state update */
	  endsc = ESL_MAX(endsc, mmx[cp9_cur][k] + CP9TSC(cp9O_ME,k));
	  
	  /*insert state*/
	  sc = ESL_MAX(    mmx[cp9_prv][k] + CP9TSC(cp9O_MI,k),
			   imx[cp9_prv][k] + CP9TSC(cp9O_II,k));
	  sc = ESL_MAX(sc, dmx[cp9_prv][k] + CP9TSC(cp9O_DI,k));
	  imx[cp9_cur][k] = ESL_MAX(sc + isc[k], -INFTY);
	  
	  /*delete state*/
	  sc = ESL_MAX(    mmx[cp9_cur][k-1] + CP9TSC(cp9O_MD,k-1),
			   imx[cp9_cur][k-1] + CP9TSC(cp9O_ID,k-1));
	  sc = ESL_MAX(sc, dmx[cp9_cur][k-1] + CP9TSC(cp9O_DD,k-1));
	  dmx[cp9_cur][k] = sc;
	  
	  /*el state*/
	  sc = -INFTY;
	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								    HMM nodes that map to right half of a MATP_MP) */
	    {
	      sc = ESL_MAX(sc, mmx[cp9_cur][k]  + CP9TSC(cp9O_MEL,k));
	      sc = ESL_MAX(sc, elmx[cp9_prv][k] + el_selfsc);
	    }
	  elmx[cp9_cur][k] = sc;
	}
	/*********************************************************************/
	else if(hsi->k_mb[(k-1)] == MB_CM && hsi->k_mb[k] == MB_CP9) { 
	  /* previous node modelled by CM, this node modelled by CP9 HMM */

	  /*match state*/
	  kp = hsi->k_prv[k];
	  sc = ESL_MAX(    mmx[cp9_prv][kp],// + CP9TSC(cp9O_MM,k-1), /* not sure about the transition from node k-1 ? (should it be from kp?), or NO transition score? */
			   imx[cp9_prv][kp]);// + CP9TSC(cp9O_IM,k-1));
	  sc = ESL_MAX(sc, dmx[cp9_prv][kp]);// + CP9TSC(cp9O_DM,k-1));
	  sc = ESL_MAX(sc, mmx[cp9_prv][0]   + CP9TSC(cp9O_BM,k));
	  /* check possibility we came from an EL, if they're valid */
	  for(c = 0; c < cm->cp9->el_from_ct[k]; c++) { /* el_from_ct[k] is >= 0 */
	    kp = cm->cp9->el_from_idx[k][c];
	    if(hsi->k_mb[kp] == MB_CP9) /* only can come from kp's EL if kp is modelled by the CP9 */
	      sc = ESL_MAX(sc, elmx[cp9_prv][kp]);
	  }
	  /* transition penalty to EL incurred when EL was entered */
	  mmx[cp9_cur][k] = ESL_MAX(sc + msc[k], -INFTY);
	  
	  /* E state update */
	  endsc = ESL_MAX(endsc, mmx[cp9_cur][k] + CP9TSC(cp9O_ME,k));
	  
	  /*insert state*/
	  sc = ESL_MAX(    mmx[cp9_prv][k] + CP9TSC(cp9O_MI,k),
			   imx[cp9_prv][k] + CP9TSC(cp9O_II,k));
	  sc = ESL_MAX(sc, dmx[cp9_prv][k] + CP9TSC(cp9O_DI,k));
	  imx[cp9_cur][k] = ESL_MAX(sc + isc[k], -INFTY);
	  
	  /*delete state*/
	  kp = hsi->k_prv[k];
	  sc = ESL_MAX(    mmx[cp9_cur][kp],// + CP9TSC(cp9O_MD,k-1),
			   imx[cp9_cur][kp]);// + CP9TSC(cp9O_ID,k-1));
	  sc = ESL_MAX(sc, dmx[cp9_cur][kp]);
	  dmx[cp9_cur][k] = sc;
	  
	  /*el state*/
	  sc = -INFTY;
	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								    HMM nodes that map to right half of a MATP_MP) */
	    {
	      sc = ESL_MAX(sc, mmx[cp9_cur][k]  + CP9TSC(cp9O_MEL,k));
	      sc = ESL_MAX(sc, elmx[cp9_prv][k] + el_selfsc);
	    }
	  elmx[cp9_cur][k] = sc;
	}
	/***********************************************************************/
	else if(hsi->k_mb[k-1] == MB_CP9 && hsi->k_mb[k] == MB_CM) { 
	  /* previous node modelled by CP9 HMM, this node modelled by CM */

	  int m_v1 = cm->cp9map->hns2cs[k][HMMMATCH][0];
	  int m_v2 = cm->cp9map->hns2cs[k][HMMMATCH][1];
	  int i_v1 = cm->cp9map->hns2cs[k][HMMINSERT][0];
	  ESL_DASSERT1((cm->cp9map->hns2cs[k][HMMINSERT][1] == -1)); /* right? */
	  int d_v1 = cm->cp9map->hns2cs[k][HMMDELETE][0];
	  int d_v2 = cm->cp9map->hns2cs[k][HMMDELETE][1];
	  
	  //printf("j: %d k: %d m_v1: %d m_v2: %d i_v1: %d i_v2: %d d_v1: %d d_v2: %d\n", j, k, m_v1, m_v2, i_v1, i_v2, d_v1, d_v2);
	  
	  /*match state*/
	  int m_max = -INFTY;
	  int i_max = -INFTY;
	  int d_max = -INFTY;
	  
	  /* determine m_max */
	  i = j - dnA[m_v1] - 1;
	  dx = ESL_MIN(dxA[m_v1], (j-1));
	  for(d = dnA[m_v1]; d <= dx; d++, i--) {
	    ip = i % (W+1);
	    m_max = ESL_MAX(m_max, mmx[ip][k-1] + alpha[cm_cur][m_v1][d]);
	  }
	  if(m_v2 != -1) { 
	    i = j - dnA[m_v2] - 1;
	    dx = ESL_MIN(dxA[m_v2], (j-1));
	    for(d = dnA[m_v2]; d <= dx; d++, i--) {
	      ip = i % (W+1);
	      m_max = ESL_MAX(m_max, mmx[ip][k-1] + alpha[cm_cur][m_v2][d]);
	    }
	  }
	  
	  /* determine i_max (assumes i_v2 == -1, as asserted above) */
	  i = j - dnA[i_v1] - 1;
	  dx = ESL_MIN(dxA[i_v1], (j-1));
	  for(d = dnA[i_v1]; d <= dx; d++, i--) {
	    ip = i % (W+1);
	    i_max = ESL_MAX(i_max, imx[ip][k-1] + alpha[cm_cur][i_v1][d]);
	  }
	  
	  /* determine d_max */
	  i = j - dnA[d_v1] - 1;
	  dx = ESL_MIN(dxA[d_v1], (j-1));
	  for(d = dnA[d_v1]; d <= dx; d++, i--) {
	    ip = i % (W+1);
	    d_max = ESL_MAX(d_max, dmx[ip][k-1] + alpha[cm_cur][d_v1][d]);
	  }
	  if(d_v2 != -1) { 
	    i = j - dnA[d_v2] - 1;
	    dx = ESL_MIN(dxA[d_v2], (j-1));
	    for(d = dnA[d_v2]; d <= dx; d++, i--) {
	      ip = i % (W+1);
	      d_max = ESL_MAX(d_max, dmx[ip][k-1] + alpha[cm_cur][d_v2][d]);
	    }
	  }
	  sc = ESL_MAX(    m_max + CP9TSC(cp9O_MM,k-1), 
			   i_max + CP9TSC(cp9O_IM,k-1));
	  sc = ESL_MAX(sc, d_max + CP9TSC(cp9O_DM,k-1));
	  //sc = ESL_MAX(sc, mmx[cp9_prv][0]   + CP9TSC(cp9O_BM,k)); /* local begin into M_k */
	  /* ? IS THIS RIGHT ? */ sc = ESL_MAX(sc, m_max + CP9TSC(cp9O_BM,k)); /* local begin into M_k */
	  mmx[cp9_cur][k] = ESL_MAX(sc, -INFTY);
	  
	  /* check possibility we came from an EL, if they're valid */
	  for(c = 0; c < cm->cp9->el_from_ct[k]; c++) { /* el_from_ct[k] is >= 0 */
	    kp = cm->cp9->el_from_idx[k][c];
	    if(hsi->k_mb[kp] == MB_CP9) /* only can come from kp's EL if kp is modelled by the CP9 */
	      sc = ESL_MAX(sc, elmx[cp9_prv][kp]);
	  }
	  /* transition penalty to EL incurred when EL was entered */
	  /* DO NOT ADD CONTRIBUTION OF EMITTING POSITION j FROM HMM NODE k, 
	   * it was emitted by the CM */
	  
	  /* insert state, comes from this node, not from k-1, so it's independent of alpha, and the recursion is normal Viterbi */
	  sc = ESL_MAX(    mmx[cp9_prv][k] + CP9TSC(cp9O_MI,k),
			   imx[cp9_prv][k] + CP9TSC(cp9O_II,k));
	  sc = ESL_MAX(sc, dmx[cp9_prv][k] + CP9TSC(cp9O_DI,k));
	  imx[cp9_cur][k] = ESL_MAX(sc, -INFTY);
	  /* DO NOT ADD CONTRIBUTION OF EMITTING POSITION j FROM HMM NODE k, 
	   * it was emitted by the CM */
	  
	  /* delete state */
	  /* I *think* we can use m_max, i_max, d_max calc'ed from m_v1, m_v2 */
	  /*match state*/
	  sc = ESL_MAX(    m_max + CP9TSC(cp9O_MD,k-1), 
			   i_max + CP9TSC(cp9O_ID,k-1));
	  sc = ESL_MAX(sc, d_max + CP9TSC(cp9O_DD,k-1));
	  dmx[cp9_cur][k] = ESL_MAX(sc, -INFTY);
	  
	  /* NO EL state update, node k is modelled by the CM, the CM could've modelled this by EL,
	   * and if that was the highest scoring subparse, it would've already been handled, rooted at 
	   * a CM match state (or some CM state that has legal local ends)
	   */
	  elmx[cp9_cur][k] = -INFTY;
	}
	/*printf("mmx[cp9_cur: %d][k: %d] %d\n", cp9_cur, k, mmx[cp9_cur][k]);
	  printf("imx[cp9_cur: %d][k: %d] %d\n", cp9_cur, k, imx[cp9_cur][k]);
	  printf("dmx[cp9_cur: %d][k: %d] %d\n", cp9_cur, k, dmx[cp9_cur][k]);*/
	/***********************************************************************/
      }
      endsc = ESL_MAX(endsc, dmx[cp9_cur][M] + CP9TSC(cp9O_DM,M)); /* transition from D_M -> end */
      endsc = ESL_MAX(endsc, imx[cp9_cur][M] + CP9TSC(cp9O_IM,M)); /* transition from I_M -> end */
      for(c = 0; c < cm->cp9->el_from_ct[M+1]; c++) { /* el_from_ct[k] is >= 0 */
	/* transition penalty to EL incurred when EL was entered */
	kp = cm->cp9->el_from_idx[M+1][c];
	if(hsi->k_mb[kp] == MB_CP9) /* only can come from kp's EL if kp is modelled by the CP9 */
	  endsc = ESL_MAX(endsc, elmx[cp9_cur][cm->cp9->el_from_idx[M+1][c]]);
      }
      erow[cp9_cur] = endsc;
      scA[jp]   = endsc;
      fsc = Scorify(endsc);
     
      /* printf("j: %d fsc: %f\n", j, fsc); */


      if(fsc > best_sc) { best_sc = fsc; best_pos= j; }

      /* determine safe start point, max of j-W+1 and i0 */
      i = ((j-W+1)> i0) ? (j-W+1) : i0;
      ip = i-i0+1;
      if(results != NULL) UpdateGammaHitMxCP9Forward(gamma, ip, jp, fsc, results);
    } /* end loop over end positions j */
  
  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, then we've reported hits in UpdateGammaHitMxCP9Forward() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) TBackGammaHitMxForward(gamma, results, i0, j0);

  /* clean up and exit */
  if(gamma != NULL) FreeGammaHitMx(gamma);

  free(jp_wA);
  free(sc_v);
  free(init_scAA[0]);
  free(init_scAA);
  if(ret_sc != NULL)     *ret_sc     = best_sc;
  if(ret_maxres != NULL) *ret_maxres = best_pos;
  if(ret_psc != NULL)    *ret_psc    = scA;
  else                    free(scA);
  ESL_DPRINTF1(("cm_cp9_FastViterbi() return score: %10.4f\n", best_sc));

  return eslOK;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* NEVERREACHED */
}

/* Function: cm_CalcAvgHitLength()
 * Date:     EPN, Mon Sep 10 10:04:52 2007
 *
 * Purpose:  Calculate the average hit length for each state of a CM using the 
 *           QDB calculation engine, and provided beta.
 *
 * Returns:  void
 */
void
cm_CalcAvgHitLength(CM_t *cm, double beta, float **ret_avglen)
{
  float *avglen = NULL;
  int safe_windowlen;

  safe_windowlen = cm->W * 2;
  while(!(BandCalculationEngine(cm, safe_windowlen, beta, TRUE, NULL, NULL, NULL, &avglen))) {
    safe_windowlen *= 2;
    if(safe_windowlen > (cm->clen * 1000)) cm_Fail("safe_windowlen big: %d\n", safe_windowlen);
  }
  /* int v; for(v = 0; v < cm->M; v++) printf("AVG LEN v: %4d d: %10.4f\n", v, avglen[v]); */

  *ret_avglen = avglen;
}


/* Function: cm_FreeHybridScanInfo()
 * Date:     EPN, Thu Nov  1 14:16:18 2007
 *
 * Purpose:  Free a HybridScanInfo_t object corresponding
 *           to CM <cm>.            
 */
void
cm_FreeHybridScanInfo(HybridScanInfo_t *hsi, CM_t *cm)
{
  int i;
  if(hsi->cm_vcalcs != NULL) free(hsi->cm_vcalcs);
  if(hsi->cp9_vcalcs != NULL)free(hsi->cp9_vcalcs);
  if(hsi->k_mb != NULL)      free(hsi->k_mb);
  if(hsi->k_nxt != NULL)     free(hsi->k_nxt);
  if(hsi->k_prv != NULL)     free(hsi->k_prv);
  if(hsi->k_nxtr != NULL)    free(hsi->k_nxtr);
  if(hsi->k_prvr != NULL)    free(hsi->k_prvr);
  if(hsi->v_mb != NULL)      free(hsi->v_mb);
  if(hsi->v_nxt != NULL)     free(hsi->v_nxt);
  if(hsi->v_prv != NULL)     free(hsi->v_prv);
  if(hsi->v_isroot != NULL)  free(hsi->v_isroot);
  if(hsi->iscandA != NULL)   free(hsi->iscandA);
  if(hsi->avglenA != NULL)   free(hsi->avglenA);
  if(hsi->startA != NULL)    free(hsi->startA);
  if(hsi->firstA != NULL)    free(hsi->firstA);
  if(hsi->lastA != NULL)     free(hsi->lastA);
  for(i = 0; i < hsi->nstarts; i++) 
    if(hsi->withinAA[i] != NULL) free(hsi->withinAA[i]);
  if(hsi->smx != NULL)       cm_FreeScanMatrix(cm, hsi->smx);
  free(hsi->withinAA);
  free(hsi);
}

/* Function: cm_CreateHybridScanInfo()
 * Date:     EPN, Wed Oct 31 05:29:22 2007
 *
 * Purpose:  Given a CM, allocate and fill a HybridScanInfo_t object
 *           for that CM. 
 *           
 *           Step 1: Determine which states are possible sub-CM filter roots:
 *           Fill an array of length cm->M with TRUE, FALSE for each
 *           state. TRUE if we should fit a Gumbel to the state b/c it may be
 *           a good sub-CM filter root state, or FALSE if not. Criteria is that
 *           a state must be a possible local entry state AND must have an average
 *           subseq length of cfg->minlen.
 * 
 *           Step 2: Find which 'start group' each state v belongs to:
 *           Each start group is defined by a start/end state pair. 
 *           Each state v belongs to the group defined by start'/end'
 *           where start' is the maximum start state index that is 
 *           less than v. The assignments of states to groups is most
 *           efficiently done using a push down stack. 
 *            
 * Returns:  Newly allocated HybridScanInfo_t object:
 */
HybridScanInfo_t *
cm_CreateHybridScanInfo(CM_t *cm, double hsi_beta, float full_cm_ncalcs)
{
  int               status;
  HybridScanInfo_t *hsi;             /* the HybridScanInfo_t object we're creating */
  int               nd;              /* counter over nodes of CM */
  int               v;               /* counter over states of CM */
  int               k;               /* counter over node of HMM */
  ESL_STACK        *pda;             /* push down stack for traversing CM */
  int               j_popped;        /* used when determining start groups */
  int               i,j;             /* used when determining start groups */
  int               on_right;        /* used when determining start groups */
  int               safe_windowlen;  /* for calc'ing qdbs */
  int               cp9_ntrans;      /* number of transitions for currently config'ed HMM */
  int               nd_clen;         /* consensus length for a subtree */
  CMEmitMap_t      *emap;            /* emit map for the cm */
  int               vx;              /* a temporary, max v */
  int              *dmin;            /* will become hsi->smx->dmin */
  int              *dmax;            /* will become hsi->smx->dmax */

  /* contract check */
  if(cm->cp9 == NULL) cm_Fail("cm_CreateHybridScanInfo(), cm->cp9 is NULL.\n");

  ESL_ALLOC(hsi, sizeof(HybridScanInfo_t));

  hsi->cm_M   = cm->M;
  hsi->cp9_M  = cm->cp9->M;
  hsi->minlen = DEFAULT_HS_MINLEN;   /* currently 7., not user changeable, should it be? is 7 not good? */
  
  /* determine average length for each subtree (state) */
  hsi->avglen_beta = DEFAULT_HS_BETA; /* currently 1E-15 */
  cm_CalcAvgHitLength(cm, hsi->avglen_beta, &(hsi->avglenA));

  /* get dmin, dmax for the hybrid scanner */
  hsi->beta = hsi_beta;
  safe_windowlen = cm->W * 2;
  while(!(BandCalculationEngine(cm, safe_windowlen, hsi_beta, FALSE, &(dmin), &(dmax), NULL, NULL))) {
    free(dmin);
    free(dmax);
    dmin = NULL;
    dmax = NULL;
    safe_windowlen *= 2;
    if(safe_windowlen > (cm->clen * 1000)) cm_Fail("safe_windowlen big: %d\n", safe_windowlen);
  }
  hsi->W = dmax[0];

  /* create the scan matrix, this stores the matrix, dmin, dmax, etc. */
  hsi->smx = cm_CreateScanMatrix(cm, hsi->W, dmin, dmax, TRUE, TRUE, FALSE); 

  /* determine number of millions of DP calculations per residue for CM and CP9 */
  /* first the full CM, using cm->dmin and cm->dmax, this value is passed in 
   * (different from # calcs using hsi->smx->dmin and hsi->smx->dmax b/c beta used to get 
   * hsi->smx->dmin hsi->smx->dmax may be different than that used for cm->dmin, cm->dmax)
   */
  hsi->full_cm_ncalcs = full_cm_ncalcs; /* this is passed in */
  /* get counts of dp calcs for each subtree in the cm, using hsi->smx->dmin, hsi->smx->dmax */
  if((status = cm_CountSearchDPCalcs(cm, NULL, 1000, hsi->smx->dmin, hsi->smx->dmax, hsi->W, &(hsi->cm_vcalcs), NULL)) != eslOK) cm_Fail("cm_CreateHybridScanInfo(), error counting DP cells.");

  /* we can calc the number of CP9 DP calcs */
  cp9_ntrans = NHMMSTATETYPES * NHMMSTATETYPES; /* 3*3 = 9 transitions in global mode */
  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) cp9_ntrans++;
  if(cm->cp9->flags & CPLAN9_LOCAL_END)   cp9_ntrans++;
  if(cm->cp9->flags & CPLAN9_EL)          cp9_ntrans++;
  hsi->full_cp9_ncalcs = (cp9_ntrans * cm->cp9->M) / 1000000.; /* convert to millions of calcs per residue */
  /* now for each possible 'subtree' of the CP9 */
  emap = CreateEmitMap(cm);
  ESL_ALLOC(hsi->cp9_vcalcs, sizeof(float) * cm->M);
  for(nd = 0; nd < cm->nodes; nd++) { 
    nd_clen   = emap->rpos[nd] - emap->lpos[nd] + 1;
    vx = (nd < (cm->nodes-1)) ? cm->nodemap[nd+1]-1 : cm->M-1;
    for(v = cm->nodemap[nd]; v <= vx; v++) { 
      hsi->cp9_vcalcs[v] = (cp9_ntrans * nd_clen) / 1000000.; /* convert to millions of calcs per residue */
    }
  }
  FreeEmitMap(emap);
  hsi->hybrid_ncalcs = hsi->full_cp9_ncalcs;

  /* determine which states are candidate sub-CM root states 
   * state v is a candidate if a local begin into v is legal and it's avg hit len exceeds our minimum */
  ESL_ALLOC(hsi->iscandA, sizeof(int) * cm->M);
  esl_vec_ISet(hsi->iscandA, cm->M, FALSE);
  hsi->ncands = 1;
  hsi->iscandA[0] = TRUE; /* ROOT_S: has to be TRUE */
  for (nd = 1; nd < cm->nodes; nd++) { /* note: nd 1 must be MATP, MATL, MATR, or BIF */
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) {
      v = cm->nodemap[nd];
      if(hsi->avglenA[v] > hsi->minlen) {
	hsi->iscandA[v] = TRUE;
	hsi->ncands++;
      }
    }
  }

  /* allocate and initialize info on which parts of model will be modelled by cm/cp9,
   * no sub CM roots have been assigned yet, so full model is modelled by cp9 hmm.
   */
  assert(cm->clen == cm->cp9->M);
  ESL_ALLOC(hsi->k_mb,  sizeof(int) * (cm->cp9->M+1));
  ESL_ALLOC(hsi->k_nxt, sizeof(int) * (cm->cp9->M+1));
  ESL_ALLOC(hsi->k_prv, sizeof(int) * (cm->cp9->M+1));
  ESL_ALLOC(hsi->k_nxtr,sizeof(int) * (cm->cp9->M+1));
  ESL_ALLOC(hsi->k_prvr,sizeof(int) * (cm->cp9->M+1));
  for(k = 0; k <= cm->cp9->M; k++) hsi->k_mb[k]  = MB_CP9;

  for(k = 0; k <= cm->cp9->M; k++) hsi->k_prv[k] = k-1; /* k_prv[0] will be -1 */
  for(k = 0; k <  cm->cp9->M; k++) hsi->k_nxt[k] = k+1; 
  hsi->k_first = 0;
  hsi->k_last  = cm->cp9->M;

  for(k = 0; k <= cm->cp9->M; k++) hsi->k_prvr[k] = k-1; /* k_prvr[0] will be -1 */
  for(k = 0; k <  cm->cp9->M; k++) hsi->k_nxtr[k] = k+1; 
  hsi->k_nxt[cm->cp9->M] = -1;
  hsi->k_nxtr[cm->cp9->M] = -1;
  hsi->k_firstr= 0;
  hsi->k_lastr = cm->cp9->M;
  

  ESL_ALLOC(hsi->v_mb,  sizeof(int) * hsi->cm_M);
  ESL_ALLOC(hsi->v_nxt, sizeof(int) * hsi->cm_M);
  ESL_ALLOC(hsi->v_prv, sizeof(int) * hsi->cm_M);
  for(v = 0; v < cm->M; v++) hsi->v_mb[v]  = MB_CP9;
  for(v = 0; v < cm->M; v++) hsi->v_prv[v] = -1; /* will become valid if/when v_mb[v] is changed to MB_CM */
  for(v = 0; v < cm->M; v++) hsi->v_nxt[v] = -1; /* will become valid if/when v_mb[v] is changed to MB_CM */
  hsi->v_first = -1;
  hsi->v_last  = -1;

  hsi->n_v_roots = 0;
  ESL_ALLOC(hsi->v_isroot, sizeof(int) * hsi->cm_M);
  for(v = 0; v < cm->M; v++) hsi->v_isroot[v] = FALSE;

  /* determine each end state's 'start group' 
   * first, allocate data structures for this info 
   */
  hsi->nstarts = CMCountStatetype(cm, S_st);
  ESL_ALLOC(hsi->startA,   sizeof(int) *   hsi->cm_M);
  esl_vec_ISet(hsi->startA, hsi->cm_M, -1);

  ESL_ALLOC(hsi->firstA,    sizeof(int) *  hsi->nstarts);
  ESL_ALLOC(hsi->lastA,     sizeof(int) *  hsi->nstarts);

  ESL_ALLOC(hsi->withinAA, sizeof(int *) * hsi->nstarts);
  for(i = 0; i < hsi->nstarts; i++) 
    {
      ESL_ALLOC(hsi->withinAA[i], sizeof(int) * hsi->nstarts);
      esl_vec_ISet(hsi->withinAA[i], hsi->nstarts, FALSE);
    }
  
  /* traverse the CM using a pda, code stolen and modified from Sean's
   *  cmemit.c:CreateEmitMap() 
   */
  nd   = 0;
  pda  = esl_stack_ICreate();
  j    = 0;
  hsi->firstA[0] = 0;

  esl_stack_IPush(pda, 0);		/* 0 = left side. 1 would = right side. */
  esl_stack_IPush(pda, j);		/* 0 = left side. 1 would = right side. */
  esl_stack_IPush(pda, nd);
  while (esl_stack_IPop(pda, &nd) != eslEOD)
    {
      esl_stack_IPop(pda, &j_popped);
      esl_stack_IPop(pda, &on_right);
      /* printf("nd: %3d j_popped: %3d on_right: %d\n", nd, j_popped, on_right); */
      
      if (on_right) 
	{
	  if(cm->ndtype[nd] == BEGL_nd ||  cm->ndtype[nd] == BEGR_nd) 
	    { 
	      /* we're done with start group j_popped, it is within
	       * all start groups currently in the stack.
	       */
	      for(i = 1; i < pda->n; i += 3) /* += 3 b/c only every third element is a state group index */
		hsi->withinAA[pda->idata[i]][j_popped] = TRUE;
	    }
	}
      else
	{
	  if (cm->ndtype[nd] == BEGL_nd || cm->ndtype[nd] == BEGR_nd) 
	    {
	      j++;
	      hsi->firstA[j] = cm->nodemap[nd];
	    }
	  if (cm->ndtype[nd] == END_nd) 
	    {
	      hsi->lastA[j] = cm->nodemap[nd];
	      for(v = hsi->firstA[j]; v <= hsi->lastA[j]; v++) hsi->startA[v] = j;
	    }
	  if (cm->ndtype[nd] == BIF_nd) 
	    {
	      hsi->lastA[j] = cm->nodemap[nd];
	      for(v = hsi->firstA[j]; v <= hsi->lastA[j]; v++) hsi->startA[v] = j;

				/* push the BIF back on for its right side  */
	      esl_stack_IPush(pda, 1);
	      esl_stack_IPush(pda, j);
	      esl_stack_IPush(pda, nd);
                            /* push node index for right child */
	      esl_stack_IPush(pda, 0);
	      esl_stack_IPush(pda, j);
	      esl_stack_IPush(pda, cm->ndidx[cm->cnum[cm->nodemap[nd]]]);   
                            /* push node index for left child */
	      esl_stack_IPush(pda, 0);
	      esl_stack_IPush(pda, j);
	      esl_stack_IPush(pda, cm->ndidx[cm->cfirst[cm->nodemap[nd]]]); 
	    }
	  else
	    {
	      /* push the node back on for right side */
	      esl_stack_IPush(pda, 1);
	      esl_stack_IPush(pda, j);
	      esl_stack_IPush(pda, nd);
	      /* push next BIF, END node on */
	      if (cm->ndtype[nd] != END_nd) {
		while(cm->ndtype[nd] != BIF_nd && cm->ndtype[nd] != END_nd) nd++;
		esl_stack_IPush(pda, 0);
		esl_stack_IPush(pda, j);
		esl_stack_IPush(pda, nd);
	      }
	    }
	}
    }
  /* printf("\n"); */
  for(v = 0; v < cm->M; v++) { /*printf("startA[%4d]: %d\n", v, hsi->startA[v]);*/ assert(hsi->startA[v] >= 0); }
  /*  for(i = 0; i < hsi->nstarts; i++) {
    printf("firstA[%2d]: %4d\nlastA [%2d]: %4d\n", i, hsi->firstA[i], i, hsi->lastA[i]);
    for(j = 0; j < hsi->nstarts; j++) 
      printf("\twithinAA[%2d][%2d] %d\n", i, j, hsi->withinAA[i][j]);
      printf("\n");
      } */
  
  /* temporary check of withinAA */
  int ileft, iright, jleft, jright;
  emap = CreateEmitMap(cm);
  for(i = 0; i < hsi->nstarts; i++)
    {
      ileft  = emap->lpos[cm->ndidx[hsi->firstA[i]]];
      iright = emap->rpos[cm->ndidx[hsi->firstA[i]]];
      for(j = 0; j < hsi->nstarts; j++)
	{
	  if(i == j) continue;
	  jleft  = emap->lpos[cm->ndidx[hsi->firstA[j]]];
	  jright = emap->rpos[cm->ndidx[hsi->firstA[j]]];

	  if(hsi->withinAA[i][j]) {
	    if(! ((ileft <= jleft) && (iright >= jright)))
	      { printf("Crap."); }
	  }
	  else {
	    if((ileft <= jleft) && (iright >= jright))
	      printf("Crapola.");
	  }
	}
    }

  cm_ValidateHybridScanInfo(cm, hsi);

  return hsi;

  FreeEmitMap(emap);

 ERROR:
  cm_Fail("memory allocation error somewhere in cm_CreateHybridScanInfo().\n");
  return NULL;
}


/* Function: cm_AddRootToHybridScanInfo()
 * Date:     EPN, Wed Oct 31 07:36:23 2007
 *
 * Purpose:  Given a CM and hybrid scan info <hsi>, add a root <v_root_to_add> to <hsi>.
 *           
 *            
 * Returns:  eslOK; dies immediately if <v_root_to_add> is incompatible with an existing
 *           v_root in hsi.
 */
int
cm_AddRootToHybridScanInfo(CM_t *cm, HybridScanInfo_t *hsi, int v_root_to_add)
{
  int nd;
  int v;
  int i,j;
  int v_left, v_right;
  int k;
  int M;
  int lpos, rpos;

  i = hsi->startA[v_root_to_add];
  /* contract check */
  if(!(hsi->iscandA[v_root_to_add])) cm_Fail("cm_AddRootToHybridScanInfo, trying to add v_root %d, but hsi->iscandA[v] is FALSE, v is not a valid local begin point.\n", v_root_to_add);
  if(cm->ndidx[v_root_to_add] == 0)  cm_Fail("cm_AddRootToHybridScanInfo, trying to add v_root %d which is is node 0, this is illegal.\n", v_root_to_add);
  /* check if v_root_to_add we want to add is incompatible with existing v_roots, 
   * this is true if start group of existing v_root is within (not independent of) 
   * start group of v_root_to_add 
   */
  if(hsi->v_isroot[v_root_to_add]) cm_Fail("Trying to add v_root %d, but it is already a root.\n", v_root_to_add);
  for(v = 0; v < cm->M; v++) { 
    if(hsi->v_isroot[v]) {
      j = hsi->startA[v];
      if(i==j)           cm_Fail("Trying to add v_root %d, conflicts with existing v_root %d, both in start group: %d.\n", v_root_to_add, v, i);
      if(hsi->withinAA[i][j]) cm_Fail("Trying to add v_root %d in start group %d, conflicts with existing v_root %d of start group %d, because groups %d and %d are not independent.\n", v_root_to_add, i, v, j, i, j);
    }
  }

  /* if we get here, v_root_to_add does not conflict with existing v_roots in hsi, so we add it */
  CMEmitMap_t *emap;
  emap = CreateEmitMap(cm);
  nd = cm->ndidx[v_root_to_add];
  
  /* update k_* data */
  M    = hsi->cp9_M;
  lpos = emap->lpos[nd];
  rpos = emap->rpos[nd];
  if(lpos == 0 && rpos == M) cm_Fail("Adding v_root %d makes all positions modelled by CM, no need for hybrid search.\n", v_root_to_add);

  /* k_nxt, k_prv */
  if(lpos == 0) hsi->k_first                 = rpos + 1;
  else          hsi->k_nxt[hsi->k_prv[lpos]] = rpos + 1;

  if(rpos == M) hsi->k_last          = lpos - 1;
  else          hsi->k_prv[rpos + 1] = hsi->k_prv[lpos];

  for(k = lpos; k <= rpos; k++) { 
    hsi->k_mb[k] = MB_CM;
    hsi->k_nxt[k] = hsi->k_prv[k] = -1; 
  }

  /* k_nxtr, k_prvr */
  if(lpos == 0) hsi->k_firstr = rpos; 
  if(rpos == M) hsi->k_lastr  = lpos; 

  if(lpos > 0) assert(hsi->k_nxtr[lpos-1] == lpos);
  if(lpos > 0) assert(hsi->k_prvr[lpos]   == lpos-1);
  if(rpos < M) { 
    hsi->k_prvr[rpos+1] = lpos;
    hsi->k_nxtr[lpos]   = rpos+1;
  }
  else hsi->k_nxtr[lpos] = -1;

  for(k = lpos+1; k <= rpos;   k++) hsi->k_nxtr[k] = -1;
  for(k = rpos;   k >= lpos+1; k--) hsi->k_prvr[k] = -1;

		  
  /* WRITE A VALIDATE FUNCTION THAT CHECKS THAT ALL VALUES IN A HYBRIDINFO_T OBJECT ARE VALID */

  /* update v_* data */
  for(v = v_root_to_add;   v <= hsi->lastA[i]; v++)     hsi->v_mb[v]  = MB_CM;
  for(v = v_root_to_add+1; v <= hsi->lastA[i]; v++)   { hsi->v_prv[v] = v-1; hsi->v_nxt[(v-1)] = v; }
  
  /* 4 possible cases, we're adding first (1), leftmost (2), rightmost (3), or a middle chunk 
   * of states (indices v_root_to_add..hsi->lastA[i]) to be modelled by CM 
   */
  if(hsi->n_v_roots == 0) { /* this is first v_root */
    hsi->v_first = v_root_to_add;
    hsi->v_last  = hsi->lastA[i];
  }
  else if(hsi->v_first > hsi->lastA[i]) { /* we're adding leftmost (smallest v) v_root */
    hsi->v_prv[hsi->v_first]  = hsi->lastA[i];
    hsi->v_nxt[hsi->lastA[i]] = hsi->v_first;
    hsi->v_first              = v_root_to_add;
    hsi->v_prv[v_root_to_add]  = -1; /* it's the first v modelled by the CM */
  }
  else if(hsi->v_last  < v_root_to_add)  { /* we're adding right most (largest v) v_root */
    hsi->v_nxt[hsi->v_last]    = v_root_to_add;
    hsi->v_prv[v_root_to_add]   = hsi->v_last;
    hsi->v_last                = hsi->lastA[i];
    hsi->v_nxt[hsi->lastA[i]]  = -1; /* it's the last v modelled by the CM */
  }
  else { /* we're adding v_root that is in the middle (not smallest v_root, not largest v_root) */
    /* find the first MB_CM chunks before and after the one we're adding (v_root_to_add..lastA[i]) */
    v_left = v_right = -1;
    for(v = 0;       v <  v_root_to_add;  v++) if(hsi->v_mb[v] == MB_CM) v_left  = v;
    for(v = cm->M-1; v >  hsi->lastA[i]; v--) if(hsi->v_mb[v] == MB_CM) v_right = v;
    assert(v_left  != -1);
    assert(v_right != -1);
    hsi->v_prv[v_right]       = hsi->lastA[i]; 
    hsi->v_nxt[hsi->lastA[i]] = v_right; 
    hsi->v_nxt[v_left]        = v_root_to_add;
    hsi->v_prv[v_root_to_add]  = v_left;
  }
  hsi->n_v_roots++;
  hsi->v_isroot[v_root_to_add] = TRUE;
  printf("adding vroot: %4d, OLD speedup vs full CM:  %10.6f\n", v_root_to_add, hsi->full_cm_ncalcs / hsi->hybrid_ncalcs);
  printf("adding vroot: %4d, OLD speedup vs full CP9: %10.6f\n", v_root_to_add, hsi->full_cp9_ncalcs / hsi->hybrid_ncalcs);
  /* update predicted number of millions of dp calcs per residue for a hybrid scan */
  hsi->hybrid_ncalcs -= hsi->cp9_vcalcs[v_root_to_add];
  hsi->hybrid_ncalcs += hsi->cm_vcalcs[v_root_to_add];
  printf("added  vroot: %4d, NEW speedup vs full CM:  %10.6f\n", v_root_to_add, hsi->full_cm_ncalcs / hsi->hybrid_ncalcs);
  printf("added  vroot: %4d, NEW speedup vs full CP9: %10.6f\n", v_root_to_add, hsi->full_cp9_ncalcs / hsi->hybrid_ncalcs);

  cm_ValidateHybridScanInfo(cm, hsi);
  FreeEmitMap(emap);

  return eslOK;
}


/* Function: cm_ValidateHybridScanInfo()
 * Date:     EPN, Thu Nov  1 07:55:38 2007
 *
 * Purpose:  Given a CM and hybrid scan info <hsi>, validate it, by making sure it's 
 *           data is consistent. (this is likely not an exhaustive check of ALL the data in hsi).
 *            
 * Returns:  eslOK; dies immediately if <hsi> is invalid.
 */
int
cm_ValidateHybridScanInfo(CM_t *cm, HybridScanInfo_t *hsi)

{
  int nd;
  int v, v2;
  int k;
  int M;
  int lpos, rpos;
  int nfound = 0;
  int i, j;

  if(hsi->cm_M != cm->M)       cm_Fail("cm_ValidateHybridScanInfo, hsi->cm_M: (%d) != cm->M (%d)\n", hsi->cm_M, cm->M);
  if(cm->cp9 == NULL)          cm_Fail("cm_ValidateHybridScanInfo, cm->cp9 is NULL");
  if(hsi->cp9_M != cm->cp9->M) cm_Fail("cm_ValidateHybridScanInfo, hsi->cp9_M: (%d) != cm->cp9->M (%d)\n", hsi->cp9_M, cm->cp9->M);
  if(hsi->cp9_M != cm->clen)   cm_Fail("cm_ValidateHybridScanInfo, hsi->cp9_M: (%d) != cm->clen (%d)\n", hsi->cp9_M, cm->clen);

  CMEmitMap_t *emap;
  emap = CreateEmitMap(cm);

  for (v = 0; v < cm->M; v++) { 
    if((hsi->v_mb[v] != MB_CM) && (hsi->v_mb[v] != MB_CP9)) cm_Fail("cm_ValidateHybridScanInfo, hsi->v_mb[v=%d] not MB_CP9(%d) nor MB_CM(%d), but %d.\n", v, MB_CP9, MB_CM, hsi->v_mb[v]);
    if(hsi->v_mb[v] == MB_CM) { 
      nd = cm->ndidx[v];
      lpos = emap->lpos[nd];
      rpos = emap->rpos[nd];
      for(k = lpos; k <= rpos; k++)  
	if(hsi->k_mb[k] != MB_CM) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_mb[k=%d] MB_CP9, but within v=%d subtree (hsi->v_mb[%d] is MB_CM)\n", k, v, v);
    }
  }

  M = hsi->cp9_M;
  for(k = 0; k <= M; k++) if(hsi->k_mb[k] != MB_CM && hsi->k_mb[k] != MB_CP9) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_mb[k=%d] not MB_CP9(%d) nor MB_CM(%d), but %d.\n", k, MB_CP9, MB_CM, hsi->k_mb[k]);
  for(k = 0; k < M; k++) { 
    if(hsi->k_nxt[k] != -1)
      if(k  != hsi->k_prv[hsi->k_nxt[k]])  cm_Fail("cm_ValidateHybridScanInfo, k: %d != hsi->k_prv[hsi->k_nxt[k]]: %d\n", k, hsi->k_prv[hsi->k_nxt[k]]);
    if(hsi->k_nxtr[k] != -1)
      if(k != hsi->k_prvr[hsi->k_nxtr[k]]) cm_Fail("cm_ValidateHybridScanInfo, k: %d != hsi->k_prvr[hsi->k_nxtr[k]]: %d\n", k, hsi->k_prvr[hsi->k_nxtr[k]]);
  }
  for(k = 1; k <= M; k++) { 
    if(hsi->k_prv[k] != -1)  
      if(k != hsi->k_nxt[hsi->k_prv[k]])  cm_Fail("cm_ValidateHybridScanInfo, hsi->k_prv[k]: %d != hsi->k_prv[hsi->k_nxt[k]]: %d\n", k, hsi->k_prv[hsi->k_nxt[k]]);
    if(hsi->k_prvr[k] != -1) 
      if(k != hsi->k_nxtr[hsi->k_prvr[k]]) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_prvr[k]: %d != hsi->k_nxtr[hsi->k_nxtr[k]]: %d\n", k, hsi->k_prvr[hsi->k_nxtr[k]]);
  }
  if(hsi->k_nxt[M]  != -1) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_nxt[M]: %d != -1\n", hsi->k_nxt[M]);
  if(hsi->k_nxtr[M] != -1) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_nxtr[M]: %d != -1\n", hsi->k_nxtr[M]);
  if(hsi->k_prv[0]  != -1) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_prv[0]: %d != -1\n", hsi->k_prv[0]);
  if(hsi->k_prvr[0] != -1) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_prvr[0]: %d != -1\n", hsi->k_prvr[0]);
  
  if(hsi->k_first < 0 || hsi->k_first > M) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_first out of 0..M=%d range (%d)\n", M, hsi->k_first);
  if(hsi->k_last  < 0 || hsi->k_last > M)  cm_Fail("cm_ValidateHybridScanInfo, hsi->k_last out of 0..M=%d range (%d)\n", M, hsi->k_last);

  if(hsi->k_firstr < 0 || hsi->k_firstr > M) cm_Fail("cm_ValidateHybridScanInfo, hsi->k_firstr out of 0..M=%d range (%d)\n", M, hsi->k_firstr);
  if(hsi->k_lastr  < 0 || hsi->k_lastr > M)  cm_Fail("cm_ValidateHybridScanInfo, hsi->k_lastr out of 0..M=%d range (%d)\n", M, hsi->k_lastr);

  if(hsi->n_v_roots == 0) { 
    if(hsi->v_first != -1) cm_Fail("cm_ValidateHybridScanInfo, hsi->v_first not -1, but hsi->n_v_roots is 0\n", hsi->v_first);
    if(hsi->v_last != -1)  cm_Fail("cm_ValidateHybridScanInfo, hsi->v_last not -1, but hsi->n_v_roots is 0\n", hsi->v_last);
    for(v = 0; v < cm->M; v++) if(hsi->v_mb[v] == MB_CM) cm_Fail("cm_ValidateHybridScanInfo, no v_roots but hsi->v_mb[v:%d] is MB_CM\n", v);
    for(v = 0; v < cm->M; v++) if(hsi->v_isroot[v])      cm_Fail("cm_ValidateHybridScanInfo, no v_roots but hsi->v_isroot[v:%d] is TRUE\n", v);
    for(k = 0; k <= M; k++)    if(hsi->k_mb[k] == MB_CM) cm_Fail("cm_ValidateHybridScanInfo, no v_roots but hsi->k_mb[k:%d] is MB_CM\n", k);

  }
  else if(hsi->n_v_roots > 0) { 
    if(hsi->v_first < 0 || hsi->k_first > cm->M) cm_Fail("cm_ValidateHybridScanInfo, hsi->v_first out of 0..M=%d range (%d)\n", cm->M, hsi->v_first);
    if(hsi->v_last  < 0 || hsi->k_last > cm->M)  cm_Fail("cm_ValidateHybridScanInfo, hsi->v_last out of 0..M=%d range (%d)\n", cm->M, hsi->v_last);
    for(v = 0; v < cm->M; v++) { 
      if(hsi->v_isroot[v]) {
	ESL_DPRINTF1(("v: %d is a sub CM root\n", v));
	nfound++;
	i = hsi->startA[v];
	for(v2 = 0; v2 < cm->M; v2++) { 
	  if(v != v2 && hsi->v_isroot[v2]) { 
	    j = hsi->startA[v2];
	    if(hsi->withinAA[i][j]) cm_Fail("cm_ValidateHybridScanInfo, v: %d v2: %d both roots, but v2 is within v's subtree (according to hsi->withinAA)\n", v, v2);
	  /* this assumes withinAA is correct */
	  }
	}
	for(v2 = v; v2 <= hsi->lastA[hsi->startA[v]]; v2++) 
	  if(hsi->v_mb[v2] != MB_CM) cm_Fail("cm_ValidateHybridScanInfo, hsi->v_mb[v=%d] MB_CP9, but it is within root v=%d subtree\n", v2, v);
      }
    }
  }
  if(nfound != hsi->n_v_roots) cm_Fail("cm_ValidateHybridScanInfo, hsi->n_v_roots: %d != number of TRUE values in hsi->v_isroot array: %d\n", hsi->n_v_roots, nfound);

  FreeEmitMap(emap);
  return eslOK;
}

/***********************************************************************
 * Function: predict_hybrid_speedups
 * 
 */
int
predict_xsub(CM_t *cm, float *cm_vcalcs, float *cm_expsc, float *cp9_expsc)
{
  int status;
  CMEmitMap_t *emap;
  double *cp9_vcalcs;
  double *xsub;
  int cp9_ntrans;
  int v, vx;
  int nd, nd_clen;
  double cp9_filter_calcs;
  double cm_filter_calcs;
  double cp9_survivor_calcs;
  double cm_survivor_calcs;

  ESL_ALLOC(xsub,       sizeof(double) * (cm->M));
  ESL_ALLOC(cp9_vcalcs, sizeof(double) * (cm->M));

  emap = CreateEmitMap(cm);
  cp9_ntrans = NHMMSTATETYPES * NHMMSTATETYPES; /* 3*3 = 9 transitions in global mode */
  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) cp9_ntrans++;
  if(cm->cp9->flags & CPLAN9_LOCAL_END)   cp9_ntrans++;
  if(cm->cp9->flags & CPLAN9_EL)          cp9_ntrans++;

  for(nd = 0; nd < cm->nodes; nd++) { 
    nd_clen   = emap->rpos[nd] - emap->lpos[nd] + 1;
    vx = (nd < (cm->nodes-1)) ? cm->nodemap[nd+1]-1 : cm->M-1;
    for(v = cm->nodemap[nd]; v <= vx; v++) { 
      cp9_vcalcs[v] = (cp9_ntrans * nd_clen) / 1000000.; /* millions of calcs per residue */
    }
  }

  cp9_filter_calcs   = cp9_vcalcs[0];
  cp9_survivor_calcs = cm_vcalcs[0] / sreEXP2(cp9_expsc[0]);

  for(v = 0; v < cm->M; v++) { 
    cm_filter_calcs  = cp9_filter_calcs - cp9_vcalcs[v] + cm_vcalcs[v];
    cm_survivor_calcs = cm_vcalcs[0] / sreEXP2(cp9_expsc[0] - cp9_expsc[v] + cm_expsc[v]);

    xsub[v] = (cp9_filter_calcs + cp9_survivor_calcs) / (cm_filter_calcs + cm_survivor_calcs);
    printf("xsub[v:%4d] %9.5f (cp9 f: %9.5f s: %12.10f) (cm f: %9.5f s: %12.10f)\n", v, xsub[v], cp9_filter_calcs, cp9_survivor_calcs, cm_filter_calcs, cm_survivor_calcs);
  }
  return eslOK;

 ERROR: 
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}

/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_HYBRIDSEARCH_BENCHMARK
/* gcc -o benchmark-hybridsearch -g -O2 -I. -L. -I../easel -L../easel -DIMPL_HYBRIDSEARCH_BENCHMARK cm_cp9_hybridsearch.c -linfernal -leasel -lm
 * ./benchmark-hybridsearch <cmfile> 
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
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "33", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,  "10000", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search in local mode [default: glocal]", 0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search w/Viterbi also", 0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search w/CYK also", 0 },
  { "--beta",    eslARG_REAL,  "1e-7", NULL, "x>0", NULL,  NULL, NULL, "set tail loss prob for hybrid scanning QDB to <x>", 5 },
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
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float          sc;
  char           *cmfile = esl_opt_GetArg(go, 1);
  CMFILE         *cmfp;	/* open input CM file stream */
  float          *vcalcs;
  double         hsi_beta;
  char           errbuf[cmERRBUFSIZE];

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);
  
  if(esl_opt_GetBoolean(go, "-l")) 
    cm->config_opts  |= CM_CONFIG_LOCAL;
  cm->config_opts |= CM_CONFIG_QDB;
  ConfigCM(cm, NULL, NULL);
  
  if((status = cm_CountSearchDPCalcs(cm, NULL, 1000, cm->dmin, cm->dmax, cm->W, &vcalcs, NULL)) != eslOK) cm_Fail("Error counting search dp calcs.");
  
  HybridScanInfo_t *hsi;
  hsi_beta = esl_opt_GetReal(go, "--beta");
  hsi = cm_CreateHybridScanInfo(cm, hsi_beta, vcalcs[0]);
  if(esl_opt_GetBoolean(go, "-v")) { 
    cm_CreateScanMatrixForCM(cm, TRUE, TRUE); /* impt to do this after QDBs set up in ConfigCM() */
  }
     
  /* for se.cm 
     cm_AddRootToHybridScanInfo(cm, hsi, 14);
     printf("added 14 to hybrid scan info\n");
  */
  /* for 5.cm (trna) */
  cm_AddRootToHybridScanInfo(cm, hsi, 69);
  printf("added  69 to hybrid scan info\n");
  cm_AddRootToHybridScanInfo(cm, hsi, 125);
  printf("added 125 to hybrid scan info\n");
  
  /*float *cm_expsc;
    float *cp9_expsc;
    cm_CalcExpSc(cm, &cm_expsc, &cp9_expsc);
    cm_CountSearchDPCalcs(cm, 1000, cm->dmin, cm->dmax, cm->W, &vcalcs);
    predict_xsub(cm, vcalcs, cm_expsc, cp9_expsc);*/


  for (i = 0; i < N; i++) {
    esl_rnd_xfIID(r, cm->null, abc->K, L, dsq);
    
    esl_stopwatch_Start(w);
    if((status = cm_cp9_HybridScan(cm, errbuf, cm->cp9_mx, dsq, hsi, 1, L, hsi->W, 0., 
				   NULL, NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
    
    printf("%4d %-30s %10.4f bits ", (i+1), "cm_cp9_HybridScan(): ", sc);
    esl_stopwatch_Stop(w);
    esl_stopwatch_Display(stdout, w, " CPU time: ");
    
    if(esl_opt_GetBoolean(go, "-v")) { 
      esl_stopwatch_Start(w);
      if((status = cp9_FastViterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL,
				   TRUE, FALSE, FALSE, NULL, NULL, NULL,
				   &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "cm_FastViterbi(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
    }
    
    if(esl_opt_GetBoolean(go, "-c")) { 
      esl_stopwatch_Start(w);
      if((status = FastCYKScan(cm, errbuf, dsq, 1, L, cm->W, 0., NULL, NULL, &sc)) != eslOK) cm_Fail(errbuf);
      printf("%4d %-30s %10.4f bits ", (i+1), "cm_FastCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
    }
    
  }
  
  cm_FreeHybridScanInfo(hsi, cm);
  FreeCM(cm);
  free(dsq);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IMPL_HYBRIDSEARCH_BENCHMARK*/

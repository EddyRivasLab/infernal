/* cm_fastsearch.c
 * DP functions for CYK and Inside CM similarity search, includes
 * fast (optimized) and reference versions. 
 * 
 * All CYK/Inside scanning functions were rewritten between
 * versions 0.81 and 1.0 Here's a list of the 1.0 functions
 * and their 0.81 analogs. All the 1.0 functions listed are in
 * this file (cm_fastsearch.c).
 *
 * 1.0 fast version   1.0 slow version   0.81 version            
 * ----------------   ----------------   -------------
 * FastCYKScan()      RefCYKScan()       scancyk.c:CYKScan()
 *                                       bandcyk.c:CYKBandedScan()
 * FastIInsideScan()  RefIInsideScan()   scaninside.c:InsideScan()
 *                                       scaninside.c:InsideBandedScan()
 * FastFInsideScan()  RefFInsideScan()   NONE
 *
 * The 1.0 functions can be run with QDB on or off, which 
 * is implicit in the cm->si ScanInfo_t data structure,
 * which includes min/max d values for each state.
 *
 * EPN, Wed Sep 12 16:53:32 2007
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

#define TSC(s,k)   (tsc[(v) * MAXCONNECT + (s)])
#define AMX(j,v,d) (alphap[(j * cm->M * (W+1)) + ((v) * (W+1) + d)])

/* Function: FastCYKScan()
 * Date:     EPN, Wed Sep 12 16:55:28 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           an optimized CYK scanning algorithm. Query-dependent 
 *           bands are used or not used as specified in ScanInfo_t <cm->si>.
 *
 * Args:     cm              - the covariance model, must have valid scaninfo (si)
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 * 
 * Note:     This function is heavily synchronized with FastFInsideScan() and FastIInsideScan(),
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
FastCYKScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
	    search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  float    *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  float   **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  float ***alpha      = si->falpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  float ***alpha_begl = si->falpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *bestr       = si->bestr;       /* [0..d..W] best root state (for local begins or 0) for this d */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  float **esc_vAA     = cm->oesc;        /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */

  /* Contract check */
  if(! cm->flags & CMH_BITS)             cm_Fail("ERROR in FastCYKScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)         cm_Fail("ERROR in FastCYKScan, cm->si is not valid.\n");
  if(j0 < i0)                            cm_Fail("ERROR in FastCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                        cm_Fail("ERROR in FastCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in FastCYKScan, CM_SEARCH_INSIDE flag raised");
  if(! (cm->si->flags & cmSI_HAS_FLOAT)) cm_Fail("ERROR in FastCYKScan, ScanInfo's cmSI_HAS_FLOAT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in FastCYKScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* precalculate the initial scores for all cells */
  init_scAA = FCalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  float const *esc_v = esc_vAA[v]; 
	  float const *tsc_v = cm->tsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  sd   = StateDelta(cm->sttype[v]);
	  cnum = cm->cnum[v];
	  dn   = dnA[v];
	  dx   = dxA[v];
	  /* if we emit right, precalc score of emitting res j from state v */
	  float esc_j = IMPOSSIBLE;
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

	    float const *arow0;
	    float const *arow1;
	    float const *arow2;
	    float const *arow3;
	    float const *arow4;
	    float const *arow5;

	    /* Note: order of cnum cases in switch and cases in each
	     * nested emitmode switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
	      arow5 = (float * const) alpha[jp_y][y+5];
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
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

	    float const *arow0;
	    float const *arow1;
	    float const *arow2;
	    float const *arow3;
	    float const *arow4;
	    float const *arow5;

	    /* Note: order of cnum cases in switch and cases in each
	     * nested emitmode switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = ESL_MAX(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ESL_MAX(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;

	    case 6: /* necessarily 2 inserts */
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
	      arow5 = (float * const) alpha[jp_y][y+5];
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = ESL_MAX(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = ESL_MAX(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = ESL_MAX(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;

	    case 5: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
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
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
	  }
	  /* if(cm->stid[v] != BEGL_S)
	     for (d = dn; d <= dx; d++) { printf("alpha[j:%4d][v:%4d][d:%4d]: %.5f\n", j, v, d, alpha[jp_v][v][d]); }*/
	} /*loop over decks v>0 */
      
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
	bestr[d] = 0;	/* root of the traceback = root state 0 */
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ESL_MAX (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(NOT_IMPOSSIBLE(cm->beginsc[y])) {
	    if(cm->stid[y] == BEGL_S)
	      {
		jp_y = j % (W+1);
		for (d = dnA[y]; d <= dxA[y]; d++) {
		  /* Is this more efficient:? 
		     bestr[d]          = (alpha[jp_v][0][d] > (alpha_begl[jp_y][y][d] + cm->beginsc[y])) ? bestr[d] : y;
		     alpha[jp_v][0][d] = ESL_MAX(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->beginsc[y]); */
		  if(alpha[jp_v][0][d] < (alpha_begl[jp_y][y][d] + cm->beginsc[y])) {
		    alpha[jp_v][0][d] = alpha_begl[jp_y][y][d] + cm->beginsc[y];
		    bestr[d] = y;
		  }
		}
	      }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		{
		  /* Is this more efficient:? 
		     bestr[d]          = (alpha[jp_v][0][d] > (alpha[jp_y][y][d] + cm->beginsc[y])) ? bestr[d] : y;
		     alpha[jp_v][0][d] = ESL_MAX(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->beginsc[y]); */
		  if(alpha[jp_v][0][d] < (alpha[jp_y][y][d] + cm->beginsc[y])) {
		    alpha[jp_v][0][d] = alpha[jp_y][y][d] + cm->beginsc[y];
		    bestr[d] = y;
		  }
		}
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, alpha[jp_v][0][d]);
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateFloatGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, FALSE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, TRUE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateFloatGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(sc_v);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc         = vsc;
  
  ESL_DPRINTF1(("FastCYKScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/* Function: FastIInsideScan()
 * Date:     EPN, Tue Nov  6 05:42:44 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           an optimized Inside scanning algorithm that uses integer scores. 
 *           Query-dependent bands are used or not used as specified in 
 *           ScanInfo_t <si>.
 *
 * Args:     cm              - the covariance model
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 * 
 * Note:     This function is heavily synchronized with FastCYKScan() and FastFInsideScan(),
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best Inside hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
FastIInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
		search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  int      *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  int     **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  int   ***alpha      = si->ialpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  int   ***alpha_begl = si->ialpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  int   **esc_vAA     = cm->ioesc;       /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */
  /* Contract check */
  if(! cm->flags & CMH_BITS)                 cm_Fail("ERROR in FastIInsideScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)             cm_Fail("ERROR in FastIInsideScan, cm->si is not valid.\n");
  if(j0 < i0)                                cm_Fail("ERROR in FastIInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                            cm_Fail("ERROR in FastIInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in FastIInsideScan, CM_SEARCH_INSIDE flag not raised");
  if(! (si->flags & cmSI_HAS_INT))           cm_Fail("ERROR in FastIInsideScan, ScanInfo's cmSI_HAS_INT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in FastIInsideScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_ISet(sc_v, (W+1), -INFTY);

  /* precalculate the initial scores for all cells */
  init_scAA = ICalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      int sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  int const *esc_v = esc_vAA[v]; 
	  int const *tsc_v = cm->itsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
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

	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for B st */
	      for (k = kmin; k <= kmax; k++) 
		sc = ILogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ILogsum (sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
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
		sc = ILogsum(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
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
		sc = ILogsum(arow5[dp_y] + tsc_v[5],
			     init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow4[dp_y] + tsc_v[4]);		
		sc = ILogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		/* add in emission score, if any */
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
		sc = ILogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
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
		sc = ILogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		

		/* add in emission score, if any */
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
		sc = ILogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
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
		sc_v[d] = ILogsum(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
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
		sc_v[d] = ILogsum(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc_v[d] = ILogsum(sc_v[d], arow4[dp_y] + tsc_v[4]);		
		sc_v[d] = ILogsum(sc_v[d], arow3[dp_y] + tsc_v[3]);		
		sc_v[d] = ILogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;

	    case 4: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      arow3 = (int * const) alpha[jp_y][y+3];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = ILogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = ILogsum(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;

	    case 5: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      arow3 = (int * const) alpha[jp_y][y+3];
	      arow4 = (int * const) alpha[jp_y][y+4];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = ILogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc_v[d] = ILogsum(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for (d = dn; d <= dx; d++, dp_y++) */
	      break;

	    case 2: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = ILogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
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
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha[jp_v][v][d]));
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha_begl[jp_v][v][d]));
	  }
	  /*if(cm->stid[v] != BEGL_S)
	    for (d = dn; d <= dx; d++) { printf("alpha[j:%4d][v:%4d][d:%4d]: %.5f\n", j, v, d, alpha[jp_v][v][d]); }*/
	} /*loop over decks v>0 */
      
      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      int const *tsc_v = cm->itsc[0];
      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = dnA[0]; d <= dxA[0]; d++) {
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(-INFTY, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ILogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, Scorify(alpha[jp_v][0][d]));
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateIntGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, TRUE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, FALSE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateIntGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(sc_v);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("FastIInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/* Function: XFastIInsideScan()
 * Date:     EPN, Tue Nov  6 05:42:44 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           an optimized Inside scanning algorithm that uses integer scores. 
 *           Query-dependent bands are used or not used as specified in 
 *           ScanInfo_t <si>.
 *
 * Args:     cm              - the covariance model
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 * 
 * Note:     This function is heavily synchronized with FastCYKScan() and FastFInsideScan(),
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best Inside hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
XFastIInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
		search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  int      *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  int     **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  int   ***alpha      = si->ialpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  int   ***alpha_begl = si->ialpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  int   **esc_vAA     = cm->ioesc;       /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */
  /* Contract check */
  if(! cm->flags & CMH_BITS)                 cm_Fail("ERROR in XFastIInsideScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)             cm_Fail("ERROR in XFastIInsideScan, cm->si is not valid.\n");
  if(j0 < i0)                                cm_Fail("ERROR in XFastIInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                            cm_Fail("ERROR in XFastIInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in XFastIInsideScan, CM_SEARCH_INSIDE flag not raised");
  if(! (si->flags & cmSI_HAS_INT))           cm_Fail("ERROR in XFastIInsideScan, ScanInfo's cmSI_HAS_INT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in XFastIInsideScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(int) * (W+1)));
  esl_vec_ISet(sc_v, (W+1), -INFTY);

  /* precalculate the initial scores for all cells */
  init_scAA = ICalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      int sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  int const *esc_v = esc_vAA[v]; 
	  int const *tsc_v = cm->itsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
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

	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for B st */
	      for (k = kmin; k <= kmax; k++) 
		sc = ILogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ILogsum (sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
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
		sc = ILogsum(init_scAA[v][dp_y], 
			     arow0[dp_y] + tsc_v[0]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		
		/* add in emission score, if any */
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
		sc = ILogsum(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow4[dp_y] + tsc_v[4]);		
		sc = ILogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		/* add in emission score, if any */
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
		sc = ILogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
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
		sc = ILogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		

		/* add in emission score, if any */
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
		sc = ILogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
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

	    for (d = dn; d <= dx; d++, dp_y++) {
	      sc_v[d] = ILogsum(init_scAA[v][dp_y],
				alpha[jp_y][y][dp_y] + tsc_v[0]);
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		sc_v[d] = ILogsum(sc_v[d], alpha[jp_y][y+yoffset][dp_y] + tsc_v[yoffset]);
	    }
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
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha[jp_v][v][d]));
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha_begl[jp_v][v][d]));
	  }
	  /* if(cm->stid[v] != BEGL_S)
	     for (d = dn; d <= dx; d++) { printf("alpha[j:%4d][v:%4d][d:%4d]: %.5f\n", j, v, d, alpha[jp_v][v][d]); }*/
	} /*loop over decks v>0 */
      
      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      int const *tsc_v = cm->itsc[0];
      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = dnA[0]; d <= dxA[0]; d++) {
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(-INFTY, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ILogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, Scorify(alpha[jp_v][0][d]));
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateIntGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, TRUE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, FALSE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateIntGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(sc_v);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("XFastIInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/* Function: X2FastIInsideScan()
 * Date:     EPN, Tue Nov  6 05:42:44 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           an optimized Inside scanning algorithm that uses integer scores. 
 *           Query-dependent bands are used or not used as specified in 
 *           ScanInfo_t <si>.
 *
 * Args:     cm              - the covariance model
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 * 
 * Note:     This function is heavily synchronized with FastCYKScan() and FastFInsideScan(),
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best Inside hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
X2FastIInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
		search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  int      *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  int     **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  int   ***alpha      = si->ialpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  int   ***alpha_begl = si->ialpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  int   **esc_vAA     = cm->ioesc;       /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */
  /* Contract check */
  if(! cm->flags & CMH_BITS)                 cm_Fail("ERROR in FastIInsideScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)             cm_Fail("ERROR in FastIInsideScan, cm->si is not valid.\n");
  if(j0 < i0)                                cm_Fail("ERROR in FastIInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                            cm_Fail("ERROR in FastIInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in FastIInsideScan, CM_SEARCH_INSIDE flag not raised");
  if(! (si->flags & cmSI_HAS_INT))           cm_Fail("ERROR in FastIInsideScan, ScanInfo's cmSI_HAS_INT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in FastIInsideScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_ISet(sc_v, (W+1), -INFTY);

  /* precalculate the initial scores for all cells */
  init_scAA = ICalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      int sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  int const *esc_v = esc_vAA[v]; 
	  int const *tsc_v = cm->itsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
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
		sc = ILogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ILogsum (sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
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
	    case 2: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = init_scAA[v][dp_y];
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		
		/* add in emission score, if any */
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

	    case 3: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = init_scAA[v][dp_y];
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		
		/* add in emission score, if any */
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
		sc = init_scAA[v][dp_y];
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow3[dp_y] + tsc_v[3]);		
		
		/* add in emission score, if any */
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
		sc = init_scAA[v][dp_y];
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, arow4[dp_y] + tsc_v[4]);		
		
		/* add in emission score, if any */
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

	    case 6: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      arow3 = (int * const) alpha[jp_y][y+3];
	      arow4 = (int * const) alpha[jp_y][y+4];
	      arow5 = (int * const) alpha[jp_y][y+5];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = init_scAA[v][dp_y];
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, arow4[dp_y] + tsc_v[4]);		
		sc = ILogsum(sc, arow5[dp_y] + tsc_v[5]);		
		
		/* add in emission score, if any */
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
	    case 2:
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = init_scAA[v][dp_y];
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
	      }
	      break;

	    case 3:
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = init_scAA[v][dp_y];
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
	      }
	      break;

	    case 4:
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      arow3 = (int * const) alpha[jp_y][y+3];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = init_scAA[v][dp_y];
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = ILogsum(sc_v[d], arow3[dp_y] + tsc_v[3]);		
	      }
	      break;

	    case 5:
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      arow3 = (int * const) alpha[jp_y][y+3];
	      arow4 = (int * const) alpha[jp_y][y+4];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = init_scAA[v][dp_y];
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = ILogsum(sc_v[d], arow3[dp_y] + tsc_v[3]);		
		sc_v[d] = ILogsum(sc_v[d], arow4[dp_y] + tsc_v[4]);		
	      }
	      break;

	    case 6:
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
	      arow3 = (int * const) alpha[jp_y][y+3];
	      arow4 = (int * const) alpha[jp_y][y+4];
	      arow5 = (int * const) alpha[jp_y][y+5];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = init_scAA[v][dp_y];
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
		sc_v[d] = ILogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = ILogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = ILogsum(sc_v[d], arow3[dp_y] + tsc_v[3]);		
		sc_v[d] = ILogsum(sc_v[d], arow4[dp_y] + tsc_v[4]);		
		sc_v[d] = ILogsum(sc_v[d], arow5[dp_y] + tsc_v[5]);		
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
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha[jp_v][v][d]));
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha_begl[jp_v][v][d]));
	  }
	  /*if(cm->stid[v] != BEGL_S)
	    for (d = dn; d <= dx; d++) { printf("alpha[j:%4d][v:%4d][d:%4d]: %.5f\n", j, v, d, alpha[jp_v][v][d]); }*/
	} /*loop over decks v>0 */
      
      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      int const *tsc_v = cm->itsc[0];
      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = dnA[0]; d <= dxA[0]; d++) {
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(-INFTY, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ILogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, Scorify(alpha[jp_v][0][d]));
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateIntGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, TRUE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, FALSE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateIntGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(sc_v);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("XFastIInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/* Function: FastFInsideScan()
 * Date:     EPN, Wed Sep 12 16:55:28 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           an optimized Inside scanning algorithm that uses float scores. 
 *           Query-dependent bands are used or not used as specified in 
 *           ScanInfo_t <si>.
 *
 * Args:     cm              - the covariance model
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 * 
 * Note:     This function is heavily synchronized with FastCYKScan() and FastIInsideScan(),
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best Inside hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
FastFInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
		search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  float    *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  float   **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  float ***alpha      = si->falpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  float ***alpha_begl = si->falpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  float **esc_vAA     = cm->oesc;        /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */
  /* Contract check */
  if(! cm->flags & CMH_BITS)                cm_Fail("ERROR in FastFInsideScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)            cm_Fail("ERROR in FastFInsideScan, cm->si is not valid.\n");
  if(j0 < i0)                               cm_Fail("ERROR in FastFInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                           cm_Fail("ERROR in FastFInsideScan, dsq is NULL\n");
  if(!(cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in FastFInsideScan, CM_SEARCH_INSIDE flag not raised");
  if(! (cm->si->flags & cmSI_HAS_FLOAT))    cm_Fail("ERROR in FastFInsideScan, ScanInfo's cmSI_HAS_FLOAT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in FastFInsideScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* precalculate the initial scores for all cells */
  init_scAA = FCalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  float const *esc_v = esc_vAA[v]; 
	  float const *tsc_v = cm->tsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  sd   = StateDelta(cm->sttype[v]);
	  cnum = cm->cnum[v];
	  dn   = dnA[v];
	  dx   = dxA[v];
	  /* if we emit right, precalc score of emitting res j from state v */
	  float esc_j = IMPOSSIBLE;
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
		sc = FLogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = FLogsum (sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
	    y    = cm->cfirst[v];
	    dp_y = dn - sd; /* initial dp_y, we increment it at end of 'for(d = ...' loop */
	    i    = j-dn+1;  /* initial i,    we decrement it when we access it, inside each possible case of the switch (cnum) below */

	    float const *arow0;
	    float const *arow1;
	    float const *arow2;
	    float const *arow3;
	    float const *arow4;
	    float const *arow5;

	    /* Note: order of cnum cases in switch and cases in each
	     * nested emitmode switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
	      arow5 = (float * const) alpha[jp_y][y+5];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow4[dp_y] + tsc_v[4]);		
		sc = FLogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = FLogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		/* add in emission score, if any */
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		

		/* add in emission score, if any */
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
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
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

	    float const *arow0;
	    float const *arow1;
	    float const *arow2;
	    float const *arow3;
	    float const *arow4;
	    float const *arow5;

	    /* Note: order of cnum cases in switch and cases in each
	     * nested emitmode switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = FLogsum(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = FLogsum(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = FLogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;

	    case 6: /* necessarily 2 inserts */
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
	      arow5 = (float * const) alpha[jp_y][y+5];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = FLogsum(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc_v[d] = FLogsum(sc_v[d], arow4[dp_y] + tsc_v[4]);		
		sc_v[d] = FLogsum(sc_v[d], arow3[dp_y] + tsc_v[3]);		
		sc_v[d] = FLogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = FLogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = FLogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;

	    case 4: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = FLogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = FLogsum(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = FLogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = FLogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;

	    case 5: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
	      arow3 = (float * const) alpha[jp_y][y+3];
	      arow4 = (float * const) alpha[jp_y][y+4];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = FLogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc_v[d] = FLogsum(sc_v[d], init_scAA[v][dp_y]);
		sc_v[d] = FLogsum(sc_v[d], arow1[dp_y] + tsc_v[1]);		
		sc_v[d] = FLogsum(sc_v[d], arow2[dp_y] + tsc_v[2]);		
		sc_v[d] = FLogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      } /* end of for (d = dn; d <= dx; d++, dp_y++) */
	      break;

	    case 2: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = FLogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc_v[d] = FLogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
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
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
	  }
	  /*if(cm->stid[v] != BEGL_S)
	     for (d = dn; d <= dx; d++) { printf("alpha[j:%4d][v:%4d][d:%4d]: %.5f\n", j, v, d, alpha[jp_v][v][d]); }*/
	} /*loop over decks v>0 */
      
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
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = FLogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(NOT_IMPOSSIBLE(cm->beginsc[y])) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) 
		alpha[jp_v][0][d] = FLogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->beginsc[y]);
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) 
		alpha[jp_v][0][d] = FLogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->beginsc[y]);
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, alpha[jp_v][0][d]);
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateFloatGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, TRUE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, FALSE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateFloatGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(sc_v);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("FastFInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/* Function: RefCYKScan()
 * Date:     EPN, Wed Sep 12 16:55:28 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           a reference CYK scanning algorithm. Query-dependent 
 *           bands are used or not used as specified in ScanInfo_t <si>.
 *
 *           This function is slower, but easier to understand than the
 *           FastCYKScan() version.
 *
 * Args:     cm              - the covariance model
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *
 * Note:     This function is heavily synchronized with RefIInsideScan() and RefCYKScan()
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
RefCYKScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
	   search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  float   **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  float ***alpha      = si->falpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  float ***alpha_begl = si->falpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *bestr       = si->bestr;       /* [0..d..W] best root state (for local begins or 0) for this d */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  float **esc_vAA     = cm->oesc;        /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */

  /* Contract check */
  if(! cm->flags & CMH_BITS)             cm_Fail("ERROR in RefCYKScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)         cm_Fail("ERROR in RefCYKScan, cm->si is not valid.\n");
  if(j0 < i0)                            cm_Fail("ERROR in RefCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                        cm_Fail("ERROR in RefCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in RefCYKScan, CM_SEARCH_INSIDE flag raised");
  if(! (cm->si->flags & cmSI_HAS_FLOAT)) cm_Fail("ERROR in RefCYKScan, ScanInfo's cmSI_HAS_FLOAT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in RefCYKScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* precalculate the initial scores for all cells */
  init_scAA = FCalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  float const *esc_v = esc_vAA[v]; 
	  float const *tsc_v = cm->tsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  sd   = StateDelta(cm->sttype[v]);
	  cnum = cm->cnum[v];
	  dn   = dnA[v];
	  dx   = dxA[v];
	  /* if we emit right, precalc score of emitting res j from state v */
	  float esc_j = IMPOSSIBLE;
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
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ESL_MAX(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* ! B_st, ! BEGL_S st */
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; 
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ESL_MAX(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);

	      switch (emitmode) {
	      case EMITLEFT:
		alpha[jp_v][v][d] = sc + esc_v[dsq[i--]];
		break;
	      case EMITNONE:
		alpha[jp_v][v][d] = sc;
		break;
	      case EMITRIGHT:
		alpha[jp_v][v][d] = sc + esc_j;
		break;		
	      case EMITPAIR:
		alpha[jp_v][v][d] = sc + esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
		break;
	      } /* end of switch emitmode */
	    } /* end of for d loop */
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
	  }
	} /*loop over decks v>0 */
      
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
	bestr[d] = 0;	/* root of the traceback = root state 0 */
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ESL_MAX (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(NOT_IMPOSSIBLE(cm->beginsc[y])) {
	    if(cm->stid[y] == BEGL_S)
	      {
		jp_y = j % (W+1);
		for (d = dnA[y]; d <= dxA[y]; d++) {
		  /* Is this more efficient:? 
		     bestr[d]          = (alpha[jp_v][0][d] > (alpha_begl[jp_y][y][d] + cm->beginsc[y])) ? bestr[d] : y;
		     alpha[jp_v][0][d] = ESL_MAX(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->beginsc[y]); */
		  if(alpha[jp_v][0][d] < (alpha_begl[jp_y][y][d] + cm->beginsc[y])) {
		    alpha[jp_v][0][d] = alpha_begl[jp_y][y][d] + cm->beginsc[y];
		    bestr[d] = y;
		  }
		}
	      }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		{
		  /* Is this more efficient:? 
		     bestr[d]          = (alpha[jp_v][0][d] > (alpha[jp_y][y][d] + cm->beginsc[y])) ? bestr[d] : y;
		     alpha[jp_v][0][d] = ESL_MAX(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->beginsc[y]); */
		  if(alpha[jp_v][0][d] < (alpha[jp_y][y][d] + cm->beginsc[y])) {
		    alpha[jp_v][0][d] = alpha[jp_y][y][d] + cm->beginsc[y];
		    bestr[d] = y;
		  }
		}
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, alpha[jp_v][0][d]);
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateFloatGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, FALSE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, TRUE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateFloatGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc         = vsc;
  
  ESL_DPRINTF1(("RefCYKScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

/* Function: RefIInsideScan()
 * Date:     EPN, Tue Nov  6 06:13:35 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           a reference CYK scanning algorithm. Query-dependent 
 *           bands are used or not used as specified in ScanInfo_t <si>.
 *
 *           This function is slower, but easier to understand than the
 *           FastCYKScan() version.
 *
 * Args:     cm              - the covariance model
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *
 * Note:     This function is heavily synchronized with RefCYKScan() and RefFInsideScan(), 
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
RefIInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
	       search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  int     **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  int   ***alpha      = si->ialpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  int   ***alpha_begl = si->ialpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  int   **esc_vAA     = cm->ioesc;       /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */
  /* Contract check */
  if(! cm->flags & CMH_BITS)                 cm_Fail("ERROR in RefIInsideScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)             cm_Fail("ERROR in RefIInsideScan, cm->si is not valid.\n");
  if(j0 < i0)                                cm_Fail("ERROR in RefIInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                            cm_Fail("ERROR in RefIInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in RefIInsideScan, CM_SEARCH_INSIDE flag not raised");
  if(! (si->flags & cmSI_HAS_INT))           cm_Fail("ERROR in RefIInsideScan, ScanInfo's cmSI_HAS_INT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in RefIInsideScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* precalculate the initial scores for all cells */
  init_scAA = ICalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      int sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  int const *esc_v = esc_vAA[v]; 
	  int const *tsc_v = cm->itsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
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
		sc = ILogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* ! B_st, ! BEGL_S st */
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; 
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);

	      switch (emitmode) {
	      case EMITLEFT:
		alpha[jp_v][v][d] = sc + esc_v[dsq[i--]];
		break;
	      case EMITNONE:
		alpha[jp_v][v][d] = sc;
		break;
	      case EMITRIGHT:
		alpha[jp_v][v][d] = sc + esc_j;
		break;		
	      case EMITPAIR:
		alpha[jp_v][v][d] = sc + esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
		break;
	      } /* end of switch emitmode */
	    } /* end of for d loop */
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha[jp_v][v][d]));
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha_begl[jp_v][v][d]));
	  }
	} /*loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      int const *tsc_v = cm->itsc[0];
      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = dnA[0]; d <= dxA[0]; d++) {
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(-INFTY, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ILogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = j % (W+1);
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, Scorify(alpha[jp_v][0][d]));
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateIntGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, TRUE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, FALSE);*/
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateIntGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("RefIInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

/* Function: XRefIInsideScan()
 * Date:     EPN, Tue Nov  6 06:13:35 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           a reference CYK scanning algorithm. Query-dependent 
 *           bands are used or not used as specified in ScanInfo_t <si>.
 *
 *           This function is slower, but easier to understand than the
 *           FastCYKScan() version.
 *
 * Args:     cm              - the covariance model
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *
 * Note:     This function is heavily synchronized with RefCYKScan() and RefFInsideScan(), 
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
XRefIInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
	       search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  int      *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  int     **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */

  /* make pointers to the ScanInfo/CM data for convenience */
  int   ***alpha      = si->ialpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  int   ***alpha_begl = si->ialpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  int   **esc_vAA     = cm->ioesc;       /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */
  /* Contract check */
  if(! cm->flags & CMH_BITS)                 cm_Fail("ERROR in XRefIInsideScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)             cm_Fail("ERROR in XRefIInsideScan, cm->si is not valid.\n");
  if(j0 < i0)                                cm_Fail("ERROR in XRefIInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                            cm_Fail("ERROR in XRefIInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in XRefIInsideScan, CM_SEARCH_INSIDE flag not raised");
  if(! (si->flags & cmSI_HAS_INT))           cm_Fail("ERROR in XRefIInsideScan, ScanInfo's cmSI_HAS_INT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in XRefIInsideScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_ISet(sc_v, (W+1), -INFTY);

  /* precalculate the initial scores for all cells */
  init_scAA = ICalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      int sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  int const *esc_v = esc_vAA[v]; 
	  int const *tsc_v = cm->itsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
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
		sc = ILogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* ! B_st, ! BEGL_S st */
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    switch (emitmode) {
	      case EMITLEFT:
		for (d = dnA[v]; d <= dxA[v]; d++) {
		  sc = init_scAA[v][d-sd]; 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		  alpha[jp_v][v][d] = sc + esc_v[dsq[i--]];
		} /* end of for d loop */
		break;
	      case EMITRIGHT:
		for (d = dnA[v]; d <= dxA[v]; d++) {
		  sc = init_scAA[v][d-sd]; 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		  alpha[jp_v][v][d] = sc + esc_j;
		} /* end of for d loop */
		break;
	      case EMITNONE:
		for (d = dnA[v]; d <= dxA[v]; d++) {
		  sc = init_scAA[v][d-sd]; 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		  alpha[jp_v][v][d] = sc;
		} /* end of for d loop */
		break;
	      case EMITPAIR:
		for (d = dnA[v]; d <= dxA[v]; d++) {
		  sc = init_scAA[v][d-sd]; 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		alpha[jp_v][v][d] = sc + esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
		} /* end of for d loop */
		break;
	    } /* end of switch emitmode */
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha[jp_v][v][d]));
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha_begl[jp_v][v][d]));
	  }
	} /*loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */

      int const *tsc_v = cm->itsc[0];
      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = dnA[0]; d <= dxA[0]; d++) {
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(-INFTY, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ILogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) { 
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->ibeginsc[y]);
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->ibeginsc[y]); 
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, Scorify(alpha[jp_v][0][d]));
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateIntGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, TRUE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, FALSE);*/
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateIntGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(sc_v);
  free(init_scAA[0]);
  free(init_scAA);
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("XRefIInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

/* Function: RefFInsideScan()
 * Date:     EPN, Sun Nov  4 16:02:17 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using
 *           a reference CYK scanning algorithm. Query-dependent 
 *           bands are used or not used as specified in ScanInfo_t <si>.
 *
 *           This function is slower, but easier to understand than the
 *           FastCYKScan() version.
 *
 * Args:     cm              - the covariance model
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
 *           dsq             - the digitized sequence
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *
 * Note:     This function is heavily synchronized with RefCYKScan() and RefIInsideScan()
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
RefFInsideScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
	       search_results_t *results, float **ret_vsc)
{
  int       status;
  cm_GammaHitMx_t *gamma;       /* semi-HMM for hit resoultion */
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
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */
  ScanInfo_t *si = cm->si;      /* for convenience */
  float    *sc_v;               /* [0..d..W] temporary score vec for each d for current j & v */
  float   **init_scAA;          /* [0..v..cm->M-1][0..d..W] initial score for each v, d for all j */
  
  /* make pointers to the ScanInfo/CM data for convenience */
  float ***alpha      = si->falpha;      /* [0..j..1][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v == BEGL_S */
  float ***alpha_begl = si->falpha_begl; /* [0..j..W][0..v..cm->M-1][0..d..W] alpha DP matrix, NULL for v != BEGL_S */
  int   **dnAA        = si->dnAA;        /* [0..v..cm->M-1][0..j..W] minimum d for v, j (for j > W use [v][W]) */
  int   **dxAA        = si->dxAA;        /* [0..v..cm->M-1][0..j..W] maximum d for v, j (for j > W use [v][W]) */
  int    *dmin        = si->dmin;        /* [0..v..cm->M-1] minimum d allowed for this state */
  int    *dmax        = si->dmax;        /* [0..v..cm->M-1] maximum d allowed for this state */
  float **esc_vAA     = cm->oesc;        /* [0..v..cm->M-1][0..a..(cm->abc->Kp | cm->abc->Kp**2)] optimized emission scores for v 
					  * and all possible emissions a (including ambiguities) */

  /* Contract check */
  if(! cm->flags & CMH_BITS)                cm_Fail("ERROR in RefFInsideScan, CMH_BITS flag is not raised.\n");
  if(! cm->flags & CMH_SCANINFO)            cm_Fail("ERROR in RefFInsideScan, cm->si is not valid.\n");
  if(j0 < i0)                               cm_Fail("ERROR in RefFInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)                           cm_Fail("ERROR in RefFInsideScan, dsq is NULL\n");
  if(!(cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in RefFInsideScan, CM_SEARCH_INSIDE flag not raised");
  if(! (cm->si->flags & cmSI_HAS_FLOAT))    cm_Fail("ERROR in RefFInsideScan, ScanInfo's cmSI_HAS_FLOAT flag is not raised");

  /* determine if we're doing banded/non-banded */
  if(si->dmin != NULL && si->dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 
  if (W > si->W) cm_Fail("ERROR in RefFInsideScan, W: %d greater than si->W: %d\n", W, si->W);

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
  if(results != NULL) gamma = cm_CreateGammaHitMx(L, i0, (cm->search_opts & CM_SEARCH_CMGREEDY), cutoff);
  else                gamma = NULL;

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* precalculate the initial scores for all cells */
  init_scAA = FCalcInitDPScores(cm);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else          { dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  float const *esc_v = esc_vAA[v]; 
	  float const *tsc_v = cm->tsc[v];
	  int emitmode = Emitmode(cm->sttype[v]);

	  /* float sc; */
	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  sd   = StateDelta(cm->sttype[v]);
	  cnum = cm->cnum[v];
	  dn   = dnA[v];
	  dx   = dxA[v];
	  /* if we emit right, precalc score of emitting res j from state v */
	  float esc_j = IMPOSSIBLE;
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
		sc = FLogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; /* state delta (sd) is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = FLogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* ! B_st, ! BEGL_S st */
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    /* printf("B BEGL j: %d, v: %d\n", j, v);*/
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d-sd]; 
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		{
		  /*printf("sc: %f\n", sc);
		  printf("d: %d\n", d);
		  printf("sd: %d\n", sd);
		  printf("tsc: %d\n", tsc_v[yoffset]);
		  printf("alpha: %f\n", alpha[jp_y][y+yoffset][d - sd]);*/
		  sc = FLogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);
		}
	      switch (emitmode) {
	      case EMITLEFT:
		alpha[jp_v][v][d] = sc + esc_v[dsq[i--]];
		break;
	      case EMITNONE:
		alpha[jp_v][v][d] = sc;
		break;
	      case EMITRIGHT:
		alpha[jp_v][v][d] = sc + esc_j;
		break;		
	      case EMITPAIR:
		alpha[jp_v][v][d] = sc + esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
		break;
	      } /* end of switch emitmode */
	    } /* end of for d loop */
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) {
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
	  }
	} /*loop over decks v>0 */
      
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
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = FLogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CMH_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(NOT_IMPOSSIBLE(cm->beginsc[y])) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = FLogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->beginsc[y]);
	      }
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) {
		alpha[jp_v][0][d] = FLogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->beginsc[y]);
	      }
	    }
	  }
	}
      }
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, alpha[jp_v][0][d]);
      /* update gamma, but only if we're reporting hits to results */
      if(results != NULL) cm_UpdateFloatGammaHitMx(gamma, jp_g, alpha[jp_v][0], dnA[0], dxA[0], si->bestr, cm->sc_boost, TRUE, results);
      /* cm_DumpScanInfoAlpha(cm, si, j, i0, FALSE); */
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;

  /* If recovering hits in a non-greedy manner, do the traceback.
   * If we were greedy, they were reported in cm_UpdateFloatGammaHitMx() for each position j */
  if(results != NULL && gamma->iamgreedy == FALSE) 
    cm_TBackGammaHitMx(gamma, results, i0, j0);

  /* clean up and return */
  if(gamma != NULL) cm_FreeGammaHitMx(gamma);
  free(jp_wA);
  free(init_scAA[0]);
  free(init_scAA);
  free(sc_v);
  if (ret_vsc != NULL) *ret_vsc         = vsc;
  
  ESL_DPRINTF1(("RefFInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

/* Function: rsearch_CYKScan()
 * Date:     EPN, Thu Oct 18 05:09:13 2007 [updated to Easel, Infernal]
 *           RJK, Sun Mar 24, 2002 [STL->DCA]
 *           SRE, Mon Aug  7 13:15:37 2000 [St. Louis] 
 *                   (from inside() in smallcyk.c)
 *
 * Purpose:  Run the inside phase of a CYK alignment algorithm, on 
 *           a complete sequence of length L.
 *
 * The following changes were made from inside() in smallyck.c:
 * 1.  Removed jp, i0, j0, W, vroot, and vend because full sequence
 * 2.  Added d<=D constraint on d for loops
 * 3.  Removed shadow matrices; we don't do traceback
 * 4.  Replace alpha with gamma[j][d][v] -- makes more efficient use of cache
 * 5.  Explicitiy define End state rather than just assigning the pre-computed
 *     "end deck" 
 * 6.  Use gamma_begl_s and gamma_begr_s [v][j][d] for optimizing bifurcation
 *     states
 * 7.  Local alignment now.
 * 8.  If passed a histogram, fills in with best j for every D+1 place
 * 9.  Only reports best hit at each j to reduce reporting complexity
 * 10. Modified so minimize dereferencing for improved speed.
 * 11. Inner loops rewritten to allow vectorization.
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           cutoff    - cutoff score to report
 *           D         - maximum size of hit
 *           results   - search_results_t to fill in; if NULL, nothing
 *                       filled in
 *
 * Returns: Score of best hit overall
 *
 */
float rsearch_CYKScan (CM_t *cm, ESL_DSQ *dsq, int L, float cutoff, int D,
		       search_results_t *results) {

  int     *bestr;               /* Best root state for d at current j */
  int      v,y,z;		/* indices for states  */
  int      j,d,i,k;		/* indices in sequence dimensions */
  int      jmod2, jmin1mod2;    /* For indices into the actual j dimension */
  int      dmax;                /* D of best hit at j */
  float    sc;  	       	/* a temporary variable holding a score */
  int      yoffset;		/* y=base+offset -- counter in child 
                                   states that v can transit to */
  int      M;                   /* Stores cm->M for loop limits */
  int      cnum;                /* Stores cm->cnum[v] for loop limits */
  int      minDj;               /* Minimum cutoff for d between j and D */
  float ***gamma;               /* The main DP matrix [j][d][v] */
  float ***gamma_begl_s;        /* For BEGL_S states -- [v][i][d] */
  float ***gamma_begr_s;        /* For BEGR_S states -- [v][j][d] */ 
  float  **gamma_jmod2;         /* gamma[jmod2] */
  float  **gamma_jmin1mod2;     /* gamma[jmin1mod2] */
  float   *gammap;              /* Pointer to last dimension of gamma to use */
  float   *gamma_begl_s_p;     
  float   *gamma_begr_s_p;
  int       sc_v_size;
  float    *sc_v;               /* Vector of possible scores to maximize over */
  float    *tsc;                /* Points to cm->tsc[v] to make pointer operation simpler in loop I want to vectorize */
  float     endsc;              /* endsc for current state [v] -- set at
				   beginning of each v loop */
  float     beginsc;            /* beginsc for current state[y] */
  char      sttype;             /* Holds cm->sttype[v] */
  float     best_score = IMPOSSIBLE;     /* Best overall score to return */
  int       status;

  /* Set M */
M = cm->M;

  if (M>D+1) {
    sc_v_size = M;
  } else {
    sc_v_size = D+1;
  }
/* #ifdef INTEL_COMPILER*/
#if 1
  printf("intel compiler baby!\n");
#endif

  /* Initialize sc_v to size of M */
  ESL_ALLOC(sc_v, (sizeof(float) * sc_v_size));


  ESL_ALLOC(gamma, (sizeof(float **) * 2));
  ESL_ALLOC(gamma[0], (sizeof(float *)*2*M));
  gamma[1] = gamma[0] + M;
  ESL_ALLOC(gamma[0][0], (sizeof(float)*2*(D+1)*M));
  gamma[1][0] = gamma[0][0] + ((D+1)*M);
  for (v=1; v<M; v++) {
    gamma[0][v] = gamma[0][v-1] + (D+1);
    gamma[1][v] = gamma[1][v-1] + (D+1);
  }

  ESL_ALLOC(gamma_begl_s, sizeof(float **)*M);
  for (v=0; v<M; v++) {
    if (cm->stid[v] == BEGL_S) {
      /* For Bifurcatoins, we may need up to D+1 */
      ESL_ALLOC(gamma_begl_s[v], sizeof(float *) * (D+1));
      for (j=0; j<D+1; j++) 
	ESL_ALLOC(gamma_begl_s[v][j], sizeof(float)*(D+1));
    } else {
      gamma_begl_s[v] = NULL;
    }
  }
  ESL_ALLOC(gamma_begr_s, sizeof(float **)*M);
  for (v=0; v<M; v++) {
    if (cm->stid[v] == BEGR_S) {
      ESL_ALLOC(gamma_begr_s[v], sizeof(float *)*2);
      for (j=0; j<2; j++)
	ESL_ALLOC(gamma_begr_s[v][j], sizeof(float)*(D+1));
    } else {
      gamma_begr_s[v] = NULL;
    }
  }

  ESL_ALLOC(bestr, sizeof(int)*(D+1));

  /* Main recursion */
  for (j=0; j<=L; j++) {
    jmod2 = j % 2;
    if (j == 0)	
      jmin1mod2 = 1;
    else 
      jmin1mod2 = (j-1)%2;
    gamma_jmod2 = gamma[jmod2];
    gamma_jmin1mod2 = gamma[jmin1mod2];
    if (j < D) {
      minDj = j;
    } else {
      minDj = D;
    }
    for (v = M-1; v > 0; v--) {          /* Handle ROOT specially */ 
      endsc = cm->endsc[v];              /* It shouldn't change in this loop */
      sttype = cm->sttype[v];            /* This also shouldn't change */
      if (sttype == E_st) {
	gammap = gamma_jmod2[v];
	*gammap = 0.;
	for (d=1; d<=minDj; d++) {
	  gammap++;
	  *gammap = IMPOSSIBLE;  /* gamma[jmod2][v][d] */
	}
      } 
      else if (sttype == D_st || sttype == S_st) {
	y = cm->cfirst[v];
	cnum = cm->cnum[v];
	tsc = cm->tsc[v];
	for (d = 0; d <= minDj; d++)
	  sc_v[d] = endsc;
	for (yoffset= 0; yoffset < cnum; yoffset++) {
	  gammap = gamma_jmod2[y+yoffset];
#if 1
	  /*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
	  for (d = 0; d <= minDj; d++) {
	    sc = gammap[d] + tsc[yoffset];
	    if (sc > sc_v[d])
	      sc_v[d] = sc;
	  }
	}
	gammap = gamma_jmod2[v];
	for (d = 0; d<= minDj; d++) {
	  sc = sc_v[d];
	  if (sc<IMPOSSIBLE) sc = IMPOSSIBLE;
	  gammap[d] = sc;
	  if (cm->stid[v] == BEGL_S) gamma_begl_s[v][(j-d+1)%(D+1)][d] = sc;
	  if (cm->stid[v] == BEGR_S) gamma_begr_s[v][jmod2][d] = sc;
	}
      }
      else if (sttype == B_st) {
	y = cm->cfirst[v];
	z = cm->cnum[v];
	for (d = 0; d <= minDj; d++) {
	  sc = endsc;
	  gamma_begl_s_p = gamma_begl_s[y][(j-d+1)%(D+1)];
	  gamma_begr_s_p = gamma_begr_s[z][jmod2];
	  for (k=0; k<=d; k++)
	    sc_v[k] = gamma_begl_s_p[d-k];
#if 1
	  /*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
	  for (k = 0; k <= d; k++) {
	    sc_v[k] += gamma_begr_s_p[k];
	    if (sc_v[k] > sc) sc = sc_v[k];
	  }
	  if (sc<IMPOSSIBLE) sc = IMPOSSIBLE;
	  gamma_jmod2[v][d] = sc;
	}
      }
      else if (sttype == MP_st) {
	gamma_jmod2[v][0] = IMPOSSIBLE;
	if (j>0) gamma_jmod2[v][1] = IMPOSSIBLE;
	y = cm->cfirst[v];
	cnum = cm->cnum[v];
	tsc = cm->tsc[v];
	for (d = 2; d <= minDj; d++) 
	  sc_v[d] = endsc;
	for (yoffset = 0; yoffset < cnum; yoffset++) {
	  gammap = gamma_jmin1mod2[y+yoffset];
#if 1
	  /*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
	  for (d = 2; d <= minDj; d++) {
	    sc = gammap[d-2] + tsc[yoffset];
	    if (sc > sc_v[d]) {
	      sc_v[d] = sc;
	    }
	  }
	}
	for (d = 2; d <= minDj; d++) {
	  i = j-d+1;
	  sc = sc_v[d];
	  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
	    sc += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
	  else
	    sc += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
	  if (sc < IMPOSSIBLE) sc = IMPOSSIBLE;
	  gamma_jmod2[v][d] = sc;
	}
      }
      else if (sttype == ML_st) {                /* IL_st done below
						    because it points to
						    itself so gamma[j][v][d]
						    depends on gamma[j][v][d-1]
						 */
	gamma_jmod2[v][0] = IMPOSSIBLE;
	y = cm->cfirst[v];
	cnum = cm->cnum[v];
	tsc = cm->tsc[v];
	for (d = 1; d <= minDj; d++) 
	  sc_v[d] = endsc;
	for (yoffset=0; yoffset<cnum; yoffset++) {
	  gammap = gamma_jmod2[y+yoffset];
#if 1
	  /*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
	  for (d = 1; d <= minDj; d++) {
	    sc = gammap[d-1] + tsc[yoffset];
	    if (sc > sc_v[d]) {
	      sc_v[d] = sc;
	    }
	  }
	}
	for (d = 1; d <= minDj; d++) {
	  i = j-d+1;
	  sc = sc_v[d];
	  if (dsq[i] < cm->abc->K)
	    sc += cm->esc[v][(int) dsq[i]];
	  else
	    sc += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	  if (sc<IMPOSSIBLE) sc = IMPOSSIBLE;
	  gamma_jmod2[v][d] = sc;
	}
      }
      else if (sttype == IL_st) {         /* ML dealt with above, iterating
					     yoffset before d.  Can't do that
					     here because gamma[j][v][d] 
					     depends on gamma[j][v][d-1] */
	gamma_jmod2[v][0] = IMPOSSIBLE;
	y = cm->cfirst[v];
	cnum = cm->cnum[v];
	tsc = cm->tsc[v];
	for (d = 1; d <= minDj; d++) {
	  sc = endsc;
	  for (yoffset=0; yoffset<cnum; yoffset++) {
	    sc_v[yoffset] = gamma_jmod2[y+yoffset][d-1] + tsc[yoffset];
	    if (sc_v[yoffset] > sc) {
	      sc = sc_v[yoffset];
	    }
	  }
	  i = j-d+1;
	  if (dsq[i] < cm->abc->K)
	    sc += cm->esc[v][(int) dsq[i]];
	  else
	    sc += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	  if (sc<IMPOSSIBLE) sc = IMPOSSIBLE;
	  gamma_jmod2[v][d] = sc;
	}
      }
      else if (sttype == IR_st || sttype == MR_st) {
	gamma_jmod2[v][0] = IMPOSSIBLE;
	y = cm->cfirst[v];
	cnum = cm->cnum[v];
	tsc = cm->tsc[v];
	for (d = 1; d <= minDj; d++) 
	  sc_v[d] = endsc;
	for (yoffset = 0; yoffset < cnum; yoffset++) {
	  gammap = gamma_jmin1mod2[y+yoffset];
#if 1
	  /*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
	  for (d = 1; d <= minDj; d++) {
	    sc = gammap[d-1] + tsc[yoffset];
	    if (sc > sc_v[d]) {
	      sc_v[d] = sc;
	    }
	  }
	}
	for (d = 1; d <= minDj; d++) {
	  sc = sc_v[d];
	  if (dsq[j] < cm->abc->K)
	    sc += cm->esc[v][(int) dsq[j]];
	  else
	    sc += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	  if (sc < IMPOSSIBLE) sc = IMPOSSIBLE;
	  gamma_jmod2[v][d] = sc;
	}
      }
    }

    /* Now do ROOT_S (v=0) -- local begins */
    /* First do standard states to transition to */
    y = cm->cfirst[0];
    cnum=cm->cnum[0];
    tsc = cm->tsc[0];
    for (d = 0; d <= minDj; d++)
      sc_v[d] = IMPOSSIBLE;
    for (yoffset = 0; yoffset < cnum; yoffset++) {
      gammap = gamma_jmod2[y+yoffset];
#if 1
      /*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
      for (d = 0; d <= minDj; d++) {
	sc = gammap[d] + tsc[yoffset];
	if (sc > sc_v[d]) {
	  sc_v[d] = sc;
	  bestr[d] = y+yoffset;
	}
      }
    }
    /* Now, if doing local BEGINS, try that */
    if (cm->flags & CMH_LOCAL_BEGIN) {
      tsc = cm->beginsc;         /* Really cm->beginsc, not tsc */
      for (y = 1; y < M; y++) {
	beginsc = tsc[y];
	gammap = gamma_jmod2[y];
#if 1
	/*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
	for (d = 0; d <= minDj; d++) {
	  sc = gammap[d] + beginsc;
	  if (sc > sc_v[d]) {
	    sc_v[d] = sc;
	    bestr[d] = y;
	  }
	}
      }
    }
    gammap = gamma_jmod2[0];
#if 1
    /*#ifdef INTEL_COMPILER*/
#pragma ivdep
#endif
    for (d = 0; d <= minDj; d++) {
      sc = sc_v[d];
      if (sc<IMPOSSIBLE) sc = IMPOSSIBLE;
      gammap[d] = sc;
    }
  
    if (results != NULL) {
      /* Now, report the hit.  At least one hit is sent back for each j here.
	 However, some hits can already be removed for the greedy overlap
	 resolution algorithm.  Specifically, at the given j, any hit with a
	 d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
      /* First, report hit with d of 1 if > cutoff */
      if (j > 0 && gamma_jmod2[0][1] >= cutoff) 
	report_hit (j, j, bestr[1], gamma_jmod2[0][1], results);

      dmax = 1;
      /* Now, if current score is greater than maximum seen previous, report
	 it if >= cutoff and set new max */
      for (d=2; d<=minDj; d++) {
	if (gamma_jmod2[0][d] > gamma_jmod2[0][dmax]) {
	  if (j > 0 && gamma_jmod2[0][d] >= cutoff)
	    report_hit (j-d+1, j, bestr[d], gamma_jmod2[0][d], results);
	  dmax = d;
	}
      }
    }
    for (d=1; d<=minDj; d++) {
      if (j > 0 && gamma_jmod2[0][d] > best_score) {
	best_score = gamma_jmod2[0][d];
      }
    }

  }
  free(gamma[0][0]);
  free(gamma[0]);
  free(gamma);

  for (v=0; v<M; v++) {
    if (gamma_begl_s[v] != NULL) {
      for (d=0; d<=D; d++) {
	free(gamma_begl_s[v][d]);
      }
      free(gamma_begl_s[v]);
    }
  }
  free (gamma_begl_s);

  for (v=0; v<M; v++) {
    if (gamma_begr_s[v] != NULL) {
      free(gamma_begr_s[v][0]);
      free(gamma_begr_s[v][1]);
      free(gamma_begr_s[v]);
    }
  }
  free(gamma_begr_s);

  free(bestr);

  free(sc_v);

  return (best_score);

 ERROR:
  cm_Fail("Memory allocation error.");
  return IMPOSSIBLE; /* never reached */
}

/* Function: cm_CountSearchDPCalcs()
 * Date:     EPN, Tue Oct 30 14:48:00 2007
 *
 * Purpose:  Determine the number of millions of DP calcs needed to scan a seq of length <L>
 *           for the subtree rooted at each state, either using bands <dmin> and <dmax>, 
 *           or not (if <dmin> == <dmax> == NULL). 
 *           <ret_vcalcs[0]> = number of dp calcs for entire model.
 *
 * Args:     cm        - the covariance model
 *           L         - length of the sequence to search 
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           W         - max d: max size of a hit
 *           ret_vcalcs- RETURN: [0..v..M-1] number of DP calcs for scanning with sub-CM at v
 *
 * Returns:  number of calcs to search L residues with full model (ret_vcalcs[0]).
 *           dies immediately if some error occurs.
 */
float 
cm_CountSearchDPCalcs(CM_t *cm, int L, int *dmin, int *dmax, int W, float **ret_vcalcs)
{
  int       status;
  float     *vcalcs;            /* [0..v..cm->M-1] # of calcs for subtree rooted at v */
  int       d;			/* a subsequence length, 0..W */
  int       j;                  /* seq index */
  int       v, w, y;            /* state indices */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       dn;                 /* temporary value for min d in for loops */
  int       dx;                 /* temporary value for max d in for loops */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */

  if(dmin != NULL && dmax != NULL) do_banded = TRUE;
  if (W > L) W = L; 

  ESL_ALLOC(vcalcs, sizeof(float) * cm->M);
  esl_vec_FSet(vcalcs, cm->M, 0.);

  /* we ignore initialization and band imposition, a little imprecise */
  /* Recursion. */
  for (j = 1; j <= L; j++) {
    for (v = cm->M-1; v > 0; v--) { /* ...almost to ROOT; we handle ROOT specially... */
      if(do_banded) { 
	dn = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	dx = ESL_MIN(j, dmax[v]); 
	dx = ESL_MIN(dx, W);
      }
      else { 
	dn = (cm->sttype[v] == MP_st) ? 2 : 1;
	dx = ESL_MIN(j, W); 
      }
      if(cm->sttype[v] == E_st) continue;
      
      if(cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S */
	y = cm->cnum[v];   /* BEGR_S */
	for (d = dn; d <= dx; d++) {
	  if(do_banded) {
	    kmin = ESL_MAX(dmin[y], (d-dmax[w]));
	    kmin = ESL_MAX(kmin, 0);
	    kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	  }
	  else { kmin = 0; kmax = d; }
	  vcalcs[v] += 1 + (kmax - kmin + 1); /* initial '1 +' is for initialization calc */
	} /* ! B_st */
      }
      else  { /* if cm->sttype[v] != B_st */
	vcalcs[v] += 1 + (cm->cnum[v] + 1) * (dx - dn + 1); /* 1 is for initialization calc */
	if(StateDelta(cm->sttype[v]) > 0) vcalcs[v] += (dx-dn+1);
      } /* end of else (v != B_st) */
    } /*loop over decks v>0 */
    
    /* determine min/max d we're allowing for the root state and this position j */
    if(do_banded) { 
      dn = ESL_MAX(dmin[0], 1); 
      dx = ESL_MIN(j, dmax[0]); 
      dx = ESL_MIN(dx, W);
    }
    else { 
      dn = 1; 
      dx = ESL_MIN(j, W); 
    }
    vcalcs[0] += (cm->cnum[v] + 1) * (dx - dn + 1);
    
    if (cm->flags & CMH_LOCAL_BEGIN) {
      for (y = 1; y < cm->M; y++) {
	if(do_banded) {
	  dn = (cm->sttype[y] == MP_st) ? ESL_MAX(dmin[y], 2) : ESL_MAX(dmin[y], 1); 
	  dn = ESL_MAX(dn, dmin[0]);
	  dx = ESL_MIN(j, dmax[y]); 
	  dx = ESL_MIN(dx, W);
	}
	else { 
	  dn = 1; 
	  dx = ESL_MIN(j, W); 
	}
	vcalcs[0] += (dx - dn + 1);
      }
    }
  } /* end loop over end positions j */
  
  /* sum up the cells for all states under each v */
  for (v = cm->M-1; v >= 0; v--) {
    if     (cm->sttype[v] == B_st) vcalcs[v] += vcalcs[cm->cnum[v]] + vcalcs[cm->cfirst[v]];
    else if(cm->sttype[v] != E_st) vcalcs[v] += vcalcs[v+1];
  }
  
  /* convert to millions of calcs */
  for (v = cm->M-1; v >= 0; v--) vcalcs[v] /= 1000000.;
  /* convert to per residue */
  for (v = cm->M-1; v >= 0; v--) vcalcs[v] /= L;

  /* for (v = cm->M-1; v >= 0; v--) printf("vcalcs[%4d]: %.3f\n", v, vcalcs[v]); */
  
  *ret_vcalcs = vcalcs;
  return vcalcs[0];

 ERROR:
  cm_Fail("count_search_target_calcs(): memory error.");
  return 0.;
}

/* EPN 03.29.06
 * Function: FastCYKScanHB()
 * Incept:   EPN, Mon Nov 12 17:45:57 2007
 *
 * Purpose:  An HMM banded version of a scanning CYK algorithm. Takes
 *           a CM_FHB_MX data structure which is indexed [v][j][d] with
 *           only cells within the bands allocated.
 *           (different than other scanning functions convention of [j][v][d]).
 *
 * Args:     cm        - the model    [0..M-1]
 *           sq        - the sequence [1..L]   
 *                     - length of the dsq
 *           vroot     - first start state of subtree (0, for whole model)
 *           vend      - last end state of subtree (cm->M-1, for whole model)
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           mx        - the dp matrix, only cells within bands in cp9b will 
 *                       be valid. 
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           cp9b      - HMM bands for <dsq> and <cm>
 *                       
 * Returns: Score of the best hit.
 */
float 
FastCYKScanHB(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
	      CM_FHB_MX *mx, int allow_begin, int *ret_b, float *ret_bsc, CP9Bands_t *cp9b)
{
  int      status;
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */
  int     *yvalidA;     /* [0..MAXCONNECT-1] TRUE if v->yoffset is legal transition (within bands) */
  float   *el_scA;      /* [0..d..W-1] probability of local end emissions of length d */
  /* variables used for memory efficient bands */
  /* ptrs to cp9b info, for convenience */
  int     *jmin  = cp9b->jmin;  
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* the DP matrix */
  float ***alpha = mx->dp; /* pointer to the alpha DP matrix */
  /* indices used for handling band-offset issues, and in the depths of the DP recursion */
  int      sd;                 /* StateDelta(cm->sttype[v]) */
  int      sdr;                /* StateRightDelta(cm->sttype[v] */
  int      jp_v, jp_y, jp_z;   /* offset j index for states v, y, z */
  int      jp_y_sdr;           /* jp_y - sdr */
  int      j_sdr;              /* j - sdr */
  int      jn, jx;             /* current minimum/maximum j allowed */
  int      jpn, jpx;           /* minimum/maximum jp_v */
  int      dp_v, dp_y;         /* d index for state v/y in alpha w/mem eff bands */
  int      dn, dx;             /* current minimum/maximum d allowed */
  int      dp_y_sd;            /* dp_y - sd */
  int      dpn, dpx;           /* minimum/maximum dp_v */
  int      kp_z;               /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      kn, kx;             /* current minimum/maximum k value */
  int      Wp;                 /* W oalso changes depending on state */
  float    tsc;                /* a transition score */
  int      yvalid_idx;         /* for keeping track of which children are valid */
  int      yvalid_ct;          /* for keeping track of which children are valid */

  /* Contract check */
  if(dsq == NULL) cm_Fail("fast_cyk_inside_align_hb(), dsq is NULL.\n");
  if (mx == NULL) cm_Fail("fast_cyk_inside_align_hb(), mx is NULL.\n");

  /* Allocations and initializations  */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  /* grow the matrix based on the current sequence and bands */
  cm_fhb_mx_GrowTo(mx, cp9b);

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* yvalidA[0..cnum[v]] will hold TRUE for states y for which a transition is legal 
   * (some transitions are impossible due to the bands) */
  ESL_ALLOC(yvalidA, sizeof(int) * MAXCONNECT);
  esl_vec_ISet(yvalidA, MAXCONNECT, FALSE);

  /* initialize all cells of the matrix to IMPOSSIBLE */
  //esl_vec_FSet(alpha[0][0], mx->ncells_valid, IMPOSSIBLE);

  /* Main recursion */
  for (v = vend; v >= vroot; v--) {
    float const *esc_v = cm->oesc[v]; /* emission scores for state v */
    float const *tsc_v = cm->tsc[v];  /* transition scores for state v */
    sd   = StateDelta(cm->sttype[v]);
    sdr  = StateRightDelta(cm->sttype[v]);
    jn   = jmin[v];
    jx   = jmax[v];
    /* Get a shadow deck to fill in and initialize all valid cells for state v */
    if (cm->sttype[v] != E_st) {
      if (cm->sttype[v] == B_st) {
	/* initialize all valid cells for state v to IMPOSSIBLE (local ends are impossible for B states) */
	assert(! (NOT_IMPOSSIBLE(cm->endsc[v])));
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	    alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	  }
	}
      } else { /* ! B_st && ! E_st */
	/* initialize all valid cells for state v */
	if(NOT_IMPOSSIBLE(cm->endsc[v])) {
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++) {
	      alpha[v][jp_v][dp_v] = el_scA[d-sd] + cm->endsc[v];
	    }
	  }
	}
	else { /* cm->endsc[v] == IMPOSSIBLE */
	  for (j = jmin[v]; j <= jmax[v]; j++) { 
	    jp_v  = j - jmin[v];
	    for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++) {
	      alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	    }
	  }
	}
      }
    }

    if(cm->sttype[v] == E_st) { 
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v = j-jmin[v];
	ESL_DASSERT1((hdmin[v][jp_v] == 0));
	ESL_DASSERT1((hdmax[v][jp_v] == 0));
	alpha[v][jp_v][0] = 0.; /* for End states, d must be 0 */
      }
    }
    else if(cm->sttype[v] == IL_st) {
      /* update alpha[v][jp_v][dp_v] cells, for IL states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
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
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > alpha[v][jp_v][dp_v])
		alpha[v][jp_v][dp_v] = sc; 
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] == IR_st) { 
      /* update alpha[v][jp_v][dp_v] cells, for IR states, loop nesting order is:
       * for j { for d { for y { } } } because they can self transit, and a 
       * alpha[v][j][d] cell must be complete (that is we must have looked at all children y) 
       * before can start calc'ing for alpha[v][j][d+1] */
      for (j = jmin[v]; j <= jmax[v]; j++) {
	jp_v = j - jmin[v];
	yvalid_ct = 0;
	j_sdr = j - sdr;
	
	/* determine which children y we can legally transit to for v, j */
	for (y = cm->cfirst[v], yoffset = 0; y < (cm->cfirst[v] + cm->cnum[v]); y++, yoffset++) 
	  if((j_sdr) >= jmin[y] && ((j_sdr) <= jmax[y])) yvalidA[yvalid_ct++] = yoffset; /* is j-sdr is valid for state y? */
	
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) { /* for each valid d for v, j */
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha */
	  for (yvalid_idx = 0; yvalid_idx < yvalid_ct; yvalid_idx++) { /* for each valid child y, for v, j */
	    yoffset = yvalidA[yvalid_idx];
	    y = cm->cfirst[v] + yoffset;
	    jp_y_sdr = j - jmin[y] - sdr;
	    
	    if((d-sd) >= hdmin[y][jp_y_sdr] && (d-sd) <= hdmax[y][jp_y_sdr]) { /* make sure d is valid for this v, j and y */
	      dp_y_sd = d - sd - hdmin[y][jp_y_sdr];
	      ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	      ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	      if ((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc_v[yoffset]) > alpha[v][jp_v][dp_v])
		{
		  alpha[v][jp_v][dp_v] = sc; 
		}
	    }
	  }
	  alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
	}
      }
    }
    else if(cm->sttype[v] != B_st) { /* entered if state v is (! IL && ! IR && ! B) */
      /* ML, MP, MR, D, S, E states cannot self transit, this means that all cells
       * in beta[v] are independent of each other, only depending on beta[y] for previously calc'ed y.
       * We can do the for loops in any nesting order, this implementation does what I think is most efficient:
       * for y { for j { for d { } } } 
       */
      for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) {
	yoffset = y - cm->cfirst[v];
	tsc = tsc_v[yoffset];
	
	jn = ESL_MAX(jmin[v], ESL_MIN(jmin[y] + sdr, jmax[y]));
	jx = ESL_MIN(jmax[v], ESL_MAX(jmax[y] + sdr, jmin[y]));
	jx = ESL_MIN(jx, jmax[y] + sdr);
	jpn = jn - jmin[v];
	jpx = jx - jmin[v];
	jp_y_sdr = jn - jmin[y] - sdr;
	
	for (jp_v = jpn; jp_v <= jpx; jp_v++, jp_y_sdr++) {
	  ESL_DASSERT1((jp_v >= 0 && jp_v <= (jmax[v]-jmin[v])));
	  ESL_DASSERT1((jp_y_sdr >= 0 && jp_y_sdr <= (jmax[y]-jmin[y])));
	  
	  dn = ESL_MAX(hdmin[v][jp_v], hdmin[y][jp_y_sdr] + sd);
	  dx = ESL_MIN(hdmax[v][jp_v], hdmax[y][jp_y_sdr] + sd);
	  dpn     = dn - hdmin[v][jp_v];
	  dpx     = dx - hdmin[v][jp_v];
	  dp_y_sd = dn - hdmin[y][jp_y_sdr] - sd;
	  	  
	  for (dp_v = dpn; dp_v <= dpx; dp_v++, dp_y_sd++) { 
	    ESL_DASSERT1((dp_v    >= 0 && dp_v     <= (hdmax[v][jp_v]     - hdmin[v][jp_v])));
	    ESL_DASSERT1((dp_y_sd >= 0 && dp_y_sd  <= (hdmax[y][jp_y_sdr] - hdmin[y][jp_y_sdr])));
	    if((sc = alpha[y][jp_y_sdr][dp_y_sd] + tsc) > alpha[v][jp_v][dp_v]) {
	      alpha[v][jp_v][dp_v] = sc;
	    }
	  }
	}
      }
      /* add in emission score, if any */
      switch(cm->sttype[v]) { 
      case ML_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]];
	}
	break;
      case MR_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[j]];
	}
	break;
      case MP_st:
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	    alpha[v][jp_v][dp_v] += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
	}
      default:
	break;
      }
      /* ensure all cells are >= IMPOSSIBLE */
      for (j = jmin[v]; j <= jmax[v]; j++) { 
	jp_v  = j - jmin[v];
	for (dp_v = 0; dp_v <= (hdmax[v][jp_v] - hdmin[v][jp_v]); dp_v++)
	  alpha[v][jp_v][dp_v] = ESL_MAX(alpha[v][jp_v][dp_v], IMPOSSIBLE);
      }
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
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++) {
	  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	  /* Find the first k value that implies a valid cell in the y and z decks.
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

	      if ((sc = alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]) 
		  > alpha[v][jp_v][dp_v]) { 
		alpha[v][jp_v][dp_v] = sc;
	      }
	    }
	  }
	}
      }
    }				/* finished calculating deck v. */
         
    /* The following loops originally access alpha[v][j0][W] but the index W will be
       in different positions due to the bands */
    
    if(allow_begin) { 
      if(j0 >= jmin[v] && j0 <= jmax[v]) { 
	jp_v = j0 - jmin[v];
	Wp = W - hdmin[v][jp_v];
	if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	/* If we get here alpha[v][jp_v][Wp] is a valid cell
	 * in the banded alpha matrix, corresponding to 
	 * alpha[v][j0][W] in the platonic matrix.
	 * Check for local begin getting us to the root. */
	  if (alpha[v][jp_v][Wp] + cm->beginsc[v] > bsc) { 
	    b   = v;
	    bsc = alpha[v][jp_v][Wp] + cm->beginsc[v];
	  }
	}
      }
      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (v == 0) { 
	if(j0 >= jmin[0] && j0 <= jmax[0]) {
	  jp_v = j0 - jmin[v];
	  Wp   = W - hdmin[v][jp_v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v]) { 
	    if (bsc > alpha[0][jp_v][Wp]) {
	      alpha[0][jp_v][Wp] = bsc;
	    }
	  }
	}
      }
    }
  } /* end loop over all v */
  /*FILE *fp; fp = fopen("cyk.mx", "w"); cm_fhb_mx_Dump(fp, mx); fclose(fp);*/
  
  Wp = W - hdmin[vroot][j0-jmin[vroot]];
  sc =     alpha[vroot][j0-jmin[vroot]][Wp];

  if (ret_b != NULL)      *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL)    *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */
  free(el_scA);
  free(yvalidA);

  ESL_DPRINTF1(("FastCYKScanHB() return sc: %f\n", sc));
  return sc;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}



/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_FASTSEARCH_BENCHMARK
/* gcc -o benchmark-fastsearch -g -O2 -I. -L. -I../easel -L../easel -DIMPL_FASTSEARCH_BENCHMARK cm_fastsearch.c -linfernal -leasel -lm
 * icc -o benchmark-fastsearch -O3 -static -I. -L. -I../easel -L../easel -DIMPL_FASTSEARCH_BENCHMARK cm_fastsearch.c -linfernal -leasel -lm
 * ./benchmark-fastsearch <cmfile>
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
  { "-e",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "emit sequences from CM, don't randomly create them", 0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search in glocal mode [default: local]", 0 },
  { "-L",        eslARG_INT,  "10000", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-o",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute old reference CYK scan implementation", 0 },
  { "-w",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute new reference CYK scan implementation", 0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute experimental CYK scan implementation", 0 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute non-banded optimized CYK scan implementation", 0 },
  { "--rsearch", eslARG_NONE,   FALSE, NULL, NULL,  NULL,"--noqdb", NULL, "also execute ported RSEARCH's CYK scan implementation", 0 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute optimized int inside scan implementation", 0 },
  { "--riins",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute reference int inside scan implementation", 0 },
  { "--oiins",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute old int inside scan implementation", 0 },
  { "--fins",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute optimized float inside scan implementation", 0 },
  { "--rfins",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute reference float inside scan implementation", 0 },
  { "--ofins",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,  "also execute old float inside scan implementation", 0 },
  { "--hbanded", eslARG_NONE,   FALSE, NULL, NULL,  NULL,"-e",   NULL, "calculate and use HMM bands in CM search", 6 },
  { "--tau",     eslARG_REAL,   "1e-7",NULL, "0<x<1",NULL,"--hbanded",  NULL, "set tail loss prob for --hbanded to <x>", 6 },
  { "--scan2bands",eslARG_NONE, FALSE, NULL, NULL,  NULL,"--hbanded",   NULL, "derive HMM bands from scanning Forward/Backward", 6 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,"--hbanded",   NULL, "use posterior sums during HMM band calculation (widens bands)", 6 },
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

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  do_random = TRUE;
  if(esl_opt_GetBoolean(go, "-e")) do_random = FALSE; 

  if(! esl_opt_GetBoolean(go, "-g")) cm->config_opts  |= CM_CONFIG_LOCAL;
  if(  esl_opt_GetBoolean(go, "--sums"))        cm->search_opts |= CM_SEARCH_SUMS;
  if(  esl_opt_GetBoolean(go, "--scan2bands"))  cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts |= CM_SEARCH_HBANDED;
  cm->tau    = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */
  cm->config_opts |= CM_CONFIG_QDB;
  ConfigCM(cm, NULL, NULL);

  if (esl_opt_GetBoolean(go, "--noqdb")) { 
    dmin = NULL; dmax = NULL;
  }
  else { dmin = cm->dmin; dmax = cm->dmax; }

  cm_CreateScanInfo(cm, TRUE, TRUE); /* impt to do this after QDBs set up in ConfigCM() */

  CM_FHB_MX  *hbmx = NULL;
  CP9Bands_t *cp9b = NULL;
  if (esl_opt_GetBoolean(go, "--hbanded")) { 
    cp9b = AllocCP9Bands(cm, cm->cp9);
    /* create the matrix, it'll be empty initially */
    hbmx = cm_fhb_mx_Create(cm->M);
  }

  /* get sequences */
  if(do_random) {
    double *dnull;
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);
    seqs_to_aln = RandomEmitSeqsToAln(r, cm->abc, dnull, 1, N, L, FALSE);
    free(dnull);
  }
  else /* don't randomly generate seqs, emit them from the CM */
    seqs_to_aln = CMEmitSeqsToAln(r, cm, 1, N, FALSE);

  for (i = 0; i < N; i++)
    {
      L = seqs_to_aln->sq[i]->n;
      dsq = seqs_to_aln->sq[i]->dsq;
      cm->search_opts  &= ~CM_SEARCH_INSIDE;

      esl_stopwatch_Start(w);
      sc = FastCYKScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      /*if (esl_opt_GetBoolean(go, "-x")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = XFastCYKScan(cm, dsq, dmin, dmax, 1, L, cm->W, 0., NULL, NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "XFastCYKScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	  }*/

      if (esl_opt_GetBoolean(go, "-w")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = RefCYKScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "RefCYKScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "-o")) 
	{ 
	  esl_stopwatch_Start(w);
	  if (esl_opt_GetBoolean(go, "--noqdb")) sc = CYKScan (cm, dsq, 1, L, cm->W, 0., NULL); 
	  else                                   sc = CYKBandedScan (cm, dsq, dmin, dmax, 1, L, cm->W, 0., NULL); 
	  printf("%4d %-30s %10.4f bits ", (i+1), "CYKBandedScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--rsearch")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = rsearch_CYKScan (cm, dsq, L, 0., cm->W, NULL); 
	  printf("%4d %-30s %10.4f bits ", (i+1), "rsearch_CYKScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      /* integer inside implementations */
      if (esl_opt_GetBoolean(go, "--iins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  sc = FastIInsideScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "FastIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");

	  esl_stopwatch_Start(w);
	  sc = XFastIInsideScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "XFastIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");

	  esl_stopwatch_Start(w);
	  sc = X2FastIInsideScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "X2FastIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--riins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  sc = RefIInsideScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "RefIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");

	  esl_stopwatch_Start(w);
	  sc = XRefIInsideScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "XRefIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--oiins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  if(esl_opt_GetBoolean(go, "--noqdb")) { 
	    esl_stopwatch_Start(w);
	    sc = iInsideScan(cm, dsq, 1, L, cm->W, 0., NULL);
	    printf("%4d %-30s %10.4f bits ", (i+1), "iInsideScan() (int no-qdb): ", sc);
	    esl_stopwatch_Stop(w);
	    esl_stopwatch_Display(stdout, w, " CPU time: ");
	  }
	  else { 
	    esl_stopwatch_Start(w);
	    sc = iInsideBandedScan(cm, dsq, dmin, dmax, 1, L, cm->W, 0., NULL);
	    printf("%4d %-30s %10.4f bits ", (i+1), "iInsideBandedScan() (int): ", sc);
	    esl_stopwatch_Stop(w);
	    esl_stopwatch_Display(stdout, w, " CPU time: ");
	  }
	}

      /* float inside implementations */
      if (esl_opt_GetBoolean(go, "--fins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  sc = FastFInsideScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "FastFInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--rfins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  sc = RefFInsideScan(cm, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "RefFInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--ofins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  if(esl_opt_GetBoolean(go, "--noqdb")) { 
	    esl_stopwatch_Start(w);
	    sc = InsideScan(cm, dsq, 1, L, cm->W, 0., NULL);
	    printf("%4d %-30s %10.4f bits ", (i+1), "InsideScan() (float no-qdb): ", sc);
	    esl_stopwatch_Stop(w);
	    esl_stopwatch_Display(stdout, w, " CPU time: ");
	  }
	  else { 
	    esl_stopwatch_Start(w);
	    sc = InsideBandedScan(cm, dsq, dmin, dmax, 1, L, cm->W, 0., NULL);
	    printf("%4d %-30s %10.4f bits ", (i+1), "InsideBandedScan() (float): ", sc);
	    esl_stopwatch_Stop(w);
	    esl_stopwatch_Display(stdout, w, " CPU time: ");
	  }
	}
      if (esl_opt_GetBoolean(go, "--hbanded")) 
	{ 
	  esl_stopwatch_Start(w);
	  cp9_Seq2Bands(cm, dsq, 1, L, cp9b, 0);   /* debug level */
	  sc = CYKBandedScan_jd(cm, dsq, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
				1, L, cm->W, 0., NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "CYKBandedScan_jd() : ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");

	  esl_stopwatch_Start(w);
	  cp9_Seq2Bands(cm, dsq, 1, L, cp9b, 0);   /* debug level */
	  sc = FastCYKScanHB(cm, dsq, L, 0, cm->M-1, 1, L, hbmx, (cm->flags & CMH_LOCAL_BEGIN), NULL, NULL, cp9b);
	  printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScanHB() : ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      printf("\n");
    }
  if(cp9b != NULL) FreeCP9Bands(cp9b);
  if(hbmx != NULL) cm_fhb_mx_Destroy(hbmx);
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
#endif /*IMPL_FASTSEARCH_BENCHMARK*/


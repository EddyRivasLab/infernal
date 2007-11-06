/* cm_fastsearch.c
 * EPN, Wed Sep 12 16:53:32 2007
 * 
 * Fast versions of CYK and Inside search functions.
 * 
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
 *           bands are used or not used as specified in ScanInfo_t <si>.
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
 * Note:     This function is heavily synchronized with FastFInsideScan() and FastIInsideScan(),
 *           any change to this function should be mirrored in those functions. 
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
FastCYKScan(CM_t *cm, ScanInfo_t *si, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
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

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in FastCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in FastCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in FastCYKScan, CM_SEARCH_INSIDE flag raised");

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
  float *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* make pointers to the ScanInfo data for convenience */
  float ***alpha      = si->alpha;
  float ***alpha_begl = si->alpha_begl;
  int   **dnAA        = si->dnAA;
  int   **dxAA        = si->dxAA;
  int    *emitmodeA   = si->emitmodeA;
  float  **esc_vAA    = si->esc_vAA;
  float **init_scAA   = si->init_scAA;
  int    *bestr       = si->bestr;
  int    *dmin        = si->dmin;
  int    *dmax        = si->dmax;

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
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	      }
	      else { kmin = 0; kmax = d; }

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = ESL_MAX(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
		sc = ESL_MAX(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc = ESL_MAX(sc, init_scAA[v][dp_y]);
		sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ESL_MAX(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc = ESL_MAX(sc, arow4[dp_y] + tsc_v[4]);		
		sc = ESL_MAX(sc, arow3[dp_y] + tsc_v[3]);		
		sc = ESL_MAX(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ESL_MAX(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ESL_MAX(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc = ESL_MAX(sc, init_scAA[v][dp_y]);
		sc = ESL_MAX(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ESL_MAX(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc = ESL_MAX(sc, init_scAA[v][dp_y]);
		sc = ESL_MAX(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ESL_MAX(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		

		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ESL_MAX(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
		switch (emitmodeA[v]) {
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep 
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = ESL_MAX(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc_v[d] = ESL_MAX(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      }
	      break; 
	    } /* end of switch(cnum) */
	    /* add in emission score (if any), and set alpha[jp_v][v][d] cell */
	    switch (emitmodeA[v]) {
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
	    } /* end of switch (emitmodeA[v]) */
	  } /* end of else (cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st) */
	  if(vsc != NULL) 
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
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
	
      if (cm->flags & CM_LOCAL_BEGIN) {
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
FastIInsideScan(CM_t *cm, ScanInfo_t *si, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
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
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in FastIInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in FastIInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in FastIInsideScan, CM_SEARCH_INSIDE flag not raised");

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
  int *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_ISet(sc_v, (W+1), -INFTY);

  /* make pointers to the ScanInfo data for convenience */
  int  ***alpha      = si->ialpha;
  int  ***alpha_begl = si->ialpha_begl;
  int   **dnAA       = si->dnAA;
  int   **dxAA       = si->dxAA;
  int    *emitmodeA  = si->emitmodeA;
  int   **esc_vAA    = si->iesc_vAA;
  int   **init_scAA  = si->iinit_scAA;
  int    *dmin       = si->dmin;
  int    *dmax       = si->dmax;

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
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	      }
	      else { kmin = 0; kmax = d; }

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = ILogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
		sc = ILogsum(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ILogsum(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow4[dp_y] + tsc_v[4]);		
		sc = ILogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ILogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ILogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc = ILogsum(sc, init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = ILogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		

		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = ILogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc = ILogsum(sc, arow0[dp_y] + tsc_v[0]);		
		switch (emitmodeA[v]) {
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (int * const) alpha[jp_y][y];
	      arow1 = (int * const) alpha[jp_y][y+1];
	      arow2 = (int * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep 
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = ILogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc_v[d] = ILogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      }
	      break; 
	    } /* end of switch(cnum) */
	    /* add in emission score (if any), and set alpha[jp_v][v][d] cell */
	    switch (emitmodeA[v]) {
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
	    } /* end of switch (emitmodeA[v]) */
	  } /* end of else (cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st) */
	  if(vsc != NULL) 
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha[jp_v][v][d]));
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha_begl[jp_v][v][d]));
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
	
      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) 
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->ibeginsc[y]);
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) 
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->ibeginsc[y]);
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
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("FastIInsideScan() return score: %10.4f\n", vsc_root)); 
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
RefIInsideScan(CM_t *cm, ScanInfo_t *si, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
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
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in FastIInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in FastIInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in FastIInsideScan, CM_SEARCH_INSIDE flag not raised");

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
  int *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_ISet(sc_v, (W+1), -INFTY);

  /* make pointers to the ScanInfo data for convenience */
  int  ***alpha      = si->ialpha;
  int  ***alpha_begl = si->ialpha_begl;
  int   **dnAA       = si->dnAA;
  int   **dxAA       = si->dxAA;
  int    *emitmodeA  = si->emitmodeA;
  int   **esc_vAA    = si->iesc_vAA;
  int   **init_scAA  = si->iinit_scAA;
  int    *dmin       = si->dmin;
  int    *dmax       = si->dmax;

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

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = ILogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
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
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ILogsum(sc, alpha[jp_y][y+yoffset][d - sd] + tsc_v[yoffset]);

	      switch (emitmodeA[v]) {
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
	      } /* end of switch emitmodeA[v] */
	    } /* end of for d loop */
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) 
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha[jp_v][v][d]));
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], Scorify(alpha_begl[jp_v][v][d]));
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
	
      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(cm->ibeginsc[y] != -INFTY) {
	    if(cm->stid[y] == BEGL_S) {
	      jp_y = jp_wA[0];
	      for (d = dnA[y]; d <= dxA[y]; d++) 
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha_begl[jp_y][y][d] + cm->ibeginsc[y]);
	    }
	    else { /* y != BEGL_S */
	      jp_y = cur;
	      for (d = dnA[y]; d <= dxA[y]; d++) 
		alpha[jp_v][0][d] = ILogsum(alpha[jp_v][0][d], alpha[jp_y][y][d] + cm->ibeginsc[y]);
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
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("FastCYKScan() return score: %10.4f\n", vsc_root)); 
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
 *           si              - ScanInfo_t for this model (includes alpha DP matrix, qdbands etc.) 
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
FastFInsideScan(CM_t *cm, ScanInfo_t *si, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
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
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in FastFInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in FastFInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in FastFInsideScan, CM_SEARCH_INSIDE flag not raised");

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
  float *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* make pointers to the ScanInfo data for convenience */
  float ***alpha      = si->alpha;
  float ***alpha_begl = si->alpha_begl;
  int   **dnAA        = si->dnAA;
  int   **dxAA        = si->dxAA;
  int    *emitmodeA   = si->emitmodeA;
  float  **esc_vAA    = si->esc_vAA;
  float **init_scAA   = si->init_scAA;
  int    *dmin        = si->dmin;
  int    *dmax        = si->dmax;

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
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	      }
	      else { kmin = 0; kmax = d; }

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = FLogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
		sc = FLogsum(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow4[dp_y] + tsc_v[4]);		
		sc = FLogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = FLogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		

		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		switch (emitmodeA[v]) {
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep 
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = FLogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc_v[d] = FLogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      }
	      break; 
	    } /* end of switch(cnum) */
	    /* add in emission score (if any), and set alpha[jp_v][v][d] cell */
	    switch (emitmodeA[v]) {
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
	    } /* end of switch (emitmodeA[v]) */
	  } /* end of else (cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st) */
	  if(vsc != NULL) 
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
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
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = FLogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CM_LOCAL_BEGIN) {
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
  if (ret_vsc != NULL) *ret_vsc = vsc;
  
  ESL_DPRINTF1(("FastCYKScan() return score: %10.4f\n", vsc_root)); 
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
RefCYKScan(CM_t *cm, ScanInfo_t *si, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
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
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in RefCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in RefCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in RefCYKScan, CM_SEARCH_INSIDE flag raised");

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

  /* Initialize sc_v to size of M */
  float *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* make pointers to the ScanInfo data for convenience */
  float ***alpha      = si->alpha;
  float ***alpha_begl = si->alpha_begl;
  int   **dnAA        = si->dnAA;
  int   **dxAA        = si->dxAA;
  int    *emitmodeA   = si->emitmodeA;
  float  **esc_vAA    = si->esc_vAA;
  float **init_scAA   = si->init_scAA;
  int    *bestr       = si->bestr;
  int    *dmin        = si->dmin;
  int    *dmax        = si->dmax;

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

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = ESL_MAX(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ESL_MAX(sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* ! B_st, ! BEGL_S st */
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ESL_MAX(sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);

	      switch (emitmodeA[v]) {
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
	      } /* end of switch emitmodeA[v] */
	    } /* end of for d loop */
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) 
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
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
	
      if (cm->flags & CM_LOCAL_BEGIN) {
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
  if (ret_vsc != NULL) *ret_vsc         = vsc;
  
  ESL_DPRINTF1(("RefCYKScan() return score: %10.4f\n", vsc_root)); 
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
RefFInsideScan(CM_t *cm, ScanInfo_t *si, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, 
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
  int       ip_g;               /* offset i for gamma (i-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int      *dnA, *dxA;          /* tmp ptr to 1 row of dnAA, dxAA */
  int       dn,   dx;           /* minimum/maximum valid d for current state */
  int       cnum;               /* number of children for current state */
  int      *jp_wA;              /* rolling pointer index for B states, gets precalc'ed */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in RefFInsideScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in RefFInsideScan, dsq is NULL\n");
  if(! (cm->search_opts & CM_SEARCH_INSIDE)) cm_Fail("ERROR in RefFInsideScan, CM_SEARCH_INSIDE flag not raised");

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
  float *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* make pointers to the ScanInfo data for convenience */
  float ***alpha      = si->alpha;
  float ***alpha_begl = si->alpha_begl;
  int   **dnAA        = si->dnAA;
  int   **dxAA        = si->dxAA;
  int    *emitmodeA   = si->emitmodeA;
  float  **esc_vAA    = si->esc_vAA;
  float **init_scAA   = si->init_scAA;
  float  *el_scA      = si->el_scA;
  int    *dmin        = si->dmin;
  int    *dmax        = si->dmax;

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

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = FLogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = FLogsum(sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* ! B_st, ! BEGL_S st */
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = FLogsum(sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);

	      switch (emitmodeA[v]) {
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
	      } /* end of switch emitmodeA[v] */
	    } /* end of for d loop */
	  } /* end of else (which was entered if ! B_st && ! BEGL_S st) */
	  if(vsc != NULL) 
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
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
	
      if (cm->flags & CM_LOCAL_BEGIN) {
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
  if (ret_vsc != NULL) *ret_vsc         = vsc;
  
  ESL_DPRINTF1(("RefFInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}


/* Function: XFastCYKScan()
 * Date:     EPN, Sat Oct 20 07:48:05 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm. If bands are NULL, reverts to non-banded
 *           (scancyk.c:CYKScan()). 
 *
 * Args:     cm              - the covariance model
 *           dsq             - the digitized sequence
 *           dmin            - minimum bound on d for state v; 0..M
 *           dmax            - maximum bound on d for state v; 0..M          
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *           ret_best_hit_sc - RETURN score of best hit (reported to results) NULL if not-wanted
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
XFastCYKScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float cutoff, 
	    search_results_t *results, float **ret_vsc, float *ret_best_hit_sc)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [j][v][d] */
  float  ***alpha_begl; 
  float    *vsc;                /* best score for each state (float) */
  float     vsc_root;           /* best overall score (score at ROOT_S) */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       ip_g;               /* offset i for gamma (i-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  float     best_hit_sc;        /* best hit score found */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int     **dnAA, **dxAA;       /* [1..j..W][0..v..M-1] min,max d value allowed for posn j, state v */
  int      *dnA,   *dxA;        /* tmp ptr to 1 row of dnAA, dxAA */

  int cnum;
  int dn;
  int dx;
  float *el_scA;
  int *jp_wA;
  float **init_scAA;
  int  ctr = 0;
  /*int yidx;*/
  /*float const *tsc = cm->tsc[0]; */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in XFastCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in XFastCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in XFastCYKScan, CM_SEARCH_INSIDE flag raised");

  /* determine if we're doing banded/non-banded */
  if(dmin != NULL && dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 

  vsc = NULL;
  if(ret_vsc != NULL) { 
    ESL_ALLOC(vsc, sizeof(float) * cm->M);
    esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);
  }
  best_hit_sc = IMPOSSIBLE;
  vsc_root    = IMPOSSIBLE;

  /*
   * alpha and alpha_begl allocations.
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
   * Note if v is NOT a BEGL_S alpha[j][v] == NULL
   */

  /* allocate alpha */
  ESL_ALLOC(alpha, (sizeof(float **) * 2));
  ESL_ALLOC(alpha[0], sizeof(float *) * cm->M);
  ESL_ALLOC(alpha[1], sizeof(float *) * cm->M);
  ESL_ALLOC(alpha[0][0], (sizeof(float) * 2 * (cm->M) * (W+1)));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      alpha[0][v] = alpha[0][0] + (v           * (W+1));
      alpha[1][v] = alpha[0][0] + ((v + cm->M) * (W+1));
    }
    else { /* BEGL_S, this is wasteful */
      alpha[0][v] = NULL;
      alpha[1][v] = NULL;
    }
  }
  /* float const *alphap = alpha[0][0]; */

  /* allocate alpha_begl */
  ESL_ALLOC(alpha_begl, (sizeof(float **) * (W+1)));
  for (j = 0; j <= W; j++) {
    ESL_ALLOC(alpha_begl[j], (sizeof(float *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	ESL_ALLOC(alpha_begl[j][v], (sizeof(float) * (W+1)));
      }
      else /* non-BEGL */
	alpha_begl[j][v] = NULL;
    }
  }
  ESL_ALLOC(bestr, (sizeof(int) * (W+1)));

  /*
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   */
  /* initialize alpha and alpha_begl */
  for(v = cm->M-1; v >= 0; v--) 
    {
      if(cm->stid[v] != BEGL_S) 
	{
	  alpha[0][v][0] = IMPOSSIBLE;
	  if      (cm->sttype[v] == E_st)  { 
	    alpha[0][v][0] = alpha[1][v][0] = 0.;
	    /* rest of E deck is IMPOSSIBLE, this rewritten if QDB is on, (slightly wasteful). */
	    for (d = 1; d <= W; d++) alpha[0][v][d] = alpha[1][v][d] = IMPOSSIBLE;
	  }
	  else if (cm->sttype[v] == MP_st) alpha[0][v][1] = alpha[1][v][1] = IMPOSSIBLE;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[0][v][0] = cm->endsc[v];
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[0][v][0] = ESL_MAX(alpha[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	      alpha[0][v][0] = ESL_MAX(alpha[0][v][0], IMPOSSIBLE);
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* BEGL_S, left child state */
	      y = cm->cnum[v];
	      alpha[0][v][0] = alpha_begl[0][w][0] + alpha[0][y][0]; 
	    }

	  alpha[1][v][0] = alpha[0][v][0];
	}
      else /* v == BEGL_S */
	{
	  y = cm->cfirst[v];
	  alpha_begl[0][v][0] = cm->endsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    alpha_begl[0][v][0] = ESL_MAX(alpha_begl[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in alpha */
	  alpha_begl[0][v][0] = ESL_MAX(alpha_begl[0][v][0], IMPOSSIBLE);
	  for (j = 1; j <= W; j++) 
	    alpha_begl[j][v][0] = alpha_begl[0][v][0];
	}
    }
      bestr[0] = -1;

  /*
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   */
  if(results != NULL) { 
    ESL_ALLOC(gamma,  sizeof(float) * (L+1));
    gamma[0] = 0;
    ESL_ALLOC(gback,  sizeof(int)   * (L+1));
    gback[0] = -1;
    ESL_ALLOC(savesc, sizeof(float) * (L+1));
    ESL_ALLOC(saver,  sizeof(int)   * (L+1));
  }
  /*
   * query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) {
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < dmin[v] && d <=W; d++) 
	  for(j = 0; j < 2; j++)
	    alpha[j][v][d] = IMPOSSIBLE;
	for (d = dmax[v]+1; d <= W;      d++) 
	  for(j = 0; j < 2; j++)
	    alpha[j][v][d] = IMPOSSIBLE;
      }
      else
	{
	  for (d = 0; d < dmin[v] && d <=W; d++) 
	    for(j = 0; j <= W; j++)
	      alpha_begl[j][v][d] = IMPOSSIBLE;
	  for (d = dmax[v]+1; d <= W;      d++) 
	    for(j = 0; j <= W; j++)
	      alpha_begl[j][v][d] = IMPOSSIBLE;
	}
    }
  }

  /* precalculate minimum and maximum d for each state and each sequence index (1..j..W). 
   * this is not always just dmin, dmax, (for ex. if j < W).
   */
  ESL_ALLOC(dnAA, sizeof(int *) * (W+1));
  ESL_ALLOC(dxAA, sizeof(int *) * (W+1));
  
  dnAA[0] = dxAA[0] = NULL; /* corresponds to j == 0, which is out of bounds */
  for(j = 1; j <= W; j++) {
    ESL_ALLOC(dnAA[j], sizeof(int) * cm->M);
    ESL_ALLOC(dxAA[j], sizeof(int) * cm->M);

    for(v = 0; v < cm->M; v++) {
      if(do_banded) { 
	dnAA[j][v] = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	dxAA[j][v] = ESL_MIN(j, dmax[v]); 
	dxAA[j][v] = ESL_MIN(dxAA[j][v], W);
      }
      else { 
	dnAA[j][v] = (cm->sttype[v] == MP_st) ? 2 : 1;
	dxAA[j][v] = ESL_MIN(j, W); 
      }
    }
  }

  /* precalculate possible emission scores for each state */
  float **esc_vAA;
  int a,b;
  ESL_ALLOC(esc_vAA, sizeof(float *) * (cm->M));
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      ESL_ALLOC(esc_vAA[v], sizeof(float) * cm->abc->Kp);
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v], cm->abc->Kp, IMPOSSIBLE);
      for(a = 0; a < cm->abc->K; a++)
	esc_vAA[v][a] = cm->esc[v][a];
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++)
	esc_vAA[v][a] = esl_abc_FAvgScore(cm->abc, a, cm->esc[v]);
      break;
    case MP_st:
      ESL_ALLOC(esc_vAA[v], sizeof(float) * (cm->abc->Kp * cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v], (cm->abc->Kp * cm->abc->Kp), IMPOSSIBLE);
      for(a = 0; a < (cm->abc->Kp-1); a++)
	for(b = 0; b < (cm->abc->Kp-1); b++)
	  if(a < cm->abc->K && b < cm->abc->K)
	    esc_vAA[v][a * cm->abc->Kp + b] = cm->esc[v][(a * cm->abc->K) + b];
	  else
	    esc_vAA[v][a * cm->abc->Kp + b] = DegeneratePairScore(cm->abc, cm->esc[v], a, b);
      break;
    default:
      esc_vAA[v] = NULL;
      break;
    }
  }

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* precalculate the initial score for all alpha[v][j][d] cells, it's independent
   * of j, so we do it here, outside the for(j...) loop */
  ESL_ALLOC(init_scAA, sizeof(float *) * (cm->M));
  for (v = 0; v < cm->M; v++) 
    {
      ESL_ALLOC(init_scAA[v], sizeof(float) * (W+1));
      if(NOT_IMPOSSIBLE(cm->endsc[v]))
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = el_scA[d] + cm->endsc[v];
      else
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = IMPOSSIBLE;
    }

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Precalculate the 'emit mode' of each state to speed up the addition of emission 
   * scores, all states are either EMITLEFT, EMITRIGHT, EMITPAIR, or EMITNONE, this
   * collapses ILs and MLs into 1 value (for example) for the switch() statement inside the for(d...) loop
   * in the heart of the recursion which saves us time.
   */
  int *emitmodeA;
  ESL_ALLOC(emitmodeA, sizeof(int) * cm->M);
  for(v = 0; v < cm->M; v++) {
    switch (cm->sttype[v]) {
    case IL_st:
    case ML_st:
      emitmodeA[v] = EMITLEFT;
      break;
    case IR_st:
    case MR_st:
      emitmodeA[v] = EMITRIGHT;
      break;		
    case MP_st:
      emitmodeA[v] = EMITPAIR;
      break;		
    default:
      emitmodeA[v] = EMITNONE;
      break;
    }
  }
  /* Initialize sc_v to size of M */
  float *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);
  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      float tsc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else {          dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  /* float const *esc_v = cm->esc[v]; */
	  float const *esc_v = esc_vAA[v]; 
	  float const *tsc_v = cm->tsc[v];
	  //float sc;
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

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = ESL_MAX(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = ESL_MAX (sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
	    /******************************************************************************/
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
	     * nested emitmodeA[v] switch is based on empirical
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
		/* ctr++; */
		sc = ESL_MAX(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc = ESL_MAX(sc, init_scAA[v][dp_y]);
		sc = ESL_MAX(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
		case EMITLEFT:
		  sc += esc_v[dsq[i--]];
		  break;
		case EMITNONE:
		  break;
		case EMITRIGHT:
		  sc += esc_j;
		  break;		
		case EMITPAIR:
 		  sc += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
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
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
		case EMITLEFT:
		  sc += esc_v[dsq[i--]];
		  break;
		case EMITNONE:
		  break;
		case EMITRIGHT:
		  sc += esc_j;
		  break;		
		case EMITPAIR: 
		  sc += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
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
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
		case EMITLEFT:
		  sc += esc_v[dsq[i--]];
		  break;
		case EMITNONE:
		  break;
		case EMITRIGHT:
		  sc += esc_j;
		  break;		
		case EMITPAIR: /* OPTIMIZE? */
		  sc += esc_v[dsq[i--]*cm->abc->Kp+dsq[j]];
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

		/* add in emission score, if any */
		switch (emitmodeA[v]) {
		case EMITRIGHT:
		  sc += esc_j;
		  break;		
		case EMITNONE:
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
		switch (emitmodeA[v]) {
		case EMITLEFT:
		  sc += esc_v[dsq[i--]];
		  break;
		case EMITNONE:
		  break;
		case EMITRIGHT:
		  sc += esc_j;
		  break;		
		} /* end of switch (cm->sttype[v]) */
		alpha[jp_v][v][d] = sc;
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;
	    }
	  }
      else { /* enter else if cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st */
	    y    = cm->cfirst[v];
	    dp_y = dn - sd; /* initial dp_y, we increment it at end of 'for(d = ...' loop */
	    i    = j-dn+1;  /* initial i,    we decrement it when we access it, inside each possible case of the switch (cnum) below */
	    for (d = dn; d <= dx; d++, dp_y++) 
	      sc_v[d] = init_scAA[v][dp_y];
	    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++, y++) {
	      tsc  = cm->tsc[v][yoffset];
	      for (d = dn, dp_y = dn-sd; d <= dx; d++, dp_y++) {
		sc_v[d] = ESL_MAX (sc_v[d], alpha[jp_y][y][dp_y] + tsc);
	      }
	    }
	    /* add in emission score (if any), and set alpha[jp_v][v][d] cell */
	    switch (emitmodeA[v]) {
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
	    } /* end of switch (emitmodeA[v]) */
	  } /* end of else (cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st) */
	  if(vsc != NULL) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
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
	bestr[d] = 0;	/* root of the traceback = root state 0 */
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + tsc_v[0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = ESL_MAX (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + tsc_v[yoffset]));
      }
	
      if (cm->flags & CM_LOCAL_BEGIN) {
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
      if(results != NULL) 
	{
	  float sc;
	  /* get information on hits */
	  if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	    {
	      /* The little semi-Markov model that deals with multihit parsing:
	       */
	      gamma[jp_g]  = gamma[jp_g-1] + 0; 
	      gback[jp_g]  = -1;
	      savesc[jp_g] = IMPOSSIBLE;
	      saver[jp_g]  = -1;
	      for (d = dnA[0]; d <= dxA[0]; d++) {
		i    = j-d+1;
		ip_g = i-i0+1;
		sc = gamma[ip_g-1] + alpha[jp_v][0][d] + cm->sc_boost; 
		/* sc_boost is experimental technique for finding hits < 0 bits. value is 0.0 if technique not used. */
		if (sc > gamma[jp_g]) {
		  gamma[jp_g]  = sc;
		  gback[jp_g]  = i;
		  savesc[jp_g] = alpha[jp_v][0][d]; 
		  saver[jp_g]  = bestr[d];
		}
	      }
	    }
	  else {
	    /* Resolving overlaps greedily (RSEARCH style),  
	     * At least one hit is sent back for each j here.
	     * However, some hits can already be removed for the greedy overlap
	     * resolution algorithm.  Specifically, at the given j, any hit with a
	     * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
	    /* First, report hit with d of 1 if > cutoff */
	    if (alpha[jp_v][0][1] >= cutoff) 
	      if(results != NULL) 
		report_hit (j, j, bestr[1], alpha[jp_v][0][1], results);
	    bestd = 1;
	    if (alpha[jp_v][0][1] > best_hit_sc)
	      best_hit_sc = alpha[jp_v][0][1];
	    
	    /* Now, if current score is greater than maximum seen previous, report
	       it if >= cutoff and set new max */
	    for (d = 2; d <= W && d <= jp_g; d++) 
	      {
		if (alpha[jp_v][0][d] > best_hit_sc) best_hit_sc = alpha[jp_v][0][d];
		if (alpha[jp_v][0][d] > alpha[jp_v][0][bestd]) {
		  if (alpha[jp_v][0][d] >= cutoff)
		    if(results != NULL) 
		      report_hit (j-d+1, j, bestr[d], alpha[jp_v][0][d], results);
		  bestd = d;
		}
	      }
	  }
	}
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;
  
  /* free alpha and alpha_begl, everything we need is in gamma.
   */

  free(alpha[0][0]);
  free(alpha[1]);
  free(alpha[0]);
  free(alpha);
  for (j = 0; j <= W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S)
	free(alpha_begl[j][v]);
    free(alpha_begl[j]);
  }
  free(alpha_begl);
  free(bestr);

  if(results != NULL) 
    {
      /*
       * Traceback stage.
       * Recover all hits: an (i,j,sc) triple for each one.
       */
      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	{
	  j     = j0;
	  while (j >= i0) 
	    {
	      jp_g = j-i0+1;
	      if (gback[jp_g] == -1) /* no hit */
		j--; 
	      else                /* a hit, a palpable hit */
		{
		  if(savesc[jp_g] > best_hit_sc) best_hit_sc = savesc[jp_g];
		  if(savesc[jp_g] >= cutoff && results != NULL) /* report the hit */
		    report_hit(gback[jp_g], j, saver[jp_g], savesc[jp_g], results);
		  j = gback[jp_g]-1;
		}
	    }
	}
      free(gback);
      free(gamma);
      free(savesc);
      free(saver);
    }

  for(v = 0; v < cm->M; v++) {
    free(init_scAA[v]);
    if(esc_vAA[v] != NULL) free(esc_vAA[v]);
  }
  free(init_scAA);
  free(esc_vAA);

  for(j = 1; j <= W; j++) {
    free(dnAA[j]);
    free(dxAA[j]);
  }
  free(dnAA);
  free(dxAA);
  free(jp_wA);
  free(el_scA);
  free(emitmodeA);
  if (ret_best_hit_sc != NULL) *ret_best_hit_sc = best_hit_sc;
  if (ret_vsc         != NULL) *ret_vsc         = vsc;
  
  /* printf("XFastCYKScan() return score: %10.4f\n", vsc_root); 
     printf("ctr: %d\n", ctr);*/
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
    if (cm->flags & CM_LOCAL_BEGIN) {
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
    
    if (cm->flags & CM_LOCAL_BEGIN) {
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


/* Function: cm_CreateScanInfo()
 * Date:     EPN, Sun Nov  4 19:56:58 2007
 *
 * Purpose:  Given a CM, allocate and initialize ScanInfo_t object for that CM. 
 *            
 * Returns:  Newly allocated HybridScanInfo_t object:
 */
ScanInfo_t *
cm_CreateScanInfo(CM_t *cm)
{
  int status;
  int d, y, yoffset, w, j;
  int v;
  int do_banded;

  if(cm->dmin == NULL && cm->dmax != NULL) cm_Fail("cm_CreateScanInfo(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmin != NULL && cm->dmax == NULL) cm_Fail("cm_CreateScanInfo(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmax != NULL && cm->W != cm->dmax[0]) cm_Fail("cm_CreateScanInfo(), cm->W: %d != cm->dmax[0]: %d\n", cm->W, cm->dmax[0]); 
  if((! cm->search_opts & CM_SEARCH_NOQDB) && (cm->dmin == NULL || cm->dmax == NULL))
    cm_Fail("cm_CreateScanInfo(), cm->dmin == NULL || cm->dmax == NULL, but !(cm->search_opts & CM_SEARCH_NOQDB)\n");

  ScanInfo_t *si;
  ESL_ALLOC(si, sizeof(ScanInfo_t));

  si->cm_M  = cm->M;
  si->W     = cm->W;
  si->dmin  = cm->dmin; /* could be NULL */
  si->dmax  = cm->dmax; /* could be NULL */
  do_banded = (cm->search_opts & CM_SEARCH_NOQDB) ? FALSE : TRUE;

  /* fill emitmodeA */
  ESL_ALLOC(si->emitmodeA, sizeof(int) * cm->M);
  for(v = 0; v < cm->M; v++) {
    switch (cm->sttype[v]) {
    case IL_st:
    case ML_st:
      si->emitmodeA[v] = EMITLEFT;
      break;
    case IR_st:
    case MR_st:
      si->emitmodeA[v] = EMITRIGHT;
      break;		
    case MP_st:
      si->emitmodeA[v] = EMITPAIR;
      break;		
    default:
      si->emitmodeA[v] = EMITNONE;
      break;
    }
  }

  /* precalculate minimum and maximum d for each state and each sequence index (1..j..W). 
   * this is not always just dmin, dmax, (for ex. if j < W). */
  int W = si->W;
  ESL_ALLOC(si->dnAA, sizeof(int *) * (W+1));
  ESL_ALLOC(si->dxAA, sizeof(int *) * (W+1));
  si->dnAA[0] = si->dxAA[0] = NULL; /* corresponds to j == 0, which is out of bounds */
  for(j = 1; j <= W; j++) {
    ESL_ALLOC(si->dnAA[j], sizeof(int) * cm->M);
    ESL_ALLOC(si->dxAA[j], sizeof(int) * cm->M);
    for(v = 0; v < cm->M; v++) {
      if(do_banded) { 
	si->dnAA[j][v] = (cm->sttype[v] == MP_st) ? ESL_MAX(si->dmin[v], 2) : ESL_MAX(si->dmin[v], 1); 
	si->dxAA[j][v] = ESL_MIN(j, si->dmax[v]); 
	si->dxAA[j][v] = ESL_MIN(si->dxAA[j][v], W);
      }
      else { 
	si->dnAA[j][v] = (cm->sttype[v] == MP_st) ? 2 : 1;
	si->dxAA[j][v] = ESL_MIN(j, W); 
      }
    }
  }

  /* precalculate possible emission scores for each state */
  int a,b;
  ESL_ALLOC(si->esc_vAA,  sizeof(float *) * (cm->M));
  ESL_ALLOC(si->iesc_vAA, sizeof(int *)   * (cm->M));
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      ESL_ALLOC(si->esc_vAA[v],  sizeof(float) * (cm->abc->Kp));
      ESL_ALLOC(si->iesc_vAA[v], sizeof(int)   * (cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(si->esc_vAA[v],  cm->abc->Kp, IMPOSSIBLE);
      esl_vec_ISet(si->iesc_vAA[v], cm->abc->Kp, -INFTY);
      for(a = 0; a < cm->abc->K; a++) {
	si->esc_vAA[v][a]  = cm->esc[v][a];
	si->iesc_vAA[v][a] = cm->iesc[v][a];
      }
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) { /* note boundary conditions, gap, missing data symbols stay IMPOSSIBLE */
	si->esc_vAA[v][a]  = esl_abc_FAvgScore(cm->abc, a, cm->esc[v]);
	si->iesc_vAA[v][a] = esl_abc_IAvgScore(cm->abc, a, cm->iesc[v]);
      }
      break;
    case MP_st:
      ESL_ALLOC(si->esc_vAA[v], sizeof(float) * (cm->abc->Kp * cm->abc->Kp));
      ESL_ALLOC(si->iesc_vAA[v], sizeof(int) * (cm->abc->Kp * cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(si->esc_vAA[v],  (cm->abc->Kp * cm->abc->Kp), IMPOSSIBLE);
      esl_vec_ISet(si->iesc_vAA[v], (cm->abc->Kp * cm->abc->Kp), -INFTY);
      for(a = 0; a < (cm->abc->Kp-1); a++)
	for(b = 0; b < (cm->abc->Kp-1); b++)
	  if(a < cm->abc->K && b < cm->abc->K) {
	    si->esc_vAA[v][a * cm->abc->Kp + b]  = cm->esc[v][(a * cm->abc->K) + b];
	    si->iesc_vAA[v][a * cm->abc->Kp + b] = cm->iesc[v][(a * cm->abc->K) + b];
	  }
	  else {
	    si->esc_vAA[v][a  * cm->abc->Kp + b]  = DegeneratePairScore(cm->abc, cm->esc[v], a, b);
	    si->iesc_vAA[v][a * cm->abc->Kp + b] = iDegeneratePairScore(cm->abc, cm->iesc[v], a, b);
	  }
      break;
    default:
      si->esc_vAA[v] = NULL;
      si->iesc_vAA[v] = NULL;
      break;
    }
  }

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(si->el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) si->el_scA[d] = cm->el_selfsc * d;
  ESL_ALLOC(si->iel_scA, sizeof(int) * (W+1));
  for(d = 0; d <= W; d++) si->iel_scA[d] = cm->iel_selfsc * d;

  /* precalculate the initial score for all alpha[v][j][d] cells, it's independent of j */
  ESL_ALLOC(si->init_scAA, sizeof(float *) * (cm->M));
  for (v = 0; v < cm->M; v++) {
    ESL_ALLOC(si->init_scAA[v], sizeof(float) * (W+1));
    if(NOT_IMPOSSIBLE(cm->endsc[v]))
      for(d = 0; d <= W; d++)
	si->init_scAA[v][d] = si->el_scA[d] + cm->endsc[v];
    else
      for(d = 0; d <= W; d++)
	si->init_scAA[v][d] = IMPOSSIBLE;
  }
  ESL_ALLOC(si->iinit_scAA, sizeof(int *) * (cm->M));
  for (v = 0; v < cm->M; v++) {
    ESL_ALLOC(si->iinit_scAA[v], sizeof(int) * (W+1));
    if(NOT_IMPOSSIBLE(cm->endsc[v])) /* we use endsc[v], not iendsc[v] */
      for(d = 0; d <= W; d++)
	si->iinit_scAA[v][d] = si->iel_scA[d] + cm->iendsc[v];
    else
      for(d = 0; d <= W; d++)
	si->iinit_scAA[v][d] = -INFTY;
  }
  /* allocate bestr, which holds best root state at alpha[0][cur][d] */
  ESL_ALLOC(si->bestr, (sizeof(int) * (W+1)));

  /* alpha, alpha_begl, ialpha, and ialpha_begl allocations, only difference is:
   *  alpha,  alpha_begl are floats
   * ialpha, ialpha_begl are ints  
   *
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
   */

  /* allocate alpha */
  ESL_ALLOC(si->alpha,        sizeof(float **) * 2);
  ESL_ALLOC(si->alpha[0],     sizeof(float *) * cm->M);
  ESL_ALLOC(si->alpha[1],     sizeof(float *) * cm->M);
  ESL_ALLOC(si->alpha[0][0],  sizeof(float) * 2 * (cm->M) * (W+1));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      si->alpha[0][v] = si->alpha[0][0] + (v           * (W+1));
      si->alpha[1][v] = si->alpha[0][0] + ((v + cm->M) * (W+1));
    }
    else si->alpha[0][v] = si->alpha[1][v] = NULL; /* BEGL_S */
  }
  /* allocate ialpha */
  ESL_ALLOC(si->ialpha,        sizeof(int **) * 2);
  ESL_ALLOC(si->ialpha[0],     sizeof(int *) * cm->M);
  ESL_ALLOC(si->ialpha[1],     sizeof(int *) * cm->M);
  ESL_ALLOC(si->ialpha[0][0],  sizeof(int) * 2 * (cm->M) * (W+1));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      si->ialpha[0][v] = si->ialpha[0][0] + (v           * (W+1));
      si->ialpha[1][v] = si->ialpha[0][0] + ((v + cm->M) * (W+1));
    }
    else si->ialpha[0][v] = si->ialpha[1][v] = NULL; /* BEGL_S */
  }
  /* allocate alpha_begl */
  ESL_ALLOC(si->alpha_begl, (sizeof(float **) * (W+1)));
  for (j = 0; j <= W; j++) {
    ESL_ALLOC(si->alpha_begl[j], (sizeof(float *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) ESL_ALLOC(si->alpha_begl[j][v], (sizeof(float) * (W+1)));
      else si->alpha_begl[j][v] = NULL; /* non-BEGL */
    }
  }
  /* allocate ialpha_begl */
  ESL_ALLOC(si->ialpha_begl, (sizeof(int **) * (W+1)));
  for (j = 0; j <= W; j++) {
    ESL_ALLOC(si->ialpha_begl[j], (sizeof(int *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) ESL_ALLOC(si->ialpha_begl[j][v], (sizeof(int) * (W+1)));
      else si->ialpha_begl[j][v] = NULL; /* non-BEGL */
    }
  }

  /* alpha initializations.
   * First initialize alpha and  alpha_begl (float matrices)
   * then initialize ialpha and ialpha_begl (int matrices)
   *
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   */

  float ***alpha      = si->alpha;
  float ***alpha_begl = si->alpha_begl;
  int  ***ialpha      = si->ialpha;
  int  ***ialpha_begl = si->ialpha_begl;

  /* initialize alpha and alpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      alpha[0][v][0] = IMPOSSIBLE;
      if (cm->sttype[v] == E_st) { 
	alpha[0][v][0] = alpha[1][v][0] = 0.;
	/* rest of E deck is IMPOSSIBLE, this is rewritten if QDB is on, (slightly wasteful). */
	for (d = 1; d <= W; d++) alpha[0][v][d] = alpha[1][v][d] = IMPOSSIBLE;
      }
      else if (cm->sttype[v] == MP_st) alpha[0][v][1] = alpha[1][v][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	alpha[0][v][0] = cm->endsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  alpha[0][v][0] = ESL_MAX(alpha[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	alpha[0][v][0] = ESL_MAX(alpha[0][v][0], IMPOSSIBLE);
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	alpha[0][v][0] = alpha_begl[0][w][0] + alpha[0][y][0]; 
      }
      alpha[1][v][0] = alpha[0][v][0];
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      alpha_begl[0][v][0] = cm->endsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	alpha_begl[0][v][0] = ESL_MAX(alpha_begl[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in alpha */
      alpha_begl[0][v][0] = ESL_MAX(alpha_begl[0][v][0], IMPOSSIBLE);
      for (j = 1; j <= W; j++) 
	alpha_begl[j][v][0] = alpha_begl[0][v][0];
    }
  }

  /* initialize ialpha and ialpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      ialpha[0][v][0] = -INFTY;
      if (cm->sttype[v] == E_st) { 
	ialpha[0][v][0] = ialpha[1][v][0] = 0.;
	/* rest of E deck is -INFTY, this is rewritten if QDB is on, (slightly wasteful). */
	for (d = 1; d <= W; d++) ialpha[0][v][d] = ialpha[1][v][d] = -INFTY;
      }
      else if (cm->sttype[v] == MP_st) ialpha[0][v][1] = ialpha[1][v][1] = -INFTY;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	ialpha[0][v][0] = cm->iendsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  ialpha[0][v][0] = ESL_MAX(ialpha[0][v][0], (ialpha[0][y+yoffset][0] + cm->itsc[v][yoffset]));
	ialpha[0][v][0] = ESL_MAX(ialpha[0][v][0], -INFTY);
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	ialpha[0][v][0] = ialpha_begl[0][w][0] + ialpha[0][y][0]; 
      }
      ialpha[1][v][0] = ialpha[0][v][0];
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      ialpha_begl[0][v][0] = cm->iendsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	ialpha_begl[0][v][0] = ESL_MAX(ialpha_begl[0][v][0], (ialpha[0][y+yoffset][0] + cm->itsc[v][yoffset])); /* careful: y is in alpha */
      ialpha_begl[0][v][0] = ESL_MAX(ialpha_begl[0][v][0], -INFTY);
      for (j = 1; j <= W; j++) 
	ialpha_begl[j][v][0] = ialpha_begl[0][v][0];
    }
  }

  /* query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) {
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < cm->dmin[v] && d <=W; d++) {
	  for(j = 0; j < 2; j++) {
	    alpha[j][v][d]  = IMPOSSIBLE;
	    ialpha[j][v][d] = -INFTY;
	  }
	}
	for (d = cm->dmax[v]+1; d <= W;      d++) { 
	  for(j = 0; j < 2; j++) {
	    alpha[j][v][d] = IMPOSSIBLE;
	    ialpha[j][v][d] = -INFTY;
	  }
	}
      }
      else {
	for (d = 0; d < cm->dmin[v] && d <=W; d++) {
	  for(j = 0; j <= W; j++) {
	       alpha_begl[j][v][d] = IMPOSSIBLE;
	      ialpha_begl[j][v][d] = -INFTY;
	  }
	}
	for (d = cm->dmax[v]+1; d <= W;      d++) {
	  for(j = 0; j <= W; j++) {
	    alpha_begl[j][v][d] = IMPOSSIBLE;
	    ialpha_begl[j][v][d] = -INFTY;
	  }
	}
      }
    }
  }
  return si;

 ERROR:
  cm_Fail("memory allocation error in cm_CreateScanInfo().\n");
  return NULL;
}


/* Function: cm_FreeScanInfo()
 * Date:     EPN, Sun Nov  4 20:57:32 2007
 *
 * Purpose:  Free a ScanInfo_t object for <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeScanInfo(CM_t *cm, ScanInfo_t *si)
{
  int j, v;
  for(v = 0; v < cm->M; v++) {
    free(si->init_scAA[v]);
    free(si->iinit_scAA[v]);
    if(si->esc_vAA[v] != NULL)  free(si->esc_vAA[v]);
    if(si->iesc_vAA[v] != NULL) free(si->iesc_vAA[v]);
  }
  free(si->init_scAA);
  free(si->esc_vAA);
  free(si->iinit_scAA);
  free(si->iesc_vAA);

  for(j = 1; j <= si->W; j++) {
    free(si->dnAA[j]);
    free(si->dxAA[j]);
  }
  free(si->dnAA);
  free(si->dxAA);
  free(si->el_scA);
  free(si->iel_scA);
  free(si->emitmodeA);

  free(si->bestr);

  /* free alpha and alpha_begl */
  free(si->alpha[0][0]);
  free(si->alpha[1]);
  free(si->alpha[0]);
  free(si->alpha);
  for (j = 0; j <= si->W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S)
	free(si->alpha_begl[j][v]);
    free(si->alpha_begl[j]);
  }
  free(si->alpha_begl);

  free(si->ialpha[0][0]);
  free(si->ialpha[1]);
  free(si->ialpha[0]);
  free(si->ialpha);
  for (j = 0; j <= si->W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S)
	free(si->ialpha_begl[j][v]);
    free(si->ialpha_begl[j]);
  }
  free(si->ialpha_begl);
  free(si);

  return eslOK;
}


/* Function: cm_DumpScanInfoAlpha()
 * Date:     EPN, Tue Nov  6 05:11:26 2007
 *
 * Purpose:  Dump current alpha matrix (either float or int).
 *            
 * Returns:  void.
 */
void
cm_DumpScanInfoAlpha(CM_t *cm, ScanInfo_t *si, int j, int i0, int doing_float)
{
  int d, v;
  int jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
  int cur = j%2;
  int prv = (j-1)%2;
  int *dnA, *dxA;
  int begl_prv = j-1 % (si->W+1);
  int begl_cur = j   % (si->W+1);
  printf("Dumping Alpha: j: %d\n", j);
  if(jp_g >= si->W) { dnA = si->dnAA[si->W]; dxA = si->dxAA[si->W]; }
  else              { dnA = si->dnAA[jp_g];  dxA = si->dxAA[jp_g]; }
  if(doing_float) { 
    for (v = si->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) { 
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10.4f\n", (j-1), v, d, si->alpha_begl[begl_prv][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10.4f\n", (j-1), v, d, si->alpha[prv][v][d]); 
      }
      printf("\n");
    }
    for (v = si->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10.4f\n", j,     v, d, si->alpha_begl[begl_cur][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10.4f\n", j,     v, d, si->alpha[cur][v][d]); 
      }
      printf("\n");
    }
  }
  else { /* doing int */
    for (v = si->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10d\n", (j-1), v, d, si->ialpha_begl[begl_prv][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10d\n", (j-1), v, d, si->ialpha[prv][v][d]); 
      }
      printf("\n\n");
    }
    for (v = si->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10d\n", j,     v, d, si->ialpha_begl[begl_cur][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10d\n", j,     v, d, si->ialpha[cur][v][d]); 
      }
      printf("\n\n");
    }
  }
  return;
}
  
/* Function: cm_CreateGammaHitMx()
 * Date:     EPN, Mon Nov  5 05:22:56 2007
 *
 * Purpose:  Allocate and initialize a gamma semi-HMM for 
 *           optimal hit resolution of a CM based scan.
 *            
 * Returns:  Newly allocated CMGammaHitMx_t object:
 */
cm_GammaHitMx_t *
cm_CreateGammaHitMx(int L, int i0, int be_greedy, float cutoff)
{
  int status;

  cm_GammaHitMx_t *gamma;
  ESL_ALLOC(gamma, sizeof(cm_GammaHitMx_t));

  gamma->L  = L;
  gamma->i0 = i0;
  gamma->iamgreedy = be_greedy;
  gamma->cutoff    = cutoff;
  /* allocate/initialize for CYK/Inside */
  ESL_ALLOC(gamma->mx,     sizeof(float) * (L+1));
  ESL_ALLOC(gamma->gback,  sizeof(int)   * (L+1));
  ESL_ALLOC(gamma->savesc, sizeof(float) * (L+1));
  ESL_ALLOC(gamma->saver,  sizeof(int)   * (L+1));
    
  gamma->mx[0]    = 0;
  gamma->gback[0] = -1;
  return gamma;

 ERROR:
  cm_Fail("memory allocation error in cm_CreateGammaHitMx().\n");
  return NULL;
}

/* Function: cm_FreeGammaHitMx()
 * Date:     EPN, Mon Nov  5 05:32:00 2007
 *
 * Purpose:  Free a gamma semi-HMM.
 *            
 * Returns:  void;
 */
void
cm_FreeGammaHitMx(cm_GammaHitMx_t *gamma)
{
  free(gamma->mx);
  free(gamma->gback);
  free(gamma->savesc);
  free(gamma->saver);
  free(gamma);

  return;
}

/* Function: cm_UpdateFloatGammaHitMx()
 * Date:     EPN, Mon Nov  5 05:41:14 2007
 *
 * Purpose:  Update a gamma semi-HMM for CM hits that end at gamma-relative position <j>.
 *            
 * Returns:  void;
 */
void
cm_UpdateFloatGammaHitMx(cm_GammaHitMx_t *gamma, int j, float *alpha_row, int dn, int dx, int *bestr, float sc_boost, 
			 int doing_inside, search_results_t *results)
{
  int i, d;
  float sc;
  int bestd;
  int r;

  /* mode 1: non-greedy  */
  if(! gamma->iamgreedy) { 
    gamma->mx[j]     = gamma->mx[j-1] + 0; 
    gamma->gback[j]  = -1;
    gamma->savesc[j] = IMPOSSIBLE;
    gamma->saver[j]  = -1;
    for (d = dn; d <= dx; d++) {
      i  = j-d+1;
      sc = gamma->mx[i-1] + alpha_row[d] + sc_boost; 
      /* sc_boost is experimental technique for finding hits < 0 bits. value is 0.0 if technique not used. */
      if (sc > gamma->mx[j]) {
	gamma->mx[j]     = sc;
	gamma->gback[j]  = i;
	gamma->savesc[j] = alpha_row[d]; 
	gamma->saver[j]  = doing_inside ? -1 : bestr[d]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
      }
    }
  }
  /* mode 2: greedy */
  if(gamma->iamgreedy) { 
    /* Resolving overlaps greedily (RSEARCH style),  
     * At least one hit is sent back for each j here.
     * However, some hits can already be removed for the greedy overlap
     * resolution algorithm.  Specifically, at the given j, any hit with a
     * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
    /* First, report hit with d of dn (min valid d) if >= cutoff */
    if (alpha_row[dn] >= gamma->cutoff) {
      r = doing_inside ? -1 : bestr[dn]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
      report_hit (j-dn+gamma->i0, j-dn+gamma->i0, r, alpha_row[dn], results);
    }
    bestd = dn;
    /* Now, if current score is greater than maximum seen previous, report
       it if >= cutoff and set new max */
    for (d = dn+1; dn <= dx; d++) {
      if (alpha_row[d] > alpha_row[bestd]) {
	if (alpha_row[d] >= gamma->cutoff) { 
	  r = doing_inside ? -1 : bestr[d]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
	  report_hit (j-d+gamma->i0, j-d+gamma->i0, r, alpha_row[d], results);
	}
	bestd = d;
      }
    }
  }
  return;
}

/* Function: cm_UpdateIntGammaHitMx()
 * Date:     EPN, Tue Nov  6 05:54:51 2007
 *
 * Purpose:  Update a gamma semi-HMM for CM hits that end at gamma-relative position <j>.
 *            
 * Returns:  void;
 */
void
cm_UpdateIntGammaHitMx(cm_GammaHitMx_t *gamma, int j, int *alpha_row, int dn, int dx, int *bestr, float sc_boost, 
		       int doing_inside, search_results_t *results)
{
  int i, d;
  float sc;
  float bestsc;
  int r;

  /* mode 1: non-greedy  */
  if(! gamma->iamgreedy) { 
    gamma->mx[j]     = gamma->mx[j-1] + 0; 
    gamma->gback[j]  = -1;
    gamma->savesc[j] = IMPOSSIBLE;
    gamma->saver[j]  = -1;
    for (d = dn; d <= dx; d++) {
      i  = j-d+1;
      sc = gamma->mx[i-1] + Scorify(alpha_row[d]) + sc_boost; 
      /* sc_boost is experimental technique for finding hits < 0 bits. value is 0.0 if technique not used. */
      if (sc > gamma->mx[j]) {
	gamma->mx[j]     = sc;
	gamma->gback[j]  = i;
	gamma->savesc[j] = Scorify(alpha_row[d]); 
	gamma->saver[j]  = doing_inside ? -1 : bestr[d]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
      }
    }
  }
  /* mode 2: greedy */
  if(gamma->iamgreedy) { 
    /* Resolving overlaps greedily (RSEARCH style),  
     * At least one hit is sent back for each j here.
     * However, some hits can already be removed for the greedy overlap
     * resolution algorithm.  Specifically, at the given j, any hit with a
     * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
    /* First, report hit with d of dn (min valid d) if >= cutoff */
    if (Scorify(alpha_row[dn]) >= gamma->cutoff) { 
      r = doing_inside ? -1 : bestr[dn]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
      report_hit (j-dn+gamma->i0, j-dn+gamma->i0, r, Scorify(alpha_row[dn]), results);
    }
    bestsc = Scorify(alpha_row[dn]);
    /* Now, if current score is greater than maximum seen previous, report
       it if >= cutoff and set new max */
    for (d = dn+1; dn <= dx; d++) {
      sc = Scorify(alpha_row[d]);
      if (sc > bestsc) {
	if (sc >= gamma->cutoff) { 
	  r = doing_inside ? -1 : bestr[dn]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
	  report_hit (j-d+gamma->i0, j-d+gamma->i0, r, sc, results);
	}
	bestsc = sc;
      }
    }
  }
  return;
}


/* Function: cm_TBackGammaHitMx()
 * Date:     EPN, Mon Nov  5 10:14:30 2007
 *
 * Purpose:  Traceback with a gamma semi-HMM for CM hits. gamma->iamgreedy should be FALSE.
 *            
 * Returns:  void; dies immediately upon an error.
 */
void
cm_TBackGammaHitMx(cm_GammaHitMx_t *gamma, search_results_t *results, int i0, int j0)
{
  int j, jp_g;

  if(gamma->iamgreedy) cm_Fail("cm_TBackGammaHitMx(), gamma->iamgreedy is TRUE.\n");   
  if(results == NULL) cm_Fail("cm_TBackGammaHitMx(), results == NULL");
  /*
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   */
  j = j0;
  while (j >= i0) {
    jp_g = j-i0+1;
    if (gamma->gback[jp_g] == -1) /* no hit */
      j--; 
    else {              /* a hit, a palpable hit */
      if(gamma->savesc[jp_g] >= gamma->cutoff) /* report the hit */
	report_hit(gamma->gback[jp_g], j, gamma->saver[jp_g], gamma->savesc[jp_g], results);
      j = gamma->gback[jp_g]-1;
    }
  }
  return;
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
  { "-L",        eslARG_INT,  "10000", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "search in glocal mode [default: local]", 0 },
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
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "benchmark driver for an optimized scanning CYK implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  CM_t            *cm;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  char            *cmfile = esl_opt_GetArg(go, 1);
  CMFILE          *cmfp;	/* open input CM file stream */
  int            *dmin;
  int            *dmax;
  ScanInfo_t     *si;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  if(! esl_opt_GetBoolean(go, "-g")) cm->config_opts  |= CM_CONFIG_LOCAL;
  cm->config_opts |= CM_CONFIG_QDB;
  ConfigCM(cm, NULL, NULL);

  if (esl_opt_GetBoolean(go, "--noqdb")) { 
    dmin = NULL; dmax = NULL;
  }
  else { dmin = cm->dmin; dmax = cm->dmax; }

  si = cm_CreateScanInfo(cm); /* impt to do this after QDBs set up in ConfigCM() */

  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, cm->null, abc->K, L, dsq);

      esl_stopwatch_Start(w);
      sc = FastCYKScan(cm, si, dsq, 1, L, cm->W, 0., NULL, NULL);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      if (esl_opt_GetBoolean(go, "-x")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = XFastCYKScan(cm, dsq, dmin, dmax, 1, L, cm->W, 0., NULL, NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "XFastCYKScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "-w")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = RefCYKScan(cm, si, dsq, 1, L, cm->W, 0., NULL, NULL);
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
	  sc = FastIInsideScan(cm, si, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "FastIInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--riins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  sc = RefIInsideScan(cm, si, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "RefIInsideScan(): ", sc);
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
	  sc = FastFInsideScan(cm, si, dsq, 1, L, cm->W, 0., NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "FastFInsideScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "--rfins")) 
	{ 
	  cm->search_opts  |= CM_SEARCH_INSIDE;
	  esl_stopwatch_Start(w);
	  sc = RefFInsideScan(cm, si, dsq, 1, L, cm->W, 0., NULL, NULL);
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
    }
  cm_FreeScanInfo(cm, si);
  FreeCM(cm);
  free(dsq);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IMPL_FASTSEARCH_BENCHMARK*/


#if 0
/* Function: FastFInsideScan()
 * Date:     EPN, Thu Nov  1 18:48:00 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded Inside algorithm. If bands are NULL, reverts to non-banded.
 *
 * Args:     cm              - the covariance model
 *           dsq             - the digitized sequence
 *           dmin            - minimum bound on d for state v; 0..M
 *           dmax            - maximum bound on d for state v; 0..M          
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *           ret_best_hit_sc - RETURN score of best hit (reported to results) NULL if not-wanted
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
FastFInsideScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float cutoff, 
		search_results_t *results, float **ret_vsc, float *ret_best_hit_sc)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [j][v][d] */
  float  ***alpha_begl; 
  float    *vsc;                /* best score for each state (float) */
  float     vsc_root;           /* best overall score (score at ROOT_S) */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       ip_g;               /* offset i for gamma (i-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  float     best_hit_sc;        /* best hit score found */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int     **dnAA, **dxAA;       /* [1..j..W][0..v..M-1] min,max d value allowed for posn j, state v */
  int      *dnA,   *dxA;        /* tmp ptr to 1 row of dnAA, dxAA */

  int cnum;
  int dn;
  int dx;
  float *el_scA;
  int *jp_wA;
  float **init_scAA;
  int  ctr = 0;
  /*int yidx;*/
  /*float const *tsc = cm->tsc[0]; */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in FastCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in FastCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in FastCYKScan, CM_SEARCH_INSIDE flag raised");

  /* determine if we're doing banded/non-banded */
  if(dmin != NULL && dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 

  vsc = NULL;
  if(ret_vsc != NULL) { 
    ESL_ALLOC(vsc, sizeof(float) * cm->M);
    esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);
  }
  best_hit_sc = IMPOSSIBLE;
  vsc_root    = IMPOSSIBLE;

  /*
   * alpha and alpha_begl allocations.
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
   * Note if v is NOT a BEGL_S alpha[j][v] == NULL
   */

  /* allocate alpha */
  ESL_ALLOC(alpha, (sizeof(float **) * 2));
  ESL_ALLOC(alpha[0], sizeof(float *) * cm->M);
  ESL_ALLOC(alpha[1], sizeof(float *) * cm->M);
  ESL_ALLOC(alpha[0][0], (sizeof(float) * 2 * (cm->M) * (W+1)));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      alpha[0][v] = alpha[0][0] + (v           * (W+1));
      alpha[1][v] = alpha[0][0] + ((v + cm->M) * (W+1));
    }
    else { /* BEGL_S, this is wasteful */
      alpha[0][v] = NULL;
      alpha[1][v] = NULL;
    }
  }
  /* float const *alphap = alpha[0][0]; */

  /* allocate alpha_begl */
  ESL_ALLOC(alpha_begl, (sizeof(float **) * (W+1)));
  for (j = 0; j <= W; j++) {
    ESL_ALLOC(alpha_begl[j], (sizeof(float *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	ESL_ALLOC(alpha_begl[j][v], (sizeof(float) * (W+1)));
      }
      else /* non-BEGL */
	alpha_begl[j][v] = NULL;
    }
  }
  ESL_ALLOC(bestr, (sizeof(int) * (W+1)));

  /*
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   */
  /* initialize alpha and alpha_begl */
  for(v = cm->M-1; v >= 0; v--) 
    {
      if(cm->stid[v] != BEGL_S) 
	{
	  alpha[0][v][0] = IMPOSSIBLE;
	  if      (cm->sttype[v] == E_st)  { 
	    alpha[0][v][0] = alpha[1][v][0] = 0.;
	    /* rest of E deck is IMPOSSIBLE, this rewritten if QDB is on, (slightly wasteful). */
	    for (d = 1; d <= W; d++) alpha[0][v][d] = alpha[1][v][d] = IMPOSSIBLE;
	  }
	  else if (cm->sttype[v] == MP_st) alpha[0][v][1] = alpha[1][v][1] = IMPOSSIBLE;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[0][v][0] = cm->endsc[v];
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[0][v][0] = FLogsum(alpha[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	      alpha[0][v][0] = ESL_MAX(alpha[0][v][0], IMPOSSIBLE);
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* BEGL_S, left child state */
	      y = cm->cnum[v];
	      alpha[0][v][0] = alpha_begl[0][w][0] + alpha[0][y][0]; 
	    }

	  alpha[1][v][0] = alpha[0][v][0];
	}
      else /* v == BEGL_S */
	{
	  y = cm->cfirst[v];
	  alpha_begl[0][v][0] = cm->endsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    alpha_begl[0][v][0] = FLogsum(alpha_begl[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in alpha */
	  alpha_begl[0][v][0] = ESL_MAX(alpha_begl[0][v][0], IMPOSSIBLE);
	  for (j = 1; j <= W; j++) 
	    alpha_begl[j][v][0] = alpha_begl[0][v][0];
	}
    }
  bestr[0] = -1;

  /*
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   */
  if(results != NULL) { 
    ESL_ALLOC(gamma,  sizeof(float) * (L+1));
    gamma[0] = 0;
    ESL_ALLOC(gback,  sizeof(int)   * (L+1));
    gback[0] = -1;
    ESL_ALLOC(savesc, sizeof(float) * (L+1));
    ESL_ALLOC(saver,  sizeof(int)   * (L+1));
  }
  /*
   * query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) {
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < dmin[v] && d <=W; d++) 
	  for(j = 0; j < 2; j++)
	    alpha[j][v][d] = IMPOSSIBLE;
	for (d = dmax[v]+1; d <= W;      d++) 
	  for(j = 0; j < 2; j++)
	    alpha[j][v][d] = IMPOSSIBLE;
      }
      else
	{
	  for (d = 0; d < dmin[v] && d <=W; d++) 
	    for(j = 0; j <= W; j++)
	      alpha_begl[j][v][d] = IMPOSSIBLE;
	  for (d = dmax[v]+1; d <= W;      d++) 
	    for(j = 0; j <= W; j++)
	      alpha_begl[j][v][d] = IMPOSSIBLE;
	}
    }
  }

  /* precalculate minimum and maximum d for each state and each sequence index (1..j..W). 
   * this is not always just dmin, dmax, (for ex. if j < W).
   */
  ESL_ALLOC(dnAA, sizeof(int *) * (W+1));
  ESL_ALLOC(dxAA, sizeof(int *) * (W+1));
  
  dnAA[0] = dxAA[0] = NULL; /* corresponds to j == 0, which is out of bounds */
  for(j = 1; j <= W; j++) {
    ESL_ALLOC(dnAA[j], sizeof(int) * cm->M);
    ESL_ALLOC(dxAA[j], sizeof(int) * cm->M);

    for(v = 0; v < cm->M; v++) {
      if(do_banded) { 
	dnAA[j][v] = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	dxAA[j][v] = ESL_MIN(j, dmax[v]); 
	dxAA[j][v] = ESL_MIN(dxAA[j][v], W);
      }
      else { 
	dnAA[j][v] = (cm->sttype[v] == MP_st) ? 2 : 1;
	dxAA[j][v] = ESL_MIN(j, W); 
      }
    }
  }

  /* precalculate possible emission scores for each state */
  float **esc_vAA;
  int a,b;
  ESL_ALLOC(esc_vAA, sizeof(float *) * (cm->M));
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      ESL_ALLOC(esc_vAA[v], sizeof(float) * cm->abc->Kp);
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v], cm->abc->Kp, IMPOSSIBLE);
      for(a = 0; a < cm->abc->K; a++)
	esc_vAA[v][a] = cm->esc[v][a];
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++)
	esc_vAA[v][a] = esl_abc_FAvgScore(cm->abc, a, cm->esc[v]);
      break;
    case MP_st:
      ESL_ALLOC(esc_vAA[v], sizeof(float) * (cm->abc->Kp * cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v], (cm->abc->Kp * cm->abc->Kp), IMPOSSIBLE);
      for(a = 0; a < (cm->abc->Kp-1); a++)
	for(b = 0; b < (cm->abc->Kp-1); b++)
	  if(a < cm->abc->K && b < cm->abc->K)
	    esc_vAA[v][a * cm->abc->Kp + b] = cm->esc[v][(a * cm->abc->K) + b];
	  else
	    esc_vAA[v][a * cm->abc->Kp + b] = DegeneratePairScore(cm->abc, cm->esc[v], a, b);
      break;
    default:
      esc_vAA[v] = NULL;
      break;
    }
  }

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* precalculate the initial score for all alpha[v][j][d] cells, it's independent
   * of j, so we do it here, outside the for(j...) loop */
  ESL_ALLOC(init_scAA, sizeof(float *) * (cm->M));
  for (v = 0; v < cm->M; v++) 
    {
      ESL_ALLOC(init_scAA[v], sizeof(float) * (W+1));

      if(NOT_IMPOSSIBLE(cm->endsc[v]))
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = el_scA[d] + cm->endsc[v];
      else
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = IMPOSSIBLE;
    }

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Precalculate the 'emit mode' of each state to speed up the addition of emission 
   * scores, all states are either EMITLEFT, EMITRIGHT, EMITPAIR, or EMITNONE, this
   * collapses ILs and MLs into 1 value (for example) for the switch() statement inside the for(d...) loop
   * in the heart of the recursion which saves us time.
   */
  int *emitmodeA;
  ESL_ALLOC(emitmodeA, sizeof(int) * cm->M);
  for(v = 0; v < cm->M; v++) {
    switch (cm->sttype[v]) {
    case IL_st:
    case ML_st:
      emitmodeA[v] = EMITLEFT;
      break;
    case IR_st:
    case MR_st:
      emitmodeA[v] = EMITRIGHT;
      break;		
    case MP_st:
      emitmodeA[v] = EMITPAIR;
      break;		
    default:
      emitmodeA[v] = EMITNONE;
      break;
    }
  }

  /* Initialize sc_v to size of M */
  float *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else {          dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;
	  /* float const *esc_v = cm->esc[v]; */
	  float const *esc_v = esc_vAA[v]; 
	  float const *tsc_v = cm->tsc[v];
	  //float sc;
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
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	      }
	      else { kmin = 0; kmax = d; }

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = FLogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
#pragma ivdep
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = FLogsum(sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
		sc = FLogsum(arow2[dp_y] + tsc_v[2],
			     arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow5[dp_y] + tsc_v[5],
			      init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow4[dp_y] + tsc_v[4]);		
		sc = FLogsum(sc, arow3[dp_y] + tsc_v[3]);		
		sc = FLogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow3[dp_y] + tsc_v[3],
			     arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		
		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow4[dp_y] + tsc_v[4],
			     arow3[dp_y] + tsc_v[3]);		
		sc = FLogsum(sc, init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow1[dp_y] + tsc_v[1]);		
		sc = FLogsum(sc, arow2[dp_y] + tsc_v[2]);		
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		

		/* add in emission score, if any */
		switch (emitmodeA[v]) {
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc = FLogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc = FLogsum(sc, arow0[dp_y] + tsc_v[0]);		
		switch (emitmodeA[v]) {
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
	     * nested emitmodeA[v] switch is based on empirical
	     * frequency in large test set, more frequent guys come
	     * earlier, so average num calcs in each switch is
	     * minimized.
	     */

	    switch (cnum) {
	    case 3: 
	      arow0 = (float * const) alpha[jp_y][y];
	      arow1 = (float * const) alpha[jp_y][y+1];
	      arow2 = (float * const) alpha[jp_y][y+2];
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		/* ctr++; */
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep
	      for (d = dn; d <= dx; d++, dp_y++) {
		sc_v[d] = FLogsum(arow1[dp_y] + tsc_v[1],
			     init_scAA[v][dp_y]);
		sc_v[d] = FLogsum(sc_v[d], arow0[dp_y] + tsc_v[0]);		
	      }
	      break; 
	    } /* end of switch(cnum) */
	    /* add in emission score (if any), and set alpha[jp_v][v][d] cell */
	    switch (emitmodeA[v]) {
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
	    } /* end of switch (emitmodeA[v]) */
	  } /* end of else (cm->sttype[v] != B_st && cm->stid[v] !=  BEGL_S st && cm->sttype[v] != IL_st && cm->sttype[v] != IR_st) */
	  if(vsc != NULL) 
	    if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	    else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
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

      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = dnA[0]; d <= dxA[0]; d++) {
	bestr[d] = 0;	/* root of the traceback = root state 0 */
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + cm->tsc[0][0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = FLogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + cm->tsc[0][yoffset]));
      }
	
      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
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
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, alpha[jp_v][0][d]);
      if(results != NULL) 
	{
	  float sc;
	  /* get information on hits */
	  if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	    {
	      /* The little semi-Markov model that deals with multihit parsing:
	       */
	      gamma[jp_g]  = gamma[jp_g-1] + 0; 
	      gback[jp_g]  = -1;
	      savesc[jp_g] = IMPOSSIBLE;
	      saver[jp_g]  = -1;
	      for (d = dnA[0]; d <= dxA[0]; d++) {
		i    = j-d+1;
		ip_g = i-i0+1;
		sc = gamma[ip_g-1] + alpha[jp_v][0][d] + cm->sc_boost; 
		/* sc_boost is experimental technique for finding hits < 0 bits. value is 0.0 if technique not used. */
		if (sc > gamma[jp_g]) {
		  gamma[jp_g]  = sc;
		  gback[jp_g]  = i;
		  savesc[jp_g] = alpha[jp_v][0][d]; 
		  saver[jp_g]  = bestr[d];
		}
	      }
	    }
	  else {
	    /* Resolving overlaps greedily (RSEARCH style),  
	     * At least one hit is sent back for each j here.
	     * However, some hits can already be removed for the greedy overlap
	     * resolution algorithm.  Specifically, at the given j, any hit with a
	     * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
	    /* First, report hit with d of 1 if > cutoff */
	    if (alpha[jp_v][0][1] >= cutoff) 
	      if(results != NULL) 
		report_hit (j, j, bestr[1], alpha[jp_v][0][1], results);
	    bestd = 1;
	    if (alpha[jp_v][0][1] > best_hit_sc)
	      best_hit_sc = alpha[jp_v][0][1];
	    
	    /* Now, if current score is greater than maximum seen previous, report
	       it if >= cutoff and set new max */
	    for (d = 2; d <= W && d <= jp_g; d++) 
	      {
		if (alpha[jp_v][0][d] > best_hit_sc) best_hit_sc = alpha[jp_v][0][d];
		if (alpha[jp_v][0][d] > alpha[jp_v][0][bestd]) {
		  if (alpha[jp_v][0][d] >= cutoff)
		    if(results != NULL) 
		      report_hit (j-d+1, j, bestr[d], alpha[jp_v][0][d], results);
		  bestd = d;
		}
	      }
	  }
	}
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;
  
  /* free alpha and alpha_begl, everything we need is in gamma.
   */

  free(alpha[0][0]);
  free(alpha[1]);
  free(alpha[0]);
  free(alpha);
  for (j = 0; j <= W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S)
	free(alpha_begl[j][v]);
    free(alpha_begl[j]);
  }
  free(alpha_begl);
  free(bestr);

  if(results != NULL) 
    {
      /*
       * Traceback stage.
       * Recover all hits: an (i,j,sc) triple for each one.
       */
      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	{
	  j     = j0;
	  while (j >= i0) 
	    {
	      jp_g = j-i0+1;
	      if (gback[jp_g] == -1) /* no hit */
		j--; 
	      else                /* a hit, a palpable hit */
		{
		  if(savesc[jp_g] > best_hit_sc) best_hit_sc = savesc[jp_g];
		  if(savesc[jp_g] >= cutoff && results != NULL) /* report the hit */
		    report_hit(gback[jp_g], j, saver[jp_g], savesc[jp_g], results);
		  j = gback[jp_g]-1;
		}
	    }
	}
      free(gback);
      free(gamma);
      free(savesc);
      free(saver);
    }

  for(v = 0; v < cm->M; v++) {
    free(init_scAA[v]);
    if(esc_vAA[v] != NULL) free(esc_vAA[v]);
  }
  free(init_scAA);
  free(esc_vAA);

  for(j = 1; j <= W; j++) {
    free(dnAA[j]);
    free(dxAA[j]);
  }
  free(dnAA);
  free(dxAA);
  free(jp_wA);
  free(el_scA);
  free(emitmodeA);
  if (ret_best_hit_sc != NULL) *ret_best_hit_sc = best_hit_sc;
  if (ret_vsc         != NULL) *ret_vsc         = vsc;
  
  ESL_DPRINTF1(("FastFInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

/* Function: RefFInsideScan()
 * Date:     EPN, Sun Nov  4 16:02:17 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded Inside algorithm. If bands are NULL, reverts to non-banded.
 *
 * Args:     cm              - the covariance model
 *           dsq             - the digitized sequence
 *           dmin            - minimum bound on d for state v; 0..M
 *           dmax            - maximum bound on d for state v; 0..M          
 *           i0              - start of target subsequence (1 for full seq)
 *           j0              - end of target subsequence (L for full seq)
 *           W               - max d: max size of a hit
 *           cutoff          - minimum score to report
 *           results         - search_results_t to add to; if NULL, don't add to it
 *           ret_vsc         - RETURN: [0..v..M-1] best score at each state v, NULL if not-wanted
 *           ret_best_hit_sc - RETURN score of best hit (reported to results) NULL if not-wanted
 *
 * Returns:  Score of best overall hit (vsc[0]). Information on hits added to <results>.
 *           <ret_vsc> is filled with an array of the best hit to each state v (if non-NULL).
 *           Dies immediately if some error occurs.
 */
float 
RefFInsideScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float cutoff, 
	       search_results_t *results, float **ret_vsc, float *ret_best_hit_sc)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [j][v][d] */
  float  ***alpha_begl; 
  float    *vsc;                /* best score for each state (float) */
  float     vsc_root;           /* best overall score (score at ROOT_S) */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  int       v, w, y;            /* state indices */
  int       jp_v;  	        /* offset j for state v */
  int       jp_y;  	        /* offset j for state y */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       ip_g;               /* offset i for gamma (i-i0+1) */
  int       dp_y;               /* offset d for state y */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  float     best_hit_sc;        /* best hit score found */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  int     **dnAA, **dxAA;       /* [1..j..W][0..v..M-1] min,max d value allowed for posn j, state v */
  int      *dnA,   *dxA;        /* tmp ptr to 1 row of dnAA, dxAA */

  int cnum;
  int dn;
  int dx;
  float *el_scA;
  int *jp_wA;
  float **init_scAA;
  int  ctr = 0;
  /*int yidx;*/
  /*float const *tsc = cm->tsc[0]; */

  /* Contract check */
  if(j0 < i0)     cm_Fail("ERROR in FastCYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL) cm_Fail("ERROR in FastCYKScan, dsq is NULL\n");
  if(cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("ERROR in FastCYKScan, CM_SEARCH_INSIDE flag raised");

  /* determine if we're doing banded/non-banded */
  if(dmin != NULL && dmax != NULL) do_banded = TRUE;

  L = j0-i0+1;
  if (W > L) W = L; 

  vsc = NULL;
  if(ret_vsc != NULL) { 
    ESL_ALLOC(vsc, sizeof(float) * cm->M);
    esl_vec_FSet(vsc, cm->M, IMPOSSIBLE);
  }
  best_hit_sc = IMPOSSIBLE;
  vsc_root    = IMPOSSIBLE;

  /*
   * alpha and alpha_begl allocations.
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
   * Note if v is NOT a BEGL_S alpha[j][v] == NULL
   */

  /* allocate alpha */
  ESL_ALLOC(alpha, (sizeof(float **) * 2));
  ESL_ALLOC(alpha[0], sizeof(float *) * cm->M);
  ESL_ALLOC(alpha[1], sizeof(float *) * cm->M);
  ESL_ALLOC(alpha[0][0], (sizeof(float) * 2 * (cm->M) * (W+1)));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      alpha[0][v] = alpha[0][0] + (v           * (W+1));
      alpha[1][v] = alpha[0][0] + ((v + cm->M) * (W+1));
    }
    else { /* BEGL_S, this is wasteful */
      alpha[0][v] = NULL;
      alpha[1][v] = NULL;
    }
  }
  /* float const *alphap = alpha[0][0]; */

  /* allocate alpha_begl */
  ESL_ALLOC(alpha_begl, (sizeof(float **) * (W+1)));
  for (j = 0; j <= W; j++) {
    ESL_ALLOC(alpha_begl[j], (sizeof(float *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	ESL_ALLOC(alpha_begl[j][v], (sizeof(float) * (W+1)));
      }
      else /* non-BEGL */
	alpha_begl[j][v] = NULL;
    }
  }
  ESL_ALLOC(bestr, (sizeof(int) * (W+1)));

  /*
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
   */
  /* initialize alpha and alpha_begl */
  for(v = cm->M-1; v >= 0; v--) 
    {
      if(cm->stid[v] != BEGL_S) 
	{
	  alpha[0][v][0] = IMPOSSIBLE;
	  if      (cm->sttype[v] == E_st)  { 
	    alpha[0][v][0] = alpha[1][v][0] = 0.;
	    /* rest of E deck is IMPOSSIBLE, this rewritten if QDB is on, (slightly wasteful). */
	    for (d = 1; d <= W; d++) alpha[0][v][d] = alpha[1][v][d] = IMPOSSIBLE;
	  }
	  else if (cm->sttype[v] == MP_st) alpha[0][v][1] = alpha[1][v][1] = IMPOSSIBLE;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[0][v][0] = cm->endsc[v];
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[0][v][0] = FLogsum(alpha[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	      alpha[0][v][0] = FLogsum(alpha[0][v][0], IMPOSSIBLE);
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* BEGL_S, left child state */
	      y = cm->cnum[v];
	      alpha[0][v][0] = alpha_begl[0][w][0] + alpha[0][y][0]; 
	    }

	  alpha[1][v][0] = alpha[0][v][0];
	}
      else /* v == BEGL_S */
	{
	  y = cm->cfirst[v];
	  alpha_begl[0][v][0] = cm->endsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    alpha_begl[0][v][0] = FLogsum(alpha_begl[0][v][0], (alpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in alpha */
	  alpha_begl[0][v][0] = FLogsum(alpha_begl[0][v][0], IMPOSSIBLE);
	  for (j = 1; j <= W; j++) 
	    alpha_begl[j][v][0] = alpha_begl[0][v][0];
	}
    }
      bestr[0] = -1;

  /*
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   */
  if(results != NULL) { 
    ESL_ALLOC(gamma,  sizeof(float) * (L+1));
    gamma[0] = 0;
    ESL_ALLOC(gback,  sizeof(int)   * (L+1));
    gback[0] = -1;
    ESL_ALLOC(savesc, sizeof(float) * (L+1));
    ESL_ALLOC(saver,  sizeof(int)   * (L+1));
  }
  /*
   * query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) {
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < dmin[v] && d <=W; d++) 
	  for(j = 0; j < 2; j++)
	    alpha[j][v][d] = IMPOSSIBLE;
	for (d = dmax[v]+1; d <= W;      d++) 
	  for(j = 0; j < 2; j++)
	    alpha[j][v][d] = IMPOSSIBLE;
      }
      else
	{
	  for (d = 0; d < dmin[v] && d <=W; d++) 
	    for(j = 0; j <= W; j++)
	      alpha_begl[j][v][d] = IMPOSSIBLE;
	  for (d = dmax[v]+1; d <= W;      d++) 
	    for(j = 0; j <= W; j++)
	      alpha_begl[j][v][d] = IMPOSSIBLE;
	}
    }
  }

  /* precalculate minimum and maximum d for each state and each sequence index (1..j..W). 
   * this is not always just dmin, dmax, (for ex. if j < W).
   */
  ESL_ALLOC(dnAA, sizeof(int *) * (W+1));
  ESL_ALLOC(dxAA, sizeof(int *) * (W+1));
  
  dnAA[0] = dxAA[0] = NULL; /* corresponds to j == 0, which is out of bounds */
  for(j = 1; j <= W; j++) {
    ESL_ALLOC(dnAA[j], sizeof(int) * cm->M);
    ESL_ALLOC(dxAA[j], sizeof(int) * cm->M);

    for(v = 0; v < cm->M; v++) {
      if(do_banded) { 
	dnAA[j][v] = (cm->sttype[v] == MP_st) ? ESL_MAX(dmin[v], 2) : ESL_MAX(dmin[v], 1); 
	dxAA[j][v] = ESL_MIN(j, dmax[v]); 
	dxAA[j][v] = ESL_MIN(dxAA[j][v], W);
      }
      else { 
	dnAA[j][v] = (cm->sttype[v] == MP_st) ? 2 : 1;
	dxAA[j][v] = ESL_MIN(j, W); 
      }
    }
  }

  /* precalculate possible emission scores for each state */
  float **esc_vAA;
  int a,b;
  ESL_ALLOC(esc_vAA, sizeof(float *) * (cm->M));
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      ESL_ALLOC(esc_vAA[v], sizeof(float) * cm->abc->Kp);
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v], cm->abc->Kp, IMPOSSIBLE);
      for(a = 0; a < cm->abc->K; a++)
	esc_vAA[v][a] = cm->esc[v][a];
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++)
	esc_vAA[v][a] = esl_abc_FAvgScore(cm->abc, a, cm->esc[v]);
      break;
    case MP_st:
      ESL_ALLOC(esc_vAA[v], sizeof(float) * (cm->abc->Kp * cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v], (cm->abc->Kp * cm->abc->Kp), IMPOSSIBLE);
      for(a = 0; a < (cm->abc->Kp-1); a++)
	for(b = 0; b < (cm->abc->Kp-1); b++)
	  if(a < cm->abc->K && b < cm->abc->K)
	    esc_vAA[v][a * cm->abc->Kp + b] = cm->esc[v][(a * cm->abc->K) + b];
	  else
	    esc_vAA[v][a * cm->abc->Kp + b] = DegeneratePairScore(cm->abc, cm->esc[v], a, b);
      break;
    default:
      esc_vAA[v] = NULL;
      break;
    }
  }

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (W+1));
  for(d = 0; d <= W; d++) el_scA[d] = cm->el_selfsc * d;

  /* precalculate the initial score for all alpha[v][j][d] cells, it's independent
   * of j, so we do it here, outside the for(j...) loop */
  ESL_ALLOC(init_scAA, sizeof(float *) * (cm->M));
  for (v = 0; v < cm->M; v++) 
    {
      ESL_ALLOC(init_scAA[v], sizeof(float) * (W+1));

      if(NOT_IMPOSSIBLE(cm->endsc[v]))
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = el_scA[d] + cm->endsc[v];
      else
	for(d = 0; d <= W; d++)
	  init_scAA[v][d] = IMPOSSIBLE;
    }

  /* allocate array for precalc'ed rolling ptrs into BEGL deck, filled inside 'for(j...' loop */
  ESL_ALLOC(jp_wA, sizeof(float) * (W+1));

  /* Precalculate the 'emit mode' of each state to speed up the addition of emission 
   * scores, all states are either EMITLEFT, EMITRIGHT, EMITPAIR, or EMITNONE, this
   * collapses ILs and MLs into 1 value (for example) for the switch() statement inside the for(d...) loop
   * in the heart of the recursion which saves us time.
   */
  int *emitmodeA;
  ESL_ALLOC(emitmodeA, sizeof(int) * cm->M);
  for(v = 0; v < cm->M; v++) {
    switch (cm->sttype[v]) {
    case IL_st:
    case ML_st:
      emitmodeA[v] = EMITLEFT;
      break;
    case IR_st:
    case MR_st:
      emitmodeA[v] = EMITRIGHT;
      break;		
    case MP_st:
      emitmodeA[v] = EMITPAIR;
      break;		
    default:
      emitmodeA[v] = EMITNONE;
      break;
    }
  }

  /* Initialize sc_v to size of M */
  float *sc_v;
  ESL_ALLOC(sc_v, (sizeof(float) * (W+1)));
  esl_vec_FSet(sc_v, (W+1), IMPOSSIBLE);

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      float sc;
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else {          dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }
      /* precalcuate all possible rolling ptrs into the BEGL deck, so we don't wastefully recalc them inside inner DP loop */
      for(d = 0; d <= W; d++) jp_wA[d] = (j-d)%(W+1);

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */
	  if(cm->sttype[v] == E_st) continue;

	  /* float const *esc_v = cm->esc[v]; */
	  float const *esc_v = esc_vAA[v]; 
	  float const *tsc_v = cm->tsc[v];
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

	      sc = init_scAA[v][d]; /* state delta is 0 for B_st */
	      for (k = kmin; k <= kmax; k++) 
		sc = FLogsum(sc, (alpha_begl[jp_wA[k]][w][d-k] + alpha[jp_y][y][k]));
	      alpha[jp_v][v][d] = sc;
	      /* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    y = cm->cfirst[v]; 
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = FLogsum(sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);
	      alpha_begl[jp_v][v][d] = sc;
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* ! B_st, ! BEGL_S st */
	    y = cm->cfirst[v]; 
	    i = j - dnA[v] + 1;
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      sc = init_scAA[v][d]; /* state delta is 0 for BEGL_S st */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		sc = FLogsum(sc, alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]);

	      switch (emitmodeA[v]) {
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
	      }
	    }
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

      /* determine min/max d we're allowing for the root state and this position j */
      jp_v = cur;
      for (d = dnA[0]; d <= dxA[0]; d++) {
	bestr[d] = 0;	/* root of the traceback = root state 0 */
	y = cm->cfirst[0];
	alpha[jp_v][0][d] = ESL_MAX(IMPOSSIBLE, alpha[cur][y][d] + cm->tsc[0][0]);
	for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++) 
	  alpha[jp_v][0][d] = FLogsum (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + cm->tsc[0][yoffset]));
      }
	
      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
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
      /* find the best score */
      for (d = dnA[0]; d <= dxA[0]; d++) 
	vsc_root = ESL_MAX(vsc_root, alpha[jp_v][0][d]);
      if(results != NULL) 
	{
	  float sc;
	  /* get information on hits */
	  if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	    {
	      /* The little semi-Markov model that deals with multihit parsing:
	       */
	      gamma[jp_g]  = gamma[jp_g-1] + 0; 
	      gback[jp_g]  = -1;
	      savesc[jp_g] = IMPOSSIBLE;
	      saver[jp_g]  = -1;
	      for (d = dnA[0]; d <= dxA[0]; d++) {
		i    = j-d+1;
		ip_g = i-i0+1;
		sc = gamma[ip_g-1] + alpha[jp_v][0][d] + cm->sc_boost; 
		/* sc_boost is experimental technique for finding hits < 0 bits. value is 0.0 if technique not used. */
		if (sc > gamma[jp_g]) {
		  gamma[jp_g]  = sc;
		  gback[jp_g]  = i;
		  savesc[jp_g] = alpha[jp_v][0][d]; 
		  saver[jp_g]  = bestr[d];
		}
	      }
	    }
	  else {
	    /* Resolving overlaps greedily (RSEARCH style),  
	     * At least one hit is sent back for each j here.
	     * However, some hits can already be removed for the greedy overlap
	     * resolution algorithm.  Specifically, at the given j, any hit with a
	     * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
	    /* First, report hit with d of 1 if > cutoff */
	    if (alpha[jp_v][0][1] >= cutoff) 
	      if(results != NULL) 
		report_hit (j, j, bestr[1], alpha[jp_v][0][1], results);
	    bestd = 1;
	    if (alpha[jp_v][0][1] > best_hit_sc)
	      best_hit_sc = alpha[jp_v][0][1];
	    
	    /* Now, if current score is greater than maximum seen previous, report
	       it if >= cutoff and set new max */
	    for (d = 2; d <= W && d <= jp_g; d++) 
	      {
		if (alpha[jp_v][0][d] > best_hit_sc) best_hit_sc = alpha[jp_v][0][d];
		if (alpha[jp_v][0][d] > alpha[jp_v][0][bestd]) {
		  if (alpha[jp_v][0][d] >= cutoff)
		    if(results != NULL) 
		      report_hit (j-d+1, j, bestr[d], alpha[jp_v][0][d], results);
		  bestd = d;
		}
	      }
	  }
	}
      if(vsc != NULL) 
	if(cm->stid[v] != BEGL_S) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
	else                      for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[jp_v][v][d]);
    } /* end loop over end positions j */
  if(vsc != NULL) vsc[0] = vsc_root;
  
  /* free alpha and alpha_begl, everything we need is in gamma.
   */

  free(alpha[0][0]);
  free(alpha[1]);
  free(alpha[0]);
  free(alpha);
  for (j = 0; j <= W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S)
	free(alpha_begl[j][v]);
    free(alpha_begl[j]);
  }
  free(alpha_begl);
  free(bestr);

  if(results != NULL) 
    {
      /*
       * Traceback stage.
       * Recover all hits: an (i,j,sc) triple for each one.
       */
      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	{
	  j     = j0;
	  while (j >= i0) 
	    {
	      jp_g = j-i0+1;
	      if (gback[jp_g] == -1) /* no hit */
		j--; 
	      else                /* a hit, a palpable hit */
		{
		  if(savesc[jp_g] > best_hit_sc) best_hit_sc = savesc[jp_g];
		  if(savesc[jp_g] >= cutoff && results != NULL) /* report the hit */
		    report_hit(gback[jp_g], j, saver[jp_g], savesc[jp_g], results);
		  j = gback[jp_g]-1;
		}
	    }
	}
      free(gback);
      free(gamma);
      free(savesc);
      free(saver);
    }

  for(v = 0; v < cm->M; v++) {
    free(init_scAA[v]);
    if(esc_vAA[v] != NULL) free(esc_vAA[v]);
  }
  free(init_scAA);
  free(esc_vAA);

  for(j = 1; j <= W; j++) {
    free(dnAA[j]);
    free(dxAA[j]);
  }
  free(dnAA);
  free(dxAA);
  free(jp_wA);
  free(el_scA);
  free(emitmodeA);
  if (ret_best_hit_sc != NULL) *ret_best_hit_sc = best_hit_sc;
  if (ret_vsc         != NULL) *ret_vsc         = vsc;
  
  ESL_DPRINTF1(("FastFInsideScan() return score: %10.4f\n", vsc_root)); 
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

#endif

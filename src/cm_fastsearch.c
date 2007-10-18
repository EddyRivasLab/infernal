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
FastCYKScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float cutoff, 
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
	  /* make this first if() statement */
	  else { /* if cm->sttype[v] != B_st || BEGL_S st */
	    y    = cm->cfirst[v];
	    dp_y = dn - sd; /* initial dp_y, we increment it at end of 'for(d = ...' loop */
	    i    = j-dn+1;  /* initial i,    we decrement it when we access it, inside each possible case of the switch (cnum) below */

	    float const *arow0;
	    float const *arow1;
	    float const *arow2;
	    float const *arow3;
	    float const *arow4;
	    float const *arow5;

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
		switch (cm->sttype[v]) {
		case IL_st:
		case ML_st:
		  sc += esc_v[dsq[i--]];
		  break;
		case IR_st:
		case MR_st:
		  sc += esc_j;
		  break;		
		case MP_st:
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
		switch (cm->sttype[v]) {
		case IL_st:
		case ML_st:
		  sc += esc_v[dsq[i--]];
		  break;
		/* can't be IR_st, IR's don't have 6 children */
		case MR_st:
		  sc += esc_j;
		  break;		
		case MP_st: 
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
		switch (cm->sttype[v]) {
		case IL_st:
		case ML_st:
		  sc += esc_v[dsq[i--]];
		  break;
		  /* can't be an IR */
		case MR_st:
		  sc += esc_j;
		  break;		
		case MP_st: /* OPTIMIZE */
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
		switch (cm->sttype[v]) {
		case IL_st:
		case ML_st:
		  sc += esc_v[dsq[i--]];
		  break;
		case IR_st:
		case MR_st:
		  sc += esc_j;
		  break;		
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
		switch (cm->sttype[v]) {
		case IL_st:
		case ML_st:
		  sc += esc_v[dsq[i--]];
		  break;
		case IR_st:
		case MR_st:
		  sc += esc_j;
		  break;		
		} /* end of switch (cm->sttype[v]) */
		alpha[jp_v][v][d] = sc;
	      } /* end of for(d = dn; d <= dx; d++) */
	      break;
	    } /* end of switch(cnum) */
	    /* for (d = dn; d <= dx; d++) 
	       printf("alpha[j:%d][v:%d][d:%d]: %10.4f\n", j, v, d, alpha[jp_v][v][d]); */
	  } /* end of else (v != B_st && v != BEGL_st) */
	  if(vsc != NULL) for (d = dn; d <= dx; d++) vsc[v] = ESL_MAX(vsc[v], alpha[jp_v][v][d]);
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
	  alpha[jp_v][0][d] = ESL_MAX (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + cm->tsc[0][yoffset]));
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
  
  if (ret_best_hit_sc != NULL) *ret_best_hit_sc = best_hit_sc;
  if (ret_vsc         != NULL) *ret_vsc         = vsc;
  
  /* printf("FastCYKScan() return score: %10.4f\n", vsc_root); 
     printf("ctr: %d\n", ctr);*/
  return vsc_root;
  
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* NEVERREACHED */
}

/* Function: OLDFastCYKScan()
 * Date:     EPN, Wed Sep 12 16:55:28 2007
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
OLDFastCYKScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, float cutoff, 
	       search_results_t *results, float **ret_vsc, float *ret_best_hit_sc)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [j][v][d] (for non-BEGL states) */
  float  ***alpha_begl;         /* CYK DP score matrix, [j][v][d] (for     BEGL states) */
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
  int       jp_w;  	        /* offset j for state w */
  int       jp_g;               /* offset j for gamma (j-i0+1) */
  int       ip_g;               /* offset i for gamma (i-i0+1) */
  int       kmin, kmax;         /* for B_st's, min/max value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       sd;                 /* StateDelta(cm->sttype[v]), # emissions from v */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */
  float     best_hit_sc;        /* best hit score found */
  int       do_banded = FALSE;  /* TRUE: use QDBs, FALSE: don't   */
  float     sc;                 /* temp score for reporting hits */
  int     **dnAA, **dxAA;       /* [1..j..W][0..v..M-1] min,max d value allowed for posn j, state v */
  int      *dnA,   *dxA;        /* tmp ptr to 1 row of dnAA, dxAA */

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
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   */
  /*ESL_ALLOC(alpha, (sizeof(float **) * cm->M));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] == BEGL_S) {
	ESL_ALLOC(alpha[v], (sizeof(float *) * (W+1)));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) alpha[v] = alpha[cm->M-1];
    else {
      ESL_ALLOC(alpha[v], sizeof(float *) * 2);
      for (j = 0; j < 2; j++) ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
    }
    }*/

  /* allocate alpha */
  ESL_ALLOC(alpha, (sizeof(float **) * 2));
  for (j = 0; j < 2; j++) {
    ESL_ALLOC(alpha[j], (sizeof(float *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] != BEGL_S) {
	ESL_ALLOC(alpha[j][v], (sizeof(float) * (W+1)));
      }
      else /* BEGL */
	alpha[j][v] = NULL;
    }
  }

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
	  if      (cm->sttype[v] == E_st)  alpha[0][v][0] = 0;
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
  /* ...we don't bother to look at local alignment starts here... */
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

  /* The main loop: scan the sequence from position i0 to j0.
   */
  for (j = i0; j <= j0; j++) 
    {
      jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
      cur  = j%2;
      prv  = (j-1)%2;
      if(jp_g >= W) { dnA = dnAA[W];     dxA = dxAA[W];    }
      else {          dnA = dnAA[jp_g];  dxA = dxAA[jp_g]; }

      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* printf("dnA[v:%d]: %d\ndxA[v:%d]: %d\n", v, dnA[v], v, dxA[v]); */

	  jp_v = (cm->stid[v] == BEGL_S) ? (j % (W+1)) : cur;
	  jp_y = (StateRightDelta(cm->sttype[v]) > 0) ? prv : cur;
	  sd   = StateDelta(cm->sttype[v]);
	  float const *esc_v = cm->esc[v];

	  if(cm->sttype[v] == B_st) {
	    w = cm->cfirst[v]; /* BEGL_S left child */
	    y = cm->cnum[v];
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      /* k is the length of the right fragment */
	      /* Careful, make sure k is consistent with bands in state w and state y. */
	      if(do_banded) {
		kmin = ESL_MAX(dmin[y], (d-dmax[w]));
		kmin = ESL_MAX(kmin, 0);
		kmax = ESL_MIN(dmax[y], (d-dmin[w]));
	      }
	      else { kmin = 0; kmax = d; }

	      alpha[jp_v][v][d] = ESL_MAX(IMPOSSIBLE, cm->endsc[v] + (cm->el_selfsc * (d - sd)));
	      for (k = kmin; k <= kmax; k++) { 
		jp_w = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		alpha[jp_v][v][d] = ESL_MAX(alpha[jp_v][v][d], (alpha_begl[jp_w][w][d-k] + alpha[jp_y][y][k]));
		/* careful: scores for w, the BEGL_S child of v, are in alpha_begl, not alpha */
	      }
	    }
	  }
	  else if (cm->stid[v] == BEGL_S) {
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      alpha_begl[jp_v][v][d] = ESL_MAX (IMPOSSIBLE, (cm->endsc[v] + (cm->el_selfsc * (d - sd))));
	      y = cm->cfirst[v]; 
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha_begl[jp_v][v][d] = ESL_MAX (alpha_begl[jp_v][v][d], (alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]));
	      /* careful: y is in alpha (all children of a BEGL_S must be non BEGL_S) */
	    }
	  }
	  else { /* if cm->sttype[v] != B_st && cm->stid[v] != BEGL_S */
	    for (d = dnA[v]; d <= dxA[v]; d++) {
	      alpha[jp_v][v][d] = ESL_MAX (IMPOSSIBLE, (cm->endsc[v] + (cm->el_selfsc * (d - sd))));
	      y = cm->cfirst[v]; 
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[jp_v][v][d] = ESL_MAX (alpha[jp_v][v][d], (alpha[jp_y][y+yoffset][d - sd] + cm->tsc[v][yoffset]));
		
	      /* add in emission score, if any */
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st: 
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[jp_v][v][d] += esc_v[(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[cur][v][d] += DegeneratePairScore(cm->abc, esc_v, dsq[i], dsq[j]);
		break;
	      case ML_st:
	      case IL_st:
		alpha[cur][v][d] += esl_abc_FAvgScore(cm->abc, dsq[i], esc_v);
		break;
	      case MR_st:
	      case IR_st:
		alpha[cur][v][d] += esl_abc_FAvgScore(cm->abc, dsq[j], esc_v);
		break;
	      } /* end of switch */
	    } /* end of d = dnA[v]; d <= dxA[v]; d++ for B_st */
	    /* for (d = dnA[v]; d <= dxA[v]; d++) 
	       printf("alpha[j:%d][v:%d][d:%d]: %10.4f\n", j, v, d, alpha[cur][v][d]); */
	  } /* end of else (v != B_st && v != BEGL_S */
	  if(vsc != NULL) { 
	    if(cm->stid[v] == BEGL_S) for (d = dnA[v]; d <= dxA[v]; d++) vsc[v] = ESL_MAX(vsc[v], alpha_begl[v][jp_v][d]); 
	    else                      for (d = dnA[v]; d <= dxA[v]; d++) vsc[v] = ESL_MAX(vsc[v], alpha[v][jp_v][d]); 
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
	  alpha[jp_v][0][d] = ESL_MAX (alpha[jp_v][0][d], (alpha[cur][y+yoffset][d] + cm->tsc[0][yoffset]));
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
  for (j = 0; j < 2; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] != BEGL_S)
	free(alpha[j][v]);
    free(alpha[j]);
  }
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

  for(j = 1; j <= W; j++) {
    free(dnAA[j]);
    free(dxAA[j]);
  }
  free(dnAA);
  free(dxAA);
  
  if (ret_best_hit_sc != NULL) *ret_best_hit_sc = best_hit_sc;
  if (ret_vsc         != NULL) *ret_vsc         = vsc;
  
  /* printf("OLDFastCYKScan() return score: %10.4f\n", vsc_root); */
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
#ifdef INTEL_COMPILER
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
#ifdef INTEL_COMPILER
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
#ifdef INTEL_COMPILER
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
#ifdef INTEL_COMPILER
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
#ifdef INTEL_COMPILER
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
#ifdef INTEL_COMPILER
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
#ifdef INTEL_COMPILER
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
#ifdef INTEL_COMPILER
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

/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_FASTSEARCH_BENCHMARK
/* gcc -o benchmark-fastsearch -g -O2 -I. -L. -I../easel -L../easel -DIMPL_FASTSEARCH_BENCHMARK cm_fastsearch.c -linfernal -leasel -lm
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
  { "-w",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute slow, reference CYK scan implementation", 0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute experimental CYK scan implementation", 0 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute non-banded optimized CYK scan implementation", 0 },
  { "--rsearch", eslARG_NONE,   FALSE, NULL, NULL,  NULL,"--noqdb", NULL, "also execute ported RSEARCH's CYK scan implementation", 0 },
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

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  cm->config_opts |= CM_CONFIG_QDB;
  ConfigCM(cm, NULL, NULL);

  if (esl_opt_GetBoolean(go, "--noqdb")) { 
    dmin = NULL; dmax = NULL;
  }
  else { dmin = cm->dmin; dmax = cm->dmax; }

  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, cm->null, abc->K, L, dsq);

      esl_stopwatch_Start(w);
      sc = FastCYKScan(cm, dsq, dmin, dmax, 1, L, cm->W, 0., NULL, NULL, NULL);
      printf("%4d %-30s %10.4f bits ", (i+1), "FastCYKScan(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");

      if (esl_opt_GetBoolean(go, "-x")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = OLDFastCYKScan(cm, dsq, dmin, dmax, 1, L, cm->W, 0., NULL, NULL, NULL);
	  printf("%4d %-30s %10.4f bits ", (i+1), "OLDFastCYKScan(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}

      if (esl_opt_GetBoolean(go, "-w")) 
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
    }

  FreeCM(cm);
  free(dsq);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IMPL_FASTSEARCH_BENCHMARK*/

/* fastsearch.c
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
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/*****************************************************************
 * 1. The CP9_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  cp9_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Wed Jun 27 08:07:01 2007 [Janelia]
 *
 * Purpose:   Create a profile of <M> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */
CP9_OPROFILE *
cp9_oprofile_Create(int M, const ESL_ALPHABET *abc)
{
  int          status;
  CP9_OPROFILE *om = NULL;
  int          x;

  /* level 0 */
  ESL_ALLOC(om, sizeof(CP9_OPROFILE));
  om->tsc = NULL;
  om->rsc = NULL;

  /* level 1 */
  ESL_ALLOC(om->tsc,   sizeof(int)   * (M+1)  * cp9O_NTRANS);
  ESL_ALLOC(om->rsc,   sizeof(int *) * abc->Kp);                     
  om->rsc[0] = NULL;

  ESL_ALLOC(om->rsc[0],sizeof(int)   * abc->Kp * (M+1) * cp9O_NR);
  for (x = 1; x < abc->Kp; x++)
    om->rsc[x] = om->rsc[0] + (x * (M+1) * cp9O_NR);

  om->M     = M;
  om->abc   = abc;

  /* NOT SURE WHAT TO DO HERE! Initializations of unused model transitions 
   * tsc[0][0..11] are unused: set to -infty. [7=TBM is used, though]. These values are touched by DP algs.
   */
  /*for (x = 0; x < cp9O_NTRANS; x++) om->tsc[x] = cp9_IMPOSSIBLE;*/
  return om;

 ERROR:
  cp9_oprofile_Destroy(om);
  return NULL;
}

/* Function:  cp9_oprofile_Destroy()
 * Synopsis:  Free an optimized profile structure.
 * Incept:    SRE, Wed Jun 27 08:07:08 2007 [Janelia]
 *
 * Purpose:   Frees profile <om>.
 */
void
cp9_oprofile_Destroy(CP9_OPROFILE *om)
{
  if (om != NULL) {
    if (om->rsc   != NULL && om->rsc[0] != NULL) free(om->rsc[0]);
    if (om->tsc   != NULL) free(om->tsc);
    if (om->rsc   != NULL) free(om->rsc);
    free(om);
  }
  return;
}


/*****************************************************************
 * 2. The CP9_OMX structure: a dynamic programming matrix.
 *****************************************************************/

/* Function:  cp9_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Wed Jun 27 08:07:29 2007 [Janelia]
 *
 * Purpose:   Allocate a reusable, resizeable <CP9_OMX> for models up to
 *            size <allocM> and sequences up to length <allocL>.
 *
 * Returns:   a pointer to the new <CP9_OMX>.
 *
 * Throws:    <NULL> on allocation error.
 */
CP9_OMX *
cp9_omx_Create(int allocM, int allocL)
{
  int      status;
  CP9_OMX  *ox;
  int      i;

  ESL_ALLOC(ox, sizeof(CP9_OMX));
  ox->dp  = NULL;

  ESL_ALLOC(ox->dp,  sizeof(int *) * (allocL+1));
  ox->dp[0] = NULL;
  
  ESL_ALLOC(ox->dp[0], sizeof(int) * (allocL+1) * (allocM+1) * cp9X_NSCELLS);
  for (i = 1; i <= allocL; i++)
    ox->dp[i] = ox->dp[0] + i*(allocM+1)*cp9X_NSCELLS;
  
  /* Initialize memory that's allocated but unused, only to keep
   * valgrind and friends happy.
   */
  /* NOT SURE WHAT TO DO HERE! 
     for (i = 0; i <= allocL; i++) 
     {*/
  //ox->dp[i][0      * cp9X_NSCELLS + cp9X_M] = cp9_IMPOSSIBLE; /* M_0 */
  //ox->dp[i][0      * cp9X_NSCELLS + cp9X_I] = cp9_IMPOSSIBLE; /* I_0 */      
  //ox->dp[i][0      * cp9X_NSCELLS + cp9X_D] = cp9_IMPOSSIBLE; /* D_0 */
  //ox->dp[i][1      * cp9X_NSCELLS + cp9X_D] = cp9_IMPOSSIBLE; /* D_1 */
  //ox->dp[i][allocM * cp9X_NSCELLS + cp9X_I] = cp9_IMPOSSIBLE; /* I_M */
  //}

  ox->nrows  = allocL+1;
  ox->ncells = (allocM+1)*(allocL+1);
  ox->M      = allocM;
  ox->L      = allocL;
  return ox;

 ERROR:
  cp9_omx_Destroy(ox);
  return NULL;
}

/* Function:  cp9_omx_Destroy()
 * Synopsis:  Frees an optimized DP matrix.
 * Incept:    SRE, Wed Jun 27 08:07:37 2007 [Janelia]
 *
 * Purpose:   Frees optimized DP matrix <ox>.
 *
 * Returns:   (void)
 */
void
cp9_omx_Destroy(CP9_OMX *ox)
{
  if (ox == NULL) return;
  if (ox->dp != NULL) {
    if (ox->dp[0] != NULL) free(ox->dp[0]);
    free(ox->dp);
  }
  free(ox);
}


/* Function:  cp9_omx_GrowTo()
 * Incept:    SRE, Wed Jul 18 09:12:34 2007 [Janelia]
 *
 * Purpose:   Assures that a DP matrix <ox> is allocated
 *            for a model of size up to <allocM> and a sequence of
 *            length up to <allocL>; reallocates if necessary.
 *            
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Returns:   <eslOK> on success, and <ox> may be reallocated upon
 *            return; any data that may have been in <ox> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <ox> must be assumed to be invalidated.
 */
int
cp9_omx_GrowTo(CP9_OMX *ox, int allocM, int allocL)
{
  int     status;
  void   *p;
  int     i;
  size_t  ncells, nrows;

  ncells = (allocM+1)*(allocL+1);
  nrows  = allocL+1;
  if (ncells <= ox->ncells && nrows <= ox->nrows) { ox->M = allocM; ox->L = allocL; return eslOK; }

  /* must we realloc the 2D matrices? (or can we get away with just
   * jiggering the row pointers, if we are growing in one dimension
   * while shrinking in another?)
   */
  if (ncells > ox->ncells) 
    {
      ESL_RALLOC(ox->dp[0], p, sizeof(int) * ncells * cp9X_NSCELLS);
      ox->ncells = ncells;
    }

  /* must we reallocate the row pointers? */
  if (nrows > ox->nrows)
    {
      ESL_RALLOC(ox->dp,  p, sizeof(int *) * nrows);
      ox->nrows  = nrows;
    }

  /* reset all the row pointers (even though we may not have to: if we
   * increased allocL without reallocating cells or changing M,
   * there's some wasted cycles here - but I'm not sure this case
   * arises, and I'm not sure it's worth thinking hard about. )
   */
  for (i = 0; i <= allocL; i++) 
    ox->dp[i] = ox->dp[0] + i * (allocM+1) * cp9X_NSCELLS;

  ox->M      = allocM;
  ox->L      = allocL;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  cp9_omx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 * Incept:    SRE, Wed Jul 18 08:35:17 2007 [Janelia]
 *
 * Purpose:   Dump matrix <ox> to stream <fp> for diagnostics.
 */
int
cp9_omx_Dump(FILE *ofp, CP9_OMX *ox)
{
  int i, k, x;

  /* Header */
  fprintf(ofp, "     ");
  for (k = 0; k <= ox->M;  k++) fprintf(ofp, "%8d ", k);
  fprintf(ofp, "%8s ", "E");
  fprintf(ofp, "%8s ", "N");
  fprintf(ofp, "%8s ", "J");
  fprintf(ofp, "%8s ", "B");
  fprintf(ofp, "%8s\n","C");
  for (k = 0; k <= ox->M+5; k++) fprintf(ofp, "%8s ", "--------");
  fprintf(ofp, "\n");
  
  /* DP matrix data */
  for (i = 0; i <= ox->L; i++)
    {
      fprintf(ofp, "%3d M ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8d ", ox->dp[i][k * cp9X_NSCELLS + cp9X_M]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d I ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8d ", ox->dp[i][k * cp9X_NSCELLS + cp9X_I]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d D ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8d ", ox->dp[i][k * cp9X_NSCELLS + cp9X_D]);
      fprintf(ofp, "\n\n");

      fprintf(ofp, "%3d EL ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8d ", ox->dp[i][k * cp9X_NSCELLS + cp9X_EL]);
      fprintf(ofp, "\n\n");
    }
  return eslOK;
}


/* Function:  cp9_oprofile_Convert()
 * Synopsis:  Convert generic profile to optimized profile.
 * Incept:    SRE, Wed Jun 27 08:09:02 2007 [Janelia]
 *
 * Purpose:   Use the <cm->cp9> scores to fill in the scores of <om>,
 *            rearranging as needed.
 * 
 * Returns:   <eslOK> on success.
 */
int
cp9_oprofile_Convert(CM_t *cm, CP9_OPROFILE *om)
{
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in cp9_oprofile_Convert, cm->cp9 is NULL.\n");

  int  k, x;
  
  /* Transition scores */
  for (k = 0 ; k <= cm->cp9->M; k++) {
    int *tsc = om->tsc + k*cp9O_NTRANS;
    
    tsc[cp9O_MM] = cm->cp9->tsc[CTMM][k];
    tsc[cp9O_MI] = cm->cp9->tsc[CTMI][k];
    tsc[cp9O_MD] = cm->cp9->tsc[CTMD][k];
    tsc[cp9O_IM] = cm->cp9->tsc[CTIM][k];
    tsc[cp9O_II] = cm->cp9->tsc[CTII][k];
    tsc[cp9O_DM] = cm->cp9->tsc[CTDM][k];
    tsc[cp9O_DD] = cm->cp9->tsc[CTDD][k];
    tsc[cp9O_ID] = cm->cp9->tsc[CTID][k];
    tsc[cp9O_DI] = cm->cp9->tsc[CTDI][k];
    tsc[cp9O_BM] = cm->cp9->bsc[k];
    tsc[cp9O_MEL]= cm->cp9->tsc[CTME][k];
    tsc[cp9O_ME] = cm->cp9->esc[k];
  }

  /* Match scores */
  for (x = 0; x < cm->abc->Kp; x++)
    for (k = 1; k <= cm->cp9->M; k++)
      {
	int *rsc = om->rsc[x] + k*cp9O_NR;
	rsc[cp9O_MSC] = cm->cp9->msc[x][k];
      }

  /* Insert scores */
  for (x = 0; x < cm->abc->Kp; x++)
    for (k = 0; k <= cm->cp9->M; k++)
      {
	int *rsc = om->rsc[x] + k*cp9O_NR;
	rsc[cp9O_ISC] = cm->cp9->isc[x][k];
      }
  om->abc   = cm->abc;
  return eslOK;
}

#define MMX(i,k)  (dp[(i)][(k) * cp9X_NSCELLS + cp9X_M])
#define IMX(i,k)  (dp[(i)][(k) * cp9X_NSCELLS + cp9X_I])
#define DMX(i,k)  (dp[(i)][(k) * cp9X_NSCELLS + cp9X_D])
#define ELMX(i,k) (dp[(i)][(k) * cp9X_NSCELLS + cp9X_EL])
   
#define TSC(s,k) (tsc[(k) * cp9O_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * cp9O_NR     + cp9O_MSC])
#define ISC(k)   (rsc[(k) * cp9O_NR     + cp9O_ISC])

/***********************************************************************
 * Function: cp9_FastViterbi()
 * 
 * Purpose:  Runs the Viterbi dynamic programming algorithm on an
 *           input subsequence (i0-j0). 
 *           Somewhat flexible based on input options.
 *    
 * Note:     IDENTICAL to cp9_FastForward() below with maxes replacing
 *           sums in DP recursion, and the possibility of returning a
 *           CP9 trace if doing_align == TRUE. See cp9_FastForward() for more info,
 *           including more verbose 'Purpose' and description of arguments.
 *
 * Returns:  if(!do_scan) log P(S,tr|M)/P(S,tr|R), as a bit score
 *           else         max log P(S,tr|M)/P(S,tr|R), for argmax subseq S of input seq i0..j0,
 */
float
cp9_FastViterbi(CM_t *cm, const CP9_OPROFILE *om, const CP9_OMX *ox, ESL_DSQ *dsq, 
		int i0, int j0, int W, float cutoff, int **ret_sc, 
		int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
		int be_efficient, CP9_dpmatrix_t **ret_mx, CP9trace_t **ret_tr)
{
  int          status;
  int          j;           /*     actual   position in the subsequence                     */
  int          jp;          /* j': relative position in the subsequence                     */
  int          i;           /* j-W: position in the subsequence                             */
  int          ip;          /* i': relative position in the subsequence                     */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  int          L;           /* j0-i0+1: subsequence length                                  */
  //CP9_dpmatrix_t *mx;       /* the CP9 DP matrix                                            */
  //int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  //int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  //  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  //int        **elmx;        /* DP matrix for EL state scores [0..1][0..cm->cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int         *scA;         /* prob (seq from j0..jp | HMM) [0..jp..cm->cp9->M]             */
  int          tmp_sc;      /* temporary int log odds score                                 */
  float        fsc;         /* float log odds score                                         */
  float        curr_sc;     /* temporary score used for filling in gamma                    */
  float       *gamma;       /* SHMM DP matrix for optimum nonoverlap resolution [0..j0-i0+1]*/
  int         *gback;       /* traceback pointers for SHMM                      [0..j0-i0+1]*/ 
  float       *savesc;      /* saves score of hit added to best parse at j      [0..j0-i0+1]*/ 
  float        best_hmm_sc; /* Best overall score from semi-HMM to return if do_scan        */
  float        best_hmm_pos;/* residue giving best_hmm_sc                                   */
  float        best_sc;     /* Best score overall, returned if 0 hits found by HMM & do_scan*/
  float        best_pos;    /* residue giving best_sc                                       */
  float        return_sc;   /* score to return, if (!do_scan) return overall Viterbi sc,    *
			     * else return best_hmm_sc if # HMM hits>0, else return best_sc */
  int          nrows;       /* num rows for DP matrix, 2 or L+1 depending on be_efficient   */
  int          c;           /* counter for EL states */
  CP9trace_t  *tr;

  int const *tsc = om->tsc;
  int      **dp  = ox->dp;
  int        M   = om->M;

  /*debug_print_cp9_params(stdout, cm->cp9, TRUE);*/

  /*printf("in cp9_FastViterbi() i0: %d j0: %d\n", i0, j0);  */
  /* Contract checks */
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in cp9_FastViterbi, cm->cp9 is NULL.\n");
  if(be_efficient && (ret_mx != NULL))
    esl_fatal("ERROR in cp9_FastViterbi, be_efficient is TRUE, but ret_mx is non-NULL\n");
  if(results != NULL && !do_scan)
    esl_fatal("ERROR in cp9_FastViterbi, passing in results data structure, but not in scanning mode.\n");
  if((cm->search_opts & CM_SEARCH_HMMGREEDY) && 
     (cm->search_opts & CM_SEARCH_HMMRESCAN))
    esl_fatal("ERROR in cp9_FastViterbi, CM_SEARCH_HMMGREEDY and CM_SEARCH_HMMRESCAN flags up, this combo not yet implemented. Implement it!\n");
  if(dsq == NULL)
    esl_fatal("ERROR in cp9_FastViterbi, dsq is NULL.");

  /* TEMPORARY */
  if(ret_tr != NULL) 
    cm_Fail("ERROR in cp9_FastViterbi, asking for ret_tr, but you haven't implemented a ViterbiTrace function that handles an optimized dp matrix yet.");

  best_sc     = IMPOSSIBLE;
  best_pos    = -1;
  best_hmm_sc = IMPOSSIBLE;
  best_hmm_pos= -1;
  L = j0-i0+1;
  if (W > L) W = L; 
  /* temporary ? */
  if(be_efficient) ESL_ALLOC(erow, sizeof(int) * 2);
  else             ESL_ALLOC(erow, sizeof(int) * (L+1));

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0.;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));

  /* Allocate DP matrix, either 2 rows or L+1 rows (depending on be_efficient),
   * M+1 columns */ 
  if(be_efficient) nrows = 2;
  else             nrows = L+1;
  /*mx = AllocCPlan9Matrix(nrows, M, &mmx, &imx, &dmx, &elmx, &erow); */

  /* scA will hold P(seq up to j | Model) in int log odds form */
  ESL_ALLOC(scA, sizeof(int) * (j0-i0+2));
			
  /* Initialization of the zero row. */
  MMX(0,0) = 0;      /* M_0 is state B, and everything starts in B */
  IMX(0,0) = -INFTY; /* I_0 is state N, can't get here without emitting*/
  DMX(0,0) = -INFTY; /* D_0 doesn't exist. */
  ELMX(0,0)= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY;   
  /*printf("mmx[jp:%d][%d]: %d %f\n", 0, 0, mmx[0][0], Score2Prob(mmx[0][0], 1.));
    printf("imx[jp:%d][%d]: %d %f\n", 0, 0, imx[0][0], Score2Prob(imx[0][0], 1.));
    printf("dmx[jp:%d][%d]: %d %f\n", 0, 0, dmx[0][0], Score2Prob(dmx[0][0], 1.));
    printf("elmx[jp:%d][%d]: %d %f\n", 0, 0, elmx[0][0], Score2Prob(elmx[0][0], 1.));*/

  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  for (k = 1; k <= M; k++)
    {
      MMX(0,k) = IMX(0,k) = ELMX(0,k) = -INFTY;      /* need seq to get here */
      DMX(0,k) = ESL_MAX(MMX(0,k-1) + TSC(cp9O_MD,k-1),
		   IMX(0,k-1) + TSC(cp9O_ID,k-1));
      DMX(0,k) = ESL_MAX(DMX(0,k), DMX(0,k-1) + TSC(cp9O_DD,k-1));
      /*printf("mmx[jp:%d][%d]: %d %f\n", 0, k, mmx[0][k], Score2Prob(mmx[0][k], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", 0, k, imx[0][k], Score2Prob(imx[0][k], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));
	printf("elmx[jp:%d][%d]: %d %f\n", 0, k, dmx[0][k], Score2Prob(dmx[0][k], 1.));*/
    }
  /* We can do a full parse through all delete states. */
  erow[0]  = DMX(0,M) + TSC(cp9O_DM,M); 
  scA[0]   = erow[0];
  fsc      = Scorify(scA[0]);
  /*printf("jp: %d j: %d fsc: %f isc: %d\n", jp, j, fsc, isc[jp]);*/
  if(fsc > best_sc) 
    {
      best_sc = fsc;
      best_pos= i0-1;
    }
  /*printf("jp: %d j: %d fsc: %f sc: %d\n", 0, i0-1, Scorify(sc[0]), sc[0]);*/

  /*****************************************************************
   * The main loop: scan the sequence from position i0 to j0.
   *****************************************************************/
  /* Recursion. */
  int ctr = 0;

  for (j = i0; j <= j0; j++)
    {
      int endsc = -INFTY;
      int const *rsc = om->rsc[dsq[j]];
      int el_selfsc = cm->cp9->el_selfsc;
      int sc;

      jp = j-i0+1;     /* jp is relative position in the sequence 1..L */
      cur = (j-i0+1); 
      prv = (j-i0);
      if(be_efficient) {
	cur %= 2;
	prv %= 2;
      }	

      /* The 1 difference between a Viterbi scanner and the 
       * regular Viterbi. In non-scanner parse must begin in B at
       * position 0 (i0-1), in scanner we can start at any position 
       * in the seq. */
      MMX(cur,0)  = (do_scan == TRUE) ? 0 : -INFTY;
      DMX(cur,0)  = -INFTY;  /*D_0 is non-existent*/
      ELMX(cur,0)  = -INFTY; /*no EL state for node 0 */

      sc = ESL_MAX(MMX(prv,0) + TSC(cp9O_MI,0),
		   IMX(prv,0) + TSC(cp9O_II,0));
      sc = ESL_MAX(sc, DMX(prv,0) + TSC(cp9O_DI,0));
      IMX(cur,0) = ESL_MAX(sc + ISC(0), -INFTY);

      /*printf("mmx[jp:%d][%d]: %d %f\n", jp, 0, mmx[cur][0], Score2Prob(mmx[cur][0], 1.));
	printf("imx[jp:%d][%d]: %d %f\n", jp, 0, imx[cur][0], Score2Prob(imx[cur][0], 1.));
	printf("dmx[jp:%d][%d]: %d %f\n", jp, 0, dmx[cur][0], Score2Prob(dmx[cur][0], 1.));*/

      for (k = 1; k <= M; k++)
	{
	  /*match state*/
	  /* ctr++; */
	  sc = ESL_MAX(    MMX(prv,k-1)  + TSC(cp9O_MM,k-1),
			   IMX(prv,k-1)  + TSC(cp9O_IM,k-1));
	  sc = ESL_MAX(sc, DMX(prv,k-1)  + TSC(cp9O_DM,k-1));
	  sc = ESL_MAX(sc, MMX(prv,0)    + TSC(cp9O_BM,k));
	  for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	    /* transition penalty to EL incurred when EL was entered */
	    sc = ESL_MAX(sc, ELMX(prv,cm->cp9->el_from_idx[k][c]));
	  MMX(cur,k) = ESL_MAX(sc + MSC(k), -INFTY);

	  /* E state update */
	  endsc = ESL_MAX(endsc, MMX(cur,k) + TSC(cp9O_ME,k));

	  /* insert state */
	  sc = ESL_MAX(MMX(prv,k) + TSC(cp9O_MI,k),
		       IMX(prv,k) + TSC(cp9O_II,k));
	  sc = ESL_MAX(sc, DMX(prv,k) + TSC(cp9O_DI,k));
	  IMX(cur,k) = ESL_MAX(sc + ISC(k), -INFTY);

	  /*delete state*/
	  sc = ESL_MAX(MMX(cur,k-1) + TSC(cp9O_MD,k-1),
		       IMX(cur,k-1) + TSC(cp9O_ID,k-1));
	  sc = ESL_MAX(sc, DMX(cur,k-1) + TSC(cp9O_DD, k-1));
	  DMX(cur,k) = ESL_MAX(sc, -INFTY);

	  /*el state*/
	  sc = -INFTY;
	  if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								    HMM nodes that map to right half of a MATP_MP) */
	    {
	      sc = ESL_MAX(sc, MMX(cur,k)   + TSC(cp9O_MEL,k));
	      sc = ESL_MAX(sc, ELMX(prv,k)  + el_selfsc);
	    }
	  ELMX(cur,k) = sc;

	  /*printf("mmx[jp:%d][%d]: %d %f\n", jp, k, mmx[cur][k], Score2Prob(mmx[cur][k], 1.));
	    printf("imx[jp:%d][%d]: %d %f\n", jp, k, imx[cur][k], Score2Prob(imx[cur][k], 1.));
	    printf("dmx[jp:%d][%d]: %d %f\n", jp, k, dmx[cur][k], Score2Prob(dmx[cur][k], 1.));*/
	}

      /* determine endsc = erow[cur] == sc[jp], the int score of the best parse ending at the current
       * position (j) of the target sequence. */
      endsc = ESL_MAX(endsc, DMX(cur,M) + TSC(cp9O_DM,M)); /* transition from D_M -> end */
      endsc = ESL_MAX(endsc, IMX(cur,M) + TSC(cp9O_IM,M)); /* transition from I_M -> end */
      for(c = 0; c < cm->cp9->el_from_ct[k]; c++) /* el_from_ct[k] is >= 0 */
	/* transition penalty to EL incurred when EL was entered */
	endsc = ESL_MAX(endsc, ELMX(cur,cm->cp9->el_from_idx[M+1][c]));

      erow[cur] = endsc;
      scA[jp]   = endsc;
      fsc = Scorify(endsc);
      /*printf("jp: %d j: %d fsc: %f sc: %d\n", jp, j, fsc, scA[jp]);*/
      if(fsc > best_sc)
	{
	  best_sc = fsc;
	  best_pos= j;
	}
      if (fsc > best_hmm_sc)
	{
	  best_hmm_sc = fsc;
	  best_hmm_pos= j;
	}
      if(!(cm->search_opts & CM_SEARCH_HMMGREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   */
	  gamma[jp]  = gamma[jp-1] + 0; /* extend without adding a new hit */
	  gback[jp]  = -1;
	  savesc[jp] = IMPOSSIBLE;
	  i = ((j-W+1)> i0) ? (j-W+1) : i0;
	  ip = i-i0+1;
	  curr_sc = gamma[ip-1] + fsc + cm->cp9_sc_boost;
	  /* cp9_sc_boost is experimental technique for finding hits < 0 bits. 
	   * value is 0.0 if technique not used. */
	  if (curr_sc > gamma[jp])
	    {
	      gamma[jp]  = curr_sc;
	      gback[jp]  = i;
	      savesc[jp] = fsc;
	    }
	}
      else
	{
	  /* Resolving overlaps greedily (RSEARCH style),  
	   * Return best hit for each j, IFF it's above threshold */
	  if (fsc >= cutoff) 
	    {
	      if(results != NULL) 
		{
		  i = ((j-W+1)> i0) ? (j-W+1) : i0;
		  /*printf("VIT greedy REPORTING HIT: i: %d j: %d fsc: %f\n", i, j, fsc);*/
		  report_hit (i, j, 0, fsc, results);
		  /* 0 is for saver, which is irrelevant for HMM hits */
		}
	    }
	}
    } /* end loop over end positions j */
      
  if((!(cm->search_opts & CM_SEARCH_HMMGREEDY)) && /* resolve overlaps optimally */
     (!doing_align || do_scan)) /* else we can save time by skipping traceback */
    {
      /*****************************************************************
       * Traceback stage.
       * Recover all hits: an (i,j,sc) triple for each one and report them.
       *****************************************************************/ 
      j           = j0;
      while (j >= i0) {
	jp = j-i0+1;
	if (gback[jp] == -1) /* no hit */
	  j--; 
	else                /* a hit, a palpable hit */
	  {
	    if(savesc[jp] > best_hmm_sc) 
	      {
		best_hmm_sc = savesc[jp];
		best_hmm_pos= j;
	      }
	    if(savesc[jp] >= cutoff)
	      {
		if(results != NULL) /* report the hit */
		  {
		    /*printf("VIT reporting hit: i: %d j: %d sc: %f\n", gback[jp], j, savesc[jp]);*/
		    report_hit(gback[jp], j, 0, savesc[jp], results); 
		    /* 0 is for saver, which is irrelevant for HMM hits */
		  }
	      }
	    j = gback[jp]-1;
	  }
      }
    }
  /* clean up and exit */
  free(gback);
  free(gamma);
  free(savesc);

  /* determine score to return: (I know, too complex) */
  if(doing_align)
    {
      return_sc = Scorify(scA[(j0-i0+1)]); /* L = j0-i0+1 */
      if(ret_bestpos != NULL) *ret_bestpos = i0;
      if(ret_tr != NULL) 
	{
	  /*CP9ViterbiTrace(cm->cp9, dsq, i0, j0, mx, &tr);*/
	  /* CP9PrintTrace(stdout, tr, cm->cp9, dsq); */
	  *ret_tr = tr;
	}
    }
  else if(best_hmm_sc <= 0.) /* scanning and there were no hits found by the 
			      * semi-HMM above 0.0 bits */
    {
      return_sc = best_sc;
      if(ret_bestpos != NULL) *ret_bestpos = best_pos;
    }
  else
    {
      return_sc = best_hmm_sc;
      if(ret_bestpos != NULL) *ret_bestpos = best_hmm_pos;
    }
  if(ret_sc != NULL) *ret_sc = scA;
  else free(scA);
  /*  if (ret_mx != NULL) *ret_mx = mx;*/
  /*else                FreeCPlan9Matrix(mx);*/
  /*printf("Viterbi return_sc: %f\n", return_sc);*/

  printf("ctr: %d\n", ctr);
  return return_sc;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}



/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef IMPL_CP9_FASTSEARCH_BENCHMARK
/* gcc -o benchmark-cp9_fastsearch -g -O2 -I. -L. -I../easel -L../easel -DIMPL_CP9_FASTSEARCH_BENCHMARK cp9_fastsearch.c -linfernal -leasel -lm
 * ./benchmark-cp9_fastsearch <cmfile>
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
  { "-L",        eslARG_INT, "500000", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also execute slow Viterbi scan implementation",  0 },
  { "-f",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "calculate full matrix, not just 2 rows",         0 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "configure HMM for glocal alignment [default: local]", 0 },
  { "--noel",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "turn local ends off [default: on, unless -g]", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "benchmark driver for the fast scanning CM plan 9 HMM Viterbi implementation";

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

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cmfile);
  if (!(CMFileRead(cmfp, &abc, &cm)))            cm_Fail("Failed to read CM");
  CMFileClose(cmfp);

  cm->config_opts |= CM_CONFIG_QDB;
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    if(! esl_opt_GetBoolean(go, "--noel")) cm->config_opts |= CM_CONFIG_HMMEL;
  }
  ConfigCM(cm, NULL, NULL);

  CP9_OPROFILE    *om      = NULL;
  CP9_OMX         *ox      = NULL;
  om = cp9_oprofile_Create(cm->cp9->M, cm->abc);
  cp9_oprofile_Convert(cm, om);
  if(esl_opt_GetBoolean(go, "-f")) 
    ox = cp9_omx_Create(om->M, L);
  else
    ox = cp9_omx_Create(om->M, 2);

  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, cm->null, abc->K, L, dsq);

      esl_stopwatch_Start(w);
      sc = cp9_FastViterbi(cm, om, ox, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			   TRUE,   /* we're scanning */
			   FALSE,  /* we're not ultimately aligning */
			   (! esl_opt_GetBoolean(go, "-f")),  /* memory efficient ? */
			   NULL,   /* don't want the DP matrix back */
			   NULL);  /* don't want traces back */
      printf("%4d %-30s %10.4f bits ", (i+1), "cp9_FastViterbi(): ", sc);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, " CPU time: ");
      

      if (esl_opt_GetBoolean(go, "-x")) 
	{ 
	  esl_stopwatch_Start(w);
	  sc = CP9Viterbi(cm, dsq, 1, L, cm->W, 0., NULL, NULL, NULL,
			  TRUE,   /* we're scanning */
			  FALSE,  /* we're not ultimately aligning */
			  (! esl_opt_GetBoolean(go, "-f")),  /* memory efficient ? */
			  NULL,   /* don't want the DP matrix back */
			  NULL);  /* don't want traces back */
	  printf("%4d %-30s %10.4f bits ", (i+1), "CP9Viterbi(): ", sc);
	  esl_stopwatch_Stop(w);
	  esl_stopwatch_Display(stdout, w, " CPU time: ");
	}
    }

  cp9_omx_Destroy(ox);
  cp9_oprofile_Destroy(om);
  FreeCM(cm);
  free(dsq);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IMPL_FASTSEARCH_BENCHMARK*/

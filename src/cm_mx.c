/* 
 * CM_HB_MX, ScanMatrix_t, and cm_GammaHitMx_t implementations: 
 * dynamic programming matrices for CMs
 * 
 * CM_HB_MX is based heavily on HMMER 3's p7_gmx.c module.
 *
 * Table of contents:
 *   1. CM_HB_MX data structure functions,
 *      matrix of float scores for HMM banded CM alignment/search
 *   2. ScanMatrix_t data structure functions,
 *      auxiliary info and matrix of float and/or int scores for 
 *      query dependent banded or non-banded CM DP search functions
 *   3. cm_GammaHitMx_t data structure functions,
 *      semi-HMM data structure for optimal resolution of overlapping
 *      hits for CM DP search functions
 *
 * EPN, Fri Oct 26 05:04:34 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/*****************************************************************
 *   1. CM_HB_MX data structure functions,
 *      matrix of float scores for HMM banded CM alignment/search
 *****************************************************************/

/* Function:  cm_hb_mx_Create()
 * Incept:    EPN, Fri Oct 26 05:05:07 2007
 *
 * Purpose:   Allocate a reusable, resizeable <CM_HB_MX> for a CM
 *            given a CP9Bands_t object that defines the bands. 
 *            
 *            We've set this up so it should be easy to allocate
 *            aligned memory, though we're not doing this yet.
 *
 * Returns:   a pointer to the new <CM_HB_MX>.
 *
 * Throws:    <NULL> on allocation error.
 */
CM_HB_MX *
cm_hb_mx_Create(int M)
{
  int     status;
  CM_HB_MX *mx = NULL;
  int     v;

  /* level 1: the structure itself */
  ESL_ALLOC(mx, sizeof(CM_HB_MX));
  mx->dp     = NULL;
  mx->dp_mem = NULL;
  mx->cp9b   = NULL;

  /* level 2: deck (state) pointers, 0.1..M, go all the way to M
   *          remember deck M is special, as it has no bands, we allocate
   *          it only if nec (if local ends are on) in cm_hb_mx_GrowTo()
   */
  ESL_ALLOC(mx->dp,  sizeof(float **) * (M+1));
 
  /* level 3: dp cell memory, when creating only allocate 1 cells per state, for j = 0, d = 0 */
  int allocL = 1;
  int allocW = 1;
  ESL_ALLOC(mx->dp_mem,  sizeof(float) * (M+1) * (allocL) * (allocW));
  ESL_ALLOC(mx->nrowsA, sizeof(int)      * (M+1));
  for (v = 0; v <= M; v++) {
    ESL_ALLOC(mx->dp[v], sizeof(float *) * (allocL));
    mx->nrowsA[v] = allocL;
    mx->dp[v][0]  = mx->dp_mem + v * (allocL) * (allocW);
  }
  mx->M            = M;
  mx->ncells_alloc = (M+1)*(allocL)*(allocW);
  mx->ncells_valid = 0;
  mx->L            = allocL; /* allocL = 1 */

  mx->size_Mb = 
    (float) (mx->M+1) * (float) sizeof(int *) +    /* nrowsA ptrs */
    (float) (mx->M+1) * (float) sizeof(float **) + /* mx->dp[] ptrs */
    (float) (mx->M+1) * (float) sizeof(float *) +  /* mx->dp[v][] ptrs */
    (float) mx->ncells_alloc * (float) sizeof(float); /* mx->dp_mem */
  mx->size_Mb *= 0.000001; /* convert to Mb */

  return mx;

 ERROR:
  if (mx != NULL) cm_hb_mx_Destroy(mx);
  return NULL;
}

/* Function:  cm_hb_mx_GrowTo()
 * Incept:    EPN, Fri Oct 26 05:19:49 2007
 *
 * Purpose:   Assures that a DP matrix <mx> is allocated
 *            for a model of exactly <mx->M> states and required number of 
 *            total cells. Determines new required size from 
 *            the CP9Bands_t object passed in, and reallocates if 
 *            necessary.
 *            
 *            If local ends are on (cm->flags & CMH_LOCAL_END), allocates
 *            a full non-banded EL deck.
 *
 *            Checks to make sure desired matrix isn't too big (see throws).
 *
 * Args:      cm     - the CM the matrix is for
 *            mx     - the matrix to grow
 *            errbuf - char buffer for reporting errors
 *            cp9b   - the bands for the current target sequence
 *            L      - the length of the current target sequence we're aligning
 *
 * Returns:   <eslOK> on success, and <mx> may be reallocated upon
 *            return; any data that may have been in <mx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslERANGE> if required size to grow to exceeds 
 *            CM_HB_MX_LIMIT (defined in structs.h). This
 *            should be caught and appropriately handled by caller. 
 *            <eslEINCOMPAT> on contract violation
 *            <eslEMEM> on memory allocation error.
 */
int
cm_hb_mx_GrowTo(CM_t *cm, CM_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L)
{
  int     status;
  void   *p;
  int     v, jp;
  int     cur_size = 0;
  size_t  ncells;
  int     jbw;
  double  Mb_needed;
  int     have_el;
  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;

  /* contract check, number of states (M) is something we don't change
   * so check this matrix has same number of 1st dim state ptrs that
   * cp9b has */
  if(cp9b == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_mx_GrowTo() entered with cp9b == NULL.\n");
  if(cp9b->cm_M != mx->M) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_mx_GrowTo() entered with mx->M: (%d) != cp9b->M (%d)\n", mx->M, cp9b->cm_M);
  
  ncells = 0;
  Mb_needed = ((float) (sizeof(int *)) * ((float) mx->M + 1)) + /* nrowsA ptrs */
    (float) (sizeof(float **)) * (float) (mx->M);               /* mx->dp[] ptrs */
  for(v = 0; v < mx->M; v++) { 
    jbw = cp9b->jmax[v] - cp9b->jmin[v]; 
    Mb_needed += (float) (sizeof(float *) * (jbw+1)); /* mx->dp[v][] ptrs */
    for(jp = 0; jp <= jbw; jp++) 
      ncells += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
  }
  if(have_el) ncells += (int) ((L+2) * (L+1) * 0.5); /* space for EL deck */

  Mb_needed += ESL_MAX(((float) (sizeof(float) * mx->ncells_alloc)), ((float) (sizeof(float) * ncells))); /* mx->dp_mem */
  Mb_needed *= 0.000001; /* convert to megabytes */
  ESL_DPRINTF1(("HMM banded matrix requested size: %.2f Mb\n", Mb_needed));
  if(Mb_needed > CM_HB_MX_MB_LIMIT) ESL_FAIL(eslERANGE, errbuf, "cm_hb_mx_GrowTo(), requested size of HMM banded DP matrix %.2f Mb > %.2f Mb limit (CM_HB_MX_MB_LIMIT from structs.h).", Mb_needed, (float) CM_HB_MX_MB_LIMIT);

  /* must we realloc the full matrix? or can we get away with just
   * jiggering the pointers, if total required num cells is less
   * than or equal to what we already have alloc'ed?
   */
  if (ncells > mx->ncells_alloc) {
      ESL_RALLOC(mx->dp_mem, p, sizeof(float) * ncells);
      mx->ncells_alloc = ncells;
  }
  mx->ncells_valid = ncells;

  for(v = 0; v < mx->M; v++) {
    jbw = cp9b->jmax[v] - cp9b->jmin[v] + 1;
    if(jbw > mx->nrowsA[v]) {
      ESL_RALLOC(mx->dp[v], p, sizeof(float *) * jbw);
      mx->nrowsA[v] = jbw;
    }
  }
  if(have_el) {
    jbw = L+1;
    if(jbw > mx->nrowsA[mx->M]) {
      ESL_RALLOC(mx->dp[mx->M], p, sizeof(float *) * jbw);
      mx->nrowsA[mx->M] = jbw;
    }
  }

  /* reset the pointers, we keep a tally of cur_size as we go,
   * we could precalc it and store it for each v,j, but that 
   * would be wasteful, as we'll only use the matrix configured
   * this way once, in a banded CYK run.
   */
  cur_size = 0;
  for(v = 0; v < mx->M; v++) { 
    /* mx->dp[v][0] = mx->dp_mem + cur_size; unnec, right? */
    for(jp = 0; jp <= (cp9b->jmax[v] - cp9b->jmin[v]); jp++) { 
      mx->dp[v][jp] = mx->dp_mem + cur_size;
      cur_size     += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
    }
  }
  if(have_el) {
    for(jp = 0; jp <= L; jp++) { 
      mx->dp[mx->M][jp] = mx->dp_mem + cur_size;
      cur_size     += jp + 1;
    }      
  }
  ESL_DASSERT1((cur_size == mx->ncells_valid));

  mx->cp9b = cp9b; /* just a reference */
  
  /* now update L and size_Mb */
  mx->L       = L;    /* length of current seq we're valid for */
  mx->size_Mb = Mb_needed;
  
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_hb_mx_Destroy()
 * Synopsis:  Frees a DP matrix.
 * Incept:    EPN, Fri Oct 26 09:04:04 2007
 *
 * Purpose:   Frees a <CM_HB_MX>.
 *
 * Returns:   (void)
 */
void
cm_hb_mx_Destroy(CM_HB_MX *mx)
{
  if (mx == NULL) return;
  int v;

  if (mx->dp      != NULL) { 
    for (v = 0; v <= mx->M; v++) 
      if(mx->dp[v] != NULL) free(mx->dp[v]);  
  }
  free(mx->dp);

  if (mx->nrowsA  != NULL)  free(mx->nrowsA);
  if (mx->dp_mem  != NULL)  free(mx->dp_mem);
  free(mx);
  return;
}

/* Function:  cm_hb_mx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 * Incept:    EPN, Fri Oct 26 09:04:46 2007
 *
 * Purpose:   Dump matrix <mx> to stream <fp> for diagnostics.
 */
int
cm_hb_mx_Dump(FILE *ofp, CM_HB_MX *mx)
{
  int v, jp, j, dp, d;

  fprintf(ofp, "M: %d\nL: %d\ncells_alloc: %d\nncells_valid: %d\n", mx->M, mx->L, mx->ncells_alloc, mx->ncells_valid);
  
  /* DP matrix data */
  for (v = 0; v < mx->M; v++) {
    for(jp = 0; jp <= mx->cp9b->jmax[v] - mx->cp9b->jmin[v]; jp++) {
      j = jp + mx->cp9b->jmin[v];
      for(dp = 0; dp <= mx->cp9b->hdmax[v][jp] - mx->cp9b->hdmin[v][jp]; dp++) {
	d = dp + mx->cp9b->hdmin[v][jp];
	fprintf(ofp, "dp[v:%5d][j:%5d][d:%5d] %8.4f\n", v, j, d, mx->dp[v][jp][dp]);
      }
      fprintf(ofp, "\n");
    }
    fprintf(ofp, "\n\n");
  }
  /* print EL deck, if it's valid */
  if(mx->nrowsA[mx->M] == (mx->L+1)) {
    for(j = 0; j <= mx->L; j++) {
      for(d = 0; d <= jp; d++) {
	fprintf(ofp, "dp[v:%5d][j:%5d][d:%5d] %8.4f\n", v, j, d, mx->dp[v][jp][dp]);
      }
      fprintf(ofp, "\n");
    }
    fprintf(ofp, "\n\n");
  }
  return eslOK;
}

/*****************************************************************
 *   2. ScanMatrix_t data structure functions,
 *      auxiliary info and matrix of float and/or int scores for 
 *      query dependent banded or non-banded CM DP search functions
 *****************************************************************/


/* Function: cm_CreateScanMatrix()
 * Date:     EPN, Sun Nov  4 19:56:58 2007
 *
 * Purpose:  Given a CM, allocate and initialize ScanMatrix_t object for that CM. 
 *            
 * Returns:  eslOK on success, dies immediately on some error
 */
int
cm_CreateScanMatrix(CM_t *cm, int do_float, int do_int)
{
  int status;
  int j, v;
  int do_banded;
  int use_hmmonly;
  use_hmmonly = ((cm->search_opts & CM_SEARCH_HMMVITERBI) ||  (cm->search_opts & CM_SEARCH_HMMFORWARD)) ? TRUE : FALSE;

  if(! (cm->flags & CMH_BITS))                 cm_Fail("cm_CreateScanMatrix(), the CM flag for valid bit scores is down.");
  if(cm->flags & CMH_SCANMATRIX)               cm_Fail("cm_CreateScanMatrix(), the CM flag for valid scan info is already up.");
  if(cm->smx != NULL)                           cm_Fail("cm_CreateScanMatrix, the cm already points to a ScanMatrix_t object.\n");
  if((! do_float) && (!do_int))                cm_Fail("cm_CreateScanMatrix, do_float and do_int both FALSE.\n");
  if(cm->dmin == NULL && cm->dmax != NULL)     cm_Fail("cm_CreateScanMatrix(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmin != NULL && cm->dmax == NULL)     cm_Fail("cm_CreateScanMatrix(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmax != NULL && cm->W != cm->dmax[0]) cm_Fail("cm_CreateScanMatrix(), cm->W: %d != cm->dmax[0]: %d\n", cm->W, cm->dmax[0]); 
  if((! cm->search_opts & CM_SEARCH_NOQDB) && (cm->dmin == NULL || cm->dmax == NULL))
    cm_Fail("cm_CreateScanMatrix(), cm->dmin == NULL || cm->dmax == NULL, but !(cm->search_opts & CM_SEARCH_NOQDB)\n");
  if(use_hmmonly && (cm->dmin != NULL && cm->dmax != NULL))
    cm_Fail("cm_CreateScanMatrix(), CM_SEARCH_HMMVITERBI or CM_SEARCH_HMMFORWARD but cm->dmin != NULL || cm->dmax != NULL\n");

  ScanMatrix_t *smx;
  ESL_ALLOC(smx, sizeof(ScanMatrix_t));

  smx->flags = 0;
  smx->cm_M  = cm->M;
  smx->W     = cm->W;
  smx->dmin  = cm->dmin; /* could be NULL */
  smx->dmax  = cm->dmax; /* could be NULL */
  do_banded = ((cm->search_opts & CM_SEARCH_NOQDB) || use_hmmonly) ? FALSE : TRUE;

  /* precalculate minimum and maximum d for each state and each sequence index (1..j..W). 
   * this is not always just dmin, dmax, (for ex. if j < W). */
  ESL_ALLOC(smx->dnAA, sizeof(int *) * (smx->W+1));
  ESL_ALLOC(smx->dxAA, sizeof(int *) * (smx->W+1));
  smx->dnAA[0] = smx->dxAA[0] = NULL; /* corresponds to j == 0, which is out of bounds */
  for(j = 1; j <= smx->W; j++) {
    ESL_ALLOC(smx->dnAA[j], sizeof(int) * cm->M);
    ESL_ALLOC(smx->dxAA[j], sizeof(int) * cm->M);
    for(v = 0; v < cm->M; v++) {
      if(do_banded) { 
	smx->dnAA[j][v] = (cm->sttype[v] == MP_st) ? ESL_MAX(smx->dmin[v], 2) : ESL_MAX(smx->dmin[v], 1); 
	smx->dxAA[j][v] = ESL_MIN(j, smx->dmax[v]); 
	smx->dxAA[j][v] = ESL_MIN(smx->dxAA[j][v], smx->W);
      }
      else { 
	smx->dnAA[j][v] = (cm->sttype[v] == MP_st) ? 2 : 1;
	smx->dxAA[j][v] = ESL_MIN(j, smx->W); 
      }
    }
  }
  /* allocate bestr, which holds best root state at alpha[0][cur][d] */
  ESL_ALLOC(smx->bestr, (sizeof(int) * (smx->W+1)));

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
   */

  smx->falpha       = NULL;
  smx->falpha_begl  = NULL;

  smx->ialpha      = NULL;
  smx->ialpha_begl = NULL;

  cm->smx = smx;
  if(do_float) /* allocate float mx and scores */
    cm_FloatizeScanMatrix(cm);
  if(do_int)   /* allocate int mx and scores */
    cm_IntizeScanMatrix(cm);
  cm->flags |= CMH_SCANMATRIX; /* raise the flag for valid CMH_SCANMATRIX */
  return eslOK;

 ERROR:
  cm_Fail("memory allocation error in cm_CreateScanMatrix().\n");
  return status; /* NEVERREACHED */
}


/* Function: cm_UpdateScanMatrix()
 * Date:     EPN, Wed Nov  7 12:49:36 2007
 *
 * Purpose:  Free, reallocate and recalculate the ScanMatrix for a CM.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_UpdateScanMatrix(CM_t *cm)
{
  /* contract check */
  if(cm->flags & CMH_SCANMATRIX)    cm_Fail("cm_UpdateScanMatrix(), the CM flag for valid scan info is already up.");
  if(cm->smx->flags & cmSMX_HAS_FLOAT) {
    cm_FreeFloatsFromScanMatrix(cm);
    cm_FloatizeScanMatrix(cm);
  }
  if(cm->smx->flags & cmSMX_HAS_INT) {
    cm_FreeIntsFromScanMatrix(cm);
    cm_IntizeScanMatrix(cm);
  }
  cm->flags |= CMH_SCANMATRIX; /* ScanMatrix is valid now */
  return eslOK;
}

/* Function: cm_FloatizeScanMatrix()
 * Date:     EPN, Wed Nov  7 10:05:55 2007
 *
 * Purpose:  Allocate and initialize float data structures in a ScanMatrix_t object for <cm>.
 *           This initializes a scanning float DP matrix for CYK/Inside, for details on that
 *           matrix see the notes by the cm_FloatizeScanMatrix() function call in 
 *           cm_CreateScanMatrix().
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FloatizeScanMatrix(CM_t *cm)
{
  int status;
  int j, v;
  int d, y, yoffset, w;
  int use_hmmonly;
  use_hmmonly = ((cm->search_opts & CM_SEARCH_HMMVITERBI) ||  (cm->search_opts & CM_SEARCH_HMMFORWARD)) ? TRUE : FALSE;
  int do_banded = ((cm->search_opts & CM_SEARCH_NOQDB) || use_hmmonly) ? FALSE : TRUE;

  /* contract check */
  if(cm->smx == NULL) cm_Fail("cm_FloatizeScanMatrix(), cm->smx is NULL.\n");
  ScanMatrix_t *smx = cm->smx;
  if(smx->flags & cmSMX_HAS_FLOAT) cm_Fail("cm_FloatizeScanMatrix(), si's cmSMX_HAS_FLOAT flag is already up.");
  if(smx->falpha != NULL)       cm_Fail("cm_FloatizeScanMatrix(), smx->falpha is not NULL.");
  if(smx->falpha_begl != NULL)  cm_Fail("cm_FloatizeScanMatrix(), smx->falpha_begl is not NULL.");
  
  /* allocate alpha */
  ESL_ALLOC(smx->falpha,        sizeof(float **) * 2);
  ESL_ALLOC(smx->falpha[0],     sizeof(float *) * cm->M);
  ESL_ALLOC(smx->falpha[1],     sizeof(float *) * cm->M);
  ESL_ALLOC(smx->falpha[0][0],  sizeof(float) * 2 * (cm->M) * (smx->W+1));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      smx->falpha[0][v] = smx->falpha[0][0] + (v           * (smx->W+1));
      smx->falpha[1][v] = smx->falpha[0][0] + ((v + cm->M) * (smx->W+1));
    }
    else smx->falpha[0][v] = smx->falpha[1][v] = NULL; /* BEGL_S */
  }
  /* allocate falpha_begl */
  ESL_ALLOC(smx->falpha_begl, (sizeof(float **) * (smx->W+1)));
  for (j = 0; j <= smx->W; j++) {
    ESL_ALLOC(smx->falpha_begl[j], (sizeof(float *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) ESL_ALLOC(smx->falpha_begl[j][v], (sizeof(float) * (smx->W+1)));
      else smx->falpha_begl[j][v] = NULL; /* non-BEGL */
    }
  }
  /* initialize falpha and falpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      smx->falpha[0][v][0] = IMPOSSIBLE;
      if (cm->sttype[v] == E_st) { 
	smx->falpha[0][v][0] = smx->falpha[1][v][0] = 0.;
	/* rest of E deck is IMPOSSIBLE, this is rewritten if QDB is on, (slightly wasteful). */
	for (d = 1; d <= smx->W; d++) smx->falpha[0][v][d] = smx->falpha[1][v][d] = IMPOSSIBLE;
      }
      else if (cm->sttype[v] == MP_st) smx->falpha[0][v][1] = smx->falpha[1][v][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	smx->falpha[0][v][0] = cm->endsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  smx->falpha[0][v][0] = ESL_MAX(smx->falpha[0][v][0], (smx->falpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	smx->falpha[0][v][0] = ESL_MAX(smx->falpha[0][v][0], IMPOSSIBLE);
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	smx->falpha[0][v][0] = smx->falpha_begl[0][w][0] + smx->falpha[0][y][0]; 
      }
      smx->falpha[1][v][0] = smx->falpha[0][v][0];
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      smx->falpha_begl[0][v][0] = cm->endsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	smx->falpha_begl[0][v][0] = ESL_MAX(smx->falpha_begl[0][v][0], (smx->falpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in smx->falpha */
      smx->falpha_begl[0][v][0] = ESL_MAX(smx->falpha_begl[0][v][0], IMPOSSIBLE);
      for (j = 1; j <= smx->W; j++) 
	smx->falpha_begl[j][v][0] = smx->falpha_begl[0][v][0];
    }
  }
  /* query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) { 
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < cm->dmin[v] && d <= smx->W; d++)
	  for(j = 0; j < 2; j++) 
	    smx->falpha[j][v][d]  = IMPOSSIBLE;
	for (d = cm->dmax[v]+1; d <= smx->W; d++)
	  for(j = 0; j < 2; j++)
	    smx->falpha[j][v][d] = IMPOSSIBLE;
      }
      else { /* not BEGL_S state */
	for (d = 0; d < cm->dmin[v] && d <= smx->W; d++)
	  for(j = 0; j <= smx->W; j++)
	    smx->falpha_begl[j][v][d] = IMPOSSIBLE;
	for (d = cm->dmax[v]+1; d <= smx->W; d++)
	  for(j = 0; j <= smx->W; j++)
	    smx->falpha_begl[j][v][d] = IMPOSSIBLE;
      }
    }
  }
  /* set the flag that tells us we've got valid floats */
  smx->flags |= cmSMX_HAS_FLOAT;
  return eslOK;

 ERROR: 
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: cm_IntizeScanMatrix()
 * Date:     EPN, Wed Nov  7 10:10:39 2007
 *
 * Purpose:  Allocate and initialize int data structures in a ScanMatrix_t object for <cm>.
 *           This initializes a scanning float DP matrix for CYK/Inside, for details on that
 *           matrix see the notes by the cm_FloatizeScanMatrix() function call in 
 *           cm_CreateScanMatrix().
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_IntizeScanMatrix(CM_t *cm)
{
  int status;
  int v, j, d, y, yoffset, w;
  int use_hmmonly;
  use_hmmonly = ((cm->search_opts & CM_SEARCH_HMMVITERBI) ||  (cm->search_opts & CM_SEARCH_HMMFORWARD)) ? TRUE : FALSE;
  int do_banded = ((cm->search_opts & CM_SEARCH_NOQDB) || use_hmmonly) ? FALSE : TRUE;

  /* contract check */
  if(cm->smx == NULL) cm_Fail("cm_IntizeScanMatrix(), cm->smx is NULL.\n");
  ScanMatrix_t *smx = cm->smx;
  if(smx->flags & cmSMX_HAS_INT) cm_Fail("cm_IntizeScanMatrix(), si's cmSMX_HAS_INT flag is already up.");
  if(smx->ialpha != NULL)       cm_Fail("cm_IntizeScanMatrix(), smx->ialpha is not NULL.");
  if(smx->ialpha_begl != NULL)  cm_Fail("cm_IntizeScanMatrix(), smx->ialpha_begl is not NULL.");

  /* allocate ialpha */
  ESL_ALLOC(smx->ialpha,        sizeof(int **) * 2);
  ESL_ALLOC(smx->ialpha[0],     sizeof(int *) * cm->M);
  ESL_ALLOC(smx->ialpha[1],     sizeof(int *) * cm->M);
  ESL_ALLOC(smx->ialpha[0][0],  sizeof(int) * 2 * (cm->M) * (smx->W+1));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      smx->ialpha[0][v] = smx->ialpha[0][0] + (v           * (smx->W+1));
      smx->ialpha[1][v] = smx->ialpha[0][0] + ((v + cm->M) * (smx->W+1));
    }
    else smx->ialpha[0][v] = smx->ialpha[1][v] = NULL; /* BEGL_S */
  }
  /* allocate ialpha_begl */
  ESL_ALLOC(smx->ialpha_begl, (sizeof(int **) * (smx->W+1)));
  for (j = 0; j <= cm->W; j++) {
    ESL_ALLOC(smx->ialpha_begl[j], (sizeof(int *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) ESL_ALLOC(smx->ialpha_begl[j][v], (sizeof(int) * (smx->W+1)));
      else smx->ialpha_begl[j][v] = NULL; /* non-BEGL */
    }
  }
  /* initialize ialpha and ialpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      smx->ialpha[0][v][0] = -INFTY;
      if (cm->sttype[v] == E_st) { 
	smx->ialpha[0][v][0] = smx->ialpha[1][v][0] = 0.;
	/* rest of E deck is -INFTY, this is rewritten if QDB is on, (slightly wasteful). */
	for (d = 1; d <= smx->W; d++) smx->ialpha[0][v][d] = smx->ialpha[1][v][d] = -INFTY;
      }
      else if (cm->sttype[v] == MP_st) smx->ialpha[0][v][1] = smx->ialpha[1][v][1] = -INFTY;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	smx->ialpha[0][v][0] = cm->iendsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  smx->ialpha[0][v][0] = ESL_MAX(smx->ialpha[0][v][0], (smx->ialpha[0][y+yoffset][0] + cm->itsc[v][yoffset]));
	smx->ialpha[0][v][0] = ESL_MAX(smx->ialpha[0][v][0], -INFTY);
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	smx->ialpha[0][v][0] = smx->ialpha_begl[0][w][0] + smx->ialpha[0][y][0]; 
      }
      smx->ialpha[1][v][0] = smx->ialpha[0][v][0];
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      smx->ialpha_begl[0][v][0] = cm->iendsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	smx->ialpha_begl[0][v][0] = ESL_MAX(smx->ialpha_begl[0][v][0], (smx->ialpha[0][y+yoffset][0] + cm->itsc[v][yoffset])); /* careful: y is in alpha */
      smx->ialpha_begl[0][v][0] = ESL_MAX(smx->ialpha_begl[0][v][0], -INFTY);
      for (j = 1; j <= smx->W; j++) 
	smx->ialpha_begl[j][v][0] = smx->ialpha_begl[0][v][0];
    }
  }
  /* query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) { 
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < cm->dmin[v] && d <= smx->W; d++)
	  for(j = 0; j < 2; j++) 
	    smx->ialpha[j][v][d] = -INFTY;
	for (d = cm->dmax[v]+1; d <= smx->W; d++)
	  for(j = 0; j < 2; j++)
	    smx->ialpha[j][v][d] = -INFTY;
      }
      else { /* not BEGL_S state */
	for (d = 0; d < cm->dmin[v] && d <= smx->W; d++)
	  for(j = 0; j <= smx->W; j++)
	    smx->ialpha_begl[j][v][d] = -INFTY;
	for (d = cm->dmax[v]+1; d <= smx->W; d++)
	  for(j = 0; j <= smx->W; j++)
	    smx->ialpha_begl[j][v][d] = -INFTY;
      }
    }
  }
  /* set the flag that tells us we've got valid ints */
  smx->flags |= cmSMX_HAS_INT;
  return eslOK;

 ERROR: 
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}  


/* Function: cm_FreeFloatsFromScanMatrix()
 * Date:     EPN, Wed Nov  7 10:03:55 2007
 *
 * Purpose:  Free float data structures in a ScanMatrix_t object for <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeFloatsFromScanMatrix(CM_t *cm)
{
  int j, v;

  /* contract check */
  if(cm->smx == NULL) cm_Fail("cm_FreeFloatsFromScanMatrix(), cm->smx is NULL.\n");
  ScanMatrix_t *smx = cm->smx;
  if(! smx->flags & cmSMX_HAS_FLOAT)    cm_Fail("cm_FreeFloatsFromScanMatrix(), si's cmSMX_HAS_FLOAT flag is down.");
  if(smx->falpha == NULL)       cm_Fail("cm_FreeFloatsFromScanMatrix(), smx->falpha is already NULL.");
  if(smx->falpha_begl == NULL)  cm_Fail("cm_FreeFloatsFromScanMatrix(), smx->falpha_begl is already NULL.");

  free(smx->falpha[0][0]);
  free(smx->falpha[1]);
  free(smx->falpha[0]);
  free(smx->falpha);
  smx->falpha = NULL;
  for (j = 0; j <= smx->W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S) {
	free(smx->falpha_begl[j][v]);
      }
    free(smx->falpha_begl[j]);
  }
  free(smx->falpha_begl);
  smx->falpha_begl = NULL;
  smx->flags &= ~cmSMX_HAS_FLOAT;
  return eslOK;
}

/* Function: cm_FreeIntsFromScanMatrix()
 * Date:     EPN, Wed Nov  7 09:56:01 2007
 *
 * Purpose:  Free int data structures in a ScanMatrix_t object for <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeIntsFromScanMatrix(CM_t *cm)
{
  int j, v;

  /* contract check */
  if(cm->smx == NULL) cm_Fail("cm_FreeIntsFromScanMatrix(), cm->smx is NULL.\n");
  ScanMatrix_t *smx = cm->smx;
  if(! smx->flags & cmSMX_HAS_INT)    cm_Fail("cm_FreeIntsFromScanMatrix(), si's cmSMX_HAS_INT flag is down.");
  if(smx->ialpha == NULL)       cm_Fail("cm_FreeIntsFromScanMatrix(), smx->ialpha is already NULL.");
  if(smx->ialpha_begl == NULL)  cm_Fail("cm_FreeIntsFromScanMatrix(), smx->ialpha_begl is already NULL.");

  free(smx->ialpha[0][0]);
  free(smx->ialpha[1]);
  free(smx->ialpha[0]);
  free(smx->ialpha);
  smx->ialpha = NULL;
  for (j = 0; j <= smx->W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S) {
	free(smx->ialpha_begl[j][v]);
      }
    free(smx->ialpha_begl[j]);
  }
  free(smx->ialpha_begl);
  smx->ialpha_begl = NULL;
  smx->flags &= ~cmSMX_HAS_INT;
  return eslOK;
}

/* Function: cm_FreeScanMatrix()
 * Date:     EPN, Sun Nov  4 20:57:32 2007
 *
 * Purpose:  Free a ScanMatrix_t object for <cm>.
 *            
 * Returns:  void
 */
void
cm_FreeScanMatrix(CM_t *cm)
{
  int j;
  /* contract check */
  if(cm->smx == NULL) cm_Fail("cm_FreeScanMatrix(), cm->smx is NULL.\n");
  ScanMatrix_t *smx = cm->smx;
  for(j = 1; j <= smx->W; j++) {
    free(smx->dnAA[j]);
    free(smx->dxAA[j]);
  }
  free(smx->dnAA);
  free(smx->dxAA);
  free(smx->bestr);
  
  if(smx->flags & cmSMX_HAS_FLOAT) cm_FreeFloatsFromScanMatrix(cm);
  if(smx->flags & cmSMX_HAS_INT)   cm_FreeIntsFromScanMatrix(cm);
  free(smx);
  cm->smx = NULL;
  cm->flags &= ~CMH_SCANMATRIX; /* drop the 'cm has valid scanmatrix flag */
  return;
}

/* Function: cm_DumpScanMatrixAlpha()
 * Date:     EPN, Tue Nov  6 05:11:26 2007
 *
 * Purpose:  Dump current alpha matrix (either float or int).
 *            
 * Returns:  void.
 */
void
cm_DumpScanMatrixAlpha(CM_t *cm, int j, int i0, int doing_float)
{
  int d, v;
  int jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
  int cur = j%2;
  int prv = (j-1)%2;
  int *dnA, *dxA;

  if(cm->smx == NULL) cm_Fail("cm_DumpScanMatrixAlpha(), cm->smx is NULL.\n");
  ScanMatrix_t *smx = cm->smx;
  if(doing_float && (! smx->flags & cmSMX_HAS_FLOAT)) cm_Fail("cm_DumpScanMatrixAlpha(), trying to print float alpha, but cmSMX_HAS_FLOAT flag is down.\n");
  if((! doing_float) && (! smx->flags & cmSMX_HAS_INT)) cm_Fail("cm_DumpScanMatrixAlpha(), trying to print int alpha, but cmSMX_HAS_INT flag is down.\n");

  int begl_prv = j-1 % (smx->W+1);
  int begl_cur = j   % (smx->W+1);

  printf("Dumping Alpha: j: %d\n", j);
  if(jp_g >= smx->W) { dnA = smx->dnAA[smx->W]; dxA = smx->dxAA[smx->W]; }
  else              { dnA = smx->dnAA[jp_g];  dxA = smx->dxAA[jp_g]; }
  if(doing_float) {
    for (v = smx->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) { 
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10.4f\n", (j-1), v, d, smx->falpha_begl[begl_prv][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10.4f\n", (j-1), v, d, smx->falpha[prv][v][d]); 
      }
      printf("\n");
    }
    for (v = smx->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10.4f\n", j,     v, d, smx->falpha_begl[begl_cur][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10.4f\n", j,     v, d, smx->falpha[cur][v][d]); 
      }
      printf("\n");
    }
  }
  else { /* doing int */
    for (v = smx->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10d\n", (j-1), v, d, smx->ialpha_begl[begl_prv][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10d\n", (j-1), v, d, smx->ialpha[prv][v][d]); 
      }
      printf("\n\n");
    }
    for (v = smx->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10d\n", j,     v, d, smx->ialpha_begl[begl_cur][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10d\n", j,     v, d, smx->ialpha[cur][v][d]); 
      }
      printf("\n\n");
    }
  }
  return;
}

/* Function: FCalcOptimizedEmitScores()
 * Date:     EPN, Tue Nov  6 17:24:45 2007
 *
 * Purpose:  Allocate, fill and return an optimized emission score vector
 *           of float scores for fast search/alignment.
 *            
 * Returns:  the 2D float emission score vector on success,
 *           dies immediately on memory allocation error.
 */
float **
FCalcOptimizedEmitScores(CM_t *cm)
{
  int status; 
  float **esc_vAA;
  int a,b;
  int v;

  /* precalculate possible emission scores for each state */
  ESL_ALLOC(esc_vAA,  sizeof(float *) * (cm->M));
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      ESL_ALLOC(esc_vAA[v],  sizeof(float) * (cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v],  cm->abc->Kp, IMPOSSIBLE);
      for(a = 0; a < cm->abc->K; a++)
	esc_vAA[v][a]  = cm->esc[v][a];
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) /* note boundary conditions, gap, missing data symbols stay IMPOSSIBLE */
	esc_vAA[v][a]  = esl_abc_FAvgScore(cm->abc, a, cm->esc[v]);
      break;
    case MP_st:
      ESL_ALLOC(esc_vAA[v], sizeof(float) * (cm->abc->Kp * cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_FSet(esc_vAA[v],  (cm->abc->Kp * cm->abc->Kp), IMPOSSIBLE);
      for(a = 0; a < (cm->abc->Kp-1); a++)
	for(b = 0; b < (cm->abc->Kp-1); b++)
	  if(a < cm->abc->K && b < cm->abc->K) 
	    esc_vAA[v][a * cm->abc->Kp + b]  = cm->esc[v][(a * cm->abc->K) + b];
	  else
	    esc_vAA[v][a  * cm->abc->Kp + b]  = DegeneratePairScore(cm->abc, cm->esc[v], a, b);
      break;
    default:
      esc_vAA[v] = NULL;
      break;
    }
  }
  return esc_vAA;

 ERROR:
  cm_Fail("memory allocation error.");
  return NULL; /* NEVERREACHED */
}


/* Function: ICalcOptimizedEmitScores()
 * Date:     EPN, Tue Nov  6 17:27:34 2007
 *
 * Purpose:  Allocate, fill and return an optimized emission score vector
 *           of integer scores for fast search/alignment.
 *            
 * Returns:  the 2D integer emission score vector on success,
 *           dies immediately on memory allocation error.
 */
int **
ICalcOptimizedEmitScores(CM_t *cm)
{
  int status;
  int **iesc_vAA;
  int v;
  int a,b;

  /* precalculate possible emission scores for each state */
  ESL_ALLOC(iesc_vAA, sizeof(int *)   * (cm->M));
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      ESL_ALLOC(iesc_vAA[v], sizeof(int)   * (cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_ISet(iesc_vAA[v], cm->abc->Kp, -INFTY);
      for(a = 0; a < cm->abc->K; a++) 
	iesc_vAA[v][a] = cm->iesc[v][a];
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) /* note boundary conditions, gap, missing data symbols stay IMPOSSIBLE */
	iesc_vAA[v][a] = esl_abc_IAvgScore(cm->abc, a, cm->iesc[v]);
      break;
    case MP_st:
      ESL_ALLOC(iesc_vAA[v], sizeof(int) * (cm->abc->Kp * cm->abc->Kp));
      /* ALLOCATE SIZE = POWER OF 2? */
      esl_vec_ISet(iesc_vAA[v], (cm->abc->Kp * cm->abc->Kp), -INFTY);
      for(a = 0; a < (cm->abc->Kp-1); a++)
	for(b = 0; b < (cm->abc->Kp-1); b++)
	  if(a < cm->abc->K && b < cm->abc->K) 
	    iesc_vAA[v][a * cm->abc->Kp + b] = cm->iesc[v][(a * cm->abc->K) + b];
	  else 
	    iesc_vAA[v][a * cm->abc->Kp + b] = iDegeneratePairScore(cm->abc, cm->iesc[v], a, b);

      break;
    default:
      iesc_vAA[v] = NULL;
      break;
    }
  }
  return iesc_vAA;

 ERROR:
  cm_Fail("memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: FreeOptimizedEmitScores()
 * Date:     EPN, Fri Nov  9 08:44:06 2007
 *
 * Purpose:  Free 2D vectors of optimized emissions scores.
 *           Either fesc_vAA or iesc_vAA (or both) must be non-NULL.
 *            
 * Returns:  void
 */
void
FreeOptimizedEmitScores(float **fesc_vAA, int **iesc_vAA, int M)
{
  int v;
  if(fesc_vAA == NULL && iesc_vAA == NULL) cm_Fail("FreeOptimizedEmitScores() but fesc and iesc are NULL.\n");

  if(fesc_vAA != NULL) { 
    for(v = 0; v < M; v++) 
      if(fesc_vAA[v] != NULL) free(fesc_vAA[v]);
    free(fesc_vAA);
    fesc_vAA = NULL;
  }

  if(iesc_vAA != NULL) { 
    for(v = 0; v < M; v++) 
      if(iesc_vAA[v] != NULL) free(iesc_vAA[v]);
    free(iesc_vAA);
    iesc_vAA = NULL;
  }
  return;
}


/* Function: FCalcInitDPScores()
 * Date:     EPN, Fri Nov  9 09:18:07 2007
 *
 * Purpose:  Allocate, fill and return the initial float scores
 *           for a scanning DP matrix for CM <cm> as it's
 *           currently configured. All [0..v..M-1][0..d..W]
 *           cells are allocated and filled, it's up to 
 *           the DP function to ignore cells outside bands.
 *            
 * Returns:  the 2D float init sc vector on success,
 *           dies immediately on memory error.
 */    
float **
FCalcInitDPScores(CM_t *cm)
{
  int status;
  float *el_scA;
  float **init_scAA;
  int v, d;

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(float) * (cm->W+1));
  for(d = 0; d <= cm->W; d++) el_scA[d] = cm->el_selfsc * d;
  /* precalculate the initial score for all alpha[v][j][d] cells, it's independent of j 
   * these scores ignore bands, (that is cells outside bands still have initsc's calc'ed)
   * it's up to the DP function to skip these cells. */
  ESL_ALLOC(init_scAA,    sizeof(float *) * (cm->M));
  ESL_ALLOC(init_scAA[0], sizeof(float)   * (cm->M) * (cm->W+1));
  for (v = 0; v < cm->M; v++) {
    init_scAA[v] = init_scAA[0] + (v * (cm->W+1));
    if(NOT_IMPOSSIBLE(cm->endsc[v])) {
      for(d = 0; d <= cm->W; d++) init_scAA[v][d] = el_scA[d] + cm->endsc[v];
    }
    else {
      for(d = 0; d <= cm->W; d++) init_scAA[v][d] = IMPOSSIBLE;
    }
  }
  free(el_scA);
  return init_scAA;

 ERROR:
  cm_Fail("memory allocation error.");
  return NULL; /* NEVERREACHED */
}


/* Function: ICalcInitDPScores()
 * Date:     EPN, Fri Nov  9 09:10:33 2007
 *
 * Purpose:  Allocate, fill and return the initial int scores
 *           for a scanning DP matrix for CM <cm> as it's
 *           currently configured. All [0..v..M-1][0..d..W]
 *           cells are allocated and filled, it's up to 
 *           the DP function to ignore cells outside bands.
 *            
 * Returns:  the 2D integer init sc vector on success,
 *           dies immediately on memory error.
 */    
int **
ICalcInitDPScores(CM_t *cm)
{
  int status;
  int *el_scA;
  int **init_scAA;
  int v, d;

  /* precalcuate all possible local end scores, for local end emits of 1..W residues */
  ESL_ALLOC(el_scA, sizeof(int) * (cm->W+1));
  for(d = 0; d <= cm->W; d++) el_scA[d] = cm->iel_selfsc * d;
  /* precalculate the initial score for all ialpha[v][j][d] cells, it's independent of j 
   * these scores ignore bands, (that is cells outside bands still have initsc's calc'ed)
   * it's up to the DP function to skip these cells. */
  ESL_ALLOC(init_scAA,    sizeof(int *) * (cm->M));
  ESL_ALLOC(init_scAA[0], sizeof(int)   * (cm->M) * (cm->W+1)); 
  for (v = 0; v < cm->M; v++) {
    init_scAA[v] = init_scAA[0] + (v * (cm->W+1));
    if(cm->iendsc[v] != -INFTY) {
      for(d = 0; d <= cm->W; d++) init_scAA[v][d] = el_scA[d] + cm->iendsc[v];
    }
    else {
      for(d = 0; d <= cm->W; d++) init_scAA[v][d] = -INFTY;
    }
  }


  free(el_scA);
  return init_scAA;

 ERROR:
  cm_Fail("memory allocation error.");
  return NULL; /* NEVERREACHED */
}


/*****************************************************************
 *   3. cm_GammaHitMx_t data structure functions,
 *      Semi HMM data structure for optimal resolution of overlapping
 *      hits for CM DP search functions.
 *****************************************************************/
  
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

/* Function: cm_UpdateGammaHitMx()
 * Date:     EPN, Mon Nov  5 05:41:14 2007
 *
 * Purpose:  Update a gamma semi-HMM for CM hits that end at gamma-relative position <j>.
 * Note:     Only difference with cm_UpdateIntGammaHitMx() is alpha_row is floats, not ints.
 *
 * Args:     gamma     - the gamma data structure
 *           j         - offset j for gamma must be between 0 and gamma->L
 *           alpha_row - row of DP matrix to examine, we look at [dn..dx], NULL if we want to report
 *                       this j is IMPOSSIBLE end point of a hit (only possible if using_hmm_bands == TRUE)
 *           dn        - minimum d to look at 
 *           dx        - maximum d to look at
 *           using_hmm_bands - if TRUE, alpha_row is offset by dn, so we look at [0..dx-dn]
 *           bestr     - [dn..dx] root state (0 or local entry) corresponding to hit stored in alpha_row
 *           doing_inside - if TRUE, we don't store bestr, we've summed over all possible starts
 *           results   - results to add to, only used in this function if gamma->iamgreedy 
 *
 * Returns:  void;

 */
void
cm_UpdateGammaHitMx(cm_GammaHitMx_t *gamma, int j, float *alpha_row, int dn, int dx, int using_hmm_bands, 
		    int *bestr, int doing_inside, search_results_t *results)
{
  int i, d;
  float sc;
  int bestd;
  int r;
  int dmin, dmax;
  int ip, jp;

  if(alpha_row == NULL && (!using_hmm_bands)) cm_Fail("cm_UpdateGammaHitMx(), alpha_row is NULL, but using_hmm_bands is FALSE.\n");
  dmin = (using_hmm_bands) ? 0     : dn;
  dmax = (using_hmm_bands) ? dx-dn : dx;

  /* mode 1: non-greedy  */
  if(! gamma->iamgreedy || alpha_row == NULL) { 
    gamma->mx[j]     = gamma->mx[j-1] + 0; 
    gamma->gback[j]  = -1;
    gamma->savesc[j] = IMPOSSIBLE;
    gamma->saver[j]  = -1;

    if(alpha_row != NULL) { 
      for (d = dmin; d <= dmax; d++) {
	i = using_hmm_bands ? j-d+1-dn  : j-d+1;
	sc = gamma->mx[i-1] + alpha_row[d];
	if (sc > gamma->mx[j]) {
	  gamma->mx[j]     = sc;
	  gamma->gback[j]  = i;
	  gamma->savesc[j] = alpha_row[d]; 
	  gamma->saver[j]  = doing_inside ? -1 : bestr[d]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
	}
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
    /* First, report hit with d of dmin (min valid d) if >= cutoff */
    if (alpha_row[dmin] >= gamma->cutoff) {
      r = doing_inside ? -1 : bestr[dmin]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
      ip = using_hmm_bands ? j-(dmin+dn)+gamma->i0 : j-dmin+gamma->i0;
      jp = j-1+gamma->i0;
      report_hit (ip, jp, r, alpha_row[dmin], results);
    }
    bestd = dmin;
    /* Now, if current score is greater than maximum seen previous, report
       it if >= cutoff and set new max */
    for (d = dmin+1; d <= dmax; d++) {
      if (alpha_row[d] > alpha_row[bestd]) {
	if (alpha_row[d] >= gamma->cutoff) { 
	  r = doing_inside ? -1 : bestr[d]; /* saver/bestr is invalid for Inside, we've summed all parses, none of this single parse crap */
	  ip = using_hmm_bands ? j-(d+dn)+gamma->i0 : j-d+gamma->i0;
	  jp = j-1+gamma->i0;
	  report_hit (ip, jp, r, alpha_row[d], results);
	}
	bestd = d;
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

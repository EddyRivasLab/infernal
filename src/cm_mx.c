/* CM_HB_MX, ScanMatrix_t, and GammaHitMx_t implementations: 
 * dynamic programming matrices for CMs
 * 
 * CM_HB_MX is based heavily on HMMER 3's p7_gmx.c module.
 *
 * Table of contents:
 *   1. CM_HB_MX data structure functions,
 *      matrix of float scores for HMM banded CM alignment/search
 *   2. CM_HB_SHADOW_MX data structure functions
 *      HMM banded shadow matrix for tracing back HMM banded CM parses
 *   3. ScanMatrix_t data structure functions,
 *      auxiliary info and matrix of float and/or int scores for 
 *      query dependent banded or non-banded CM DP search functions
 *   4. GammaHitMx_t data structure functions,
 *      semi-HMM data structure for optimal resolution of overlapping
 *      hits for CM and CP9 DP search functions
 *
 * EPN, Fri Oct 26 05:04:34 2007
 * SVN $Id$
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

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
  int allocL = 1;
  int allocW = 1;

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
 
  /* level 3: dp cell memory, when creating only allocate 1 cell per state, for j = 0, d = 0 */
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
 *            size_limit- max number of Mb for DP matrix, if matrix is bigger -> return eslERANGE
 *
 * Returns:   <eslOK> on success, and <mx> may be reallocated upon
 *            return; any data that may have been in <mx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslERANGE> if required size to grow to exceeds <size_limit>.
 *            This should be caught and appropriately handled by caller. 
 *            <eslEINCOMPAT> on contract violation
 *            <eslEMEM> on memory allocation error.
 */
int
cm_hb_mx_GrowTo(CM_t *cm, CM_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit)
{
  int     status;
  void   *p;
  int     v, jp;
  int     cur_size = 0;
  int64_t ncells;
  int     jbw;
  float   Mb_needed;
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
  ESL_DPRINTF2(("HMM banded matrix requested size: %.2f Mb\n", Mb_needed));
  if(Mb_needed > size_limit) ESL_FAIL(eslERANGE, errbuf, "requested HMM banded DP mx of %.2f Mb > %.2f Mb limit.\nIncrease limit with --mxsize or tau with --tau.", Mb_needed, (float) size_limit);

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

  fprintf(ofp, "M: %d\nL: %d\ncells_alloc: %" PRId64 "\nncells_valid: %" PRId64 "\n", mx->M, mx->L, mx->ncells_alloc, mx->ncells_valid);
  
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


/* Function:  cm_hb_mx_NumCellsNeeded()
 * Incept:    EPN, Fri Oct 29 13:03:05 2010
 *
 * Purpose:   Given a model and CP9_bands_t object with pre-calced bands for 
 *            a target, determine the number of cells required in a 
 *            CM_HB_MX for the target given the bands.
 *            
 *
 * Args:      cm     - the CM the matrix is for
 *            mx     - the matrix to grow
 *            errbuf - char buffer for reporting errors
 *            cp9b   - the bands for the current target sequence
 *            L      - the length of the current target sequence we're aligning
 *            size_limit- max number of Mb for DP matrix, if matrix is bigger -> return eslERANGE
 *
 * Returns:   <eslOK> on success
 *
 */
int
cm_hb_mx_NumCellsNeeded(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_ncells)
{
  int     v, jp;
  int64_t ncells;
  int     jbw;
  int     have_el;

  have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;

  /* contract check */
  if(cp9b == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_mx_GrowTo() entered with cp9b == NULL.\n");

  ncells = 0;
  for(v = 0; v < cp9b->cm_M; v++) { 
    jbw = cp9b->jmax[v] - cp9b->jmin[v]; 
    for(jp = 0; jp <= jbw; jp++) 
      ncells += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
  }
  if(have_el) ncells += (int) ((L+2) * (L+1) * 0.5); /* space for EL deck */

  *ret_ncells = ncells;
  return eslOK;
}


/*****************************************************************
 *   2. CM_HB_SHADOOW_MX data structure functions,
 *      HMM banded shadow matrix for tracing back HMM banded CM parses
 *****************************************************************/

/* Function:  cm_hb_shadow_mx_Create()
 * Incept:    EPN, Fri Oct 26 05:05:07 2007
 *
 * Purpose:   Allocate a reusable, resizeable <CM_HB_SHADOW_MX> for a CM <cm>
 *            given a CP9Bands_t object that defines the bands. The CM
 *            is needed so we know which decks need to be int's (BIF_B states)
 *            and which need to be char's (all other states).
 *            
 *            We've set this up so it should be easy to allocate
 *            aligned memory, though we're not doing this yet.
 *
 * Returns:   a pointer to the new <CM_HB_SHADOW_MX>.
 *
 * Throws:    <NULL> on allocation error.
 */
CM_HB_SHADOW_MX *
cm_hb_shadow_mx_Create(CM_t *cm, int M)
{
  int     status;
  CM_HB_SHADOW_MX *mx = NULL;
  int     v;
  int     nbifs, nb;
  int allocL = 1;
  int allocW = 1;

  /* level 1: the structure itself */
  ESL_ALLOC(mx, sizeof(CM_HB_SHADOW_MX));
  mx->yshadow     = NULL;
  mx->yshadow_mem = NULL;
  mx->kshadow     = NULL;
  mx->kshadow_mem = NULL;
  mx->cp9b        = NULL;

  /* level 2: deck (state) pointers, 0.1..M-1, M (EL deck) is irrelevant for the
   *          shadow matrix.
   */
  ESL_ALLOC(mx->yshadow,  sizeof(char **) * M);
  ESL_ALLOC(mx->kshadow,  sizeof(int **)  * M);
 
  /* level 3: matrix cell memory, when creating only allocate 1 cell per state, for j = 0, d = 0 */
  nbifs = CMCountStatetype(cm, B_st);

  ESL_ALLOC(mx->nrowsA, sizeof(int)      * M);
  ESL_ALLOC(mx->yshadow_mem, (sizeof(char) * (M-nbifs)));
  ESL_ALLOC(mx->kshadow_mem, (sizeof(int)  * (nbifs)));

  nb = 0;
  for (v = 0; v < M; v++) {
    if(cm->sttype[v] == B_st) { 
      ESL_ALLOC(mx->kshadow[v], sizeof(int *) * (allocL));
      mx->kshadow[v][0] = mx->kshadow_mem + nb * (allocL) * (allocW);
      mx->yshadow[v] = NULL;
      nb++;
    }
    else { 
      ESL_ALLOC(mx->yshadow[v], sizeof(char *) * (allocL));
      mx->yshadow[v][0] = mx->yshadow_mem + (v-nb) * (allocL) * (allocW);
      mx->kshadow[v] = NULL;
    }
    mx->nrowsA[v] = allocL;
  }
  mx->M            = M;
  mx->nbifs        = nbifs;
  mx->y_ncells_alloc = (M-nbifs)*(allocL)*(allocW);
  mx->y_ncells_valid = 0;
  mx->k_ncells_alloc = (nbifs)*(allocL)*(allocW);
  mx->k_ncells_valid = 0;
  mx->L            = allocL; /* allocL = 1 */

  mx->size_Mb = 
    (float) (mx->M) * (float) sizeof(int *) +    /* nrowsA ptrs */
    (float) (mx->M - nbifs) * (float) sizeof(char **) + /* mx->yshadow[] ptrs */
    (float) (mx->M - nbifs) * (float) sizeof(char *) +  /* mx->yshadow[v][] ptrs */
    (float) (nbifs)         * (float) sizeof(int **) + /* mx->yshadow[] ptrs */
    (float) (nbifs)         * (float) sizeof(int *) +  /* mx->yshadow[v][] ptrs */
    (float) mx->y_ncells_alloc * (float) sizeof(char) + /* mx->yshadow_mem */
    (float) mx->k_ncells_alloc * (float) sizeof(int);  /* mx->kshadow_mem */
  mx->size_Mb *= 0.000001; /* convert to Mb */

  return mx;

 ERROR:
  if (mx != NULL) cm_hb_shadow_mx_Destroy(mx);
  return NULL;
}

/* Function:  cm_hb_shadow_mx_GrowTo()
 * Incept:    EPN, Fri Oct 26 05:19:49 2007
 *
 * Purpose:   Assures that a shadow matrix <mx> is allocated
 *            for a model of exactly <mx->M> states and required number of 
 *            total cells. Determines new required size from 
 *            the CP9Bands_t object passed in, and reallocates if 
 *            necessary.
 *            
 *            Checks that the matrix has been created for the current CM.
 *            Check is that  mx->yshadow[v] == NULL when v is a B_st and
 *                           mx->kshadow[v] != NULL when v is a B_st.
 *
 *            DOES NOT check size of desired matrix, because an
 *            actual DP matrix should be much bigger (dp mx is floats,
 *            whereas a shadow mx is mostly chars with some ints).
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
 * Throws:    <eslEINCOMPAT> if mx does not appeared to be created for this cm
 *                           This should be caught and appropriately handled by caller. 
 *            <eslEINCOMPAT> on contract violation
 *            <eslEMEM> on memory allocation error.
 */
int
cm_hb_shadow_mx_GrowTo(CM_t *cm, CM_HB_SHADOW_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L)
{
  int     status;
  void   *p;
  int     v, jp;
  int     y_cur_size, k_cur_size = 0;
  size_t  y_ncells, k_ncells;
  int     jbw;
  double  Mb_needed;
  /* contract check, number of states (M) is something we don't change
   * so check this matrix has same number of 1st dim state ptrs that
   * cp9b has */
  if(cp9b == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_shadow_mx_GrowTo() entered with cp9b == NULL.\n");
  if(cp9b->cm_M != mx->M) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_shadow_mx_GrowTo() entered with mx->M: (%d) != cp9b->M (%d)\n", mx->M, cp9b->cm_M);
  
  y_ncells = k_ncells = 0;
  Mb_needed = ((float) (sizeof(int *)) * ((float) mx->M + 1)) + /* nrowsA ptrs */
    (float) (sizeof(char **)) * (float) (mx->M-mx->nbifs) +     /* mx->yshadow[] ptrs */
    (float) (sizeof(int  **)) * (float) (mx->nbifs);            /* mx->kshadow[] ptrs */
  for(v = 0; v < mx->M; v++) { 
    jbw = cp9b->jmax[v] - cp9b->jmin[v]; 
    if(cm->sttype[v] == B_st) { 
      if(mx->kshadow[v] == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_shadow_mx_GrowTo() v is a B st %d, but mx->kshadow[v] == NULL.\n", v);
      if(mx->yshadow[v] != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_shadow_mx_GrowTo() v is a B st %d, but mx->yshadow[v] != NULL.\n", v);
      Mb_needed += (float) (sizeof(int *) * (jbw+1)); /* mx->kshadow[v][] ptrs */
      for(jp = 0; jp <= jbw; jp++) 
	k_ncells += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
    }
    else { 
      if(mx->kshadow[v] != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_shadow_mx_GrowTo() v is not a B st %d, but mx->kshadow[v] != NULL.\n", v);
      if(mx->yshadow[v] == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_shadow_mx_GrowTo() v is not a B st %d, but mx->kshadow[v] == NULL.\n", v);
      Mb_needed += (float) (sizeof(char *) * (jbw+1)); /* mx->yshadow[v][] ptrs */
      for(jp = 0; jp <= jbw; jp++) 
	y_ncells += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
    }
  }
  Mb_needed += ESL_MAX(((float) (sizeof(char) * mx->y_ncells_alloc)), ((float) (sizeof(char) * y_ncells))); /* mx->yshadow_mem */
  Mb_needed += ESL_MAX(((float) (sizeof(char) * mx->k_ncells_alloc)), ((float) (sizeof(char) * k_ncells))); /* mx->kshadow_mem */
  Mb_needed *= 0.000001; /* convert to megabytes */

  /* must we realloc the full yshadow and kshadow matrices? 
   * or can we get away with just jiggering the pointers, if 
   * total required num cells is less than or equal to what 
   * we already have alloc'ed?
   */

  /* handle yshadow */
  if (y_ncells > mx->y_ncells_alloc) {
      ESL_RALLOC(mx->yshadow_mem, p, sizeof(char) * y_ncells);
      mx->y_ncells_alloc = y_ncells;
  }
  mx->y_ncells_valid = y_ncells;

  /* handle kshadow */
  if (k_ncells > mx->k_ncells_alloc) {
      ESL_RALLOC(mx->kshadow_mem, p, sizeof(int) * k_ncells);
      mx->k_ncells_alloc = k_ncells;
  }
  mx->k_ncells_valid = k_ncells;

  /* make sure each row is big enough */
  for(v = 0; v < mx->M; v++) {
    jbw = cp9b->jmax[v] - cp9b->jmin[v] + 1;
    if(jbw > mx->nrowsA[v]) {
      if(cm->sttype[v] == B_st) { 
	ESL_RALLOC(mx->kshadow[v], p, sizeof(int *) * jbw);
	mx->nrowsA[v] = jbw;
      }
      else { 
	ESL_RALLOC(mx->yshadow[v], p, sizeof(char *) * jbw);
	mx->nrowsA[v] = jbw;
      }
    }
  }

  /* reset the pointers, we keep a tally of number of cells
   * we've seen in each matrix (y_cur_size and k_cur_size) as we go,
   * we could precalc it and store it for each v,j, but that 
   * would be wasteful, as we'll only use the matrix configured
   * this way once, in a banded CYK run.
   */
  y_cur_size = k_cur_size = 0;
  for(v = 0; v < mx->M; v++) { 
    if(cm->sttype[v] == B_st) { 
      for(jp = 0; jp <= (cp9b->jmax[v] - cp9b->jmin[v]); jp++) { 
	mx->kshadow[v][jp] = mx->kshadow_mem + k_cur_size;
	k_cur_size        += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
      }
    }
    else { 
      for(jp = 0; jp <= (cp9b->jmax[v] - cp9b->jmin[v]); jp++) { 
	mx->yshadow[v][jp] = mx->yshadow_mem + y_cur_size;
	y_cur_size        += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
      }
    }
  }
  /* TEMP */
  assert(mx->y_ncells_valid == y_cur_size);
  assert(mx->k_ncells_valid == k_cur_size);
  ESL_DASSERT1((y_cur_size == mx->y_ncells_valid));
  ESL_DASSERT1((k_cur_size == mx->k_ncells_valid));

  mx->cp9b = cp9b; /* just a reference */
  
  /* now update L and size_Mb */
  mx->L       = L;    /* length of current seq we're valid for */
  mx->size_Mb = Mb_needed;
  
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_hb_shadow_mx_Destroy()
 * Synopsis:  Frees a DP matrix.
 * Incept:    EPN, Fri Oct 26 09:04:04 2007
 *
 * Purpose:   Frees a <CM_HB_SHADOW_MX>.
 *
 * Returns:   (void)
 */
void
cm_hb_shadow_mx_Destroy(CM_HB_SHADOW_MX *mx)
{
  if (mx == NULL) return;
  int v;

  if (mx->yshadow      != NULL) { 
    for (v = 0; v < mx->M; v++) 
      if(mx->yshadow[v] != NULL) free(mx->yshadow[v]);  
  }
  free(mx->yshadow);

  if (mx->kshadow      != NULL) { 
    for (v = 0; v < mx->M; v++) 
      if(mx->kshadow[v] != NULL) free(mx->kshadow[v]);  
  }
  free(mx->kshadow);

  if (mx->nrowsA  != NULL)       free(mx->nrowsA);
  if (mx->yshadow_mem  != NULL)  free(mx->yshadow_mem);
  if (mx->kshadow_mem  != NULL)  free(mx->kshadow_mem);
  free(mx);
  return;
}

/* Function:  cm_hb_shadow_mx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 * Incept:    EPN, Fri Oct 26 09:04:46 2007
 *
 * Purpose:   Dump matrix <mx> to stream <fp> for diagnostics.
 */
int
cm_hb_shadow_mx_Dump(FILE *ofp, CM_t *cm, CM_HB_SHADOW_MX *mx)
{
  int v, jp, j, dp, d;

  fprintf(ofp, "M: %d\nnbifs: %d\nL: %d\ny_ncells_alloc: %d\ny_ncells_valid: %d\nk_ncells_alloc: %d\nk_ncells_valid: %d\n", mx->M, mx->nbifs, mx->L, mx->y_ncells_alloc, mx->y_ncells_valid, mx->k_ncells_alloc, mx->k_ncells_valid);
  
  /* yshadow/kshadow matrix data */
  for (v = 0; v < mx->M; v++) {
    if(cm->sttype[v] == B_st) { 
      for(jp = 0; jp <= mx->cp9b->jmax[v] - mx->cp9b->jmin[v]; jp++) {
	j = jp + mx->cp9b->jmin[v];
	for(dp = 0; dp <= mx->cp9b->hdmax[v][jp] - mx->cp9b->hdmin[v][jp]; dp++) {
	  d = dp + mx->cp9b->hdmin[v][jp];
	  fprintf(ofp, "kshad[v:%5d][j:%5d][d:%5d] %8d\n", v, j, d, mx->kshadow[v][jp][dp]);
	}
	fprintf(ofp, "\n");
      }
      fprintf(ofp, "\n\n");
    }
    else { 
      for(jp = 0; jp <= mx->cp9b->jmax[v] - mx->cp9b->jmin[v]; jp++) {
	j = jp + mx->cp9b->jmin[v];
	for(dp = 0; dp <= mx->cp9b->hdmax[v][jp] - mx->cp9b->hdmin[v][jp]; dp++) {
	  d = dp + mx->cp9b->hdmin[v][jp];
	  fprintf(ofp, "kshad[v:%5d][j:%5d][d:%5d] %8c\n", v, j, d, mx->yshadow[v][jp][dp]);
	}
	fprintf(ofp, "\n");
      }
      fprintf(ofp, "\n\n");
    }
  }
  return eslOK;
}

/*****************************************************************
 *   3. ScanMatrix_t data structure functions,
 *      auxiliary info and matrix of float and/or int scores for 
 *      query dependent banded or non-banded CM DP search functions
 *****************************************************************/

/* Function: cm_CreateScanMatrix()
 * Date:     EPN, Sun Nov  4 19:56:58 2007
 *
 * Purpose:  Given relevant info, allocate and initialize ScanMatrix_t object.
 *            
 * Returns:  eslOK on success, dies immediately on some error
 */
ScanMatrix_t *
cm_CreateScanMatrix(CM_t *cm, int W, int *dmin, int *dmax, double beta_W, double beta_qdb, int do_banded, int do_float, int do_int)
{ 
  int status;
  ScanMatrix_t *smx;
  int v,j;

  if((!do_float) && (!do_int)) cm_Fail("cm_CreateScanMatrix(), do_float and do_int both FALSE.\n");
  if(do_banded && (dmin == NULL || dmax == NULL)) cm_Fail("cm_CreateScanMatrix(), do_banded is TRUE, but dmin or dmax is NULL.\n");

  ESL_ALLOC(smx, sizeof(ScanMatrix_t));

  smx->flags    = 0;
  smx->cm_M     = cm->M;
  smx->W        = W;
  smx->dmin     = dmin; /* could be NULL */
  smx->dmax     = dmax; /* could be NULL */
  smx->beta_W   = beta_W; 
  smx->beta_qdb = beta_qdb; 

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

  smx->falpha          = NULL;
  smx->falpha_begl     = NULL;
  smx->falpha_mem      = NULL;
  smx->falpha_begl_mem = NULL;

  smx->ialpha          = NULL;
  smx->ialpha_begl     = NULL;
  smx->ialpha_mem      = NULL;
  smx->ialpha_begl_mem = NULL;

  smx->ncells_alpha      = 0;
  smx->ncells_alpha_begl = 0;

  if(do_float) /* allocate float mx and scores */
    cm_FloatizeScanMatrix(cm, smx);
  if(do_int)   /* allocate int mx and scores */
    cm_IntizeScanMatrix(cm, smx);
  return smx;

 ERROR:
  cm_Fail("memory allocation error in cm_CreateScanMatrix().\n");
  return NULL; /* NEVERREACHED */
}

/* Function: cm_CreateScanMatrixForCM()
 * Date:     EPN, Fri Nov 30 06:07:23 2007
 *
 * Purpose:  Given a CM, allocate and initialize ScanMatrix_t object for that CM. 
 *           Most of work is done by cm_CreateScanMatrix(). 
 *
 * Returns:  eslOK on success, dies immediately on some error
 */
int
cm_CreateScanMatrixForCM(CM_t *cm, int do_float, int do_int)
{
  int do_banded;
  double beta_W;

  if(cm->flags & CMH_SCANMATRIX)               cm_Fail("cm_CreateScanMatrixForCM(), the CM flag for valid scan info is already up.");
  if(cm->smx != NULL)                          cm_Fail("cm_CreateScanMatrixForCM(), the cm already points to a ScanMatrix_t object.\n");
  if(cm->dmin == NULL && cm->dmax != NULL)     cm_Fail("cm_CreateScanMatrixForCM(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmin != NULL && cm->dmax == NULL)     cm_Fail("cm_CreateScanMatrixForCM(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmax != NULL && cm->W > cm->dmax[0])  cm_Fail("cm_CreateScanMatrixForCM(), cm->W: %d > cm->dmax[0]: %d (cm->W can be less than cm->dmax[0] but not greater)\n", cm->W, cm->dmax[0]); 
  if((cm->dmin != NULL && cm->dmax != NULL) && (! (cm->flags & CMH_QDB))) 
     cm_Fail("cm_CreateScanMatrixForCM(), cm->dmin != NULL, cm->dmax != NULL, but CMH_QDB flag down, bands are invalid\n"); 
  if((! cm->search_opts & CM_SEARCH_NOQDB) && (cm->dmin == NULL || cm->dmax == NULL))
    cm_Fail("cm_CreateScanMatrixForCM(), cm->dmin == NULL || cm->dmax == NULL, but !(cm->search_opts & CM_SEARCH_NOQDB)\n");

  do_banded   = (cm->search_opts & CM_SEARCH_NOQDB) ? FALSE : TRUE;
  
  if(do_banded) beta_W = ESL_MAX(cm->beta_W, cm->beta_qdb);
  else beta_W = cm->beta_W;
  /*printf("TEMP cm_CreateScanMatrix(), do_banded: %d beta_W: %g beta_qdb: %g\n", do_banded, beta_W, cm->beta_qdb);*/
  cm->smx = cm_CreateScanMatrix(cm, cm->W, cm->dmin, cm->dmax, beta_W, cm->beta_qdb, do_banded, do_float, do_int);
  cm->flags |= CMH_SCANMATRIX; /* raise the flag for valid CMH_SCANMATRIX */
  if(cm->si != NULL) { /* CM has searchinfo, update the final round so it points at cm->smx */
    UpdateSearchInfoForNewSMX(cm);
  }

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
cm_FloatizeScanMatrix(CM_t *cm, ScanMatrix_t *smx)
{
  int status;
  int j, v;
  int y, yoffset, w;
  int use_hmmonly;
  use_hmmonly = ((cm->search_opts & CM_SEARCH_HMMVITERBI) ||  (cm->search_opts & CM_SEARCH_HMMFORWARD)) ? TRUE : FALSE;
  int n_begl;
  int n_non_begl;
  int cur_cell;

  /* contract check */
  if(smx->flags & cmSMX_HAS_FLOAT) cm_Fail("cm_FloatizeScanMatrix(), si's cmSMX_HAS_FLOAT flag is already up.");
  if(smx->falpha != NULL)       cm_Fail("cm_FloatizeScanMatrix(), smx->falpha is not NULL.");
  if(smx->falpha_begl != NULL)  cm_Fail("cm_FloatizeScanMatrix(), smx->falpha_begl is not NULL.");
  
  /* allocate alpha 
   * we allocate only as many cells as necessary,
   * for falpha,      we only allocate for non-BEGL_S states,
   * for falpha_begl, we only allocate for     BEGL_S states
   *
   * note: deck for the EL state, cm->M is never used for scanners
   */
  n_begl = 0;
  for (v = 0; v < cm->M; v++) if (cm->stid[v] == BEGL_S) n_begl++;
  n_non_begl = cm->M - n_begl;

  /* allocate falpha */
  /* j == 0 v == 0 cells, followed by j == 1 v == 0, then j == 0 v == 1 etc.. */
  ESL_ALLOC(smx->falpha,        sizeof(float **) * 2);
  ESL_ALLOC(smx->falpha[0],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, falpha[0][v] will be NULL */
  ESL_ALLOC(smx->falpha[1],     sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, falpha[0][v] will be NULL */
  ESL_ALLOC(smx->falpha_mem,    sizeof(float) * 2 * n_non_begl * (smx->W+1));
  if((smx->flags & cmSMX_HAS_INT) && ((2 * n_non_begl * (smx->W+1)) != smx->ncells_alpha)) 
    cm_Fail("cm_FloatizeScanMatrix(), cmSMX_HAS_INT flag raised, but smx->ncells_alpha %d != %d (predicted num float cells size)\n", smx->ncells_alpha, (2 * n_non_begl * (smx->W+1)));
  smx->ncells_alpha = 2 * n_non_begl * (smx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    if (cm->stid[v] != BEGL_S) {
      smx->falpha[0][v] = smx->falpha_mem + cur_cell;
      cur_cell += smx->W+1;
      smx->falpha[1][v] = smx->falpha_mem + cur_cell;
      cur_cell += smx->W+1;
    }
    else { 
      smx->falpha[0][v] = NULL;
      smx->falpha[1][v] = NULL;
    }
  }
  if(cur_cell != smx->ncells_alpha) cm_Fail("cm_FloatizeScanMatrix(), error allocating falpha, cell cts differ %d != %d\n", cur_cell, smx->ncells_alpha);

  /* allocate falpha_begl */
  /* j == d, v == 0 cells, followed by j == d+1, v == 0, etc. */
  ESL_ALLOC(smx->falpha_begl, sizeof(float **) * (smx->W+1));
  for (j = 0; j <= smx->W; j++) 
    ESL_ALLOC(smx->falpha_begl[j],  sizeof(float *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, falpha_begl[0][v] will be NULL */
  ESL_ALLOC(smx->falpha_begl_mem,   sizeof(float) * (smx->W+1) * n_begl * (smx->W+1));
  if((smx->flags & cmSMX_HAS_INT) && (((smx->W+1) * n_begl * (smx->W+1)) != smx->ncells_alpha_begl)) 
    cm_Fail("cm_IntizeScanMatrix(), cmSMX_HAS_FLOAT flag raised, but smx->ncells_alpha_begl %d != %d (predicted num float cells size)\n", smx->ncells_alpha_begl, ((smx->W+1) * n_begl * (smx->W+1)));
  smx->ncells_alpha_begl = (smx->W+1) * n_begl * (smx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    for (j = 0; j <= smx->W; j++) { 
      if (cm->stid[v] == BEGL_S) {
	smx->falpha_begl[j][v] = smx->falpha_begl_mem + cur_cell;
	cur_cell += smx->W+1;
      }
      else smx->falpha_begl[j][v] = NULL;
    }
  }
  if(cur_cell != smx->ncells_alpha_begl) cm_Fail("cm_FloatizeScanMatrix(), error allocating falpha_begl, cell cts differ %d != %d\n", cur_cell, smx->ncells_alpha_begl);

  /* Initialize matrix */
  /* First, init entire matrix to IMPOSSIBLE */
  esl_vec_FSet(smx->falpha_mem,      smx->ncells_alpha,      IMPOSSIBLE);
  esl_vec_FSet(smx->falpha_begl_mem, smx->ncells_alpha_begl, IMPOSSIBLE);
  /* Now, initialize cells that should not be IMPOSSIBLE in falpha and falpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      if (cm->sttype[v] == E_st) { 
	smx->falpha[0][v][0] = smx->falpha[1][v][0] = 0.;
	/* rest of E deck is IMPOSSIBLE, it's already set */
      }
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
cm_IntizeScanMatrix(CM_t *cm, ScanMatrix_t *smx)
{
  int status;
  int v, j, y, yoffset, w;
  int use_hmmonly;
  use_hmmonly = ((cm->search_opts & CM_SEARCH_HMMVITERBI) ||  (cm->search_opts & CM_SEARCH_HMMFORWARD)) ? TRUE : FALSE;
  int n_begl;
  int n_non_begl;
  int cur_cell;

  /* contract check */
  if(smx->flags & cmSMX_HAS_INT) cm_Fail("cm_IntizeScanMatrix(), si's cmSMX_HAS_INT flag is already up.");
  if(smx->ialpha != NULL)       cm_Fail("cm_IntizeScanMatrix(), smx->ialpha is not NULL.");
  if(smx->ialpha_begl != NULL)  cm_Fail("cm_IntizeScanMatrix(), smx->ialpha_begl is not NULL.");

  /* allocate alpha 
   * we allocate only as many cells as necessary,
   * for ialpha,      we only allocate for non-BEGL_S states,
   * for ialpha_begl, we only allocate for     BEGL_S states
   *
   * note: deck for the EL state, cm->M is never used for scanners
   */
  n_begl = 0;
  for (v = 0; v < cm->M; v++) if (cm->stid[v] == BEGL_S) n_begl++;
  n_non_begl = cm->M - n_begl;

  /* allocate ialpha */
  /* j == 0 v == 0 cells, followed by j == 1 v == 0, then j == 0 v == 1 etc.. */
  ESL_ALLOC(smx->ialpha,        sizeof(int **) * 2);
  ESL_ALLOC(smx->ialpha[0],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, ialpha[0][v] will be NULL */
  ESL_ALLOC(smx->ialpha[1],     sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v == BEGL_S, ialpha[0][v] will be NULL */
  ESL_ALLOC(smx->ialpha_mem,    sizeof(int) * 2 * n_non_begl * (smx->W+1));
  if((smx->flags & cmSMX_HAS_FLOAT) && ((2 * n_non_begl * (smx->W+1)) != smx->ncells_alpha)) 
    cm_Fail("cm_IntizeScanMatrix(), cmSMX_HAS_INT flag raised, but smx->ncells_alpha %d != %d (predicted num int cells size)\n", smx->ncells_alpha, (2 * n_non_begl * (smx->W+1)));
  smx->ncells_alpha = 2 * n_non_begl * (smx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    if (cm->stid[v] != BEGL_S) {
      smx->ialpha[0][v] = smx->ialpha_mem + cur_cell;
      cur_cell += smx->W+1;
      smx->ialpha[1][v] = smx->ialpha_mem + cur_cell;
      cur_cell += smx->W+1;
    }
    else { 
      smx->ialpha[0][v] = NULL;
      smx->ialpha[1][v] = NULL;
    }
  }
  if(cur_cell != smx->ncells_alpha) cm_Fail("cm_IntizeScanMatrix(), error allocating ialpha, cell cts differ %d != %d\n", cur_cell, smx->ncells_alpha);
  
  /* allocate ialpha_begl */
  /* j == d, v == 0 cells, followed by j == d+1, v == 0, etc. */
  ESL_ALLOC(smx->ialpha_begl, sizeof(int **) * (smx->W+1));
  for (j = 0; j <= smx->W; j++) 
    ESL_ALLOC(smx->ialpha_begl[j],  sizeof(int *) * (cm->M)); /* we still allocate cm->M ptrs, if v != BEGL_S, ialpha_begl[0][v] will be NULL */
  ESL_ALLOC(smx->ialpha_begl_mem,   sizeof(int) * (smx->W+1) * n_begl * (smx->W+1));
  if((smx->flags & cmSMX_HAS_FLOAT) && (((smx->W+1) * n_begl * (smx->W+1)) != smx->ncells_alpha_begl)) 
    cm_Fail("cm_IntizeScanMatrix(), cmSMX_HAS_INT flag raised, but smx->ncells_alpha_begl %d != %d (predicted num int cells size)\n", smx->ncells_alpha_begl, ((smx->W+1) * n_begl * (smx->W+1)));
  smx->ncells_alpha_begl = (smx->W+1) * n_begl * (smx->W+1);

  cur_cell = 0;
  for (v = 0; v < cm->M; v++) {	
    for (j = 0; j <= smx->W; j++) { 
      if (cm->stid[v] == BEGL_S) {
	smx->ialpha_begl[j][v] = smx->ialpha_begl_mem + cur_cell;
	cur_cell += smx->W+1;
      }
      else smx->ialpha_begl[j][v] = NULL;
    }
  }
  if(cur_cell != smx->ncells_alpha_begl) cm_Fail("cm_IntizeScanMatrix(), error allocating ialpha_begl, cell cts differ %d != %d\n", cur_cell, smx->ncells_alpha_begl);

  /* Initialize matrix */
  /* First, init entire matrix to -INFTY */
  esl_vec_ISet(smx->ialpha_mem,      smx->ncells_alpha,      -INFTY);
  esl_vec_ISet(smx->ialpha_begl_mem, smx->ncells_alpha_begl, -INFTY);
  /* Now, initialize cells that should not be -INFTY in ialpha and ialpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      if (cm->sttype[v] == E_st) { 
	smx->ialpha[0][v][0] = smx->ialpha[1][v][0] = 0.;
	/* rest of E deck is -INFTY, it's already set */
      }
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
 * Purpose:  Free float data structures in a ScanMatrix_t object 
 *           corresponding to <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeFloatsFromScanMatrix(CM_t *cm, ScanMatrix_t *smx)
{
  int j;

  /* contract check */
  if(! smx->flags & cmSMX_HAS_FLOAT)    cm_Fail("cm_FreeFloatsFromScanMatrix(), si's cmSMX_HAS_FLOAT flag is down.");
  if(smx->falpha == NULL)       cm_Fail("cm_FreeFloatsFromScanMatrix(), smx->falpha is already NULL.");
  if(smx->falpha_begl == NULL)  cm_Fail("cm_FreeFloatsFromScanMatrix(), smx->falpha_begl is already NULL.");

  free(smx->falpha_mem);
  free(smx->falpha[1]);
  free(smx->falpha[0]);
  free(smx->falpha);
  smx->falpha = NULL;

  free(smx->falpha_begl_mem);
  for (j = 0; j <= smx->W; j++) free(smx->falpha_begl[j]);
  free(smx->falpha_begl);
  smx->falpha_begl = NULL;

  smx->flags &= ~cmSMX_HAS_FLOAT;
  return eslOK;
}

/* Function: cm_FreeIntsFromScanMatrix()
 * Date:     EPN, Wed Nov  7 09:56:01 2007
 *
 * Purpose:  Free int data structures in a ScanMatrix_t object 
 *           corresponding to <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeIntsFromScanMatrix(CM_t *cm, ScanMatrix_t *smx)
{
  int j;

  /* contract check */
  if(! smx->flags & cmSMX_HAS_INT)    cm_Fail("cm_FreeIntsFromScanMatrix(), si's cmSMX_HAS_INT flag is down.");
  if(smx->ialpha == NULL)       cm_Fail("cm_FreeIntsFromScanMatrix(), smx->ialpha is already NULL.");
  if(smx->ialpha_begl == NULL)  cm_Fail("cm_FreeIntsFromScanMatrix(), smx->ialpha_begl is already NULL.");

  free(smx->ialpha_mem);
  free(smx->ialpha[1]);
  free(smx->ialpha[0]);
  free(smx->ialpha);
  smx->ialpha = NULL;

  free(smx->ialpha_begl_mem);
  for (j = 0; j <= smx->W; j++) free(smx->ialpha_begl[j]);
  free(smx->ialpha_begl);
  smx->ialpha_begl = NULL;

  smx->flags &= ~cmSMX_HAS_INT;
  return eslOK;
}

/* Function: cm_FreeScanMatrix()
 * Date:     EPN, Sun Nov  4 20:57:32 2007
 *
 * Purpose:  Free a ScanMatrix_t object corresponding
 *           to CM <cm>.
 *            
 * Returns:  void
 */
void
cm_FreeScanMatrix(CM_t *cm, ScanMatrix_t *smx)
{
  int j;
  if(! ((cm->flags & CMH_SCANMATRIX) && (smx == cm->smx))) { /* don't free the cm->smx's dmin, dmax */
    if(smx->dmin != cm->dmin && smx->dmin != NULL) { free(smx->dmin); smx->dmin = NULL; }
    if(smx->dmax != cm->dmax && smx->dmax != NULL) { free(smx->dmax); smx->dmax = NULL; }
  }

  for(j = 1; j <= smx->W; j++) {
    free(smx->dnAA[j]);
    free(smx->dxAA[j]);
  }
  free(smx->dnAA);
  free(smx->dxAA);
  free(smx->bestr);
  
  if(smx->flags & cmSMX_HAS_FLOAT) cm_FreeFloatsFromScanMatrix(cm, smx);
  if(smx->flags & cmSMX_HAS_INT)   cm_FreeIntsFromScanMatrix(cm, smx);
  free(smx);
  return;
}


/* Function: cm_FreeScanMatrixForCM()
 * Date:     EPN, Fri Nov 30 06:47:12 2007
 *
 * Purpose:  Free a ScanMatrix_t object cm->smx for <cm>.
 *            
 * Returns:  void
 */
void
cm_FreeScanMatrixForCM(CM_t *cm)
{
  /* contract check */
  if(cm->smx == NULL) cm_Fail("cm_FreeScanMatrixForCM(), cm->smx is NULL.\n");
  cm_FreeScanMatrix(cm, cm->smx);
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
  ESL_DSQ a,b;
  int v;
  int cur_cell;
  int npairs = 0;
  int nsinglets = 0;
  float *ptr_to_start; /* points to block allocated to esc_vAA[0], nec b/c esc_vAA[0] gets set to NULL, because v == 0 is non-emitter */
  float **leftAA;
  float **rightAA;

  /* count pairs, singlets */
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      nsinglets++;
      break;
    case MP_st:
      npairs++;
      break;
    }
  }

  /* set up our left and right vectors for all possible non-canonical residues,
   * these are calc'ed once and passed to FastPairScore*() functions to minimize
   * run time. 
   */
  ESL_ALLOC(leftAA,  sizeof(float *) * cm->abc->Kp);
  ESL_ALLOC(rightAA, sizeof(float *) * cm->abc->Kp);
  for(a = 0; a <= cm->abc->K; a++) leftAA[a] = rightAA[a] = NULL; /* canonicals and gap, left/right unnec */
  for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) {
    ESL_ALLOC(leftAA[a],  sizeof(float) * cm->abc->K);
    ESL_ALLOC(rightAA[a], sizeof(float) * cm->abc->K);
    esl_vec_FSet(leftAA[a],  cm->abc->K, 0.);
    esl_vec_FSet(rightAA[a], cm->abc->K, 0.);
    esl_abc_FCount(cm->abc, leftAA[a],  a, 1.);
    esl_abc_FCount(cm->abc, rightAA[a], a, 1.);
  }
  leftAA[cm->abc->Kp-1] = rightAA[cm->abc->Kp-1] = NULL; /* missing data, left/right unnec */

  /* precalculate possible emission scores for each state */
  ESL_ALLOC(esc_vAA,     sizeof(float *) * (cm->M));
  ESL_ALLOC(esc_vAA[0],  sizeof(float)   * ((cm->abc->Kp * nsinglets) + (cm->abc->Kp * cm->abc->Kp * npairs)));
  ptr_to_start = esc_vAA[0];
  cur_cell = 0;
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      esc_vAA[v] = ptr_to_start + cur_cell;
      cur_cell += cm->abc->Kp;
      for(a = 0; a < cm->abc->K; a++) /* all canonical residues */
	esc_vAA[v][a]  = cm->esc[v][a]; 
      esc_vAA[v][cm->abc->K] = IMPOSSIBLE; /* gap symbol is impossible */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) /* all ambiguous residues */
	esc_vAA[v][a]  = esl_abc_FAvgScore(cm->abc, a, cm->esc[v]);
      esc_vAA[v][cm->abc->Kp-1] = IMPOSSIBLE; /* missing data is IMPOSSIBLE */
      break;
    case MP_st:
      esc_vAA[v] = ptr_to_start + cur_cell;
      esl_vec_FSet(esc_vAA[v], cm->abc->Kp * cm->abc->Kp, IMPOSSIBLE); /* init all cells to IMPOSSIBLE */
      cur_cell += cm->abc->Kp * cm->abc->Kp;
      /* a is canonical, b is canonical */
      for(a = 0; a < cm->abc->K; a++) { 
	for(b = 0; b < cm->abc->K; b++) { 
	  esc_vAA[v][(a * cm->abc->Kp) + b]  = cm->esc[v][(a * cm->abc->K) + b];
	}
      }
      /* a is not canonical, b is canonical */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) { 
	for(b = 0; b < cm->abc->K; b++) { 
	  esc_vAA[v][(a * cm->abc->Kp) + b]  = FastPairScoreLeftOnlyDegenerate(cm->abc->K, cm->esc[v], leftAA[a], b);
	}
      }	  
      /* a is canonical, b is not canonical */
      for(a = 0; a < cm->abc->K; a++) { 
	for(b = cm->abc->K+1; b < cm->abc->Kp-1; b++) { 
	  esc_vAA[v][(a * cm->abc->Kp) + b]  = FastPairScoreRightOnlyDegenerate(cm->abc->K, cm->esc[v], rightAA[b], a);
	}
      }	  
      /* a is not canonical, b is not canonical */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) { 
	for(b = cm->abc->K+1; b < cm->abc->Kp-1; b++) { 
	  esc_vAA[v][(a * cm->abc->Kp) + b]  = FastPairScoreBothDegenerate(cm->abc->K, cm->esc[v], leftAA[a], rightAA[b]);
	}
      }	  
      /* everything else, when either a or b is gap or missing data, stays IMPOSSIBLE */
      break;
    default:
      esc_vAA[v] = NULL;
      break;
    }
  }
  for(a = 0; a < cm->abc->Kp; a++) { 
    if(leftAA[a] != NULL)  free(leftAA[a]);
    if(rightAA[a] != NULL) free(rightAA[a]);
  }
  free(leftAA);
  free(rightAA);
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
  ESL_DSQ a,b;
  int v;
  int cur_cell;
  int npairs = 0;
  int nsinglets = 0;
  int *ptr_to_start; /* points to block allocated to iesc_vAA[0], nec b/c esc_vAA[0] gets set to NULL, because v == 0 is non-emitter */
  float **leftAA;
  float **rightAA;

  /* count pairs, singlets */
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      nsinglets++;
      break;
    case MP_st:
      npairs++;
      break;
    }
  }

  /* set up our left and right vectors for all possible non-canonical residues,
   * these are calc'ed once and passed to FastPairScore*() functions to minimize
   * run time. 
   */
  ESL_ALLOC(leftAA,  sizeof(float *) * cm->abc->Kp);
  ESL_ALLOC(rightAA, sizeof(float *) * cm->abc->Kp);
  for(a = 0; a <= cm->abc->K; a++) leftAA[a] = rightAA[a] = NULL; /* canonicals and gap, left/right unnec */
  for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) {
    ESL_ALLOC(leftAA[a],  sizeof(float) * cm->abc->K);
    ESL_ALLOC(rightAA[a], sizeof(float) * cm->abc->K);
    esl_vec_FSet(leftAA[a],  cm->abc->K, 0.);
    esl_vec_FSet(rightAA[a], cm->abc->K, 0.);
    esl_abc_FCount(cm->abc, leftAA[a],  a, 1.);
    esl_abc_FCount(cm->abc, rightAA[a], a, 1.);
  }
  leftAA[cm->abc->Kp-1] = rightAA[cm->abc->Kp-1] = NULL; /* missing data, left/right unnec */

  /* precalculate possible emission scores for each state */
  ESL_ALLOC(iesc_vAA,     sizeof(int *) * (cm->M));
  ESL_ALLOC(iesc_vAA[0],  sizeof(int)   * ((cm->abc->Kp * nsinglets) + (cm->abc->Kp * cm->abc->Kp * npairs)));
  ptr_to_start = iesc_vAA[0];
  cur_cell = 0;
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      iesc_vAA[v] = ptr_to_start + cur_cell;
      cur_cell += cm->abc->Kp;
      for(a = 0; a < cm->abc->K; a++) /* all canonical residues */
	iesc_vAA[v][a]  = cm->iesc[v][a]; 
      iesc_vAA[v][cm->abc->K] = -INFTY; /* gap symbol is impossible */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) /* all ambiguous residues */
	iesc_vAA[v][a]  = esl_abc_IAvgScore(cm->abc, a, cm->iesc[v]);
      iesc_vAA[v][cm->abc->Kp-1] = -INFTY; /* missing data is IMPOSSIBLE */
      break;
    case MP_st:
      iesc_vAA[v] = ptr_to_start + cur_cell;
      esl_vec_ISet(iesc_vAA[v], cm->abc->Kp * cm->abc->Kp, -INFTY); /* init all cells to -INFTY */
      cur_cell += cm->abc->Kp * cm->abc->Kp;
      /* a is canonical, b is canonical */
      for(a = 0; a < cm->abc->K; a++) { 
	for(b = 0; b < cm->abc->K; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = cm->iesc[v][(a * cm->abc->K) + b];
	}
      }
      /* a is not canonical, b is canonical */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) { 
	for(b = 0; b < cm->abc->K; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = iFastPairScoreLeftOnlyDegenerate(cm->abc->K, cm->iesc[v], leftAA[a], b);
	}
      }	  
      /* a is canonical, b is not canonical */
      for(a = 0; a < cm->abc->K; a++) { 
	for(b = cm->abc->K+1; b < cm->abc->Kp-1; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = iFastPairScoreRightOnlyDegenerate(cm->abc->K, cm->iesc[v], rightAA[b], a);
	}
      }	  
      /* a is not canonical, b is not canonical */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) { 
	for(b = cm->abc->K+1; b < cm->abc->Kp-1; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = iFastPairScoreBothDegenerate(cm->abc->K, cm->iesc[v], leftAA[a], rightAA[b]);
	}
      }	  
      /* everything else, when either a or b is gap or missing data, stays -INFTY */
      break;
    default:
      iesc_vAA[v] = NULL;
      break;
    }
  }
  for(a = 0; a < cm->abc->Kp; a++) { 
    if(leftAA[a] != NULL)  free(leftAA[a]);
    if(rightAA[a] != NULL) free(rightAA[a]);
  }
  free(leftAA);
  free(rightAA);
  return iesc_vAA;

 ERROR:
  cm_Fail("memory allocation error.");
  return NULL; /* NEVERREACHED */
}


/* Function: ICopyOptimizedEmitScoresFromFloats()
 * Date:     EPN, Wed Aug 20 14:47:13 2008
 *
 * Purpose:  Allocate and fill an optimized emission score
 *           vector of integer scores. For degenerate
 *           residues calculating the scores is somewhat
 *           compute-intensive, so don't calc them, but copy/integerize
 *           scores from the CM's pre-calculated float
 *           optimized emission score vector.
 *           This is done only because is fast.
 *            
 * Returns:  the 2D integer emission score vector on success,
 *           dies immediately on memory allocation error.
 */
int **
ICopyOptimizedEmitScoresFromFloats(CM_t *cm, float **oesc)
{
  int status; 
  int **iesc_vAA;
  ESL_DSQ a,b;
  int v;
  int cur_cell;
  int npairs = 0;
  int nsinglets = 0;
  int *ptr_to_start; /* points to block allocated to iesc_vAA[0], nec b/c esc_vAA[0] gets set to NULL, because v == 0 is non-emitter */

  /* count pairs, singlets */
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      nsinglets++;
      break;
    case MP_st:
      npairs++;
      break;
    }
  }

  /* fill emission scores for each state */
  ESL_ALLOC(iesc_vAA,     sizeof(int *) * (cm->M));
  ESL_ALLOC(iesc_vAA[0],  sizeof(int)   * ((cm->abc->Kp * nsinglets) + (cm->abc->Kp * cm->abc->Kp * npairs)));
  ptr_to_start = iesc_vAA[0];
  cur_cell = 0;
  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      iesc_vAA[v] = ptr_to_start + cur_cell;
      cur_cell += cm->abc->Kp;
      for(a = 0; a < cm->abc->K; a++) /* all canonical residues */
	iesc_vAA[v][a]  = cm->iesc[v][a]; 
      iesc_vAA[v][cm->abc->K] = -INFTY; /* gap symbol is impossible */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) /* all ambiguous residues */
	iesc_vAA[v][a]  = (int) floor(0.5 + INTSCALE * oesc[v][a]); /* COPY, don't calc */
      iesc_vAA[v][cm->abc->Kp-1] = -INFTY; /* missing data is IMPOSSIBLE */
      break;
    case MP_st:
      iesc_vAA[v] = ptr_to_start + cur_cell;
      esl_vec_ISet(iesc_vAA[v], cm->abc->Kp * cm->abc->Kp, -INFTY); /* init all cells to -INFTY */
      cur_cell += cm->abc->Kp * cm->abc->Kp;
      /* a is canonical, b is canonical */
      for(a = 0; a < cm->abc->K; a++) { 
	for(b = 0; b < cm->abc->K; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = cm->iesc[v][(a * cm->abc->K) + b];
	}
      }
      /* a is not canonical, b is canonical */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) { 
	for(b = 0; b < cm->abc->K; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = (int) floor(0.5 + INTSCALE * oesc[v][(a * cm->abc->Kp)+b]); /* COPY, don't calc */
	}
      }	  
      /* a is canonical, b is not canonical */
      for(a = 0; a < cm->abc->K; a++) { 
	for(b = cm->abc->K+1; b < cm->abc->Kp-1; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = (int) floor(0.5 + INTSCALE * oesc[v][(a * cm->abc->Kp)+b]); /* COPY, don't calc */
	}
      }	  
      /* a is not canonical, b is not canonical */
      for(a = cm->abc->K+1; a < cm->abc->Kp-1; a++) { 
	for(b = cm->abc->K+1; b < cm->abc->Kp-1; b++) { 
	  iesc_vAA[v][(a * cm->abc->Kp) + b]  = (int) floor(0.5 + INTSCALE * oesc[v][(a * cm->abc->Kp)+b]); /* COPY, don't calc */
	}
      }	  
      /* everything else, when either a or b is gap or missing data, stays -INFTY */
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

/* Function: DumpOptimizedEmitScores()
 */
void
DumpOptimizedEmitScores(CM_t *cm, FILE *fp)
{
  int v;

  for(v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      if(cm->oesc  != NULL) { 
	fprintf(fp, "SF v: %5d\t", v);
	esl_vec_FDump(fp, cm->oesc[v], cm->abc->Kp, NULL);
      }
      if(cm->ioesc != NULL) { 
	fprintf(fp, "SI v: %5d\t", v);
	esl_vec_IDump(fp, cm->ioesc[v], cm->abc->Kp, NULL);
      }
      break;
    case MP_st:
      if(cm->oesc  != NULL) { 
	fprintf(fp, "PF v: %5d\t", v);
	esl_vec_FDump(fp, cm->oesc[v], cm->abc->Kp * cm->abc->Kp, NULL);
      }
      if(cm->ioesc != NULL) { 
	esl_vec_IDump(fp, cm->ioesc[v], cm->abc->Kp * cm->abc->Kp, NULL);
	fprintf(fp, "PI v: %5d\t", v);
      }
      break;
    }
  }
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
  if(fesc_vAA == NULL && iesc_vAA == NULL) cm_Fail("FreeOptimizedEmitScores() but fesc and iesc are NULL.\n");

  if(fesc_vAA != NULL) { 
    if(fesc_vAA[1] != NULL) { 
      free(fesc_vAA[1]); /* note: we free [1], but we alloc'ed to [0], why? b/c fesc_vAA[0] is set to NULL after it's
			  *       used for allocation b/c it's the ROOT_S state, a non-emitter, then fesc_vAA[1] is set
			  *       to point where it used to point (it's the ROOT_IL state, an emitter).
			  */
    }
    free(fesc_vAA);
    fesc_vAA = NULL;
  }

  if(iesc_vAA != NULL) { 
    if(iesc_vAA[1] != NULL) { 
      free(iesc_vAA[1]); /* note: we free [1], but we alloc'ed to [0], why? b/c iesc_vAA[0] is set to NULL after it's
			  *       used for allocation b/c it's the ROOT_S state, a non-emitter, then iesc_vAA[1] is set
			  *       to point where it used to point (it's the ROOT_IL state, an emitter).
			  */
    }
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
 *   4. GammaHitMx_t data structure functions,
 *      Semi HMM data structure for optimal resolution of overlapping
 *      hits for CM DP search functions.
 *****************************************************************/
  
/* Function: CreateGammaHitMx()
 * Date:     EPN, Mon Nov  5 05:22:56 2007
 *
 * Purpose:  Allocate and initialize a gamma semi-HMM for 
 *           optimal hit resolution of a CM based scan.
 *           If(do_backward), L position is init'ed instead of
 *           0th position, for Backward HMM scans.
 * 
 * Returns:  Newly allocated GammaHitMx_t object:
 */
GammaHitMx_t *
CreateGammaHitMx(int L, int i0, int be_greedy, float cutoff, int do_backward)
{
  int status;
  GammaHitMx_t *gamma;
  ESL_ALLOC(gamma, sizeof(GammaHitMx_t));

  gamma->L  = L;
  gamma->i0 = i0;
  gamma->iamgreedy = be_greedy;
  gamma->cutoff    = cutoff;
  /* allocate/initialize for CYK/Inside */
  ESL_ALLOC(gamma->mx,     sizeof(float) * (L+1));
  ESL_ALLOC(gamma->gback,  sizeof(int)   * (L+1));
  ESL_ALLOC(gamma->savesc, sizeof(float) * (L+1));
  ESL_ALLOC(gamma->saver,  sizeof(int)   * (L+1));
    
  if(do_backward) { 
    gamma->mx[L]    = 0;
    gamma->gback[L] = -1;
  } 
  else { 
    gamma->mx[0]    = 0;
    gamma->gback[0] = -1;
  }
  return gamma;

 ERROR:
  cm_Fail("memory allocation error in cm_CreateGammaHitMx().\n");
  return NULL;
}

/* Function: FreeGammaHitMx()
 * Date:     EPN, Mon Nov  5 05:32:00 2007
 *
 * Purpose:  Free a gamma semi-HMM.
 *            
 * Returns:  void;
 */
void
FreeGammaHitMx(GammaHitMx_t *gamma)
{
  free(gamma->mx);
  free(gamma->gback);
  free(gamma->savesc);
  free(gamma->saver);
  free(gamma);

  return;
}

/* Function: UpdateGammaHitMxCM()
 * Date:     EPN, Mon Nov  5 05:41:14 2007
 *
 * Purpose:  Update a gamma semi-HMM for CM hits that end at gamma-relative position <j>.
 *
 * Args:     cm        - the model, used only for it's alphabet and null model
 *           errbuf    - for reporting errors
 *           gamma     - the gamma data structure
 *           j         - offset j for gamma must be between 0 and gamma->L
 *           alpha_row - row of DP matrix to examine, we look at [dn..dx], NULL if we want to report
 *                       this j is IMPOSSIBLE end point of a hit (only possible if using_hmm_bands == TRUE)
 *           dn        - minimum d to look at 
 *           dx        - maximum d to look at
 *           using_hmm_bands - if TRUE, alpha_row is offset by dn, so we look at [0..dx-dn]
 *           bestr     - [dn..dx] root state (0 or local entry) corresponding to hit stored in alpha_row
 *           results   - results to add to, only used in this function if gamma->iamgreedy 
 *           W         - window size, max size of a hit, only used if we're doing a NULL3 correction (act != NULL)
 *           act       - [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1)
 *
 * Returns:  eslOK on succes; eslEMEM on memory allocation error;
 *
 */
int
UpdateGammaHitMxCM(CM_t *cm, char *errbuf, GammaHitMx_t *gamma, int j, float *alpha_row, int dn, int dx, int using_hmm_bands, 
		   int *bestr, search_results_t *results, int W, double **act)
{
  int status;
  int i, d;
  int bestd;
  int r;
  int dmin, dmax;
  int ip, jp;
  float *comp = NULL;    /* 0..a..cm->abc-K-1, the composition of residue a within the hit being reported */
  int a;
  float null3_correction;
  int do_report_hit;
  float hit_sc, cumulative_sc, bestd_sc;

  if(alpha_row == NULL && (!using_hmm_bands)) cm_Fail("UpdateGammaHitMxCM(), alpha_row is NULL, but using_hmm_bands is FALSE.\n");
  dmin = (using_hmm_bands) ? 0                 : dn; 
  dmax = (using_hmm_bands) ? ESL_MIN(dx-dn, j) : dx;
  if(act != NULL) ESL_ALLOC(comp, sizeof(float) * cm->abc->K);

  /* mode 1: non-greedy  */
  if((! gamma->iamgreedy) || alpha_row == NULL) { 
    gamma->mx[j]     = gamma->mx[j-1] + 0; 
    gamma->gback[j]  = -1;
    gamma->savesc[j] = IMPOSSIBLE;
    gamma->saver[j]  = -1;

    if(alpha_row != NULL) { 
      for (d = dmin; d <= dmax; d++) {
	i = using_hmm_bands ? j-d+1-dn  : j-d+1;
	hit_sc = alpha_row[d];
	cumulative_sc = gamma->mx[i-1] + hit_sc;
	/*printf("CAND hit %3d..%3d: %8.2f\n", i, j, hit_sc);*/
	if (cumulative_sc > gamma->mx[j]) {
	  do_report_hit = TRUE;
	  if(act != NULL && NOT_IMPOSSIBLE(hit_sc)) { /* do NULL3 score correction */
	    for(a = 0; a < cm->abc->K; a++) { 
	      comp[a] = act[j%(W+1)][a] - act[(i-1)%(W+1)][a]; 
	      /*printf("a: %5d j/W: %5d i-1/W: %5d j[a]: %.3f i-1[a]: %.3f c[a]: %.3f\n", a, j%(W+1), (i-1%W), act[(j%(W+1))][a], act[((i-1)%(W+1))][a], comp[a]);*/
	    }
	    esl_vec_FNorm(comp, cm->abc->K);
	    ScoreCorrectionNull3(cm->abc, cm->null, comp, j-i+1, cm->null3_omega, &null3_correction);
	    hit_sc -= null3_correction;
	    cumulative_sc -= null3_correction;
	    do_report_hit = (cumulative_sc > gamma->mx[j]) ? TRUE : FALSE;
	    /*printf("GOOD hit %3d..%3d: %8.2f  %10.6f  %8.2f\n", i, j, hit_sc+null3_correction, null3_correction, hit_sc);*/
	  }
	  if(do_report_hit) { 
	    gamma->mx[j]     = cumulative_sc;
	    gamma->gback[j]  = i + (gamma->i0-1);
	    gamma->savesc[j] = hit_sc;
	    gamma->saver[j]  = bestr[d]; 
	  }
	}
      }
    }
  }
  /* mode 2: greedy */
  else if(gamma->iamgreedy && dmin <= dmax) { /* if dmin >= dmax, no valid d for this j exists, don't report any hits */
    /* Resolving overlaps greedily (RSEARCH style),  
     * At least one hit is sent back for each j here.
     * However, some hits can already be removed for the greedy overlap
     * resolution algorithm.  Specifically, at the given j, any hit with a
     * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
    /* First, report hit with d of dmin (min valid d) if >= cutoff */
    hit_sc = alpha_row[dmin];
    if (hit_sc >= gamma->cutoff && NOT_IMPOSSIBLE(hit_sc)) {
      do_report_hit = TRUE;
      r = bestr[dmin]; 
      ip = using_hmm_bands ? j-(dmin+dn)+gamma->i0 : j-dmin+gamma->i0;
      jp = j-1+gamma->i0;
      assert(ip >= gamma->i0);
      assert(jp >= gamma->i0);
      if(act != NULL) { /* do NULL3 score correction */
	for(a = 0; a < cm->abc->K; a++) comp[a] = act[j%(W+1)][a] - act[(j-dmin+1-1)%(W+1)][a];
	esl_vec_FNorm(comp, cm->abc->K);
	ScoreCorrectionNull3(cm->abc, cm->null, comp, jp-ip+1, cm->null3_omega, &null3_correction);
	hit_sc -= null3_correction;
	do_report_hit = (hit_sc >= gamma->cutoff) ? TRUE : FALSE;
      }
      if(do_report_hit) { 
	/*printf("\t0 %.3f %.3f ip: %d jp: %d r: %d\n", hit_sc+null3_correction, hit_sc, ip, jp, r);*/
	ReportHit (ip, jp, r, hit_sc, results);
      }
    }
    bestd    = dmin;
    bestd_sc = hit_sc;
    /* Now, if current score is greater than maximum seen previous, report
     * it if >= cutoff and set new max */
    for (d = dmin+1; d <= dmax; d++) {
      hit_sc = alpha_row[d];
      if (hit_sc > bestd_sc) {
	if (hit_sc >= gamma->cutoff && NOT_IMPOSSIBLE(hit_sc)) { 
	  do_report_hit = TRUE;
	  r = bestr[d]; 
	  ip = using_hmm_bands ? j-(d+dn)+gamma->i0 : j-d+gamma->i0;
	  jp = j-1+gamma->i0;
	  assert(ip >= gamma->i0);
	  assert(jp >= gamma->i0);
	  if(act != NULL) { /* do NULL3 score correction */
	    for(a = 0; a < cm->abc->K; a++) comp[a] = act[j%(W+1)][a] - act[(j-d+1-1)%(W+1)][a];
	    esl_vec_FNorm(comp, cm->abc->K);
	    ScoreCorrectionNull3(cm->abc, cm->null, comp, jp-ip+1, cm->null3_omega, &null3_correction);
	    hit_sc -= null3_correction;
	    do_report_hit = ((hit_sc > bestd_sc) && (hit_sc >= gamma->cutoff)) ? TRUE : FALSE;
	  }
	  if(do_report_hit) { 
	    /*printf("\t1 %.3f %.3f ip: %d jp: %d r: %d\n", hit_sc+null3_correction, hit_sc, ip, jp, r);*/
	    ReportHit (ip, jp, r, hit_sc, results);
	  }
	}
	if(hit_sc > bestd_sc) { bestd = d; bestd_sc = hit_sc; } /* we need to check again b/c if null3, hit_sc -= null3_correction */
      }
    }
  }
  if(comp != NULL) free(comp);
  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return eslEMEM; /* NEVERREACHED */
}


/* Function: UpdateGammaHitMxCP9Forward()
 * Date:     EPN, Wed Nov 28 16:55:18 2007
 *
 * Purpose:  Update a gamma semi-HMM for forward-direction CP9 hits that span i..j.
 * 
 * Args:     cp9       - the model, used only for it's alphabet and null model
 *           errbuf    - for reporting errors
 *           gamma     - the gamma data structure
 *           i         - offset i for gamma must be between 0 and gamma->L, this is usually max(j-W+1, i0)
 *           j         - offset j for gamma must be between 0 and gamma->L
 *           sc        - score of best hit for i..j
 *           results   - results to add to, only used in this function if gamma->iamgreedy 
 *           W         - window size, max size of a hit, only used if we're doing a NULL3 correction (act != NULL)
 *           act       - [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1)
 * Returns:  eslOK on succes; eslEMEM on memory allocation error;
 */
int
UpdateGammaHitMxCP9Forward(CP9_t *cp9, char *errbuf, GammaHitMx_t *gamma, int i, int j, float hit_sc, search_results_t *results, int W, double **act)
{
  int status;
  float cumulative_sc;
  int a;
  int ip, jp;
  float *comp = NULL;    /* 0..a..cm->abc-K-1, the composition of residue a within the hit being reported */
  float null3_correction;

  /* update hit_sc if do_null3 == TRUE (do_null3 == TRUE iff act != NULL) */
  if(act != NULL) { /* do NULL3 score correction */
    /*printf("\ti: %4d i: %4d j: %4d hit_sc: %.2f ", i, i, j, hit_sc);*/
    ESL_ALLOC(comp, sizeof(float) * cp9->abc->K);
    for(a = 0; a < cp9->abc->K; a++) comp[a] = act[j%(W+1)][a] - act[(i-1)%(W+1)][a]; /* comp[a] is act[j][a] - act[i-1][a] */
    esl_vec_FNorm(comp, cp9->abc->K);
    /*esl_vec_FDump(stdout, comp, cp9->abc->K, NULL);*/
    /*printf("(A:%.2f C:%.2f G:%.2f U:%.2f)  ", comp[0], comp[1], comp[2], comp[3]);*/
    ScoreCorrectionNull3(cp9->abc, cp9->null, comp, j-i+1, cp9->null3_omega, &null3_correction);
    hit_sc -= null3_correction;
  }
  /*printf("%.2f (%.2f)\n", hit_sc, null3_correction);*/
  /* mode 1: non-greedy  */
  if(! gamma->iamgreedy) {
    gamma->mx[j]     = gamma->mx[j-1] + 0; 
    gamma->gback[j]  = -1;
    gamma->savesc[j] = IMPOSSIBLE;
    gamma->saver[j]  = -1;

    cumulative_sc = gamma->mx[i-1] + hit_sc;
    if (cumulative_sc > gamma->mx[j]) {
      gamma->mx[j]     = cumulative_sc;
      gamma->gback[j]  = i + (gamma->i0-1);
      gamma->savesc[j] = hit_sc;
      gamma->saver[j]  = 0; /* saver is invalid for CP9 hits */
    }
  }
  /* mode 2: greedy */
  if(gamma->iamgreedy) { 
    /* Return best hit for each j, IFF it's above threshold */
    if (hit_sc >= gamma->cutoff) { 
      ip = i-1 + gamma->i0;
      jp = j-1 + gamma->i0;
      ReportHit (ip, jp, 0, hit_sc, results); /* 0 is for saver, which is irrelevant for HMM hits */
    }
  }
  if(comp != NULL) free(comp);
  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return eslEMEM; /* NEVERREACHED */
}


/* Function: UpdateGammaHitMxCP9Backward()
 * Date:     EPN, Wed Nov 28 16:55:09 2007
 *
 * Purpose:  Update a gamma semi-HMM for backward-direction CP9 hits that span i..j.
 *
 *          Note: There's a 'backwards-specific' off-by-one here, that
 *          only occurs b/c we're going backwards, this is probably
 *          implementation specific (meaning getting rid of it is possible,
 *          but I haven't figured out how to), but we deal
 *          with it (albeit confusingly) as follows:
 * 
 *          '*off-by-one*' marked comments below refers to this issue:
 *          All Backward hits are rooted in M_O, the B (begin) state,
 *          which is * a non-emitter.  let i = ip+i0-1 => ip = i-i0+1;
 *          so sc[ip] = * backward->mmx[ip][0] = summed log prob of
 *          all parses that end at * j0, and start at position i+1 of
 *          the sequence (because i+1 is the * last residue whose
 *          emission has been accounted for). As a result, * gamma
 *          indexing is off-by-one with respect to sequence position,
 *          * hence the i+1 or i-1 in the following code blocks, each
 *          marked by * "*off-by-one*" comment below.  for example:
 *          let i0 = 2 gamma[ip=4], * normally this means ip=4
 *          corresponds to i=5 but due to this * off-by-one sc[ip=4]
 *          corresponds to hits that start at i=6
 *
 *          NOTE: (EPN, Fri Sep 25 10:23:01 2009) I introduced a 
 *          'bug' fix for UpdateGammaHitMxCP9Forward() (above) that is specific
 *          to NULL3 (see the comments in that function above). 
 *          The analogous bug fixing change was not made to *this* function b/c
 *          *in the current implementation* cp9_Backward (or cp9_ViterbiBackward)
 *          is only called after Forward has been run, so 'j' (the analog of i
 *          in Backward) is already known, and the problem of overestimating
 *          the hit length and inflating the null3 penalty goes away. 
 *          IF THIS IS CHANGED, if Backward were run as a scanner by itself, 
 *          without a j from Forward (though I can't imagine why you'd do this)
 *          one would want to think about adding the analogous bug fix in this 
 *          function. The reason NOT to do the bug fix here is that
 *          sometimes the hit length is longer than clen, in which
 *          case you want the null3 correction to be larger than it would be
 *          if the hit length were clen. See 
 *          /groups/eddy/home/nawrockie/notebook/9_0924_inf_rfam_question/00LOG
 *          for more details.
 *            
 * Args:     cp9       - the model, used only for it's alphabet and null model
 *           errbuf    - for reporting errors
 *           gamma     - the gamma data structure
 *           i         - offset i for gamma must be between 0 and gamma->L
 *           j         - offset j for gamma must be between 0 and gamma->L, this is usually i+W-1
 *           sc        - score of best hit for i..j
 *           results   - results to add to, only used in this function if gamma->iamgreedy 
 *           W         - window size, max size of a hit, only used if we're doing a NULL3 correction (act != NULL)
 *           act       - [0..j..W-1][0..a..abc->K-1], alphabet count, count of residue a in dsq from 1..jp where j = jp%(W+1)
 *
 * Returns:  eslOK on succes; eslEMEM on memory allocation error;
 */
int
UpdateGammaHitMxCP9Backward(CP9_t *cp9, char *errbuf, GammaHitMx_t *gamma, int i, int j, float hit_sc, search_results_t *results, int W, double **act)
{
  int status;
  float cumulative_sc;
  int a;
  int ip, jp;
  float *comp = NULL;    /* 0..a..abc-K-1, the composition of residue a within the hit being reported */
  float null3_correction;

  ESL_DASSERT1((i <= j));

  /* update hit_sc if do_null3 == TRUE (do_null3 == TRUE iff act != NULL) */
  if(act != NULL) { /* do NULL3 score correction */
    ESL_ALLOC(comp, sizeof(float) * cp9->abc->K);
    for(a = 0; a < cp9->abc->K; a++) comp[a] = act[(i+1)%(W+1)][a] - act[(j+1)%(W+1)][a]; /* *off-by-one*  we're going backwards: comp[a] is act[i][a] - act[j+1][a], but there's an off-by-one (see above) */
    esl_vec_FNorm(comp, cp9->abc->K);
    ScoreCorrectionNull3(cp9->abc, cp9->null, comp, j-(i+1)+1, cp9->null3_omega, &null3_correction);
    hit_sc -= null3_correction;
    /*esl_vec_FDump(stdout, comp, cp9->abc->K, NULL);*/
    /*printf("%.3f %.3f\n", hit_sc + null3_correction, hit_sc);*/
  }

  /* mode 1: non-greedy  */
  if(! gamma->iamgreedy) {
    gamma->mx[i]     = gamma->mx[i+1] + 0; 
    gamma->gback[i]  = -1;
    gamma->savesc[i] = IMPOSSIBLE;
    gamma->saver[i]  = -1;

    cumulative_sc = gamma->mx[j+1-1] + hit_sc; /* *off-by-one */
    if (cumulative_sc > gamma->mx[i]) {
      gamma->mx[i]     = cumulative_sc;
      gamma->gback[i]  = j + (gamma->i0-1);
      gamma->savesc[i] = hit_sc;
      gamma->saver[i]  = 0; /* saver is invalid for CP9 hits */
    }
  }

  /* mode 2: greedy */
  if(gamma->iamgreedy) { 
    /* Return best hit for each j, IFF it's above threshold */
    if (hit_sc >= gamma->cutoff) { 
      ip = i-1 + gamma->i0;
      jp = j-1 + gamma->i0;
      ReportHit (ip+1, jp, 0, hit_sc, results); /* 0 is for saver, which is irrelevant for HMM hits */
    }
  }
  if(comp != NULL) free(comp);
  return eslOK;

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Memory allocation error.\n");
  return eslEMEM; /* NEVERREACHED */
}


/* Function: TBackGammaHitMxForward()
 * Date:     EPN, Mon Nov  5 10:14:30 2007
 *
 * Purpose:  Traceback with a gamma semi-HMM for CM/CP9 hits in the forward
 *           direction. See TBackGammaHitMxBackward() for backward direction.
 *           gamma->iamgreedy should be FALSE.
 *            
 * Returns:  void; dies immediately upon an error.
 */
void
TBackGammaHitMxForward(GammaHitMx_t *gamma, search_results_t *results, int i0, int j0)
{
  int j, jp_g;

  if(gamma->iamgreedy) cm_Fail("cm_TBackGammaHitMx(), gamma->iamgreedy is TRUE.\n");   
  if(results == NULL)  cm_Fail("cm_TBackGammaHitMx(), results == NULL");
  /* Recover all hits: an (i,j,sc) triple for each one.
   */
  j = j0;
  while (j >= i0) {
    jp_g = j-i0+1;
    /*printf("TBACK j: %d sc: %.2f\n", j, gamma->savesc[jp_g]);*/
    if (gamma->gback[jp_g] == -1) j--; /* no hit */
    else {              /* a hit, a palpable hit */
      if(gamma->savesc[jp_g] >= gamma->cutoff) /* report the hit */
	ReportHit(gamma->gback[jp_g], j, gamma->saver[jp_g], gamma->savesc[jp_g], results);
      j = gamma->gback[jp_g]-1;
    }
  }
  return;
}


/* Function: TBackGammaHitMxBackward()
 * Date:     EPN, Wed Nov 28 16:53:48 2007
 *
 * Purpose:  Traceback with a gamma semi-HMM for CP9 hits in the backward
 *           direction. See TBackGammaHitMxForward() for forward direction.
 *           gamma->iamgreedy should be FALSE.
 *
 *           Note: Remember the 'backward-specific' off-by-one b/t seq 
 *           index and gamma (see *off-by-one* in 'Purpose' section 
 *           of UpdateGammaHitMxCP9Backward, above)
 *            
 * Returns:  void; dies immediately upon an error.
 */
void
TBackGammaHitMxBackward(GammaHitMx_t *gamma, search_results_t *results, int i0, int j0)
{
  int i, ip_g;

  if(gamma->iamgreedy) cm_Fail("cm_TBackGammaHitMx(), gamma->iamgreedy is TRUE.\n");   
  if(results == NULL)  cm_Fail("cm_TBackGammaHitMx(), results == NULL");
  /* Recover all hits: an (i,j,sc) triple for each one.
   */
  i = i0;
  while (i <= j0) {
    ip_g   = (i-1)-i0+1; /* *off-by-one*, i-1, ip corresponds to i+1 
			  * (yes: i *-* 1 =>   ip corresponds to i *+* 1) */
    if (gamma->gback[ip_g] == -1) i++; /* no hit */ 
    else {              /* a hit, a palpable hit */
      if(gamma->savesc[ip_g] >= gamma->cutoff) { /* report the hit */
	ESL_DASSERT1((i <= gamma->gback[ip_g]));
	ReportHit(i, gamma->gback[ip_g], gamma->saver[ip_g], gamma->savesc[ip_g], results);
      }
      i = gamma->gback[ip_g]+1;
    }
  }
  return;
}

#if 0
/* Function: cm_UpdateScanMatrix()
 * Date:     EPN, Wed Nov  7 12:49:36 2007
 *
 * Purpose:  Free, reallocate and recalculate ScanMatrix cm->smx>
 *           for CM <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_UpdateScanMatrixForCM(CM_t *cm)
{
  /* contract check */
  if(cm->flags & CMH_SCANMATRIX)    cm_Fail("cm_UpdateScanMatrix(), the CM flag for valid scan info is already up.");
  if(cm->smx->flags & cmSMX_HAS_FLOAT) {
    cm_FreeFloatsFromScanMatrix(cm, cm->smx);
    cm_FloatizeScanMatrix(cm, cm->smx);
  }
  if(cm->smx->flags & cmSMX_HAS_INT) {
    cm_FreeIntsFromScanMatrix(cm, cm->smx);
    cm_IntizeScanMatrix(cm, cm->smx);
  }
  cm->flags |= CMH_SCANMATRIX; /* ScanMatrix is valid now */
  return eslOK;
}
#endif


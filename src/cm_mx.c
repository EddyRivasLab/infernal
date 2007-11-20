/* CM_*MX* implementations: dynamic programming matrices for CMs
 * 
 * This is based heavily on HMMER 3's p7_gmx.c module.
 *
 * Table of contents:
 *   1. CM_HB_MX data structure functions,
 *      matrix of float scores for HMM banded CM alignment/search
 * EPN, Fri Oct 26 05:04:34 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"

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
   *          for deck M: only allocate the pointer, leave it pointing to NULL,
   *          deck M is special for 2 reasons: 
   *            1) only outside uses it, but it allocates it itself (cyk, inside don't use it)
   *            2) there are no bands on state M, the EL state
   */
  ESL_ALLOC(mx->dp,  sizeof(float **) * (M+1));
  mx->dp[M] = NULL;
 
  /* level 3: dp cell memory, when creating only allocate 1 cells per state, for j = 0, d = 0 */
  int allocL = 1;
  int allocW = 1;
  ESL_ALLOC(mx->dp_mem,  sizeof(float) * (M) * (allocL) * (allocW));
  ESL_ALLOC(mx->nrowsA, sizeof(int)      * (M));
  for (v = 0; v < M; v++) {
    ESL_ALLOC(mx->dp[v], sizeof(float *) * (allocL));
    mx->nrowsA[v] = allocL;
    mx->dp[v][0]  = mx->dp_mem + v * (allocL) * (allocW);
  }
  mx->M            = M;
  mx->ncells_alloc = M*(allocL)*(allocW);
  mx->ncells_valid = 0;
  mx->L            = allocL; /* allocL = 1 */

  mx->size_Mb = 
    (float) mx->M * (float) sizeof(int *) +    /* nrowsA ptrs */
    (float) mx->M * (float) sizeof(float **) + /* mx->dp[] ptrs */
    (float) mx->M * (float) sizeof(float *) +  /* mx->dp[v][] ptrs */
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
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Args:      mx     - the matrix to grow
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
cm_hb_mx_GrowTo(CM_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L)
{
  int     status;
  void   *p;
  int     v, jp, j;
  int     cur_size = 0;
  size_t  ncells;
  int     jbw;
  float   Mb_needed;

  /* contract check, number of states (M) is something we don't change
   * so check this matrix has same number of 1st dim state ptrs that
   * cp9b has */
  if(cp9b == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_mx_GrowTo() entered with cp9b == NULL.\n");
  if(cp9b->cm_M != mx->M) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_hb_mx_GrowTo() entered with mx->M: (%d) != cp9b->M (%d)\n", mx->M, cp9b->cm_M);
  
  ncells = 0;
  Mb_needed = ((float) (sizeof(int *)) * ((float) mx->M)) + /* nrowsA ptrs */
    (float) (sizeof(float **)) * (float) (mx->M); /* mx->dp[] ptrs */
  for(v = 0; v < mx->M; v++) { 
    jbw = (cp9b->jmax[v] - cp9b->jmin[v]); 
    Mb_needed += (float) (sizeof(float *) * jbw); /* mx->dp[v][] ptrs */
    for(jp = 0; jp <= jbw; jp++) 
      ncells += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
  }
  Mb_needed += ESL_MAX((float) (sizeof(float) * mx->ncells_alloc), (float) (sizeof(float) * ncells)); /* mx->dp_mem */
  Mb_needed *= 0.000001;  /* convert to Mb */
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
      cur_size    += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
    }
  }
  /* printf("ncells:   %d\n", ncells);
     printf("cur_size: %d\n", cur_size);*/

  mx->dp[mx->M] = NULL;
  mx->cp9b = cp9b; /* just a reference */
  
  /* finally, free the EL deck if it exists, this is the only reason
   * we keep track of mx->L, so we can free this if nec, 
   * the cm->M EL deck will be alloc'ed only if needed in FastOutsideAlignHB(): */
  if (mx->dp[mx->M] != NULL) { 
    for(j = 0; j <= mx->L; j++) { 
      if(mx->dp[mx->M][j] != NULL) free(mx->dp[mx->M][j]);
    }
    free(mx->dp[mx->M]);
    mx->dp[mx->M] = NULL;
  }	 
  /* now update L and size_Mb */
  mx->L       = L;    /* length of current seq we're valid for */
  mx->size_Mb = Mb_needed;
  
  ESL_DPRINTF1(("HMM banded mx->size_Mb: %.2f\n", mx->size_Mb));
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
  int v, j;

  if (mx->dp      != NULL) { 
    for (v = 0; v <= mx->M; v++) 
      if(mx->dp[v] != NULL) free(mx->dp[v]);  
  }
  /* free the mx->M deck if it exists */
  if (mx->dp[mx->M] != NULL) { 
    for(j = 0; j <= mx->L; j++) { 
      if(mx->dp[mx->M][j] != NULL) free(mx->dp[mx->M][j]);
    }
    free(mx->dp[mx->M]);
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
  return eslOK;
}

#if 0
/* Function:  cm_hb_mx_Size()
 * Incept:    EPN, Mon Nov 19 15:37:25 2007
 *
 * Purpose:   Return actual and needed size of a CM_HB_MX.
 *
 * Returns:   (void)
 *            Mb allocated in <ret_alloc>
 *            Mb currently needed (valid) in <ret_needed>
 *            Mb needed for an EL deck in <ret_el>
 */
void
cm_hb_mx_Size(CM_HB_MX *mx, float *ret_alloc, float *ret_needed, float *ret_el)
{
  float Mb_alloc;
  float Mb_needed;
  float Mb_el;
  int   v;
  int   jbw;

  Mb_alloc  = (float) (sizeof(float) * mx->ncells_alloc);
  Mb_needed = (float) (sizeof(float) * mx->ncells_valid);

  if(mx->ncells_valid != 0 && mx->cp9b != NULL) {  
    Mb_alloc  += (float) (sizeof(int *) * mx->M);   /* nrowsA ptrs */
    Mb_needed += (float) (sizeof(int *) * mx->M);   /* nrowsA ptrs */
    Mb_alloc  += (float) (sizeof(float **) * mx->M); /* dp 1st dim pts */
    Mb_needed += (float) (sizeof(float **) * mx->M); /* dp 1st dim pts */
    for(v = 0; v < mx->M; v++) {
      jbw = mx->cp9b->jmax[v] - mx->cp9b->jmin[v] + 1;
      Mb_alloc  += (float) (sizeof(float *) * jbw);  /* dp 2nd dim ptrs */
      Mb_needed += (float) (sizeof(float *) * jbw);  /* dp 2nd dim ptrs */
    }
  }

  Mb_el  = (float) (sizeof(float *) * (mx->L + 1));  /* EL deck pointers 0..j..L */
  Mb_el += (float) (sizeof(float) * (mx->L + 2) * (mx->L + 1) * 0.5);  /* EL deck cells */

  /* convert to Mb */
  Mb_alloc  *= 0.000001;
  Mb_needed *= 0.000001;
  Mb_el     *= 0.000001;

  if(ret_alloc  != NULL) *ret_alloc  = Mb_alloc;
  if(ret_needed != NULL) *ret_needed = Mb_needed;
  if(ret_el     != NULL) *ret_el     = Mb_el;

#if eslDEBUGLEVEL >= 1
  printf("\ncm_hb_mx_Size(): Mb_alloc:  %10.5f\n", Mb_alloc);
  printf("cm_hb_mx_Size(): Mb_needed: %10.5f\n", Mb_needed);
  printf("cm_hb_mx_Size(): Mb_el:     %10.5f\n", Mb_el);
#endif
  return;
}
#endif
/*****************************************************************
 * @LICENSE@
 *****************************************************************/


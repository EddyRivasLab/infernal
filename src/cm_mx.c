/* CM_*MX* implementations: dynamic programming matrices for CMs
 * 
 * This is based heavily on HMMER 3's p7_gmx.c module.
 *
 * EPN, Fri Oct 26 05:04:34 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"

#include "funcs.h"
#include "structs.h"

/* Function:  cm_fhb_mx_Create()
 * Incept:    EPN, Fri Oct 26 05:05:07 2007
 *
 * Purpose:   Allocate a reusable, resizeable <CM_FHBMX_VJD> for a CM
 *            given a CP9Bands_t object that defines the bands. 
 *            
 *            We've set this up so it should be easy to allocate
 *            aligned memory, though we're not doing this yet.
 *
 * Returns:   a pointer to the new <CM_FHB_MX>.
 *
 * Throws:    <NULL> on allocation error.
 */
CM_FHB_MX *
cm_fhb_mx_Create(int M)
{
  int     status;
  CM_FHB_MX *mx = NULL;
  int     v;

  /* level 1: the structure itself */
  ESL_ALLOC(mx, sizeof(CM_FHB_MX));
  mx->dp     = NULL;
  mx->dp_mem = NULL;
  mx->cp9b   = NULL;

  /* level 2: deck (state) pointers, 0.1..M, go all the way to M, cyk, inside don't use it, but outside does */
  ESL_ALLOC(mx->dp,  sizeof(float **) * (M+1));
 
  /* level 3: dp cell memory, when creating only allocate 1 cells per state, for j = 0, d = 0 */
  int allocL = 1;
  int allocW = 1;
  ESL_ALLOC(mx->dp_mem,  sizeof(float) * (M+1) * (allocL) * (allocW));
  
  ESL_ALLOC(mx->dp,     sizeof(float **) * (M+1));
  ESL_ALLOC(mx->nrowsA, sizeof(int)      * (M+1));
  for (v = 0; v <= M; v++) {
    ESL_ALLOC(mx->dp[v], sizeof(float *) * (allocL));
    mx->nrowsA[v] = allocL;
    mx->dp[v][0]  = mx->dp_mem + v * (allocL) * (allocW);
  }
  mx->M      = M;
  mx->ncells = (M+1)*(allocL)*(allocW);
  return mx;

 ERROR:
  if (mx != NULL) cm_fhb_mx_Destroy(mx);
  return NULL;
}

/* Function:  cm_fhb_mx_GrowTo()
 * Incept:    EPN, Fri Oct 26 05:19:49 2007
 *
 * Purpose:   Assures that a DP matrix <mx> is allocated
 *            for a model of exactly <mx->M> states and <mx->ncells>
 *            total cells. Determines new required size from 
 *            the CP9Bands_t object passed in, and reallocates if 
 *            necessary.
 *            
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Returns:   <eslOK> on success, and <mx> may be reallocated upon
 *            return; any data that may have been in <mx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <mx> must be assumed to be invalidated.
 */
int
cm_fhb_mx_GrowTo(CM_FHB_MX *mx, CP9Bands_t *cp9b)
{
  int     status;
  void   *p;
  int     i, v, jp;
  int     cur_size = 0;
  size_t  ncells;
  int     jbw;

  /* contract check, number of states (M) is something we don't change
   * so check this matrix has same number of 1st dim state ptrs that
   * cp9b has */
  if(cp9b == NULL)     cm_Fail("cm_fhb_mx_Growto() entered with cp9b == NULL.\n");
  if(cp9b->cm_M != mx->M) cm_Fail("cm_fhb_mx_Growto() entered with mx->M: (%d) != cp9b->M (%d)\n", mx->M, cp9b->cm_M);

  ncells = 0;
  for(v = 0; v < mx->M; v++)
    for(jp = 0; jp <= (cp9b->jmax[v] - cp9b->jmin[v]); jp++)
      ncells += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;

  /* must we realloc the full matrix? or can we get away with just
   * jiggering the pointers, if total required num cells is less
   * than or equal to what we already have alloc'ed?
   */
  if (ncells > mx->ncells) 
    {
      ESL_RALLOC(mx->dp_mem, p, sizeof(float) * ncells);
      mx->ncells = ncells;
    }

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

  /* TO DO: add mx->dp[M], the EL deck, somehow,
   * maybe do this separately, in a different data structure
   * b/c M doesn't have bands? or you could pass in a flag
   * to allocate this mem or not.
   */
  mx->dp[mx->M] = NULL;
  mx->cp9b = cp9b; /* just a reference */

  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_fhb_mx_Destroy()
 * Synopsis:  Frees a DP matrix.
 * Incept:    EPN, Fri Oct 26 09:04:04 2007
 *
 * Purpose:   Frees a <CM_FHB_MX>.
 *
 * Returns:   (void)
 */
void
cm_fhb_mx_Destroy(CM_FHB_MX *mx)
{
  if (mx == NULL) return;

  if (mx->dp      != NULL)  free(mx->dp);
  if (mx->dp_mem  != NULL)  free(mx->dp_mem);
  free(mx);
  return;
}

/* Function:  cm_fhb_mx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 * Incept:    EPN, Fri Oct 26 09:04:46 2007
 *
 * Purpose:   Dump matrix <mx> to stream <fp> for diagnostics.
 */
int
cm_fhb_mx_Dump(FILE *ofp, CM_FHB_MX *mx)
{
  int v, jp, j, dp, d;

  fprintf(ofp, "M: %d\nncells: %d\n", mx->M, mx->ncells);
  
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




/*****************************************************************
 * @LICENSE@
 *****************************************************************/

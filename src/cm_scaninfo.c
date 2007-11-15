/* ScanInfo_t implementations: information for scanning CYK/Inside DP functions.
 * 
 * EPN, Wed Nov  7 10:14:46 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"


/* Function: cm_CreateScanInfo()
 * Date:     EPN, Sun Nov  4 19:56:58 2007
 *
 * Purpose:  Given a CM, allocate and initialize ScanInfo_t object for that CM. 
 *            
 * Returns:  eslOK on success, dies immediately on some error
 */
int
cm_CreateScanInfo(CM_t *cm, int do_float, int do_int)
{
  int status;
  int j, v;
  int do_banded;

  if(! (cm->flags & CMH_BITS)) cm_Fail("cm_CreateScanInfo(), the CM flag for valid bit scores is down.");
  if(cm->flags & CMH_SCANINFO) cm_Fail("cm_CreateScanInfo(), the CM flag for valid scan info is already up.");
  if(cm->si != NULL) cm_Fail("cm_CreateScanInfo, the cm already points to a ScanInfo_t object.\n");
  if((! do_float) && (!do_int)) cm_Fail("cm_CreateScanInfo, do_float and do_int both FALSE.\n");
  if(cm->dmin == NULL && cm->dmax != NULL) cm_Fail("cm_CreateScanInfo(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmin != NULL && cm->dmax == NULL) cm_Fail("cm_CreateScanInfo(), cm->dmin == NULL, cm->dmax != NULL\n"); 
  if(cm->dmax != NULL && cm->W != cm->dmax[0]) cm_Fail("cm_CreateScanInfo(), cm->W: %d != cm->dmax[0]: %d\n", cm->W, cm->dmax[0]); 
  if((! cm->search_opts & CM_SEARCH_NOQDB) && (cm->dmin == NULL || cm->dmax == NULL))
    cm_Fail("cm_CreateScanInfo(), cm->dmin == NULL || cm->dmax == NULL, but !(cm->search_opts & CM_SEARCH_NOQDB)\n");

  ScanInfo_t *si;
  ESL_ALLOC(si, sizeof(ScanInfo_t));

  si->flags = 0;
  si->cm_M  = cm->M;
  si->W     = cm->W;
  si->dmin  = cm->dmin; /* could be NULL */
  si->dmax  = cm->dmax; /* could be NULL */
  do_banded = (cm->search_opts & CM_SEARCH_NOQDB) ? FALSE : TRUE;

  /* precalculate minimum and maximum d for each state and each sequence index (1..j..W). 
   * this is not always just dmin, dmax, (for ex. if j < W). */
  ESL_ALLOC(si->dnAA, sizeof(int *) * (si->W+1));
  ESL_ALLOC(si->dxAA, sizeof(int *) * (si->W+1));
  si->dnAA[0] = si->dxAA[0] = NULL; /* corresponds to j == 0, which is out of bounds */
  for(j = 1; j <= si->W; j++) {
    ESL_ALLOC(si->dnAA[j], sizeof(int) * cm->M);
    ESL_ALLOC(si->dxAA[j], sizeof(int) * cm->M);
    for(v = 0; v < cm->M; v++) {
      if(do_banded) { 
	si->dnAA[j][v] = (cm->sttype[v] == MP_st) ? ESL_MAX(si->dmin[v], 2) : ESL_MAX(si->dmin[v], 1); 
	si->dxAA[j][v] = ESL_MIN(j, si->dmax[v]); 
	si->dxAA[j][v] = ESL_MIN(si->dxAA[j][v], si->W);
      }
      else { 
	si->dnAA[j][v] = (cm->sttype[v] == MP_st) ? 2 : 1;
	si->dxAA[j][v] = ESL_MIN(j, si->W); 
      }
    }
  }
  /* allocate bestr, which holds best root state at alpha[0][cur][d] */
  ESL_ALLOC(si->bestr, (sizeof(int) * (si->W+1)));

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

  si->falpha       = NULL;
  si->falpha_begl  = NULL;

  si->ialpha      = NULL;
  si->ialpha_begl = NULL;

  cm->si = si;
  if(do_float) /* allocate float mx and scores */
    cm_FloatizeScanInfo(cm);
  if(do_int)   /* allocate int mx and scores */
    cm_IntizeScanInfo(cm);
  cm->flags |= CMH_SCANINFO; /* raise the flag for valid CMH_SCANINFO */
  return eslOK;

 ERROR:
  cm_Fail("memory allocation error in cm_CreateScanInfo().\n");
  return status; /* NEVERREACHED */
}


/* Function: cm_UpdateScanInfo()
 * Date:     EPN, Wed Nov  7 12:49:36 2007
 *
 * Purpose:  Free, reallocate and recalculate the ScanInfo for a CM.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_UpdateScanInfo(CM_t *cm)
{
  /* contract check */
  if(cm->flags & CMH_SCANINFO)    cm_Fail("cm_UpdateScanInfo(), the CM flag for valid scan info is already up.");
  if(cm->si->flags & cmSI_HAS_FLOAT) {
    cm_FreeFloatsFromScanInfo(cm);
    cm_FloatizeScanInfo(cm);
  }
  if(cm->si->flags & cmSI_HAS_INT) {
    cm_FreeIntsFromScanInfo(cm);
    cm_IntizeScanInfo(cm);
  }
  cm->flags |= CMH_SCANINFO; /* ScanInfo is valid now */
  return eslOK;
}

/* Function: cm_FloatizeScanInfo()
 * Date:     EPN, Wed Nov  7 10:05:55 2007
 *
 * Purpose:  Allocate and initialize float data structures in a ScanInfo_t object for <cm>.
 *           This initializes a scanning float DP matrix for CYK/Inside, for details on that
 *           matrix see the notes by the cm_FloatizeScanInfo() function call in 
 *           cm_CreateScanInfo().
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FloatizeScanInfo(CM_t *cm)
{
  int status;
  int j, v;
  int d, y, yoffset, w;
  int do_banded = (cm->search_opts & CM_SEARCH_NOQDB) ? FALSE : TRUE;

  /* contract check */
  if(cm->si == NULL) cm_Fail("cm_FloatizeScanInfo(), cm->si is NULL.\n");
  ScanInfo_t *si = cm->si;
  if(si->flags & cmSI_HAS_FLOAT) cm_Fail("cm_FloatizeScanInfo(), si's cmSI_HAS_FLOAT flag is already up.");
  if(si->falpha != NULL)       cm_Fail("cm_FloatizeScanInfo(), si->falpha is not NULL.");
  if(si->falpha_begl != NULL)  cm_Fail("cm_FloatizeScanInfo(), si->falpha_begl is not NULL.");
  
  /* allocate alpha */
  ESL_ALLOC(si->falpha,        sizeof(float **) * 2);
  ESL_ALLOC(si->falpha[0],     sizeof(float *) * cm->M);
  ESL_ALLOC(si->falpha[1],     sizeof(float *) * cm->M);
  ESL_ALLOC(si->falpha[0][0],  sizeof(float) * 2 * (cm->M) * (si->W+1));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      si->falpha[0][v] = si->falpha[0][0] + (v           * (si->W+1));
      si->falpha[1][v] = si->falpha[0][0] + ((v + cm->M) * (si->W+1));
    }
    else si->falpha[0][v] = si->falpha[1][v] = NULL; /* BEGL_S */
  }
  /* allocate falpha_begl */
  ESL_ALLOC(si->falpha_begl, (sizeof(float **) * (si->W+1)));
  for (j = 0; j <= si->W; j++) {
    ESL_ALLOC(si->falpha_begl[j], (sizeof(float *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) ESL_ALLOC(si->falpha_begl[j][v], (sizeof(float) * (si->W+1)));
      else si->falpha_begl[j][v] = NULL; /* non-BEGL */
    }
  }
  /* initialize falpha and falpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      si->falpha[0][v][0] = IMPOSSIBLE;
      if (cm->sttype[v] == E_st) { 
	si->falpha[0][v][0] = si->falpha[1][v][0] = 0.;
	/* rest of E deck is IMPOSSIBLE, this is rewritten if QDB is on, (slightly wasteful). */
	for (d = 1; d <= si->W; d++) si->falpha[0][v][d] = si->falpha[1][v][d] = IMPOSSIBLE;
      }
      else if (cm->sttype[v] == MP_st) si->falpha[0][v][1] = si->falpha[1][v][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	si->falpha[0][v][0] = cm->endsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  si->falpha[0][v][0] = ESL_MAX(si->falpha[0][v][0], (si->falpha[0][y+yoffset][0] + cm->tsc[v][yoffset]));
	si->falpha[0][v][0] = ESL_MAX(si->falpha[0][v][0], IMPOSSIBLE);
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	si->falpha[0][v][0] = si->falpha_begl[0][w][0] + si->falpha[0][y][0]; 
      }
      si->falpha[1][v][0] = si->falpha[0][v][0];
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      si->falpha_begl[0][v][0] = cm->endsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	si->falpha_begl[0][v][0] = ESL_MAX(si->falpha_begl[0][v][0], (si->falpha[0][y+yoffset][0] + cm->tsc[v][yoffset])); /* careful: y is in si->falpha */
      si->falpha_begl[0][v][0] = ESL_MAX(si->falpha_begl[0][v][0], IMPOSSIBLE);
      for (j = 1; j <= si->W; j++) 
	si->falpha_begl[j][v][0] = si->falpha_begl[0][v][0];
    }
  }
  /* query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) { 
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < cm->dmin[v] && d <= si->W; d++)
	  for(j = 0; j < 2; j++) 
	    si->falpha[j][v][d]  = IMPOSSIBLE;
	for (d = cm->dmax[v]+1; d <= si->W; d++)
	  for(j = 0; j < 2; j++)
	    si->falpha[j][v][d] = IMPOSSIBLE;
      }
      else { /* not BEGL_S state */
	for (d = 0; d < cm->dmin[v] && d <= si->W; d++)
	  for(j = 0; j <= si->W; j++)
	    si->falpha_begl[j][v][d] = IMPOSSIBLE;
	for (d = cm->dmax[v]+1; d <= si->W; d++)
	  for(j = 0; j <= si->W; j++)
	    si->falpha_begl[j][v][d] = IMPOSSIBLE;
      }
    }
  }
  /* set the flag that tells us we've got valid floats */
  si->flags |= cmSI_HAS_FLOAT;
  return eslOK;

 ERROR: 
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: cm_IntizeScanInfo()
 * Date:     EPN, Wed Nov  7 10:10:39 2007
 *
 * Purpose:  Allocate and initialize int data structures in a ScanInfo_t object for <cm>.
 *           This initializes a scanning float DP matrix for CYK/Inside, for details on that
 *           matrix see the notes by the cm_FloatizeScanInfo() function call in 
 *           cm_CreateScanInfo().
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_IntizeScanInfo(CM_t *cm)
{
  int status;
  int v, j, d, y, yoffset, w;
  int do_banded = (cm->search_opts & CM_SEARCH_NOQDB) ? FALSE : TRUE;

  /* contract check */
  if(cm->si == NULL) cm_Fail("cm_IntizeScanInfo(), cm->si is NULL.\n");
  ScanInfo_t *si = cm->si;
  if(si->flags & cmSI_HAS_INT) cm_Fail("cm_IntizeScanInfo(), si's cmSI_HAS_INT flag is already up.");
  if(si->ialpha != NULL)       cm_Fail("cm_IntizeScanInfo(), si->ialpha is not NULL.");
  if(si->ialpha_begl != NULL)  cm_Fail("cm_IntizeScanInfo(), si->ialpha_begl is not NULL.");

  /* allocate ialpha */
  ESL_ALLOC(si->ialpha,        sizeof(int **) * 2);
  ESL_ALLOC(si->ialpha[0],     sizeof(int *) * cm->M);
  ESL_ALLOC(si->ialpha[1],     sizeof(int *) * cm->M);
  ESL_ALLOC(si->ialpha[0][0],  sizeof(int) * 2 * (cm->M) * (si->W+1));
  for (v = cm->M-1; v >= 0; v--) {	
    if (cm->stid[v] != BEGL_S) {
      si->ialpha[0][v] = si->ialpha[0][0] + (v           * (si->W+1));
      si->ialpha[1][v] = si->ialpha[0][0] + ((v + cm->M) * (si->W+1));
    }
    else si->ialpha[0][v] = si->ialpha[1][v] = NULL; /* BEGL_S */
  }
  /* allocate ialpha_begl */
  ESL_ALLOC(si->ialpha_begl, (sizeof(int **) * (si->W+1)));
  for (j = 0; j <= cm->W; j++) {
    ESL_ALLOC(si->ialpha_begl[j], (sizeof(int *) * (cm->M)));
    for (v = cm->M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) ESL_ALLOC(si->ialpha_begl[j][v], (sizeof(int) * (si->W+1)));
      else si->ialpha_begl[j][v] = NULL; /* non-BEGL */
    }
  }
  /* initialize ialpha and ialpha_begl */
  for(v = cm->M-1; v >= 0; v--) {
    if(cm->stid[v] != BEGL_S) {
      si->ialpha[0][v][0] = -INFTY;
      if (cm->sttype[v] == E_st) { 
	si->ialpha[0][v][0] = si->ialpha[1][v][0] = 0.;
	/* rest of E deck is -INFTY, this is rewritten if QDB is on, (slightly wasteful). */
	for (d = 1; d <= si->W; d++) si->ialpha[0][v][d] = si->ialpha[1][v][d] = -INFTY;
      }
      else if (cm->sttype[v] == MP_st) si->ialpha[0][v][1] = si->ialpha[1][v][1] = -INFTY;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) {
	y = cm->cfirst[v];
	si->ialpha[0][v][0] = cm->iendsc[v];
	for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	  si->ialpha[0][v][0] = ESL_MAX(si->ialpha[0][v][0], (si->ialpha[0][y+yoffset][0] + cm->itsc[v][yoffset]));
	si->ialpha[0][v][0] = ESL_MAX(si->ialpha[0][v][0], -INFTY);
      }
      else if (cm->sttype[v] == B_st) {
	w = cm->cfirst[v]; /* BEGL_S, left child state */
	y = cm->cnum[v];
	si->ialpha[0][v][0] = si->ialpha_begl[0][w][0] + si->ialpha[0][y][0]; 
      }
      si->ialpha[1][v][0] = si->ialpha[0][v][0];
    }
    else { /* v == BEGL_S */
      y = cm->cfirst[v];
      si->ialpha_begl[0][v][0] = cm->iendsc[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	si->ialpha_begl[0][v][0] = ESL_MAX(si->ialpha_begl[0][v][0], (si->ialpha[0][y+yoffset][0] + cm->itsc[v][yoffset])); /* careful: y is in alpha */
      si->ialpha_begl[0][v][0] = ESL_MAX(si->ialpha_begl[0][v][0], -INFTY);
      for (j = 1; j <= si->W; j++) 
	si->ialpha_begl[j][v][0] = si->ialpha_begl[0][v][0];
    }
  }
  /* query-dependent band imposition.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  if(do_banded) { 
    for (v = 0; v < cm->M; v++) { 
      if(cm->stid[v] != BEGL_S) {
	for (d = 0; d < cm->dmin[v] && d <= si->W; d++)
	  for(j = 0; j < 2; j++) 
	    si->ialpha[j][v][d] = -INFTY;
	for (d = cm->dmax[v]+1; d <= si->W; d++)
	  for(j = 0; j < 2; j++)
	    si->ialpha[j][v][d] = -INFTY;
      }
      else { /* not BEGL_S state */
	for (d = 0; d < cm->dmin[v] && d <= si->W; d++)
	  for(j = 0; j <= si->W; j++)
	    si->ialpha_begl[j][v][d] = -INFTY;
	for (d = cm->dmax[v]+1; d <= si->W; d++)
	  for(j = 0; j <= si->W; j++)
	    si->ialpha_begl[j][v][d] = -INFTY;
      }
    }
  }
  /* set the flag that tells us we've got valid ints */
  si->flags |= cmSI_HAS_INT;
  return eslOK;

 ERROR: 
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}  


/* Function: cm_FreeFloatsFromScanInfo()
 * Date:     EPN, Wed Nov  7 10:03:55 2007
 *
 * Purpose:  Free float data structures in a ScanInfo_t object for <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeFloatsFromScanInfo(CM_t *cm)
{
  int j, v;

  /* contract check */
  if(cm->si == NULL) cm_Fail("cm_FreeFloatsFromScanInfo(), cm->si is NULL.\n");
  ScanInfo_t *si = cm->si;
  if(! si->flags & cmSI_HAS_FLOAT)    cm_Fail("cm_FreeFloatsFromScanInfo(), si's cmSI_HAS_FLOAT flag is down.");
  if(si->falpha == NULL)       cm_Fail("cm_FreeFloatsFromScanInfo(), si->falpha is already NULL.");
  if(si->falpha_begl == NULL)  cm_Fail("cm_FreeFloatsFromScanInfo(), si->falpha_begl is already NULL.");

  free(si->falpha[0][0]);
  free(si->falpha[1]);
  free(si->falpha[0]);
  free(si->falpha);
  si->falpha = NULL;
  for (j = 0; j <= si->W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S) {
	free(si->falpha_begl[j][v]);
      }
    free(si->falpha_begl[j]);
  }
  free(si->falpha_begl);
  si->falpha_begl = NULL;
  si->flags &= ~cmSI_HAS_FLOAT;
  return eslOK;
}

/* Function: cm_FreeIntsFromScanInfo()
 * Date:     EPN, Wed Nov  7 09:56:01 2007
 *
 * Purpose:  Free int data structures in a ScanInfo_t object for <cm>.
 *            
 * Returns:  eslOK on success, dies immediately on an error.
 */
int
cm_FreeIntsFromScanInfo(CM_t *cm)
{
  int j, v;

  /* contract check */
  if(cm->si == NULL) cm_Fail("cm_FreeIntsFromScanInfo(), cm->si is NULL.\n");
  ScanInfo_t *si = cm->si;
  if(! si->flags & cmSI_HAS_INT)    cm_Fail("cm_FreeIntsFromScanInfo(), si's cmSI_HAS_INT flag is down.");
  if(si->ialpha == NULL)       cm_Fail("cm_FreeIntsFromScanInfo(), si->ialpha is already NULL.");
  if(si->ialpha_begl == NULL)  cm_Fail("cm_FreeIntsFromScanInfo(), si->ialpha_begl is already NULL.");

  free(si->ialpha[0][0]);
  free(si->ialpha[1]);
  free(si->ialpha[0]);
  free(si->ialpha);
  si->ialpha = NULL;
  for (j = 0; j <= si->W; j++) {
    for (v = 0; v < cm->M; v++) 
      if (cm->stid[v] == BEGL_S) {
	free(si->ialpha_begl[j][v]);
      }
    free(si->ialpha_begl[j]);
  }
  free(si->ialpha_begl);
  si->ialpha_begl = NULL;
  si->flags &= ~cmSI_HAS_INT;
  return eslOK;
}

/* Function: cm_FreeScanInfo()
 * Date:     EPN, Sun Nov  4 20:57:32 2007
 *
 * Purpose:  Free a ScanInfo_t object for <cm>.
 *            
 * Returns:  void
 */
void
cm_FreeScanInfo(CM_t *cm)
{
  int j;
  /* contract check */
  if(cm->si == NULL) cm_Fail("cm_FreeScanInfo(), cm->si is NULL.\n");
  ScanInfo_t *si = cm->si;
  for(j = 1; j <= si->W; j++) {
    free(si->dnAA[j]);
    free(si->dxAA[j]);
  }
  free(si->dnAA);
  free(si->dxAA);
  free(si->bestr);
  
  if(si->flags & cmSI_HAS_FLOAT) cm_FreeFloatsFromScanInfo(cm);
  if(si->flags & cmSI_HAS_INT)   cm_FreeIntsFromScanInfo(cm);
  free(si);
  cm->si = NULL;
  cm->flags &= ~CMH_SCANINFO; /* drop the 'cm has valid scaninfo' flag */
  return;
}

/* Function: cm_DumpScanInfoAlpha()
 * Date:     EPN, Tue Nov  6 05:11:26 2007
 *
 * Purpose:  Dump current alpha matrix (either float or int).
 *            
 * Returns:  void.
 */
void
cm_DumpScanInfoAlpha(CM_t *cm, int j, int i0, int doing_float)
{
  int d, v;
  int jp_g = j-i0+1; /* j is actual index in j, jp_g is offset j relative to start i0 (index in gamma* data structures) */
  int cur = j%2;
  int prv = (j-1)%2;
  int *dnA, *dxA;

  if(cm->si == NULL) cm_Fail("cm_DumpScanInfoAlpha(), cm->si is NULL.\n");
  ScanInfo_t *si = cm->si;
  if(doing_float && (! si->flags & cmSI_HAS_FLOAT)) cm_Fail("cm_DumpScanInfoAlpha(), trying to print float alpha, but cmSI_HAS_FLOAT flag is down.\n");
  if((! doing_float) && (! si->flags & cmSI_HAS_INT)) cm_Fail("cm_DumpScanInfoAlpha(), trying to print int alpha, but cmSI_HAS_INT flag is down.\n");

  int begl_prv = j-1 % (si->W+1);
  int begl_cur = j   % (si->W+1);

  printf("Dumping Alpha: j: %d\n", j);
  if(jp_g >= si->W) { dnA = si->dnAA[si->W]; dxA = si->dxAA[si->W]; }
  else              { dnA = si->dnAA[jp_g];  dxA = si->dxAA[jp_g]; }
  if(doing_float) {
    for (v = si->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) { 
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10.4f\n", (j-1), v, d, si->falpha_begl[begl_prv][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j-1:%4d][%4d][%4d]: %10.4f\n", (j-1), v, d, si->falpha[prv][v][d]); 
      }
      printf("\n");
    }
    for (v = si->cm_M-1; v >= 0; v--) {	
      if (cm->stid[v] == BEGL_S) {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10.4f\n", j,     v, d, si->falpha_begl[begl_cur][v][d]); 
      }
      else {
	for(d = dnA[v]; d <= dxA[v]; d++) printf("A[j  :%4d][%4d][%4d]: %10.4f\n", j,     v, d, si->falpha[cur][v][d]); 
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

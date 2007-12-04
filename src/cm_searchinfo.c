/* SearchInfo_t implementations: information for CM/CP9 
 * filters and scans.
 *
 * EPN, Tue Nov 27 08:42:08 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"


/* Function: cm_CreateSearchInfo()
 * Date:     EPN, Tue Nov 27 12:57:24 2007
 *
 * Purpose:  Allocate and initialize a search info object that
 *           specifies that no filtering should be done. 
 *            
 * Returns:  cm->si points to a new SearchInfo_t object
 *           eslOK on success, dies immediately on some error
 */
int
cm_CreateSearchInfo(CM_t *cm, int cutoff_type, float cutoff)
{
  int status;

  if(cm->si != NULL)  cm_Fail("cm_CreateSearchInfo(), the cm already points to a SearchInfo_t object.\n");
  
  SearchInfo_t *si;
  ESL_ALLOC(si, sizeof(SearchInfo_t));
  
  si->nrounds = 0;
  ESL_ALLOC(si->search_opts, sizeof(int)                * (si->nrounds+1));
  ESL_ALLOC(si->cutoff_type, sizeof(int)                * (si->nrounds+1));
  ESL_ALLOC(si->cutoff,      sizeof(float)              * (si->nrounds+1));
  ESL_ALLOC(si->stype,       sizeof(int)                * (si->nrounds+1));
  ESL_ALLOC(si->smx,         sizeof(ScanMatrix_t *)     * (si->nrounds+1));
  ESL_ALLOC(si->hsi,         sizeof(HybridScanInfo_t *) * (si->nrounds+1));

  si->search_opts[0] = cm->search_opts;
  si->cutoff_type[0] = cutoff_type;
  si->cutoff[0]      = cutoff;
  si->stype[0]       = SEARCH_WITH_CM;
  si->smx[0]         = cm->smx;  /* could be NULL */
  si->hsi[0]         = NULL;

  cm->si = si;
  return eslOK;

 ERROR:
  cm_Fail("cm_CreateSearchInfo(), memory allocate error.");
  return status; /* NEVERREACHED */
}  


/* Function: cm_AddFilterToSearchInfo()
 * Date:     EPN, Tue Nov 27 13:00:23 2007
 *
 * Purpose:  Add an HMM filter as the 1st round filter for CM <cm>.
 *           A new SearchInfo_t object <fi> is created, and the existing
 *           information from cm->si is copied into it. cm->si is then
 *           freed and cm->si is set to point at fi.            
 * 
 * Returns:  cm->si points to a new SearchInfo_t object
 *           eslOK on success, dies immediately on some error
 */
int
cm_AddFilterToSearchInfo(CM_t *cm, int cyk_filter, int inside_filter, int viterbi_filter, int forward_filter, 
			 int hybrid_filter, ScanMatrix_t *smx, HybridScanInfo_t *hsi, int cutoff_type, float cutoff)
{
  int status;
  int n;
  int orig_nrounds;

  if(cm->si == NULL)                 cm_Fail("cm_AddFilterToSearchInfo(), the cm does not point to a SearchInfo_t object.\n");
  if((cyk_filter + inside_filter + viterbi_filter + forward_filter + hybrid_filter) != 1)
    cm_Fail("cm_AddFilterToSearchInfo(), cyk_filter: %d\ninside_filter: %d\nviterbi_filter: %d\nforward_filter: %d\nhybrid_filter: %d. Exactly 1 of these must be 1, the rest 0s.\n", cyk_filter, inside_filter, viterbi_filter, forward_filter, hybrid_filter);
  if(cyk_filter    && smx == NULL) cm_Fail("cm_AddFilterToSearchInfo(), cyk_filter: %d but smx == NULL\n", cyk_filter);
  if(inside_filter && smx == NULL) cm_Fail("cm_AddFilterToSearchInfo(), inside_filter: %d but smx == NULL\n", inside_filter);
  if(hybrid_filter && hsi == NULL) cm_Fail("cm_AddFilterToSearchInfo(), hybrid_filter: %d but hsi == NULL\n", hybrid_filter);

  orig_nrounds = cm->si->nrounds;

  /* allocate new si object, with 1 more round than cm->si, and set round 0 */
  SearchInfo_t *si;
  ESL_ALLOC(si, sizeof(SearchInfo_t));
  si->nrounds = orig_nrounds+1;
  ESL_ALLOC(si->search_opts, sizeof(int)   * (si->nrounds+1));
  ESL_ALLOC(si->cutoff_type, sizeof(int)   * (si->nrounds+1));
  ESL_ALLOC(si->cutoff,      sizeof(float) * (si->nrounds+1));
  ESL_ALLOC(si->stype,       sizeof(int)   * (si->nrounds+1));
  ESL_ALLOC(si->smx,         sizeof(ScanMatrix_t     *) * (si->nrounds+1));
  ESL_ALLOC(si->hsi,         sizeof(HybridScanInfo_t *) * (si->nrounds+1));

  si->search_opts[0] = 0;
  si->search_opts[0] |= CM_SEARCH_NOALIGN;
  if(cyk_filter)     ;/* do nothing, CYK is default */
  if(inside_filter)  si->search_opts[0] |= CM_SEARCH_INSIDE;
  if(viterbi_filter) si->search_opts[0] |= CM_SEARCH_HMMVITERBI;
  if(forward_filter) si->search_opts[0] |= CM_SEARCH_HMMFORWARD;
  
  si->cutoff_type[0] = cutoff_type;
  si->cutoff[0]      = cutoff;
  if(viterbi_filter || forward_filter) { 
    si->stype[0] = SEARCH_WITH_HMM;
    si->smx[0]   = NULL;
    si->hsi[0]   = NULL;
  }
  if(cyk_filter || inside_filter)  { 
    si->stype[0] = SEARCH_WITH_CM;
    si->smx[0]   = smx;
    si->hsi[0]   = NULL;
  }
  si->hsi[0] = NULL;
  if(hybrid_filter) { 
    si->stype[0] = SEARCH_WITH_HYBRID;
    si->smx[0]   = NULL;
    si->hsi[0]   = hsi;
  }
       
  /* copy existing information for other rounds from old cm->si */
  for(n = 0; n <= orig_nrounds; n++) { 
    si->search_opts[(n+1)] = cm->si->search_opts[n];
    si->cutoff_type[(n+1)] = cm->si->cutoff_type[n];
    si->cutoff[(n+1)]      = cm->si->cutoff[n];
    si->stype[(n+1)]       = cm->si->stype[n];
    /* and copy the ptr to smx and hsi */
    si->smx[(n+1)] = cm->si->smx[n];
    si->hsi[(n+1)] = cm->si->hsi[n];
  }    

  /* free old cm->si, but be careful not to free smx[] and hsi[], we're still pointing to those */
  free(cm->si->search_opts);
  free(cm->si->cutoff_type);
  free(cm->si->cutoff);
  free(cm->si->stype);
  free(cm->si->smx);
  free(cm->si->hsi);
  free(cm->si);

  cm->si = si;
  return eslOK;

 ERROR:
  cm_Fail("cm_AddFilterToSearchInfo(), memory allocate error.");
  return status; /* NEVERREACHED */
}  

/* Function: cm_FreeSearchInfo()
 * Date:     EPN, Tue Nov 27 08:48:49 2007
 *
 * Purpose:  Free a SearchInfo_t object corresponding to 
 *           CM <cm>.
 *            
 * Returns:  void
 */
void
cm_FreeSearchInfo(SearchInfo_t *si, CM_t *cm)
{
  int n;

  for(n = 0; n <  si->nrounds; n++) if(si->smx[n] != NULL) cm_FreeScanMatrix(cm, si->smx[n]);
  /* we don't free si->smx[nrounds] b/c it == cm->smx */
  for(n = 0; n <= si->nrounds; n++) if(si->hsi[n] != NULL) cm_FreeHybridScanInfo(si->hsi[n], cm); 
  free(si->search_opts);
  free(si->cutoff_type);
  free(si->cutoff);
  free(si->stype);
  free(si->smx);
  free(si->hsi);

  free(si);
  return;
}

/* Function: cm_DumpSearchInfo()
 * Date:     EPN, Tue Nov 27 08:50:54 2007
 *
 * Purpose:  Dump a CM's search info (except scan matrix and hybrid scan info) to stdout. 
 *            
 * Returns:  void.
 */
void
cm_DumpSearchInfo(SearchInfo_t *si)
{
  int n, v;
  printf("\nSearchInfo summary:\n");
  printf("nrounds: %d\n", si->nrounds);
  for(n = 0; n <= si->nrounds; n++) { 
    printf("\nround: %d\n", n);
    if(si->stype[n] == SEARCH_WITH_HMM)    printf("\ttype: HMM\n"); 
    if(si->stype[n] == SEARCH_WITH_HYBRID) printf("\ttype: Hybrid\n"); 
    if(si->stype[n] == SEARCH_WITH_CM)     printf("\ttype: CM\n"); 
    DumpSearchOpts(si->search_opts[n]);
    if(si->cutoff_type[n] == SCORE_CUTOFF) printf("\tcutoff     : %10.4f bits\n", si->cutoff[n]);
    else                                   printf("\tcutoff     : %10.4f E-value\n", si->cutoff[n]);
    if(si->hsi[n] != NULL) { 
      printf("\tHybrid info:\n");
      printf("\t\tNumber of sub CM roots: %d\n", si->hsi[n]->n_v_roots);
      for(v = 0; v < si->hsi[n]->cm_M; v++) 
	if(si->hsi[n]->v_isroot[v]) printf("\t\tstate %d is a root\n", v);
    }
  }
  return;
}


/* Function: cm_ValidateSearchInfo()
 * Date:     EPN, Tue Nov 27 09:25:10 2007
 *
 * Purpose:  Validate a Search Info <si> object for CM <cm>.
 *            
 * Returns:  void.
 */
void
cm_ValidateSearchInfo(CM_t *cm, SearchInfo_t *si)
{
  int n, sum;
  int do_noqdb;
  int do_hbanded;
  int do_hmmscanbands;
  int do_sums;
  int do_inside;
  int do_toponly;
  int do_noalign;
  int do_null2;
  int do_rsearch;
  int do_cmgreedy;
  int do_hmmgreedy;
  int do_hmmviterbi;
  int do_hmmforward;

  for(n = 0; n <= si->nrounds; n++) { 
    do_noqdb       = (si->search_opts[n] & CM_SEARCH_NOQDB)        ? TRUE : FALSE;
    do_hbanded     = (si->search_opts[n] & CM_SEARCH_HBANDED)      ? TRUE : FALSE;
    do_hmmscanbands= (si->search_opts[n] & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;
    do_sums        = (si->search_opts[n] & CM_SEARCH_SUMS)         ? TRUE : FALSE;
    do_inside      = (si->search_opts[n] & CM_SEARCH_INSIDE)       ? TRUE : FALSE;
    do_toponly     = (si->search_opts[n] & CM_SEARCH_TOPONLY)      ? TRUE : FALSE;
    do_noalign     = (si->search_opts[n] & CM_SEARCH_NOALIGN)      ? TRUE : FALSE;
    do_null2       = (si->search_opts[n] & CM_SEARCH_NULL2)        ? TRUE : FALSE;
    do_rsearch     = (si->search_opts[n] & CM_SEARCH_RSEARCH)      ? TRUE : FALSE;
    do_cmgreedy    = (si->search_opts[n] & CM_SEARCH_CMGREEDY)     ? TRUE : FALSE;
    do_hmmgreedy   = (si->search_opts[n] & CM_SEARCH_HMMGREEDY)    ? TRUE : FALSE;
    do_hmmviterbi  = (si->search_opts[n] & CM_SEARCH_HMMVITERBI)   ? TRUE : FALSE;
    do_hmmforward  = (si->search_opts[n] & CM_SEARCH_HMMFORWARD)   ? TRUE : FALSE;

    if(n < si->nrounds) { 
      if(!do_noalign) cm_Fail("cm_ValidateSearchInfo(), round %d has CM_SEARCH_NOALIGN flag down.\n", n);
      if(si->stype[n] == SEARCH_WITH_HMM) {
	sum = do_noqdb + do_hbanded + do_hmmscanbands + do_sums + do_inside + do_toponly + do_null2 + do_rsearch + do_cmgreedy;
	if(sum != 0 || (do_hmmviterbi + do_hmmforward != 1)) { 
	  printf("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HMM but search opts are invalid\n", n);
	  DumpSearchOpts(si->search_opts[n]);
	  cm_Fail("This is fatal.");
	}
	if(si->smx[n] != NULL) cm_Fail("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HMM but smx[%d] is non-NULL\n", n, n);
	if(si->hsi[n] != NULL) cm_Fail("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HMM but hsi[%d] is non-NULL\n", n, n);
      }
      else if (si->stype[n] == SEARCH_WITH_HYBRID) {
	if(si->hsi[n] == NULL)      cm_Fail("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but hsi[%d] is NULL\n", n, n);
	if(si->hsi[n]->smx == NULL) cm_Fail("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but hsi[%d]->smx is NULL\n", n, n);
	if(si->smx[n] != NULL)      cm_Fail("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but smx[%d] is not NULL\n", n, n);
	if(si->hsi[n]->v_isroot[0]) cm_Fail("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID and hsi->vi_isroot[0] is TRUE, this shouldn't happen, we might as well filter with a SEARCH_WITH_CM filter.");
	sum = do_hbanded + do_hmmscanbands + do_sums + do_toponly + do_null2 + do_rsearch + do_cmgreedy;	
	if(sum != 0 || (do_hmmviterbi + do_hmmforward != 1)) { 
	  printf("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but search opts are invalid\n", n);
	  DumpSearchOpts(si->search_opts[n]);
	  cm_Fail("This is fatal.");
	}
      }
      else if (si->stype[n] == SEARCH_WITH_CM) {
	if(si->smx[n] == NULL) cm_Fail("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_CM but smx[%d] is NULL\n", n, n);
	sum = do_hbanded + do_hmmscanbands + do_sums + do_toponly + do_null2 + do_rsearch + do_hmmviterbi + do_hmmforward;	
	if(sum != 0) {
	  printf("cm_ValidateSearchInfo(), round %d is SEARCH_WITH_CM but search opts are invalid\n", n);
	    DumpSearchOpts(si->search_opts[n]);
	    cm_Fail("This is fatal.");
	}
      }
      else cm_Fail("cm_ValidateSearchInfo(), round %d is neither type SEARCH_WITH_HMM, SEARCH_WITH_HYBRID, nor SEARCH_WITH_CM\n", n);
    }
    else { /* round n == si->nrounds */
      /* check final round, in which we're done filtering */
      if(si->stype[si->nrounds] == SEARCH_WITH_HYBRID) cm_Fail("cm_ValidateSearchInfo(), final round %d is SEARCH_WITH_HYBRID.", n);
      if(si->stype[si->nrounds] == SEARCH_WITH_CM && si->smx[n] == NULL)    cm_Fail("cm_ValidateSearchInfo(), final round %d is SEARCH_WITH_CM but smx is NULL.", n);
      if(si->stype[si->nrounds] == SEARCH_WITH_CM && si->smx[n] != cm->smx) cm_Fail("cm_ValidateSearchInfo(), final round %d is SEARCH_WITH_CM but smx != cm->smx.", n);
      if(si->hsi[si->nrounds]   != NULL)      cm_Fail("cm_ValidateSearchInfo(), final round hsi non-NULL.");
      if(do_hmmviterbi || do_hmmforward) { /* searching with only an HMM */
	if(do_hmmviterbi + do_hmmforward != 1) cm_Fail("cm_ValidateSearchInfo(), final round %d specifies HMM viterbi and HMM forward.\n");
	if(do_inside)                          cm_Fail("cm_ValidateSearchInfo(), final round %d specifies HMM viterbi or HMM forward but also Inside.\n");
      }
    }
  }
  ESL_DPRINTF1(("SearchInfo validated.\n"));
  return;
}


/* Function: cm_UpdateSearchInfoCutoff()
 * Date:     EPN, Tue Nov 27 13:43:21 2007
 *
 * Purpose:  Update the cutoff value for a specified round of filtering
 *            
 * Returns:  void, dies if some error
 */
void
cm_UpdateSearchInfoCutoff(CM_t *cm, int nround, int cutoff_type, float cutoff)
{
  if(cm->si == NULL)           cm_Fail("cm_UpdateSearchInfoCutoff(), cm->si is NULL.");
  if(nround > cm->si->nrounds) cm_Fail("cm_UpdateSearchInfoCutoff(), requested round %d is > cm->si->nrounds\n", nround, cm->si->nrounds);
  cm->si->cutoff_type[nround] = cutoff_type;
  cm->si->cutoff[nround]      = cutoff;
  return;
}

/* Function: DumpSearchOpts()
 * Date:     EPN, Tue Nov 27 09:02:21 2007
 *
 * Purpose:  Print search options that are turned on in a search_opts integer.
 *            
 * Returns:  void.
 */
void
DumpSearchOpts(int search_opts)
{
  if(search_opts & CM_SEARCH_NOQDB)        printf("\tCM_SEARCH_NOQDB\n");
  if(search_opts & CM_SEARCH_HBANDED)      printf("\tCM_SEARCH_HBANDED\n");
  if(search_opts & CM_SEARCH_HMMSCANBANDS) printf("\tCM_SEARCH_HMMSCANBANDS\n");
  if(search_opts & CM_SEARCH_SUMS)         printf("\tCM_SEARCH_SUMS\n");
  if(search_opts & CM_SEARCH_INSIDE)       printf("\tCM_SEARCH_INSIDE\n");
  if(search_opts & CM_SEARCH_TOPONLY)      printf("\tCM_SEARCH_TOPONLY\n");
  if(search_opts & CM_SEARCH_NOALIGN)      printf("\tCM_SEARCH_NOALIGN\n");
  if(search_opts & CM_SEARCH_NULL2)        printf("\tCM_SEARCH_NULL2\n");
  if(search_opts & CM_SEARCH_RSEARCH)      printf("\tCM_SEARCH_RSEARCH\n");
  if(search_opts & CM_SEARCH_CMGREEDY)     printf("\tCM_SEARCH_CMGREEDY\n");
  if(search_opts & CM_SEARCH_HMMGREEDY)    printf("\tCM_SEARCH_HMMGREEDY\n");
  if(search_opts & CM_SEARCH_HMMVITERBI)   printf("\tCM_SEARCH_HMMVITERBI\n");
  if(search_opts & CM_SEARCH_HMMFORWARD)   printf("\tCM_SEARCH_HMMFORWARD\n");
  return;
}

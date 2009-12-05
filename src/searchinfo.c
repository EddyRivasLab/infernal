/* searchinfo.c
 * SearchInfo_t implementations: information for CM/CP9 
 * filters and scans.
 *
 * EPN, Tue Nov 27 08:42:08 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_exponential.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"


/* Function: CreateSearchInfo()
 * Date:     EPN, Tue Nov 27 12:57:24 2007
 *
 * Purpose:  Allocate and initialize a search info object that
 *           specifies that no filtering should be done. 
 *            
 * Returns:  cm->si points to a new SearchInfo_t object
 *           eslOK on success, dies immediately on some error
 */
int
CreateSearchInfo(CM_t *cm, int cutoff_type, float sc_cutoff, float e_cutoff)
{
  int status;
  int use_hmmonly;

  if(cm->si != NULL)  cm_Fail("CreateSearchInfo(), the cm already points to a SearchInfo_t object.\n");
  use_hmmonly = ((cm->search_opts & CM_SEARCH_HMMVITERBI) ||  (cm->search_opts & CM_SEARCH_HMMFORWARD)) ? TRUE : FALSE;

  SearchInfo_t *si;
  ESL_ALLOC(si, sizeof(SearchInfo_t));
  
  si->nrounds = 0;
  ESL_ALLOC(si->search_opts, sizeof(int)                * (si->nrounds+1));
  ESL_ALLOC(si->cutoff_type, sizeof(int)                * (si->nrounds+1));
  ESL_ALLOC(si->sc_cutoff,   sizeof(float)              * (si->nrounds+1));
  ESL_ALLOC(si->e_cutoff,    sizeof(float)              * (si->nrounds+1));
  ESL_ALLOC(si->stype,       sizeof(int)                * (si->nrounds+1));
  ESL_ALLOC(si->smx,         sizeof(ScanMatrix_t *)     * (si->nrounds+1));
  ESL_ALLOC(si->hsi,         sizeof(HybridScanInfo_t *) * (si->nrounds+1));

  si->search_opts[0] = cm->search_opts;
  si->cutoff_type[0] = cutoff_type;
  si->sc_cutoff[0]   = sc_cutoff;
  si->e_cutoff[0]    = e_cutoff;
  si->stype[0]       = use_hmmonly ? SEARCH_WITH_HMM : SEARCH_WITH_CM;
  si->smx[0]         = cm->smx;  /* could be NULL */
  si->hsi[0]         = NULL;

  cm->si = si;
  return eslOK;

 ERROR:
  cm_Fail("CreateSearchInfo(), memory allocate error.");
  return status; /* NEVERREACHED */
}  


/* Function: AddFilterToSearchInfo()
 * Date:     EPN, Tue Nov 27 13:00:23 2007
 *
 * Purpose:  Add a filter as the 1st round filter for CM <cm>.
 *           A new SearchInfo_t object <fi> is created, and the existing
 *           information from cm->si is copied into it. cm->si is then
 *           freed and cm->si is set to point at fi.            
 * 
 * Returns:  cm->si points to a new SearchInfo_t object
 *           eslOK on success, dies immediately on some error
 */
int
AddFilterToSearchInfo(CM_t *cm, int cyk_filter, int inside_filter, int viterbi_filter, int forward_filter, int hybrid_filter, 
		      ScanMatrix_t *smx, HybridScanInfo_t *hsi, int cutoff_type, float sc_cutoff, float e_cutoff, int do_null3)
{
  int status;
  int n;
  int orig_nrounds;

  if(cm->si == NULL)                 cm_Fail("AddFilterToSearchInfo(), the cm does not point to a SearchInfo_t object.\n");
  if((cyk_filter + inside_filter + viterbi_filter + forward_filter + hybrid_filter) != 1)
    cm_Fail("AddFilterToSearchInfo(), cyk_filter: %d\ninside_filter: %d\nviterbi_filter: %d\nforward_filter: %d\nhybrid_filter: %d. Exactly 1 of these must be 1, the rest 0s.\n", cyk_filter, inside_filter, viterbi_filter, forward_filter, hybrid_filter);
  if(cyk_filter    && smx == NULL) cm_Fail("AddFilterToSearchInfo(), cyk_filter: %d but smx == NULL\n", cyk_filter);
  if(inside_filter && smx == NULL) cm_Fail("AddFilterToSearchInfo(), inside_filter: %d but smx == NULL\n", inside_filter);
  if(hybrid_filter && hsi == NULL) cm_Fail("AddFilterToSearchInfo(), hybrid_filter: %d but hsi == NULL\n", hybrid_filter);

  orig_nrounds = cm->si->nrounds;

  /* allocate new si object, with 1 more round than cm->si, and set round 0 */
  SearchInfo_t *si;
  ESL_ALLOC(si, sizeof(SearchInfo_t));
  si->nrounds = orig_nrounds+1;
  ESL_ALLOC(si->search_opts, sizeof(int)   * (si->nrounds+1));
  ESL_ALLOC(si->cutoff_type, sizeof(int)   * (si->nrounds+1));
  ESL_ALLOC(si->sc_cutoff,   sizeof(float) * (si->nrounds+1));
  ESL_ALLOC(si->e_cutoff,    sizeof(float) * (si->nrounds+1));
  ESL_ALLOC(si->stype,       sizeof(int)   * (si->nrounds+1));
  ESL_ALLOC(si->smx,         sizeof(ScanMatrix_t     *) * (si->nrounds+1));
  ESL_ALLOC(si->hsi,         sizeof(HybridScanInfo_t *) * (si->nrounds+1));

  si->search_opts[0] = 0;
  si->search_opts[0] |= CM_SEARCH_NOALIGN;
  if(cyk_filter)     ;/* do nothing, CYK is default */
  if(inside_filter)  si->search_opts[0] |= CM_SEARCH_INSIDE;
  if(viterbi_filter) si->search_opts[0] |= CM_SEARCH_HMMVITERBI;
  if(forward_filter) si->search_opts[0] |= CM_SEARCH_HMMFORWARD;
  if(sc_cutoff < -eslSMALLX1) { /* if we're asking to return negative scores, turn on the greedy hit resolution algorithm (that's the only way we can return negative scoring hits */
    if(viterbi_filter || forward_filter) si->search_opts[0] |= CM_SEARCH_HMMGREEDY;
    if(cyk_filter     || inside_filter)  si->search_opts[0] |= CM_SEARCH_CMGREEDY;
  }
  else { /* turn greedy options off, (they may already be off) */
    if(cm->si->stype[0] == SEARCH_WITH_HMM) si->search_opts[0] &= ~CM_SEARCH_HMMGREEDY;
    if(cm->si->stype[0] == SEARCH_WITH_CM)  si->search_opts[0] &= ~CM_SEARCH_CMGREEDY;
  }
  
  si->cutoff_type[0] = cutoff_type;
  si->sc_cutoff[0]   = sc_cutoff;
  si->e_cutoff[0]    = e_cutoff;
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
  if(do_null3) si->search_opts[0] |= CM_SEARCH_NULL3;
       
  /* copy existing information for other rounds from old cm->si */
  for(n = 0; n <= orig_nrounds; n++) { 
    si->search_opts[(n+1)] = cm->si->search_opts[n];
    si->cutoff_type[(n+1)] = cm->si->cutoff_type[n];
    si->sc_cutoff[(n+1)]   = cm->si->sc_cutoff[n];
    si->e_cutoff[(n+1)]    = cm->si->e_cutoff[n];
    si->stype[(n+1)]       = cm->si->stype[n];
    /* and copy the ptr to smx and hsi */
    si->smx[(n+1)] = cm->si->smx[n];
    si->hsi[(n+1)] = cm->si->hsi[n];
  }    

  /* free old cm->si, but be careful not to free smx[] and hsi[], we're still pointing to those */
  free(cm->si->search_opts);
  free(cm->si->cutoff_type);
  free(cm->si->sc_cutoff);
  free(cm->si->e_cutoff);
  free(cm->si->stype);
  free(cm->si->smx);
  free(cm->si->hsi);
  free(cm->si);

  cm->si = si;
  return eslOK;

 ERROR:
  cm_Fail("AddFilterToSearchInfo(), memory allocate error.");
  return status; /* NEVERREACHED */
}  

/* Function: FreeSearchInfo()
 * Date:     EPN, Tue Nov 27 08:48:49 2007
 *
 * Purpose:  Free a SearchInfo_t object corresponding to 
 *           CM <cm>.
 *            
 * Returns:  void
 */
void
FreeSearchInfo(SearchInfo_t *si, CM_t *cm)
{
  int n;

  for(n = 0; n <  si->nrounds; n++) if(si->smx[n] != NULL) cm_FreeScanMatrix(cm, si->smx[n]);
  /* we don't free si->smx[nrounds] b/c it == cm->smx */
  for(n = 0; n <= si->nrounds; n++) if(si->hsi[n] != NULL) cm_FreeHybridScanInfo(si->hsi[n], cm); 
  free(si->search_opts);
  free(si->cutoff_type);
  free(si->sc_cutoff);
  free(si->e_cutoff);
  free(si->stype);
  free(si->smx);
  free(si->hsi);

  free(si);
  return;
}

/* Function: DumpSearchInfo()
 * Date:     EPN, Tue Nov 27 08:50:54 2007
 *
 * Purpose:  Dump a CM's search info (except scan matrix and hybrid scan info) to stdout. 
 *            
 * Returns:  void.
 */
void
DumpSearchInfo(SearchInfo_t *si)
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
    if(si->cutoff_type[n] == SCORE_CUTOFF) printf("\tcutoff     : %10.4f bits\n",    si->sc_cutoff[n]);
    else                                   printf("\tcutoff     : %10.4f E-value\n", si->e_cutoff[n]);
    if(si->hsi[n] != NULL) { 
      printf("\tHybrid info:\n");
      printf("\t\tNumber of sub CM roots: %d\n", si->hsi[n]->n_v_roots);
      for(v = 0; v < si->hsi[n]->cm_M; v++) 
	if(si->hsi[n]->v_isroot[v]) printf("\t\tstate %d is a root\n", v);
    }
  }
  return;
}


/* Function: ValidateSearchInfo()
 * Date:     EPN, Tue Nov 27 09:25:10 2007
 *
 * Purpose:  Validate a Search Info <si> object for CM <cm>.
 *            
 * Returns:  void.
 */
void
ValidateSearchInfo(CM_t *cm, SearchInfo_t *si)
{
  int n, sum;
  int do_noqdb;
  int do_hbanded;
  int do_hmmalnbands;
  int do_sums;
  int do_inside;
  int do_noalign;
  int do_rsearch;
  int do_cmgreedy;
  int do_hmmgreedy;
  int do_hmmviterbi;
  int do_hmmforward;

  for(n = 0; n <= si->nrounds; n++) { 
    do_noqdb       = (si->search_opts[n] & CM_SEARCH_NOQDB)        ? TRUE : FALSE;
    do_hbanded     = (si->search_opts[n] & CM_SEARCH_HBANDED)      ? TRUE : FALSE;
    do_hmmalnbands = (si->search_opts[n] & CM_SEARCH_HMMALNBANDS)  ? TRUE : FALSE;
    do_sums        = (si->search_opts[n] & CM_SEARCH_SUMS)         ? TRUE : FALSE;
    do_inside      = (si->search_opts[n] & CM_SEARCH_INSIDE)       ? TRUE : FALSE;
    do_noalign     = (si->search_opts[n] & CM_SEARCH_NOALIGN)      ? TRUE : FALSE;
    do_rsearch     = (si->search_opts[n] & CM_SEARCH_RSEARCH)      ? TRUE : FALSE;
    do_cmgreedy    = (si->search_opts[n] & CM_SEARCH_CMGREEDY)     ? TRUE : FALSE;
    do_hmmgreedy   = (si->search_opts[n] & CM_SEARCH_HMMGREEDY)    ? TRUE : FALSE;
    do_hmmviterbi  = (si->search_opts[n] & CM_SEARCH_HMMVITERBI)   ? TRUE : FALSE;
    do_hmmforward  = (si->search_opts[n] & CM_SEARCH_HMMFORWARD)   ? TRUE : FALSE;

    if(n < si->nrounds) { 
      if(!do_noalign) cm_Fail("ValidateSearchInfo(), round %d has CM_SEARCH_NOALIGN flag down.\n", n);
      if(si->stype[n] == SEARCH_WITH_HMM) {
	sum = do_noqdb + do_hbanded + do_hmmalnbands + do_sums + do_inside + do_rsearch + do_cmgreedy;
	if(sum != 0 || (do_hmmviterbi + do_hmmforward != 1)) { 
	  printf("ValidateSearchInfo(), round %d is SEARCH_WITH_HMM but search opts are invalid\n", n);
	  DumpSearchOpts(si->search_opts[n]);
	  cm_Fail("This is fatal.");
	}
	if(si->smx[n] != NULL) cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_HMM but smx[%d] is non-NULL\n", n, n);
	if(si->hsi[n] != NULL) cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_HMM but hsi[%d] is non-NULL\n", n, n);
      }
      else if (si->stype[n] == SEARCH_WITH_HYBRID) {
	if(si->hsi[n] == NULL)      cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but hsi[%d] is NULL\n", n, n);
	if(si->hsi[n]->smx == NULL) cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but hsi[%d]->smx is NULL\n", n, n);
	if(si->smx[n] != NULL)      cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but smx[%d] is not NULL\n", n, n);
	if(si->hsi[n]->v_isroot[0]) cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID and hsi->vi_isroot[0] is TRUE, this shouldn't happen, we might as well filter with a SEARCH_WITH_CM filter.");
	sum = do_hbanded + do_hmmalnbands + do_sums + do_rsearch + do_cmgreedy;	
	if(sum != 0 || (do_hmmviterbi + do_hmmforward != 1)) { 
	  printf("ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but search opts are invalid\n", n);
	  DumpSearchOpts(si->search_opts[n]);
	  cm_Fail("This is fatal.");
	}
      }
      else if (si->stype[n] == SEARCH_WITH_CM) {
	if(si->smx[n] == NULL) cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_CM but smx[%d] is NULL\n", n, n);
	sum = do_hbanded + do_hmmalnbands + do_sums + do_rsearch + do_hmmviterbi + do_hmmforward;	
	if(sum != 0) {
	  printf("ValidateSearchInfo(), round %d is SEARCH_WITH_CM but search opts are invalid\n", n);
	    DumpSearchOpts(si->search_opts[n]);
	    cm_Fail("This is fatal.");
	}
      }
      else cm_Fail("ValidateSearchInfo(), round %d is neither type SEARCH_WITH_HMM, SEARCH_WITH_HYBRID, nor SEARCH_WITH_CM\n", n);
    }
    else { /* round n == si->nrounds */
      /* check final round, in which we're done filtering */
      if(si->stype[si->nrounds] == SEARCH_WITH_HYBRID) cm_Fail("ValidateSearchInfo(), final round %d is SEARCH_WITH_HYBRID.", n);
      if(si->stype[si->nrounds] == SEARCH_WITH_CM && si->smx[n] == NULL)    cm_Fail("ValidateSearchInfo(), final round %d is SEARCH_WITH_CM but smx is NULL.", n);
      if(si->stype[si->nrounds] == SEARCH_WITH_CM && si->smx[n] != cm->smx) cm_Fail("ValidateSearchInfo(), final round %d is SEARCH_WITH_CM but smx != cm->smx.", n);
      if(si->hsi[si->nrounds]   != NULL)      cm_Fail("ValidateSearchInfo(), final round hsi non-NULL.");
      if(do_hmmviterbi || do_hmmforward) { /* searching with only an HMM */
	if(do_hmmviterbi + do_hmmforward != 1) cm_Fail("ValidateSearchInfo(), final round %d specifies HMM viterbi and HMM forward.\n");
	if(do_inside)                          cm_Fail("ValidateSearchInfo(), final round %d specifies HMM viterbi or HMM forward but also Inside.\n");
      }
    }
  }
  ESL_DPRINTF1(("SearchInfo validated.\n"));
  return;
}


/* Function: UpdateSearchInfoCutoff()
 * Date:     EPN, Tue Nov 27 13:43:21 2007
 *
 * Purpose:  Update the cutoff value for a specified round of filtering
 *            
 * Returns:  void, dies if some error
 */
void
UpdateSearchInfoCutoff(CM_t *cm, int nround, int cutoff_type, float sc_cutoff, float e_cutoff)
{
  if(cutoff_type == E_CUTOFF)  
    if(! (cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("UpdateSearchInfoCutoff(), cm->si is NULL.");
  if(cm->si == NULL)           cm_Fail("UpdateSearchInfoCutoff(), cm->si is NULL.");
  if(nround > cm->si->nrounds) cm_Fail("UpdateSearchInfoCutoff(), requested round %d is > cm->si->nrounds\n", nround, cm->si->nrounds);
  cm->si->cutoff_type[nround] = cutoff_type;
  cm->si->sc_cutoff[nround]   = sc_cutoff;
  cm->si->e_cutoff[nround]    = e_cutoff;
  if(sc_cutoff < -eslSMALLX1) { /* if we're asking to return negative scores, turn on the greedy hit resolution algorithm (that's the only way we can return negative scoring hits */
    if(cm->si->stype[nround] == SEARCH_WITH_HMM) cm->si->search_opts[nround] |= CM_SEARCH_HMMGREEDY;
    if(cm->si->stype[nround] == SEARCH_WITH_CM)  cm->si->search_opts[nround] |= CM_SEARCH_CMGREEDY;
  }
  else { /* turn greedy options off, (they may already be off) */
    if(cm->si->stype[nround] == SEARCH_WITH_HMM) cm->si->search_opts[nround] &= ~CM_SEARCH_HMMGREEDY;
    if(cm->si->stype[nround] == SEARCH_WITH_CM)  cm->si->search_opts[nround] &= ~CM_SEARCH_CMGREEDY;
  }

  return;
}

/* Function: UpdateSearchInfoForExpMode
 * Date:     EPN, Thu Jan 24 11:57:20 2008
 * Purpose:  Given a exp tail mode and a search round <round>, update cm->si
 *           SearchInfo_t for that round.
 *           that exp tail mode. 
 *
 *           0. EXP_CM_GC : !cm->si->search_opts[round] & CM_SEARCH_INSIDE
 *           1. EXP_CM_GI : !cm->si->search_opts[round] & CM_SEARCH_INSIDE
 *           4. EXP_CP9_GV:  cm->si->search_opts[round] & CM_SEARCH_HMMVITERBI
 *                          !cm->si->search_opts[round] & CM_SEARCH_HMMFORWARD
 *           5. EXP_CP9_GF:  cm->si->search_opts[round] & CM_SEARCH_HMMVITERBI
 *                          !cm->si->search_opts[round] & CM_SEARCH_HMMFORWARD
 *           3. EXP_CM_LC :  cm->si->search_opts[round] & CM_SEARCH_INSIDE
 *           2. EXP_CM_LI :  cm->si->search_opts[round] & CM_SEARCH_INSIDE
 *           6. EXP_CP9_LV: !cm->si->search_opts[round] & CM_SEARCH_HMMVITERBI
 *                           cm->si->search_opts[round] & CM_SEARCH_HMMFORWARD
 *           7. EXP_CP9_LF: !cm->si->search_opts[round] & CM_SEARCH_HMMVITERBI
 *                           cm->si->search_opts[round] & CM_SEARCH_HMMFORWARD
 * 
 * Args:
 *           CM           - the covariance model
 *           exp_mode     - the mode 0..EXP_NMODES-1
 */
void
UpdateSearchInfoForExpMode(CM_t *cm, int round, int exp_mode)
{
  if(cm->si == NULL)           cm_Fail("UpdateSearchInfoForExpMode(), cm->si is NULL.");
  if(round > cm->si->nrounds)  cm_Fail("UpdateSearchInfoForExpMode(), requested round %d is > cm->si->nrounds\n", round, cm->si->nrounds);

  ESL_DASSERT1((exp_mode >= 0 && exp_mode < EXP_NMODES));

  switch (exp_mode) {
  case EXP_CM_GC: 
  case EXP_CM_LC: /* CYK, local or glocal */
    cm->si->search_opts[round] &= ~CM_SEARCH_INSIDE;
    cm->si->search_opts[round] &= ~CM_SEARCH_HMMVITERBI;
    cm->si->search_opts[round] &= ~CM_SEARCH_HMMFORWARD;
    cm->si->stype[round] = SEARCH_WITH_CM;
    break;
  case EXP_CM_GI: 
  case EXP_CM_LI: /* Inside, local or glocal */
    cm->si->search_opts[round] |= CM_SEARCH_INSIDE;
    cm->si->search_opts[round] &= ~CM_SEARCH_HMMVITERBI;
    cm->si->search_opts[round] &= ~CM_SEARCH_HMMFORWARD;
    cm->si->stype[round] = SEARCH_WITH_CM;
    break;
  case EXP_CP9_GV: 
  case EXP_CP9_LV: /* Viterbi, local or glocal */
    cm->si->search_opts[round] &= ~CM_SEARCH_INSIDE;
    cm->si->search_opts[round] |= CM_SEARCH_HMMVITERBI;
    cm->si->search_opts[round] &= ~CM_SEARCH_HMMFORWARD;
    cm->si->stype[round] = SEARCH_WITH_HMM;
    break;
  case EXP_CP9_GF: 
  case EXP_CP9_LF: /* Forward, local or glocal */
    cm->si->search_opts[round] &= ~CM_SEARCH_INSIDE;
    cm->si->search_opts[round] &= ~CM_SEARCH_HMMVITERBI;
    cm->si->search_opts[round] |= CM_SEARCH_HMMFORWARD;
    cm->si->stype[round] = SEARCH_WITH_HMM;
    break;
  default: 
    cm_Fail("UpdateSearchInfoForExpMode(): bogus exp_mode: %d\n", exp_mode);
  }
  return;
}


/* Function: UpdateSearchInfoForNewSMX()
 * Date:     EPN, Thu Feb 14 11:38:23 2008
 * Purpose:  Point the ScanMatrix_t for the final round of searching to the 
 *           cm->smx.
 *
 * 
 *           CM           - the covariance model
 */
void
UpdateSearchInfoForNewSMX(CM_t *cm)
{
  if(cm->si == NULL)                cm_Fail("UpdateSearchInfoForNewSMX(), cm->si is NULL.");
  if(cm->smx == NULL)               cm_Fail("UpdateSearchInfoForNewSMX(), cm->smx is NULL.");
  if(!(cm->flags & CMH_SCANMATRIX)) cm_Fail("UpdateSearchInfoForNewSMX(), CMH_SCAN_MATRIX flag is down.");

  cm->si->smx[cm->si->nrounds] = cm->smx;
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
  if(search_opts & CM_SEARCH_HMMALNBANDS)  printf("\tCM_SEARCH_HMMALNBANDS\n");
  if(search_opts & CM_SEARCH_SUMS)         printf("\tCM_SEARCH_SUMS\n");
  if(search_opts & CM_SEARCH_INSIDE)       printf("\tCM_SEARCH_INSIDE\n");
  if(search_opts & CM_SEARCH_NOALIGN)      printf("\tCM_SEARCH_NOALIGN\n");
  if(search_opts & CM_SEARCH_RSEARCH)      printf("\tCM_SEARCH_RSEARCH\n");
  if(search_opts & CM_SEARCH_CMGREEDY)     printf("\tCM_SEARCH_CMGREEDY\n");
  if(search_opts & CM_SEARCH_HMMGREEDY)    printf("\tCM_SEARCH_HMMGREEDY\n");
  if(search_opts & CM_SEARCH_HMMVITERBI)   printf("\tCM_SEARCH_HMMVITERBI\n");
  if(search_opts & CM_SEARCH_HMMFORWARD)   printf("\tCM_SEARCH_HMMFORWARD\n");
  return;
}


/*
 * Function: CreateResults ()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Creates a results type of specified size
 */
search_results_t *CreateResults (int size) {
  int status;
  search_results_t *results;
  int i;

  if(size == 0) return NULL;

  ESL_ALLOC(results, sizeof(search_results_t));
  results->num_results = 0;
  results->num_allocated = size;
  ESL_ALLOC(results->data, sizeof(search_result_node_t)*size);
  for(i = 0; i < size; i++) {
    results->data[i].start  = -1;    
    results->data[i].stop   = -1;    
    results->data[i].bestr  = -1;    
    results->data[i].score  = IMPOSSIBLE;    
    results->data[i].tr     = NULL;
    results->data[i].pcode  = NULL;
  }
  return (results);
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

/* Function: ExpandResults ()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Expands a results structure by specified amount
 */
void ExpandResults (search_results_t *results, int additional) {
  int status;
  void *tmp;
  int i;
  ESL_RALLOC(results->data, tmp, sizeof(search_result_node_t) * (results->num_allocated+additional));

  for(i = results->num_allocated; i < (results->num_allocated + additional); i++) {
    results->data[i].start  = -1;    
    results->data[i].stop   = -1;    
    results->data[i].bestr  = -1;    
    results->data[i].score  = IMPOSSIBLE;    
    results->data[i].tr     = NULL;
    results->data[i].pcode  = NULL;
  }

  results->num_allocated+=additional;
  return;
 ERROR:
  cm_Fail("Memory reallocation error.");
}

/* Function: AppendResults()
 * Date:     EPN, Wed Aug 29 08:58:28 2007
 *
 * Purpose:  Add result nodes from one results structure onto
 *           another by copying data and manipulating pointers. 
 *           Originally written to add results returned from 
 *           an MPI worker to a growing 'master' results structure 
 *           in the MPI master. 
 *
 *           The search_results_node_t's (results->data)
 *           must have their start, stop, bestr, and score data
 *           copied because the results->data is not a set 
 *           of pointers (the whole reason for this is so 
 *           we can call quicksort() on the hits, which requires
 *           we have a fixed distance between them in memory).
 *           The parsetree, however, can have it's pointer 
 *           simply switched.
 *
 *           i0 is an offset in the sequence, (i0-1) is added to
 *           data[i].start and data[i].stop for hits in src_results. 
 *           i0 == 1 means no offset. i0 != 1 usually useful for hits 
 *           returned by MPI workers who were searching a database
 *           subsequence.
 *
 * Note:     Because the dest_results->data now points to some of the
 *           src_results->data, be careful not to free src_results with
 *           FreeResults() if you don't want to lose dest_results.
 */
void AppendResults (search_results_t *src_results, search_results_t *dest_results, int i0) {
  int i, ip;
  for(i = 0; i < src_results->num_results; i++) 
    {
      ip = dest_results->num_results;
      ReportHit (src_results->data[i].start+i0-1, src_results->data[i].stop+i0-1, 
		 src_results->data[i].bestr,      src_results->data[i].score,
		 dest_results);
      if(src_results->data[i].tr != NULL)
	(*dest_results).data[ip].tr = (*src_results).data[i].tr;
      if(src_results->data[i].pcode != NULL)
	(*dest_results).data[ip].pcode = (*src_results).data[i].pcode;
    }
  return;
}

/* Function: FreeResults()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Frees a results structure
 */
void FreeResults (search_results_t *r) {
  int i;
  if (r != NULL) {
    for (i=0; i < r->num_allocated; i++) {
      if (r->data[i].tr     != NULL) FreeParsetree(r->data[i].tr);
      if (r->data[i].pcode  != NULL) free(r->data[i].pcode);
    }
    free (r->data);
    free(r);
  }
}


/* Function: CompareResultsByScore()
 * Date:     RJK, Wed Apr 10, 2002 [St. Louis]
 * Purpose:  Compares two search_result_node_ts based on score and returns -1
 *           if first is higher score than second, 0 if equal, 1 if first
 *           score is lower.  This results in sorting by score, highest
 *           first.
 */
int CompareResultsByScore (const void *a_void, const void *b_void) {
  search_result_node_t *a, *b;
 
  a = (search_result_node_t *)a_void;
  b = (search_result_node_t *)b_void;

  if (a->score < b->score)
    return (1);
  else if (a->score > b->score)
    return (-1);
  else if (a->start < b->start)
    return (1);
  else if (a->start > b->start)
    return (-1);
  else
    return (0);
}

/* Function: CompareResultsByEndPoint()
 * Date:     RJK, Wed Apr 10, 2002 [St. Louis]
 * Purpose:  Compares two search_result_node_ts based on end point and returns -1
 *           if first is higher score than second, 0 if equal, 1 if first
 *           score is lower.  This results in sorting by end points j, highest
 *           first.
 */
int CompareResultsByEndPoint (const void *a_void, const void *b_void) {
  search_result_node_t *a, *b;
 
  a = (search_result_node_t *)a_void;
  b = (search_result_node_t *)b_void;

  if      (a->stop  < b->stop)  return ( 1);
  else if (a->stop  > b->stop)  return (-1);
  else if (a->start < b->start) return ( 1);
  else if (a->start > b->start) return (-1);
  else                          return ( 0);
}

/* Function: UpdateHitScoresWithNull2Or3()
 * Date:     EPN, Tue Apr  1 06:32:34 2008
 * Purpose:  Update scores for each hit based on NULL2 or NULL3 correction.
 *           Remove any hits that are now below <sc_cutoff>. To not remove
 *           any hits, pass <sc_cutoff> == IMPOSSIBLE.
 * 
 * Args:    
 *           cm      - the covariance model
 *           errbuf  - for printing error messages
 *           si      - SearchInfo, relevant round is final one, si->nrounds
 *           results - the hits data structure
 *           dsq     - digitized seq hits lie within, needed to determine gc content
 *           first_result   - index of first result in search_results to align (if dsq_mode)
 *           sc_cutoff - we'll remove hits that have new bit scores < <sc_cutoff> (IMPOSSIBLE => don't remove any hits)
 *           do_null2 - TRUE to do NULL2 score correction
 *           do_null3 - TRUE to do NULL3 score correction
 *           sort_by_score    - TRUE to sort list by score at the end of the function IFF we remove any hits
 *           sort_by_endpoint - TRUE to sort list by end point at the end of the function IFF we remove any hits
 */
int UpdateHitScoresWithNull2Or3(CM_t *cm, char *errbuf, SearchInfo_t *si, search_results_t *results, ESL_DSQ *dsq, int first_result, float sc_cutoff, int do_null2, int do_null3, int sort_by_score, int sort_by_endpoint)
{
  int status;
  int i, x;
  search_result_node_t swap;
  float corr_sc;             /* score correction */
  int orig_num_results; /* number of results when function entered */

  /*printf("in UpdateHitScoresWithNull2Or3()\n");*/

  /* Check contract */
  if(dsq == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), dsq is NULL\n");
  if(si == NULL)      ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), si is NULL\n");
  if(si->stype[si->nrounds] != SEARCH_WITH_HMM && si->stype[si->nrounds] != SEARCH_WITH_CM) ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), final search round is neither SEARCH_WITH_HMM nor SEARCH_WITH_CM.\n");
  if(do_null2 && do_null3)    ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), do_null2 is TRUE and do_null3 is TRUE. (only one should be true)\n");
  if(!do_null2 && !do_null3)  ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), do_null2 is FALSE and do_null3 is FALSE. (one should be true)\n");
  if(sort_by_score && sort_by_endpoint)   ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), sort_by_score is TRUE and sort_by_endpoint is TRUE (only one should be true).\n");
  if(!sort_by_score && !sort_by_endpoint) ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), sort_by_score is FALSE and sort_by_endpoint is FALSE (one should be true).\n");
  if(first_result > results->num_results) ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), first_result %d > results->num_results %d\n", first_result, results->num_results);

  if (results == NULL) return eslOK;

  orig_num_results = results->num_results;
  for (i=first_result; i<results->num_results; i++) {
    if(results->data[i].tr == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateHitScoresWithNull2Or3(), parsetree for hit i: %d is NULL.\n", i);
    if(do_null2) { 
      if((status = ParsetreeScoreCorrectionNull2(cm, errbuf, results->data[i].tr, dsq, results->data[i].start, &corr_sc)) != eslOK) return status;
    }
    else if(do_null3) { 
      /*printf("\norig score: %.3f\n", results->data[i].score);*/
      if((status = ParsetreeScoreCorrectionNull3(cm, errbuf, results->data[i].tr, dsq, results->data[i].start, &corr_sc)) != eslOK) return status;
    }
    results->data[i].score -= corr_sc;
    /*printf("new  score: %.3f\n", results->data[i].score);*/
    if(results->data[i].score < sc_cutoff) { 
      results->data[i].start = -1;
    }
  }

  for (x=0; x<results->num_results; x++) {
    while (results->num_results > 0 && 
	   results->data[results->num_results-1].start == -1)
      results->num_results--;
    if (x<results->num_results && results->data[x].start == -1) {
      swap = results->data[x];
      results->data[x] = results->data[results->num_results-1];
      results->data[results->num_results-1] = swap;
      results->num_results--;
    }
  }
  while (results->num_results > 0 && results->data[results->num_results-1].start == -1)
    results->num_results--;

  if(results->num_results != orig_num_results) { 
    if(sort_by_score)    SortResultsByScore(results);
    if(sort_by_endpoint) SortResultsByEndPoint(results);
  }
  return eslOK;
}  


/* Function: RemoveOverlappingHits ()
 * Date:     EPN, Tue Apr  3 14:36:38 2007
 * Plucked verbatim out of RSEARCH.
 * RSEARCH date: RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Purpose:  Given a list of hits, removes overlapping hits to produce
 * a list consisting of at most one hit covering each nucleotide in the
 * sequence.  Works as follows:
 * 1.  quicksort hits 
 * 2.  For each hit, sees if any nucleotide covered yet
 *     If yes, remove hit
 *     If no, mark each nt as covered
 * 
 * Args: 
 *        i0    - first position of subsequence results work for (1 for full seq)
 *        j0    - last  position of subsequence results work for (L for full seq)
 */
void RemoveOverlappingHits (search_results_t *results, int i0, int j0)
{
  int status;
  char *covered_yet;
  int x,y;
  int covered;
  int L;
  int yp;          /* offset position, yp = y-i0+1 */
  search_result_node_t swap;

  if (results == NULL)
    return;

  if (results->num_results == 0)
    return;

  L = j0-i0+1;
  ESL_ALLOC(covered_yet, sizeof(char)*(L+1));
  for (x=0; x<=L; x++)
    covered_yet[x] = 0;

  SortResultsByScore(results);

  for (x=0; x<results->num_results; x++) {
    covered = 0;
    for (y=results->data[x].start; y<=results->data[x].stop && !covered; y++) {
      {
	yp = y-i0+1; 
	if(yp <= 0 || yp > L) { 
	  printf("whoa\n");
	}
	assert(yp > 0 && yp <= L);
	if (covered_yet[yp] != 0) {
	  covered = 1;
	} 
      }
    }
    if (covered == 1) {
      results->data[x].start = -1;        /* Flag -- remove later to keep sorted */
    } else {
      for (y=results->data[x].start; y<=results->data[x].stop; y++) {
	yp = y-i0+1; 
	covered_yet[yp] = 1;
      }
    }
  }
  free (covered_yet);

  for (x=0; x < results->num_results; x++) {
    while (results->num_results > 0 &&
	   results->data[results->num_results-1].start == -1)
      results->num_results--;
    if (x < results->num_results && results->data[x].start == -1) {
      swap = results->data[x];
      results->data[x] = results->data[results->num_results-1];
      results->data[results->num_results-1] = swap;
      results->num_results--;
    }
  }
  while (results->num_results > 0 &&
	 results->data[results->num_results-1].start == -1)
    results->num_results--;

  SortResultsByScore(results);
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
}

/* Function: RemoveHitsOverECutoff
 * Date:     RJK, Tue Oct 8, 2002 [St. Louis]
 * Purpose:  Given an E-value cutoff, lambdas, mus, a sequence, and
 *           a list of results, calculates GC content for each hit, 
 *           calculates E-value, and decides whether to keep hit or not.
 * 
 * Args:    
 *           cm      - the covariance model
 *           si      - SearchInfo, relevant round is <sround>
 *           sround  - round of search we're removing hits for
 *           results - the hits data structure
 *           dsq     - seq hits lie within, needed to determine gc content
 *           first_result   - index of first result in search_results to align (if dsq_mode)
 *           sort_by_score    - TRUE to sort list by score at the end of the function IFF we remove any hits
 *           sort_by_endpoint - TRUE to sort list by end point at the end of the function IFF we remove any hits
 *
 * Returns: eslOK on success, eslEINCOMPAT on contract violation
 */
int RemoveHitsOverECutoff (CM_t *cm, char *errbuf, SearchInfo_t *si, int sround, search_results_t *results, ESL_DSQ *dsq, int first_result, int sort_by_score, int sort_by_endpoint)
{
  int gc_comp;
  int i, x;
  search_result_node_t swap;
  float score_for_Eval; /* the score we'll determine the statistical signifance of. */
  int cm_exp_mode;      /* exp tail mode if we're using CM hits */
  int cp9_exp_mode;     /* exp tail mode if we're using HMM hits */
  int p;                /* relevant partition */
  ExpInfo_t **exp;      /* pointer to exp tail to use */
  float cutoff;         /* the max E-value we want to keep */
  int orig_num_results; /* number of results when function entered */

  /* Check contract */
  if(!(cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "remove_hits_over_e_cutoff(), but CM has no exp tail stats");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "remove_hits_over_e_cutoff(), dsq == NULL.");
  if(si == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "remove_hits_over_e_cutoff(), si == NULL.");
  if(si->stype[si->nrounds] != SEARCH_WITH_HMM && si->stype[si->nrounds] != SEARCH_WITH_CM) ESL_FAIL(eslEINCOMPAT, errbuf, "remove_hits_over_e_cutoff(), final search round is neither SEARCH_WITH_HMM nor SEARCH_WITH_CM.");
  if(sort_by_score == TRUE  && sort_by_endpoint == TRUE)  ESL_FAIL(eslEINCOMPAT, errbuf, "remove_hits_over_e_cutoff(), sort_by_score and sort_by_endpoint both TRUE.");
  if(sort_by_score == FALSE && sort_by_endpoint == FALSE) ESL_FAIL(eslEINCOMPAT, errbuf, "remove_hits_over_e_cutoff(), sort_by_score and sort_by_endpoint both FALSE (disallowed, but shouldn't be a problem to allow).");
  if(first_result > results->num_results) ESL_FAIL(eslEINCOMPAT, errbuf, "RemoveHitsOverECutoff(), first_result %d > results->num_results %d\n", first_result, results->num_results);
  
  if (results == NULL) return eslOK;

  /* Determine exp tail mode to use */
  orig_num_results = results->num_results;
  CM2ExpMode(cm, si->search_opts[sround], &cm_exp_mode, &cp9_exp_mode);
  exp = (si->stype[sround] == SEARCH_WITH_HMM) ? cm->stats->expAA[cp9_exp_mode] : cm->stats->expAA[cm_exp_mode];
  
  ESL_DASSERT1((si->cutoff_type[sround] == E_CUTOFF));
  cutoff = si->e_cutoff[sround];
  
  for (i=first_result; i<results->num_results; i++) {
    gc_comp = get_gc_comp (cm->abc, dsq, results->data[i].start, results->data[i].stop);
    p = cm->stats->gc2p[gc_comp];
    score_for_Eval = results->data[i].score;
    if (Score2E(score_for_Eval, exp[p]->mu_extrap, exp[p]->lambda, exp[p]->cur_eff_dbsize) > cutoff)  
      results->data[i].start = -1;
  }
  
  for (x=0; x<results->num_results; x++) {
    while (results->num_results > 0 && 
	   results->data[results->num_results-1].start == -1)
      results->num_results--;
    if (x<results->num_results && results->data[x].start == -1) {
      swap = results->data[x];
      results->data[x] = results->data[results->num_results-1];
      results->data[results->num_results-1] = swap;
      results->num_results--;
    }
  }
  while (results->num_results > 0 && results->data[results->num_results-1].start == -1)
    results->num_results--;
  /* only sort if we removed anything */
  if(results->num_results != orig_num_results) { 
    if(sort_by_score)    SortResultsByScore(results);
    if(sort_by_endpoint) SortResultsByEndPoint(results);
  }
  return eslOK;
}  

/* Function: SortResults()
 * Date:    RJK,  Sun Mar 31, 2002 [AA Flight 2869 LGA->STL]
 * Purpose: Given a results array, sorts it by score with a call to qsort
 *
 */
void SortResultsByScore (search_results_t *results) 
{
  qsort (results->data, results->num_results, sizeof(search_result_node_t), CompareResultsByScore);
}

/* Function: SortResultsByEndPoint()
 * Date:     EPN, Thu Jan 24 13:06:50 2008
 * Purpose:  Given a results array, sorts it by end point with a call to qsort
 *
 */
void SortResultsByEndPoint (search_results_t *results) 
{
  qsort (results->data, results->num_results, sizeof(search_result_node_t), CompareResultsByEndPoint);
}

/*
 * Function: ReportHit()
 * Date:     RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Given j,d, coordinates, a score, and a search_results_t data type,
 * adds result into the set of reportable results.  Naively adds hit.
 *
 * Non-overlap algorithm is now done in the scanning routine by Sean's
 * Semi-HMM code.  I've just kept the hit report structure for convenience.
 */
void ReportHit (int i, int j, int bestr, float score, search_results_t *results) 
{
  /*  printf("ReportHit() n:%5d i:%5d j:%5d r:%5d sc:%8.3f\n", results->num_results, i, j, bestr, score);*/

  if(results == NULL) cm_Fail("in ReportHit, but results is NULL\n");
  if (results->num_results == results->num_allocated) 
    ExpandResults (results, INIT_RESULTS);

  results->data[results->num_results].score = score;
  results->data[results->num_results].start = i;
  results->data[results->num_results].stop = j;
  results->data[results->num_results].bestr = bestr;
  results->data[results->num_results].tr = NULL;
  results->data[results->num_results].pcode = NULL;
  results->num_results++;
}


/* Function: PrintResults
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easelfied: EPN, Fri Dec  8 08:29:05 2006 
 * Purpose:  Given the needed information, prints the results.
 *
 *           cm                  the model
 *           fp                  open file ptr to print to
 *           tabfp               open file ptr to print optional tabular output to
 *           si                  SearchInfo, relevant round is final one, si->nrounds
 *           abc                 alphabet to use for output
 *           cons                consensus seq for model (query seq)
 *           dbseq               the database seq
 *           name                sequence name
 *           len                 length of the sequence
 *           do_top              are we doing the plus  (top)    strand?
 *           do_bottom           are we doing the minus (bottom) strand?
 *           do_noncompensatory  are we printing the top line that marks all negative scoring and half bps?
 *           do_noncanonical     are we printing the top line that marks all negative scoring non-canonical and half bps?
 *           namewidth           max length of a name in the target file, for pretty formatting 
 */
void PrintResults (CM_t *cm, FILE *fp, FILE *tabfp, SearchInfo_t *si, const ESL_ALPHABET *abc, CMConsensus_t *cons, dbseq_t *dbseq, 
		   int do_top, int do_bottom, int do_noncompensatory, int do_noncanonical, int namewidth)
{
  int i;
  char *name;
  int len;
  search_results_t *results;
  Fancyali_t *ali;
  int in_revcomp;
  int header_printed = 0;
  int gc_comp;
  float score_for_Eval; /* the score we'll determine the statistical significance of */
  CMEmitMap_t *emap;    /* consensus emit map for the CM */
  int do_stats;        
  ExpInfo_t **exp;      /* pointer to exp tail to use */
  int cm_exp_mode;      /* exp tail mode if we're using CM hits */
  int cp9_exp_mode;     /* exp tail mode if we're using HMM hits */
  int p;                /* relevant partition */
  int offset;         
  int init_rci;         /* initial strand that's been searched, 0 if do_top, else 1 */
  double Eval, Pval;    /* E value and P value of a hit */
  int cm_namewidth;     /* length for printing model name field to tab file */

  /* Contract check: we allow the caller to specify the alphabet they want the 
   * resulting MSA in, but it has to make sense (see next few lines). */
  if(cm->abc->type == eslRNA) { 
      if(abc->type != eslRNA && abc->type != eslDNA)
	cm_Fail("PrintResults(), cm alphabet is RNA, but requested output alphabet is neither DNA nor RNA.");
  }
  else if(cm->abc->K != abc->K) cm_Fail("PrintResults(), cm alphabet size is %d, but requested output alphabet size is %d.", cm->abc->K, abc->K);
  if(si == NULL) cm_Fail("PrintResults(), si == NULL.\n");
  if(si->stype[si->nrounds] != SEARCH_WITH_HMM && si->stype[si->nrounds] != SEARCH_WITH_CM) cm_Fail("PrintResults(), final search round is neither SEARCH_WITH_HMM nor SEARCH_WITH_CM.\n");
  if((!do_top) && (!do_bottom)) cm_Fail("PrintResults(), do_top FALSE, and do_bottom FALSE, what's the point?\n");

  if((si->cutoff_type[si->nrounds] == E_CUTOFF)  && !(cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("PrintResults(), stats wanted but CM has no exp tail stats\n");
  do_stats = (cm->flags & CMH_EXPTAIL_STATS) ? TRUE : FALSE;

  if(do_noncompensatory && do_noncanonical) { 
    cm_Fail("PrintResults(), incompatible flags do_noncompensatory and do_noncanonical both set as TRUE, caller should pick one only.\n");
  }

  if(do_stats) { /* determine exp tail mode to use */
    CM2ExpMode(cm, si->search_opts[si->nrounds], &cm_exp_mode, &cp9_exp_mode);
    exp = (si->stype[si->nrounds] == SEARCH_WITH_HMM) ? cm->stats->expAA[cp9_exp_mode] : cm->stats->expAA[cm_exp_mode];
  }
  emap = CreateEmitMap(cm);
  name = dbseq->sq[0]->name;
  len  = dbseq->sq[0]->n;

  init_rci = do_top ? 0 : 1; 
  for (in_revcomp = init_rci; in_revcomp <= do_bottom; in_revcomp++) {
    results = dbseq->results[in_revcomp];
    if (results == NULL || results->num_results == 0) continue;
      
    if (!header_printed) {
      header_printed = 1;
      fprintf(fp, ">%s\n\n", name);
    }
    fprintf(fp, "  %s strand results:\n\n", in_revcomp ? "Minus" : "Plus");
    
    /* sort hits by bit score */
    SortResultsByScore(results);

    cm_namewidth = ESL_MAX(strlen(cm->name), strlen("model name"));
  
    for (i=0; i<results->num_results; i++) {
      gc_comp = get_gc_comp (cm->abc, dbseq->sq[in_revcomp]->dsq, 
			     results->data[i].start, results->data[i].stop);
      fprintf(fp, " Query = %d - %d, Target = %d - %d\n", 
	      (emap->lpos[cm->ndidx[results->data[i].bestr]] + 1 
	       - StateLeftDelta(cm->sttype[results->data[i].bestr])),
	      (emap->rpos[cm->ndidx[results->data[i].bestr]] - 1 
	       + StateRightDelta(cm->sttype[results->data[i].bestr])),
	      COORDINATE(in_revcomp, results->data[i].start, len), 
	      COORDINATE(in_revcomp, results->data[i].stop, len));

      if(tabfp != NULL) { /* print tabular output also */
	fprintf(tabfp, "  %-*s  %-*s  %10d  %10d  %5d  %5d  %8.2f  ", 
		cm_namewidth, /* max(length of "model name", length(cm->name)) */
		cm->name,     /* model name */
		namewidth,    /* max length of target seq in file, passed in by caller */
		name,         /* target seq name */
		COORDINATE(in_revcomp, results->data[i].start, len), 
		COORDINATE(in_revcomp, results->data[i].stop, len),
		(emap->lpos[cm->ndidx[results->data[i].bestr]] + 1 
		 - StateLeftDelta(cm->sttype[results->data[i].bestr])),
		(emap->rpos[cm->ndidx[results->data[i].bestr]] - 1 
		 + StateRightDelta(cm->sttype[results->data[i].bestr])),
		results->data[i].score);
      }

      if (do_stats) {
	p = cm->stats->gc2p[gc_comp];
	score_for_Eval = results->data[i].score;
	Pval = esl_exp_surv((double) score_for_Eval, exp[p]->mu_extrap, exp[p]->lambda);
	Eval = Pval * exp[p]->cur_eff_dbsize;
	fprintf(fp, " Score = %.2f, E = %.4g, P = %.4g, GC = %3d\n", results->data[i].score, Eval, Pval, gc_comp);
	if(tabfp != NULL) { 
	  fprintf(tabfp, "%4.2e  %3d\n", Eval, gc_comp);
	}
      } 
      else { /* don't print E-value stats */
	fprintf(fp, " Score = %.2f, GC = %3d\n", results->data[i].score, gc_comp);
	if(tabfp != NULL) fprintf(tabfp, "%8s  %3d\n", "-", gc_comp);
      }
      fprintf(fp, "\n");
      if (results->data[i].tr != NULL) {
	/* careful here, all parsetrees have emitl/emitr sequence indices
	 * relative to the hit subsequence of the dsq (i.e. emitl[0] always = 1),
	 * so we pass dsq + start-1.
	 */
	ali = CreateFancyAli (abc, results->data[i].tr, cm, cons, 
			      dbseq->sq[in_revcomp]->dsq + 
			      (results->data[i].start-1), do_noncanonical,
			      results->data[i].pcode);
	
	if(in_revcomp) offset = len - 1;
	else           offset = 0;
	PrintFancyAli(fp, ali,
		      (COORDINATE(in_revcomp, results->data[i].start, len)-1), /* offset in sq index */
		      in_revcomp, (do_noncanonical || do_noncompensatory));
	FreeFancyAli(ali);
	fprintf(fp, "\n");
      }
    }
  }
  fflush(stdout);
  FreeEmitMap(emap);
}

/*
 * Function: ScoresFromResults()
 * Date:     EPN, Thu Jan 24 05:36:54 2008
 * Purpose:  Given a search_results_t results data structure, return the scores within 
 *           it as a vector of floats in <ret_scA> and the number of scores you've
 *           returned in <ret_scN>.
 *
 * Returns:  eslOK on success, <ret_scA> alloc'ed and filled with scores, <ret_scN> set
 *           as number of scores in <ret_scA>. If <ret_scN> == 0, <ret_scA> will be NULL.
 */
int ScoresFromResults(search_results_t *results, char *errbuf, float **ret_scA, int *ret_scN) 
{
  int status;
  float *scA = NULL;
  int    scN = results->num_results;
  int    i;

  if(ret_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "ScoresFromResults(), ret_scA is NULL.");

  ESL_ALLOC(scA, sizeof(float) * scN);
  for(i = 0; i < scN; i++) scA[i] = results->data[i].score;
  *ret_scA = scA;
  *ret_scN = scN;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "ScoresFromResults(), memory allocation error.");
}


/* Function: CountScanDPCalcs()
 * Date:     EPN, Wed Aug 22 09:08:03 2007
 *
 * Purpose:  Count all DP calcs for a CM scan against a 
 *           sequence of length L. Similar to smallcyk.c's
 *           CYKDemands() but takes into account number of
 *           transitions from each state, and is concerned
 *           with a scanning dp matrix, not an alignment matrix.
 *
 * Args:     cm     - the model
 *           L      - length of sequence
 *           use_qdb- TRUE to enforce cm->dmin and cm->dmax for calculation
 *
 * Returns: (float) the total number of DP calculations, either using QDB or not.
 */
float
CountScanDPCalcs(CM_t *cm, int L, int use_qdb)
{
  int v, j;
  float dpcalcs = 0.;
  float dpcalcs_bif = 0.;
  
  /* Contract check */
  if(cm->W > L) cm_Fail("ERROR in CountScanDPCalcs(), cm->W: %d exceeds L: %d\n", cm->W, L);

  float  dpcells     = 0.;
  float  dpcells_bif = 0.;
  int d,w,y,kmin,kmax, bw;

  if(! use_qdb) 
    {
      dpcells = (cm->W * L) - (cm->W * (cm->W-1) * .5); /* fillable dp cells per state (deck) */
      for (j = 1; j < cm->W; j++)
	dpcells_bif += ((j    +2) * (j    +1) * .5) - 1;
      for (j = cm->W; j <= L; j++)
	dpcells_bif += ((cm->W+2) * (cm->W+1) * .5) - 1; 
      dpcalcs_bif = CMCountStatetype(cm, B_st) * dpcells_bif; /* no choice of transition */
      for(v = 0; v < cm->M; v++)
	{
	  if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	    {
	      dpcalcs += dpcells * cm->cnum[v]; /* cnum choices of transitions */

	      /* non-obvious subtractions that match implementation in scancyk.c::CYKScan() */
	      if(v == 0) dpcalcs  -= dpcells; 
	      if(cm->sttype[v] == MP_st) dpcalcs  -= L * cm->cnum[v];
	    }
	}
    }
  else /* use_qdb */
    {
      for(v = 0; v < cm->M; v++)
	{
	  if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	    {
	      bw = cm->dmax[v] - cm->dmin[v] + 1; /* band width */
	      if(cm->dmin[v] == 0) bw--;
	      dpcalcs += ((bw * L) - (bw * (bw-1) * 0.5)) * cm->cnum[v];

	      /* non-obvious subtractions that match implementation in cm_qdband.c::CYKBandedScan() */
	      if(v == 0) dpcalcs  -= ((bw * L) - (bw * (bw-1) * 0.5)); 
	      if(cm->sttype[v] == MP_st) dpcalcs  -= bw * cm->cnum[v];
	    }

	  else if(cm->sttype[v] == B_st)
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      for (j = 1; j <= L; j++)
		{
		  d = (cm->dmin[v] > 0) ? cm->dmin[v] : 1;
		  for (; d <= cm->dmax[v] && d <= j; d++)
		    {
		      if(cm->dmin[y] > (d-cm->dmax[w])) kmin = cm->dmin[y];
		      else kmin = d-cm->dmax[w];
		      if(kmin < 0) kmin = 0;
		      if(cm->dmax[y] < (d-cm->dmin[w])) kmax = cm->dmax[y];
		      else kmax = d-cm->dmin[w];
		      if(kmin <= kmax)
			{
			  bw = (kmax - kmin + 1);
			  dpcalcs_bif += bw;
			}
		    }
		}
	    }
	}
    }
  /*printf("%d CountScanDPCalcs dpc     %.0f\n", use_qdb, dpcalcs);
    printf("%d CountScanDPCalcs dpc_bif %.0f\n", use_qdb, dpcalcs_bif);
    printf("%d CountScanDPCalcs total   %.0f\n", use_qdb, dpcalcs + dpcalcs_bif);*/
  return dpcalcs + dpcalcs_bif;
}


/* Function: CreateBestFilterInfo()
 * Date:     EPN, Mon Dec 10 12:16:22 2007
 *
 * Purpose:  Allocate and initialize a best filter info object. 
 *            
 * Returns:  Newly allocated BestFilterInfo_t object on success, NULL if some error occurs
 */
BestFilterInfo_t *
CreateBestFilterInfo()
{
  int status;

  BestFilterInfo_t *bf = NULL;
  ESL_ALLOC(bf, sizeof(BestFilterInfo_t));

  bf->cm_M      = 0;
  bf->ftype     = FILTER_NOTYETSET;
  bf->cm_eval   = 0.;
  bf->F         = 0.;
  bf->N         = 0;
  bf->db_size   = 0;
  bf->full_cm_ncalcs       = 1.; /* so calc'ed spdup = fil_plus_surv_ncalcs / full_cm_ncalcs = eslINFINITY */
  bf->fil_ncalcs           = eslINFINITY;
  bf->fil_plus_surv_ncalcs = eslINFINITY;
  bf->e_cutoff  = 0.;
  bf->hbeta     = 0.;
  bf->v_isroot  = NULL;   
  bf->np        = 0;
  bf->hexpA     = NULL;   
  bf->is_valid  = FALSE;
  return bf;

 ERROR:
  return NULL; /* reached if memory error */
}  

/* Function: SetBestFilterInfoHMM()
 * Date:     EPN, Mon Dec 10 13:30:34 2007
 *
 * Purpose:  Update a BestFilterInfo_t object to type
 *           FILTER_WITH_HMM_VITERBI or FILTER_WITH_HMM_FORWARD.
 *            
 * Returns:  
 *           eslOK on success, eslEINCOMPAT if contract violated
 */
int 
SetBestFilterInfoHMM(BestFilterInfo_t *bf, char *errbuf, int cm_M, float cm_eval, float F, int N, int db_size, float full_cm_ncalcs, int ftype, float e_cutoff, float fil_ncalcs, float fil_plus_surv_ncalcs)
{
  if(ftype != FILTER_WITH_HMM_VITERBI && ftype != FILTER_WITH_HMM_FORWARD) ESL_FAIL(eslEINCOMPAT, errbuf, "SetBestFilterInfoHMM(), ftype is neither FILTER_WITH_HMM_VITERBI nor FILTER_WITH_HMM_FORWARD.\n");
  if(bf->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "SetBestFilterInfoHMM(), bf->is_valid=TRUE (shouldn't happen, only time to set filter=HMM is when initializing BestFilter.\n");

  bf->cm_M      = cm_M;
  bf->cm_eval   = cm_eval;
  bf->F         = F;
  bf->N         = N;
  bf->db_size   = db_size;
  bf->full_cm_ncalcs       = full_cm_ncalcs;

  bf->ftype      = ftype;
  bf->e_cutoff   = e_cutoff;
  bf->fil_ncalcs = fil_ncalcs;
  bf->fil_plus_surv_ncalcs = fil_plus_surv_ncalcs;
  bf->is_valid   = TRUE;
  return eslOK;
}  

/* Function: SetBestFilterInfoHybrid()
 * Date:     EPN, Mon Dec 10 13:36:13 2007
 *
 * Purpose:  Update a BestFilterInfo_t object to type FILTER_WITH_HYBRID.
 *            
 * Returns:  eslOK on success, eslEINCOMPAT if contract violated, eslEMEM if memory error.
 */
int 
SetBestFilterInfoHybrid(BestFilterInfo_t *bf, char *errbuf, int cm_M, float cm_eval, float F, int N, int db_size, float full_cm_ncalcs, float e_cutoff, float fil_ncalcs, float fil_plus_surv_ncalcs, HybridScanInfo_t *hsi, int np, ExpInfo_t **hexpA)
{
  int status;
  int p;

  if(hsi   == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SetBestFilterInfoHybrid(), hsi is NULL.\n");
  if(hexpA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SetBestFilterInfoHybrid(), hexpA is NULL.\n");

  bf->cm_M      = cm_M;
  bf->cm_eval   = cm_eval;
  bf->F         = F;
  bf->N         = N;
  bf->db_size   = db_size;
  bf->full_cm_ncalcs = full_cm_ncalcs;

  bf->ftype     = FILTER_WITH_HYBRID;
  bf->e_cutoff  = e_cutoff;
  bf->fil_ncalcs = fil_ncalcs;
  bf->fil_plus_surv_ncalcs = fil_plus_surv_ncalcs;
  bf->is_valid  = TRUE;

  if(bf->v_isroot != NULL) free(bf->v_isroot); /* probably unnec, but safe */
  ESL_ALLOC(bf->v_isroot, sizeof(int) * bf->cm_M);
  esl_vec_ICopy(hsi->v_isroot, bf->cm_M, bf->v_isroot);

  /* if bf->hexp exists, free it */
  if(bf->np != 0) for(p = 0; p < bf->np; p++) free(bf->hexpA[p]);
  free(bf->hexpA);

  ESL_ALLOC(bf->hexpA, sizeof(ExpInfo_t *) * np);
  bf->np = np;
  for(p = 0; p < bf->np; p++) {
    bf->hexpA[p] = DuplicateExpInfo(hexpA[p]);
    if(bf->hexpA[p] == NULL) goto ERROR;
  }
  return eslOK;

 ERROR:
  return status;
}  


/* Function: FreeBestFilterInfo()
 * Date:     EPN, Mon Dec 10 13:40:43 2007
 *
 * Purpose:  Free a BestFilterInfo_t object, but not the hsi it points to (if any)
 *            
 * Returns:  void
 */
void
FreeBestFilterInfo(BestFilterInfo_t *bf)
{
  int p;
  if(bf->np != 0) for(p = 0; p < bf->np; p++) free(bf->hexpA[p]);
  free(bf->hexpA);
  if(bf->v_isroot != NULL) free(bf->v_isroot);
  free(bf);
  return;
}  


/* Function: DumpBestFilterInfo()
 * Date:     EPN, Mon Dec 10 12:22:10 2007
 *
 * Purpose:  Print out relevant info in a best filter info object.
 *            
 * Returns:  
 *           eslOK on success, dies immediately on some error
 */
void
DumpBestFilterInfo(BestFilterInfo_t *bf)
{
  int v, p;

  if(! (bf->is_valid)) {
    printf("BestFilterInfo_t not yet valid.\n");
    return;
  }

  printf("BestFilterInfo_t:\n");
  if(bf->ftype == FILTER_WITH_HMM_VITERBI)
    printf("type: FILTER_WITH_HMM_VITERBI\n");
  else if(bf->ftype == FILTER_WITH_HMM_FORWARD)
    printf("type: FILTER_WITH_HMM_FORWARD\n");
  else if(bf->ftype == FILTER_WITH_HYBRID) {
    printf("type: FILTER_WITH_HYBRID\n");
    printf("sub CM roots:\n");
    for(v = 0; v < bf->cm_M; v++) { 
      if(bf->v_isroot[v]) printf("\tv: %d\n", v);
    }
    if(bf->hexpA != NULL) 
      for(p = 0; p < bf->np; p++) { 
	printf("\nHybrid Exp tail, partition: %d\n", p);
	debug_print_expinfo(bf->hexpA[p]);
      }
  }
  printf("CM E value cutoff:     %10.4f\n", bf->cm_eval);
  printf("F:                     %10.4f\n", bf->F);
  printf("N:                     %10d\n",   bf->N);
  printf("DB size (for E-vals):  %10d\n",   bf->db_size);
  printf("Filter E value cutoff: %10.4f\n", bf->e_cutoff);
  printf("Full CM scan DP calcs: %10.4f\n", bf->full_cm_ncalcs);
  printf("Filter scan DP calcs:  %10.4f\n", bf->fil_ncalcs);
  printf("Fil + survivor calcs:  %10.4f\n", bf->fil_plus_surv_ncalcs);
  printf("Speedup:               %10.4f\n", (bf->full_cm_ncalcs / bf->fil_plus_surv_ncalcs));
  return;
}  

/* Function: CreateHMMFilterInfo()
 * Date:     EPN, Tue Jan 15 18:04:56 2008
 *
 * Purpose:  Allocate and initialize a HMM filter info object. 
 *            
 * Returns:  Newly allocated HMMFilterInfo_t object on success, NULL if some error occurs
 */
HMMFilterInfo_t *
CreateHMMFilterInfo()
{
  int status;

  HMMFilterInfo_t *hfi = NULL;
  ESL_ALLOC(hfi, sizeof(HMMFilterInfo_t));

  hfi->is_valid  = FALSE;
  hfi->F         = 0.;
  hfi->N         = 0;
  hfi->dbsize    = 0;
  hfi->ncut      = 0;
  hfi->cm_E_cut  = NULL;
  hfi->fwd_E_cut = NULL;
  hfi->always_better_than_Smax = FALSE;
  return hfi;

 ERROR:
  return NULL; /* reached if memory error */
}  

/* Function: SetHMMFilterInfo()
 * Date:     EPN, Tue Jan 15 18:08:31 2008
 *
 * Purpose:  Fill data for a HMMFilterInfo_t object.
 *            
 * Returns:  eslOK on success, eslEINCOMPAT if contract violated
 */
int 
SetHMMFilterInfoHMM(HMMFilterInfo_t *hfi, char *errbuf, float F, int N, int dbsize, int ncut, float *cm_E_cut, float *fwd_E_cut, int always_better_than_Smax)
{
  int status;
  int i;
  if(hfi->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "SetHMMFilterInfoHMM(), hfi->is_valid=TRUE (shouldn't happen, only time to set filter=HMM is when initializing HMMFilter.\n");
  hfi->ncut      = ncut;
  hfi->F         = F;
  hfi->N         = N;
  hfi->dbsize    = dbsize;
  ESL_ALLOC(hfi->cm_E_cut,  sizeof(float) * ncut);
  ESL_ALLOC(hfi->fwd_E_cut, sizeof(float) * ncut);
  for(i = 0; i < ncut; i++) { 
    hfi->cm_E_cut[i]  = cm_E_cut[i];
    hfi->fwd_E_cut[i] = fwd_E_cut[i];
  }
  hfi->always_better_than_Smax = always_better_than_Smax;

  hfi->is_valid   = TRUE;
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "SetHMMFilterInfoHMM(), memory allocation error.");
}  

/* Function: FreeHMMFilterInfo()
 * Date:     EPN, Tue Jan 15 18:13:15 2008
 *
 * Purpose:  Free a HMMFilterInfo_t object
 *            
 * Returns:  void
 */
void
FreeHMMFilterInfo(HMMFilterInfo_t *hfi)
{
  if(hfi->cm_E_cut  != NULL) free(hfi->cm_E_cut);
  if(hfi->fwd_E_cut != NULL) free(hfi->fwd_E_cut);
  free(hfi);
  return;
}  

/* Function: DumpHMMFilterInfo()
 * Date:     EPN, Mon Dec 10 12:22:10 2007
 *
 * Purpose:  Print out relevant info in a hmm filter info object.
 *           Does some expensive calculations (like QDB calc to
 *           get average length of hits) to determine predicted
 *           speedups, etc, when using the filters.
 *            
 * Returns:  eslOK on success, other Easel status code on some error
 */
int
DumpHMMFilterInfo(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int cmi, int namewidth, char *namedashes)
{
  int i;
  int status;
  float avg_hit_len;
  float cm_ncalcs_per_res;
  int   W; /* window size calculated using cm->beta_W */
  float hmm_ncalcs_per_res;
  float cm_bit_sc;
  float hmm_bit_sc;
  float cm_E;
  float hmm_E;

  /* contract checks */
  if(! (cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "DumpHMMFilterInfo(), cm does not have exp tail stats.");

  /* update the CM's exp tail params for current dbsize, this may have already been done by caller, but we do it again to be safe */
  if((status = UpdateExpsForDBSize(cm, errbuf, dbsize)) != eslOK) return status;

  if(! (hfi->is_valid)) {
    fprintf(fp, "HMMFilterInfo_t not yet valid.\n");
    return eslOK;
  }

  if((status = cm_GetAvgHitLen        (cm,      errbuf, &avg_hit_len))        != eslOK) return status;
  if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &hmm_ncalcs_per_res)) != eslOK) return status;
  if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, cm->beta_qdb, &cm_ncalcs_per_res, &W))  != eslOK) return status;

  fprintf(fp, "#\n");
  fprintf(fp, "# %4s  %-*s  %5s  %6s  %7s  %7s  %7s\n", "idx",  namewidth, "name",     "clen",   "F",      "nseq",    "db (Mb)", "always?");
  fprintf(fp, "# %4s  %-*s  %5s  %6s  %7s  %7s  %7s\n", "----", namewidth, namedashes, "-----",  "------", "-------", "-------", "-------");
  fprintf(fp, "%6d  %-*s  %5d  %6.4f  %7d  %7.1f  %7s\n",
	 cmi, namewidth, cm->name, cm->clen, hfi->F, hfi->N, (double) dbsize / 1000000.,
	 hfi->always_better_than_Smax ? "yes" : "no");
  fprintf(fp, "#\n");
  fprintf(fp, "#\n");
  fprintf(fp, "#%5s  %s", "", "CM E-value cutoff / HMM Forward E-value filter cutoff pairs:\n#\n");
  fprintf(fp, "#%5s  %-4s  %10s  %6s  %10s  %6s  %6s  %7s  %7s\n", "", "idx",  "cm E",       "cm bit", "hmm E",      "hmmbit", "surv",   "xhmm",    "speedup");
  fprintf(fp, "#%5s  %-4s  %10s  %6s  %10s  %6s  %6s  %7s  %7s\n", "", "----", "----------", "------", "----------", "------", "------", "-------", "-------");
  for(i = 0; i < hfi->ncut; i++) {
    cm_E  = hfi->cm_E_cut[i]  * ((double) dbsize / (double) hfi->dbsize);
    hmm_E = hfi->fwd_E_cut[i] * ((double) dbsize / (double) hfi->dbsize);
    if((status = E2MinScore(cm, errbuf, cm_mode,  cm_E,  &cm_bit_sc))  != eslOK) return status;
    if((status = E2MinScore(cm, errbuf, hmm_mode, hmm_E, &hmm_bit_sc)) != eslOK) return status;
    fprintf(fp, "%4s  %6d  ", "", (i+1));
    if(cm_E < 0.01)  fprintf(fp, "%10.2e  ", cm_E);
    else             fprintf(fp, "%10.3f  ", cm_E);
    fprintf(fp, "%6.1f  ", cm_bit_sc);
    if(hmm_E < 0.01) fprintf(fp, "%10.2e  ", hmm_E);
    else             fprintf(fp, "%10.3f  ", hmm_E);
    fprintf(fp, "%6.1f  ", hmm_bit_sc);
    fprintf(fp, "%6.4f  %7.1f  %7.1f\n", 
	   GetHMMFilterS      (hfi, i, W, avg_hit_len),
	   GetHMMFilterXHMM   (hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res),
	   GetHMMFilterSpeedup(hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res));
  }
  return eslOK;
}  


/* Function: DumpHMMFilterInfoForCME()
 * Date:     EPN, Fri Jan 18 10:48:01 2008
 *
 * Purpose:  Given a CM E-value cutoff <cm_E> for a specified
 *           database size <dbsize>, print info on the HMM filter cutoff 
 *           (if any) that would be used for that CM E-value cutoff
 *           in cmsearch.
 *           Does some expensive calculations (like QDB calc to
 *           get average length of hits) to determine predicted
 *           speedups, etc, when using the filters.
 *            
 * Returns:  eslOK on success, other Easel status code on some error
 *           optionally (if non-NULL):
 *           <ret_cm_bit_sc>:          cm bit score cutoff
 *           <ret_hmm_E>:              HMM filter threshold E-value
 *           <ret_hmm_bit_sc>:         HMM filter threshold bit score
 *           <ret_S>:                  predicted filter survival fraction 
 *           <ret_xhmm>:               predicted xhmm factor (predicted speed * hmm speed) 
 *           <ret_spdup>:              predicted speedup from using filter versus only CM search with cm->beta_qdb QDBs
 *           <ret_cm_ncalcs_per_res>:  millions of dp calcs for CM search of 1 residue
 *           <ret_hmm_ncalcs_per_res>: millions of dp calcs for HMM search of 1 residue
 *           <ret_do_filter>:          TRUE if filtering predicted to save time, FALSE if not
 */
int
DumpHMMFilterInfoForCME(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int cmi, float cm_E, int do_header, int namewidth, char *namedashes,
			float *ret_cm_bit_sc, float *ret_hmm_E, float *ret_hmm_bit_sc, float *ret_S, float *ret_xhmm, float *ret_spdup, float *ret_cm_ncalcs_per_res, float *ret_hmm_ncalcs_per_res, int *ret_do_filter)
{
  int i;
  int status;
  float avg_hit_len;
  float cm_ncalcs_per_res;
  int   W; /* window size calculated using cm->beta_qdb */
  float hmm_ncalcs_per_res;
  float cm_bit_sc;
  float hmm_bit_sc;
  float hmm_E;
  float S;     /* predicted survival fraction */
  float xhmm;  /* filtered scan will take <xhmm> times as long as HMM only scan */
  float spdup; /* predicted speedup of filtered scan relative to CM scan */

  /* contract checks */
  if(! (cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "DumpHMMFilterInfoForCME(), cm does not have Exp Tail stats.");

  /* update the CM's exp tail params for current dbsize, this may have already been done by caller, but we do it again to be safe */
  if((status = UpdateExpsForDBSize(cm, errbuf, dbsize)) != eslOK) return status;

  if(! (hfi->is_valid)) {
    fprintf(fp, "HMMFilterInfo_t not yet valid.\n");
    return eslOK;
  }

  if((status = cm_GetAvgHitLen        (cm,      errbuf, &avg_hit_len))        != eslOK) return status;
  if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &hmm_ncalcs_per_res)) != eslOK) return status;
  if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, cm->beta_qdb, &cm_ncalcs_per_res, &W))  != eslOK) return status;

  if(do_header) { 
    fprintf(fp, "# %4s  %-*s  %5s  %8s  %6s  %6s  %6s  %7s  %7s\n", "idx",  namewidth, "name",            "clen",  "cm E",     "cm bit", "hmmbit", "surv",   "xhmm",    "speedup");
    fprintf(fp, "# %4s  %-*s  %5s  %8s  %6s  %6s  %6s  %7s  %7s\n", "----", namewidth, namedashes,       "-----", "--------", "------", "------", "------", "-------", "-------");
  }
  fprintf(fp, "%6d  %-*s  %5d  ", cmi, namewidth, cm->name, cm->clen);
  if((status = GetHMMFilterFwdECutGivenCME(hfi, errbuf, cm_E, dbsize, &i)) != eslOK) return status;
  if((status = E2MinScore(cm, errbuf, cm_mode,  cm_E,  &cm_bit_sc))  != eslOK) return status;
  if(cm_E < 0.01)  fprintf(fp, "%4.2e  ", cm_E);
  else             fprintf(fp, "%8.3f  ", cm_E);
  fprintf(fp, "%6.1f  ", cm_bit_sc);

  if(i != -1) { 
    hmm_E = hfi->fwd_E_cut[i] * ((double) dbsize / (double) hfi->dbsize);
    if((status = E2MinScore(cm, errbuf, hmm_mode, hmm_E, &hmm_bit_sc)) != eslOK) return status;
    fprintf(fp, "%6.1f  ", hmm_bit_sc);
    S     = GetHMMFilterS      (hfi, i, W, avg_hit_len);
    xhmm  = GetHMMFilterXHMM   (hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res);
    spdup = GetHMMFilterSpeedup(hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res);
  }
  else { /* i == -1, special case: filter isn't worth it */
    fprintf(fp, "%6s  ", "-");
    S    = 1.0;
    xhmm = cm_ncalcs_per_res / hmm_ncalcs_per_res;
    spdup = 1.0;
  }
  fprintf(fp, "%6.4f  %7.1f  %7.1f\n", S, xhmm, spdup);

  if(ret_cm_bit_sc != NULL)          *ret_cm_bit_sc          = cm_bit_sc;
  if(ret_hmm_E      != NULL)         *ret_hmm_E      = (i == -1) ? -1. : hmm_E;
  if(ret_hmm_bit_sc != NULL)         *ret_hmm_bit_sc = (i == -1) ? -1. : hmm_bit_sc;
  if(ret_S     != NULL)              *ret_S                  = S;
  if(ret_xhmm  != NULL)              *ret_xhmm               = xhmm;
  if(ret_spdup != NULL)              *ret_spdup              = spdup;
  if(ret_cm_ncalcs_per_res  != NULL) *ret_cm_ncalcs_per_res  = cm_ncalcs_per_res;
  if(ret_hmm_ncalcs_per_res != NULL) *ret_hmm_ncalcs_per_res = hmm_ncalcs_per_res;
  if(ret_do_filter != NULL)          *ret_do_filter          = (i == -1) ? FALSE : TRUE;
  return eslOK;
}  

/* Function: DumpHMMFilterInfoForCMBitScore()
 * Date:     EPN, Fri Jan 18 10:47:54 2008
 *
 * Purpose:  Given a CM bit score cutoff <cm_bit_sc> for a specified
 *           database size <dbsize>, print info on the HMM filter cutoff 
 *           (if any) that would be used for that CM bit score cutoff
 *           in cmsearch.
 *           Simply determines CM E value that corresponds to <cm_bit_sc> 
 *           then calls that DumpHMMFilterInfoForCME().
 *             
 * Returns:  eslOK on success, other Easel status code on some error
 *           optionally (if non-NULL):
 *           <ret_cm_E>:               cm E value cutoff
 *           <ret_hmm_bit_sc>:         HMM filter threshold bit score
 *           <ret_S>:                  predicted filter survival fraction 
 *           <ret_xhmm>:               predicted xhmm factor (predicted speed * hmm speed) 
 *           <ret_spdup>:              predicted speedup from using filter versus only CM search with cm->beta_qdb QDBs
 *           <ret_cm_ncalcs_per_res>:  millions of dp calcs for CM search OF 1 RESIDUE
 *           <ret_hmm_ncalcs_per_res>: millions of dp calcs for HMM search OF 1 RESIDUE
 *           <ret_do_filter>:          TRUE if filtering predicted to save time, FALSE if not
 */
int
DumpHMMFilterInfoForCMBitScore(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int cmi, float cm_bit_sc, int do_header, int namewidth, char *namespaces,
			       float *ret_cm_E, float *ret_hmm_E, float *ret_hmm_bit_sc, float *ret_S, float *ret_xhmm, float *ret_spdup, float *ret_cm_ncalcs_per_res, float *ret_hmm_ncalcs_per_res, int *ret_do_filter)
{
  int status;
  float cm_E;

  if((status = Score2MaxE(cm, errbuf, cm_mode, cm_bit_sc, &cm_E)) != eslOK)  return status;
  if((status = DumpHMMFilterInfoForCME(fp, hfi, errbuf, cm, cm_mode, hmm_mode, dbsize, cmi, cm_E, do_header, namewidth, namespaces, 
				       NULL, ret_hmm_E, ret_hmm_bit_sc, ret_S, ret_xhmm, ret_spdup, ret_cm_ncalcs_per_res, ret_hmm_ncalcs_per_res, ret_do_filter)) != eslOK) return status;

  if(ret_cm_E != NULL) *ret_cm_E = cm_E;
  return eslOK;
}  


/* Function: PlotHMMFilterInfo()
 * Date:     EPN, Mon Dec 10 12:22:10 2007
 *
 * Purpose:  Print out HMM filter stats for each cut point in an 
 *           HMMFilterInfo_t object in xmgrace format. Run in 
 *           8 possible modes, depending on value of <mode>:
 * 
 *           mode  #define (structs.h)  x values            y values
 *           ----  -------------------  ------------------  --------------------
 *              0  FTHR_PLOT_CME_HMME   CM E-value cutoffs  HMM E-value cutoffs 
 *              1  FTHR_PLOT_CME_S      CM E-value cutoffs  predicted survival fraction (S)
 *              2  FTHR_PLOT_CME_XHMM   CM E-value cutoffs  predicted xhmm (factor slower than HMM only scan) 
 *              3  FTHR_PLOT_CME_SPDUP  CM E-value cutoffs  predicted speedup using filter 
 *              4  FTHR_PLOT_CMB_HMMB   CM bit sc cutoffs   HMM bit score cutoffs 
 *              5  FTHR_PLOT_CMB_S      CM bit sc cutoffs   predicted survival fraction (S)
 *              6  FTHR_PLOT_CMB_XHMM   CM bit sc cutoffs   predicted xhmm (factor slower than HMM only scan) 
 *              7  FTHR_PLOT_CMB_SPDUP  CM bit sc cutoffs   predicted speedup using filter 
 *            
 * Returns:  eslOK on success, other Easel status code on some error
 */
int
PlotHMMFilterInfo(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int mode)
{
  int i;
  int status;
  float avg_hit_len;
  float cm_ncalcs_per_res;
  int   W; /* window size calculated using cm->beta_qdb */
  float hmm_ncalcs_per_res;
  float cm_bit_sc;
  float hmm_bit_sc;
  float cm_E;
  float hmm_E;

  /* contract checks */
  if(mode < 0 || mode >= FTHR_NPLOT) ESL_FAIL(eslEINCOMPAT, errbuf, "PlotHMMFilterInfo(), mode: %d is outside allowed range of [%d-%d]", mode, 0, (FTHR_NPLOT-1));

  /* updated the CM's exp tail params for current dbsize, this may have already been done by caller, but we do it again to be safe */
  if((status = UpdateExpsForDBSize(cm, errbuf, dbsize)) != eslOK) return status;

  if(! (hfi->is_valid)) {
    fprintf(fp, "HMMFilterInfo_t not yet valid.\n");
    return eslOK;
  }

  if(mode != FTHR_PLOT_CME_HMME) { /* these calculations are unnec for FTHR_PLOT_CME_HMME */
    if((status = cm_GetAvgHitLen        (cm,      errbuf, &avg_hit_len))        != eslOK) return status;
    if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &hmm_ncalcs_per_res)) != eslOK) return status;
    if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, cm->beta_qdb, &cm_ncalcs_per_res, &W))  != eslOK) return status;
  }
  for(i = 0; i < hfi->ncut; i++) {
    cm_E  = hfi->cm_E_cut[i]  * ((double) dbsize / (double) hfi->dbsize);
    hmm_E = hfi->fwd_E_cut[i] * ((double) dbsize / (double) hfi->dbsize);
    if((status = E2MinScore(cm, errbuf, cm_mode,  cm_E,  &cm_bit_sc))  != eslOK) return status;
    if((status = E2MinScore(cm, errbuf, hmm_mode, hmm_E, &hmm_bit_sc)) != eslOK) return status;

    switch (mode) { 
    case FTHR_PLOT_CME_HMME:  fprintf(fp, "%g\t%g\n", cm_E, hmm_E); break;
    case FTHR_PLOT_CME_S:     fprintf(fp, "%g\t%f\n", cm_E, GetHMMFilterS(hfi, i, W, avg_hit_len)); break;
    case FTHR_PLOT_CME_XHMM:  fprintf(fp, "%g\t%f\n", cm_E, GetHMMFilterXHMM(hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res)); break; 
    case FTHR_PLOT_CME_SPDUP: fprintf(fp, "%g\t%f\n", cm_E, GetHMMFilterSpeedup(hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res)); break;
    case FTHR_PLOT_CMB_HMMB:  fprintf(fp, "%f\t%f\n", cm_bit_sc, hmm_bit_sc); break; 
    case FTHR_PLOT_CMB_S:     fprintf(fp, "%f\t%f\n", cm_bit_sc, GetHMMFilterS(hfi, i, W, avg_hit_len)); break;
    case FTHR_PLOT_CMB_XHMM:  fprintf(fp, "%f\t%f\n", cm_bit_sc, GetHMMFilterXHMM(hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res)); break;
    case FTHR_PLOT_CMB_SPDUP: fprintf(fp, "%f\t%f\n", cm_bit_sc, GetHMMFilterSpeedup(hfi, i, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res)); break;
    }
  }
  fprintf(fp, "&\n");
  return eslOK;
}    

/* Function: GetHMMFilterS()
 * Date:     EPN, Wed Jan 16 21:21:55 2008
 *
 * Purpose:  Return the survival fraction S for a given
 *           cut point in an HMM filter object.
 *            
 * Returns:  Survival fraction for cut point 'cut',
 *           dies if cut >= hfi->ncut (out of bounds)
 */
float 
GetHMMFilterS(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len)
{
  if(cut >= hfi->ncut) cm_Fail("HMMFilterS() request cut point %d, when only %d exist.", cut, hfi->ncut);
  float surv_res_per_hit = ((float) 2 * W) - avg_hit_len;
  return((hfi->fwd_E_cut[cut] * surv_res_per_hit) / (float) hfi->dbsize);
}

/* Function: GetHMMFilterTotalCalcs()
 * Date:     EPN, Wed Jan 16 21:37:54 2008
 *
 * Purpose:  Returns the predicted number of millions
 *           of DP calcs for a HMM filter scan plus
 *           the CM search of the survivors for
 *           a sequence of length <hfi->dbsize>.

 * Returns:  number of millions of dp calcs
 *           dies if cut >= hfi->ncut (out of bounds)
 */
float 
GetHMMFilterTotalCalcs(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len, float cm_ncalcs_per_res, float hmm_ncalcs_per_res)
{
  if(cut >= hfi->ncut) cm_Fail("HMMFilterS() request cut point %d, when only %d exist.", cut, hfi->ncut);
  float S = GetHMMFilterS(hfi, cut, W, avg_hit_len);
  float f = hmm_ncalcs_per_res * (float) hfi->dbsize; 
  float c = cm_ncalcs_per_res  * (float) hfi->dbsize; 
  return(f + (S * c));
}

/* Function: GetHMMFilterXHMM()
 * Date:     EPN, Wed Jan 16 21:29:37 2008
 *
 * Purpose:  Return the <xhmm> factor for a given
 *           cut point in an HMM filter object.
 *           If <xhmm> = 2.0, the HMM filter plus
 *           the CM search of the survivors should
 *           take twice as long as the HMM ONLY scan.
 *            
 * Returns:  <xhmm>
 *           dies if cut >= hfi->ncut (out of bounds)
 */
float 
GetHMMFilterXHMM(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len, float cm_ncalcs_per_res, float hmm_ncalcs_per_res)
{
  if(cut >= hfi->ncut) cm_Fail("HMMFilterS() request cut point %d, when only %d exist.", cut, hfi->ncut);
  float total_calcs = GetHMMFilterTotalCalcs(hfi, cut, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res);
  float f = hmm_ncalcs_per_res * (float) hfi->dbsize; 
  return(total_calcs / f);
}

/* Function: GetHMMFilterSpeedup()
 * EPN, Wed Jan 16 21:44:13 2008
 *
 * Purpose:  Return the predicted speedup for 
 *           a HMM filter scan plus CM search 
 *           of survivors versus a non-filtered
 *           CM search of hfi->dbsize residues.
 *            
 * Returns:  predicted speedup
 *           dies if cut >= hfi->ncut (out of bounds)
 */
float 
GetHMMFilterSpeedup(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len, float cm_ncalcs_per_res, float hmm_ncalcs_per_res)
{
  if(cut >= hfi->ncut) cm_Fail("HMMFilterS() request cut point %d, when only %d exist.", cut, hfi->ncut);
  float total_calcs = GetHMMFilterTotalCalcs(hfi, cut, W, avg_hit_len, cm_ncalcs_per_res, hmm_ncalcs_per_res);
  float c = cm_ncalcs_per_res  * (float) hfi->dbsize; 
  return(c / total_calcs);
}

/* Function: GetHMMFilterFwdECutGivenCME()
 * EPN, Fri Jan 18 09:48:55 2008
 *
 * Purpose:  Given a CM E value cutoff <cm_E> to use in 
 *           cmsearch for a search of a database of <dbsize> 
 *           size, determine the appropriate HMM filter cutoff 
 *           to use for a HMM forward filter, if any.
 *            
 * Returns:  <ret_cut_pt>: index of forward filter E value cutoff to use in hfi->fwd_E_cut[].
 *           if we shouldn't bother filtering, we set *ret_cut_pt to -1.
 *           eslOK on success, eslEINCOMPAT on contract violation.
 */
int
GetHMMFilterFwdECutGivenCME(HMMFilterInfo_t *hfi, char *errbuf, float cm_E, long dbsize, int *ret_cut_pt)
{
  int c = -1;
  double dbsize_factor = (double) dbsize / (double) hfi->dbsize;

  /* contract check */
  if(ret_cut_pt == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "GetHMMFilterFwdECutGivenCME(), ret_cut_pt == NULL");

  while(((c+1) < hfi->ncut)  /* we're not at final cut point */
	&& ((hfi->cm_E_cut[(c+1)] * dbsize_factor) > cm_E))  /* next cut point is still > E_cut */
    { c++; }
  /* now c is the max c for which cm_E_cut[c] * dbsize_factor < cm_E 
   * boundary cases: 
   *
   * if c == hfi->ncut-1: cm_E_cut[c] * dbsize_factor > cm_E for all c == 0..hfi->ncut-1
   *                      This means cm_E is better than the ((1. - FTHR_MAXFRACT) * hfi->N)th
   *                      best CM E value we observed after sampling and finding the optimal CM hit in 
   *                      hfi->N seqs in cmcalibrate. In this case, we use the fwd_E_cut[hfi->ncut-1] 
   *                      for the forward filter. This should be safe b/c this fwd_E value has empirically
   *                      shown it can find hfi->F fraction for a worse CM E-value cutoff than cm_E.
   *
   * if c == -1:          cm_E_cut[0] * dbsize_factor < cm_E
   *                      What this means depends on the value of the hfi->always_better_than_Smax parameter:
   *
   *                      (A) if hfi->always_better_than_Smax == FALSE: 
   *                          cm_E is greater than the maximum cm E value cutoff we determined it was 
   *                          worth filtering for in cmcalibrate. In this case we turn filtering off 
   *                          by returning -1.
   *
   *                      (B) if hfi->always_better_than_Smax == TRUE: 
   *                          cm_E is greater than the maximum cm E value we observed after sampling 
   *                          and finding the optimal CM hit in hfi->N seqs in cmcalibrate. In this 
   *                          case, it's unclear what the best stratgey is. The current strategy is 
   *                          to use the fwd_E_cut[0] E cutoff by returning 0. This is potentially 
   *                          unsafe because it means there are sequences that have E values E such 
   *                          that 
   *                             cm_E > E > cm_E_cut[0]
   *                          that our HMM filter will likely miss, but there's really nothing we 
   *                          can do about it because we have no way of sampling these sequences
   *                          w.r.t the true CM distribution over sequences. That's exactly what 
   *                          we tried to do in cmcalibrate but we did not sample deeply enough.
   *                          (To get around this case, I tried some tricks with exponentiating 
   *                          the CM's parameters by <x> < 1 to try to sample sequences with lower
   *                          scores but this mucks with the CM distribution in a way I don't understand, 
   *                          so it was abandonded.)
   */
  if(c == -1 && hfi->always_better_than_Smax) c = 0; /* case B above */
  *ret_cut_pt = c;
  return eslOK;
}

/* Function: GetHMMFilterFwdECutGivenCMBitScore()
 * EPN, Fri Jan 18 10:37:20 2008
 *
 * Purpose:  Given a CM bit score cutoff <cm_bit_sc> to use in 
 *           cmsearch for a search of a database of <dbsize> 
 *           size, determine the appropriate HMM filter cutoff 
 *           to use for a HMM forward filter, if any.
 *            
 * Returns:  <ret_cut_pt>: index of forward filter E value cutoff to use in hfi->fwd_E_cut[].
 *           if we shouldn't bother filtering, we set *ret_cut_pt to -1.
 *           eslOK on success, eslEINCOMPAT on contract violation.
 */
int
GetHMMFilterFwdECutGivenCMBitScore(HMMFilterInfo_t *hfi, char *errbuf, float cm_bit_sc, long dbsize, int *ret_cut_pt, CM_t *cm, int cm_mode)
{
  int status;
  float cm_E; /* the max 'safe' E value cm_bit_sc corresponds to, across all partitions */
  int cut_pt;

  /* contract check */
  if(ret_cut_pt == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "GetHMMFilterFwdECutGivenCMBitScore(), ret_cut_pt == NULL");

  if((status = Score2MaxE(cm, errbuf, cm_mode, cm_bit_sc, &cm_E)) != eslOK)          return status;
  if((status = GetHMMFilterFwdECutGivenCME(hfi, errbuf, cm_E, dbsize, &cut_pt)) != eslOK) return status; 

  *ret_cut_pt = cut_pt;

  return eslOK;
}


/* Function: SurvFract2E()
 * Date:     EPN, Mon Jan 21 09:27:07 2008
 *
 * Purpose:  Return the E-value that will give 
 *           a specific survival fraction S given
 *           W, avg_hit_len, and dbsize.
 */
float 
SurvFract2E(float S, int W, float avg_hit_len, long dbsize)
{
  float surv_res_per_hit = ((float) 2 * W) - avg_hit_len;
  return((S * (double) dbsize) / surv_res_per_hit);
}


/* Function: E2SurvFract()
 * Date:     EPN, Mon Jan 21 09:27:07 2008
 *
 * Purpose:  Return the survival fraction S given 
 *           an E-value cutoff, W, avg_hit_len, and dbsize.
 */
float 
E2SurvFract(float E, int W, float avg_hit_len, long dbsize, int do_pad)
{
  float surv_res_per_hit;
  float surv_fract;
  surv_res_per_hit = (do_pad) ? (((float) 2 * W) - avg_hit_len) : avg_hit_len;
  surv_fract = (E * surv_res_per_hit) / ((double) dbsize);
  surv_fract = ESL_MIN(surv_fract, 1.0); /* max survival fraction is 1.0 */
  return surv_fract;
}


/* Function: Results2SurvFract()
 * Date:     EPN, Mon Mar  3 14:07:59 2008
 *
 * Purpose:  Return the survival fraction from a round
 *           of filtering, after possibly adding a W pad to both sides
 *           if (do_pad) (i.e. if i..j is a hit, j-(W-1)..i+(W-1) survives)
 *           and possibly collapsing overlapping hits together (if do_collapse).
 */
int
Results2SurvFract(CM_t *cm, char *errbuf, int i0, int j0, search_results_t *results, int do_pad, int do_collapse, float *ret_survfract)
{
  int prev_j, nhits, h, next_j, i, j;
  int surv_res = 0;
  float survfract;

  prev_j = j0;
  nhits  = results->num_results;
  /* To be safe, we only trust that i..j of our filter-passing hit is within the real hit,
   * so we add (W-1) to start point i and subtract (W-1) from j, and treat this region j-(W-1)..i+(W-1)
   * as having survived the filter.
   */
  for(h = 0; h < nhits; h++) {
    if(results->data[h].stop > prev_j) ESL_FAIL(eslEINCOMPAT, errbuf, "Results2SurvFract(): j's not in descending order");
    prev_j = results->data[h].stop;
    
    if(do_pad) { 
      i = ((results->data[h].stop  - (cm->W-1)) >= 1)    ? (results->data[h].stop  - (cm->W-1)) : 1;
      j = ((results->data[h].start + (cm->W-1)) <= j0)   ? (results->data[h].start + (cm->W-1)) : j0;
      
      if((h+1) < nhits) next_j = ((results->data[h+1].start + (cm->W-1)) <= j0) ? (results->data[h+1].start + (cm->W-1)) : j0;
      else              next_j = -1;
      
      /* possibly collapse multiple overlapping hits together into a single hit. 
       */
      if(do_collapse) { 
	while(((h+1) < nhits) && (next_j >= i)) { /* suck in hit */
	  h++;
	  i = ((results->data[h].stop - (cm->W-1)) >= 1) ? (results->data[h].stop - (cm->W-1)) : 1;
	  if((h+1) < nhits) next_j = ((results->data[h+1].start + (cm->W-1)) <= j0) ? (results->data[h+1].start + (cm->W-1)) : j0;
	  else              next_j = -1;
	}
      }
    } /* end of if(do_pad) */
    else { /* do_pad == FALSE */
      i = results->data[h].start;
      j = results->data[h].stop;
      
      if((h+1) < nhits) next_j = results->data[h+1].stop;
      else              next_j = -1;
      
      /* possibly collapse multiple overlapping hits together into a single hit. 
       */
      if(do_collapse) { 
	while(((h+1) < nhits) && (next_j >= i)) { /* suck in hit */
	  h++;
	  i = results->data[h].start;
	  if((h+1) < nhits) next_j = results->data[h+1].stop;
	  else              next_j = -1;
	}
      }
    } 
    surv_res += j-i+1;
  }
  survfract = ((float) surv_res) / ((float) (j0-i0+1));
  *ret_survfract = survfract;
  return eslOK;
}

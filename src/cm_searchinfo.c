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
    results->data[i].pcode1 = NULL;
    results->data[i].pcode2 = NULL;
  }
  return (results);
 ERROR:
  esl_fatal("Memory allocation error.");
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
    results->data[i].pcode1 = NULL;
    results->data[i].pcode2 = NULL;
  }

  results->num_allocated+=additional;
  return;
 ERROR:
  esl_fatal("Memory reallocation error.");
}

/* Function: AppendResults()
 * Date:     EPN, Wed Aug 29 08:58:28 2007
 *
 * Purpose:  Add result nodes from one results structure onto
 *           another by copying data and manipulating pointers. 
 *           Originally written to add results returned from 
 *           an MPI worker to a growing 'master' results structure 
 *           in the MPI master. 

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
      report_hit (src_results->data[i].start+i0-1, src_results->data[i].stop+i0-1, 
		  src_results->data[i].bestr,      src_results->data[i].score,
		  dest_results);
      if(src_results->data[i].tr != NULL)
	(*dest_results).data[ip].tr = (*src_results).data[i].tr;
      if(src_results->data[i].pcode1 != NULL)
	(*dest_results).data[ip].pcode1 = (*src_results).data[i].pcode1;
      if(src_results->data[i].pcode2 != NULL)
	(*dest_results).data[ip].pcode2 = (*src_results).data[i].pcode2;
    }
  return;
}

/*
 * Function: FreeResults()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Frees a results structure
 */
void FreeResults (search_results_t *r) {
  int i;
  if (r != NULL) {
    for (i=0; i < r->num_allocated; i++) {
      if (r->data[i].tr     != NULL) FreeParsetree(r->data[i].tr);
      if (r->data[i].pcode1 != NULL) free(r->data[i].pcode1);
      if (r->data[i].pcode2 != NULL) free(r->data[i].pcode2);
    }
    free (r->data);
    free(r);
  }
}


/*
 * Function: compare_results()
 * Date:     RJK, Wed Apr 10, 2002 [St. Louis]
 * Purpose:  Compares two search_result_node_ts based on score and returns -1
 *           if first is higher score than second, 0 if equal, 1 if first
 *           score is lower.  This results in sorting by score, highest
 *           first.
 */
int compare_results (const void *a_void, const void *b_void) {
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

/*
 * Function: remove_overlapping_hits ()
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
void remove_overlapping_hits (search_results_t *results, int i0, int j0)
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

  sort_results (results);

  for (x=0; x<results->num_results; x++) {
    covered = 0;
    for (y=results->data[x].start; y<=results->data[x].stop && !covered; y++) {
      {
	yp = y-i0+1; 
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

  sort_results(results);
  return;

 ERROR:
  esl_fatal("Memory allocation error.");
}

/*
 * Function: sort_results()
 * Date:    RJK,  Sun Mar 31, 2002 [AA Flight 2869 LGA->STL]
 * Purpose: Given a results array, sorts it with a call to qsort
 *
 */
void sort_results (search_results_t *results) 
{
  qsort (results->data, results->num_results, sizeof(search_result_node_t), compare_results);
}

/*
 * Function: report_hit()
 * Date:     RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Given j,d, coordinates, a score, and a search_results_t data type,
 * adds result into the set of reportable results.  Naively adds hit.
 *
 * Non-overlap algorithm is now done in the scanning routine by Sean's
 * Semi-HMM code.  I've just kept the hit report structure for convenience.
 */
void report_hit (int i, int j, int bestr, float score, search_results_t *results) 
{

  if(results == NULL) 
    esl_fatal("in report_hit, but results is NULL\n");
  if (results->num_results == results->num_allocated) 
    ExpandResults (results, INIT_RESULTS);

  results->data[results->num_results].score = score;
  results->data[results->num_results].start = i;
  results->data[results->num_results].stop = j;
  results->data[results->num_results].bestr = bestr;
  results->data[results->num_results].tr = NULL;
  results->data[results->num_results].pcode1 = NULL;
  results->data[results->num_results].pcode2 = NULL;
  results->num_results++;
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
  if(cm->W > L) esl_fatal("ERROR in CountScanDPCalcs(), cm->W: %d exceeds L: %d\n", cm->W, L);

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

	      /* non-obvious subtractions that match implementation in bandcyk.c::CYKBandedScan() */
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

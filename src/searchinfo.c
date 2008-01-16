/* searchinfo.c
 * SearchInfo_t implementations: information for CM/CP9 
 * filters and scans.
 *
 * EPN, Tue Nov 27 08:42:08 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_gumbel.h"
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
		      ScanMatrix_t *smx, HybridScanInfo_t *hsi, int cutoff_type, float sc_cutoff, float e_cutoff)
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
      if(!do_noalign) cm_Fail("ValidateSearchInfo(), round %d has CM_SEARCH_NOALIGN flag down.\n", n);
      if(si->stype[n] == SEARCH_WITH_HMM) {
	sum = do_noqdb + do_hbanded + do_hmmscanbands + do_sums + do_inside + do_toponly + do_null2 + do_rsearch + do_cmgreedy;
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
	sum = do_hbanded + do_hmmscanbands + do_sums + do_toponly + do_null2 + do_rsearch + do_cmgreedy;	
	if(sum != 0 || (do_hmmviterbi + do_hmmforward != 1)) { 
	  printf("ValidateSearchInfo(), round %d is SEARCH_WITH_HYBRID but search opts are invalid\n", n);
	  DumpSearchOpts(si->search_opts[n]);
	  cm_Fail("This is fatal.");
	}
      }
      else if (si->stype[n] == SEARCH_WITH_CM) {
	if(si->smx[n] == NULL) cm_Fail("ValidateSearchInfo(), round %d is SEARCH_WITH_CM but smx[%d] is NULL\n", n, n);
	sum = do_hbanded + do_hmmscanbands + do_sums + do_toponly + do_null2 + do_rsearch + do_hmmviterbi + do_hmmforward;	
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
  if(cm->si == NULL)           cm_Fail("UpdateSearchInfoCutoff(), cm->si is NULL.");
  if(nround > cm->si->nrounds) cm_Fail("UpdateSearchInfoCutoff(), requested round %d is > cm->si->nrounds\n", nround, cm->si->nrounds);
  cm->si->cutoff_type[nround] = cutoff_type;
  cm->si->sc_cutoff[nround]   = sc_cutoff;
  cm->si->e_cutoff[nround]    = e_cutoff;
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
    results->data[i].pcode1 = NULL;
    results->data[i].pcode2 = NULL;
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
  cm_Fail("Memory allocation error.");
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
    cm_Fail("in report_hit, but results is NULL\n");
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


/* Function: print_results
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easelfied: EPN, Fri Dec  8 08:29:05 2006 
 * Purpose:  Given the needed information, prints the results.
 *
 *           cm                  the model
 *           si                  SearchInfo, relevant round is final one, si->nrounds
 *           abc                 alphabet to use for output
 *           cons                consensus seq for model (query seq)
 *           dbseq               the database seq
 *           name                sequence name
 *           len                 length of the sequence
 *           do_top              are we doing the plus  (top)    strand?
 *           do_bottom           are we doing the minus (bottom) strand?
 *           do_noncompensatory  are we printing the optional non-compensatory line?
 */
void print_results (CM_t *cm, SearchInfo_t *si, const ESL_ALPHABET *abc, CMConsensus_t *cons, dbseq_t *dbseq, int do_top, int do_bottom, int do_noncompensatory)
{
  int i;
  char *name;
  int len;
  search_results_t *results;
  Fancyali_t *ali;
  int in_revcomp;
  int header_printed = 0;
  int gc_comp;
  float score_for_Eval; /* the score we'll determine the statistical significance
			 * of. This will be the bit score stored in
			 * dbseq unless (cm->flags & CM_ENFORCED)
			 * in which case we subtract cm->enf_scdiff first. */
  CMEmitMap_t *emap;    /* consensus emit map for the CM */
  int do_stats;        
  GumbelInfo_t **gum;      /* pointer to gum to use */
  int cm_gum_mode;      /* Gumbel mode if we're using CM hits */
  int cp9_gum_mode;     /* Gumbel mode if we're using HMM hits */
  int p;                /* relevant partition */
  int offset;         
  int init_rci;         /* initial strand that's been searched, 0 if do_top, else 1 */

  /* Contract check: we allow the caller to specify the alphabet they want the 
   * resulting MSA in, but it has to make sense (see next few lines). */
  if(cm->abc->type == eslRNA) { 
      if(abc->type != eslRNA && abc->type != eslDNA)
	cm_Fail("print_results(), cm alphabet is RNA, but requested output alphabet is neither DNA nor RNA.");
  }
  else if(cm->abc->K != abc->K) cm_Fail("print_results(), cm alphabet size is %d, but requested output alphabet size is %d.", cm->abc->K, abc->K);
  if(si == NULL) cm_Fail("print_results(), si == NULL.\n");
  if(si->stype[si->nrounds] != SEARCH_WITH_HMM && si->stype[si->nrounds] != SEARCH_WITH_CM) cm_Fail("print_results(), final search round is neither SEARCH_WITH_HMM nor SEARCH_WITH_CM.\n");
  if((!do_top) && (!do_bottom)) cm_Fail("print_results(), do_top FALSE, and do_bottom FALSE, what's the point?\n");

  do_stats = (si->cutoff_type[si->nrounds] == E_CUTOFF) ? TRUE : FALSE;
  if(do_stats  && !(cm->flags & CMH_GUMBEL_STATS)) cm_Fail("print_results(), stats wanted but CM has no Gumbel stats\n");

  if(do_stats) { /* determine Gumbel mode to use */
    CM2Gumbel_mode(cm, cm->search_opts, &cm_gum_mode, &cp9_gum_mode);
    gum = (si->stype[si->nrounds] == SEARCH_WITH_HMM) ? cm->stats->gumAA[cp9_gum_mode] : cm->stats->gumAA[cm_gum_mode];
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
      printf (">%s\n\n", name);
    }
    printf ("  %s strand results:\n\n", in_revcomp ? "Minus" : "Plus");
    
    /*for (i=0; i<results->num_results; i++) 
      printf("hit: %5d start: %5d stop: %5d len: %5d emitl[0]: %5d emitr[0]: %5d score: %9.3f\n", i, results->data[i].start, results->data[i].stop, len, results->data[i].tr->emitl[0], results->data[i].tr->emitr[0], results->data[i].score);*/
    for (i=0; i<results->num_results; i++) {
      gc_comp = get_gc_comp (dbseq->sq[in_revcomp], 
			     results->data[i].start, results->data[i].stop);
      printf (" Query = %d - %d, Target = %d - %d\n", 
	      (emap->lpos[cm->ndidx[results->data[i].bestr]] + 1 
	       - StateLeftDelta(cm->sttype[results->data[i].bestr])),
	      (emap->rpos[cm->ndidx[results->data[i].bestr]] - 1 
	       + StateRightDelta(cm->sttype[results->data[i].bestr])),
	      COORDINATE(in_revcomp, results->data[i].start, len), 
	      COORDINATE(in_revcomp, results->data[i].stop, len));
      if (do_stats) {
	p = cm->stats->gc2p[gc_comp];
	score_for_Eval = results->data[i].score;
	if(cm->flags & CM_ENFORCED) {
	  printf("\n\torig sc: %.3f", score_for_Eval);
	  score_for_Eval -= cm->enf_scdiff;
	  printf(" new sc: %.3f (diff: %.3f\n\n", score_for_Eval, cm->enf_scdiff);
	}
	printf (" Score = %.2f, E = %.4g, P = %.4g, GC = %3d\n", results->data[i].score,
		RJK_ExtremeValueE(score_for_Eval, gum[p]->mu, 
				  gum[p]->lambda),
		esl_gumbel_surv((double) score_for_Eval, gum[p]->mu, 
				gum[p]->lambda), gc_comp);
	/*printf("  Mu[gc=%d]: %f, Lambda[gc=%d]: %f\n", gc_comp, mu[gc_comp], gc_comp,
	  lambda[gc_comp]);
	  ExtremeValueP(results->data[i].score, mu[gc_comp], 
	  lambda[gc_comp]));*/
      } 
      else { /* don't print E-value stats */
	printf (" Score = %.2f, GC = %3d\n", results->data[i].score, gc_comp);
      }
      printf ("\n");
      if (results->data[i].tr != NULL) {
	/* careful here, all parsetrees have emitl/emitr sequence indices
	 * relative to the hit subsequence of the dsq (i.e. emitl[0] always = 1),
	 * so we pass dsq + start-1.
	 */
	ali = CreateFancyAli (abc, results->data[i].tr, cm, cons, 
			      dbseq->sq[in_revcomp]->dsq + 
			      (results->data[i].start-1), 
			      results->data[i].pcode1, results->data[i].pcode2);
	
	if(in_revcomp) offset = len - 1;
	else           offset = 0;
	PrintFancyAli(stdout, ali,
		      (COORDINATE(in_revcomp, results->data[i].start, len)-1), /* offset in sq index */
		      in_revcomp, do_noncompensatory);
	FreeFancyAli(ali);
	printf ("\n");
      }
    }
  }
  fflush(stdout);
  FreeEmitMap(emap);
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
  bf->hgumA     = NULL;   
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
  if(bf->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "SetBestFilterInfoHMM(), bf->is_valid is TRUE (shouldn't happen, only time to set filter as HMM is when initializing BestFilter object.\n");

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
SetBestFilterInfoHybrid(BestFilterInfo_t *bf, char *errbuf, int cm_M, float cm_eval, float F, int N, int db_size, float full_cm_ncalcs, float e_cutoff, float fil_ncalcs, float fil_plus_surv_ncalcs, HybridScanInfo_t *hsi, int np, GumbelInfo_t **hgumA)
{
  int status;
  int p;

  if(hsi   == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SetBestFilterInfoHybrid(), hsi is NULL.\n");
  if(hgumA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "SetBestFilterInfoHybrid(), hgumA is NULL.\n");

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

  /* if bf->hgum exists, free it */
  if(bf->np != 0) for(p = 0; p < bf->np; p++) free(bf->hgumA[p]);
  free(bf->hgumA);

  ESL_ALLOC(bf->hgumA, sizeof(GumbelInfo_t *) * np);
  bf->np = np;
  for(p = 0; p < bf->np; p++) {
    bf->hgumA[p] = DuplicateGumbelInfo(hgumA[p]);
    if(bf->hgumA[p] == NULL) goto ERROR;
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
  if(bf->np != 0) for(p = 0; p < bf->np; p++) free(bf->hgumA[p]);
  free(bf->hgumA);
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
    if(bf->hgumA != NULL) 
      for(p = 0; p < bf->np; p++) { 
	printf("\nHybrid Gumbel, partition: %d\n", p);
	debug_print_gumbelinfo(bf->hgumA[p]);
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
  if(hfi->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "SetHMMFilterInfoHMM(), hfi->is_valid is TRUE (shouldn't happen, only time to set filter as HMM is when initializing HMMFilter object.\n");
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
 *           
 *            
 * Returns:  
 *           eslOK on success, dies immediately on some error
 */
void
DumpHMMFilterInfo(HMMFilterInfo_t *hfi)
{
  int i;

  if(! (hfi->is_valid)) {
    printf("HMMFilterInfo_t not yet valid.\n");
    return;
  }

  printf("HMMFilterInfo_t:\n");

  printf("F:                     %10.4f\n", hfi->F);
  printf("N:                     %10d\n",   hfi->N);
  printf("DB size (for E-vals):  %10d\n",   hfi->dbsize);
  printf("ncut:                  %10d\n",   hfi->ncut);

  for(i = 0; i < hfi->ncut; i++) 
    printf("%5d cm_E: %15.7f fwd_E: %15.7f\n", i, hfi->cm_E_cut[i], hfi->fwd_E_cut[i]);

  return;
}  


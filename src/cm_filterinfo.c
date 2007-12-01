/* FilterInfo_t implementations: information for filters of CM scans.
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


/* Function: cm_CreateFilterInfo
 * Date:     EPN, Tue Nov 27 12:57:24 2007
 *
 * Purpose:  Allocate and initialize a filter info object that
 *           specifies that no filtering should be done. 
 *            
 * Returns:  cm->fi points to a new FilterInfo_t object
 *           eslOK on success, dies immediately on some error
 */
int
cm_CreateFilterInfo(CM_t *cm, float cutoff)
{
  int status;

  if(cm->fi != NULL)              cm_Fail("cm_CreateFilterInfo(), the cm already points to a FilterInfo_t object.\n");
  
  
  FilterInfo_t *fi;
  ESL_ALLOC(fi, sizeof(FilterInfo_t));
  
  fi->nrounds = 0;
  ESL_ALLOC(fi->search_opts, sizeof(int)   * (fi->nrounds+1));
  ESL_ALLOC(fi->cutoff,      sizeof(float) * (fi->nrounds+1));
  ESL_ALLOC(fi->ftype,       sizeof(int)   * (fi->nrounds+1));
  ESL_ALLOC(fi->hsi,         sizeof(HybridScanInfo_t *) * (fi->nrounds+1));

  fi->search_opts[0] = cm->search_opts;
  fi->cutoff[0] = cutoff;
  fi->ftype[0]  = NO_FILTER;
  fi->hsi[0]    = NULL;

  cm->fi = fi;
  return eslOK;

 ERROR:
  cm_Fail("cm_CreateFilterInfo(), memory allocate error.");
  return status; /* NEVERREACHED */
}  


/* Function: cm_AddHMMFilterInfo()
 * Date:     EPN, Tue Nov 27 13:00:23 2007
 *
 * Purpose:  Add an HMM filter as the 1st round filter for CM <cm>.
 *           A new FilterInfo_t object <fi> is created, and the existing
 *           information from cm->fi is copied into it. cm->fi is then
 *           freed and cm->fi is set to point at fi.            
 * 
 * Returns:  cm->fi points to a new FilterInfo_t object
 *           eslOK on success, dies immediately on some error
 */
int
cm_AddHMMFilterInfo(CM_t *cm, int do_viterbi, float cutoff)
{
  int status;
  int n;
  int orig_nrounds;

  if(cm->fi == NULL)                 cm_Fail("cm_AddHMMFilterInfo(), the cm does not point to a FilterInfo_t object.\n");

  orig_nrounds = cm->fi->nrounds;

  FilterInfo_t *fi;
  ESL_ALLOC(fi, sizeof(FilterInfo_t));

  /* allocate new fi object, with 1 more round than cm->fi, and set round 0 */
  fi->nrounds = orig_nrounds+1;
  ESL_ALLOC(fi->search_opts, sizeof(int)   * (fi->nrounds+1));
  ESL_ALLOC(fi->cutoff,      sizeof(float) * (fi->nrounds+1));
  ESL_ALLOC(fi->ftype,       sizeof(int)   * (fi->nrounds+1));
  ESL_ALLOC(fi->search_opts, sizeof(HybridScanInfo_t *) * (fi->nrounds+1));

  fi->search_opts[0] = 0;
  fi->search_opts[0] |= CM_SEARCH_NOALIGN;
  if(do_viterbi) fi->search_opts[0] |= CM_SEARCH_HMMVITERBI;
  else           fi->search_opts[0] |= CM_SEARCH_HMMFORWARD;
  
  fi->cutoff[0] = cutoff;
  fi->ftype[0]  = FILTER_WITH_HMM;
  fi->hsi[0]    = NULL;

  /* copy existing information for other rounds from old cm->fi */
  for(n = 0; n <= orig_nrounds; n++) { 
    fi->search_opts[(n+1)] = cm->fi->search_opts[n];
    fi->cutoff[(n+1)] = cm->fi->cutoff[n];
    fi->ftype[(n+1)] = cm->fi->ftype[n];
    /* and copy the ptr to hsi */
    fi->hsi[(n+1)] = cm->fi->hsi[n];
  }    

  /* free old cm->fi, but be careful not to free hsi[], we're still pointing to those */
  free(cm->fi->search_opts);
  free(cm->fi->cutoff);
  free(cm->fi->ftype);
  free(cm->fi->hsi);

  cm->fi = fi;
  return eslOK;

 ERROR:
  cm_Fail("cm_AddHMMFilterInfo(), memory allocate error.");
  return status; /* NEVERREACHED */
}  


/* Function: cm_FreeFilterInfo()
 * Date:     EPN, Tue Nov 27 08:48:49 2007
 *
 * Purpose:  Free a FilterInfo_t object corresponding to 
 *           CM <cm>.
 *            
 * Returns:  void
 */
void
cm_FreeFilterInfo(FilterInfo_t *fi, CM_t *cm)
{
  int n;

  for(n = 0; n <= fi->nrounds; n++) if(fi->hsi[n] != NULL) cm_FreeHybridScanInfo(fi->hsi[n], cm); 
  free(fi->search_opts);
  free(fi->cutoff);
  free(fi->ftype);
  free(fi->hsi);

  free(fi);
  return;
}

/* Function: cm_DumpFilterInfo()
 * Date:     EPN, Tue Nov 27 08:50:54 2007
 *
 * Purpose:  Dump a CM's filter info (except hybrid scan info) to stdout. 
 *            
 * Returns:  void.
 */
void
cm_DumpFilterInfo(FilterInfo_t *fi)
{
  int n, v;
  printf("\nFilterInfo summary:\n");
  printf("nrounds: %d\n", fi->nrounds);
  for(n = 0; n <= fi->nrounds; n++) { 
    printf("\nround: %d\n", n);
    if(fi->ftype[n] == FILTER_WITH_HMM)    printf("\ttype: HMM\n"); 
    if(fi->ftype[n] == FILTER_WITH_HYBRID) printf("\ttype: Hybrid\n"); 
    if(fi->ftype[n] == NO_FILTER)          printf("\ttype: CM (no filter)\n"); 
    DumpSearchOpts(fi->search_opts[n]);
    printf("\tcutoff: %10.4f bits\n", fi->cutoff[n]);
    if(fi->hsi[n] != NULL) { 
      printf("\tHybrid info:\n");
      printf("\t\tNumber of sub CM roots: %d\n", fi->hsi[n]->n_v_roots);
      for(v = 0; v < fi->hsi[n]->cm_M; v++) 
	if(fi->hsi[n]->v_isroot[v]) printf("\t\tstate %d is a root\n", v);
    }
  }
  return;
}


/* Function: cm_ValidateFilterInfo()
 * Date:     EPN, Tue Nov 27 09:25:10 2007
 *
 * Purpose:  Validate a Filter Info object.
 *            
 * Returns:  void.
 */
void
cm_ValidateFilterInfo(FilterInfo_t *fi)
{
  int n, sum;
  int do_noqdb;
  int do_hmmonly;
  int do_hmmfilter;
  int do_hmmrescan;
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
  int use_full_cm;

  for(n = 0; n <= fi->nrounds; n++) { 
    do_noqdb       = (fi->search_opts[n] & CM_SEARCH_NOQDB)        ? TRUE : FALSE;
    do_hmmonly     = (fi->search_opts[n] & CM_SEARCH_HMMONLY)      ? TRUE : FALSE;
    do_hmmfilter   = (fi->search_opts[n] & CM_SEARCH_HMMFILTER)    ? TRUE : FALSE;
    do_hbanded     = (fi->search_opts[n] & CM_SEARCH_HBANDED)      ? TRUE : FALSE;
    do_hmmscanbands= (fi->search_opts[n] & CM_SEARCH_HMMSCANBANDS) ? TRUE : FALSE;
    do_sums        = (fi->search_opts[n] & CM_SEARCH_SUMS)         ? TRUE : FALSE;
    do_inside      = (fi->search_opts[n] & CM_SEARCH_INSIDE)       ? TRUE : FALSE;
    do_toponly     = (fi->search_opts[n] & CM_SEARCH_TOPONLY)      ? TRUE : FALSE;
    do_noalign     = (fi->search_opts[n] & CM_SEARCH_NOALIGN)      ? TRUE : FALSE;
    do_null2       = (fi->search_opts[n] & CM_SEARCH_NULL2)        ? TRUE : FALSE;
    do_rsearch     = (fi->search_opts[n] & CM_SEARCH_RSEARCH)      ? TRUE : FALSE;
    do_cmgreedy    = (fi->search_opts[n] & CM_SEARCH_CMGREEDY)     ? TRUE : FALSE;
    do_hmmgreedy   = (fi->search_opts[n] & CM_SEARCH_HMMGREEDY)    ? TRUE : FALSE;
    do_hmmviterbi  = (fi->search_opts[n] & CM_SEARCH_HMMVITERBI)   ? TRUE : FALSE;
    do_hmmforward  = (fi->search_opts[n] & CM_SEARCH_HMMFORWARD)   ? TRUE : FALSE;

    if(n < fi->nrounds) { 
      if(!do_noalign) cm_Fail("cm_ValidateFilterInfo(), round %d has CM_SEARCH_NOALIGN flag down.\n", n);
      if(fi->ftype[n] == FILTER_WITH_HMM) {
	sum = do_noqdb + do_hmmonly + do_hmmfilter + do_hmmrescan + do_hbanded + do_hmmscanbands + do_sums + do_inside + do_toponly + do_null2 + do_rsearch + do_cmgreedy;
	if(sum != 0 || (do_hmmviterbi + do_hmmforward != 1)) { 
	  printf("cm_ValidateFilterInfo(), round %d is FILTER_WITH_HMM but search opts are invalid\n", n);
	  DumpSearchOpts(fi->search_opts[n]);
	  cm_Fail("This is fatal.");
	}
	if(fi->hsi[n] != NULL) cm_Fail("cm_ValidateFilterInfo(), round %d is FILTER_WITH_HMM but hsi[%d] is non-NULL\n", n, n);
      }
      else if (fi->ftype[n] == FILTER_WITH_HYBRID) {
	if(fi->hsi[n] == NULL) cm_Fail("cm_ValidateFilterInfo(), round %d is FILTER_WITH_HYBRID but hsi[%d] is NULL\n", n, n);
	use_full_cm = (fi->hsi[n]->v_isroot[0]) ? TRUE : FALSE;
	if(use_full_cm) { 
	  sum = do_hmmonly + do_hmmfilter + do_hmmrescan + do_hbanded + do_hmmscanbands + do_sums + do_toponly + do_null2 + do_rsearch + do_hmmviterbi + do_hmmforward;	
	  if(sum != 0) {
	    printf("cm_ValidateFilterInfo(), round %d is FILTER_WITH_HYBRID with FULL CM but search opts are invalid\n", n);
	    DumpSearchOpts(fi->search_opts[n]);
	    cm_Fail("This is fatal.");
	  }
	}
	else {
	  sum = do_hmmonly + do_hmmfilter + do_hmmrescan + do_hbanded + do_hmmscanbands + do_sums + do_toponly + do_null2 + do_rsearch + do_cmgreedy;	
	  if(sum != 0 || (do_hmmviterbi + do_hmmforward != 1)) { 
	    printf("cm_ValidateFilterInfo(), round %d is FILTER_WITH_HYBRID with non-full CM but search opts are invalid\n", n);
	    DumpSearchOpts(fi->search_opts[n]);
	    cm_Fail("This is fatal.");
	  }
	}
      }
      else cm_Fail("cm_ValidateFilterInfo(), round %d is neither type FILTER_WITH_HMM nor FILTER_WITH_HYBRID\n", n);
    }
    else { /* round n == fi->nrounds */
      /* check final round, in which we're done filtering */
      if(fi->ftype[fi->nrounds] != NO_FILTER) cm_Fail("cm_ValidateFilterInfo(), final round %d is not NO_FILTER.", n);
      if(fi->hsi[fi->nrounds]   != NULL)      cm_Fail("cm_ValidateFilterInfo(), final round hsi non-NULL.");
      if(do_hmmviterbi || do_hmmforward) { /* searching with only an HMM */
	if(do_hmmviterbi + do_hmmforward != 1) cm_Fail("cm_ValidateFilterInfo(), final round %d specifies HMM viterbi and HMM forward.\n");
	if(do_inside)                          cm_Fail("cm_ValidateFilterInfo(), final round %d specifies HMM viterbi or HMM forward but also Inside.\n");
      }
    }
  }
  ESL_DPRINTF1(("FilterInfo validated.\n"));
  return;
}


/* Function: cm_UpdateFilterInfoCutoff()
 * Date:     EPN, Tue Nov 27 13:43:21 2007
 *
 * Purpose:  Update the cutoff value for a specified round of filtering
 *            
 * Returns:  void, dies if some error
 */
void
cm_UpdateFilterInfoCutoff(CM_t *cm, int nround, float cutoff)
{
  if(cm->fi == NULL)           cm_Fail("cm_UpdateFilterInfoCutoff(), cm->fi is NULL.");
  if(nround > cm->fi->nrounds) cm_Fail("cm_UpdateFilterInfoCutoff(), requested round %d is > cm->fi->nrounds\n", nround, cm->fi->nrounds);
  cm->fi->cutoff[nround] = cutoff;
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
  if(search_opts & CM_SEARCH_HMMONLY)      printf("\tCM_SEARCH_HMMONLY\n");
  if(search_opts & CM_SEARCH_HMMFILTER)    printf("\tCM_SEARCH_HMMFILTER\n");
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

/* CM_TOPHITS: implementation of ranked list of top-scoring hits based
 * closely on HMMER3's P7_TOPHITS but with some complexity removed
 * since we don't have to worry about domains nor with iterative
 * searches (yet) with CM hits, and some extra complexity for bridging
 * the gap with older versions of Infernal.
 *
 * Contents:
 *    1. The CM_TOPHITS object.
 *    2. Standard (human-readable) output of pipeline results.
 *    3. Tabular (parsable) output of pipeline results.
 *    4. Debugging/dev code.
 *    5. Benchmark driver.
 *    6. Test driver.
 *    7. Copyright and license information.
 * 
 * EPN, Tue May 24 13:03:31 2011
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "easel.h"
#include "hmmer.h"
#include "funcs.h"
#include "structs.h"


/*****************************************************************
 * 1. The CM_TOPHITS object
 *****************************************************************/

/* Function:  cm_tophits_Create()
 * Synopsis:  Allocate a hit list.
 * Incept:    EPN, Tue May 24 13:11:06 2011
 *            SRE, Fri Dec 28 07:17:51 2007 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Allocates a new <CM_TOPHITS> hit list and return a pointer
 *            to it.
 *
 * Throws:    <NULL> on allocation failure.
 */
CM_TOPHITS *
cm_tophits_Create(void)
{
  CM_TOPHITS *h = NULL;
  int         default_nalloc = 256;
  int         status;

  ESL_ALLOC(h, sizeof(CM_TOPHITS));
  h->hit    = NULL;
  h->unsrt  = NULL;

  ESL_ALLOC(h->hit,   sizeof(CM_HIT *) * default_nalloc);
  ESL_ALLOC(h->unsrt, sizeof(CM_HIT)   * default_nalloc);
  h->Nalloc    = default_nalloc;
  h->N         = 0;
  h->nreported = 0;
  h->nincluded = 0;
  h->is_sorted_by_score            = TRUE;  /* but only because there's 0 hits */
  h->is_sorted_for_overlap_removal = FALSE; /* actually this is true with 0 hits, but for safety, 
				             * we don't want both sorted_* fields as TRUE */
  h->is_sorted_by_position         = FALSE; /* ditto */
  h->hit[0]    = h->unsrt;  /* if you're going to call it "sorted" when it contains just one hit, you need this */
  return h;

 ERROR:
  cm_tophits_Destroy(h);
  return NULL;
}

/* Function:  cm_tophits_Grow()
 * Synopsis:  Reallocates a larger hit list, if needed.
 * Incept:    EPN, Tue May 24 13:13:52 2011
 *            SRE, Fri Dec 28 07:37:27 2007 [Janelia] (p7_tophits.c)
 *
 * Purpose:   If list <h> cannot hold another hit, doubles
 *            the internal allocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case,
 *            the data in <h> are unchanged.
 */
int
cm_tophits_Grow(CM_TOPHITS *h)
{
  void   *p;
  CM_HIT *ori    = h->unsrt;
  int     Nalloc = h->Nalloc * 2;    /* grow by doubling */
  int     i;
  int     status;

  if (h->N < h->Nalloc) return eslOK; /* we have enough room for another hit */

  ESL_RALLOC(h->hit,   p, sizeof(CM_HIT *) * Nalloc);
  ESL_RALLOC(h->unsrt, p, sizeof(CM_HIT)   * Nalloc);

  /* If we grow a sorted list, we have to translate the pointers
   * in h->hit, because h->unsrt might have just moved in memory. 
   */
  if (h->is_sorted_by_score || h->is_sorted_for_overlap_removal || h->is_sorted_by_position) { 
    for (i = 0; i < h->N; i++)
      h->hit[i] = h->unsrt + (h->hit[i] - ori);
  }

  h->Nalloc = Nalloc;
  return eslOK;

 ERROR:
  return eslEMEM;
}

/* Function:  cm_tophits_CreateNextHit()
 * Synopsis:  Get pointer to new structure for recording a hit.
 * Incept:    EPN, Tue May 24 13:14:12 2011
 *            SRE, Tue Mar 11 08:44:53 2008 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Ask the top hits object <h> to do any necessary
 *            internal allocation and bookkeeping to add a new,
 *            empty hit to its list; return a pointer to 
 *            this new <CM_HIT> structure for data to be filled
 *            in by the caller.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
cm_tophits_CreateNextHit(CM_TOPHITS *h, CM_HIT **ret_hit)
{
  CM_HIT *hit = NULL;
  int     status;

  if ((status = cm_tophits_Grow(h)) != eslOK) goto ERROR;
  
  hit = &(h->unsrt[h->N]);
  h->N++;
  if (h->N >= 2) { 
    h->is_sorted_by_score            = FALSE;
    h->is_sorted_for_overlap_removal = FALSE;
    h->is_sorted_by_position         = FALSE;
  }

  hit->name             = NULL;
  hit->acc              = NULL;
  hit->desc             = NULL;

  hit->start            = 1;
  hit->stop             = 1;
  hit->in_rc            = FALSE;
  hit->score            = 0.0;
  hit->pvalue           = 0.0;
  hit->evalue           = 0.0;

  hit->cm_idx           = -1;
  hit->seq_idx          = -1;
  hit->pass_idx         = -1;

  hit->srcL             = -1;
  hit->maxW             = -1;

  hit->ad               = NULL;
  hit->flags            = CM_HIT_FLAGS_DEFAULT;

  *ret_hit = hit;
  return eslOK;

 ERROR:
  *ret_hit = NULL;
  return status;
}

/* hit_sorter_by_score(), hit_sorter_for_overlap_removal and hit_sorter_by_position: qsort's pawns, below */
static int
hit_sorter_by_score(const void *vh1, const void *vh2)
{
  CM_HIT *h1 = *((CM_HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  CM_HIT *h2 = *((CM_HIT **) vh2);

  if      (h1->score < h2->score) return  1; /* first key, bit score, high to low */
  else if (h1->score > h2->score) return -1;
  else {
    if      (h1->seq_idx > h2->seq_idx) return  1; /* second key, seq_idx (unique id for sequences), low to high */
    else if (h1->seq_idx < h2->seq_idx) return -1;
    else {
      if      (h1->start > h2->start) return  1; /* third key, start position, low to high */
      else if (h1->start < h2->start) return -1;
      else                            return  (h1->pass_idx < h2->pass_idx ? 1 : -1 ); /* fourth key, pass_idx, high to low */
    }
  }
}

static int
hit_sorter_for_overlap_removal(const void *vh1, const void *vh2)
{
  CM_HIT *h1 = *((CM_HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  CM_HIT *h2 = *((CM_HIT **) vh2);

  if      (h1->cm_idx > h2->cm_idx)       return  1; /* first key, cm_idx (unique id for models), low to high */
  else if (h1->cm_idx < h2->cm_idx)       return -1; /* first key, cm_idx (unique id for models), low to high */
  else { 
    if      (h1->seq_idx > h2->seq_idx)   return  1; /* second key, seq_idx (unique id for sequences), low to high */
    else if (h1->seq_idx < h2->seq_idx)   return -1;
    else { 
      /* same sequence, sort by strand, stop position then start position (if revcomp) or start position then stop position (if !revcomp) */
      if     (h1->in_rc > h2->in_rc)      return  1; /* third key, strand (h1->in_rc = 1, h1->in_rc = 0), forward, then reverse */
      else if(h1->in_rc < h2->in_rc)      return -1; /*                   (h1->in_rc = 0, h2->in_rc = 1), forward, then reverse */
      else {
	if     (h1->score < h2->score)    return  1; /* fourth key is bit score, high to low */
	else if(h1->score > h2->score)    return -1; 
	else { 
	  if     (h1->start > h2->start)  return  1; /* fifth key is start position, low to high (irregardless of in_rc value) */
	  else if(h1->start < h2->start)  return -1; 
	  else                            return  0;
	}
      }
    }
  }
}

static int
hit_sorter_by_position(const void *vh1, const void *vh2)
{
  CM_HIT *h1 = *((CM_HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  CM_HIT *h2 = *((CM_HIT **) vh2);

  if      (h1->cm_idx > h2->cm_idx)     return  1; /* first key, cm_idx (unique id for models), low to high */
  else if (h1->cm_idx < h2->cm_idx)     return -1; /* first key, cm_idx (unique id for models), low to high */
  else { 
    if      (h1->seq_idx > h2->seq_idx) return  1; /* second key, seq_idx (unique id for sequences), low to high */
    else if (h1->seq_idx < h2->seq_idx) return -1;
    else { 
      /* same sequence, sort by strand, stop position then start position (if revcomp) or start position then stop position (if !revcomp) */
      if     (h1->in_rc > h2->in_rc)    return  1; /* third key, strand (h1->in_rc = 1, h1->in_rc = 0), forward, then reverse */
      else if(h1->in_rc < h2->in_rc)    return -1; /*                   (h1->in_rc = 0, h2->in_rc = 1), forward, then reverse */
      else if(h1->in_rc) { 
	if     (h1->stop > h2->stop)    return  1; /* both revcomp:     fourth key is stop  position, low to high */
	else if(h1->stop < h2->stop)    return -1; 
	else                            return (h1->start  > h2->start  ? 1 : -1 ); /* both revcomp, same stop position, fourth key is start position, low to high */
      }
      else               {
	if     (h1->start > h2->start)  return  1; /* both !revcomp:    fourth key is start position, low to high */
	else if(h1->start < h2->start)  return -1; 
	else                            return (h1->stop  > h2->stop    ? 1 : -1 ); /* both !revcomp, same start position, fourth key is stop position, low to high */
      }
    }
  }
}

/* Function:  cm_tophits_SortByScore()
 * Synopsis:  Sorts a hit list by bit score.
 * Incept:    EPN, Tue May 24 13:30:23 2011
 *            SRE, Fri Dec 28 07:51:56 2007 (p7_tophits_Sort())
 *
 * Purpose:   Sorts a top hit list by score. After this call,
 *            <h->hit[i]> points to the i'th ranked <CM_HIT> for all
 *            <h->N> hits. First sort key is score (high to low), 
 *            second is seq_idx (low to high), third is 
 *            start position (low to high).

 * Returns:   <eslOK> on success.
 */
int
cm_tophits_SortByScore(CM_TOPHITS *h)
{
  int i;

  if (h->is_sorted_by_score) { 
    h->is_sorted_for_overlap_removal = FALSE;
    h->is_sorted_by_position         = FALSE;
    return eslOK;
  }
  /* initialize hit ptrs, this also unsorts if already sorted by seq_idx */
  for (i = 0; i < h->N; i++) h->hit[i] = h->unsrt + i;
  if (h->N > 1)  qsort(h->hit, h->N, sizeof(CM_HIT *), hit_sorter_by_score);
  h->is_sorted_for_overlap_removal = FALSE;
  h->is_sorted_by_position         = FALSE;
  h->is_sorted_by_score            = TRUE;
  return eslOK;
}

/* Function:  cm_tophits_SortForOverlapRemoval()
 * Synopsis:  Sorts a hit list by cm index, sequence index, strand, then score.
 * Incept:    EPN, Tue Dec 20 09:17:47 2011
 *
 * Purpose:   Sorts a top hit list to ease removal of duplicates.
 *            After this call, <h->hit[i]> points to the i'th ranked
 *            <CM_HIT> for all <h->N> hits. First sort key is cm_idx
 *            (low to high), second is seq_idx (low to high) position
 *            (low to high), third is strand (forward then reverse),
 *            fourth key is score. 
 *
 * Returns:   <eslOK> on success.
 */
int
cm_tophits_SortForOverlapRemoval(CM_TOPHITS *h)
{
  int i;

  if (h->is_sorted_for_overlap_removal) { 
    h->is_sorted_by_score    = FALSE;
    h->is_sorted_by_position = FALSE;
    return eslOK;
  }
  /* initialize hit ptrs, this also unsorts if already sorted by score */
  for (i = 0; i < h->N; i++) h->hit[i] = h->unsrt + i;
  if (h->N > 1)  qsort(h->hit, h->N, sizeof(CM_HIT *), hit_sorter_for_overlap_removal);
  h->is_sorted_by_score            = FALSE;
  h->is_sorted_by_position         = FALSE;
  h->is_sorted_for_overlap_removal = TRUE;

  return eslOK;
}

/* Function:  cm_tophits_SortByPosition()
 * Synopsis:  Sorts a hit list by cm index, sequence index, strand, then position.
 * Incept:    EPN, Tue Dec 20 09:17:47 2011
 *
 * Purpose:   Sorts a top hit list to ease merging of nearby hits after
 *            padding start and stop (like we do at the end of 
 *            cm_pli_SeqCYKFilter()).
 *            After this call, <h->hit[i]> points to the i'th ranked
 *            <CM_HIT> for all <h->N> hits. First sort key is cm_idx
 *            (low to high), second is seq_idx (low to high) position
 *            (low to high), third is strand (forward then reverse),
 *            fourth key is start (if ! in_rc) or stop (if in_rc) 
 *            (low to high for both strands).
 *
 * Returns:   <eslOK> on success.
 */
int
cm_tophits_SortByPosition(CM_TOPHITS *h)
{
  int i;

  if (h->is_sorted_by_position) { 
    h->is_sorted_by_score            = FALSE;
    h->is_sorted_for_overlap_removal = FALSE;
    return eslOK;
  }
  /* initialize hit ptrs, this also unsorts if already sorted by score */
  for (i = 0; i < h->N; i++) h->hit[i] = h->unsrt + i;
  if (h->N > 1)  qsort(h->hit, h->N, sizeof(CM_HIT *), hit_sorter_by_position);
  h->is_sorted_by_score            = FALSE;
  h->is_sorted_for_overlap_removal = FALSE;
  h->is_sorted_by_position         = TRUE;

  return eslOK;
}

/* Function:  cm_tophits_Merge()
 * Synopsis:  Merge two top hits lists.
 * Incept:    EPN, Tue May 24 13:30:39 2011
 *            SRE, Fri Dec 28 09:32:12 2007 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Merge <h2> into <h1>. Upon return, <h1>
 *            contains the unsorted, merged list. <h2>
 *            is effectively destroyed; caller should
 *            not access it further, and may as well free
 *            it immediately.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, and
 *            both <h1> and <h2> remain valid.
 */
int
cm_tophits_Merge(CM_TOPHITS *h1, CM_TOPHITS *h2)
{
  void    *p;
  CM_HIT  *new2;
  int      i;
  int      Nalloc = h1->Nalloc + h2->Nalloc;
  int      status;

  /* Attempt our allocations, so we fail early if we fail. 
   * Reallocating h1->unsrt screws up h1->hit, so fix it.
   */
  ESL_RALLOC(h1->unsrt, p, sizeof(CM_HIT)   * Nalloc);
  ESL_RALLOC(h1->hit,   p, sizeof(CM_HIT *) * Nalloc);

  /* Append h2's unsorted data array to h1. h2's data begin at <new2> */
  new2 = h1->unsrt + h1->N;
  memcpy(new2, h2->unsrt, sizeof(CM_HIT) * h2->N);

  /* h2 now turns over management of name, acc, desc memory to h1;
   * nullify its pointers, to prevent double free.  */
  for (i = 0; i < h2->N; i++)
    {
      h2->unsrt[i].name = NULL;
      h2->unsrt[i].acc  = NULL;
      h2->unsrt[i].desc = NULL;
      h2->unsrt[i].ad   = NULL;
    }

  /* Construct the new grown h1 */
  h1->Nalloc = Nalloc;
  h1->N     += h2->N;
  h1->is_sorted_by_score            = FALSE;
  h1->is_sorted_for_overlap_removal = FALSE;
  h1->is_sorted_by_position         = FALSE;

  /* reset pointers in sorted list (not really nec because we're not sorted) */
  for (i = 0; i < h1->N; i++) h1->hit[i] = h1->unsrt + i;

  return eslOK;
  
 ERROR:
  return status;
}

static int
integer_textwidth(long n)
{
  int w = (n < 0)? 1 : 0;
  while (n != 0) { n /= 10; w++; }
  return w;
}

/* Function:  cm_tophits_GetMaxPositionLength()
 * Synopsis:  Returns maximum position length in hit list (targets).
 * Incept:    EPN, Tue May 24 13:31:13 2011
 *            TJW, Mon May 24 14:16:16 EDT 2010 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Returns the length of the longest hit location (start/end)
 *            of all the registered hits, in chars. This is useful when
 *            deciding how to format output.
 *
 *            The maximum is taken over all registered hits. This
 *            opens a possible side effect: caller might print only
 *            the top hits, and the max name length in these top hits
 *            may be different than the max length over all the hits.
 *
 *            Used specifically for nhmmer output, so expects only one
 *            domain per hit
 *
 *            If there are no hits in <h>, or none of the
 *            hits have names, returns 0.
 */
int
cm_tophits_GetMaxPositionLength(CM_TOPHITS *h)
{
  int i, max;

  max = 0;
  for (i = 0; i < h->N; i++) {
    max = ESL_MAX(max, 
		  ESL_MAX(integer_textwidth(h->unsrt[i].start), integer_textwidth(h->unsrt[i].stop)));
  }
  return max;
}

/* Function:  cm_tophits_GetMaxNameLength()
 * Synopsis:  Returns maximum name length in hit list (targets).
 * Incept:    EPN, Tue May 24 13:34:12 2011
 *            SRE, Fri Dec 28 09:00:13 2007 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Returns the maximum name length of all the registered
 *            hits, in chars. This is useful when deciding how to
 *            format output.
 *            
 *            The maximum is taken over all registered hits. This
 *            opens a possible side effect: caller might print only
 *            the top hits, and the max name length in these top hits
 *            may be different than the max length over all the hits.
 *            
 *            If there are no hits in <h>, or none of the
 *            hits have names, returns 0.
 */
int
cm_tophits_GetMaxNameLength(CM_TOPHITS *h)
{
  int i, max;
  for (max = 0, i = 0; i < h->N; i++)
    if (h->unsrt[i].name != NULL) {
      max = ESL_MAX(max, strlen(h->unsrt[i].name));
    }
  return max;
}

/* Function:  cm_tophits_GetMaxAccessionLength()
 * Synopsis:  Returns maximum accession length in hit list (targets).
 * Incept:    EPN, Tue May 24 13:35:23 2011
 *            SRE, Tue Aug 25 09:18:33 2009 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Same as <cm_tophits_GetMaxNameLength()>, but for
 *            accessions. If there are no hits in <h>, or none
 *            of the hits have accessions, returns 0.
 */
int
cm_tophits_GetMaxAccessionLength(CM_TOPHITS *h)
{
  int i, max, n;
  for (max = 0, i = 0; i < h->N; i++)
    if (h->unsrt[i].acc != NULL) {
      n   = strlen(h->unsrt[i].acc);
      max = ESL_MAX(n, max);
    }
  return max;
}

/* Function:  cm_tophits_GetMaxShownLength()
 * Synopsis:  Returns max shown name/accession length in hit list.
 * Incept:    EPN, Tue May 24 13:35:50 2011
 *            SRE, Tue Aug 25 09:30:43 2009 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Same as <cm_tophits_GetMaxNameLength()>, but 
 *            for the case when --acc is on, where
 *            we show accession if one is available, and 
 *            fall back to showing the name if it is not.
 *            Returns the max length of whatever is being
 *            shown as the reported "name".
 */
int
cm_tophits_GetMaxShownLength(CM_TOPHITS *h)
{
  int i, max, n;
  for (max = 0, i = 0; i < h->N; i++)
  {
      if (h->unsrt[i].acc != NULL && h->unsrt[i].acc[0] != '\0')
    {
      n   = strlen(h->unsrt[i].acc);
      max = ESL_MAX(n, max);
    }
      else if (h->unsrt[i].name != NULL)
    {
      n   = strlen(h->unsrt[i].name);
      max = ESL_MAX(n, max);
    }
  }
  return max;
}

/* Function:  cm_tophits_Reuse()
 * Synopsis:  Reuse a hit list, freeing internals.
 * Incept:    EPN, Tue May 24 13:36:15 2011
 *            SRE, Fri Jun  6 15:39:05 2008 [Janelia] (p7_tophits.c)
 * 
 * Purpose:   Reuse the tophits list <h>; save as 
 *            many malloc/free cycles as possible,
 *            as opposed to <Destroy()>'ing it and
 *            <Create>'ing a new one.
 */
int
cm_tophits_Reuse(CM_TOPHITS *h)
{
  int i;

  if (h == NULL) return eslOK;
  if (h->unsrt != NULL) 
  {
      for (i = 0; i < h->N; i++)
    {
      if (h->unsrt[i].name != NULL) free(h->unsrt[i].name);
      if (h->unsrt[i].acc  != NULL) free(h->unsrt[i].acc);
      if (h->unsrt[i].desc != NULL) free(h->unsrt[i].desc);
      if (h->unsrt[i].ad != NULL)   cm_alidisplay_Destroy(h->unsrt[i].ad);
    }
  }
  h->N         = 0;
  h->is_sorted_by_score            = TRUE;  /* because there's no hits */
  h->is_sorted_for_overlap_removal = FALSE; /* actually this is true with 0 hits, but for safety, 
			           	     * we don't want multiple sorted_* fields as TRUE */
  h->is_sorted_by_position         = FALSE; /* ditto */

  h->hit[0]    = h->unsrt;
  return eslOK;
}

/* Function:  cm_tophits_Destroy()
 * Synopsis:  Frees a hit list.
 * Incept:    EPN, Tue May 24 13:36:49 2011
 *            SRE, Fri Dec 28 07:33:21 2007 [Janelia] (p7_tophits.c)
 */
void
cm_tophits_Destroy(CM_TOPHITS *h)
{
  int i;
  if (h == NULL) return;
  if (h->hit   != NULL) free(h->hit);
  if (h->unsrt != NULL) {
    for (i = 0; i < h->N; i++) {
      if (h->unsrt[i].name != NULL) free(h->unsrt[i].name);
      if (h->unsrt[i].acc  != NULL) free(h->unsrt[i].acc);
      if (h->unsrt[i].desc != NULL) free(h->unsrt[i].desc);
      if (h->unsrt[i].ad != NULL)   cm_alidisplay_Destroy(h->unsrt[i].ad);
    }
    free(h->unsrt);
  }
  free(h);
  return;
}

/* Function:  cm_tophits_CloneHitMostly()
 * Synopsis:  Add a new hit to <dest_th>, a clone of hit <h> in <src_th>.
 * Incept:    EPN, Wed May 25 08:17:35 2011
 *
 * Purpose:   Create a new hit in the CM_TOPHITS object <dest_th>
 *            and copy the information from hit <h> in the sorted
 *            hitlist <src_th> into it.
 * 
 *            NOTE: we do not copy name, acc, desc, or ad.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
cm_tophits_CloneHitMostly(CM_TOPHITS *src_th, int h, CM_TOPHITS *dest_th)
{
  CM_HIT *hit = NULL;
  int     status;

  if ((status = cm_tophits_CreateNextHit(dest_th, &hit)) != eslOK) goto ERROR;
  hit->cm_idx   = src_th->hit[h]->cm_idx;
  hit->seq_idx  = src_th->hit[h]->seq_idx;
  hit->pass_idx = src_th->hit[h]->pass_idx;
  hit->start    = src_th->hit[h]->start;
  hit->stop     = src_th->hit[h]->stop;
  hit->in_rc    = src_th->hit[h]->in_rc;
  hit->root     = src_th->hit[h]->root;
  hit->mode     = src_th->hit[h]->mode;
  hit->score    = src_th->hit[h]->score;
  hit->pvalue   = src_th->hit[h]->pvalue;
  hit->evalue   = src_th->hit[h]->evalue;
  hit->srcL     = src_th->hit[h]->srcL;
  hit->maxW     = src_th->hit[h]->maxW;
  hit->flags    = src_th->hit[h]->flags;
  hit->ad       = NULL;

  return eslOK;

 ERROR:
  return status;
}  

/* Function:  cm_tophits_ComputeEvalues()
 * Synopsis:  Compute E-values given effective database size.
 * Incept:    EPN, Tue Sep 28 05:26:20 2010
 *
 * Purpose:   After the cm pipeline has completed, the CM_TOPHITS object
 *            contains objects with p-values that haven't yet been
 *            converted to e-values.
 *
 * Returns:   <eslOK> on success.
 */
int
cm_tophits_ComputeEvalues(CM_TOPHITS *th, double eZ, int istart)
{
  int i; 

  for (i = istart; i < th->N ; i++) { 
    th->unsrt[i].evalue = th->unsrt[i].pvalue * eZ;
  }
  return eslOK;
}

/* Helper functions for removing overlapping hits from a set of hits
 * that all use the same model, are to the same sequence and are on
 * the same strand. The hits in the set are sorted by score. 
 * 
 * remove_overlaps_one_seq_fast(): O(N) with number of hits N, but
 * less memory efficient (O(L), requiring an allocation of a char
 * array the length of the sequence). This is how RSEARCH and infernal
 * 1.0.2 did overlap removal. 
 * 
 * remove_overlaps_one_seq_memeff(): O(N^2), but O(1) memory (not
 * counting the O(N) memory hit list). This is fast for small numbers
 * of sequences; up to about 10k seqs (10,000 hits takes about ~0.3
 * seconds with benchmark_cm_tophits, but 100,000 hits takes 43s).
 * 
 * cm_tophits_RemoveOverlaps() picks a function to call based on 
 * length of the sequence L and number of hits N.
 *
 * Both functions: 
 * Returns eslOK on success
 * Returns eslEINVAL if not all hits in the set have equal srcL, cm_idx, seq_idx, or in_rc values
 * Returns eslERANGE if a hit includes positions outside of 1..srcL
 * Returns eslEMEM if out of memory (remove_overlaps_one_seq_fast() only)
 * errbuf is filled if not returning eslOK
 */
int remove_overlaps_one_seq_fast(CM_TOPHITS *th, int64_t idx1, int64_t idx2, char *errbuf)
{ 
  int     status;
  int64_t i;
  char   *covered = NULL;  /* [1..pos..srcL] is position pos covered by a hit we've examined and kept? */
  int     remove_flag = FALSE;
  int64_t srcL    = th->hit[idx1]->srcL;
  int64_t cm_idx  = th->hit[idx1]->cm_idx;
  int64_t seq_idx = th->hit[idx1]->seq_idx;
  int     in_rc   = th->hit[idx1]->in_rc;
  int64_t pos, min, max; /* position indices in the sequence */

  /*printf("in remove_overlaps_one_seq_fast() i: %" PRId64 " j: %" PRId64 "\n", idx1, idx2);*/

  ESL_ALLOC(covered, sizeof(char) * (srcL+1));
  for(i = 0; i <= srcL; i++) covered[i] = FALSE;

  for(i = idx1; i <= idx2; i++) { 
    /* verify that what we think is true is true */
    if(th->hit[i]->srcL    != srcL)                       ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, srcL inconsistent, hit %" PRId64, i);
    if(th->hit[i]->cm_idx  != cm_idx)                     ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, cm_idx is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->seq_idx != seq_idx)                    ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, seq_idx is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->in_rc   != in_rc)                      ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, in_rc is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->start < 1 || th->hit[i]->start > srcL) ESL_FAIL(eslERANGE, errbuf, "removing overlapping hits, start posn is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->stop  < 1 || th->hit[i]->stop  > srcL) ESL_FAIL(eslERANGE, errbuf, "removing overlapping hits, stop posn is inconsistent, hit %" PRId64, i);

    if(! (th->hit[i]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) { 
      /* i is not a duplicate that's already been removed */
      min = ESL_MIN(th->hit[i]->start, th->hit[i]->stop); 
      max = ESL_MAX(th->hit[i]->start, th->hit[i]->stop); 
      remove_flag = FALSE;
      for(pos = min; pos <= max; pos++) { 
	if(covered[pos] == TRUE) { remove_flag = TRUE; break; }
      }

      if(remove_flag) { 
	th->hit[i]->flags |=  CM_HIT_IS_REMOVED_DUPLICATE;
	th->hit[i]->flags &= ~CM_HIT_IS_REPORTED;  /* could be set if pli->use_bit_cutoffs (--cut_ga, --cut_nc, --cut_tc) */
	th->hit[i]->flags &= ~CM_HIT_IS_INCLUDED;  /* could be set if pli->use_bit_cutoffs (--cut_ga, --cut_nc, --cut_tc) */
      }
      else { 
	for(pos = min; pos <= max; pos++) covered[pos] = TRUE;
      }
    }
  }

  if(covered != NULL) free(covered);
  return eslOK;

 ERROR:
  if(covered != NULL) free(covered);
  ESL_FAIL(status, errbuf, "removing overlapping hits, out of memory");
  return status; /* NOT REACHED */
}

int remove_overlaps_one_seq_memeff(CM_TOPHITS *th, int64_t idx1, int64_t idx2, char *errbuf)
{
  int64_t i, j;            
  int64_t srcL    = th->hit[idx1]->srcL;
  int64_t cm_idx  = th->hit[idx1]->cm_idx;
  int64_t seq_idx = th->hit[idx1]->seq_idx;
  int     in_rc   = th->hit[idx1]->in_rc;

  /*printf("in remove_overlaps_one_seq_memeff() i: %" PRId64 " j: %" PRId64 "\n", idx1, idx2);*/

  for(i = idx1; i <= idx2; i++) { 
    /* verify that what we think is true is true */
    if(th->hit[i]->srcL    != srcL)                       ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, srcL inconsistent, hit %" PRId64, i);
    if(th->hit[i]->cm_idx  != cm_idx)                     ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, cm_idx is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->seq_idx != seq_idx)                    ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, seq_idx is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->in_rc   != in_rc)                      ESL_FAIL(eslEINVAL, errbuf, "removing overlapping hits, in_rc is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->start < 1 || th->hit[i]->start > srcL) ESL_FAIL(eslERANGE, errbuf, "removing overlapping hits, start posn is inconsistent, hit %" PRId64, i);
    if(th->hit[i]->stop  < 1 || th->hit[i]->stop  > srcL) ESL_FAIL(eslERANGE, errbuf, "removing overlapping hits, stop posn is inconsistent, hit %" PRId64, i);

    if(! (th->hit[i]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) { 
      /* i is not a duplicate that's already been removed */
      for(j = i+1; j <= idx2; j++) { 
	if(! (th->hit[j]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) { 
	  /* j has not already been removed */
	  /*printf("comparing %" PRId64 " and %" PRId64 "\n", i, j);*/
	  if(th->hit[j]->in_rc == FALSE &&                 /* both i and j are on forward strand */
	     (! (th->hit[j]->stop < th->hit[i]->start)) && /* one of two ways in which i and j DO NOT overlap */
	     (! (th->hit[i]->stop < th->hit[j]->start))) { /* the other way   in which i and j DO NOT overlap */
	    /* i and j overlap, i is better scoring so remove j */
	    th->hit[j]->flags |=  CM_HIT_IS_REMOVED_DUPLICATE;
	    th->hit[j]->flags &= ~CM_HIT_IS_REPORTED;  /* could be set if pli->use_bit_cutoffs (--cut_ga, --cut_nc, --cut_tc) */
	    th->hit[j]->flags &= ~CM_HIT_IS_INCLUDED;  /* could be set if pli->use_bit_cutoffs (--cut_ga, --cut_nc, --cut_tc) */
	    /*printf("\tremoved j (FWD) %5" PRId64 "..%5" PRId64 " overlaps with %5" PRId64 "..%5" PRId64 "\n", th->hit[i]->start, th->hit[i]->stop, th->hit[j]->start, th->hit[j]->stop);*/
	  }
	  else if(th->hit[j]->in_rc == TRUE &&                  /* both i and j are on reverse strand */
		  (! (th->hit[j]->start < th->hit[i]->stop)) && /* one of two ways in which i and j DO NOT overlap */
		  (! (th->hit[i]->start < th->hit[j]->stop))) { /* the other way   in which i and j DO NOT overlap */
	    /* i and j overlap, i is better scoring so remove j */
	    th->hit[j]->flags |=  CM_HIT_IS_REMOVED_DUPLICATE;
	    th->hit[j]->flags &= ~CM_HIT_IS_REPORTED;  /* could be set if pli->use_bit_cutoffs (--cut_ga, --cut_nc, --cut_tc) */
	    th->hit[j]->flags &= ~CM_HIT_IS_INCLUDED;  /* could be set if pli->use_bit_cutoffs (--cut_ga, --cut_nc, --cut_tc) */
	    /*printf("\tremoved j (REV) %5" PRId64 "..%5" PRId64 " overlaps with %5" PRId64 "..%5" PRId64 "\n", th->hit[i]->start, th->hit[i]->stop, th->hit[j]->start, th->hit[j]->stop);*/
	  }
	}
      }
    }
  }
  return eslOK;
}
/* Function:  cm_tophits_RemoveOverlaps()
 * Synopsis:  Remove overlapping hits from a tophits object sorted by seq_idx.
 * Incept:    EPN, Tue Jun 14 05:42:31 2011
 *
 * Purpose:   After the CM pipeline has completed, the CM_TOPHITS
 *            object may contain overlapping hits if the target was
 *            broken into overlapping windows. Scan through the
 *            tophits object, which, upon entering, is sorted by cm
 *            index, sequence index, strand, and score and remove hits
 *            greedily. For each hit for the same model, sequence and
 *            strand, remove all lower scoring hits that overlap with
 *            it. This is done by one of two helper functions
 *            (remove_overlaps_one_seq_fast() or
 *            remove_overlaps_one_seq_memeff() described near their
 *            definition above) depending on the length of the 
 *            sequence and the number of hits.
 *         
 * Returns:   eslOK on success. 
 *            eslEINVAL if th is not sorted appropriately, errbuf filled
 *            eslERANGE if a hit includes positions outside of 1..srcL, errbuf filled
 *            eslEMEM if we run out of memory, errbuf filled
 */
int
cm_tophits_RemoveOverlaps(CM_TOPHITS *th, char *errbuf)
{
  int status;
  int64_t i, j;            
  int64_t nhits  = 0;

  if (! th->is_sorted_for_overlap_removal) ESL_FAIL(status, errbuf, "cm_tophits_RemoveOverlaps() list is not sorted appropriately");
  if (th->N<2) return eslOK;

  i = 0;
  while(i < th->N) { 
    j = i+1;
    while(j < th->N && 
	  th->hit[j]->cm_idx  == th->hit[i]->cm_idx  &&
	  th->hit[j]->seq_idx == th->hit[i]->seq_idx &&
	  th->hit[j]->in_rc   == th->hit[i]->in_rc) { 
      j++;
    }
    if(j != (i+1)) { 
      /* Hits i to j-1 form a set of hits that all share cm_idx, seq_idx
       * and in_rc. Remove overlaps from this set in 1 of 2 ways,
       * depending on length of source sequence and number of hits.
       * remove_overlaps_one_seq_fast() will need to allocate a char
       * array of size srcL, but it is significantly faster when nhits
       * is big.
       */ 
      nhits = (j-1)-i+1;
      if(nhits < 5000 || th->hit[i]->srcL > 256000000) { 
	if((status = remove_overlaps_one_seq_memeff(th, i, j-1, errbuf)) != eslOK) return status;
      }
      else { /* use fast, non mem-efficient way if >= 5000 hits and L is < 256 Mb */
	if((status = remove_overlaps_one_seq_fast  (th, i, j-1, errbuf)) != eslOK) return status;
      }
    }
    i = j; /* skip ahead to begin next set */
  }
  /*printf("Leaving cm_tophits_RemoveOverlaps()\n");
    cm_tophits_Dump(stdout, th);*/

  return eslOK;
}

/* Function:  cm_tophits_UpdateHitPositions()
 * Synopsis:  Update sequence positions in a hit list.
 * Incept:    EPN, Wed May 25 09:12:52 2011
 *
 * Purpose: For hits <hit_start..th->N-1> in a hit list, update the
 *          positions of start, stop, ad->sqto and ad->sqfrom.  This
 *          is necessary when we've searched a chunk of subsequence
 *          that originated somewhere within a larger sequence.
 *          
 *
 * Returns: <eslOK> on success.
 */
 int
 cm_tophits_UpdateHitPositions(CM_TOPHITS *th, int hit_start, int64_t seq_start, int in_revcomp)
{ 
  int i;
  if(in_revcomp) { 
    for (i = hit_start; i < th->N ; i++) {
      th->unsrt[i].start = seq_start - th->unsrt[i].start + 1;
      th->unsrt[i].stop  = seq_start - th->unsrt[i].stop  + 1;
      th->unsrt[i].in_rc = TRUE;
      if(th->unsrt[i].ad != NULL) { 
	th->unsrt[i].ad->sqfrom = seq_start - th->unsrt[i].ad->sqfrom + 1;
	th->unsrt[i].ad->sqto   = seq_start - th->unsrt[i].ad->sqto + 1;
      }
    }

  }
  else { 
    for (i = hit_start; i < th->N ; i++) {
      th->unsrt[i].start += seq_start-1;
      th->unsrt[i].stop  += seq_start-1;
      th->unsrt[i].in_rc = FALSE;
      if(th->unsrt[i].ad != NULL) { 
	th->unsrt[i].ad->sqfrom += seq_start-1;
	th->unsrt[i].ad->sqto   += seq_start-1;
      }
    }
  }
  return eslOK;
}

/* Function:  cm_tophits_SetSourceLengths()
 * Synopsis:  Update hit->srcL for all hits in a hitlist.
 * Incept:    EPN, Wed Nov 23 14:57:53 2011
 *
 * Purpose: For all hits in a hitlist, set the
 *          srcL value (full length of source sequence the hit originated
 *          in) given an array of <nseqs> source lengths indexed by seq_idx.
 *
 * Returns: <eslOK> on success.
 *          <eslEINVAL> if hit->seq_idx >= nseqs (we have a hit to a sequence
 *          we don't know the length of - something went wrong)
 */
 int
 cm_tophits_SetSourceLengths(CM_TOPHITS *th, int64_t *srcL, uint64_t nseqs)
{ 
  int i;
  for (i = 0; i < th->N ; i++) { 
    if(th->unsrt[i].seq_idx >= nseqs) return eslEINVAL;
    th->unsrt[i].srcL = srcL[th->unsrt[i].seq_idx];
  }  

  return eslOK;
}

/*---------------- end, CM_TOPHITS object -----------------------*/

/*****************************************************************
 * 2. Standard (human-readable) output of pipeline results
 *****************************************************************/


/* Function:  cm_tophits_Threshold()
 * Synopsis:  Apply score and E-value thresholds to a hitlist before output.
 * Incept:    EPN, Tue May 24 13:44:13 2011
 *            SRE, Tue Dec  9 09:04:55 2008 [Janelia] (p7_tophits.c)
 *
 * Purpose:   After a pipeline has completed, go through it and mark all
 *            the targets and domains that are "significant" (satisfying
 *            the reporting thresholds set for the pipeline). 
 *            
 *            Also sets the final total number of reported and
 *            included targets, the number of reported and included
 *            targets in each target, and the size of the search space
 *            for per-domain conditional E-value calculations,
 *            <pli->domZ>. By default, <pli->domZ> is the number of
 *            significant targets reported.
 *
 *            If model-specific thresholds were used in the pipeline,
 *            we cannot apply those thresholds now. They were already
 *            applied in the pipeline. In this case all we're
 *            responsible for here is counting them (setting
 *            nreported, nincluded counters).
 *            
 * Returns:   <eslOK> on success.
 */
int
cm_tophits_Threshold(CM_TOPHITS *th, CM_PIPELINE *pli)
{
  int h;   /* counter */
  
  /* Flag reported, included targets (if we're using general thresholds) */
  if (! pli->use_bit_cutoffs) {
    for (h = 0; h < th->N; h++) {
      if (! (th->hit[h]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) { 
	if (cm_pli_TargetReportable(pli, th->hit[h]->score, th->hit[h]->evalue)) {
	  th->hit[h]->flags |= CM_HIT_IS_REPORTED;
	  if (cm_pli_TargetIncludable(pli, th->hit[h]->score, th->hit[h]->evalue))
	    th->hit[h]->flags |= CM_HIT_IS_INCLUDED;
	}
      }
    }
  }

  /* Count reported, included targets */
  th->nreported = 0;
  th->nincluded = 0;
  for (h = 0; h < th->N; h++) { 
    if (! (th->hit[h]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) { 
      if (th->hit[h]->flags & CM_HIT_IS_REPORTED)  th->nreported++;
      if (th->hit[h]->flags & CM_HIT_IS_INCLUDED)  th->nincluded++;
    }
  }
  
  return eslOK;
}


/* Function:  cm_tophits_Targets()
 * Synopsis:  Standard output format for a top target hits list.
 * Incept:    EPN, Tue May 24 13:48:29 2011
 *            SRE, Tue Dec  9 09:10:43 2008 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Output a list of the reportable top target hits in <th> 
 *            in human-readable ASCII text format to stream <ofp>, using
 *            final pipeline accounting stored in <pli>. 
 * 
 *            Hits must have a valid alignment display in order
 *            to determine if they're truncated or not. 
 *
 *            The tophits list <th> should already be sorted (see
 *            <cm_tophits_Sort()> and thresholded (see
 *            <cm_tophits_Threshold>).
 *
 * Returns:   <eslOK> on success.
 */
int
cm_tophits_Targets(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli, int textw)
{
  int    status;
  int    h,i;
  int    namew;
  int    posw;
  int    descw;
  int    rankw;
  char  *showname;
  int    have_printed_incthresh = FALSE;
  int    nprinted = 0;

  char *namestr     = NULL;
  char *posstr      = NULL;
  char *rankstr     = NULL;
  char *cur_rankstr = NULL;

  /* when --acc is on, we'll show accession if available, and fall back to name */
  if (pli->show_accessions) namew = ESL_MAX(8, cm_tophits_GetMaxShownLength(th));
  else                      namew = ESL_MAX(8, cm_tophits_GetMaxNameLength(th));

  posw = ESL_MAX(6, cm_tophits_GetMaxPositionLength(th));

  if (textw >  0) descw = ESL_MAX(32, textw - namew - (2*posw) - 31); /* 31 chars excluding desc is from the format: 1 + 2 + 9+2(E-value) + 6+2(score) + 1+2(strand) + 2+2+2+2(spacing) */
  else            descw = 0;                                          /* unlimited desc length is handled separately */

  rankw = ESL_MAX(4, (integer_textwidth(th->nreported)+2));

  ESL_ALLOC(namestr,     sizeof(char) * (namew+1));
  ESL_ALLOC(posstr,      sizeof(char) * (posw+1));
  ESL_ALLOC(rankstr,     sizeof(char) * (rankw+1));
  ESL_ALLOC(cur_rankstr, sizeof(char) * (rankw+1));

  for(i = 0; i < namew; i++) { namestr[i] = '-'; } namestr[namew] = '\0';
  for(i = 0; i < posw;  i++) { posstr[i]  = '-'; } posstr[posw]   = '\0';
  for(i = 0; i < rankw; i++) { rankstr[i] = '-'; } rankstr[rankw] = '\0';
  cur_rankstr[rankw] = '\0';

  fprintf(ofp, "Hit scores:\n");
  fprintf(ofp, " %*s  %9s  %6s  %-*s  %*s  %*s  %1s  %5s  %4s  %s\n", rankw, "rank",  "E-value",   " score", namew, (pli->mode == CM_SEARCH_SEQS ? "sequence":"model"), posw, "start", posw, "end", "", "trunc", "pass", "description");
  fprintf(ofp, " %*s  %9s  %6s  %-*s  %*s  %*s  %1s  %5s  %4s  %s\n", rankw, rankstr, "---------", "------", namew, namestr, posw, posstr, posw, posstr, "", "-----", "----", "-----------");
  
  nprinted = 0;
  for (h = 0; h < th->N; h++) { 
    if (th->hit[h]->flags & CM_HIT_IS_REPORTED) { 

      if (! (th->hit[h]->flags & CM_HIT_IS_INCLUDED) && ! have_printed_incthresh) {
	fprintf(ofp, "  ------ inclusion threshold ------\n");
	have_printed_incthresh = TRUE;
      }
      
      if (pli->show_accessions) { 
	/* the --acc option: report accessions rather than names if possible */
	if (th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0') showname = th->hit[h]->acc;
	else                                                       showname = th->hit[h]->name;
      }
      else {
	showname = th->hit[h]->name;
      }
      
      sprintf(cur_rankstr, "(%d)", nprinted+1);

      fprintf(ofp, " %*s  %9.2g  %6.1f  %-*s  %*" PRId64 "  %*" PRId64 "  %c  %5s  %4d  ",
	      rankw, cur_rankstr,
	      th->hit[h]->evalue,
	      th->hit[h]->score,
	      namew, showname,
	      posw, th->hit[h]->start,
	      posw, th->hit[h]->stop,
	      (th->hit[h]->start < th->hit[h]->stop ? '+' : '-'), 
	      cm_alidisplay_TruncString(th->hit[h]->ad), 
	      th->hit[h]->pass_idx);
      
      if (textw > 0) fprintf(ofp, "%-.*s\n", descw, th->hit[h]->desc == NULL ? "" : th->hit[h]->desc);
      else           fprintf(ofp, "%s\n",           th->hit[h]->desc == NULL ? "" : th->hit[h]->desc);
      /* do NOT use *s with unlimited (INT_MAX) line length. Some systems
       * have an fprintf() bug here (we found one on an Opteron/SUSE Linux
       * system (#h66))
       */
      nprinted++;
    }
  }
  if (th->nreported == 0) fprintf(ofp, "\n   [No hits detected that satisfy reporting thresholds]\n");

  if(namestr     != NULL) free(namestr);
  if(posstr      != NULL) free(posstr);
  if(rankstr     != NULL) free(rankstr);
  if(cur_rankstr != NULL) free(cur_rankstr);

  return eslOK;

 ERROR:
  if(namestr != NULL) free(namestr);
  if(posstr  != NULL) free(posstr);

  return status;
}

/* Function:  cm_tophits_HitAlignments()
 * Synopsis:  Standard output format for alignments.
 * Incept:    EPN, Tue May 24 13:50:53 2011
 *            SRE, Tue Dec  9 09:32:32 2008 [Janelia] (p7_tophits.c::p7_tophits_Domains())
 *
 * Purpose:   For each reportable target sequence, output a tabular summary
 *            of each hit followed by its alignment.
 * 
 *            Similar to <cm_tophits_Targets()>; see additional notes there.
 */
int
cm_tophits_HitAlignments(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli, int textw)
{
  int status;
  int h, i, nprinted;
  int namew, descw, rankw;
  char *showname;
  char *rankstr     = NULL;
  char *cur_rankstr = NULL;

  /* next 4 characters indicate whether alignment ends internal to model/sequence
   * and if it is in a truncated alignment mode. */
  char lmod;
  char rmod;
  char lseq; 
  char rseq;

  rankw = ESL_MAX(4, (integer_textwidth(th->nreported)+2));
  ESL_ALLOC(rankstr,     sizeof(char) * (rankw+1));
  ESL_ALLOC(cur_rankstr, sizeof(char) * (rankw+1));
  for(i = 0; i < rankw; i++) { rankstr[i] = '-'; } rankstr[rankw] = '\0';
  cur_rankstr[rankw] = '\0';

  fprintf(ofp, "Hit alignments:\n");

  nprinted = 0;
  for (h = 0; h < th->N; h++) {
    if (th->hit[h]->flags & CM_HIT_IS_REPORTED) {
      if (pli->show_accessions && th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0') {
	showname = th->hit[h]->acc;
	namew    = strlen(th->hit[h]->acc);
      }
      else {
	showname = th->hit[h]->name;
	namew = strlen(th->hit[h]->name);
      }

      if (textw > 0) {
	descw = ESL_MAX(32, textw - namew - 5);
	fprintf(ofp, ">> %s  %-.*s\n", showname, descw, (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc));
      }
      else {
	fprintf(ofp, ">> %s  %s\n",    showname,        (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc));
      }

      /* The hit info display is 101+rankw char wide:, where rankw is the maximum of 4 and 2 plus the number of digits in th->N. 
       * If (! pli->show_alignments) the width drops to 72+rankw.
       *     rank  score   E-value cm from   cm to       seq from      seq to      trunc  acc  bands   mx Mb aln secs  
       *     ---- ------ --------- ------- -------    ----------- -----------      ----- ---- ------ ------- --------
       *      (1)  123.4    6.8e-9       3      72 []         412         492 + ..       0.98    yes    1.30     0.04  
       *     (12)  123.4    1.8e-3       1      72 []         180         103 - ..       0.90    yes    0.65     2.23  
       *    rankw 123456 123456789 1234567 1234567 12 12345678901 12345678901 1 12 12345 1234  12345 1234567 12345678
       *        0         1         2         3         4         5         6         7         8        9        10
       *        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345789012345678901
       * 
       * In rare cases, when computing posteriors is not feasible in allowable memory, 
       * the "acc" column will be replaced by a "cyksc" colum which is 6 characters wide 
       * instead of 4. 
       */

      fprintf(ofp, " %*s %1s %6s %9s %7s %7s %2s %11s %11s %1s %2s %5s",  rankw, "rank", "", "score", "E-value", "cm from", "cm to", "", "seq from", "seq to", "", "", "trunc");
      if(th->hit[h]->ad->ppline) { fprintf(ofp, " %4s %5s %7s %7s", "acc",   "bands", "mx Mb",   "seconds"); }
      else                       { fprintf(ofp, " %6s %5s %7s %7s", "cyksc", "bands", "mx Mb",   "seconds"); }

      fprintf(ofp, "\n");

      fprintf(ofp, " %*s %1s %6s %9s %7s %7s %2s %11s %11s %1s %2s %5s",  rankw, rankstr,  "", "------", "---------", "-------", "-------", "", "-----------", "-----------", "", "", "-----");
      if(th->hit[h]->ad->ppline) { fprintf(ofp, " %4s %5s %7s %7s", "----",   "-----", "-------", "-------"); }
      else                       { fprintf(ofp, " %6s %5s %7s %7s", "------", "-----", "-------", "-------"); }
      fprintf(ofp, "\n");

      if(cm_alidisplay_Is5PTrunc(th->hit[h]->ad)) { /* 5' truncated */
	lmod = '~';
	if(th->hit[h]->in_rc) { lseq = th->hit[h]->ad->sqfrom == th->hit[h]->srcL ? '{' : '~'; }
	else                  { lseq = th->hit[h]->ad->sqfrom == 1                ? '{' : '~'; }
      }
      else { /* not 5' truncated */
	lmod = th->hit[h]->ad->cfrom_emit == 1 ? '[' : '.';
	if(th->hit[h]->in_rc) { lseq = th->hit[h]->ad->sqfrom == th->hit[h]->srcL ? '[' : '.'; }
	else                  { lseq = th->hit[h]->ad->sqfrom == 1                ? '[' : '.'; }
      }
      if(cm_alidisplay_Is3PTrunc(th->hit[h]->ad)) { /* 3' truncated */
	rmod = '~';
	if(th->hit[h]->in_rc) { rseq = th->hit[h]->ad->sqto == 1                ? '}' : '~'; }
	else                  { rseq = th->hit[h]->ad->sqto == th->hit[h]->srcL ? '}' : '~'; }
      }	
      else { /* not 3' truncated */
	rmod = th->hit[h]->ad->cto_emit == th->hit[h]->ad->clen ? ']' : '.';
	if(th->hit[h]->in_rc) { rseq = th->hit[h]->ad->sqto == 1                ? ']' : '.'; }
	else                  { rseq = th->hit[h]->ad->sqto == th->hit[h]->srcL ? ']' : '.'; }
      }

      sprintf(cur_rankstr, "(%d)", nprinted+1);

      fprintf(ofp, " %*s %c %6.1f %9.2g %7d %7d %c%c %11" PRId64 " %11" PRId64 " %c %c%c %5s",
	      rankw, cur_rankstr,
	      (th->hit[h]->flags & CM_HIT_IS_INCLUDED ? '!' : '?'),
	      th->hit[h]->score,
	      th->hit[h]->evalue,
	      th->hit[h]->ad->cfrom_emit,
	      th->hit[h]->ad->cto_emit,
	      lmod, rmod, 
	      th->hit[h]->start,
	      th->hit[h]->stop,
	      (th->hit[h]->start < th->hit[h]->stop ? '+' : '-'),
	      lseq, rseq, 
	      cm_alidisplay_TruncString(th->hit[h]->ad));
      if(th->hit[h]->ad->ppline) { fprintf(ofp, " %4.2f", th->hit[h]->ad->avgpp); }
      else                       { fprintf(ofp, " %6.2f", th->hit[h]->ad->sc); }
      fprintf(ofp, " %5s %7.2f %7.2f\n\n",	
	      (th->hit[h]->ad->used_hbands ? "yes" : "no"),
	      th->hit[h]->ad->matrix_Mb,
	      th->hit[h]->ad->elapsed_secs);
      /*cm_alidisplay_Dump(ofp, th->hit[h]->ad);*/
      cm_alidisplay_Print(ofp, th->hit[h]->ad, 40, textw, pli->show_accessions, TRUE);
      fprintf(ofp, "\n");
      nprinted++;
    }
  }
  if (th->nreported == 0) { 
    fprintf(ofp, "\n   [No hits detected that satisfy reporting thresholds]\n"); 
  }

  if(rankstr     != NULL) free(rankstr);
  if(cur_rankstr != NULL) free(cur_rankstr);
  return eslOK;

ERROR: 
  if(rankstr     != NULL) free(rankstr);
  if(cur_rankstr != NULL) free(cur_rankstr);
  return status;

}
  
/* Function:  cm_tophits_HitAlignmentStatistics()
 * Synopsis:  Final alignment statistics output from a list of hits.
 * Incept:    EPN, Thu Jun 16 04:22:11 2011
 *
 * Purpose:   Print summary statistics on alignments performed 
 *            for all hits in tophits object <th> to stream <ofp>.
 *            Statistics are compiled and reported for all 
 *            reported and included hits as well as for each
 *            type independently.
 * 
 *            <used_cyk> is TRUE if pli->align_cyk was TRUE, if 
 *            HMM banded CYK was used instead of HMM banded 
 *            optimal accuracy for alignment.
 *
 * Returns:   <eslOK> on success.
 */
int
cm_tophits_HitAlignmentStatistics(FILE *ofp, CM_TOPHITS *th, int used_cyk)
{
  uint64_t h;
  int is_reported;
  int is_included;

  int64_t rep_naln_hb, rep_naln_dccyk;
  double  tot_rep_matrix_Mb_hb, tot_rep_matrix_Mb_dccyk;
  double  min_rep_matrix_Mb_hb, min_rep_matrix_Mb_dccyk;
  double  max_rep_matrix_Mb_hb, max_rep_matrix_Mb_dccyk;
  double  avg_rep_matrix_Mb_hb, avg_rep_matrix_Mb_dccyk;
  double  tot_rep_elapsed_secs_hb, tot_rep_elapsed_secs_dccyk;
  double  min_rep_elapsed_secs_hb, min_rep_elapsed_secs_dccyk;
  double  max_rep_elapsed_secs_hb, max_rep_elapsed_secs_dccyk;
  double  avg_rep_elapsed_secs_hb, avg_rep_elapsed_secs_dccyk;

  int64_t inc_naln_hb, inc_naln_dccyk;
  double  tot_inc_matrix_Mb_hb, tot_inc_matrix_Mb_dccyk;
  double  min_inc_matrix_Mb_hb, min_inc_matrix_Mb_dccyk;
  double  max_inc_matrix_Mb_hb, max_inc_matrix_Mb_dccyk;
  double  avg_inc_matrix_Mb_hb, avg_inc_matrix_Mb_dccyk;
  double  tot_inc_elapsed_secs_hb, tot_inc_elapsed_secs_dccyk;
  double  min_inc_elapsed_secs_hb, min_inc_elapsed_secs_dccyk;
  double  max_inc_elapsed_secs_hb, max_inc_elapsed_secs_dccyk;
  double  avg_inc_elapsed_secs_hb, avg_inc_elapsed_secs_dccyk;

  /* initalize */
  rep_naln_hb = rep_naln_dccyk = 0;
  tot_rep_matrix_Mb_hb = tot_rep_matrix_Mb_dccyk = 0.;
  min_rep_matrix_Mb_hb = min_rep_matrix_Mb_dccyk = 0.;
  max_rep_matrix_Mb_hb = max_rep_matrix_Mb_dccyk = 0.;
  tot_rep_elapsed_secs_hb = tot_rep_elapsed_secs_dccyk = 0.;
  min_rep_elapsed_secs_hb = min_rep_elapsed_secs_dccyk = 0.;
  max_rep_elapsed_secs_hb = max_rep_elapsed_secs_dccyk = 0.;

  inc_naln_hb = inc_naln_dccyk = 0;
  tot_inc_matrix_Mb_hb = tot_inc_matrix_Mb_dccyk = 0.;
  min_inc_matrix_Mb_hb = min_inc_matrix_Mb_dccyk = 0.;
  max_inc_matrix_Mb_hb = max_inc_matrix_Mb_dccyk = 0.;
  tot_inc_elapsed_secs_hb = tot_inc_elapsed_secs_dccyk = 0.;
  min_inc_elapsed_secs_hb = min_inc_elapsed_secs_dccyk = 0.;
  max_inc_elapsed_secs_hb = max_inc_elapsed_secs_dccyk = 0.;
 
  for(h = 0; h < th->N; h++) { 
    is_reported = (th->unsrt[h].flags & CM_HIT_IS_REPORTED) ? TRUE : FALSE;
    is_included = (th->unsrt[h].flags & CM_HIT_IS_INCLUDED) ? TRUE : FALSE;
    if(th->unsrt[h].ad != NULL) { 
      /* update reported stats */
      if(is_reported) { 
	if(th->unsrt[h].ad->used_hbands) { 
	  rep_naln_hb++;
	  tot_rep_matrix_Mb_hb    += th->unsrt[h].ad->matrix_Mb;
	  tot_rep_elapsed_secs_hb += th->unsrt[h].ad->elapsed_secs;
	  if(rep_naln_hb == 1) { 
	    min_rep_matrix_Mb_hb    = th->unsrt[h].ad->matrix_Mb;
	    max_rep_matrix_Mb_hb    = th->unsrt[h].ad->matrix_Mb;
	    min_rep_elapsed_secs_hb = th->unsrt[h].ad->elapsed_secs;
	    max_rep_elapsed_secs_hb = th->unsrt[h].ad->elapsed_secs;
	  }
	  else { 
	    min_rep_matrix_Mb_hb    = ESL_MIN(min_rep_matrix_Mb_hb,    th->unsrt[h].ad->matrix_Mb);
	    max_rep_matrix_Mb_hb    = ESL_MAX(max_rep_matrix_Mb_hb,    th->unsrt[h].ad->matrix_Mb);
	    min_rep_elapsed_secs_hb = ESL_MIN(min_rep_elapsed_secs_hb, th->unsrt[h].ad->elapsed_secs);
	    max_rep_elapsed_secs_hb = ESL_MAX(max_rep_elapsed_secs_hb, th->unsrt[h].ad->elapsed_secs);
	  }
	}
	else { 
	  rep_naln_dccyk++;
	  tot_rep_matrix_Mb_dccyk    += th->unsrt[h].ad->matrix_Mb;
	  tot_rep_elapsed_secs_dccyk += th->unsrt[h].ad->elapsed_secs;
	  if(rep_naln_dccyk == 1) { 
	    min_rep_matrix_Mb_dccyk    = th->unsrt[h].ad->matrix_Mb;
	    max_rep_matrix_Mb_dccyk    = th->unsrt[h].ad->matrix_Mb;
	    min_rep_elapsed_secs_dccyk = th->unsrt[h].ad->elapsed_secs;
	    max_rep_elapsed_secs_dccyk = th->unsrt[h].ad->elapsed_secs;
	  }
	  else { 
	    min_rep_matrix_Mb_dccyk    = ESL_MIN(min_rep_matrix_Mb_dccyk,    th->unsrt[h].ad->matrix_Mb);
	    max_rep_matrix_Mb_dccyk    = ESL_MAX(max_rep_matrix_Mb_dccyk,    th->unsrt[h].ad->matrix_Mb);
	    min_rep_elapsed_secs_dccyk = ESL_MIN(min_rep_elapsed_secs_dccyk, th->unsrt[h].ad->elapsed_secs);
	    max_rep_elapsed_secs_dccyk = ESL_MAX(max_rep_elapsed_secs_dccyk, th->unsrt[h].ad->elapsed_secs);
	  }
	}
      }
      /* update included stats */
      if(is_included) { 
	if(th->unsrt[h].ad->used_hbands) { 
	  inc_naln_hb++;
	  tot_inc_matrix_Mb_hb    += th->unsrt[h].ad->matrix_Mb;
	  tot_inc_elapsed_secs_hb += th->unsrt[h].ad->elapsed_secs;
	  if(inc_naln_hb == 1) { 
	    min_inc_matrix_Mb_hb    = th->unsrt[h].ad->matrix_Mb;
	    max_inc_matrix_Mb_hb    = th->unsrt[h].ad->matrix_Mb;
	    min_inc_elapsed_secs_hb = th->unsrt[h].ad->elapsed_secs;
	    max_inc_elapsed_secs_hb = th->unsrt[h].ad->elapsed_secs;
	  }
	  else { 
	    min_inc_matrix_Mb_hb    = ESL_MIN(min_inc_matrix_Mb_hb,    th->unsrt[h].ad->matrix_Mb);
	    max_inc_matrix_Mb_hb    = ESL_MAX(max_inc_matrix_Mb_hb,    th->unsrt[h].ad->matrix_Mb);
	    min_inc_elapsed_secs_hb = ESL_MIN(min_inc_elapsed_secs_hb, th->unsrt[h].ad->elapsed_secs);
	    max_inc_elapsed_secs_hb = ESL_MAX(max_inc_elapsed_secs_hb, th->unsrt[h].ad->elapsed_secs);
	  }
	}
	else { 
	  inc_naln_dccyk++;
	  tot_inc_matrix_Mb_dccyk    += th->unsrt[h].ad->matrix_Mb;
	  tot_inc_elapsed_secs_dccyk += th->unsrt[h].ad->elapsed_secs;
	  if(inc_naln_dccyk == 1) { 
	    min_inc_matrix_Mb_dccyk    = th->unsrt[h].ad->matrix_Mb;
	    max_inc_matrix_Mb_dccyk    = th->unsrt[h].ad->matrix_Mb;
	    min_inc_elapsed_secs_dccyk = th->unsrt[h].ad->elapsed_secs;
	    max_inc_elapsed_secs_dccyk = th->unsrt[h].ad->elapsed_secs;
	  }
	  else { 
	    min_inc_matrix_Mb_dccyk    = ESL_MIN(min_inc_matrix_Mb_dccyk,    th->unsrt[h].ad->matrix_Mb);
	    max_inc_matrix_Mb_dccyk    = ESL_MAX(max_inc_matrix_Mb_dccyk,    th->unsrt[h].ad->matrix_Mb);
	    min_inc_elapsed_secs_dccyk = ESL_MIN(min_inc_elapsed_secs_dccyk, th->unsrt[h].ad->elapsed_secs);
	    max_inc_elapsed_secs_dccyk = ESL_MAX(max_inc_elapsed_secs_dccyk, th->unsrt[h].ad->elapsed_secs);
	  }
	}
      }
    }
  }
  /* Output */
  if(rep_naln_hb > 0) { 
    avg_rep_matrix_Mb_hb    = tot_rep_matrix_Mb_hb    / (double) rep_naln_hb;
    avg_rep_elapsed_secs_hb = tot_rep_elapsed_secs_hb / (double) rep_naln_hb;
  }
  if(rep_naln_dccyk > 0) { 
    avg_rep_matrix_Mb_dccyk    = tot_rep_matrix_Mb_dccyk    / (double) rep_naln_dccyk;
    avg_rep_elapsed_secs_dccyk = tot_rep_elapsed_secs_dccyk / (double) rep_naln_dccyk;
  }
  if(inc_naln_hb > 0) { 
    avg_inc_matrix_Mb_hb    = tot_inc_matrix_Mb_hb    / (double) inc_naln_hb;
    avg_inc_elapsed_secs_hb = tot_inc_elapsed_secs_hb / (double) inc_naln_hb;
  }
  if(inc_naln_dccyk > 0) { 
    avg_inc_matrix_Mb_dccyk    = tot_inc_matrix_Mb_dccyk    / (double) inc_naln_dccyk;
    avg_inc_elapsed_secs_dccyk = tot_inc_elapsed_secs_dccyk / (double) inc_naln_dccyk;
  }

  fprintf(ofp, "Hit alignment statistics summary:\n");
  fprintf(ofp, "---------------------------------\n");
  //fprintf(ofp, "Query:       %s  [CLEN=%d]\n", cm->name, cm->clen);
  //if (cm->acc)  fprintf(ofp, "Accession:   %s\n", cm->acc);
  //if (cm->desc) fprintf(ofp, "Description: %s\n", cm->desc);
  //fprintf(ofp, "\n");
  if((rep_naln_hb + rep_naln_dccyk) > 0) { 
    fprintf(ofp, "%8s  %-22s  %9s  %25s  %34s\n", "", "", "", "    matrix size (Mb)     ", "      alignment time (secs)       ");
    fprintf(ofp, "%8s  %-22s  %9s  %25s  %34s\n", "", "", "", "-------------------------", "----------------------------------");
    fprintf(ofp, "%8s  %-22s  %9s  %7s  %7s  %7s  %7s  %7s  %7s  %7s\n", "category", "      algorithm       ", "# alns", "minimum", "average", "maximum", "minimum", "average", "maximum", "total");
    fprintf(ofp, "%8s  %-22s  %9s  %7s  %7s  %7s  %7s  %7s  %7s  %7s\n", "--------", "----------------------", "---------", "-------", "-------", "-------", "-------", "-------", "-------", "-------");
    /* reported */
    if(rep_naln_hb > 0) { 
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n", "reported", 
	      (used_cyk) ? "HMM banded CYK" : "HMM banded optimal acc", 
	      rep_naln_hb, min_rep_matrix_Mb_hb, avg_rep_matrix_Mb_hb, max_rep_matrix_Mb_hb, 
	      min_rep_elapsed_secs_hb, avg_rep_elapsed_secs_hb, max_rep_elapsed_secs_hb, tot_rep_elapsed_secs_hb);
    }
    else {
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7s  %7s  %7s  %7s  %7s  %7s  %7s\n", "reported", 
	      (used_cyk) ? "HMM banded CYK" : "HMM banded optimal acc", 
	      rep_naln_hb, "-", "-", "-", "-", "-", "-", "-");
    }
    if(rep_naln_dccyk > 0) { 
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n", "reported", "nonbanded D&C CYK", rep_naln_dccyk, 
	      min_rep_matrix_Mb_dccyk,    avg_rep_matrix_Mb_dccyk,    max_rep_matrix_Mb_dccyk, 
	      min_rep_elapsed_secs_dccyk, avg_rep_elapsed_secs_dccyk, max_rep_elapsed_secs_dccyk, tot_rep_elapsed_secs_dccyk);
    }
    else {
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7s  %7s  %7s  %7s  %7s  %7s  %7s\n", "reported", "nonbanded D&C CYK", rep_naln_dccyk, "-", "-", "-", "-", "-", "-", "-");
    }
    /* included */
    if(inc_naln_hb > 0) { 
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n", "included", 
	      (used_cyk) ? "HMM banded CYK" : "HMM banded optimal acc", 
	      inc_naln_hb, min_inc_matrix_Mb_hb,    avg_inc_matrix_Mb_hb,    max_inc_matrix_Mb_hb, 
	      min_inc_elapsed_secs_hb, avg_inc_elapsed_secs_hb, max_inc_elapsed_secs_hb, tot_inc_elapsed_secs_hb);
    }
    else {
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7s  %7s  %7s  %7s  %7s  %7s  %7s\n", "included", 
	      (used_cyk) ? "HMM banded CYK" : "HMM banded optimal acc", 
	      inc_naln_hb, "-", "-", "-", "-", "-", "-", "-");
    }
    if(inc_naln_dccyk > 0) { 
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n", "included", "nonbanded D&C CYK", inc_naln_dccyk, 
	      min_inc_matrix_Mb_dccyk,    avg_inc_matrix_Mb_dccyk,    max_inc_matrix_Mb_dccyk, 
	      min_inc_elapsed_secs_dccyk, avg_inc_elapsed_secs_dccyk, max_inc_elapsed_secs_dccyk, tot_inc_elapsed_secs_dccyk);
    }
    else {
      fprintf(ofp, "%8s  %-22s  %9" PRId64 "  %7s  %7s  %7s  %7s  %7s  %7s  %7s\n", "included", "nonbanded D&C CYK", inc_naln_dccyk, "-", "-", "-", "-", "-", "-", "-");
    }
  }
  else { 
    fprintf(ofp, "\n   [No hits detected that satisfy reporting thresholds]\n");
  }
  return eslOK;
}
/*---------------- end, standard output format ------------------*/

/*****************************************************************
 * 3. Tabular (parsable) output of pipeline results.
 *****************************************************************/

/* Function:  cm_tophits_TabularTargets()
 * Synopsis:  Output parsable table of per-sequence hits.
 * Incept:    EPN, Tue May 24 14:24:06 2011
 *            SRE, Wed Mar 18 15:26:17 2009 [Janelia] (p7_tophits.c)
 *
 * Purpose:   Output a parseable table of reportable per-sequence hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 */
int
cm_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, CM_TOPHITS *th, CM_PIPELINE *pli, int show_header)
{
  int status;
  int i;
  int tnamew = ESL_MAX(20, cm_tophits_GetMaxNameLength(th));
  int qnamew = ESL_MAX(20, strlen(qname));
  int qaccw  = ((qacc != NULL) ? ESL_MAX(9, strlen(qacc)) : 9);
  int taccw  = ESL_MAX(9, cm_tophits_GetMaxAccessionLength(th));
  int posw   = ESL_MAX(8, cm_tophits_GetMaxPositionLength(th));

  char *qnamestr = NULL;
  char *tnamestr = NULL;
  char *qaccstr  = NULL;
  char *taccstr  = NULL;
  char *posstr   = NULL;

  ESL_ALLOC(tnamestr, sizeof(char) * (tnamew));
  ESL_ALLOC(taccstr,  sizeof(char) * (taccw+1));
  ESL_ALLOC(qnamestr, sizeof(char) * (qnamew+1));
  ESL_ALLOC(qaccstr,  sizeof(char) * (qaccw+1));
  ESL_ALLOC(posstr,   sizeof(char) * (posw+1));

  for(i = 0; i < tnamew-1; i++) { tnamestr[i] = '-'; } tnamestr[tnamew-1] = '\0'; /* need to account for single '#' */
  for(i = 0; i < taccw;    i++) { taccstr[i]  = '-'; } taccstr[taccw]     = '\0';
  for(i = 0; i < qnamew;   i++) { qnamestr[i] = '-'; } qnamestr[qnamew]   = '\0';
  for(i = 0; i < qaccw;    i++) { qaccstr[i]  = '-'; } qaccstr[qaccw]     = '\0';
  for(i = 0; i < posw;     i++) { posstr[i]   = '-'; } posstr[posw]       = '\0';

  int h;

  if (show_header) { 
    fprintf(ofp, "#%-*s %-*s %-*s %-*s %7s %7s %*s %*s %6s %5s %4s %9s %6s %-s\n",
	    tnamew-1, "target name", taccw, "accession",  qnamew, "query name", qaccw, "accession", 
	    "cm from", "cm to", 
	    posw, "hit from", posw, "hit to", "strand", "trunc", "pass", "E-value", "score", "description of target");
    fprintf(ofp, "#%-*s %-*s %-*s %-*s %-7s %-7s %*s %*s %6s %5s %4s %9s %6s %s\n",
	    tnamew-1, tnamestr, taccw, taccstr, qnamew, qnamestr, qaccw, qaccstr, 
	    "-------", "-------", 
	    posw, posstr, posw, posstr, "------", "-----", "----", "---------", "-----", "---------------------");
  }
  for (h = 0; h < th->N; h++) { 
    if (th->hit[h]->flags & CM_HIT_IS_REPORTED)    {
      /* print occurs in three statements, b/c cfrom/cto can only be printed if we have a alignment display computed */
      fprintf(ofp, "%-*s %-*s %-*s %-*s ",
	      tnamew, th->hit[h]->name,
	      taccw,  ((th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0') ? th->hit[h]->acc : "-"),
	      qnamew, qname,
	      qaccw,  ((qacc != NULL && qacc[0] != '\0') ? qacc : "-"));

      if(th->hit[h]->ad == NULL) { fprintf(ofp, "%7s %7s ", "-", "-"); }
      else                       { fprintf(ofp, "%7d %7d ", th->hit[h]->ad->cfrom_emit, th->hit[h]->ad->cto_emit); }

      fprintf(ofp, "%*" PRId64 " %*" PRId64 " %6s %5s %4d %9.2g %6.1f %s\n",
	      posw, th->hit[h]->start,
	      posw, th->hit[h]->stop,
	      (th->hit[h]->in_rc == TRUE) ? "-" : "+",
	      cm_alidisplay_TruncString(th->hit[h]->ad), 
	      th->hit[h]->pass_idx, 
	      th->hit[h]->evalue,
	      th->hit[h]->score,
	      (th->hit[h]->desc != NULL) ? th->hit[h]->desc : "-");
    }
  }
  if(qnamestr != NULL) free(qnamestr);
  if(tnamestr != NULL) free(tnamestr);
  if(qaccstr  != NULL) free(qaccstr);
  if(qaccstr  != NULL) free(taccstr);
  if(posstr   != NULL) free(posstr);

  return eslOK;

 ERROR:
  if(qnamestr != NULL) free(qnamestr);
  if(tnamestr != NULL) free(tnamestr);
  if(qaccstr  != NULL) free(qaccstr);
  if(qaccstr  != NULL) free(taccstr);
  if(posstr   != NULL) free(posstr);

  return status;
}

/* Function:  cm_tophits_TabularTail()
 * Synopsis:  Print a trailer on a tabular output file.
 * Incept:    SRE, Tue Jan 11 16:13:30 2011 [Janelia]
 *
 * Purpose:   Print some metadata as a trailer on a tabular output file:
 *            date/time, the program, INFERNAL version info, the
 *            pipeline mode (SCAN or SEARCH), the query and target
 *            filenames, a spoof commandline recording the entire
 *            program configuration, and a "fini!" that's useful for
 *            detecting successful output completion.
 *
 * Args:      ofp       - open tabular output file (--tblout)
 *            progname  - "cmscan", for example
 *            pipemode  - CM_SEARCH_SEQS | CM_SCAN_MODELS
 *            qfile     - name of query file, or '-' for stdin, or '[none]' if NULL
 *            tfile     - name of target file, or '-' for stdin, or '[none]' if NULL
 *            go        - program configuration; used to generate spoofed command line
 *
 * Returns:   <eslOK>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslESYS> if time() or ctime_r() system calls fail.
 *                        
 * Xref:      J7/54
 */
int
cm_tophits_TabularTail(FILE *ofp, const char *progname, enum cm_pipemodes_e pipemode, const char *qfile, const char *tfile, const ESL_GETOPTS *go)
{
   time_t date           = time(NULL);
   char  *spoof_cmd      = NULL;
   char  *cwd            = NULL;
   char   timestamp[32];
   char   modestamp[16];
   int    status;

  if ((status = esl_opt_SpoofCmdline(go, &spoof_cmd)) != eslOK) goto ERROR;
  if (date == -1)                                               ESL_XEXCEPTION(eslESYS, "time() failed");
  if ((ctime_r(&date, timestamp)) == NULL)                      ESL_XEXCEPTION(eslESYS, "ctime_r() failed");
  switch (pipemode) {
  case CM_SEARCH_SEQS: strcpy(modestamp, "SEARCH"); break;
  case CM_SCAN_MODELS: strcpy(modestamp, "SCAN");   break;
  default:             ESL_EXCEPTION(eslEINCONCEIVABLE, "wait, what? no such pipemode");
  }
  esl_getcwd(&cwd);

  fprintf(ofp, "#\n");
  fprintf(ofp, "# Program:         %s\n",      (progname == NULL) ? "[none]" : progname);
  fprintf(ofp, "# Version:         %s (%s)\n", INFERNAL_VERSION, INFERNAL_DATE);
  fprintf(ofp, "# Pipeline mode:   %s\n",      modestamp);
  fprintf(ofp, "# Query file:      %s\n",      (qfile    == NULL) ? "[none]" : qfile);
  fprintf(ofp, "# Target file:     %s\n",      (tfile    == NULL) ? "[none]" : tfile);
  fprintf(ofp, "# Option settings: %s\n",      spoof_cmd);
  fprintf(ofp, "# Current dir:     %s\n",      (cwd      == NULL) ? "[unknown]" : cwd);
  fprintf(ofp, "# Date:            %s",        timestamp); /* timestamp ends in \n */
  fprintf(ofp, "# [ok]\n");

  free(spoof_cmd);
  if (cwd) free(cwd);
  return eslOK;

 ERROR:
  if (spoof_cmd) free(spoof_cmd);
  if (cwd)       free(cwd);
  return status;
}
/*------------------- end, tabular output -----------------------*/

/*****************************************************************
 * 4. Debugging/dev code
 *****************************************************************/

/* Function:  cm_tophits_Dump()
 * Synopsis:  Print contents of CM_TOPHITS for inspection.
 *
 * Purpose:   Print contents of the <CM_TOPHITS> <th> to
 *            stream <fp> for inspection. 
 *
 * Returns:   <eslOK>
 */
int
cm_tophits_Dump(FILE *fp, const CM_TOPHITS *th)
{
  uint64_t i;

  fprintf(fp, "CM_TOPHITS dump\n");
  fprintf(fp, "------------------\n");
  
  fprintf(fp, "N                             = %" PRId64 "\n", th->N);
  fprintf(fp, "Nalloc                        = %" PRId64 "\n", th->Nalloc);
  fprintf(fp, "nreported                     = %" PRId64 "\n", th->nreported);
  fprintf(fp, "nincluded                     = %" PRId64 "\n", th->nincluded);
  fprintf(fp, "is_sorted_by_score            = %s\n",  th->is_sorted_by_score            ? "TRUE" : "FALSE");
  fprintf(fp, "is_sorted_for_overlap_removal = %s\n",  th->is_sorted_for_overlap_removal ? "TRUE" : "FALSE");
  fprintf(fp, "is_sorted_by_position         = %s\n",  th->is_sorted_by_position         ? "TRUE" : "FALSE");
  if(th->is_sorted_by_score) { 
    for (i = 0; i < th->N; i++) {
      fprintf(fp, "SCORE SORTED HIT %" PRId64 ":\n", i);
      cm_hit_Dump(fp, th->hit[i]);
    }
  }
  else if(th->is_sorted_for_overlap_removal) { 
    for (i = 0; i < th->N; i++) {
      fprintf(fp, "OVERLAP REMOVAL SORTED HIT %" PRId64 ":\n", i);
      cm_hit_Dump(fp, th->hit[i]);
    }
  }
  else if(th->is_sorted_by_position) { 
    for (i = 0; i < th->N; i++) {
      fprintf(fp, "POSITION SORTED HIT %" PRId64 ":\n", i);
      cm_hit_Dump(fp, th->hit[i]);
    }
  }
  else { 
    for (i = 0; i < th->N; i++) {
      fprintf(fp, "UNSORTED HIT %" PRId64 ":\n", i);
      cm_hit_Dump(fp, &(th->unsrt[i]));
    }
  }
  return eslOK;
}

/* Function:  cm_hit_AllowTruncation()
 * Synopsis:  Determine if a hit returned from a truncated DP scanner
 *            is allowed, and return TRUE if it is and should be
 *            reported. Otherwise return FALSE. 
 *            
 *            The decision depends on the pipeline pass index, the
 *            mode of the hit, the locality mode of the CM and whether
 *            or not the hit contains the first and/or final residue
 *            of its source sequence.
 * 
 * Args:      cm       - the model, we need its emitmap and clen
 *            pass_idx - pass index hit was found in, if 
 *                       PLI_PASS_STD_ANY or PLI_PASS_5P_AND_3P_ANY
 *                       we allow all hits.
 *            start    - start position of the hit
 *            stop     - end position of the hit
 *            i0       - start position of source sequence
 *            j0       - end position of source sequence
 *            mode     - marginal mode of hit
 *            b        - entry state of the hit, we used a 
 *                       truncated begin into this state.
 * Returns:   informative string
 */
int 
cm_hit_AllowTruncation(CM_t *cm, int pass_idx, int64_t start, int64_t stop, int64_t i0, int64_t j0, char mode, int b)
{
  int in_local_mode = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;
  int nd = cm->ndidx[b];
  int lpos = (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd) ? cm->emap->lpos[nd] : cm->emap->lpos[nd] + 1;
  int rpos = (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATR_nd) ? cm->emap->rpos[nd] : cm->emap->rpos[nd] - 1;

  /* if our pass index allows 'any' hit, return TRUE */
  if(pass_idx == PLI_PASS_STD_ANY || pass_idx == PLI_PASS_5P_AND_3P_ANY) return TRUE;

  /* always allow full sequence hits i0..j0 */
  if(start == i0 && stop == j0) return TRUE; 

  /* if we get here, hit does not include full seq i0..j0 */
  if(in_local_mode) { 
    switch(mode) { 
    case TRMODE_J: /* local J hit that doesn't span full seq: always allow it */
      return TRUE; break;
    case TRMODE_L: /* local L hit that doesn't span full seq: j0 must be included */
      return (stop  == j0) ? TRUE : FALSE; break;
    case TRMODE_R: /* local R hit that doesn't span full seq: i0 must be included */
      return (start == i0) ? TRUE : FALSE; break;
    case TRMODE_T: /* local T hit that doesn't span full seq: don't allow it */
      return FALSE; break;
    default:
      return FALSE; break;
    }
  }
  else { /* local begins are off */
    switch(mode) { 
    case TRMODE_J: /* global J hit that doesn't span full seq: b must span 1..clen to allow it*/
      return (lpos == 1 && rpos == cm->clen) ? TRUE : FALSE; break;
    case TRMODE_L: /* global L hit that doesn't span full seq: j0 must be included and b must span 1 to allow it */
      return (stop  == j0   && lpos == 1)        ? TRUE : FALSE; break;
    case TRMODE_R: /* global R hit that doesn't span full seq: i0 must be included and b must span clen to allow it */
      return (start == i0   && rpos == cm->clen) ? TRUE : FALSE; break;
    case TRMODE_T: /* global T hit that doesn't span full seq: don't allow it */
      return FALSE; break;
    default:
      return FALSE; break;
    }
  }
  return FALSE; /* not reached */
}

/* Function:  cm_hit_Dump()
 * Synopsis:  Print contents of a CM_HIT for inspection.
 *
 * Purpose:   Print contents of the <CM_HIT> <h> to
 *            stream <fp> for inspection. 
 *
 * Returns:   <eslOK>
 */
int
cm_hit_Dump(FILE *fp, const CM_HIT *h)
{
  fprintf(fp, "CM_HIT dump\n");
  fprintf(fp, "------------------\n");
  fprintf(fp, "name      = %s\n",  h->name);
  fprintf(fp, "acc       = %s\n",  (h->acc  != NULL) ? h->acc  : "NULL");
  fprintf(fp, "desc      = %s\n",  (h->desc != NULL) ? h->desc : "NULL");
  fprintf(fp, "cm_idx    = %" PRId64 "\n", h->cm_idx);
  fprintf(fp, "seq_idx   = %" PRId64 "\n", h->seq_idx);
  fprintf(fp, "pass_idx  = %d\n",          h->pass_idx);
  fprintf(fp, "start     = %" PRId64 "\n", h->start);
  fprintf(fp, "stop      = %" PRId64 "\n", h->stop);
  fprintf(fp, "srcL      = %" PRId64 "\n", h->srcL);
  fprintf(fp, "maxW      = %d\n",  h->maxW);
  fprintf(fp, "in_rc     = %s\n",  h->in_rc ? "TRUE" : "FALSE");
  fprintf(fp, "root      = %d\n",  h->root);
  fprintf(fp, "mode      = %s\n",  MarginalMode(h->mode));
  fprintf(fp, "score     = %f\n",  h->score);
  fprintf(fp, "pvalue    = %f\n",  h->pvalue);
  fprintf(fp, "evalue    = %f\n",  h->evalue);
  if(h->flags == 0) { 
    fprintf(fp, "flags     = NONE\n");
  }
  else { 
    fprintf(fp, "flags:\n");
    if(h->flags & CM_HIT_IS_REPORTED)          fprintf(fp, "\tCM_HIT_IS_REPORTED\n");
    if(h->flags & CM_HIT_IS_INCLUDED)          fprintf(fp, "\tCM_HIT_IS_INCLUDED\n");
    if(h->flags & CM_HIT_IS_REMOVED_DUPLICATE) fprintf(fp, "\tCM_HIT_IS_REMOVED_DUPLICATE\n");
  }
  if(h->ad == NULL) { 
    fprintf(fp, "ad        = NULL\n");
  }
  else { 
    cm_alidisplay_Dump(fp, h->ad);
  }
  return eslOK;
}


/*****************************************************************
 * 4. Benchmark driver
 *****************************************************************/
#ifdef CM_TOPHITS_BENCHMARK
/* 
  gcc -o benchmark-cm-tophits -std=gnu99 -g -O2 -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DCM_TOPHITS_BENCHMARK cm_tophits.c -linfernal -lhmmer -leasel -lm 
  ./benchmark-cm-tophits

  As of 28 Dec 07, shows 0.20u for 10 lists of 10,000 hits each (at least ~100x normal expectation),
  so we expect top hits list time to be negligible for typical hmmsearch/hmmscan runs.
  
  If needed, we do have opportunity for optimization, however - especially in memory handling.
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stopwatch.h"
#include "esl_random.h"

#include "hmmer.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,    "181", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-M",        eslARG_INT,     "10", NULL, NULL,  NULL,  NULL, NULL, "number of top hits lists to simulate and merge",   0 },
  { "-N",        eslARG_INT,  "10000", NULL, NULL,  NULL,  NULL, NULL, "number of top hits to simulate per list",          0 },
  { "-X",        eslARG_INT,   "1000", NULL, NULL,  NULL,  NULL, NULL, "number of target sequences hits can come from",    0 },
  { "-Y",        eslARG_INT,     "10", NULL, NULL,  NULL,  NULL, NULL, "number of models hits can be to",                  0 },
  { "-Z",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "length of hits (fixed)",                           0 },
  { "-L",        eslARG_INT, "100000", NULL, NULL,  NULL,  NULL, NULL, "length of target sequences",                       0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump hit list after removing overlaps",            0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "benchmark driver for CM_TOPHITS";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = cm_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_STOPWATCH  *w        = esl_stopwatch_Create();
  ESL_RANDOMNESS *r        = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             N        = esl_opt_GetInteger(go, "-N");
  int             M        = esl_opt_GetInteger(go, "-M");
  int             X        = esl_opt_GetInteger(go, "-X");
  int             Y        = esl_opt_GetInteger(go, "-Y");
  int             L        = esl_opt_GetInteger(go, "-L");
  int             hitlen   = esl_opt_GetInteger(go, "-Z");
  CM_TOPHITS    **h        = NULL;
  CM_HIT         *hit      = NULL;
  float          *scores   = NULL;
  int            *seq_idxes= NULL;
  int            *cm_idxes = NULL;
  int            *starts   = NULL;
  int            *stops    = NULL;
  char            name[]   = "not_unique_name";
  char            acc[]    = "not_unique_acc";
  char            desc[]   = "Test description for the purposes of making the benchmark allocate space";
  int             i,j;
  int             nhits;
  int             status;
  char            errbuf[eslERRBUFSIZE];

  /* prep work: generate our sort keys before starting to time anything    */
  ESL_ALLOC(h,         sizeof(CM_TOPHITS *) * M); /* allocate pointers for M lists */
  ESL_ALLOC(scores,    sizeof(double) * N * M);   
  ESL_ALLOC(seq_idxes, sizeof(int) * N * M);   
  ESL_ALLOC(cm_idxes,  sizeof(int) * N * M);   
  ESL_ALLOC(starts,    sizeof(int) * N * M);   
  ESL_ALLOC(stops,     sizeof(int) * N * M);   
  for (i = 0; i < N*M; i++) scores[i]    = esl_random(r);
  for (i = 0; i < N*M; i++) seq_idxes[i] = esl_rnd_Roll(r, X);
  for (i = 0; i < N*M; i++) cm_idxes[i]  = esl_rnd_Roll(r, Y);
  for (i = 0; i < N*M; i++) starts[i]    = esl_rnd_Roll(r, L) + 1;
  for (i = 0; i < N*M; i++) stops[i]     = (esl_rnd_Roll(r, 2) == 0) ? ESL_MIN(starts[i] + hitlen - 1, L) : ESL_MAX(starts[i] - hitlen + 1, 1); 

  esl_stopwatch_Start(w);

  /* generate M "random" lists and sort them */
  for (j = 0; j < M; j++)
    {
      h[j] = cm_tophits_Create();
      for (i = 0; i < N; i++) 
	{ 
	  cm_tophits_CreateNextHit(h[j], &hit);
	  esl_strdup(name, -1, &(hit->name));
	  esl_strdup(acc,  -1, &(hit->acc));
	  esl_strdup(desc, -1, &(hit->desc));
	  hit->start   = starts[j*N+i];
	  hit->stop    = stops[j*N+i];
	  hit->score   = scores[j*N+i];
	  hit->seq_idx = seq_idxes[j*N+i];
	  hit->cm_idx  = cm_idxes[j*N+i];
	  hit->srcL    = L;
	  if     (hit->start < hit->stop) { hit->in_rc = FALSE; }
	  else if(hit->start > hit->stop) { hit->in_rc = TRUE;  }
	  else                            { hit->in_rc = (esl_rnd_Roll(r, 2) == 0) ? FALSE : TRUE; }
	}
      if(esl_opt_GetBoolean(go, "-v")) cm_tophits_Dump(stdout, h[j]);
    }

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time hit creation:                            ");
  esl_stopwatch_Start(w);

  /* then merge them into one big list in h[0] */
  for (j = 1; j < M; j++)
    {
      cm_tophits_Merge(h[0], h[j]);
      cm_tophits_Destroy(h[j]);
    }      
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time cm_tophits_Merge():                      ");
  esl_stopwatch_Start(w);
  
  cm_tophits_SortForOverlapRemoval(h[0]);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time cm_tophits_SortForOverlapRemoval():      ");
  esl_stopwatch_Start(w);

  if(esl_opt_GetBoolean(go, "-v")) cm_tophits_Dump(stdout, h[0]);

  if((status = cm_tophits_RemoveOverlaps(h[0], errbuf)) != eslOK) cm_Fail(errbuf);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time cm_tophits_RemoveOverlaps():             ");
  esl_stopwatch_Start(w);

  cm_tophits_SortByScore(h[0]);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time cm_tophits_SortByScore():                ");
  
  if(esl_opt_GetBoolean(go, "-v")) cm_tophits_Dump(stdout, h[0]);

  /* determine number of valid (not removed) hits */
  nhits = 0;
  for(i = 0; i < h[0]->N; i++) { 
    if(! (h[0]->hit[i]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) { 
      nhits++; 
    }
  }

  printf("# number of lists:               %d\n", M);
  printf("# sequence length                %d\n", L);
  printf("# number of sequences            %d\n", X);
  printf("# number of models               %d\n", Y);
  printf("# initial number of hits         %d\n", N*M);
  printf("# hit length                     %d\n", hitlen);
  printf("# number of non-overlapping hits %d\n", nhits);

  cm_tophits_Destroy(h[0]);
  status = eslOK;

 ERROR:
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  if (scores    != NULL) free(scores);
  if (seq_idxes != NULL) free(seq_idxes);
  if (cm_idxes  != NULL) free(cm_idxes);
  if (starts    != NULL) free(starts);
  if (stops     != NULL) free(stops);
  if (h         != NULL) free(h);
  return status;
}
#endif /*CM_TOPHITS_BENCHMARK*/
/*****************************************************************
 * 5. Test driver
 *****************************************************************/

#ifdef CM_TOPHITS_TESTDRIVE
/*
  gcc -o cm_tophits_utest -std=gnu99 -g -O2 -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel -DCM_TOPHITS_TESTDRIVE cm_tophits.c -linfernal -lhmmer -leasel -lm 
  ./tophits_test
*/
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stopwatch.h"
#include "esl_random.h"

#include "hmmer.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,    "181", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of top hits to simulate",                   0 },
  { "-X",        eslARG_INT,   "1000", NULL, NULL,  NULL,  NULL, NULL, "number of target sequences hits can come from",    0 },
  { "-Y",        eslARG_INT,     "10", NULL, NULL,  NULL,  NULL, NULL, "number of CMs hits can be to",                     0 },
  { "-L",        eslARG_INT, "100000", NULL, NULL,  NULL,  NULL, NULL, "length of target sequences",                       0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options]";
static char banner[] = "test driver for CM_TOPHITS";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = cm_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r        = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             N        = esl_opt_GetInteger(go, "-N");
  int             L        = esl_opt_GetInteger(go, "-L");
  int             X        = esl_opt_GetInteger(go, "-X");
  int             Y        = esl_opt_GetInteger(go, "-Y");
  CM_TOPHITS     *h1       = NULL;
  CM_TOPHITS     *h2       = NULL;
  CM_TOPHITS     *h3       = NULL;
  char            name[]   = "not_unique_name";
  char            acc[]    = "not_unique_acc";
  char            desc[]   = "Test description for the purposes of making the test driver allocate space";
  CM_HIT         *hit = NULL;
  int             i;
  int             status;
  char            errbuf[eslERRBUFSIZE];

  h1 = cm_tophits_Create();
  h2 = cm_tophits_Create();
  h3 = cm_tophits_Create();
  
  for (i = 0; i < N; i++) 
    {
      /* add hit to h1 */
      cm_tophits_CreateNextHit(h1, &hit);
      esl_strdup(name, -1, &(hit->name));
      esl_strdup(acc,  -1, &(hit->acc));
      esl_strdup(desc, -1, &(hit->desc));
      hit->start   = esl_rnd_Roll(r, L);
      hit->stop    = esl_rnd_Roll(r, L);
      hit->score   = esl_random(r);
      hit->seq_idx = esl_rnd_Roll(r, X);
      hit->cm_idx  = esl_rnd_Roll(r, Y);
      hit->srcL    = L;

      /* add hit to h2 */
      cm_tophits_CreateNextHit(h2, &hit);
      esl_strdup(name, -1, &(hit->name));
      esl_strdup(acc,  -1, &(hit->acc));
      esl_strdup(desc, -1, &(hit->desc));
      hit->start   = esl_rnd_Roll(r, L);
      hit->stop    = esl_rnd_Roll(r, L);
      hit->score   = 10.0 * esl_random(r);
      hit->seq_idx = esl_rnd_Roll(r, X);
      hit->cm_idx  = esl_rnd_Roll(r, Y);
      hit->srcL    = L;

      /* add hit to h3 */
      cm_tophits_CreateNextHit(h3, &hit);
      esl_strdup(name, -1, &(hit->name));
      esl_strdup(acc,  -1, &(hit->acc));
      esl_strdup(desc, -1, &(hit->desc));
      hit->start   = esl_rnd_Roll(r, L);
      hit->stop    = esl_rnd_Roll(r, L);
      hit->score   = 0.1 * esl_random(r);
      hit->seq_idx = esl_rnd_Roll(r, X);
      hit->cm_idx  = esl_rnd_Roll(r, Y);
      hit->srcL    = L;
    }

  /* add a few more hits */
  cm_tophits_CreateNextHit(h1, &hit);
  esl_strdup("third", -1, &(hit->name));
  hit->start   = L;
  hit->stop    = L;
  hit->score   = 20.0;
  hit->cm_idx  = 0;
  hit->seq_idx = 0;
  hit->srcL    = L+1;

  cm_tophits_CreateNextHit(h2, &hit);
  esl_strdup("second", -1, &(hit->name));
  hit->start   = L;
  hit->stop    = L;
  hit->score   = 30.0;
  hit->cm_idx  = 0;
  hit->seq_idx = 0;
  hit->srcL    = L+1;

  cm_tophits_CreateNextHit(h3, &hit);
  esl_strdup("first", -1, &(hit->name));
  hit->start   = L;
  hit->stop    = L;
  hit->score   = 40.0;
  hit->cm_idx  = 0;
  hit->seq_idx = 0;
  hit->srcL    = L+1;

  cm_tophits_CreateNextHit(h1, &hit);
  esl_strdup("thirdtolast", -1, &(hit->name));
  hit->start   = L+1;
  hit->stop    = L+1;
  hit->score   = -1.0;
  hit->cm_idx  = 0;
  hit->seq_idx = 0;
  hit->srcL    = L+1;

  cm_tophits_CreateNextHit(h2, &hit);
  esl_strdup("secondtolast", -1, &(hit->name));
  hit->start   = L+1;
  hit->stop    = L+1;
  hit->score   = -2.0;
  hit->cm_idx  = 0;
  hit->seq_idx = 0;
  hit->srcL    = L+1;

  cm_tophits_CreateNextHit(h3, &hit);
  esl_strdup("last", -1, &(hit->name));
  hit->start   = L+1;
  hit->stop    = L+1;
  hit->score   = -3.0;
  hit->cm_idx  = 0;
  hit->seq_idx = 0;
  hit->srcL    = L+1;
  
  cm_tophits_SortForOverlapRemoval(h1);
  if((status = cm_tophits_RemoveOverlaps(h1, errbuf)) != eslOK) cm_Fail(errbuf);
  cm_tophits_SortByScore(h1);
  if (strcmp(h1->hit[0]->name,   "third")        != 0)   esl_fatal("sort 1 failed (top is %s = %f)",  h1->hit[0]->name,   h1->hit[0]->score);
  if (strcmp(h1->hit[N+1]->name, "thirdtolast")  != 0)   esl_fatal("sort 1 failed (last is %s = %f)", h1->hit[N+1]->name, h1->hit[N+1]->score);

  cm_tophits_Merge(h1, h2);
  cm_tophits_SortForOverlapRemoval(h1);
  if((status = cm_tophits_RemoveOverlaps(h1, errbuf)) != eslOK) cm_Fail(errbuf);
  cm_tophits_SortByScore(h1);
  if (strcmp(h1->hit[0]->name,     "second")        != 0)   esl_fatal("sort 2 failed (top is %s = %f)",            h1->hit[0]->name,     h1->hit[0]->score);
  if (strcmp(h1->hit[1]->name,     "third")         != 0)   esl_fatal("sort 2 failed (second is %s = %f)",         h1->hit[1]->name,     h1->hit[1]->score);
  if (strcmp(h1->hit[2*N+2]->name, "thirdtolast")   != 0)   esl_fatal("sort 2 failed (second to last is %s = %f)", h1->hit[2*N+2]->name, h1->hit[2*N+2]->score);
  if (strcmp(h1->hit[2*N+3]->name, "secondtolast")  != 0)   esl_fatal("sort 2 failed (last is %s = %f)",           h1->hit[2*N+3]->name, h1->hit[2*N+3]->score);
  if (   h1->hit[0]->flags     & CM_HIT_IS_REMOVED_DUPLICATE)  esl_fatal("RemoveOverlaps failed 1");
  if (! (h1->hit[1]->flags     & CM_HIT_IS_REMOVED_DUPLICATE)) esl_fatal("RemoveOverlaps failed 2");
  if (   h1->hit[2*N+2]->flags & CM_HIT_IS_REMOVED_DUPLICATE)  esl_fatal("RemoveOverlaps failed 3");
  if (! (h1->hit[2*N+3]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) esl_fatal("RemoveOverlaps failed 4");

  cm_tophits_Merge(h1, h3);
  cm_tophits_SortForOverlapRemoval(h1);
  if((status = cm_tophits_RemoveOverlaps(h1, errbuf)) != eslOK) cm_Fail(errbuf);
  cm_tophits_SortByScore(h1);
  if (strcmp(h1->hit[0]->name,     "first")         != 0)   esl_fatal("sort 3 failed (top    is %s = %f)",         h1->hit[0]->name,     h1->hit[0]->score);
  if (strcmp(h1->hit[1]->name,     "second")        != 0)   esl_fatal("sort 3 failed (second is %s = %f)",         h1->hit[1]->name,     h1->hit[1]->score);
  if (strcmp(h1->hit[2]->name,     "third")         != 0)   esl_fatal("sort 3 failed (third  is %s = %f)",         h1->hit[2]->name,     h1->hit[2]->score);
  if (strcmp(h1->hit[3*N+3]->name, "thirdtolast")   != 0)   esl_fatal("sort 3 failed (third to last is %s = %f)",  h1->hit[3*N+3]->name, h1->hit[3*N+3]->score);
  if (strcmp(h1->hit[3*N+4]->name, "secondtolast")  != 0)   esl_fatal("sort 3 failed (second to last is %s = %f)", h1->hit[3*N+4]->name, h1->hit[3*N+4]->score);
  if (strcmp(h1->hit[3*N+5]->name, "last")          != 0)   esl_fatal("sort 3 failed (last is %s = %f)",           h1->hit[3*N+5]->name, h1->hit[3*N+5]->score);
  if (   h1->hit[0]->flags     & CM_HIT_IS_REMOVED_DUPLICATE)  esl_fatal("RemoveOverlaps failed 5");
  if (! (h1->hit[1]->flags     & CM_HIT_IS_REMOVED_DUPLICATE)) esl_fatal("RemoveOverlaps failed 6");
  if (! (h1->hit[2]->flags     & CM_HIT_IS_REMOVED_DUPLICATE)) esl_fatal("RemoveOverlaps failed 7");
  if (   h1->hit[3*N+3]->flags & CM_HIT_IS_REMOVED_DUPLICATE)  esl_fatal("RemoveOverlaps failed 8");
  if (! (h1->hit[3*N+4]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) esl_fatal("RemoveOverlaps failed 9");
  if (! (h1->hit[3*N+5]->flags & CM_HIT_IS_REMOVED_DUPLICATE)) esl_fatal("RemoveOverlaps failed 10");

  if (cm_tophits_GetMaxNameLength(h1) != strlen(name)) esl_fatal("GetMaxNameLength() failed");

  cm_tophits_Destroy(h1);
  cm_tophits_Destroy(h2);
  cm_tophits_Destroy(h3);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*CM_TOPHITS_TESTDRIVE*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/






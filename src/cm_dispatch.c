/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_dispatch.c
 * EPN, Wed Dec  6 06:11:46 2006
 * 
 * Dispatch functions for aligning and searching seqs
 * with a CM.
 * 
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "easel.h"
#include "esl_gumbel.h"
#include "esl_msa.h"         
#include "esl_stack.h"
#include "esl_stopwatch.h"   

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* 
 * Function: ActuallySearchTarget()
 * Incept:   EPN, Wed Nov 14 10:43:16 2007
 *
 * Purpose:  Given a CM and a sequence, call the correct search algorithm
 *           based on search_opts. 
 * 
 * Args:     cm              - the covariance model
 *           errbuf          - char buffer for reporting errors
 *           sround          - filtering round we're currently on, 
 *                             if sround == cm->fi->nrounds, we're done filtering (and possibly never filtered)
 *           dsq             - the target sequence (digitized)
 *           i0              - start of target subsequence (often 1, beginning of dsq)
 *           j0              - end of target subsequence (often L, end of dsq)
 *           sround          - filter round, if > 0, we're filtering
 *           results         - [0..cm-fi->nrounds] search_results_t to keep results for each round in, must be non NULL and empty
 *           ret_flen        - RETURN: subseq len that survived filter (NULL if not filtering)
 *           ret_sc          - RETURN: Highest scoring hit from search (even if below cutoff).
 *
 * Returns: eslOK on success.
 */
int ActuallySearchTarget(CM_t *cm, char *errbuf, int sround, ESL_DSQ *dsq, int i0, int j0, search_results_t **results, int *ret_flen, float *ret_sc)
{
  int               status;          /* easel status code */
  float             sc;              /* score of best hit in seq */
  float             bwd_sc;          /* score of best hit from Backward HMM algs */
  int               flen;            /* filter length, length of i0..j0 that survives filter */
  int               h;               /* counter over hits */
  int               i, j;            /* subseq start/end points */
  int               do_collapse;     /* TRUE to collapse multiple overlapping hits together into a single hit */
  int               next_j;          /* for collapsing overlapping hits together */
  int               min_i;           /* a start point, used if we're scanning with HMM */
  int               h_existing;      /* number of hits in round_results that exist when this function is entered */
  SearchInfo_t     *si = cm->si;     /* the SearchInfo */

  /* convenience pointers to cm->si for this 'filter round' of searching */
  float             cutoff;          /* cutoff for this round, HMM or CM, whichever is relevant for this round */
  int               stype;           /* search type for this round SEARCH_WITH_HMM, SEARCH_WITH_HYBRID, or SEARCH_WITH_CM */
  ScanMatrix_t     *smx;             /* scan matrix for this round, != NULL only if SEARCH_WITH_CM, and must == cm->smx if we're in the final round */
  HybridScanInfo_t *hsi;             /* hybrid scan info for this round, NULL unless stype is SEARCH_WITH_HYBRID */
  search_results_t *round_results;   /* search_results for this round */

  /* Contract checks */
  if(!(cm->flags & CMH_BITS))          ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), CMH_BITS flag down.\n");
  if(si == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(): search info cm->si is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(): dsq is NULL.");
  if(!(cm->flags & CMH_BITS))          ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(): CMH_BITS flag down.\n");
  if(sround > si->nrounds)             ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(): current round %d is greater than cm->si->nrounds: %d\n", sround, si->nrounds);
  if(results[sround] == NULL)          ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(): results for current round %d are NULL\n", sround);

  ESL_DPRINTF1(("In ActuallySearchTarget(), round: %d\n", sround));

  flen = (j0-i0+1);

  /* TEMPORARY */
  if(si->stype[sround] == SEARCH_WITH_HYBRID) cm_Fail("ActuallySearchTarget, hybrid filtering not yet implemented.\n");

  /* copy info for this round from SearchInfo fi */
  cm->search_opts = si->search_opts[sround]; 
  cutoff          = si->cutoff[sround];
  stype           = si->stype[sround];
  smx             = si->smx[sround]; /* may be NULL */
  hsi             = si->hsi[sround]; /* may be NULL */
  round_results   = results[sround]; /* may not be NULL, contract enforced this */
  h_existing      = round_results->num_results; /* remember this, b/c we only want rescan survivors found in *this* function call */

  /* SEARCH_WITH_HMM section */
  if(stype == SEARCH_WITH_HMM) { 
    /* some SEARCH_WITH_HMM specific contract checks */
    if(cm->cp9 == NULL)                    ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), trying to use CP9 HMM that is NULL.\n");
    if(!(cm->cp9->flags & CPLAN9_HASBITS)) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");
    if(smx != NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), round %d, SEARCH_WITH_HMM but smx != NULL.\n", sround);
    if(hsi != NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), round %d, SEARCH_WITH_HMM but hsi != NULL.\n", sround);
    if(! (cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD))
      ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), search type for this round is SEARCH_WITH_HMM, but CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags are both down.");

    search_results_t *fwd_results;
    /* Scan the (sub)seq in forward direction w/Viterbi or Forward, getting j end points of hits above cutoff */
    fwd_results = CreateResults(INIT_RESULTS);
    if(cm->search_opts & CM_SEARCH_HMMVITERBI) { 
      if((status = cp9_FastViterbi(cm, errbuf, cm->cp9_mx, dsq, i0, j0, cm->W, cutoff, fwd_results, 
				   TRUE,   /* we're scanning */
				   FALSE,  /* we're not ultimately aligning */
				   TRUE,   /* be memory efficient */
				   NULL, NULL, NULL,  /* don't return best score at each posn, best scoring posn, or traces */
				   &sc)) != eslOK) return status;
    }
    else if(cm->search_opts & CM_SEARCH_HMMFORWARD) { 
      if((status = cp9_FastForward(cm, errbuf, cm->cp9_mx, dsq, i0, j0, cm->W, cutoff, fwd_results,
				   TRUE,   /* we're scanning */
				   FALSE,  /* we're not ultimately aligning */
				   TRUE,   /* be memory efficient */
				   NULL, NULL, /* don't return best score at each posn, or best scoring posn */
				   &sc)) != eslOK) return status;
    }
    /* Remove overlapping hits, if we're being greedy */
    if(cm->search_opts & CM_SEARCH_HMMGREEDY) { /* resolve overlaps by being greedy */
      ESL_DASSERT1((i0 == 1)); /* EPN, Tue Nov 27 13:59:31 2007 not sure why this is here */
      remove_overlapping_hits (fwd_results, i0, j0);
    }

    /* determine start points (i) of the hits based on backward direction (Viterbi or Backward) scan starting at j */
    for(h = 0; h < fwd_results->num_results; h++) {
      min_i = (fwd_results->data[h].stop - cm->W + 1) >= 1 ? (fwd_results->data[h].stop - cm->W + 1) : 1;
      if(cm->search_opts & CM_SEARCH_HMMVITERBI) { 
	if((status = cp9_FastViterbiBackward(cm, errbuf, cm->cp9_mx, dsq, min_i, fwd_results->data[h].stop, cm->W, cutoff, 
					     round_results, /* report hits to this round's results */
					     TRUE,   /* we're scanning */
					     FALSE,  /* we're not ultimately aligning */
					     TRUE,   /* be memory efficient */
					     NULL, NULL, NULL,  /* don't return best score at each posn, best scoring posn, or traces */
					     &bwd_sc)) != eslOK) return status;
      }
      else { 
	if((status = Xcp9_FastBackward(cm, errbuf, cm->cp9_mx, dsq, min_i, fwd_results->data[h].stop, cm->W, cutoff, 
				       round_results, /* report hits to this round's results */
				       TRUE,   /* we're scanning */
				       FALSE,  /* we're not ultimately aligning */
				       TRUE,   /* be memory efficient */
				       NULL, NULL,   /* don't return best score at each posn, best scoring posn */
				       &bwd_sc)) != eslOK) return status;
      }
      /* this only works if we've saved the matrices, and didn't do scan mode for both Forward and Backward:
       * debug_check_CP9_FB(fmx, bmx, cm->cp9, cur_best_hmm_bsc, i0, j0, dsq); */
      if(bwd_sc > sc) sc = bwd_sc;
    }	  
    FreeResults(fwd_results);
  }
  /* end of SEARCH_WITH_HMM section */
  else if(stype == SEARCH_WITH_HYBRID) { 
    /* contract check */
    if(smx != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), round %d, SEARCH_WITH_HYBRID but smx != NULL.\n", sround);
    if(hsi == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(): current round %d is type SEARCH_WITH_HYBRID, but hsi is NULL\n", sround);
    if((status = cm_cp9_HybridScan(cm, errbuf, cm->cp9_mx, dsq, hsi, i0, j0, hsi->W, cutoff, round_results, 
				   NULL, NULL, /* don't return best score at each posn, and best scoring posn */
				   &sc)) != eslOK) return status;
  }  
  else { /* stype == SEARCH_WITH_CM */
    ESL_DASSERT1((stype == SEARCH_WITH_CM));
    if(smx == NULL)                             ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), round %d, SEARCH_WITH_CM but smx == NULL.\n", sround);
    if(sround == si->nrounds && smx != cm->smx) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(), final round %d, SEARCH_WITH_CM but smx != cm->smx.\n", sround);
    if(hsi != NULL)                             ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallySearchTarget(): round %d is type SEARCH_WITH_CM, but hsi is NULL\n", sround);

    if(cm->search_opts & CM_SEARCH_HBANDED) {
      if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, i0, j0, cm->cp9b, TRUE, 0)) != eslOK) return status; 
      if(cm->search_opts & CM_SEARCH_INSIDE) { if((status = FastFInsideScanHB(cm, errbuf, dsq, i0, j0, cutoff, round_results, cm->hbmx, &sc)) != eslOK) return status; }
      else                                   { if((status = FastCYKScanHB    (cm, errbuf, dsq, i0, j0, cutoff, round_results, cm->hbmx, &sc)) != eslOK) return status; }
    }
    else { /* don't do HMM banded search */
      if(cm->search_opts & CM_SEARCH_INSIDE) { if((status = FastIInsideScan(cm, errbuf, smx, dsq, i0, j0, cutoff, round_results, NULL, &sc)) != eslOK) return status; }
      else                                   { if((status = FastCYKScan    (cm, errbuf, smx, dsq, i0, j0, cutoff, round_results, NULL, &sc)) != eslOK) return status; }
    }    
  }

  if(sround < si->nrounds) { /* we're filtering */
    int   prev_j = j0;
    int   nhits  = round_results->num_results;
    
    /* To be safe, we only trust that i..j of our filter-passing hit is within the real hit,
     * so we add (W-1) to start point i and subtract (W-1) from j, and treat this region j-(W-1)..i+(W-1)
     * as having survived the filter.
     */
    for(h = h_existing; h < nhits; h++) {
      if(round_results->data[h].stop > prev_j) ESL_EXCEPTION(eslEINCOMPAT, "j's not in descending order");
      prev_j = round_results->data[h].stop;

      i = ((round_results->data[h].stop  - (cm->W-1)) >= 1)    ? (round_results->data[h].stop  - (cm->W-1)) : 1;
      j = ((round_results->data[h].start + (cm->W-1)) <= j0)   ? (round_results->data[h].start + (cm->W-1)) : j0;

      if((h+1) < nhits) next_j = ((round_results->data[h+1].start + (cm->W-1)) <= j0) ? (round_results->data[h+1].start + (cm->W-1)) : j0;
      else              next_j = -1;
      
      /* Collapse multiple overlapping hits together into a single hit. 
       * *Unless* our next round of searching is the final one, and we're going to do HMM banded search,
       * in which case we want to treat each hit separately, so we get more reasonable bands.
       */
      do_collapse = (((sround+1) == si->nrounds) && (si->search_opts[si->nrounds] & CM_SEARCH_HBANDED)) ? FALSE : TRUE;
      if(do_collapse) { 
	while(((h+1) < nhits) && (next_j >= i)) { /* suck in hit */
	  h++;
	  i = ((round_results->data[h].stop - (cm->W-1)) >= 1) ? (round_results->data[h].stop - (cm->W-1)) : 1;
	  if((h+1) < nhits) next_j = ((round_results->data[h+1].start + (cm->W-1)) <= j0) ? (round_results->data[h+1].start + (cm->W-1)) : j0;
	  else              next_j = -1;
	  printf("\tsucked in subseq: hit %d new_i: %d j (still): %d\n", h, i, j);
	}
      }
      /* next round: research this chunk that survived the filter */
      if((status = ActuallySearchTarget(cm, errbuf, (sround+1), dsq, i, j, results, NULL, NULL)) != eslOK) return status;
    }
  }
  else { /* we're done filtering, and we're done searching, get alignments if nec */
    if((round_results->num_results > 0) && (! (cm->search_opts & CM_SEARCH_NOALIGN))) {
      if((status = ActuallyAlignTargets(cm, errbuf, NULL, 
					dsq, round_results, h_existing,     /* put function into dsq_mode, designed for aligning search hits */
					0, 0, 0, NULL)) != eslOK) return status;
    }
  }
  if(ret_sc != NULL) *ret_sc = sc;
  return eslOK;
}  


/* 
 * Function: ActuallyAlignTargets
 * Incept:   EPN, Thu Nov 15 11:35:23 2007
 *
 * Purpose:  Given a CM and sequences, do preliminaries, call the correct 
 *           alignment function and return parsetrees and optionally postal codes 
 *           (if cm->align_opts & CM_ALIGN_POST).
 *
 *           Two different modes are possible dependent on input args. Mode
 *           is checked for during contract enforcement.
 *
 *           sq_mode: seqs_to_aln != NULL; dsq == NULL; results == NULL.
 *                    align the seqs_to_aln->nseq ESL_SQ sq sequences store
 *                    parsetrees or CP9 traces and/or postal codes in
 *                    seqs_to_aln.
 *
 *          dsq_mode: seqs_to_aln == NULL; dsq != NULL, results != NULL.
 *                    align the search results (hits) in results, which
 *                    are all subsequences of a single sequence (dsq).
 *                    parstrees are stored in seacrh_results.
 *
 * Args:     CM             - the covariance model
 *           errbuf         - char buffer for reporting errors
 *           seqs_to_aln    - the sequences (if sq_mode)
 *           dsq            - a single digitized sequence (if dsq_mode)
 *           search_results - search results with subsequence indices of dsq to align (if dsq_mode)
 *           first_result   - index of first result in search_results to align (if dsq_mode)
 *           bdump_level    - verbosity level for band related print statements
 *           debug_level    - verbosity level for debugging print statements
 *           silent_mode    - TRUE to not print anything, FALSE to print scores 
 *           r              - source of randomness (NULL unless CM_ALIGN_SAMPLE)
 * 
 * Returns:  eslOK on success;
 *           eslERANGE if required memory for a DP matrix is too big;
 *           eslEINCOMPAT if input parameters violate contract;
 *           eslEMEM on memory allocation error;
 *           eslFAIL if some other error;
 *           if(!eslOK) errbuf is filled with informative error message
 */
int
ActuallyAlignTargets(CM_t *cm, char *errbuf, seqs_to_aln_t *seqs_to_aln, ESL_DSQ *dsq, search_results_t *search_results,
		     int first_result, int bdump_level, int debug_level, int silent_mode, ESL_RANDOMNESS *r)
{
  int status;
  ESL_STOPWATCH *watch;         /* for timings */
  int sq_mode  = FALSE;         /* we're aligning nseq seqs in sq */
  int dsq_mode = FALSE;         /* we're aligning search_results->num_results seqs, all subseqs of dsq */
  int nalign   = 0;             /* number of sequences we're aligning */
  ESL_DSQ *cur_dsq;             /* ptr to digitized sequence we're currently aligning */
  Parsetree_t **cur_tr;         /* pointer to the pointer to the parsetree we're currently creating */
  int L;                        /* length of sequence/subseq we're currently aligning */
  int i;                        /* counter over sequences */
  int ip;                       /* offset index in search_results */
  int v;                        /* state counter */
  char        **postcode1 = NULL;/* posterior decode array of strings, tens place ('9' for 93) */
  char        **postcode2 = NULL;/* posterior decode array of strings, ones place ('3' for 93) */
  Parsetree_t **tr       = NULL;/* parse trees for the sequences */
  CP9trace_t  **cp9_tr   = NULL;/* CP9 traces for the sequences */
  float         sc;		/* score for one sequence alignment */
  float         maxsc;	        /* max score in all seqs */
  float         minsc;	        /* min score in all seqs */
  float         avgsc;      	/* avg score over all seqs */
  float         tmpsc;          /* temporary score */

  /* variables related to CM Plan 9 HMMs */
  CP9_t       *hmm;             /* constructed CP9 HMM */
  CP9Bands_t  *cp9b;            /* data structure for hmm bands (bands on the hmm states) 
				 * and arrays for CM state bands, derived from HMM bands */
  CP9Map_t       *cp9map;       /* maps the hmm to the cm and vice versa */
  float           swentry;	/* S/W aggregate entry probability       */
  float           swexit;       /* S/W aggregate exit probability        */

  /* variables related to the do_sub option */
  int                spos;         /* HMM node most likely to have emitted posn 1 of target seq */
  int                spos_state;   /* HMM state type for curr spos 0=match or 1=insert */
  int                epos;         /* HMM node most likely to have emitted posn L of target seq */
  int                epos_state;   /* HMM state type for curr epos 0=match or  1=insert */

  CMSubMap_t        *submap;
  CM_t              *sub_cm;       /* sub covariance model                      */
  CP9_t             *sub_hmm;      /* constructed CP9 HMM */
  CP9Map_t          *sub_cp9map;   /* maps the sub_hmm to the sub_cm and vice versa */
  CP9Bands_t        *sub_cp9b;     /* data structure for hmm bands (bands on the hmm states) 
				    * and arrays for CM state bands, derived from HMM bands */

  CM_t              *orig_cm;      /* the original, template covariance model the sub CM was built from */
  CP9_t             *orig_hmm;     /* original CP9 HMM built from orig_cm */
  CP9Map_t          *orig_cp9map;  /* original CP9 map */
  CP9Bands_t        *orig_cp9b;    /* original CP9Bands */
  Parsetree_t       *orig_tr;      /* parsetree for the orig_cm; created from the sub_cm parsetree */

  /* variables related to query dependent banding (qdb) */
  int    expand_flag;           /* TRUE if the dmin and dmax vectors have just been 
				 * expanded (in which case we want to recalculate them 
				 * before we align a new sequence), and FALSE if not*/
  int *orig_dmin;               /* original dmin values passed in */
  int *orig_dmax;               /* original dmax values passed in */

  /* variables related to inside/outside */
  float           ***alpha = NULL;    /* alpha DP matrix for non-banded Inside() */
  float           ***beta  = NULL;    /* beta DP matrix for non-baned Outside() */
  CM_HB_MX           *out_mx;         /* outside matrix for HMM banded Outside() */

  float             *parsesc; /* parsetree scores of each sequence */

  /* declare and initialize options */
  int do_small     = FALSE;   /* TRUE to use D&C small alignment algs */
  int do_local     = FALSE;   /* TRUE to do local alignment */
  int do_qdb       = FALSE;   /* TRUE to do QDB alignment */
  int do_hbanded   = FALSE;   /* TRUE to do HMM banded alignment */
  int use_sums     = FALSE;   /* TRUE to use posterior sums for HMM banded alignment */
  int do_sub       = FALSE;   /* TRUE to align to a sub CM */
  int do_hmmonly   = FALSE;   /* TRUE to align with an HMM only */
  int do_scoreonly = FALSE;   /* TRUE to only calculate the score */
  int do_inside    = FALSE;   /* TRUE to do Inside also */
  int do_post      = FALSE;   /* TRUE to calculate posterior probabilities */
  int do_timings   = FALSE;   /* TRUE to report timings */
  int do_check     = FALSE;   /* TRUE to check posteriors from Inside/Outside */
  int do_sample    = FALSE;   /* TRUE to sample from an Inside matrix */
  int do_optacc    = FALSE;   /* TRUE to find optimally accurate alignment instead of CYK */
  int do_hmmsafe   = FALSE;   /* TRUE to realign seqs with HMM banded parses < 0. bits (only works if !do_optacc && !do_post && do_hbanded)*/

  /* Contract check */
  if(!(cm->flags & CMH_BITS))                            ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), CMH_BITS flag down.\n");
  if(r == NULL && (cm->align_opts & CM_ALIGN_SAMPLE))    ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), no source of randomness, but CM_ALIGN_SAMPLE alignment option on.\n");
  if(r != NULL && (!(cm->align_opts & CM_ALIGN_SAMPLE))) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), we have a source of randomness, but CM_ALIGN_SAMPLE alignment option off.\n");
  if((cm->align_opts & CM_ALIGN_POST)      && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), CM_ALIGN_POST and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), CM_ALIGN_SCOREONLY and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_POST))       ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), CM_ALIGN_SCOREONLY and CM_ALIGN_POST options are incompatible.\n");

  /* determine mode */
  if     (seqs_to_aln != NULL && (dsq == NULL && search_results == NULL))  sq_mode = TRUE;
  else if(seqs_to_aln == NULL && (dsq != NULL && search_results != NULL)) dsq_mode = TRUE;
  else   ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), can't determine mode (sq_mode or dsq_mode).\n");

  if( sq_mode && (seqs_to_aln->sq        == NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in sq_mode, seqs_to_aln->sq is NULL.\n");
  if( sq_mode && (seqs_to_aln->tr        != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in sq_mode, seqs_to_aln->tr is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->cp9_tr    != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in sq_mode, seqs_to_aln->cp9_tr is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->postcode1 != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in sq_mode, seqs_to_aln->postcode1 is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->postcode2 != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in sq_mode, seqs_to_aln->postcode2 is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->sc        != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in sq_mode, seqs_to_aln->sc is non-NULL.\n");
  
  if(dsq_mode && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in dsq_mode, CM_ALIGN_HMMVITERBI option on.\n");
  if(dsq_mode && (cm->align_opts & CM_ALIGN_INSIDE))     ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in dsq_mode, CM_ALIGN_INSIDE option on.\n");
  if(dsq_mode && (cm->align_opts & CM_ALIGN_SAMPLE))     ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in dsq_mode, CM_ALIGN_SAMPLE option on.\n");
  if(dsq_mode && search_results == NULL)                 ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in dsq_mode, search_results are NULL.\n");
  if(dsq_mode && (first_result > search_results->num_results)) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), in dsq_mode, first_result: %d > search_results->num_results: %d\n", first_result, search_results->num_results);

  /* save a copy of the align_opts we entered function with, we may change some of these for
   * individual target sequences, and we want to be able to change them back
   */
  if(cm->align_opts  & CM_ALIGN_SMALL)      do_small     = TRUE;
  if(cm->config_opts & CM_CONFIG_LOCAL)     do_local     = TRUE;
  if(cm->align_opts  & CM_ALIGN_QDB)        do_qdb       = TRUE;
  if(cm->align_opts  & CM_ALIGN_HBANDED)    do_hbanded   = TRUE;
  if(cm->align_opts  & CM_ALIGN_SUMS)       use_sums     = TRUE;
  if(cm->align_opts  & CM_ALIGN_SUB)        do_sub       = TRUE;
  if(cm->align_opts  & CM_ALIGN_HMMVITERBI) do_hmmonly   = TRUE;
  if(cm->align_opts  & CM_ALIGN_INSIDE)     do_inside    = TRUE;
  if(cm->align_opts  & CM_ALIGN_POST)       do_post      = TRUE;
  if(cm->align_opts  & CM_ALIGN_TIME)       do_timings   = TRUE;
  if(cm->align_opts  & CM_ALIGN_CHECKINOUT) do_check     = TRUE;
  if(cm->align_opts  & CM_ALIGN_SCOREONLY)  do_scoreonly = TRUE;
  if(cm->align_opts  & CM_ALIGN_SAMPLE)     do_sample    = TRUE;
  if(cm->align_opts  & CM_ALIGN_OPTACC)     do_optacc    = TRUE;
  if(cm->align_opts  & CM_ALIGN_HMMSAFE)    do_hmmsafe   = TRUE;

  /* another contract check */
  if((do_sample + do_inside + do_post + do_hmmonly + do_scoreonly) > 1) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), exactly 0 or 1 of the following must be TRUE (== 1):\n\tdo_sample = %d\n\tdo_inside = %d\n\t do_post = %d\n\tdo_hmmonly = %d\n\tdo_scoreonly = %d\n", do_sample, do_inside, do_post, do_hmmonly, do_scoreonly);

  if(debug_level > 0) {
    printf("do_local    : %d\n", do_local);
    printf("do_qdb      : %d\n", do_qdb);
    printf("do_hbanded  : %d\n", do_hbanded);
    printf("use_sums    : %d\n", use_sums);
    printf("do_sub      : %d\n", do_sub);
    printf("do_hmmonly  : %d\n", do_hmmonly);
    printf("do_inside   : %d\n", do_inside);
    printf("do_small    : %d\n", do_small);
    printf("do_post     : %d\n", do_post);
    printf("do_timings  : %d\n", do_timings);
    printf("do_check    : %d\n", do_check);
    printf("do_scoreonly: %d\n", do_scoreonly);
    printf("do_sample   : %d\n", do_sample);
    printf("do_optacc   : %d\n", do_optacc);
    printf("do_hmmsafe  : %d\n", do_hmmsafe);
  }

  /* allocate out_mx, if needed, only if !do_sub, if do_sub each sub CM will need to allocate a new out_mx */
  out_mx = NULL;
  if((!do_sub) && (do_hbanded && (do_optacc || (do_post)))) out_mx = cm_hb_mx_Create(cm->M);

  if      (sq_mode)   nalign = seqs_to_aln->nseq;
  else if(dsq_mode) { nalign = search_results->num_results - first_result; silent_mode = TRUE; }

  /* If sqmode: potentially allocate tr, cp9_tr, and postcodes. We'll set
   * seqs_to_aln->tr, seqs_to_aln->cp9_tr, seqs_to_aln->postcode1, 
   * and seqs_to_aln->postcode2 to these guys at end of function.
   * 
   * If dsqmode: do not allocate parsetree pointers, they already exist 
   * in search_results.
   */
  tr       = NULL;
  cp9_tr   = NULL;
  postcode1= NULL;
  postcode2= NULL;
  if(sq_mode) {
    if(!do_hmmonly && !do_scoreonly && !do_inside)
      ESL_ALLOC(tr, sizeof(Parsetree_t *) * nalign);
    else if(do_hmmonly) /* do_hmmonly */
      ESL_ALLOC(cp9_tr, sizeof(CP9trace_t *) * nalign);
  }   
  ESL_ALLOC(parsesc, sizeof(float) * nalign);
  if(do_post) {
    ESL_ALLOC(postcode1, sizeof(char **) * nalign);
    ESL_ALLOC(postcode2, sizeof(char **) * nalign);
  }
  minsc =  FLT_MAX;
  maxsc = -FLT_MAX;
  avgsc = 0;
  watch = esl_stopwatch_Create();

  if(do_hbanded || do_sub) { /* We need a CP9 HMM to build sub_cms */
    if(cm->cp9 == NULL)                    ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets, trying to use CP9 HMM that is NULL\n");
    if(cm->cp9b == NULL)                   ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets, cm->cp9b is NULL\n");
    if(!(cm->cp9->flags & CPLAN9_HASBITS)) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets, trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");
    
    /* Keep data for the original CM safe; we'll be doing
     * pointer swapping to ease the sub_cm alignment implementation. */
    hmm         = cm->cp9;
    cp9b        = cm->cp9b;
    cp9map      = cm->cp9map;
    orig_hmm    = hmm;
    orig_cp9b   = cp9b;
    orig_cp9map = cp9map;
  }
  /* Copy the QD bands in case we expand them. */
  if(do_qdb) {
    if(bdump_level > 1) debug_print_bands(stdout, cm, cm->dmin, cm->dmax);
    expand_flag = FALSE;
    /* Copy dmin and dmax, so we can replace them after expansion */
    ESL_ALLOC(orig_dmin, sizeof(int) * cm->M);
    ESL_ALLOC(orig_dmax, sizeof(int) * cm->M);
    for(v = 0; v < cm->M; v++) {
      orig_dmin[v] = cm->dmin[v];
      orig_dmax[v] = cm->dmax[v];
    }	  
  }
  if(do_sub) { /* to get spos and epos for the sub_cm, 
	        * we config the HMM to local mode with equiprobable start/end points.*/
    swentry = ((hmm->M)-1.)/hmm->M; /* all start pts equiprobable, including 1 */
    swexit  = ((hmm->M)-1.)/hmm->M; /* all end   pts equiprobable, including M */
    CPlan9SWConfig(hmm, swentry, swexit);
    CP9Logoddsify(hmm);
  }
  orig_cm = cm;
  
  /*****************************************************************
   *  Collect parse trees for each sequence
   *****************************************************************/
  for (i = 0; i < nalign; i++) {
    if(do_timings) esl_stopwatch_Start(watch);
    if (sq_mode) { 
      cur_dsq = seqs_to_aln->sq[i]->dsq;
      cur_tr  = &(tr[i]);
      L       = seqs_to_aln->sq[i]->n;
    }
    else if (dsq_mode) {
      ip      = i + first_result;
      cur_dsq = dsq + search_results->data[ip].start - 1;
      cur_tr  = &(search_results->data[ip].tr);
      L       = search_results->data[ip].stop - search_results->data[ip].start + 1;
      ESL_DASSERT1((L >= 0));
    }
    if (L == 0) continue; /* silently skip zero length seqs */

    /* Special case, if do_hmmonly, align seq with Viterbi, print score and move on to next seq */
    if(sq_mode && do_hmmonly) {
      if(sq_mode && !silent_mode) printf("Aligning (to a CP9 HMM w/viterbi) %-20s", seqs_to_aln->sq[i]->name);
      if((status = cp9_FastViterbi(cm, errbuf, cm->cp9_mx, cur_dsq, 1, L, L, 0., NULL,
				   FALSE,  /* we are not scanning */
				   TRUE,   /* we are aligning */
				   FALSE,  /* don't be memory efficient */
				   NULL, NULL, /* don't return best sc at each posn, or best scoring posn */
				   &(cp9_tr[i]), /* return the trace */
				   &sc)) != eslOK) return status;
      if(sq_mode && !silent_mode) printf(" score: %10.2f bits\n", sc);
      parsesc[i] = sc;
      continue;
    }
    /* Special case, if do_scoreonly, align seq with full CYK inside, just to 
     * get the score. For testing, probably in cmscore. */
    if(sq_mode && do_scoreonly) {
      if(sq_mode && !silent_mode) printf("Aligning (w/full CYK score only) %-30s", seqs_to_aln->sq[i]->name);
      sc = CYKInsideScore(cm, cur_dsq, L, 0, 1, L, NULL, NULL); /* don't do QDB mode */
      if(sq_mode && !silent_mode) printf("    score: %10.2f bits\n", sc);
      parsesc[i] = sc;
      continue;
    }

    /* Potentially, do HMM calculations. */
    if((!do_sub) && do_hbanded) {
      if((status = cp9_Seq2Bands(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, cur_dsq, 1, L, orig_cp9b, FALSE, debug_level)) != eslOK) return status;
    }
    else if(do_sub) { 
      /* If we're in sub mode:
       * (1) Get HMM posteriors 
       * (2) Infer the start (spos) and end (epos) HMM states by 
       *     looking at the posterior matrix.
       * (3) Build the sub_cm from the original CM.
       *
       * If we're also doing HMM banded alignment to sub CM:
       * (4) Build a new CP9 HMM from the sub CM.
       * (5) Do Forward/Backward again, and get HMM bands 
       */
      
      /* (1) Get HMM posteriors */
      if((status = cp9_Seq2Posteriors(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, cur_dsq, 1, L, debug_level)) != eslOK) return status; 
      
      /* (2) infer the start and end HMM nodes (consensus cols) from posterior matrix.
       * Remember: we're necessarily in CP9 local mode, the --sub option turns local mode on. 
       */
      CP9NodeForPosn(orig_hmm, 1, L, 1, orig_cm->cp9_bmx, &spos, &spos_state, 0., TRUE, debug_level);
      CP9NodeForPosn(orig_hmm, 1, L, L, orig_cm->cp9_bmx, &epos, &epos_state, 0., FALSE, debug_level);
      /* Deal with special cases for sub-CM alignment:
       * If the most likely state to have emitted the first or last residue
       * is the insert state in node 0, it only makes sense to start modelling
       * at consensus column 1. */
      if(spos == 0 && spos_state == 1) spos = 1;
      if(epos == 0 && epos_state == 1) epos = 1;
      /* If most-likely HMM node to emit final position comes BEFORE most-likely HMM node to emit first position,
       * our HMM alignment is crap, default to using the full CM. */
      if(epos < spos) { spos = 1; epos = cm->cp9->M; } 
	  
      /* (3) Build the sub_cm from the original CM. */
      if(!(build_sub_cm(orig_cm, &sub_cm, 
			spos, epos,         /* first and last col of structure kept in the sub_cm  */
			&submap,            /* maps from the sub_cm to cm and vice versa           */
			debug_level)))      /* print or don't print debugging info                 */
	ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), unexpected error building a sub CM for seq %d.", i);
      /* Configure the sub_cm, the same as the cm, this will build a CP9 HMM if (do_hbanded), this will also:  */
      /* (4) Build a new CP9 HMM from the sub CM. */
      ConfigCM(sub_cm, NULL, NULL);
      cm    = sub_cm; /* orig_cm still points to the original CM */
      if(do_hbanded) { /* we're doing HMM banded alignment to the sub_cm */
	/* Get the HMM bands for the sub_cm */
	sub_hmm    = sub_cm->cp9;
	sub_cp9b   = sub_cm->cp9b;
	sub_cp9map = sub_cm->cp9map;
	/* (5) Do Forward/Backward again, and get HMM bands */
	if((status = cp9_Seq2Bands(sub_cm, errbuf, sub_cm->cp9_mx, sub_cm->cp9_bmx, sub_cm->cp9_bmx, cur_dsq, 1, L, sub_cp9b, FALSE, debug_level)) != eslOK) return status;
	hmm           = sub_hmm;    
	cp9b          = sub_cp9b;
	cp9map        = sub_cp9map;

	/* Create the out_mx if needed, cm == sub_cm */
	if(do_optacc || do_post) out_mx = cm_hb_mx_Create(cm->M);
      }
    }

    if(do_qdb) {
      /*Check if we need to reset the query dependent bands b/c they're currently expanded. */
      if(expand_flag) {
	for(v = 0; v < cm->M; v++) { cm->dmin[v] = orig_dmin[v]; cm->dmax[v] = orig_dmax[v]; }
	expand_flag = FALSE;
      }
      if((L < cm->dmin[0]) || (L > cm->dmax[0])) { 
	/* the seq we're aligning is outside the root band, so we expand.*/
	ExpandBands(cm, L, cm->dmin, cm->dmax);
	if(sq_mode && debug_level > 0) printf("Expanded bands for seq : %s\n", seqs_to_aln->sq[i]->name);
	if(bdump_level > 2) { printf("printing expanded bands :\n"); debug_print_bands(stdout, cm, cm->dmin, cm->dmax); }
	expand_flag = TRUE;
      }
    }

    if(sq_mode && !silent_mode) { 
      if(do_sub) printf("Aligning (to a sub CM) %-20s", seqs_to_aln->sq[i]->name);
      else       printf("Aligning %-30s", seqs_to_aln->sq[i]->name);
    }

    /* beginning of large if() else if() else if() ... statement */
    if(do_sample) { 
      if(do_hbanded) { /* sampling from inside HMM banded matrix */
	if((status = FastInsideAlignHB(cm, errbuf, cur_dsq, 1, L, cm->hbmx, NULL)) != eslOK) return status; /* errbuf will have been filled by FastInsideAlignHB() */
	if((status = SampleFromInsideHB(r, cm, errbuf, cur_dsq, L, cm->hbmx, cur_tr, &sc)) != eslOK) return status; /* errbuf will have been filled by SampleFromInsideHB() */
      }
      else { /* sampling from inside matrix, but not HMM banded */
	if((status = FastInsideAlign(cm, errbuf, cur_dsq, 1, L, &alpha, NULL)) != eslOK) return status; /* errbuf will have been filled by FastInsideAlign() */
	if((status = SampleFromInside(r, cm, errbuf, cur_dsq, L,  alpha, cur_tr, &sc)) != eslOK) return status; /* errbuf will have been filled by SampleFromInside() */
      }
    }
    else if(do_inside) { 
      if(do_hbanded) { /* HMM banded inside only */
	if((status = FastInsideAlignHB(cm, errbuf, cur_dsq, 1, L, cm->hbmx, &sc)) != eslOK) return status; /* errbuf will have been filled by FastInsideAlignHB() */
      }
      else { /* non-banded inside only */
	if((status = FastInsideAlign(cm, errbuf, cur_dsq, 1, L, NULL, &sc)) != eslOK) return status; 
      }
    }
    else if (do_small) { /* small D&C CYK alignment */
      if(do_qdb) { /* use QDBs when doing D&C CYK */
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, cm->dmin, cm->dmax);
	if(bdump_level > 0) qdb_trace_info_dump(cm, *cur_tr, cm->dmin, cm->dmax, bdump_level);
      }
      else { /* small D&C CYK non-banded alignment */
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL); /* we're not in QDB mode */
	if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level); /* informative as to where optimal parse goes outside bands */
      }
    }
    else if(do_qdb) { /* non-small, QDB banded CYK alignment */
      sc = CYKInside(cm, cur_dsq, L, 0, 1, L, cur_tr, cm->dmin, cm->dmax);
      if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
    }
    else { /* non-small, non-QDB CYK or optimal accuracy alignment */
      if(do_hbanded) { 
	if(do_post) { /* HMM banded CYK or optimal accuracy, posterior annotated */
	  if((status = FastAlignHB(cm, errbuf, cur_dsq, L, 1, L, cm->hbmx, do_optacc, out_mx, cur_tr, &(postcode1[i]), &(postcode2[i]), &sc)) != eslOK) return status;
	}
	else { /* HMM banded CYK or optimal accuracy, no posteriors */
	  if((status = FastAlignHB(cm, errbuf, cur_dsq, L, 1, L, cm->hbmx, do_optacc, out_mx, cur_tr, NULL, NULL, &sc)) != eslOK) {
	    if(do_optacc) return status; /* we can't handle a memory overload if we're trying to do optimal accuracy */
	    else if (status == eslERANGE) { /* we can still do CYK D&C alignment with QDBs derived from the HMM bands */
	      hd2safe_hd_bands(cm->M, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);
	      printf("Doing QDB D&C because HMM banded parse of seq %d was too memory intensive.\n", i);
	      sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, cp9b->safe_hdmin, cp9b->safe_hdmax);
	    }
	    else return status; /* get here (!do_optacc) && FastAlignHB() returned status other than eslOK and eslERANGE */
	  }
	  /* SHOULD I UNCOMMENT THIS?: 
	   * if we're aligning search hits, and we're !do_optacc, we realign if the HMM banded parse was > 0.01 bits less than the search score for this hit */
	  /* if(dsq_mode && (! cm->search_opts & CM_SEARCH_INSIDE)) {
	    if(!do_optacc && (sc - search_results->data[i].score) < 0.01) {
	      ESL_DPRINTF1(("Realigning hit: %d with D&C b/c HMM banded parse (%.3f bits) too-far-off search score (%.3f bits).\n", sc, search_results->data[i].score));
	      FreeParsetree(*cur_tr);
	      sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL);
	    }
	    }*/
	  if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level); /* allows you to see where the non-banded parse went outside the bands. */
	} 
	/* if CM_ALIGN_HMMSAFE option is enabled, realign seqs w/HMM banded parses < 0 bits,
	 * this should never happen in we're doing optimal accuracy or appending posteriors, due to option checking in cmalign, cmscore,
	 * but we check here to be safe */
	if(cm->align_opts & CM_ALIGN_HMMSAFE && sc < 0.) { 
	  if(do_post)   ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets() cm->align_opts option CM_ALIGN_HMMSAFE is ON at same time as incompatible option CM_ALIGN_POST.\n");
	  if(do_optacc) ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets() cm->align_opts option CM_ALIGN_HMMSAFE is ON at same time as incompatible option CM_ALIGN_OPTACC.\n");
	  tmpsc = sc;
	  if(!silent_mode) printf("\n%s HMM banded parse had a negative score, realigning with non-banded CYK.\n", seqs_to_aln->sq[i]->name);
	  FreeParsetree(*cur_tr);
	  sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL); /* we're not in QDB mode */
	  if(!silent_mode && fabs(sc-tmpsc) < 0.01) printf("HMM banded parse was the optimal parse.\n\n");
	  else if (!silent_mode) printf("HMM banded parse was non-optimal, it was %.2f bits below the optimal.\n\n", (fabs(sc-tmpsc)));
	}
      }
      else { 
	if(do_post) { /* non-banded CYK or optimal accuracy, posterior annotated */
	  if((status = FastAlign(cm, errbuf, cur_dsq, L, 1, L, &alpha, do_optacc, &beta, cur_tr, &(postcode1[i]), &(postcode2[i]), &sc)) != eslOK) return status;
	}
	else { /* non-banded CYK or optimal accuracy, no posteriors */
	  if((status = FastAlign(cm, errbuf, cur_dsq, L, 1, L, &alpha, do_optacc, &beta, cur_tr, NULL, NULL, &sc)) != eslOK) return status;
	}
	if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level); /* allows you to see where the non-banded parse went outside the bands. */
      } 
    }
    /* end of large if() else if() else if() else statement */
    /* done alignment for this seq */

    avgsc += sc;
    if (sc > maxsc) maxsc = sc;
    if (sc < minsc) minsc = sc;
      
    if(!silent_mode) { 
      if(!do_optacc) printf("    score: %10.2f bits\n", sc);
      else           printf("    score: %14.6f average posterior probability\n", sc);
    }
    parsesc[i] = sc;

    /* check parsetree score if cm->align_opts & CM_ALIGN_CHECKPARSESC */
    if((cm->align_opts & CM_ALIGN_CHECKPARSESC) && (!(cm->flags & CM_IS_SUB))) { 
      if(do_optacc) 
	ESL_FAIL(eslEINCOMPAT, errbuf, "ActuallyAlignTargets(), cm->align_opts CM_ALIGN_CHECKPARSESC, is on, but incompatible with another enabled option: CM_ALIGN_OPTACC.\n");
      if (fabs(sc - ParsetreeScore(cm, tr[i], cur_dsq, FALSE)) >= 0.01)
	ESL_FAIL(eslFAIL, errbuf, "ActuallyAlignTargets(), seq: %d alignment score %.3f differs from its parse tree's score: %.3f", i, sc, ParsetreeScore(cm, tr[i], cur_dsq, FALSE));
    }

    /* If requested, or if debug level high enough, print out the parse tree */
    if((cm->align_opts & CM_ALIGN_PRINTTREES) || (debug_level > 2)) { 
      fprintf(stdout, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], cur_dsq, FALSE));;
      ParsetreeDump(stdout, tr[i], cm, cur_dsq, NULL, NULL);
      fprintf(stdout, "//\n");
    }
    /* Dump the trace with info on i, j and d bands
     * if bdump_level is high enough */
    if(bdump_level > 0 && do_hbanded)
      ijdBandedTraceInfoDump(cm, tr[i], cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 1);

    if(do_sub) {
      if(! do_inside) { 
	/* Convert the sub_cm parsetree to a full CM parsetree */
	if(debug_level > 0) ParsetreeDump(stdout, *cur_tr, cm, cur_dsq, NULL, NULL);
	if(!(sub_cm2cm_parsetree(orig_cm, sub_cm, &orig_tr, *cur_tr, submap, debug_level))) { 
	  /* ParsetreeDump(stdout, orig_tr, orig_cm, cur_dsq, NULL, NULL); */
	  ESL_FAIL(eslFAIL, errbuf, "ActuallyAlignTargets(), Unable to convert sub CM parsetree to original CM parsetree. This shouldn't happen.");
	}
	if(debug_level > 0) { 
	  printf("\n\nConverted original trace:\n");
	  ParsetreeDump(stdout, orig_tr, orig_cm, cur_dsq, NULL, NULL);
	}
	/* Replace the sub_cm trace with the converted orig_cm trace. */
	FreeParsetree(*cur_tr);
	*cur_tr = orig_tr;
      }
      /* free sub_cm variables, we build a new sub CM for each seq */
      if(out_mx != NULL) { cm_hb_mx_Destroy(out_mx); out_mx = NULL; }

      FreeSubMap(submap);
      FreeCM(sub_cm); /* cm and sub_cm now point to NULL */
    }
    /* free alpha and beta, we need to allocate new ones for each seq */
    if(alpha != NULL)  { free_vjd_matrix(alpha, cm->M, 1, L); alpha = NULL; }
    if(beta  != NULL)  { free_vjd_matrix(beta,  cm->M, 1, L); beta = NULL; }
    if(do_timings) { 
      esl_stopwatch_Stop(watch); 
      esl_stopwatch_Display(stdout, watch, "seq alignment CPU time: ");
      printf("\n");
    }
  }
  /* done aligning all nalign seqs. */
  /* Clean up. */
  if(out_mx != NULL) cm_hb_mx_Destroy(out_mx);
  if(alpha != NULL)  { free_vjd_matrix(alpha, cm->M, 1, L); alpha = NULL; }
  if(beta  != NULL)  { free_vjd_matrix(beta,  cm->M, 1, L); beta = NULL; }
  if (do_qdb) {
    free(orig_dmin);
    free(orig_dmax);
  }
  esl_stopwatch_Destroy(watch);
  
  if(sq_mode) {
    seqs_to_aln->tr       = tr;       /* could be NULL */
    seqs_to_aln->cp9_tr   = cp9_tr;   /* could be NULL */
    seqs_to_aln->postcode1= postcode1;/* could be NULL */
    seqs_to_aln->postcode2= postcode2;/* could be NULL */
    seqs_to_aln->sc       = parsesc;  /* shouldn't be NULL */
  }
  else { /* dsq mode */
    if(do_post) { 
      for (i = 0; i < nalign; i++) {
	search_results->data[i].pcode1 = postcode1[i];
	search_results->data[i].pcode2 = postcode2[i];
      }
      /* we've copied the 1D postcode ptrs, free the 2D, ptr to the ptrs */
      free(postcode1);
      free(postcode2);
    }
    free(parsesc);
  }
  
  return eslOK;
  ERROR:
  ESL_FAIL(eslEMEM, errbuf, "ActuallyAlignTargets(), Memory allocation error.");
  return status; /* NEVERREACHED */
}

/*
 * Function: read_next_search_seq
 *
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easeled: EPN, Fri Dec  8 11:40:20 2006
 *
 * Purpose:  Given a dbfp and whether or not to take the reverse complement,
 *           reads in the next sequence and prepares reverse complement.
 *
 * Returns:  eslOK on success; eslEOF if end of file, 
 *           some other status code from esl_sqio_Read() if an error occurs.
 */
int read_next_search_seq (const ESL_ALPHABET *abc, ESL_SQFILE *dbfp, int do_revcomp, dbseq_t **ret_dbseq) 
{
  int status;
  dbseq_t *dbseq = NULL;

  ESL_ALLOC(dbseq, sizeof(dbseq_t));
  dbseq->sq[0] = NULL;
  dbseq->sq[1] = NULL;

  dbseq->sq[0] = esl_sq_CreateDigital(abc);

  while((status = esl_sqio_Read(dbfp, dbseq->sq[0])) == eslOK && (dbseq->sq[0]->n == 0)) /* skip zero length seqs */
    esl_sq_Reuse(dbseq->sq[0]);

  if(status != eslOK) goto ERROR;

  if (do_revcomp)
    {
      /* make a new ESL_SQ object, to store the reverse complement */
      if((dbseq->sq[1] = esl_sq_CreateDigitalFrom(abc, dbseq->sq[0]->name, dbseq->sq[0]->dsq, 
						  dbseq->sq[0]->n, dbseq->sq[0]->desc, 
						  dbseq->sq[0]->acc, dbseq->sq[0]->ss)) == NULL) goto ERROR;
      /* reverse complement it in place */
      revcomp(dbseq->sq[1]->abc, dbseq->sq[1], dbseq->sq[1]);
    }
  dbseq->results[0] = NULL;
  dbseq->results[1] = NULL;

  *ret_dbseq = dbseq;
  return eslOK;

 ERROR:
  if(dbseq->sq[0] != NULL) esl_sq_Destroy(dbseq->sq[0]);
  if(dbseq->sq[1] != NULL) esl_sq_Destroy(dbseq->sq[1]);
  if(dbseq != NULL) free(dbseq);
  return status;
}

/*
 * Function: CreateSeqsToAln()
 * Date:     EPN, Sat Sep  1 10:51:28 2007
 *
 * Purpose:  Allocate and return a seqs_to_aln_t data structure.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  An initialized and allocated (for nalloc seqs) 
 *           seqs_to_aln_t object.
 *           Dies immediately on a memory error.
 */
seqs_to_aln_t *CreateSeqsToAln(int size, int i_am_mpi_master)
{
  int status;
  int i;
  seqs_to_aln_t *seqs_to_aln;

  ESL_ALLOC(seqs_to_aln, sizeof(seqs_to_aln_t));
  ESL_ALLOC(seqs_to_aln->sq,      sizeof(ESL_SQ *)      * size);
  seqs_to_aln->tr       = NULL;
  seqs_to_aln->cp9_tr   = NULL;
  seqs_to_aln->postcode1= NULL;
  seqs_to_aln->postcode2= NULL;
  seqs_to_aln->sc       = NULL;
  seqs_to_aln->nalloc   = size;
  seqs_to_aln->nseq     = 0;

  if(i_am_mpi_master) {
    ESL_ALLOC(seqs_to_aln->tr,       sizeof(Parsetree_t *) * size);
    ESL_ALLOC(seqs_to_aln->cp9_tr,   sizeof(CP9trace_t)    * size);
    ESL_ALLOC(seqs_to_aln->postcode1,sizeof(char **)       * size);
    ESL_ALLOC(seqs_to_aln->postcode2,sizeof(char **)       * size);
    ESL_ALLOC(seqs_to_aln->sc,       sizeof(float)         * size);
    for(i = 0; i < size; i++) {
      seqs_to_aln->tr[i]       = NULL;
      seqs_to_aln->cp9_tr[i]   = NULL;
      seqs_to_aln->postcode1[i]= NULL;
      seqs_to_aln->postcode2[i]= NULL;
      seqs_to_aln->sc[i]       = IMPOSSIBLE;
    }
  }
  return seqs_to_aln;

 ERROR:
  cm_Fail("Memory error.");
  return NULL; /* NEVERREACHED */
}

/*
 * Function: CreateSeqsToAlnFromSq()
 * Date:     EPN, Wed Sep  5 18:13:11 2007
 *
 * Purpose:  Allocate and return a seqs_to_aln_t data structure, setting
 *           seqs_to_aln->sq ptr to an input ptr to ESL_SQs.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  An initialized and allocated (for nalloc seqs) 
 *           seqs_to_aln_t object.
 *           Dies immediately on a memory error.
 */
seqs_to_aln_t *CreateSeqsToAlnFromSq(ESL_SQ **sq, int size, int i_am_mpi_master)
{
  int status;
  int i;
  seqs_to_aln_t *seqs_to_aln;

  ESL_ALLOC(seqs_to_aln, sizeof(seqs_to_aln_t));
  seqs_to_aln->sq       = sq;
  seqs_to_aln->tr       = NULL;
  seqs_to_aln->cp9_tr   = NULL;
  seqs_to_aln->postcode1= NULL;
  seqs_to_aln->postcode2= NULL;
  seqs_to_aln->sc       = NULL;
  seqs_to_aln->nalloc   = size;
  seqs_to_aln->nseq     = size;

  if(i_am_mpi_master) {
    ESL_ALLOC(seqs_to_aln->tr,       sizeof(Parsetree_t *) * size);
    ESL_ALLOC(seqs_to_aln->cp9_tr,   sizeof(CP9trace_t)    * size);
    ESL_ALLOC(seqs_to_aln->postcode1,sizeof(char **)       * size);
    ESL_ALLOC(seqs_to_aln->postcode2,sizeof(char **)       * size);
    ESL_ALLOC(seqs_to_aln->sc,       sizeof(float)         * size);
    for(i = 0; i < size; i++) {
      seqs_to_aln->tr[i]       = NULL;
      seqs_to_aln->cp9_tr[i]   = NULL;
      seqs_to_aln->postcode1[i]= NULL;
      seqs_to_aln->postcode2[i]= NULL;
      seqs_to_aln->sc[i]       = IMPOSSIBLE;
    }
  }
  return seqs_to_aln;

 ERROR:
  cm_Fail("Memory error.");
  return NULL; /* NEVERREACHED */
}

/*
 * Function: GrowSeqsToAln()
 * Date:     EPN, Sat Sep  1 11:10:22 2007
 *
 * Purpose:  Grow a seqs_to_aln_t object by <new_alloc>.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  eslOK;
 */
int GrowSeqsToAln(seqs_to_aln_t *seqs_to_aln, int new_alloc, int i_am_mpi_master)
{
  int status;
  void *tmp;
  int i;

  ESL_RALLOC(seqs_to_aln->sq, tmp, sizeof(ESL_SQ *) * (seqs_to_aln->nalloc + new_alloc)); 

  if(i_am_mpi_master) {
    ESL_RALLOC(seqs_to_aln->tr,       tmp, sizeof(Parsetree_t *) * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->cp9_tr,   tmp, sizeof(CP9trace_t)    * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->postcode1,tmp, sizeof(char **)       * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->postcode2,tmp, sizeof(char **)       * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->sc,       tmp, sizeof(float)         * (seqs_to_aln->nalloc + new_alloc));
    for(i = seqs_to_aln->nalloc; i < (seqs_to_aln->nalloc + new_alloc); i++) {
      seqs_to_aln->tr[i]       = NULL;
      seqs_to_aln->cp9_tr[i]   = NULL;
      seqs_to_aln->postcode1[i]= NULL;
      seqs_to_aln->postcode2[i]= NULL;
      seqs_to_aln->sc[i]       = IMPOSSIBLE;
    }
  }
  
  seqs_to_aln->nalloc += new_alloc;
  return eslOK;

 ERROR:
  cm_Fail("Memory reallocation error.");
  return status; /* NEVERREACHED */
}


/*
 * Function: FreeSeqsToAln()
 *
 * Date:     EPN, Sat Sep  1 11:18:39 2007
 *
 * Purpose:  Free a seqs_to_aln_t object.
 *
 * Returns:  void
 *
 */
void FreeSeqsToAln(seqs_to_aln_t *s) 
{
  int i;
  
  if(s->sq != NULL) /* with MPI workers, we sometimes free the sequences outside this function */
    {
      for (i=0; i < s->nseq; i++) 
	if(s->sq[i] != NULL) esl_sq_Destroy(s->sq[i]);
      free(s->sq);
    }

  if(s->tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->tr[i] != NULL) FreeParsetree(s->tr[i]);
    free(s->tr);
  }

  if(s->cp9_tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->cp9_tr[i] != NULL) CP9FreeTrace(s->cp9_tr[i]);
    free(s->cp9_tr);
  }
 
  if(s->postcode1 != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->postcode1[i] != NULL) free(s->postcode1[i]);
    free(s->postcode1);
  }

  if(s->postcode2 != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->postcode2[i] != NULL) free(s->postcode2[i]);
    free(s->postcode2);
  }

  if(s->sc != NULL) free(s->sc);
  
  free(s);
}

/*
 * Function: FreePartialSeqsToAln()
 *
 * Date:     EPN, Wed Sep  5 06:58:39 2007
 *
 * Purpose:  Free specified parts of a seqs_to_aln_t object. 
 *
 * Returns:  void
 *
 */
void FreePartialSeqsToAln(seqs_to_aln_t *s, int do_free_sq, int do_free_tr, int do_free_cp9_tr, int do_free_post, int do_free_sc) 
{
  int i;
  
  if(do_free_sq && s->sq != NULL) {
    for (i=0; i < s->nseq; i++) 
      if(s->sq[i] != NULL) esl_sq_Destroy(s->sq[i]);
    free(s->sq);
    s->sq = NULL;
  }

  if(do_free_tr && s->tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->tr[i] != NULL) FreeParsetree(s->tr[i]);
    free(s->tr);
    s->tr = NULL;
  }

  if(do_free_cp9_tr && s->cp9_tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->cp9_tr[i] != NULL) CP9FreeTrace(s->cp9_tr[i]);
    free(s->cp9_tr);
    s->cp9_tr = NULL;
  }
 
  if(do_free_post && s->postcode1 != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->postcode1[i] != NULL) free(s->postcode1[i]);
    free(s->postcode1);
    s->postcode1 = NULL;
  }

  if(do_free_post && s->postcode2 != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->postcode2[i] != NULL) free(s->postcode2[i]);
    free(s->postcode2);
    s->postcode2 = NULL;
  }

  if(do_free_sc && s->sc != NULL) {
    free(s->sc);
    s->sc = NULL;
  }
}

/*
 * Function: ReadSeqsToAln()
 * Date:     EPN, Fri Aug 31 15:20:37 2007
 *
 * Purpose:  Given a pointer to a seq file we're reading seqs to align
 *           from, read in nseq seqs from the seq file, or 
 *           if nseq == 0 && do_real_all == TRUE, read all the seqs.
 *           Add the sequences to a growing seqs_to_aln_t object,
 *           a pointer to which is passed in.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  <eslOK> on success with <*ret_seqs_to_aln_t> filled with 
 *           seqs to align, *ret_seqs_to_aln_t->nseq gives number of seqs.
 *           Dies immediately on failure with informative error message.
 */
int ReadSeqsToAln(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int nseq, int do_read_all, seqs_to_aln_t *seqs_to_aln, int i_am_mpi_master) 
{
  int status;
  int keep_reading = TRUE;
  int i;
  int nseq_orig;

  /* contract check */
  if(  do_read_all && nseq != 0) cm_Fail("if do_read_all is TRUE,  nseq must be zero.");
  if(! do_read_all && nseq <= 0) cm_Fail("if do_read_all is FALSE, nseq must be at least 1.");

  nseq_orig = seqs_to_aln->nseq;
  i         = seqs_to_aln->nseq;
  if(i == seqs_to_aln->nalloc) GrowSeqsToAln(seqs_to_aln, 100, i_am_mpi_master);

  seqs_to_aln->sq[i] = esl_sq_CreateDigital(abc);
  while (keep_reading && (status = esl_sqio_Read(seqfp, seqs_to_aln->sq[i])) == eslOK) {
    if(seqs_to_aln->sq[i]->n == 0) { esl_sq_Reuse(seqs_to_aln->sq[i]); continue; }
    i++;
    if(i == seqs_to_aln->nalloc) GrowSeqsToAln(seqs_to_aln, 100, i_am_mpi_master);
    if(! do_read_all && (i - nseq_orig) == nseq)   keep_reading = FALSE; 
    seqs_to_aln->sq[i] = esl_sq_CreateDigital(abc);
  }
  /* destroy the last sequence that was alloc'ed but not filled */
  esl_sq_Destroy(seqs_to_aln->sq[i]);
  if ((  do_read_all && status  != eslEOF) || 
      (! do_read_all && (status != eslEOF && status != eslOK)))
    cm_Fail("Parse failed, line %d, file %s:\n%s", 
	    seqfp->linenumber, seqfp->filename, seqfp->errbuf);

  seqs_to_aln->nseq = i;
  return status;

}

/*
 * Function: CMEmitSeqsToAln()
 * Date:     EPN, Tue Sep  4 13:22:11 2007   
 *
 * Purpose:  Create a seqs_to_aln object by generating sequences
 *           from a CM.
 *
 * Note:     Sequences are allocated slightly different if the MPI master
 *           calls this function, to allow us to store them after receiving
 *           them back from workers in any order.
 *
 * Args:     r               - source of randomness
 *           cm              - CM to emit from
 *           ncm             - number for CM (only for naming seqs if cm->name == NULL)
 *           nseq            - number of seqs to emit
 *           i_am_mpi_master - TRUE if called from MPI master (see Note)
 *
 * Returns:  Ptr to a newly allocated seqs_to_aln object with nseq sequences to align.
 *           Dies immediately on failure with informative error message.
 */
seqs_to_aln_t *CMEmitSeqsToAln(ESL_RANDOMNESS *r, CM_t *cm, int ncm, int nseq, int i_am_mpi_master) 
{
  int status;
  seqs_to_aln_t *seqs_to_aln = NULL;
  char *name = NULL;
  int namelen;
  int L;
  int i;

  seqs_to_aln = CreateSeqsToAln(nseq, i_am_mpi_master);

  namelen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
  if(cm->name != NULL) namelen += strlen(cm->name) + 1;
  ESL_ALLOC(name, sizeof(char) * namelen);

  for(i = 0; i < nseq; i++)
    {
      if(cm->name != NULL) sprintf(name, "%s-%d", cm->name, i+1);
      else                 sprintf(name, "%d-%d", ncm, i+1);
      EmitParsetree(cm, r, name, TRUE, NULL, &(seqs_to_aln->sq[i]), &L);
    }
  seqs_to_aln->nseq = nseq;

  free(name);
  return seqs_to_aln;


 ERROR:
  cm_Fail("memory allocation error");
  return NULL;
}

/*
 * Function: RandomEmitSeqsToAln()
 * Date:     EPN, Tue Sep  4 13:42:16 2007
 *
 * Purpose:  Create a seqs_to_aln object by generating sequences
 *           randomly from a background distro.
 *
 * Note:     Sequences are allocated slightly different if the MPI master
 *           calls this function, to allow us to store them after receiving
 *           them back from workers in any order.
 *
 * Args:     r               - source of randomness
 *           abc             - alphabet 
 *           pdist           - probability distribution to use for emitting
 *           extranum        - use this as first part of sequence name (could be ncm)
 *           nseq            - number of seqs to emit
 *           L               - length to make sequences 
 *           i_am_mpi_master - TRUE if called from MPI master (see Note)
 *
 * Returns:  Ptr to a newly allocated seqs_to_aln object with nseq sequences to align.
 *           Dies immediately on failure with informative error message.
 */
seqs_to_aln_t *RandomEmitSeqsToAln(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, double *pdist, int extranum, int nseq, int L, int i_am_mpi_master) 
{
  int status;
  seqs_to_aln_t *seqs_to_aln = NULL;
  char *name = NULL;
  int namelen;
  int i;
  ESL_DSQ *randdsq = NULL;

  seqs_to_aln = CreateSeqsToAln(nseq, i_am_mpi_master);
  ESL_ALLOC(randdsq,      sizeof(ESL_DSQ)* (L+2));

  namelen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
  ESL_ALLOC(name, sizeof(char) * namelen);

  for(i = 0; i < nseq; i++)
    {
      sprintf(name, "randseq%d-%d", extranum, i+1);
      if (esl_rnd_xIID(r, pdist, abc->K, L, randdsq)  != eslOK) cm_Fail("RandomEmitSeqsToAln(): failure creating random sequence.");
      if((seqs_to_aln->sq[i] = esl_sq_CreateDigitalFrom(abc, name, randdsq, L, NULL, NULL, NULL)) == NULL) 
	 cm_Fail("RandomEmitSeqsToAln() error.");
    }
  seqs_to_aln->nseq = nseq;

  free(name);
  free(randdsq);
  return seqs_to_aln;

 ERROR:
  cm_Fail("memory allocation error");
  return NULL;
}

/*
 * Function: print_results
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
 *           do_complement       are we doing the minus strand
 */
void print_results (CM_t *cm, SearchInfo_t *si, const ESL_ALPHABET *abc, CMConsensus_t *cons, dbseq_t *dbseq, int do_complement)
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

  /* Contract check: we allow the caller to specify the alphabet they want the 
   * resulting MSA in, but it has to make sense (see next few lines). */
  if(cm->abc->type == eslRNA) { 
      if(abc->type != eslRNA && abc->type != eslDNA)
	cm_Fail("print_results(), cm alphabet is RNA, but requested output alphabet is neither DNA nor RNA.");
  }
  else if(cm->abc->K != abc->K) cm_Fail("print_results(), cm alphabet size is %d, but requested output alphabet size is %d.", cm->abc->K, abc->K);
  if(si == NULL) cm_Fail("print_results(), si == NULL.\n");
  if(si->stype[si->nrounds] != SEARCH_WITH_HMM && si->stype[si->nrounds] != SEARCH_WITH_CM) cm_Fail("print_results(), final search round is neither SEARCH_WITH_HMM nor SEARCH_WITH_CM.\n");

  do_stats = (si->cutoff_type[si->nrounds] == E_CUTOFF) ? TRUE : FALSE;
  if(do_stats  && !(cm->flags & CMH_GUMBEL_STATS)) cm_Fail("print_results(), stats wanted but CM has no Gumbel stats\n");

  if(do_stats) { /* Determine Gumbel mode to use */
    CM2Gumbel_mode(cm, cm->search_opts, &cm_gum_mode, &cp9_gum_mode);
    gum = (si->stype[si->nrounds] == SEARCH_WITH_HMM) ? cm->stats->gumAA[cp9_gum_mode] : cm->stats->gumAA[cm_gum_mode];
  }
  emap = CreateEmitMap(cm);
  name = dbseq->sq[0]->name;
  len  = dbseq->sq[0]->n;

  for (in_revcomp = 0; in_revcomp <= do_complement; in_revcomp++) {
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
		      in_revcomp);
	FreeFancyAli(ali);
	printf ("\n");
      }
    }
  }
  fflush(stdout);
  FreeEmitMap(emap);
}

/*
 * Function: remove_hits_over_e_cutoff
 * Date:     RJK, Tue Oct 8, 2002 [St. Louis]
 * Purpose:  Given an E-value cutoff, lambdas, mus, a sequence, and
 *           a list of results, calculates GC content for each hit, 
 *           calculates E-value, and decides whether to keep hit or not.
 * 
 * Args:    
 *           cm      - the covariance model
 *           si      - SearchInfo, relevant round is final one, si->nrounds
 *           results - the hits data structure
 *           seq     - seq hits lie within, needed to determine gc content
 *           used_HMM- TRUE if hits are to the CM's CP9 HMM, not the CMa
 */
void remove_hits_over_e_cutoff (CM_t *cm, SearchInfo_t *si, search_results_t *results, ESL_SQ *sq)
{
  int gc_comp;
  int i, x;
  search_result_node_t swap;
  float score_for_Eval; /* the score we'll determine the statistical signifance
			 * of. This will be the bit score stored in
			 * dbseq unless (cm->flags & CM_ENFORCED)
			 * in which case we subtract cm->enf_scdiff first. */
  int cm_gum_mode;      /* Gumbel mode if we're using CM hits */
  int cp9_gum_mode;     /* Gumbel mode if we're using HMM hits */
  int p;                /* relevant partition */
  GumbelInfo_t **gum;      /* pointer to gum to use */
  float cutoff;         /* the max E-value we want to keep */

  /* Check contract */
  if(!(cm->flags & CMH_GUMBEL_STATS)) cm_Fail("remove_hits_over_e_cutoff(), but CM has no gumbel stats\n");
  if(!(sq->flags & eslSQ_DIGITAL))    cm_Fail("remove_hits_over_e_cutoff(), sequences is not digitized.\n");
  if(si == NULL) cm_Fail("remove_hits_over_e_cutoff(), si == NULL.\n");
  if(si->stype[si->nrounds] != SEARCH_WITH_HMM && si->stype[si->nrounds] != SEARCH_WITH_CM) cm_Fail("remove_hits_over_e_cutoff(), final search round is neither SEARCH_WITH_HMM nor SEARCH_WITH_CM.\n");

  if (results == NULL) return;

  /* Determine Gumbel mode to use */
  CM2Gumbel_mode(cm, si->search_opts[si->nrounds], &cm_gum_mode, &cp9_gum_mode);
  gum = (si->stype[si->nrounds] == SEARCH_WITH_HMM) ? cm->stats->gumAA[cp9_gum_mode] : cm->stats->gumAA[cm_gum_mode];
  
  ESL_DASSERT1((si->cutoff_type[si->nrounds] == E_CUTOFF));
  cutoff = si->cutoff[si->nrounds];
  
  for (i=0; i<results->num_results; i++) {
    gc_comp = get_gc_comp (sq, results->data[i].start, results->data[i].stop);
    p = cm->stats->gc2p[gc_comp];
    score_for_Eval = results->data[i].score;
    if(cm->flags & CM_ENFORCED) {
      /*printf("\n\tRM orig sc: %.3f", score_for_Eval);*/
      score_for_Eval -= cm->enf_scdiff;
      /*printf(" new sc: %.3f (diff: %.3f\n", score_for_Eval, cm->enf_scdiff);*/
    }
    /*printf("score_for_Eval: %f \n", score_for_Eval);*/
    if (RJK_ExtremeValueE(score_for_Eval, gum[p]->mu, gum[p]->lambda) > cutoff)  
      results->data[i].start = -1;
    /*printf("Eval: %f, start: %d\n", RJK_ExtremeValueE(score_for_Eval,
      mu[gc_comp], lambda[gc_comp]),
      results->data[i].start);*/
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
  sort_results(results);
}  

/* Function: revcomp()
 * Incept:   EPN, Tue Aug  7 10:05:14 2007
 *           based on Squid's revcomp()
 *
 * Purpose:  Reverse complement ESL_SQ seq; store in comp.
 *           Can revcomp "in place" (revcomp(seq, seq)).
 *           sq can be in digital or text form.
 *
 * Args:     comp  - destination for reverse complement of sq
 *           seq   - sequence to reverse complement
 *
 * Returns:  eslOK on success;
 *           Dies immediately if any error occurs.
 */
int
revcomp(const ESL_ALPHABET *abc, ESL_SQ *comp, ESL_SQ *sq)
{
  int status;
  int do_digital = FALSE;
  int i;

  /* contract checks */
  if (comp == NULL)
    cm_Fail("ERROR in revcomp, comp is NULL.");
  if(sq == NULL)
    cm_Fail("ERROR in revcomp, sq is NULL.");
  if(    sq->flags & eslSQ_DIGITAL   &&  (! (comp->flags & eslSQ_DIGITAL)))
    cm_Fail("ERROR in revcomp, sq is digital, comp is not.");
  if((! (sq->flags & eslSQ_DIGITAL)) &&      comp->flags & eslSQ_DIGITAL)
    cm_Fail("ERROR in revcomp, comp is digital, sq is not.");
  if(abc->type != eslRNA && abc->type != eslDNA)
    cm_Fail("ERROR in revcomp, alphabet type must be RNA or DNA.");
  if(comp->n < sq->n)
    cm_Fail("ERROR in revcomp, comp->n is smaller than sq->n.");

  if(sq->flags & eslSQ_DIGITAL) do_digital = TRUE;

  if(do_digital) {
    if((status = esl_rnd_XReverse(sq->dsq, sq->n, comp->dsq)) != eslOK) 
      goto ERROR; 
  }
  else {
    if((status = esl_rnd_CReverse(sq->seq, comp->seq) != eslOK))
      goto ERROR; 
  } 

  if(do_digital)
    {
      for(i = 1; i <= sq->n; i++)
	{
	  if(sq->dsq[i] >= abc->Kp) { status = eslEINVAL; goto ERROR; }
	  switch (abc->sym[sq->dsq[i]]) {
	  case 'A': 
	    if(abc->type == eslRNA) 
	      comp->dsq[i] = abc->inmap[(int) 'U']; 
	    else
	      comp->dsq[i] = abc->inmap[(int) 'T']; 
	    break;
	  case 'C': comp->dsq[i] = abc->inmap[(int) 'G']; break;
	  case 'G': comp->dsq[i] = abc->inmap[(int) 'C']; break;
	  case 'T': comp->dsq[i] = abc->inmap[(int) 'A']; break;
	  case 'U': comp->dsq[i] = abc->inmap[(int) 'A']; break;
	  case 'R': comp->dsq[i] = abc->inmap[(int) 'Y']; break;
	  case 'Y': comp->dsq[i] = abc->inmap[(int) 'R']; break;
	  case 'M': comp->dsq[i] = abc->inmap[(int) 'K']; break;
	  case 'K': comp->dsq[i] = abc->inmap[(int) 'M']; break;
	  case 'S': comp->dsq[i] = abc->inmap[(int) 'S']; break;
	  case 'W': comp->dsq[i] = abc->inmap[(int) 'W']; break;
	  case 'H': comp->dsq[i] = abc->inmap[(int) 'D']; break;
	  case 'D': comp->dsq[i] = abc->inmap[(int) 'H']; break;
	  case 'B': comp->dsq[i] = abc->inmap[(int) 'V']; break;
	  case 'V': comp->dsq[i] = abc->inmap[(int) 'B']; break;
	  default:  break;		/* anything else? leave it; it's prob a gap or an X */
	  }
	}
    }
  else
    {
      for(i = 0; i < sq->n; i++)
	{
	  if(islower(sq->seq[i])) { status = eslEINVAL; goto ERROR; }
	     switch (sq->seq[i]) {
	     case 'A': 
	       if(abc->type == eslRNA) 
		 comp->seq[i] = 'U';
	       else
		 comp->seq[i] = 'T'; 
	       break;
	     case 'C': comp->seq[i] = 'G'; break;
	     case 'G': comp->seq[i] = 'C'; break;
	     case 'T': comp->seq[i] = 'A'; break;
	     case 'U': comp->seq[i] = 'A'; break;
	     case 'R': comp->seq[i] = 'Y'; break;
	     case 'Y': comp->seq[i] = 'R'; break;
	     case 'M': comp->seq[i] = 'K'; break;
	     case 'K': comp->seq[i] = 'M'; break;
	     case 'S': comp->seq[i] = 'S'; break;
	     case 'W': comp->seq[i] = 'W'; break;
	     case 'H': comp->seq[i] = 'D'; break;
	     case 'D': comp->seq[i] = 'H'; break;
	     case 'B': comp->seq[i] = 'V'; break;
	     case 'V': comp->seq[i] = 'B'; break;
	     default:  break;		/* anything else? leave it; it's prob a gap or an X */
	     }
	}
    }
  return eslOK;

 ERROR: 
  cm_Fail("Unexpected error code: %d in revcomp().", status);
  return status; /* NOTREACHED */
}
    


/************************************************************
 * @LICENSE@
 ************************************************************/
/* dispatch.c
 * 
 * The two all-important dispatch functions:
 * DispatchSearch()     calls appropriate DP search functions.
 * DispatchAlignments() calls appropriate DP alignment functions.
 * 
 * EPN, Wed Dec  6 06:11:46 2006
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
#include "esl_msa.h"         
#include "esl_stopwatch.h"   

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* 
 * Function: DispatchSearch()
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
 *           size_limit      - max number of Mb for DP matrix, if matrix is bigger return eslERANGE 
 *           ret_flen        - RETURN: subseq len that survived filter (NULL if not filtering)
 *           ret_sc          - RETURN: Highest scoring hit from search (even if below cutoff).
 *
 * Returns: eslOK on success. eslERANGE if we're doing HMM banded alignment and requested matrix is too big.
 */
int DispatchSearch(CM_t *cm, char *errbuf, int sround, ESL_DSQ *dsq, int i0, int j0, search_results_t **results, float size_limit, int *ret_flen, float *ret_sc)
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
  if(!(cm->flags & CMH_BITS))          ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), CMH_BITS flag down.\n");
  if(si == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(): search info cm->si is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(): dsq is NULL.");
  if(!(cm->flags & CMH_BITS))          ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(): CMH_BITS flag down.\n");
  if(sround > si->nrounds)             ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(): current round %d is greater than cm->si->nrounds: %d\n", sround, si->nrounds);
  if(results[sround] == NULL)          ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(): results for current round %d are NULL\n", sround);

  ESL_DPRINTF1(("In DispatchSearch(), round: %d\n", sround));

  flen = (j0-i0+1);

  /* TEMPORARY */
  if(si->stype[sround] == SEARCH_WITH_HYBRID) cm_Fail("DispatchSearch, hybrid filtering not yet implemented.\n");

  /* copy info for this round from SearchInfo fi */
  cm->search_opts = si->search_opts[sround]; 
  cutoff          = si->sc_cutoff[sround]; /* this will be a bit score regardless of whether the cutoff_type == E_CUTOFF */
  stype           = si->stype[sround];
  smx             = si->smx[sround]; /* may be NULL */
  hsi             = si->hsi[sround]; /* may be NULL */
  round_results   = results[sround]; /* may not be NULL, contract enforced this */
  h_existing      = round_results->num_results; /* remember this, b/c we only want rescan survivors found in *this* function call */

  /* SEARCH_WITH_HMM section */
  if(stype == SEARCH_WITH_HMM) { 
    /* some SEARCH_WITH_HMM specific contract checks */
    if(cm->cp9 == NULL)                    ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), trying to use CP9 HMM that is NULL.\n");
    if(!(cm->cp9->flags & CPLAN9_HASBITS)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");
    if(smx != NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), round %d, SEARCH_WITH_HMM but smx != NULL.\n", sround);
    if(hsi != NULL)                        ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), round %d, SEARCH_WITH_HMM but hsi != NULL.\n", sround);
    if(! ((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD)))
      ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), search type for this round is SEARCH_WITH_HMM, but CM_SEARCH_HMMVITERBI and CM_SEARCH_HMMFORWARD flags are both down.");

    search_results_t *fwd_results;
    /* Scan the (sub)seq in forward direction w/Viterbi or Forward, getting j end points of hits above cutoff */
    fwd_results = CreateResults(INIT_RESULTS);
    if(cm->search_opts & CM_SEARCH_HMMVITERBI) { 
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, i0, j0, cm->W, cutoff, fwd_results, 
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       NULL, NULL, NULL,  /* don't return best score at each posn, best scoring posn, or traces */
			       &sc)) != eslOK) return status;
    }
    else if(cm->search_opts & CM_SEARCH_HMMFORWARD) { 
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, i0, j0, cm->W, cutoff, fwd_results,
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       NULL, NULL, /* don't return best score at each posn, or best scoring posn */
			       &sc)) != eslOK) return status;
    }
#if 0 /* EPN, Tue Jan 22 14:26:57 2008, we don't need to remove overlaps, they'll be collapsed below
       * (unless do_hbanded), and even if they weren't it's okay b/c we remove overlaps before
       * we print results in cmsearch.
       */
    /* Remove overlapping hits, if we're being greedy */
    if(cm->search_opts & CM_SEARCH_HMMGREEDY) { /* resolve overlaps by being greedy */
      ESL_DASSERT1((i0 == 1)); /* EPN, Tue Nov 27 13:59:31 2007 not sure why this is here */
      remove_overlapping_hits (fwd_results, i0, j0);
    }
#endif

    /* determine start points (i) of the hits based on backward direction (Viterbi or Backward) scan starting at j */
    for(h = 0; h < fwd_results->num_results; h++) {
      min_i = (fwd_results->data[h].stop - cm->W + 1) >= 1 ? (fwd_results->data[h].stop - cm->W + 1) : 1;
      if(cm->search_opts & CM_SEARCH_HMMVITERBI) { 
	if((status = cp9_ViterbiBackward(cm, errbuf, cm->cp9_mx, dsq, min_i, fwd_results->data[h].stop, cm->W, cutoff, 
					 round_results, /* report hits to this round's results */
					 TRUE,   /* we're scanning */
					 FALSE,  /* we're not ultimately aligning */
					 TRUE,   /* be memory efficient */
					 NULL, NULL, NULL,  /* don't return best score at each posn, best scoring posn, or traces */
					 &bwd_sc)) != eslOK) return status;
      }
      else { 
	if((status = cp9_Backward(cm, errbuf, cm->cp9_mx, dsq, min_i, fwd_results->data[h].stop, cm->W, cutoff, 
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
    if(smx != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), round %d, SEARCH_WITH_HYBRID but smx != NULL.\n", sround);
    if(hsi == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(): current round %d is type SEARCH_WITH_HYBRID, but hsi is NULL\n", sround);
    if((status = cm_cp9_HybridScan(cm, errbuf, cm->cp9_mx, dsq, hsi, i0, j0, hsi->W, cutoff, round_results, 
				   NULL, NULL, /* don't return best score at each posn, and best scoring posn */
				   &sc)) != eslOK) return status;
  }  
  else { /* stype == SEARCH_WITH_CM */
    ESL_DASSERT1((stype == SEARCH_WITH_CM));
    if(smx == NULL)                             ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), round %d, SEARCH_WITH_CM but smx == NULL.\n", sround);
    if(sround == si->nrounds && smx != cm->smx) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(), final round %d, SEARCH_WITH_CM but smx != cm->smx.\n", sround);
    if(hsi != NULL)                             ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchSearch(): round %d is type SEARCH_WITH_CM, but hsi is NULL\n", sround);

    if(cm->search_opts & CM_SEARCH_HBANDED) {
      if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, i0, j0, cm->cp9b, TRUE, 0)) != eslOK) return status; 
      if(cm->search_opts & CM_SEARCH_INSIDE) { if((status = FastFInsideScanHB(cm, errbuf, dsq, i0, j0, cutoff, round_results, cm->hbmx, size_limit, &sc)) != eslOK) return status; }
      else                                   { if((status = FastCYKScanHB    (cm, errbuf, dsq, i0, j0, cutoff, round_results, cm->hbmx, size_limit, &sc)) != eslOK) return status; }
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
	  ESL_DPRINTF1(("\tsucked in subseq: hit %d new_i: %d j (still): %d\n", h, i, j));
	}
      }
      /* next round: research this chunk that survived the filter */
      if((status = DispatchSearch(cm, errbuf, (sround+1), dsq, i, j, results, size_limit, NULL, NULL)) != eslOK) return status;
    }
  }
  else { /* we're done filtering, and we're done searching, get alignments if nec */
    if((round_results->num_results > 0) && (! (cm->search_opts & CM_SEARCH_NOALIGN))) {
      if((status = DispatchAlignments(cm, errbuf, NULL, 
					dsq, round_results, h_existing,     /* put function into dsq_mode, designed for aligning search hits */
					0, 0, 0, NULL, size_limit)) != eslOK) return status;
    }
  }
  if(ret_sc != NULL) *ret_sc = sc;
  return eslOK;
}  

/* 
 * Function: DispatchAlignments
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
 *           size_limit     - max number of Mb for a DP matrix, if requestd matrix is bigger return eslERANGE 
 * 
 * Returns:  eslOK on success;
 *           eslERANGE if required memory for a DP matrix is too big;
 *           eslEINCOMPAT if input parameters violate contract;
 *           eslEMEM on memory allocation error;
 *           eslFAIL if some other error;
 *           if(!eslOK) errbuf is filled with informative error message
 */
int
DispatchAlignments(CM_t *cm, char *errbuf, seqs_to_aln_t *seqs_to_aln, ESL_DSQ *dsq, search_results_t *search_results,
		   int first_result, int bdump_level, int debug_level, int silent_mode, ESL_RANDOMNESS *r, float size_limit)
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
  if(!(cm->flags & CMH_BITS))                            ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CMH_BITS flag down.\n");
  if(r == NULL && (cm->align_opts & CM_ALIGN_SAMPLE))    ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), no source of randomness, but CM_ALIGN_SAMPLE alignment option on.\n");
  if(r != NULL && (!(cm->align_opts & CM_ALIGN_SAMPLE))) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), we have a source of randomness, but CM_ALIGN_SAMPLE alignment option off.\n");
  if((cm->align_opts & CM_ALIGN_POST)      && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CM_ALIGN_POST and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CM_ALIGN_SCOREONLY and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_POST))       ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CM_ALIGN_SCOREONLY and CM_ALIGN_POST options are incompatible.\n");

  /* determine mode */
  if     (seqs_to_aln != NULL && (dsq == NULL && search_results == NULL))  sq_mode = TRUE;
  else if(seqs_to_aln == NULL && (dsq != NULL && search_results != NULL)) dsq_mode = TRUE;
  else   ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), can't determine mode (sq_mode or dsq_mode).\n");

  if( sq_mode && (seqs_to_aln->sq        == NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in sq_mode, seqs_to_aln->sq is NULL.\n");
  if( sq_mode && (seqs_to_aln->tr        != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in sq_mode, seqs_to_aln->tr is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->cp9_tr    != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in sq_mode, seqs_to_aln->cp9_tr is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->postcode1 != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in sq_mode, seqs_to_aln->postcode1 is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->postcode2 != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in sq_mode, seqs_to_aln->postcode2 is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->sc        != NULL))  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in sq_mode, seqs_to_aln->sc is non-NULL.\n");
  
  if(dsq_mode && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in dsq_mode, CM_ALIGN_HMMVITERBI option on.\n");
  if(dsq_mode && (cm->align_opts & CM_ALIGN_INSIDE))     ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in dsq_mode, CM_ALIGN_INSIDE option on.\n");
  if(dsq_mode && (cm->align_opts & CM_ALIGN_SAMPLE))     ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in dsq_mode, CM_ALIGN_SAMPLE option on.\n");
  if(dsq_mode && search_results == NULL)                 ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in dsq_mode, search_results are NULL.\n");
  if(dsq_mode && (first_result > search_results->num_results)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), in dsq_mode, first_result: %d > search_results->num_results: %d\n", first_result, search_results->num_results);

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
  if((do_sample + do_inside + do_post + do_hmmonly + do_scoreonly) > 1) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), exactly 0 or 1 of the following must be TRUE (== 1):\n\tdo_sample = %d\n\tdo_inside = %d\n\t do_post = %d\n\tdo_hmmonly = %d\n\tdo_scoreonly = %d\n", do_sample, do_inside, do_post, do_hmmonly, do_scoreonly);

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
  if((watch = esl_stopwatch_Create()) == NULL) goto ERROR;

  if(do_hbanded || do_sub) { /* We need a CP9 HMM to build sub_cms */
    if(cm->cp9 == NULL)                    ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments, trying to use CP9 HMM that is NULL\n");
    if(cm->cp9b == NULL)                   ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments, cm->cp9b is NULL\n");
    if(!(cm->cp9->flags & CPLAN9_HASBITS)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments, trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");
    
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
    CPlan9SWConfig(cm->cp9, swentry, swexit, FALSE); /* FALSE means don't make I_0, D_1, I_M unreachable (like a local CM, undesirable for sub CM strategy)) */
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
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, cur_dsq, 1, L, L, 0., NULL,
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
	ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), unexpected error building a sub CM for seq %d.", i);
      /* Configure the sub_cm, the same as the cm, this will build a CP9 HMM if (do_hbanded), this will also:  */
      /* (4) Build a new CP9 HMM from the sub CM. */
      ConfigCM(sub_cm, NULL, NULL, FALSE);
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
	if((status = FastInsideAlignHB(cm, errbuf, cur_dsq, 1, L, size_limit, cm->hbmx, NULL)) != eslOK) return status; /* errbuf will have been filled by FastInsideAlignHB() */
	if((status = SampleFromInsideHB(r, cm, errbuf, cur_dsq, L, cm->hbmx, cur_tr, &sc)) != eslOK) return status; /* errbuf will have been filled by SampleFromInsideHB() */
      }
      else { /* sampling from inside matrix, but not HMM banded */
	if((status = FastInsideAlign(cm, errbuf, cur_dsq, 1, L, size_limit, &alpha, NULL)) != eslOK) return status; /* errbuf will have been filled by FastInsideAlign() */
	if((status = SampleFromInside(r, cm, errbuf, cur_dsq, L, alpha, cur_tr, &sc)) != eslOK) return status; /* errbuf will have been filled by SampleFromInside() */
      }
    }
    else if(do_inside) { 
      if(do_hbanded) { /* HMM banded inside only */
	if((status = FastInsideAlignHB(cm, errbuf, cur_dsq, 1, L, size_limit, cm->hbmx, &sc)) != eslOK) return status; /* errbuf will have been filled by FastInsideAlignHB() */
      }
      else { /* non-banded inside only */
	if((status = FastInsideAlign(cm, errbuf, cur_dsq, 1, L, size_limit, NULL, &sc)) != eslOK) return status; 
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
	  if((status = FastAlignHB(cm, errbuf, cur_dsq, L, 1, L, size_limit, cm->hbmx, do_optacc, out_mx, cur_tr, &(postcode1[i]), &(postcode2[i]), &sc)) != eslOK) return status;
	}
	else { /* HMM banded CYK or optimal accuracy, no posteriors */
	  if((status = FastAlignHB(cm, errbuf, cur_dsq, L, 1, L, size_limit, cm->hbmx, do_optacc, out_mx, cur_tr, NULL, NULL, &sc)) != eslOK) {
	    if(do_optacc) return status; /* we can't handle a memory overload if we're trying to do optimal accuracy */
	    else if (status == eslERANGE) { /* we can still do CYK D&C alignment with QDBs derived from the HMM bands */
	      hd2safe_hd_bands(cm->M, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);
	      ESL_DPRINTF1(("Doing D&C because HMM banded parse of seq %d was too memory intensive.\n", i));
	      printf("Doing D&C because HMM banded parse of seq %d was too memory intensive.\n", i); 
	      sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL); /* we're not in QDB mode */
	    }
	    else return status; /* get here (!do_optacc) && FastAlignHB() returned status other than eslOK and eslERANGE */
	  }
	  /* if we're aligning search hits, and we're !do_optacc, we realign if the HMM banded parse was > 0.01 bits less than the search score for this hit */
	  if(dsq_mode && (! (cm->search_opts & CM_SEARCH_INSIDE))) {
	    if((!do_optacc) && ((fabs(sc - search_results->data[i].score)) > 0.01)) {
	      ESL_DPRINTF1(("Realigning hit: %d with D&C b/c HMM banded parse (%.3f bits) too-far-off search score (%.3f bits).\n", i, sc, search_results->data[i].score));
	      printf("Realigning hit: %d with D&C b/c HMM banded parse (%.3f bits) too-far-off search score (%.3f bits).\n", i, sc, search_results->data[i].score);
	      FreeParsetree(*cur_tr);
	      sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL);
	    }
	  }
	} 
	/* if CM_ALIGN_HMMSAFE option is enabled, realign seqs w/HMM banded parses < 0 bits,
	 * this should never happen in we're doing optimal accuracy or appending posteriors, due to option checking in cmalign, cmscore,
	 * but we check here to be safe */
	if(cm->align_opts & CM_ALIGN_HMMSAFE && sc < 0.) { 
	  if(do_post)   ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments() cm->align_opts option CM_ALIGN_HMMSAFE is ON at same time as incompatible option CM_ALIGN_POST.\n");
	  if(do_optacc) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments() cm->align_opts option CM_ALIGN_HMMSAFE is ON at same time as incompatible option CM_ALIGN_OPTACC.\n");
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
	  if((status = FastAlign(cm, errbuf, cur_dsq, L, 1, L, size_limit, &alpha, do_optacc, &beta, cur_tr, &(postcode1[i]), &(postcode2[i]), &sc)) != eslOK) return status;
	}
	else { /* non-banded CYK or optimal accuracy, no posteriors */
	  if((status = FastAlign(cm, errbuf, cur_dsq, L, 1, L, size_limit, &alpha, do_optacc, &beta, cur_tr, NULL, NULL, &sc)) != eslOK) return status;
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
	ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), cm->align_opts CM_ALIGN_CHECKPARSESC, is on, but incompatible with another enabled option: CM_ALIGN_OPTACC.\n");
      if (fabs(sc - ParsetreeScore(cm, tr[i], cur_dsq, FALSE)) >= 0.01)
	ESL_FAIL(eslFAIL, errbuf, "DispatchAlignments(), seq: %d alignment score %.3f differs from its parse tree's score: %.3f", i, sc, ParsetreeScore(cm, tr[i], cur_dsq, FALSE));
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
	  ESL_FAIL(eslFAIL, errbuf, "DispatchAlignments(), Unable to convert sub CM parsetree to original CM parsetree. This shouldn't happen.");
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
  ESL_FAIL(eslEMEM, errbuf, "DispatchAlignments(), Memory allocation error.");
  return status; /* NEVERREACHED */
}

/***********************************************************
 * @LICENSE@
 ************************************************************/
/* cm_postprob.c
 * EPN 05.08.06
 * 
 * Functions for working with posterior probabilities for CMs.
 * Includes non-banded functions as well as banded ones (bands 
 * in the j and d dimensions). 
 *
 * Most of the alignment functions (ex: Inside) have two versions:
 * FInside() and IInside(), the former uses float log odds scores, the
 * latter uses scaled int log odds scores. Floats are more precise but
 * about 10X slower than ints, the difference is b/c it's necessary to go
 * into and out of log space to add floats in FLogsum, while ILogsum
 * uses a precomputed lookup table with ints.
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "easel.h"		/* general sequence analysis library    */
#include "esl_msa.h"   
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dirichlet.h"
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_fileparser.h"
#include "esl_gamma.h"
#include "esl_getopts.h"
#include "esl_gev.h"
#include "esl_gumbel.h"
#include "esl_histogram.h"
#include "esl_hyperexp.h"
#include "esl_keyhash.h"
#include "esl_minimizer.h"
#include "esl_mixgev.h"
/*#include "esl_mpi.h"*/
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msaweight.h"
#include "esl_normal.h"
#include "esl_paml.h"
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_regexp.h"
#include "esl_rootfinder.h"
#include "esl_scorematrix.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_stopwatch.h"
#include "esl_stretchexp.h"
/*#include "esl_swat.h"*/
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_weibull.h"
#include "esl_wuss.h"

#include "funcs.h"		/* external functions                   */
#include "old_funcs.h"		/* old external functions               */
#include "structs.h"		/* data structures, macros, #define's   */

static int   get_iemission_score(CM_t *cm, ESL_DSQ *dsq, int v, int i, int j);
static float get_femission_score(CM_t *cm, ESL_DSQ *dsq, int v, int i, int j);

/* 
 * Function: OldActuallyAlignTargets
 * Incept:   EPN, Tue Dec  5 14:25:02 2006
 *
 * Purpose:  Given a CM and sequences, do preliminaries, call the correct 
 *           alignment function and return parsetrees and optionally postal codes 
 *           (if cm->align_opts & CM_ALIGN_POST).
 *           Uses version 0.81 (old) DP functions, as opposed to the fast version
 *           1.0 DP functions called by ActuallyAlignTargets().
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
 *           Dies if there's an error (not good for MPI).
 */
int
OldActuallyAlignTargets(CM_t *cm, seqs_to_aln_t *seqs_to_aln, ESL_DSQ *dsq, search_results_t *search_results,
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
  CM_t              *sub_cm;       /* sub covariance model                      */
  CMSubMap_t        *submap;
  CP9Bands_t        *sub_cp9b;     /* data structure for hmm bands (bands on the hmm states) 
				    * and arrays for CM state bands, derived from HMM bands */
  CM_t              *orig_cm;      /* the original, template covariance model the sub CM was built from */
  int                spos;         /* HMM node most likely to have emitted posn 1 of target seq */
  int                spos_state;   /* HMM state type for curr spos 0=match or 1=insert */
  int                epos;         /* HMM node most likely to have emitted posn L of target seq */
  int                epos_state;   /* HMM state type for curr epos 0=match or  1=insert */
  Parsetree_t       *orig_tr;      /* parsetree for the orig_cm; created from the sub_cm parsetree */

  CP9_t             *sub_hmm;      /* constructed CP9 HMM */
  CP9Map_t          *sub_cp9map;   /* maps the sub_hmm to the sub_cm and vice versa */
  CP9_t             *orig_hmm;     /* original CP9 HMM built from orig_cm */
  CP9Map_t          *orig_cp9map;  /* original CP9 map */
  CP9Bands_t        *orig_cp9b;    /* original CP9Bands */

  /* variables related to query dependent banding (qdb) */
  int    expand_flag;           /* TRUE if the dmin and dmax vectors have just been 
				 * expanded (in which case we want to recalculate them 
				 * before we align a new sequence), and FALSE if not*/
  int *orig_dmin;               /* original dmin values passed in */
  int *orig_dmax;               /* original dmax values passed in */

  /* variables related to inside/outside */
  /*float           ***alpha;*/     /* alpha DP matrix for Inside() */
  /*float           ***beta; */     /* beta DP matrix for Inside() */
  /*float           ***post; */     /* post DP matrix for Inside() */
  int             ***alpha;    /* alpha DP matrix for Inside() */
  int             ***beta;     /* beta DP matrix for Inside() */
  int             ***post;     /* post DP matrix for Inside() */

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

  /* Contract check */
  if(!(cm->flags & CMH_BITS))                            cm_Fail("OldActuallyAlignTargets(), CMH_BITS flag down.\n");
  if(r == NULL && (cm->align_opts & CM_ALIGN_SAMPLE))    cm_Fail("OldActuallyAlignTargets(), no source of randomness, but CM_ALIGN_SAMPLE alignment option on.\n");
  if(r != NULL && (!(cm->align_opts & CM_ALIGN_SAMPLE))) cm_Fail("OldActuallyAlignTargets(), we have a source of randomness, but CM_ALIGN_SAMPLE alignment option off.\n");
  if((cm->align_opts & CM_ALIGN_POST)      && (cm->align_opts & CM_ALIGN_HMMVITERBI)) cm_Fail("OldActuallyAlignTargets(), CM_ALIGN_POST and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_HMMVITERBI)) cm_Fail("OldActuallyAlignTargets(), CM_ALIGN_SCOREONLY and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_POST))       cm_Fail("OldActuallyAlignTargets(), CM_ALIGN_SCOREONLY and CM_ALIGN_POST options are incompatible.\n");

  /* determine mode */
  if     (seqs_to_aln != NULL && (dsq == NULL && search_results == NULL))  sq_mode = TRUE;
  else if(seqs_to_aln == NULL && (dsq != NULL && search_results != NULL)) dsq_mode = TRUE;
  else   cm_Fail("OldActuallyAlignTargets(), can't determine mode (sq_mode or dsq_mode).\n");

  if( sq_mode && (seqs_to_aln->sq       == NULL)) cm_Fail("OldActuallyAlignTargets(), in sq_mode, seqs_to_aln->sq is NULL.\n");
  if( sq_mode && (seqs_to_aln->tr       != NULL)) cm_Fail("OldActuallyAlignTargets(), in sq_mode, seqs_to_aln->tr is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->cp9_tr   != NULL)) cm_Fail("OldActuallyAlignTargets(), in sq_mode, seqs_to_aln->cp9_tr is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->postcode1!= NULL)) cm_Fail("OldActuallyAlignTargets(), in sq_mode, seqs_to_aln->postcode1 is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->postcode2!= NULL)) cm_Fail("OldActuallyAlignTargets(), in sq_mode, seqs_to_aln->postcode2 is non-NULL.\n");
  if( sq_mode && (seqs_to_aln->sc       != NULL)) cm_Fail("OldActuallyAlignTargets(), in sq_mode, seqs_to_aln->sc is non-NULL.\n");
  
  if(dsq_mode && (cm->align_opts & CM_ALIGN_HMMVITERBI)) cm_Fail("OldActuallyAlignTargets(), in dsq_mode, CM_ALIGN_HMMVITERBI option on.\n");
  if(dsq_mode && (cm->align_opts & CM_ALIGN_POST))       cm_Fail("OldActuallyAlignTargets(), in dsq_mode, CM_ALIGN_POST option on.\n");
  if(dsq_mode && (cm->align_opts & CM_ALIGN_INSIDE))     cm_Fail("OldActuallyAlignTargets(), in dsq_mode, CM_ALIGN_INSIDE option on.\n");
  if(dsq_mode && (cm->align_opts & CM_ALIGN_SAMPLE))     cm_Fail("OldActuallyAlignTargets(), in dsq_mode, CM_ALIGN_SAMPLE option on.\n");
  if(dsq_mode && search_results == NULL)                 cm_Fail("OldActuallyAlignTargets(), in dsq_mode, search_results are NULL.\n");
  if(dsq_mode && (first_result > search_results->num_results)) cm_Fail("OldActuallyAlignTargets(), in dsq_mode, first_result: %d > search_results->num_results: %d\n", first_result, search_results->num_results);

  /* set the options based on cm->align_opts */
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

  /* another contract check */
  if((do_sample + do_inside + do_post) > 1) cm_Fail("OldActuallyAlignTargets(), exactly 0 or 1 of the following must be TRUE (== 1):\n\tdo_sample = %d\n\tdo_inside = %d\n\tdo_post%d\n\tdo_hmmonly: %d\n\tdo_scoreonly: %d\n", do_sample, do_inside, do_post, do_hmmonly, do_scoreonly);

  if(debug_level > 0) {
    printf("do_local  : %d\n", do_local);
    printf("do_qdb    : %d\n", do_qdb);
    printf("do_hbanded: %d\n", do_hbanded);
    printf("use_sums  : %d\n", use_sums);
    printf("do_sub    : %d\n", do_sub);
    printf("do_hmmonly: %d\n", do_hmmonly);
    printf("do_inside : %d\n", do_inside);
    printf("do_small  : %d\n", do_small);
    printf("do_post   : %d\n", do_post);
    printf("do_timings: %d\n", do_timings);
  }

  if      (sq_mode)   nalign = seqs_to_aln->nseq;
  else if(dsq_mode) { nalign = search_results->num_results - first_result; silent_mode = TRUE; }

  /* If sqmode: potentially allocate tr, cp9_tr, and postcodes. We'll set
   * seqs_to_aln->tr, seqs_to_aln->cp9_tr, seqs_to_aln->postcode1, and 
   * seqs_to_aln->postcode2 to these guys at end of function.
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
    if(do_post) {
      ESL_ALLOC(postcode1, sizeof(char **) * nalign);
      ESL_ALLOC(postcode2, sizeof(char **) * nalign);
    }
  }   
  ESL_ALLOC(parsesc, sizeof(float) * nalign);

  minsc = FLT_MAX;
  maxsc = -FLT_MAX;
  avgsc = 0;
  watch = esl_stopwatch_Create();

  if(do_hbanded || do_sub) /* We need a CP9 HMM to build sub_cms */
    {
      if(cm->cp9 == NULL)                    cm_Fail("OldActuallyAlignTargets, trying to use CP9 HMM that is NULL\n");
      if(cm->cp9b == NULL)                   cm_Fail("OldActuallyAlignTargets, cm->cp9b is NULL\n");
      if(!(cm->cp9->flags & CPLAN9_HASBITS)) cm_Fail("OldActuallyAlignTargets, trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");

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
      ip      = i + first_result; /* offset index in search_results structures */
      cur_dsq = dsq + search_results->data[ip].start - 1;
      cur_tr  = &(search_results->data[ip].tr);
      L       = search_results->data[ip].stop - search_results->data[ip].start + 1;
    }
    if (L == 0) continue; /* silently skip zero length seqs */

    /* Special case, if do_hmmonly, align seq with Viterbi, print score and move on to next seq */
    if(sq_mode && do_hmmonly) {
      if(sq_mode && !silent_mode) printf("Aligning (to a CP9 HMM w/viterbi) %-20s", seqs_to_aln->sq[i]->name);
      sc = CP9ViterbiAlign(cur_dsq, 1, L, cm->cp9, cm->cp9_mx, &(cp9_tr[i]));
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
      if((status = cp9_Seq2Bands(orig_cm, NULL, orig_cm->cp9_mx, orig_cm->cp9_bmx, cm->cp9_bmx, cur_dsq, 1, L, orig_cp9b, FALSE, debug_level)) != eslOK) cm_Fail("OldActuallyAlignTargets(), unrecoverable error in cp9_Seq2Bands().");
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
      if((status = cp9_Seq2Posteriors(orig_cm, NULL, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, cur_dsq, 1, L, debug_level)) != eslOK) cm_Fail("OldActuallyAlignTargets(), unrecoverable error in cp9_Seq2Posteriors().");
      
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
	cm_Fail("ERROR OldActuallyAlignTargets(), building sub CM.");
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
	if((status = cp9_Seq2Bands(sub_cm, NULL, sub_cm->cp9_mx, sub_cm->cp9_bmx, sub_cm->cp9_bmx, cur_dsq, 1, L, sub_cp9b, FALSE, debug_level)) != eslOK)  cm_Fail("OldActuallyAlignTargets(), unrecoverable error in cp9_Seq2Bands().");
	hmm           = sub_hmm;    
	cp9b          = sub_cp9b;
	cp9map        = sub_cp9map;
      }
    }
  
    /* Determine which CYK alignment algorithm to use, based
     * on command-line options AND memory requirements.
     */
    if(do_hbanded) {
      /* write a function to determine size of jd banded memory
       * req'd, and set do_small to true if its > thresh.
       if(do_small) * We're only going to band on d in memory, but 
       * we need to calculate safe_hd bands on the d dimension. 
       {
      */
    }
    if(do_qdb) {
      /*Check if we need to reset the query dependent bands b/c they're currently expanded. */
      if(expand_flag) {
	for(v = 0; v < cm->M; v++) {
	  cm->dmin[v] = orig_dmin[v];
	  cm->dmax[v] = orig_dmax[v];
	}
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
	IInside_b_jd_me(cm, cur_dsq, 1, L,
			TRUE,	   /* non memory-saving mode, we sample from alpha mx */
			NULL, &alpha,/* fill alpha, and return it, we'll sample a parsetree from it */
			NULL, NULL,  /* manage your own deckpool, I don't want it */
			do_local,    /* TRUE to allow local begins */
			cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	sc = ParsetreeSampleFromIInside_b_jd_me(r, cm, cur_dsq, L, alpha, cp9b, cur_tr, NULL);
      }
      else { /* sampling from inside matrix, but not HMM banded */
	IInside(cm, cur_dsq, 1, L,
		TRUE,            /* save full alpha, so we can sample from it,  */
		NULL, &alpha,    /* fill alpha, and return it, we'll sample a parsetree from it */
		NULL, NULL,      /* manage your own deckpool, I don't want it */
		do_local);       /* TRUE to allow local begins */
	sc = ParsetreeSampleFromIInside(r, cm, cur_dsq, L, alpha, cur_tr, NULL);
      }
    }
    else if(do_inside) { 
      if(do_hbanded) { /* HMM banded inside only */
	sc = IInside_b_jd_me(cm, cur_dsq, 1, L,
			     TRUE,	/* non memory-saving mode */
			     NULL, NULL,	/* manage your own matrix, I don't want it */
			     NULL, NULL,	/* manage your own deckpool, I don't want it */
			     do_local,    /* TRUE to allow local begins */
			     cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
      }
      else { /* non-banded inside only */
	sc = IInside(cm, cur_dsq, 1, L,
		     FALSE,       /* memory-saving mode */
		     NULL, NULL,	/* manage your own matrix, I don't want it */
		     NULL, NULL,	/* manage your own deckpool, I don't want it */
		     do_local);       /* TRUE to allow local begins */
      }
    }
    else if (do_small) { /* small D&C CYK alignment */
      if(do_qdb) { /* use QDBs when doing D&C CYK */
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, 
				 cur_tr, cm->dmin, cm->dmax);
	if(bdump_level > 0) qdb_trace_info_dump(cm, *cur_tr, cm->dmin, cm->dmax, bdump_level);
      }
      else if(do_hbanded) { /* use QDBs (safe d bands) derived from HMM bands when doing D&C CYK, HMM bands were not tight enough to allow HMM banded full CYK*/
	/* Calc the safe d bands */
	hd2safe_hd_bands(cm->M, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
			 cp9b->safe_hdmin, cp9b->safe_hdmax);
	if(debug_level > 3) { printf("\nprinting hd bands\n\n"); debug_print_hd_bands(cm, cp9b->hdmin, cp9b->hdmax, cp9b->jmin, cp9b->jmax); printf("\ndone printing hd bands\n\n"); }
	/* Note the following CYK call will not enforce j bands, even though user specified --hbanded. */
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, cp9b->safe_hdmin, cp9b->safe_hdmax);
	if(bdump_level > 0) qdb_trace_info_dump(cm, *cur_tr, cm->dmin, cm->dmax, bdump_level);
      }
      else { /* small D&C CYK non-banded alignment */
	/*printf("DEBUG PRINTING CM PARAMS BEFORE D&C CALL\n"); debug_print_cm_params(cm); printf("DONE DEBUG PRINTING CM PARAMS BEFORE D&C CALL\n");*/
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL); /* we're not in QDB mode */
	if(bdump_level > 0) { 
	  /* We want QDB info but QDBs weren't used.  Useful if you're curious why a QDB banded parse is crappy 
	   * relative to non-banded parse, e.g. allows you to see where the non-banded parse went outside the bands. */
	  qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
	}
      }
    }
    else if(do_qdb) { /* non-small, QDB banded CYK alignment */
      sc = CYKInside(cm, cur_dsq, L, 0, 1, L, cur_tr, cm->dmin, cm->dmax);
      if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
    }
    else if(do_hbanded) { /* non-small, HMM banded CYK alignment */
      sc = CYKInside_b_jd(cm, cur_dsq, L, 0, 1, L, cur_tr, cp9b->jmin, 
			  cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);
      /* if CM_ALIGN_HMMSAFE option is enabled, realign seqs w/HMM banded parses < 0 bits */
      if(cm->align_opts & CM_ALIGN_HMMSAFE && sc < 0.) { 
	tmpsc = sc;
	if(!silent_mode) printf("\n%s HMM banded parse had a negative score, realigning with non-banded CYK.\n", seqs_to_aln->sq[i]->name);
	FreeParsetree(*cur_tr);
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL); /* we're not in QDB mode */
	if(!silent_mode && fabs(sc-tmpsc) < 0.01) printf("HMM banded parse was the optimal parse.\n\n");
	else if (!silent_mode) printf("HMM banded parse was non-optimal, it was %.2f bits below the optimal.\n\n", (fabs(sc-tmpsc)));
      }	      
    }
    else { /* non-small, non-banded CYK alignment */
      sc = CYKInside(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL);
      if(bdump_level > 0) { 
	/* We want band info but --hbanded wasn't used.  Useful if you're curious why a banded parse is crappy 
	 * relative to non-banded parse, e.g. allows you to see where the non-banded parse went outside the bands. */
	qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
      }
    }
    /* end of large if() else if() else if() else statement */
  
    if(do_post) { /* do Inside() and Outside() runs and use alpha and beta to get posteriors */
      ESL_ALLOC(post, sizeof(int **) * (cm->M+1));
      if(do_hbanded) { /* HMM banded Inside/Outside --> posteriors */
	for (v = 0; v < cm->M; v++) post[v] = Ialloc_jdbanded_vjd_deck(L, 1, L, cp9b->jmin[v], cp9b->jmax[v], cp9b->hdmin[v], cp9b->hdmax[v]);
	post[cm->M] = NULL;
	post[cm->M] = Ialloc_vjd_deck(L, 1, L);
	sc = IInside_b_jd_me(cm, cur_dsq, 1, L,
			     TRUE,	/* save full alpha so we can run outside */
			     NULL, &alpha,	/* fill alpha, and return it, needed for IOutside() */
			     NULL, NULL,	/* manage your own deckpool, I don't want it */
			     do_local,       /* TRUE to allow local begins */
			     cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	sc = IOutside_b_jd_me(cm, cur_dsq, 1, L,
			      TRUE,	/* save full beta */
			      NULL, &beta,	/* fill beta, and return it, needed for ICMPosterior() */
			      NULL, NULL,	/* manage your own deckpool, I don't want it */
			      do_local,       /* TRUE to allow local begins */
			      alpha, &alpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
			      do_check,      /* TRUE to check Outside probs agree with Inside */
			      cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	ICMPosterior_b_jd_me(L, cm, alpha, NULL, beta, NULL, post, &post,
			     cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);
	ICMPostalCode_b_jd_me(cm, L, post, tr[i], cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
			      &(postcode1[i]), &(postcode2[i]));
      }
      else { /* non-HMM banded Inside/Outside --> posteriors */
	for (v = 0; v < cm->M+1; v++) post[v] = Ialloc_vjd_deck(L, 1, L);
	sc = IInside(cm, cur_dsq, 1, L,
		     TRUE,	/* save full alpha so we can run outside */
		     NULL, &alpha,	/* fill alpha, and return it, needed for IOutside() */
		     NULL, NULL,	/* manage your own deckpool, I don't want it */
		     do_local);       /* TRUE to allow local begins */
	sc = IOutside(cm, cur_dsq, 1, L,
		      TRUE,	/* save full beta */
		      NULL, &beta,	/* fill beta, and return it, needed for CMPosterior() */
		      NULL, NULL,	/* manage your own deckpool, I don't want it */
		      do_local,       /* TRUE to allow local begins */
		      alpha, &alpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
		      do_check);      /* TRUE to check Outside probs agree with Inside */
	ICMPosterior(L, cm, alpha, NULL, beta, NULL, post, &post); /* this frees alpha, beta */
	if(do_check) { 
	  ICMCheckPosterior(L, cm, post);
	  printf("\nPosteriors checked (I).\n\n");
	}
	ICMPostalCode(cm, L, post, tr[i], &(postcode1[i]), &(postcode2[i]));
      }
      /* free post  */
      if(post != NULL) Ifree_vjd_matrix(post, cm->M, 1, L);
    }
    /* done alignment for this seq */

    avgsc += sc;
    if (sc > maxsc) maxsc = sc;
    if (sc < minsc) minsc = sc;
      
    if(!silent_mode) printf("    score: %10.2f bits\n", sc);
    parsesc[i] = sc;
      
    /* check parsetree score if cm->align_opts & CM_ALIGN_CHECKPARSESC */
    if((cm->align_opts & CM_ALIGN_CHECKPARSESC) &&
       (!(cm->flags & CM_IS_SUB))) { 
      if (fabs(sc - ParsetreeScore(cm, tr[i], cur_dsq, FALSE)) >= 0.01)
	cm_Fail("ERROR in actually_align_target(), alignment score differs from its parse tree's score");
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
	  printf("\n\nIncorrectly converted original trace:\n");
	  ParsetreeDump(stdout, orig_tr, orig_cm, cur_dsq, NULL, NULL);
	  cm_Fail("this shouldn't happen.");
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
      FreeSubMap(submap);
      FreeCM(sub_cm); /* cm and sub_cm now point to NULL */
    }
    if(do_timings) { 
      esl_stopwatch_Stop(watch); 
      esl_stopwatch_Display(stdout, watch, "seq alignment CPU time: ");
      printf("\n");
    }
  }
  /* done aligning all nalign seqs. */
  /* Clean up. */
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
    free(parsesc);
  }

  return eslOK;
 ERROR:
  cm_Fail("Memory allocation error.");
  return status; /* NEVERREACHED */
}

/*****************************************************************
 * CM {F,I}Inside() & {F,I}Outside() functions.
 * Important: smallcyk.c has functions named Inside(), Outside() or some 
 *            variation, but these are implementations of variations of
 *            the CYK-Inside functions, as mentioned in the 2002 D&C
 *            BMC Bioinformatics publication.
 *****************************************************************/ 

/*
 * Function: FInside() 
 *           IInside() 
 * 
 * Purpose:  The Inside dynamic programming algorithm for CMs.
 *           Works directly with floats, stepping into and out of
 *           log space as necessary.
 * 
 *           Based on inside() in smallcyk.c, with the following 
 *           differences: necessarily align the sequence to the 
 *           full model (not possible to align to subtrees as in
 *           smallcyk.c's inside()), also no shadow matrix is
 *           kept because we're interested in ALL paths, finally
 *           we don't care about the best local begin state for the
 *           same reason.
 *  
 * Purpose:  Run the Inside() alignment algorithm, on a 
 *           subsequence from i0..j0, using the entire model.
 * 
 *           A note on the loop conventions. We're going to keep the
 *           sequence (sq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: i0-1..j in the rows, 
 *           and up to j0-i0+1 in the columns (d dimension). Where this is 
 *           handled the most is in two variables: L, which is the length of 
 *           the subsequence (j0-i0+1), and jp (read: j'), which is the 
 *           *relative* j w.r.t. the subsequence, ranging from 0..L, and 
 *           then d ranges from 0 to jp, and j is calculated from jp (i0-1+jp).
 *           
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in a dpool, 
 *           the decks *must* be sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the digitized sequence 
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           alpha     - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_alpha - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *                       
 *
 * Returns: log P(S|M)/P(S|R), as a bit score.
 */
float 
FInside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
       float ***alpha, float ****ret_alpha, 
       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
       int allow_begin)
{
  int      status;
  float  **end;         /* we re-use the end deck. */
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  float    bsc;		/* total score for using local begin states */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in FInside(), dsq is NULL.\n");

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */

  if (dpool == NULL) dpool = deckpool_create();
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, 0, E_st);
  for (jp = 0; jp <= L; jp++) {
    j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
    end[j][0] = 0.;
    for (d = 1; d <= jp; d++) end[j][d] = IMPOSSIBLE;
  }

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v <= (cm->M-1); v++) touch[v] = cm->pnum[v];

  /* EPN: now that the EL self loop has a transition score, its
   *      necessary to keep track of the alpha EL deck (alpha[cm->M])
   */
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      if (! deckpool_pop(dpool, &(alpha[cm->M]))) 
	alpha[cm->M] = alloc_vjd_deck(L, i0, j0);
      for (jp = 0; jp <= L; jp++) {
	j = i0-1+jp;
	alpha[cm->M][j][0] = 0.;
	/*alpha[cm->M][j][0] = IMPOSSIBLE;*/
	for (d = 1; d <= jp; d++)
	  {
	    alpha[cm->M][j][d] = (cm->el_selfsc * (d));
	  }
      }
    }

  /* Main recursion
   */
  for (v = (cm->M - 1); v >= 0; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      } 
      if (! deckpool_pop(dpool, &(alpha[v]))) 
	alpha[v] = alloc_vjd_deck(L, i0, j0);

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = FLogsum(alpha[v][j][d], (alpha[y][j][d] 
							      + cm->tsc[v][yoffset]));
		  }
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == B_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		z = cm->cnum[v];
		  
		alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0];
		for (k = 1; k <= d; k++)
		  alpha[v][j][d] = FLogsum(alpha[v][j][d], (alpha[y][j-k][d-k] 
							   + alpha[z][j][k]));
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    if (jp > 0) alpha[v][j][1] = IMPOSSIBLE;
	    for (d = 2; d <= jp; d++) 
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = FLogsum(alpha[v][j][d], (alpha[y][j-1][d-2] 
							      + cm->tsc[v][yoffset]));
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    for (d = 1; d <= jp; d++)
	      {
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = FLogsum(alpha[v][j][d], (alpha[y][j][d-1] 
							    + cm->tsc[v][yoffset]));
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][dsq[i]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    for (d = 1; d <= jp; d++)
	      {
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = FLogsum(alpha[v][j][d], (alpha[y][j-1][d-1] 
							      + cm->tsc[v][yoffset]));
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][dsq[j]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);		
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}				/* finished calculating deck v. */
      
      /* Keep track of contributions of local begins */
      if (allow_begin)
	{
	  bsc = FLogsum(bsc, (alpha[v][j0][L] + cm->beginsc[v]));
	}

      /* If we're at the root state, record contribution of local begins */
      if (allow_begin && v == 0)
	{
	  alpha[0][j0][L] = FLogsum(alpha[0][j0][L], bsc);
	}	  

      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { /* we can definitely release the S children of a bifurc. */
	    y = cm->cfirst[v]; deckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	    z = cm->cnum[v];   deckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) 
		  {
		    if (cm->sttype[y] == E_st) { 
		      nends--; 
		      if (nends == 0) { deckpool_push(dpool, end); end = NULL;}
		    } else 
		      deckpool_push(dpool, alpha[y]);
		    alpha[y] = NULL;
		  }
	      }
	  }
      }
  } /* end loop over all v */

  /*
    printf("INSIDE, printing alpha\n");
    debug_print_alpha(alpha, cm, L);
    printf("INSIDE, done printing alpha\n");
  */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  sc = alpha[0][j0][L];

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { deckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { deckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (deckpool_pop(dpool, &end)) free_vjd_deck(end, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  /*printf("\n\tFInside() sc  : %f\n", sc);*/
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* IInside() is the same as FInside, but uses scaled int log odds scores
 * instead of float log odds scores.
 */
float
IInside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
       int   ***alpha, int   ****ret_alpha, 
       Ideckpool_t *dpool, Ideckpool_t **ret_dpool,
       int allow_begin)
{
  int      status;
  int    **end;         /* we re-use the end deck. */
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    return_sc;   /* the return score, converted to bits (Scorified) */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  int      bsc;		/* total score for using local begin states */
  
  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in IInside(), dsq is NULL.\n");

  /* Allocations and initializations
   */
  bsc = -INFTY;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */

  if (dpool == NULL) dpool = Ideckpool_create();
  if (! Ideckpool_pop(dpool, &end))
    end = Ialloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, 0, E_st);
  for (jp = 0; jp <= L; jp++) {
    j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
    end[j][0] = 0.;
    for (d = 1; d <= jp; d++) end[j][d] = -INFTY;
  }

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(int   **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v <= (cm->M-1); v++) touch[v] = cm->pnum[v];

  /* EPN: now that the EL self loop has a transition score, its
   *      necessary to keep track of the alpha EL deck (alpha[cm->M])
   */
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      if (! Ideckpool_pop(dpool, &(alpha[cm->M]))) 
	alpha[cm->M] = Ialloc_vjd_deck(L, i0, j0);
      for (jp = 0; jp <= L; jp++) {
	j = i0-1+jp;
	alpha[cm->M][j][0] = Prob2Score(1.0, 1.0);
	/*alpha[cm->M][j][0] = IMPOSSIBLE;*/
	for (d = 1; d <= jp; d++)
	  {
	    alpha[cm->M][j][d] = (cm->iel_selfsc * (d));
	  }
      }
    }

  /* Main recursion
   */
  for (v = (cm->M - 1); v >= 0; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      } 
      if (! Ideckpool_pop(dpool, &(alpha[v]))) 
	alpha[v] = Ialloc_vjd_deck(L, i0, j0);

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = ILogsum(alpha[v][j][d], (alpha[y][j][d] 
							      + cm->itsc[v][yoffset]));
		  }
		if (alpha[v][j][d] < -INFTY) alpha[v][j][d] = -INFTY;
	      }
	  }
	}
      else if (cm->sttype[v] == B_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		z = cm->cnum[v];
		  
		alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0];
		for (k = 1; k <= d; k++)
		  alpha[v][j][d] = ILogsum(alpha[v][j][d], (alpha[y][j-k][d-k] 
							   + alpha[z][j][k]));
		if (alpha[v][j][d] < -INFTY) alpha[v][j][d] = -INFTY;
	      }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = -INFTY;
	    if (jp > 0) alpha[v][j][1] = -INFTY;
	    for (d = 2; d <= jp; d++) 
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = ILogsum(alpha[v][j][d], (alpha[y][j-1][d-2] 
							      + cm->itsc[v][yoffset]));
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->iesc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][d] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < -INFTY) alpha[v][j][d] = -INFTY;
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = -INFTY;
	    for (d = 1; d <= jp; d++)
	      {
		alpha[v][j][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = ILogsum(alpha[v][j][d], (alpha[y][j][d-1] 
							    + cm->itsc[v][yoffset]));
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][d] += cm->iesc[v][dsq[i]];
		else
		  alpha[v][j][d] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);		
		
		if (alpha[v][j][d] < -INFTY) alpha[v][j][d] = -INFTY;
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (jp = 0; jp <= L; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = -INFTY;
	    for (d = 1; d <= jp; d++)
	      {
		alpha[v][j][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = ILogsum(alpha[v][j][d], (alpha[y][j-1][d-1] 
							      + cm->itsc[v][yoffset]));
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->iesc[v][dsq[j]];
		else
		  alpha[v][j][d] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);		
		
		if (alpha[v][j][d] < -INFTY) alpha[v][j][d] = -INFTY;
	      }
	  }
	}				/* finished calculating deck v. */
      
      /* Keep track of contributions of local begins */
      if (allow_begin)
	{
	  bsc = ILogsum(bsc, (alpha[v][j0][L] + cm->ibeginsc[v]));
	}

      /* If we're at the root state, record contribution of local begins */
      if (allow_begin && v == 0)
	{
	  alpha[0][j0][L] = ILogsum(alpha[0][j0][L], bsc);
	}	  

      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { /* we can definitely release the S children of a bifurc. */
	    y = cm->cfirst[v]; Ideckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	    z = cm->cnum[v];   Ideckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) 
		  {
		    if (cm->sttype[y] == E_st) { 
		      nends--; 
		      if (nends == 0) { Ideckpool_push(dpool, end); end = NULL;}
		    } else 
		      Ideckpool_push(dpool, alpha[y]);
		    alpha[y] = NULL;
		  }
	      }
	  }
      }
  } /* end loop over all v */

  /*
    printf("INSIDE, printing alpha\n");
    debug_print_alpha(alpha, cm, L);
    printf("INSIDE, done printing alpha\n");
  */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  return_sc = Scorify(alpha[0][j0][L]);

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ideckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { Ideckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (Ideckpool_pop(dpool, &end)) Ifree_vjd_deck(end, i0, j0);
    Ideckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  ESL_DPRINTF1(("\tIInside() sc  : %f\n", return_sc));
  return return_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/***********************************************************************
 * Function: FOutside()
 * Date:     EPN 05.08.06
 * SRE, Tue Aug  8 10:42:52 2000 [St. Louis]
 *
 * Purpose:  The Outside dynamic programming algorithm for CMs.
 *           Works directly with floats, stepping into and out of 
 *           log space as necessary.
 *  
 *           Derived from smallcyk.c::CYKOutside() and smallcyk.c::outsdie(). 
 * 
 *           Align a subsequence to the full model, i.e. we're given
 *           i0 and j0, beginning and end positions of the subseq we're
 *           considering.
 * 
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha and beta, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           beta       - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_beta  - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           alpha     - the alpha matrix from a Inside run, if do_check is FALSE
 *                        only decks for S states must be valid, else all must be
 *                        valid.
 *           ret_alpha - if non-NULL, return the alpha matrix with one or more
 *                       decks available for examination (see "do_full")
 *           do_check  - TRUE to do time-consuming check to make sure
 *                       beta and alpha are consistent (only if NON-LOCAL mode)
 * 
 * Returns: log P(S|M)/P(S|R), as a bit score.
 ***********************************************************************/
float 
FOutside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
	 float ***beta, float ****ret_beta, 
	 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	 int allow_begin, float ***alpha, float ****ret_alpha, 
	 int do_check)
{
  int      status;
  int      v,y,z;	/* indices for states */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  float    escore;	/* an emission score, tmp variable */
  int      voffset;	/* index of v in t_v(y) transition scores */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  float    bsc;		/* total score for using local begin states */
  float    return_sc;   /* P(S|M)/P(S|R) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  float  **end;         /* we re-use the end deck. */
  int      fail_flag = FALSE; /* set to TRUE if do_check and we see a problem */

  /* Contract check */
  if(dsq == NULL) cm_Fail("ERROR in FOutside(), dsq is NULL.\n");
  if (cm->flags & CMH_LOCAL_END) { do_check = FALSE; } 
  /* Code for checking doesn't apply in local mode. See below. */

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* if caller didn't give us a matrix, make one.
   * Allocate room for M+1 decks because we might need the EL deck (M)
   * if we're doing local alignment.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck. Root is necessarily the ROOT_S state 0.
   */
  if (! deckpool_pop(dpool, &(beta[0])))
    beta[0] = alloc_vjd_deck(L, i0, j0);
  for (jp = 0; jp <= L; jp++) {
    j = i0-1+jp;
    for (d = 0; d <= jp; d++)
      beta[0][j][d] = IMPOSSIBLE;
  }
  beta[0][j0][L] = 0;		

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   */
  if (cm->flags & CMH_LOCAL_END) {
    if (! deckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = alloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= L; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[cm->M][j][d] = IMPOSSIBLE;
    }
    
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  /* Main loop down through the decks
   */
  /*for (v = 2; v < cm->M; v++) */ /*EPN is this 2 b/c Durbin p.287 
				     has state 2 in the algorithm? b/c state 1 is root*/
  for (v = 1; v < cm->M; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_vjd_deck(L, i0, j0);

      /* Init the whole deck to IMPOSSIBLE
       */
      for (jp = L; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  beta[v][j][d] = IMPOSSIBLE;
      }

      /* If we can do a local begin into v, also init with that. 
       * By definition, beta[0][j0][L] == 0.
       */ 
      if (i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN))
	beta[v][j0][L] = cm->beginsc[v];

      /* main recursion:
       */
      for (jp = L; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  {
	    if (cm->stid[v] == BEGL_S) 
	      {
		y = cm->plast[v];	/* the parent bifurcation    */
		z = cm->cnum[y];	/* the other (right) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j][0]; /* init on k=0 */
		for (k = 1; k <= L-j; k++)
		  beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+k][d+k] + alpha[z][j+k][k]));
	      }
	    else if (cm->stid[v] == BEGR_S) 
	      {
		y = cm->plast[v];	        /* the parent bifurcation    */
		z = cm->cfirst[y];	/* the other (left) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j-d][0];	/* init on k=0 */
		for (k = 1; k <= j-d; k++) 
		  beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j][d+k] + alpha[z][j-d][k]));
	      }
	    else
	      {
		i = j-d+1;
		for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
		  voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
		  switch(cm->sttype[y]) {
		  case MP_st: 
		    if (j == j0 || d == jp) continue; /* boundary condition */
		    
		    if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		      escore = cm->esc[y][(dsq[i-1]*cm->abc->K+dsq[j+1])];
		    else
		      escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
		    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+1][d+2] + cm->tsc[y][voffset]
							   + escore));
		    break;
		    
		  case ML_st:
		  case IL_st: 
		    if (d == jp) continue;	/* boundary condition (note when j=0, d=0)*/
		    
		    if (dsq[i-1] < cm->abc->K) 
		      escore = cm->esc[y][dsq[i-1]];
		    else
		      escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j][d+1] + cm->tsc[y][voffset] 
							   + escore));
		    break;
		    
		  case MR_st:
		  case IR_st:
		    if (j == j0) continue;
		    
		    if (dsq[j+1] < cm->abc->K) 
		      escore = cm->esc[y][dsq[j+1]];
		    else
		      escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);
		    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j+1][d+1] + cm->tsc[y][voffset] 
							   + escore));
		    break;
		    
		  case S_st:
		  case E_st:
		  case D_st:
		    beta[v][j][d] = FLogsum(beta[v][j][d], (beta[y][j][d] + cm->tsc[y][voffset])); 
		      break;
		    
		  default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
		  }/* end switch over states*/
		} /* ends for loop over parent states. we now know beta[v][j][d] for this d */
		if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;
	      }	/* ends else entered for non-BEGL/BEGR states*/	
	  } /* ends loop over d. We know all beta[v][j][d] in this row j*/
      }/* end loop over jp. We know the beta's for the whole deck.*/
	
      /* Deal with local alignment end transitions v->EL
       * (EL = deck at M.)
       */
      if (NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (jp = 0; jp <= L; jp++) { 
	  j = i0-1+jp;
	  for (d = 0; d <= jp; d++) 
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st: 
		if (j == j0 || d == jp) continue; /* boundary condition */
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+1][d+2] + cm->endsc[v] 
								+ escore));
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;	
		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[v][dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j][d+1] + cm->endsc[v] 
								+ escore));
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[v][dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j+1][d+1] + cm->endsc[v]
								+ escore));
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][j][d] + cm->endsc[v]));
		break;

	      case B_st:  
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
		/* note that although B is a valid vend for a segment we'd do
                   outside on, B->EL is set to be impossible, by the local alignment
                   config. There's no point in having a B->EL because B is a nonemitter
                   (indeed, it would introduce an alignment ambiguity). The same
		   alignment case is handled by the X->EL transition where X is the
		   parent consensus state (S, MP, ML, or MR) above the B. Thus,
		   this code is relying on the NOT_IMPOSSIBLE() test, above,
		   to make sure the sttype[vend]=B case gets into this switch.
		*/
	      } /* end switch over parent state type v */
	    } /* end inner loop over d */
	} /* end outer loop over jp */
      } /* end conditional section for dealing w/ v->EL local end transitions */

      /* Look at v's parents; if we're reusing memory (! do_full)
       * push the parents that we don't need any more into the pool.
       */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { deckpool_push(dpool, beta[y]); beta[y] = NULL; }
	}
      }
    } /* end loop over decks v. */

  /* EPN this code block is not superfluous for Inside() */
  /*     below are SRE's notes from a CYK inside() function */
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = L; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->el_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE, printing beta\n");
    debug_print_alpha(beta, cm, L);
    printf("OUTSIDE, done printing beta\n");
  */

  if(do_check && (!(cm->flags & CMH_LOCAL_END))) 
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    {
      /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
       * using both the Outside (beta) and Inside (alpha) matrices,
       * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
       * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
       *                                                = P(S|M) / P(S|R)  
       * in v,j,d coordinates this is:
       * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
       *                                                = P(S|M) / P(S|R)
       */
	 
      for(n = 0; n < cm->nodes; n++)
	{
	  sc = IMPOSSIBLE;
	  num_split_states = SplitStatesInNode(cm->ndtype[n]);
	  for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++)
	    {
	      for (jp = 0; jp <= L; jp++) 
		{ 
		  j = i0-1+jp;
		  for (d = 0; d <= jp; d++) 
		    {
		      sc = FLogsum(sc, (alpha[v][j][d] + beta[v][j][d]));
		    }
		}
	      }
	  diff = alpha[0][j0][L] - sc;
	  if(diff > 0.001 || diff < -0.001)
	    {
	      fail_flag = TRUE;
	      printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
		     n, sc, alpha[0][j0][L], diff);
	    }
	}

    }

  /* IF not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state with 
   * d = 0.
   * We pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * IF local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */

  if(!(cm->flags & CMH_LOCAL_END))
    {
      return_sc = IMPOSSIBLE;
      for (jp = 0; jp <= L; jp++) 
	{ 
	  j = i0-1+jp;
	  /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	  return_sc = FLogsum(return_sc, (beta[cm->M-1][j][0]));
	}
    }
  else /* return_sc = P(S|M) / P(S|R) from Inside() */
    {
      return_sc = alpha[0][j0][L];
    }
  /* If the caller doesn't want the beta matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
      deckpool_push(dpool, beta[cm->M]);
      beta[cm->M] = NULL; 
    }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the alpha matrix, free it 
   * Else, pass it back to him.
   * EPN - if we free the alpha and beta matrix the deck pool has all the 
   *       decks from alpha and beta, not sure if this is desirable.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { deckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { deckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    float **a;
    while (deckpool_pop(dpool, &a)) free_vjd_deck(a, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }
  free(touch);

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tFOutside() sc : %f\n", return_sc));
  else
    ESL_DPRINTF1(("\tFOutside() sc : %f (LOCAL mode; sc is from Inside)\n", return_sc));

  return return_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* IOutside() is the same as FOutside, but uses scaled int log odds scores
 * instead of float log odds scores.
 */
float
IOutside(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
	 int   ***beta, int   ****ret_beta, 
	 Ideckpool_t *dpool, Ideckpool_t **ret_dpool,
	 int allow_begin, int   ***alpha, int   ****ret_alpha, 
	 int do_check)
{
  int      status;
  int      v,y,z;	/* indices for states */
  int      j,d,i,k;	/* indices in sequence dimensions */
  int      isc;     	/* a temporary variable holding an int score */
  float    fsc;     	/* a temporary variable holding a float score */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      iescore;	/* an emission score, tmp variable */
  int      voffset;	/* index of v in t_v(y) transition scores */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  int      bsc;		/* total score for using local begin states */
  int      ireturn_sc;  /* P(S|M)/P(S|R), a scaled int*/
  float    freturn_sc;  /* P(S|M)/P(S|R), a float (Scorified ireturn_sc) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  int    **end;         /* we re-use the end deck. */
  int      fail_flag = FALSE; /* set to TRUE if do_check and we see a problem */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in IOutside(), dsq is NULL.\n");
  if (cm->flags & CMH_LOCAL_END) { do_check = FALSE; } 
  /* Code for checking doesn't apply in local mode. See below. */

  /* Allocations and initializations
   */
  bsc = -INFTY;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = Ideckpool_create();

  /* if caller didn't give us a matrix, make one.
   * Allocate room for M+1 decks because we might need the EL deck (M)
   * if we're doing local alignment.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(int   **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck. Root is necessarily the ROOT_S state 0.
   */
  if (! Ideckpool_pop(dpool, &(beta[0])))
    beta[0] = Ialloc_vjd_deck(L, i0, j0);
  for (jp = 0; jp <= L; jp++) {
    j = i0-1+jp;
    for (d = 0; d <= jp; d++)
      beta[0][j][d] = -INFTY;
  }
  beta[0][j0][L] = 0;		

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   */
  if (cm->flags & CMH_LOCAL_END) {
    if (! Ideckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = Ialloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= L; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[cm->M][j][d] = -INFTY;
    }
    
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  /* Main loop down through the decks
   */
  /*for (v = 2; v < cm->M; v++) */ /*EPN is this 2 b/c Durbin p.287 
				     has state 2 in the algorithm? b/c state 1 is root*/
  for (v = 1; v < cm->M; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! Ideckpool_pop(dpool, &(beta[v])))
	beta[v] = Ialloc_vjd_deck(L, i0, j0);

      /* Init the whole deck to -INFTY
       */
      for (jp = L; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  beta[v][j][d] = -INFTY;
      }

      /* If we can do a local begin into v, also init with that. 
       * By definition, beta[0][j0][L] == 0.
       */ 
      if (i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN))
	beta[v][j0][L] = cm->ibeginsc[v];

      /* main recursion:
       */
      for (jp = L; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  {
	    if (cm->stid[v] == BEGL_S) 
	      {
		y = cm->plast[v];	/* the parent bifurcation    */
		z = cm->cnum[y];	/* the other (right) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j][0]; /* init on k=0 */
		for (k = 1; k <= L-j; k++)
		  beta[v][j][d] = ILogsum(beta[v][j][d], (beta[y][j+k][d+k] + alpha[z][j+k][k]));
	      }
	    else if (cm->stid[v] == BEGR_S) 
	      {
		y = cm->plast[v];	        /* the parent bifurcation    */
		z = cm->cfirst[y];	/* the other (left) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j-d][0];	/* init on k=0 */
		for (k = 1; k <= j-d; k++) 
		  beta[v][j][d] = ILogsum(beta[v][j][d], (beta[y][j][d+k] + alpha[z][j-d][k]));
	      }
	    else
	      {
		i = j-d+1;
		for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
		  voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
		  switch(cm->sttype[y]) {
		  case MP_st: 
		    if (j == j0 || d == jp) continue; /* boundary condition */
		    
		    if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		      iescore = cm->iesc[y][(dsq[i-1]*cm->abc->K+dsq[j+1])];
		    else
		      iescore = iDegeneratePairScore(cm->abc, cm->iesc[y], dsq[i-1], dsq[j+1]);
		    beta[v][j][d] = ILogsum(beta[v][j][d], (beta[y][j+1][d+2] + cm->itsc[y][voffset]
							   + iescore));
		    break;
		    
		  case ML_st:
		  case IL_st: 
		    if (d == jp) continue;	/* boundary condition (note when j=0, d=0)*/
		    
		    if (dsq[i-1] < cm->abc->K) 
		      iescore = cm->iesc[y][dsq[i-1]];
		    else
		      iescore = esl_abc_IAvgScore(cm->abc, dsq[i-1], cm->iesc[y]);
		    beta[v][j][d] = ILogsum(beta[v][j][d], (beta[y][j][d+1] + cm->itsc[y][voffset] 
							   + iescore));
		    break;
		    
		  case MR_st:
		  case IR_st:
		    if (j == j0) continue;
		    
		    if (dsq[j+1] < cm->abc->K) 
		      iescore = cm->iesc[y][dsq[j+1]];
		    else
		      iescore = esl_abc_IAvgScore(cm->abc, dsq[j+1], cm->iesc[y]);
		    beta[v][j][d] = ILogsum(beta[v][j][d], (beta[y][j+1][d+1] + cm->itsc[y][voffset] 
							   + iescore));
		    break;
		    
		  case S_st:
		  case E_st:
		  case D_st:
		    beta[v][j][d] = ILogsum(beta[v][j][d], (beta[y][j][d] + cm->itsc[y][voffset])); 
		      break;
		    
		  default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
		  }/* end switch over states*/
		} /* ends for loop over parent states. we now know beta[v][j][d] for this d */
		if (beta[v][j][d] < -INFTY) beta[v][j][d] = -INFTY;
	      }	/* ends else entered for non-BEGL/BEGR states*/	
	  } /* ends loop over d. We know all beta[v][j][d] in this row j*/
      }/* end loop over jp. We know the beta's for the whole deck.*/
	
      /* Deal with local alignment end transitions v->EL
       * (EL = deck at M.)
       */
      if (cm->iendsc[v] != -INFTY) {
	for (jp = 0; jp <= L; jp++) { 
	  j = i0-1+jp;
	  for (d = 0; d <= jp; d++) 
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st: 
		if (j == j0 || d == jp) continue; /* boundary condition */
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  iescore = cm->iesc[v][(dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  iescore = iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i-1], dsq[j+1]);
		beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][j+1][d+2] + cm->iendsc[v] 
								+ iescore));
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;	
		if (dsq[i-1] < cm->abc->K) 
		  iescore = cm->iesc[v][dsq[i-1]];
		else
		  iescore = esl_abc_IAvgScore(cm->abc, dsq[i-1], cm->iesc[v]);
		beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][j][d+1] + cm->iendsc[v] 
								+ iescore));
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < cm->abc->K) 
		  iescore = cm->iesc[v][dsq[j+1]];
		else
		  iescore = esl_abc_IAvgScore(cm->abc, dsq[j+1], cm->iesc[v]);
		beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][j+1][d+1] + cm->iendsc[v]
								+ iescore));
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][j][d] + cm->iendsc[v]));
		break;

	      case B_st:  
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
		/* note that although B is a valid vend for a segment we'd do
                   outside on, B->EL is set to be impossible, by the local alignment
                   config. There's no point in having a B->EL because B is a nonemitter
                   (indeed, it would introduce an alignment ambiguity). The same
		   alignment case is handled by the X->EL transition where X is the
		   parent consensus state (S, MP, ML, or MR) above the B. Thus,
		   this code is relying on the (cm->iendsc[v] != -INFTY) test, above,
		   to make sure the sttype[vend]=B case gets into this switch.
		*/
	      } /* end switch over parent state type v */
	    } /* end inner loop over d */
	} /* end outer loop over jp */
      } /* end conditional section for dealing w/ v->EL local end transitions */

      /* Look at v's parents; if we're reusing memory (! do_full)
       * push the parents that we don't need any more into the pool.
       */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { Ideckpool_push(dpool, beta[y]); beta[y] = NULL; }
	}
      }
    } /* end loop over decks v. */

  /* EPN this code block is not superfluous for Inside() */
  /*     below are SRE's notes from a CYK inside() function */
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = L; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->iel_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE, printing beta\n");
    debug_print_alpha(beta, cm, L);
    printf("OUTSIDE, done printing beta\n");
  */

  /* EPN, Mon Apr 23 16:34:37 2007 
   * Found problems with the following test (I think they are due to
   * precision) when testing prior to release 0.8. I'm taking it
   * out. It's redundant anyway, the cmalign --checkpost option calls
   * the function CMCheckPosterior() which checks that Inside and
   * Outside matrices agree, albeit in a different way. I think this
   * test is valid, it's the precision issue with int log odds scores,
   * though I can't be sure right now.
   */
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) 
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    {
      /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
       * using both the Outside (beta) and Inside (alpha) matrices,
       * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
       * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
       *                                                = P(S|M) / P(S|R)  
       * in v,j,d coordinates this is:
       * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
       *                                                = P(S|M) / P(S|R)
       */

      for(n = 0; n < cm->nodes; n++)
	{
	  isc = -INFTY;
	  num_split_states = SplitStatesInNode(cm->ndtype[n]);
	  for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++)
	    {
	      for (jp = 0; jp <= L; jp++) 
		{ 
		  j = i0-1+jp;
		  for (d = 0; d <= jp; d++) 
		    isc = ILogsum(isc, (alpha[v][j][d] + beta[v][j][d]));
		}
	      }
	  fsc = Scorify(isc);
	  diff = fsc - (Scorify(alpha[0][j0][L]));
	  if(diff > 0.01 || diff < -0.01)
	    {
	      fail_flag = TRUE;
	      printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
		     n, fsc, Scorify(alpha[0][j0][L]), diff);
	    }
	}

    }

  /* IF not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state with 
   * d = 0.
   * We pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * IF local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */

  if(!(cm->flags & CMH_LOCAL_END))
    {
      ireturn_sc = -INFTY;
      for (jp = 0; jp <= L; jp++) 
	{ 
	  j = i0-1+jp;
	  /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	  ireturn_sc = ILogsum(ireturn_sc, (beta[cm->M-1][j][0]));
	}
      freturn_sc = Scorify(ireturn_sc);
    }
  else /* ireturn_sc = P(S|M) / P(S|R) from Inside() */
    {
      freturn_sc = Scorify(alpha[0][j0][L]);
    }
  /* If the caller doesn't want the beta matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (beta[v] != NULL) { Ideckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
      Ideckpool_push(dpool, beta[cm->M]);
      beta[cm->M] = NULL; 
    }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the alpha matrix, free it 
   * Else, pass it back to him.
   * EPN - if we free the alpha and beta matrix the deck pool has all the 
   *       decks from alpha and beta, not sure if this is desirable.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ideckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { Ideckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    int   **a;
    while (Ideckpool_pop(dpool, &a)) Ifree_vjd_deck(a, i0, j0);
    Ideckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }
  free(touch);

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tIOutside() sc : %f\n", freturn_sc));
  else
    ESL_DPRINTF1(("\tIOutside() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  return freturn_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/***************************************************************/
/* Function: FCMPosterior() 
 *           
 * EPN 05.25.06 based on IHH's P7EmitterPosterior() from HMMER's postprob.c
 *
 * Purpose:  Combines Inside and Outside cubes into a posterior
 *           probability cube.
 *           The entry in post[v][j][d] is the log of the
 *           posterior probability of a parse subtree rooted at v 
 *           emitting the subsequence i..j (i=j-d+1).
 *           The caller must allocate space for the cube, although the
 *           beta matrix from Outside can be used instead (overwriting it will not
 *           compromise the algorithm).
 *           
 * Args:     L        - length of sequence
 *           cm       - the model
 *           alpha    - pre-calculated Inside matrix 
 *           ret_alpha - if non-NULL, return the matrix as it was passed in,
 *                       else free it.
 *           beta     - pre-calculated Outside matrix
 *           ret_beta - if non-NULL, return the matrix as it was passed in,
 *                       else free it.
 *           post     - pre-allocated dynamic programming cube
 *           ret_post - if non-NULL, return the posterior matrix,
 *                       else free it.
 *           
 * Return:   void
 */
void
FCMPosterior(int L, CM_t *cm, float ***alpha, float ****ret_alpha, float ***beta, float ****ret_beta, 
	     float ***post, float ****ret_post)
{
  int   v, j, d;
  float sc;
  float  **end;         /* used for freeing alpha b/c we re-use the end deck. */
  int vmax;
  sc = alpha[0][L][L];
  
  /* If local ends are on, start with the EL state (cm->M), otherwise
   * its not a valid deck.
   */
  vmax = cm->M-1;
  if (cm->flags & CMH_LOCAL_END) vmax = cm->M;

  for (v = vmax; v >= 0; v--) 
    for (j = 0; j <= L; j++) 
      for (d = 0; d <= j; d++)
	{
	  post[v][j][d] = alpha[v][j][d] + beta[v][j][d] - sc;
	  if(v == vmax)
	    {
	      /*printf("v: %3d | j: %3d | d: %3d | alpha: %5.2f | beta: %5.2f\n", v, j, d, alpha[v][j][d], beta[v][j][d]);*/
	      /*printf("post[%d][%d][%d]: %f\n", cm->M, j, d, post[cm->M][j][d]);*/
	    }
	}
  /* If the caller doesn't want the matrix, free it and free the decks in the pool
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= cm->M; v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { free_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { free_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  }
  else *ret_alpha = alpha;

  /* If the caller doesn't want the beta matrix, free it along with the decks.
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= cm->M; v++) 
      if (beta[v] != NULL) { free_vjd_deck(beta[v], 1, L); beta[v] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the post matrix, free it, though
   * it would be *stupid* for the caller not to want it in current implementation.
   */
  if (ret_post == NULL) {
    for (v = 0; v <= cm->M; v++) 
      if (post[v] != NULL) { free_vjd_deck(post[v], 1, L); post[v] = NULL; }
    free(post);
  } else *ret_post = post;

}

/* ICMPosterior() is the same as CMPosterior, but uses scaled int log odds scores
 * instead of float log odds scores.
 */
void
ICMPosterior(int L, CM_t *cm, int   ***alpha, int   ****ret_alpha, int   ***beta, int   ****ret_beta, 
	    int   ***post, int   ****ret_post)
{
  int   v, j, d;
  int   sc;
  int    **end;         /* used for freeing alpha b/c we re-use the end deck. */
  int vmax;

  /*printf("in ICMPosterior()\n");*/
  sc = alpha[0][L][L];
  
  /* If local ends are on, start with the EL state (cm->M), otherwise
   * its not a valid deck.
   */
  vmax = cm->M-1;
  if (cm->flags & CMH_LOCAL_END) vmax = cm->M;

  for (v = vmax; v >= 0; v--) 
    for (j = 0; j <= L; j++) 
      for (d = 0; d <= j; d++)
	{
	  /* printf("v: %2d | j: %2d | d: %2d | alpha[%d][%d][%d]: %d | beta[%d][%d][%d]: %d\n",  v, j, d, v, j, d, alpha[v][j][d], v,j,d, beta[v][j][d]);  */
	  post[v][j][d] = alpha[v][j][d] + beta[v][j][d] - sc;
	  /*printf("v: %3d | j: %3d | d: %3d | alpha: %10d | beta: %10d | post: %10d\n", v, j, d, alpha[v][j][d], beta[v][j][d], post[v][j][d]);*/
	  if(v == vmax)
	    {
	      /*printf("v: %3d | j: %3d | d: %3d | alpha: %5.2f | beta: %5.2f\n", v, j, d, alpha[v][j][d], beta[v][j][d]);*/
	      /*printf("post[%d][%d][%d]: %f\n", cm->M, j, d, post[cm->M][j][d]);*/
	    }
	}
  /* If the caller doesn't want the matrix, free it and free the decks in the pool
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= cm->M; v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ifree_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { Ifree_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  }
  else *ret_alpha = alpha;

  /* If the caller doesn't want the beta matrix, free it along with the decks.
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= cm->M; v++) 
      if (beta[v] != NULL) { Ifree_vjd_deck(beta[v], 1, L); beta[v] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the post matrix, free it, though
   * it would be *stupid* for the caller not to want it in current implementation.
   */
  if (ret_post == NULL) {
    for (v = 0; v <= cm->M; v++) 
      if (post[v] != NULL) { Ifree_vjd_deck(post[v], 1, L); post[v] = NULL; }
    free(post);
  } else *ret_post = post;

}


/* ICMPostalCode() is the same as CMPostalCode, but uses scaled int log odds scores
 * instead of float log odds scores.
 */
void
ICMPostalCode(CM_t *cm, int L, int ***post, Parsetree_t *tr, char **ret_pcode1, char **ret_pcode2)
{
  int status;
  int x, v, i, j, d, r, p;
  char *pcode1;
  char *pcode2;

  ESL_ALLOC(pcode1, (L+1) * sizeof(char)); 
  ESL_ALLOC(pcode2, (L+1) * sizeof(char)); 

  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];
      d = j-i+1;
      /* printf("x: %2d | v: %2d | i: %2d | j: %2d | d: %2d | post[%d][%d][%d]: %d\n", x, v, i, j, d, v, j, d, post[v][j][d]);*/
      /*
       * Only P, L, R states have emissions.
       */
      if (cm->sttype[v] == MP_st) {
	p = Iscore2postcode(post[v][j][d]);
	if(p == 100) { 
	  pcode1[i-1] = pcode1[j-1] = '*';
	  pcode2[i-1] = pcode2[j-1] = '*';
	}
	else {
	  pcode1[i-1] = pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	p = Iscore2postcode(post[v][j][d]);
	if(p == 100) { 
	  pcode1[i-1] = '*';
	  pcode2[i-1] = '*';
	}
	else {
	  pcode1[i-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	p = Iscore2postcode(post[v][j][d]);
	if(p == 100) { 
	  pcode1[j-1] = '*';
	  pcode2[j-1] = '*';
	}
	else {
	  pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == EL_st) /*special case*/ {
	for(r = (i-1); r <= (j-1); r++)
	  {
	    d = j - (r+1) + 1;
	    p = Iscore2postcode(post[v][j][d]);
	    if(p == 100) { 
	      pcode1[r] = '*';
	      pcode2[r] = '*';
	    }
	    else {
	      pcode1[r] = '0' + (char) (p / 10);
	      pcode2[r] = '0' + (char) (p % 10);
	    }
	    /*printf("r: %d | post[%d][%d][%d]: %f | sc: %c\n", r, v, j, d, post[v][j][d], postcode[r]);*/
	  }
      }
    }
  pcode1[L] = '\0';
  pcode2[L] = '\0';
  *ret_pcode1 = pcode1;
  *ret_pcode2 = pcode2;
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
  return; /* never reached */
}


void
CMPostalCode_b_jd_me(CM_t *cm, int L, float ***post, Parsetree_t *tr,
		    int *jmin, int *jmax, int **hdmin, int **hdmax, char **ret_pcode1, char **ret_pcode2)
{
  int status;
  int x, v, i, j, d, r, p;
  char *pcode1;
  char *pcode2;
  int jp_v, dp_v;
  ESL_ALLOC(pcode1, (L+1) * sizeof(char)); 
  ESL_ALLOC(pcode2, (L+1) * sizeof(char)); 

  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];
      d = j-i+1;
      /*
       * Only P, L, R states have emissions.
       */
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];

      dp_v = d - hdmin[v][jp_v];
      if (cm->sttype[v] == MP_st) {
	p = Fscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[i-1] = pcode1[j-1] = '*';
	  pcode2[i-1] = pcode2[j-1] = '*';
	}
	else {
	  pcode1[i-1] = pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	p = Fscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[i-1] = '*';
	  pcode2[i-1] = '*';
	}
	else {
	  pcode1[i-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	p = Fscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[j-1] = '*';
	  pcode2[j-1] = '*';
	}
	else {
	  pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == EL_st) /*special case*/ {
	for(r = (i-1); r <= (j-1); r++) {
	  d = j - (r+1) + 1;
	  p = Fscore2postcode(post[v][j][d]);
	    if(p == 100) { 
	      pcode1[r] = '*';
	      pcode2[r] = '*';
	    }
	    else {
	      pcode1[r] = '0' + (char) (p / 10);
	      pcode2[r] = '0' + (char) (p % 10);
	    }
	    /*printf("r: %d | post[%d][%d][%d]: %f | sc: %c\n", r, v, j, d, post[v][j][d], postcode[r]);*/
	  }
      }
    }
  pcode1[L] = '\0';
  pcode2[L] = '\0';
  *ret_pcode1 = pcode1;
  *ret_pcode2 = pcode2;
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
  return; /* never reached */
}


int 
Iscore2postcode(int sc)
{
  int i;
  i = (int) (Score2Prob(sc, 1.) * 100.);
  ESL_DASSERT1((i >= 0 && i <= 100)); 
  return i;
}

void
ICMPostalCode_b_jd_me(CM_t *cm, int L, int ***post, Parsetree_t *tr,
		      int *jmin, int *jmax, int **hdmin, int **hdmax, char **ret_pcode1, char **ret_pcode2)
{
  int status;
  int x, v, i, j, d, r, p;
  char *pcode1;
  char *pcode2;
  int jp_v, dp_v;

  ESL_ALLOC(pcode1, (L+1) * sizeof(char)); 
  ESL_ALLOC(pcode2, (L+1) * sizeof(char)); 

  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];
      d = j-i+1;
      /*
       * Only P, L, R states have emissions.
       */
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
      if (cm->sttype[v] == MP_st) {
	p = Iscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[i-1] = pcode1[j-1] = '*';
	  pcode2[i-1] = pcode2[j-1] = '*';
	}
	else {
	  pcode1[i-1] = pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	p = Iscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[i-1] = '*';
	  pcode2[i-1] = '*';
	}
	else {
	  pcode1[i-1] = '0' + (char) (p / 10);
	  pcode2[i-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	p = Iscore2postcode(post[v][jp_v][dp_v]);
	if(p == 100) { 
	  pcode1[j-1] = '*';
	  pcode2[j-1] = '*';
	}
	else {
	  pcode1[j-1] = '0' + (char) (p / 10);
	  pcode2[j-1] = '0' + (char) (p % 10);
	}
      } else if (cm->sttype[v] == EL_st) /*special case*/ {
	for(r = (i-1); r <= (j-1); r++) {
	  d = j - (r+1) + 1;
	  p = Iscore2postcode(post[v][j][d]);
	    if(p == 100) { 
	      pcode1[r] = '*';
	      pcode2[r] = '*';
	    }
	    else {
	      pcode1[r] = '0' + (char) (p / 10);
	      pcode2[r] = '0' + (char) (p % 10);
	    }
	    /*printf("r: %d | post[%d][%d][%d]: %f | sc: %c\n", r, v, j, d, post[v][j][d], postcode[r]);*/
	  }
      }
    }
  pcode1[L] = '\0';
  pcode2[L] = '\0';
  *ret_pcode1 = pcode1;
  *ret_pcode2 = pcode2;
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
  return; /* never reached */
}

/* ICMCheckPosterior() is the same as CMCheckPosterior(), but uses scaled int log odds scores
 * instead of float log odds scores.
 */
void
ICMCheckPosterior(int L, CM_t *cm, int ***post)
{
  int   v, j, d, k;
  int   sc;
  float fsc;
  for (k = 1; k <= L; k++) 
    {
      sc = -INFTY;
      for (v = (cm->M - 1); v >= 0; v--) 
	{
	  if((cm->sttype[v] == MP_st) ||
	     (cm->sttype[v] == ML_st) ||
	     (cm->sttype[v] == IL_st))
	    {
	      for (j = k; j <= L; j++)
		{
		  d = j-k+1;
		  sc = ILogsum(sc, (post[v][j][d]));
		  /* printf("adding L v: %d | i: %d | j: %d | d: %d | sc: %10d\n", v, (j-d+1), j, d, sc); */
		}
	    }
	  if((cm->sttype[v] == MP_st) ||
	     (cm->sttype[v] == MR_st) ||
	     (cm->sttype[v] == IR_st))
	    {
	      for (d = 1; d <= k; d++)
		{
		  sc = ILogsum(sc, (post[v][k][d]));
		  /* printf("adding R v: %d | i: %d | j: %d | d: %d | sc: %10d\n", v, (k-d+1), k, d, sc); */
		}
	    }
	}
      /* Finally factor in possibility of a local end, i.e. that the EL state
       * may have "emitted" this residue.
       */
      if (cm->flags & CMH_LOCAL_END) {
	for (j = k; j <= L; j++)
	  {
	    d = j-k+1;
	    /*printf("EL adding L v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", cm->M, (j-d+1), j, d, post[cm->M][j][d]);*/
	    sc = ILogsum(sc, (post[cm->M][j][d]));
	  }
      }

      fsc = Scorify(sc);
      /* printf("k: %4d sc: %10.7f\n", k, fsc); */
      if(((fsc - 0.) > 0.01) || ((fsc - 0.) < -0.01))
	{
	  cm_Fail("residue position %d has summed prob of %5.4f (2^%5.4f) in posterior cube.\n", k, (sreEXP2(fsc)), fsc);
	}
      /*printf("k: %d | total: %10.2f\n", k, (sreEXP2(sc)));*/
    }  
  ESL_DPRINTF1(("ICMCheckPosterior() passed, all residues have summed probability of emission of 1.0.\n"));
}


/*****************************************************************
 * CM {F,I}Inside_b_jd_me() & {F,I}Outside_b_jd_me() functions.
 * Banded versions of {F,I}Inside() and {F,I}Outside() that only 
 * allocate cells within state dependent j bands, and 
 * state & j dependent d bands.
 *****************************************************************/ 

/* Function: FInside_b_jd_me()
 * EPN 05.26.06
 * 
 * Purpose:  The banded Inside dynamic programming algorithm for CMs.
 *           Banded in both the j and d dimensions, and works
 *           with transformed coordinates for memory efficiency, 
 *           only alpha cells within bands are allocated.
 *           Works directly with floats, stepping into and out of
 *           log space as necessary.
 * 
 *           Based on inside() in smallcyk.c, with the following 
 *           differences: necessarily align the sequence to the 
 *           full model (not possible to align to subtrees as in
 *           smallcyk.c's inside()), also no shadow matrix is
 *           kept because we're interested in ALL paths, finally
 *           we don't care about the best local begin state for the
 *           same reason.
 *  
 *           Run the Inside() alignment algorithm, on a 
 *           subsequence from i0..j0, using the entire model, enforcing
 *           bands in both the j and d dimensions.
 * 
 *           A note on the loop conventions. We're going to keep the
 *           sequence (dsq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: i0-1..j in the rows, 
 *           and up to j0-i0+1 in the columns (d dimension). Where this is 
 *           handled the most is in two variables: L, which is the length of 
 *           the subsequence (j0-i0+1), and jp (read: j'), which is the 
 *           *relative* j w.r.t. the subsequence, ranging from 0..W, and 
 *           then d ranges from 0 to jp, and j is calculated from jp (i0-1+jp).
 *           
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in a dpool, 
 *           the decks *must* be sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           alpha     - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_alpha - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *
 * Returns: log P(S|M)/P(S|R), as a bit score.
 */
float 
FInside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		float ***alpha, float ****ret_alpha, 
		struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int      status;
  float  **end;         /* we re-use the end deck. */
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  float    bsc;		/* total score for using local begin states */
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Lp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in FInside_b_jd_me(), dsq is NULL.\n");
  if(i0 != 1) 
    cm_Fail("ERROR: FInside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */

  if (dpool == NULL) dpool = deckpool_create();
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, 0, E_st);
  for (jp = 0; jp <= L; jp++) {
    j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
    end[j][0] = 0.;
    for (d = 1; d <= jp; d++) end[j][d] = IMPOSSIBLE;
  }

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v <= (cm->M-1); v++) touch[v] = cm->pnum[v];

  /* EPN: now that the EL self loop has a transition score, its
   *      necessary to keep track of the alpha EL deck (alpha[cm->M]).
   *      There's no bands on the EL state. 
   */
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      if (! deckpool_pop(dpool, &(alpha[cm->M]))) 
	alpha[cm->M] = alloc_vjd_deck(L, i0, j0);
      for (jp = 0; jp <= L; jp++) {
	j = i0-1+jp;
	/*alpha[cm->M][j][0] = IMPOSSIBLE;*/
	alpha[cm->M][j][0] = 0.;
	for (d = 1; d <= jp; d++)
	  {
	    alpha[cm->M][j][d] = (cm->el_selfsc * (d));
	  }
      }
    }

  /* Main recursion
   */
  for (v = (cm->M - 1); v >= 0; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      } 
      if (! deckpool_pop(dpool, &(alpha[v]))) 
	alpha[v] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);

      /* We've only allocated alpha cells that are within the bands
       * on the j and d dimensions. This means we have to deal
       * with all sorts of offset issues, but we don't have to 
       * waste time setting cells outside the bands to IMPOSSIBLE.
       */

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		  alpha[v][jp_v][dp_v]  = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  /* treat EL as emitting only on self transition */
		  for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		    {
		      yoffset = y - cm->cfirst[v];
		      if(j >= jmin[y] && j <= jmax[y]) 
			/* Enforces j is valid for state y */
			{
			  jp_y = j - jmin[y];
			  if(d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y])
			    {
			      dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							     in alpha w/mem eff bands */
			      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			       * corresponding to alpha[y][j][d] in the platonic matrix.
			       */
			      alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y][dp_y] 
										    + cm->tsc[v][yoffset]));
			    }
			}
		    }
		  if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
		}
	    }
	}
      else if (cm->sttype[v] == B_st)
	{
	  y = cm->cfirst[v];
	  z = cm->cnum[v];
	  /* Any valid j must be within both state v and state z's j band 
	   * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
	   * constructed, but we'll check anyway. 
	   */
	  tmp_jmin = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	  tmp_jmax = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];

	  /* For any values of j within v's j band but outside of z's j band,
	   * we have to set the corresponding alpha cells to IMPOSSIBLE.
	   * This is done be the following two ugly for loops, 
	   * which will only be looked at once for each B state, and
	   * even then only *very* rarely entered. This
	   * is why they're here, seemingly out of place before the 
	   * main j loop below, where similar performing code would be 
	   * looked at on the order of j times, instead of just once.
	   */
	  for(j = jmin[v]; j < tmp_jmin; j++)
	    {
	      jp_v = j-jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d-hdmin[v][jp_v];
		  alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
		}
	    }
	  if(tmp_jmax < jmax[v])
	    for(j = (tmp_jmax+1); j <= jmax[v]; j++)
	      {
		jp_v = j-jmin[v];
		for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		  {
		    dp_v = d-hdmin[v][jp_v];
		    alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
		  }
	      }
	  /* the main j loop */
	  for (j = tmp_jmin; j <= tmp_jmax; j++)
	    {
	      jp_v = j - jmin[v];
	      jp_y = j - jmin[y];
	      jp_z = j - jmin[z];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */

		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 6 inequalities (some may be redundant):
		   * (1) k >= j-jmax[y];
		   * (2) k <= j-jmin[y]; 
		   *     1 and 2 guarantee (j-k) is within state y's j band
		   *
		   * (3) k >= hdmin[z][j-jmin[z]];
		   * (4) k <= hdmax[z][j-jmin[z]]; 
		   *     3 and 4 guarantee k is within z's j=(j), d band
		   *
		   * (5) k >= d-hdmax[y][j-jmin[y]-k];
		   * (6) k <= d-hdmin[y][j-jmin[y]-k]; 
		   *     5 and 6 guarantee (d-k) is within state y's j=(j-k) d band
		   */
		  kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
		  /* kmin satisfies inequalities (1) and (3) */
		  kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
		  /* kmax satisfies inequalities (2) and (4) */
		  /* RHS of inequalities 5 and 6 are dependent on k, so we check
		   * for these within the next for loop.
		   */
		  alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* initialize */
		  for(k = kmin; k <= kmax; k++)
		    {
		      if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
			{
			  /* for current k, all 6 inequalities have been satisified 
			   * so we know the cells corresponding to the platonic 
			   * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
			   * alpha[z][j][k] are all within the bands. These
			   * cells correspond to alpha[v][jp_v][dp_v], 
			   * alpha[y][jp_y-k][d-hdmin[y][jp_y-k]-k],
			   * and alpha[z][jp_z][k-hdmin[z][jp_z]];
			   */
			  kp_z = k-hdmin[z][jp_z];
			  dp_y = d-hdmin[y][jp_y-k];

			  alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y-k][dp_y-k] 
										+ alpha[z][jp_z][kp_z]));
			}
		    }
		  if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
		  /* CYK Full ME Bands used 5 end block */
		}
	    }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces (j-1) is valid for state y */
		      {
			jp_y = j - jmin[y];
			if((d-2) >= hdmin[y][jp_y-1] && (d-2) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-2] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-2] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y-1][dp_y-2] 
										  + cm->tsc[v][yoffset]));
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][jp_v][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if(j >= jmin[y] && j <= jmax[y]) /* Enforces j is valid for state y */
		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y] && (d-1) <= hdmax[y][jp_y])
			  {
			    dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							   in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j][d-1] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y][dp_y-1] 
										  + cm->tsc[v][yoffset]));
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][dsq[i]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces j-1 is valid for state y */
		      
		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y-1] && (d-1) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-1] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = FLogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y-1][dp_y-1] 
										  + cm->tsc[v][yoffset]));
			  }
		      }
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][dsq[j]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}				/* finished calculating deck v. */
      
      /* Keep track of contributions of local begins */
      /* The following loops originally access alpha[v][j0][L] but the index L will be
	 in different positions due to the bands */
      if(j0 >= jmin[v] && j0 <= jmax[v])
	{
	  jp_v = j0 - jmin[v];
	  if(L >= hdmin[v][jp_v] && L <= hdmax[v][jp_v])
	    {
	      Lp = L - hdmin[v][jp_v];
	      /* If we get here alpha[v][jp_v][Lp] is a valid cell
	       * in the banded alpha matrix, corresponding to 
	       * alpha[v][j0][L] in the platonic matrix.
	       */
	      /* Check for local begin getting us to the root.
	       * This is "off-shadow": if/when we trace back, we'll handle this
	       * case separately (and we'll know to do it because we'll immediately
	       * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	       * to jump right to state b; see below)
	       */
	      if (allow_begin)
		{
		  bsc = FLogsum(bsc, (alpha[v][jp_v][Lp] + cm->beginsc[v]));
		}
	    }
	}

      /* If we're at the root state, record contribution of local begins */
      if (allow_begin && v == 0)
	{
	  if(j0 >= jmin[0] && j0 <= jmax[0])
	    {
	      jp_v = j0 - jmin[v];
	      if(L >= hdmin[v][jp_v] && L <= hdmax[v][jp_v])
		{
		  if (allow_begin && v == 0)
		    alpha[0][jp_v][Lp] = FLogsum(alpha[0][jp_v][Lp], bsc);
		}
	    }
	}	  

      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { 
	    /* Original code block : */
	    /* we can definitely release the S children of a bifurc. 
	       y = cm->cfirst[v]; deckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	       z = cm->cnum[v];   deckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	     End of original code block */
	    /* New ME code : */
	    y = cm->cfirst[v]; deckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	    z = cm->cnum[v];   deckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	    free_vjd_deck(alpha[y], i0, j0);
	    alpha[y] = NULL;
	    free_vjd_deck(alpha[z], i0, j0);
	    alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) 
		  {
		    if (cm->sttype[y] == E_st) { 
		      nends--; 
		      /* Original code : if (nends == 0) { deckpool_push(dpool, end); end = NULL;} */
		      /* ME code deletes the previous line, we don't mess with end, because
			 it is used later */
		    } else 
		      /* original code (deck reuse) deckpool_push(dpool, alpha[y]);*/
		      /* new ME code : */
		      {
			//printf("calling free vjd deck for alpha[y=%d]\n", y);
			free_vjd_deck(alpha[y], i0, j0);
		      }
		      alpha[y] = NULL;
		  }
	      }
	  }
      }
  } /* end loop over all v */

  /*
    printf("INSIDE JD, printing alpha\n");
    debug_print_alpha_banded_jd(alpha, cm, L, jmin, jmax, hdmin, hdmax);
    printf("INSIDE JD, done printing alpha\n");
  */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  Lp = L - hdmin[0][j0-jmin[0]];
  sc = alpha[0][j0-jmin[0]][Lp];

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { deckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { deckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (deckpool_pop(dpool, &end)) free_vjd_deck(end, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  ESL_DPRINTF1(("\n\tFInside_b_jd_me() sc  : %f\n", sc));
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

float 
IInside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		int   ***alpha, int   ****ret_alpha, 
		Ideckpool_t *dpool, Ideckpool_t **ret_dpool,
		int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int      status;       
  int    **end;         /* we re-use the end deck. */
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    return_sc;   /* the return score, converted to bits (Scorified) */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  int      bsc;		/* total score for using local begin states */
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Lp;             /* L also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in IInside_b_jd_me(), dsq is NULL.\n");
  if(i0 != 1) 
    cm_Fail("ERROR: IInside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");

  /* Allocations and initializations
   */
  bsc = -INFTY;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */

  if (dpool == NULL) dpool = Ideckpool_create();
  if (! Ideckpool_pop(dpool, &end))
    end = Ialloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, 0, E_st);
  for (jp = 0; jp <= L; jp++) {
    j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
    end[j][0] = 0.;
    for (d = 1; d <= jp; d++) end[j][d] = -INFTY;
  }

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(int   **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v <= (cm->M-1); v++) touch[v] = cm->pnum[v];

  /* EPN: now that the EL self loop has a transition score, its
   *      necessary to keep track of the alpha EL deck (alpha[cm->M]).
   *      There's no bands on the EL state. 
   */
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      if (! Ideckpool_pop(dpool, &(alpha[cm->M]))) 
	alpha[cm->M] = Ialloc_vjd_deck(L, i0, j0);
      for (jp = 0; jp <= L; jp++) {
	j = i0-1+jp;
	/*alpha[cm->M][j][0] = -INFTY;*/
	alpha[cm->M][j][0] = 0.;
	for (d = 1; d <= jp; d++)
	  {
	    alpha[cm->M][j][d] = (cm->iel_selfsc * (d));
	  }
      }
    }

  /* Main recursion
   */
  for (v = (cm->M - 1); v >= 0; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      } 
      if (! Ideckpool_pop(dpool, &(alpha[v]))) 
	alpha[v] = Ialloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);

      /* We've only allocated alpha cells that are within the bands
       * on the j and d dimensions. This means we have to deal
       * with all sorts of offset issues, but we don't have to 
       * waste time setting cells outside the bands to -INFTY.
       */

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		  alpha[v][jp_v][dp_v]  = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  /* treat EL as emitting only on self transition */
		  for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		    {
		      yoffset = y - cm->cfirst[v];
		      if(j >= jmin[y] && j <= jmax[y]) 
			/* Enforces j is valid for state y */
			{
			  jp_y = j - jmin[y];
			  if(d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y])
			    {
			      dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							     in alpha w/mem eff bands */
			      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			       * corresponding to alpha[y][j][d] in the platonic matrix.
			       */
			      alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y][dp_y] 
										    + cm->itsc[v][yoffset]));
			    }
			}
		    }
		  if (alpha[v][jp_v][dp_v] < -INFTY) alpha[v][jp_v][dp_v] = -INFTY;
		}
	    }
	}
      else if (cm->sttype[v] == B_st)
	{
	  y = cm->cfirst[v];
	  z = cm->cnum[v];
	  /* Any valid j must be within both state v and state z's j band 
	   * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
	   * constructed, but we'll check anyway. 
	   */
	  tmp_jmin = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	  tmp_jmax = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];

	  /* For any values of j within v's j band but outside of z's j band,
	   * we have to set the corresponding alpha cells to -INFTY.
	   * This is done be the following two ugly for loops, 
	   * which will only be looked at once for each B state, and
	   * even then only *very* rarely entered. This
	   * is why they're here, seemingly out of place before the 
	   * main j loop below, where similar performing code would be 
	   * looked at on the order of j times, instead of just once.
	   */
	  for(j = jmin[v]; j < tmp_jmin; j++)
	    {
	      jp_v = j-jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d-hdmin[v][jp_v];
		  alpha[v][jp_v][dp_v] = -INFTY; /* this won't be changed */
		}
	    }
	  if(tmp_jmax < jmax[v])
	    for(j = (tmp_jmax+1); j <= jmax[v]; j++)
	      {
		jp_v = j-jmin[v];
		for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		  {
		    dp_v = d-hdmin[v][jp_v];
		    alpha[v][jp_v][dp_v] = -INFTY; /* this won't be changed */
		  }
	      }
	  /* the main j loop */
	  for (j = tmp_jmin; j <= tmp_jmax; j++)
	    {
	      jp_v = j - jmin[v];
	      jp_y = j - jmin[y];
	      jp_z = j - jmin[z];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */

		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 6 inequalities (some may be redundant):
		   * (1) k >= j-jmax[y];
		   * (2) k <= j-jmin[y]; 
		   *     1 and 2 guarantee (j-k) is within state y's j band
		   *
		   * (3) k >= hdmin[z][j-jmin[z]];
		   * (4) k <= hdmax[z][j-jmin[z]]; 
		   *     3 and 4 guarantee k is within z's j=(j), d band
		   *
		   * (5) k >= d-hdmax[y][j-jmin[y]-k];
		   * (6) k <= d-hdmin[y][j-jmin[y]-k]; 
		   *     5 and 6 guarantee (d-k) is within state y's j=(j-k) d band
		   */
		  kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
		  /* kmin satisfies inequalities (1) and (3) */
		  kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
		  /* kmax satisfies inequalities (2) and (4) */
		  /* RHS of inequalities 5 and 6 are dependent on k, so we check
		   * for these within the next for loop.
		   */
		  alpha[v][jp_v][dp_v] = -INFTY; /* initialize */
		  for(k = kmin; k <= kmax; k++)
		    {
		      if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
			{
			  /* for current k, all 6 inequalities have been satisified 
			   * so we know the cells corresponding to the platonic 
			   * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
			   * alpha[z][j][k] are all within the bands. These
			   * cells correspond to alpha[v][jp_v][dp_v], 
			   * alpha[y][jp_y-k][d-hdmin[y][jp_y-k]-k],
			   * and alpha[z][jp_z][k-hdmin[z][jp_z]];
			   */
			  kp_z = k-hdmin[z][jp_z];
			  dp_y = d-hdmin[y][jp_y-k];

			  alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y-k][dp_y-k] 
										+ alpha[z][jp_z][kp_z]));
			}
		    }
		  if (alpha[v][jp_v][dp_v] < -INFTY) alpha[v][jp_v][dp_v] = -INFTY;
		  /* CYK Full ME Bands used 5 end block */
		}
	    }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces (j-1) is valid for state y */
		      {
			jp_y = j - jmin[y];
			if((d-2) >= hdmin[y][jp_y-1] && (d-2) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-2] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-2] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y-1][dp_y-2] 
										  + cm->itsc[v][yoffset]));
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->iesc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][jp_v][dp_v] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);
		
		if (alpha[v][jp_v][dp_v] < -INFTY) alpha[v][jp_v][dp_v] = -INFTY;
	      }
	    }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if(j >= jmin[y] && j <= jmax[y]) /* Enforces j is valid for state y */
		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y] && (d-1) <= hdmax[y][jp_y])
			  {
			    dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							   in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j][d-1] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y][dp_y-1] 
										  + cm->itsc[v][yoffset]));
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->iesc[v][dsq[i]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);
		if (alpha[v][jp_v][dp_v] < -INFTY) alpha[v][jp_v][dp_v] = -INFTY;
	      }
	    }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces j-1 is valid for state y */
		      
		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y-1] && (d-1) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-1] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = ILogsum(alpha[v][jp_v][dp_v], (alpha[y][jp_y-1][dp_y-1] 
										  + cm->itsc[v][yoffset]));
			  }
		      }
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->iesc[v][dsq[j]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);
		
		if (alpha[v][jp_v][dp_v] < -INFTY) alpha[v][jp_v][dp_v] = -INFTY;
	      }
	    }
	}				/* finished calculating deck v. */
      
      /* Keep track of contributions of local begins */
      /* The following loops originally access alpha[v][j0][L] but the index L will be
	 in different positions due to the bands */
      if(j0 >= jmin[v] && j0 <= jmax[v])
	{
	  jp_v = j0 - jmin[v];
	  if(L >= hdmin[v][jp_v] && L <= hdmax[v][jp_v])
	    {
	      Lp = L - hdmin[v][jp_v];
	      /* If we get here alpha[v][jp_v][Lp] is a valid cell
	       * in the banded alpha matrix, corresponding to 
	       * alpha[v][j0][L] in the platonic matrix.
	       */
	      /* Check for local begin getting us to the root.
	       * This is "off-shadow": if/when we trace back, we'll handle this
	       * case separately (and we'll know to do it because we'll immediately
	       * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	       * to jump right to state b; see below)
	       */
	      if (allow_begin)
		{
		  bsc = ILogsum(bsc, (alpha[v][jp_v][Lp] + cm->ibeginsc[v]));
		}
	    }
	}

      /* If we're at the root state, record contribution of local begins */
      if (allow_begin && v == 0)
	{
	  if(j0 >= jmin[0] && j0 <= jmax[0])
	    {
	      jp_v = j0 - jmin[v];
	      if(L >= hdmin[v][jp_v] && L <= hdmax[v][jp_v])
		{
		  if (allow_begin && v == 0)
		    alpha[0][jp_v][Lp] = ILogsum(alpha[0][jp_v][Lp], bsc);
		}
	    }
	}	  

      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { 
	    /* Original code block : */
	    /* we can definitely release the S children of a bifurc. 
	       y = cm->cfirst[v]; Ideckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	       z = cm->cnum[v];   Ideckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	     End of original code block */
	    /* New ME code : */
	    y = cm->cfirst[v]; Ideckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	    z = cm->cnum[v];   Ideckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	    Ifree_vjd_deck(alpha[y], i0, j0);
	    alpha[y] = NULL;
	    Ifree_vjd_deck(alpha[z], i0, j0);
	    alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) 
		  {
		    if (cm->sttype[y] == E_st) { 
		      nends--; 
		      /* Original code : if (nends == 0) { Ideckpool_push(dpool, end); end = NULL;} */
		      /* ME code deletes the previous line, we don't mess with end, because
			 it is used later */
		    } else 
		      /* original code (deck reuse) Ideckpool_push(dpool, alpha[y]);*/
		      /* new ME code : */
		      {
			/* printf("calling free vjd deck for alpha[y=%d]\n", y); */
			Ifree_vjd_deck(alpha[y], i0, j0);
		      }
		      alpha[y] = NULL;
		  }
	      }
	  }
      }
  } /* end loop over all v */

  /*
    printf("INSIDE JD, printing alpha\n");
    debug_print_alpha_banded_jd(alpha, cm, L, jmin, jmax, hdmin, hdmax);
    printf("INSIDE JD, done printing alpha\n");
  */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  Lp = L - hdmin[0][j0-jmin[0]];
  return_sc = Scorify(alpha[0][j0-jmin[0]][Lp]);

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ideckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { Ideckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (Ideckpool_pop(dpool, &end)) Ifree_vjd_deck(end, i0, j0);
    Ideckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  ESL_DPRINTF1(("\tIInside_b_jd_me() sc  : %f\n", return_sc));
  return return_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/***********************************************************************
 * Function: FOutside_b_jd_me()
 * Date:     EPN 05.26.06
 * SRE, Tue Aug  8 10:42:52 2000 [St. Louis]
 *
 * Purpose:  The Outside dynamic programming algorithm for CMs.
 *           Works directly with floats, stepping into and out of 
 *           log space as necessary.
 *  
 *           Derived from smallcyk.c::CYKOutside() and smallcyk.c::outsdie(). 
 * 
 *           Align a subsequence to the full model, i.e. we're given
 *           i0 and j0, beginning and end positions of the subseq we're
 *           considering.
 * 
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha and beta, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           beta       - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_beta  - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           alpha     - the alpha matrix from a Inside run, if do_check is FALSE
 *                        only decks for S states must be valid, else all must be
 *                        valid.
 *           ret_alpha - if non-NULL, return the alpha matrix with one or more
 *                       decks available for examination (see "do_full")
 *           do_check  - TRUE to do time-consuming check to make sure
 *                       beta and alpha are consistent (only if NON-LOCAL mode)
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * 
 * Returns: log P(S|M)/P(S|R), as a bit score.
 ***********************************************************************/
float 
FOutside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		 float ***beta, float ****ret_beta, 
		 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		 int allow_begin, float ***alpha, float ****ret_alpha, 
		 int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int      status;
  int      v,y,z;	/* indices for states */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  float    escore;	/* an emission score, tmp variable */
  int      voffset;	/* index of v in t_v(y) transition scores */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  float    bsc;		/* total score for using local begin states */
  float    return_sc;   /* P(S|M)/P(S|R) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  float  **end;         /* we re-use the end deck. */
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Lp;             /* L also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;
  int      fail_flag = FALSE; /* set to TRUE if do_check and we see a problem */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in FOutside_b_jd_me(), dsq is NULL.\n");
  if(i0 != 1) 
    cm_Fail("ERROR: FOutside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");
  if (cm->flags & CMH_LOCAL_END) { do_check = FALSE; } 
  /* Code for checking doesn't apply in local mode. See below. */

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* if caller didn't give us a matrix, make one.
   * Allocate room for M+1 decks because we might need the EL deck (M)
   * if we're doing local alignment.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck. Root is necessarily the ROOT_S state 0.
   */
  if (! deckpool_pop(dpool, &(beta[0])))
    beta[0] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[0], jmax[0], hdmin[0], hdmax[0]);
  for (j = jmin[0]; j <= jmax[0]; j++)
    {
      jp_v = j - jmin[0];
      for (d = hdmin[0][jp_v]; d <= hdmax[0][jp_v]; d++)
	{
	  dp_v = d - hdmin[0][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  beta[0][jp_v][dp_v] = IMPOSSIBLE;
	}
    }
  /* non banded line: beta[0][j0][L] = 0; */
  jp_v = j0 - jmin[0];
  Lp = L - hdmin[0][jp_v];
  assert(L >= hdmin[0][jp_v]);
  assert(L <= hdmax[0][jp_v]);
  beta[0][jp_v][Lp] = 0.;

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   * EL deck has no bands as currently implemented.
   */
  if (cm->flags & CMH_LOCAL_END) {
    if (! deckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = alloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= L; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[cm->M][j][d] = IMPOSSIBLE;
    }
    
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  /* Main loop down through the decks
   */
  /*for (v = 2; v < cm->M; v++) */ /*EPN is this 2 b/c Durbin p.287 
				     has state 2 in the algorithm? b/c state 1 is root*/
  for (v = 1; v < cm->M; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);

      /* Init the whole deck to IMPOSSIBLE
       */
      for (j = jmin[v]; j <= jmax[v]; j++)
	{
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	    {
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      beta[v][jp_v][dp_v] = IMPOSSIBLE;
	    }
	}

      /* If we can do a local begin into v, also init with that. 
       * By definition, beta[0][j0][L] == 0.
       */ 
      if (cm->flags & CMH_LOCAL_BEGIN)
	{
	  if((j0 >= jmin[v]) && (j0 <= jmax[v]))
	    {
	      jp_v = j0 - jmin[v];
	      if((L >= hdmin[v][jp_v]) && L <= hdmax[v][jp_v])
		{
		  Lp = L - hdmin[v][jp_v];
		  beta[v][jp_v][Lp] = cm->beginsc[v];
		}
	    }
	}
      /* main recursion: reorganized relative to FOutside() for simplification of
       * band-related issues.
       */
      for (j = jmax[v]; j >= jmin[v]; j--)
	{
	  jp_v = j - jmin[v];
	  for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--)
	    {
	      i = j-d+1;
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	      if (cm->stid[v] == BEGL_S) 
		{
		  y = cm->plast[v];	/* the parent bifurcation    */
		  z = cm->cnum[y];	/* the other (right) S state */
		  jp_y = j - jmin[y];
		  jp_z = j - jmin[z];
		  
		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 8 inequalities (some may be redundant):
		   * NOTE: these are different from those in Inside() (for one thing, v and y
		   *       (BEGL_S and BIF_B here respectively) are switched relative to Inside.
		   *
		   * (1) k <= jmax[y] - j;
		   * (2) k >= jmin[y] - j;
		   * (3) k <= jmax[z] - j;
		   * (4) k >= jmin[z] - j;
		   *     1 and 2 guarantee (j+k) is within state y's j band
		   *     3 and 4 guarantee (j+k) is within state z's j band
		   *
		   * (5) k >= hdmin[y][j-jmin[y]+k] - d;
		   * (6) k <= hdmax[y][j-jmin[y]+k] - d; 
		   *     5 and 6 guarantee k+d is within y's j=(j+k), d band
		   *
		   * (7) k >= hdmin[z][j-jmin[z]+k];
		   * (8) k <= hdmax[z][j-jmin[z]+k]; 
		   *     5 and 6 guarantee k is within state z's j=(j+k) d band
		   * 
		   * Note, below:
		   * tmp_jmin = MAX(jmin[y], jmin[z];
		   * tmp_jmax = MIN(jmax[y], jmax[z];
		   */
		  tmp_jmin = (jmin[y] > jmin[z]) ? jmin[y] : jmin[z];
		  tmp_jmax = (jmax[y] < jmax[z]) ? jmax[y] : jmax[z];

		  kmin = tmp_jmin - j;
		  kmax = tmp_jmax - j;
		  /* kmin and kmax satisfy inequalities (1-4) */
		  /* RHS of inequalities 5-8 are dependent on k, so we check
		   * for these within the next for loop.
		   */
		  for(k = kmin; k <= kmax; k++)
		    {
		      if(k < (hdmin[y][jp_y+k] - d) || k > (hdmax[y][jp_y+k] - d)) continue; 
		      /* above line continues if inequality 5 or 6 is violated */
		      if(k < (hdmin[z][jp_z+k]) || k > (hdmax[z][jp_z+k])) continue; 
		      /* above line continues if inequality 7 or 8 is violated */
		      
		      /* if we get here for current k, all 8 inequalities have been satisified 
		       * so we know the cells corresponding to the platonic 
		       * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
		       * alpha[z][j+k][k] are all within the bands. These
		       * cells correspond to beta[v][jp_v][dp_v], 
		       * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
		       * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
		       */
		      kp_z = k-hdmin[z][jp_z+k];
		      dp_y = d-hdmin[y][jp_y+k];
		      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] 
									  + alpha[z][jp_z+k][kp_z]));
		    }
		}
	      else if (cm->stid[v] == BEGR_S) 
		{
		  y = cm->plast[v];	/* the parent bifurcation    */
		  z = cm->cfirst[y];	/* the other (left) S state */

		  jp_y = j - jmin[y];
		  jp_z = j - jmin[z];
		  
		  /* For j to be valid for state y: *
		   * jmin[y] >= j >= jmax[y]
		   * These are independent of k so we check outside of k loop below
		   * For j to be valid for state z: *
		   * (jmin[z] + d) >= j >= (jmax[z] + d)
		   */
		  if(j < jmin[y] || j > jmax[y]) continue;
		  if((j < (jmin[z] + d)) || (j > (jmax[z]+d))) continue;

		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 4 inequalities (some may be redundant):
		   * NOTE: these are different from those in Inside() (for one thing, v and y
		   *       (BEGR_S and BIF_B here respectively) are switched relative to Inside.
		   *
		   * (1) k >= hdmin[y][j-jmin[y]] - d;
		   * (2) k <= hdmax[y][j-jmin[y]] - d;
		   *     1 and 2 guarantee (d+k) is within state y's j=(j) d band
		   *
		   * (3) k >= hdmin[z][j-jmin[z]-d];
		   * (4) k <= hdmax[z][j-jmin[z]-d];
		   *     3 and 4 guarantee k is within z's j=(j-d), d band
		   *
		   */
		  kmin = ((hdmin[y][jp_y]-d) > (hdmin[z][jp_z-d])) ? (hdmin[y][jp_y]-d) : (hdmin[z][jp_z-d]);
		  /* kmin satisfies inequalities (1) and (3) */
		  kmax = ((hdmax[y][jp_y]-d) < (hdmax[z][jp_z-d])) ? (hdmax[y][jp_y]-d) : (hdmax[z][jp_z-d]);
		  /* kmax satisfies inequalities (2) and (4) */

		  for(k = kmin; k <= kmax; k++)
		    {
		      /* for current k, all 4 inequalities have been satisified 
		       * so we know the cells corresponding to the platonic 
		       * matrix cells beta[v][j][d], beta[y][j][d+k], and
		       * alpha[z][j-d][k] are all within the bands. These
		       * cells correspond to beta[v][jp_v][dp_v], 
		       * beta[y][jp_y+k][d-hdmin[y][jp_y]+k],
		       * and alpha[z][jp_z-d][k-hdmin[z][jp_z-d]];
		       */
		      kp_z = k-hdmin[z][jp_z-d];
		      dp_y = d-hdmin[y][jp_y];
		      beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
									  + alpha[z][jp_z-d][kp_z]));
		    }
		}
	      else
		{
		  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) 
		    {
		      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

		      switch(cm->sttype[y]) {
		      case MP_st: 
			if (j == j0 || d == j) continue; /* boundary condition */
			if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+2) < hdmin[y][(jp_y+1)] || (d+2) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
			 * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y+1];  /* d index for state y */
			
			if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
			  escore = cm->esc[y][(dsq[i-1]*cm->abc->K+dsq[j+1])];
			else
			  escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
			beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
									    + cm->tsc[y][voffset] + escore));
			break;

		      case ML_st:
		      case IL_st: 
			if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
			if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+1) < hdmin[y][jp_y] || (d+1) > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
			 * corresponding to alpha[y][j][d+1] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y];  /* d index for state y */
			if (dsq[i-1] < cm->abc->K) 
			  escore = cm->esc[y][dsq[i-1]];
			else
			  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
			beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
									    + cm->tsc[y][voffset] + escore));
			break;
		    
		      case MR_st:
		      case IR_st:
			if (j == j0) continue;
			if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+1) < hdmin[y][(jp_y+1)] || (d+1) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
			 * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][(jp_y+1)];  /* d index for state y */
			if (dsq[j+1] < cm->abc->K) 
			  escore = cm->esc[y][dsq[j+1]];
			else
			  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);
			/*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
			beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
									    + cm->tsc[y][voffset] + escore));
			break;

		      case S_st:
		      case E_st:
		      case D_st:
			if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if (d < hdmin[y][jp_y] || d > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			 * corresponding to alpha[y][j][d] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y];  /* d index for state y */
			beta[v][jp_v][dp_v] = FLogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->tsc[y][voffset])); 
			break;

		      default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
		      }/* end switch over states*/
		  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
		if (beta[v][jp_v][dp_v] < IMPOSSIBLE) beta[v][jp_v][dp_v] = IMPOSSIBLE;
		}	/* ends else entered for non-BEGL/BEGR states*/	
	    } /* ends loop over d. We know all beta[v][j][d] in this row j*/
	}    /* end loop over jp. We know the beta's for the whole deck.*/
      /* Deal with local alignment end transitions v->EL
       * (EL = deck at M.)
       */
      if (NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (j = jmin[v]; j <= jmax[v]; j++)
	  {
	    jp_v = j - jmin[v];
	    for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		i = j-d+1;
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		switch (cm->sttype[v]) {
		case MP_st: 
		  if (j == j0 || d == j) continue; /* boundary condition */
		  if (((j+1) > jmax[v]) || ((d+2) > hdmax[v][jp_v])) continue; /*boundary condition*/
		  if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		    escore = cm->esc[v][(dsq[i-1]*cm->abc->K+dsq[j+1])];
		  else
		    escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		  beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+2] + cm->endsc[v] 
								  + escore));
		break;
	      case ML_st:
	      case IL_st:
		if (d == j) continue;	
		if ((d+1) > hdmax[v][jp_v]) continue; /*boundary condition*/
		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[v][dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v+1] + cm->endsc[v] 
								+ escore));
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (((j+1) > jmax[v]) || ((d+1) > hdmax[v][jp_v])) continue; /*boundary condition*/
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[v][dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+1] + cm->endsc[v] 
								+ escore));
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v] + cm->endsc[v] 
								+ escore));
		break;
	      case B_st:  
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
		/* note that although B is a valid vend for a segment we'd do
                   outside on, B->EL is set to be impossible, by the local alignment
                   config. There's no point in having a B->EL because B is a nonemitter
                   (indeed, it would introduce an alignment ambiguity). The same
		   alignment case is handled by the X->EL transition where X is the
		   parent consensus state (S, MP, ML, or MR) above the B. Thus,
		   this code is relying on the NOT_IMPOSSIBLE() test, above,
		   to make sure the sttype[vend]=B case gets into this switch.
		*/
	      } /* end switch over parent state type v */
	    } /* end inner loop over d */
	} /* end outer loop over jp */
      } /* end conditional section for dealing w/ v->EL local end transitions */

      /* Look at v's parents; if we're reusing memory (! do_full)
       * push the parents that we don't need any more into the pool.
       */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { deckpool_push(dpool, beta[y]); beta[y] = NULL; }
	}
      }
    } /* end loop over decks v. */

  /* EPN left the below SRE block out */
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = L; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = FLogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->el_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE JD, printing beta\n");
    debug_print_alpha_banded_jd(beta, cm, L, jmin, jmax, hdmin, hdmax);
    printf("OUTSIDE JD, done printing beta\n");
  */

  Lp = L - hdmin[0][j0-jmin[0]];
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) 
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    {
      /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
       * using both the Outside (beta) and Inside (alpha) matrices,
       * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
       * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
       *                                                = P(S|M) / P(S|R)  
       * in v,j,d coordinates this is:
       * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
       *                                                = P(S|M) / P(S|R)
       */
	 
      for(n = 0; n < cm->nodes; n++)
	{
	  sc = IMPOSSIBLE;
	  num_split_states = SplitStatesInNode(cm->ndtype[n]);
	  for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++)
	    {
	      for (j = jmin[v]; j <= jmax[v]; j++)
		{
		  jp_v = j - jmin[v];
		  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		    {
		      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		      /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
			printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
		      sc = FLogsum(sc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
		    }
		}
	      }
	  /*printf("checking node: %d | sc: %.6f\n", n, sc);*/
	  diff = alpha[0][j0-jmin[0]][Lp] - sc;
	  if(diff > 0.001 || diff < -0.001)
	    {
	      fail_flag = TRUE;
	      printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
		     n, sc, alpha[0][(j0-jmin[0])][Lp], diff);
	    }
	}

    }

  /* IF not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * IF local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END))
    {
      return_sc = IMPOSSIBLE;
      v = cm->M-1;
      for (j = jmin[v]; j <= jmax[v]; j++)
	{
	  jp_v = j - jmin[v];
	  assert(hdmin[v][jp_v] == 0);
	  /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	  return_sc = FLogsum(return_sc, (beta[v][jp_v][0]));
	}
    }
  else /* return_sc = P(S|M) / P(S|R) from Inside() */
    {
      return_sc = alpha[0][(j0-jmin[0])][Lp];
    }

  /* If the caller doesn't want the beta matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
      deckpool_push(dpool, beta[cm->M]);
      beta[cm->M] = NULL; 
    }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the alpha matrix, free it 
   * Else, pass it back to him.
   * EPN - if we free the alpha and beta matrix the deck pool has all the 
   *       decks from alpha and beta, not sure if this is desirable.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { deckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { deckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    float **a;
    while (deckpool_pop(dpool, &a)) free_vjd_deck(a, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }
  free(touch);

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");
  
  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tFOutside_b_jd_me() sc : %f\n", return_sc));
  else
    ESL_DPRINTF1(("\tFOutside_b_jd_me() sc : %f (LOCAL mode; sc is from Inside)\n", return_sc));

  return return_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

float 
IOutside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int do_full,
		 int   ***beta, int   ****ret_beta, 
		 Ideckpool_t *dpool, Ideckpool_t **ret_dpool,
		 int allow_begin, int   ***alpha, int   ****ret_alpha, 
		 int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int      status;
  int      v,y,z;	/* indices for states */
  int      j,d,i,k;	/* indices in sequence dimensions */
  int      isc;     	/* a temporary variable holding an int score */
  float    fsc;     	/* a temporary variable holding a float score */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      iescore;	/* an emission score, tmp variable */
  int      voffset;	/* index of v in t_v(y) transition scores */
  int      L;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  int      bsc;		/* total score for using local begin states */
  int      ireturn_sc;  /* P(S|M)/P(S|R), a scaled int*/
  float    freturn_sc;  /* P(S|M)/P(S|R), a float (Scorified ireturn_sc) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  int    **end;         /* we re-use the end deck. */
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Lp;             /* L also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;
  int      fail_flag = FALSE; /* set to TRUE if do_check and we see a problem */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in IOutside_b_jd_me(), dsq is NULL.\n");
  if(i0 != 1) 
    cm_Fail("ERROR: IOutside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");

  if (cm->flags & CMH_LOCAL_END) { do_check = FALSE; } 
  /* Code for checking doesn't apply in local mode. See below. */

  /* Allocations and initializations
   */
  bsc = -INFTY;
  L   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = Ideckpool_create();

  /* if caller didn't give us a matrix, make one.
   * Allocate room for M+1 decks because we might need the EL deck (M)
   * if we're doing local alignment.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(int   **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck. Root is necessarily the ROOT_S state 0.
   */
  if (! Ideckpool_pop(dpool, &(beta[0])))
    beta[0] = Ialloc_jdbanded_vjd_deck(L, i0, j0, jmin[0], jmax[0], hdmin[0], hdmax[0]);
  for (j = jmin[0]; j <= jmax[0]; j++)
    {
      jp_v = j - jmin[0];
      for (d = hdmin[0][jp_v]; d <= hdmax[0][jp_v]; d++)
	{
	  dp_v = d - hdmin[0][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  beta[0][jp_v][dp_v] = -INFTY;
	}
    }
  /* non banded line: beta[0][j0][L] = 0; */
  jp_v = j0 - jmin[0];
  Lp = L - hdmin[0][jp_v];
  assert(L >= hdmin[0][jp_v]);
  assert(L <= hdmax[0][jp_v]);
  beta[0][jp_v][Lp] = 0.;

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   * EL deck has no bands as currently implemented.
   */
  if (cm->flags & CMH_LOCAL_END) {
    if (! Ideckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = Ialloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= L; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[cm->M][j][d] = -INFTY;
    }
    
    /* Le don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  /* Main loop down through the decks
   */
  /*for (v = 2; v < cm->M; v++) */ /*EPN is this 2 b/c Durbin p.287 
				     has state 2 in the algorithm? b/c state 1 is root? */
  for (v = 1; v < cm->M; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! Ideckpool_pop(dpool, &(beta[v])))
	beta[v] = Ialloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);

      /* Init the whole deck to -INFTY
       */
      for (j = jmin[v]; j <= jmax[v]; j++)
	{
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	    {
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      beta[v][jp_v][dp_v] = -INFTY;
	    }
	}

      /* If we can do a local begin into v, also init with that. 
       * By definition, beta[0][j0][L] == 0.
       */ 
      if (cm->flags & CMH_LOCAL_BEGIN)
	{
	  if((j0 >= jmin[v]) && (j0 <= jmax[v]))
	    {
	      jp_v = j0 - jmin[v];
	      if((L >= hdmin[v][jp_v]) && L <= hdmax[v][jp_v])
		{
		  Lp = L - hdmin[v][jp_v];
		  beta[v][jp_v][Lp] = cm->ibeginsc[v];
		}
	    }
	}
      /* main recursion: reorganized relative to FOutside() for simplification of
       * band-related issues.
       */
      for (j = jmax[v]; j >= jmin[v]; j--)
	{
	  jp_v = j - jmin[v];
	  for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--)
	    {
	      i = j-d+1;
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	      if (cm->stid[v] == BEGL_S) 
		{
		  y = cm->plast[v];	/* the parent bifurcation    */
		  z = cm->cnum[y];	/* the other (right) S state */
		  jp_y = j - jmin[y];
		  jp_z = j - jmin[z];
		  
		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 8 inequalities (some may be redundant):
		   * NOTE: these are different from those in Inside() (for one thing, v and y
		   *       (BEGL_S and BIF_B here respectively) are switched relative to Inside.
		   *
		   * (1) k <= jmax[y] - j;
		   * (2) k >= jmin[y] - j;
		   * (3) k <= jmax[z] - j;
		   * (4) k >= jmin[z] - j;
		   *     1 and 2 guarantee (j+k) is within state y's j band
		   *     3 and 4 guarantee (j+k) is within state z's j band
		   *
		   * (5) k >= hdmin[y][j-jmin[y]+k] - d;
		   * (6) k <= hdmax[y][j-jmin[y]+k] - d; 
		   *     5 and 6 guarantee k+d is within y's j=(j+k), d band
		   *
		   * (7) k >= hdmin[z][j-jmin[z]+k];
		   * (8) k <= hdmax[z][j-jmin[z]+k]; 
		   *     5 and 6 guarantee k is within state z's j=(j+k) d band
		   * 
		   * Note, below:
		   * tmp_jmin = MAX(jmin[y], jmin[z];
		   * tmp_jmax = MIN(jmax[y], jmax[z];
		   */
		  tmp_jmin = (jmin[y] > jmin[z]) ? jmin[y] : jmin[z];
		  tmp_jmax = (jmax[y] < jmax[z]) ? jmax[y] : jmax[z];

		  kmin = tmp_jmin - j;
		  kmax = tmp_jmax - j;
		  /* kmin and kmax satisfy inequalities (1-4) */
		  /* RHS of inequalities 5-8 are dependent on k, so we check
		   * for these within the next for loop.
		   */
		  for(k = kmin; k <= kmax; k++)
		    {
		      if(k < (hdmin[y][jp_y+k] - d) || k > (hdmax[y][jp_y+k] - d)) continue; 
		      /* above line continues if inequality 5 or 6 is violated */
		      if(k < (hdmin[z][jp_z+k]) || k > (hdmax[z][jp_z+k])) continue; 
		      /* above line continues if inequality 7 or 8 is violated */
		      
		      /* if we get here for current k, all 8 inequalities have been satisified 
		       * so we know the cells corresponding to the platonic 
		       * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
		       * alpha[z][j+k][k] are all within the bands. These
		       * cells correspond to beta[v][jp_v][dp_v], 
		       * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
		       * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
		       */
		      kp_z = k-hdmin[z][jp_z+k];
		      dp_y = d-hdmin[y][jp_y+k];
		      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] 
									  + alpha[z][jp_z+k][kp_z]));
		    }
		}
	      else if (cm->stid[v] == BEGR_S) 
		{
		  y = cm->plast[v];	/* the parent bifurcation    */
		  z = cm->cfirst[y];	/* the other (left) S state */

		  jp_y = j - jmin[y];
		  jp_z = j - jmin[z];
		  
		  /* For j to be valid for state y: *
		   * jmin[y] >= j >= jmax[y]
		   * These are independent of k so we check outside of k loop below
		   * For j to be valid for state z: *
		   * (jmin[z] + d) >= j >= (jmax[z] + d)
		   */
		  if(j < jmin[y] || j > jmax[y]) continue;
		  if((j < (jmin[z] + d)) || (j > (jmax[z]+d))) continue;

		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 4 inequalities (some may be redundant):
		   * NOTE: these are different from those in Inside() (for one thing, v and y
		   *       (BEGR_S and BIF_B here respectively) are switched relative to Inside.
		   *
		   * (1) k >= hdmin[y][j-jmin[y]] - d;
		   * (2) k <= hdmax[y][j-jmin[y]] - d;
		   *     1 and 2 guarantee (d+k) is within state y's j=(j) d band
		   *
		   * (3) k >= hdmin[z][j-jmin[z]-d];
		   * (4) k <= hdmax[z][j-jmin[z]-d];
		   *     3 and 4 guarantee k is within z's j=(j-d), d band
		   *
		   */
		  kmin = ((hdmin[y][jp_y]-d) > (hdmin[z][jp_z-d])) ? (hdmin[y][jp_y]-d) : (hdmin[z][jp_z-d]);
		  /* kmin satisfies inequalities (1) and (3) */
		  kmax = ((hdmax[y][jp_y]-d) < (hdmax[z][jp_z-d])) ? (hdmax[y][jp_y]-d) : (hdmax[z][jp_z-d]);
		  /* kmax satisfies inequalities (2) and (4) */

		  for(k = kmin; k <= kmax; k++)
		    {
		      /* for current k, all 4 inequalities have been satisified 
		       * so we know the cells corresponding to the platonic 
		       * matrix cells beta[v][j][d], beta[y][j][d+k], and
		       * alpha[z][j-d][k] are all within the bands. These
		       * cells correspond to beta[v][jp_v][dp_v], 
		       * beta[y][jp_y+k][d-hdmin[y][jp_y]+k],
		       * and alpha[z][jp_z-d][k-hdmin[z][jp_z-d]];
		       */
		      kp_z = k-hdmin[z][jp_z-d];
		      dp_y = d-hdmin[y][jp_y];
		      beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
									  + alpha[z][jp_z-d][kp_z]));
		    }
		}
	      else
		{
		  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) 
		    {
		      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

		      switch(cm->sttype[y]) {
		      case MP_st: 
			if (j == j0 || d == j) continue; /* boundary condition */
			if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+2) < hdmin[y][(jp_y+1)] || (d+2) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
			 * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y+1];  /* d index for state y */
			
			if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
			  iescore = cm->iesc[y][(dsq[i-1]*cm->abc->K+dsq[j+1])];
			else
			  iescore = iDegeneratePairScore(cm->abc, cm->iesc[y], dsq[i-1], dsq[j+1]);
			beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
									    + cm->itsc[y][voffset] + iescore));
			break;

		      case ML_st:
		      case IL_st: 
			if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
			if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+1) < hdmin[y][jp_y] || (d+1) > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
			 * corresponding to alpha[y][j][d+1] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y];  /* d index for state y */
			if (dsq[i-1] < cm->abc->K) 
			  iescore = cm->iesc[y][dsq[i-1]];
			else
			  iescore = esl_abc_IAvgScore(cm->abc, dsq[i-1], cm->iesc[y]);
			beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
									    + cm->itsc[y][voffset] + iescore));
			break;
		    
		      case MR_st:
		      case IR_st:
			if (j == j0) continue;
			if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+1) < hdmin[y][(jp_y+1)] || (d+1) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
			 * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][(jp_y+1)];  /* d index for state y */
			if (dsq[j+1] < cm->abc->K) 
			  iescore = cm->iesc[y][dsq[j+1]];
			else
			  iescore = esl_abc_IAvgScore(cm->abc, dsq[j+1], cm->iesc[y]);
			/*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
			beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
									    + cm->itsc[y][voffset] + iescore));
			break;

		      case S_st:
		      case E_st:
		      case D_st:
			if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if (d < hdmin[y][jp_y] || d > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			 * corresponding to alpha[y][j][d] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y];  /* d index for state y */
			beta[v][jp_v][dp_v] = ILogsum(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->itsc[y][voffset])); 
			break;

		      default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
		      }/* end switch over states*/
		  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
		if (beta[v][jp_v][dp_v] < -INFTY) beta[v][jp_v][dp_v] = -INFTY;
		}	/* ends else entered for non-BEGL/BEGR states*/	
	    } /* ends loop over d. We know all beta[v][j][d] in this row j*/
	}    /* end loop over jp. We know the beta's for the whole deck.*/
      /* Deal with local alignment end transitions v->EL
       * (EL = deck at M.)
       */
      if (cm->iendsc[v] != -INFTY) {
	for (j = jmin[v]; j <= jmax[v]; j++)
	  {
	    jp_v = j - jmin[v];
	    for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		i = j-d+1;
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		switch (cm->sttype[v]) {
		case MP_st: 
		  if (j == j0 || d == j) continue; /* boundary condition */
		  if (((j+1) > jmax[v]) || ((d+2) > hdmax[v][jp_v])) continue; /*boundary condition*/
		  if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		    iescore = cm->iesc[v][(dsq[i-1]*cm->abc->K+dsq[j+1])];
		  else
		    iescore = iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i-1], dsq[j+1]);
		  beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+2] + cm->iendsc[v] 
								  + iescore));
		break;
	      case ML_st:
	      case IL_st:
		if (d == j) continue;	
		if ((d+1) > hdmax[v][jp_v]) continue; /*boundary condition*/
		if (dsq[i-1] < cm->abc->K) 
		  iescore = cm->iesc[v][dsq[i-1]];
		else
		  iescore = esl_abc_IAvgScore(cm->abc, dsq[i-1], cm->iesc[v]);
		beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v+1] + cm->iendsc[v] 
								+ iescore));
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (((j+1) > jmax[v]) || ((d+1) > hdmax[v][jp_v])) continue; /*boundary condition*/
		if (dsq[j+1] < cm->abc->K) 
		  iescore = cm->iesc[v][dsq[j+1]];
		else
		  iescore = esl_abc_IAvgScore(cm->abc, dsq[j+1], cm->iesc[v]);
		beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+1] + cm->iendsc[v] 
								+ iescore));
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[v][jp_v][dp_v] + cm->iendsc[v] 
								+ iescore));
		break;
	      case B_st:  
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
		/* note that although B is a valid vend for a segment we'd do
                   outside on, B->EL is set to be impossible, by the local alignment
                   config. There's no point in having a B->EL because B is a nonemitter
                   (indeed, it would introduce an alignment ambiguity). The same
		   alignment case is handled by the X->EL transition where X is the
		   parent consensus state (S, MP, ML, or MR) above the B. Thus,
		   this code is relying on the (cm->iendsc[v] != -INFTY) test, above,
		   to make sure the sttype[vend]=B case gets into this switch.
		*/
	      } /* end switch over parent state type v */
	    } /* end inner loop over d */
	} /* end outer loop over jp */
      } /* end conditional section for dealing w/ v->EL local end transitions */

      /* Look at v's parents; if we're reusing memory (! do_full)
       * push the parents that we don't need any more into the pool.
       */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { Ideckpool_push(dpool, beta[y]); beta[y] = NULL; }
	}
      }
    } /* end loop over decks v. */

  /* EPN left the below SRE block out */
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = L; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = ILogsum(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->iel_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE JD, printing beta\n");
    debug_print_alpha_banded_jd(beta, cm, L, jmin, jmax, hdmin, hdmax);
    printf("OUTSIDE JD, done printing beta\n");
  */

  Lp = L - hdmin[0][j0-jmin[0]];
  if(do_check && (!(cm->flags & CMH_LOCAL_END))) 
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    {
      /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
       * using both the Outside (beta) and Inside (alpha) matrices,
       * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
       * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
       *                                                = P(S|M) / P(S|R)  
       * in v,j,d coordinates this is:
       * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
       *                                                = P(S|M) / P(S|R)
       */
	 
      for(n = 0; n < cm->nodes; n++)
	{
	  isc = -INFTY;
	  num_split_states = SplitStatesInNode(cm->ndtype[n]);
	  for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++)
	    {
	      for (j = jmin[v]; j <= jmax[v]; j++)
		{
		  jp_v = j - jmin[v];
		  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		    {
		      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		      /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
			printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
		      isc = ILogsum(isc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
		    }
		}
	      }
	  fsc = Scorify(isc);
	  /*printf("checking node: %d | sc: %.6f\n", n, fsc);*/
	  diff = fsc - (Scorify(alpha[0][j0-jmin[0]][Lp]));
	  if(diff > 0.01 || diff < -0.01)
	    {
	      fail_flag = TRUE;
	      printf("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
		     n, fsc, Scorify(alpha[0][(j0-jmin[0])][Lp]), diff);
	    }
	}

    }

  /* IF not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * IF local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CMH_LOCAL_END))
    {
      ireturn_sc = -INFTY;
      v = cm->M-1;
      for (j = jmin[v]; j <= jmax[v]; j++)
	{
	  jp_v = j - jmin[v];
	  assert(hdmin[v][jp_v] == 0);
	  /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	  ireturn_sc = ILogsum(ireturn_sc, (beta[v][jp_v][0]));
	}
      freturn_sc = Scorify(ireturn_sc);
    }
  else /* return_sc = P(S|M) / P(S|R) from Inside() */
    {
      freturn_sc = Scorify(alpha[0][(j0-jmin[0])][Lp]);
    }

  /* If the caller doesn't want the beta matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (beta[v] != NULL) { Ideckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
      Ideckpool_push(dpool, beta[cm->M]);
      beta[cm->M] = NULL; 
    }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the alpha matrix, free it 
   * Else, pass it back to him.
   * EPN - if we free the alpha and beta matrix the deck pool has all the 
   *       decks from alpha and beta, not sure if this is desirable.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ideckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { Ideckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    int   **a;
    while (Ideckpool_pop(dpool, &a)) Ifree_vjd_deck(a, i0, j0);
    Ideckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }
  free(touch);

  if(fail_flag) cm_Fail("Not all nodes passed posterior check.");

  if(!(cm->flags & CMH_LOCAL_END))
    ESL_DPRINTF1(("\tIOutside_b_jd_me() sc : %f\n", freturn_sc));
  else
    ESL_DPRINTF1(("\tIOutside_b_jd_me() sc : %f (LOCAL mode; sc is from Inside)\n", freturn_sc));

  return freturn_sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/***************************************************************/
/* Function: CMPosterior_b_jd_me()
 *           
 * EPN 05.27.06 based on IHH's P7EmitterPosterior() from HMMER's postprob.c
 *
 * Purpose:  Combines banded Inside and Outside cubes into a posterior
 *           probability cube. The Inside and Outside cubes are banded
 *           in both the j and d dimensions, any cells outside of
 *           bands do not exist in memory
 *           The entry in post[v][jp_v][dp_v] is the log of the
 *           posterior probability of a parse subtree rooted at v 
 *           emitting the subsequence i..j (i=j-d+1). Where j = jp_v + jmin[v],
 *           and d = dp_v + hdmin[v][jp_v].
 *           The caller must allocate space for the cube, although the
 *           beta matrix from Outside can be used instead (overwriting it will not
 *           compromise the algorithm).
 *           
 * Args:     L        - length of sequence
 *           cm       - the model
 *           alpha    - pre-calculated Inside matrix 
 *           beta     - pre-calculated Outside matrix
 *           post     - pre-allocated dynamic programming cube
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * Return:   void
 */
void
CMPosterior_b_jd_me(int L, CM_t *cm, float ***alpha, float ****ret_alpha,
		    float ***beta, float ****ret_beta, float ***post,
		    float ****ret_post, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int   v, j, d;
  float sc;
  int      jp_v; /* j index for state v in alpha/beta with mem eff bands */
  int      dp_v; /* d index for state v in alpha/beta w/mem eff bands */
  int      Lp;
  float  **end;         /* used for freeing alpha b/c we re-use the end deck. */
  
  Lp = L - hdmin[0][L-jmin[0]];
  sc = alpha[0][L-jmin[0]][Lp];
  
  /* If local ends are on, start with the EL state (cm->M), otherwise
   * its not a valid deck.
   */
  if (cm->flags & CMH_LOCAL_END)
    {
      for(j = 0; j <= L; j++) 
	for (d = 0; d <= j; d++)
	  {
	    post[cm->M][j][d] = alpha[cm->M][j][d] + beta[cm->M][j][d] - sc;
	    /*printf("v: %3d | j: %3d | d: %3d | alpha : %5.2f | beta : %5.2f\n", cm->M, j, d, alpha[cm->M][j][d], beta[cm->M][j][d]);*/
	    /*printf("post[%d][%d][%d]: %f\n", cm->M, j, d, post[cm->M][j][d]);*/
	  }  
    }

  for (v = (cm->M-1); v >= 0; v--) 
    for (j = jmin[v]; j <= jmax[v]; j++) 
      {
	jp_v = j - jmin[v];
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	{
	  dp_v = d - hdmin[v][jp_v];
	  /*printf("v: %3d | jp_v: %3d | dp_v: %3d | alpha: %5.2f | beta: %5.2f\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	  post[v][jp_v][dp_v] = alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v] - sc;
	}  
      }

  /* If the caller doesn't want the matrix, free it and free the decks in the pool
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { free_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { free_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  }
  else *ret_alpha = alpha;

  /* If the caller doesn't want the beta matrix, free it along with the decks.
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M); v++) 
      if (beta[v] != NULL) { free_vjd_deck(beta[v], 1, L); beta[v] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the post matrix, free it, though
   * it would be *stupid* for the caller not to want it in current implementation.
   */
  if (ret_post == NULL) {
    for (v = 0; v <= (cm->M); v++) 
      if (post[v] != NULL) { free_vjd_deck(post[v], 1, L); post[v] = NULL; }
    free(post);
  } else *ret_post = post;
}

void
ICMPosterior_b_jd_me(int L, CM_t *cm, int   ***alpha, int   ****ret_alpha,
		    int   ***beta, int   ****ret_beta, int   ***post,
		    int ****ret_post, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int   v, j, d;
  int   sc;
  int      jp_v; /* j index for state v in alpha/beta with mem eff bands */
  int      dp_v; /* d index for state v in alpha/beta w/mem eff bands */
  int      Lp;
  int    **end;         /* used for freeing alpha b/c we re-use the end deck. */
  
  Lp = L - hdmin[0][L-jmin[0]];
  sc = alpha[0][L-jmin[0]][Lp];
  
  /* If local ends are on, start with the EL state (cm->M), otherwise
   * it's not a valid deck.
   */
  if (cm->flags & CMH_LOCAL_END)
    {
      for(j = 0; j <= L; j++) 
	for (d = 0; d <= j; d++)
	  {
	    post[cm->M][j][d] = alpha[cm->M][j][d] + beta[cm->M][j][d] - sc;
	    /*printf("v: %3d | j: %3d | d: %3d | alpha : %5.2f | beta : %d\n", cm->M, j, d, alpha[cm->M][j][d], beta[cm->M][j][d]);
	      printf("post[%d][%d][%d]: %d\n", cm->M, j, d, post[cm->M][j][d]);*/
	  }  
    }

  for (v = (cm->M-1); v >= 0; v--) 
    for (j = jmin[v]; j <= jmax[v]; j++) 
      {
	jp_v = j - jmin[v];
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	{
	  dp_v = d - hdmin[v][jp_v];
	  post[v][jp_v][dp_v] = alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v] - sc;
	  /*printf("v: %3d | jp_v: %3d | dp_v: %3d | alpha: %10d | beta: %10d | post: %10d\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v], post[v][jp_v][dp_v]);*/
	  /* if(Score2Prob(post[v][jp_v][dp_v], 1.) > 1.03) printf("%f v: %d jp_v: %d dp_v: %d isc: %d\n", Score2Prob(post[v][jp_v][dp_v], 1.), v, jp_v, dp_v, post[v][jp_v][dp_v]); */
	  /* assert(Score2Prob(post[v][jp_v][dp_v], 1.) <= 1.001); */
	}  
      }

  /* If the caller doesn't want the matrix, free it and free the decks in the pool
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ifree_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { Ifree_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  }
  else *ret_alpha = alpha;

  /* If the caller doesn't want the beta matrix, free it along with the decks.
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M); v++) 
      if (beta[v] != NULL) { Ifree_vjd_deck(beta[v], 1, L); beta[v] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the post matrix, free it, though
   * it would be *stupid* for the caller not to want it in current implementation.
   */
  if (ret_post == NULL) {
    for (v = 0; v <= (cm->M); v++) 
      if (post[v] != NULL) { Ifree_vjd_deck(post[v], 1, L); post[v] = NULL; }
    free(post);
  } else *ret_post = post;
}

/*################################################################*/
/* Integer versions of functions necessary for implementation
 * of cm_postprob.c functions using scaled integer log odds scores. Copied
 * and minimally changed from smallcyk.c.
 */

/*################################################################*/
/* Functions: Ideckpool_*()
 * Date:      EPN, Fri Dec 29 06:50:51 2006
 *            adapted from deckpool_*() functions in smallcyk.c by:
 *            SRE, Wed Aug  2 10:43:17 2000 [St. Louis]
 *
 * Purpose:   Implementation of a pushdown stack for storing decks 
 *            of the inside or outside dynamic programming matrices
 *            with ints (not floats - this is the *I*deckpool part),
 *            with the usual _create, _push, _pop, and _free API. 
 *            
 *            The deck pool allows us to efficiently reuse memory,
 *            so long as our DP algorithms step through the decks
 *            as their outermost loop.
 *            
 *            Works for either coordinate system (vjd or vji) 
 *            and subseq variants, because it's simply managing
 *            a deck as a float **.
 */
Ideckpool_t *
Ideckpool_create(void)
{
  int status;
  Ideckpool_t *dpool;

  ESL_ALLOC(dpool, sizeof(Ideckpool_t));
  dpool->block  = 10;		/* configurable if you want */
  ESL_ALLOC(dpool->pool, sizeof(int **) * dpool->block);
  dpool->nalloc = dpool->block;
  dpool->n      = 0;
  return dpool;
  
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
void 
Ideckpool_push(Ideckpool_t *dpool, int **deck)
{
  int status;
  void *tmp;
  if (dpool->n == dpool->nalloc) {
    dpool->nalloc += dpool->block;
    ESL_RALLOC(dpool->pool, tmp, sizeof(int **) * dpool->nalloc);
  }
  dpool->pool[dpool->n] = deck;
  dpool->n++;
  ESL_DPRINTF3(("Ideckpool_push\n"));
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
}
int
Ideckpool_pop(Ideckpool_t *d, int ***ret_deck)
{
  if (d->n == 0) { *ret_deck = NULL; return 0;}
  d->n--;
  *ret_deck = d->pool[d->n];
  ESL_DPRINTF3(("Ideckpool_pop\n"));
  return 1;
}
void
Ideckpool_free(Ideckpool_t *d)
{
  free(d->pool);
  free(d);
}
/*================================================================*/

/*################################################################*/
/* Functions: I*_vjd_*
 * Date:     EPN, Fri Dec 29 06:55:29 2006
 *           adapted from *_vjd_*() functions in smallcyk.c by:
 *           SRE, Sat Aug 12 16:27:37 2000 [Titusville]
 * 
 * 
 * Purpose:  Allocation and freeing of 3D matrices and 2D decks
 *           of ints (not floats - this is the *I*_vjd_ part)
 *           in the vjd coord system. These can be called on
 *           subsequences i..j, not just the full sequence 1..L,
 *           so they need i,j... if you're doing the full sequence
 *           just pass 1,L.
 *           
 *           We don't need to deal with shadow matrices differently
 *           between float and int log odds scores, we can use
 *           the original functions in smallcyk.c
 */
int **
Ialloc_vjd_deck(int L, int i, int j)
{
  int status;
  int **a;
  int     jp;

  ESL_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  ESL_ALLOC(a, sizeof(int *) * (L+1)); /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0;   jp < i-1;    jp++) a[jp]     = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp]     = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) ESL_ALLOC(a[jp+i-1], sizeof(int) * (jp+1));
  return a;

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
int
Isize_vjd_deck(int L, int i, int j)
{
  float Mb;
  int   jp;
  Mb = (float) (sizeof(int *) * (L+1));
  for (jp = 0; jp <= j-i+1; jp++)
    Mb += (float) (sizeof(int) * (jp+1));
  return (Mb / 1000000.);
}

void
Ifree_vjd_deck(int **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp+i-1]);
  free(a);
}
void
Ifree_vjd_matrix(int ***a, int M, int i, int j)
{
  int v;
  for (v = 0; v <= M; v++)
    if (a[v] != NULL)		/* protect against double free's of reused decks (ends) */
      { Ifree_vjd_deck(a[v], i, j); a[v] = NULL; }
  free(a);
}
/*================================================================*/

/* A single vjd deck allocation function: Ialloc_jdbanded_vjd_deck()
 * is needed to support scaled int log odds scores in HMM banded (in
 * the j and d dimensions) Inside and Outside functions. This is copied
 * from hbandcyk.c and modified to handle int log odds scores, not floats.
 */
int **
Ialloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  int     status;
  int   **a;
  int     jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;
  /*printf("in alloc JD banded vjd deck, L : %d, i : %d, j : %d, jmin : %d, jmax : %d\n", L, i, j, jmin, jmax);*/

  if(j < jmin || i > jmax)
    cm_Fail("ERROR called alloc_jdbanded_vjd_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);

  ESL_ALLOC(a, sizeof(int *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0; jp <= L;     jp++) a[jp]     = NULL;

  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jfirst: %d | jlast: %d\n", jfirst, jlast);
      printf("jp: %d | max : %d\n", jp, (jlast)); 
      printf("hdmax[%d]: %d\n", (jp-jmin), hdmax[jp-jmin]);
      */
      if(hdmax[jp-jmin] > (jp+1))
	{
	  /* Based on my current understanding this shouldn't happen, it means there's a valid d
	   * in the hd band that is invalid because its > j. I check, or ensure, that this
	   * doesn't happen when I'm constructing the d bands.
	   */
	  cm_Fail("jd banded error 0.\n");
	}
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*a is offset only the first (jlast-jfirst+1) elements will be non-NULL*/
      ESL_ALLOC(a[jp-jfirst], sizeof(int) * bw);
      /*printf("\tallocated a[%d] | bw: %d\n", (jp-jfirst), bw);*/
    }
  return a;

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}


/*
 * Function: ParsetreeSampleFromIInside()
 * Incept:   EPN, Thu Sep  6 09:58:26 2007
 *          
 * Purpose:  Sample a parsetree from an non-banded integer Inside matrix.
 *           
 * Args:     r        - source of randomness
 *           cm       - the model
 *           dsq      - digitized sequence
 *           L        - length of dsq, alpha *must* go from 1..L
 *           alpha    - pre-calculated Inside matrix (ints)
 *           ret_tr   - ptr to parsetree we'll return (*must* be non-NULL)
 *           ret_alpha- pass NULL to free input alpha, otherwise it's passed back here
 *
 * Return:   score of sampled parsetree; dies immediately with cm_Fail if an error occurs.
 */
float
ParsetreeSampleFromIInside(ESL_RANDOMNESS *r, CM_t *cm, ESL_DSQ *dsq, int L, int ***alpha, Parsetree_t **ret_tr, int ****ret_alpha)
{
  int          status;             /* easel status code */
  int          v, y, z, b;         /* state indices */
  int          yoffset;            /* transition offset in a states transition vector */
  int          i, j;               /* sequence position indices */
  int          d;                  /* j - i + 1; the current subseq length */
  int          k;                  /* right subseq fragment length for bifurcs */
  int          nd;                 /* node index */
  int          bifparent;          /* for connecting bifurcs */
  Parsetree_t *tr;                 /* trace we're building */
  ESL_STACK   *pda;                /* the stack */
  float        pvec[MAXCONNECT+1]; /* prob vector of possible paths to take, (max num children + 1 for possibility of EL) */
  float       *bifvec;             /* pvec for choosing transition out of BIF_B states */
  float       *rootvec;            /* pvec for choosing transition out of ROOT_S if local begins are on */
  float        maxsc;              /* max score in our vector of scores of possible subparses */
  int          el_is_possible;     /* TRUE if we can jump to EL from current state (and we're in local mode) FALSE if not */
  int          ntrans;             /* number of transitions for current state */
  int          isc = 0;            /* score of the parsetree we're sampling */
  int        **end;                /* used for freeing alpha b/c we re-use the end deck. */

  /* contract check */
  if(ret_tr == NULL) cm_Fail("ParsetreeSampleFromIInside(), ret_tr is NULL.");
  if(r      == NULL) cm_Fail("ParsetreeSampleFromIInside(), source of randomness r is NULL.");
  
  /* initialize pvec */
  esl_vec_FSet(pvec, (MAXCONNECT+1), 0.);

  /* Create a parse tree structure and initialize it by adding the root state.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Stochastically traceback through the Inside matrix 
   * this section of code is stolen and adapted from smallcyk.c:insideT() 
   */
  pda = esl_stack_ICreate();
  v = 0;

  j = d = L;
  i = 1;
  isc = 0;
  while (1) {
    if (cm->sttype[v] == B_st) {
      y = cm->cfirst[v];
      z = cm->cnum[v];

      ESL_ALLOC(bifvec, sizeof(float) * (d+1));
      /* set bifvec[] as (float-ized) log odds scores for each valid left fragment length */
      for(k = 0; k <= d; k++) 
	bifvec[k] = Scorify(alpha[y][j-k][d-k] + alpha[z][j][k]);
      maxsc = esl_vec_FMax(bifvec, (d+1));
      esl_vec_FIncrement(bifvec, (d+1), (-1. * maxsc));
      esl_vec_FScale(bifvec, (d+1), log(2));
      esl_vec_FLogNorm(bifvec, (d+1));
      k = esl_rnd_FChoose(r, bifvec, (d+1));
      free(bifvec);

      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */

      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);

      v = y;
    } else {
      if((v > 0) || (! (cm->flags & CMH_LOCAL_BEGIN))) /* ROOT_S with local begins on is a special case that we handle below */
	{ 
	  /* choose which transition we take */
	  esl_vec_FSet(pvec, (MAXCONNECT+1), IMPOSSIBLE); /* not really necessary */
	  isc += get_iemission_score(cm, dsq, v, i, j); 
	  
	  /* set pvec[] as (float-ized) log odds scores for each child we can transit to, 
	   * plus a local end (if possible) */
	  ntrans = cm->cnum[v];
	  el_is_possible = FALSE;
	  if((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) { 
	    el_is_possible = TRUE; 
	    ntrans++; 
	  }
	  for(yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = yoffset + cm->cfirst[v];
	    pvec[yoffset] = Scorify(cm->itsc[v][yoffset] + 
				    alpha[y][j - StateRightDelta(cm->sttype[v])][d - StateDelta(cm->sttype[v])]);
	  }
	  if(el_is_possible) pvec[cm->cnum[v]] = Scorify(cm->iendsc[v] + 
						 alpha[cm->M][j][d]); /* EL is silent when we transition into it from non-EL */
	  /* note: we can treat the log odds scores as log probs, because
	   * the log probability of the null model is the same for each,
	   * so essentially we've divided each score by the same constant, so 
	   * the *relative* proportion of the log odds scores is the
	   * same as the relative proportion of the log probabilities (seq | model) */
	  
	  maxsc = esl_vec_FMax(pvec, ntrans);
	  esl_vec_FIncrement(pvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (pvec, ntrans, log(2));
	  esl_vec_FLogNorm(pvec, ntrans);
	  yoffset = esl_rnd_FChoose(r, pvec, ntrans);
	  if(yoffset < cm->cnum[v]) isc += cm->itsc[v][yoffset]; 
	  else {
	    isc += cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
	    yoffset = USED_EL; /* we chose EL */
	  }
	}
      else /* v == 0 && (cm->flags && CMH_LOCAL_BEGIN) ( local begins are on )*/
	{
	  ntrans = cm->M; /* pretend all states are possible to begin into, but they're not as some will remain IMPOSSIBLE */
	  ESL_ALLOC(rootvec, sizeof(float) * (ntrans));
	  esl_vec_FSet(rootvec, ntrans, IMPOSSIBLE);
	  rootvec[cm->nodemap[1]] = Scorify(cm->ibeginsc[cm->nodemap[1]] + 
					    alpha[cm->nodemap[1]][j][d]); /* ROOT_S is silent */
	  for (nd = 2; nd < cm->nodes; nd++) {
	    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
		cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
	      {
		rootvec[cm->nodemap[nd]] = Scorify(cm->ibeginsc[cm->nodemap[nd]] + 
						   alpha[cm->nodemap[nd]][j][d]); /* ROOT_S is silent */
	      }
	  }
	  /* this block is shared with v > 0 block, but we repeat it here so we don't need another if statement */
	  maxsc = esl_vec_FMax(rootvec, ntrans);
	  esl_vec_FIncrement(rootvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (rootvec, ntrans, log(2));
	  esl_vec_FLogNorm(rootvec, ntrans);
	  b = esl_rnd_FChoose(r, rootvec, ntrans);
	  /* end of similar block with v > 0 */
	  isc += cm->ibeginsc[b];
	  yoffset = USED_LOCAL_BEGIN; 
	  free(rootvec); /* we will not need this again */
	}

      /*printf("v : %d | r : %d | z : %d | 1 : %d | \n", v, r, z, 1);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */

  /* If the caller doesn't want the alpha matrix, free it
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= cm->M; v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ifree_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
	}
    if (end != NULL) { Ifree_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  *ret_tr = tr; /* contract checked ret_tr was non-NULL */
  return Scorify(isc);

 ERROR:
  cm_Fail("memory error.");
  return 0.; /* NEVERREACHED */
}

/*
 * Function: ParsetreeSampleFromIInside_b_jd_me()
 * Incept:   EPN, Fri Sep  7 11:02:15 2007
 *          
 * Purpose:  Sample a parsetree from a integer Inside matrix banded in the j and d dimensions.
 *           
 * Args:     r        - source of randomness
 *           cm       - the model
 *           dsq      - digitized sequence
 *           L        - length of dsq, alpha *must* go from 1..L
 *           alpha    - pre-calculated Inside matrix (ints)
 *           cp9b     - CP9Bands data structure giving bands on j and d dimensions
 *           ret_tr   - ptr to parsetree we'll return (*must* be non-NULL)
 *           ret_alpha- pass NULL to free input alpha, otherwise it's passed back here
 *
 * Return:   score of sampled parsetree; dies immediately with cm_Fail if an error occurs.
 */
float
ParsetreeSampleFromIInside_b_jd_me(ESL_RANDOMNESS *r, CM_t *cm, ESL_DSQ *dsq, int L, int ***alpha, CP9Bands_t *cp9b, Parsetree_t **ret_tr, int ****ret_alpha)
{
  int          status;             /* easel status code */
  int          v, y, z, b;         /* state indices */
  int          yoffset;            /* transition offset in a states transition vector */
  int          i, j;               /* sequence position indices */
  int          jp_v, jp_y, jp_z;   /* positions, offset inside j band */
  int          kmin, kmax;         /* min/max k in current d band */
  int          d;                  /* j - i + 1; the current subseq length */
  int          dp_v, dp_y;         /* length, offset inside a d band */
  int          k;                  /* right subseq fragment length for bifurcs */
  int          kp_z;               /* right fragment length, offset inside a d band */
  int          nd;                 /* node index */
  int          bifparent;          /* for connecting bifurcs */
  Parsetree_t *tr;                 /* trace we're building */
  ESL_STACK   *pda;                /* the stack */
  float        pvec[MAXCONNECT+1]; /* prob vector of possible paths to take, (max num children + 1 for possibility of EL) */
  float       *bifvec;             /* pvec for choosing transition out of BIF_B states */
  float       *rootvec;            /* pvec for choosing transition out of ROOT_S if local begins are on */
  float        maxsc;              /* max score in our vector of scores of possible subparses */
  int          el_is_possible;     /* TRUE if we can jump to EL from current state (and we're in local mode) FALSE if not */
  int          ntrans;             /* number of transitions for current state */
  int          isc = 0;            /* score of the parsetree we're sampling */
  int        **end;                /* used for freeing alpha b/c we re-use the end deck. */
  int          seen_valid;         /* for checking we have at least one valid path to take  */
  int          sd;                 /* state delta for current state, residues emitted left + residues emitted right */
  int          sdr;                /* state right delta for current state, residues emitted right */
  int         *jmin;               /* ptr to CP9Bands cp9b's j min band, for convenience */
  int         *jmax;               /* ptr to CP9Bands cp9b's j max band, for convenience */
  int        **hdmin;              /* ptr to CP9Bands cp9b's d min band, for convenience */
  int        **hdmax;              /* ptr to CP9Bands cp9b's d max band, for convenience */

  /* for convenience */
  jmin  = cp9b->jmin;
  jmax  = cp9b->jmax;
  hdmin = cp9b->hdmin;
  hdmax = cp9b->hdmax;

  /* contract check */
  if(ret_tr == NULL) cm_Fail("ParsetreeSampleFromIInside_b_jd_me(), ret_tr is NULL.");
  if(r      == NULL) cm_Fail("ParsetreeSampleFromIInside_b_jd_me(), source of randomness r is NULL.");
  
  /* initialize pvec */
  esl_vec_FSet(pvec, (MAXCONNECT+1), 0.);

  /* Create a parse tree structure and initialize it by adding the root state.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Stochastically traceback through the Inside matrix 
   * this section of code is stolen and adapted from hbandcyk.c:insideT_b_jd_me() 
   */
  pda = esl_stack_ICreate();
  v = 0;

  j = d = L;
  i = 1;

  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];
  isc  = 0;
  while (1) {
    if(cm->sttype[v] != EL_st && d > hdmax[v][jp_v])
      cm_Fail("ERROR in ParsetreeSampleFromIInside_b_jd_me(). d : %d > hdmax[%d] (%d)\n", d, v, hdmax[v][jp_v]);
    if(cm->sttype[v] != EL_st && d < hdmin[v][jp_v])
      cm_Fail("ERROR in ParsetreeSampleFromIInside_b_jd_me(). d : %d < hdmin[%d] (%d)\n", d, v, hdmin[v][jp_v]);

    if (cm->sttype[v] == B_st) {
      y = cm->cfirst[v];
      z = cm->cnum[v];
      jp_z = j-jmin[z];
      k = kp_z + hdmin[z][jp_z];  /* k = offset len of right fragment */

      ESL_ALLOC(bifvec, sizeof(float) * (d+1));
      /* set bifvec[] as (float-ized) log odds scores for each valid left fragment length,
       * we have to be careful to check that the corresponding alpha cell for each length is valid
       */
      esl_vec_FSet(bifvec, (d+1), IMPOSSIBLE); /* only valid d's will be reset to a non-IMPOSSIBLE score */

      /* This search for valid k's is complex, and uncommented. It was taken from
       * hbandcyk.c:inside_b_jd_me(), the B_st case. The code there is commented somewhat
       * extensively. I'm pretty sure this is the most efficient (or at least close to it) 
       * way to find the valid cells in the DP matrix we're looking for. 
       */
      jp_v = j - jmin[v];
      jp_y = j - jmin[y];
      jp_z = j - jmin[z];
      if(j < jmin[v] || j > jmax[v]) cm_Fail("ParsetreeSampleFromIInside_b_jd_me() B_st v: %d j: %d outside band jmin: %d jmax: %d\n", v, j, jmin[v], jmax[v]);
      if(d < hdmin[v][jp_v] || d > hdmax[v][jp_v]) cm_Fail("ParsetreeSampleFromIInside_b_jd_me() B_st v: %d j: %d d: %d outside band dmin: %d dmax: %d\n", v, j, d, hdmin[v][jp_v], hdmax[v][jp_v]);
      seen_valid = FALSE;
      kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
      kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
      for(k = kmin; k <= kmax; k++)
	{
	  if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
	    {
	      kp_z = k-hdmin[z][jp_z];
	      dp_y = d-hdmin[y][jp_y-k];
	      bifvec[k] = alpha[y][jp_y-k][dp_y-k] + alpha[z][jp_z][kp_z]; 
	      seen_valid = TRUE;
	    }
	}
      if(!seen_valid) cm_Fail("ParsetreeSampleFromIInside_b_jd_me() number of valid transitions (for a B_st) is 0. You thought this was impossible.");
      maxsc = esl_vec_FMax(bifvec, (d+1));
      esl_vec_FIncrement(bifvec, (d+1), (-1. * maxsc));
      esl_vec_FScale(bifvec, (d+1), log(2));
      esl_vec_FLogNorm(bifvec, (d+1));
      k = esl_rnd_FChoose(r, bifvec, (d+1));
      free(bifvec);

      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */

      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);

      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else {
      if((v > 0) || (! (cm->flags & CMH_LOCAL_BEGIN))) /* ROOT_S with local begins on is a special case that we handle below */
	{ 
	  /* Choose which transition we take.
	   * Set pvec[] as (float-ized) log odds scores for each child we can transit to, 
	   * plus a local end (if possible). We only want to look at valid transitions, that
	   * is, those that do not violate the bands (correspond to accessing cells that actually
	   * exist in the DP matrix). 
	   */
	  seen_valid = FALSE;
	  esl_vec_FSet(pvec, (MAXCONNECT+1), IMPOSSIBLE); /* only transitions that correspond to valid cells will be reset to a non-IMPOSSIBLE score */
	  isc += get_iemission_score(cm, dsq, v, i, j); 
	  sdr = StateRightDelta(cm->sttype[v]);
	  sd  = StateDelta(cm->sttype[v]);
	  for(yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
	    {
	      y = yoffset + cm->cfirst[v];
	      if((j - sdr) >= jmin[y] && (j - sdr) <= jmax[y]) 
		{ /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if((d - sd) >= hdmin[y][jp_y-sdr] && (d - sd) <= hdmax[y][jp_y-sdr])
		    {
		      dp_y = d - hdmin[y][(jp_y - sdr)];  /* d index for state y 
							     in alpha w/mem eff bands */
		      /* if we get here alpha[y][jp_y-sdr][dp_y-sd] is a valid alpha cell
		       * corresponding to alpha[y][j-sdr][d-sd] in the platonic matrix.
		       */
		      pvec[yoffset] = Scorify(cm->itsc[v][yoffset] + 
					      alpha[y][jp_y - sdr][dp_y - sd]);
		      seen_valid = TRUE;
		    }
		}		
	    }
	  if(!seen_valid) {
	    cm_Fail("ParsetreeSampleFromIInside_b_jd_me() number of valid transitions is 0. You thought this was impossible.");
	  }
	  if((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) 
	    el_is_possible = TRUE; 
	  else 
	    el_is_possible = FALSE;
	  if(el_is_possible) pvec[cm->cnum[v]] = Scorify(cm->iendsc[v] + 
							 alpha[cm->M][j][d]); /* EL is silent when we transition into it from non-EL */
	  ntrans = cm->cnum[v] + el_is_possible;
	  maxsc = esl_vec_FMax(pvec, ntrans);
	  esl_vec_FIncrement(pvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (pvec, ntrans, log(2));
	  esl_vec_FLogNorm(pvec, ntrans);
	  yoffset = esl_rnd_FChoose(r, pvec, ntrans);
	  if(yoffset < cm->cnum[v]) isc += cm->itsc[v][yoffset]; 
	  else {
	    isc += cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
	    yoffset = USED_EL; /* we chose EL */
	  }
	}
      else /* v == 0 && (cm->flags && CMH_LOCAL_BEGIN) ( local begins are on )*/
	{
	  seen_valid = FALSE;
	  ntrans = cm->M; /* pretend all states are possible to begin into, but they're not as some will remain IMPOSSIBLE */
	  ESL_ALLOC(rootvec, sizeof(float) * (ntrans));
	  esl_vec_FSet(rootvec, ntrans, IMPOSSIBLE);

	  /* Set all the legal states that we can local begin into to appropriate scores.
	   * Only states y that have a non-zero cm->beginsc[y] AND have alpha[y][j][d]
	   * within their bands are legal.
	   */
	  for (nd = 1; nd < cm->nodes; nd++) {
	    if ((nd == 1) || /* we can transit into node 1 no matter what */
		(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
		 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd))
	      {
		y = cm->nodemap[nd];
		if(j >= jmin[y] && j <= jmax[y]) 
		  { /* enforces j is valid for state y */
		    jp_y = j - jmin[y];
		    if(d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y])
		      {
			dp_y = d - hdmin[y][jp_y];
			rootvec[y] = Scorify(cm->ibeginsc[y] +
					     alpha[y][jp_y][dp_y]); /* ROOT_S is silent */
			seen_valid = TRUE;
		      }
		  }
	      }
	  }
	  if(!seen_valid) cm_Fail("ParsetreeSampleFromIInside_b_jd_me() number of valid transitions (from ROOT_S!) is 0. You thought this was impossible.");
	  /* this block is shared with v > 0 block, but we repeat it here so we don't need another if statement */
	  maxsc = esl_vec_FMax(rootvec, ntrans);
	  esl_vec_FIncrement(rootvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (rootvec, ntrans, log(2));
	  esl_vec_FLogNorm(rootvec, ntrans);
	  b = esl_rnd_FChoose(r, rootvec, ntrans);
	  /* end of similar block with v > 0 */
	  isc += cm->ibeginsc[b];
	  yoffset = USED_LOCAL_BEGIN; 
	  free(rootvec); /* we will not need this again */
	}

      /*printf("v : %d | r : %d | z : %d | 1 : %d | \n", v, r, z, 1);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	  jp_v = j;
	  dp_v = d;
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */

  /* If the caller doesn't want the alpha matrix, free it
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= cm->M; v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { Ifree_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
	}
    if (end != NULL) { Ifree_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  *ret_tr = tr; /* contract checked ret_tr was non-NULL */
  return Scorify(isc);

 ERROR:
  cm_Fail("memory error.");
  return 0.; /* NEVERREACHED */
}


/*
 * Function: GetIEmissionScore()
 * Incept:   EPN, Thu Sep  6 13:35:36 2007
 *          
 * Purpose:  Given a CM, dsq, state index and coordinates return the integer emission
 *           score.
 *           
 * Args:     cm       - the model
 *           dsq      - digitized sequence
 *           v        - state index
 *           i        - dsq index for first position of subseq for subtree at v
 *           j        - dsq index for last position of subseq for subtree at v
 *
 * Return:   integer emission score, 0 if state is non-emitter.
 */
int
get_iemission_score(CM_t *cm, ESL_DSQ *dsq, int v, int i, int j)
{
  if     (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) return cm->ioesc[v][dsq[i]];
  else if(cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) return cm->ioesc[v][dsq[j]];
  else if(cm->sttype[v] == MP_st)                           return cm->ioesc[v][dsq[i]*cm->abc->Kp+dsq[j]];
  else return 0;
}


/*
 * Function: ParsetreeSampleFromFInside()
 * Incept:   EPN, Thu Nov 15 16:45:32 2007
 *          
 * Purpose:  Sample a parsetree from an non-banded float Inside matrix.
 *           
 * Args:     r        - source of randomness
 *           cm       - the model
 *           dsq      - digitized sequence
 *           L        - length of dsq, alpha *must* go from 1..L
 *           alpha    - pre-calculated Inside matrix (floats)
 *           ret_tr   - ptr to parsetree we'll return (*must* be non-NULL)
 *           ret_alpha- pass NULL to free input alpha, otherwise it's passed back here
 *
 * Return:   score of sampled parsetree; dies immediately with cm_Fail if an error occurs.
 */
float
ParsetreeSampleFromFInside(ESL_RANDOMNESS *r, CM_t *cm, ESL_DSQ *dsq, int L, float ***alpha, Parsetree_t **ret_tr, float ****ret_alpha)
{
  int          status;             /* easel status code */
  int          v, y, z, b;         /* state indices */
  int          yoffset;            /* transition offset in a states transition vector */
  int          i, j;               /* sequence position indices */
  int          d;                  /* j - i + 1; the current subseq length */
  int          k;                  /* right subseq fragment length for bifurcs */
  int          nd;                 /* node index */
  int          bifparent;          /* for connecting bifurcs */
  Parsetree_t *tr;                 /* trace we're building */
  ESL_STACK   *pda;                /* the stack */
  float        pvec[MAXCONNECT+1]; /* prob vector of possible paths to take, (max num children + 1 for possibility of EL) */
  float       *bifvec;             /* pvec for choosing transition out of BIF_B states */
  float       *rootvec;            /* pvec for choosing transition out of ROOT_S if local begins are on */
  float        maxsc;              /* max score in our vector of scores of possible subparses */
  int          el_is_possible;     /* TRUE if we can jump to EL from current state (and we're in local mode) FALSE if not */
  int          ntrans;             /* number of transitions for current state */
  int          fsc = 0.;           /* score of the parsetree we're sampling */

  /* contract check */
  if(ret_tr == NULL) cm_Fail("ParsetreeSampleFromFInside(), ret_tr is NULL.");
  if(r      == NULL) cm_Fail("ParsetreeSampleFromFInside(), source of randomness r is NULL.");
  
  /* initialize pvec */
  esl_vec_FSet(pvec, (MAXCONNECT+1), 0.);

  /* Create a parse tree structure and initialize it by adding the root state.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Stochastically traceback through the Inside matrix 
   * this section of code is stolen and adapted from smallcyk.c:insideT() 
   */
  pda = esl_stack_ICreate();
  v = 0;

  j = d = L;
  i = 1;
  fsc = 0.;
  while (1) {
    if (cm->sttype[v] == B_st) {
      y = cm->cfirst[v];
      z = cm->cnum[v];

      ESL_ALLOC(bifvec, sizeof(float) * (d+1));
      /* set bifvec[] as (float-ized) log odds scores for each valid left fragment length */
      for(k = 0; k <= d; k++) 
	bifvec[k] = alpha[y][j-k][d-k] + alpha[z][j][k];
      maxsc = esl_vec_FMax(bifvec, (d+1));
      esl_vec_FIncrement(bifvec, (d+1), (-1. * maxsc));
      esl_vec_FScale(bifvec, (d+1), log(2));
      esl_vec_FLogNorm(bifvec, (d+1));
      k = esl_rnd_FChoose(r, bifvec, (d+1));
      free(bifvec);

      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */

      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);

      v = y;
    } else {
      if((v > 0) || (! (cm->flags & CMH_LOCAL_BEGIN))) /* ROOT_S with local begins on is a special case that we handle below */
	{ 
	  /* choose which transition we take */
	  esl_vec_FSet(pvec, (MAXCONNECT+1), IMPOSSIBLE); /* not really necessary */
	  fsc += get_femission_score(cm, dsq, v, i, j); 
	  
	  /* set pvec[] as (float-ized) log odds scores for each child we can transit to, 
	   * plus a local end (if possible) */
	  ntrans = cm->cnum[v];
	  el_is_possible = FALSE;
	  if((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) { 
	    el_is_possible = TRUE; 
	    ntrans++; 
	  }
	  for(yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
	    y = yoffset + cm->cfirst[v];
	    pvec[yoffset] = cm->tsc[v][yoffset] + 
	      alpha[y][j - StateRightDelta(cm->sttype[v])][d - StateDelta(cm->sttype[v])];
	  }
	  if(el_is_possible) pvec[cm->cnum[v]] = cm->endsc[v] + 
			       alpha[cm->M][j][d]; /* EL is silent when we transition into it from non-EL */
	  /* note: we can treat the log odds scores as log probs, because
	   * the log probability of the null model is the same for each,
	   * so essentially we've divided each score by the same constant, so 
	   * the *relative* proportion of the log odds scores is the
	   * same as the relative proportion of the log probabilities (seq | model) */
	  
	  maxsc = esl_vec_FMax(pvec, ntrans);
	  esl_vec_FIncrement(pvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (pvec, ntrans, log(2));
	  esl_vec_FLogNorm(pvec, ntrans);
	  yoffset = esl_rnd_FChoose(r, pvec, ntrans);
	  if(yoffset < cm->cnum[v]) fsc += cm->tsc[v][yoffset]; 
	  else {
	    fsc += cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
	    yoffset = USED_EL; /* we chose EL */
	  }
	}
      else /* v == 0 && (cm->flags && CMH_LOCAL_BEGIN) ( local begins are on )*/
	{
	  ntrans = cm->M; /* pretend all states are possible to begin into, but they're not as some will remain IMPOSSIBLE */
	  ESL_ALLOC(rootvec, sizeof(float) * (ntrans));
	  esl_vec_FSet(rootvec, ntrans, IMPOSSIBLE);
	  rootvec[cm->nodemap[1]] = cm->beginsc[cm->nodemap[1]] + alpha[cm->nodemap[1]][j][d]; /* ROOT_S is silent */
	  for (nd = 2; nd < cm->nodes; nd++) {
	    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
		cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
	      {
		rootvec[cm->nodemap[nd]] = cm->beginsc[cm->nodemap[nd]] + alpha[cm->nodemap[nd]][j][d]; /* ROOT_S is silent */
	      }
	  }
	  /* this block is shared with v > 0 block, but we repeat it here so we don't need another if statement */
	  maxsc = esl_vec_FMax(rootvec, ntrans);
	  esl_vec_FIncrement(rootvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (rootvec, ntrans, log(2));
	  esl_vec_FLogNorm(rootvec, ntrans);
	  b = esl_rnd_FChoose(r, rootvec, ntrans);
	  /* end of similar block with v > 0 */
	  fsc += cm->beginsc[b];
	  yoffset = USED_LOCAL_BEGIN; 
	  free(rootvec); /* we will not need this again */
	}

      /*printf("v : %d | r : %d | z : %d | 1 : %d | \n", v, r, z, 1);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */

  /* If the caller doesn't want the alpha matrix, free it
   */
  if (ret_alpha == NULL) free_vjd_matrix(alpha, cm->M, 1, L);
  else *ret_alpha = alpha;

  *ret_tr = tr; /* contract checked ret_tr was non-NULL */
  return fsc;

 ERROR:
  cm_Fail("memory error.");
  return 0.; /* NEVERREACHED */
}

/*
 * Function: ParsetreeSampleFromFInside_b_jd_me()
 * Incept:   EPN, Fri Sep  7 11:02:15 2007
 *          
 * Purpose:  Sample a parsetree from a integer Inside matrix banded in the j and d dimensions.
 *           
 * Args:     r        - source of randomness
 *           cm       - the model
 *           dsq      - digitized sequence
 *           L        - length of dsq, alpha *must* go from 1..L
 *           alpha    - pre-calculated Inside matrix (ints)
 *           cp9b     - CP9Bands data structure giving bands on j and d dimensions
 *           ret_tr   - ptr to parsetree we'll return (*must* be non-NULL)
 *           ret_alpha- pass NULL to free input alpha, otherwise it's passed back here
 *
 * Return:   score of sampled parsetree; dies immediately with cm_Fail if an error occurs.
 */
float
ParsetreeSampleFromFInside_b_jd_me(ESL_RANDOMNESS *r, CM_t *cm, ESL_DSQ *dsq, int L, float ***alpha, CP9Bands_t *cp9b, Parsetree_t **ret_tr, float ****ret_alpha)
{
  int          status;             /* easel status code */
  int          v, y, z, b;         /* state indices */
  int          yoffset;            /* transition offset in a states transition vector */
  int          i, j;               /* sequence position indices */
  int          jp_v, jp_y, jp_z;   /* positions, offset inside j band */
  int          kmin, kmax;         /* min/max k in current d band */
  int          d;                  /* j - i + 1; the current subseq length */
  int          dp_v, dp_y;         /* length, offset inside a d band */
  int          k;                  /* right subseq fragment length for bifurcs */
  int          kp_z;               /* right fragment length, offset inside a d band */
  int          nd;                 /* node index */
  int          bifparent;          /* for connecting bifurcs */
  Parsetree_t *tr;                 /* trace we're building */
  ESL_STACK   *pda;                /* the stack */
  float        pvec[MAXCONNECT+1]; /* prob vector of possible paths to take, (max num children + 1 for possibility of EL) */
  float       *bifvec;             /* pvec for choosing transition out of BIF_B states */
  float       *rootvec;            /* pvec for choosing transition out of ROOT_S if local begins are on */
  float        maxsc;              /* max score in our vector of scores of possible subparses */
  int          el_is_possible;     /* TRUE if we can jump to EL from current state (and we're in local mode) FALSE if not */
  int          ntrans;             /* number of transitions for current state */
  int          fsc = 0.;           /* score of the parsetree we're sampling */
  int          seen_valid;         /* for checking we have at least one valid path to take  */
  int          sd;                 /* state delta for current state, residues emitted left + residues emitted right */
  int          sdr;                /* state right delta for current state, residues emitted right */
  int         *jmin;               /* ptr to CP9Bands cp9b's j min band, for convenience */
  int         *jmax;               /* ptr to CP9Bands cp9b's j max band, for convenience */
  int        **hdmin;              /* ptr to CP9Bands cp9b's d min band, for convenience */
  int        **hdmax;              /* ptr to CP9Bands cp9b's d max band, for convenience */

  /* for convenience */
  jmin  = cp9b->jmin;
  jmax  = cp9b->jmax;
  hdmin = cp9b->hdmin;
  hdmax = cp9b->hdmax;

  /* contract check */
  if(ret_tr == NULL) cm_Fail("ParsetreeSampleFromFInside_b_jd_me(), ret_tr is NULL.");
  if(r      == NULL) cm_Fail("ParsetreeSampleFromFInside_b_jd_me(), source of randomness r is NULL.");
  
  /* initialize pvec */
  esl_vec_FSet(pvec, (MAXCONNECT+1), 0.);

  /* Create a parse tree structure and initialize it by adding the root state.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Stochastically traceback through the Inside matrix 
   * this section of code is stolen and adapted from hbandcyk.c:insideT_b_jd_me() 
   */
  pda = esl_stack_ICreate();
  v = 0;

  j = d = L;
  i = 1;

  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];
  fsc  = 0.;
  while (1) {
    if(cm->sttype[v] != EL_st && d > hdmax[v][jp_v])
      cm_Fail("ERROR in ParsetreeSampleFromFInside_b_jd_me(). d : %d > hdmax[%d] (%d)\n", d, v, hdmax[v][jp_v]);
    if(cm->sttype[v] != EL_st && d < hdmin[v][jp_v])
      cm_Fail("ERROR in ParsetreeSampleFromFInside_b_jd_me(). d : %d < hdmin[%d] (%d)\n", d, v, hdmin[v][jp_v]);

    if (cm->sttype[v] == B_st) {
      y = cm->cfirst[v];
      z = cm->cnum[v];
      jp_z = j-jmin[z];
      k = kp_z + hdmin[z][jp_z];  /* k = offset len of right fragment */

      ESL_ALLOC(bifvec, sizeof(float) * (d+1));
      /* set bifvec[] as (float-ized) log odds scores for each valid left fragment length,
       * we have to be careful to check that the corresponding alpha cell for each length is valid
       */
      esl_vec_FSet(bifvec, (d+1), IMPOSSIBLE); /* only valid d's will be reset to a non-IMPOSSIBLE score */

      /* This search for valid k's is complex, and uncommented. It was taken from
       * hbandcyk.c:inside_b_jd_me(), the B_st case. The code there is commented somewhat
       * extensively. I'm pretty sure this is the most efficient (or at least close to it) 
       * way to find the valid cells in the DP matrix we're looking for. 
       */
      jp_v = j - jmin[v];
      jp_y = j - jmin[y];
      jp_z = j - jmin[z];
      if(j < jmin[v] || j > jmax[v]) cm_Fail("ParsetreeSampleFromFInside_b_jd_me() B_st v: %d j: %d outside band jmin: %d jmax: %d\n", v, j, jmin[v], jmax[v]);
      if(d < hdmin[v][jp_v] || d > hdmax[v][jp_v]) cm_Fail("ParsetreeSampleFromFInside_b_jd_me() B_st v: %d j: %d d: %d outside band dmin: %d dmax: %d\n", v, j, d, hdmin[v][jp_v], hdmax[v][jp_v]);
      seen_valid = FALSE;
      kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
      kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
      for(k = kmin; k <= kmax; k++)
	{
	  if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
	    {
	      kp_z = k-hdmin[z][jp_z];
	      dp_y = d-hdmin[y][jp_y-k];
	      bifvec[k] = alpha[y][jp_y-k][dp_y-k] + alpha[z][jp_z][kp_z]; 
	      seen_valid = TRUE;
	    }
	}
      if(!seen_valid) cm_Fail("ParsetreeSampleFromFInside_b_jd_me() number of valid transitions (for a B_st) is 0. You thought this was impossible.");
      maxsc = esl_vec_FMax(bifvec, (d+1));
      esl_vec_FIncrement(bifvec, (d+1), (-1. * maxsc));
      esl_vec_FScale(bifvec, (d+1), log(2));
      esl_vec_FLogNorm(bifvec, (d+1));
      k = esl_rnd_FChoose(r, bifvec, (d+1));
      free(bifvec);

      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */

      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);

      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else {
      if((v > 0) || (! (cm->flags & CMH_LOCAL_BEGIN))) /* ROOT_S with local begins on is a special case that we handle below */
	{ 
	  /* Choose which transition we take.
	   * Set pvec[] as (float-ized) log odds scores for each child we can transit to, 
	   * plus a local end (if possible). We only want to look at valid transitions, that
	   * is, those that do not violate the bands (correspond to accessing cells that actually
	   * exist in the DP matrix). 
	   */
	  seen_valid = FALSE;
	  esl_vec_FSet(pvec, (MAXCONNECT+1), IMPOSSIBLE); /* only transitions that correspond to valid cells will be reset to a non-IMPOSSIBLE score */
	  fsc += get_femission_score(cm, dsq, v, i, j); 
	  sdr = StateRightDelta(cm->sttype[v]);
	  sd  = StateDelta(cm->sttype[v]);
	  for(yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
	    {
	      y = yoffset + cm->cfirst[v];
	      if((j - sdr) >= jmin[y] && (j - sdr) <= jmax[y]) 
		{ /* enforces j is valid for state y */
		  jp_y = j - jmin[y];
		  if((d - sd) >= hdmin[y][jp_y-sdr] && (d - sd) <= hdmax[y][jp_y-sdr])
		    {
		      dp_y = d - hdmin[y][(jp_y - sdr)];  /* d index for state y 
							     in alpha w/mem eff bands */
		      /* if we get here alpha[y][jp_y-sdr][dp_y-sd] is a valid alpha cell
		       * corresponding to alpha[y][j-sdr][d-sd] in the platonic matrix.
		       */
		      pvec[yoffset] = cm->tsc[v][yoffset] + alpha[y][jp_y - sdr][dp_y - sd];
		      seen_valid = TRUE;
		    }
		}		
	    }
	  if(!seen_valid) {
	    cm_Fail("ParsetreeSampleFromFInside_b_jd_me() number of valid transitions is 0. You thought this was impossible.");
	  }
	  if((cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v])) 
	    el_is_possible = TRUE; 
	  else 
	    el_is_possible = FALSE;
	  if(el_is_possible) pvec[cm->cnum[v]] = cm->endsc[v] + alpha[cm->M][j][d]; /* EL is silent when we transition into it from non-EL */
	  ntrans = cm->cnum[v] + el_is_possible;
	  maxsc = esl_vec_FMax(pvec, ntrans);
	  esl_vec_FIncrement(pvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (pvec, ntrans, log(2));
	  esl_vec_FLogNorm(pvec, ntrans);
	  yoffset = esl_rnd_FChoose(r, pvec, ntrans);
	  if(yoffset < cm->cnum[v]) fsc += cm->tsc[v][yoffset]; 
	  else {
	    fsc += cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
	    yoffset = USED_EL; /* we chose EL */
	  }
	}
      else /* v == 0 && (cm->flags && CMH_LOCAL_BEGIN) ( local begins are on )*/
	{
	  seen_valid = FALSE;
	  ntrans = cm->M; /* pretend all states are possible to begin into, but they're not as some will remain IMPOSSIBLE */
	  ESL_ALLOC(rootvec, sizeof(float) * (ntrans));
	  esl_vec_FSet(rootvec, ntrans, IMPOSSIBLE);

	  /* Set all the legal states that we can local begin into to appropriate scores.
	   * Only states y that have a non-zero cm->beginsc[y] AND have alpha[y][j][d]
	   * within their bands are legal.
	   */
	  for (nd = 1; nd < cm->nodes; nd++) {
	    if ((nd == 1) || /* we can transit into node 1 no matter what */
		(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
		 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd))
	      {
		y = cm->nodemap[nd];
		if(j >= jmin[y] && j <= jmax[y]) 
		  { /* enforces j is valid for state y */
		    jp_y = j - jmin[y];
		    if(d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y])
		      {
			dp_y = d - hdmin[y][jp_y];
			rootvec[y] = cm->beginsc[y] + alpha[y][jp_y][dp_y]; /* ROOT_S is silent */
			seen_valid = TRUE;
		      }
		  }
	      }
	  }
	  if(!seen_valid) cm_Fail("ParsetreeSampleFromFInside_b_jd_me() number of valid transitions (from ROOT_S!) is 0. You thought this was impossible.");
	  /* this block is shared with v > 0 block, but we repeat it here so we don't need another if statement */
	  maxsc = esl_vec_FMax(rootvec, ntrans);
	  esl_vec_FIncrement(rootvec, ntrans, (-1. * maxsc));
	  /* get from log_2 to log_e, so we can use easel's log vec ops */
	  esl_vec_FScale  (rootvec, ntrans, log(2));
	  esl_vec_FLogNorm(rootvec, ntrans);
	  b = esl_rnd_FChoose(r, rootvec, ntrans);
	  /* end of similar block with v > 0 */
	  fsc += cm->beginsc[b];
	  yoffset = USED_LOCAL_BEGIN; 
	  free(rootvec); /* we will not need this again */
	}

      /*printf("v : %d | r : %d | z : %d | 1 : %d | \n", v, r, z, 1);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	  jp_v = j;
	  dp_v = d;
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */

  /* If the caller doesn't want the alpha matrix, free it
   */
  if (ret_alpha == NULL) free_vjd_matrix(alpha, cm->M, 1, L);
  else *ret_alpha = alpha;

  *ret_tr = tr; /* contract checked ret_tr was non-NULL */
  return fsc;

 ERROR:
  cm_Fail("memory error.");
  return 0.; /* NEVERREACHED */
}



/*
 * Function: get_femission_score()
 * Incept:   EPN, Thu Nov 15 16:48:56 2007
 *          
 * Purpose:  Given a CM, dsq, state index and coordinates return the float emission
 *           score.
 *           
 * Args:     cm       - the model
 *           dsq      - digitized sequence
 *           v        - state index
 *           i        - dsq index for first position of subseq for subtree at v
 *           j        - dsq index for last position of subseq for subtree at v
 *
 * Return:   float emission score, 0 if state is non-emitter.
 */
float
get_femission_score(CM_t *cm, ESL_DSQ *dsq, int v, int i, int j)
{
  if     (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) return cm->oesc[v][dsq[i]];
  else if(cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) return cm->oesc[v][dsq[j]];
  else if(cm->sttype[v] == MP_st)                           return cm->oesc[v][dsq[i]*cm->abc->Kp+dsq[j]];
  else return 0.;
}



/*******************************************************************************
 * 11.04.05
 * EPN 
 * Memory efficient banded versions of selected smallcyk.c functions that
 * enforce bands in the d and j dimensions informed by an HMM Forward/Backward
 * posterior decode of the target sequence.
 * 
 * These functions are modified from their originals in smallcyk.c to make 
 * HMM banded FULL (not D&C) CYK alignment memory efficient. The starting
 * point for CYKInside_b_jd() was CYKInside_b_me() in smallcyk.c. The main
 * difference is that bands in the j dimension are enforced, and the d
 * bands have j dependence. 
 * 
 * CYK_Inside_b_jd() only allocates cells within the j AND d bands.
 * 
 * Comments from smallcyk.c pertaining to CYK_Inside_b_me():
 * .........................................................................
 * The only real difficulty implementing memory efficient
 * bands is in being able to determine what cell alpha[v][j][d] from the 
 * non-memory efficient code corresponds to in the memory-efficient code (we'll
 * call the corresponding cell a[v'][j'][d'] or a[vp][jp][dp]).  The reason
 * v != v'; j != j' and d != d' is because the primes are offset due to the
 * fact that some of the original alpha matrix deck (a[v]) has not been allocated
 * due to the bands.  Therefore all of the differences between the *_b_me() functions
 * and their *_b() versions is to deal with the offset issue.
 * .........................................................................
 * 
 *******************************************************************************/

/* The alignment engine. 
 */
static float inside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, 
			    int r, int z, int i0, int j0, 
			    int do_full,
			    float ***alpha, float ****ret_alpha, 
			    struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			    void ****ret_shadow, 
			    int allow_begin, int *ret_b, float *ret_bsc,
			    int *jmin, int *jmax,
			    int **hdmin, int **hdmax,
			    int *safe_hdmin, int *safe_hdmax);

/* The traceback routine.
 */

static float insideT_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			     int r, int z, int i0, int j0, int allow_begin,
			     int *jmin, int *jax,
			     int **hdmin, int **hdmax,
			     int *safe_hdmin, int *safe_hdmax);

#define BE_EFFICIENT  0		/* setting for do_full: small memory mode */
#define BE_PARANOID   1		/* setting for do_full: keep whole matrix, perhaps for debugging */

/* Special flags for use in shadow (traceback) matrices, instead of
 * offsets to connected states. When yshad[0][][] is USED_LOCAL_BEGIN,
 * the b value returned by inside() is the best connected state (a 0->b
 * local entry). When yshad[v][][] is USED_EL, there is a v->EL transition
 * and the remaining subsequence is aligned to the EL state. 
 */
#define USED_LOCAL_BEGIN 101
#define USED_EL          102


/* Function: CYKInside_b_jd()
 *           EPN 11.04.05
 * based on CYKInside_b() which was based on CYKInside()
 *
 * Only difference is bands are used in d and j dimesions: 
 *
 * Date:     SRE, Sun Jun  3 19:48:33 2001 [St. Louis]
 *
 * Purpose:  Wrapper for the insideT_b_jd_me() routine - solve
 *           a full alignment problem, return the traceback
 *           and the score, without dividing & conquering, using bands.
 *           
 *           Analogous to CYKDivideAndConquer() in many respects;
 *           see the more extensive comments in that function for
 *           more details on shared aspects.
 *           
 * Args:     cm     - the covariance model
 *           sq     - the sequence, 1..L
 *           r      - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *           dmin   - minimum d bound for each state v; [0..v..M-1]
 *           dmax   - maximum d bound for each state v; [0..v..M-1]
 *
 * Returns:  score of the alignment in bits.
 */
float
CYKInside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, 
	       int *jmin, int *jmax, int **hdmin, int **hdmax, int *dmin, int *dmax)
{
  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in CYKInside_b_jd(), dsq is NULL.\n");

  Parsetree_t *tr;
  int          z;
  float        sc;

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));
    printf("alignment strategy:CYKInside_b_jd:b:nosmall\n"); 
    printf("L: %d\n", L);*/

  /* Trust, but verify.
   * Check out input parameters.
   */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CMH_LOCAL_BEGIN)) cm_Fail("internal error: we're not in local mode, but r is not root");
    if (cm->stid[r] != MATP_MP && cm->stid[r] != MATL_ML &&
	cm->stid[r] != MATR_MR && cm->stid[r] != BIF_B)
      cm_Fail("internal error: trying to do a local begin at a non-mainline start");
  }

  /* Create the parse tree, and initialize.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */
  z  = cm->M-1;
  sc = 0.;

  /* Deal with case where we already know a local entry transition 0->r
   */
  if (r != 0)
    {
      InsertTraceNode(tr, 0,  TRACE_LEFT_CHILD, i0, j0, r);
      z  =  CMSubtreeFindEnd(cm, r);
      sc =  cm->beginsc[r];
    }

  /* Solve the whole thing with one call to insideT_b_jd.  This calls
     a memory efficient insideT function, which only allocates cells
     in alpha within the bands. 
   */
  sc += insideT_b_jd_me(cm, dsq, L, tr, r, z, i0, j0, (r==0), jmin, jmax, hdmin, hdmax, 
			dmin, dmax);

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  ESL_DPRINTF1(("returning from CYKInside_b_jd() sc : %f\n", sc));

  return sc;
}


 
/* EPN 03.29.06
 * Function: inside_b_jd_me()
 * based on inside_b_me() which was ...
 * based on inside()
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phase of a CYK alignment algorithm
 *           using bands in the j and d dimension from obtained
 *           from an HMM forwards-backwards run. This function
 *           is memory efficient in the j AND d dimension.
 * 
 *           To be able to consistently handle end states, the
 *           original SRE behavior of reusing the end deck was
 *           abandoned. Now each end state has its own deck, which
 *           makes this implementation easier because each state
 *           has its own bands on j, and thus has a state specific
 *           offset with alpha[end][jp][dp] in the banded mem eff 
 *           matrix corresponding to alpha[end][jp+jmin[end]][dp+hdmin[v][jp_v]]
 *           in the platonic matrix.
 *           
 *           The deck re-use strategy in general does not work with
 *           this implementation b/c each state has it's own j-specific
 *           bands. 
 *
 *           Notes from inside():
 *           A note on the loop conventions. We're going to keep the
 *           sequence (dsq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: only vroot..vend
 *           in the decks, i0-1..j in the rows, and up to j0-i0+1 in
 *           the columns (d dimension). Where this is handled the most
 *           is in two variables: W, which is the length of the subsequence
 *           (j0-i0+1), and is oft used in place of L in the usual CYK;
 *           and jp (read: j'), which is the *relative* j w.r.t. the
 *           subsequence, ranging from 0..W, and then d ranges from 
 *           0 to jp, and j is calculated from jp (i0-1+jp).
 *           In this banded version, there are more offset issues,
 *           these are detailed with comments in the code.
 *
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in an
 *           alpha, we expect that alpha[vroot..vend] are all NULL
 *           decks already; any other decks <vroot and >vend will
 *           be preserved. If we pass in a dpool, the decks *must* be
 *           sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *           We also deal with local begins, by keeping track of the optimal
 *           state that we could enter and account for the whole target 
 *           sequence: b = argmax_v  alpha_v(i0,j0) + log t_0(v),
 *           and bsc is the score for that. 
 *
 *           If vroot==0, i0==1, and j0==L (e.g. a complete alignment),
 *           the optimal alignment might use a local begin transition, 0->b,
 *           and we'd have to be able to trace that back. For any
 *           problem where the caller sets allow_begin, we return a valid b 
 *           (the optimal 0->b choice) and bsc (the score if 0->b is used).
 *           If a local begin is part of the optimal parse tree, the optimal
 *           alignment score returned by inside() will be bsc and yshad[0][L][L] 
 *           will be USE_LOCAL_BEGIN, telling insideT() to check b and
 *           start with a local 0->b entry transition. When inside()
 *           is called on smaller subproblems (v != 0 || i0 > 1 || j0
 *           < L), we're using inside() as an engine in divide &
 *           conquer, and we don't use the overall return score nor
 *           shadow matrices, but we do need allow_begin, b, and bsc for
 *           divide&conquer to sort out where a local begin might be used.
 *
 * Args:     cm        - the model    [0..M-1]
 *           sq        - the sequence [1..L]   
 *                     - length of the dsq
 *           vroot     - first start state of subtree (0, for whole model)
 *           vend      - last end state of subtree (cm->M-1, for whole model)
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           alpha     - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_alpha - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * int *safe_hdmin     - safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0)
 * int *safe_hdmax     - safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0)
 *                       
 * Returns: Score of the optimal alignment.  
 */
static float 
inside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	       float ***alpha, float ****ret_alpha, 
	       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	       void ****ret_shadow, 
	       int allow_begin, int *ret_b, float *ret_bsc,
	       int *jmin, int *jmax, int **hdmin, int **hdmax,
	       int *safe_hdmin, int *safe_hdmax)
{
  int      status;
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  void  ***shadow;      /* shadow matrix for tracebacks */
  int    **kshad;       /* a shadow deck for bifurcations */
  char   **yshad;       /* a shadow deck for every other kind of state */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */

  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;
  float  **tmp_deck;       /* temp variable, used only to free deckpool at end */

  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in inside_b_jd_me(), dsq is NULL.\n");

  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;     v < vroot; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) touch[v] = cm->pnum[v];
  for (v = vend+1;v < cm->M; v++) touch[v] = 0;

  /* The shadow matrix, if caller wants a traceback.
   * We do some pointer tricks here to save memory. The shadow matrix
   * is a void ***. Decks may either be char ** (usually) or
   * int ** (for bifurcation decks). Watch out for the casts.
   * For most states we only need
   * to keep y as traceback info, and y <= 6. For bifurcations,
   * we need to keep k, and k <= L, and L might be fairly big.
   * (We could probably limit k to an unsigned short ... anyone
   * aligning an RNA > 65536 would need a big computer... but
   * we'll hold off on that for now. We could also pack more
   * traceback pointers into a smaller space since we only really
   * need 3 bits, not 8.)
   */
  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* Main recursion
   */
  for (v = vend; v >= vroot; v--) 
    {
      /* First we need a deck to fill in. With memory efficient bands 
       * we don't reuse decks b/c each state has different bands and therefore
       * different deck sizes, so we ALWAYS allocate a deck here.
       */
      alpha[v] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
      //printf("allocated 
      if (cm->sttype[v] != E_st) {
	if (ret_shadow != NULL) {
	  if (cm->sttype[v] == B_st) {
	    kshad     = alloc_jdbanded_vjd_kshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	    shadow[v] = (void **) kshad;
	  } else {
	    yshad     = alloc_jdbanded_vjd_yshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	    shadow[v] = (void **) yshad;
	  }
	}
      }

      /* We've only allocated alpha cells that are within the bands
       * on the j and d dimensions. This means we have to deal
       * with all sorts of offset issues, but we don't have to 
       * waste time setting cells outside the bands to IMPOSSIBLE.
       */
      if (cm->sttype[v] == E_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  if(d != 0)
		    cm_Fail("band on E state %d has a non-zero d value within its j band for j:%d\n", v, j);
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v
						 in alpha w/mem eff bands */
		  alpha[v][jp_v][dp_v] = 0.; /* for End states, d must be 0 */
		}		    
	    }
	  continue;
	}  
      else if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		  alpha[v][jp_v][dp_v]  = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  /* treat EL as emitting only on self transition */
		  if (ret_shadow != NULL) yshad[jp_v][dp_v]  = USED_EL; 
		  for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		    {
		      yoffset = y - cm->cfirst[v];
		      if(j >= jmin[y] && j <= jmax[y]) 
			/* Enforces j is valid for state y */
			{
			  jp_y = j - jmin[y];
			  if(d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y])
			    {
			      dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							     in alpha w/mem eff bands */
			      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			       * corresponding to alpha[y][j][d] in the platonic matrix.
			       */
			      if ((sc = alpha[y][jp_y][dp_y] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
				{
				  alpha[v][jp_v][dp_v] = sc; 
				  if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
				}
			    }
			}
		    }
		  if (alpha[v][jp_v][dp_v] < IMPOSSIBLE)
		    alpha[v][jp_v][dp_v] = IMPOSSIBLE;
		}
	    }
	}
      else if (cm->sttype[v] == B_st)
	{
	  y = cm->cfirst[v];
	  z = cm->cnum[v];
	  /* Any valid j must be within both state v and state z's j band 
	   * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
	   * constructed, but we'll check anyway. 
	   */
	  tmp_jmin = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	  tmp_jmax = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];

	  /* For any values of j within v's j band but outside of z's j band,
	   * we have to set the corresponding alpha cells to IMPOSSIBLE.
	   * This is done be the following two ugly for loops, 
	   * which will only be looked at once for each B state, and
	   * even then only *very* rarely entered. This
	   * is why they're here, seemingly out of place before the 
	   * main j loop below, where similar performing code would be 
	   * looked at on the order of j times, instead of just once.
	   */
	  for(j = jmin[v]; j < tmp_jmin; j++)
	    {
	      jp_v = j-jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d-hdmin[v][jp_v];
		  alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
		}
	    }
	  if(tmp_jmax < jmax[v])
	    for(j = (tmp_jmax+1); j <= jmax[v]; j++)
	      {
		jp_v = j-jmin[v];
		for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		  {
		    dp_v = d-hdmin[v][jp_v];
		    alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
		  }
	      }
	  /* the main j loop */
	  for (j = tmp_jmin; j <= tmp_jmax; j++)
	    {
	      jp_v = j - jmin[v];
	      jp_y = j - jmin[y];
	      jp_z = j - jmin[z];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */

		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 6 inequalities (some may be redundant):
		   * (1) k >= j-jmax[y];
		   * (2) k <= j-jmin[y]; 
		   *     1 and 2 guarantee (j-k) is within state y's j band
		   *
		   * (3) k >= hdmin[z][j-jmin[z]];
		   * (4) k <= hdmax[z][j-jmin[z]]; 
		   *     3 and 4 guarantee k is within z's j=(j), d band
		   *
		   * (5) k >= d-hdmax[y][j-jmin[y]-k];
		   * (6) k <= d-hdmin[y][j-jmin[y]-k]; 
		   *     5 and 6 guarantee (d-k) is within state y's j=(j-k) d band
		   */
		  kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
		  /* kmin satisfies inequalities (1) and (3) */
		  kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
		  /* kmax satisfies inequalities (2) and (4) */
		  /* RHS of inequalities 5 and 6 are dependent on k, so we check
		   * for these within the next for loop.
		   */
		  alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* initialize */
		  for(k = kmin; k <= kmax; k++)
		    {
		      if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
			{
			  /* for current k, all 6 inequalities have been satisified 
			   * so we know the cells corresponding to the platonic 
			   * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
			   * alpha[z][j][k] are all within the bands. These
			   * cells correspond to alpha[v][jp_v][dp_v], 
			   * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
			   * and alpha[z][jp_z][k-hdmin[jp_z]];
			   */
			  kp_z = k-hdmin[z][jp_z];
			  dp_y = d-hdmin[y][jp_y-k];

			  if ((sc = alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]) 
			      > alpha[v][jp_v][dp_v])
			    {
			      alpha[v][jp_v][dp_v] = sc;
			      if (ret_shadow != NULL) kshad[jp_v][dp_v] = kp_z;
			    }
			}
		    }
		  if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
		}
	    }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if(ret_shadow != NULL) yshad[jp_v][dp_v] = USED_EL;
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces (j-1) is valid for state y */

		      {
			jp_y = j - jmin[y];
			if((d-2) >= hdmin[y][jp_y-1] && (d-2) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-2] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-2] in the platonic matrix.
			     */
			    if ((sc = alpha[y][jp_y-1][dp_y-2] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
			      {
				alpha[v][jp_v][dp_v] = sc; 
				if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
			      }
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][jp_v][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if(ret_shadow != NULL) yshad[jp_v][dp_v] = USED_EL;
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if(j >= jmin[y] && j <= jmax[y]) /* Enforces j is valid for state y */
		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y] && (d-1) <= hdmax[y][jp_y])
			  {
			    dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							   in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j][d-1] in the platonic matrix.
			     */
			    if ((sc = alpha[y][jp_y][dp_y-1] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
			      {
				alpha[v][jp_v][dp_v] = sc; 
				if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
			      }
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if(ret_shadow != NULL) yshad[jp_v][dp_v] = USED_EL;
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces j-1 is valid for state y */

		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y-1] && (d-1) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-1] in the platonic matrix.
			     */
			    if ((sc = alpha[y][jp_y-1][dp_y-1] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
			      {
				alpha[v][jp_v][dp_v] = sc; 
				if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
			      }
			  }
		      }
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      /*if((cm->sttype[v] != IL_st) && (cm->sttype[v] != IR_st) && (cm->sttype[v] != B_st)) {
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++, i--) {
	    printf("alpha[v: %4d][jp_v: %4d][dp_v: %4d]: %.4f\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v]);
	    
	  }
	  printf("\n");
	}
	printf("\n\n");
	}*/
  
      /* The following loops originally access alpha[v][j0][W] but the index W will be
	 in different positions due to the bands */
      if(j0 >= jmin[v] && j0 <= jmax[v])
	{
	  jp_v = j0 - jmin[v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v])
	    {
	      Wp = W - hdmin[v][jp_v];
	      /* If we get here alpha[v][jp_v][Wp] is a valid cell
	       * in the banded alpha matrix, corresponding to 
	       * alpha[v][j0][W] in the platonic matrix.
	       */
	      /* Check for local begin getting us to the root.
	       * This is "off-shadow": if/when we trace back, we'll handle this
	       * case separately (and we'll know to do it because we'll immediately
	       * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	       * to jump right to state b; see below)
	       */
	      if (allow_begin && alpha[v][jp_v][Wp] + cm->beginsc[v] > bsc) 
		{
		  b   = v;
		  bsc = alpha[v][jp_v][Wp] + cm->beginsc[v];
		}
	    }
	}
      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (v == 0)
	{
	  if(j0 >= jmin[0] && j0 <= jmax[0])
	    {
	      jp_v = j0 - jmin[v];
	      if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v])
		{
		  if (allow_begin && v == 0 && bsc > alpha[0][jp_v][Wp]) {
		    alpha[0][jp_v][Wp] = bsc;
		    if (ret_shadow != NULL) yshad[jp_v][Wp] = USED_LOCAL_BEGIN;
		  }
		}
	    }
	}
      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, free them, we don't 
       * reuse decks with bands b/c each state has different deck size.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { 
	    /* we can definitely release the S children of a bifurc. */
	    y = cm->cfirst[v];
	    z = cm->cnum[v];  
	    free_vjd_deck(alpha[y], i0, j0);
	    alpha[y] = NULL;
	    free_vjd_deck(alpha[z], i0, j0);
	    alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) 
		  {
		    free_vjd_deck(alpha[y], i0, j0);
		    alpha[y] = NULL;
		  }
	      }
	  }
      }
    } /* end loop over all v */
  /*debug_print_alpha_banded_jd(alpha, cm, L, jmin, jmax, hdmin, hdmax);*/

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid 
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  
  Wp = W - hdmin[vroot][j0-jmin[vroot]];
  sc =     alpha[vroot][j0-jmin[vroot]][Wp];

  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = vroot; v <= vend; v++) 
      if (alpha[v] != NULL) { 
	deckpool_push(dpool, alpha[v]); alpha[v] = NULL;
      }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (deckpool_pop(dpool, &tmp_deck)) free_vjd_deck(tmp_deck, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  /*printf("inside jd me returning sc: %f\n", sc);*/

  return sc;

 ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

/* Function: insideT_b_jd_me()
 *           EPN 03.29.06
 * *based on insideT(), only difference is memory efficient bands on the j and d dimensions
 *  are used : 
 *
 * Date:     SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]
 *
 * Purpose:  Call inside, get vjd shadow matrix;
 *           then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 */
static float
insideT_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
		int r, int z, int i0, int j0, 
		int allow_begin, int *jmin, int *jmax,
		int **hdmin, int **hdmax,
		int *safe_hdmin, int *safe_hdmax)
{
  /* Contract check */
  if(dsq == NULL)
    cm_Fail("ERROR in insideT_b_jd_me(), dsq is NULL.");

  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;
  int       jp_v;               /* j-jmin[v] for current j, and current v */
  int       dp_v;               /* d-hdmin[v][jp_v] for current j, current v, current d*/
  int       jp_z;               /* j-jmin[z] for current j, and current z */
  int       kp_z;               /* the k value (d dim) from the shadow matrix
				 * giving the len of right fragment offset in deck z,
				 * k = kp_z + hdmin[z][jp_z]*/

  sc = inside_b_jd_me(cm, dsq, L, r, z, i0, j0, 
		      BE_EFFICIENT,	/* memory-saving mode */
		      /*BE_PARANOID,*/	/* non-memory-saving mode */
		      NULL, NULL,	/* manage your own matrix, I don't want it */
		      NULL, NULL,	/* manage your own deckpool, I don't want it */
		      &shadow,		/* return a shadow matrix to me. */
		      allow_begin,      /* TRUE to allow local begins */
		      &b, &bsc,	/* if allow_begin is TRUE, gives info on optimal b */
		      jmin, jmax,    /* bands on j */
		      hdmin, hdmax,  /* j dependent bands on d */
		      safe_hdmin, safe_hdmax);

  pda = esl_stack_ICreate();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];

  while (1) {
    if(cm->sttype[v] != EL_st && d > hdmax[v][jp_v])
      cm_Fail("ERROR in insideT_b_jd(). d : %d > hdmax[%d] (%d)\n", d, v, hdmax[v]);
    if(cm->sttype[v] != EL_st && d < hdmin[v][jp_v])
      cm_Fail("ERROR in insideT_b_jd(). d : %d < hdmin[%d] (%d)\n", d, v, hdmin[v]);
    
    if (cm->sttype[v] == B_st) {
      kp_z = ((int **) shadow[v])[jp_v][dp_v];   /* kp = offset len of right fragment */
      z = cm->cnum[v];
      jp_z = j-jmin[z];
      k = kp_z + hdmin[z][jp_z];  /* k = offset len of right fragment */
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */
      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else {
      yoffset = ((char **) shadow[v])[jp_v][dp_v];
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	  jp_v = j;
	  dp_v = d;
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  return sc;
}

  
/* Functions: *jdbanded_*_vjd_*
 * EPN 03.29.06 these functions were derived from their 
 *              *_vjd_* analogs from SRE's smallcyk.c
 * notes from smallcyk.c:
 * Date:     SRE, Sat Aug 12 16:27:37 2000 [Titusville]
 *
 * Purpose:  Allocation and freeing of 3D matrices and 2D decks
 *           in the vjd coord system. These can be called on
 *           subsequences i..j, not just the full sequence 1..L,
 *           so they need i,j... if you're doing the full sequence
 *           just pass 1,L.
 *           
 *           Also deal with shadow matrices and shadow decks in the
 *           vjd coordinate system. Note that bifurcation shadow decks
 *           need more dynamic range than other shadow decks, hence
 *           a separation into "kshadow" (BIFURC) and "yshadow" (other
 *           states) decks, and some casting shenanigans in
 *           a full ***shadow matrix.
 *           
 *           Values in yshad are offsets to the next connected state,
 *           or a flag for local alignment. Possible offsets range from
 *           0..5 (maximum of 6 connected states). The flags are
 *           USED_LOCAL_BEGIN (101) and USED_EL (102), defined at
 *           the top of this file. Only yshad[0][L][L] (e.g. root state 0,
 *           aligned to the whole sequence) may be set to USED_LOCAL_BEGIN.
 *           (Remember that the dynamic range of yshad, as a char, is 
 *           0..127, in ANSI C; we don't know if a machine will make it
 *           signed or unsigned.)
 */
float **
alloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  int     status;
  float **a;
  int     jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;
  /*printf("in alloc JD banded vjd deck, L : %d, i : %d, j : %d, jmin : %d, jmax : %d\n", L, i, j, jmin, jmax);*/

  if(j < jmin || i > jmax)
    cm_Fail("ERROR called alloc_jdbanded_vjd_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);

  ESL_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  ESL_ALLOC(a, sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0; jp <= L;     jp++) a[jp]     = NULL;

  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jfirst: %d | jlast: %d\n", jfirst, jlast);
      printf("jp: %d | max : %d\n", jp, (jlast)); 
      printf("hdmax[%d]: %d\n", (jp-jmin), hdmax[jp-jmin]);
      */
      ESL_DASSERT2(hdmax[jp-jmin] <= (jp+1))
      /* Based on my current understanding the above line should never be false, if it is means there's a valid d
       * in the hd band that is invalid because its > j. I think I check, or ensure, that this
       * doesn't happen when I'm constructing the d bands.
       */
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*a is offset only the first (jlast-jfirst+1) elements will be non-NULL*/
      ESL_ALLOC(a[jp-jfirst], sizeof(float) * bw);
      /*printf("\tallocated a[%d] | bw: %d\n", (jp-jfirst), bw);*/
    }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

#if 0
/******************************************************************/
/* The below functions were written during debugging, and print
   out either the shadow or alpha matrix.  They are kept
   here just in case they're needed again.  Note : the functions
   that print out the entire matrix are really only useful
   when the BE_PARANOID flag is set, meaning that decks are
   never freed until the end.
*/
/*================================================================*/

/* Debugging functions that print info to STDOUT */
static void debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
					int **hdmin, int **hdmax);
static void debug_print_shadow_banded_jd(void ***shadow, CM_t *cm, int L, int *jmin, int *jmax, 
					 int **hdmin, int **hdmax);
static void debug_print_shadow_banded_deck_jd(int v, void ***shadow, CM_t *cm, int L, int *jmin, int *jmax,
					      int **hdmin, int **hdmax);

/* EPN 03.29.06
   debug_print_alpha_banded_jd()
 * Function: debug_print_alpha_banded_jd
 *
 * Purpose:  Print alpha matrix banded in j and d dimensions
 */
void
debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
			    int **hdmin, int **hdmax)
{
  int v, j, d, dp, jp, max_v;

  printf("\nPrinting banded alpha matrix :\n");
  printf("************************************\n");
  max_v = cm->M-1;
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      max_v = cm->M;
    }
  for(v = 0; v < max_v; v++)
    {
      printf("====================================\n");
      for(j = jmin[v]; j <= jmax[v]; j++)
	{
	  printf("------------------------------------\n");
	  for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	    {
	      jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	      dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	      printf("alpha_jd[%2d][%2d][%2d] : %6.2f | j: %4d | d: %4d\n", v, jp, dp, alpha[v][jp][dp], j, d);
	    }
	}
    }
  printf("****************\n\n");
}

/* EPN 05.16.05
   debug_print_shadow_banded()
 * Function: debug_print_shadow_banded
 *
 * Purpose:  Print banded shadow matrix 
 */
static void
debug_print_shadow_banded_jd(void ***shadow, CM_t *cm, int L, int *jmin, int *jmax, 
			     int **hdmin, int **hdmax)
{
  int v, j, d, dp, jp;
  int yoffset;
  char yoffset_c;

  printf("\nPrinting banded shadow matrix :\n");
  printf("************************************\n");
  for(v = 0; v < cm->M; v++)
    {
      printf("====================================\n");
      for(j = jmin[v]; j <= jmax[v]; j++)
	{
	  printf("------------------------------------\n");
	  for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	    {
	      jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	      dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	      if(cm->sttype[v] == E_st)
		{
		  printf("END state\n");
		}
	      else
		{
		  if(cm->sttype[v] == B_st)
		    {
		      yoffset = ((int **) shadow[v])[jp][dp];
		      printf("INT  shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
		    }
		  else
		    {
		      yoffset_c = ((char **) shadow[v])[jp][dp];
		      printf("CHAR shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
		    }
		}
	    }
	}
    }
  printf("****************\n\n");
}

/* EPN 05.16.05
   debug_print_shadow_banded_deck_jd()
 * Function: debug_print_shadow_banded_deck_jd
 *
 * Purpose:  Print banded (in j and d dimensions) shadow matrix deck
 */

static void
debug_print_shadow_banded_deck_jd(int v, void ***shadow, CM_t *cm, int L, int *jmin, int *jmax,
				  int **hdmin, int **hdmax)
{
  int j, d, dp, jp;
  int yoffset;

  printf("\nPrinting banded shadow matrix deck for v : %d:\n", v);
  printf("====================================\n");
  for(j = jmin[v]; j <= jmax[v]; j++)
    {
      printf("------------------------------------\n");
      for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	{
	  jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	  dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	  if(cm->sttype[v] == E_st)
	    {
	      printf("END state\n");
	    }
	  else
	    {
	      yoffset = ((char **) shadow[v])[jp][dp];
	      printf("shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
	    }
	}
    }
}
#endif

#if 0
/* Here are the non-memory efficient functions, kept around for reference */
/* The alignment engine (not memory efficient) */
static float inside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, 
			 int r, int z, int i0, int j0, 
			 int do_full,
			 float ***alpha, float ****ret_alpha, 
			 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			 void ****ret_shadow, 
			 int allow_begin, int *ret_b, float *ret_bsc,
			 int *jmin, int *jmax,
			 int **hdmin, int **hdmax,
			 int *safe_hdmin, int *safe_hdmax);

/* The traceback routine (not memory efficient) */
static float insideT_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			  int r, int z, int i0, int j0, int allow_begin,
			  int *jmin, int *jax, 
			  int **hdmin, int **hdmax,
			  int *safe_hdmin, int *safe_hdmax);

#endif


#if 0
/******************************************************************/
/* The below functions were written during debugging, and print
   out either the shadow or alpha matrix.  They are kept
   here just in case they're needed again.  Note : the functions
   that print out the entire matrix are really only useful
   when the BE_PARANOID flag is set, meaning that decks are
   never freed until the end.
*/
/*================================================================*/

/* Debugging functions that print info to STDOUT */
static void debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
					int **hdmin, int **hdmax);
static void debug_print_shadow_banded_jd(void ***shadow, CM_t *cm, int L, int *jmin, int *jmax, 
					 int **hdmin, int **hdmax);
static void debug_print_shadow_banded_deck_jd(int v, void ***shadow, CM_t *cm, int L, int *jmin, int *jmax,
					      int **hdmin, int **hdmax);

/* EPN 03.29.06
   debug_print_alpha_banded_jd()
 * Function: debug_print_alpha_banded_jd
 *
 * Purpose:  Print alpha matrix banded in j and d dimensions
 */
void
debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
			    int **hdmin, int **hdmax)
{
  int v, j, d, dp, jp, max_v;

  printf("\nPrinting banded alpha matrix :\n");
  printf("************************************\n");
  max_v = cm->M-1;
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      max_v = cm->M;
    }
  for(v = 0; v < max_v; v++)
    {
      printf("====================================\n");
      for(j = jmin[v]; j <= jmax[v]; j++)
	{
	  printf("------------------------------------\n");
	  for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	    {
	      jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	      dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	      printf("alpha_jd[%2d][%2d][%2d] : %6.2f | j: %4d | d: %4d\n", v, jp, dp, alpha[v][jp][dp], j, d);
	    }
	}
    }
  printf("****************\n\n");
}

/* EPN 05.16.05
   debug_print_shadow_banded()
 * Function: debug_print_shadow_banded
 *
 * Purpose:  Print banded shadow matrix 
 */
static void
debug_print_shadow_banded_jd(void ***shadow, CM_t *cm, int L, int *jmin, int *jmax, 
			     int **hdmin, int **hdmax)
{
  int v, j, d, dp, jp;
  int yoffset;
  char yoffset_c;

  printf("\nPrinting banded shadow matrix :\n");
  printf("************************************\n");
  for(v = 0; v < cm->M; v++)
    {
      printf("====================================\n");
      for(j = jmin[v]; j <= jmax[v]; j++)
	{
	  printf("------------------------------------\n");
	  for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	    {
	      jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	      dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	      if(cm->sttype[v] == E_st)
		{
		  printf("END state\n");
		}
	      else
		{
		  if(cm->sttype[v] == B_st)
		    {
		      yoffset = ((int **) shadow[v])[jp][dp];
		      printf("INT  shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
		    }
		  else
		    {
		      yoffset_c = ((char **) shadow[v])[jp][dp];
		      printf("CHAR shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
		    }
		}
	    }
	}
    }
  printf("****************\n\n");
}

/* EPN 05.16.05
   debug_print_shadow_banded_deck_jd()
 * Function: debug_print_shadow_banded_deck_jd
 *
 * Purpose:  Print banded (in j and d dimensions) shadow matrix deck
 */

static void
debug_print_shadow_banded_deck_jd(int v, void ***shadow, CM_t *cm, int L, int *jmin, int *jmax,
				  int **hdmin, int **hdmax)
{
  int j, d, dp, jp;
  int yoffset;

  printf("\nPrinting banded shadow matrix deck for v : %d:\n", v);
  printf("====================================\n");
  for(j = jmin[v]; j <= jmax[v]; j++)
    {
      printf("------------------------------------\n");
      for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	{
	  jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	  dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	  if(cm->sttype[v] == E_st)
	    {
	      printf("END state\n");
	    }
	  else
	    {
	      yoffset = ((char **) shadow[v])[jp][dp];
	      printf("shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
	    }
	}
    }
}
#endif

#if 0
/* Here are the non-memory efficient functions, kept around for reference */
/* The alignment engine (not memory efficient) */
static float inside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, 
			 int r, int z, int i0, int j0, 
			 int do_full,
			 float ***alpha, float ****ret_alpha, 
			 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			 void ****ret_shadow, 
			 int allow_begin, int *ret_b, float *ret_bsc,
			 int *jmin, int *jmax,
			 int **hdmin, int **hdmax,
			 int *safe_hdmin, int *safe_hdmax);

/* The traceback routine (not memory efficient) */
static float insideT_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			  int r, int z, int i0, int j0, int allow_begin,
			  int *jmin, int *jax, 
			  int **hdmin, int **hdmax,
			  int *safe_hdmin, int *safe_hdmax);


/* EPN 
 * Function: inside_b_jd()
 * based on inside_b_me() which was ...
 * based on inside()
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phased of a CYK alignment algorithm
 *           using bands in the j and d dimension from an 
 *           HMM forwards backwards run. This function
 *           is memory efficient in the d dimension (only
 *           allocate within "safe" d bands). 
 *           Further assumes we're aligning a full sequence 
 *           (1..L) NOT a subsequence (i0..j0), and aligns it
 *           to the full model. This is very different from inside()
 *           which aligns a subsequence to a subtree of the model.
 *           
 *           Notes from inside():
 *           A note on the loop conventions. We're going to keep the
 *           sequence (dsq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: only vroot..vend
 *           in the decks, i0-1..j in the rows, and up to j0-i0+1 in
 *           the columns (d dimension). Where this is handled the most
 *           is in two variables: W, which is the length of the subsequence
 *           (j0-i0+1), and is oft used in place of L in the usual CYK;
 *           and jp (read: j'), which is the *relative* j w.r.t. the
 *           subsequence, ranging from 0..W, and then d ranges from 
 *           0 to jp, and j is calculated from jp (i0-1+jp).
 *           
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in an
 *           alpha, we expect that alpha[vroot..vend] are all NULL
 *           decks already; any other decks <vroot and >vend will
 *           be preserved. If we pass in a dpool, the decks *must* be
 *           sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *           We also deal with local begins, by keeping track of the optimal
 *           state that we could enter and account for the whole target 
 *           sequence: b = argmax_v  alpha_v(i0,j0) + log t_0(v),
 *           and bsc is the score for that. 
 *
 *           If vroot==0, i0==1, and j0==L (e.g. a complete alignment),
 *           the optimal alignment might use a local begin transition, 0->b,
 *           and we'd have to be able to trace that back. For any
 *           problem where the caller sets allow_begin, we return a valid b 
 *           (the optimal 0->b choice) and bsc (the score if 0->b is used).
 *           If a local begin is part of the optimal parse tree, the optimal
 *           alignment score returned by inside() will be bsc and yshad[0][L][L] 
 *           will be USE_LOCAL_BEGIN, telling insideT() to check b and
 *           start with a local 0->b entry transition. When inside()
 *           is called on smaller subproblems (v != 0 || i0 > 1 || j0
 *           < L), we're using inside() as an engine in divide &
 *           conquer, and we don't use the overall return score nor
 *           shadow matrices, but we do need allow_begin, b, and bsc for
 *           divide&conquer to sort out where a local begin might be used.
 *
 * Args:     cm        - the model    [0..M-1]
 *           sq        - the sequence [1..L]   
 *           vroot     - first start state of subtree (0, for whole model)
 *           vend      - last end state of subtree (cm->M-1, for whole model)
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           alpha     - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_alpha - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * int *safe_hdmin     - safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0)
 * int *safe_hdmax     - safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0)
 *                       
 * Returns: Score of the optimal alignment.  
 */
static float 
inside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	    float ***alpha, float ****ret_alpha, 
	    struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	    void ****ret_shadow, 
	    int allow_begin, int *ret_b, float *ret_bsc,
	    int *jmin, int *jmax, int **hdmin, int **hdmax,
	    int *safe_hdmin, int *safe_hdmax)
{
  float  **end;         /* we re-use the end deck. */
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */

  void  ***shadow;      /* shadow matrix for tracebacks */
  int    **kshad;       /* a shadow deck for bifurcations */
  char   **yshad;       /* a shadow deck for every other kind of state */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */

  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      dp_z;           /* d index for state z in alpha w/mem eff bands */
  int      kp;             /* k prime - keeps track of what k should be now
			     that we're using memory efficient bands */
  int      Wp;             /* W also changes depending on state */

  if(dsq == NULL)
    esl_fatal("ERROR, dsq is NULL.");

  /* 11.04.05 jd addition: */
  if(i0 != 1)
    {
      printf("inside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");
      exit(1);
    }
  if(j0 != L)
    {
      printf("inside_b_jd requires that j0 be L. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
  if(vroot != 0)
    {
      printf("inside_b_jd requires that vroot be 0. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
  if(vend != cm->M-1)
    {
      printf("inside_b_jd requires that vend be cm->M-1. This function is not set up for subsequence alignment.\n");
      exit(1);
    }

  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops 
				 * This must be L because i0 must be 1 and j0 must be L
				 */
  
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, vroot, E_st);
  for (j = 0; j <= W; j++) {
    end[j][0] = 0.;
    for (d = 1; d <= j; d++) end[j][d] = IMPOSSIBLE;
  }

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;     v < vroot; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) touch[v] = cm->pnum[v];
  for (v = vend+1;v < cm->M; v++) touch[v] = 0;

  /* The shadow matrix, if caller wants a traceback.
   * We do some pointer tricks here to save memory. The shadow matrix
   * is a void ***. Decks may either be char ** (usually) or
   * int ** (for bifurcation decks). Watch out for the casts.
   * For most states we only need
   * to keep y as traceback info, and y <= 6. For bifurcations,
   * we need to keep k, and k <= L, and L might be fairly big.
   * (We could probably limit k to an unsigned short ... anyone
   * aligning an RNA > 65536 would need a big computer... but
   * we'll hold off on that for now. We could also pack more
   * traceback pointers into a smaller space since we only really
   * need 3 bits, not 8.)
   */
  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* Main recursion
   */
  for (v = vend; v >= vroot; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      } 
      if (! deckpool_pop(dpool, &(alpha[v]))) 
	/* CYK Full ME Bands used 1 */
	/* original line : alpha[v] = alloc_vjd_deck(L, i0, j0);*/
	alpha[v] = alloc_banded_vjd_deck(L, i0, j0, safe_hdmin[v], safe_hdmax[v]);
      
      if (ret_shadow != NULL) {
	if (cm->sttype[v] == B_st) {
	  /* CYK Full ME Bands used 2 */
	  /* original line : kshad     = alloc_vjd_kshadow_deck(L, i0, j0); */
	  kshad     = alloc_banded_vjd_kshadow_deck(L, i0, j0, safe_hdmin[v], safe_hdmax[v]);
	  shadow[v] = (void **) kshad;
	} else {
	  /* CYK Full ME Bands used 3 */
	  /* original line : yshad     = alloc_vjd_yshadow_deck(L, i0, j0); */
	  yshad     = alloc_banded_vjd_yshadow_deck(L, i0, j0, safe_hdmin[v], safe_hdmax[v]);
	  shadow[v] = (void **) yshad;
	}
      }

      /* 11.05.05
       * One strategy is to set all cells OUTSIDE bands to IMPOSSIBLE.
       * I think I'll run into problems doing this because some cells
       * are inside the j bands and inside the safe_hd bands, but not
       * inside the j dependent d bands. These cells though allocated, 
       * will potentially never get filled. 
       * One way to deal with this (though inefficient) is to set
       * ALL cells to impossible. Below is the the first strategy, only 
       * setting some cells to impossible. 
       
       **************************************************************
       * Strategy 1: only set cells outside j bands to IMPOSSIBLE:
       **************************************************************
       */
      /* j bands used 1.
       * Set all cells for j's outside of bands to IMPOSSIBLE 
       * Further, set any cells outside of hd j specific band for valid j's 
       * to IMPOSSIBLE. Remember, we've allocated only 
       * (safe_hdmax[v] - safe_hdmin[v] +1) cells for each vj deck.
       * Take advantage of fact that we know we're aligning the full sequence 1..L.
       */
      /* Following loop starts at safe_hdmin[v] because the j section
       * of a banded vjd deck is not allocated if j < dmin[v] because
       * there's no way that j can be used.
       */
      /*
      for (j = safe_hdmin[v]; j < jmin[v]; j++)
      {
      */
	  /* this j is outside the j band, set all d to IMPOSSIBLE */
      /*
	  for (d = safe_hdmin[v]; d <= safe_hdmax[v] && d <= j; d++)
	    {
	      alpha[v][j][d-safe_hdmin[v]] = IMPOSSIBLE;
	    }
	}

      for (j = jmax[v] + 1; j <= W; j++)
	{
      */
	  /* this j is outside the j band, set all d to IMPOSSIBLE */
      /*	  for (d = safe_hdmin[v]; d <= safe_hdmax[v] && d <= j; d++)
	    {
	      alpha[v][j][d-safe_hdmin[v]] = IMPOSSIBLE;
	    }
	}
	*************************************************************
	*/

      /**************************************************************
	* Strategy 2: set all allocated cells to IMPOSSIBLE at first.
	**************************************************************
	*/

      for (j = safe_hdmin[v]; j <= W; j++)
	{
	  for (d = safe_hdmin[v]; d <= safe_hdmax[v] && d <= j; d++)
	    {
	      alpha[v][j][d-safe_hdmin[v]] = IMPOSSIBLE;
	    }
	}

      //printf("2 v: %d\n", v);
	/*************************************************************/
      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  /* j bands used 2. */

	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      //printf("3 j: %d\n", j);
	      for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
		{
		  assert(d >= safe_hdmin[v]);
		  assert(d <= safe_hdmax[v]);
		  
		  y = cm->cfirst[v];
		  /* CYK Full ME Bands used 4 begin block */
		  /* original block */
		  /* alpha[v][j][d]  = cm->endsc[v];*/	/* init w/ local end */ 
		  /*if (ret_shadow != NULL) yshad[j][d]  = USED_EL; 
		    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		    if ((sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc; 
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		    }
		    if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		  */
		  /* new ME block */
		  dp_v = d - safe_hdmin[v];  /* d index for state v in alpha w/mem eff bands */
		  alpha[v][j][dp_v]  = cm->endsc[v];	/* init w/ local end */
		  if (ret_shadow != NULL) yshad[j][dp_v]  = USED_EL; 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		    {
		      dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
							    in alpha w/mem eff bands */
		      /* check to make sure the cell we're about to query is within the
			 bands for state y; this might be more complex than necessary */
		      if((dp_y >= 0) && ((dp_y < (j - (safe_hdmin[y+yoffset]) + 1))
					 && (dp_y < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
			{
			  if ((sc = alpha[y+yoffset][j][dp_y] + cm->tsc[v][yoffset]) >  alpha[v][j][dp_v]) {
			    alpha[v][j][dp_v] = sc; 
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  }
			}
		    }
		  if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		  /* CYK Full ME Bands used 4 end block */
		}
	    }
	}
      else if (cm->sttype[v] == B_st)
	{
	  /* j bands used 3. */
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //printf("3 j: %d\n", j);
	    /* Bands used */
	    /* old line :	for (d = 0; d <= jp; d++) */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);
		y = cm->cfirst[v];
		z = cm->cnum[v];

		/* CYK Full ME Bands used 5 begin block */
		/* original block */
		/*
		alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0];
		if (ret_shadow != NULL) kshad[j][d] = 0;
		for (k = 1; k <= d; k++)
		  if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) kshad[j][d] = k;
		  }
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/

		/* 11.04.05 Left this comment block here (from inside_b_me()) */

		/* new ME block : */
		/* 05.30.05 Fixed a small bug here */
		/* The changes made to this section of code in the memory efficient
		 * banded implementation are the most complex changes necessary to 
		 * get memory efficiency.  The reason is because there are indices in 
		 * two other states for a B_st, y and z (instead of just y).  This
		 * means that when we're dealing with a dp_v that is d minus a v-state
		 * specific offset, we also have to worry about the y-state offset
		 * and z-state offset.
		 * Let's set kp as the equivalent of k from the old code, but
		 * now we have to take into account the offsets.  To remain as
		 * consistent as possible with the old code, we will keep the
		 * indexing in z the same in the recursion, and figure out what
		 * the corresponding indices involving state y are.  
		 * So the old recursion code is : 
		 *
		 * for (jp = 0; jp <= W; jp++) {
		 * j = i0-1+jp;
		 * for (d = 0; d <= jp; d++) 
		 * {
		 *   alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0]; *INIT*
		 *   if (ret_shadow != NULL) kshad[j][d] = 0;
		 *   for (k = 1; k <= d; k++)
		 *   *RECURSION*
		 *   if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) {
		 *     alpha[v][j][d] = sc;
		 *     if (ret_shadow != NULL) kshad[j][d] = k; }
		 * 
		 * So we'll minimally change alpha[z][j][k] to alpha[z][j][kp]
		 * The INIT may change because although alpha[z][j][0] MUST be
		 * within the bands (because dmin[z] >= 0), the corresponding
		 * cell in alpha[y] might not be within the bands for y.  
		 * That cell is alpha[y][j-dmin[z]-kp][d-dmin[y]-dmin[z]-kp]
		 * because k = kp + dmin[z] (it probably takes some time writing
		 * down the new and old equations, and staring and thinking for a 
		 * while - I would write down more here - but this is already pretty
		 * verbose ... ).
		 * 
		 * Therefore we can't just start with k (or kp)  = 0 
		 * (like the old code did), because that might not be valid.
		 *
		 * First we need to determine the smallest kp for which we can 
		 * do a valid traceback, which means the alpha cell for both the y
		 * state and z state are within the bands.  For a kp to be valid given
		 * the following code, the following three inequalities have to be
		 * true.
		 *
		 * (1) d-dmin[z]-kp <= dmax[y]  
		 * (2) d-dmin[z]-kp >= dmin[y]
		 * (3) kp <= dmax[z]-dmin[z]
		 *
		 * (1) and (2) need to be satisified to guarantee that the cell we
		 * are going to access in the alpha[y] deck is within the bands for
		 * state y.  (3) is necessary to guarantee that the cell we are
		 * going to access in the alpha[z] deck is within the bands for 
		 * state z.
		 * We can rearrange 1 and 2 : 
		 *
		 * (1) kp >= d-dmax[y]-dmin[z]
		 * (2) kp <= d-dmin[y]-dmin[z]
		 * 
		 * First to check to see if ANY kp is valid, we can first
		 * check to make sure that (d-dmin[y]-dmin[z]) (RHS of (2))
		 * is >= 0.  If not, then kp can never be 0 or greater. 
		 * So it can never be valid. So we check for this at
		 * the beginning.
		 * 
		 * So, to find the minimal kp that satisfies (1), (2) and (3)
		 * I set kp = d-dmax[y]-dmin[z], and then check that it kp >= 0
		 * If kp < 0, we set it to 0.  Then we check to make sure kp
		 * satisfies (3) (It has to satisfy (2) if it satisfies (1)
		 * because dmax[y] >= dmin[y]).  This is our *INIT* assignment.
		 * Next we incrementally step through all valid kp values, we'll need 
		 * a for loop with two conditions to check in the 'while' portion.  
		 * Namely, that kp satisfies inequalities (2) and (3), that is
		 * kp <= (d-dmin[y]-dmin[z]) and kp <= (dmax[z]-dmin[z])
		 * This is marked in the code by *RECUR*
		 *
		 * Also, we want to make sure the while statement from the 
		 * original for loop (non-banded) is also satisfied.  This
		 * statement is k <= d.  We're dealing with kp, and k = kp+dmin[z]
		 * so this statement becomes kp <= d-dmin[z].  However, inequality
		 * (2) (kp <= d-dmin[y]-dmin[z]) takes care of this because dmin[y] >= 0
		 * 
		 */

		dp_v = d - safe_hdmin[v];  /* d index for state v in alpha w/mem eff bands */
		dp_y = d - safe_hdmin[y];  /* d index for state y in alpha w/mem eff bands */
		dp_z = d - safe_hdmin[z];  /* d index for state z in alpha w/mem eff bands */

		/* First make sure we have any valid kp, we know from inequality (2)
		   that kp <= d-safe_hdmin[y]-safe_hdmin[z] so if this is < 0 then no kp
		   is valid (see notes above) */

		if((d-safe_hdmin[y]-safe_hdmin[z]) >= 0)
		{
		  if(j < safe_hdmax[y]) kp = d-safe_hdmin[z]-j;
		  else kp = d-safe_hdmin[z]-safe_hdmax[y];
		  if(kp < 0) kp = 0;
		  if(kp <= safe_hdmax[z] - safe_hdmin[z]) /* make sure its valid in deck alpha[z] */
		    {
		      alpha[v][j][dp_v] = alpha[y][j-safe_hdmin[z]-kp][d-safe_hdmin[y]-safe_hdmin[z]-kp] 
			+ alpha[z][j][kp];
		      if (ret_shadow != NULL) kshad[j][dp_v] = kp;
		      for (kp = kp+1; kp <= (d-safe_hdmin[y]-safe_hdmin[z]) && kp <= (safe_hdmax[z]-safe_hdmin[z]);
			   kp++)
			{
			  /*printf("v is %d | checking y : %d z : %d\n", v, y, z);
			  printf("y comp          : alpha[%d][%d][%d] is %f\n", y, (j-safe_hdmin[z]-kp),(d-safe_hdmin[y]-safe_hdmin[z]-kp), 
				 alpha[y][j-safe_hdmin[z]-kp][d-safe_hdmin[y]-safe_hdmin[z]-kp]);
			  printf("z comp          : alpha[%d][%d][%d] is %f\n", z, j, kp, alpha[z][j][kp]);
			  printf("existing v comp : alpha[%d][%d][%d] is %f\n", v, j, dp_v, alpha[v][j][dp_v]);
			  printf("\n");*/
			  /* the following if statement ensures that the alpha cell for 
			     state y and the cell for state z that we are about to query 
			     is in fact within the bands for state y and state z respectively*/
			  if ((sc = alpha[y][j-safe_hdmin[z]-kp][d-safe_hdmin[y]-safe_hdmin[z]-kp] 
			       + alpha[z][j][kp]) > alpha[v][j][dp_v]) 
			    {
			      alpha[v][j][dp_v] = sc;
			      if (ret_shadow != NULL) kshad[j][dp_v] = kp;
			    }
			}
		    }
		}
		else alpha[v][j][dp_v] = IMPOSSIBLE;
		/*else esl_fatal("cell in alpha matrix was not filled in due to bands.\n");*/
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 5 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //printf("3 j: %d\n", j);

	    /* CYK Full ME Bands used 6 */
	    /* Deleted because I realized this was no longer needed */

	    /* Bands used */
	    /* old line :	for (d = 2; d <= jp; d++) */
	    /* we assume hdmin[v][j-jmin[v]] >= 2 */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);

		y = cm->cfirst[v];
		/* CYK Full ME Bands used 7 block */
		/* original code block below : */
		/*
		alpha[v][j][d] = cm->endsc[v];  init w/ local end 
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-2] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/
		/* new ME code block : */
		dp_v = d - safe_hdmin[v]; /* d index for state v in alpha w/mem eff bands */
		/*printf("j: %d | v: %d | d: %d | dp_v: %d | safe_hdmin[v]: %d\n", j, v, d, dp_v, safe_hdmin[v]);*/
				
		alpha[v][j][dp_v] = cm->endsc[v]; /* init w/ local end */
		if(ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
							  in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-2) >= 0) && (((dp_y-2) < (j - (safe_hdmin[y+yoffset]) + 1))
					   && ((dp_y-2) < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j-1][dp_y-2] + cm->tsc[v][yoffset]) >  alpha[v][j][dp_v])
			  {
			    alpha[v][j][dp_v] = sc;
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 7 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //	    printf("3 j: %d\n", j);
	    /* CYK Full ME Bands used 8 */
	    /* Deleted because I realized this was no longer needed */

	    /* Bands used */
	    /* old line :	for (d = 1; d <= jp; d++) */
	    /* we assume safe_hdmin[v] >= 1 */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);

		y = cm->cfirst[v];
		/* CYK Full ME Bands used 9 block */
		/* original code block below : */
		/*
		alpha[v][j][d] = cm->endsc[v];
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  } 
		
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/
		/* new ME code block : */
		dp_v = d - safe_hdmin[v]; /* d index for state v in alpha w/mem eff bands */
		/*printf("v: %d | j: %d | dp_v: %d | j-jmin[v]: %d | safe_hdmin[%d]: %d | d: %d\n", v, j, dp_v, (j-jmin[v]), v, safe_hdmin[v], d);*/
		alpha[v][j][dp_v] = cm->endsc[v]; /* init w/ local end */
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-1) >= 0) && (((dp_y-1) < (j - (safe_hdmin[y+yoffset]) + 1))
					   && ((dp_y-1) < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j][dp_y-1] + cm->tsc[v][yoffset]) >  alpha[v][j][dp_v]) 
			  {
			    alpha[v][j][dp_v] = sc;
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  } 
		      }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 9 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //printf("3 j: %d\n", j);
	    /* CYK Full ME Bands used 10 */
	    /* Deleted because I realized this was no longer needed */

	    /* Bands used */
	    /* old line :	for (d = 1; d <= jp; d++) */
	    /* we assume safe_hdmin[v] >= 1 */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);

		y = cm->cfirst[v];
		/* CYK Full ME Bands used 11 block */
		/* original code block below : */
		/*
		alpha[v][j][d] = cm->endsc[v];
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/
		/* new ME code block : */
		dp_v = d - safe_hdmin[v]; /* d index for state v in alpha w/mem eff bands */
		alpha[v][j][dp_v] = cm->endsc[v]; /* init w/ local end */
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-1) >= 0) && (((dp_y-1) < (j - (safe_hdmin[y+yoffset]) + 1))
				      && ((dp_y-1) < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j-1][dp_y-1] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v])
			  {
			    alpha[v][j][dp_v] = sc;
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  }
		      }
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 11 end block */
	      }
	  }
	}				/* finished calculating deck v. */
      
      /* CYK Full ME Bands used 12 block */
      /* The following loops originally access alpha[v][j0][W] but the index W will be
	 in different positions due to the bands */

      /* ME Added the following two lines */
      Wp = W - safe_hdmin[v];
      /* We need to make sure that Wp is within the bands */
      if(Wp >= 0 && Wp <= (safe_hdmax[v] - safe_hdmin[v]))
	{
	  /* ME all subsequent changes in this block simply replace
	     W with Wp (so wherever Wp is, there used to be a W) */

	  /* Check for local begin getting us to the root.
	   * This is "off-shadow": if/when we trace back, we'll handle this
	   * case separately (and we'll know to do it because we'll immediately
	   * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	   * to jump right to state b; see below)
	   */
	  if (allow_begin && alpha[v][j0][Wp] + cm->beginsc[v] > bsc) 
	    {
	      b   = v;
	      bsc = alpha[v][j0][Wp] + cm->beginsc[v];
	    }

	  /* Check for whether we need to store an optimal local begin score
	   * as the optimal overall score, and if we need to put a flag
	   * in the shadow matrix telling insideT() to use the b we return.
	   */
	  if (allow_begin && v == 0 && bsc > alpha[0][j0][Wp]) {
	    alpha[0][j0][Wp] = bsc;
	    if (ret_shadow != NULL) yshad[j0][Wp] = USED_LOCAL_BEGIN;
	  }
	}
      /* CYK Full ME Bands used 12 end block */

      /* CYK Full ME Bands used 13 block */
      /* The following block implements the deck reuse strategy, however, here
	 we can't do that, because for each state, the bands are different, so 
	 we can't use old Decks, but rather must allocate a new one, and free
	 the old one. Lines specific to ME are indicated, and original lines
	 are commented out */

      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { 
	    /* Original code block : */
	    /* we can definitely release the S children of a bifurc. 
	       y = cm->cfirst[v]; deckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	       z = cm->cnum[v];   deckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	     End of original code block */
	    /* New ME code : */
	    y = cm->cfirst[v];
	    z = cm->cnum[v];  
	    free_vjd_deck(alpha[y], i0, j0);
	    alpha[y] = NULL;
	    free_vjd_deck(alpha[z], i0, j0);
	    alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) 
		  {
		    if (cm->sttype[y] == E_st) { 
		      nends--; 
		      /* Original code : if (nends == 0) { deckpool_push(dpool, end); end = NULL;} */
		      /* ME code deletes the previous line, we don't mess with end, because
			 it is used later */
		    } else 
		      /* original code (deck reuse) deckpool_push(dpool, alpha[y]);*/
		      /* new ME code : */
		      {
			//printf("calling free vjd deck for alpha[y=%d]\n", y);
			free_vjd_deck(alpha[y], i0, j0);
		      }
		      alpha[y] = NULL;
		  }
	      }
	  }
	/* CYK Full ME Bands used 13 end block */
      }
    } /* end loop over all v */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  
  /* CYK Full ME Bands used 14 */
  /* original line :  sc       = alpha[vroot][j0][W];*/
  Wp = W - safe_hdmin[vroot];
  sc       = alpha[vroot][j0][Wp];

  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = vroot; v <= vend; v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { deckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { deckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (deckpool_pop(dpool, &end)) free_vjd_deck(end, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  return sc;
}

/* Function: insideT_b_jd()
 *           EPN 05.24.05
 * *based on insideT(), only difference is memory efficient bands are used : 
 *
 * Date:     SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]
 *
 * Purpose:  Call inside, get vjd shadow matrix;
 *           then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 */
static float
insideT_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
	     int r, int z, int i0, int j0, 
	     int allow_begin, int *jmin, int *jmax,
	     int **hdmin, int **hdmax,
	     int *safe_hdmin, int *safe_hdmax)
{
  if(dsq == NULL)
    esl_fatal("ERROR(), dsq is NULL.\n");

  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;
  int       dp;                 /* add explanation */
  int       kp;                 /* add explanation */

  sc = inside_b_jd(cm, dsq, L, r, z, i0, j0, 
		   BE_EFFICIENT,	/* memory-saving mode */
		   NULL, NULL,	/* manage your own matrix, I don't want it */
		   NULL, NULL,	/* manage your own deckpool, I don't want it */
		   &shadow,		/* return a shadow matrix to me. */
		   allow_begin,      /* TRUE to allow local begins */
		   &b, &bsc,	/* if allow_begin is TRUE, gives info on optimal b */
		   jmin, jmax,    /* bands on j */
		   hdmin, hdmax,  /* j dependent bands on d */
		   safe_hdmin, safe_hdmax);

  pda = esl_stack_ICreate();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  /*printf("Starting traceback in insideT_b_me()\n");*/
  while (1) {
    /* CYK Full ME Bands used 15 */
    /* 2 lines below added */
    /*
    assert(j >= jmin[v]);
    assert(j <= jmax[v]);
    assert(d >= hdmin[v][j0]);
    assert(d <= hdmax[v][j0]);

    assert(d >= safe_hdmin[v]);
    assert(d <= safe_hdmax[v]);
    */
    if(cm->sttype[v] != EL_st && d > safe_hdmax[v])
      {
	printf("ERROR in insideT_b_jd(). d : %d > safe_hdmax[%d] (%d)\n", d, v, safe_hdmax[v]);
	exit(1);
      }
    if(cm->sttype[v] != EL_st && d < safe_hdmin[v])
      {
	printf("ERROR in insideT_b_jd(). d : %d < safe_hdmin[%d] (%d)\n", d, v, safe_hdmin[v]);
	exit(1);
      }

    /* Deal with end local states */
    if(cm->sttype[v] != EL_st)
      dp = d - safe_hdmin[v];
    else
      dp = d;
    
    if (cm->sttype[v] == B_st) {
      /* CYK Full ME Bands used 16 */
      /* original line : k = ((int **) shadow[v])[j][d];  */
      /* new 3 lines below replace it */
      assert(v >= 0);
      kp = ((int **) shadow[v])[j][dp];   /* kp = offset len of right fragment */
      z = cm->cnum[v];
      k = kp + safe_hdmin[z];  /* k = offset len of right fragment */
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */
      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      /* CYK Full ME Bands used 17 */
      /* original line : esl_stack_IPop(pda, &d); */
      esl_stack_IPop(pda, &dp);
      /* CYK Full ME Bands used 18 */
      /* line below added */

      /* Deal with end local states */
      if(cm->sttype[v] != EL_st)
	d = dp + safe_hdmin[y];
      else
      d = dp;

      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
    } else {
      yoffset = ((char **) shadow[v])[j][dp];

      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    esl_fatal("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	}
    }
  }
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  return sc;
}
#endif


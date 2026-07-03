/* CM_ALNDATA: data structure containing information relevant to the
 * alignment of a sequence and its output. Currently, mainly used by
 * cmalign, but also used in cmbuild with the --refine option.
 * 
 * Contents:
 *    1. The CM_ALNDATA object.
 *    2. Alignment workunit processing functions, which create and 
 *       fill CM_ALNDATA objects.
 *
 * EPN, Fri Jan  6 09:00:31 2012
 */
#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "easel.h"

#include "hmmer.h"

#include "infernal.h"

static int sub_alignment_prep(CM_t *orig_cm, char *errbuf, ESL_SQ *sq, CMSubMap_t **ret_submap, CM_t **ret_sub_cm);

/*****************************************************************
 * 1. The CM_ALNDATA object
 *****************************************************************/

/* Function:  cm_alndata_Create()
 * Synopsis:  Allocate a CM_ALNDATA object.
 * Incept:    EPN, Fri Jan  6 09:01:33 2012
 *            
 * Purpose:   Allocates a new <CM_ALNDATA> and returns a pointer
 *            to it.
 *
 * Throws:    <NULL> on allocation failure.
 */
CM_ALNDATA *
cm_alndata_Create(void)
{
  int status;
  CM_ALNDATA *data = NULL;

  ESL_ALLOC(data, sizeof(CM_ALNDATA));
  data->sq         = NULL;
  data->idx        = -1;
  data->tr         = NULL;
  data->sc         = 0.;
  data->pp         = 0.;
  data->ppstr      = NULL;
  data->spos       = -1;
  data->epos       = -1;
  data->secs_bands = 0.;
  data->secs_aln   = 0.;
  data->mb_tot     = 0.;
  data->tau        = -1.;
  
  return data;

 ERROR: 
  return NULL;
}

/* Function:  cm_alndata_Destroy()
 * Synopsis:  Free a CM_ALNDATA object.
 * Incept:    EPN, Fri Jan  6 09:10:55 2012
 *            
 * Purpose:  Frees a <CM_ALNDATA> object, but only frees the
 *           ESL_SQ <sq> if <free_sq> is TRUE. Often this is
 *           only a pointer to a sequence in another data 
 *           structure that will be free'd with that structure.
 *
 * Returns:  void.
 */
void
cm_alndata_Destroy(CM_ALNDATA *data, int free_sq)
{ 
  if(data == NULL) return;

  if(free_sq && data->sq != NULL) esl_sq_Destroy(data->sq);
  if(data->tr    != NULL)         FreeParsetree(data->tr);
  if(data->ppstr != NULL)         free(data->ppstr);
  free(data);

  return;
}

/*****************************************************************
 * 2. Alignment workunit processing functions
 *****************************************************************/

/* Function: sub_alignment_prep()
 * Date:     EPN, Mon Jan  9 05:25:26 2012
 *
 * Purpose:  Prepare for an alignment workunit in sub-mode.
 *
 * Args:     orig_cm    - the covariance model
 *           errbuf     - char buffer for reporting errors
 *           sq         - the sequence we're creating the sub CM for
 *           ret_submap - RETURN: the sub CM to original CM map, created here
 *           ret_sub_cm - RETURN: the sub CM, created here
 *
 * Returns:  eslOK on success;
 *           eslEMEM if we run out of memory;
 *           eslEINVAL on other error, errbuf is filled;
 *           <ret_dataA> is alloc'ed and filled with sq_block->count CM_ALNDATA objects.
 */
int
sub_alignment_prep(CM_t *orig_cm, char *errbuf, ESL_SQ *sq, CMSubMap_t **ret_submap, CM_t **ret_sub_cm)
{
  int          status;            /* easel status */
  CM_t        *sub_cm  = NULL;    /* the sub CM */
  CMSubMap_t  *submap  = NULL;    /* map from mother CM to sub CM, and vice versa */
  int          spos;              /* HMM node most likely to have emitted posn 1 of target seq */
  int          spos_state;        /* HMM state type for curr spos 0=match or 1=insert */
  int          epos;              /* HMM node most likely to have emitted posn L of target seq */
  int          epos_state;        /* HMM state type for curr epos 0=match or 1=insert */

  /* step 1. predict start and end positions (HMM nodes) from posterior matrix */
  if((status = cp9_Seq2Posteriors(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, sq->dsq, 1, sq->L, 0)) != eslOK) return status; 
  CP9NodeForPosn(orig_cm->cp9, 1, sq->L,     1, orig_cm->cp9_bmx, &spos, &spos_state, 0., TRUE,  0);
  CP9NodeForPosn(orig_cm->cp9, 1, sq->L, sq->L, orig_cm->cp9_bmx, &epos, &epos_state, 0., FALSE, 0);
  /* Deal with special cases for sub-CM alignment: If the most
   * likely state to have emitted the first or last residue is the
   * insert state in node 0, it only makes sense to start modelling
   * at consensus column 1. */
  if(spos == 0 && spos_state == 1) spos = 1;
  if(epos == 0 && epos_state == 1) epos = 1;
  /* If most-likely HMM node to emit final position comes BEFORE or
   * EQUALS the most-likely HMM node to emit first position, our HMM
   * alignment is crap, default to using the full CM. (note: If
   * EQUALS we could be right, but we can't build a CM from a single
   * consensus column (see notes in cm_modelmaker.c::cm_from_guide),
   * and I would argue we don't really care about getting single
   * residue alignments correct anyway. 
   */
  if(epos <= spos) { spos = 1; epos = orig_cm->cp9->M; } 
  
  /* step 2. build the sub_cm from the original CM. */
  if((status = build_sub_cm(orig_cm, errbuf, &sub_cm, 
			    spos, epos,                /* first and last col of structure kept in the sub_cm  */
			    &submap,                   /* this maps from the sub_cm to cm and vice versa      */
			    0)) != eslOK)              /* don't print debugging info */
    return status;

  /* step 3. configure the sub_cm */
  if((status = cm_ConfigureSub(sub_cm, errbuf, -1, orig_cm, submap)) != eslOK) return status; 

  *ret_sub_cm = sub_cm;
  *ret_submap = submap;

  return eslOK;
}

/* Function: DispatchSqBlockAlignment()
 * Date:     EPN, Fri Dec 30 14:59:43 2011
 *
 * Purpose:  Given a CM and a block of sequences, align the
 *           sequence(s) using the appropriate alignment function and
 *           return relevant data for eventual output in <ret_dataA>.
 *           This function simply calls DispatchSqAlignment() serially
 *           for each sequence in the block, and creates an array
 *           of the <ret_data> DispatchSqAlignment() returns.
 *
 *           Currently <mode>, <cp9b_valid> and <pass_idx> values sent
 *           to DispatchSqAlignment() are hard-coded to
 *           TRMODE_UNKNOWN, FALSE, and PLI_PASS_5P_AND_3P_FORCE (if
 *           cm->align_opts & CM_ALIGN_TRUNC) or PLI_PASS_STD_ANY (if
 *           (! cm->align_opts & CM_ALIGN_TRUNC)). This is because
 *           this function is only used by the alignment pipeline, in
 *           which these values are correct. If this changes, we may
 *           want caller to pass in an array of modes, cp9b_valids and
 *           pass_idx values, one per sq.
 *
 *           If (cm->flags & CM_ALIGN_XTAU) we'll potentially tighten
 *           HMM bands until the required DP matrices are below out
 *           limit (<mxsize>). cm->maxtau is the max allowed tau value
 *           during this iterative band tightening, and cm->xtau is
 *           the factor by which we multiply cm->tau at each iteration
 *           during band tightening.
 *
 * Args:     cm        - the covariance model
 *           errbuf    - char buffer for reporting errors
 *           sq_block  - block of sequences to align
 *           mxsize    - max size in Mb of allowable DP mx
 *           w         - stopwatch for timing individual stages
 *           w_tot     - stopwatch for timing total time per seq
 *           r         - RNG, req'd if CM_ALIGN_SAMPLE, can be NULL otherwise
 *           ret_dataA - RETURN: newly created array of CM_ALNDATA objects
 *
 * Returns:  eslOK on success;
 *           eslEINCOMPAT on contract violation, errbuf is filled;
 *           eslEMEM if we run out of memory;
 *           <ret_dataA> is alloc'ed and filled with sq_block->count CM_ALNDATA objects.
 */
int
DispatchSqBlockAlignment(CM_t *cm, char *errbuf, ESL_SQ_BLOCK *sq_block, float mxsize, ESL_STOPWATCH *w, 
			 ESL_STOPWATCH *w_tot, ESL_RANDOMNESS *r, CM_ALNDATA ***ret_dataA)
{
  int           status;          /* easel status */
  int           j;               /* counter over parsetrees */
  CM_ALNDATA  **dataA = NULL;    /* CM_ALNDATA array we'll create and return */
  ESL_SQ       *sqp;             /* ptr to a ESL_SQ */
  int           pass_idx;        /* pass_idx passed to DispatchSqAlignment() */
  char          mode;            /* mode passed to DispatchSqAlignment() */
  int           cp9b_valid;      /* passed to DispatchSqAlignment() */
  CM_P7_OM_HOLDER om_holder;     /* reusable LOCAL p7 profile/OPROFILE for the
				  * --p7pinbridge SW scan, built once and reused
				  * across this block (brief 090). Safe because
				  * this call owns the whole block. */

  ESL_ALLOC(dataA, sizeof(CM_ALNDATA *) * ESL_MAX(1, sq_block->count)); // avoid 0 malloc
  for(j = 0; j < sq_block->count; j++) dataA[j] = NULL;

  /* DispatchSqAligment() needs a mode, pipeline pass index, and
   * knowledge of whether cm->cp9b are valid for sequence to align
   * (see note in 'Purpose' above). Currently the relevant values 
   * for these are as follows:
   */
  mode       = TRMODE_UNKNOWN;
  pass_idx   = (cm->align_opts & CM_ALIGN_TRUNC) ? PLI_PASS_5P_AND_3P_FORCE : PLI_PASS_STD_ANY; 
  cp9b_valid = FALSE;

  /* main loop: for each sequence, call DispatchSqAlignment() to do the work */
  cm_p7_om_holder_Init(&om_holder);
  for(j = 0; j < sq_block->count; j++) {
    sqp = sq_block->list + j;
    if((status = DispatchSqAlignment(cm, errbuf, sqp, sq_block->first_seqidx + j, mxsize, mode, pass_idx, cp9b_valid, w, w_tot, r, &om_holder, &(dataA[j]))) != eslOK) { cm_p7_om_holder_Reset(&om_holder); goto ERROR; }
  }
  cm_p7_om_holder_Reset(&om_holder);
  *ret_dataA = dataA;

  return eslOK;

 ERROR: 
  if(dataA != NULL) { 
    for(j = 0; j < sq_block->count; j++) { 
      if(dataA[j] != NULL) cm_alndata_Destroy(dataA[j], FALSE);
    }
    free(dataA);
  }
  *ret_dataA = NULL;
  fprintf(stderr, "Problem during alignment of sequence %s\n", sqp->name);
  if(status == eslEMEM) ESL_FAIL(status, errbuf, "DispatchSqBlockAlignment(), out of memory");
  else return status; /* errbuf was filled by DispatchSqAlignment() */
}

/* Function: DispatchSqAlignment()
 * Date:     EPN, Thu Jan 12 14:47:26 2012
 *
 * Purpose:  Given a CM and a sequence, align the sequence(s) using
 *           the appropriate alignment function and return relevant
 *           data for eventual output in <ret_data>. 
 *
 *           This function can be called from either an alignment
 *           pipeline (i.e. cmalign) or a search/scan pipeline
 *           (i.e. cmsearch or cmscan). <idx> is the (overloaded) flag
 *           for determining which, if -1, we're a search/scan
 *           pipeline. This is only relevant because in a search/scan
 *           pipeline we don't care about determining spos/epos so we
 *           don't call ParsetreeToCMBounds().
 *                        
 *           If (cm->flags & CM_ALIGN_XTAU) we'll potentially tighten
 *           HMM bands until the required DP matrices are below out
 *           limit (<mxsize>). cm->maxtau is the max allowed tau value
 *           during this iterative band tightening, and cm->xtau is
 *           the factor by which we multiply cm->tau at each iteration
 *           during band tightening.
 *
 * Args:     cm         - the covariance model
 *           errbuf     - char buffer for reporting errors
 *           sq         - sequence to align
 *           idx        - index of sequence (may be used to reorder data later)
 *           mxsize     - max size in Mb of allowable DP mx 
 *           mode       - preset mode of alignment (TRMODE_UNKNOWN if unknown)
 *           pass_idx   - pipeline pass index, determines trunc penalty
 *           cp9b_valid - TRUE if cm->cp9b are valid, don't compute HMM bands
 *           w          - stopwatch for timing individual stages, can be NULL
 *           w_tot      - stopwatch for timing total time per seq, can be NULL
 *           r          - RNG, req'd if CM_ALIGN_SAMPLE, can be NULL otherwise
 *           om_holder  - reusable LOCAL p7 profile/OPROFILE holder for the
 *                        --p7pinbridge SW scan (brief 090); built once per
 *                        worker/block and reused across sequences. Can be NULL
 *                        (then the pinbridge wrapper builds/frees its own per
 *                        call). MUST be per-thread (not shared across threads).
 *           ret_data   - RETURN: newly created CM_ALNDATA object
 *
 * Returns:  eslOK on success;
 *           eslEINCOMPAT on contract violation, errbuf is filled;
 *           eslEMEM if we run out of memory;
 *           <ret_data> is alloc'ed and filled.
 */
int
DispatchSqAlignment(CM_t *cm, char *errbuf, ESL_SQ *sq, int64_t idx, float mxsize, char mode, int pass_idx,
		    int cp9b_valid, ESL_STOPWATCH *w, ESL_STOPWATCH *w_tot, ESL_RANDOMNESS *r,
		    CM_P7_OM_HOLDER *om_holder, CM_ALNDATA **ret_data)
{
  int           status;            /* easel status */
  CM_ALNDATA   *data         = NULL; /* CM_ALNDATA we'll create and fill */
  float         sc           = 0.;   /* score from alignment function */
  float         pp           = 0.;   /* average PP from alignment function */
  Parsetree_t  *tr           = NULL; /* ptr to a parsetree */
  char         *ppstr        = NULL; /* ptr to a PP string */
  float         secs_bands   = 0.;   /* seconds elapsed for band calculation */
  float         secs_aln     = 0.;   /* seconds elapsed for alignment calculation */
  float         mb_tot       = 0.;   /* size of all DP matrices used for alignment */
  double        tau          = -1.;  /* tau used for calculating bands */
  float         thresh1      = -1.;  /* cp9b->thresh1 used for calculating bands */
  float         thresh2      = -1.;  /* cp9b->thresh2 used for calculating bands */
  int           spos         = -1;   /* start posn: first non-gap CM consensus position */
  int           epos         = -1;   /* end   posn: final non-gap CM consensus position */
  double        save_tau     = cm->tau; /* cm->tau upon entrance, we restore before leaving */
  float         save_thresh1 = (cm->cp9b == NULL) ? -1. : cm->cp9b->thresh1;
  float         save_thresh2 = (cm->cp9b == NULL) ? -1. : cm->cp9b->thresh2;

  /* alignment options */
  int do_nonbanded = (cm->align_opts & CM_ALIGN_NONBANDED) ? TRUE  : FALSE;
  int do_qdb       = (cm->align_opts & CM_ALIGN_QDB)       ? TRUE  : FALSE;
  int do_hbanded   = (do_nonbanded || do_qdb)              ? FALSE : TRUE;
  int do_optacc    = (cm->align_opts & CM_ALIGN_OPTACC)    ? TRUE  : FALSE;
  int do_sample    = (cm->align_opts & CM_ALIGN_SAMPLE)    ? TRUE  : FALSE;
  int do_post      = (cm->align_opts & CM_ALIGN_POST)      ? TRUE  : FALSE;
  int do_sub       = (cm->align_opts & CM_ALIGN_SUB)       ? TRUE  : FALSE;
  int do_small     = (cm->align_opts & CM_ALIGN_SMALL)     ? TRUE  : FALSE;
  int do_trunc     = (cm->align_opts & CM_ALIGN_TRUNC)     ? TRUE  : FALSE;
  int do_xtau      = (cm->align_opts & CM_ALIGN_XTAU)      ? TRUE  : FALSE;
  int do_p7band    = (cm->align_opts & CM_ALIGN_P7BANDED)  ? TRUE  : FALSE;
  int doing_search = FALSE;
  /* Brief 120: IBV HMM-divergence fallback. Set when cm_TrAlignHB / cm_AlignHB
   * fails on IBV-derived bands and we've already rebuilt with vitband for
   * this sequence; prevents infinite retry.
   */
  int ibv_fallback_used = FALSE;

#if eslDEBUGLEVEL >= 1
  printf("#DEBUG: in DispatchSqAlignment() %s\n", sq->name);
  printf("#DEBUG: \tdo_nonbanded: %d\n", do_nonbanded);
  printf("#DEBUG: \tdo_optacc:    %d\n", do_optacc);
  printf("#DEBUG: \tdo_sample:    %d\n", do_sample);
  printf("#DEBUG: \tdo_post:      %d\n", do_post);
  printf("#DEBUG: \tdo_sub:       %d\n", do_sub);
  printf("#DEBUG: \tdo_small:     %d\n", do_small);
  printf("#DEBUG: \tdo_trunc:     %d\n", do_trunc);
  printf("#DEBUG: \tdo_qdb:       %d\n", do_qdb);
  printf("#DEBUG: \tdoing_search: %d\n", doing_search);
#endif
  
  /* sub-mode specific variables (wouldn't be needed if sub mode were not supported) */
  CM_t        *orig_cm = cm;      /* pointer to the original CM */
  CM_t        *sub_cm  = NULL;    /* the sub CM */
  CMSubMap_t  *submap  = NULL;    /* map from mother CM to sub CM, and vice versa */
  Parsetree_t *full_tr = NULL;    /* converted parsetree to full CM */

  /* contract check */
  if(do_small  && do_hbanded)       ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do small and HMM banded alignment");
  if(do_small  && do_optacc)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do small and opt acc alignment");
  if(do_post   && do_small)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do PP and small alignment");
  if(do_optacc && do_sample)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to sample and do optacc alignment");
  if(do_sub    && do_small)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do sub and small alignment");
  if(do_sub    && do_trunc)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do sub and truncated alignment");
  if(do_sample && r == NULL)        ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to sample but RNG r == NULL");
  if(do_xtau   && ! do_hbanded)     ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to multiply tau without HMM banded alignment");
  if(do_xtau   && cp9b_valid)       ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to multiply tau but HMM bands already valid");
  if(do_qdb    && do_nonbanded)     ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do qdb and nonbanded alignment");
  if(do_qdb    && do_trunc)         ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs and truncated alignment");
  /* qdb + trunc combo disallowed only b/c no function exists for it yet */
  if(do_qdb    && (! do_small))     ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs but not divide and conquer");
  /* qdb + small combo disallowed b/c only non-HMM banded non-small alignment functions are not set up to use QDBs */
  if(do_qdb && cm->qdbinfo == NULL) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs but cm->qdbinfo is NULL");
  }
  if(do_qdb && (cm->qdbinfo->dmin2 == NULL || cm->qdbinfo->dmax2 == NULL)) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to use qdbs but cm->qdbinfo is NULL");
  }
  if(do_trunc && (! cm_pli_PassAllowsTruncation(pass_idx))) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() trying to do truncated alignment, but pass_idx doesn't allow truncation (PLI_PASS_STD_ANY)");
  }
  if(pass_idx == PLI_PASS_STD_ANY && (mode == TRMODE_L || mode == TRMODE_R || mode == TRMODE_T)) { 
    ESL_XFAIL(eslEINCOMPAT, errbuf, "DispatchSqAlignment() mode is L, R, or T, but pass_idx is PLI_PASS_STD_ANY");
  }

  if(w_tot != NULL) esl_stopwatch_Start(w_tot);

  /* do sub-mode specific pre-alignment steps, if nec */
  if(do_sub) { 
    if((status = sub_alignment_prep(cm, errbuf, sq, &submap, &sub_cm)) != eslOK) goto ERROR;
    cm = sub_cm;
  }

  if(w != NULL) esl_stopwatch_Start(w);
  /* do small D&C alignment, if nec */
  if(do_small) { 
    if(do_trunc) { 
      sc = TrCYK_DnC(cm, sq->dsq, sq->L, 0, 1, sq->L, pass_idx, FALSE, &tr); /* FALSE: don't reproduce 1.0 behavior */
      mb_tot = 4. * CYKNonQDBSmallMbNeeded(cm, sq->L); /* not sure how accurate this is */
    }
    else { 
      /* with QDB, always use dmin2/dmax2, the looser of the two sets of QDBs in cm->qdbinfo */
      sc = CYKDivideAndConquer(cm, sq->dsq, sq->L, 0, 1, sq->L, &tr, 
			       (do_qdb) ? cm->qdbinfo->dmin2 : NULL, 
			       (do_qdb) ? cm->qdbinfo->dmax2 : NULL);
      mb_tot = CYKNonQDBSmallMbNeeded(cm, sq->L);
    }
  }
  else { /* do_small is FALSE */
    if(do_nonbanded || do_qdb) { /* do not use HMM bands */
      if(do_trunc) { 
	if((status = cm_TrAlignSizeNeeded(cm, errbuf, sq->L, mxsize, do_sample, do_post, 
					  NULL, NULL, NULL, &mb_tot)) != eslOK) goto ERROR;
	if((status = cm_TrAlign(cm, errbuf, sq->dsq, sq->L, mxsize, mode, pass_idx, 
				do_optacc, do_sample, cm->trnb_mx, cm->trnb_shmx, cm->trnb_omx, 
				cm->trnb_emx, r, do_post ? &ppstr : NULL, &tr, NULL, &pp, &sc)) != eslOK) goto ERROR;
      }
      else {
	if((status = cm_AlignSizeNeeded(cm, errbuf, sq->L, mxsize, do_sample, do_post, 
					NULL, NULL, NULL, &mb_tot)) != eslOK) goto ERROR;
	if((status = cm_Align(cm, errbuf, sq->dsq, sq->L, mxsize, do_optacc, do_sample, cm->nb_mx, cm->nb_shmx, 
			      cm->nb_omx, cm->nb_emx, r, do_post ? &ppstr : NULL, &tr, &pp, &sc)) != eslOK) goto ERROR;
      }
    }
    else { /* use HMM bands */
      if(! cp9b_valid) {
	/* TODO #9 mitigation: --p7band produces too-narrow k-envelopes for small-M
	 * models (M < ~200), causing accuracy regression on rmark4e (Lacto-usp, atoC,
	 * snoZ152, ar45, SNORA47). Fall back to unbanded CP9 F/B for small CMs; the
	 * absolute wall savings from --p7band on tiny CMs is negligible. See brief 059. */
#define P7BAND_MIN_M 0   /* was 200; gate was introduced in brief 062 (session 19) for
                            the pinbridge-era band-derivation path. The current F+B IBV
                            (briefs 116-126) is a different algorithm; gate removed per
                            brief 141. */
	if(do_p7band && cm->fp7 != NULL && cm->fp7->M < P7BAND_MIN_M) {
	  fprintf(stderr, "#P7BAND_SKIP M=%d reason=small_M_acc_gap (threshold=%d)\n",
		  cm->fp7->M, P7BAND_MIN_M);
	  do_p7band = FALSE;
	}
	if(do_p7band && cm->fp7 != NULL) {
	  /* p7-derived bands: Viterbi trace -> kmin/kmax -> banded CP9 F/B -> CM bands */
	  P7_PROFILE *gm_p7b  = NULL;
	  P7_GMX     *gx_p7b  = NULL;
	  P7_BG      *bg_p7b  = NULL;
	  P7_TRACE   *tr_p7b  = NULL;
	  int        *p7_kmin  = NULL;
	  int        *p7_kmax  = NULL;
	  int        *p7_i2k   = NULL;
	  int         p7_ncells = 0;

	  /* Create P7 objects */
	  bg_p7b = p7_bg_Create(cm->abc);
	  gm_p7b = p7_profile_Create(cm->fp7->M, cm->abc);
	  /* Configure profile for band derivation:
	   * - Standard (non-truncated) alignment: GLOCAL (must align full model)
	   * - Truncated alignment: T-profile (LOCAL + forced full-sequence parse).
	   *   The T-profile allows entry/exit at any model node but forces the parse
	   *   to include the first and last residue (N->N and C->C set to -inf).
	   */
	  if(do_trunc) {
	    p7_ProfileConfig(cm->fp7, bg_p7b, gm_p7b, sq->L, p7_LOCAL);
	    p7_ProfileConfig5PrimeAnd3PrimeTrunc(gm_p7b, sq->L);
	  } else {
	    p7_ProfileConfig(cm->fp7, bg_p7b, gm_p7b, sq->L, p7_GLOCAL);
	  }
	  /* gx_p7b (full O(M*L) p7 matrix) is allocated lazily only where the
	   * unbanded p7_Seq2BandsVit path actually needs it (brief 094). The
	   * --p7pinbridge success path never touches it, so we avoid the eager
	   * full-matrix alloc that (a) defeats pinbridge's large-M memory win and
	   * (b) overflows int32 in p7_gmx_Create at M=L ~ 1.5e5 (e.g. HSV). */
	  tr_p7b = p7_trace_Create();

	  /* Build local nodepad copy with p7bpad (p7padplus) added */
	  int *local_nodepad = NULL;
	  if(cm->flags & CMH_P7NODEPAD) {
	    int k;
	    ESL_ALLOC(local_nodepad, sizeof(int) * (cm->fp7->M + 1));
	    for(k = 0; k <= cm->fp7->M; k++) local_nodepad[k] = cm->p7_cm_nodepad[k] + cm->p7bpad;
	  }

	  /* Derive p7 bands: either SW-pinbridge prefilter + banded Viterbi (--p7pinbridge)
	   * or full unbanded p7_GViterbi (default). Pinbridge wrapper falls back to
	   * full p7_Seq2BandsVit if the banded trace fails inside the prefilter band
	   * (signaled by ncells=0). Pinbridge is correct for both truncated and
	   * non-truncated alignment.
	   *
	   * M-threshold gate: pinbridge has per-sequence overhead (build OPROFILE,
	   * SSE scan setup, LSIS allocation, banded traceback) that on small models
	   * exceeds the cost of full p7_GViterbi at O(LM). For M < 200 the full
	   * Viterbi wins; gate pinbridge on M >= 200 to capture the big-M speedup
	   * without the tiny-M tail regressions. */
	  struct timespec _ta_p7b, _tb_p7b;
	  const char *_p7b_kind = NULL;
	  clock_gettime(CLOCK_MONOTONIC, &_ta_p7b);
	  if (cm->p7_use_ibv) {
	    _p7b_kind = "p7ibv";
	    /* F+B direct-band band derivation (brief 120). Does NOT take gm/gx
	     * because it extracts transitions directly from cm->fp7. We still
	     * built gm/gx above for the vitband fallback path.
	     *
	     * brief 124: with --p7ibv-mem, dispatch to the divide-and-conquer
	     * O(M*logL) band deriver, byte-identical to the flat path but with
	     * dramatically lower peak memory at large M/L.
	     */
	    if (cm->p7_ibv_wv) {
	      /* Brief 169: windowed-Viterbi band = MAP-trace i2k +/- F+B-halfwidth
	       * pad (cm->p7_wv_nodepad, calibrated once at align-time setup). */
	      _p7b_kind = "p7ibv-wv";
	      int *wv_nodepad = NULL;
	      int  wk;
	      if (cm->p7_wv_nodepad == NULL)
		ESL_FAIL(eslEINVAL, errbuf, "--p7ibv-wv: cm->p7_wv_nodepad not calibrated");
	      ESL_ALLOC(wv_nodepad, sizeof(int) * (cm->fp7->M + 1));
	      for (wk = 0; wk <= cm->fp7->M; wk++) wv_nodepad[wk] = cm->p7_wv_nodepad[wk] + cm->p7bpad;
	      status = p7_Seq2BandsWV(cm, errbuf, sq->dsq, sq->L, wv_nodepad,
				      do_trunc, /* brief 171 */
				      &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	      free(wv_nodepad);
	    } else if (cm->p7_ibv_mem) {
	      _p7b_kind = "p7ibv-dnc";
	      status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, sq->L,
					   cm->p7_ibv_delta, cm->p7_ibv_base_slab,
					   TRUE, /* do_boundary_widen: CM-side preserves current behavior */
					   FALSE, /* brief 172: do_kband (unbanded; --p7ibv-mem keeps exact delta band) */
					   do_trunc, /* brief 171 */
					   cm->p7_ibv_mode, cm->p7_ibv_width, /* brief 140 */
					   &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    } else {
	      status = p7_Seq2BandsIBV(cm, errbuf, sq->dsq, sq->L,
				       cm->p7_ibv_delta, do_trunc, /* brief 171 */
				       cm->p7_ibv_mode, cm->p7_ibv_width, /* brief 140 */
				       &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    }
	    /* No internal ncells==0 fallback here: IBV always produces a band.
	     * Empty rows default to [1, M] inside the kernel.
	     */
	  } else if (cm->p7_use_kmeranchor) {
	    /* Brief 026: k-mer best-window anchor. Blind diagonal-dominance guide
	     * deriver feeding the unmodified p7_pins2bands_nodepad. Opt-in. */
	    _p7b_kind = "kmeranchor";
	    status = p7_Seq2BandsKmerAnchor(cm, errbuf, sq->dsq, sq->L, local_nodepad,
	                                    &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    /* ncells==0 => no usable anchor; fall back to full unbanded Viterbi band
	     * derivation (same shape as the pinbridge ncells==0 fallback below). */
	    if (status == eslOK && p7_ncells == 0) {
	      _p7b_kind = "kmeranchor->vitband";
	      if (gx_p7b == NULL) gx_p7b = p7_gmx_Create(cm->fp7->M, sq->L);
	      status = p7_Seq2BandsVit(errbuf, gm_p7b, gx_p7b, bg_p7b, tr_p7b,
	                               sq->dsq, sq->L, cm->p7bpad, local_nodepad,
	                               0, 0, &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    }
	  } else if (cm->p7_use_pinbridge) {
	    _p7b_kind = "pinbridge";
	    status = p7_Seq2BandsPinBridgeWrap(cm, errbuf, gm_p7b, bg_p7b, tr_p7b,
					       sq->dsq, sq->L, cm->p7bpad,
					       local_nodepad,
					       0, 0, /* hopback=0, vitend=0 */
					       om_holder,
					       &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    /* If pinbridge couldn't produce a trace (rare; band missed the trace
	     * entirely), fall back to full unbanded Viterbi for this sequence. */
	    if (status == eslOK && p7_ncells == 0) {
	      if (gx_p7b == NULL) gx_p7b = p7_gmx_Create(cm->fp7->M, sq->L);
	      status = p7_Seq2BandsVit(errbuf, gm_p7b, gx_p7b, bg_p7b, tr_p7b,
				       sq->dsq, sq->L, cm->p7bpad,
				       local_nodepad,
				       0, 0,
				       &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	    }
	  } else {
	    _p7b_kind = "vitband";
	    if (gx_p7b == NULL) gx_p7b = p7_gmx_Create(cm->fp7->M, sq->L);
	    status = p7_Seq2BandsVit(errbuf, gm_p7b, gx_p7b, bg_p7b, tr_p7b,
				     sq->dsq, sq->L, cm->p7bpad,
				     local_nodepad,
				     0, 0, /* hopback=0, vitend=0 */
				     &p7_i2k, &p7_kmin, &p7_kmax, &p7_ncells);
	  }
	  clock_gettime(CLOCK_MONOTONIC, &_tb_p7b);
	  {
	    double _p7b_s = (_tb_p7b.tv_sec - _ta_p7b.tv_sec) +
	                    (_tb_p7b.tv_nsec - _ta_p7b.tv_nsec) / 1e9;
	    fprintf(stderr, "#P7BAND_TIME %s kind=%s L=%d M=%d t=%.6f ncells=%d ibvmode=%d ibvwidth=%d ibvdelta=%d\n",
	            sq->name, _p7b_kind, (int)sq->L, cm->fp7->M, _p7b_s,
	            p7_ncells, cm->p7_ibv_mode, cm->p7_ibv_width, cm->p7_ibv_delta); /* brief 140: ncells = band-size discriminator */
	  }

	  /* Debug: report Viterbi band stats */
	  if(status == eslOK && p7_ncells > 0) {
	    int dbg_npin = 0, dbg_minpin = sq->L+1, dbg_maxpin = 0;
	    int dbg_i;
	    for(dbg_i = 1; dbg_i <= sq->L; dbg_i++) {
	      if(p7_i2k[dbg_i] != -1) {
		dbg_npin++;
		if(dbg_i < dbg_minpin) dbg_minpin = dbg_i;
		if(dbg_i > dbg_maxpin) dbg_maxpin = dbg_i;
	      }
	    }
	    float dbg_avgbw = (float)p7_ncells / (float)sq->L;
	    fprintf(stderr, "#P7BAND %s L=%d M=%d npins=%d pin_range=[%d,%d] ncells=%d avg_bw=%.1f\n",
		    sq->name, (int)sq->L, cm->fp7->M, dbg_npin, dbg_minpin, dbg_maxpin, p7_ncells, dbg_avgbw);
	  }
	  else {
	    fprintf(stderr, "#P7BAND %s L=%d M=%d FAILED (status=%d ncells=%d)\n",
		    sq->name, (int)sq->L, cm->fp7->M, status, p7_ncells);
	  }

	  if(status == eslOK && p7_ncells > 0) {
	    /* Use p7 bands to derive CM bands via p7-banded CP9 F/B with tau-ratcheting */
	    struct timespec _ta_cp9, _tb_cp9;
	    clock_gettime(CLOCK_MONOTONIC, &_ta_cp9);
	    status = cp9_IterateSeq2BandsP7B(cm, errbuf, sq->dsq, sq->L, p7_kmin, p7_kmax,
					     1, sq->L, pass_idx, mxsize,
					     doing_search, do_sample, do_post,
					     cm->maxtau, 0, 0, NULL);
	    clock_gettime(CLOCK_MONOTONIC, &_tb_cp9);
	    double _cp9_s = (_tb_cp9.tv_sec - _ta_cp9.tv_sec) + (_tb_cp9.tv_nsec - _ta_cp9.tv_nsec)/1e9;
	    fprintf(stderr, "#P7PB_POST M=%d L=%d cp9_iterate=%.4f\n", cm->fp7->M, (int)sq->L, _cp9_s);
	  }
	  else {
	    /* Viterbi found no path or error; fall back to standard cp9 bands */
	    status = eslERANGE;
	  }

	  /* Free p7 objects */
	  if(local_nodepad) free(local_nodepad);
	  if(p7_i2k)  free(p7_i2k);
	  if(p7_kmin) free(p7_kmin);
	  if(p7_kmax) free(p7_kmax);
	  p7_trace_Destroy(tr_p7b);
	  p7_gmx_Destroy(gx_p7b);
	  p7_profile_Destroy(gm_p7b);
	  p7_bg_Destroy(bg_p7b);

	  if(status != eslOK) {
	    /* P7B bands too wide even at maxtau; fall back to standard cp9 band derivation */
	    /* Brief 159: the standard fallback below builds a *non-banded* full CP9
	     * F/B matrix, (L+1)*(M+1) cells. At genome scale (M~L~2e5) that alone is
	     * hundreds of GB, far over --mxsize, and unlike the CM DP matrix there is
	     * no banding/tau lever to shrink the CP9 F/B. The p7-banded path already
	     * returned !eslOK (matrix too big), so attempting the even-larger
	     * non-banded CP9 F/B would only OOM (it did: brief-157 SIGSEGV via int32
	     * overflow, now int64-correct but still 155 GB -> SIGKILL). If the
	     * non-banded CP9 F/B (fwd+bck) would itself exceed --mxsize, refuse with
	     * eslERANGE up front -- the eslERANGE-before-alloc convention used for
	     * every CM DP matrix (cm_mx.c). Only fires in this already-over-budget
	     * p7band fallback, so plain non-p7band cmalign is unaffected (its CP9 F/B
	     * is intentionally not mxsize-gated, matching cm_*AlignSizeNeededHB). */
	    float cp9fb_Mb = 2.0 * (float) SizeNeededCP9Matrix(sq->L, cm->cp9->M, NULL, NULL);
	    if(cp9fb_Mb > mxsize)
	      ESL_XFAIL(eslERANGE, errbuf,
			"non-banded CP9 F/B band derivation needs %.1f > %.1f Mb limit.\nUse --mxsize, --maxtau or --tau (this seq needs a p7-banded/--ckpt path).",
			cp9fb_Mb, (float) mxsize);
	    if(do_xtau) {
	      if((status = cp9_IterateSeq2Bands(cm, errbuf, sq->dsq, 1, sq->L, pass_idx, mxsize, doing_search, do_sample, do_post, 1,
						cm->maxtau, NULL)) != eslOK) goto ERROR;
	    }
	    else {
	      if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq,
					 1, sq->L, cm->cp9b, doing_search, pass_idx, 0)) != eslOK) goto ERROR;
	    }
	  }
	}
	else if(do_xtau) { /* multiply tau (if nec) until required mx is below Mb limit (mxsize) */
	  if((status = cp9_IterateSeq2Bands(cm, errbuf, sq->dsq, 1, sq->L, pass_idx, mxsize, doing_search, do_sample, do_post, 1 /*do_iterate*/,
					    cm->maxtau, NULL)) != eslOK) goto ERROR;
	}
	else {
	  if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq, 
				     1, sq->L, cm->cp9b, doing_search, pass_idx, 0)) != eslOK) goto ERROR;
	}
	if(w != NULL) esl_stopwatch_Stop(w);
	secs_bands = (w == NULL) ? 0. : w->elapsed;
	tau     = cm->tau;
	thresh1 = cm->cp9b->thresh1;
	thresh2 = cm->cp9b->thresh2;
	/* note: we don't set these three if cp9b_valid is TRUE */
      }

      /* DIAGNOSTIC: always print CP9 banded matrix cell count (whether
       * cykbands is on or off) so we can compare CP9-only vs cykbands runs. */
      {
        double _cp9_cells = 0.;
        CP9Bands_t *_cp9b = cm->cp9b;
        int _v, _jp;
        for(_v = 0; _v < cm->M; _v++)
          for(_jp = 0; _jp <= _cp9b->jmax[_v] - _cp9b->jmin[_v]; _jp++)
            if(hd_min(_cp9b, _v, _jp) <= hd_max(_cp9b, _v, _jp))
              _cp9_cells += hd_max(_cp9b, _v, _jp) - hd_min(_cp9b, _v, _jp) + 1;
        fprintf(stderr, "#P7PB_POST M=%d L=%d cp9_band_cells=%.0f\n",
                (cm->fp7 ? cm->fp7->M : 0), (int)sq->L, _cp9_cells);
      }

      /* CYK pre-pass: run CYK on CP9 bands, derive tighter per-state bands before Inside/Outside */
      if(cm->p7_use_cykbands) {
	struct timespec _ta_cyk, _tb_cyk;
	clock_gettime(CLOCK_MONOTONIC, &_ta_cyk);

	/* measure original band area (total d-band cells across all states) */
	double _orig_cells = 0.;
	{
	  CP9Bands_t *_cp9b = cm->cp9b;
	  int _v, _jp;
	  for(_v = 0; _v < cm->M; _v++)
	    for(_jp = 0; _jp <= _cp9b->jmax[_v] - _cp9b->jmin[_v]; _jp++)
	      if(hd_min(_cp9b, _v, _jp) <= hd_max(_cp9b, _v, _jp))
		_orig_cells += hd_max(_cp9b, _v, _jp) - hd_min(_cp9b, _v, _jp) + 1;
	}

	Parsetree_t *_cyk_tr = NULL;
	float        _cyk_sc  = 0.;
	int          _cyk_ok  = FALSE;

	if(do_trunc) {
	  char _cyk_mode = TRMODE_UNKNOWN;
	  float _cyk_avgpp = 0.;
	  if(cm_TrAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize,
			  TRMODE_UNKNOWN, pass_idx, FALSE/*do_optacc*/, FALSE/*do_sample*/,
			  cm->trhb_mx, cm->trhb_shmx, NULL/*post_mx*/, NULL/*emit_mx*/,
			  NULL/*r*/, NULL/*ret_ppstr*/,
			  &_cyk_tr, &_cyk_mode, &_cyk_avgpp, &_cyk_sc) == eslOK && _cyk_tr != NULL) {
	    _cyk_ok = TRUE;
	  }
	} else {
	  if(cm_alignT_hb(cm, errbuf, sq->dsq, sq->L, mxsize, FALSE/*do_optacc*/,
			  cm->hb_mx, cm->hb_shmx, NULL/*emit_mx*/,
			  &_cyk_tr, &_cyk_sc) == eslOK && _cyk_tr != NULL) {
	    _cyk_ok = TRUE;
	  }
	}

	double _tight_cells = _orig_cells;
	if(_cyk_ok) {
	  /* Per-state pad was archived 2026-05-19 (see cm_CYKPerstatePadCompute
	   * doc comment for failure analysis). Production uses uniform pad. */
	  if(cm_BandsFromCYKParsetree(cm, errbuf, _cyk_tr,
				      1, sq->L, cm->p7_cykbands_pad, NULL, FALSE,
				      cm->cp9b, pass_idx, 0) == eslOK) {
	    _tight_cells = 0.;
	    CP9Bands_t *_cp9b = cm->cp9b;
	    int _v, _jp;
	    for(_v = 0; _v < cm->M; _v++)
	      for(_jp = 0; _jp <= _cp9b->jmax[_v] - _cp9b->jmin[_v]; _jp++)
		if(hd_min(_cp9b, _v, _jp) <= hd_max(_cp9b, _v, _jp))
		  _tight_cells += hd_max(_cp9b, _v, _jp) - hd_min(_cp9b, _v, _jp) + 1;
	  }
	  FreeParsetree(_cyk_tr);
	}

	clock_gettime(CLOCK_MONOTONIC, &_tb_cyk);
	double _cyk_s   = (_tb_cyk.tv_sec  - _ta_cyk.tv_sec)  + (_tb_cyk.tv_nsec  - _ta_cyk.tv_nsec) /1e9;
	double _b_ratio = (_orig_cells > 0.) ? _tight_cells / _orig_cells : 1.0;
	fprintf(stderr, "#P7PB_POST M=%d L=%d cyk_prepass=%.4f band_area_ratio=%.4f orig_cells=%.0f tight_cells=%.0f\n",
		(cm->fp7 ? cm->fp7->M : 0), (int)sq->L, _cyk_s, _b_ratio, _orig_cells, _tight_cells);
      }

      /* --dump-bands: write per-(v,j) band TSV before cm_AlignHB/cm_TrAlignHB */
      if(cm->p7_dump_bands_file != NULL && cm->cp9b != NULL) {
        FILE *_dbfp = fopen(cm->p7_dump_bands_file, "w");
        if(_dbfp != NULL) {
          CP9Bands_t *_dbands = cm->cp9b;
          int _dv, _dj, _djp;
          fprintf(_dbfp, "v\tstate_type\tj\thdmin\thdmax\tdwidth\n");
          for(_dv = 0; _dv < cm->M; _dv++) {
            for(_djp = 0; _djp <= _dbands->jmax[_dv] - _dbands->jmin[_dv]; _djp++) {
              _dj = _dbands->jmin[_dv] + _djp;
              int _dmin = hd_min(_dbands, _dv, _djp);
              int _dmax = hd_max(_dbands, _dv, _djp);
              int _dwidth = (_dmax >= _dmin) ? (_dmax - _dmin + 1) : 0;
              fprintf(_dbfp, "%d\t%s\t%d\t%d\t%d\t%d\n",
                      _dv, Statetype(cm->sttype[_dv]), _dj, _dmin, _dmax, _dwidth);
            }
          }
          fclose(_dbfp);
        }
      }

      if(w != NULL) esl_stopwatch_Start(w);
	  struct timespec _ta_cm, _tb_cm;
	  clock_gettime(CLOCK_MONOTONIC, &_ta_cm);
	CM_ALIGN_HB_RETRY:
	  if(do_trunc) {
		/* brief 126 merge: keep cd577024's #DBG-009 instrumentation, but route
		 * SizeNeededHB failure to CM_ALIGN_HB_CHECK_FB (IBV vitband fallback)
		 * instead of directly to ERROR, so the brief-120 IBV fallback stays live
		 * in the trunc path. For non-IBV runs CHECK_FB falls through to ERROR. */
		status = cm_TrAlignSizeNeededHB(cm, errbuf, sq->L, mxsize, do_sample, do_post,
					    NULL, NULL, NULL, NULL, NULL, &mb_tot);
		fprintf(stderr, "#DBG-009 trunc SizeNeededHB status=%d mb_tot=%.2f mxsize=%.2f do_post=%d errbuf=[%s]\n",
			status, mb_tot, (float) mxsize, do_post, errbuf);
		if(status != eslOK) goto CM_ALIGN_HB_CHECK_FB;
		/* checkpointed sqrt(M)-memory TRUNCATED OptAcc path: engaged by --ckpt
		 * (CM_ALIGN_CHECKPT) for the global, pure-MATL-chain (bps=0) OptAcc case it
		 * supports (marginal modes J/L/R, T absent); stock cm_TrAlignHB() otherwise
		 * (byte-identical output, but full-cube memory). On failure fall through to
		 * CM_ALIGN_HB_CHECK_FB (IBV vitband fallback), not ERROR. */
		int do_trckpt = ((cm->align_opts & CM_ALIGN_CHECKPT) && do_optacc && (! do_sample) &&
				 cm_CheckptTrAlignHB_Qualifies(cm)) ? TRUE : FALSE;
		if(do_trckpt) {
		  status = cm_CheckptTrAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, mode, pass_idx,
					       cm->trhb_emx, do_post ? &ppstr : NULL, &tr, NULL, &pp, &sc);
		}
		else {
	      	  status = cm_TrAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, mode, pass_idx,
				    do_optacc, do_sample, cm->trhb_mx, cm->trhb_shmx, cm->trhb_omx,
				    cm->trhb_emx, r, do_post ? &ppstr : NULL, &tr, NULL, &pp, &sc);
		}
      }
      else {
	if((status = cm_AlignSizeNeededHB(cm, errbuf, sq->L, mxsize, do_sample, do_post,
					  NULL, NULL, NULL, NULL, NULL, &mb_tot)) != eslOK) goto CM_ALIGN_HB_CHECK_FB;
	/* checkpointed sqrt(M)-memory OptAcc path: engaged by --ckpt (CM_ALIGN_CHECKPT)
	 * only for the non-truncated, global, pure-MATL-chain OptAcc case it supports;
	 * stock cm_AlignHB() otherwise (byte-identical output, but full-cube memory).
	 * On failure fall through to CM_ALIGN_HB_CHECK_FB (IBV vitband fallback), not ERROR. */
	int do_checkpt = ((cm->align_opts & CM_ALIGN_CHECKPT) && do_optacc && (! do_sample) &&
			  cm_CheckptAlignHB_Qualifies(cm)) ? TRUE : FALSE;
	if(do_checkpt) {
	  status = cm_CheckptAlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, cm->hb_emx,
				     do_post ? &ppstr : NULL, &tr, &pp, &sc);
	}
	else {
	  status = cm_AlignHB(cm, errbuf, sq->dsq, sq->L, mxsize, do_optacc, do_sample, cm->hb_mx, cm->hb_shmx,
			      cm->hb_omx, cm->hb_emx, r, do_post ? &ppstr : NULL, &tr, &pp, &sc);
	}
      }
    CM_ALIGN_HB_CHECK_FB:
      /* Brief 120 IBV HMM-divergence fallback: cm_TrInsideAlignHB returns
       * eslEAMBIGUOUS "no valid parsetree found" on the 3/14 brief-117 seqs
       * where HMM-Viterbi disagrees with the CM's preferred parse (brief 117
       * §5). Re-derive bands using the unbanded p7_Seq2BandsVit path (PAD=20
       * around the p7-Viterbi trace) and retry the CM alignment once.
       */
      if (status != eslOK && cm->p7_use_ibv && !ibv_fallback_used) {
	fprintf(stderr, "#P7IBV_FALLBACK %s L=%d M=%d status=%d errbuf='%s'\n",
		sq->name, (int)sq->L, (cm->fp7 ? cm->fp7->M : 0), status, errbuf);
	ibv_fallback_used = TRUE;
	errbuf[0] = '\0';
	status    = eslOK;

	/* Self-contained vitband redo: build temporary p7 profile/gx/bg/trace,
	 * call p7_Seq2BandsVit, push through cp9_IterateSeq2BandsP7B, free.
	 */
	{
	  P7_PROFILE *fb_gm = p7_profile_Create(cm->fp7->M, cm->abc);
	  P7_BG      *fb_bg = p7_bg_Create(cm->abc);
	  P7_GMX     *fb_gx = p7_gmx_Create(cm->fp7->M, sq->L);
	  P7_TRACE   *fb_tr = p7_trace_Create();
	  int        *fb_i2k = NULL, *fb_kmin = NULL, *fb_kmax = NULL;
	  int         fb_ncells = 0;
	  int        *fb_nodepad = NULL;
	  int         fbk;

	  if (do_trunc) {
	    p7_ProfileConfig(cm->fp7, fb_bg, fb_gm, sq->L, p7_LOCAL);
	    p7_ProfileConfig5PrimeAnd3PrimeTrunc(fb_gm, sq->L);
	  } else {
	    p7_ProfileConfig(cm->fp7, fb_bg, fb_gm, sq->L, p7_GLOCAL);
	  }
	  if (cm->flags & CMH_P7NODEPAD) {
	    fb_nodepad = (int *) malloc(sizeof(int) * (cm->fp7->M + 1));
	    if (fb_nodepad) {
	      for (fbk = 0; fbk <= cm->fp7->M; fbk++)
		fb_nodepad[fbk] = cm->p7_cm_nodepad[fbk] + cm->p7bpad;
	    }
	  }
	  status = p7_Seq2BandsVit(errbuf, fb_gm, fb_gx, fb_bg, fb_tr,
				   sq->dsq, sq->L, cm->p7bpad, fb_nodepad,
				   0, 0,
				   &fb_i2k, &fb_kmin, &fb_kmax, &fb_ncells);
	  if (status == eslOK && fb_ncells > 0) {
	    status = cp9_IterateSeq2BandsP7B(cm, errbuf, sq->dsq, sq->L,
					     fb_kmin, fb_kmax,
					     1, sq->L, pass_idx, mxsize,
					     doing_search, do_sample, do_post,
					     cm->maxtau, 0, 0, NULL);
	  }
	  if (fb_i2k) free(fb_i2k);
	  if (fb_kmin) free(fb_kmin);
	  if (fb_kmax) free(fb_kmax);
	  if (fb_nodepad) free(fb_nodepad);
	  p7_trace_Destroy(fb_tr);
	  p7_gmx_Destroy(fb_gx);
	  p7_profile_Destroy(fb_gm);
	  p7_bg_Destroy(fb_bg);
	}
	if (status != eslOK) goto ERROR;
	goto CM_ALIGN_HB_RETRY;
      }
      if (status != eslOK) goto ERROR;
      clock_gettime(CLOCK_MONOTONIC, &_tb_cm);
      double _cm_s = (_tb_cm.tv_sec - _ta_cm.tv_sec) + (_tb_cm.tv_nsec - _ta_cm.tv_nsec)/1e9;
      fprintf(stderr, "#P7PB_POST M=%d L=%d cm_align_hb=%.4f%s\n",
              (cm->fp7 ? cm->fp7->M : 0), (int)sq->L, _cm_s,
              ibv_fallback_used ? " ibv_fallback=1" : "");
    }
  }

  if(do_sub) { /* add size of original CM's CP9 matrices used for calculating start/end position */
    mb_tot += orig_cm->cp9_mx->size_Mb;
    mb_tot += orig_cm->cp9_bmx->size_Mb;
  }

  if(w != NULL) esl_stopwatch_Stop(w);
  secs_aln = (w == NULL) ? 0. : w->elapsed;

  if(do_sub) { 
    /* convert sub cm parsetree to a full CM parsetree */
    if((status = sub_cm2cm_parsetree(orig_cm, cm, &full_tr, tr, submap, 0)) != eslOK) ESL_XFAIL(status, errbuf, "out of memory, converting sub parsetree to full parsetree");
    /* free sub data structures, we're done with them */
    FreeParsetree(tr);   tr     = full_tr;
    FreeCM(cm);          cm     = orig_cm;
    FreeSubMap(submap);  submap = NULL;
  }
  
  /* determine start and end points of the parsetree, 
   * but only if we're not in a search/scan pipeline 
   */
  if(idx != -1) { /* we're not in a search/scan pipeline */
    if((status = ParsetreeToCMBounds(cm, tr, TRUE, TRUE, errbuf, NULL, NULL, NULL, NULL, &spos, &epos)) != eslOK) goto ERROR;
  }
  
  /* create and fill data */
  ESL_ALLOC(data, sizeof(CM_ALNDATA));
  data->sq         = sq;
  data->idx        = idx;
  data->tr         = tr;
  data->sc         = sc;
  data->pp         = (do_post)      ? pp     : 0.;
  data->ppstr      = (do_post)      ? ppstr  : NULL;
  data->spos       = spos;
  data->epos       = epos;
  data->secs_bands = (do_nonbanded) ? 0.     : secs_bands;
  data->secs_aln   = secs_aln;
  data->mb_tot     = mb_tot;
  data->tau        = tau;
  data->thresh1    = thresh1;
  data->thresh2    = thresh2;
  if(w_tot != NULL) esl_stopwatch_Stop(w_tot);
  data->secs_tot   = (w_tot == NULL) ? 0. : w_tot->elapsed;

  *ret_data = data;

  cm->tau = save_tau;
  if(cm->cp9b != NULL) { 
    cm->cp9b->thresh1 = save_thresh1;
    cm->cp9b->thresh2 = save_thresh2;
  }
  return eslOK;

 ERROR: 
  cm->tau = save_tau;
  if(cm->cp9b != NULL) { 
    cm->cp9b->thresh1 = save_thresh1;
    cm->cp9b->thresh2 = save_thresh2;
  }
  if(data != NULL) cm_alndata_Destroy(data, FALSE);
  *ret_data = NULL;

  if(status == eslEMEM) ESL_FAIL(status, errbuf, "DispatchSqAlignment(), out of memory");

  return status; 
}

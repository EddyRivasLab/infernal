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
  data->sqp        = NULL;
  data->idx        = -1;
  data->tr         = NULL;
  data->sc         = 0.;
  data->pp         = 0.;
  data->ppstr      = NULL;
  data->cm_from    = -1;
  data->cm_to      = -1;
  data->secs_bands = 0.;
  data->secs_aln   = 0.;
  
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
cm_alndata_Destroy(CM_ALNDATA *data, int free_sqp)
{ 
  if(free_sqp && data->sqp != NULL) esl_sq_Destroy(data->sqp);
  if(data->tr    != NULL)           FreeParsetree(data->tr);
  if(data->ppstr != NULL)           free(data->ppstr);
  free(data);
}

/*****************************************************************
 * 2. Alignment workunit processing functions
 *****************************************************************/

/* Function: sub_alignment_workunit_prep()
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
sub_alignment_workunit_prep(CM_t *orig_cm, char *errbuf, ESL_SQ *sq, CMSubMap_t **ret_submap, CM_t **ret_sub_cm)
{
  int          status;            /* easel status */
  CM_t        *sub_cm  = NULL;    /* the sub CM */
  CMSubMap_t  *submap  = NULL;    /* map from mother CM to sub CM, and vice versa */
  int          spos;              /* HMM node most likely to have emitted posn 1 of target seq */
  int          spos_state;        /* HMM state type for curr spos 0=match or 1=insert */
  int          epos;              /* HMM node most likely to have emitted posn L of target seq */
  int          epos_state;        /* HMM state type for curr epos 0=match or 1=insert */

  /* sub-specific step 1. predict start and end positions (HMM nodes) from posterior matrix */
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
  
  /* sub-specific step 2. build the sub_cm from the original CM. */
  if((status = build_sub_cm(orig_cm, errbuf, &sub_cm, 
			    spos, epos,                /* first and last col of structure kept in the sub_cm  */
			    &submap,                   /* this maps from the sub_cm to cm and vice versa      */
			    0)) != eslOK)              /* don't print debugging info */
    return status;

  /* sub-specific step 3. configure the sub_cm */
  if((status = cm_ConfigureSub(sub_cm, errbuf, -1, orig_cm, submap)) != eslOK) return status; /* FALSE says: don't calculate W, we won't need it */

  *ret_sub_cm = sub_cm;
  *ret_submap = submap;

  return eslOK;
}

/* Function: ProcessAlignmentWorkunit()
 * Date:     EPN, Fri Dec 30 14:59:43 2011
 *
 * Purpose:  Given a CM and a block of sequences, align the sequences 
 *           using the appropriate alignment function and return 
 *           relevant data for eventual output in <ret_dataA>.
 *
 * Args:     cm              - the covariance model
 *           errbuf          - char buffer for reporting errors
 *           sq_block        - block of sequences to align
 *           mxsize          - max size in Mb of allowable DP mx
 *           w               - stopwatch for timing individual stages
 *           w_tot           - stopwatch for timing total time per seq
 *           ret_dataA       - array of CM_ALNDATA objects to return
 *
 * Returns:  eslOK on success;
 *           eslEINCOMPAT on contract violation, errbuf is filled;
 *           eslEMEM if we run out of memory;
 *           <ret_dataA> is alloc'ed and filled with sq_block->count CM_ALNDATA objects.
 */
int
ProcessAlignmentWorkunit(CM_t *cm, char *errbuf, ESL_SQ_BLOCK *sq_block, float mxsize, ESL_STOPWATCH *w, ESL_STOPWATCH *w_tot, CM_ALNDATA ***ret_dataA)
{
  int           status;          /* easel status */
  int           j;               /* counter over parsetrees */
  CM_ALNDATA **dataA = NULL;    /* CM_ALNDATA array we'll create */
  ESL_SQ       *sqp;             /* ptr to a ESL_SQ */
  float         sc, pp;          /* score, average PP, from alignment function */
  Parsetree_t  *tr = NULL;       /* ptr to a parsetree */
  char         *ppstr = NULL;    /* ptr to a PP string */
  TruncOpts_t  *tro = NULL;      /* truncated alignment info, remains NULL if --notrunc */
  float         secs_bands = 0.; /* seconds elapsed for band calculation */
  float         secs_aln;        /* seconds elapsed for alignment calculation */
  int   cfrom_emit, cto_emit;    /* CM boundaries of aligned sequence */
  int   cfrom_span, cto_span;    /* CM boundaries of parsetree */

  /* alignment options */
  int do_nonbanded = (cm->align_opts & CM_ALIGN_NONBANDED) ? TRUE : FALSE;
  int do_optacc    = (cm->align_opts & CM_ALIGN_OPTACC)    ? TRUE : FALSE;
  int do_post      = (cm->align_opts & CM_ALIGN_POST)      ? TRUE : FALSE;
  int do_sample    = (cm->align_opts & CM_ALIGN_SAMPLE)    ? TRUE : FALSE;
  int do_sub       = (cm->align_opts & CM_ALIGN_SUB)       ? TRUE : FALSE;
  int do_small     = (cm->align_opts & CM_ALIGN_SMALL)     ? TRUE : FALSE;
  int do_trunc     = (cm->align_opts & CM_ALIGN_TRUNC)     ? TRUE : FALSE;
  int doing_search = FALSE;

  /* non-banded truncated matrices, used only if --nonbanded */
  CM_TR_MX        *trmx   = NULL;
  CM_TR_SHADOW_MX *trshmx = NULL;
  CM_TR_MX        *tromx  = NULL;
  CM_TR_EMIT_MX   *tremx  = NULL;

  /* non-banded matrices, used only if --nonbanded & --notrunc */
  CM_MX        *mx   = NULL;
  CM_SHADOW_MX *shmx = NULL;
  CM_MX        *omx  = NULL;
  CM_EMIT_MX   *emx  = NULL;

  /* sub-mode specific variables (wouldn't be needed if sub mode were not supported) */
  CM_t        *orig_cm = cm;      /* pointer to the original CM */
  CM_t        *sub_cm  = NULL;    /* the sub CM */
  CMSubMap_t  *submap  = NULL;    /* map from mother CM to sub CM, and vice versa */
  Parsetree_t *full_tr  = NULL;   /* converted parsetree to full CM */

  /* contract check */
  if(do_small && (! do_nonbanded)) ESL_FAIL(eslEINCOMPAT, errbuf, "ProcessAlignmentWorkunit() trying to do small and HMM banded alignment");
  if(do_post  && do_small)         ESL_FAIL(eslEINCOMPAT, errbuf, "ProcessAlignmentWorkunit() trying to do PP and small alignment");
  if(do_sub   && do_small)         ESL_FAIL(eslEINCOMPAT, errbuf, "ProcessAlignmentWorkunit() trying to do sub and small alignment");
  if(do_sub   && do_trunc)         ESL_FAIL(eslEINCOMPAT, errbuf, "ProcessAlignmentWorkunit() trying to do sub and truncated alignment");
  if(sq_block->count <= 0)         ESL_FAIL(eslEINCOMPAT, errbuf, "ProcessAlignmentWorkunit() received empty block");

  if(do_trunc) { 
    tro = CreateTruncOpts();
    tro->allowL = tro->allowR = TRUE;
  }
  /* create non-banded matrices, if nec (if do_sub: matrices will be built per-sequence) */
  if(do_nonbanded && (! do_sub) && (! do_small)) { 
    if(do_trunc) { 
      if((trmx   = cm_tr_mx_Create(cm))        == NULL) goto ERROR;
      if((tromx  = cm_tr_mx_Create(cm))        == NULL) goto ERROR;
      if((trshmx = cm_tr_shadow_mx_Create(cm)) == NULL) goto ERROR;
      if((tremx  = cm_tr_emit_mx_Create(cm))   == NULL) goto ERROR;
    }
    else { 
      if((mx     = cm_mx_Create(cm))           == NULL) goto ERROR;
      if((omx    = cm_mx_Create(cm))           == NULL) goto ERROR;
      if((shmx   = cm_shadow_mx_Create(cm))    == NULL) goto ERROR;
      if((emx    = cm_emit_mx_Create(cm))      == NULL) goto ERROR;
    }
  }

  ESL_ALLOC(dataA, sizeof(CM_ALNDATA *) * sq_block->count);
  for(j = 0; j < sq_block->count; j++) dataA[j] = NULL;

  /* main loop: for each sequence, create a parsetree */
  for(j = 0; j < sq_block->count; j++) { 
    if(w_tot != NULL) esl_stopwatch_Start(w_tot);
    sqp = sq_block->list + j;

    /* do sub-mode specific pre-alignment steps, if nec */
    if(do_sub) { 
      if((status = sub_alignment_workunit_prep(cm, errbuf, sqp, &submap, &sub_cm)) != eslOK) return status;
      if(do_nonbanded) { /* create the non-banded matrices */
	if((mx   = cm_mx_Create(sub_cm))        == NULL) goto ERROR;
	if((omx  = cm_mx_Create(sub_cm))        == NULL) goto ERROR;
	if((shmx = cm_shadow_mx_Create(sub_cm)) == NULL) goto ERROR;
	if((emx  = cm_emit_mx_Create(sub_cm))   == NULL) goto ERROR;
      }
      cm = sub_cm;
    }

    if(w != NULL) esl_stopwatch_Start(w);
    /* do small D&C alignment, if nec */
    if(do_small) { 
      if(do_trunc) sc = TrCYK_DnC          (cm, sqp->dsq, sqp->L, 0, 1, sqp->L, &tr);
      else         sc = CYKDivideAndConquer(cm, sqp->dsq, sqp->L, 0, 1, sqp->L, &tr, NULL, NULL);
    }
    else { /* do_small is FALSE */
      if(do_nonbanded) { /* do not use HMM bands */
	if(do_trunc) { 
	  status = cm_TrAlign(cm, errbuf, sqp->dsq, sqp->L, mxsize, TRMODE_UNKNOWN, do_optacc, do_sample, 
			      trmx, trshmx, tromx, tremx, NULL, &ppstr, &tr, NULL, &pp, &sc);
	}
	else {
	  status = cm_Align(cm, errbuf, sqp->dsq, sqp->L, mxsize, do_optacc, do_sample, 
			    mx, shmx, omx, emx, NULL, &ppstr, &tr, &pp, &sc);
	}
	if(status != eslOK) return status;
      }
      else { /* use HMM bands */
	if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sqp->dsq, 
				   1, sqp->L, cm->cp9b, doing_search, tro, 0)) != eslOK) return status;
	if(w != NULL) esl_stopwatch_Stop(w);
	secs_bands = (w == NULL) ? 0. : w->elapsed;
	
	if(w != NULL) esl_stopwatch_Start(w);
	if(do_trunc) { 
	  status = cm_TrAlignHB(cm, errbuf, sqp->dsq, sqp->L, mxsize, TRMODE_UNKNOWN, do_optacc, do_sample, 
				cm->trhbmx, cm->trshhbmx, cm->trohbmx, cm->trehbmx, NULL, &ppstr, &tr, NULL, &pp, &sc);
	}
	else { 
	  status = cm_AlignHB(cm, errbuf, sqp->dsq, sqp->L, mxsize, do_optacc, do_sample, 
			      cm->hbmx, cm->shhbmx, cm->ohbmx, cm->ehbmx, NULL, &ppstr, &tr, &pp, &sc);
	}
	if(status != eslOK) return status;
      }
    }
    if(w != NULL) esl_stopwatch_Stop(w);
    secs_aln = (w == NULL) ? 0. : w->elapsed;

    if(do_sub) { 
      /* convert sub cm parsetree to a full CM parsetree */
      if((status = sub_cm2cm_parsetree(orig_cm, cm, &full_tr, tr, submap, 0)) != eslOK) ESL_FAIL(status, errbuf, "out of memory, converting sub parsetree to full parsetree");
      /* free sub data structures, we're done with them */
      FreeParsetree(tr);   tr     = full_tr;
      FreeCM(cm);          cm     = orig_cm;
      FreeSubMap(submap);  submap = NULL;
      if(mx     != NULL) cm_mx_Destroy(mx);
      if(shmx   != NULL) cm_shadow_mx_Destroy(shmx);
      if(omx    != NULL) cm_mx_Destroy(omx);
      if(emx    != NULL) cm_emit_mx_Destroy(emx);
    }

    /* determine start and end points of the parsetree */
    if((status = ParsetreeToCMBounds(cm, tr, errbuf, &cfrom_span, &cto_span, &cfrom_emit, &cto_emit)) != eslOK) return status;

    /* create and fill dataA[j] */
    ESL_ALLOC(dataA[j], sizeof(CM_ALNDATA));
    dataA[j]->sqp        = sqp;
    dataA[j]->idx        = sq_block->first_seqidx + j;
    dataA[j]->tr         = tr;
    dataA[j]->sc         = sc;
    dataA[j]->pp         = (do_post)      ? pp     : 0.;
    dataA[j]->ppstr      = (do_post)      ? ppstr  : NULL;
    dataA[j]->cm_from    = cfrom_emit;
    dataA[j]->cm_to      = cto_emit;
    dataA[j]->secs_bands = (do_nonbanded) ? 0.     : secs_bands;
    dataA[j]->secs_aln   = secs_aln;
    if(w_tot != NULL) esl_stopwatch_Stop(w_tot);
    dataA[j]->secs_tot   = (w_tot == NULL) ? 0. : w_tot->elapsed;
  }
  *ret_dataA = dataA;

  /* clean up */
  if(! do_sub) { 
    if(mx     != NULL) cm_mx_Destroy(mx);
    if(shmx   != NULL) cm_shadow_mx_Destroy(shmx);
    if(omx    != NULL) cm_mx_Destroy(omx);
    if(emx    != NULL) cm_emit_mx_Destroy(emx);
  }
  if(trmx   != NULL) cm_tr_mx_Destroy(trmx);
  if(trshmx != NULL) cm_tr_shadow_mx_Destroy(trshmx);
  if(tromx  != NULL) cm_tr_mx_Destroy(tromx);
  if(tremx  != NULL) cm_tr_emit_mx_Destroy(tremx);
  return eslOK;

 ERROR: 
  if(! do_sub) { 
    if(mx    != NULL) cm_mx_Destroy(mx);
    if(shmx  != NULL) cm_shadow_mx_Destroy(shmx);
    if(omx   != NULL) cm_mx_Destroy(omx);
    if(emx   != NULL) cm_emit_mx_Destroy(emx);
  }
  if(trmx   != NULL) cm_tr_mx_Destroy(trmx);
  if(trshmx != NULL) cm_tr_shadow_mx_Destroy(trshmx);
  if(tromx  != NULL) cm_tr_mx_Destroy(tromx);
  if(tremx  != NULL) cm_tr_emit_mx_Destroy(tremx);
  if(dataA != NULL) { 
    for(j = 0; j < sq_block->count; j++) cm_alndata_Destroy(dataA[j], FALSE);
    free(dataA);
  }

  ESL_FAIL(status, errbuf, "ProcessAlignmentWorkunit(), out of memory");
  return status; /* NOT REACHED */
}

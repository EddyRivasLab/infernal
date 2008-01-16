/*
 * EPN, Wed Nov 22 12:00:20 2006
 * Ported from RSEARCH-1.1
 *
 * stats.c
 *
 * Routines for calculating statistics in Infernal.  
 *
 * Robert J. Klein
 * Separate file started May 17, 2002
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <string.h>

#include "easel.h"
#include "esl_gumbel.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* Function: AllocCMStats()
 * Date:     EPN, Wed May  2 14:10:25 2007      
 *
 * Purpose:  Allocate a CMStats_t data structure given
 *           the number of partitions and relevant info needed
 *           for creating the HMMFilter_t objects.
 * 
 * Args:     np - number of partitions 
 */
CMStats_t *
AllocCMStats(int np)
{
  int status;
  CMStats_t  *cmstats;
  int i, p;

  ESL_ALLOC(cmstats, sizeof(struct cmstats_s));

  cmstats->np = np;
  ESL_ALLOC(cmstats->ps, sizeof(int) * cmstats->np);
  ESL_ALLOC(cmstats->pe, sizeof(int) * cmstats->np);
  ESL_ALLOC(cmstats->gumAA, sizeof(GumbelInfo_t **) * GUM_NMODES);
  ESL_ALLOC(cmstats->hfiA,  sizeof(HMMFilterInfo_t *) * FTHR_NMODES);
  for(i = 0; i < GUM_NMODES; i++) {
    ESL_ALLOC(cmstats->gumAA[i], sizeof(GumbelInfo_t *) * cmstats->np);
    for(p = 0; p < cmstats->np; p++) {
      cmstats->gumAA[i][p] = CreateGumbelInfo();
      if(cmstats->gumAA[i][p] == NULL) goto ERROR; /* memory error */
    }
  }
  for(i = 0; i < FTHR_NMODES; i++) {
    cmstats->hfiA[i] = CreateHMMFilterInfo();
    if(cmstats->hfiA[i] == NULL) goto ERROR; /* memory error */
  }
  return cmstats;

 ERROR:
  cm_Fail("AllocCMStats() memory allocation error.");
  return NULL; /* never reached */
}

/* Function: FreeCMStats()
 * Returns: (void) 
 */
void 
FreeCMStats(CMStats_t *cmstats)
{
  int i, p;
  for(i = 0; i < GUM_NMODES; i++)
    {
      for(p = 0; p < cmstats->np; p++)
	if(cmstats->gumAA[i][p] != NULL) free(cmstats->gumAA[i][p]);
      free(cmstats->gumAA[i]);
    }
  free(cmstats->gumAA);
  for(i = 0; i < FTHR_NMODES; i++)
    FreeHMMFilterInfo(cmstats->hfiA[i]);
  free(cmstats->hfiA);
  free(cmstats->ps);
  free(cmstats->pe);
  free(cmstats);
}  


/* Function: debug_print_cmstats
 */
int debug_print_cmstats(CMStats_t *cmstats, int has_fthr)
{
  int p;
  printf("Num partitions: %d\n", cmstats->np);
  for (p = 0; p < cmstats->np; p++)
    {
      printf("Partition %d: start: %d end: %d\n", p, cmstats->ps[p], cmstats->pe[p]);
      printf("cm_lc Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CM_LC][p]);
      printf("cm_gc Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CM_GC][p]);
      printf("cm_li Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CM_LI][p]);
      printf("cm_gi Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CM_GI][p]);
      printf("cp9_lv Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CP9_LV][p]);
      printf("cp9_gv Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CP9_GV][p]);
      printf("cp9_lf Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CP9_LF][p]);
      printf("cp9_gf Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[GUM_CP9_GF][p]);
      printf("\n\n");
    }

  if(has_fthr)
    {
      printf("Filter CM_LC info:\n");
      DumpHMMFilterInfo(cmstats->hfiA[FTHR_CM_LC]);
      printf("Filter CM_LI info:\n");
      DumpHMMFilterInfo(cmstats->hfiA[FTHR_CM_LI]);
      printf("Filter CM_GC info:\n");
      DumpHMMFilterInfo(cmstats->hfiA[FTHR_CM_GC]);
      printf("Filter CM_GI info:\n");
      DumpHMMFilterInfo(cmstats->hfiA[FTHR_CM_GI]);
      printf("\n\n");
    }
  return eslOK;
}

/* Function: debug_print_gumbelinfo
 */
int debug_print_gumbelinfo(GumbelInfo_t *gum)
{
  if(gum->is_valid) printf("N: %d L: %d lambda: %.5f mu: %.5f\n", gum->N, gum->L, gum->lambda, gum->mu);
  else             printf("invalid (not yet set)\n");
  return eslOK;
}

/* Function: get_gc_comp
 * Date:     EPN, Tue Aug  7 08:51:06 2007
 * Purpose:  Given a sequence and start and stop coordinates, returns 
 *           integerized GC composition of the region. This Easelfied
 *           version replaces RSEARCH's version.
 */
int get_gc_comp(ESL_SQ *sq, int start, int stop) 
{
  /* contract check */
  if(sq->abc == NULL)
    cm_Fail("get_gc_comp expects sq to have a valid alphabet.");
  if(! (sq->flags & eslSQ_DIGITAL))
    cm_Fail("get_gc_comp expects sq to be digitized");
  if(sq->abc->type != eslRNA && sq->abc->type != eslDNA)
    cm_Fail("get_gc_comp expects alphabet of RNA or DNA");

  int status;
  int i;
  float *ct = NULL;
  float  gc;

  if (start > stop) {
    i = start;
    start = stop;
    stop = i;
  }
  ESL_ALLOC(ct, sizeof(float) * sq->abc->K);
  esl_vec_FSet(ct, sq->abc->K, 0.);
  for (i = start; i <= stop; i++)
  {
    if(esl_abc_XIsGap(sq->abc, sq->dsq[i])) cm_Fail("in get_gc_comp, res %d is a gap!\n", i);
    esl_abc_FCount(sq->abc, ct, sq->dsq[i], 1.);
  }
  gc = ct[sq->abc->inmap[(int) 'G']] + ct[sq->abc->inmap[(int) 'C']];
  gc /= ((float) (stop-start+1));
  gc *= 100.;
  free(ct);
  return ((int) gc);

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0; /* never reached */
}

/*
 * Function: GetDBInfo()
 *
 * Date:     Easelification: EPN, Thu Dec  7 06:07:58 2006
 *           (initial - RSEARCH::get_dbinfo()) RJK, Thu Apr 11, 2002 [St. Louis]
 *
 * Purpose:  Given a sequence file name, determine the total size of the
 *           seqs in the file (DB) and GC content information.
 *
 * Args:     abc      - alphabet for sq file
 *           sqfp     - open sequence file
 *           ret_N    - RETURN: total length (residues) or all seqs in seqfile
 *           gc_ct    - RETURN: gc_ct[x] observed 100-nt segments with GC% of x [0..100] 
 *
 * Dies on parse error of sqfile.
 */
void GetDBInfo (const ESL_ALPHABET *abc, ESL_SQFILE *sqfp, long *ret_N, double **ret_gc_ct) 
{
  int               status;
  ESL_SQ           *sq;
  int               i, j, jp;  
  long              N = 0;
  double           *gc_ct;
  int               gc;
  int               all_ambig_flag; /* used to check if curr DB chunk is all ambiguous characters 
				     * usually Ns, if it is, we don't count it towards the GC content info */
  /*printf("in GetDBInfo\n");*/

  ESL_ALLOC(gc_ct, sizeof(double) * GC_SEGMENTS);
  for (i=0; i<GC_SEGMENTS; i++)
    gc_ct[i] = 0.;

  if(ret_gc_ct != NULL) sq = esl_sq_CreateDigital(abc); 
  else                  sq = esl_sq_Create(); /* allows abc to be NULL if we don't care about gc stats */
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) 
    { 
      N += sq->n;
      /*printf("new N: %d\n", N);*/
      if(ret_gc_ct != NULL) 
	{
	  for(i = 1; i <= sq->n; i += 100)
	    {
	      j = (i+99 <= sq->n) ? i+99 : sq->n;
	      gc = get_gc_comp(sq, i, j);
	      /*printf(">%d.raw\n", i);*/
	      all_ambig_flag = TRUE;
	      for(jp = 0; jp < 100 && (jp+i) < sq->n; jp++) {
		if(sq->dsq[i+jp] < abc->K) {
		  all_ambig_flag = FALSE; 
		  break; 
		}
	      }
	      /*printf("N: %d i: %d gc: %d\n", N, i, gc);*/
	      /* scale gc for chunks < 100 nt */
	      if(j < 100) gc *= 100. / (float) j;
	      /* don't count GC content of chunks < 20 nt, very hacky;
	       * don't count GC content of chunks that are all N, this
	       * will be common in RepeatMasked genomes where poly-Ns could
	       * skew the base composition stats of the genome */
	      if(j > 20 && !all_ambig_flag)
		{
		  /*printf("j: %d i: %d adding 1 to gc_ct[%d]\n", j, i, ((int) gc));*/
		  gc_ct[(int) gc] += 1.;
		}
	    }
	}
      esl_sq_Reuse(sq); 
    } 
  if (status != eslEOF) 
    cm_Fail("Parse failed, line %d, file %s:\n%s", 
	      sqfp->linenumber, sqfp->filename, sqfp->errbuf); 
  esl_sq_Destroy(sq); 
  esl_sqio_Rewind(sqfp);

#ifdef PRINT_GC_COUNTS
  for (i=0; i<GC_SEGMENTS; i++) 
    printf ("%d\t%.4f\n", i, gc_ct[i]);
#endif

  if(ret_N != NULL)      *ret_N     = N;
  if(ret_gc_ct != NULL)  *ret_gc_ct = gc_ct;
  else free(gc_ct);
  return; 

 ERROR:
  cm_Fail("Memory allocation error.");
}

/* Function: E2Score()
 * Date:     EPN, Tue Dec 11 15:40:25 2007 
 *           (morphed from RSEARCH: RJK, Mon April 15 2002 [St. Louis])
 *
 * Purpose:  Given an E-value <E> and a CM with valid Gumbel stats 
 *           determine the minimal bit score that will give an E-value 
 *           of <E> across all partitions. This will be a safe bit score
 *           cutoff to use when returning hits in DispatchSearch().
 *           Because no score less than this will have an E-value 
 *           less than E.
 *
 * Returns:  eslOK on success, <ret_sc> filled with bit score
 *           error code on failure, errbuf filled with message
 */
int E2Score (CM_t *cm, char *errbuf, int gum_mode, float E, float *ret_sc)
{
  float low_sc, sc;
  int p;

  /* contract checks */
  if(!(cm->flags & CMH_GUMBEL_STATS))        ESL_FAIL(eslEINCOMPAT, errbuf, "E2Score, CM's CMH_GUMBEL_STATS flag is down.");
  if(ret_sc == NULL)                         ESL_FAIL(eslEINCOMPAT, errbuf, "E2Score, ret_sc is NULL");

  low_sc = cm->stats->gumAA[gum_mode][0]->mu - (log(E) / cm->stats->gumAA[gum_mode][0]->lambda);
  if(! cm->stats->gumAA[gum_mode][0]->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "E2Score, CM's gumbel stats for mode: %d partition: %d are invalid.", gum_mode, p);
  for(p = 1; p < cm->stats->np; p++) {
    if(! cm->stats->gumAA[gum_mode][p]->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "E2Score, CM's gumbel stats for mode: %d partition: %d are invalid.", gum_mode, p);
    sc = cm->stats->gumAA[gum_mode][p]->mu - (log(E) / cm->stats->gumAA[gum_mode][p]->lambda);
    low_sc = ESL_MIN(low_sc, sc);
  }
  *ret_sc = low_sc;
  return eslOK;
}


/* Function: Score2E()
 * Date:     EPN, Wed Jan 16 14:25:44 2008
 *
 * Purpose:  Given a bit score <sc> and a CM with valid Gumbel stats
 *           determine the maximal E-value that will be assigned to 
 *           bit score of <sc> across all partitions. 
 *           This will be a 'safe' E-value cutoff to use 
 *           to always return hits with bit score of <sc> or greater.
 *           because no E-value higher than this will be assigned to
 *           a bit score greater than <sc>.
 *
 * Returns:  eslOK on success, <ret_E> filled with E value
 *           error code on failure, errbuf filled with message
 */
int Score2E (CM_t *cm, char *errbuf, int gum_mode, float sc, float *ret_E)
{
  float high_E, E;
  int p;

  /* contract checks */
  if(!(cm->flags & CMH_GUMBEL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), CM's CMH_GUMBEL_STATS flag is down.");
  if(ret_E == NULL)                   ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), ret_E is NULL");

  high_E = RJK_ExtremeValueE(sc, cm->stats->gumAA[gum_mode][0]->mu, cm->stats->gumAA[gum_mode][0]->lambda);

  if(! cm->stats->gumAA[gum_mode][0]->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), CM's gumbel stats for mode: %d partition: %d are invalid.", gum_mode, p);
  for(p = 1; p < cm->stats->np; p++) {
    if(! cm->stats->gumAA[gum_mode][p]->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), CM's gumbel stats for mode: %d partition: %d are invalid.", gum_mode, p);
    E = RJK_ExtremeValueE(sc, cm->stats->gumAA[gum_mode][p]->mu, cm->stats->gumAA[gum_mode][p]->lambda);
    high_E = ESL_MAX(high_E, E);
  }
  *ret_E = high_E;
  return eslOK;
}


/*
 * Function: RJK_ExtremeValueE
 * Date:     RJK, Mon Sep 30, 2002 [St. Louis]
 * Purpose:  Given a score (x), mu, and lambda, calculates 
 *           E=exp(-1*lambda(x-mu)) using first part of code from Sean's
 *           ExtremeValueP
 */
double RJK_ExtremeValueE (float x, double mu, double lambda) {
                        /* avoid underflow fp exceptions near P=0.0*/
  if ((lambda * (x - mu)) >= 2.3 * (double) DBL_MAX_10_EXP) 
    return 0.0;
  else 
    return(exp(-1. * lambda * (x - mu)));
}

/*
 * Function: MinScCutoff
 * Date:     EPN, Mon May  7 17:36:56 2007
 *
 * Purpose:  Return the minimum bit score cutoff for CM
 *           for round n in SearchInfo_t si.
 *           Trivial if si->cutoff_type[n] == SCORE_CUTOFF,
 *           but if E_CUTOFF return minimal bit score across 
 *           all partitions for the E cutoff in the 
 *           appropriate search algorithm.
 *
 */
float MinScCutoff (CM_t *cm, SearchInfo_t *si, int n)
{
  float E, low_sc, sc;
  int cm_mode, cp9_mode, gum_mode;
  int p; 

  /* contract check */
  if(si == NULL)      cm_Fail("MinCMScCutoff(), si == NULL.\n");
  if(n > si->nrounds) cm_Fail("MinCMScCutoff(), n (%d) > si->nrounds\n", n, si->nrounds);

  if(si->cutoff_type[n] == SCORE_CUTOFF) return si->sc_cutoff[n];
  
  /* if we get here, cutoff_type is E_CUTOFF we better have stats */
  ESL_DASSERT1((si->cutoff_type[n] == E_CUTOFF));
  if(!(cm->flags & CMH_GUMBEL_STATS)) cm_Fail("ERROR in MinScCutoff, cutoff type E value, but no stats.\n");

  /* Determine appropriate Gumbel mode */
  CM2Gumbel_mode(cm, si->search_opts[n], &cm_mode, &cp9_mode);
  E = si->e_cutoff[n];

  if(si->stype[n] == SEARCH_WITH_HMM) {
    ESL_DASSERT1(((si->search_opts[n] & CM_SEARCH_HMMVITERBI) || (si->search_opts[n] & CM_SEARCH_HMMFORWARD)));
    gum_mode = cp9_mode; 
  }
  else if (si->stype[n] == SEARCH_WITH_CM) {
    ESL_DASSERT1((! ((si->search_opts[n] & CM_SEARCH_HMMVITERBI) || (si->search_opts[n] & CM_SEARCH_HMMFORWARD))));
    gum_mode = cm_mode; 
  }
  else cm_Fail("MinScCutoff(), asking for E-value cutoff for SEARCH_WITH_HYBRID search round.\n");

  low_sc = cm->stats->gumAA[cm_mode][0]->mu - 
    (log(E) / cm->stats->gumAA[gum_mode][0]->lambda);
  for (p = 1; p < cm->stats->np; p++) {
    sc = cm->stats->gumAA[gum_mode][p]->mu - 
      (log(E) / cm->stats->gumAA[gum_mode][p]->lambda);
    if (sc < low_sc) low_sc = sc;
  }
  return (low_sc);
}


/*
 * Function: CM2Gumbel_mode
 * Date:     EPN, Mon May  7 17:43:28 2007
 * Purpose:  Return the gum_mode for the CM and HMM
 *           given the flags and search options in the
 *           CM data structure.
 */
int CM2Gumbel_mode(CM_t *cm, int search_opts, int *ret_cm_gum_mode, int *ret_cp9_gum_mode)
{
  int cm_gum_mode;
  int cp9_gum_mode;

  /* check contract */
  if(!(cm->flags & CMH_CP9) || cm->cp9 == NULL) cm_Fail("ERROR no CP9 in CM2Gumbel_mode()\n");
  if(cm->flags & CMH_LOCAL_BEGIN) {
    if(search_opts & CM_SEARCH_INSIDE) cm_gum_mode = GUM_CM_LI;
    else               	               cm_gum_mode = GUM_CM_LC;
  }
  else {
    if(search_opts & CM_SEARCH_INSIDE) cm_gum_mode = GUM_CM_GI;
    else        	               cm_gum_mode = GUM_CM_GC;
  }

  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) {
    if(search_opts & CM_SEARCH_HMMFORWARD) cp9_gum_mode = GUM_CP9_LF;
    else                                   cp9_gum_mode = GUM_CP9_LV;
  }
  else {
    if(search_opts & CM_SEARCH_HMMFORWARD) cp9_gum_mode = GUM_CP9_GF;
    else                                   cp9_gum_mode = GUM_CP9_GV;
  }

  if(ret_cm_gum_mode  != NULL) *ret_cm_gum_mode  = cm_gum_mode;
  if(ret_cp9_gum_mode != NULL) *ret_cp9_gum_mode = cp9_gum_mode;
  return eslOK;
}


/*
 * Function: CM2FthrMode
 * Date:     EPN, Tue Dec 11 13:16:35 2007
 * Purpose:  Return the filter threshold mode for the CM 
 *           given CM's flags and a passed in search options
 *           int.
 * 
 * Returns: eslOK on success, eslEINCOMPAT if search_opts indicate
 *          we're doing HMM search, errbuf is filled with error message.
 */
int CM2FthrMode(CM_t *cm, char *errbuf, int search_opts, int *ret_fthr_mode)
{
  int fthr_mode;

  /* check contract */
  if(search_opts & CM_SEARCH_HMMVITERBI) ESL_FAIL(eslEINCOMPAT, errbuf, "CM2FThrMode(), search_opts CM_SEARCH_HMMVITERBI flag raised.\n");
  if(search_opts & CM_SEARCH_HMMFORWARD) ESL_FAIL(eslEINCOMPAT, errbuf, "CM2FThrMode(), search_opts CM_SEARCH_HMMFORWARD flag raised.\n");
  if(ret_fthr_mode == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "CM2FThrMode(), ret_fthr_mode is NULL.");

  if(cm->flags & CMH_LOCAL_BEGIN) {
    if(search_opts & CM_SEARCH_INSIDE) fthr_mode = FTHR_CM_LI;
    else               	               fthr_mode = FTHR_CM_LC;
  }
  else {
    if(search_opts & CM_SEARCH_INSIDE) fthr_mode = FTHR_CM_GI;
    else        	               fthr_mode = FTHR_CM_GC;
  }
  if(ret_fthr_mode  != NULL) *ret_fthr_mode  = fthr_mode;
  return eslOK;
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
  cutoff = si->e_cutoff[si->nrounds];
  
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


/* Function: GumModeIsLocal()
 * Date:     EPN, Mon Dec 10 09:07:59 2007
 * Purpose:  Given a gumbel mode, return TRUE if it corresponds to 
 *           a local model configuration.
 *
 * Args:     gum_mode     - the mode 0..GUM_NMODES-1
 */
int
GumModeIsLocal(int gum_mode)
{
  ESL_DASSERT1((gum_mode >= 0 && gum_mode < GUM_NMODES));

  switch (gum_mode) {
  case GUM_CM_LC: 
  case GUM_CM_LI: 
  case GUM_CP9_LV: 
  case GUM_CP9_LF: 
    return TRUE;
    break;
  case GUM_CM_GC:
  case GUM_CM_GI: 
  case GUM_CP9_GV: 
  case GUM_CP9_GF: 
    return FALSE;
    break;
  default: 
    cm_Fail("GumModeIsLocal(): bogus gum_mode: %d\n", gum_mode);
  }
  return FALSE; /* never reached */
}


/* Function: GumModeIsForCM()
 * Date:     EPN, Mon Dec 10 09:11:55 2007
 * Purpose:  Given a gumbel mode, return TRUE if it corresponds to 
 *           a CM (not a CP9 HMM).
 *
 * Args:     gum_mode     - the mode 0..GUM_NMODES-1
 */
int
GumModeIsForCM(int gum_mode)
{
  ESL_DASSERT1((gum_mode >= 0 && gum_mode < GUM_NMODES));

  switch (gum_mode) {
  case GUM_CM_LC: 
  case GUM_CM_LI: 
  case GUM_CM_GC:
  case GUM_CM_GI: 
    return TRUE;
    break;
  case GUM_CP9_LV: 
  case GUM_CP9_LF: 
  case GUM_CP9_GV: 
  case GUM_CP9_GF: 
    return FALSE;
    break;
  default: 
    cm_Fail("GumModeIsForCM(): bogus gum_mode: %d\n", gum_mode);
  }
  return FALSE; /* never reached */
}


/* Function: GumModeToSearchOpts
 * Date:     EPN, Mon Dec 10 09:14:12 2007
 * Purpose:  Given a gumbel mode, update cm->search_opts to reflect
 *           that gumbel mode. 
 *
 *           0. GUM_CM_GC : !cm->search_opts & CM_SEARCH_INSIDE
 *           1. GUM_CM_GI : !cm->search_opts & CM_SEARCH_INSIDE
 *           4. GUM_CP9_GV:  cm->search_opts & CM_SEARCH_HMMVITERBI
 *                          !cm->search_opts & CM_SEARCH_HMMFORWARD
 *           5. GUM_CP9_GF:  cm->search_opts & CM_SEARCH_HMMVITERBI
 *                          !cm->search_opts & CM_SEARCH_HMMFORWARD
 *           3. GUM_CM_LC :  cm->search_opts & CM_SEARCH_INSIDE
 *           2. GUM_CM_LI :  cm->search_opts & CM_SEARCH_INSIDE
 *           6. GUM_CP9_LV: !cm->search_opts & CM_SEARCH_HMMVITERBI
 *                           cm->search_opts & CM_SEARCH_HMMFORWARD
 *           7. GUM_CP9_LF: !cm->search_opts & CM_SEARCH_HMMVITERBI
 *                           cm->search_opts & CM_SEARCH_HMMFORWARD
 * 
 * Args:
 *           CM           - the covariance model
 *           gum_mode     - the mode 0..GUM_NMODES-1
 */
int
GumModeToSearchOpts(CM_t *cm, int gum_mode)
{
  ESL_DASSERT1((gum_mode >= 0 && gum_mode < GUM_NMODES));

  switch (gum_mode) {
  case GUM_CM_GC: 
  case GUM_CM_LC: /* CYK, local or glocal */
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    break;
  case GUM_CM_GI: 
  case GUM_CM_LI: /* Inside, local or glocal */
    cm->search_opts |= CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    break;
  case GUM_CP9_GV: 
  case GUM_CP9_LV: /* Viterbi, local or glocal */
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts |= CM_SEARCH_HMMVITERBI;
    cm->search_opts &= ~CM_SEARCH_HMMFORWARD;
    break;
  case GUM_CP9_GF: 
  case GUM_CP9_LF: /* Forward, local or glocal */
    cm->search_opts &= ~CM_SEARCH_INSIDE;
    cm->search_opts &= ~CM_SEARCH_HMMVITERBI;
    cm->search_opts |= CM_SEARCH_HMMFORWARD;
    break;
  default: 
    cm_Fail("GumModeToSearchOpts(): bogus gum_mode: %d\n", gum_mode);
  }
  return eslOK;
}


/* Function: GumModeToFthrMode()
 * Date:     EPN, Mon Dec 10 09:31:48 2007
 * Purpose:  Given a gumbel mode, return it's corresponding
 *           filter threshold mode, or -1 if there is none
 *
 * Args:     gum_mode     - the mode 0..GUM_NMODES-1
 */
int
GumModeToFthrMode(int gum_mode)
{
  ESL_DASSERT1((gum_mode >= 0 && gum_mode < GUM_NMODES));

  switch (gum_mode) {
  case GUM_CM_GC: 
    return FTHR_CM_GC;
    break;
  case GUM_CM_GI: 
    return FTHR_CM_GI;
    break;
  case GUM_CM_LC: 
    return FTHR_CM_LC;
    break;
  case GUM_CM_LI: 
    return FTHR_CM_LI;
    break;
  case GUM_CP9_LV: 
  case GUM_CP9_LF: 
  case GUM_CP9_GV: 
  case GUM_CP9_GF: 
    return -1;
    break;
  default: 
    cm_Fail("GumModeToFthrMode(): bogus gum_mode: %d\n", gum_mode);
  }
  return FALSE; /* never reached */
}

/* Function: CreateGumbelInfo()
 * Date:     EPN, Tue Dec 11 05:25:06 2007
 *
 * Purpose:  Allocate and minimally initialize a gumbel info object.
 *            
 * Returns:  Newly allocated GumbelInfo_t object on success, NULL if some error occurs
 */
GumbelInfo_t *
CreateGumbelInfo()
{
  int status;

  GumbelInfo_t *gum = NULL;
  ESL_ALLOC(gum, sizeof(GumbelInfo_t));

  gum->N = 0;
  gum->L = 0;
  gum->mu = 0.;
  gum->lambda = 0.;
  gum->is_valid = FALSE;
  return gum;

 ERROR:
  return NULL; /* reached if memory error */
}  


/* Function: SetGumbelInfo()
 * Date:     EPN, Wed Dec 12 13:43:36 2007
 *
 * Purpose:  Set parameters of a gumbel info object and raise it's is_valid 'flag'.
 *            
 * Returns:  void
 */
void 
SetGumbelInfo(GumbelInfo_t *gum, double mu, double lambda, int L, int N)
{
  gum->N = N;
  gum->L = L;
  gum->mu = mu;
  gum->lambda = lambda;
  gum->is_valid = TRUE;
  return;
}  


/* Function: DuplicateGumbelInfo()
 * Date:     EPN, Tue Dec 11 05:28:13 2007
 *
 * Purpose:  Duplicate a gumbel info object.
 *            
 * Returns:  Newly allocated GumbelInfo_t object on success, NULL if some error occurs
 */
GumbelInfo_t *
DuplicateGumbelInfo(GumbelInfo_t *src)
{
  int status;

  GumbelInfo_t *dest = NULL;
  ESL_ALLOC(dest, sizeof(GumbelInfo_t));

  dest->N = src->N;
  dest->L = src->L;
  dest->mu = src->mu;
  dest->lambda = src->lambda;
  dest->is_valid = src->is_valid;
  return dest;

 ERROR:
  return NULL; /* reached if memory error */
}  


/* Function:  DescribeGumMode()
 * Incept:    EPN, Mon Jan  7 18:04:31 2008
 *
 * Purpose:   Returns the Gumbel mode in text.
 *            For example, <DescribeGumMode(GUM_CM_GC)>
 *            returns "glocal CM  CYK".
 */
char *
DescribeGumMode(int gum_mode)
{
  switch (gum_mode) {
  case GUM_CP9_GV: return "HMM  glc  vit";
  case GUM_CP9_GF: return "HMM  glc  fwd";
  case GUM_CM_GC:  return " CM  glc  cyk";
  case GUM_CM_GI:  return " CM  glc  ins";
  case GUM_CP9_LV: return "HMM  loc  vit";
  case GUM_CP9_LF: return "HMM  loc  fwd";
  case GUM_CM_LC:  return " CM  loc  cyk";
  case GUM_CM_LI:  return " CM  loc  ins";
  default:     return "?";
  }
}


/* Function:  DescribeFthrMode()
 * Incept:    EPN, Mon Jan  7 18:41:41 2008
 *
 * Purpose:   Returns the Filter thresold mode in text.
 *            For example, <DescribeFThrMode(GUM_CM_GC)>
 *            returns "glocal CM CYK".
 */
char *
DescribeFthrMode(int fthr_mode)
{
  switch (fthr_mode) {
  case FTHR_CM_GC:  return "glc  cyk";
  case FTHR_CM_GI:  return "glc  ins";
  case FTHR_CM_LC:  return "loc  cyk";
  case FTHR_CM_LI:  return "loc  ins";
  default:     return "?";
  }
}

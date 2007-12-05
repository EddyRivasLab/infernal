/*
 * EPN, Wed Nov 22 12:00:20 2006
 * Ported from RSEARCH-1.1
 *
 * stats.c
 *
 * Routines for calculating statistics in rsearch.  
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

/*
 * Function: AllocCMStats()
 * Date:     EPN, Wed May  2 14:10:25 2007      
 *
 * Purpose:  Allocate a CMStats_t data structure given
 *           the number of partitions.
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
  ESL_ALLOC(cmstats->gumAA, sizeof(struct gumbelinfo_s **) * NGUMBELMODES);
  ESL_ALLOC(cmstats->fthrA, sizeof(struct cp9filterthr_s *) * NFTHRMODES);
  for(i = 0; i < NGUMBELMODES; i++)
    {
      ESL_ALLOC(cmstats->gumAA[i], sizeof(struct gumbelinfo_s *));
      for(p = 0; p < cmstats->np; p++) { 
	ESL_ALLOC(cmstats->gumAA[i][p], sizeof(struct gumbelinfo_s));
	cmstats->gumAA[i][p]->isvalid = FALSE;
      }
    }
  for(i = 0; i < NFTHRMODES; i++) {
    ESL_ALLOC(cmstats->fthrA[i], sizeof(struct cp9filterthr_s));
    cmstats->fthrA[i]->isvalid = FALSE;
  }
  return cmstats;

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

/* Function: FreeCMStats()
 * Returns: (void) 
 */
void 
FreeCMStats(CMStats_t *cmstats)
{
  int i, p;
  for(i = 0; i < NGUMBELMODES; i++)
    {
      for(p = 0; p < cmstats->np; p++)
	free(cmstats->gumAA[i][p]);
      free(cmstats->gumAA[i]);
    }
  free(cmstats->gumAA);
  for(i = 0; i < NFTHRMODES; i++)
    free(cmstats->fthrA[i]);
  free(cmstats->fthrA);
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
      debug_print_gumbelinfo(cmstats->gumAA[CM_LC][p]);
      printf("cm_gc Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CM_GC][p]);
      printf("cm_li Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CM_LI][p]);
      printf("cm_gi Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CM_GI][p]);
      printf("cp9_lv Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CP9_LV][p]);
      printf("cp9_gv Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CP9_GV][p]);
      printf("cp9_lv Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CP9_LF][p]);
      printf("cp9_gv Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CP9_GF][p]);
      printf("\n\n");
    }

  if(has_fthr)
    {
      printf("fthr lc filter threshold:\n");
      debug_print_filterthrinfo(cmstats, cmstats->fthrA[CM_LC]);
      printf("fthr gc filter threshold:\n");
      debug_print_filterthrinfo(cmstats, cmstats->fthrA[CM_GC]);
      printf("fthr li filter threshold:\n");
      debug_print_filterthrinfo(cmstats, cmstats->fthrA[CM_LI]);
      printf("fthr gi filter threshold:\n");
      debug_print_filterthrinfo(cmstats, cmstats->fthrA[CM_GI]);
      printf("\n\n");
    }
  return eslOK;
}

/* Function: debug_print_gumbelinfo
 */
int debug_print_gumbelinfo(GumbelInfo_t *gum)
{
  if(gum->isvalid) printf("N: %d L: %d lambda: %.5f mu: %.5f\n", gum->N, gum->L, gum->lambda, gum->mu);
  else             printf("invalid (not yet set)\n");
  return eslOK;
}

/* Function: debug_print_filterthrinfo
 */
int debug_print_filterthrinfo(CMStats_t *cmstats, CP9FilterThr_t *fthr)
{
  if(! fthr->isvalid) { printf("invalid (not yet set)\n"); return eslOK; }

  double l_x;
  double g_x;
  double tmp_K, tmp_mu;
  tmp_K = exp(cmstats->gumAA[CP9_GF][0]->mu * cmstats->gumAA[CP9_GF][0]->lambda) / 
    cmstats->gumAA[CP9_GF][0]->L;
  tmp_mu = log(tmp_K * ((double) fthr->db_size)) / cmstats->gumAA[CP9_GF][0]->lambda;
  g_x = tmp_mu - (log(fthr->g_eval) / cmstats->gumAA[CP9_GF][0]->lambda);

  tmp_K = exp(cmstats->gumAA[CP9_LF][0]->mu * cmstats->gumAA[CP9_LF][0]->lambda) / 
    cmstats->gumAA[CP9_LF][0]->L;
  tmp_mu = log(tmp_K * ((double) fthr->db_size)) / cmstats->gumAA[CP9_LF][0]->lambda;
  l_x = tmp_mu - (log(fthr->l_eval) / cmstats->gumAA[CP9_LF][0]->lambda);
  printf("\tN: %d gsc: %.5f F: %.5f (%.5f bits) lsc: %.5f F: %.5f (%.5f bits)\n\tcmsc: %.5f db_size: %d was_fast: %d\n",
	 fthr->N, fthr->g_eval, fthr->g_F, g_x, fthr->l_eval, fthr->l_F, l_x, fthr->cm_eval, fthr->db_size, fthr->was_fast);
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
  int               i, j;  
  long              N = 0;
  double           *gc_ct;
  int               gc;
  char              allN[100]; /* used to check if curr DB chunk is all N's, if it is,
				* we don't count it towards the GC content info */
  int               allN_flag; /* stays up if curr DB chunk is all Ns */
  /*printf("in GetDBInfo\n");*/

  ESL_ALLOC(gc_ct, sizeof(double) * GC_SEGMENTS);
  for (i=0; i<GC_SEGMENTS; i++)
    gc_ct[i] = 0.;

  for (j=0; j<100; j++)
    allN[j] = 'N';
  
  if(ret_gc_ct != NULL) sq = esl_sq_CreateDigital(abc); 
  else                  sq = esl_sq_Create(); /* allows abc to be NULL if we don't care about gc stats */
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) 
    { 
      N += sq->n;
      /*printf("new N: %d\n", N);*/
      if(ret_gc_ct != NULL) 
	{
	  for(i = 0; i < sq->n; i += 100)
	    {
	      j = (i+99 <= sq->n) ? i+99 : sq->n;
	      gc = get_gc_comp(sq, i, j);
	      /*printf(">%d.raw\n", i);*/
	      allN_flag = TRUE;
	      for(j = 0; j < 100 && (j+i) < sq->n; j++)
		{
		  if(abc->sym[sq->dsq[(i+j)]] != 'N')
		    {
		      allN_flag = FALSE;
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
	      if(j > 20 && !allN_flag)
		{
		  /*printf("j: %d i: %d N: %d adding 1 to gc_ct[%d]\n", j, i, N, ((int) gc));*/
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

  if(ret_N != NULL)      *ret_N     = N;
  if(ret_gc_ct != NULL)  *ret_gc_ct = gc_ct;
  else free(gc_ct);
#ifdef PRINT_GC_COUNTS
  for (i=0; i<GC_SEGMENTS; i++) 
    printf ("%d\t%d\n", i, gc_ct[i]);
#endif
  
  return; 

 ERROR:
  cm_Fail("Memory allocation error.");
}

/*
 * Function: e_to_score()
 * Date:     RJK, Mon April 15 2002 [St. Louis]
 * Purpose:  Given an E-value and mu, lambda, and N, returns a value for S
 *           that will give such an E-value.  Basically the inverse of 
 *           ExtremeValueE
 */
float e_to_score (float E, double *mu, double *lambda) {
  float lowest_score, result;
  int i;

  lowest_score = mu[0] - (log(E)/lambda[0]);
  for (i=1; i<GC_SEGMENTS; i++) {
    result = mu[i] - (log(E)/lambda[i]);
    if (result < lowest_score)
      lowest_score = result;
  }
  return (lowest_score);
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

  if(si->cutoff_type[n] == SCORE_CUTOFF) return si->cutoff[n];
  
  /* if we get here, cutoff_type is E_CUTOFF we better have stats */
  ESL_DASSERT1((si->cutoff_type[n] == E_CUTOFF));
  if(!(cm->flags & CMH_GUMBEL_STATS)) cm_Fail("ERROR in MinScCutoff, cutoff type E value, but no stats.\n");

  /* Determine appropriate Gumbel mode */
  CM2Gumbel_mode(cm, si->search_opts[n], &cm_mode, &cp9_mode);
  E = si->cutoff[n];

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
    if(search_opts & CM_SEARCH_INSIDE) cm_gum_mode = CM_LI;
    else               	               cm_gum_mode = CM_LC;
  }
  else {
    if(search_opts & CM_SEARCH_INSIDE) cm_gum_mode = CM_GI;
    else        	               cm_gum_mode = CM_GC;
  }

  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) {
    if(search_opts & CM_SEARCH_HMMFORWARD) cp9_gum_mode = CP9_LF;
    else                                   cp9_gum_mode = CP9_LV;
  }
  else {
    if(search_opts & CM_SEARCH_HMMFORWARD) cp9_gum_mode = CP9_GF;
    else                                   cp9_gum_mode = CP9_GV;
  }

  if(ret_cm_gum_mode  != NULL) *ret_cm_gum_mode  = cm_gum_mode;
  if(ret_cp9_gum_mode != NULL) *ret_cp9_gum_mode = cp9_gum_mode;
  return eslOK;
}

/* Function: CopyFThrInfo()
 * Incept:   EPN, Fri May  4 15:54:51 2007
 */
int CopyFThrInfo(CP9FilterThr_t *src, CP9FilterThr_t *dest)
{
  dest->N           = src->N;
  dest->cm_eval     = src->cm_eval;
  dest->l_eval      = src->l_eval;
  dest->l_F         = src->l_F;
  dest->g_eval      = src->g_eval;
  dest->g_F         = src->g_F;
  dest->db_size     = src->db_size;
  dest->was_fast    = src->was_fast;
  return eslOK;
}

/* Function: CopyCMStatsGumbel()
 * Incept:   EPN, Mon May  7 06:04:58 2007
 * 
 * Purpose:  Copy the Gumbel stats in a source CMStats_t object into
 *           a pre-alloc'ed destination CMStats_t object.
 */
int CopyCMStatsGumbel(CMStats_t *src, CMStats_t *dest)
{
  int i, p;

  /* Check contract */
  if(src->np != dest->np)
    cm_Fail("ERROR in CopyCMStatsGumbel() src->np: %d not equal to alloc'ed dest->np: %d\n", src->np, dest->np);

  for(p = 0; p < src->np; p++)
    {
      dest->ps[p] = src->ps[p];
      dest->pe[p] = src->pe[p];
    }
  for(i = 0; i < GC_SEGMENTS; i++)
    dest->gc2p[i] = src->gc2p[i]; 

  for(i = 0; i < NGUMBELMODES; i++)
    {
      for(p = 0; p < src->np; p++)
	{
	  dest->gumAA[i][p]->N      = src->gumAA[i][p]->N;
	  dest->gumAA[i][p]->L      = src->gumAA[i][p]->L;
	  dest->gumAA[i][p]->mu     = src->gumAA[i][p]->mu;
	  dest->gumAA[i][p]->lambda = src->gumAA[i][p]->lambda;
	}
    }
  return eslOK;
}

/* Function: CopyCMStats()
 * Incept:   EPN, Tue May 29 06:00:41 2007
 * 
 * Purpose:  Copy the Gumbel and possibly CP9 Filter
 *           stats in a source CMStats_t object into
 *           a pre-alloc'ed destination CMStats_t object.
 */
int CopyCMStats(CMStats_t *src, CMStats_t *dest)
{
  /* Check contract */
  if(src->np != dest->np)
    cm_Fail("ERROR in CopyCMStats() src->np: %d not equal to alloc'ed dest->np: %d\n", src->np, dest->np);
  
  CopyCMStatsGumbel(src, dest);
  CopyFThrInfo(src->fthrA[CM_LC], dest->fthrA[CM_LC]);
  CopyFThrInfo(src->fthrA[CM_GC], dest->fthrA[CM_GC]);
  CopyFThrInfo(src->fthrA[CM_LI], dest->fthrA[CM_LI]);
  CopyFThrInfo(src->fthrA[CM_GI], dest->fthrA[CM_GI]);
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

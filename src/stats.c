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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <string.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stats.h"
#include "mpifuncs.h"
#include "histogram.h"
#include "easel.h"
#include "esl_histogram.h"
#include "esl_gumbel.h"
#include "esl_sqio.h"
#include "cm_dispatch.h"

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
  CMStats_t  *cmstats;
  int i, p;

  cmstats = (struct cmstats_s *) MallocOrDie(sizeof(struct cmstats_s));

  cmstats->np = np;
  cmstats->ps = MallocOrDie(sizeof(int) * cmstats->np);
  cmstats->pe = MallocOrDie(sizeof(int) * cmstats->np);
  cmstats->gumAA = MallocOrDie(sizeof(struct gumbelinfo_s **) * NGUMBELMODES);
  cmstats->fthrA = MallocOrDie(sizeof(struct cp9filterthr_s *) * NFTHRMODES);
  for(i = 0; i < NGUMBELMODES; i++)
    {
      cmstats->gumAA[i] = MallocOrDie(sizeof(struct gumbelinfo_s *));
      for(p = 0; p < cmstats->np; p++)
	cmstats->gumAA[i][p] = MallocOrDie(sizeof(struct gumbelinfo_s));
    }
  for(i = 0; i < NFTHRMODES; i++)
    cmstats->fthrA[i]  = MallocOrDie(sizeof(struct cp9filterthr_s));
  return cmstats;
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


/*
 * Function: SetCMCutoff
 * Date:     EPN, Thu May 17 13:30:41 2007
 * Purpose:  Fill cm->cutoff and cm->cutoff_type.
 */
int SetCMCutoff(CM_t *cm, int cm_cutoff_type, float cm_sc_cutoff, float cm_e_cutoff)
{
  if(cm->search_opts & CM_SEARCH_HMMONLY) /* CM score cutoff won't be used */
    {
      cm_cutoff_type = SCORE_CUTOFF;
      cm_sc_cutoff   = 0.;
    }
  else
    {
      cm->cutoff_type = cm_cutoff_type;
      if(cm->cutoff_type == SCORE_CUTOFF)
	cm->cutoff = cm_sc_cutoff;
      else 
	{
	  cm->cutoff = cm_e_cutoff;
	  if(!(cm->flags & CM_GUMBEL_STATS) && (!(cm->search_opts & CM_SEARCH_HMMONLY)))
	    Die("ERROR trying to use E-values but none in CM file.\nUse cmcalibrate or try -T.\n");
	}
    }
  return eslOK;
}

/*
 * Function: SetCP9Cutoff
 * Date:     EPN, Thu May 17 13:33:14 2007
 * Purpose:  Fill cm->cp9_cutoff and cm->cp9_cutoff_type.
 */
int SetCP9Cutoff(CM_t *cm, int cp9_cutoff_type, float cp9_sc_cutoff, float cp9_e_cutoff,
		 float cm_e_cutoff)
{
  if(cm->search_opts & CM_SEARCH_HMMONLY || 
     cm->search_opts & CM_SEARCH_HMMFILTER)
    {
      cm->cp9_cutoff_type = cp9_cutoff_type;  
      if(cm->cp9_cutoff_type == SCORE_CUTOFF)
	cm->cp9_cutoff = cp9_sc_cutoff;
      else 
	{

	  if(!(cm->flags & CM_GUMBEL_STATS))
	    Die("ERROR trying to use E-values but none in CM file.\nUse cmcalibrate or try --hmmT.\n");
	  if(cp9_e_cutoff < DEFAULT_MIN_CP9_E_CUTOFF) cp9_e_cutoff = DEFAULT_MIN_CP9_E_CUTOFF;
	  if(cm->cutoff_type == E_CUTOFF && cp9_e_cutoff < cm_e_cutoff) cp9_e_cutoff = cm_e_cutoff;
	  cm->cp9_cutoff = cp9_e_cutoff;
	}
    }
  else /* we won't use the CP9 at all, set score cutoff with 0 bit cutoff */
    {
      cm->cp9_cutoff_type = SCORE_CUTOFF;
      cm->cp9_cutoff      = 0.;
    }
  
  return eslOK;
}

/*
 * Function: PrintSearchInfo
 * Date:     EPN, Thu May 17 14:47:36 2007
 * Purpose:  Print info about search (cutoffs, algorithm, etc.) to file or stdout 
 */
int PrintSearchInfo(FILE *fp, CM_t *cm, int cm_mode, int cp9_mode, long N)
{
  int p, n;
  int clen = 0;
  float surv_fract;
  float avg_hit_len;

  for(n = 0; n < cm->nodes; n++) 
    if     (cm->ndtype[n] == MATP_nd) clen += 2;
    else if(cm->ndtype[n] == MATL_nd || cm->ndtype[n] == MATR_nd) clen += 1;

  if(!(cm->search_opts & CM_SEARCH_HMMONLY))
    {
      if(cm->cutoff_type == E_CUTOFF)
	{
	  fprintf(fp, "CM cutoff (E value):  %.2f\n", cm->cutoff);
	  for(p = 0; p < cm->stats->np; p++)
	    fprintf(fp, "   GC %2d-%3d bit sc:  %.2f mu: %.5f lambda: %.5f\n", cm->stats->ps[p], cm->stats->pe[p], 
		    (cm->stats->gumAA[cm_mode][p]->mu - 
		     (log(cm->cutoff) / cm->stats->gumAA[cm_mode][p]->lambda)), 
		    cm->stats->gumAA[cm_mode][p]->mu, cm->stats->gumAA[cm_mode][p]->lambda);
	}		       
      else if (cm->cutoff_type == SCORE_CUTOFF) 
	fprintf(fp, "CM cutoff (bit sc):   %.2f\n", cm->cutoff);
      printf ("CM search algorithm:  ");
      if(cm->search_opts & CM_SEARCH_INSIDE) fprintf(fp, "Inside\n");
      else fprintf(fp, "CYK\n");
      printf ("CM configuration:     ");
      if(cm->flags & CM_LOCAL_BEGIN) fprintf(fp, "Local\n");
      else fprintf(fp, "Glocal\n");
    }
  else 
    fprintf(fp, "Scanning with CP9 HMM only\n");
  if (cm->search_opts & CM_SEARCH_HMMFILTER)
    fprintf(fp, "Filtering with a CP9 HMM\n");
  
  if(cm->search_opts & CM_SEARCH_HMMONLY || 
     cm->search_opts & CM_SEARCH_HMMFILTER)
    {
      if(cm->cp9_cutoff_type == E_CUTOFF)
	{
	  if(!(cm->flags & CM_GUMBEL_STATS))
	    Die("ERROR trying to use E-values but none in CM file.\nUse cmcalibrate or try -T and/or --hmmT.\n");

	  /* Predict survival fraction from filter based on E-value, consensus length, W and N */
	  if(cp9_mode == CP9_G) avg_hit_len = clen;       /* should be weighted sum of gamma[0] from QDB calc */
	  if(cp9_mode == CP9_L) avg_hit_len = clen * 0.5; /* should be weighted sum of gamma[0] from QDB calc */
	  surv_fract = (cm->cp9_cutoff * ((2. * cm->W) - avg_hit_len)) / ((double) N); 
	  /* HMM filtering sends j-W..i+W to be researched with CM for HMM hits i..j */
	  fprintf(fp, "CP9 cutoff (E value): %.2f\n", cm->cp9_cutoff);
	  fprintf(fp, "   Predicted survival fraction: %.5f (1/%.3f)\n", surv_fract, (1./surv_fract));
	  for(p = 0; p < cm->stats->np; p++)
	    fprintf(fp, "   GC %2d-%3d bit sc:  %.2f mu: %.5f lambda: %.5f\n", cm->stats->ps[p], cm->stats->pe[p], 
		    (cm->stats->gumAA[cp9_mode][p]->mu - 
		     (log(cm->cp9_cutoff) / cm->stats->gumAA[cp9_mode][p]->lambda)), 
		    cm->stats->gumAA[cp9_mode][p]->mu, cm->stats->gumAA[cp9_mode][p]->lambda);
	}
      else if (cm->cp9_cutoff_type == SCORE_CUTOFF) 
	fprintf(fp, "CP9 cutoff (bit sc):  %.2f\n", cm->cp9_cutoff);
      printf ("CP9 search algorithm: Forward/Backward\n");
      printf ("CP9 configuration:    ");
      if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) fprintf(fp, "Local\n");
      else fprintf(fp, "Glocal\n");
    }
  printf     ("N (db size, nt):      %ld\n\n", N);
  fflush(stdout);
  return eslOK;
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
      printf("cp9_l Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CP9_L][p]);
      printf("cp9_g Gumbel:\t");
      debug_print_gumbelinfo(cmstats->gumAA[CP9_G][p]);
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
  printf("N: %d L: %d lambda: %.5f mu: %.5f\n", gum->N, gum->L, gum->lambda, gum->mu);
  return eslOK;
}

/* Function: debug_print_filterthrinfo
 */
int debug_print_filterthrinfo(CMStats_t *cmstats, CP9FilterThr_t *fthr)
{
  double l_x;
  double g_x;
  double tmp_K, tmp_mu;
  tmp_K = exp(cmstats->gumAA[CP9_G][0]->mu * cmstats->gumAA[CP9_G][0]->lambda) / 
    cmstats->gumAA[CP9_G][0]->L;
  tmp_mu = log(tmp_K * ((double) fthr->db_size)) / cmstats->gumAA[CP9_G][0]->lambda;
  g_x = tmp_mu - (log(fthr->g_eval) / cmstats->gumAA[CP9_G][0]->lambda);

  tmp_K = exp(cmstats->gumAA[CP9_L][0]->mu * cmstats->gumAA[CP9_L][0]->lambda) / 
    cmstats->gumAA[CP9_L][0]->L;
  tmp_mu = log(tmp_K * ((double) fthr->db_size)) / cmstats->gumAA[CP9_L][0]->lambda;
  l_x = tmp_mu - (log(fthr->l_eval) / cmstats->gumAA[CP9_L][0]->lambda);
  printf("\tN: %d gsc: %.5f (%.5f bits) lsc: %.5f (%.5f bits)\n\tcmsc: %.5f fraction: %.3f db_size: %d was_fast: %d\n",
	 fthr->N, fthr->g_eval, g_x, fthr->l_eval, l_x, fthr->cm_eval, fthr->fraction, fthr->db_size, fthr->was_fast);
  return eslOK;
}

/* Function: get_gc_comp
 * Date:     RJK, Mon Oct 7, 2002 [St. Louis]
 * Purpose:  Given a sequence and start and stop coordinates, returns 
 *           integer GC composition of the region 
 */
int get_gc_comp(char *seq, int start, int stop) {
  int i;
  int gc_ct;
  char c;

  if (start > stop) {
    i = start;
    start = stop;
    stop = i;
  }
  gc_ct = 0;
  /* Careful: seq is indexed 0..n-1 so start and
   * stop are off-by-one. This is a bug in RSEARCH-1.1 */
  for (i=(start-1); i<=(stop-1); i++) {
    c = resolve_degenerate(seq[i]);
    if (c=='G' || c == 'C')
      gc_ct++;
  }
  return ((int)(100.*gc_ct/(stop-start+1)));
}
 


/*
 * Function: random_from_string
 * Date:     September, 1998 (approx.) -- from hmmgcc
 * This function returns a character randomly chosen from the string.
 * Used in conjunction with the function below that resolves degenerate code
 * nucleotides.
 */
char random_from_string (char *s) {
  int i;
  do 
    {
      /*i = (int) ((float)(strlen(s)-1)*sre_random()/(RAND_MAX+1.0));*/
      i = (int) ((float)(strlen(s))*sre_random());
    } while (i<0 || i>=strlen(s));
  return(s[i]);
}

/*
 * Function: resolve_degenerate
 * Date:     September, 1998 (from hmmgcc)
 * This function resolves "degnerate" nucleotides by selecting a random 
 * A, C, G, or T as appropriate by the code present there.  Returns
 * the character passed in if that character does not represent a
 * non-degnerate nucleotide (either A, C, G, or T or not representative
 * at all of a nucleotide.
 *
 * The degenerate code used here is:
 * (taken from http://www.neb.com/neb/products/REs/RE_code.html
 *
 *                         R = G or A
 *                         K = G or T
 *                         B = not A (C or G or T)
 *                         V = not T (A or C or G)
 *                         Y = C or T
 *                         S = G or C
 *                         D = not C (A or G or T)
 *                         N = A or C or G or T
 *                         M = A or C
 *                         W = A or T
 *                         H = not G (A or C or T)
 *
 * This function assumes all letters are already uppercased via toupper
 * before calling.  In other words, it will return a "n" if passed an "n"
 * because it will assume that the symbol for all nucleotides will be passed
 * in as "N".
 */
char resolve_degenerate (char c) {
  c = toupper(c);
  switch (c) {
    case 'R' : return(random_from_string("GA"));
    case 'K' : return(random_from_string("GT"));
    case 'B' : return(random_from_string("CGT"));
    case 'V' : return(random_from_string("ACG"));
    case 'Y' : return(random_from_string("CT"));
    case 'S' : return(random_from_string("GC"));
    case 'D' : return(random_from_string("AGT"));
    case 'N' : return(random_from_string("ACGT"));
    case 'M' : return(random_from_string("AC"));
    case 'W' : return(random_from_string("AT"));
    case 'H' : return(random_from_string("ACT"));
  }
  return(c);
}

/*
 * Function: GetDBInfo()
 * Date:     Easelification: EPN, Thu Dec  7 06:07:58 2006
 *           (initial - RSEARCH::get_dbinfo()) RJK, Thu Apr 11, 2002 [St. Louis]
 * Purpose:  Given a sequence file name, determine the total size of the
 *           seqs in the file (DB) and GC content information.
 * Args:     seqfile  - name of sequence file
 *           ret_N    - RETURN: total length (residues) or all seqs in seqfile
 *           gc_ct    - RETURN: gc_ct[x] observed 100-nt segments with GC% of x [0..100] 
 * seqfile   
 */
void GetDBInfo (ESL_SQFILE *sqfp, long *ret_N, double **ret_gc_ct) 
{
  ESL_SQ           *sq;
  int               i, j;  
  long              N = 0;
  double           *gc_ct = MallocOrDie(sizeof(double) * GC_SEGMENTS);
  int               status;
  int               gc;
  char              c;
  char              allN[100]; /* used to check if curr DB chunk is all N's, if it is,
				* we don't count it towards the GC content info */
  int               allN_flag; /* stays up if curr DB chunk is all Ns */
  /*printf("in GetDBInfo\n");*/

  for (i=0; i<GC_SEGMENTS; i++)
    gc_ct[i] = 0.;

  for (j=0; j<100; j++)
    allN[j] = 'N';
  
  sq = esl_sq_Create(); 
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) 
    { 
      N += sq->n;
      /*printf("new N: %d\n", N);*/
      for(i = 0; i < sq->n; i += 100)
	{
	  gc = 0;
	  /*printf(">%d.raw\n", i);*/
	  allN_flag = TRUE;
	  for(j = 0; j < 100 && (j+i) < sq->n; j++)
	    {
	      if(allN_flag && sq->seq[(j+i)] != 'N')
		allN_flag = FALSE;
	      /*printf("%c", sq->seq[(j+i)]);*/
	      c = resolve_degenerate(sq->seq[(j+i)]);
	      if (c == 'G' || c == 'C') gc++;
	    }
	  /*printf("\n>%d.resolved\n", i);*/
	  for(j = 0; j < 100 && (j+i) < sq->n; j++)
	    {
	      c = resolve_degenerate(sq->seq[(j+i)]);
	      /*printf("%c", c);*/
	    }
	  /*printf("N: %d i: %d gc: %d\n", N, i, gc);*/
	  /* scale gc for chunks < 100 nt */
	  if(j < 100)
	    {
	      gc *= 100. / (float) j;
	    }
	  /*if(allN_flag)
	    printf("allN_flag UP!\n");*/
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
      esl_sq_Reuse(sq); 
    } 
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %d, file %s:\n%s", 
	      sqfp->linenumber, sqfp->filename, sqfp->errbuf); 
  esl_sq_Destroy(sq); 
  esl_sqio_Rewind(sqfp);

  if(ret_N != NULL)      *ret_N     = N;
  if(ret_gc_ct != NULL)  *ret_gc_ct = gc_ct;
#ifdef PRINT_GC_COUNTS
  for (i=0; i<GC_SEGMENTS; i++) 
    printf ("%d\t%d\n", i, gc_ct[i]);
#endif
  
  return; 
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
 * Function: MinCMScCutoff
 * Date:     EPN, Mon May  7 17:36:56 2007
 * Purpose:  Return the minimum bit score cutoff for a CM.
 *           Trivial if cm->cutoff_type == SCORE_CUTOFF,
 *           if E_CUTOFF return minimal bit score across 
 *           all partitions for the E cutoff in the 
 *           appropriate search algorithm (local/glocal
 *           CYK/Inside combo)
 */
float MinCMScCutoff (CM_t *cm)
{
  float E, low_sc, sc;
  int gum_mode;
  int p; 

  if(cm->cutoff_type == SCORE_CUTOFF)
    return cm->cutoff;
  
  /* we better have stats */
  if(!(cm->flags & CM_GUMBEL_STATS))
    Die("ERROR in MinCMScCutoff, cutoff type E value, but no stats.\n");

  /* Determine appropriate Gumbel mode */
  CM2Gumbel_mode(cm, &gum_mode, 
	      NULL); /* don't care about CP9 Gumbel mode */
  /*printf("in MinCMScCutoff, gum_mode: %d\n", gum_mode);*/
  E = cm->cutoff;

  low_sc = cm->stats->gumAA[gum_mode][0]->mu - 
    (log(E) / cm->stats->gumAA[gum_mode][0]->lambda);
  for (p = 1; p < cm->stats->np; p++) 
    {
      sc = cm->stats->gumAA[gum_mode][p]->mu - 
	(log(E) / cm->stats->gumAA[gum_mode][p]->lambda);
      if (sc < low_sc)
	low_sc = sc;
  }
  return (low_sc);
}

/*
 * Function: MinCP9ScCutoff
 * Date:     EPN, Mon May  7 17:36:56 2007
 * Purpose:  Return the minimum bit score cutoff for a CM's
 *           CP9 HMM. Trivial if cm->cp9_cutoff_type == SCORE_CUTOFF,
 *           if E_CUTOFF return minimal bit score across 
 *           all partitions for the E cutoff in the 
 *           appropriate search algorithm (local/glocal)
 */
float MinCP9ScCutoff (CM_t *cm)
{
  float E, low_sc, sc;
  int gum_mode;
  int p;

  if(cm->cp9_cutoff_type == SCORE_CUTOFF)
    return cm->cp9_cutoff;
  
  /* we better have stats */
  if(!(cm->flags & CM_GUMBEL_STATS))
    Die("ERROR in MinCMScCutoff, cutoff type E value, but no stats.\n");

  /* Determine appropriate Gumbel mode */
  CM2Gumbel_mode(cm, NULL,  /* don't care about CM Gumbel mode */
	      &gum_mode);
  E = cm->cp9_cutoff;

  low_sc = cm->stats->gumAA[gum_mode][0]->mu - 
    (log(E) / cm->stats->gumAA[gum_mode][0]->lambda);
  for (p=1; p < cm->stats->np; p++) 
    {
      sc = cm->stats->gumAA[gum_mode][p]->mu - 
	(log(E) / cm->stats->gumAA[gum_mode][p]->lambda);
      if (sc < low_sc)
	low_sc = sc;
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
int CM2Gumbel_mode(CM_t *cm, int *ret_cm_gum_mode, 
		   int *ret_cp9_gum_mode)
{
  int cm_gum_mode;
  int cp9_gum_mode;

  /* check contract */
  if(!(cm->flags & CM_CP9) || cm->cp9 == NULL)
    Die("ERROR no CP9 in CM2Gumbel_mode()\n");

  if(cm->flags & CM_LOCAL_BEGIN)
    {
      if(cm->search_opts & CM_SEARCH_INSIDE)
	cm_gum_mode = CM_LI;
      else
	cm_gum_mode = CM_LC;
    }
  else
    {
      if(cm->search_opts & CM_SEARCH_INSIDE)
	cm_gum_mode = CM_GI;
      else
	cm_gum_mode = CM_GC;
    }

  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN)
    cp9_gum_mode = CP9_L;
  else
    cp9_gum_mode = CP9_G;

  if(ret_cm_gum_mode  != NULL) *ret_cm_gum_mode  = cm_gum_mode;
  if(ret_cp9_gum_mode != NULL) *ret_cp9_gum_mode = cp9_gum_mode;
  return eslOK;
}



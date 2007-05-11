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
  cmstats->evdAA = MallocOrDie(sizeof(struct evdinfo_s **) * NEVDMODES);
  cmstats->fthrA = MallocOrDie(sizeof(struct cp9filterthr_s *) * NFTHRMODES);
  for(i = 0; i < NEVDMODES; i++)
    {
      cmstats->evdAA[i] = MallocOrDie(sizeof(struct evdinfo_s *));
      for(p = 0; p < cmstats->np; p++)
	cmstats->evdAA[i][p] = MallocOrDie(sizeof(struct evdinfo_s));
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
  for(i = 0; i < NEVDMODES; i++)
    {
      for(p = 0; p < cmstats->np; p++)
	free(cmstats->evdAA[i][p]);
      free(cmstats->evdAA[i]);
    }
  free(cmstats->evdAA);
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
      printf("cm_lc EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_LC][p]);
      printf("cm_gc EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_GC][p]);
      printf("cm_li EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_LI][p]);
      printf("cm_gi EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CM_GI][p]);
      printf("cp9_l EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CP9_L][p]);
      printf("cp9_g EVD:\t");
      debug_print_evdinfo(cmstats->evdAA[CP9_G][p]);
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

/* Function: debug_print_evdinfo
 */
int debug_print_evdinfo(EVDInfo_t *evd)
{
  printf("N: %d L: %d lambda: %.5f mu: %.5f\n", evd->N, evd->L, evd->lambda, evd->mu);
  return eslOK;
}

/* Function: debug_print_filterthrinfo
 */
int debug_print_filterthrinfo(CMStats_t *cmstats, CP9FilterThr_t *fthr)
{
  double l_x;
  double g_x;
  
  g_x = cmstats->evdAA[CP9_G][0]->mu - 
    (log(fthr->g_eval) / cmstats->evdAA[CP9_G][0]->lambda);
  l_x = cmstats->evdAA[CP9_L][0]->mu - 
    (log(fthr->l_eval) / cmstats->evdAA[CP9_L][0]->lambda);
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
void GetDBInfo (ESL_SQFILE *sqfp, long *ret_N, int **ret_gc_ct) 
{
  ESL_SQ           *sq;
  int               i, j;  
  long              N = 0;
  int              *gc_ct = MallocOrDie(sizeof(int) * GC_SEGMENTS);
  int               status;
  int               gc;
  char              c;
  char              allN[100]; /* used to check if curr DB chunk is all N's, if it is,
				* we don't count it towards the GC content info */
  int               allN_flag; /* stays up if curr DB chunk is all Ns */
  /*printf("in GetDBInfo\n");*/

  for (i=0; i<GC_SEGMENTS; i++)
    gc_ct[i] = 0;

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
	      gc_ct[(int) gc]++;
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
  int evd_mode;
  int p; 

  if(cm->cutoff_type == SCORE_CUTOFF)
    return cm->cutoff;
  
  /* we better have stats */
  if(!(cm->flags & CM_EVD_STATS))
    Die("ERROR in MinCMScCutoff, cutoff type E value, but no stats.\n");

  /* Determine appropriate EVD mode */
  CM2EVD_mode(cm, &evd_mode, 
	      NULL); /* don't care about CP9 EVD mode */
  E = cm->cutoff;

  low_sc = cm->stats->evdAA[evd_mode][0]->mu - 
    (log(E) / cm->stats->evdAA[evd_mode][0]->lambda);
  for (p = 1; p < cm->stats->np; p++) 
    {
      sc = cm->stats->evdAA[evd_mode][p]->mu - 
	(log(E) / cm->stats->evdAA[evd_mode][p]->lambda);
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
  int evd_mode;
  int p;

  if(cm->cp9_cutoff_type == SCORE_CUTOFF)
    return cm->cp9_cutoff;
  
  /* we better have stats */
  if(!(cm->flags & CM_EVD_STATS))
    Die("ERROR in MinCMScCutoff, cutoff type E value, but no stats.\n");

  /* Determine appropriate EVD mode */
  CM2EVD_mode(cm, NULL,  /* don't care about CM EVD mode */
	      &evd_mode);
  E = cm->cp9_cutoff;

  low_sc = cm->stats->evdAA[evd_mode][0]->mu - 
    (log(E) / cm->stats->evdAA[evd_mode][0]->lambda);
  for (p=1; p < cm->stats->np; p++) 
    {
      sc = cm->stats->evdAA[evd_mode][p]->mu - 
	(log(E) / cm->stats->evdAA[evd_mode][p]->lambda);
      if (sc < low_sc)
	low_sc = sc;
  }
  return (low_sc);
}


/*
 * Function: CM2EVD_mode
 * Date:     EPN, Mon May  7 17:43:28 2007
 * Purpose:  Return the EVD_mode for the CM and HMM
 *           given the flags and search options in the
 *           CM data structure.
 */
int CM2EVD_mode(CM_t *cm, int *ret_cm_evd_mode, 
		int *ret_cp9_evd_mode)
{
  int cm_evd_mode;
  int cp9_evd_mode;

  /* check contract */
  if(!(cm->flags & CM_CP9) || cm->cp9 == NULL)
    Die("ERROR no CP9 in CM2EVD_mode()\n");

  if(cm->flags & CM_LOCAL_BEGIN)
    {
      if(cm->search_opts & CM_SEARCH_INSIDE)
	cm_evd_mode = CM_LI;
      else
	cm_evd_mode = CM_LC;
    }
  else
    {
      if(cm->search_opts & CM_SEARCH_INSIDE)
	cm_evd_mode = CM_GI;
      else
	cm_evd_mode = CM_GC;
    }

  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN)
    cp9_evd_mode = CP9_L;
  else
    cp9_evd_mode = CP9_G;

  if(ret_cm_evd_mode  != NULL) *ret_cm_evd_mode  = cm_evd_mode;
  if(ret_cp9_evd_mode != NULL) *ret_cp9_evd_mode = cp9_evd_mode;
  return eslOK;
}


#if 0
/*
 * Function: serial_make_histogram()
 * Date:     Mon Apr 1 2002 [St. Louis]
 * Purpose:  Makes a histogram using random sequences.  Returns mu and lambda.
 *           Makes random sequences of length dblen, finds best hit at 
 *           arbitrary j's every D nucleotides along database.
 *           Also returns K (how much to scale N in calculating E-value)
 *
 * Inputs:   gc_comp     %GC of random seq
 *           cm          the model
 *           num_samples number of samples to take
 *           sample_length  length of each sample
 *           use_easel - TRUE to use Easel's histogram and EVD fitters, false to
 *                       use RSEARCH versions.
 *
 */  
void serial_make_histogram (int *gc_count, int *partitions, int num_partitions,
			    CM_t *cm, int num_samples, int sample_length, 
			    int doing_cp9_stats, int use_easel)
{
  int i;
  char *randseq;
  char *dsq;
  /*struct histogram_s *h_old;*/
  ESL_HISTOGRAM *h;
  float *nt_p;                /* Distribution for random sequences */
  float score;
  int cur_partition;
  float cur_gc_freq[GC_SEGMENTS];
  float gc_comp;
  double curr_lambda;         /* lambda for current partition */
  double curr_mu;             /* mu for current partition */
  double *xv;
  int z;
  int n;

  /*printf("in serial_make_histogram, nparts: %d sample_len: %d cp9_stats: %d do_ins: %d do_enf: %d\n", num_partitions, sample_length, doing_cp9_stats, (cm->search_opts & CM_SEARCH_INSIDE), (cm->config_opts & CM_CONFIG_ENFORCE));*/

  /* Allocate for random distribution */
  nt_p = MallocOrDie(sizeof(float)*Alphabet_size); 

  /* For each partition */
  for (cur_partition = 0; cur_partition < num_partitions; cur_partition++) 
    {
      /*printf("cur_partition: %d\n", cur_partition);*/
      
      /* Initialize histogram; these numbers are guesses */
      /*if(!use_easel) h_old = AllocHistogram (0, 100, 100);*/
      h     = esl_histogram_CreateFull(0., 100., 1.);    
      
      /* Set up cur_gc_freq */
      for (i=0; i<GC_SEGMENTS; i++) 
	{
	  if (partitions[i] == cur_partition) 
	    {
	      cur_gc_freq[i] = (float)gc_count[i];
	      /*printf("cur_gc_freq(cur_partition:%d)[i:%d]: %f\n", cur_partition, i, cur_gc_freq[i]);*/
	    } 
	  else
	    cur_gc_freq[i] = 0.;
	}

      FNorm(cur_gc_freq, GC_SEGMENTS);

      /* Take num_samples samples */
      for (i=0; i<num_samples; i++) 
	{
	  /* Get random GC content */
	  gc_comp = 0.01*FChoose (cur_gc_freq, GC_SEGMENTS);
	  /*printf("SH GC: %f part: %d randseq%d\n", gc_comp, cur_partition, i);*/

	  nt_p[1] = nt_p[2] = 0.5*gc_comp;
	  nt_p[0] = nt_p[3] = 0.5*(1. - gc_comp);
	  
	  /* Get random sequence */
	  randseq = RandomSequence (Alphabet, nt_p, Alphabet_size, sample_length);

	  /* Digitize the sequence, parse it, and add to histogram */
	  dsq = DigitizeSequence (randseq, sample_length);
	  
	  /* Do the scan */
	  score = 
	    actually_search_target(cm, dsq, 1, sample_length, 
				   0.,    /* cutoff is 0 bits (actually we'll find highest
					   * negative score if it's < 0.0) */
				   0.,    /* CP9 cutoff is 0 bits */
				   NULL,  /* don't keep results */
				   FALSE, /* don't filter with a CP9 HMM */
				   (!doing_cp9_stats), /* TRUE if we're calc'ing CM stats */
				   doing_cp9_stats,    /* TRUE if we're calc'ing CP9 stats */
				   NULL);          /* filter fraction N/A */
	  /*if(i % 100 == 0)
	    printf("(%4d) SCORE: %f\n", i, score);*/
	  /* Add best score to histogram */
	  /*if(!use_easel) AddToHistogram (h_old, score); */
	  esl_histogram_Add(h, score);
	  /*printf("\n\nSH RANDOMSEQ BIT SC: %f\n", score);
	    printf(">randseq\n");
	    printf("%s\n", randseq);*/

	}

      /* Fit the histogram.  */
      /*if(!use_easel)
	ExtremeValueFitHistogram (h_old, TRUE, 9999);*/

      /* Fit the scores to a Gumbel */
      
      /* If the esl_histogram example for 'complete data, high scores fit as
       *  censored Gumbel' is the correct approach we do this: */
      esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
      
      /* If the esl_histogram example for 'censored data, fit as censored Gumbel
       *  is the correct approach, we would do this (BUT IT'S NOT!): */
      /*esl_histogram_GetData(h, &xv, &n);*/
      
      esl_gumbel_FitCensored(xv, n, z, xv[0], &curr_mu, &curr_lambda);
      
    for (i=0; i<GC_SEGMENTS; i++) 
      {
	if (partitions[i] == cur_partition) 
	{
	  if(use_easel)
	    {
	      if(doing_cp9_stats)
		{
		  cm->cp9_lambda[i] = curr_lambda;
		  cm->cp9_K[i] = exp(curr_mu * curr_lambda)/sample_length;
		}
	      else /* we're calc'ing stats for the CM */
		{
		  cm->lambda[i] = curr_lambda;
		  cm->K[i] = exp(curr_mu * curr_lambda)/sample_length;
		}
	      /*printf("ESL i: %d lambda: %f K: %f\n", i, lambda[i], K[i]);*/
	      /*printf("OLD i: %d lambda: %f K: %f\n\n", i, ((double) h_old->param[EVD_LAMBDA]),
		((double) exp(h_old->param[EVD_MU]*h_old->param[EVD_LAMBDA])/sample_length));*/
	    }
	  else /* use RSEARCH's histogram code */
	    {
	      ;/*lambda[i] = (double) h_old->param[EVD_LAMBDA];
		K[i] = (double) exp(h_old->param[EVD_MU]*h_old->param[EVD_LAMBDA])/sample_length;*/
	      /*printf("OLD i: %d lambda: %f K: %f\n", i, lambda[i], K[i]);
		printf("ESL i: %d lambda: %f K: %f\n\n", i, curr_lambda, exp(curr_mu * curr_lambda)/sample_length);*/
	    }
	}
      }
    /*if(!use_easel) FreeHistogram(h_old);*/
    esl_histogram_Destroy(h);
  }
  free(nt_p);
}

#ifdef USE_MPI
void parallel_make_histogram (int *gc_count, int *partitions, int num_partitions, 
			      CM_t *cm, int num_samples, int sample_length,
			      int doing_cp9_stats,
			      int mpi_my_rank, int mpi_num_procs, 
			      int mpi_master_rank) 
{
  /*struct histogram_s **h;*/
  ESL_HISTOGRAM **h;
  int z;
  int n;
  double *xv;
  double curr_lambda;         /* lambda for current partition */
  double curr_mu;             /* mu for current partition */

  db_seq_t **randseqs;          /* The random sequences */
  int randseq_index;
  int num_seqs_made;
  float *nt_p;                  /* Distribution for random sequences */
  int proc_index;
  job_t **process_status;
  job_t *job_queue = NULL;
  char *dsq;
  int seqlen;
  char job_type;
  float score;
  int dummy;              /* To hold bestr, which isn't used in these jobs */
  int i = 0;
  int cur_partition;
  float **cur_gc_freqs;
  float gc_comp;
  char *tmp_name;
  float *enf_vec;             /* vector for FChoose to pick starting point for enf_seq */
  int enf_start;           /* starting point for enf_seq */
  int x;

  /* Infernal specific variables (not in RSEARCH's stats.c) */
  int    nhits;			/* number of hits in a seq */
  int   *hitr;			/* initial states for hits */
  int   *hiti;                  /* start positions of hits */
  int   *hitj;                  /* end positions of hits */
  float *hitsc;			/* scores of hits */

  /*printf("in parallel_make_histogram, nparts: %d sample_len: %d cp9_stats: %d do_ins: %d do_enf: %d\n", num_partitions, sample_length, doing_cp9_stats, (cm->search_opts & CM_SEARCH_INSIDE), (cm->config_opts & CM_CONFIG_ENFORCE));*/

  tmp_name = sre_strdup("random", -1);

  if (mpi_my_rank == mpi_master_rank) 
    {
      /* Allocate random distribution */
      nt_p = MallocOrDie(sizeof(float)*Alphabet_size); 
      
      /* Allocate histograms and set up cur_gc_freq */
      h = MallocOrDie(sizeof (ESL_HISTOGRAM *)*num_partitions);
      cur_gc_freqs = MallocOrDie(sizeof(float *)*num_partitions);
      for (cur_partition = 0; cur_partition<num_partitions; cur_partition++) {
	/* Initialize histogram; these numbers are guesses */
	h[cur_partition] = esl_histogram_CreateFull(0., 100., 1.);    
	/*h[cur_partition] = AllocHistogram (0, 100, 100);*/
	
	/* Set up cur_gc_freq */
	cur_gc_freqs[cur_partition] = MallocOrDie(sizeof(float)*GC_SEGMENTS);
	for (i=0; i<GC_SEGMENTS; i++) {
	  if (partitions[i] == cur_partition) {
	    cur_gc_freqs[cur_partition][i] = (float)gc_count[i];
	    /*printf("cur_gc_freqs[cur_partition:%d][i:%d]: %f\n", cur_partition, i, cur_gc_freqs[cur_partition][i]);*/
	  } else {
	    cur_gc_freqs[cur_partition][i] = 0.;
	  }
	}
	FNorm (cur_gc_freqs[cur_partition], GC_SEGMENTS);
	/*printf("\n\n");*/
      }
      /* Set up arrays to hold pointers to active seqs and jobs on
	 processes */
      randseqs = MallocOrDie(sizeof(db_seq_t *) * mpi_num_procs);
      process_status = MallocOrDie(sizeof(job_t *)*mpi_num_procs);
      for (randseq_index=0; randseq_index<mpi_num_procs; randseq_index++)
	randseqs[randseq_index] = NULL;
      for (proc_index = 0; proc_index < mpi_num_procs; proc_index++)
	process_status[proc_index] = NULL;
      
      num_seqs_made = 0;
      cur_partition = 0;
      
      do {
	/* Check for idle processes.  Send jobs */
	for (proc_index=0; proc_index<mpi_num_procs; proc_index++) 
	  {
	    if (proc_index == mpi_master_rank) continue; /* Skip master process */
	    if (process_status[proc_index] == NULL) 
	      {         
		/* I'm idle -- need a job */
		if (job_queue == NULL) {
		  for (randseq_index = 0; randseq_index <mpi_num_procs; 
		       randseq_index++) {
		    if (randseqs[randseq_index] == NULL) break;
		  }
		  if (randseq_index == mpi_num_procs)
		    Die ("Tried to read more than %d seqs at once\n", mpi_num_procs);
		  /* Get random sequence, digitize it */
		  if (num_seqs_made == num_samples * num_partitions)
		    /* We've got all we need.  Stop */
		    break;
		  
		  cur_partition = num_seqs_made / num_samples;
		  randseqs[randseq_index] = MallocOrDie(sizeof(db_seq_t));
		  
		  /* Get random GC content */
		  gc_comp = 0.01*FChoose(cur_gc_freqs[cur_partition], GC_SEGMENTS);
		  nt_p[1] = nt_p[2] = 0.5*gc_comp;
		  nt_p[0] = nt_p[3] = 0.5*(1. - gc_comp);
		  /*printf("PH GC: %f part: %d %s\n", gc_comp, cur_partition, tmp_name);*/
		  randseqs[randseq_index]->sq[0] = 
		    esl_sq_CreateFrom(tmp_name, 
				      RandomSequence (Alphabet, nt_p, Alphabet_size, sample_length),
				      NULL, NULL, NULL);
		  
		  randseqs[randseq_index]->sq[0]->dsq = 
		    DigitizeSequence(randseqs[randseq_index]->sq[0]->seq, 
				     randseqs[randseq_index]->sq[0]->n);

		  randseqs[randseq_index]->best_score = IMPOSSIBLE;
		  
		  randseqs[randseq_index]->partition = cur_partition;
		  
		  job_queue = search_enqueue(randseqs[randseq_index], randseq_index, cm->W,
					     FALSE, SEARCH_HIST_SCAN_WORK);
		  num_seqs_made++;
		  
		}
		if (job_queue != NULL)
		  {
		    fflush(stdout);
		    search_send_next_job (&job_queue, process_status + proc_index, 
					  proc_index);
		  }
	      }
	  }
	/* Wait for next reply */
	if (search_procs_working(process_status, mpi_num_procs, mpi_master_rank)) {
	  randseq_index = search_check_hist_results (randseqs, process_status, cm->W);
	  /* If the sequence is done */
	  if (randseqs[randseq_index]->chunks_sent == 0) {
	    /* Get best score at cm->W and add */
	    /*AddToHistogram (h[randseqs[randseq_index]->partition], 
	      randseqs[randseq_index]->best_score);*/
	    esl_histogram_Add(h[randseqs[randseq_index]->partition], randseqs[randseq_index]->best_score);

	    /*printf("\n\nPH RANDOMSEQ BIT SC: %f\n", randseqs[randseq_index]->best_score);
	      printf(">randseq\n");
	      printf("%s\n", randseqs[randseq_index]->sq[0]->seq);*/
	    esl_sq_Destroy(randseqs[randseq_index]->sq[0]);
	    free(randseqs[randseq_index]);
	    randseqs[randseq_index] = NULL;
	  }
	}
      } while (num_seqs_made < num_samples*num_partitions || job_queue != NULL ||
	       search_procs_working(process_status, mpi_num_procs, mpi_master_rank));
      
      /* Terminate the processes */
      for (proc_index=0; proc_index<mpi_num_procs; proc_index++) {
	if (proc_index != mpi_master_rank) {
	  search_send_terminate (proc_index);
	}
      }
      
      /* Fit the histogram.  */
      for (cur_partition=0; cur_partition<num_partitions; cur_partition++) {
	/*ExtremeValueFitHistogram (h[cur_partition], TRUE, 9999);*/
	
	esl_histogram_GetTailByMass(h[cur_partition], 0.5, &xv, &n, &z); /* fit to right 50% */
	esl_gumbel_FitCensored(xv, n, z, xv[0], &curr_mu, &curr_lambda);
	
	for (i=0; i<GC_SEGMENTS; i++) 
	  {
	    if (partitions[i] == cur_partition) 
	      {
		if(doing_cp9_stats) 
		  {
		    cm->cp9_lambda[i] = curr_lambda;
		    cm->cp9_K[i] = exp(curr_mu * curr_lambda)/sample_length;
		    /*printf("P ESL i: %d lambda: %f K: %f\n", i, lambda[i], K[i]);*/
		  }
		else /* we're calcing stats for the CM */
		  {
		    cm->lambda[i] = curr_lambda;
		    cm->K[i] = exp(curr_mu * curr_lambda)/sample_length;
		    /*printf("P ESL i: %d lambda: %f K: %f\n", i, lambda[i], K[i]);*/
		  }
	      }
	  }
      }
      free(randseqs);
      free(process_status);
      free(nt_p);
      for (i=0; i<cur_partition; i++) {
	esl_histogram_Destroy(h[i]);
	free(cur_gc_freqs[i]);
      }
      free(cur_gc_freqs);
      free(h);
    }
  else 
    {
      dsq = NULL;
      do 
	{
	  job_type = search_receive_job(&seqlen, &dsq, &dummy, mpi_master_rank);
	  if (job_type == SEARCH_HIST_SCAN_WORK) 
	    {
	      score = 
		actually_search_target(cm, dsq, 1, sample_length, 
				       0.,    /* CM cutoff is 0 bits (actually we'll find highest
					       * negative score if it's < 0.0) */
				       0.,    /* CP9 cutoff is 0 bits */
				       NULL,  /* don't keep results */
				       FALSE, /* don't filter with a CP9 HMM */
				       (!doing_cp9_stats), /* TRUE if we're calc'ing CM stats */
				       doing_cp9_stats,    /* TRUE if we're calc'ing CP9 stats */
				       NULL);          /* filter fraction N/A */
	      /*printf("score: %f\n", score);*/
	      search_send_hist_scan_results (score, mpi_master_rank);
	    }
	  if (dsq != NULL)
	    free(dsq);
	  dsq = NULL;
	} while (job_type != TERMINATE_WORK);
    }
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif
#endif

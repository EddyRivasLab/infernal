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
#include "esl_histogram.h"
#include "esl_gumbel.h"
/*#include "mpifuncs.h"*/

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
 *           D           maximum value for D
 *           num_samples number of samples to take
 *           sample_length  length of each sample
 *           mu_p        pointer to final mu
 *           lambda_p    pointer to final lambda
 *           dmin   - minimum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *           dmax   - maximum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *           hmm    - CP 9 HMM to scan with, NULL if we're scanning with the CM
 */  
void serial_make_histogram (int *gc_count, int *partitions, int num_partitions,
			    CM_t *cm, int D, int num_samples,
			    int sample_length, double *lambda, double *K,
			    int *dmin, int *dmax, struct cplan9_s *hmm)
{
  int i;
  char *randseq;
  char *dsq;
  /*struct histogram_s *h;*/
  ESL_HISTOGRAM *h;
  double *scores;             /* [0..num_samples-1] best score from each sample */
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

  /* Infernal specific variables (not in RSEARCH's stats.c) */
  int    nhits;			/* number of hits in a seq */
  int   *hitr;			/* initial states for hits */
  int   *hiti;                  /* start positions of hits */
  int   *hitj;                  /* end positions of hits */
  float *hitsc;			/* scores of hits */

  printf("in serial_make_histogram, nparts: %d D: %d sample_len: %d\n", num_partitions, D, sample_length);

  scores = MallocOrDie(sizeof (double) * num_samples);

  if (num_samples == 0) {
    for (i=0; i<GC_SEGMENTS; i++) {
      K[i] = 0.;
      lambda[i] = 0.;
    }
    return;        /* No sequence to analyse */
  }

  /* Allocate for random distribution */
  nt_p = MallocOrDie(sizeof(float)*Alphabet_size); 

  /* For each partition */
  for (cur_partition = 0; cur_partition < num_partitions; cur_partition++) {
    printf("cur_partition: %d\n", cur_partition);
   
    /* Initialize histogram; these numbers are guesses */
    /*h = AllocHistogram (0, 100, 100);*/
    h = esl_histogram_CreateFull(0., 100., 1.);    

    /* Set up cur_gc_freq */
    for (i=0; i<GC_SEGMENTS; i++) {
      if (partitions[i] == cur_partition) {
	cur_gc_freq[i] = (float)gc_count[i];
	/*printf("cur_gc_freq[%d] set to %f\n", i, cur_gc_freq[i]);*/
      } else {
	cur_gc_freq[i] = 0.;
      }
    }
    FNorm(cur_gc_freq, GC_SEGMENTS);

    /* Take num_samples samples */
    for (i=0; i<num_samples; i++) {

      /* Get random GC content */
      gc_comp = 0.01*FChoose (cur_gc_freq, GC_SEGMENTS);
      nt_p[1] = nt_p[2] = 0.5*gc_comp;
      nt_p[0] = nt_p[3] = 0.5*(1. - gc_comp);

      /* Get random sequence */
      randseq = RandomSequence (Alphabet, nt_p, Alphabet_size, sample_length);
    
      /* Digitize the sequence, parse it, and add to histogram */
      dsq = DigitizeSequence (randseq, sample_length);

      /* Do the scan */
      if(hmm != NULL)
	score = CP9ForwardScan(dsq, 1, sample_length, D, hmm, NULL,  /* don't save the DP mx */
			       &nhits, &hitr, &hiti, &hitj, &hitsc, 
			       0);
      else if(dmin == NULL && dmax == NULL)
	score = CYKScan(cm, dsq, 1, sample_length, D,
			&nhits, &hitr, &hiti, &hitj, &hitsc, 0);
      else /* use QDB */
	score = CYKBandedScan(cm, dsq, dmin, dmax, 1, sample_length, D,
			      &nhits, &hitr, &hiti, &hitj, &hitsc, 0);

      if(i % 100 == 0)
	printf("(%4d) SCORE: %f\n", i, score);
      /* Add best score to histogram */
      /* AddToHistogram (h, score); */
      esl_histogram_Add(h, score);

      /* Free stuff */
      free (randseq);
      free (dsq);
      free(hitr);
      free(hiti);
      free(hitj);
      free(hitsc);
    }

    /*esl_gumbel_FitComplete(scores, num_samples, &curr_mu, &curr_lambda);*/
    /*esl_gumbel_FitCensored(scores, num_samples, &curr_mu, &curr_lambda);*/
    /* Fit the scores to a Gumbel */

    /* If the esl_histogram example for 'complete data, high scores fit as
     *  censored Gumbel' is the correct approach we do this: */
    esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */

    /* If the esl_histogram example for 'censored data, fit as censored Gumbel
     *  is the correct approach, we would do this (BUT IT'S NOT!): */
      /*esl_histogram_GetData(h, &xv, &n);*/
    
    esl_gumbel_FitCensored(xv, n, z, xv[0], &curr_mu, &curr_lambda);
    
    for (i=0; i<GC_SEGMENTS; i++) {
      if (partitions[i] == cur_partition) {
	lambda[i] = curr_lambda;
	K[i] = exp(curr_mu * curr_lambda)/sample_length;
	printf("i: %d lambda: %f K: %f\n", i, lambda[i], K[i]);
      }
    }
    esl_histogram_Destroy(h);
  }
  free(nt_p);
}

#ifdef USE_MPI
void parallel_make_histogram (int *gc_count, int *partitions, int num_partitions, 
			      CM_t *cm, int D, int num_samples,
			      int sample_length,float *lambda, float *K, 
			      int mpi_my_rank, int mpi_num_procs, 
			      int mpi_master_rank) {
  struct histogram_s **h;
  db_seq_t **randseqs;          /* The random sequences */
  int randseq_index;
  int num_seqs_made;
  float *nt_p;                  /* Distribution for random sequences */
  int proc_index;
  job_t **process_status;
  job_t *job_queue = NULL;
  char *seq;
  int seqlen;
  char job_type;
  float score;
  int dummy;              /* To hold bestr, which isn't used in these jobs */
  int i = 0;
  int cur_partition;
  float **cur_gc_freqs;
  float gc_comp;

  if (num_samples == 0) {
    for (i = 0; i<GC_SEGMENTS; i++) {
      K[i] = 0.;
      lambda[i] = 0.;
    }
  }

  if (mpi_my_rank == mpi_master_rank) {

    /* Allocate random distribution */
    nt_p = MallocOrDie(sizeof(float)*Alphabet_size); 

    /* Allocate histograms and set up cur_gc_freq */
    h = MallocOrDie(sizeof (struct histogram_s *)*num_partitions);
    cur_gc_freqs = MallocOrDie(sizeof(float *)*num_partitions);
    for (cur_partition = 0; cur_partition<num_partitions; cur_partition++) {
      /* Initialize histogram; these numbers are guesses */
      h[cur_partition] = AllocHistogram (0, 100, 100);
 
      /* Set up cur_gc_freq */
      cur_gc_freqs[cur_partition] = MallocOrDie(sizeof(float)*GC_SEGMENTS);
      for (i=0; i<GC_SEGMENTS; i++) {
	if (partitions[i] == cur_partition) {
	  cur_gc_freqs[cur_partition][i] = (float)gc_count[i];
	} else {
	  cur_gc_freqs[cur_partition][i] = 0.;
	}
      }
      FNorm (cur_gc_freqs[cur_partition], GC_SEGMENTS);
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
      for (proc_index=0; proc_index<mpi_num_procs; proc_index++) {
	if (proc_index == mpi_master_rank) continue; /* Skip master process */
	if (process_status[proc_index] == NULL) {         
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
	    randseqs[randseq_index]->sqinfo.len = sample_length;
	    /* Get random GC content */
	    gc_comp = 0.01*FChoose(cur_gc_freqs[cur_partition], GC_SEGMENTS);
	    nt_p[1] = nt_p[2] = 0.5*gc_comp;
	    nt_p[0] = nt_p[3] = 0.5*(1. - gc_comp);
	    randseqs[randseq_index]->seq[0] = 
	      RandomSequence (Alphabet, nt_p, Alphabet_size, 
			      randseqs[randseq_index]->sqinfo.len);
	    randseqs[randseq_index]->dsq[0] = 
	      DigitizeSequence (randseqs[randseq_index]->seq[0], 
				randseqs[randseq_index]->sqinfo.len);
	    randseqs[randseq_index]->best_score = 0;
	    randseqs[randseq_index]->partition = cur_partition;
	    job_queue = enqueue (randseqs[randseq_index], randseq_index, D,
				 FALSE, HIST_SCAN_WORK);
	    num_seqs_made++;
	  }
	  if (job_queue != NULL)
	    send_next_job (&job_queue, process_status + proc_index, 
			   proc_index);
	} 
      } 
      /* Wait for next reply */
      if (procs_working(process_status, mpi_num_procs, mpi_master_rank)) {
	randseq_index = check_hist_results (randseqs, process_status, D);
	/* If the sequence is done */
	if (randseqs[randseq_index]->chunks_sent == 0) {
	  /* Get best score at D and add */
	  AddToHistogram (h[randseqs[randseq_index]->partition], 
			  randseqs[randseq_index]->best_score);
	  free(randseqs[randseq_index]->dsq[0]);
	  free(randseqs[randseq_index]->seq[0]);
	  free(randseqs[randseq_index]);
	  randseqs[randseq_index] = NULL;
	}
      }
    } while (num_seqs_made < num_samples*num_partitions || job_queue != NULL ||
	     procs_working(process_status, mpi_num_procs, mpi_master_rank));
      
    /* Terminate the processes */
    for (proc_index=0; proc_index<mpi_num_procs; proc_index++) {
      if (proc_index != mpi_master_rank) {
	send_terminate (proc_index);
      }
    }
    
    /* Fit the histogram.  */
    for (cur_partition=0; cur_partition<num_partitions; cur_partition++) {
      ExtremeValueFitHistogram (h[cur_partition], TRUE, 9999);
      for (i=0; i<GC_SEGMENTS; i++) {
	if (partitions[i] == cur_partition) {
	  lambda[i] = h[cur_partition]->param[EVD_LAMBDA];
	  K[i] = exp(h[cur_partition]->param[EVD_MU]*h[cur_partition]->param[EVD_LAMBDA])/sample_length;
	}
      }
    }
    free(randseqs);
    free(process_status);
    free(nt_p);
    for (i=0; i<cur_partition; i++) {
      FreeHistogram(h[i]);
      free(cur_gc_freqs[i]);
    }
    free(cur_gc_freqs);
    free(h);
  } else {
    seq = NULL;
    do {
      job_type = recieve_job (&seqlen, &seq, &dummy, mpi_master_rank);
      if (job_type == HIST_SCAN_WORK) {
	score = CYKScan (cm, seq, seqlen, IMPOSSIBLE, D, NULL);
	send_hist_scan_results (score, mpi_master_rank);
      } 
      if (seq != NULL)
	free(seq);
      seq = NULL;
    } while (job_type != TERMINATE_WORK);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

/*
 * Function: random_from_string
 * Date:     September, 1998 (approx.) -- from hmmgcc
 * This function returns a character randomly chosen from the string.
 * Used in conjunction with the function below that resolves degenerate code
 * nucleotides.
 */
char random_from_string (char *s) {
  int i;
  do {
    i = (int) ((float)(strlen(s)-1)*sre_random()/(RAND_MAX+1.0));
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
 * Function: get_dbinfo()
 * Date:     RJK, Thu Apr 11, 2002 [St. Louis]
 * Purpose:  Given a pointer to an open but as yet unread sequence database 
 *           file, place to put the size, partition info, and place to put
 *           GC content, calculates total size of database and GC info 
 */
void get_dbinfo (SQFILE *dbfp, long *N, int *gc_count) {
  long counter = 0;
  char *seq;
  char c;
  int i,j;
  int gc;
  SQINFO sqinfo;

  for (i=0; i<GC_SEGMENTS; i++)
    gc_count[i] = 0;

  while (ReadSeq (dbfp, dbfp->format, &seq, &sqinfo)) {
    counter += sqinfo.len;
    for (i = 0; i<sqinfo.len; i+=100) {
      gc = 0;
      for (j=0; j<100; j++) {
	c = resolve_degenerate(seq[i+j]);
	if (c == 'G' || c == 'C')
	  gc++;
      }
      gc_count[gc]++;
    }
    FreeSequence (seq, &sqinfo);
  }
#ifdef PRINT_GC_COUNTS
  for (i=0; i<GC_SEGMENTS; i++) 
    printf ("%d\t%d\n", i, gc_count[i]);
#endif
  *N = counter;
  SeqfileRewind(dbfp);
}

/*
 * Function: e_to_score()
 * Date:     RJK, Mon April 15 2002 [St. Louis]
 * Purpose:  Given an E-value and mu, lambda, and N, returns a value for S
 *           that will give such an E-value.  Basically the inverse of 
 *           ExtremeValueE
 */
float e_to_score (float E, float *mu, float *lambda) {
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
  if ((lambda * (x - mu)) >= 2.3 * (double) DBL_MAX_10_EXP) {
    return 0.0;
  } else {
    return(exp(-1. * lambda * (x - mu)));
  }
}

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

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "easel.h"
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "sre_stack.h"
#include "hmmband.h"         
#include "cm_postprob.h"
#include "stats.h"
#include "esl_gumbel.h"
#include "mpifuncs.h"
#include "cm_dispatch.h"

/* Helper functions called by the main functions (main functions
 * declared in cm_dispatch.h) 
 */
static db_seq_t *read_next_seq (ESL_SQFILE *dbfp, int do_revcomp);
static void print_results (CM_t *cm, CMConsensus_t *cons, db_seq_t *dbseq,
			   int do_complement, int do_stats, double *mu, 
			   double *lambda) ;
static int get_gc_comp(char *seq, int start, int stop);
static void remove_hits_over_e_cutoff (CM_t *cm, scan_results_t *results, char *seq,
				       float cutoff, double *lambda, double *mu);
static seqs_to_aln_t *read_next_aln_seqs(ESL_SQFILE *seqfp, int nseq, int index);

/* coordinate -- macro that checks if it's reverse complement and if so 
   returns coordinate in original strand
   a = true if revcomp, false if not
   b = the position in current seq
   c = length of the seq
*/
#define coordinate(a,b,c) ( a ? -1*b+c+1 : b)

/*
 * Function: serial_search_database
 * Date:     EPN, Thu Dec  7 18:07:24 2006
 *           original code by RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Given an open database file, a model, and various parameters, does
 *           the search using CYKScan and then determines and prints out the 
 *           alignments.
 *
 * Parameters:        dbfp         the database
 *                    cm           the model
 *                    cons         precalc'ed consensus info for display
 */
void serial_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons)
{
  int reversed;                /* Am I currently doing reverse complement? */
  int i,a;
  db_seq_t *dbseq;
  float min_cm_cutoff;
  float min_cp9_cutoff;
  int do_revcomp;
  int do_align;
  /*printf("in serial_search database do_align: %d do_revcomp: %d\n", do_align, do_revcomp);*/
  
  /* Determine minimum cutoff for CM and for CP9 */
  if (cm->cutoff_type == SCORE_CUTOFF) 
    min_cm_cutoff = cm->cutoff;
  else 
    min_cm_cutoff = e_to_score (cm->cutoff, cm->mu, cm->lambda);

  if (cm->cp9_cutoff_type == SCORE_CUTOFF) 
    min_cp9_cutoff = cm->cp9_cutoff;
  else 
    min_cp9_cutoff = e_to_score (cm->cp9_cutoff, cm->cp9_mu, cm->cp9_lambda);
  
  do_revcomp = (!(cm->search_opts & CM_SEARCH_TOPONLY));
  do_align   = (!(cm->search_opts & CM_SEARCH_NOALIGN));

  while ((dbseq = read_next_seq(dbfp, do_revcomp)))
    {
      for (reversed = 0; reversed <= do_revcomp; reversed++) 
	{
	  /*printf("SEARCHING >%s %d\n", dbseq->sq[reversed]->name, reversed);*/
	  /* Scan */
	  dbseq->results[reversed] = CreateResults(INIT_RESULTS);
	  actually_search_target(cm, dbseq->sq[reversed]->dsq, 1, dbseq->sq[reversed]->n, 
				 min_cm_cutoff, min_cp9_cutoff, 
				 dbseq->results[reversed], 
				 TRUE,  /* filter with a CP9 HMM if appropriate */
				 FALSE, /* we're not building a histogram for CM stats  */
				 FALSE, /* we're not building a histogram for CP9 stats */
				 NULL); /* filter fraction, TEMPORARILY NULL            */
	  if(cm->search_opts & CM_SEARCH_GREEDY)
	    {
	      remove_overlapping_hits (dbseq->results[reversed],
				       dbseq->sq[reversed]->n);
	    }
	  if (cm->cutoff_type == E_CUTOFF) 
	    remove_hits_over_e_cutoff (cm, dbseq->results[reversed],
				       dbseq->sq[reversed]->seq,
				       cm->cutoff, cm->lambda, cm->mu);

	  /* Align results */
	  if (do_align) 
	    {
	    for (i=0; i<dbseq->results[reversed]->num_results; i++) 
	      {
		CYKDivideAndConquer
		  (cm, dbseq->sq[reversed]->dsq, dbseq->sq[reversed]->n,
		   dbseq->results[reversed]->data[i].bestr,
		   dbseq->results[reversed]->data[i].start, 
		   dbseq->results[reversed]->data[i].stop, 
		   &(dbseq->results[reversed]->data[i].tr),
		   cm->dmin, cm->dmax); /* dmin and dmax will be NULL if non-banded 
					 * alternatively, could always pass NULL to 
					 * always do non-banded alignment. */
		
		/* Now, subtract out the starting point of the result so 
		   that it can be added in later.  This makes the print_alignment
		   routine compatible with the parallel version, while not needing
		   to send the entire database seq over for each alignment. */
		for (a=0; a<dbseq->results[reversed]->data[i].tr->n; a++) 
		  {
		    dbseq->results[reversed]->data[i].tr->emitl[a] -= 
		      (dbseq->results[reversed]->data[i].start - 1);
		    dbseq->results[reversed]->data[i].tr->emitr[a] -= 
		      (dbseq->results[reversed]->data[i].start - 1);
		  }
	      }
	    }
	}
      /* Print results */
      if(cm->search_opts & CM_SEARCH_HMMONLY) /* pass CP9 EVD params */
	print_results (cm, cons, dbseq, do_revcomp, (cm->search_opts & CM_SEARCH_CP9STATS),
		       cm->cp9_mu, cm->cp9_lambda);
      else /* pass CM EVD params */
	print_results (cm, cons, dbseq, do_revcomp, (cm->search_opts & CM_SEARCH_CMSTATS),
		       cm->mu, cm->lambda);

      fflush (stdout);
      
      FreeResults(dbseq->results[0]);
      esl_sq_Destroy(dbseq->sq[0]);
      if (do_revcomp) 
	{
	  esl_sq_Destroy(dbseq->sq[1]);
	  FreeResults(dbseq->results[1]);
	}
      free(dbseq);
    }
}

#ifdef USE_MPI
/*
 * Function: parallel_search_database 
 * Date:     original (RSEARCH): RJK, Tue May 28, 2002 [St. Louis] 
 *           infernalized: EPN, Sat Dec  9 12:57:42 2006
 * Purpose:  Given the same parameters as serial_search_database, does
 *           the database search with alignments and printing results, but
 *           in a parallel fashion.
 *     
 *           It the master node performs the following tasks:
 *           while (!eof(dbfp) && slaves_working) {
 *             For each empty process, send next job
 *             Wait for a result
 *             For each empty process, send next job
 *             Process result
 *           }
 *            
 *
 *           The slave processes do the following:
 *           while (1) {
 *             receive_job
 *             if (terminate code) return;
 *             do the job
 *             send the results
 *           }
 */
void parallel_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons,
			       int mpi_my_rank, int mpi_master_rank, 
			       int mpi_num_procs) 
{
  char job_type;
  int seqlen;
  char *seq;
  scan_results_t *results;
  db_seq_t **active_seqs;
  job_t **process_status;
  int eof = FALSE;
  job_t *job_queue = NULL;
  int proc_index, active_seq_index;
  int bestr;       /* Best root state -- for alignments */
  Parsetree_t *tr;
  float min_cm_cutoff;
  float min_cp9_cutoff;
  int do_revcomp;
  int do_align;


  do_revcomp = (!(cm->search_opts & CM_SEARCH_TOPONLY));
  do_align   = (!(cm->search_opts & CM_SEARCH_NOALIGN));

  /* Determine minimum cutoff for CM and for CP9 */
  if (cm->cutoff_type == SCORE_CUTOFF) 
    min_cm_cutoff = cm->cutoff;
  else 
    min_cm_cutoff = e_to_score (cm->cutoff, cm->mu, cm->lambda);

  if (cm->cp9_cutoff_type == SCORE_CUTOFF) 
    min_cp9_cutoff = cm->cp9_cutoff;
  else 
    min_cp9_cutoff = e_to_score (cm->cp9_cutoff, cm->cp9_mu, cm->cp9_lambda);

  /*printf("B PSD rank: %4d mast: %4d cm_cutoff: %f cp9_cutoff: %f\n", mpi_my_rank, mpi_master_rank, min_cm_cutoff, min_cp9_cutoff);*/

  if (mpi_my_rank == mpi_master_rank) 
    {
      /* Set up arrays to hold pointers to active seqs and jobs on
	 processes */
      active_seqs = MallocOrDie(sizeof(db_seq_t *) * mpi_num_procs);
      process_status = MallocOrDie(sizeof(job_t *) * mpi_num_procs);
      for (active_seq_index=0; active_seq_index<mpi_num_procs; active_seq_index++) 
	active_seqs[active_seq_index] = NULL;
      for (proc_index = 0; proc_index < mpi_num_procs; proc_index++)
	process_status[proc_index] = NULL;
      
      do
	{
	  /* Check for idle processes.  Send jobs */
	  for (proc_index=0; proc_index<mpi_num_procs; proc_index++) 
	    {
	      if (proc_index == mpi_master_rank) continue;  /* Skip master process */
	      if (process_status[proc_index] == NULL) 
		{         
		  /* I'm idle -- need a job */
		  if (job_queue == NULL) 
		    { 
		      /* Queue is empty */
		      /* Find next non-master open process */
		      for (active_seq_index=0; active_seq_index<mpi_num_procs; active_seq_index++) 
			if (active_seqs[active_seq_index] == NULL) break;
		      if (active_seq_index == mpi_num_procs) 
			Die ("Tried to read more than %d seqs at once\n", mpi_num_procs);
		      active_seqs[active_seq_index] = read_next_seq(dbfp, do_revcomp);
		      if (active_seqs[active_seq_index] == NULL) 
			{
			  eof = TRUE;
			  break;            /* Queue is empty and no more seqs */
			}
		      else
			{
			  job_queue = search_enqueue (active_seqs[active_seq_index], 
						      active_seq_index, cm->W, do_revcomp, 
						      SEARCH_STD_SCAN_WORK);
			}
		    }
		  if (job_queue != NULL)
		    search_send_next_job (&job_queue, process_status + proc_index, proc_index);
		} 
	    }
	  /* Wait for next reply */
	  if (search_procs_working(process_status, mpi_num_procs, mpi_master_rank)) 
	    {
	      active_seq_index = search_check_results (active_seqs, process_status, cm->W);
	      if (active_seqs[active_seq_index]->chunks_sent == 0) 
		{
		  remove_overlapping_hits
		    (active_seqs[active_seq_index]->results[0], 
		     active_seqs[active_seq_index]->sq[0]->n);
		  if (cm->cutoff_type == E_CUTOFF)
		    remove_hits_over_e_cutoff 
		      (cm, active_seqs[active_seq_index]->results[0],
		       active_seqs[active_seq_index]->sq[0]->seq,
		       cm->cutoff, cm->lambda, cm->mu);
		  if (do_revcomp) 
		    {
		      remove_overlapping_hits 
			(active_seqs[active_seq_index]->results[1],
			 active_seqs[active_seq_index]->sq[1]->n);
		      if (cm->cutoff_type == E_CUTOFF)
			remove_hits_over_e_cutoff 
			  (cm, active_seqs[active_seq_index]->results[1],
			   active_seqs[active_seq_index]->sq[1]->seq,
			   cm->cutoff, cm->lambda, cm->mu);
		    }
		  /* Check here if doing alignments and queue them or check if 
		     all done */
		  if (do_align && 
		      active_seqs[active_seq_index]->alignments_sent == -1) 
		    {
		      search_enqueue_alignments (&job_queue, active_seqs[active_seq_index],
						 active_seq_index, do_revcomp, ALN_WORK);
		    }
		  if (!do_align || 
		      active_seqs[active_seq_index]->alignments_sent == 0) 
		    {
		      /* Print results */
		      if(cm->search_opts & CM_SEARCH_HMMONLY) /* pass CP9 EVD params */
			print_results (cm, cons, active_seqs[active_seq_index], 
				       do_revcomp, (cm->search_opts & CM_SEARCH_CP9STATS), 
				       cm->cp9_mu, cm->cp9_lambda);
		      else /* pass CM EVD params */
			print_results (cm, cons, active_seqs[active_seq_index], 
				       do_revcomp, (cm->search_opts & CM_SEARCH_CMSTATS), cm->mu, cm->lambda);
		      
		    if (do_revcomp) 
		      {
			FreeResults(active_seqs[active_seq_index]->results[1]);
			esl_sq_Destroy(active_seqs[active_seq_index]->sq[1]);
		      }
		    FreeResults(active_seqs[active_seq_index]->results[0]);
		    esl_sq_Destroy(active_seqs[active_seq_index]->sq[0]);
		    active_seqs[active_seq_index] = NULL;
		  }
		}
	    }
	} while (!eof || job_queue != NULL || 
		 (search_procs_working(process_status, mpi_num_procs, mpi_master_rank)));
      for (proc_index=0; proc_index<mpi_num_procs; proc_index++) 
	if (proc_index != mpi_master_rank) 
	  search_send_terminate (proc_index);

      free(active_seqs);
      free(process_status);
    } 
  else  /* not the master node */
    {
      seq = NULL;
      do 
	{
	  job_type = search_receive_job (&seqlen, &seq, &bestr, mpi_master_rank);
	  if (job_type == SEARCH_STD_SCAN_WORK) 
	    {
	      /* Do the scan */
	      results = CreateResults(INIT_RESULTS);
	      actually_search_target(cm, seq, 1, seqlen, min_cm_cutoff, min_cp9_cutoff,
				     results, /* keep results */
				     TRUE,  /* filter with a CP9 HMM if appropriate */
				     FALSE, /* we're not building a histogram for CM stats  */
				     FALSE, /* we're not building a histogram for CP9 stats */
				     NULL); /* filter fraction, TEMPORARILY NULL            */

	      search_send_scan_results (results, mpi_master_rank);
	      FreeResults(results);
	    } 
	  else if (job_type == ALN_WORK && do_align) 
	    {
	      CYKDivideAndConquer(cm, seq, seqlen, bestr, 1, seqlen, &tr,
				  cm->dmin, cm->dmax);
	      search_send_align_results (tr, mpi_master_rank);
	      FreeParsetree(tr);
	    }
	  if (seq != NULL)
	    free(seq);
	  seq = NULL;
	} while (job_type != TERMINATE_WORK);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  /*printf("E PSD rank: %4d mast: %4d\n", mpi_my_rank, mpi_master_rank);*/
}
#endif


/*
 * Function: read_next_seq
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easeled: EPN, Fri Dec  8 11:40:20 2006
 * Purpose:  Given a dbfp and whether or not to take the reverse complement,
 *           reads in the sequence and prepares reverse complement.
 */
db_seq_t *read_next_seq (ESL_SQFILE *dbfp, int do_revcomp) 
{
  db_seq_t *ret_dbseq;
  int status;
  char *tmp_seq;

  ret_dbseq = MallocOrDie(sizeof(db_seq_t));

  ret_dbseq->sq[0] = esl_sq_Create();
  status = (esl_sqio_Read(dbfp, ret_dbseq->sq[0]) == eslOK);

  while(status && ret_dbseq->sq[0]->n == 0) /* skip zero length seqs */
    {
      esl_sq_Reuse(ret_dbseq->sq[0]);
      status = (esl_sqio_Read(dbfp, ret_dbseq->sq[0]) == eslOK);
    }
  if(!status)
    return NULL;

  s2upper(ret_dbseq->sq[0]->seq);                      /* Makes uppercase */
  /* Following line will be unnecessary once ESL_SQ objects have dsq's implemented
   * (i.e. allocated and filled within a esl_sqio_Read() call */
  ret_dbseq->sq[0]->dsq = DigitizeSequence (ret_dbseq->sq[0]->seq, ret_dbseq->sq[0]->n);

  if (do_revcomp)
    {
      /* make a new SQ object, to store the reverse complement */
      tmp_seq = MallocOrDie(sizeof(char) * (ret_dbseq->sq[0]->n+1));
      revcomp(tmp_seq, ret_dbseq->sq[0]->seq);
      ret_dbseq->sq[1] = esl_sq_CreateFrom(ret_dbseq->sq[0]->name, tmp_seq, 
					   ret_dbseq->sq[0]->desc, ret_dbseq->sq[0]->acc, 
					   ret_dbseq->sq[0]->ss);
      free(tmp_seq);
      /* Following line will be unnecessary once ESL_SQ objects have dsq's implemented
       * (i.e. allocated and filled within a esl_sq_CreateFrom() call */
      ret_dbseq->sq[1]->dsq = DigitizeSequence (ret_dbseq->sq[1]->seq, ret_dbseq->sq[1]->n);
    }
  ret_dbseq->results[0] = NULL;
  ret_dbseq->results[1] = NULL;

  return(ret_dbseq);
}

/* EPN, Mon Jan  8 06:42:59 2007
 * 
 * Function: actually_search_target
 * 
 * Purpose:  Given a CM and a sequence, call the correct search algorithm
 *           based on cm->search_opts.
 * 
 * Args:     CM         - the covariance model
 *           dsq        - the digitized target sequence
 *           i0         - start of target subsequence (often 1, beginning of dsq)
 *           j0         - end of target subsequence (often L, end of dsq)
 *           cm_cutoff  - minimum CM  score to report 
 *           cp9_cutoff - minimum CP9 score to report (or keep if filtering)
 *           results    - scan_results_t to add to; if NULL, don't add to it
 *           do_filter  - TRUE if we should filter, but only if cm->search_opts tells us to 
 *           doing_cm_stats - TRUE if the reason we're scanning this seq is to build
 *                            a histogram to calculate EVDs for the CM, in this
 *                            case we don't filter regardless of what cm->search_opts says.
 *           doing_cp9_stats- TRUE if we're calc'ing stats for the CP9, in this 
 *                            case we always run CP9ForwardScan()
 *           ret_flen   - RETURN: subseq len that survived filter (NULL if not filtering)
 * Returns: Highest scoring hit from search (even if below cutoff).
 */
float actually_search_target(CM_t *cm, char *dsq, int i0, int j0, float cm_cutoff, 
			     float cp9_cutoff, scan_results_t *results, int do_filter, 
			     int doing_cm_stats, int doing_cp9_stats, int *ret_flen)
{
  float sc;
  int flen;

  /*printf("in actually_search_target: i0: %d j0: %d do_filter: %d doing_cm_stats: %d doing_cp9_stats: %d\n", i0, j0, do_filter, doing_cm_stats, doing_cp9_stats);
    printf("\ti0: %d j0: %d filter: %d\n", i0, j0, do_filter);*/
  flen = (j0-i0+1);

  if(doing_cm_stats && doing_cp9_stats)
    Die("ERROR in actually_search_target doing_cm_stats and doing_cp9_stats both TRUE.\n");
  
  /* check for CP9 related (either filtering or HMMONLY) options first */

  if(doing_cp9_stats || 
     ((do_filter) && 
     (!doing_cm_stats) &&  /* if we're doing CM stats, don't filter. */
     ((cm->search_opts & CM_SEARCH_HMMONLY) ||
      (cm->search_opts & CM_SEARCH_HMMFB) ||
      (cm->search_opts & CM_SEARCH_HMMWEINBERG))))
    {
      sc = CP9Scan_dispatch(cm, dsq, i0, j0, cm->W, cm_cutoff, cp9_cutoff, results, doing_cp9_stats, ret_flen);;
    }
  else
    {
      if(cm->search_opts & CM_SEARCH_NOQDB)
	if(cm->search_opts & CM_SEARCH_INSIDE)
	sc = iInsideScan(cm, dsq, i0, j0, cm->W, cm_cutoff, results);
      else /* don't do inside */
	sc = CYKScan (cm, dsq, i0, j0, cm->W, cm_cutoff, results);
    else /* use QDB */
      if(cm->search_opts & CM_SEARCH_INSIDE)
	sc = iInsideBandedScan(cm, dsq, cm->dmin, cm->dmax, i0, j0, cm->W, cm_cutoff, results);
      else /* do't do inside */
	sc = CYKBandedScan (cm, dsq, cm->dmin, cm->dmax, i0, j0, cm->W, cm_cutoff, results);
    }    
  /*printf("returning from actually_search_target, sc: %f\n", sc);*/
  return sc;
}  

/*
 * Function: print_results
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easelfied: EPN, Fri Dec  8 08:29:05 2006 
 * Purpose:  Given the needed information, prints the results.
 *
 *           cm                  the model
 *           cons                consensus seq for model (query seq)
 *           dbseq               the database seq
 *           name                sequence name
 *           len                 length of the sequence
 *           in_revcomp          are we doing the minus strand
 *           do_stats            should we calculate stats?
 *           mu, lambda          for statistics
 */
void print_results (CM_t *cm, CMConsensus_t *cons, db_seq_t *dbseq,
		    int do_complement, int do_stats, double *mu, 
		    double *lambda) 
{
  int i;
  char *name;
  int len;
  scan_results_t *results;
  Fancyali_t *ali;
  int in_revcomp;
  int header_printed = 0;
  int gc_comp;
  float score_for_Eval; /* the score we'll determine the statistical signifance
			 * of. This will be the bit score stored in
			 * dbseq unless (cm->flags & CM_ENFORCED)
			 * in which case we subtract cm->enf_scdiff first. */

  name = dbseq->sq[0]->name;
  len = dbseq->sq[0]->n;

  for (in_revcomp = 0; in_revcomp <= do_complement; in_revcomp++) 
    {
      results = dbseq->results[in_revcomp];
      if (results == NULL || results->num_results == 0) continue;
      
      if (!header_printed) 
	{
	  header_printed = 1;
	  printf (">%s\n\n", name);
	}
      printf ("  %s strand results:\n\n", in_revcomp ? "Minus" : "Plus");
      
      for (i=0; i<results->num_results; i++) 
	{
	  gc_comp = get_gc_comp (dbseq->sq[in_revcomp]->seq, 
				 results->data[i].start, results->data[i].stop);
	  printf (" Query = %d - %d, Target = %d - %d\n", 
		  cons->lpos[cm->ndidx[results->data[i].bestr]]+1,
		  cons->rpos[cm->ndidx[results->data[i].bestr]]+1,
		  coordinate(in_revcomp, results->data[i].start, len), 
		  coordinate(in_revcomp, results->data[i].stop, len));
	  if (do_stats) 
	    {
	      score_for_Eval = results->data[i].score;
	      if(cm->flags & CM_ENFORCED)
		{
		  printf("\n\torig sc: %.3f", score_for_Eval);
		  score_for_Eval -= cm->enf_scdiff;
		  printf(" new sc: %.3f (diff: %.3f\n\n", score_for_Eval, cm->enf_scdiff);
		}
	      printf (" Score = %.2f, E = %.4g, P = %.4g, GC = %3d\n", results->data[i].score,
		      RJK_ExtremeValueE(score_for_Eval, mu[gc_comp], 
					lambda[gc_comp]),
		      esl_gumbel_surv((double) score_for_Eval, mu[gc_comp], 
				      lambda[gc_comp]), gc_comp);
	      /*printf("  Mu[gc=%d]: %f, Lambda[gc=%d]: %f\n", gc_comp, mu[gc_comp], gc_comp,
		lambda[gc_comp]);*/
	      /*ExtremeValueP(results->data[i].score, mu[gc_comp], 
		lambda[gc_comp]));*/
	    } 
	  else 
	    {
	      printf (" Score = %.2f, GC = %3d\n", results->data[i].score, gc_comp);
	    }
	  printf ("\n");
	  if (results->data[i].tr != NULL) 
	    {
	      ali = CreateFancyAli (results->data[i].tr, cm, cons, 
				    dbseq->sq[in_revcomp]->dsq +
				    (results->data[i].start-1));
	      PrintFancyAli(stdout, ali);
	      FreeFancyAli(ali);
	      printf ("\n");
	    }
	}
    }
  fflush(stdout);
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
 * Function: remove_hits_over_e_cutoff
 * Date:     RJK, Tue Oct 8, 2002 [St. Louis]
 * Purpose:  Given an E-value cutoff, lambdas, mus, a sequence, and
 *           a list of results, calculates GC content for each hit, 
 *           calculates E-value, and decides wheter to keep hit or not.
 */
void remove_hits_over_e_cutoff (CM_t *cm, scan_results_t *results, char *seq,
				float cutoff, double *lambda, double *mu) 
{
  int gc_comp;
  int i, x;
  scan_result_node_t swap;
  float score_for_Eval; /* the score we'll determine the statistical signifance
			 * of. This will be the bit score stored in
			 * dbseq unless (cm->flags & CM_ENFORCED)
			 * in which case we subtract cm->enf_scdiff first. */

  /*printf("in remove_hits_over_e_cutoff()\n");*/
  if (results == NULL)
    return;

  for (i=0; i<results->num_results; i++) 
    {
      gc_comp = get_gc_comp (seq, results->data[i].start, results->data[i].stop);
      score_for_Eval = results->data[i].score;
      if(cm->flags & CM_ENFORCED)
	{
	  /*printf("\n\tRM orig sc: %.3f", score_for_Eval);*/
	  score_for_Eval -= cm->enf_scdiff;
	  /*printf(" new sc: %.3f (diff: %.3f\n", score_for_Eval, cm->enf_scdiff);*/
	}
      /*printf("score_for_Eval: %f \n", score_for_Eval);*/
      if (RJK_ExtremeValueE(score_for_Eval,
			    mu[gc_comp], lambda[gc_comp]) > cutoff) 
	{
	  results->data[i].start = -1;
	}
      /*printf("Eval: %f, start: %d\n", RJK_ExtremeValueE(score_for_Eval,
	mu[gc_comp], lambda[gc_comp]),
	results->data[i].start);*/
    }
  
  for (x=0; x<results->num_results; x++) 
    {
      while (results->num_results > 0 && 
	     results->data[results->num_results-1].start == -1)
	results->num_results--;
      if (x<results->num_results && results->data[x].start == -1) 
	{
	  swap = results->data[x];
	  results->data[x] = results->data[results->num_results-1];
	  results->data[results->num_results-1] = swap;
	  results->num_results--;
	}
    }
  while (results->num_results > 0 &&
	 results->data[results->num_results-1].start == -1)
    results->num_results--;
  sort_results(results);
}  


/*
 * Function: read_next_aln_seqs
 * Date:     EPN, Fri Dec 29 17:44:25 2006
 *
 * Purpose:  Given a pointer to a seq file we're reading seqs to align
 *           from, read in nseq seqs from the seq file. 
 */
seqs_to_aln_t * read_next_aln_seqs(ESL_SQFILE *seqfp, int nseq, int index) 
{
  seqs_to_aln_t *ret_seqs_to_aln;
  int status;
  int i;

  ret_seqs_to_aln = MallocOrDie(sizeof(seqs_to_aln_t));
  ret_seqs_to_aln->sq = MallocOrDie(sizeof(ESL_SQ *) * nseq);
  /*ret_seqs_to_aln->sq = MallocOrDie(sizeof(ESL_SQ *) * nseq);*/
  
  for(i=0; i < nseq; i++)
  {
    ret_seqs_to_aln->sq[i] = esl_sq_Create();
    status = (esl_sqio_Read(seqfp, ret_seqs_to_aln->sq[i]) == eslOK);
    while(status && ret_seqs_to_aln->sq[i]->n == 0) /* skip zero length seqs */
      {
	esl_sq_Reuse(ret_seqs_to_aln->sq[i]);
	status = (esl_sqio_Read(seqfp, ret_seqs_to_aln->sq[i]) == eslOK);
      }
    if(!status)
      {
	if(i == 0) return NULL; /* we're at the end of the file and we didn't read any sequences. */
	else break; /* we're at the end of the file, but we read some sequences */
      }
    /* Following line will be unnecessary once ESL_SQ objects have dsq's implemented
     * (i.e. allocated and filled within a esl_sqio_Read() call) */
    ret_seqs_to_aln->sq[i]->dsq = DigitizeSequence (ret_seqs_to_aln->sq[i]->seq, ret_seqs_to_aln->sq[i]->n);
  }

  ret_seqs_to_aln->nseq = i; /* however many seqs we read, up to nseq */
  ret_seqs_to_aln->tr   = NULL;
  ret_seqs_to_aln->postcode = NULL;
  ret_seqs_to_aln->index = index;
  return(ret_seqs_to_aln);
}


/* EPN, Tue Dec  5 14:25:02 2006
 * 
 * Function: serial_align_targets()
 * 
 * Purpose:  Given a CM and a sequence file name, do preliminaries, align w/the correct 
 *           alignment function and print out the alignment.
 * 
 * Args:     seqfp        - the open sequence file
 *           CM           - the covariance model
 *           ret_sq       - RETURN: the sequences (EASEL)
 *           ret_tr       - RETURN: the parsetrees for seqs in seqfp
 *           ret_postcode - RETURN: the postal codes (NULL if not doing posteriors)
 *           ret_nseq     - RETURN: the number of seqs in seqfp
 *           bdump_level  - verbosity level for band related print statements
 *           debug_level  - verbosity level for debugging print statements
 *           silent_mode  - TRUE to not print anything, FALSE to print scores 
 * 
 */
void
serial_align_targets(ESL_SQFILE *seqfp, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, 
		     char ***ret_postcode, int *ret_nseq, int bdump_level, int debug_level, 
		     int silent_mode)
{
  Parsetree_t    **tr;          /* parse trees for the sequences */
  ESL_SQ         **sq;          /* the sequences */
  int              nalloc;      /* seqs allocated thus far */
  int              i;           /* seq index */
  int              status;
  char           **postcode;
  int              nseq;

  /*printf("in serial_align_targets\n");*/

  /*****************************************************************
   * Read the sequences from the open sequence file. 
   *****************************************************************/
  i = 0;
  nalloc = 10;
  sq = MallocOrDie(sizeof(ESL_SQ *) * nalloc);
  sq[i] = esl_sq_Create();
  while ((status = esl_sqio_Read(seqfp, sq[i])) == eslOK)
  {
    /*printf("Read %12s: length %d\n", sq[i]->name, sq[i]->n);*/
    /* Following line will be unnecessary once ESL_SQ objects have dsq's implemented
     * (i.e. allocated and filled within a esl_sqio_Read() call */
    sq[i]->dsq = DigitizeSequence (sq[i]->seq, sq[i]->n);

    if(++i == nalloc)
    {
      nalloc += 10;
      sq = ReallocOrDie(sq, (sizeof(ESL_SQ *) * nalloc));
    }
    sq[i] = esl_sq_Create();
  }
  /* destroy the last sequence that was alloc'ed but not filled */
  esl_sq_Destroy(sq[i]);
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %d, file %s:\n%s", 
	      seqfp->linenumber, seqfp->filename, seqfp->errbuf);
  nseq = i;

  /***************************************************************
   * Align all the sequences and collect parsetrees with 1 call
   * to align_target_seqs()
   ****************************************************************/
  actually_align_targets(cm, sq, nseq, &tr, &postcode, bdump_level, debug_level, silent_mode);

  /* Clean up and return */
  *ret_tr = tr;
  if((cm->align_opts & CM_ALIGN_POST) && ret_postcode != NULL) *ret_postcode = postcode;
  *ret_nseq = nseq;
  *ret_sq   = sq;
  /*printf("leaving serial_align_targets\n");*/
  return;
}

#ifdef USE_MPI
/* EPN, Tue Dec  5 14:25:02 2006
 * 
 * Function: parallel_align_targets
 * 
 * Purpose:  Given a CM, an open ESL_SQ sequence file, do preliminaries, and align
 *           seqs in parallel. After finishing, print out the alignment.
 * 
 * Args:     seqfp        - the sequence file with target seqs to align
 *           CM           - the covariance model
 *           seqfp        - the sequence file with target seqs to align
 *           ret_sq       - RETURN: the sequences (EASEL)
 *           ret_tr       - RETURN: the parsetrees for seqs in seqfp
 *           ret_postcode - RETURN: the postal codes (NULL if not doing posteriors)
 *           ret_nseq     - RETURN: the number of seqs in seqfp
 *           bdump_level  - verbosity level for band related print statements
 *           debug_level  - verbosity level for debugging print statements
 *           silent_mode  - TRUE to not print anything, FALSE to print scores 
 *           
 *           mpi_my_rank  - rank of current processor 
 *           mpi_master_rank - the master's rank
 *           mpi_num_procs - number of processes 
 */
void
parallel_align_targets(ESL_SQFILE *seqfp, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr,
		       char ***ret_postcode, int *ret_nseq, int bdump_level, int debug_level,
		       int silent_mode, int mpi_my_rank, int mpi_master_rank, int mpi_num_procs)
{
  char job_type;
  seqs_to_aln_t **active_seqs;
  int *process_status;
  int eof = FALSE;
  int proc_index, active_seq_index;
  int nseq_per_job = 1;
  seqs_to_aln_t *seqs_to_aln;
  int i;
  Parsetree_t    **all_tr;          /* parse trees for the all the sequences in order they were read */
  ESL_SQ         **all_sq;          /* all sequences in the order they were read                     */
  char           **all_postcode;    /* post codes for all the sequences in order they were read      */
  int              nseq_read;       /* number of sequences read overall */
  int              nalloc;
  int              alloc_chunk = 2; 
  int              do_post;

  if(cm->align_opts & CM_ALIGN_POST)
    do_post = TRUE;
  else
    do_post = FALSE;

  /*printf("in parallel_align_targets rank: %d master: %d nprocs: %d do_post: %d\n", mpi_my_rank, mpi_master_rank, mpi_num_procs, do_post);*/
  if (mpi_my_rank == mpi_master_rank) 
    {
      nalloc = alloc_chunk;
      all_tr = MallocOrDie(sizeof(Parsetree_t *) * nalloc);
      all_sq = MallocOrDie(sizeof(ESL_SQ *)      * nalloc);
      all_postcode = MallocOrDie(sizeof(char *)  * nalloc);
      nseq_read  = 0;

      /* Set up arrays to hold pointers to active seqs and jobs on
	 processes */
      active_seqs = MallocOrDie(sizeof(seqs_to_aln_t *) * mpi_num_procs);
      process_status = MallocOrDie(sizeof(int) * mpi_num_procs);
      for (active_seq_index=0; active_seq_index<mpi_num_procs; 
	   active_seq_index++) 
	active_seqs[active_seq_index] = NULL;
      for (proc_index = 0; proc_index < mpi_num_procs; proc_index++)
	process_status[proc_index] = IDLE;
      
      do
	{
	  /* Check for idle processes.  Send jobs */
	  for (proc_index=0; proc_index<mpi_num_procs; proc_index++) 
	    {
	      if (proc_index == mpi_master_rank) continue;  /* Skip master process */
	      if (process_status[proc_index] == IDLE) 
		{         
		  /* I'm idle -- need a job */
		  /* Find next non-master open process */
		  /*for (active_seq_index=0; active_seq_index<mpi_num_procs; active_seq_index++) 
		    if (active_seqs[active_seq_index] == NULL) break;
		    if (active_seq_index == mpi_num_procs) 
		    Die ("Tried to read more than %d seqs at once\n", mpi_num_procs);*/

		  active_seqs[proc_index] = read_next_aln_seqs(seqfp, nseq_per_job, nseq_read);

		  if (active_seqs[proc_index] == NULL) 
		    {
		      eof = TRUE;
		      break;            /* No more seqs */
		    }
		  else
		    {
		      for(i = nseq_read; i < (nseq_read + active_seqs[proc_index]->nseq); i++)
			{
			  if(i == nalloc)
			    {
			      nalloc += alloc_chunk;
			      all_sq       = ReallocOrDie(all_sq, sizeof(ESL_SQ *) * nalloc);
			      all_tr       = ReallocOrDie(all_tr, sizeof(Parsetree_t *) * nalloc);
			      all_postcode = ReallocOrDie(all_postcode, sizeof(char *) * nalloc);
			    }
			  all_sq[i] = active_seqs[proc_index]->sq[(i - nseq_read)];
			}
		      nseq_read += active_seqs[proc_index]->nseq;
		      aln_send_next_job(active_seqs[proc_index], proc_index);
		      /* Set this process as having the job */
		      process_status[proc_index] = BUSY;
		    }
		} 
	    }
	  /* Wait for next reply */
	  if (aln_procs_working(process_status, mpi_num_procs, mpi_master_rank)) 
	    {
	      aln_check_results (all_tr, all_postcode, &process_status);
	      /* Free active seq seqs_to_aln_t somehow */
	    }
	} while (!eof || (aln_procs_working(process_status, mpi_num_procs, mpi_master_rank)));
      for (proc_index=0; proc_index<mpi_num_procs; proc_index++) 
	if (proc_index != mpi_master_rank)
	  aln_send_terminate (proc_index);
      free(active_seqs);
      free(process_status);
      *ret_tr = all_tr;
      *ret_postcode = all_postcode;
      *ret_sq = all_sq;
      *ret_nseq = nseq_read;
    } 
  else  /* not the master node */
    {
      do 
	{
	  job_type = aln_receive_job (&seqs_to_aln, mpi_master_rank);
	  if (job_type == ALN_WORK) 
	    {
	      silent_mode = FALSE;
	      debug_level = 0;
	      bdump_level = 0;
	      /* align the targets */
	      actually_align_targets(cm, seqs_to_aln->sq, seqs_to_aln->nseq, &(seqs_to_aln->tr), 
				     &(seqs_to_aln->postcode), bdump_level, debug_level, silent_mode);

	      /*printf("done actually_align_targets\n");*/
	      aln_send_results(seqs_to_aln, do_post, mpi_master_rank);
	      /*Free_seqs_to_aln(seqs_to_aln); */
	    }		  
	} while (job_type != TERMINATE_WORK);
      ret_tr = NULL;
      ret_sq = NULL;
      ret_nseq = NULL;
      ret_postcode = NULL;
    }
  MPI_Barrier(MPI_COMM_WORLD);

  /*printf("leaving parallel_align_targets rank: %4d mast: %4d\n", mpi_my_rank, mpi_master_rank);*/
}
#endif

/* EPN, Tue Dec  5 14:25:02 2006
 * 
 * Function: actually_align_targets
 * 
 * Purpose:  Given a CM and sequences, do preliminaries, call the correct 
 *           CYK function and return parsetrees and optionally postal codes 
 *           (if cm->align_opts & CM_ALIGN_POST)
 * 
 * Args:     CM           - the covariance model
 *           sq           - the sequences
 *           nseq         - number of seqs we're aligning
 *           ret_tr       - RETURN: parsetrees (pass NULL if trace isn't wanted)
 *           ret_postcode - RETURN: postal code string
 *           bdump_level  - verbosity level for band related print statements
 *           debug_level  - verbosity level for debugging print statements
 *           silent_mode  - TRUE to not print anything, FALSE to print scores 
 */
void
actually_align_targets(CM_t *cm, ESL_SQ **sq, int nseq, Parsetree_t ***ret_tr, char ***ret_postcode,
		       int bdump_level, int debug_level, int silent_mode)
{
  Stopwatch_t      *watch;      /* for timings */
  int i;                        /* counter over sequences */
  int v;                        /* state counter */
  char           **postcode;    /* posterior decode array of strings        */
  Parsetree_t    **tr;          /* parse trees for the sequences */
  float            sc;		/* score for one sequence alignment */
  float            maxsc;	/* max score in all seqs */
  float            minsc;	/* min score in all seqs */
  float            avgsc;	/* avg score over all seqs */
  float            tmpsc;       /* temporary score */

  /* variables related to CM Plan 9 HMMs */
  CP9_t       *hmm;                     /* constructed CP9 HMM */
  CP9Bands_t  *cp9b;                    /* data structure for hmm bands (bands on the hmm states) 
				         * and arrays for CM state bands, derived from HMM bands*/
  CP9Map_t       *cp9map;        /* maps the hmm to the cm and vice versa */
  CP9_dpmatrix_t *cp9_post;      /* growable DP matrix for posterior decode              */
  CP9_dpmatrix_t *cp9_mx;        /* growable DP matrix for viterbi                       */
  float           swentry;	 /* S/W aggregate entry probability       */
  float           swexit;        /* S/W aggregate exit probability        */

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
  Parsetree_t     *orig_tr;        /* parsetree for the orig_cm; created from the sub_cm parsetree */

  CP9_t           *sub_hmm;        /* constructed CP9 HMM; written to hmmfile              */
  CP9Map_t        *sub_cp9map;     /* maps the sub_hmm to the sub_cm and vice versa */
  CP9_t           *orig_hmm;       /* original CP9 HMM built from orig_cm */
  CP9Map_t        *orig_cp9map;    
  CP9Bands_t      *orig_cp9b; 

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

  int do_local   = FALSE;
  int do_qdb     = FALSE;
  int do_hbanded = FALSE;
  int use_sums   = FALSE;
  int do_sub     = FALSE;
  int do_fullsub = FALSE;
  int do_hmmonly = FALSE;
  int do_scoreonly = FALSE;
  int do_inside  = FALSE;
  int do_outside = FALSE;
  int do_small   = TRUE;
  int do_post    = FALSE;
  int do_timings = FALSE;
  int do_check   = FALSE;

  /*printf("in actually_align_targets\n");*/

  /* set the options based on cm->align_opts */
  if(cm->config_opts & CM_CONFIG_LOCAL)     do_local   = TRUE;
  if(cm->align_opts  & CM_ALIGN_QDB)        do_qdb     = TRUE;
  if(cm->align_opts  & CM_ALIGN_HBANDED)    do_hbanded = TRUE;
  if(cm->align_opts  & CM_ALIGN_SUMS)       use_sums   = TRUE;
  if(cm->align_opts  & CM_ALIGN_SUB)        do_sub     = TRUE;
  if(cm->align_opts  & CM_ALIGN_FSUB)       do_fullsub = TRUE;
  if(cm->align_opts  & CM_ALIGN_HMMONLY)    do_hmmonly = TRUE;
  if(cm->align_opts  & CM_ALIGN_INSIDE)     do_inside  = TRUE;
  if(cm->align_opts  & CM_ALIGN_OUTSIDE)    do_outside = TRUE;
  if(cm->align_opts  & CM_ALIGN_NOSMALL)    do_small   = FALSE;
  if(cm->align_opts  & CM_ALIGN_POST)       do_post    = TRUE;
  if(cm->align_opts  & CM_ALIGN_TIME)       do_timings = TRUE;
  if(cm->align_opts  & CM_ALIGN_CHECKINOUT) do_check   = TRUE;
  if(cm->align_opts  & CM_ALIGN_SCOREONLY)  do_scoreonly = TRUE;

  if(do_fullsub)
    {
      do_sub = TRUE;
      cm->align_opts |= CM_ALIGN_SUB;
    }
  /*printf("do_local  : %d\n", do_local);
    printf("do_qdb    : %d\n", do_qdb);
    printf("do_hbanded: %d\n", do_hbanded);
    printf("use_sums  : %d\n", use_sums);
    printf("do_sub    : %d\n", do_sub);
    printf("do_fsub   : %d\n", do_fullsub);
    printf("do_hmmonly: %d\n", do_hmmonly);
    printf("do_inside : %d\n", do_inside);
    printf("do_outside: %d\n", do_outside);
    printf("do_small  : %d\n", do_small);
    printf("do_post   : %d\n", do_post);
    printf("do_timings: %d\n", do_timings);*/
    
  tr    = MallocOrDie(sizeof(Parsetree_t) * nseq);
  minsc = FLT_MAX;
  maxsc = -FLT_MAX;
  avgsc = 0;
  watch = StopwatchCreate(); 

  if(do_hbanded || do_sub) /* We need a CP9 HMM to build sub_cms */
    {
      /* Keep this data for the original CM safe; we'll be doing
       * pointer swapping to ease the sub_cm alignment implementation. */
      hmm         = cm->cp9;
      cp9map      = cm->cp9map;

      orig_hmm    = hmm;
      orig_cp9map = cp9map;
    }

  /* Copy the QD bands in case we expand them. */
  if(do_qdb)
    {
      if(bdump_level > 1) 
	  /*printf("cm->beta:%f\n", cm->beta);*/
	  debug_print_bands(cm, cm->dmin, cm->dmax);
      expand_flag = FALSE;
      /* Copy dmin and dmax, so we can replace them after expansion */
      orig_dmin = MallocOrDie(sizeof(int) * cm->M);
      orig_dmax = MallocOrDie(sizeof(int) * cm->M);
      for(v = 0; v < cm->M; v++)
	{
	  orig_dmin[v] = cm->dmin[v];
	  orig_dmax[v] = cm->dmax[v];
	}
    }	  

  if(do_sub) /* to get spos and epos for the sub_cm, 
	      * we config the HMM to local mode with equiprobably start/end points.*/
      {
	/*printf("configuring the CM plan 9 HMM for local alignment.\n");*/
	swentry= ((hmm->M)-1.)/hmm->M; /* all start pts equiprobable, including 1 */
	swexit = ((hmm->M)-1.)/hmm->M; /* all end   pts equiprobable, including M */
	CPlan9SWConfig(hmm, swentry, swexit);
	CP9Logoddsify(hmm);
      }

  if(do_post)
    postcode = malloc(sizeof(char *) * nseq);
  if(do_hbanded)
    {
      cp9b = AllocCP9Bands(cm, cm->cp9);
      orig_cp9b = cp9b; 
    }
  orig_cm = cm;
  /*****************************************************************
   *  Collect parse trees for each sequence
   *****************************************************************/

  for (i = 0; i < nseq; i++)
    {
      StopwatchZero(watch);
      StopwatchStart(watch);

      if (sq[i]->n == 0) continue;

      /* Special case, if do_hmmonly, align seq with Viterbi, print score and move 
       * on to next seq */
      if(do_hmmonly)
	{
	  cp9_mx  = CreateCPlan9Matrix(1, cm->cp9->M, 25, 0);
	  if(!silent_mode) printf("Aligning (to a CP9 HMM w/viterbi) %-20s", sq[i]->name);
	  sc = CP9Viterbi(sq[i]->dsq, 1, sq[i]->n, cm->cp9, cp9_mx);
	  if(!silent_mode) printf(" score: %10.2f bits\n", sc);
	  FreeCPlan9Matrix(cp9_mx);
	  continue;
	}
      /* Special case, if do_scoreonly, align seq with full CYK inside, just to 
       * get the score. For testing, probably in cmscore. */
      if(do_scoreonly)
	{
	  if(!silent_mode) printf("Aligning (w/full CYK score only) %-30s", sq[i]->name);
	  sc = CYKInsideScore(cm, sq[i]->dsq, 0, 1, sq[i]->n, sq[i]->n,
			      NULL, NULL); /* don't do QDB mode */
	  if(!silent_mode) printf("    score: %10.2f bits\n", sc);
	  continue;
	}

      /* Potentially, do HMM calculations. */
      if(do_hbanded)
	{
	  if(do_sub)
	    CP9_seq2bands(orig_cm, sq[i]->dsq, 1, sq[i]->n, orig_cp9b, 
			  &cp9_post, /* we DO want the posterior matrix back */
			  debug_level);
	  else
	    CP9_seq2bands(orig_cm, sq[i]->dsq, 1, sq[i]->n, orig_cp9b, 
			  NULL, /* we don't want the posterior matrix back */
			  debug_level);
	}
      /* If we're in sub mode:
       * (1) Get HMM posteriors (if do_hbanded, we already have them) 
       * (2) Infer the start (spos) and end (epos) HMM states by 
       *     looking at the posterior matrix.
       * (3) Build the sub_cm from the original CM.
       *
       * If we're also doing HMM banded alignment:
       * (4) Build a new CP9 HMM from the sub CM.
       * (5) Do Forward/Backward again, and get a new posterior matrix.
       */
      if(do_sub)
	{
	  /* (1) Get HMM posteriors (if do_hbanded, we already have them) */
	  if(!do_hbanded) 
	    CP9_seq2posteriors(orig_cm, sq[i]->dsq, 1, sq[i]->n, &cp9_post, debug_level); 
	  
	  /* (2) infer the start and end HMM nodes (consensus cols) from posterior matrix.
	   * Remember: we're necessarily in CP9 local mode, the --sub option turns local mode on. 
	   */
	  CP9NodeForPosn(orig_hmm, 1, sq[i]->n, 1,             cp9_post, &spos, &spos_state, 
			 do_fullsub, 0., TRUE, debug_level);
	  CP9NodeForPosn(orig_hmm, 1, sq[i]->n, sq[i]->n, cp9_post, &epos, &epos_state, 
			 do_fullsub, 0., FALSE, debug_level);
	  /* If the most likely state to have emitted the first or last residue
	   * is the insert state in node 0, it only makes sense to start modelling
	   * at consensus column 1. */
	  if(spos == 0 && spos_state == 1) 
	      spos = 1;
	  if(epos == 0 && epos_state == 1) 
	      epos = 1;
	  if(epos < spos) /* This is a possible but hopefully rarely encountered situation. */
	    epos = spos;
	  
	  /* (3) Build the sub_cm from the original CM. */
	  if(!(build_sub_cm(orig_cm, &sub_cm, 
			    spos, epos,         /* first and last col of structure kept in the sub_cm  */
			    &submap,            /* maps from the sub_cm to cm and vice versa           */
			    do_fullsub,         /* build or not build a sub CM that models all columns */
			    debug_level)))      /* print or don't print debugging info                 */
	    Die("Couldn't build a sub CM from the CM\n");
	  /* Configure the sub_cm, the same as the cm, this will build a CP9 HMM if (do_hbanded) */
	  ConfigCM(sub_cm, NULL, NULL);

	  cm    = sub_cm; /* orig_cm still points to the original CM */
	  if(do_hbanded) /* we're doing HMM banded alignment to the sub_cm */
	    {
	      /* Get the HMM bands for the sub_cm */
	      sub_hmm    = sub_cm->cp9;
	      sub_cp9map = sub_cm->cp9map;
	      sub_cp9b   = AllocCP9Bands(sub_cm, sub_cm->cp9);
	      CP9_seq2bands(sub_cm, sq[i]->dsq, 1, sq[i]->n, sub_cp9b, 
			    NULL, /* we don't want posterior matrix back */
			    debug_level);
	      hmm           = sub_hmm;    
	      cp9map        = sub_cp9map;
	      cp9b          = sub_cp9b;
	    }
	}
      /* Determine which CYK alignment algorithm to use, based
       * on command-line options AND memory requirements.
       */
      if(do_hbanded)
	{
	  /* write a function to determine size of jd banded memory
	   * req'd, and set do_small to true if its > thresh.
	   if(do_small) * We're only going to band on d in memory, but 
	   * we need to calculate safe_hd bands on the d dimension. 
	   {
	  */
	}
      if(do_qdb)
	{
	  /*Check if we need to reset the query dependent bands b/c they're currently expanded. */
	  if(expand_flag)
	    {
	      for(v = 0; v < cm->M; v++)
		{
		  cm->dmin[v] = orig_dmin[v];
		  cm->dmax[v] = orig_dmax[v];
		}
	      expand_flag = FALSE;
	    }
	  if((sq[i]->n < cm->dmin[0]) || (sq[i]->n > cm->dmax[0]))
	    {
	      /* the seq we're aligning is outside the root band, so we expand.*/
	      ExpandBands(cm, sq[i]->n, cm->dmin, cm->dmax);
	      if(debug_level > 0) printf("Expanded bands for seq : %s\n", sq[i]->name);
	      if(bdump_level > 2) 
		{
		  printf("printing expanded bands :\n");
		  debug_print_bands(cm, cm->dmin, cm->dmax);
		}
	      expand_flag = TRUE;
	    }
	}

      if(!silent_mode) 
	{
	  if(do_sub) 
	    printf("Aligning (to a sub CM) %-20s", sq[i]->name);
	  else
	    printf("Aligning %-30s", sq[i]->name);
	}
      if (do_inside)
	{
	  if(do_hbanded)
	    {
	      sc = IInside_b_jd_me(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
				   BE_PARANOID,	/* non memory-saving mode */
				   NULL, NULL,	/* manage your own matrix, I don't want it */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,        /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	    }
	  else
	    {
	      sc = IInside(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
			   BE_EFFICIENT,	/* memory-saving mode */
			   NULL, NULL,	/* manage your own matrix, I don't want it */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	    }

	}
      else if(do_outside)
	{	
	  if(do_hbanded)
	    {
	      
	      sc = IInside_b_jd_me(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,        /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      /*do_check = TRUE;*/
	      sc = IOutside_b_jd_me(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
				    BE_PARANOID,	/* save full beta */
				    NULL, NULL,	/* manage your own matrix, I don't want it */
				    NULL, NULL,	/* manage your own deckpool, I don't want it */
				    do_local,       /* TRUE to allow local begins */
				    alpha,          /* alpha matrix from FInside_b_jd_me() */
				    NULL,           /* don't save alpha */
				    do_check,       /* TRUE to check Outside probs agree with Inside */
				    cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	    }
	  else
	    {
	      sc = IInside(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = IOutside(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
			    BE_PARANOID,	/* save full beta */
			    NULL, NULL,	/* manage your own matrix, I don't want it */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha,         /* alpha matrix from IInside() */
			    NULL,           /* don't save alpha */
			    do_check);      /* TRUE to check Outside probs agree with Inside */
	    }
	}
      else if (do_small) 
	{
	  if(do_qdb)
	    {
	      sc = CYKDivideAndConquer(cm, sq[i]->dsq, sq[i]->n, 0, 1, sq[i]->n, 
				       &(tr[i]), cm->dmin, cm->dmax);
	      if(bdump_level > 0)
 		qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
	    }
	  else if(do_hbanded) /*j and d bands not tight enough to allow HMM banded full CYK*/
	    {
	      /* Calc the safe d bands */
	      hd2safe_hd_bands(cm->M, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
			       cp9b->safe_hdmin, cp9b->safe_hdmax);
	      if(debug_level > 3)
		{
		  printf("\nprinting hd bands\n\n");
		  debug_print_hd_bands(cm, cp9b->hdmin, cp9b->hdmax, cp9b->jmin, cp9b->jmax);
		  printf("\ndone printing hd bands\n\n");
		}
	      /* Note the following CYK call will not enforce j bands, even
	       * though user specified --hbanded. */
	      sc = CYKDivideAndConquer(cm, sq[i]->dsq, sq[i]->n, 0, 1, sq[i]->n, 
				       &(tr[i]), cp9b->safe_hdmin, cp9b->safe_hdmax);
	      if(bdump_level > 0)
		qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
	    }
	  else
	    {
	      /*printf("DEBUG PRINTING CM PARAMS BEFORE D&C CALL\n");
		debug_print_cm_params(cm);
		printf("DONE DEBUG PRINTING CM PARAMS BEFORE D&C CALL\n");*/

	      sc = CYKDivideAndConquer(cm, sq[i]->dsq, sq[i]->n, 0, 1, sq[i]->n, &(tr[i]),
				       NULL, NULL); /* we're not in QDB mode */
	      if(bdump_level > 0)
		{
		  /* We want band info but --banded wasn't used.  Useful if you're curious
		   * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
		   * to see where the non-banded parse went outside the bands.
		   */
		  qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
		}
	    }
	}
      else if(do_qdb)
	{
	  sc = CYKInside(cm, sq[i]->dsq, sq[i]->n, 0, 1, sq[i]->n, &(tr[i]), cm->dmin, cm->dmax);
	  if(bdump_level > 0)
	    qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
	}
      else if(do_hbanded)
	{
	  sc = CYKInside_b_jd(cm, sq[i]->dsq, sq[i]->n, 0, 1, sq[i]->n, &(tr[i]), cp9b->jmin, 
			      cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);
	  if(bdump_level > 0)
	    qdb_trace_info_dump(cm, tr[i], cp9b->safe_hdmin, cp9b->safe_hdmax, bdump_level);
	  /* if CM_ALIGN_HMMSAFE option is enabled, realign seqs w/HMM banded parses < 0 bits */
	  if(cm->align_opts & CM_ALIGN_HMMSAFE && sc < 0.)
	    {
	      tmpsc = sc;
	      printf("\n%s HMM banded parse had a negative score, realigning with non-banded CYK.\n", sq[i]->name);
	      FreeParsetree(tr[i]);
	      sc = CYKDivideAndConquer(cm, sq[i]->dsq, sq[i]->n, 0, 1, sq[i]->n, &(tr[i]),
				       NULL, NULL); /* we're not in QDB mode */
	      if(fabs(sc-tmpsc) < 0.01)
		printf("HMM banded parse was the optimal parse.\n\n");
	      else
		printf("HMM banded parse was non-optimal, it was %.2f bits below the optimal.\n\n", (fabs(sc-tmpsc)));
	    }	      
	}
      else
	{
	  sc = CYKInside(cm, sq[i]->dsq, sq[i]->n, 0, 1, sq[i]->n, &(tr[i]), NULL, NULL);
	  if(bdump_level > 0)
	    {
	      /* We want band info but --hbanded wasn't used.  Useful if you're curious
	       * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
	       * to see where the non-banded parse went outside the bands.
	       */
	      qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
	    }
	}
      if(do_post) /* Do Inside() and Outside() runs and use alpha and beta to get posteriors */
	{	
	  /*alpha = MallocOrDie(sizeof(float **) * (cm->M));
	  beta  = MallocOrDie(sizeof(float **) * (cm->M+1));
	  */
	  post  = MallocOrDie(sizeof(int **) * (cm->M+1));
	  /*
	  for (v = 0; v < cm->M; v++) alpha[v] = NULL;
	  for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
	  */
	  if(do_hbanded)
	    {
	      for (v = 0; v < cm->M; v++)
		{
		  post[v] = NULL;
		  post[v] = Ialloc_jdbanded_vjd_deck(sq[i]->n, 1, sq[i]->n, cp9b->jmin[v], 
						      cp9b->jmax[v], cp9b->hdmin[v], cp9b->hdmax[v]);
		}
	      post[cm->M] = NULL;
	      post[cm->M] = Ialloc_vjd_deck(sq[i]->n, 1, sq[i]->n);
	      sc = IInside_b_jd_me(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for IOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,       /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      sc = IOutside_b_jd_me(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
				    BE_PARANOID,	/* save full beta */
				    NULL, &beta,	/* fill beta, and return it, needed for ICMPosterior() */
				    NULL, NULL,	/* manage your own deckpool, I don't want it */
				    do_local,       /* TRUE to allow local begins */
				    alpha, &alpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
				    do_check,      /* TRUE to check Outside probs agree with Inside */
				    cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      ICMPosterior_b_jd_me(sq[i]->n, cm, alpha, NULL, beta, NULL, post, &post,
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);
	      postcode[i] = ICMPostalCode_b_jd_me(cm, sq[i]->n, post, tr[i],
						  cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);
	      /*postcode[i] = CMPostalCode_b_jd_me(cm, sq[i]->n, post, tr[i],
		cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);*/
	    }
	  else
	    {
	      for (v = 0; v < cm->M+1; v++)
		{
		  post[v] = NULL;
		  post[v] = Ialloc_vjd_deck(sq[i]->n, 1, sq[i]->n);
		}
	      sc = IInside(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for IOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = IOutside(cm, sq[i]->dsq, sq[i]->n, 1, sq[i]->n,
			    BE_PARANOID,	/* save full beta */
			    NULL, &beta,	/* fill beta, and return it, needed for CMPosterior() */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha, &alpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
			    do_check);      /* TRUE to check Outside probs agree with Inside */
	      ICMPosterior(sq[i]->n, cm, alpha, NULL, beta, NULL, post, &post);
	      if(do_check || TRUE)
		{
		  ICMCheckPosterior(sq[i]->n, cm, post);
		  printf("\nPosteriors checked (I).\n\n");
		}
	      postcode[i] = ICMPostalCode(cm, sq[i]->n, post, tr[i]);
	      /*postcode[i] = CMPostalCode(cm, sq[i]->n, post, tr[i]);*/
	    }

	  /* free post */
	  if(post != NULL)
	    {
	      for (v = 0; v <= (cm->M); v++)
		if (post[v] != NULL) { Ifree_vjd_deck(post[v], 1, sq[i]->n); post[v] = NULL;}
	      free(post);
	    }
	}
      avgsc += sc;
      if (sc > maxsc) maxsc = sc;
      if (sc < minsc) minsc = sc;
      
      if(!silent_mode) printf("    score: %10.2f bits\n", sc);
      
      /* check parsetree score if cm->align_opts & CM_ALIGN_CHECKPARSESC */
      if((cm->align_opts & CM_ALIGN_CHECKPARSESC) &&
	 (!(cm->flags & CM_IS_SUB) && (!(cm->flags & CM_IS_FSUB))))
	{
	  if (fabs(sc - ParsetreeScore(cm, tr[i], sq[i]->dsq, FALSE)) >= 0.01)
	    Die("ERROR in actually_align_target(), alignment score differs from its parse tree's score");
	}

      /* If debug level high enough, print out the parse tree */
      if((cm->align_opts & CM_ALIGN_PRINTTREES) || (debug_level > 2))
	{
	  fprintf(stdout, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], sq[i]->dsq, FALSE));;
	  ParsetreeDump(stdout, tr[i], cm, sq[i]->dsq);
	  fprintf(stdout, "//\n");
	}
      /* Dump the trace with info on i, j and d bands
       * if bdump_level is high enough */
      if(bdump_level > 0 && do_hbanded)
	ijd_banded_trace_info_dump(cm, tr[i], cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, 
				   cp9b->hdmin, cp9b->hdmax, 1);
      
      /* Clean up the structures we use calculating HMM bands, that are allocated
       * differently for each sequence. 
       */
      if(do_hbanded)
	{
	  for(v = 0; v < cm->M; v++)
	    { 
	      free(cp9b->hdmin[v]); 
	      free(cp9b->hdmax[v]);
	    }
	}
      if(do_sub && !(do_inside || do_outside))
	{
	  /* Convert the sub_cm parsetree to a full CM parsetree */
	  if(debug_level > 0)
	    ParsetreeDump(stdout, tr[i], cm, sq[i]->dsq);
	  if(!(sub_cm2cm_parsetree(orig_cm, sub_cm, &orig_tr, tr[i], submap, do_fullsub, debug_level)))
	    {
	      printf("\n\nIncorrectly converted original trace:\n");
	      ParsetreeDump(stdout, orig_tr, orig_cm, sq[i]->dsq);
	      exit(1);
	    }
	  if(debug_level > 0)
	    {
	      printf("\n\nConverted original trace:\n");
	      ParsetreeDump(stdout, orig_tr, orig_cm, sq[i]->dsq);
	    }
	  /* Replace the sub_cm trace with the converted orig_cm trace. */
	  FreeParsetree(tr[i]);
	  tr[i] = orig_tr;
	  
	  FreeSubMap(submap);
	  FreeCM(sub_cm); /* cm and sub_cm now point to NULL */
	  if(do_hbanded)
	    FreeCP9Bands(sub_cp9b);
	}
      StopwatchStop(watch);
      if(do_timings) 
	{ 
	  StopwatchDisplay(stdout, "seq alignment CPU time: ", watch);
	  printf("\n");
	}
    }
  /* Clean up. */
  if(do_hbanded)
    FreeCP9Bands(orig_cp9b);
  if (do_qdb)
    {
      free(orig_dmin);
      free(orig_dmax);
    }
  StopwatchFree(watch);
  
  *ret_tr = tr; 
  if (do_post) *ret_postcode = postcode; 
  else ret_postcode = NULL;
  
  /*  if(ret_post_spos != NULL)
   *ret_post_spos = post_spos;
   if(ret_post_epos != NULL)
   *ret_post_epos = post_epos;
   if(ret_dist_spos != NULL)
   *ret_dist_spos = dist_spos;
   if(ret_dist_epos != NULL)
   *ret_dist_epos = dist_epos;*/
  
  /*printf("leaving actually_align_targets()\n");*/
}


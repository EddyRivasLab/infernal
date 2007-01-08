/*
 * mpifuncs.h
 
 * Basic functions for using MPI in infernal.
 * Functions can be organized into 2 groups, 
 * 1 group starts with "search_", these
 * were copied from Robbie Klein's RSEARCH
 * and in many cases untouched (besides 
 * renaming). The other group starts with
 * "aln_", these were made for MPI alignment
 * (which is easier to implement than search)
 * by copying and morphing Robbie's functions.
 * Functions that don't start with either
 * "search_" or "aln_" are general.
 *
 * Robert J. Klein
 * May 28, 2002
 * 
 * EPN, Thu Jan  4 14:17:06 2007
 */

#ifndef _MPIFUNCS_H
#define _MPIFUNCS_H

#ifdef USE_MPI

#include "mpi.h"
#include "structs.h"

/* Work types */
#define TERMINATE_WORK            0
#define SEARCH_STD_SCAN_WORK      1
#define SEARCH_ALIGN_WORK         2
#define SEARCH_HIST_SCAN_WORK     3
#define ALN_WORK                  4

/* Results types */
#define SEARCH_STD_SCAN_RESULTS   0
#define SEARCH_ALIGN_RESULTS      1
#define SEARCH_HIST_SCAN_RESULTS  2
#define ALN_RESULTS               3

/* Communication tags */
#define JOB_PACKET_TAG                   1
#define SEQ_TAG                          2
#define SEARCH_STD_SCAN_RESULTS_SIZE_TAG 3
#define SEARCH_STD_SCAN_RESULTS_TAG      4
#define SEARCH_HIST_RESULTS_TAG          5
#define SEARCH_ALIGN_RESULTS_SIZE_TAG    6
#define SEARCH_ALIGN_RESULTS_TAG         7
#define ALN_JOB_SIZE_TAG                 8
#define ALN_JOB_TAG                      9
#define ALN_RESULTS_SIZE_TAG             10
#define ALN_RESULTS_TAG                  11
#define ALN_JOB_PACKET_TAG               12

typedef struct _job_t {
  char job_type;
  int seqlen;
  char *dsq;
  int db_seq_index;
  char in_revcomp;
  int index;            /* result index for alignment, seq index for scan */
  int bestr;            /* Best root state -- only used for alignment jobs,
			   0 for scanning jobs */
  struct _job_t *next;
} job_t;


/***************************************************************************
 * General functions, which can be used for either MPI search or alignment *
 ***************************************************************************/
/* Get rank of master process (lowest ranked one that can do I/O 
 * Also checks the version string */
extern int  get_master_rank (MPI_Comm comm, int mpi_my_rank);
extern void broadcast_cm    (CM_t **cm, int mpi_my_rank, int mpi_master_rank);

/**************************************************************
 * MPI search functions ("search_*") Originally from RSEARCH. *
 **************************************************************/

/* First broadcast of information */
extern void search_first_broadcast (int *num_samples, float *W_scale,
				    int mpi_my_rank, int mpi_master_rank);

/* Second broadcast of information */
extern void search_second_broadcast (float *sc_cutoff, float *e_cutoff, int *cutoff_type, 
				     double *mu, double *lambda, 
				     double *K, long *N, int mpi_my_rank, int mpi_master_rank);

/* Get job from master process */
extern char search_receive_job (int *seqlen_p, char **seq_p, int *bestr_p, int mpi_master_rank);

/* Send results of a scan (scan_results_t) */
extern void search_send_scan_results (scan_results_t *results, int mpi_master_node);

/* Send results of an alignment (Parsetree_t *) */
extern void search_send_align_results (Parsetree_t *tr, int mpi_master_node);

/* TRUE if a slave is still working, FALSE if all done */
extern int search_procs_working (job_t **process_status, int mpi_num_procs, int mpi_master_rank);

/* Chunk the seq */
extern job_t *search_enqueue (db_seq_t *active_seq, int db_seq_index, int D,
		       int do_revcomp, int job_type);

/* Put alignments into the queue */
extern void search_enqueue_alignments (job_t **queue, db_seq_t *active_seq, int db_seq_index,
				int do_revcomp, int job_type);

/* Send the next job from the queue */
extern void search_send_next_job (job_t **queue, job_t **process_status, int rank_to_send_to); 

/* Do a blocking call to MPI_Recv until a process finishes, then process results */
extern int search_check_results (db_seq_t **active_seqs, job_t **process_status, int D);

/* Check for results from a histogram scan */
extern int search_check_hist_results (db_seq_t **seq, job_t **process_status, int D);

/* Send histogram scan results */
extern void search_send_hist_scan_results (float score, int mpi_master_node);

/* Send the termination code to rank i */
extern void search_send_terminate (int i);

/**************************************************************
 * MPI alignment functions ("aln_*") (EPN 01.04.07)
 **************************************************************/
extern int  aln_procs_working   (int *process_status, int mpi_num_procs, int mpi_master_rank);
extern void aln_send_next_job   (seqs_to_aln_t *seqs_to_aln, int rank_to_send_to);
extern char aln_receive_job     (seqs_to_aln_t **ret_seqs_to_aln, int mpi_master_rank);
extern void aln_send_results    (seqs_to_aln_t *seqs_to_aln, int do_post, int mpi_master_node);
extern int  aln_check_results   (Parsetree_t **all_parsetrees, char **all_postcodes, 
				 int **process_status);
extern void aln_send_terminate  (int rank_to_send_to);

#endif
#endif

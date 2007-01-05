/*
 * mpifuncs.c
 *
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

#ifdef USE_MPI

#define BUFSIZE 16384
#define MIN_CHUNK_D_MULTIPLIER 10
#define MAX_CHUNK_SIZE 1000000

#include <string.h>

#include "mpi.h"
#include "mpifuncs.h"
#include "structs.h"
#include "funcs.h"

#define VERSION_STRING "INFERNAL 0.71"
#define VERSION_STRING_SIZE 100

/***************************************************************************
 * General functions, which can be used for either MPI search or alignment *
 ***************************************************************************/
/*
 * Function: get_master_rank
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Given a communicator, returns the lowest ranked one that can
 *           do I/O.
 *           Also checks the version string -- makes sure that it's the same
 *           on all procs.
 */
int get_master_rank (MPI_Comm comm, int mpi_my_rank) {
  int *io_proc_rank_p;
  int i;
  char versionbuf[VERSION_STRING_SIZE];
  
  MPI_Attr_get (comm, MPI_IO, &io_proc_rank_p, &i);
  
  if (i == 0)                 /* Not MPI compliant */
    return(MPI_PROC_NULL);
  if (*io_proc_rank_p == MPI_PROC_NULL)
    return (MPI_PROC_NULL);
  if (*io_proc_rank_p == MPI_ANY_SOURCE)
    return (0);

  /* Take min of procs allowed to do I/O.  */
  MPI_Allreduce (io_proc_rank_p, &i, 1, MPI_INT, MPI_MIN, comm);

  /* i is now master rank */
  /* Broadcast the version from master rank */
  if (i == mpi_my_rank)
    strncpy (versionbuf, VERSION_STRING, VERSION_STRING_SIZE-1);
  MPI_Bcast (versionbuf, VERSION_STRING_SIZE, MPI_CHAR, i, comm);
  
  if (strncmp (versionbuf, VERSION_STRING, VERSION_STRING_SIZE))
    Die ("Version strings %s and %s don't match\n", versionbuf, VERSION_STRING);
  return (i);
}

/**************************************************************
 * MPI search functions ("search_*") Originally from RSEARCH. *
 **************************************************************/
/*
 * Function: search_first_broadcast()
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Broadcasts the first set of parameters needed in all processes
 *
 * dmin and dmax will be NULL if QDB is not being used
 */
void search_first_broadcast (int *num_samples, int *windowlen, float *W_scale,
			     CM_t **cm, int *do_qdb, int **dmin, int **dmax, 
			     int *do_inside, int mpi_my_rank, 
			     int mpi_master_rank) {
  char buf[BUFSIZE];      /* Buffer for packing it all but the bulk of the CM */
  int position = 0;         /* Where I am in the buffer */
  int nstates, nnodes;
  float el_selfsc;
  int  iel_selfsc;
  int W;
  /*printf("entered first_broadcast: do_qdb: %d do_inside: %d my: %d master: %d\n", *do_qdb, *do_inside, 
    mpi_my_rank, mpi_master_rank);*/

  position = 0;
  if (mpi_my_rank == mpi_master_rank) 
    {   /* I'm in charge */
      nstates = (*cm)->M;
      nnodes = (*cm)->nodes;
      el_selfsc  = (*cm)->el_selfsc;
      iel_selfsc = (*cm)->iel_selfsc;
      W = (*cm)->W;
      
      /* Some ints */
      MPI_Pack (num_samples, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (windowlen, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      /* A float */
      MPI_Pack (W_scale, 1, MPI_FLOAT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      /* Basics of the model */
      MPI_Pack (&nstates, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD); 
      MPI_Pack (&nnodes, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->flags), 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&el_selfsc, 1, MPI_FLOAT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&iel_selfsc, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&W, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (do_qdb, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (do_inside, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    }
  /* Broadcast to everyone */
  MPI_Bcast (buf, BUFSIZE, MPI_PACKED, mpi_master_rank, MPI_COMM_WORLD);

  /* Decode this first set */
  position = 0;
  if (mpi_my_rank != mpi_master_rank) 
    {
      MPI_Unpack (buf, BUFSIZE, &position, num_samples, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, windowlen, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, W_scale, 1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &nstates, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &nnodes, 1, MPI_INT, MPI_COMM_WORLD);
      *cm = CreateCM (nnodes, nstates);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->flags), 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &el_selfsc, 1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &iel_selfsc, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &W, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, do_qdb, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, do_inside, 1, MPI_INT, MPI_COMM_WORLD);
      (*cm)->W = W;
      (*cm)->el_selfsc  =  el_selfsc;
      (*cm)->iel_selfsc = iel_selfsc;
    }
  /* Now we broadcast the rest of the model using many calls to MPI_Bcast.  
     This is inefficient, but is probably negligible compared to the actual 
     searches */
  MPI_Bcast ((*cm)->null, Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->sttype, nstates, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ndidx, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->stid, nstates, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->cfirst, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->cnum, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->plast, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->pnum, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
 
  MPI_Bcast ((*cm)->nodemap, nnodes, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ndtype, nnodes, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);

  MPI_Bcast ((*cm)->begin, nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->end, nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->beginsc, nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->endsc, nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ibeginsc, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->iendsc, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);

  /* if in QDB mode: broadcast dmin and dmax */
  if(*do_qdb) 
    {
      if(mpi_my_rank != mpi_master_rank)
	{
	  *dmin = MallocOrDie(sizeof(int) * nstates);
	  *dmax = MallocOrDie(sizeof(int) * nstates);
	}
      MPI_Bcast (*dmin, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
      MPI_Bcast (*dmax, nstates, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
    }

  /* These next calls depend on Sean's FMX2Alloc to be what CreateCM calls, and to allocate one large
     memory chunk at x[0] (where x is float **) and then fill in x[1]..x[n] with the appropriate offsets into
     this chunk of memory */
  MPI_Bcast ((*cm)->t[0], nstates*MAXCONNECT, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->e[0], nstates*Alphabet_size*Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->tsc[0], nstates*MAXCONNECT, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->esc[0], nstates*Alphabet_size*Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->itsc[0], nstates*MAXCONNECT, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->iesc[0], nstates*Alphabet_size*Alphabet_size, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);

  /*printf("leaving search_first_broadcast: do_qdb: %d do_inside: %d\n", (*do_qdb), (*do_inside));*/
}

/*
 * Function: search_second_broadcast
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Second broadcast of information to all processes.
 */
void search_second_broadcast (float *sc_cutoff, float *e_cutoff, int *cutoff_type, int *do_revcomp, 
			      int *do_align,
			      double *mu, double *lambda, double *K, long *N,
			      int mpi_my_rank, int mpi_master_rank) {
  char buf[BUFSIZE];      /* Buffer for packing it all but the bulk of the CM */
  int position = 0;         /* Where I am in the buffer */

  position = 0;
  if (mpi_my_rank == mpi_master_rank) {   /* I'm in charge */
    MPI_Pack (sc_cutoff, 1, MPI_FLOAT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (e_cutoff, 1, MPI_FLOAT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (cutoff_type, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (do_revcomp, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (do_align, 1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (N, 1, MPI_LONG, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (mu, GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (lambda, GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
    MPI_Pack (K, GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
  }
  MPI_Bcast (buf, BUFSIZE, MPI_PACKED, mpi_master_rank, MPI_COMM_WORLD);

  /* Decode this first set */
  position = 0;
  if (mpi_my_rank != mpi_master_rank) {
    MPI_Unpack (buf, BUFSIZE, &position, sc_cutoff, 1, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack (buf, BUFSIZE, &position, e_cutoff, 1, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack (buf, BUFSIZE, &position, cutoff_type, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, BUFSIZE, &position, do_revcomp, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, BUFSIZE, &position, do_align, 1, MPI_INT, MPI_COMM_WORLD) ;
    MPI_Unpack (buf, BUFSIZE, &position, N, 1, MPI_LONG, MPI_COMM_WORLD);
    MPI_Unpack (buf, BUFSIZE, &position, mu, GC_SEGMENTS, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack (buf, BUFSIZE, &position, lambda, GC_SEGMENTS, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack (buf, BUFSIZE, &position, K, GC_SEGMENTS, MPI_FLOAT, MPI_COMM_WORLD);
  }
}

/*
 * Function: search_receive_job
 * Date:     RJK, Mon May 28, 2002 [St. Louis]
 * Purpose:  Calls MPI_Recv to receive a work packet, and unpacks the packet
 */
char search_receive_job (int *seqlen_p, char **seq_p, int *bestr_p,
		  int mpi_master_rank) {
  char buf[32];
  MPI_Status status;
  int position = 0;
  char job_type;
  char *seq = NULL;

  /* Get the job type and the sequence length (0 if no sequence) */
  MPI_Recv (buf, 32, MPI_PACKED, mpi_master_rank, JOB_PACKET_TAG, 
	    MPI_COMM_WORLD, &status);
  MPI_Unpack (buf, 32, &position, &job_type, 1, MPI_CHAR, MPI_COMM_WORLD);
  MPI_Unpack (buf, 32, &position, seqlen_p, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack (buf, 32, &position, bestr_p, 1, MPI_INT, MPI_COMM_WORLD);

  if (*seqlen_p > 0) {
    /* Receive a partial sequence and convert to digitized sequence format
       by placing sentinels at end */
    seq = MallocOrDie(sizeof(char)*(*seqlen_p+2));
    seq[0] = seq[*seqlen_p+1] = DIGITAL_SENTINEL;
    MPI_Recv (seq+1, *seqlen_p, MPI_CHAR, mpi_master_rank, SEQ_TAG, 
	      MPI_COMM_WORLD, &status);
  }
  *seq_p = seq;
  return(job_type);
}

/*
 * Function: search_send_scan_results
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Given a list of scan results, sends them to master node using
 *           MPI_Send
 */
void search_send_scan_results (scan_results_t *results, int mpi_master_node) {
  char *buf;
  int pos = 0;
  int i;
  int bufsize;
  char results_type = SEARCH_STD_SCAN_RESULTS;

  bufsize = sizeof(char)+sizeof(int)+(sizeof(scan_result_node_t)*(results->num_results+1));
  buf = MallocOrDie(bufsize);
  /* Send the size of the results */
  MPI_Send (&bufsize, 1, MPI_INT, mpi_master_node, SEARCH_STD_SCAN_RESULTS_SIZE_TAG, MPI_COMM_WORLD);

  /* Send the results */
  MPI_Pack (&results_type, 1, MPI_CHAR, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (&(results->num_results), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  for (i=0; i<results->num_results; i++) {
    MPI_Pack(&(results->data[i].start), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
    MPI_Pack(&(results->data[i].stop), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
    MPI_Pack(&(results->data[i].bestr), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
    MPI_Pack(&(results->data[i].score), 1, MPI_FLOAT, buf, bufsize, &pos, MPI_COMM_WORLD);
  }
  MPI_Send (buf, bufsize, MPI_PACKED, mpi_master_node, SEARCH_STD_SCAN_RESULTS_TAG, MPI_COMM_WORLD);
}

/*
 * Function: search_send_align_results
 * Date:     RJK, Tue Jun 25, 2002 [St. Louis]
 * Purpose:  Sends a parse tree back to the master process 
 */
void search_send_align_results (Parsetree_t *tr, int mpi_master_node) {
  char *buf;
  int pos = 0;
  int bufsize;
  char results_type = ALN_RESULTS;

  bufsize = sizeof(char)+sizeof(int)+sizeof(int)+6*((tr->n)+1)*sizeof(int);
  buf = MallocOrDie(bufsize);
  /* Send the size of the results */
  MPI_Send (&bufsize, 1, MPI_INT, mpi_master_node, SEARCH_STD_SCAN_RESULTS_SIZE_TAG, MPI_COMM_WORLD);

  /* Send the results */
  MPI_Pack (&results_type, 1, MPI_CHAR, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (&(tr->memblock), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (&(tr->n), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (tr->emitl, tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (tr->emitr, tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (tr->state, tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (tr->nxtl, tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (tr->nxtr, tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (tr->prv, tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);

  MPI_Send (buf, bufsize, MPI_PACKED, mpi_master_node, SEARCH_STD_SCAN_RESULTS_TAG, MPI_COMM_WORLD);
}

/*
 * Function: search_procs_working
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 * Purpose:  Returns TRUE if a slave is still working (status not NULL), false
 *           otherwise.
 */
int search_procs_working (job_t **process_status, int mpi_num_procs, int mpi_master_rank) { 
  int i;
  for (i=0; i<mpi_num_procs; i++) {
    if (i == mpi_master_rank) continue;
    if (process_status[i] != NULL) 
      return(TRUE);
  }
  return (FALSE);
}

/*
 * Function: search_enqueue
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 * Purpose:  Chunks the given database sequence and enqueues all chunks for
 *           scanning.  Queue is a linked list.
 */
job_t *search_enqueue (db_seq_t *active_seq, int db_seq_index, 
		int D, int do_revcomp, int job_type) 
{
  job_t *queue = NULL;
  job_t *cur_tail = NULL;
  job_t *new_entry;
  int chunksize;
  int chunkoffset;
  int curpos;
  int in_revcomp;
  int mpi_num_procs;

  MPI_Comm_size (MPI_COMM_WORLD, &mpi_num_procs);

  /* 
   * Set the chunk size as follows:
   * 1.  Ideally take smallest multiple of D that gives result greater than
   *     (seqlen+D*(num_procs-2))/(num_procs-1)
   *     This should put one chunk on each processor.
   * 2.  If this is less than MIN_CHUNK_D_MULTIPLIER * D, use that value.
   * 3.  If this is greater than MAX_CHUNK_SIZE, use that.
   */
  chunksize = ((active_seq->sq[0]->n+D*(mpi_num_procs-2))/(mpi_num_procs-1))+1;
  chunksize = ((chunksize/D)+1)*D;

  /*printf("seq len: %d chunksize: %d procs: %d\n", active_seq->sq[0]->n, chunksize, mpi_num_procs);*/

  if (do_revcomp)
    chunksize *= 2;
  if (chunksize < MIN_CHUNK_D_MULTIPLIER * D) 
    chunksize = MIN_CHUNK_D_MULTIPLIER * D;
  if (chunksize > MAX_CHUNK_SIZE)
    chunksize = MAX_CHUNK_SIZE;

  chunkoffset = chunksize - D;
  active_seq->chunks_sent = 0;
  active_seq->alignments_sent = -1;     /* None sent yet */
  for (in_revcomp = 0; in_revcomp <= do_revcomp; in_revcomp++) {
    for (curpos = 1; curpos <= active_seq->sq[0]->n; curpos += chunkoffset) {
      new_entry = MallocOrDie (sizeof(db_seq_t));
      new_entry->next = NULL;
      if (cur_tail != NULL)
	cur_tail->next = new_entry;
      if (queue == NULL)
	queue = new_entry;
      cur_tail = new_entry;
      new_entry -> job_type = job_type;
      if (chunksize < active_seq->sq[0]->n - curpos + 1)
	new_entry -> seqlen = chunksize;
      else
	new_entry -> seqlen = active_seq->sq[0]->n - curpos + 1;
      new_entry -> dsq = active_seq->sq[in_revcomp]->dsq + curpos;
      new_entry -> db_seq_index = db_seq_index;
      new_entry->in_revcomp = (char)in_revcomp;
      new_entry->index = curpos;
      new_entry->bestr = 0;          /* Not used here */
      /*printf("new_entry len: %d i: %d j: %d\n", new_entry->seqlen, curpos, (curpos+ new_entry->seqlen));*/
      active_seq->chunks_sent++;
    }
  }
  return (queue);
}

/*
 * Function: search_enqueue_alignments
 * Date:     RJK, Tue Jun 25, 2002 [St. Louis]
 * Purpose:  Given a sequence and a job type, takes each result for sequence
 *           and enqueues an alignment.
 */
void search_enqueue_alignments (job_t **queue, db_seq_t *active_seq, int db_seq_index,
			 int do_revcomp, int job_type) {
  job_t *cur_tail = NULL;
  job_t *new_entry;
  int result_index;
  int in_revcomp;
  
  /* First, set cur tail to end if already stuff in queue */
  if (*queue != NULL) {
    cur_tail = *queue;
    while (cur_tail->next != NULL) {
      cur_tail = cur_tail -> next;
    }
  }

  active_seq->alignments_sent = 0;

  for (in_revcomp = 0; in_revcomp <= do_revcomp; in_revcomp++) {
    if (active_seq->results[in_revcomp] == NULL) continue;    /* No results here */
    for (result_index = 0; 
	 result_index < active_seq->results[in_revcomp]->num_results;
	 result_index++) {
      new_entry = MallocOrDie (sizeof(db_seq_t));
      new_entry->next = NULL;
      if (cur_tail != NULL)
	cur_tail->next = new_entry;
      if (*queue == NULL)
	*queue = new_entry;
      cur_tail = new_entry;
      new_entry -> job_type = job_type;
      new_entry->seqlen = 
	active_seq->results[in_revcomp]->data[result_index].stop - 
	active_seq->results[in_revcomp]->data[result_index].start + 1;
      new_entry -> dsq = active_seq->sq[in_revcomp]->dsq + 
	active_seq->results[in_revcomp]->data[result_index].start;
      new_entry -> db_seq_index = db_seq_index;
      new_entry->in_revcomp = (char)in_revcomp;
      new_entry->index = result_index;
      new_entry->bestr = 
	active_seq->results[in_revcomp]->data[result_index].bestr;
      active_seq->alignments_sent++;
    }
  }
}

/*
 * Function: search_send_next_job
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 * Purpose:  Takes first job in queue, sends it to the given rank, and updates
 *           queue and process_status to reflect this.
 */
void search_send_next_job (job_t **queue, job_t **process_status, int rank_to_send_to) {
  char buf[32];
  int position = 0;

  /* Set this process as having the job */
  *process_status = *queue;

  /* Remove job from queue */
  *queue = (*queue)->next;

  MPI_Pack (&((*process_status)->job_type), 1, MPI_CHAR, buf, 32, &position, MPI_COMM_WORLD);
  MPI_Pack (&((*process_status)->seqlen), 1, MPI_INT, buf, 32, &position, MPI_COMM_WORLD);
  MPI_Pack (&((*process_status)->bestr), 1, MPI_INT, buf, 32, &position, MPI_COMM_WORLD);
  MPI_Send (buf, position, MPI_PACKED, rank_to_send_to, JOB_PACKET_TAG, MPI_COMM_WORLD);

  if ((*process_status)->seqlen > 0) {
    MPI_Send ((*process_status)->dsq, (*process_status)->seqlen, MPI_CHAR, rank_to_send_to, SEQ_TAG, MPI_COMM_WORLD);
  }
}

/*
 * Function: search_check_results
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 * Purpose:  Does a blocking call to MPI_Recv until a process finishes, then
 *           processes results and returns.
 */
int search_check_results (db_seq_t **active_seqs, job_t **process_status, int D) {
  char *buf;
  int bufsize;
  MPI_Status status;
  int data_from;        /* Who's sending us data */
  char results_type;
  int position = 0;
  int num_results;
  int start, stop, bestr;
  float score;
  db_seq_t *cur_seq;
  char in_revcomp;
  int index;
  int i;
  int cur_seq_index;
  Parsetree_t *tr;

  /* Get the size of the buffer */
  MPI_Recv (&bufsize, 1, MPI_INT, MPI_ANY_SOURCE, 
	    SEARCH_STD_SCAN_RESULTS_SIZE_TAG, MPI_COMM_WORLD, &status);
  data_from = status.MPI_SOURCE;
  buf = MallocOrDie(sizeof(char)*bufsize);

  /* Figure out the sequence it belongs to */
  cur_seq_index = process_status[data_from]->db_seq_index;
  cur_seq = active_seqs[cur_seq_index];
  index = process_status[data_from]->index;
  in_revcomp = process_status[data_from]->in_revcomp;

  /* Clear this job -- it's done */
  free(process_status[data_from]);
  process_status[data_from] = NULL;

  /* Now get the results */
  MPI_Recv (buf, bufsize, MPI_PACKED, data_from, SEARCH_STD_SCAN_RESULTS_TAG, 
	    MPI_COMM_WORLD, &status);
  MPI_Unpack (buf, bufsize, &position, &results_type, 1, 
	      MPI_CHAR, MPI_COMM_WORLD);

  if (results_type == SEARCH_STD_SCAN_RESULTS) {
    MPI_Unpack (buf, bufsize, &position, &num_results, 1, 
		MPI_INT, MPI_COMM_WORLD);
    if (num_results > 0 && cur_seq->results[(int)in_revcomp] == NULL) {
      cur_seq->results[(int)in_revcomp] = CreateResults (INIT_RESULTS);
    }
    for (i=0; i<num_results; i++) {
      MPI_Unpack(buf, bufsize, &position, &start, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buf, bufsize, &position, &stop, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buf, bufsize, &position, &bestr, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buf, bufsize, &position, &score,1, MPI_FLOAT, MPI_COMM_WORLD);
      /* Don't report hits from first D nucleotides in second overlapping seq
	 because it wasn't a full analysis for seqs ending there -- longer
	 seqs missed */
      if (index == 1 || stop > D)
	report_hit (index+start-1, index+stop-1, bestr, score, 
		    cur_seq->results[(int)in_revcomp]);
    }
    cur_seq->chunks_sent--;
  } else if (results_type == ALN_RESULTS) {
    tr = MallocOrDie(sizeof(Parsetree_t));
    /* Get size of the tree */
    MPI_Unpack (buf, bufsize, &position, &(tr->memblock), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, &(tr->n), 1, MPI_INT, MPI_COMM_WORLD);
    /* Allocate it */
    tr->emitl = MallocOrDie(sizeof(int)*tr->n);
    tr->emitr = MallocOrDie(sizeof(int)*tr->n);
    tr->state = MallocOrDie(sizeof(int)*tr->n);
    tr->nxtl = MallocOrDie(sizeof(int)*tr->n);
    tr->nxtr = MallocOrDie(sizeof(int)*tr->n);
    tr->prv = MallocOrDie(sizeof(int)*tr->n);
    tr->nalloc = tr->n;
    /* Unpack it */
    MPI_Unpack (buf, bufsize, &position, tr->emitl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->emitr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->state, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->prv, tr->n, MPI_INT, MPI_COMM_WORLD);
    cur_seq->results[(int)in_revcomp]->data[index].tr = tr;
    cur_seq->alignments_sent--;
  } else {
    Die ("Got result type %d when expecting SEARCH_STD_SCAN_RESULTS (%d) or ALN_RESULTS (%d)\n", results_type, SEARCH_STD_SCAN_RESULTS, ALN_RESULTS);
  }
  free(buf);
  return(cur_seq_index);
}

/*
 * Function: search_check_hist_results
 * Date:     RJK, Sat Jun 01, 2002 [St. Louis]
 * Purpose:  Does a blocking call to MPI_Recv until a histogram scan process 
 *           finishes, then processes results and returns.  Based on
 *           check_results.
 */
int search_check_hist_results (db_seq_t **seqs, job_t **process_status, int D) {

  char buf[32];
  MPI_Status status;
  int data_from;        /* Who's sending us data */
  char results_type;
  int position = 0;
  int db_seq_index;
  float score;

  /* Get the hist results -- score to add and whether this score was
     masked by better score or not */
  MPI_Recv (buf, 32, MPI_PACKED, MPI_ANY_SOURCE, SEARCH_HIST_RESULTS_TAG,
	    MPI_COMM_WORLD, &status);
  data_from = status.MPI_SOURCE;

  /* Figure out the sequence it belongs to */
  db_seq_index = process_status[data_from]->db_seq_index;

  /* Clear this job -- it's done */
  free(process_status[data_from]);
  process_status[data_from] = NULL;


  MPI_Unpack (buf, 32, &position, &results_type, 1, 
	      MPI_CHAR, MPI_COMM_WORLD);
  if (results_type == SEARCH_HIST_SCAN_RESULTS) {
    MPI_Unpack (buf, 32, &position, &score, 1, MPI_FLOAT, 
		MPI_COMM_WORLD);
    seqs[db_seq_index]->chunks_sent--;
    if (seqs[db_seq_index]->best_score < score) {
      seqs[db_seq_index]->best_score = score;
    }
  } else {
    Die ("Got result type %d when expecting SEARCH_HIST_SCAN_RESULTS (%d)\n", 
	 results_type, SEARCH_HIST_SCAN_RESULTS);
  }
  return(db_seq_index);
}

/*
 * Function: search_send_hist_scan_results
 * Date:     RJK, Sat Jun 01, 2002 [St. Louis]
 * Purpose:  Given best score for building a histogram, sends to master
 *           node using MPI_Send
 */
void search_send_hist_scan_results (float score, int mpi_master_node) {
  char *buf;
  int pos = 0;
  int bufsize;
  char results_type = SEARCH_HIST_SCAN_RESULTS;

  bufsize = sizeof(char)+sizeof(float);
  buf = MallocOrDie(bufsize);

  /* Send the results */
  MPI_Pack (&results_type, 1, MPI_CHAR, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (&score, 1, MPI_FLOAT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Send (buf, bufsize, MPI_PACKED, mpi_master_node, SEARCH_HIST_RESULTS_TAG, MPI_COMM_WORLD);
}


/*
 * Function: search_send_terminate
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 * Purpose:  Sends termination code to specified process 
 */
void search_send_terminate (int i) {
  char buf[16];
  int position = 0;
  char job_type = TERMINATE_WORK;
  int len = 0;

  MPI_Pack (&job_type, 1, MPI_CHAR, buf, 16, &position, MPI_COMM_WORLD);
  MPI_Pack (&len, 1, MPI_INT, buf, 16, &position, MPI_COMM_WORLD);
  MPI_Send (buf, position, MPI_PACKED, i, JOB_PACKET_TAG, MPI_COMM_WORLD);
}

/**************************************************************
 * MPI alignment functions ("aln_*") (EPN 01.04.07)
 **************************************************************/
/*
 * Function: aln_broadcast()
 * Date:     EPN, Thu Jan  4 14:33:07 2007
 * Purpose:  Broadcasts the CM for alignment (in cmalign).
 *
 */
void aln_broadcast (CM_t **cm, int mpi_my_rank, int mpi_master_rank) 
{
  char buf[BUFSIZE];      /* Buffer for packing it all but the bulk of the CM */
  int position = 0;         /* Where I am in the buffer */
  int nstates, nnodes;
  int enf_len;
  /*printf("entered aln_broadcast my: %d master: %d\n", mpi_my_rank, mpi_master_rank);*/
  
  position = 0;
  if (mpi_my_rank == mpi_master_rank) 
    {   /* I'm in charge */
      nstates = (*cm)->M;
      nnodes = (*cm)->nodes;
      
      /* Basics of the model */
      MPI_Pack (&nstates,              1, MPI_INT,   buf, BUFSIZE, &position, MPI_COMM_WORLD); 
      MPI_Pack (&nnodes,               1, MPI_INT,   buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->flags),       1, MPI_INT,   buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->align_flags), 1, MPI_INT,   buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->el_selfsc),   1, MPI_FLOAT, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->iel_selfsc),  1, MPI_INT,   buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->W),           1, MPI_INT,   buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->enf_start),   1, MPI_INT,   buf, BUFSIZE, &position, MPI_COMM_WORLD);
      if((*cm)->enf_start != 0)
	{
	  enf_len = strlen((*cm)->enf_seq);
	  MPI_Pack (&enf_len,          1,       MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
	  MPI_Pack (&((*cm)->enf_seq), enf_len, MPI_CHAR,   buf, BUFSIZE, &position, MPI_COMM_WORLD);
	}
    }
  /* Broadcast to everyone */
  MPI_Bcast (buf, BUFSIZE, MPI_PACKED, mpi_master_rank, MPI_COMM_WORLD);

  /* Decode this first set */
  position = 0;
  if (mpi_my_rank != mpi_master_rank) 
    {
      MPI_Unpack (buf, BUFSIZE, &position, &nstates, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &nnodes, 1, MPI_INT, MPI_COMM_WORLD);
      *cm = CreateCM (nnodes, nstates);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->flags),       1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->align_flags), 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->el_selfsc),   1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->iel_selfsc),  1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->W),           1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_start),   1, MPI_INT, MPI_COMM_WORLD);
      if((*cm)->enf_start != 0)
	{
	  MPI_Unpack (buf, BUFSIZE, &position, &enf_len,                  1, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_seq),   enf_len, MPI_INT, MPI_COMM_WORLD);
	}
    }
  /* Now we broadcast the rest of the model using many calls to MPI_Bcast.  
     This is inefficient, but is probably negligible compared to the actual 
     searches */
  MPI_Bcast ((*cm)->null,   Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->sttype, nstates,       MPI_CHAR,  mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ndidx,  nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->stid,   nstates,       MPI_CHAR,  mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->cfirst, nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->cnum,   nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->plast,  nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->pnum,   nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
 
  MPI_Bcast ((*cm)->nodemap, nnodes, MPI_INT,  mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ndtype,  nnodes, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);

  MPI_Bcast ((*cm)->begin,    nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->end,      nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->beginsc,  nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->endsc,    nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ibeginsc, nstates, MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->iendsc,   nstates, MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);

  /* These next calls depend on Sean's FMX2Alloc to be what CreateCM calls, and to allocate one large
     memory chunk at x[0] (where x is float **) and then fill in x[1]..x[n] with the appropriate offsets into
     this chunk of memory */
  MPI_Bcast ((*cm)->t[0], nstates*MAXCONNECT, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->e[0], nstates*Alphabet_size*Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->tsc[0], nstates*MAXCONNECT, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->esc[0], nstates*Alphabet_size*Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->itsc[0], nstates*MAXCONNECT, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->iesc[0], nstates*Alphabet_size*Alphabet_size, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  /*printf("leaving aln_broadcast\n");*/
  return;
}


/*
 * Function: aln_procs_working
 * Date:     EPN, Sat Dec 30 21:38:54 2006
 * Purpose:  Returns TRUE if a slave is still working (status BUSY), false
 *           otherwise.
 */
int aln_procs_working (int *process_status, int mpi_num_procs, int mpi_master_rank) { 
  int i;
  for (i=0; i<mpi_num_procs; i++) 
    {
      if (i == mpi_master_rank) continue;
      if (process_status[i] == BUSY) return(TRUE); 
    }
  return (FALSE);
}

/*
 * Function: aln_send_next_job
 * Date:     EPN, Mon Jan  1 10:37:31 2007
 * Purpose:  Given a seqs_to_aln_t data structure with sequences to send to 
 *           a slave, send the important bits to the slave.
 * Args:
 *         seqs_to_aln - the seqs to send to the slave for alignment
 */
void aln_send_next_job (seqs_to_aln_t *seqs_to_aln, int rank_to_send_to) 
{
  char *buf;
  int bufsize;
  int position = 0;
  char job_type = ALN_WORK;
  int i;
  int *namelen;
  /*printf("in aln_send_next_job, rank_to_send_to: %d\n", rank_to_send_to);*/

  bufsize = sizeof(char) + ((3 + (2 * seqs_to_aln->nseq)) * sizeof(int));
  buf = MallocOrDie(bufsize);
  
  /* Send the size of the job */
  MPI_Send (&bufsize, 1, MPI_INT, rank_to_send_to, ALN_JOB_SIZE_TAG, MPI_COMM_WORLD);

  MPI_Pack (&job_type, 1, MPI_CHAR, buf, bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack (&(seqs_to_aln->index), 1, MPI_INT, buf, bufsize, &position, MPI_COMM_WORLD);

  namelen = MallocOrDie(sizeof(int) * seqs_to_aln->nseq);

  /* Send the seqs_to_aln_t data structure, piecewise */
  /* First send the number of sequences */
  MPI_Pack (&(seqs_to_aln->nseq), 1, MPI_INT, buf, bufsize, &position, MPI_COMM_WORLD);
  /* For each sequence, send the length of each sequence, and the length of each sequence name */
  for(i = 0; i < seqs_to_aln->nseq; i++)
    {
      namelen[i] = strlen(seqs_to_aln->sq[i]->name)+1;
      MPI_Pack (&(seqs_to_aln->sq[i]->n), 1, MPI_INT, buf, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack (&(namelen[i]), 1, MPI_INT, buf, bufsize, &position, MPI_COMM_WORLD);
    }
  MPI_Send (buf, position, MPI_PACKED, rank_to_send_to, ALN_JOB_PACKET_TAG, MPI_COMM_WORLD);

  /* Send the names of the sequences */
  for(i = 0; i < seqs_to_aln->nseq; i++)
    if(seqs_to_aln->sq[i]->n > 0)
      MPI_Send (seqs_to_aln->sq[i]->name, (namelen[i]), MPI_CHAR, rank_to_send_to, 
		SEQ_TAG, MPI_COMM_WORLD);

  /* send the digitized seqs only, (not the full ESL_SQ data structure) */
  for(i = 0; i < seqs_to_aln->nseq; i++)
    {
      if(seqs_to_aln->sq[i]->n > 0)
	MPI_Send (seqs_to_aln->sq[i]->dsq, (seqs_to_aln->sq[i]->n+2), MPI_CHAR, rank_to_send_to, SEQ_TAG, MPI_COMM_WORLD);
    }

  /*printf("leaving aln_send_next_job, rank_to_send_to: %d\n", rank_to_send_to);*/
  free(namelen);
}

/*
 * Function: aln_receive_job
 * Date:     EPN, Sat Dec 30 22:13:24 2006
 * Purpose:  Calls MPI_Recv to receive a aln work packet, and unpacks the packet
 */
char aln_receive_job (seqs_to_aln_t **ret_seqs_to_aln, int mpi_master_rank) 
{
  char *buf;
  int bufsize;
  MPI_Status status;
  int position = 0;
  char job_type;
  seqs_to_aln_t *seqs_to_aln;
  seqs_to_aln = MallocOrDie(sizeof(seqs_to_aln_t));
  int i;
  int *namelen;

  /* Get the size of the buffer */
  MPI_Recv (&bufsize, 1, MPI_INT, MPI_ANY_SOURCE, ALN_JOB_SIZE_TAG, MPI_COMM_WORLD, &status);
  buf = MallocOrDie(sizeof(char)*bufsize);

  printf("buf size: %d\n", bufsize);
  
  /* Get the job */
  MPI_Recv (buf, bufsize, MPI_PACKED, mpi_master_rank, ALN_JOB_PACKET_TAG, 
	    MPI_COMM_WORLD, &status);
  MPI_Unpack (buf, bufsize, &position, &job_type, 1, MPI_CHAR, MPI_COMM_WORLD);

  if(job_type == ALN_WORK)
    {
      MPI_Unpack (buf, bufsize, &position, &(seqs_to_aln->index), 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, bufsize, &position, &(seqs_to_aln->nseq),  1, MPI_INT, MPI_COMM_WORLD);
      seqs_to_aln->sq = MallocOrDie(sizeof(ESL_SQ *) * seqs_to_aln->nseq);
      namelen         = MallocOrDie(sizeof(int)      * seqs_to_aln->nseq);

      for(i = 0; i < seqs_to_aln->nseq; i++)
	{
	  seqs_to_aln->sq[i] = esl_sq_Create();
	  MPI_Unpack (buf, bufsize, &position, &(seqs_to_aln->sq[i]->n), 1, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, &(namelen[i]), 1, MPI_INT, MPI_COMM_WORLD);
	}
      
      /* Receive the names of the sequences */
      for(i = 0; i < seqs_to_aln->nseq; i++)
	{
	  if(seqs_to_aln->sq[i]->n > 0)
	    MPI_Recv (seqs_to_aln->sq[i]->name, namelen[i], MPI_CHAR, 
		      mpi_master_rank, SEQ_TAG, MPI_COMM_WORLD, &status);
	}

      /* Receive the digitized sequences */
      for(i = 0; i < seqs_to_aln->nseq; i++)
	{
	  if(seqs_to_aln->sq[i]->n > 0)
	    {
	      seqs_to_aln->sq[i]->dsq = MallocOrDie(sizeof(char)*((seqs_to_aln->sq[i]->n)+2));
	      MPI_Recv (seqs_to_aln->sq[i]->dsq, (seqs_to_aln->sq[i]->n+2), MPI_CHAR, 
			mpi_master_rank, SEQ_TAG, MPI_COMM_WORLD, &status);
	    }
	}
      *ret_seqs_to_aln = seqs_to_aln;
    }
  else if(job_type != TERMINATE_WORK)
    Die("ERROR in aln_receive_job did not receive ALN_WORK or TERMINATE_WORK signal\n");

  return(job_type);
}

/*
 * Function: aln_send_results
 * Date:     EPN, Sat Dec 30 22:05:30 2006
 * Purpose:  Sends parse tree(s) and potentially postcodes, back to the master process. 
 */
void aln_send_results (seqs_to_aln_t *seqs_to_aln, int do_post, int mpi_master_node) 
{
  char *buf;
  int pos = 0;
  int bufsize;
  int i;
  char results_type = ALN_RESULTS;

  bufsize = sizeof(char) + sizeof(int) + sizeof(int) + sizeof(int);
  for(i = 0; i < seqs_to_aln->nseq; i++)
    bufsize += (6 * (seqs_to_aln->tr[i]->n + 1) * sizeof(int));

  if(do_post) /* add size of postcodes */
    for(i = 0; i < seqs_to_aln->nseq; i++)
      bufsize += (seqs_to_aln->sq[i]->n + 1) * sizeof(char) + sizeof(int);

  buf = MallocOrDie(bufsize);

  /* Send the size of the results */
  MPI_Send (&bufsize, 1, MPI_INT, mpi_master_node, ALN_RESULTS_SIZE_TAG, MPI_COMM_WORLD);

  /* Send the results */
  MPI_Pack (&results_type,         1, MPI_CHAR, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (&(seqs_to_aln->index), 1, MPI_INT , buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (&(seqs_to_aln->nseq),  1, MPI_INT , buf, bufsize, &pos, MPI_COMM_WORLD);
  for(i = 0; i < seqs_to_aln->nseq; i++)
    {
      MPI_Pack (&(seqs_to_aln->tr[i]->memblock),1,                     MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
      MPI_Pack (&(seqs_to_aln->tr[i]->n),       1,                     MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
      MPI_Pack (seqs_to_aln->tr[i]->emitl,      seqs_to_aln->tr[i]->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
      MPI_Pack (seqs_to_aln->tr[i]->emitr,      seqs_to_aln->tr[i]->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
      MPI_Pack (seqs_to_aln->tr[i]->state,      seqs_to_aln->tr[i]->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
      MPI_Pack (seqs_to_aln->tr[i]->nxtl,       seqs_to_aln->tr[i]->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
      MPI_Pack (seqs_to_aln->tr[i]->nxtr,       seqs_to_aln->tr[i]->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
      MPI_Pack (seqs_to_aln->tr[i]->prv,        seqs_to_aln->tr[i]->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
    }
  MPI_Pack (&do_post,  1, MPI_INT , buf, bufsize, &pos, MPI_COMM_WORLD);

  if(do_post)
    for(i = 0; i < seqs_to_aln->nseq; i++)
      {
	MPI_Pack (&(seqs_to_aln->sq[i]->n), 1                        , MPI_INT,  buf, bufsize, &pos, MPI_COMM_WORLD);
	MPI_Pack (seqs_to_aln->postcode[i],   (seqs_to_aln->sq[i]->n+1), MPI_CHAR, buf, bufsize, &pos, MPI_COMM_WORLD);
      }
  MPI_Send (buf, bufsize, MPI_PACKED, mpi_master_node, ALN_RESULTS_TAG, MPI_COMM_WORLD);
}

/*
 * Function: aln_check_results
 * Date:     EPN, Sat Dec 30 21:49:42 2006
 * Purpose:  Does a blocking call to MPI_Recv until a process finishes, then
 *           processes results and returns.
 */
int aln_check_results (Parsetree_t **all_parsetrees, char **all_postcodes, int **process_status)
{
  char *buf;
  int bufsize;
  MPI_Status status;
  int data_from;        /* Who's sending us data */
  char results_type;
  int position = 0;
  int i;
  Parsetree_t *tr;
  int nseq;
  int index;
  int postsize;
  int do_post;
  char *postcode;

  /*printf("in aln_check_results\n");*/

  /* Get the size of the buffer */
  MPI_Recv (&bufsize, 1, MPI_INT, MPI_ANY_SOURCE, 
	    ALN_RESULTS_SIZE_TAG, MPI_COMM_WORLD, &status);
  data_from = status.MPI_SOURCE;
  buf = MallocOrDie(sizeof(char)*bufsize);

  /* Clear this job -- it's done */
  (*process_status)[data_from] = IDLE;

  /* Now get the results */
  MPI_Recv (buf, bufsize, MPI_PACKED, data_from, ALN_RESULTS_TAG, 
	    MPI_COMM_WORLD, &status);
  MPI_Unpack (buf, bufsize, &position, &results_type, 1, 
	      MPI_CHAR, MPI_COMM_WORLD);

  if (results_type == ALN_RESULTS) 
    {
      MPI_Unpack (buf, bufsize, &position, &index, 1, 
		  MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, bufsize, &position, &nseq, 1, 
		  MPI_INT, MPI_COMM_WORLD);
      for (i=0; i < nseq; i++) 
	{
	  /* Get a parsetree for each sequence */
	  tr = MallocOrDie(sizeof(Parsetree_t));
	  /* Get size of the tree */
	  MPI_Unpack (buf, bufsize, &position, &(tr->memblock), 1, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, &(tr->n), 1, MPI_INT, MPI_COMM_WORLD);
	  /* Allocate it */
	  tr->emitl = MallocOrDie(sizeof(int)*tr->n);
	  tr->emitr = MallocOrDie(sizeof(int)*tr->n);
	  tr->state = MallocOrDie(sizeof(int)*tr->n);
	  tr->nxtl  = MallocOrDie(sizeof(int)*tr->n);
	  tr->nxtr  = MallocOrDie(sizeof(int)*tr->n);
	  tr->prv   = MallocOrDie(sizeof(int)*tr->n);
	  tr->nalloc = tr->n;
	  /* Unpack it */
	  MPI_Unpack (buf, bufsize, &position, tr->emitl, tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->emitr, tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->state, tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->nxtl,  tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->nxtr,  tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->prv,   tr->n, MPI_INT, MPI_COMM_WORLD);
	  /* add the parsetree onto the master array of parsetrees */
	  all_parsetrees[index++] = tr;
	}

      MPI_Unpack (buf, bufsize, &position, &do_post,  1, MPI_INT, MPI_COMM_WORLD);

      if(do_post)
	{
	  index -= nseq;
	  for (i=0; i < nseq; i++) 
	    {
	      /* unpack the postcodes */
	      MPI_Unpack (buf, bufsize, &position, &postsize,   1       , MPI_INT,  MPI_COMM_WORLD);
	      postcode = MallocOrDie(sizeof(char) * (postsize+1));
	      MPI_Unpack (buf, bufsize, &position, postcode, (postsize+1), MPI_CHAR, MPI_COMM_WORLD);
	      all_postcodes[index++] = postcode;
	    }
	}
    }
  else 
    Die ("Got result type %d when expecting ALN_RESULTS (%d)\n", results_type, ALN_RESULTS);
      
  free(buf);
  return(data_from);
}

/*
 * Function: aln_send_terminate
 * Date:     EPN, Thu Jan  4 16:16:45 2007
 * Purpose:  Sends termination code to specified process, specific for cmalign (not cmsearch).
 */
void aln_send_terminate (int rank_to_send_to) 
{
  int bufsize;
  char *buf;
  int position = 0;
  char job_type = TERMINATE_WORK;

  bufsize = sizeof(char);
  buf = MallocOrDie(bufsize);
  
  /* Send the bufsize (only b/c aln_receive_job expects it) */
  MPI_Send (&bufsize, 1, MPI_INT, rank_to_send_to, ALN_JOB_SIZE_TAG, MPI_COMM_WORLD);

  MPI_Pack (&job_type, 1, MPI_CHAR, buf, bufsize, &position, MPI_COMM_WORLD);
  MPI_Send (buf, position, MPI_PACKED, rank_to_send_to, ALN_JOB_PACKET_TAG, MPI_COMM_WORLD);
}


#endif



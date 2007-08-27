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

#ifdef HAVE_MPI

#include "esl_config.h"
#include "config.h"

#include <string.h>

#include "mpi.h"
#include "mpifuncs.h"
#include "structs.h"
#include "funcs.h"
#include "cm_dispatch.h"
#include "stats.h"

#if 0

#define BUFSIZE 16384
#define MIN_CHUNK_D_MULTIPLIER 10
#define MAX_CHUNK_SIZE 1000000

#define VERSION_STRING "INFERNAL 1.0"
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

/*
 * Function: broadcast_cm()
 * Date:     EPN, Thu Jan  4 14:33:07 2007
 * Purpose:  Broadcasts the CM for alignment or search.
 *
 */
void broadcast_cm (CM_t **cm, int mpi_my_rank, int mpi_master_rank) 
{
  char buf[BUFSIZE];      /* Buffer for packing it all but the bulk of the CM */
  int position = 0;         /* Where I am in the buffer */
  int nstates, nnodes;
  int enf_len;
  int nparts;
  int i;
  int p;

  position = 0;
  if (mpi_my_rank == mpi_master_rank) 
    {   /* I'm in charge */
      /* contract check, if we claim to have Gumbel stats, we better have them */
      if((*cm)->flags & CM_GUMBEL_STATS && (*cm)->stats == NULL)
	esl_fatal("ERROR in broadcast_cm() master node claims to have Gumbel stats but cm->stats is NULL!\n");
      nstates = (*cm)->M;
      nnodes = (*cm)->nodes;
      
      /* Basics of the model */
      MPI_Pack (&nstates,                  1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD); 
      MPI_Pack (&nnodes,                   1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->flags),           1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->config_opts),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->align_opts),      1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->search_opts),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->el_selfsc),       1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->iel_selfsc),      1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->W),               1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->enf_start),       1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->enf_scdiff),      1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->sc_boost),        1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_sc_boost),    1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->ffract),          1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cutoff_type),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cutoff),          1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_cutoff_type), 1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_cutoff),      1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->beta),            1, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->tau),             1, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->hmmpad),          1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->pbegin),          1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->pend),            1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);

      /* Take special care with enf_len, this is used later to get cm->enf_seq if nec */
      if((*cm)->enf_start != 0) enf_len = strlen((*cm)->enf_seq);
      else enf_len = 0;
      MPI_Pack (&enf_len,                  1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);

      /* Take special care with number of partitions, used later to get cm->stats if nec */
      if((*cm)->flags & CM_GUMBEL_STATS) nparts = (*cm)->stats->np;
      else nparts = 0;
      MPI_Pack (&nparts,                  1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);

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
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->flags),           1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->config_opts),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->align_opts),      1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->search_opts),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->el_selfsc),       1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->iel_selfsc),      1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->W),               1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_start),       1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_scdiff),      1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->sc_boost),        1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_sc_boost),    1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->ffract),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cutoff_type),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cutoff),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_cutoff_type), 1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_cutoff),      1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->beta),            1, MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->tau),             1, MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->hmmpad),          1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->pbegin),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->pend),            1, MPI_FLOAT,   MPI_COMM_WORLD);

      MPI_Unpack (buf, BUFSIZE, &position, &enf_len,                  1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &nparts,                   1, MPI_INT,   MPI_COMM_WORLD);
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

  /* Broadcast the enf_seq, if it's NULL (enf_start == 0) we don't */
  if((*cm)->enf_start != 0)
    {
      if (mpi_my_rank != mpi_master_rank) 
	ESL_ALLOC((*cm)->enf_seq, sizeof(char) * (enf_len+1));
      MPI_Bcast((*cm)->enf_seq, enf_len, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);
      if (mpi_my_rank != mpi_master_rank) 
	(*cm)->enf_seq[enf_len] = '\0';
    }

  /* Broadcast the Gumbel stats if they exist 
   * IMPT: currently filter threshold stats are NOT broadcasted as they're only
   * used to get cm->cp9_cutoff, which is broadcasted separately. We could get 
   * away with not broadcasting these stats too - though we'd have to modify 
   * parallel_search_database() to be independent on Gumbel params */
  if((*cm)->flags & CM_GUMBEL_STATS) /* flags were already sent/received */
    {
      if (mpi_my_rank != mpi_master_rank) 
	(*cm)->stats = AllocCMStats(nparts); /* nparts was already sent/recd */
      for(i = 0; i < NGUMBELMODES; i++)
	for(p = 0; p < nparts; p++)
	  {	    
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->N),      1, MPI_INT,    mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->L),      1, MPI_INT,    mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->mu),     1, MPI_DOUBLE, mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->lambda), 1, MPI_DOUBLE, mpi_master_rank, MPI_COMM_WORLD);
	  }
    }
  return;
}

/**************************************************************
 * MPI search functions ("search_*") Originally from RSEARCH. *
 **************************************************************/
/*
 * Function: search_first_broadcast()
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Broadcasts only num_samples and W_scale for parallel cmsearch.
 *
 */
void search_first_broadcast (int *num_samples, float *W_scale,
			     int mpi_my_rank, int mpi_master_rank) 
{
  int  status;
  int   bufsize = sizeof(int) + sizeof(float);
  char *buf;                /* Buffer for packing */
  int   position = 0;         /* Where I am in the buffer */
  ESL_ALLOC(buf, bufsize);

  /*printf("entered search_first_broadcast: my: %d master: %d\n", mpi_my_rank, mpi_master_rank);*/

  position = 0;
  if (mpi_my_rank == mpi_master_rank) 
    {   /* I'm in charge */
      MPI_Pack (num_samples, 1, MPI_INT,   buf, bufsize, &position, MPI_COMM_WORLD);
      MPI_Pack (W_scale,     1, MPI_FLOAT, buf, bufsize, &position, MPI_COMM_WORLD);
    }
  /* Broadcast to everyone */
  MPI_Bcast (buf, bufsize, MPI_PACKED, mpi_master_rank, MPI_COMM_WORLD);

  /* Decode it */
  position = 0;
  if (mpi_my_rank != mpi_master_rank) 
    {
      MPI_Unpack (buf, BUFSIZE, &position, num_samples, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, W_scale, 1, MPI_FLOAT, MPI_COMM_WORLD);
    }
  /*printf("leaving search_first_broadcast: do_qdb: %d do_inside: %d\n", (*do_qdb), (*do_inside));*/

  return;
 ERROR:
  esl_fatal("Memory allocation error.");
}

/*
 * Function: search_second_broadcast
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Second broadcast of information to all processes.
 */
void search_second_broadcast (CM_t **cm, long *N, int mpi_my_rank, int mpi_master_rank) 
{
  char buf[BUFSIZE];      /* Buffer for packing it all but the bulk of the CM */
  int position = 0;         /* Where I am in the buffer */
  double lambda[GC_SEGMENTS];
  double mu[GC_SEGMENTS];
  double K[GC_SEGMENTS];
  double cp9_lambda[GC_SEGMENTS];
  double cp9_mu[GC_SEGMENTS];
  double cp9_K[GC_SEGMENTS];
  int     i;

  position = 0;
  if (mpi_my_rank == mpi_master_rank)    /* I'm in charge */
    { 
      /* we always send N */
      MPI_Pack (&N, 1, MPI_LONG, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      if((*cm)->search_opts & CM_SEARCH_CMSTATS) /* pack the CM Gumbel parameters */
	{
	  //MPI_Pack ((*cm)->lambda,GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
	  //MPI_Pack ((*cm)->K,     GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
	  //MPI_Pack ((*cm)->mu,    GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
	}
      if((*cm)->search_opts & CM_SEARCH_CP9STATS) /* pack the CP9 Gumbel parameters */
	{
	  //MPI_Pack ((*cm)->cp9_lambda, GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
	  //MPI_Pack ((*cm)->cp9_mu,     GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
	  //MPI_Pack ((*cm)->cp9_K,      GC_SEGMENTS, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
	}
    }
  MPI_Bcast (buf, BUFSIZE, MPI_PACKED, mpi_master_rank, MPI_COMM_WORLD);
  
  /* Decode the packet */
  position = 0;
  if (mpi_my_rank != mpi_master_rank) 
    {
      MPI_Unpack (buf, BUFSIZE, &position, N, 1, MPI_LONG, MPI_COMM_WORLD);
      /* this is fragile, we rely on the fact that we broadcasted the CM
       * in broadcast_cm() earlier, so it will have the same cm->*opts as
       * the one we're sending from the master node. */
      if((*cm)->search_opts & CM_SEARCH_CMSTATS) /* unpack the CM Gumbel parameters */
	{
	  //MPI_Unpack (buf, BUFSIZE, &position, ((*cm)->lambda), GC_SEGMENTS, MPI_DOUBLE, MPI_COMM_WORLD);
	  //MPI_Unpack (buf, BUFSIZE, &position, ((*cm)->mu),     GC_SEGMENTS, MPI_DOUBLE, MPI_COMM_WORLD);
	  //MPI_Unpack (buf, BUFSIZE, &position, ((*cm)->K),      GC_SEGMENTS, MPI_DOUBLE, MPI_COMM_WORLD);
	  for(i = 0; i < GC_SEGMENTS; i++)
	    {
	      //(*cm)->lambda[i] = lambda[i];
	      //(*cm)->mu[i]     = mu[i];
	      //(*cm)->K[i]      = K[i];
	    }
	}
      if((*cm)->search_opts & CM_SEARCH_CP9STATS) /* unpack the CP9 Gumbel parameters */
	{
	  //MPI_Unpack (buf, BUFSIZE, &position, cp9_lambda, GC_SEGMENTS, MPI_DOUBLE, MPI_COMM_WORLD);
	  //MPI_Unpack (buf, BUFSIZE, &position, cp9_mu,     GC_SEGMENTS, MPI_DOUBLE, MPI_COMM_WORLD);
	  //MPI_Unpack (buf, BUFSIZE, &position, cp9_K,      GC_SEGMENTS, MPI_DOUBLE, MPI_COMM_WORLD);
	  for(i = 0; i < GC_SEGMENTS; i++)
	    {
	      //(*cm)->cp9_lambda[i] = cp9_lambda[i];
	      //(*cm)->cp9_mu[i]     = cp9_mu[i];
	      //(*cm)->cp9_K[i]      = cp9_K[i];
	    }
	}
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

  /*  printf("in receive job seqlen: %d\n", *seqlen_p);*/
  if (*seqlen_p > 0) {
    /* Receive a partial sequence and convert to digitized sequence format
       by placing sentinels at end */
    ESL_ALLOC(seq = MallocOrDie(sizeof(char)*(*seqlen_p+2)));
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
  ESL_ALLOC(buf, bufsize);
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

  bufsize = sizeof(char)+sizeof(int)+sizeof(int)+7*((tr->n)+1)*sizeof(int);
  ESL_ALLOC(buf, bufsize);
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
  MPI_Pack (tr->prv,  tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
  MPI_Pack (tr->mode, tr->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);

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

  /*printf("in search_enqueue, D: %d do_revcomp: %d\n", D, do_revcomp);
    printf("seq len: %d chunksize: %d procs: %d\n", active_seq->sq[0]->n, chunksize, mpi_num_procs);*/

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
      ESL_ALLOC(new_entry, sizeof(db_seq_t));
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
      ESL_ALLOC(new_entry, sizeof(db_seq_t));
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

  /*printf("in search_send_next_job sent len: %d\n", (*process_status)->seqlen);*/
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
  ESL_ALLOC(buf, (char)*bufsize);

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
    ESL_ALLOC(tr, sizeof(Parsetree_t));
    /* Get size of the tree */
    MPI_Unpack (buf, bufsize, &position, &(tr->memblock), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, &(tr->n), 1, MPI_INT, MPI_COMM_WORLD);
    /* Allocate it */
    ESL_ALLOC(tr->emitl, sizeof(int)*tr->n);
    ESL_ALLOC(tr->emitr, sizeof(int)*tr->n);
    ESL_ALLOC(tr->state, sizeof(int)*tr->n);
    ESL_ALLOC(tr->nxtl,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->nxtr,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->prv,   sizeof(int)*tr->n);
    ESL_ALLOC(tr->mode,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->nalloc,sizeof(int)*tr->n);
    tr->nalloc = tr->n;

    /* Unpack it */
    MPI_Unpack (buf, bufsize, &position, tr->emitl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->emitr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->state, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->prv, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->mode, tr->n, MPI_INT, MPI_COMM_WORLD);
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
  ESL_ALLOC(buf, bufsize);

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
  ESL_ALLOC(buf, bufsize);
  
  /* Send the size of the job */
  MPI_Send (&bufsize, 1, MPI_INT, rank_to_send_to, ALN_JOB_SIZE_TAG, MPI_COMM_WORLD);

  MPI_Pack (&job_type, 1, MPI_CHAR, buf, bufsize, &position, MPI_COMM_WORLD);
  MPI_Pack (&(seqs_to_aln->index), 1, MPI_INT, buf, bufsize, &position, MPI_COMM_WORLD);

  ESL_ALLOC(namelen, sizeof(int) * seqs_to_aln->nseq);

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
  ESL_ALLOC(seqs_to_aln, sizeof(seqs_to_aln_t));
  int i;
  int *namelen;

  /* Get the size of the buffer */
  MPI_Recv (&bufsize, 1, MPI_INT, MPI_ANY_SOURCE, ALN_JOB_SIZE_TAG, MPI_COMM_WORLD, &status);
  ESL_ALLOC(buf, sizeof(char)*bufsize);

  /* Get the job */
  MPI_Recv (buf, bufsize, MPI_PACKED, mpi_master_rank, ALN_JOB_PACKET_TAG, 
	    MPI_COMM_WORLD, &status);
  MPI_Unpack (buf, bufsize, &position, &job_type, 1, MPI_CHAR, MPI_COMM_WORLD);

  if(job_type == ALN_WORK)
    {
      MPI_Unpack (buf, bufsize, &position, &(seqs_to_aln->index), 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, bufsize, &position, &(seqs_to_aln->nseq),  1, MPI_INT, MPI_COMM_WORLD);
      ESL_ALLOC(seqs_to_aln->sq, sizeof(ESL_SQ *) * seqs_to_aln->nseq);
      ESL_ALLOC(namelen,         sizeof(int)      * seqs_to_aln->nseq);

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
	      ESL_ALLOC(seqs_to_aln->sq[i]->dsq, sizeof(char)*((seqs_to_aln->sq[i]->n)+2));
	      MPI_Recv (seqs_to_aln->sq[i]->dsq, (seqs_to_aln->sq[i]->n+2), MPI_CHAR, 
			mpi_master_rank, SEQ_TAG, MPI_COMM_WORLD, &status);
	    }
	}
      *ret_seqs_to_aln = seqs_to_aln;
    }
  else if(job_type != TERMINATE_WORK)
    esl_fatal("ERROR in aln_receive_job did not receive ALN_WORK or TERMINATE_WORK signal\n");

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
    bufsize += (7 * (seqs_to_aln->tr[i]->n + 1) * sizeof(int));

  if(do_post) /* add size of postcodes */
    for(i = 0; i < seqs_to_aln->nseq; i++)
      bufsize += (seqs_to_aln->sq[i]->n + 1) * sizeof(char) + sizeof(int);

  ESL_ALLOC(buf, bufsize);

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
      MPI_Pack (seqs_to_aln->tr[i]->mode,       seqs_to_aln->tr[i]->n, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
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
  ESL_ALLOC(buf, sizeof(char)*bufsize);

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
	  ESL_ALLOC(tr, sizeof(Parsetree_t));
	  /* Get size of the tree */
	  MPI_Unpack (buf, bufsize, &position, &(tr->memblock), 1, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, &(tr->n), 1, MPI_INT, MPI_COMM_WORLD);
	  /* Allocate it */
	  ESL_ALLOC(tr->emitl, sizeof(int)*tr->n);
	  ESL_ALLOC(tr->emitr, sizeof(int)*tr->n);
	  ESL_ALLOC(tr->state, sizeof(int)*tr->n);
	  ESL_ALLOC(tr->nxtl,  sizeof(int)*tr->n);
	  ESL_ALLOC(tr->nxtr,  sizeof(int)*tr->n);
	  ESL_ALLOC(tr->prv,   sizeof(int)*tr->n);
	  ESL_ALLOC(tr->mode,  sizeof(int)*tr->n);
	  ESL_ALLOC(tr->nalloc,sizeof(int)*tr->n);
	  tr->nalloc = tr->n;

	  /* Unpack it */
	  MPI_Unpack (buf, bufsize, &position, tr->emitl, tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->emitr, tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->state, tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->nxtl,  tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->nxtr,  tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->prv,   tr->n, MPI_INT, MPI_COMM_WORLD);
	  MPI_Unpack (buf, bufsize, &position, tr->mode,  tr->n, MPI_INT, MPI_COMM_WORLD);
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
	      ESL_ALLOC(postcode, sizeof(char) * (postsize+1));
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
  ESL_ALLOC(buf, bufsize);
  
  /* Send the bufsize (only b/c aln_receive_job expects it) */
  MPI_Send (&bufsize, 1, MPI_INT, rank_to_send_to, ALN_JOB_SIZE_TAG, MPI_COMM_WORLD);

  MPI_Pack (&job_type, 1, MPI_CHAR, buf, bufsize, &position, MPI_COMM_WORLD);
  MPI_Send (buf, position, MPI_PACKED, rank_to_send_to, ALN_JOB_PACKET_TAG, MPI_COMM_WORLD);
}

/**************************************************************************************/
/* EPN, Thu May 10 10:11:18 2007 New functions roughly following Easel/H3 conventions */
/* Function: mpi_worker_search_target()
 * Incept:   EPN, Wed May  9 17:07:48 2007
 * Purpose:  The main control for an MPI worker process for searching sequences. 
 *           Worker receives CM, then loops over receipt of sequences, returning
 *           best score and results data structure for each.
 *           Never do revcomp, we'll call this function twice once with 
 *           plus once with minus strand.
 */
void
mpi_worker_search_target(CM_t *cm, int my_rank)
{
  int status;
  char *dsq = NULL;
  int   L;
  float best_sc;

  int doing_cm_stats  = FALSE;
  int doing_cp9_stats = FALSE;

  int ctr = 0;
  if(cm->search_opts & CM_SEARCH_HMMONLY) doing_cp9_stats = TRUE;
  else doing_cm_stats = TRUE;
  /* Main loop */
  while (dsq_MPIRecv(&dsq, &L) == eslOK)
    {
      best_sc = actually_search_target(cm, dsq, 1, L, 
				       0.,    /* minimum CM bit cutoff, irrelevant (?) */
				       0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
				       NULL,  /* do not keep results */
				       FALSE, /* do not filter with a CP9 HMM */
				       doing_cm_stats, doing_cp9_stats,
				       NULL); /* filter fraction, nobody cares */
      
      MPI_Send(&(best_sc), 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
      free(dsq);
    }
  return;

 ERROR:
  if (dsq != NULL) free(dsq);
  return;
}

/* Function:  dsq_MPISend()
 * Incept:    EPN, Wed May  9 17:30:14 2007
 *
 * Purpose:   Send sequence <dsq> to processor <dest>.
 *            
 *            If <dsq> is NULL, sends a end-of-data signal to <dest>, to
 *            tell it to shut down.
 */
int
dsq_MPISend(char *dsq, int L, int dest)
{
  if(dsq == NULL) 
    {
      int eod = -1;
      MPI_Send(&eod, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      return eslOK;
    }
  MPI_Send(&(L), 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
  /* receiver will now allocate storage, before reading on...*/
  MPI_Send(dsq, (L+2), MPI_CHAR, dest, 0, MPI_COMM_WORLD);
  return eslOK;
}

/* Function:  dsq_MPISend()
 * Incept:    EPN, Thu Jun  7 15:02:34 2007
 *
 * Purpose:   Send sequence <dsq> and max bit sc <maxsc> to processor <dest>.
 *            
 *            If <dsq> is NULL, sends a end-of-data signal to <dest>, to
 *            tell it to shut down.
 */
int
dsq_maxsc_MPISend(char *dsq, int L, float maxsc, int dest)
{
  if(dsq == NULL) 
    {
      int eod = -1;
      MPI_Send(&eod, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      return eslOK;
    }
  MPI_Send(&(L), 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
  /* receiver will now allocate storage, before reading on...*/
  MPI_Send(dsq, (L+2), MPI_CHAR, dest, 0, MPI_COMM_WORLD);
  MPI_Send(&(maxsc), 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
  return eslOK;
}

/* Function:  dsq_MPIRecv()
 * Incept:    EPN, Wed May  9 17:34:43 2007
 *
 * Purpose:   Receive a sequence sent from the master MPI process (src=0)
 *            on a worker MPI process. 
 *            
 *            If it receives an end-of-data signal, returns <eslEOD>.
 */
int
dsq_MPIRecv(char **ret_dsq, int *ret_L)
{
  int status;
  char *dsq = NULL;
  MPI_Status mpistatus;
  int L;
  
  MPI_Recv(&L, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  if (L == -1) return eslEOD;
  ESL_ALLOC(dsq, sizeof(char) * (L+2));
  MPI_Recv(dsq, (L+2), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &mpistatus);
  *ret_L   = L;
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  dsq_maxsc_MPIRecv()
 * Incept:    EPN, Thu Jun  7 15:00:29 2007    
 *
 * Purpose:   Receive a sequence and maximum score 
 *            sent from the master MPI process (src=0)
 *            on a worker MPI process. 
 *            
 *            If it receives an end-of-data signal, returns <eslEOD>.
 */
int
dsq_maxsc_MPIRecv(char **ret_dsq, int *ret_L, float *ret_maxsc)
{
  int status;
  char *dsq = NULL;
  MPI_Status mpistatus;
  int L;
  float maxsc;
  
  MPI_Recv(&L, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  if (L == -1) return eslEOD;
  ESL_ALLOC(dsq, sizeof(char) * (L+2));
  MPI_Recv(dsq, (L+2), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(&maxsc, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  *ret_L   = L;
  *ret_dsq = dsq;
  *ret_maxsc = maxsc;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_MPIBroadcast()
 * Incept:    EPN, Wed May  9 17:24:53 2007
 *
 * Purpose:   Sends CM <cm> to processor <dest>.
 *            
 *            If <cm> is NULL, sends a end-of-data signal to <dest>, to
 *            tell it to shut down.
 *            
 */
int
cm_MPIBroadcast(CM_t *cm)
{
  return eslOK;
}  

/* Function: mpi_worker_cm_and_cp9_search()
 * Incept:   EPN, Thu May 10 10:04:02 2007
 * Purpose:  The main control for an MPI worker process for searching sequences
 *           twice, once with a CM and once with a CP9, both scores are returned.
 *           Called in mpi_FindCP9FilterThreshold9).
 * Args:
 *           cm       - the covariance model
 *           do_fast  - don't search with CM, only do CP9 search
 *           my_rank  - my MPI rank
 */
void
mpi_worker_cm_and_cp9_search(CM_t *cm, int do_fast, int my_rank)
{
  int status;
  char *dsq = NULL;
  int   L;
  float *scores = NULL;
  ESL_ALLOC(scores, sizeof(float) * 2);
  int was_hmmonly;
  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  /* Main loop */
  while (dsq_MPIRecv(&dsq, &L) == eslOK)
    {
      /* Do the CM search first */
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      if(do_fast)
	scores[0] = 0.;
      else
	scores[0] = actually_search_target(cm, dsq, 1, L, 
					   0.,    /* minimum CM bit cutoff, irrelevant (?) */
					   0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
					   NULL,  /* do not keep results */
					   FALSE, /* do not filter with a CP9 HMM */
					   FALSE, FALSE, /* not doing CM or CP9 Gumbel calcs */
					   NULL); /* filter fraction, nobody cares */
      /* DO NOT CALL actually_search_target b/c that will run Forward then 
       * Backward to get score of best hit, but we'll be detecting by a
       * Forward scan (then running Backward only on hits above our threshold),
       * since we're calc'ing the threshold here it's impt we only do Forward.
       */
      scores[1] =  CP9Forward(cm, dsq, 1, L, cm->W, 0., 
			      NULL,   /* don't return scores of hits */
			      NULL,   /* don't return posns of hits */
			      NULL,   /* don't keep track of hits */
			      TRUE,   /* we're scanning */
			      FALSE,  /* we're not ultimately aligning */
			      FALSE,  /* we're not rescanning */
			      TRUE,   /* be memory efficient */
			      NULL);  /* don't want the DP matrix back */
      
      MPI_Send(scores, 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); /* send together so results don't interleave */
      free(dsq);
    }
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;
  if(scores != NULL) free(scores);
  /*printf("\trank: %d RETURNING!\n", my_rank);*/
  return;

 ERROR:
  if (dsq != NULL) free(dsq);
  if (scores != NULL) free(scores);
  return;
}

/* Function: mpi_worker_cm_and_cp9_search_maxsc()
 * Incept:   EPN, Thu Jun  7 15:02:54 2007   
 * Purpose:  The main control for an MPI worker process for searching sequences
 *           with decreasingly fast techniques, quitting if any technique 
 *           returns a score greater than some specified maximum bit score. 
 *           The goal is to determine if the optimal parse is within a 
 *           given range during a empirical HMM filter threshold calculation. 
 *           Called in mpi_FindCP9FilterThreshold().
 * Args:
 *           cm       - the covariance model
 *           do_fast  - don't search with CM, only do CP9 search
 *           my_rank  - my MPI rank
 */
void
mpi_worker_cm_and_cp9_search_maxsc(CM_t *cm, int do_fast, int do_minmax, int my_rank)
{
  int status;
  char *dsq = NULL;
  int   L;
  float maxsc;
  float *scores = NULL;
  ESL_ALLOC(scores, sizeof(float) * 2);
  int was_hmmonly;
  int was_hbanded;
  float orig_tau;
  int passed_flag = FALSE;
  orig_tau = cm->tau;

  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  if(cm->search_opts & CM_SEARCH_HBANDED) was_hbanded = TRUE;
  else was_hbanded = FALSE;
  /* Main loop */
  while (dsq_maxsc_MPIRecv(&dsq, &L, &maxsc) == eslOK)
    {
      /* Do the CM search first */
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      if(do_fast)
	scores[0] = 0.;
      else if(do_minmax)
	{
	  cm->search_opts |= CM_SEARCH_HBANDED;
	  cm->tau = 0.1;
	  scores[0] = actually_search_target(cm, dsq, 1, L,
					 0.,    /* cutoff is 0 bits (actually we'll find highest
						 * negative score if it's < 0.0) */
					 0.,    /* CP9 cutoff is 0 bits */
					 NULL,  /* don't keep results */
					 FALSE, /* don't filter with a CP9 HMM */
					 FALSE, /* we're not calcing CM  stats */
					 FALSE, /* we're not calcing CP9 stats */
					 NULL); /* filter fraction N/A */
	  if(scores[0] < maxsc) /* search with another, less strict (lower tau)  HMM banded parse */
	    {
	      cm->tau = 1e-10;
	      scores[0] = actually_search_target(cm, dsq, 1, L,
					     0.,    /* cutoff is 0 bits (actually we'll find highest
						     * negative score if it's < 0.0) */
					     0.,    /* CP9 cutoff is 0 bits */
					     NULL,  /* don't keep results */
					     FALSE, /* don't filter with a CP9 HMM */
					     FALSE, /* we're not calcing CM  stats */
					     FALSE, /* we're not calcing CP9 stats */
					     NULL); /* filter fraction N/A */
	    }
	}
      else
	scores[0] = actually_search_target(cm, dsq, 1, L, 
					   0.,    /* minimum CM bit cutoff, irrelevant (?) */
					   0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
					   NULL,  /* do not keep results */
					   FALSE, /* do not filter with a CP9 HMM */
					   FALSE, FALSE, /* not doing CM or CP9 Gumbel calcs */
					   NULL); /* filter fraction, nobody cares */

      /* Now do HMM search, but if do_minmax, only do HMM search 
       * if our CM score hasn't exceeded the max */
      if(do_minmax && scores[0] > maxsc)
	scores[1] = 0.;
      else
	/* DO NOT CALL actually_search_target b/c that will run Forward then 
	 * Backward to get score of best hit, but we'll be detecting by a
	 * Forward scan (then running Backward only on hits above our threshold),
	 * since we're calc'ing the threshold here it's impt we only do Forward.
	 */
	scores[1] =  CP9Forward(cm, dsq, 1, L, cm->W, 0., 
				NULL,   /* don't return scores of hits */
				NULL,   /* don't return posns of hits */
				NULL,   /* don't keep track of hits */
				TRUE,   /* we're scanning */
				FALSE,  /* we're not ultimately aligning */
				FALSE,  /* we're not rescanning */
				TRUE,   /* be memory efficient */
				NULL);  /* don't want the DP matrix back */
      MPI_Send(scores, 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); /* send together so results don't interleave */
      free(dsq);
    }
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;
  if(was_hbanded) cm->search_opts |= CM_SEARCH_HBANDED;
  else cm->search_opts &= ~CM_SEARCH_HBANDED;

  if(scores != NULL) free(scores);
  /*printf("\trank: %d RETURNING!\n", my_rank);*/
  return;

 ERROR:
  if (dsq != NULL) free(dsq);
  if (scores != NULL) free(scores);
  return;
}
#endif

/* Function:  cm_MPIBroadcast()
 * Incept:    EPN, Wed May  9 17:24:53 2007
 *
 * Purpose:   Sends CM <cm> to processor <dest>.
 *            
 *            If <cm> is NULL, sends a end-of-data signal to <dest>, to
 *            tell it to shut down.
 *            
 */
int
cm_MPIBroadcast(CM_t *cm)
{
  int   status;
  int   code;
  int   sz, n, pos;

  /* Figure out size */
  if (MPI_Pack_size(1, MPI_INT, comm, &n) != 0) ESL_XEXCEPTION(eslESYS, "mpi pack size failed");
  if (cm != NULL) {
    if ((status = cm_MPIPackSize(cm, comm, &sz)) != eslOK) return status;
    n += sz;
  }

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Pack the status code and HMM into the buffer */
  pos  = 0;
  code = (hmm == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &pos, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed");
  if (hmm != NULL) {
    if ((status = p7_hmm_MPIPack(hmm, *buf, n, &pos, comm)) != eslOK) return status;
  }

  /* Send the packed HMM to the destination. */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0)  ESL_EXCEPTION(eslESYS, "mpi send failed");
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_MPIPackSize()
 * Synopsis:  Calculates size needed to pack a CM.
 * Incept:    EPN, Mon Aug 27 10:34:15 2007
 *            based on p7_hmm_MPIPackSize() from HMMER3.
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_MPIPack()> will need to pack a CM
 *            <cm> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
cm_MPIPackSize(CM_t *cm, MPI_Comm comm, int *ret_n)
{
  int   status;
  int   n = 0;
  int   K = cm->abc->K;
  int   M = cm->M;
  int   nnodes = cm->nodes;
  int   sz;

  if (MPI_Pack_size(1,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   
  n += 15*sz; 
  /* M, nodes, flags, config_opts, search_opts, align_opts, nseq, clen, iel_selfsc, W (10)
   * enf_start, cutoff_type, cp9_cutoff_type, hmmpad, np (5) 
   */

  if (MPI_Pack_size(1,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += (12+K)*sz; 
  /* el_selfsc, enf_scdiff, sc_boost, cp9_sc_boost, cutoff, cp9_cutoff, pbegin, pend, eff_nseq (9) 
   * ga, tc, nc, null (3+K) */

  if (MPI_Pack_size(1,      MPI_DOUBLE, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* beta, tau */

  if (MPI_Pack_size(M,        MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* sttype, stid */

  if (MPI_Pack_size(M,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 7*sz; 
  /* ndidx, cfirst, cnum, plast, pnum, ibeginsc, iendsc */

  if (MPI_Pack_size(M,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 4*sz; 
  /* begin, end, beginsc, endsc */

  if (MPI_Pack_size(nnodes,    MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* nodemap, ndtype */

  if (MPI_Pack_size(M*K*K,   MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* e, esc */

  if (MPI_Pack_size(M*MAXCONNECT,MPI_FLOAT,comm,&sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* t, tsc */

  if (MPI_Pack_size(M*K*K,    MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* iesc */

  if (MPI_Pack_size(M*MAXCONNECT,MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* itsc */

  if ((status = esl_mpi_PackOptSize(cm->enfseq, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; 
  n += sz; /* enfseq */

  if ((status = esl_mpi_PackOptSize(cm->name, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; 
  n += sz; /* name */

  if (cm->flags & CMH_ACC)  { if ((status = esl_mpi_PackOptSize(cm->acc, -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  if (cm->flags & CMH_DESC) { if ((status = esl_mpi_PackOptSize(cm->desc,-1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  
  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;

}

/* Function:  cm_justread_MPIPackSize()
 *
 * Synopsis:  Calculates size needed to pack a CM that has
 *            just been read from a CM file by a CMFileRead()
 *            call, we'll need to pack far less than a fully
 *            configure CM in this case.
 *
 * Incept:    EPN, Mon Aug 27 10:34:15 2007
 *            based on p7_hmm_MPIPackSize() from HMMER3.
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_MPIPack()> will need to pack a CM
 *            <cm> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
cm_justread_MPIPackSize(CM_t *cm, MPI_Comm comm, int *ret_n)
{
  int   status;
  int   n = 0;
  int   K = cm->abc->K;
  int   M = cm->M;
  int   nnodes = cm->nodes;
  int   sz;

  if (MPI_Pack_size(1,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   
  n += 4*sz; /* M, nodes, nseq, clen */ 

  if (MPI_Pack_size(1,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += (5+K)*sz; 
  /* el_selfsc, eff_nseq, ga, tc, nc, null (5+K) */

  if (MPI_Pack_size(M,        MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* sttype, stid */

  if (MPI_Pack_size(M,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 7*sz; 
  /* ndidx, cfirst, cnum, plast, pnum, ibeginsc, iendsc */

  if (MPI_Pack_size(nnodes,    MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* nodemap, ndtype */

  if (MPI_Pack_size(M*K*K,   MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* e */

  if (MPI_Pack_size(M*MAXCONNECT,MPI_FLOAT,comm,&sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* t */

  if ((status = esl_mpi_PackOptSize(cm->name, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; 
  n += sz; /* name */

  if (cm->flags & CMH_ACC)  { if ((status = esl_mpi_PackOptSize(cm->acc, -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  if (cm->flags & CMH_DESC) { if ((status = esl_mpi_PackOptSize(cm->desc,-1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  
  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;

}
#endif
void broadcast_cm (CM_t **cm, int mpi_my_rank, int mpi_master_rank) 
{
  char buf[BUFSIZE];      /* Buffer for packing it all but the bulk of the CM */
  int position = 0;         /* Where I am in the buffer */
  int nstates, nnodes;
  int enf_len;
  int nparts;
  int i;
  int p;

  position = 0;
  if (mpi_my_rank == mpi_master_rank) 
    {   /* I'm in charge */
      /* contract check, if we claim to have Gumbel stats, we better have them */
      if((*cm)->flags & CM_GUMBEL_STATS && (*cm)->stats == NULL)
	esl_fatal("ERROR in broadcast_cm() master node claims to have Gumbel stats but cm->stats is NULL!\n");
      nstates = (*cm)->M;
      nnodes = (*cm)->nodes;
      
      /* Basics of the model */
      MPI_Pack (&nstates,                  1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD); 
      MPI_Pack (&nnodes,                   1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->flags),           1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->config_opts),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->align_opts),      1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->search_opts),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->el_selfsc),       1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->iel_selfsc),      1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->W),               1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->enf_start),       1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->enf_scdiff),      1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->sc_boost),        1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_sc_boost),    1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cutoff_type),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cutoff),          1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_cutoff_type), 1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_cutoff),      1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->beta),            1, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->tau),             1, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->hmmpad),          1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->pbegin),          1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->pend),            1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);

      /* Take special care with enf_len, this is used later to get cm->enf_seq if nec */
      if((*cm)->enf_start != 0) enf_len = strlen((*cm)->enf_seq);
      else enf_len = 0;
      MPI_Pack (&enf_len,                  1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);

      /* Take special care with number of partitions, used later to get cm->stats if nec */
      if((*cm)->flags & CM_GUMBEL_STATS) nparts = (*cm)->stats->np;
      else nparts = 0;
      MPI_Pack (&nparts,                  1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);

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
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->flags),           1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->config_opts),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->align_opts),      1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->search_opts),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->el_selfsc),       1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->iel_selfsc),      1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->W),               1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_start),       1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_scdiff),      1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->sc_boost),        1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_sc_boost),    1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->ffract),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cutoff_type),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cutoff),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_cutoff_type), 1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_cutoff),      1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->beta),            1, MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->tau),             1, MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->hmmpad),          1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->pbegin),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->pend),            1, MPI_FLOAT,   MPI_COMM_WORLD);

      MPI_Unpack (buf, BUFSIZE, &position, &enf_len,                  1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &nparts,                   1, MPI_INT,   MPI_COMM_WORLD);
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

  /* Broadcast the enf_seq, if it's NULL (enf_start == 0) we don't */
  if((*cm)->enf_start != 0)
    {
      if (mpi_my_rank != mpi_master_rank) 
	ESL_ALLOC((*cm)->enf_seq, sizeof(char) * (enf_len+1));
      MPI_Bcast((*cm)->enf_seq, enf_len, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);
      if (mpi_my_rank != mpi_master_rank) 
	(*cm)->enf_seq[enf_len] = '\0';
    }

  /* Broadcast the Gumbel stats if they exist 
   * IMPT: currently filter threshold stats are NOT broadcasted as they're only
   * used to get cm->cp9_cutoff, which is broadcasted separately. We could get 
   * away with not broadcasting these stats too - though we'd have to modify 
   * parallel_search_database() to be independent on Gumbel params */
  if((*cm)->flags & CM_GUMBEL_STATS) /* flags were already sent/received */
    {
      if (mpi_my_rank != mpi_master_rank) 
	(*cm)->stats = AllocCMStats(nparts); /* nparts was already sent/recd */
      for(i = 0; i < NGUMBELMODES; i++)
	for(p = 0; p < nparts; p++)
	  {	    
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->N),      1, MPI_INT,    mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->L),      1, MPI_INT,    mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->mu),     1, MPI_DOUBLE, mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->lambda), 1, MPI_DOUBLE, mpi_master_rank, MPI_COMM_WORLD);
	  }
    }
  return eslOK;
}  


/* BEGIN NEW 1.0 CODE HEREHEREHEREHERE */
/*
 * Function: search_results_MPIUnpack() 
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 * Purpose:  Does a blocking call to MPI_Recv until a process finishes, then
 *           processes results and returns.
 */
#if 0
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
  ESL_ALLOC(buf, (char)*bufsize);

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
    ESL_ALLOC(tr, sizeof(Parsetree_t));
    /* Get size of the tree */
    MPI_Unpack (buf, bufsize, &position, &(tr->memblock), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, &(tr->n), 1, MPI_INT, MPI_COMM_WORLD);
    /* Allocate it */
    ESL_ALLOC(tr->emitl, sizeof(int)*tr->n);
    ESL_ALLOC(tr->emitr, sizeof(int)*tr->n);
    ESL_ALLOC(tr->state, sizeof(int)*tr->n);
    ESL_ALLOC(tr->nxtl,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->nxtr,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->prv,   sizeof(int)*tr->n);
    ESL_ALLOC(tr->mode,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->nalloc,sizeof(int)*tr->n);
    tr->nalloc = tr->n;

    /* Unpack it */
    MPI_Unpack (buf, bufsize, &position, tr->emitl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->emitr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->state, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->prv, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->mode, tr->n, MPI_INT, MPI_COMM_WORLD);
    cur_seq->results[(int)in_revcomp]->data[index].tr = tr;
    cur_seq->alignments_sent--;
  } else {
    Die ("Got result type %d when expecting SEARCH_STD_SCAN_RESULTS (%d) or ALN_RESULTS (%d)\n", results_type, SEARCH_STD_SCAN_RESULTS, ALN_RESULTS);
  }
  free(buf);
  return(cur_seq_index);
}

#endif


#endif

/*
 * stats.h
 * 
 * Header file for stats.c
 */

#ifndef _stats_h
#define _stats_h
#include "esl_sqio.h"

extern void serial_make_histogram (int *gc_count, int *partitions, int num_partitions,
				   CM_t *cm, int D, int num_samples, 
				   int sample_length, double *lambda, double *K,
				   int *dmin, int *dmax, struct cplan9_s *hmm,
				   int use_easel);

#ifdef USE_MPI
extern void parallel_make_histogram (int *gc_count, int *partitions, int num_partitions,
				     CM_t *cm, int D, int num_samples,
				     int sample_length, double *lambda, double *K, 
				     int *dmin, int *dmax, 
				     int mpi_my_rank, int num_procs, 
				     int mpi_master_rank);
#endif
extern void GetDBInfo(ESL_SQFILE *sqfp, long *ret_N, int **ret_gc_ct);

extern float e_to_score (float E, double *mu, double *lambda);

extern double RJK_ExtremeValueE (float x, double mu, double lambda);

extern char resolve_degenerate (char c);


#endif

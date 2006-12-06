/*
 * stats.h
 * 
 * Header file for stats.c
 */

#ifndef _stats_h
#define _stats_h

extern void serial_make_histogram (int *gc_count, int *partitions, int num_partitions,
				   CM_t *cm, int D, int num_samples, 
				   int sample_length, double *lambda, double *K,
				   int *dmin, int *dmax, struct cplan9_s *hmm,
				   int use_easel);

#ifdef USE_MPI
extern void parallel_make_histogram (int *gc_count, int *partitions, int num_partitions,
				     CM_t *cm, int D, int num_samples,
				     int sample_length, double *lambda, double *K, 
				     int mpi_my_rank, int num_procs, 
				     int mpi_master_rank);
#endif

extern void get_dbinfo (SQFILE *dbfp, long *N, int *gc_count);

extern float e_to_score (float E, float *mu, float *lambda);

extern double RJK_ExtremeValueE (float x, double mu, double lambda);

extern char resolve_degenerate (char c);


#endif

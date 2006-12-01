/*
 * stats.h
 * 
 * Header file for stats.c
 */

#ifndef _stats_h
#define _stats_h

void serial_make_histogram (int *gc_count, int *partitions, int num_partitions,
                            CM_t *cm, int D, int num_samples, 
			    int sample_length, double *lambda, double *K,
			    int *dmin, int *dmax, struct cplan9_s *hmm);

#ifdef USE_MPI
void parallel_make_histogram (int *gc_count, int *partitions, int num_partitions,
			      CM_t *cm, int D, int num_samples,
			      int sample_length, float *lambda, float *K, 
			      int mpi_my_rank, int num_procs, 
			      int mpi_master_rank);
#endif

void get_dbinfo (SQFILE *dbfp, long *N, int *gc_count);

float e_to_score (float E, float *mu, float *lambda);

double RJK_ExtremeValueE (float x, double mu, double lambda);

char resolve_degenerate (char c);


#endif

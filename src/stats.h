/*
 * stats.h
 * 
 * Header file for stats.c
 */

#ifndef _stats_h
#define _stats_h
#include "esl_sqio.h"
#include "config.h"

extern CMStats_t *AllocCMStats(int np);
extern void FreeCMStats(CMStats_t *cmstats);
extern int debug_print_cmstats(CMStats_t *cmstats);
extern int debug_print_evdinfo(EVDInfo_t *evd);
extern int debug_print_filterthrinfo(CMStats_t *cmstats, CP9FilterThr_t *fthr);

extern int  get_gc_comp(char *seq, int start, int stop);

extern void serial_make_histogram (int *gc_count, int *partitions, int num_partitions,
				   CM_t *cm, int num_samples, 
				   int sample_length, int doing_cp9_stats,
				   int use_easel);

#ifdef USE_MPI
void parallel_make_histogram (int *gc_count, int *partitions, int num_partitions, 
			      CM_t *cm, int num_samples, int sample_length,
			      int doing_cp9_stats,
			      int mpi_my_rank, int mpi_num_procs, 
			      int mpi_master_rank);
#endif

extern void GetDBInfo(ESL_SQFILE *sqfp, long *ret_N, int **ret_gc_ct);

extern float e_to_score (float E, double *mu, double *lambda);

extern double RJK_ExtremeValueE (float x, double mu, double lambda);

extern char resolve_degenerate (char c);

extern float MinCMScCutoff (CM_t *cm);
extern float MinCP9ScCutoff (CM_t *cm);
extern int   CM2EVD_mode(CM_t *cm, int *ret_cm_evd_mode, int *ret_cp9_evd_mode);

#endif

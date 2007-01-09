/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_wrappers.h
 * 
 * Functions that actually do the work for cmalign and cmsearch,
 * in parallel and serial varieties.
 *
 * Eric Nawrocki
 */

#ifndef CMWRAPPERS_INCLUDED
#define CMWRAPPERS_INCLUDED

#include "config.h"
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "squid.h"

extern void 
serial_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons, int cutoff_type, 
			float cutoff, double *mu, double *lambda);
extern void 
parallel_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons,
			  int cutoff_type, float cutoff, 
			  double *mu, double *lambda, 
			  int mpi_my_rank, int mpi_master_rank, 
			  int mpi_num_procs) ;
extern float
actually_search_target(CM_t *cm, char *dsq, int i0, int j0, float cutoff, 
		       scan_results_t *results, int do_filter, int *ret_flen);

extern void
serial_align_targets(ESL_SQFILE *seqfp, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, 
		     char ***ret_postcode, int *ret_nseq, int bdump_level, int debug_level, 
		     int silent_mode);
extern void
parallel_align_targets(ESL_SQFILE *seqfp, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr,
		       char ***ret_postcode, int *ret_nseq, int bdump_level, int debug_level,
		       int silent_mode, int mpi_my_rank, int mpi_master_rank, int mpi_num_procs);
extern void
actually_align_targets(CM_t *cm, ESL_SQ **sq, int nseq, Parsetree_t ***ret_tr, char ***ret_postcode,
		       int bdump_level, int debug_level, int silent_mode);


#endif

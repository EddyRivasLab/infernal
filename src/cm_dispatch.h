/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_dispatch.h
 * 
 * Functions that actually do the work for cmalign and cmsearch,
 * in parallel and serial varieties.
 *
 * Eric Nawrocki
 */

#ifndef CMDISPATCH_INCLUDED
#define CMDISPATCH_INCLUDED

#include "esl_config.h"
#include "config.h"

#include "easel.h"
#include "esl_sqio.h"

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */

extern void 
serial_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons);
extern void 
parallel_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons,
			  int mpi_my_rank, int mpi_master_rank, int mpi_num_procs) ;
extern float 
actually_search_target(CM_t *cm, char *dsq, int i0, int j0, float cm_cutoff, 
		       float cp9_cutoff, scan_results_t *results, int do_filter, 
		       int doing_cm_stats, int doing_cp9_stats, int *ret_flen);
extern void
serial_align_targets(ESL_SQFILE *seqfp, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, 
		     char ***ret_postcode, CP9trace_t ***ret_cp9_tr, int *ret_nseq,
		     int bdump_level, int debug_level, int silent_mode);
extern void
parallel_align_targets(ESL_SQFILE *seqfp, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr,
		       char ***ret_postcode, CP9trace_t ***ret_cp9_tr, int *ret_nseq,
		       int bdump_level, int debug_level,
		       int silent_mode, int mpi_my_rank, int mpi_master_rank, int mpi_num_procs);
extern void
actually_align_targets(CM_t *cm, ESL_SQ **sq, int nseq, Parsetree_t ***ret_tr, char ***ret_postcode,
		       CP9trace_t ***ret_cp9_tr, int bdump_level, int debug_level, int silent_mode);

extern int PrintSearchInfo(FILE *fp, CM_t *cm, int cm_mode, int cp9_mode, long N);

#endif

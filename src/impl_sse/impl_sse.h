/* SSE optimized implementation for Infernal
 */

#ifndef CM_IMPL_SSE_INCLUDED
#define CM_IMPL_SSE_INCLUDED

#include <xmmintrin.h>
#include <emmintrin.h>

#include "config.h"
#include "structs.h"
#include "funcs.h"

/*****************************************************************
 * 1. CM_OPTIMIZED: an optimized score profile
 *****************************************************************/

typedef struct cm_optimized_s {
  /* basic model information */
  int       M;
  char     *sttype;
  int      *cfirst;
  int      *cnum;

  const ESL_ALPHABET *abc;

  /* epi16-specific information */
  float     scale_w;
  int16_t **tsc;
  int16_t **oesc;
  int16_t  *beginsc;
  int16_t  *endsc;
  int16_t   el_selfsc;
} CM_OPTIMIZED;

/*****************************************************************
 * 2. Declarations of the external API.
 *****************************************************************/

/* cm_optimized.c */
int cm_optimized_Convert(const CM_t *cm, CM_OPTIMIZED *ocm);

/* sse_cm_dpsearch.c */
int SSECYKScan(CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq,
	int i0, int j0, float cutoff, search_results_t *results,
	int do_null3, float **ret_vsc, float *ret_sc);

/* sse_cm_dpsmall.c */
float SSE_CYKInside(CM_t *cm, ESL_DSQ *dsq, int L, int r,
	int i0, int j0, Parsetree_t **ret_tr);
float SSE_CYKInsideScore(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0);
float SSE_CYKDemands(CM_t *cm, int L, int be_quiet);
float SSE_CYKDivideAndConquer(CM_t *cm, ESL_DSQ *dsq, int L, int r,
	int i0, int j0, Parsetree_t **ret_tr);

/* sse_util.c */
inline __m128  alt_rightshift_ps(__m128 a, __m128 b);
inline __m128i sse_leftshift_epi16(__m128i a, __m128i b);
inline __m128i sse_rightshift_epi16(__m128i a, __m128i b);

#endif /* CM_IMPL_SSE_INCLUDED */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/


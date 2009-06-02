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
 * 2. CM_CONSENSUS: a match-only profile
 *****************************************************************/

typedef struct cm_consensus_s {
  int       M;
  char     *sttype;
  int      *next; /* Overloaded: we'll usually fill this in, even
                     though next[i] = i+1 for most cases.  For
                     B_st, next will be the _right_ child (and
                     the left will always be i+1                */

  ESL_ALPHABET *abc;

  /* floating point scores, if necessary */

  /* Reduced-precision uchar scores */
  float     scale_b;
  uint8_t   base_b;
//uint8_t   bias_b;
  uint8_t **oesc;
} CM_CONSENSUS;

/*****************************************************************
 * 3. Declarations of the external API.
 *****************************************************************/

/* cm_optimized.c */
CM_OPTIMIZED* cm_optimized_Convert(const CM_t *cm);
void cm_optimized_Free(CM_OPTIMIZED *ocm);

CM_CONSENSUS* cm_consensus_Convert(CM_t *cm);
void cm_consensus_Free(CM_CONSENSUS *ccm);

/* sse_cm_dpsearch.c */
int SSE_CYKScan(CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq,
	int i0, int j0, float cutoff, search_results_t *results,
	int do_null3, float **ret_vsc, float *ret_sc);

/* sse_cm_dpsmall.c */
float SSE_CYKInside(CM_t *cm, ESL_DSQ *dsq, int L, int r,
	int i0, int j0, Parsetree_t **ret_tr);
float SSE_CYKInsideScore(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0);
float SSE_CYKDemands(CM_t *cm, int L, int be_quiet);
float SSE_CYKDivideAndConquer(CM_t *cm, ESL_DSQ *dsq, int L, int r,
	int i0, int j0, Parsetree_t **ret_tr);
int SSE_CYKFilter_epi16(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
        int allow_begin, int *ret_b, int *ret_bsc);


/* sse_util.c */
inline __m128  alt_rightshift_ps(__m128 a, __m128 b);
inline __m128i sse_setlw_neginfv(__m128i a);

#endif /* CM_IMPL_SSE_INCLUDED */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/


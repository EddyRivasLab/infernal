/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_postprob.h
 * 
 * Functions for working with posterior probabilities for CMs.
 *
 */

#ifndef CMPOSTPROB_INCLUDED
#define CMPOSTPROB_INCLUDED

#include "config.h"
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "squid.h"

extern float FInside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
		     float ***alpha, float ****ret_alpha, 
		     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		     int allow_begin);
extern float FOutside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
		      float ***beta, float ****ret_beta, 
		      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		      int allow_begin, float ***alpha, int do_check);
#endif

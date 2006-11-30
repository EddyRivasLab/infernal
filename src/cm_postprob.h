/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_postprob.h
 * 
 * Functions for working with posterior probabilities for CMs.
 * Eric Nawrocki
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
		      int allow_begin, float ***alpha, float ****ret_alpha, int do_check);
extern void  CMPosterior(int L, CM_t *cm, float ***alpha, float ****ret_alpha, float ***beta, 
			 float ****ret_beta, float ***post, float ****ret_post);
extern char *CMPostalCode(CM_t *cm, int L, float ***post, Parsetree_t *tr);
extern char Fscore2postcode(float sc);
extern float FScore2Prob(float sc, float null);
extern void CMCheckPosterior(int L, CM_t *cm, float ***post);
extern float FInside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
			     float ***alpha, float ****ret_alpha, 
			     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			     int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern float FOutside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
			      float ***beta, float ****ret_beta, 
			      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			      int allow_begin, float ***alpha, float ****ret_alpha, 
			      int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void CMPosterior_b_jd_me(int L, CM_t *cm, float ***alpha, float ****ret_alpha, 
				float ***beta, float ****ret_beta, float ***post, float ****ret_post,
				int *jmin, int *jmax, int **hdmin, int **hdmax);
extern char *CMPostalCode_b_jd_me(CM_t *cm, int L, float ***post, Parsetree_t *tr,
				int *jmin, int *jmax, int **hdmin, int **hdmax);
     
#endif

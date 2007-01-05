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
extern float IInside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
		     int ***alpha, int ****ret_alpha, 
		     struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
		     int allow_begin);
extern float FOutside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
		      float ***beta, float ****ret_beta, 
		      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		      int allow_begin, float ***alpha, float ****ret_alpha, int do_check);
extern float IOutside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
		      int ***beta, int ****ret_beta, 
		      struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
		      int allow_begin, int ***alpha, int ****ret_alpha, int do_check);
extern void   CMPosterior(int L, CM_t *cm, float ***alpha, float ****ret_alpha, float ***beta, 
			  float ****ret_beta, float ***post, float ****ret_post);
extern void  ICMPosterior(int L, CM_t *cm, int ***alpha, int ****ret_alpha, int ***beta, 
			 int ****ret_beta, int ***post, int ****ret_post);
extern char  *CMPostalCode(CM_t *cm, int L, float ***post, Parsetree_t *tr);
extern char *ICMPostalCode(CM_t *cm, int L, int ***post, Parsetree_t *tr);
extern char Fscore2postcode(float sc);
extern char Iscore2postcode(int sc);
extern float FScore2Prob(float sc, float null);
extern void  CMCheckPosterior(int L, CM_t *cm, float ***post);
extern void ICMCheckPosterior(int L, CM_t *cm, int ***post);

extern float FInside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
			     float ***alpha, float ****ret_alpha, 
			     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			     int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern float IInside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
			     int ***alpha, int ****ret_alpha, 
			     struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
			     int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern float FOutside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
			      float ***beta, float ****ret_beta, 
			      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			      int allow_begin, float ***alpha, float ****ret_alpha, 
			      int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern float IOutside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
			      int ***beta, int ****ret_beta, 
			      struct Ideckpool_s *dpool, struct Ideckpool_s **ret_dpool,
			      int allow_begin, int ***alpha, int ****ret_alpha, 
			      int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void  CMPosterior_b_jd_me(int L, CM_t *cm, float ***alpha, float ****ret_alpha, 
				 float ***beta, float ****ret_beta, float ***post, float ****ret_post,
				 int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void ICMPosterior_b_jd_me(int L, CM_t *cm, int ***alpha, int ****ret_alpha, 
				 int ***beta, int ****ret_beta, int ***post, int ****ret_post,
				 int *jmin, int *jmax, int **hdmin, int **hdmax);
extern char  *CMPostalCode_b_jd_me(CM_t *cm, int L, float ***post, Parsetree_t *tr,
				   int *jmin, int *jmax, int **hdmin, int **hdmax);
extern char *ICMPostalCode_b_jd_me(CM_t *cm, int L, int ***post, Parsetree_t *tr,
				   int *jmin, int *jmax, int **hdmin, int **hdmax);
     
/* And new memory management routines analogous to those in smallcyk.c for
 * handling scaled int log odds scores instead of floats. 
 */
extern Ideckpool_t *Ideckpool_create(void);
extern void    Ideckpool_push(struct Ideckpool_s *dpool, int **deck);
extern int     Ideckpool_pop(struct Ideckpool_s *d, int ***ret_deck);
extern void    Ideckpool_free(struct Ideckpool_s *d);
extern int   **Ialloc_vjd_deck(int L, int i, int j);
extern int     Isize_vjd_deck(int L, int i, int j);
extern void    Ifree_vjd_deck(int **a, int i, int j);
extern void    Ifree_vjd_matrix(int ***a, int M, int i, int j);

extern int ** Ialloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
#endif

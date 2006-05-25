/************************************************************
 * @LICENSE@
 ************************************************************/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "hmmer_funcs.h"
#include "hmmer_structs.h"

/* cm_postprob.c
 * 
 * Functions for working with posterior probabilities for CMs.
 * EPN 05.08.06
 */

/*****************************************************************
 * CM Inside() & Outside() functions.
 * Important: smallcyk.c has functions named Inside(), Outside() or some 
 *            variation, but these are implementations of variations of
 *            the CYK-Inside functions, as mentioned in the 2002 D&C
 *            BMC Bioinformatics publication.
 *****************************************************************/ 

/* Function: FInside()
 * 
 * Purpose:  The Inside dynamic programming algorithm for CMs.
 *           Works directly with floats, stepping into and out of
 *           log space as necessary.
 * 
 *           Based on inside() in smallcyk.c, with the following 
 *           differences: necessarily align the sequence to the 
 *           full model (not possible to align to subtrees as in
 *           smallcyk.c's inside()), also no shadow matrix is
 *           kept because we're interested in ALL paths, finally
 *           we don't care about the best local begin state for the
 *           same reason.
 *  
 * Purpose:  Run the Inside() alignment algorithm, on a 
 *           subsequence from i0..j0, using the entire model.
 * 
 *           A note on the loop conventions. We're going to keep the
 *           sequence (dsq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: i0-1..j in the rows, 
 *           and up to j0-i0+1 in the columns (d dimension). Where this is 
 *           handled the most is in two variables: W, which is the length of 
 *           the subsequence (j0-i0+1), and is oft used in place of L in the 
 *           usual CYK; and jp (read: j'), which is the *relative* j w.r.t. the
 *           subsequence, ranging from 0..W, and then d ranges from 
 *           0 to jp, and j is calculated from jp (i0-1+jp).
 *           
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in a dpool, 
 *           the decks *must* be sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           alpha     - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_alpha - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *                       
 *
 * Returns: log P(S|M)/P(S|R), as a bit score.
 */
float 
FInside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
       float ***alpha, float ****ret_alpha, 
       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
       int allow_begin)
{
  float  **end;         /* we re-use the end deck. */
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  float    bsc;		/* total score for using local begin states */
  
  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */

  if (dpool == NULL) dpool = deckpool_create();
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, 0, E_st);
  for (jp = 0; jp <= W; jp++) {
    j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
    end[j][0] = 0.;
    for (d = 1; d <= jp; d++) end[j][d] = IMPOSSIBLE;
  }

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    alpha = MallocOrDie(sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v <= (cm->M-1); v++) touch[v] = cm->pnum[v];

  /* Main recursion
   */
  for (v = (cm->M - 1); v >= 0; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      } 
      if (! deckpool_pop(dpool, &(alpha[v]))) 
	alpha[v] = alloc_vjd_deck(L, i0, j0);

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y+yoffset][j][d] 
							   + cm->tsc[v][yoffset]));
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == B_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		z = cm->cnum[v];
		  
		alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0];
		for (k = 1; k <= d; k++)
		  alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y][j-k][d-k] 
							   + alpha[z][j][k]));
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    if (jp > 0) alpha[v][j][1] = IMPOSSIBLE;
	    for (d = 2; d <= jp; d++) 
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y+yoffset][j-1][d-2] 
							   + cm->tsc[v][yoffset]));
		i = j-d+1;
		if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    for (d = 1; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y+yoffset][j][d-1] 
							   + cm->tsc[v][yoffset]));
		i = j-d+1;
		if (dsq[i] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    for (d = 1; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y+yoffset][j-1][d-1] 
							   + cm->tsc[v][yoffset]));
		if (dsq[j] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}				/* finished calculating deck v. */
      
      /* Keep track of contributions of local begins */
      if (allow_begin)
	{
	  bsc = LogSum2(bsc, (alpha[v][j0][W] + cm->beginsc[v]));
	}

      /* If we're at the root state, record contribution of local begins */
      if (allow_begin && v == 0)
	{
	  alpha[0][j0][W] = LogSum2(alpha[0][j0][W], bsc);
	}	  

      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { /* we can definitely release the S children of a bifurc. */
	    y = cm->cfirst[v]; deckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	    z = cm->cnum[v];   deckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) 
		  {
		    if (cm->sttype[y] == E_st) { 
		      nends--; 
		      if (nends == 0) { deckpool_push(dpool, end); end = NULL;}
		    } else 
		      deckpool_push(dpool, alpha[y]);
		    alpha[y] = NULL;
		  }
	      }
	  }
      }
  } /* end loop over all v */

  /* debug_print_alpha(alpha, cm, L);*/

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  sc = alpha[0][j0][W];

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { deckpool_push(dpool, alpha[v]); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { deckpool_push(dpool, end); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (deckpool_pop(dpool, &end)) free_vjd_deck(end, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  printf("***returning from FInside() sc : %f\n", sc);
  return sc;
}

/***********************************************************************
 * Function: FOutside()
 * Date:     EPN 05.08.06
 * SRE, Tue Aug  8 10:42:52 2000 [St. Louis]
 *
 * Purpose:  The Outside dynamic programming algorithm for CMs.
 *           Works directly with floats, stepping into and out of 
 *           log space as necessary.
 *  
 *           Derived from smallcyk.c::CYKOutside() and smallcyk.c::outsdie(). 
 * 
 *           Align a subsequence to the full model, i.e. we're given
 *           i0 and j0, beginning and end positions of the subseq we're
 *           considering.
 * 
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in alpha, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           alpha     - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_alpha - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           alpha      - the alpha matrix from a Inside run, if do_check is FALSE
 *                        only decks for S states must be valid, else all must be
 *                        valid.
 *           do_check   - TRUE to check that probabilities are correctly calc'ed
 *                        that is they give correct P(S|M) (at least same as in
 *                        alpha[0][L][L] from the Inside run.
 * Returns: log P(S|M)/P(S|R), as a bit score.
 ***********************************************************************/
float 
FOutside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
	 float ***beta, float ****ret_beta, 
	 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	 int allow_begin, float ***alpha, int do_check)
{
  int      v,y,z;	/* indices for states */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  float    escore;	/* an emission score, tmp variable */
  int      voffset;	/* index of v in t_v(y) transition scores */
  int      W;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  float    bsc;		/* total score for using local begin states */
  float    return_sc;   /* P(S|M) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */

  /* Allocations and initializations
   */
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* if caller didn't give us a matrix, make one.
   * Allocate room for M+1 decks because we might need the EL deck (M)
   * if we're doing local alignment.
   */
  if (beta == NULL) {
    beta = MallocOrDie(sizeof(float **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck. Root is necessarily the ROOT_S state 0.
   */
  if (! deckpool_pop(dpool, &(beta[0])))
    beta[0] = alloc_vjd_deck(L, i0, j0);
  for (jp = 0; jp <= W; jp++) {
    j = i0-1+jp;
    for (d = 0; d <= jp; d++)
      beta[0][j][d] = IMPOSSIBLE;
  }
  beta[0][j0][W] = 0;		

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   */
  if (cm->flags & CM_LOCAL_END) {
    if (! deckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = alloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[cm->M][j][d] = IMPOSSIBLE;
    }
    
    /* We don't have to worry about vroot -> EL transitions the way 
     * smallcyk.c::outside() does, because vroot = 0.
     */
  }

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  /* Main loop down through the decks
   */
  /*for (v = 2; v < cm->M; v++) */ /*EPN is this 2 b/c Durbin p.287 
				     has state 2 in the algorithm? b/c state 1 is root*/
  for (v = 1; v < cm->M; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_vjd_deck(L, i0, j0);

      /* Init the whole deck to IMPOSSIBLE
       */
      for (jp = W; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  beta[v][j][d] = IMPOSSIBLE;
      }

      /* If we can do a local begin into v, also init with that. 
       * By definition, beta[0][j0][W] == 0.
       */ 
      if (i0 == 1 && j0 == L && (cm->flags & CM_LOCAL_BEGIN))
	beta[v][j0][W] = cm->beginsc[v];

      /* main recursion:
       */
      for (jp = W; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  {
	    if (cm->stid[v] == BEGL_S) 
	      {
		y = cm->plast[v];	/* the parent bifurcation    */
		z = cm->cnum[y];	/* the other (right) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j][0]; /* init on k=0 */
		for (k = 1; k <= L-j; k++)
		  beta[v][j][d] = LogSum2(beta[v][j][d], (beta[y][j+k][d+k] + alpha[z][j+k][k]));
	      }
	    else if (cm->stid[v] == BEGR_S) 
	      {
		y = cm->plast[v];	        /* the parent bifurcation    */
		z = cm->cfirst[y];	/* the other (left) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j-d][0];	/* init on k=0 */
		for (k = 1; k <= j-d; k++) 
		  beta[v][j][d] = LogSum2(beta[v][j][d], (beta[y][j][d+k] + alpha[z][j-d][k]));
	      }
	    else
	      {
		if(!(do_check)) alpha[v][j][d] = IMPOSSIBLE;
		i = j-d+1;
		for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
		  voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */
		  switch(cm->sttype[y]) {
		  case MP_st: 
		    if (j == j0 || d == jp) continue; /* boundary condition */
		    
		    if (dsq[i-1] < Alphabet_size && dsq[j+1] > Alphabet_size)
		      escore = cm->esc[y][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
		    else
		      escore = DegeneratePairScore(cm->esc[y], dsq[i-1], dsq[j+1]);
		    beta[v][j][d] = LogSum2(beta[v][j][d], (beta[y][j+1][d+2] + cm->tsc[y][voffset]
							   + escore));
		    break;
		    
		  case ML_st:
		  case IL_st: 
		    if (d == jp) continue;	/* boundary condition (note when j=0, d=0)*/
		    
		    if (dsq[i-1] < Alphabet_size) 
		      escore = cm->esc[y][(int) dsq[i-1]];
		    else
		      escore = DegenerateSingletScore(cm->esc[y], dsq[i-1]);
		    beta[v][j][d] = LogSum2(beta[v][j][d], (beta[y][j][d+1] + cm->tsc[y][voffset] 
							   + escore));
		    break;
		    
		  case MR_st:
		  case IR_st:
		    if (j == j0) continue;
		    
		    if (dsq[j+1] < Alphabet_size) 
		      escore = cm->esc[y][(int) dsq[j+1]];
		    else
		      escore = DegenerateSingletScore(cm->esc[y], dsq[j+1]);
		    beta[v][j][d] = LogSum2(beta[v][j][d], (beta[y][j+1][d+1] + cm->tsc[y][voffset] 
							   + escore));
		    break;
		    
		  case S_st:
		  case E_st:
		  case D_st:
		    beta[v][j][d] = LogSum2(beta[v][j][d], (beta[y][j][d] + cm->tsc[y][voffset])); 
		      break;
		    
		  default: Die("bogus child state %d\n", cm->sttype[y]);
		  }/* end switch over states*/
		} /* ends for loop over parent states. we now know beta[v][j][d] for this d */
		if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;
	      }	/* ends else entered for non-BEGL/BEGR states*/	
	  } /* ends loop over d. We know all beta[v][j][d] in this row j*/
      }/* end loop over jp. We know the beta's for the whole deck.*/
	
      /* Deal with local alignment end transitions v->EL
       * (EL = deck at M.)
       */
      if (NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (jp = 0; jp <= W; jp++) { 
	  j = i0-1+jp;
	  for (d = 0; d <= jp; d++) 
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st: 
		if (j == j0 || d == jp) continue; /* boundary condition */
		if (dsq[i-1] < Alphabet_size && dsq[j+1] > Alphabet_size)
		  escore = cm->esc[v][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->esc[v], dsq[i-1], dsq[j+1]);
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][j+1][d+2] + cm->endsc[v] 
							       + (cm->el_selfsc * d) + escore));
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;	
		if (dsq[i-1] < Alphabet_size) 
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = DegenerateSingletScore(cm->esc[v], dsq[i-1]);
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][j][d+1] + cm->endsc[v] 
							       + (cm->el_selfsc * d) + escore));
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < Alphabet_size) 
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = DegenerateSingletScore(cm->esc[v], dsq[j+1]);
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][j+1][d+1] + cm->endsc[v] 
							       + (cm->el_selfsc * d) + escore));
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][j][d] + cm->endsc[v] 
							       + (cm->el_selfsc * d) + escore));
		break;
	      case B_st:  
	      default: Die("bogus parent state %d\n", cm->sttype[v]);
		/* note that although B is a valid vend for a segment we'd do
                   outside on, B->EL is set to be impossible, by the local alignment
                   config. There's no point in having a B->EL because B is a nonemitter
                   (indeed, it would introduce an alignment ambiguity). The same
		   alignment case is handled by the X->EL transition where X is the
		   parent consensus state (S, MP, ML, or MR) above the B. Thus,
		   this code is relying on the NOT_IMPOSSIBLE() test, above,
		   to make sure the sttype[vend]=B case gets into this switch.
		*/
	      } /* end switch over parent state type v */
	    } /* end inner loop over d */
	} /* end outer loop over jp */
      } /* end conditional section for dealing w/ v->EL local end transitions */

      /* Look at v's parents; if we're reusing memory (! do_full)
       * push the parents that we don't need any more into the pool.
       */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { deckpool_push(dpool, beta[y]); beta[y] = NULL; }
	}
      }
    } /* end loop over decks v. */

  /* EPN left the below SRE block out */
#if 0
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CM_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	if ((sc = beta[cm->M][j][d+1]) > beta[cm->M][j][d])
	  beta[cm->M][j][d] = sc;
    }
  }
#endif

  /*debug_print_alpha(beta, cm, L);*/

  if(do_check)
    {
      /* Determine P(S|M) (probability of the sequence given the model) 
       * using both the Outside (beta) and Inside (alpha) matrices,
       * and ensure they're consistent with P(S|M) from the Inside calculation.
       * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
       *                                                = P(S|M)  
       * in v,j,d coordinates this is:
       * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
       *                                                = P(S|M) 
       */
	 
      for(n = 0; n < cm->nodes; n++)
	{
	  sc = IMPOSSIBLE;
	  num_split_states = SplitStatesInNode(cm->ndtype[n]);
	  for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++)
	    {
	      for (jp = 0; jp <= W; jp++) 
		{ 
		  j = i0-1+jp;
		  for (d = 0; d <= jp; d++) 
		    {
		      sc = LogSum2(sc, (alpha[v][j][d] + beta[v][j][d]));
		    }
		}
	    }
	  diff = alpha[0][j0][W] - sc;
	  if(diff > 0.001 || diff < -0.001)
	    {
	      Die("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
		  n, sc, alpha[0][j0][W], diff);
	    }
	}

    }

  printf("        inside score: %10.2f\n", alpha[0][j0][W]);
  /* Check 2: Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) */

  /* Calculation P(S|M) given only the beta matrix, we use the deck for the final END_E state,
   * but the deck for any END_E state could be used. 
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M)
   * NOTE: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M)
   */
  return_sc = IMPOSSIBLE;
  for (jp = 0; jp <= W; jp++) 
    { 
      j = i0-1+jp;
      /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
      return_sc = LogSum2(return_sc, (beta[cm->M-1][j][0]));
    }
  printf("        check2 score: %10.2f\n", return_sc);

  /* If the caller doesn't want the matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CM_LOCAL_END) {
      deckpool_push(dpool, beta[cm->M]);
      beta[cm->M] = NULL; 
    }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    float **a;
    while (deckpool_pop(dpool, &a)) free_vjd_deck(a, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }
  free(touch);

  printf("***returning from FOutside() sc : %f\n", return_sc);
  return return_sc;
}
/***************************************************************/

/* Function: LogSum2()
 * 
 * Purpose:  Returns the log_2 of the sum of two log_2 probabilities.
 *           log(exp(p1)+exp(p2)) = p1 + log(1 + exp(p2-p1)) for p1 > p2
 *           Note that this is in log_2 space.
 */
float 
LogSum2(float p1, float p2)
{
  if (p1 > p2)
    return (p1-p2 > 50.) ? p1 : p1 + sreLOG2(1. + pow(2.,(p2-p1)));
  else
    return (p2-p1 > 50.) ? p2 : p2 + sreLOG2(1. + pow(2.,(p1-p2)));
}



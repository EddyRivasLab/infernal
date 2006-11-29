/***********************************************************
 * @LICENSE@
 ************************************************************/
/* cm_postprob.c
 * EPN 05.08.06
 * 
 * Functions for working with posterior probabilities for CMs.
 * Includes non-banded functions as well as banded ones (bands 
 * in the j and d dimensions)

 */

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


/*****************************************************************
 * CM FInside() & FOutside() functions.
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

  /* EPN: now that the EL self loop has a transition score, its
   *      necessary to keep track of the alpha EL deck (alpha[cm->M])
   */
  if(cm->flags & CM_LOCAL_BEGIN)
    {
      if (! deckpool_pop(dpool, &(alpha[cm->M]))) 
	alpha[cm->M] = alloc_vjd_deck(L, i0, j0);
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	alpha[cm->M][j][0] = 0.;
	/*alpha[cm->M][j][0] = IMPOSSIBLE;*/
	for (d = 1; d <= jp; d++)
	  {
	    alpha[cm->M][j][d] = (cm->el_selfsc * (d));
	  }
      }
    }

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
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y][j][d] 
							      + cm->tsc[v][yoffset]));
		  }
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
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y][j-1][d-2] 
							      + cm->tsc[v][yoffset]));
		  }
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
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y][j][d-1] 
							    + cm->tsc[v][yoffset]));
		  }
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
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    alpha[v][j][d] = LogSum2(alpha[v][j][d], (alpha[y][j-1][d-1] 
							      + cm->tsc[v][yoffset]));
		  }
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

  /*
    printf("INSIDE, printing alpha\n");
    debug_print_alpha(alpha, cm, L);
    printf("INSIDE, done printing alpha\n");
  */

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
  printf("***returning from FInside() sc  : %f\n", sc);
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
 *           do_full   - if TRUE, we save all the decks in alpha and beta, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           beta       - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_beta  - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           alpha     - the alpha matrix from a Inside run, if do_check is FALSE
 *                        only decks for S states must be valid, else all must be
 *                        valid.
 *           ret_alpha - if non-NULL, return the alpha matrix with one or more
 *                       decks available for examination (see "do_full")
 *           do_check  - TRUE to do time-consuming check to make sure
 *                       beta and alpha are consistent (only if NON-LOCAL mode)
 * 
 * Returns: log P(S|M)/P(S|R), as a bit score.
 ***********************************************************************/
float 
FOutside(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
	 float ***beta, float ****ret_beta, 
	 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	 int allow_begin, float ***alpha, float ****ret_alpha, 
	 int do_check)
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
  float    return_sc;   /* P(S|M)/P(S|R) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  float  **end;         /* we re-use the end deck. */

  if (cm->flags & CM_LOCAL_END) { do_check = FALSE; } 
  /* Code for checking doesn't apply in local mode. See below. */

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
								+ escore));
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;	
		if (dsq[i-1] < Alphabet_size) 
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = DegenerateSingletScore(cm->esc[v], dsq[i-1]);
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][j][d+1] + cm->endsc[v] 
								+ escore));
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < Alphabet_size) 
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = DegenerateSingletScore(cm->esc[v], dsq[j+1]);
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][j+1][d+1] + cm->endsc[v]
								+ escore));
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][j][d] + cm->endsc[v]));
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

  /* EPN this code block is not superfluous for Inside() */
  /*     below are SRE's notes from a CYK inside() function */
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CM_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->el_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE, printing beta\n");
    debug_print_alpha(beta, cm, L);
    printf("OUTSIDE, done printing beta\n");
  */

  if(do_check && (!(cm->flags & CM_LOCAL_END))) 
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    {
      /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
       * using both the Outside (beta) and Inside (alpha) matrices,
       * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
       * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
       *                                                = P(S|M) / P(S|R)  
       * in v,j,d coordinates this is:
       * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
       *                                                = P(S|M) / P(S|R)
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

  /* IF not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state with 
   * d = 0.
   * We pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * IF local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */

  if(!(cm->flags & CM_LOCAL_END))
    {
      return_sc = IMPOSSIBLE;
      for (jp = 0; jp <= W; jp++) 
	{ 
	  j = i0-1+jp;
	  /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	  return_sc = LogSum2(return_sc, (beta[cm->M-1][j][0]));
	}
    }
  else /* return_sc = P(S|M) / P(S|R) from Inside() */
    {
      return_sc = alpha[0][j0][W];
    }
  /* If the caller doesn't want the beta matrix, free it.
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

  /* If the caller doesn't want the alpha matrix, free it 
   * Else, pass it back to him.
   * EPN - if we free the alpha and beta matrix the deck pool has all the 
   *       decks from alpha and beta, not sure if this is desirable.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
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
    float **a;
    while (deckpool_pop(dpool, &a)) free_vjd_deck(a, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }
  free(touch);

  if(!(cm->flags & CM_LOCAL_END))
    printf("***returning from FOutside() sc : %f\n", return_sc);
  else
    printf("***returning from FOutside() sc : %f (LOCAL mode; sc is from Inside)\n", return_sc);
  return return_sc;
}

/* Function: FScore2Prob()
 * 
 * Purpose:  Convert a float log_2 odds score back to a probability;
 *           needs the null model probability, if any, to do the conversion.
 */
float 
FScore2Prob(float sc, float null)
{
  /*printf("in FScore2Prob: %10.2f sreEXP2: %10.2f\n", sc, (sreEXP2(sc)));*/
  if (sc == IMPOSSIBLE) return 0.;
  else              return (null * sreEXP2(sc));
}

/***************************************************************/
/* Function: CMPosterior()
 *           
 * EPN 05.25.06 based on IHH's P7EmitterPosterior() from HMMER's postprob.c
 *
 * Purpose:  Combines Inside and Outside cubes into a posterior
 *           probability cube.
 *           The entry in post[v][j][d] is the log of the
 *           posterior probability of a parse subtree rooted at v 
 *           emitting the subsequence i..j (i=j-d+1).
 *           The caller must allocate space for the cube, although the
 *           beta matrix from Outside can be used instead (overwriting it will not
 *           compromise the algorithm).
 *           
 * Args:     L        - length of sequence
 *           cm       - the model
 *           alpha    - pre-calculated Inside matrix 
 *           ret_alpha - if non-NULL, return the matrix as it was passed in,
 *                       else free it.
 *           beta     - pre-calculated Outside matrix
 *           ret_beta - if non-NULL, return the matrix as it was passed in,
 *                       else free it.
 *           post     - pre-allocated dynamic programming cube
 *           ret_post - if non-NULL, return the posterior matrix,
 *                       else free it.
 *           
 * Return:   void
 */
void
CMPosterior(int L, CM_t *cm, float ***alpha, float ****ret_alpha, float ***beta, float ****ret_beta, 
	    float ***post, float ****ret_post)
{
  int   v, j, d;
  float sc;
  float  **end;         /* used for freeing alpha b/c we re-use the end deck. */
  int vmax;
  sc = alpha[0][L][L];
  
  /* If local ends are on, start with the EL state (cm->M), otherwise
   * its not a valid deck.
   */
  vmax = cm->M-1;
  if (cm->flags & CM_LOCAL_END) vmax = cm->M;

  for (v = vmax; v >= 0; v--) 
    for (j = 0; j <= L; j++) 
      for (d = 0; d <= j; d++)
	{
	  post[v][j][d] = alpha[v][j][d] + beta[v][j][d] - sc;
	  if(v == vmax)
	    {
	      /*printf("v: %3d | j: %3d | d: %3d | alpha: %5.2f | beta: %5.2f\n", v, j, d, alpha[v][j][d], beta[v][j][d]);*/
	      /*printf("post[%d][%d][%d]: %f\n", cm->M, j, d, post[cm->M][j][d]);*/
	    }
	}
  /* If the caller doesn't want the matrix, free it and free the decks in the pool
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= cm->M; v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { free_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { free_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  }
  else *ret_alpha = alpha;

  /* If the caller doesn't want the beta matrix, free it along with the decks.
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= cm->M; v++) 
      if (beta[v] != NULL) { free_vjd_deck(beta[v], 1, L); beta[v] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the post matrix, free it, though
   * it would be *stupid* for the caller not to want it in current implementation.
   */
  if (ret_post == NULL) {
    for (v = 0; v <= cm->M; v++) 
      if (post[v] != NULL) { free_vjd_deck(post[v], 1, L); post[v] = NULL; }
    free(post);
  } else *ret_post = post;

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


/* Function: CMPostalCode()
 * Date:     EPN 05.25.06 based on SRE's PostalCode() 
 *           from HMMER's postprob.c
 *
 * Purpose:  Given a parse tree and a posterior
 *           probability cube, calculate a string that
 *           represents the confidence values on each 
 *           residue in the sequence.
 *           
 *           The code string is 0..L-1  (L = len of target seq),
 *           so it's in the coordinate system of the sequence string;
 *           off by one from dsq; and convertible to the coordinate
 *           system of aseq using MakeAlignedString().
 *           
 *           Values are 0-9,*  
 *           for example, 9 means with >=90% posterior probabiility,
 *           residue i is aligned to the state k that it
 *           is assigned to in the given trace.
 *
 * Args:     L    - length of seq
 *           post - posterior prob cube: see CMPosterior()
 *           *tr  - parsetree to get a Postal code string for.   
 *
 * Returns:  char * array of codes, 0..L-1
 *           Caller is responsible for free'ing it.
 */
char
Fscore2postcode(float sc)
{
  char i;
  /*printf("sc: %10.2f\n", sc);*/
  i = (char) (FScore2Prob(sc, 1.) * 10.);
  return ((i > 9) ? '*' : '0'+i);
}
char *
CMPostalCode(CM_t *cm, int L, float ***post, Parsetree_t *tr)
{
  int x, v, i, j, d, r;
  char *postcode;

  postcode = MallocOrDie((L+1) * sizeof(char)); 

  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];
      d = j-i+1;
      /*printf("x: %2d | v: %2d | i: %2d | j: %2d | d: %2d | post[%d][%d][%d]: %f\n", x, v, i, j, d, v, j, d, post[v][j][d]);*/
      /*
       * Only P, L, R states have emissions.
       */
      if (cm->sttype[v] == MP_st) {
	postcode[i-1] = Fscore2postcode(post[v][j][d]);
	postcode[j-1] = Fscore2postcode(post[v][j][d]);
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	postcode[i-1] = Fscore2postcode(post[v][j][d]);
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	postcode[j-1] = Fscore2postcode(post[v][j][d]);
      } else if (cm->sttype[v] == EL_st) /*special case*/ {
	for(r = (i-1); r <= (j-1); r++)
	  {
	    d = j - (r+1) + 1;
	    postcode[r] = Fscore2postcode(post[v][j][d]);
	    /*printf("r: %d | post[%d][%d][%d]: %f | sc: %c\n", r, v, j, d, post[v][j][d], postcode[r]);*/
	  }
      }
    }
  postcode[L] = '\0';
  return(postcode);
}

/***************************************************************/
/* Function: CMCheckPosterior()
 *           
 * EPN 05.25.06 
 *
 * Purpose:  Given a posterior probability cube, check to make
 *           sure that for each residue k of the sequence:
 *           \sum_v p(v | k emitted from v) = 1.0
 *           To check this, we have to allow possibility that 
 *           the res at posn k was emitted from a left 
 *           emitter or a right emitter.
 *           
 * Args:     L        - length of sequence
 *           cm       - the model
 *           post     - pre-allocated dynamic programming cube
 *           
 * Return:   void
 */
void
CMCheckPosterior(int L, CM_t *cm, float ***post)
{
  int   v, j, d, k;
  float sc;

  for (k = 1; k <= L; k++) 
    {
      sc = IMPOSSIBLE;
      for (v = (cm->M - 1); v >= 0; v--) 
	{
	  if((cm->sttype[v] == MP_st) ||
	     (cm->sttype[v] == ML_st) ||
	     (cm->sttype[v] == IL_st))
	    {
	      for (j = k; j <= L; j++)
		{
		  /*printf("adding L v: %d | i: %d | j: %d | d: %d\n", v, (j-d+1), j, d);*/
		  d = j-k+1;
		  sc = LogSum2(sc, (post[v][j][d]));
		}
	    }
	  if((cm->sttype[v] == MP_st) ||
	     (cm->sttype[v] == MR_st) ||
	     (cm->sttype[v] == IR_st))
	    {
	      for (d = 1; d <= k; d++)
		{
		  /*printf("adding R v: %d | i: %d | j: %d | d: %d\n", v, (k-d+1), k, d);*/
		  sc = LogSum2(sc, (post[v][k][d]));
		}
	    }
	}
      /* Finally factor in possibility of a local end, i.e. that the EL state
       * may have "emitted" this residue.
       */
      if (cm->flags & CM_LOCAL_END) {
	for (j = k; j <= L; j++)
	  {
	    d = j-k+1;
	    /*printf("EL adding L v: %d | i: %d | j: %d | d: %d post[v][j][d]: %5.2f\n", cm->M, (j-d+1), j, d, post[cm->M][j][d]);*/
	    sc = LogSum2(sc, (post[cm->M][j][d]));
	  }
      }
      
      if(((sc - 0.) > 0.0001) || ((sc - 0.) < -0.0001))
	{
	  Die("residue position %d has summed prob of %5.4f (2^%5.4f) in posterior cube.\n", k, (sreEXP2(sc)), sc);
	}
      /*printf("k: %d | total: %10.2f\n", k, (sreEXP2(sc)));*/
    }  
}

/*****************************************************************
 * CM FInside_b_jd_me() & FOutside_b_jd_me() functions.
 * Banded versions of FInside() and FOutside() that only 
 * allocate cells within state dependent j bands, and 
 * state & j dependent d bands.
 *****************************************************************/ 

/* Function: FInside_b_jd_me()
 * EPN 05.26.06
 * 
 * Purpose:  The banded Inside dynamic programming algorithm for CMs.
 *           Banded in both the j and d dimensions, and works
 *           with transformed coordinates for memory efficiency, 
 *           only alpha cells within bands are allocated.
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
 *           subsequence from i0..j0, using the entire model, enforcing
 *           bands in both the j and d dimensions.
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
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *
 * Returns: log P(S|M)/P(S|R), as a bit score.
 */
float 
FInside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
		float ***alpha, float ****ret_alpha, 
		struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		int allow_begin, int *jmin, int *jmax, int **hdmin, int **hdmax)
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
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;

  if(i0 != 1)
    {
      printf("FInside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");
      exit(1);
    }
  if(j0 != L)
    {
      printf("FInside_b_jd requires that j0 be L. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
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

  /* EPN: now that the EL self loop has a transition score, its
   *      necessary to keep track of the alpha EL deck (alpha[cm->M]).
   *      There's no bands on the EL state. 
   */
  if(cm->flags & CM_LOCAL_BEGIN)
    {
      if (! deckpool_pop(dpool, &(alpha[cm->M]))) 
	alpha[cm->M] = alloc_vjd_deck(L, i0, j0);
      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	/*alpha[cm->M][j][0] = IMPOSSIBLE;*/
	alpha[cm->M][j][0] = 0.;
	for (d = 1; d <= jp; d++)
	  {
	    alpha[cm->M][j][d] = (cm->el_selfsc * (d));
	  }
      }
    }

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
	alpha[v] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);

      /* We've only allocated alpha cells that are within the bands
       * on the j and d dimensions. This means we have to deal
       * with all sorts of offset issues, but we don't have to 
       * waste time setting cells outside the bands to IMPOSSIBLE.
       */

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		  alpha[v][jp_v][dp_v]  = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  /* treat EL as emitting only on self transition */
		  for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		    {
		      yoffset = y - cm->cfirst[v];
		      if(j >= jmin[y] && j <= jmax[y]) 
			/* Enforces j is valid for state y */
			{
			  jp_y = j - jmin[y];
			  if(d >= hdmin[y][jp_y] && d <= hdmax[y][jp_y])
			    {
			      dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							     in alpha w/mem eff bands */
			      /* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			       * corresponding to alpha[y][j][d] in the platonic matrix.
			       */
			      alpha[v][jp_v][dp_v] = LogSum2(alpha[v][jp_v][dp_v], (alpha[y][jp_y][dp_y] 
										    + cm->tsc[v][yoffset]));
			    }
			}
		    }
		  if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
		}
	    }
	}
      else if (cm->sttype[v] == B_st)
	{
	  y = cm->cfirst[v];
	  z = cm->cnum[v];
	  /* Any valid j must be within both state v and state z's j band 
	   * I think jmin[v] <= jmin[z] is guaranteed by the way bands are 
	   * constructed, but we'll check anyway. 
	   */
	  tmp_jmin = (jmin[v] > jmin[z]) ? jmin[v] : jmin[z];
	  tmp_jmax = (jmax[v] < jmax[z]) ? jmax[v] : jmax[z];

	  /* For any values of j within v's j band but outside of z's j band,
	   * we have to set the corresponding alpha cells to IMPOSSIBLE.
	   * This is done be the following two ugly for loops, 
	   * which will only be looked at once for each B state, and
	   * even then only *very* rarely entered. This
	   * is why they're here, seemingly out of place before the 
	   * main j loop below, where similar performing code would be 
	   * looked at on the order of j times, instead of just once.
	   */
	  for(j = jmin[v]; j < tmp_jmin; j++)
	    {
	      jp_v = j-jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d-hdmin[v][jp_v];
		  alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
		}
	    }
	  if(tmp_jmax < jmax[v])
	    for(j = (tmp_jmax+1); j <= jmax[v]; j++)
	      {
		jp_v = j-jmin[v];
		for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		  {
		    dp_v = d-hdmin[v][jp_v];
		    alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* this won't be changed */
		  }
	      }
	  /* the main j loop */
	  for (j = tmp_jmin; j <= tmp_jmax; j++)
	    {
	      jp_v = j - jmin[v];
	      jp_y = j - jmin[y];
	      jp_z = j - jmin[z];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */

		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 6 inequalities (some may be redundant):
		   * (1) k >= j-jmax[y];
		   * (2) k <= j-jmin[y]; 
		   *     1 and 2 guarantee (j-k) is within state y's j band
		   *
		   * (3) k >= hdmin[z][j-jmin[z]];
		   * (4) k <= hdmax[z][j-jmin[z]]; 
		   *     3 and 4 guarantee k is within z's j=(j), d band
		   *
		   * (5) k >= d-hdmax[y][j-jmin[y]-k];
		   * (6) k <= d-hdmin[y][j-jmin[y]-k]; 
		   *     5 and 6 guarantee (d-k) is within state y's j=(j-k) d band
		   */
		  kmin = ((j-jmax[y]) > (hdmin[z][jp_z])) ? (j-jmax[y]) : hdmin[z][jp_z];
		  /* kmin satisfies inequalities (1) and (3) */
		  kmax = ( jp_y       < (hdmax[z][jp_z])) ?  jp_y       : hdmax[z][jp_z];
		  /* kmax satisfies inequalities (2) and (4) */
		  /* RHS of inequalities 5 and 6 are dependent on k, so we check
		   * for these within the next for loop.
		   */
		  alpha[v][jp_v][dp_v] = IMPOSSIBLE; /* initialize */
		  for(k = kmin; k <= kmax; k++)
		    {
		      if((k >= d - hdmax[y][jp_y-k]) && k <= d - hdmin[y][jp_y-k])
			{
			  /* for current k, all 6 inequalities have been satisified 
			   * so we know the cells corresponding to the platonic 
			   * matrix cells alpha[v][j][d], alpha[y][j-k][d-k], and
			   * alpha[z][j][k] are all within the bands. These
			   * cells correspond to alpha[v][jp_v][dp_v], 
			   * alpha[y][jp_y-k][d-hdmin[y][jp_y-k]-k],
			   * and alpha[z][jp_z][k-hdmin[z][jp_z]];
			   */
			  kp_z = k-hdmin[z][jp_z];
			  dp_y = d-hdmin[y][jp_y-k];

			  alpha[v][jp_v][dp_v] = LogSum2(alpha[v][jp_v][dp_v], (alpha[y][jp_y-k][dp_y-k] 
										+ alpha[z][jp_z][kp_z]));
			}
		    }
		  if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
		  /* CYK Full ME Bands used 5 end block */
		}
	    }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces (j-1) is valid for state y */
		      {
			jp_y = j - jmin[y];
			if((d-2) >= hdmin[y][jp_y-1] && (d-2) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-2] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-2] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = LogSum2(alpha[v][jp_v][dp_v], (alpha[y][jp_y-1][dp_y-2] 
										  + cm->tsc[v][yoffset]));
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		else
		  alpha[v][jp_v][dp_v] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);
		
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if(j >= jmin[y] && j <= jmax[y]) /* Enforces j is valid for state y */
		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y] && (d-1) <= hdmax[y][jp_y])
			  {
			    dp_y = d - hdmin[y][jp_y];  /* d index for state y 
							   in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j][d-1] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = LogSum2(alpha[v][jp_v][dp_v], (alpha[y][jp_y][dp_y-1] 
										  + cm->tsc[v][yoffset]));
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < Alphabet_size)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][jp_v][dp_v] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		alpha[v][jp_v][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		for (y = cm->cfirst[v]; y < (cm->cfirst[v] + cm->cnum[v]); y++) 
		  {
		    yoffset = y - cm->cfirst[v];
		    if((j-1) >= jmin[y] && (j-1) <= jmax[y]) /* Enforces j-1 is valid for state y */
		      
		      {
			jp_y = j - jmin[y];
			if((d-1) >= hdmin[y][jp_y-1] && (d-1) <= hdmax[y][jp_y-1])
			  {
			    dp_y = d - hdmin[y][jp_y-1];  /* d index for state y 
							     in alpha w/mem eff bands */
			    /* if we get here alpha[y][jp_y-1][dp_y-1] is a valid alpha cell
			     * corresponding to alpha[y][j-1][d-1] in the platonic matrix.
			     */
			    alpha[v][jp_v][dp_v] = LogSum2(alpha[v][jp_v][dp_v], (alpha[y][jp_y-1][dp_y-1] 
										  + cm->tsc[v][yoffset]));
			  }
		      }
		  }
		if (dsq[j] < Alphabet_size)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][jp_v][dp_v] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}				/* finished calculating deck v. */
      
      /* Keep track of contributions of local begins */
      /* The following loops originally access alpha[v][j0][W] but the index W will be
	 in different positions due to the bands */
      if(j0 >= jmin[v] && j0 <= jmax[v])
	{
	  jp_v = j0 - jmin[v];
	  if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v])
	    {
	      Wp = W - hdmin[v][jp_v];
	      /* If we get here alpha[v][jp_v][Wp] is a valid cell
	       * in the banded alpha matrix, corresponding to 
	       * alpha[v][j0][W] in the platonic matrix.
	       */
	      /* Check for local begin getting us to the root.
	       * This is "off-shadow": if/when we trace back, we'll handle this
	       * case separately (and we'll know to do it because we'll immediately
	       * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	       * to jump right to state b; see below)
	       */
	      if (allow_begin)
		{
		  bsc = LogSum2(bsc, (alpha[v][jp_v][Wp] + cm->beginsc[v]));
		}
	    }
	}

      /* If we're at the root state, record contribution of local begins */
      if (allow_begin && v == 0)
	{
	  if(j0 >= jmin[0] && j0 <= jmax[0])
	    {
	      jp_v = j0 - jmin[v];
	      if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v])
		{
		  if (allow_begin && v == 0)
		    alpha[0][jp_v][Wp] = LogSum2(alpha[0][jp_v][Wp], bsc);
		}
	    }
	}	  

      /* Now, if we're trying to reuse memory in our normal mode (e.g. ! do_full):
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { 
	    /* Original code block : */
	    /* we can definitely release the S children of a bifurc. 
	       y = cm->cfirst[v]; deckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	       z = cm->cnum[v];   deckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	     End of original code block */
	    /* New ME code : */
	    y = cm->cfirst[v]; deckpool_push(dpool, alpha[y]); alpha[y] = NULL;
	    z = cm->cnum[v];   deckpool_push(dpool, alpha[z]); alpha[z] = NULL;
	    free_vjd_deck(alpha[y], i0, j0);
	    alpha[y] = NULL;
	    free_vjd_deck(alpha[z], i0, j0);
	    alpha[z] = NULL;
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
		      /* Original code : if (nends == 0) { deckpool_push(dpool, end); end = NULL;} */
		      /* ME code deletes the previous line, we don't mess with end, because
			 it is used later */
		    } else 
		      /* original code (deck reuse) deckpool_push(dpool, alpha[y]);*/
		      /* new ME code : */
		      {
			//printf("calling free vjd deck for alpha[y=%d]\n", y);
			free_vjd_deck(alpha[y], i0, j0);
		      }
		      alpha[y] = NULL;
		  }
	      }
	  }
      }
  } /* end loop over all v */

  /*
    printf("INSIDE JD, printing alpha\n");
    debug_print_alpha_banded_jd(alpha, cm, L, jmin, jmax, hdmin, hdmax);
    printf("INSIDE JD, done printing alpha\n");
  */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  Wp = W - hdmin[0][j0-jmin[0]];
  sc = alpha[0][j0-jmin[0]][Wp];

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
  printf("***returning from FInside_b_jd_me() sc  : %f\n", sc);
  return sc;
}

/***********************************************************************
 * Function: FOutside_b_jd_me()
 * Date:     EPN 05.26.06
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
 *           do_full   - if TRUE, we save all the decks in alpha and beta, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (vroot) is valid
 *                       at the end.
 *           beta       - if non-NULL, this is an existing matrix, with NULL
 *                       decks for vroot..vend, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_beta  - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated decks sized
 *                       for this subsequence i0..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..j0 subseq,
 *                       because of the size of the subseq decks.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           alpha     - the alpha matrix from a Inside run, if do_check is FALSE
 *                        only decks for S states must be valid, else all must be
 *                        valid.
 *           ret_alpha - if non-NULL, return the alpha matrix with one or more
 *                       decks available for examination (see "do_full")
 *           do_check  - TRUE to do time-consuming check to make sure
 *                       beta and alpha are consistent (only if NON-LOCAL mode)
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * 
 * Returns: log P(S|M)/P(S|R), as a bit score.
 ***********************************************************************/
float 
FOutside_b_jd_me(CM_t *cm, char *dsq, int L, int i0, int j0, int do_full,
		 float ***beta, float ****ret_beta, 
		 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		 int allow_begin, float ***alpha, float ****ret_alpha, 
		 int do_check, int *jmin, int *jmax, int **hdmin, int **hdmax)
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
  float    return_sc;   /* P(S|M)/P(S|R) */
  int      n;           /* counter over nodes, used only if do_check = TRUE */
  int      num_split_states; /* temp variable used only if do_check = TRUE */
  float    diff;        /* temp variable used only if do_check = TRUE */
  float  **end;         /* we re-use the end deck. */
  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;

  if(i0 != 1)
    {
      printf("FOutside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");
      exit(1);
    }
  if(j0 != L)
    {
      printf("FOutside_b_jd requires that j0 be L. This function is not set up for subsequence alignment.\n");
      exit(1);
    }

  if (cm->flags & CM_LOCAL_END) { do_check = FALSE; } 
  /* Code for checking doesn't apply in local mode. See below. */

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
    beta[0] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[0], jmax[0], hdmin[0], hdmax[0]);
  for (j = jmin[0]; j <= jmax[0]; j++)
    {
      jp_v = j - jmin[0];
      for (d = hdmin[0][jp_v]; d <= hdmax[0][jp_v]; d++)
	{
	  dp_v = d - hdmin[0][jp_v];  /* d index for state v in alpha w/mem eff bands */
	  beta[0][jp_v][dp_v] = IMPOSSIBLE;
	}
    }
  /* non banded line: beta[0][j0][W] = 0; */
  jp_v = j0 - jmin[0];
  Wp = W - hdmin[0][jp_v];
  assert(W >= hdmin[0][jp_v]);
  assert(W <= hdmax[0][jp_v]);
  beta[0][jp_v][Wp] = 0.;

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   * EL deck has no bands as currently implemented.
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
	beta[v] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);

      /* Init the whole deck to IMPOSSIBLE
       */
      for (j = jmin[v]; j <= jmax[v]; j++)
	{
	  jp_v = j - jmin[v];
	  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	    {
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      beta[v][jp_v][dp_v] = IMPOSSIBLE;
	    }
	}

      /* If we can do a local begin into v, also init with that. 
       * By definition, beta[0][j0][W] == 0.
       */ 
      if (i0 == 1 && j0 == L && (cm->flags & CM_LOCAL_BEGIN))
	{
	  if((j0 >= jmin[v]) && (j0 <= jmax[v]))
	    {
	      jp_v = j0 - jmin[v];
	      if((W >= hdmin[v][jp_v]) && W <= hdmax[v][jp_v])
		{
		  Wp = W - hdmin[v][jp_v];
		  beta[v][jp_v][Wp] = cm->beginsc[v];
		}
	    }
	}
      /* main recursion: reorganized relative to FOutside() for simplification of
       * band-related issues.
       */
      for (j = jmax[v]; j >= jmin[v]; j--)
	{
	  jp_v = j - jmin[v];
	  for (d = hdmax[v][jp_v]; d >= hdmin[v][jp_v]; d--)
	    {
	      i = j-d+1;
	      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
	      
	      if (cm->stid[v] == BEGL_S) 
		{
		  y = cm->plast[v];	/* the parent bifurcation    */
		  z = cm->cnum[y];	/* the other (right) S state */
		  jp_y = j - jmin[y];
		  jp_z = j - jmin[z];
		  
		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 8 inequalities (some may be redundant):
		   * NOTE: these are different from those in Inside() (for one thing, v and y
		   *       (BEGL_S and BIF_B here respectively) are switched relative to Inside.
		   *
		   * (1) k <= jmax[y] - j;
		   * (2) k >= jmin[y] - j;
		   * (3) k <= jmax[z] - j;
		   * (4) k >= jmin[z] - j;
		   *     1 and 2 guarantee (j+k) is within state y's j band
		   *     3 and 4 guarantee (j+k) is within state z's j band
		   *
		   * (5) k >= hdmin[y][j-jmin[y]+k] - d;
		   * (6) k <= hdmax[y][j-jmin[y]+k] - d; 
		   *     5 and 6 guarantee k+d is within y's j=(j+k), d band
		   *
		   * (7) k >= hdmin[z][j-jmin[z]+k];
		   * (8) k <= hdmax[z][j-jmin[z]+k]; 
		   *     5 and 6 guarantee k is within state z's j=(j+k) d band
		   * 
		   * Note, below:
		   * tmp_jmin = MAX(jmin[y], jmin[z];
		   * tmp_jmax = MIN(jmax[y], jmax[z];
		   */
		  tmp_jmin = (jmin[y] > jmin[z]) ? jmin[y] : jmin[z];
		  tmp_jmax = (jmax[y] < jmax[z]) ? jmax[y] : jmax[z];

		  kmin = tmp_jmin - j;
		  kmax = tmp_jmax - j;
		  /* kmin and kmax satisfy inequalities (1-4) */
		  /* RHS of inequalities 5-8 are dependent on k, so we check
		   * for these within the next for loop.
		   */
		  for(k = kmin; k <= kmax; k++)
		    {
		      if(k < (hdmin[y][jp_y+k] - d) || k > (hdmax[y][jp_y+k] - d)) continue; 
		      /* above line continues if inequality 5 or 6 is violated */
		      if(k < (hdmin[z][jp_z+k]) || k > (hdmax[z][jp_z+k])) continue; 
		      /* above line continues if inequality 7 or 8 is violated */
		      
		      /* if we get here for current k, all 8 inequalities have been satisified 
		       * so we know the cells corresponding to the platonic 
		       * matrix cells alpha[v][j][d], alpha[y][j+k][d+k], and
		       * alpha[z][j+k][k] are all within the bands. These
		       * cells correspond to beta[v][jp_v][dp_v], 
		       * beta[y][jp_y+k][d-hdmin[y][jp_y+k]+k],
		       * and alpha[z][jp_z][k-hdmin[z][jp_z+k]];
		       */
		      kp_z = k-hdmin[z][jp_z+k];
		      dp_y = d-hdmin[y][jp_y+k];
		      beta[v][jp_v][dp_v] = LogSum2(beta[v][jp_v][dp_v], (beta[y][jp_y+k][dp_y+k] 
									  + alpha[z][jp_z+k][kp_z]));
		    }
		}
	      else if (cm->stid[v] == BEGR_S) 
		{
		  y = cm->plast[v];	/* the parent bifurcation    */
		  z = cm->cfirst[y];	/* the other (left) S state */

		  jp_y = j - jmin[y];
		  jp_z = j - jmin[z];
		  
		  /* For j to be valid for state y: *
		   * jmin[y] >= j >= jmax[y]
		   * These are independent of k so we check outside of k loop below
		   * For j to be valid for state z: *
		   * (jmin[z] + d) >= j >= (jmax[z] + d)
		   */
		  if(j < jmin[y] || j > jmax[y]) continue;
		  if((j < (jmin[z] + d)) || (j > (jmax[z]+d))) continue;

		  /* Find the first k value that implies a valid cell in the y and z decks.
		   * This k must satisfy the following 4 inequalities (some may be redundant):
		   * NOTE: these are different from those in Inside() (for one thing, v and y
		   *       (BEGR_S and BIF_B here respectively) are switched relative to Inside.
		   *
		   * (1) k >= hdmin[y][j-jmin[y]] - d;
		   * (2) k <= hdmax[y][j-jmin[y]] - d;
		   *     1 and 2 guarantee (d+k) is within state y's j=(j) d band
		   *
		   * (3) k >= hdmin[z][j-jmin[z]-d];
		   * (4) k <= hdmax[z][j-jmin[z]-d];
		   *     3 and 4 guarantee k is within z's j=(j-d), d band
		   *
		   */
		  kmin = ((hdmin[y][jp_y]-d) > (hdmin[z][jp_z-d])) ? (hdmin[y][jp_y]-d) : (hdmin[z][jp_z-d]);
		  /* kmin satisfies inequalities (1) and (3) */
		  kmax = ((hdmax[y][jp_y]-d) < (hdmax[z][jp_z-d])) ? (hdmax[y][jp_y]-d) : (hdmax[z][jp_z-d]);
		  /* kmax satisfies inequalities (2) and (4) */

		  for(k = kmin; k <= kmax; k++)
		    {
		      /* for current k, all 4 inequalities have been satisified 
		       * so we know the cells corresponding to the platonic 
		       * matrix cells beta[v][j][d], beta[y][j][d+k], and
		       * alpha[z][j-d][k] are all within the bands. These
		       * cells correspond to beta[v][jp_v][dp_v], 
		       * beta[y][jp_y+k][d-hdmin[y][jp_y]+k],
		       * and alpha[z][jp_z-d][k-hdmin[z][jp_z-d]];
		       */
		      kp_z = k-hdmin[z][jp_z-d];
		      dp_y = d-hdmin[y][jp_y];
		      beta[v][jp_v][dp_v] = LogSum2(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+k] 
									  + alpha[z][jp_z-d][kp_z]));
		    }
		}
	      else
		{
		  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) 
		    {
		      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

		      switch(cm->sttype[y]) {
		      case MP_st: 
			if (j == j0 || d == j) continue; /* boundary condition */
			if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+2) < hdmin[y][(jp_y+1)] || (d+2) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y+1][dp_y+2] is a valid alpha cell
			 * corresponding to alpha[y][j+1][d+2] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y+1];  /* d index for state y */
			
			if (dsq[i-1] < Alphabet_size && dsq[j+1] > Alphabet_size)
			  escore = cm->esc[y][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
			else
			  escore = DegeneratePairScore(cm->esc[y], dsq[i-1], dsq[j+1]);
			beta[v][jp_v][dp_v] = LogSum2(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+2] 
									    + cm->tsc[y][voffset] + escore));
			break;

		      case ML_st:
		      case IL_st: 
			if (d == j) continue;	/* boundary condition (note when j=0, d=0)*/
			if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+1) < hdmin[y][jp_y] || (d+1) > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y][dp_y+1] is a valid alpha cell
			 * corresponding to alpha[y][j][d+1] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y];  /* d index for state y */
			if (dsq[i-1] < Alphabet_size) 
			  escore = cm->esc[y][(int) dsq[i-1]];
			else
			  escore = DegenerateSingletScore(cm->esc[y], dsq[i-1]);
			beta[v][jp_v][dp_v] = LogSum2(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y+1] 
									    + cm->tsc[y][voffset] + escore));
			break;
		    
		      case MR_st:
		      case IR_st:
			if (j == j0) continue;
			if ((j+1) < jmin[y] || (j+1) > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if ((d+1) < hdmin[y][(jp_y+1)] || (d+1) > hdmax[y][(jp_y+1)]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y+1][dp_y+1] is a valid alpha cell
			 * corresponding to alpha[y][j+1][d+1] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][(jp_y+1)];  /* d index for state y */
			if (dsq[j+1] < Alphabet_size) 
			  escore = cm->esc[y][(int) dsq[j+1]];
			else
			  escore = DegenerateSingletScore(cm->esc[y], dsq[j+1]);
			/*printf("j: %d | jmin[y]: %d | jmax[y]: %d | jp_v: %d | dp_v: %d | jp_y: %d | dp_y: %d\n", j, jmin[y], jmax[y], jp_v, dp_v, jp_y, dp_y);*/
			beta[v][jp_v][dp_v] = LogSum2(beta[v][jp_v][dp_v], (beta[y][jp_y+1][dp_y+1] 
									    + cm->tsc[y][voffset] + escore));
			break;

		      case S_st:
		      case E_st:
		      case D_st:
			if (j < jmin[y] || j > jmax[y]) continue; /* enforces j is valid for state y */
			jp_y = j - jmin[y];
			if (d < hdmin[y][jp_y] || d > hdmax[y][jp_y]) continue; /* enforces d is valid for state y */
			/* if we get here alpha[y][jp_y][dp_y] is a valid alpha cell
			 * corresponding to alpha[y][j][d] in the platonic matrix.
			 */
			dp_y = d - hdmin[y][jp_y];  /* d index for state y */
			beta[v][jp_v][dp_v] = LogSum2(beta[v][jp_v][dp_v], (beta[y][jp_y][dp_y] + cm->tsc[y][voffset])); 
			break;

		      default: Die("bogus child state %d\n", cm->sttype[y]);
		      }/* end switch over states*/
		  } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
		if (beta[v][jp_v][dp_v] < IMPOSSIBLE) beta[v][jp_v][dp_v] = IMPOSSIBLE;
		}	/* ends else entered for non-BEGL/BEGR states*/	
	    } /* ends loop over d. We know all beta[v][j][d] in this row j*/
	}    /* end loop over jp. We know the beta's for the whole deck.*/
      /* Deal with local alignment end transitions v->EL
       * (EL = deck at M.)
       */
      if (NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (j = jmin[v]; j <= jmax[v]; j++)
	  {
	    jp_v = j - jmin[v];
	    for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	      {
		i = j-d+1;
		dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		switch (cm->sttype[v]) {
		case MP_st: 
		  if (j == j0 || d == j) continue; /* boundary condition */
		  if (((j+1) > jmax[v]) || ((d+2) > hdmax[v][jp_v])) continue; /*boundary condition*/
		  if (dsq[i-1] < Alphabet_size && dsq[j+1] > Alphabet_size)
		    escore = cm->esc[v][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
		  else
		    escore = DegeneratePairScore(cm->esc[v], dsq[i-1], dsq[j+1]);
		  beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+2] + cm->endsc[v] 
								  + escore));
		break;
	      case ML_st:
	      case IL_st:
		if (d == j) continue;	
		if ((d+1) > hdmax[v][jp_v]) continue; /*boundary condition*/
		if (dsq[i-1] < Alphabet_size) 
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = DegenerateSingletScore(cm->esc[v], dsq[i-1]);
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][jp_v][dp_v+1] + cm->endsc[v] 
								+ escore));
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (((j+1) > jmax[v]) || ((d+1) > hdmax[v][jp_v])) continue; /*boundary condition*/
		if (dsq[j+1] < Alphabet_size) 
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = DegenerateSingletScore(cm->esc[v], dsq[j+1]);
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][jp_v+1][dp_v+1] + cm->endsc[v] 
								+ escore));
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[v][jp_v][dp_v] + cm->endsc[v] 
								+ escore));
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
  /*#if 0*/
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CM_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	beta[cm->M][j][d] = LogSum2(beta[cm->M][j][d], (beta[cm->M][j][d+1]
							+ cm->el_selfsc));
    }
  }
  /*#endif*/

  /*
    printf("OUTSIDE JD, printing beta\n");
    debug_print_alpha_banded_jd(beta, cm, L, jmin, jmax, hdmin, hdmax);
    printf("OUTSIDE JD, done printing beta\n");
  */

  Wp = W - hdmin[0][j0-jmin[0]];
  if(do_check && (!(cm->flags & CM_LOCAL_END))) 
    /* Local ends make the following test invalid because it is not true that
     * exactly 1 state in each node's split set must be visited in each parse. 
     */
    {
      /* Determine P(S|M) / P(S|R) (probability of the sequence given the model) 
       * using both the Outside (beta) and Inside (alpha) matrices,
       * and ensure they're consistent with P(S|M) / P(S|R) from the Inside calculation.
       * For all v in each split set: Sum_v [ Sum_i,j ( alpha[v][i][j] * beta[v][i][j] ) ] 
       *                                                = P(S|M) / P(S|R)  
       * in v,j,d coordinates this is:
       * For all v in each split set: Sum_v [ Sum_j,(d<=j) ( alpha[v][j][d] * beta[v][j][d] ) ]
       *                                                = P(S|M) / P(S|R)
       */
	 
      for(n = 0; n < cm->nodes; n++)
	{
	  sc = IMPOSSIBLE;
	  num_split_states = SplitStatesInNode(cm->ndtype[n]);
	  for(v = cm->nodemap[n]; v < cm->nodemap[n] + num_split_states; v++)
	    {
	      for (j = jmin[v]; j <= jmax[v]; j++)
		{
		  jp_v = j - jmin[v];
		  for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		    {
		      dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		      /*printf("node %d | adding alpha beta: v: %d | jp_v: %d | dp_v: %d| j: %d | d: %d\n", n, v, jp_v, dp_v, j, d);
			printf("\talpha: %f | beta: %f\n", alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
		      sc = LogSum2(sc, (alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v]));
		    }
		}
	      }
	  printf("checking node: %d | sc: %.6f\n", n, sc);
	  diff = alpha[0][j0-jmin[0]][Wp] - sc;
	  if(diff > 0.001 || diff < -0.001)
	    {
	      Die("ERROR: node %d P(S|M): %.5f inconsistent with Inside P(S|M): %.5f (diff: %.5f)\n", 
		  n, sc, alpha[0][(j0-jmin[0])][Wp], diff);
	    }
	}

    }

  /* IF not in local mode, we can calculate P(S|M) / P(S|R) given only the 
   * beta matrix as follows:
   * 
   * IF local ends are off, we know each parse MUST visit each END_E state,
   * we pick final END_E state state cm->M-1 (though any END_E could be used here):
   *
   * Sum_j=0 to L (alpha[M-1][j][0] * beta[M-1][j][0]) = P(S|M) / P(S|R)
   *
   * Note: alpha[M-1][j][0] = 0.0 for all j 
   *       because all parse subtrees rooted at an END_E must have d=0, (2^0 = 1.0)
   * therefore: 
   * Sum_j=0 to L (beta[M-1][j][0]) = P(S|M) / P(S|R)
   * 
   * IF local ends are on, each parse MUST visit either each END_E state with d=0
   * or the EL state but d can vary, so we can't use this test (believe me I tried
   * to get a similar test working, but I'm convinced you need alpha to get P(S|M)
   * in local mode).
   */
  if(!(cm->flags & CM_LOCAL_END))
    {
      return_sc = IMPOSSIBLE;
      v = cm->M-1;
      for (j = jmin[v]; j <= jmax[v]; j++)
	{
	  jp_v = j - jmin[v];
	  assert(hdmin[v][jp_v] == 0);
	  /* printf("\talpha[%3d][%3d][%3d]: %5.2f | beta[%3d][%3d][%3d]: %5.2f\n", (cm->M-1), (j), 0, alpha[(cm->M-1)][j][0], (cm->M-1), (j), 0, beta[(cm->M-1)][j][0]);*/
	  return_sc = LogSum2(return_sc, (beta[v][jp_v][0]));
	}
    }
  else /* return_sc = P(S|M) / P(S|R) from Inside() */
    {
      return_sc = alpha[0][(j0-jmin[0])][Wp];
    }

  /* If the caller doesn't want the beta matrix, free it.
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

  /* If the caller doesn't want the alpha matrix, free it 
   * Else, pass it back to him.
   * EPN - if we free the alpha and beta matrix the deck pool has all the 
   *       decks from alpha and beta, not sure if this is desirable.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M-1); v++) 
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
    float **a;
    while (deckpool_pop(dpool, &a)) free_vjd_deck(a, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }
  free(touch);

  if(!(cm->flags & CM_LOCAL_END))
    printf("***returning from FOutside_b_jd_me() sc : %f\n", return_sc);
  else
    printf("***returning from FOutside_b_jd_me() sc : %f (LOCAL mode; sc is from Inside)\n", return_sc);
  return return_sc;
}

/***************************************************************/
/* Function: CMPosterior_b_jd_me()
 *           
 * EPN 05.27.06 based on IHH's P7EmitterPosterior() from HMMER's postprob.c
 *
 * Purpose:  Combines banded Inside and Outside cubes into a posterior
 *           probability cube. The Inside and Outside cubes are banded
 *           in both the j and d dimensions, any cells outside of
 *           bands do not exist in memory
 *           The entry in post[v][jp_v][dp_v] is the log of the
 *           posterior probability of a parse subtree rooted at v 
 *           emitting the subsequence i..j (i=j-d+1). Where j = jp_v + jmin[v],
 *           and d = dp_v + hdmin[v][jp_v].
 *           The caller must allocate space for the cube, although the
 *           beta matrix from Outside can be used instead (overwriting it will not
 *           compromise the algorithm).
 *           
 * Args:     L        - length of sequence
 *           cm       - the model
 *           alpha    - pre-calculated Inside matrix 
 *           beta     - pre-calculated Outside matrix
 *           post     - pre-allocated dynamic programming cube
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * Return:   void
 */
void
CMPosterior_b_jd_me(int L, CM_t *cm, float ***alpha, float ****ret_alpha,
		    float ***beta, float ****ret_beta, float ***post,
		    float ****ret_post, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int   v, j, d;
  float sc;
  int      jp_v; /* j index for state v in alpha/beta with mem eff bands */
  int      dp_v; /* d index for state v in alpha/beta w/mem eff bands */
  int      Lp;
  float  **end;         /* used for freeing alpha b/c we re-use the end deck. */
  
  Lp = L - hdmin[0][L-jmin[0]];
  sc = alpha[0][L-jmin[0]][Lp];
  
  /* If local ends are on, start with the EL state (cm->M), otherwise
   * its not a valid deck.
   */
  if (cm->flags & CM_LOCAL_END)
    {
      for(j = 0; j <= L; j++) 
	for (d = 0; d <= j; d++)
	  {
	    post[cm->M][j][d] = alpha[cm->M][j][d] + beta[cm->M][j][d] - sc;
	    /*printf("v: %3d | j: %3d | d: %3d | alpha : %5.2f | beta : %5.2f\n", cm->M, j, d, alpha[cm->M][j][d], beta[cm->M][j][d]);*/
	    /*printf("post[%d][%d][%d]: %f\n", cm->M, j, d, post[cm->M][j][d]);*/
	  }  
    }

  for (v = (cm->M-1); v >= 0; v--) 
    for (j = jmin[v]; j <= jmax[v]; j++) 
      {
	jp_v = j - jmin[v];
	for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
	{
	  dp_v = d - hdmin[v][jp_v];
	  /*printf("v: %3d | jp_v: %3d | dp_v: %3d | alpha: %5.2f | beta: %5.2f\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v], beta[v][jp_v][dp_v]);*/
	  post[v][jp_v][dp_v] = alpha[v][jp_v][dp_v] + beta[v][jp_v][dp_v] - sc;
	}  
      }

  /* If the caller doesn't want the matrix, free it and free the decks in the pool
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = 0; v <= (cm->M); v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { free_vjd_deck(alpha[v], 1, L); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { free_vjd_deck(end, 1, L); end = NULL; }
    free(alpha);
  }
  else *ret_alpha = alpha;

  /* If the caller doesn't want the beta matrix, free it along with the decks.
   */
  if (ret_beta == NULL) {
    for (v = 0; v <= (cm->M); v++) 
      if (beta[v] != NULL) { free_vjd_deck(beta[v], 1, L); beta[v] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the post matrix, free it, though
   * it would be *stupid* for the caller not to want it in current implementation.
   */
  if (ret_post == NULL) {
    for (v = 0; v <= (cm->M); v++) 
      if (post[v] != NULL) { free_vjd_deck(post[v], 1, L); post[v] = NULL; }
    free(post);
  } else *ret_post = post;
}

char *
CMPostalCode_b_jd_me(CM_t *cm, int L, float ***post, Parsetree_t *tr,
		    int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int x, v, i, j, d, r;
  char *postcode;
  int jp_v, dp_v;
  postcode = MallocOrDie((L+1) * sizeof(char)); 

  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];
      i = tr->emitl[x];
      j = tr->emitr[x];
      d = j-i+1;
      /*
       * Only P, L, R states have emissions.
       */
      if(cm->sttype[v] != EL_st)
	{
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
      else
	{
	  jp_v = j;
	  dp_v = d;
	}
      if (cm->sttype[v] == MP_st) {
	postcode[i-1] = Fscore2postcode(post[v][jp_v][dp_v]);
	postcode[j-1] = Fscore2postcode(post[v][jp_v][dp_v]);
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	postcode[i-1] = Fscore2postcode(post[v][jp_v][dp_v]);
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	postcode[j-1] = Fscore2postcode(post[v][jp_v][dp_v]);
      } else if (cm->sttype[v] == EL_st) /*special case*/ {
	for(r = (i-1); r <= (j-1); r++)
	  {
	    d = j - (r+1) + 1;
	    postcode[r] = Fscore2postcode(post[v][j][d]);
	    /*printf("r: %d | post[%d][%d][%d]: %f | sc: %c\n", r, v, j, d, post[v][j][d], postcode[r]);*/
	  }
      }
    }
  postcode[L] = '\0';
  return(postcode);
}

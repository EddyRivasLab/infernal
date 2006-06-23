/* hbandcyk.c (created using smallcyk.c as a template EPN 11.04.05)
 * 
 * Alignment of a CM to a target sequence using HMM imposed bands
 * on the j and d dimensions of the CYK matrix.
 * These algorithms align to the entire target sequence
 * (e.g. global alignment). For sequence-local alignment, see
 * scancyk.c.
 *
 *################################################################
 * CYKInside_b_jd()       - Align model to sequence, using normal 
 *                          CYK/Inside algorithm with the CYK
 *                          matrix banded in the j and d dimension.
 *                          The mother of all banded CYK functions.
 *################################################################
 * 
 * Also included are some functions for determining bands,
 * namely: 
 *  ij2d_bands()
 *  hd2safe_hd_bands()
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

#include "stopwatch.h"          /* squid's process timing module        */

/*******************************************************************************
 * 11.04.05
 * EPN 
 * Memory efficient banded versions of selected smallcyk.c functions that
 * enforce bands in the d and j dimensions informed by an HMM forwards/backwards
 * alignment of the target sequence.
 * 
 * These functions are modified from their originals in smallcyk.c to make 
 * FULL (not D&C) CYK alignment somewhat memory efficient. The starting 
 * point for CYKInside_b_jd() was CYKInside_b_me() in smallcyk.c. The main
 * difference is that bands in the j dimension are enforced, and the d
 * bands have j dependence. 
 * The original CYK_Inside_b_me() function was dubbed 'memory efficient' (the _me part)
 * because it only allocates cells of the alpha or shadow matrix which are within the 
 * bands in the d dimension. Here, cells outside the j bands are still 
 * allocated (though they don't need to be).
 * 
 * CYK_Inside_b_jd() only allocates cells within the j AND d bands.
 * 
 * Comments from smallcyk.c pertaining to CYK_Inside_b_me():
 * .........................................................................
 * The only real difficulty implementing memory efficient
 * bands is in being able to determine what cell alpha[v][j][d] from the 
 * non-memory efficient code corresponds to in the memory-efficient code (we'll
 * call the corresponding cell a[v'][j'][d'] or a[vp][jp][dp]).  The reason
 * v != v'; j != j' and d != d' is because the primes are offset due to the
 * fact that some of the original alpha matrix deck (a[v]) has not been allocated
 * due to the bands.  Therefore all of the differences between the *_b_me() functions
 * and their *_b() versions is to deal with the offset issue.
 * .........................................................................
 * 
 *******************************************************************************/

/* The alignment engines. 
 */
static float inside_b_jd(CM_t *cm, char *dsq, int L, 
			 int r, int z, int i0, int j0, 
			 int do_full,
			 float ***alpha, float ****ret_alpha, 
			 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			 void ****ret_shadow, 
			 int allow_begin, int *ret_b, float *ret_bsc,
			 int *jmin, int *jmax,
			 int **hdmin, int **hdmax,
			 int *safe_hdmin, int *safe_hdmax);

static float inside_b_jd_me(CM_t *cm, char *dsq, int L, 
			    int r, int z, int i0, int j0, 
			    int do_full,
			    float ***alpha, float ****ret_alpha, 
			    struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			    void ****ret_shadow, 
			    int allow_begin, int *ret_b, float *ret_bsc,
			    int *jmin, int *jmax,
			    int **hdmin, int **hdmax,
			    int *safe_hdmin, int *safe_hdmax);

/* The traceback routines.
 */

static float insideT_b_jd(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
			  int r, int z, int i0, int j0, int allow_begin,
			  int *jmin, int *jax, 
			  int **hdmin, int **hdmax,
			  int *safe_hdmin, int *safe_hdmax);

static float insideT_b_jd_me(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
			     int r, int z, int i0, int j0, int allow_begin,
			     int *jmin, int *jax,
			     int **hdmin, int **hdmax,
			     int *safe_hdmin, int *safe_hdmax);

/* Allocation functions. Derived from analogous functions in smallcyk.c */
static int   ** alloc_jdbanded_vjd_kshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
static char  ** alloc_jdbanded_vjd_yshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);

/* Debugging functions that print info to STDOUT */
static void debug_print_shadow_banded_jd(void ***shadow, CM_t *cm, int L, int *jmin, int *jmax, 
					 int **hdmin, int **hdmax);
static void debug_print_shadow_banded_deck_jd(int v, void ***shadow, CM_t *cm, int L, int *jmin, int *jmax,
					      int **hdmin, int **hdmax);


#define BE_EFFICIENT  0		/* setting for do_full: small memory mode */
#define BE_PARANOID   1		/* setting for do_full: keep whole matrix, perhaps for debugging */

/* Special flags for use in shadow (traceback) matrices, instead of
 * offsets to connected states. When yshad[0][][] is USED_LOCAL_BEGIN,
 * the b value returned by inside() is the best connected state (a 0->b
 * local entry). When yshad[v][][] is USED_EL, there is a v->EL transition
 * and the remaining subsequence is aligned to the EL state. 
 */
#define USED_LOCAL_BEGIN 101
#define USED_EL          102

/* Function: CYKInside_b_jd()
 *           EPN 11.04.05
 * based on CYKInside_b()
 * which was:
 * based on CYKInside()
 *
 * Only difference is bands are used in d and j dimesions: 
 *
 * Date:     SRE, Sun Jun  3 19:48:33 2001 [St. Louis]
 *
 * Purpose:  Wrapper for the insideT_b_jd() routine - solve
 *           a full alignment problem, return the traceback
 *           and the score, without dividing & conquering, using bands.
 *           
 *           Analogous to CYKDivideAndConquer() in many respects;
 *           see the more extensive comments in that function for
 *           more details on shared aspects.
 *           
 * Args:     cm     - the covariance model
 *           dsq    - the sequence, 1..L
 *           L      - length of sequence
 *           r      - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *           dmin   - minimum d bound for each state v; [0..v..M-1]
 *           dmax   - maximum d bound for each state v; [0..v..M-1]
 *
 * Returns:  score of the alignment in bits.
 */
float
CYKInside_b_jd(CM_t *cm, char *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, 
	       int *jmin, int *jmax, int **hdmin, int **hdmax, int *dmin, int *dmax)
{
  Parsetree_t *tr;
  int          z;
  float        sc;
  Stopwatch_t     *watch;
  
  watch = StopwatchCreate();
  StopwatchZero(watch);
  StopwatchStart(watch);

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));*/
  printf("alignment strategy:CYKInside_b_jd:b:nosmall\n"); 
  /*  printf("L: %d\n", L);*/
  
  /* Trust, but verify.
   * Check out input parameters.
   */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CM_LOCAL_BEGIN)) Die("internal error: we're not in local mode, but r is not root");
    if (cm->stid[r] != MATP_MP && cm->stid[r] != MATL_ML &&
	cm->stid[r] != MATR_MR && cm->stid[r] != BIF_B)
      Die("internal error: trying to do a local begin at a non-mainline start");
  }

  /* Create the parse tree, and initialize.
   */
  tr = CreateParsetree();
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */
  z  = cm->M-1;
  sc = 0.;

  /* Deal with case where we already know a local entry transition 0->r
   */
  if (r != 0)
    {
      InsertTraceNode(tr, 0,  TRACE_LEFT_CHILD, i0, j0, r);
      z  =  CMSubtreeFindEnd(cm, r);
      sc =  cm->beginsc[r];
    }

  /* Solve the whole thing with one call to insideT_b_jd.
     This calls a memory efficient (in the d-dimension only, not in the j dimension)
     insideT function, which only allocates cells in alpha within the bands.
     To use the non-memory efficient implementation call insideT_b()
     with the same arguments
   */
  sc += insideT_b_jd_me(cm, dsq, L, tr, r, z, i0, j0, (r==0), jmin, jmax, hdmin, hdmax, 
			dmin, dmax);

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  printf("returning from CYKInside_b_jd() sc : %f\n", sc); 

  StopwatchStop(watch);
  StopwatchDisplay(stdout, "banded CYK (j and d) CPU time: ", watch);
  StopwatchFree(watch);

  return sc;
}

 
/* EPN 
 * Function: inside_b_jd()
 * based on inside_b_me() which was ...
 * based on inside()
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phased of a CYK alignment algorithm
 *           using bands in the j and d dimension from an 
 *           HMM forwards backwards run. This function
 *           is memory efficient in the d dimension (only
 *           allocate within "safe" d bands). 
 *           Further assumes we're aligning a full sequence 
 *           (1..L) NOT a subsequence (i0..j0), and aligns it
 *           to the full model. This is very different from inside()
 *           which aligns a subsequence to a subtree of the model.
 *           
 *           Notes from inside():
 *           A note on the loop conventions. We're going to keep the
 *           sequence (dsq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: only vroot..vend
 *           in the decks, i0-1..j in the rows, and up to j0-i0+1 in
 *           the columns (d dimension). Where this is handled the most
 *           is in two variables: W, which is the length of the subsequence
 *           (j0-i0+1), and is oft used in place of L in the usual CYK;
 *           and jp (read: j'), which is the *relative* j w.r.t. the
 *           subsequence, ranging from 0..W, and then d ranges from 
 *           0 to jp, and j is calculated from jp (i0-1+jp).
 *           
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in an
 *           alpha, we expect that alpha[vroot..vend] are all NULL
 *           decks already; any other decks <vroot and >vend will
 *           be preserved. If we pass in a dpool, the decks *must* be
 *           sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *           We also deal with local begins, by keeping track of the optimal
 *           state that we could enter and account for the whole target 
 *           sequence: b = argmax_v  alpha_v(i0,j0) + log t_0(v),
 *           and bsc is the score for that. 
 *
 *           If vroot==0, i0==1, and j0==L (e.g. a complete alignment),
 *           the optimal alignment might use a local begin transition, 0->b,
 *           and we'd have to be able to trace that back. For any
 *           problem where the caller sets allow_begin, we return a valid b 
 *           (the optimal 0->b choice) and bsc (the score if 0->b is used).
 *           If a local begin is part of the optimal parse tree, the optimal
 *           alignment score returned by inside() will be bsc and yshad[0][L][L] 
 *           will be USE_LOCAL_BEGIN, telling insideT() to check b and
 *           start with a local 0->b entry transition. When inside()
 *           is called on smaller subproblems (v != 0 || i0 > 1 || j0
 *           < L), we're using inside() as an engine in divide &
 *           conquer, and we don't use the overall return score nor
 *           shadow matrices, but we do need allow_begin, b, and bsc for
 *           divide&conquer to sort out where a local begin might be used.
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           vroot     - first start state of subtree (0, for whole model)
 *           vend      - last end state of subtree (cm->M-1, for whole model)
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
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * int *safe_hdmin     - safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0)
 * int *safe_hdmax     - safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0)
 *                       
 * Returns: Score of the optimal alignment.  
 */
static float 
inside_b_jd(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	    float ***alpha, float ****ret_alpha, 
	    struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	    void ****ret_shadow, 
	    int allow_begin, int *ret_b, float *ret_bsc,
	    int *jmin, int *jmax, int **hdmin, int **hdmax,
	    int *safe_hdmin, int *safe_hdmax)
{
  float  **end;         /* we re-use the end deck. */
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */

  void  ***shadow;      /* shadow matrix for tracebacks */
  int    **kshad;       /* a shadow deck for bifurcations */
  char   **yshad;       /* a shadow deck for every other kind of state */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */

  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      dp_z;           /* d index for state z in alpha w/mem eff bands */
  int      kp;             /* k prime - keeps track of what k should be now
			     that we're using memory efficient bands */
  int      Wp;             /* W also changes depending on state */

  /* 11.04.05 jd addition: */
  if(i0 != 1)
    {
      printf("inside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");
      exit(1);
    }
  if(j0 != L)
    {
      printf("inside_b_jd requires that j0 be L. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
  if(vroot != 0)
    {
      printf("inside_b_jd requires that vroot be 0. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
  if(vend != cm->M-1)
    {
      printf("inside_b_jd requires that vend be cm->M-1. This function is not set up for subsequence alignment.\n");
      exit(1);
    }

  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops 
				 * This must be L because i0 must be 1 and j0 must be L
				 */
  
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, vroot, E_st);
  for (j = 0; j <= W; j++) {
    end[j][0] = 0.;
    for (d = 1; d <= j; d++) end[j][d] = IMPOSSIBLE;
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
  for (v = 0;     v < vroot; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) touch[v] = cm->pnum[v];
  for (v = vend+1;v < cm->M; v++) touch[v] = 0;

  /* The shadow matrix, if caller wants a traceback.
   * We do some pointer tricks here to save memory. The shadow matrix
   * is a void ***. Decks may either be char ** (usually) or
   * int ** (for bifurcation decks). Watch out for the casts.
   * For most states we only need
   * to keep y as traceback info, and y <= 6. For bifurcations,
   * we need to keep k, and k <= L, and L might be fairly big.
   * (We could probably limit k to an unsigned short ... anyone
   * aligning an RNA > 65536 would need a big computer... but
   * we'll hold off on that for now. We could also pack more
   * traceback pointers into a smaller space since we only really
   * need 3 bits, not 8.)
   */
  if (ret_shadow != NULL) {
    shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* Main recursion
   */
  for (v = vend; v >= vroot; v--) 
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
	/* CYK Full ME Bands used 1 */
	/* original line : alpha[v] = alloc_vjd_deck(L, i0, j0);*/
	alpha[v] = alloc_banded_vjd_deck(L, i0, j0, safe_hdmin[v], safe_hdmax[v]);
      
      if (ret_shadow != NULL) {
	if (cm->sttype[v] == B_st) {
	  /* CYK Full ME Bands used 2 */
	  /* original line : kshad     = alloc_vjd_kshadow_deck(L, i0, j0); */
	  kshad     = alloc_banded_vjd_kshadow_deck(L, i0, j0, safe_hdmin[v], safe_hdmax[v]);
	  shadow[v] = (void **) kshad;
	} else {
	  /* CYK Full ME Bands used 3 */
	  /* original line : yshad     = alloc_vjd_yshadow_deck(L, i0, j0); */
	  yshad     = alloc_banded_vjd_yshadow_deck(L, i0, j0, safe_hdmin[v], safe_hdmax[v]);
	  shadow[v] = (void **) yshad;
	}
      }

      /* 11.05.05
       * One strategy is to set all cells OUTSIDE bands to IMPOSSIBLE.
       * I think I'll run into problems doing this because some cells
       * are inside the j bands and inside the safe_hd bands, but not
       * inside the j dependent d bands. These cells though allocated, 
       * will potentially never get filled. 
       * One way to deal with this (though inefficient) is to set
       * ALL cells to impossible. Below is the the first strategy, only 
       * setting some cells to impossible. 
       
       **************************************************************
       * Strategy 1: only set cells outside j bands to IMPOSSIBLE:
       **************************************************************
       */
      /* j bands used 1.
       * Set all cells for j's outside of bands to IMPOSSIBLE 
       * Further, set any cells outside of hd j specific band for valid j's 
       * to IMPOSSIBLE. Remember, we've allocated only 
       * (safe_hdmax[v] - safe_hdmin[v] +1) cells for each vj deck.
       * Take advantage of fact that we know we're aligning the full sequence 1..L.
       */
      /* Following loop starts at safe_hdmin[v] because the j section
       * of a banded vjd deck is not allocated if j < dmin[v] because
       * there's no way that j can be used.
       */
      /*
      for (j = safe_hdmin[v]; j < jmin[v]; j++)
      {
      */
	  /* this j is outside the j band, set all d to IMPOSSIBLE */
      /*
	  for (d = safe_hdmin[v]; d <= safe_hdmax[v] && d <= j; d++)
	    {
	      alpha[v][j][d-safe_hdmin[v]] = IMPOSSIBLE;
	    }
	}

      for (j = jmax[v] + 1; j <= W; j++)
	{
      */
	  /* this j is outside the j band, set all d to IMPOSSIBLE */
      /*	  for (d = safe_hdmin[v]; d <= safe_hdmax[v] && d <= j; d++)
	    {
	      alpha[v][j][d-safe_hdmin[v]] = IMPOSSIBLE;
	    }
	}
	*************************************************************
	*/

      /**************************************************************
	* Strategy 2: set all allocated cells to IMPOSSIBLE at first.
	**************************************************************
	*/

      for (j = safe_hdmin[v]; j <= W; j++)
	{
	  for (d = safe_hdmin[v]; d <= safe_hdmax[v] && d <= j; d++)
	    {
	      alpha[v][j][d-safe_hdmin[v]] = IMPOSSIBLE;
	    }
	}

      //printf("2 v: %d\n", v);
	/*************************************************************/
      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  /* j bands used 2. */

	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      //printf("3 j: %d\n", j);
	      for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
		{
		  assert(d >= safe_hdmin[v]);
		  assert(d <= safe_hdmax[v]);
		  
		  y = cm->cfirst[v];
		  /* CYK Full ME Bands used 4 begin block */
		  /* original block */
		  /* alpha[v][j][d]  = cm->endsc[v];*/	/* init w/ local end */ 
		  /*if (ret_shadow != NULL) yshad[j][d]  = USED_EL; 
		    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		    if ((sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc; 
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		    }
		    if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		  */
		  /* new ME block */
		  dp_v = d - safe_hdmin[v];  /* d index for state v in alpha w/mem eff bands */
		  alpha[v][j][dp_v]  = cm->endsc[v];	/* init w/ local end */
		  if (ret_shadow != NULL) yshad[j][dp_v]  = USED_EL; 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		    {
		      dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
							    in alpha w/mem eff bands */
		      /* check to make sure the cell we're about to query is within the
			 bands for state y; this might be more complex than necessary */
		      if((dp_y >= 0) && ((dp_y < (j - (safe_hdmin[y+yoffset]) + 1))
					 && (dp_y < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
			{
			  if ((sc = alpha[y+yoffset][j][dp_y] + cm->tsc[v][yoffset]) >  alpha[v][j][dp_v]) {
			    alpha[v][j][dp_v] = sc; 
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  }
			}
		    }
		  if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		  /* CYK Full ME Bands used 4 end block */
		}
	    }
	}
      else if (cm->sttype[v] == B_st)
	{
	  /* j bands used 3. */
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //printf("3 j: %d\n", j);
	    /* Bands used */
	    /* old line :	for (d = 0; d <= jp; d++) */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);
		y = cm->cfirst[v];
		z = cm->cnum[v];

		/* CYK Full ME Bands used 5 begin block */
		/* original block */
		/*
		alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0];
		if (ret_shadow != NULL) kshad[j][d] = 0;
		for (k = 1; k <= d; k++)
		  if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) kshad[j][d] = k;
		  }
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/

		/* 11.04.05 Left this comment block here (from inside_b_me()) */

		/* new ME block : */
		/* 05.30.05 Fixed a small bug here */
		/* The changes made to this section of code in the memory efficient
		 * banded implementation are the most complex changes necessary to 
		 * get memory efficiency.  The reason is because there are indices in 
		 * two other states for a B_st, y and z (instead of just y).  This
		 * means that when we're dealing with a dp_v that is d minus a v-state
		 * specific offset, we also have to worry about the y-state offset
		 * and z-state offset.
		 * Let's set kp as the equivalent of k from the old code, but
		 * now we have to take into account the offsets.  To remain as
		 * consistent as possible with the old code, we will keep the
		 * indexing in z the same in the recursion, and figure out what
		 * the corresponding indices involving state y are.  
		 * So the old recursion code is : 
		 *
		 * for (jp = 0; jp <= W; jp++) {
		 * j = i0-1+jp;
		 * for (d = 0; d <= jp; d++) 
		 * {
		 *   alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0]; *INIT*
		 *   if (ret_shadow != NULL) kshad[j][d] = 0;
		 *   for (k = 1; k <= d; k++)
		 *   *RECURSION*
		 *   if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) {
		 *     alpha[v][j][d] = sc;
		 *     if (ret_shadow != NULL) kshad[j][d] = k; }
		 * 
		 * So we'll minimally change alpha[z][j][k] to alpha[z][j][kp]
		 * The INIT may change because although alpha[z][j][0] MUST be
		 * within the bands (because dmin[z] >= 0), the corresponding
		 * cell in alpha[y] might not be within the bands for y.  
		 * That cell is alpha[y][j-dmin[z]-kp][d-dmin[y]-dmin[z]-kp]
		 * because k = kp + dmin[z] (it probably takes some time writing
		 * down the new and old equations, and staring and thinking for a 
		 * while - I would write down more here - but this is already pretty
		 * verbose ... ).
		 * 
		 * Therefore we can't just start with k (or kp)  = 0 
		 * (like the old code did), because that might not be valid.
		 *
		 * First we need to determine the smallest kp for which we can 
		 * do a valid traceback, which means the alpha cell for both the y
		 * state and z state are within the bands.  For a kp to be valid given
		 * the following code, the following three inequalities have to be
		 * true.
		 *
		 * (1) d-dmin[z]-kp <= dmax[y]  
		 * (2) d-dmin[z]-kp >= dmin[y]
		 * (3) kp <= dmax[z]-dmin[z]
		 *
		 * (1) and (2) need to be satisified to guarantee that the cell we
		 * are going to access in the alpha[y] deck is within the bands for
		 * state y.  (3) is necessary to guarantee that the cell we are
		 * going to access in the alpha[z] deck is within the bands for 
		 * state z.
		 * We can rearrange 1 and 2 : 
		 *
		 * (1) kp >= d-dmax[y]-dmin[z]
		 * (2) kp <= d-dmin[y]-dmin[z]
		 * 
		 * First to check to see if ANY kp is valid, we can first
		 * check to make sure that (d-dmin[y]-dmin[z]) (RHS of (2))
		 * is >= 0.  If not, then kp can never be 0 or greater. 
		 * So it can never be valid. So we check for this at
		 * the beginning.
		 * 
		 * So, to find the minimal kp that satisfies (1), (2) and (3)
		 * I set kp = d-dmax[y]-dmin[z], and then check that it kp >= 0
		 * If kp < 0, we set it to 0.  Then we check to make sure kp
		 * satisfies (3) (It has to satisfy (2) if it satisfies (1)
		 * because dmax[y] >= dmin[y]).  This is our *INIT* assignment.
		 * Next we incrementally step through all valid kp values, we'll need 
		 * a for loop with two conditions to check in the 'while' portion.  
		 * Namely, that kp satisfies inequalities (2) and (3), that is
		 * kp <= (d-dmin[y]-dmin[z]) and kp <= (dmax[z]-dmin[z])
		 * This is marked in the code by *RECUR*
		 *
		 * Also, we want to make sure the while statement from the 
		 * original for loop (non-banded) is also satisfied.  This
		 * statement is k <= d.  We're dealing with kp, and k = kp+dmin[z]
		 * so this statement becomes kp <= d-dmin[z].  However, inequality
		 * (2) (kp <= d-dmin[y]-dmin[z]) takes care of this because dmin[y] >= 0
		 * 
		 */

		dp_v = d - safe_hdmin[v];  /* d index for state v in alpha w/mem eff bands */
		dp_y = d - safe_hdmin[y];  /* d index for state y in alpha w/mem eff bands */
		dp_z = d - safe_hdmin[z];  /* d index for state z in alpha w/mem eff bands */

		/* First make sure we have any valid kp, we know from inequality (2)
		   that kp <= d-safe_hdmin[y]-safe_hdmin[z] so if this is < 0 then no kp
		   is valid (see notes above) */

		if((d-safe_hdmin[y]-safe_hdmin[z]) >= 0)
		{
		  if(j < safe_hdmax[y]) kp = d-safe_hdmin[z]-j;
		  else kp = d-safe_hdmin[z]-safe_hdmax[y];
		  if(kp < 0) kp = 0;
		  if(kp <= safe_hdmax[z] - safe_hdmin[z]) /* make sure its valid in deck alpha[z] */
		    {
		      alpha[v][j][dp_v] = alpha[y][j-safe_hdmin[z]-kp][d-safe_hdmin[y]-safe_hdmin[z]-kp] 
			+ alpha[z][j][kp];
		      if (ret_shadow != NULL) kshad[j][dp_v] = kp;
		      for (kp = kp+1; kp <= (d-safe_hdmin[y]-safe_hdmin[z]) && kp <= (safe_hdmax[z]-safe_hdmin[z]);
			   kp++)
			{
			  /*printf("v is %d | checking y : %d z : %d\n", v, y, z);
			  printf("y comp          : alpha[%d][%d][%d] is %f\n", y, (j-safe_hdmin[z]-kp),(d-safe_hdmin[y]-safe_hdmin[z]-kp), 
				 alpha[y][j-safe_hdmin[z]-kp][d-safe_hdmin[y]-safe_hdmin[z]-kp]);
			  printf("z comp          : alpha[%d][%d][%d] is %f\n", z, j, kp, alpha[z][j][kp]);
			  printf("existing v comp : alpha[%d][%d][%d] is %f\n", v, j, dp_v, alpha[v][j][dp_v]);
			  printf("\n");*/
			  /* the following if statement ensures that the alpha cell for 
			     state y and the cell for state z that we are about to query 
			     is in fact within the bands for state y and state z respectively*/
			  if ((sc = alpha[y][j-safe_hdmin[z]-kp][d-safe_hdmin[y]-safe_hdmin[z]-kp] 
			       + alpha[z][j][kp]) > alpha[v][j][dp_v]) 
			    {
			      alpha[v][j][dp_v] = sc;
			      if (ret_shadow != NULL) kshad[j][dp_v] = kp;
			    }
			}
		    }
		}
		else alpha[v][j][dp_v] = IMPOSSIBLE;
		/*else Die("cell in alpha matrix was not filled in due to bands.\n");*/
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 5 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //printf("3 j: %d\n", j);

	    /* CYK Full ME Bands used 6 */
	    /* Deleted because I realized this was no longer needed */

	    /* Bands used */
	    /* old line :	for (d = 2; d <= jp; d++) */
	    /* we assume hdmin[v][j-jmin[v]] >= 2 */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);

		y = cm->cfirst[v];
		/* CYK Full ME Bands used 7 block */
		/* original code block below : */
		/*
		alpha[v][j][d] = cm->endsc[v];  init w/ local end 
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-2] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		
		i = j-d+1;
		if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/
		/* new ME code block : */
		dp_v = d - safe_hdmin[v]; /* d index for state v in alpha w/mem eff bands */
		/*printf("j: %d | v: %d | d: %d | dp_v: %d | safe_hdmin[v]: %d\n", j, v, d, dp_v, safe_hdmin[v]);*/
				
		alpha[v][j][dp_v] = cm->endsc[v]; /* init w/ local end */
		if(ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
							  in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-2) >= 0) && (((dp_y-2) < (j - (safe_hdmin[y+yoffset]) + 1))
					   && ((dp_y-2) < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j-1][dp_y-2] + cm->tsc[v][yoffset]) >  alpha[v][j][dp_v])
			  {
			    alpha[v][j][dp_v] = sc;
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  }
		      }
		  }
		i = j-d+1;
		if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		  alpha[v][j][dp_v] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		else
		  alpha[v][j][dp_v] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);
		
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 7 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //	    printf("3 j: %d\n", j);
	    /* CYK Full ME Bands used 8 */
	    /* Deleted because I realized this was no longer needed */

	    /* Bands used */
	    /* old line :	for (d = 1; d <= jp; d++) */
	    /* we assume safe_hdmin[v] >= 1 */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);

		y = cm->cfirst[v];
		/* CYK Full ME Bands used 9 block */
		/* original code block below : */
		/*
		alpha[v][j][d] = cm->endsc[v];
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  } 
		
		i = j-d+1;
		if (dsq[i] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/
		/* new ME code block : */
		dp_v = d - safe_hdmin[v]; /* d index for state v in alpha w/mem eff bands */
		/*printf("v: %d | j: %d | dp_v: %d | j-jmin[v]: %d | safe_hdmin[%d]: %d | d: %d\n", v, j, dp_v, (j-jmin[v]), v, safe_hdmin[v], d);*/
		alpha[v][j][dp_v] = cm->endsc[v]; /* init w/ local end */
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-1) >= 0) && (((dp_y-1) < (j - (safe_hdmin[y+yoffset]) + 1))
					   && ((dp_y-1) < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j][dp_y-1] + cm->tsc[v][yoffset]) >  alpha[v][j][dp_v]) 
			  {
			    alpha[v][j][dp_v] = sc;
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  } 
		      }
		  }
		i = j-d+1;
		if (dsq[i] < Alphabet_size)
		  alpha[v][j][dp_v] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][dp_v] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 9 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++) {
	    //printf("3 j: %d\n", j);
	    /* CYK Full ME Bands used 10 */
	    /* Deleted because I realized this was no longer needed */

	    /* Bands used */
	    /* old line :	for (d = 1; d <= jp; d++) */
	    /* we assume safe_hdmin[v] >= 1 */
	    for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++)
	      {
		assert(d >= safe_hdmin[v]);
		assert(d <= safe_hdmax[v]);

		y = cm->cfirst[v];
		/* CYK Full ME Bands used 11 block */
		/* original code block below : */
		/*
		alpha[v][j][d] = cm->endsc[v];
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		if (dsq[j] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
		*/
		/* new ME code block : */
		dp_v = d - safe_hdmin[v]; /* d index for state v in alpha w/mem eff bands */
		alpha[v][j][dp_v] = cm->endsc[v]; /* init w/ local end */
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - safe_hdmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-1) >= 0) && (((dp_y-1) < (j - (safe_hdmin[y+yoffset]) + 1))
				      && ((dp_y-1) < (safe_hdmax[y+yoffset] - safe_hdmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j-1][dp_y-1] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v])
			  {
			    alpha[v][j][dp_v] = sc;
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  }
		      }
		  }
		if (dsq[j] < Alphabet_size)
		  alpha[v][j][dp_v] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][dp_v] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 11 end block */
	      }
	  }
	}				/* finished calculating deck v. */
      
      /* CYK Full ME Bands used 12 block */
      /* The following loops originally access alpha[v][j0][W] but the index W will be
	 in different positions due to the bands */

      /* ME Added the following two lines */
      Wp = W - safe_hdmin[v];
      /* We need to make sure that Wp is within the bands */
      if(Wp >= 0 && Wp <= (safe_hdmax[v] - safe_hdmin[v]))
	{
	  /* ME all subsequent changes in this block simply replace
	     W with Wp (so wherever Wp is, there used to be a W) */

	  /* Check for local begin getting us to the root.
	   * This is "off-shadow": if/when we trace back, we'll handle this
	   * case separately (and we'll know to do it because we'll immediately
	   * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
	   * to jump right to state b; see below)
	   */
	  if (allow_begin && alpha[v][j0][Wp] + cm->beginsc[v] > bsc) 
	    {
	      b   = v;
	      bsc = alpha[v][j0][Wp] + cm->beginsc[v];
	    }

	  /* Check for whether we need to store an optimal local begin score
	   * as the optimal overall score, and if we need to put a flag
	   * in the shadow matrix telling insideT() to use the b we return.
	   */
	  if (allow_begin && v == 0 && bsc > alpha[0][j0][Wp]) {
	    alpha[0][j0][Wp] = bsc;
	    if (ret_shadow != NULL) yshad[j0][Wp] = USED_LOCAL_BEGIN;
	  }
	}
      /* CYK Full ME Bands used 12 end block */

      /* CYK Full ME Bands used 13 block */
      /* The following block implements the deck reuse strategy, however, here
	 we can't do that, because for each state, the bands are different, so 
	 we can't use old Decks, but rather must allocate a new one, and free
	 the old one. Lines specific to ME are indicated, and original lines
	 are commented out */

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
	    y = cm->cfirst[v];
	    z = cm->cnum[v];  
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
	/* CYK Full ME Bands used 13 end block */
      }
    } /* end loop over all v */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  
  /* CYK Full ME Bands used 14 */
  /* original line :  sc       = alpha[vroot][j0][W];*/
  Wp = W - safe_hdmin[vroot];
  sc       = alpha[vroot][j0][Wp];

  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = vroot; v <= vend; v++) /* be careful of our reuse of the end deck -- free it only once */
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
  if (ret_shadow != NULL) *ret_shadow = shadow;
  return sc;
}

/* Function: insideT_b_jd()
 *           EPN 05.24.05
 * *based on insideT(), only difference is memory efficient bands are used : 
 *
 * Date:     SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]
 *
 * Purpose:  Call inside, get vjd shadow matrix;
 *           then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 */
static float
insideT_b_jd(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
	     int r, int z, int i0, int j0, 
	     int allow_begin, int *jmin, int *jmax,
	     int **hdmin, int **hdmax,
	     int *safe_hdmin, int *safe_hdmax)
{
  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  Nstack_t *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;
  int       dp;                 /* add explanation */
  int       kp;                 /* add explanation */

  sc = inside_b_jd(cm, dsq, L, r, z, i0, j0, 
		   BE_EFFICIENT,	/* memory-saving mode */
		   NULL, NULL,	/* manage your own matrix, I don't want it */
		   NULL, NULL,	/* manage your own deckpool, I don't want it */
		   &shadow,		/* return a shadow matrix to me. */
		   allow_begin,      /* TRUE to allow local begins */
		   &b, &bsc,	/* if allow_begin is TRUE, gives info on optimal b */
		   jmin, jmax,    /* bands on j */
		   hdmin, hdmax,  /* j dependent bands on d */
		   safe_hdmin, safe_hdmax);

  pda = CreateNstack();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  /*printf("Starting traceback in insideT_b_me()\n");*/
  while (1) {
    /* CYK Full ME Bands used 15 */
    /* 2 lines below added */
    /*
    assert(j >= jmin[v]);
    assert(j <= jmax[v]);
    assert(d >= hdmin[v][j0]);
    assert(d <= hdmax[v][j0]);

    assert(d >= safe_hdmin[v]);
    assert(d <= safe_hdmax[v]);
    */
    if(cm->sttype[v] != EL_st && d > safe_hdmax[v])
      {
	printf("ERROR in insideT_b_jd(). d : %d > safe_hdmax[%d] (%d)\n", d, v, safe_hdmax[v]);
	exit(1);
      }
    if(cm->sttype[v] != EL_st && d < safe_hdmin[v])
      {
	printf("ERROR in insideT_b_jd(). d : %d < safe_hdmin[%d] (%d)\n", d, v, safe_hdmin[v]);
	exit(1);
      }

    /* Deal with end local states */
    if(cm->sttype[v] != EL_st)
      dp = d - safe_hdmin[v];
    else
      dp = d;
    
    if (cm->sttype[v] == B_st) {
      /* CYK Full ME Bands used 16 */
      /* original line : k = ((int **) shadow[v])[j][d];  */
      /* new 3 lines below replace it */
      assert(v >= 0);
      kp = ((int **) shadow[v])[j][dp];   /* kp = offset len of right fragment */
      z = cm->cnum[v];
      k = kp + safe_hdmin[z];  /* k = offset len of right fragment */
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      PushNstack(pda, j);	/* remember the end j    */
      PushNstack(pda, k);	/* remember the subseq length k */
      PushNstack(pda, tr->n-1);	/* remember the trace index of the parent B state */
      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (! PopNstack(pda, &bifparent)) break;
      /* CYK Full ME Bands used 17 */
      /* original line : PopNstack(pda, &d); */
      PopNstack(pda, &dp);
      /* CYK Full ME Bands used 18 */
      /* line below added */

      /* Deal with end local states */
      if(cm->sttype[v] != EL_st)
	d = dp + safe_hdmin[y];
      else
      d = dp;

      PopNstack(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
    } else {
      yoffset = ((char **) shadow[v])[j][dp];

      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    Die("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	}
    }
  }
  FreeNstack(pda);  /* it should be empty; we could check; naaah. */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  return sc;
}

/****************************************************************************/
void
PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax,
		     int W)
{
  int v;
  int j;
  int max;
  double after, before;

  printf("Printing DP cells saved using j and d bands:\n");
  before = after = 0;
  for (v = 0; v < cm->M; v++) 
    {
      for(j = 0; j <= W; j++)
	{
	  if (cm->sttype[v] != E_st) 
	    {
	      before += j + 1;
	    }
	}
      for(j = jmin[v]; j <= jmax[v]; j++)
	{
	  if (cm->sttype[v] != E_st) 
	    {
	      max = (j < hdmax[v][j-jmin[v]]) ? j : hdmax[v][j-jmin[v]];
	      after += max - hdmin[v][j-jmin[v]] + 1;
	    }
	}
    }
  printf("Before:  something like %.0f\n", before);
  printf("After:   something like %.0f\n", after);
  printf("Speedup: maybe %.2f fold\n\n", (float) before / (float) after);
}

 
/* EPN 03.29.06
 * Function: inside_b_jd_me()
 * based on inside_b_me() which was ...
 * based on inside()
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phase of a CYK alignment algorithm
 *           using bands in the j and d dimension from obtained
 *           from an HMM forwards-backwards run. This function
 *           is memory efficient in the j AND d dimension.
 * 
 *           To be able to consistently handle end states, the
 *           original SRE behavior of reusing the end deck was
 *           abandoned. Now each end state has its own deck, which
 *           makes this implementation easier because each state
 *           has its own bands on j, and thus has a state specific
 *           offset with alpha[end][jp][dp] in the banded mem eff 
 *           matrix corresponding to alpha[end][jp+jmin[end]][dp+hdmin[v][jp_v]]
 *           in the platonic matrix.
 *           
 *           Notes from inside():
 *           A note on the loop conventions. We're going to keep the
 *           sequence (dsq) and the matrix (alpha) in the full coordinate
 *           system: [0..v..M-1][0..j..L][0..d..j]. However, we're
 *           only calculating a part of that matrix: only vroot..vend
 *           in the decks, i0-1..j in the rows, and up to j0-i0+1 in
 *           the columns (d dimension). Where this is handled the most
 *           is in two variables: W, which is the length of the subsequence
 *           (j0-i0+1), and is oft used in place of L in the usual CYK;
 *           and jp (read: j'), which is the *relative* j w.r.t. the
 *           subsequence, ranging from 0..W, and then d ranges from 
 *           0 to jp, and j is calculated from jp (i0-1+jp).
 *           
 *           The caller is allowed to provide us with a preexisting
 *           matrix and/or deckpool (thru "alpha" and "dpool"), or
 *           have them newly created by passing NULL. If we pass in an
 *           alpha, we expect that alpha[vroot..vend] are all NULL
 *           decks already; any other decks <vroot and >vend will
 *           be preserved. If we pass in a dpool, the decks *must* be
 *           sized for the same subsequence i0,j0.
 *           
 *           Note that the (alpha, ret_alpha) calling idiom allows the
 *           caller to provide an existing matrix or not, and to
 *           retrieve the calculated matrix or not, in any combination.
 *           
 *           We also deal with local begins, by keeping track of the optimal
 *           state that we could enter and account for the whole target 
 *           sequence: b = argmax_v  alpha_v(i0,j0) + log t_0(v),
 *           and bsc is the score for that. 
 *
 *           If vroot==0, i0==1, and j0==L (e.g. a complete alignment),
 *           the optimal alignment might use a local begin transition, 0->b,
 *           and we'd have to be able to trace that back. For any
 *           problem where the caller sets allow_begin, we return a valid b 
 *           (the optimal 0->b choice) and bsc (the score if 0->b is used).
 *           If a local begin is part of the optimal parse tree, the optimal
 *           alignment score returned by inside() will be bsc and yshad[0][L][L] 
 *           will be USE_LOCAL_BEGIN, telling insideT() to check b and
 *           start with a local 0->b entry transition. When inside()
 *           is called on smaller subproblems (v != 0 || i0 > 1 || j0
 *           < L), we're using inside() as an engine in divide &
 *           conquer, and we don't use the overall return score nor
 *           shadow matrices, but we do need allow_begin, b, and bsc for
 *           divide&conquer to sort out where a local begin might be used.
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           vroot     - first start state of subtree (0, for whole model)
 *           vend      - last end state of subtree (cm->M-1, for whole model)
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
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           jmin      - minimum j bound for each state v; [0..v..M-1]
 *           jmax      - maximum j bound for each state v; [0..v..M-1]
 *           hdmin     - minimum d bound for each state v and valid j; 
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 *           hdmax     - maximum d bound for each state v and valid j;
 *                       [0..v..M-1][0..j0..(jmax[v]-jmin[v])]
 *                       careful: j dimension offset. j0-jmin[v] = j;
 * int *safe_hdmin     - safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0)
 * int *safe_hdmax     - safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0)
 *                       
 * Returns: Score of the optimal alignment.  
 */
static float 
inside_b_jd_me(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	       float ***alpha, float ****ret_alpha, 
	       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	       void ****ret_shadow, 
	       int allow_begin, int *ret_b, float *ret_bsc,
	       int *jmin, int *jmax, int **hdmin, int **hdmax,
	       int *safe_hdmin, int *safe_hdmax)
{
  int      nends;       /* counter that tracks when we can release end deck to the pool */
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  void  ***shadow;      /* shadow matrix for tracebacks */
  int    **kshad;       /* a shadow deck for bifurcations */
  char   **yshad;       /* a shadow deck for every other kind of state */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */

  /* variables used for memory efficient bands */
  int      dp_v;           /* d index for state v in alpha w/mem eff bands */
  int      dp_y;           /* d index for state y in alpha w/mem eff bands */
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;
  float  **tmp_deck;       /* temp variable, used only to free deckpool at end */

  /* 11.04.05 jd addition: */
  /*  if(i0 != 1)
    {
      printf("inside_b_jd requires that i0 be 1. This function is not set up for subsequence alignment\n");
      exit(1);
    }
  if(j0 != L)
    {
      printf("inside_b_jd requires that j0 be L. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
  if(vroot != 0)
    {
      printf("inside_b_jd requires that vroot be 0. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
  if(vend != cm->M-1)
    {
      printf("inside_b_jd requires that vend be cm->M-1. This function is not set up for subsequence alignment.\n");
      exit(1);
    }
  */
  
  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();
  nends = CMSubtreeCountStatetype(cm, vroot, E_st);

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
  for (v = 0;     v < vroot; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) touch[v] = cm->pnum[v];
  for (v = vend+1;v < cm->M; v++) touch[v] = 0;

  /* The shadow matrix, if caller wants a traceback.
   * We do some pointer tricks here to save memory. The shadow matrix
   * is a void ***. Decks may either be char ** (usually) or
   * int ** (for bifurcation decks). Watch out for the casts.
   * For most states we only need
   * to keep y as traceback info, and y <= 6. For bifurcations,
   * we need to keep k, and k <= L, and L might be fairly big.
   * (We could probably limit k to an unsigned short ... anyone
   * aligning an RNA > 65536 would need a big computer... but
   * we'll hold off on that for now. We could also pack more
   * traceback pointers into a smaller space since we only really
   * need 3 bits, not 8.)
   */
  if (ret_shadow != NULL) {
    shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* Main recursion
   */
  for (v = vend; v >= vroot; v--) 
    {
      /*printf("v: %d\n", v);*/
      /* First we need a deck to fill in.
       * 1. see if we can take something from the pool
       *    if we're an end state, fill it
       * 2. else, allocate a new deck.
       */
      if (! deckpool_pop(dpool, &(alpha[v]))) 
	alpha[v] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);

      if (cm->sttype[v] != E_st) {
	if (ret_shadow != NULL) {
	  if (cm->sttype[v] == B_st) {
	    /* CYK Full ME Bands used 2 */
	    /* original line : kshad     = alloc_vjd_kshadow_deck(L, i0, j0); */
	    kshad     = alloc_jdbanded_vjd_kshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	    shadow[v] = (void **) kshad;
	  } else {
	    /* CYK Full ME Bands used 3 */
	    /* original line : yshad     = alloc_vjd_yshadow_deck(L, i0, j0); */
	    yshad     = alloc_jdbanded_vjd_yshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	    shadow[v] = (void **) yshad;
	  }
	}
      }

      /* We've only allocated alpha cells that are within the bands
       * on the j and d dimensions. This means we have to deal
       * with all sorts of offset issues, but we don't have to 
       * waste time setting cells outside the bands to IMPOSSIBLE.
       */
      if (cm->sttype[v] == E_st)
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  if(d != 0)
		    Die("band on E state %d has a non-zero d value within its j band for j:%d\n", v, j);
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v
						 in alpha w/mem eff bands */
		  alpha[v][jp_v][dp_v] = 0.; /* for End states, d must be 0 */
		}		    
	    }
	  continue;
	}  
      else if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (j = jmin[v]; j <= jmax[v]; j++)
	    {
	      jp_v = j - jmin[v];
	      for (d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; d++)
		{
		  dp_v = d - hdmin[v][jp_v];  /* d index for state v in alpha w/mem eff bands */
		  alpha[v][jp_v][dp_v]  = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  /* treat EL as emitting only on self transition */
		  if (ret_shadow != NULL) yshad[jp_v][dp_v]  = USED_EL; 
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
			      if ((sc = alpha[y][jp_y][dp_y] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
				{
				  alpha[v][jp_v][dp_v] = sc; 
				  if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
				}
			    }
			}
		    }
		  if (alpha[v][jp_v][dp_v] < IMPOSSIBLE)
		    alpha[v][jp_v][dp_v] = IMPOSSIBLE;
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
			   * alpha[y][jp_y-k][d-hdmin[jp_y-k]-k],
			   * and alpha[z][jp_z][k-hdmin[jp_z]];
			   */
			  kp_z = k-hdmin[z][jp_z];
			  dp_y = d-hdmin[y][jp_y-k];

			  if ((sc = alpha[y][jp_y-k][dp_y - k] + alpha[z][jp_z][kp_z]) 
			      > alpha[v][jp_v][dp_v])
			    {
			      alpha[v][jp_v][dp_v] = sc;
			      if (ret_shadow != NULL) kshad[jp_v][dp_v] = kp_z;
			    }
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
		if(ret_shadow != NULL) yshad[jp_v][dp_v] = USED_EL;
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
			    if ((sc = alpha[y][jp_y-1][dp_y-2] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
			      {
				alpha[v][jp_v][dp_v] = sc; 
				if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
			      }
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
		if(ret_shadow != NULL) yshad[jp_v][dp_v] = USED_EL;
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
			    if ((sc = alpha[y][jp_y][dp_y-1] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
			      {
				alpha[v][jp_v][dp_v] = sc; 
				if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
			      }
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
		if(ret_shadow != NULL) yshad[jp_v][dp_v] = USED_EL;
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
			    if ((sc = alpha[y][jp_y-1][dp_y-1] + cm->tsc[v][yoffset]) > alpha[v][jp_v][dp_v])
			      {
				alpha[v][jp_v][dp_v] = sc; 
				if (ret_shadow != NULL) yshad[jp_v][dp_v] = yoffset;
			      }
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
	      if (allow_begin && alpha[v][jp_v][Wp] + cm->beginsc[v] > bsc) 
		{
		  b   = v;
		  bsc = alpha[v][jp_v][Wp] + cm->beginsc[v];
		}
	    }
	}
      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (v == 0)
	{
	  if(j0 >= jmin[0] && j0 <= jmax[0])
	    {
	      jp_v = j0 - jmin[v];
	      if(W >= hdmin[v][jp_v] && W <= hdmax[v][jp_v])
		{
		  if (allow_begin && v == 0 && bsc > alpha[0][jp_v][Wp]) {
		    alpha[0][jp_v][Wp] = bsc;
		    if (ret_shadow != NULL) yshad[jp_v][Wp] = USED_LOCAL_BEGIN;
		  }
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
	    y = cm->cfirst[v];
	    z = cm->cnum[v];  
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
			free_vjd_deck(alpha[y], i0, j0);
		      }
		      alpha[y] = NULL;
		  }
	      }
	  }
      }
    } /* end loop over all v */
  /*debug_print_alpha_banded_jd(alpha, cm, L, jmin, jmax, hdmin, hdmax);*/

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid 
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  
  Wp = W - hdmin[vroot][j0-jmin[vroot]];
  sc =     alpha[vroot][j0-jmin[vroot]][Wp];

  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = vroot; v <= vend; v++) 
      if (alpha[v] != NULL) { 
	deckpool_push(dpool, alpha[v]); alpha[v] = NULL;
      }
    free(alpha);
  } else *ret_alpha = alpha;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    while (deckpool_pop(dpool, &tmp_deck)) free_vjd_deck(tmp_deck, i0, j0);
    deckpool_free(dpool);
  } else {
    *ret_dpool = dpool;
  }

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  /*printf("inside jd me returning sc: %f\n", sc);*/

  return sc;
}

/* Function: insideT_b_jd_me()
 *           EPN 03.29.06
 * *based on insideT(), only difference is memory efficient bands on the j and d dimensions
 *  are used : 
 *
 * Date:     SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]
 *
 * Purpose:  Call inside, get vjd shadow matrix;
 *           then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 */
static float
insideT_b_jd_me(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
		int r, int z, int i0, int j0, 
		int allow_begin, int *jmin, int *jmax,
		int **hdmin, int **hdmax,
		int *safe_hdmin, int *safe_hdmax)
{
  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  Nstack_t *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;
  int       jp_v;               /* j-jmin[v] for current j, and current v */
  int       dp_v;               /* d-hdmin[v][jp_v] for current j, current v, current d*/
  int       jp_z;               /* j-jmin[z] for current j, and current z */
  int       kp_z;               /* the k value (d dim) from the shadow matrix
				 * giving the len of right fragment offset in deck z,
				 * k = kp_z + hdmin[z][jp_z]*/

  sc = inside_b_jd_me(cm, dsq, L, r, z, i0, j0, 
		      BE_EFFICIENT,	/* memory-saving mode */
		      /*BE_PARANOID,*/	/* non-memory-saving mode */
		      NULL, NULL,	/* manage your own matrix, I don't want it */
		      NULL, NULL,	/* manage your own deckpool, I don't want it */
		      &shadow,		/* return a shadow matrix to me. */
		      allow_begin,      /* TRUE to allow local begins */
		      &b, &bsc,	/* if allow_begin is TRUE, gives info on optimal b */
		      jmin, jmax,    /* bands on j */
		      hdmin, hdmax,  /* j dependent bands on d */
		      safe_hdmin, safe_hdmax);

  pda = CreateNstack();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];

  while (1) {
    if(cm->sttype[v] != EL_st && d > hdmax[v][jp_v])
      {
	printf("ERROR in insideT_b_jd(). d : %d > safe_hdmax[%d] (%d)\n", d, v, safe_hdmax[v]);
	exit(1);
      }
    if(cm->sttype[v] != EL_st && d < hdmin[v][jp_v])
      {
	printf("ERROR in insideT_b_jd(). d : %d < safe_hdmin[%d] (%d)\n", d, v, safe_hdmin[v]);
	exit(1);
      }
    
    if (cm->sttype[v] == B_st) {
      kp_z = ((int **) shadow[v])[jp_v][dp_v];   /* kp = offset len of right fragment */
      z = cm->cnum[v];
      jp_z = j-jmin[z];
      k = kp_z + hdmin[z][jp_z];  /* k = offset len of right fragment */
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      PushNstack(pda, j);	/* remember the end j    */
      PushNstack(pda, k);	/* remember the subseq length k */
      PushNstack(pda, tr->n-1);	/* remember the trace index of the parent B state */
      /* Deal with attaching left start state.
       */
      j = j-k;
      d = d-k;
      i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      /* We don't trace back from an E or EL. Instead, we're done with the
       * left branch of the tree, and we try to swing over to the right
       * branch by popping a right start off the stack and attaching
       * it. If the stack is empty, then we're done with the
       * traceback altogether. This is the only way to break the
       * while (1) loop.
       */
      if (! PopNstack(pda, &bifparent)) break;
      PopNstack(pda, &d);
      PopNstack(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
      jp_v = j - jmin[v];
      dp_v = d - hdmin[v][jp_v];
    } else {
      yoffset = ((char **) shadow[v])[jp_v][dp_v];
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    Die("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL) 
	{	/* a local alignment end */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;		/* now we're in EL. */
	  jp_v = j;
	  dp_v = d;
	}
      else if (yoffset == USED_LOCAL_BEGIN) 
	{ /* local begin; can only happen once, from root */
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	  v = b;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
      else 
	{
	  y = cm->cfirst[v] + yoffset;
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	  v = y;
	  jp_v = j - jmin[v];
	  dp_v = d - hdmin[v][jp_v];
	}
    }
  }
  FreeNstack(pda);  /* it should be empty; we could check; naaah. */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  return sc;
}

  
/* Functions: *jdbanded_*_vjd_*
 * EPN 03.29.06 these functions were derived from their 
 *              *_vjd_* analogs from SRE's smallcyk.c
 * notes from smallcyk.c:
 * Date:     SRE, Sat Aug 12 16:27:37 2000 [Titusville]
 *
 * Purpose:  Allocation and freeing of 3D matrices and 2D decks
 *           in the vjd coord system. These can be called on
 *           subsequences i..j, not just the full sequence 1..L,
 *           so they need i,j... if you're doing the full sequence
 *           just pass 1,L.
 *           
 *           Also deal with shadow matrices and shadow decks in the
 *           vjd coordinate system. Note that bifurcation shadow decks
 *           need more dynamic range than other shadow decks, hence
 *           a separation into "kshadow" (BIFURC) and "yshadow" (other
 *           states) decks, and some casting shenanigans in
 *           a full ***shadow matrix.
 *           
 *           Values in yshad are offsets to the next connected state,
 *           or a flag for local alignment. Possible offsets range from
 *           0..5 (maximum of 6 connected states). The flags are
 *           USED_LOCAL_BEGIN (101) and USED_EL (102), defined at
 *           the top of this file. Only yshad[0][L][L] (e.g. root state 0,
 *           aligned to the whole sequence) may be set to USED_LOCAL_BEGIN.
 *           (Remember that the dynamic range of yshad, as a char, is 
 *           0..127, in ANSI C; we don't know if a machine will make it
 *           signed or unsigned.)
 */
float **
alloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  float **a;
  int     jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;
  /*printf("in alloc JD banded vjd deck, L : %d, i : %d, j : %d, jmin : %d, jmax : %d\n", L, i, j, jmin, jmax);*/

  if(j < jmin || i > jmax)
    {
      printf("ERROR called alloc_jdbanded_vjd_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);
      exit(1);
    }

  SQD_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  a = MallocOrDie(sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0; jp <= L;     jp++) a[jp]     = NULL;

  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jfirst: %d | jlast: %d\n", jfirst, jlast);
      printf("jp: %d | max : %d\n", jp, (jlast)); 
      printf("hdmax[%d]: %d\n", (jp-jmin), hdmax[jp-jmin]);
      */
      if(hdmax[jp-jmin] > (jp+1))
	{
	  /* Based on my current understanding this shouldn't happen, it means there's a valid d
	   * in the hd band that is invalid because its > j. I check, or ensure, that this
	   * doesn't happen when I'm constructing the d bands.
	   */
	  printf("jd banded error 0.\n");
	  exit(1);
	}
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*a is offset only the first (jlast-jfirst+1) elements will be non-NULL*/
      a[jp-jfirst] = MallocOrDie(sizeof(float) * bw);
      /*printf("\tallocated a[%d] | bw: %d\n", (jp-jfirst), bw);*/
    }
  return a;
}

static char  **
alloc_jdbanded_vjd_yshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  char **a;
  int    jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;

  if(j < jmin || i > jmax)
    {
      printf("ERROR called alloc_jdbanded_vjd_yshadow_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);
      exit(1);
    }

  SQD_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  a = MallocOrDie(sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  for (jp = (jlast-jfirst+1); jp <= L;     jp++) a[jp]     = NULL;

  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jp: %d | max : %d\n", jp, (jlast)); */
      if(hdmax[jp-jmin] > (jp+1))
	{
	  /* Based on my current understanding this shouldn't happen, it means there's a valid d
	   * in the hd band that is invalid because its > j. I think I check, or ensure, that this
	   * doesn't happen when I'm constructing the d bands.
	   */
	  printf("jd banded error 1.\n");
	  exit(1);
	}
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*printf("\tallocated a[%d]\n", (jp-jfirst));*/
      a[jp-jfirst] = MallocOrDie(sizeof(char) * bw);
    }
  return a;
}
static int** 
alloc_jdbanded_vjd_kshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  int **a;
  int   jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;

  if(j < jmin || i > jmax)
    {
      printf("ERROR called alloc_jdbanded_vjd_kshadow_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);
      exit(1);
    }

  SQD_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  a = MallocOrDie(sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  for (jp = (jlast-jfirst+1); jp <= L;     jp++) a[jp]     = NULL;

  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jp: %d | max : %d\n", jp, (jlast)); */
      if(hdmax[jp-jmin] > (jp+1))
	{
	  /* Based on my current understanding this shouldn't happen, it means there's a valid d
	   * in the hd band that is invalid because its > j. I think I check, or ensure, that this
	   * doesn't happen when I'm constructing the d bands.
	   */
	  printf("jd banded error 0.\n");
	  exit(1);
	}
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*printf("\tallocated a[%d]\n", (jp-jfirst));*/
      a[jp-jfirst] = MallocOrDie(sizeof(int) * bw);
    }
  return a;

}

/******************************************************************/
/* The below functions were written during debugging, and print
   out either the shadow or alpha matrix.  They are kept
   here just in case they're needed again.  Note : the functions
   that print out the entire matrix are really only useful
   when the BE_PARANOID flag is set, meaning that decks are
   never freed until the end.
*/
/*================================================================*/

/* EPN 05.16.05
   debug_print_shadow_banded()
 * Function: debug_print_shadow_banded
 *
 * Purpose:  Print banded shadow matrix 
 */

static void
debug_print_shadow_banded_jd(void ***shadow, CM_t *cm, int L, int *jmin, int *jmax, 
			     int **hdmin, int **hdmax)
{
  int v, j, d, dp, jp;
  int yoffset;
  char yoffset_c;

  printf("\nPrinting banded shadow matrix :\n");
  printf("************************************\n");
  for(v = 0; v < cm->M; v++)
    {
      printf("====================================\n");
      for(j = jmin[v]; j <= jmax[v]; j++)
	{
	  printf("------------------------------------\n");
	  for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	    {
	      jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	      dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	      if(cm->sttype[v] == E_st)
		{
		  printf("END state\n");
		}
	      else
		{
		  if(cm->sttype[v] == B_st)
		    {
		      yoffset = ((int **) shadow[v])[jp][dp];
		      printf("INT  shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
		    }
		  else
		    {
		      yoffset_c = ((char **) shadow[v])[jp][dp];
		      printf("CHAR shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
		    }
		}
	    }
	}
    }
  printf("****************\n\n");
}

/* EPN 05.16.05
   debug_print_shadow_banded_deck_jd()
 * Function: debug_print_shadow_banded_deck_jd
 *
 * Purpose:  Print banded (in j and d dimensions) shadow matrix deck
 */

static void
debug_print_shadow_banded_deck_jd(int v, void ***shadow, CM_t *cm, int L, int *jmin, int *jmax,
				  int **hdmin, int **hdmax)
{
  int j, d, dp, jp;
  int yoffset;

  printf("\nPrinting banded shadow matrix deck for v : %d:\n", v);
  printf("====================================\n");
  for(j = jmin[v]; j <= jmax[v]; j++)
    {
      printf("------------------------------------\n");
      for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	{
	  jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	  dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	  if(cm->sttype[v] == E_st)
	    {
	      printf("END state\n");
	    }
	  else
	    {
	      yoffset = ((char **) shadow[v])[jp][dp];
	      printf("shadow_banded_jd[%2d][%2d][%2d] : %d| j: %d | d: %d\n", v, jp, dp, yoffset, jp, dp);
	    }
	}
    }
}

/* EPN 03.29.06
   debug_print_alpha_banded_jd()
 * Function: debug_print_alpha_banded_jd
 *
 * Purpose:  Print alpha matrix banded in j and d dimensions
 */
void
debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
			    int **hdmin, int **hdmax)
{
  int v, j, d, dp, jp, max_v;

  printf("\nPrinting banded alpha matrix :\n");
  printf("************************************\n");
  max_v = cm->M-1;
  if(cm->flags & CM_LOCAL_BEGIN)
    {
      max_v = cm->M;
    }
  for(v = 0; v < max_v; v++)
    {
      printf("====================================\n");
      for(j = jmin[v]; j <= jmax[v]; j++)
	{
	  printf("------------------------------------\n");
	  for (d = hdmin[v][j-jmin[v]]; d <= hdmax[v][j-jmin[v]]; d++) 
	    {
	      jp = j - jmin[v]; // j index for state v in alpha w/mem eff bands
	      dp = d - hdmin[v][j-jmin[v]]; // d index for state v in alpha w/mem eff bands
	      printf("alpha_jd[%2d][%2d][%2d] : %6.2f | j: %4d | d: %4d\n", v, jp, dp, alpha[v][jp][dp], j, d);
	    }
	}
    }
  printf("****************\n\n");
}


/*****************************************************************************
 * EPN 11.03.05
 * Function: ij2d_bands()
 *
 * Purpose:  Determine the band for each cm state v on d (the band on the 
 *           length of the subsequence emitted from the subtree rooted
 *           at state v). These are easily calculated given the bands on i
 *           and j.
 * 
 * arguments:
 *
 * CM_t *cm         the CM 
 * int  W           length of sequence we're aligning
 * int *imin        imin[v] = first position in band on i for state v
 * int *imax        imax[v] = last position in band on i for state v
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 * int **hdmin      hdmin[v][j0] = first position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Filled in this function.
 * int **hdmax      hdmax[v][j0] = last position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Filled in this function.
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 *****************************************************************************/
void
ij2d_bands(CM_t *cm, int W, int *imin, int *imax, int *jmin, int *jmax,
	   int **hdmin, int **hdmax, int debug_level)
{
  int v;            /* counter over states of the CM */
  int j0;           /* counter over valid j's, but offset. j0+jmin[v] = actual j */
  int state_min_d;  /* minimum d allowed for a state, ex: MP_st = 2, ML_st = 1. etc. */
  for(v = 0; v < cm->M; v++)
    {
      if(cm->sttype[v] == E_st)
	{
	  for(j0 = 0; j0 <= (jmax[v]-jmin[v]); j0++)
	    {
	      hdmin[v][j0] = 0;
	      hdmax[v][j0] = 0;
	    }
	}
      else
	{
	  if((cm->sttype[v] == ML_st) ||
	     (cm->sttype[v] == MR_st) ||
	     (cm->sttype[v] == IL_st) ||
	     (cm->sttype[v] == IR_st))
	    state_min_d = 1;
	  else if(cm->sttype[v] == MP_st)
	    state_min_d = 2;
	  else
	    state_min_d = 0;

	  for(j0 = 0; j0 <= (jmax[v]-jmin[v]); j0++)
	    {
	      hdmin[v][j0] = (j0+jmin[v]) - imax[v] + 1;
	      hdmax[v][j0] = (j0+jmin[v]) - imin[v] + 1;
	      if(hdmin[v][j0] < state_min_d)
		{
		  hdmin[v][j0] = state_min_d;
		  /*printf("ERROR ij2dbands: v: %d | hdmin[act j: %d] : %d\n", v, (j0+jmin[v]), hdmin[v][j0]);*/
		}
	      if(hdmax[v][j0] < state_min_d)
		{
		  hdmax[v][j0] = state_min_d;
		  /*printf("ERROR ij2dbands: v: %d | hdmin[act j: %d] : %d\n", v, (j0+jmin[v]), hdmin[v][j0]);*/
		}
	      if(hdmax[v][j0] > W)
		{
		  hdmax[v][j0] = W;
		  /*printf("ERROR ij2dbands: v: %d | hdmax[act j: %d] : %d | W : %d\n", v, (j0+jmin[v]), hdmax[v][j0], W);*/
		}
	      if(debug_level == 2)
		{
		  printf("hd[%d][j=%d]: min: %d | max: %d\n", v, (j0+jmin[v]), hdmin[v][j0], hdmax[v][j0]);
		}
	    }
	}
    }
}


/*****************************************************************************
 * EPN 11.04.05
 * Function: hd2safe_hd_bands
 *
 * Purpose:  HMMERNAL milestone 4 function. Given 
 *           hdmin and hdmax 2D arrays, simply
 *           fill safe_hdmin and safe_hdmax (1D arrays):
 *           safe_hdmin[v] = min_d (hdmin[v][j0])
 *           safe_hdmax[v] = max_d (hdmax[v][j0])
 * 
 * arguments:
 * int M            num states in the CM.
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 * int **hdmin      hdmin[v][j0] = first position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 * int **hdmax      hdmax[v][j0] = last position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 * int *safe_hdmin  safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0)
 *                  filled in this function.
 * int *safe_hdmax  safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0)
 *                  filled in this function.
 *****************************************************************************/
void
hd2safe_hd_bands(int M, int *jmin, int *jmax, int **hdmin, int **hdmax,
		 int *safe_hdmin, int *safe_hdmax)

{
  int v;            /* counter over states of the CM */
  int j0;           /* counter over valid j's, but offset. j0+jmin[v] = actual j */

  for(v = 0; v < M; v++)
    {
      safe_hdmin[v] = hdmin[v][0];
      safe_hdmax[v] = hdmax[v][0];
      /*printf("j0: %2d | j: %2d | v: %3d | smin %d | smax : %d\n", 0, (jmin[v]), v, safe_hdmin[v], safe_hdmax[v]);*/
      for(j0 = 1; j0 <= (jmax[v]-jmin[v]); j0++)
	{
	  if(hdmin[v][j0] < safe_hdmin[v])
	    safe_hdmin[v] = hdmin[v][j0];
	  if(hdmax[v][j0] > safe_hdmax[v])
	    safe_hdmax[v] = hdmax[v][j0];
	  /*printf("j0: %2d | j: %2d | v: %3d | smin %d | smax : %d\n", j0, (j0+jmin[v]), v, safe_hdmin[v], safe_hdmax[v]);*/
	}
    }

  for(v = 0; v < M; v++)
    {
      for(j0 = 1; j0 <= (jmax[v]-jmin[v]); j0++)
	{
	  /* check to make sure that hdmax[v][j0] doesn't exceed j */
	  if(hdmax[v][j0] > (j0+jmin[v]))
	    {
	      printf("ERROR: hd2safe_bands.1\n");
	      /*exit(1);*/
	    }
	}
    }
}

/*****************************************************************************
 * EPN 01.18.06
 * Function: debug_print_hd_bands
 *
 * Purpose:  Print out the v and j dependent hd bands.
 * 
 *****************************************************************************/
void
debug_print_hd_bands(CM_t *cm, int **hdmin, int **hdmax, int *jmin, int *jmax)
{
  int v, j;
  char **sttypes;
  char **nodetypes;

  sttypes = malloc(sizeof(char *) * 10);
  sttypes[0] = "D";
  sttypes[1] = "MP";
  sttypes[2] = "ML";
  sttypes[3] = "MR";
  sttypes[4] = "IL";
  sttypes[5] = "IR";
  sttypes[6] = "S";
  sttypes[7] = "E";
  sttypes[8] = "B";
  sttypes[9] = "EL";

  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  printf("\nPrinting hd bands :\n");
  printf("****************\n");
  for(v = 0; v < cm->M; v++)
   {
     for(j = jmin[v]; j <= jmax[v]; j++) 
       {
	 printf("band v:%d j:%d n:%d %-4s %-2s min:%d max:%d\n", v, j, cm->ndidx[v], nodetypes[(int) cm->ndtype[(int) cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], hdmin[v][j-jmin[v]], hdmax[v][j-jmin[v]]);
       }
     printf("\n");
   }
  printf("****************\n\n");

  free(sttypes);
  free(nodetypes);

}


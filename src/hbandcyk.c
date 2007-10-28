/* 
 * hbandcyk.c (created using smallcyk.c as a template)
 * EPN 11.04.05
 * 
 * Alignment of a CM to a target sequence using HMM imposed bands
 * on the j and d dimensions of the CYK matrix.
 * These algorithms align to the entire target sequence
 * (e.g. global alignment). For sequence-local alignment, see
 * scancyk.c.
 * 
 * Also included here is CYKBandedScan_jd() which is 
 * a still experimental CYK scanning function that enforces
 * bands in the j AND d dimension.
 * 
 * NOTE: All of the alignment functions (not the scanner) are
 * memory-efficient in that they only allocate DP cells within bands
 * on j and d.  Initially less memory efficent functions were written,
 * they are kept at the end of this file for reference.
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

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_stack.h"

#include "funcs.h"
#include "structs.h"


/*******************************************************************************
 * 11.04.05
 * EPN 
 * Memory efficient banded versions of selected smallcyk.c functions that
 * enforce bands in the d and j dimensions informed by an HMM Forward/Backward
 * posterior decode of the target sequence.
 * 
 * These functions are modified from their originals in smallcyk.c to make 
 * HMM banded FULL (not D&C) CYK alignment memory efficient. The starting
 * point for CYKInside_b_jd() was CYKInside_b_me() in smallcyk.c. The main
 * difference is that bands in the j dimension are enforced, and the d
 * bands have j dependence. 
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

/* The alignment engine. 
 */
static float inside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, 
			    int r, int z, int i0, int j0, 
			    int do_full,
			    float ***alpha, float ****ret_alpha, 
			    struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			    void ****ret_shadow, 
			    int allow_begin, int *ret_b, float *ret_bsc,
			    int *jmin, int *jmax,
			    int **hdmin, int **hdmax,
			    int *safe_hdmin, int *safe_hdmax);

/* The traceback routine.
 */

static float insideT_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			     int r, int z, int i0, int j0, int allow_begin,
			     int *jmin, int *jax,
			     int **hdmin, int **hdmax,
			     int *safe_hdmin, int *safe_hdmax);

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


/* Function: PrintDPCellsSaved_jd()
 * Prints out an estimate of the speed up due to j and d bands */
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
	if (cm->sttype[v] != E_st) 
	  before += j + 1;
      for(j = jmin[v]; j <= jmax[v]; j++)
	if (cm->sttype[v] != E_st) 
	  {
	    max = (j < hdmax[v][j-jmin[v]]) ? j : hdmax[v][j-jmin[v]];
	    after += max - hdmin[v][j-jmin[v]] + 1;
	  }
    }
  printf("Before:  something like %.0f\n", before);
  printf("After:   something like %.0f\n", after);
  printf("Speedup: maybe %.2f fold\n\n", (float) before / (float) after);
}

/* Function: CYKInside_b_jd()
 *           EPN 11.04.05
 * based on CYKInside_b() which was based on CYKInside()
 *
 * Only difference is bands are used in d and j dimesions: 
 *
 * Date:     SRE, Sun Jun  3 19:48:33 2001 [St. Louis]
 *
 * Purpose:  Wrapper for the insideT_b_jd_me() routine - solve
 *           a full alignment problem, return the traceback
 *           and the score, without dividing & conquering, using bands.
 *           
 *           Analogous to CYKDivideAndConquer() in many respects;
 *           see the more extensive comments in that function for
 *           more details on shared aspects.
 *           
 * Args:     cm     - the covariance model
 *           sq     - the sequence, 1..L
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
CYKInside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, 
	       int *jmin, int *jmax, int **hdmin, int **hdmax, int *dmin, int *dmax)
{
  /* Contract check */
  if(dsq == NULL)
    esl_fatal("ERROR in CYKInside_b_jd(), dsq is NULL.\n");

  Parsetree_t *tr;
  int          z;
  float        sc;

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));
    printf("alignment strategy:CYKInside_b_jd:b:nosmall\n"); 
    printf("L: %d\n", L);*/

  /* Trust, but verify.
   * Check out input parameters.
   */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CM_LOCAL_BEGIN)) esl_fatal("internal error: we're not in local mode, but r is not root");
    if (cm->stid[r] != MATP_MP && cm->stid[r] != MATL_ML &&
	cm->stid[r] != MATR_MR && cm->stid[r] != BIF_B)
      esl_fatal("internal error: trying to do a local begin at a non-mainline start");
  }

  /* Create the parse tree, and initialize.
   */
  tr = CreateParsetree(100);
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

  /* Solve the whole thing with one call to insideT_b_jd.  This calls
     a memory efficient insideT function, which only allocates cells
     in alpha within the bands. 
   */
  sc += insideT_b_jd_me(cm, dsq, L, tr, r, z, i0, j0, (r==0), jmin, jmax, hdmin, hdmax, 
			dmin, dmax);

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  /*printf("returning from CYKInside_b_jd() sc : %f\n", sc); */

  return sc;
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
 *           The deck re-use strategy in general does not work with
 *           this implementation b/c each state has it's own j-specific
 *           bands. 
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
 *           In this banded version, there are more offset issues,
 *           these are detailed with comments in the code.
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
 *           sq        - the sequence [1..L]   
 *                     - length of the dsq
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
inside_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	       float ***alpha, float ****ret_alpha, 
	       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	       void ****ret_shadow, 
	       int allow_begin, int *ret_b, float *ret_bsc,
	       int *jmin, int *jmax, int **hdmin, int **hdmax,
	       int *safe_hdmin, int *safe_hdmax)
{
  /* Contract check */
  if(dsq == NULL)
    esl_fatal("ERROR in inside_b_jd_me(), dsq is NULL.\n");

  int      status;
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
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
  int      kp_z;           /* k (in the d dim) index for state z in alpha w/mem eff bands */
  int      Wp;             /* W also changes depending on state */
  int      jp_v, jp_y, jp_z;
  int      kmin, kmax;
  int      tmp_jmin, tmp_jmax;
  float  **tmp_deck;       /* temp variable, used only to free deckpool at end */

  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the sequence -- used in many loops */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* if caller didn't give us a matrix, make one.
   * It's important to allocate for M+1 decks (deck M is for EL, local
   * alignment) - even though Inside doesn't need EL, Outside does,
   * and we might reuse this memory in a call to Outside.  
   */
  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
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
    ESL_ALLOC(shadow, sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* Main recursion
   */
  for (v = vend; v >= vroot; v--) 
    {
      /* First we need a deck to fill in. With memory efficient bands 
       * we don't reuse decks b/c each state has different bands and therefore
       * different deck sizes, so we ALWAYS allocate a deck here.
       */
      alpha[v] = alloc_jdbanded_vjd_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
      //printf("allocated 
      if (cm->sttype[v] != E_st) {
	if (ret_shadow != NULL) {
	  if (cm->sttype[v] == B_st) {
	    kshad     = alloc_jdbanded_vjd_kshadow_deck(L, i0, j0, jmin[v], jmax[v], hdmin[v], hdmax[v]);
	    shadow[v] = (void **) kshad;
	  } else {
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
		    esl_fatal("band on E state %d has a non-zero d value within its j band for j:%d\n", v, j);
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
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][jp_v][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
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
		if (dsq[i] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
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
		if (dsq[j] < cm->abc->K)
		  alpha[v][jp_v][dp_v] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][jp_v][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		if (alpha[v][jp_v][dp_v] < IMPOSSIBLE) alpha[v][jp_v][dp_v] = IMPOSSIBLE;
	      }
	    }
	}
      /*if((cm->sttype[v] != IL_st) && (cm->sttype[v] != IR_st) && (cm->sttype[v] != B_st)) {
	for (j = jmin[v]; j <= jmax[v]; j++) { 
	  jp_v  = j - jmin[v];
	  i     = j - hdmin[v][jp_v] + 1;
	  for (dp_v = 0, d = hdmin[v][jp_v]; d <= hdmax[v][jp_v]; dp_v++, d++, i--) {
	    printf("alpha[v: %4d][jp_v: %4d][dp_v: %4d]: %.4f\n", v, jp_v, dp_v, alpha[v][jp_v][dp_v]);
	    
	  }
	  printf("\n");
	}
	printf("\n\n");
	}*/
  
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
       * Look at our children; if they're fully released, free them, we don't 
       * reuse decks with bands b/c each state has different deck size.
       */
      if (! do_full) {
	if (cm->sttype[v] == B_st) 
	  { 
	    /* we can definitely release the S children of a bifurc. */
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
		    free_vjd_deck(alpha[y], i0, j0);
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

 ERROR: 
  esl_fatal("Memory allocation error.\n");
  return 0.; /* never reached */
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
insideT_b_jd_me(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
		int r, int z, int i0, int j0, 
		int allow_begin, int *jmin, int *jmax,
		int **hdmin, int **hdmax,
		int *safe_hdmin, int *safe_hdmax)
{
  /* Contract check */
  if(dsq == NULL)
    esl_fatal("ERROR in insideT_b_jd_me(), dsq is NULL.");

  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;                /* stack that tracks bifurc parent of a right start */
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

  pda = esl_stack_ICreate();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  jp_v = j - jmin[v];
  dp_v = d - hdmin[v][jp_v];

  while (1) {
    if(cm->sttype[v] != EL_st && d > hdmax[v][jp_v])
      esl_fatal("ERROR in insideT_b_jd(). d : %d > hdmax[%d] (%d)\n", d, v, hdmax[v]);
    if(cm->sttype[v] != EL_st && d < hdmin[v][jp_v])
      esl_fatal("ERROR in insideT_b_jd(). d : %d < hdmin[%d] (%d)\n", d, v, hdmin[v]);
    
    if (cm->sttype[v] == B_st) {
      kp_z = ((int **) shadow[v])[jp_v][dp_v];   /* kp = offset len of right fragment */
      z = cm->cnum[v];
      jp_z = j-jmin[z];
      k = kp_z + hdmin[z][jp_z];  /* k = offset len of right fragment */
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */
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
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
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
      default:    esl_fatal("'Inconceivable!'\n'You keep using that word...'");
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
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */
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
  int     status;
  float **a;
  int     jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;
  /*printf("in alloc JD banded vjd deck, L : %d, i : %d, j : %d, jmin : %d, jmax : %d\n", L, i, j, jmin, jmax);*/

  if(j < jmin || i > jmax)
    esl_fatal("ERROR called alloc_jdbanded_vjd_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);

  ESL_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  ESL_ALLOC(a, sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
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
      ESL_DASSERT2(hdmax[jp-jmin] <= (jp+1))
      /* Based on my current understanding the above line should never be false, if it is means there's a valid d
       * in the hd band that is invalid because its > j. I think I check, or ensure, that this
       * doesn't happen when I'm constructing the d bands.
       */
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*a is offset only the first (jlast-jfirst+1) elements will be non-NULL*/
      ESL_ALLOC(a[jp-jfirst], sizeof(float) * bw);
      /*printf("\tallocated a[%d] | bw: %d\n", (jp-jfirst), bw);*/
    }
  return a;
 ERROR:
  esl_fatal("Memory allocation error.");
  return NULL; /* never reached */
}

char  **
alloc_jdbanded_vjd_yshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  int    status;
  char **a;
  int    jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;

  if(j < jmin || i > jmax)
    esl_fatal("ERROR called alloc_jdbanded_vjd_yshadow_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);

  ESL_ALLOC(a, sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  for (jp = (jlast-jfirst+1); jp <= L;     jp++) a[jp]     = NULL;

  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jp: %d | max : %d\n", jp, (jlast)); */
      ESL_DASSERT2(hdmax[jp-jmin] <= (jp+1))
      /* Based on my current understanding the above line should never be false, if it is means there's a valid d
       * in the hd band that is invalid because its > j. I think I check, or ensure, that this
       * doesn't happen when I'm constructing the d bands.
       */
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*printf("\tallocated a[%d]\n", (jp-jfirst));*/
      ESL_ALLOC(a[jp-jfirst], sizeof(char) * bw);
    }
  return a;
 ERROR:
  esl_fatal("Memory allocation error.");
  return NULL; /* never reached */
}
int** 
alloc_jdbanded_vjd_kshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax)
{
  int   status;
  int **a;
  int   jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/
  int     jfirst, jlast;

  if(j < jmin || i > jmax)
    esl_fatal("ERROR called alloc_jdbanded_vjd_kshadow_deck for i: %d j: %d which is outside the band on j, jmin: %d | jmax: %d\n", i, j, jmin, jmax);

  ESL_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  ESL_ALLOC(a, sizeof(float *) * (L+1));  /* always alloc 0..L rows, some of which are NULL */
  jfirst = ((i-1) > jmin) ? (i-1) : jmin;
  jlast = (j < jmax) ? j : jmax;
  for (jp = (jlast-jfirst+1); jp <= L;     jp++) a[jp]     = NULL;

  /* jfirst is the first valid j, jlast is the last */
  for (jp = jfirst; jp <= jlast; jp++)
    {
      /*printf("jp: %d | max : %d\n", jp, (jlast)); */
      ESL_DASSERT2(hdmax[jp-jmin] <= (jp+1))
      /* Based on my current understanding the above line should never be false, if it is means there's a valid d
       * in the hd band that is invalid because its > j. I think I check, or ensure, that this
       * doesn't happen when I'm constructing the d bands.
       */
      bw = hdmax[jp-jmin] - hdmin[jp-jmin] +1;

      /*printf("\tallocated a[%d]\n", (jp-jfirst));*/
      ESL_ALLOC(a[jp-jfirst], sizeof(int) * bw);
    }
  return a;
 ERROR:
  esl_fatal("Memory allocation error.");
  return NULL; /* never reached */
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
 * EPN, Thu Apr 26 13:27:16 2007
 * Function: combine_qdb_hmm_d_bands()
 *
 * Purpose:  Given hdmin and hdmax bands, and query dependent bands (QDBs)
 *           in cm->dmin and cm->dmax, combine them by redefining the 
 *           hdmin and hdmax bands where necessary:
 *           hdmin[v][j] = max(hdmin[v][j], dmin[v])
 *           hdmax[v][j] = min(hdmin[v][j], dmin[v])
 * 
 * arguments:
 *
 * CM_t *cm         the CM 
 * int *jmin        jmin[v] = first position in band on j for state v
 * int *jmax        jmax[v] = last position in band on j for state v
 * int **hdmin      hdmin[v][j0] = first position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Redefined in this function.
 * int **hdmax      hdmax[v][j0] = last position in band on d for state v
 *                                 and j position: j = j0+jmin[v].
 *                  Redefined in this function.
 *****************************************************************************/
void
combine_qdb_hmm_d_bands(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax)
{
  int v;            /* counter over states of the CM */
  int jp;           /* counter over valid j's, but offset. jp+jmin[v] = actual j */

  /* Contract check */
  if(!(cm->flags & CM_QDB))
    esl_fatal("ERROR, in combine_qdb_hmm_d_bands(), CM QDBs invalid.\n");
  if(cm->dmin == NULL || cm->dmax == NULL)
    esl_fatal("ERROR, in combine_qdb_hmm_d_bands() but cm->dmin and/or cm->dmax is NULL.\n");

  for(v = 0; v < cm->M; v++)
    {
      for(jp = 0; jp <= (jmax[v]-jmin[v]); jp++)
	{
	  hdmin[v][jp] = hdmin[v][jp] > cm->dmin[v] ? hdmin[v][jp] : cm->dmin[v];
	  hdmax[v][jp] = hdmax[v][jp] < cm->dmax[v] ? hdmax[v][jp] : cm->dmax[v];
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
  int status;
  int v, j;
  char **sttypes;
  char **nodetypes;

  ESL_ALLOC(sttypes, sizeof(char *) * 10);
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

  ESL_ALLOC(nodetypes, sizeof(char *) * 8);
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
  return;

 ERROR:
  esl_fatal("Memory allocation error.");
}

/* Function: CYKBandedScan_jd() 
 * Date:     EPN, 06.22.06
 * based on: CYKBandedScan() by SRE in bandcyk.c
 *
 * Purpose:  Scan a (sub)sequence for matches to a covariance model, using HMM
 *           derived bands on the j and d dimensions. Intended for use on 
 *           subsequences with endpoints i and j, where i and j were determined
 *           using a HMM scan of the full sequence. This function then refines
 *           the positions of i and j, as well as deriving a CYK score that is
 *           more informative than an HMM based score. 
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *
 *           jmin, jmax set the state specific bounds on the j dimension. 0..v..cm->M-1.
 *           hdmin, hdmax set the state and j specific bounds on the d dimension, indexed
 *           [0..v..cm-M-1][0..(jmax[v]-jmin[v]+1)].
 *           
 *           The j band for v is jmin[v]..jmax[v], inclusive; that is,
 *           jmin[v] is the minimum allowed j for state v;
 *           jmax[v] is the maximum; 
 * 
 *           The d bands are v and j specific (in contrast to the a priori d bands
 *           which are only v specific), the d band for v and j is 
 *           hdmin[v][j-jmin[v]]..hdmax[v][j-jmin[v]] inclusive;
 *           hdmin[v][j-jmin[v]] is the minimum allowed d for state v and end point j
 *           hdmax[v][j-jmin[v]] is the maximum; 
 *
 * Args:     cm        - the covariance model
 *           sq        - digitized sequence to search; i0..j0
 *           jmin      - minimum bound on j for state v; 0..M
 *           jmax      - maximum bound on j for state v; 0..M
 *           hdmin     - minimum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           hdmax     - maximum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results   - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float
CYKBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
		 int j0, int W, float cutoff, search_results_t *results)
{
  int       status;
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  float  ***alpha_mem;          /* pointers to original alpha memory */
  float    *imp_row;           /* an IMPOSSIBLE deck (full of IMPOSSIBLE scores), 
				 * pointed to when j is outside j band for v */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int      gamma_j;             /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int      gamma_i;             /* i index in the gamma matrix */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp_roll;   	        /* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       tmp_dmin, tmp_dmax; /* temp variables for ensuring we stay within d bands within loops */
  int       tmp_kmin, tmp_kmax; /* temp vars for B_st's, min/max k values consistent with bands*/

  int      jp_v, jp_y, jp_w;    /* mem eff banded j index in states v, y, and z 
				 * jp_x = j-jmin[x] */
  int      L;                   /* length of subsequence (j0-i0+1) */
  int      x;
  int      tmp_y;
  float    tmp_alpha_w, tmp_alpha_y;
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  
  /* Contract checks */
  if((!(cm->search_opts & CM_SEARCH_NOQDB)) && (cm->dmin == NULL || cm->dmax == NULL))
    esl_fatal("ERROR in CYKBandedScan_jd(), trying to use QDB, but cm->dmin, cm->dmax are NULL.\n");
  if(dsq == NULL)
    esl_fatal("ERROR in CYKBandedScan_jd(), dsq is NULL.");

  if(!(cm->search_opts & CM_SEARCH_NOQDB)) /* we're doing qdb */
    combine_qdb_hmm_d_bands(cm, jmin, jmax, hdmin, hdmax);

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  /*printf("in CYKBandedScan_jd i0: %5d | j0: %5d | L: %5d | W: %5d\n", i0, j0, L, W);*/
  if (W > L) W = L; /* shouldn't look longer than seq length L */

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));*/
  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that unlike the other CYK scan functions, E memory is not shared: 
   * this is because the E deck will be different for different j values
   * due to the j bands. 
   * 
   *****************************************************************/
  ESL_ALLOC(alpha, sizeof(float **) * cm->M);
  ESL_ALLOC(alpha_mem, sizeof(float **) * cm->M);
  /* we use alpha_mem to remember where each alpha row (alpha[v][j]) is
   * in case we've set alpha[v][cur] to imp_row (the precalc'ed IMPOSSIBLE row)
   * in a prior iteration, and we are about to overwrite it, and we don't
   * want to overwrite our special IMPOSSIBLE row.
   */
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], sizeof(float *)  * (W+1));
	ESL_ALLOC(alpha_mem[v], sizeof(float *)  * (W+1));
	for (j = 0; j <= W; j++)
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(float) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	ESL_ALLOC(alpha_mem[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(float) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

  /* Initialize the impossible deck, which we'll point to for 
   * j positions that are outside of the j band on v */
  ESL_ALLOC(imp_row, sizeof(float) * (W+1));
  for (d = 0; d <= W; d++) imp_row[d] = IMPOSSIBLE;
    
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {

      gamma_j = j-i0+1;
      cur = (j-i0+1)%2; /* cur == 1 when j == i0 */
      prv = (j-i0)  %2; /* prv == 0 when j == i0 */

      /*****************************************************************
       * alpha initializations.
       * For the jd (HMM) banded strategy, we initialize inside the j loop,
       * because no cells are j-independent: for j's outside
       * the bands for a state v, should have ALL cells = IMPOSSIBLE.
       *****************************************************************/ 
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* Check to see if we're outside the bounds on j */
	  if(j < jmin[v] || j > jmax[v])
	    {
	      if (cm->stid[v] == BEGL_S) 
		{
		  jp_roll = j % (W+1); 
		  for (d = 0; d <= W; d++) 
		    alpha[v][jp_roll][d] = IMPOSSIBLE;
		}
	      else
		{
		  jp_roll = cur;
		  alpha[v][jp_roll] = imp_row;
		}
	      /* Special boundary case: have to initialize alpha[v][prv] also */
	      if (j == i0)
		alpha[v][prv] = imp_row;
	      /*for (d = 0; d <= W; d++) 
		alpha[v][prv][d] = IMPOSSIBLE;*/

	      continue;
	    }
	  if(alpha[v][cur] == imp_row)
	    {
	      /* we don't want to overwrite imp_row */
	      alpha[v][cur] = alpha_mem[v][cur];
	      for (d = 0; d <= W; d++) 
		alpha[v][cur][d] = IMPOSSIBLE;
	    }
	  
	  /* else we initialize on d = 0 */
	  alpha[v][cur][0] = IMPOSSIBLE;
	  
	  if      (cm->sttype[v] == E_st)  alpha[v][cur][0] = 0;
	  else if (cm->sttype[v] == MP_st) alpha[v][cur][1] = alpha[v][prv][1] = IMPOSSIBLE;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[v][cur][0] = cm->endsc[v];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		if ((sc = alpha[y+yoffset][cur][0] + cm->tsc[v][yoffset]) > alpha[v][cur][0]) 
		  alpha[v][cur][0] = sc;
	      /* ...we don't bother to look at local alignment starts here... */
	      bestr[cur] = -1;
	      if (alpha[v][cur][0] < IMPOSSIBLE) alpha[v][cur][0] = IMPOSSIBLE;	
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* w is BEGL_S */
	      y = cm->cnum[v];   /* y is BEGR_S */
	      /* original line: 
	       * alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	       * we can't use that because alpha[w][0][0] and alpha[y][0][0] 
	       * may have been reset to IMPOSSIBLE, so we recalculate what they
	       * should be (this is wasteful):
	       */
	      tmp_y = cm->cfirst[w];
	      tmp_alpha_w = cm->endsc[w];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[w]; yoffset++)
		{
		  if ((sc = alpha[tmp_y+yoffset][cur][0] + cm->tsc[w][yoffset]) > tmp_alpha_w)
		    tmp_alpha_w = sc;
		}
	      tmp_y = cm->cfirst[y];
	      tmp_alpha_y = cm->endsc[y];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[y]; yoffset++)
		if ((sc = alpha[tmp_y+yoffset][cur][0] + cm->tsc[y][yoffset]) > tmp_alpha_y)
		  tmp_alpha_y = sc;
	      alpha[v][cur][0] = tmp_alpha_w + tmp_alpha_y;
	      if (alpha[v][cur][0] < IMPOSSIBLE) alpha[v][cur][0] = IMPOSSIBLE;	
	    }
	  
	  /* Special boundary case: have to initialize alpha[v][prv] also */
	  if(j == i0) alpha[v][prv][0] = alpha[v][cur][0];

	  if (cm->stid[v] == BEGL_S) 
	    {
	      alpha[v][prv][0] = alpha[v][cur][0];
	      for (x = 2; x <= W; x++) 
		alpha[v][x][0] = alpha[v][cur][0];
	    }
	  /* done initialization */
	  
	  jp_v = j - jmin[v];
	  /* Impose the bands.
	   *   We have to do this inside the main loop because d bands are
	   *   dependent on v AND j. 
	   */
	  if (cm->stid[v] == BEGL_S) jp_roll = j % (W+1); else jp_roll = cur;

	  for (d =0; d < hdmin[v][jp_v] && d <=W; d++) 
	    alpha[v][jp_roll][d] = IMPOSSIBLE;
	  for (d = hdmax[v][jp_v]+1; d <= W;      d++) 
	    alpha[v][jp_roll][d] = IMPOSSIBLE;
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp_roll][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[v][yoffset]) > alpha[v][jp_roll][d]) 
		      alpha[v][jp_roll][d] = sc;
		  if (alpha[v][jp_roll][d] < IMPROBABLE) alpha[v][jp_roll][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      /* Five inequalities must be satisfied to ensure that j and k 
	       * and k combinations correspond with alpha cells within the bands 
	       * on states y and w. 
	       * Below: jp_y = j - jmin[y] & jp_w = j - jmin[w]
	       *
	       * (1) j   >= jmin[y]          && j   <= jmax[y]
	       * (2) j-k >= jmin[w]          && j-k <= jmax[w]
	       * (3) k   >= hdmin[y][jp_y]   && k   <= hdmax[y][jp_y]
	       * (4) d-k >= hdmin[w][jp_w-k] && d-k <= hdmax[w][jp_w-k]
	       * (5) d   >= hdmin[v][jp_v]   && d   <= hdmax[v][jp_v]
	       */

	      /* initialize with endsc for all valid d for B state v */
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		}
	      /* Following code is careful, and not 'efficient' */
	      if(j >= jmin[y] && j <= jmax[y]) /* ensures (1): that j is valid for state y */
		{
		  jp_y = j - jmin[y];
		  jp_w = j - jmin[w]; 
		  i = j-d+1;
		  /*
		    printf("valid j: %d | jp_y: %d | jp_w: %d\n", j, jp_y, jp_w);
		    printf("hdmin[v][jp_v]: %d | hdmin[y][jp_y]: %d\n", hdmin[v][jp_v], hdmin[y][jp_y]);
		    printf("hdmax[v][jp_v]: %d | hdmax[y][jp_y]: %d\n", hdmax[v][jp_v], hdmax[y][jp_y]);
		  */
		  for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		    {
		      /* k is the length of the right fragment */
		      tmp_kmin = ((j-jmax[w]) > hdmin[y][jp_y]) ? (j-jmax[w]) : hdmin[y][jp_y];
		      if(tmp_kmin < 0) tmp_kmin = 0;
		      /* HEREHEREHEREHEREHEREHEREHERE, ensure that hdmax[w][jp_w-tmp_kmin] is valid
		       * before accessing it! */
		      if(tmp_kmin < d-hdmax[w][jp_w-tmp_kmin]) tmp_kmin = d-hdmax[w][jp_w-tmp_kmin];
		      /* tmp_kmin is now smallest k that satisfies (2), (3), and (4) */

		      tmp_kmax = ((j-jmin[w]) < hdmax[y][jp_y]) ? (j-jmin[w]) : hdmax[y][jp_y];
		      /* HEREHEREHEREHEREHEREHEREHERE, ensure that hdmin[w][jp_w-tmp_kmax] is valid
		       * before accessing it! */
		      if(tmp_kmax > d-hdmin[w][jp_w-tmp_kmax]) tmp_kmax = d-hdmin[w][jp_w-tmp_kmax];
		      /* tmp_kmax is now largest k that satisfies (2), (3), and (4) */
		      /*printf("tmp_kmin: %d | tmp_kmax: %d\n", tmp_kmin, tmp_kmax);*/
		      for (k = tmp_kmin; k <= tmp_kmax; k++)
			{
			  jp_roll = (j-k)%(W+1); /* jp_roll is rolling index into BEGL_S (state w) 
						  * deck j dimension */
			  if ((sc = alpha[w][jp_roll][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			    alpha[v][cur][d] = sc;
			}
		      if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		      /*printf("B alpha[%d][%d][%d]: %f\n", v, cur, d, alpha[v][cur][d]);*/
		    }
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      /* Check to see if we're within bounds on j */
      if(j < jmin[0] || j > jmax[0])
	{
	  /*printf("j: %d (gamma_j: %d) IMPOSSIBLE BABY! min: %d max: %d\n", j, gamma_j, jmin[0], jmax[0]);*/
	  for (d = 0; d <= W; d++) 
	    alpha[0][jp_roll][d] = IMPOSSIBLE; 
	  /* Inform the little semi-Markov model that deals with multihit parsing
	   * that a hit is impossible, j is outside root band on j:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  continue;
	}
      /* if we get here, j is within ROOT_S state 0's band */

      /* first initialize on d = 0 */
      alpha[0][0][0] = IMPOSSIBLE;
      y = cm->cfirst[v];
      alpha[0][0][0] = cm->endsc[v];
      /* treat EL as emitting only on self transition */
      for (yoffset = 0; yoffset < cm->cnum[0]; yoffset++)
	if ((sc = alpha[y+yoffset][0][0] + cm->tsc[0][yoffset]) > alpha[0][0][0]) 
	  alpha[0][0][0] = sc;
      /* ...we don't bother to look at local alignment starts here... */
      bestr[0] = -1;
      if (alpha[0][0][0] < IMPOSSIBLE) alpha[0][0][0] = IMPOSSIBLE;	
      alpha[0][1][0] = alpha[0][0][0];
      /* done initialization on d = 0 */

      jp_v = j - jmin[0];
      /* Impose the bands.
       *   We have to do this here because d bands are
       *   dependent on v AND j. 
       */
      for (d =0; d < hdmin[0][jp_v] && d <=W; d++) 
	alpha[0][cur][d] = IMPOSSIBLE;
      for (d = hdmax[0][jp_v]+1; d <= W;      d++) 
	alpha[0][cur][d] = IMPOSSIBLE;
      
      for (d = hdmin[0][jp_v]; ((d <= hdmax[0][jp_v] && d <= gamma_j) && d <= W); d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    {
	      if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
		{
		  alpha[0][cur][d] = sc;
		}
	    }
	  /*printf("j: %d | alpha[0][cur][%d]: %f\n", j, d, alpha[0][cur][d]);*/
	  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	}

      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(j >= jmin[y] && j <= jmax[y]) 
	    {
	      jp_y = j - jmin[y];
	      tmp_dmin = (hdmin[y][jp_y] > hdmin[0][jp_v]) ? hdmin[y][jp_y] : hdmin[0][jp_v];
	      tmp_dmax = (hdmax[y][jp_y] < hdmax[0][jp_v]) ? hdmax[y][jp_y] : hdmax[0][jp_v];
	      if(tmp_dmax > j) tmp_dmax = j;
	      for (d = tmp_dmin; d <= tmp_dmax; d++)
		{
		  if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
		  else                       sc = alpha[y][cur][d]     + cm->beginsc[y]; /* BUG! alpha[y][cur][d] outside bands */
		  if (sc > alpha[0][cur][d]) {
		    alpha[0][cur][d] = sc;
		    bestr[d]         = y;
		  }
		  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
		  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
		}
	    }
	}
      }
	      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = hdmin[0][jp_v]; (d <= hdmax[0][jp_v] && d <= gamma_j) && d <= W; d++) 
	{
	  i       = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  assert(i > 0);
	  /*printf("v: %d gamma_i: %d d: %d\n", v, gamma_i, d);*/
	  /*printf("alpha[0][j:%3d][d:%3d]: %f\n", j, d, alpha[0][cur][d]);*/
	  sc = gamma[gamma_i-1] + alpha[0][cur][d] + cm->sc_boost;
	  /* sc_boost is experimental technique for finding hits < 0 bits. 
	   * value is 0.0 if technique not used. */
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = alpha[0][cur][d]; 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */

  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   * be careful about only freeing our impossible deck, imp_row, once
   *****************************************************************/ 
for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) 
	  if(alpha_mem[v][j] != imp_row) free(alpha_mem[v][j]); 
	free(alpha[v]);
	free(alpha_mem[v]);
      } else {
	if(alpha_mem[v][0] != imp_row) free(alpha_mem[v][0]);
	if(alpha_mem[v][1] != imp_row) free(alpha_mem[v][1]); 
	free(alpha[v]);
	free(alpha_mem[v]);
      }
    }
  free(imp_row);
  free(alpha);
  free(alpha_mem);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  return best_score;
 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}

/* Function: iInsideBandedScan_jd() 
 * Date    : EPN, Fri Apr 27 09:32:38 2007
 *           
 *           
 * Purpose:  Identical to CYKBandedScan_jd(), but sums replaces maxes.
 *           Scan a (sub)sequence for matches to a covariance model, using HMM
 *           derived bands on the j and d dimensions. Intended for use on 
 *           subsequences with endpoints i and j, where i and j were determined
 *           using a HMM scan of the full sequence. This function then refines
 *           the positions of i and j, as well as deriving a CYK score that is
 *           more informative than an HMM based score. 
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *           Log sums are performed using scaled ints with ILogsum().
 *
 *           jmin, jmax set the state specific bounds on the j dimension. 0..v..cm->M-1.
 *           hdmin, hdmax set the state and j specific bounds on the d dimension, indexed
 *           [0..v..cm-M-1][0..(jmax[v]-jmin[v]+1)].
 *           
 *           The j band for v is jmin[v]..jmax[v], inclusive; that is,
 *           jmin[v] is the minimum allowed j for state v;
 *           jmax[v] is the maximum; 
 * 
 *           The d bands are v and j specific (in contrast to the a priori d bands
 *           which are only v specific), the d band for v and j is 
 *           hdmin[v][j-jmin[v]]..hdmax[v][j-jmin[v]] inclusive;
 *           hdmin[v][j-jmin[v]] is the minimum allowed d for state v and end point j
 *           hdmax[v][j-jmin[v]] is the maximum; 
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; i0..j0
 *           jmin      - minimum bound on j for state v; 0..M
 *           jmax      - maximum bound on j for state v; 0..M
 *           hdmin     - minimum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           hdmax     - maximum bound on j for state v and end posn j;
 *                       [0..M-1][0..(jmax[v]-jmin[v]+1)          
 *           i0        - start of target subsequence (1 for beginning of dsq)
 *           j0        - end of target subsequence (L for end of dsq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results   - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float
iInsideBandedScan_jd(CM_t *cm, ESL_DSQ *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
		    int j0, int W, float cutoff, search_results_t *results)
{
  int        status;
  int     ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int     ***alpha_mem;          /* pointers to original alpha memory */
  int       *imp_row;           /* an IMPOSSIBLE deck (full of -INFTY scores), 
				 * pointed to when j is outside j band for v */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int      gamma_j;             /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int      gamma_i;             /* i index in the gamma matrix */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp_roll;   	        /* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       tmp_dmin, tmp_dmax; /* temp variables for ensuring we stay within d bands within loops */
  int       tmp_kmin, tmp_kmax; /* temp vars for B_st's, min/max k values consistent with bands*/

  int      jp_v, jp_y, jp_w;    /* mem eff banded j index in states v, y, and z 
				 * jp_x = j-jmin[x] */
  int      L;                   /* length of subsequence (j0-i0+1) */
  int      x;
  int      tmp_y;
  int      tmp_alpha_w, tmp_alpha_y;
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  
  /* Contract checks */
  if((!(cm->search_opts & CM_SEARCH_NOQDB)) && (cm->dmin == NULL || cm->dmax == NULL))
    esl_fatal("ERROR in iInsideBandedScan_jd(), trying to use QDB, but cm->dmin, cm->dmax are NULL.\n");
  if(dsq == NULL)
    esl_fatal("ERROR in iInsideBandedScan_jd(), dsq is NULL.");

  if(!(cm->search_opts & CM_SEARCH_NOQDB)) /* we're doing qdb */
    combine_qdb_hmm_d_bands(cm, jmin, jmax, hdmin, hdmax);

  best_score     = -INFTY;
  best_neg_score = -INFTY;
  L = j0-i0+1;
  /*printf("in iInsideBandedScan_jd i0: %5d | j0: %5d | L: %5d | W: %5d\n", i0, j0, L, W);*/
  if (W > L) W = L; /* shouldn't look longer than seq length L */

  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (j0-i0+1));*/
  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that unlike the other CYK scan functions, E memory is not shared: 
   * this is because the E deck will be different for different j values
   * due to the j bands. 
   * 
   *****************************************************************/
  ESL_ALLOC(alpha,     sizeof(int **) * cm->M);
  ESL_ALLOC(alpha_mem, sizeof(int **) * cm->M);
  /* we use alpha_mem to remember where each alpha row (alpha[v][j]) is
   * in case we've set alpha[v][cur] to imp_row (the precalc'ed -INFTY row)
   * in a prior iteration, and we are about to overwrite it, and we don't
   * want to overwrite our special -INFTY row.
   */
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v],     sizeof(int *)  * (W+1));
	ESL_ALLOC(alpha_mem[v], sizeof(int *)  * (W+1));
	for (j = 0; j <= W; j++)
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(int) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
    else 
      {
	ESL_ALLOC(alpha[v],     sizeof(int *) * 2);
	ESL_ALLOC(alpha_mem[v], sizeof(int *) * 2);
	for (j = 0; j < 2; j++) 
	  {
	    ESL_ALLOC(alpha_mem[v][j], sizeof(int) * (W+1));
	    alpha[v][j]     = alpha_mem[v][j];
	  }
      }
  }
  ESL_ALLOC(bestr, sizeof(int) * (W+1));

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  ESL_ALLOC(gamma, sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback, sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc,sizeof(float) * (L+1));
  ESL_ALLOC(saver, sizeof(int)   * (L+1));

  /* Initialize the impossible deck, which we'll point to for 
   * j positions that are outside of the j band on v */
  ESL_ALLOC(imp_row, sizeof(int) * (W+1));
  for (d = 0; d <= W; d++) imp_row[d] = -INFTY;
    
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {

      gamma_j = j-i0+1;
      cur = (j-i0+1)%2; /* cur == 1 when j == i0 */
      prv = (j-i0)  %2; /* prv == 0 when j == i0 */

      /*****************************************************************
       * alpha initializations.
       * For the jd (HMM) banded strategy, we initialize inside the j loop,
       * because no cells are j-independent: for j's outside
       * the bands for a state v, should have ALL cells = -INFTY.
       *****************************************************************/ 
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  /* Check to see if we're outside the bounds on j */
	  if(j < jmin[v] || j > jmax[v])
	    {
	      if (cm->stid[v] == BEGL_S) 
		{
		  jp_roll = j % (W+1); 
		  for (d = 0; d <= W; d++) 
		    alpha[v][jp_roll][d] = -INFTY;
		}
	      else
		{
		  jp_roll = cur;
		  alpha[v][jp_roll] = imp_row;
		}
	      /* Special boundary case: have to initialize alpha[v][prv] also */
	      if (j == i0)
		alpha[v][prv] = imp_row;
	      /*for (d = 0; d <= W; d++) 
		alpha[v][prv][d] = -INFTY;*/

	      continue;
	    }
	  if(alpha[v][cur] == imp_row)
	    {
	      /* we don't want to overwrite imp_row */
	      alpha[v][cur] = alpha_mem[v][cur];
	      for (d = 0; d <= W; d++) 
		alpha[v][cur][d] = -INFTY;
	    }
	  
	  /* else we initialize on d = 0 */
	  alpha[v][cur][0] = -INFTY;
	  
	  if      (cm->sttype[v] == E_st)  alpha[v][cur][0] = 0;
	  else if (cm->sttype[v] == MP_st) alpha[v][cur][1] = alpha[v][prv][1] = -INFTY;
	  else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	    {
	      y = cm->cfirst[v];
	      alpha[v][cur][0] = cm->iendsc[v];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		alpha[v][cur][0] = ILogsum(alpha[v][cur][0], (alpha[y+yoffset][cur][0] 
							      + cm->itsc[v][yoffset]));
	      /* ...we don't bother to look at local alignment starts here... */
	      bestr[cur] = -1;
	      if (alpha[v][cur][0] < -INFTY) alpha[v][cur][0] = -INFTY;	
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v]; /* w is BEGL_S */
	      y = cm->cnum[v];   /* y is BEGR_S */
	      /* original line: 
	       * alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	       * we can't use that because alpha[w][0][0] and alpha[y][0][0] 
	       * may have been reset to -INFTY, so we recalculate what they
	       * should be (this is wasteful):
	       */
	      tmp_y = cm->cfirst[w];
	      tmp_alpha_w = cm->iendsc[w];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[w]; yoffset++)
		{
		  tmp_alpha_w = ILogsum(tmp_alpha_w, (alpha[tmp_y+yoffset][cur][0] 
						      + cm->itsc[w][yoffset]));
		}
	      tmp_y = cm->cfirst[y];
	      tmp_alpha_y = cm->iendsc[y];
	      /* treat EL as emitting only on self transition */
	      for (yoffset = 0; yoffset < cm->cnum[y]; yoffset++)
		tmp_alpha_y = ILogsum(tmp_alpha_y, (alpha[tmp_y+yoffset][cur][0]
						    + cm->itsc[y][yoffset]));
	      alpha[v][cur][0] = tmp_alpha_w + tmp_alpha_y;
	      /*printf("! alpha[v][j:%3d][d:%3d]: %f\n", v, j, 0, Scorify(alpha[v][cur][0]));*/
	      if (alpha[v][cur][0] < -INFTY) alpha[v][cur][0] = -INFTY;	
	    }
	  
	  /* Special boundary case: have to initialize alpha[v][prv] also */
	  if(j == i0) alpha[v][prv][0] = alpha[v][cur][0];

	  if (cm->stid[v] == BEGL_S) 
	    {
	      alpha[v][prv][0] = alpha[v][cur][0];
	      for (x = 2; x <= W; x++) 
		alpha[v][x][0] = alpha[v][cur][0];
	    }
	  /* done initialization */
	  
	  jp_v = j - jmin[v];
	  /* Impose the bands.
	   *   We have to do this inside the main loop because d bands are
	   *   dependent on v AND j. 
	   */
	  if (cm->stid[v] == BEGL_S) jp_roll = j % (W+1); else jp_roll = cur;

	  for (d =0; d < hdmin[v][jp_v] && d <=W; d++) 
	    alpha[v][jp_roll][d] = -INFTY;
	  for (d = hdmax[v][jp_v]+1; d <= W;      d++) 
	    alpha[v][jp_roll][d] = -INFTY;
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp_roll][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][jp_roll][d] = ILogsum(alpha[v][jp_roll][d], (alpha[y+yoffset][cur][d] + 
									  cm->itsc[v][yoffset]));
		  if (alpha[v][jp_roll][d] < -INFTY) alpha[v][jp_roll][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-2] + 
								  cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += iDegeneratePairScore(cm->abc, cm->iesc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][cur][d-1] + 
								  cm->itsc[v][yoffset]));
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[i], cm->iesc[v]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		  /*printf("alpha[v][j:%3d][d:%3d]: %f\n", v, j, d, Scorify(alpha[v][cur][d]));*/

		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[y+yoffset][prv][d-1] + 
								  cm->itsc[v][yoffset]));
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->iesc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_IAvgScore(cm->abc, dsq[j], cm->iesc[v]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      /* Five inequalities must be satisfied to ensure that j and k 
	       * and k combinations correspond with alpha cells within the bands 
	       * on states y and w. 
	       * Below: jp_y = j - jmin[y] & jp_w = j - jmin[w]
	       *
	       * (1) j   >= jmin[y]          && j   <= jmax[y]
	       * (2) j-k >= jmin[w]          && j-k <= jmax[w]
	       * (3) k   >= hdmin[y][jp_y]   && k   <= hdmax[y][jp_y]
	       * (4) d-k >= hdmin[w][jp_w-k] && d-k <= hdmax[w][jp_w-k]
	       * (5) d   >= hdmin[v][jp_v]   && d   <= hdmax[v][jp_v]
	       */

	      /* initialize with endsc for all valid d for B state v */
	      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		{
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		}
	      /* Following code is careful, and not 'efficient' */
	      if(j >= jmin[y] && j <= jmax[y]) /* ensures (1): that j is valid for state y */
		{
		  jp_y = j - jmin[y];
		  jp_w = j - jmin[w]; 
		  i = j-d+1;
		  /*
		    printf("valid j: %d | jp_y: %d | jp_w: %d\n", j, jp_y, jp_w);
		    printf("hdmin[v][jp_v]: %d | hdmin[y][jp_y]: %d\n", hdmin[v][jp_v], hdmin[y][jp_y]);
		    printf("hdmax[v][jp_v]: %d | hdmax[y][jp_y]: %d\n", hdmax[v][jp_v], hdmax[y][jp_y]);
		  */
		  for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) /* ensures (5) above */
		    {
		      /* k is the length of the right fragment */
		      tmp_kmin = ((j-jmax[w]) > hdmin[y][jp_y]) ? (j-jmax[w]) : hdmin[y][jp_y];
		      if(tmp_kmin < 0) tmp_kmin = 0;
		      if(tmp_kmin < d-hdmax[w][jp_w-tmp_kmin]) tmp_kmin = d-hdmax[w][jp_w-tmp_kmin];
		      /* tmp_kmin is now smallest k that satisfies (2), (3), and (4) */

		      tmp_kmax = ((j-jmin[w]) < hdmax[y][jp_y]) ? (j-jmin[w]) : hdmax[y][jp_y];
		      if(tmp_kmax > d-hdmin[w][jp_w-tmp_kmax]) tmp_kmax = d-hdmin[w][jp_w-tmp_kmax];
		      /* tmp_kmax is now largest k that satisfies (2), (3), and (4) */
		      /*printf("tmp_kmin: %d | tmp_kmax: %d\n", tmp_kmin, tmp_kmax);*/
		      for (k = tmp_kmin; k <= tmp_kmax; k++)
			{
			  jp_roll = (j-k)%(W+1); /* jp_roll is rolling index into BEGL_S (state w) 
						  * deck j dimension */
			  alpha[v][cur][d] = ILogsum(alpha[v][cur][d], (alpha[w][jp_roll][d-k] + 
									alpha[y][cur][k]));
			}
		      if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		      /*printf("B alpha[%d][%d][%d]: %f\n", v, cur, d, alpha[v][cur][d]);*/
		    }
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * IMPOSSIBLE), which is probably a little too fragile of a method. 
       */
      /* Check to see if we're within bounds on j */
      if(j < jmin[0] || j > jmax[0])
	{
	  /*printf("j: %d (gamma_j: %d) -INFTY BABY! min: %d max: %d\n", j, gamma_j, jmin[0], jmax[0]);*/
	  for (d = 0; d <= W; d++) 
	    alpha[0][jp_roll][d] = -INFTY; 
	  /* Inform the little semi-Markov model that deals with multihit parsing
	   * that a hit is impossible, j is outside root band on j:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  continue;
	}
      /* if we get here, j is within ROOT_S state 0's band */

      /* first initialize on d = 0 */
      alpha[0][0][0] = -INFTY;
      y = cm->cfirst[v];
      alpha[0][0][0] = cm->iendsc[v];
      /* treat EL as emitting only on self transition */
      for (yoffset = 0; yoffset < cm->cnum[0]; yoffset++)
	alpha[0][0][0] = ILogsum(alpha[0][0][0], (alpha[y+yoffset][0][0]
						  + cm->itsc[0][yoffset]));
      /* ...we don't bother to look at local alignment starts here... */
      bestr[0] = -1;
      if (alpha[0][0][0] < -INFTY) alpha[0][0][0] = -INFTY;	
      alpha[0][1][0] = alpha[0][0][0];
      /* done initialization on d = 0 */

      jp_v = j - jmin[0];
      /* Impose the bands.
       *   We have to do this here because d bands are
       *   dependent on v AND j. 
       */
      for (d =0; d < hdmin[0][jp_v] && d <=W; d++) 
	alpha[0][cur][d] = -INFTY;
      for (d = hdmax[0][jp_v]+1; d <= W;      d++) 
	alpha[0][cur][d] = -INFTY;
      
      for (d = hdmin[v][jp_v]; ((d <= hdmax[v][jp_v] && d <= gamma_j) && d <= W); d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->itsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    {
	      alpha[0][cur][d] = ILogsum(alpha[0][cur][d], (alpha[y+yoffset][cur][d]
							    + cm->itsc[0][yoffset]));
	    }
	  /*printf("j: %d | alpha[0][cur][%d]: %f\n", j, d, alpha[0][cur][d]);*/
	  if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
	  if (Scorify(alpha[0][cur][d]) > best_neg_score) best_neg_score = Scorify(alpha[0][cur][d]);
	}

      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  if(j >= jmin[y] && j <= jmax[y]) 
	    {
	      jp_y = j - jmin[y];
	      tmp_dmin = (hdmin[y][jp_y] > hdmin[0][jp_v]) ? hdmin[y][jp_y] : hdmin[0][jp_v];
	      tmp_dmax = (hdmax[y][jp_y] < hdmax[0][jp_v]) ? hdmax[y][jp_y] : hdmax[0][jp_v];
	      if(tmp_dmax > j) tmp_dmax = j;
	      for (d = tmp_dmin; d <= tmp_dmax; d++)
		{
		  if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->ibeginsc[y];
		  else                       sc = alpha[y][cur][d]     + cm->ibeginsc[y];
		  if (sc > alpha[0][cur][d]) {
		    alpha[0][cur][d] = sc;
		    bestr[d]         = y;
		  }
		  if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
		  if (Scorify(alpha[0][cur][d]) > best_neg_score) best_neg_score = Scorify(alpha[0][cur][d]);
		}
	    }
	}
      }
	      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = hdmin[0][jp_v]; (d <= hdmax[0][jp_v] && d <= gamma_j) && d <= W; d++) 
	{
	  i       = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  assert(i > 0);
	  /*printf("v: %d gamma_i: %d d: %d\n", v, gamma_i, d);*/
	  /*printf("alpha[0][j:%3d][d:%3d]: %f\n", j, d, Scorify(alpha[0][cur][d]));*/
	  sc = gamma[gamma_i-1] + Scorify(alpha[0][cur][d]) + cm->sc_boost;
	  /* sc_boost is experimental technique for finding hits < 0 bits. 
	   * value is 0.0 if technique not used. */
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = Scorify(alpha[0][cur][d]); 
	      saver[gamma_j]  = bestr[d];
	    }
	}
    } /* end loop over end positions j */
  
  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   * be careful about only freeing our impossible deck, imp_row, once
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) 
	  if(alpha[v][j] != imp_row) free(alpha[v][j]); 
	free(alpha[v]);
      } else {
	if(alpha[v][0] != imp_row) free(alpha[v][0]);
	if(alpha[v][1] != imp_row) free(alpha[v][1]); 
	free(alpha[v]);
      }
    }
  free(imp_row);
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  j     = j0;
  while (j >= i0) 
    {
      gamma_j = j-i0+1;
      if (gback[gamma_j] == -1) /* no hit */
	j--; 
      else                /* a hit, a palpable hit */
	{
	  if(savesc[gamma_j] > best_score) 
	    best_score = savesc[gamma_j];
	  if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
	    report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	  j = gback[gamma_j]-1;
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  return best_score;

 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* never reached */
}

#if 0
/******************************************************************/
/* The below functions were written during debugging, and print
   out either the shadow or alpha matrix.  They are kept
   here just in case they're needed again.  Note : the functions
   that print out the entire matrix are really only useful
   when the BE_PARANOID flag is set, meaning that decks are
   never freed until the end.
*/
/*================================================================*/

/* Debugging functions that print info to STDOUT */
static void debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
					int **hdmin, int **hdmax);
static void debug_print_shadow_banded_jd(void ***shadow, CM_t *cm, int L, int *jmin, int *jmax, 
					 int **hdmin, int **hdmax);
static void debug_print_shadow_banded_deck_jd(int v, void ***shadow, CM_t *cm, int L, int *jmin, int *jmax,
					      int **hdmin, int **hdmax);

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
#endif

#if 0
/* Here are the non-memory efficient functions, kept around for reference */
/* The alignment engine (not memory efficient) */
static float inside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, 
			 int r, int z, int i0, int j0, 
			 int do_full,
			 float ***alpha, float ****ret_alpha, 
			 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			 void ****ret_shadow, 
			 int allow_begin, int *ret_b, float *ret_bsc,
			 int *jmin, int *jmax,
			 int **hdmin, int **hdmax,
			 int *safe_hdmin, int *safe_hdmax);

/* The traceback routine (not memory efficient) */
static float insideT_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			  int r, int z, int i0, int j0, int allow_begin,
			  int *jmin, int *jax, 
			  int **hdmin, int **hdmax,
			  int *safe_hdmin, int *safe_hdmax);


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
 *           sq        - the sequence [1..L]   
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
inside_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
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

  if(dsq == NULL)
    esl_fatal("ERROR, dsq is NULL.");

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
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
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
    ESL_ALLOC(shadow, sizeof(void **) * cm->M);
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
		/*else esl_fatal("cell in alpha matrix was not filled in due to bands.\n");*/
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
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);

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
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		
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
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		
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
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
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
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
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
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
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
insideT_b_jd(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
	     int r, int z, int i0, int j0, 
	     int allow_begin, int *jmin, int *jmax,
	     int **hdmin, int **hdmax,
	     int *safe_hdmin, int *safe_hdmax)
{
  if(dsq == NULL)
    esl_fatal("ERROR(), dsq is NULL.\n");

  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;                /* stack that tracks bifurc parent of a right start */
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

  pda = esl_stack_ICreate();
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
      esl_stack_IPush(pda, j);	/* remember the end j    */
      esl_stack_IPush(pda, k);	/* remember the subseq length k */
      esl_stack_IPush(pda, tr->n-1);	/* remember the trace index of the parent B state */
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
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      /* CYK Full ME Bands used 17 */
      /* original line : esl_stack_IPop(pda, &d); */
      esl_stack_IPop(pda, &dp);
      /* CYK Full ME Bands used 18 */
      /* line below added */

      /* Deal with end local states */
      if(cm->sttype[v] != EL_st)
	d = dp + safe_hdmin[y];
      else
      d = dp;

      esl_stack_IPop(pda, &j);
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
      default:    esl_fatal("'Inconceivable!'\n'You keep using that word...'");
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
  esl_stack_Destroy(pda);  /* it should be empty; we could check; naaah. */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  return sc;
}
#endif


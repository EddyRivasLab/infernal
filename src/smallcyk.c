/* smallcyk.c
 * SRE, Wed Aug  2 08:42:49 2000 [St. Louis]
 * CVS $Id$
 * 
 * Divide and conquer CYK alignment.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include <stdio.h>

#include "squid.h"

#include "structs.h"
#include "nstack.h"
#include "funcs.h"

/*################################################################*/
/* Functions: deckpool_*()
 * Date:      SRE, Wed Aug  2 10:43:17 2000 [St. Louis]
 *
 * Purpose:   Implementation of a pushdown stack for storing decks
 *            of the inside or outside dynamic programming matrices, with the
 *            usual _create, _push, _pop, and _free API. 
 *            
 *            The deck pool allows us to efficiently reuse memory,
 *            so long as our DP algorithms step through the decks
 *            as their outermost loop.
 *            
 *            Works for either coordinate system (vjd or vji) 
 *            and subseq variants, because it's simply managing
 *            a deck as a float **.
 */
struct deckpool_s {
  float ***pool;
  int      n;
  int      nalloc;
  int      block;
};
static struct deckpool_s *
deckpool_create(void)
{
  struct deckpool_s *dpool;

  dpool = MallocOrDie(sizeof(struct deckpool_s));
  dpool->block  = 10;		/* configurable if you want */
  dpool->pool   = MallocOrDie(sizeof(float **) * dpool->block);
  dpool->nalloc = dpool->block;;
  dpool->n      = 0;
  return dpool;
}
static void 
deckpool_push(struct deckpool_s *dpool, float **deck)
{
  if (dpool->n == dpool->nalloc) {
    dpool->nalloc += dpool->block;
    dpool->pool = ReallocOrDie(dpool->pool, sizeof(float **) * dpool->nalloc);
  }
  dpool->pool[dpool->n] = deck;
  dpool->n++;
}
static int
deckpool_pop(struct deckpool_s *d, float ***ret_deck)
{
  if (d->n == 0) { *ret_deck = NULL; return 0;}
  d->n--;
  *ret_deck = d->pool[d->n];
  return 1;
}
static void
deckpool_free(struct deckpool_s *d)
{
  free(d->pool);
  free(d);
}
/*================================================================*/


/*################################################################*/
/* Functions: *_vjd_*
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
 */
static float **
alloc_vjd_deck(int L, int i, int j)
{
  float **a;
  int     jp;
  
  a = MallocOrDie(sizeof(float *) * (L+1)); /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0;   jp < i-1;    jp++) a[jp]     = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp]     = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) a[jp+i-1] = MallocOrDie(sizeof(float) * (jp+1));

  return a;
}
static void
free_vjd_deck(float **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp]);
  free(a);
}
static void
free_vjd_matrix(float ***a, int vroot, int vend, int i, int j)
{
  int v;
  for (v = vroot; v <= vend; v++)
    if (a[v] != NULL)		/* protect against double free's of reused decks (ends) */
      { free_vjd_deck(a[v], i, j); a[v] = NULL; }
  free(a);
}
static char **
alloc_vjd_yshadow_deck(int L, int i, int j)
{
  char **a;
  int    jp;
  
  a = MallocOrDie(sizeof(char *) * (L+1)); /* always alloc 0..L rows, same as alloc_deck */
  for (jp = 0;   jp < i-1;    jp++) a[jp] = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp] = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) a[jp+i-1] = MallocOrDie(sizeof(char) * (jp+1));
  return a;
}
static void
free_vjd_yshadow_deck(char **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp]);
  free(a);
}
static int **
alloc_vjd_kshadow_deck(int L, int i, int j)
{
  int **a;
  int   jp;
  
  a = MallocOrDie(sizeof(int *) * (L+1)); /* always alloc 0..L rows, same as alloc_deck */
  for (jp = 0;   jp <  i-1;   jp++) a[jp] = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) a[jp+i-1] = MallocOrDie(sizeof(int) * (jp+1));
  for (jp = j+1; jp <= L;     jp++) a[jp] = NULL;
  return a;
}
static void
free_vjd_kshadow_deck(int **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp]);
  free(a);
}
static void
free_vjd_shadow_matrix(void ***shadow, CM_t *cm, int vroot, int vend, int i, int j)
{
  int v;
  for (v = vroot; v <= vend; v++)
    if (shadow[v] != NULL) {
      if (cm->sttype[v] == B_st) free_vjd_kshadow_deck((int **)  shadow[v], i, j);
      else                       free_vjd_yshadow_deck((char **) shadow[v], i, j);
    }
  free(shadow);
}
/*================================================================*/


/*################################################################*/
/* Functions: *_vji_*
 * Date:     SRE, Sat Aug 12 16:44:55 2000 [Titusville]
 *
 * Purpose:  Allocation and freeing of 3D matrices and 2D decks
 *           in the vji coordinate system. Since these are used
 *           by the segment_alignment_engine() special case, they work only
 *           on a defined segment in the 3D matrix: they need
 *           two triplets (vroot, i0, j0), (vend, i1, j1) 
 *           defining the known optimal endpoints of a segment from
 *           an S state to a B state.
 *
 *           Further, by our definition, we know there's no B states
 *           in between, so the shadow matrix doesn't need any
 *           special casting tricks the way the more generally
 *           used vjd system does.
 */
static float **                 /* allocation of a score deck. */
alloc_vji_deck(int i0, int j0, int i1, int j1)
{
  float **a;
  int     jp;
  a = MallocOrDie(sizeof(float *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    a[jp] = MallocOrDie(sizeof(float)*(i1-i0+1));
  return a;
}
static void			/* free'ing a score deck */
free_vji_deck(float **a, int j0, int j1)
{
  int jp;
  for (jp = 0; jp <= j0-j1; jp++) 
    if (a[jp] != NULL) free(a[jp]);
  free(a);
}
static void
free_vji_matrix(float ***a, int vroot, int j0, int vend, int j1)
{
  int v;
  for (v = vroot; v <= vend; v++) 
    if (a[v] != NULL) free_vji_deck(a[v], j0, j1);
  free(a);
}
static char **		        /* allocation of a traceback ptr (shadow matrix) deck */
alloc_vji_shadow_deck(int i0, int j0, int i1, int j1)
{
  char **a;
  int     jp;
  a = MallocOrDie(sizeof(char *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    a[jp] = MallocOrDie(sizeof(char)*(i1-i0+1));
  return a;
}
static void	                /* free'ing a shadow deck */
free_vji_shadow_deck(char **a, int j0, int j1)
{
  int jp;
  for (jp = 0; jp <= j0-j1; jp++) 
    if (a[jp] != NULL) free(a[jp]);
  free(a);
}
static void
free_vji_shadow_matrix(char ***a, int vroot, int j0, int vend, int j1)
{
  int v;
  for (v = vroot; v <= vend; v++) 
    if (a[v] != NULL) free_vji_shadow_deck(a[v], j0, j1);
  free(a);
}
/*================================================================*/



static float splitting_engine(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
			      int r, int vend, int i0, int j0);
static float cyk_inside_engine(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, 
			       int do_full, float ***alpha, float ****ret_alpha, 
			       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			       void ****ret_shadow);
static void  cyk_outside_engine(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0,
				int do_full, float ***beta, float ****ret_beta,
				struct deckpool_s *dpool, struct deckpool_s **ret_dpool);
static void  cyk_vjd_trace(CM_t *cm, void ***shadow, Parsetree_t *tr, int r, int i0, int j0);
static void  segment_alignment_engine(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
				      int r, int i0, int j0,
				      int s, int i1, int j1);

float
CYKInside(CM_t *cm, char *dsq, int L, Parsetree_t **ret_tr)
{
  Parsetree_t *tr;
  float        sc;
  void      ***shadow;

  /* Create a parse tree structure.
   * The traceback machinery expects to build on a start state already
   * in the parsetree, so initialize by adding the root state.
   */
  tr = CreateParsetree();
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Call the Inside() engine; the fill stage of the DP algorithm,
   * for the whole matrix.
   */
  sc = cyk_inside_engine(cm, dsq, L, 
			 0, cm->M-1, /* the whole model, from state 0 to M-1      */
			 1, L,       /* the whole sequence, from 1 to L           */
			 FALSE,      /* in memory saving mode                     */
			 NULL, NULL, /* create and destroy your own alpha matrix  */
			 NULL, NULL, /* create and destory your own deck pool     */
			 &shadow);   /* give me a shadow matrix to trace          */

  /* Trace back and recover a parse tree, to get the alignment
   * of the model to the sequence.
   */
  cyk_vjd_trace(cm, shadow, tr, 0, 1, L);
  
  /* Free memory and return
   */
  free_vjd_shadow_matrix(shadow, cm, 0, cm->M-1, 1, L);
  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  return sc;
}

float
CYKDivideAndConquer(CM_t *cm, char *dsq, int L, Parsetree_t **ret_tr)
{
  Parsetree_t *tr;
  float        sc;

  /* Create a parse tree structure.
   * The traceback machinery expects to build on a start state already
   * in the parsetree, so initialize by adding the root state.
   */
  tr = CreateParsetree();
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Start the divide and conquer recursion: call the splitting
   * engine on the whole DP cube.
   */
  sc = splitting_engine(cm, dsq, L, tr, 0, cm->M-1, 1, L);

  /* Free memory and return
   */
  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  return sc;
}


void
CYKDemands(CM_t *cm, int L)
{
  float Mb_per_deck;		/* megabytes per deck */
  int   bif_decks;		/* bifurcation decks  */
  int   nends;			/* end decks (only need 1, even for multiple E's */
  int   maxdecks;		/* maximum # of decks needed by CYKInside() */
  float smallmemory;		/* how much memory small version of CYKInside() needs */
  float bigmemory;		/* how much memory a full CYKInside() would take */
  float dpcells;		/* # of dp cells */
  float bifcalcs;		/* # of inner loops executed for bifurcation calculations */
  float dpcalcs;		/* # of inner loops executed for non-bif calculations */
  int   j;

  Mb_per_deck = (float)(L+2)*(float)(L+1)*2. / 1000000.;
  bif_decks   = CMCountStatetype(cm, B_st);
  nends       = CMCountStatetype(cm, E_st);
  maxdecks    = CYKDeckCount(cm);
  smallmemory = (float) maxdecks * Mb_per_deck;
  bigmemory   = (float) (cm->M - nends +1) * Mb_per_deck;
  dpcells     = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - nends +1);
  
  bifcalcs = 0.;
  for (j = 0; j <= L; j++)
    bifcalcs += (float)(j+1)*(float)(j+2)/2.;
  bifcalcs *= (float) bif_decks;
  
  dpcalcs = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - bif_decks - nends +1);

  printf("CYK cpu/memory demand estimates:\n");
  printf("Mb per cyk deck:                %.2f\n", Mb_per_deck);
  printf("# of decks (M):                 %d\n",   cm->M);
  printf("# of decks needed in small CYK: %d\n",   maxdecks);
  printf("RAM needed for full CYK, Mb:    %.2f\n", bigmemory);
  printf("RAM needed for small CYK, Mb:   %.2f\n", smallmemory);
  printf("# of dp cells, total:           %.0f\n", dpcells);
  printf("# of non-bifurc dp cells:       %.0f\n", dpcalcs);
  printf("# of bifurcations:              %d\n",   bif_decks);
  printf("# of bifurc dp inner loop calcs:%.0f\n", bifcalcs);
  printf("# of dp inner loops:            %.0f\n", dpcalcs+bifcalcs);
}

int
CYKDeckCount(CM_t *cm)
{
  Nstack_t *pda;		/* pushdown stack simulating the deck pool */
  int       v,y,z;		/* state indices */
  int       nends;
  int       ndecks;
  int      *touch;		/* keeps track of how many higher decks still need this deck */

  /* Initializations, mirroring key parts of CYKInside()
   */
  ndecks = 1;			/* the end deck, which we always need. */
  nends  = CMCountStatetype(cm, E_st);
  pda    = CreateNstack();
  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    touch[v] = cm->pnum[v];

  for (v = cm->M-1; v >= 0; v--)
    {
      if (cm->sttype[v] != E_st) {
	if (! PopNstack(pda, &y)) ndecks++; /* simulated allocation of a new deck */
      }
      
      if (cm->sttype[v] == B_st) { /* release both S children of a bifurc */
	y = cm->cfirst[v];
	z = cm->cnum[v];
	PushNstack(pda, y);
	PushNstack(pda, z);
      } else {
	for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	  {
	    touch[y]--;
	    if (touch[y] == 0) 
	      {
		if (cm->sttype[y] == E_st) { 
		  nends--; 
		  if (nends == 0) { PushNstack(pda, cm->M-1); }
		} else 
		  PushNstack(pda, y);
	      }
	  }
      }
    }
  FreeNstack(pda);
  free(touch);
  printf("Maximum # of decks = %d\n", ndecks);
  return ndecks;
}







/* Function: cyk_inside_engine()
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phase of a CYK alignment algorithm, on a 
 *           subsequence from i0..j0, using a subtree of a model
 *           anchored at a start state vroot, and ending at an end
 *           state vend. (It is a feature of the model layout in
 *           a CM structure that all subtrees are contiguous in the
 *           model.)
 *           
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
 *
 * Returns:  
 */
static float 
cyk_inside_engine(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
		  float ***alpha, float ****ret_alpha, 
		  struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		  void ****ret_shadow)
{
  float  **end;                 /* we re-use the end deck. */
  int      nends;               /* counter that tracks when we can release end deck to the pool */
  int     *touch;               /* keeps track of how many higher decks still need this deck */
  int      v,y,z;		/* indices for states  */
  int      j,d,i,k;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  int      yoffset;		/* y=base+offset -- counter in child states that v can transit to */
  int      W;			/* subsequence length */
  int      jp;			/* j': relative position in the subsequence  */
  void  ***shadow;              /* shadow matrix for tracebacks */
  int    **kshad;               /* a shadow deck for bifurcations */
  char   **yshad;               /* a shadow deck for every other kind of state */

  /* Allocations and initializations
   */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  W     = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, vroot, E_st);
  for (jp = 0; jp <= W; jp++) {
    j = i0+jp-1;		/* e.g. j runs from 0..L on whole seq */
    end[j][0] = 0.;
    for (d = 1; d <= jp; d++) end[j][d] = IMPOSSIBLE;
  }
				/* if caller didn't give us a matrix, make one */
  if (alpha == NULL) {
    alpha = MallocOrDie(sizeof(float **) * cm->M);
    for (v = 0; v < cm->M; v++) alpha[v] = NULL;
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
	alpha[v] = alloc_vjd_deck(L, i0, j0);

      if (ret_shadow != NULL) {
	if (cm->sttype[v] == B_st) {
	  kshad     = alloc_vjd_kshadow_deck(L, i0, j0); 
	  shadow[v] = (void **) kshad;
	} else {
	  yshad     = alloc_vjd_yshadow_deck(L, i0, j0); 
	  shadow[v] = (void **) yshad;
	}
      }

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = 0; d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d]  = alpha[y][j][d] + cm->tsc[v][0];
		if (ret_shadow != NULL) yshad[j][d]  = 0;
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc; 
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
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
		if (ret_shadow != NULL) kshad[j][d] = 0;
		for (k = 1; k <= d; k++)
		  if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) kshad[j][d] = k;
		  }
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
		alpha[v][j][d] = alpha[y][j-1][d-2] + cm->tsc[v][0];
		if (ret_shadow != NULL) yshad[j][d] = 0;
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
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
		alpha[v][j][d] = alpha[y][j][d-1] + cm->tsc[v][0];
		if (ret_shadow != NULL) yshad[j][d] = 0;
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
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
		alpha[v][j][d] = alpha[y][j-1][d-1] + cm->tsc[v][0];
		if (ret_shadow != NULL) yshad[j][d] = 0;
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		if (dsq[j] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}				/* finished calculating deck v. */
      
      
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
  } /* end loop over v */

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  sc = alpha[vroot][j0][W];

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

/* Function: cyk_vjd_trace()
 * Date:     SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]
 *
 * Purpose:  Traceback in a vjd shadow matrix produced by the
 *           cyk_inside_engine. Append the trace to a given
 *           traceback, which has an S_state at tr->n-1.
 */
static void
cyk_vjd_trace(CM_t *cm, void ***shadow, Parsetree_t *tr, int r, int i0, int j0)
{
  Nstack_t *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			/* reco */
  int       y, yoffset;
  int       bifparent;

  pda = CreateNstack();
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;
  while (1) {
    if (cm->sttype[v] == B_st) {
      k = ((int **) shadow[v])[j][d];   /* k = len of right fragment */

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
    } else if (cm->sttype[v] == E_st) {
      /* We don't trace back from an E. Instead, we're done with the
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
    } else {
      yoffset = ((char **)shadow[v])[j][d];
      y = cm->cfirst[v] + yoffset;

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

      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    }
  }
  FreeNstack(pda);  /* it should be empty; we could check; naaah. */
}



/* Function: cyk_outside_engine()
 * Date:     SRE, Tue Aug  8 10:42:52 2000 [St. Louis]
 *
 * Purpose:  Run the outside version of a CYK alignment algorithm,
 *           on a subsequence i0..j0 of a digitized sequence dsq [1..L],
 *           using a linear segment of a model anchored at a start state vroot
 *           (possibly the absolute root, 0) and ending at an end
 *           state or bifurcation state vend. There must be no
 *           start, end, or bifurcation states in the path other than 
 *           these termini: this is not a full Outside implementation,
 *           it is only the bit that's necessary in the divide
 *           and conquer alignment algorithm.
 *           
 *           Much of the behavior in calling conventions, etc., is
 *           analogous to the cyk_inside_engine(); see its preface
 *           for more info.
 *           
 *           At the end of the routine, the bottom deck (vend) is valid.
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           vroot     - first state of linear model segment (always a start)
 *           vend      - last state of linear model segment 
 *           i0        - first position in subseq to align (1, for whole seq)
 *           j0        - last position in subseq to align (L, for whole seq)
 *           do_full   - if TRUE, we save all the decks in beta, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the lowermost deck (vend) is valid
 *                       at the end.
 *           beta      - if non-NULL, this is an existing matrix, with NULL
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
 */
static void
cyk_outside_engine(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0,
		   int do_full, float ***beta, float ****ret_beta,
		   struct deckpool_s *dpool, struct deckpool_s **ret_dpool)
{
  int      v,y;			/* indices for states */
  int      j,d,i;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */
  int      W;			/* subsequence length */
  int      jp;			/* j': relative position in the subsequence, 0..W */
  int      voffset;		/* index of v in t_v(y) transition scores */

  /* Allocations and initializations
   */
  W = j0-i0+1;		/* the length of the subsequence: used in many loops */

  			/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

                        /* if caller didn't give us a matrix, make one */
  if (beta == NULL) {
    beta = MallocOrDie(sizeof(float **) * cm->M);
    for (v = 0; v < cm->M; v++) beta[v] = NULL;
  }

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0;      v < vroot; v++) touch[v] = 0;
  for (v = vend+1; v < cm->M; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) {
    if (cm->sttype[v] == B_st) touch[v] = 2; /* well, we'll never use this, but set it anyway. */
    else                       touch[v] = cm->cnum[v];
  }
				
  /* Initialize the root deck. This probably isn't the most efficient way to do it.
   */
  if (! deckpool_pop(dpool, &(beta[vroot])))
    beta[vroot] = alloc_vjd_deck(L, i0, j0);
  for (jp = 0; jp <= W; jp++) {
    j = i0-1+jp;
    for (d = 0; d <= jp; d++)
      beta[vroot][j][d] = IMPOSSIBLE;
  }
  beta[vroot][j0][W] = 0;		
  
  /* Main loop down through the decks
   */
  for (v = vroot+1; v <= vend; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_vjd_deck(L, i0, j0);

      /* main recursion:
       */
      for (jp = W; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  {
	    beta[v][j][d] = IMPOSSIBLE;
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
		
		if ((sc = beta[y][j+1][d+2] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;

	      case ML_st:
	      case IL_st: 
		if (d == jp) continue;	/* boundary condition (note when j=0, d=0*/

		if (dsq[i-1] < Alphabet_size) 
		  escore = cm->esc[y][(int) dsq[i-1]];
		else
		  escore = DegenerateSingletScore(cm->esc[y], dsq[i-1]);
		  
		if ((sc = beta[y][j][d+1] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		  
		if (dsq[j+1] < Alphabet_size) 
		  escore = cm->esc[y][(int) dsq[j+1]];
		else
		  escore = DegenerateSingletScore(cm->esc[y], dsq[j+1]);

		if ((sc = beta[y][j+1][d+1] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		if ((sc = beta[y][j][d] + cm->tsc[y][voffset]) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;

	      default: Die("bogus parent state %d\n", cm->sttype[y]);
	      }/* end switch over states*/
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;

	  } /* ends loop over d. We know all beta[v][j][d] in this row j*/

      }/* end loop over jp. We know the beta's for the whole deck.*/

      /* Finished deck v.
       * now look at its parents; if we're reusing memory (! do_full)
       * push the parents that we don't need any more into the pool.
       */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { deckpool_push(dpool, beta[y]); beta[y] = NULL; }
	}
      }

    } /* end loop over decks v. */

  /* If the caller doesn't want the matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = vroot; v <= vend; v++)
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
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
}


static float
splitting_engine(CM_t *cm, char *dsq, int L, Parsetree_t *tr, int r, int vend, int i0, int j0)
{
  float ***alpha;
  float ***beta;
  struct deckpool_s *pool;
  int      v,y,z;		/* state indices */
  int      jp;			/* j': relative position in subseq, 0..W */
  int      W;			/* length of subseq i0..j0 */
  float    sc;			/* tmp variable for a score */
  int      j,d,k;		/* sequence indices */
  float    best_sc;		/* optimal score at the optimal split point */
  int      best_k;		/* optimal k for the optimal split */
  int      best_d;		/* optimal d for the optimal split */
  int      best_j;		/* optimal j for the optimal split */
  int      tv;			/* remember the position of a bifurc in the trace. */

  /* 1. Traverse down from r, find first end or bifurc;
   *    set up our state quartet r,v,y,z
   */
  for (v = r; v <= vend; v++)
    {
      if (cm->sttype[v] == B_st) break;
      if (cm->sttype[v] == E_st) {
	void ***shadow;
	/* we're clear to an END! finish this off. This is how the recursion terminates.
         * Have the inside engine do everything and free its memory -- all we want it
         * to do is add to the trace.
	 */
	sc = cyk_inside_engine(cm, dsq, L, r, v, i0, j0, FALSE, NULL, NULL, NULL, NULL, &shadow);
	cyk_vjd_trace(cm, shadow, tr, r, i0, j0);
	free_vjd_shadow_matrix(shadow, cm, r, v, i0, j0);
	return sc; /* no split to be done. more code needed here eventually. */
      }
    }
  y = cm->cfirst[v];		/* left S  */
  z = cm->cnum[v];		/* right S */

  /* 2. Calculate alpha[y] deck and alpha[z] deck.
   */
  cyk_inside_engine(cm, dsq, L, y, z-1,  i0, j0, FALSE, NULL,  &alpha, NULL, &pool, NULL);
  cyk_inside_engine(cm, dsq, L, z, vend, i0, j0, FALSE, alpha, &alpha, pool, &pool, NULL);

  /* 3. Calculate beta[v] deck (stick it in alpha). Let the pool get free'd.
   */
  cyk_outside_engine(cm, dsq, L, r, v, i0, j0, TRUE, alpha, &beta, pool, NULL);

  /* 4. Find the optimal split.
   */
  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  for (jp = 0; jp <= W; jp++) 
    {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	for (k = 0; k <= d; k++)
	  if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k] + beta[v][j][d]) > best_sc) 
	    {
	      best_sc = sc;
	      best_k  = k;
	      best_j  = j;
	      best_d  = d;
	    }
    }

  /* now, the optimal split:
   * left fragment: i1 = j-d+1, j1 = j-k, vroot = y, vend = z-1
   * right frag:    i2 = j-k+1, j2 = j,   vroot = z, vend = caller's vend
   * Execute the following steps:
   *    1. do the segment_alignment_engine() on the segment we've just defined,
   *       which gives us the traceback down to the bifurc (inclusive)
   *    2. add a left start to the trace.
   *    3. recursively call the splitting engine on the left side
   *    4. add a right start to the bifurc
   */
  segment_alignment_engine(cm, dsq, L, tr, r, i0, j0, v, best_j-best_d+1, best_j);
  tv = tr->n-1;

  InsertTraceNode(tr, tv, TRACE_LEFT_CHILD, best_j-best_d+1, best_j-best_k, y);
  splitting_engine(cm, dsq, L, tr, y, z-1, best_j-best_d+1, best_j-best_k);
  InsertTraceNode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j, z);
  splitting_engine(cm, dsq, L, tr, z, vend, best_j-best_k+1, best_j);

  /* The two alpha matrices and the beta matrix
   * actually all point to the same memory, since no
   * decks in Inside and Outside needed to overlap. 
   * Free 'em all in one call.
   */
  free_vjd_matrix(alpha, r, vend, i0, j0);

  return best_sc;
}


/* Function: segment_alignment_engine()
 * Date:     SRE, Thu Aug 10 12:05:48 2000 [St. Louis]
 *
 * Purpose:  Given that the splitting_engine has determined
 *           that the optimal alignment is bounded by
 *           r,i0,j0 on the outside and s,i1,j1 on the inside;
 *           find the optimal alignment in between.
 *           
 *           r must be a start state.
 *           s must be a bifurc state. (In principle it could
 *             also be an end state, but in that case, we would
 *             have to know the exact position i1,j1 for the end,
 *             which is not normally known.)
 *             
 *           The states r+1..s-1 must be a continuous
 *           linear model segment -- no bifurcations, no ends.
 *
 * Args:     
 *
 * Returns:  
 *
 * Example:  
 */
static void
segment_alignment_engine(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
			 int r, int i0, int j0,
			 int s, int i1, int j1)
{
  float ***a;                   /* the score matrix -- memory is reused */
  char  ***shadow;              /* the shadow matrix -- traceback ptrs -- memory is kept */
  struct deckpool_s *dpool;     /* stack for reuse of score matrix decks */
  int     v,i,j;
  int     jp, ip;		/* j' and i' -- in the matrix coords */
  int    *touch;                /* keeps track of whether we can free a deck yet or not */
  int     y, yoffset;
  float   sc;			/* tmp variable holding a score */

  /* Allocations, initializations
   */
  dpool = deckpool_create();
  
  a = MallocOrDie(sizeof(float **) * cm->M);
  for (v = 0; v < cm->M; v++) a[v] = NULL;

  shadow = MallocOrDie(sizeof(char **) * cm->M);
  for (v = 0; v < cm->M; v++) shadow[v] = NULL; 

				/* init the last deck, s; we don't need a shadow for it. */
  a[s] = alloc_vji_deck(i0, j0, i1, j1);
  for (jp = 0; jp <= j0-j1; jp++) 
    for (ip = 0; ip <= i1-i0; ip++) 
      a[s][jp][ip] = IMPOSSIBLE;
  a[s][0][i1-i0] = 0.;

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0;   v < r;  v++) touch[v] = 0;
  for (v = r;   v <= s; v++) touch[v] = cm->pnum[v];
  for (v = s+1; v < cm->M; v++) touch[v] = 0;

  /* Main recursion
   */
  for (v = s-1; v >= r; v--)
    {
      /* Get a deck and a shadow deck.
       */
      if (! deckpool_pop(dpool, &(a[v]))) 
	a[v] = alloc_vji_deck(i0, j0, i1, j1);
      shadow[v] = alloc_vji_shadow_deck(i0,j0,i1,j1);

				/* reassert our definition of a model segment */
      if (cm->sttype[v] == E_st || cm->sttype[v] == B_st || (cm->sttype[v] == S_st && v > r))
	Die("you told me you wouldn't ever do that again.");
      
      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= j0-j1; jp++) 
	    for (ip = i1-i0; ip >= 0; ip--) {
	      y = cm->cfirst[v];
	      a[v][jp][ip]      = a[y][jp][ip] + cm->tsc[v][0];
	      shadow[v][jp][ip] = (char) 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { a[v][jp][ip] = sc; shadow[v][jp][ip] = (char) yoffset; }
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;
	    }
	} else if (cm->sttype[v] == MP_st) {
	  for (ip = i1-i0; ip >= 0; ip--) a[v][0][ip] = IMPOSSIBLE; /* boundary condition */

	  for (jp = 1; jp <= j0-j1; jp++) { 
	    j = jp+j1;
	    a[v][jp][i1-i0] = IMPOSSIBLE; /* boundary condition */
	    for (ip = i1-i0-1; ip >= 0; ip--) {
	      i = ip+i0;
	      y = cm->cfirst[v];
	      a[v][jp][ip] = a[y][jp-1][ip+1] + cm->tsc[v][0];
	      shadow[v][jp][ip] = 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		   { a[v][jp][ip] = sc; shadow[v][jp][ip] = (char) yoffset; }

	      if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		a[v][jp][ip] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
	      else
		a[v][jp][ip] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
	  
	  for (jp = 0; jp <= j0-j1; jp++) { 
	    a[v][jp][i1-i0] = IMPOSSIBLE; /* boundary condition */
	    for (ip = i1-i0-1; ip >= 0; ip--) {
	      i = ip+i0;
	      y = cm->cfirst[v];
	      a[v][jp][ip] = a[y][jp][ip+1] + cm->tsc[v][0];
	      shadow[v][jp][ip] = 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { a[v][jp][ip] = sc; shadow[v][jp][ip] = (char) yoffset; }
	      
	      if (dsq[i] < Alphabet_size)
		a[v][jp][ip] += cm->esc[v][(int) dsq[ip]];
	      else
		a[v][jp][ip] += DegenerateSingletScore(cm->esc[v], dsq[i]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
	  for (ip = i1-i0; ip >= 0; ip--) a[v][0][ip] = IMPOSSIBLE; /* boundary condition */

	  for (jp = 1; jp <= j0-j1; jp++) { 
	    j = jp+j1;
	    for (ip = i1-i0; ip >= 0; ip--) {
	      y = cm->cfirst[v];
	      a[v][jp][ip]      = a[y][jp-1][ip] + cm->tsc[v][0];
	      shadow[v][jp][ip] = 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { a[v][jp][ip] = sc; shadow[v][jp][ip] = (char) yoffset; }
	      
	      if (dsq[j] < Alphabet_size)
		a[v][jp][ip] += cm->esc[v][(int) dsq[j]];
	      else
		a[v][jp][ip] += DegenerateSingletScore(cm->esc[v], dsq[j]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} /* finished calculating deck v */

      /* Now, try to reuse memory under v
       */
      for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	{
	  touch[y]--;
	  if (touch[y] == 0) { 
	    deckpool_push(dpool, a[y]);
	    a[y] = NULL;
	  }
	}
    } /* end loop over v; we now have a complete matrix */
	
  /* We're done with the score matrix, so blow it away, along with the
   * other machinery -- everything but the shadow matrix.
   * deck r is active in a[][][], and we have an unknown number of decks in the pool.
   */
  free_vji_deck(a[r], j0, j1);
  while (deckpool_pop(dpool, &(a[r]))) free_vji_deck(a[r], j0, j1);
  deckpool_free(dpool);
  free(a);
  free(touch);

  /* We've got a complete shadow matrix. Trace it back. We know
   * that the trace will begin with the start state r, at i0,j0
   * (e.g. jp=j0-j1, ip=0)
   */
  v  = r;
  j = j0;
  i = i0;
  while (v < s) {
    jp = j-j1;
    ip = i-i0;

    /* 1. figure out the next state (deck) in the shadow matrix.
     */ 
    yoffset = shadow[v][jp][ip];
    y = cm->cfirst[v] + yoffset;

    /* 2. figure out the i,j for state y, which is dependent 
     *    on what v emits (if anything)
     */
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
    
    /* 3. Attach y,i,j to the trace. This new node always attaches
     *    to the end of the growing trace -- e.g. trace node
     *    tr->n-1.
     */
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
    v = y;
  }
  
  /* We're done. Our traceback has just ended. We have just attached
   * the bifurcation at s,i1,j1; it is in the traceback at node tr->n-1.
   * Free the shadow matrix and return.
   */
  for (v = r; v < s; v++) free_vji_shadow_deck(shadow[v], j0, j1);
  free(shadow);
}


/* Function: CYKOutside()
 * Date:     SRE, Mon Aug  7 07:45:37 2000 [St. Louis]
 */
void
CYKOutside(CM_t *cm, char *dsq, int L, float ***alpha)
{
  float ***beta;		/* the scoring cube [v=0..M-1][j=0..L][d=0..j]*/
  int      v,y,z;		/* indices for states */
  int      j,d,i,k;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  struct deckpool_s *dpool;     /* a pool of decks for beta that we can reuse */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */

  /* Allocations and initializations
   */
  beta = MallocOrDie(sizeof(float **) * cm->M);
  for (v = 0; v < cm->M; v++) beta[v] = NULL;

  dpool = deckpool_create();

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; j++)
      beta[0][j][d] = IMPOSSIBLE; /* can prob speed this initialization up */
  beta[0][L][L] = 0;		
  
  /* Main loop down through the decks
   */
  for (v = 2; v < cm->M; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_vjd_deck(L, 1, L);

      /* main recursion:
       */
      for (j = L; j >= 0; j--)
	for (d = j; d >= 0; d--) 
	  {
	    if (cm->stid[v] == BEGL_S) 
	      {
		y = cm->plast[v];	/* the parent bifurcation    */
		z = cm->cnum[y];	/* the other (right) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j][0]; /* init on k=0 */
		for (k = 1; k <= L-j; k++)
		  if ((sc = beta[y][j+k][d+k] + alpha[z][j+k][k]) > beta[v][j][d])
		    beta[v][j][d] = sc;
	      }
	    else if (cm->stid[v] == BEGR_S) 
	      {
		y = cm->plast[v];	        /* the parent bifurcation    */
		z = cm->cfirst[y];	/* the other (left) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j-d][0];	/* init on k=0 */
		for (k = 1; k <= j-d; k++) 
		  if ((sc = beta[y][j][d+k] + alpha[z][j-d][k]) > beta[v][j][d])
		    beta[v][j][d] = sc;
	      }
	    else
	      {
		alpha[v][j][d] = IMPOSSIBLE;
		i = j-d+1;
		for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
		  switch(cm->sttype[j]) {
		  case MP_st: 
		    if (d == j || d == j-1) continue; /* boundary condition */

		    if (dsq[i-1] < Alphabet_size && dsq[j+1] > Alphabet_size)
		      escore = cm->esc[y][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
		    else
		      escore = DegeneratePairScore(cm->esc[y], dsq[i-1], dsq[j+1]);

		    if ((sc = beta[y][j+1][d+2] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;

		  case ML_st:
		  case IL_st: 
		    if (d == j) continue;	/* boundary condition (note when j=0, d=0*/

		    if (dsq[i-1] < Alphabet_size) 
		      escore = cm->esc[y][(int) dsq[i-1]];
		    else
		      escore = DegenerateSingletScore(cm->esc[y], dsq[i-1]);
		  
		    if ((sc = beta[y][j][d+1] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;
		  
		  case MR_st:
		  case IR_st:
		    if (d == j || j == L) continue;
		  
		    if (dsq[j+1] < Alphabet_size) 
		      escore = cm->esc[y][(int) dsq[j+1]];
		    else
		      escore = DegenerateSingletScore(cm->esc[y], dsq[j+1]);

		    if ((sc = beta[y][j+1][d+1] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;
		  
		  case B_st:
		  case E_st:
		  case D_st:
		    if ((sc = beta[y][j][d] + cm->tsc[y][v]) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;

		  default: Die("bogus parent state %d\n", cm->sttype[y]);
		  }/* end switch over states*/
		}
	      }/*ends our handling of beta[v][j][d] */
	    if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;
	  }

      /* Finished deck v.
       * now worry about reuse of memory in beta:
       */
      for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--)
	{
	  touch[y]--;
	  if (touch[y] == 0) {
	    deckpool_push(dpool, beta[y]);
	    beta[y] = NULL;
	  }
	}
    } /* end loop over decks v. */

  free(touch);
  /*dpool*/
  /*beta*/
}

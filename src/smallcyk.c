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

#include "config.h"
#include "structs.h"
#include "nstack.h"
#include "funcs.h"

struct deckpool_s {
  float ***pool;
  int      n;
  int      nalloc;
  int      block;
};

/* The dividers and conquerors.
 */
static float generic_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
			      int r, int vend, int i0, int j0);
static float wedge_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
			    int r, int z, int i0, int j0);
static void  v_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
			int r, int z, int i0, int i1, int j1, int j0);

/* The alignment engines. 
 */
static float inside(CM_t *cm, char *dsq, int L, 
		    int r, int z, int i0, int j0, 
		    int do_full,
		    float ***alpha, float ****ret_alpha, 
		    struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		    void ****ret_shadow);
static void  outside(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0,
		     int do_full, float ***beta, float ****ret_beta,
		     struct deckpool_s *dpool, struct deckpool_s **ret_dpool);
static float vinside(CM_t *cm, char *dsq, int L, 
		     int r, int z, int i0, int i1, int j1, int j0,
		     int do_full, float ***a, float ****ret_a,
		     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		     char ****ret_shadow);
static void  voutside(CM_t *cm, char *dsq, int L, 
		      int r, int z, int i0, int i1, int j1, int j0,
		      int do_full, float ***beta, float ****ret_beta,
		      struct deckpool_s *dpool, struct deckpool_s **ret_dpool);

/* The traceback routines.
 */
static float insideT(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
		     int r, int z, int i0, int j0);
static float vinsideT(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
		      int r, int z, int i0, int i1, int j1, int j0);

/* The size calculators.
 */
static float insideT_size(CM_t *cm, int L, int r, int z, int i0, int j0);
static float vinsideT_size(CM_t *cm, int r, int z, int i0, int i1, int j1, int j0);
static int   cyk_deck_count(CM_t *cm, int r, int z);

/* The memory management routines.
 */
static struct  deckpool_s *deckpool_create(void);
static void    deckpool_push(struct deckpool_s *dpool, float **deck);
static int     deckpool_pop(struct deckpool_s *d, float ***ret_deck);
static void    deckpool_free(struct deckpool_s *d);
static float **alloc_vjd_deck(int L, int i, int j);
static float   size_vjd_deck(int L, int i, int j);
static void    free_vjd_deck(float **a, int i, int j);
static void    free_vjd_matrix(float ***a, int r, int z, int i, int j);
static char  **alloc_vjd_yshadow_deck(int L, int i, int j);
static float   size_vjd_yshadow_deck(int L, int i, int j);
static void    free_vjd_yshadow_deck(char **a, int i, int j);
static int   **alloc_vjd_kshadow_deck(int L, int i, int j);
static float   size_vjd_kshadow_deck(int L, int i, int j);
static void    free_vjd_kshadow_deck(int **a, int i, int j);
static void    free_vjd_shadow_matrix(void ***shadow, CM_t *cm, int r, int z, int i, int j);
static float **alloc_vji_deck(int i0, int i1, int j1, int j0);
static float   size_vji_deck(int i0, int i1, int j1, int j0);
static void    free_vji_deck(float **a, int j1, int j0);
static void    free_vji_matrix(float ***a, int r, int z, int j1, int j0);
static char  **alloc_vji_shadow_deck(int i0, int i1, int j1, int j0);
static float   size_vji_shadow_deck(int i0, int i1, int j1, int j0);
static void    free_vji_shadow_deck(char **a, int j1, int j0);
static void    free_vji_shadow_matrix(char ***a, int r, int z, int j1, int j0);

#define BE_EFFICIENT  0		/* setting for do_full: small memory mode */
#define BE_PARANOID   1		/* setting for do_full: keep whole matrix */

/*################################################################
 * The API.
 * 
 * CYKDivideAndConquer() - the main routine. Align a model to a sequence.
 *################################################################
 */  


/* Function: CYKDivideAndConquer()
 * Date:     SRE, Sun Jun  3 19:32:14 2001 [St. Louis]
 *
 * Purpose:  Align a CM to a sequence using the divide and conquer
 *           algorithm. Return the score (in bits) and a traceback
 *           structure.
 *
 * Args:     cm     - the covariance model
 *           dsq    - the sequence, 1..L
 *           L      - length of sequence
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *
 * Returns:  score of the alignment in bits.
 */
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

  /* Start the divide and conquer recursion: call the generic_splitter()
   * on the whole DP cube.
   */
  sc = generic_splitter(cm, dsq, L, tr, 0, cm->M-1, 1, L);

  /* Free memory and return
   */
  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  return sc;
}

/* Function: CYKInside()
 * Date:     SRE, Sun Jun  3 19:48:33 2001 [St. Louis]
 *
 * Purpose:  Wrapper for the insideT() routine - solve
 *           a full alignment problem, return the traceback
 *           and the score.
 *           
 * Args:     cm     - the covariance model
 *           dsq    - the sequence, 1..L
 *           L      - length of sequence
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *
 * Returns:  score of the alignment in bits.
 */
float
CYKInside(CM_t *cm, char *dsq, int L, Parsetree_t **ret_tr)
{
  Parsetree_t *tr;
  float        sc;

  /* Create a parse tree structure.
   * The traceback machinery expects to build on a start state already
   * in the parsetree, so initialize by adding the root state.
   */
  tr = CreateParsetree();
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, 1, L, 0); /* init: attach the root S */

  /* Solve the whole thing with one call to insideT.
   */
  sc = insideT(cm, dsq, L, tr, 0, cm->M-1, 1, L);

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  return sc;
}

/* Function: CYKDemands()
 * Date:     SRE, Sun Jun  3 20:00:54 2001 [St. Louis]
 *
 * Purpose:  Print out information on the computational
 *           complexity of an alignment problem for divide
 *           and conquer versus the full CYK.
 *
 * Args:     cm - the model
 *           L  - length of sequence.
 */
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

  Mb_per_deck = size_vjd_deck(L, 1, L);
  bif_decks   = CMCountStatetype(cm, B_st);
  nends       = CMCountStatetype(cm, E_st);
  maxdecks    = cyk_deck_count(cm, 0, cm->M-1);
  smallmemory = (float) maxdecks * Mb_per_deck;
  bigmemory   = (float) (cm->M - nends +1) * Mb_per_deck;
  dpcells     = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - nends +1);
  
  bifcalcs = 0.;
  for (j = 0; j <= L; j++)
    bifcalcs += (float)(j+1)*(float)(j+2)/2.;
  bifcalcs *= (float) bif_decks;
  
  dpcalcs = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - bif_decks - nends +1);

  printf("CYK cpu/memory demand estimates:\n");
  printf("Mb per cyk deck:                 %.4f\n", Mb_per_deck);
  printf("# of decks (M):                  %d\n",   cm->M);
  printf("# of decks needed in small CYK:  %d\n",   maxdecks);
  printf("RAM needed for full CYK, Mb:     %.2f\n", bigmemory);
  printf("RAM needed for small CYK, Mb:    %.2f\n", smallmemory);
  printf("# of dp cells, total:            %.3g\n", dpcells);
  printf("# of non-bifurc dp cells:        %.3g\n", dpcalcs);
  printf("# of bifurcations:               %d\n",   bif_decks);
  printf("# of bifurc dp inner loop calcs: %.3g\n", bifcalcs);
  printf("# of dp inner loops:             %.3g\n", dpcalcs+bifcalcs);
}


/*################################################################
 * The dividers and conquerors. 
 *################################################################*/  

/* Function: generic_splitter()
 * Date:     SRE, Sat May 12 15:08:38 2001 [CSHL]
 *
 * Purpose:  Solve a "generic problem": best parse of
 *           a possibly bifurcated subgraph cm^r_vend to
 *           a substring dsq[i0..j0]. r is always a start
 *           state (S_st) and vend is always an end 
 *           state (E_st).
 *           Given: a cm subgraph from r..vend
 *                  a subsequence from i0..j0
 *           Attaches the optimal trace T{r..vend}, exclusive of r
 *           and inclusive of vend, to tr.
 *           
 *           A full divide & conquer never terminates
 *           in generic_splitter; the recursion must
 *           terminate in v_splitter and wedge_splitter;
 *           so we don't test an end-of-recursion boundary.
 *           
 * Args:     cm -  model
 *           dsq - digitized sequence 1..L
 *           L   - length of dsq
 *           tr  - the traceback we're adding on to.
 *           r   - index of a start state (S_st) in the model       
 *           vend- index of an end state (E_st) in the model
 *           i0  - start in the sequence (1..L)
 *           j0  - end in the sequence (1..L)
 *
 * Returns:  score of the optimal parse of dsq(i0..j0) with cm^r_z 
 */
static float
generic_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
		 int r, int vend, int i0, int j0)
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

  /* 1. If the generic problem is small enough, solve it with inside^T,
   *    and append the trace to tr.
   */
  if (insideT_size(cm, L, r, vend, i0, j0) < RAMLIMIT) {
    SQD_DPRINTF1(("Solving a generic w/ insideT - G%d[%s]..%d[%s], %d..%d\n",
		  r, UniqueStatetype(cm->stid[r]),
		  vend, UniqueStatetype(cm->stid[vend]),
		  i0, j0));
    sc = insideT(cm, dsq, L, tr, r, vend, i0, j0);
    return sc;
  }

  /* 2. Traverse down from r, find first bifurc.
   *    The lowest a bifurc could be: B-S-E/S-IL-E = vend-5
   *                                   
   */
  for (v = r; v <= vend-5; v++)
    if (cm->sttype[v] == B_st) break; /* found the first bifurcation, now v */

  /* 3. If there was no bifurcation, this is a wedge problem; solve it
   *    with wedge_splitter. 
   */
  if (v > vend-5) {		/* no bifurc? it's a wedge problem  */
    if (cm->sttype[vend] != E_st) Die("inconceivable.");
    sc = wedge_splitter(cm, dsq, L, tr, r, vend, i0, j0);
    return sc;
  }

  /* Set up the state quartet r,v,y,z, for a divide and conquer
   * solution of the generic problem.
   */
  y = cm->cfirst[v];		/* left S  (v+1) */
  z = cm->cnum[v];		/* right S */

  /* Calculate alpha[y] deck and alpha[z] deck.
   */
  inside(cm, dsq, L, y, z-1,  i0, j0, BE_EFFICIENT, NULL,  &alpha, NULL, &pool, NULL);
  inside(cm, dsq, L, z, vend, i0, j0, BE_EFFICIENT, alpha, &alpha, pool, &pool, NULL);

  /* Calculate beta[v] deck (stick it in alpha). Let the pool get free'd.
   */
  outside(cm, dsq, L, r, v, i0, j0, BE_EFFICIENT, alpha, &beta, pool, NULL);

  /* Find the optimal split.
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

  /* Free now, before recursing.
   * The two alpha matrices and the beta matrix
   * actually all point to the same memory, since no
   * decks in Inside and Outside needed to overlap. 
   * Free 'em all in one call.
   */
  free_vjd_matrix(alpha, r, vend, i0, j0);

  /* now, the optimal split, into a V problem and two generic problems.
   * left fragment: i1 = j-d+1, j1 = j-k, vroot = y, vend = z-1
   * right frag:    i2 = j-k+1, j2 = j,   vroot = z, vend = caller's vend
   * 
   * The problems must be solved in a particular order, since we're
   * constructing the trace in a postorder traversal.
   */
  SQD_DPRINTF1(("Generic splitter:\n"));
  SQD_DPRINTF1(("   V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		v, UniqueStatetype(cm->stid[v]),
		i0, best_j-best_d+1, best_j, j0));
  SQD_DPRINTF1(("   generic: G%d[%s]..%d[%s], %d..%d\n", 
		y, UniqueStatetype(cm->stid[y]),
		z-1, UniqueStatetype(cm->stid[z-1]),
		best_j-best_d+1, best_j-best_k));
  SQD_DPRINTF1(("   generic: G%d[%s]..%d[%s], %d..%d\n", 
		z, UniqueStatetype(cm->stid[z]),
		vend, UniqueStatetype(cm->stid[vend]),
		best_j-best_k+1, best_j));

  v_splitter(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0);
  tv = tr->n-1;

  InsertTraceNode(tr, tv, TRACE_LEFT_CHILD, best_j-best_d+1, best_j-best_k, y);
  generic_splitter(cm, dsq, L, tr, y, z-1, best_j-best_d+1, best_j-best_k);
  InsertTraceNode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j, z);
  generic_splitter(cm, dsq, L, tr, z, vend, best_j-best_k+1, best_j);

  return best_sc;
}

/* Function: wedge_splitter()
 * Date:     SRE, Sun May 13 08:44:15 2001 [CSHL genome mtg]
 *
 * Purpose:  Solve a "wedge problem": best parse of an 
 *           unbifurcated subgraph cm^r..z to a substring
 *           dsq[i0..j0]. r may be a start state (when
 *           the wedge problem comes from being a special case
 *           of a generic problem) or a non-insert state
 *           (D, MP, ML, MR) (when the wedge comes from a
 *           previous wedge_splitter). z is always an end state.
 *           
 *           Attaches the optimal trace T(r..z), exclusive
 *           of r and inclusive of z, to the growing trace tr.
 *           
 *           Deal with a divide and conquer boundary condition:
 *           the next non-insert state after r is the end state z.
 *           All remaining sequence of i0..j0 that r doesn't emit
 *           must be dealt with by insert states.
 *
 * Args:     cm -  model
 *           dsq - digitized sequence 1..L
 *           L   - length of dsq
 *           tr  - the traceback we're adding on to.
 *           r   - index of the first state in the subgraph (S; or D,ML,MR,or MP)
 *           z   - index of an end state (E_st) in the model
 *           i0  - start in the sequence (1..L)
 *           j0  - end in the sequence (1..L)
 *
 * Returns:  The score of the best parse in bits.
 */
static float 
wedge_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int j0)
{
  float ***alpha;
  float ***beta;
  struct deckpool_s *pool;
  float sc;
  float best_sc;
  int   v,w,y;
  int   W;
  int   d, jp, j;
  int   best_v, best_d, best_j;
  int   midnode;
  
  /* 1. If the wedge problem is either a boundary condition,
   *    or small enough, solve it with inside^T and append
   *    the trace to tr. 
   *    It's formally possible that someone could set RAMLIMIT
   *    to something so small that even the boundary condition
   *    couldn't be done with inside^T - but that'd be a silly
   *    thing to do, so we ignore RAMLIMIT in that case.
   */
  if (cm->ndidx[z] == cm->ndidx[r] + 1 ||
      insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT) 
    {
      SQD_DPRINTF1(("Solving a wedge:   G%d[%s]..%d[%s], %d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		z, UniqueStatetype(cm->stid[z]),
		i0,j0));
      sc = insideT(cm, dsq, L, tr, r, z, i0, j0);
      return sc;
    }

  /* 2. Find our split set, w..y
   *    We choose the node in the middle.
   *    This can't be a BIF_nd (we're a wedge), or an END_nd (midnode
   *    can't be z) but it could be any other node including
   *    begin nodes (i.e. it might be that w==y).
   */
  midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r]) / 2);
  w = cm->nodemap[midnode];
  y = cm->cfirst[w]-1;

  /* 3. Calculate inside up to w, and outside down to y.
   *    We rely on a side effect of how deallocation works
   *    in these routines; the w..y decks are guaranteed
   *    to be retained.
   */
  inside(cm, dsq, L, w, z, i0, j0, BE_EFFICIENT, NULL, &alpha, NULL, &pool, NULL);
  outside(cm, dsq, L, r, y, i0, j0, BE_EFFICIENT, NULL, &beta, pool, NULL);

  /* 4. Find the optimal split: best_v, best_d, best_j
   */
  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  for (v = w; v <= y; v++)
    for (jp = 0; jp <= W; jp++) 
      {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++) 
	  if ((sc = alpha[v][j][d] + beta[v][j][d]) > best_sc)
	    {
	      best_sc = sc;
	      best_v  = v;
	      best_d  = d;
	      best_j  = j;
	    }
      }

  /* free now, before recursing!
   */
  free_vjd_matrix(alpha, w, z, i0, j0);
  free_vjd_matrix(beta, r, v, i0, j0);

  /* 5. The optimal split into a V problem and a wedge problem:
   *    i1 = best_j-best_d+1, j1 = best_j
   *    the V problem:     r..v, i0..i1, j1..j0
   *    the wedge problem: v..z, i1..j1
   *    
   *    These have to solved in the order given because we're
   *    constructing the trace in postorder traversal.
   */
  SQD_DPRINTF1(("Wedge splitter:\n"));
  SQD_DPRINTF1(("   V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		best_v, UniqueStatetype(cm->stid[best_v]),
		i0, best_j-best_d+1, best_j, j0));
  SQD_DPRINTF1(("   wedge:   G%d[%s]..%d[%s], %d..%d\n", 
		best_v, UniqueStatetype(cm->stid[best_v]),
		z, UniqueStatetype(cm->stid[z]),
		best_j-best_d+1, best_j));

  v_splitter(cm, dsq, L, tr, r, best_v, i0, best_j-best_d+1, best_j, j0);
  wedge_splitter(cm, dsq, L, tr, best_v, z, best_j-best_d+1, best_j);

  return best_sc;
}



/* Function: v_splitter()
 * Date:     SRE, Thu May 31 19:47:57 2001 [Kaldi's]
 *
 * Purpose:  Solve a "V problem": best parse of an unbifurcated
 *           subgraph cm^r..z to a one-hole subsequence
 *           i0..i1 // j1..j0. 
 *           
 *           Attaches the optimal trace T(r..z), exclusive of
 *           r, inclusive of z, to the growing trace tr.
 *           
 *           r and z can be any non-insert state. 
 *
 * Args:     cm -  model
 *           dsq - digitized sequence 1..L
 *           L   - length of dsq
 *           tr  - the traceback we're adding on to.
 *           r   - index of the first state in the subgraph 
 *           z   - index of the last state in the subgraph
 *           i0,i1 - first part of the subsequence (1..L)
 *           j1,j0 - second part of the subsequence (1..L)
 * 
 * Returns:  (void)
 */
static void
v_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
	   int r, int z, int i0, int i1, int j1, int j0)
{
  float ***alpha, ***beta;      /* inside and outside matrices */
  struct deckpool_s *pool;      /* pool for holding alloced decks */
  float sc;			/* tmp variable holding a score */
  int   v,w,y;			/* state indexes */
  int   ip,jp;
  int   best_v;
  int   best_i, best_j;		/* optimal i', j' split point */
  float best_sc;		/* score at optimal split point */
  int   midnode;

  /* 1. If the V problem is either a boundary condition, or small
   *    enough, solve it with v_inside^T and append the trace to tr.
   */
   if (cm->ndidx[z] == cm->ndidx[r] + 1 ||
      vinsideT_size(cm, r, z, i0, i1, j1, j0) < RAMLIMIT)
    {
      SQD_DPRINTF1(("Solving a V:   G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		z, UniqueStatetype(cm->stid[z]),
		i0,j1,j1,j0));
      vinsideT(cm, dsq, L, tr, r, z, i0, i1, j1, j0);
      return;
    }

  /* 2. Find our split set, w..y.
   *    Choose the node in the middle.
   */
  midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r]) / 2);
  w = cm->nodemap[midnode];
  y = cm->cfirst[w]-1;

  /* 3. Calculate v_inside up to w, and v_outside down to y.
   *    As with wedge_splitter(), we rely on a side effect of how
   *    deallocation works, so the w..y decks are retained
   *    in alpha and beta even though we're in small memory mode.
   */
  vinside (cm, dsq, L, w, z, i0, i1, j1, j0, BE_EFFICIENT, NULL, &alpha, NULL, &pool, NULL);
  voutside(cm, dsq, L, r, y, i0, i1, j1, j0, BE_EFFICIENT, NULL, &beta,  pool, NULL);

  /* 4. Find the optimal split: v, ip, jp. 
   */
  best_sc = IMPOSSIBLE;
  for (v = w; v <= y; v++)
    for (ip = 0; ip <= i1-i0; ip++)
      for (jp = 0; jp <= j0-j1; jp++)
	if ((sc = alpha[v][jp][ip] + beta[v][jp][ip]) > best_sc)
	  {
	    best_sc = sc;
	    best_v  = v;
	    best_i  = ip + i0;
	    best_j  = jp + j1;
	  }

  /* Free now, before recursing!
   */
  free_vji_matrix(alpha, r, z, j1, j0);
  free_vji_matrix(beta,  r, z, j1, j0);

  /* The optimal split into two V problems:
   *    V:   r..v, i0..i', j'..j0
   *    V:   v..z, i'..i1, j1..j'
   * Solve in this order, because we're constructing the
   * trace in postorder traversal.
   */
  SQD_DPRINTF1(("V splitter:\n"));
  SQD_DPRINTF1(("   V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		best_v, UniqueStatetype(cm->stid[best_v]),
		i0, best_i, best_j, j0));
  SQD_DPRINTF1(("   V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		best_v, UniqueStatetype(cm->stid[best_v]),
		z, UniqueStatetype(cm->stid[z]),
		best_i, i1, j1, best_j));

  v_splitter(cm, dsq, L, tr, r,      best_v, i0,     best_i, best_j, j0);
  v_splitter(cm, dsq, L, tr, best_v, z,      best_i, i1,     j1,     best_j);
  return;
}


/*****************************************************************
 * The alignment engines:
 *     inside   - given generic or wedge problem G^r_z to i0..j0, return score and matrix
 *     outside  - given unbifurcated G^r_z to i0..j0, return matrix
 *     
 *     vinside  - given V problem G^r_z to i0..i1//j1..j0, return score and matrix
 *     voutside - given unbifurcated G^r_z to i0..i1//j1..j0, return matrix
 ******************************************************************/

/* Function: inside()
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
inside(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
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


/* Function: outside()
 * Date:     SRE, Tue Aug  8 10:42:52 2000 [St. Louis]
 *
 * Purpose:  Run the outside version of a CYK alignment algorithm,
 *           on a subsequence i0..j0 of a digitized sequence dsq [1..L],
 *           using a linear segment of a model anchored at a start state
 *           (possibly the absolute root, 0) or (MP,ML,MR,D) and ending at an end
 *           state, bifurcation state, or (MP|ML|MR|D) vend. There must be no
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
 *           vroot     - first state of linear model segment (S; MP|ML|MR|D)
 *           vend      - last state of linear model segment  (B; E; MP|ML|MR|D)
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
outside(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0,
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
  int      w1,w2;		/* bounds of split set */

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

  /* Initialize the root deck.
   * If the root is in a split set, initialize the whole split set.
   */
  w1 = cm->nodemap[cm->ndidx[vroot]]; /* first state in split set */
  w2 = cm->cfirst[w1]-1;	      /* last state in split set w1<=vroot<=w2 */

  for (v = w1; v <= w2; v++) {
    if (! deckpool_pop(dpool, &(beta[v])))
      beta[v] = alloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[v][j][d] = IMPOSSIBLE;
    }
  }
  beta[vroot][j0][W] = 0;		

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0;      v < w1; v++) touch[v] = 0; /* note: top of split set w1, not vroot */
  for (v = vend+1; v < cm->M; v++) touch[v] = 0;
  for (v = w1; v <= vend; v++) {
    if (cm->sttype[v] == B_st) touch[v] = 2; /* well, we'll never use this, but set it anyway. */
    else                       touch[v] = cm->cnum[v];
  }
				

  
  /* Main loop down through the decks
   */
  for (v = w2+1; v <= vend; v++)
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
	      if (y < vroot) continue; /* deal with split sets */
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
    for (v = w1; v <= vend; v++) /* start at w1 - top of split set - not vroot */
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


/* Function: vinside()
 * Date:     SRE, Sat Jun  2 09:24:51 2001 [Kaldi's]
 *
 * Purpose:  Run the inside phase of the CYK alignment algorithm for
 *           a V problem: an unbifurcated CM subgraph from
 *           r..z, aligned to a one-hole subsequence
 *           i0..i1 // j1..j0, exclusive of z,i1,j1.
 *           
 *           This is done in the vji coord system, where
 *           both our j and i coordinates are transformed.
 *           The Platonic matrix runs [j1..j0][i0..i1].
 *           The actual matrix runs [0..j0-j1][0..i1-i0].
 *           To transform a sequence coord i to a transformed
 *           coord i', subtract i0; to transform i' to i,
 *           add i0.
 *           
 *           The conventions for alpha and dpool are the
 *           same as cyk_inside_engine().
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           r         - first start state of subtree (0, for whole model)
 *           z         - last end state of subtree (cm->M-1, for whole model)
 *           i0,i1     - first subseq part of the V problem
 *           j1,j0     - second subseq part 
 *           do_full   - if TRUE, we save all the decks in alpha, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the uppermost deck (r) is valid
 *                       at the end.
 *           a         - if non-NULL, this is an existing matrix, with NULL
 *                       decks for r..z, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_a     - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated vji decks sized
 *                       for this subsequence i0..i1//j0..j1.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..i1//j0..j1 subseq
 *                       because of the size of the subseq decks.
 *           ret_shadow- if non-NULL, the caller wants a shadow matrix, because
 *                       he intends to do a traceback.
 * 
 * Returns:  score.
 */
static float
vinside(CM_t *cm, char *dsq, int L, 
	int r, int z, int i0, int i1, int j1, int j0,
	int do_full, float ***a, float ****ret_a,
	struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	char ****ret_shadow)
{
  char  ***shadow;              /* the shadow matrix -- traceback ptrs -- memory is kept */
  int     v,i,j;
  int     w1,w2;		/* bounds of the split set */
  int     jp, ip;		/* j' and i' -- in the matrix coords */
  int    *touch;                /* keeps track of whether we can free a deck yet or not */
  int     y, yoffset;
  float   sc;			/* tmp variable holding a score */

  /* Allocations, initializations
   */
  if (dpool == NULL) dpool = deckpool_create();
  
  if (a == NULL) {
    a = MallocOrDie(sizeof(float **) * cm->M);
    for (v = 0; v < cm->M; v++) a[v] = NULL;
  }
				/* the whole split set w<=z<=y must be initialized */
  w1 = cm->nodemap[cm->ndidx[z]];
  w2 = cm->cfirst[w1]-1;
  for (v = w1; v <= w2; v++) { 
    if (! deckpool_pop(dpool, &(a[v]))) 
      a[v] = alloc_vji_deck(i0, i1, j1, j0);
    for (jp = 0; jp <= j0-j1; jp++) 
      for (ip = 0; ip <= i1-i0; ip++) 
	a[v][jp][ip] = IMPOSSIBLE;
  }
  a[z][0][i1-i0] = 0.;

  if (ret_shadow != NULL) {
    shadow = MallocOrDie(sizeof(char **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL; 
  }

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0;   v < r;  v++) touch[v] = 0;
  for (v = r;   v <= w2; v++) touch[v] = cm->pnum[v]; /* note w2 not z: to bottom of split set */
  for (v = w2+1; v < cm->M; v++) touch[v] = 0;

  /* Main recursion
   */
  for (v = w1-1; v >= r; v--)
    {
      /* Get a deck and a shadow deck.
       */
      if (! deckpool_pop(dpool, &(a[v]))) 
	a[v] = alloc_vji_deck(i0,i1,j1,j0);
      if (ret_shadow != NULL) 
	shadow[v] = alloc_vji_shadow_deck(i0,i1,j1,j0);

				/* reassert our definition of a V problem */
      if (cm->sttype[v] == E_st || cm->sttype[v] == B_st || (cm->sttype[v] == S_st && v > r))
	Die("you told me you wouldn't ever do that again.");
      
      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= j0-j1; jp++) 
	    for (ip = i1-i0; ip >= 0; ip--) {
	      y = cm->cfirst[v];
	      a[v][jp][ip]      = a[y][jp][ip] + cm->tsc[v][0];
	      if (ret_shadow != NULL) shadow[v][jp][ip] = (char) 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc;
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
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
	      if (ret_shadow != NULL) shadow[v][jp][ip] = (char) 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		   { 
		     a[v][jp][ip] = sc; 
		     if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		   }
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
	      if (ret_shadow != NULL) shadow[v][jp][ip] = 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc; 
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
	      
	      if (dsq[i] < Alphabet_size)
		a[v][jp][ip] += cm->esc[v][(int) dsq[i]];
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
	      if (ret_shadow != NULL) shadow[v][jp][ip] = 0;
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc; 
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
	      
	      if (dsq[j] < Alphabet_size)
		a[v][jp][ip] += cm->esc[v][(int) dsq[j]];
	      else
		a[v][jp][ip] += DegenerateSingletScore(cm->esc[v], dsq[j]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} /* finished calculating deck v */

      /* Now, try to reuse memory under v.
       */
      if (! do_full) {
	for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	  {
	    touch[y]--;
	    if (touch[y] == 0) { 
	      deckpool_push(dpool, a[y]);
	      a[y] = NULL;
	    }
	  }
      }
    } /* end loop over v; we now have a complete matrix */
	
  /* Keep the score.
   */
  sc = a[r][j0-j1][0];

  /* If the caller doesn't want the score matrix back, blow
   * it away (saving decks in the pool). Else, pass it back.
   */
  if (ret_a == NULL) {
    for (v = r; v <= w2; v++)	/* note: go all the way to the bottom of the split set */
      if (a[v] != NULL) {
	deckpool_push(dpool, a[v]);
	a[v] = NULL;
      }
    free(a);
  } else *ret_a = a;
    
  /* If caller doesn't want the deck pool, blow it away.
   * Else, pass it back.
   */
  if (ret_dpool == NULL) {
    float **foo;
    while (deckpool_pop(dpool, &foo)) 
      free_vji_deck(foo, j1,j0);
    deckpool_free(dpool);
  } else *ret_dpool = dpool;

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  return sc;
}


/* Function: voutside()
 * Date:     SRE, Sun Jun  3 15:44:41 2001 [St. Louis]
 *
 * Purpose:  Run the outside version of a CYK alignment algorithm for
 *           a V problem: an unbifurcated CM subgraph from r..z, aligned
 *           to a one-whole subsequence i0..i1//j1..j0, exclusive of
 *           z, i1, j1.
 *           
 *           This is done in the vji coordinate system, where both
 *           our j and i coordinates are transformed. The Platonic
 *           ideal matrix runs [j1..j0][i0..i1]. The implemented
 *           matrix runs [0..j0-j1][0..i1-i0].
 *           
 *           Much of the behavior in calling conventions, etc., is
 *           analogous to the cyk_inside_engine(); see its preface
 *           for more info. Unlike the inside engines, we never 
 *           need to calculate a shadow matrix - outside engines are
 *           only used for divide and conquer steps.
 *
 * Args:     cm        - the model    [0..M-1]
 *           dsq       - the sequence [1..L]   
 *           L         - length of the dsq
 *           r         - first state of linear model segment (S; MP, ML, MR, or D)
 *           z         - last state of linear model segment (B; MP, ML, MR, or D)
 *           i0,i1     - subsequence before the hole  (1..L)
 *           j1,j0     - subsequence after the hole (1..L)
 *           do_full   - if TRUE, we save all the decks in beta, instead of
 *                       working in our default memory-efficient mode where 
 *                       we reuse decks and only the lowermost decks (inc. z) are valid
 *                       at the end.
 *           beta      - if non-NULL, this is an existing matrix, with NULL
 *                       decks for r..z, and we'll fill in those decks
 *                       appropriately instead of creating a new matrix
 *           ret_beta  - if non-NULL, return the matrix with one or more
 *                       decks available for examination (see "do_full")
 *           dpool     - if non-NULL, this is an existing deck pool, possibly empty,
 *                       but usually containing one or more allocated vji decks sized
 *                       for this subsequence i0..i1//j1..j0.
 *           ret_dpool - if non-NULL, return the deck pool for reuse -- these will
 *                       *only* be valid on exactly the same i0..i1//j1..j0 subseq,
 *                       because of the size of the subseq decks.
 */
static void
voutside(CM_t *cm, char *dsq, int L, 
	 int r, int z, int i0, int i1, int j1, int j0,
	 int do_full, float ***beta, float ****ret_beta,
	 struct deckpool_s *dpool, struct deckpool_s **ret_dpool)
{
  int      v,y;			/* indices for states */
  int      i,j;			/* indices in sequence dimensions */
  int      ip, jp;		/* transformed sequence indices */
  float    sc;			/* a temporary variable holding a score */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */
  int      voffset;		/* index of v in t_v(y) transition scores */

  /* Allocations and initializations
   */
  			/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

                        /* if caller didn't give us a matrix, make one */
  if (beta == NULL) {
    beta = MallocOrDie(sizeof(float **) * cm->M);
    for (v = 0; v < cm->M; v++) beta[v] = NULL;
  }
  /* Initialize the root deck. This probably isn't the most efficient way to do it.
   */
  if (! deckpool_pop(dpool, &(beta[r])))
    beta[r] = alloc_vji_deck(i0,i1,j1,j0);
  for (jp = 0; jp <= j0-j1; jp++) {
    for (ip = 0; ip <= i1-i0; ip++)
      beta[r][jp][ip] = IMPOSSIBLE;
  }
  beta[r][j0-j1][0] = 0;		

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0;   v < r;     v++) touch[v] = 0;
  for (v = z+1; v < cm->M; v++) touch[v] = 0;
  for (v = r;   v <= z;    v++) {
    if (cm->sttype[v] == B_st) touch[v] = 2; /* well, we never use this, but be complete */
    else                       touch[v] = cm->cnum[v];
  }
  
  /* Main loop down through the decks
   */
  for (v = r+1; v <= z; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_vji_deck(i0,i1,j1,j0);

      /* main recursion:
       */
      for (jp = j0-j1; jp >= 0; jp--) {
	j = jp+j1;
	for (ip = 0; ip <= i1-i0; ip++) 
	  {
	    i = ip+i0;
	    beta[v][jp][ip] = IMPOSSIBLE;

	    for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	      if (y < r) continue; /* deal with split sets */
	      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

	      switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == j0 || i == i0) continue; /* boundary condition */

		if (dsq[i-1] < Alphabet_size && dsq[j+1] > Alphabet_size)
		  escore = cm->esc[y][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->esc[y], dsq[i-1], dsq[j+1]);
		
		if ((sc = beta[y][jp+1][ip-1]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;

	      case ML_st:
	      case IL_st: 
		if (i == i0) continue;	/* boundary condition */

		if (dsq[i-1] < Alphabet_size) 
		  escore = cm->esc[y][(int) dsq[i-1]];
		else
		  escore = DegenerateSingletScore(cm->esc[y], dsq[i-1]);
		  
		if ((sc = beta[y][jp][ip-1]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		  
		if (dsq[j+1] < Alphabet_size) 
		  escore = cm->esc[y][(int) dsq[j+1]];
		else
		  escore = DegenerateSingletScore(cm->esc[y], dsq[j+1]);

		if ((sc = beta[y][jp+1][ip]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		if ((sc = beta[y][jp][ip] + cm->tsc[y][voffset]) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;

	      default: Die("bogus parent state %d\n", cm->sttype[y]);
	      }/* end switch over states*/
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (beta[v][jp][ip] < IMPOSSIBLE) beta[v][jp][ip] = IMPOSSIBLE;

	  } /* ends loop over ip. We know all beta[v][jp][ip] in this row jp */

      }/* end loop over jp. We know the beta's for the whole deck.*/

      /* Finished deck v.
       * now look at its parents; if we're reusing memory (! do_full)
       * push the parents that we don't need any more into the pool.
       */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { 
	    deckpool_push(dpool, beta[y]); 
	    beta[y] = NULL; 
	  }
	}
      }

    } /* end loop over decks v. */

  /* If the caller doesn't want the matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation!)
   */
  if (ret_beta == NULL) {
    for (v = r; v <= z; v++)
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  /* If the caller doesn't want the deck pool, free it. 
   * Else, pass it back to him.
   */
  if (ret_dpool == NULL) {
    float **a;
    while (deckpool_pop(dpool, &a)) 
      free_vji_deck(a,j1,j0);
    deckpool_free(dpool);
  } else *ret_dpool = dpool;

  free(touch);
}

/*****************************************************************
 * The traceback routines
 *   insideT  - run inside(), append trace in postorder traversal
 *   vinsideT - run vinside(), append trace in postorder traversal
 *****************************************************************/

/* Function: insideT()
 * Date:     SRE, Fri Aug 11 12:08:18 2000 [Pittsburgh]
 *
 * Purpose:  Call inside, get vjd shadow matrix;
 *           then trace back. Append the trace to a given
 *           traceback, which already has state r at tr->n-1.
 */
static float
insideT(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
	int r, int z, int i0, int j0)
{
  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  Nstack_t *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;

  sc = inside(cm, dsq, L, r, z, i0, j0, 
	      BE_EFFICIENT,	/* memory-saving mode */
	      NULL, NULL,	/* manage your own matrix, I don't want it */
	      NULL, NULL,	/* manage your own deckpool, I don't want it */
	      &shadow);		/* return a shadow matrix to me. */

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
  free_vjd_shadow_matrix(shadow, cm, r, z, i0, j0);
  return sc;
}

/* Function: vinsideT()
 * Date:     SRE, Sat Jun  2 14:40:13 2001 [St. Louis]
 *
 * Purpose:  Call vinside(), get vji shadow matrix for a V problem;
 *           then trace back. Append the trace to a
 *           given traceback, which has state r already at
 *           t->n-1.
 */
static float
vinsideT(CM_t *cm, char *dsq, int L, Parsetree_t *tr, 
	 int r, int z, int i0, int i1, int j1, int j0)
{
  char ***shadow;
  float   sc;
  int     v,y;
  int     j,i;
  int     jp,ip;
  int     yoffset;

  sc = vinside(cm, dsq, L, r, z, i0, i1, j1, j0,
	       BE_EFFICIENT,	/* memory-saving mode */
	       NULL, NULL,	/* manage your own matrix, I don't want it */
	       NULL, NULL,	/* manage your own deckpool, I don't want it */
	       &shadow);	/* return a shadow matrix to me. */

  /* We've got a complete shadow matrix. Trace it back. We know
   * that the trace will begin with the start state r, at i0,j0
   * (e.g. jp=j0-j1, ip=0)
   */
  v = r;
  j = j0;
  i = i0;
  while (v < z) {
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
   * state z for i1,j1; it is in the traceback at node tr->n-1.
   */
  free_vji_shadow_matrix(shadow, r, z, j1, j0);
  return sc;
}

/*****************************************************************
 * The size calculators:
 *    insideT_size()   - Mb required by insideT
 *    vinsideT_size()  - Mb required by vinsideT
 *****************************************************************/ 

/* Function: insideT_size()
 * Date:     SRE, Sun Jun  3 17:56:08 2001 [St. Louis]
 *
 * Purpose:  Calculate the # of Mb required to run insideT()
 *           and solve a generic or wedge problem without any
 *           more divide/conquer.
 */
static float
insideT_size(CM_t *cm, int L, int r, int z, int i0, int j0)
{
  float Mb;
  int   maxdecks;
  int   nends;
  int   nbif;

  nends = CMSegmentCountStatetype(cm, r, z, E_st);
  nbif  = CMSegmentCountStatetype(cm, r, z, B_st);
  maxdecks = cyk_deck_count(cm, r, z);

  Mb = (float) (sizeof(float **) * cm->M) / 1000000.;  /* the score matrix */
  Mb += (float) maxdecks * size_vjd_deck(L, i0, j0);
  Mb += (float) (sizeof(int) * cm->M) / 1000000.;      /* the touch array */

  Mb += (float) (sizeof(void **) * cm->M) / 1000000.;
  Mb += (float) (z-r+1-nends-nbif) * size_vjd_yshadow_deck(L, i0, j0);
  Mb += (float) nbif * size_vjd_kshadow_deck(L, i0, j0);

  return Mb;
}

static float
vinsideT_size(CM_t *cm, int r, int z, int i0, int i1, int j1, int j0)
{
  float Mb;
  int   maxdecks;

  Mb = (float) (sizeof(float **) * cm->M) / 1000000.;
  maxdecks = cyk_deck_count(cm, r, z);
  Mb += maxdecks * size_vji_deck(i0,i1,j1,j0);
  Mb += (float)(z-r) * size_vji_shadow_deck(i0,i1,j1,j0);
  return Mb;
}

/* Function: cyk_deck_count()
 * Date:     SRE, Sun Jun  3 20:05:18 2001 [St. Louis]
 *
 * Purpose:  calculate and return the maximum number of
 *           decks that would be required in memory to
 *           solve an alignment problem involving a CM
 *           subgraph from r..z.
 */
static int
cyk_deck_count(CM_t *cm, int r, int z)
{
  Nstack_t *pda;		/* pushdown stack simulating the deck pool */
  int       v,w,y;		/* state indices */
  int       nends;
  int       ndecks;
  int      *touch;		/* keeps track of how many higher decks still need this deck */

  /* Initializations, mirroring key parts of CYKInside()
   */
  ndecks = 1;			/* deck z, which we always need to start with. */
  nends  = CMSegmentCountStatetype(cm, r, z, E_st);
  pda    = CreateNstack();
  touch  = MallocOrDie(sizeof(int) * cm->M);
  for (v = r; v < z; v++)
    touch[v] = cm->pnum[v];

  for (v = z; v >= r; v--)
    {
      if (cm->sttype[v] != E_st) {
	if (! PopNstack(pda, &y)) ndecks++; /* simulated allocation of a new deck */
      }
      
      if (cm->sttype[v] == B_st) { /* release both S children of a bifurc */
	w = cm->cfirst[v];
	y = cm->cnum[v];
	PushNstack(pda, w);
	PushNstack(pda, y);
      } else {
	for (w = cm->cfirst[v]; w < cm->cfirst[v]+cm->cnum[v]; w++)
	  {
	    touch[w]--;
	    if (touch[w] == 0) 
	      {
		if (cm->sttype[w] == E_st) { 
		  nends--; 
		  if (nends == 0) { PushNstack(pda, cm->M-1); }
		} else 
		  PushNstack(pda, w);
	      }
	  }
      }
    }
  FreeNstack(pda);
  free(touch);
  return ndecks;
}

/*################################################################
 * The memory management routines.
 ################################################################*/

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

  SQD_DPRINTF3(("alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  a = MallocOrDie(sizeof(float *) * (L+1)); /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0;   jp < i-1;    jp++) a[jp]     = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp]     = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) a[jp+i-1] = MallocOrDie(sizeof(float) * (jp+1));
  return a;
}
static float
size_vjd_deck(int L, int i, int j)
{
  float Mb;
  int   jp;
  Mb = (float) (sizeof(float *) * (L+1));
  for (jp = 0; jp <= j-i+1; jp++)
    Mb += (float) (sizeof(float) * (jp+1));
  return (Mb / 1000000.);
}
static void
free_vjd_deck(float **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp+i-1]);
  free(a);
}
static void
free_vjd_matrix(float ***a, int r, int z, int i, int j)
{
  int v;
  for (v = r; v <= z; v++)
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
static float
size_vjd_yshadow_deck(int L, int i, int j)
{
  float  Mb;
  int    jp;
  Mb = (float) (sizeof(char *) * (L+1));
  for (jp = 0; jp <= j-i+1; jp++) 
    Mb += (float) (sizeof(char) * (jp+1));
  return Mb / 1000000.;
}
static void
free_vjd_yshadow_deck(char **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp+i-1]);
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
static float
size_vjd_kshadow_deck(int L, int i, int j)
{
  float Mb;
  int   jp;
  
  Mb = (float)(sizeof(int *) * (L+1)); 
  for (jp = 0;   jp <= j-i+1; jp++)
    Mb += (float) (sizeof(int) * (jp+1));
  return Mb / 1000000.;
}
static void
free_vjd_kshadow_deck(int **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp]);
  free(a);
}
static void
free_vjd_shadow_matrix(void ***shadow, CM_t *cm, int r, int z, int i, int j)
{
  int v;
  for (v = r; v <= z; v++)
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
 *           only for solving V problems, they work only
 *           on a defined cube in the 3D matrix: they need
 *           two triplets (r, i0, j0), (z, i1, j1) 
 *           defining the known optimal endpoints of a segment from
 *           an S state to a B state.
 *
 *           By definition of V problems, there's no B states
 *           in between, so the shadow matrix doesn't need any
 *           special casting tricks the way the more generally
 *           used vjd system does.
 */
static float **                 /* allocation of a score deck. */
alloc_vji_deck(int i0, int i1, int j1, int j0)
{
  float **a;
  int     jp;
  a = MallocOrDie(sizeof(float *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    a[jp] = MallocOrDie(sizeof(float)*(i1-i0+1));
  return a;
}
static float
size_vji_deck(int i0, int i1, int j1, int j0)
{
  float Mb;
  int   jp;
  Mb = (float)(sizeof(float *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    Mb += (float)(sizeof(float)*(i1-i0+1));
  return Mb / 1000000.;
}
static void			/* free'ing a score deck */
free_vji_deck(float **a, int j1, int j0)
{
  int jp;
  for (jp = 0; jp <= j0-j1; jp++) 
    if (a[jp] != NULL) free(a[jp]);
  free(a);
}
static void
free_vji_matrix(float ***a, int r, int z, int j1, int j0)
{
  int v;
  for (v = r; v <= z; v++) 
    if (a[v] != NULL) free_vji_deck(a[v], j1, j0);
  free(a);
}
static char **		        /* allocation of a traceback ptr (shadow matrix) deck */
alloc_vji_shadow_deck(int i0, int i1, int j1, int j0)
{
  char **a;
  int     jp;
  a = MallocOrDie(sizeof(char *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    a[jp] = MallocOrDie(sizeof(char)*(i1-i0+1));
  return a;
}
static float		        /* allocation of a traceback ptr (shadow matrix) deck */
size_vji_shadow_deck(int i0, int i1, int j1, int j0)
{
  float   Mb;
  int     jp;
  Mb = (float)(sizeof(char *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    Mb += (float)(sizeof(char)*(i1-i0+1));
  return Mb / 1000000;
}
static void	                /* free'ing a shadow deck */
free_vji_shadow_deck(char **a, int j1, int j0)
{
  int jp;
  for (jp = 0; jp <= j0-j1; jp++) 
    if (a[jp] != NULL) free(a[jp]);
  free(a);
}
static void
free_vji_shadow_matrix(char ***a, int r, int z, int j1, int j0)
{
  int v;
  for (v = r; v <= z; v++) 
    if (a[v] != NULL) free_vji_shadow_deck(a[v], j1, j0);
  free(a);
}


/*################################################################
 * Unused code - 
 *     a reference implementation of the real Outside() algorithm,
 *     including bifurcations. 
 *################################################################*/     
#if 0
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
#endif 

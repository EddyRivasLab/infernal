/* cm_dpsmall.c  (formerly smallcyk.c)
 * SRE, Wed Aug  2 08:42:49 2000 [St. Louis]
 * 
 * Alignment of a CM to a target (sub)sequence.
 *
 * Implementation of the CM divide and conquer alignment algorithm 
 * described in [Eddy02]. Also implements standard CYK/Inside 
 * optimal alignment by dynamic programming [Durbin98]. 
 *
 * These algorithms align to the entire target (sub)sequence
 * (e.g. global alignment). For sequence-local alignment, see
 * scancyk.c.
 */

/*################################################################
 * smallcyk's external API:
 * 
 * CYKDivideAndConquer()    - The divide and conquer algorithm. Align
 *                            a model to a (sub)sequence.
 * CYKInside()              - Align model to (sub)sequence, using normal 
 *                            CYK/Inside algorithm.
 * CYKInsideScore()         - Calculate the CYK/Inside score of optimal 
 *                            alignment, without recovering the alignment; 
 *                            allows timing CYK/Inside without blowing
 *                            out memory, for large target RNAs.
 *                          
 * CYKDemands()             - Print a bunch of info comparing predicted d&c
 *                            time/memory requirements to standard CYK/inside
 *                            time/memory requirements.
 * 
 * All of these functions can take query dependent bands (dmin
 * and dmax) or have them passed as NULL.				
 *################################################################
 */  

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

/* The dividers and conquerors.
 */
static float generic_splitter(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			      int r, int vend, int i0, int j0);
static float wedge_splitter(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			    int r, int z, int i0, int j0);
static void  v_splitter(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			int r, int z, int i0, int i1, int j1, int j0, int useEL);

/* The alignment engines. 
 */
static float inside(CM_t *cm, ESL_DSQ *dsq, int L,
		    int r, int z, int i0, int j0, int do_full,
		    float ***alpha, float ****ret_alpha, 
		    struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		    void ****ret_shadow, int allow_begin, int *ret_b, float *ret_bsc);
static void  outside(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
		     int do_full, float ***beta, float ****ret_beta,
		     struct deckpool_s *dpool, struct deckpool_s **ret_dpool);
static float vinside(CM_t *cm, ESL_DSQ *dsq, int L,
		     int r, int z, int i0, int i1, int j1, int j0, int useEL,
		     int do_full, float ***a, float ****ret_a,
		     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		     char ****ret_shadow,
		     int allow_begin, int *ret_b, float *ret_bsc);
static void  voutside(CM_t *cm, ESL_DSQ *dsq, int L, 
		      int r, int z, int i0, int i1, int j1, int j0, int useEL,
		      int do_full, float ***beta, float ****ret_beta,
		      struct deckpool_s *dpool, struct deckpool_s **ret_dpool);

/* The traceback routines.
 */
static float insideT(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
		     int r, int z, int i0, int j0, int allow_begin, 
		     int *dmin, int *dmax);
static float vinsideT(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
		      int r, int z, int i0, int i1, int j1, int j0, int useEL, 
		      int allow_begin, int *dmin, int *dmax);

/* The size calculators.
 */
float insideT_size(CM_t *cm, int L, int r, int z, int i0, int j0);
float vinsideT_size(CM_t *cm, int r, int z, int i0, int i1, int j1, int j0);
static int   cyk_deck_count(CM_t *cm, int r, int z);
static int   cyk_extra_decks(CM_t *cm);

/* The memory management routines are in infernal.h so hmmband.c can access them 
 */

/*******************************************************************************
 * EPN: QDB-banded D&C engines are named *_qdb() (renamed from *_b() 2026).
 * HMM-banded D&C engines (to be added) will be named *_hb().
 * Functions that I don't think need a banded version are indicated with a U
 * before their names.
 * 
 * To change *most* of the following code from banded to normal versions, two
 * 'replace-string's would be done : 
 * (1) replace '_b(' with '(' : to replace all banded function calls with calls
 *     to their non-banded versions.
 * (2) replace ', dmin, dmax)' with ')' : all banded functions have exactly
 *     two extra variables passed in, dmin a pointer to an int array with minimum
 *     bands, and dmax, a pointer to an int array with maximum bands.  Further,
 *     these are always the last two variables passed into a function.
 *
 * There are two classes of changes that were made to the original functions
 * to make (what I think are) functioning banded versions (*_qdb()).
 *
 * Class 1 : vjd deck changes - using dmin and dmax as bands
 * Class 2 : vji deck changes - using imin and imax (derived from dmin and dmax)
 *
 * Class 2 changes occur only within v problems, only functions : v_splitter_qdb(),
 * vinside_qdb(), and voutside_qdb().
 * 
 * The class 1 changes are more straightforward relative to the class 2 changes.
 * This is completely due to the fact that the vjd coordinate system directly
 * uses d (distance of subsequence in parse tree rooted at state v) which 
 * corresponds conveniently with dmin and dmax.
 * 
 * Class 1 changes are usually involved with a for loop that involves
 * the d index in either the alpha or the beta matrix.  The original for loops
 * are simply replaced with a new for loop that enforces the bands.
 *
 * Class 2 changes that involve the vji decks involve several offset variables
 * because the implicit d value for a given vji cell has to be calculated.  The 
 * formula for that conversion is simple :   d = j-i+1
 * in the code however, jp and ip are used where jp = j-j1 and ip = i-i0.
 * so we have :  d = (jp+j1) - (ip+i0) + 1
 *
 * The way this is handled is only one possible way (and not necessarily the best way) 
 * but saves some calculations from being repeated and is somewhat consistent with
 * analagous code elsewhere.  Also the way it's handled here is somewhat general
 * and could be easily changed. 
 *
 * That approach is to use an imin[] and imax[] vector, somewhat analagous to
 * dmin[] and dmax[], indexed by states where states in the imin
 * and imax vectors are offset (usually by r or w1) because v problems don't involve
 * the entire set of 0..M-1 states.  Because determining a d for a given vji
 * cell depends on both jp and ip, we can't calculate the bands for a given
 * state (vji deck) independent of jp.  Therefore, imin[] and imax[] are calculated
 * independent of jp, and jp must be added within a for(jp...) loop to determine
 * the actual band in the i dimension.  
 *
 * So imin[v-r] = j1-i0-dmax[v]+1;
 *    imax[v-r] = j1-i0-dmin[v]+1;
 *
 * Here's an example of using imin and imax within a for(jp ... ) loop : 
 *	  for (jp = 0; jp <= j0-j1; jp++) 
 *	    {
 * 	      if((imax[v-r]+jp) > (i1-i0)) ip = (i1-i0);
 *	      else ip = imax[v-r] + jp;
 * 	      for(; ip >= imin[v-r]+jp && ip >= 0; ip--) {
 * 
 * Code where bands are used in the vji deck are marked with "Bands used ip X" where X
 * is a number (1-19).  Some of these sections have been commented out as I slowly
 * realized they were mistakes or unnecessary.  There are, admittedly scattered, notes
 * on how I arrived at each of these in :
 * ~nawrocki/lab/rRNA/inf/infernal_0426/banded_testing_0207/00LOG
 * 
 * Other changes of both class 1 and 2 involves imposing the bands during 
 * the initialization step of the alpha matrix.  These changes
 * add additional code that sets all cells outside the bands to IMPOSSIBLE.
 * 
 *******************************************************************************/

/* The QDB-banded dividers and conquerors.
 */
static float generic_splitter_qdb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
				int r, int vend, int i0, int j0, int *dmin, int *dmax);
static float wedge_splitter_qdb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			      int r, int z, int i0, int j0, int *dmin, int *dmax);
static void  v_splitter_qdb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			  int r, int z, int i0, int i1, int j1, int j0, int useEL,
			  int *dmin, int *dmax);

/* The banded alignment engines. 
 */
static float inside_qdb(CM_t *cm, ESL_DSQ *dsq, int L,
		      int r, int z, int i0, int j0, 
		      int do_full,
		      float ***alpha, float ****ret_alpha, 
		      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		      void ****ret_shadow, 
		      int allow_begin, int *ret_b, float *ret_bsc,
		      int *dmin, int *dmax);
static void  outside_qdb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
		       int do_full, float ***beta, float ****ret_beta,
		       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		       int *dmin, int *dmax);
static float vinside_qdb(CM_t *cm, ESL_DSQ *dsq, int L,
		       int r, int z, int i0, int i1, int j1, int j0, int useEL,
		       int do_full, float ***a, float ****ret_a,
		       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		       char ****ret_shadow,
		       int allow_begin, int *ret_b, float *ret_bsc,
		       int *dmin, int *dmax);
static void  voutside_qdb(CM_t *cm, ESL_DSQ *dsq, int L,
			int r, int z, int i0, int i1, int j1, int j0, int useEL,
			int do_full, float ***beta, float ****ret_beta,
			struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			int *dmin, int *dmax);

/* No banded versions of the traceback routines because the non-banded
 * functions can be used.*/

/* No banded size calculators right now. */

/*******************************************************************************
 * EPN 2026: HMM-banded (CP9) divide-and-conquer CYK engines, named *_hb().
 *
 * Stage 1a.1 (brief 26_0610-007): class-1 (vjd-deck) engines are BANDED here, using the
 * CP9/HMM bands in <cp9b> (per-v j-band jmin[v]..jmax[v], and per-(v,j) d-band
 * hdmin[v][j-jmin[v]]..hdmax[v][j-jmin[v]]) -- the SAME bands cm_CYKInsideAlignHB
 * uses. Memory model is option (a): FULL sub-problem-rectangle deck allocation
 * (identical to the NULL-band/qdb D&C, via alloc_vjd_deck()); the bands are
 * enforced ONLY in the DP recursion loops (which (j,d) cells get computed),
 * with all out-of-band cells left IMPOSSIBLE by a full-deck init. Because decks
 * stay full and in platonic [v][j][d] coords, no jp_v/dp_v band-offset
 * translation (as in the mem-efficient CM_HB_MX) is needed: reads of out-of-band
 * child cells simply return IMPOSSIBLE, exactly mirroring cm_CYKInsideAlignHB's
 * skip-out-of-band behavior.
 *
 * Class-2 (vji-deck / v-problem) work is left EXACT (unbanded): the *_hb
 * splitters hand V problems to the EXISTING exact v_splitter() (NULL bands).
 * This makes the *_hb D&C byte-identical to cm_CYKInsideAlignHB only on
 * sequences whose optimum does not clip the bands (banded opt == unbanded opt).
 * Banding the v-problems is Stage 1a.2.
 *******************************************************************************/
static float generic_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
				 int r, int vend, int i0, int j0, CP9Bands_t *cp9b);
static float wedge_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			       int r, int z, int i0, int j0, CP9Bands_t *cp9b);
static float inside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
		       int r, int z, int i0, int j0,
		       int do_full,
		       float ***alpha, float ****ret_alpha,
		       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
		       void ****ret_shadow,
		       int allow_begin, int *ret_b, float *ret_bsc,
		       CP9Bands_t *cp9b);
static void  outside_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
			int do_full, float ***beta, float ****ret_beta,
			struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			CP9Bands_t *cp9b, int **ret_eldmax);
/* brief 26_0430-227: banded EL (local-end) outside deck helpers (defined below,
 * near outside_hb); free_el_banded_vjd_deck is used by the splitters above it. */
static void  free_el_banded_vjd_deck(float **a, int i, int j, int *eldmax);
static float insideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			int r, int z, int i0, int j0, int allow_begin, CP9Bands_t *cp9b);
/* Stage 1a.2 (brief 26_0610-009): banded V-problem (class-2 vji) engines. */
static void  v_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			   int r, int z, int i0, int i1, int j1, int j0, int useEL,
			   CP9Bands_t *cp9b);
static float vinside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
			int r, int z, int i0, int i1, int j1, int j0, int useEL,
			int do_full, float ***a, float ****ret_a,
			char ****ret_shadow,
			int allow_begin, int *ret_b, float *ret_bsc, CP9Bands_t *cp9b);
static void  voutside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
			 int r, int z, int i0, int i1, int j1, int j0, int useEL,
			 int do_full, float ***beta, float ****ret_beta, CP9Bands_t *cp9b);
static float vinsideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			 int r, int z, int i0, int i1, int j1, int j0, int useEL,
			 int allow_begin, CP9Bands_t *cp9b);

/* High-water-mark instrumentation for the D&C-HB live-deck working set.
 * Defined just below; incremented on alloc_vjd_deck()/pop and decremented on
 * push/free within the *_hb path, to measure peak simultaneously-live deck
 * bytes (the capacity-win evidence for brief 26_0610-007). */

/*******************************************************************************
 * 05.24.05
 * EPN MEMORY EFFICIENT BANDED VERSIONS OF SELECTED FUNCTIONS
 * Memory efficient banded functions are named *_b_me()
 * 
 * These functions are modified from their originals to make the memory
 * efficient banded FULL (not D&C) CYK implementation work.  These functions
 * are dubbed 'memory efficient' because they only allocate cells of the
 * alpha or shadow matrix which are within the bands.  The non-memory efficient
 * functions (*_b()) still allocate the same memory as the non-banded functions,
 * but only use the cells within the bands, here we actually don't even allocate
 * unnecessary cells.  The only real difficulty implementing memory efficient
 * bands is in being able to determine what cell alpha[v][j][d] from the 
 * non-memory efficient code corresponds to in the memory-efficient code (we'll
 * call the corresponding cell a[v'][j'][d'] or a[vp][jp][dp]).  The reason
 * v != v'; j != j' and d != d' is because the primes are offset due to the
 * fact that some of the original alpha matrix deck (a[v]) has not been allocated
 * due to the bands.  Therefore all of the differences between the *_b_me() functions
 * and their *_b() versions is to deal with the offset issue.
 * 
 * All changes from the original (non-memory efficient) banded code have been
 * marked with comments beginning 'CYK Full ME Bands Used'.
 *  
 * There are only two functions that need seperate _qdb_me() versions, because
 * the non D&C alignment algorithm only involves three functions, CYKInside(),
 * inside(), and insideT(), and the CYKInside() is really only a wrapper, 
 * for which the memory efficient implementation has no effect, so all we
 * need is inside_qdb_me() and insideT_qdb_me().
 * 
 *******************************************************************************/

/* The alignment engines. 
 */
static float inside_qdb_me(CM_t *cm, ESL_DSQ *dsq, int L,
			 int r, int z, int i0, int j0, 
			 int do_full,
			 float ***alpha, float ****ret_alpha, 
			 void ****ret_shadow, 
			 int allow_begin, int *ret_b, float *ret_bsc,
			 int *dmin, int *dmax);

/* The traceback routines.
 * At first, it wasn't immediately obvious that a *_me version of  
 * this function was needed, but there's some crazy offset issues. [EPN]
 */

static float insideT_qdb_me(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
			  int r, int z, int i0, int j0, int allow_begin,
			  int *dmin, int *dmax);


/* Function: CYKDivideAndConquer()
 * Date:     SRE, Sun Jun  3 19:32:14 2001 [St. Louis]
 *
 * Purpose:  Align a CM to a (sub)sequence using the divide and conquer
 *           algorithm. Return the score (in bits) and a traceback
 *           structure.
 *           
 *           The simplest call to this, for a model cm and a sequence
 *           dsq of length L and no bands on d:
 *               CYKDivideAndConquer(cm, dsq, L, 0, 1, &tr, NULL, NULL);
 *           which will align the model to the entire sequence. (The alignment
 *           will be global w.r.t the sequence.) 
 *           
 *           Sometimes we already know the second state in the traceback:
 *           a CYKScan() will tell us r, for a 0->r local begin transition.
 *           (It also tells us i0, j0: the bounds of a high-scoring subsequence
 *           hit in the target sequence.)  We take all this information in
 *           as a shortcut. The 0->r transition is still counted
 *           towards the score. That is, CYKDivideAndConquer() always
 *           gives a parsetree rooted at state 0, the root, and the sc
 *           we return is the score for that complete parse tree.
 *
 * Args:     cm     - the covariance model
 *           dsq    - the digitized sequence, 1..L
 *           L      - length of sequence
 *           r      - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0     - start of target subsequence (often 1, beginning of sq)
 *           j0     - end of target subsequence (often L, end of sq)
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *           dmin   - minimum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *           dmax   - maximum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *
 * Returns: score of the alignment in bits.  
 */
float
CYKDivideAndConquer(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, 
		    int *dmin, int *dmax)
{
  Parsetree_t *tr;
  float        sc;
  int          z;

  /*printf("alignment strategy:CYKDivideAndConquer:nb:small\n");*/
  /* Trust, but verify.
   * Check out input parameters.
   */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CMH_LOCAL_BEGIN)) cm_Fail("internal error: we're not in local mode, but r is not root");
    if (cm->stid[r] != MATP_MP && cm->stid[r] != MATL_ML &&
	cm->stid[r] != MATR_MR && cm->stid[r] != BIF_B)
      cm_Fail("internal error: trying to do a local begin at a non-mainline start");
  }

  /* Create a parse tree structure.
   * The traceback machinery expects to build on a start state already
   * in the parsetree, so initialize by adding the root state.
   */
  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, i0, j0, 0); /* init: attach the root S */
  z  = cm->M-1;
  sc = 0.;

  /* If r != 0, we already know we're starting with a local entry transition 0->r;
   * add that node too, and count the begin transition towards the score. We have
   * just done our one allowed local begin, so allow_begin becomes FALSE.
   */
  if (r != 0) 
    {
      InsertTraceNode(tr, 0,  TRACE_LEFT_CHILD, i0, j0, r);
      z  =  CMSubtreeFindEnd(cm, r);
      sc =  cm->beginsc[r];
    }

  /* Start the divide and conquer recursion: call the generic_splitter()
   * or generic_splitter_qdb() on the whole DP cube.
   */
  if(dmin == NULL && dmax == NULL)
    sc += generic_splitter(cm, dsq, L, tr, r, z, i0, j0);
  else
    sc += generic_splitter_qdb(cm, dsq, L, tr, r, z, i0, j0, dmin, dmax);
    
  /* Free memory and return
   */
  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  ESL_DPRINTF1(("#DEBUG: returning from CYKDivideAndConquer() sc : %f\n", sc));
  return sc;
}

/*****************************************************************
 * High-water-mark instrumentation for D&C live-deck working set
 * (brief 26_0610-007 capacity-win evidence).
 *
 * When cyk_dnc_track is TRUE, alloc_vjd_deck()/free_vjd_deck() (and the shadow
 * deck allocators/freers) add/subtract the deck's byte size to cyk_dnc_cur_bytes
 * and bump cyk_dnc_max_bytes. This measures the peak number of bytes held in
 * simultaneously-allocated vjd decks -- i.e. the divide-and-conquer working set.
 * Decks parked in a deckpool still count (their memory is still held), which is
 * exactly what we want to compare against the single monolithic CM_HB_MX.
 *****************************************************************/
int    cyk_dnc_track     = FALSE; /* set TRUE to enable byte accounting          */
double cyk_dnc_cur_bytes = 0.0;   /* bytes currently held in class-1 vjd decks   */
double cyk_dnc_vji_bytes = 0.0;   /* bytes currently held in class-2 vji decks   */
double cyk_dnc_max_bytes = 0.0;   /* high-water mark of (vjd + vji) total bytes   */
double cyk_dnc_max_vjd   = 0.0;   /* vjd component at the moment of the high-water */
double cyk_dnc_max_vji   = 0.0;   /* vji component at the moment of the high-water */
double cyk_dnc_shad_cur  = 0.0;   /* brief 26_0610-050: bytes in the current leaf's shadow */
double cyk_dnc_shad_max  = 0.0;   /* brief 26_0610-050: shadow high-water (largest leaf)   */
int    cyk_dnc_vji_row_floats = 0;/* current V-problem row width (i1-i0+1), set at  */
                                  /* vinside()/voutside() entry so free_vji_deck()   */
                                  /* (which lacks i0/i1) can decrement correctly.    */

/* Recompute the live total and, if it sets a new high-water, snapshot the
 * class-1 (banded vjd) vs class-2 (full vji V-problem) split. Called from every
 * deck allocator on the growth path (brief 26_0610-008 measurement). */
void
cyk_dnc_note(void)
{
  double tot = cyk_dnc_cur_bytes + cyk_dnc_vji_bytes;
  if (tot > cyk_dnc_max_bytes) {
    cyk_dnc_max_bytes = tot;
    cyk_dnc_max_vjd   = cyk_dnc_cur_bytes;
    cyk_dnc_max_vji   = cyk_dnc_vji_bytes;
  }
}

void
CYKDeckTrackReset(void)
{
  cyk_dnc_track     = TRUE;
  cyk_dnc_cur_bytes = 0.0;
  cyk_dnc_vji_bytes = 0.0;
  cyk_dnc_max_bytes = 0.0;
  cyk_dnc_max_vjd   = 0.0;
  cyk_dnc_max_vji   = 0.0;
  cyk_dnc_shad_cur  = 0.0;   /* brief 26_0610-050 shadow high-water (see below) */
  cyk_dnc_shad_max  = 0.0;
}
double
CYKDeckTrackMaxMb(void)
{
  cyk_dnc_track = FALSE;
  return cyk_dnc_max_bytes / 1000000.;
}
/* class-1 (banded vjd) bytes held at the moment of the peak total */
double
CYKDeckTrackVjdAtPeakMb(void)
{
  return cyk_dnc_max_vjd / 1000000.;
}
/* class-2 (full vji V-problem) bytes held at the moment of the peak total */
double
CYKDeckTrackVjiAtPeakMb(void)
{
  return cyk_dnc_max_vji / 1000000.;
}

/* brief 26_0610-050 shadow accounting: the SHADOW high-water (yshad/kshad/Lshad/Rshad/
 * Lkmode/Rkmode), separate from the score-deck frontier above. Shadows live only
 * inside a single leaf tr_insideT_hb / tr_vinsideT_hb call (built by the engine,
 * consumed by that call's traceback, then freed); those calls do not nest, so the
 * counter is bumped by the shadow allocators and zeroed by each leaf after it frees
 * its shadow. The max therefore = the largest single-subproblem shadow, NOT the
 * full O(N) parse-tree shadow (the thing 045/046/047 could not bank). */
static void cyk_dnc_shad_add(double bytes)
{ if (! cyk_dnc_track) return; cyk_dnc_shad_cur += bytes; if (cyk_dnc_shad_cur > cyk_dnc_shad_max) cyk_dnc_shad_max = cyk_dnc_shad_cur; }
void   CYKShadowTrackReset(void)      { cyk_dnc_shad_cur = 0.0; cyk_dnc_shad_max = 0.0; }
void   CYKShadowTrackLeafDone(void)   { cyk_dnc_shad_cur = 0.0; }  /* leaf freed its shadow */
double CYKShadowTrackMaxMb(void)      { return cyk_dnc_shad_max / 1000000.; }

/* Function: CYKDivideAndConquerHB()
 * Date:     EPN 2026 [brief 26_0610-007]
 *
 * Purpose:  HMM-banded (CP9) divide-and-conquer CYK alignment. Mirrors
 *           CYKDivideAndConquer() but carries the CP9 bands <cp9b> down the
 *           recursion via generic_splitter_hb(). Class-1 (vjd) sub-problems are
 *           solved with the banded *_hb engines; class-2 (V) sub-problems are
 *           solved EXACTLY with the existing v_splitter() (Stage 1a.1 scope).
 *
 *           Same calling convention as CYKDivideAndConquer() with r usually 0,
 *           i0=1, j0=L for a full global alignment. <cp9b> must already have
 *           been filled for this <dsq>/<L> (e.g. by cp9_Seq2Bands()).
 *
 * Args:     cm     - the covariance model (with cm->cp9b == cp9b typically)
 *           dsq    - digitized sequence 1..L
 *           L      - length of sequence
 *           r      - root state of subgraph (usually 0)
 *           i0     - start of target subsequence (usually 1)
 *           j0     - end of target subsequence (usually L)
 *           ret_tr - RETURN: traceback (NULL if not wanted)
 *           cp9b   - CP9 bands (jmin/jmax + hdmin/hdmax) for dsq/L
 *
 * Returns:  score of the alignment in bits.
 */
float
CYKDivideAndConquerHB(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr,
		      CP9Bands_t *cp9b)
{
  Parsetree_t *tr;
  float        sc;
  int          z;

  if (cp9b == NULL) cm_Fail("CYKDivideAndConquerHB(): cp9b is NULL");

  /* Trust, but verify. (Same checks as CYKDivideAndConquer().) */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CMH_LOCAL_BEGIN)) cm_Fail("internal error: we're not in local mode, but r is not root");
    if (cm->stid[r] != MATP_MP && cm->stid[r] != MATL_ML &&
	cm->stid[r] != MATR_MR && cm->stid[r] != BIF_B)
      cm_Fail("internal error: trying to do a local begin at a non-mainline start");
  }

  tr = CreateParsetree(100);
  InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, i0, j0, 0); /* init: attach the root S */
  z  = cm->M-1;
  sc = 0.;

  if (r != 0)
    {
      InsertTraceNode(tr, 0,  TRACE_LEFT_CHILD, i0, j0, r);
      z  =  CMSubtreeFindEnd(cm, r);
      sc =  cm->beginsc[r];
    }

  /* brief 26_0430-225: DNC_MEM_VERBOSE ground-truth peak-memory diagnostic,
   * mirroring INFERNAL_CKPT_VERBOSE's (cm_dpalign.c) env-gated fprintf style.
   * Reuses the existing (previously dormant -- no caller anywhere before this
   * brief) high-water-mark instrumentation built for brief 26_0610-007/008/050
   * (CYKDeckTrackReset/MaxMb/VjdAtPeakMb/VjiAtPeakMb + the shadow tracker
   * above): that infrastructure already accounts every live banded-vjd deck,
   * class-2 (V-problem) vji deck, and shadow deck byte-for-byte at alloc/free
   * time, so no new byte-counting is added here -- just the env gate + report. */
  int dnc_mem_verbose = (getenv("DNC_MEM_VERBOSE") != NULL) ? TRUE : FALSE;
  if (dnc_mem_verbose) { CYKDeckTrackReset(); CYKShadowTrackReset(); }

  /* Start the banded divide and conquer recursion. */
  sc += generic_splitter_hb(cm, dsq, L, tr, r, z, i0, j0, cp9b);

  if (dnc_mem_verbose) {
    double peak_mb = CYKDeckTrackMaxMb();  /* NOTE: this call also sets cyk_dnc_track = FALSE */
    double vjd_mb  = CYKDeckTrackVjdAtPeakMb();
    double vji_mb  = CYKDeckTrackVjiAtPeakMb();
    double shad_mb = CYKShadowTrackMaxMb();
    fprintf(stderr, "# CYKDivideAndConquerHB (D&C) engaged: M=%d L=%d sc=%.5f  DnC-DP peak=%.2f Mb (vjd=%.2f Mb vji=%.2f Mb)  shadow-peak=%.2f Mb\n",
            cm->M, L, sc, peak_mb, vjd_mb, vji_mb, shad_mb);
  }

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  ESL_DPRINTF1(("#DEBUG: returning from CYKDivideAndConquerHB() sc : %f\n", sc));
  return sc;
}

/* Function: CYKInside()
 * Date:     SRE, Sun Jun  3 19:48:33 2001 [St. Louis]
 *
 * Purpose:  Wrapper for the insideT() routine - solve
 *           a full alignment problem, return the traceback
 *           and the score, without dividing & conquering.
 *           
 *           Analogous to CYKDivideAndConquer() in many respects;
 *           see the more extensive comments in that function for
 *           more details on shared aspects.
 *           
 * Args:     cm     - the covariance model
 *           sq    - the sequence, 1..L
 *           r      - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0     - start of target subsequence (often 1, beginning of sq)
 *           j0     - end of target subsequence (often L, end of sq)
 *           ret_tr - RETURN: traceback (pass NULL if trace isn't wanted)
 *           dmin   - minimum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *           dmax   - maximum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *
 * Returns:  score of the alignment in bits.
 */
float
CYKInside(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr,
	  int *dmin, int *dmax)
{
  Parsetree_t *tr;
  int          z;
  float        sc;

  /* Trust, but verify.
   * Check out input parameters.
   */
  if (cm->stid[r] != ROOT_S) {
    if (! (cm->flags & CMH_LOCAL_BEGIN)) cm_Fail("internal error: we're not in local mode, but r is not root");
    if (cm->stid[r] != MATP_MP && cm->stid[r] != MATL_ML &&
	cm->stid[r] != MATR_MR && cm->stid[r] != BIF_B)
      cm_Fail("internal error: trying to do a local begin at a non-mainline start");
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

  /* Solve the whole thing with one call to insideT.
   */
  /* if we're non-banded use the original function */
  if(dmin == NULL && dmax == NULL)
    sc += insideT(cm, dsq, L, tr, r, z, i0, j0, (r==0), 
		  dmin, dmax);
  /* if we're using query dependent bands, call the 
   * memory efficient QDB alignment version.
   */
  else
    sc += insideT_qdb_me(cm, dsq, L, tr, r, z, i0, j0, (r==0),
      dmin, dmax);
  /* To call the non-memory efficient version, uncomment
   * the following line: */
  /*sc += insideT(cm, dsq, L, tr, r, z, i0, j0, (r==0),    dmin, dmax);*/

  if (ret_tr != NULL) *ret_tr = tr; else FreeParsetree(tr);
  return sc;
}

/* Function: CYKInsideScore()
 * Date:     SRE, Tue Apr  9 05:21:22 2002 [St. Louis]
 *
 * Purpose:  Wrapper for the inside() routine. Solve
 *           a full alignment problem in one pass of inside,
 *           in memory-saving mode, returning only the score.
 *           
 *           Fairly useless. Written just to obtain timings
 *           for SSU and LSU alignments, for comparison to
 *           divide and conquer.
 *
 *           Analogous to CYKDivideAndConquer() in many respects;
 *           see the more extensive comments in that function for
 *           more details on shared aspects.
 *           
 * Args:     cm     - the covariance model
 *           dsq    - the sequence, 1..L
 *           L      - length of sequence
 *           r      - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0     - start of target subsequence (often 1, beginning of sq)
 *           j0     - end of target subsequence (often L, end of sq)
 *           dmin   - minimum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *           dmax   - maximum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *
 * Returns:  score of the alignment in bits.
 */
float
CYKInsideScore(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, int *dmin, int *dmax)
{
  int    z;
  float  sc;

  z           = cm->M-1;
  sc          = 0.;

  if (r != 0) 
    {
      z  =  CMSubtreeFindEnd(cm, r);
      sc =  cm->beginsc[r];
    }

  if(dmin == NULL && dmax == NULL)
    sc +=  inside(cm, dsq, L, r, z, i0, j0, FALSE, 
		  NULL, NULL, NULL, NULL, NULL,
		  (r==0), NULL, NULL);
  else
    sc +=  inside_qdb(cm, dsq, L, r, z, i0, j0, FALSE, 
		    NULL, NULL, NULL, NULL, NULL,
		    (r==0), NULL, NULL, dmin, dmax);

  return sc;
}


/* Function: CYKDemands()
 * Date:     SRE, Sun Jun  3 20:00:54 2001 [St. Louis]
 *
 * Purpose:  Print out information on the computational
 *           complexity of an alignment problem for divide
 *           and conquer versus the full CYK.
 *
 * Args:     cm     - the model
 *           L      - length of sequence.
 *           dmin   - minimum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *           dmax   - maximum d bound for each state v; [0..v..M-1] (NULL if non-banded)
 *           be_quiet - TRUE to not print info, just return number of DP calcs
 * 
 * Returns: (float) the total number of DP calculations, either using QDB (if
 *                  dmin & dmax are non-NULL) or not using QDB.
 */
float
CYKDemands(CM_t *cm, int L, int *dmin, int *dmax, int be_quiet)
{
  float Mb_per_deck;    /* megabytes per deck */
  int   bif_decks;	/* bifurcation decks  */
  int   nends;		/* end decks (only need 1, even for multiple E's */
  int   maxdecks;	/* maximum # of decks needed by CYKInside() */
  int   extradecks;     /* max # of extra decks needed for bifurcs */
  float smallmemory;	/* how much memory small version of CYKInside() needs */
  float bigmemory;	/* how much memory a full CYKInside() would take */
  float dpcells;	/* # of dp cells */
  float bifcalcs;	/* # of inner loops executed for bifurcation calculations */
  float bifcalcs_b;	/* # of inner loops executed for bifurcation calculations in QDB */
  float dpcalcs;	/* # of inner loops executed for non-bif calculations */
  float dpcalcs_b;	/* # of inner loops executed for bifurcation calculations in QDB */
  int   j;
  float avg_Mb_per_banded_deck;    /* average megabytes per deck in mem efficient big mode */
  int   v, y, z, d, kmin, kmax; /* for QDB calculations */

  Mb_per_deck = size_vjd_deck(L, 1, L);
  bif_decks   = CMCountStatetype(cm, B_st);
  nends       = CMCountStatetype(cm, E_st);
  maxdecks    = cyk_deck_count(cm, 0, cm->M-1);
  extradecks  = cyk_extra_decks(cm);
  smallmemory = (float) maxdecks * Mb_per_deck;
  bifcalcs = 0.;
  for (j = 0; j <= L; j++)
    bifcalcs += (float)(j+1)*(float)(j+2)/2.;
  bifcalcs *= (float) bif_decks;
  dpcalcs = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - bif_decks - nends +1);
  if(dmin == NULL && dmax == NULL)
    {
      bigmemory   = (float) (cm->M - nends +1) * Mb_per_deck;
      dpcells     = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - nends +1);
      avg_Mb_per_banded_deck = 0.; /* irrelevant */
    }
  else
    {
      dpcells = 0.;
      dpcalcs_b = 0.;
      for(v = 0; v < cm->M; v++)
	{
	  dpcells   += (float) (L+1) * (float) (dmax[v] - dmin[v] + 1.);
	  if(cm->sttype[v] != B_st)
	    dpcalcs_b   += (float) (L+1) * (float) (dmax[v] - dmin[v] + 1.);
	  for(d = dmin[v]; d <= dmax[v]; d++)
	    {
	      dpcells -= (float) d; /* subtract out cells for which d <= j */
	      if(cm->sttype[v] != B_st)
		dpcalcs_b   -= (float) d; 
	    }
	}
      bigmemory   = (sizeof(float) * dpcells) / 1000000.;
      avg_Mb_per_banded_deck = bigmemory / ((float) cm->M -nends + 1);
      /* bigmemory and avg_Mb_per_banded_deck should be treated as approximates,
       * I'm not sure if they're exactly correct. EPN, Mon Nov  6 07:56:13 2006 */

      /* for QDB, to get bifcalcs, we need to count all the cells within the bands on
       * left and right childs y and z of v, that are consistent with band on v 
       * there's probably a more efficient way of doing this. */
      bifcalcs_b = 0.;
      for (v = 0; v < cm->M; v++)
	{
	  if(cm->sttype[v] == B_st)
	    {
	      y = cm->cfirst[v];
	      z = cm->cnum[v];
	      for (j = 0; j <= L; j++)
		{
		  for (d = dmin[v]; d <= dmax[v] && d <= j; d++)
		    {
		      if(dmin[z] > (d-dmax[y])) kmin = dmin[z];
		      else kmin = d-dmax[y];
		      if(kmin < 0) kmin = 0;
		      if(dmax[z] < (d-dmin[y])) kmax = dmax[z];
		      else kmax = d-dmin[y];
		      if(kmin <= kmax)
			bifcalcs_b += (float)(kmax - kmin + 1);
		    }
		}
	    }
	}
    }

  if(dmin == NULL && dmax == NULL)
    {
      if(!be_quiet)
	{
	  printf("CYK cpu/memory demand estimates:\n");
	  printf("Mb per cyk deck:                  %.4f\n", Mb_per_deck);
	  printf("# of decks (M):                   %d\n",   cm->M);
	  printf("# of decks needed in small CYK:   %d\n",   maxdecks);
	  printf("# of extra decks needed:          %d\n",   extradecks);
	  printf("RAM needed for full CYK, Mb:      %.2f\n", bigmemory);
	  printf("RAM needed for small CYK, Mb:     %.2f\n", smallmemory);
	  printf("# of dp cells, total:             %.3g\n", dpcells);
	  printf("# of non-bifurc dp cells:         %.3g\n", dpcalcs);
	  printf("# of bifurcations:                %d\n",   bif_decks);
	  printf("# of bifurc dp inner loop calcs:  %.3g\n", bifcalcs);
	  printf("# of dp inner loops:              %.3g\n", dpcalcs+bifcalcs);
	}
      return (dpcalcs + bifcalcs);
    }
  else /* QDB */
    {
      if(!be_quiet)
	{
	  printf("QDB CYK cpu/memory demand estimates:\n");
	  printf("Mb per cyk deck:                     %.4f\n", Mb_per_deck);
	  printf("Avg Mb per QDB cyk deck:             %.4f\n", avg_Mb_per_banded_deck);
	  printf("# of decks (M):                      %d\n",   cm->M);
	  printf("# of decks needed in small QDB CYK:  %d\n",   maxdecks);
	  printf("# of extra decks needed:             %d\n",   extradecks);
	  printf("RAM needed for full QDB CYK, Mb:     %.2f\n", bigmemory);
	  printf("RAM needed for small QDB CYK, Mb:    %.2f\n", smallmemory);
	  printf("# of QDB dp cells, total:            %.3g\n", dpcells);
	  printf("# of QDB non-bifurc dp cells:        %.3g\n", dpcalcs_b);
	  printf("# of bifurcations:                   %d\n",   bif_decks);
	  printf("# of QDB bifurc dp inner loop calcs: %.3g\n", bifcalcs_b);
	  printf("# of QDB dp inner loops:             %.3g\n", dpcalcs_b+bifcalcs_b);
	  printf("Estimated small CYK QDB aln speedup: %.4f\n", ((dpcalcs+bifcalcs)/(dpcalcs_b+bifcalcs_b)));
	}
      return (dpcalcs_b + bifcalcs_b);
    }
}


/* Function: CYKNonQDBSmallMbNeeded()
 * Date:     EPN, Fri May 27 11:43:56 2011
 *
 * Purpose:  Return number of Mb needed for non-QDB
 *           divide and conquer CYK.
 *
 * Args:     cm     - the model
 *           L      - length of sequence.
 * 
 * Returns: Number of Mb required.
 */
float
CYKNonQDBSmallMbNeeded(CM_t *cm, int L)
{
  float Mb_per_deck;    /* megabytes per deck */
  int   maxdecks;	/* maximum # of decks needed by CYKInside() */
  float smallmemory;	/* how much memory small version of CYKInside() needs */

  Mb_per_deck = size_vjd_deck(L, 1, L);
  maxdecks    = cyk_deck_count(cm, 0, cm->M-1);
  smallmemory = (float) maxdecks * Mb_per_deck;
  return smallmemory;
}

/*################################################################
 * The dividers and conquerors. 
 *################################################################*/  

/* Function: generic_splitter()
 * Date:     SRE, Sat May 12 15:08:38 2001 [CSHL]
 *
 * Purpose:  Solve a "generic problem": best parse of
 *           a possibly bifurcated subgraph cm^r_z to
 *           a substring sq->sq[i0..j0]. r is usually a start
 *           state (S_st) but may be any non-end state type in 
 *           the case of local alignment begins (ROOT 0->r).
 *           z is always an end state (E_st).
 *
 *           Given: a cm subgraph from r..z
 *                  a subsequence from i0..j0
 *           Attaches the optimal trace T{r..z}, exclusive of r
 *           and inclusive of z, to tr.
 *           
 *           A full divide & conquer never terminates
 *           in generic_splitter; the recursion must
 *           terminate in v_splitter and wedge_splitter;
 *           so we don't test an end-of-recursion boundary.
 *           
 * Args:     cm          - model
 *           sq          - sequence, digitized, 1..L
 *           tr          - the traceback we're adding on to.
 *           r           - index of the root state of this problem in the model       
 *           z           - index of an end state (E_st) in the model
 *           i0          - start in the sequence (1..L)
 *           j0          - end in the sequence (1..L)
 *
 * Returns:  score of the optimal parse of sq(i0..j0) with cm^r_z 
 */
static float
generic_splitter(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
		 int r, int z, int i0, int j0)
{
  float ***alpha;
  float ***beta;
  struct deckpool_s *pool;
  int      v,w,y;		/* state indices */
  int      wend, yend;		/* indices for end of subgraphs rooted at w,y */
  int      jp;			/* j': relative position in subseq, 0..W */
  int      W;			/* length of subseq i0..j0 */
  float    sc;			/* tmp variable for a score */
  int      j,d,k;		/* sequence indices */
  float    best_sc;		/* optimal score at the optimal split point */
  int      best_k;		/* optimal k for the optimal split */
  int      best_d;		/* optimal d for the optimal split */
  int      best_j;		/* optimal j for the optimal split */
  int      tv;			/* remember the position of a bifurc in the trace. */
  int      b1,b2;		/* argmax_v for 0->v local begin transitions */
  float    b1_sc, b2_sc;	/* max_v scores for 0->v local begin transitions */

  /* 1. If the generic problem is small enough, solve it with inside^T,
   *    and append the trace to tr.
   */
  if (insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT) {
    ESL_DPRINTF2(("#DEBUG: Solving a generic w/ insideT - G%d[%s]..%d[%s], %d..%d\n",
		  r, UniqueStatetype(cm->stid[r]),
		  z, UniqueStatetype(cm->stid[z]),
		  i0, j0));
    sc = insideT(cm, dsq, L, tr, r, z, i0, j0, (r==0), 
		NULL, NULL); /* two NULLs mean 'don't use bands' */

    return sc;
  }

  /* 2. Traverse down from r, find first bifurc.
   *    The lowest a bifurc could be: B-S-E/S-IL-E = vend-5
   *                                   
   */
  for (v = r; v <= z-5; v++)
    if (cm->sttype[v] == B_st) break; /* found the first bifurcation, now v */

  /* 3. If there was no bifurcation, this is a wedge problem; solve it
   *    with wedge_splitter. 
   */
  if (v > z-5) {		/* no bifurc? it's a wedge problem  */
    if (cm->sttype[z] != E_st) cm_Fail("inconceivable.");
    sc = wedge_splitter(cm, dsq, L, tr, r, z, i0, j0);
    return sc;
  }

  /* Set up the state quartet r,v,w,y for a divide and conquer
   * solution of the generic problem.
   */
  w = cm->cfirst[v];		/* index of left S  */
  y = cm->cnum[v];		/* index right S    */
  if (w < y) { wend = y-1; yend = z; }
  else       { yend = w-1; wend = z; }

  /* Calculate alpha[w] deck and alpha[y] deck.
   * We also get b1: best choice for 0->b local begin. b1_sc is the score if we do this.
   * Analogous for b2, b2_sc on the other side.
   */
  inside(cm, dsq, L, w, wend, i0, j0, BE_EFFICIENT, NULL,  &alpha, NULL, &pool, NULL, 
	 (r==0), &b1, &b1_sc);
  inside(cm, dsq, L, y, yend, i0, j0, BE_EFFICIENT, alpha, &alpha, pool, &pool, NULL,
	 (r==0), &b2, &b2_sc);

  /* Calculate beta[v] deck (stick it in alpha). Let the pool get free'd.
   * (If we're doing local alignment, deck M is the beta[EL] deck.)
   */
  outside(cm, dsq, L, r, v, i0, j0, BE_EFFICIENT, alpha, &beta, pool, NULL);

  /* Find the optimal split at the B.
   */
  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  for (jp = 0; jp <= W; jp++) 
    {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	for (k = 0; k <= d; k++)
	  if ((sc = alpha[w][j-k][d-k] + alpha[y][j][k] + beta[v][j][d]) > best_sc) 
	    {
	      best_sc = sc;
	      best_k  = k;
	      best_j  = j;
	      best_d  = d;
	    }
    }

  /* Local alignment only: maybe we're better off in EL?
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = 0; jp <= W; jp++) 
      {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--)
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc;
	    best_k  = -1;	/* special flag for local end, EL. */
	    best_j  = j;
	    best_d  = d;
	  }
      }
  }

  /* Local alignment only: maybe we're better off in ROOT?
   */
  if (r == 0 && cm->flags & CMH_LOCAL_BEGIN) {
    if (b1_sc > best_sc) {
      best_sc = b1_sc;
      best_k  = -2;		/* flag for using local begin into left wedge w..wend */
      best_j  = j0;		
      best_d  = W;
    }
    if (b2_sc > best_sc) {
      best_sc = b2_sc;
      best_k  = -3;		/* flag for using local begin into right wedge y..yend */
      best_j  = j0;		
      best_d  = W;
    }
  }

  /* Free now, before recursing.
   * The two alpha matrices and the beta matrix
   * actually all point to the same memory, since no
   * decks in Inside and Outside needed to overlap. 
   * Free 'em all in one call.
   */
  free_vjd_matrix(alpha, cm->M, i0, j0);

  /* If we're in EL, instead of B, the optimal alignment is entirely
   * in a V problem that's still above us. The TRUE flag sets useEL.
   */
  if (best_k == -1) {	
    v_splitter(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, TRUE);    
    return best_sc;
  } 

  /* Else: if we're in the root 0, we know which r we did our local begin into.
   * We have a generic problem rooted there. The FALSE flag disallows
   * any further local begins.
   */
  if (best_k == -2) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b1);
    z = CMSubtreeFindEnd(cm, b1);
    generic_splitter(cm, dsq, L, tr, b1, z, i0, j0);
    return best_sc;
  }
  if (best_k == -3) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b2);
    z = CMSubtreeFindEnd(cm, b2);
    generic_splitter(cm, dsq, L, tr, b2, z, i0, j0);
    return best_sc;
  }

  /* Else (the usual case), ok, we did use B in the optimal split.
   * Split now into a V problem and two generic problems, and recurse
   * left fragment: i1 = j-d+1, j1 = j-k, vroot = w, vend = wend
   * right frag:    i2 = j-k+1, j2 = j,   vroot = y, vend = yend
   * 
   * The problems must be solved in a particular order, since we're
   * constructing the trace in a postorder traversal.
   */
  ESL_DPRINTF2(("#DEBUG: Generic splitter:\n"));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		v, UniqueStatetype(cm->stid[v]),
		i0, best_j-best_d+1, best_j, j0));
  ESL_DPRINTF2(("#DEBUG:    generic: G%d[%s]..%d[%s], %d..%d\n", 
		w,    UniqueStatetype(cm->stid[w]),
		wend, UniqueStatetype(cm->stid[wend]),
		best_j-best_d+1, best_j-best_k));
  ESL_DPRINTF2(("#DEBUG:    generic: G%d[%s]..%d[%s], %d..%d\n", 
		y,    UniqueStatetype(cm->stid[y]),
		yend, UniqueStatetype(cm->stid[yend]),
		best_j-best_k+1, best_j));

  v_splitter(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, FALSE);
  tv = tr->n-1;

  InsertTraceNode(tr, tv, TRACE_LEFT_CHILD, best_j-best_d+1, best_j-best_k, w);
  generic_splitter(cm, dsq, L, tr, w, wend, best_j-best_d+1, best_j-best_k);
  InsertTraceNode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j, y);
  generic_splitter(cm, dsq, L, tr, y, yend, best_j-best_k+1, best_j);

  return best_sc;
}

/* Function: wedge_splitter()
 * Date:     SRE, Sun May 13 08:44:15 2001 [CSHL genome mtg]
 *
 * Purpose:  Solve a "wedge problem": best parse of an 
 *           unbifurcated subgraph cm^r..z to a substring
 *           sq->sq[i0..j0]. r may be a start state (when
 *           the wedge problem comes from being a special case
 *           of a generic problem) or a non-insert state
 *           (D, MP, ML, MR) (when the wedge comes from a
 *           previous wedge_splitter), or indeed, any non-end
 *           state (when wedge comes from a local begin).
 *           z, however, is always an end state.
 *           
 *           Attaches the optimal trace T(r..z), exclusive
 *           of r and inclusive of z, to the growing trace tr.
 *           
 *           Deal with a divide and conquer boundary condition:
 *           the next non-insert state after r is the end state z.
 *           All remaining sequence of i0..j0 that r doesn't emit
 *           must be dealt with by insert states.
 *
 * Args:     cm          - model
 *           sq          - digitized sequence 1..L
 *           tr          - the traceback we're adding on to.
 *           r           - index of the first state in the subgraph
 *           z           - index of an end state (E_st) in the model
 *           i0          - start in the sequence (1..L)
 *           j0          - end in the sequence (1..L)
 *
 * Returns:  The score of the best parse in bits.
 */
static float 
wedge_splitter(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int j0)
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
  int   b;	/* optimal local begin: b = argmax_v alpha_v(i0,j0) + t_0(v) */
  float bsc;	/* score for optimal local begin      */
  
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
      ESL_DPRINTF2(("#DEBUG: Solving a wedge:   G%d[%s]..%d[%s], %d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		z, UniqueStatetype(cm->stid[z]),
		i0,j0));
      sc = insideT(cm, dsq, L, tr, r, z, i0, j0, (r==0),
		   NULL, NULL); /* two NULLs mean 'don't use bands' */

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
   *    b will contain the optimal 0->v state for a local begin, and bsc
   *    is the score for using it.
   *    beta[cm->M] will contain the EL deck, if needed for local ends.
   */
  inside(cm, dsq, L, w, z, i0, j0, BE_EFFICIENT, 
	 NULL, &alpha, NULL, &pool, NULL, 
	 (r==0), &b, &bsc);
  outside(cm, dsq, L, r, y, i0, j0, BE_EFFICIENT, NULL, &beta, pool, NULL);

  /* 4. Find the optimal split at the split set: best_v, best_d, best_j
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

  /* Local alignment ends only: maybe we're better off in EL, 
   * not in the split set?
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = 0; jp <= W; jp++) 
      {
	j = i0-1+jp;
	for (d = 0; d <= jp; d++)
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc;
	    best_v  = -1;	/* flag for local alignment. */
	    best_j  = j;
	    best_d  = d;
	  }
      }
  }

  /* Local alignment begins only: maybe we're better off in the root.
   */
  if (r==0 && (cm->flags & CMH_LOCAL_BEGIN)) {
    if (bsc > best_sc) {
      best_sc = bsc;
      best_v  = -2;		/* flag for local alignment */
      best_j  = j0;
      best_d  = W;
    }
  }

  /* free now, before recursing!
   */
  free_vjd_matrix(alpha, cm->M, i0, j0);
  free_vjd_matrix(beta,  cm->M, i0, j0);

  /* If we're in EL, instead of the split set, the optimal alignment
   * is entirely in a V problem that's still above us. The TRUE
   * flag sets useEL. It doesn't matter which state in the split
   * set w..y we use as the end of the graph; vinside() will have to
   * initialize the whole thing to IMPOSSIBLE anyway.
   */  
  if (best_v == -1) {
    v_splitter(cm, dsq, L, tr, r, w, i0, best_j-best_d+1, best_j, j0, TRUE);    
    return best_sc;
  }

  /* If we're in the root because of a local begin, the local alignment
   * is entirely in a wedge problem that's still below us, rooted at b.
   * The FALSE flag prohibits any more local begins in this and subsequent
   * problems. 
   */
  if (best_v == -2) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b);
    wedge_splitter(cm, dsq, L, tr, b, z, i0, j0);
    return best_sc; 
  }

  /* Else (usual case): the optimal split into a V problem and a wedge problem:
   *    i1 = best_j-best_d+1, j1 = best_j
   *    the V problem:     r..v, i0..i1, j1..j0
   *    the wedge problem: v..z, i1..j1
   *    
   *    These have to solved in the order given because we're
   *    constructing the trace in postorder traversal.
   */
  ESL_DPRINTF2(("#DEBUG: Wedge splitter:\n"));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		best_v, UniqueStatetype(cm->stid[best_v]),
		i0, best_j-best_d+1, best_j, j0));
  ESL_DPRINTF2(("#DEBUG:    wedge:   G%d[%s]..%d[%s], %d..%d\n", 
		best_v, UniqueStatetype(cm->stid[best_v]),
		z, UniqueStatetype(cm->stid[z]),
		best_j-best_d+1, best_j));

  v_splitter(cm, dsq, L, tr, r, best_v, i0, best_j-best_d+1, best_j, j0, FALSE);
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
 * Args:     cm          -  model
 *           sq          - digitized sequence 1..L
 *           tr          - the traceback we're adding on to.
 *           r           - index of the first state in the subgraph 
 *           z           - index of the last state in the subgraph
 *           i0,i1       - first part of the subsequence (1..L)
 *           j1,j0       - second part of the subsequence (1..L)
 *           useEL       - TRUE if i1,j1 aligned to EL, not z
 * 
 * Returns:  (void)
 */
static void
v_splitter(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
	   int r, int z, int i0, int i1, int j1, int j0, 
	   int useEL)
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
  int   b;			/* optimal choice for a 0->b local begin  */
  float bsc;			/* score if we use the local begin */

  /* 1. If the V problem is either a boundary condition, or small
   *    enough, solve it with v_inside^T and append the trace to tr.
   *    (With local alignment, we might even see a lone B state
   *     get handed to v_splitter(); hence the r==z case.)
   */
   if (cm->ndidx[z] == cm->ndidx[r] + 1 || r == z || 
      vinsideT_size(cm, r, z, i0, i1, j1, j0) < RAMLIMIT)
    {
      ESL_DPRINTF2(("#DEBUG: Solving a V:   G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		z, UniqueStatetype(cm->stid[z]),
		i0,j1,j1,j0));
      vinsideT(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, (r==0),
		NULL, NULL); /* two NULLs mean 'don't use bands' */
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
   *    beta[cm->M] is the EL deck, needed for local ends.
   */
  vinside (cm, dsq, L, w, z, i0, i1, j1, j0, useEL, BE_EFFICIENT, 
	   NULL, &alpha, NULL, &pool, NULL, (r==0), &b, &bsc);
  voutside(cm, dsq, L, r, y, i0, i1, j1, j0, useEL, BE_EFFICIENT, 
	   NULL, &beta,  pool, NULL);

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

  /* Local alignment ends: maybe we're better off in EL, not
   * the split set?
   */
  if (useEL && (cm->flags & CMH_LOCAL_END)) {
    for (ip = 0; ip <= i1-i0; ip++)
      for (jp = 0; jp <= j0-j1; jp++)
	if ((sc = beta[cm->M][jp][ip]) > best_sc) {
	  best_sc = sc;
	  best_v  = -1;
	  best_i  = ip + i0;
	  best_j  = jp + j1;
	}
  }
	
  /* Local alignment begins: maybe we're better off in root...
   */
  if (r==0 && (cm->flags & CMH_LOCAL_BEGIN)) {
    if (bsc > best_sc) {
      best_sc = bsc;
      best_v  = -2;
      best_i  = i0;
      best_j  = j0;
    }
  }

  /* Free now, before recursing!
   */
  free_vji_matrix(alpha, cm->M, j1, j0);
  free_vji_matrix(beta,  cm->M, j1, j0);

  /* If we're in EL, instead of the split set, the optimal
   * alignment is entirely in a V problem that's still above us.
   * The TRUE flag sets useEL; we propagate allow_begin. 
   */
  if (best_v == -1) {
    v_splitter(cm, dsq, L, tr, r, w, i0, best_i, best_j, j0, TRUE);    
    return;
  }

  /* If we used a local begin, the optimal alignment is
   * entirely in a V problem that's still below us, rooted
   * at b, for the entire one-hole sequence. The FALSE
   * flag prohibits more local begin transitions; we propagate
   * useEL.
   */
  if (best_v == -2) {
    if (b != z) 
      {
	InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b);
      }
    v_splitter(cm, dsq, L, tr, b, z, i0, i1, j1, j0, useEL);    
    return;
  }

  /* The optimal split into two V problems:
   *    V:   r..v, i0..i', j'..j0
   *    V:   v..z, i'..i1, j1..j'
   * Solve in this order, because we're constructing the
   * trace in postorder traversal.
   */
  ESL_DPRINTF2(("#DEBUG: V splitter:\n"));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		best_v, UniqueStatetype(cm->stid[best_v]),
		i0, best_i, best_j, j0));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		best_v, UniqueStatetype(cm->stid[best_v]),
		z, UniqueStatetype(cm->stid[z]),
		best_i, i1, j1, best_j));

  v_splitter(cm, dsq, L, tr, r,      best_v, i0,     best_i, best_j, j0, FALSE);
  v_splitter(cm, dsq, L, tr, best_v, z,      best_i, i1,     j1,     best_j, useEL);
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
 *           sequence (sq) and the matrix (alpha) in the full coordinate
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
 *                       
 *
 * Returns: Score of the optimal alignment.  
 */
static float 
inside(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
       float ***alpha, float ****ret_alpha, 
       struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
       void ****ret_shadow, 
       int allow_begin, int *ret_b, float *ret_bsc)
{
  int      status;
  float  **end;         /* we re-use the end deck. */
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

  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, vroot, E_st);
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
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * (cm->M+1));
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
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][d]  = USED_EL; 
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc; 
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
//printf("j%2d v%2d ",j,v);
//for (d = 0; d <= W && d <= j; d++) { printf("%10.2e ",alpha[v][j][d]); }
//printf("\n");
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
//printf("j%2d v%2d ",j,v);
//for (d = 0; d <= W && d <= j; d++) { printf("%10.2e ",alpha[v][j][d]); }
//printf("\n");
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
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-2] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
//printf("j%2d v%2d ",j,v);
//for (d = 0; d <= W && d <= j; d++) { printf("%10.2e ",alpha[v][j][d]); }
//printf("\n");
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
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  } 
		
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][dsq[i]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
//printf("j%2d v%2d ",j,v);
//for (d = 0; d <= W && d <= j; d++) { printf("%10.2e ",alpha[v][j][d]); }
//printf("\n");
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
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][dsq[j]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
//printf("j%2d v%2d ",j,v);
//for (d = 0; d <= W && d <= j; d++) { printf("%10.2e ",alpha[v][j][d]); }
//printf("\n");
	  }
	}				/* finished calculating deck v. */
      
      /* Check for local begin getting us to the root.
       * This is "off-shadow": if/when we trace back, we'll handle this
       * case separately (and we'll know to do it because we'll immediately
       * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
       * to jump right to state b; see below)
       */
      if (allow_begin && alpha[v][j0][W] + cm->beginsc[v] > bsc) 
	{
	  b   = v;
	  bsc = alpha[v][j0][W] + cm->beginsc[v];
	}

      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (allow_begin && v == 0 && bsc > alpha[0][j0][W]) {
	alpha[0][j0][W] = bsc;
	if (ret_shadow != NULL) yshad[j0][W] = USED_LOCAL_BEGIN;
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
  sc       = alpha[vroot][j0][W];
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

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}


/* Function: outside()
 * Date:     SRE, Tue Aug  8 10:42:52 2000 [St. Louis]
 *
 * Purpose:  Run the outside version of a CYK alignment algorithm,
 *           on a subsequence i0..j0 of a digitized sequence sq [1..L],
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
outside(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
	int do_full, float ***beta, float ****ret_beta,
	struct deckpool_s *dpool, struct deckpool_s **ret_dpool)
{
  int      status;
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

  /* if caller didn't give us a matrix, make one.
   * Allocate room for M+1 decks because we might need the EL deck (M)
   * if we're doing local alignment.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck.
   * If the root is in a split set, initialize the whole split set.
   */
  w1 = cm->nodemap[cm->ndidx[vroot]]; /* first state in split set */
  if (cm->sttype[vroot] == B_st) {    /* special boundary case of Outside on a single B state. */
    w2 = w1;
    if (vend != vroot) cm_Fail("oh no. not again.");
  } else
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

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   */
  if (cm->flags & CMH_LOCAL_END) {
    if (! deckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = alloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[cm->M][j][d] = IMPOSSIBLE;
    }
    
    /* We have to worry about vroot -> EL transitions.
     * since we start the main recursion at w2+1. This requires a 
     * laborious partial unroll of the main recursion, grabbing
     * the stuff relevant to a beta[EL] calculation for just the
     * vroot->EL transition.
     */
    if (NOT_IMPOSSIBLE(cm->endsc[vroot])) {
      switch (cm->sttype[vroot]) {
      case MP_st:
	if (W < 2) break;
	if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
	else
	  escore = DegeneratePairScore(cm->abc, cm->esc[vroot], dsq[i0], dsq[j0]);
	beta[cm->M][j0-1][W-2] = cm->endsc[vroot] + 
	  (cm->el_selfsc * (W-2)) + escore;

	if (beta[cm->M][j0-1][W-2] < IMPOSSIBLE) beta[cm->M][j0-1][W-2] = IMPOSSIBLE;
	break;
      case ML_st:
      case IL_st:
	if (W < 1) break;
	if (dsq[i0] < cm->abc->K) 
	  escore = cm->esc[vroot][(int) dsq[i0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[vroot]);
	beta[cm->M][j0][W-1] = cm->endsc[vroot] + 
	  (cm->el_selfsc * (W-1)) + escore;

	if (beta[cm->M][j0][W-1] < IMPOSSIBLE) beta[cm->M][j0][W-1] = IMPOSSIBLE;
	break;
      case MR_st:
      case IR_st:
	if (W < 1) break;
	if (dsq[j0] < cm->abc->K) 
	  escore = cm->esc[vroot][(int) dsq[j0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[vroot]);
	beta[cm->M][j0-1][W-1] = cm->endsc[vroot] + 
	  (cm->el_selfsc * (W-1)) + escore;
	
	if (beta[cm->M][j0-1][W-1] < IMPOSSIBLE) beta[cm->M][j0-1][W-1] = IMPOSSIBLE;
	break;
      case S_st:
      case D_st:
	beta[cm->M][j0][W] = cm->endsc[vroot] + 
	  (cm->el_selfsc * W);
	if (beta[cm->M][j0][W] < IMPOSSIBLE) beta[cm->M][j0][W] = IMPOSSIBLE;
	break;
      case B_st:		/* can't start w/ bifurcation at vroot. */
      default: cm_Fail("bogus parent state %d\n", cm->sttype[vroot]);
      }
    }
  }
  
  ESL_ALLOC(touch, sizeof(int) * cm->M);
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
      if (vroot == 0 && i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN))
	beta[v][j0][W] = cm->beginsc[v];

      /* main recursion:
       */
      for (jp = W; jp >= 0; jp--) {
	j = i0-1+jp;
	for (d = jp; d >= 0; d--) 
	  {
	    i = j-d+1;
	    for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	      if (y < vroot) continue; /* deal with split sets */
	      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

	      switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == j0 || d == jp) continue; /* boundary condition */

		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
		
		if ((sc = beta[y][j+1][d+2] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;

	      case ML_st:
	      case IL_st: 
		if (d == jp) continue;	/* boundary condition (note when j=0, d=0*/

		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		  
		if ((sc = beta[y][j][d+1] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		  
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);

		if ((sc = beta[y][j+1][d+1] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		if ((sc = beta[y][j][d] + cm->tsc[y][voffset]) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;

	      default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
	      }/* end switch over states*/
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;


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
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		if ((sc = beta[v][j+1][d+2] + cm->endsc[v] + 
		     (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;	
		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		if ((sc = beta[v][j][d+1] + cm->endsc[v] + 
		     (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  /*(cm->el_selfsc * (d+1)) + escore) > beta[cm->M][j][d])*/
		  beta[cm->M][j][d] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		if ((sc = beta[v][j+1][d+1] + cm->endsc[v] + 
		     (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		     /*(cm->el_selfsc * (d+1)) + escore) > beta[cm->M][j][d])*/
		  beta[cm->M][j][d] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if ((sc = beta[v][j][d] + cm->endsc[v] +
		     (cm->el_selfsc * d)) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case B_st:  
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
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

#if 0
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      for (d = jp-1; d >= 0; d--) /* careful w/ boundary here */
	if ((sc = beta[cm->M][j][d+1]) > beta[cm->M][j][d])
	  beta[cm->M][j][d] = sc;
    }
  }
#endif

  /* If the caller doesn't want the matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = w1; v <= vend; v++) /* start at w1 - top of split set - not vroot */
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
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
  return;

 ERROR:
  cm_Fail("Memory allocation error.\n");
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
 *           useEL     - if TRUE, V problem ends at EL/i1/j1, not z/i1/j1
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
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        

 * 
 * Returns:  score.
 */
static float
vinside(CM_t *cm, ESL_DSQ *dsq, int L, 
	int r, int z, int i0, int i1, int j1, int j0, int useEL,
	int do_full, float ***a, float ****ret_a,
	struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	char ****ret_shadow,
	int allow_begin, int *ret_b, float *ret_bsc)
{
  int      status;
  char  ***shadow;              /* the shadow matrix -- traceback ptrs -- memory is kept */
  int     v,i,j;
  int     w1,w2;		/* bounds of the split set */
  int     jp, ip;		/* j' and i' -- in the matrix coords */
  int    *touch;                /* keeps track of whether we can free a deck yet or not */
  int     y, yoffset;
  float   sc;			/* tmp variable holding a score */
  int      b;			/* best local begin state */
  float    bsc;			/* score for using the best local begin state */

  /*printf("***in vinside()****\n");
    printf("\tr  : %d\n", r);
    printf("\tz  : %d\n", z);
    printf("\ti0 : %d\n", i0);
    printf("\ti1 : %d\n", i1);
    printf("\tj1 : %d\n", j1);
    printf("\tj0 : %d\n", j0);
  */

  /* Allocations, initializations.
   * Remember to allocate for M+1 decks, in case we reuse this
   * memory for a local alignment voutside() calculation.
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  if (cyk_dnc_track) cyk_dnc_vji_row_floats = i1 - i0 + 1; /* brief 26_0610-008 vji accounting */
  if (dpool == NULL) dpool = deckpool_create();
  if (a == NULL) {
    ESL_ALLOC(a, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) a[v] = NULL;
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

  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(char **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL; 
  }

  /* Initialize the one non-IMPOSSIBLE cell as a boundary
   * condition.
   * If local alignment (useEL=1), we must connect z to EL;
   * we would init a[EL][0][i1-i0] = 0. But, we're not explicitly
   * keeping an EL deck, we're swallowing it into the recursion.
   * So, we unroll a chunk of the main recursion;
   * we have to laboriously figure out from the statetype z 
   * and our position where and what our initialization is.
   * Else, for global alignments, we simply connect to z,0,i1-i0.
   */
  ip = i1-i0;
  jp = 0;
  if (! useEL) 
    a[z][jp][ip] = 0.;
  else 
    {
      if (ret_shadow != NULL) 
	shadow[z] = alloc_vji_shadow_deck(i0,i1,j1,j0); 

      switch (cm->sttype[z]) {
      case D_st:
      case S_st:
	/*a[z][jp][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));
	if (ret_shadow != NULL) shadow[z][jp][ip] = USED_EL;
	break;
      case MP_st:
	if (i0 == i1 || j1 == j0) break;
	/*a[z][jp+1][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp+1][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));
	if (dsq[i1-1] < cm->abc->K && dsq[j1+1] < cm->abc->K)
	  a[z][jp+1][ip-1] += cm->esc[z][(int) (dsq[i1-1]*cm->abc->K+dsq[j1+1])];
	else
	  a[z][jp+1][ip-1] += DegeneratePairScore(cm->abc, cm->esc[z], dsq[i1-1], dsq[j1+1]);
	if (ret_shadow != NULL) shadow[z][jp+1][ip-1] = USED_EL;
	if (a[z][jp+1][ip-1] < IMPOSSIBLE) a[z][jp+1][ip-1] = IMPOSSIBLE;
	break;
      case ML_st:
      case IL_st:
	if (i0==i1) break;
	/*a[z][jp][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));
	if (dsq[i1-1] < cm->abc->K)
	  a[z][jp][ip-1] += cm->esc[z][(int) dsq[i1-1]];
	else
	  a[z][jp][ip-1] += esl_abc_FAvgScore(cm->abc, dsq[i1-1], cm->esc[z]);
	if (ret_shadow != NULL) shadow[z][jp][ip-1] = USED_EL;
	if (a[z][jp][ip-1] < IMPOSSIBLE) a[z][jp][ip-1] = IMPOSSIBLE;
	break;
      case MR_st:
      case IR_st:
	if (j1==j0) break;
	/*a[z][jp+1][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp+1][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));
	if (dsq[j1+1] < cm->abc->K)
	  a[z][jp+1][ip] += cm->esc[z][(int) dsq[j1+1]];
	else
	  a[z][jp+1][ip] += esl_abc_FAvgScore(cm->abc, dsq[j1+1], cm->esc[z]);
	if (ret_shadow != NULL) shadow[z][jp+1][ip] = USED_EL;
	if (a[z][jp+1][ip] < IMPOSSIBLE) a[z][jp+1][ip] = IMPOSSIBLE;
	break;
      }

    } /* done initializing the appropriate cell for useEL=TRUE */

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;   v < r;  v++) touch[v] = 0;
  for (v = r;   v <= w2; v++) touch[v] = cm->pnum[v]; /* note w2 not z: to bottom of split set */
  for (v = w2+1; v < cm->M; v++) touch[v] = 0;

  /* A special case. If vinside() is called on empty sequences,
   * we might do a begin transition right into z.
   */ 
  if (allow_begin && j0-j1 == 0 && i1-i0 == 0)
    {
      b   = z;
      bsc = a[z][0][0] + cm->beginsc[z];
      if (z == 0) { 
	a[0][0][0] = bsc;
	if (ret_shadow != NULL) shadow[0][0][0] = USED_LOCAL_BEGIN;
      }
    }

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
	cm_Fail("you told me you wouldn't ever do that again.");
      
      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= j0-j1; jp++) 
	    for (ip = i1-i0; ip >= 0; ip--) {
	      /*printf("D S jp : %d | ip : %d\n", jp, ip);*/
	      y = cm->cfirst[v];
	      a[v][jp][ip]      = a[y][jp][ip] + cm->tsc[v][0];
	      /*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
	      if (ret_shadow != NULL) shadow[v][jp][ip] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		  > a[v][jp][ip])) {
		a[v][jp][ip]      = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc;
		    /*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
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
	      /*printf("MP jp : %d | ip : %d\n", jp, ip);*/
	      i = ip+i0;
	      y = cm->cfirst[v];
	      a[v][jp][ip] = a[y][jp-1][ip+1] + cm->tsc[v][0];
	      /*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
	      if (ret_shadow != NULL) shadow[v][jp][ip] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		  > a[v][jp][ip])) {
		a[v][jp][ip]      = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		   { 
		     a[v][jp][ip] = sc; 
		     /*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
		     if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		   }
	      if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		a[v][jp][ip] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
	      else
		a[v][jp][ip] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
	  
	  for (jp = 0; jp <= j0-j1; jp++) { 
	    a[v][jp][i1-i0] = IMPOSSIBLE; /* boundary condition */
	    for (ip = i1-i0-1; ip >= 0; ip--) {
	      /*printf("ML IL jp : %d | ip : %d\n", jp, ip);*/
	      i = ip+i0;
	      y = cm->cfirst[v];
	      a[v][jp][ip] = a[y][jp][ip+1] + cm->tsc[v][0];
	      if (ret_shadow != NULL) shadow[v][jp][ip] = 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		  > a[v][jp][ip])) {
		a[v][jp][ip]      = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		/*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc; 
		    /*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
	      
	      if (dsq[i] < cm->abc->K)
		a[v][jp][ip] += cm->esc[v][dsq[i]];
	      else
		a[v][jp][ip] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
	  for (ip = i1-i0; ip >= 0; ip--) a[v][0][ip] = IMPOSSIBLE; /* boundary condition */

	  for (jp = 1; jp <= j0-j1; jp++) { 
	    j = jp+j1;
	    for (ip = i1-i0; ip >= 0; ip--) {
	      /*printf("MR IR jp : %d | ip : %d\n", jp, ip);*/
	      y = cm->cfirst[v];
	      a[v][jp][ip]      = a[y][jp-1][ip] + cm->tsc[v][0];
	      /*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
	      if (ret_shadow != NULL) shadow[v][jp][ip] = 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		  > a[v][jp][ip])) {
		a[v][jp][ip] = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc; 
		    /*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
	      
	      if (dsq[j] < cm->abc->K)
		a[v][jp][ip] += cm->esc[v][dsq[j]];
	      else
		a[v][jp][ip] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} /* finished calculating deck v */

      /* Check for local begin getting us to the root.
       */
      if (allow_begin && a[v][j0-j1][0] + cm->beginsc[v] > bsc) 
	{
	  b   = v;
	  bsc = a[v][j0-j1][0] + cm->beginsc[v];
	}

      /* Check whether we need to store the local begin score
       * for a possible traceback.
       */
      if (allow_begin && v == 0 && bsc > a[0][j0-j1][0]) 
	{
	  a[0][j0-j1][0] = bsc;
	  if (ret_shadow != NULL) shadow[v][j0-j1][0] = USED_LOCAL_BEGIN;
	}


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
  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */


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

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
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
 *           analogous to inside() and vinside(); see their prefaces
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
 *           useEL     - if TRUE, worry about local alignment.
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
voutside(CM_t *cm, ESL_DSQ *dsq, int L, 
	 int r, int z, int i0, int i1, int j1, int j0, int useEL,
	 int do_full, float ***beta, float ****ret_beta,
	 struct deckpool_s *dpool, struct deckpool_s **ret_dpool)
{
  int      status;
  int      v,y;			/* indices for states */
  int      i,j;			/* indices in sequence dimensions */
  int      ip, jp;		/* transformed sequence indices */
  float    sc;			/* a temporary variable holding a score */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */
  int      voffset;		/* index of v in t_v(y) transition scores */


  /* Allocations and initializations
   */
  if (cyk_dnc_track) cyk_dnc_vji_row_floats = i1 - i0 + 1; /* brief 26_0610-008 vji accounting */
  			/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* If caller didn't give us a matrix, make one.
   * Remember to allow for deck M, the EL deck, for local alignments.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) beta[v] = NULL;
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

  /* Initialize the EL deck, if we're in local mode w.r.t. ends.
   * Deal with the special initialization case of the root state r
   * immediately transitioning to EL, if we're supposed to use EL.
   */
  if (useEL && cm->flags & CMH_LOCAL_END) {
    if (! deckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = alloc_vji_deck(i0,i1,j1,j0);
    for (jp = 0; jp <= j0-j1; jp++) {
      for (ip = 0; ip <= i1-i0; ip++)
	beta[cm->M][jp][ip] = IMPOSSIBLE;
    }
  }
  if (useEL && NOT_IMPOSSIBLE(cm->endsc[r])) {
    switch(cm->sttype[r]) {
    case MP_st:
      if (i0 == i1 || j1 == j0) break;
      if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	escore = cm->esc[r][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
      else
	escore = DegeneratePairScore(cm->abc, cm->esc[r], dsq[i0], dsq[j0]);
      beta[cm->M][j0-j1-1][1] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0-1)-(i0+1)+1)) + escore;
      break;
    case ML_st:
    case IL_st:
      if (i0 == i1) break;
      if (dsq[i0] < cm->abc->K) 
	escore = cm->esc[r][(int) dsq[i0]];
      else
	escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[r]);
      beta[cm->M][j0-j1][1] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0)-(i0+1)+1)) + escore;
      break;
    case MR_st:
    case IR_st:
      if (j0==j1) break;
      if (dsq[j0] < cm->abc->K) 
	escore = cm->esc[r][(int) dsq[j0]];
      else
	escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[r]);
      beta[cm->M][j0-j1-1][0] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0-1)-(i0)+1)) + escore;
      break;
    case S_st:
    case D_st:
      beta[cm->M][j0-j1][0] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0)-(i0)+1));
      break;
    default:  cm_Fail("bogus parent state %d\n", cm->sttype[r]);
    }
  }
      
  /* Initialize the "touch" array, used for figuring out
   * when a deck is no longer touched, so it can be free'd.
   */
  ESL_ALLOC(touch, sizeof(int) * cm->M);
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

      /* Init the whole deck to IMPOSSIBLE.
       */
      for (jp = j0-j1; jp >= 0; jp--) 
	for (ip = 0; ip <= i1-i0; ip++) 
	  beta[v][jp][ip] = IMPOSSIBLE;

      /* If we can get into deck v by a local begin transition, do an init
       * with that.
       */
      if (r == 0 && i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN))
	{
	  if (cm->beginsc[v] > beta[v][j0-j1][0]) 
	    beta[v][j0-j1][0] = cm->beginsc[v];
	}

      /* main recursion:
       */
      for (jp = j0-j1; jp >= 0; jp--) {
	j = jp+j1;
	for (ip = 0; ip <= i1-i0; ip++) 
	  {
	    i = ip+i0;

	    for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	      if (y < r) continue; /* deal with split sets */
	      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

	      switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == j0 || i == i0) continue; /* boundary condition */
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
		
		if ((sc = beta[y][jp+1][ip-1]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;

	      case ML_st:
	      case IL_st: 
		if (i == i0) continue;	/* boundary condition */

		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		  
		if ((sc = beta[y][jp][ip-1]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		  
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);

		if ((sc = beta[y][jp+1][ip]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		if ((sc = beta[y][jp][ip] + cm->tsc[y][voffset]) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;

	      default: cm_Fail("bogus parent state %d\n", cm->sttype[y]);
	      }/* end switch over states*/
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (beta[v][jp][ip] < IMPOSSIBLE) beta[v][jp][ip] = IMPOSSIBLE;

	  } /* ends loop over ip. We know all beta[v][jp][ip] in this row jp */

      }/* end loop over jp. We know the beta's for the whole deck.*/

      /* Deal with local alignment
       * transitions v->EL, if we're doing local alignment and there's a 
       * possible transition.
       */
      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (jp = j0-j1; jp >= 0; jp--) {
	  j = jp+j1;
	  for (ip = 0; ip <= i1-i0; ip++) 
	    {
	      i = ip+i0;
	      switch (cm->sttype[v]) {
	      case MP_st:
		if (j == j0 || i == i0) continue; /* boundary condition */
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		if ((sc = beta[v][jp+1][ip-1] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		  beta[cm->M][jp][ip] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (i == i0) continue;
		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		if ((sc = beta[v][jp][ip-1] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		  beta[cm->M][jp][ip] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		if ((sc = beta[v][jp+1][ip] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		  beta[cm->M][jp][ip] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if ((sc = beta[v][jp][ip] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1))) > beta[cm->M][jp][ip])
		    beta[cm->M][jp][ip] = sc;
		break;
	      default:  cm_Fail("bogus parent state %d\n", cm->sttype[y]);
	      } /* end switch over parent v state type */
	    } /* end loop over ip */
	} /* end loop over jp */
      }
	
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

#if 0 
  /* superfluous code, I think...*/
  /* Deal with the last step needed for local alignment
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   */
  if (useEL && cm->flags & CMH_LOCAL_END) {
    for (jp = j0-j1; jp >= 0; jp--) 
      for (ip = 1; ip <= i1-i0; ip++) /* careful w/ boundary here */
	if ((sc = beta[cm->M][jp][ip-1]) > beta[cm->M][jp][ip]) 
	  beta[cm->M][jp][ip] = sc;
  }
#endif

  /* If the caller doesn't want the matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation!)
   */
  if (ret_beta == NULL) {
    for (v = r; v <= z; v++)
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
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
    while (deckpool_pop(dpool, &a)) 
      free_vji_deck(a,j1,j0);
    deckpool_free(dpool);
  } else *ret_dpool = dpool;

  free(touch);
  return;

 ERROR:
  cm_Fail("Memory allocation error.\n");
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
 *
 *           If we're not in banded mode, dmin and dmax should
 *           be passed in as NULL.
 */
static float
insideT(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
	int r, int z, int i0, int j0, 
	int allow_begin, int *dmin, int *dmax)
{

  int       status;
  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;

  if(dmin == NULL && dmax == NULL)
    {
      sc = inside(cm, dsq, L, r, z, i0, j0, 
		  BE_EFFICIENT,	/* memory-saving mode */
		  NULL, NULL,	/* manage your own matrix, I don't want it */
		  NULL, NULL,	/* manage your own deckpool, I don't want it */
		  &shadow,	/* return a shadow matrix to me. */
		  allow_begin,  /* TRUE to allow local begins */
		  &b, &bsc);	/* if allow_begin is TRUE, gives info on optimal b */
    }
  else
    {
      sc = inside_qdb(cm, dsq, L, r, z, i0, j0, 
		    BE_EFFICIENT,/* memory-saving mode */
		    NULL, NULL,	 /* manage your own matrix, I don't want it */
		    NULL, NULL,	 /* manage your own deckpool, I don't want it */
		    &shadow,	 /* return a shadow matrix to me. */
		    allow_begin, /* TRUE to allow local begins */
		    &b, &bsc,	 /* if allow_begin is TRUE, gives info on optimal b */
		    dmin, dmax); /* the bands */
    }      
  
  pda = esl_stack_ICreate();
  if(pda == NULL) goto ERROR;
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  /*printf("Starting traceback in insideT()\n");*/
  while (1) {
    if (cm->sttype[v] == B_st) {
      k = ((int **) shadow[v])[j][d];   /* k = len of right fragment */

      /* Store info about the right fragment that we'll retrieve later:
       */
      if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;	/* remember the end j    */
      if((status = esl_stack_IPush(pda, k)) != eslOK) goto ERROR;	/* remember the subseq length k */
      if((status = esl_stack_IPush(pda, tr->n-1)) != eslOK) goto ERROR;	/* remember the trace index of the parent B state */

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
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
    } else {
      yoffset = ((char **) shadow[v])[j][d];

      /*printf("v : %d | r : %d | z : %d | i0 : %d | \n", v, r, z, i0);*/
      /*printf("\tyoffset : %d\n", yoffset);*/
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
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

 ERROR: 
  cm_Fail("Memory allocation error.");
  return 0.; /* NEVERREACHED */
}

/* Function: vinsideT()
 * Date:     SRE, Sat Jun  2 14:40:13 2001 [St. Louis]
 *
 * Purpose:  Call vinside(), get vji shadow matrix for a V problem;
 *           then trace back. Append the trace to a
 *           given traceback, which has state r already at
 *           t->n-1.
 *
 *           If we're not in banded mode, dmin and dmax should
 *           be passed in as NULL.
 */
static float
vinsideT(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
	 int r, int z, int i0, int i1, int j1, int j0, int useEL, 
	 int allow_begin, int *dmin, int *dmax)
{
  char ***shadow;
  float   sc;
  int     v,y;
  int     j,i;
  int     jp,ip;
  int     yoffset;
  int     b;
  float   bsc;

  /* If we can deduce the traceback unambiguously without
   * doing any DP... do it.
   */
  if (r == z) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, r);
    return 0.;
  }

  if(dmin == NULL && dmax == NULL)
    {
      sc = vinside(cm, dsq, L, r, z, i0, i1, j1, j0, useEL,
		   BE_EFFICIENT,	/* memory-saving mode */
		   NULL, NULL,	/* manage your own matrix, I don't want it */
		   NULL, NULL,	/* manage your own deckpool, I don't want it */
		   &shadow,      	/* return a shadow matrix to me. */
		   allow_begin,     /* TRUE to allow local begin transitions */
		   &b, &bsc);       /* info on optimal local begin */
    }
  else
    {
      sc = vinside_qdb(cm, dsq, L, r, z, i0, i1, j1, j0, useEL,
		     BE_EFFICIENT,	/* memory-saving mode */
		     NULL, NULL,	/* manage your own matrix, I don't want it */
		     NULL, NULL,	/* manage your own deckpool, I don't want it */
		     &shadow,      	/* return a shadow matrix to me. */
		     allow_begin,       /* TRUE to allow local begin transitions */
		     &b, &bsc,          /* info on optimal local begin */
		     dmin, dmax);
    }
  /* We've got a complete shadow matrix. Trace it back. We know
   * that the trace will begin with the start state r, at i0,j0
   * (e.g. jp=j0-j1, ip=0)
   */
  v = r;
  j = j0;
  i = i0;

  /*printf("Starting traceback in vinsideT()\n");*/
  while (1) {
    jp = j-j1;
    ip = i-i0;

    /* 1. figure out the next state (deck) in the shadow matrix.
     */ 
    /*printf("v : %d | jp : %d | ip : %d | i0 : %d | \n", v, jp, ip, i0);*/
    yoffset = shadow[v][jp][ip];
    /*printf("\tyoffset : %d\n", yoffset);*/

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
    default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
    }

    /* If the traceback pointer (yoffset) is -1, that's a special
     * flag for a local alignment end, e.g. transition to EL (state "M").
     */
    if (yoffset == USED_EL) 
      {
	InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	break;			/* one way out of the while loop */
      }
    else if (yoffset == USED_LOCAL_BEGIN) 
      {
	InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
	v = b;
	if (! useEL && v == z) break; /* the other way out of the while loop */
      }
    else
      {
	/*    Attach y,i,j to the trace. This new node always attaches
	 *    to the end of the growing trace -- e.g. trace node
	 *    tr->n-1.
	 */
	y = cm->cfirst[v] + yoffset;
	InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
	v = y;
	if (! useEL && v == z) break; /* the other way out of the while loop */
      }
  }
  
  /* We're done. Our traceback has just ended. We have just attached
   * state z for i1,j1; it is in the traceback at node tr->n-1.
   */
  free_vji_shadow_matrix(shadow, cm->M, j1, j0);
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
float
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

float
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
 *           
 *           For a whole model, except for trivially small models with no
 *           stacked base pairs, this is almost invariably 
 *           10+1+cyk_extra_decks(): MATP-MATP connections require
 *           10 decks (6 states in current node, 4 states in connected
 *           split set of next node). We share 1 end state deck. All
 *           other decks are retained S decks, needed for bifurcation
 *           calculations.  
 */
static int
cyk_deck_count(CM_t *cm, int r, int z)
{
  int       status;
  ESL_STACK *pda;	/* pushdown stack simulating the deck pool */
  int       v,w,y;	/* state indices */
  int       nends;
  int       ndecks;
  int      *touch;	/* keeps track of how many higher decks still need this deck */

  /* Initializations, mirroring key parts of CYKInside()
   */
  ndecks = 1;			/* deck z, which we always need to start with. */
  nends  = CMSegmentCountStatetype(cm, r, z, E_st);
  pda    = esl_stack_ICreate();
  if(pda == NULL) goto ERROR;

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v < r;     v++) touch[v] = 0;
  for (v = r; v < z;     v++) touch[v] = cm->pnum[v];
  for (v = z; v < cm->M; v++) touch[v] = 0;

  for (v = z; v >= r; v--)
    {
      if (cm->sttype[v] != E_st) {
	if (esl_stack_IPop(pda, &y) == eslEOD) ndecks++; /* simulated allocation of a new deck */
      }
      
      if (cm->sttype[v] == B_st) { /* release both S children of a bifurc */
	w = cm->cfirst[v];
	y = cm->cnum[v];
	if((status =esl_stack_IPush(pda, w)) != eslOK) goto ERROR;
	if((status = esl_stack_IPush(pda, y)) != eslOK) goto ERROR;
      } else {
	for (w = cm->cfirst[v]; w < cm->cfirst[v]+cm->cnum[v]; w++)
	  {
	    touch[w]--;
	    if (touch[w] == 0) 
	      {
		if (cm->sttype[w] == E_st) { 
		  nends--; 
		  if (nends == 0) { if((status = esl_stack_IPush(pda, cm->M-1)) != eslOK) goto ERROR; }
		} else 
		  if((status = esl_stack_IPush(pda, w)) != eslOK) goto ERROR;
	      }
	  }
      }
    }
  free(touch);
  esl_stack_Destroy(pda);
  return ndecks;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0; /* never reached */
}

/* Function: cyk_extra_decks()
 * Date:     SRE, Sun Apr  7 14:42:48 2002 [St. Louis]
 *
 * Purpose:  Calculate the number of extra
 *           decks that will be needed to accommodate bifurc
 *           calculations.
 *
 * Args:     cm - the model.
 *
 * Returns:  # of extra decks.
 */
static int
cyk_extra_decks(CM_t *cm)
{
  int  max;
  int  x;
  int  v;

  max = x = 0;
  for (v = cm->M-1; v >= 0; v--) 
    {
      if      (cm->sttype[v] == S_st) x++;
      else if (cm->sttype[v] == B_st) x-=2;
      if (x > max) max = x;
    }
  return max-1;			/* discount ROOT S */
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
struct deckpool_s *
deckpool_create(void)
{
  int status;
  struct deckpool_s *dpool;

  ESL_ALLOC(dpool, sizeof(struct deckpool_s));
  dpool->block  = 10;		/* configurable if you want */
  ESL_ALLOC(dpool->pool, sizeof(float **) * dpool->block);
  dpool->nalloc = dpool->block;;
  dpool->n      = 0;
  return dpool;
 ERROR:
  cm_Fail("Memory allocation error.\n");
  return NULL; /* never reached */
}
void 
deckpool_push(struct deckpool_s *dpool, float **deck)
{
  int   status;
  void *tmp;
  if (dpool->n == dpool->nalloc) {
    dpool->nalloc += dpool->block;
    ESL_RALLOC(dpool->pool, tmp, sizeof(float **) * dpool->nalloc);
  }
  dpool->pool[dpool->n] = deck;
  dpool->n++;
  ESL_DPRINTF3(("#DEBUG: deckpool_push\n"));
  return;
 ERROR:
  cm_Fail("Memory reallocation error.\n");
}
int
deckpool_pop(struct deckpool_s *d, float ***ret_deck)
{
  if (d->n == 0) { *ret_deck = NULL; return 0;}
  d->n--;
  *ret_deck = d->pool[d->n];
  ESL_DPRINTF3(("#DEBUG: deckpool_pop\n"));
  return 1;
}
void
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
/* data (per-row float arrays) bytes of a vjd deck spanning rows i-1..j; the
 * O(W^2) bulk that dominates the working set. Excludes the (L+1)-pointer array
 * (small, and free_vjd_deck() doesn't know L). Used by the D&C-HB high-water
 * instrumentation (cyk_dnc_track). */
double
vjd_deck_data_bytes(int i, int j)
{
  long W = (long) (j - i + 1);                 /* max d for this deck */
  return (double) sizeof(float) * (double) ((W+1)*(W+2)/2);
}

float **
alloc_vjd_deck(int L, int i, int j)
{
  int status;
  float **a;
  int     jp;
  ESL_DPRINTF3(("#DEBUG: alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  ESL_ALLOC(a, sizeof(float *) * (L+1)); /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0;   jp < i-1;    jp++) a[jp]     = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp]     = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) ESL_ALLOC(a[jp+i-1], sizeof(float) * (jp+1));
  if (cyk_dnc_track) {
    cyk_dnc_cur_bytes += vjd_deck_data_bytes(i, j);
    cyk_dnc_note();
  }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
float
size_vjd_deck(int L, int i, int j)
{
  float Mb;
  int   jp;
  Mb = (float) (sizeof(float *) * (L+1));
  for (jp = 0; jp <= j-i+1; jp++)
    Mb += (float) (sizeof(float) * (jp+1));
  return (Mb / 1000000.);
}
void
free_vjd_deck(float **a, int i, int j)
{
  int jp;
  if (cyk_dnc_track) cyk_dnc_cur_bytes -= vjd_deck_data_bytes(i, j);
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp+i-1]);
  free(a);
}
void
free_vjd_matrix(float ***a, int M, int i, int j)
{
  int v;
  for (v = 0; v <= M; v++)
    if (a[v] != NULL)		/* protect against double free's of reused decks (ends) */
      { free_vjd_deck(a[v], i, j); a[v] = NULL; }
  free(a);
}
char **
alloc_vjd_yshadow_deck(int L, int i, int j)
{
  int status;
  char **a;
  int    jp;
  ESL_ALLOC(a, sizeof(char *) * (L+1)); /* always alloc 0..L rows, same as alloc_deck */
  for (jp = 0;   jp < i-1;    jp++) a[jp] = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp] = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) ESL_ALLOC(a[jp+i-1], sizeof(char) * (jp+1));
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
float
size_vjd_yshadow_deck(int L, int i, int j)
{
  float  Mb;
  int    jp;
  Mb = (float) (sizeof(char *) * (L+1));
  for (jp = 0; jp <= j-i+1; jp++) 
    Mb += (float) (sizeof(char) * (jp+1));
  return Mb / 1000000.;
}
void
free_vjd_yshadow_deck(char **a, int i, int j)
{
  int jp;
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp+i-1]);
  free(a);
}
int **
alloc_vjd_kshadow_deck(int L, int i, int j)
{
  int status;
  int **a;
  int   jp;
  ESL_ALLOC(a, sizeof(int *) * (L+1)); /* always alloc 0..L rows, same as alloc_deck */
  for (jp = 0;   jp <  i-1;   jp++) a[jp] = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) ESL_ALLOC(a[jp+i-1], sizeof(int) * (jp+1));
  for (jp = j+1; jp <= L;     jp++) a[jp] = NULL;
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
float
size_vjd_kshadow_deck(int L, int i, int j)
{
  float Mb;
  int   jp;
  
  Mb = (float)(sizeof(int *) * (L+1)); 
  for (jp = 0;   jp <= j-i+1; jp++)
    Mb += (float) (sizeof(int) * (jp+1));
  return Mb / 1000000.;
}
void
free_vjd_kshadow_deck(int **a, int i, int j)
{
  int jp;
  /*11.14.05 old line: for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp]);*/
  for (jp = 0; jp <= j-i+1; jp++) if (a[jp+i-1] != NULL) free(a[jp-i+1]);
  free(a);
}
void
free_vjd_shadow_matrix(void ***shadow, CM_t *cm, int i, int j)
{
  int v;
  for (v = 0; v < cm->M; v++)
    if (shadow[v] != NULL) {
      if (cm->sttype[v] == B_st) free_vjd_kshadow_deck((int **)  shadow[v], i, j);
      else                       free_vjd_yshadow_deck((char **) shadow[v], i, j);
      shadow[v] = NULL;
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
float **                 /* allocation of a score deck. */
alloc_vji_deck(int i0, int i1, int j1, int j0)
{
  int status; 
  float **a;
  int     jp;
  ESL_DPRINTF3(("#DEBUG: alloc_vji_deck : %.4f\n", size_vji_deck(i0,i1,j1,j0)));
  ESL_ALLOC(a, sizeof(float *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    ESL_ALLOC(a[jp], sizeof(float)*(i1-i0+1));
  if (cyk_dnc_track) {  /* brief 26_0610-008: count the (full, unbanded) V-problem decks */
    cyk_dnc_vji_bytes += (double) sizeof(float) * (double)(j0-j1+1) * (double)(i1-i0+1);
    cyk_dnc_note();
  }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
float
size_vji_deck(int i0, int i1, int j1, int j0)
{
  float Mb;
  int   jp;
  Mb = (float)(sizeof(float *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    Mb += (float)(sizeof(float)*(i1-i0+1));
  return Mb / 1000000.;
}
void			/* free'ing a score deck */
free_vji_deck(float **a, int j1, int j0)
{
  int jp;
  ESL_DPRINTF3(("#DEBUG: free_vji_deck called\n"));
  if (cyk_dnc_track)  /* brief 26_0610-008: width (i1-i0+1) stashed by vinside/voutside entry */
    cyk_dnc_vji_bytes -= (double) sizeof(float) * (double)(j0-j1+1) * (double) cyk_dnc_vji_row_floats;
  for (jp = 0; jp <= j0-j1; jp++) 
    if (a[jp] != NULL) free(a[jp]);
  free(a);
}
void
free_vji_matrix(float ***a, int M, int j1, int j0)
{
  int v;
  /* Free the whole matrix - even if we used only a subset of
   * the decks, all initialization routines init all decks 0..M
   * to NULL, so this is safe. (see bug #i2).
   */                         
  for (v = 0; v <= M; v++) 
    if (a[v] != NULL) { free_vji_deck(a[v], j1, j0); a[v] = NULL; }
  free(a);
}
char **		        /* allocation of a traceback ptr (shadow matrix) deck */
alloc_vji_shadow_deck(int i0, int i1, int j1, int j0)
{
  int status; 
  char **a;
  int     jp;
  ESL_ALLOC(a, sizeof(char *) * (j0-j1+1)); 
  for (jp = 0; jp <= j0-j1; jp++)
    ESL_ALLOC(a[jp], sizeof(char)*(i1-i0+1));
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
float		        /* allocation of a traceback ptr (shadow matrix) deck */
size_vji_shadow_deck(int i0, int i1, int j1, int j0)
{
  float   Mb;
  int     jp;
  Mb = (float)(sizeof(char *) * (j0-j1+1));
  for (jp = 0; jp <= j0-j1; jp++)
    Mb += (float)(sizeof(char)*(i1-i0+1));
  return Mb / 1000000;
}
void	                /* free'ing a shadow deck */
free_vji_shadow_deck(char **a, int j1, int j0)
{
  int jp;
  for (jp = 0; jp <= j0-j1; jp++) 
    if (a[jp] != NULL) free(a[jp]);
  free(a);
}
void
free_vji_shadow_matrix(char ***a, int M, int j1, int j0)
{
  int v;
  for (v = 0; v < M; v++) 
    if (a[v] != NULL) { free_vji_shadow_deck(a[v], j1, j0); a[v] = NULL; }
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
CYKOutside(CM_t *cm, ESL_DSQ *dsq, int L, float ***alpha)
{
  int      status;
  float ***beta;		/* the scoring cube [v=0..M-1][j=0..L][d=0..j]*/
  int      v,y,z;		/* indices for states */
  int      j,d,i,k;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  struct deckpool_s *dpool;     /* a pool of decks for beta that we can reuse */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */

  /* Allocations and initializations
   */
  ESL_ALLOC(beta, (sizeof(float **) * cm->M));
  for (v = 0; v < cm->M; v++) beta[v] = NULL;

  dpool = deckpool_create();

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; j++)
      beta[0][j][d] = IMPOSSIBLE; /* can prob speed this initialization up */
  beta[0][L][L] = 0;		
  
  /* Main loop down through the decks
   */
  /* EPN bug fix 05.25.06. Durbin et. al. p.287 CM Outside alg uses state
   * indices 1..M, with state 1 = ROOT_S, so there's an off-by-one 
   * w.r.t this implementation. Following loop followed Durbin convention,
   * but should follow implemented convention:
   * OLD LINE: for (v = 2; v < cm->M; v++)
   */
  for (v = 1; v < cm->M; v++)  
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

		    if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		      escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		    else
		      escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);

		    if ((sc = beta[y][j+1][d+2] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;

		  case ML_st:
		  case IL_st: 
		    if (d == j) continue;	/* boundary condition (note when j=0, d=0*/

		    if (dsq[i-1] < cm->abc->K) 
		      escore = cm->esc[y][(int) dsq[i-1]];
		    else
		      escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		  
		    if ((sc = beta[y][j][d+1] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;
		  
		  case MR_st:
		  case IR_st:
		    if (d == j || j == L) continue;
		  
		    if (dsq[j+1] < cm->abc->K) 
		      escore = cm->esc[y][(int) dsq[j+1]];
		    else
		      escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);

		    if ((sc = beta[y][j+1][d+1] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;
		  
		  case B_st:
		  case E_st:
		  case D_st:
		    if ((sc = beta[y][j][d] + cm->tsc[y][v]) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;

		  default: cm_Fail("bogus parent state %d\n", cm->sttype[y]);
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
  return;
 ERROR:
  cm_Fail("Memory allocation error.");
}
#endif 

/*################################################################
 * The banded dividers and conquerors. 
 *################################################################*/  

/* Function: generic_splitter_qdb()
 *           EPN 05.19.05
 * *based on generic_splitter(), only difference is bands are used : 
 * Date:     SRE, Sat May 12 15:08:38 2001 [CSHL]
 *
 * Purpose:  Solve a "generic problem": best parse of
 *           a possibly bifurcated subgraph cm^r_z to
 *           a substring dsq[i0..j0]. r is usually a start
 *           state (S_st) but may be any non-end state type in 
 *           the case of local alignment begins (ROOT 0->r).
 *           z is always an end state (E_st).
 *
 *           Given: a cm subgraph from r..z
 *                  a subsequence from i0..j0
 *           Attaches the optimal trace T{r..z}, exclusive of r
 *           and inclusive of z, to tr.
 *           
 *           A full divide & conquer never terminates
 *           in generic_splitter; the recursion must
 *           terminate in v_splitter and wedge_splitter;
 *           so we don't test an end-of-recursion boundary.
 *           
 * Args:     cm          - model
 *           sq          - digitized sequence 1..L
 *           tr          - the traceback we're adding on to.
 *           r           - index of the root state of this problem in the model       
 *           z           - index of an end state (E_st) in the model
 *           i0          - start in the sequence (1..L)
 *           j0          - end in the sequence (1..L)
 *           dmin   - minimum d bound for each state v; [0..v..M-1]
 *           dmax   - maximum d bound for each state v; [0..v..M-1]
 *
 * Returns:  score of the optimal parse of dsq(i0..j0) with cm^r_z 
 */
static float
generic_splitter_qdb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
		 int r, int z, int i0, int j0, int *dmin, int *dmax)
{
  float ***alpha;
  float ***beta;
  struct deckpool_s *pool;
  int      v,w,y;		/* state indices */
  int      wend, yend;		/* indices for end of subgraphs rooted at w,y */
  int      jp;			/* j': relative position in subseq, 0..W */
  int      W;			/* length of subseq i0..j0 */
  float    sc;			/* tmp variable for a score */
  int      j,d,k;		/* sequence indices */
  float    best_sc;		/* optimal score at the optimal split point */
  int      best_k;		/* optimal k for the optimal split */
  int      best_d;		/* optimal d for the optimal split */
  int      best_j;		/* optimal j for the optimal split */
  int      tv;			/* remember the position of a bifurc in the trace. */
  int      b1,b2;		/* argmax_v for 0->v local begin transitions */
  float    b1_sc, b2_sc;	/* max_v scores for 0->v local begin transitions */

  /* 1. If the generic problem is small enough, solve it with insideT,
   *    and append the trace to tr.
   */
  if (insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT) {
    ESL_DPRINTF2(("#DEBUG: Solving a generic w/ insideT - G%d[%s]..%d[%s], %d..%d\n",
		  r, UniqueStatetype(cm->stid[r]),
		  z, UniqueStatetype(cm->stid[z]),
		  i0, j0));
    sc = insideT(cm, dsq, L, tr, r, z, i0, j0, (r==0), dmin, dmax);
    return sc;
  }

  /* 2. Traverse down from r, find first bifurc.
   *    The lowest a bifurc could be: B-S-E/S-IL-E = vend-5
   *                                   
   */
  for (v = r; v <= z-5; v++)
    if (cm->sttype[v] == B_st) break; /* found the first bifurcation, now v */

  /* 3. If there was no bifurcation, this is a wedge problem; solve it
   *    with wedge_splitter. 
   */
  if (v > z-5) {		/* no bifurc? it's a wedge problem  */
    if (cm->sttype[z] != E_st) cm_Fail("inconceivable.");
    sc = wedge_splitter_qdb(cm, dsq, L, tr, r, z, i0, j0, dmin, dmax);
    return sc;
  }

  /* Set up the state quartet r,v,w,y for a divide and conquer
   * solution of the generic problem.
   */
  w = cm->cfirst[v];		/* index of left S  */
  y = cm->cnum[v];		/* index right S    */
  if (w < y) { wend = y-1; yend = z; }
  else       { yend = w-1; wend = z; }

  /* Calculate alpha[w] deck and alpha[y] deck.
   * We also get b1: best choice for 0->b local begin. b1_sc is the score if we do this.
   * Analogous for b2, b2_sc on the other side.
   */
  inside_qdb(cm, dsq, L, w, wend, i0, j0, BE_EFFICIENT, NULL,  &alpha, NULL, &pool, NULL, 
	   (r==0), &b1, &b1_sc, dmin, dmax);
  inside_qdb(cm, dsq, L, y, yend, i0, j0, BE_EFFICIENT, alpha, &alpha, pool, &pool, NULL,
	   (r==0), &b2, &b2_sc, dmin, dmax);

  /* Calculate beta[v] deck (stick it in alpha). Let the pool get free'd.
   * (If we're doing local alignment, deck M is the beta[EL] deck.)
   */
  outside_qdb(cm, dsq, L, r, v, i0, j0, BE_EFFICIENT, alpha, &beta, pool, NULL, dmin, dmax);

  /* Find the optimal split at the B.
   */
  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  for (jp = 0; jp <= W; jp++) 
    {
      j = i0-1+jp;
      /* Bands used */
      /* old line : for (d = 0; d <= jp; d++) */
      for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	for (k = 0; k <= d; k++)
	  if ((sc = alpha[w][j-k][d-k] + alpha[y][j][k] + beta[v][j][d]) > best_sc) 
	    {
	      best_sc = sc;
	      best_k  = k;
	      best_j  = j;
	      best_d  = d;
	    }
    }

  /* Local alignment only: maybe we're better off in EL?
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = 0; jp <= W; jp++) 
      {
	j = i0-1+jp;
	/* There is no band on the EL state */
	for (d = 0; d <= jp; d++) 
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc;
	    best_k  = -1;	/* special flag for local end, EL. */
	    best_j  = j;
	    best_d  = d;
	  }
      }
  }
  
  /* Local alignment only: maybe we're better off in ROOT?
   */
  if (r == 0 && cm->flags & CMH_LOCAL_BEGIN) {
    if (b1_sc > best_sc) {
      best_sc = b1_sc;
      best_k  = -2;		/* flag for using local begin into left wedge w..wend */
      best_j  = j0;		
      best_d  = W;
    }
    if (b2_sc > best_sc) {
      best_sc = b2_sc;
      best_k  = -3;		/* flag for using local begin into right wedge y..yend */
      best_j  = j0;		
      best_d  = W;
    }
  }

  /* Free now, before recursing.
   * The two alpha matrices and the beta matrix
   * actually all point to the same memory, since no
   * decks in Inside and Outside needed to overlap. 
   * Free 'em all in one call.
   */
  free_vjd_matrix(alpha, cm->M, i0, j0);

  /* If we're in EL, instead of B, the optimal alignment is entirely
   * in a V problem that's still above us. The TRUE flag sets useEL.
   */
  if (best_k == -1) {	
    v_splitter_qdb(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, TRUE, dmin, dmax);    
    return best_sc;
  } 

  /* Else: if we're in the root 0, we know which r we did our local begin into.
   * We have a generic problem rooted there. The FALSE flag disallows
   * any further local begins.
   */
  if (best_k == -2) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b1);
    z = CMSubtreeFindEnd(cm, b1);
    generic_splitter_qdb(cm, dsq, L, tr, b1, z, i0, j0, dmin, dmax);
    return best_sc;
  }
  if (best_k == -3) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b2);
    z = CMSubtreeFindEnd(cm, b2);
    generic_splitter_qdb(cm, dsq, L, tr, b2, z, i0, j0, dmin, dmax);
    return best_sc;
  }

  /* Else (the usual case), ok, we did use B in the optimal split.
   * Split now into a V problem and two generic problems, and recurse
   * left fragment: i1 = j-d+1, j1 = j-k, vroot = w, vend = wend
   * right frag:    i2 = j-k+1, j2 = j,   vroot = y, vend = yend
   * 
   * The problems must be solved in a particular order, since we're
   * constructing the trace in a postorder traversal.
   */
  ESL_DPRINTF2(("#DEBUG: Generic splitter:\n"));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		v, UniqueStatetype(cm->stid[v]),
		i0, best_j-best_d+1, best_j, j0));
  ESL_DPRINTF2(("#DEBUG:    generic: G%d[%s]..%d[%s], %d..%d\n", 
		w,    UniqueStatetype(cm->stid[w]),
		wend, UniqueStatetype(cm->stid[wend]),
		best_j-best_d+1, best_j-best_k));
  ESL_DPRINTF2(("#DEBUG:    generic: G%d[%s]..%d[%s], %d..%d\n", 
		y,    UniqueStatetype(cm->stid[y]),
		yend, UniqueStatetype(cm->stid[yend]),
		best_j-best_k+1, best_j));

  v_splitter_qdb(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, FALSE, dmin, dmax);
  tv = tr->n-1;

  InsertTraceNode(tr, tv, TRACE_LEFT_CHILD, best_j-best_d+1, best_j-best_k, w);
  generic_splitter_qdb(cm, dsq, L, tr, w, wend, best_j-best_d+1, best_j-best_k, dmin, dmax);
  InsertTraceNode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j, y);
  generic_splitter_qdb(cm, dsq, L, tr, y, yend, best_j-best_k+1, best_j, dmin, dmax);

  return best_sc;
}

/* Function: wedge_splitter_qdb()
 *           EPN 05.19.05
 * *based on wedge_splitter(), only difference is bands are used : 
 * Date:     SRE, Sun May 13 08:44:15 2001 [CSHL genome mtg]
 *
 * Purpose:  Solve a "wedge problem": best parse of an 
 *           unbifurcated subgraph cm^r..z to a substring
 *           dsq[i0..j0]. r may be a start state (when
 *           the wedge problem comes from being a special case
 *           of a generic problem) or a non-insert state
 *           (D, MP, ML, MR) (when the wedge comes from a
 *           previous wedge_splitter), or indeed, any non-end
 *           state (when wedge comes from a local begin).
 *           z, however, is always an end state.
 *           
 *           Attaches the optimal trace T(r..z), exclusive
 *           of r and inclusive of z, to the growing trace tr.
 *           
 *           Deal with a divide and conquer boundary condition:
 *           the next non-insert state after r is the end state z.
 *           All remaining sequence of i0..j0 that r doesn't emit
 *           must be dealt with by insert states.
 *
 * Args:     cm          - model
 *           dsq         - digitized sequence 1..L
 *           L           - length of dsq
 *           tr          - the traceback we're adding on to.
 *           r           - index of the first state in the subgraph
 *           z           - index of an end state (E_st) in the model
 *           i0          - start in the sequence (1..L)
 *           j0          - end in the sequence (1..L)
 *           dmin   - minimum d bound for each state v; [0..v..M-1]
 *           dmax   - maximum d bound for each state v; [0..v..M-1]
 *
 * Returns:  The score of the best parse in bits.
 */
static float 
wedge_splitter_qdb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int j0,
		 int *dmin, int *dmax)
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
  int   b;	/* optimal local begin: b = argmax_v alpha_v(i0,j0) + t_0(v) */
  float bsc;	/* score for optimal local begin      */
  
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
      ESL_DPRINTF2(("#DEBUG: Solving a wedge:   G%d[%s]..%d[%s], %d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		z, UniqueStatetype(cm->stid[z]),
		i0,j0));
      sc = insideT(cm, dsq, L, tr, r, z, i0, j0, (r==0), dmin, dmax);
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
   *    b will contain the optimal 0->v state for a local begin, and bsc
   *    is the score for using it.
   *    beta[cm->M] will contain the EL deck, if needed for local ends.
   */
  inside_qdb(cm, dsq, L, w, z, i0, j0, BE_EFFICIENT, 
	   NULL, &alpha, NULL, &pool, NULL, 
	   (r==0), &b, &bsc, dmin, dmax);
  outside_qdb(cm, dsq, L, r, y, i0, j0, BE_EFFICIENT, NULL, &beta, pool, NULL,
  dmin, dmax);

  /* 4. Find the optimal split at the split set: best_v, best_d, best_j
   */
  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  for (v = w; v <= y; v++)
    for (jp = 0; jp <= W; jp++) 
      {
	j = i0-1+jp;
	for (d = dmin[v]; d <= dmax[v] && d <= jp; d++) 
	  if ((sc = alpha[v][j][d] + beta[v][j][d]) > best_sc)
	    {
	      best_sc = sc;
	      best_v  = v;
	      best_d  = d;
	      best_j  = j;
	    }
      }

  /* Local alignment ends only: maybe we're better off in EL, 
   * not in the split set?
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = 0; jp <= W; jp++) 
      {
	j = i0-1+jp;
	/* There is no band on the EL state */
	for (d = 0; d <= jp; d++) 
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc;
	    best_v  = -1;	/* flag for local alignment. */
	    best_j  = j;
	    best_d  = d;
	  }
      }
  }

  /* Local alignment begins only: maybe we're better off in the root.
   */
  if (r==0 && (cm->flags & CMH_LOCAL_BEGIN)) {
    if (bsc > best_sc) {
      best_sc = bsc;
      best_v  = -2;		/* flag for local alignment */
      best_j  = j0;
      best_d  = W;
    }
  }

  /* free now, before recursing!
   */
  free_vjd_matrix(alpha, cm->M, i0, j0);
  free_vjd_matrix(beta,  cm->M, i0, j0);

  /* If we're in EL, instead of the split set, the optimal alignment
   * is entirely in a V problem that's still above us. The TRUE
   * flag sets useEL. It doesn't matter which state in the split
   * set w..y we use as the end of the graph; vinside() will have to
   * initialize the whole thing to IMPOSSIBLE anyway.
   */  
  if (best_v == -1) {
    v_splitter_qdb(cm, dsq, L, tr, r, w, i0, best_j-best_d+1, best_j, j0, TRUE, dmin, dmax);    
    return best_sc;
  }

  /* If we're in the root because of a local begin, the local alignment
   * is entirely in a wedge problem that's still below us, rooted at b.
   * The FALSE flag prohibits any more local begins in this and subsequent
   * problems. 
   */
  if (best_v == -2) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b);
    wedge_splitter_qdb(cm, dsq, L, tr, b, z, i0, j0, dmin, dmax);
    return best_sc; 
  }

  /* Else (usual case): the optimal split into a V problem and a wedge problem:
   *    i1 = best_j-best_d+1, j1 = best_j
   *    the V problem:     r..v, i0..i1, j1..j0
   *    the wedge problem: v..z, i1..j1
   *    
   *    These have to solved in the order given because we're
   *    constructing the trace in postorder traversal.
   */
  ESL_DPRINTF2(("#DEBUG: Wedge splitter:\n"));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		best_v, UniqueStatetype(cm->stid[best_v]),
		i0, best_j-best_d+1, best_j, j0));
  ESL_DPRINTF2(("#DEBUG:    wedge:   G%d[%s]..%d[%s], %d..%d\n", 
		best_v, UniqueStatetype(cm->stid[best_v]),
		z, UniqueStatetype(cm->stid[z]),
		best_j-best_d+1, best_j));

  v_splitter_qdb(cm, dsq, L, tr, r, best_v, i0, best_j-best_d+1, best_j, j0, FALSE,
	       dmin, dmax);
  wedge_splitter_qdb(cm, dsq, L, tr, best_v, z, best_j-best_d+1, best_j, dmin, dmax);
  return best_sc;
}



/*################################################################
 * The HMM-banded (CP9) dividers, conquerors, and engines (*_hb).
 * See the block comment near the top of this file for the design.
 *################################################################*/

/*****************************************************************
 * Banded-cell vjd deck allocation for the HMM-banded D&C (brief 26_0610-008, Stage 1b).
 *
 * Unlike Stage 1a.1 (full O(W^2)-triangle decks via alloc_vjd_deck()), a Stage
 * 1b deck for state v stores ONLY the in-band cells -- exactly the cell set
 * CM_HB_MX holds. Row indexing stays PLATONIC (a[j], j in subproblem [i0-1..j0]
 * intersected with v's j-band [jmin[v]..jmax[v]]); only the d dimension is
 * offset-indexed: a[j][d - hdmin[v][j-jmin[v]]], per-row width
 * hdmax[v][jp_v]-hdmin[v][jp_v]+1. Out-of-band rows are NULL; out-of-band cells
 * simply don't exist. Cross-state reads must therefore be guarded + offset by
 * hb_inband() (the analogue of the full-deck IMPOSSIBLE init 1a.1 relied on).
 *
 * Banded decks have per-state shapes, so they are NOT pooled (a popped deck of
 * the wrong shape would corrupt memory); the *_hb engines allocate fresh and
 * free immediately. The CYKDeckTrack high-water counts the BANDED bytes.
 *****************************************************************/

/* Storage offset of cell (v,j,d) in a banded-HB deck spanning subproblem
 * [i0-1..j0], or "not present" (return 0). Platonic j; d offset by hdmin. */
static int
hb_inband(CP9Bands_t *cp9b, int v, int j, int d, int i0, int j0, int *ret_dp)
{
  int jp_v;
  if (j < i0-1 || j > j0)                       return 0; /* outside this subproblem deck */
  if (j < cp9b->jmin[v] || j > cp9b->jmax[v])   return 0; /* outside v's j-band            */
  jp_v = j - cp9b->jmin[v];
  if (d < hd_min(cp9b, v, jp_v) || d > hd_max(cp9b, v, jp_v)) return 0; /* outside d-band     */
  *ret_dp = d - hd_min(cp9b, v, jp_v);
  return 1;
}

/* data bytes (per-row float arrays) of a banded-HB deck for state v over
 * subproblem [i0-1..j0]; matches what alloc/free add/subtract from the counter. */
static double
banded_hb_vjd_deck_bytes(int i0, int j0, int v, CP9Bands_t *cp9b)
{
  int    j, jlo, jhi, jp_v, w;
  double tot = 0.;
  jlo = ESL_MAX(i0-1, cp9b->jmin[v]);
  jhi = ESL_MIN(j0,   cp9b->jmax[v]);
  for (j = jlo; j <= jhi; j++) {
    jp_v = j - cp9b->jmin[v];
    w = hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v) + 1;
    if (w > 0) tot += (double) sizeof(float) * (double) w;
  }
  return tot;
}

static float **
alloc_banded_hb_vjd_deck(int L, int i0, int j0, int v, CP9Bands_t *cp9b)
{
  int     status;
  float **a;
  int     j, jlo, jhi, jp_v, w;
  ESL_ALLOC(a, sizeof(float *) * (L+1));
  for (j = 0; j <= L; j++) a[j] = NULL;
  jlo = ESL_MAX(i0-1, cp9b->jmin[v]);
  jhi = ESL_MIN(j0,   cp9b->jmax[v]);
  for (j = jlo; j <= jhi; j++) {
    jp_v = j - cp9b->jmin[v];
    w = hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v) + 1;
    if (w > 0) ESL_ALLOC(a[j], sizeof(float) * w);
    else       a[j] = NULL;
  }
  if (cyk_dnc_track) {
    cyk_dnc_cur_bytes += banded_hb_vjd_deck_bytes(i0, j0, v, cp9b);
    cyk_dnc_note();
  }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

static void
free_banded_hb_vjd_deck(float **a, int i0, int j0, int v, CP9Bands_t *cp9b)
{
  int j, jlo, jhi;
  if (a == NULL) return;
  if (cyk_dnc_track) cyk_dnc_cur_bytes -= banded_hb_vjd_deck_bytes(i0, j0, v, cp9b);
  jlo = ESL_MAX(i0-1, cp9b->jmin[v]);
  jhi = ESL_MIN(j0,   cp9b->jmax[v]);
  for (j = jlo; j <= jhi; j++) if (a[j] != NULL) free(a[j]);
  free(a);
}

static char **
alloc_banded_hb_vjd_yshadow_deck(int L, int i0, int j0, int v, CP9Bands_t *cp9b)
{
  int    status;
  char **a;
  int    j, jlo, jhi, jp_v, w;
  ESL_ALLOC(a, sizeof(char *) * (L+1));
  for (j = 0; j <= L; j++) a[j] = NULL;
  jlo = ESL_MAX(i0-1, cp9b->jmin[v]);
  jhi = ESL_MIN(j0,   cp9b->jmax[v]);
  { double nb = 0.0;
  for (j = jlo; j <= jhi; j++) {
    jp_v = j - cp9b->jmin[v];
    w = hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v) + 1;
    if (w > 0) { ESL_ALLOC(a[j], sizeof(char) * w); nb += (double) w * sizeof(char); }
    else       a[j] = NULL;
  }
  cyk_dnc_shad_add(nb); }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

static int **
alloc_banded_hb_vjd_kshadow_deck(int L, int i0, int j0, int v, CP9Bands_t *cp9b)
{
  int    status;
  int  **a;
  int    j, jlo, jhi, jp_v, w;
  ESL_ALLOC(a, sizeof(int *) * (L+1));
  for (j = 0; j <= L; j++) a[j] = NULL;
  jlo = ESL_MAX(i0-1, cp9b->jmin[v]);
  jhi = ESL_MIN(j0,   cp9b->jmax[v]);
  { double nb = 0.0;
  for (j = jlo; j <= jhi; j++) {
    jp_v = j - cp9b->jmin[v];
    w = hd_max(cp9b, v, jp_v) - hd_min(cp9b, v, jp_v) + 1;
    if (w > 0) { ESL_ALLOC(a[j], sizeof(int) * w); nb += (double) w * sizeof(int); }
    else       a[j] = NULL;
  }
  cyk_dnc_shad_add(nb); }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

/* free a banded-HB score matrix (alpha/beta): each deck v is banded, EXCEPT the
 * EL deck M (full, unbanded -- "no band on EL"). NULL decks skipped. The banded
 * shadow matrix is freed by the existing free_vjd_shadow_matrix() (NULL-safe). */
static void
free_banded_hb_vjd_matrix(float ***a, CM_t *cm, int i0, int j0, CP9Bands_t *cp9b)
{
  int v;
  for (v = 0; v <= cm->M; v++) {
    if (a[v] == NULL) continue;
    if (v == cm->M) free_vjd_deck(a[v], i0, j0);                 /* EL deck: full */
    else            free_banded_hb_vjd_deck(a[v], i0, j0, v, cp9b);
    a[v] = NULL;
  }
  free(a);
}

/*****************************************************************
 * Banded-cell vji deck allocation for the HMM-banded D&C (brief 26_0610-009, Stage 1a.2).
 *
 * V problems are solved in the vji coordinate system: the Platonic deck for state
 * v over a one-hole subseq [i0..i1]//[j1..j0] is a[jp][ip], jp=j-j1 in [0..j0-j1],
 * ip=i-i0 in [0..i1-i0], with d = j-i+1 (the full subsequence length subtended by
 * v, INCLUDING the hole -- exactly what hdmin/hdmax band). A Stage 1a.2 deck
 * stores ONLY the in-band cells: row jp (j=jp+j1) is present iff j is in v's
 * j-band [jmin[v]..jmax[v]], and within that row only i in [ilo..ihi], where the
 * d-band [hdmin..hdmax] at jpb=j-jmin[v] maps to i = j-d+1, i.e.
 *   ilo = max(i0, j - hdmax[v][jpb] + 1),  ihi = min(i1, j - hdmin[v][jpb] + 1).
 * Cell access is a[jp][i - ilo(jp)] via vji_inband(). Row indexing stays platonic
 * (outer dim j0-j1+1, out-of-band rows NULL) so the existing NULL-safe vji shadow
 * free routines apply; only the i dimension is offset, mirroring the vjd Stage 1b.
 *
 * THE HAZARD (brief 26_0610-003): the QDB vji path crashed with a negative malloc when a
 * collapsed band left v_splitter_qdb with degenerate corners. Here, with tight
 * per-(v,j) CP9 bands, empty rows are NORMAL, so the allocator NEVER allocs a
 * non-positive row (w<=0 -> NULL row), and the outer dim is clamped >= 1. The
 * degenerate-recursion hazard itself is designed out in v_splitter_hb (the
 * explicit empty-band path), not here. Banded decks are NOT pooled (per-(v,j)
 * shapes); engines alloc fresh + free immediately. CYKDeckTrack counts the banded
 * bytes via cyk_dnc_vji_bytes (the same counter the full vji decks used).
 *****************************************************************/

/* Storage offset of cell (v, platonic j, platonic i) in a banded-HB vji deck for
 * a V-problem [i0..i1]//[j1..j0], or "not present" (return 0). vji analogue of
 * hb_inband(): row jp=j-j1 (caller indexes a[jp]); the i dimension is offset by
 * the per-row in-band minimum ilo. */
static int
vji_inband(CP9Bands_t *cp9b, int v, int j, int i, int i0, int i1, int j1, int j0, int *ret_ip)
{
  int jpb, d, ilo;
  if (j < j1 || j > j0)                         return 0; /* outside V-problem j-range */
  if (j < cp9b->jmin[v] || j > cp9b->jmax[v])   return 0; /* outside v's j-band        */
  if (i < i0 || i > i1)                         return 0; /* outside V-problem i-range */
  jpb = j - cp9b->jmin[v];
  d   = j - i + 1;
  if (d < hd_min(cp9b, v, jpb) || d > hd_max(cp9b, v, jpb)) return 0; /* outside d-band  */
  ilo = j - hd_max(cp9b, v, jpb) + 1;            /* smallest in-band i (before clamp)   */
  if (ilo < i0) ilo = i0;
  *ret_ip = i - ilo;
  return 1;
}

/* data bytes (per-row float arrays) of a banded-HB vji deck for state v; matches
 * what alloc/free add/subtract from cyk_dnc_vji_bytes. */
static double
banded_hb_vji_deck_bytes(int i0, int i1, int j1, int j0, int v, CP9Bands_t *cp9b)
{
  int    j, jpb, ilo, ihi, w;
  double tot = 0.;
  for (j = ESL_MAX(j1, cp9b->jmin[v]); j <= ESL_MIN(j0, cp9b->jmax[v]); j++) {
    jpb = j - cp9b->jmin[v];
    ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
    ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
    w   = ihi - ilo + 1;
    if (w > 0) tot += (double) sizeof(float) * (double) w;
  }
  return tot;
}

static float **
alloc_banded_hb_vji_deck(int i0, int i1, int j1, int j0, int v, CP9Bands_t *cp9b)
{
  int     status;
  float **a;
  int     jp, j, jpb, ilo, ihi, w;
  int     njp = j0 - j1 + 1;
  if (njp < 1) njp = 1;                         /* defensive: never alloc < 1 row ptr  */
  ESL_ALLOC(a, sizeof(float *) * njp);
  for (jp = 0; jp < njp; jp++) a[jp] = NULL;
  for (j = ESL_MAX(j1, cp9b->jmin[v]); j <= ESL_MIN(j0, cp9b->jmax[v]); j++) {
    jp  = j - j1;
    jpb = j - cp9b->jmin[v];
    ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
    ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
    w   = ihi - ilo + 1;
    if (w > 0) ESL_ALLOC(a[jp], sizeof(float) * w);  /* w>0 guaranteed by guard        */
    else       a[jp] = NULL;
  }
  if (cyk_dnc_track) {
    cyk_dnc_vji_bytes += banded_hb_vji_deck_bytes(i0, i1, j1, j0, v, cp9b);
    cyk_dnc_note();
  }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

static void
free_banded_hb_vji_deck(float **a, int i0, int i1, int j1, int j0, int v, CP9Bands_t *cp9b)
{
  int jp;
  if (a == NULL) return;
  if (cyk_dnc_track) cyk_dnc_vji_bytes -= banded_hb_vji_deck_bytes(i0, i1, j1, j0, v, cp9b);
  for (jp = 0; jp <= j0-j1; jp++) if (a[jp] != NULL) free(a[jp]);
  free(a);
}

static char **
alloc_banded_hb_vji_shadow_deck(int i0, int i1, int j1, int j0, int v, CP9Bands_t *cp9b)
{
  int    status;
  char **a;
  int    jp, j, jpb, ilo, ihi, w;
  int    njp = j0 - j1 + 1;
  if (njp < 1) njp = 1;
  ESL_ALLOC(a, sizeof(char *) * njp);
  for (jp = 0; jp < njp; jp++) a[jp] = NULL;
  { double nb = 0.0;
  for (j = ESL_MAX(j1, cp9b->jmin[v]); j <= ESL_MIN(j0, cp9b->jmax[v]); j++) {
    jp  = j - j1;
    jpb = j - cp9b->jmin[v];
    ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
    ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
    w   = ihi - ilo + 1;
    if (w > 0) { ESL_ALLOC(a[jp], sizeof(char) * w); nb += (double) w * sizeof(char); }
    else       a[jp] = NULL;
  }
  cyk_dnc_shad_add(nb); }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

/* Set every in-band cell of a banded-HB vji deck for state v to IMPOSSIBLE. */
static void
banded_hb_vji_init_impossible(float **a, int i0, int i1, int j1, int j0, int v, CP9Bands_t *cp9b)
{
  int j, jp, jpb, ilo, ihi, i;
  for (j = ESL_MAX(j1, cp9b->jmin[v]); j <= ESL_MIN(j0, cp9b->jmax[v]); j++) {
    jp  = j - j1;
    jpb = j - cp9b->jmin[v];
    ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
    ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
    for (i = ilo; i <= ihi; i++) a[jp][i - ilo] = IMPOSSIBLE;
  }
}

/* free a banded-HB vji score matrix (alpha/beta): each deck v is banded, EXCEPT
 * the EL deck M (full vji, unbanded -- "no band on EL"). NULL decks skipped. The
 * banded vji shadow matrix is freed by the existing free_vji_shadow_matrix()
 * (NULL-safe; out-of-band rows are NULL). */
static void
free_banded_hb_vji_matrix(float ***a, CM_t *cm, int i0, int i1, int j1, int j0, CP9Bands_t *cp9b)
{
  int v;
  for (v = 0; v <= cm->M; v++) {
    if (a[v] == NULL) continue;
    if (v == cm->M) free_vji_deck(a[v], j1, j0);              /* EL deck: full vji */
    else            free_banded_hb_vji_deck(a[v], i0, i1, j1, j0, v, cp9b);
    a[v] = NULL;
  }
  free(a);
}

/* Function: generic_splitter_hb()
 *           EPN 2026 [brief 26_0610-007]
 *
 * Purpose:  HMM-banded analogue of generic_splitter_qdb(). Identical control
 *           flow, but uses the *_hb engines (inside_hb/outside_hb/insideT_hb)
 *           and the CP9 bands in <cp9b> (per-v j-band + per-(v,j) d-band) when
 *           searching for the optimal bifurcation split. V problems are handed
 *           to the EXACT (unbanded) v_splitter() (Stage 1a.1 class-2 scope).
 */
static float
generic_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
		    int r, int z, int i0, int j0, CP9Bands_t *cp9b)
{
  float ***alpha;
  float ***beta;
  int     *eldmax = NULL;	/* brief 26_0430-227: banded EL deck per-row d-band edges */
  struct deckpool_s *pool;
  int      v,w,y;		/* state indices */
  int      wend, yend;		/* indices for end of subgraphs rooted at w,y */
  int      jp;			/* j': relative position in subseq, 0..W */
  int      jp_v;                /* j index in v's band: j - jmin[v] */
  int      W;			/* length of subseq i0..j0 */
  float    sc;			/* tmp variable for a score */
  int      j,d,k;		/* sequence indices */
  float    best_sc;		/* optimal score at the optimal split point */
  int      best_k;		/* optimal k for the optimal split */
  int      best_d;		/* optimal d for the optimal split */
  int      best_j;		/* optimal j for the optimal split */
  int      tv;			/* remember the position of a bifurc in the trace. */
  int      b1,b2;		/* argmax_v for 0->v local begin transitions */
  float    b1_sc, b2_sc;	/* max_v scores for 0->v local begin transitions */
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  /* 1. If the generic problem is small enough, solve it with insideT_hb. */
  if (insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT) {
    sc = insideT_hb(cm, dsq, L, tr, r, z, i0, j0, (r==0), cp9b);
    return sc;
  }

  /* 2. Traverse down from r, find first bifurc. */
  for (v = r; v <= z-5; v++)
    if (cm->sttype[v] == B_st) break;

  /* 3. No bifurcation -> wedge problem. */
  if (v > z-5) {
    if (cm->sttype[z] != E_st) cm_Fail("inconceivable.");
    sc = wedge_splitter_hb(cm, dsq, L, tr, r, z, i0, j0, cp9b);
    return sc;
  }

  /* Set up the state quartet r,v,w,y. */
  w = cm->cfirst[v];
  y = cm->cnum[v];
  if (w < y) { wend = y-1; yend = z; }
  else       { yend = w-1; wend = z; }

  /* Calculate alpha[w] and alpha[y] decks, and beta[v]. */
  inside_hb(cm, dsq, L, w, wend, i0, j0, BE_EFFICIENT, NULL,  &alpha, NULL, &pool, NULL,
	    (r==0), &b1, &b1_sc, cp9b);
  inside_hb(cm, dsq, L, y, yend, i0, j0, BE_EFFICIENT, alpha, &alpha, pool, &pool, NULL,
	    (r==0), &b2, &b2_sc, cp9b);
  outside_hb(cm, dsq, L, r, v, i0, j0, BE_EFFICIENT, alpha, &beta, pool, NULL, cp9b, &eldmax);

  /* Find the optimal split at the B, within v's bands. */
  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++)
    {
      jp   = j - (i0-1);
      jp_v = j - jmin[v];
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	{
	  int dp_v = d - hd_min(cp9b, v, jp_v);           /* beta[v] offset (v,j,d in band) */
	  for (k = 0; k <= d; k++)
	    {
	      int dp_w, dp_y; /* the w (at j-k,d-k) and y (at j,k) child cells must
			       * both be in-band; else this split is IMPOSSIBLE */
	      if (hb_inband(cp9b, w, j-k, d-k, i0, j0, &dp_w) &&
		  hb_inband(cp9b, y, j,   k,   i0, j0, &dp_y) &&
		  (sc = alpha[w][j-k][dp_w] + alpha[y][j][dp_y] + beta[v][j][dp_v]) > best_sc)
		{
		  best_sc = sc;
		  best_k  = k;
		  best_j  = j;
		  best_d  = d;
		}
	    }
	}
    }

  /* Local alignment only: maybe we're better off in EL?  brief 26_0430-227: the
   * EL deck beta[cm->M] is now banded (row j holds d=0..eldmax[j]); cells above
   * eldmax[j] were IMPOSSIBLE in the old full deck, so clamping the read to
   * eldmax[j] leaves best_sc byte-identical.  eldmax[j]<0 => empty (NULL) row. */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = 0; jp <= W; jp++)
      {
	int dhi;
	j = i0-1+jp;
	if (eldmax[j] < 0) continue;
	dhi = ESL_MIN(jp, eldmax[j]);
	for (d = 0; d <= dhi; d++)
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc;
	    best_k  = -1;
	    best_j  = j;
	    best_d  = d;
	  }
      }
  }

  /* Local alignment only: maybe we're better off in ROOT? */
  if (r == 0 && cm->flags & CMH_LOCAL_BEGIN) {
    if (b1_sc > best_sc) { best_sc = b1_sc; best_k = -2; best_j = j0; best_d = W; }
    if (b2_sc > best_sc) { best_sc = b2_sc; best_k = -3; best_j = j0; best_d = W; }
  }

  /* brief 26_0430-227: here outside_hb() was handed the inside `alpha` deck-
   * array as its working beta array (the call above passes `alpha` as the beta
   * INPUT arg), so the returned `beta` IS `alpha` -- the SAME float*** pointer,
   * with the EL deck stored in the shared alpha[cm->M]==beta[cm->M] slot.  So
   * `beta` must NOT be freed separately (that double-frees; this is exactly what
   * made brief 26_0430-226's attempted "beta leak" fix crash -- there was never
   * a separate beta array to leak, only the EL deck, which the alpha free below
   * already reclaimed).  Free the banded EL deck explicitly first (eldmax-driven
   * width + correct cyk_dnc_track accounting, vs the full-triangle free_vjd_deck
   * that free_banded_hb_vjd_matrix's v==cm->M branch would use), NULL the shared
   * slot, then free the shared array once. */
  if ((cm->flags & CMH_LOCAL_END) && beta[cm->M] != NULL) {
    free_el_banded_vjd_deck(beta[cm->M], i0, j0, eldmax);
    beta[cm->M] = NULL;      /* == alpha[cm->M]: alpha's matrix-free then skips it */
  }
  free_banded_hb_vjd_matrix(alpha, cm, i0, j0, cp9b);
  free(eldmax); eldmax = NULL;

  /* EL case: V problem above us; solve with banded v_splitter_hb (Stage 1a.2). */
  if (best_k == -1) {
    v_splitter_hb(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, TRUE, cp9b);
    return best_sc;
  }
  if (best_k == -2) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b1);
    z = CMSubtreeFindEnd(cm, b1);
    generic_splitter_hb(cm, dsq, L, tr, b1, z, i0, j0, cp9b);
    return best_sc;
  }
  if (best_k == -3) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b2);
    z = CMSubtreeFindEnd(cm, b2);
    generic_splitter_hb(cm, dsq, L, tr, b2, z, i0, j0, cp9b);
    return best_sc;
  }

  /* Usual case: banded V problem (Stage 1a.2) + two banded generic problems. */
  v_splitter_hb(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, FALSE, cp9b);
  tv = tr->n-1;

  InsertTraceNode(tr, tv, TRACE_LEFT_CHILD, best_j-best_d+1, best_j-best_k, w);
  generic_splitter_hb(cm, dsq, L, tr, w, wend, best_j-best_d+1, best_j-best_k, cp9b);
  InsertTraceNode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j, y);
  generic_splitter_hb(cm, dsq, L, tr, y, yend, best_j-best_k+1, best_j, cp9b);

  return best_sc;
}

/* Function: wedge_splitter_hb()
 *           EPN 2026 [brief 26_0610-007]
 *
 * Purpose:  HMM-banded analogue of wedge_splitter_qdb(). V problems solved with
 *           the EXACT (unbanded) v_splitter().
 */
static float
wedge_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int j0,
		  CP9Bands_t *cp9b)
{
  float ***alpha;
  float ***beta;
  int     *eldmax = NULL;	/* brief 26_0430-227: banded EL deck per-row d-band edges */
  struct deckpool_s *pool;
  float sc;
  float best_sc;
  int   v,w,y;
  int   W;
  int   d, jp, j, jp_v;
  int   best_v, best_d, best_j;
  int   midnode;
  int   b;
  float bsc;
  int  *jmin  = cp9b->jmin;
  int  *jmax  = cp9b->jmax;
  int **hdmin = cp9b->hdmin;
  int **hdmax = cp9b->hdmax;

  /* 1. Boundary condition or small enough -> insideT_hb. */
  if (cm->ndidx[z] == cm->ndidx[r] + 1 ||
      insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT)
    {
      sc = insideT_hb(cm, dsq, L, tr, r, z, i0, j0, (r==0), cp9b);
      return sc;
    }

  /* 2. Find our split set, w..y. */
  midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r]) / 2);
  w = cm->nodemap[midnode];
  y = cm->cfirst[w]-1;

  /* 3. Banded inside up to w, banded outside down to y. */
  inside_hb(cm, dsq, L, w, z, i0, j0, BE_EFFICIENT,
	    NULL, &alpha, NULL, &pool, NULL,
	    (r==0), &b, &bsc, cp9b);
  outside_hb(cm, dsq, L, r, y, i0, j0, BE_EFFICIENT, NULL, &beta, pool, NULL, cp9b, &eldmax);

  /* 4. Find the optimal split at the split set, within each v's bands. */
  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  for (v = w; v <= y; v++)
    for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++)
      {
	jp   = j - (i0-1);
	jp_v = j - jmin[v];
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	  {
	    int dp_v = d - hd_min(cp9b, v, jp_v);          /* v,j,d in band by loop bounds */
	    if ((sc = alpha[v][j][dp_v] + beta[v][j][dp_v]) > best_sc)
	      {
		best_sc = sc;
		best_v  = v;
		best_d  = d;
		best_j  = j;
	      }
	  }
      }

  /* Local ends: maybe better in EL?  brief 26_0430-227: EL deck now banded on
   * the upper d-edge (eldmax[j]); cells above were IMPOSSIBLE, so clamping the
   * read to eldmax[j] is byte-identical.  eldmax[j]<0 => empty (NULL) row. */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = 0; jp <= W; jp++)
      {
	int dhi;
	j = i0-1+jp;
	if (eldmax[j] < 0) continue;
	dhi = ESL_MIN(jp, eldmax[j]);
	for (d = 0; d <= dhi; d++)
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc;
	    best_v  = -1;
	    best_j  = j;
	    best_d  = d;
	  }
      }
  }

  /* Local begins: maybe better in root? */
  if (r==0 && (cm->flags & CMH_LOCAL_BEGIN)) {
    if (bsc > best_sc) { best_sc = bsc; best_v = -2; best_j = j0; best_d = W; }
  }

#ifdef HBDNC_DEBUG
  if (! NOT_IMPOSSIBLE(best_sc)) {
    fprintf(stderr, "## WEDGE_HB no split: r=%d z=%d i0=%d j0=%d W=%d  splitset w=%d y=%d\n",
	    r, z, i0, j0, W, w, y);
    { int vv; for (vv=w; vv<=y; vv++) {
        int na=0, nb=0, jj, dd; int ajlo=99999,ajhi=-1,bjlo=99999,bjhi=-1; int both=0;
        for (jj=ESL_MAX(i0-1,jmin[vv]); jj<=ESL_MIN(j0,jmax[vv]); jj++) {
          int jpv=jj-jmin[vv], jpp=jj-(i0-1);
          for (dd=hd_min(cp9b, vv, jpv); dd<=hd_max(cp9b, vv, jpv) && dd<=jpp; dd++) {
            int ddp=dd-hd_min(cp9b, vv, jpv);
            int va=NOT_IMPOSSIBLE(alpha[vv][jj][ddp]), vb=NOT_IMPOSSIBLE(beta[vv][jj][ddp]);
            if(va){na++; if(jj<ajlo)ajlo=jj; if(jj>ajhi)ajhi=jj;}
            if(vb){nb++; if(jj<bjlo)bjlo=jj; if(jj>bjhi)bjhi=jj;}
            if(va&&vb)both++;
          } }
        fprintf(stderr, "##   v=%d type=%d jband[%d..%d]: alpha-valid=%d(j%d..%d) beta-valid=%d(j%d..%d) both=%d\n",
                vv, cm->sttype[vv], jmin[vv], jmax[vv], na,ajlo,ajhi, nb,bjlo,bjhi, both);
        { int jj2=j0, jpv=jj2-jmin[vv], jpp=jj2-(i0-1); int adlo=99999,adhi=-1,bdlo=99999,bdhi=-1;
          for (dd=hd_min(cp9b, vv, jpv); dd<=hd_max(cp9b, vv, jpv) && dd<=jpp; dd++) {
            int ddp=dd-hd_min(cp9b, vv, jpv);
            if(NOT_IMPOSSIBLE(alpha[vv][jj2][ddp])){if(dd<adlo)adlo=dd;if(dd>adhi)adhi=dd;}
            if(NOT_IMPOSSIBLE(beta[vv][jj2][ddp])){if(dd<bdlo)bdlo=dd;if(dd>bdhi)bdhi=dd;} }
          fprintf(stderr, "##      @j=j0=%d band-d[%d..%d](cap jp=%d): alpha-d[%d..%d] beta-d[%d..%d]\n",
                  jj2, hd_min(cp9b, vv, jpv), hd_max(cp9b, vv, jpv), jpp, adlo,adhi, bdlo,bdhi); } } }
    cm_Fail("wedge_splitter_hb: band-infeasible subproblem (see stderr dump)");
  }
#endif

  free_banded_hb_vjd_matrix(alpha, cm, i0, j0, cp9b);
  /* brief 26_0430-227: free the banded EL deck explicitly (eldmax-driven width,
   * not the full-triangle free_vjd_deck that free_banded_hb_vjd_matrix assumes
   * for v==cm->M), then the rest of beta. */
  if ((cm->flags & CMH_LOCAL_END) && beta[cm->M] != NULL) {
    free_el_banded_vjd_deck(beta[cm->M], i0, j0, eldmax);
    beta[cm->M] = NULL;
  }
  free_banded_hb_vjd_matrix(beta,  cm, i0, j0, cp9b);
  free(eldmax); eldmax = NULL;

  if (best_v == -1) {
    v_splitter_hb(cm, dsq, L, tr, r, w, i0, best_j-best_d+1, best_j, j0, TRUE, cp9b);
    return best_sc;
  }
  if (best_v == -2) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b);
    wedge_splitter_hb(cm, dsq, L, tr, b, z, i0, j0, cp9b);
    return best_sc;
  }

  /* Usual case: banded V problem (Stage 1a.2) + banded wedge problem. */
  v_splitter_hb(cm, dsq, L, tr, r, best_v, i0, best_j-best_d+1, best_j, j0, FALSE, cp9b);
  wedge_splitter_hb(cm, dsq, L, tr, best_v, z, best_j-best_d+1, best_j, cp9b);
  return best_sc;
}

/* Function: inside_hb()
 *           EPN 2026 [brief 26_0610-007]
 *
 * Purpose:  HMM-banded (CP9) inside engine, full-deck (option-a) allocation.
 *           Mirrors inside_qdb() but enforces the CP9 j-band (jmin/jmax) and the
 *           per-(v,j) d-band (hdmin/hdmax indexed by jp_v = j-jmin[v]) in the DP
 *           loops, leaving out-of-band cells IMPOSSIBLE (full-deck init). The
 *           band arithmetic (esp. the bifurcation k-loop intersection) mirrors
 *           cm_CYKInsideAlignHB() (cm_dpalign.c ~L1102), but in platonic
 *           [v][j][d] coordinates (no jp_v/dp_v offset translation needed).
 *
 *           Unlike inside_qdb(), E states get their own deck (set to 0 only at
 *           in-band (j, d=0) cells) rather than sharing a single all-j end deck;
 *           with a j-band a shared all-j end deck would over-permit E at
 *           out-of-band j and break byte-exactness vs cm_CYKInsideAlignHB().
 */
static float
inside_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	  float ***alpha, float ****ret_alpha,
	  struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	  void ****ret_shadow,
	  int allow_begin, int *ret_b, float *ret_bsc,
	  CP9Bands_t *cp9b)
{
  int      status;
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i,k;	/* indices in sequence dimensions */
  float    sc;		/* a temporary variable holding a score */
  int      yoffset;	/* y=base+offset -- counter in child states that v can transit to */
  int      W;		/* subsequence length */
  int      jp;		/* j': relative position in the subsequence  */
  int      jp_v;        /* j index in v's band: j - jmin[v] */
  void  ***shadow;      /* shadow matrix for tracebacks */
  int    **kshad;       /* a shadow deck for bifurcations */
  char   **yshad;       /* a shadow deck for every other kind of state */
  int      b;		/* best local begin state */
  float    bsc;		/* score for using the best local begin state */
  int      jn, jx;      /* min/max j for state v in this subproblem */
  int      yy, zz;      /* bifurcation children */
  int      jp_y, jp_z;  /* j index in children's bands */
  int      kn, kx;      /* min/max k for bifurcation */
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;
  if (dpool == NULL) dpool = deckpool_create();

  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * (cm->M+1));
  for (v = 0;     v < vroot; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) touch[v] = cm->pnum[v];
  for (v = vend+1;v < cm->M; v++) touch[v] = 0;

  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* Main recursion */
  for (v = vend; v >= vroot; v--)
    {
      /* Allocate a fresh banded deck (no pooling: banded decks have per-state
       * shapes). E states get their own (banded) deck too, unlike inside_qdb. */
      alpha[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);

      if (ret_shadow != NULL && cm->sttype[v] != E_st) {
	if (cm->sttype[v] == B_st) {
	  kshad     = alloc_banded_hb_vjd_kshadow_deck(L, i0, j0, v, cp9b);
	  shadow[v] = (void **) kshad;
	} else {
	  yshad     = alloc_banded_hb_vjd_yshadow_deck(L, i0, j0, v, cp9b);
	  shadow[v] = (void **) yshad;
	}
      }

      jn = ESL_MAX(i0-1, jmin[v]);
      jx = ESL_MIN(j0,   jmax[v]);

      /* Banded IMPOSSIBLE init: every allocated (in-band) cell -> IMPOSSIBLE.
       * (Out-of-band cells don't exist; cross-state reads treat them as such.) */
      for (j = jn; j <= jx; j++) {
	jp_v = j - jmin[v];
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++)
	  alpha[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE;
      }

      if (cm->sttype[v] == E_st)
	{
	  int dpe;
	  for (j = jn; j <= jx; j++)             /* d must be 0 for E */
	    if (hb_inband(cp9b, v, j, 0, i0, j0, &dpe)) alpha[v][j][dpe] = 0.;
	}
      else if (cm->sttype[v] == D_st || cm->sttype[v] == S_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo; /* child (y+yoffset) at same (j,d); skip if out-of-band */
		    if (hb_inband(cp9b, y+yoffset, j, d, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == B_st)
	{
	  yy = cm->cfirst[v];	/* left  child */
	  zz = cm->cnum[v];	/* right child */
	  jn = ESL_MAX(jn, jmin[zz]); /* j must be in both v and z j-bands */
	  jx = ESL_MIN(jx, jmax[zz]);
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    jp_y = j - jmin[yy];
	    jp_z = j - jmin[zz];
	    kn = ESL_MAX(j - jmax[yy], hd_min(cp9b, zz, jp_z));
	    kn = ESL_MAX(kn, 0);
	    kx = ESL_MIN(jp_y, hd_max(cp9b, zz, jp_z));      /* jp_y == j - jmin[yy] */
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		for (k = kn; k <= kx; k++)
		  if ((k >= d - hd_max(cp9b, yy, jp_y-k)) && (k <= d - hd_min(cp9b, yy, jp_y-k)))
		    {
		      /* children in-band by the kn/kx + guard construction: offset directly.
		       * yy cell (j-k,d-k): row index in yy's band is jp_y-k. zz cell (j,k). */
		      int dp_yk = (d-k) - hd_min(cp9b, yy, jp_y-k);
		      int dp_zk =  k    - hd_min(cp9b, zz, jp_z);
		      if ((sc = alpha[yy][j-k][dp_yk] + alpha[zz][j][dp_zk]) > alpha[v][j][dp_v]) {
			alpha[v][j][dp_v] = sc;
			if (ret_shadow != NULL) kshad[j][dp_v] = k;
		      }
		    }
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo; /* child (y+yoffset) at (j-1,d-2) */
		    if (hb_inband(cp9b, y+yoffset, j-1, d-2, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j-1][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo; /* child (y+yoffset) at (j,d-1) */
		    if (hb_inband(cp9b, y+yoffset, j, d-1, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][dsq[i]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo; /* child (y+yoffset) at (j-1,d-1) */
		    if (hb_inband(cp9b, y+yoffset, j-1, d-1, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j-1][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][dsq[j]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	}				/* finished calculating deck v. */

      /* local begin bookkeeping. The (j0,W) cell may be out of v's band, in
       * which case it's conceptually IMPOSSIBLE and contributes nothing. */
      if (allow_begin) {
	int dpb;
	if (hb_inband(cp9b, v, j0, W, i0, j0, &dpb) &&
	    alpha[v][j0][dpb] + cm->beginsc[v] > bsc)
	  {
	    b   = v;
	    bsc = alpha[v][j0][dpb] + cm->beginsc[v];
	  }
      }
      if (allow_begin && v == 0) {
	int dpb0;
	if (hb_inband(cp9b, 0, j0, W, i0, j0, &dpb0) && bsc > alpha[0][j0][dpb0]) {
	  alpha[0][j0][dpb0] = bsc;
	  if (ret_shadow != NULL) yshad[j0][dpb0] = USED_LOCAL_BEGIN;
	}
      }

      /* reuse memory: release fully-used children (banded decks are not pooled,
       * so free them outright -- the live set is exactly what's still needed). */
      if (! do_full) {
	if (cm->sttype[v] == B_st)
	  {
	    y = cm->cfirst[v]; free_banded_hb_vjd_deck(alpha[y], i0, j0, y, cp9b); alpha[y] = NULL;
	    z = cm->cnum[v];   free_banded_hb_vjd_deck(alpha[z], i0, j0, z, cp9b); alpha[z] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) { free_banded_hb_vjd_deck(alpha[y], i0, j0, y, cp9b); alpha[y] = NULL; }
	      }
	  }
      }
  } /* end loop over all v */

  { int dpr; sc = hb_inband(cp9b, vroot, j0, W, i0, j0, &dpr) ? alpha[vroot][j0][dpr] : IMPOSSIBLE; }
  if (ret_b != NULL)   *ret_b   = b;
  if (ret_bsc != NULL) *ret_bsc = bsc;

  if (ret_alpha == NULL) {
    for (v = vroot; v <= vend; v++)
      if (alpha[v] != NULL) { free_banded_hb_vjd_deck(alpha[v], i0, j0, v, cp9b); alpha[v] = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* The deckpool holds no banded decks (we never pool them); just dispose of it. */
  if (ret_dpool == NULL) deckpool_free(dpool);
  else                   *ret_dpool = dpool;

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/*****************************************************************
 * brief 26_0430-227: banded local-end (EL) outside deck for the HMM-banded D&C.
 *
 * outside_hb()/tr_outside_hb() store the EL state (cm->M) outside scores in a
 * vjd deck that stock allocates as a FULL O(L^2/2) lower triangle every call
 * (alloc_vjd_deck()), regardless of the CP9 bands.  At genome scale this single
 * deck IS the D&C memory peak (brief 26_0430-226: sarscov2 ~1763 Mb, ~all of it
 * this deck).  But brief 26_0610-106 already proved that every write to
 * beta[cm->M][j][d] lands in a per-v-band-derived region (the v->EL feed
 * iterates only v's own banded footprint, shifted by the per-type (elsj=sdr,
 * elsd=sd) offset).  So the storage above the per-row upper d-edge is provably
 * IMPOSSIBLE and can be dropped -- exactly what cm_dpalign.c's checkpointed
 * engine already does (ckpt_el_compute_dmax / ckpt_el_deck_alloc, R-L.2).
 *
 * Here we mirror that, adapted to D&C's per-sub-call [i0,j0] window: an eldmax[]
 * array (absolute row j, -1 == empty row) bounds the deck; row j stores only
 * d=0..eldmax[j].  Unlike stock's EL deck this one is NEVER pooled (deckpool
 * reuse assumes a uniform full-triangle size; a band-dependent size would hand
 * back a wrong-sized deck) -- it is always freshly allocated and freed by its
 * owner, which also removes the aliasing that made the beta leak in
 * generic_splitter_hb() unsafe to free (see brief 26_0430-226 note there).
 *****************************************************************/

/* Allocate a banded EL vjd deck: (L+1) row pointers (rows outside [i-1..j] and
 * rows with eldmax<0 are NULL), row r holding eldmax[r]+1 floats.  Cells are
 * left uninitialized here (caller IMPOSSIBLE-inits its own band, mirroring the
 * stock init loop).  Accounts banded bytes in the cyk_dnc_track high-water. */
static float **
alloc_el_banded_vjd_deck(int L, int i, int j, int *eldmax)
{
  int     status;
  float **a;
  int     r;
  double  nb = 0.;
  ESL_ALLOC(a, sizeof(float *) * (L+1));
  for (r = 0;   r <= L; r++) a[r] = NULL;
  for (r = i-1; r <= j; r++) {
    if (eldmax[r] < 0) continue;
    ESL_ALLOC(a[r], sizeof(float) * (eldmax[r]+1));
    nb += (double) (eldmax[r]+1) * sizeof(float);
  }
  if (cyk_dnc_track) { cyk_dnc_cur_bytes += nb; cyk_dnc_note(); }
  return a;
 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
/* Free a banded EL vjd deck.  <eldmax> (may be NULL) is used only to subtract
 * the banded byte count from the cyk_dnc_track high-water; it MUST be the same
 * array passed to alloc_el_banded_vjd_deck() so alloc/free stay balanced. */
static void
free_el_banded_vjd_deck(float **a, int i, int j, int *eldmax)
{
  int r;
  if (a == NULL) return;
  if (cyk_dnc_track && eldmax != NULL) {
    double nb = 0.;
    for (r = i-1; r <= j; r++) if (eldmax[r] >= 0) nb += (double) (eldmax[r]+1) * sizeof(float);
    cyk_dnc_cur_bytes -= nb;
  }
  for (r = i-1; r <= j; r++) if (a[r] != NULL) free(a[r]);
  free(a);
}

/* Compute eldmax[0..L] (absolute row j; -1 = empty) for outside_hb()'s EL deck
 * over sub-problem (vroot..vend, [i0,j0]).  eldmax[j] is a provable superset of
 * the max d ever WRITTEN to beta[cm->M][j][.], = the union of:
 *   (a) vroot's boundary v->EL unroll   (cm_dpsmall.c ~5340-5379), and
 *   (b) each main-loop state v's banded v->EL feed (~5491-5558),
 * using the SAME per-type (elsj=sdr, elsd=sd) read-cell shifts that feed uses.
 * The extra per-cell `continue` guards inside the feed only REMOVE writes, so
 * bounding by the feed's d-loop maximum (eldhi, itself already clamped to jp)
 * is a safe upper bound.  Mirrors ckpt_el_compute_dmax() but is per-sub-call. */
static void
outside_hb_el_dmax(CM_t *cm, int L, int vroot, int vend, int i0, int j0,
		   CP9Bands_t *cp9b, int *eldmax)
{
  int *jmin = cp9b->jmin;
  int *jmax = cp9b->jmax;
  int  v, elJ, r;
  int  W = j0 - i0 + 1;
  int  w1, w2;

  for (r = 0; r <= L; r++) eldmax[r] = -1;

  /* (a) vroot boundary v->EL unroll */
  if (NOT_IMPOSSIBLE(cm->endsc[vroot])) {
    switch (cm->sttype[vroot]) {
    case MP_st:             if (W >= 2 && W-2 > eldmax[j0-1]) eldmax[j0-1] = W-2; break;
    case ML_st: case IL_st: if (W >= 1 && W-1 > eldmax[j0])   eldmax[j0]   = W-1; break;
    case MR_st: case IR_st: if (W >= 1 && W-1 > eldmax[j0-1]) eldmax[j0-1] = W-1; break;
    case S_st:  case D_st:  if (W   > eldmax[j0]) eldmax[j0] = W;   break;
    default: break;
    }
  }

  /* (b) main-loop states' banded v->EL feed.  The main loop is v = w2+1..vend,
   * where w1/w2 delimit vroot's split set (copied from outside_hb ~5315-5320). */
  w1 = cm->nodemap[cm->ndidx[vroot]];
  if (cm->sttype[vroot] == B_st) w2 = w1;
  else                           w2 = cm->cfirst[w1]-1;
  for (v = w2+1; v <= vend; v++) {
    int elsj = 0, elsd = 0, elJlo, elJhi;
    if (! NOT_IMPOSSIBLE(cm->endsc[v])) continue;
    switch (cm->sttype[v]) {
    case MP_st:                       elsj = 1; elsd = 2; break;
    case ML_st: case IL_st:           elsj = 0; elsd = 1; break;
    case MR_st: case IR_st:           elsj = 1; elsd = 1; break;
    case S_st:  case D_st: case E_st: elsj = 0; elsd = 0; break;
    default: continue;
    }
    elJlo = ESL_MAX(i0-1, jmin[v]);
    elJhi = ESL_MIN(j0,   jmax[v]);
    for (elJ = elJlo; elJ <= elJhi; elJ++) {
      int eljpv = elJ - jmin[v];
      int j     = elJ - elsj;
      int jp    = j - (i0-1);
      int eldhi, eldlo;
      if (jp < 0) continue;
      eldhi = hd_max(cp9b, v, eljpv) - elsd; if (eldhi > jp) eldhi = jp;
      eldlo = hd_min(cp9b, v, eljpv) - elsd; if (eldlo < 0)  eldlo = 0;
      if (eldhi < eldlo) continue;                 /* feed d-loop empty: no write */
      if (eldhi > eldmax[j]) eldmax[j] = eldhi;
    }
  }
}

/* brief 26_0430-227: eldmax[0..L] for tr_outside_hb()'s three EL decks
 * (plane 0=J beta[cm->M], 1=L betaL[cm->M], 2=R betaR[cm->M]).  Same per-sub-call
 * structure as outside_hb_el_dmax(), with per-plane v->EL feed shifts (elsj,elsd)
 * and per-plane feed gating (L feed only Lvalid[v], R feed only Rvalid[v]) that
 * mirror the feed loops at cm_dpsmall.c ~7891 (J) / ~7972 (L) / ~8044 (R).  The
 * vroot boundary v->EL seeds touch only rows j0 and j0-1, and are covered
 * CONSERVATIVELY (full width on those two rows) when the boundary is active for
 * this plane -- a safe superset (a v-feed write at row j0/j0-1 also has d<=jp, so
 * full width on those two rows dominates every write there; the whole O(L^2) win
 * lives in the interior v-feed rows).  Boundary is inactive when vroot==0 (the
 * truncated-begin top-level call), so the biggest deck stays precisely banded.
 * This helper is deterministic in (cm,cp9b,vroot,vend,i0,j0,plane), so the
 * writer (tr_outside_hb) and the reader (tr_*_splitter_hb) recompute the SAME
 * eldmax independently -- no need to thread it through the signature. */
static void
tr_outside_hb_el_dmax(CM_t *cm, int L, int vroot, int vend, int i0, int j0,
		      int plane, CP9Bands_t *cp9b, int *eldmax)
{
  int *jmin = cp9b->jmin;
  int *jmax = cp9b->jmax;
  int  v, elJ, r;
  int  W = j0 - i0 + 1;
  int  w1, w2;

  for (r = 0; r <= L; r++) eldmax[r] = -1;

  /* (a) vroot boundary v->EL seeds (rows j0, j0-1 only), conservative */
  if (vroot != 0 && NOT_IMPOSSIBLE(cm->endsc[vroot])) {
    int active = (plane == 0) ||
                 (plane == 1 && cp9b->Lvalid[vroot]) ||
                 (plane == 2 && cp9b->Rvalid[vroot]);
    if (active) {
      if (W   > eldmax[j0])   eldmax[j0]   = W;      /* boundary d up to W   (S/D) */
      if (j0-1 >= 0 && W-1 > eldmax[j0-1]) eldmax[j0-1] = W-1;  /* d up to W-1     */
    }
  }

  /* (b) main-loop states' banded v->EL feed, per-plane */
  w1 = cm->nodemap[cm->ndidx[vroot]];
  if (cm->sttype[vroot] == B_st) w2 = w1;
  else                           w2 = cm->cfirst[w1]-1;
  for (v = w2+1; v <= vend; v++) {
    int elsj = 0, elsd = 0, elJlo, elJhi;
    if (! NOT_IMPOSSIBLE(cm->endsc[v])) continue;
    if (plane == 1 && ! cp9b->Lvalid[v]) continue;   /* L feed gated on Lvalid[v] */
    if (plane == 2 && ! cp9b->Rvalid[v]) continue;   /* R feed gated on Rvalid[v] */
    switch (cm->sttype[v]) {
    case MP_st:             if (plane==0) { elsj=1; elsd=2; } else if (plane==1) { elsj=0; elsd=1; } else { elsj=1; elsd=1; } break;
    case ML_st: case IL_st: if (plane==0) { elsj=0; elsd=1; } else if (plane==1) { elsj=0; elsd=1; } else { elsj=0; elsd=0; } break;
    case MR_st: case IR_st: if (plane==0) { elsj=1; elsd=1; } else if (plane==1) { elsj=0; elsd=0; } else { elsj=1; elsd=1; } break;
    case S_st:  case D_st: case E_st: elsj = 0; elsd = 0; break;
    default: continue;
    }
    elJlo = ESL_MAX(i0-1, jmin[v]);
    elJhi = ESL_MIN(j0,   jmax[v]);
    for (elJ = elJlo; elJ <= elJhi; elJ++) {
      int eljpv = elJ - jmin[v];
      int j     = elJ - elsj;
      int jp    = j - (i0-1);
      int eldhi, eldlo;
      if (jp < 0) continue;
      eldhi = hd_max(cp9b, v, eljpv) - elsd; if (eldhi > jp) eldhi = jp;
      eldlo = hd_min(cp9b, v, eljpv) - elsd; if (eldlo < 0)  eldlo = 0;
      if (eldhi < eldlo) continue;
      if (eldhi > eldmax[j]) eldmax[j] = eldhi;
    }
  }
}

/* Function: outside_hb()
 *           EPN 2026 [brief 26_0610-007]
 *
 * Purpose:  HMM-banded (CP9) outside engine, full-deck (option-a) allocation.
 *           Mirrors outside_qdb() but computes beta[v][j][d] only at v's in-band
 *           (j,d) cells; out-of-band cells stay IMPOSSIBLE (full-deck init).
 *           Parent contributions read full decks, so out-of-band parent cells
 *           return IMPOSSIBLE and don't contribute. The v->EL section loops the
 *           full d range (relying on out-of-band v cells being IMPOSSIBLE).
 */
static void
outside_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
	   int do_full, float ***beta, float ****ret_beta,
	   struct deckpool_s *dpool, struct deckpool_s **ret_dpool, CP9Bands_t *cp9b,
	   int **ret_eldmax)
{
  int      status;
  int     *eldmax = NULL;       /* brief 26_0430-227: per-row upper d-edge of the banded EL deck */
  int      v,y;			/* indices for states */
  int      j,d,i;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */
  int      W;			/* subsequence length */
  int      jp;			/* j': relative position in the subsequence, 0..W */
  int      jp_v;                /* j - jmin[v] */
  int      voffset;		/* index of v in t_v(y) transition scores */
  int      w1,w2;		/* bounds of split set */
  int      jn, jx;              /* min/max j for state v in this subproblem */
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  W = j0-i0+1;
  if (dpool == NULL) dpool = deckpool_create();

  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck (and the whole split set it lives in). */
  w1 = cm->nodemap[cm->ndidx[vroot]];
  if (cm->sttype[vroot] == B_st) {
    w2 = w1;
    if (vend != vroot) cm_Fail("oh no. not again.");
  } else
    w2 = cm->cfirst[w1]-1;

  for (v = w1; v <= w2; v++) {
    beta[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
    for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
      jp_v = j - jmin[v];
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) beta[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE;
    }
  }
  { int dpr; if (hb_inband(cp9b, vroot, j0, W, i0, j0, &dpr)) beta[vroot][j0][dpr] = 0; }

  /* Initialize the EL deck at M, if local ends.  brief 26_0430-227: allocate it
   * BANDED on the upper d-edge (eldmax[j]) instead of as a full O(W^2) triangle,
   * and never pool it (a band-dependent size is incompatible with the shared
   * deckpool's uniform-size reuse).  Only the banded (d=0..eldmax[j]) cells are
   * IMPOSSIBLE-inited; cells above eldmax[j] are provably never written/read. */
  if (cm->flags & CMH_LOCAL_END) {
    ESL_ALLOC(eldmax, sizeof(int) * (L+1));
    outside_hb_el_dmax(cm, L, vroot, vend, i0, j0, cp9b, eldmax);
    beta[cm->M] = alloc_el_banded_vjd_deck(L, i0, j0, eldmax);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      if (eldmax[j] < 0) continue;
      for (d = 0; d <= eldmax[j]; d++) beta[cm->M][j][d] = IMPOSSIBLE;
    }
    /* vroot -> EL boundary unroll (no band on EL). */
    if (NOT_IMPOSSIBLE(cm->endsc[vroot])) {
      switch (cm->sttype[vroot]) {
      case MP_st:
	if (W < 2) break;
	if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
	else
	  escore = DegeneratePairScore(cm->abc, cm->esc[vroot], dsq[i0], dsq[j0]);
	beta[cm->M][j0-1][W-2] = cm->endsc[vroot] + (cm->el_selfsc * (W-2)) + escore;
	if (beta[cm->M][j0-1][W-2] < IMPOSSIBLE) beta[cm->M][j0-1][W-2] = IMPOSSIBLE;
	break;
      case ML_st:
      case IL_st:
	if (W < 1) break;
	if (dsq[i0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) dsq[i0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[vroot]);
	beta[cm->M][j0][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + escore;
	if (beta[cm->M][j0][W-1] < IMPOSSIBLE) beta[cm->M][j0][W-1] = IMPOSSIBLE;
	break;
      case MR_st:
      case IR_st:
	if (W < 1) break;
	if (dsq[j0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) dsq[j0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[vroot]);
	beta[cm->M][j0-1][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + escore;
	if (beta[cm->M][j0-1][W-1] < IMPOSSIBLE) beta[cm->M][j0-1][W-1] = IMPOSSIBLE;
	break;
      case S_st:
      case D_st:
	beta[cm->M][j0][W] = cm->endsc[vroot] + (cm->el_selfsc * W);
	if (beta[cm->M][j0][W] < IMPOSSIBLE) beta[cm->M][j0][W] = IMPOSSIBLE;
	break;
      case B_st:
      default: cm_Fail("bogus parent state %d\n", cm->sttype[vroot]);
      }
    }
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;      v < w1; v++) touch[v] = 0;
  for (v = vend+1; v < cm->M; v++) touch[v] = 0;
  for (v = w1; v <= vend; v++) {
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
  }

  /* Main loop down through the decks */
  for (v = w2+1; v <= vend; v++)
    {
      beta[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);

      /* banded IMPOSSIBLE init (every allocated in-band cell) */
      for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
	jp_v = j - jmin[v];
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) beta[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE;
      }

      /* local begin into v, if the (j0,W) cell is in v's band */
      if ((vroot == 0 && i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN))
	  && (jmin[v] <= j0 && jmax[v] >= j0)
	  && (hd_min(cp9b, v, j0-jmin[v]) <= W && hd_max(cp9b, v, j0-jmin[v]) >= W))
	beta[v][j0][W - hd_min(cp9b, v, j0-jmin[v])] = cm->beginsc[v];

      /* main recursion: only v's in-band (j,d) cells.
       * j and d are iterated in DECREASING order: insert (IL/IR) self-transitions
       * make beta[v][j][d] depend on beta[v][j+1][d+1] (IR) or beta[v][j][d+1] (IL),
       * so the larger-(j,d) cells must be computed first (mirrors outside_qdb). */
      jn = ESL_MAX(i0-1, jmin[v]);
      jx = ESL_MIN(j0,   jmax[v]);
      for (j = jx; j >= jn; j--) {
	jp   = j - (i0-1);
	jp_v = j - jmin[v];
	for (d = ESL_MIN(hd_max(cp9b, v, jp_v), jp); d >= hd_min(cp9b, v, jp_v); d--)
	  {
	    int dp_v = d - hd_min(cp9b, v, jp_v);   /* v,j,d in band by loop bounds */
	    int dp_y;                        /* parent cell offset (when in-band) */
	    i = j-d+1;
	    for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	      if (y < vroot) continue;
	      voffset = v - cm->cfirst[y];

	      switch(cm->sttype[y]) {
	      case MP_st:
		if (j == j0 || d == jp) continue;
		if (! hb_inband(cp9b, y, j+1, d+2, i0, j0, &dp_y)) continue;
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
		if ((sc = beta[y][j+1][dp_y] + cm->tsc[y][voffset] + escore) > beta[v][j][dp_v])
		  beta[v][j][dp_v] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;
		if (! hb_inband(cp9b, y, j, d+1, i0, j0, &dp_y)) continue;
		if (dsq[i-1] < cm->abc->K)
		  escore = cm->esc[y][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		if ((sc = beta[y][j][dp_y] + cm->tsc[y][voffset] + escore) > beta[v][j][dp_v])
		  beta[v][j][dp_v] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (! hb_inband(cp9b, y, j+1, d+1, i0, j0, &dp_y)) continue;
		if (dsq[j+1] < cm->abc->K)
		  escore = cm->esc[y][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);
		if ((sc = beta[y][j+1][dp_y] + cm->tsc[y][voffset] + escore) > beta[v][j][dp_v])
		  beta[v][j][dp_v] = sc;
		break;
	      case S_st:
	      case E_st:
	      case D_st:
		if (! hb_inband(cp9b, y, j, d, i0, j0, &dp_y)) continue;
		if ((sc = beta[y][j][dp_y] + cm->tsc[y][voffset]) > beta[v][j][dp_v])
		  beta[v][j][dp_v] = sc;
		break;
	      default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
	      }
	    } /* end loop over parents */
	    if (beta[v][j][dp_v] < IMPOSSIBLE) beta[v][j][dp_v] = IMPOSSIBLE;
	  } /* end loop over d */
      } /* end loop over j */

      /* v->EL local end transitions (EL = deck M, full/unbanded). The beta[v]
       * source is banded: read it through hb_inband (out-of-band -> IMPOSSIBLE,
       * skip). beta[cm->M] is a full deck, so its [j][d] index is unoffset.
       *
       * brief 26_0610-106: band-limit the v->EL feed to v's OWN banded footprint,
       * porting brief 26_0610-101's fix (tr_outside_hb's J-feed) to this
       * non-truncated sibling. The old loop swept the full 0..W x 0..jp triangle
       * (O(W^2)) for EVERY state v with a local end -- confirmed by measurement
       * (Task A) to be 99.94% of all EL-deck work here, the dominant cost of
       * local-mode D&C at genome scale. But a cell (j,d) can only update
       * beta[cm->M] when v's read cell -- shifted by the per-type (elsj,elsd)
       * below -- is IN v's band; every out-of-band read hit `!hb_inband ->
       * continue` and did nothing. So iterate ONLY v's band rows (elJ = read
       * row) and their shifted d-range. The body (switch, hb_inband guard,
       * boundary gates, all index math) is BYTE-IDENTICAL to the old sweep; the
       * tighter bounds are a provable superset of the old update set, so the
       * result is byte-exact. Read-cell shifts (identical to 101's J-feed
       * table, since this is the same J-only recursion): MP=(j+1,d+2)
       * ML/IL=(j,d+1) MR/IR=(j+1,d+1) S/D/E=(j,d). */
      if (NOT_IMPOSSIBLE(cm->endsc[v])) {
	int dp_v;
	int elsj = 0, elsd = 0, elJ, elJlo, elJhi;
	switch (cm->sttype[v]) {
	case MP_st:                       elsj = 1; elsd = 2; break;
	case ML_st: case IL_st:           elsj = 0; elsd = 1; break;
	case MR_st: case IR_st:           elsj = 1; elsd = 1; break;
	case S_st:  case D_st: case E_st: elsj = 0; elsd = 0; break;
	case B_st:
	default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	}
	elJlo = ESL_MAX(i0-1, jmin[v]); elJhi = ESL_MIN(j0, jmax[v]);
	for (elJ = elJlo; elJ <= elJhi; elJ++) {
	  int eljpv = elJ - jmin[v], eldlo, eldhi;
	  j  = elJ - elsj;
	  jp = j - (i0-1);
	  if (jp < 0) continue;
	  eldlo = hd_min(cp9b, v, eljpv) - elsd; if (eldlo < 0)  eldlo = 0;
	  eldhi = hd_max(cp9b, v, eljpv) - elsd; if (eldhi > jp) eldhi = jp;
	  for (d = eldlo; d <= eldhi; d++)
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st:
		if (j == j0 || d == jp) continue;
		if (! hb_inband(cp9b, v, j+1, d+2, i0, j0, &dp_v)) continue;
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		if ((sc = beta[v][j+1][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;
		if (! hb_inband(cp9b, v, j, d+1, i0, j0, &dp_v)) continue;
		if (dsq[i-1] < cm->abc->K)
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		if ((sc = beta[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (! hb_inband(cp9b, v, j+1, d+1, i0, j0, &dp_v)) continue;
		if (dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		if ((sc = beta[v][j+1][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if (! hb_inband(cp9b, v, j, d, i0, j0, &dp_v)) continue;
		if ((sc = beta[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d)) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case B_st:
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	      } /* end switch */
	    } /* end loop over d */
	} /* end loop over jp */
      } /* end v->EL */

      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { free_banded_hb_vjd_deck(beta[y], i0, j0, y, cp9b); beta[y] = NULL; }
	}
      }
    } /* end loop over decks v. */

  if (ret_beta == NULL) {
    for (v = w1; v <= vend; v++)
      if (beta[v] != NULL) { free_banded_hb_vjd_deck(beta[v], i0, j0, v, cp9b); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
      free_el_banded_vjd_deck(beta[cm->M], i0, j0, eldmax);   /* brief 26_0430-227: banded EL deck */
      beta[cm->M] = NULL;
    }
    free(beta);
    free(eldmax); eldmax = NULL;   /* deck freed; nothing left for eldmax to describe */
  } else *ret_beta = beta;

  /* brief 26_0430-227: hand the banded EL deck's per-row d-band edges to the
   * caller (it owns the deck when ret_beta!=NULL, and needs eldmax to clamp its
   * beta[cm->M] reads and to free the deck).  If the caller took the deck but
   * doesn't want eldmax, or there is no EL deck, free/pass NULL as appropriate. */
  if (ret_eldmax != NULL) *ret_eldmax = eldmax;
  else                    free(eldmax);

  /* The deckpool holds no banded decks (we never pool them); just dispose of it. */
  if (ret_dpool == NULL) deckpool_free(dpool);
  else                   *ret_dpool = dpool;
  free(touch);
  return;
 ERROR:
  cm_Fail("Memory allocation error.");
}

/* Function: insideT_hb()
 *           EPN 2026 [brief 26_0610-007]
 *
 * Purpose:  HMM-banded analogue of insideT(): run inside_hb() with a shadow
 *           matrix, then trace back. The traceback loop is band-agnostic (it
 *           just follows shadow pointers), so it is identical to insideT().
 */
static float
insideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
	   int r, int z, int i0, int j0, int allow_begin, CP9Bands_t *cp9b)
{
  int       status;
  void   ***shadow;
  float     sc;
  ESL_STACK *pda;
  int       v,j,d,i;
  int       k;
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;

  sc = inside_hb(cm, dsq, L, r, z, i0, j0,
		 BE_EFFICIENT,
		 NULL, NULL,
		 NULL, NULL,
		 &shadow,
		 allow_begin,
		 &b, &bsc,
		 cp9b);

  pda = esl_stack_ICreate();
  if(pda == NULL) goto ERROR;
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  while (1) {
    if (cm->sttype[v] == B_st) {
      /* shadow decks are banded: (v,j,d) is on the optimal path, hence in-band */
      k = ((int **) shadow[v])[j][d - hd_min(cp9b, v, j - cp9b->jmin[v])];
      if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;
      if((status = esl_stack_IPush(pda, k)) != eslOK) goto ERROR;
      if((status = esl_stack_IPush(pda, tr->n-1)) != eslOK) goto ERROR;
      j = j-k;
      d = d-k;
      i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];
      y = cm->cnum[v];
      i = j-d+1;
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
    } else {
      yoffset = ((char **) shadow[v])[j][d - hd_min(cp9b, v, j - cp9b->jmin[v])];
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL)
	{
	  InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
	  v = cm->M;
	}
      else if (yoffset == USED_LOCAL_BEGIN)
	{
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
  esl_stack_Destroy(pda);
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* NEVERREACHED */
}

/*################################################################
 * Stage 1a.2 (brief 26_0610-009): banded V-problem (class-2 vji) engines.
 *
 * These are the banded analogues of the EXACT vinside()/voutside()/vinsideT()/
 * v_splitter(), NOT ports of the hazardous QDB vinside_qdb()/.../v_splitter_qdb().
 * They mirror the exact engines' structure and tie-breaking line-for-line (so the
 * non-clip optimum is byte-identical to Stage 1b, which used the exact engines),
 * but (a) store/iterate ONLY the in-band cells of each vji deck (offset-indexed
 * via vji_inband()), and (b) design out the QDB negative-malloc crash: the deck
 * allocator never makes a non-positive row, and v_splitter_hb has an EXPLICIT
 * empty-band path that terminates cleanly instead of recursing on degenerate
 * corners. The bands are the SAME per-(v,j) CP9 bands cm_CYKInsideAlignHB() uses,
 * so the whole D&C now explores exactly the oracle's banded space -- the basis
 * for byte-exactness vs the oracle on CLIP cases (the Stage 1a.2 gate).
 *################################################################*/

/* Function: vinside_hb()  [brief 26_0610-009] -- banded analogue of vinside(). */
static float
vinside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
	   int r, int z, int i0, int i1, int j1, int j0, int useEL,
	   int do_full, float ***a, float ****ret_a,
	   char ****ret_shadow,
	   int allow_begin, int *ret_b, float *ret_bsc, CP9Bands_t *cp9b)
{
  int      status;
  char  ***shadow = NULL;
  int      v, i, j;
  int      w1, w2;
  int      jp, op, op_y;
  int     *touch;
  int      y, yoffset;
  float    sc, val;
  int      b;
  float    bsc;
  int      jpb, ilo, ihi;
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  b   = -1;
  bsc = IMPOSSIBLE;
  if (cyk_dnc_track) cyk_dnc_vji_row_floats = i1 - i0 + 1; /* for the full EL deck's free */

  if (a == NULL) {
    ESL_ALLOC(a, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) a[v] = NULL;
  }

  /* the whole split set w1<=v<=w2 (z's node) is allocated + IMPOSSIBLE-init */
  w1 = cm->nodemap[cm->ndidx[z]];
  w2 = cm->cfirst[w1]-1;
  for (v = w1; v <= w2; v++) {
    a[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b);
    banded_hb_vji_init_impossible(a[v], i0, i1, j1, j0, v, cp9b);
  }

  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(char **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* Boundary condition: seed the one non-IMPOSSIBLE cell at z. Global: connect z
   * at (j1,i1). Local (useEL): connect z to EL (unroll one recursion step). Every
   * write is guarded -- an out-of-band seed means the V-problem is unreachable. */
  if (! useEL) {
    if (vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) a[z][0][op] = 0.;
  } else {
    if (ret_shadow != NULL) shadow[z] = alloc_banded_hb_vji_shadow_deck(i0,i1,j1,j0,z,cp9b);
    switch (cm->sttype[z]) {
    case D_st:
    case S_st:
      if (vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) {
	a[z][0][op] = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (ret_shadow != NULL) shadow[z][0][op] = USED_EL;
      }
      break;
    case MP_st:
      if (i0 == i1 || j1 == j0) break;
      if (vji_inband(cp9b, z, j1+1, i1-1, i0,i1,j1,j0, &op)) {
	val = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (dsq[i1-1] < cm->abc->K && dsq[j1+1] < cm->abc->K)
	  val += cm->esc[z][(int) (dsq[i1-1]*cm->abc->K+dsq[j1+1])];
	else
	  val += DegeneratePairScore(cm->abc, cm->esc[z], dsq[i1-1], dsq[j1+1]);
	if (val < IMPOSSIBLE) val = IMPOSSIBLE;
	a[z][1][op] = val;
	if (ret_shadow != NULL) shadow[z][1][op] = USED_EL;
      }
      break;
    case ML_st:
    case IL_st:
      if (i0 == i1) break;
      if (vji_inband(cp9b, z, j1, i1-1, i0,i1,j1,j0, &op)) {
	val = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (dsq[i1-1] < cm->abc->K)
	  val += cm->esc[z][(int) dsq[i1-1]];
	else
	  val += esl_abc_FAvgScore(cm->abc, dsq[i1-1], cm->esc[z]);
	if (val < IMPOSSIBLE) val = IMPOSSIBLE;
	a[z][0][op] = val;
	if (ret_shadow != NULL) shadow[z][0][op] = USED_EL;
      }
      break;
    case MR_st:
    case IR_st:
      if (j1 == j0) break;
      if (vji_inband(cp9b, z, j1+1, i1, i0,i1,j1,j0, &op)) {
	val = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (dsq[j1+1] < cm->abc->K)
	  val += cm->esc[z][(int) dsq[j1+1]];
	else
	  val += esl_abc_FAvgScore(cm->abc, dsq[j1+1], cm->esc[z]);
	if (val < IMPOSSIBLE) val = IMPOSSIBLE;
	a[z][1][op] = val;
	if (ret_shadow != NULL) shadow[z][1][op] = USED_EL;
      }
      break;
    }
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;    v < r;      v++) touch[v] = 0;
  for (v = r;    v <= w2;    v++) touch[v] = cm->pnum[v]; /* to bottom of split set */
  for (v = w2+1; v < cm->M;  v++) touch[v] = 0;

  /* Special case: begin transition straight into z on empty subsequences. */
  if (allow_begin && j0-j1 == 0 && i1-i0 == 0) {
    if (vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) {
      b   = z;
      bsc = a[z][0][op] + cm->beginsc[z];
      if (z == 0) {
	a[0][0][op] = bsc;
	if (ret_shadow != NULL) shadow[0][0][op] = USED_LOCAL_BEGIN;
      }
    }
  }

  /* Main recursion: v = w1-1 downto r. Only in-band cells of each deck. */
  for (v = w1-1; v >= r; v--)
    {
      a[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b);
      banded_hb_vji_init_impossible(a[v], i0, i1, j1, j0, v, cp9b);
      if (ret_shadow != NULL)
	shadow[v] = alloc_banded_hb_vji_shadow_deck(i0, i1, j1, j0, v, cp9b);

      if (cm->sttype[v] == E_st || cm->sttype[v] == B_st || (cm->sttype[v] == S_st && v > r))
	cm_Fail("you told me you wouldn't ever do that again.");

      for (j = ESL_MAX(j1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
	jp  = j - j1;
	jpb = j - jmin[v];
	ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
	ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
	/* i DECREASING (matches exact vinside): IL self-transitions read child
	 * (j,i+1), so the larger-i cell must be computed first. */
	for (i = ihi; i >= ilo; i--)
	  {
	    int d = j - i + 1;
	    op = i - ilo;
	    y  = cm->cfirst[v];

	    if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) {
	      /* child y at (j,i) */
	      if (vji_inband(cp9b, y, j, i, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j, i, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }
	    else if (cm->sttype[v] == MP_st) {
	      /* child y at (j-1,i+1) */
	      if (vji_inband(cp9b, y, j-1, i+1, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp-1][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j-1, i+1, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp-1][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		a[v][jp][op] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
	      else
		a[v][jp][op] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }
	    else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
	      /* child y at (j,i+1) */
	      if (vji_inband(cp9b, y, j, i+1, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j, i+1, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (dsq[i] < cm->abc->K)
		a[v][jp][op] += cm->esc[v][dsq[i]];
	      else
		a[v][jp][op] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }
	    else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
	      /* child y at (j-1,i) */
	      if (vji_inband(cp9b, y, j-1, i, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp-1][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j-1, i, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp-1][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (dsq[j] < cm->abc->K)
		a[v][jp][op] += cm->esc[v][dsq[j]];
	      else
		a[v][jp][op] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }
	  } /* end loop over i */
      } /* end loop over j; finished deck v */

      /* local begin bookkeeping (cell (j0,i0)) */
      if (allow_begin && vji_inband(cp9b, v, j0, i0, i0,i1,j1,j0, &op) &&
	  a[v][j0-j1][op] + cm->beginsc[v] > bsc) {
	b   = v;
	bsc = a[v][j0-j1][op] + cm->beginsc[v];
      }
      if (allow_begin && v == 0 && vji_inband(cp9b, 0, j0, i0, i0,i1,j1,j0, &op) &&
	  bsc > a[0][j0-j1][op]) {
	a[0][j0-j1][op] = bsc;
	if (ret_shadow != NULL) shadow[v][j0-j1][op] = USED_LOCAL_BEGIN;
      }

      /* reuse memory: free children whose touch count hits 0 (banded: no pool). */
      if (! do_full) {
	for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++) {
	  touch[y]--;
	  if (touch[y] == 0) { free_banded_hb_vji_deck(a[y], i0,i1,j1,j0, y, cp9b); a[y] = NULL; }
	}
      }
    } /* end loop over v */

  { int dpr; sc = vji_inband(cp9b, r, j0, i0, i0,i1,j1,j0, &dpr) ? a[r][j0-j1][dpr] : IMPOSSIBLE; }
  if (ret_b   != NULL) *ret_b   = b;
  if (ret_bsc != NULL) *ret_bsc = bsc;

  if (ret_a == NULL) {
    for (v = r; v <= w2; v++)
      if (a[v] != NULL) { free_banded_hb_vji_deck(a[v], i0,i1,j1,j0, v, cp9b); a[v] = NULL; }
    free(a);
  } else *ret_a = a;

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.; /* never reached */
}

/* Function: voutside_hb()  [brief 26_0610-009] -- banded analogue of voutside(). */
static void
voutside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
	    int r, int z, int i0, int i1, int j1, int j0, int useEL,
	    int do_full, float ***beta, float ****ret_beta, CP9Bands_t *cp9b)
{
  int      status;
  int      v, y;
  int      i, j;
  int      jp, ip, op, op_y;
  float    sc;
  int     *touch;
  float    escore;
  int      voffset;
  int      jpb, ilo, ihi;
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  if (cyk_dnc_track) cyk_dnc_vji_row_floats = i1 - i0 + 1; /* for the full EL deck's free */

  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) beta[v] = NULL;
  }

  /* root deck */
  beta[r] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, r, cp9b);
  banded_hb_vji_init_impossible(beta[r], i0, i1, j1, j0, r, cp9b);
  if (vji_inband(cp9b, r, j0, i0, i0,i1,j1,j0, &op)) beta[r][j0-j1][op] = 0;

  /* EL deck: FULL/unbanded ("no band on EL"). */
  if (useEL && cm->flags & CMH_LOCAL_END) {
    beta[cm->M] = alloc_vji_deck(i0, i1, j1, j0);
    for (jp = 0; jp <= j0-j1; jp++)
      for (ip = 0; ip <= i1-i0; ip++)
	beta[cm->M][jp][ip] = IMPOSSIBLE;
  }
  if (useEL && NOT_IMPOSSIBLE(cm->endsc[r])) {
    switch (cm->sttype[r]) {
    case MP_st:
      if (i0 == i1 || j1 == j0) break;
      if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	escore = cm->esc[r][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
      else
	escore = DegeneratePairScore(cm->abc, cm->esc[r], dsq[i0], dsq[j0]);
      beta[cm->M][j0-j1-1][1] = cm->endsc[r] + (cm->el_selfsc * ((j0-1)-(i0+1)+1)) + escore;
      break;
    case ML_st:
    case IL_st:
      if (i0 == i1) break;
      if (dsq[i0] < cm->abc->K) escore = cm->esc[r][(int) dsq[i0]];
      else                      escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[r]);
      beta[cm->M][j0-j1][1] = cm->endsc[r] + (cm->el_selfsc * ((j0)-(i0+1)+1)) + escore;
      break;
    case MR_st:
    case IR_st:
      if (j0 == j1) break;
      if (dsq[j0] < cm->abc->K) escore = cm->esc[r][(int) dsq[j0]];
      else                      escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[r]);
      beta[cm->M][j0-j1-1][0] = cm->endsc[r] + (cm->el_selfsc * ((j0-1)-(i0)+1)) + escore;
      break;
    case S_st:
    case D_st:
      beta[cm->M][j0-j1][0] = cm->endsc[r] + (cm->el_selfsc * ((j0)-(i0)+1));
      break;
    default:  cm_Fail("bogus parent state %d\n", cm->sttype[r]);
    }
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;   v < r;     v++) touch[v] = 0;
  for (v = z+1; v < cm->M; v++) touch[v] = 0;
  for (v = r;   v <= z;    v++) {
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
  }

  /* Main loop down through the decks. */
  for (v = r+1; v <= z; v++)
    {
      beta[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b);
      banded_hb_vji_init_impossible(beta[v], i0, i1, j1, j0, v, cp9b);

      /* local begin into v at (j0,i0) */
      if (r == 0 && i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN)) {
	if (vji_inband(cp9b, v, j0, i0, i0,i1,j1,j0, &op) && cm->beginsc[v] > beta[v][j0-j1][op])
	  beta[v][j0-j1][op] = cm->beginsc[v];
      }

      /* main recursion: in-band cells, j DECREASING, i INCREASING (insert
       * self-transitions: beta[v][j][i] depends on beta[v][j][i-1] (IL) /
       * beta[v][j+1][i] (IR), so larger-(j),smaller-(i) cells must precede). */
      for (j = ESL_MIN(j0, jmax[v]); j >= ESL_MAX(j1, jmin[v]); j--) {
	jp  = j - j1;
	jpb = j - jmin[v];
	ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
	ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
	for (i = ilo; i <= ihi; i++) {
	  op = i - ilo;
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    if (y < r) continue;
	    voffset = v - cm->cfirst[y];
	    switch (cm->sttype[y]) {
	    case MP_st:
	      if (j == j0 || i == i0) continue;
	      if (! vji_inband(cp9b, y, j+1, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
	      else
		escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
	      if ((sc = beta[y][jp+1][op_y] + cm->tsc[y][voffset] + escore) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    case ML_st:
	    case IL_st:
	      if (i == i0) continue;
	      if (! vji_inband(cp9b, y, j, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K) escore = cm->esc[y][(int) dsq[i-1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
	      if ((sc = beta[y][jp][op_y] + cm->tsc[y][voffset] + escore) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    case MR_st:
	    case IR_st:
	      if (j == j0) continue;
	      if (! vji_inband(cp9b, y, j+1, i, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[j+1] < cm->abc->K) escore = cm->esc[y][(int) dsq[j+1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);
	      if ((sc = beta[y][jp+1][op_y] + cm->tsc[y][voffset] + escore) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    case S_st:
	    case E_st:
	    case D_st:
	      if (! vji_inband(cp9b, y, j, i, i0,i1,j1,j0, &op_y)) continue;
	      if ((sc = beta[y][jp][op_y] + cm->tsc[y][voffset]) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    default: cm_Fail("bogus parent state %d\n", cm->sttype[y]);
	    }
	  } /* end loop over parents */
	  if (beta[v][jp][op] < IMPOSSIBLE) beta[v][jp][op] = IMPOSSIBLE;
	} /* end loop over i */
      } /* end loop over j */

      /* v->EL local-end transitions (EL deck M is full/unbanded; read beta[v]
       * banded via vji_inband, write beta[M] at the full [jp][ip] index). */
      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (jp = j0-j1; jp >= 0; jp--) {
	  j = jp + j1;
	  for (ip = 0; ip <= i1-i0; ip++) {
	    i = ip + i0;
	    switch (cm->sttype[v]) {
	    case MP_st:
	      if (j == j0 || i == i0) continue;
	      if (! vji_inband(cp9b, v, j+1, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
	      else
		escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
	      if ((sc = beta[v][jp+1][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    case ML_st:
	    case IL_st:
	      if (i == i0) continue;
	      if (! vji_inband(cp9b, v, j, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K) escore = cm->esc[v][(int) dsq[i-1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
	      if ((sc = beta[v][jp][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    case MR_st:
	    case IR_st:
	      if (j == j0) continue;
	      if (! vji_inband(cp9b, v, j+1, i, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[j+1] < cm->abc->K) escore = cm->esc[v][(int) dsq[j+1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
	      if ((sc = beta[v][jp+1][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    case S_st:
	    case D_st:
	    case E_st:
	      if (! vji_inband(cp9b, v, j, i, i0,i1,j1,j0, &op_y)) continue;
	      if ((sc = beta[v][jp][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1))) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    default:  cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	    }
	  }
	}
      }

      /* reuse memory: push parents no longer needed (banded: free outright). */
      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { free_banded_hb_vji_deck(beta[y], i0,i1,j1,j0, y, cp9b); beta[y] = NULL; }
	}
      }
    } /* end loop over decks v */

  if (ret_beta == NULL) {
    for (v = r; v <= z; v++)
      if (beta[v] != NULL) { free_banded_hb_vji_deck(beta[v], i0,i1,j1,j0, v, cp9b); beta[v] = NULL; }
    if (useEL && cm->flags & CMH_LOCAL_END) { free_vji_deck(beta[cm->M], j1, j0); beta[cm->M] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  free(touch);
  return;

 ERROR:
  cm_Fail("Memory allocation error.\n");
}

/* Function: vinsideT_hb()  [brief 26_0610-009] -- banded analogue of vinsideT(). */
static float
vinsideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
	    int r, int z, int i0, int i1, int j1, int j0, int useEL,
	    int allow_begin, CP9Bands_t *cp9b)
{
  char ***shadow;
  float   sc;
  int     v, y;
  int     j, i;
  int     op;
  int     yoffset;
  int     b;
  float   bsc;

  if (r == z) {
    InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, r);
    return 0.;
  }

  sc = vinside_hb(cm, dsq, L, r, z, i0, i1, j1, j0, useEL,
		  BE_EFFICIENT, NULL, NULL, &shadow,
		  allow_begin, &b, &bsc, cp9b);

  v = r;
  j = j0;
  i = i0;
  while (1) {
    /* Defensive (brief 26_0610-003 hazard): if the on-path cell fell out of band (empty/
     * collapsed V-problem), terminate cleanly by attaching z -- never crash. */
    if (! vji_inband(cp9b, v, j, i, i0,i1,j1,j0, &op)) {
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, z);
      break;
    }
    yoffset = shadow[v][j-j1][op];

    switch (cm->sttype[v]) {
    case D_st:            break;
    case MP_st: i++; j--; break;
    case ML_st: i++;      break;
    case MR_st:      j--; break;
    case IL_st: i++;      break;
    case IR_st:      j--; break;
    case S_st:            break;
    default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
    }

    if (yoffset == USED_EL) {
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M);
      break;
    }
    else if (yoffset == USED_LOCAL_BEGIN) {
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, b);
      v = b;
      if (! useEL && v == z) break;
    }
    else {
      y = cm->cfirst[v] + yoffset;
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y);
      v = y;
      if (! useEL && v == z) break;
    }
  }

  free_vji_shadow_matrix(shadow, cm->M, j1, j0);
  return sc;
}

/* Function: v_splitter_hb()  [brief 26_0610-009] -- banded analogue of v_splitter(). */
static void
v_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
	      int r, int z, int i0, int i1, int j1, int j0, int useEL, CP9Bands_t *cp9b)
{
  float ***alpha, ***beta;
  float    sc;
  int      v, w, y;
  int      i, j, jp, ip, op;
  int      best_v, best_i, best_j;
  float    best_sc;
  int      midnode;
  int      b;
  float    bsc;

  /* 1. Boundary condition or small enough -> banded vinsideT_hb. */
  if (cm->ndidx[z] == cm->ndidx[r] + 1 || r == z ||
      vinsideT_size(cm, r, z, i0, i1, j1, j0) < RAMLIMIT) {
    vinsideT_hb(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, (r==0), cp9b);
    return;
  }

  /* 2. Find the split set, w..y. */
  midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r]) / 2);
  w = cm->nodemap[midnode];
  y = cm->cfirst[w]-1;

  /* 3. Banded vinside up to w, banded voutside down to y. */
  vinside_hb (cm, dsq, L, w, z, i0, i1, j1, j0, useEL, BE_EFFICIENT,
	      NULL, &alpha, NULL, (r==0), &b, &bsc, cp9b);
  voutside_hb(cm, dsq, L, r, y, i0, i1, j1, j0, useEL, BE_EFFICIENT,
	      NULL, &beta, cp9b);

  /* 4. Find the optimal split: v, i, j -- over in-band cells only. The loop
   * NESTING (v, then ip, then jp) and strict '>' MUST match the exact v_splitter:
   * with ties broken by first-to-reach-the-max, a different nesting would pick a
   * different equal-scoring split and diverge from the cm_CYKInsideAlignHB
   * parsetree (seen on LSU, where V-problems are big enough to actually split). */
  best_sc = IMPOSSIBLE;
  best_v  = -99;
  for (v = w; v <= y; v++)
    for (ip = 0; ip <= i1-i0; ip++) {
      i = ip + i0;
      for (jp = 0; jp <= j0-j1; jp++) {
	j = jp + j1;
	if (vji_inband(cp9b, v, j, i, i0,i1,j1,j0, &op) &&
	    (sc = alpha[v][jp][op] + beta[v][jp][op]) > best_sc) {
	  best_sc = sc;
	  best_v  = v;
	  best_i  = i;
	  best_j  = j;
	}
      }
    }

  /* Local ends: maybe better off in EL (full deck, no band). */
  if (useEL && (cm->flags & CMH_LOCAL_END)) {
    for (ip = 0; ip <= i1-i0; ip++)
      for (jp = 0; jp <= j0-j1; jp++)
	if ((sc = beta[cm->M][jp][ip]) > best_sc) {
	  best_sc = sc;
	  best_v  = -1;
	  best_i  = ip + i0;
	  best_j  = jp + j1;
	}
  }

  /* Local begins: maybe better off in root. */
  if (r == 0 && (cm->flags & CMH_LOCAL_BEGIN)) {
    if (bsc > best_sc) {
      best_sc = bsc;
      best_v  = -2;
      best_i  = i0;
      best_j  = j0;
    }
  }

  free_banded_hb_vji_matrix(alpha, cm, i0, i1, j1, j0, cp9b);
  free_banded_hb_vji_matrix(beta,  cm, i0, i1, j1, j0, cp9b);

  /* EXPLICIT empty-band path (THE brief-26_0610-003 hazard, designed out): if no in-band
   * split point beat IMPOSSIBLE, this V-problem is unreachable under the bands.
   * Do NOT recurse with degenerate best_i/best_j (that is the QDB negative-malloc
   * crash). Terminate cleanly via the banded base case (which defensively
   * attaches z). */
  if (best_v == -99 || ! NOT_IMPOSSIBLE(best_sc)) {
    vinsideT_hb(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, (r==0), cp9b);
    return;
  }

  if (best_v == -1) {
    v_splitter_hb(cm, dsq, L, tr, r, w, i0, best_i, best_j, j0, TRUE, cp9b);
    return;
  }
  if (best_v == -2) {
    if (b != z)
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b);
    v_splitter_hb(cm, dsq, L, tr, b, z, i0, i1, j1, j0, useEL, cp9b);
    return;
  }

  /* The optimal split into two V problems:
   *    V:  r..v,      i0..best_i, best_j..j0
   *    V:  best_v..z, best_i..i1, j1..best_j
   */
  v_splitter_hb(cm, dsq, L, tr, r,      best_v, i0,     best_i, best_j, j0, FALSE, cp9b);
  v_splitter_hb(cm, dsq, L, tr, best_v, z,      best_i, i1,     j1,     best_j, useEL, cp9b);
  return;
}

/*################################################################
 * Truncated banded HMM-banded D&C CYK (brief 26_0610-044, rung R4.4a: J-plane).
 *
 * These tr_*_hb functions are the TRUNCATED analogues of the non-truncated
 * banded D&C CYK family above (CYKDivideAndConquerHB + *_hb engines, briefs
 * 007-009). For rung R4.4a we fill/search the J (Joint) plane ONLY: the J-plane
 * interior recurrence is IDENTICAL to standard CYK (so the engine bodies are
 * near-verbatim copies of inside_hb/outside_hb/etc.), and the only genuinely-new
 * truncated machinery is the *truncated-begin ROOT*:
 *
 *   In truncated alignment the ONLY exit from ROOT_S (state 0) is a 'truncated
 *   begin' into ANY state v (with cp9b->Jvalid[v]) spanning the full subsequence,
 *   carrying a penalty cm->trp->{l,g}_ptyAA[pty_idx][v] (replaces cm->beginsc[v]).
 *   This REPLACES the normal ROOT_S->children descent entirely. We mirror the
 *   oracle cm_TrCYKInsideAlignHB()'s root block (cm_dpalign_trunc.c ~2787-2846):
 *   state 0 is excluded from the normal recurrence, and Jalpha[0][L][L] =
 *   max_v (Jalpha[v][L][L] + trpenalty[v]), recorded with a USED_TRUNC_BEGIN
 *   shadow. The begin search is penalty-AWARE (penalty folded into the score),
 *   so the D&C score == the oracle's Jalpha[0][L][L] byte-for-byte. Every
 *   parsetree node is tagged with its marginal mode (all TRMODE_J this rung) via
 *   InsertTraceNodewithMode(). The L/R/T planes are deferred to R4.4b/c; the
 *   r_allow_L/r_allow_R/v_allow_T flags are threaded through (FALSE this rung) so
 *   those rungs just turn on the planes.
 *
 * The truncated begin (penalty + Jvalid gating + USED_TRUNC_BEGIN + skip-state-0-
 * normal-recurrence) is the surgical change applied to inside/outside/vinside/
 * voutside; the splitters are copies that tag modes. The pty_idx and local/global
 * choice are held in file-static tr_dnc_pty_idx / tr_dnc_local, set once at the
 * TrCYKDivideAndConquerHB() entry (mirrors the cyk_dnc_track static pattern).
 *################################################################*/

static int tr_dnc_pty_idx = -1;   /* truncation-penalty index (TRPENALTY_*), set at entry */
static int tr_dnc_local   = 0;    /* TRUE if CMH_LOCAL_BEGIN (use l_ptyAA), else g_ptyAA  */

/* The (penalty-aware) truncated-begin score for entering state v, or IMPOSSIBLE. */
static float
tr_trpenalty(CM_t *cm, int v)
{
  if (tr_dnc_pty_idx < 0) return IMPOSSIBLE;
  return tr_dnc_local ? cm->trp->l_ptyAA[tr_dnc_pty_idx][v]
                      : cm->trp->g_ptyAA[tr_dnc_pty_idx][v];
}

/* TR_LR [brief 26_0610-045, R4.4b]: the L and R marginal planes for the truncated D&C
 * inside engine. The J plane stays in the existing 'alpha'/'shadow' (so the J
 * path is byte-identical to R4.4a); this bundle holds the additional L/R score
 * planes, their y/k shadows, the B-state child-mode shadows (Lkmode/Rkmode), and
 * the per-mode truncated-begin entry (Lb/Lbsc, Rb/Rbsc). fill_L/fill_R select
 * which planes to compute. When lr==NULL or neither fill flag is set, tr_inside_hb
 * behaves exactly as R4.4a (J only). The score planes are freed as the engine
 * walks (like J alpha); only the shadows are retained (for the traceback). */
typedef struct tr_lr_s {
  int      fill_L, fill_R;
  int      fill_T;               /* brief 26_0610-051 (R4.4c): also compute the T plane (B states only) */
  int      ret_planes;           /* if TRUE, return Lalpha/Ralpha instead of freeing */
  float ***Lalpha,  ***Ralpha;   /* L/R score planes (per-state banded vjd decks); also
				  * serve as the in/out array for splitter chaining   */
  void  ***Lshad,   ***Rshad;    /* L/R yshad(char**)/kshad(int**) shadows, by state */
  char  ***Lkmode,  ***Rkmode;   /* L/R B-state child-mode shadows (B states only)   */
  int      Lb, Rb;               /* L/R truncated-begin entry states                 */
  float    Lbsc, Rbsc;           /* their scores (penalty folded in)                 */
  /* T marginal (brief 26_0610-051): T mode appears ONLY at B states, reached ONLY by a
   * truncated begin into the B at full span -> children are R-left (BEGL) + L-right
   * (BEGR). So the T data is just the full-span T-combine per B state (Tfull) and
   * its k* split (Tfullk), plus the T begin (Tb/Tbsc). No T plane/outside is needed:
   * the begin scalar is the only T outside (oracle cm_TrCYKOutsideAlignHB:6743), and
   * the parent region above the B is degenerate (the begin jumps straight in). */
  float   *Tfull;                /* Talpha[v][j0][W] per state (IMPOSSIBLE off B/in-band) */
  int     *Tfullk;               /* the T-combine k* at the full-span cell, per B state   */
  int      Tb;                   /* T truncated-begin entry B state                       */
  float    Tbsc;                 /* its score (penalty folded in)                        */
} TR_LR;

/* TR_VLR [brief 26_0610-046]: the L/R marginal planes for the V-problem (vji) inside
 * engine, used by the whole-solve path for marginal/mixed-mode V-problems. The J
 * plane stays in the existing 'a'/'shadow'; this bundle holds the L/R vji score
 * planes, their shadows, and the B-... (no B in V-problems) child-mode shadows
 * Lmode/Rmode (the marginal->child-mode decode for the vji traceback). */
typedef struct tr_vlr_s {
  int      fill_L, fill_R;
  float ***La,    ***Ra;         /* L/R vji score planes (per-state banded vji decks) */
  char  ***Lsh,   ***Rsh;        /* L/R vji yshadows (yoffset, by state)              */
  char  ***Lmode, ***Rmode;      /* L/R marginal child-mode shadows                   */
  int      b_v, b_i, b_j, b_mode;/* best truncated-begin (root V-problem) entry        */
  float    b_sc;
  int      Lb, Rb;               /* L/R truncated-begin entry states (root V-problem)  */
  float    Lbsc, Rbsc;           /* their scores (penalty folded in)                  */
} TR_VLR;

/* free an L/R y/k shadow matrix bundle (void*** of char or int decks) allocated
 * by tr_inside_hb; mirrors free_vjd_shadow_matrix but for the banded L/R shadows.
 * NULL-safe per state and overall. */
static void
tr_free_lr_shadow(void ***shadow, char ***kmode, CM_t *cm, int i0, int j0)
{
  int v, j;
  if (shadow != NULL) {
    for (v = 0; v < cm->M; v++) {
      if (shadow[v] == NULL) continue;
      if (cm->sttype[v] == B_st) { for (j = i0-1; j <= j0; j++) if (((int **)shadow[v])[j]) free(((int **)shadow[v])[j]); }
      else                       { for (j = i0-1; j <= j0; j++) if (((char**)shadow[v])[j]) free(((char**)shadow[v])[j]); }
      free(shadow[v]);
    }
    free(shadow);
  }
  if (kmode != NULL) {
    for (v = 0; v < cm->M; v++) {
      if (kmode[v] == NULL) continue;
      for (j = i0-1; j <= j0; j++) if (kmode[v][j]) free(kmode[v][j]);
      free(kmode[v]);
    }
    free(kmode);
  }
}


/* forward declarations (mutually recursive) */
static float tr_generic_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
				    int r, int z, int i0, int j0,
				    int r_allow_J, int r_allow_L, int r_allow_R, int r_allow_T, CP9Bands_t *cp9b);
static float tr_wedge_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
				  int r, int z, int i0, int j0,
				  int r_allow_J, int r_allow_L, int r_allow_R, CP9Bands_t *cp9b);
static void  tr_v_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			      int r, int z, int i0, int i1, int j1, int j0, int useEL,
			      int r_allow_J, int r_allow_L, int r_allow_R,
			      int z_allow_J, int z_allow_L, int z_allow_R, CP9Bands_t *cp9b);
static float tr_inside_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
			  float ***alpha, float ****ret_alpha,
			  struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			  void ****ret_shadow, int allow_begin, int *ret_b, float *ret_bsc,
			  TR_LR *lr, CP9Bands_t *cp9b);
static void  tr_outside_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
			   int do_full, int fill_L, int fill_R,
			   float ***beta, float ****ret_beta,
			   float ***betaL, float ****ret_betaL,
			   float ***betaR, float ****ret_betaR,
			   struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
			   float *ret_bsc, int *ret_bv, int *ret_bj, int *ret_bmode, CP9Bands_t *cp9b);
static float tr_insideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			   int r, int z, int i0, int j0, int allow_begin,
			   int r_allow_J, int r_allow_L, int r_allow_R, int r_allow_T, CP9Bands_t *cp9b);
static float tr_vinside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
			   int r, int z, int i0, int i1, int j1, int j0, int useEL,
			   int do_full, float ***a, float ****ret_a, char ****ret_shadow,
			   int allow_begin, int *ret_b, float *ret_bsc,
			   int z_allow_J, int z_allow_L, int z_allow_R,
			   TR_VLR *vlr, CP9Bands_t *cp9b);
static void  tr_voutside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
			    int r, int z, int i0, int i1, int j1, int j0, int useEL,
			    int do_full, float ***beta, float ****ret_beta, CP9Bands_t *cp9b);
static float tr_vinsideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
			    int r, int z, int i0, int i1, int j1, int j0, int useEL,
			    int allow_begin, int r_allow_J, int r_allow_L, int r_allow_R,
			    int z_allow_J, int z_allow_L, int z_allow_R, CP9Bands_t *cp9b);

/* Function: tr_inside_hb()  [brief 26_0610-044, R4.4a]
 *
 * Purpose:  J-plane truncated analogue of inside_hb(). Identical banded CYK
 *           recurrence, EXCEPT: (1) state 0 (ROOT_S) gets NO normal transitions
 *           when allow_begin (its alpha is set only by the truncated begin), and
 *           (2) the begin bookkeeping uses the penalty-aware truncated begin
 *           (tr_trpenalty(), over all cp9b->Jvalid[] states, marked
 *           USED_TRUNC_BEGIN) rather than cm->beginsc[]/USED_LOCAL_BEGIN.
 */
static float
tr_inside_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	     float ***alpha, float ****ret_alpha,
	     struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	     void ****ret_shadow, int allow_begin, int *ret_b, float *ret_bsc,
	     TR_LR *lr, CP9Bands_t *cp9b)
{
  int      status;
  int     *touch;
  int      v,y,z;
  int      j,d,i,k;
  float    sc;
  int      yoffset;
  int      W;
  int      jp;
  int      jp_v;
  void  ***shadow;
  int    **kshad;
  char   **yshad;
  int      b;
  float    bsc;
  int      jn, jx;
  int      yy, zz;
  int      jp_y, jp_z;
  int      kn, kx;
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  /* L/R marginal-plane state (brief 26_0610-045, R4.4b). All NULL/off when lr==NULL. */
  int      fill_L = (lr != NULL) ? lr->fill_L : FALSE;
  int      fill_R = (lr != NULL) ? lr->fill_R : FALSE;
  int      fill_T = (lr != NULL) ? lr->fill_T : FALSE;  /* brief 26_0610-051 (R4.4c) */
  int      ret_planes = (lr != NULL) ? lr->ret_planes : FALSE;
  /* T (brief 26_0610-051) needs the L+R child planes to form the B-state T-combine
   * (Ralpha[BEGL] + Lalpha[BEGR]); force them on whenever fill_T. */
  float   *Tfull  = NULL;     /* Talpha[v][j0][W] per state (begin scan + traceback) */
  int     *Tfullk = NULL;     /* the T-combine k* at the full-span cell, per B state  */
  int      Tb   = -1;
  float    Tbsc = IMPOSSIBLE;
  float ***Lalpha = NULL, ***Ralpha = NULL;
  void  ***Lshad  = NULL, ***Rshad  = NULL;  /* yshad(char**)/kshad(int**) per state */
  char  ***Lkmode = NULL, ***Rkmode = NULL;  /* B-state child-mode shadow            */
  int      Lb = -1, Rb = -1;
  float    Lbsc = IMPOSSIBLE, Rbsc = IMPOSSIBLE;
  int      sd, sdl, sdr;
  int    **Lkshad = NULL, **Rkshad = NULL;
  char   **Lyshad = NULL, **Ryshad = NULL;

  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;
  if (fill_T) { fill_L = TRUE; fill_R = TRUE; }   /* T-combine reads R(BEGL)+L(BEGR) */
  if (dpool == NULL) dpool = deckpool_create();

  if (fill_T) {
    ESL_ALLOC(Tfull,  sizeof(float) * cm->M);
    ESL_ALLOC(Tfullk, sizeof(int)   * cm->M);
    for (v = 0; v < cm->M; v++) { Tfull[v] = IMPOSSIBLE; Tfullk[v] = 0; }
  }

  if (alpha == NULL) {
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, sizeof(int) * (cm->M+1));
  for (v = 0;     v < vroot; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) touch[v] = cm->pnum[v];
  for (v = vend+1;v < cm->M; v++) touch[v] = 0;

  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(void **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
  }

  /* L/R score planes: reuse the caller's array if chaining (lr->Lalpha != NULL,
   * the generic-splitter w-then-y pattern), else allocate fresh. Freed as we walk
   * (like J) unless ret_planes, in which case the kept decks are returned. */
  if (fill_L) {
    Lalpha = (lr != NULL) ? lr->Lalpha : NULL;
    if (Lalpha == NULL) { ESL_ALLOC(Lalpha, sizeof(float **) * (cm->M+1)); for (v = 0; v <= cm->M; v++) Lalpha[v] = NULL; }
  }
  if (fill_R) {
    Ralpha = (lr != NULL) ? lr->Ralpha : NULL;
    if (Ralpha == NULL) { ESL_ALLOC(Ralpha, sizeof(float **) * (cm->M+1)); for (v = 0; v <= cm->M; v++) Ralpha[v] = NULL; }
  }
  if (ret_shadow != NULL && fill_L) {
    ESL_ALLOC(Lshad,  sizeof(void **) * cm->M); for (v = 0; v < cm->M; v++) Lshad[v]  = NULL;
    ESL_ALLOC(Lkmode, sizeof(char **) * cm->M); for (v = 0; v < cm->M; v++) Lkmode[v] = NULL;
  }
  if (ret_shadow != NULL && fill_R) {
    ESL_ALLOC(Rshad,  sizeof(void **) * cm->M); for (v = 0; v < cm->M; v++) Rshad[v]  = NULL;
    ESL_ALLOC(Rkmode, sizeof(char **) * cm->M); for (v = 0; v < cm->M; v++) Rkmode[v] = NULL;
  }

  for (v = vend; v >= vroot; v--)
    {
      sd  = StateDelta(cm->sttype[v]);
      sdl = StateLeftDelta(cm->sttype[v]);
      sdr = StateRightDelta(cm->sttype[v]);
      alpha[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
      if (fill_L) Lalpha[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
      if (fill_R) Ralpha[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);

      if (ret_shadow != NULL && cm->sttype[v] != E_st) {
	if (cm->sttype[v] == B_st) {
	  kshad     = alloc_banded_hb_vjd_kshadow_deck(L, i0, j0, v, cp9b);
	  shadow[v] = (void **) kshad;
	  if (fill_L) { Lkshad = alloc_banded_hb_vjd_kshadow_deck(L, i0, j0, v, cp9b); Lshad[v] = (void **) Lkshad;
	                Lkmode[v] = alloc_banded_hb_vjd_yshadow_deck(L, i0, j0, v, cp9b); }
	  if (fill_R) { Rkshad = alloc_banded_hb_vjd_kshadow_deck(L, i0, j0, v, cp9b); Rshad[v] = (void **) Rkshad;
	                Rkmode[v] = alloc_banded_hb_vjd_yshadow_deck(L, i0, j0, v, cp9b); }
	} else {
	  yshad     = alloc_banded_hb_vjd_yshadow_deck(L, i0, j0, v, cp9b);
	  shadow[v] = (void **) yshad;
	  if (fill_L) { Lyshad = alloc_banded_hb_vjd_yshadow_deck(L, i0, j0, v, cp9b); Lshad[v] = (void **) Lyshad; }
	  if (fill_R) { Ryshad = alloc_banded_hb_vjd_yshadow_deck(L, i0, j0, v, cp9b); Rshad[v] = (void **) Ryshad; }
	}
      } else { Lkshad = Rkshad = NULL; Lyshad = Ryshad = NULL; }

      /* per-state per-mode validity gating (mirrors the HB oracle's do_{L,R}_v).
       * Jvalid[v] is always TRUE for in-band states (so the J path needs no gate;
       * that is why R4.4a matched byte-exact), but Lvalid/Rvalid are selective. */
      int do_L_v = fill_L && cp9b->Lvalid[v];
      int do_R_v = fill_R && cp9b->Rvalid[v];

      jn = ESL_MAX(i0-1, jmin[v]);
      jx = ESL_MIN(j0,   jmax[v]);

      for (j = jn; j <= jx; j++) {
	jp_v = j - jmin[v];
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) {
	  int dpi = d - hd_min(cp9b, v, jp_v);
	  alpha[v][j][dpi] = IMPOSSIBLE;
	  if (fill_L) Lalpha[v][j][dpi] = IMPOSSIBLE;
	  if (fill_R) Ralpha[v][j][dpi] = IMPOSSIBLE;
	  /* shadow init: B states use k-shadows (init 0) + kmode (init J); other
	   * states use y-shadows (init USED_EL, mirroring the oracle). */
	  if (cm->sttype[v] == B_st) {
	    if (fill_L && Lkshad) { Lkshad[j][dpi] = 0; Lkmode[v][j][dpi] = TRMODE_J; }
	    if (fill_R && Rkshad) { Rkshad[j][dpi] = 0; Rkmode[v][j][dpi] = TRMODE_J; }
	  } else {
	    if (fill_L && Lyshad) Lyshad[j][dpi] = USED_EL;
	    if (fill_R && Ryshad) Ryshad[j][dpi] = USED_EL;
	  }
	}
      }

      if (allow_begin && v == 0)
	{
	  /* TRUNCATED: ROOT_S has NO normal transitions; alpha[0] is set only by
	   * the truncated begin (handled in the begin bookkeeping below). */
	}
      else if (cm->sttype[v] == E_st)
	{
	  int dpe;
	  for (j = jn; j <= jx; j++)
	    if (hb_inband(cp9b, v, j, 0, i0, j0, &dpe)) {
	      alpha[v][j][dpe] = 0.;
	      if (do_L_v) Lalpha[v][j][dpe] = 0.;
	      if (do_R_v) Ralpha[v][j][dpe] = 0.;
	    }
	}
      else if (cm->sttype[v] == D_st || cm->sttype[v] == S_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo;
		    if (hb_inband(cp9b, y+yoffset, j, d, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* L/R marginal (brief 26_0610-045): D/S do not emit, so the marginal mode
		 * passes through unchanged (L child stays L, R child stays R). */
		if (do_L_v || do_R_v) {
		  if (d == 0) {
		    if (do_L_v) { Lalpha[v][j][dp_v] = IMPOSSIBLE; if (ret_shadow != NULL && cm->sttype[v] == S_st) Lyshad[j][dp_v] = USED_TRUNC_END; }
		    if (do_R_v) { Ralpha[v][j][dp_v] = IMPOSSIBLE; if (ret_shadow != NULL && cm->sttype[v] == S_st) Ryshad[j][dp_v] = USED_TRUNC_END; }
		  } else {
		    int dp_yo;
		    /* brief 26_0610-096: marginal EL for a D/S state (d>0). Mirrors the J-plane
		     * EL init at the top of this D/S branch (6653) and the emitting-state
		     * marginal EL below, and the oracle's state-type-INDEPENDENT EL reinit
		     * (cm_TrCYKInsideAlignHB ~7011-7069). Without it, the root-split candidate
		     * scoring (which reads Ralpha[BEGL]) cannot credit a BEGL_S->EL R-mode
		     * terminus, so on CsrB-sample13_5ptr it picked a suboptimal split (the
		     * wedge-level OUTSIDE betaR[cm->M] found the EL, but the INSIDE used for
		     * split selection did not -> a 0.73-bit inside/outside inconsistency).
		     * sdl==sdr==0 for D/S; endsc[v]==IMPOSSIBLE for non-EL states leaves the
		     * cell effectively IMPOSSIBLE (clamped below), so this is a no-op there.
		     * Gated on Lvalid[cm->M]/Rvalid[cm->M] EXACTLY as the oracle (7031/7048):
		     * the marginal EL *deck* must be valid for this mode, else no marginal EL
		     * (omitting this gate over-credited a phantom L-EL on ar45-sample26_p60,
		     * flipping its resolved mode L->R for a 0.6-bit-worse parse). */
		    if (do_L_v && cp9b->Lvalid[cm->M]) { Lalpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdl)); if (ret_shadow != NULL) Lyshad[j][dp_v] = USED_EL; }
		    if (do_R_v && cp9b->Rvalid[cm->M]) { Ralpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdr)); if (ret_shadow != NULL) Ryshad[j][dp_v] = USED_EL; }
		    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		      int yy2 = cm->cfirst[v] + yoffset;
		      if (do_L_v && cp9b->Lvalid[yy2] && hb_inband(cp9b, yy2, j, d, i0, j0, &dp_yo) &&
			  (sc = Lalpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Lalpha[v][j][dp_v]) {
			Lalpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Lyshad[j][dp_v] = yoffset + TRMODE_L_OFFSET;
		      }
		      if (do_R_v && cp9b->Rvalid[yy2] && hb_inband(cp9b, yy2, j, d, i0, j0, &dp_yo) &&
			  (sc = Ralpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Ralpha[v][j][dp_v]) {
			Ralpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Ryshad[j][dp_v] = yoffset + TRMODE_R_OFFSET;
		      }
		    }
		    if (do_L_v && Lalpha[v][j][dp_v] < IMPOSSIBLE) Lalpha[v][j][dp_v] = IMPOSSIBLE;
		    if (do_R_v && Ralpha[v][j][dp_v] < IMPOSSIBLE) Ralpha[v][j][dp_v] = IMPOSSIBLE;
		  }
		}
	      }
	  }
	}
      else if (cm->sttype[v] == B_st)
	{
	  int do_L_y, do_L_z, do_R_y, do_R_z;
	  yy = cm->cfirst[v];
	  zz = cm->cnum[v];
	  do_L_y = fill_L && cp9b->Lvalid[yy]; do_L_z = fill_L && cp9b->Lvalid[zz];
	  do_R_y = fill_R && cp9b->Rvalid[yy]; do_R_z = fill_R && cp9b->Rvalid[zz];
	  jn = ESL_MAX(jn, jmin[zz]);
	  jx = ESL_MIN(jx, jmax[zz]);
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    jp_y = j - jmin[yy];
	    jp_z = j - jmin[zz];
	    kn = ESL_MAX(j - jmax[yy], hd_min(cp9b, zz, jp_z));
	    kn = ESL_MAX(kn, 0);
	    kx = ESL_MIN(jp_y, hd_max(cp9b, zz, jp_z));
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		for (k = kn; k <= kx; k++)
		  if ((k >= d - hd_max(cp9b, yy, jp_y-k)) && (k <= d - hd_min(cp9b, yy, jp_y-k)))
		    {
		      int dp_yk = (d-k) - hd_min(cp9b, yy, jp_y-k);
		      int dp_zk =  k    - hd_min(cp9b, zz, jp_z);
		      if ((sc = alpha[yy][j-k][dp_yk] + alpha[zz][j][dp_zk]) > alpha[v][j][dp_v]) {
			alpha[v][j][dp_v] = sc;
			if (ret_shadow != NULL) kshad[j][dp_v] = k;
		      }
		      /* L: left child full (J), right child L marginal.  yy's J-plane
		       * cross-term is only a valid truncated contribution if yy's own
		       * node lies entirely within the observed (non-truncated) sequence
		       * (brief 26_0610-070: cp9b->Jvalid[yy], mirrors the do_J_y convention used
		       * throughout cm_dpalign_trunc.c's stock truncated recursions). */
		      if (do_L_v && do_L_z && cp9b->Jvalid[yy] &&
			  (sc = alpha[yy][j-k][dp_yk] + Lalpha[zz][j][dp_zk]) > Lalpha[v][j][dp_v]) {
			Lalpha[v][j][dp_v] = sc;
			if (ret_shadow != NULL) { Lkshad[j][dp_v] = k; Lkmode[v][j][dp_v] = TRMODE_J; }
		      }
		      /* R: left child R marginal, right child full (J).  Same Jvalid[zz]
		       * gate on the right child's J-plane cross-term (brief 26_0610-070). */
		      if (do_R_v && do_R_y && cp9b->Jvalid[zz] &&
			  (sc = Ralpha[yy][j-k][dp_yk] + alpha[zz][j][dp_zk]) > Ralpha[v][j][dp_v]) {
			Ralpha[v][j][dp_v] = sc;
			if (ret_shadow != NULL) { Rkshad[j][dp_v] = k; Rkmode[v][j][dp_v] = TRMODE_J; }
		      }
		      /* T (brief 26_0610-051): both ends truncated -> left child R, right child L;
		       * k != 0, k != d (both children non-empty). Only the full-span cell
		       * (j0,W) is needed -- that is the only place T occurs (the begin enters
		       * the B at full span). Transcribed from oracle cm_TrCYKInsideAlignHB
		       * (cm_dpalign_trunc.c:2692-2698). */
		      if (fill_T && j == j0 && d == W && k != 0 && k != d &&
			  do_R_y && do_L_z &&
			  (sc = Ralpha[yy][j-k][dp_yk] + Lalpha[zz][j][dp_zk]) > Tfull[v]) {
			Tfull[v] = sc; Tfullk[v] = k;
		      }
		    }
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	  /* B special case 1: k==0, full sequence on left child (right child empty). */
	  if (do_L_v) {
	    int jnn = ESL_MAX(jmin[v], jmin[yy]);
	    int jxx = ESL_MIN(jmax[v], jmax[yy]);
	    for (j = ESL_MAX(i0-1, jnn); j <= ESL_MIN(j0, jxx); j++) {
	      int jpp = j - (i0-1);
	      int dnn = ESL_MAX(hd_min(cp9b, v, j-jmin[v]), hd_min(cp9b, yy, j-jmin[yy]));
	      int dxx = ESL_MIN(hd_max(cp9b, v, j-jmin[v]), hd_max(cp9b, yy, j-jmin[yy]));
	      for (d = dnn; d <= dxx && d <= jpp; d++) {
		int dp_v = d - hd_min(cp9b, v, j-jmin[v]);
		int dp_y = d - hd_min(cp9b, yy, j-jmin[yy]);
		/* brief 26_0610-070: yy's J-plane is only a valid contribution here if yy's own
		 * node lies entirely within the observed sequence (cp9b->Jvalid[yy]). */
		if (cp9b->Jvalid[yy] && (sc = alpha[yy][j][dp_y]) > Lalpha[v][j][dp_v]) {
		  Lalpha[v][j][dp_v] = sc;
		  if (ret_shadow != NULL) { Lkshad[j][dp_v] = 0; Lkmode[v][j][dp_v] = TRMODE_J; }
		}
		if (do_L_y && (sc = Lalpha[yy][j][dp_y]) > Lalpha[v][j][dp_v]) {
		  Lalpha[v][j][dp_v] = sc;
		  if (ret_shadow != NULL) { Lkshad[j][dp_v] = 0; Lkmode[v][j][dp_v] = TRMODE_L; }
		}
	      }
	    }
	  }
	  /* B special case 2: k==d, full sequence on right child (left child empty). */
	  if (do_R_v) {
	    int jnn = ESL_MAX(jmin[v], jmin[zz]);
	    int jxx = ESL_MIN(jmax[v], jmax[zz]);
	    for (j = ESL_MAX(i0-1, jnn); j <= ESL_MIN(j0, jxx); j++) {
	      int jpp = j - (i0-1);
	      int dnn = ESL_MAX(hd_min(cp9b, v, j-jmin[v]), hd_min(cp9b, zz, j-jmin[zz]));
	      int dxx = ESL_MIN(hd_max(cp9b, v, j-jmin[v]), hd_max(cp9b, zz, j-jmin[zz]));
	      for (d = dnn; d <= dxx && d <= jpp; d++) {
		int dp_v = d - hd_min(cp9b, v, j-jmin[v]);
		int dp_z = d - hd_min(cp9b, zz, j-jmin[zz]);
		/* brief 26_0610-070: zz's J-plane is only a valid contribution here if zz's own
		 * node lies entirely within the observed sequence (cp9b->Jvalid[zz]).
		 * Without this gate, a RIGHT_FULL bifurcation can pull in zz's plain
		 * (truncation-unaware) classical Inside value even when zz's node
		 * structurally spans past the truncation boundary -- exactly the
		 * crash reproduced in brief 26_0610-070 (frag5p_s12 L24, v=518/z=566). */
		if (cp9b->Jvalid[zz] && (sc = alpha[zz][j][dp_z]) > Ralpha[v][j][dp_v]) {
		  Ralpha[v][j][dp_v] = sc;
		  if (ret_shadow != NULL) { Rkshad[j][dp_v] = d; Rkmode[v][j][dp_v] = TRMODE_J; }
		}
		if (do_R_z && (sc = Ralpha[zz][j][dp_z]) > Ralpha[v][j][dp_v]) {
		  Ralpha[v][j][dp_v] = sc;
		  if (ret_shadow != NULL) { Rkshad[j][dp_v] = d; Rkmode[v][j][dp_v] = TRMODE_R; }
		}
	      }
	    }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo;
		    if (hb_inband(cp9b, y+yoffset, j-1, d-2, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j-1][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* L marginal: MP emits only the LEFT residue; transit child uses (y,j,d-sdl). */
		if (do_L_v) {
		  int dp_yo;
		  Lalpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdl));
		  if (ret_shadow != NULL) Lyshad[j][dp_v] = USED_EL;
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		    int yy2 = cm->cfirst[v] + yoffset;
		    if (hb_inband(cp9b, yy2, j, d-sdl, i0, j0, &dp_yo)) {
		      /* brief 26_0610-093: gate the MP L<-J cross-term on cp9b->Jvalid[yy2]
		       * (oracle cm_TrCYKInsideAlignHB MP L-recursion, do_J_y). Same phantom-J
		       * over-score as the ML/IL and IR/MR cross fixes -- fires for structured
		       * (MATP) subtrees in L mode. */
		      if (cp9b->Jvalid[yy2] && (sc = alpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Lalpha[v][j][dp_v]) {
			Lalpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Lyshad[j][dp_v] = yoffset + TRMODE_J_OFFSET;
		      }
		      if (cp9b->Lvalid[yy2] && (sc = Lalpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Lalpha[v][j][dp_v]) {
			Lalpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Lyshad[j][dp_v] = yoffset + TRMODE_L_OFFSET;
		      }
		    }
		  }
		  if (d >= 2) Lalpha[v][j][dp_v] += cm->lmesc[v][dsq[i]];
		  else { Lalpha[v][j][dp_v] = cm->lmesc[v][dsq[i]]; if (ret_shadow != NULL) Lyshad[j][dp_v] = USED_TRUNC_END; }
		  if (Lalpha[v][j][dp_v] < IMPOSSIBLE) Lalpha[v][j][dp_v] = IMPOSSIBLE;
		}
		/* R marginal: MP emits only the RIGHT residue; transit child uses (y,j-sdr,d-sdr). */
		if (do_R_v) {
		  int dp_yo;
		  Ralpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdr));
		  if (ret_shadow != NULL) Ryshad[j][dp_v] = USED_EL;
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		    int yy2 = cm->cfirst[v] + yoffset;
		    if (hb_inband(cp9b, yy2, j-sdr, d-sdr, i0, j0, &dp_yo)) {
		      /* brief 26_0610-093: gate the MP R<-J cross-term on cp9b->Jvalid[yy2]
		       * (oracle MP R-recursion, do_J_y). R-mirror of the MP L fix above. */
		      if (cp9b->Jvalid[yy2] && (sc = alpha[yy2][j-sdr][dp_yo] + cm->tsc[v][yoffset]) > Ralpha[v][j][dp_v]) {
			Ralpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Ryshad[j][dp_v] = yoffset + TRMODE_J_OFFSET;
		      }
		      if (cp9b->Rvalid[yy2] && (sc = Ralpha[yy2][j-sdr][dp_yo] + cm->tsc[v][yoffset]) > Ralpha[v][j][dp_v]) {
			Ralpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Ryshad[j][dp_v] = yoffset + TRMODE_R_OFFSET;
		      }
		    }
		  }
		  if (d >= 2) Ralpha[v][j][dp_v] += cm->rmesc[v][dsq[j]];
		  else { Ralpha[v][j][dp_v] = cm->rmesc[v][dsq[j]]; if (ret_shadow != NULL) Ryshad[j][dp_v] = USED_TRUNC_END; }
		  if (Ralpha[v][j][dp_v] < IMPOSSIBLE) Ralpha[v][j][dp_v] = IMPOSSIBLE;
		}
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo;
		    if (hb_inband(cp9b, y+yoffset, j, d-1, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][dsq[i]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* L marginal: IL/ML emit left; child uses (y,j,d-sd). brief 26_0610-092:
		 * an L-mode LEFT-emitter transits ONLY to its child's L-plane, NEVER to the
		 * child's J-plane -- matching the trusted oracle cm_TrCYKInsideAlignHB ML/IL
		 * L-recursion (cm_dpalign_trunc.c:7050-7054, which gates on do_L_y only and
		 * has no Jalpha[y] term). The prior extra alpha[yy2] (J-child) term was an
		 * inside-begin inflation: where a child is not Lvalid (or Jalpha[y]>Lalpha[y]),
		 * it fabricated an L-parse via a joint child, over-scoring the L-begin scan
		 * (contrast the MP L-plane and the ML/IL R-plane below, where the oracle DOES
		 * read the J child -- that asymmetry is exactly the marginal-mode rule). */
		if (do_L_v) {
		  int dp_yo;
		  Lalpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdl));
		  if (ret_shadow != NULL) Lyshad[j][dp_v] = USED_EL;
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		    int yy2 = cm->cfirst[v] + yoffset;
		    if (cp9b->Lvalid[yy2] && hb_inband(cp9b, yy2, j, d-sd, i0, j0, &dp_yo)) {
		      if ((sc = Lalpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Lalpha[v][j][dp_v]) {
			Lalpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Lyshad[j][dp_v] = yoffset + TRMODE_L_OFFSET;
		      }
		    }
		  }
		  if (d >= 2) { if (dsq[i] < cm->abc->K) Lalpha[v][j][dp_v] += cm->esc[v][dsq[i]];
		               else Lalpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]); }
		  else { if (dsq[i] < cm->abc->K) Lalpha[v][j][dp_v] = cm->esc[v][dsq[i]];
		         else Lalpha[v][j][dp_v] = esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		         if (ret_shadow != NULL) Lyshad[j][dp_v] = USED_TRUNC_END; }
		  if (Lalpha[v][j][dp_v] < IMPOSSIBLE) Lalpha[v][j][dp_v] = IMPOSSIBLE;
		}
		/* R marginal: IL/ML emit nothing in R mode; child uses (y,j,d); no IL self-transit. */
		if (do_R_v) {
		  int dp_yo, Ryoffset0 = (cm->sttype[v] == IL_st) ? 1 : 0;
		  Ralpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdr));
		  if (ret_shadow != NULL) Ryshad[j][dp_v] = USED_EL;
		  for (yoffset = Ryoffset0; yoffset < cm->cnum[v]; yoffset++) {
		    int yy2 = cm->cfirst[v] + yoffset;
		    if (hb_inband(cp9b, yy2, j, d, i0, j0, &dp_yo)) {
		      /* brief 26_0610-093: gate the R<-J cross-term (R-mode left-emitter
		       * converts to a joint child) on cp9b->Jvalid[yy2], matching the trusted
		       * oracle cm_TrCYKInsideAlignHB ML/IL R-recursion (cm_dpalign_trunc.c:7091,
		       * do_J_y). Without it the D&C's always-filled J deck alpha[yy2] for a
		       * !Jvalid child fabricates a phantom R->J conversion, over-scoring Ralpha[v]
		       * at off-optimal cells -> the structured-case generic-splitter score
		       * inflation (parse still correct; the traceback ignores the phantom). */
		      if (cp9b->Jvalid[yy2] && (sc = alpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Ralpha[v][j][dp_v]) {
			Ralpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Ryshad[j][dp_v] = yoffset + TRMODE_J_OFFSET;
		      }
		      if (cp9b->Rvalid[yy2] && (sc = Ralpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Ralpha[v][j][dp_v]) {
			Ralpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Ryshad[j][dp_v] = yoffset + TRMODE_R_OFFSET;
		      }
		    }
		  }
		  if (Ralpha[v][j][dp_v] < IMPOSSIBLE) Ralpha[v][j][dp_v] = IMPOSSIBLE;
		}
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = jn; j <= jx; j++) {
	    jp   = j - (i0-1);
	    jp_v = j - jmin[v];
	    for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	      {
		int dp_v = d - hd_min(cp9b, v, jp_v);
		y = cm->cfirst[v];
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		  {
		    int dp_yo;
		    if (hb_inband(cp9b, y+yoffset, j-1, d-1, i0, j0, &dp_yo) &&
			(sc = alpha[y+yoffset][j-1][dp_yo] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v]) {
		      alpha[v][j][dp_v] = sc;
		      if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
		    }
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][dsq[j]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* R marginal: IR/MR emit right; child uses (y,j-sdr,d-sd). brief 26_0610-092:
		 * symmetric to the ML/IL L-plane fix above -- an R-mode RIGHT-emitter transits
		 * ONLY to its child's R-plane, NEVER the J-plane (oracle cm_TrCYKInsideAlignHB
		 * IR/MR R-recursion, cm_dpalign_trunc.c:7096-7100, gates do_R_y only). The prior
		 * extra alpha[yy2] (J-child) term was the R-mode mirror of the inside-begin
		 * inflation. */
		if (do_R_v) {
		  int dp_yo;
		  Ralpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdr));
		  if (ret_shadow != NULL) Ryshad[j][dp_v] = USED_EL;
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		    int yy2 = cm->cfirst[v] + yoffset;
		    if (cp9b->Rvalid[yy2] && hb_inband(cp9b, yy2, j-sdr, d-sd, i0, j0, &dp_yo)) {
		      if ((sc = Ralpha[yy2][j-sdr][dp_yo] + cm->tsc[v][yoffset]) > Ralpha[v][j][dp_v]) {
			Ralpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Ryshad[j][dp_v] = yoffset + TRMODE_R_OFFSET;
		      }
		    }
		  }
		  if (d >= 2) { if (dsq[j] < cm->abc->K) Ralpha[v][j][dp_v] += cm->esc[v][dsq[j]];
		               else Ralpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]); }
		  else { if (dsq[j] < cm->abc->K) Ralpha[v][j][dp_v] = cm->esc[v][dsq[j]];
		         else Ralpha[v][j][dp_v] = esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		         if (ret_shadow != NULL) Ryshad[j][dp_v] = USED_TRUNC_END; }
		  if (Ralpha[v][j][dp_v] < IMPOSSIBLE) Ralpha[v][j][dp_v] = IMPOSSIBLE;
		}
		/* L marginal: IR/MR emit nothing in L mode; child uses (y,j,d); no IR self-transit. */
		if (do_L_v) {
		  int dp_yo, Lyoffset0 = (cm->sttype[v] == IR_st) ? 1 : 0;
		  Lalpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d - sdl));
		  if (ret_shadow != NULL) Lyshad[j][dp_v] = USED_EL;
		  for (yoffset = Lyoffset0; yoffset < cm->cnum[v]; yoffset++) {
		    int yy2 = cm->cfirst[v] + yoffset;
		    if (hb_inband(cp9b, yy2, j, d, i0, j0, &dp_yo)) {
		      /* brief 26_0610-093: gate the L<-J cross-term on cp9b->Jvalid[yy2] (the
		       * R-mirror of the ML/IL R-plane fix above; oracle IR/MR L-recursion
		       * cm_dpalign_trunc.c:7216, do_J_y). Same phantom-J-child inflation. */
		      if (cp9b->Jvalid[yy2] && (sc = alpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Lalpha[v][j][dp_v]) {
			Lalpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Lyshad[j][dp_v] = yoffset + TRMODE_J_OFFSET;
		      }
		      if (cp9b->Lvalid[yy2] && (sc = Lalpha[yy2][j][dp_yo] + cm->tsc[v][yoffset]) > Lalpha[v][j][dp_v]) {
			Lalpha[v][j][dp_v] = sc; if (ret_shadow != NULL) Lyshad[j][dp_v] = yoffset + TRMODE_L_OFFSET;
		      }
		    }
		  }
		  if (Lalpha[v][j][dp_v] < IMPOSSIBLE) Lalpha[v][j][dp_v] = IMPOSSIBLE;
		}
	      }
	  }
	}

      /* TRUNCATED-BEGIN bookkeeping (penalty-aware; replaces non-trunc local begin).
       * Per mode: gated on {J,L,R}valid[v]; penalty folded into the score. */
      if (allow_begin && v != 0 && cp9b->Jvalid[v]) {
	int dpb;
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) &&
	    hb_inband(cp9b, v, j0, W, i0, j0, &dpb) &&
	    alpha[v][j0][dpb] + trpen > bsc)
	  {
	    b   = v;
	    bsc = alpha[v][j0][dpb] + trpen;
	  }
      }
      if (allow_begin && v != 0 && do_L_v && cp9b->Lvalid[v]) {
	int dpb;
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) &&
	    hb_inband(cp9b, v, j0, W, i0, j0, &dpb) &&
	    Lalpha[v][j0][dpb] + trpen > Lbsc)
	  { Lb = v; Lbsc = Lalpha[v][j0][dpb] + trpen; }
      }
      if (allow_begin && v != 0 && do_R_v && cp9b->Rvalid[v]) {
	int dpb;
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) &&
	    hb_inband(cp9b, v, j0, W, i0, j0, &dpb) &&
	    Ralpha[v][j0][dpb] + trpen > Rbsc)
	  { Rb = v; Rbsc = Ralpha[v][j0][dpb] + trpen; }
      }
      /* T begin (brief 26_0610-051): only B states, gated cp9b->Tvalid[0] && Tvalid[v]; the
       * begin enters the B in T mode at full span (oracle cm_TrCYKInsideAlignHB:2837). */
      if (allow_begin && v != 0 && fill_T && cm->sttype[v] == B_st &&
	  cp9b->Tvalid[0] && cp9b->Tvalid[v] && NOT_IMPOSSIBLE(Tfull[v])) {
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) && Tfull[v] + trpen > Tbsc)
	  { Tb = v; Tbsc = Tfull[v] + trpen; }
      }
      if (allow_begin && v == 0 && cp9b->Jvalid[0]) {
	int dpb0;
	if (hb_inband(cp9b, 0, j0, W, i0, j0, &dpb0)) {
	  alpha[0][j0][dpb0] = bsc;
	  if (ret_shadow != NULL) yshad[j0][dpb0] = USED_TRUNC_BEGIN;
	}
      }
      if (allow_begin && v == 0 && fill_L && cp9b->Lvalid[0]) {
	int dpb0;
	if (hb_inband(cp9b, 0, j0, W, i0, j0, &dpb0)) {
	  Lalpha[0][j0][dpb0] = Lbsc;
	  if (ret_shadow != NULL) Lyshad[j0][dpb0] = USED_TRUNC_BEGIN;
	}
      }
      if (allow_begin && v == 0 && fill_R && cp9b->Rvalid[0]) {
	int dpb0;
	if (hb_inband(cp9b, 0, j0, W, i0, j0, &dpb0)) {
	  Ralpha[0][j0][dpb0] = Rbsc;
	  if (ret_shadow != NULL) Ryshad[j0][dpb0] = USED_TRUNC_BEGIN;
	}
      }

      if (! do_full) {
	if (cm->sttype[v] == B_st)
	  {
	    y = cm->cfirst[v]; free_banded_hb_vjd_deck(alpha[y], i0, j0, y, cp9b); alpha[y] = NULL;
	    z = cm->cnum[v];   free_banded_hb_vjd_deck(alpha[z], i0, j0, z, cp9b); alpha[z] = NULL;
	    if (fill_L) { free_banded_hb_vjd_deck(Lalpha[y], i0, j0, y, cp9b); Lalpha[y] = NULL;
	                  free_banded_hb_vjd_deck(Lalpha[z], i0, j0, z, cp9b); Lalpha[z] = NULL; }
	    if (fill_R) { free_banded_hb_vjd_deck(Ralpha[y], i0, j0, y, cp9b); Ralpha[y] = NULL;
	                  free_banded_hb_vjd_deck(Ralpha[z], i0, j0, z, cp9b); Ralpha[z] = NULL; }
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0) {
		  free_banded_hb_vjd_deck(alpha[y], i0, j0, y, cp9b); alpha[y] = NULL;
		  if (fill_L) { free_banded_hb_vjd_deck(Lalpha[y], i0, j0, y, cp9b); Lalpha[y] = NULL; }
		  if (fill_R) { free_banded_hb_vjd_deck(Ralpha[y], i0, j0, y, cp9b); Ralpha[y] = NULL; }
		}
	      }
	  }
      }
  } /* end loop over all v */

  { int dpr; sc = hb_inband(cp9b, vroot, j0, W, i0, j0, &dpr) ? alpha[vroot][j0][dpr] : IMPOSSIBLE; }
  if (ret_b != NULL)   *ret_b   = b;
  if (ret_bsc != NULL) *ret_bsc = bsc;

  if (ret_alpha == NULL) {
    for (v = vroot; v <= vend; v++)
      if (alpha[v] != NULL) { free_banded_hb_vjd_deck(alpha[v], i0, j0, v, cp9b); alpha[v] = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  /* L/R score planes: when ret_planes (splitter), KEEP them and return via the lr
   * bundle; otherwise (insideT base case, which uses only the shadows) free any
   * decks still live (do_full kept them, or this was the top problem). */
  if (fill_L && ! ret_planes) {
    for (v = vroot; v <= vend; v++) if (Lalpha[v] != NULL) { free_banded_hb_vjd_deck(Lalpha[v], i0, j0, v, cp9b); Lalpha[v] = NULL; }
    free(Lalpha); Lalpha = NULL;
  }
  if (fill_R && ! ret_planes) {
    for (v = vroot; v <= vend; v++) if (Ralpha[v] != NULL) { free_banded_hb_vjd_deck(Ralpha[v], i0, j0, v, cp9b); Ralpha[v] = NULL; }
    free(Ralpha); Ralpha = NULL;
  }

  if (ret_dpool == NULL) deckpool_free(dpool);
  else                   *ret_dpool = dpool;

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  /* T (brief 26_0610-051): Tfull was consumed into Tbsc above; free it always. Tfullk (the
   * full-span k* per B state) is needed by tr_insideT_hb's traceback, so return it
   * when ret_shadow != NULL (the traceback path); otherwise free it. */
  if (Tfull != NULL) { free(Tfull); Tfull = NULL; }
  if (Tfullk != NULL && ret_shadow == NULL) { free(Tfullk); Tfullk = NULL; }
  /* hand the L/R shadows (+ B-state kmode) and/or kept score planes back via lr */
  if (lr != NULL) {
    lr->Lshad = Lshad; lr->Rshad = Rshad; lr->Lkmode = Lkmode; lr->Rkmode = Rkmode;
    lr->Lb = Lb; lr->Lbsc = Lbsc; lr->Rb = Rb; lr->Rbsc = Rbsc;
    lr->Lalpha = Lalpha; lr->Ralpha = Ralpha;
    lr->Tb = Tb; lr->Tbsc = Tbsc; lr->Tfullk = (ret_shadow != NULL) ? Tfullk : NULL;
  }
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.;
}

/* Function: tr_outside_hb()  [brief 26_0610-044, R4.4a J; brief 26_0610-049, R4.4b-pt2c-i: 2-D L/R]
 *
 * Purpose:  J-plane truncated analogue of outside_hb(). Identical banded outside
 *           recurrence, EXCEPT the root deck handling implements the truncated
 *           begin: when vroot==0 the normal root seed beta[0][j0][W]=0 is OMITTED
 *           (no normal ROOT_S descent), and every cp9b->Jvalid[] state in the
 *           subgraph gets a penalty-aware truncated-begin entry beta[v][j0][W] =
 *           tr_trpenalty(v) (replacing the local-begin cm->beginsc[v] injection),
 *           which then propagates down. This makes beta[v] carry "begin into an
 *           ancestor -> descend to v", exactly the oracle's semantics.
 *
 *           Brief 26_0610-049 (R4.4b-part-2c-i): the L/R marginal outside is now stored as
 *           full 2-D banded vjd decks betaL[v][j][d] / betaR[v][j][d] (NOT the 1-D
 *           rows brief 26_0610-046 used). The L/R recurrence is transcribed cell-for-cell
 *           from the canonical banded oracle cm_TrCYKOutsideAlignHB()
 *           (cm_dpalign_trunc.c:6952-7102): the marginal-mode boundary gates are the
 *           subproblem boundaries j==j0 (L) / i==i0 (R) -- the D&C sub-region
 *           restriction of the oracle's global j==L / i==1. The 2-D form is
 *           band-consistent by construction (it can represent "j==j0", which the
 *           1-D betaL[v][j][dp_v] could not -> 047's mixed-split over-score is gone).
 *           Decks are banded vjd and free on the same touch-count loop as the J
 *           beta, so the live frontier stays O(log N). T (brief 26_0610-051) needs NO
 *           outside here: its only "outside" is the begin scalar at a B's full-span
 *           cell (handled in the splitter), matching oracle cm_TrCYKOutsideAlignHB
 *           where Tbeta[v][L][L]=trpenalty with no T recurrence.
 *           MARGINAL LOCAL-END (Lbeta/Rbeta at deck M): brief 26_0610-090 BUILT it
 *           (was the brief 26_0610-049 deferral; brief 26_0610-088 root-caused the 2/1800
 *           matl300 L-mode undershoots 086 flagged as this gap). betaL[cm->M]/betaR[cm->M]
 *           are now allocated/seeded/propagated below with per-mode marginal emission
 *           transcribed cell-for-cell from the monolithic Lelbeta/Relbeta
 *           (cm_dpalign_trunc.c:2152-2226) + inline el_selfsc*d to match the J feed; the
 *           wedge/generic splitters carry L/R EL candidates, and the best_v=-1 EL
 *           traceback routes through tr_v_splitter_hb->tr_vinsideT_hb in the winning
 *           parent mode (already EL+mode aware). RESULT: sample5_3ptr now byte-exact vs
 *           the oracle (153-node parse; the brief-088 proof betaL[420][j][d+1]+endsc[420]+
 *           el_selfsc*d+esc = -6.18375 is exactly what this deck now supplies). NOTE:
 *           sample9_5ptr is NOT this gap -- it is a separate PRE-EXISTING inside-begin
 *           inflation (D&C inside-L[606]=-13.436 vs oracle Ldp[606]=-15.175), unchanged
 *           by this fix (own follow-on brief). Valid for CMH_LOCAL_END on and off.
 */
static void
tr_outside_hb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
	      int do_full, int fill_L, int fill_R,
	      float ***beta, float ****ret_beta,
	      float ***betaL, float ****ret_betaL,
	      float ***betaR, float ****ret_betaR,
	      struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	      float *ret_bsc, int *ret_bv, int *ret_bj, int *ret_bmode, CP9Bands_t *cp9b)
{
  int      status;
  int      v,y;
  int      j,d,i;
  float    sc;
  int     *touch;
  float    escore;
  int      W;
  int      jp;
  int      jp_v;
  int      voffset;
  int      w1,w2;
  int      jn, jx;
  float    b_sc   = IMPOSSIBLE;  /* best marginal terminus ("local hit in parent")    */
  int      b_v    = -1;
  int      b_j    = -1;
  int      b_mode = -1;
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  int     *eldmaxJ = NULL, *eldmaxL = NULL, *eldmaxR = NULL;  /* brief 26_0430-227: banded EL d-band edges (J/L/R) */

  W = j0-i0+1;
  if (dpool == NULL) dpool = deckpool_create();

  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }
  /* L/R marginal beta pointer arrays (2-D banded vjd decks per state, brief 26_0610-049).
   * Reuse caller's if any (chaining), else allocate fresh; decks recycle with the
   * J decks (touch counting). */
  if (fill_L && betaL == NULL) { ESL_ALLOC(betaL, sizeof(float **) * (cm->M+1)); for (v = 0; v <= cm->M; v++) betaL[v] = NULL; }
  if (fill_R && betaR == NULL) { ESL_ALLOC(betaR, sizeof(float **) * (cm->M+1)); for (v = 0; v <= cm->M; v++) betaR[v] = NULL; }

  w1 = cm->nodemap[cm->ndidx[vroot]];
  if (cm->sttype[vroot] == B_st) {
    w2 = w1;
    if (vend != vroot) cm_Fail("oh no. not again.");
  } else if (cm->sttype[vroot] == IL_st || cm->sttype[vroot] == IR_st) {
    /* TRUNCATED begins (unlike classical CMH_LOCAL_BEGIN) may enter directly at
     * an insert state -- a fragment can start/end mid-insertion (brief 26_0610-069 fix).
     * An insert state is not part of any node's split set (w1 above is its
     * node's first SPLIT state, not vroot itself), so the split-set grouping
     * below doesn't apply; vroot has no split-set siblings to allocate here,
     * only its own single deck. */
    w1 = w2 = vroot;
  } else
    w2 = cm->cfirst[w1]-1;

  for (v = w1; v <= w2; v++) {
    beta[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
    for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
      jp_v = j - jmin[v];
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) beta[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE;
    }
    if (fill_L) { betaL[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
      for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) { jp_v = j - jmin[v];
        for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) betaL[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE; } }
    if (fill_R) { betaR[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
      for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) { jp_v = j - jmin[v];
        for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) betaR[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE; } }
  }
  /* TRUNCATED: omit the normal root seed (beta[0][j0][W]=0) when vroot==0 -- in
   * truncated mode there is no normal ROOT_S descent. For interior subproblems
   * (vroot!=0) keep the standard root seed (J + the marginal-mode root seeds). The
   * L/R marginal root cell is (j0,W) -- same cell as J (the full-span subproblem
   * corner), the 2-D analogue of the 1-D betaL[vroot][i0]/betaR[vroot][j0]. */
  if (vroot != 0) {
    int dpr; if (hb_inband(cp9b, vroot, j0, W, i0, j0, &dpr)) {
      beta[vroot][j0][dpr] = 0;
      if (fill_L && cp9b->Lvalid[vroot]) betaL[vroot][j0][dpr] = 0.;
      if (fill_R && cp9b->Rvalid[vroot]) betaR[vroot][j0][dpr] = 0.;
    }
  }

  if (cm->flags & CMH_LOCAL_END) {
    /* brief 26_0430-227: allocate the three EL decks (J/L/R) BANDED on the upper
     * d-edge (eldmax*[j]) instead of full O(W^2) triangles, and never pool them
     * (band-dependent size is incompatible with the shared deckpool).  Only the
     * banded (d=0..eldmax*[j]) cells are IMPOSSIBLE-inited.  The three eldmax
     * arrays are recomputed (deterministically) by the splitter readers below. */
    ESL_ALLOC(eldmaxJ, sizeof(int) * (L+1));
    tr_outside_hb_el_dmax(cm, L, vroot, vend, i0, j0, 0, cp9b, eldmaxJ);
    beta[cm->M] = alloc_el_banded_vjd_deck(L, i0, j0, eldmaxJ);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      if (eldmaxJ[j] < 0) continue;
      for (d = 0; d <= eldmaxJ[j]; d++) beta[cm->M][j][d] = IMPOSSIBLE;
    }
    /* brief 26_0610-090: marginal L/R local-end (EL) outside decks, now banded
     * (brief 26_0430-227).  Built only when the corresponding marginal plane is
     * requested (fill_L/fill_R). */
    if (fill_L) {
      ESL_ALLOC(eldmaxL, sizeof(int) * (L+1));
      tr_outside_hb_el_dmax(cm, L, vroot, vend, i0, j0, 1, cp9b, eldmaxL);
      betaL[cm->M] = alloc_el_banded_vjd_deck(L, i0, j0, eldmaxL);
      for (jp = 0; jp <= W; jp++) { j = i0-1+jp; if (eldmaxL[j] < 0) continue; for (d = 0; d <= eldmaxL[j]; d++) betaL[cm->M][j][d] = IMPOSSIBLE; }
    }
    if (fill_R) {
      ESL_ALLOC(eldmaxR, sizeof(int) * (L+1));
      tr_outside_hb_el_dmax(cm, L, vroot, vend, i0, j0, 2, cp9b, eldmaxR);
      betaR[cm->M] = alloc_el_banded_vjd_deck(L, i0, j0, eldmaxR);
      for (jp = 0; jp <= W; jp++) { j = i0-1+jp; if (eldmaxR[j] < 0) continue; for (d = 0; d <= eldmaxR[j]; d++) betaR[cm->M][j][d] = IMPOSSIBLE; }
    }
    if (vroot != 0 && NOT_IMPOSSIBLE(cm->endsc[vroot])) {
      switch (cm->sttype[vroot]) {
      case MP_st:
	if (W < 2) break;
	if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
	else
	  escore = DegeneratePairScore(cm->abc, cm->esc[vroot], dsq[i0], dsq[j0]);
	beta[cm->M][j0-1][W-2] = cm->endsc[vroot] + (cm->el_selfsc * (W-2)) + escore;
	if (beta[cm->M][j0-1][W-2] < IMPOSSIBLE) beta[cm->M][j0-1][W-2] = IMPOSSIBLE;
	break;
      case ML_st:
      case IL_st:
	if (W < 1) break;
	if (dsq[i0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) dsq[i0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[vroot]);
	beta[cm->M][j0][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + escore;
	if (beta[cm->M][j0][W-1] < IMPOSSIBLE) beta[cm->M][j0][W-1] = IMPOSSIBLE;
	break;
      case MR_st:
      case IR_st:
	if (W < 1) break;
	if (dsq[j0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) dsq[j0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[vroot]);
	beta[cm->M][j0-1][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + escore;
	if (beta[cm->M][j0-1][W-1] < IMPOSSIBLE) beta[cm->M][j0-1][W-1] = IMPOSSIBLE;
	break;
      case S_st:
      case D_st:
	beta[cm->M][j0][W] = cm->endsc[vroot] + (cm->el_selfsc * W);
	if (beta[cm->M][j0][W] < IMPOSSIBLE) beta[cm->M][j0][W] = IMPOSSIBLE;
	break;
      case B_st:
      default: cm_Fail("bogus parent state %d\n", cm->sttype[vroot]);
      }
      /* brief 26_0610-090: marginal vroot EL seeds -- mathematically the L/R EL
       * propagate feed (below) applied to vroot at its full-span root cell
       * (betaL/betaR[vroot][j0][W]=0, seeded above), which the propagate loop
       * (v=w2+1..vend) skips for vroot. Per-mode marginal emission, derived by the
       * same rule as the propagate: L emits only the LEFT residue (MP->lmesc, ML/IL
       * ->esc), gated to the L boundary j==j0 (full span); R emits only the RIGHT
       * residue (MP->rmesc, MR/IR->esc), gated to the R boundary i==i0 (full span).
       * MR/IR in L (and ML/IL in R) emit nothing (their residue is truncated away). */
      if (fill_L && cp9b->Lvalid[vroot]) {
	switch (cm->sttype[vroot]) {
	case MP_st:
	  if (W < 1) break;
	  betaL[cm->M][j0][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + cm->lmesc[vroot][dsq[i0]];
	  if (betaL[cm->M][j0][W-1] < IMPOSSIBLE) betaL[cm->M][j0][W-1] = IMPOSSIBLE;
	  break;
	case ML_st:
	case IL_st:
	  if (W < 1) break;
	  if (dsq[i0] < cm->abc->K) escore = cm->esc[vroot][(int) dsq[i0]];
	  else                      escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[vroot]);
	  betaL[cm->M][j0][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + escore;
	  if (betaL[cm->M][j0][W-1] < IMPOSSIBLE) betaL[cm->M][j0][W-1] = IMPOSSIBLE;
	  break;
	case MR_st:
	case IR_st:
	case S_st:
	case D_st:
	  betaL[cm->M][j0][W] = cm->endsc[vroot] + (cm->el_selfsc * W);
	  if (betaL[cm->M][j0][W] < IMPOSSIBLE) betaL[cm->M][j0][W] = IMPOSSIBLE;
	  break;
	case B_st:
	default: cm_Fail("bogus parent state %d\n", cm->sttype[vroot]);
	}
      }
      if (fill_R && cp9b->Rvalid[vroot]) {
	switch (cm->sttype[vroot]) {
	case MP_st:
	  if (W < 1) break;
	  betaR[cm->M][j0-1][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + cm->rmesc[vroot][dsq[j0]];
	  if (betaR[cm->M][j0-1][W-1] < IMPOSSIBLE) betaR[cm->M][j0-1][W-1] = IMPOSSIBLE;
	  break;
	case MR_st:
	case IR_st:
	  if (W < 1) break;
	  if (dsq[j0] < cm->abc->K) escore = cm->esc[vroot][(int) dsq[j0]];
	  else                      escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[vroot]);
	  betaR[cm->M][j0-1][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + escore;
	  if (betaR[cm->M][j0-1][W-1] < IMPOSSIBLE) betaR[cm->M][j0-1][W-1] = IMPOSSIBLE;
	  break;
	case ML_st:
	case IL_st:
	case S_st:
	case D_st:
	  betaR[cm->M][j0][W] = cm->endsc[vroot] + (cm->el_selfsc * W);
	  if (betaR[cm->M][j0][W] < IMPOSSIBLE) betaR[cm->M][j0][W] = IMPOSSIBLE;
	  break;
	case B_st:
	default: cm_Fail("bogus parent state %d\n", cm->sttype[vroot]);
	}
      }
    }
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;      v < w1; v++) touch[v] = 0;
  for (v = vend+1; v < cm->M; v++) touch[v] = 0;
  for (v = w1; v <= vend; v++) {
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
  }

  for (v = w2+1; v <= vend; v++)
    {
      beta[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);

      for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
	jp_v = j - jmin[v];
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) beta[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE;
      }
      if (fill_L) { betaL[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
	for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) { jp_v = j - jmin[v];
	  for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) betaL[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE; } }
      if (fill_R) { betaR[v] = alloc_banded_hb_vjd_deck(L, i0, j0, v, cp9b);
	for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) { jp_v = j - jmin[v];
	  for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v); d++) betaR[v][j][d - hd_min(cp9b, v, jp_v)] = IMPOSSIBLE; } }

      /* TRUNCATED-BEGIN injection: at the top (vroot==0, full span i0=1,j0=L) any
       * Jvalid state v may be entered via a truncated begin with penalty
       * tr_trpenalty(v). (Replaces the non-trunc CMH_LOCAL_BEGIN / cm->beginsc[v]
       * injection; unconditional on local mode -- truncated begins always apply.)
       * L/R (brief 26_0610-049): inject into the 2-D marginal root cell (j0,W) -- the same
       * full-span corner as J (oracle cm_TrCYKOutsideAlignHB:6739-6741) -- for each
       * {L,R}valid v.
       *
       * The injection itself is correct (byte-exact vs the oracle's own
       * beta[v][j0][W]=trpen; brief 26_0610-086). The 081-reported LOCAL-BEGIN
       * inflation was NOT here: it was the downward-propagation loop just below
       * chaining this trpen THROUGH a non-J-valid (cp9b->Jvalid[v]==FALSE) wedge
       * state (the injection was already Jvalid-gated, but the J-plane
       * propagation was not) -- fixed in brief 26_0610-086 by the do_J_v gate
       * below (and its traceback twin in tr_vinside_hb). Do NOT "fix" by simply
       * disabling this injection -- that regresses legitimate wedge-entry
       * candidates (081 confirmed experimentally). */
      if (vroot == 0 && i0 == 1 && j0 == L &&
	  (jmin[v] <= j0 && jmax[v] >= j0)
	  && (hd_min(cp9b, v, j0-jmin[v]) <= W && hd_max(cp9b, v, j0-jmin[v]) >= W)) {
	int dpW = W - hd_min(cp9b, v, j0-jmin[v]);
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen)) {
	  if (cp9b->Jvalid[v]            && trpen > beta[v][j0][dpW])  beta[v][j0][dpW]  = trpen;
	  if (fill_L && cp9b->Lvalid[v]  && trpen > betaL[v][j0][dpW]) betaL[v][j0][dpW] = trpen;
	  if (fill_R && cp9b->Rvalid[v]  && trpen > betaR[v][j0][dpW]) betaR[v][j0][dpW] = trpen;
	}
      }

      jn = ESL_MAX(i0-1, jmin[v]);
      jx = ESL_MIN(j0,   jmax[v]);
      for (j = jx; j >= jn; j--) {
	jp   = j - (i0-1);
	jp_v = j - jmin[v];
	for (d = ESL_MIN(hd_max(cp9b, v, jp_v), jp); d >= hd_min(cp9b, v, jp_v); d--)
	  {
	    int dp_v = d - hd_min(cp9b, v, jp_v);
	    int dp_y;
	    /* brief 26_0610-086: gate the pure-J outside propagation on cp9b->Jvalid[v],
	     * mirroring the do_L_v/do_R_v gates below (and the Jvalid[v] gates already
	     * present on this block's marginal->J cross-terms + the truncated-begin
	     * injection above). Without it, beta[v]'s J-plane is filled for a state v
	     * whose node structurally cannot appear in a pure-J parse (cp9b->Jvalid[v]
	     * ==FALSE), letting a truncated-begin penalty propagate a chain of (good)
	     * emissions THROUGH that non-J-valid state -- inflating tr_generic_splitter_hb's
	     * "all-J split" (beta[v]-based) candidate to an unreachable score. The oracle
	     * (cm_TrCYKInsideAlignHB/cm_TrCYKOutsideAlignHB) never allocates a J deck for
	     * a !Jvalid state, so its beta stays IMPOSSIBLE; do_J_v reproduces that mask. */
	    int do_J_v = cp9b->Jvalid[v];
	    int do_L_v = fill_L && cp9b->Lvalid[v];
	    int do_R_v = fill_R && cp9b->Rvalid[v];
	    i = j-d+1;
	    for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	      if (y < vroot) continue;
	      voffset = v - cm->cfirst[y];

	      /* Brief 26_0610-049: J keeps the exact R4.4a recurrence; L/R now fill the
	       * 2-D banded decks betaL[v][j][d]/betaR[v][j][d], transcribed cell-
	       * for-cell from cm_TrCYKOutsideAlignHB (cm_dpalign_trunc.c:7000-7093).
	       * Marginal boundary gates are j==j0 (L) / i==i0 (R) -- the D&C sub-
	       * region restriction of the oracle's global j==L / i==1. */
	      switch(cm->sttype[y]) {
	      case MP_st:
		/* J <- J : gate j!=j0 && d!=jp ; y at (j+1,d+2) ; pair esc */
		if (j != j0 && d != jp) {
		  if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		    escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		  else
		    escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
		  if (do_J_v && hb_inband(cp9b, y, j+1, d+2, i0, j0, &dp_y) &&
		    (sc = beta[y][j+1][dp_y] + cm->tsc[y][voffset] + escore) > beta[v][j][dp_v])
		  beta[v][j][dp_v] = sc;
		}
		/* L branch (MP parent emits LEFT in L mode): gate j==j0 && d!=jp ; y at (j,d+1) ; lmesc */
		if (fill_L && j == j0 && d != jp) {
		  int dpyL;
		  if (hb_inband(cp9b, y, j, d+1, i0, j0, &dpyL) && NOT_IMPOSSIBLE(betaL[y][j][dpyL])) {
		    sc = betaL[y][j][dpyL] + cm->tsc[y][voffset] + cm->lmesc[y][dsq[i-1]];
		    if (cp9b->Jvalid[v] && sc > beta[v][j][dp_v])  beta[v][j][dp_v]  = sc;
		    if (do_L_v          && sc > betaL[v][j][dp_v]) betaL[v][j][dp_v] = sc;
		  }
		}
		/* R branch (MP parent emits RIGHT in R mode): gate i==i0 && j!=j0 ; y at (j+1,d+1) ; rmesc */
		if (fill_R && i == i0 && j != j0) {
		  int dpyR;
		  if (hb_inband(cp9b, y, j+1, d+1, i0, j0, &dpyR) && NOT_IMPOSSIBLE(betaR[y][j+1][dpyR])) {
		    sc = betaR[y][j+1][dpyR] + cm->tsc[y][voffset] + cm->rmesc[y][dsq[j+1]];
		    if (cp9b->Jvalid[v] && sc > beta[v][j][dp_v])  beta[v][j][dp_v]  = sc;
		    if (do_R_v          && sc > betaR[v][j][dp_v]) betaR[v][j][dp_v] = sc;
		  }
		}
		break;
	      case ML_st:
	      case IL_st:
		/* J<-J and L<-L : gate d!=jp ; y at (j,d+1) ; left esc(i-1) */
		if (d != jp) {
		  if (dsq[i-1] < cm->abc->K)
		    escore = cm->esc[y][(int) dsq[i-1]];
		  else
		    escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		  if (hb_inband(cp9b, y, j, d+1, i0, j0, &dp_y)) {
		    if (do_J_v && (sc = beta[y][j][dp_y] + cm->tsc[y][voffset] + escore) > beta[v][j][dp_v])
		      beta[v][j][dp_v] = sc;
		    if (do_L_v && NOT_IMPOSSIBLE(betaL[y][j][dp_y]) &&
		      (sc = betaL[y][j][dp_y] + cm->tsc[y][voffset] + escore) > betaL[v][j][dp_v])
		      betaL[v][j][dp_v] = sc;
		  }
		}
		/* R<-R (no emit): gate i==i0 && v!=y (IL self-transit guard) ; y at (j,d) */
		if (fill_R && i == i0 && v != y) {
		  int dpyR;
		  if (hb_inband(cp9b, y, j, d, i0, j0, &dpyR) && NOT_IMPOSSIBLE(betaR[y][j][dpyR])) {
		    sc = betaR[y][j][dpyR] + cm->tsc[y][voffset];
		    if (cp9b->Jvalid[v] && sc > beta[v][j][dp_v])  beta[v][j][dp_v]  = sc;
		    if (do_R_v          && sc > betaR[v][j][dp_v]) betaR[v][j][dp_v] = sc;
		  }
		}
		break;
	      case MR_st:
	      case IR_st:
		/* J<-J and R<-R : gate j!=j0 ; y at (j+1,d+1) ; right esc(j+1) */
		if (j != j0) {
		  if (dsq[j+1] < cm->abc->K)
		    escore = cm->esc[y][(int) dsq[j+1]];
		  else
		    escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);
		  if (hb_inband(cp9b, y, j+1, d+1, i0, j0, &dp_y)) {
		    if (do_J_v && (sc = beta[y][j+1][dp_y] + cm->tsc[y][voffset] + escore) > beta[v][j][dp_v])
		      beta[v][j][dp_v] = sc;
		    if (do_R_v && NOT_IMPOSSIBLE(betaR[y][j+1][dp_y]) &&
		      (sc = betaR[y][j+1][dp_y] + cm->tsc[y][voffset] + escore) > betaR[v][j][dp_v])
		      betaR[v][j][dp_v] = sc;
		  }
		}
		/* L<-L (no emit): gate j==j0 && v!=y (IR self-transit guard) ; y at (j,d) */
		if (fill_L && j == j0 && v != y) {
		  int dpyL;
		  if (hb_inband(cp9b, y, j, d, i0, j0, &dpyL) && NOT_IMPOSSIBLE(betaL[y][j][dpyL])) {
		    sc = betaL[y][j][dpyL] + cm->tsc[y][voffset];
		    if (cp9b->Jvalid[v] && sc > beta[v][j][dp_v])  beta[v][j][dp_v]  = sc;
		    if (do_L_v          && sc > betaL[v][j][dp_v]) betaL[v][j][dp_v] = sc;
		  }
		}
		break;
	      case S_st:
	      case E_st:
	      case D_st:
		if (! hb_inband(cp9b, y, j, d, i0, j0, &dp_y)) continue;
		if (do_J_v && (sc = beta[y][j][dp_y] + cm->tsc[y][voffset]) > beta[v][j][dp_v])
		  beta[v][j][dp_v] = sc;
		if (do_L_v && NOT_IMPOSSIBLE(betaL[y][j][dp_y]) &&
		  (sc = betaL[y][j][dp_y] + cm->tsc[y][voffset]) > betaL[v][j][dp_v])
		  betaL[v][j][dp_v] = sc;
		if (do_R_v && NOT_IMPOSSIBLE(betaR[y][j][dp_y]) &&
		  (sc = betaR[y][j][dp_y] + cm->tsc[y][voffset]) > betaR[v][j][dp_v])
		  betaR[v][j][dp_v] = sc;
		break;
	      default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
	      }
	    }
	    if (beta[v][j][dp_v] < IMPOSSIBLE) beta[v][j][dp_v] = IMPOSSIBLE;
	    if (do_L_v && betaL[v][j][dp_v] < IMPOSSIBLE) betaL[v][j][dp_v] = IMPOSSIBLE;
	    if (do_R_v && betaR[v][j][dp_v] < IMPOSSIBLE) betaR[v][j][dp_v] = IMPOSSIBLE;
	  }
      }

      /* Marginal terminus b_sc (2-D, brief 26_0610-049): the "local hit in parent" -- the
       * optimal marginal parse stays within r..v and ends at v via a TRUNCATED END
       * (the child bifurcation/wedge is empty). The truncated-end inside base is the
       * d==1 single-emission cell (tr_inside_hb USED_TRUNC_END, :6775/:6795): for an
       * emitting v in the marginal mode, inside = v's marginal emission. So the
       * terminus is beta{L,R}[v][j][d==1] + that emission, scanned over j, for
       * left/right-emitting v only (non-emitters get IMPOSSIBLE truncated-end). The
       * d==1 restriction is what the 1-D row implicitly encoded; an earlier all-d
       * scan over-claimed -> bifurcated-5S R over-score. */
      if (fill_L && cp9b->Lvalid[v] && (cm->sttype[v]==MP_st || cm->sttype[v]==ML_st || cm->sttype[v]==IL_st))
	for (j = ESL_MAX(i0, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
	  int dpv1;
	  if (! hb_inband(cp9b, v, j, 1, i0, j0, &dpv1)) continue;
	  if (! NOT_IMPOSSIBLE(betaL[v][j][dpv1])) continue;
	  escore = (cm->sttype[v]==MP_st) ? cm->lmesc[v][dsq[j]]
	         : ((dsq[j] < cm->abc->K) ? cm->esc[v][(int) dsq[j]] : esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]));
	  if (betaL[v][j][dpv1] + escore > b_sc) { b_sc = betaL[v][j][dpv1] + escore; b_v = v; b_j = j; b_mode = TRMODE_L; }
	}
      if (fill_R && cp9b->Rvalid[v] && (cm->sttype[v]==MP_st || cm->sttype[v]==MR_st || cm->sttype[v]==IR_st))
	for (j = ESL_MAX(i0, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
	  int dpv1;
	  if (! hb_inband(cp9b, v, j, 1, i0, j0, &dpv1)) continue;
	  if (! NOT_IMPOSSIBLE(betaR[v][j][dpv1])) continue;
	  escore = (cm->sttype[v]==MP_st) ? cm->rmesc[v][dsq[j]]
	         : ((dsq[j] < cm->abc->K) ? cm->esc[v][(int) dsq[j]] : esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]));
	  if (betaR[v][j][dpv1] + escore > b_sc) { b_sc = betaR[v][j][dpv1] + escore; b_v = v; b_j = j; b_mode = TRMODE_R; }
	}
      if (NOT_IMPOSSIBLE(cm->endsc[v])) {
	int dp_v;
	/* brief 26_0610-101: band-limit the v->EL feed to v's OWN banded footprint.
	 * The old loop swept the full 0..W x 0..jp triangle (O(W^2)) for EVERY state v
	 * with a local end -- the dominant cost of local-mode D&C at genome scale
	 * (~88-98% of all EL-deck work, brief 26_0610-099/100/101). But a cell (j,d) can
	 * only update beta[cm->M] when v's read cell -- shifted by the per-type (elsj,elsd)
	 * below -- is IN v's band; every out-of-band read hit `!hb_inband -> continue`
	 * and did nothing. So we iterate ONLY v's band rows (elJ = read row) and their
	 * shifted d-range. The body (switch, hb_inband guard, boundary gates, all index
	 * math) is BYTE-IDENTICAL to the old sweep; the tighter bounds are a provable
	 * superset of the old update set, so the result is byte-exact (Task C). Read-cell
	 * shifts: MP=(j+1,d+2) ML/IL=(j,d+1) MR/IR=(j+1,d+1) S/D/E=(j,d). */
	int elsj = 0, elsd = 0, elJ, elJlo, elJhi;
	switch (cm->sttype[v]) {
	case MP_st:                       elsj = 1; elsd = 2; break;
	case ML_st: case IL_st:           elsj = 0; elsd = 1; break;
	case MR_st: case IR_st:           elsj = 1; elsd = 1; break;
	case S_st:  case D_st: case E_st: elsj = 0; elsd = 0; break;
	case B_st:
	default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	}
	elJlo = ESL_MAX(i0-1, jmin[v]); elJhi = ESL_MIN(j0, jmax[v]);
	for (elJ = elJlo; elJ <= elJhi; elJ++) {
	  int eljpv = elJ - jmin[v], eldlo, eldhi;
	  j  = elJ - elsj;
	  jp = j - (i0-1);
	  if (jp < 0) continue;
	  eldlo = hd_min(cp9b, v, eljpv) - elsd; if (eldlo < 0)  eldlo = 0;
	  eldhi = hd_max(cp9b, v, eljpv) - elsd; if (eldhi > jp) eldhi = jp;
	  for (d = eldlo; d <= eldhi; d++)
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st:
		if (j == j0 || d == jp) continue;
		if (! hb_inband(cp9b, v, j+1, d+2, i0, j0, &dp_v)) continue;
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		if ((sc = beta[v][j+1][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;
		if (! hb_inband(cp9b, v, j, d+1, i0, j0, &dp_v)) continue;
		if (dsq[i-1] < cm->abc->K)
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		if ((sc = beta[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (! hb_inband(cp9b, v, j+1, d+1, i0, j0, &dp_v)) continue;
		if (dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		if ((sc = beta[v][j+1][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if (! hb_inband(cp9b, v, j, d, i0, j0, &dp_v)) continue;
		if ((sc = beta[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d)) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case B_st:
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	      }
	    }
	}
      }

      /* brief 26_0610-090: L-marginal v->EL feed into betaL[cm->M]. Transcribed
       * from the trusted monolithic Lelbeta recurrence (cm_dpalign_trunc.c:2152-2188,
       * itself cm_TrOutsideAlignHB's v->EL block) into the D&C's own (j,d)/hb_inband
       * frame, with el_selfsc*d applied INLINE (matching the J feed above; the
       * monolithic defers this self-loop fold). L mode is right-truncated: a
       * left-emitting parent (MP/ML/IL) emits only its LEFT residue and feeds betaL[v]
       * at (j,d+1); a right-emitting parent (MR/IR) emits nothing and feeds betaL[v]
       * at (j,d) gated to the L boundary j==j0; MP also carries j==j0. brief 26_0610-088
       * proved the MP-absent left-emitting case byte-exact. */
      if (fill_L && cp9b->Lvalid[v] && NOT_IMPOSSIBLE(cm->endsc[v])) {
	int dp_v;
	/* brief 26_0610-101: band-limit the L-marginal v->EL feed (see the J-feed note
	 * above). Byte-identical body; read-cell shifts here: MP=(j,d+1) ML/IL=(j,d+1)
	 * MR/IR=(j,d) S/D/E=(j,d). */
	int elsj = 0, elsd = 0, elJ, elJlo, elJhi;
	switch (cm->sttype[v]) {
	case MP_st:                       elsj = 0; elsd = 1; break;
	case ML_st: case IL_st:           elsj = 0; elsd = 1; break;
	case MR_st: case IR_st:           elsj = 0; elsd = 0; break;
	case S_st:  case D_st: case E_st: elsj = 0; elsd = 0; break;
	case B_st:
	default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	}
	elJlo = ESL_MAX(i0-1, jmin[v]); elJhi = ESL_MIN(j0, jmax[v]);
	for (elJ = elJlo; elJ <= elJhi; elJ++) {
	  int eljpv = elJ - jmin[v], eldlo, eldhi;
	  j  = elJ - elsj;
	  jp = j - (i0-1);
	  if (jp < 0) continue;
	  eldlo = hd_min(cp9b, v, eljpv) - elsd; if (eldlo < 0)  eldlo = 0;
	  eldhi = hd_max(cp9b, v, eljpv) - elsd; if (eldhi > jp) eldhi = jp;
	  for (d = eldlo; d <= eldhi; d++)
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st:
		if (j != j0 || d == jp) continue;
		if (! hb_inband(cp9b, v, j, d+1, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaL[v][j][dp_v])) continue;
		if ((sc = betaL[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + cm->lmesc[v][dsq[i-1]]) > betaL[cm->M][j][d])
		  betaL[cm->M][j][d] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;
		if (! hb_inband(cp9b, v, j, d+1, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaL[v][j][dp_v])) continue;
		if (dsq[i-1] < cm->abc->K) escore = cm->esc[v][(int) dsq[i-1]];
		else                       escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		if ((sc = betaL[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > betaL[cm->M][j][d])
		  betaL[cm->M][j][d] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j != j0) continue;
		if (! hb_inband(cp9b, v, j, d, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaL[v][j][dp_v])) continue;
		if ((sc = betaL[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d)) > betaL[cm->M][j][d])
		  betaL[cm->M][j][d] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if (! hb_inband(cp9b, v, j, d, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaL[v][j][dp_v])) continue;
		if ((sc = betaL[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d)) > betaL[cm->M][j][d])
		  betaL[cm->M][j][d] = sc;
		break;
	      case B_st:
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	      }
	    }
	}
      }

      /* brief 26_0610-090: R-marginal v->EL feed into betaR[cm->M]. Transcribed from
       * the trusted monolithic Relbeta recurrence (cm_dpalign_trunc.c:2190-2226).
       * R mode is left-truncated: a right-emitting parent (MP/MR/IR) emits only its
       * RIGHT residue and feeds betaR[v] at (j+1,d+1); a left-emitting parent (ML/IL)
       * emits nothing and feeds betaR[v] at (j,d) gated to the R boundary i==i0; MP
       * also carries i==i0. */
      if (fill_R && cp9b->Rvalid[v] && NOT_IMPOSSIBLE(cm->endsc[v])) {
	int dp_v;
	/* brief 26_0610-101: band-limit the R-marginal v->EL feed (see the J-feed note
	 * above). Byte-identical body; read-cell shifts here: MP=(j+1,d+1) ML/IL=(j,d)
	 * MR/IR=(j+1,d+1) S/D/E=(j,d). */
	int elsj = 0, elsd = 0, elJ, elJlo, elJhi;
	switch (cm->sttype[v]) {
	case MP_st:                       elsj = 1; elsd = 1; break;
	case ML_st: case IL_st:           elsj = 0; elsd = 0; break;
	case MR_st: case IR_st:           elsj = 1; elsd = 1; break;
	case S_st:  case D_st: case E_st: elsj = 0; elsd = 0; break;
	case B_st:
	default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	}
	elJlo = ESL_MAX(i0-1, jmin[v]); elJhi = ESL_MIN(j0, jmax[v]);
	for (elJ = elJlo; elJ <= elJhi; elJ++) {
	  int eljpv = elJ - jmin[v], eldlo, eldhi;
	  j  = elJ - elsj;
	  jp = j - (i0-1);
	  if (jp < 0) continue;
	  eldlo = hd_min(cp9b, v, eljpv) - elsd; if (eldlo < 0)  eldlo = 0;
	  eldhi = hd_max(cp9b, v, eljpv) - elsd; if (eldhi > jp) eldhi = jp;
	  for (d = eldlo; d <= eldhi; d++)
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st:
		if (i != i0 || j == j0) continue;
		if (! hb_inband(cp9b, v, j+1, d+1, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaR[v][j+1][dp_v])) continue;
		if ((sc = betaR[v][j+1][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + cm->rmesc[v][dsq[j+1]]) > betaR[cm->M][j][d])
		  betaR[cm->M][j][d] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (i != i0) continue;
		if (! hb_inband(cp9b, v, j, d, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaR[v][j][dp_v])) continue;
		if ((sc = betaR[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d)) > betaR[cm->M][j][d])
		  betaR[cm->M][j][d] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (! hb_inband(cp9b, v, j+1, d+1, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaR[v][j+1][dp_v])) continue;
		if (dsq[j+1] < cm->abc->K) escore = cm->esc[v][(int) dsq[j+1]];
		else                       escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		if ((sc = betaR[v][j+1][dp_v] + cm->endsc[v] + (cm->el_selfsc * d) + escore) > betaR[cm->M][j][d])
		  betaR[cm->M][j][d] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if (! hb_inband(cp9b, v, j, d, i0, j0, &dp_v)) continue;
		if (! NOT_IMPOSSIBLE(betaR[v][j][dp_v])) continue;
		if ((sc = betaR[v][j][dp_v] + cm->endsc[v] + (cm->el_selfsc * d)) > betaR[cm->M][j][d])
		  betaR[cm->M][j][d] = sc;
		break;
	      case B_st:
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	      }
	    }
	}
      }

      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) {
	    free_banded_hb_vjd_deck(beta[y], i0, j0, y, cp9b); beta[y] = NULL;
	    if (fill_L) { free_banded_hb_vjd_deck(betaL[y], i0, j0, y, cp9b); betaL[y] = NULL; }
	    if (fill_R) { free_banded_hb_vjd_deck(betaR[y], i0, j0, y, cp9b); betaR[y] = NULL; }
	  }
	}
      }
    }

  if (ret_beta == NULL) {
    for (v = w1; v <= vend; v++)
      if (beta[v] != NULL) { free_banded_hb_vjd_deck(beta[v], i0, j0, v, cp9b); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
      free_el_banded_vjd_deck(beta[cm->M], i0, j0, eldmaxJ);   /* brief 26_0430-227: banded EL deck */
      beta[cm->M] = NULL;
    }
    free(beta);
  } else *ret_beta = beta;

  /* L/R marginal 2-D banded decks (brief 26_0610-049): return them (splitter reads them,
   * and frees via free_banded_hb_vjd_matrix which handles the deck M) or free
   * everything here. brief 26_0610-090: deck M (EL) is now allocated for L/R too;
   * brief 26_0430-227: banded (free with free_el_banded_vjd_deck). */
  if (fill_L) {
    if (ret_betaL != NULL) *ret_betaL = betaL;
    else { for (v = w1; v <= vend; v++) if (betaL[v] != NULL) { free_banded_hb_vjd_deck(betaL[v], i0, j0, v, cp9b); betaL[v] = NULL; }
           if ((cm->flags & CMH_LOCAL_END) && betaL[cm->M] != NULL) { free_el_banded_vjd_deck(betaL[cm->M], i0, j0, eldmaxL); betaL[cm->M] = NULL; }
           free(betaL); }
  }
  if (fill_R) {
    if (ret_betaR != NULL) *ret_betaR = betaR;
    else { for (v = w1; v <= vend; v++) if (betaR[v] != NULL) { free_banded_hb_vjd_deck(betaR[v], i0, j0, v, cp9b); betaR[v] = NULL; }
           if ((cm->flags & CMH_LOCAL_END) && betaR[cm->M] != NULL) { free_el_banded_vjd_deck(betaR[cm->M], i0, j0, eldmaxR); betaR[cm->M] = NULL; }
           free(betaR); }
  }

  if (ret_dpool == NULL) deckpool_free(dpool);
  else                   *ret_dpool = dpool;
  free(touch);
  free(eldmaxJ); free(eldmaxL); free(eldmaxR);   /* brief 26_0430-227: NULL-safe */
  if (ret_bsc   != NULL) *ret_bsc   = b_sc;
  if (ret_bv    != NULL) *ret_bv    = b_v;
  if (ret_bj    != NULL) *ret_bj    = b_j;
  if (ret_bmode != NULL) *ret_bmode = b_mode;
  return;
 ERROR:
  cm_Fail("Memory allocation error.");
}

/* Function: tr_insideT_hb()  [brief 26_0610-045, R4.4b]
 *
 * Purpose:  Truncated analogue of insideT_hb(): run tr_inside_hb() filling the
 *           J plane plus whichever marginal planes (L/R) the root permits
 *           (r_allow_L/r_allow_R), then trace back with FULL marginal-mode
 *           tracking (mirrors cm_tr_alignT_hb): each node carries its on-path
 *           mode, B states consult the current mode's kshadow + Lkmode/Rkmode
 *           child-mode shadows, marginal states advance only the emitting end,
 *           and USED_TRUNC_END / USED_TRUNC_BEGIN terminate / root the marginal.
 *           For a single preset root mode (the byte-exact gate) exactly one of
 *           r_allow_{J,L,R} is set.
 */
static float
tr_insideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
	      int r, int z, int i0, int j0, int allow_begin,
	      int r_allow_J, int r_allow_L, int r_allow_R, int r_allow_T, CP9Bands_t *cp9b)
{
  int       status;
  void   ***shadow;
  TR_LR     lr;
  float     sc, retsc;
  ESL_STACK *pda_i, *pda_c;
  int       v,j,d,i;
  int       k;
  int       y, yoffset;
  int       bifparent;
  int       b, bb;
  float     bsc;
  char      mode, nxtmode, prvmode;
  int      *jmin = cp9b->jmin, *jmax = cp9b->jmax;
  int     **hdmin = cp9b->hdmin, **hdmax = cp9b->hdmax;

  /* T (brief 26_0610-051) needs the L+R child planes/shadows for the BEGL(R)+BEGR(L)
   * subtree traceback; fill them whenever r_allow_T (tr_inside_hb forces them too). */
  lr.fill_L = r_allow_L || r_allow_T; lr.fill_R = r_allow_R || r_allow_T; lr.fill_T = r_allow_T;
  lr.ret_planes = FALSE;
  lr.Lalpha = lr.Ralpha = NULL; lr.Lshad = lr.Rshad = NULL; lr.Lkmode = lr.Rkmode = NULL;
  lr.Lb = lr.Rb = -1; lr.Lbsc = lr.Rbsc = IMPOSSIBLE;
  lr.Tfull = NULL; lr.Tfullk = NULL; lr.Tb = -1; lr.Tbsc = IMPOSSIBLE;

  sc = tr_inside_hb(cm, dsq, L, r, z, i0, j0, BE_EFFICIENT,
		    NULL, NULL, NULL, NULL, &shadow, allow_begin, &b, &bsc, &lr, cp9b);

  /* resolve the root mode and the begin-target state / returned score */
  if      (r_allow_L) { mode = TRMODE_L; retsc = lr.Lbsc; bb = lr.Lb; }
  else if (r_allow_R) { mode = TRMODE_R; retsc = lr.Rbsc; bb = lr.Rb; }
  else if (r_allow_T) { mode = TRMODE_T; retsc = lr.Tbsc; bb = lr.Tb; }
  else                { mode = TRMODE_J; retsc = sc;      bb = b;     }

  pda_i = esl_stack_ICreate();
  pda_c = esl_stack_CCreate();
  if(pda_i == NULL || pda_c == NULL) goto ERROR;
  v = r; j = j0; i = i0; d = j0-i0+1;

  while (1) {
    int jp_v = 0, dp_v = 0, oob = 0;
    if (cm->sttype[v] == EL_st) { oob = 0; }
    else {
      oob = (j < jmin[v] || j > jmax[v]);
      if (!oob) { jp_v = j - jmin[v]; oob = (d < hd_min(cp9b, v, jp_v) || d > hd_max(cp9b, v, jp_v)); if (!oob) dp_v = d - hd_min(cp9b, v, jp_v); }
    }

    if (cm->sttype[v] == B_st) {
      if      (mode == TRMODE_J) k = ((int**)shadow[v])[j][dp_v];
      else if (mode == TRMODE_L) k = ((int**)lr.Lshad[v])[j][dp_v];
      else if (mode == TRMODE_R) k = ((int**)lr.Rshad[v])[j][dp_v];
      else                       k = lr.Tfullk[v];          /* T: full-span k* (brief 26_0610-051) */
      prvmode = mode;
      if      (mode == TRMODE_J) nxtmode = TRMODE_J;
      else if (mode == TRMODE_L) nxtmode = TRMODE_L;       /* L: right child stays L */
      else if (mode == TRMODE_R) nxtmode = lr.Rkmode[v][j][dp_v]; /* R: right child J or R */
      else                       nxtmode = TRMODE_L;       /* T: right child always L (oracle :615) */
      if((status = esl_stack_CPush(pda_c, nxtmode))    != eslOK) goto ERROR;
      if((status = esl_stack_IPush(pda_i, j))          != eslOK) goto ERROR;
      if((status = esl_stack_IPush(pda_i, k))          != eslOK) goto ERROR;
      if((status = esl_stack_IPush(pda_i, tr->n-1))    != eslOK) goto ERROR;
      if      (prvmode == TRMODE_J) mode = TRMODE_J;
      else if (prvmode == TRMODE_L) mode = lr.Lkmode[v][j][dp_v]; /* L: left child J or L */
      else if (prvmode == TRMODE_R) mode = TRMODE_R;             /* R: left child stays R */
      else                          mode = TRMODE_R;             /* T: left child always R (oracle :628) */
      j = j-k; d = d-k; i = j-d+1;
      y = cm->cfirst[v];
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y, mode);
      v = y;
    } else if (cm->sttype[v] == E_st || cm->sttype[v] == EL_st) {
      if (esl_stack_IPop(pda_i, &bifparent) == eslEOD) break;
      esl_stack_IPop(pda_i, &d);
      esl_stack_IPop(pda_i, &j);
      esl_stack_CPop(pda_c, &mode);
      v = tr->state[bifparent];
      y = cm->cnum[v];
      i = j-d+1;
      InsertTraceNodewithMode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y, mode);
      v = y;
    } else {
      int allow_S_trunc_end = 0;
      /* super-special: a BEGL_S/BEGR_S at d==0 in a disallowed/out-of-band mode is
       * a truncated end (the sister child emitted the whole subsequence). */
      if ((cm->stid[v]==BEGL_S || cm->stid[v]==BEGR_S) && d==0 &&
	  ((mode==TRMODE_J && !cp9b->Jvalid[v]) || (mode==TRMODE_L && !cp9b->Lvalid[v]) ||
	   (mode==TRMODE_R && !cp9b->Rvalid[v]) || oob))
	allow_S_trunc_end = 1;

      if      (allow_S_trunc_end)  yoffset = USED_TRUNC_END;
      else if (mode == TRMODE_J)   yoffset = ((char**)shadow[v])[j][dp_v];
      else if (mode == TRMODE_L)   yoffset = ((char**)lr.Lshad[v])[j][dp_v];
      else if (mode == TRMODE_R)   yoffset = ((char**)lr.Rshad[v])[j][dp_v];
      else                         yoffset = USED_TRUNC_BEGIN;  /* T at root v==0: begin into the B (brief 26_0610-051; oracle :683-685) */

      if      (yoffset == USED_TRUNC_BEGIN) { nxtmode = mode; }
      else if (yoffset == USED_TRUNC_END)   { }
      else if (yoffset == USED_EL)          { }
      else if (yoffset >= TRMODE_R_OFFSET)  { nxtmode = TRMODE_R; yoffset -= TRMODE_R_OFFSET; }
      else if (yoffset >= TRMODE_L_OFFSET)  { nxtmode = TRMODE_L; yoffset -= TRMODE_L_OFFSET; }
      else                                  { nxtmode = TRMODE_J; yoffset -= TRMODE_J_OFFSET; }

      switch (cm->sttype[v]) {
      case D_st: break;
      case MP_st: if (mode==TRMODE_J) { i++; j--; } else if (mode==TRMODE_L && d>0) i++; else if (mode==TRMODE_R && d>0) j--; break;
      case ML_st: if (mode==TRMODE_J || (mode==TRMODE_L && d>0)) i++; break;
      case MR_st: if (mode==TRMODE_J || (mode==TRMODE_R && d>0)) j--; break;
      case IL_st: if (mode==TRMODE_J || (mode==TRMODE_L && d>0)) i++; break;
      case IR_st: if (mode==TRMODE_J || (mode==TRMODE_R && d>0)) j--; break;
      case S_st:  break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
      }
      d = j-i+1;

      if (yoffset == USED_EL || yoffset == USED_TRUNC_END) {
	if (yoffset == USED_EL) InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M, mode);
	v = cm->M;
      } else if (yoffset == USED_TRUNC_BEGIN) {
	InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, bb, mode);
	v = bb;
      } else {
	mode = nxtmode;
	y = cm->cfirst[v] + yoffset;
	InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y, mode);
	v = y;
      }
    }
  }
  esl_stack_Destroy(pda_i);
  esl_stack_Destroy(pda_c);
  tr_free_lr_shadow(lr.Lshad, lr.Lkmode, cm, i0, j0);
  tr_free_lr_shadow(lr.Rshad, lr.Rkmode, cm, i0, j0);
  if (lr.Tfullk != NULL) free(lr.Tfullk);   /* T full-span k* array (brief 26_0610-051) */
  free_vjd_shadow_matrix(shadow, cm, i0, j0);
  CYKShadowTrackLeafDone();   /* brief 26_0610-050: this leaf's shadow is now freed */
  return retsc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.;
}

/* Function: tr_vinside_hb()  [brief 26_0610-044 R4.4a (J); brief 26_0610-046 R4.4b-pt2 (L/R)]
 *
 * J-plane unchanged from R4.4a. When vlr != NULL the L/R marginal vji planes are
 * also computed (mirrors truncyk.c::tr_vinside, in the banded vji idiom): La/Ra
 * score planes, Lsh/Rsh yshadows, and Lmode/Rmode marginal->child-mode shadows.
 * The boundary seed at z is gated by z_allow_{J,L,R}; the recurrences fill the
 * planes selected by vlr->fill_{L,R}. T is OFF (R4.4c).
 */
static float
tr_vinside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
	      int r, int z, int i0, int i1, int j1, int j0, int useEL,
	      int do_full, float ***a, float ****ret_a, char ****ret_shadow,
	      int allow_begin, int *ret_b, float *ret_bsc,
	      int z_allow_J, int z_allow_L, int z_allow_R, TR_VLR *vlr, CP9Bands_t *cp9b)
{
  int      status;
  char  ***shadow = NULL;
  int      v, i, j;
  int      w1, w2;
  int      jp, op, op_y;
  int     *touch;
  int      y, yoffset;
  float    sc, val;
  int      b;
  float    bsc;
  int      jpb, ilo, ihi;
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  /* L/R marginal vji state (brief 26_0610-046). All off when vlr==NULL. */
  int      fill_L = (vlr != NULL) ? vlr->fill_L : FALSE;
  int      fill_R = (vlr != NULL) ? vlr->fill_R : FALSE;
  float ***La = NULL, ***Ra = NULL;
  char  ***Lsh = NULL, ***Rsh = NULL, ***Lmode = NULL, ***Rmode = NULL;
  int      Lb = -1, Rb = -1;
  float    Lbsc = IMPOSSIBLE, Rbsc = IMPOSSIBLE;

  b   = -1;
  bsc = IMPOSSIBLE;
  if (cyk_dnc_track) cyk_dnc_vji_row_floats = i1 - i0 + 1;

  if (a == NULL) {
    ESL_ALLOC(a, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) a[v] = NULL;
  }
  if (fill_L) { ESL_ALLOC(La, sizeof(float **) * (cm->M+1)); for (v = 0; v <= cm->M; v++) La[v] = NULL; }
  if (fill_R) { ESL_ALLOC(Ra, sizeof(float **) * (cm->M+1)); for (v = 0; v <= cm->M; v++) Ra[v] = NULL; }

  w1 = cm->nodemap[cm->ndidx[z]];
  /* TRUNCATED marginal termini / begins (unlike classical CMH_LOCAL_BEGIN) may
   * land directly on an insert state -- a fragment can begin/end mid-insertion
   * (brief 26_0610-069 fix). Normally z is the split set's own last state (cfirst[w1]-1)
   * and w1..w2 is exactly that split set; the main recursion below (v = w1-1
   * downto r) never expands w1..w2 itself; it only reads them as already-seeded
   * boundary values. When z is instead a same-node insert state (index beyond
   * the split set, since inserts are numbered as split states' first children),
   * z needs the SAME treatment: seeded, not recursively expanded (its own
   * children lie beyond z, outside this V-problem's domain and never allocated
   * here). So the seeded range must extend through z, not stop at cfirst[w1]-1. */
  w2 = ESL_MAX(cm->cfirst[w1]-1, z);
  for (v = w1; v <= w2; v++) {
    a[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b);
    banded_hb_vji_init_impossible(a[v], i0, i1, j1, j0, v, cp9b);
    if (fill_L) { La[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b); banded_hb_vji_init_impossible(La[v], i0, i1, j1, j0, v, cp9b); }
    if (fill_R) { Ra[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b); banded_hb_vji_init_impossible(Ra[v], i0, i1, j1, j0, v, cp9b); }
  }

  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(char **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL;
    if (fill_L) { ESL_ALLOC(Lsh, sizeof(char **)*cm->M); ESL_ALLOC(Lmode, sizeof(char **)*cm->M);
                  for (v=0;v<cm->M;v++){Lsh[v]=NULL;Lmode[v]=NULL;} }
    if (fill_R) { ESL_ALLOC(Rsh, sizeof(char **)*cm->M); ESL_ALLOC(Rmode, sizeof(char **)*cm->M);
                  for (v=0;v<cm->M;v++){Rsh[v]=NULL;Rmode[v]=NULL;} }
  }

  if (! useEL) {
    if (z_allow_J && vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) a[z][0][op] = 0.;
    if (fill_L && z_allow_L && vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) La[z][0][op] = 0.;
    if (fill_R && z_allow_R && vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) Ra[z][0][op] = 0.;
  } else {
    if (ret_shadow != NULL) shadow[z] = alloc_banded_hb_vji_shadow_deck(i0,i1,j1,j0,z,cp9b);
    switch (cm->sttype[z]) {
    case D_st:
    case S_st:
      if (vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) {
	a[z][0][op] = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (ret_shadow != NULL) shadow[z][0][op] = USED_EL;
      }
      break;
    case MP_st:
      if (i0 == i1 || j1 == j0) break;
      if (vji_inband(cp9b, z, j1+1, i1-1, i0,i1,j1,j0, &op)) {
	val = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (dsq[i1-1] < cm->abc->K && dsq[j1+1] < cm->abc->K)
	  val += cm->esc[z][(int) (dsq[i1-1]*cm->abc->K+dsq[j1+1])];
	else
	  val += DegeneratePairScore(cm->abc, cm->esc[z], dsq[i1-1], dsq[j1+1]);
	if (val < IMPOSSIBLE) val = IMPOSSIBLE;
	a[z][1][op] = val;
	if (ret_shadow != NULL) shadow[z][1][op] = USED_EL;
      }
      break;
    case ML_st:
    case IL_st:
      if (i0 == i1) break;
      if (vji_inband(cp9b, z, j1, i1-1, i0,i1,j1,j0, &op)) {
	val = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (dsq[i1-1] < cm->abc->K)
	  val += cm->esc[z][(int) dsq[i1-1]];
	else
	  val += esl_abc_FAvgScore(cm->abc, dsq[i1-1], cm->esc[z]);
	if (val < IMPOSSIBLE) val = IMPOSSIBLE;
	a[z][0][op] = val;
	if (ret_shadow != NULL) shadow[z][0][op] = USED_EL;
      }
      break;
    case MR_st:
    case IR_st:
      if (j1 == j0) break;
      if (vji_inband(cp9b, z, j1+1, i1, i0,i1,j1,j0, &op)) {
	val = cm->endsc[z] + (cm->el_selfsc * ((j1)-(i1)+1));
	if (dsq[j1+1] < cm->abc->K)
	  val += cm->esc[z][(int) dsq[j1+1]];
	else
	  val += esl_abc_FAvgScore(cm->abc, dsq[j1+1], cm->esc[z]);
	if (val < IMPOSSIBLE) val = IMPOSSIBLE;
	a[z][1][op] = val;
	if (ret_shadow != NULL) shadow[z][1][op] = USED_EL;
      }
      break;
    }
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;    v < r;      v++) touch[v] = 0;
  for (v = r;    v <= w2;    v++) touch[v] = cm->pnum[v];
  for (v = w2+1; v < cm->M;  v++) touch[v] = 0;

  /* TRUNCATED begin straight into z on empty subsequences. */
  if (allow_begin && j0-j1 == 0 && i1-i0 == 0 && z != 0 && cp9b->Jvalid[z]) {
    float trpen = tr_trpenalty(cm, z);
    if (NOT_IMPOSSIBLE(trpen) && vji_inband(cp9b, z, j1, i1, i0,i1,j1,j0, &op)) {
      b   = z;
      bsc = a[z][0][op] + trpen;
      if (z == 0) {
	a[0][0][op] = bsc;
	if (ret_shadow != NULL) shadow[0][0][op] = USED_TRUNC_BEGIN;
      }
    }
  }

  for (v = w1-1; v >= r; v--)
    {
      a[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b);
      banded_hb_vji_init_impossible(a[v], i0, i1, j1, j0, v, cp9b);
      if (fill_L) { La[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b); banded_hb_vji_init_impossible(La[v], i0, i1, j1, j0, v, cp9b); }
      if (fill_R) { Ra[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b); banded_hb_vji_init_impossible(Ra[v], i0, i1, j1, j0, v, cp9b); }
      if (ret_shadow != NULL) {
	shadow[v] = alloc_banded_hb_vji_shadow_deck(i0, i1, j1, j0, v, cp9b);
	if (fill_L) { Lsh[v] = alloc_banded_hb_vji_shadow_deck(i0,i1,j1,j0,v,cp9b); Lmode[v] = alloc_banded_hb_vji_shadow_deck(i0,i1,j1,j0,v,cp9b); }
	if (fill_R) { Rsh[v] = alloc_banded_hb_vji_shadow_deck(i0,i1,j1,j0,v,cp9b); Rmode[v] = alloc_banded_hb_vji_shadow_deck(i0,i1,j1,j0,v,cp9b); }
      }

      if (cm->sttype[v] == E_st || cm->sttype[v] == B_st || (cm->sttype[v] == S_st && v > r))
	cm_Fail("you told me you wouldn't ever do that again.");

      /* TRUNCATED: ROOT_S (state 0 == r) gets no normal transitions; its a[] stays
       * IMPOSSIBLE except for the truncated begin (handled in begin bookkeeping). */
      if (allow_begin && v == 0) {
	/* nothing */
      }
      else
      for (j = ESL_MAX(j1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++) {
	jp  = j - j1;
	jpb = j - jmin[v];
	ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
	ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
	for (i = ihi; i >= ilo; i--)
	  {
	    int d = j - i + 1;
	    op = i - ilo;
	    y  = cm->cfirst[v];

	    /* brief 26_0610-086: mask the pure-J plane for a state v whose node
	     * structurally cannot appear in a pure-J parse (cp9b->Jvalid[v]==FALSE).
	     * The oracle (cm_TrCYKInsideAlignHB) never allocates a J deck for such a
	     * state, so it can be neither a J entry nor a J pass-through; here a[v] is
	     * read directly by v's parent's J recurrence (ungated), so an unmasked a[v]
	     * would let a truncated-begin J parse route THROUGH the non-J-valid state,
	     * reconstructing the same unreachable good-emission wedge the tr_outside_hb
	     * fix eliminates from the SCORE (this is its traceback twin). Mirrors the
	     * brief 26_0610-075 Jvalid[yy2] gate already on the L/R marginal cross-terms
	     * below, which the pure-J plane was missing. */
	    if (! cp9b->Jvalid[v]) {
	      a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	    }
	    else if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) {
	      if (vji_inband(cp9b, y, j, i, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j, i, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }
	    else if (cm->sttype[v] == MP_st) {
	      if (vji_inband(cp9b, y, j-1, i+1, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp-1][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j-1, i+1, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp-1][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		a[v][jp][op] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
	      else
		a[v][jp][op] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }
	    else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
	      if (vji_inband(cp9b, y, j, i+1, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j, i+1, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (dsq[i] < cm->abc->K)
		a[v][jp][op] += cm->esc[v][dsq[i]];
	      else
		a[v][jp][op] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }
	    else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
	      if (vji_inband(cp9b, y, j-1, i, i0,i1,j1,j0, &op_y))
		a[v][jp][op] = a[y][jp-1][op_y] + cm->tsc[v][0];
	      else a[v][jp][op] = IMPOSSIBLE;
	      if (ret_shadow != NULL) shadow[v][jp][op] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) &&
		  ((cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])))) > a[v][jp][op])) {
		a[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][op] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++)
		if (vji_inband(cp9b, y+yoffset, j-1, i, i0,i1,j1,j0, &op_y) &&
		    (sc = a[y+yoffset][jp-1][op_y] + cm->tsc[v][yoffset]) > a[v][jp][op]) {
		  a[v][jp][op] = sc;
		  if (ret_shadow != NULL) shadow[v][jp][op] = (char) yoffset;
		}
	      if (dsq[j] < cm->abc->K)
		a[v][jp][op] += cm->esc[v][dsq[j]];
	      else
		a[v][jp][op] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	      if (a[v][jp][op] < IMPOSSIBLE) a[v][jp][op] = IMPOSSIBLE;
	    }

	    /* ---- L/R marginal vji recurrences (brief 26_0610-047): a cell-for-cell
	     * translation of the vjd tr_inside_hb L/R recurrences into vji
	     * coords (d = j-i+1). MP uses lmesc/rmesc (NOT pair-esc); ML/IL & MR/IR
	     * try BOTH a J child (the J->marginal heal) and the same-mode child;
	     * D/S pass the mode through. Lmode/Rmode record the child mode for the
	     * mode-tracking traceback. ---- */
	    {
	      int   styp = cm->sttype[v];
	      int   sdl  = StateLeftDelta(styp);
	      int   sdr  = StateRightDelta(styp);
	      int   op_y;
	      float esc_i, esc_j;
	      esc_i = (dsq[i] < cm->abc->K) ? cm->esc[v][dsq[i]] : esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	      esc_j = (dsq[j] < cm->abc->K) ? cm->esc[v][dsq[j]] : esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);

	      /* ----- L marginal ----- */
	      if (fill_L && cp9b->Lvalid[v]) {
		if (styp == D_st || styp == S_st) {
		  if (d == 0) {
		    La[v][jp][op] = IMPOSSIBLE;
		    if (ret_shadow != NULL && styp == S_st) { Lsh[v][jp][op] = USED_TRUNC_END; Lmode[v][jp][op] = TRMODE_L; }
		    else if (ret_shadow != NULL) { Lsh[v][jp][op] = (char) 0; Lmode[v][jp][op] = TRMODE_L; }
		  } else {
		    La[v][jp][op] = IMPOSSIBLE;
		    if (ret_shadow != NULL) { Lsh[v][jp][op] = (char) 0; Lmode[v][jp][op] = TRMODE_L; }
		    /* brief 26_0610-096: marginal EL for a D/S state (d>0), L mirror of the R
		     * fix below -- see the full rationale there. sdl==0 for D/S. */
		    if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
		      La[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - sdl));
		      if (ret_shadow != NULL) Lsh[v][jp][op] = USED_EL;
		    }
		    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		      int yy2 = cm->cfirst[v] + yoffset;
		      if (cp9b->Lvalid[yy2] && vji_inband(cp9b, yy2, j, i, i0,i1,j1,j0, &op_y) &&
			  (sc = La[yy2][jp][op_y] + cm->tsc[v][yoffset]) > La[v][jp][op]) {
			La[v][jp][op] = sc; if (ret_shadow != NULL) { Lsh[v][jp][op] = (char) yoffset; Lmode[v][jp][op] = TRMODE_L; }
		      }
		    }
		    if (La[v][jp][op] < IMPOSSIBLE) La[v][jp][op] = IMPOSSIBLE;
		  }
		}
		else { /* emitting states: MP / ML / IL / MR / IR */
		  La[v][jp][op] = IMPOSSIBLE;
		  if (ret_shadow != NULL) { Lsh[v][jp][op] = USED_EL; Lmode[v][jp][op] = TRMODE_J; }
		  if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
		    La[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - sdl));
		    if (ret_shadow != NULL) Lsh[v][jp][op] = USED_EL;
		  }
		  if (styp == MP_st || styp == ML_st || styp == IL_st) {
		    /* L emits the left residue -> child drops it (j, i+1) */
		    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		      int yy2 = cm->cfirst[v] + yoffset;
		      if (vji_inband(cp9b, yy2, j, i+1, i0,i1,j1,j0, &op_y)) {
			/* brief 26_0610-092: only an MP L-plane transits to the child's J-plane;
			 * a plain LEFT-emitter (ML/IL) in L mode transits ONLY to the child's
			 * L-plane (oracle cm_TrCYKInsideAlignHB ML/IL L-recursion,
			 * cm_dpalign_trunc.c:7050-7054). The prior unconditional J-child term was
			 * the tr_vinside_hb copy of the inside-begin inflation this brief fixes in
			 * tr_inside_hb -- it is what made the byte-exact-score traceback still
			 * reconstruct an ML->EL phantom. (brief 26_0610-075's Jvalid[yy2] gate stays
			 * for the MP case.) */
			if (styp == MP_st && cp9b->Jvalid[yy2] && (sc = a[yy2][jp][op_y] + cm->tsc[v][yoffset]) > La[v][jp][op]) {
			  La[v][jp][op] = sc; if (ret_shadow != NULL) { Lsh[v][jp][op] = (char) yoffset; Lmode[v][jp][op] = TRMODE_J; }
			}
			if (cp9b->Lvalid[yy2] && (sc = La[yy2][jp][op_y] + cm->tsc[v][yoffset]) > La[v][jp][op]) {
			  La[v][jp][op] = sc; if (ret_shadow != NULL) { Lsh[v][jp][op] = (char) yoffset; Lmode[v][jp][op] = TRMODE_L; }
			}
		      }
		    }
		    if (d >= 2) La[v][jp][op] += (styp == MP_st) ? cm->lmesc[v][dsq[i]] : esc_i;
		    else { La[v][jp][op] = (styp == MP_st) ? cm->lmesc[v][dsq[i]] : esc_i;
			   if (ret_shadow != NULL) Lsh[v][jp][op] = USED_TRUNC_END; }
		  }
		  else { /* MR / IR : L mode emits nothing here -> child keeps d (j, i) */
		    int Lyoffset0 = (styp == IR_st) ? 1 : 0;
		    for (yoffset = Lyoffset0; yoffset < cm->cnum[v]; yoffset++) {
		      int yy2 = cm->cfirst[v] + yoffset;
		      if (vji_inband(cp9b, yy2, j, i, i0,i1,j1,j0, &op_y)) {
			/* brief 26_0610-075: same Jvalid[yy2] gate as the MP/ML/IL branch above. */
			if (cp9b->Jvalid[yy2] && (sc = a[yy2][jp][op_y] + cm->tsc[v][yoffset]) > La[v][jp][op]) {
			  La[v][jp][op] = sc; if (ret_shadow != NULL) { Lsh[v][jp][op] = (char) yoffset; Lmode[v][jp][op] = TRMODE_J; }
			}
			if (cp9b->Lvalid[yy2] && (sc = La[yy2][jp][op_y] + cm->tsc[v][yoffset]) > La[v][jp][op]) {
			  La[v][jp][op] = sc; if (ret_shadow != NULL) { Lsh[v][jp][op] = (char) yoffset; Lmode[v][jp][op] = TRMODE_L; }
			}
		      }
		    }
		  }
		  if (La[v][jp][op] < IMPOSSIBLE) La[v][jp][op] = IMPOSSIBLE;
		}
	      }

	      /* ----- R marginal ----- */
	      if (fill_R && cp9b->Rvalid[v]) {
		if (styp == D_st || styp == S_st) {
		  if (d == 0) {
		    Ra[v][jp][op] = IMPOSSIBLE;
		    if (ret_shadow != NULL && styp == S_st) { Rsh[v][jp][op] = USED_TRUNC_END; Rmode[v][jp][op] = TRMODE_R; }
		    else if (ret_shadow != NULL) { Rsh[v][jp][op] = (char) 0; Rmode[v][jp][op] = TRMODE_R; }
		  } else {
		    Ra[v][jp][op] = IMPOSSIBLE;
		    if (ret_shadow != NULL) { Rsh[v][jp][op] = (char) 0; Rmode[v][jp][op] = TRMODE_R; }
		    /* brief 26_0610-096: marginal EL for a D/S state (d>0). The oracle's EL
		     * reinit (cm_TrCYKInsideAlignHB ~7011-7069) is state-type-INDEPENDENT --
		     * every state with a valid endsc[v] seeds its L/R (and J) decks with the
		     * local-end score before the recurrence. This D/S branch was the one place
		     * that omitted it (the emitting-state branch below already does it), so an
		     * R-mode BEGL_S->EL terminus was scorable in the outside (betaR[cm->M]) but
		     * NOT reconstructable in the inside -> the CsrB-sample13_5ptr NULL-deref.
		     * sdr==0 for D/S, so this matches Ralpha[cm->M][j][d]+endsc[v] = el*d+endsc. */
		    if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
		      Ra[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - sdr));
		      if (ret_shadow != NULL) Rsh[v][jp][op] = USED_EL;
		    }
		    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		      int yy2 = cm->cfirst[v] + yoffset;
		      if (cp9b->Rvalid[yy2] && vji_inband(cp9b, yy2, j, i, i0,i1,j1,j0, &op_y) &&
			  (sc = Ra[yy2][jp][op_y] + cm->tsc[v][yoffset]) > Ra[v][jp][op]) {
			Ra[v][jp][op] = sc; if (ret_shadow != NULL) { Rsh[v][jp][op] = (char) yoffset; Rmode[v][jp][op] = TRMODE_R; }
		      }
		    }
		    if (Ra[v][jp][op] < IMPOSSIBLE) Ra[v][jp][op] = IMPOSSIBLE;
		  }
		}
		else { /* emitting states */
		  Ra[v][jp][op] = IMPOSSIBLE;
		  if (ret_shadow != NULL) { Rsh[v][jp][op] = USED_EL; Rmode[v][jp][op] = TRMODE_J; }
		  if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
		    Ra[v][jp][op] = cm->endsc[v] + (cm->el_selfsc * (d - sdr));
		    if (ret_shadow != NULL) Rsh[v][jp][op] = USED_EL;
		  }
		  if (styp == MP_st || styp == MR_st || styp == IR_st) {
		    /* R emits the right residue -> child drops it (j-1, i) */
		    for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
		      int yy2 = cm->cfirst[v] + yoffset;
		      if (vji_inband(cp9b, yy2, j-1, i, i0,i1,j1,j0, &op_y)) {
			/* brief 26_0610-092: R-mode mirror -- only an MP R-plane transits to the
			 * child's J-plane; a plain RIGHT-emitter (MR/IR) in R mode transits ONLY to
			 * the child's R-plane (oracle IR/MR R-recursion, cm_dpalign_trunc.c:7158-7162). */
			if (styp == MP_st && cp9b->Jvalid[yy2] && (sc = a[yy2][jp-1][op_y] + cm->tsc[v][yoffset]) > Ra[v][jp][op]) {
			  Ra[v][jp][op] = sc; if (ret_shadow != NULL) { Rsh[v][jp][op] = (char) yoffset; Rmode[v][jp][op] = TRMODE_J; }
			}
			if (cp9b->Rvalid[yy2] && (sc = Ra[yy2][jp-1][op_y] + cm->tsc[v][yoffset]) > Ra[v][jp][op]) {
			  Ra[v][jp][op] = sc; if (ret_shadow != NULL) { Rsh[v][jp][op] = (char) yoffset; Rmode[v][jp][op] = TRMODE_R; }
			}
		      }
		    }
		    if (d >= 2) Ra[v][jp][op] += (styp == MP_st) ? cm->rmesc[v][dsq[j]] : esc_j;
		    else { Ra[v][jp][op] = (styp == MP_st) ? cm->rmesc[v][dsq[j]] : esc_j;
			   if (ret_shadow != NULL) Rsh[v][jp][op] = USED_TRUNC_END; }
		  }
		  else { /* ML / IL : R mode emits nothing here -> child keeps d (j, i) */
		    int Ryoffset0 = (styp == IL_st) ? 1 : 0;
		    for (yoffset = Ryoffset0; yoffset < cm->cnum[v]; yoffset++) {
		      int yy2 = cm->cfirst[v] + yoffset;
		      if (vji_inband(cp9b, yy2, j, i, i0,i1,j1,j0, &op_y)) {
			/* brief 26_0610-075: same Jvalid[yy2] gate as the MP/MR/IR branch above. */
			if (cp9b->Jvalid[yy2] && (sc = a[yy2][jp][op_y] + cm->tsc[v][yoffset]) > Ra[v][jp][op]) {
			  Ra[v][jp][op] = sc; if (ret_shadow != NULL) { Rsh[v][jp][op] = (char) yoffset; Rmode[v][jp][op] = TRMODE_J; }
			}
			if (cp9b->Rvalid[yy2] && (sc = Ra[yy2][jp][op_y] + cm->tsc[v][yoffset]) > Ra[v][jp][op]) {
			  Ra[v][jp][op] = sc; if (ret_shadow != NULL) { Rsh[v][jp][op] = (char) yoffset; Rmode[v][jp][op] = TRMODE_R; }
			}
		      }
		    }
		  }
		  if (Ra[v][jp][op] < IMPOSSIBLE) Ra[v][jp][op] = IMPOSSIBLE;
		}
	      }
	    }
	  }
      }

      /* TRUNCATED begin bookkeeping (cell (j0,i0)) */
      if (allow_begin && v != 0 && cp9b->Jvalid[v] && vji_inband(cp9b, v, j0, i0, i0,i1,j1,j0, &op)) {
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) && a[v][j0-j1][op] + trpen > bsc) {
	  b   = v;
	  bsc = a[v][j0-j1][op] + trpen;
	}
      }
      if (allow_begin && v == 0 && cp9b->Jvalid[0] && vji_inband(cp9b, 0, j0, i0, i0,i1,j1,j0, &op)) {
	a[0][j0-j1][op] = bsc;
	if (ret_shadow != NULL) shadow[v][j0-j1][op] = USED_TRUNC_BEGIN;
      }
      /* marginal truncated-begin bookkeeping (cell (j0,i0)), per mode */
      if (allow_begin && v != 0 && fill_L && cp9b->Lvalid[v] && vji_inband(cp9b, v, j0, i0, i0,i1,j1,j0, &op)) {
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) && La[v][j0-j1][op] + trpen > Lbsc) { Lb = v; Lbsc = La[v][j0-j1][op] + trpen; }
      }
      if (allow_begin && v == 0 && fill_L && cp9b->Lvalid[0] && vji_inband(cp9b, 0, j0, i0, i0,i1,j1,j0, &op)) {
	La[0][j0-j1][op] = Lbsc;
	if (ret_shadow != NULL) Lsh[v][j0-j1][op] = USED_TRUNC_BEGIN;
      }
      if (allow_begin && v != 0 && fill_R && cp9b->Rvalid[v] && vji_inband(cp9b, v, j0, i0, i0,i1,j1,j0, &op)) {
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) && Ra[v][j0-j1][op] + trpen > Rbsc) { Rb = v; Rbsc = Ra[v][j0-j1][op] + trpen; }
      }
      if (allow_begin && v == 0 && fill_R && cp9b->Rvalid[0] && vji_inband(cp9b, 0, j0, i0, i0,i1,j1,j0, &op)) {
	Ra[0][j0-j1][op] = Rbsc;
	if (ret_shadow != NULL) Rsh[v][j0-j1][op] = USED_TRUNC_BEGIN;
      }

      if (! do_full) {
	for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++) {
	  touch[y]--;
	  if (touch[y] == 0) { free_banded_hb_vji_deck(a[y], i0,i1,j1,j0, y, cp9b); a[y] = NULL;
	    if (fill_L) { free_banded_hb_vji_deck(La[y], i0,i1,j1,j0, y, cp9b); La[y] = NULL; }
	    if (fill_R) { free_banded_hb_vji_deck(Ra[y], i0,i1,j1,j0, y, cp9b); Ra[y] = NULL; }
	  }
	}
      }
    }

  { int dpr; sc = vji_inband(cp9b, r, j0, i0, i0,i1,j1,j0, &dpr) ? a[r][j0-j1][dpr] : IMPOSSIBLE; }
  if (ret_b   != NULL) *ret_b   = b;
  if (ret_bsc != NULL) *ret_bsc = bsc;

  if (ret_a == NULL) {
    for (v = r; v <= w2; v++)
      if (a[v] != NULL) { free_banded_hb_vji_deck(a[v], i0,i1,j1,j0, v, cp9b); a[v] = NULL; }
    free(a);
  } else *ret_a = a;

  /* marginal L/R score planes are never returned (only the shadows feed the
   * mode-tracking traceback); free any decks still live. */
  if (fill_L) {
    for (v = r; v <= w2; v++) if (La[v] != NULL) { free_banded_hb_vji_deck(La[v], i0,i1,j1,j0, v, cp9b); La[v] = NULL; }
    free(La); La = NULL;
  }
  if (fill_R) {
    for (v = r; v <= w2; v++) if (Ra[v] != NULL) { free_banded_hb_vji_deck(Ra[v], i0,i1,j1,j0, v, cp9b); Ra[v] = NULL; }
    free(Ra); Ra = NULL;
  }

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  /* hand the marginal shadows (+ child-mode shadows) and per-mode begins to the
   * caller (tr_vinsideT_hb) via the vlr bundle. */
  if (vlr != NULL) {
    vlr->Lsh = Lsh; vlr->Rsh = Rsh; vlr->Lmode = Lmode; vlr->Rmode = Rmode;
    vlr->Lb = Lb; vlr->Lbsc = Lbsc; vlr->Rb = Rb; vlr->Rbsc = Rbsc;
  }
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return 0.;
}

/* Function: tr_voutside_hb()  [brief 26_0610-044, R4.4a] -- J-plane truncated analogue of voutside_hb(). */
static void
tr_voutside_hb(CM_t *cm, ESL_DSQ *dsq, int L,
	       int r, int z, int i0, int i1, int j1, int j0, int useEL,
	       int do_full, float ***beta, float ****ret_beta, CP9Bands_t *cp9b)
{
  int      status;
  int      v, y;
  int      i, j;
  int      jp, ip, op, op_y;
  float    sc;
  int     *touch;
  float    escore;
  int      voffset;
  int      jpb, ilo, ihi;
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;

  if (cyk_dnc_track) cyk_dnc_vji_row_floats = i1 - i0 + 1;

  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) beta[v] = NULL;
  }

  beta[r] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, r, cp9b);
  banded_hb_vji_init_impossible(beta[r], i0, i1, j1, j0, r, cp9b);
  /* TRUNCATED: omit normal root seed when r==0 (no normal ROOT_S descent). */
  if (r != 0 && vji_inband(cp9b, r, j0, i0, i0,i1,j1,j0, &op)) beta[r][j0-j1][op] = 0;

  if (useEL && cm->flags & CMH_LOCAL_END) {
    beta[cm->M] = alloc_vji_deck(i0, i1, j1, j0);
    for (jp = 0; jp <= j0-j1; jp++)
      for (ip = 0; ip <= i1-i0; ip++)
	beta[cm->M][jp][ip] = IMPOSSIBLE;
  }
  if (r != 0 && useEL && NOT_IMPOSSIBLE(cm->endsc[r])) {
    switch (cm->sttype[r]) {
    case MP_st:
      if (i0 == i1 || j1 == j0) break;
      if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	escore = cm->esc[r][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
      else
	escore = DegeneratePairScore(cm->abc, cm->esc[r], dsq[i0], dsq[j0]);
      beta[cm->M][j0-j1-1][1] = cm->endsc[r] + (cm->el_selfsc * ((j0-1)-(i0+1)+1)) + escore;
      break;
    case ML_st:
    case IL_st:
      if (i0 == i1) break;
      if (dsq[i0] < cm->abc->K) escore = cm->esc[r][(int) dsq[i0]];
      else                      escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[r]);
      beta[cm->M][j0-j1][1] = cm->endsc[r] + (cm->el_selfsc * ((j0)-(i0+1)+1)) + escore;
      break;
    case MR_st:
    case IR_st:
      if (j0 == j1) break;
      if (dsq[j0] < cm->abc->K) escore = cm->esc[r][(int) dsq[j0]];
      else                      escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[r]);
      beta[cm->M][j0-j1-1][0] = cm->endsc[r] + (cm->el_selfsc * ((j0-1)-(i0)+1)) + escore;
      break;
    case S_st:
    case D_st:
      beta[cm->M][j0-j1][0] = cm->endsc[r] + (cm->el_selfsc * ((j0)-(i0)+1));
      break;
    default:  cm_Fail("bogus parent state %d\n", cm->sttype[r]);
    }
  }

  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;   v < r;     v++) touch[v] = 0;
  for (v = z+1; v < cm->M; v++) touch[v] = 0;
  for (v = r;   v <= z;    v++) {
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
  }

  for (v = r+1; v <= z; v++)
    {
      beta[v] = alloc_banded_hb_vji_deck(i0, i1, j1, j0, v, cp9b);
      banded_hb_vji_init_impossible(beta[v], i0, i1, j1, j0, v, cp9b);

      /* TRUNCATED begin into v at (j0,i0) -- penalty-aware, all Jvalid states. */
      if (r == 0 && i0 == 1 && j0 == L && cp9b->Jvalid[v]) {
	float trpen = tr_trpenalty(cm, v);
	if (NOT_IMPOSSIBLE(trpen) && vji_inband(cp9b, v, j0, i0, i0,i1,j1,j0, &op) && trpen > beta[v][j0-j1][op])
	  beta[v][j0-j1][op] = trpen;
      }

      for (j = ESL_MIN(j0, jmax[v]); j >= ESL_MAX(j1, jmin[v]); j--) {
	jp  = j - j1;
	jpb = j - jmin[v];
	ilo = j - hd_max(cp9b, v, jpb) + 1;  if (ilo < i0) ilo = i0;
	ihi = j - hd_min(cp9b, v, jpb) + 1;  if (ihi > i1) ihi = i1;
	for (i = ilo; i <= ihi; i++) {
	  op = i - ilo;
	  for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	    if (y < r) continue;
	    voffset = v - cm->cfirst[y];
	    switch (cm->sttype[y]) {
	    case MP_st:
	      if (j == j0 || i == i0) continue;
	      if (! vji_inband(cp9b, y, j+1, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
	      else
		escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
	      if ((sc = beta[y][jp+1][op_y] + cm->tsc[y][voffset] + escore) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    case ML_st:
	    case IL_st:
	      if (i == i0) continue;
	      if (! vji_inband(cp9b, y, j, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K) escore = cm->esc[y][(int) dsq[i-1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
	      if ((sc = beta[y][jp][op_y] + cm->tsc[y][voffset] + escore) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    case MR_st:
	    case IR_st:
	      if (j == j0) continue;
	      if (! vji_inband(cp9b, y, j+1, i, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[j+1] < cm->abc->K) escore = cm->esc[y][(int) dsq[j+1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);
	      if ((sc = beta[y][jp+1][op_y] + cm->tsc[y][voffset] + escore) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    case S_st:
	    case E_st:
	    case D_st:
	      if (! vji_inband(cp9b, y, j, i, i0,i1,j1,j0, &op_y)) continue;
	      if ((sc = beta[y][jp][op_y] + cm->tsc[y][voffset]) > beta[v][jp][op])
		beta[v][jp][op] = sc;
	      break;
	    default: cm_Fail("bogus parent state %d\n", cm->sttype[y]);
	    }
	  }
	  if (beta[v][jp][op] < IMPOSSIBLE) beta[v][jp][op] = IMPOSSIBLE;
	}
      }

      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (jp = j0-j1; jp >= 0; jp--) {
	  j = jp + j1;
	  for (ip = 0; ip <= i1-i0; ip++) {
	    i = ip + i0;
	    switch (cm->sttype[v]) {
	    case MP_st:
	      if (j == j0 || i == i0) continue;
	      if (! vji_inband(cp9b, v, j+1, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
	      else
		escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
	      if ((sc = beta[v][jp+1][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    case ML_st:
	    case IL_st:
	      if (i == i0) continue;
	      if (! vji_inband(cp9b, v, j, i-1, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[i-1] < cm->abc->K) escore = cm->esc[v][(int) dsq[i-1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
	      if ((sc = beta[v][jp][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    case MR_st:
	    case IR_st:
	      if (j == j0) continue;
	      if (! vji_inband(cp9b, v, j+1, i, i0,i1,j1,j0, &op_y)) continue;
	      if (dsq[j+1] < cm->abc->K) escore = cm->esc[v][(int) dsq[j+1]];
	      else                       escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
	      if ((sc = beta[v][jp+1][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1)) + escore) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    case S_st:
	    case D_st:
	    case E_st:
	      if (! vji_inband(cp9b, v, j, i, i0,i1,j1,j0, &op_y)) continue;
	      if ((sc = beta[v][jp][op_y] + cm->endsc[v] + (cm->el_selfsc * (j-i+1))) > beta[cm->M][jp][ip])
		beta[cm->M][jp][ip] = sc;
	      break;
	    default:  cm_Fail("bogus parent state %d\n", cm->sttype[v]);
	    }
	  }
	}
      }

      if (! do_full) {
	for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	  touch[y]--;
	  if (touch[y] == 0) { free_banded_hb_vji_deck(beta[y], i0,i1,j1,j0, y, cp9b); beta[y] = NULL; }
	}
      }
    }

  if (ret_beta == NULL) {
    for (v = r; v <= z; v++)
      if (beta[v] != NULL) { free_banded_hb_vji_deck(beta[v], i0,i1,j1,j0, v, cp9b); beta[v] = NULL; }
    if (useEL && cm->flags & CMH_LOCAL_END) { free_vji_deck(beta[cm->M], j1, j0); beta[cm->M] = NULL; }
    free(beta);
  } else *ret_beta = beta;

  free(touch);
  return;

 ERROR:
  cm_Fail("Memory allocation error.\n");
}

/* Function: tr_vinsideT_hb()  [brief 26_0610-044 R4.4a (J); brief 26_0610-047 R4.4b-pt2b (L/R)]
 *
 * Whole-solve a (possibly marginal) V-problem and trace it back with full
 * marginal-mode tracking. J (r_allow_J & z_allow_J) is the R4.4a path; L/R route
 * the root mode through the now-implemented marginal vji recurrences. Each node
 * carries its on-path mode; the per-mode shadow (J=shadow, L=Lsh/Lmode,
 * R=Rsh/Rmode) supplies the child yoffset + child mode; marginal states advance
 * only the emitting end; USED_TRUNC_END terminates the marginal chain (the
 * emitting leaf), USED_EL/USED_TRUNC_BEGIN as in J.
 */
static float
tr_vinsideT_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
	       int r, int z, int i0, int i1, int j1, int j0, int useEL,
	       int allow_begin, int r_allow_J, int r_allow_L, int r_allow_R,
	       int z_allow_J, int z_allow_L, int z_allow_R, CP9Bands_t *cp9b)
{
  char ***shadow;
  TR_VLR  vlr;
  float   sc;
  int     v, y;
  int     j, i, d, op;
  int     yoffset;
  int     b, bb;
  float   bsc;
  char    mode, nxtmode;

  /* the V-problem's root mode (exactly one of r_allow_{J,L,R} is set) */
  mode = r_allow_L ? TRMODE_L : (r_allow_R ? TRMODE_R : TRMODE_J);

  if (r == z) {
    InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, r, mode);
    return 0.;
  }

  vlr.fill_L = r_allow_L; vlr.fill_R = r_allow_R;
  vlr.La = vlr.Ra = NULL; vlr.Lsh = vlr.Rsh = NULL; vlr.Lmode = vlr.Rmode = NULL;
  vlr.Lb = vlr.Rb = -1; vlr.Lbsc = vlr.Rbsc = IMPOSSIBLE;

  sc = tr_vinside_hb(cm, dsq, L, r, z, i0, i1, j1, j0, useEL,
		     BE_EFFICIENT, NULL, NULL, &shadow,
		     allow_begin, &b, &bsc,
		     z_allow_J, z_allow_L, z_allow_R, &vlr, cp9b);

  bb = (mode == TRMODE_L) ? vlr.Lb : (mode == TRMODE_R) ? vlr.Rb : b;

  v = r;
  j = j0;
  i = i0;
  while (1) {
    if (! vji_inband(cp9b, v, j, i, i0,i1,j1,j0, &op)) {
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, z, mode);
      break;
    }
    d = j - i + 1;
    if      (mode == TRMODE_J) { yoffset = shadow[v][j-j1][op];   nxtmode = TRMODE_J; }
    else if (mode == TRMODE_L) { yoffset = vlr.Lsh[v][j-j1][op];  nxtmode = vlr.Lmode[v][j-j1][op]; }
    else                       { yoffset = vlr.Rsh[v][j-j1][op];  nxtmode = vlr.Rmode[v][j-j1][op]; }

    /* advance coords past v's emission (mode-dependent; only the emitting end) */
    switch (cm->sttype[v]) {
    case D_st:  break;
    case S_st:  break;
    case MP_st: if (mode==TRMODE_J) { i++; j--; } else if (mode==TRMODE_L && d>0) i++; else if (mode==TRMODE_R && d>0) j--; break;
    case ML_st: if (mode==TRMODE_J || (mode==TRMODE_L && d>0)) i++; break;
    case IL_st: if (mode==TRMODE_J || (mode==TRMODE_L && d>0)) i++; break;
    case MR_st: if (mode==TRMODE_J || (mode==TRMODE_R && d>0)) j--; break;
    case IR_st: if (mode==TRMODE_J || (mode==TRMODE_R && d>0)) j--; break;
    default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
    }

    if (yoffset == USED_TRUNC_END) {
      break;                          /* marginal leaf: the chain ends here */
    }
    else if (yoffset == USED_EL) {
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, cm->M, mode);
      break;
    }
    else if (yoffset == USED_TRUNC_BEGIN) {
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, bb, mode);
      v = bb;
      if (! useEL && v == z) break;
    }
    else {
      mode = nxtmode;
      y = cm->cfirst[v] + yoffset;
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, y, mode);
      v = y;
      if (! useEL && v == z) break;
    }
  }

  free_vji_shadow_matrix(shadow, cm->M, j1, j0);
  if (r_allow_L) { free_vji_shadow_matrix(vlr.Lsh, cm->M, j1, j0); free_vji_shadow_matrix(vlr.Lmode, cm->M, j1, j0); }
  if (r_allow_R) { free_vji_shadow_matrix(vlr.Rsh, cm->M, j1, j0); free_vji_shadow_matrix(vlr.Rmode, cm->M, j1, j0); }
  CYKShadowTrackLeafDone();   /* brief 26_0610-050: this V-problem's shadow is now freed */
  return sc;
}

/* Function: tr_v_splitter_hb()  [brief 26_0610-044 R4.4a (J); brief 26_0610-046 R4.4b-pt2 (L/R)]
 *
 * J V-problems (r_allow_J && z_allow_J) keep the exact R4.4a D&C split. Marginal
 * or mixed-mode V-problems (any L/R in r_allow_* or z_allow_*) are solved WHOLE
 * via the mode-aware tr_vinsideT_hb (the marginal vji cube is tiny -- 009 banded
 * it to ~0 even on LSU -- so the whole-solve is memory-safe and avoids a separate
 * marginal voutside + marginal v-split). T is OFF (R4.4c).
 */
static void
tr_v_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
		 int r, int z, int i0, int i1, int j1, int j0, int useEL,
		 int r_allow_J, int r_allow_L, int r_allow_R,
		 int z_allow_J, int z_allow_L, int z_allow_R, CP9Bands_t *cp9b)
{
  float ***alpha, ***beta;
  float    sc;
  int      v, w, y;
  int      i, j, jp, ip, op;
  int      best_v, best_i, best_j;
  float    best_sc;
  int      midnode;
  int      b;
  float    bsc;

  /* Marginal / mixed-mode V-problem: solve whole (mode-aware traceback). */
  if (r_allow_L || r_allow_R || z_allow_L || z_allow_R) {
    tr_vinsideT_hb(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, (r==0),
		   r_allow_J, r_allow_L, r_allow_R, z_allow_J, z_allow_L, z_allow_R, cp9b);
    return;
  }

  /* ---- pure-J V-problem: the R4.4a memory-efficient D&C split (unchanged) ---- */
  if (cm->ndidx[z] == cm->ndidx[r] + 1 || r == z ||
      vinsideT_size(cm, r, z, i0, i1, j1, j0) < RAMLIMIT) {
    tr_vinsideT_hb(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, (r==0),
		   r_allow_J, r_allow_L, r_allow_R, z_allow_J, z_allow_L, z_allow_R, cp9b);
    return;
  }

  midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r]) / 2);
  w = cm->nodemap[midnode];
  y = cm->cfirst[w]-1;

  tr_vinside_hb (cm, dsq, L, w, z, i0, i1, j1, j0, useEL, BE_EFFICIENT,
		 NULL, &alpha, NULL, (r==0), &b, &bsc, TRUE, FALSE, FALSE, NULL, cp9b);
  tr_voutside_hb(cm, dsq, L, r, y, i0, i1, j1, j0, useEL, BE_EFFICIENT,
		 NULL, &beta, cp9b);

  best_sc = IMPOSSIBLE;
  best_v  = -99;
  for (v = w; v <= y; v++)
    for (ip = 0; ip <= i1-i0; ip++) {
      i = ip + i0;
      for (jp = 0; jp <= j0-j1; jp++) {
	j = jp + j1;
	if (vji_inband(cp9b, v, j, i, i0,i1,j1,j0, &op) &&
	    (sc = alpha[v][jp][op] + beta[v][jp][op]) > best_sc) {
	  best_sc = sc;
	  best_v  = v;
	  best_i  = i;
	  best_j  = j;
	}
      }
    }

  if (useEL && (cm->flags & CMH_LOCAL_END)) {
    for (ip = 0; ip <= i1-i0; ip++)
      for (jp = 0; jp <= j0-j1; jp++)
	if ((sc = beta[cm->M][jp][ip]) > best_sc) {
	  best_sc = sc;
	  best_v  = -1;
	  best_i  = ip + i0;
	  best_j  = jp + j1;
	}
  }

  if (r == 0) {
    if (bsc > best_sc) {
      best_sc = bsc;
      best_v  = -2;
      best_i  = i0;
      best_j  = j0;
    }
  }

  free_banded_hb_vji_matrix(alpha, cm, i0, i1, j1, j0, cp9b);
  free_banded_hb_vji_matrix(beta,  cm, i0, i1, j1, j0, cp9b);

  if (best_v == -99 || ! NOT_IMPOSSIBLE(best_sc)) {
    tr_vinsideT_hb(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, (r==0),
		   r_allow_J, r_allow_L, r_allow_R, z_allow_J, z_allow_L, z_allow_R, cp9b);
    return;
  }

  if (best_v == -1) {
    tr_v_splitter_hb(cm, dsq, L, tr, r, w, i0, best_i, best_j, j0, TRUE,
		     r_allow_J, r_allow_L, r_allow_R, TRUE, FALSE, FALSE, cp9b);
    return;
  }
  if (best_v == -2) {
    if (b != z)
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b, TRMODE_J);
    tr_v_splitter_hb(cm, dsq, L, tr, b, z, i0, i1, j1, j0, useEL,
		     r_allow_J, r_allow_L, r_allow_R, z_allow_J, z_allow_L, z_allow_R, cp9b);
    return;
  }

  tr_v_splitter_hb(cm, dsq, L, tr, r,      best_v, i0,     best_i, best_j, j0, FALSE,
		   r_allow_J, r_allow_L, r_allow_R, TRUE, FALSE, FALSE, cp9b);
  tr_v_splitter_hb(cm, dsq, L, tr, best_v, z,      best_i, i1,     j1,     best_j, useEL,
		   TRUE, FALSE, FALSE, z_allow_J, z_allow_L, z_allow_R, cp9b);
  return;
}

/* Function: tr_wedge_splitter_hb()  [brief 26_0610-044, R4.4a] -- J-plane truncated analogue of wedge_splitter_hb(). */
static float
tr_wedge_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int j0,
		     int r_allow_J, int r_allow_L, int r_allow_R, CP9Bands_t *cp9b)
{
  float ***alpha = NULL, ***Lalpha = NULL, ***Ralpha = NULL;
  float ***beta  = NULL;
  float ***betaL = NULL, ***betaR = NULL;
  int      fill_L = r_allow_L, fill_R = r_allow_R;
  float sc;
  float best_sc;
  int   v,w,y;
  int   W;
  int   d, jp, j, jp_v;
  int   best_v, best_d, best_j;
  int   p_mode, c_mode;          /* parent (V-problem) mode / child (wedge) mode */
  int   midnode;
  TR_LR lr;
  int   binJ;  float binJsc;     int binv, binmode;  float binsc;
  float boutsc; int boutv, boutj, boutmode;
  int  *jmin  = cp9b->jmin;
  int  *jmax  = cp9b->jmax;
  int **hdmin = cp9b->hdmin;
  int **hdmax = cp9b->hdmax;
  int  *eldmaxJ = NULL, *eldmaxL = NULL, *eldmaxR = NULL;  /* brief 26_0430-227: banded EL d-band edges (recomputed) */

  if (cm->ndidx[z] == cm->ndidx[r] + 1 ||
      insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT)
    {
      sc = tr_insideT_hb(cm, dsq, L, tr, r, z, i0, j0, (r==0), r_allow_J, r_allow_L, r_allow_R, FALSE, cp9b);
      return sc;
    }

  midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r]) / 2);
  w = cm->nodemap[midnode];
  y = cm->cfirst[w]-1;

  lr.fill_L = fill_L; lr.fill_R = fill_R; lr.fill_T = FALSE; lr.ret_planes = TRUE; lr.Lalpha = NULL; lr.Ralpha = NULL;  /* T never reaches the wedge (brief 26_0610-051) */
  tr_inside_hb(cm, dsq, L, w, z, i0, j0, BE_EFFICIENT,
	       NULL, &alpha, NULL, NULL, NULL,
	       (r==0), &binJ, &binJsc, &lr, cp9b);
  Lalpha = lr.Lalpha; Ralpha = lr.Ralpha;
  tr_outside_hb(cm, dsq, L, r, y, i0, j0, BE_EFFICIENT, fill_L, fill_R,
		NULL, &beta, NULL, &betaL, NULL, &betaR, NULL, NULL,
		&boutsc, &boutv, &boutj, &boutmode, cp9b);

  /* brief 26_0430-227: recompute the banded EL decks' per-row d-band edges (same
   * vroot=r/vend=y that tr_outside_hb used) to clamp the EL-terminus reads below
   * and free the banded decks; the helper is deterministic so this matches the
   * writer's allocation exactly. */
  if (cm->flags & CMH_LOCAL_END) {
    if ((eldmaxJ = malloc(sizeof(int) * (L+1))) == NULL) cm_Fail("brief 26_0430-227: eldmaxJ OOM");
    tr_outside_hb_el_dmax(cm, L, r, y, i0, j0, 0, cp9b, eldmaxJ);
    if (fill_L) { if ((eldmaxL = malloc(sizeof(int) * (L+1))) == NULL) cm_Fail("eldmaxL OOM"); tr_outside_hb_el_dmax(cm, L, r, y, i0, j0, 1, cp9b, eldmaxL); }
    if (fill_R) { if ((eldmaxR = malloc(sizeof(int) * (L+1))) == NULL) cm_Fail("eldmaxR OOM"); tr_outside_hb_el_dmax(cm, L, r, y, i0, j0, 2, cp9b, eldmaxR); }
  }

  /* best inside begin (child local hit), across J/L/R */
  binsc = binJsc; binv = binJ; binmode = TRMODE_J;
  if (fill_L && lr.Lbsc > binsc) { binsc = lr.Lbsc; binv = lr.Lb; binmode = TRMODE_L; }
  if (fill_R && lr.Rbsc > binsc) { binsc = lr.Rbsc; binv = lr.Rb; binmode = TRMODE_R; }

  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  best_v  = -99;
  p_mode  = c_mode = TRMODE_J;
  for (v = w; v <= y; v++)
    for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++)
      {
	jp   = j - (i0-1);
	jp_v = j - jmin[v];
	for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	  {
	    int dp_v = d - hd_min(cp9b, v, jp_v);
	    int haveL = fill_L && NOT_IMPOSSIBLE(betaL[v][j][dp_v]);
	    int haveR = fill_R && NOT_IMPOSSIBLE(betaR[v][j][dp_v]);
	    /* J */
	    if ((sc = alpha[v][j][dp_v] + beta[v][j][dp_v]) > best_sc)
	      { best_sc=sc; best_v=v; best_d=d; best_j=j; p_mode=TRMODE_J; c_mode=TRMODE_J; }
	    /* L: parent v in L, child L (plane-consistent).  brief 26_0610-093: the
	     * former mismatched candidate (alpha[v]+betaL[v], p_mode=L/c_mode=J) is
	     * REMOVED. In a linear wedge stretch a state v carries a single mode, so its
	     * inside*outside product is meaningful only plane-consistently (Lalpha[v]+
	     * betaL[v] for v-in-L; alpha[v]+beta[v] for v-in-J). A marginal->J truncation
	     * conversion needs no mismatched split: at-or-below v it is already inside
	     * Lalpha[v] (tr_inside_hb ML/IL own-mode L<-J term), and above v it is already
	     * inside beta[v] (tr_outside_hb feeds the child J-plane from an L/R-mode parent,
	     * e.g. :7495/:7530/:7556). alpha[v]+betaL[v] pairs a J-inside subtree with an
	     * L-outside path -- no single-mode parse -- so it can only over-score (betaL[v]
	     * drops right-side emission mass beta[v] keeps) and its c_mode=J routing rebuilds
	     * v's subtree in J -> the sample9 4-node phantom (092 root-caused this). */
	    if (haveL && cp9b->Lvalid[v] && (sc = Lalpha[v][j][dp_v] + betaL[v][j][dp_v]) > best_sc)
	      { best_sc=sc; best_v=v; best_d=d; best_j=j; p_mode=TRMODE_L; c_mode=TRMODE_L; }
	    /* R: parent v in R, child R (plane-consistent; R mirror -- the mismatched
	     * alpha[v]+betaR[v] candidate is REMOVED for the same reason, brief 26_0610-093). */
	    if (haveR && cp9b->Rvalid[v] && (sc = Ralpha[v][j][dp_v] + betaR[v][j][dp_v]) > best_sc)
	      { best_sc=sc; best_v=v; best_d=d; best_j=j; p_mode=TRMODE_R; c_mode=TRMODE_R; }
	  }
      }

  if (cm->flags & CMH_LOCAL_END) {
    /* J EL terminus, NOT gated on r_allow_J (brief 26_0610-096 -- the mirror of the
     * tr_generic_splitter_hb note: a J EL terminus is a legal winner in any mode; the
     * sample13 fix lives in tr_vinside_hb's D/S marginal-EL reconstruction, not a gate). */
    /* brief 26_0430-227: EL decks are banded (row j: d=0..eldmax*[j], NULL rows
     * where eldmax*[j]<0); clamp the read to eldmax*[j] (cells above were
     * IMPOSSIBLE, so best_sc is byte-identical). */
    for (jp = 0; jp <= W; jp++)
      {
	int dhi;
	j = i0-1+jp;
	if (eldmaxJ[j] < 0) continue;
	dhi = ESL_MIN(jp, eldmaxJ[j]);
	for (d = 0; d <= dhi; d++)
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc; best_v = -1; best_j = j; best_d = d; p_mode = TRMODE_J; c_mode = TRMODE_T;
	  }
      }
    /* brief 26_0610-090: L/R marginal EL candidates (parent v in L/R descends to EL,
     * child empty). betaL/betaR[cm->M] now built above. best_v=-1 with p_mode=L/R;
     * traceback routes through tr_v_splitter_hb (parent-mode) -> tr_vinsideT_hb, which
     * already reconstructs the marginal EL leaf (USED_EL in the active mode). */
    if (fill_L)
      for (jp = 0; jp <= W; jp++)
	{
	  int dhi;
	  j = i0-1+jp;
	  if (eldmaxL[j] < 0) continue;
	  dhi = ESL_MIN(jp, eldmaxL[j]);
	  for (d = 0; d <= dhi; d++)
	    if ((sc = betaL[cm->M][j][d]) > best_sc) {
	      best_sc = sc; best_v = -1; best_j = j; best_d = d; p_mode = TRMODE_L; c_mode = TRMODE_T;
	    }
	}
    if (fill_R)
      for (jp = 0; jp <= W; jp++)
	{
	  int dhi;
	  j = i0-1+jp;
	  if (eldmaxR[j] < 0) continue;
	  dhi = ESL_MIN(jp, eldmaxR[j]);
	  for (d = 0; d <= dhi; d++)
	    if ((sc = betaR[cm->M][j][d]) > best_sc) {
	      best_sc = sc; best_v = -1; best_j = j; best_d = d; p_mode = TRMODE_R; c_mode = TRMODE_T;
	    }
	}
  }

  if (r==0) {
    if (binsc > best_sc) { best_sc = binsc; best_v = -2; best_j = j0; best_d = W; p_mode = TRMODE_T; c_mode = binmode; }
  }
  /* outside marginal terminus (parent local hit, child empty). Use >= so the
   * terminus WINS ties vs a standard split that lands on the same marginal-end
   * cell: the oracle inside DP forces the marginal-end at d<2 (it OVERWRITES the
   * transit), so a tied "split-and-continue-below" must collapse to the terminus
   * (else the wedge below emits spurious all-delete-to-end nodes -- brief 26_0610-050). */
  if (boutsc >= best_sc) { best_sc = boutsc; best_v = -3; best_j = boutj; best_d = 1; p_mode = boutmode; c_mode = TRMODE_T; }

  /* brief 26_0430-227: free the three banded EL decks explicitly (eldmax-driven
   * width + correct cyk_dnc accounting) before the generic matrix free, then the
   * rest.  betaL/betaR here are fresh arrays (tr_outside_hb got NULL inputs), so
   * unlike outside_hb's generic caller there is no alpha aliasing to avoid. */
  if (cm->flags & CMH_LOCAL_END) {
    if (beta[cm->M]  != NULL) { free_el_banded_vjd_deck(beta[cm->M],  i0, j0, eldmaxJ); beta[cm->M]  = NULL; }
    if (fill_L && betaL[cm->M] != NULL) { free_el_banded_vjd_deck(betaL[cm->M], i0, j0, eldmaxL); betaL[cm->M] = NULL; }
    if (fill_R && betaR[cm->M] != NULL) { free_el_banded_vjd_deck(betaR[cm->M], i0, j0, eldmaxR); betaR[cm->M] = NULL; }
  }
  free_banded_hb_vjd_matrix(alpha, cm, i0, j0, cp9b);
  free_banded_hb_vjd_matrix(beta,  cm, i0, j0, cp9b);
  if (fill_L) { free_banded_hb_vjd_matrix(Lalpha, cm, i0, j0, cp9b); free_banded_hb_vjd_matrix(betaL, cm, i0, j0, cp9b); }
  if (fill_R) { free_banded_hb_vjd_matrix(Ralpha, cm, i0, j0, cp9b); free_banded_hb_vjd_matrix(betaR, cm, i0, j0, cp9b); }
  free(eldmaxJ); free(eldmaxL); free(eldmaxR);

  /* TRUNCATED infeasibility guard: nothing in-band beat IMPOSSIBLE. */
  if (best_v == -99) return best_sc;

  if (best_v == -1) {           /* EL (parent p_mode to EL, child empty) */
    /* brief 26_0610-090: route in the winning parent mode. z_allow=p_mode mirrors the
     * -3 marginal-terminus routing; for p_mode==J this is (TRUE,FALSE,FALSE) -- exactly
     * the prior J-only behavior, so no J regression. */
    tr_v_splitter_hb(cm, dsq, L, tr, r, w, i0, best_j-best_d+1, best_j, j0, TRUE,
		     r_allow_J, r_allow_L, r_allow_R,
		     (p_mode==TRMODE_J), (p_mode==TRMODE_L), (p_mode==TRMODE_R), cp9b);
    return best_sc;
  }
  if (best_v == -2) {           /* inside begin (parent empty), child local hit */
    InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, binv, binmode);
    tr_wedge_splitter_hb(cm, dsq, L, tr, binv, z, i0, j0,
			 (binmode==TRMODE_J), (binmode==TRMODE_L), (binmode==TRMODE_R), cp9b);
    return best_sc;
  }
  if (best_v == -3) {           /* outside terminus (parent marginal, child empty) */
    tr_v_splitter_hb(cm, dsq, L, tr, r, boutv, i0, best_j, best_j, j0, FALSE,
		     r_allow_J, r_allow_L, r_allow_R,
		     (p_mode==TRMODE_J), (p_mode==TRMODE_L), (p_mode==TRMODE_R), cp9b);
    return best_sc;
  }

  /* standard split: V-problem r..v (bottom mode c_mode), wedge v..z (mode c_mode) */
  tr_v_splitter_hb(cm, dsq, L, tr, r, best_v, i0, best_j-best_d+1, best_j, j0, FALSE,
		   r_allow_J, r_allow_L, r_allow_R,
		   (c_mode==TRMODE_J), (c_mode==TRMODE_L), (c_mode==TRMODE_R), cp9b);
  tr_wedge_splitter_hb(cm, dsq, L, tr, best_v, z, best_j-best_d+1, best_j,
		       (c_mode==TRMODE_J), (c_mode==TRMODE_L), (c_mode==TRMODE_R), cp9b);
  return best_sc;
}

/* Function: tr_generic_splitter_hb()  [brief 26_0610-044 R4.4a (J); brief 26_0610-046 R4.4b-pt2 (L/R)]
 *
 * The mode-aware banded truncated generic splitter. J keeps the exact R4.4a logic
 * (best_k>=0 bifurcation / -1 EL / -2,-3 child begins). L/R add: marginal split
 * cases (parent v in L or R, children J or marginal) scored against the 1-D
 * marginal outside betas (betaL/betaR), the marginal child begins, and -4 (the
 * "local hit in parent" marginal terminus from the outside b_sc). v_mode/w_mode/
 * y_mode track the resolved per-node modes (an empty marginal child is carried as
 * the parent mode v_mode, brief 26_0610-050). T (brief 26_0610-051, R4.4c): T appears only at B
 * states, reached only by a truncated begin into the B at full span -> a T split
 * (v_mode=T, w=R-left BEGL, y=L-right BEGR) whose "outside" is just the begin
 * scalar (no T outside recurrence), plus deeper-B T begins via b1/b2. The T split
 * attaches root->B(T) + children directly (no V-problem: the parent region above a
 * T B is degenerate). Ported from truncyk.c::tr_generic_splitter, gated by the bands.
 */
static float
tr_generic_splitter_hb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
		       int r, int z, int i0, int j0,
		       int r_allow_J, int r_allow_L, int r_allow_R, int r_allow_T, CP9Bands_t *cp9b)
{
  float ***alpha = NULL, ***Lalpha = NULL, ***Ralpha = NULL;
  float ***beta  = NULL;
  float ***betaL = NULL, ***betaR = NULL;
  /* T (brief 26_0610-051) forces the L+R planes on (the T-combine reads R(BEGL)+L(BEGR)). */
  int      fill_L = r_allow_L || r_allow_T, fill_R = r_allow_R || r_allow_T;
  int      fill_T = r_allow_T;
  int      v,w,y;
  int      wend, yend;
  int      jp;
  int      jp_v;
  int      W;
  float    sc;
  int      j,d,k;
  float    best_sc;
  int      best_k, best_d, best_j;
  int      v_mode, w_mode, y_mode;
  int      tv;
  TR_LR    lr1, lr2;
  int      b1J, b2J;       float b1Jsc, b2Jsc;
  int      b1_v, b2_v, b1_mode, b2_mode;   float b1_sc, b2_sc;
  float    b3_sc;  int b3_v, b3_j, b3_mode;
  int     *jmin  = cp9b->jmin;
  int     *jmax  = cp9b->jmax;
  int    **hdmin = cp9b->hdmin;
  int    **hdmax = cp9b->hdmax;
  int     *eldmaxJ = NULL, *eldmaxL = NULL, *eldmaxR = NULL;  /* brief 26_0430-227: banded EL d-band edges (recomputed) */

  if (insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT) {
    sc = tr_insideT_hb(cm, dsq, L, tr, r, z, i0, j0, (r==0), r_allow_J, r_allow_L, r_allow_R, r_allow_T, cp9b);
    return sc;
  }

  for (v = r; v <= z-5; v++)
    if (cm->sttype[v] == B_st) break;

  if (v > z-5) {
    /* No bifurcation in r..z. T is only valid at B states, so a T problem always
     * has a bifurcation -> this branch is reached only for J/L/R (the R-left/L-right
     * subtrees of a T split recurse here as pure R/L). */
    if (cm->sttype[z] != E_st) cm_Fail("inconceivable.");
    sc = tr_wedge_splitter_hb(cm, dsq, L, tr, r, z, i0, j0, r_allow_J, r_allow_L, r_allow_R, cp9b);
    return sc;
  }

  w = cm->cfirst[v];
  y = cm->cnum[v];
  if (w < y) { wend = y-1; yend = z; }
  else       { yend = w-1; wend = z; }

  /* inside for both child subtrees (J + L/R planes [+ T begins]), chaining the arrays */
  lr1.fill_L = fill_L; lr1.fill_R = fill_R; lr1.fill_T = fill_T; lr1.ret_planes = TRUE; lr1.Lalpha = NULL; lr1.Ralpha = NULL;
  tr_inside_hb(cm, dsq, L, w, wend, i0, j0, BE_EFFICIENT, NULL,  &alpha, NULL, NULL, NULL,
	       (r==0), &b1J, &b1Jsc, &lr1, cp9b);
  lr2.fill_L = fill_L; lr2.fill_R = fill_R; lr2.fill_T = fill_T; lr2.ret_planes = TRUE; lr2.Lalpha = lr1.Lalpha; lr2.Ralpha = lr1.Ralpha;
  tr_inside_hb(cm, dsq, L, y, yend, i0, j0, BE_EFFICIENT, alpha, &alpha, NULL, NULL, NULL,
	       (r==0), &b2J, &b2Jsc, &lr2, cp9b);
  Lalpha = lr2.Lalpha; Ralpha = lr2.Ralpha;
  /* outside for the parent (J reuses the alpha array; 1-D L/R betas + terminus b3) */
  tr_outside_hb(cm, dsq, L, r, v, i0, j0, BE_EFFICIENT, fill_L, fill_R,
		alpha, &beta, NULL, &betaL, NULL, &betaR, NULL, NULL,
		&b3_sc, &b3_v, &b3_j, &b3_mode, cp9b);

  /* brief 26_0430-227: recompute the banded EL decks' per-row d-band edges (same
   * vroot=r/vend=v tr_outside_hb used) to clamp the EL-terminus reads and free the
   * banded decks; deterministic helper, so it matches the writer's allocation. */
  if (cm->flags & CMH_LOCAL_END) {
    if ((eldmaxJ = malloc(sizeof(int) * (L+1))) == NULL) cm_Fail("brief 26_0430-227: eldmaxJ OOM");
    tr_outside_hb_el_dmax(cm, L, r, v, i0, j0, 0, cp9b, eldmaxJ);
    if (fill_L) { if ((eldmaxL = malloc(sizeof(int) * (L+1))) == NULL) cm_Fail("eldmaxL OOM"); tr_outside_hb_el_dmax(cm, L, r, v, i0, j0, 1, cp9b, eldmaxL); }
    if (fill_R) { if ((eldmaxR = malloc(sizeof(int) * (L+1))) == NULL) cm_Fail("eldmaxR OOM"); tr_outside_hb_el_dmax(cm, L, r, v, i0, j0, 2, cp9b, eldmaxR); }
  }

  /* best inside begin per child subtree (across J/L/R/T). A T begin targets a
   * deeper bifurcation inside the subtree (brief 26_0610-051); it competes with J/L/R. */
  b1_sc = b1Jsc; b1_v = b1J; b1_mode = TRMODE_J;
  if (fill_L && lr1.Lbsc > b1_sc) { b1_sc = lr1.Lbsc; b1_v = lr1.Lb; b1_mode = TRMODE_L; }
  if (fill_R && lr1.Rbsc > b1_sc) { b1_sc = lr1.Rbsc; b1_v = lr1.Rb; b1_mode = TRMODE_R; }
  if (fill_T && lr1.Tbsc > b1_sc) { b1_sc = lr1.Tbsc; b1_v = lr1.Tb; b1_mode = TRMODE_T; }
  b2_sc = b2Jsc; b2_v = b2J; b2_mode = TRMODE_J;
  if (fill_L && lr2.Lbsc > b2_sc) { b2_sc = lr2.Lbsc; b2_v = lr2.Lb; b2_mode = TRMODE_L; }
  if (fill_R && lr2.Rbsc > b2_sc) { b2_sc = lr2.Rbsc; b2_v = lr2.Rb; b2_mode = TRMODE_R; }
  if (fill_T && lr2.Tbsc > b2_sc) { b2_sc = lr2.Tbsc; b2_v = lr2.Tb; b2_mode = TRMODE_T; }

  W = j0-i0+1;
  best_sc = IMPOSSIBLE;
  best_k  = -99;
  v_mode  = w_mode = y_mode = TRMODE_J;
  for (j = ESL_MAX(i0-1, jmin[v]); j <= ESL_MIN(j0, jmax[v]); j++)
    {
      jp   = j - (i0-1);
      jp_v = j - jmin[v];
      for (d = hd_min(cp9b, v, jp_v); d <= hd_max(cp9b, v, jp_v) && d <= jp; d++)
	{
	  int dp_v = d - hd_min(cp9b, v, jp_v);
	  int haveL = fill_L && NOT_IMPOSSIBLE(betaL[v][j][dp_v]);
	  int haveR = fill_R && NOT_IMPOSSIBLE(betaR[v][j][dp_v]);
	  /* T (brief 26_0610-051): the only T "outside" is the begin scalar at the B's
	   * full-span cell (oracle cm_TrCYKOutsideAlignHB:6743): trpenalty at the
	   * root (r==0), else 0 (deeper-B recursion; penalty already in b1/b2).
	   * No T outside recurrence -> T contributes only at (j0,W). */
	  float Tbeta_v = IMPOSSIBLE;
	  if (fill_T && cp9b->Tvalid[v] && j == j0 && d == W) {
	    if (r == 0) { if (cp9b->Tvalid[0]) Tbeta_v = tr_trpenalty(cm, v); }
	    else        Tbeta_v = 0.;
	  }
	  int haveT = NOT_IMPOSSIBLE(Tbeta_v);
	  for (k = 0; k <= d; k++)
	    {
	      int dp_w, dp_y;
	      int inw = hb_inband(cp9b, w, j-k, d-k, i0, j0, &dp_w);
	      int iny = hb_inband(cp9b, y, j,   k,   i0, j0, &dp_y);
	      if (! (inw && iny)) continue;
	      /* all-J split.  brief 26_0610-075: this candidate combines v's own classical
	       * Outside (beta[v]) with both children's classical Inside (alpha[w], alpha[y]) --
	       * exactly the 3-participant classical/J cross-term 070 gated everywhere else in
	       * this function, but 070's audit did not reach this specific (unmarked) candidate.
	       * Without the gate, a B state whose own node structurally cannot be reached in
	       * pure J mode (cp9b->Jvalid[v]==FALSE) can still win here purely because alpha[]/
	       * beta[] are always filled (never masked to IMPOSSIBLE) regardless of Jvalid, and
	       * best_sc starts at IMPOSSIBLE -- producing a (bkind,J) pin the downstream pinned
	       * checkpoint posterior engine (pin_tr_inside_B) correctly refuses to allocate a
	       * deck for, causing the NULL deref one layer downstream (mirrors 070's own
	       * diagnosis of the analogous bug at its other sites). */
	      if (cp9b->Jvalid[v] && cp9b->Jvalid[w] && cp9b->Jvalid[y] &&
		  (sc = alpha[w][j-k][dp_w] + alpha[y][j][dp_y] + beta[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=k; best_j=j; best_d=d; v_mode=TRMODE_J; w_mode=TRMODE_J; y_mode=TRMODE_J; }
	      /* L: v in L; w=J, y=L  (k>0).  brief 26_0610-070: w's J-plane cross-term needs
	       * cp9b->Jvalid[w] -- w's node must lie entirely within the observed
	       * sequence for its plain (truncation-unaware) alpha[w] to be a valid
	       * contribution (mirrors the do_J_y convention in cm_dpalign_trunc.c). */
	      if (haveL && k > 0 && cp9b->Lvalid[y] && cp9b->Jvalid[w] &&
		  (sc = alpha[w][j-k][dp_w] + Lalpha[y][j][dp_y] + betaL[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=k; best_j=j; best_d=d; v_mode=TRMODE_L; w_mode=TRMODE_J; y_mode=TRMODE_L; }
	      /* L: v in L; w=J, y=J  (k>0).  brief 26_0610-070: gate both children's J-planes. */
	      if (haveL && k > 0 && cp9b->Jvalid[w] && cp9b->Jvalid[y] &&
		  (sc = alpha[w][j-k][dp_w] + alpha[y][j][dp_y] + betaL[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=k; best_j=j; best_d=d; v_mode=TRMODE_L; w_mode=TRMODE_J; y_mode=TRMODE_J; }
	      /* R: v in R; w=R, y=J  (k<d).  brief 26_0610-070: gate y's J-plane. */
	      if (haveR && k < d && cp9b->Rvalid[w] && cp9b->Jvalid[y] &&
		  (sc = Ralpha[w][j-k][dp_w] + alpha[y][j][dp_y] + betaR[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=k; best_j=j; best_d=d; v_mode=TRMODE_R; w_mode=TRMODE_R; y_mode=TRMODE_J; }
	      /* R: v in R; w=J, y=J  (k<d).  brief 26_0610-070: gate both children's J-planes. */
	      if (haveR && k < d && cp9b->Jvalid[w] && cp9b->Jvalid[y] &&
		  (sc = alpha[w][j-k][dp_w] + alpha[y][j][dp_y] + betaR[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=k; best_j=j; best_d=d; v_mode=TRMODE_R; w_mode=TRMODE_J; y_mode=TRMODE_J; }
	      /* T: v in T; w=R (left BEGL), y=L (right BEGR); 1<=k<=d-1 (both children
	       * non-empty). Tbeta_v is the begin scalar (only at j0,W). (brief 26_0610-051) */
	      if (haveT && k > 0 && k < d && cp9b->Rvalid[w] && cp9b->Lvalid[y] &&
		  (sc = Ralpha[w][j-k][dp_w] + Lalpha[y][j][dp_y] + Tbeta_v) > best_sc)
		{ best_sc=sc; best_k=k; best_j=j; best_d=d; v_mode=TRMODE_T; w_mode=TRMODE_R; y_mode=TRMODE_L; }
	    }
	  /* k=0: left child covers all d, right child empty (L marginal) */
	  if (haveL) {
	    int dp_w;
	    if (hb_inband(cp9b, w, j, d, i0, j0, &dp_w)) {
	      if (cp9b->Lvalid[w] && (sc = Lalpha[w][j][dp_w] + betaL[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=0; best_j=j; best_d=d; v_mode=TRMODE_L; w_mode=TRMODE_L; y_mode=TRMODE_T; }
	      /* brief 26_0610-070: w's J-plane cross-term needs cp9b->Jvalid[w] (LEFT_FULL
	       * mirror of the RIGHT_FULL crash site below). */
	      if (cp9b->Jvalid[w] && (sc = alpha[w][j][dp_w] + betaL[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=0; best_j=j; best_d=d; v_mode=TRMODE_L; w_mode=TRMODE_J; y_mode=TRMODE_T; }
	    }
	  }
	  /* k=d: right child covers all d, left child empty (R marginal) */
	  if (haveR) {
	    int dp_y;
	    if (hb_inband(cp9b, y, j, d, i0, j0, &dp_y)) {
	      if (cp9b->Rvalid[y] && (sc = Ralpha[y][j][dp_y] + betaR[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=d; best_j=j; best_d=d; v_mode=TRMODE_R; w_mode=TRMODE_T; y_mode=TRMODE_R; }
	      /* brief 26_0610-070: y's J-plane cross-term needs cp9b->Jvalid[y] -- the
	       * tr_generic_splitter_hb (true D&C, large-subtree) mirror of the
	       * tr_inside_hb "B special case 2" crash site (frag5p_s12 L24). */
	      if (cp9b->Jvalid[y] && (sc = alpha[y][j][dp_y] + betaR[v][j][dp_v]) > best_sc)
		{ best_sc=sc; best_k=d; best_j=j; best_d=d; v_mode=TRMODE_R; w_mode=TRMODE_T; y_mode=TRMODE_J; }
	    }
	  }
	}
    }

  if (cm->flags & CMH_LOCAL_END) {
    /* J EL terminus (classical beta[cm->M]). brief 26_0610-096: this candidate is NOT
     * gated on r_allow_J -- a J-mode EL terminus is a legitimate winner even in an
     * L/R solve (the oracle's bulk-init USED_EL admits it in every mode; ar45-sample26_p60
     * L-solve's true optimum uses exactly this candidate). An earlier 096 draft gated it
     * on r_allow_J to keep the CsrB-sample13_5ptr V-problem's dispatch mode consistent,
     * but that dropped ar45-sample26_p60 by 0.6 bits (flipped L->R). The real sample13
     * fix is the D/S marginal-EL reconstruction in tr_vinside_hb (see there): the R-mode
     * traceback terminates at the BEGL_S->EL leaf regardless of which (tied) EL candidate
     * won here, so no gate is needed and none is correct. */
    /* brief 26_0430-227: EL decks banded; clamp reads to eldmax*[j] (byte-identical). */
    for (jp = 0; jp <= W; jp++)
      {
	int dhi;
	j = i0-1+jp;
	if (eldmaxJ[j] < 0) continue;
	dhi = ESL_MIN(jp, eldmaxJ[j]);
	for (d = 0; d <= dhi; d++)
	  if ((sc = beta[cm->M][j][d]) > best_sc) {
	    best_sc = sc; best_k = -1; best_j = j; best_d = d;
	    v_mode = TRMODE_J; w_mode = TRMODE_T; y_mode = TRMODE_T;
	  }
      }
    /* brief 26_0610-090: L/R marginal EL candidates (parent v in L/R descends to EL,
     * child empty). Gated on r_allow_L/r_allow_R -- NOT fill_L/fill_R (which are also
     * ON for T, since T needs the L+R planes for its combine): an L/R EL terminus is
     * only a valid winner in an actual L/R solve, never in a T solve (whose optimum is
     * the bifurcation begin). This keeps J and T solves byte-identical. */
    if (r_allow_L)
      for (jp = 0; jp <= W; jp++)
	{
	  int dhi;
	  j = i0-1+jp;
	  if (eldmaxL[j] < 0) continue;
	  dhi = ESL_MIN(jp, eldmaxL[j]);
	  for (d = 0; d <= dhi; d++)
	    if ((sc = betaL[cm->M][j][d]) > best_sc) {
	      best_sc = sc; best_k = -1; best_j = j; best_d = d;
	      v_mode = TRMODE_L; w_mode = TRMODE_T; y_mode = TRMODE_T;
	    }
	}
    if (r_allow_R)
      for (jp = 0; jp <= W; jp++)
	{
	  int dhi;
	  j = i0-1+jp;
	  if (eldmaxR[j] < 0) continue;
	  dhi = ESL_MIN(jp, eldmaxR[j]);
	  for (d = 0; d <= dhi; d++)
	    if ((sc = betaR[cm->M][j][d]) > best_sc) {
	      best_sc = sc; best_k = -1; best_j = j; best_d = d;
	      v_mode = TRMODE_R; w_mode = TRMODE_T; y_mode = TRMODE_T;
	    }
	}
  }

  if (r == 0) {
    if (b1_sc > best_sc) { best_sc = b1_sc; best_k = -2; best_j = j0; best_d = W; }
    if (b2_sc > best_sc) { best_sc = b2_sc; best_k = -3; best_j = j0; best_d = W; }
  }
  /* "local hit in parent" marginal terminus (outside b_sc): the parse stays in
   * the parent r..v region marginally and never reaches the bifurcation. */
  if (b3_sc > best_sc) { best_sc = b3_sc; best_k = -4; best_j = b3_j; v_mode = b3_mode; }

  /* brief 26_0430-227: free the banded EL decks explicitly (eldmax-driven width +
   * correct cyk_dnc accounting).  beta IS alpha here (tr_outside_hb reused the
   * inside alpha array in place, cf. the J-feed call above), so the J EL deck
   * lives in the shared alpha[cm->M]==beta[cm->M] slot -- free it there, NULL the
   * slot, and never free beta separately (that would double-free alpha).
   * betaL/betaR are fresh arrays (NULL inputs), freed normally. */
  if (cm->flags & CMH_LOCAL_END) {
    if (beta[cm->M] != NULL) { free_el_banded_vjd_deck(beta[cm->M], i0, j0, eldmaxJ); beta[cm->M] = NULL; }  /* == alpha[cm->M] */
    if (fill_L && betaL[cm->M] != NULL) { free_el_banded_vjd_deck(betaL[cm->M], i0, j0, eldmaxL); betaL[cm->M] = NULL; }
    if (fill_R && betaR[cm->M] != NULL) { free_el_banded_vjd_deck(betaR[cm->M], i0, j0, eldmaxR); betaR[cm->M] = NULL; }
  }
  free_banded_hb_vjd_matrix(alpha, cm, i0, j0, cp9b);
  if (fill_L) { free_banded_hb_vjd_matrix(Lalpha, cm, i0, j0, cp9b); free_banded_hb_vjd_matrix(betaL, cm, i0, j0, cp9b); }
  if (fill_R) { free_banded_hb_vjd_matrix(Ralpha, cm, i0, j0, cp9b); free_banded_hb_vjd_matrix(betaR, cm, i0, j0, cp9b); }
  free(eldmaxJ); free(eldmaxL); free(eldmaxR);

  /* TRUNCATED infeasibility guard: no in-band split/EL/begin/terminus beat
   * IMPOSSIBLE -> no valid parse under the bands; return cleanly (no garbage recurse). */
  if (best_k == -99) return best_sc;

  if (best_k == -1) {
    /* brief 26_0610-090: route in the winning parent mode v_mode (z_allow=v_mode).
     * For v_mode==J this is (TRUE,FALSE,FALSE) -- exactly the prior J-only behavior. */
    tr_v_splitter_hb(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, TRUE,
		     r_allow_J, r_allow_L, r_allow_R,
		     (v_mode==TRMODE_J), (v_mode==TRMODE_L), (v_mode==TRMODE_R), cp9b);
    return best_sc;
  }
  if (best_k == -2) {
    InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b1_v, b1_mode);
    z = CMSubtreeFindEnd(cm, b1_v);
    tr_generic_splitter_hb(cm, dsq, L, tr, b1_v, z, i0, j0,
			   (b1_mode==TRMODE_J), (b1_mode==TRMODE_L), (b1_mode==TRMODE_R), (b1_mode==TRMODE_T), cp9b);
    return best_sc;
  }
  if (best_k == -3) {
    InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b2_v, b2_mode);
    z = CMSubtreeFindEnd(cm, b2_v);
    tr_generic_splitter_hb(cm, dsq, L, tr, b2_v, z, i0, j0,
			   (b2_mode==TRMODE_J), (b2_mode==TRMODE_L), (b2_mode==TRMODE_R), (b2_mode==TRMODE_T), cp9b);
    return best_sc;
  }
  if (best_k == -4) {
    /* marginal terminus: solve the parent V-problem r..b3_v with a degenerate hole */
    tr_v_splitter_hb(cm, dsq, L, tr, r, b3_v, i0, best_j, best_j, j0, FALSE,
		     r_allow_J, r_allow_L, r_allow_R,
		     (v_mode==TRMODE_J), (v_mode==TRMODE_L), (v_mode==TRMODE_R), cp9b);
    return best_sc;
  }

  /* T bifurcation (brief 26_0610-051): the parse begins (T) directly into the B state v,
   * splitting into left child BEGL in R mode + right child BEGR in L mode. The
   * parent region r..v is degenerate (the truncated begin jumps straight to v), so
   * we attach v + children directly -- NO V-problem (the begin has no descent above
   * a B). best_k in [1,d-1]. At r==0 attach the root->v begin here; for a deeper-B
   * recursion (r!=0, v==r) v was already attached by the b1/b2 begin pre-attach. */
  if (v_mode == TRMODE_T) {
    if (r == 0)
      InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, v, TRMODE_T);
    tv = tr->n-1;
    InsertTraceNodewithMode(tr, tv, TRACE_LEFT_CHILD,  best_j-best_d+1, best_j-best_k, w, TRMODE_R);
    tr_generic_splitter_hb(cm, dsq, L, tr, w, wend, best_j-best_d+1, best_j-best_k,
			   FALSE, FALSE, TRUE, FALSE, cp9b);
    InsertTraceNodewithMode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j,       y, TRMODE_L);
    tr_generic_splitter_hb(cm, dsq, L, tr, y, yend, best_j-best_k+1, best_j,
			   FALSE, TRUE, FALSE, FALSE, cp9b);
    return best_sc;
  }

  /* bifurcation split (best_k in [0,d]); one child may be empty (mode TRMODE_T) */
  tr_v_splitter_hb(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, FALSE,
		   r_allow_J, r_allow_L, r_allow_R,
		   (v_mode==TRMODE_J), (v_mode==TRMODE_L), (v_mode==TRMODE_R), cp9b);
  tv = tr->n-1;

  if (w_mode != TRMODE_T) {
    InsertTraceNodewithMode(tr, tv, TRACE_LEFT_CHILD, best_j-best_d+1, best_j-best_k, w, w_mode);
    tr_generic_splitter_hb(cm, dsq, L, tr, w, wend, best_j-best_d+1, best_j-best_k,
			   (w_mode==TRMODE_J), (w_mode==TRMODE_L), (w_mode==TRMODE_R), FALSE, cp9b);
  } else
    /* empty (truncated-away) marginal child: the oracle (tr_insideT_hb) labels it
     * with the PARENT marginal mode (v_mode), not TRMODE_T (brief 26_0610-050). */
    InsertTraceNodewithMode(tr, tv, TRACE_LEFT_CHILD, best_j-best_d+1, best_j-best_d, w, v_mode);

  if (y_mode != TRMODE_T) {
    InsertTraceNodewithMode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j, y, y_mode);
    tr_generic_splitter_hb(cm, dsq, L, tr, y, yend, best_j-best_k+1, best_j,
			   (y_mode==TRMODE_J), (y_mode==TRMODE_L), (y_mode==TRMODE_R), FALSE, cp9b);
  } else
    InsertTraceNodewithMode(tr, tv, TRACE_RIGHT_CHILD, best_j-best_k+1, best_j, y, v_mode);

  return best_sc;
}

/* Function: TrCYKDivideAndConquerHB()  [brief 26_0610-044 R4.4a; 045 R4.4b L/R; 051 R4.4c T]
 *
 * Purpose:  HMM-banded truncated divide-and-conquer CYK. <preset_mode> selects
 *           the marginal mode to solve (TRMODE_J / TRMODE_L / TRMODE_R / TRMODE_T;
 *           all four modes are now supported -> the truncated D&C CYK is complete).
 *           Returns that mode in <ret_mode> and the optimal parsetree
 *           (is_std=FALSE, pass_idx + trpenalty set, every node mode-tagged) in
 *           <ret_tr>. The returned score includes the (penalty-folded) truncated
 *           begin so it equals the oracle cm_TrCYKInsideAlignHB()'s
 *           {J,L,R,T}alpha[0][L][L] byte-for-byte for the chosen mode.
 *
 * LOCAL-BEGIN inflation bug (brief 26_0610-081 root-caused; FIXED brief 26_0610-086):
 * in LOCAL config (CMH_LOCAL_BEGIN) this function previously could return an
 * unreachable, too-high score -- e.g. a tRNA case reported 29.4 bits (forced J)
 * when the trusted oracle cm_TrAlignHB() proves the true forced-J optimum is 6.5
 * bits. Root cause (081, confirmed live via gdb; 086, confirmed vs the oracle's
 * own cm_TrCYKOutsideAlignHB beta): tr_outside_hb()'s pure-J outside propagation
 * (and its traceback twin, tr_vinside_hb's pure-J recurrence) was NOT gated on
 * cp9b->Jvalid[v] -- unlike the truncated-begin injection, the marginal->J
 * cross-terms (brief 070), and the L/R planes, which all already were. So a
 * truncated-begin penalty could chain a run of (good) emissions THROUGH a state
 * whose node structurally cannot appear in a pure-J parse (Jvalid[v]==FALSE; the
 * oracle never even allocates its J deck), inflating tr_generic_splitter_hb's
 * "all-J split" (beta[v]-based) candidate. Concretely, for the tRNA repro the
 * oracle's own beta[54][69][60] = -22.23 (duality: +28.75 inside = 6.52) but the
 * ungated D&C propagated +0.66 -- exactly the ~23-bit inflation. The fix is the
 * do_J_v = cp9b->Jvalid[v] gate on the four pure-J beta[v] writes in
 * tr_outside_hb (~line 7363) plus the Jvalid[v] mask on tr_vinside_hb's pure-J
 * plane. Validated byte-exact vs the oracle (argmax ncand-loop) across a
 * 1000-seq tRNA local panel + 5S + a bps=0 matl300 panel; global (-g) unaffected
 * (Jvalid all-TRUE there -> the gate is a no-op). NOTE (081/086): a blanket
 * removal of the injection is NOT the fix -- it regresses legitimate wedge-entry
 * candidates (e.g. the same case's forced-J optimum, 6.5 bits at entry state
 * v=12). See subagent-summaries/081_2026-07-11_* and 086_2026-07-13_* in the
 * 26_0610 notebook dir. (Separately, forcing a NON-optimal preset_mode can still
 * "leak" the J winner via the all-J-split candidate not being gated on
 * r_allow_J -- a pre-existing, out-of-086-scope issue that does not affect the
 * argmax ncand-loop production/fallback path.)
 */
float
TrCYKDivideAndConquerHB(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, int pass_idx,
			char preset_mode, char *ret_mode, Parsetree_t **ret_tr, CP9Bands_t *cp9b)
{
  Parsetree_t *tr;
  float        sc;
  int          z;
  int          b;
  int          r_allow_J = (preset_mode == TRMODE_J);
  int          r_allow_L = (preset_mode == TRMODE_L);
  int          r_allow_R = (preset_mode == TRMODE_R);
  int          r_allow_T = (preset_mode == TRMODE_T);   /* brief 26_0610-051 (R4.4c) */

  if (cp9b == NULL) cm_Fail("TrCYKDivideAndConquerHB(): cp9b is NULL");
  if (r != 0)       cm_Fail("TrCYKDivideAndConquerHB(): r must be 0 (truncated begins enter from ROOT_S)");
  if (! (r_allow_J || r_allow_L || r_allow_R || r_allow_T))
    cm_Fail("TrCYKDivideAndConquerHB(): preset_mode must be J, L, R, or T");

  /* set the file-static truncation-penalty context (read by tr_trpenalty()) */
  if ((tr_dnc_pty_idx = cm_tr_penalties_IdxForPass(pass_idx)) == -1)
    cm_Fail("TrCYKDivideAndConquerHB(): unexpected pass idx: %d", pass_idx);
  tr_dnc_local = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;

  tr = CreateParsetree(100);
  tr->is_std   = FALSE;
  tr->pass_idx = pass_idx;
  /* root node, in the preset mode. The truncated begin target <b> is recovered
   * from the parse's first non-root node after D&C completes. */
  InsertTraceNodewithMode(tr, -1, TRACE_LEFT_CHILD, i0, j0, 0, preset_mode);
  z = cm->M-1;

  /* J and L/R all route through the mode-aware HMM-banded D&C generic splitter.
   * J is R4.4a (mode J throughout). L/R (R4.4b): the 2-D marginal outside + combine
   * (049, SCORE byte-exact) drives the split, and the FULL mode-tracking traceback
   * is now byte-exact too (brief 26_0610-050, R4.4b-part-2c-ii): the part-2c-i SCORE-only
   * guard is removed. The two part-2c-ii parse fixes were (a) empty (truncated-away)
   * marginal bifurcation children must carry the PARENT marginal mode (v_mode), not
   * TRMODE_T, and (b) the marginal terminus (outside b_sc) must WIN ties vs a
   * standard split landing on the same marginal-end cell (oracle forces the
   * marginal-end at d<2). T is OFF (R4.4c). Driver-only entry (rung-3/4 not yet in
   * the cmalign dispatch on this branch). */
  /* brief 26_0430-225: DNC_MEM_VERBOSE ground-truth peak-memory diagnostic,
   * truncated D&C entry -- same reused CYKDeckTrack + CYKShadowTrack
   * instrumentation as CYKDivideAndConquerHB() above (that infra already
   * covers this trunc code path too: tr_inside_hb/tr_outside_hb/
   * tr_generic_splitter_hb allocate via the same alloc_banded_hb_vjd_deck-
   * family calls, already instrumented). */
  int dnc_mem_verbose = (getenv("DNC_MEM_VERBOSE") != NULL) ? TRUE : FALSE;
  if (dnc_mem_verbose) { CYKDeckTrackReset(); CYKShadowTrackReset(); }

  sc = tr_generic_splitter_hb(cm, dsq, L, tr, 0, z, i0, j0, r_allow_J, r_allow_L, r_allow_R, r_allow_T, cp9b);

  if (dnc_mem_verbose) {
    double peak_mb = CYKDeckTrackMaxMb();
    double vjd_mb  = CYKDeckTrackVjdAtPeakMb();
    double vji_mb  = CYKDeckTrackVjiAtPeakMb();
    double shad_mb = CYKShadowTrackMaxMb();
    fprintf(stderr, "# TrCYKDivideAndConquerHB (D&C, trunc) engaged: M=%d L=%d mode=%c sc=%.5f  DnC-DP peak=%.2f Mb (vjd=%.2f Mb vji=%.2f Mb)  shadow-peak=%.2f Mb\n",
            cm->M, L, preset_mode, sc, peak_mb, vjd_mb, vji_mb, shad_mb);
  }

  /* the truncated-begin entry state is the first state attached below ROOT_S */
  b = (tr->n > 1) ? tr->state[1] : 0;
  tr->trpenalty = tr_trpenalty(cm, b);

  if (ret_mode != NULL) *ret_mode = preset_mode;
  if (ret_tr   != NULL) *ret_tr = tr; else FreeParsetree(tr);
  return sc;
}

/*****************************************************************
 * Brief 26_0430-225: cm_DnCAlignSizeNeededHB() / cm_TrDnCAlignSizeNeededHB()
 * -- pre-alignment memory estimators for the HMM-banded divide-and-conquer
 * engine (CYKDivideAndConquerHB() / TrCYKDivideAndConquerHB()).
 *
 * UNLIKE the --ckpt estimators (cm_dpalign.c / cm_dpalign_trunc.c), this is
 * NOT an exact structural replay: D&C's actual recursion tree shape is
 * DATA-DEPENDENT.  generic_splitter_hb() (cm_dpsmall.c:4656-4799) picks its
 * split point (best_j/best_k) by DP SCORE at each bifurcation
 * (cm_dpsmall.c:4721-4746) -- which residues end up left vs. right of a
 * given bifurcation, and therefore how big each recursive child's banded
 * sub-matrix (alpha[w..wend] + alpha[y..yend], restricted to the chosen
 * i0..j0 subrange) actually is, cannot be known without running the DP.
 * So this function computes a data-INDEPENDENT UPPER BOUND instead: for
 * EVERY bifurcation state v in the CM, the byte cost if its two children's
 * node ranges were each allocated over the FULL sequence range [1,L] (the
 * real recursion always uses some SUBRANGE of [1,L] for any nested call,
 * so this is >= the true cost for that bifurcation); take the max over all
 * bifurcations.  This is the right kind of estimate for an --mxsize gate
 * (which must not under-promise), but is expected to overestimate,
 * possibly by a lot, on deeply unbalanced trees -- brief 226's job is to
 * measure the actual overestimation factor via DNC_MEM_VERBOSE.
 *
 * KNOWN GAP (verified, not assumed -- flagged for brief 226): beta[v] (the
 * bifurcation-state Outside deck, cm_dpsmall.c:4719/4737) is allocated by
 * outside_hb() and is NEVER FREED anywhere in generic_splitter_hb() (grep-
 * confirmed: only `alpha` gets a free_banded_hb_vjd_matrix() call, at
 * cm_dpsmall.c:4769; `beta` has no matching free in the function body).
 * This means beta[v] decks accumulate, unfreed, across EVERY bifurcation
 * visited over the WHOLE recursion -- not just the single widest live
 * sub-problem this estimator bounds.  For a CM with many bifurcations
 * (e.g. an rRNA model), the true peak measured by DNC_MEM_VERBOSE could
 * exceed this bound by O(#bifurcations-visited x avg-beta-deck-size),
 * which is NOT captured below (beta[v] is a single-state deck, individually
 * small, but the accumulation is unbounded across a long alignment run).
 * This is an existing property of CYKDivideAndConquerHB()/
 * generic_splitter_hb() uncovered while writing this estimator, not
 * something this brief introduces or fixes -- out of scope to fix here,
 * but the self-check below directly tests whether it matters in practice.
 *
 * Also NOT modeled: the class-2 "V-problem" (wedge) contribution
 * (v_splitter()'s EXACT unbanded triangular decks, tracked separately by
 * the existing cyk_dnc_vji_bytes counter) -- computing its size requires
 * modeling insideT_size()'s RAMLIMIT cutoff and the V-problem's own d-range,
 * which this estimator does not attempt.  DNC_MEM_VERBOSE's vji=... report
 * captures it in the ground truth; this estimator's vjd-only bound may
 * therefore UNDERESTIMATE on CMs/sequences where the V-problem dominates.
 *****************************************************************/

/* sum of nc(v) = per-(v,jp) hd band-width, v in [lo..hi], full L bands
 * (same per-cell formula cm_hb_mx_SizeNeeded_ex uses, cm_mx.c:1149-1183). */
static int64_t
mxest_dnc_range_nc(CP9Bands_t *cp9b, int lo, int hi)
{
  int v, jp;
  int64_t nc = 0;
  for (v = lo; v <= hi; v++) {
    int njr = cp9b->jmax[v] - cp9b->jmin[v] + 1; if (njr < 0) njr = 0;
    for (jp = 0; jp < njr; jp++) { int w = hd_max(cp9b, v, jp) - hd_min(cp9b, v, jp) + 1; if (w > 0) nc += w; }
  }
  return nc;
}

/* Brief 26_0430-228: banded upper bound on outside_hb()'s single EL deck
 * (state cm->M), replacing the pre-227 full-triangle size_vjd_deck(L,1,L)
 * term now that outside_hb_el_dmax()/alloc_el_banded_vjd_deck() (cm_dpsmall.c
 * :5329/5375, brief 227) band the real engine's EL allocation to eldmax[]
 * rows.  Calling outside_hb_el_dmax() with (vroot=0,vend=cm->M-1,i0=1,j0=L)
 * -- the same (vroot,vend,i0,j0) CYKDivideAndConquerHB()'s top-level call
 * uses (cm_dpsmall.c:544-566: z=cm->M-1, i0/j0 passed straight through,
 * r defaults to 0) and also the widest possible sub-problem span -- gives a
 * per-row eldmax[] that safely bounds every narrower recursive sub-call's
 * eldmax too: a smaller vroot..vend main-loop range only removes v's that
 * can contribute to eldmax[j], and a narrower [i0,j0] only shrinks
 * jp = j-(i0-1), and both changes can only lower eldhi, never raise it. So
 * this single data-independent call gives a safe upper bound, using the
 * exact per-row cost function the real allocator uses (sum over eldmax[r]>=0
 * rows of (eldmax[r]+1) floats). */
static float
mxest_dnc_el_mb(CM_t *cm, int L, CP9Bands_t *cp9b)
{
  int      status;
  int     *eldmax = NULL;
  int      r;
  int64_t  nfloats = 0;
  ESL_ALLOC(eldmax, sizeof(int) * (L+1));
  outside_hb_el_dmax(cm, L, 0, cm->M-1, 1, L, cp9b, eldmax);
  for (r = 0; r <= L; r++) if (eldmax[r] >= 0) nfloats += (int64_t) (eldmax[r]+1);
  free(eldmax);
  return (float) (nfloats * sizeof(float) / 1000000.);
 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* Truncated analogue of mxest_dnc_el_mb(), for tr_outside_hb()'s up-to-three
 * concurrently-live EL decks (J always; L/R per fill_L/fill_R, mirroring
 * tr_outside_hb()'s per-plane allocation at cm_dpsmall.c ~7617-7636).  Each
 * plane's eldmax[] is computed independently by tr_outside_hb_el_dmax() with
 * the same (vroot=0,vend=cm->M-1,i0=1,j0=L) top-level/widest-span argument
 * as mxest_dnc_el_mb() uses, for the same reason (a safe upper bound over
 * every narrower recursive sub-call). */
static float
mxest_dnc_tr_el_mb(CM_t *cm, int L, CP9Bands_t *cp9b, int fill_L, int fill_R)
{
  int      status;
  int     *eldmax = NULL;
  int      r, plane;
  int64_t  nfloats;
  float    totmb = 0.;
  ESL_ALLOC(eldmax, sizeof(int) * (L+1));
  for (plane = 0; plane <= 2; plane++) {
    if (plane == 1 && ! fill_L) continue;
    if (plane == 2 && ! fill_R) continue;
    tr_outside_hb_el_dmax(cm, L, 0, cm->M-1, 1, L, plane, cp9b, eldmax);
    nfloats = 0;
    for (r = 0; r <= L; r++) if (eldmax[r] >= 0) nfloats += (int64_t) (eldmax[r]+1);
    totmb += (float) (nfloats * sizeof(float) / 1000000.);
  }
  free(eldmax);
  return totmb;
 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* Function: cm_DnCAlignSizeNeededHB()
 * Incept:   Brief 26_0430-225
 *
 * Purpose:  Predict an UPPER BOUND (not an exact prediction -- see the file
 *           header above) on the peak Mb CYKDivideAndConquerHB() will need
 *           to align a length-<L> sequence to <cm> under its current
 *           cm->cp9b bands, without running any alignment.
 *
 * Args:     cm, errbuf, L - usual
 *           ret_vjdmb - RETURN: bound on class-1 (banded vjd alpha+beta) peak Mb
 *           ret_shmb  - RETURN: bound on shadow-deck (yshadow/kshadow) Mb at that same peak
 *           ret_totmb - RETURN: ret_vjdmb + ret_shmb
 *
 * Returns:  <eslOK> on success; <eslEINCOMPAT> if cm->cp9b is NULL.
 */
int
cm_DnCAlignSizeNeededHB(CM_t *cm, char *errbuf, int L, float *ret_vjdmb, float *ret_shmb, float *ret_totmb)
{
  int status;
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_DnCAlignSizeNeededHB(): cm->cp9b is NULL");
  CP9Bands_t *cp9b = cm->cp9b;
  int v, w, y, wend, yend;
  int64_t best_vjd_nc = 0, best_sh_bytes = 0;
  /* brief 26_0430-228: outside_hb()'s EL deck (state cm->M) is now BANDED
   * (brief 227: outside_hb_el_dmax()/alloc_el_banded_vjd_deck(), cm_dpsmall.c
   * :5329/5375) rather than the old full-triangle alloc_vjd_deck() this
   * comment used to describe -- see mxest_dnc_el_mb() above for the banded
   * upper-bound replacement (was size_vjd_deck(L,1,L) pre-228). */
  int have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  float elmb  = have_el ? mxest_dnc_el_mb(cm, L, cp9b) : 0.;

  int has_bif = FALSE;
  for (v = 0; v < cm->M; v++) {
    if (cm->sttype[v] != B_st) continue;
    has_bif = TRUE;
    w = cm->cfirst[v]; y = cm->cnum[v];
    if (w < y) { wend = y-1; yend = cm->M-1; } else { yend = w-1; wend = cm->M-1; }
    /* NOTE: wend/yend above use cm->M-1 as a safe (over-)estimate of each
     * subtree's true end state (CMSubtreeFindEnd(cm,w)/(cm,y) would be
     * tighter, but the wider range is still a valid upper bound and avoids
     * a second full state-array walk per bifurcation). */
    int64_t nc_w = mxest_dnc_range_nc(cp9b, w, ESL_MIN(wend, CMSubtreeFindEnd(cm, w)));
    int64_t nc_y = mxest_dnc_range_nc(cp9b, y, ESL_MIN(yend, CMSubtreeFindEnd(cm, y)));
    int64_t nc_v = mxest_dnc_range_nc(cp9b, v, v); /* beta[v]: single-state deck */
    int64_t vjd_nc = nc_w + nc_y + nc_v;
    if (vjd_nc > best_vjd_nc) {
      best_vjd_nc = vjd_nc;
      /* shadow cost at this same split: yshadow (char) for every non-B state
       * in [w..wend]+[y..yend], kshadow (int) for any B states among them
       * (alloc_banded_hb_vjd_yshadow_deck/kshadow_deck, cm_dpsmall.c:4433-4479) */
      int64_t sh_bytes = 0, vv;
      int wend_t = ESL_MIN(wend, CMSubtreeFindEnd(cm, w)), yend_t = ESL_MIN(yend, CMSubtreeFindEnd(cm, y));
      for (vv = w; vv <= wend_t; vv++) sh_bytes += mxest_dnc_range_nc(cp9b, vv, vv) * (cm->sttype[vv] == B_st ? sizeof(int) : sizeof(char));
      for (vv = y; vv <= yend_t; vv++) sh_bytes += mxest_dnc_range_nc(cp9b, vv, vv) * (cm->sttype[vv] == B_st ? sizeof(int) : sizeof(char));
      best_sh_bytes = sh_bytes;
    }
  }
  if (! has_bif) {
    /* no bifurcations: D&C never splits the node range, so its peak equals
     * the ordinary full (non-checkpointed) banded matrix -- same formula
     * cm_hb_mx_SizeNeeded() uses (cm_mx.c:1149-1183). */
    best_vjd_nc = mxest_dnc_range_nc(cp9b, 0, cm->M-1);
    best_sh_bytes = mxest_dnc_range_nc(cp9b, 0, cm->M-1) * sizeof(char);
  }

  float vjdmb = (float) (best_vjd_nc * sizeof(float) / 1000000.) + elmb;
  float shmb  = (float) (best_sh_bytes / 1000000.);
  if (ret_vjdmb != NULL) *ret_vjdmb = vjdmb;
  if (ret_shmb  != NULL) *ret_shmb  = shmb;
  if (ret_totmb != NULL) *ret_totmb = vjdmb + shmb;
  return eslOK;
}

/* Function: cm_TrDnCAlignSizeNeededHB()
 * Incept:   Brief 26_0430-225
 *
 * Purpose:  Truncated analogue of cm_DnCAlignSizeNeededHB(), for
 *           TrCYKDivideAndConquerHB() / tr_generic_splitter_hb().  Same
 *           upper-bound strategy, J/L/R-tripled per bifurcation child
 *           exactly as tr_generic_splitter_hb() allocates alpha+Lalpha+
 *           Ralpha together (cm_dpsmall.c:9305-9312: tr_inside_hb() with
 *           lr1.ret_planes=TRUE returns all three requested planes at
 *           once for the SAME [w..wend]/[y..yend] node range).  <preset_mode>
 *           selects fill_L/fill_R exactly as cm_TrFillFromMode() does.
 *
 * Args:     cm, errbuf, L, preset_mode - usual (preset_mode: TRMODE_J/L/R/T)
 *           ret_vjdmb, ret_shmb, ret_totmb - as cm_DnCAlignSizeNeededHB()
 *
 * Returns:  <eslOK> on success; <eslEINCOMPAT>/other on bad preset_mode or NULL cp9b.
 */
int
cm_TrDnCAlignSizeNeededHB(CM_t *cm, char *errbuf, int L, char preset_mode, float *ret_vjdmb, float *ret_shmb, float *ret_totmb)
{
  int status;
  int fill_L, fill_R, fill_T;
  if (cm->cp9b == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm_TrDnCAlignSizeNeededHB(): cm->cp9b is NULL");
  if ((status = cm_TrFillFromMode(preset_mode, &fill_L, &fill_R, &fill_T)) != eslOK)
    ESL_FAIL(status, errbuf, "cm_TrDnCAlignSizeNeededHB(): bad preset_mode");
  CP9Bands_t *cp9b = cm->cp9b;
  int planes = 1 + (fill_L?1:0) + (fill_R?1:0); /* J always; L/R per mode -- cm_dpsmall.c:9264 */
  int v, w, y, wend, yend;
  int64_t best_vjd_nc = 0, best_sh_bytes = 0;
  /* brief 26_0430-228: tr_outside_hb() now allocates one BANDED EL deck per
   * active plane (J always, L/R per fill_L/fill_R -- brief 227:
   * tr_outside_hb_el_dmax()/alloc_el_banded_vjd_deck(), cm_dpsmall.c
   * :5443/7617-7636) rather than the old full-triangle decks this comment
   * used to describe -- see mxest_dnc_tr_el_mb() above for the banded
   * upper-bound replacement (was (float)planes * size_vjd_deck(L,1,L)
   * pre-228). */
  int have_el = (cm->flags & CMH_LOCAL_END) ? TRUE : FALSE;
  float elmb  = have_el ? mxest_dnc_tr_el_mb(cm, L, cp9b, fill_L, fill_R) : 0.;

  int has_bif = FALSE;
  for (v = 0; v < cm->M; v++) {
    if (cm->sttype[v] != B_st) continue;
    has_bif = TRUE;
    w = cm->cfirst[v]; y = cm->cnum[v];
    if (w < y) { wend = y-1; yend = cm->M-1; } else { yend = w-1; wend = cm->M-1; }
    int wend_t = ESL_MIN(wend, CMSubtreeFindEnd(cm, w)), yend_t = ESL_MIN(yend, CMSubtreeFindEnd(cm, y));
    int64_t nc_w = mxest_dnc_range_nc(cp9b, w, wend_t);
    int64_t nc_y = mxest_dnc_range_nc(cp9b, y, yend_t);
    int64_t nc_v = mxest_dnc_range_nc(cp9b, v, v);
    int64_t vjd_nc = (int64_t)planes * (nc_w + nc_y) + nc_v; /* beta[v] (J-only 1-D outside) not tripled */
    if (vjd_nc > best_vjd_nc) {
      best_vjd_nc = vjd_nc;
      int64_t sh_bytes = 0, vv;
      for (vv = w; vv <= wend_t; vv++) sh_bytes += (int64_t)planes * mxest_dnc_range_nc(cp9b, vv, vv) * (cm->sttype[vv] == B_st ? sizeof(int) : sizeof(char));
      for (vv = y; vv <= yend_t; vv++) sh_bytes += (int64_t)planes * mxest_dnc_range_nc(cp9b, vv, vv) * (cm->sttype[vv] == B_st ? sizeof(int) : sizeof(char));
      best_sh_bytes = sh_bytes;
    }
  }
  if (! has_bif) {
    best_vjd_nc = (int64_t)planes * mxest_dnc_range_nc(cp9b, 0, cm->M-1);
    best_sh_bytes = (int64_t)planes * mxest_dnc_range_nc(cp9b, 0, cm->M-1) * sizeof(char);
  }

  float vjdmb = (float) (best_vjd_nc * sizeof(float) / 1000000.) + elmb;
  float shmb  = (float) (best_sh_bytes / 1000000.);
  if (ret_vjdmb != NULL) *ret_vjdmb = vjdmb;
  if (ret_shmb  != NULL) *ret_shmb  = shmb;
  if (ret_totmb != NULL) *ret_totmb = vjdmb + shmb;
  return eslOK;
}

/* Function: vsplitter_b()
 *           EPN 05.19.05
 * *based on vsplitter(), only difference is bands are used :
 *
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
 * Args:     cm          -  model
 *           dsq         - digitized sequence 1..L
 *           L           - length of dsq
 *           tr          - the traceback we're adding on to.
 *           r           - index of the first state in the subgraph 
 *           z           - index of the last state in the subgraph
 *           i0,i1       - first part of the subsequence (1..L)
 *           j1,j0       - second part of the subsequence (1..L)
 *           useEL       - TRUE if i1,j1 aligned to EL, not z
 *           dmin        - minimum d bound for each state v; [0..v..M-1]
 *           dmax        - maximum d bound for each state v; [0..v..M-1]
 * 
 * Returns:  (void)
 */
static void
v_splitter_qdb(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr,
	   int r, int z, int i0, int i1, int j1, int j0, 
	   int useEL, int *dmin, int *dmax)
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
  int   b;			/* optimal choice for a 0->b local begin  */
  float bsc;			/* score if we use the local begin */
  int  *imin = NULL;            /* minimum i bound for each state v; [0..y-w] 
                                 * calculated using *dmin; offset from v, the
				 * band that corresponds to state v, is imin[v-w] */
  int  *imax = NULL;            /* maximum i bound for each state v; [0..y-w] 
                                 * calculated using *dmax; offset from v, the
				 * band that corresponds to state v, is imax[v-w] */


  /* 1. If the V problem is either a boundary condition, or small
   *    enough, solve it with v_inside^T and append the trace to tr.
   *    (With local alignment, we might even see a lone B state
   *     get handed to v_splitter(); hence the r==z case.)
   */
   if (cm->ndidx[z] == cm->ndidx[r] + 1 || r == z || 
      vinsideT_size(cm, r, z, i0, i1, j1, j0) < RAMLIMIT)
    {
      ESL_DPRINTF2(("#DEBUG: Solving a V:   G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		z, UniqueStatetype(cm->stid[z]),
		i0,j1,j1,j0));
      vinsideT(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, (r==0), dmin, dmax);
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
   *    beta[cm->M] is the EL deck, needed for local ends.
   */
  vinside_qdb (cm, dsq, L, w, z, i0, i1, j1, j0, useEL, BE_EFFICIENT, 
	     NULL, &alpha, NULL, &pool, NULL, (r==0), &b, &bsc,
	     dmin, dmax);
  voutside_qdb(cm, dsq, L, r, y, i0, i1, j1, j0, useEL, BE_EFFICIENT, 
	     NULL, &beta,  pool, NULL, dmin, dmax);

  /* 4. Find the optimal split: v, ip, jp. 
   */
  /* Bands used ip 1A */
  imin = malloc(sizeof (int) * (y-w+1));
  imax = malloc(sizeof (int) * (y-w+1));

  best_sc = IMPOSSIBLE;
  for (v = w; v <= y; v++)
    {
      /* Bands used ip 1B */

      /* Fill imin[v-w] and imax[v-w] as we go, one of many ways to handle imin and imax 
       * Remember state indices in imin and imax are offset from v because imin and 
       * imax run [0..y-w], ==> dmin[v] corresponds to imin[v-w] 
       */

      imin[v-w] = j1-i0-dmax[v]+1;
      imax[v-w] = j1-i0-dmin[v]+1;

      /*orig lines : for (ip = 0; ip <= i1-i0; ip++) 
       *                    for (jp = 0; jp <= j0-j1; jp++)
       *the order is switched here because the band on ip depends
       *on jp.
       */
      for (jp = 0; jp <= j0-j1; jp++)
	{
	  if((imin[v-w]+jp) < 0) ip = 0;
	  else ip = imin[v-w]+jp;
	  for (; (ip <= imax[v-w]+jp) && ip <= (i1-i0); ip++) 
	    if ((sc = alpha[v][jp][ip] + beta[v][jp][ip]) > best_sc)
	      {
		best_sc = sc;
		best_v  = v;
		best_i  = ip + i0;
		best_j  = jp + j1;
	      }
	}
    }
  /* Local alignment ends: maybe we're better off in EL, not
   * the split set?
   */
  if (useEL && (cm->flags & CMH_LOCAL_END)) {
    /* There is no band on the EL state */
    for (ip = 0; ip <= i1-i0; ip++)
      for (jp = 0; jp <= j0-j1; jp++)
	if ((sc = beta[cm->M][jp][ip]) > best_sc) {
	  best_sc = sc;
	  best_v  = -1;
	  best_i  = ip + i0;
	  best_j  = jp + j1;
	}
  }
	
  /* Local alignment begins: maybe we're better off in root...
   */
  if (r==0 && (cm->flags & CMH_LOCAL_BEGIN)) {
    if (bsc > best_sc) {
      best_sc = bsc;
      best_v  = -2;
      best_i  = i0;
      best_j  = j0;
    }
  }

  /* Free now, before recursing!
   */
  free_vji_matrix(alpha, cm->M, j1, j0);
  free_vji_matrix(beta,  cm->M, j1, j0);

  /* If we're in EL, instead of the split set, the optimal
   * alignment is entirely in a V problem that's still above us.
   * The TRUE flag sets useEL; we propagate allow_begin. 
   */
  if (best_v == -1) {
    v_splitter_qdb(cm, dsq, L, tr, r, w, i0, best_i, best_j, j0, TRUE, dmin, dmax);    
    if(imin != NULL) free(imin);
    if(imax != NULL) free(imax);
    return;
  }

  /* If we used a local begin, the optimal alignment is
   * entirely in a V problem that's still below us, rooted
   * at b, for the entire one-hole sequence. The FALSE
   * flag prohibits more local begin transitions; we propagate
   * useEL.
   */
  if (best_v == -2) {
    if (b != z) 
      InsertTraceNode(tr, tr->n-1, TRACE_LEFT_CHILD, i0, j0, b);
    v_splitter_qdb(cm, dsq, L, tr, b, z, i0, i1, j1, j0, useEL, dmin, dmax);    
    return;
  }

  /* The optimal split into two V problems:
   *    V:   r..v, i0..i', j'..j0
   *    V:   v..z, i'..i1, j1..j'
   * Solve in this order, because we're constructing the
   * trace in postorder traversal.
   */
  ESL_DPRINTF2(("#DEBUG: V splitter:\n"));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		r, UniqueStatetype(cm->stid[r]),
		best_v, UniqueStatetype(cm->stid[best_v]),
		i0, best_i, best_j, j0));
  ESL_DPRINTF2(("#DEBUG:    V:       G%d[%s]..%d[%s], %d..%d//%d..%d\n", 
		best_v, UniqueStatetype(cm->stid[best_v]),
		z, UniqueStatetype(cm->stid[z]),
		best_i, i1, j1, best_j));

  v_splitter_qdb(cm, dsq, L, tr, r,      best_v, i0,     best_i, best_j, j0, FALSE,
	       dmin, dmax);
  v_splitter_qdb(cm, dsq, L, tr, best_v, z,      best_i, i1,     j1,     best_j, useEL,
	       dmin, dmax);
  
  free(imax);
  free(imin);
  return;
}


/*****************************************************************
 * The alignment engines, using bands:
 *     inside_qdb   - given generic or wedge problem G^r_z to i0..j0, return score and matrix
 *     outside_qdb  - given unbifurcated G^r_z to i0..j0, return matrix
 *     
 *     vinside_qdb  - given V problem G^r_z to i0..i1//j1..j0, return score and matrix
 *     voutside_qdb - given unbifurcated G^r_z to i0..i1//j1..j0, return matrix
 ******************************************************************/

/* Function: inside_qdb()
 *           EPN 05.19.05
 * *based on inside(), only difference is bands are used : 
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  (See inside())
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
 *           dmin   - minimum d bound for each state v; [0..v..M-1]
 *           dmax   - maximum d bound for each state v; [0..v..M-1]
 *                       
 * Returns: Score of the optimal alignment.  
 */
static float 
inside_qdb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	 float ***alpha, float ****ret_alpha, 
	 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	 void ****ret_shadow, 
	 int allow_begin, int *ret_b, float *ret_bsc,
	 int *dmin, int *dmax)
{
  int      status;
  float  **end;         /* we re-use the end deck. */
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
  int      kmax;        /* for B_st's, maximum k value consistent with bands*/
  
  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();
  if (! deckpool_pop(dpool, &end))
    end = alloc_vjd_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, vroot, E_st);
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
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch,  sizeof(int) * cm->M);
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

      /* Impose bands by setting all cells outside the bands to 0 
       * This is independent of state type so we do it outside
       * the following set of if then statements. 
       */

      for (jp = 0; jp <= W; jp++) {
	j = i0-1+jp;
	for (d = 0; d < dmin[v] && d <= jp; d++)
	  alpha[v][j][d] = IMPOSSIBLE;
	for (d = dmax[v]+1; d <= jp;     d++) 
	  alpha[v][j][d] = IMPOSSIBLE;
      }
      
      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][d]  = USED_EL; 
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
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
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		z = cm->cnum[v];
		/* Careful, in qdb, we only want to look at alpha cells that are
		 * within the bands for all states involved (v, y and z) */
		/* k is the length of the right fragment */
		if(dmin[z] > (d-dmax[y])) k = dmin[z];
		else k = d-dmax[y];
		if(k < 0) k = 0;
		
		if(dmax[z] < (d-dmin[y])) kmax = dmax[z];
		else kmax = d-dmin[y];
		
		if(k <= kmax)
		  {
		    alpha[v][j][d] = alpha[y][j-k][d-k] + alpha[z][j][k];
		    if (ret_shadow != NULL) kshad[j][d] = k;
		    for (k=k+1; k <= kmax; k++)
		      {
			if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d]) {
			  alpha[v][j][d] = sc;
			  if (ret_shadow != NULL) kshad[j][d] = k;
			}
		      }
		  }
		else alpha[v][j][d] = IMPOSSIBLE;
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
	    /* dmin[v] must be >= 2 */
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-2] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		
		i = j-d+1;
		if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    /* dmin[v] must be >= 1 */
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) >  alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  } 
		
		i = j-d+1;
		if (dsq[i] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][dsq[i]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    alpha[v][j][0] = IMPOSSIBLE;
	    /* dmin[v] must be >= 1 */
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][d] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d]) {
		    alpha[v][j][d] = sc;
		    if (ret_shadow != NULL) yshad[j][d] = yoffset;
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][d] += cm->esc[v][dsq[j]];
		else
		  alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}				/* finished calculating deck v. */
      
      /* Check for local begin getting us to the root.
       * This is "off-shadow": if/when we trace back, we'll handle this
       * case separately (and we'll know to do it because we'll immediately
       * see a USED_LOCAL_BEGIN flag in the shadow matrix, telling us
       * to jump right to state b; see below)
       */
      if (allow_begin && alpha[v][j0][W] + cm->beginsc[v] > bsc) 
	{
	  b   = v;
	  bsc = alpha[v][j0][W] + cm->beginsc[v];
	}

      /* Check for whether we need to store an optimal local begin score
       * as the optimal overall score, and if we need to put a flag
       * in the shadow matrix telling insideT() to use the b we return.
       */
      if (allow_begin && v == 0 && bsc > alpha[0][j0][W]) {
	alpha[0][j0][W] = bsc;
	if (ret_shadow != NULL) yshad[j0][W] = USED_LOCAL_BEGIN;
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

  /* Now we free our memory. 
   * if we've got do_full set, all decks vroot..vend are now valid (end is shared).
   * else, only vroot deck is valid now and all others vroot+1..vend are NULL, 
   * and end is NULL.
   * We could check this status to be sure (and we used to) but now we trust. 
   */
  sc       = alpha[vroot][j0][W];
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

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}


/* Function: outside_qdb()
 *           EPN 05.19.05
 * *based on outside(), only difference is bands are used : 
 *
 * Date:     SRE, Tue Aug  8 10:42:52 2000 [St. Louis]
 * Purpose:  (See outside())
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
 *           dmin      - minimum d bound for each state v; [0..v..M-1]
 *           dmax      - maximum d bound for each state v; [0..v..M-1]
 */
static void
outside_qdb(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0,
	  int do_full, float ***beta, float ****ret_beta,
	  struct deckpool_s *dpool, struct deckpool_s **ret_dpool, int *dmin, int *dmax)
{
  int      status;
  int      v,y;			/* indices for states */
  int      j,d,i;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */
  int      W;			/* subsequence length */
  int      jp;			/* j': relative position in the subsequence, 0..W */
  int      voffset;		/* index of v in t_v(y) transition scores */
  int      w1,w2;		/* bounds of split set */
  int      dv;                  /* StateDelta() for state v */

  /* Allocations and initializations
   */
  W = j0-i0+1;		/* the length of the subsequence: used in many loops */

  			/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* if caller didn't give us a matrix, make one.
   * Allocate room for M+1 decks because we might need the EL deck (M)
   * if we're doing local alignment.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
  }

  /* Initialize the root deck.
   * If the root is in a split set, initialize the whole split set.
   */
  w1 = cm->nodemap[cm->ndidx[vroot]]; /* first state in split set */
  if (cm->sttype[vroot] == B_st) {    /* special boundary case of Outside on a single B state. */
    w2 = w1;
    if (vend != vroot) cm_Fail("oh no. not again.");
  } else
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

  /* Initialize the EL deck at M, if we're doing local alignment w.r.t. ends.
   */
  if (cm->flags & CMH_LOCAL_END) {
    if (! deckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = alloc_vjd_deck(L, i0, j0);
    for (jp = 0; jp <= W; jp++) {
      j = i0-1+jp;
      for (d = 0; d <= jp; d++)
	beta[cm->M][j][d] = IMPOSSIBLE;
    }
    
    /* We have to worry about vroot -> EL transitions.
     * since we start the main recursion at w2+1. This requires a 
     * laborious partial unroll of the main recursion, grabbing
     * the stuff relevant to a beta[EL] calculation for just the
     * vroot->EL transition.
     */
    if (NOT_IMPOSSIBLE(cm->endsc[vroot])) {
      switch (cm->sttype[vroot]) {
      case MP_st:
	if (W < 2) break;
	if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	  escore = cm->esc[vroot][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
	else
	  escore = DegeneratePairScore(cm->abc, cm->esc[vroot], dsq[i0], dsq[j0]);
	beta[cm->M][j0-1][W-2] = cm->endsc[vroot] + 
	  (cm->el_selfsc * (W-2)) + escore;
	if (beta[cm->M][j0-1][W-2] < IMPOSSIBLE) beta[cm->M][j0-1][W-2] = IMPOSSIBLE;
	break;
      case ML_st:
      case IL_st:
	if (W < 1) break;
	if (dsq[i0] < cm->abc->K) 
	  escore = cm->esc[vroot][(int) dsq[i0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[vroot]);
	beta[cm->M][j0][W-1] = cm->endsc[vroot] + 
	  (cm->el_selfsc * (W-1)) + escore;
	if (beta[cm->M][j0][W-1] < IMPOSSIBLE) beta[cm->M][j0][W-1] = IMPOSSIBLE;
	break;
      case MR_st:
      case IR_st:
	if (W < 1) break;
	if (dsq[j0] < cm->abc->K) 
	  escore = cm->esc[vroot][(int) dsq[j0]];
	else
	  escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[vroot]);
	beta[cm->M][j0-1][W-1] = cm->endsc[vroot] + 
	  (cm->el_selfsc * (W-1)) + escore;
	if (beta[cm->M][j0-1][W-1] < IMPOSSIBLE) beta[cm->M][j0-1][W-1] = IMPOSSIBLE;
	break;
      case S_st:
      case D_st:
	beta[cm->M][j0][W] = cm->endsc[vroot] + 
	  (cm->el_selfsc * W);
	if (beta[cm->M][j0][W] < IMPOSSIBLE) beta[cm->M][j0][W] = IMPOSSIBLE;
	break;
      case B_st:		/* can't start w/ bifurcation at vroot. */
      default: cm_Fail("bogus parent state %d\n", cm->sttype[vroot]);
      }
    }
  }
  
  ESL_ALLOC(touch, sizeof(int) * cm->M);
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
      if ((vroot == 0 && i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN))
	  && (dmin[v] <= W && dmax[v] >= W))
	  beta[v][j0][W] = cm->beginsc[v];

      /* main recursion:
       */
      for (jp = W; jp >= 0; jp--) {
	j = i0-1+jp;
	if((dmax[v]) > jp) d = jp;
	else d = (dmax[v]);
	for (; d >= (dmin[v]); d--)
	  {
	    i = j-d+1;
	    for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	      if (y < vroot) continue; /* deal with split sets */
	      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

	      switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == j0 || d == jp) continue; /* boundary condition */

		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
		
		if ((sc = beta[y][j+1][d+2] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;

	      case ML_st:
	      case IL_st: 
		if (d == jp) continue;	/* boundary condition (note when j=0, d=0*/

		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		  
		if ((sc = beta[y][j][d+1] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		  
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);

		if ((sc = beta[y][j+1][d+1] + cm->tsc[y][voffset] + escore) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		if ((sc = beta[y][j][d] + cm->tsc[y][voffset]) > beta[v][j][d])
		  beta[v][j][d] = sc;
		break;

	      default: cm_Fail("bogus child state %d\n", cm->sttype[y]);
	      }/* end switch over states*/
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] = IMPOSSIBLE;
	    
	    
	  } /* ends loop over d. We know all beta[v][j][d] in this row j*/
      }/* end loop over jp. We know the beta's for the whole deck.*/
      
      /* Deal with local alignment end transitions v->EL
       * (EL = deck at M.)
       */
      if (NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (jp = 0; jp <= W; jp++) { 
	  j = i0-1+jp;
	  /* Careful here, we're filling in beta[cm->M][j][d] which is unbanded
	   * by adding beta[v][j+{0,1}][d+dv] to endsc[v], and we know there's a 
	   * band on v, so we can save time here as follows:
	   */
	  dv = StateDelta(cm->sttype[v]);
	  for (d = (dmin[v]-dv); d <= (dmax[v]-dv) && d <= jp; d++)
	    {
	      i = j-d+1;
	      switch (cm->sttype[v]) {
	      case MP_st: 
		if (j == j0 || d == jp) continue; /* boundary condition */
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		if ((sc = beta[v][j+1][d+2] + cm->endsc[v] + 
		     (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (d == jp) continue;	
		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		if ((sc = beta[v][j][d+1] + cm->endsc[v] + 
		     (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		if ((sc = beta[v][j+1][d+1] + cm->endsc[v] + 
		     (cm->el_selfsc * d) + escore) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if ((sc = beta[v][j][d] + cm->endsc[v] +
		     (cm->el_selfsc * d)) > beta[cm->M][j][d])
		  beta[cm->M][j][d] = sc;
		break;
	      case B_st:  
	      default: cm_Fail("bogus parent state %d\n", cm->sttype[v]);
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

#if 0
  /* SRE: this code is superfluous, yes??? */
  /* Deal with last step needed for local alignment 
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   * (EL = deck at M.)
   */
  if (cm->flags & CMH_LOCAL_END) {
    for (jp = W; jp > 0; jp--) { /* careful w/ boundary here */
      j = i0-1+jp;
      /* There is no band on the EL state */
      for (d = jp-1; d >= 0; d--)
	if ((sc = beta[cm->M][j][d+1]) > beta[cm->M][j][d])
	  beta[cm->M][j][d] = sc;
    }
  }
#endif

  /* If the caller doesn't want the matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation...)
   */
  if (ret_beta == NULL) {
    for (v = w1; v <= vend; v++) /* start at w1 - top of split set - not vroot */
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
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
  return;
 ERROR:
  cm_Fail("Memory allocation error.");
}


/* Function: vinside_qdb()
 *           EPN 05.19.05
 * *based on vinside(), only difference is bands are used : 
 * 
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
 *           useEL     - if TRUE, V problem ends at EL/i1/j1, not z/i1/j1
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
 *           allow_begin- TRUE to allow 0->b local alignment begin transitions. 
 *           ret_b     - best local begin state, or NULL if unwanted
 *           ret_bsc   - score for using ret_b, or NULL if unwanted                        
 *           dmin      - minimum d bound for each state v; [0..v..M-1]
 *           dmax      - maximum d bound for each state v; [0..v..M-1]
 * 
 * Returns:  score.
 */
static float
vinside_qdb(CM_t *cm, ESL_DSQ *dsq, int L, 
	int r, int z, int i0, int i1, int j1, int j0, int useEL,
	int do_full, float ***a, float ****ret_a,
	struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	char ****ret_shadow,
	int allow_begin, int *ret_b, float *ret_bsc, int *dmin, int *dmax)
{
  int     status;
  char  ***shadow;              /* the shadow matrix -- traceback ptrs -- memory is kept */
  int     v,i,j;
  int     w1,w2;		/* bounds of the split set */
  int     jp, ip;		/* j' and i' -- in the matrix coords */
  int    *touch;                /* keeps track of whether we can free a deck yet or not */
  int     y, yoffset;
  float   sc;			/* tmp variable holding a score */
  int      b;			/* best local begin state */
  float    bsc;			/* score for using the best local begin state */
  int     *imin;                /* minimum i bound for each state v; [0..w1-r] 
                                 * calculated using *dmin; offset from v, the
				 * band that corresponds to state v, is imin[v-r] */
  int     *imax;                /* maximum i bound for each state v; [0..w1-r] 
                                 * calculated using *dmax; offset from v, the
				 * band that corresponds to state v, is imax[v-r] */ 

  /*debugging block*/
  /*printf("***in vinside_qdb()****\n");
  printf("\tr  : %d\n", r);
  printf("\tz  : %d\n", z);
  printf("\ti0 : %d\n", i0);
  printf("\ti1 : %d\n", i1);
  printf("\tj1 : %d\n", j1);
  printf("\tj0 : %d\n", j0);
  */

  /* Allocations, initializations.
   * Remember to allocate for M+1 decks, in case we reuse this 
   * memorry for a local alignment voutside() calculation.
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  if (dpool == NULL) dpool = deckpool_create();
  if (a == NULL) {
    ESL_ALLOC(a, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) a[v] = NULL;
  }
				/* the whole split set w<=z<=y must be initialized */
  w1 = cm->nodemap[cm->ndidx[z]];
  w2 = cm->cfirst[w1]-1;

  /* Bands used ip 3 */
  /* Allocate imin and imax */

  imin = malloc(sizeof (int) * (w1-r+1));
  imax = malloc(sizeof (int) * (w1-r+1));

  for (v = w1; v <= w2; v++) { 
    if (! deckpool_pop(dpool, &(a[v]))) 
      a[v] = alloc_vji_deck(i0, i1, j1, j0);
    for (jp = 0; jp <= j0-j1; jp++) 
      for (ip = 0; ip <= i1-i0; ip++) 
	a[v][jp][ip] = IMPOSSIBLE;
  }

  if (ret_shadow != NULL) {
    ESL_ALLOC(shadow, sizeof(char **) * cm->M);
    for (v = 0; v < cm->M; v++) shadow[v] = NULL; 
  }

  /* Initialize the one non-IMPOSSIBLE cell as a boundary
   * condition.
   * If local alignment (useEL=1), we must connect z to EL;
   * we would init a[EL][0][i1-i0] = 0. But, we're not explicitly
   * keeping an EL deck, we're swallowing it into the recursion.
   * So, we unroll a chunk of the main recursion;
   * we have to laboriously figure out from the statetype z 
   * and our position where and what our initialization is.
   * Else, for global alignments, we simply connect to z,0,i1-i0.
   */

  ip = i1-i0;
  jp = 0;

  if (! useEL) 
    a[z][jp][ip] = 0.;
  else 
    {
      if (ret_shadow != NULL) 
	shadow[z] = alloc_vji_shadow_deck(i0,i1,j1,j0); 

      switch (cm->sttype[z]) {
      case D_st:
      case S_st:
	/*a[z][jp][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));
	if (ret_shadow != NULL) shadow[z][jp][ip] = USED_EL;
	break;
      case MP_st:
	if (i0 == i1 || j1 == j0) break;
	/*a[z][jp+1][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp+1][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));

	if (dsq[i1-1] < cm->abc->K && dsq[j1+1] < cm->abc->K)
	  a[z][jp+1][ip-1] += cm->esc[z][(int) (dsq[i1-1]*cm->abc->K+dsq[j1+1])];
	else
	  a[z][jp+1][ip-1] += DegeneratePairScore(cm->abc, cm->esc[z], dsq[i1-1], dsq[j1+1]);
	if (ret_shadow != NULL) shadow[z][jp+1][ip-1] = USED_EL;
	if (a[z][jp+1][ip-1] < IMPOSSIBLE) a[z][jp+1][ip-1] = IMPOSSIBLE;
	break;
      case ML_st:
      case IL_st:
	if (i0==i1) break;
	/*a[z][jp][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp][ip-1] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));

	if (dsq[i1-1] < cm->abc->K)
	  a[z][jp][ip-1] += cm->esc[z][(int) dsq[i1-1]];
	else
	  a[z][jp][ip-1] += esl_abc_FAvgScore(cm->abc, dsq[i1-1], cm->esc[z]);
	if (ret_shadow != NULL) shadow[z][jp][ip-1] = USED_EL;
	if (a[z][jp][ip-1] < IMPOSSIBLE) a[z][jp][ip-1] = IMPOSSIBLE;
	break;
      case MR_st:
      case IR_st:
	if (j1==j0) break;
	/*a[z][jp+1][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1 - StateDelta(cm->sttype[z])));*/
	a[z][jp+1][ip] = cm->endsc[z] + (cm->el_selfsc * ((jp+j1)-(ip+i0)+1));
	
	if (dsq[j1+1] < cm->abc->K)
	  a[z][jp+1][ip] += cm->esc[z][(int) dsq[j1+1]];
	else
	  a[z][jp+1][ip] += esl_abc_FAvgScore(cm->abc, dsq[j1+1], cm->esc[z]);
	if (ret_shadow != NULL) shadow[z][jp+1][ip] = USED_EL;
	if (a[z][jp+1][ip] < IMPOSSIBLE) a[z][jp+1][ip] = IMPOSSIBLE;
	break;
      }
    } /* done initializing the appropriate cell for useEL=TRUE */
  
  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0;   v < r;  v++) touch[v] = 0;
  for (v = r;   v <= w2; v++) touch[v] = cm->pnum[v]; /* note w2 not z: to bottom of split set */
  for (v = w2+1; v < cm->M; v++) touch[v] = 0;

  /* A special case. If vinside() is called on empty sequences,
   * we might do a begin transition right into z.
   */ 
  
  /* EPN 05.19.05 
     We are setting alpha cells in the following block, we should make
     sure they're within the bands */
  
  if (allow_begin && j0-j1 == 0 && i1-i0 == 0)
    {
      b   = z;
      bsc = a[z][0][0] + cm->beginsc[z];
      if (z == 0) { 
	a[0][0][0] = bsc;
	if (ret_shadow != NULL) shadow[0][0][0] = USED_LOCAL_BEGIN;
      }
    }
  
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

      /* Bands used ip 8 */
      /* First fill imin[v] and imax[v] */
      
      /* debugging block */
      /*
      if((dmin[v] > (j0-i0+1)) || (dmax[v] < (j1-i1+1)))
	{
	  printf("ERROR vinside_qdb() whole deck is outside bands\n");
	  printf("v : %d\n", v);
	  printf("dmin[v] : %d\n", dmin[v]);
	  printf("dmax[v] : %d\n", dmax[v]);
	  printf("i0 : %d\n", i0);
	  printf("i1 : %d\n", i1);
	  printf("j1 : %d\n", j1);
	  printf("j0 : %d\n", j0);
	}
      */
  
      imin[v-r] = j1-i0-dmax[v]+1;
      imax[v-r] = j1-i0-dmin[v]+1;

      /* Bands used ip 8 continued */
      /* Impose bands by setting all cells outside the bands to IMPOSSIBLE 
       * This is independent of state type so we do it outside
       * the following set of if then statements. 
       * Alternatively, it could be done within each of the following
       * if(cm->sttype[v] == *) statements - matter of style I suppose.
       */

      for (jp = 0; jp <= j0-j1; jp++) {
	for (ip = 0; ip < (imin[v-r]+jp) && ip<=(i1-i0); ip++)
	  {
	    a[v][jp][ip] = IMPOSSIBLE;
	  }
	if((imax[v-r]+jp) > (i1-i0)) ip = (i1-i0+1);
	else ip = imax[v-r]+jp+1;
	if(ip < 0) ip = 0;
	for (; ip <= (i1-i0); ip++) 
	  {
	    a[v][jp][ip] = IMPOSSIBLE;
	  }
      }      
      /* reassert our definition of a V problem */
      if (cm->sttype[v] == E_st || cm->sttype[v] == B_st || (cm->sttype[v] == S_st && v > r))
	cm_Fail("you told me you wouldn't ever do that again.");
      
      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= j0-j1; jp++) 
	    {
	      /* Bands used ip 9 */
	      /* old line :  for (ip = i1-i0; ip >= 0; ip--) { */
	      /* Use the imin[v-r] and imax[v-r] we have already set (see Bands used ip 3B) */
	      /* Remember 'state' indices in imin and imax are offset from v because imin and 
		 imax run [0..z-r], ==> dmin[v] corresponds to imin[v-r] */
	      if((imax[v-r]+jp) > (i1-i0)) ip = (i1-i0);
	      else ip = imax[v-r] + jp;
	      for(; ip >= imin[v-r]+jp && ip >= 0; ip--) {
		y = cm->cfirst[v];
		a[v][jp][ip]      = a[y][jp][ip] + cm->tsc[v][0];
		if (ret_shadow != NULL) shadow[v][jp][ip] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		   > a[v][jp][ip])) {
		a[v][jp][ip]      = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = a[y+yoffset][jp][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		    { 
		      a[v][jp][ip] = sc;
		      if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		    }
		if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;
	      }
	    }
	} else if (cm->sttype[v] == MP_st) {
	  
	  /* EPN following line redundant? are these cells already IMPOSSIBLE
	     due to band imposition */
	  
	  for (ip = i1-i0; ip >= 0; ip--) a[v][0][ip] = IMPOSSIBLE; /* boundary condition */
	  
	  for (jp = 1; jp <= j0-j1; jp++) { 
	    j = jp+j1;
	    a[v][jp][i1-i0] = IMPOSSIBLE; /* boundary condition */
	    /* Bands used ip 10 */
	    /* old line :  for (ip = i1-i0-1; ip >= 0; ip--) { */
	    /* Use the imin[v-w1] and imax[v-w1] we have already set (see Bands used ip 3B) */
	    /* Remember 'state' indices in imin and imax are offset from v because imin and 
	       imax run [0..z-r], ==> dmin[v] corresponds to imin[v-r] */
	    if((imax[v-r]+jp) > (i1-i0-1)) ip = (i1-i0-1);
	    else ip = imax[v-r] + jp;
	    for(; ip >= imin[v-r]+jp && ip >= 0; ip--) {
	      i = ip+i0;
	      y = cm->cfirst[v];
	      a[v][jp][ip] = a[y][jp-1][ip+1] + cm->tsc[v][0];
	      if (ret_shadow != NULL) shadow[v][jp][ip] = (char) 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		  > a[v][jp][ip])) {
		a[v][jp][ip]      = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc; 
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
	      if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		a[v][jp][ip] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
	      else
		a[v][jp][ip] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
	  
	  for (jp = 0; jp <= j0-j1; jp++) { 
	    a[v][jp][i1-i0] = IMPOSSIBLE; /* boundary condition */
	    /* Bands used ip 11 */
	    /* old line :  for (ip = i1-i0-1; ip >= 0; ip--) { */
	    /* Use the imin[v-w1] and imax[v-w1] we have already set (see Bands used ip 3B) */
	    /* Remember 'state' indices in imin and imax are offset from v because imin and 
	       imax run [0..z-r], ==> dmin[v] corresponds to imin[v-r] */
	    if((imax[v-r]+jp) > (i1-i0-1)) ip = (i1-i0-1);
	    else ip = imax[v-r] + jp;
	    for(; ip >= imin[v-r]+jp && ip >= 0; ip--) {
	      i = ip+i0;
	      y = cm->cfirst[v];
	      a[v][jp][ip] = a[y][jp][ip+1] + cm->tsc[v][0];
	      if (ret_shadow != NULL) shadow[v][jp][ip] = 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		  > a[v][jp][ip])) {
		a[v][jp][ip]      = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		/*printf("set a[%d][%d][%d] to %f\n", v, jp, ip, sc);*/
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp][ip+1] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc; 
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
	      
	      if (dsq[i] < cm->abc->K)
		a[v][jp][ip] += cm->esc[v][dsq[i]];
	      else
		a[v][jp][ip] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
	  /* EPN following line redundant? are these cells already IMPOSSIBLE
	     due to band imposition */
	  for (ip = i1-i0; ip >= 0; ip--) a[v][0][ip] = IMPOSSIBLE; /* boundary condition */
	  
	  for (jp = 1; jp <= j0-j1; jp++) { 
	    j = jp+j1;
	    /* Bands used ip 12 */
	    /* old line :  for (ip = i1-i0; ip >= 0; ip--) { */
	    /* Use the imin[v-w1] and imax[v] we have already set (see Bands used ip 3B) */
	    /* Remember 'state' indices in imin and imax are offset from v because imin and 
	       imax run [0..z-r], ==> dmin[v] corresponds to imin[v-r] */
	    /*05.20 for (ip = imax[v-r]; ip >= imin[v-r]; ip--) {		*/
	    if((imax[v-r]+jp) > (i1-i0)) ip = (i1-i0);
	    else ip = imax[v-r] + jp;
	    for(; ip >= imin[v-r]+jp && ip >= 0; ip--) {
	      y = cm->cfirst[v];
	      a[v][jp][ip]      = a[y][jp-1][ip] + cm->tsc[v][0];
	      if (ret_shadow != NULL) shadow[v][jp][ip] = 0;
	      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v]) && 
		  ((cm->endsc[v] + (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v]))))
		  > a[v][jp][ip])) {
		a[v][jp][ip] = cm->endsc[v] + 
		  (cm->el_selfsc * (((jp+j1)-(ip+i0)+1) - StateDelta(cm->sttype[v])));
		if (ret_shadow != NULL) shadow[v][jp][ip] = USED_EL;
	      }
	      for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		if ((sc = a[y+yoffset][jp-1][ip] + cm->tsc[v][yoffset]) >  a[v][jp][ip])
		  { 
		    a[v][jp][ip] = sc; 
		    if (ret_shadow != NULL) shadow[v][jp][ip] = (char) yoffset; 
		  }
	      
	      if (dsq[j] < cm->abc->K)
		a[v][jp][ip] += cm->esc[v][dsq[j]];
	      else
		a[v][jp][ip] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
	      if (a[v][jp][ip] < IMPOSSIBLE) a[v][jp][ip] = IMPOSSIBLE;  
	    }
	  }
	} /* finished calculating deck v */
      
      /* Check for local begin getting us to the root.
       */
      if (allow_begin && a[v][j0-j1][0] + cm->beginsc[v] > bsc) 
	{
	  b   = v;
	  bsc = a[v][j0-j1][0] + cm->beginsc[v];
	}

      /* Check whether we need to store the local begin score
       * for a possible traceback.
       */
      if (allow_begin && v == 0 && bsc > a[0][j0-j1][0]) 
	{
	  a[0][j0-j1][0] = bsc;
	  if (ret_shadow != NULL) shadow[v][j0-j1][0] = USED_LOCAL_BEGIN;
	}


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
  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */
  
  
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
  free(imax);
  free(imin);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}


/* Function: voutside_qdb()
 *           EPN 05.19.05
 * *based on voutside(), only difference is bands are used : 
 *
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
 *           analogous to inside() and vinside(); see their prefaces
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
 *           useEL     - if TRUE, worry about local alignment.
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
 *           dmin      - minimum d bound for each state v; [0..v..M-1]
 *           dmax      - maximum d bound for each state v; [0..v..M-1]
 * 
 */
static void
voutside_qdb(CM_t *cm, ESL_DSQ *dsq, int L, 
	   int r, int z, int i0, int i1, int j1, int j0, int useEL,
	   int do_full, float ***beta, float ****ret_beta,
	   struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
	   int *dmin, int *dmax)
{
  int      status;
  int      v,y;			/* indices for states */
  int      i,j;			/* indices in sequence dimensions */
  int      ip, jp;		/* transformed sequence indices */
  float    sc;			/* a temporary variable holding a score */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */
  int      voffset;		/* index of v in t_v(y) transition scores */
  int     *imin;                /* minimum i bound for each state v; [0..r-z] 
                                 * calculated using *dmin; offset from v, the
				 * band that corresponds to state v, is imin[v-r] */
  int     *imax;                /* maximum i bound for each state v; [0..r-z] 
                                 * calculated using *dmax; offset from v, the
				 * band that corresponds to state v, is imax[v-r] */
  int      dv;                  /* state delta */				   

  /* Allocations and initializations
   */
  			/* if caller didn't give us a deck pool, make one */
  if (dpool == NULL) dpool = deckpool_create();

  /* If caller didn't give us a matrix, make one.
   * Remember to allow for deck M, the EL deck, for local alignments.
   */
  if (beta == NULL) {
    ESL_ALLOC(beta, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) beta[v] = NULL;
  }
  /* Initialize the root deck. This probably isn't the most efficient way to do it.
   */
  if (! deckpool_pop(dpool, &(beta[r])))
    beta[r] = alloc_vji_deck(i0,i1,j1,j0);
  for (jp = 0; jp <= j0-j1; jp++) {
    for (ip = 0; ip <= i1-i0; ip++)
      beta[r][jp][ip] = IMPOSSIBLE;
  }
  /* Bands used ip 15 */
  /* We want to make sure that imin[0] <= 0; but we don't have imin[0] */
  /* First calculate imin[0], then assert its less than 0, not sure
     if this is necessary, imin[0] == 0 may be guaranteed, I'll use
     the assert here to be safe*/
  /* Note imin[0] corresponds to state r */

  imin = malloc(sizeof (int) * (z-r+1));
  imax = malloc(sizeof (int) * (z-r+1));

  /* debugging block */
  /*
  if((dmin[r] > (j0-i0)) || (dmax[r] < (j1-i1)))
    {
      printf("ERROR voutside_qdb()\n");
      printf("v : %d\n", r);
      printf("dmin[v] : %d\n", dmin[r]);
      printf("dmax[v] : %d\n", dmax[r]);
      printf("i0 : %d\n", i0);
      printf("i1 : %d\n", i1);
      printf("j1 : %d\n", j1);
      printf("j0 : %d\n", j0);
    }
  */

  assert(dmin[r] <= (j0-i0)+1); 
  assert(dmax[r] >= (j1-i1)+1); 

  imin[0] = j1-i0-dmax[r]+1;
  imax[0] = j1-i0-dmin[r]+1;

  assert(imin[0] <= 0);

  beta[r][j0-j1][0] = 0;		
  
  /* Initialize the EL deck, if we're in local mode w.r.t. ends.
   * Deal with the special initialization case of the root state r
   * immediately transitioning to EL, if we're supposed to use EL.
   */
  
  if (useEL && cm->flags & CMH_LOCAL_END) {
    if (! deckpool_pop(dpool, &(beta[cm->M])))
      beta[cm->M] = alloc_vji_deck(i0,i1,j1,j0);
    for (jp = 0; jp <= j0-j1; jp++) {
      for (ip = 0; ip <= i1-i0; ip++)
	beta[cm->M][jp][ip] = IMPOSSIBLE;
    }
  }
  if (useEL && NOT_IMPOSSIBLE(cm->endsc[r])) {
    switch(cm->sttype[r]) {
    case MP_st:
      if (i0 == i1 || j1 == j0) break;
      if (dsq[i0] < cm->abc->K && dsq[j0] < cm->abc->K)
	escore = cm->esc[r][(int) (dsq[i0]*cm->abc->K+dsq[j0])];
      else
	escore = DegeneratePairScore(cm->abc, cm->esc[r], dsq[i0], dsq[j0]);
      beta[cm->M][j0-j1-1][1] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0-1)-(i0+1)+1)) + escore;
      break;
    case ML_st:
    case IL_st:
      if (i0 == i1) break;
      if (dsq[i0] < cm->abc->K) 
	escore = cm->esc[r][(int) dsq[i0]];
      else
	escore = esl_abc_FAvgScore(cm->abc, dsq[i0], cm->esc[r]);
      beta[cm->M][j0-j1][1] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0)-(i0+1)+1)) + escore;
      break;
    case MR_st:
    case IR_st:
      if (j0==j1) break;
      if (dsq[j0] < cm->abc->K) 
	escore = cm->esc[r][(int) dsq[j0]];
      else
	escore = esl_abc_FAvgScore(cm->abc, dsq[j0], cm->esc[r]);
      beta[cm->M][j0-j1-1][0] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0-1)-(i0)+1)) + escore;
      break;
    case S_st:
    case D_st:
      beta[cm->M][j0-j1][0] = cm->endsc[r] + 
	(cm->el_selfsc * ((j0)-(i0)+1));
      break;
    default:  cm_Fail("bogus parent state %d\n", cm->sttype[r]);
    }
  }
      
  /* Initialize the "touch" array, used for figuring out
   * when a deck is no longer touched, so it can be free'd.
   */
  ESL_ALLOC(touch, sizeof(int) * cm->M);
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
      /* Bands used ip 16 */
      /* Fill imin[v-r+1] and imax[v-r+1] as we go, one of many ways to handle imin and imax */
      /* Remember 'state' indices in imin and imax are offset from v because imin and 
	 imax run [0..z-r+1], ==> dmin[v] corresponds to imin[v-r] */

      imin[v-r] = j1-i0-dmax[v]+1;
      imax[v-r] = j1-i0-dmin[v]+1;

      /* An awkward situation here.  If dmin[v] > i1, imin[v-r] will be 0
	 however, we don't want to query ANY cells (in other words
	 none of the following for(ip*) loops should ever be entered)
	 because in this case the whole vji deck is outside the bands, so
	 the bestsc we want is IMPOSSIBLE (which was set before the
	 for (v = w; v <= y; v++) loop).  There is probably a better
         way to do this but I'll explicitly check for this situation.
         Note - it's okay if dmax < i0 (which also means the entire
         deck is outside the bands) because this will make the
         for(ip*) loops always evaluate to false because imin[v-r] will
         be 0 and imax[v-r] will be < 0.*/
      /* This situation is recapitulated in v_splitter_qdb() */

      /* unnecssary 05.22
	 05.20 code : if(dmin[v] > i1) imin[v-r] = imax[v-r]+1;  */
	 /* now the for(ip) loops
						    will never be entered
						    (see above comments) */

      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_vji_deck(i0,i1,j1,j0);

      /* Init the whole deck to IMPOSSIBLE.
       */
      for (jp = j0-j1; jp >= 0; jp--) 
	for (ip = 0; ip <= i1-i0; ip++) 
	  beta[v][jp][ip] = IMPOSSIBLE;

      /* We've set the whole matrix to impossible, everything outside bands must be impossible */
	 
      /* If we can get into deck v by a local begin transition, do an init
       * with that.
       */
      if (r == 0 && i0 == 1 && j0 == L && (cm->flags & CMH_LOCAL_BEGIN))
	{
	  if (cm->beginsc[v] > beta[v][j0-j1][0]) 
	    beta[v][j0-j1][0] = cm->beginsc[v];
	}

      /* main recursion:
       */
      for (jp = j0-j1; jp >= 0; jp--) {
	j = jp+j1;
	/* Bands used ip 17 */
	/* old line :	for (ip = 0; ip <= i1-i0; ip++) */
	/* Remember 'state' indices in imin and imax are offset from v because imin and 
	   imax run [0..z-r+1], ==> dmin[v] corresponds to imin[v-r] */
	/* 05.20 for (ip = imin[v-r]; ip <= imax[v-r]; ip++) */

	if((imin[v-r]+jp) < 0) ip = 0;
        else ip = imin[v-r]+jp;
        for(; ip <= imax[v-r] + jp && ip <= (i1-i0); ip++)
	  {
	    i = ip+i0;

	    for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
	      if (y < r) continue; /* deal with split sets */
	      voffset = v - cm->cfirst[y]; /* gotta calculate the transition score index for t_y(v) */

	      switch(cm->sttype[y]) {
	      case MP_st: 
		if (j == j0 || i == i0) continue; /* boundary condition */

		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[y][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[y], dsq[i-1], dsq[j+1]);
		
		if ((sc = beta[y][jp+1][ip-1]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;

	      case ML_st:
	      case IL_st: 
		if (i == i0) continue;	/* boundary condition */

		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[y]);
		  
		if ((sc = beta[y][jp][ip-1]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;
		  
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		  
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[y][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[y]);

		if ((sc = beta[y][jp+1][ip]+cm->tsc[y][voffset]+escore) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;
		  
	      case S_st:
	      case E_st:
	      case D_st:
		if ((sc = beta[y][jp][ip] + cm->tsc[y][voffset]) > beta[v][jp][ip])
		  beta[v][jp][ip] = sc;
		break;

	      default: cm_Fail("bogus parent state %d\n", cm->sttype[y]);
	      }/* end switch over states*/
	    } /* ends for loop over parent states. we now know beta[v][j][d] for this d */
	    if (beta[v][jp][ip] < IMPOSSIBLE) beta[v][jp][ip] = IMPOSSIBLE;

	  } /* ends loop over ip. We know all beta[v][jp][ip] in this row jp */
	
      }/* end loop over jp. We know the beta's for the whole deck.*/
      
      /* Deal with local alignment
       * transitions v->EL, if we're doing local alignment and there's a 
       * possible transition.
       */
      if (useEL && NOT_IMPOSSIBLE(cm->endsc[v])) {
	for (jp = j0-j1; jp >= 0; jp--) {
	  j = jp+j1;
	  /* Careful here, we're filling in beta[cm->M][jp][ip] which is unbanded
	   * by adding beta[v][jp+{0,1}][ip-{0,1}] to endsc[v], and we know there's a 
	   * i band on v (imin[v-r]..imax[v-r], so we can save time here as follows:
	   */
	  dv = StateDelta(cm->sttype[v]);
	  if((imin[v-r]+jp+dv) < 0) ip = 0;
	  else ip = imin[v-r]+jp+dv;
	  for(; (ip<=imax[v-r]+jp+dv) && ip <= (i1-i0); ip++)
	    {
	      i = ip+i0;
	      switch (cm->sttype[v]) {
	      case MP_st:
		if (j == j0 || i == i0) continue; /* boundary condition */
		if (dsq[i-1] < cm->abc->K && dsq[j+1] < cm->abc->K)
		  escore = cm->esc[v][(int) (dsq[i-1]*cm->abc->K+dsq[j+1])];
		else
		  escore = DegeneratePairScore(cm->abc, cm->esc[v], dsq[i-1], dsq[j+1]);
		if ((sc = beta[v][jp+1][ip-1] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1))
		     + escore) > beta[cm->M][jp][ip])
		  beta[cm->M][jp][ip] = sc;
		break;
	      case ML_st:
	      case IL_st:
		if (i == i0) continue;
		if (dsq[i-1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[i-1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[i-1], cm->esc[v]);
		if ((sc = beta[v][jp][ip-1] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1))
		     + escore) > beta[cm->M][jp][ip])
		  beta[cm->M][jp][ip] = sc;
		break;
	      case MR_st:
	      case IR_st:
		if (j == j0) continue;
		if (dsq[j+1] < cm->abc->K) 
		  escore = cm->esc[v][(int) dsq[j+1]];
		else
		  escore = esl_abc_FAvgScore(cm->abc, dsq[j+1], cm->esc[v]);
		if ((sc = beta[v][jp+1][ip] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1))
		     + escore) > beta[cm->M][jp][ip])
		  beta[cm->M][jp][ip] = sc;
		break;
	      case S_st:
	      case D_st:
	      case E_st:
		if ((sc = beta[v][jp][ip] + cm->endsc[v] + 
		     (cm->el_selfsc * (j-i+1)))
		     > beta[cm->M][jp][ip])
		  beta[cm->M][jp][ip] = sc;
		break;
	      default:  cm_Fail("bogus parent state %d\n", cm->sttype[y]);
	      } /* end switch over parent v state type */
	    } /* end loop over ip */
	} /* end loop over jp */
      }
      
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
  
#if 0 
  /* superfluous code, I think...*/
  /* Deal with the last step needed for local alignment
   * w.r.t. ends: left-emitting, zero-scoring EL->EL transitions.
   */
  if (useEL && cm->flags & CMH_LOCAL_END) {
    for (jp = j0-j1; jp >= 0; jp--) 
      {
	/* Bands used ip 19 */
	/* Actually the bands are not used here, because there are no bands for 
	   state cm->M.  I'll just leave the unbanded code alone here.  Not sure
	   how to think about bands in terms of local alignment??? */
	for (ip = 1; ip <= i1-i0; ip++) /* careful w/boundary here */
	  if ((sc = beta[cm->M][jp][ip-1]) > beta[cm->M][jp][ip]) 
	    beta[cm->M][jp][ip] = sc;
      }
  }
#endif
  
  /* If the caller doesn't want the matrix, free it.
   * (though it would be *stupid* for the caller not to want the
   * matrix in the current implementation!)
   */
  if (ret_beta == NULL) {
    for (v = r; v <= z; v++)
      if (beta[v] != NULL) { deckpool_push(dpool, beta[v]); beta[v] = NULL; }
    if (cm->flags & CMH_LOCAL_END) {
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
    while (deckpool_pop(dpool, &a)) 
      free_vji_deck(a,j1,j0);
    deckpool_free(dpool);
  } else *ret_dpool = dpool;
  
  free(touch);
  free(imax);
  free(imin);
  return;
 ERROR:
  cm_Fail("Memory allocation error.");
}


/* For the Full CYK memory efficient banded implementation we need 
 *  banded versions of some of the memory management routines 
 *
 * The D&C banded implementation is not memory efficient, in that
 * it requires the same amount of memory as the non-banded D&C implementation.
 * This means that we still allocate the same memory as we would without bands, 
 * we just set all cells of alpha or beta that are outside of the bands to 
 * IMPOSSIBLE.  Because of this we should be able to use the same memory management 
 * routines as the non-banded implementation.
 *
 * Therefore we can use the D&C memory routines for banded D&C.
 */

/*################################################################*/
/* EPN *_banded_vjd_* 
   adapted from *_vjd_* from SRE*/

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
alloc_banded_vjd_deck(int L, int i, int j, int min, int max)
{
  int     status;
  float **a;
  int     jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/

  /*printf("in alloc banded vjd deck, L : %d, i : %d, j : %d, min : %d, max : %d\n", L, i, j, min, max);*/

  ESL_DPRINTF3(("#DEBUG: alloc_vjd_deck : %.4f\n", size_vjd_deck(L,i,j)));
  ESL_ALLOC(a, sizeof(float *) * (L+1)); /* always alloc 0..L rows, some of which are NULL */
  for (jp = 0;   jp < i-1;    jp++) a[jp]     = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp]     = NULL;
  for (jp = 0; jp <= j-i+1; jp++) 
    {
      if(jp > max)
	bw = max - min + 1;
      else
	bw = jp - (min) + 1;

      if(bw > 0)
	{
	  /*printf("\tallocated a[%d]\n", jp+i-1);*/
	  ESL_ALLOC(a[jp+i-1], sizeof(float) * bw);
	}
      else
	{
	  a[jp+i-1] = NULL;
	  /*printf("\tdid not allocate a[%d]\n", jp+i-1);*/
	}
    }
  return a;

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}

char **
alloc_banded_vjd_yshadow_deck(int L, int i, int j, int min, int max)
{
  int    status;
  char **a;
  int    jp;
  int    bw; /* width of band, depends on jp, so we need to calculate
	        this inside the jp loop*/

  ESL_ALLOC(a, sizeof(char *) * (L+1)); /* always alloc 0..L rows, same as alloc_deck */
  for (jp = 0;   jp < i-1;    jp++) a[jp] = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp] = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) 
    {
      if(jp > max)
	bw = max - min + 1;
      else
	bw = jp - min + 1;
      if(bw > 0)
	{
	  ESL_ALLOC(a[jp+i-1], sizeof(char) * (bw));
	}
      else a[jp+i-1] = NULL;
    }
  return a;

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}
int **
alloc_banded_vjd_kshadow_deck(int L, int i, int j, int min, int max)
{
  int   status;
  int **a;
  int   jp;
  int     bw; /* width of band, depends on jp, so we need to calculate
	         this inside the jp loop*/

  ESL_ALLOC(a, sizeof(int *) * (L+1)); /* always alloc 0..L rows, same as alloc_deck */
  for (jp = 0;   jp <  i-1;   jp++) a[jp] = NULL;
  for (jp = j+1; jp <= L;     jp++) a[jp] = NULL;
  for (jp = 0;   jp <= j-i+1; jp++) 
    {
      if(jp > max) bw = max - min + 1;
      else bw = jp - min + 1;
      if(bw > 0)
	{
	  ESL_ALLOC(a[jp+i-1], sizeof(int) * bw);
	}
      else a[jp+i-1] = NULL;
    }
  
  return a;

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
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
/* EPN 05.09.05
   debug_print_shadow()
 * Function: debug_print_shadow
 *
 * Purpose:  Print shadow matrix 
 */

void
debug_print_shadow(void ***shadow, CM_t *cm, int L)
{
  int v, j, d;
  int yoffset;
  
  printf("\nPrinting alpha matrix :\n");
  printf("************************************\n");
  for(v = 0; v < cm->M; v++)
    {
      printf("====================================\n");
      for(j = 0; j <= L; j++)
	{
	  printf("------------------------------------\n");
	  for(d = 0; d <= j; d++)
	    {
	      if(cm->sttype[v] == E_st)
		{
		  printf("END state\n");
		}
	      else
		{
		  if(cm->sttype[v] == B_st)
		    {
		      yoffset = ((int **) shadow[v])[j][d];
		      printf("INT  shadow[%2d][%2d][%2d] : %d\n", v, j, d, yoffset);
		    }
		  else
		    {
		      yoffset = ((int **) shadow[v])[j][d];
		      printf("CHAR shadow[%2d][%2d][%2d] : %d\n", v, j, d, yoffset);
		    }
		}
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

void
debug_print_shadow_banded(void ***shadow, CM_t *cm, int L, int *dmin, int *dmax)
{
  int v, j, d, vdp;
  int yoffset;

  printf("\nPrinting banded shadow matrix :\n");
  printf("************************************\n");
  for(v = 0; v < cm->M; v++)
    {
      printf("====================================\n");
      for(j = 0; j <= L; j++)
	{
	  printf("------------------------------------\n");
	  /* there may be a problem with using j and not jp */
	  for (d = dmin[v]; d <= dmax[v] && d <= j; d++) 
	    {
	      vdp = d - dmin[v]; /* d index for state v in alpha w/mem eff bands */
	      if(cm->sttype[v] == E_st)
		{
		  printf("END state\n");
		}
	      else
		{
		  if(cm->sttype[v] == B_st)
		    {
		      yoffset = ((int **) shadow[v])[j][vdp];
		      printf("INT  shadow[%2d][%2d][%2d] : %d | d is %d\n", v, j, vdp, yoffset, d);
		    }
		  else
		    {
		      yoffset = ((int **) shadow[v])[j][vdp];
		      printf("CHAR shadow[%2d][%2d][%2d] : %d | d is %d\n", v, j, vdp, yoffset, d);
		    }
		}
	    }
	}
    }
  printf("****************\n\n");
}

/* EPN 05.16.05
   debug_print_shadow_banded_deck()
 * Function: debug_print_shadow_banded_deck
 *
 * Purpose:  Print banded shadow matrix deck
 */

void
debug_print_shadow_banded_deck(int v, void ***shadow, CM_t *cm, int L, int *dmin, int *dmax)
{
  int j, d, vdp;
  int yoffset;

  printf("\nPrinting banded shadow matrix deck for v : %d:\n", v);
  printf("====================================\n");
  for(j = 0; j <= L; j++)
    {
      printf("------------------------------------\n");
      /* there may be a problem with using j and not jp*/
      for (d = dmin[v]; d <= dmax[v] && d <= j; d++) 
	{
	  vdp = d - dmin[v]; /* d index for state v in alpha w/mem eff bands */

	  if(cm->sttype[v] == E_st)
	    {
	      printf("END state\n");
	    }
	  else
	    {
	      yoffset = ((char **) shadow[v])[j][vdp];
	      printf("shadow_banded[%2d][%2d][%2d] : %d| d is %d\n", v, j, vdp, yoffset, d);
	    }
	}
    }
}



/* EPN 05.09.05
   debug_print_alpha_banded()
 * Function: debug_print_alpha_banded
 *
 * Purpose:  Print alpha matrix 
 */
void
debug_print_alpha_banded(float ***alpha, CM_t *cm, int L, int *dmin, int *dmax)
{
  int v, j, d, vdp, max_v;

  printf("\nPrinting banded alpha matrix :\n");
  printf("************************************\n");
  max_v = cm->M-1;
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      max_v = cm->M;
    }
  for(v = 0; v <= max_v; v++)
    {
      printf("====================================\n");
      for(j = 0; j <= L; j++)
	{
	  printf("------------------------------------\n");
	  for (d = dmin[v]; d <= dmax[v] && d <= j; d++) 
	    {
	      vdp = d - dmin[v]; /* d index for state v in alpha w/mem eff bands */
	      printf("alpha[%2d][%2d][%2d] : %6.2f | d is %d\n", v, j, vdp, alpha[v][j][vdp], d);
	    }
	}
    }
  printf("****************\n\n");
}

/* EPN 05.09.05
   debug_print_alpha()
 * Function: debug_print_alpha
 *
 * Purpose:  Print alpha matrix 
 */

void
debug_print_alpha(float ***alpha, CM_t *cm, int L)
{
  int v, j, d, max_v;

  printf("\nPrinting alpha matrix :\n");
  printf("************************************\n");
  max_v = cm->M-1;
  if(cm->flags & CMH_LOCAL_BEGIN)
    {
      max_v = cm->M;
    }
  for(v = 0; v <= max_v; v++)
    {
      printf("====================================\n");
      for(j = 0; j <= L; j++)
	{
	  printf("------------------------------------\n");
	  for(d = 0; d <= j; d++)
	    {
	      printf("alpha[%2d][%2d][%2d] : %6.2f\n", v, j, d, alpha[v][j][d]);
	    }
	}
    }
  printf("****************\n\n");
}


/* EPN Memory efficient banded functions */
/* Function: inside_qdb_me()
 *
 * Based on inside(), only difference is bands are used : 
 * further the bands are used in a memory-efficient way
 * Another big difference is that we can't employ the deck
 * reuse strategy because the size of each deck depends
 * on the band for that state, so each deck can be different.
 *
 * Comments below are from inside(): 
 * 
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
 *           dmin      - minimum d bound for each state v; [0..v..M-1]
 *           dmax      - maximum d bound for each state v; [0..v..M-1]
 *                       
 * Returns: Score of the optimal alignment.  
 */
static float 
inside_qdb_me(CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
	    float ***alpha, float ****ret_alpha, 
	    void ****ret_shadow, 
	    int allow_begin, int *ret_b, float *ret_bsc,
	    int *dmin, int *dmax)
{
  int      status;
  float  **end;         /* we re-use the end deck. */
  //int      nends;     /* counter that tracks when we can release end deck to the pool, not needed in me version*/
  int     *touch;       /* keeps track of how many higher decks still need this deck */
  int      v,y,z;	/* indices for states  */
  int      j,d,i;	/* indices in sequence dimensions */
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
  int      kp;             /* k' - what k should be, now that we're banded */
  int      Wp;             /* W also changes depending on state */

  /* Allocations and initializations
   */
  b   = -1;
  bsc = IMPOSSIBLE;
  W   = j0-i0+1;		/* the length of the subsequence -- used in many loops  */
				/* if caller didn't give us a deck pool, make one */
  end = alloc_vjd_deck(L, i0, j0);
  //nends = CMSubtreeCountStatetype(cm, vroot, E_st);
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
    ESL_ALLOC(alpha, sizeof(float **) * (cm->M+1));
    for (v = 0; v <= cm->M; v++) alpha[v] = NULL;
  }

  ESL_ALLOC(touch, (sizeof(int) * cm->M));
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
      alpha[v] = alloc_banded_vjd_deck(L, i0, j0, dmin[v], dmax[v]);
      
      if (ret_shadow != NULL) {
	if (cm->sttype[v] == B_st) {
	  kshad     = alloc_banded_vjd_kshadow_deck(L, i0, j0, dmin[v], dmax[v]);
	  shadow[v] = (void **) kshad;
	} else {
	  yshad     = alloc_banded_vjd_yshadow_deck(L, i0, j0, dmin[v], dmax[v]);
	  shadow[v] = (void **) yshad;
	}
      }

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		dp_v = d - dmin[v];  /* d index for state v in alpha w/mem eff bands */

		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][dp_v]  = USED_EL; 
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - dmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* check to make sure the cell we're about to query is within the
		       bands for state y; this might be more complex than necessary */
		    if((dp_y >= 0) && ((dp_y < (jp - (dmin[y+yoffset]) + 1))
				       && (dp_y < (dmax[y+yoffset] - dmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j][dp_y] + cm->tsc[v][yoffset]) >  alpha[v][j][dp_v]) {
			  alpha[v][j][dp_v] = sc; 
			  if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			}
		      }
		  }
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == B_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		z = cm->cnum[v];

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
		dp_v = d - dmin[v];  /* d index for state v in alpha w/mem eff bands */
		dp_y = d - dmin[y];  /* d index for state y in alpha w/mem eff bands */

		/* First make sure we have any valid kp, we know from inequality (2)
		   that kp <= d-dmin[y]-dmin[z] so if this is < 0 then no kp
		   is valid (see notes above) */

		if((d-dmin[y]-dmin[z]) >= 0)
		{
		  if(jp < dmax[y]) kp = d-dmin[z]-jp;
		  else kp = d-dmin[z]-dmax[y];
		  if(kp < 0) kp = 0;
		  if(kp <= dmax[z] - dmin[z]) /* make sure its valid in deck alpha[z] */
		    {
		      alpha[v][j][dp_v] = alpha[y][j-dmin[z]-kp][d-dmin[y]-dmin[z]-kp] 
			+ alpha[z][j][kp];
		      if (ret_shadow != NULL) kshad[j][dp_v] = kp;
		      for (kp = kp+1; kp <= (d-dmin[y]-dmin[z]) && kp <= (dmax[z]-dmin[z]);
			   kp++)
			{
			  /* the following if statement ensures that the alpha cell for 
			     state y and the cell for state z that we are about to query 
			     is in fact within the bands for state y and state z respectively*/
			  if ((sc = alpha[y][j-dmin[z]-kp][d-dmin[y]-dmin[z]-kp] 
			       + alpha[z][j][kp]) > alpha[v][j][dp_v]) 
			    {
			      alpha[v][j][dp_v] = sc;
			      if (ret_shadow != NULL) kshad[j][dp_v] = kp;
			    }
			}
		    }
		}
		else alpha[v][j][dp_v] = IMPOSSIBLE;
		/*else cm_Fail("cell in alpha matrix was not filled in due to bands.\n");*/
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    /* We assume dmin[v] >= 2 (it has to be) */
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		dp_v = d - dmin[v]; /* d index for state v in alpha w/mem eff bands */
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if(ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - dmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-2) >= 0) && (((dp_y-2) < (jp - (dmin[y+yoffset]) + 1))
					   && ((dp_y-2) < (dmax[y+yoffset] - dmin[y+yoffset] + 1))))
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
		  alpha[v][j][dp_v] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		else
		  alpha[v][j][dp_v] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 7 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;

	    /* we assume dmin[v] >= 1, it has to be */
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		dp_v = d - dmin[v]; /* d index for state v in alpha w/mem eff bands */
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - dmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-1) >= 0) && (((dp_y-1) < (jp - (dmin[y+yoffset]) + 1))
				      && ((dp_y-1) < (dmax[y+yoffset] - dmin[y+yoffset] + 1))))
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
		  alpha[v][j][dp_v] += cm->esc[v][dsq[i]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 9 end block */
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (jp = 0; jp <= W; jp++) {
	    j = i0-1+jp;
	    for (d = dmin[v]; d <= dmax[v] && d <= jp; d++)
	      {
		y = cm->cfirst[v];
		dp_v = d - dmin[v]; /* d index for state v in alpha w/mem eff bands */
		alpha[v][j][dp_v] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		/* treat EL as emitting only on self transition */
		if (ret_shadow != NULL) yshad[j][dp_v] = USED_EL;
		for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) 
		  {
		    dp_y = d - dmin[y+yoffset];  /* d index for state (y+yoffset) 
						   in alpha w/mem eff bands */
		    /* the following if statement ensures that the alpha cell for 
		       state y that we are about to query is in fact within the
		       bands for state y */
		    if(((dp_y-1) >= 0) && (((dp_y-1) < (jp - (dmin[y+yoffset]) + 1))
				      && ((dp_y-1) < (dmax[y+yoffset] - dmin[y+yoffset] + 1))))
		      {
			if ((sc = alpha[y+yoffset][j-1][dp_y-1] + cm->tsc[v][yoffset]) > alpha[v][j][dp_v])
			  {
			    alpha[v][j][dp_v] = sc;
			    if (ret_shadow != NULL) yshad[j][dp_v] = yoffset;
			  }
		      }
		  }
		if (dsq[j] < cm->abc->K)
		  alpha[v][j][dp_v] += cm->esc[v][dsq[j]];
		else
		  alpha[v][j][dp_v] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		
		if (alpha[v][j][dp_v] < IMPOSSIBLE) alpha[v][j][dp_v] = IMPOSSIBLE;
		/* CYK Full ME Bands used 11 end block */
	      }
	  }
	}				/* finished calculating deck v. */
      
      /* The following loops originally access alpha[v][j0][W] but the index W will be
	 in different positions due to the bands */

      Wp = W - dmin[v];
      /* We need to make sure that Wp is within the bands */
      if(Wp >= 0 && Wp <= (dmax[v] - dmin[v]))
	{
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
      /* In the non-banded code, we used the deck reuse strategy, however, here
	 we can't do that, because for each state, the bands are different, so 
	 we can't use old decks, but rather must allocate a new one, and free
	 the old one. */

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
		    if (cm->sttype[y] == E_st) { 
		      //nends--; 
		      /* Original code : if (nends == 0) { deckpool_push(dpool, end); end = NULL;} */
		      /* ME code deletes the previous line, we don't mess with end, because
			 it is used later */
		    } else 
		      free_vjd_deck(alpha[y], i0, j0);
		    alpha[y] = NULL;
		  }
	      }
	  }
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
  Wp = W - dmin[vroot];
  sc       = alpha[vroot][j0][Wp];

  if (ret_b != NULL)   *ret_b   = b;    /* b is -1 if allow_begin is FALSE. */
  if (ret_bsc != NULL) *ret_bsc = bsc;  /* bsc is IMPOSSIBLE if allow_begin is FALSE */

  /* If the caller doesn't want the matrix, free it (saving the decks in the pool!)
   * Else, pass it back to him.
   */
  if (ret_alpha == NULL) {
    for (v = vroot; v <= vend; v++) /* be careful of our reuse of the end deck -- free it only once */
      if (alpha[v] != NULL) { 
	if (cm->sttype[v] != E_st) { free_vjd_deck(alpha[v], i0, j0); alpha[v] = NULL; }
	else end = alpha[v]; 
      }
    if (end != NULL) { free_vjd_deck(end, i0, j0); end = NULL; }
    free(alpha);
  } else *ret_alpha = alpha;

  free(touch);
  if (ret_shadow != NULL) *ret_shadow = shadow;
  return sc;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* never reached */
}

/* Function: insideT_qdb_me()
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
insideT_qdb_me(CM_t *cm, ESL_DSQ *dsq, int L, Parsetree_t *tr, 
	     int r, int z, int i0, int j0, 
	     int allow_begin, int *dmin, int *dmax)
{
  int       status;
  void   ***shadow;             /* the traceback shadow matrix */
  float     sc;			/* the score of the CYK alignment */
  ESL_STACK *pda;                /* stack that tracks bifurc parent of a right start */
  int       v,j,d,i;		/* indices for state, j, subseq len */
  int       k;			
  int       y, yoffset;
  int       bifparent;
  int       b;
  float     bsc;
  int       dp;                 /* dp: d' d offset in current state v's band; dp = d - dmin[v] */
  int       kp;                 /* dp: k' k offset in current state v's band; kp = k - dmin[v] */

  sc = inside_qdb_me(cm, dsq, L, r, z, i0, j0, 
		   BE_EFFICIENT,	/* memory-saving mode */
		   NULL, NULL,	        /* manage your own matrix, I don't want it */
		   &shadow,		/* return a shadow matrix to me. */
		   allow_begin,         /* TRUE to allow local begins */
		   &b, &bsc,	        /* if allow_begin is TRUE, gives info on optimal b */
		   dmin, dmax);

  pda = esl_stack_ICreate();
  if(pda == NULL) goto ERROR;
  v = r;
  j = j0;
  i = i0;
  d = j0-i0+1;

  while (1) {
    if(v == cm->M)
      dp = d;
    else
      dp = d - dmin[v];
    if(v != cm->M)
      {
	assert(d <= dmax[v]);
	assert(d >= dmin[v]);
      }
    if (cm->sttype[v] == B_st) {
      assert(v >= 0);
      kp = ((int **) shadow[v])[j][dp];   /* kp = offset len of right fragment */
      z = cm->cnum[v];
      k = kp + dmin[z];  /* k = len of right fragment */
      
      /* Store info about the right fragment that we'll retrieve later:
       */
      if((status = esl_stack_IPush(pda, j)) != eslOK) goto ERROR;	/* remember the end j    */
      if((status = esl_stack_IPush(pda, k)) != eslOK) goto ERROR;	/* remember the subseq length k */
      if((status = esl_stack_IPush(pda, tr->n-1)) != eslOK) goto ERROR;	/* remember the trace index of the parent B state */
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
      /* Note: we don't pop dp below, but d, because we're either in an E state
       * in which case d must be 0, or the EL state, which has no
       * dmin and dmax band, so if we pop dp and add dmin[v] to get d,
       * we'll f*** everything up, as Sam Griffiths-Jones found
       * when preparing Rfam 8.0 on 08.04.06.
       */
      esl_stack_IPop(pda, &d);
      esl_stack_IPop(pda, &j);
      v = tr->state[bifparent];	/* recover state index of B */
      y = cm->cnum[v];		/* find state index of right S */
      i = j-d+1;
				/* attach the S to the right */
      InsertTraceNode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, y);
      v = y;
    } else {
      yoffset = ((char **) shadow[v])[j][dp];
      if((((int) yoffset) != USED_LOCAL_BEGIN) && (((int) yoffset) != USED_EL))
	{
	  if(!((yoffset >= 0) && yoffset <= cm->M))
	    y = cm->cfirst[v] + yoffset;
	}
      if((yoffset != USED_LOCAL_BEGIN) && (yoffset != USED_EL))
	assert(yoffset >= 0 &&  yoffset <= cm->M);
      switch (cm->sttype[v]) {
      case D_st:            break;
      case MP_st: i++; j--; break;
      case ML_st: i++;      break;
      case MR_st:      j--; break;
      case IL_st: i++;      break;
      case IR_st:      j--; break;
      case S_st:            break;
      default:    cm_Fail("'Inconceivable!'\n'You keep using that word...'");
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

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0.; /* NEVERREACHED */
}


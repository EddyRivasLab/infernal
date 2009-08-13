/* cplan9.h 
 *
 Support for a CM plan 9 HMM architecture.
 * Data structures were morphed from their plan 7 analogs.
 * 
 * EPN 02.26.06 (based on SRE's HMMER's plan7.h)
 *
 * A few functions were hijacked from HMMER 2.4 and placed
 * here without modification. These first 4 are all from 
 * HMMER::mathsupport.c
 * 
 * ILogSum() (and auxiliary funcs associated with it)
 * Score2Prob()
 * Prob2Score()
 * Scorify()
 *
 * And one from HMMER:plan7.c
 * DegenerateSymbolScore()
 */  

#include "structs.h"
#include "config.h"

/* Functions stolen from HMMER-2.4::mathsupport.c */
extern int   ILogsum(int p1, int p2);
extern int   Prob2Score(float p, float null);
extern float Score2Prob(int sc, float null);
extern float Scorify(int sc);
extern int   DegenerateSymbolScore(float *p, float *null, int ambig);

#ifndef CPLAN9_INCLUDED
#define CPLAN9_INCLUDED

/* Structure: cm_plan9_s
 * 
 * 03.10.06 EPN: Original intended use of CM plan 9 structure is to read a CM
 * file, and build CM plan 9 HMM based on the CM, first by determining the 
 * probabilities for each state of the HMM, and then logoddsifying the model. 
 *
 * Declaration of a CM Plan 9 profile-HMM structure.
 * Modified from a plan 7 (with (hopefully) minimal change) to mirror a CM 
 * as closely as possible.
 * 
 * The model has two forms:
 * 1. The "core" model has 0..M nodes, node 0 is special, its "match" state
 *    is really state B (which is forced silent by having hmm->mat[0] = NULL and
 *    hmm->msc[0] = NULL), its "insert" state is really state N (with emission
 *    probs hmm->ins[0]), and it has NO DELETE STATE. 
 * 
 *    hmm->t[0][CTMM]: 0. (B->M_1 transition is hmm->begin[1])
 *    hmm->t[0][CTMI]: transition from B to N (I_0); 
 *    hmm->t[0][CTMD]: transition from B to D_1;
 *    hmm->t[0][CTME]: null (transition from B to an EL state is impossible)
 *    hmm->t[0][CTIM]: transition from N to M_1;
 *    hmm->t[0][CTII]: N self transition; 
 *    hmm->t[0][CTID]: N -> D_1
 *    hmm->t[0][CTDM]: null
 *    hmm->t[0][CTDI]: null
 *    hmm->t[0][CTDD]: null
 *    
 *    t[0..M] are the state transition probs. t[k][CTME] is an
 *    end-local probability, the EL states can only be reached by a
 *    subset of match states, this probability is -INFTY for states
 *    that can't reach the EL. 
 *
 *    t[M] are special, because this node transits to the end (E
 *    state). The E state is (sort-of) treated as match state M+1, as
 *    t[M][CTIM] is the transition from I_M to E, t[M][CTDM] is the
 *    transition from D_M to E. However, t[M][CTMM] is always 0.0,
 *    the transition from M_M to E is end[hmm->M]; t[M][CTMD],
 *    t[M][CTDD], t[M][CTDI] are set as 0.0.
 *    
 *    mat[1..M] are match emission probs.
 *    ins[0..M] are insert emission probs.  (ins[0] is state N emission probs)
 *
 *    The CM_PLAN9_HASPROB flag is up when these all correspond to a fully normalized
 *    profile HMM.
 *    
 * 2. The "logoddsified" model is the configured model, converted to
 *    integer score form and ready for alignment algorithms. 
 *    bsc, esc scores correspond to begin, and end probabilities.
 *    
 *    The CPLAN9_HASBITS flag is up when both of these are ready for
 *    alignment.
 *    
 */
typedef struct cplan9_s {
  /* The main model in probability form: data-dependent probabilities.
   * Transition probabilities are usually accessed as a
   *   two-D array: hmm->t[k][CTMM], for instance. They are allocated
   *   such that they can also be stepped through in 1D by pointer
   *   manipulations, for efficiency in DP algorithms.
   * CPLAN9_HASPROBS flag is raised when these probs are all valid.
   */
  const ESL_ALPHABET *abc;      /* pointer to the alphabet, usually points to cm->abc */
  int     M;                    /* length of the model (# nodes)        +*/
  float **t;                    /* transition prob's. t[0..M][0..9]   +*/
  float **mat;                  /* match emissions.  mat[1..M][0..3]   +*/ 
  float **ins;                  /* insert emissions. ins[0..M][0..3] +*/

  /* The unique states of CM Plan 9 in probability form.
   * These are the algorithm-dependent, data-independent probabilities.
   * Some parts of the code may briefly use a trick of copying tbd1
   *   into begin[0]; this makes it easy to call FChoose() or FNorm()
   *   on the resulting vector. However, in general begin[0] is not
   *   a valid number.
   * CPLAN9_HASALG flag is up when these probs are all valid.
   */
  float  *begin;                /* 1..M B->M state transitions                +*/
  float  *end;                  /* 1..M M->E state transitions (!= a dist!)   +*/

  /* The model's log-odds score form.
   * These are created from the probabilities by LogoddsifyHMM_cp9().
   * By definition, null[] emission scores are all zero.
   * Note that emission distributions are over possible alphabet symbols,
   * not just the unambiguous protein or DNA alphabet: we
   * precalculate the scores for all IUPAC degenerate symbols we
   * may see. 
   *
   * Note the reversed indexing on msc, isc, tsc -- for efficiency reasons.
   * They're not probability vectors any more so we can reorder them
   * without wildly complicating our life.
   * 
   * The _mem ptrs are where the real memory is alloc'ed and free'd,
   * as opposed to where it is accessed.
   * This came in with Erik Lindahl's altivec port; it allows alignment on
   * 16-byte boundaries. In the non-altivec code, this is just a little
   * redundancy; tsc and tsc_mem point to the same thing, for example.
   * 
   * CPLAN9_HASBITS flag is up when these scores are valid.
   */
  int  **tsc;                   /* transition scores     [0.9][0.M]       +*/
  int  **msc;                   /* match emission scores [0.MAXDEGEN-1][1.M] +*/
  int  **isc;                   /* ins emission scores   [0.MAXDEGEN-1][0.M] +*/
  int   *bsc;                   /* begin transitions     [1.M]              +*/
  int   *esc;			/* end transitions       [1.M]              +*/
  int   *tsc_mem, *msc_mem, *isc_mem, *bsc_mem, *esc_mem;

  /* The null model probabilities.
   */
  float *null;         /* "random sequence" emission prob's     +*/
  float  p1;                       /* null model loop probability           +*/
  float  el_self;                  /* EL transition self loop probability    */
  int    el_selfsc;                /* EL transition self loop score          */
  int   *has_el;                   /* has_el[k] is TRUE if node k has an EL state */
  int   *el_from_ct;               /* el_from_ct[k] is the number of HMM nodes kp
				    * where a transition from kp's EL state to k's
				    * match state is valid. */
  int  **el_from_idx;              /* [0..M+1][] el_from_idx[k] is an array of 
				    * size el_from_idx[k] each element is a node 
				    * kp where a transition from kp's EL state 
				    * to k's match state is allowed */
  int  **el_from_cmnd;             /* [0..M+1][] el_from_cmnd[k] is an array of 
				    * size el_from_idx[k] element i is the CM
				    * node that the EL transition to k to 
				    * el_from_idx[k][i] corresponds with, used
				    * only for building alignments from traces. */
  int flags;                       /* bit flags indicating state of HMM, valid data +*/
} CP9_t;

/* Flag codes for cplan9->flags.
 */
#define CPLAN9_HASBITS     (1<<0)    /* raised if model has log-odds scores      */
#define CPLAN9_HASPROB     (1<<1)    /* raised if model has probabilities        */
#define CPLAN9_LOCAL_BEGIN (1<<2)    /* raised if model has local begins turned on */
#define CPLAN9_LOCAL_END   (1<<3)    /* raised if model has S/W local ends turned on */
#define CPLAN9_EL          (1<<4)    /* raised if model has EL local ends turned on */

/* Indices for CM Plan9 main model state transitions.
 * Used for indexing hmm->t[k][]
 * mnemonic: Cm plan 9 Transition from Match to Match = CTMM
 */
#define CTMM  0
#define CTMI  1
#define CTMD  2
#define CTME  3
#define CTIM  4
#define CTII  5
#define CTID  6
#define CTDM  7
#define CTDI  8
#define CTDD  9

/* Declaration of CM Plan9 dynamic programming matrix structure.
 */
typedef struct cp9_dpmatrix_s {
  int **mmx;			/* match scores  [0.1..N][0..M] */
  int **imx;			/* insert scores [0.1..N][0..M] */
  int **dmx;			/* delete scores [0.1..N][0..M] */
  int **elmx;			/* end local scores [0.1..N][0..M] */
  int  *erow;                   /* score for E state [0.1..N] */
  /* Hidden ptrs where the real memory is kept; this trick was
   * introduced by Erik Lindahl with the Altivec port; it's used to
   * align xmx, etc. on 16-byte boundaries for cache optimization.
   */
  void *mmx_mem, *imx_mem, *dmx_mem, *elmx_mem;

  int *  workspace;      /* Workspace for altivec (aligned ptr)    */
  int *  workspace_mem;  /* Actual allocated pointer for workspace */
  
  /* The other trick brought in w/ the Lindahl Altivec port; dp matrix
   * is retained and grown, rather than reallocated for every HMM or sequence.
   * Keep track of current allocated-for size in rows (sequence length N)
   * and columns (HMM length M). Also keep track of pad sizes: how much
   * we should overallocate rows or columns when we reallocate. If pad = 0,
   * then we're not growable in this dimension.
   */
  int maxN;			/* alloc'ed for seq of length N; N+1 rows */
  int maxM;			/* alloc'ed for HMM of length M; M+1 cols */

  int padN;			/* extra pad in sequence length/rows */
  int padM;			/* extra pad in HMM length/columns   */
} CP9_dpmatrix_t;


/* CM Plan 9 model state types
 * used in traceback structure
 */
#define CSTBOGUS 0
#define CSTM     1
#define CSTD     2
#define CSTI     3
#define CSTB     4  /* M_0 the B state */
#define CSTE     5  /* the end state, M_(k+1) */
#define CSTEL    6  /* an EL (end local) state */
/* Structure: cp9trace_s
 * 
 * Traceback structure for alignments of model to sequence.
 * Each array in a trace_s is 0..tlen-1.
 * Element 0 is always to M_0 (match state of node 0)
 * Element tlen-1 is always to the E_st
 */
typedef struct cp9trace_s {
  int   tlen;                   /* length of traceback                           */
  char *statetype;              /* state type used for alignment                 */
  int  *nodeidx;                /* idx of aligned node, 0..M if M or I 1..M if D */
  int  *pos;                    /* position in dsq, 1..L, or 0 if none           */ 
} CP9trace_t;


#endif /* CPLAN9_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/


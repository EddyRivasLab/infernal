#ifndef STRUCTSH_INCLUDED
#define STRUCTSH_INCLUDED

/* structs.h
 * SRE, 28 Feb 2000 [Seattle]
 * SVN $Id$
 * 
 * Declarations of structures and global variables;
 * definitions of constants; and macros.
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include "esl_config.h"
#include "config.h"

#include "easel.h"
#include "esl_dirichlet.h"
#include "esl_random.h"
#include "esl_sqio.h"

#define USE_NEWLOGSUM 1
#define USE_OLDLOGSUM 0

/* various default parameters for CMs and CP9 HMMs */ 
#define DEFAULT_CM_CUTOFF 0.1
#define DEFAULT_CM_CUTOFF_TYPE E_CUTOFF
#define DEFAULT_CP9_CUTOFF 0.0
#define DEFAULT_CP9_CUTOFF_TYPE SCORE_CUTOFF
#define DEFAULT_MIN_CP9_E_CUTOFF 1.0
#define DEFAULT_BETA   0.0000001
#define DEFAULT_TAU    0.0000001
#define DEFAULT_HMMPAD 0
#define DEFAULT_PBEGIN 0.05  /* EPN 06.29.07 (formerly 0.5) */
#define DEFAULT_PEND   0.05  /* EPN 06.29.07 (formerly 0.5) */
#define DEFAULT_ETARGET 0.54 /* EPN 07.10.07 (formerly (v0.7->v0.8)= 2.-0.54 = 1.46 */

/* default num samples for CM and CP9 E-values */
#define DEFAULT_NUM_SAMPLES 1000

#define GC_SEGMENTS 101                   /* Possible integer GC contents */

#define BE_EFFICIENT  0		/* setting for do_full: small memory mode */
#define BE_PARANOID   1		/* setting for do_full: keep whole matrix, perhaps for debugging */

/* Constants for type of cutoff */
#define SCORE_CUTOFF 0
#define E_CUTOFF     1

/* Alphabet information is declared here, and defined in globals.c.
 */
#define MAXABET     4
#define CP9MAXABET  4 /* should be same as MAXABET */
#define MAXDEGEN   17

/* We're moderately paranoid about underflow and overflow errors, so
 * we do some checking on the magnitude of the scores.
 * 
 * IMPOSSIBLE, the "-infinity" value in a DP matrix must be > -FLT_MAX/3, so that 
 *   we can add three of them together (alpha + tsc + esc) and not get an
 *   underflow error. ANSI guarantees us FLT_MAX >= 1e37. 
 * 
 * MAXSCOREVAL is the maximum value we allow in any alpha cell, tsc, or esc.
 *
 * IMPROBABLE must be > IMPOSSIBLE + 2* MAXSCOREVAL; such that adding
 *  any two valid score values (say, an alpha and a tsc) to IMPOSSIBLE
 *  gives us a number < IMPROBABLE, and we can reset it to IMPOSSIBLE.
 *  
 * NOT_IMPOSSIBLE() exists because we can't compare floating point by 
 * equality.
 *  
 * sreLOG2() exists because we want to work in bits, and we will need
 * to take log(0).
*/
#define IMPOSSIBLE  -1e36
#define MAXSCOREVAL  1e35
#define IMPROBABLE  -5e35
#define NOT_IMPOSSIBLE(x)  ((x) > -9.999e35) 
#define sreLOG2(x)  ((x) > 0 ? log(x) * 1.44269504 : IMPOSSIBLE)
#define sreEXP2(x)  (exp((x) * 0.69314718 )) 
#define epnEXP10(x) (exp((x) * 2.30258509 ))
#define NOTZERO(x)  (fabs(x - 0.) > -1e6)
#define INFTY       987654321   /* infinity for purposes of integer DP cells       */

/***********************************************************************************
 * CM Plan 9 HMM information                                                       */

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
   * CPLAN9_HASPROB flag is raised when these probs are all valid.
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
  float  null[CP9MAXABET];         /* "random sequence" emission prob's     +*/
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

/************************************************************************************
 * End of CM Plan 9 HMM information.
 ************************************************************************************/

/* CM State types. (cm->sttype[])
 */
#define MAXCONNECT 6            /* maximum number of states per node */

#define D_st     0		/* delete       */
#define MP_st    1		/* match-pair   */
#define ML_st    2		/* match-left   */
#define MR_st    3		/* match-right  */
#define IL_st    4		/* insert-left  */
#define IR_st    5		/* insert-right */
#define S_st     6		/* start        */
#define E_st     7		/* end          */
#define B_st     8		/* bifurcation  */
#define EL_st    9              /* local end    */

/* CM Node types (8) (cm->ndtype[])
 */
#define NODETYPES 8		

#define DUMMY_nd -1
#define BIF_nd    0
#define MATP_nd   1
#define MATL_nd   2
#define MATR_nd   3
#define BEGL_nd   4		
#define BEGR_nd   5
#define ROOT_nd   6		
#define END_nd    7

/* CM Unique state identifiers  (cm->stid[])
 */
#define UNIQUESTATES 21

#define DUMMY   -1
#define ROOT_S  0
#define ROOT_IL 1
#define ROOT_IR 2
#define BEGL_S  3
#define BEGR_S  4
#define BEGR_IL 5
#define MATP_MP 6
#define MATP_ML 7
#define MATP_MR 8
#define MATP_D  9
#define MATP_IL 10
#define MATP_IR 11
#define MATL_ML 12
#define MATL_D  13
#define MATL_IL 14
#define MATR_MR 15
#define MATR_D  16
#define MATR_IR 17
#define END_E   18
#define BIF_B   19
#define END_EL  20

/* Flags used in InsertTraceNode()
 */
#define TRACE_LEFT_CHILD  1
#define TRACE_RIGHT_CHILD 2

/* Flags used to define PDA moves, 
 * in display.c and emit.c (if not elsewhere)
 *
 */
#define PDA_RESIDUE 0
#define PDA_STATE   1
#define PDA_MARKER  2

/* Structure: CP9Map_t
 * Incept:    EPN, 10.23.06
 *
 * Maps a CM to a CM plan 9 HMM and vice versa. 
 *    Consensus positions are indexed 1..hmm_M.
 *
 * The lpos and rpos arrays are somewhat redundant with 
 * CMEmitMap_t, but they're not identical. I didn't realize emitmap existed
 * prior to implementing CP9 HMMs - if I had I would have used emitmap's, but
 * it's difficult to go back and use emitmap's now.
 * 
 * See emitmap.c for implementation and more documentation.
 */
typedef struct cp9map_s {
  int   *nd2lpos;   /* [0..cm_nodes] left bound of consensus for subtree under this nd,
		     *               -1 for non-{MATL|MATR|MATP} nodes */
  int   *nd2rpos;   /* [0..cm_nodes] right bound of consensus for subtree under this nd  
		     *               -1 for non-{MATL|MATR|MATP} nodes */
  int   *pos2nd;    /* [1..clen], the MATL, MATR or MATP CM node that maps to this 
		     *            consensus column */
  int  **cs2hn;     /* [0..cm_M][0..1], 1 or 2 HMM nodes that maps to this CM state 
		     *                  [x][1] is -1 if state x maps to only 1 HMM node */
  int  **cs2hs;     /* [0..cm_M][0..1], 1 or 2 HMM states (0=MATCH, 1=INSERT, 2=DELETE)
		     *            that maps to this CM state 
		     *            [x][0] corresponds to the HMM node in cs2hn[x][0],
		     *            [x][1] corresponds to the HMM node in cs2hn[x][1] (or -1) */
  int ***hns2cs;    /* [0..clen][0..2][0..1]
		     * hns2cs[x][y][0] 1st CM state that maps to HMM node x state type y (0,1,2)
		     * hns2cs[x][y][1] 2nd CM state that maps to HMM node x state type y (0,1,2)
		     *                 -1 if only 1 CM state maps to HMM node x state y */
  int    hmm_M;     /* consensus length, the number of HMM nodes */
  int    cm_M;      /* number of states in the CM this HMM maps to */
  int    cm_nodes;  /* number of nodes in the CM this HMM maps to */
} CP9Map_t;

/* Structure GumbelInfo_t
 */
typedef struct gumbelinfo_s {
  int    N;             /* number of samples stats calc'ed from        */
  int    L;             /* length of samples stats calc'ed from        */
  double mu;		/* location param for gumbel, calced w/K,lambda*/
  double lambda;	/* scale param gumbel                          */
} GumbelInfo_t;

/* Structure CP9FThresh_t: CP9 HMM filter thresholds, determined empirically
 * by sampling from the CM
 */
typedef struct cp9filterthr_s {
  int   N;             /* number of CM hits used to get threshold ((N*fraction) passed)*/
  float cm_eval;       /* CM E-value threshold, we rejected worse than   */
  float l_eval;        /*  local CP9 scanning E-value threshold    */
  float g_eval;        /* glocal CP9 scanning E-value threshold    */
  float l_F;           /* fraction of empirical CM hits survive filter for l_eval cutoff */
  float g_F;           /* fraction of empirical CM hits survive filter for g_eval cutoff */
  int   db_size;       /* db size used to calculate Gum mu for *_eval calculations */
  int   was_fast;      /* TRUE if hacky fast method for calcing thresholds was used */
} CP9FilterThr_t;


/* Structure SubFilterInfo_t: Information on possible sub CM filters for a CM.                           
 * States of a CM are grouped into 'start groups'. There is one start group                              
 * for each start state of the CM. A 'start group' begins with a start state and ends                    
 * with a E or B state, and includes all states in between.                                              
 */                                                                                                      
typedef struct subfilterinfo_s {
  int    M;            /* # states in the CM */                                                          
  int    nstarts;      /* # start states (and start groups) in the CM */                                 
  int    ncands;       /* number of candidate states, these *could* be sub CM roots */                   
  double beta;         /* beta used for calculating avglenA */                                           
  float  minlen;       /* minimum average length (avglen) a candidate state must have */                 
  int   *iscandA;      /* [0..v..cm->M-1] TRUE if state v is a candidate sub CM root, FALSE otherwise */   
  float *avglenA;      /* [0..v..cm->M-1] average length of a hit rooted at v (from QDB) */                
  int   *startA;       /* [0..i..cm->M-1] start group this state belongs to */                               
  int   *firstA;       /* [0..i..nstarts-1], first state in start state i's group */                     
  int   *lastA;        /* [0..i..nstarts-1], last state in start state i's group */                      
  int  **withinAA;     /* [0..i..nstarts-1][0..j..nstarts-1] = TRUE if start state j's group             
                        * is within start state i's group.                                               
                        *  emap->startA[cm->nodemap[i]]->lpos < emap->startA[cm->nodemap[j]]->lpos  &&   
                        *  emap->endA  [cm->nodemap[i]]->rpos > emap->endA  [cm->nodemap[j]]->rpos       
                        */			
} SubFilterInfo_t;

/* Structure CMStats_t
 */
typedef struct cmstats_s {
  int np;                    /* number of partitions, df: 1 */
  int *ps;                   /* start GC content [0..100] of each partition */
  int *pe;                   /* end   GC content [0..100] of each partition */
  int gc2p[GC_SEGMENTS];     /* map from GC content to partition number     */
  GumbelInfo_t ***gumAA;     /* [0..NSTATMODES-1][0..np-1] */
  CP9FilterThr_t **fthrA;    /* [0..NFTHRMODES-1] */
} CMStats_t;


/* Stat modes, 
 * 0..NSTATMODES-1 are first dimension of cmstats->gumAA 
 * 0..NFTHRMODES-1 are only dimension cmstats->fthrA 
 */
#define CM_LC 0  
#define CM_GC 1
#define CM_LI 2
#define CM_GI 3
#define CP9_L 4
#define CP9_G 5
#define NGUMBELMODES 6
#define NFTHRMODES   4
#define NCMMODES     4
#define NCP9MODES    2

/* Structure: CM_t
 * Incept:    SRE, 9 Mar 2000 [San Carlos CA]
 * 
 * A covariance model. M states, arranged logically as a directed graph
 * (on a binary tree backbone); arranged physically as a set of arrays 0..M-1.
 *
 * State 0 is always the root state. State M-1 is always an end state.
 * 
 * EPN 12.19.06: added arrays to hold integer log-odds scores for faster 
 * inside/outside
 * 
 */
typedef struct cm_s {			
			/* General information about the model:            */
  char *name;		/*   name of the model                             */
  char *acc;		/*   optional accession number for model, or NULL  */
  char *desc;		/*   optional description of the model, or NULL    */
  char *annote;         /*   consensus column annotation line, or NULL     */ /* ONLY PARTIALLY IMPLEMENTED, BEWARE */

  /* new for 1.0 */
  char  *comlog;	/*   command line(s) that built model      (mandatory) */ /* String, \0-terminated */
  int    nseq;		/*   number of training sequences          (mandatory) */
  float  eff_nseq;	/*   effective number of seqs (<= nseq)    (mandatory) */
  char  *ctime;		/*   creation date                         (mandatory) */
  float  ga;	        /*   per-seq/per-domain gathering thresholds (bits) (CMH_GA) */
  float  tc;            /*   per-seq/per-domain trusted cutoff (bits)       (CMH_TC) */
  float  nc;	        /*   per-seq/per-domain noise cutoff (bits)         (CMH_NC) */


			/* Information about the null model:               */
  float *null;          /*   residue probabilities [0..3]                  */

			/* Information about the state type:               */
  int   M;		/*   number of states in the model                 */
  int   clen;		/*   consensus length (2*MATP+MATL+MATR)           */
  char *sttype;		/*   type of state this is; e.g. MP_st             */
  int  *ndidx;		/*   index of node this state belongs to           */
  char *stid;		/*   unique state identifier; e.g. MATP_MP         */

			/* Information about its connectivity in CM:       */
  int  *cfirst;		/*   index of left child state                     */
  int  *cnum;		/*   overloaded: for non-BIF: # connections;       */
			/*               for BIF: right child S_st         */
  int  *plast;          /*   index to first parent state                   */
  int  *pnum;           /*   number of parent connections                  */

			/* Information mapping nodes->states               */
  int   nodes;		/*   number of nodes in the model                  */
  int  *nodemap;        /*   nodemap[5] = idx first state, node 5          */
  char *ndtype;		/*   type of node, e.g. MATP_nd                    */

                        /* Parameters of the probabilistic model:          */
  float **t;		/*   Transition prob's [0..M-1][0..MAXCONNECT-1]   */
  float **e;		/*   Emission probabilities.  [0..M-1][0..15]      */
  float  *begin;	/*   Local alignment start probabilities [0..M-1]  */
  float  *end;		/*   Local alignment ending probabilities [0..M-1] */

			/* Parameters of the log odds model:               */
  float **tsc;		/*   Transition score vector, log odds             */
  float **esc;		/*   Emission score vector, log odds               */
  float *beginsc;	/*   Score for ROOT_S -> state v (local alignment) */
  float *endsc;   	/*   Score for state_v -> EL (local alignment)     */

			/* Scaled int parameters of the log odds model:    */
  int  **itsc;		/*   Transition score vector, scaled log odds int  */
  int  **iesc;		/*   Emission score vector, scaled log odds int    */
  int   *ibeginsc;      /*   Score for ROOT_S -> state v (local alignment) */
  int   *iendsc;  	/*   Score for state_v -> EL (local alignment)     */

  int    flags;		/* status flags                                    */

  /* query dependent bands (QDB) on subsequence lengths at each state                         */
  int   *dmin;          /* minimum d bound for each state v; [0..v..M-1] (NULL if non-banded) */
  int   *dmax;          /* maximum d bound for each state v; [0..v..M-1] (NULL if non-banded) */
  double beta;          /* tail loss probability for QDB                                      */
  double tau;           /* tail loss probability for HMM target dependent banding             */

  /* added by EPN, Tue Jan  2 14:24:08 2007 */
  int       config_opts;/* model configuration options                                        */
  int       align_opts; /* alignment options                                                  */
  int       search_opts;/* search options                                                     */
  CP9_t    *cp9;        /* a CM Plan 9 HMM, always built when the model is read from a file   */
  CP9Map_t *cp9map;     /* the map from the Plan 9 HMM to the CM and vice versa               */
  int       enf_start;  /* if(cm->config_opts & CM_CONFIG_ENFORCE) the first posn to enforce, else 0 */
  char     *enf_seq;    /* if(cm->config_opts & CM_CONFIG_ENFORCE) the subseq to enforce, else NULL  */
  float     enf_scdiff; /* if(cm->config_opts & CM_CONFIG_ENFORCE) the difference in scoring  *
			 * cm->enfseq b/t the non-enforced & enforced CMs, this is subtracted *
			 * from bit scores in cmsearch before assigned E-value stats which    *
			 * are always calc'ed (histograms built) using non-enforced CMs/CP9s  */
  float     sc_boost;   /* value added to CYK bit scores during search (usually 0.)           */
  float cp9_sc_boost;   /* value added to Forward bit scores during CP9 search (usually 0.)   */
  float     ffract;     /* desired filter fraction (0.99 -> filter out 99% of db), default: 0.*/
  float    *root_trans; /* transition probs from state 0, saved IFF zeroed in ConfigLocal()   */
  int       hmmpad;     /* if(cm->search_opts & CM_SEARCH_HMMPAD) # of res to -/+ from i/j    */
  float     pbegin;     /* local begin prob to spread across internal nodes for local mode    */
  float     pend;       /* local end prob to spread across internal nodes for local mode      */
  
  /* search cutoffs */
  int       cutoff_type;/* either SC_CUTOFF or E_CUTOFF                                       */
  float     cutoff;     /* min bit score or max E val to keep in a scan (depending on cutoff_type) */
  int   cp9_cutoff_type;/* either SC_CUTOFF or E_CUTOFF                                       */
  float cp9_cutoff;     /* min bit score or max E val to keep from a CP9 scan                 */
  
  int    W;             /* max d: max size of a hit (EPN 08.18.05) */
  float  el_selfsc;     /* score of a self transition in the EL state
			 * the EL state emits only on self transition (EPN 11.15.05)*/
  int   iel_selfsc;     /* scaled int version of el_selfsc         */

  CMStats_t *stats;     /* holds Gumbel stats and HMM filtering thresholds */

  /* From 1.0-ification, based on HMMER3 */
  const  ESL_ALPHABET *abc;     /* ptr to alphabet info (cm->abc->K is alphabet size)*/
} CM_t;

/* status flags, cm->flags */
#define CMH_BITS               (1<<0)  /* CM has valid log odds scores             */
#define CMH_ACC                (1<<1)  /* accession number is available            */
#define CMH_DESC               (1<<2)  /* description exists                       */
#define CMH_GA                 (1<<3)  /* gathering threshold exists               */
#define CMH_TC                 (1<<4)  /* trusted cutoff exists                    */
#define CMH_NC                 (1<<5)  /* noise cutoff exists                      */

#define CM_LOCAL_BEGIN         (1<<6)  /* Begin distribution is active (local ali) */
#define CM_LOCAL_END           (1<<7)  /* End distribution is active (local ali)   */
#define CM_GUMBEL_STATS        (1<<8)  /* Gumbel stats for local/glocal CYK/Ins set*/
#define CM_FTHR_STATS          (1<<9)  /* CP9 HMM filter threshold stats are set   */
#define CM_QDB                 (1<<10) /* query-dependent bands, QDB valid         */
#define CM_CP9                 (1<<11) /* CP9 HMM is valid in cm->cp9              */
#define CM_CP9STATS            (1<<12) /* CP9 HMM has Gumbel stats                 */
#define CM_IS_SUB              (1<<13) /* the CM is a sub CM                       */
#define CM_ENFORCED            (1<<14) /* CM is reparam'ized to enforce a subseq   */
#define CM_IS_RSEARCH          (1<<15) /* the CM was parameterized a la RSEARCH    */
#define CM_RSEARCHTRANS        (1<<16) /* CM has/will have RSEARCH transitions     */
#define CM_RSEARCHEMIT         (1<<17) /* CM has/will have RSEARCH emissions       */

/* model configuration options, cm->config_opts */
#define CM_CONFIG_LOCAL        (1<<0)  /* configure the model for local alignment  */
#define CM_CONFIG_HMMLOCAL     (1<<1)  /* configure the CP9   for local alignment  */
#define CM_CONFIG_HMMEL        (1<<2)  /* configure the CP9   for local alignment  */
#define CM_CONFIG_ENFORCE      (1<<3)  /* enforce a subseq be incl. in each parse  */
#define CM_CONFIG_ENFORCEHMM   (1<<4)  /* build CP9 HMM to only enforce subseq     */
#define CM_CONFIG_ZEROINSERTS  (1<<5)  /* make all insert emissions equiprobable   */
#define CM_CONFIG_QDB          (1<<6)  /* calculate query dependent bands          */

/* alignment options, cm->align_opts */
#define CM_ALIGN_NOSMALL       (1<<0)  /* DO NOT use small CYK D&C                 */
#define CM_ALIGN_QDB           (1<<1)  /* use QD bands                             */
#define CM_ALIGN_HBANDED       (1<<2)  /* use HMM bands                            */
#define CM_ALIGN_SUMS          (1<<3)  /* if using HMM bands, use posterior sums   */
#define CM_ALIGN_SUB           (1<<4)  /* build a sub CM for each seq to align     */
#define CM_ALIGN_FSUB          (1<<5)  /* build a 'full sub' CM for each seq       */
#define CM_ALIGN_HMMONLY       (1<<6)  /* use a CP9 HMM only to align              */
#define CM_ALIGN_INSIDE        (1<<7)  /* use Inside, not CYK                      */
#define CM_ALIGN_OUTSIDE       (1<<8)  /* use Outside, not CYK (for testing)       */
#define CM_ALIGN_POST          (1<<9)  /* do inside/outside and append posteriors  */
#define CM_ALIGN_TIME          (1<<10) /* print out alignment timings              */
#define CM_ALIGN_CHECKINOUT    (1<<11) /* check inside/outside calculations        */
#define CM_ALIGN_CHECKPARSESC  (1<<12) /* check parsetree score against aln alg sc */
#define CM_ALIGN_PRINTTREES    (1<<13) /* print parsetrees to stdout               */
#define CM_ALIGN_HMMSAFE       (1<<14) /* realign seqs w/HMM banded CYK bit sc < 0 */
#define CM_ALIGN_SCOREONLY     (1<<15) /* do full CYK/inside to get score only     */
#define CM_ALIGN_SAMPLE        (1<<16) /* sample parsetrees from the inside matrix */
#define CM_ALIGN_FLUSHINSERTS  (1<<17) /* flush inserts L/R like pre 1.0 infernal  */

/* search options, cm->search_opts */
#define CM_SEARCH_NOQDB        (1<<0)  /* DO NOT use QDB to search (QDB is default)*/
#define CM_SEARCH_HMMONLY      (1<<1)  /* use a CP9 HMM only to search             */
#define CM_SEARCH_HMMFILTER    (1<<2)  /* filter w/CP9 HMM, using forward/backward */
#define CM_SEARCH_HMMPAD       (1<<3)  /* +/- cm->hmmpad residues from j/i in HMMFB*/
#define CM_SEARCH_HMMRESCAN    (1<<4)  /* rescan HMM hits b/c Forward is inf length*/
#define CM_SEARCH_HBANDED      (1<<5)  /* use HMM bands for search                 */
#define CM_SEARCH_HMMSCANBANDS (1<<6)  /* filter w/CP9 HMM, and derive HMM bands   */
#define CM_SEARCH_SUMS         (1<<7)  /* if using HMM bands, use posterior sums   */
#define CM_SEARCH_INSIDE       (1<<8)  /* scan with Inside, not CYK                */
#define CM_SEARCH_TOPONLY      (1<<9)  /* don't search reverse complement          */
#define CM_SEARCH_NOALIGN      (1<<10) /* don't align hits, just report locations  */
#define CM_SEARCH_NULL2        (1<<11) /* use post hoc second null model           */
#define CM_SEARCH_CMSTATS      (1<<12) /* calculate E-value statistics for CM      */
#define CM_SEARCH_CP9STATS     (1<<13) /* calculate E-value stats for CP9 HMM      */
#define CM_SEARCH_FFRACT       (1<<14) /* filter to filter fraction cm->ffract     */
#define CM_SEARCH_RSEARCH      (1<<15) /* use RSEARCH parameterized CM             */
#define CM_SEARCH_CMGREEDY     (1<<16) /* use greedy alg to resolve CM overlaps    */
#define CM_SEARCH_HMMGREEDY    (1<<17) /* use greedy alg to resolve HMM overlaps   */

/* Structure: CMFILE
 * Incept:    SRE, Tue Aug 13 10:16:39 2002 [St. Louis]
 *
 * An open CM database for reading. 
 * (When writing, we just use a normal stdio.h FILE.)
 * API is implemented in cmio.c
 */
typedef struct cmfile_s {
  FILE     *f;                  /* open file for reading */
  char     *fname;              /* name of the CM file; [STDIN] if -           */
  ESL_SSI  *ssi;                /* ptr to open SSI index, or NULL if unavailable */
  int       is_binary;		/* TRUE if file is in binary format */
  int       byteswap;		/* TRUE if binary and we need to swap byte order */
  int       mode;		/* type of SSI offset (part of SSI API) */
  off_t     offset;             /* disk offset of the CM that was read last */
} CMFILE;


/* Structure: Parsetree_t
 * Incept:    SRE 29 Feb 2000 [Seattle]
 * 
 * Binary tree structure for storing a traceback of an alignment.
 * 
 * Also used for tracebacks of model constructions. Then, 
 * "state" is misused for a node (not state) index. 
 * 
 * Example of a traceback (from ParsetreeDump(), from a tRNA
 * model:
 * 
 * > DF6280
 * idx   emitl  emitr   state  nxtl  nxtr  prv   tsc   esc
 * ----- ------ ------ ------- ----- ----- ----- ----- -----
 *    0     1     74      0S      1    -1    -1 -0.58  0.00
 *    1     1     74A     3MR     2    -1     0 -0.74  0.41
 *    2     1G    73C     6MP     3    -1     1 -0.87  1.58
 * ...<snip>...
 *   11    10     66     54B     12    43    10  0.00  0.00
 *   12    10     44    124S     13    -1    11  0.00  0.00
 *   13    10     44    125B     14    28    12  0.00  0.00
 * ...<snip>...
 *   60    61U    61    120ML    61    -1    59 -0.22  0.87
 *   61    -1     -1    123E     -1    -1    60  0.00  0.00
 * ----- ------ ------ ------- ----- ----- ----- ----- -----
 *    
 * That is, emitl and emitr are always valid and always represent
 * the bounds of the subsequence accounted for by the parse
 * subtree rooted at this state. (Except for end states, which
 * are -1,-1). nxtl is always a valid state (again except for E
 * states, which are -1. nxtr is only != -1 for bifurcation states.
 *    
 * For reasons of malloc() efficiency, the binary tree is organized
 * in a set of arrays. 
 */
typedef struct parsetree_s {
  int *emitl;		/* i position in seq or ali (1..L or alen) */
  int *emitr;		/* j position in seq or ali (1..L or alen) */
  int *state;		/* y of state (0..M-1)                     */
  int *mode;		/* mode of state (used in marginal         *
                         * alignment), (0,1,2,3)                   */

  int *nxtl;		/* index in trace of left child            */
  int *nxtr;		/* index in trace of right child           */
  int *prv;		/* index in trace of parent                */

  int  n;		/* number of elements in use so far        */
  int  nalloc;		/* number of elements allocated for        */
  int  memblock;	/* size of malloc() chunk, # of elems      */
} Parsetree_t;


/* Structure: CMConsensus_t
 * Incept:    SRE, Thu May 23 16:55:04 2002 [St. Louis]
 * 
 * Created by display.c:CreateCMConsensus(). 
 * Preprocesses a CM into consensus information that is needed by
 * display.c:CreateFancyAli().
 *
 *   ct[x]:  Zuker-style ct map.
 *           indicates the pairing partner for consensus position ct[x].
 *           x can be 0..clen-1
 *           ct[x] is -1 (no partner) or 0..clen-1 (coord of partner)
 *           
 *   (lpos, rpos may be redundant w/ CMEmitMap_t now.)
 *   (off-by-one w.r.t. CMEmitMap_t; 1..clen is better)
 */
typedef struct consensus_s {
  char *cseq;           /* consensus sequence display string; 0..clen-1     */
  char *cstr;		/* consensus structure display string; 0..clen-1    */
  int  *ct;             /* Zuker-style ct pairing map; [0..clen-1]          */
  int  *lpos;		/* maps node->consensus position; 0..nodes-1        */
  int  *rpos;		/* maps node->consensus position; 0..nodes-1        */
  int   clen;		/* length of cseq, cstr                             */
} CMConsensus_t;

/* Structure: Fancyali_t
 * Incept:    SRE, May 2002 [St. Louis]
 * 
 * See display.c:CreateFancyAli(). 
 * An alignment of a CM to a target sequence, formatted for display.
 */
typedef struct fancyali_s {
  char *annote;         /* reference annotation line (NULL if unavail) */
  char *cstr;		/* CM consensus structure line                 */
  char *cseq;		/* CM consensus sequence line                  */
  char *mid;		/* alignment identity middle line              */
  char *aseq;		/* aligned target sequence                     */
  int  *scoord;		/* coords 1..L for aligned dsq chars           */
  int  *ccoord;		/* coords 1..clen for aligned consensus chars  */
  int   len;		/* len of the strings above                    */
  int   cfrom, cto;	/* max bounds in ccoord                        */
  int   sqfrom, sqto;	/* max bounds in scoord                        */
} Fancyali_t;

/* Structure: CMEmitMap_t
 * Incept:    SRE, Thu Aug  8 12:47:49 2002 [St. Louis]
 *
 * Maps model nodes to consensus positions.
 *    Consensus positions are indexed 1..clen.
 *    Each array (lpos, rpos, epos) is 0..nodes-1.
 *    Residues from an MP go into lpos and rpos in the consensus.
 *    Residues from an IL follow lpos.
 *    Residues from an IR precede rpos.
 *    Residues from an EL follow epos[nd] for the nd that went to EL.
 *    For nonemitters, rpos and lpos are a non-inclusive bound: for
 *      example, rpos[0], lpos[0] are 0,clen+1.
 *    There are no dummy values; all rpos, lpos, epos are valid coords
 *      0..clen+1 in the consensus.
 *
 * See emitmap.c for implementation and more documentation.
 */
typedef struct emitmap_s {
  int *lpos;           /* left bound of consensus for subtree under nd   */
  int *rpos;           /* right bound of consensus for subtree under nd  */
  int *epos;           /* EL inserts come after this consensus pos */
  int  clen;           /* consensus length */
} CMEmitMap_t;

/* Structure: CMSubMap_t
 * Incept:    EPN, 10.23.06
 *
 * Maps a template CM to a sub CM and vice versa. 
 *    Consensus positions are indexed 1..clen.
 *
 * The *node_cc_left and *node_cc_right arrays are redundant with CMEmitMap_t.
 * See emitmap.c for implementation and more documentation.
 */
typedef struct submap_s {
  int spos;            /* first consensus column this sub_cm models */
  int epos;            /* final consensus column this sub_cm models */
  int sstruct;         /* first consensus column this sub_cm models structure of */
  int estruct;         /* final consensus column this sub_cm models structure of */

  int **s2o_smap;      /* v = [0..sub_M-1] [0..1], orig_cm state(s) that maps to v */
  int **o2s_smap;      /* v = [0..orig_M-1][0..1], sub_cm  state(s) that maps to v */
  int  *s2o_id;        /* v = [0..sub_M-1] TRUE if sub_cm state v maps identically *
                        * to a orig_cm state (this will be s2o_smap[v][0])         */

  int  sub_clen;       /* consensus length orig_cm */
  int  orig_clen;      /* consensus length sub_cm  */

  int sub_M;           /* number of states in the sub CM */
  int orig_M;          /* number of states in the original CM */
} CMSubMap_t;


/* Structure: CMSubInfo_t
 * Incept:    EPN, 10.23.06
 *
 * Information on a sub CM, used for checking the sub CM 
 * construction procedure works.
 *    Consensus positions are indexed 1..clen.
 *
 */
typedef struct subinfo_s {
  int  *imp_cc;         /* [0..(epos-spos+1)] ret_imp_cc[k] != 0 if it is 
			* impossible for CP9 node (consensus column) k to be
			* calculated in the sub_cm to have distros to match the
			* corresponding CP9 node in the original CM - due to
			* topological differences in the architecture of the
			* sub_cm and orig_cm.
			*/
  int  *apredict_ct;   /* For an analytical test, the number of times we 
			* predict we'll fail the test for an HMM node for each 
			* of 6 cases - 6 different reasons we predict we'll fail.
			*/
  int  *spredict_ct;   /* For as sampling test, the number of times we 
			* predict we'll fail the test for an HMM node for each 
			* of 6 cases - 6 different reasons we predict we'll fail.
			*/
  int  *awrong_ct;     /* Subset of cases in apredict_ct that were incorrectly 
			* predicted. */
  int  *swrong_ct;     /* Subset of cases in spredict_ct that were incorrectly 
			* predicted. */
} CMSubInfo_t;

/* Structure: CP9Bands_t
 * Incept:    EPN, 10.27.06
 *
 * CP9 HMM bands and the resulting CM bands for HMM banded
 * CYK (or Inside or Outside) algorithms.
 *
 */
typedef struct cp9bands_s {
  int      hmm_M;             /* Number of nodes in CP9 HMM */
  int      cm_M;              /* Number of nodes in CM, the CP9 HMM was built from */

  /* data structures for hmm bands (bands on the hmm states) */
  int     *pn_min_m;          /* HMM band: minimum position node k match state will emit  */
  int     *pn_max_m;          /* HMM band: maximum position node k match state will emit  */
  int     *pn_min_i;          /* HMM band: minimum position node k insert state will emit */
  int     *pn_max_i;          /* HMM band: maximum position node k insert state will emit */
  int     *pn_min_d;          /* HMM band: minimum position that was emitted prior to entering
			       * node k delete state */
  int     *pn_max_d;          /* HMM band: maximum position that was emitted prior to entering
			       * node k delete state */
  int     *isum_pn_m;         /* [1..k..M] sum over i of log post probs from post->mmx[i][k]*/
  int     *isum_pn_i;         /* [1..k..M] sum over i of log post probs from post->imx[i][k]*/
  int     *isum_pn_d;         /* [1..k..M] sum over i of log post probs from post->dmx[i][k]*/

  /* arrays for CM state bands, derived from HMM bands */
  int *imin;                  /* [1..M] imin[v] = first position in band on i for state v*/
  int *imax;                  /* [1..M] imax[v] = last position in band on i for state v*/
  int *jmin;                  /* [1..M] jmin[v] = first position in band on j for state v*/
  int *jmax;                  /* [1..M] jmax[v] = last position in band on j for state v*/
  int **hdmin;                /* [v=1..M][0..(jmax[v]-jmin[v])] 
			       * hdmin[v][j0] = first position in band on d for state v, and position
			       * j = jmin[v] + j0.*/
  int **hdmax;                /* [v=1..M][0..(jmax[v]-jmin[v])] 
			       * hdmin[v][j0] = last position in band on d for state v, and position
			       * j = jmin[v] + j0.*/
  int *safe_hdmin;            /* [1..M] safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0) */
  int *safe_hdmax;            /* [1..M] safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0) */
} CP9Bands_t;



/* used by CM Plan 9 HMM structures */
#define HMMMATCH  0
#define HMMINSERT 1
#define HMMDELETE 2

/* structures from RSEARCH */
#define INIT_RESULTS 100
typedef struct _search_result_node_t {
  int start;
  int stop;
  int bestr;   /* Best root state */
  float score;
  Parsetree_t *tr;
} search_result_node_t;

typedef struct _search_results_t {
  search_result_node_t *data;
  int num_results;
  int num_allocated;
} search_results_t;

typedef struct _dbseq_t {
  ESL_SQ *sq[2];
  search_results_t *results[2];
  int chunks_sent;
  int alignments_sent;           /* -1 is flag for none queued yet */
  float best_score;              /* Best score for scan of this sequence */
  int partition;                 /* For histogram building */
} dbseq_t;

/* sequences to align, for cmalign and cmscore (implemented to ease MPI) */
typedef struct _seqs_to_aln_t {
  ESL_SQ  **sq;                  /* the sequences */
  int nseq;                      /* number of sequences */
  int nalloc;                    /* number of sequences alloc'ed */
  Parsetree_t **tr;              /* parsetrees */
  CP9trace_t **cp9_tr;           /* CP9 traces, usually NULL unless tr is NULL */
  char **postcode;               /* postal codes, left NULL unless do_post */
  float *sc;                     /* score for each seq, can be parsetree score (usually if tr != NULL),
				  * CP9 trace score (usually if cp9_tr != NULL), but could also be
				  * score for the sub parsetree (in case of sub CM alignment)
				  */
} seqs_to_aln_t;

/* The integer log odds score deckpool for integer versions of 
 * Inside and Outside, see cm_postprob.c */
typedef struct Ideckpool_s {
  int   ***pool;
  int      n;
  int      nalloc;
  int      block;
} Ideckpool_t;

/* Structure: Prior_t
 * 
 * Dirichlet priors on all model parameters. 
 */
typedef struct {
  /* transition priors */
  int    tsetnum;                           /* number of transition sets to read in */
  int    tsetmap[UNIQUESTATES][NODETYPES];  /* tsetmap[a][b] is for transition set from ustate a to node b */
  ESL_MIXDCHLET **t;	                    /* array of transition priors, 0..tsetnum-1 */

  /* emission priors */
  ESL_MIXDCHLET *mbp;		/* consensus base pair emission prior */
  ESL_MIXDCHLET *mnt;		/* consensus singlet emission prior */
  ESL_MIXDCHLET *i;		/* nonconsensus singlet emission prior */

  /* bookkeeping */
  int  maxnq;			/* maximum # of components in any prior */
  int  maxnalpha;		/* maximum # of parameters in any prior */
} Prior_t;

#define BUSY 1
#define IDLE 0

/* RSEARCH defaults defined here */
#define DEFAULT_RMATRIX "RIBOSUM85-60"
#define DEFAULT_RALPHA 10.
#define DEFAULT_RBETA 5.
#define DEFAULT_RALPHAP 0.
#define DEFAULT_RBETAP 15.
/* the RSEARCH default is below, it was changed b/c 
 * with no local begin penalty, a glocal hit is ALWAYS
 * going to be decomposed into it's subtrees.
 *#define DEFAULT_RBEGINSC 0.  */
#define DEFAULT_RBEGINSC -0.01
#define DEFAULT_RENDSC -15.

/* The six classes of states in RSEARCH */
#define M_cl 0
#define IL_cl 1
#define DL_cl 2
#define IR_cl 3
#define DR_cl 4
#define DB_cl 5

/* Two modes for padding residues to HMM hits */
#define PAD_SUBI_ADDJ 1
#define PAD_ADDI_SUBJ 2

/* MPI tags */
#define MPI_WORK_EOD    0
#define MPI_WORK_SEARCH 1

#define MPI_RESULTS_SEARCH 2
/* MPI worker time constraints */
#define MPI_WORKER_MIN_SEC 5.
#define MPI_WORKER_MAX_SEC 60.
#define MDPC_SEC 50.

/* RSEARCH macros/#defines etc. (from rnamat.h) */

#define RNAPAIR_ALPHABET "AAAACCCCGGGGUUUU"
#define RNAPAIR_ALPHABET2 "ACGUACGUACGUACGU"
/* Returns true if pos. C of seq B of msa A is a gap */
#define is_rna_gap(A, B, C) (esl_abc_CIsGap(A->abc, A->aseq[B][C]))
/* Returns true if position C of digitized sequence B of msa A is a canonical */
#define is_defined_rna_nucleotide(A, B, C) (esl_abc_CIsCanonical(A->abc, A->aseq[B][C]))
#define unpairedmat_size (matrix_index(3,3) + 1)
#define pairedmat_size (matrix_index (15,15) + 1)
/* Maps to index of matrix, using binary representation of
 * nucleotides (unsorted).
 * See lab book 7, p. 3-4 for details of mapping function (RJK) */
#define matrix_index(X,Y) ((X>Y) ? X*(X+1)/2+Y: Y*(Y+1)/2+X)
/* Matrix type
 * Contains array in one dimension (to be indexed later), matrix size,
 * H, and E. 
 */
typedef struct _matrix_t {
  double *matrix;
  int edge_size;         /* Size of one edge, e.g. 4 for 4x4 matrix */
  int full_size;         /* Num of elements, e.g. 10 for 4x4 matirx */
  double H;
  double E;
} matrix_t;

/* Full matrix definition, includes the g background freq vector (g added by EPN). */
typedef struct _fullmat_t {
  const ESL_ALPHABET *abc;/* alphabet, we enforce it's eslRNA */
  matrix_t *unpaired;
  matrix_t *paired;
  char     *name;
  float    *g;           /* EPN: the background distro, g vector in RSEARCH paper
			  * this now appears in the RIBOSUM matrix files */
  int       scores_flag; /* TRUE if matrix values are log odds scores, FALSE if 
			  * they're target probs, or unfilled */
  int       probs_flag;  /* TRUE if matrix values are target probs, FALSE if 
			  * they're log odds scores, or unfilled */
} fullmat_t;

/* BE_EFFICIENT and BE_PARANOID are alternative (exclusive) settings
 * for the do_full? argument to the alignment engines.
 */
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

/* EPN, Fri Sep  7 16:49:43 2007
 * From HMMER3's p7_config.h:
 *
 * Sean's notes (verbatim):
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * In Forward algorithm implementations, we use a table lookup in
 * p7_FLogsum() to calculate summed probabilities in log
 * space. p7_INTSCALE defines the precision of the calculation; the
 * default of 1000.0 means rounding differences to the nearest 0.001
 * nat. p7_LOGSUM_TBL defines the size of the lookup table; the
 * default of 16000 means entries are calculated for differences of 0
 * to 16.000 nats (when p7_INTSCALE is 1000.0).  e^{-p7_LOGSUM_TBL /
 * p7_INTSCALE} should be on the order of the machine FLT_EPSILON,
 * typically 1.2e-7.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * EPN: Infernal uses bits, not nats. 1.2e-7 =~ 2^-23 =~ e^-16. 
 *      And I've removed the p7_ prefixes.
 */
#if USE_NEWLOGSUM
#define INTSCALE     1000.0f
#define LOGSUM_TBL   23000
#endif
/* Below is infernal -->v0.81 constants for logsums */
#if USE_OLDLOGSUM
#define INTSCALE    1000.0      /* scaling constant for floats to integer scores */
#define LOGSUM_TBL  20000       /* controls precision of ILogsum()            */
#endif


/* from cp9_fastsearch.c, based on HMMER 3's impl_jb.c */
/*****************************************************************
 * 1. CP9_OPROFILE: a scoring profile
 *****************************************************************/

/* Indices for transition scores gm->tsc[k][] */
/* order is optimized for dynamic programming */
enum cp9o_tsc_e { 
  /*  cp9O_MM = 0,
  cp9O_IM = 1, 
  cp9O_DM = 2, 
  cp9O_MD = 3,
  cp9O_DD = 4, 
  cp9O_MI = 5, 
  cp9O_II = 6, 
  cp9O_ID = 7, 
  cp9O_DI = 8, 
  cp9O_BM = 9, 
  cp9O_MEL=10, 
  cp9O_ME =11 */
  cp9O_MM = 0,
  cp9O_IM = 1, 
  cp9O_DM = 2, 
  cp9O_BM = 3, 
  cp9O_MI = 4, 
  cp9O_II = 5, 
  cp9O_DI = 6, 
  cp9O_MD = 7,
  cp9O_ID = 8, 
  cp9O_DD = 9, 
  cp9O_ME =10,
  cp9O_MEL=11 
};
#define cp9O_NTRANS 12

/* Indices for residue emission score vectors
 */
enum cp9o_rsc_e {
  cp9O_MSC = 0,
  cp9O_ISC = 1 
};
#define cp9O_NR 2

typedef struct cp9_oprofile_s {
  int     M;
  int    *tsc;	/* [0.1..M-1][0..cp9X_NTSC-1] */
  int   **rsc;	/* [0..Kp-1][0.1..M][cp9X_NR] */
  const ESL_ALPHABET *abc;
} CP9_OPROFILE;


/*****************************************************************
 * 2. CP9_OMX: a dynamic programming matrix
 *****************************************************************/
enum cp9x_scells_e {
  cp9X_M = 0, 
  cp9X_I = 1,
  cp9X_D = 2, 
  cp9X_EL= 3
};
#define cp9X_NSCELLS 4

typedef struct cp9_omx_s {
  int M;
  int L;

  size_t ncells;
  size_t nrows;
  
  int **dp;			/* [0.1..L][0.1..M][0..cp9X_NSCELLS-1] */
} CP9_OMX;



#endif /*STRUCTSH_INCLUDED*/


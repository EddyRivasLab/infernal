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
#include "p7_config.h"
#include "config.h"

#include "easel.h"
#include "esl_sq.h"
#include "esl_dirichlet.h"
#include "esl_random.h"
#include "esl_sqio.h"

#include "hmmer.h"
#if 0
#include "impl_sse.h"
#endif
#define cmERRBUFSIZE 1024

/* TEMP EPN, Tue Aug 19 17:40:32 2008 */
#define PRINTNOW 0
#define PRINTNOW2 0
/* TEMP */

/* various default parameters for CMs and CP9 HMMs */ 
#define DEFAULT_MIN_CP9_E_CUTOFF 1.0
#define DEFAULT_BETA             0.0000001
#define DEFAULT_TAU              0.0000001
#define DEFAULT_PBEGIN           0.05   /* EPN 06.29.07 (formerly 0.5) */
#define DEFAULT_PEND             0.05   /* EPN 06.29.07 (formerly 0.5) */
#define DEFAULT_ETARGET          0.59   /* EPN 07.10.07 (formerly (v0.7->v0.8)= 2.-0.54 = 1.46 */
#define DEFAULT_NULL2_OMEGA      0.000015258791 /* 1/(2^16), the hard-coded prior probability of the null2 model */
#define DEFAULT_NULL3_OMEGA      0.000015258791 /* 1/(2^16), the hard-coded prior probability of the null3 model */
#define V1P0_NULL2_OMEGA         0.03125        /* 1/(2^5),  the prior probability of the null2 model for infernal versions 0.56 through 1.0.2 */
#define V1P0_NULL3_OMEGA         0.03125        /* 1/(2^5),  the prior probability of the null3 model for infernal versions 0.56 through 1.0.2 */

/* max number of parititons for cmcalibrate */
#define MAX_PARTITIONS 20

/* database size for E-values in Mb for HMM filter thresholds */
#define FTHR_DBSIZE 1000000

/* number of possible integer GC contents, example 40 = 0.40 GC */
#define GC_SEGMENTS 101

#define BE_EFFICIENT  0		/* setting for do_full: small memory mode */
#define BE_PARANOID   1		/* setting for do_full: keep whole matrix, perhaps for debugging */

/* Constants for type of cutoff */
#define SCORE_CUTOFF 0
#define E_CUTOFF     1

/* P7 HMM E-value parameters, borrowed from HMMER, but with 2 new additions: CM_p7_GMU and CM_p7_GLAMBDA */
#define CM_p7_NEVPARAM 8  /* number of statistical parameters stored in cm->p7_evparam */
enum cm_p7_evparams_e {CM_p7_LMMU  = 0, CM_p7_LMLAMBDA = 1, CM_p7_LVMU = 2,  CM_p7_LVLAMBDA = 3, CM_p7_LFTAU = 4, CM_p7_LFLAMBDA = 5, CM_p7_GFMU = 6, CM_p7_GFLAMBDA = 7 };
#define CM_p7_EVPARAM_UNSET -99999.0f  /* if p7_evparam[0] is unset, then all unset                         */

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

/* constants used in cm_pipeline for tallying up number of residues surviving each stage */
#define p7_SURV_F1  0
#define p7_SURV_F1b 1
#define p7_SURV_F2  2
#define p7_SURV_F2b 3
#define p7_SURV_F3  4
#define p7_SURV_F3b 5
#define Np7_SURV    6

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
   * These are created from the probabilities by CP9Logoddsify().
   * By definition, null[] emission scores are all zero.
   * Note that emission distributions are over possible alphabet symbols,
   * not just the unambiguous DNA alphabet: we
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
   * The otsc are reordered transition scores [0..k..M][0..cp9O_NTRANS-1], 
   * for efficiency in DP calculations. This reordering is based on HMMER3. 
   * 
   * CPLAN9_HASBITS flag is up when these scores are valid.
   */
  int  **tsc;                   /* transition scores     [0.9][0.M]       +*/
  int  **msc;                   /* match emission scores [0.Kp-1][1.M] +*/
  int  **isc;                   /* ins emission scores   [0.Kp-1][0.M] +*/
  int   *bsc;                   /* begin transitions     [1.M]              +*/
  int   *esc;			/* end transitions       [1.M]              +*/
  int   *tsc_mem, *msc_mem, *isc_mem, *bsc_mem, *esc_mem;
  int   *otsc;                  /* transition scores [0..M][0..9], special ordering */

  /* The null model probabilities.
   */
  float  *null;                    /* "random sequence" emission prob's     +*/
  float  null2_omega;              /* prior probability of null2 model, copied from CM */
  float  null3_omega;              /* prior probability of null3 model, copied from CM */
  float  p1;                       /* null model loop probability           +*/
  /* local end, EL state parameters */
  float  el_self;                  /* EL transition self loop probability    */
  int    el_selfsc;                /* EL transition self loop score          */
  int   *has_el;                   /* has_el[k] is TRUE if node k has an EL state */
  int   *el_from_ct;               /* [0..M+1] el_from_ct[k] is the number of HMM nodes kp
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
#define CTMEL 3
#define CTIM  4
#define CTII  5
#define CTID  6
#define CTDM  7
#define CTDI  8
#define CTDD  9
#define cp9_NTRANS 10

#define cp9_TRANS_MATCH_OFFSET  0 /* hmm->t[k][0] is first transition out of match */
#define cp9_TRANS_INSERT_OFFSET 4 /* hmm->t[k][4] is first transition out of insert */
#define cp9_TRANS_DELETE_OFFSET 7 /* hmm->t[k][7] is first transition out of delete */
#define cp9_TRANS_NMATCH        4 /* there are 4 transitions out of match */
#define cp9_TRANS_NINSERT       3 /* there are 3 transitions out of insert */
#define cp9_TRANS_NDELETE       3 /* there are 3 transitions out of delete */

/* Declaration of CM Plan9 dynamic programming matrix structure.
 */
typedef struct cp9_mx_s {
  int **mmx;			/* match scores  [0.1..N][0..M] */
  int **imx;			/* insert scores [0.1..N][0..M] */
  int **dmx;			/* delete scores [0.1..N][0..M] */
  int **elmx;			/* end local scores [0.1..N][0..M] */
  int  *erow;                   /* score for E state [0.1..N] */
  /* Hidden ptrs where the real memory is kept */
  int *mmx_mem, *imx_mem, *dmx_mem, *elmx_mem;

  int    M;             /* number of nodes in HMM this mx corresponds to, never changes */
  int    rows;          /* generally L or 2, # of DP rows in seq dimension, where L is length of seq,
			 * == 2 if we're scanning in mem efficient mode, 
			 * never shrinks, but can increase to 'grow' the matrix
			 */
  float  size_Mb;       /* current size of matrix in Megabytes */
  
  /* variables added for HMMER3 p7 HMM banding of CP9 HMM dp algorithms */
  int *kmin;            /* OPTIONAL (can be null) [0.1..i..rows] = k, minimum node for residue i is k */
  int *kmax;            /* OPTIONAL (can be null) [0.1..i..rows] = k, maximum node for residue i is k */
  int  ncells_allocated; /* number of cells allocated in matrix */
  int  ncells_valid;     /* number of cells currently valid in the matrix */

} CP9_MX;

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

/* status flags, cm->flags */
#define CMH_BITS                (1<<0)  /* CM has valid log odds scores             */
#define CMH_ACC                 (1<<1)  /* accession number is available            */
#define CMH_DESC                (1<<2)  /* description exists                       */
#define CMH_GA                  (1<<3)  /* gathering threshold exists               */
#define CMH_TC                  (1<<4)  /* trusted cutoff exists                    */
#define CMH_NC                  (1<<5)  /* noise cutoff exists                      */
#define CMH_LOCAL_BEGIN         (1<<6)  /* Begin distribution is active (local ali) */
#define CMH_LOCAL_END           (1<<7)  /* End distribution is active (local ali)   */
#define CMH_EXPTAIL_STATS       (1<<8)  /* exponential tail stats set               */
#define CMH_FILTER_STATS        (1<<9)  /* filter threshold stats are set           */
#define CMH_QDB                 (1<<10) /* query-dependent bands, QDBs valid        */
#define CMH_CP9                 (1<<11) /* CP9 HMM is valid in cm->cp9              */
#define CMH_CP9STATS            (1<<12) /* CP9 HMM has exp tail stats               */
#define CMH_SCANMATRIX          (1<<13) /* ScanMatrix smx is valid                  */
#define CMH_MLP7                (1<<14) /* 'maximum likelihood' p7 is valid in cm->mlp7 */
#define CMH_MLP7_STATS          (1<<15) /* ml p7 HMM exponential tail stats set    */
#define CMH_AP7                 (1<<16) /* at least 1 additional p7 is valid in cm->ap7 */
#define CMH_AP7_STATS           (1<<17) /* additional p7 HMM exponential tail stats set */

#define CM_IS_SUB               (1<<18) /* the CM is a sub CM                       */
#define CM_IS_RSEARCH           (1<<19) /* the CM was parameterized a la RSEARCH    */
#define CM_RSEARCHTRANS         (1<<20) /* CM has/will have RSEARCH transitions     */
#define CM_RSEARCHEMIT          (1<<21) /* CM has/will have RSEARCH emissions       */
#define CM_EMIT_NO_LOCAL_BEGINS (1<<22) /* emitted parsetrees will never have local begins */
#define CM_EMIT_NO_LOCAL_ENDS   (1<<23) /* emitted parsetrees will never have local ends   */

/* model configuration options, cm->config_opts */
#define CM_CONFIG_LOCAL        (1<<0)  /* configure the model for local alignment  */
#define CM_CONFIG_HMMLOCAL     (1<<1)  /* configure the CP9   for local alignment  */
#define CM_CONFIG_HMMEL        (1<<2)  /* configure the CP9 for EL local aln       */
#define CM_CONFIG_QDB          (1<<3)  /* calculate query dependent bands          */

/* alignment options, cm->align_opts */
#define CM_ALIGN_SMALL         (1<<0)  /* use small CYK D&C                        */
#define CM_ALIGN_QDB           (1<<1)  /* use QD bands                             */
#define CM_ALIGN_HBANDED       (1<<2)  /* use HMM bands                            */
#define CM_ALIGN_SUMS          (1<<3)  /* if using HMM bands, use posterior sums   */
#define CM_ALIGN_SUB           (1<<4)  /* build a sub CM for each seq to align     */
#define CM_ALIGN_HMMVITERBI    (1<<5)  /* use a CP9 HMM only to align, w/viterbi   */
#define CM_ALIGN_INSIDE        (1<<6)  /* use Inside, not CYK                      */
#define CM_ALIGN_POST          (1<<7)  /* do inside/outside and append posteriors  */
#define CM_ALIGN_CHECKINOUT    (1<<8)  /* check inside/outside calculations        */
#define CM_ALIGN_CHECKPARSESC  (1<<9)  /* check parsetree score against aln alg sc */
#define CM_ALIGN_PRINTTREES    (1<<10) /* print parsetrees to stdout               */
#define CM_ALIGN_HMMSAFE       (1<<11) /* realign seqs w/HMM banded CYK bit sc < 0 */
#define CM_ALIGN_SCOREONLY     (1<<12) /* do full CYK/inside to get score only     */
#define CM_ALIGN_SAMPLE        (1<<13) /* sample parsetrees from the inside matrix */
#define CM_ALIGN_FLUSHINSERTS  (1<<14) /* flush inserts L/R like pre 1.0 infernal  */
#define CM_ALIGN_CHECKFB       (1<<15) /* check forward/backward CP9 HMM calcs     */
#define CM_ALIGN_OPTACC        (1<<16) /* no CYK, aln w/Holmes/Durbin opt accuracy */
#define CM_ALIGN_HMM2IJOLD     (1<<17) /* use old hmm2ij band calculation alg      */
#define CM_ALIGN_P7BANDED      (1<<18) /* use p7 HMM bands to band the CP9 HMM     */

/* search options, cm->search_opts */
#define CM_SEARCH_NOQDB        (1<<0)  /* DO NOT use QDB to search (QDB is default)*/
#define CM_SEARCH_HBANDED      (1<<1)  /* use HMM bands for search                 */
#define CM_SEARCH_HMMALNBANDS  (1<<2)  /* force full aln when deriving HMM bands   */
#define CM_SEARCH_SUMS         (1<<3)  /* if using HMM bands, use posterior sums   */
#define CM_SEARCH_INSIDE       (1<<4)  /* scan with Inside, not CYK                */
#define CM_SEARCH_NOALIGN      (1<<5)  /* don't align hits, just report locations  */
#define CM_SEARCH_RSEARCH      (1<<6)  /* use RSEARCH parameterized CM             */
#define CM_SEARCH_CMGREEDY     (1<<7)  /* use greedy alg to resolve CM overlaps    */
#define CM_SEARCH_HMMGREEDY    (1<<8)  /* use greedy alg to resolve HMM overlaps   */
#define CM_SEARCH_HMMVITERBI   (1<<9)  /* search with CP9 HMM Viterbi              */
#define CM_SEARCH_HMMFORWARD   (1<<10) /* search with CP9 HMM Forward              */
#define CM_SEARCH_HMM2IJOLD    (1<<11) /* use old hmm2ij band calculation alg      */
#define CM_SEARCH_NULL2        (1<<12) /* use NULL2 score correction               */
#define CM_SEARCH_NULL3        (1<<13) /* use NULL3 score correction               */

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
  char *rf;             /* reference annotation line (NULL if unavail) */
  char *cstr;		/* CM consensus structure line                 */
  char *cseq;		/* CM consensus sequence line                  */
  char *mid;		/* alignment identity middle line              */
  char *top;		/* optional, non-compensatory 'x' top line     */
  char *aseq;		/* aligned target sequence                     */
  char *pcode;          /* aligned posteriors 'ones' place (9 in 93)   */
  int  *scoord;		/* coords 1..L for aligned dsq chars           */
  int  *ccoord;	        /* coords 1..clen for aligned consensus chars  */
  int   len;		/* len of the strings above                    */
  int   cfrom, cto;	/* max bounds in ccoord                        */
  int   sqfrom, sqto;	/* max bounds in scoord                        */

  char *hmmname;		/* name of HMM                          */
  char *hmmacc;			/* accession of HMM; or [0]='\0'        */
  char *hmmdesc;		/* description of HMM; or [0]='\0'      */
  
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
  int     *isum_pn_m;         /* [0..k..hmm_M] sum over i of log post probs from post->mmx[i][k]*/
  int     *isum_pn_i;         /* [0..k..hmm_M] sum over i of log post probs from post->imx[i][k]*/
  int     *isum_pn_d;         /* [0..k..hmm_M] sum over i of log post probs from post->dmx[i][k]*/

  /* arrays for CM state bands, derived from HMM bands */
  int *imin;                  /* [0..cm_M-1] imin[v] = first position in band on i for state v*/
  int *imax;                  /* [0..cm_M-1] imax[v] = last position in band on i for state v*/
  int *jmin;                  /* [0..cm_M-1] jmin[v] = first position in band on j for state v*/
  int *jmax;                  /* [0..cm_M-1] jmax[v] = last position in band on j for state v*/
  int **hdmin;                /* [0..cm_M-1][0..(jmax[v]-jmin[v])] 
			       * hdmin[v][j0] = first position in band on d for state v, and position
			       * j = jmin[v] + j0. */
  int **hdmax;                /* [0..cm_M-1][0..(jmax[v]-jmin[v])] 
			       * hdmin[v][j0] = last position in band on d for state v, and position
			       * j = jmin[v] + j0.*/
  int *hdmin_mem;             /* actual memory for hdmin */
  int *hdmax_mem;             /* actual memory for hdmax */
  int *safe_hdmin;            /* [0..cm_M-1] safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0) */
  int *safe_hdmax;            /* [0..cm_M-1] safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0) */

  /* info on size of bands */
  int hd_needed;              /* Sum_v cp9b->jmax[v] - cp9b->jmin[v] + 1, number of hd arrays needed */
  int hd_alloced;             /* number of hd arrays currently alloc'ed */

} CP9Bands_t;

/* used by CM Plan 9 HMM structures */
#define HMMMATCH  0
#define HMMINSERT 1
#define HMMDELETE 2
#define NHMMSTATETYPES 3

/* structures from RSEARCH */
#define INIT_RESULTS 100
typedef struct _search_result_node_t {
  int start;
  int stop;
  int bestr;   /* Best root state */
  float score;
  Parsetree_t *tr;
  char *pcode;           /* posterior code string, left NULL unless cm->search_opts & CM_SEARCH_POST */
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
  char **postcode;               /* posterior code string, sometimes NULL */
  float *sc;                     /* score for each seq, can be parsetree score (usually if tr != NULL),
				  * CP9 trace score (usually if cp9_tr != NULL), but could also be
				  * score for the sub parsetree (in case of sub CM alignment)
				  */
  float *pp;                     /* average posterior probability for each seq, if applicable, IMPOSSIBLE if not */
  float *struct_sc;              /* contribution of (MATP emission scores minus marginalized scores) for each tr */ 
} seqs_to_aln_t;

struct deckpool_s {
  float ***pool;
  int      n;
  int      nalloc;
  int      block;
};

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

#define MPI_MIN_CHUNK_W_MULTIPLIER 10
#define MPI_MAX_CHUNK_SIZE 1000000

/* RSEARCH macros/#defines etc. (from rnamat.h) */

/* coordinate -- macro that checks if it's reverse complement and if so 
   returns coordinate in original strand
   a = true if revcomp, false if not
   b = the position in current seq
   c = length of the seq
*/
#define COORDINATE(a,b,c) ( a ? -1*b+c+1 : b)
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
#if 1
#define INTSCALE     1000.0f
#define LOGSUM_TBL   23000
#endif
#if 0
#define INTSCALE    1000.0      /* scaling constant for floats to integer scores */
#define LOGSUM_TBL  20000       /* controls precision of ILogsum()            */
#endif


/* from cp9_fastsearch.c, based on HMMER 3's impl_jb.c */

/* Indices for transition scores tsc[k][] */
/* order is optimized for dynamic programming */
enum cp9o_tsc_e { 
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

enum cp9_locality_e {
  CP9_LOCAL_BEGIN_END_OFF_AND_EL_OFF = 0,
  CP9_LOCAL_BEGIN_END_OFF_AND_EL_ON  = 1,
  CP9_LOCAL_BEGIN_END_ON_AND_EL_OFF  = 2,
  CP9_LOCAL_BEGIN_END_ON_AND_EL_ON   = 3
};
#define nCP9_LOCALITIES 4


enum cm_locality_e {
  CM_LOCAL_MODE = 0,
  CM_GLOCAL_MODE = 1
};
#define nCM_LOCALITIES 2

enum emitmode_e {
  EMITLEFT  = 0,
  EMITRIGHT = 1,
  EMITPAIR  = 2,
  EMITNONE  = 3
};
#define nEMITMODES 4 

/* Declaration of CM dynamic programming matrix structure for 
 * alignment with float scores in vjd (state idx, aln posn,
 * subseq len) coordinates. May be banded in j and/or d dimensions.
 */
typedef struct cm_hb_mx_s {
  int  M;		/* number of states (1st dim ptrs) in current mx */
  int  L;               /* length of sequence the matrix currently corresponds to */

  int64_t    ncells_alloc;	/* current cell allocation limit */
  int64_t    ncells_valid;	/* current number of valid cells */
  float      size_Mb;       /* current size of matrix in Megabytes */

  int   *nrowsA;        /* [0..v..M] current number allocated rows for deck v */

  float ***dp;          /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
  float   *dp_mem;      /* the actual mem, points to dp[0][0][0] */

  CP9Bands_t *cp9b;     /* the CP9Bands_t object associated with this
			 * matrix, which defines j, d, bands for each
			 * state, only a reference, so don't free
			 * it when mx is freed. */
} CM_HB_MX;


/* Declaration of CM dynamic programming matrix structure for 
 * alignment with float scores in vjd (state idx, aln posn,
 * subseq len) coordinates. May be banded in j and/or d dimensions.
 */
typedef struct cm_hb_shadow_mx_s {
  int  M;		/* number of states (1st dim ptrs) in current mx */
  int  L;               /* length of sequence the matrix currently corresponds to */

  int    y_ncells_alloc;  /* current cell allocation limit in yshadow*/
  int    y_ncells_valid;  /* current number of valid cells in yshadow */
  int    k_ncells_alloc;  /* current cell allocation limit in kshadow*/
  int    k_ncells_valid;  /* current number of valid cells in kshadow */
  float  size_Mb;         /* current size of matrix in Megabytes */

  int   *nrowsA;          /* [0..v..M] current number allocated rows for deck v */
  int    nbifs;           /* number of B_st's in the cm this mx was created for */

  /* yshadow holds the shadow matrix for all non-BIF_B states, yshadow[v] == NULL if cm->sttype[v] == B_st */
  char ***yshadow;       /* [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
  char   *yshadow_mem;   /* the actual mem, points to yshadow[0][0][0] */

  /* kshadow holds the shadow matrix for all BIF_B states, kshadow[v] == NULL if cm->sttype[v] != B_st */
  int ***kshadow;       /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
  int   *kshadow_mem;   /* the actual mem, points to kshadow[0][0][0] */

  CP9Bands_t *cp9b;     /* the CP9Bands_t object associated with this
			 * matrix, which defines j, d, bands for each
			 * state, only a reference, so don't free
			 * it when mx is freed. */
} CM_HB_SHADOW_MX;

/* Structure ScanMatrix_t: Information used by all CYK/Inside scanning functions,
 * compiled together into one data structure for convenience. 
 */
#define cmSMX_HAS_FLOAT (1 << 0)  /* if float versions of alpha and precalc'ed scores are valid */
#define cmSMX_HAS_INT   (1 << 1)  /* if int   versions of alpha and precalc'ed scores are valid */
typedef struct scanmx_s {
  /* general info about the model/search */
  int     cm_M;        /* # states in the CM */
  int     W;           /* max hit size */
  int    *dmin;        /* [0..v..cm->M-1] min subtree length for v using beta, just a ref, NULL for non-banded */
  int    *dmax;        /* [0..v..cm->M-1] max subtree length for v using beta, just a ref, NULL for non-banded */
  int   **dnAA;        /* [1..j..W][0..v..M-1] max d value allowed for posn j, state v */
  int   **dxAA;        /* [1..j..W][0..v..M-1] max d value allowed for posn j, state v */
  int    *bestr;       /* auxil info: best root state at alpha[0][cur][d] */
  int     flags;       /* flags for what info has been set (can be float and/or int versions of alpha) */
  double  beta_qdb;    /* tail loss prob used for calc'ing dmin/dmax, invalid if dmin==dmax==NULL */
  double  beta_W;      /* tail loss prob used for calc'ing W, often == beta_qdb, may be greater, can't be less */

  /* falpha dp matrices [0..j..1][0..v..cm->M-1][0..d..W] for float implementations of CYK/Inside */
  float ***falpha;          /* non-BEGL_S states for float versions of CYK/Inside */
  float ***falpha_begl;     /*     BEGL_S states for float versions of CYK/Inside */
  float   *falpha_mem;      /* ptr to the actual memory for falpha */
  float   *falpha_begl_mem; /* ptr to the actual memory for falpha_begl */

  /* ialpha dp matrices [0..j..1][0..v..cm->M-1][0..d..W] for integer implementations of CYK/Inside */
  int   ***ialpha;          /* non-BEGL_S states for int   versions of CYK/Inside */
  int   ***ialpha_begl;     /*     BEGL_S states for int   versions of CYK/Inside */
  int     *ialpha_mem;      /* ptr to the actual memory for ialpha */
  int     *ialpha_begl_mem; /* ptr to the actual memory for ialpha_begl */

  int      ncells_alpha;      /* number of alloc'ed, valid cells for falpha and ialpha matrices, alloc'ed as contiguous block */
  int      ncells_alpha_begl; /* number of alloc'ed, valid cells for falpha_begl and ialpha_begl matrices, alloc'ed as contiguous block */
} ScanMatrix_t;

/* Structure GammaHitMx_t: gamma semi-HMM used for optimal hit resolution
 * of a CM or CP9 scan. All arrays are 0..L.
 */
typedef struct gammahitmx_s {
  int       L;                  /* length of sequence, arrays are size L+1 */
  float    *mx;                 /* [0..L] SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* [0..L] traceback pointers for SHMM */ 
  float    *savesc;             /* [0..L] saves score of hit added to best parse at j */
  int      *saver;		/* [0..L] saves initial non-ROOT state of best parse ended at j */
  float     cutoff;             /* minimum score to report */
  int       i0;                 /* position of first residue in sequence (gamma->mx[0] corresponds to this residue) */
  int       iamgreedy;          /* TRUE to use RSEARCH's greedy overlap resolution alg, FALSE to use optimal alg */
} GammaHitMx_t;

/* Structure Theta_t: probability a parsetree of score <= x will be emitted from
 *                    the subtree rooted at v.           
 * of a CM scan. All arrays are 0..L.
 */
typedef struct theta_s {
  int       L;                  /* length of sequence, arrays are size L+1 (or L+2 if iambackward = TRUE) */
  float    *mx;                 /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  float     cutoff;             /* minimum score to report */
  int       i0;                 /* position of first residue in sequence (gamma->mx[0] corresponds to this residue) */
} Theta_t;

/* Structure SearchInfo_t: 
 * 
 * Information for CM searches, including info on filters.  
 * <nrounds> holds number of rounds of filtering.  
 * <search_opts>, <sc_cutoff>, <e_cutoff> <stype>, and <hsi>
 * are all arrays of length <nrounds + 1>, running [0..nrounds].  
 * The final value in all those arrays (index <nrounds>) corresponds to
 * the final scan, when filtering is finished.  A special case is when
 * <nrounds> == 0, in this case we're not filtering.
 *
 * A note about the mandatory use of two cutoffs <sc_cutoff> and <e_cutoff>
 * If <cutoff_type> == E_CUTOFF, <sc_cutoff> is the minimal bit score
 * that satisfies the E_CUTOFF for all partitions, this is used for
 * reporting hits in SearchDispath(), but the final cutoff used is still
 * <e_cutoff>.
 * If <cutoff_type> == SCORE_CUTOFF, <e_cutoff> is set to -1., and never
 * used. It should be considered invalid.
 *
 */                                                                                                      
typedef struct searchinfo_s {
  int    nrounds;            /* number of rounds of filtering, if 0, we're not filtering */
  int   *stype;              /* [0..n..nrounds] search 'type' "SEARCH_WITH_HMM", "SEARCH_WITH_CM" */
  int   *search_opts;        /* [0..n..nrounds] search options for each round of filtering, including the final round */
  int   *cutoff_type;        /* [0..n..nrounds] SCORE_CUTOFF or E_CUTOFF */
  float *sc_cutoff;          /* [0..n..nrounds] bit score cutoff threshold for each round, always valid, 
			      * if cutoff_type[n] == E_CUTOFF this is minimal bit score across all partitions for e_cutoff */
  float *e_cutoff;           /* [0..n..nrounds] E-value cutoff threshold for each round, ONLY valid if if cutoff_type[n] == E_CUTOFF */
  ScanMatrix_t     **smx;    /* [0..n..nrounds] scanning DP matrix for each round, for final round (n==nrounds) si->smx[nrounds] == cm->smx */
} SearchInfo_t;

/* possible values for stype[] array in SearchInfo_t objects */
#define SEARCH_WITH_HMM    0  
#define SEARCH_WITH_CM     1


/* Structure ExpInfo_t:
 *
 * Info on an exponential tail that describes score distribution in random sequence 
 * of a given algorithm, model configuration (can be 1 of 8 modes, local/glocal of
 * each  viterbi, forward, cyk, or inside). Fit in cmcalibrate and stored in the
 * cm file. All values with the sole exception of <cur_eff_dbsize> are never 
 * changed once read, and some of the params are actually unnecessary downstream 
 * of cmcalibrate, but are potentially informative to the user.
 */
typedef struct expinfo_s {
  long   cur_eff_dbsize;/* the total number of possible hits we expect for current database search, 
			 * cur_eff_dbsize = ceiling((current dbsize) / <dbsize> * <nrandhits>) */
  double lambda;	/* scale param exponential tail */
  double mu_extrap;	/* offset/location param for exponential tail extrapolated to include all <nrandhits> from cmcalibrate, 
			 * mu_corrected = mu_orig - log(1./tailp) / lambda */
  double mu_orig;	/* offset/location param for exponential tail's original fit to tailp of rand seq score histogram in cmcalibrate */
  long   dbsize;        /* db size in residues that was used in cmcalibrate */
  int    nrandhits;     /* total number of hits in random sequence database in cmcalibrate */ 
  double tailp;         /* fractional tail mass threshold for hit histogram in random sequences in cmcalibrate */
  int    is_valid;      /* TRUE if this expinfo_s object is valid (it's parameters have been set), FALSE if not */
} ExpInfo_t;

/* Structure HMMFilterInfo_t: 
 * 
 * Information for HMM filters of CM searches as determined in cmcalibrate. 
 */                                                                                                      
typedef struct hmmfilterinfo_s {
  int    is_valid;                /* TRUE if values have been set, FALSE if not */
  float  F;                       /* fraction of empirical CM hits required to survive filter */
  int    N;                       /* number of CM hits used to get threshold ((N*F) passed)*/
  long   dbsize;                  /* db size used to calculate E-values in fwd_E_cut and cm_E_cut, SHOULD ALWAYS BE FTHR_DBSIZE_MB 1 */
  int    ncut;                    /* number of filter threshold cutoff pairs we have, size of fwd_E_cut and cm_E_cut arrays */
  float *cm_E_cut;                /* [0..i..ncut-1] CM E-value cutoff used to determine fwd_E_cut[i], fwd_E_cut[i], 
				   * at least F fraction of hits with optimal hits with CM E-value < cm_E_cut[i] were able
				   * to be recognized with by a forward HMM scan with E-value cutoff fwd_E_cut[i].
				   * These are sorted in decreasing order, from worst, highest E-value to best, lowest.
				   */
  float *fwd_E_cut;               /* [0..i..ncut-1] cutoff E-value threshold for HMM forward filter, using this cutoff
				   * we were able to find at least F fraction of CYK hits with E-value of cm_E_cut[i] or better.
				   * (we can use this E-value and db_size and exponential tailto get bit score for each partition) 
				   */
  int    always_better_than_Smax; /* tells us what we should do if given a CM E-value cutoff worse (higher) than cm_E_cut[0]:
				   * If TRUE  we should use fwd_E_cut[0] as the HMM filter cutoff. In this case cm_E_cut[0] was
				   *          the worst CM E-value we observed. And we could still recognize F fraction of the hits
				   *          with an expected survival fraction < Smax.
				   * If FALSE we should turn HMM filtering off in this scenario. In this case cm_E_cut[0] was
				   *          not the worst E-value we observed, so there were some for which we couldn't
				   *          filter effectively (cm_E_cut[0] is the worst E-value cutoff we did observe for
				   *          which we can filter with expected survival fraction < Smax.
				   */
} HMMFilterInfo_t;

/* different possible x/y values for HMM filter threshold plots in cmstat */
#define FTHR_PLOT_CME_HMME  0 /* HMM filter E-value cutoffs versus CM E-value cutoffs */
#define FTHR_PLOT_CME_S     1 /* predicted survival fraction versus CM E-value cutoffs */
#define FTHR_PLOT_CME_XHMM  2 /* predicted xhmm (factor slower than HMM only scan) versus CM E-value cutoffs */
#define FTHR_PLOT_CME_SPDUP 3 /* predicted speedup with filter versus CM E-value cutoffs */
#define FTHR_PLOT_CMB_HMMB  4 /* HMM filter bit score cutoffs versus CM bit score cutoffs */
#define FTHR_PLOT_CMB_S     5 /* predicted survival fraction versus CM bit score cutoffs */
#define FTHR_PLOT_CMB_XHMM  6 /* predicted xhmm (factor slower than HMM only scan) versus CM bit score cutoffs */
#define FTHR_PLOT_CMB_SPDUP 7 /* predicted speedup with filter versus CM bit score cutoffs */
#define FTHR_NPLOT          8

/* Structure BestFilterInfo_t: 
 * 
 * Information for the predicted best filter for CM searches
 * as determined in cmcalibrate (only used in cmcalibrate).
 */                                                                                                      
typedef struct bestfilterinfo_s {
  int           cm_M;                /* # states in the CM */
  float         cm_eval;             /* CM E-value threshold, we rejected worse than   */
  float         F;                   /* fraction of empirical CM hits required to survive filter */
  int           N;                   /* number of CM hits used to get threshold ((N*F) passed)*/
  int           db_size;             /* db size used to calculate exponential mu for *_eval calculations */
  int           is_valid;            /* TRUE if values have been set, FALSE if not */
  int           ftype;               /* FILTER_WITH_HMM_VITERBI, FILTER_WITH_HMM_FORWARD or FILTER_NOTYETSET */
  float         e_cutoff;            /* cutoff E-value threshold for filter (we can use this and db_size and exponential tail to get bit score for each partition) */
  float         full_cm_ncalcs;      /* millions of DP calcs for full CM scan of length db_size */
  float         fil_ncalcs;          /* millions of DP calcs for filter scan of length db_size */
  float         fil_plus_surv_ncalcs;/* millions of DP calcs for filter scan + full CM scan of survivors of length db_size */
} BestFilterInfo_t;



/* possible values for ftype[] array in FilterInfo_t objects */
#define FILTER_WITH_HMM_VITERBI 0  
#define FILTER_WITH_HMM_FORWARD 1  
#define FILTER_NOTYETSET        2

/* Structure CMStats_t
 */
typedef struct cmstats_s {
  int np;                    /* number of partitions, default: 1 */
  int *ps;                   /* start GC content [0..100] of each partition */
  int *pe;                   /* end   GC content [0..100] of each partition */
  int gc2p[GC_SEGMENTS];     /* map from GC content to partition number     */
  ExpInfo_t ***expAA;        /* [0..EXP_NMODES-1][0..np-1] */
  HMMFilterInfo_t **hfiA;    /* [0..FTHR_NMODES-1] */
} CMStats_t;


/* Exponential tail statistics modes, a different exp tail fit exists for each mode
 * 0..EXP_NMODES-1 are first dimension of cmstats->expAA 
 * order is important, it's exploited by cmcalibrate 
 */
/* 
As of 1.0.2:
*/
#define EXP_CP9_GV 0
#define EXP_CP9_GF 1
#define EXP_CM_GC  2  
#define EXP_CM_GI  3
#define EXP_CP9_LV 4
#define EXP_CP9_LF 5
#define EXP_CM_LC  6
#define EXP_CM_LI  7
#define EXP_NMODES 8

/* Filter threshold modes, used in cmcalibrate 
 * 0..FTHR_MODES-1 are only dimension cmstats->fthrA 
 */
#define FTHR_CM_GC 0
#define FTHR_CM_GI 1
#define FTHR_CM_LC 2
#define FTHR_CM_LI 3
#define FTHR_NMODES 4

/* Structure ComLog_t: command line info used to build/calibrate a CM.
 * 
 * bcom, bdate must be non-NULL. 
 * ccom1, cdate1 is non-NULL only if at least 1 cmcalibrate call was performed for this cm
 * ccom2, cdate2 is non-NULL only if at > 1 cmcalibrate calls were performed for this cm AND
 * the most recent cmcalibrate call had --filonly enabled (meaning only filter thresholds were rewritten).
 * 
 */
typedef struct comlog_s {
  char     *bcom;           /* command line used for cmbuild, if --gibbs used w/o --seed, --seed will be artificially appended */
  char     *bdate;          /* date of cmbuild call */
  char     *ccom;           /* command line used for first of up to two cmcalibrate calls, if -s not used, -s will be artificially appended */
  char     *cdate;          /* date of first of up to two cmcalibrate call */
} ComLog_t;

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
  char *rf;             /*   consensus column annotation line, or NULL     */ /* ONLY PARTIALLY IMPLEMENTED, BEWARE */

  /* new as of v1.0 */
  ComLog_t *comlog;	/*   creation dates and command line(s) that built/calibrated the model (mandatory) */
  int    nseq;		/*   number of training sequences          (mandatory) */
  float  eff_nseq;	/*   effective number of seqs (<= nseq)    (mandatory) */
  float  ga;	        /*   per-seq gathering thresholds (bits) (CMH_GA) */
  float  tc;            /*   per-seq trusted cutoff (bits)       (CMH_TC) */
  float  nc;	        /*   per-seq noise cutoff (bits)         (CMH_NC) */


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
  float **oesc;         /*   Optimized emission score log odds float vec   */
  float *beginsc;	/*   Score for ROOT_S -> state v (local alignment) */
  float *endsc;   	/*   Score for state_v -> EL (local alignment)     */

                        /* Parameters used in marginal alignments          */
  float **lmesc;        /*   Left marginal emission scores (log odds)      */
  float **rmesc;        /*   Right marginal emission scores (log odds)     */                

			/* Scaled int parameters of the log odds model:    */
  int  **itsc;		/*   Transition score vector, scaled log odds int  */
  int  **iesc;		/*   Emission score vector, scaled log odds int    */
  int  **ioesc;         /*   Optimized emission score log odds int vector  */
  int   *ibeginsc;      /*   Score for ROOT_S -> state v (local alignment) */
  int   *iendsc;  	/*   Score for state_v -> EL (local alignment)     */

  float   null2_omega;  /* prior probability of the null2 model (if it is used) */
  float   null3_omega;  /* prior probability of the null3 model (if it is used) */

  int    flags;		/* status flags                                    */

  /* W and query dependent bands (QDB) on subsequence lengths at each state */
  int   *dmin;          /* minimum d bound for each state v; [0..v..M-1] (NULL if non-banded) */
  int   *dmax;          /* maximum d bound for each state v; [0..v..M-1] (NULL if non-banded) */
  int    W;             /* max d: max size of a hit (EPN 08.18.05)                            */
  double beta_qdb;      /* tail loss probability for QDB calculation used to set dmin/dmax    */
  double beta_W;        /* tail loss probability for QDB calculation used to set W, often     *
			 * equal to beta_qdb, but not always. beta_W >= beta_qdb ALWAYS.      *
			 * If beta_W > beta_qdb, dmax[0] > W, d values > W are not allowed    *
			 * in the DP algorithms though (enforced sneakily when the            *
			 * ScanMatrix_t is built). However, if beta_W > beta_qdb, we still    *
			 * can get less sensitivity loss w.r.t non-banded than if bands were  * 
			 * tighter with beta_W == beta_qdb; because some subtrees (think      * 
			 * BEGL's and BEGRs) still have wider bands, it's just the nodes near * 
			 * the root that will have their dmax values truncated to <= cm->W.   */
  double tau;           /* tail loss probability for HMM target dependent banding             */

  /* added by EPN, Tue Jan  2 14:24:08 2007 */
  int        config_opts;/* model configuration options                                        */
  int        align_opts; /* alignment options                                                  */
  int        search_opts;/* search options                                                     */
  CP9_t     *cp9;        /* a CM Plan 9 HMM, always built when the model is read from a file   */
  CP9Map_t  *cp9map;     /* the map from the Plan 9 HMM to the CM and vice versa               */
  CP9Bands_t *cp9b;      /* the CP9 bands                                                      */
  float     *root_trans; /* transition probs from state 0, saved IFF zeroed in ConfigLocal()   */
  float      pbegin;     /* local begin prob to spread across internal nodes for local mode    */
  float      pend;       /* local end prob to spread across internal nodes for local mode      */
  
  float  el_selfsc;     /* score of a self transition in the EL state
			 * the EL state emits only on self transition (EPN 11.15.05)*/
  int   iel_selfsc;     /* scaled int version of el_selfsc         */

  /* DP matrices and some auxiliary info for DP algorithms */
  ScanMatrix_t *smx;     /* matrices, info for CYK/Inside scans with this CM */
  CM_HB_MX     *hbmx;    /* growable HMM banded float matrix */
  CP9_MX       *cp9_mx;  /* growable CP9 DP matrix */
  CP9_MX       *cp9_bmx; /* another growable CP9 DP matrix, 'b' is for backward,
			  * only alloc'ed to any significant size if we do Forward,Backward->Posteriors */

  /* search info describing the cmsearch filtering strategy, NULL unless created in cmsearch */
  SearchInfo_t *si;      /* describes each round of filtering, and final round of searching */

  /* statistics */
  CMStats_t *stats;      /* holds exponential tail stats and HMM filtering thresholds */

  /* p7 hmms, added 08.05.08 */
  P7_HMM       *mlp7;         /* the maximum likelihood p7 HMM, built from the CM  */
  float         mlp7_evparam[CM_p7_NEVPARAM]; /* E-value params (CMH_MLP7_STATS) */

  /* p7 hmms, added 08.05.08 */
  int          nap7;           /* number of additional p7 HMMs read from file */
  P7_HMM      **ap7A;          /* query p7 HMM */
  float       **ap7_evparamAA; /* E-value params (CMH_AP7_STATS) */

  const  ESL_ALPHABET *abc; /* ptr to alphabet info (cm->abc->K is alphabet size)*/
} CM_t;


enum cm_pipemodes_e { CM_SEARCH_SEQS = 0, CM_SCAN_MODELS = 1 };
enum cm_zsetby_e    { CM_ZSETBY_DBSIZE = 0, CM_ZSETBY_OPTION = 1, CM_ZSETBY_FILEINFO = 2 };

typedef struct cm_pipeline_s {
  /* Dynamic programming matrices                                           */
  P7_OMX       *oxf;		/* one-row Forward matrix, accel pipe       */
  P7_OMX       *oxb;		/* one-row Backward matrix, accel pipe      */
  P7_OMX       *fwd;		/* full Fwd matrix for envelopes            */
  P7_OMX       *bck;		/* full Bck matrix for envelopes            */
  P7_GMX       *gxf;		/* generic Forward matrix                   */
  P7_GMX       *gxb;		/* generic Backward matrix                  */
  P7_GMX       *gfwd;		/* generic full Fwd matrix for envelopes    */
  P7_GMX       *gbck;		/* generic full Bck matrix for envelopes    */
  ScanMatrix_t *fsmx;           /* scan matrix for CYK filter stage         */
  ScanMatrix_t *smx;            /* scan matrix for final stage              */

  /* Domain/envelope postprocessing                                         */
  ESL_RANDOMNESS *r;		/* random number generator                  */
  int             do_reseeding; /* TRUE: reseed for reproducible results    */
  P7_DOMAINDEF   *ddef;		/* domain definition workflow               */

  /* Reporting threshold settings                                           */
  int     by_E;		        /* TRUE to cut per-target report off by E   */
  double  E;	                /* per-target E-value threshold             */
  double  T;	                /* per-target bit score threshold           */
  int     use_bit_cutoffs;      /* (FALSE | CMH_GA | CMH_TC | CMH_NC)       */

  /* Inclusion threshold settings                                           */
  int     inc_by_E;		/* TRUE to threshold inclusion by E-values  */
  double  incE;			/* per-target inclusion E-value threshold   */
  double  incT;			/* per-target inclusion score threshold     */

  /* Tracking search space sizes for E value calculations                   */
  double  Z;			/* eff # targs searched (per-target E-val)  */
  enum cm_zsetby_e Z_setby;   	/* how Z was set                            */
  
  /* Threshold settings for pipeline                                        */
  int     do_max;	        /* TRUE to run in slow/max mode             */
  int     do_rfam;	        /* TRUE to run in Rfam pipeline (fast) mode */
  double  F1;		        /* MSV filter threshold                     */
  double  F2;		        /* Viterbi filter threshold                 */
  double  F3;		        /* uncorrected Forward filter threshold     */
  double  F4;		        /* glocal Forward filter thr                */
  double  F5;		        /* glocal env def filter thr                */
  double  F6;		        /* CYK filter thr                           */
  double  F1b;		        /* bias-corrected MSV filter threshold      */
  double  F2b;		        /* bias-corrected Viterbi filter threshold  */
  double  F3b;		        /* bias-corrected Forward filter threshold  */
  double  F4b;		        /* bias-corrected gloc Forward filter threshold  */
  double  F5b;		        /* bias-corrected env def filter threshold  */
  double  orig_F4;		/* glocal Forward filter thr                */
  double  orig_F4b;		/* bias-corrected glocal Fwd threshold      */
  double  orig_F5;		/* per-envelope Forward filter thr          */
  double  orig_F5b;	        /* bias-corrected per-envelope threshold    */
  int     do_cykenv;	        /* TRUE to redefine envelopes after CYK     */
  double  F6env;	        /* CYK envelope P-value threshold           */
  int     do_F3env;	        /* TRUE to redefine envelopes after local fwd */
  int     do_F4env;	        /* TRUE to redefine envelopes after glocal fwd */
  int     do_F4env_strict;      /* TRUE to be strict with glocal fwd env redefn */
  double  F3envX;               /* max avg number of passes through model for local fwd env redefn */
  double  F4envX;               /* max avg number of passes through model for glocal fwd env redefn */
  int     F3ns;                 /* number of samples for local fwd env redfn */
  int     F4ns;                 /* number of samples for glocal fwd env redfn */
  int     do_cm;		/* TRUE to use CM for at least one stage    */
  int     do_hmm;		/* TRUE to use HMM for at least one stage   */
  int     do_envelopes;		/* TRUE to find envelopes in windows prior to CM stages */
  int     do_pad;		/* TRUE to pad hits based on cm->W          */
  int     do_msvmerge;		/* TRUE to merge MSV hits, FALSE not to     */
  int     do_msv;		/* TRUE to filter with MSV, FALSE not to    */
  int     do_shortmsv;		/* TRUE to filter with standard MSV, not Longtarget variant */
  int     do_msvbias;	        /* TRUE to use biased comp HMM filter w/MSV */
  int     do_vitbias;      	/* TRUE to use biased comp HMM filter w/Vit */
  int     do_fwdbias;     	/* TRUE to use biased comp HMM filter w/Fwd */
  int     do_gfwdbias;     	/* TRUE to use biased comp HMM filter w/gFwd*/
  int     do_envbias;     	/* TRUE to use biased comp HMM filter w/ddef*/
  int     do_vit;		/* TRUE to filter with Vit, FALSE not to    */
  int     do_fwd;		/* TRUE to filter with Fwd, FALSE not to    */
  int     do_gfwd;		/* TRUE to filter w/glocal Fwd, FALSE not to*/
  int     do_cyk;		/* TRUE to filter with CYK, FALSE not to    */
  int     do_null2;		/* TRUE to use null2 score corrections      */
  int     do_null3;		/* TRUE to use null3 score corrections      */
  int     do_localenv;          /* TRUE to define envelopes in local mode   */
  int     do_wsplit;            /* TRUE to split MSV windows > wmult*W      */
  int     do_wcorr;             /* TRUE to correct for window size          */
  double  wmult;                /* scalar * W, for do_wsplit                */
  int     do_oldcorr;           /* TRUE to use old correction for env def   */
  int     do_nocorr;            /* TRUE to use no correction for env def    */
  int     do_envwinbias;        /* TRUE to calc env bias for entire window  */
  int     do_fwdbias_sampling;  /* TRUE to calculate Fwd bias (F3b) based on sampled traces */
  int     do_gmsv;              /* TRUE to use generic MSV */
  int     do_filcmW;            /* TRUE to use CM's window length for all HMM filters */
  int     fwdbias_ns;           /* number of samples for do_fwdbias_sampling */
  int     do_glen;              /* TRUE to use len-dependent glc p7 thresholds */
  int     glen_min;             /* min clen for len-dependent glc p7 thr    */
  int     glen_max;             /* max clen for len-dependent glc p7 thr    */
  int     glen_step;            /* step size for halving glc p7 thr if do_glen */

  /* Parameters controlling p7 domain/envelope defintion */
  float  rt1;   	/* controls when regions are called. mocc[i] post prob >= dt1 : triggers a region around i */
  float  rt2;		/* controls extent of regions. regions extended until mocc[i]-{b,e}occ[i] < dt2            */
  float  rt3;		/* controls when regions are flagged for split: if expected # of E preceding B is >= dt3   */
  int    ns;            /* number of traceback samples for domain/envelope def */

  /* CM search options for fourth filter and final stage */
  int     fcyk_cm_search_opts;  /* CYK filter stage search opts             */
  int     final_cm_search_opts; /* final stage search opts                  */
  int     fcyk_cm_exp_mode;     /* CYK filter exp mode                      */
  int     final_cm_exp_mode;    /* final stage exp mode   e                 */
  double  fcyk_beta;            /* QDB beta for CYK filter stage            */
  double  final_beta;           /* QDB beta for final stage                 */
  double  fcyk_tau;             /* HMM bands tau for CYK filter stage       */
  double  final_tau;            /* HMM bands tau for final stage            */

  /* Accounting. (reduceable in threaded/MPI parallel version)              */
  uint64_t      nmodels;        /* # of HMMs searched                       */
  uint64_t      nseqs;	        /* # of sequences searched                  */
  uint64_t      nres;	        /* # of residues searched                   */
  uint64_t      nnodes;	        /* # of model nodes searched                */
  uint64_t      n_past_msv;	/* # windows that pass MSVFilter()          */
  uint64_t      n_past_vit;	/* # windows that pass ViterbiFilter()      */
  uint64_t      n_past_fwd;	/* # windows that pass ForwardFilter()      */
  uint64_t      n_past_gfwd;	/* # windows that pass glocal GForward()    */
  uint64_t      n_past_edef;	/* # envelopess that pass envelope definition */
  uint64_t      n_past_cyk;	/* # windows that pass CYK filter           */
  uint64_t      n_past_ins;	/* # windows that pass Inside               */
  uint64_t      n_output;	/* # alignments that make it to the final output */
  uint64_t      n_past_msvbias;	/* # windows that pass MSV bias filter      */
  uint64_t      n_past_vitbias;	/* # windows that pass Vit bias filter      */
  uint64_t      n_past_fwdbias;	/* # windows that pass Fwd bias filter      */
  uint64_t      n_past_gfwdbias;/* # windows that pass gFwd bias filter     */
  uint64_t      n_past_edefbias;/* # envelopes that pass env bias filter    */
  uint64_t      pos_past_msv;	/* # positions that pass MSVFilter()        */
  uint64_t      pos_past_vit;	/* # positions that pass ViterbiFilter()    */
  uint64_t      pos_past_fwd;	/* # positions that pass ForwardFilter()    */
  uint64_t      pos_past_gfwd;	/* # positions that pass glocal GForward()  */
  uint64_t      pos_past_edef;	/* # positions that pass env definition     */
  uint64_t      pos_past_cyk;	/* # positions that pass CYK filter         */
  uint64_t      pos_past_ins;	/* # positions that pass Inside             */
  uint64_t      pos_output;	/* # positions that make it to the final output */
  uint64_t      pos_past_msvbias;/* # positions that pass MSV bias filter */
  uint64_t      pos_past_vitbias;/* # positions that pass Vit bias filter */
  uint64_t      pos_past_fwdbias;/* # positions that pass Fwd bias filter */
  uint64_t      pos_past_gfwdbias;/*# positions that pass gFwd bias filter*/
  uint64_t      pos_past_edefbias;/* # positions that pass dom def bias filter */

  int           do_time_F1;      /* TRUE to abort after Stage 1 MSV */
  int           do_time_F2;      /* TRUE to abort after Stage 2 Vit */
  int           do_time_F3;      /* TRUE to abort after Stage 3 Fwd */
  int           do_time_F4;      /* TRUE to abort after Stage 4 glocal Fwd */
  int           do_time_F5;      /* TRUE to abort after Stage 5 env def */
  int           do_time_F6;      /* TRUE to abort after Stage 6 CYK */

  enum cm_pipemodes_e mode;    	/* CM_SCAN_MODELS | CM_SEARCH_SEQS          */
  int           do_top;         /* TRUE to do top    strand (usually TRUE) */
  int           do_bot;         /* TRUE to do bottom strand (usually TRUE) */
  int 		W;              /* window length */
  int 		clen;           /* consensus length of model */
  //float         p7_evparam[CM_p7_NEVPARAM]; /* E-value params */

  int           show_accessions;/* TRUE to output accessions not names      */
  int           show_alignments;/* TRUE to output alignments (default)      */

  int           use_cyk;        /* TRUE to use CYK instead of optimal accuracy    */
  int           align_hbanded;  /* TRUE to do HMM banded alignment, when possible */
  float         hb_size_limit;  /* maximum size in Mb allowed for HB alignment    */

  CMFILE       *cmfp;		/* COPY of open CM database (if scan mode) */
  char          errbuf[eslERRBUFSIZE];
} CM_PIPELINE;

/* Structure: CM_ALIDISPLAY
 * 
 * Alignment of a sequence to a CM, formatted for printing.
 * Based on HMMER's P7_ALIDISPLAY.
 *
 * For an alignment of L residues and names C chars long, requires
 * 7L + 2C + 30 bytes; for typical case of L=100,C=10, that's
 * <0.8 Kb.
 */
typedef struct cm_alidisplay_s {
  char *rfline;                 /* reference coord info; or NULL        */
  char *nline;                  /* negative scoring noncanonicals       */
  char *csline;                 /* consensus structure info             */
  char *model;                  /* aligned query consensus sequence     */
  char *mline;                  /* "identities", conservation +'s, etc. */
  char *aseq;                   /* aligned target sequence              */
  char *ppline;			/* posterior prob annotation; or NULL   */
  int   N;			/* length of strings                    */
  
  char *cmname;	    	        /* name of HMM                          */
  char *cmacc;			/* accession of HMM; or [0]='\0'        */
  char *cmdesc;		        /* description of HMM; or [0]='\0'      */
  int   cfrom;		        /* min bound in ccoord, start position in CM */
  int   cto;			/* max bound in ccoord, end position in CM   */
  int   clen;			/* consensus length of model            */
  
  char *sqname;			/* name of target sequence              */
  char *sqacc;			/* accession of target seq; or [0]='\0' */
  char *sqdesc;			/* description of targ seq; or [0]='\0' */
  long  sqfrom;			/* min bound in scoord, start position in sequence (1..L) */
  long  sqto;		        /* max bound in scoord, end position in sequence   (1..L) */
  long  L;			/* length of sequence                   */

  int    used_optacc;           /* TRUE if aln alg was optacc, FALSE if CYK */
  float  aln_sc;		/* if(used_optacc) avg PP of all aligned residues, else CYK score */
  int    used_hbands;           /* TRUE if aln used HMM bands, FALSE if not */
  float  matrix_Mb;             /* size of DP matrix used in Mb, either HMM banded CYK/OA or D&C CYK */
  double elapsed_secs;          /* number of seconds required for alignment */

  int   memsize;                /* size of allocated block of char memory */
  char *mem;		        /* memory used for the char data above  */
} CM_ALIDISPLAY;

#define CM_HIT_FLAGS_DEFAULT 0
#define CM_HIT_IS_INCLUDED      (1<<0)
#define CM_HIT_IS_REPORTED      (1<<1)
#define CM_HIT_IS_NEW           (1<<2)
#define CM_HIT_IS_DROPPED       (1<<3)

/* Structure: CM_HIT
 * 
 * Info about a high-scoring database hit, kept so we can output a
 * sorted list of high hits at the end.
 *
 * sqfrom and sqto are the coordinates that will be shown
 * in the results, not coords in arrays... therefore, reverse
 * complements have sqfrom > sqto
 */
typedef struct cm_hit_s {
  char          *name;		/* name of the target               (mandatory)           */
  char          *acc;		/* accession of the target          (optional; else NULL) */
  char          *desc;		/* description of the target        (optional; else NULL) */
  double         sortkey;       /* number to sort by; big is better                       */

  int64_t        start, stop;   /* start/end points of hit */
  float          score;		/* bit score of the hit (with corrections) */
  double         pvalue;	/* P-value of the hit   (with corrections) */
  double         evalue;	/* E-value of the hit   (with corrections) */
  CM_ALIDISPLAY *ad;            /* alignment display */

  uint32_t       flags;         /* CM_HIT_IS_REPORTED | CM_HIT_IS_INCLUDED | CM_HIT_IS_NEW | CM_HIT_IS_DROPPED */
} CM_HIT;

/* Structure: CM_TOPHITS
 * merging when we prepare to output results. "hit" list is NULL and
 * unavailable until after we do a sort.  
 */
typedef struct cm_tophits_s {
  CM_HIT **hit;         /* sorted pointer array                     */
  CM_HIT  *unsrt;	/* unsorted data storage                    */
  uint64_t Nalloc;	/* current allocation size                  */
  uint64_t N;		/* number of hits in list now               */
  uint64_t nreported;	/* number of hits that are reportable       */
  uint64_t nincluded;	/* number of hits that are includable       */
  int      is_sorted;	/* TRUE when h->hit valid for all N hits    */
} CM_TOPHITS;

#endif /*STRUCTSH_INCLUDED*/


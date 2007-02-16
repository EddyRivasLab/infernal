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

#include "squid.h"
#include "ssi.h"               /* CMFILE supports SSI indexes */
#include "easel.h"
#include "esl_sqio.h"
#include "cplan9.h"

/* various default parameters for CMs and CP9 HMMs */ 
#define DEFAULT_CM_CUTOFF 0.0
#define DEFAULT_CM_CUTOFF_TYPE SCORE_CUTOFF
#define DEFAULT_CP9_CUTOFF 0.0
#define DEFAULT_CP9_CUTOFF_TYPE SCORE_CUTOFF
#define DEFAULT_BETA   0.0000001
#define DEFAULT_TAU    0.0000001

#define GC_SEGMENTS 101                   /* Possible integer GC contents */

#define BE_EFFICIENT  0		/* setting for do_full: small memory mode */
#define BE_PARANOID   1		/* setting for do_full: keep whole matrix, perhaps for debugging */

/* Constants for type of cutoff */
#define SCORE_CUTOFF 0
#define E_CUTOFF     1

/* Alphabet information is declared here, and defined in globals.c.
 */
#define MAXABET           4	/* max for Alphabet_size             */
#define MAXDEGEN         17     /* maximum for Alphabet_iupac        */
#define DIGITAL_GAP      126	/* see alphabet.c:DigitizeSequence() */
#define DIGITAL_SENTINEL 127    
#define INTSCALE    1000.0      /* scaling constant for floats to integer scores */
#define LOGSUM_TBL  20000       /* controls precision of ILogsum()            */

extern int   Alphabet_type;
extern int   Alphabet_size;
extern int   Alphabet_iupac;
extern char *Alphabet;
extern char  Degenerate[MAXDEGEN][MAXABET];
extern int   DegenCount[MAXDEGEN];

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

/* For CM Plan 9 HMMs which has scores as integers */
#define INFTY       987654321   /* infinity for purposes of integer DP cells       */

/* State types. (cm->sttype[])
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

/* Node types (8) (cm->ndtype[])
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

/* Unique state identifiers  (cm->stid[])
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

			/* Information about the null model:               */
  float *null;          /*   residue probabilities [0..3]                  */

			/* Information about the state type:               */
  int   M;		/*   number of states in the model                 */
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

			/* Scaled int parameters of the log odds model:       */
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
  /* search cutoffs */
  int       cutoff_type;/* either SC_CUTOFF or E_CUTOFF                                       */
  float     cutoff;     /* min bit score or max E val to keep in a scan (depending on cutoff_type) */
  int   cp9_cutoff_type;/* either SC_CUTOFF or E_CUTOFF                                       */
  float cp9_cutoff;     /* min bit score or max E val to keep from a CP9 scan                 */
  
  /* EVD statistics for the CM */
  double     lambda[GC_SEGMENTS];   /* EVD lambda, one for each GC segment   */
  double     K     [GC_SEGMENTS];   /* EVD K, one for each GC segment        */
  double     mu    [GC_SEGMENTS];   /* EVD mu, one for each GC segment       */
  double cp9_lambda[GC_SEGMENTS];   /* CP9's EVD lambda, one for each GC segment   */
  double cp9_K     [GC_SEGMENTS];   /* CP9's EVD K, one for each GC segment        */
  double cp9_mu    [GC_SEGMENTS];   /* CP9's EVD mu, one for each GC segment       */

  int    W;             /* max d: max size of a hit (EPN 08.18.05) */
  float  el_selfsc;     /* score of a self transition in the EL state
			 * the EL state emits only on self transition (EPN 11.15.05)*/
  int   iel_selfsc;     /* scaled int version of el_selfsc         */


} CM_t;

/* status flags, cm->flags */
#define CM_LOCAL_BEGIN         (1<<0)  /* Begin distribution is active (local ali) */
#define CM_LOCAL_END           (1<<1)  /* End distribution is active (local ali)   */
#define CM_STATS               (1<<2)  /* EVD stats, mu, lambda, K are set         */
#define CM_CP9                 (1<<3)  /* CP9 HMM is valid in cm->cp9              */
#define CM_CP9STATS            (1<<4)  /* CP9 HMM has EVD stats                    */
#define CM_IS_SUB              (1<<5)  /* the CM is a sub CM                       */
#define CM_IS_FSUB             (1<<6)  /* the CM is a fullsub CM                   */
#define CM_IS_RSEARCH          (1<<7)  /* the CM was parameterized a la RSEARCH    */
#define CM_ENFORCED            (1<<8)  /* CM is reparam'ized to enforce a subseq   */

/* model configuration options, cm->config_opts */
#define CM_CONFIG_LOCAL        (1<<0)  /* configure the model for local alignment  */
#define CM_CONFIG_ENFORCE      (1<<1)  /* enforce a subseq be incl. in each parse  */
#define CM_CONFIG_ENFORCEHMM   (1<<2)  /* build CP9 HMM to only enforce subseq     */
#define CM_CONFIG_ELSILENT     (1<<3)  /* disallow EL state emissions              */
#define CM_CONFIG_ZEROINSERTS  (1<<4)  /* make all insert emissions equiprobable   */
#define CM_CONFIG_QDB          (1<<5)  /* calculate query dependent bands          */

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

/* search options, cm->search_opts */
#define CM_SEARCH_NOQDB        (1<<0)  /* DO NOT use QDB to search (QDB is default)*/
#define CM_SEARCH_HMMONLY      (1<<1)  /* use a CP9 HMM only to search             */
#define CM_SEARCH_HMMFB        (1<<2)  /* filter w/CP9 HMM, forward/backward mode  */
#define CM_SEARCH_HMMWEINBERG  (1<<3)  /* filter w/CP9 HMM, Zasha Weinberg mode    */
#define CM_SEARCH_HMMRESCAN    (1<<4)  /* rescan HMM hits b/c Forward is inf length*/
#define CM_SEARCH_HMMSCANBANDS (1<<5)  /* filter w/CP9 HMM, and derive HMM bands   */
#define CM_SEARCH_SUMS         (1<<6)  /* if using HMM bands, use posterior sums   */
#define CM_SEARCH_INSIDE       (1<<7)  /* scan with Inside, not CYK                */
#define CM_SEARCH_TOPONLY      (1<<8)  /* don't search reverse complement          */
#define CM_SEARCH_NOALIGN      (1<<9)  /* don't align hits, just report locations  */
#define CM_SEARCH_NULL2        (1<<10) /* use post hoc second null model           */
#define CM_SEARCH_CMSTATS      (1<<11) /* calculate E-value statistics for CM      */
#define CM_SEARCH_CP9STATS     (1<<12) /* calculate E-value stats for CP9 HMM      */
#define CM_SEARCH_FFRACT       (1<<13) /* filter to filter fraction cm->ffract     */
#define CM_SEARCH_RSEARCH      (1<<14) /* use RSEARCH parameterized CM             */

/* Structure: CMFILE
 * Incept:    SRE, Tue Aug 13 10:16:39 2002 [St. Louis]
 *
 * An open CM database for reading. 
 * (When writing, we just use a normal stdio.h FILE.)
 * API is implemented in cmio.c
 */
typedef struct cmfile_s {
  FILE     *f;                  /* open file for reading */
  SSIFILE  *ssi;                /* ptr to open SSI index, or NULL if unavailable */
  int       is_binary;		/* TRUE if file is in binary format */
  int       byteswap;		/* TRUE if binary and we need to swap byte order */
  SSIOFFSET offset;		/* disk offset of the CM that was read last */
  int       mode;		/* type of SSI offset (part of SSI API) */
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
typedef struct _scan_result_node_t {
  int start;
  int stop;
  int bestr;   /* Best root state */
  float score;
  Parsetree_t *tr;
} scan_result_node_t;

typedef struct _scan_results_t {
  scan_result_node_t *data;
  int num_results;
  int num_allocated;
} scan_results_t;

typedef struct _db_seq_t {
  ESL_SQ  *sq[2];
  scan_results_t *results[2];
  int chunks_sent;
  int alignments_sent;           /* -1 is flag for none queued yet */
  float best_score;              /* Best score for scan of this sequence */
  int partition;                 /* For histogram building */
} db_seq_t;

/* structure for MPI cmalign */
typedef struct _seqs_to_aln_t {
  ESL_SQ  **sq;                  /* the sequences */
  int nseq;                      /* number of sequences */
  Parsetree_t **tr;              /* parsetrees */
  char **postcode;               /* postal codes, left NULL unless do_post */
  int index;                     /* the index of the first sq (sq[0]) in master structure */
} seqs_to_aln_t;

/* The integer log odds score deckpool for integer versions of 
 * Inside and Outside, see cm_postprob.c */
typedef struct Ideckpool_s {
  int   ***pool;
  int      n;
  int      nalloc;
  int      block;
} Ideckpool_t;


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

#endif /*STRUCTSH_INCLUDED*/



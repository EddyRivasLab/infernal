#ifndef STRUCTSH_INCLUDED
#define STRUCTSH_INCLUDED

/* structs.h
 * SRE, 28 Feb 2000 [Seattle]
 * CVS $Id$
 * 
 * Declarations of structures and global variables;
 * definitions of constants; and macros.
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include "ssi.h"               /* CMFILE supports SSI indexes */

/* Alphabet information is declared here, and defined in globals.c.
 */
#define MAXABET          4	/* max for Alphabet_size             */
#define MAXDEGEN         17     /* maximum for Alphabet_iupac        */
#define DIGITAL_GAP      126	/* see alphabet.c:DigitizeSequence() */
#define DIGITAL_SENTINEL 127    
extern int   Alphabet_type;
extern int   Alphabet_size;
extern int   Alphabet_iupac;
extern char *Alphabet;
extern char  Degenerate[MAXDEGEN][MAXABET];
extern int   DegenCount[MAXDEGEN];
extern int   CMTransitionIndex[20][20];

/* The "-infinity" value in a DP matrix must be > -FLT_MAX/3, so that 
 * we can add three of them together (alpha + tsc + esc) and not get an
 * underflow error. ANSI guarantees us FLT_MAX >= 1e37. 
 * NOT_IMPOSSIBLE() exists because we can't compare floating point by 
 * equality.
 * sreLOG2() exists because we want to work in bits, and we will need
 * to take log(0).
 */
#define IMPOSSIBLE -1e36
#define NOT_IMPOSSIBLE(x)  ((x) > -9.999e35) 
#define sreLOG2(x)  ((x) > 0 ? log(x) * 1.44269504 : IMPOSSIBLE)

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


/* Structure: CM_t
 * Incept:    SRE, 9 Mar 2000 [San Carlos CA]
 * 
 * A covariance model. M states, arranged logically as a directed graph
 * (on a binary tree backbone); arranged physically as a set of arrays 0..M-1.
 *
 * State 0 is always the root state. State M-1 is always an end state.
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
  char *stid;		/*   unique state identifier; e.g. MATP-MP         */

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

  int    flags;		/* status flags                                    */
} CM_t;

/* Status flags for a CM, cm->flags
 */
#define CM_LOCAL_BEGIN  (1<<0)	/* Begin distribution is active (local ali) */
#define CM_LOCAL_END    (1<<1)  /* End distribution is active (local ali)   */


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
 * For reasons of malloc() efficiency, the binary tree is organized
 * in a set of arrays. 
 */
typedef struct parsetree_s {
  int *emitl;		/* i position in seq or ali (1..L or alen) */
  int *emitr;		/* j position in seq or ali (1..L or alen) */
  int *state;		/* y of state (0..M-1)                     */

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

#endif /*STRUCTSH_INCLUDED*/

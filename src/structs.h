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

/* Alphabet information is declared here, and defined in globals.c.
 */
#define MAXABET          4	/* maximum for Alphabet_size (used for static arrays) */
#define MAXDEGEN         17     /* maximum for Alphabet_iupac */
#define DIGITAL_GAP      126	/* see alphabet.c:DigitizeSequence() */
#define DIGITAL_SENTINEL 127    
extern int   Alphabet_type;
extern int   Alphabet_size;
extern int   Alphabet_iupac;
extern char *Alphabet;
extern char  Degenerate[MAXDEGEN][MAXABET];
extern int   DegenCount[MAXDEGEN];
extern int   CMTransitionIndex[20][20];

/* The "-infinity" value in a DP matrix should be > -FLT_MAX/3, so that 
 * we can add three of them together (alpha + tsc + esc) and not get an
 * underflow error. ANSI guarantees us FLT_MAX >= 1e37. 
 */
#define IMPOSSIBLE -1e36
#define NOT_IMPOSSIBLE(x)  ((x) > -9.999e35) /* can't compare floating point by equality */
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
				/* General information about the model:                  */
  char *name;			/*   name of the model                                   */
  char *acc;			/*   optional accession number for model                 */
  char *desc;			/*   optional description of the model                   */
  char *annote;                 /*   consensus column annotation line                    */

				/* Information about the null model:                     */
  float *null;                  /*   residue probabilities [0..3]                        */

				/* Information about the state type:                     */
  int   M;			/*   number of states in the model                       */
  char *sttype;			/*   type of state this is; e.g. MP_st                   */
  int  *ndidx;			/*   index of node this state belongs to                 */
  char *stid;			/*   unique identifier for this state type; e.g. MATP-MP */

				/* Information about its connectivity in the CM:         */
  int  *cfirst;			/*   index to first child state we connect to            */
  int  *cnum;			/*   overloaded: for non-BIF, = number of connections;   */
				/*               for BIF,     = second child S_st        */
  int  *plast;                  /*   index to first parent state we connect to           */
  int  *pnum;                   /*   number of parent connections                        */

				/* Information mapping nodes->states                     */
  int   nodes;			/*   number of nodes in the model                        */
  int  *nodemap;                /* i.e. nodemap[5] = index of first state of node 5      */
  char *ndtype;			/* type of node, e.g. MATP_nd                            */

				/* Parameters of the probabilistic model:                */
  float **t;			/*   Transition probabilities [0..M-1][0..MAXCONNECT-1]  */
  float **e;			/*   Emission probabilities.  [0..M-1][0..15]            */
  float  *begin;		/*   Local alignment start probabilities [0..M-1]        */
  float  *end;			/*   Local alignment ending probabilities [0..M-1]       */

				/* Parameters of the log odds model:                     */
  float **tsc;			/*   Transition score vector, log odds                   */
  float **esc;			/*   Emission score vector, log odds                     */
  float *beginsc;		/*   Score for ROOT_S -> state v (local alignment)       */
  float *endsc;   		/*   Score for state_v -> EL (local alignment)           */

  int    flags;			/* status flags */
} CM_t;

/* Status flags for a CM, cm->flags
 */
#define CM_LOCAL_BEGIN  (1<<0)	/* Begin distribution is active (local alignment) */
#define CM_LOCAL_END    (1<<1)  /* End distribution is active (local alignment)   */


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
  int *emitl;			/* i position in sequence or alignment (1..L or alen) */
  int *emitr;			/* j position in sequence or alignment (1..L or alen) */
  int *state;			/* y of state (0..M-1)                  */

  int *nxtl;			/* index in trace of left child  */
  int *nxtr;			/* index in trace of right child */
  int *prv;			/* index in trace of parent      */

  int  n;			/* number of elements in use so far     */
  int  nalloc;			/* number of elements allocated for     */
  int  memblock;		/* size of malloc() chunk, # of elems   */
} Parsetree_t;


/* Structure: CMConsensus_t
 * Incept:    SRE, Thu May 23 16:55:04 2002 [St. Louis]
 * 
 * Created by display.c:CreateCMConsensus(). 
 * Preprocesses a CM into consensus information that is needed by
 * display.c:CreateFancyAli().
 */
typedef struct consensus_s {
  char *cseq;                   /* consensus sequence display string; 0..clen-1     */
  char *cstr;			/* consensus structure display string; 0..clen-1    */
  int  *lpos;			/* maps node->consensus position; 0..nodes-1        */
  int  *rpos;			/* maps node->consensus position; 0..nodes-1        */
  int   clen;			/* length of cseq, cstr                             */
} CMConsensus_t;

typedef struct fancyali_s {
  char *annote;                 /* reference annotation line (NULL if unavail) */
  char *cstr;			/* CM consensus structure line */
  char *cseq;			/* CM consensus sequence line  */
  char *mid;			/* alignment identity middle line */
  char *aseq;			/* aligned target sequence */
  int   len;
} Fancyali_t;



#endif /*STRUCTSH_INCLUDED*/

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
extern int   Alphabet_type;
extern int   Alphabet_size;
extern int   Alphabet_iupac;
extern char *Alphabet;


/* State types. (9)
 */
#define STATETYPES 9            /* different kinds of states         */
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

/* Node types. (8)
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

/* Unique state identifiers. (20)
 */
#define UNIQUESTATES 20

#define DUMMY   -1
#define ROOT_S  0
#define ROOT_IL 1
#define ROOT_IR 2
#define BEGL_S  3
#define BEGR_S  4
#define BEGR_IL 5
#define MATP_D  6
#define MATP_MP 7
#define MATP_ML 8
#define MATP_MR 9
#define MATP_IL 10
#define MATP_IR 11
#define MATL_D  12
#define MATL_ML 13
#define MATL_IL 14
#define MATR_D  15
#define MATR_MR 16
#define MATR_IR 17
#define END_E   18
#define BIF_B   19

/* Flags used in InsertTraceNode()
 */
#define TRACE_LEFT_CHILD  1
#define TRACE_RIGHT_CHILD 2


/* Structure: CM_t
 * Incept:    SRE, 9 Mar 2000 [San Carlos CA]
 * 
 * A covariance model. M states, arranged logically as a directed acyclic graph
 * (on a binary tree backbone); arranged physically as a set of arrays 0..M-1.
 *
 * State 0 is always the root state. State M-1 is always an end state.
 */
typedef struct cm_s {			
				/* General information about the model:                  */
  char *name;			/*   name of the model                                   */
  char *acc;			/*   optional accession number for model                 */
  char *desc;			/*   optional description of the model                   */
  int   nodes;			/*   number of nodes in the model                        */
  int   M;			/*   number of states in the model                       */

				/* Information about the state type:                     */
  char *sttype;			/*   type of state this is; e.g. MP_st                   */
  int  *ndidx;			/*   index of node this state belongs to                 */
  char *ndtype;			/*   type of node this state belongs to; e.g. MATP_nd    */
  char *stid;			/*   unique identifier for this state type; e.g. MATP-MP */

				/* Information about its connectivity in the CM:         */
  int  *cfirst;			/*   index to first child state we connect to            */
  int  *cnum;			/*   overloaded: for non-BIF, = number of connections;   */
				/*               for BIF,     = second child S_st        */
  int  *plast;                  /*   index to first parent state we connect to           */
  int  *pnum;                   /*   number of parent connections                        */

				/* Parameters of the probabilistic model:                */
  float **t;			/*   Transition probabilities [0..M-1][0..ynum-1]        */
  float **e;			/*   Emission probabilities.  [0..M-1][0..15]            */

				/* Parameters of the log odds model:                     */
  float **tsc;			/*   Transition score vector, log odds                   */
  float **esc;			/*   Emission score vector, log odds                     */
} CM_t;


/* Structure: Parsetree_t
 * Incept:    SRE 29 Feb 2000 [Seattle]
 * 
 * Binary tree structure for storing a traceback of an alignment.
 * 
 * Also used for tracebacks of model constructions. Then, some
 * fields are misused: "state" is used for a node (not state) index,
 * and "type" is used for a node (not state) type.
 *
 * For reasons of malloc() efficiency, the binary tree is organized
 * in a set of arrays. 
 */
typedef struct parsetree_s {
  int *emitl;			/* i position (0..N-1) or -1 if nothing */
  int *emitr;			/* j position (0..N-1) or -1 if nothing */
  int *state;			/* y of state (0..M-1)                  */
  int *type;			/* type of state                        */

  int *nxtl;			/* index in trace of left child  */
  int *nxtr;			/* index in trace of right child */
  int *prv;			/* index in trace of parent      */

  int  n;			/* number of elements in use so far     */
  int  nalloc;			/* number of elements allocated for     */
  int  memblock;		/* size of malloc() chunk, # of elems   */
} Parsetree_t;




#endif /*STRUCTSH_INCLUDED*/

/* The all-encompassing include file for INFERNAL.
 *
 *    1. CM_CORE:      a core model.
 *    2. CM_PROFILE: a scoring profile, and its implicit model.
 *    3. CM_BG:      a null (background) model.
 *    4. CM_TRACE:   a traceback path (alignment of seq to profile).
 *    5. CM_CMFILE:  a CM save file or database, open for reading.
 *    6. CM_GMX:     a "generic" dynamic programming matrix
 *    7. CM_DPRIOR:  mixture Dirichlet prior for CMs
 *    8. Other routines in Infernal's exposed API.
 * 
 * EPN, Wed Aug  1 07:09:21 2007 [DC]
 * SVN $Id$
 */
#ifndef CM_INFERNALH_INCLUDED
#define CM_INFERNALH_INCLUDED

#include "cm_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"	/* ESL_DSQ, ESL_ALPHABET */
#include "esl_dmatrix.h"	/* ESL_DMATRIX           */
#include "esl_msa.h"		/* ESL_MSA               */
#include "esl_random.h"		/* ESL_RANDOMNESS        */
#include "esl_sqio.h"		/* ESL_SQ                */
#include "esl_histogram.h"      /* ESL_HISTOGRAM         */
#include "esl_dirichlet.h"	/* ESL_MIXDCHLET         */

/* Search modes.
 */
#define cm_NO_MODE   0
#define cm_LOCAL     1		/*  local */
#define cm_GLOCAL    2		/* glocal */

#define cm_IsLocal(mode)  (mode == cm_LOCAL)


/*****************************************************************
 * 1. CM_CORE: a core model.
 *****************************************************************/

/* Flag codes for cm->flags.
 * Flags marked with ! may not be changed nor used for other meanings;
 * such flags were stored in old CM files, and we must preserve their
 * meaning to preserve reverse compatibility.
 */
#define cmH_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define cmH_DESC    (1<<1)    /* description exists                              !*/
#define cmH_RF      (1<<2)    /* #RF annotation available                        !*/
#define cmH_HASPROB (1<<3)    /* obsolete (was: model in probability form)       !*/
#define cmH_STATS   (1<<4)    /* obsolete (was: model has EVD stats calibrated)  !*/
#define cmH_ACC     (1<<5)    /* accession number is available                   !*/
#define cmH_GA      (1<<6)    /* gathering threshold available                   !*/
#define cmH_TC      (1<<7)    /* trusted cutoff available                        !*/
#define cmH_NC      (1<<8)    /* noise cutoff available                          !*/

/* State types. (cm->sttype[])
 */
#define cmH_MAXTRANSITIONS 6            /* maximum number transitions from a state */

/*
enum p7t_statetype_e {
  D_st    =  0,
  MP_st   =  1,
  ML_st   =  2,
  MR_st   =  3,
  IL_st   =  4,
  IR_st   =  5,
  S_st    =  6,
  E_st    =  7,
  B_st    =  8,
  EL_st   =  9,
};
#define cmT_NSTATETYPES 10
*/

#define  D_st   0
#define  MP_st  1
#define  ML_st  2
#define  MR_st  3
#define  IL_st  4
#define  IR_st  5
#define  S_st   6
#define  E_st   7
#define  B_st   8
#define  EL_st  9

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

/* Structure: CM_CORE
 * Incept:    SRE, 9 Mar 2000 [San Carlos CA]
 * 
 * A covariance model. M states, arranged logically as a directed graph
 * (on a binary tree backbone); arranged physically as a set of arrays 0..M-1.
 *
 * State 0 is always the root state. State M-1 is always an end state.
 * 
 */
typedef struct cm_core_s {			
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
  float   el;           /*   EL->EL self transition probability            */

  /* Annotation. Everything but <name> and <cs> is optional. Flags are set when
   * optional values are set. All the char *'s are proper nul-terminated
   * strings, not just arrays. (hmm->map is an int array).
   */
  char  *name;                  /* name of the model                     (mandatory) */ /* String, \0-terminated */
  char  *cs;                    /* consensus structure line      1..M    (mandatory) */ /* String; 0=' ', M+1='\0' */
  char  *acc;			/* accession number of model (Rfam)      (cmH_ACC)   */ /* String, \0-terminated */
  char  *desc;                  /* brief (1-line) description of model   (cmH_DESC)  */ /* String, \0-terminated */
  char  *rf;                    /* reference line from alignment 1..M    (cmH_RF)    */ /* String; 0=' ', M+1='\0' */
  char  *comlog;		/* command line(s) that built model      (mandatory) */ /* String, \0-terminated */
  int    nseq;			/* number of training sequences          (mandatory) */
  float  eff_nseq;		/* effective number of seqs (<= nseq)    (mandatory) */
  char  *ctime;			/* creation date                         (mandatory) */
  int    checksum;              /* checksum of training sequences        (mandatory) */
  /* emit map goes here */

  float  ga;	                /* gathering threshold (bits)  (cmH_GA) */
  float  tc;                    /* trusted cutoff (bits)       (cmH_TC) */
  float  nc;	                /* noise cutoff (bits)         (cmH_NC) */
  off_t  offset;                /* CM record offset on disk */
  int    flags;               /* status flags */

  const ESL_ALPHABET *abc; /* ptr to alphabet info (cm->abc->K is alphabet size) */
  const CM_BG *bg;         /* ptr to the background (null) model                 */
} CM_CORE;



/*****************************************************************
 * 2. CM_PROFILE: a scoring profile, and its implicit model.
 *****************************************************************/

/* Indices for residue emission score vectors
 */
enum cmp_rsc_e {
  cmP_MSC = 0, 
  cmP_ISC = 1
};
#define cmP_NR 2

typedef struct cm_profile_s {
  int     mode;        	/* configured algorithm mode (e.g. cm_LOCAL)              */ 
  int     M;		/* number of states in the model                          */
  float  *tsc;          /* transitions  [0.1..M-1][0..cmP_MAXCONNECT-1], hand-indexed  */
  float **rsc;          /* emissions    [0.1..M-1][cm_MAXABET*cm_MAXABET], hand-indexed       */
  /* Below should be temporary, HMMER3 doesn't have them, all begin/ends equiprobable  */
  float  *beginsc;	/*   Local alignment start probabilities [0..M-1]  */
  float  *endsc;		/*   Local alignment ending probabilities [0..M-1] */

  /* Objects we keep references to */
  const ESL_ALPHABET    *abc_r;	/* copy of pointer to appropriate alphabet     */
  const struct cm_core_s  *cm_r;	/* who's your daddy                            */
  const struct cm_bg_s  *bg_r;	/* background null model                       */
} CM_PROFILE;




/*****************************************************************
 * 3. CM_BG: a null (background) model.
 *****************************************************************/

typedef struct cm_bg_s {
  ESL_ALPHABET *abc;		/* reference to alphabet in use       */
  float *f;			/* residue frequencies [0..K-1] */
} CM_BG;

/*****************************************************************
 * 4. CM_PARSETREE:  a parsetree (alignment of seq to profile).
 *****************************************************************/
/* Structure: CM_PARSETREE
 *
 * Incept:    SRE 29 Feb 2000 [Seattle]
 * 
 * Binary tree structure for storing a traceback of an alignment.
 * 
 * Also used for tracebacks of model constructions. Then, 
 * "state" is misused for a node (not state) index. 
 * 
 * Example of a traceback (from cm_parsetree_Dump(), from a tRNA
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
typedef struct cm_parsetree_s {
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
} CM_PARSETREE;


/*****************************************************************
 * 5. CM_CMFILE:  an CM save file or database, open for reading.
 *****************************************************************/

typedef struct cm_cmfile_s {
  FILE     *f;                  /* open file for reading */
  char     *fname;	        /* name of the CM file; [STDIN] if -           */
  int       is_binary;		/* TRUE if file is in binary format */
  int       byteswap;		/* TRUE if binary and we need to swap byte order */
  int (*parser)(struct cm_cmfile_s *, ESL_ALPHABET **, CM_CORE **);  /* parsing function */
  int       mode;		/* type of SSI offset (part of SSI API) */
  int           do_gzip;	/* TRUE if f is "gzip -dc |" (will pclose(f))    */ 
  int           do_stdin;       /* TRUE if f is stdin (won't close f)            */
  ESL_SSI  *ssi;                /* open SSI index; or NULL if none. */
} CM_CMFILE;




/*****************************************************************
 * 6. CM_GMX: a "generic" dynamic programming matrix
 * EPN, Wed Aug  1 07:52:44 2007 UNTOUCHED, not needed for cmbuild
 *****************************************************************/

enum cmg_scells_e {
  cmG_M = 0,
  cmG_I = 1,
  cmG_D = 2,
};
#define cmG_NSCELLS 3

enum cmg_xcells_e {
  cmG_E  = 0,
  cmG_N  = 1,
  cmG_J  = 2,
  cmG_B  = 3,
  cmG_C  = 4
};
#define cmG_NXCELLS 5


typedef struct cm_gmx_s {
  int  M;		/* actual model dimension (model 1..M)    */
  int  L;		/* actual sequence dimension (seq 1..L)   */
  
  size_t ncells;	/* current cell allocation limit: >= (M+1)*(L+1) */
  size_t nrows;    	/* current row allocation limit:  >= L+1  */

  float **dp;           /*  [0.1..L][0.1..M][0..cmG_NSCELLS-1] */
  float  *xmx;          /*  [0.1..L][0..cmG_NXCELLS-1]         */

  float  *xmx_mem;	
  float  *dp_mem;
} CM_GMX;



/*****************************************************************
 * 7. CM_DPRIOR: mixture Dirichlet prior for profile CMs
 *****************************************************************/

/* Structure: Prior_t
 * 
 * Dirichlet priors on all model parameters. 
 */
typedef struct cm_dprior_s {
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
} CM_DPRIOR;


/*****************************************************************
 * 8. Other routines in INFERNAL's exposed API.
 *****************************************************************/

/* build.c */
extern int cm_Handmodelmaker(ESL_MSA *msa,                CM_CORE **ret_cm, CM_TRACE ***ret_tr);
extern int cm_Fastmodelmaker(ESL_MSA *msa, float symfrac, CM_CORE **ret_cm, CM_TRACE ***ret_tr);

/* errors.c */
extern void cm_Die (char *format, ...);
extern void cm_Fail(char *format, ...);

/* eweight.c */
extern int  cm_EntropyWeight(const CM_CORE *cm, const CM_BG *bg, const CM_DPRIOR *pri, double infotarget, double *ret_Neff);

/* infernal.c */
extern void  cm_banner(FILE *fp, char *progname, char *banner);
extern float cm_SILO2Lod(int silo);
extern int   cm_AminoFrequencies(float *f);

/* mpisupport.c */
#ifdef HAVE_MPI
extern int cm_core_MPISend(CM_CORE *cm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_core_MPIPackSize(CM_CORE *cm, MPI_Comm comm, int *ret_n);
extern int cm_core_MPIPack(CM_CORE *cm, char *buf, int n, int *position, MPI_Comm comm);
extern int cm_core_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, CM_CORE **ret_cm);
extern int cm_core_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, CM_CORE **ret_cm);

extern int cm_profile_MPISend(CM_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_profile_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, const CM_BG *bg,
			      char **buf, int *nalloc,  CM_PROFILE **ret_gm);
#endif /*HAVE_MPI*/

/* cm_bg.c */
extern CM_BG *cm_bg_Create(const ESL_ALPHABET *abc);
extern CM_BG *cm_bg_CreateUniform(const ESL_ALPHABET *abc);
extern int    cm_bg_Dump(FILE *ofp, CM_BG *bg);
extern void   cm_bg_Destroy(CM_BG *bg);
extern int    cm_bg_SetLength(CM_BG *bg, int L);
extern int    cm_bg_NullOne(const CM_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

/* cm_core.c */
/*      1. The CM_CORE object: allocation, initialization, destruction. */
extern CM_CORE *cm_core_Create(int M, const ESL_ALPHABET *abc);
extern CM_CORE *cm_core_CreateShell(void);
extern int     cm_core_CreateBody(CM_CORE *cm, int M, const ESL_ALPHABET *abc);
extern void    cm_core_Destroy(CM_CORE *cm);
extern int     cm_core_CopyParameters(const CM_CORE *src, CM_CORE *dest);
extern CM_CORE *cm_core_Duplicate(const CM_CORE *cm);
extern int     cm_core_Scale(CM_CORE *cm, double scale);
extern int     cm_core_Zero(CM_CORE *cm);
extern char   *cm_core_DescribeStatetype(char st);
/*      2. Convenience routines for setting fields in a CM. */
extern int     cm_core_SetName(CM_CORE *cm, char *name);
extern int     cm_core_SetAccession(CM_CORE *cm, char *acc);
extern int     cm_core_SetDescription(CM_CORE *cm, char *desc);
extern int     cm_core_AppendComlog(CM_CORE *cm, int argc, char **argv);
extern int     cm_core_SetCtime(CM_CORE *cm);
/*      3. Convenience routines for getting info from a CM. */
extern int   cm_core_CountStatetype(CM_t *cm, char type);
extern int   cm_core_SegmentCountStatetype(CM_t *cm, int r, int z, char type);
extern int   cm_core_SubtreeCountStatetype(CM_t *cm, int v, char type);
extern int   cm_core_SubtreeFindEnd(CM_t *cm, int v);
extern int   cm_core_CalculateStateIndex(CM_t *cm, int node, char utype);
extern int   cm_core_TotalStatesInNode(int ndtype);
extern int   cm_core_SplitStatesInNode(int ndtype);
extern int   cm_core_InsertStatesInNode(int ndtype);
extern int   cm_core_StateDelta(int sttype);
extern int   cm_core_StateLeftDelta(int sttype);
extern int   cm_core_StateRightDelta(int sttype);
extern int   cm_core_NEmitAlph(const ESL_ALPHABET *abc, int sttype);

/*      3. Renormalization and rescaling counts in core CMs. */
extern int     cm_core_Rescale(CM_CORE *cm, float scale);
extern int     cm_core_Renormalize(CM_CORE *cm);
/*      4. Debugging and development code. */
extern int     cm_core_Dump(FILE *fp, CM_CORE *cm);
extern int     cm_core_Sample        (ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, CM_CORE **ret_cm);
extern int     cm_core_SampleUngapped(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, CM_CORE **ret_cm);
extern int     cm_core_SampleEnumerable(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, CM_CORE **ret_cm);
extern int     cm_core_SampleUniform (ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, 
				     float tmi, float tii, float tmd, float tdd,  CM_CORE **ret_cm);
extern int     cm_core_Compare(CM_CORE *h1, CM_CORE *h2, float tol);
extern int     cm_core_Validate(CM_CORE *cm, float tol, char *errbuf);
/*      5. Other routines in the API */
extern int     cm_core_CalculateOccupancy(const CM_CORE *cm, float *occ);



/* cm_cmfile.c */
extern int  cm_cmfile_Open(char *filename, char *env, CM_CMFILE **ret_hfp);
extern void cm_cmfile_Close(CM_CMFILE *hfp);
extern int  cm_cmfile_Write(FILE *fp, CM_CORE *cm);
extern int  cm_cmfile_Read(CM_CMFILE *hfp, ESL_ALPHABET **ret_abc,  CM_CORE **ret_cm);
extern int  cm_cmfile_PositionByKey(CM_CMFILE *hfp, const char *key);

/* cm_prior.c */
extern CM_DPRIOR *cm_dprior_CreateAmino(void);
extern CM_DPRIOR *cm_dprior_CreateNucleic(void);
extern CM_DPRIOR *cm_dprior_CreateLaplace(ESL_ALPHABET *abc);
extern void       cm_dprior_Destroy(CM_DPRIOR *pri);
extern int        cm_ParameterEstimation(CM_CORE *cm, const CM_DPRIOR *pri);

/* cm_profile.c */
extern CM_PROFILE *cm_profile_Create(int M, const ESL_ALPHABET *abc);
extern CM_PROFILE *cm_profile_Clone(const CM_PROFILE *gm);
extern int         cm_profile_SetNullEmissions(CM_PROFILE *gm);
extern void        cm_profile_Destroy(CM_PROFILE *gm);
extern int         cm_profile_IsLocal(const CM_PROFILE *gm);
extern int         cm_profile_IsMultihit(const CM_PROFILE *gm);
extern int         cm_profile_GetT(const CM_PROFILE *gm, char st1, int k1, 
				   char st2, int k2, float *ret_tsc);
extern int         cm_profile_Validate(const CM_PROFILE *gm, float tol);
extern int         cm_profile_Compare(CM_PROFILE *gm1, CM_PROFILE *gm2, float tol);

#endif /*CM_INFERNALH_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/

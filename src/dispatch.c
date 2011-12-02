/************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 ************************************************************/
/* dispatch.c
 * 
 * The two all-important dispatch functions:
 * DispatchSearch()     calls appropriate DP search functions.
 * DispatchAlignments() calls appropriate DP alignment functions.
 * 
 * EPN, Wed Dec  6 06:11:46 2006
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "easel.h"
#include "esl_msa.h"         
#include "esl_sse.h"         
#include "esl_stopwatch.h"   
#include "esl_vectorops.h"

#include "hmmer.h"
#if 0
#include "impl_sse.h"
#endif

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* Function: DispatchAlignments()
 * Incept:   EPN, Thu Nov 15 11:35:23 2007
 *
 * Purpose:  Given a CM and sequences, do preliminaries, call the
 *           correct alignment function and return parsetrees and
 *           optionally posterior code strings (if cm->align_opts &
 *           CM_ALIGN_POST).  We align the seqs_to_aln->nseq ESL_SQ sq
 *           sequences and store parsetrees or CP9 traces and/or
 *           posterior codes in seqs_to_aln.
 *
 * Args:     CM             - the covariance model
 *           errbuf         - char buffer for reporting errors
 *           seqs_to_aln    - the sequences
 *           bdump_level    - verbosity level for band related print statements
 *           debug_level    - verbosity level for debugging print statements
 *           silent_mode    - TRUE to not print anything, FALSE to print scores 
 *           do_null3       - TRUE to apply null3 correction to scores, FALSE not to
 *           do_trunc       - TRUE if we're aligning with truncated CYK/Inside/Outside
 *           r              - source of randomness (NULL unless CM_ALIGN_SAMPLE)
 *           size_limit     - max number of Mb for a DP matrix, if requestd matrix is bigger return eslERANGE 
 *           ofp            - output file to print scores to as we're aligning
 *           sfp            - special output file to print extra score information to, NULL if none 
 *           iidx           - initial index to use for descriptive output lines (often '1')
 * 
 * Returns:  eslOK on success;
 *           eslERANGE if required memory for a DP matrix is too big;
 *           eslEINCOMPAT if input parameters violate contract;
 *           eslEMEM on memory allocation error;
 *           eslFAIL if some other error;
 *           if(!eslOK) errbuf is filled with informative error message
 */
int
DispatchAlignments(CM_t *cm, char *errbuf, seqs_to_aln_t *seqs_to_aln, 
		   int bdump_level, int debug_level, int silent_mode, int do_null3, TruncOpts_t *tro,
		   ESL_RANDOMNESS *r, 
		   float size_limit, FILE *ofp, FILE *sfp, int iidx,
		   int pad7, int len7, float sc7, int end7, 
		   float mprob7, float mcprob7, float iprob7, float ilprob7)
{
  int status;
  ESL_STOPWATCH *watch;         /* for timings */
  ESL_STOPWATCH *watch2;        /* for timings */
  int nalign   = 0;             /* number of sequences we're aligning */
  ESL_DSQ *cur_dsq;             /* ptr to digitized sequence we're currently aligning */
  Parsetree_t **cur_tr;         /* pointer to the pointer to the parsetree we're currently creating */
  int L;                        /* length of sequence/subseq we're currently aligning */
  int i;                        /* counter over sequences */
  int v;                        /* state counter */
  char        **postcode = NULL;/* posterior decode array of strings */
  Parsetree_t **tr       = NULL;/* parse trees for the sequences */
  CP9trace_t  **cp9_tr   = NULL;/* CP9 traces for the sequences */
  float         sc;		/* score for one sequence alignment */
  float         ins_sc;		/* inside score for one sequence */
  float         maxsc;	        /* max score in all seqs */
  float         minsc;	        /* min score in all seqs */
  float         avgsc;      	/* avg score over all seqs */
  float         tmpsc;          /* temporary score */
  float         struct_sc;      /* structure component of the score */
  float         primary_sc;     /* primary sequence component of the score */
  float         null3_correction; /* correction in bits due to NULL3 */
  int           namewidth;      /* for dynamic width of name strings for nice tab delimited formatting */
  CMEmitMap_t  *cur_emap;       /* for passing to ParsetreeScore() to get spos epos */
  char          time_buf[128];  /* string for printing timings (safely holds up to 10^14 years) */

  /* variables related to CM Plan 9 HMMs */
  CP9_t       *hmm;             /* constructed CP9 HMM */
  CP9Bands_t  *cp9b;            /* data structure for hmm bands (bands on the hmm states) 
				 * and arrays for CM state bands, derived from HMM bands */
  CP9Map_t       *cp9map;       /* maps the hmm to the cm and vice versa */
  float           swentry;	/* S/W aggregate entry probability       */
  float           swexit;       /* S/W aggregate exit probability        */

  /* variables related to the do_sub option */
  int                spos;         /* HMM node most likely to have emitted posn 1 of target seq */
  int                spos_state;   /* HMM state type for curr spos 0=match or 1=insert */
  int                epos;         /* HMM node most likely to have emitted posn L of target seq */
  int                epos_state;   /* HMM state type for curr epos 0=match or  1=insert */
  int                cur_spos, cur_epos;
  
  CMSubMap_t        *submap;
  CM_t              *sub_cm;       /* sub covariance model                      */
  CP9_t             *sub_hmm;      /* constructed CP9 HMM */
  CP9Map_t          *sub_cp9map;   /* maps the sub_hmm to the sub_cm and vice versa */
  CP9Bands_t        *sub_cp9b;     /* data structure for hmm bands (bands on the hmm states) 
				    * and arrays for CM state bands, derived from HMM bands */

  CM_t              *orig_cm;      /* the original, template covariance model the sub CM was built from */
  CP9_t             *orig_hmm;     /* original CP9 HMM built from orig_cm */
  CP9Map_t          *orig_cp9map;  /* original CP9 map */
  CP9Bands_t        *orig_cp9b;    /* original CP9Bands */
  Parsetree_t       *orig_tr;      /* parsetree for the orig_cm; created from the sub_cm parsetree */

  /* variables related to query dependent banding (qdb) */
  int    expand_flag;           /* TRUE if the dmin and dmax vectors have just been 
				 * expanded (in which case we want to recalculate them 
				 * before we align a new sequence), and FALSE if not*/
  int *orig_dmin;               /* original dmin values passed in */
  int *orig_dmax;               /* original dmax values passed in */

  /* variables related to non-banded cyk/inside/outside */
  CM_MX             *mx      = NULL;  /* alpha DP matrix for non-banded CYK/Inside() */
  CM_MX             *out_mx  = NULL;  /* outside matrix for HMM banded Outside() */
  CM_SHADOW_MX      *shmx    = NULL;  /* shadow matrix for non-banded tracebacks */
  CM_EMIT_MX        *emit_mx = NULL;  /* emit matrix for OptAccAlign() */

  float             *parsesc; /* parsetree scores of each sequence */
  float             *parsepp; /* optimal parse posterior probability of each sequence, if any */
  float             *parse_struct_sc; /* contribution of MATP emissions - marginalized emissions to parse score, approximation of 'structural contribution' to score */
  int have_parsetrees = FALSE; /* TRUE if we'll be creating parsetrees for each seq, TRUE if (!do_hmmonly && !do_scoreonly && !do_inside) */

  /* declare and initialize options */
  int do_small     = FALSE;   /* TRUE to use D&C small alignment algs */
  int do_local     = FALSE;   /* TRUE to do local alignment */
  int do_qdb       = FALSE;   /* TRUE to do QDB alignment */
  int do_hbanded   = FALSE;   /* TRUE to do HMM banded alignment */
  int use_sums     = FALSE;   /* TRUE to use posterior sums for HMM banded alignment */
  int do_sub       = FALSE;   /* TRUE to align to a sub CM */
  int do_hmmonly   = FALSE;   /* TRUE to align with an HMM only */
  int do_scoreonly = FALSE;   /* TRUE to only calculate the score */
  int do_inside    = FALSE;   /* TRUE to do Inside also */
  int do_post      = FALSE;   /* TRUE to calculate posterior probabilities */
  int do_check     = FALSE;   /* TRUE to check posteriors from Inside/Outside */
  int do_sample    = FALSE;   /* TRUE to sample from an Inside matrix */
  int do_optacc    = FALSE;   /* TRUE to find optimally accurate alignment instead of CYK */
  int do_hmmsafe   = FALSE;   /* TRUE to realign seqs with HMM banded parses < 0. bits (only works if !do_optacc && !do_post && do_hbanded)*/
  int do_p7banded  = FALSE;   /* TRUE to use a HMMER3 plan 7 HMM to get pins prior to CP9 banding, requires do_hbanded == TRUE */

  /* TEMP */
  P7_GMX          *gx      = NULL;     /* DP matrix                               */
#if 0
  P7_OMX          *ox      = NULL;     /* optimized DP matrix                     */
#endif
  P7_PROFILE *gm = NULL;
  P7_BG *bg = NULL;
  P7_TRACE *p7_tr;
  /*P7_ALIDISPLAY  *ad      = NULL;*/
  double **phi;       /* phi array, phi[k][v] is expected number of times
			 state v (0 = match, 1 insert, 2 = delete) in 
			 cp9 hmm node k is visited. Node 0 is special, 
			 state 0 = B state, state 1 = N_state, state 2 = NULL */
  int *kmin, *kmax;
  CMEmitMap_t  *emap;           
  int k;
#if 0
  int ipos
  int cm_k;
  double i_ncells_total = 0.;
  double ncells_banded = 0.;
  double ncells_total = 0.;
  int i_npins, i_ncorrect;
  int i_nodes_n, i_nodes_ncorrect;
  int *cm_i2k;
#endif
  int i_ncells_banded = 0;
  int nodes_n, nodes_ncorrect;
  int npins, ncorrect;
  npins = ncorrect = 0;
  nodes_n = nodes_ncorrect = 0;
#if 0 
  ox  = p7_omx_Create(200, 0, 0);       /* ox is a one-row matrix for M=200 */
#endif
  bg = p7_bg_Create(cm->abc);
  gx = p7_gmx_Create(200, 400);	/* initial alloc is for M=200, L=400; will grow as needed */
  gm = p7_profile_Create (cm->mlp7->M, cm->mlp7->abc);
  p7_ProfileConfig(cm->mlp7, bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; MSVFilter requires local mode */
#if 0 
  om = p7_oprofile_Create(cm->mlp7_om, cm->mlp7->M, cm->mlp7->abc);
  p7_oprofile_Convert(gm, cm->mlp7_om); 
#endif
  if ((p7_tr = p7_trace_Create()) == NULL)  cm_Fail("trace creation failed");

  /*p7_omx_SetDumpMode(stdout, ox, TRUE);*/
  int *p7_i2k;
  /* TEMP */

  /* Contract check */
  if(!(cm->flags & CMH_BITS))                            ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CMH_BITS flag down.\n");
  if(r == NULL && (cm->align_opts & CM_ALIGN_SAMPLE))    ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), no source of randomness, but CM_ALIGN_SAMPLE alignment option on.\n");
  if((cm->align_opts & CM_ALIGN_POST)      && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CM_ALIGN_POST and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_HMMVITERBI)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CM_ALIGN_SCOREONLY and CM_ALIGN_HMMVITERBI options are incompatible.\n");
  if((cm->align_opts & CM_ALIGN_SCOREONLY) && (cm->align_opts & CM_ALIGN_POST))       ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), CM_ALIGN_SCOREONLY and CM_ALIGN_POST options are incompatible.\n");
  if(!silent_mode && ofp == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), not silent mode, but ofp is NULL\n");

  if(seqs_to_aln->sq        == NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), seqs_to_aln->sq is NULL.\n");
  if(seqs_to_aln->tr        != NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), seqs_to_aln->tr is non-NULL.\n");
  if(seqs_to_aln->cp9_tr    != NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), seqs_to_aln->cp9_tr is non-NULL.\n");
  if(seqs_to_aln->postcode  != NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), seqs_to_aln->postcode is non-NULL.\n");
  if(seqs_to_aln->sc        != NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), seqs_to_aln->sc is non-NULL.\n");
  
  /* save a copy of the align_opts we entered function with, we may change some of these for
   * individual target sequences, and we want to be able to change them back
   */
  if(cm->align_opts  & CM_ALIGN_SMALL)      do_small     = TRUE;
  if(cm->config_opts & CM_CONFIG_LOCAL)     do_local     = TRUE;
  if(cm->align_opts  & CM_ALIGN_QDB)        do_qdb       = TRUE;
  if(cm->align_opts  & CM_ALIGN_HBANDED)    do_hbanded   = TRUE;
  if(cm->align_opts  & CM_ALIGN_SUMS)       use_sums     = TRUE;
  if(cm->align_opts  & CM_ALIGN_SUB)        do_sub       = TRUE;
  if(cm->align_opts  & CM_ALIGN_HMMVITERBI) do_hmmonly   = TRUE;
  if(cm->align_opts  & CM_ALIGN_INSIDE)     do_inside    = TRUE;
  if(cm->align_opts  & CM_ALIGN_POST)       do_post      = TRUE;
  if(cm->align_opts  & CM_ALIGN_CHECKINOUT) do_check     = TRUE;
  if(cm->align_opts  & CM_ALIGN_SCOREONLY)  do_scoreonly = TRUE;
  if(cm->align_opts  & CM_ALIGN_SAMPLE)     do_sample    = TRUE;
  if(cm->align_opts  & CM_ALIGN_OPTACC)     do_optacc    = TRUE;
  if(cm->align_opts  & CM_ALIGN_HMMSAFE)    do_hmmsafe   = TRUE;
  if(cm->align_opts  & CM_ALIGN_P7BANDED)   do_p7banded  = TRUE;

  /* another contract check */

  if((do_inside + do_post + do_hmmonly + do_scoreonly) > 1) { 
    printf("\tdo_inside = %d\n\tdo_post = %d\n\tdo_hmmonly = %d\n\tdo_scoreonly = %d\n", do_inside, do_post, do_hmmonly, do_scoreonly);
    ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), exactly 0 or 1 of the above must be TRUE (== 1).");
  }
  if((do_qdb + do_optacc) > 1) { 
    printf("\tdo_qdb = %d\n\tdo_optacc = %d\n", do_qdb, do_optacc);
    ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), exactly 0 or 1 of the above must be TRUE (== 1).");
  }
  if((do_qdb + do_post) > 1) { 
    printf("\tdo_qdb = %d\n\tdo_post = %d\n", do_qdb, do_post);
    ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), exactly 0 or 1 of the above must be TRUE (== 1).");
  }
  if(do_p7banded && (!do_hbanded)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), do_p7banded is TRUE but do_hbanded is FALSE.");

  if(debug_level > 0) {
    printf("do_local    : %d\n", do_local);
    printf("do_qdb      : %d\n", do_qdb);
    printf("do_hbanded  : %d\n", do_hbanded);
    printf("do_p7banded : %d\n", do_p7banded);
    printf("use_sums    : %d\n", use_sums);
    printf("do_sub      : %d\n", do_sub);
    printf("do_hmmonly  : %d\n", do_hmmonly);
    printf("do_inside   : %d\n", do_inside);
    printf("do_small    : %d\n", do_small);
    printf("do_post     : %d\n", do_post);
    printf("do_check    : %d\n", do_check);
    printf("do_scoreonly: %d\n", do_scoreonly);
    printf("do_sample   : %d\n", do_sample);
    printf("do_optacc   : %d\n", do_optacc);
    printf("do_hmmsafe  : %d\n", do_hmmsafe);
  }
#if PRINTNOW
  printf("pad7:      %4d\n", pad7);
  printf("len7:      %4d\n", len7);
  printf("sc7:       %.2f\n", sc7);
  printf("end7:      %4d\n", end7);
  printf("mprob7:  %.4f\n",  mprob7);
  printf("mcprob7: %.4f\n",  mcprob7);
  printf("iprob7:  %.4f\n",  iprob7);
  printf("ilprob7: %.4f\n",  ilprob7);
  printf("\n");
#endif

  nalign = seqs_to_aln->nseq;

  /* If sqmode: potentially allocate tr, cp9_tr, and postcodes. We'll set
   * seqs_to_aln->tr, seqs_to_aln->cp9_tr, seqs_to_aln->postcode, 
   * to these guys at end of function.
   * 
   * If dsqmode: do not allocate parsetree pointers, they already exist 
   * in search_results.
   */
  tr       = NULL;
  cp9_tr   = NULL;
  postcode = NULL;

  have_parsetrees = (!do_hmmonly && !do_scoreonly && !do_inside) ? TRUE : FALSE;
  if(have_parsetrees)
    ESL_ALLOC(tr, sizeof(Parsetree_t *) * nalign);
  else if(do_hmmonly) /* do_hmmonly */
    ESL_ALLOC(cp9_tr, sizeof(CP9trace_t *) * nalign);
  ESL_ALLOC(parsesc, sizeof(float) * nalign);
  if(do_post) {
    ESL_ALLOC(postcode, sizeof(char **) * nalign);
  }
  if(do_optacc) ESL_ALLOC(parsepp,   sizeof(float) * nalign);
  else          parsepp = NULL;
  if(have_parsetrees) ESL_ALLOC(parse_struct_sc, sizeof(float) * nalign);
  else          parse_struct_sc = NULL;

  minsc =  FLT_MAX;
  maxsc = -FLT_MAX;
  avgsc = 0;
  if((watch = esl_stopwatch_Create()) == NULL) goto ERROR;
  if((watch2 = esl_stopwatch_Create()) == NULL) goto ERROR;

  if(do_hbanded || do_sub) { /* We need a localized CP9 HMM to build sub_cms */
    if(cm->cp9loc == NULL)                 ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments, trying to use CP9 HMM that is NULL\n");
    if(cm->cp9b   == NULL)                 ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments, cm->cp9b is NULL\n");
    if(!(cm->cp9loc->flags & CPLAN9_HASBITS)) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments, trying to use CP9 HMM with CPLAN9_HASBITS flag down.\n");
    
    /* Keep data for the original CM safe; we'll be doing
     * pointer swapping to ease the sub_cm alignment implementation. */
    hmm         = cm->cp9loc;
    cp9b        = cm->cp9b;
    cp9map      = cm->cp9map;
    orig_hmm    = hmm;
    orig_cp9b   = cp9b;
    orig_cp9map = cp9map;
  }
  /* Copy the QD bands in case we expand them. */
  if(do_qdb) {
    if(bdump_level > 1) debug_print_bands(stdout, cm, cm->dmin, cm->dmax);
    expand_flag = FALSE;
    /* Copy dmin and dmax, so we can replace them after expansion */
    ESL_ALLOC(orig_dmin, sizeof(int) * cm->M);
    ESL_ALLOC(orig_dmax, sizeof(int) * cm->M);
    for(v = 0; v < cm->M; v++) {
      orig_dmin[v] = cm->dmin[v];
      orig_dmax[v] = cm->dmax[v];
    }	  
  }
  if(do_sub) { /* to get spos and epos for the sub_cm, 
	        * we config the HMM to local mode with equiprobable start/end points.*/
    swentry = ((hmm->M)-1.)/hmm->M; /* all start pts equiprobable, including 1 */
    swexit  = ((hmm->M)-1.)/hmm->M; /* all end   pts equiprobable, including M */
    CPlan9SWConfig(cm->cp9loc, swentry, swexit, FALSE, cm->ndtype[1]); /* FALSE means don't make I_0, D_1, I_M unreachable (like a local CM, undesirable for sub CM strategy)) */
    CP9Logoddsify(hmm);
  }
  if(! do_hbanded) { /* we need non-banded matrices for alignment */
    mx      = cm_mx_Create(cm);
    out_mx  = cm_mx_Create(cm);
    shmx    = cm_shadow_mx_Create(cm);
    emit_mx = cm_emit_mx_Create(cm);
  }

  orig_cm = cm;
  emap = CreateEmitMap(cm);
  fill_phi_cp9((cm->flags & CMH_LOCAL_BEGIN) ? cm->cp9loc : cm->cp9glb, &phi, 1, TRUE);
  
  /* if not in silent mode, print the header for the sequence info */
  if(!silent_mode) { 
    char *namedashes;
    int ni;
    namewidth = 8; /* length of 'seq name' */
    /* determine the longest name in seqs_to_aln */
    for(ni = 0; ni < seqs_to_aln->nseq; ni++) namewidth = ESL_MAX(namewidth, strlen(seqs_to_aln->sq[ni]->name));
    ESL_ALLOC(namedashes, sizeof(char) * namewidth+1);
    namedashes[namewidth] = '\0';
    for(ni = 0; ni < namewidth; ni++) namedashes[ni] = '-';

    if(cm->align_opts & CM_ALIGN_OPTACC) { 
      if(have_parsetrees) { 
	fprintf(ofp, "#\n");
	fprintf(ofp, "# %9s  %-*s  %5s  %18s  %8s  %11s\n", "",         namewidth, "", "",                  "   bit scores     ",   "",         "");
	fprintf(ofp, "# %9s  %-*s  %5s  %18s  %8s  %11s\n", "",         namewidth, "", "",                  "------------------",   "",         "");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %8s  %8s  %11s\n", "seq idx",    namewidth, "seq name",   "len", "total",    "struct",   "avg prob", "elapsed");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %8s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "--------", "--------", "--------", "-----------");
	if(sfp != NULL) { 
	  fprintf(sfp, "#\n");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %48s  %8s  %11s\n", "",         namewidth, "", "", "",               "",       "                   bit scores                   ",   "",         "");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %48s  %8s  %11s\n", "",         namewidth, "", "", "",               "",       "------------------------------------------------",   "",         "");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %8s  %8s  %8s  %8s  %8s  %8s  %11s\n", "seq idx",    namewidth, "seq name",   "len", "spos",  "epos",  "total",    "primary",  "struct",   "trans",    "null3",    "avg prob", "elapsed");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %8s  %8s  %8s  %8s  %8s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "-----", "-----", "--------", "--------", "--------", "--------", "--------", "--------", "-----------");
	}
      }
      else { 
	fprintf(ofp, "#\n");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %8s  %11s\n", "seq idx",    namewidth, "seq name",   "len",  "bit sc",   "avg prob", "elapsed");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "--------", "--------", "-----------");
	if(sfp != NULL) { 
	  fprintf(sfp, "#\n");
	  fprintf(sfp, "# %9s  %-*s  %5s  %8s  %8s  %11s\n", "seq idx",    namewidth, "seq name",   "len",  "bit sc",   "avg prob", "elapsed");
	  fprintf(sfp, "# %9s  %-*s  %5s  %8s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "--------", "--------", "-----------");
	}
      }
    }
    else { 
      if(have_parsetrees) { 
	fprintf(ofp, "#\n");
	fprintf(ofp, "# %9s  %-*s  %5s  %18s  %11s\n", "",         namewidth, "", "",                  "   bit scores     ",   "");
	fprintf(ofp, "# %9s  %-*s  %5s  %18s  %11s\n", "",         namewidth, "", "",                  "------------------",   "");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %8s  %11s\n", "seq idx",    namewidth, "seq name",   "len", "total",    "struct",   "elapsed");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "--------", "--------", "-----------");
	if(sfp != NULL) { 
	  fprintf(sfp, "#\n");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %48s  %11s\n", "",         namewidth, "", "", "",               "",       "                   bit scores                   ",   "");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %48s  %11s\n", "",         namewidth, "", "", "",               "",       "------------------------------------------------",   "");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %8s  %8s  %8s  %8s  %8s  %11s\n", "seq idx",    namewidth, "seq name",   "len", "spos",  "epos",  "total",    "primary",  "struct",   "trans",    "null3",    "elapsed");
	  fprintf(sfp, "# %9s  %-*s  %5s  %5s  %5s  %8s  %8s  %8s  %8s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "-----", "-----", "--------", "--------", "--------", "--------", "--------", "-----------");
	}
      }
      else { 
	fprintf(ofp, "#\n");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %11s\n",  "seq idx",   namewidth, "seq name",   "len", "bit sc",    "elapsed");
	fprintf(ofp, "# %9s  %-*s  %5s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "--------", "-----------");
	if(sfp != NULL) { 
	  fprintf(sfp, "#\n");
	  fprintf(sfp, "# %9s  %-*s  %5s  %8s  %11s\n",  "seq idx",   namewidth, "seq name",   "len", "bit sc",    "elapsed");
	  fprintf(sfp, "# %9s  %-*s  %5s  %8s  %11s\n",  "---------", namewidth, namedashes, "-----", "--------", "-----------");
	}
      }
    }
    free(namedashes);
  }
  
  /*****************************************************************
   *  Collect parse trees for each sequence
   *****************************************************************/
  for (i = 0; i < nalign; i++) {
    if(!silent_mode) esl_stopwatch_Start(watch);
    cur_dsq = seqs_to_aln->sq[i]->dsq;
    cur_tr  = &(tr[i]);
    L       = seqs_to_aln->sq[i]->n;
    if (L == 0) continue; /* silently skip zero length seqs */


    /* Special case, if do_hmmonly, align seq with Viterbi, print score and move on to next seq */
    if(do_hmmonly) {
      if(!silent_mode) {
	fprintf(ofp, "  %9d  %-*s  %5" PRId64, (iidx+i), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n);
	if(sfp != NULL) fprintf(sfp, "  %9d  %-*s  %5" PRId64, (iidx+i), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n);
      }
      if((status = cp9_Viterbi((cm->flags & CMH_LOCAL_BEGIN) ? cm->cp9loc : cm->cp9glb,
			       errbuf, cm->cp9_mx, cur_dsq, 1, L, 
			       FALSE,  /* we are not scanning */
			       TRUE,   /* we are aligning */
			       FALSE,  /* don't be memory efficient */
			       NULL, NULL, /* don't return best sc at each posn, or best scoring posn */
			       &(cp9_tr[i]), /* return the trace */
			       &sc)) != eslOK) return status;
      null3_correction = 0.;
      if(do_null3) ScoreCorrectionNull3CompUnknown(cm->abc, cm->null, cur_dsq, 1, L, cm->null3_omega, &null3_correction);
      if(!silent_mode) { 
	esl_stopwatch_Stop(watch); 
	FormatTimeString(time_buf, watch->user, TRUE);
	fprintf(ofp, "  %8.2f  %11s\n", sc - null3_correction, time_buf);
      }
      parsesc[i] = sc;
      continue;
    }
    /* Special case, if do_scoreonly, align seq with full CYK inside, just to 
     * get the score. For testing, probably in cmscore. */
    if(do_scoreonly) {
      if(!silent_mode) {
	fprintf(ofp, "  %9d  %-*s  %5" PRId64, (iidx+i), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n);
	if(sfp != NULL) { 
	  fprintf(ofp, "  %9d  %-*s  %5" PRId64, (iidx+i), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n);
	}
      }
      sc = CYKInsideScore(cm, cur_dsq, L, 0, 1, L, NULL, NULL); /* don't do QDB mode */
      if(!silent_mode) fprintf(ofp, "  %8.2f  ", sc);
      parsesc[i] = sc;
      continue;
    }

    /* Potentially, do HMM calculations. */
    if((!do_sub) && do_hbanded) {
      if(do_p7banded) { 
	if((status =  p7_Seq2Bands   (orig_cm, errbuf, gm, gx, bg, p7_tr, cur_dsq, L, phi, sc7, len7, end7, mprob7, mcprob7, iprob7, ilprob7, pad7, &p7_i2k, &kmin, &kmax, &i_ncells_banded)) != eslOK) return status;
	if((status = cp9_Seq2BandsP7B(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, cur_dsq, L, orig_cp9b, kmin, kmax, debug_level)) != eslOK) return status; 
      }
      else { 
	if((status = cp9_Seq2Bands(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, cur_dsq, 1, L, orig_cp9b, FALSE, tro, debug_level)) != eslOK) return status; 
      }
    }
    else if(do_sub) { 
      /* If we're in sub mode:
       * (0) Get p7 bands (only if do_p7banded) 
       * (1) Get HMM posteriors 
       * (2) Infer the start (spos) and end (epos) HMM states by 
       *     looking at the posterior matrix.
       * (3) Build the sub_cm from the original CM.
       *
       * If we're also doing HMM banded alignment to sub CM:
       * (4) Build a new CP9 HMM from the sub CM.
       * (5) Do Forward/Backward again, and get HMM bands 
       */
      
      if(do_p7banded) { /* use P7 bands to constrain CP9 dp calculations */
	/* (0) Get p7 bands */
	if((status = p7_Seq2Bands(orig_cm, errbuf, gm, gx, bg, p7_tr, cur_dsq, L, phi, sc7, len7, end7, mprob7, mcprob7, iprob7, ilprob7, pad7, &p7_i2k, &kmin, &kmax, &i_ncells_banded)) != eslOK) return status;
	/* (1) Get HMM posteriors */
	if((status = cp9_Seq2PosteriorsP7B(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, cur_dsq, L, kmin, kmax, debug_level)) != eslOK) return status; 
	/* (2) infer the start and end HMM nodes (consensus cols) from posterior matrix.
	 * Remember: we're necessarily in CP9 local mode, the --sub option turns local mode on. 
	 */
	if((status = CP9NodeForPosnP7B(orig_hmm, errbuf, 1, orig_cm->cp9_bmx, kmin[1], kmax[1], &spos, &spos_state, debug_level)) != eslOK) return status;
	if((status = CP9NodeForPosnP7B(orig_hmm, errbuf, L, orig_cm->cp9_bmx, kmin[L], kmax[L], &epos, &epos_state, debug_level)) != eslOK) return status;
      }
      else { /* ! do_p7banded */ 
	/* (1) Get HMM posteriors */
	if((status = cp9_Seq2Posteriors(orig_cm, errbuf, orig_cm->cp9_mx, orig_cm->cp9_bmx, orig_cm->cp9_bmx, cur_dsq, 1, L, debug_level)) != eslOK) return status; 
	
	/* (2) infer the start and end HMM nodes (consensus cols) from posterior matrix.
	 * Remember: we're necessarily in CP9 local mode, the --sub option turns local mode on. 
	 */
	CP9NodeForPosn(orig_hmm, 1, L, 1, orig_cm->cp9_bmx, &spos, &spos_state, 0., TRUE, debug_level);
	CP9NodeForPosn(orig_hmm, 1, L, L, orig_cm->cp9_bmx, &epos, &epos_state, 0., FALSE, debug_level);
      }

      /* Deal with special cases for sub-CM alignment:
       * If the most likely state to have emitted the first or last residue
       * is the insert state in node 0, it only makes sense to start modelling
       * at consensus column 1. */
      if(spos == 0 && spos_state == 1) spos = 1;
      if(epos == 0 && epos_state == 1) epos = 1;
      /* If most-likely HMM node to emit final position comes BEFORE or EQUALS the most-likely HMM node to emit first position,
       * our HMM alignment is crap, default to using the full CM. (note: If EQUALS we could be right, but we can't build a
       * CM from a single consensus column (see notes in cm_modelmaker.c::cm_from_guide), and I would argue we don't really care about
       * getting single residue alignments correct anyway. */
      if(epos <= spos) { spos = 1; epos = orig_hmm->M; } 

      /* (3) Build the sub_cm from the original CM. */
      esl_stopwatch_Start(watch2);
      if((status = build_sub_cm(orig_cm, errbuf, &sub_cm, 
				spos, epos,               /* first and last col of structure kept in the sub_cm  */
				&submap,                  /* maps from the sub_cm to cm and vice versa           */
				debug_level)) != eslOK)    /* print or don't print debugging info                 */
	return status;
      esl_stopwatch_Stop(watch2);
      FormatTimeString(time_buf, watch2->user, TRUE);
#if PRINTNOW
      fprintf(stdout, "build sub cm      %11s\n", time_buf);
#endif

      /* Configure the sub_cm, the same as the cm, this will build a CP9 HMM if (do_hbanded), this will also:  */
      /* (4) Build a new CP9 HMM from the sub CM. */
      esl_stopwatch_Start(watch2);
      if((status = ConfigCM(sub_cm, errbuf, FALSE, orig_cm, submap)) != eslOK) return status; /* FALSE says: don't calculate W, we won't need it */
      esl_stopwatch_Stop(watch2);
      FormatTimeString(time_buf, watch2->user, TRUE);
#if PRINTNOW
      fprintf(stdout, "config sub cm      %11s\n", time_buf);
#endif
      cm    = sub_cm; /* orig_cm still points to the original CM */
      if(do_hbanded) { /* we're doing HMM banded alignment to the sub_cm */
	/* Get the HMM bands for the sub_cm */
	sub_hmm    = (sub_cm->flags & CMH_LOCAL_BEGIN) ? sub_cm->cp9loc : sub_cm->cp9glb;
	sub_cp9b   = sub_cm->cp9b;
	sub_cp9map = sub_cm->cp9map;

	if(do_p7banded) { 
	  P7BandsAdjustForSubCM(kmin, kmax, L, spos, epos);
	  if((status = cp9_Seq2BandsP7B(sub_cm, errbuf, sub_cm->cp9_mx, sub_cm->cp9_bmx, sub_cm->cp9_bmx, cur_dsq, L, sub_cp9b, kmin, kmax, debug_level)) != eslOK) return status;
	}
	else { /* ! do_p7banded */ 
	  /* (5) Do Forward/Backward again, and get HMM bands */
	  if((status = cp9_Seq2Bands(sub_cm, errbuf, sub_cm->cp9_mx, sub_cm->cp9_bmx, sub_cm->cp9_bmx, cur_dsq, 1, L, sub_cp9b, FALSE, tro, debug_level)) != eslOK) return status;
	}
	hmm           = sub_hmm;    
	cp9b          = sub_cp9b;
	cp9map        = sub_cp9map;
      }
    }

    if(do_qdb) {
      /*Check if we need to reset the query dependent bands b/c they're currently expanded. */
      if(expand_flag) {
	for(v = 0; v < cm->M; v++) { cm->dmin[v] = orig_dmin[v]; cm->dmax[v] = orig_dmax[v]; }
	expand_flag = FALSE;
      }
      if((L < cm->dmin[0]) || (L > cm->dmax[0])) { 
	/* the seq we're aligning is outside the root band, so we expand.*/
	ExpandBands(cm, L, cm->dmin, cm->dmax);
	if(debug_level > 0) fprintf(ofp, "# Expanded bands for seq : %s\n", seqs_to_aln->sq[i]->name);
	if(bdump_level > 2) { fprintf(ofp, "printing expanded bands :\n"); debug_print_bands(ofp, cm, cm->dmin, cm->dmax); }
	expand_flag = TRUE;
      }
    }

    if(!silent_mode) { 
      fprintf(ofp, "  %9d  %-*s  %5" PRId64, (iidx+i), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n);
      if(sfp != NULL) fprintf(sfp, "  %9d  %-*s  %5" PRId64, (iidx+i), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n);
    }

    esl_stopwatch_Start(watch2);  

    /* beginning of large if() else if() else if() ... statement */
    if(do_inside) { 
      if(do_hbanded) { /* HMM banded inside only */
	if((status = cm_InsideAlignHB(cm, errbuf, cur_dsq, L, size_limit, cm->hbmx, &sc)) != eslOK) return status; /* errbuf will have been filled by cm_InsideAlignHB() */
      }
      else { /* non-banded inside only */
	if((status = cm_InsideAlign(cm, errbuf, cur_dsq, L, size_limit, mx, &sc)) != eslOK) return status; 
      }
    }
    else if (do_small) { /* small D&C CYK alignment */
      if(do_qdb) { /* use QDBs when doing D&C CYK */
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, cm->dmin, cm->dmax);
	if(bdump_level > 0) qdb_trace_info_dump(cm, *cur_tr, cm->dmin, cm->dmax, bdump_level);
      }
      else { /* small D&C CYK non-banded alignment */
	sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL); /* we're not in QDB mode */
	if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level); /* informative as to where optimal parse goes outside bands */
      }
    }
    else if(do_qdb) { /* non-small, QDB banded CYK alignment */
      sc = CYKInside(cm, cur_dsq, L, 0, 1, L, cur_tr, cm->dmin, cm->dmax);
      if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level);
    }
    else { /* non-small, non-QDB CYK or optimal accuracy alignment or sample an alignment from Inside matrix */
      if(do_hbanded) { /* HMM banded CYK, optimal accuracy or sample, either with or without posteriors annotated */
	status = cm_AlignHB(cm, errbuf, cur_dsq, L, size_limit, do_optacc, do_sample, cm->hbmx, cm->shhbmx, cm->ohbmx, cm->ehbmx,
			    ( do_sample             ? r              : NULL), 
			    ( do_post               ? &(postcode[i]) : NULL), 
			    ((do_post || do_optacc) ? &ins_sc        : NULL), 
			    cur_tr, &sc);
	if (status == eslERANGE && (! do_optacc) && (! do_sample) && (! do_post)) { 
	  ESL_DPRINTF1(("# Doing D&C because HMM banded parse of seq %d was too memory intensive.\n", i));
	  fprintf(ofp, "# Doing D&C because HMM banded parse of seq %d was too memory intensive.\n", i); 
	  sc = CYKDivideAndConquer(cm, cur_dsq, L, 0, 1, L, cur_tr, NULL, NULL); /* we're not in QDB mode */
	}
	else if(status != eslOK) return status;
      } /* end of if(do_hbanded) */
      else { /* non-banded */
	status = cm_Align(cm, errbuf, cur_dsq, L, size_limit, do_optacc, do_sample, mx, shmx, out_mx, emit_mx,
			  ( do_sample             ? r              : NULL), 
			  ( do_post               ? &(postcode[i]) : NULL), 
			  ((do_post || do_optacc) ? &ins_sc        : NULL), 
			  cur_tr, &sc);
	if(status != eslOK) return status;
	if(bdump_level > 0) qdb_trace_info_dump(cm, tr[i], cm->dmin, cm->dmax, bdump_level); /* allows you to see where the non-banded parse went outside the QD bands. */
      } 
    }
    /* end of large if() else if() else if() else statement */
    /* done alignment for this seq */

    struct_sc = IMPOSSIBLE;
    if(have_parsetrees) { 
      if(*cur_tr == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments() should have parsetrees, but don't (coding error).");
      cur_emap = CreateEmitMap(cm); /* TEMPORARY!, we don't want to build this separately for each seq, but with sub CMs I'm not sure how to accomplish that */
      if((status = ParsetreeScore(cm, cur_emap, errbuf, *cur_tr, cur_dsq, FALSE, &tmpsc, &struct_sc, &primary_sc, &cur_spos, &cur_epos)) != eslOK) return status;
      /* if(do_sub) { cur_spos += (spos-1); cur_epos += (spos-1); } */
      /* currently primary_sc, cur_spos, cur_epos are not output, but they may be someday */
      FreeEmitMap(cur_emap);
    }
    /* determine NULL3 score correction, which is independent of the parsetree */
    null3_correction = 0.;
    if(do_null3) ScoreCorrectionNull3CompUnknown(cm->abc, cm->null, cur_dsq, 1, L, cm->null3_omega, &null3_correction);
    /* EPN, Thu Dec 17 14:07:04 2009
     * Following code would correct for null3 in struct_sc and primary_sc,
     * I'm undecided as to whether that's a good idea. 
     * Currently no correction is performed, but this code is left here for reference.
     *
     * nbp_emits   = (float) ParsetreeCountMPEmissions(cm, *cur_tr);
     * adjust struct_sc  for null3 (this is inexact) 
     * struct_sc  -= (nbp_emits / float (L)) * null3_correction;                
     * adjust primary_sc  for null3 (this is inexact) 
     * primary_sc -= (((float) L - nbp_emits) / float (L)) * null3_correction;
     */
    if(!silent_mode) { 
      if(have_parsetrees) { 
 
	if(do_optacc)  fprintf(ofp, "  %8.2f  %8.2f  %8.3f  ", ins_sc - null3_correction, struct_sc, sc);
	else           fprintf(ofp, "  %8.2f  %8.2f  ",            sc - null3_correction, struct_sc);
	if(sfp != NULL) { 
	  if(do_optacc)  fprintf(sfp, "  %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.3f  ", cur_spos, cur_epos, ins_sc - null3_correction, primary_sc, struct_sc, (ins_sc - primary_sc - struct_sc), -1. * null3_correction, sc);
	  else           fprintf(sfp, "  %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  ",        cur_spos, cur_epos,     sc - null3_correction, primary_sc, struct_sc, (    sc - primary_sc - struct_sc), -1. * null3_correction);
	}
      }	
      else { 
	if(do_optacc)  fprintf(ofp, "  %8.2f  %8.3f  ", ins_sc - null3_correction, sc);
	else           fprintf(ofp, "  %8.2f  ", sc - null3_correction);
	if(sfp != NULL) { 
	  if(do_optacc)  fprintf(sfp, "  %8.2f  %8.3f  ", ins_sc - null3_correction, sc);
	  else           fprintf(sfp, "  %8.2f  ", sc - null3_correction);
	}
      }
    }
    /* if we did optimally accurate alignment:
     *  - <sc> is average posterior probability of optimally accurate parse, 
     *  - <ins_sc> is the Inside score for the target sequence,
     * else (non-optimally accurate alignment), for ex, CYK:
     *  - <sc> is parse score in bits
     *  - <ins_sc> is irrelevant 
     */
    if(do_optacc) { 
      parsesc[i] = ins_sc; 
      parsepp[i] = sc; /* sc will be an average posterior probability if do_optacc */
      if(parse_struct_sc != NULL) parse_struct_sc[i] = struct_sc;
      avgsc += ins_sc;
      if (ins_sc > maxsc) maxsc = ins_sc;
      if (ins_sc < minsc) minsc = ins_sc;
    }
    else { 
      parsesc[i] = sc; /* parsepp == NULL if !do_optacc */
      if(parse_struct_sc != NULL) parse_struct_sc[i] = struct_sc;
      avgsc += sc;
      if (sc > maxsc) maxsc = sc;
      if (sc < minsc) minsc = sc;
    }
  
    /* check parsetree score if cm->align_opts & CM_ALIGN_CHECKPARSESC */
    if((cm->align_opts & CM_ALIGN_CHECKPARSESC) && (!(cm->flags & CM_IS_SUB))) { 
      if(do_optacc) 
	ESL_FAIL(eslEINCOMPAT, errbuf, "DispatchAlignments(), align_opts CM_ALIGN_CHECKPARSESC is on, but incompatible with raised flag CM_ALIGN_OPTACC.\n");
      if((status = ParsetreeScore(cm, NULL, errbuf, tr[i], cur_dsq, FALSE, &tmpsc, NULL, NULL, NULL, NULL)) != eslOK) return status;
      if (fabs(sc - tmpsc) >= 0.03)
	ESL_FAIL(eslFAIL, errbuf, "DispatchAlignments(), seq: %d alignment score %.3f differs from its parse tree's score: %.3f", i, sc, tmpsc);
    }

    /* If requested, or if debug level high enough, print out the parse tree */
    if((cm->align_opts & CM_ALIGN_PRINTTREES) || (debug_level > 2)) { 
      if((status = ParsetreeScore(cm, NULL, errbuf, tr[i], cur_dsq, FALSE, &tmpsc, &struct_sc, &primary_sc, NULL, NULL)) != eslOK) return status;
      fprintf(ofp, "  %16s %.2f bits\n", "SCORE:", tmpsc);
      fprintf(ofp, "  %16s %.2f bits\n", "STRUCTURE SCORE:", struct_sc);
      fprintf(ofp, "  %16s %.2f bits\n", "PRIMARY   SCORE:", primary_sc);
      ParsetreeDump(stdout, tr[i], cm, cur_dsq, NULL, NULL);
      fprintf(ofp, "//\n");
    }
    /* Dump the trace with info on i, j and d bands
     * if bdump_level is high enough */
    if(bdump_level > 0 && do_hbanded)
      ijdBandedTraceInfoDump(cm, tr[i], cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 1);

    if(do_sub) {
      if(! do_inside) { 
	/* Convert the sub_cm parsetree to a full CM parsetree */
	if(debug_level > 0) ParsetreeDump(stdout, *cur_tr, cm, cur_dsq, NULL, NULL);
	if(!(sub_cm2cm_parsetree(orig_cm, sub_cm, &orig_tr, *cur_tr, submap, debug_level))) { 
	  /* ParsetreeDump(stdout, orig_tr, orig_cm, cur_dsq, NULL, NULL); */
	  ESL_FAIL(eslFAIL, errbuf, "DispatchAlignments(), Unable to convert sub CM parsetree to original CM parsetree. This shouldn't happen.");
	}
	if(debug_level > 0) { 
	  fprintf(ofp, "\n\nConverted original trace:\n");
	  ParsetreeDump(stdout, orig_tr, orig_cm, cur_dsq, NULL, NULL);
	}
	/* Replace the sub_cm trace with the converted orig_cm trace. */
	FreeParsetree(*cur_tr);
	*cur_tr = orig_tr;
      }

      FreeSubMap(submap);
      FreeCM(sub_cm); /* cm and sub_cm now point to NULL */
    }
    if(!silent_mode) { 
      esl_stopwatch_Stop(watch); 
      FormatTimeString(time_buf, watch->user, TRUE);
      fprintf(ofp, "%11s\n", time_buf);
      if(sfp != NULL) fprintf(sfp, "%11s\n", time_buf);
    }
    esl_stopwatch_Stop(watch2); 
    FormatTimeString(time_buf, watch2->user, TRUE);
#if PRINTNOW
    fprintf(stdout, "CM parse        %11s\n", time_buf);
#endif

#if 0    
    if(do_p7banded) { 
      /* compare parsetree with p7 msv alignment */
      if((status = Parsetree2i_to_k(cm, emap, seqs_to_aln->sq[i]->n, errbuf, *cur_tr, &cm_i2k)) != eslOK) return status;;
      i_ncells_total = (double) (seqs_to_aln->sq[i]->n * cm->clen);
      ncells_total  += i_ncells_total;
      ncells_banded += (double) i_ncells_banded;
      
      i_ncorrect = 0;
      i_npins = 0;
      i_nodes_n = i_nodes_ncorrect = 0;
      float mprob = 1.;
      int posn = 0;
      for(ipos = 1; ipos <= seqs_to_aln->sq[i]->n; ipos++) { 
	if(p7_i2k[ipos] != -1) { 
	  /*printf("\t%4d\t%4d\t%4d\t%1s\n", ipos, p7_i2k[ipos], cm_i2k[ipos], ((abs(p7_i2k[ipos]-cm_i2k[ipos])) == 0) ? "+" : "-");*/
	  mprob *= phi[p7_i2k[ipos]][HMMMATCH];
	  posn++;
	  if(p7_i2k[ipos] == cm_i2k[ipos]) { 
	    i_ncorrect++;
	    ;/*printf("\t%4d\t%4d\t%4d\t\t%3d\t\t%.3f\t\t%.3f\t%.3f\t%.3f\t\t%.3f\t%.3f\t%.3f\t%1s\n", ipos, p7_i2k[ipos], cm_i2k[ipos], posn, mprob, phi[(ipos-1)][HMMMATCH], phi[(ipos-1)][HMMDELETE], phi[(ipos-1)][HMMINSERT], phi[ipos][HMMMATCH], phi[ipos][HMMDELETE], phi[ipos][HMMINSERT], ((abs(p7_i2k[ipos]-cm_i2k[ipos])) == 0) ? "+" : "-");*/
	  }
	  else { 
	    ;printf("\t%4d\t%4d\t%4d\t\t%3d\t\t%.3f\t\t%.3f\t%.3f\t%.3f\t\t%.3f\t%.3f\t%.3f\t%1s\n", ipos, p7_i2k[ipos], cm_i2k[ipos], posn, mprob, phi[(ipos-1)][HMMMATCH], phi[(ipos-1)][HMMDELETE], phi[(ipos-1)][HMMINSERT], phi[ipos][HMMMATCH], phi[ipos][HMMDELETE], phi[ipos][HMMINSERT], ((abs(p7_i2k[ipos]-cm_i2k[ipos])) == 0) ? "+" : "-");
	  }
	  i_npins++;
	}
	else { mprob = 1.; posn = 0; }
	cm_k = cm_i2k[ipos];
	if(cm_k < 0) cm_k *= -1;
	if(kmin[ipos] <= cm_k && kmax[ipos] >= cm_k) i_nodes_ncorrect++;
	else { 
	  ;printf("\tNODE OFF %4d\t%4d\t%4d\t\t%4d\t\t%4d\n", ipos, p7_i2k[ipos], cm_i2k[ipos], kmin[ipos], kmax[ipos]);
	}
	i_nodes_n++;
    }
      
      ncorrect += i_ncorrect;
      npins    += i_npins;
      nodes_ncorrect += i_nodes_ncorrect;
      nodes_n        += i_nodes_n;
      
      printf("> %5d %d/%d (%.3f) MSV p7 pins correct.\n", i, i_ncorrect, i_npins, (float) i_ncorrect / (float) i_npins);
      printf("> %5d %d/%d (%.3f) MSV p7 nodes potentially correct.\n", i, i_nodes_ncorrect, i_nodes_n, (float) i_nodes_ncorrect / (float) i_nodes_n);
      printf("> %5d %.0f/%.0f (%.8f) matrix pruned away.\n", i, (float) i_ncells_banded, i_ncells_total, 1. - (i_ncells_banded / i_ncells_total));
      free(cm_i2k);
    }
    printf("\n");
    printf("$ TOTAL %d/%d (%.3f) MSV p7 pins correct.\n", ncorrect, npins, (float) ncorrect / (float) npins);
    printf("$ TOTAL %d/%d (%.3f) consensus nodes potentially correct.\n", nodes_ncorrect, nodes_n, (float) nodes_ncorrect / (float) nodes_n);
    printf("$ TOTAL %.0f/%.0f (%.8f) matrix pruned away.\n", ncells_banded, ncells_total, 1. - (ncells_banded / ncells_total));
#endif
    if(do_p7banded) { 
      free(p7_i2k);
      free(kmin);
      free(kmax);
    }
  }

  /* done aligning all nalign seqs. */
  /* Clean up. */
  if(mx != NULL)      cm_mx_Destroy(mx);
  if(out_mx != NULL)  cm_mx_Destroy(out_mx);
  if(shmx != NULL)    cm_shadow_mx_Destroy(shmx);
  if(emit_mx != NULL) cm_emit_mx_Destroy(emit_mx);
  if (do_qdb) {
    free(orig_dmin);
    free(orig_dmax);
  }
  esl_stopwatch_Destroy(watch);
  esl_stopwatch_Destroy(watch2);
  
  seqs_to_aln->tr       = tr;       /* could be NULL */
  seqs_to_aln->cp9_tr   = cp9_tr;   /* could be NULL */
  seqs_to_aln->postcode = postcode; /* could be NULL */
  seqs_to_aln->sc       = parsesc;  /* shouldn't be NULL */
  seqs_to_aln->pp       = parsepp;  /* could be NULL */
  seqs_to_aln->struct_sc= parse_struct_sc; /* could be NULL */

  p7_trace_Destroy(p7_tr);
  p7_bg_Destroy(bg);
  p7_profile_Destroy(gm);
#if 0 
  p7_omx_Destroy(ox);
#endif
  p7_gmx_Destroy(gx);
  FreeEmitMap(emap);
  if(phi != NULL) { 
    for(k = 0; k <= orig_cm->cp9loc->M; k++) free(phi[k]);
    free(phi);
  }

  return eslOK;
  ERROR:
  ESL_FAIL(eslEMEM, errbuf, "DispatchAlignments(), Memory allocation error.");
  return status; /* NEVERREACHED */
}


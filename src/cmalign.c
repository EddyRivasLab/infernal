/* cmalign.c
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 * SVN $Id$
 * 
 * Align sequences to a CM.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "hmmer_funcs.h"
#include "hmmer_structs.h"
#include "sre_stack.h"
#include "hmmband.h"         
#include "cm_postprob.h"

#define BE_EFFICIENT  0		/* setting for do_full: small memory mode */
#define BE_PARANOID   1		/* setting for do_full: keep whole matrix, perhaps for debugging */

static void ExpandBands(CM_t *cm, int qlen, int *dmin, int *dmax);
static void banded_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, int *dmax, int bdump_level);

static char banner[] = "cmalign - align sequences to an RNA CM";

static char usage[]  = "\
Usage: cmalign [-options] <cmfile> <sequence file>\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -l     : local; align locally w.r.t. the model\n\
   -o <f> : output the alignment file to file <f>\n\
   -q     : quiet; suppress verbose banner\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input alignment is in format <s>\n\
   --nosmall     : use normal alignment algorithm, not d&c\n\
   --regress <f> : save regression test data to file <f>\n\
   --full        : include all match columns in output alignment\n\
   --tfile <f>   : dump individual sequence tracebacks to file <f>\n\
   --banddump <n>: set verbosity of band info print statements to <n> [1..3]\n\
   --dlev <n>    : set verbosity of debugging print statements to <n> [1..3]\n\
   --time        : print timings for alignment, band calculation, etc.\n\
   --inside      : don't align; return scores from the Inside algorithm \n\
   --outside     : don't align; return scores from the Outside algorithm\n\
   --post        : align with CYK and append posterior probabilities\n\
   --checkpost   : check that posteriors are correctly calc'ed\n\
   --sub         : build sub CM for columns b/t HMM predicted start/end points\n\
   --fsub        : build sub CM for structure b/t HMM predicted start/end points\n\
\n\
  * HMM banded alignment related options:\n\
   --hbanded     : use exptl CM plan 9 HMM banded CYK aln algorithm\n\
   --hbandp <f>  : tail loss prob for --hbanded [default: 0.0001]\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
   --checkcp9    : check the CP9 empirically by generating sequences\n\
\n\
  * Query dependent banded (qdb) alignment related options:\n\
   --qdb         : use query dependent banded CYK alignment algorithm\n\
   --beta <f>    : tail loss prob for --qdb [default:1E-7]\n\
   --expand      : naively expand qd bands if target seq is outside root band\n\
    -W <n>       : window size for calc'ing qd bands (df: precalc'd in cmbuild)\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-l", TRUE, sqdARG_NONE },
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--nosmall",    FALSE, sqdARG_NONE },
  { "--regress",    FALSE, sqdARG_STRING },
  { "--qdb",     FALSE, sqdARG_NONE },
  { "--beta",      FALSE, sqdARG_FLOAT},
  { "--expand", FALSE, sqdARG_NONE},
  { "--tfile",     FALSE, sqdARG_STRING },
  { "--banddump"  , FALSE, sqdARG_INT},
  { "-W", TRUE, sqdARG_INT },
  { "--full", FALSE, sqdARG_NONE },
  { "--dlev",       FALSE, sqdARG_INT },
  { "--hbanded",    FALSE, sqdARG_NONE },
  { "--hbandp",     FALSE, sqdARG_FLOAT},
  { "--sums",       FALSE, sqdARG_NONE},
  { "--time",       FALSE, sqdARG_NONE},
  { "--inside",     FALSE, sqdARG_NONE},
  { "--outside",    FALSE, sqdARG_NONE},
  { "--post",       FALSE, sqdARG_NONE},
  { "--checkpost",  FALSE, sqdARG_NONE},
  { "--sub",        FALSE, sqdARG_NONE},
  { "--fsub",       FALSE, sqdARG_NONE},
  { "--checkcp9",   FALSE, sqdARG_NONE},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  CMFILE          *cmfp;	/* open CM file */
  char            *seqfile;     /* file to read sequences from */
  int              format;      /* format of sequence file */
  CM_t            *cm;          /* a covariance model       */
  char           **rseq;        /* RNA sequences to align */
  int              nseq;	/* number of rseq */
  SQINFO          *sqinfo;      /* optional info attached to sequences */
  char           **dsq;         /* digitized RNA sequences */
  Parsetree_t    **tr;          /* parse trees for the sequences */
  MSA             *msa;         /* alignment that's created */
  char            *outfile;	/* optional output file name */
  FILE            *ofp;         /* an open output file */
  int              i;		/* counter over sequences/parsetrees */
  int              v;           /* counter over states of the CM */
  int              k;           /* counter over HMM nodes */

  float            sc;		/* score for one sequence alignment */
  float            maxsc;	/* max score in all seqs */
  float            minsc;	/* min score in all seqs */
  float            avgsc;	/* avg score over all seqs */
  
  char  *regressfile;           /* regression test data file */
  char  *tracefile;		/* file to dump debugging traces to        */
  int    be_quiet;		/* TRUE to suppress verbose output & banner */
  int    do_local;		/* TRUE to config the model in local mode   */
  int    do_small;		/* TRUE to do divide and conquer alignments */
  int    windowlen;             /* window length for calculating bands */
  int    do_full;               /* TRUE to output all match columns in output alignment */
  int    do_qdb;                /* TRUE to do qdb CYK (either d&c or full)  */
  int    bdump_level;           /* verbosity level for --banddump option, 0 is OFF */

  char  *optname;               /* name of option found by Getopt()        */
  char  *optarg;                /* argument found by Getopt()              */
  int    optind;                /* index in argv[]                         */
  int    safe_windowlen;        /* initial windowlen (W) used for calculating bands
				 * in BandCalculationEngine(). For truncation error 
				 * handling this should be > than expected dmax[0]*/
  int    debug_level;           /* verbosity level for debugging printf() statements,
			         * passed to many functions. */
  Stopwatch_t  *watch1;         /* for overall timings */
  Stopwatch_t  *watch2;         /* for HMM band calc timings */
  int    time_flag;             /* TRUE to print timings, FALSE not to */

  /* query dependent bands data structures */
  double   qdb_beta;	        /* tail loss probability for query dependent banding */
  int    do_expand;             /* TRUE to naively expand query dependent bands when necessary */
  int    expand_flag;           /* TRUE if the dmin and dmax vectors have just been 
				 * expanded (in which case we want to recalculate them 
				 * before we align a new sequence), and FALSE if not*/
  int      set_window;          /* TRUE to set window length due to -W option*/
  double **gamma;               /* P(subseq length = n) for each state v    */
  int     *dmin;                /* minimum d bound for state v, [0..v..M-1] */
  int     *dmax;                /* maximum d bound for state v, [0..v..M-1] */

  /* HMMERNAL!: hmm banded alignment data structures */
  /* data structures for hmm bands (bands on the hmm states) */
  int     *pn_min_m;          /* HMM band: minimum position node k match state will emit  */
  int     *pn_max_m;          /* HMM band: maximum position node k match state will emit  */
  int     *pn_min_i;          /* HMM band: minimum position node k insert state will emit */
  int     *pn_max_i;          /* HMM band: maximum position node k insert state will emit */
  int     *pn_min_d;          /* HMM band: minimum position node k delete state will emit */
  int     *pn_max_d;          /* HMM band: maximum position node k delete state will emit */
  double   hbandp;            /* tail loss probability for hmm bands */
  int      use_sums;          /* TRUE to fill and use the posterior sums, false not to. */
  int    **isum_pn_m;         /* [1..k..M] sum over i of log post probs from post->mmx[i][k]*/
  int    **isum_pn_i;         /* [1..k..M] sum over i of log post probs from post->imx[i][k]*/
  int    **isum_pn_d;         /* [1..k..M] sum over i of log post probs from post->dmx[i][k]*/

  /* data structures for mapping HMM to CM */
  int *node_cc_left;    /* consensus column each CM node's left emission maps to
			 * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *node_cc_right;   /* consensus column each CM node's right emission maps to
			 * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *cc_node_map;     /* node that each consensus column maps to (is modelled by)
			 * [1..hmm_nmc] */
  int **cs2hn_map;      /* 2D CM state to HMM node map, 1st D - CM state index
		         * 2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
		         * that contains state that maps to CM state, -1 if none.*/
  int **cs2hs_map;      /* 2D CM state to HMM node map, 1st D - CM state index
		         * 2nd D - 2 elements for up to 2 matching HMM states, 
		         * value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state,
		         * -1 if none.
		         * For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
                         * maps to CM state v.*/
  int ***hns2cs_map;    /* 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
		         * HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
		         * 2 matching CM states, value: CM states that map, -1 if none.
		         * For example: CM states hsn2cs_map[k][0][0] and hsn2cs_map[k][0][1]
		         * map to HMM node k's match state.*/

  /* arrays for CM state bands, derived from HMM bands */
  int  do_hbanded;      /* TRUE to do CM Plan 9 HMM banded CYKInside_b_jd() using bands on d and j dim*/
  int *imin;            /* [1..M] imin[v] = first position in band on i for state v*/
  int *imax;            /* [1..M] imax[v] = last position in band on i for state v*/
  int *jmin;            /* [1..M] jmin[v] = first position in band on j for state v*/
  int *jmax;            /* [1..M] jmax[v] = last position in band on j for state v*/
  int **hdmin;          /* [v=1..M][0..(jmax[v]-jmin[v])] 
			 * hdmin[v][j0] = first position in band on d for state v, and position
			 * j = jmin[v] + j0.*/
  int **hdmax;          /* [v=1..M][0..(jmax[v]-jmin[v])] 
			 * hdmin[v][j0] = last position in band on d for state v, and position
			 * j = jmin[v] + j0.*/
  int *safe_hdmin;      /* [1..M] safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0) */
  int *safe_hdmax;      /* [1..M] safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0) */

  /* data structures for HMMs */
  char  *hmmfile;      /* file to read HMMs from                  */
  int   ks;            /* Counter over HMM state types (0 (match), 1(ins) or 2 (del))*/
  float forward_sc; 
  float backward_sc; 
  int   hmm_M;         /* Number of nodes in either the Plan 7 HMM or CM Plan 9 HMM */
  
  /* CM Plan 9 */
  CP9HMMFILE            *hmmfp;     /* opened CP9 hmmfile for reading                       */
  struct cplan9_s       *hmm;       /* constructed CP9 HMM; written to hmmfile              */
  int                    cp9_M;         /* number of nodes in CP9 HMM (MATL+MATR+2*MATP)        */
  struct cp9_dpmatrix_s *cp9_mx;        /* growable DP matrix for viterbi                       */
  struct cp9_dpmatrix_s *cp9_fwd;       /* growable DP matrix for forward                       */
  struct cp9_dpmatrix_s *cp9_bck;       /* growable DP matrix for backward                      */
  struct cp9_dpmatrix_s *cp9_posterior; /* growable DP matrix for posterior decode              */
  float                  swentry;	/* S/W aggregate entry probability       */
  float                  swexit;        /* S/W aggregate exit probability        */
  int                    do_checkcp9;   /* TRUE to check the CP9 HMM by generating sequences */
  int                    seed;	        /* random number generator seed (only used if(do_checkcp9) */

  unsigned char    **p7dsq;     /* digitized RNA sequences for the HMM */

  /* Alternatives to CYK */
  int                do_inside; /* TRUE to use the Inside algorithm instead of CYK */
  int                do_outside;/* TRUE to use the Outside algorithm instead of CYK */
  int                do_check;  /* TRUE to check Inside and Outside probabilities */
  int                do_post;   /* TRUE to use the Outside algorithm instead of CYK */
  char             **postcode;  /* posterior decode array of strings        */
  char              *apostcode; /* aligned posterior decode array           */
  float           ***alpha;     /* alpha DP matrix for Inside() */
  float           ***beta;      /* beta DP matrix for Inside() */
  float           ***post;      /* post DP matrix for Inside() */
  int j;

  /* the --sub option */
  int                do_sub;       /* TRUE to use HMM to infer start and end point of each seq
				    * and build a separate sub CM for alignment of that seq */
  int                do_fullsub;   /* TRUE to only remove structure outside HMM predicted start
				    * (spos) and end points (epos) */
  CM_t              *sub_cm;       /* sub covariance model                      */
  CM_t              *orig_cm;      /* the original, template covariance model the sub CM was built from */
  int                spos;         /* HMM node most likely to have emitted posn 1 of target seq */
  int                spos_state;   /* HMM state type for curr spos 0=match or 1=insert */
  int                epos;         /* HMM node most likely to have emitted posn L of target seq */
  int                epos_state;   /* HMM state type for curr epos 0=match or  1=insert */

  int              **orig2sub_smap;/* 2D state map from orig_cm (template) to sub_cm.
				    * 1st dimension - state index in orig_cm 
				    * 2nd D - 2 elements for up to 2 matching sub_cm states */
  int              **sub2orig_smap;/* 2D state map from orig_cm (template) to sub_cm.
				    * 1st dimension - state index in sub_cm (0..sub_cm->M-1)
				    * 2nd D - 2 elements for up to 2 matching orig_cm states */
  char            *check_outfile;  /* output file name for subCM stk file*/
  FILE            *check_ofp;      /* an open output file */
  Parsetree_t     *orig_tr;        /* parsetree for the orig_cm; created from the sub_cm parsetree */
  
  int *sub_node_cc_left; 
  int *sub_node_cc_right;
  int *sub_cc_node_map;
  int *orig_node_cc_left; 
  int *orig_node_cc_right;
  int *orig_cc_node_map;
  int **sub_cs2hn_map;  
  int **sub_cs2hs_map;  
  int ***sub_hns2cs_map;
  int **orig_cs2hn_map;  
  struct cplan9_s       *sub_hmm;       /* constructed CP9 HMM; written to hmmfile              */
  struct cplan9_s       *orig_hmm;      /* original CP9 HMM built from orig_cm */
  int                    sub_cp9_M;         /* number of nodes in CP9 HMM (MATL+MATR+2*MATP)        */
  int                    orig_hmm_M;        /* number of nodes in the original, template HMM */
  double **orig_phi;
  double **sub_phi_trad;    
  double **sub_phi_trunc;    
  int *imp_cc;                  /* imp_cc[k] = 1 if CP9 node k is an impossible case to get 
				 * the right transition distros for the sub_cm. */

  /*********************************************** 
   * Parse command line
   ***********************************************/
  format      = SQFILE_UNKNOWN;
  be_quiet    = FALSE;
  do_local    = FALSE;
  do_small    = TRUE;
  outfile     = NULL;
  regressfile = NULL;
  tracefile   = NULL;
  windowlen   = 200;
  set_window  = FALSE;
  do_qdb = FALSE;
  qdb_beta     = 0.0000001;
  do_full     = FALSE;
  do_expand = FALSE;
  bdump_level = 0;
  debug_level = 0;
  do_hbanded  = FALSE;
  hbandp      = 0.0001;
  use_sums    = FALSE;
  time_flag   = TRUE;
  do_inside   = FALSE;
  do_outside  = FALSE;
  do_check    = FALSE;
  do_post     = FALSE;
  do_sub      = FALSE;
  do_fullsub  = FALSE;
  do_checkcp9 = FALSE;
  check_outfile = "check.stk";
  seed         = time ((time_t *) NULL);

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-l")          == 0) do_local    = TRUE;
    else if (strcmp(optname, "-o")          == 0) outfile     = optarg;
    else if (strcmp(optname, "-q")          == 0) be_quiet    = TRUE;
    else if (strcmp(optname, "--nosmall")   == 0) do_small    = FALSE;
    else if (strcmp(optname, "--regress")   == 0) regressfile = optarg;
    else if (strcmp(optname, "--qdb")       == 0) do_qdb       = TRUE;
    else if (strcmp(optname, "--beta")      == 0) qdb_beta     = atof(optarg);
    else if (strcmp(optname, "-W")          == 0) {
      windowlen    = atoi(optarg); 
      set_window = TRUE; } 
    else if (strcmp(optname, "--full")      == 0) do_full      = TRUE;
    else if (strcmp(optname, "--expand")  == 0) do_expand  = TRUE;
    else if (strcmp(optname, "--banddump")  == 0) bdump_level  = atoi(optarg);
    else if (strcmp(optname, "--tfile")     == 0) tracefile    = optarg;
    else if (strcmp(optname, "--dlev")      == 0) debug_level  = atoi(optarg);
    else if (strcmp(optname, "--time")      == 0) time_flag    = TRUE;
    else if (strcmp(optname, "--inside")    == 0) do_inside    = TRUE;
    else if (strcmp(optname, "--outside")   == 0) { do_outside = TRUE; do_check = TRUE; }
    else if (strcmp(optname, "--post")      == 0) do_post      = TRUE;
    else if (strcmp(optname, "--checkpost") == 0) do_check     = TRUE;
    else if (strcmp(optname, "--sub")       == 0) do_sub       = TRUE; 
    else if (strcmp(optname, "--fsub")      == 0) { do_sub = TRUE; do_fullsub = TRUE; }
    else if (strcmp(optname, "--hbanded")   == 0) { do_hbanded = TRUE; do_small = FALSE; }
    else if (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--checkcp9")  == 0) do_checkcp9  = TRUE;
    else if (strcmp(optname, "--informat")  == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
    }
    else if (strcmp(optname, "-h") == 0) {
      MainBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

  if(do_inside && do_outside)
    Die("Please pick either --inside or --outside (--outside will run Inside()\nalso and check to make sure Inside() and Outside() scores agree).\n");
  if(do_checkcp9 && do_hbanded == FALSE)
    Die("--checkcp9 only makes sense with --hbanded\n");
  if(do_sub && do_local)
    Die("--sub and --local combination not yet supported.\n");
  if(do_sub && do_qdb)
    Die("Please pick either --sub or --qdb.\n");

  /* currently not set up for local alignment and posterior decode or outside run */
  /*  if(do_local && (do_outside || do_post))
    {
      Die("Can't do -l and either --post or --outside.\n");
    }
  */

  if (bdump_level > 3) Die("Highest available --banddump verbosity level is 3\n%s", usage);
  if (do_expand && (!(do_qdb))) Die("Doesn't make sense to use --expand option without --qdb option\n", usage);
  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  seqfile = argv[optind++]; 


  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (format == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    format = SQFILE_FASTA;
  
  /*****************************************************************
   * Input and configure the CM
   *****************************************************************/

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);
  orig_cm = cm;
  
  if (! (set_window)) windowlen = cm->W;
  /* Now that we know what windowlen is, we need to ensure that
   * cm->el_selfsc * W >= IMPOSSIBLE (cm->el_selfsc is the score for an EL self transition)
   * This is done because we potentially multiply cm->el_selfsc * W, and add
   * that to IMPOSSIBLE. To avoid underflow issues this value must be less than
   * 3 * IMPOSSIBLE. Here we guarantee its less than 2 * IMPOSSIBLE (to be safe).
   */
  if((cm->el_selfsc * windowlen) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (windowlen+1));

  if (do_local && do_hbanded)
    {
      printf("Warning: banding with an HMM (--hbanded) and allowing\nlocal alignment (-l). There's no telling what will happen.\n");
    } 

  CMLogoddsify(cm);
  /*CMHackInsertScores(cm);*/	/* "TEMPORARY" fix for bad priors */

  /*****************************************************************
   * Optionally, input and configure a Plan 7 or CM Plan 9 HMM
   * for banded alignment.
   *****************************************************************/
  if(do_hbanded || do_sub)
    {
      Alphabet_type = hmmNOTSETYET;
      SetAlphabet(hmmNUCLEIC); /* Set up the hmmer_alphabet global variable */
    }

  /*****************************************************************
   * Input and digitize the unaligned sequences
   *****************************************************************/

  if (! ReadMultipleRseqs(seqfile, format, &rseq, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);
  dsq = MallocOrDie(sizeof(char *) * nseq);
  for (i = 0; i < nseq; i++) 
    dsq[i] = DigitizeSequence(rseq[i], sqinfo[i].len);
  if(do_hbanded || do_sub)
    {
      p7dsq = MallocOrDie(sizeof(unsigned char *) * nseq);
      for(i = 0; i < nseq; i++)
	p7dsq[i] = hmmer_DigitizeSequence(rseq[i], sqinfo[i].len);
    }      
  /*********************************************** 
   * Show the banner
   ***********************************************/

  if (! be_quiet) 
    {
      MainBanner(stdout, banner);
      printf(   "CM file:              %s\n", cmfile);
      printf(   "Sequence file:        %s\n", seqfile);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /*****************************************************************
   * Do the work: 
   *        collect parse trees for each sequence
   *        run them thru Parsetrees2Alignment().
   *****************************************************************/
  
  tr    = MallocOrDie(sizeof(Parsetree_t) * nseq);
  minsc = FLT_MAX;
  maxsc = -FLT_MAX;
  avgsc = 0;

  if(do_post)
    {
      postcode = malloc(sizeof(char *) * nseq);
    }      

  watch1 = StopwatchCreate(); /*watch1 is used to time each step individually*/
  watch2 = StopwatchCreate(); /*watch2 times the full alignment (including band calc)
				for each seq*/

  if(do_hbanded || do_sub) /* We need a CP9 HMM to build sub_cms */
    {
      cp9_M = 0;
      for(v = 0; v <= cm->M; v++)
	{
	  if(cm->stid[v] ==  MATP_MP)
	    cp9_M += 2;
	  else if(cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
	    cp9_M ++;
	}
      /* build the HMM data structure */
      hmm = AllocCPlan9(cp9_M);
      ZeroCPlan9(hmm);

      /* Get information mapping the HMM to the CM and vice versa, used
       * for mapping bands. */
      map_consensus_columns(cm, hmm->M, &node_cc_left, &node_cc_right,
			    &cc_node_map, debug_level);
      
      CP9_map_cm2hmm_and_hmm2cm(cm, hmm, node_cc_left, node_cc_right, 
				cc_node_map, &cs2hn_map, &cs2hs_map, 
				&hns2cs_map, debug_level);
      
      /* fill in parameters of HMM using the CM and some ideas/formulas/tricks
       * from Zasha Weinberg's thesis (~p.123) (see CP9_cm2wrhmm.c code) */
      if(!(CP9_cm2wrhmm(cm, hmm, node_cc_left, node_cc_right, cc_node_map, cs2hn_map,
			cs2hs_map, hns2cs_map, debug_level)))
	Die("Couldn't build a CM Plan 9 HMM from the CM.\n");
      if(do_checkcp9)
	{
	  sre_srandom(seed);
	  //if(!(CP9_check_wrhmm_by_sampling(cm, hmm, 1, hmm->M, hns2cs_map, 0.05, 100000)))
	  if(!(CP9_check_wrhmm_by_sampling(cm, hmm, 1, hmm->M, hns2cs_map, 0.01, 100000, NULL,
					   NULL, NULL, debug_level)))
	    Die("CM Plan 9 fails sampling check!\n");
	  else
	    printf("CM Plan 9 passed sampling check.\n");
	}

      cp9_mx  = CreateCPlan9Matrix(1, hmm->M, 25, 0);
      hmm_M = hmm->M;

      /* Keep this data for the original CM safe; we'll be doing
       * pointer swapping to ease the sub_cm alignment implementation. */
      orig_hmm = hmm;
      orig_hmm_M = hmm_M;
      orig_node_cc_left =  node_cc_left;
      orig_node_cc_right = node_cc_right;
      orig_cc_node_map   = cc_node_map;
      orig_cs2hn_map     = cs2hn_map;

      StopwatchZero(watch2);
      StopwatchStart(watch2);
    }

  /* Relocated ConfigLocal() call to here, below the CM Plan 9 construction.
   * Otherwise its impossible to make a CM Plan 9 HMM from the local CM
   * that passes the current tests to ensure the HMM is "close enough" to
   * the CM. This is something to look into later.
   */
  if (do_local)
    { 
      ConfigLocal(cm, 0.5, 0.5);
      CMLogoddsify(cm);
      /*CMHackInsertScores(cm);*/	/* "TEMPORARY" fix for bad priors */
    }

  if((do_local && do_hbanded) || do_sub) /* to get spos and epos for the sub_cm, 
					  * we config the HMM to local mode.*/
      {
	printf("configuring the CM plan 9 HMM for local alignment.\n");
	swentry           = 0.5;
	swexit            = 0.5;
	CPlan9SWConfig(hmm, swentry, swexit);
	CP9Logoddsify(hmm);

	if(do_sub)
	  orig_tr    = MallocOrDie(sizeof(Parsetree_t));
      }

  /* set up the query dependent bands, this has to be done after the ConfigLocal() call */
  if(do_qdb || bdump_level > 0)
    {
      safe_windowlen = windowlen * 2;
      while(!(BandCalculationEngine(cm, safe_windowlen, qdb_beta, 0, &dmin, &dmax, &gamma, do_local)))
	{
	  FreeBandDensities(cm, gamma);
	  free(dmin);
	  free(dmax);
	  safe_windowlen *= 2;
	  printf("ERROR BandCalculationEngine returned false, windowlen adjusted to %d\n", safe_windowlen);
	}

      /* EPN 11.13.05 
       * An important design decision.
       * We're changing the windowlen value here. By default,
       * windowlen is read from the cm file (set to cm->W). 
       * Here we're doing a banded alignment though. Its pointless to allow
       * a windowlen that's greater than the largest possible banded hit 
       * (which is dmax[0]). So we reset windowlen to dmax[0].
       * Its also possible that BandCalculationEngine() returns a dmax[0] that 
       * is > cm->W. This should only happen if the qdb_beta we're using now is < 1E-7 
       * (1E-7 is the qdb_beta value used to determine cm->W in cmbuild). If this 
       * happens, the current implementation reassigns windowlen to this larger value.
       * NOTE: if W was set at the command line, the command line value is 
       *       always used.
       */
      if(!(set_window))
	windowlen = dmax[0];
      if(bdump_level > 1) 
	{
	  printf("qdb_beta:%f\n", qdb_beta);
	  debug_print_bands(cm, dmin, dmax);
	  PrintDPCellsSaved(cm, dmin, dmax, windowlen);
	}
      expand_flag = FALSE;
    }

/* Allocate data structures for use with HMM banding strategy */
  if(do_hbanded)
    {
      pn_min_m    = malloc(sizeof(int) * (hmm_M+1));
      pn_max_m    = malloc(sizeof(int) * (hmm_M+1));
      pn_min_i    = malloc(sizeof(int) * (hmm_M+1));
      pn_max_i    = malloc(sizeof(int) * (hmm_M+1));
      pn_min_d    = malloc(sizeof(int) * (hmm_M+1));
      pn_max_d    = malloc(sizeof(int) * (hmm_M+1));
      imin        = malloc(sizeof(int) * cm->M);
      imax        = malloc(sizeof(int) * cm->M);
      jmin        = malloc(sizeof(int) * cm->M);
      jmax        = malloc(sizeof(int) * cm->M);
      hdmin       = malloc(sizeof(int *) * cm->M);
      hdmax       = malloc(sizeof(int *) * cm->M);
      safe_hdmin  = malloc(sizeof(int) * cm->M);
      safe_hdmax  = malloc(sizeof(int) * cm->M);
      isum_pn_m   = malloc(sizeof(int *) * nseq);
      isum_pn_i   = malloc(sizeof(int *) * nseq);
      isum_pn_d   = malloc(sizeof(int *) * nseq);
      
      for (i = 0; i < nseq; i++)
	{
	  isum_pn_m[i] = malloc(sizeof(int) * (hmm_M+1));
	  isum_pn_i[i] = malloc(sizeof(int) * (hmm_M+1));
	  isum_pn_d[i] = malloc(sizeof(int) * (hmm_M+1));
	}
    }

  for (i = 0; i < nseq; i++)
    {
      StopwatchZero(watch2);
      StopwatchStart(watch2);

      /* Potentially, do HMM calculations. */
      if(do_hbanded || do_sub)
	{
	  /* We want HMM posteriors for this sequence to the full length (non-sub) HMM */
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);

	  /* Step 1: Get HMM posteriors.*/
	  /*sc = CP9Viterbi(p7dsq[i], 1, sqinfo[i].len, hmm, cp9_mx);*/
	  forward_sc = CP9Forward(p7dsq[i], 1, sqinfo[i].len, orig_hmm, &cp9_fwd);
	  printf("CP9 i: %d | forward_sc : %f\n", i, forward_sc);
	  backward_sc = CP9Backward(p7dsq[i], 1, sqinfo[i].len, orig_hmm, &cp9_bck);
	  printf("CP9 i: %d | backward_sc: %f\n", i, backward_sc);

	  /*debug_check_CP9_FB(cp9_fwd, cp9_bck, hmm, forward_sc, 1, sqinfo[i].len, p7dsq[i]);*/
	  cp9_posterior = cp9_bck;
	  CP9FullPosterior(p7dsq[i], 1, sqinfo[i].len, orig_hmm, cp9_fwd, cp9_bck, cp9_posterior);
	}
      /* If we're in sub mode:
       * (1) Get HMM posteriors. (we already did this above)
       * (2) Infer the start (spos) and end (epos) HMM states by 
       *     looking at the posterior matrix.
       * (3) Build the sub_cm from the original CM.
       *
       * If we're also doing HMM banded alignment:
       * (4) Build a new CP9 HMM from the sub CM.
       * (4) Build a new CP9 HMM from the sub CM.
       * (4) Build a new CP9 HMM from the sub CM.
       * (5) Do Forward/Backward again, and get a new posterior matrix.
       */
      if(do_sub)
	{
	  /* (2) infer the start and end HMM states by looking at the posterior matrix.
	   * Remember: we're necessarily in local mode, the --sub option turns local mode on. 
	   */
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, 1,             cp9_posterior, &spos, &spos_state);
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, sqinfo[i].len, cp9_posterior, &epos, &epos_state);
	  
	  if(spos == 0 && spos_state == 1) 
	    {
	      /* The most likely state to have emitted first residue
	       * is the insert state in node 0, in this case the first consensus 
	       * column we want to model is 1, not 0.
	       */
	      spos = 1;
	    }
	  if(epos == 0 && epos_state == 1) 
	    {
	      /* The most likely state to have emitted first residue
	       * is the insert state in node 0, in this case the first consensus 
	       * column we want to model is 1, not 0.
	       */
	      epos = 1;
	    }
	  if(epos < spos) /* This can happen, if it does our alignment will be crap, still we 
			     shouldn't crash. */
	    epos = spos;
	  

	  /* (3) Build the sub_cm from the original CM. */
	  build_sub_cm(orig_cm, &sub_cm, spos, epos, &orig2sub_smap, &sub2orig_smap, 
		       NULL,                   /* Don't return imp_cc, I'm not checking the sub_cm */
		       NULL, NULL, NULL, NULL, /* Don't return *predict* info, I'm not checking the sub_cm */
		       0.00001,                /* threshold for tests irrelevant b/c next is FALSE */
		       do_fullsub,             /* Says to build or not build a sub CM that models all columns */
		       FALSE, FALSE,           /* don't do analytical or sampling check */
		       0.01, 10000,            /* chi-square threshold and num samples irrelevant b/c prev is FALSE */
		       debug_level);           /* TRUE to print debugging info */
	  
	  cm    = sub_cm; /* orig_cm still points to the original CM */

	  if(do_hbanded) /* we're doing HMM banded alignment to the sub_cm */
	    {
	      /* (4) Build a new CP9 HMM from the sub CM. */
	      /* Eventually, I think we can do this by just adjusting the parameters of the original HMM 
		 CP9_2sub_cp9(hmm, &sub_hmm2, spos, epos, orig_phi);
	      */
	      sub_hmm = AllocCPlan9(epos-spos+1);
	      ZeroCPlan9(sub_hmm);
	      /* Get information mapping the HMM to the CM and vice versa, used
	       * for mapping bands. 
	       */

	      map_consensus_columns(sub_cm, sub_hmm->M, &sub_node_cc_left, &sub_node_cc_right,
				    &sub_cc_node_map, debug_level);
	      
	      CP9_map_cm2hmm_and_hmm2cm(sub_cm, sub_hmm, sub_node_cc_left, sub_node_cc_right,
					sub_cc_node_map, &sub_cs2hn_map, &sub_cs2hs_map, 
					&sub_hns2cs_map, debug_level);
	      if(!(CP9_cm2wrhmm(sub_cm, sub_hmm, sub_node_cc_left, sub_node_cc_right, 
				sub_cc_node_map, sub_cs2hn_map,
				sub_cs2hs_map, sub_hns2cs_map, debug_level)))
		Die("Couldn't build a CM Plan 9 HMM from the sub CM.\n");
	      
	      /* (5) Do Forward/Backward again, and get a new posterior matrix. 
	       * We have to free cp9_fwd and cp9_posterior because we used them 
	       * to find spos and epos. */

	      FreeCPlan9Matrix(cp9_fwd);
	      FreeCPlan9Matrix(cp9_posterior);
	      forward_sc = CP9Forward(p7dsq[i], 1, sqinfo[i].len, sub_hmm, &cp9_fwd);
	      printf("CP9 i: %d | forward_sc : %f\n", i, forward_sc);
	      backward_sc = CP9Backward(p7dsq[i], 1, sqinfo[i].len, sub_hmm, &cp9_bck);
	      printf("CP9 i: %d | backward_sc: %f\n", i, backward_sc);
	      /*debug_check_CP9_FB(cp9_fwd, cp9_bck, hmm, forward_sc, 1, sqinfo[i].len, p7dsq[i]);*/
	      cp9_posterior = cp9_bck;
	      CP9FullPosterior(p7dsq[i], 1, sqinfo[i].len, sub_hmm, cp9_fwd, cp9_bck, cp9_posterior);
	      /* cp9_posterior has the posteriors for the sub_hmm */
	      
	      /* Change some pointers so that the functions that create bands use the
	       * sub_* data structures. The orig_* data structures will still point
	       * to the original CM versions. */
	      hmm_M         = sub_hmm->M;
	      hmm           = sub_hmm;    
	      node_cc_left  = sub_node_cc_left;
	      node_cc_right = sub_node_cc_right;
	      cc_node_map   = sub_cc_node_map;  
	      cs2hn_map     = sub_cs2hn_map;
	    }
	}
      if(do_hbanded)
	{
	  StopwatchStop(watch1);
	  if(time_flag) StopwatchDisplay(stdout, "CP9 Forward/Backward CPU time: ", watch1);
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);
      
	  /* Align the current seq to the cp9 HMM, we don't care
	   * about the trace, just the posteriors.
	   * Step 1: Get HMM posteriors. (we already did this above, if do_sub,
	   *                              the posteriors are for the sub_hmm)
	   * Step 2: posteriors -> HMM bands.
	   * Step 3: HMM bands  ->  CM bands.
	   */
	  
	  /* Step 2: posteriors -> HMM bands.*/
	  if(!(use_sums))
	    {
	      /* match states */
	      CP9_hmm_band_bounds(cp9_posterior->mmx, 1, sqinfo[i].len, hmm_M,
				  NULL, pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
	      /* insert states */
	      CP9_hmm_band_bounds(cp9_posterior->imx, 1, sqinfo[i].len, hmm_M,
				  NULL, pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
	      /* delete states (note: delete_flag set to TRUE) */
	      CP9_hmm_band_bounds(cp9_posterior->dmx, 1, sqinfo[i].len, hmm_M,
				  NULL, pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);
	    }
	  else
	    {
	      CP9_ifill_post_sums(cp9_posterior, 1, sqinfo[i].len, hmm_M,
				  isum_pn_m[i], isum_pn_i[i], 
				  isum_pn_d[i]);
	      /* match states */
	      CP9_hmm_band_bounds(cp9_posterior->mmx, 1, sqinfo[i].len, hmm_M,
				  isum_pn_m[i], pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
	      /* insert states */
	      CP9_hmm_band_bounds(cp9_posterior->imx, 1, sqinfo[i].len, hmm_M,
				  isum_pn_i[i], pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
	      /* delete states */
	      CP9_hmm_band_bounds(cp9_posterior->dmx, 1, sqinfo[i].len, hmm_M,
				  isum_pn_d[i], pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);
	    }
	  if(debug_level != 0)
	    {
	      printf("printing hmm bands\n");
	      print_hmm_bands(stdout, sqinfo[i].len, hmm_M, pn_min_m, pn_max_m, pn_min_i,
			      pn_max_i, pn_min_d, pn_max_d, hbandp, debug_level);
	    }
	  
	  /* Step 3: HMM bands  ->  CM bands. */
	  printf("10.24.06 cm->nodes: %d\n", cm->nodes);
	  printf("10.24.06 hmm_M    : %d\n", hmm_M);
	  hmm2ij_bands(cm, hmm_M, node_cc_left, node_cc_right, cc_node_map, 
		       1, sqinfo[i].len, pn_min_m, pn_max_m, pn_min_i, pn_max_i, 
		       pn_min_d, pn_max_d, imin, imax, jmin, jmax, cs2hn_map,
		       debug_level);
	  
	  StopwatchStop(watch1);
	  if(time_flag) StopwatchDisplay(stdout, "CP9 Band calculation CPU time: ", watch1);
	  /* Use the CM bands on i and j to get bands on d, specific to j. */
	  for(v = 0; v < cm->M; v++)
	    {
	      hdmin[v] = malloc(sizeof(int) * (jmax[v] - jmin[v] + 1));
	      hdmax[v] = malloc(sizeof(int) * (jmax[v] - jmin[v] + 1));
	    }
	  ij2d_bands(cm, sqinfo[i].len, imin, imax, jmin, jmax,
		     hdmin, hdmax, -1);
	  
	  if(debug_level != 0)
	    PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, sqinfo[i].len);
	  
	  FreeCPlan9Matrix(cp9_fwd);
	  FreeCPlan9Matrix(cp9_posterior);
	  /* Done with the HMM. On to the CM. */
	}
      
      /* Determine which CYK alignment algorithm to use, based
       * on command-line options AND memory requirements.
       */
      if(do_hbanded)
	{
	  /* write a function to determine size of jd banded memory
	   * req'd, and set do_small to true if its > thresh.
	   if(do_small) * We're only going to band on d in memory, but 
	   * we need to calculate safe_hd bands on the d dimension. 
	   {
	  */
	}
      
      if(do_expand)
	{
	  /*First, check to see if we need to reset the apriori bands b/c 
	   * they're currently expanded. */
	  if(expand_flag)
	    {
	      FreeBandDensities(cm, gamma);
	      free(dmin);
	      free(dmax);
	      while(!(BandCalculationEngine(cm, safe_windowlen, qdb_beta, 0, &dmin, &dmax, &gamma, do_local)))
		{
		  FreeBandDensities(cm, gamma);
		  free(dmin);
		  free(dmax);
		  safe_windowlen *= 2;
		}
	      expand_flag = FALSE;
	    }
	  if((sqinfo[i].len < dmin[0]) || (sqinfo[i].len > dmax[0]))
	    {
	      /* the seq we're aligning is outside the root band, so we expand.*/
	      ExpandBands(cm, sqinfo[i].len, dmin, dmax);
	      printf("Expanded bands for seq : %s\n", sqinfo[i].name);
	      if(bdump_level > 2) 
		{
		  printf("printing expanded bands :\n");
		  debug_print_bands(cm, dmin, dmax);
		}
	      expand_flag = TRUE;
	    }
	}
      else 
	{
	  if (do_qdb && (sqinfo[i].len < dmin[0] || sqinfo[i].len > dmax[0]))
	    {
	      /* the seq we're aligning is outside the root band, but
	       * --expand was not enabled, so we die.*/
	      Die("Length of sequence to align (%d nt) lies outside the root band.\ndmin[0]: %d and dmax[0]: %d\nImpossible to align with query dependent banded CYK unless you try --expand.\n%s", sqinfo[i].len, dmin[0], dmax[0], usage);
	    }
	}
      printf("Aligning %s\n", sqinfo[i].name);
      if (do_inside)
	{
	  if(do_hbanded)
	    sc = FInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				 BE_PARANOID,	/* non memory-saving mode */
				 NULL, NULL,	/* manage your own matrix, I don't want it */
				 NULL, NULL,	/* manage your own deckpool, I don't want it */
				 do_local,        /* TRUE to allow local begins */
				 jmin, jmax, hdmin, hdmax); /* j and d bands */
	  else
	    sc = FInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			 BE_EFFICIENT,	/* memory-saving mode */
			 NULL, NULL,	/* manage your own matrix, I don't want it */
			 NULL, NULL,	/* manage your own deckpool, I don't want it */
			 do_local);       /* TRUE to allow local begins */
	}
      else if(do_outside)
	{	
	  if(do_hbanded)
	    {
	      sc = FInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,        /* TRUE to allow local begins */
				   jmin, jmax, hdmin, hdmax); /* j and d bands */
	      /*do_check = TRUE;*/
	      sc = FOutside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				    BE_PARANOID,	/* save full beta */
				    NULL, NULL,	/* manage your own matrix, I don't want it */
				    NULL, NULL,	/* manage your own deckpool, I don't want it */
				    do_local,       /* TRUE to allow local begins */
				    alpha,          /* alpha matrix from FInside_b_jd_me() */
				    NULL,           /* don't save alpha */
				    do_check,       /* TRUE to check Outside probs agree with Inside */
				    jmin, jmax, hdmin, hdmax); /* j and d bands */
	    }
	  else
	    {
	      sc = FInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = FOutside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			    BE_PARANOID,	/* save full beta */
			    NULL, NULL,	/* manage your own matrix, I don't want it */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha,          /* alpha matrix from FInside() */
			    NULL,           /* don't save alpha */
			    do_check);      /* TRUE to check Outside probs agree with Inside */
	    }
	}
      else if (do_small) 
	{
	  if(do_qdb)
	    {
	      sc = CYKDivideAndConquer_b(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, 
					      &(tr[i]), dmin, dmax);
	      if(bdump_level > 0)
		banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	  else if(do_hbanded) /*j and d bands not tight enough to allow HMM banded full CYK*/
	    {
	      /* Calc the safe d bands */
	      hd2safe_hd_bands(cm->M, jmin, jmax, hdmin, hdmax, safe_hdmin, safe_hdmax);
	      if(debug_level > 3)
		{
		  printf("\nprinting hd bands\n\n");
		  debug_print_hd_bands(cm, hdmin, hdmax, jmin, jmax);
		  printf("\ndone printing hd bands\n\n");
		}
	      /* Note the following CYK call will not enforce j bands, even
	       * though user specified --hbanded. */
	      sc = CYKDivideAndConquer_b(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, 
					      &(tr[i]), safe_hdmin, safe_hdmax);
	      if(bdump_level > 0)
		banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	  else
	    {
	      sc = CYKDivideAndConquer(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]));
	      if(bdump_level > 0)
		{
		  /* We want band info but --banded wasn't used.  Useful if you're curious
		   * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
		   * to see where the non-banded parse went outside the bands.
		   */
		  banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
		}
	    }
	}
      else if(do_qdb)
	{
	  sc = CYKInside(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), dmin, dmax);
	  if(bdump_level > 0)
	    banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	}
      else if(do_hbanded)
	{
	  sc = CYKInside_b_jd(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), jmin, jmax,
	    hdmin, hdmax, safe_hdmin, safe_hdmax);
	  if(bdump_level > 0)
	    banded_trace_info_dump(cm, tr[i], safe_hdmin, safe_hdmax, bdump_level);
	}
      else
	{
	  sc = CYKInside(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), NULL, NULL);
	  if(bdump_level > 0)
	    {
	      /* We want band info but --banded wasn't used.  Useful if you're curious
	       * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
	       * to see where the non-banded parse went outside the bands.
	       */
	      banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	}
      if(do_post) /* Do Inside() and Outside() runs and use alpha and beta to get posteriors */
	{	
	  /*alpha = MallocOrDie(sizeof(float **) * (cm->M));
	  beta  = MallocOrDie(sizeof(float **) * (cm->M+1));
	  */
	  post  = MallocOrDie(sizeof(float **) * (cm->M+1));
	  /*
	  for (v = 0; v < cm->M; v++) alpha[v] = NULL;
	  for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
	  */
	  if(do_hbanded)
	    {
	      for (v = 0; v < cm->M; v++)
		{
		  post[v] = NULL;
		  post[v] = alloc_jdbanded_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len, jmin[v], jmax[v], hdmin[v], hdmax[v]);
		}
	      post[cm->M] = NULL;
	      post[cm->M] = alloc_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len);
	      sc = FInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,       /* TRUE to allow local begins */
				   jmin, jmax, hdmin, hdmax); /* j and d bands */
	      sc = FOutside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				    BE_PARANOID,	/* save full beta */
				    NULL, &beta,	/* fill beta, and return it, needed for CMPosterior() */
				    NULL, NULL,	/* manage your own deckpool, I don't want it */
				    do_local,       /* TRUE to allow local begins */
				    alpha, &alpha,  /* alpha matrix from FInside(), and save it for CMPosterior*/
				    do_check,      /* TRUE to check Outside probs agree with Inside */
				    jmin, jmax, hdmin, hdmax); /* j and d bands */
	      CMPosterior_b_jd_me(sqinfo[i].len, cm, alpha, NULL, beta, NULL, post, &post,
				    jmin, jmax, hdmin, hdmax);
	      postcode[i] = CMPostalCode_b_jd_me(cm, sqinfo[i].len, post, tr[i],
						 jmin, jmax, hdmin, hdmax);
	    }
	  else
	    {
	      for (v = 0; v < cm->M+1; v++)
		{
		  post[v] = NULL;
		  post[v] = alloc_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len);
		}
	      sc = FInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = FOutside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			    BE_PARANOID,	/* save full beta */
			    NULL, &beta,	/* fill beta, and return it, needed for CMPosterior() */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha, &alpha,  /* alpha matrix from FInside(), and save it for CMPosterior*/
			    do_check);      /* TRUE to check Outside probs agree with Inside */
	      CMPosterior(sqinfo[i].len, cm, alpha, NULL, beta, NULL, post, &post);
	      if(do_check)
		{
		  CMCheckPosterior(sqinfo[i].len, cm, post);
		  printf("\nPosteriors checked.\n\n");
		}
	      postcode[i] = CMPostalCode(cm, sqinfo[i].len, post, tr[i]);
	    }

	  /* free post */
	  if(post != NULL)
	    {
	      for (v = 0; v <= (cm->M); v++)
		if (post[v] != NULL) { free_vjd_deck(post[v], 1, sqinfo[i].len); post[v] = NULL;}
	      free(post);
	    }
	}

      avgsc += sc;
      if (sc > maxsc) maxsc = sc;
      if (sc < minsc) minsc = sc;

      printf("Alignment score for %30s: %8.2f bits\n", sqinfo[i].name, sc);

      /* If debug level high enough, print out the parse tree */
      if(debug_level > 2)
	{
	  fprintf(stdout, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], dsq[i], FALSE));;
	  ParsetreeDump(stdout, tr[i], cm, dsq[i]);
	  fprintf(stdout, "//\n");
	}
      /* Dump the trace with info on i, j and d bands
       * if bdump_level is high enough */
      if(bdump_level > 0)
	ijd_banded_trace_info_dump(cm, tr[i], imin, imax, jmin, jmax, hdmin, hdmax, 1);
      
      /* Clean up the structures we use calculating HMM bands, that are allocated
       * differently for each sequence. 
       */
      if(do_hbanded)
	{
	  for(v = 0; v < cm->M; v++)
	    { 
	      free(hdmin[v]); 
	      free(hdmax[v]);
	    }
	  StopwatchStop(watch2);
	  if(time_flag) StopwatchDisplay(stdout, "band calc and jd CYK CPU time: ", watch2);
	}
      if(do_sub)
	{
	  /* Convert the sub_cm parsetree to a full CM parsetree */
	  if(debug_level > 0)
	    ParsetreeDump(stdout, tr[i], cm, dsq[i]);
	  if(!(sub_cm2cm_parsetree(orig_cm, sub_cm, &orig_tr, tr[i], spos, epos, orig2sub_smap, sub2orig_smap, FALSE, debug_level)))
	    {
	      printf("\n\nIncorrectly converted original trace:\n");
	      ParsetreeDump(stdout, orig_tr, orig_cm, dsq[i]);
	      exit(1);
	    }
	  if(debug_level > 0)
	    {
	      printf("\n\nConverted original trace:\n");
	      ParsetreeDump(stdout, orig_tr, orig_cm, dsq[i]);
	    }
	  /* Replace the sub_cm trace with the converted orig_cm trace. */
	  FreeParsetree(tr[i]);
	  tr[i] = orig_tr;

	  FreeCM(sub_cm); /* cm and sub_cm now point to NULL */
	  if(do_hbanded)
	    FreeCPlan9(sub_hmm);
	}
    }
  avgsc /= nseq;

  if(!(do_inside || do_outside))
    {
      msa = Parsetrees2Alignment(orig_cm, dsq, sqinfo, NULL, tr, nseq, do_full);
      
      if(do_post)
	{
	  for (i = 0; i < nseq; i++)
	    {
	      MakeAlignedString(msa->aseq[i], msa->alen, postcode[i], &apostcode); 
	      MSAAppendGR(msa, "POST", i, apostcode);
	      free(apostcode);
	      free(postcode[i]);
	    }
	  free(postcode);
	}

      /*****************************************************************
       * Output the alignment.
       *****************************************************************/
      
      if (outfile != NULL && (ofp = fopen(outfile, "w")) != NULL) 
	{
	  WriteStockholm(ofp, msa);
	  printf("Alignment saved in file %s\n", outfile);
	  fclose(ofp);
	}
      else
	WriteStockholm(stdout, msa);
      
      /* Detailed traces for debugging training set.
       */
      if (tracefile != NULL)       
	{
	  printf("%-40s ... ", "Saving parsetrees"); fflush(stdout);
	  if ((ofp = fopen(tracefile,"w")) == NULL)
	    Die("failed to open trace file %s", tracefile);
	  for (i = 0; i < msa->nseq; i++) 
	    {
	      fprintf(ofp, "> %s\n", msa->sqname[i]);
	      fprintf(ofp, "  SCORE : %.2f bits\n", ParsetreeScore(orig_cm, tr[i], dsq[i], FALSE));;
	      ParsetreeDump(ofp, tr[i], orig_cm, dsq[i]);
	      fprintf(ofp, "//\n");
	    }
	  fclose(ofp);
	  printf("done. [%s]\n", tracefile);
	}


      if (regressfile != NULL && (ofp = fopen(regressfile, "w")) != NULL) 
	{
	  /* Must delete author info from msa, because it contains version
	   * and won't diff clean in regression tests.
	   */
	  free(msa->au); msa->au = NULL;
	  WriteStockholm(ofp, msa);
	  fclose(ofp);
	}
    }

  /*****************************************************************
   * Clean up and exit.
   *****************************************************************/
  
  for (i = 0; i < nseq; i++) 
    {
      if(!(do_inside || do_outside)) FreeParsetree(tr[i]);
      free(dsq[i]);
      FreeSequence(rseq[i], &(sqinfo[i]));
    }
  if (do_qdb)
    {
      FreeBandDensities(cm, gamma);
      free(dmin);
      free(dmax);
    }

  if(do_hbanded || do_sub)
    {
      for(v = 0; v <= orig_cm->M; v++)
	{
	  free(orig_cs2hn_map[v]);
	  free(cs2hs_map[v]);
	}
      free(cs2hn_map);
      free(cs2hs_map);
      for(k = 0; k <= hmm_M; k++)
	{
	  for(ks = 0; ks < 3; ks++)
	    free(hns2cs_map[k][ks]);
	  free(hns2cs_map[k]);
	}
      free(hns2cs_map);
      free(node_cc_left);
      free(node_cc_right);
      free(cc_node_map);
      for (i = 0; i < nseq; i++) 
	free(p7dsq[i]);
      free(p7dsq);
    }
  if(do_hbanded)
    {
      /* HMMERNAL */
      free(imin);
      free(imax);
      free(jmin);
      free(jmax);
      free(hdmin);
      free(hdmax);
      free(pn_min_m);
      free(pn_max_m);
      free(pn_min_i);
      free(pn_max_i);
      free(pn_min_d);
      free(pn_max_d);
      free(safe_hdmin);
      free(safe_hdmax);
      for (i = 0; i < nseq; i++) 
	{ 
	  free(isum_pn_m[i]); 
	  free(isum_pn_i[i]);
	  free(isum_pn_d[i]);
	  /* Uncomment if calc'ing a postcode */
	}
      free(isum_pn_m);
      free(isum_pn_i);
      free(isum_pn_d);
    }

  if(do_hbanded || do_sub)
    {
      FreeCPlan9Matrix(cp9_mx);
      FreeCPlan9(orig_hmm);
    }

  if(!(do_inside || do_outside)) MSAFree(msa);
  FreeCM(orig_cm);
  free(rseq);
  free(sqinfo);
  free(dsq);
  free(tr);
  
  SqdClean();
  StopwatchFree(watch1);
  StopwatchFree(watch2);
  return 0;
}




/* EPN 07.22.05
 * ExpandBands()
 * Function: ExpandBands
 *
 * Purpose:  Called when the sequence we are about to align 
 *           using bands is either shorter in length than
 *           the dmin on the root state, or longer in length
 *           than the dmax on the root state.
 *            
 *           This function expands the bands on ALL states
 *           v=1..cm->M-1 in the following manner :
 *           
 *           case 1 : target len < dmin[0]
 *                    subtract (dmin[0]-target len) from
 *                    dmin of all states, and ensure
 *                    dmin[v]>=0 for all v.
 *                    Further :
 *                    if cm->sttype[v] == MP_st ensure dmin[v]>=2;
 *                    if cm->sttype[v] == IL_st || ML_st ensure dmin[v]>=1;
 *                    if cm->sttype[v] == IR_st || MR_st ensure dmin[v]>=1;
 *                        
 *           case 2 : target len > dmax[0]
 *                    add (target len-dmax[0] to dmax
 *                    of all states.
 *
 *           Prior to handling such situtations with this
 *           hack, the program would choke and die.  This
 *           hacky approach is used as a simple, inefficient
 *           not well thought out, but effective way to 
 *           solve this problem.
 * 
 * Args:    cm       - the CM
 *          tlen     - length of target sequence about to be aligned
 *          dmin     - minimum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          dmax     - maximum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *
 * Returns: (void) 
 */

static void
ExpandBands(CM_t *cm, int tlen, int *dmin, int *dmax)
{
  int v;
  int diff;
  int root_min;
  int root_max;
  int M = cm->M;
  root_min = dmin[0];
  root_max = dmax[0];

  if(tlen < root_min)
    {
      diff = root_min - tlen;
      for(v=0; v<M; v++)
	{
	  dmin[v] -= diff;
	  if((cm->sttype[v] == MP_st) && (dmin[v] < 2)) 
	    dmin[v] = 2;
	  else if(((cm->sttype[v] == IL_st) || (cm->sttype[v] == ML_st)) 
		  && (dmin[v] < 1)) 
	    dmin[v] = 1;
	  else if(((cm->sttype[v] == IR_st) || (cm->sttype[v] == MR_st)) 
		  && (dmin[v] < 1)) 
	    dmin[v] = 1;
	  else if(dmin[v] < 0) 
	    dmin[v] = 0;
	}
    }
  else if(tlen > root_max)
    {
      diff = tlen - root_min;
      for(v=0; v<M; v++)
	{
	  dmax[v] += diff;
	}
    }
  printf("Expanded bands : \n");
}

/* EPN 08.15.05
 * banded_trace_info_dump()
 * Function: banded_trace_info_dump
 *
 * Purpose:  Called when the user has enabled the --banddump
 *           options.  This function determines how close the
 *           trace was to the bands at each state in the trace,
 *           and prints out that information in differing levels
 *           of verbosity depending on an input parameter 
 *           (bdump_level).
 * 
 * Args:    tr       - the parsetree (trace)
 *          dmin     - minimum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          dmax     - maximum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          bdump_level - level of verbosity
 * Returns: (void) 
 */

static void
banded_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, int *dmax, int bdump_level)
{
  char **sttypes;
  char **nodetypes;
  int v, i, j, d, tpos;
  int mindiff;            /* d - dmin[v] */
  int maxdiff;            /* dmax[v] - d */

  sttypes = malloc(sizeof(char *) * 10);
  sttypes[0] = "D";
  sttypes[1] = "MP";
  sttypes[2] = "ML";
  sttypes[3] = "MR";
  sttypes[4] = "IL";
  sttypes[5] = "IR";
  sttypes[6] = "S";
  sttypes[7] = "E";
  sttypes[8] = "B";
  sttypes[9] = "EL";

  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  for (tpos = 0; tpos < tr->n; tpos++)
    {
      v  = tr->state[tpos];
      i = tr->emitl[tpos];
      j = tr->emitr[tpos];
      d = j-i+1;

      if(cm->sttype[v] != EL_st)
	{
	  mindiff = d-dmin[v];
	  maxdiff = dmax[v]-d;
	  if(bdump_level > 1 || ((mindiff < 0) || (maxdiff < 0)))
	    printf("%-4s %-3s v: %4d | d: %4d | dmin: %4d | dmax: %4d | %3d | %3d |\n", nodetypes[(int) cm->ndtype[(int) cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], v, d, dmin[v], dmax[v], mindiff, maxdiff);
	}
      else
	{
	  if(bdump_level > 1)
	    printf("%-8s v: %4d | d: %4d |\n", sttypes[(int) cm->sttype[v]], v, d);
	}
    }
  free(sttypes);
  free(nodetypes);
}


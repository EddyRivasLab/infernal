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
   --banddump <n>: set verbosity of band info print statements to <n> [1..3]\n\
   --dlev <n>    : set verbosity of debugging print statements to <n> [1..3]\n\
   --time        : print timings for alignment, band calculation, etc.\n\
   --inside      : don't align; return scores from the Inside algorithm \n\
   --outside     : don't align; return scores from the Outside algorithm\n\
   --post        : align with CYK and append posterior probabilities\n\
   --checkpost   : check that posteriors are correctly calc'ed\n\
   --sub         : use HMM predicted start and end points to build a sub CM\n\
\n\
  * HMM banded alignment related options:\n\
   --hbanded     : use exptl HMM banded CYK aln algorithm (df: builds CP9 HMM) \n\
   --hbandp <f>  : tail loss prob for --hbanded (default:0.0001)\n\
   --cp9 <f>     : use the CM plan 9 HMM in file <f> for band calculation\n\
   --p7 <f>      : use the plan 7 HMMER 2.4 HMM in file <f> for band calculation\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
   --checkcp9    : check the CP9 empirically by generating sequences\n\
\n\
  * A priori banded alignment related options:\n\
   --apbanded    : use experimental a priori (ap) banded CYK alignment algorithm\n\
   --apbandp <f> : tail loss prob for --apbanded (default:0.0001)\n\
   --apexpand    : naively expand ap bands if target seq is outside root band\n\
    -W <n>       : window size for calc'ing ap bands (df: precalc'd in cmbuild)\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-l", TRUE, sqdARG_NONE },
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--nosmall",    FALSE, sqdARG_NONE },
  { "--regress",    FALSE, sqdARG_STRING },
  { "--apbanded",     FALSE, sqdARG_NONE },
  { "--apbandp",      FALSE, sqdARG_FLOAT},
  { "--apexpand", FALSE, sqdARG_NONE},
  { "--banddump"  , FALSE, sqdARG_INT},
  { "-W", TRUE, sqdARG_INT },
  { "--full", FALSE, sqdARG_NONE },
  { "--dlev",       FALSE, sqdARG_INT },
  { "--hbanded",    FALSE, sqdARG_NONE },
  { "--hbandp",     FALSE, sqdARG_FLOAT},
  { "--sums",       FALSE, sqdARG_NONE},
  { "--cp9",        FALSE, sqdARG_STRING},
  { "--p7",         FALSE, sqdARG_STRING},
  { "--time",       FALSE, sqdARG_NONE},
  { "--inside",     FALSE, sqdARG_NONE},
  { "--outside",    FALSE, sqdARG_NONE},
  { "--post",       FALSE, sqdARG_NONE},
  { "--checkpost",  FALSE, sqdARG_NONE},
  { "--sub",        FALSE, sqdARG_NONE},
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
  int    be_quiet;		/* TRUE to suppress verbose output & banner */
  int    do_local;		/* TRUE to config the model in local mode   */
  int    do_small;		/* TRUE to do divide and conquer alignments */
  int    windowlen;             /* window length for calculating bands */
  int    do_full;               /* TRUE to output all match columns in output alignment */
  int    do_apbanded;           /* TRUE to do a priori banded CYK (either d&c or full)*/
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

  /* a priori bands data structures */
  double   apbandp;	        /* tail loss probability for a priori banding */
  int    do_apexpand;           /* TRUE to naively expand a priori bands when necessary */
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
  int  do_hbanded;      /* TRUE to do HMM banded CYKInside_b_jd() using bands on d and j dim*/
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
  enum {		/* Type of HMM to use for banding */
    NONE,
    HMM_CP9,
    HMM_P7
  } hmm_type;
  char  *hmmfile;      /* file to read HMMs from                  */
  int   ks;            /* Counter over HMM state types (0 (match), 1(ins) or 2 (del))*/
  float forward_sc; 
  float backward_sc; 
  int   hmm_M;         /* Number of nodes in either the Plan 7 HMM or CM Plan 9 HMM */
  
  /* CM Plan 9 */
  CP9HMMFILE            *cp9_hmmfp;     /* opened CP9 hmmfile for reading                       */
  struct cplan9_s       *cp9_hmm;       /* constructed CP9 HMM; written to hmmfile              */
  int                    cp9_M;         /* number of nodes in CP9 HMM (MATL+MATR+2*MATP)        */
  int                    read_cp9_flag; /* TRUE to read a CP9 HMM from file, FALSE to build one */
  struct cp9_dpmatrix_s *cp9_mx;        /* growable DP matrix for viterbi                       */
  struct cp9_dpmatrix_s *cp9_fwd;       /* growable DP matrix for forward                       */
  struct cp9_dpmatrix_s *cp9_bck;       /* growable DP matrix for backward                      */
  struct cp9_dpmatrix_s *cp9_posterior; /* growable DP matrix for posterior decode              */
  float                  swentry;	/* S/W aggregate entry probability       */
  float                  swexit;        /* S/W aggregate exit probability        */
  int                    do_checkcp9;   /* TRUE to check the CP9 HMM by generating sequences */
  int                    seed;	        /* random number generator seed (only used if(do_checkcp9) */

  /* Plan 7 */
  HMMFILE           *hmmfp;     /* opened hmmfile for reading              */
  struct plan7_s    *hmm;       /* HMM to align to                         */ 
  struct dpmatrix_s *mx;        /* growable DP matrix                      */
  struct dpmatrix_s *fwd;       /* growable DP matrix for forwards         */
  struct dpmatrix_s *bck;       /* growable DP matrix for backwards        */
  struct dpmatrix_s *posterior; /* growable DP matrix for posterior decode */
  struct p7trace_s **p7tr;      /* traces for aligned sequences            */
  /*char             **p7postcode;*//* posterior decode array of strings       */
  unsigned char    **p7dsq;     /* digitized RNA sequences (plan 7 version)*/

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
  int                do_sub;      /* TRUE to use HMM to infer start and end point of hit 
				   * and remove structure before start and after end */
  int                hmm_start_node; /* HMM node most likely to have emitted posn 1 of target seq */
  int                hmm_start_state;/* HMM state type for hmm_start_node 0=match, 1=insert */
  int                hmm_end_node;   /* HMM node most likely to have emitted posn L of target seq */
  int                hmm_end_state;  /* HMM state type for hmm_end_node 0=match, 1=insert */
  int                **orig2sub_smap;/* 2D state map from orig_cm (template) to sub_cm.
				      * 1st dimension - state index in orig_cm 
				      * 2nd D - 2 elements for up to 2 matching sub_cm states */
  int                **sub2orig_smap;/* 2D state map from orig_cm (template) to sub_cm.
				      * 1st dimension - state index in sub_cm (0..sub_cm->M-1)
				      * 2nd D - 2 elements for up to 2 matching orig_cm states */
  char            *check_outfile;	/* output file name for subCM stk file*/
  FILE            *check_ofp;         /* an open output file */
  CM_t            *sub_cm;       /* sub covariance model                      */
  int *sub_node_cc_left; 
  int *sub_node_cc_right;
  int *sub_cc_node_map;
  int **sub_cs2hn_map;  
  int **sub_cs2hs_map;  
  int ***sub_hns2cs_map;
  struct cplan9_s       *sub_cp9_hmm;       /* constructed CP9 HMM; written to hmmfile              */
  struct cplan9_s       *sub_cp9_hmm2; 
  int                    sub_cp9_M;         /* number of nodes in CP9 HMM (MATL+MATR+2*MATP)        */
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
  windowlen   = 200;
  set_window  = FALSE;
  do_apbanded = FALSE;
  apbandp     = 0.0001;
  do_full     = FALSE;
  do_apexpand = FALSE;
  bdump_level = 0;
  debug_level = 0;
  do_hbanded  = FALSE;
  hbandp      = 0.0001;
  use_sums    = FALSE;
  hmm_type    = NONE;
  read_cp9_flag = FALSE;
  time_flag   = TRUE;
  do_inside   = FALSE;
  do_outside  = FALSE;
  do_check    = FALSE;
  do_post     = FALSE;
  do_sub      = FALSE;
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
    else if (strcmp(optname, "--apbanded")  == 0) do_apbanded = TRUE;
    else if (strcmp(optname, "--apbandp")   == 0) apbandp     = atof(optarg);
    else if (strcmp(optname, "-W")          == 0) {
      windowlen    = atoi(optarg); 
      set_window = TRUE; } 
    else if (strcmp(optname, "--full")      == 0) do_full      = TRUE;
    else if (strcmp(optname, "--apexpand")  == 0) do_apexpand  = TRUE;
    else if (strcmp(optname, "--banddump")  == 0) bdump_level  = atoi(optarg);
    else if (strcmp(optname, "--dlev")      == 0) debug_level  = atoi(optarg);
    else if (strcmp(optname, "--time")      == 0) time_flag    = TRUE;
    else if (strcmp(optname, "--inside")    == 0) do_inside    = TRUE;
    else if (strcmp(optname, "--outside")   == 0) { do_outside = TRUE; do_check = TRUE; }
    else if (strcmp(optname, "--post")      == 0) do_post      = TRUE;
    else if (strcmp(optname, "--checkpost") == 0) do_check     = TRUE;
    else if (strcmp(optname, "--sub")       == 0) 
      { 
	do_sub  = TRUE;
	hmm_type = HMM_CP9; 
      }

    else if (strcmp(optname, "--hbanded")   == 0) {
      do_hbanded = TRUE; do_small = FALSE; 
      if(hmm_type == NONE) hmm_type = HMM_CP9; 
      } /* default HMM: CP9, but let --p7 override */
    else if (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--checkcp9")  == 0) do_checkcp9  = TRUE;
    else if (strcmp(optname, "--cp9")       == 0) { 
      if(hmm_type == HMM_P7) Die("Can't do --cp9 and --p7, pick one HMM type.\n");
      hmm_type  = HMM_CP9; hmmfile = optarg; read_cp9_flag = TRUE; 
    }
    else if (strcmp(optname, "--p7")        == 0) {
      if(read_cp9_flag == TRUE) Die("Can't do --cp9 and --p7, pick one HMM type.\n");
      hmm_type  = HMM_P7;  hmmfile = optarg;  
    }
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
    {
      Die("Please pick either --inside or --outside (--outside will run Inside()\nalso and check to make sure Inside() and Outside() scores agree).\n");
    }

  if(do_sub && do_local)
    {
      Die("--sub and --local combination not yet supported.\n");
    }
  if(do_checkcp9 && (hmm_type != HMM_CP9 || do_hbanded == FALSE))
    {
      Die("--checkcp9 only makes sense with --hbanded (and not --p7).\n");
    }

  /* currently not set up for local alignment and posterior decode or outside run */
  /*  if(do_local && (do_outside || do_post))
    {
      Die("Can't do -l and either --post or --outside.\n");
    }
  */
  /* default to HMM banded alignment if --p7 or --cp9 was enabled */
  /*if(hmm_type != NONE)
    {
      do_hbanded = TRUE;
      do_small = FALSE;
    }
  */

  if (bdump_level > 3) Die("Highest available --banddump verbosity level is 3\n%s", usage);
  if (do_apexpand && (!(do_apbanded))) Die("Doesn't make sense to use --apexpand option without --apbanded option\n", usage);
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
  if(hmm_type != NONE)
    {
      Alphabet_type = hmmNOTSETYET;
      if(hmm_type == HMM_P7)
	{
	  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
	    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
	  if (!HMMFileRead(hmmfp, &hmm)) 
	    Die("Failed to read any HMMs from %s\n", hmmfile);
	  HMMFileClose(hmmfp);
	  hmm->xt[XTE][MOVE] = 1.;     /* only 1 domain/sequence ("global" alignment) */
	  hmm->xt[XTE][LOOP] = 0.;
	  P7Logoddsify(hmm, TRUE);
	  if (hmm == NULL) 
	    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
	}
      if(hmm_type == HMM_CP9)
	{
	  if(read_cp9_flag)
	    {
	      if ((cp9_hmmfp = CP9_HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
		Die("Failed to open HMM file %s\n%s", hmmfile, usage);
	      if (!CP9_HMMFileRead(cp9_hmmfp, &cp9_hmm)) 
		Die("Failed to read any CP9 HMMs from %s\n", hmmfile);
	      CP9_HMMFileClose(cp9_hmmfp);
	      CP9Logoddsify(cp9_hmm);
	      if (cp9_hmm == NULL) 
		Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
	    }
	  else /* build (don't read in) a CM Plan 9 from the CM (default --hbanded) */
	    SetAlphabet(hmmNUCLEIC); /* Set up the hmmer_alphabet global variable */
	}
    }

  /*****************************************************************
   * Input and digitize the unaligned sequences
   *****************************************************************/

  if (! ReadMultipleRseqs(seqfile, format, &rseq, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);
  dsq = MallocOrDie(sizeof(char *) * nseq);
  for (i = 0; i < nseq; i++) 
    dsq[i] = DigitizeSequence(rseq[i], sqinfo[i].len);
  if(hmm_type != NONE)
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

  if(hmm_type == HMM_P7)
    {
      /*p7tr  = MallocOrDie(sizeof(struct p7trace_s *) * nseq);*/
      mx  = CreatePlan7Matrix(1, hmm->M, 25, 0);
      /*p7postcode = malloc(sizeof(char *) * nseq);*/
      map_consensus_columns(cm, hmm->M, &node_cc_left, &node_cc_right,
			    &cc_node_map, debug_level);
      /*printf("HMM type is HMM_P7, number of nodes: %d\n", hmm->M);*/
      /* Get information mapping the HMM to the CM and vice versa. */
      
      P7_map_cm2hmm_and_hmm2cm(cm, hmm, node_cc_left, node_cc_right, 
			       cc_node_map, &cs2hn_map, &cs2hs_map, 
			       &hns2cs_map, debug_level);
      hmm_M = hmm->M;
    }      
  
  if(hmm_type == HMM_CP9)
    {
      if(!(read_cp9_flag)) /* we're building a CP9 HMM, not reading one in */
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
	  cp9_hmm = AllocCPlan9(cp9_M);
	  ZeroCPlan9(cp9_hmm);
	}

      /* Get information mapping the HMM to the CM and vice versa, used
       * for mapping bands. 
       * This is relevant if we're building (read_cp9_flag == FALSE)
       * or reading from a file (read_cp9_flag == TRUE)
       */
      map_consensus_columns(cm, cp9_hmm->M, &node_cc_left, &node_cc_right,
			    &cc_node_map, debug_level);
      
      CP9_map_cm2hmm_and_hmm2cm(cm, cp9_hmm, node_cc_left, node_cc_right, 
				cc_node_map, &cs2hn_map, &cs2hs_map, 
				&hns2cs_map, debug_level);
      
      if(!(read_cp9_flag)) /* actually build the CM Plan 9 HMM */
      {
	/* fill in parameters of HMM using the CM and some ideas/formulas/tricks
	 * from Zasha Weinberg's thesis (~p.123) (see CP9_cm2wrhmm.c code) */
	if(!(CP9_cm2wrhmm(cm, cp9_hmm, node_cc_left, node_cc_right, cc_node_map, cs2hn_map,
		     cs2hs_map, hns2cs_map, debug_level)))
	  Die("Couldn't build a CM Plan 9 HMM from the CM.\n");
	if(do_checkcp9)
	  {
	    sre_srandom(seed);
	    //if(!(CP9_check_wrhmm_by_sampling(cm, cp9_hmm, 1, cp9_hmm->M, hns2cs_map, 0.05, 100000)))
	    if(!(CP9_check_wrhmm_by_sampling(cm, cp9_hmm, 1, cp9_hmm->M, hns2cs_map, 0.01, 100000, NULL,
					     NULL, NULL, debug_level)))
	      Die("CM Plan 9 fails sampling check!\n");
	    else
	      printf("CM Plan 9 passed sampling check.\n");
	  }
      }
      else /* a CP9 HMM has already been read from a file. */
	{
	/* Check to make sure the HMM we read from the input file is
	 * "close enough" to the CM, based on psi and phi values (see CP9_cm2wrhmm.c).
	 * If not, we build a new one from the CM. This is a design decision, 
	 * and I'm not sure if its the best or desired behavior (EPN)
	 */
	  if(!(CP9_check_wrhmm(cm, cp9_hmm, hns2cs_map, cc_node_map, debug_level)))
	    {
	      printf("\nCM Plan 9 HMM read from %s not similar enough to the CM.\n", hmmfile);
	      printf("Building a new CM plan 9 HMM directly from the CM.\n\n");
	      /* build a new CM plan 9 HMM */
	      if(!(CP9_cm2wrhmm(cm, cp9_hmm, node_cc_left, node_cc_right, cc_node_map, cs2hn_map,
				cs2hs_map, hns2cs_map, debug_level)))
		Die("Couldn't build a CM Plan 9 HMM from the CM\n");
	    }
	}
      CP9Logoddsify(cp9_hmm);
      fill_phi_cp9(cp9_hmm, &orig_phi, 1);
      /*debug_print_cp9_params(cp9_hmm);*/

      cp9_mx  = CreateCPlan9Matrix(1, cp9_hmm->M, 25, 0);
      hmm_M = cp9_hmm->M;

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
      if(hmm_type == HMM_CP9)
	{
	  printf("configuring the CM plan 9 HMM for local alignment.\n");
	  swentry           = 0.5;
	  swexit            = 0.5;
	  CPlan9SWConfig(cp9_hmm, swentry, swexit);
	  CP9Logoddsify(cp9_hmm);
	}
    }

  /* set up the a priori bands, this has to be done after the ConfigLocal() call */
  if(do_apbanded || bdump_level > 0)
    {
      safe_windowlen = windowlen * 2;
      while(!(BandCalculationEngine(cm, safe_windowlen, apbandp, 0, &dmin, &dmax, &gamma, do_local)))
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
       * is > cm->W. This should only happen if the apbandp we're using now is < 1E-7 
       * (1E-7 is the apbandp value used to determine cm->W in cmbuild). If this 
       * happens, the current implementation reassigns windowlen to this larger value.
       * NOTE: if W was set at the command line, the command line value is 
       *       always used.
       */
      if(!(set_window))
	{
	  windowlen = dmax[0];
	}
      if(bdump_level > 1) 
	{
	  printf("apbandp:%f\n", apbandp);
	  debug_print_bands(cm, dmin, dmax);
	}
      expand_flag = FALSE;
      PrintDPCellsSaved(cm, dmin, dmax, windowlen);
    }

  if (do_sub && hmm_type == HMM_CP9)
    {
      printf("configuring the CM plan 9 HMM for local alignment.\n");
      swentry           = 0.5;
      swexit            = 0.5;
      printf("PRINTING ORIG SUB CP9 HMM PARAMETERS GLOCAL 0\n");
      debug_print_cp9_params(cp9_hmm);
      printf("DONE PRINTING ORIG SUB CP9 HMM PARAMETERS GLOCAL 0\n\n");
      CPlan9SWConfig(cp9_hmm, swentry, swexit);
      printf("PRINTING ORIG SUB CP9 HMM PARAMETERS LOCAL 0\n");
      debug_print_cp9_params(cp9_hmm);
      printf("DONE PRINTING ORIG SUB CP9 HMM PARAMETERS LOCAL 0\n\n");
      CP9Logoddsify(cp9_hmm);
    }
  /* Allocate data structures for use with HMM banding strategy (either P7 or CP9) */
  if(hmm_type != NONE)
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
      if(hmm_type == HMM_CP9)
	{
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);
	  /* Align the current seq to the cp9 HMM, we don't care
	   * about the trace, just the posteriors.
	   * Step 1: Get HMM posteriors.
	   * Step 2: posteriors -> HMM bands.
	   * Step 3: HMM bands  ->  CM bands.
	   */
	  
	  /* Step 1: Get HMM posteriors.*/
	  /*sc = CP9Viterbi(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, cp9_mx);*/
	  forward_sc = CP9Forward(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, &cp9_fwd);
	  printf("CP9 i: %d | forward_sc : %f\n", i, forward_sc);
	  backward_sc = CP9Backward(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, &cp9_bck);
	  printf("CP9 i: %d | backward_sc: %f\n", i, backward_sc);

	  /*debug_check_CP9_FB(cp9_fwd, cp9_bck, cp9_hmm, forward_sc, 1, sqinfo[i].len, p7dsq[i]);*/
	  cp9_posterior = cp9_bck;
	  CP9FullPosterior(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, cp9_fwd, cp9_bck, cp9_posterior);

	  /* If we're in sub mode:
	   * (1) infer the start and end HMM states by looking at the posterior matrix.
	   * (2) modify the CM model by removing unnecessary structure and marginalizing
	   * (3) (NECESSARY?) build a new CP9 HMM from the modified CM
	   * (4) (NECESSARY?) do Forward/Backward again, and get a new posterior matrix.
	   */
	  if(do_sub)
	    {
	      /* we're necessarily in local mode, the --sub option turns local mode on */
	      CP9NodeForPosn(cp9_hmm, 1, sqinfo[i].len, 1,             cp9_posterior, &hmm_start_node, &hmm_start_state);
	      CP9NodeForPosn(cp9_hmm, 1, sqinfo[i].len, sqinfo[i].len, cp9_posterior, &hmm_end_node,   &hmm_end_state);
	      
	      /* Given the original CM, and the start and end HMM nodes, build a new CM by removing
	       * structure outside the start and end HMM nodes, and marginalizing.
	       */
	      /* Uncomment below to build a sub_cm with same number of consensus columns, just structure
	       * removed. */
	      /*build_sub_cm(cm, &sub_cm, hmm_start_node, hmm_end_node, 1, cp9_hmm->M, orig2sub_smap, sub2orig_smap);*/

	      /* Uncomment below to build a sub_cm that only models consensus columns between HMM start and end
	       * node. */
	      build_sub_cm(cm, &sub_cm, hmm_start_node, hmm_end_node, orig2sub_smap, sub2orig_smap, 
			   &imp_cc, NULL, NULL, NULL, NULL, 0.00001, FALSE, FALSE, FALSE, 0.01, 
			   10000, debug_level);

	      exit(1);

	      /* Following function call samples for cm and sub_cm and builds CP9 HMMs from each set of samples,
	       * then prints out the parameters of those HMMs.
	      check_sub_cm_by_sampling2(cm, sub_cm, hmm_start_node, hmm_end_node, 10000);
	      exit(1);*/

	      /**************************************************************************/
	      /* For checking if our CP9->subCP9 truncation modification method is working we build a
	       * CP9 HMM from the sub_cm using the full blown method, and check its parameters
	       * against the subCP9 built from the truncation-built one.
	       *****************BEGIN EVENTUALLY ERASABLE BLOCK*************************/
	      sub_cp9_M = 0;
	      for(v = 0; v <= sub_cm->M; v++)
		{
		  if(sub_cm->stid[v] ==  MATP_MP)
		    sub_cp9_M += 2;
		  else if(sub_cm->stid[v] == MATL_ML || sub_cm->stid[v] == MATR_MR)
		    sub_cp9_M ++;
		}

	      StopwatchZero(watch1);
	      StopwatchStart(watch1);
	      /* build the HMM data structure */
	      sub_cp9_hmm = AllocCPlan9(sub_cp9_M);
	      printf("sub_cp9_M: %d\n", sub_cp9_M);
	      ZeroCPlan9(sub_cp9_hmm);

	      /* Get information mapping the HMM to the CM and vice versa, used
	       * for mapping bands. */
	      map_consensus_columns(sub_cm, sub_cp9_hmm->M, &sub_node_cc_left, &sub_node_cc_right,
				    &sub_cc_node_map, debug_level);
	      
	      printf("check 1\n");
	      CP9_map_cm2hmm_and_hmm2cm(sub_cm, sub_cp9_hmm, sub_node_cc_left, sub_node_cc_right, 
					sub_cc_node_map, &sub_cs2hn_map, &sub_cs2hs_map, 
					&sub_hns2cs_map, debug_level);
	      if(!(CP9_cm2wrhmm(sub_cm, sub_cp9_hmm, sub_node_cc_left, sub_node_cc_right, sub_cc_node_map, 
				sub_cs2hn_map, sub_cs2hs_map, sub_hns2cs_map, debug_level)))
		Die("Couldn't build a CM Plan 9 HMM from the sub CM.\n");
	      StopwatchStop(watch1);
	      StopwatchDisplay(stdout, "CP9 BUILDING TRADITIONAL TIME: ", watch1);

	      printf("check 2\n");
	      printf("sub_cp9_hmm->M: %d\n", sub_cp9_hmm->M);
	      printf("PRINTING ORIG SUB CP9 HMM PARAMETERS LOCAL\n");
	      debug_print_cp9_params(cp9_hmm);
	      printf("DONE PRINTING ORIG SUB CP9 HMM PARAMETERS LOCAL\n\n");

	      /* Reconfig the cp9_hmm as glocal */
	      CPlan9GlobalConfig(cp9_hmm);

	      printf("PRINTING ORIG SUB CP9 HMM PARAMETERS GLOCAL 1\n");
	      debug_print_cp9_params(cp9_hmm);
	      printf("DONE PRINTING ORIG SUB CP9 HMM PARAMETERS GLOCAL 1\n\n");

	      /*sub_CPlan9GlobalConfig(cp9_hmm, hmm_start_node, hmm_end_node, orig_phi);*/

	      StopwatchZero(watch1);
	      StopwatchStart(watch1);
	      /* Build the sub CP9 HMM by copying as much of the original cp9_hmm as possible */
	      CP9_2sub_cp9(cp9_hmm, &sub_cp9_hmm2, hmm_start_node, hmm_end_node, orig_phi);
	      StopwatchStop(watch1);
	      StopwatchDisplay(stdout, "CP9 BUILDING TRUNCATION TIME: ", watch1);

	      /* fill phi arrays for each CP9 HMM */
	      fill_phi_cp9(sub_cp9_hmm, &sub_phi_trad, 1);
	      fill_phi_cp9(sub_cp9_hmm2, &sub_phi_trunc, 1);
	      
	      printf("sub_cp9_hmm->M: %d\n", sub_cp9_hmm->M);

	      printf("MODEL START: %d\n", hmm_start_node);
	      printf("MODEL END  : %d\n", hmm_end_node);

	      printf("PRINTING TRAD SUB CP9 HMM PARAMETERS\n");
	      debug_print_cp9_params(sub_cp9_hmm);
	      printf("DONE PRINTING TRAD SUB CP9 HMM PARAMETERS\n\n");

	      printf("PRINTING TRUNC SUB CP9 HMM PARAMETERS\n");
	      debug_print_cp9_params(sub_cp9_hmm2);
	      printf("DONE PRINTING TRUNC SUB CP9 HMM PARAMETERS\n\n");
	      
	      printf("PRINTING ORIG CP9 HMM PHI PARAMETERS\n");
	      debug_print_phi_cp9(cp9_hmm, orig_phi);
	      printf("DONE PRINTING ORIG CP9 HMM PHI PARAMETERS\n");

	      printf("PRINTING TRAD SUB CP9 HMM PHI PARAMETERS\n");
	      debug_print_phi_cp9(sub_cp9_hmm, sub_phi_trad);
	      printf("DONE PRINTING TRAD SUB CP9 HMM PHI PARAMETERS\n\n");

	      printf("PRINTING TRUNC SUB CP9 HMM PHI PARAMETERS\n");
	      debug_print_phi_cp9(sub_cp9_hmm2, sub_phi_trunc);
	      printf("DONE PRINTING TRUNC SUB CP9 HMM PHI PARAMETERS\n\n");

	      exit(1);
	      
	      /* Derive bands for the sub_cm using the newly configured CP9 HMM */
	      /* Step 1: Get HMM posteriors.*/
	      /*sc = CP9Viterbi(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, cp9_mx);*/
	      forward_sc = CP9Forward(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, &cp9_fwd);
	      printf("CP9 i: %d | forward_sc : %f\n", i, forward_sc);
	      backward_sc = CP9Backward(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, &cp9_bck);
	      printf("CP9 i: %d | backward_sc: %f\n", i, backward_sc);
	      
	      /*debug_check_CP9_FB(cp9_fwd, cp9_bck, cp9_hmm, forward_sc, 1, sqinfo[i].len, p7dsq[i]);*/
	      cp9_posterior = cp9_bck;
	      CP9FullPosterior(p7dsq[i], 1, sqinfo[i].len, cp9_hmm, cp9_fwd, cp9_bck, cp9_posterior);
	      
	      exit(1);
	      /*cm = sub_cm;*/
	      /* Build a new CM Plan 9 HMM */

	    }
	  StopwatchStop(watch1);
	  if(time_flag) StopwatchDisplay(stdout, "CP9 Forward/Backward CPU time: ", watch1);
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);
	  
	  /* Step 2: posteriors -> HMM bands.*/
	  if(!(use_sums))
	    {
	      /* match states */
	      CP9_hmm_band_bounds(cp9_posterior->mmx, 1, sqinfo[i].len, cp9_hmm->M,
				  NULL, pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
	      /* insert states */
	      CP9_hmm_band_bounds(cp9_posterior->imx, 1, sqinfo[i].len, cp9_hmm->M,
				  NULL, pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
	      /* delete states (note: delete_flag set to TRUE) */
	      CP9_hmm_band_bounds(cp9_posterior->dmx, 1, sqinfo[i].len, cp9_hmm->M,
				  NULL, pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);
	    }
	  else
	    {
	      CP9_ifill_post_sums(cp9_posterior, 1, sqinfo[i].len, cp9_hmm->M,
				  isum_pn_m[i], isum_pn_i[i], 
				  isum_pn_d[i]);
	      /* match states */
	      CP9_hmm_band_bounds(cp9_posterior->mmx, 1, sqinfo[i].len, cp9_hmm->M,
				  isum_pn_m[i], pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
	      /* insert states */
	      CP9_hmm_band_bounds(cp9_posterior->imx, 1, sqinfo[i].len, cp9_hmm->M,
				  isum_pn_i[i], pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
	      /* delete states */
	      CP9_hmm_band_bounds(cp9_posterior->dmx, 1, sqinfo[i].len, cp9_hmm->M,
				  isum_pn_d[i], pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);
	    }
	  if(debug_level != 0)
	    {
	      printf("printing hmm bands\n");
	      print_hmm_bands(stdout, sqinfo[i].len, cp9_hmm->M, pn_min_m, pn_max_m, pn_min_i,
			      pn_max_i, pn_min_d, pn_max_d, hbandp, debug_level);
	    }
	  
	  /* Step 3: HMM bands  ->  CM bands. */
	  hmm2ij_bands(cm, cp9_hmm->M, node_cc_left, node_cc_right, cc_node_map, 
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

      else if(hmm_type == HMM_P7)
	{
	  /* Align the current seq to the p7 HMM, we don't care
	   * about the trace, just the posteriors.
	   * Step 1: Get HMM posteriors.
	   * Step 2: posteriors -> HMM bands.
	   * Step 3: HMM bands  ->  CM bands.
	   */
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);
	  
	  /* Step 1: Get HMM posteriors.*/
	  /*
	  if (P7ViterbiSpaceOK(sqinfo[i].len, hmm->M, mx))
	    (void) P7Viterbi(p7dsq[i], sqinfo[i].len, hmm, mx, &(p7tr[i]));
	  else
	    (void) P7SmallViterbi(p7dsq[i], sqinfo[i].len, hmm, mx, &(p7tr[i]));
	  */
	  forward_sc = P7Forward(p7dsq[i], sqinfo[i].len, hmm, &fwd);
	  /*printf("P7 i: %d | forward_sc : %f\n", i, forward_sc);*/
	  backward_sc = P7Backward(p7dsq[i], sqinfo[i].len, hmm, &bck);
	  /*printf("P7 i: %d | backward_sc: %f\n", i, backward_sc);*/
	  posterior = bck;
	  P7FullPosterior(sqinfo[i].len, hmm, fwd, bck, posterior);
	  if(debug_level > 2)
	    P7_debug_print_post_decode(sqinfo[i].len, hmm->M, posterior);
	  
	  StopwatchStop(watch1);
	  if(time_flag) StopwatchDisplay(stdout, "P7 Forward/Backward CPU time: ", watch1);
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);
	  
	  /* Uncomment this to assign a postal code (IHH's postprob.c) */
	  /*p7postcode[i] = PostalCode(sqinfo[i].len, posterior, p7tr[i]);
	  printf("p7postcode[i]\n%s\n", p7postcode[i]);
	  */

	  /* Step 2: posteriors -> HMM bands.*/
	  if(!(use_sums))
	    {
	      /* match states */
	      P7_hmm_band_bounds(posterior->mmx, sqinfo[i].len, hmm->M, NULL,
				 pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
	      /* insert states */
	      P7_hmm_band_bounds(posterior->imx, sqinfo[i].len, hmm->M, NULL,
				 pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
	      /* delete states */
	      P7_hmm_band_bounds(posterior->dmx, sqinfo[i].len, hmm->M, NULL,
				 pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);

	    }
	  else
	    {
	      P7_ifill_post_sums(posterior, sqinfo[i].len, hmm->M,
				 isum_pn_m[i], isum_pn_i[i], 
				 isum_pn_d[i]);
	      /* match states */
	      P7_hmm_band_bounds(posterior->mmx, sqinfo[i].len, hmm->M, isum_pn_m[i],
				 pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
	      /* insert states */
	      P7_hmm_band_bounds(posterior->imx, sqinfo[i].len, hmm->M, isum_pn_i[i],
				 pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
	      /* delete states (note: delete_flag set to TRUE) */
	      P7_hmm_band_bounds(posterior->dmx, sqinfo[i].len, hmm->M, isum_pn_d[i],
				 pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);
	    }

	  /* The pn_min_i and pn_max_i arrays are not complete.
	   * HMMER Plan 7 doesn't allow inserts after the final match state,
	   * or deletes before the first match, or after the final match state,
	   * Therefore pn_min_i and pn_max_i are only filled from 1..M-1.
	   * But we need the band for HMM node M as well. We use a hack
	   * and assign it the same as the pn_max_m[M] and pn_min_m[M].
	   */
	  P7_last_hmm_insert_state_hack(hmm->M, pn_min_m, pn_max_m, pn_min_i, pn_max_i);
	  P7_last_and_first_hmm_delete_state_hack(hmm->M, pn_min_m, pn_max_m, pn_min_d, 
						  pn_max_d, sqinfo[i].len);

	  if(debug_level != 0)
	    {
	      printf("printing hmm bands\n");
	      print_hmm_bands(stdout, sqinfo[i].len, hmm->M, pn_min_m, pn_max_m, pn_min_i,
			      pn_max_i, pn_min_d, pn_max_d, hbandp, debug_level);
	      /* Get the bands imin, imax, jmin, and jmax for the CM */
	    }

	  /* Step 3: HMM bands  ->  CM bands. */
	  hmm2ij_bands(cm, hmm->M, node_cc_left, node_cc_right, cc_node_map, 
		       1, sqinfo[i].len, pn_min_m, pn_max_m, pn_min_i, pn_max_i, 
		       pn_min_d, pn_max_d, imin, imax, jmin, jmax, cs2hn_map,
		       debug_level);
	  
	  StopwatchStop(watch1);
	  if(time_flag) StopwatchDisplay(stdout, "P7 Band calculation CPU time: ", watch1);

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
	
	  FreePlan7Matrix(fwd);
	  FreePlan7Matrix(posterior);
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
      
      if(do_apexpand)
	{
	  /*First, check to see if we need to reset the apriori bands b/c 
	   * they're currently expanded. */
	  if(expand_flag)
	    {
	      FreeBandDensities(cm, gamma);
	      free(dmin);
	      free(dmax);
	      while(!(BandCalculationEngine(cm, safe_windowlen, apbandp, 0, &dmin, &dmax, &gamma, do_local)))
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
	  if (do_apbanded && (sqinfo[i].len < dmin[0] || sqinfo[i].len > dmax[0]))
	    {
	      /* the seq we're aligning is outside the root band, but
	       * --apexpand was not enabled, so we die.*/
	      Die("Length of sequence to align (%d nt) lies outside the root band.\ndmin[0]: %d and dmax[0]: %d\nImpossible to align with a priori banded CYK unless you try --apexpand.\n%s", sqinfo[i].len, dmin[0], dmax[0], usage);
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
	  if(do_apbanded)
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
      else if(do_apbanded)
	{
	  sc = CYKInside_b(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), dmin, dmax);
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
	  sc = CYKInside(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]));
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

      /* If debug level high enough, print out the parse tree */
      if(debug_level > 2)
	{
	  fprintf(stdout, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], dsq[i], FALSE));;
	  ParsetreeDump(stdout, tr[i], cm, dsq[i]);
	  fprintf(stdout, "//\n");
	}
      /* Dump the trace with info on i, j and d bands
       * if bdump_level is high enough or HMM bands were calc'ed but not 
       * used (useful in the latter case to determine if the optimal trace 
       * is within our bands on i and j and d*/
      if((hmm_type != NONE && !do_hbanded) || bdump_level > 0)
	{
	  /*printf("bands on i\n");
	  debug_print_bands(cm, imin, imax);
	  printf("bands on j\n");
	  debug_print_bands(cm, jmin, jmax);
	  */
	  ijd_banded_trace_info_dump(cm, tr[i], imin, imax, jmin, jmax, hdmin, hdmax, 1);
	}
      
      /* Clean up the structures we use calculating HMM bands, that are allocated
       * differently for each sequence. 
       */
      if(hmm_type != NONE)
	for(v = 0; v < cm->M; v++)
	  { free(hdmin[v]); free(hdmax[v]); }

      StopwatchStop(watch2);
      if(time_flag) StopwatchDisplay(stdout, "band calc and jd CYK CPU time: ", watch2);
    }
  avgsc /= nseq;

  if(!(do_inside || do_outside))
    {
      msa = Parsetrees2Alignment(cm, dsq, sqinfo, NULL, tr, nseq, do_full);
      
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
  if (do_apbanded)
    {
      FreeBandDensities(cm, gamma);
      free(dmin);
      free(dmax);
    }

  if(hmm_type != NONE)
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
      free(node_cc_left);
      free(node_cc_right);
      free(cc_node_map);
      free(safe_hdmin);
      free(safe_hdmax);
      for (i = 0; i < nseq; i++) 
	{ 
	  free(p7dsq[i]);
	  free(isum_pn_m[i]); 
	  free(isum_pn_i[i]);
	  free(isum_pn_d[i]);
	  /* Uncomment if calc'ing a postcode */
	  if(hmm_type == HMM_P7) 
	    {
	      /*P7FreeTrace(p7tr[i]);*/
	      /*free(p7postcode[i]);*/
	    }
	}
      free(p7dsq);
      free(isum_pn_m);
      free(isum_pn_i);
      free(isum_pn_d);
      for(v = 0; v <= cm->M; v++)
	{
	  free(cs2hn_map[v]);
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
    }
  if(hmm_type == HMM_P7)
    {
      /*free(p7tr);*/
      /*free(p7postcode);*/
      FreePlan7(hmm);
      FreePlan7Matrix(mx);
    }      
  if(hmm_type == HMM_CP9)
    {
      FreeCPlan9Matrix(cp9_mx);
      FreeCPlan9(cp9_hmm);
    }

  if(!(do_inside || do_outside)) MSAFree(msa);
  FreeCM(cm);
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


/* cmsearch.c
 * SRE, Fri May  3 13:58:18 2002
 * SVN $Id$
 * 
 * Search sequences with a CM.
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

#include "squid.h"		/* general sequence analysis library    */
#include "easel.h"              /* better general sequence analysis library */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process tcp9b->iming module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "hmmband.h"
#include "stats.h"
#include "esl_gumbel.h"
#include "esl_sqio.h"
#include "mpifuncs.h"

/* From rsearch-1.1 various defaults defined here */
#define DEFAULT_NUM_SAMPLES 0
#define DEFAULT_CUTOFF 0.0
#define DEFAULT_CUTOFF_TYPE SCORE_CUTOFF

static int in_mpi;

#ifdef USE_MPI
/*
 * Function: exit_from_mpi
 * Date:     RJK, Thu Jun 6, 2002 [St. Louis]
 * Purpose:  Calls MPI_Abort on exit if in_mpi flag is 1, otherwise
 *           returns
 */
void exit_from_mpi () {
  if (in_mpi)
    MPI_Abort (MPI_COMM_WORLD, -1);
}
#endif

static int QDBFileRead(FILE *fp, CM_t *cm, int **ret_dmin, int **ret_dmax);
static int set_partitions(int **ret_partitions, int *num_partitions, char *list);
static int debug_print_stats(int *partitions, int num_partitions, double *lambda, double *mu);

static char banner[] = "cmsearch - search a sequence database with an RNA covariance model";

static char usage[]  = "\
Usage: cmsearch [-options] <cmfile> <sequence file>\n\
The sequence file is expected to be in FASTA format.\n\
  Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -s <f> : use cutoff bit score of <f> [default: 0]\n\
   -E <f> : use cutoff E-value of <f> (default ignored; not-calc'ed)\n\
   -n <n> : determine EVD with <n> samples (default with -E: 1000)\n\
   -W <n> : set scanning window size to <n> (default: precalc'd in cmbuild)\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input alignment is in format <s>, not FASTA\n\
   --toponly     : only search the top strand\n\
   --local       : do local alignment\n\
   --noalign     : find start/stop only; don't do alignments\n\
   --dumptrees   : dump verbose parse tree information for each hit\n\
   --partition <n>[,<n>]... : partition points for different GC content EVDs\n\
   --inside      : scan with Inside, not CYK (caution much slower(!))\n\
   --null2       : turn on the post hoc second null model [df:OFF]\n\
   --learninserts: do not set insert emission scores to 0\n\
   --negscore <f>: set min bit score to report as <f> < 0 (experimental)\n\
   --enfstart <n>: enforce MATL stretch starting at CM node <n>\n\
   --enfseq   <s>: enforce MATL stretch starting at --enfstart <n> emits seq <s>\n\
\n\
  * Filtering options using a CM plan 9 HMM (*in development*):\n\
   --hmmfb        : use Forward to get end points & Backward to get start points\n\
   --hmmweinberg  : use Forward to get end points, subtract W for start points\n\
   --hmmpad <n>   : subtract/add <n> residues from start/end [df:0]\n\
   --hmmonly      : don't use CM at all, just scan with HMM (Forward + Backward)\n\
   --hthresh <f>  : HMM reporting bit score threshold [df: 0]\n\
\n\
  * Options for accelerating CM search/alignment (*in development*):\n\
   --beta <f>    : tail loss prob for QBD (default:0.0000001)\n\
   --noqdb       : DO NOT use query dependent bands (QDB) to accelerate CYK\n\
   --qdbfile <f> : read QDBs from file <f> (outputted from cmbuild)\n\
   --hbanded     : use HMM bands from a CM plan 9 HMM scan for CYK\n\
   --hbandp <f>  : tail loss prob for --hbanded (default:0.0001)\n\
   --banddump    : print bands for each state\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
   --scan2bands  : use scanning Forward and Backward to get bands (EXPTL!)\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-s", TRUE, sqdARG_FLOAT }, 
  { "-E", TRUE, sqdARG_FLOAT }, 
  { "-n", TRUE, sqdARG_INT }, 
  { "-W", TRUE, sqdARG_INT }, 
  { "--dumptrees",  FALSE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--local",      FALSE, sqdARG_NONE },
  { "--noalign",    FALSE, sqdARG_NONE },
  { "--toponly",    FALSE, sqdARG_NONE },
  { "--inside",     FALSE, sqdARG_NONE },
  { "--null2",      FALSE, sqdARG_NONE },
  { "--learninserts",FALSE, sqdARG_NONE},
  { "--negscore",   FALSE, sqdARG_FLOAT},
  { "--hmmfb",      FALSE, sqdARG_NONE },
  { "--hmmweinberg",FALSE, sqdARG_NONE},
  { "--hmmpad",     FALSE, sqdARG_INT },
  { "--hmmonly",    FALSE, sqdARG_NONE },
  { "--hthresh",    FALSE, sqdARG_FLOAT},
  { "--noqdb",      FALSE, sqdARG_NONE },
  { "--qdbfile",    FALSE, sqdARG_STRING},
  { "--beta",       FALSE, sqdARG_FLOAT},
  { "--hbanded",    FALSE, sqdARG_NONE },
  { "--hbandp",     FALSE, sqdARG_FLOAT},
  { "--banddump",   FALSE, sqdARG_NONE},
  { "--sums",       FALSE, sqdARG_NONE},
  { "--scan2hbands",FALSE, sqdARG_NONE},
  { "--partition",  FALSE, sqdARG_STRING},
  { "--enfstart",   FALSE, sqdARG_INT},
  { "--enfseq",     FALSE, sqdARG_STRING}
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  char            *seqfile;     /* file to read sequences from */
  int              format = eslSQFILE_UNKNOWN;   /* format of sequence file */
  int              status;
  CMFILE          *cmfp;        /* open CM file for reading */
  int             querylen;     /* consensus length of the CM query */
  ESL_SQFILE	  *dbfp;        /* open seqfile for reading */
  CM_t            *cm;          /* a covariance model       */
  ESL_SQ           *sq;
  char            *dsq;         /* digitized RNA sequence */
  Stopwatch_t     *watch;       /* times band calc, then search time */
  int              i, ip;
  int              reversed;    /* TRUE when we're doing the reverse complement strand */
  int              maxlen;
  Parsetree_t     *tr;		/* parse of an individual hit */
  CMConsensus_t   *cons;	/* precalculated consensus info for display */
  Fancyali_t      *ali;         /* alignment, formatted for display */

  double  **gamma;              /* P(subseq length = n) for each state v    */
  int     *dmin = NULL;		/* minimum d bound for state v, [0..v..M-1] */
  int     *dmax = NULL; 	/* maximum d bound for state v, [0..v..M-1] */
  double   beta;		/* tail loss probability for a priori banding */

  /* information on hits found with the CM */
  int    nhits;			/* number of hits in a seq */
  int   *hitr;			/* initial states for hits */
  int   *hiti;                  /* start positions of hits */
  int   *hitj;                  /* end positions of hits */
  float *hitsc;			/* scores of hits */

  /* information on hits found with the CP9 HMM derived from the CM */
  int    hmm_nhits;		/* number of hits in a seq */
  int   *hmm_hitr;		/* initial states for hits */
  int   *hmm_hiti;              /* start positions of hits */
  int   *hmm_hitj;              /* end positions of hits */
  float *hmm_hitsc;		/* scores of hits */

  /* temp info on hits found within a CYKBandedScan_jd() call 
   * these are copied to the 'master' hit* structures, free'd
   * and reallocated as nec
   */
  int    tmp_nhits;		 /* number of hits in a seq */
  int   *tmp_hitr;		 /* initial states for hits */
  int   *tmp_hiti;               /* start positions of hits */
  int   *tmp_hitj;               /* end positions of hits */
  float *tmp_hitsc;		 /* scores of hits */

  int    W;   		        /* maximum len of hit; scanning window size */
  int    do_revcomp;		/* true to do reverse complement too */
  int    do_local;		/* TRUE to do local alignment */
  int    do_align;              /* TRUE to calculate and show alignments */
  int    do_dumptrees;		/* TRUE to dump parse trees */
  int    do_qdb=0;		/* TRUE to do a priori banded CYK */
  int    read_qdb;		/* TRUE to read QDBs from a cmbuild outputted file */
  char  *qdb_file;              /* file to read QDBs from (output from cmbuild) */
  FILE  *qdb_fp;
  int    do_projectx;           /* TRUE to activate special in-progress testing code */
  int    do_bdump;              /* TRUE to print out bands */
  /*EPN 08.18.05*/
  int    set_window;            /* TRUE to set window length due to -W option*/

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */
  int      safe_windowlen;	/* initial windowlen (W) used for calc'ing bands */

  /* CM Plan 9 HMM data structures */
  struct cplan9_s       *cp9_hmm;       /* constructed CP9 HMM; written to hmmfile              */
  CP9Map_t              *cp9map;        /* maps the hmm to the cm and vice versa                */
  /*struct cp9_dpmatrix_s *cp9_mx;*/        /* growable DP matrix for viterbi                       */
  struct cp9_dpmatrix_s *cp9_fwd;       /* growable DP matrix for forward                       */
  struct cp9_dpmatrix_s *cp9_bck;       /* growable DP matrix for backward                      */
  struct cp9_dpmatrix_s *cp9_posterior; /* growable DP matrix for posterior decode              */
  float                  swentry;	/* S/W aggregate entry probability       */
  float                  swexit;        /* S/W aggregate exit probability        */
  float                  fwd_sc;        /* score for Forward() run */
  float                  fb_sc;         /* score from Forward or Backward */
  float                  sc;            /* score from CYK */

  /* HMMERNAL!: hmm banded alignment data structures */
  CP9Bands_t *cp9b;             /* data structure for hmm bands (bands on the hmm states) 
				 * and arrays for CM state bands, derived from HMM bands */
  int         do_hbanded;       /* TRUE to first scan with a CP9 HMM to derive bands for a CYK scan */
  double      hbandp;           /* tail loss probability for hmm bands */
  int         use_sums;         /* TRUE to fill and use the posterior sums, false not to. */

  int    v;             /* counter over states of the CM */
  int    x;
  int    debug_level;   /* verbosity level for debugging printf() statements,
			 * passed to many functions. */
  float hmm_thresh;     /* bit score threshold for reporting hits to HMM */

  int do_inside=0;      /* TRUE to use scanning Inside algorithm instead of CYK */
  int alloc_nhits;	/* used to grow the hit arrays */
  int do_scan2hbands;   /* TRUE to use scanning Forward and Backward algs instead of traditional
			 * FB algs to get bands on seqs surviving the filter */
  int   do_filter;              /* TRUE to scan with a CM Plan 9 HMM */
  int   filter_fb;              /* TRUE to use Forward to get start points and Backward
				 * to get end points of promising subsequences*/
  int   filter_weinberg;        /* TRUE to use Forward to get start points and subtract W
				 * to get start points of promising subsequences*/
  int   do_hmmonly;             /* TRUE to scan with a CM Plan 9 HMM ONLY!*/
  int   hmm_pad;                /* number of residues to add to and subtract from end and 
				 * start points of HMM hits prior to CM reevauation, 
				 * respectively. */
  float filter_fraction;        /* fraction of sequence not included in any HMM hit. 
				 * Roughly the fraction of the database filtered out,
				 * but doens't check for overlap: could overcount some bases.
				 */
  int   do_null2;		/* TRUE to adjust scores with null model #2 */
  int   do_zero_inserts;        /* TRUE to zero insert emission scores */

  /* The enforce option (--enfstart and --enfseq), added specifically for enforcing the template 
   * region for telomerase RNA searches */
  int   do_enforce;             /* TRUE to read .enforce file and enforce MATL stretch */
  int   enf_start;              /* if (do_enforce), first MATL node to enforce each parse enter */
  int   enf_end;                /* if (do_enforce), last  MATL node to enforce each parse enter */
  char *enf_seq;                /* if (do_enforce), the subsequence to enforce in nodes from enf_start to 
				 * enf_end */
  int   nd;

  /* E-value statistics (ported from rsearch-1.1) */
  int sample_length= 0;         /* length of samples to use for calc'ing stats (2*W) */
  int num_samples;              /* number of samples to use to calculate EVDs */
  int cutoff_type;              /* either E_CUTOFF for e-values or SCORE_CUTOFF for bit scores */
  float sc_cutoff;              /* bit score cutoff, min bit score to report */
  float e_cutoff;               /* E-value cutoff, min E value to report */
  float cutoff;                 /* either E-value cutoff or bit score cutoff */
  float e_value;                
  int do_stats;                 /* TRUE to calculate E-value statistics */
  double lambda[GC_SEGMENTS];
  double K[GC_SEGMENTS];
  double mu[GC_SEGMENTS];       /* Set from lambda, K, N */
  int seed;                     /* Random seed */
  float score_boost = 0.0;      /* add this value to CYK scores to try to return negative scoring
				 * hits above -1 * score_boost (experimental) */
  long N;                        /* Effective number of sequences for this search */
  float W_scale = 2.0;           /* scale we'll multiply W by to get sample_length */
  int defined_N = FALSE;

  int do_partitions;            /* TRUE if --partition enabled */
  int *partitions;              /* What partition each percentage point goes to */
  int num_partitions = 1;
  int *gc_ct;                   /* gc_ct[x] observed 100-nt segs in DB with GC% of x [0..100] */
  int gc_comp; 
#ifdef USE_MPI
  int mpi_my_rank;              /* My rank in MPI */
  int mpi_num_procs;            /* Total number of processes */
  int mpi_master_rank;          /* Rank of master process */

  /* Initailize MPI, get values for rank and num procs */
  MPI_Init (&argc, &argv);
  
  atexit (exit_from_mpi);
  in_mpi = 1;                /* Flag for exit_from_mpi() */

  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_num_procs);

  /*
   * Determine master process.  This is the lowest ranking one that can do I/O
   */
  mpi_master_rank = get_master_rank (MPI_COMM_WORLD, mpi_my_rank);

  /*printf("B CS rank: %4d master: %4d num: %4d\n", mpi_my_rank, mpi_master_rank, mpi_num_procs);*/

  /* If I'm the master, do the following set up code -- parse arguments, read
     in matrix and query, build model */
  if (mpi_my_rank == mpi_master_rank) {
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = SQFILE_UNKNOWN;
  set_window        = FALSE;
  do_revcomp        = TRUE;
  do_local          = FALSE;
  do_align          = TRUE;
  do_dumptrees      = FALSE;
  do_qdb            = TRUE;      /* QDB is default */
  read_qdb          = FALSE;
  qdb_file          = NULL;
  qdb_fp            = NULL;
  beta              = 0.0000001;
  do_projectx       = FALSE;
  do_bdump          = FALSE;
  do_hmmonly        = FALSE;
  do_filter         = FALSE;
  filter_fb         = FALSE;
  filter_weinberg   = FALSE;
  do_inside         = FALSE;
  hmm_thresh        = 0.;
  do_hbanded        = FALSE;
  hbandp            = 0.0001;
  use_sums          = FALSE;
  do_scan2hbands    = FALSE;
  hmm_pad           = 0;
  do_null2          = FALSE;
  do_zero_inserts   = TRUE;
  sample_length     = 0;
  sc_cutoff         = 0;
  e_cutoff          = 10;
  cutoff_type       = DEFAULT_CUTOFF_TYPE;
  do_stats          = FALSE;
  score_boost       = 0.;
  debug_level       = 0;
  do_partitions     = FALSE;
  num_samples       = 0;
  do_enforce        = FALSE;
  enf_start         = 0;
  enf_end           = 0;
  enf_seq           = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if       (strcmp(optname, "-W")          == 0) 
      { 
	W    = atoi(optarg); set_window = TRUE; 
	if(W < 2) Die("-W <f>, W must be at least 2.\n");
      }
    else if       (strcmp(optname, "-s")          == 0) 
      { 
	sc_cutoff    = atof(optarg); 
	cutoff_type  = SCORE_CUTOFF;
      }
    else if       (strcmp(optname, "-E")          == 0) 
      { 
	e_cutoff = atof(optarg); 
	cutoff_type  = E_CUTOFF;
	do_stats = TRUE;
	num_samples = 1000;
      }
    else if       (strcmp(optname, "-n")          == 0) 
      { 
	num_samples = atoi(optarg);
      }
    else if  (strcmp(optname, "--dumptrees") == 0) do_dumptrees = TRUE;
    else if  (strcmp(optname, "--dumptrees") == 0) do_dumptrees = TRUE;
    else if  (strcmp(optname, "--local")     == 0) do_local     = TRUE;
    else if  (strcmp(optname, "--noalign")   == 0) do_align     = FALSE;
    else if  (strcmp(optname, "--toponly")   == 0) do_revcomp   = FALSE;
    else if  (strcmp(optname, "--inside")    == 0) do_inside    = TRUE;
    else if  (strcmp(optname, "--null2")     == 0) do_null2     = TRUE;
    else if  (strcmp(optname, "--learninserts")== 0) do_zero_inserts = FALSE;
    else if  (strcmp(optname, "--negscore")   == 0) score_boost = -1. * atof(optarg);
    else if  (strcmp(optname, "--enfstart")  == 0) { do_enforce = TRUE; enf_start = atoi(optarg); }
    else if  (strcmp(optname, "--enfseq")    == 0) { do_enforce = TRUE; enf_seq = optarg; } 
    else if  (strcmp(optname, "--hmmfb")   == 0)   { do_filter = TRUE; filter_fb  = TRUE; }
    else if  (strcmp(optname, "--hmmweinberg")   == 0)   
      {
	do_filter = TRUE; filter_weinberg  = TRUE;
	printf("--hmmweinberg not yet supported.\n"); exit(1);
      }
    else if  (strcmp(optname, "--hmmpad")    == 0) { hmm_pad = atoi(optarg); }
    else if  (strcmp(optname, "--hmmonly")   == 0) { do_hmmonly = TRUE; do_align = FALSE; } 
    else if  (strcmp(optname, "--hthresh")   == 0) hmm_thresh   = atof(optarg);
    else if  (strcmp(optname, "--beta")   == 0) beta      = atof(optarg);
    else if  (strcmp(optname, "--noqdb")  == 0) do_qdb    = FALSE;
    else if  (strcmp(optname, "--qdbfile")== 0) { read_qdb  = TRUE; qdb_file = optarg; }
    else if  (strcmp(optname, "--hbanded")   == 0) do_hbanded   = TRUE; 
    else if  (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if  (strcmp(optname, "--banddump")  == 0) do_bdump     = TRUE;
    else if  (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if  (strcmp(optname, "--scan2hbands")== 0) do_scan2hbands= TRUE;
    else if (strcmp (optname, "--partition") == 0) 
      {
	do_partitions = TRUE;
	if (!(set_partitions (&partitions, &num_partitions, optarg)))
	  Die("Specify partitions separated by commas, no spaces, range 1-100, integers only\n");
      }
    else if  (strcmp(optname, "--informat")  == 0) {
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

  if(do_bdump && !do_qdb)
    Die("The --banddump option is incompatible with the --noqdb option.\n");
  if(num_samples != 0 && !do_stats)
    Die("The -n option only makes sense with -E also.\n");
  if(do_enforce && enf_seq == NULL)
    Die("--enfstart only makes sense with --enfseq also.\n");
  if(do_enforce && enf_start == 0)
    Die("--enfseq only makes sense with --enfstart (which can't be 0) also.\n");

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  seqfile = argv[optind++]; 

  /**********************************************
   * Seed random number generator
   **********************************************/
  /* Seed the random number generator with the time */
  /* This is the seed used in the HMMER code */
  seed = (time ((time_t *) NULL));   
  /*seed = 33;*/
  sre_srandom (seed);
  printf ("Random seed: %d\n", seed);

  printf ("W scale of %.1f\n", W_scale);

  /* Initialize partition array, only if we haven't already in set_partitions,
   * in this case we don't use partitions, a single Gumbel is used.
   */
  if(!(do_partitions))
    {
      partitions = MallocOrDie(sizeof(int) * GC_SEGMENTS+1);
      for (i = 0; i < GC_SEGMENTS; i++) 
	partitions[i] = 0;
    }
  /*********************************************** 
   * Preliminaries: open our files for i/o; get a CM
   ***********************************************/

  watch = StopwatchCreate();

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);

  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);

  if(do_qdb && do_hbanded) 
    Die("Can't do --qdb and --hbanded. Pick one.\n");
  if (do_scan2hbands && !(do_hbanded))
    Die("Can't pick --scan2hbands without --hbanded option.\n");
  if (do_hbanded && !(do_filter))
    Die("Can't pick --hbanded without --hmmfb or --hmmweinberg filtering option.\n");
  if (read_qdb && !(do_qdb))
    Die("--qdbfile and --noqdb don't make sense together.\n");
  if (score_boost < 0)
    Die("for --negscore <f>, <f> must be negative.\n");
  if (score_boost != 0. && do_stats)
    Die("--negscore and -E combination not yet implemented.\n");

  /* Open the sequence (db) file */
  status = esl_sqfile_Open(seqfile, format, NULL, &dbfp);
  if (status == eslENOTFOUND) esl_fatal("No such file."); 
  else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
  else if (status == eslEINVAL) esl_fatal("Can’t autodetect stdin or .gz."); 
  else if (status != eslOK) esl_fatal("Failed to open sequence database file, code %d.", status); 

  /* EPN 08.18.05 */
  if (! (set_window)) W = cm->W;
  /*printf("***cm->W : %d***\n", cm->W);*/

  if (hmm_pad >= (W/2)) 
    Die("Value for --hmmpad is too high (must be less than W/2=%d).\n", (int) (W/2));

  CMLogoddsify(cm);
  if(do_zero_inserts)
    CMHackInsertScores(cm);	/* "TEMPORARY" fix for bad priors */
      
  if(do_enforce)
    {
      enf_end = enf_start + strlen(enf_seq) - 1;
      EnforceSubsequence(cm, enf_start, enf_seq);
    }

  cons = CreateCMConsensus(cm, 3.0, 1.0); 
  querylen = cons->clen;

  if (do_filter || do_hmmonly || do_hbanded)
    {
      /* build a CM Plan 9 HMM, and use it to scan. */
      if(!build_cp9_hmm(cm, &cp9_hmm, &cp9map, FALSE, 0.0001, debug_level))
	Die("Couldn't build a CP9 HMM from the CM\n");
      /*debug_print_cp9_params(cp9_hmm); */
    }
  if(do_hbanded)
    cp9b = AllocCP9Bands(cm, cp9_hmm);
  
  /* Relocated ConfigLocal() call to here, AFTER the CM Plan 9 construction.
   * Otherwise its impossible to make a CM Plan 9 HMM from the local CM
   * that passes the current tests to ensure the HMM is "close enough" to
   * the CM. 
   */
  if (do_local)
    { 
      if(do_enforce)
	{
	  ConfigLocalEnforce(cm, 0.5, 0.5, enf_start, enf_end);
	  CMLogoddsify(cm);
	}
      else
	{
	  ConfigLocal(cm, 0.5, 0.5);
	  CMLogoddsify(cm);
	}
      if(do_zero_inserts)
	CMHackInsertScores(cm);	/* "TEMPORARY" fix for bad priors */
      if(do_filter || do_hmmonly || do_hbanded)
	{
	  printf("configuring the CM plan 9 HMM for local alignment.\n");
	  swentry           = 0.5;
	  swexit            = 0.5;
	  CPlan9SWConfig(cp9_hmm, swentry, swexit);
	  CP9Logoddsify(cp9_hmm);
	}
    }

  if (do_qdb || do_bdump)
    {
      /*StopwatchZero(watch);
	StopwatchStart(watch);*/
      if(read_qdb)
	{
	  /* read the bands from a file */
	  if ((qdb_fp = fopen(qdb_file, "r")) == NULL)
	    Die("failed to open QDB file %s", qdb_file);
	  if(!(QDBFileRead(qdb_fp, cm, &dmin, &dmax)))
	    {
	      Die("ERROR reading QDB file: %s.\nDoes it correspond (same number of states) to this model?\n", qdb_file);
	    }
	  fclose(qdb_fp);
	}
      else /* calculate the bands */
	{
	  /* start stopwatch for timing the band calculation */
	  safe_windowlen = W * 2;
	  while(!(BandCalculationEngine(cm, safe_windowlen, beta, 0, &dmin, &dmax, &gamma, do_local)))
	    {
	      /*Die("BandCalculationEngine() failed.\n");*/
	      FreeBandDensities(cm, gamma);
	      free(dmin);
	      free(dmax);
	      safe_windowlen *= 2;
	      /*printf("ERROR BandCalculationEngine returned false, windowlen adjusted to %d\n", safe_windowlen);*/
	    }
	}	  
      /*
       * Change W (which is cm->W as read from the CM file) to dmax[0].
       * dmax[0] could be > cm->W if the beta we're using now is < 1E-7 
       * (1E-7 is the beta value used to determine cm->W in cmbuild). 
       * If W was set at the command line, that value is always used.
       */
      if(!(set_window))
	W = dmax[0];
      if(do_bdump) 
	{
	  printf("beta:%f\n", beta);
	  debug_print_bands(cm, dmin, dmax);
	  PrintDPCellsSaved(cm, dmin, dmax, W);
	}
      /*StopwatchStop(watch);
	StopwatchDisplay(stdout, "\nCPU time (band calc): ", watch);*/
    }

  /* EPN 11.18.05 Now that we know what W is, we need to ensure that
   * cm->el_selfsc * W >= IMPOSSIBLE (cm->el_selfsc is the score for an EL self transition)
   * because we will potentially multiply cm->el_selfsc * W, and add that to 
   * 2 * IMPOSSIBLE, and IMPOSSIBLE must be > -FLT_MAX/3 so we can add it together 3 
   * times (see structs.h). 
   */
  if((cm->el_selfsc * W) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (W+1) * 3);

#ifdef USE_MPI
  }   /* End of first block that is only done by master process */
  /* Barrier for debugging */
  MPI_Barrier(MPI_COMM_WORLD);

  printf("do_qdb: %d\n", do_qdb);
  /* Here we need to broadcast the following parameters:
     num_samples, W, W_scale, cm and dmin and dmax if do_qdb == TRUE*/
    first_broadcast(&num_samples, &W, &W_scale, &cm,  
		    &do_qdb, &dmin, &dmax, &do_inside, mpi_my_rank, mpi_master_rank);
    
#endif
  /**************************************************
   * Make the histogram
   *************************************************/
  
#ifdef USE_MPI
  if (mpi_my_rank == mpi_master_rank) 
    {
#endif
      GetDBInfo(dbfp, &N, &gc_ct);
      if (do_revcomp) N*=2;
#ifdef USE_MPI
    }
#endif
  /* Set sample_length to 2*W if not yet set */
  if (sample_length == 0) sample_length = W_scale * W;
  /*printf("0 W: %d W_scale: %f sample_length: %d\n", W, W_scale, sample_length);*/
  if (num_samples > 0) /* num_samples will be 0 unless do_stats is TRUE */
    {
#ifdef USE_MPI
    /*StopwatchZero(watch);
      StopwatchStart(watch);*/
    if (mpi_num_procs > 1)
	parallel_make_histogram(gc_ct, partitions, num_partitions,
				cm, W, num_samples, sample_length, lambda, K, 
				dmin, dmax, do_inside, mpi_my_rank, mpi_num_procs, 
				mpi_master_rank);
    else 
#endif
      if(do_hmmonly)
	serial_make_histogram (gc_ct, partitions, num_partitions,
			       cm, W, 
			       num_samples, sample_length, lambda, K, 
			       dmin, dmax, cp9_hmm, do_inside, TRUE);
      else
	serial_make_histogram (gc_ct, partitions, num_partitions,
			       cm, W, 
			       num_samples, sample_length, lambda, K, 
			       dmin, dmax, NULL, do_inside, TRUE);

    /*StopwatchStop(watch);
      StopwatchDisplay(stdout, "\nCPU time (histogram): ", watch);*/
    
#ifdef USE_MPI
    if (mpi_my_rank == mpi_master_rank) {
#endif
      
      /* Set mu from K, lambda, N */
      if (num_samples > 0)
	for (i=0; i<GC_SEGMENTS; i++) {
	  mu[i] = log(K[i]*N)/lambda[i];
	}
      
      debug_print_stats(partitions, num_partitions, lambda, mu);

      if (cutoff_type == E_CUTOFF) {
	if (num_samples < 1) 
	  Die ("Cannot use -E option without defined lambda and K\n");
      }
      
      if (num_samples > 0) {
	printf ("Statistics calculated with simulation of %d samples of length %d\n", num_samples, sample_length);
	if (num_partitions == 1) {
	  printf ("No partition points\n");
	} else {
	    printf ("Partition points are: ");
	    for (i=0; i<=GC_SEGMENTS; i++) {
	      if (partitions[i] != partitions[i-1]) {
		printf ("%d ", i);
	      }
	    }
	}
	for (i=0; i<GC_SEGMENTS; i++) {
	  /*printf ("GC = %d\tlambda = %.4f\tmu = %.4f\n", i, lambda[i], mu[i]);*/
	}
	printf ("N = %ld\n", N);
	if (cutoff_type == SCORE_CUTOFF) {
	  printf ("Using score cutoff of %.2f\n", sc_cutoff);
	} else {
	  printf ("Using E cutoff of %.2f\n", e_cutoff);
	}
      } else {
	printf ("lambda and K undefined -- no statistics\n");
	printf ("Using score cutoff of %.2f\n", sc_cutoff);
      }
      fflush(stdout);
    }
  /*******************************************************************************
   * End of making the histogram.
   *******************************************************************************/
  /*************************************************
   *    Do the search
   *************************************************/
#ifdef USE_MPI
  }                   /* Done with second master-only block */
  
  /* Now I need to broadcast the following parameters:
     cutoff, cutoff_type, do_revcomp, do_align, mu, lambda, K, N */
  second_broadcast(&sc_cutoff, &e_cutoff, &cutoff_type, &do_revcomp, &do_align, mu, lambda, 
		   K, &N, mpi_my_rank, mpi_master_rank);
#endif
  if(cutoff_type == E_CUTOFF) cutoff = e_cutoff;
  else cutoff = sc_cutoff;

#ifdef USE_MPI
  if (mpi_num_procs > 1)
    parallel_search_database (dbfp, cm, cons, W, cutoff_type, cutoff, do_revcomp, do_align,
			      do_stats, mu, lambda, dmin, dmax, do_inside,
			      mpi_my_rank, mpi_master_rank, mpi_num_procs);
  else
#endif
    serial_search_database (dbfp, cm, cons, W, cutoff_type, cutoff, do_revcomp, do_align,
			    do_stats, mu, lambda, dmin, dmax, do_inside);
  
  FreeCM(cm);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  in_mpi = 0;
  if (mpi_my_rank == mpi_master_rank) {
#endif
    printf ("Fin\n");
#ifdef USE_MPI
  }
#endif
  return EXIT_SUCCESS;
}

/* OLD CMSEARCH CODE TEMPORARILY LEFT FOR REFERENCE */
#if 0
  /* start stopwatch for timing the search */
  StopwatchZero(watch);
  StopwatchStart(watch);

  maxlen   = 0;
  reversed = FALSE;
  sq = esl_sq_Create(); 
  while (reversed || esl_sqio_Read(dbfp, sq) == eslOK)
    {
      if (sq->n == 0) continue; 	/* silently skip len 0 seqs */
      if (sq->n > maxlen) maxlen = sq->n;
      /* Following line will be unnecessary once ESL_SQ objects have dsq's implemented
       * (i.e. allocated and filled within a esl_sqio_Read() call) once that 
       * happens we can set dsq = sq->dsq or just just sq->dsq instead of dsq.*/      
      dsq = DigitizeSequence(sq->seq, sq->n);

      if (do_filter || do_hmmonly || do_hbanded) 
	/* either scan only with CP9 HMM, or use it to infer bands for CYK.
	 * Information on hits found with the HMM are in hmm_hit* arrays. */
	{
	  hmm_nhits = 0;

	  if(do_hmmonly) /* scan only with CP9 HMM */
	    {
	      fwd_sc = CP9ForwardScan(dsq, 1, sq->n, W, cp9_hmm, &cp9_fwd, 
				      &hmm_nhits, &hmm_hitr, &hmm_hiti, &hmm_hitj, &hmm_hitsc, 
				      hmm_thresh);
	      /*printf("forward  sc: %f\n", fwd_sc);*/
	      /*printf("hmmer_Alphabet: %s\n", hmmer_Alphabet);
		sc = CP9Forward(dsq, sq->n, cp9_hmm, &cp9_fwd);*/
	      /* We're only using the HMM */
	      nhits = hmm_nhits;
	      hitr  = hmm_hitr;
	      hiti  = hmm_hiti;
	      hitj  = hmm_hitj;
	      hitsc = hmm_hitsc;
	    }
	  else /* use CP9 HMM to infer bands for CYK */
	    {
	      /* Scan the sequence with "scanning" Forward and Backward algorithms,
	       * Forward gives likely endpoints (j's), Backward gives likely start points
	       * (i's). Then combine them to get likely i and j pairs (see 
	       * code in CP9_scan.c.)
	       */
	      fb_sc = CP9ForwardBackwardScan(dsq, 1, sq->n, W, cp9_hmm, &cp9_fwd, &cp9_bck,
					     &hmm_nhits, &hmm_hitr, &hmm_hiti, &hmm_hitj, &hmm_hitsc, hmm_thresh,
					     hmm_pad);
	      /*printf("forward/backward sc: %f\n", fb_sc);*/
	      FreeCPlan9Matrix(cp9_fwd);
	      FreeCPlan9Matrix(cp9_bck);

	      /* Print HMM hits */
	      filter_fraction = 0.;
	      for (i = 0; i < hmm_nhits; i++)
		{
	      	  printf("HMM hit %-4d: %6d %6d %8.2f bits\n", i, 
			 reversed ? sq->n - hmm_hiti[i] + 1 : hmm_hiti[i], 
			 reversed ? sq->n - hmm_hitj[i] + 1 : hmm_hitj[i],
			 hmm_hitsc[i]);
		  filter_fraction += hmm_hitj[i] - hmm_hiti[i] + 1;
		}
	      filter_fraction /= sq->n;
	      printf("Fraction removed by filter is about: %5.2f (%6.2f speed-up)\n", (1.-filter_fraction), (1.0 / (filter_fraction)));

	      /* For each hit defined by an i-j pair, use non-scanning, traditional
	       * Forward and Backward algorithms to get hmm bands, we'll use to 
	       * get j and d bands (a la HMM banded cmalign).
	       * Step 1: Get HMM posteriors.
	       * Step 2: posteriors -> HMM bands.
	       * Step 3: HMM bands  ->  CM bands.
	       */
	      nhits = 0; /* number of CM hits is set to 0, but we'll check each HMM hit with the CM,
			  * and this number may grow */
	      alloc_nhits = 10;
	      hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
	      hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
	      hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
	      hitsc = MallocOrDie(sizeof(float) * alloc_nhits);
	      for (i = 0; i < hmm_nhits; i++)
		{
		  if(!do_hbanded) 
		    {
		      if (do_qdb)
			if (do_inside)
			  InsideBandedScan(cm, dsq, dmin, dmax, hmm_hiti[i], hmm_hitj[i], W, 
					   &tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, sc_cutoff);
			else
			  CYKBandedScan(cm, dsq, dmin, dmax, hmm_hiti[i], hmm_hitj[i], W, 
					&tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, sc_cutoff);
		      else if (do_inside)
			InsideScan(cm, dsq, hmm_hiti[i], hmm_hitj[i], W, 
				   &tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, sc_cutoff);
		      else
			CYKScan(cm, dsq, hmm_hiti[i], hmm_hitj[i], W, 
				&tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, sc_cutoff);
		      for (x = 0; x < tmp_nhits; x++)
			{
			  /* copy the new hit info the 'master' data structures */
			  /* Inefficient and ugly but other strategies were confounded by mysterious
			   * memory errors.
			   */
			  hitr[nhits] = tmp_hitr[x];
			  hiti[nhits] = tmp_hiti[x];
			  hitj[nhits] = tmp_hitj[x];
			  hitsc[nhits] = tmp_hitsc[x];
			  nhits++;
			  if (nhits == alloc_nhits) {
			    hitr  = ReallocOrDie(hitr,  sizeof(int)   * (alloc_nhits + 10));
			    hitj  = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits + 10));
			    hiti  = ReallocOrDie(hiti,  sizeof(int)   * (alloc_nhits + 10));
			    hitsc = ReallocOrDie(hitsc, sizeof(float) * (alloc_nhits + 10));
			    alloc_nhits += 10;
			  }
			}
		      free(tmp_hitr);
		      free(tmp_hiti);
		      free(tmp_hitj);
		      free(tmp_hitsc);
		    }
		  else /* do_hbanded */
		    {
		      /*********************************************************/
		      /* TO DO: write function that encapsulates this block... */
		      
		      /* Default behavior: use traditional Forward and Backward algs */
		      if(!(do_scan2hbands))
			{
			  /* Step 1: Get HMM posteriors.*/
			  fb_sc = CP9Forward(dsq, hmm_hiti[i], hmm_hitj[i], cp9_hmm, &cp9_fwd);
			  /*printf("hit: %d i: %d j: %d forward_sc : %f\n", i, hiti[i], hitj[i], fb_sc);*/
			  fb_sc = CP9Backward(dsq, hmm_hiti[i], hmm_hitj[i], cp9_hmm, &cp9_bck);
			  /*printf("CP9 i: %d | backward_sc: %f\n", i, fb_sc);*/
			}
		      else /* scan the subsequence to get the bands (VERY EXPERIMENTAL) */
			{
			  /* Step 1: Get HMM posteriors.*/
			  fb_sc = CP9ForwardScan(dsq, hmm_hiti[i], hmm_hitj[i], W, cp9_hmm, &cp9_fwd,
						 &tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, hmm_thresh);
			  /*printf("SCAN F hit: %d i: %d j: %d forward_sc : %f\n", i, hmm_hiti[i], hmm_hitj[i], fb_sc);*/
			  free(tmp_hitr);
			  free(tmp_hiti);
			  free(tmp_hitj);
			  free(tmp_hitsc);
			  fb_sc = CP9BackwardScan(dsq, hmm_hiti[i], hmm_hitj[i], W, cp9_hmm, &cp9_bck,
						  &tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, hmm_thresh);
			  /*printf("SCAN B hit: %d i: %d j: %d bckward_sc : %f\n", i, hmm_hiti[i], hmm_hitj[i], fb_sc);*/
			  free(tmp_hitr);
			  free(tmp_hiti);
			  free(tmp_hitj);
			  free(tmp_hitsc);
			}
		      /*debug_check_CP9_FB(cp9_fwd, cp9_bck, cp9_hmm, fb_sc, hiti[i], hitj[i], dsq);*/
		      cp9_posterior = cp9_bck;
		      CP9FullPosterior(dsq, hmm_hiti[i], hmm_hitj[i], cp9_hmm, cp9_fwd, cp9_bck, cp9_posterior);
		      /* Step 2: posteriors -> HMM bands.
		       * NOTE: HMM bands have offset sequence indices, from 1 to W, 
		       * with W = hitj[i] - hiti[i] + 1.
		       */
		      if(use_sums)
			CP9_ifill_post_sums(cp9_posterior, hiti[i], hitj[i], cp9b->hmm_M,
					    cp9b->isum_pn_m, cp9b->isum_pn_i, cp9b->isum_pn_d);
		      
		      /* match states */
		      CP9_hmm_band_bounds(cp9_posterior->mmx, hiti[i], hitj[i], cp9b->hmm_M,
					  cp9b->isum_pn_m, cp9b->pn_min_m, cp9b->pn_max_m, 
					  (1.-hbandp), HMMMATCH, use_sums, debug_level);
		      /* insert states */
		      CP9_hmm_band_bounds(cp9_posterior->imx, hiti[i], hitj[i], cp9b->hmm_M,
					  cp9b->isum_pn_i, cp9b->pn_min_i, cp9b->pn_max_i, 
				      (1.-hbandp), HMMINSERT, use_sums, debug_level);
		      /* delete states */
		      CP9_hmm_band_bounds(cp9_posterior->dmx, hiti[i], hitj[i], cp9b->hmm_M,
					  cp9b->isum_pn_d, cp9b->pn_min_d, cp9b->pn_max_d, 
					  (1.-hbandp), HMMDELETE, use_sums, debug_level);

		      if(debug_level != 0)
			{
			  printf("printing hmm bands\n");
			  print_hmm_bands(stdout, sq->n, cp9b->hmm_M, cp9b->pn_min_m, 
					  cp9b->pn_max_m, cp9b->pn_min_i, cp9b->pn_max_i, 
					  cp9b->pn_min_d, cp9b->pn_max_d, hbandp, debug_level);
			}
		      
		      /* Step 3: HMM bands  ->  CM bands. */
		      hmm2ij_bands(cm, cp9map, hmm_hiti[i], hmm_hitj[i], cp9b->pn_min_m, cp9b->pn_max_m, 
				   cp9b->pn_min_i, cp9b->pn_max_i, cp9b->pn_min_d, cp9b->pn_max_d, 
				   cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, debug_level);
		      /* we're going to use bands to search, so we relax ROOT_S bands,
		       * otherwise search is forced to return optimal hit from hmm_hiti[i]
		       * to hmm_hitj[i] b/c imin[0]=imax[0]=hmm_hiti[i] and 
		       * jmin[0]=jmax[0]=hmm_hitj[i].
		       */
		      relax_root_bands(cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax);
		      
		      /* Use the CM bands on i and j to get bands on d, specific to j. */
		      for(v = 0; v < cm->M; v++)
			{
			  cp9b->hdmin[v] = malloc(sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
			  cp9b->hdmax[v] = malloc(sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
			}
		      ij2d_bands(cm, (hmm_hitj[i] - hmm_hiti[i] + 1), cp9b->imin, cp9b->imax, 
				 cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, -1);
		      
		      /*if(debug_level != 0)*/
		      /*PrintDPCellsSaved_jd(cm, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->
			hdmax, (hmm_hitj[i] - hmm_hiti[i] + 1));*/
		      
		      /*debug_print_hd_bands(cm, cp9b->hdmin, cp9b->hdmax, cp9b->jmin, cp9b->jmax);*/
		      
		      FreeCPlan9Matrix(cp9_fwd);
		      FreeCPlan9Matrix(cp9_bck);
		      
		      /* End of block to be encapsulated in function.          */
		      /*********************************************************/
		      
		      if(!do_inside)
			{
			  /* Scan the sequence with a scanning CYK constrained by the HMM bands (experimental) */
			  /*printf("CYK banded jd scanning hit: %d | i: %d | j: %d\n", i, hmm_hiti[i], hmm_hitj[i]);*/
			  /*debug_print_hd_bands(cm, hdmin, hdmax, jmin, jmax);*/
			  CYKBandedScan_jd(cm, dsq, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
					   hmm_hiti[i], hmm_hitj[i], W, 
					   &tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, sc_cutoff);
			}
		      else /* do_inside */
			{			  
			  /* Scan the sequence with a scanning Inside constrained by the HMM bands (experimental) */
			  /*printf("Inside banded jd scanning hit: %d | i: %d | j: %d\n", i, hmm_hiti[i], hmm_hitj[i]);*/
			  InsideBandedScan_jd(cm, dsq, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
					      hmm_hiti[i], hmm_hitj[i], W, 
					      &tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, sc_cutoff);
			}
		      /* copy the new hit info the 'master' data structures */
		      /* Inefficient and ugly but other strategies were confounded by mysterious
		       * memory errors.
		       */
		      for (x = 0; x < tmp_nhits; x++)
			{
			  hitr[nhits] = tmp_hitr[x];
			  hiti[nhits] = tmp_hiti[x];
			  hitj[nhits] = tmp_hitj[x];
			  hitsc[nhits] = tmp_hitsc[x];
			  nhits++;
			  if (nhits == alloc_nhits) {
			    hitr  = ReallocOrDie(hitr,  sizeof(int)   * (alloc_nhits + 10));
			    hitj  = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits + 10));
			    hiti  = ReallocOrDie(hiti,  sizeof(int)   * (alloc_nhits + 10));
			    hitsc = ReallocOrDie(hitsc, sizeof(float) * (alloc_nhits + 10));
			    alloc_nhits += 10;
			  }
			}
		      free(tmp_hitr);
		      free(tmp_hiti);
		      free(tmp_hitj);
		      free(tmp_hitsc);
		      
		      /* If we're done with them, free hdmin and hdmax, these are allocated
		       * differently for each sequence. */
		      if(!(do_align))
			{		 
			  for(v = 0; v < cm->M; v++)
			    {
			      free(cp9b->hdmin[v]);
			      free(cp9b->hdmax[v]);
			    }
			}
		    }
		}
	    }
	}
      else if (do_qdb)
	if (do_inside)
	  InsideBandedScan(cm, dsq, dmin, dmax, 1, sq->n, W, 
			   &nhits, &hitr, &hiti, &hitj, &hitsc, sc_cutoff);
	else
	  CYKBandedScan(cm, dsq, dmin, dmax, 1, sq->n, W, 
		      &nhits, &hitr, &hiti, &hitj, &hitsc, sc_cutoff);
      else if (do_inside)
	InsideScan(cm, dsq, 1, sq->n, W, 
		   &nhits, &hitr, &hiti, &hitj, &hitsc, sc_cutoff);
      else
	CYKScan(cm, dsq, 1, sq->n, W, 
		&nhits, &hitr, &hiti, &hitj, &hitsc, sc_cutoff);
      if (! reversed) printf("sequence: %s\n", sq->name);

      ip = 0;
      for (i = 0; i < nhits; i++)
	{
	  if(!do_null2)
	    {
	      if (do_stats) 
		{
		  gc_comp = get_gc_comp (sq->seq, hiti[i], hitj[i]);
		  e_value = RJK_ExtremeValueE(hitsc[i], mu[gc_comp], 
					      lambda[gc_comp]);
		  printf("\tgc: %d sc: %.2f E: %g P: %g\n", gc_comp, hitsc[i], e_value,
			 esl_gumbel_surv((double) hitsc[i], mu[gc_comp], 
					 lambda[gc_comp]));
		  
		  if(e_value <= e_cutoff)
		    {
		      printf("hit %-4d: %6d %6d %8.2f bits   E = %8.3g, P = %8.3g\n", ip, 
			     reversed ? sq->n - hiti[i] + 1 : hiti[i], 
			     reversed ? sq->n - hitj[i] + 1 : hitj[i],
			     hitsc[i], 
			     e_value,
			     esl_gumbel_surv((double) hitsc[i], mu[gc_comp], 
					     lambda[gc_comp]));
		      ip++;
		    }
		}
	      else 
		{
		  printf("hit %-4d: %6d %6d %8.2f bits\n", ip, 
			 reversed ? sq->n - hiti[i] + 1 : hiti[i], 
			 reversed ? sq->n - hitj[i] + 1 : hitj[i],
			 hitsc[i]);
		  ip++;
		}
	    }
	  if (do_null2 || do_align)
	    {
	      /* For the null2 score correction we need a trace, so we have to do 
	       * the alignment.
	       */
	      if(!(do_hbanded))
		{
		  CYKDivideAndConquer(cm, dsq, sq->n, 
				      hitr[i], hiti[i], hitj[i], &tr, 
				      NULL, NULL); /* don't use qd bands */

		  if (do_null2) sc = hitsc[i] - CM_TraceScoreCorrection(cm, tr, dsq);
		  else          sc = hitsc[i];

		  if(do_stats) 
		    {
		      gc_comp = get_gc_comp (sq->seq, hiti[i], hitj[i]);
		      e_value = RJK_ExtremeValueE(hitsc[i], mu[gc_comp], 
						  lambda[gc_comp]);
		    }
		  if ((!do_stats && sc >= sc_cutoff) ||
		      (do_stats  && e_value <= e_cutoff))
		    {
		      printf("hit %-4d: %6d %6d %8.2f bits\n", ip, 
			     reversed ? sq->n - hiti[i] + 1 : hiti[i], 
			     reversed ? sq->n - hitj[i] + 1 : hitj[i],
			     sc);
		      ip++;

		      if (do_align)
			{
			  ali = CreateFancyAli(tr, cm, cons, dsq);
			  PrintFancyAli(stdout, ali);
			  FreeFancyAli(ali);
			}
		    }

		  if (do_dumptrees)
		    if(!do_qdb) 
		      ParsetreeDump(stdout, tr, cm, dsq);
		    else
		      BandedParsetreeDump(stdout, tr, cm, dsq, gamma, W, dmin, dmax);
		  printf("\tscore = %.2f\n\n", ParsetreeScore(cm,tr,dsq, do_null2));
		  
		  FreeParsetree(tr);
		}
	      else /* do_hbanded==TRUE */
		{
		  /* Derive HMM bands, following block is too long... */
		  /*********************************************************/
		  /* TO DO: write function that encapsulates this block... */
		  /* Step 1: Get HMM posteriors.*/
		  fb_sc = CP9Forward(dsq, hiti[i], hitj[i], cp9_hmm, &cp9_fwd);
		  /*printf("hit: %d i: %d j: %d forward_sc : %f\n", i, hiti[i], hitj[i], fb_sc);*/
		  fb_sc = CP9Backward(dsq, hiti[i], hitj[i], cp9_hmm, &cp9_bck);
		  /*printf("CP9 i: %d | backward_sc: %f\n", i, fb_sc);*/

		  /*debug_check_CP9_FB(cp9_fwd, cp9_bck, cp9_hmm, fb_sc, hiti[i], hitj[i], dsq);*/
		  cp9_posterior = cp9_bck;
		  CP9FullPosterior(dsq, hiti[i], hitj[i], cp9_hmm, cp9_fwd, cp9_bck, cp9_posterior);

		  /* Step 2: posteriors -> HMM bands.
		   * NOTE: HMM bands have offset sequence indices, from 1 to W, 
		   * with W = hitj[i] - hiti[i] + 1.
		   */
		  if(use_sums)
		    CP9_ifill_post_sums(cp9_posterior, hiti[i], hitj[i], cp9b->hmm_M,
					cp9b->isum_pn_m, cp9b->isum_pn_i, cp9b->isum_pn_d);

		  /* match states */
		  CP9_hmm_band_bounds(cp9_posterior->mmx, hiti[i], hitj[i], cp9b->hmm_M,
				      cp9b->isum_pn_m, cp9b->pn_min_m, cp9b->pn_max_m, 
				      (1.-hbandp), HMMMATCH, use_sums, debug_level);
		  /* insert states */
		  CP9_hmm_band_bounds(cp9_posterior->imx, hiti[i], hitj[i], cp9b->hmm_M,
				      cp9b->isum_pn_i, cp9b->pn_min_i, cp9b->pn_max_i, 
				      (1.-hbandp), HMMINSERT, use_sums, debug_level);
		  /* delete states */
		  CP9_hmm_band_bounds(cp9_posterior->dmx, hiti[i], hitj[i], cp9b->hmm_M,
				      cp9b->isum_pn_d, cp9b->pn_min_d, cp9b->pn_max_d, 
				      (1.-hbandp), HMMDELETE, use_sums, debug_level);

		  if(debug_level != 0)
		    {
		      printf("printing hmm bands\n");
		      print_hmm_bands(stdout, sq->n, cp9b->hmm_M, cp9b->pn_min_m, 
				      cp9b->pn_max_m, cp9b->pn_min_i, cp9b->pn_max_i, 
				      cp9b->pn_min_d, cp9b->pn_max_d, hbandp, debug_level);
		    }
		  
		  /* Step 3: HMM bands  ->  CM bands. */
		  hmm2ij_bands(cm, cp9map, hiti[i], hitj[i], cp9b->pn_min_m, cp9b->pn_max_m, 
			       cp9b->pn_min_i, cp9b->pn_max_i, cp9b->pn_min_d, cp9b->pn_max_d, 
			       cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, debug_level);
	  
		  /* Use the CM bands on i and j to get bands on d, specific to j. */
		  for(v = 0; v < cm->M; v++)
		    {
		      cp9b->hdmin[v] = malloc(sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
		      cp9b->hdmax[v] = malloc(sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
		    }
		  ij2d_bands(cm, (hitj[i] - hiti[i] + 1), cp9b->imin, cp9b->imax, 
			     cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, -1);

		  /*if(debug_level != 0)*/
		  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (hitj[i] - hiti[i] + 1));*/

		  /*debug_print_hd_bands(cm, hdmin, hdmax, jmin, jmax);*/

		  FreeCPlan9Matrix(cp9_fwd);
		  FreeCPlan9Matrix(cp9_bck);

		  /* End of block to be encapsulated in function.          */
		  /*********************************************************/

		  printf("Aligning subseq from i: %d to j: %d using HMM bands\n", hiti[i], hitj[i]);
		  /* we don't overwrite hitsc[i], the optimal score,
		   * which may be missed by the HMM banded alignment 
		   */
		  sc = CYKInside_b_jd(cm, dsq, sq->n, 0, hiti[i], hitj[i], &tr, 
				      cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
				      cp9b->safe_hdmin, cp9b->safe_hdmax);
		  if(do_null2)
		    {
		      sc -= CM_TraceScoreCorrection(cm, tr, dsq);
		      if(sc >= sc_cutoff) /* only print alignments with
					  correct CYK scores > reporting sc_cutoff */
			{
			  printf("hit %-4d: %6d %6d %8.2f bits\n", ip, 
				 reversed ? sq->n - hiti[i] + 1 : hiti[i], 
				 reversed ? sq->n - hitj[i] + 1 : hitj[i],
				 hitsc[i]);
			  printf("\tCYK i: %5d j: %5d sc: %10.2f bits\n", hiti[i], hitj[i], sc);
			  ip++;
			}
		    }
		  if(sc >= sc_cutoff)
		    {
		      ali = CreateFancyAli(tr, cm, cons, dsq);
		      PrintFancyAli(stdout, ali);
		      FreeFancyAli(ali);
		    }
		  FreeParsetree(tr);

		  /* Free hdmin and hdmax, these are allocated
		   * differently for each sequence. 
		   */
		  for(v = 0; v < cm->M; v++)
		    {
		      free(cp9b->hdmin[v]);
		      free(cp9b->hdmax[v]);
		    }
		  /*if(bdump_level > 0)
		    banded_trace_info_dump(cm, tr[i], safe_hdmin, safe_hdmax, bdump_level);
		  */
		}
	    }
	}
      free(hitr);
      free(hiti);
      free(hitj);
      free(hitsc);
      if((do_hbanded || do_filter) && !do_hmmonly)
	{
	  free(hmm_hitr);
	  free(hmm_hiti);
	  free(hmm_hitj);
	  free(hmm_hitsc);
	}
      free(dsq);
      if (! reversed && do_revcomp) {
	revcomp(sq->seq,sq->seq);
	reversed = TRUE;
      } else {
	reversed = FALSE;
	esl_sq_Reuse(sq);
      }
    }

  StopwatchStop(watch);
  StopwatchDisplay(stdout, "\nCPU time (search)   : ", watch);
  if(do_filter || do_hmmonly) printf("CP9 Forward memory  :   %8.2f MB\n", CP9ForwardScanRequires(cp9_hmm, maxlen, W));
  printf("CYK memory          :   %8.2f MB\n\n", CYKScanRequires(cm, maxlen, W));

  if (do_qdb)
    {
      if(!read_qdb)
	FreeBandDensities(cm, gamma);
      free(dmin);
      free(dmax);
    }
    
  if(do_hbanded)
    FreeCP9Bands(cp9b);
  
  FreeCMConsensus(cons);
  FreeCM(cm);
  esl_sqfile_Close(dbfp); 
  esl_sq_Destroy(sq); 
  StopwatchFree(watch);
  SqdClean();
  return EXIT_SUCCESS;
}
#endif
  
static int  
QDBFileRead(FILE *fp, CM_t *cm, int **ret_dmin, int **ret_dmax)
{
  char   *buf;
  int     n;			/* length of buf */
  char   *s;
  int     M;			/* number of states in model */
  int     v;		        /* counter for states */
  char   *tok;
  int     toklen;
  int    *dmin;
  int    *dmax;
  int     read_v;

  /* format of QDB file: 
   * line  1        :<cm->M>
   * lines 2 -> M+1 :<v> <dmin> <dmax> */

  buf = NULL;
  n   = 0;
  if (feof(fp) || sre_fgets(&buf, &n, fp) == NULL) return 0;

  s   = buf;
  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
  if (! IsInt(tok))                                     goto FAILURE;
  M = atoi(tok);
  if(M != cm->M) goto FAILURE;

  dmin = MallocOrDie(sizeof(int) * cm->M);
  dmax = MallocOrDie(sizeof(int) * cm->M);  
  *ret_dmin = NULL;
  *ret_dmax = NULL;

  v = 0;
  while (sre_fgets(&buf, &n, fp) != NULL) 
    {
      if (strncmp(buf, "//", 2) == 0) 
	break;
      s   = buf;
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      read_v = atoi(tok);
      if(v != read_v) goto FAILURE;

      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      dmin[v] = atoi(tok);

      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      dmax[v] = atoi(tok);

      v++;
    }
  if(v != M) goto FAILURE;

  *ret_dmin = dmin;
  *ret_dmax = dmax;
  
  if (buf != NULL) free(buf);
  return 1;

 FAILURE:
  if (cm != NULL)  FreeCM(cm);
  if (buf != NULL) free(buf);
  return 0;
}

/* Function: set_partitions
 * Date:     RJK, Mon, Oct 7, 2002 [St. Louis]
 *           Infernalification by EPN 
 * Purpose:  Set partitions array given a 'list' from the command line.
 */
int set_partitions (int **ret_partitions, int *num_partitions, char *list) 
{
  int *partitions;   /* What partition each percentage point goes to */
  int *partition_pt; /* partition_pt[i] is TRUE if i is to be a partition point */
  int i;
  int cur_point;
  int cur_partition;

  printf("in set partitions\n");
  partitions   = MallocOrDie(sizeof(int) * GC_SEGMENTS+1);
  partition_pt = MallocOrDie(sizeof(int) * GC_SEGMENTS+1);
  
  /* initialize partition_pt array */
  for (i = 0; i < GC_SEGMENTS; i++) 
    partition_pt[i] = FALSE;

  /* Read the partition points */
  cur_point = 0;
  while (*list != '\0') 
    {
      if (isdigit(*list)) 
	cur_point = (cur_point * 10) + (*list - '0');
      else if (*list == ',') 
	{
	  if (cur_point <= 0 || cur_point >= GC_SEGMENTS)
	    return(0); 
	  else
	    partition_pt[cur_point] = TRUE;
	  cur_point = 0;
	}
      else 
	return(0);
      list++;
    }
  /* keep track of the last partition point, which wasn't followed by a comma */
  if (cur_point <= 0 || cur_point >= GC_SEGMENTS)
    return(0);
  else
    partition_pt[cur_point] = TRUE;
  
  /* Set the partitions */
  cur_partition = 0;
  partitions[0] = 0;
  /* first possible point for 2nd partition is 1 */
  for(i=1; i < GC_SEGMENTS; i++)
    {
      if(partition_pt[i]) 
	cur_partition++;
      partitions[i] = cur_partition;
    }

  *ret_partitions = partitions;
  (*num_partitions) = cur_partition + 1;
  free(partition_pt);

  /*  for(i=0; i < GC_SEGMENTS; i++)
  printf("Partition[%d]: %d\n", i, partitions[i]);*/
  return(1);
}

/* Function: debug_print_stats
 */
int debug_print_stats(int *partitions, int num_partitions, double *lambda, double *mu)
{
  int i;
  int cur_partition;
  float sc;

  printf("in debug_print_stats num_partitions: %d\n", num_partitions);
  cur_partition = 0;
  printf("cur_partition: %d\n", cur_partition);
  for (i=0; i<GC_SEGMENTS; i++) 
    {
      printf("i: %d\n", i);
      if (partitions[i] == cur_partition) 
	{
	  printf("partition i:%d starts at: %d\n", cur_partition, i);
	  for(sc = 0.0; sc < 100.0; sc +=1.)
	    {
	      printf (" DEBUG Score = %.2f, E = %.4g, P = %.4g\n", sc,
		      RJK_ExtremeValueE(sc, mu[i], lambda[i]),
		      esl_gumbel_surv((double) sc, mu[i], lambda[i]));
	      printf("\tmu[%d]: %f lambda[%d]: %f\n", i, mu[i], i, lambda[i]);
	    }
	  printf("\n");
	  cur_partition++;
	} 
    }
  printf("end of debug_print_stats\n");
  return 1;
}


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
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "hmmband.h"
#include "stats.h"
#include "esl_gumbel.h"
#include "esl_sqio.h"
#include "mpifuncs.h"
#include "cm_wrappers.h"

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
   -E <f> : use cutoff E-value of <f> (default ignored; not-calc'ed)\n\
   -n <n> : determine EVD with <n> samples (default with -E: 1000)\n\
   -S <f> : use cutoff bit score of <f> [default: 0]\n\
   -R     : build RSEARCH CM from 1 seq <aln file> (replaces <cmfile>)\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input alignment is in format <s>, not FASTA\n\
   --toponly     : only search the top strand\n\
   --local       : do local alignment\n\
   --noalign     : find start/stop only; don't do alignments\n\
   --window <n>  : set scanning window size to <n> (default: precalc'd in cmbuild)\n\
   --dumptrees   : dump verbose parse tree information for each hit\n\
   --partition <n>[,<n>]... : partition points for different GC content EVDs\n\
   --inside      : scan with Inside, not CYK (caution much slower(!))\n\
   --null2       : turn on the post hoc second null model [df:OFF]\n\
   --learninserts: do not set insert emission scores to 0\n\
   --negsc    <f>: set min bit score to report as <f> < 0 (experimental)\n\
   --enfstart <n>: enforce MATL stretch starting at CM node <n>\n\
   --enfseq <s>  : enforce MATL stretch starting at --enfstart <n> emits seq <s>\n\
\n\
  * Filtering options using a CM plan 9 HMM (*in development*):\n\
   --hmmfb        : use Forward to get end points & Backward to get start points\n\
   --hmmweinberg  : use Forward to get end points, subtract W for start points\n\
   --hmmpad <n>   : subtract/add <n> residues from start/end [df:0]\n\
   --hmmonly      : don't use CM at all, just scan with HMM (Forward + Backward)\n\
   --hmmE <f>     : use cutoff E-value of <f> for CP9 (possibly filtered) scan\n\
   --hmmS <f>     : use cutoff bit score of <f> for CP9 (possibly filtered) scan\n\
   --hmmnegsc <f> : set min bit score to report as <f> < 0 (experimental)\n\
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
  { "-S", TRUE, sqdARG_FLOAT }, 
  { "-E", TRUE, sqdARG_FLOAT }, 
  { "-n", TRUE, sqdARG_INT }, 
  { "-R", TRUE, sqdARG_NONE }, 
  { "--dumptrees",  FALSE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--local",      FALSE, sqdARG_NONE },
  { "--noalign",    FALSE, sqdARG_NONE },
  { "--toponly",    FALSE, sqdARG_NONE },
  { "--window",     FALSE, sqdARG_INT }, 
  { "--inside",     FALSE, sqdARG_NONE },
  { "--null2",      FALSE, sqdARG_NONE },
  { "--learninserts",FALSE, sqdARG_NONE},
  { "--negsc",      FALSE, sqdARG_FLOAT},
  { "--hmmfb",      FALSE, sqdARG_NONE },
  { "--hmmweinberg",FALSE, sqdARG_NONE},
  { "--hmmpad",     FALSE, sqdARG_INT },
  { "--hmmonly",    FALSE, sqdARG_NONE },
  { "--hmmE",       FALSE, sqdARG_FLOAT},
  { "--hmmS",       FALSE, sqdARG_FLOAT},
  { "--hmmnegsc",   FALSE, sqdARG_FLOAT},
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
  char            *cmfile;      /* file to read CM from                     */	
  CMFILE          *cmfp;        /* open CM file for reading                 */
  char            *seqfile;     /* file to read sequences from              */
  int              format;      /* format of sequence file                  */
  int              status;      /* status of the sequence file              */
  char            *optname;     /* name of option found by Getopt()         */
  char            *optarg;      /* argument found by Getopt()               */
  int              optind;      /* index in argv[]                          */
  ESL_SQFILE	  *dbfp;        /* open seqfile for reading                 */
  CM_t            *cm;          /* a covariance model                       */
  Stopwatch_t     *watch;       /* times band calc, then search time        */
  int              i;           /* counter                                  */
  CMConsensus_t   *cons;	/* precalculated consensus info for display */
  double            beta;       /* tail loss prob for query dependent bands */
  int             *preset_dmin = NULL; /* if --qdbfile, dmins read from file*/
  int             *preset_dmax = NULL; /* if --qdbfile, dmaxs read from file*/
  int              do_revcomp;  /* true to do reverse complement also       */
  int              do_local;    /* TRUE to do local alignment               */
  int              do_align;    /* TRUE to calculate and show alignments    */
  int              do_dumptrees;/* TRUE to dump parse trees                 */
  int              do_qdb;      /* TRUE to do query dependent banded CYK    */
  int              read_qdb;    /* TRUE to read QDBs from cmbuild file      */
  char            *qdb_file;    /* cmbuild output file to read QDBs from    */
  FILE            *qdb_fp;      /* the open qdb_file                        */
  int              do_bdump;    /* TRUE to print out bands                  */
  int              set_window;  /* TRUE to set window len b/c of -W option  */
  int              set_W;	/* W set at command line, only works --noqdb*/
  int              seed;        /* Random seed                              */
  float            sc_boost;    /* value added to CYK bit scores, allows    *
				 * hits > (-1 * sc_boost) (EXPERIMENTAL)    */
  float            cp9_sc_boost;/* value added to Forward bit scores, allows*
				 * hits > (-1 * sc_boost) (EXPERIMENTAL)    */

  /* E-value statistics (ported from rsearch-1.1). 
   * Some variables set to a value as they're defined to ease 
   * MPI implementation. 
   */
  int   do_cm_stats;            /* TRUE to calculate E-value stats for CM   */
  int   sample_length= 0;       /* sample len used for calc'ing stats (2*W) */
  int   num_samples;            /* # samples used to calculate EVDs         */
  int   cm_cutoff_type;         /* E_CUTOFF or SCORE_CUTOFF for CM          */
  float cm_sc_cutoff;           /* min CM bit score to report               */
  float cm_e_cutoff;            /* max CM E value to report                 */
  int   cp9_cutoff_type;        /* E_CUTOFF or SCORE_CUTOFF for CP9         */
  float cp9_sc_cutoff;          /* min CP9 bit score to report              */
  float cp9_e_cutoff;           /* max CP9 E value to report                */
  long  N;                      /* effective number of seqs for this search */
  float W_scale = 2.0;          /* W_scale * W= sample_length               */
  int   do_partitions;          /* TRUE if --partition enabled              */
  int  *partitions;             /* partition each GC % point seg goes to    */
  int   num_partitions = 1;     /* number of partitions                     */
  int  *gc_ct;                  /* gc_ct[x] observed 100-nt segs in DB with *
				 * GC% of x [0..100]                        */

  /* The enforce option (--enfstart and --enfseq), added specifically for   *
   * enforcing the template region for telomerase RNA searches              */
  int   do_enforce;             /* TRUE to enforce a MATL stretch is used   */
  int   enf_start;              /* first MATL node to enforce each parse use*/
  int   enf_end;                /* last  MATL node to enforce each parse use*/
  char *enf_seq;                /* the subsequence to enforce emitted by    *
                                 * in nodes from enf_start to enf_end       */

  /* HMMERNAL!: hmm banded alignment data structures */
  int         do_hbanded;  /* TRUE to scan 1st with a CP9 to derive bands for CYK scan */
  double      hbandp;      /* tail loss probability for hmm bands */
  int         use_sums;    /* TRUE to fill and use the posterior sums, false not to. */
  int         debug_level; /* verbosity level for debugging printf() statements,
			    * passed to many functions. */

  int   do_inside;         /* TRUE to use scanning Inside algorithm instead of CYK */
  int   do_scan2hbands;    /* TRUE to use scanning Forward/Backward algs instead 
			    * of traditional FB algs to get bands on seqs surviving 
			    * the filter */
  int   do_hmmfb;          /* TRUE to use Forward to get start points and Backward
			    * to get end points of promising subsequences*/
  int   do_hmmweinberg;    /* TRUE to use Forward to get start points and subtract W
			    * to get start points of promising subsequences*/
  int   do_hmmonly;        /* TRUE to scan with a CM Plan 9 HMM ONLY!*/
  int   do_cp9_stats;      /* TRUE to calculate CP9 stats for HMM */
  int   hmm_pad;           /* number of residues to add to and subtract from end and 
			    * start points of HMM hits prior to CM reevauation, 
			    * respectively. */
  int   do_null2;	   /* TRUE to adjust scores with null model #2 */
  int   do_zero_inserts;   /* TRUE to zero insert emission scores */
  
  /* variables related to the -R rsearch option */
  int      do_rsearch;	   /* TRUE to build RSEARCH CM                    */
  char    *alifile;        /* seqfile to read single seq alignment from   */ 
  int      aliformat;	   /* format of seqfile                       */
  MSAFILE *afp;            /* open alignment file                     */
  MSA     *msa;            /* a multiple sequence alignment           */
  fullmat_t       *fullmat;     /* The full matrix */
  char matrixname[256];         /* Name of the matrix, from -m */
  FILE            *matfp;       /* open matrix file for reading */
  float           ralpha =   DEFAULT_RALPHA; 
  float           rbeta =    DEFAULT_RBETA;
  float           ralphap =  DEFAULT_RALPHAP;
  float           rbetap =   DEFAULT_RBETAP;
  float           rbeginsc = DEFAULT_RBEGINSC;
  float           rendsc =   DEFAULT_RENDSC;
  float           rmin_alpha_beta_sum;
  int             rquerylen;     /* length of the query file */

#ifdef USE_MPI
  int mpi_my_rank;              /* My rank in MPI */
  int mpi_num_procs;            /* Total number of processes */
  int mpi_master_rank;          /* Rank of master process */

  /* Initailize MPI, get values for rank and naum procs */
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
  do_revcomp        = TRUE;
  do_local          = FALSE;
  do_align          = TRUE;
  do_dumptrees      = FALSE;
  do_qdb            = TRUE;      /* QDB is default */
  read_qdb          = FALSE;
  qdb_file          = NULL;
  qdb_fp            = NULL;
  beta              = DEFAULT_BETA;
  do_bdump          = FALSE;
  do_hmmonly        = FALSE;
  set_window        = FALSE;
  do_hmmfb          = FALSE;
  do_hmmweinberg    = FALSE;
  do_cp9_stats      = FALSE;
  do_inside         = FALSE;
  do_hbanded        = FALSE;
  hbandp            = DEFAULT_HBANDP;
  use_sums          = FALSE;
  do_scan2hbands    = FALSE;
  hmm_pad           = 0;
  do_null2          = FALSE;
  do_zero_inserts   = TRUE;
  sample_length     = 0;
  cm_cutoff_type    = DEFAULT_CM_CUTOFF_TYPE;
  cm_sc_cutoff      = DEFAULT_CM_CUTOFF;
  cm_e_cutoff       = DEFAULT_CM_CUTOFF;
  cp9_cutoff_type   = DEFAULT_CP9_CUTOFF_TYPE;
  cp9_sc_cutoff     = DEFAULT_CP9_CUTOFF;
  cp9_e_cutoff      = DEFAULT_CP9_CUTOFF;
  do_cm_stats       = FALSE;
  do_cp9_stats      = FALSE;
  sc_boost          = 0.;
  cp9_sc_boost      = 0.;
  debug_level       = 0;
  do_partitions     = FALSE;
  num_samples       = 0;
  do_enforce        = FALSE;
  enf_start         = 0;
  enf_end           = 0;
  enf_seq           = NULL;
  do_rsearch        = FALSE;
  aliformat         = MSAFILE_UNKNOWN;   /* autodetect by default */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if       (strcmp(optname, "-W")          == 0) 
      { 
	set_W = atoi(optarg); set_window = TRUE; 
	if(set_W < 2) Die("-W <f>, W must be at least 2.\n");
      }
    else if       (strcmp(optname, "-S")          == 0) 
      { 
	cm_sc_cutoff    = atof(optarg); 
	cm_cutoff_type  = SCORE_CUTOFF;
      }
    else if       (strcmp(optname, "-E")          == 0) 
      { 
	cm_e_cutoff = atof(optarg); 
	cm_cutoff_type  = E_CUTOFF;
	do_cm_stats = TRUE;
	num_samples = 1000;
      }
    else if  (strcmp(optname, "-R")          == 0) do_rsearch   = TRUE;
    else if  (strcmp(optname, "-n")          == 0) num_samples  = atoi(optarg);
    else if  (strcmp(optname, "--dumptrees") == 0) do_dumptrees = TRUE;
    else if  (strcmp(optname, "--local")     == 0) do_local     = TRUE;
    else if  (strcmp(optname, "--noalign")   == 0) do_align     = FALSE;
    else if  (strcmp(optname, "--toponly")   == 0) do_revcomp   = FALSE;
    else if  (strcmp(optname, "--inside")    == 0) do_inside    = TRUE;
    else if  (strcmp(optname, "--null2")     == 0) do_null2     = TRUE;
    else if  (strcmp(optname, "--learninserts")== 0) do_zero_inserts = FALSE;
    else if  (strcmp(optname, "--negsc")       == 0) sc_boost = -1. * atof(optarg);
    else if  (strcmp(optname, "--enfstart")    == 0) { do_enforce = TRUE; enf_start = atoi(optarg); }
    else if  (strcmp(optname, "--enfseq")      == 0) { do_enforce = TRUE; enf_seq = optarg; } 
    else if  (strcmp(optname, "--hmmfb")       == 0)   do_hmmfb = TRUE;
    else if  (strcmp(optname, "--hmmweinberg") == 0)   do_hmmweinberg = TRUE;
    else if  (strcmp(optname, "--hmmnegsc")    == 0) cp9_sc_boost = -1. * atof(optarg);
    else if  (strcmp(optname, "--hmmE")        == 0)   
      { 
	cp9_e_cutoff = atof(optarg); 
	cp9_cutoff_type  = E_CUTOFF;
	do_cp9_stats = TRUE;
	num_samples = 1000;
      }
    else if  (strcmp(optname, "--hmmS")        == 0)   
      { 
	cp9_sc_cutoff = atof(optarg); 
	cp9_cutoff_type  = SCORE_CUTOFF;
      }
    else if  (strcmp(optname, "--hmmpad")    == 0) { hmm_pad = atoi(optarg); }
    else if  (strcmp(optname, "--hmmonly")   == 0) { do_hmmonly = TRUE; do_align = FALSE; } 
    else if  (strcmp(optname, "--beta")   == 0) beta      = atof(optarg);
    else if  (strcmp(optname, "--noqdb")  == 0) do_qdb    = FALSE;
    else if  (strcmp(optname, "--qdbfile")== 0) { read_qdb  = TRUE; qdb_file = optarg; }
    else if  (strcmp(optname, "--hbanded")   == 0) do_hbanded   = TRUE; 
    else if  (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if  (strcmp(optname, "--banddump")  == 0) do_bdump     = TRUE;
    else if  (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if  (strcmp(optname, "--scan2hbands")== 0) do_scan2hbands= TRUE;
    else if  (strcmp (optname, "--partition") == 0) 
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
  
  /* Check for incompatible option combos. */
  if(do_cp9_stats && 
     ((!do_hmmonly) && (!do_hmmweinberg)))
    Die("--hmmE only makes sense with --hmmonly or --hmmweinberg.\n");
  if(do_cm_stats && do_hmmonly)
    Die("-E and --hmmonly combo doesn't make sense, did you mean --hmmE and --hmmonly?\n");
  if(do_bdump && !do_qdb)
    Die("The --banddump option is incompatible with the --noqdb option.\n");
  if(num_samples != 0 && (!do_cm_stats && !do_cp9_stats))
    Die("The -n option only makes sense with -E or --hmmE also.\n");
  if(do_enforce && enf_seq == NULL)
    Die("--enfstart only makes sense with --enfseq also.\n");
  if(do_enforce && enf_start == 0)
    Die("--enfseq only makes sense with --enfstart (which can't be 0) also.\n");
  if(do_qdb && do_hbanded) 
    Die("Can't do --qdb and --hbanded. Pick one.\n");
  if (do_scan2hbands && !(do_hbanded))
    Die("Can't pick --scan2hbands without --hbanded option.\n");
  if (do_hbanded && !(do_hmmweinberg || do_hmmfb))
    Die("Can't pick --hbanded without --hmmfb or --hmmweinberg filtering option.\n");
  if (read_qdb && !(do_qdb))
    Die("--qdbfile and --noqdb don't make sense together.\n");
  if (sc_boost < 0)
    Die("for --negsc <f>, <f> must be negative.\n");
  if (cp9_sc_boost < 0)
    Die("for --hmmnegsc <f>, <f> must be negative.\n");
  if(do_bdump && !(do_qdb))
    Die("--banddump and --noqdb combination not supported.\n");
  if(set_window && do_qdb)
    Die("--window only works with --noqdb.\n");
#if USE_MPI
  if(read_qdb && ((mpi_num_procs > 1) && (mpi_my_rank == mpi_master_rank)))
    Die("Sorry, you can't read in bands with --qdbfile in MPI mode.\n");
#endif

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  if (!do_rsearch) cmfile = argv[optind++];
  else alifile = argv[optind++];
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
  /**************************************************
   * Preliminaries: open our files for i/o; get a CM.
   * We configure the CM with ConfigCM() later.
   ************************************************/
  watch = StopwatchCreate();
  
  if(!do_rsearch)
    {
      if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
	Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
      if (! CMFileRead(cmfp, &cm))
	Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
      if (cm == NULL) 
	Die("%s empty?\n", cmfile);
      CMFileClose(cmfp);
    }
  else  /* do RSEARCH */
    {
      /* Read in the alignment (with 1 seq) to build RSEARCH CM from */
      if ((afp = MSAFileOpen(alifile, format, NULL)) == NULL)
	Die("Alignment file %s could not be opened for reading", alifile);
      /* read file */
      if ((msa = MSAFileRead(afp)) == NULL) 
	Die ("Could not read in query\n");
      MSAFileClose(afp);

      /* Set up matrix */
      matrixname[0]          = '\0';
      if ((matfp = MatFileOpen (DEFAULT_RMATRIX, getenv("RNAMAT"), matrixname)) == NULL) 
	Die ("Failed to open matrix file\n%s\n", usage);
      if (! (fullmat = ReadMatrix(matfp)))
	Die ("Failed to read matrix file \n%s\n", usage);
      printf ("Matrix: %s\n", fullmat->name);
      /* print alphas and betas */
      printf ("Alpha: %.2f\n", ralpha);
      printf ("Beta: %.2f\n", rbeta);
      printf ("Alpha': %.2f\n", ralphap);
      printf ("Beta': %.2f\n", rbetap);
      /* check if alpha+beta sum too low  */
      rmin_alpha_beta_sum = get_min_alpha_beta_sum(fullmat);
      if (ralpha + rbeta < rmin_alpha_beta_sum) 
	Die ("alpha + beta must sum to less than %.2f\n", rmin_alpha_beta_sum);
      if (ralphap + rbetap < rmin_alpha_beta_sum) 
	Die ("alpha' + beta' must sum to less than %.2f\n", rmin_alpha_beta_sum);
  
      /* Build the CM */
      cm = build_cm (msa, fullmat, &rquerylen, ralpha, rbeta, ralphap, rbetap,
		     rbeginsc, rendsc);
      /* done with the matrix, free it */
      FreeMat(fullmat);

      cm->W = 2 * MSAAverageSequenceLength(msa);
      /* NO NULL MODEL SET! */
      cm->el_selfsc  = 0.; /* temporary to match RSEARCH */
      cm->iel_selfsc = Prob2Score(sreEXP2(cm->el_selfsc), 1.0);
    } 
  
  /* Set CM and CP9 parameters that can be changed at command line */
  cm->beta         = beta;     /* this will be DEFAULT_BETA unless set at command line */
  cm->hbandp       = hbandp;   /* this will be DEFAULT_HBANDP unless set at command line */
  cm->sc_boost     = sc_boost; /* this will be 0.0 unless set at command line */
  cm->cp9_sc_boost = sc_boost; /* this will be 0.0 unless set at command line */
  if (set_window) cm->W = set_W;

  /* Set the cutoffs */
  cm->cutoff_type = cm_cutoff_type;  /* this will be DEFAULT_CM_CUTOFF_TYPE unless set at command line */
  if(cm->cutoff_type == SCORE_CUTOFF)
    cm->cutoff = cm_sc_cutoff;
  else
    cm->cutoff = cm_e_cutoff;
  cm->cp9_cutoff_type = cp9_cutoff_type; /* this will be DEFAULT_CP9_CUTOFF_TYPE unless set at command line */
  if(cm->cp9_cutoff_type == SCORE_CUTOFF)
    cm->cp9_cutoff = cp9_sc_cutoff;
  else
    cm->cp9_cutoff = cp9_e_cutoff;
  
  /* Update cm->config_opts and cm->search_opts based on command line options */
  if(do_local)        cm->config_opts |= CM_CONFIG_LOCAL;
  if(do_zero_inserts) cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  if(!(do_qdb))       cm->search_opts |= CM_SEARCH_NOQDB;
  if(do_hmmonly)      cm->search_opts |= CM_SEARCH_HMMONLY;
  if(do_hmmfb)        cm->search_opts |= CM_SEARCH_HMMFB;
  if(do_hmmweinberg)  cm->search_opts |= CM_SEARCH_HMMWEINBERG;
  if(do_scan2hbands)  cm->search_opts |= CM_SEARCH_SCANBANDS;
  if(use_sums)        cm->search_opts |= CM_SEARCH_SUMS;
  if(do_inside)       cm->search_opts |= CM_SEARCH_INSIDE;
  if(!do_revcomp)     cm->search_opts |= CM_SEARCH_TOPONLY;
  if(!do_align)       cm->search_opts |= CM_SEARCH_NOALIGN;
  if(do_null2)        cm->search_opts |= CM_SEARCH_NULL2;
  if(do_cm_stats)     cm->search_opts |= CM_SEARCH_CMSTATS;
  if(do_cp9_stats)    cm->search_opts |= CM_SEARCH_CP9STATS;
  if(do_rsearch) 
    {
      /* TEMPORARILY: we don't do anything (no QDB, no CP9s, no CMLogoddisfy calls etc.) 
       * if in RSEARCH mode */
      /* RSEARCH: always local, never QDB, never null2 (?), never anything CP9 related */
      cm->search_opts |= CM_SEARCH_NOQDB;
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->search_opts |= CM_SEARCH_RSEARCH;
      cm->search_opts &= ~CM_SEARCH_NULL2;
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      cm->search_opts &= ~CM_SEARCH_HMMFB;
      cm->search_opts &= ~CM_SEARCH_HMMWEINBERG;
      cm->search_opts &= ~CM_SEARCH_SCANBANDS;
      cm->search_opts &= ~CM_SEARCH_CP9STATS;
      cm->flags       |= CM_IS_RSEARCH;
    }
  if(do_enforce)
    {
      cm->config_opts |= CM_CONFIG_ENFORCE;
      cm->enf_start    = enf_start; 
      cm->enf_seq      = enf_seq;
    }

  if(do_qdb) cm->config_opts |= CM_CONFIG_QDB;
  if(read_qdb)
    {
      /* read the bands from a file */
      if ((qdb_fp = fopen(qdb_file, "r")) == NULL)
	Die("failed to open QDB file %s", qdb_file);
      if(!(QDBFileRead(qdb_fp, cm, &preset_dmin, &preset_dmax)))
	{
	  Die("ERROR reading QDB file: %s.\nDoes it correspond (same number of states) to this model?\n", qdb_file);
	}
      fclose(qdb_fp);
    }
  else
    preset_dmin = preset_dmax = NULL;
  
  /*******************************
   * Open the sequence (db) file *
   *******************************/
  status = esl_sqfile_Open(seqfile, format, NULL, &dbfp);
  if (status == eslENOTFOUND) esl_fatal("No such file."); 
  else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
  else if (status == eslEINVAL) esl_fatal("Can’t autodetect stdin or .gz."); 
  else if (status != eslOK) esl_fatal("Failed to open sequence database file, code %d.", status); 

#ifdef USE_MPI
}   /* End of first block that is only done by master process */
  /* Barrier for debugging */
  MPI_Barrier(MPI_COMM_WORLD);

  /* Here we need to broadcast the following parameters:
     num_samples, W, W_scale, and the CM */
  broadcast_cm(&cm, mpi_my_rank, mpi_master_rank);
  search_first_broadcast(&num_samples, &W_scale, mpi_my_rank, 
			 mpi_master_rank);
  
#endif

  /* Configure the CM for search based on cm->config_opts 
   * and cm->search_opts. Set local mode, make cp9 HMM, 
   * calculate QD bands etc., preset_dmin and preset_dmax
   * are NULL unless --qdbfile, and you can't enable 
   * --qdbfile in MPI mode (we check for this and die if
   * we're trying to above). */
  ConfigCM(cm, preset_dmin, preset_dmax);

  cons = CreateCMConsensus(cm, 3.0, 1.0); 
#ifdef USE_MPI
  if(mpi_my_rank == mpi_master_rank)
    {
#endif
  if(do_bdump && (!(cm->search_opts & CM_SEARCH_NOQDB))) 
    {
      printf("beta:%f\n", cm->beta);
      debug_print_bands(cm, cm->dmin, cm->dmax);
      PrintDPCellsSaved(cm, cm->dmin, cm->dmax, cm->W);
    }
#ifdef USE_MPI
    }
#endif

  /**************************************************
   * Make the histogram(s)
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
  if (sample_length == 0) sample_length = W_scale * cm->W;
  /*printf("0 W: %d W_scale: %f sample_length: %d\n", W, W_scale, sample_length);*/
  printf("cm stats: %d (%d) cp9 stats: %d (%d)\n", (cm->search_opts & CM_SEARCH_CMSTATS),
    do_cm_stats, (cm->search_opts & CM_SEARCH_CP9STATS), do_cp9_stats);
  
  if (cm->search_opts & CM_SEARCH_CMSTATS) 
    {
      printf("CALCING CM STATS\n");
#ifdef USE_MPI
      /*StopwatchZero(watch);
	StopwatchStart(watch);*/
      if (mpi_num_procs > 1)
	parallel_make_histogram(gc_ct, partitions, num_partitions,
				cm, num_samples, sample_length, 
				FALSE, /* we're not doing CP9 stats */
				mpi_my_rank, mpi_num_procs, mpi_master_rank);
      else 
#endif
	serial_make_histogram (gc_ct, partitions, num_partitions,
			       cm, num_samples, sample_length, 
			       FALSE, /* we're not doing CP9 stats */
			       TRUE); /* use easel, discard eventually */
    }      
  /*StopwatchStop(watch);
    StopwatchDisplay(stdout, "\nCPU time (histogram): ", watch);*/

  /* If we're calculating stats for the CP9, build CP9 histograms */
  if (cm->search_opts & CM_SEARCH_CP9STATS) 
    {
      printf("CALCING CP9 STATS\n");
#ifdef USE_MPI
      /*StopwatchZero(watch);
	StopwatchStart(watch);*/
      if (mpi_num_procs > 1)
	parallel_make_histogram(gc_ct, partitions, num_partitions,
				cm, num_samples, sample_length, 
				TRUE, /* we are doing CP9 stats */
				mpi_my_rank, mpi_num_procs, mpi_master_rank);
      else 
#endif
	serial_make_histogram (gc_ct, partitions, num_partitions,
			       cm, num_samples, sample_length, 
			       TRUE, /* we are doing CP9 stats */
			       TRUE);/* use easel, discard eventually */
    }      
  printf("done calcing CP9 stats\n");
  /*StopwatchStop(watch);
    StopwatchDisplay(stdout, "\nCPU time (histogram): ", watch);*/
      
#ifdef USE_MPI
  if (mpi_my_rank == mpi_master_rank) {
#endif
    
    /* Set CM mu from K, lambda, N */
    if (cm->search_opts & CM_SEARCH_CMSTATS)
      {
	for (i=0; i<GC_SEGMENTS; i++) 
	  cm->mu[i] = log(cm->K[i]*N)/cm->lambda[i];
	//debug_print_stats(partitions, num_partitions, cm->lambda, cm->mu);
      }    
    /* else they've been set to default 0.0s in ConfigCM() */

    /* Set CP9 mu from K, lambda, N */
    if (cm->search_opts & CM_SEARCH_CP9STATS)
      {
	for (i=0; i<GC_SEGMENTS; i++) 
	  cm->cp9_mu[i] = log(cm->cp9_K[i]*N)/cm->cp9_lambda[i];
	//debug_print_stats(partitions, num_partitions, cm->cp9_lambda, cm->cp9_mu);
      }    
    /* else they've been set to default 0.0s in ConfigCM() */

    if (cm->search_opts & CM_SEARCH_CMSTATS) 
      {
	printf ("CM statistics calculated with simulation of %d samples of length %d\n", num_samples, sample_length);
	if (num_partitions == 1) 
	  printf ("No partition points\n");
	else 
	  {
	    printf ("Partition points are: ");
	    for (i=0; i<=GC_SEGMENTS; i++)
	      if (partitions[i] != partitions[i-1]) 
		printf ("%d ", i);
	  }
	for (i=0; i<GC_SEGMENTS; i++) 
	  ;/*printf ("GC = %d\tlambda = %.4f\tmu = %.4f\n", i, lambda[i], mu[i]);*/
	printf ("N = %ld\n", N);
	if (cm->cutoff_type == SCORE_CUTOFF) 
	  printf ("Using score cutoff of %.2f\n", cm->cutoff);
	else 
	  printf ("Using E cutoff of %.2f\n", cm->cutoff);
      } 
    else 
	{
	  printf ("lambda and K undefined -- no statistics\n");
	  printf ("Using score cutoff of %.2f\n", cm->cutoff);
	}
    fflush(stdout);

    if (cm->search_opts & CM_SEARCH_CP9STATS) 
      {
	printf ("CP9 statistics calculated with simulation of %d samples of length %d\n", num_samples, sample_length);
	if (num_partitions == 1) 
	  printf ("No partition points\n");
	else 
	  {
	    printf ("Partition points are: ");
	    for (i=0; i<=GC_SEGMENTS; i++)
	      if (partitions[i] != partitions[i-1]) 
		printf ("%d ", i);
	  }
	for (i=0; i<GC_SEGMENTS; i++) 
	  ;/*printf ("GC = %d\tlambda = %.4f\tmu = %.4f\n", i, lambda[i], mu[i]);*/
	printf ("N = %ld\n", N);
	if (cm->cp9_cutoff_type == SCORE_CUTOFF) 
	  printf ("Using score cutoff of %.2f\n", cm->cp9_cutoff);
	else 
	  printf ("Using E cutoff of %.2f\n", cm->cp9_cutoff);
      } 
    else if(cm->search_opts & CM_SEARCH_HMMONLY || cm->search_opts & CM_SEARCH_HMMWEINBERG)
	{
	  printf ("lambda and K undefined -- no statistics\n");
	  printf ("Using score cutoff of %.2f\n", cm->cp9_cutoff);
	}
    fflush(stdout);

    /*************************************************
     * End of making the histogram(s).
     *************************************************/
    /*************************************************
     *    Do the search
     *************************************************/
#ifdef USE_MPI
    } /* Done with second master-only block */
  
  /* Second broadcast, send N, and the EVD stats if
   * we're doing CM stats and/or CP9 stats. */
  search_second_broadcast(&cm, &N, mpi_my_rank, mpi_master_rank);

  printf("cm->cutoff: %f cp9->cutoff: %f rank: %d\n", cm->cutoff, cm->cp9_cutoff, mpi_my_rank);

  if (mpi_num_procs > 1)
    parallel_search_database (dbfp, cm, cons,
			      mpi_my_rank, mpi_master_rank, mpi_num_procs);
  else
#endif
    serial_search_database (dbfp, cm, cons);
  
  FreeCM(cm);
  free(partitions);
  StopwatchFree(watch);
  esl_sqfile_Close(dbfp);
  FreeCMConsensus(cons);
  free(gc_ct);
  if(do_rsearch)
    MSAFree(msa);
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
  for (i=0; i<GC_SEGMENTS; i++) 
    {
      /*printf("i: %d\n", i);*/
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

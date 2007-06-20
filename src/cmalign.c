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
#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "squid.h"		/* general sequence analysis library    */
#include "easel.h"		/* newer general seq analysis library   */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "sre_stack.h"
#include "hmmband.h"         
#include "cplan9.h"
#include "cm_postprob.h"
#include "mpifuncs.h"
#include "cm_dispatch.h"
#include "esl_sqio.h"		

static int in_mpi;

static void include_alignment(char *seqfile, int use_rf, float gapthresh, CM_t *cm, 
			      ESL_SQ ***ret_sq, Parsetree_t ***ret_tr,
			      int *ret_nseq, int *ret_nnewseq);
static int compare_cms(CM_t *cm1, CM_t *cm2);

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
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

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
static char banner[] = "mpicmalign - align sequences to an RNA CM";
static char usage[]  = "\
Usage: mpicmalign [-options] <cmfile> <sequence file>\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -l     : local; align locally w.r.t. the model\n\
   -o <f> : output the alignment file to file <f>\n\
   -q     : quiet; suppress verbose banner\n\
";
#else
static char banner[] = "cmalign - align sequences to an RNA CM";
static char usage[]  = "\
Usage: cmalign [-options] <cmfile> <sequence file>\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -l     : local; align locally w.r.t. the model\n\
   -o <f> : output the alignment file to file <f>\n\
   -q     : quiet; suppress verbose banner\n\
";
#endif

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
   --zeroinserts : set insert emission scores to 0\n\
   --sub         : build sub CM for columns b/t HMM predicted start/end points\n\
   --elsilent    : disallow local end (EL) emissions\n\
   --enfstart <n>: enforce MATL stretch starting at consensus position <n>\n\
   --enfseq   <s>: enforce MATL stretch starting at --enfstart <n> emits seq <s>\n\
\n\
  * HMM banded alignment related options (*in development*):\n\
   --hbanded     : accelerate using CM plan 9 HMM banded CYK aln algorithm\n\
   --tau <x>     : set tail loss prob for --hbanded to <x> [default: 1E-7]\n\
   --hsafe       : realign (non-banded) seqs with HMM banded CYK score < 0 bits\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
   --hmmonly     : align with the CM Plan 9 HMM (no alignment given, just score)\n\
\n\
  * Query dependent banded (qdb) alignment related options:\n\
   --qdb         : use query dependent banded CYK alignment algorithm\n\
   --beta <x>    : set tail loss prob for QDB to <x> [default:1E-7]\n\
\n\
  * Options for including the alignment used to build the CM in the output:\n\
   --withali <f> : incl. alignment in <f> (must be aln <cm file> was built from)\n\
   --rf          : (only with --withali) cmbuild --rf was used\n\
   --gapthresh <x>:(only with --withali) cmbuild --gapthresh <x> was used\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-l", TRUE, sqdARG_NONE },
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--nosmall",    FALSE, sqdARG_NONE },
  { "--regress",    FALSE, sqdARG_STRING },
  { "--qdb",        FALSE, sqdARG_NONE },
  { "--beta",       FALSE, sqdARG_FLOAT},
  { "--tfile",      FALSE, sqdARG_STRING },
  { "--banddump"  , FALSE, sqdARG_INT},
  { "--full",       FALSE, sqdARG_NONE },
  { "--dlev",       FALSE, sqdARG_INT },
  { "--hbanded",    FALSE, sqdARG_NONE },
  { "--tau",       FALSE, sqdARG_FLOAT},
  { "--hsafe",      FALSE, sqdARG_NONE},
  { "--hmmonly",    FALSE, sqdARG_NONE },
  /*  { "--sums",       FALSE, sqdARG_NONE},*/ 
  { "--time",       FALSE, sqdARG_NONE},
  { "--inside",     FALSE, sqdARG_NONE},
  { "--outside",    FALSE, sqdARG_NONE},
  { "--post",       FALSE, sqdARG_NONE},
  { "--checkpost",  FALSE, sqdARG_NONE},
  { "--sub",        FALSE, sqdARG_NONE},
  { "--elsilent",   FALSE, sqdARG_NONE},
  { "--enfstart",   FALSE, sqdARG_INT},
  { "--enfseq",     FALSE, sqdARG_STRING},
  { "--zeroinserts",FALSE, sqdARG_NONE},
  { "--withali",    FALSE, sqdARG_STRING },
  { "--gapthresh",  FALSE, sqdARG_FLOAT},
  { "--rf",         FALSE, sqdARG_NONE }
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from                     */	
  CMFILE          *cmfp;	/* open CM file                             */
  char            *seqfile;     /* file to read sequences from              */
  ESL_SQFILE	  *seqfp;       /* open seqfile for reading                 */
  ESL_SQ	 **sq;          /* the target sequences to align            */
  char            *optname;     /* name of option found by Getopt()         */
  char            *optarg;      /* argument found by Getopt()               */
  int              optind;      /* index in argv[]                          */

  int              status;
  int              format;      /* format of sequence file                  */
  int              nseq;	/* number of seqs                           */
  CM_t            *cm;          /* a covariance model                       */
  Parsetree_t    **tr;          /* parse trees for the sequences            */
  CP9trace_t     **cp9_tr;      /* CP9 traces for the sequences             */
  MSA             *msa;         /* alignment that's created                 */
  char            *outfile;	/* optional output file name                */
  FILE            *ofp;         /* an open output file                      */
  int              i;		/* counter over sequences/parsetrees        */
  char            *regressfile; /* regression test data file                */
  char            *tracefile;	/* file to dump debugging traces to         */
  CMConsensus_t   *con;         /* consensus information for the CM         */

  /* Options, modifiable via the command line */
  int              be_quiet;	/* TRUE to suppress verbose output & banner */
  int              do_local;	/* TRUE to config the model in local mode   */
  int              do_small;	/* TRUE to do divide and conquer alignments */
  int              do_full;     /* TRUE to include all match cols in aln    */
  int              do_qdb;      /* TRUE to do qdb CYK (either d&c or full)  */
  int              bdump_level; /* verbosity level for --banddump, 0 is OFF */
  int              debug_level; /* verbosity level for debugging printf's s */
  int              do_timings;  /* TRUE to print timings, FALSE not to      */
  double           qdb_beta;	/* tail loss prob for query dependent bands */

  /* HMM banded alignment */
  int              do_hbanded;  /* TRUE to use CP9 HMM to band CYK          */
  double           tau;         /* tail loss probability for hmm bands      */
  int              use_sums;    /* TRUE: use the posterior sums w/HMM bands */
  int              do_hmmonly;  /* TRUE: align with the HMM, not the CM     */
  int              do_hsafe;    /* TRUE: realign seqs with banded sc < 0    */
  /* Alternatives to CYK */
  int              do_inside;   /* TRUE to use Inside algorithm, not CYK    */
  int              do_outside;  /* TRUE to use Outside algorithm, not CYK   */
  int              do_check;    /* TRUE to check Inside & Outside probs     */
  int              do_post;     /* TRUE to do Inside/Outside, not CYK       */
  char           **postcode;    /* posterior decode array of strings        */
  char            *apostcode;   /* aligned posterior decode array           */

  /* the --sub option */
  int              do_sub;      /* TRUE to use HMM to infer start and end   *
				 * point of each seq and build a separate   *
				 * sub CM for alignment of that seq         */
  int              do_fullsub;  /* TRUE to only remove structure outside HMM*
				 * predicted start (spos) and end points    *
				 * (epos) (STILL IN DEVELOPMENT!)           */
  int              do_elsilent; /* TRUE to disallow EL emissions            */
  int              do_zero_inserts; /* TRUE to zero insert emission scores */

  /* The enforce option (--enfstart and --enfseq), added specifically for   *
   * enforcing the template region for telomerase RNA searches              */
  int   do_enforce;             /* TRUE to enforce a MATL stretch is used   */
  int   enf_cc_start;           /* first consensus position to enforce      */
  char *enf_seq;                /* the subsequence to enforce emitted by    *
                                 * MATL nodes starting at cm->enf_start     */

  /* variables for appending query aln CM was built from to output aln      */
  int    do_withali;            /* TRUE to incl. original query aln         */
  char  *withali;               /* name of additional aln file to align     */
  int    withali_nseq;          /* num seqs in withali ali file, 0 if none  */
  int    ip;                    /* used to  assign post codes w/--withali   */
  int    use_rf;		/* TRUE was used #=RF to define consensus   */
  float  gapthresh;		/* 0=all cols inserts; 1=all cols consensus */

  in_mpi = 0;			/* to silence overzealous compiler warnings */
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  int mpi_my_rank;              /* My rank in MPI                           */
  int mpi_num_procs;            /* Total number of processes                */
  int mpi_master_rank;          /* Rank of master process                   */
  Stopwatch_t  *mpi_watch;      /* for timings in MPI mode                  */

  /* Initialize MPI, get values for rank and num procs */
  MPI_Init (&argc, &argv);
  
  atexit (exit_from_mpi);
  in_mpi = 1;                /* Flag for exit_from_mpi() */

  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_num_procs);

  /* Determine master process.  This is the lowest ranking one that can do I/O */
  mpi_master_rank = get_master_rank (MPI_COMM_WORLD, mpi_my_rank);

  /*printf("B CS rank: %4d master: %4d num: %4d\n", mpi_my_rank, mpi_master_rank, mpi_num_procs);*/

  /* If I'm the master, do the following set up code -- parse arguments, read
     in matrix and query, build model */
  if (mpi_my_rank == mpi_master_rank) {
    mpi_watch = StopwatchCreate(); 
    StopwatchZero(mpi_watch);
    StopwatchStart(mpi_watch);
#endif

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
  do_qdb      = FALSE;
  qdb_beta    = DEFAULT_BETA;
  do_full     = FALSE;
  bdump_level = 0;
  debug_level = 0;
  do_hbanded  = FALSE;
  do_hmmonly  = FALSE;
  do_hsafe    = FALSE;
  tau         = DEFAULT_TAU;
  use_sums    = FALSE;
  do_timings  = FALSE;
  do_inside   = FALSE;
  do_outside  = FALSE;
  do_post     = FALSE;
  do_check    = FALSE;
  do_sub      = FALSE;
  do_fullsub  = FALSE;
  do_elsilent = FALSE;
  do_zero_inserts=FALSE;
  do_enforce  = FALSE;
  enf_cc_start= 0;
  enf_seq     = NULL;
  do_withali  = FALSE;
  withali     = NULL;
  use_rf      = FALSE;
  gapthresh   = 0.5;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-l")          == 0) do_local    = TRUE;
    else if (strcmp(optname, "-o")          == 0) outfile     = optarg;
    else if (strcmp(optname, "-q")          == 0) be_quiet    = TRUE;
    else if (strcmp(optname, "--nosmall")   == 0) do_small    = FALSE;
    else if (strcmp(optname, "--regress")   == 0) regressfile = optarg;
    else if (strcmp(optname, "--qdb")       == 0) do_qdb       = TRUE;
    else if (strcmp(optname, "--beta")      == 0) qdb_beta     = atof(optarg);
    else if (strcmp(optname, "--full")      == 0) do_full      = TRUE;
    else if (strcmp(optname, "--banddump")  == 0) bdump_level  = atoi(optarg);
    else if (strcmp(optname, "--tfile")     == 0) tracefile    = optarg;
    else if (strcmp(optname, "--dlev")      == 0) debug_level  = atoi(optarg);
    else if (strcmp(optname, "--time")      == 0) do_timings   = TRUE;
    else if (strcmp(optname, "--inside")    == 0) do_inside    = TRUE;
    else if (strcmp(optname, "--outside")   == 0) do_outside   = TRUE; 
    else if (strcmp(optname, "--post")      == 0) do_post      = TRUE;
    else if (strcmp(optname, "--checkpost") == 0) do_check     = TRUE;
    else if (strcmp(optname, "--sub")       == 0) do_sub       = TRUE; 
    else if (strcmp(optname, "--elsilent")  == 0) do_elsilent  = TRUE;
    else if (strcmp(optname, "--hbanded")   == 0) { do_hbanded = TRUE; do_small = FALSE; }
    else if (strcmp(optname, "--tau")       == 0) tau          = atof(optarg);
    else if (strcmp(optname, "--hsafe")     == 0) do_hsafe     = TRUE;
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--hmmonly")   == 0) { do_hmmonly   = TRUE; do_timings = TRUE; } 
    else if (strcmp(optname, "--enfstart")  == 0) { do_enforce = TRUE; enf_cc_start = atoi(optarg); }
    else if (strcmp(optname, "--enfseq")    == 0) { do_enforce = TRUE; enf_seq = optarg; } 
    else if (strcmp(optname, "--zeroinserts")== 0) do_zero_inserts = TRUE;
    else if (strcmp(optname, "--withali")   == 0) { do_withali = TRUE; withali = optarg; } 
    else if (strcmp(optname, "--rf")        == 0) { do_withali = TRUE; use_rf  = TRUE;   } 
    else if (strcmp(optname, "--gapthresh") == 0) { do_withali = TRUE; gapthresh = atof(optarg);} 
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

  /* Check for incompatible or misused options */
  if(do_inside && do_outside)
    Die("Please pick either --inside or --outside (--outside will run Inside()\nalso and check to make sure Inside() and Outside() scores agree).\n");
  if((do_inside || do_outside) && (outfile != NULL))
    Die("-o <f> cannot be used in combination with --inside or --outside, as no alignment is created.\n");
  if(do_sub && do_local)
    Die("--sub and -l combination not supported.\n");
  if(do_sub && do_qdb)
    Die("Please pick either --sub or --qdb.\n");
  /*if(do_hmmonly)
    Die("--hmmonly not yet implemented.\n");*/
  if(do_enforce && enf_seq == NULL)
    Die("--enfstart only makes sense with --enfseq also.\n");
  if(do_enforce && enf_cc_start == 0)
    Die("--enfseq only makes sense with --enfstart (which can't be 0) also.\n");
  if(bdump_level > 0 && (!do_qdb))
    Die("--bdump option does not work with --noqdb.\n");
  if(do_qdb && do_hbanded)
    Die("--qdb and --hbanded combo not supported, pick one.\n");
  if (bdump_level > 3) 
    Die("Highest available --banddump verbosity level is 3\n%s", usage);
  if (do_hsafe && !do_hbanded)
    Die("--hsafe only makes sense with --hbanded\n%s", usage);
  if(withali != NULL && (do_inside || do_outside || do_hmmonly))
    Die("--withali does not work with --hmmonly, --inside or --outside\n%s", usage);
  if(do_withali && withali == NULL)
    Die("--rf and --gapthresh only make sense with --withali\n%s", usage);
  if (do_local && do_hbanded)
    {
      printf("Warning: banding with an HMM (--hbanded) and allowing\nlocal alignment (-l). This may not work very well.\n");
    } 

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
   * Input the CM  (we configure it with ConfigCM() later).
   *****************************************************************/

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);

  cm->beta   = qdb_beta; /* this will be DEFAULT_BETA unless changed at command line */
  cm->tau    = tau;      /* this will be DEFAULT_TAU unless changed at command line */

  /* Update cm->config_opts and cm->align_opts based on command line options */
  if(do_local)        
    {
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
    }
  if(do_elsilent)     cm->config_opts |= CM_CONFIG_ELSILENT;
  if(do_zero_inserts) cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  if(do_hbanded)      cm->align_opts  |= CM_ALIGN_HBANDED;
  if(use_sums)        cm->align_opts  |= CM_ALIGN_SUMS;
  if(do_sub)          cm->align_opts  |= CM_ALIGN_SUB;
  if(do_fullsub)      cm->align_opts  |= CM_ALIGN_FSUB;
  if(do_hmmonly)      cm->align_opts  |= CM_ALIGN_HMMONLY;
  if(do_inside)       cm->align_opts  |= CM_ALIGN_INSIDE;
  if(do_outside)      cm->align_opts  |= CM_ALIGN_OUTSIDE;
  if(!do_small)       cm->align_opts  |= CM_ALIGN_NOSMALL;
  if(do_post)         cm->align_opts  |= CM_ALIGN_POST;
  if(do_timings)      cm->align_opts  |= CM_ALIGN_TIME;
  if(do_check)        cm->align_opts  |= CM_ALIGN_CHECKINOUT;
  if(do_hsafe)        cm->align_opts  |= CM_ALIGN_HMMSAFE;
  if(do_enforce)
    {
      cm->config_opts |= CM_CONFIG_ENFORCE;
      cm->enf_start    = EnforceFindEnfStart(cm, enf_cc_start); 
      cm->enf_seq      = enf_seq;
    }
  if(do_qdb)          
    { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }
  /*****************************************************************
   * Open the target sequence file
   *****************************************************************/
  status = esl_sqfile_Open(seqfile, format, NULL, &seqfp);
  if (status == eslENOTFOUND) esl_fatal("No such file."); 
  else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
  else if (status == eslEINVAL) esl_fatal("Can’t autodetect stdin or .gz."); 
  else if (status != eslOK) esl_fatal("Failed to open sequence database file, code %d.", status); 

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

  /****************************************************************************
   * Do the work: read in seqs and align them either serially or in parallel. 
   *
   * This is done within {serial, parallel}_align_targets(). The function
   * that actually aligns the targets is actually_align_targets(), The code 
   * was organized in this way to make implementing MPI parallelization easier
   * (EPN 01.04.07).
   ****************************************************************************/

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  }   /* End of first block that is only done by master process */
  /* Barrier for debugging */
  MPI_Barrier(MPI_COMM_WORLD);

  /* Broadcast the CM, before we configure it. */
  broadcast_cm(&cm, mpi_my_rank, mpi_master_rank);

  if (mpi_num_procs > 1)
    {
      /* Configure the CM for alignment based on cm->config_opts and cm->align_opts.
       * set local mode, make cp9 HMM, calculate QD bands etc. */
      ConfigCM(cm, NULL, NULL);
      if(cm->config_opts & CM_CONFIG_ENFORCE) ConfigCMEnforce(cm);

      parallel_align_targets(seqfp, cm, &sq, &tr, &postcode, NULL, &nseq,
			     bdump_level, debug_level, 
			     TRUE, /* be_quiet=TRUE we don't print scores in MPI
				    * mode yet, b/c they get jumbled due to potentially
				    * multiple simultaneous writes to stdout */
			     mpi_my_rank, mpi_master_rank, mpi_num_procs);
    }
  else
#endif /* (end of if USE_MPI) */
    {
      /* Configure the CM for alignment based on cm->config_opts and cm->align_opts.
       * set local mode, make cp9 HMM, calculate QD bands etc. */
      ConfigCM(cm, NULL, NULL);
      if(cm->config_opts & CM_CONFIG_ENFORCE) ConfigCMEnforce(cm);
      if(do_hmmonly)
	serial_align_targets(seqfp, cm, &sq, &tr, &postcode, &cp9_tr, &nseq, bdump_level, debug_level, 
			     be_quiet);
      else
	serial_align_targets(seqfp, cm, &sq, &tr, &postcode, NULL, &nseq, bdump_level, debug_level, 
			     be_quiet);
    }
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  if (mpi_my_rank == mpi_master_rank) {
#endif  
    /***************************************************************                     
     * Create the MSA.                                                                   
     ****************************************************************/                   
    msa = NULL;                                                                          
    if(!((cm->align_opts & CM_ALIGN_INSIDE) || (cm->align_opts & CM_ALIGN_OUTSIDE)))
      {                                                                                  
	/* optionally include a fixed alignment provided with --withali */
	withali_nseq = 0;
	if(withali != NULL)
	  include_alignment(withali, use_rf, gapthresh, cm, &sq, &tr, &nseq, &withali_nseq);

	if(!do_hmmonly)
	  msa = ESL_Parsetrees2Alignment(cm, sq, NULL, tr, nseq, do_full);                 
	else
	  {
	    con = CreateCMConsensus(cm, 3.0, 1.0);
	    msa = CP9Traces2Alignment(cm, sq, NULL, nseq, cp9_tr, do_full);
	    FreeCMConsensus(con);
	  }

	if(cm->align_opts & CM_ALIGN_POST)                                              
        {                                                                              
	  if(postcode == NULL)
	    Die("ERROR CM_ALIGN_POST flag is up, but {serial,parallel}_align_targets() did not return post codes.\n");
          for (i = withali_nseq; i < nseq; i++)                                                   
            {                                                                          
	      ip = i - withali_nseq;
	      MakeAlignedString(msa->aseq[i], msa->alen, postcode[ip], &apostcode);     
              MSAAppendGR(msa, "POST", i, apostcode);                                  
              free(apostcode);                                                         
              free(postcode[ip]);                                                       
            }                                                                          
          free(postcode);                                                              
        }                                                                              

	/*****************************************************************
	 * Output the alignment.
	 *****************************************************************/
	
	printf("\n");
	if (outfile != NULL && (ofp = fopen(outfile, "w")) != NULL) 
	  {
	    WriteStockholm(ofp, msa);
	    printf("Alignment saved in file %s\n", outfile);
	    fclose(ofp);
	  }
	else
	  WriteStockholm(stdout, msa);
	
	/* Detailed traces for debugging training set. */
	if (tracefile != NULL)       
	  {
	    if(do_hmmonly)
	      { printf("%-40s ... ", "Saving CP9 HMM traces"); fflush(stdout); }
	    else
	      { printf("%-40s ... ", "Saving CM parsetrees"); fflush(stdout); }
	    if ((ofp = fopen(tracefile,"w")) == NULL)
	      Die("failed to open trace file %s", tracefile);
	    for (i = 0; i < msa->nseq; i++) 
	      {
		fprintf(ofp, "> %s\n", msa->sqname[i]);
		if(do_hmmonly)
		  {
		    fprintf(ofp, "  SCORE : %.2f bits\n", CP9TraceScore(cm->cp9, sq[i]->dsq, cp9_tr[i]));
		    CP9PrintTrace(ofp, cp9_tr[i], cm->cp9, sq[i]->dsq);
		  }
		else
		  {
		    fprintf(ofp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], sq[i]->dsq, FALSE));;
		    ParsetreeDump(ofp, tr[i], cm, sq[i]->dsq);
		  }
		fprintf(ofp, "//\n");
	      }
	    fclose(ofp);
	    printf("done. [%s]\n", tracefile);
	  }
	
	if (regressfile != NULL && (ofp = fopen(regressfile, "w")) != NULL) 
	  {
	    /* Must delete author info from msa, because it contains version
	     * and won't diff clean in regression tests. */
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
	if((!(do_inside || do_outside)) && !do_hmmonly) FreeParsetree(tr[i]);
	if(do_hmmonly) CP9FreeTrace(cp9_tr[i]);
	esl_sq_Destroy(sq[i]);
      }
    esl_sqfile_Close(seqfp);

    if(!(do_inside || do_outside)) MSAFree(msa);
    free(sq);
    free(tr);
    if(do_hmmonly) free(cp9_tr);
    SqdClean();

#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  } /* end  of block to convert traces to alignment and clean up */
  if(mpi_my_rank == mpi_master_rank)
    {
      StopwatchStop(mpi_watch);
      StopwatchDisplay(stdout, "MPI Master node time: ", mpi_watch);
      StopwatchFree(mpi_watch);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  in_mpi = 0;
#endif
  FreeCM(cm);
#if defined(USE_MPI) && defined(MPI_EXECUTABLE)
  /*printf("EXITING rank:%d\n", mpi_my_rank);*/
#endif
  return EXIT_SUCCESS;
}


/* Function: include_alignment()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Given the name of a multiple alignment file,
 *           align that alignment to the CM, and add parsetrees
 *           to an existing array of parsetrees. 
 *           Based on HMMER 2.4devl's hmmalign.c::include_alignment.
 *
 * Args:     seqfile      - the alignment file name to include
 *           use_rf       - TRUE to use #=GC RF line for consensus
 *           gapthresh    - gap threshold for consensus columns
 *           CM           - the covariance model
 *           seqfp        - the alignment file with ali to include
 *           ret_sq       - RETURN: the sequences (EASEL)
 *           ret_tr       - RETURN: the parsetrees for seqs in seqfp
 *           ret_nseq     - RETURN: the new total number of seqs 
 *           ret_nnewseq  - RETURN: number of new seqs added
 * 
 * Returns:  new, realloc'ed arrays for sq, tr; nseq is
 *           increased by number of seqs in seqfp.a
 */
static void include_alignment(char *seqfile, int use_rf, float gapthresh, CM_t *cm, 
			      ESL_SQ ***ret_sq, Parsetree_t ***ret_tr,
			      int *ret_nseq, int *ret_nnewseq)
{
  /* TODO: Easelfy! */
  int           format;	  /* format of alignment file */
  MSA          *msa;      /* alignment we're including  */
  MSAFILE      *afp;      /* open alignment file      */
  int           i;	  /* counter over aseqs       */
  int           ip;	  /* offset counter over aseqs */
  CM_t         *new_cm;   /* CM built from MSA, we check it has same guide tree as 'cm' */
  CM_t         *newer_cm; /* used briefly if we rebalance new_cm */
  char        **dsq;	  /* digitized aligned sequences             */
  Parsetree_t  *mtr;      /* master structure tree from the MSA */
  char        **uaseq;    /* unaligned seqs, dealigned from the MSA */
  int           apos;     /*   aligned position index */
  int           uapos;    /* unaligned position index */
  int           x;        /* counter of parsetree nodes */
  int         **map;      /* [0..msa->nseq-1][0..msa->alen] map from aligned
			   * positions to unaligned (non-gap) positions */

  format = MSAFILE_UNKNOWN;	/* invoke Babelfish */
  if ((afp = MSAFileOpen(seqfile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", seqfile);
  if ((msa = MSAFileRead(afp)) == NULL)
    Die("Failed to read an alignment from %s\n", seqfile);
  MSAFileClose(afp);

  dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);

  /* Some input data cleaning. */
  if (msa->ss_cons == NULL) 
    Die("failed... Alignment has no consensus structure annotation.");
  if (! clean_cs(msa->ss_cons, msa->alen))
    Die("failed... Failed to parse consensus structure annotation.");

  /* Build a CM from a master guide tree built from the msa, 
   * then check to make sure this CM has same emit map as the CM
   * we've had passed in. This is fragile and hopefully temporary. 
   * Another solution would be to use a checksum, but CM files don't 
   * have checksums yet.
   */
  HandModelmaker(msa, dsq, use_rf, gapthresh, &new_cm, &mtr);
  if(!(compare_cms(cm, new_cm)))
    {
      newer_cm = CMRebalance(new_cm);
      FreeCM(new_cm);
      new_cm = newer_cm;
      if(!(compare_cms(cm, new_cm)))
	Die("ERROR, in include_alignment(), CMs differ (even after rebalancing).\nAre you sure you used %s to build this CM?\nDid you use --rf or --gapthresh <x> options with cmbuild?\nIf so, use them again them with cmalign.", seqfile);
    }
  /* If we get here, a CM built from the seqfile MSA has same node 
   * architecture as the CM that was passed in, now we don't care
   * about the new_cm anymore, free it. 
   */
  FreeCM(new_cm);

  /* For each seq in the MSA, map the aligned sequences coords to 
   * the unaligned coords, we stay in digitized seq coords (1..alen),
   * we need this for converting parsetrees from Transmogrify, which
   * have emitl and emitr in aligned coords to unaligned coords, so 
   * we can call Parsetrees2Alignment() with them. */
  map = MallocOrDie(sizeof(int *) * msa->nseq);
  for (i = 0; i < msa->nseq; i++)
    {
      map[i] = MallocOrDie(sizeof(int) * (msa->alen+1));
      map[i][0] = -1; /* invalid */
      uapos = 1;
      for(apos = 0; apos < msa->alen; apos++)
	{
	  if (!isgap(msa->aseq[i][apos]))
	    map[i][(apos+1)] = uapos++;
	  else
	    map[i][(apos+1)] = -1;
	}
    }

  DealignAseqs(msa->aseq, msa->nseq, &uaseq);
  *ret_tr  = ReallocOrDie((*ret_tr), (sizeof(Parsetree_t *) * (*ret_nseq + msa->nseq)));
  *ret_sq  = ReallocOrDie((*ret_sq), (sizeof (ESL_SQ *)     * (*ret_nseq + msa->nseq)));

  /* Swap some pointers so the included alignment appears at the top of the output 
   * alignment instead of the bottom. */
  for(i = 0; i < *ret_nseq; i++)
    {
      ip = i + msa->nseq;
      (*ret_tr)[ip] = (*ret_tr)[i];
      (*ret_sq)[ip] = (*ret_sq)[i];
    }

  /* Transmogrify each aligned seq to get a parsetree */
  for (i = 0; i < msa->nseq; i++)
    {
      /*ip = i + *ret_nseq;*/
      (*ret_tr)[i] = Transmogrify(cm, mtr, dsq[i], msa->aseq[i], msa->alen);

      /* ret_tr[i] is in alignment coords, convert it to unaligned coords, */
      for(x = 0; x < (*ret_tr)[i]->n; x++)
	{
	  if((*ret_tr)[i]->emitl[x] != -1)
	    (*ret_tr)[i]->emitl[x] = map[i][(*ret_tr)[i]->emitl[x]];
	  if((*ret_tr)[i]->emitr[x] != -1)
	    (*ret_tr)[i]->emitr[x] = map[i][(*ret_tr)[i]->emitr[x]];
	}
      (*ret_sq)[i]      = esl_sq_CreateFrom(msa->sqname[i], uaseq[i], NULL, NULL, NULL);
      (*ret_sq)[i]->dsq = DigitizeSequence ((*ret_sq)[i]->seq, ((*ret_sq)[i]->n));
    }
  *ret_nseq    += msa->nseq;
  *ret_nnewseq  = msa->nseq;

  /* Clean up and exit. */
  for(i = 0; i < msa->nseq; i++)
    {
      free(map[i]);
      free(uaseq[i]);
    }
  free(map);
  free(uaseq);
  FreeParsetree(mtr);
  Free2DArray((void**)dsq, msa->nseq);
  MSAFree(msa);
  return;
}


/* Function: compare_cms()
 * EPN, Tue Mar  6 08:32:12 2007
 *
 * Purpose:  Given two CMs, cm1 and cm2, compare them, returning TRUE 
 *           iff they have the same guide tree (same node architecture).
 *
 * Args:     cm1          - covariance model number 1
 *           cm2          - covariance model number 2
 * 
 * Returns:  TRUE if CMs have same guide tree, FALSE otherwise
 */
static int compare_cms(CM_t *cm1, CM_t *cm2)
{
  int          nd; 
  if(cm1->nodes != cm2->nodes) return FALSE;
  for(nd = 0; nd < cm1->nodes; nd++)
    if(cm1->ndtype[nd] != cm2->ndtype[nd]) return FALSE;
  return TRUE;
}

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
#include "sre_stack.h"
#include "hmmband.h"         
#include "cm_postprob.h"
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
   --zeroinserts : set insert emission scores to 0\n\
   --sub         : build sub CM for columns b/t HMM predicted start/end points\n\
   --elsilent    : disallow local end (EL) emissions\n\
   --enfstart <n>: enforce MATL stretch starting at CM node <n>\n\
   --enfseq   <s>: enforce MATL stretch starting at --enfstart <n> emits seq <s>\n\
\n\
  * HMM banded alignment related options (IN DEVELOPMENT):\n\
   --hbanded     : use experimental CM plan 9 HMM banded CYK aln algorithm\n\
   --hbandp <f>  : tail loss prob for --hbanded [default: 0.0001]\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
   --hmmonly     : align with the CM Plan 9 HMM\n\
\n\
  * Query dependent banded (qdb) alignment related options:\n\
   --qdb         : use query dependent banded CYK alignment algorithm\n\
   --beta <f>    : tail loss prob for --qdb [default:0.0000001]\n\
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
  { "--hbandp",     FALSE, sqdARG_FLOAT},
  { "--sums",       FALSE, sqdARG_NONE},
  { "--time",       FALSE, sqdARG_NONE},
  { "--inside",     FALSE, sqdARG_NONE},
  { "--outside",    FALSE, sqdARG_NONE},
  { "--post",       FALSE, sqdARG_NONE},
  { "--checkpost",  FALSE, sqdARG_NONE},
  { "--zeroinserts",FALSE, sqdARG_NONE},
  { "--sub",        FALSE, sqdARG_NONE},
  { "--elsilent",   FALSE, sqdARG_NONE},
  { "--enfstart",   FALSE, sqdARG_INT},
  { "--enfseq",     FALSE, sqdARG_STRING},
  { "--zeroinserts",FALSE, sqdARG_NONE}
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
  MSA             *msa;         /* alignment that's created                 */
  char            *outfile;	/* optional output file name                */
  FILE            *ofp;         /* an open output file                      */
  int              i;		/* counter over sequences/parsetrees        */
  char            *regressfile; /* regression test data file                */
  char            *tracefile;	/* file to dump debugging traces to         */

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

  /* HMM banded alignment data structures */
  int              do_hbanded;  /* TRUE to use CP9 HMM to band CYK          */
  double           hbandp;      /* tail loss probability for hmm bands      */
  int              use_sums;    /* TRUE: use the posterior sums w/HMM bands */
  int              do_hmmonly;  /* TRUE: align with the HMM, not the CM     */
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
  int   enf_start;              /* first MATL node to enforce each parse use*/
  int   enf_end;                /* last  MATL node to enforce each parse use*/
  char *enf_seq;                /* the subsequence to enforce emitted by    *
                                 * in nodes from enf_start to enf_end       */

#ifdef USE_MPI
  int mpi_my_rank;              /* My rank in MPI                           */
  int mpi_num_procs;            /* Total number of processes                */
  int mpi_master_rank;          /* Rank of master process                   */
  Stopwatch_t  *mpi_watch;      /* for timings in MPI mode                  */

  /* Initailize MPI, get values for rank and num procs */
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
  hbandp      = DEFAULT_HBANDP;
  use_sums    = FALSE;
  do_timings  = FALSE;
  do_inside   = FALSE;
  do_outside  = FALSE;
  do_post     = FALSE;
  do_sub      = FALSE;
  do_fullsub  = FALSE;
  do_elsilent = FALSE;
  do_zero_inserts =FALSE;
  do_enforce   = FALSE;
  enf_start    = 0;
  enf_end      = 0;
  enf_seq      = NULL;

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
    else if (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--hmmonly")   == 0) do_hmmonly   = TRUE;
    else if (strcmp(optname, "--enfstart")  == 0) { do_enforce = TRUE; enf_start = atoi(optarg); }
    else if (strcmp(optname, "--enfseq")    == 0) { do_enforce = TRUE; enf_seq = optarg; } 
    else if (strcmp(optname, "--zeroinserts")== 0) do_zero_inserts = TRUE;
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

  if(do_sub && do_local)
    Die("--sub and -l combination not supported.\n");
  if(do_sub && do_qdb)
    Die("Please pick either --sub or --qdb.\n");
  if(do_hmmonly)
    Die("--hmmonly not yet implemented.\n");
  if(do_enforce && enf_seq == NULL)
    Die("--enfstart only makes sense with --enfseq also.\n");
  if(do_enforce && enf_start == 0)
    Die("--enfseq only makes sense with --enfstart (which can't be 0) also.\n");
  if(bdump_level > 0 && (!do_qdb))
    Die("--bdump option does not work with --noqdb.\n");
  if(do_qdb && do_hbanded)
    Die("--qdb and --hbanded combo not supported, pick one.\n");
  if (bdump_level > 3) 
    Die("Highest available --banddump verbosity level is 3\n%s", usage);
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
  cm->hbandp = hbandp;   /* this will be DEFAULT_HBANDP unless changed at command line */

  /* Update cm->config_opts and cm->align_opts based on command line options */
  if(do_local)        cm->config_opts |= CM_CONFIG_LOCAL;
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
  if(do_enforce)
    {
      cm->config_opts |= CM_CONFIG_ENFORCE;
      cm->enf_start = enf_start; 
      cm->enf_seq   = enf_seq;
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

#ifdef USE_MPI
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
      parallel_align_targets(seqfp, cm, &sq, &tr, &postcode, &nseq,
			     bdump_level, debug_level, be_quiet,
			     mpi_my_rank, mpi_master_rank, mpi_num_procs);
      printf("done parallel align_seqs.\n");
    }
  else
#endif /* (end of if USE_MPI) */
    {
      /* Configure the CM for alignment based on cm->config_opts and cm->align_opts.
       * set local mode, make cp9 HMM, calculate QD bands etc. */
      ConfigCM(cm, NULL, NULL);
      serial_align_targets(seqfp, cm, &sq, &tr, &postcode, &nseq, bdump_level, debug_level, 
			   be_quiet);
    }
#ifdef USE_MPI
  if (mpi_my_rank == mpi_master_rank) {
#endif  
    /***************************************************************                     
     * Create the MSA.                                                                   
     ****************************************************************/                   
    msa = NULL;                                                                          
    if(!((cm->align_opts & CM_ALIGN_INSIDE) || (cm->align_opts & CM_ALIGN_OUTSIDE)))   
      {                                                                                  
	msa = ESL_Parsetrees2Alignment(cm, sq, NULL, tr, nseq, do_full);                 
	if(cm->align_opts & CM_ALIGN_POST)                                              
        {                                                                              
	  if(postcode == NULL)
	    Die("ERROR CM_ALIGN_POST flag is up, but {serial,parallel}_align_targets() did not return post codes.\n");
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
	    printf("%-40s ... ", "Saving parsetrees"); fflush(stdout);
	    if ((ofp = fopen(tracefile,"w")) == NULL)
	      Die("failed to open trace file %s", tracefile);
	    for (i = 0; i < msa->nseq; i++) 
	      {
		fprintf(ofp, "> %s\n", msa->sqname[i]);
		fprintf(ofp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], sq[i]->dsq, FALSE));;
		ParsetreeDump(ofp, tr[i], cm, sq[i]->dsq);
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
	if(!(do_inside || do_outside)) FreeParsetree(tr[i]);
	esl_sq_Destroy(sq[i]);
      }
    esl_sqfile_Close(seqfp);

    if(!(do_inside || do_outside)) MSAFree(msa);
    free(sq);
    free(tr);
    SqdClean();

#ifdef USE_MPI
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
#ifdef USE_MPI
  printf("EXITING rank:%d\n", mpi_my_rank);
#endif
  return EXIT_SUCCESS;
}



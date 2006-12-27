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
   --fsub <f>    : sub CM w/structure b/t HMM start/end pts w/ > <f> prob mass\n\
   --elsilent    : disallow local end (EL) emissions\n\
   --enfstart <n>: enforce MATL stretch starting at CM node <n>\n\
   --enfseq   <s>: enforce MATL stretch starting at --enfstart <n> emits seq <s>\n\
\n\
  * HMM banded alignment related options (IN DEVELOPMENT):\n\
   --hbanded     : use experimental CM plan 9 HMM banded CYK aln algorithm\n\
   --hbandp <f>  : tail loss prob for --hbanded [default: 0.0001]\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
   --hmmonly     : align with the CM Plan 9 HMM\n\
   --checkcp9    : check the CP9 empirically by generating sequences\n\
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
  { "--sub",        FALSE, sqdARG_NONE},
  { "--fsub",       FALSE, sqdARG_FLOAT},
  { "--elsilent",   FALSE, sqdARG_NONE},
  { "--checkcp9",   FALSE, sqdARG_NONE},
  { "--enfstart",   FALSE, sqdARG_INT},
  { "--enfseq",     FALSE, sqdARG_STRING}
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
  char  *regressfile;           /* regression test data file */
  char  *tracefile;		/* file to dump debugging traces to        */
  int    be_quiet;		/* TRUE to suppress verbose output & banner */
  int    do_local;		/* TRUE to config the model in local mode   */
  int    do_small;		/* TRUE to do divide and conquer alignments */
  int    do_full;               /* TRUE to output all match columns in output alignment */
  int    do_qdb;                /* TRUE to do qdb CYK (either d&c or full)  */
  int    bdump_level;           /* verbosity level for --banddump option, 0 is OFF */
  char  *optname;               /* name of option found by Getopt()        */
  char  *optarg;                /* argument found by Getopt()              */
  int    optind;                /* index in argv[]                         */
  int    debug_level;           /* verbosity level for debugging printf() statements */
  int    do_timings;             /* TRUE to print timings, FALSE not to */
  double   qdb_beta;	        /* tail loss probability for query dependent banding */

  /* HMMERNAL!: hmm banded alignment data structures */
  int         do_hbanded;       /* TRUE to do CM Plan 9 HMM banded CYKInside_b_jd() using bands on d and j dim*/
  double      hbandp;           /* tail loss probability for hmm bands */
  int         use_sums;         /* TRUE to fill and use the posterior sums, false not to. */
  int         do_hmmonly;       /* TRUE to align with the HMM, not the CM */
  
  /* CM Plan 9 */
  int                    do_checkcp9;   /* TRUE to check the CP9 HMM by generating sequences */
  int                    seed;	        /* random number generator seed (only used if(do_checkcp9)) */

  /* Alternatives to CYK */
  int                do_inside; /* TRUE to use the Inside algorithm instead of CYK */
  int                do_outside;/* TRUE to use the Outside algorithm instead of CYK */
  int                do_check;  /* TRUE to check Inside and Outside probabilities */
  int                do_post;   /* TRUE to do a posterior decode instead of CYK */
  char             **postcode;  /* posterior decode array of strings        */
  char              *apostcode; /* aligned posterior decode array           */

  /* the --sub option */
  int                do_sub;       /* TRUE to use HMM to infer start and end point of each seq
				    * and build a separate sub CM for alignment of that seq */
  int                do_fullsub;   /* TRUE to only remove structure outside HMM predicted start
				    * (spos) and end points (epos) */
  float              fsub_pmass;   /* probability mass from HMM posteriors req'd of start before sstruct
				    * and end after estruct */
  CM_t              *orig_cm;      /* the original, template covariance model the sub CM was built from */
  char            *check_outfile;  /* output file name for subCM stk file*/
  
  /* Options for checking sub_cm construction */
  int do_atest;                 /* TRUE to build 2 ML HMMs, one from the CM and one from
				 * the sub_cm, analytically, and check to make sure
				 * the corresponding parameters of these two HMMS
				 * are within 'pthresh' of each other */
  float atest_pthresh;          /* Probability threshold for atest */

  int               do_elsilent;  /* TRUE to disallow EL emissions, by setting EL self transition prob
				   * to as close to IMPOSSIBLE as we can and avoid underflow errors */
  
  /* The enforce option (--enfstart and --enfseq), added specifically for enforcing the template 
   * region for telomerase RNA searches */
  int   do_enforce;             /* TRUE to read .enforce file and enforce MATL stretch */
  int   enf_start;              /* if (do_enforce), first MATL node to enforce each parse enter */
  int   enf_end;                /* if (do_enforce), last  MATL node to enforce each parse enter */
  char *enf_seq;                /* if (do_enforce), the subsequence to enforce in nodes from enf_start to 
				 * enf_end */
  int   nd;

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
  do_qdb = FALSE;
  qdb_beta    = 0.0000001;
  do_full     = FALSE;
  bdump_level = 0;
  debug_level = 0;
  do_hbanded  = FALSE;
  do_hmmonly  = FALSE;
  hbandp      = 0.0001;
  use_sums    = FALSE;
  do_timings  = FALSE;
  do_inside   = FALSE;
  do_outside  = FALSE;
  do_check    = FALSE;
  do_post     = FALSE;
  do_sub      = FALSE;
  do_fullsub  = FALSE;
  do_elsilent = FALSE;
  do_checkcp9 = FALSE;
  do_atest    = FALSE;
  atest_pthresh = 0.00001;
  check_outfile = "check.stk";
  seed         = time ((time_t *) NULL);
  fsub_pmass   = 0.;
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
    else if (strcmp(optname, "--outside")   == 0) { do_outside = TRUE; do_check = TRUE; }
    else if (strcmp(optname, "--post")      == 0) do_post      = TRUE;
    else if (strcmp(optname, "--checkpost") == 0) do_check     = TRUE;
    else if (strcmp(optname, "--sub")       == 0) do_sub       = TRUE; 
    else if (strcmp(optname, "--fsub")      == 0) 
      { do_sub = TRUE; do_fullsub = TRUE; fsub_pmass = atof(optarg); }
    else if (strcmp(optname, "--elsilent")  == 0) do_elsilent  = TRUE;
    else if (strcmp(optname, "--hbanded")   == 0) { do_hbanded = TRUE; do_small = FALSE; }
    else if (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--hmmonly")   == 0) do_hmmonly   = TRUE;
    else if (strcmp(optname, "--checkcp9")  == 0) do_checkcp9  = TRUE;
    else if (strcmp(optname, "--enfstart")  == 0) { do_enforce = TRUE; enf_start = atoi(optarg); }
    else if (strcmp(optname, "--enfseq")    == 0) { do_enforce = TRUE; enf_seq = optarg; } 
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
  if(do_sub && do_local && !do_fullsub)
    Die("--sub and -l combination not supported.\n");
  if(do_sub && do_qdb)
    Die("Please pick either --sub or --qdb.\n");
  if(do_hmmonly)
    Die("--hmmonly not yet implemented.\n");
  if(do_enforce && enf_seq == NULL)
    Die("--enfstart only makes sense with --enfseq also.\n");
  if(do_enforce && enf_start == 0)
    Die("--enfseq only makes sense with --enfstart (which can't be 0) also.\n");

  if (bdump_level > 3) Die("Highest available --banddump verbosity level is 3\n%s", usage);
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

  /* We need to ensure that cm->el_selfsc * W >= IMPOSSIBLE
   * (cm->el_selfsc is the score for an EL self transition) This is
   * done because we potentially multiply cm->el_selfsc * W, and add
   * that to IMPOSSIBLE. To avoid underflow issues this value must be
   * less than 3 * IMPOSSIBLE. Here, to be safe, we guarantee its less
   * than 2 * IMPOSSIBLE.
   */

  if((cm->el_selfsc * cm->W) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (cm->W+1));

  if (do_local && do_hbanded)
    {
      printf("Warning: banding with an HMM (--hbanded) and allowing\nlocal alignment (-l). This may not work very well.\n");
    } 

  CMLogoddsify(cm);
  /*CMHackInsertScores(cm);*/	/* "TEMPORARY" fix for bad priors */

  if(do_enforce)
    {
      enf_end = enf_start + strlen(enf_seq) - 1;
      EnforceSubsequence(cm, enf_start, enf_seq);
    }
  /*****************************************************************
   * Input and digitize the unaligned sequences
   *****************************************************************/

  if (! ReadMultipleRseqs(seqfile, format, &rseq, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);
  dsq = MallocOrDie(sizeof(char *) * nseq);
  for (i = 0; i < nseq; i++) 
    dsq[i] = DigitizeSequence(rseq[i], sqinfo[i].len);
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

  
  AlignSeqsWrapper(cm, dsq, sqinfo, nseq, &tr, do_local, do_small, do_qdb,
		   0.0000001, do_hbanded, use_sums, hbandp,
		   do_sub, do_fullsub, fsub_pmass, do_hmmonly, do_inside, do_outside, do_check, 
		   do_post, &postcode, do_timings, bdump_level, debug_level, FALSE,
		   do_enforce, enf_start, enf_end, do_elsilent, 
		   NULL, NULL, NULL, NULL, NULL, NULL); 
  /* last 6 NULL args are specific to partial-test.c */

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
      
      printf("\n");
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
  if(!(do_inside || do_outside)) MSAFree(msa);
  FreeCM(orig_cm);
  free(rseq);
  free(sqinfo);
  free(dsq);
  free(tr);
  
  SqdClean();
  return 0;
}

  

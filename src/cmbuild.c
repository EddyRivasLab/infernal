/* cmbuild.c
 * SRE, Thu Jul 27 13:19:43 2000 [StL]
 * CVS $Id$
 * 
 * Construct a CM from a given multiple sequence alignment.
 *  
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

static char banner[] = "cmbuild - build an RNA covariance model from an alignment";

static char usage[]  = "\
Usage: cmbuild [-options] <cmfile output> <alignment file>\n\
The alignment file is expected to be in Stockholm format.\n\
  Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -A     : append; append this CM to <cmfile>\n\
   -B     : Babelfish; autodetect alternative alignment file format\n\
   -F     : force; allow overwriting of <cmfile>\n\
";

static char experts[] = "\
   --rf          : use #=RF alignment annotation to specify consensus vs. insertion\n\
   --gapthresh   : fraction of gaps to allow in a consensus column (0..1)\n\
   --binary      : save output model (cmfile) in binary, not ASCII text\n\
   --informat <s>: specify that input alignment is in format <s>, not Stockholm\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-A", TRUE, sqdARG_NONE },
  { "-B", TRUE, sqdARG_NONE },
  { "-F", TRUE, sqdARG_NONE },
  { "--binary",    FALSE, sqdARG_NONE },
  { "--gapthresh", FALSE, sqdARG_FLOAT},
  { "--informat",  FALSE, sqdARG_STRING },
  { "--rf",        FALSE, sqdARG_NONE },

};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;	/* file to write CM to                     */
  char            *alifile;     /* seqfile to read alignment from          */ 
  int              format;	/* format of seqfile                       */
  MSAFILE         *afp;         /* open alignment file                     */
  FILE            *cmfp;	/* output file for CMs                     */
  MSA             *msa;         /* a multiple sequence alignment           */
  int              nali;	/* number of alignments processed          */
  int              idx;		/* counter over seqs                       */
  Parsetree_t     *mtr;         /* master structure tree from the alignment*/
  Parsetree_t    **tr;		/* inidividual traces from alignment       */
  CM_t            *cm;          /* a covariance model                      */
  Stopwatch_t     *watch;	/* timer to run  */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  int   do_append;		/* TRUE to append CM to cmfile             */
  int   do_binary;		/* TRUE to write CM in binary not ASCII    */
  int   allow_overwrite;	/* true to allow overwriting cmfile        */
  char  fpopts[3];		/* options to open a file with, e.g. "ab"  */
  int   use_rf;			/* TRUE to use #=RF annot to define consensus cols  */
  float gapthresh;		/* 0=all cols are inserts; 1=all cols are consensus */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = MSAFILE_STOCKHOLM;
  do_append         = FALSE;
  do_binary         = FALSE;
  allow_overwrite   = FALSE;
  gapthresh         = 0.5;
  use_rf            = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-A") == 0)          do_append         = TRUE; 
    else if (strcmp(optname, "-B") == 0)          format            = MSAFILE_UNKNOWN;
    else if (strcmp(optname, "-F") == 0)          allow_overwrite   = TRUE;
    else if (strcmp(optname, "--gapthresh") == 0) gapthresh         = atof(optarg);
    else if (strcmp(optname, "--rf")        == 0) use_rf            = TRUE;

    else if (strcmp(optname, "--informat") == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == MSAFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
      if (! IsAlignmentFormat(format))
	Die("%s is an unaligned format, can't read as an alignment", optarg);
    }
    else if (strcmp(optname, "-h") == 0) {
      Banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  alifile = argv[optind++]; 

  if (!allow_overwrite && !do_append && FileExists(cmfile))
    Die("CM file %s already exists. Rename or delete it.", cmfile); 

  /*********************************************** 
   * Preliminaries: open our files for i/o
   ***********************************************/

				/* Open the alignment */
  if ((afp = MSAFileOpen(alifile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", alifile);
  
				/* Open the CM output file */
  if (do_append) strcpy(fpopts, "a");
  else           strcpy(fpopts, "w");
  if (do_binary) strcat(fpopts, "b");
  if ((cmfp = fopen(cmfile, fpopts)) == NULL)
    Die("Failed to open CM file %s for %s\n", cmfile, 
	do_append ? "appending" : "writing");

  watch = StopwatchCreate();

  /*********************************************** 
   * Show the banner
   ***********************************************/
  
  Banner(stdout, banner);
  printf("Alignment file:                    %s\n", alifile);
  printf("File format:                       %s\n", 
	 SeqfileFormat2String(afp->format));
  printf("New CM file:                       %s %s\n",
	 cmfile, do_append? "[appending]" : "");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /*********************************************** 
   * Get alignment(s), build CMs one at a time
   ***********************************************/

  nali = 0;
  while ((msa = MSAFileRead(afp)) != NULL)
    {
      /* Print some stuff about what we're about to do.
       */
      if (msa->name != NULL) printf("Alignment:           %s\n",  msa->name);
      else                   printf("Alignment:           #%d\n", nali+1);
      printf                       ("Number of sequences: %d\n",  msa->nseq);
      printf                       ("Number of columns:   %d\n",  msa->alen);
      puts("");
      fflush(stdout);
      
      /* Make alignment upper case, because some symbol counting
       * things are case-sensitive.
       */
      for (idx = 0; idx < msa->nseq; idx++)
	s2upper(msa->aseq[idx]);

      /* Construct a model
       */
      HandModelmaker(msa, use_rf, gapthresh, &cm, &mtr, &tr);
      /* PrintParsetree(stdout, mtr);  */
      PrintCM(stdout, cm); 
      SummarizeMasterTrace(stdout, mtr); 
      SummarizeCM(stdout, cm); 
      
      CMSimpleProbify(cm);
      CMSetDefaultNullModel(cm);

o      CMLogoddsify(cm);
      for (idx = 0; idx < msa->nseq; idx++)
	{
	  char *rseq, *dsq;
	  int   L;
	  
	  MakeDealignedString(msa->aseq[idx], msa->alen, msa->aseq[idx], &rseq);
	  L = strlen(rseq);
	  dsq = DigitizeSequence(rseq, L);

#if 0
	  StopwatchZero(watch);
	  StopwatchStart(watch);
	  CYKInside(cm, dsq, L);
	  StopwatchStop(watch);
	  StopwatchDisplay(stdout, "CPU time: ", watch);
#endif
	  ParsetreeDump(stdout, tr[idx], cm, dsq+1);
	  printf("trace score says: %.2f\n",
		 ParsetreeScore(cm, tr[idx], msa->aseq[idx])/ 0.693);

	  free(rseq);
	  free(dsq);
	}


      WriteBinaryCM(cmfp, cm);

      FreeParsetree(mtr);
      FreeCM(cm);
      MSAFree(msa);
    }


  /* Clean up and exit
   */
  StopwatchFree(watch);
  MSAFileClose(afp);
  fclose(cmfp);
  SqdClean();
  return 0;
}

 

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
   --tfile <s>   : dump traces to debugging file <s>\n\
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
  { "--tfile",     FALSE, sqdARG_STRING },

};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void dump_traces(char *tracefile, MSA *msa, Parsetree_t **tr, CM_t *cm, char **dsq); 

int
main(int argc, char **argv)
{
  char            *cmfile;	/* file to write CM to                     */
  char            *alifile;     /* seqfile to read alignment from          */ 
  int              format;	/* format of seqfile                       */
  MSAFILE         *afp;         /* open alignment file                     */
  FILE            *cmfp;	/* output file for CMs                     */
  MSA             *msa;         /* a multiple sequence alignment           */
  char           **dsq;		/* digitized aligned sequences             */
  int              nali;	/* number of alignments processed          */
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
  char *tracefile;		/* file to dump debugging traces to */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = MSAFILE_STOCKHOLM;
  do_append         = FALSE;
  do_binary         = FALSE;
  allow_overwrite   = FALSE;
  gapthresh         = 0.5;
  use_rf            = FALSE;
  tracefile         = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-A") == 0)          do_append         = TRUE; 
    else if (strcmp(optname, "-B") == 0)          format            = MSAFILE_UNKNOWN;
    else if (strcmp(optname, "-F") == 0)          allow_overwrite   = TRUE;
    else if (strcmp(optname, "--gapthresh") == 0) gapthresh         = atof(optarg);
    else if (strcmp(optname, "--rf")        == 0) use_rf            = TRUE;
    else if (strcmp(optname, "--tfile")     == 0) tracefile         = optarg;

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
      
      /* Digitize the alignment: this takes care of
       * case sensivitity (A vs. a), speeds all future
       * array indexing, and deals with the poor fools
       * who would give us horrid DNA (T) instead of
       * lovely RNA (U). It does cause one wee problem:
       * you need to keep in mind that a digitized seq
       * is indexed 1..alen, but msa (and its annotation!!)
       * is indexed 0..alen-1.
       */
      dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);

      /* Construct a model
       */
      HandModelmaker(msa, dsq, use_rf, gapthresh, &cm, &mtr, &tr);
      /* PrintParsetree(stdout, mtr);  */
      PrintCM(stdout, cm); 
      SummarizeMasterTrace(stdout, mtr); 
      MasterTraceDisplay(stdout, mtr, cm);
      SummarizeCM(stdout, cm); 
      
      CMSimpleProbify(cm);
      CMSetDefaultNullModel(cm);

      CMLogoddsify(cm);

      if (tracefile != NULL) 
	dump_traces(tracefile, msa, tr, cm, dsq);

      WriteBinaryCM(cmfp, cm);

      FreeParsetree(mtr);
      FreeCM(cm);
      Free2DArray((void**)dsq, msa->nseq);
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


/* Function: dump_traces()
 * Date:     SRE, Sun Aug  6 09:36:32 2000 [St. Louis]
 *
 * Purpose:  for debugging information: dump individual 
 *           transmogrified tracebacks and trace scores
 *           to a file.
 */
static void
dump_traces(char *tracefile, MSA *msa, Parsetree_t **tr, CM_t *cm, char **dsq)
{
  FILE *fp;
  int   i;

  if ((fp = fopen(tracefile,"w")) == NULL)
    Die("failed to open trace file %s", tracefile);

  for (i = 0; i < msa->nseq; i++) {
    fprintf(fp, "> %s\n", msa->sqname[i]);
    fprintf(fp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], dsq[i]) / 0.693);;
    ParsetreeDump(fp, tr[i], cm, dsq[i]);
    fprintf(fp, "//\n");
  }
  fclose(fp);
}

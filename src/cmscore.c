/* cmscore.c
 * SRE, Thu Aug  3 17:08:45 2000 [StL]
 * CVS $Id$
 * 
 * Score a CM against unaligned sequence examples.
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

static char banner[] = "cmscore - score an RNA covariance model against sequences";

static char usage[]  = "\
Usage: cmscore [-options] <cmfile> <sequence file>\n\
The sequence file is expected to be in FASTA format.\n\
  Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -B     : Babelfish; autodetect alternative sequence file format\n\
";

static char experts[] = "\
   --informat <s>: specify that input alignment is in format <s>, not FASTA\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-B", TRUE, sqdARG_NONE },
  { "--informat",  FALSE, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  char            *seqfile;     /* file to read sequences from */
  int              format;      /* format of sequence file */
  FILE            *cmfp;        /* open CM file for reading */
  SQFILE	  *sqfp;        /* open seqfile for reading */
  CM_t            *cm;          /* a covariance model       */
  char            *seq;         /* RNA sequence */
  SQINFO           sqinfo;      /* optional info attached to seq */
  char            *dsq;         /* digitized RNA sequence */
  Stopwatch_t     *watch;
  float            sc1,  sc2;	/* score of a sequence */
  Parsetree_t     *tr1, *tr2;	/* a traceback */
  
  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = SQFILE_FASTA;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-B") == 0)          format            = SQFILE_UNKNOWN;
    else if (strcmp(optname, "--informat") == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
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
  seqfile = argv[optind++]; 

  
  /*********************************************** 
   * Preliminaries: open our files for i/o; get a CM
   ***********************************************/

  watch = StopwatchCreate();

  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);
  if ((cmfp = fopen(cmfile, "rb")) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);

  if (! ReadBinaryCM(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);

  CMLogoddsify(cm);

  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
    {
      CYKDemands(cm, sqinfo.len); 

      if (sqinfo.len == 0) continue; 	/* silently skip len 0 seqs */
      
      dsq = DigitizeSequence(seq, sqinfo.len);

#if 0
      StopwatchZero(watch);
      StopwatchStart(watch);
      sc1 = CYKInside(cm, dsq, sqinfo.len, &tr1);  
      ParsetreeDump(stdout, tr1, cm, dsq);
      printf("%-12s : %.2f  %.2f\n", sqinfo.name, sc1/0.693,
	     ParsetreeScore(cm, tr1, dsq)/0.693);
      StopwatchStop(watch);
      StopwatchDisplay(stdout, "CPU time: ", watch);
#endif

      StopwatchZero(watch);
      StopwatchStart(watch);
      sc2 = CYKDivideAndConquer(cm, dsq, sqinfo.len, &tr2);  
      ParsetreeDump(stdout, tr2, cm, dsq);
      printf("%-12s : %.2f  %.2f\n", sqinfo.name, sc2/0.693,
	     ParsetreeScore(cm, tr2, dsq)/0.693);
      /* ParsetreeCompare(tr1, tr2);  */
      StopwatchStop(watch);
      StopwatchDisplay(stdout, "CPU time: ", watch);

      FreeSequence(seq, &sqinfo);
      /*      FreeParsetree(tr1);  */
      FreeParsetree(tr2); 
      free(dsq);
    }

  FreeCM(cm);
  fclose(cmfp);
  SeqfileClose(sqfp);
  StopwatchFree(watch);
  SqdClean();
  return EXIT_SUCCESS;
}

/* cmsearch.c
 * SRE, Fri May  3 13:58:18 2002
 * CVS $Id$
 * 
 * Search sequences with a CM.
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

static char banner[] = "cmsearch - search a sequence database with an RNA covariance model";

static char usage[]  = "\
Usage: cmsearch [-options] <cmfile> <sequence file>\n\
The sequence file is expected to be in FASTA format.\n\
  Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -W <n> : set scanning window size to <n> (default: 200)\n\
";

static char experts[] = "\
   --informat <s>: specify that input alignment is in format <s>, not FASTA\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-W", TRUE, sqdARG_INT }, 
  { "--informat",   FALSE, sqdARG_STRING },
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
  int              windowlen;	/* maximum len of hit; scanning window size */
  int              i;

  int    nhits;			/* number of hits in a seq */
  int   *hiti;                  /* start positions of hits */
  int   *hitj;                  /* end positions of hits */
  float *hitsc;			/* scores of hits */
  
  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = SQFILE_UNKNOWN;
  windowlen         = 200;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if       (strcmp(optname, "-W")          == 0) windowlen = atoi(optarg);
    else if  (strcmp(optname, "--informat")  == 0) {
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
      if (sqinfo.len == 0) continue; 	/* silently skip len 0 seqs */
      dsq = DigitizeSequence(seq, sqinfo.len);

      StopwatchZero(watch);
      StopwatchStart(watch);
      CYKScan(cm, dsq, sqinfo.len, windowlen, 
	      &nhits, &hiti, &hitj, &hitsc);

      printf("sequence: %s\n", sqinfo.name);
      for (i = 0; i < nhits; i++)
	printf("hit %-4d: %6d %6d %8.2f bits\n", i, hiti[i], hitj[i], hitsc[i]/0.693);
	  
      free(hiti);
      free(hitj);
      free(hitsc);

      StopwatchStop(watch);
      StopwatchDisplay(stdout, "CPU time: ", watch);
      printf("memory:   %.2f MB\n\n", CYKScanRequires(cm, sqinfo.len, windowlen));
      FreeSequence(seq, &sqinfo);
      free(dsq);
    }

  FreeCM(cm);
  fclose(cmfp);
  SeqfileClose(sqfp);
  StopwatchFree(watch);
  SqdClean();
  return EXIT_SUCCESS;
}

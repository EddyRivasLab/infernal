/* check-emit.c
 * SRE, Tue Oct 14 08:04:30 2003 [St. Louis]
 * 
 * "Unit testing" of Infernal's emit.c.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 * CVS $Id$
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

static char usage[] = "\
Usage: emit-test [-options] <cmfile>\n\
  where test options are:\n\
  (none, at present)\n\
";

static struct opt_s OPTIONS[] = {
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  CMFILE          *cmfp;        /* open CM file for reading */
  CM_t            *cm;          /* a covariance model       */
  Parsetree_t     *tr;		/* sampled parsetree from the CM */
  char            *seq;		/* emitted sequence, alphabetic */
  char            *dsq;		/* emitted sequence, digital */
  int              N;		/* sequence length */
  
  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  /*********************************************** 
   * Parse command line
   ***********************************************/

 while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
 }

 if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
 cmfile = argv[optind++];
 
  /*********************************************** 
   * Preliminaries: get our CM
   ***********************************************/
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);

  EmitParsetree(cm, &tr, &seq, &dsq, &N);
  ParsetreeDump(stdout, tr, cm, dsq);

  FreeParsetree(tr);
  free(seq);
  free(dsq);
  FreeCM(cm);
  SqdClean();
  exit(0);
}




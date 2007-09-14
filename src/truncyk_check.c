/* truncyk_check.c
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */

static char banner[] = "truncyk_check - score RNA covariance model against sequences";

static char usage[]  = "\
Usage: truncyk_check [-options] <cmfile> <sequence file>\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input sequence file is in format <s>\n\
   --regress <f> : save regression test data to file <f>\n\
   --scoreonly   : for full CYK/inside stage, do only score, save memory\n\
   --smallonly   : do only d&c, don't do full CYK/inside\n\
   --stringent   : require the two parse trees to be identical\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "--informat",   FALSE, sqdARG_STRING },
  { "--regress",    FALSE, sqdARG_STRING },
  { "--scoreonly",  FALSE, sqdARG_NONE },
  { "--smallonly",  FALSE, sqdARG_NONE },
  { "--stringent",  FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  char            *seqfile;     /* file to read sequences from */
  int              format;      /* format of sequence file */
  CMFILE          *cmfp;        /* open CM file for reading */
  SQFILE	  *sqfp;        /* open seqfile for reading */
  CM_t            *cm;          /* a covariance model       */
  char            *seq;         /* RNA sequence */
  SQINFO           sqinfo;      /* optional info attached to seq */
  char            *dsq;         /* digitized RNA sequence */
  Stopwatch_t     *watch;
  float            sc1,  sc2;	/* score of a sequence */
  Parsetree_t     *tr1, *tr2;	/* a traceback */
  float            ptsc1, ptsc2; /* scores from interpreting parsetrees */
  int              v, model_len;
  float            bsc;
  
  int   do_local;		/* TRUE to align locally w.r.t. model       */
  int   do_scoreonly;		/* TRUE for score-only (small mem) full CYK */
  int   do_smallonly;		/* TRUE to do only d&c, not full CYK/inside */
  int   compare_stringently;	/* TRUE to demand identical parse trees     */
  char *regressfile;		/* name of regression data file to save     */
  FILE *regressfp;              /* open filehandle for writing regressions  */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format              = SQFILE_UNKNOWN;
  do_local            = TRUE;
  do_scoreonly        = FALSE;
  do_smallonly        = FALSE;
  regressfile         = NULL;
  compare_stringently = FALSE;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "--regress")   == 0) regressfile         = optarg;
    else if (strcmp(optname, "--smallonly") == 0) do_smallonly        = TRUE;
    else if (strcmp(optname, "--scoreonly") == 0) do_scoreonly        = TRUE;
    else if (strcmp(optname, "--stringent") == 0) compare_stringently = TRUE;
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

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  seqfile = argv[optind++]; 
  
  /*********************************************** 
   * Preliminaries: open our files for i/o; get a CM
   ***********************************************/

  watch = StopwatchCreate();

  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);

  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);

				/* open regression test data file */
  if (regressfile != NULL) {
    if ((regressfp = fopen(regressfile, "w")) == NULL)
      Die("Failed to open regression test file %s", regressfile);
  }

  
  if (do_local) ConfigLocal(cm, 0.5, 0.5);
  CMLogoddsify(cm);
  CMHackInsertScores(cm);	/* TEMPORARY: FIXME */

  /* EPN 11.18.05 Now that know what windowlen is, we need to ensure that
   * cm->el_selfsc * W >= IMPOSSIBLE (cm->el_selfsc is the score for an EL self transition)
   * This is done because we are potentially multiply cm->el_selfsc * W, and adding
   * that to IMPOSSIBLE. To avoid underflow issues this value must be less than
   * 3 * IMPOSSIBLE. Here we guarantee its less than 2 * IMPOSSIBLE (to be safe).
   */
  if((cm->el_selfsc * cm->W) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (cm->W+1));

  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
    {
      /* CYKDemands(cm, sqinfo.len);  */

      if (sqinfo.len == 0) continue; 	/* silently skip len 0 seqs */
      
      dsq = DigitizeSequence(seq, sqinfo.len);

      tr1 = tr2 = NULL;

      /* Length correction for Parsetree scores */
      model_len = 0;
      for ( v = 0; v < cm->M; v++ )
      {
         if      ( cm->stid[v] == MATP_MP ) model_len += 2;
         else if ( cm->stid[v] == MATL_ML ) model_len += 1;
         else if ( cm->stid[v] == MATR_MR ) model_len += 1;
      }
      /* 2.0 instead of 2 to force floating point division, not integer division */
      bsc = sreLOG2(2.0/(model_len*(model_len+1)));

      if (! do_smallonly) {
	printf("Full inside algorithm:\n");
	printf("----------------------\n");
	StopwatchZero(watch);
	StopwatchStart(watch);
	if (do_scoreonly) {
	  sc1 = TrCYK_Inside(cm, dsq, sqinfo.len, 0, 1, sqinfo.len, NULL);
	  printf("%-12s : %.2f\n", sqinfo.name, sc1);
	} else {
	  sc1 = TrCYK_Inside(cm, dsq, sqinfo.len, 0, 1, sqinfo.len, &tr1);  
	  ParsetreeDump(stdout, tr1, cm, dsq);
          ptsc1 = ParsetreeScore(cm, tr1, dsq, FALSE);
          ptsc1 += bsc;
	  printf("%-12s : %.2f  %.2f\n", sqinfo.name, sc1, ptsc1);
	}
	StopwatchStop(watch);
	StopwatchDisplay(stdout, "CPU time: ", watch);
	puts("");
      }

      printf("Divide and conquer algorithm:\n");
      printf("-------------------------------\n");
      StopwatchZero(watch);
      StopwatchStart(watch);
      sc2 = TrCYK_DnC(cm, dsq, sqinfo.len, 0, 1, sqinfo.len, &tr2);  
      ParsetreeDump(stdout, tr2, cm, dsq);
      ptsc2 = ParsetreeScore(cm, tr2, dsq, FALSE);
      ptsc2 += bsc;
      printf("%-12s : %.2f  %.2f\n", sqinfo.name, sc2, ptsc2);
      StopwatchStop(watch);
      StopwatchDisplay(stdout, "CPU time: ", watch);
      puts("");

      /* Test that the two solutions are identical; or if not identical,
       * at least they're alternative solutions w/ equal score.
       * If not, fail w/ non-zero exit status, so qc protocols
       * can catch the problem.
       */
      if (tr1 != NULL && fabs(sc1 - ptsc1) >= 0.01)
	Die("TrCYKInside score differs from its parse tree's score\n");
      if (tr2 != NULL && fabs(sc2 - ptsc2) >= 0.01)
	Die("TrCYKDivideAndConquer score differs from its parse tree's score\n");
      if (!do_smallonly && fabs(sc1 - sc2) >= 0.01) 
	Die("TrCYKInside score differs from TrCYKDivideAndConquer\n");
      if (tr1 != NULL && tr2 != NULL && 
	  compare_stringently && !ParsetreeCompare(tr1, tr2))
	Die("Parse trees for TrCYKInside and TrCYKDivideAndConquer differ\n");
      
      /* Save regression test data
       */
      if (regressfile != NULL) {
	if (tr1 != NULL) ParsetreeDump(regressfp, tr1, cm, dsq);
	if (tr2 != NULL) ParsetreeDump(regressfp, tr2, cm, dsq);
      }

      FreeSequence(seq, &sqinfo);
      if (tr1 != NULL) FreeParsetree(tr1);  
      if (tr2 != NULL) FreeParsetree(tr2); 
      free(dsq);
    }

  if (regressfile != NULL) fclose(regressfp);
  FreeCM(cm);
  CMFileClose(cmfp);
  SeqfileClose(sqfp);
  StopwatchFree(watch);
  SqdClean();
  return EXIT_SUCCESS;
}

/* cmscore.c
 * SRE, Thu Aug  3 17:08:45 2000 [StL]
 * SVN $Id$
 * 
 * Score a CM against unaligned sequence examples.
 * 
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
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

static char banner[] = "cmscore - score RNA covariance model against sequences";

static char usage[]  = "\
Usage: cmscore [-options] <cmfile> <sequence file>\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input sequence file is in format <s>\n\
   --local       : do local alignment (w.r.t. model)\n\
   --regress <f> : save regression test data to file <f>\n\
   --scoreonly   : for full CYK/inside stage, do only score, save memory\n\
   --smallonly   : do only d&c, don't do full CYK/inside\n\
   --stringent   : require the two parse trees to be identical\n\
   --qdb         : use query dependent banded CYK alignment algorithm\n\
   --beta <f>    : tail loss prob for --qdb [default:0.0000001]\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "--informat",   FALSE, sqdARG_STRING },
  { "--local",      FALSE, sqdARG_NONE },
  { "--regress",    FALSE, sqdARG_STRING },
  { "--scoreonly",  FALSE, sqdARG_NONE },
  { "--smallonly",  FALSE, sqdARG_NONE },
  { "--stringent",  FALSE, sqdARG_NONE },
  { "--qdb",        FALSE, sqdARG_NONE },
  { "--beta",       FALSE, sqdARG_FLOAT},
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
  
  int   do_local;		/* TRUE to align locally w.r.t. model       */
  int   do_scoreonly;		/* TRUE for score-only (small mem) full CYK */
  int   do_smallonly;		/* TRUE to do only d&c, not full CYK/inside */
  int   compare_stringently;	/* TRUE to demand identical parse trees     */
  char *regressfile;		/* name of regression data file to save     */
  FILE *regressfp;              /* open filehandle for writing regressions  */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  int      do_qdb;              /* TRUE to do qdb CYK (either d&c or full)  */
  double   qdb_beta;	        /* tail loss probability for query dependent banding */
  double **gamma;               /* P(subseq length = n) for each state v    */
  int     *dmin;                /* minimum d bound for state v, [0..v..M-1] */
  int     *dmax;                /* maximum d bound for state v, [0..v..M-1] */
  int      safe_windowlen;      /* initial windowlen (W) used for calculating bands */
  int      expand_flag;         /* TRUE if the dmin and dmax vectors have just been 
				 * expanded (in which case we want to recalculate them 
				 * before we align a new sequence), and FALSE if not*/
  /*********************************************** 
   * Parse command line
   ***********************************************/

  format              = SQFILE_UNKNOWN;
  do_local            = FALSE;
  do_scoreonly        = FALSE;
  do_smallonly        = FALSE;
  regressfile         = NULL;
  do_qdb              = FALSE;
  qdb_beta            = 0.0000001;
  compare_stringently = FALSE;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "--local")     == 0) do_local            = TRUE;
    else if (strcmp(optname, "--regress")   == 0) regressfile         = optarg;
    else if (strcmp(optname, "--smallonly") == 0) do_smallonly        = TRUE;
    else if (strcmp(optname, "--scoreonly") == 0) do_scoreonly        = TRUE;
    else if (strcmp(optname, "--stringent") == 0) compare_stringently = TRUE;
    else if (strcmp(optname, "--qdb")       == 0) do_qdb              = TRUE;
    else if (strcmp(optname, "--qdb")       == 0) do_qdb              = TRUE;
    else if (strcmp(optname, "--beta")      == 0) qdb_beta            = atof(optarg);
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
  /*CMHackInsertScores(cm);*/	/* TEMPORARY: FIXME */

  /* Now that know what windowlen is, we need to ensure that 
   * cm->el_selfsc * W >= IMPOSSIBLE (cm->el_selfsc is the score for an EL self transition)
   * This is done because we are potentially multiply cm->el_selfsc * W, and adding
   * that to IMPOSSIBLE. To avoid underflow issues this value must be less than
   * 3 * IMPOSSIBLE. Here, to be safe, we guarantee its less than 2 * IMPOSSIBLE.
   */
  if((cm->el_selfsc * cm->W) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (cm->W+1));

  /* set up the query dependent bands, this has to be done after the ConfigLocal() call */
  if(do_qdb)
    {
      safe_windowlen = cm->W * 2;
      while(!(BandCalculationEngine(cm, safe_windowlen, qdb_beta, 0, &dmin, &dmax, &gamma, do_local)))
	{
	  FreeBandDensities(cm, gamma);
	  free(dmin);
	  free(dmax);
	  safe_windowlen *= 2;
	  /*printf("ERROR BandCalculationEngine returned false, windowlen adjusted to %d\n", safe_windowlen);*/
	}
      expand_flag = FALSE;
    }
  else
    dmin = dmax = NULL;

  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
    {
      if (sqinfo.len == 0) continue; 	/* silently skip len 0 seqs */

      CYKDemands(cm, sqinfo.len, dmin, dmax);  /* dmin and dmax will be NULL if not in QDB mode */

      if(do_qdb)
	{
	  /*Check if we need to reset the query dependent bands b/c they're currently expanded. */
	  if(expand_flag)
	    {
	      FreeBandDensities(cm, gamma);
	      free(dmin);
	      free(dmax);
	      while(!(BandCalculationEngine(cm, safe_windowlen, qdb_beta, 0, &dmin, &dmax, &gamma, do_local)))
		{
		  FreeBandDensities(cm, gamma);
		  free(dmin);
		  free(dmax);
		  safe_windowlen *= 2;
		}
	      expand_flag = FALSE;
	    }
	  if((sqinfo.len < dmin[0]) || (sqinfo.len > dmax[0]))
	    {
	      /* the seq we're aligning is outside the root band, so we expand.*/
	      ExpandBands(cm, sqinfo.len, dmin, dmax);
	      printf("Expanded bands for seq : %s\n", sqinfo.name);
	      expand_flag = TRUE;
	    }
	}
      dsq = DigitizeSequence(seq, sqinfo.len);

      tr1 = tr2 = NULL;

      if (! do_smallonly) {
	printf("Full inside algorithm:\n");
	printf("----------------------\n");
	StopwatchZero(watch);
	StopwatchStart(watch);
	if (do_scoreonly) {
	  sc1 = CYKInsideScore(cm, dsq, sqinfo.len, 0, 1, sqinfo.len, 
			       dmin, dmax); /* dmin, dmax = NULL if we're not in QDB mode */
	  printf("%-12s : %.2f\n", sqinfo.name, sc1);
	} else {
	  sc1 = CYKInside(cm, dsq, sqinfo.len, 0, 1, sqinfo.len, &tr1, 
			  dmin, dmax); /* dmin, dmax = NULL if we're not in QDB mode */
	  ParsetreeDump(stdout, tr1, cm, dsq);
	  printf("%-12s : %.2f  %.2f\n", sqinfo.name, sc1,
		 ParsetreeScore(cm, tr1, dsq, FALSE));
	}
	StopwatchStop(watch);
	StopwatchDisplay(stdout, "CPU time: ", watch);
	puts("");
      }

      printf("Divide and conquer algorithm:\n");
      printf("-------------------------------\n");
      StopwatchZero(watch);
      StopwatchStart(watch);
      sc2 = CYKDivideAndConquer(cm, dsq, sqinfo.len, 0, 1, sqinfo.len, &tr2,
				dmin, dmax); /* dmin, dmax = NULL if we're not in QDB mode */
      ParsetreeDump(stdout, tr2, cm, dsq);
      printf("%-12s : %.2f  %.2f\n", sqinfo.name, sc2,
	     ParsetreeScore(cm, tr2, dsq, FALSE));
      StopwatchStop(watch);
      StopwatchDisplay(stdout, "CPU time: ", watch);
      puts("");

      /* Test that the two solutions are identical; or if not identical,
       * at least they're alternative solutions w/ equal score.
       * If not, fail w/ non-zero exit status, so qc protocols
       * can catch the problem.
       */
      if (tr1 != NULL && fabs(sc1 - ParsetreeScore(cm, tr1, dsq, FALSE)) >= 0.01)
	Die("CYKInside score differs from its parse tree's score");
      if (tr2 != NULL && fabs(sc2 - ParsetreeScore(cm, tr2, dsq, FALSE)) >= 0.01)
	Die("CYKDivideAndConquer score differs from its parse tree's score");
      if (!do_smallonly && fabs(sc1 - sc2) >= 0.01) 
	Die("CYKInside score differs from CYKDivideAndConquer");
      if (tr1 != NULL && tr2 != NULL && 
	  compare_stringently && !ParsetreeCompare(tr1, tr2))
	Die("Parse trees for CYKInside and CYKDivideAndConquer differ");
      
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

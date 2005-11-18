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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */

static void debug_print_bands(CM_t *cm, int *dmin, int *dmax);

static char banner[] = "cmsearch - search a sequence database with an RNA covariance model";

static char usage[]  = "\
Usage: cmsearch [-options] <cmfile> <sequence file>\n\
The sequence file is expected to be in FASTA format.\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -W <n> : set scanning window size to <n> (default: precalc'd in cmbuild)\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input alignment is in format <s>, not FASTA\n\
   --toponly     : only search the top strand\n\
   --local       : do local alignment\n\
   --noalign     : find start/stop only; don't do alignments\n\
   --dumptrees   : dump verbose parse tree information for each hit\n\
   --banded      : use experimental banded CYK scanning algorithm\n\
   --bandp <f>   : tail loss prob for --banded [df: 0.0001]\n\
   --banddump    : print bands for each state\n\
   --X           : project X!\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-W", TRUE, sqdARG_INT }, 
  { "--dumptrees",  FALSE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--local",      FALSE, sqdARG_NONE },
  { "--noalign",    FALSE, sqdARG_NONE },
  { "--toponly",    FALSE, sqdARG_NONE },
  { "--banded",     FALSE, sqdARG_NONE },
  { "--bandp",      FALSE, sqdARG_FLOAT},
  { "--banddump",   FALSE, sqdARG_NONE},
  { "--X",          FALSE, sqdARG_NONE },
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
  int              i;
  int              reversed;    /* TRUE when we're doing the reverse complement strand */
  int              maxlen;
  Parsetree_t     *tr;		/* parse of an individual hit */
  CMConsensus_t   *cons;	/* precalculated consensus info for display */
  Fancyali_t      *ali;         /* alignment, formatted for display */

  double  **gamma;              /* P(subseq length = n) for each state v    */
  int     *dmin;		/* minimum d bound for state v, [0..v..M-1] */
  int     *dmax; 		/* maximum d bound for state v, [0..v..M-1] */
  double   bandp;		/* tail loss probability for banding */

  int    nhits;			/* number of hits in a seq */
  int   *hitr;			/* initial states for hits */
  int   *hiti;                  /* start positions of hits */
  int   *hitj;                  /* end positions of hits */
  float *hitsc;			/* scores of hits */

  int    windowlen;		/* maximum len of hit; scanning window size */
  int    do_revcomp;		/* true to do reverse complement too */
  int    do_local;		/* TRUE to do local alignment */
  int    do_align;              /* TRUE to calculate and show alignments */
  int    do_dumptrees;		/* TRUE to dump parse trees */
  int    do_banded;		/* TRUE to do banded CYK */
  int    do_projectx;           /* TRUE to activate special in-progress testing code */
  int    do_bdump;              /* TRUE to print out bands */
  /*EPN 08.18.05*/
  int    set_window;            /* TRUE to set window length due to -W option*/

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */


  /*EPN 11.11.05 */
  int      safe_windowlen;	/* initial windowlen (W) used for calculating bands
				 * in BandCalculationEngine().
				 * this needs to be significantly bigger than what
				 * we expect dmax[0] to be, for truncation error
				 * handling.
				 */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = SQFILE_UNKNOWN;
  windowlen         = 200;
  set_window        = FALSE;
  do_revcomp        = TRUE;
  do_local          = FALSE;
  do_align          = TRUE;
  do_dumptrees      = FALSE;
  do_banded         = FALSE;
  bandp             = 0.0001;
  do_projectx       = FALSE;
  do_bdump          = FALSE;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if       (strcmp(optname, "-W")          == 0) 
      { windowlen    = atoi(optarg); set_window = TRUE; }
    else if  (strcmp(optname, "--dumptrees") == 0) do_dumptrees = TRUE;
    else if  (strcmp(optname, "--local")     == 0) do_local     = TRUE;
    else if  (strcmp(optname, "--noalign")   == 0) do_align     = FALSE;
    else if  (strcmp(optname, "--toponly")   == 0) do_revcomp   = FALSE;
    else if  (strcmp(optname, "--banded")    == 0) do_banded    = TRUE;
    else if  (strcmp(optname, "--bandp")     == 0) bandp        = atof(optarg);
    else if  (strcmp(optname, "--X")         == 0) do_projectx  = TRUE;
    else if  (strcmp(optname, "--banddump")  == 0) do_bdump     = TRUE;
    else if  (strcmp(optname, "--informat")  == 0) {
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

  /* EPN 08.18.05 */
  if (! (set_window)) windowlen = cm->W;
  /*printf("***cm->W : %d***\n", cm->W);*/

  if (do_local) ConfigLocal(cm, 0.5, 0.5);
  CMLogoddsify(cm);
  CMHackInsertScores(cm);	/* make insert emissions score zero. "TEMPORARY" FIX. */
  cons = CreateCMConsensus(cm, 3.0, 1.0); 

  if (do_banded || do_projectx || do_bdump)
    {
      safe_windowlen = windowlen * 2;
      while(!(BandCalculationEngine(cm, safe_windowlen, bandp, 0, &dmin, &dmax, &gamma, do_local)))
	{
	  FreeBandDensities(cm, gamma);
	  free(dmin);
	  free(dmax);
	  safe_windowlen *= 2;
	}

      /* EPN 11.11.05 
       * An important design decision.
       * We're changing the windowlen value here. By default,
       * windowlen is read from the cm file (set to cm->W). 
       * Here we're doing a banded scan though. Its pointless to allow
       * a windowlen that's greater than the largest possible banded hit 
       * (which is dmax[0]). So we reset windowlen to dmax[0].
       * Its also possible that BandCalculationEngine() returns a dmax[0] that 
       * is > cm->W. This should only happen if the bandp we're using now is < 1E-7 
       * (1E-7 is the bandp value used to determine cm->W in cmbuild). If this 
       * happens, the current implementation reassigns windowlen to this larger value.
       * NOTE: if W was set at the command line, the command line value is 
       *       always used.
       */
      if(!(set_window))
	{
	  windowlen = dmax[0];
	}
      if(do_bdump) 
	{
	  printf("bandp:%f\n", bandp);
	  debug_print_bands(cm, dmin, dmax);
	}
    }
  StopwatchZero(watch);
  StopwatchStart(watch);
  
  /* EPN 11.18.05 Now that know what windowlen is, we need to ensure that
   * cm->el_selfsc * W >= IMPOSSIBLE (cm->el_selfsc is the score for an EL self transition)
   * This is done because we are potentially multiply cm->el_selfsc * W, and adding
   * that to IMPOSSIBLE. To avoid underflow issues this value must be less than
   * 3 * IMPOSSIBLE. Here we guarantee its less than 2 * IMPOSSIBLE (to be safe).
   */
  if((cm->el_selfsc * windowlen) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (windowlen+1));

  maxlen   = 0;
  reversed = FALSE;
  while (reversed || ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
    {
      if (sqinfo.len == 0) continue; 	/* silently skip len 0 seqs */
      /* EPN 08.12.05
       * Without the next if(), a strange, unresolved segmentation 
       * fault causing memory error occurs.  This fix should be okay 
       * as long as we're not trying to do --local because no glocal
       * hits can be smaller than the minimum bound of the root band.  
       * Currently, the --banded --local scenario doesn't work (at 
       * least) for this reason.
       */
      if (do_banded && (sqinfo.len < dmin[0] || sqinfo.len > dmax[0]))
	{
	  /*printf("sequence: %s\n", sqinfo.name);*/
	  /*continue;*/
	}
      if (sqinfo.len > maxlen) maxlen = sqinfo.len;
      dsq = DigitizeSequence(seq, sqinfo.len);

      if (do_banded)
	  CYKBandedScan(cm, dsq, dmin, dmax, sqinfo.len, windowlen, 
			&nhits, &hitr, &hiti, &hitj, &hitsc);
      else
	CYKScan(cm, dsq, sqinfo.len, windowlen, 
		&nhits, &hitr, &hiti, &hitj, &hitsc);
      if (! reversed) printf("sequence: %s\n", sqinfo.name);
      for (i = 0; i < nhits; i++)
	{
	  printf("hit %-4d: %6d %6d %8.2f bits\n", i, 
		 reversed ? sqinfo.len - hiti[i] + 1 : hiti[i], 
		 reversed ? sqinfo.len - hitj[i] + 1 : hitj[i],
		 hitsc[i]);
	  if (do_align) 
	    {
	      CYKDivideAndConquer(cm, dsq, sqinfo.len, 
				  hitr[i], hiti[i], hitj[i], &tr);
	      ali = CreateFancyAli(tr, cm, cons, dsq);
	      PrintFancyAli(stdout, ali);
	      
	      if (do_dumptrees) {
		ParsetreeDump(stdout, tr, cm, dsq);
		printf("\tscore = %.2f\n", ParsetreeScore(cm,tr,dsq));
	      }
	      if (do_projectx) {
		BandedParsetreeDump(stdout, tr, cm, dsq, gamma, windowlen, dmin, dmax);
	      }
	      
	      FreeFancyAli(ali);
	      FreeParsetree(tr);
	    }
	}
      
      free(hitr);
      free(hiti);
      free(hitj);
      free(hitsc);
      
      free(dsq);
      if (! reversed && do_revcomp) {
	revcomp(seq,seq);
	reversed = TRUE;
      } else {
	reversed = FALSE;
	FreeSequence(seq, &sqinfo);
      }
    }
  StopwatchStop(watch);
  StopwatchDisplay(stdout, "\nCPU time: ", watch);
  printf("memory:   %.2f MB\n\n", CYKScanRequires(cm, maxlen, windowlen));

  if (do_banded)
    {
      FreeBandDensities(cm, gamma);
      free(dmin);
      free(dmax);
    }
  FreeCMConsensus(cons);
  FreeCM(cm);
  CMFileClose(cmfp);
  SeqfileClose(sqfp);
  StopwatchFree(watch);
  SqdClean();
  return EXIT_SUCCESS;
}

/* EPN 05.09.05
  debug_print_bands()
 * Function: debug_print_bands
 *
 * Purpose:  Print bands for each state.
 */

static void
debug_print_bands(CM_t *cm, int *dmin, int *dmax)
{
  int v;

  printf("\nPrinting bands :\n");
  printf("****************\n");
  for(v = 0; v < cm->M; v++)
   {
     printf("band state:%d type:%d min:%d max:%d\n", v, cm->sttype[v], dmin[v], dmax[v]);
   }
  printf("****************\n\n");
}


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
#include <ctype.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

static char banner[] = "cmbuild - build RNA covariance model from alignment";

static char usage[]  = "\
Usage: cmbuild [-options] <cmfile output> <alignment file>\n\
The alignment file is expected to be in Stockholm format.\n\
  Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -A     : append; append this CM to <cmfile>\n\
   -F     : force; allow overwriting of <cmfile>\n\
";

static char experts[] = "\
   --nobalance   : don't rebalance the CM; number in strict preorder\n\
   --rf          : use #=RF alignment annotation to specify consensus\n\
   --gapthresh   : fraction of gaps to allow in a consensus column (0..1)\n\
   --informat <s>: specify input alignment is in format <s>, not Stockholm\n\
\n\
 Verbose output files, useful for detailed information about the CM:\n\
   --cfile <f>   : save count vectors to file <f>\n\
   --gtree <f>   : save tree description of master tree to file <f>\n\
   --gtbl  <f>   : save tabular description of master tree to file <f>\n\
   --cmtbl <f>   : save tabular description of CM topology to file <f>\n\
   --tfile <f>   : dump individual sequence tracebacks to file <f>\n\

";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-A", TRUE, sqdARG_NONE },
  { "-F", TRUE, sqdARG_NONE },
  { "--nobalance", FALSE, sqdARG_NONE },
  { "--cfile",     FALSE, sqdARG_STRING },
  { "--cmtbl",     FALSE, sqdARG_STRING },
  { "--gapthresh", FALSE, sqdARG_FLOAT},
  { "--gtbl",      FALSE, sqdARG_STRING },
  { "--gtree",     FALSE, sqdARG_STRING },
  { "--informat",  FALSE, sqdARG_STRING },
  { "--rf",        FALSE, sqdARG_NONE },
  { "--tfile",     FALSE, sqdARG_STRING },

};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void dump_traces(char *tracefile, MSA *msa, Parsetree_t **tr, 
			CM_t *cm, char **dsq); 
static int  save_countvectors(char *cfile, CM_t *cm);
static int  clean_cs(char *cs, int alen);

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
  Stopwatch_t     *watch;	/* timer to run                            */
  int              idx;         /* sequence index                          */
  int              avlen;	/* average sequence length in MSA          */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  int   do_append;		/* TRUE to append CM to cmfile             */
  int   do_balance;		/* TRUE to balance the CM                  */
  int   do_binary;		/* TRUE to use binary file format for CM   */
  int   allow_overwrite;	/* true to allow overwriting cmfile        */
  char  fpopts[3];		/* options to open a file with, e.g. "ab"  */
  int   use_rf;			/* TRUE to use #=RF to define consensus    */
  float gapthresh;		/* 0=all cols inserts; 1=all cols consensus*/

  FILE *ofp;                    /* filehandle to dump info to */
  char *cfile;                  /* file to dump count vectors to */
  char *tracefile;		/* file to dump debugging traces to        */
  char *cmtblfile;              /* file to dump CM tabular description to  */
  char *gtreefile;              /* file to dump guide tree to              */
  char *gtblfile;               /* file to dump guide tree table to        */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = MSAFILE_UNKNOWN;   /* autodetect by default */
  do_append         = FALSE;
  do_balance        = TRUE;
  do_binary         = TRUE;
  allow_overwrite   = FALSE;
  gapthresh         = 0.5;
  use_rf            = FALSE;

  cfile             = NULL;
  tracefile         = NULL;
  cmtblfile         = NULL;
  gtblfile          = NULL;
  gtreefile         = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-A") == 0)          do_append         = TRUE; 
    else if (strcmp(optname, "-B") == 0)          format            = MSAFILE_UNKNOWN;
    else if (strcmp(optname, "-F") == 0)          allow_overwrite   = TRUE;
    else if (strcmp(optname, "--balance")   == 0) do_balance        = TRUE;
    else if (strcmp(optname, "--gapthresh") == 0) gapthresh         = atof(optarg);
    else if (strcmp(optname, "--rf")        == 0) use_rf            = TRUE;

    else if (strcmp(optname, "--cfile")     == 0) cfile             = optarg;
    else if (strcmp(optname, "--gtbl")      == 0) gtblfile          = optarg;
    else if (strcmp(optname, "--gtree")     == 0) gtreefile         = optarg;
    else if (strcmp(optname, "--cmtbl")     == 0) cmtblfile         = optarg;
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

  /* If we're outputting tranmogrified traces, shut CM rebalancing
   * off; else, our numbering of states in traces will not be consistent with
   * numbering of states in CM.
   */
  if (tracefile != NULL)  do_balance = FALSE;

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
      avlen = (int) MSAAverageSequenceLength(msa);

      /* Print some stuff about what we're about to do.
       */
      if (msa->name != NULL) printf("Alignment:           %s\n",  msa->name);
      else                   printf("Alignment:           #%d\n", nali+1);
      printf                       ("Number of sequences: %d\n",  msa->nseq);
      printf                       ("Number of columns:   %d\n",  msa->alen);
      printf                       ("Average seq length:  %d\n",  avlen);
      puts("");
      fflush(stdout);
      
      /* Some input data cleaning. 
       */
      if (use_rf && msa->rf == NULL) 
	Die("Alignment has no reference coord annotation.");
      if (msa->ss_cons == NULL) 
	Die("Alignment has no consensus structure annotation.");
      if (! clean_cs(msa->ss_cons, msa->alen))
	Die("Failed to parse consensus structure annotation.");
      printf("%-40s ... [done]\n", "Alignment format checks"); fflush(stdout);

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
      printf("%-40s ... [done]\n", "Digitizing alignment"); fflush(stdout);

      /* Construct a model.
       */
      HandModelmaker(msa, dsq, use_rf, gapthresh, &cm, &mtr, &tr);

      /* Save a count vector file, if asked.
       * Used primarily for making data files for training priors.
       */
      if (cfile != NULL) {
	printf("%-40s ... ", "Saving count vector file");
	fflush(stdout);
	if (! save_countvectors(cfile, cm)) printf("[FAILED]\n");
	else                                printf("done. [%s]\n", cfile);
      }

      CMSimpleProbify(cm);
      CMSetDefaultNullModel(cm);
      CMLogoddsify(cm);

      /* Rebalance the model (default).
       * This is shut off by --nobalance option, or if traces
       * are being saved to a file using --tfile option.
       */
      if (do_balance) 
	{
	  CM_t *new;
	  new = CMRebalance(cm);
	  FreeCM(cm);
	  cm = new;
	}

      /* Dump optional information to files:
       */
      /* 1. Tabular description of CM topology */
      if (cmtblfile != NULL) 
	{
	  if ((ofp = fopen(cmtblfile, "w")) == NULL) 
	    Die("Failed to open cm table file %s", cmtblfile);
	  PrintCM(ofp, cm); 	  
	  fclose(ofp);
	}

      /* 2. Tabular description of guide tree topology */
      if (gtblfile != NULL) 
	{
	  if ((ofp = fopen(gtblfile, "w")) == NULL) 
	    Die("Failed to open guide tree table file %s", gtblfile);
	  PrintParsetree(ofp, mtr);  
	  fclose(ofp);
	}

      /* 3. Tree description of guide tree topology */
      if (gtreefile != NULL) 
	{
	  if ((ofp = fopen(gtreefile, "w")) == NULL) 
	    Die("Failed to open guide tree file %s", gtreefile);
	  MasterTraceDisplay(ofp, mtr, cm);
	  fclose(ofp);
	}

      /* SummarizeMasterTrace(stdout, mtr); */

      SummarizeCM(stdout, cm); 
      puts("");
      CYKDemands(cm, avlen);
      puts("");
      
		/* note: if do_balance is true, the indices on these
                   tracebacks won't match the model. */
      if (tracefile != NULL)       
	dump_traces(tracefile, msa, tr, cm, dsq);

      WriteBinaryCM(cmfp, cm);

      for (idx = 0; idx < msa->nseq; idx++)
	FreeParsetree(tr[idx]);
      free(tr);

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
    fprintf(fp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], dsq[i]));;
    ParsetreeDump(fp, tr[i], cm, dsq[i]);
    fprintf(fp, "//\n");
  }
  fclose(fp);
}


/* Function: save_countvectors()
 * Date:     SRE, Tue May  7 16:21:10 2002 [St. Louis]
 *
 * Purpose:  Save emission count vectors to a file.
 *           Used to gather data for training Dirichlet priors.
 *
 * Args:     cfile  - name of file to save vectors to.
 *           cm     - a model containing counts (before probify'ing)
 *
 */
static int
save_countvectors(char *cfile, CM_t *cm)
{
  FILE *fp;
  int   v,x;

  if ((fp = fopen(cfile, "w")) == NULL) return 0;
  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] == MP_st || 
	  cm->sttype[v] == ML_st || 
	  cm->sttype[v] == MR_st) 
	{
	  fprintf(fp, "%-7s ", UniqueStatetype(cm->stid[v]));
	  if (cm->sttype[v] == MP_st) {
	    for (x = 0; x < Alphabet_size*Alphabet_size; x++)
	      fprintf(fp, "%8.3f ", cm->e[v][x]);
	  } else {
	    for (x = 0; x < Alphabet_size; x++)
	      fprintf(fp, "%8.3f ", cm->e[v][x]);
	  }
	  fprintf(fp, "\n");
	}
    }
  fclose(fp);
  return 1;
}


/* Functions: clean_cs()
 * Date:      SRE, Fri May 17 14:52:42 2002 [St. Louis]
 *
 * Purpose:   Verify and (if needed) clean the consensus structure annotation.
 */
static int
clean_cs(char *cs, int alen)
{
  int   i;
  int  *ct;
  int   status;
  int   nright = 0;
  int   nleft = 0;
  int   nbad = 0;
  char  example;
  int   first;
  int   has_pseudoknots = FALSE;

  /* 1. Maybe we're ok and don't need any cleaning.
   */
  status = KHS2ct(cs, alen, FALSE, &ct);
  free(ct);
  if (status == 1) return 1;

  /* 2. Maybe we have a good CS line but it annotates one or
   *    or more pseudoknots that have to be deleted.
   */
  if ((status = KHS2ct(cs, alen, TRUE, &ct)) == 1) { 
    has_pseudoknots = TRUE; 
    printf("    [Consensus structure has annotated pseudoknots that will be ignored.]\n");
    fflush(stdout);
  }
  free(ct);

  /* 3. Delete everything we don't recognize.
   */
  for (i = 0; i < alen; i++)
    {
      if      (isgap(cs[i])) ;
      else if (cs[i] == '>') nleft++;  
      else if (cs[i] == '<') nright++; 
      else if (has_pseudoknots && isalpha(cs[i])) cs[i] = '.';
      else {	/* count bad chars; remember first one; replace w/gap */
	if (nbad == 0) { example = cs[i]; first = i; }
	nbad++;
	cs[i] = '.';
      }
    }
  printf("    [Removed %d bad chars from consensus line. Example: a %c at position %d.]\n",
	 nbad, example, first);
  fflush(stdout);

  /* Check it again.
   */
  status = KHS2ct(cs, alen, FALSE, &ct);
  free(ct);
  if (status == 1) return 1;

  printf("    [Failed to parse the consensus structure line. %d >, %d <]\n", nleft, nright);
  return 0;
}

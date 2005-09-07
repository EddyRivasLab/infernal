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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */
#include "vectorops.h"

#include "structs.h"		/* data structures, macros, #define's   */
#include "prior.h"		/* mixture Dirichlet prior */
#include "funcs.h"		/* external functions                   */

static char banner[] = "cmbuild - build RNA covariance model from alignment";

static char usage[]  = "\
Usage: cmbuild [-options] <cmfile output> <alignment file>\n\
The alignment file is expected to be in Stockholm format.\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -n <s> : name this CM <s>\n\
   -A     : append; append this CM to <cmfile>\n\
   -F     : force; allow overwriting of <cmfile>\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --binary       : save the model in binary format\n\
   --rf           : use reference coordinate annotation to specify consensus\n\
   --gapthresh <x>: fraction of gaps to allow in a consensus column (0..1)\n\
   --informat <s> : specify input alignment is in format <s>, not Stockholm\n\
\n\
 * sequence weighting options (default: GSC weighting):\n\
   --wgiven       : use weights as annotated in alignment file\n\
   --wnone        : no weighting; re-set all weights to 1.\n\
   --wgsc         : use Gerstein/Sonnhammer/Chothia tree weights (default)\n\
\n\
 * verbose output files, useful for detailed information about the CM:\n\
   --cfile <f>    : save count vectors to file <f>\n\
   --cmtbl <f>    : save tabular description of CM topology to file <f>\n\
   --emap  <f>    : save consensus emit map to file <f>\n\
   --gtree <f>    : save tree description of master tree to file <f>\n\
   --gtbl  <f>    : save tabular description of master tree to file <f>\n\
   --tfile <f>    : dump individual sequence tracebacks to file <f>\n\
\n\
 * debugging, experimentation:\n\
   --nobalance    : don't rebalance the CM; number in strict preorder\n\
   --regress <f>  : save regression test information to file <f>\n\
   --treeforce    : score first seq in alignment and show parsetree\n\
   --ignorant     : strip the structural info from input alignment\n\
\n\
 * customize the Dirichlet priors:\n\
   --priorfile <f> : read priors from file <f>\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-n", TRUE, sqdARG_STRING },
  { "-A", TRUE, sqdARG_NONE },
  { "-F", TRUE, sqdARG_NONE },
  { "--binary",    FALSE, sqdARG_NONE },
  { "--nobalance", FALSE, sqdARG_NONE },
  { "--cfile",     FALSE, sqdARG_STRING },
  { "--cmtbl",     FALSE, sqdARG_STRING },
  { "--emap",      FALSE, sqdARG_STRING },
  { "--gapthresh", FALSE, sqdARG_FLOAT},
  { "--gtbl",      FALSE, sqdARG_STRING },
  { "--gtree",     FALSE, sqdARG_STRING },
  { "--informat",  FALSE, sqdARG_STRING },
  { "--regress",   FALSE, sqdARG_STRING },
  { "--rf",        FALSE, sqdARG_NONE },
  { "--tfile",     FALSE, sqdARG_STRING },
  { "--treeforce", FALSE, sqdARG_NONE },
  { "--wgiven",    FALSE, sqdARG_NONE },
  { "--wnone",     FALSE, sqdARG_NONE },
  { "--wgsc",      FALSE, sqdARG_NONE },
  { "--priorfile", FALSE, sqdARG_STRING },
  { "--ignorant",  FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static int  save_countvectors(char *cfile, CM_t *cm);
static int  clean_cs(char *cs, int alen);
/* EPN 08.18.05 */
static int MSAMaxSequenceLength(MSA *msa);

/* EPN 09.07.05 */
static void StripWUSS(char *ss);

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
  Parsetree_t     *tr;		/* individual traces from alignment        */
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
  int   treeforce;		/* number of seqs to show parsetrees for   */
  char *setname;                /* name to give to HMM, overriding others  */
  enum  weight_flags {		/* sequence weighting strategy to use      */
    WGT_GIVEN, WGT_NONE, WGT_GSC } weight_strategy;

  FILE *ofp;                    /* filehandle to dump info to */
  char *cfile;                  /* file to dump count vectors to */
  char *emapfile;		/* file to dump emit map to */
  char *tracefile;		/* file to dump debugging traces to        */
  char *cmtblfile;              /* file to dump CM tabular description to  */
  char *gtreefile;              /* file to dump guide tree to              */
  char *gtblfile;               /* file to dump guide tree table to        */
  char *regressionfile;		/* file to dump regression test info to    */
  FILE *regressfp;		/* open file to dump regression test info  */

  /*ADDED EPN 01.31.05*/
  char *prifile;                /* file with prior data */
  Prior_t *pri;                 /* mixture Dirichlet prior structure */

  /* Added EPN 08.18.05 So we can calculate an appropriate W when building the CM */
  double  **gamma;		/* cumulative distribution p(len <= n) for state v */
  int     *dmin;		/* minimum d bound for state v, [0..v..M-1] */
  int     *dmax; 		/* maximum d bound for state v, [0..v..M-1] */
  double   bandp;		/* tail loss probability for banding */
  int      safe_windowlen;	/* initial windowlen (W) used for calculating bands,
				 * once bands are calculated we'll set cm->W to 
				 * dmax[0] */
  /* Added EPN 09.07.05 So we can force ignorant (all single-stranded) models if we
     so choose*/
  int      be_ignorant;         /* TRUE to strip all bp information from the input
				 * consensus structure.
				 */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = MSAFILE_UNKNOWN;   /* autodetect by default */
  do_append         = FALSE;
  do_balance        = TRUE;
  do_binary         = FALSE;	/* default: save CMs in ASCII flatfile format */
  allow_overwrite   = FALSE;
  gapthresh         = 0.5;
  use_rf            = FALSE;
  treeforce         = 0;
  weight_strategy   = WGT_GSC;	/* default: GSC sequence weighting */
  setname           = NULL;	/* default: get CM name from alifile, or filename */
  pri               = NULL;

  cfile             = NULL;
  emapfile          = NULL;
  tracefile         = NULL;
  cmtblfile         = NULL;
  gtblfile          = NULL;
  gtreefile         = NULL;
  regressionfile    = NULL;
  prifile           = NULL;
  
  /* EPN 08.18.05 bandp used in BandBounds() for calculating W(W=dmax[0]) 
   *              Brief expt on effect of different bandp values in 
   *              ~nawrocki/notebook/5_0818_inf_cmbuild_bands/00LOG
   */
  /* bandp             = 0.0001; */
  bandp             = 0.0000001; 
  safe_windowlen    = 0;
  
  be_ignorant       = FALSE;	/* default: leave in bp information */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-A") == 0)          do_append         = TRUE; 
    else if (strcmp(optname, "-B") == 0)          format            = MSAFILE_UNKNOWN;
    else if (strcmp(optname, "-F") == 0)          allow_overwrite   = TRUE;
    else if (strcmp(optname, "-n") == 0)          setname           = optarg;
    else if (strcmp(optname, "--binary")    == 0) do_binary         = TRUE;
    else if (strcmp(optname, "--nobalance") == 0) do_balance        = FALSE;
    else if (strcmp(optname, "--gapthresh") == 0) gapthresh         = atof(optarg);
    else if (strcmp(optname, "--rf")        == 0) use_rf            = TRUE;
    else if (strcmp(optname, "--treeforce") == 0) treeforce         = 1;
    else if (strcmp(optname, "--wgiven")    == 0) weight_strategy   = WGT_GIVEN;
    else if (strcmp(optname, "--wgsc")      == 0) weight_strategy   = WGT_GSC;
    else if (strcmp(optname, "--wnone")     == 0) weight_strategy   = WGT_NONE;
    else if (strcmp(optname, "--cfile")     == 0) cfile             = optarg;
    else if (strcmp(optname, "--emap")      == 0) emapfile          = optarg;
    else if (strcmp(optname, "--gtbl")      == 0) gtblfile          = optarg;
    else if (strcmp(optname, "--gtree")     == 0) gtreefile         = optarg;
    else if (strcmp(optname, "--cmtbl")     == 0) cmtblfile         = optarg;
    else if (strcmp(optname, "--tfile")     == 0) tracefile         = optarg;
    else if (strcmp(optname, "--regress")   == 0) regressionfile    = optarg;
    else if (strcmp(optname, "--priorfile") == 0) prifile         = optarg;
    else if (strcmp(optname, "--ignorant")  == 0) be_ignorant       = TRUE;

    else if (strcmp(optname, "--informat") == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == MSAFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
      if (! IsAlignmentFormat(format))
	Die("%s is an unaligned format, can't read as an alignment", optarg);
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
  if ((cmfp = fopen(cmfile, fpopts)) == NULL)
    Die("Failed to open CM file %s for %s\n", cmfile, 
	do_append ? "appending" : "writing");

				/* open regression test file */
  if (regressionfile != NULL) {
    if ((regressfp = fopen(regressionfile, "w")) == NULL) 
      Die("Failed to open regression test file %s", regressionfile);
  }

  if (prifile != NULL)
    {
      FILE *pfp;
      if ((pfp = fopen(prifile, "r")) == NULL)
	Die("Failed to open prior file %s\n", prifile);
      if ((pri = Prior_Read(pfp)) == NULL)
	Die("Failed to parse prior file %s\n", prifile);
      fclose(pfp);
    }
  else 
    pri = Prior_Default();

  watch = StopwatchCreate();

  /*********************************************** 
   * Show the banner
   ***********************************************/
  
  MainBanner(stdout, banner);
  printf("Alignment file:                    %s\n", alifile);
  printf("File format:                       %s\n", 
	 SeqfileFormat2String(afp->format));
  printf("Model construction strategy:       %s\n",
	 (use_rf)? "Manual, from RF annotation" : "Fast/ad-hoc");
  printf("Sequence weighting strategy:       ");
  switch (weight_strategy) {
  case WGT_GIVEN: puts("use annotation in alifile, if any"); break;
  case WGT_NONE:  puts("no weights"); break;
  case WGT_GSC:   puts("GSC tree weights"); break;
  }
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
      printf("%-40s ... ", "Alignment format checks"); fflush(stdout);
      if (use_rf && msa->rf == NULL) 
	Die("failed... Alignment has no reference coord annotation.");
      if (msa->ss_cons == NULL) 
	Die("failed... Alignment has no consensus structure annotation.");
      if (! clean_cs(msa->ss_cons, msa->alen))
	Die("failed... Failed to parse consensus structure annotation.");
      printf("done.\n");

      /* If user inputted the --ignorant option, strip the consensus
       * structure of all base pair information.
       */
      if (be_ignorant)
	StripWUSS(msa->ss_cons);

      /* Sequence weighting. Default: GSC weights. If WGT_GIVEN,
       * do nothing.
       */
      if (weight_strategy == WGT_NONE) 
	{
	  FSet(msa->wgt, msa->nseq, 1.0);
	}
      else if (weight_strategy == WGT_GSC)
	{
	  printf("%-40s ... ", "Weighting sequences by GSC rule");
	  fflush(stdout);
	  GSCWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
	  printf("done.\n");
	}

      /* Digitize the alignment: this takes care of
       * case sensivitity (A vs. a), speeds all future
       * array indexing, and deals with the poor fools
       * who would give us horrid DNA (T) instead of
       * lovely RNA (U). It does cause one wee problem:
       * you need to keep in mind that a digitized seq
       * is indexed 1..alen, but msa (and its annotation!!)
       * is indexed 0..alen-1.
       */
      printf("%-40s ... ", "Digitizing alignment"); fflush(stdout);
      dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);
      printf("done.\n");

      /* Construct a model, and collect observed counts.
       * Note on "treeforce": this is the number of sequences that we
       *   will ignore for the purposes of count-collection and parameterization
       *   of the CM. We will only use this for debugging. These seqs
       *   are then dumped as full parsetrees to stdout.
       */
      printf("%-40s ... ", "Constructing model architecture");
      fflush(stdout);
      HandModelmaker(msa, dsq, use_rf, gapthresh, &cm, &mtr);
      if (do_balance) 
	{
	  CM_t *new;
	  new = CMRebalance(cm);
	  FreeCM(cm);
	  cm = new;
	}
      for (idx = treeforce; idx < msa->nseq; idx++)
	{
	  tr = Transmogrify(cm, mtr, dsq[idx], msa->aseq[idx], msa->alen);
	  ParsetreeCount(cm, tr, dsq[idx], msa->wgt[idx]);
	  FreeParsetree(tr);
	}
      printf("done.\n");

      /* Before converting to probabilities,
       * save a count vector file, if asked.
       * Used primarily for making data files for training priors.
       */
      if (cfile != NULL) {
	printf("%-40s ... ", "Saving count vector file"); fflush(stdout);
	if (! save_countvectors(cfile, cm)) printf("[FAILED]\n");
	else                                printf("done. [%s]\n", cfile);
      }

      /* Convert to probabilities, and the global log-odds form
       * we save the model in.
       */
      printf("%-40s ... ", "Converting counts to probabilities"); fflush(stdout);

      PriorifyCM(cm, pri);
      CMSetDefaultNullModel(cm);
      CMLogoddsify(cm);
      printf("done.\n");

      /* EPN 08.18.05
       * Calculate W for the model based on the seed seqs' tracebacks.
       *              Brief expt on effect of different bandp and 
       *              safe_windowlen scaling factor values in 
       *              ~nawrocki/notebook/5_0818_inf_cmbuild_bands/00LOG
       */
      printf("%-40s ... ", "Calculating windowlen for model"); fflush(stdout);
      safe_windowlen = 2 * MSAMaxSequenceLength(msa);
      /*printf("safe_windowlen : %d\n", safe_windowlen);*/
      /*printf("bandp          : %12f\n", bandp);*/
      gamma = BandDistribution(cm, safe_windowlen);
      BandBounds(gamma, cm->M, safe_windowlen, bandp, &dmin, &dmax);
      cm->W = dmax[0];
      /*printf("cm->W : %d\n", cm->W);*/
      printf("done.\n");

      /* Give the model a name (mandatory in the CM file).
       * Order of precedence:
       *      1. -n option  (only if a single alignment in file)
       *      2. msa->name  (only in Stockholm or SELEX files)
       *      3. filename, without tail (e.g. "rnaseP.msa" becomes "rnaseP")
       * Also, add any optional annotations.     
       */
      printf("%-40s ... ", "Naming and annotating model"); fflush(stdout);
      if (nali == 0)
	{
	  if      (setname != NULL)   cm->name = Strdup(setname);
	  else if (msa->name != NULL) cm->name = Strdup(msa->name);
	  else                        cm->name = FileTail(alifile, TRUE);
	}
      else
	{
	  if (setname != NULL)
	    Die("FAILED.\nOops. Wait. You can't use -n w/ an alignment database");
	  else if (msa->name != NULL)
	    cm->name = Strdup(msa->name);
	  else
	    Die("FAILED.\nOops. Wait. I need a name annotated in each alignment");
	}
      if (msa->acc  != NULL) cm->acc  = Strdup(msa->acc);
      if (msa->desc != NULL) cm->desc = Strdup(msa->desc);
      printf("done.\n");

      /* Save the CM. 
       */
      printf("%-40s ... ", "Saving model to file"); fflush(stdout);
      CMFileWrite(cmfp, cm, do_binary);
      printf("done. [%s]\n", cmfile);

      /* Dump optional information to files:
       */
      /* Tabular description of CM topology */
      if (cmtblfile != NULL) 
	{
	  printf("%-40s ... ", "Saving CM topology table"); fflush(stdout);
	  if ((ofp = fopen(cmtblfile, "w")) == NULL) 
	    Die("Failed to open cm table file %s", cmtblfile);
	  PrintCM(ofp, cm); 	  
	  fclose(ofp);
	  printf("done. [%s]\n", cmtblfile);
	}

      /* Tabular description of guide tree topology */
      if (gtblfile != NULL) 
	{
	  printf("%-40s ... ", "Saving guide tree table"); fflush(stdout);
	  if ((ofp = fopen(gtblfile, "w")) == NULL) 
	    Die("Failed to open guide tree table file %s", gtblfile);
	  PrintParsetree(ofp, mtr);  
	  fclose(ofp);
	  printf("done. [%s]\n", gtblfile);
	}

      /* Emit map.
       */
      if (emapfile != NULL) 
	{
	  CMEmitMap_t *emap;

	  printf("%-40s ... ", "Saving emit map"); fflush(stdout);
	  if ((ofp = fopen(emapfile, "w")) == NULL) 
	    Die("Failed to open emit map file %s", emapfile);
	  emap = CreateEmitMap(cm);
	  DumpEmitMap(ofp, emap, cm);
	  FreeEmitMap(emap);
	  fclose(ofp);
	  printf("done. [%s]\n", emapfile);
	}

      /* Tree description of guide tree topology */
      if (gtreefile != NULL) 
	{
	  printf("%-40s ... ", "Saving guide tree dendrogram"); fflush(stdout);
	  if ((ofp = fopen(gtreefile, "w")) == NULL) 
	    Die("Failed to open guide tree file %s", gtreefile);
	  MasterTraceDisplay(ofp, mtr, cm);
	  fclose(ofp);
	  printf("done. [%s]\n", gtreefile);
	}

      /* Detailed traces for the training set.
       */
      if (tracefile != NULL)       
	{
	  printf("%-40s ... ", "Saving parsetrees"); fflush(stdout);
	  if ((ofp = fopen(tracefile,"w")) == NULL)
	    Die("failed to open trace file %s", tracefile);
	  for (idx = treeforce; idx < msa->nseq; idx++) 
	    {
	      tr = Transmogrify(cm, mtr, dsq[idx], msa->aseq[idx], msa->alen);
	      fprintf(ofp, "> %s\n", msa->sqname[idx]);
	      fprintf(ofp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr, dsq[idx]));;
	      ParsetreeDump(ofp, tr, cm, dsq[idx]);
	      fprintf(ofp, "//\n");
	      FreeParsetree(tr);
	    }
	  fclose(ofp);
	  printf("done. [%s]\n", tracefile);
	}

      /* Regression test info.
       */
      if (regressionfile != NULL) {
	printf("%-40s ... ", "Saving regression test data"); fflush(stdout);
	SummarizeCM(regressfp, cm);
	PrintCM(regressfp, cm);
	PrintParsetree(regressfp, mtr);
	MasterTraceDisplay(regressfp, mtr, cm);
	for (idx = treeforce; idx < msa->nseq; idx++) 
	  {
	    tr = Transmogrify(cm, mtr, dsq[idx], msa->aseq[idx], msa->alen);
	    fprintf(regressfp, "> %s\n", msa->sqname[idx]);
	    fprintf(regressfp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr, dsq[idx]));
	    ParsetreeDump(regressfp, tr, cm, dsq[idx]);
	    fprintf(regressfp, "//\n"); 
	  }
	  printf("done. [%s]\n", regressionfile);
      }

      /* Detailed parsetrees for the test set of forced parsetrees.
       * We reconfig the model into local alignment.
       */
      if (treeforce) 
	{
	  ConfigLocal(cm, 0.5, 0.5);	  
	  CMLogoddsify(cm);
	  for (idx = 0; idx < treeforce; idx++) 
	    {
	      tr = Transmogrify(cm, mtr, dsq[idx], msa->aseq[idx], msa->alen);
	      printf("> %s\n", msa->sqname[idx]);
	      printf("  SCORE : %.2f bits\n", ParsetreeScore(cm, tr, dsq[idx]));
	      ParsetreeDump(stdout, tr, cm, dsq[idx]);
	      printf("//\n");
	      FreeParsetree(tr);
	    }
	}

      puts("");
      SummarizeCM(stdout, cm);  
      puts("");
      CYKDemands(cm, avlen);     

      FreeParsetree(mtr);
      FreeCM(cm);
      Free2DArray((void**)dsq, msa->nseq);
      MSAFree(msa);
      fflush(cmfp);
      puts("//\n");
      nali++;
    }


  /* Clean up and exit
   */
  if (regressionfile != NULL) fclose(regressfp);
  StopwatchFree(watch);
  MSAFileClose(afp);
  Prior_Destroy(pri);
  fclose(cmfp);
  SqdClean();

  /* EPN 08.18.05 */
  DMX2Free(gamma);
  free(dmin);
  free(dmax);

  return 0;
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

  /* Print emission counts */
  if ((fp = fopen(cfile, "w")) == NULL) return 0;
  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] == MP_st || 
	  cm->sttype[v] == ML_st || 
	  cm->sttype[v] == MR_st) 
	{
	  fprintf(fp, "E\t%-7s ", UniqueStatetype(cm->stid[v]));
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

  /* Print transition counts */
  for (v = 0; v < cm->M; v++)
    {
      if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	{
	  fprintf(fp, "T\t%-7s : %-2d", UniqueStatetype(cm->stid[v]), cm->ndtype[(cm->ndidx[v] + 1)]);
	  for (x = 0; x < cm->cnum[v]; x++)
	    {
	      fprintf(fp, "%8.3f ", cm->t[v][x]);
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
  status = WUSS2ct(cs, alen, FALSE, &ct);
  free(ct);
  if (status == 1) return 1;

  /* 2. Maybe we have a good CS line but it annotates one or
   *    or more pseudoknots that have to be deleted.
   */
  if ((status = WUSS2ct(cs, alen, TRUE, &ct)) == 1) { 
    has_pseudoknots = TRUE; 
    printf("    [Consensus structure has annotated pseudoknots that will be ignored.]\n");
    fflush(stdout);
  }
  free(ct);

  /* 3. Delete everything we don't recognize.
   */
  for (i = 0; i < alen; i++)
    {
      if      (strchr("{[(<", cs[i]) != NULL) nleft++;  
      else if (strchr(">)]}", cs[i]) != NULL) nright++; 
      else if (strchr(":_-,.~", cs[i]) != NULL) ;
      else if (has_pseudoknots && isalpha((int) cs[i])) cs[i] = '.';
      else {	/* count bad chars; remember first one; replace w/gap */
	if (nbad == 0) { example = cs[i]; first = i; }
	nbad++;
	cs[i] = '.';
      }
    }
  if (nbad > 0) {
    printf("    [Removed %d bad chars from consensus line. Example: a %c at position %d.]\n",
	   nbad, example, first);
    fflush(stdout);
  }

  /* Check it again.
   */
  status = WUSS2ct(cs, alen, FALSE, &ct);
  free(ct);
  if (status == 1) return 1;

  printf("    [Failed to parse the consensus structure line.]\n");
  return 0;
}

/* EPN 08.18.05
 * This function really belongs in msa.c in squid (or easel I guess) 
 * but was placed here to minimize both number of modified files and
 * confusion.
 */

/* Function: MSAMaxSequenceLength()
 * based on Function: MSAAverageSequenceLength()
 * (comments below from MSAAverageSequenceLength())
 *
 * Date:     SRE, Sat Apr  6 09:41:34 2002 [St. Louis]
 *
 * Purpose:  Return the maximum length of the (unaligned) sequences
 *           in the MSA.
 *
 * Args:     msa  - the alignment
 *
 * Returns:  average length
 */
int 
MSAMaxSequenceLength(MSA *msa)
{
  int   i;
  int max;
  
  max = 0;
  for (i = 0; i < msa->nseq; i++) 
    max = MAX(DealignedLength(msa->aseq[i]), max);

  return max;
}


/* Function:  StripWUSS()
 * EPN 09.07.05
 *
 * Purpose:   Strips a secondary structure string in WUSS notation 
 *            of all base pair information, resulting in a completely single 
 *            stranded structure: the secondary structure string is modified.
 *            
 *            Characters <([{  are converted to :   (left base of base pairs)
 *            Characters >)]}  are converted to :   (right base of base pairs)
 *            Characters _-,   are converted to :   (unpaired bases)
 *            Characters  .:~  are untouched        
 *            Pseudoknot characters are converted to : as well.
 *
 * Args:      ss - the string to convert
 *
 * Returns:   (void)
 */
void
StripWUSS(char *ss)
{
  char *s;

  for (s = ss; *s != '\0'; s++)
      if ((*s != '~') && (*s != '.')) *s = ':';
  return;
}

/* cmbuild.c
 * SRE, Thu Jul 27 13:19:43 2000 [StL]
 * SVN $Id$
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
#include "cm_eweight.h"         /* LSJ's entropy weighted functions ported
				 * from HMMER2.4devl [by EPN] */
#include "hmmer_funcs.h"
#include "hmmer_structs.h"
#include "hmmband.h"

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
   --bandp <x>    : tail loss prob for calc'ing W using bands [default: 1E-7]\n\
   --elself <x>   : set EL self transition prob to <x> [df: 0.94]\n\
   --nodetach     : do not 'detach' one of two inserts that model same column\n\ 
\n\
 * sequence weighting options [default: GSC weighting]:\n\
   --wgiven       : use weights as annotated in alignment file\n\
   --wnone        : no weighting; re-set all weights to 1.\n\
   --wgsc         : use Gerstein/Sonnhammer/Chothia tree weights [default]\n\
\n\
 * alternative effective sequence number strategies:\n\
   --effent       : entropy loss target [default]\n\
   --eloss <x>    : for --effent: set target loss to <x> [default: 0.54]\n\
   --effnone      : effective sequence number is just # of seqs\n\
\n\
 * verbose output files, useful for detailed information about the CM:\n\
   --cfile <f>    : save count vectors to file <f>\n\
   --cmtbl <f>    : save tabular description of CM topology to file <f>\n\
   --emap  <f>    : save consensus emit map to file <f>\n\
   --gtree <f>    : save tree description of master tree to file <f>\n\
   --gtbl  <f>    : save tabular description of master tree to file <f>\n\
   --tfile <f>    : dump individual sequence tracebacks to file <f>\n\
   --bfile <f>    : save band data to file <f>\n\
\n\
 * debugging, experimentation:\n\
   --nobalance    : don't rebalance the CM; number in strict preorder\n\
   --regress <f>  : save regression test information to file <f>\n\
   --treeforce    : score first seq in alignment and show parsetree\n\
   --ignorant     : strip the structural info from input alignment\n\
   --dlev <x>     : set verbosity of debugging print statements to <x> [1..3]\n\
\n\
 * customization of null model and priors:\n\
   --null      <f>: read null (random sequence) model from file <f>\n\
   --priorfile <f>: read priors from file <f>\n\
\n\
 * build an HMM\n\
   model architecture types: \n\
   --p7g <f>      : build glocal (NW) plan 7 HMMER 2.4 HMM to file <f>\n\
   --p7l <f>      : build  local (SW) plan 7 HMMER 2.4 HMM to file <f>\n\
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
  { "--bfile",     FALSE, sqdARG_STRING },
  { "--bandp",     FALSE, sqdARG_FLOAT},
  { "--ignorant",  FALSE, sqdARG_NONE },
  { "--null",      FALSE, sqdARG_STRING },
  { "--effnone",   FALSE, sqdARG_NONE },
  { "--effent",    FALSE, sqdARG_NONE },
  /*{ "--effrelent",  FALSE, sqdARG_NONE },*/
  { "--eloss",     FALSE, sqdARG_FLOAT },
  { "--elself",    FALSE, sqdARG_FLOAT},
  { "--dlev",      FALSE, sqdARG_INT },
  { "--p7g",       FALSE, sqdARG_STRING},
  { "--p7l",       FALSE, sqdARG_STRING},
  { "--nodetach",  FALSE, sqdARG_NONE},
  { "--noprior",   FALSE, sqdARG_NONE}
}
;
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static int  save_countvectors(char *cfile, CM_t *cm);
static int  clean_cs(char *cs, int alen);
static int  MSAMaxSequenceLength(MSA *msa);
static void model_trace_info_dump(FILE *ofp, CM_t *cm, Parsetree_t *tr, 
				  char *aseq);
static void PrintBandDensity(FILE *fp, double **gamma, int v, int W,
			     int min, int max);
static void StripWUSS(char *ss);

int
main(int argc, char **argv)
{
  char            *cmfile;	/* file to write CM to                     */
  char            *hmmfile;     /* file to write HMM to                    */
  char            *alifile;     /* seqfile to read alignment from          */ 
  int              format;	/* format of seqfile                       */
  MSAFILE         *afp;         /* open alignment file                     */
  FILE            *cmfp;	/* output file for CMs                     */
  MSA             *msa;         /* a multiple sequence alignment           */
  char           **dsq;		/* digitized aligned sequences             */
  unsigned char  **p7dsq;	/* digitized aligned sequences             */
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

  enum {			/* Effective sequence number strategy:        */
    EFF_NOTSETYET,
    EFF_NONE, 			/* --effnone: eff_nseq is nseq                */
    EFF_ENTROPY                 /* --effent:  entropy loss target             */
    /*EFF_RELENTROPY               --effrelent:  relative entropy loss target */
  } eff_strategy;

  FILE *ofp;                    /* filehandle to dump info to */
  char *cfile;                  /* file to dump count vectors to */
  char *emapfile;		/* file to dump emit map to */
  char *tracefile;		/* file to dump debugging traces to        */
  char *cmtblfile;              /* file to dump CM tabular description to  */
  char *gtreefile;              /* file to dump guide tree to              */
  char *gtblfile;               /* file to dump guide tree table to        */
  char *regressionfile;		/* file to dump regression test info to    */
  FILE *regressfp;		/* open file to dump regression test info  */

  /* EPN modifications (beginning 01.31.05) */
  char       *prifile;          /* file with prior data */
  Prior_t    *pri;              /* mixture Dirichlet prior structure */

  /* EPN 08.18.05 So we can calculate an appropriate W when building the CM */
  double   **gamma;             /* P(subseq length = n) for each state v    */
  int       *dmin;		/* minimum d bound for state v, [0..v..M-1] */
  int       *dmax; 		/* maximum d bound for state v, [0..v..M-1] */
  double     bandp;		/* tail loss probability for banding        */
  int        safe_windowlen;	/* initial windowlen (W) used for calculating bands,
				 * once bands are calculated we'll set cm->W to dmax[0] */
  char      *bandfile;
  int        v;                 /* counter over states                      */
  int        be_ignorant;       /* TRUE to strip all bp info from the input structure */
  char      *rndfile;		/* random sequence model file to read       */

  /* EPN 11.07.05 entropy weighting */
  float      eloss;		/* target entropy loss, entropy-weights    */
  int        eloss_set;		/* TRUE if eloss was set on commandline    */
  float      etarget;		/* target entropy (background - eloss)     */
  float      eff_nseq;		/* effective sequence number               */
  int        eff_nseq_set;	/* TRUE if eff_nseq has been calculated    */
  float      randomseq[MAXABET];/* null sequence model                     */
  int        do_local;		/* always FALSE for now, used in BandCalculationEngine().
				 * Does cmbuild have a --local option in its future? */
  int        save_gamma;	/* TRUE to save the gamma matrix in
				 * BandCalculationEngine().*/
  /* EPN 11.15.05 User has the option to set the EL self transition probability
   * at the command line. This value is converted to a score (log2(prob)) and
   * stored in the CM file. By default this probability is 1.0, which has score
   * 0, so by default this self transition probability doesn't affect the performance
   * of the model at all relative to a v0.6 or earlier model.*/
  float      el_selfprob;      /* EL state's self transition probability. This is really hacky.
				* EL states are different, they don't have a transition score vector,
				* This probability is converted to a score and used to penalize very 
				* large regions 'skipped' by EL states. 
				* see ~nawrocki/notebook/5_1115_inf_el_trans_prob/00LOG for details.
				*/
  /* EPN 10.10.05 HMMERNAL (HMM banded alignment) additions */
  struct plan7_s    *p7hmm;        /* constructed p7 HMM; written to hmmfile  */
  struct p7prior_s  *p7pri;        /* Dirichlet priors to use                 */
  struct p7trace_s **p7tr;         /* fake tracebacks for aseq's              */
  int                checksum;     /* checksum of the alignment               */
  float              p1;           /* null sequence model p1 transition       */
  char              *name;         /* name of the HMM                         */
  enum p7_config {                 /* algorithm configuration strategy        */
    P7_BASE_CONFIG, P7_LS_CONFIG, P7_FS_CONFIG, P7_SW_CONFIG } cfg_strategy;
  float              swentry;      /* S/W aggregate entry probability         */
  float              swexit;       /* S/W aggregate exit probability          */
  FILE              *hmmfp;        /* HMM output file handle                  */
  int                debug_level;  /* level of debugging print statements     */
  int                build_p7hmm;  /* TRUE to build a plan 7 HMM              */
		      
  /* variables for 'detach' mode: we find columns that are modelled by two insert states (an ambiguity
   * in the CM architecture), and *handle* them (right hand in fist pounding open left palm) by 
   * picking the one that is never parameterized by any counts in the seed alignment (only the prior) 
   * and detaching it from the rest of the model by setting all transitions into it as 0.0: 
   * a remarkably ad hoc approach
   */
  int                do_detach;   /* TRUE to 'detach' off one of two insert states that 
				   * insert at the same position. */
  int                no_prior;    /* TRUE to not use a prior */

  /* Do hmmbuild.c stuff */
  p7pri = P7DefaultInfernalPrior();
  /* Set up the null/random seq model */
  P7DefaultNullModel(randomseq, &p1);
  cfg_strategy      = P7_BASE_CONFIG;  /* Our default is global (NW) style */
  swentry           = 0.5;
  swexit            = 0.5;
  Alphabet_type     = hmmNOTSETYET;	/* initially unknown */   

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
  bandfile          = NULL;
  rndfile           = NULL;
  
  eff_strategy      = EFF_ENTROPY; /* 11.28.05 EPN entropy weighting is default. */
  eloss_set         = FALSE;
  eff_nseq_set      = FALSE;
  do_local          = FALSE;
  save_gamma        = FALSE;
  eloss             = 0.54;  /* EPN: empirically determined optimal eloss using RMARK benchmark*/ 

  bandp             = 0.0000001; 
  safe_windowlen    = 0;
  be_ignorant       = FALSE; /* default: leave in bp information */
  el_selfprob       = 0.94;  /* EPN: empirically determined optimal 
			      * EL self transition prob using RMARK benchmark 
			      * (11.28.05) */
  debug_level       = 0; 
  do_detach         = TRUE;  /* default: detach 1 of 2 ambiguous inserts */
  no_prior          = FALSE; /* default: use a prior */
  build_p7hmm       = TRUE;  /* default: don't build a p7 HMM */
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
    else if (strcmp(optname, "--priorfile") == 0) prifile           = optarg;
    else if (strcmp(optname, "--bfile")     == 0) {bandfile         = optarg; save_gamma = TRUE;}
    else if (strcmp(optname, "--bandp")     == 0) bandp             = atof(optarg);
    else if (strcmp(optname, "--ignorant")  == 0) be_ignorant       = TRUE;
    else if (strcmp(optname, "--null")      == 0) rndfile           = optarg;
    else if (strcmp(optname, "--nodetach")  == 0) do_detach         = FALSE;   
    else if (strcmp(optname, "--noprior")   == 0) no_prior          = TRUE;   
    else if (strcmp(optname, "--effent")    == 0) eff_strategy      = EFF_ENTROPY;
    /*else if (strcmp(optname, "--effrelent")  == 0) eff_strategy   = EFF_RELENTROPY;*/
    else if (strcmp(optname, "--effnone")   == 0) eff_strategy      = EFF_NONE;
    else if (strcmp(optname, "--eloss")     == 0) { eloss  = atof(optarg); eloss_set  = TRUE; }

    else if (strcmp(optname, "--elself")    == 0) { 
      el_selfprob= atof(optarg); 
      if(el_selfprob < 0 || el_selfprob > 1)
	Die("EL self transition probability must be between 0 and 1.\n");
    }
    else if (strcmp(optname, "--dlev")      == 0)  debug_level  = atoi(optarg);
    else if (strcmp(optname, "--p7g")        == 0) { 
      if(build_p7hmm) Die("Can't do --p7g or --p7l, pick one HMM type.\n");
      build_p7hmm = TRUE; hmmfile = optarg; 
      cfg_strategy = P7_BASE_CONFIG;  /* NW style */
    }
    else if (strcmp(optname, "--p7l")        == 0) { 
      if(build_p7hmm) Die("Can't do --p7g or --p7l, pick one HMM type.\n");
      build_p7hmm = TRUE; hmmfile = optarg; 
      cfg_strategy = P7_SW_CONFIG;  /* SW style */
    }
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

  /* Open the HMM file */
  if(build_p7hmm)
    {
      if ((hmmfp = fopen(hmmfile, fpopts)) == NULL)
	Die("Failed to open HMM file %s for %s\n", hmmfile,
	    do_append ? "appending" : "writing");
    }
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
  printf("Null model used:                   %s\n",
	 (rndfile == NULL) ? "(default)" : rndfile);
  printf("Prior used:                        %s\n",
	 (prifile == NULL) ? "(default)" : prifile);
  printf("Effective sequence # calculation:  ");
  if (eff_strategy == EFF_NONE)      puts("none; use actual seq #");
  else if (eff_strategy == EFF_ENTROPY) {
    puts("entropy targeting");
    if (eloss_set)
      printf("  mean target entropy loss:        %.2f bits\n", eloss);
    else
      printf("  mean target entropy loss:        %.2f bits [default]\n", eloss);
  }
  
  /*EPN 11.28.05
   *Uncomment code below to use relative entropy weighting (identical 
   *to entropy weighting for equiprobable null distribution.)
   *
   * else if (eff_strategy == EFF_RELENTROPY) {
   * puts("relative entropy targeting");
   * if (eloss_set)
   *  printf("  mean target relative entropy loss:  %.2f bits\n", eloss);
   * else
   * printf("  mean target relative entropy loss:  %.2f bits [default]\n", eloss);
   * }
   */

  printf("Sequence weighting strategy:       ");
  switch (weight_strategy) {
  case WGT_GIVEN: puts("use annotation in alifile, if any"); break;
  case WGT_NONE:  puts("no weights"); break;
  case WGT_GSC:   puts("GSC tree weights"); break;
  }
  printf("New CM file:                       %s %s\n",
	 cmfile, do_append? "[appending]" : "");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /* EPN 11.07.05 - EFF_ENTROPY effective sequence number strategy
   *                ported from HMMER 2.4devl. 
   * If we're using the entropy-target strategy for effective
   * sequence number calculation, set the default target entropy loss.
   */
  if (eff_strategy == EFF_ENTROPY) 
    {
      /* 11.23.05
       * EPN Following line enforces an equiprobable 
       * background distribution for entropy weighting
       * (causing etarget to be unaffected by a null 
       * distribution read in with the --null option). 
       */
      etarget = 2.0 - eloss;
      /* 11.23.05 - Below is the old strategy. Admittedly,
       * not sure which (if either) is 'right'. 
       */
      /*etarget = FEntropy(randomseq, MAXABET) - eloss;*/
    }
  /* End of effective seq number port code block. */

  /*********************************************** 
   * Get alignment(s), build CMs one at a time
   ***********************************************/

  nali = 0;
  while ((msa = MSAFileRead(afp)) != NULL)
    {

      if (Alphabet_type == hmmNOTSETYET)
	{
	  printf("%-40s ... ", "Determining alphabet");
	  fflush(stdout);
	  DetermineAlphabet(msa->aseq, msa->nseq);
	  if      (Alphabet_type == hmmNUCLEIC) puts("done. [DNA]");
	  else if (Alphabet_type == hmmAMINO)   puts("done. [protein]");
	  else                                  puts("done.");
	}

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

      eff_nseq_set = FALSE;

      /* --- Post-alphabet initialization section ---
       * If we do this before we've seen the first alignment, then
       * Alphabet_size is uninitialized, and CMReadNullModel() won't
       * work. Not a good reason I know, assuming our Alphabet_size
       * is always 4... 
       * A consequence of stealing code from HMMER.
       */
      if(nali == 0)
	{
	  /* Set up the null/random seq model */
	  if (rndfile == NULL)  CMDefaultNullModel(randomseq);
	  else                  CMReadNullModel(rndfile, randomseq);
	}

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
      if(build_p7hmm)
	hmmer_DigitizeAlignment(msa, &p7dsq);
      printf("done.\n");

      if (eff_strategy == EFF_NONE)
	{
	  eff_nseq = (float) msa->nseq;
	  eff_nseq_set = TRUE;
	}

      /* Construct a model, and collect observed counts.
       * Note on "treeforce": this is the number of sequences that we
       *   will ignore for the purposes of count-collection and parameterization
       *   of the CM. We will only use this for debugging. These seqs
       *   are then dumped as full parsetrees to stdout.
       */
      printf("%-40s ... ", "Constructing model architecture");
      fflush(stdout);
      HandModelmaker(msa, dsq, use_rf, gapthresh, &cm, &mtr);

      printf("done.\n");

      if(do_balance)
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
      if(do_detach)
	{
	  printf("%-40s ... ", "Finding and checking dual inserts");
	  cm_find_and_detach_dual_inserts(cm, 
					  TRUE,   /* Do check (END_E-1) insert states have 0 counts */
					  FALSE); /* Don't detach the states yet, wait til CM is priorified */
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

      /* EPN 11.07.05 - EFF_ENTROPY effective sequence number strategy
       *                ported from HMMER 2.4devl. 
       * Effective sequence number calculation.
       * (if we don't have eff_nseq yet, calculate it now).
       */
      if (! eff_nseq_set) {
	if (eff_strategy == EFF_ENTROPY) {
	  printf("%-40s ... ", "Determining eff seq # by entropy target");
	  fflush(stdout);
	  eff_nseq = CM_Eweight(cm, pri, (float) msa->nseq, etarget);
	}
	/*EPN 11.28.05
	 * Uncomment this block for relative entropy weighting.
	 * else if (eff_strategy == EFF_RELENTROPY) {
	 * printf("%-40s ... ", "Determining eff seq # by relative entropy target");
	 * fflush(stdout);
	 * eff_nseq = CM_Eweight_RE(cm, pri, (float) msa->nseq, eloss, randomseq);
	 * }
	 */
	else Die("no effective seq #: shouldn't happen");

	CMRescale(cm, eff_nseq / (float) msa->nseq);
	eff_nseq_set = TRUE;
	printf("done. [%.2f]\n", eff_nseq);
      }/* End of effective seq number port code block. */
      
      /* Convert to probabilities, and the global log-odds form
       * we save the model in.
       */
      printf("%-40s ... ", "Converting counts to probabilities"); fflush(stdout);
      CMSetNullModel(cm, randomseq);
      if(!no_prior)
	PriorifyCM(cm, pri);

      if(do_detach) /* Detach dual inserts where appropriate, if
		     * we get here we've already checked these states */
	{
	  cm_find_and_detach_dual_inserts(cm, 
					  FALSE, /* Don't check states have 0 counts (they won't due to priors) */
					  TRUE); /* Detach the states by setting trans probs into them as 0.0   */
	}
      CMLogoddsify(cm);
      printf("done.\n");

      /* EPN 08.18.05
       * Calculate W for the model based on the seed seqs' tracebacks.
       *              Brief expt on effect of different bandp and 
       *              safe_windowlen scaling factor values in 
       *              ~nawrocki/notebook/5_0818_inf_cmbuild_bands/00LOG
       */
      printf("%-40s ... ", "Calculating max hit length for model"); fflush(stdout);
      safe_windowlen = 2 * MSAMaxSequenceLength(msa);

      /* EPN 11.13.05 
       * BandCalculationEngine() replaces BandDistribution() and BandBounds().
       * See ~nawrocki/notebook/5_1111_inf_banded_dist_vs_bandcalc/00LOG for details.
       */
      while(!(BandCalculationEngine(cm, safe_windowlen, bandp, save_gamma, &dmin, &dmax, &gamma, do_local)))
	{
	  free(dmin);
	  free(dmax);
	  FreeBandDensities(cm, gamma);	  
	  /*printf("Failure in BandCalculationEngine(). W:%d | bandp: %4e\n", safe_windowlen, bandp);*/
	  safe_windowlen *= 2;
	}

      /*printf("Success in BandCalculationEngine(). W:%d | bandp: %4e\n", safe_windowlen, bandp);*/
      /*debug_print_bands(cm, dmin, dmax);*/
      cm->W = dmax[0];
      printf("done. [%d]\n", cm->W);

      /*11.15.05 EPN Set the EL self transition score, by default its log2(0.94).*/
      cm->el_selfsc = sreLOG2(el_selfprob);
      /* Next line is very hacky. 
       * We want to avoid underflow errors. structs.h explains
       * how IMPOSSIBLE must be > -FLT_MAX/3 so we can add it together 3 
       * times with an underflow. Here, we may potentially be adding el_selfsc
       * together W times. (And W can change in cmsearch or cmalign). Here
       * we'll ensure we can multiply el_selfsc by 2W and still avoid underflows,
       * and we'll check in cmsearch to make sure that W * cm->el_selfsc < (IMPOSSIBLE*3)
       * and we'll die if it isn't. We shouldn't face this in cmsearch situation
       * unless the user wants to set W as something greater than twice what
       * it is set as in the .cm file.
       */
      if(cm->el_selfsc < (IMPOSSIBLE/(2 * cm->W)))
	cm->el_selfsc = (IMPOSSIBLE/(2 * cm->W));
      
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

      /* Optionally, build a plan 7 HMMER 2.4 HMM */
      if(build_p7hmm)
	{
	  /*****************************************************************/
	  /* 10.10.05 Build a plan 7 HMM from the seed alignment.  
	   * Most code torn out of HMMER 2.4devl */
	  printf("%-40s ... ", "Building HMM with hmmer's P7Fastmodelmaker()");
	  checksum = GCGMultchecksum(msa->aseq, msa->nseq);
	  P7Fastmodelmaker(msa, p7dsq, gapthresh, &p7hmm, &p7tr);
	  p7hmm->checksum = checksum;
	  printf("done.\n");
	  
	  /* HMM entropy weighting strategy is same as for the CM
	   * by using the effective sequence number (eff_nseq)
	   * from the CM.
	   */
	  
	  if (eff_strategy == EFF_ENTROPY) {
	    printf("%-40s ... ", "Using eff seq # from CM entropy target");
	    fflush(stdout);
	    Plan7Rescale(p7hmm, eff_nseq / (float) msa->nseq);
	    printf("done. [%.2f]\n", eff_nseq);
	  }
	  
	  /* Record the null model in the HMM;
	   * add prior contributions in pseudocounts and renormalize.
	   */
	  printf("%-40s ... ", "Converting counts to probabilities");
	  fflush(stdout);
	  Plan7SetNullModel(p7hmm, randomseq, p1);
	  P7PriorifyHMM(p7hmm, p7pri);
	  printf("done.\n");

	  /* Model configuration, temporary.
	   * hmmbuild assumes that it's given an alignment of single domains,
	   * and the alignment may contain fragments. So, for the purpose of
	   * scoring the sequences (or, optionally, MD/ME weighting),
	   * configure the model into hmmsw mode. Later we'll
	   * configure the model according to how the user wants to
	   * use it.
	   */
	  Plan7SWConfig(p7hmm, 0.5, 0.5);

	  /* Give the model a name.
	   * We deal with this differently depending on whether
	   * we're in an alignment database or a single alignment.
	   * 
	   * If a single alignment, priority is:
	   *      1. Use -n <name> if set.
	   *      2. Use msa->name (avail in Stockholm or SELEX formats only)
	   *      3. If all else fails, use alignment file name without
	   *         filename extension (e.g. "globins.slx" gets named "globins"
	   *         
	   * If a multiple MSA database (e.g. Stockholm/Pfam), 
	   * only msa->name is applied. -n is not allowed.
	   * if msa->name is unavailable, or -n was used,
	   * a fatal error is thrown.
	   * 
	   * Because we can't tell whether we've got more than one
	   * alignment 'til we're on the second one, these fatal errors
	   * only happen after the first HMM has already been built.
	   * Oh well.
	   */
	  printf("%-40s ... ", "Setting model name, etc.");
	  fflush(stdout);
	  if (nali == 0)		/* first (only?) HMM in file:  */
	    {
	      if      (setname != NULL)   name = Strdup(setname);
	      else if (msa->name != NULL) name = Strdup(msa->name);
	      else                        name = FileTail(alifile, TRUE);
	    }
	  else
	    {
	      if (setname != NULL) 
		Die("Oops. Wait. You can't use -n with an alignment database.");
	      else if (msa->name != NULL) name = Strdup(msa->name);
	      else
		Die("Oops. Wait. I need name annotation on each alignment.\n");
	    }
	  Plan7SetName(p7hmm, name);
	  free(name);
	  /* Transfer other information from the alignment to
	   * the HMM. This typically only works for Stockholm or SELEX format
	   * alignments, so these things are conditional/optional.
	   */
	  if (msa->acc  != NULL) Plan7SetAccession(p7hmm,   msa->acc);
	  if (msa->desc != NULL) Plan7SetDescription(p7hmm, msa->desc);
	  
	  if (msa->cutoff_is_set[MSA_CUTOFF_GA1] && msa->cutoff_is_set[MSA_CUTOFF_GA2])
	    { p7hmm->flags |= PLAN7_GA; p7hmm->ga1 = msa->cutoff[MSA_CUTOFF_GA1]; p7hmm->ga2 = msa->cutoff[MSA_CUTOFF_GA2]; }
	  if (msa->cutoff_is_set[MSA_CUTOFF_TC1] && msa->cutoff_is_set[MSA_CUTOFF_TC2])
	    { p7hmm->flags |= PLAN7_TC; p7hmm->tc1 = msa->cutoff[MSA_CUTOFF_TC1]; p7hmm->tc2 = msa->cutoff[MSA_CUTOFF_TC2]; }
	  if (msa->cutoff_is_set[MSA_CUTOFF_NC1] && msa->cutoff_is_set[MSA_CUTOFF_NC2])
	    { p7hmm->flags |= PLAN7_NC; p7hmm->nc1 = msa->cutoff[MSA_CUTOFF_NC1]; p7hmm->nc2 = msa->cutoff[MSA_CUTOFF_NC2]; }
	  /* Record some other miscellaneous information in the HMM,
	   * like how/when we built it.
	   */
	  Plan7ComlogAppend(p7hmm, argc, argv);
	  Plan7SetCtime(p7hmm);
	  p7hmm->nseq = msa->nseq;
	  printf("done. [%s]\n", p7hmm->name); 
	  /* Print information for the user
	   */
	  printf("\nConstructed a profile HMM (length %d)\n", p7hmm->M);
	  PrintPlan7Stats(stdout, p7hmm, p7dsq, msa->nseq, p7tr); 
	  printf("\n");
	  /* Configure the model for chosen algorithm
	   */
	  
	  printf("%-40s ... ", "Finalizing model configuration");
	  fflush(stdout);
	  switch (cfg_strategy) {
	  case P7_BASE_CONFIG:  Plan7GlobalConfig(p7hmm);              break;
	  case P7_SW_CONFIG:    Plan7SWConfig(p7hmm, swentry, swexit); break;
	  case P7_LS_CONFIG:    Die("whoa, we're not set up for LS HMM config.");
	  case P7_FS_CONFIG:    Die("whoa, we're not set up for FS HMM config.");
	  default:              Die("bogus configuration choice");
	  }
	  printf("done.\n");
	  
	  /* Save new HMM to disk: open a file for appending or writing.
	   */
	  printf("%-40s ... ", "Saving model to file");
	  fflush(stdout);
	  if (do_binary) WriteBinHMM(hmmfp, p7hmm);
	  else           WriteAscHMM(hmmfp, p7hmm);
	  printf("done.\n");

	} /* end of if(build_hmm) */

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
	      fprintf(ofp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr, dsq[idx], FALSE));;
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
	    fprintf(regressfp, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr, dsq[idx], FALSE));
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
	      printf("  SCORE : %.2f bits\n", ParsetreeScore(cm, tr, dsq[idx], FALSE));
	      ParsetreeDump(stdout, tr, cm, dsq[idx]);
	      printf("//\n");
	      FreeParsetree(tr);
	    }
	}

      /* EPN 08.18.05 Detailed band info for the training set (seed seqs).
       */
      if (bandfile != NULL)       
	{
	  printf("%-40s ... ", "Saving band information"); fflush(stdout);
	  if ((ofp = fopen(bandfile,"w")) == NULL)
	    Die("failed to open band file %s", bandfile);

	  /* We want band information for bands we'd use in a cmsearch
	   * with bandp as its set now (default is 1E-7, but it can
	   * be set at the command line). We already have gamma from
	   * the band calculation we used to get cm->W.
	   */
	  fprintf(ofp, "acc:%s\n", msa->acc);
	  fprintf(ofp, "bandp:%g\n", bandp);
	  for (v = 0; v < cm->M; v++)
	    if(cm->sttype[v] == S_st)
		PrintBandDensity(ofp, gamma, v, cm->W, dmin[v], dmax[v]);
	  for (idx = treeforce; idx < msa->nseq; idx++)
	    {
	      fprintf(ofp, "> %s\n", msa->sqname[idx]);
	      tr = Transmogrify(cm, mtr, dsq[idx], msa->aseq[idx], msa->alen);
	      model_trace_info_dump(ofp, cm, tr, msa->aseq[idx]); 
	      FreeParsetree(tr);
	    }
	  fprintf(ofp, "//\n");
	  fclose(ofp);
	  printf("done. [%s]\n", bandfile);
	}

      puts("");
      SummarizeCM(stdout, cm);  
      puts("");
      CYKDemands(cm, avlen);     

      /* Free aln specific HMM related data structures */
      if(build_p7hmm)
	{
	  Free2DArray((void**)p7dsq, msa->nseq);
	  for (idx = 0; idx < msa->nseq; idx++) 
	    P7FreeTrace(p7tr[idx]);
	  free(p7tr);
	  FreePlan7(p7hmm);
	}

      /* Free aln specific CM related data structures */
      FreeParsetree(mtr);
      Free2DArray((void**)dsq, msa->nseq);
      MSAFree(msa);
      fflush(cmfp);
      puts("//\n");
      nali++;

      FreeBandDensities(cm, gamma);	  
      free(dmin);
      free(dmax);
      FreeCM(cm);
    }


  /* Clean up and exit
   */
  if (regressionfile != NULL) fclose(regressfp);
  StopwatchFree(watch);
  MSAFileClose(afp);
  Prior_Destroy(pri);
  fclose(cmfp);
  SqdClean();

  P7FreePrior(p7pri);
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

/* EPN 08.18.05
 * model_trace_info_dump()
 * Function: model_trace_info_dump
 *
 * Purpose:  Given a trace from a sequence used to create the model, 
 *           print the subsequence length rooted at each start state.  
 *           Tricky because the sequence positions in a Parsetree_t tr
 *           returned from Transmogrify refer to aligned positions.
 *           We want subsequence lengths that refer to unaligned lengths.
 * 
 * Args:    ofp      - filehandle to print to
 *          cm       - the CM
 *          tr       - the parsetree (trace)
 *          aseq     - the aligned sequence the trace corresponds to
 * Returns: (void) 
 */

static void
model_trace_info_dump(FILE *ofp, CM_t *cm, Parsetree_t *tr, char *aseq)
{
  int a, i, j, tpos, d, l, r;
  int *map;

  map = MallocOrDie (sizeof(int) * strlen(aseq));
  
  a=0;
  for (i = 0; i < strlen(aseq); i++)
    if (! isgap(aseq[i])) map[i] = a++;
    else map[i] = -1;

  for (tpos = 0; tpos < tr->n; tpos++)
    if(cm->sttype[tr->state[tpos]] == S_st)
      {
	l = tr->emitl[tpos]-1;
	r = tr->emitr[tpos]-1;
	i = map[l];
	j = map[r];
	/* tr->emitl[tpos]-1 might map to a gap (root node emits the gaps
	 * also). So we look for first residue that exists in the unaligned
	 * seq.  Then we do the same for j, looking backwards.
	 */ 
	while (i == -1)
	  i = map[++l];
	while (j == -1)
	  j = map[--r];

	d = j-i+1;
	/* assume ofp is open (probably not good) */
	fprintf(ofp, "state:%d d:%d\n", tr->state[tpos], d);
	/*fprintf(ofp, "state:%d d:%d i:%d j:%d emitl:%d emitr:%d\n", tr->state[tpos], d, i, j, tr->emitl[tpos], tr->emitr[tpos]);*/
      }
  free(map);
}

/* EPN 08.18.05
 * PrintBandDensity()
 * Function: PrintBandDensity
 *
 * Purpose:  Given gamma, a state index v, and a W (maximum hit len)
 *           print out the probability that a subsequence rooted at 
 *           v will have lengths 0 to W.
 *
 * Args:    fp       - filehandle to print to
 *          gamma    - cumulative probability distribution P(length <= n) for state v;
 *                     [0..v..M-1][0..W] 
 *          v        - state index         
 *          W        - maximum subseq len in DP
 *          min      - dmin[v]
 *          max      - dmax[v]
 *
 * Returns: (void) 
 */

static void
PrintBandDensity(FILE *fp, double **gamma, int v, int W, int min, int max)
{
  int n;

  fprintf(fp, "band for state:%d min:%d max:%d\n", v, min, max);
  for (n = 0; n <= W; n++)
    fprintf(fp, "%d:%.12f\n", n, gamma[v][n]);
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

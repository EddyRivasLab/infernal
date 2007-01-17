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
#include "cm_wrappers.h"	/* alignment functions                   */

void SummarizeAlignOptions(CM_t *cm);

static char banner[] = "cmscore - score RNA covariance model against sequences";

static char usage[]  = "\
Usage: cmscore [-options] <cmfile> <sequence file>\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -i     : print individual timings & score comparisons, not just summary\n\
";

static char experts[] = "\
  * Modes:\n\
   --dc          : compare divide and conquer versus standard CYK (default)\n\
   --qdb         : compare QDB alignment versus non-banded d&c/standard\n\
   --hbanded     : compare HMM banded alignment versus non-banded d&c/standard\n\
   --hbandp <f>  : tail loss prob for --hbanded [default: 0.0001]\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
   --hmmonly     : align with the CM Plan 9 HMM\n\
   --regress <f> : save regression test data to file <f>\n\
   --scoreonly   : for full CYK/inside stage, do only score, save memory\n\
   --smallonly   : do only d&c, don't do full CYK/inside\n\
   --stringent   : require the two parse trees to be identical\n\
   --dcqdb       : compare QDB d&c versus QDB standard CYK\n\
   --beta <f>    : tail loss prob for QDB [default:0.0000001]\n\
   --trees       : print parsetrees\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-i", TRUE, sqdARG_NONE }, 
  { "--dc",         FALSE, sqdARG_NONE },
  { "--qdb",        FALSE, sqdARG_NONE },
  { "--hbanded",    FALSE, sqdARG_NONE },
  { "--hbandp",     FALSE, sqdARG_FLOAT},
  { "--sums",       FALSE, sqdARG_NONE},
  { "--hmmonly",    FALSE, sqdARG_NONE},
  { "--informat",   FALSE, sqdARG_STRING },
  { "--local",      FALSE, sqdARG_NONE },
  { "--regress",    FALSE, sqdARG_STRING },
  { "--scoreonly",  FALSE, sqdARG_NONE },
  { "--smallonly",  FALSE, sqdARG_NONE },
  { "--stringent",  FALSE, sqdARG_NONE },
  { "--dcqdb",      FALSE, sqdARG_NONE },
  { "--beta",       FALSE, sqdARG_FLOAT},
  { "--trees",       FALSE, sqdARG_FLOAT}
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  char            *seqfile;     /* file to read sequences from */
  ESL_SQFILE	  *seqfp;       /* open seqfile for reading                 */
  ESL_SQ	 **sq;          /* the target sequences to align            */
  int              format;      /* format of sequence file */
  int              status;      /* status of seqfp                          */
  CMFILE          *cmfp;        /* open CM file for reading */
  CM_t            *cm;          /* a covariance model       */
  char            *seq;         /* RNA sequence */
  SQINFO           sqinfo;      /* optional info attached to seq */
  char            *dsq;         /* digitized RNA sequence */
  Stopwatch_t     *watch1;      /* for timing of stage 1 alignment */
  Stopwatch_t     *watch2;      /* for timing of stage 2 alignment */
  float            sc1,  sc2;	/* score of a sequence */
  Parsetree_t     *tr1, *tr2;	/* a traceback */
  
  int   do_local;		/* TRUE to align locally w.r.t. model       */
  int   do_scoreonly;		/* TRUE for score-only (small mem) full CYK */
  int   do_smallonly;		/* TRUE to do only d&c, not full CYK/inside */
  int   do_trees;		/* TRUE to print parse trees to stdout      */
  int   do_individuals;         /* TRUE to print individual scores/times    */
  int   compare_stringently;	/* TRUE to demand identical parse trees     */
  char *regressfile;		/* name of regression data file to save     */
  FILE *regressfp;              /* open filehandle for writing regressions  */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */
  double   qdb_beta;	        /* tail loss probability for query dependent banding */
  double   hbandp; 	        /* tail loss probability for HMM banding */
  int      use_sums;            /* TRUE: use the posterior sums w/HMM bands */
  int     *dmin;                /* minimum d bound for state v, [0..v..M-1] */
  int     *dmax;                /* maximum d bound for state v, [0..v..M-1] */
  int      be_quiet;            /* TRUE to not print out scores during alignment */
  int      i;                   /* counter over seqs                       */
  int      bdump_level;         /* verbosity level for band-related printfs */
  int      debug_level;         /* verbosity level for debugging printfs */
  float    spdup;               /* stage 1 time elapsed / stage 2 time elapsed */

  ESL_SQ	 **s1_sq;          /* the target sequences to align            */
  Parsetree_t    **s1_tr;          /* stage 1 parse trees for the sequences */
  int              s1_nseq;        /* number of target sequences             */
  float           *s1_sc;          /* scores from stage 1 */

  ESL_SQ	 **s2_sq;          /* the target sequences to align            */
  Parsetree_t    **s2_tr;          /* stage 2 parse trees for the sequences */
  int              s2_nseq;        /* number of target sequences             */
  float           *s2_sc;          /* scores from stage 2 */

  int              diff_ct;        /* number of seqs w/ different scores b/t stage 1 and 2 */ 
  float            diff_sc;        /* total score difference b/t stage 1 and 2 */ 
  /* Stage 1 options: alignment options for first alignment */
  int      s1_do_qdb;           /* TRUE to do qdb CYK in stage 1 */
  int      s1_do_small;         /* TRUE to do d&c in stage 1 */

  /* Stage 2 options: alignment options for second alignment */
  int      do_s2;               /* TRUE to do stage 2 alignment at all */
  int      s2_do_qdb;           /* TRUE to do qdb CYK in stage 2 */
  int      s2_do_small;         /* TRUE to do qdb CYK in stage 2 */
  int      s2_do_hbanded;       /* TRUE to do HMM banded CYK in stage 2 */
  int      s2_do_hmmonly;       /* TRUE: stage 2 align with the HMM, not the CM     */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format              = SQFILE_UNKNOWN;
  do_local            = FALSE;
  do_scoreonly        = FALSE;
  do_smallonly        = FALSE;
  do_trees            = FALSE;
  do_individuals      = FALSE;
  regressfile         = NULL;
  qdb_beta            = 0.0000001;
  compare_stringently = FALSE;
  s1_do_qdb           = FALSE;
  s1_do_small         = TRUE;
  do_s2               = TRUE;
  s2_do_qdb           = FALSE;
  s2_do_small         = FALSE;
  s2_do_hbanded       = FALSE;
  s2_do_hmmonly       = FALSE;
  qdb_beta            = DEFAULT_BETA;
  hbandp              = DEFAULT_HBANDP;
  use_sums            = FALSE;
  be_quiet            = FALSE;
  bdump_level         = 0;
  debug_level         = 0;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-i")          == 0) do_individuals      = TRUE;
    else if (strcmp(optname, "--local")     == 0) do_local            = TRUE;
    else if (strcmp(optname, "--regress")   == 0) regressfile         = optarg;
    else if (strcmp(optname, "--regress")   == 0) regressfile         = optarg;
    else if (strcmp(optname, "--smallonly") == 0) do_smallonly        = TRUE;
    else if (strcmp(optname, "--scoreonly") == 0) do_scoreonly        = TRUE;
    else if (strcmp(optname, "--stringent") == 0) compare_stringently = TRUE;
    else if (strcmp(optname, "--beta")      == 0) qdb_beta            = atof(optarg);
    else if (strcmp(optname, "--dc")        == 0) { /* this is default */ }
    else if (strcmp(optname, "--hbanded")   == 0) s2_do_hbanded       = TRUE;
    else if (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--hmmonly")   == 0) { Die("--hmmonly not yet implemented.\n"); s2_do_hmmonly   = TRUE; }
    else if (strcmp(optname, "--qdb")       == 0) 
      {
	if(s2_do_qdb) Die("--qdb and --dcqdb combination doesn't make sense, pick one.\n");
	s2_do_qdb = TRUE;
	s2_do_small = TRUE;
      }
    else if (strcmp(optname, "--dcqdb")     == 0) 
      {
	if(s2_do_qdb) Die("--qdb and --dcqdb combination doesn't make sense, pick one.\n");
	s1_do_qdb = s2_do_qdb = TRUE;
      }
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
  
  /*****************************************************************
   * Input the CM for stage 1 alignment (we configure it with ConfigCM() later).
   *****************************************************************/

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);

  cm->beta   = qdb_beta; /* this will be DEFAULT_BETA unless changed at command line */
  cm->hbandp = hbandp;   /* this will be DEFAULT_HBANDP unless changed at command line */
  /*****************************************************************
   * Open the target sequence file
   *****************************************************************/
  status = esl_sqfile_Open(seqfile, format, NULL, &seqfp);
  if (status == eslENOTFOUND) esl_fatal("No such file."); 
  else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
  else if (status == eslEINVAL) esl_fatal("Can’t autodetect stdin or .gz."); 
  else if (status != eslOK) esl_fatal("Failed to open sequence database file, code %d.", status); 

				/* open regression test data file */
  if (regressfile != NULL) {
    if ((regressfp = fopen(regressfile, "w")) == NULL)
      Die("Failed to open regression test file %s", regressfile);
  }

  /*********************************************** 
   * Show the banner
   ***********************************************/

  if (! be_quiet) 
    {
      MainBanner(stdout, banner);
      printf(   "CM file:              %s\n", cmfile);
      printf(   "Sequence file:        %s\n", seqfile);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /****************************************************************************
   * Do the work: read in seqs and align them twice. 
   ****************************************************************************/

  /* Configure the CM for stage 1 alignment */
  /* Update cm->config_opts and cm->align_opts based on command line options */
  /* enable option to check parsetree score against the alignment score */
  watch1 = StopwatchCreate(); 
  StopwatchZero(watch1);
  StopwatchStart(watch1);

  cm->align_opts |= CM_ALIGN_CHECKPARSESC;
  if (do_individuals) cm->align_opts |= CM_ALIGN_TIME;
  if (do_trees) cm->align_opts |= CM_ALIGN_PRINTTREES;
  if (do_local) cm->config_opts |= CM_CONFIG_LOCAL;
  if(s1_do_qdb)          
    { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }
  if(!s1_do_small) cm->align_opts |= CM_ALIGN_NOSMALL;

  ConfigCM(cm, NULL, NULL);

  printf("Stage 1 alignment:\n");
  SummarizeAlignOptions(cm);
  serial_align_targets(seqfp, cm, &s1_sq, &s1_tr, NULL, &s1_nseq, bdump_level, debug_level, 
		       (!do_individuals));
  StopwatchStop(watch1);
  esl_sqfile_Close(seqfp);
  FreeCM(cm);

  if(do_s2)
    {
      /*****************************************************************
       * Input the CM for stage 2 alignment (we configure it with ConfigCM() later).
       *****************************************************************/

      if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
	Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
      if (! CMFileRead(cmfp, &cm))
	Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
      if (cm == NULL) 
	Die("%s empty?\n", cmfile);
      CMFileClose(cmfp);
      
      watch2 = StopwatchCreate(); 
      StopwatchZero(watch2);
      StopwatchStart(watch2);

      cm->beta   = qdb_beta; /* this will be DEFAULT_BETA unless changed at command line */
      cm->hbandp = hbandp;   /* this will be DEFAULT_HBANDP unless changed at command line */
      /* Reopen the sequence file, we'll wastefully reread the seqs. */
      status = esl_sqfile_Open(seqfile, format, NULL, &seqfp);
      if (status == eslENOTFOUND) esl_fatal("No such file."); 
      else if (status == eslEFORMAT) esl_fatal("Format unrecognized."); 
      else if (status == eslEINVAL) esl_fatal("Can’t autodetect stdin or .gz."); 
      else if (status != eslOK) esl_fatal("Failed to open sequence database file, code %d.", status); 

      /* 
      /* Configure the CM for stage 2 alignment */
      /* enable option to check parsetree score against the alignment score */
      cm->align_opts |= CM_ALIGN_CHECKPARSESC;
      if (do_individuals) cm->align_opts |= CM_ALIGN_TIME;
      if (do_trees) cm->align_opts  |= CM_ALIGN_PRINTTREES;
      if (do_local) cm->config_opts |= CM_CONFIG_LOCAL;
      if(s2_do_qdb)          
	{ 
	  cm->align_opts  |= CM_ALIGN_QDB;
	  cm->config_opts |= CM_CONFIG_QDB;
	}
      if(!s2_do_small)  cm->align_opts |= CM_ALIGN_NOSMALL;
      if(s2_do_hbanded) cm->align_opts |= CM_ALIGN_HBANDED;
      if(use_sums)      cm->align_opts |= CM_ALIGN_SUMS;
      if(s2_do_hmmonly)    cm->align_opts |= CM_ALIGN_HMMONLY;

      ConfigCM(cm, NULL, NULL);
      printf("Stage 2 alignment:\n");
      SummarizeAlignOptions(cm);
      serial_align_targets(seqfp, cm, &s2_sq, &s2_tr, NULL, &s2_nseq, bdump_level, debug_level, 
			   (!do_individuals));
      StopwatchStop(watch2);
      esl_sqfile_Close(seqfp);
    }

  if(s1_nseq != s2_nseq) Die("ERROR, in stage 1, aligned %d seqs, in stage 2, aligned %d seqs (should be equal).\n", s1_nseq, s2_nseq);

  /* Compare parsetrees from stage 1 and stage 2 and collect stats */
  s1_sc = MallocOrDie(sizeof (float) * s1_nseq);
  s2_sc = MallocOrDie(sizeof (float) * s1_nseq);

  diff_sc = 0.;
  diff_ct = 0;
  for(i = 0; i < s1_nseq; i++)
    {
      s1_sc[i] = ParsetreeScore(cm, s1_tr[i], s1_sq[i]->dsq, FALSE);
      s2_sc[i] = ParsetreeScore(cm, s2_tr[i], s2_sq[i]->dsq, FALSE);
      if(do_individuals)
	printf("%-12s S1: %.3f S2: %.3f diff: %.3f\n", s1_sq[i]->name, s1_sc[i], s2_sc[i], fabs(s1_sc[i]-s2_sc[i]));
      if(fabs(s1_sc[i] -  s2_sc[i]) > 0.0001)
	{
	  diff_ct++;
	  diff_sc += fabs(s1_sc[i] - s2_sc[i]);
	}
    }

  /* Print summary */ 

  printf("Results summary:\n");
  printf("---------------------------------\n");
  printf("Number seqs aligned:     %d\n", s1_nseq);
  StopwatchDisplay(stdout, "Stage 1 time:            ", watch1);
  if(do_s2)
    {
      StopwatchDisplay(stdout, "Stage 2 time:            ", watch2);
      spdup = watch1->user / watch2->user;
      printf("2/1 speedup (user):      %.2f\n", spdup);
    }
  printf("Total bit score diff:      %.2f\n", (diff_sc / ((float) s1_nseq)));
  printf("Avg bit score diff:      %.2f\n", (diff_sc / ((float) s1_nseq)));
  printf("Num   diff (>1e-4):      %d\n", (diff_ct));
  printf("Fract diff (>1e-4):      %.5f\n", (((float) diff_ct) / ((float) s1_nseq)));
  printf("\n");

  StopwatchFree(watch1);
  if(do_s2)
    {
      StopwatchFree(watch2);
    }

  /* 
  printf("\n\n(%4d) %4d (%.3f) %.3f (%.3f)\n", s1_nseq, diff_ct, ((float) diff_ct / s1_nseq), diff_sc, ((float) diff_sc / diff_ct));

  FreeCM(cm);
  /* Save regression test data
   */
  /*  if (regressfile != NULL) {
    if (tr1 != NULL) ParsetreeDump(regressfp, tr1, cm, dsq);
    if (tr2 != NULL) ParsetreeDump(regressfp, tr2, cm, dsq);
    }*/
  /*
      FreeSequence(seq, &sqinfo);
      if (tr1 != NULL) FreeParsetree(tr1);  
      if (tr2 != NULL) FreeParsetree(tr2); 
      free(dsq);
      }*/

  /*  if (regressfile != NULL) fclose(regressfp);

  StopwatchFree(watch);
  SqdClean();
  return EXIT_SUCCESS;*/
}

/* Function: SummarizeAlignOptions
 * Date:     EPN, Wed Jan 17 09:08:18 2007
 * Purpose:  Print out alignment options in pretty format. 
 */
void SummarizeAlignOptions(CM_t *cm)
{
  printf("---------------------------------\n");
  /* Algorithm */
  if(cm->align_opts & CM_ALIGN_INSIDE)
    printf("Algorithm:               Inside\n");
  else if(cm->align_opts & CM_ALIGN_NOSMALL)
    printf("Algorithm:               CYK Standard\n");
  else if(cm->align_opts & CM_ALIGN_HMMONLY)
    printf("Algorithm:               CP9 HMM Viterbi\n");
  else 
    printf("Algorithm:               CYK D&C\n");

  /* Bands */
  if(cm->align_opts & CM_ALIGN_HBANDED)
    {
      if(cm->align_opts & CM_ALIGN_SUMS)
	printf("Bands:                   CP9 HMM (sums)\n");
      else
	printf("Bands:                   CP9 HMM\n"); 
      printf("Tail loss:               %g\n", cm->hbandp);
    }
  else if(cm->align_opts & CM_ALIGN_QDB)
    {
      printf("Bands:                   QDB\n"); 
      printf("Tail loss:               %g\n", cm->beta);
    }
  else
    {
      printf("Bands:                   None\n"); 
      printf("Tail loss:               0.0\n", cm->beta);
    }
  /* Locality */
  if(cm->config_opts & CM_CONFIG_LOCAL)
    printf("Local:                   Yes\n");
  else
    printf("Local:                   No\n");

  /* Sub mode? */
  if(cm->align_opts & CM_ALIGN_SUB)
    printf("Sub mode.\n");
  if(cm->align_opts & CM_ALIGN_FSUB)
    printf("Full Sub mode.\n");

  printf("---------------------------------\n");
  return;
}

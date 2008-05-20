/************************************************************
 * @LICENSE@
 ************************************************************/

/* cmstat.c
 * EPN, Tue Aug 21 12:50:34 2007
 *
 * Display summary statistics for an CM or CM database 
 * (such as Rfam). 
 *
 * Based on SRE's hmmstat.c from HMMER3.
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sqio.h"
#include "esl_stats.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

#define ONELINEOPTS  "-m,--all,--le,--ge,--lfc,--gfc,--lfi,--gfi" /* exclusive choice of summary stats */
#define NOTMOPTS     "--le,--ge,--lfc,--gfc,--lfi,--gfi"    /* incompatible with -g, --beta, and --search */
#define CMCUTOPTS    "-E,-T,--ga,--tc,--nc,--range"         /* exclusive choice for CM cutoff */

static ESL_OPTIONS options[] = {
  /* name           type      default      env  range     toggles      reqs     incomp    help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-g",        eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL,    NOTMOPTS, "configure CM for glocal alignment [default: local]", 1 },
  { "-m",        eslARG_NONE,"default",    NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of model statistics", 1 },
  { "-Z",        eslARG_REAL,   FALSE,     NULL, NULL,      NULL,      NULL,        NULL, "set Z (database size in *Mb*) to <x> for E-value calculations", 1},
  { "--all",     eslARG_NONE,   FALSE,     NULL, NULL,      "-m",      NULL, ONELINEOPTS, "print model, E-value and filter thresholds stats", 1 },
  { "--le",      eslARG_NONE,   FALSE,     NULL, NULL,      "-m",      NULL, ONELINEOPTS, "only print one line summary of  local E-value statistics", 1 },
  { "--ge",      eslARG_NONE,   FALSE,     NULL, NULL,      "-m",      NULL, ONELINEOPTS, "only print one line summary of glocal E-value statistics", 1 },
  { "--beta",    eslARG_REAL,   "1E-7",    NULL, "0<x<1",   NULL,      NULL,    NOTMOPTS, "set tail loss prob for QDB stats to <x>", 1 },
  { "--search",  eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL,    NOTMOPTS, "do search timing experiments", 1 },
  { "--cmL",     eslARG_INT,    "1000",    NULL, "n>0",     NULL,"--search",        NULL, "length of sequences for CM search stats", 1 },
  { "--hmmL",    eslARG_INT,    "100000",  NULL, "n>0",     NULL,"--search",        NULL, "length of sequences for CP9 HMM search stats", 1 },
  { "--qdbfile", eslARG_OUTFILE, NULL,     NULL, NULL,      NULL,      "-m",        NULL, "save query-dependent bands (QDBs) for each state to file <f>", 1 },
  { "--lfi",     eslARG_NONE,   FALSE,     NULL, NULL,      "-m",      NULL, ONELINEOPTS, "only print summary of  local Inside filter threshold stats", 2 },
  { "--gfi",     eslARG_NONE,   FALSE,     NULL, NULL,      "-m",      NULL, ONELINEOPTS, "only print summary of glocal Inside filter threshold stats", 2 },
  { "--lfc",     eslARG_NONE,   FALSE,     NULL, NULL,      "-m",      NULL, ONELINEOPTS, "only print summary of  local CYK    filter threshold stats", 2 },
  { "--gfc",     eslARG_NONE,   FALSE,     NULL, NULL,      "-m",      NULL, ONELINEOPTS, "only print summary of glocal CYK    filter threshold stats", 2 },
  { "-E",        eslARG_REAL,   "0.1",     NULL, "x>0",     NULL,      NULL,   CMCUTOPTS, "print HMM filter stats for cmsearch E cutoff <x>", 2}, 
  { "-T",        eslARG_REAL,   NULL,      NULL, NULL,      NULL,      NULL,   CMCUTOPTS, "print HMM filter stats for cmsearch bit cutoff <x>", 2}, 
  { "--ga",      eslARG_NONE,   NULL,      NULL, NULL,      NULL,      NULL,   CMCUTOPTS, "print HMM filter stats for Rfam GA cutoff", 2}, 
  { "--tc",      eslARG_NONE,   NULL,      NULL, NULL,      NULL,      NULL,   CMCUTOPTS, "print HMM filter stats for Rfam TC cutoff", 2}, 
  { "--nc",      eslARG_NONE,   NULL,      NULL, NULL,      NULL,      NULL,   CMCUTOPTS, "print HMM filter stats for Rfam NC cutoff", 2}, 
  { "--range",   eslARG_NONE,   NULL,      NULL, NULL,      NULL,      NULL,   CMCUTOPTS, "print full range CM/HMM filter cutoff points", 2}, 
  { "--seqfile", eslARG_INFILE, FALSE,     NULL, NULL,      NULL,      NULL,        "-Z", "compute E-value cutoffs for sequence file <f>", 2 },
  { "--toponly", eslARG_NONE,   FALSE,     NULL, NULL,      NULL,"--seqfile",       NULL, "with --seqfile, only consider top-strand", 2 },
  { "--efile",   eslARG_OUTFILE,NULL,      NULL, NULL,      NULL,      NULL,        NULL, "output HMM filter E-val cutoff vs CM E-val cutoff plots to <f>", 3},
  { "--bfile",   eslARG_OUTFILE,NULL,      NULL, NULL,      NULL,      NULL,        NULL, "output HMM filter bit sc cutoff vs CM bit sc cutoff plots to <f>", 3},
  { "--sfile",   eslARG_OUTFILE,NULL,      NULL, NULL,      NULL,      NULL,        NULL, "output predicted survival fraction vs CM cutoff plots to <f>", 3},
  { "--xfile",   eslARG_OUTFILE,NULL,      NULL, NULL,      NULL,      NULL,        NULL, "output predicted xhmm (calcs * HMM) vs CM cutoff plots to <f>", 3},
  { "--afile",   eslARG_OUTFILE,NULL,      NULL, NULL,      NULL,      NULL,        NULL, "output predicted acceleration vs CM cutoff plots to <f>", 3},
  { "--bits",    eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL,        NULL, "with --{s,x,a}file, plot CM bit score cutoffs not E-values", 3},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "display summary statistics for CMs";

static int    summarize_search(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w, FILE *ofp); 
static int    initialize_cm   (CM_t *cm, int cm_mode, int hmm_mode);
static int    print_run_info(const ESL_GETOPTS *go, char *errbuf, ESL_RANDOMNESS *r);
extern int    get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing   */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory error, stopwatch not created.\n");
  esl_stopwatch_Start(w);

  ESL_ALPHABET    *abc = NULL;  /* alphabet                  */
  ESL_RANDOMNESS  *r   = NULL;  /* source of randomness      */
  ESL_STOPWATCH   *s_w = NULL;  /* for timing expts          */
  char            *cmfile;	/* name of input CM file     */ 
  CMFILE          *cmfp;	/* open input CM file stream */
  CM_t            *cm;          /* CM most recently read     */
  int              ncm;         /* CM index                  */
  char             errbuf[cmERRBUFSIZE]; /* for error messages */
  int              status;      /* easel status */
  int              p, i;        /* counters */
  int              fthr_mode;   /* filter threshold mode */
  int              cm_mode, hmm_mode; /* exp tail modes for CM, CP9 */
  long             dbsize;      /* database size E-values correspond to */
  FILE            *qdbfp = NULL;/* for --qdbfile output */
  FILE            *efp = NULL;  /* for --efile output */
  FILE            *bfp = NULL;  /* for --bfile output */
  FILE            *sfp = NULL;  /* for --sfile output */
  FILE            *xfp = NULL;  /* for --xfile output */
  FILE            *afp = NULL;  /* for --afile output */
  int              seen_exps_yet = FALSE; /* set to true if exp tails read */
  int              seen_fthr_yet    = FALSE; /* set to true if filter threshold stats read */
  /* variables for filter threshold stats */
  int              do_filter_stats;                 /* TRUE if --lfc, --gfc, --lfi or --gfi */
  int              do_avg_stats;                    /* TRUE to print average stats */
  float            avg_clen = 0.;                   /* average consensus length */
  float            cm_E, avg_cm_E = 0.;             /* cm E value cutoff */
  float            cm_bit_sc, avg_cm_bit_sc = 0.;   /* cm bit score cutoff */
  float            hmm_E, avg_hmm_E = 0.;           /* hmm E value cutoff */
  float            hmm_bit_sc, avg_hmm_bit_sc = 0.; /* hmm bit score cutoff */
  float            S, avg_S = 0.;                   /* survival fraction */
  float            xhmm, avg_xhmm = 0.;             /* filter scan takes <xhmm> times as long as hmm only scan */
  float            spdup, avg_spdup = 0.;           /* predicted speedup using filter */
  float            tot_xhmm = 0.;                   /* total xhmm for all CMs */
  float            tot_spdup = 0.;                  /* total predicted speedup using filters for all CMs */
  float            cm_ncalcs, tot_cm_ncalcs = 0.;   /* total millions of dp calcs for CM per residue for all CMs */
  float            hmm_ncalcs, tot_hmm_ncalcs = 0.; /* total millions of dp calcs for HMM per residue for all CMs */
  float            tot_cm_surv_plus_fil_calcs = 0.; /* total millions of dp calcs for HMM per residue plus CM search of survivors for all CMs */
  int              do_filter;                       /* TRUE if it's worth filtering for current model */
  int              nfilter = 0;                     /* number of models it's worth filtering for */

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  /* Process command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\n  options for printing filter threshold statistics:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      puts("\n  optional xmgrace plots for --lfc, --gfc, --lfi or --gfi:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\n  where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  do_filter_stats = ((esl_opt_GetBoolean(go, "--all")) || (esl_opt_GetBoolean(go, "--lfc")) || (esl_opt_GetBoolean(go, "--lfi")) || (esl_opt_GetBoolean(go, "--gfc")) || (esl_opt_GetBoolean(go, "--gfi")));

  /* check for incompatible options combos too convoluted for esl_getopts */
  if (esl_opt_GetBoolean(go, "--bits")) { 
    if((esl_opt_IsDefault(go, "--sfile")) && (esl_opt_IsDefault(go, "--xfile")) && (esl_opt_IsDefault(go, "--afile"))) {
      cm_Fail("--bits only works with --sfile, --xfile and/or --afile");
    }
  }
  if((! do_filter_stats) && (! ((esl_opt_IsDefault(go, "--sfile")) && (esl_opt_IsDefault(go, "--xfile")) && (esl_opt_IsDefault(go, "--afile")) && (esl_opt_IsDefault(go, "--efile")) && (esl_opt_IsDefault(go, "--bfile"))))) { 
    cm_Fail("--{e,b,s,x,a}file options only work with --lfc, --lfi, --gfc or --gfi");
  }
  if((! do_filter_stats) && (!(esl_opt_IsDefault(go, "--seqfile")))) {
    cm_Fail("--seqfile only makes sense with one of: --lfc, --lfi, --gfc, or --gfi");
  }
  if((! do_filter_stats) && (! ((esl_opt_IsDefault(go, "-E")) && (esl_opt_IsDefault(go, "-T")) && (esl_opt_IsDefault(go, "--ga")) && (esl_opt_IsDefault(go, "--nc")) && (esl_opt_IsDefault(go, "--tc")) && (esl_opt_IsDefault(go, "--range"))))) {
    cm_Fail("-E,-T,--ga,--nc, --tc, and --range options only work with --lfc, --lfi, --gfc or --gfi");
  }

  cm_banner(stdout, argv[0], banner);
  cmfile     = esl_opt_GetArg(go, 1); 
  r = esl_randomness_CreateTimeseeded();
  s_w = esl_stopwatch_Create();
  if(r   == NULL) cm_Fail("Failed to create RNG, probably out of memory.\n");
  if(s_w == NULL) cm_Fail("Failed to create stopwatch, probably out of memory.\n");

  if((status  = print_run_info (go, errbuf, r))  != eslOK) cm_Fail(errbuf);

  /* Initializations: open the CM file
   */
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    cm_Fail("Failed to open covariance model save file %s\n", cmfile);

  /* open optional output files */
  /* --qdbfile */
  if (! esl_opt_IsDefault(go, "--qdbfile")) 
    if ((qdbfp = fopen(esl_opt_GetString(go, "--qdbfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --qdbfile output file %s\n", esl_opt_GetString(go, "--qdbfile"));
  /* --efile */
  if (! esl_opt_IsDefault(go, "--efile"))
    if ((efp = fopen(esl_opt_GetString(go, "--efile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --efile output file %s\n", esl_opt_GetString(go, "--efile"));
  /* --bfile */
  if (! esl_opt_IsDefault(go, "--bfile"))
    if ((bfp = fopen(esl_opt_GetString(go, "--bfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --bfile output file %s\n", esl_opt_GetString(go, "--bfile"));
  /* --sfile */
  if (! esl_opt_IsDefault(go, "--sfile"))
    if ((sfp = fopen(esl_opt_GetString(go, "--sfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --sfile output file %s\n", esl_opt_GetString(go, "--sfile"));
  /* --xfile */
  if (! esl_opt_IsDefault(go, "--xfile"))
    if ((xfp = fopen(esl_opt_GetString(go, "--xfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --xfile output file %s\n", esl_opt_GetString(go, "--xfile"));
  /* --afile */
  if (! esl_opt_IsDefault(go, "--afile"))
    if ((afp = fopen(esl_opt_GetString(go, "--afile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --afile output file %s\n", esl_opt_GetString(go, "--afile"));

  /* if --seqfile enabled, read the sequence file to get the database size, else we use 1 Mb as db size */
  if(!(esl_opt_IsDefault(go, "--seqfile"))) { 
    /* open input sequence file */
    ESL_SQFILE      *sqfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--seqfile"), eslSQFILE_UNKNOWN, NULL, &sqfp);
    if (status == eslENOTFOUND)    cm_Fail("File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--seqfile"));
    else if (status == eslEFORMAT) cm_Fail("Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--seqfile"));
    else if (status == eslEINVAL)  cm_Fail("Can’t autodetect stdin or .gz."); 
    else if (status != eslOK)      cm_Fail("Sequence file open failed with error %d\n", status);
    /* GetDBSize() reads all sequences, rewinds seq file and returns db size */
    if((status = GetDBSize(sqfp, errbuf, &(dbsize))) != eslOK) cm_Fail(errbuf);
    esl_sqfile_Close(sqfp); 
    if (! esl_opt_GetBoolean(go, "--toponly")) dbsize *= 2;
  }
  else if (!(esl_opt_IsDefault(go, "-Z"))) { dbsize = (long) (esl_opt_GetReal(go, "-Z") * 1000000.); /* convert to Mb then to a long */ }
  else dbsize = FTHR_DBSIZE; /* 1 Mb */

  /* Main body: read CMs one at a time, print stats 
   * Options for only printing 1 line per CM: 
   * if -m (default): print general model stats, and optionally determine search stats (if --search)
   * else if --le:    print  local exp tail stats
   * else if --ge:    print glocal exp tail stats
   * else if --lfc:   print  local CYK    filter threshold stats
   * else if --gfc:   print glocal CYK    filter threshold stats
   * else if --lfc:   print  local Inside filter threshold stats
   * else if --gfc:   print glocal Inside filter threshold stats
   * else if --all:   print all stat categories, one category at a time 
   */
  int do_model    = esl_opt_GetBoolean(go, "-m");
  int do_locale   = esl_opt_GetBoolean(go, "--le");
  int do_glocale  = esl_opt_GetBoolean(go, "--ge");
  int do_localfc  = esl_opt_GetBoolean(go, "--lfc");
  int do_glocalfc = esl_opt_GetBoolean(go, "--gfc");
  int do_localfi  = esl_opt_GetBoolean(go, "--lfi");
  int do_glocalfi = esl_opt_GetBoolean(go, "--gfi");
  /* assert we only have one of our exclusive modes on, getops should've already handled this actually */
  if((do_model + do_locale + do_glocale + do_localfc + do_glocalfc + do_localfi + do_glocalfi) > 1) 
    cm_Fail("error parsing options, exactly 1 or 0 of the following should be true (1):\ndo_model: %d\ndo_locale: %d\ndo_glocale: %d\ndo_localfc: %d\ndo_glocalfc: %d\ndo_localfi: %d\ndo_glocalfi: %d\n", do_model, do_locale, do_glocale, do_localfc, do_glocalfc, do_localfi, do_glocalfi);
  if(esl_opt_GetBoolean(go, "--all"))  /* bombard the unprepared user with info */
  ///if((do_model + do_locale + do_glocale + do_localfc + do_glocalfc + do_localfi + do_glocalfi) == 0)  /* no one line options selected, do them all */
  ///do_model = do_locale = do_glocale = TRUE;
    do_model = do_locale = do_glocale = do_localfc = do_glocalfc = do_localfi = do_glocalfi = TRUE;
  int doing_locale   = FALSE;
  int doing_glocale  = FALSE;
  int doing_localfc  = FALSE;
  int doing_glocalfc = FALSE;
  int doing_localfi  = FALSE;
  int doing_glocalfi = FALSE;

  /* print general model stats (default) */
  if(do_model) { 
    ncm = 0;
    while (CMFileRead(cmfp, &abc, &cm))
      {
	if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cmfile);
	if(ncm == 0 || (esl_opt_GetBoolean(go, "--search"))) { 
	  fprintf(stdout, "#\n");
	  fprintf(stdout, "# %-4s %-20s %8s %8s %4s %5s %5s %3s %13s\n",    "",     "",                     "",         "",         "",     "",      "",      "",    " rel entropy ");
	  fprintf(stdout, "# %-4s %-20s %8s %8s %4s %5s %5s %3s %13s\n",    "",     "",                     "",         "",         "",     "",      "",      "",    "-------------");
	  fprintf(stdout, "# %-4s %-20s %8s %8s %4s %5s %5s %3s %6s %6s\n", "idx",  "name",                 "nseq",     "eff_nseq", "clen", "W",     "M",     "bif", "CM",     "HMM");
	  fprintf(stdout, "# %-4s %-20s %8s %8s %4s %5s %5s %3s %6s %6s\n", "----", "--------------------", "--------", "--------", "----", "-----", "-----", "---", "------", "------");
	}
	ncm++;

	if(! esl_opt_IsDefault(go, "--beta")) cm->beta_qdb = esl_opt_GetReal(go, "--beta");
	/* else cm->beta_qdb will be set as cm->beta_W read from CM file */
	cm->config_opts |= CM_CONFIG_QDB;
	/* update cm->config_opts */
	if(! esl_opt_GetBoolean(go, "-g"))
	  {
	    cm->config_opts |= CM_CONFIG_LOCAL;
	    cm->config_opts |= CM_CONFIG_HMMLOCAL;
	    cm->config_opts |= CM_CONFIG_HMMEL;
	  }
	ConfigCM(cm, TRUE); /* TRUE says: calculate W */

	/* print qdbs to file if nec */
	if(qdbfp != NULL) debug_print_bands(qdbfp, cm, cm->dmin, cm->dmax);

	fprintf(stdout, "%6d %-20.20s %8d %8.2f %4d %5d %5d %3d %6.2f %6.2f\n",
		ncm,
		cm->name,
		cm->nseq,
		cm->eff_nseq,
		cm->clen,
		cm->W,
		cm->M,
		CMCountStatetype(cm, B_st),
		cm_MeanMatchRelativeEntropy(cm),
		cp9_MeanMatchRelativeEntropy(cm));

	if(esl_opt_GetBoolean(go, "--search")) { if((status = summarize_search(go, errbuf, cm, r, s_w, stdout))    != eslOK) cm_Fail(errbuf); }
	FreeCM(cm);
      }    
  }
  /* print local or glocal exp tail stats if requested */
  for(i = 0; i < 2; i++) { 
    if(i == 0 && !do_locale)  continue;
    if(i == 1 && !do_glocale) continue;
    if(i == 0) { doing_locale = TRUE;  doing_glocale = FALSE; }
    if(i == 1) { doing_locale = FALSE; doing_glocale = TRUE;  }
    ncm = 0;
    seen_exps_yet = FALSE;
    CMFileRewind(cmfp);
    while (CMFileRead(cmfp, &abc, &cm)) {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cmfile);
      ncm++;
      if(cm->flags & CMH_EXPTAIL_STATS) {
	if(!seen_exps_yet) {
	  fprintf(stdout, "#\n");
	  if(doing_locale) fprintf(stdout, "# local exponential tail statistics \n");
	  else             fprintf(stdout, "# glocal exponential tail statistics \n");
	  fprintf(stdout, "#\n");
	  fprintf(stdout, "# %-4s %-15s %2s %2s %3s %11s %11s %11s %11s\n",             "",     "",                "",   "",   "",    "cyk",            "inside",         "viterbi",        "forward");
	  fprintf(stdout, "# %-4s %-15s %2s %2s %3s %11s %11s %11s %11s\n",             "",     "",                "",   "",   "",    "-----------",    "-----------",    "-----------",    "-----------");
	  fprintf(stdout, "# %-4s %-15s %2s %2s %3s %5s %5s %5s %5s %5s %5s %5s %5s\n", "idx",  "name",            "p",  "ps", "pe",  "mu",    "lmbda", "mu",    "lmbda", "mu",    "lmbda", "mu",    "lmbda");
	  fprintf(stdout, "# %-4s %-15s %2s %2s %3s %5s %5s %5s %5s %5s %5s %5s %5s\n", "----", "---------------", "--", "--", "---", "-----", "-----", "-----", "-----", "-----", "-----", "-----", "-----");
	  seen_exps_yet = TRUE;
	}
	for(p = 0; p < cm->stats->np; p++) { 
	  if(doing_locale) {
	    fprintf(stdout, "%6d %-15.15s %2d %2d %3d %5.1f %5.3f %5.1f %5.3f %5.1f %5.3f %5.1f %5.3f\n",
		   ncm,
		   cm->name,
		   p+1,
		   cm->stats->ps[p], cm->stats->pe[p],
		   cm->stats->expAA[EXP_CM_LC][p]->mu_extrap,  cm->stats->expAA[EXP_CM_LC][p]->lambda,
		   cm->stats->expAA[EXP_CM_LI][p]->mu_extrap,  cm->stats->expAA[EXP_CM_LI][p]->lambda,
		   cm->stats->expAA[EXP_CP9_LV][p]->mu_extrap, cm->stats->expAA[EXP_CP9_LV][p]->lambda,
		   cm->stats->expAA[EXP_CP9_LF][p]->mu_extrap, cm->stats->expAA[EXP_CP9_LF][p]->lambda);
	  }
	  else { /* glocal */
	    fprintf(stdout, "%6d %-15.15s %2d %2d %3d %5.1f %5.3f %5.1f %5.3f %5.1f %5.3f %5.1f %5.3f\n",
		   ncm,
		   cm->name,
		   p+1,
		   cm->stats->ps[p], cm->stats->pe[p],
		   cm->stats->expAA[EXP_CM_GC][p]->mu_extrap,  cm->stats->expAA[EXP_CM_GC][p]->lambda,
		   cm->stats->expAA[EXP_CM_GI][p]->mu_extrap,  cm->stats->expAA[EXP_CM_GI][p]->lambda,
		   cm->stats->expAA[EXP_CP9_GV][p]->mu_extrap, cm->stats->expAA[EXP_CP9_GV][p]->lambda,
		   cm->stats->expAA[EXP_CP9_GF][p]->mu_extrap, cm->stats->expAA[EXP_CP9_GF][p]->lambda);
	  }
	}
      }
      FreeCM(cm);
    }
    if(!seen_exps_yet) {
      if(doing_locale  && esl_opt_GetBoolean(go, "--le"))  cm_Fail("--le option enabled but none of the CMs in %s have exp tail stats.", cmfile);
      if(doing_glocale && esl_opt_GetBoolean(go, "--ge")) cm_Fail("--ge option enabled but none of the CMs in %s have exp tail stats.", cmfile);
      if(doing_glocale && (! esl_opt_GetBoolean(go, "--ge")))   fprintf(stdout, "# No E-value exp tail statistics.\n");
    }
  }
  /* print filter threshold stats if requested */
  for(i = 0; i < 4; i++) { /* 4 possible modes, glocal cyk, glocal inside, local cyk, local inside */
    if(i == 0 && !do_localfc)  continue;
    if(i == 1 && !do_glocalfc) continue;
    if(i == 2 && !do_localfi)  continue;
    if(i == 3 && !do_glocalfi) continue;
    if(i == 0) { doing_localfc  = TRUE;  doing_glocalfc = doing_localfi  = doing_glocalfi = FALSE; }
    if(i == 1) { doing_glocalfc = TRUE;  doing_localfc  = doing_localfi  = doing_glocalfi = FALSE; }
    if(i == 2) { doing_localfi  = TRUE;  doing_localfc  = doing_glocalfc = doing_glocalfi = FALSE; }
    if(i == 3) { doing_glocalfi = TRUE;  doing_localfc  = doing_glocalfc = doing_localfi  = FALSE; }

    ncm = 0;
    seen_fthr_yet = FALSE;
    avg_clen = avg_cm_E = avg_cm_bit_sc = avg_hmm_E = avg_hmm_bit_sc = avg_S = avg_xhmm = avg_spdup = 0.;
    tot_xhmm = tot_spdup = tot_cm_ncalcs = tot_hmm_ncalcs = tot_cm_surv_plus_fil_calcs = 0.;

    if(doing_glocalfc) { fthr_mode = FTHR_CM_GC; cm_mode = EXP_CM_GC; hmm_mode = EXP_CP9_GF; }
    if(doing_glocalfi) { fthr_mode = FTHR_CM_GI; cm_mode = EXP_CM_GI; hmm_mode = EXP_CP9_GF; }
    if(doing_localfc)  { fthr_mode = FTHR_CM_LC; cm_mode = EXP_CM_LC; hmm_mode = EXP_CP9_LF; }
    if(doing_localfi)  { fthr_mode = FTHR_CM_LI; cm_mode = EXP_CM_LI; hmm_mode = EXP_CP9_LF; }
    CMFileRewind(cmfp);
    while (CMFileRead(cmfp, &abc, &cm)) {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cmfile);
      ncm++;
      if(cm->flags & CMH_FILTER_STATS) {
	if(! (cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("cm: %d has filter threshold stats, but no exp tail stats, this shouldn't happen.");
	/* update the exp tail for the dbsize of the HMM filters */
	if((status = UpdateExpsForDBSize(cm, errbuf, dbsize)) != eslOK) cm_Fail(errbuf);
	
	/* initialize model and determine average hit length, number of CM DP calcs per residue and number of HMM DP calcs per residue */
	initialize_cm(cm, cm_mode, hmm_mode);
	if(!seen_fthr_yet) {
	  fprintf(stdout, "#\n");
	  if(doing_localfc)  fprintf(stdout, "# local CYK filter threshold stats ");
	  if(doing_glocalfc) fprintf(stdout, "# glocal CYK filter threshold stats ");
	  if(doing_localfi)  fprintf(stdout, "# local Inside filter threshold stats ");
	  if(doing_glocalfi) fprintf(stdout, "# glocal Inside filter threshold stats ");
	}  
	if(! esl_opt_GetBoolean(go, "--range")) {
	  do_avg_stats = TRUE;
	  if((! esl_opt_IsDefault(go, "-T")) || (! esl_opt_IsDefault(go, "--ga")) || (! esl_opt_IsDefault(go, "--nc")) || (! esl_opt_IsDefault(go, "--tc"))) { 
	    if(! esl_opt_IsDefault(go, "-T")) { 
	      cm_bit_sc = esl_opt_GetReal(go, "-T");
	      if (!seen_fthr_yet) fprintf(stdout, "for bit score cutoff of %4g\n#\n", esl_opt_GetReal(go, "-T"));
	    }
	    else if(! esl_opt_IsDefault(go, "--ga")) { 
	      cm_bit_sc = cm->ga;
	      if(!seen_fthr_yet) fprintf(stdout, "for Rfam GA gathering cutoff from CM file\n#\n");
	      if(! (cm->flags & CMH_GA)) ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file for cm: %d, can't use --ga.", ncm);
	    }
	    else if(! esl_opt_IsDefault(go, "--tc")) { 
	      cm_bit_sc = cm->tc;
	      if(!seen_fthr_yet) fprintf(stdout, "for Rfam TC gathering cutoff from CM file\n#\n");
	      if(! (cm->flags & CMH_TC)) ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file for cm: %d, can't use --tc.", ncm);
	    }
	    else if(! esl_opt_IsDefault(go, "--nc")) { 
	      cm_bit_sc = cm->nc;
	      if(!seen_fthr_yet) fprintf(stdout, "for Rfam NC gathering cutoff from CM file\n#\n");
	      if(! (cm->flags & CMH_NC)) ESL_FAIL(eslEINVAL, errbuf, "No NC gathering threshold in CM file for cm: %d, can't use --nc.", ncm);
	    }
	    if((status = DumpHMMFilterInfoForCMBitScore(stdout, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, ncm, cm_bit_sc, (!seen_fthr_yet),
							&cm_E, &hmm_E, &hmm_bit_sc, &S, &xhmm, &spdup, &cm_ncalcs, &hmm_ncalcs, &do_filter)) != eslOK) cm_Fail(errbuf);
	    if(do_filter) nfilter++;
	  }
	  else { /* default with --lfc,--gfc,--lfi,--gfi and ! --range is to use -E, with default value 0.1, the default cmsearch strategy */
	    cm_E = esl_opt_GetReal(go, "-E");
	    if(!seen_fthr_yet) fprintf(stdout, "for E-value cutoff of %4g per %.4f Mb\n#\n", cm_E, (dbsize / 1000000.));
	    if((status = DumpHMMFilterInfoForCME(stdout, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, ncm, cm_E, (!seen_fthr_yet), 
						 &cm_bit_sc, &hmm_E, &hmm_bit_sc, &S, &xhmm, &spdup, &cm_ncalcs, &hmm_ncalcs, &do_filter)) != eslOK) cm_Fail(errbuf);
	    if(do_filter) nfilter++;
	  }
	  avg_clen       += cm->clen;
	  avg_cm_E       += cm_E;
	  avg_cm_bit_sc  += cm_bit_sc;
	  avg_hmm_bit_sc += hmm_bit_sc;
	  avg_S          += S;
	  avg_xhmm       += xhmm;
	  avg_spdup      += spdup;
	  tot_cm_ncalcs  += cm_ncalcs;
	  tot_hmm_ncalcs += hmm_ncalcs;
	  tot_cm_surv_plus_fil_calcs += (cm_ncalcs / spdup);
	}
	else { /* --range enabled */
	  if(!seen_fthr_yet) fprintf(stdout, "for all filter cutoffs in CM file\n#\n");
	  if((status = DumpHMMFilterInfo(stdout, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, ncm)) != eslOK) cm_Fail(errbuf);
	}
	if(efp != NULL) if((status = PlotHMMFilterInfo(efp, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, FTHR_PLOT_CME_HMME)) != eslOK) cm_Fail(errbuf);
	if(bfp != NULL) if((status = PlotHMMFilterInfo(bfp, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, FTHR_PLOT_CMB_HMMB)) != eslOK) cm_Fail(errbuf);
	if(sfp != NULL) if((status = PlotHMMFilterInfo(sfp, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, (esl_opt_GetBoolean(go, "--bits") ? FTHR_PLOT_CMB_S     : FTHR_PLOT_CME_S))) != eslOK) cm_Fail(errbuf);
	if(xfp != NULL) if((status = PlotHMMFilterInfo(xfp, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, (esl_opt_GetBoolean(go, "--bits") ? FTHR_PLOT_CMB_XHMM  : FTHR_PLOT_CME_XHMM))) != eslOK) cm_Fail(errbuf);
	if(afp != NULL) if((status = PlotHMMFilterInfo(afp, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, (esl_opt_GetBoolean(go, "--bits") ? FTHR_PLOT_CMB_SPDUP : FTHR_PLOT_CME_SPDUP))) != eslOK) cm_Fail(errbuf);
	seen_fthr_yet = TRUE;
      }
      FreeCM(cm);
    }
    if(!seen_fthr_yet) { 
      if(doing_localfc  && esl_opt_GetBoolean(go, "--lfc"))  cm_Fail("--lfc option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_glocalfc && esl_opt_GetBoolean(go, "--gfc")) cm_Fail("--gfc option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_localfi  && esl_opt_GetBoolean(go, "--lfi"))  cm_Fail("--lfi option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_glocalfi && esl_opt_GetBoolean(go, "--gfi")) cm_Fail("--gfi option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_glocalfi && (! esl_opt_GetBoolean(go, "--gfi"))) fprintf(stdout, "# No filter threshold statistics.\n");
    }
    if(do_avg_stats && ncm > 1) { 
      avg_clen       /= ncm;
      avg_cm_E       /= ncm;
      avg_cm_bit_sc  /= ncm;
      avg_hmm_E      /= nfilter;
      avg_hmm_bit_sc /= nfilter;
      avg_S          /= ncm;
      avg_xhmm       /= ncm;
      avg_spdup      /= ncm;
      tot_spdup       = tot_cm_ncalcs / tot_cm_surv_plus_fil_calcs;
      tot_xhmm        = tot_cm_surv_plus_fil_calcs / tot_hmm_ncalcs;
      fprintf(stdout, "# %4s  %-15s  %5s  %8s  %6s  %6s  %6s  %7s  %7s\n", "----", "---------------", "-----", "--------", "------", "------", "------", "-------", "-------");
      fprintf(stdout, "%6s  %-15s  %5d  ", "-", "*Average*", (int) (avg_clen+0.5));
      if(avg_cm_E < 0.01)  fprintf(stdout, "%4.2e  ", avg_cm_E);
      else                 fprintf(stdout, "%8.3f  ", avg_cm_E);
      fprintf(stdout, "%6.1f  %6.1f  %6.4f  %7.1f  %7.1f\n", avg_cm_bit_sc, avg_hmm_bit_sc, avg_S, avg_xhmm, avg_spdup);

      fprintf(stdout, "%6s  %-15s  %5s  %8s  %6s  %6s  %6s  %7.1f  %7.1f\n", "-", "*Total*", "-", "-", "-", "-", "-", tot_xhmm, tot_spdup);
    }
  }
  esl_alphabet_Destroy(abc);
  printf("#\n");
  esl_stopwatch_Destroy(s_w);
  esl_randomness_Destroy(r);
  CMFileClose(cmfp);
  if(qdbfp != NULL) { printf("# Query-dependent bands saved in file %s.\n",                              esl_opt_GetString(go, "--qdbfile")); fclose(qdbfp); }
  if(efp != NULL)   { printf("# HMM filter E-val cutoff vs CM E-val cutoff plots saved in file %s.\n",   esl_opt_GetString(go, "--efile")); fclose(efp); }
  if(bfp != NULL)   { printf("# HMM filter bit sc cutoff vs CM bit sc cutoff plots saved in file %s.\n", esl_opt_GetString(go, "--bfile")); fclose(bfp); }
  if(sfp != NULL)   { printf("# Predicted survival fraction vs CM cutoff plots saved in file %s.\n",     esl_opt_GetString(go, "--sfile")); fclose(sfp); }
  if(xfp != NULL)   { printf("# Predicted xhmm (calcs * HMM) vs CM cutoff plots saved in file %s.\n",    esl_opt_GetString(go, "--xfile")); fclose(xfp); }
  if(afp != NULL)   { printf("# Predicted acceleration vs CM cutoff plots saved in file %s.\n",          esl_opt_GetString(go, "--afile")); fclose(afp); }
  esl_getopts_Destroy(go);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  esl_stopwatch_Destroy(w);
  return 0;
}

/* Function:  summarize_search()
 * Incept:    EPN, Tue Aug 21 20:00:28 2007
 *
 * Purpose:   Summarize search statistics to varying extents
 *            based on command-line options.
 */
int
summarize_search(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w, FILE *ofp) 
{
  int status;     
  int L_cm  = esl_opt_GetInteger(go, "--cmL");  /* length sequence to search with CM */
  int L_cp9 = esl_opt_GetInteger(go, "--hmmL"); /* length sequence to search with CP9 */
  float dpc;       /* number of     mega-DP calcs for search of length L */
  float dpc_q;     /* number of QDB mega-DP calcs for search of length L */
  float th_acc;    /* theoretical QDB acceleration */
  float dpc_v;     /* number of CP9 mega-DP calcs for search of length L */

  /* optional, -t related variables */
  ESL_DSQ *dsq_cm; /* digitized sequence of length L for CM  timings  */
  ESL_DSQ *dsq_cp9;/* digitized sequence of length L for CP9 timings  */
  float t_c;       /* number of seconds (w->user) for        CYK search */
  float t_i;       /* number of seconds (w->user) for     Inside search */
  float t_cq;      /* number of seconds (w->user) for QDB    CYK search */
  float t_iq;      /* number of seconds (w->user) for QDB Inside search */
  float t_v;       /* number of seconds (w->user) for CP9 Viterbi search */
  float t_f;       /* number of seconds (w->user) for CP9 Forward search */

  if(L_cm < cm->W) { L_cm = cm->W; fprintf(stdout, "\tL increased to minimum size of cm->W (%d)\n", L_cm); }
  ESL_ALLOC(dsq_cm,  sizeof(ESL_DSQ) * (L_cm +2));
  ESL_ALLOC(dsq_cp9, sizeof(ESL_DSQ) * (L_cp9+2));
  esl_rsq_xfIID(r, cm->null, cm->abc->K, L_cm,  dsq_cm);
  esl_rsq_xfIID(r, cm->null, cm->abc->K, L_cp9, dsq_cp9);

  /* estimate speedup due to QDB */
  dpc    = CountScanDPCalcs(cm, L_cm, FALSE) / 1000000.;
  dpc_q  = CountScanDPCalcs(cm, L_cm, TRUE)  / 1000000.;
  th_acc = dpc / dpc_q;

  if(esl_opt_GetBoolean(go, "-g")) dpc_v  = (float) (cm->clen+1) * L_cp9 * 9;  /*  9 transitions queried per HMM node: 9 main model no local begin, end, nor EL */
  else                             dpc_v  = (float) (cm->clen+1) * L_cp9 * 12; /* 12 transitions queried per HMM node: 9 main model, local begin, end, and EL */
  dpc_v /= 1000000.;

  /* First create scan info for non-QDB runs */
  int *tmp_dmin = cm->dmin;
  int *tmp_dmax = cm->dmax;
  cm->dmin = NULL;
  cm->dmax = NULL;
  cm->search_opts |= CM_SEARCH_NOQDB;
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("summarize_search(), CreateScanMatrixForCM() call failed.");
  
  /* cyk */
  esl_stopwatch_Start(w);
  if((status = FastCYKScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, TRUE, NULL, NULL)) != eslOK) goto ERROR;
  /*CYKScan (cm, dsq_cm, 1, L_cm, cm->W, 0., NULL);*/
  esl_stopwatch_Stop(w);
  t_c = w->user;

  /* inside */
  cm->search_opts |= CM_SEARCH_INSIDE;
  esl_stopwatch_Start(w);
  if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, TRUE, NULL, NULL)) != eslOK) goto ERROR;
  /* iInsideScan (cm, dsq_cm, 1, L_cm, cm->W, 0., NULL); */
  esl_stopwatch_Stop(w);
  t_i = w->user;

  /* reset cm->dmin, cm->dmax, recalc scanmatrix */
  cm->dmin = tmp_dmin;
  cm->dmax = tmp_dmax;
  cm->search_opts &= ~CM_SEARCH_NOQDB;
  cm->search_opts &= ~CM_SEARCH_INSIDE;
  cm_FreeScanMatrixForCM(cm);
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("summarize_search(), CreateScanMatrix() call failed.");

  /* qdb cyk */
  esl_stopwatch_Start(w);
  if((status = FastCYKScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, TRUE, NULL, NULL)) != eslOK) goto ERROR;
  /*CYKBandedScan (cm, dsq_cm, cm->dmin, cm->dmax, 1, L_cm, cm->W, 0., NULL); */
  esl_stopwatch_Stop(w);
  t_cq = w->user;

  /* qdb inside */
  cm->search_opts |= CM_SEARCH_INSIDE;
  esl_stopwatch_Start(w);
  if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, TRUE, NULL, NULL)) != eslOK) goto ERROR;
  /*iInsideBandedScan (cm, dsq_cm, cm->dmin, cm->dmax, 1, L_cm, cm->W, 0., NULL);*/
  esl_stopwatch_Stop(w);
  t_iq = w->user;
  
  /* CP9 viterbi */
  esl_stopwatch_Start(w);
  if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq_cp9, 1, L_cp9, cm->W, 0., NULL,
			   TRUE,   /* we're scanning */
			   FALSE,  /* we're not ultimately aligning */
			   TRUE,   /* be memory efficient */
			   TRUE,   /* do NULL3 score corrections, for accurate timings */
			   NULL, NULL,
			   NULL,   /* don't want traces back */
			   NULL)) != eslOK) goto ERROR;

  esl_stopwatch_Stop(w);
  t_v = w->user;

  /* CP9 forward */
  esl_stopwatch_Start(w);
  if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq_cp9, 1, L_cp9, cm->W, 0., NULL,
			   TRUE,   /* we're scanning */
			   FALSE,  /* we're not ultimately aligning */
			   TRUE,   /* be memory efficient */
			   TRUE,   /* do NULL3 score corrections, for accurate timings */
			   NULL, NULL, NULL)) != eslOK) goto ERROR;
  esl_stopwatch_Stop(w);
  t_f = w->user;
  
  /* Experimental method of accelerating Forward by determining guarantees that scores will be non -INFTY,
   * see cp9_dp.c:cp9_FastForward for more. 
   * 
   * determine the minimum length we can search safely with the optimized forward implementation. 
   * int   minL = 0;  minimum length can safely scan with optimized Forward(), -1 ==> any length 
   *if((status = cp9_WorstForward(cm, errbuf, cm->cp9_mx, -INFTY, TRUE, FALSE, &minL)) != eslOK) goto ERROR;
   *int be_safe = FALSE; should we be safe, and not use optimized Forward()? 
   * ESL_DPRINTF1(("minL: %d L: %d\n", minL, L));
   *if(minL != -1 && minL <= L) be_safe = TRUE;
   *esl_stopwatch_Start(w);
   *if((status = cp9_FastForward(cm, errbuf, cm->cp9_mx, dsq_cp9, 1, L_cp9, cm->W, 0., NULL, 
   * TRUE,   
   * FALSE,
   * TRUE,
   * be_safe,
   * NULL, NULL, NULL)) != eslOK) goto ERROR;
   */

  /* print results */
  float mc_s; /* million calcs / second */
  float kb_s; /* kilobases / second */
  float emp_acc; /* empirical acceleration from QDB */
  float L_cm_kb  = (float) L_cm / 1000.;
  float L_cp9_kb = (float) L_cp9 / 1000.;
  fprintf(stdout, "#\n");
  fprintf(stdout, "#\t\t\t search statistics:\n");
  fprintf(stdout, "#\t\t\t %7s %7s %6s %8s   %5s %5s %5s\n",           "alg",     "Mc/kb",   "Mc/s",   "kb/s",     "beta",   "qdbXt", "qdbXe");
  fprintf(stdout, "#\t\t\t %7s %7s %6s %8s   %5s %5s %5s\n",           "-------", "-------", "------", "--------", "-----",  "-----", "-----");
  /* cyk non-banded */
  float dpc_kb = dpc * (1000. / (float) L_cm); /* convert to cells per KB */
  mc_s = dpc / t_c; 
  kb_s = ((float) L_cm_kb) / t_c; 
  fprintf(stdout, " \t\t\t %7s %7.1f %6.1f %8.2f   %5s %5s %5s\n", "cyk",     dpc_kb, mc_s, kb_s, "-", "-", "-");
  mc_s = dpc / t_i; 
  kb_s = ((float) L_cm_kb) / t_i; 
  fprintf(stdout, " \t\t\t %7s %7.1f %6.1f %8.2f   %5s %5s %5s\n", "inside",  dpc_kb, mc_s, kb_s, "-", "-", "-");
  float dpc_q_kb = dpc_q * (1000. / (float) L_cm); /* convert to cells per KB */
  mc_s = dpc_q / t_cq; 
  kb_s = ((float) L_cm_kb) / t_cq; 
  emp_acc = t_c / t_cq; 
  fprintf(stdout, " \t\t\t %7s %7.1f %6.1f %8.2f   %5g %5.1f %5.1f\n", "cyk",     dpc_q_kb, mc_s, kb_s, cm->beta_qdb, th_acc, emp_acc);
  mc_s = dpc_q / t_iq; 
  kb_s = ((float) L_cm_kb) / t_iq; 
  emp_acc = t_i / t_iq; 
  fprintf(stdout, " \t\t\t %7s %7.1f %6.1f %8.2f   %5g %5.1f %5.1f\n", "inside",  dpc_q_kb, mc_s, kb_s, cm->beta_qdb, th_acc, emp_acc);
  mc_s = dpc_v / t_v; 
  kb_s = ((float) L_cp9_kb) / t_v; 
  float dpc_v_kb = dpc_v * (1000. / (float) L_cp9); /* convert to cells per KB */
  fprintf(stdout, " \t\t\t %7s %7.1f %6.1f %8.2f   %5s %5s %5s\n", "viterbi",  dpc_v_kb, mc_s, kb_s, "-", "-", "-");
  mc_s = dpc_v / t_f; 
  kb_s = ((float) L_cp9_kb) / t_f; 
  fprintf(stdout, " \t\t\t %7s %7.1f %6.1f %8.2f   %5s %5s %5s\n", "forward",  dpc_v_kb, mc_s, kb_s, "-", "-", "-");
  
  free(dsq_cm);
  free(dsq_cp9);
  return eslOK;

 ERROR:
  return status; 
}

/* initialize_cm()
 * Setup the CM based on the exp tail mode,
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(CM_t *cm, int cm_mode, int hmm_mode)
{
  /* Update cm->config_opts based on exp tail mode */
  if(ExpModeIsLocal(cm_mode))  cm->config_opts |= CM_CONFIG_LOCAL;
  if(ExpModeIsLocal(hmm_mode)) {
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  ConfigCM(cm, TRUE); /* TRUE says: calculate W */

  return eslOK;
}

/* Function: print_run_info
 * Date:     EPN, Mon Mar  3 09:47:26 2008
 *
 * Purpose:  Print information on this run of cmstat.
 *           Command used to run it, and execution date.
 *
 * Returns:  eslOK on success
 */
static int
print_run_info(const ESL_GETOPTS *go, char *errbuf, ESL_RANDOMNESS *r)
{
  int status;
  char *command;
  char *date;

  if((status = get_command(go, errbuf, &command)) != eslOK) return status;
  if((status = GetDate    (errbuf, &date))    != eslOK) return status;

  fprintf(stdout, "%-10s %s\n",  "# command:", command);
  fprintf(stdout, "%-10s %s\n",  "# date:",    date);
  if(esl_opt_GetBoolean(go, "--search")) fprintf(stdout, "%-10s %ld\n", "# seed:",    esl_randomness_GetSeed(r));

  free(command);
  free(date);
  return eslOK;
}

/* Function: get_command
 * Date:     EPN, Mon Mar  3 09:48:55 2008
 *
 * Purpose:  Return the command used to call cmstat
 *           in <ret_command>.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
int 
get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command)
{
  int status;
  int i;
  char *command = NULL;

  for (i = 0; i < go->argc; i++) { /* copy all command line options and args */
    if((status = esl_strcat(&(command),  -1, go->argv[i], -1)) != eslOK) goto ERROR;
    if(i < (go->argc-1)) if((status = esl_strcat(&(command), -1, " ", 1)) != eslOK) goto ERROR;
  }
  *ret_command = command;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_command(): memory allocation error.");
  return status;
}

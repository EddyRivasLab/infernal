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
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_sqio.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

#define ONELINEOPTS  "-m,--le,--ge,--lfc,--gfc,--lfi,--gfi" /* exclusive choice of summary stats */
#define NOTMOPTS     "--le,--ge,--lfc,--gfc,--lfi,--gfi"    /* incompatible with -g, --beta, and --search */

static ESL_OPTIONS options[] = {
  /* name           type      default      env  range     toggles      reqs     incomp    help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL,      NULL,   "show brief help on version and usage",   1 },
  { "-g",        eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, NOTMOPTS,    "configure CM for glocal alignment [default: local]", 1 },
  { "-m",        eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of model statistics", 1 },
  { "--le",      eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of  local E-value statistics", 1 },
  { "--ge",      eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of glocal E-value statistics", 1 },
  { "--lfc",     eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of  local CYK    filter threshold stats", 1 },
  { "--gfc",     eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of glocal CYK    filter threshold stats", 1 },
  { "--lfi",     eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of  local Inside filter threshold stats", 1 },
  { "--gfi",     eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, ONELINEOPTS, "only print one line summary of glocal Inside filter threshold stats", 1 },
  { "--seqfile", eslARG_INFILE, FALSE,     NULL, NULL,      NULL,      NULL,      NULL,   "for filter stats: compute E-value cutoffs for sequence file <f>", 1 },
  { "--toponly", eslARG_NONE,   FALSE,     NULL, NULL,      NULL,"--seqfile",     NULL,   "with --seqfile, only consider top-strand", 1},
  { "--beta",    eslARG_REAL,   "1E-7",    NULL, "0<x<1",   NULL,      NULL, NOTMOPTS,    "set tail loss prob for W calc and QDB stats to <x>", 1 },
  { "--search",  eslARG_NONE,   FALSE,     NULL, NULL,      NULL,      NULL, NOTMOPTS,    "do search timing experiments", 1 },
  { "--cmL",     eslARG_INT,    "1000",    NULL, "n>0",     NULL,"--search",      NULL,   "length of sequences for CM search stats", 1 },
  { "--hmmL",    eslARG_INT,    "100000",  NULL, "n>0",     NULL,"--search",      NULL,   "length of sequences for CP9 HMM search stats", 1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


static char usage[]  = "[-options] <cmfile>";
static char banner[] = "display summary statistics for CMs";

static int    summarize_search(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w); 
/*static int    summarize_alignment(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w); */
static float  count_align_dp_calcs(CM_t *cm, int L);
static int    initialize_cm(CM_t *cm, int cm_mode, int hmm_mode);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing   */
  ESL_ALPHABET    *abc = NULL;  /* alphabet                  */
  ESL_RANDOMNESS  *r   = NULL;  /* source of randomness      */
  ESL_STOPWATCH   *w   = NULL;  /* for timings               */
  char            *cmfile;	/* name of input CM file     */ 
  CMFILE          *cmfp;	/* open input CM file stream */
  CM_t            *cm;          /* CM most recently read     */
  int              ncm;         /* CM index                  */
  char             errbuf[cmERRBUFSIZE];
  int              status;
  int seen_gumbels_yet = FALSE;
  int seen_fthr_yet    = FALSE;
  int p, i;
  int   fthr_mode;
  int   cm_mode, hmm_mode;
  long dbsize;

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
  cm_banner(stdout, argv[0], banner);

  cmfile     = esl_opt_GetArg(go, 1); 
  r = esl_randomness_CreateTimeseeded();
  w = esl_stopwatch_Create();
  if(r == NULL) cm_Fail("Failed to create RNG, probably out of memory.\n");
  if(w == NULL) cm_Fail("Failed to create stopwatch, probably out of memory.\n");

  /* Initializations: open the CM file
   */
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    cm_Fail("Failed to open covariance model save file %s\n", cmfile);

  /* if --seqfile enabled, read the sequence file to get the database size, else we use 1 Mb as db size */
  if(!(esl_opt_IsDefault(go, "--seqfile"))) { 
    if((!(esl_opt_GetBoolean(go, "--lfc"))) && (!(esl_opt_GetBoolean(go, "--lfi"))) && (!(esl_opt_GetBoolean(go, "--gfc"))) && (!(esl_opt_GetBoolean(go, "--gfi"))))
      cm_Fail("--seqfile only makes sense with one of: --lfc, --lfi, --gfc, or --gfi");
    /* open input sequence file */
    ESL_SQFILE      *sqfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--seqfile"), eslSQFILE_UNKNOWN, NULL, &sqfp);
    if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", esl_opt_GetString(go, "--seqfile"));
    else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", esl_opt_GetString(go, "--seqfile"));
    else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
    else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
    /* GetDBInfo() reads all sequences, rewinds seq file and returns db size */
    GetDBInfo(NULL, sqfp, &(dbsize), NULL);  
    esl_sqfile_Close(sqfp); 
    if (! esl_opt_GetBoolean(go, "--toponly")) dbsize *= 2;
  }
  else dbsize = FTHR_DBSIZE; /* 1 Mb */

  /* Main body: read CMs one at a time, print stats 
   * Options for only printing 1 line per CM: 
   * if -m:              print general model stats, and optionally determine search stats (if --search)
   * else if  --le:  print  local gumbel stats
   * else if --ge:  print glocal gumbel stats
   * else if  --lfc: print  local CYK    filter threshold stats
   * else if --gfc: print glocal CYK    filter threshold stats
   * else if  --lfc: print  local Inside filter threshold stats
   * else if --gfc: print glocal Inside filter threshold stats
   * else (default):     print all stat categories, one category at a time 
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
  if((do_model + do_locale + do_glocale + do_localfc + do_glocalfc + do_localfi + do_glocalfi) == 0)  /* no one line options selected, do them all */
    do_model = do_locale = do_glocale = TRUE;
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
	if(ncm == 0 || (esl_opt_GetBoolean(go, "--search"))) { 
	  printf("#\n");
	  printf("# %-4s %-20s %8s %8s %4s %5s %5s %3s %13s\n",    "",     "",                     "",         "",         "",     "",      "",      "",    "    relent   ");
	  printf("# %-4s %-20s %8s %8s %4s %5s %5s %3s %13s\n",    "",     "",                     "",         "",         "",     "",      "",      "",    "-------------");
	  printf("# %-4s %-20s %8s %8s %4s %5s %5s %3s %6s %6s\n", "idx",  "name",                 "nseq",     "eff_nseq", "clen", "W",     "M",     "bif", "CM",     "HMM");
	  printf("# %-4s %-20s %8s %8s %4s %5s %5s %3s %6s %6s\n", "----", "--------------------", "--------", "--------", "----", "-----", "-----", "---", "------", "------");
	}
	ncm++;


	cm->beta = esl_opt_GetReal(go, "--beta");
	cm->config_opts |= CM_CONFIG_QDB;
	/* update cm->config_opts */
	if(! esl_opt_GetBoolean(go, "-g"))
	  {
	    cm->config_opts |= CM_CONFIG_LOCAL;
	    cm->config_opts |= CM_CONFIG_HMMLOCAL;
	    cm->config_opts |= CM_CONFIG_HMMEL;
	  }
	/*ConfigQDB(cm);*/
	ConfigCM(cm, NULL, NULL);

	printf("%6d %-20s %8d %8.2f %4d %5d %5d %3d %6.2f %6.2f\n",
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

	if(esl_opt_GetBoolean(go, "--search")) { if((status = summarize_search(go, errbuf, cm, r, w))    != eslOK) cm_Fail(errbuf); }
	FreeCM(cm);
      }    
  }
  /* print local or glocal gumbel stats if requested */
  for(i = 0; i < 2; i++) { 
    if(i == 0 && !do_locale)  continue;
    if(i == 1 && !do_glocale) continue;
    if(i == 0) { doing_locale = TRUE;  doing_glocale = FALSE; }
    if(i == 1) { doing_locale = FALSE; doing_glocale = TRUE;  }
    ncm = 0;
    seen_gumbels_yet = FALSE;
    CMFileRewind(cmfp);
    while (CMFileRead(cmfp, &abc, &cm)) {
      ncm++;
      if(cm->flags & CMH_GUMBEL_STATS) {
	if(!seen_gumbels_yet) {
	  printf("#\n");
	  if(doing_locale) printf("# local gumbel statistics \n");
	  else             printf("# glocal gumbel statistics \n");
	  printf("#\n");
	  printf("# %-4s %-15s %2s %2s %3s %11s %11s %11s %11s\n",             "",     "",                "",   "",   "",    "cyk",            "inside",         "viterbi",        "forward");
	  printf("# %-4s %-15s %2s %2s %3s %11s %11s %11s %11s\n",             "",     "",                "",   "",   "",    "-----------",    "-----------",    "-----------",    "-----------");
	  printf("# %-4s %-15s %2s %2s %3s %5s %5s %5s %5s %5s %5s %5s %5s\n", "idx",  "name",            "p",  "ps", "pe",  "mu",    "lmbda", "mu",    "lmbda", "mu",    "lmbda", "mu",    "lmbda");
	  printf("# %-4s %-15s %2s %2s %3s %5s %5s %5s %5s %5s %5s %5s %5s\n", "----", "---------------", "--", "--", "---", "-----", "-----", "-----", "-----", "-----", "-----", "-----", "-----");
	  seen_gumbels_yet = TRUE;
	}
	for(p = 0; p < cm->stats->np; p++) { 
	  if(doing_locale) {
	    printf("%6d %-15s %2d %2d %3d %5.2f %5.3f %5.2f %5.3f %5.2f %5.3f %5.2f %5.3f\n",
		   ncm,
		   cm->name,
		   p+1,
		   cm->stats->ps[p], cm->stats->pe[p],
		   cm->stats->gumAA[GUM_CM_LC][p]->mu,  cm->stats->gumAA[GUM_CM_LC][p]->lambda,
		   cm->stats->gumAA[GUM_CM_LI][p]->mu,  cm->stats->gumAA[GUM_CM_LI][p]->lambda,
		   cm->stats->gumAA[GUM_CP9_LV][p]->mu, cm->stats->gumAA[GUM_CP9_LV][p]->lambda,
		   cm->stats->gumAA[GUM_CP9_LF][p]->mu, cm->stats->gumAA[GUM_CP9_LF][p]->lambda);
	  }
	  else { /* glocal */
	    printf("%6d %-15s %2d %2d %3d %5.2f %5.3f %5.2f %5.3f %5.2f %5.3f %5.2f %5.3f\n",
		   ncm,
		   cm->name,
		   p+1,
		   cm->stats->ps[p], cm->stats->pe[p],
		   cm->stats->gumAA[GUM_CM_GC][p]->mu,  cm->stats->gumAA[GUM_CM_GC][p]->lambda,
		   cm->stats->gumAA[GUM_CM_GI][p]->mu,  cm->stats->gumAA[GUM_CM_GI][p]->lambda,
		   cm->stats->gumAA[GUM_CP9_GV][p]->mu, cm->stats->gumAA[GUM_CP9_GV][p]->lambda,
		   cm->stats->gumAA[GUM_CP9_GF][p]->mu, cm->stats->gumAA[GUM_CP9_GF][p]->lambda);
	  }
	}
      }
      FreeCM(cm);
    }
    if(!seen_gumbels_yet) {
      if(doing_locale  && esl_opt_GetBoolean(go, "--le"))  cm_Fail("--le option enabled but none of the CMs in %s have Gumbel stats.", cmfile);
      if(doing_glocale && esl_opt_GetBoolean(go, "--ge")) cm_Fail("--ge option enabled but none of the CMs in %s have Gumbel stats.", cmfile);
      if(doing_glocale && (! esl_opt_GetBoolean(go, "--ge")))   printf("# No E-value Gumbel statistics.\n");
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

    if(doing_glocalfc) { fthr_mode = FTHR_CM_GC; cm_mode = GUM_CM_GC; hmm_mode = GUM_CP9_GF; }
    if(doing_glocalfi) { fthr_mode = FTHR_CM_GI; cm_mode = GUM_CM_GI; hmm_mode = GUM_CP9_GF; }
    if(doing_localfc)  { fthr_mode = FTHR_CM_LC; cm_mode = GUM_CM_LC; hmm_mode = GUM_CP9_LF; }
    if(doing_localfi)  { fthr_mode = FTHR_CM_LI; cm_mode = GUM_CM_LI; hmm_mode = GUM_CP9_LF; }
    CMFileRewind(cmfp);
    while (CMFileRead(cmfp, &abc, &cm)) {
      ncm++;
      if(cm->flags & CMH_FILTER_STATS) {
	if(! (cm->flags & CMH_GUMBEL_STATS)) cm_Fail("cm: %d has filter threshold stats, but no Gumbel stats, this shouldn't happen.");
	/* we expect theoretical db size used to calc filter threshold stats is FTHR_DBSIZE (1 Mb) */
	if(cm->stats->hfiA[fthr_mode]->dbsize != FTHR_DBSIZE)
	  cm_Fail("Expected db size of %d in cm file for filter thr stats, but read db size of %d residues.", FTHR_DBSIZE, cm->stats->hfiA[fthr_mode]->dbsize);
	/* update the Gumbels for the dbsize of the HMM filters */
	if((status = UpdateGumbelsForDBSize(cm, errbuf, dbsize)) != eslOK) cm_Fail(errbuf);
	
	/* initialize model and determine average hit length, number of CM DP calcs per residue and number of HMM DP calcs per residue */
	initialize_cm(cm, cm_mode, hmm_mode);
	if(!seen_fthr_yet) {
	  printf("#\n");
	  if(doing_localfc)  printf("# local CYK filter threshold statistics\n");
	  if(doing_glocalfc) printf("# glocal CYK filter threshold statistics\n");
	  if(doing_localfi)  printf("# local Inside filter threshold statistics\n");
	  if(doing_glocalfi) printf("# glocal Inside filter threshold statistics\n");
	  printf("#\n");
	  seen_fthr_yet = TRUE;
	}  
	if((status = DumpHMMFilterInfo(stdout, cm->stats->hfiA[fthr_mode], errbuf, cm, cm_mode, hmm_mode, dbsize, ncm)) != eslOK) cm_Fail(errbuf);
      }
      FreeCM(cm);
    }
    if(!seen_fthr_yet) { 
      if(doing_localfc  && esl_opt_GetBoolean(go, "--lfc"))  cm_Fail("--lfc option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_glocalfc && esl_opt_GetBoolean(go, "--gfc")) cm_Fail("--gfc option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_localfi  && esl_opt_GetBoolean(go, "--lfi"))  cm_Fail("--lfi option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_glocalfi && esl_opt_GetBoolean(go, "--gfi")) cm_Fail("--gfi option enabled but none of the CMs in %s have filter threshold stats.", cmfile);
      if(doing_glocalfi && (! esl_opt_GetBoolean(go, "--gfi"))) printf("# No filter threshold statistics.\n");
    }
  }
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  CMFileClose(cmfp);
  esl_getopts_Destroy(go);
  return 0;
}

/* Function:  summarize_search()
 * Incept:    EPN, Tue Aug 21 20:00:28 2007
 *
 * Purpose:   Summarize search statistics to varying extents
 *            based on command-line options.
 */
int
summarize_search(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w) 
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

  if(L_cm < cm->W) { L_cm = cm->W; printf("\tL increased to minimum size of cm->W (%d)\n", L_cm); }
  ESL_ALLOC(dsq_cm,  sizeof(ESL_DSQ) * L_cm +2);
  ESL_ALLOC(dsq_cp9, sizeof(ESL_DSQ) * L_cp9+2);
  esl_rnd_xfIID(r, cm->null, cm->abc->K, L_cm,  dsq_cm);
  esl_rnd_xfIID(r, cm->null, cm->abc->K, L_cp9, dsq_cp9);

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
  if((status = FastCYKScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
  /*CYKScan (cm, dsq_cm, 1, L_cm, cm->W, 0., NULL);*/
  esl_stopwatch_Stop(w);
  t_c = w->user;

  /* inside */
  cm->search_opts |= CM_SEARCH_INSIDE;
  esl_stopwatch_Start(w);
  if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
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
  if((status = FastCYKScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
  /*CYKBandedScan (cm, dsq_cm, cm->dmin, cm->dmax, 1, L_cm, cm->W, 0., NULL); */
  esl_stopwatch_Stop(w);
  t_cq = w->user;

  /* qdb inside */
  cm->search_opts |= CM_SEARCH_INSIDE;
  esl_stopwatch_Start(w);
  if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq_cm, 1, L_cm, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
  /*iInsideBandedScan (cm, dsq_cm, cm->dmin, cm->dmax, 1, L_cm, cm->W, 0., NULL);*/
  esl_stopwatch_Stop(w);
  t_iq = w->user;
  
  /* CP9 viterbi */
  esl_stopwatch_Start(w);
  if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq_cp9, 1, L_cp9, cm->W, 0., NULL,
			   TRUE,   /* we're scanning */
			   FALSE,  /* we're not ultimately aligning */
			   TRUE,   /* be memory efficient */
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
  printf("#\n");
  printf("#\t\t\t search statistics:\n");
  printf("#\t\t\t %7s %6s %6s %8s   %5s %5s %5s\n",           "alg",     "Mc/kb",  "Mc/s",   "kb/s",     "beta",   "qdbXt", "qdbXe");
  printf("#\t\t\t %7s %6s %6s %8s   %5s %5s %5s\n",           "-------", "------", "------", "--------", "-----",  "-----", "-----");
  /* cyk non-banded */
  float dpc_kb = dpc * (1000. / (float) L_cm); /* convert to cells per KB */
  mc_s = dpc / t_c; 
  kb_s = ((float) L_cm_kb) / t_c; 
  printf(" \t\t\t %7s %6.1f %6.1f %8.2f   %5s %5s %5s\n", "cyk",     dpc_kb, mc_s, kb_s, "", "", "");
  mc_s = dpc / t_i; 
  kb_s = ((float) L_cm_kb) / t_i; 
  printf(" \t\t\t %7s %6.1f %6.1f %8.2f   %5s %5s %5s\n", "inside",  dpc_kb, mc_s, kb_s, "", "", "");
  float dpc_q_kb = dpc_q * (1000. / (float) L_cm); /* convert to cells per KB */
  mc_s = dpc_q / t_cq; 
  kb_s = ((float) L_cm_kb) / t_cq; 
  emp_acc = t_c / t_cq; 
  printf(" \t\t\t %7s %6.1f %6.1f %8.2f   %5g %5.1f %5.1f\n", "cyk",     dpc_q_kb, mc_s, kb_s, cm->beta, th_acc, emp_acc);
  mc_s = dpc_q / t_iq; 
  kb_s = ((float) L_cm_kb) / t_iq; 
  emp_acc = t_i / t_iq; 
  printf(" \t\t\t %7s %6.1f %6.1f %8.2f   %5g %5.1f %5.1f\n", "inside",  dpc_q_kb, mc_s, kb_s, cm->beta, th_acc, emp_acc);
  mc_s = dpc_v / t_v; 
  kb_s = ((float) L_cp9_kb) / t_v; 
  float dpc_v_kb = dpc_v * (1000. / (float) L_cp9); /* convert to cells per KB */
  printf(" \t\t\t %7s %6.1f %6.1f %8.2f   %5s %5s %5s\n", "viterbi",  dpc_v_kb, mc_s, kb_s, "", "", "");
  mc_s = dpc_v / t_f; 
  kb_s = ((float) L_cp9_kb) / t_f; 
  printf(" \t\t\t %7s %6.1f %6.1f %8.2f   %5s %5s %5s\n", "forward",  dpc_v_kb, mc_s, kb_s, "", "", "");
  
  free(dsq_cm);
  free(dsq_cp9);
  return eslOK;

 ERROR:
  return status; 
}

/* Function:  summarize_alignment()
 * Incept:    
 *
 * Purpose:   Summarize alignment statistics to varying extents
 *            based on command-line options.
 */
int
summarize_alignment(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w) 
{
  /* HERE: do HMM banded alignment stats
   * sample N=100 seqs, and calculate posteriors, determine new
   * number of CYK DP calcs AND CP9 F/B calcs to get bands. */
  int status;
  float dpc;  /* # DP calcs for non-banded alignment of consensus */
  CMConsensus_t *con = NULL;            /* consensus info for the CM */
  ESL_SQ *csq = NULL;
  float t_dc; /* user seconds time for D&C alignment */
  float t_hb; /* user seconds time for HMM banded alignment */
  float mc_s; /* million calcs/second */
  float size_limit = esl_opt_GetReal(go, "--mxsize");

  /* Create and align consensus sequence for D&C stats */
  CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &con);
  if((csq = esl_sq_CreateFrom("consensus", con->cseq, NULL, NULL, NULL)) == NULL)
    { status = eslEMEM; goto ERROR; }
  esl_sq_Digitize(cm->abc, csq);
  dpc = count_align_dp_calcs(cm, csq->n) / 1000000.;

  /* cyk inside (score only) */
  esl_stopwatch_Start(w);
  /*CYKDivideAndConquer(cm, csq->dsq, csq->n, 0, 1, csq->n, NULL, NULL, NULL);*/
  CYKInsideScore(cm, csq->dsq, csq->n, 0, 1, csq->n, 
		 NULL, NULL); /* don't do QDB mode */
  esl_stopwatch_Stop(w);
  t_dc = w->user;
  mc_s = dpc / t_dc;

  /* HMM banded */
  /* Emit N seqs, and align them, to get total time up to reasonable level,
   * and to average out tightness of bands */
  int N = esl_opt_GetInteger(go, "-N");
  seqs_to_aln_t *seqs_to_aln = NULL;
  ESL_SQ **sq = NULL;
  ESL_ALLOC(sq, sizeof(ESL_SQ *) * N);
  int L; 
  float L_avg = 0.; 
  int i;
  for(i = 0; i < N; i++)
    {
      if((status = EmitParsetree(cm, errbuf, r, "seq", TRUE, NULL, &(sq[i]), &L)) != eslOK) goto ERROR;
      /*esl_sqio_Write(stdout, sq[i], eslSQFILE_FASTA);*/
      L_avg += L;
    }
  L_avg /= (float) N;
  cm->align_opts |= CM_ALIGN_HBANDED;
  esl_stopwatch_Start(w);
  seqs_to_aln = CreateSeqsToAlnFromSq(sq, N, FALSE);
  if((status = DispatchAlignments(cm, errbuf, seqs_to_aln, NULL, NULL, 0, 0, 0, TRUE, NULL, size_limit)) != eslOK) goto ERROR;
  esl_stopwatch_Stop(w);
  t_hb = w->user / (float) N;
  FreeSeqsToAln(seqs_to_aln);

  printf("#\n");
  printf("#\t\t\t Alignment statistics:\n");
  printf("#\t\t\t %7s %6s %6s %8s %8s %8s\n",             "alg",     "Mc",     "L",     "Mc/s",    "s/seq",        "accel");
  printf("#\t\t\t %7s %6s %6s %8s %8s %8s\n",             "-------", "------", "------","--------", "--------", "--------");
  /*mc_s = dpc / t_dc; */
  printf(" \t\t\t %7s %6.1f %6d %8.1f %8.3f %8s\n",       "cyk",      dpc,      csq->n,  mc_s,      t_dc,       "-");
  printf(" \t\t\t %7s %6s %6.0f %8s %8.3f %8.2f\n",       "hb cyk",   "?",      L_avg,   "?",       t_hb,       (t_dc/t_hb));

  esl_sq_Destroy(csq);
  FreeCMConsensus(con);
  return eslOK;
 ERROR:
  cm_Fail("ERROR code %d in summarize_stats().", status);
  return status; /* NOTREACHED */
}

/* Function: count_align_dp_calcs()
 * Date:     EPN, Wed Aug 22 09:08:03 2007
 *
 * Purpose:  Count all non-d&c inside DP calcs for a CM 
 *           alignment of a seq of length L. Similar to cm_dpsmall.c's
 *           CYKDemands() but takes into account number of
 *           transitions from each state, and is concerned
 *           with a scanning dp matrix, not an alignment matrix.
 *
 * Args:     cm     - the model
 *           L      - length of sequence
 *
 * Returns: (float) the total number of DP calculations.
 */
float count_align_dp_calcs(CM_t *cm, int L)
{
  int v, j;
  float dpcalcs = 0.;
  float dpcalcs_bif = 0.;
  
  float  dpcells     = 0.;
  float  dpcells_bif = 0.;

  dpcells = (L+2) * (L+1) * 0.5; /* fillable dp cells per state (deck) */
  for (j = 0; j <= L; j++)
    dpcells_bif += (j+2) * (j+1) * .5;
  dpcalcs_bif = CMCountStatetype(cm, B_st) * dpcells_bif; /* no choice of transitions */
  for(v = 0; v < cm->M; v++)
    if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
      dpcalcs += dpcells * cm->cnum[v]; /* cnum choices of transitions */

  return dpcalcs + dpcalcs_bif;
}

/* initialize_cm()
 * Setup the CM based on the Gumbel mode,
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(CM_t *cm, int cm_mode, int hmm_mode)
{
  /* Update cm->config_opt based on gumbel mode */
  if(GumModeIsLocal(cm_mode))  cm->config_opts |= CM_CONFIG_LOCAL;
  if(GumModeIsLocal(hmm_mode)) {
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  ConfigCM(cm, NULL, NULL);

  return eslOK;
}

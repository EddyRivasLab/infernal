/* cmcalibrate.c
 * Score a CM and a CM Plan 9 HMM against random sequence 
 * data to set the statistical parameters for E-value determination,
 * and CP9 HMM filtering thresholds. 
 * 
 * EPN, Wed May  2 07:02:52 2007
 * based on HMMER-2.3.2's hmmcalibrate.c from SRE
 *
 * MPI example:  
 * qsub -N testrun -j y -R y -b y -cwd -V -pe lam-mpi-tight 32 'mpirun C ./cmcalibrate --mpi foo.cm > foo.out'
 *  
 ************************************************************
 * @LICENSE@
 ************************************************************
 */

#include "esl_config.h"
#include "config.h"	

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_mpi.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_ratematrix.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#define MPI_FINISHED_EXP    -1 /* message to send to workers */
#define MPI_FINISHED_FILTER -2 /* message to send to workers */

#include "funcs.h"		/* external functions                   */
#include "structs.h"

static ESL_OPTIONS options[] = {
  /* name                type           default env   range     toggles      reqs       incomp  help  docgroup*/
  { "-h",                eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "show brief help on version and usage",   1 },
  { "-s",                eslARG_INT,    NULL,   NULL, "n>0",    NULL,        NULL,        NULL, "set random number generator seed to <n>",  1 },
  { "--forecast",        eslARG_INT,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "don't do calibration, forecast running time with <n> processors", 1 },
  { "--devhelp",         eslARG_NONE,   NULL,   NULL, NULL,     NULL,        NULL,        NULL, "show list of undocumented developer options", 1 },
#ifdef HAVE_MPI
  { "--mpi",            eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "run as an MPI parallel program", 1 },  
#endif
  /* options for exp tail fitting */
  { "--exp-cmL-glc",    eslARG_REAL,    "1.5",  NULL, "0.1<=x<=1000.", NULL,NULL,         NULL, "set glocal  CM     Mb random seq length to <x>", 2 },
  { "--exp-cmL-loc",    eslARG_REAL,    "1.5",  NULL, "0.1<=x<=1000.", NULL,NULL,         NULL, "set  local  CM     Mb random seq length to <x>", 2 },
  { "--exp-hmmLn-glc",  eslARG_REAL,    "15.",  NULL, "2.<=x<=1000.",   NULL,NULL,        NULL, "set glocal HMM min Mb random seq length to <x>", 2 },
  { "--exp-hmmLn-loc",  eslARG_REAL,    "15.",  NULL, "2.<=x<=1000.",   NULL,NULL,        NULL, "set  local HMM min Mb random seq length to <x>", 2 },
  { "--exp-hmmLx",      eslARG_REAL,    "1000.",NULL, "10.<=x<=1001.",  NULL,NULL,        NULL, "set        HMM max Mb random seq length to <x>", 2 },
  { "--exp-fract",      eslARG_REAL,    "0.10", NULL, "0.01<=x<=1.0",   NULL,NULL,        NULL, "set min fraction of HMM vs CM DP calcs to <x>", 2 },
  { "--exp-tailn-cglc", eslARG_INT,    "250",  NULL, "n>=100",  NULL,        NULL,"--exp-tailp","fit the top <n> hits/Mb in histogram for  CM local modes", 2 },
  { "--exp-tailn-cloc", eslARG_INT,    "750",  NULL, "n>=100",  NULL,        NULL,"--exp-tailp","fit the top <n> hits/Mb in histogram for  CM glocal modes", 2 },
  { "--exp-tailn-hglc", eslARG_INT,     "25",  NULL, "n>=10",   NULL,        NULL,"--exp-tailp","fit the top <n> hits/Mb in histogram for HMM local modes", 2 },
  { "--exp-tailn-hloc", eslARG_INT,     "75",  NULL, "n>=10",   NULL,        NULL,"--exp-tailp","fit the top <n> hits/Mb in histogram for HMM glocal modes", 2 },
  { "--exp-tailp",      eslARG_REAL,    NULL,   NULL, "0.0<x<0.6",NULL,      NULL,        NULL, "set fraction of histogram tail to fit to exp tail to <x>", 2 },
  { "--exp-tailxn",     eslARG_INT,     "1000", NULL, "n>=50",  NULL,        NULL,"--exp-tailp", "w/--exp-tailp, set max num hits in tail to fit as <n>", 2 },
  { "--exp-beta",       eslARG_REAL,    "1E-15",NULL, "x>0",    NULL,        NULL,"--exp-no-qdb","set tail loss prob for QDB to <x>", 2 },
  { "--exp-no-qdb",     eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "do not use QDBs for calibrating CM search modes", 2 },
  { "--exp-hfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save fitted score histogram(s) to file <f>", 2 },
  { "--exp-sfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save survival plot to file <f>", 2 },
  { "--exp-qqfile",     eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save Q-Q plot for score histogram(s) to file <f>", 2 },
  { "--exp-ffile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,        NULL, "save lambdas for different tail fit probs to file <f>", 2 },
  /* options for HMM filter threshold calculation */
  { "--fil-N",          eslARG_INT,     "10000",NULL, "100<=n<=100000",NULL,  NULL,       NULL, "number of emitted sequences for HMM filter threshold calc",    3 },
  { "--fil-F",          eslARG_REAL,    "0.995",NULL, "0<x<=1", NULL,        NULL,        NULL, "required fraction of seqs that survive HMM filter", 3},
  { "--fil-xhmm",       eslARG_REAL,    "2.0",  NULL, "x>=1.1", NULL,        NULL,        NULL, "set target time for filtered search as <x> times HMM time", 3},
  { "--fil-tau",        eslARG_REAL,    "1e-7", NULL, "0<x<1",  NULL,        NULL,"--fil-nonbanded", "set tail loss prob for HMM banding <x>", 3 },
  { "--fil-gemit",      eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "when calc'ing filter thresholds, always emit globally from CM",  3},
  { "--fil-dfile",      eslARG_OUTFILE, NULL,   NULL, NULL,     NULL,        NULL,"--exp-pfile", "save filter threshold data (HMM and CM scores) to file <s>", 3},
  /* Other options */
  { "--mxsize",         eslARG_REAL,    "2048.0",NULL, "x>0.",  NULL,        NULL,        NULL, "set maximum allowable HMM banded DP matrix size to <x> Mb", 4 },
  /* All options below are developer options, only shown if --devhelp invoked */
  /* Developer option, print extra info */
  { "-v",                eslARG_NONE,   FALSE,  NULL, NULL,     NULL,        NULL, "--forecast", "print arguably interesting info",  101},
#ifdef HAVE_MPI
  /* Developer option, for debugging */
  { "--stall",          eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,        NULL, "arrest after start: for debugging MPI under gdb", 101 },  
#endif
  /* Developer exponential tail options the average user doesn't need to know about */
  { "--exp-random",     eslARG_NONE,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "use GC content of random null background model of CM",  102},
  { "--exp-T",          eslARG_REAL,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "set bit sc cutoff for exp tail fitting to <x> [df: -INFTY]", 102 },
  { "--exp-pfile",      eslARG_INFILE,  NULL,   NULL, NULL,     NULL,  "--exp-gc",        NULL, "read partition info for exp tails from file <f>", 102},
  { "--exp-gc",         eslARG_INFILE,  NULL,   NULL, NULL,     NULL,        NULL,        NULL, "use GC content distribution from file <f>",  102},
  /* Developer filter threshold options the average user doesn't need to know about */
  { "--fil-nonbanded",  eslARG_NONE,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "do not use HMM banded search for filter calculation", 104},
  { "--fil-aln2bands",  eslARG_NONE,    FALSE,  NULL, NULL,     NULL,        NULL,"--fil-nonbanded", "derive HMM bands w/o scanning Forward/Backward", 104 },
  /* Developer options related to experiment local begin/end modes */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,    "--pbegin", "set all local begins as equiprobable", 103 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,    NULL,    "--pend", "set all local end probs to <x>", 103 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,    NULL,        NULL, "set aggregate local begin prob to <x>", 103 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,    NULL,        NULL, "set aggregate local end prob to <x>", 103 },
  { "--no-null3", eslARG_NONE,   FALSE,  NULL, NULL,      NULL,        NULL,        NULL, "turn OFF the NULL3 post hoc additional null model", 103 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char            *cmfile;	      /* name of input CM file  */ 
  ESL_RANDOMNESS  *r;
  ESL_ALPHABET    *abc;
  ESL_STOPWATCH   *w_stage;           /* stopwatch for each exp, filter stage */
  double          *gc_freq;
  double          *pgc_freq;
  CMStats_t      **cmstatsA;          /* the CM stats data structures, 1 for each CM */
  int              ncm;                /* what number CM we're on */
  int              be_verbose;	       /* print extra info */
  int              cmalloc;            /* number of cmstats we have allocated */
  char            *tmpfile;            /* tmp file we're writing to */
  char            *mode;               /* write mode, "w" or "wb"                     */
  int              np;                 /* number of partitions, must be 1 unless --exp-pfile invoked,
					* once set, is never changed */
  int             *pstart;             /* [0..p..np-1], begin points for partitions, end pts are implicit */
  float            avg_hit_len;        /* average CM hit length, calc'ed using QDB calculation algorithm */
  double          *dnull;              /* double version of cm->null, for generating random seqs */

  /* number of sequences and the length of each seq for exp tail fitting, set such that:
   * exp_{cm,hmm}N_{loc,glc} are the number of 10 Kb seqs we'll search for CM/HMM local/glocal
   * exponential tail fitting:
   *
   * exp_cmN_loc  = (esl_opt_GetBoolean(go, "--exp-cmL-loc")  * 1,000,000) / 100,000; 
   * exp_cmN_glc  = (esl_opt_GetBoolean(go, "--exp-cmL-glc")  * 1,000,000) / 100,000; 
   * exp_hmmN_loc = (esl_opt_GetBoolean(go, "--exp-hmmL-loc") * 1,000,000) / 100,000; 
   * exp_hmmN_glc = (esl_opt_GetBoolean(go, "--exp-hmmL-glc") * 1,000,000) / 100,000; 
   *
   * We don't search just 1 long sequence (for ex of 1 Mb) b/c using sequence lengths
   * above 10 Kb for exp tail calibration can yield millions of hits 
   * (for CM scans) before overlaps are removed, which requires a lot
   * of memory. 
   */
  int              exp_cmN_loc;        /* number of 10 Kb seqs for  local CM exp tail fitting */
  int              exp_cmN_glc;        /* number of 10 Kb seqs for glocal CM exp tail fitting */
  int              exp_hmmN_loc;       /* number of 10 Kb seqs for  local HMM exp tail fitting */
  int              exp_hmmN_glc;       /* number of 10 Kb seqs for glocal HMM exp tail fitting */
  int              expL;               /* the size of seq chunks to search, set as 10,000 (10 Kb) */

  /* info for the comlog we'll add to the cmfiles */
  char            *ccom;               /* command line used in this execution of cmcalibrate */
  char            *cdate;              /* date of this execution of cmcalibrate */

  /* the following data is modified for each CM, and in some cases for each exp mode for each CM,
   * it is assumed to be 'current' in many functions.
   */
  float            fil_cm_ncalcs;   /* millions of calcs for full CM scan of 1 residue, with QDBs from beta == cm->beta_qdb */
  float            cp9_ncalcs;      /* millions of calcs for CP9 HMM scan of 1 residue, updated when model is localfied */
  /* mpi */
  int              do_mpi;
  int              my_rank;
  int              nproc;
  int              do_stall;          /* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE          *cmfp;	      /* open input CM file stream       */
  FILE            *exphfp;            /* optional output for exp tail histograms */
  FILE            *expsfp;            /* optional output for exp tail survival plot */
  FILE            *expqfp;            /* optional output for exp tail QQ file */
  FILE            *exptfitfp;         /* optional output for exp tail fit file */
  FILE            *fildfp;            /* optional output for filter scores */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit exp tails for E-values and determine HMM filter thresholds";

static int init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static int   mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int   mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int  process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
				   float **ret_cyk_scA, float **ret_ins_scA, float **ret_fwd_scA, int **ret_partA);
static int  initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);

static int  set_partition_gc_freq(struct cfg_s *cfg, int p);
static int  fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, int exp_mode, double *ret_mu, double *ret_lambda, int *ret_nrandhits, float *ret_tailp);
static int  get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, double *dnull, int L, ESL_DSQ **ret_dsq);
static int  get_cmemit_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int *ret_L, int *ret_p, ESL_DSQ **ret_dsq);
static int  read_partition_file(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  switch_global_to_local(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, char *errbuf);
extern int  update_dp_calcs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm);
static int  get_hmm_filter_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_scA, float *fwd_scA, int *partA, int cm_mode, HMMFilterInfo_t *hfi);
static int  compare_fseq_by_cm_Eval (const void *a_void, const void *b_void);
static int  compare_fseq_by_fwd_Eval(const void *a_void, const void *b_void);
static int  set_dnull(CM_t *cm, char *errbuf, double **ret_dnull);
static int  print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int  print_per_cm_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  print_per_cm_summary        (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, double psec, double asec);
static int  print_exp_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, int expN, int expL, int p, double psec);
static int  print_fil_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, double psec);
static int  print_post_calibration_info (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, double **exp_psecAA, double *fil_psecA, double **exp_asecAA, double *fil_asecA);
static int  estimate_time_for_exp_round (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_res);
static int  estimate_time_for_fil_round (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_seq);
static int  update_hmm_exp_length(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  get_genomic_sequence_from_hmm(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int L, ESL_DSQ **ret_dsq);
/*static int  predict_time_for_exp_stage(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, int cmN, int hmmN, int expL, float *ret_seconds);*/
/*static int  print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);*/

int
main(int argc, char **argv)
{
  int              status;
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory allocation error, stopwatch could not be created.");
  esl_stopwatch_Start(w);
  struct cfg_s     cfg;
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
  if (esl_opt_GetBoolean(go, "--devhelp") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nexponential tail distribution fitting options :");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nHMM filter threshold calculation options :");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\nother options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nundocumented developer options for debugging:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\nundocumented exp tail related developer options:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80);
      puts("\nundocumented filter related developer options:");
      esl_opt_DisplayHelp(stdout, go, 104, 2, 80);
      puts("\nother undocumented developer options:");
      esl_opt_DisplayHelp(stdout, go, 103, 2, 80);
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nexponential tail distribution fitting options :");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nHMM filter threshold calculation options :");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\nother options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  /* Initialize configuration shared across all kinds of masters
   * and workers in this .c file.
   */
  cfg.cmfile  = esl_opt_GetArg(go, 1);
  if (cfg.cmfile == NULL) cm_Fail("Failed to read <cmfile> argument from command line.");
  cfg.r           = NULL; 
  cfg.abc         = NULL; 
  cfg.w_stage     = NULL; 
  cfg.gc_freq     = NULL; 
  cfg.pgc_freq    = NULL; 
  cfg.be_verbose  = FALSE;
  cfg.cmalloc     = 128;
  cfg.tmpfile     = NULL;
  cfg.mode        = NULL;
  cfg.np          = 1;     /* default number of partitions is 1, changed if --exp-pfile */
  cfg.pstart      = NULL;  /* allocated (by default to size 1) in init_master_cfg() */
  cfg.avg_hit_len = 0.;
  cfg.expL        = 10000; /* 10 Kb chunks are searched */

  /* Initial allocations for results per CM;
   * we'll resize these arrays dynamically as we read more CMs.
   */
  cfg.cmalloc  = 128;
  ESL_ALLOC(cfg.cmstatsA, sizeof(CMStats_t *) * cfg.cmalloc);
  cfg.ncm      = 0;

  cfg.ccom      = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */
  cfg.cdate     = NULL;  /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */

  cfg.fil_cm_ncalcs = 0;
  cfg.cp9_ncalcs    = 0;

  cfg.do_mpi   = FALSE;
  cfg.my_rank  = 0;
  cfg.nproc    = 0;
  cfg.do_stall = FALSE;
#ifdef HAVE_MPI
  cfg.do_stall = esl_opt_GetBoolean(go, "--stall");
#endif

  /* calculate sequence lengths and quantities for exp tail fitting, */
  int cmL_total_nt_loc  = (int) (1000000. * esl_opt_GetReal(go, "--exp-cmL-loc"));
  int cmL_total_nt_glc  = (int) (1000000. * esl_opt_GetReal(go, "--exp-cmL-glc"));
  int hmmL_total_nt_loc = (int) (1000000. * esl_opt_GetReal(go, "--exp-hmmLn-loc"));
  int hmmL_total_nt_glc = (int) (1000000. * esl_opt_GetReal(go, "--exp-hmmLn-glc"));

  /* determine the number of 10 Kb chunks (cfg.expL = 10000) to search to reach the totals */
  if(cmL_total_nt_loc  < (cfg.expL + 1)) cm_Fail("with --exp-cmL-loc <x>, <x> must be at least %.3f.", cfg.expL / 1000000.);
  if(cmL_total_nt_glc  < (cfg.expL + 1)) cm_Fail("with --exp-cmL-glc <x>, <x> must be at least %.3f.", cfg.expL / 1000000.);
  if(hmmL_total_nt_loc < (cfg.expL + 1)) cm_Fail("with --exp-cmL-loc <x>, <x> must be at least %.3f.", cfg.expL / 1000000.);
  if(hmmL_total_nt_glc < (cfg.expL + 1)) cm_Fail("with --exp-cmL-glc <x>, <x> must be at least %.3f.", cfg.expL / 1000000.);
  cfg.exp_cmN_loc  = (int) (((float) cmL_total_nt_loc  / (float) cfg.expL) + 0.999999); 
  cfg.exp_cmN_glc  = (int) (((float) cmL_total_nt_glc  / (float) cfg.expL) + 0.999999); 
  cfg.exp_hmmN_loc = (int) (((float) hmmL_total_nt_loc / (float) cfg.expL) + 0.999999); 
  cfg.exp_hmmN_glc = (int) (((float) hmmL_total_nt_glc / (float) cfg.expL) + 0.999999); 

  cfg.cmfp     = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.exphfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.expsfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.expqfp   = NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.exptfitfp= NULL; /* ALWAYS remains NULL for mpi workers */
  cfg.fildfp   = NULL; /* ALWAYS remains NULL for mpi workers */

  ESL_DASSERT1((EXP_CP9_GV == 0));
  ESL_DASSERT1((EXP_CP9_GF == 1));
  ESL_DASSERT1((EXP_CM_GC  == 2));
  ESL_DASSERT1((EXP_CM_GI  == 3));
  ESL_DASSERT1((EXP_CP9_LV == 4));
  ESL_DASSERT1((EXP_CP9_LF == 5));
  ESL_DASSERT1((EXP_CM_LC  == 6));
  ESL_DASSERT1((EXP_CM_LI  == 7));
  ESL_DASSERT1((EXP_NMODES == 8));
  ESL_DASSERT1((FTHR_CM_GC == 0));
  ESL_DASSERT1((FTHR_CM_GI == 1));
  ESL_DASSERT1((FTHR_CM_LC == 2));
  ESL_DASSERT1((FTHR_CM_LI == 3));
  ESL_DASSERT1((FTHR_NMODES== 4));

  /* This is our stall point, if we need to wait until we get a
   * debugger attached to this process for debugging (especially
   * useful for MPI):
   */
  while (cfg.do_stall); 

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      char             errbuf[cmERRBUFSIZE]; /* for error messages in mpi_master() */
      if(! esl_opt_IsDefault(go, "--forecast")) cm_Fail("--forecast is incompatible with --mpi.");
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you run mpirun?)");

      if (cfg.my_rank > 0)  { status = mpi_worker(go, &cfg); }
      else { 
	cm_banner(stdout, argv[0], banner);
	status = mpi_master(go, &cfg, errbuf);
      }
      /* check status, if eslOK, we continue, else we exit. either way we call MPI_Finalize() */
      if(status == eslOK) { 
	esl_stopwatch_Stop(w);
	esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
	MPI_Finalize();
      }
      else { /* status != eslOK, master has error message in errbuf, worker does not */
	MPI_Finalize();
	if(cfg.my_rank == 0) cm_Fail(errbuf); /* master */
	else                 return 0;        /* worker */
      }
    }
  else
#endif /*HAVE_MPI*/
    {
      cm_banner(stdout, argv[0], banner);
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }

  if(esl_opt_IsDefault(go, "--forecast") && cfg.my_rank == 0) { /* master, serial or mpi */
    /* Rewind the CM file for a second pass.
     * Write a temporary CM file with new stats information in it
     */
    int   status;
    int   cmi;
    CM_t *cm;
    FILE *outfp;
    sigset_t blocksigs;  /* list of signals to protect from             */
    char     errbuf[cmERRBUFSIZE];

    CMFileRewind(cfg.cmfp);
    if (esl_FileExists(cfg.tmpfile))                    cm_Fail("Ouch. Temporary file %s appeared during the run.", cfg.tmpfile);
    if ((outfp = fopen(cfg.tmpfile, cfg.mode)) == NULL) cm_Fail("Ouch. Temporary file %s couldn't be opened for writing.", cfg.tmpfile); 
    
    for (cmi = 0; cmi < cfg.ncm; cmi++) {
      if ((status = CMFileRead(cfg.cmfp, errbuf, &(cfg.abc), &cm)) != eslOK) cm_Fail("Ran out of CMs too early in pass 2");
      if (cm == NULL)                                                        cm_Fail("CM file %s was corrupted? Parse failed in pass 2", cfg.cmfile);

      /* update the cm->comlog info */
      if((status = update_comlog(go, errbuf, cfg.ccom, cfg.cdate, cm)) != eslOK) cm_Fail(errbuf);
	
      if(cm->flags & CMH_EXPTAIL_STATS) FreeCMStats(cm->stats); 
      cm->stats = cfg.cmstatsA[cmi];
      cm->flags |= CMH_EXPTAIL_STATS; 
      cm->flags |= CMH_FILTER_STATS; 
      if((status = CMFileWrite(outfp, cm, cfg.cmfp->is_binary, errbuf)) != eslOK) cm_Fail(go->errbuf);
      FreeCM(cm);
    } /* end of from idx = 0 to ncm */
    
    /* Now, carefully remove original file and replace it
     * with the tmpfile. Note the protection from signals;
     * we wouldn't want a user to ctrl-C just as we've deleted
     * their CM file but before the new one is moved.
     */
    CMFileClose(cfg.cmfp);
    if (fclose(outfp)   != 0)                            cm_Fail("system error during rewrite of CM file");
    if (sigemptyset(&blocksigs) != 0)                    cm_Fail("system error during rewrite of CM file.");;
    if (sigaddset(&blocksigs, SIGINT) != 0)              cm_Fail("system error during rewrite of CM file.");;
    if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   cm_Fail("system error during rewrite of CM file.");;
    if (remove(cfg.cmfile) != 0)                         cm_Fail("system error during rewrite of CM file.");;
    if (rename(cfg.tmpfile, cfg.cmfile) != 0)            cm_Fail("system error during rewrite of CM file.");;
    if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) cm_Fail("system error during rewrite of CM file.");;
    free(cfg.tmpfile);
    
    /* master specific cleaning */
    if (cfg.exphfp   != NULL) { 
      fclose(cfg.exphfp);
      printf("# Histogram of high scoring hits in random seqs saved to file %s.\n", esl_opt_GetString(go, "--exp-hfile"));
    }
    if (cfg.expsfp   != NULL) { 
      fclose(cfg.expsfp);
      printf("# Survival plot for exponential tails saved to file %s.\n", esl_opt_GetString(go, "--exp-sfile"));
    }
    if (cfg.expqfp   != NULL) { 
      fclose(cfg.expqfp);
      printf("# Exponential tail QQ plots saved to file %s.\n", esl_opt_GetString(go, "--exp-qqfile"));
    }
    if (cfg.exptfitfp   != NULL) { 
      fclose(cfg.exptfitfp);
      printf("# Exponential tail fit points saved to file %s.\n", esl_opt_GetString(go, "--exp-ffile"));
    }
    if (cfg.fildfp   != NULL) { 
      fclose(cfg.fildfp);
      printf("# Filter histograms saved to file %s.\n", esl_opt_GetString(go, "--fil-dfile"));
    }
    if (cfg.ccom     != NULL) free(cfg.ccom);
    if (cfg.cdate    != NULL) free(cfg.cdate);
    if (cfg.cmstatsA != NULL) free(cfg.cmstatsA);
  }

  /* clean up */
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.pstart    != NULL) free(cfg.pstart);
  if (cfg.w_stage   != NULL) esl_stopwatch_Destroy(cfg.w_stage);
  if (cfg.r         != NULL) esl_randomness_Destroy(cfg.r);
  if (cfg.my_rank == 0) { 
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
  }
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;

 ERROR: 
  cm_Fail("Memory allocation error.");
  return 1; /* NEVERREACHED */
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Allocates/sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->exphfp      - optional output file
 *    cfg->expsfp      - optional output file
 *    cfg->expqfp      - optional output file
 *    cfg->fildfp      - optional output file
 *    cfg->gc_freq     - observed GC freqs (if --exp-gc invoked)
 *    cfg->cmstatsA    - the stats, allocated only
 *    cfg->np          - number of partitions, never changes once set 
 *    cfg->pstart      - array of partition starts 
 *    cfg->r           - source of randomness
 *    cfg->tmpfile     - temp file for rewriting cm file
 *    cfg->be_verbose  - print extra info? 
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errbuf) by the ESL_FAIL mech.
 */
static int
init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* optionally, open exp tail histogram file */
  if (esl_opt_GetString(go, "--exp-hfile") != NULL) 
    {
      if ((cfg->exphfp = fopen(esl_opt_GetString(go, "--exp-hfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail histogram save file %s for writing\n", esl_opt_GetString(go, "--exp-hfile"));
    }

  /* optionally, open survival plot */
  if (esl_opt_GetString(go, "--exp-sfile") != NULL) 
    {
      if ((cfg->expsfp = fopen(esl_opt_GetString(go, "--exp-sfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open survival plot save file %s for writing\n", esl_opt_GetString(go, "--exp-sfile"));
    }

  /* optionally, open exp tail QQ plot file */
  if (esl_opt_GetString(go, "--exp-qqfile") != NULL) 
    {
      if ((cfg->expqfp = fopen(esl_opt_GetString(go, "--exp-qqfile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail QQ plot save file %s for writing\n", esl_opt_GetString(go, "--exp-qqfile"));
    }

  /* optionally, open exp tail tail fit prob file */
  if (esl_opt_GetString(go, "--exp-ffile") != NULL) 
    {
      if ((cfg->exptfitfp = fopen(esl_opt_GetString(go, "--exp-ffile"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail save file %s for writing\n", esl_opt_GetString(go, "--exp-ffile"));
    }

  /* optionally, open filter threshold data file */
  if (esl_opt_GetString(go, "--fil-dfile") != NULL) {
    if ((cfg->fildfp = fopen(esl_opt_GetString(go, "--fil-dfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --fil-dfile output file %s\n", esl_opt_GetString(go, "--fil-dfile"));
    }

  /* optionally, get distribution of GC content from an input database (default is use cm->null for GC distro) */
  if(esl_opt_GetString(go, "--exp-gc") != NULL) {
    ESL_ALPHABET *tmp_abc = NULL;
    tmp_abc = esl_alphabet_Create(eslRNA);
    ESL_SQFILE      *dbfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--exp-gc"), eslSQFILE_UNKNOWN, NULL, &dbfp);
    if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "No such file: %s.", esl_opt_GetString(go, "--exp-gc")); 
    else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "file: %s format unrecognized.", esl_opt_GetString(go, "--exp-gc")); 
    else if (status != eslOK)      ESL_FAIL(status, errbuf, "Failed to open sequence database file %s, code %d.", esl_opt_GetString(go, "--exp-gc"), status); 
    if((status = GetDBInfo(tmp_abc, dbfp, errbuf, NULL, NULL, &(cfg->gc_freq))) != eslOK) return status; 
    esl_vec_DNorm(cfg->gc_freq, GC_SEGMENTS);
    esl_alphabet_Destroy(tmp_abc);
    esl_sqfile_Close(dbfp); 
   /* allocate pgc_freq, the gc freqs per partition, used to sample seqs for different partitions */
    ESL_ALLOC(cfg->pgc_freq, sizeof(double) * GC_SEGMENTS);
  }

  /* set up the partition data that's used for all CMs */
  if(esl_opt_IsDefault(go, "--exp-pfile")) { /* by default we have 1 partition 0..100 */
    ESL_ALLOC(cfg->pstart, sizeof(int) * 1);
    cfg->np        = 1;
    cfg->pstart[0] = 0;
  }
  else { /* setup cfg->np and cfg->pstart in read_partition_file() */
    if((status = read_partition_file(go, cfg, errbuf)) != eslOK) return status;
  }

  if (esl_opt_GetString(go, "--fil-dfile") != NULL) { if(cfg->np != 1) ESL_FAIL(eslEINVAL, errbuf, "--fil-dfile only works with a single partition\n"); }

  cfg->be_verbose = FALSE;
  if (esl_opt_GetBoolean(go, "-v")) cfg->be_verbose = TRUE;        

  /* seed master's RNG */
  if (! esl_opt_IsDefault(go, "-s")) 
    cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
  else cfg->r = esl_randomness_CreateTimeseeded();

  /* create the stopwatch */
  cfg->w_stage = esl_stopwatch_Create();

  /* From HMMER 2.4X hmmcalibrate.c:
   * Generate calibrated CM(s) in a tmp file in the current
   * directory. When we're finished, we delete the original
   * CM file and rename() this one. That way, the worst
   * effect of a catastrophic failure should be that we
   * leave a tmp file lying around, but the original CM
   * file remains uncorrupted. tmpnam() doesn't work portably here,
   * because it'll put the file in /tmp and we won't
   * necessarily be able to rename() it from there.
   */
  ESL_ALLOC(cfg->tmpfile, (sizeof(char) * (strlen(cfg->cmfile) + 5)));
  strcpy(cfg->tmpfile, cfg->cmfile);
  strcat(cfg->tmpfile, ".xxx");	/* could be more inventive here... */
  if (esl_FileExists(cfg->tmpfile))
    ESL_FAIL(eslFAIL, errbuf, "temporary file %s already exists; please delete it first", cfg->tmpfile);
  if (cfg->cmfp->is_binary) cfg->mode = "wb";
  else                      cfg->mode = "w"; 

  if(cfg->r        == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create master RNG.");
  if(cfg->w_stage  == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create stopwatch.");

  /* fill cfg->ccom, and cfg->cdate */
  if((status = get_cmcalibrate_comlog_info(go, cfg, errbuf)) != eslOK) return status;

  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "init_master_cfg(), memory allocation error."); 
  return eslEMEM; /* NEVER REACHED */
}

/* serial_master()
 * The serial version of cmcalibrate.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;                /* Easel status */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  int      p;                     /* partition index */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  long     seed;                  /* RNG seed */
  int      exp_mode;              /* ctr over exp tail modes */
  double   cm_psec;               /* predicted number of seconds for calibrating current CM */
  double   cm_asec;               /* predicted number of seconds for calibrating current CM */
  double   total_psec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   total_asec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double **exp_asecAA;            /* stores actual timings for each exp tail fit stage, for each CM, each partition */
  double  *fil_asecA;             /* stores actual timings for each filter stage for each CM */
  double **exp_psecAA;            /* stores predicted timings for each exp tail fit stage, for each CM, each partition */
  double  *fil_psecA;             /* stores predicted timings for each filter stage for each CM */
  double   psec;                  /* predicted seconds */

  /* exptail related vars */
  int               expN;                                        /* number of length <expL> sequences to search for exp tail fitting of current exp mode */
  int               exp_cm_cyk_mode;                             /* CYK    exp mode CM is in EXP_CM_LC or EXP_CM_GC */
  int               exp_cm_ins_mode;                             /* Inside exp mode CM is in EXP_CM_LI or EXP_CM_GI */
  int               exp_scN = 0;                                 /* number of hits reported thus far, for all seqs */
  float            *exp_scA = NULL;                              /* [0..exp_scN-1] hit scores for all seqs */
  ESL_DSQ          *dsq = NULL;                                  /* digitized sequence to search */
  search_results_t *results;                                     /* results (hits) from current sequence */
  double           *dnull = NULL;                                /* double version of cm->null, for generating random seqs */
  int               i;                                           /* counter over sequences */
  int               h;                                           /* counter over hits */
  double            tmp_mu, tmp_lambda;                          /* temporary mu and lambda used for setting exp tails */
  int               tmp_nrandhits;                               /* temporary number of rand hits found */
  float             tmp_tailp;                                   /* temporary tail mass probability fit to an exponential */

  /* filter threshold related vars */
  int      filN = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences to search for filter threshold calculation */
  int      fthr_mode = 0;         /* CM mode for filter threshold calculation, FTHR_CM_GC, FTHR_CM_GI, FTHR_CM_LC, FTHR_CM_LI */
  float   *fil_cyk_scA = NULL;    /* [0..filN-1] best cm cyk score for each emitted seq */
  float   *fil_ins_scA = NULL;    /* [0..filN-1] best cm insidei score for each emitted seq */
  float   *fil_fwd_scA = NULL;    /* [0..filN-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA   = NULL;    /* [0..filN-1] partition of CM emitted seq */
  int      fil_cm_cyk_mode;       /* CYK    fthr mode CM is in FTHR_CM_LC or FTHR_CM_GC */
  int      fil_cm_ins_mode;       /* Inside fthr mode CM is in FTHR_CM_LI or FTHR_CM_GI */
  

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = print_run_info (go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  
  while ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK)
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->cmstatsA, tmp, sizeof(CMStats_t *) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;
      
      if((status = initialize_cm(go, cfg, errbuf, cm))               != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf, cm))          != eslOK) cm_Fail(errbuf);
      if((status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len))) != eslOK) cm_Fail(errbuf);
      if((status = print_per_cm_column_headings(go, cfg, errbuf, cm))!= eslOK) cm_Fail(errbuf);
      if((status = update_hmm_exp_length(go, cfg, errbuf, cm))       != eslOK) cm_Fail(errbuf);
      if(esl_opt_IsDefault(go, "--exp-gc")) { /* only setup dnull if --exp-gc NOT enabled */
	if((status = set_dnull(cm, errbuf, &dnull))                    != eslOK) cm_Fail(errbuf); 
      }
      
      /* allocate the exp_{a,p}secAA and fil_{a,p}secA arrays that hold {actual,predicted} times */
      ESL_ALLOC(exp_asecAA, sizeof(double *) * (EXP_NMODES));
      ESL_ALLOC(fil_asecA,  sizeof(double) *   (EXP_NMODES));
      ESL_ALLOC(exp_psecAA, sizeof(double *) * (EXP_NMODES));
      ESL_ALLOC(fil_psecA,  sizeof(double) *   (EXP_NMODES));
      esl_vec_DSet(fil_asecA, EXP_NMODES, 0.);
      esl_vec_DSet(fil_psecA, EXP_NMODES, 0.);
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) { 
	ESL_ALLOC(exp_asecAA[exp_mode], sizeof(double) * (cfg->np));
	ESL_ALLOC(exp_psecAA[exp_mode], sizeof(double) * (cfg->np));
	esl_vec_DSet(exp_asecAA[exp_mode], cfg->np, 0.);
	esl_vec_DSet(exp_psecAA[exp_mode], cfg->np, 0.);
      }
      cm_psec = cm_asec = 0.;

      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
	if(ExpModeIsLocal(exp_mode)) { expN = ExpModeIsForCM(exp_mode) ? cfg->exp_cmN_loc : cfg->exp_hmmN_loc; }
	else                         { expN = ExpModeIsForCM(exp_mode) ? cfg->exp_cmN_glc : cfg->exp_hmmN_glc; }
	
	/* do we need to switch from glocal configuration to local? */
	if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf))      != eslOK) cm_Fail(errbuf);
	  if((status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len))) != eslOK) cm_Fail(errbuf);
	}
	/* update search info for round 0 (final round) for exp tail mode */
	UpdateSearchInfoForExpMode(cm, 0, exp_mode);

	/* We want to use the same seqs for exp tail fittings of all CM modes and HMM modes, 
	 * so we free RNG, then create a new one and reseed it with the initial seed,
	 * The following pairs of modes will have identical sequences used for each member of the pair:
	 * 1. EXP_CP9_GV and EXP_CP9_GF
	 * 2. EXP_CM_GC  and EXP_CM_GI
	 * 3. EXP_CP9_LV and EXP_CP9_LF
	 * 4. EXP_CM_LC  and EXP_CM_LI
	 */
	seed = esl_randomness_GetSeed(cfg->r);
	esl_randomness_Destroy(cfg->r);
	cfg->r = esl_randomness_Create(seed);
	
	/************************************/
	/* exponential tail fitting section */
	/************************************/
	/* calculate exp tails for this exp mode */
	/* determine length of seqs to search for exp tail fitting */
	ESL_DASSERT1((cfg->np == cfg->cmstatsA[cmi]->np));
	for (p = 0; p < cfg->np; p++) {
	  if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	  /* estimate time for this round */
	  if(p == 0) { 
	    if((status = estimate_time_for_exp_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf); 
	    psec *= expN * cfg->expL; /* psec was per residue */
	    /* with --forecast, take into account parallelization */
	    if((! esl_opt_IsDefault(go, "--forecast")) && (esl_opt_GetInteger(go, "--forecast") > 1)) psec /= (esl_opt_GetInteger(go, "--forecast") - 1);
	  }
	  else psec = exp_psecAA[exp_mode][0];
	  exp_psecAA[exp_mode][p] = psec;
	  cm_psec    += psec;
	  total_psec += psec;
	  print_exp_line(go, cfg, errbuf, exp_mode, expN, cfg->expL, p, psec);
	  if(! esl_opt_IsDefault(go, "--forecast")) continue; /* special mode, we don't do the calibration, just print the predicting timings */

	  esl_stopwatch_Start(cfg->w_stage);
	  fflush(stdout);

	  ESL_DPRINTF1(("\n\ncalling ProcessSearchWorkunit to fit exp tail for p: %d EXP mode: %d\n", p, exp_mode));
	  
	  exp_scN  = 0;
	  for(i = 0; i < expN; i++) { 
	    /* do the work, fit the histogram, update exp tail info in cmstats */
	    /* generate sequence, either randomly from background null or from hard-wired 5 state HMM that emits genome like sequence */
	    if(esl_opt_GetBoolean(go, "--exp-random")) { 
	      if((status = get_random_dsq(cfg, errbuf, cm, dnull, cfg->expL, &dsq)) != eslOK) cm_Fail(errbuf); 
	    }
	    else { 
	      if((status = get_genomic_sequence_from_hmm(cfg, errbuf, cm, cfg->expL, &dsq)) != eslOK) cm_Fail(errbuf); 
	    }
	    if((status = ProcessSearchWorkunit (cm,  errbuf, dsq, cfg->expL, &results, esl_opt_GetReal(go, "--mxsize"), cfg->my_rank, NULL, NULL)) != eslOK) cm_Fail(errbuf);
	    RemoveOverlappingHits(results, 1, cfg->expL);
	    
	    if(results->num_results > 0) { 
	      if(i == 0) ESL_ALLOC (exp_scA, sizeof(float) * (exp_scN + results->num_results));
	      else       ESL_RALLOC(exp_scA, tmp, sizeof(float) * (exp_scN + results->num_results));
	      for(h = 0; h < results->num_results; h++) exp_scA[(exp_scN+h)] = results->data[h].score;
	      exp_scN += results->num_results;
	    }
	    FreeResults(results);
	    free(dsq);
	  }
	  if(cfg->exptfitfp != NULL) { 
	    fprintf(cfg->exptfitfp, "# CM: %s\n", cm->name);
	    fprintf(cfg->exptfitfp, "# mode: %12s\n", DescribeExpMode(exp_mode));
	  }
	  if((status = fit_histogram(go, cfg, errbuf, exp_scA, exp_scN, exp_mode, &tmp_mu, &tmp_lambda, &tmp_nrandhits, &tmp_tailp)) != eslOK) cm_Fail(errbuf);
	  SetExpInfo(cfg->cmstatsA[cmi]->expAA[exp_mode][p], tmp_lambda, tmp_mu, (long) (cfg->expL * expN), tmp_nrandhits, tmp_tailp);
	  
	  esl_stopwatch_Stop(cfg->w_stage);
	  exp_asecAA[exp_mode][p] = cfg->w_stage->elapsed;
	  cm_asec += cfg->w_stage->elapsed;
	  total_asec += cfg->w_stage->elapsed;
	  FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
	  printf("  %10s\n", time_buf);
	  free(exp_scA);
	} /* end of for loop over partitions */
      
	/****************************/
	/* filter threshold section */
	/****************************/
	if(exp_mode == EXP_CM_GI || exp_mode == EXP_CM_LI) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  /* estimate time for this round */
	  if((status = estimate_time_for_fil_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf);
	  psec *= filN;
	  /* with --forecast, take into account parallelization */
	  if((! esl_opt_IsDefault(go, "--forecast")) && (esl_opt_GetInteger(go, "--forecast") > 1)) psec /= (esl_opt_GetInteger(go, "--forecast") - 1);
	  fil_psecA[exp_mode] = psec;
	  cm_psec    += psec;
	  total_psec += psec;
	  print_fil_line(go, cfg, errbuf, exp_mode, psec);
	  if(! esl_opt_IsDefault(go, "--forecast")) continue; /* special mode, we don't do the calibration, just print the predicting timings */

	  esl_stopwatch_Start(cfg->w_stage);
	  fthr_mode = ExpModeToFthrMode(exp_mode);
	  /* search emitted sequences to get filter thresholds for HMM and each candidate sub CM root state */
	  ESL_DPRINTF1(("\n\ncalling process_filter_workunit to get HMM filter thresholds for p: %d mode: %d\n", p, exp_mode));
	  
	  if((status = process_filter_workunit (go, cfg, errbuf, cm, filN, &fil_cyk_scA, &fil_ins_scA, &fil_fwd_scA, &fil_partA)) != eslOK) cm_Fail(errbuf);
	  
	  exp_cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? EXP_CM_LC  : EXP_CM_GC;
	  exp_cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? EXP_CM_LI  : EXP_CM_GI;
	  fil_cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? FTHR_CM_LC : FTHR_CM_GC;
	  fil_cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? FTHR_CM_LI : FTHR_CM_GI;
	  /* set cutoffs for forward HMM filters, first for CYK, then for Inside */
	  if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_cyk_scA, fil_fwd_scA, fil_partA, exp_cm_cyk_mode, cfg->cmstatsA[cmi]->hfiA[fil_cm_cyk_mode])) != eslOK) cm_Fail(errbuf);
	  if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_ins_scA, fil_fwd_scA, fil_partA, exp_cm_ins_mode, cfg->cmstatsA[cmi]->hfiA[fil_cm_ins_mode])) != eslOK) cm_Fail(errbuf);
	  free(fil_cyk_scA);
	  free(fil_ins_scA);
	  free(fil_fwd_scA);
	  free(fil_partA);
	  
	  esl_stopwatch_Stop(cfg->w_stage);
	  fil_asecA[exp_mode] = cfg->w_stage->elapsed;
	  cm_asec += cfg->w_stage->elapsed;
	  total_asec += cfg->w_stage->elapsed;
	  FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
	  printf("  %10s\n", time_buf);
	}
      } /* end of for(exp_mode = 0; exp_mode < NCMMODES; exp_mode++) */
      if(cfg->be_verbose) if((status = debug_print_cmstats(cm, errbuf, cfg->cmstatsA[cmi], TRUE)) != eslOK) cm_Fail(errbuf);
      print_per_cm_summary(go, cfg, errbuf, cm, cm_psec, cm_asec);
      if(esl_opt_IsDefault(go, "--forecast")) { if((status = print_post_calibration_info(go, cfg, errbuf, stdout, cm, exp_psecAA, fil_psecA, exp_asecAA, fil_asecA)) != eslOK) cm_Fail(errbuf); }
      free(dnull);
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) { 
	free(exp_asecAA[exp_mode]); 
	free(exp_psecAA[exp_mode]); 
      }
      free(exp_asecAA);
      free(exp_psecAA);
      free(fil_asecA); 
      free(fil_psecA); 
      FreeCM(cm);
      printf("//\n");
      fflush(stdout);
    }
  if(status != eslEOF) cm_Fail(errbuf);
  
  if(cfg->ncm > 1 && (! esl_opt_IsDefault(go, "--forecast"))) { 
    fprintf(stdout, "#\n");
    FormatTimeString(time_buf, total_psec, FALSE);
    fprintf(stdout, "# total predicted time for all %d CMs: %s\n", cfg->ncm, time_buf);
    fprintf(stdout, "#\n");
  }
  return;
      
 ERROR:
  cm_Fail("Memory allocation error.");
  return; /*NEVERREACHED*/
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmcalibrate
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * A master returns eslOK if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * If a recoverable error occurs, errbuf is filled with an error message
 * from the master or a worker, and it's sent back while returning a
 * non-eslOK error code.
 * 
 * Recoverable errors include most worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers. The 
 * 
 * Some worker side errors (such as ESL_ALLOCs) are likely to be 
 * unrecoverable and will almost certainly cause MPI to crash
 * uncleanly, they're only here because I couldn't find a way around
 * them without massive reimplementation. Hopefully they rarely occur.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * cm_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int 
mpi_master(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int      xstatus       = eslOK;	/* changes from OK on recoverable error */
  int      status;                /* Easel status */
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  int      p;                     /* partition index */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  int      n, i;                  /* counters */
  int      have_work     = TRUE;  /* TRUE while work remains  */
  int      nproc_working = 0;	  /* number of worker processes working, up to nproc-1 */
  int      wi;          	  /* rank of next worker to get a job to work on */
  int      wi_error = 0;          /* worker index that sent back an error message, if an error occurs */
  char    *buf           = NULL;  /* input/output buffer, for packed MPI messages */
  int      bn            = 0;     /* size of buf */
  int      pos = 1;               /* posn in buf */
  int      nseq_sent        = 0;  /* number of seqs we've told workers to work on */
  int      nseq_this_worker = 0;  /* number of seqs to tell current worker to work on */
  int      nseq_just_recv   = 0;  /* number of seqs we just received scores for from a worker */
  int      nseq_recv        = 0;  /* number of seqs we've received thus far this round from workers */
  MPI_Status mpistatus;           /* MPI status... */
  int      msg;                   /* holds integer telling workers we've finished current stage */
  int      exp_mode;              /* ctr over exp tail modes */
  int      h;                     /* ctr over hits */
  long     seed;                  /* for seeding the master's RNG */
  /* variables for predicted and actual timings */
  double   psec;                  /* predicted number of seconds */
  double   cm_psec;               /* predicted number of seconds for calibrating current CM */
  double   cm_asec;               /* predicted number of seconds for calibrating current CM */
  double   total_asec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   total_psec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double **exp_asecAA;            /* stores actual timings for each exp tail fit stage, for each CM, each partition */
  double  *fil_asecA;             /* stores actual timings for each filter stage for each CM */
  double **exp_psecAA;            /* stores predicted timings for each exp tail fit stage, for each CM, each partition */
  double  *fil_psecA;             /* stores predicted timings for each filter stage for each CM */
  /* exponential tail related vars */
  int      expN;                                        /* number of length <cfg->expL> sequences to search for exp tail fitting of current exp mode */
  int      exp_scN = 0;                                 /* number of hits reported thus far, for all seqs */
  float   *exp_scA = NULL;                              /* [0..exp_scN-1] hit scores for all seqs */
  int      exp_cm_cyk_mode;                             /* CYK    exp mode CM is in EXP_CM_LC or EXP_CM_GC */
  int      exp_cm_ins_mode;                             /* Inside exp mode CM is in EXP_CM_LI or EXP_CM_GI */
  int      si;                                          /* sequence index, the index of the last sequence generated */
  int      si_recv;                                     /* sequence index of the sequence we've just received results for from a worker */
  int      seqpos = 1;                                  /* sequence position in the current sequence */
  int      len;                                         /* length of a sequence chunk */
  int      chunksize;                                   /* size of chunks for each worker */
  double   tmp_mu, tmp_lambda;                          /* temporary mu and lambda used for setting exp tails */
  int      tmp_nrandhits;                               /* temporary number of rand hits found */
  float    tmp_tailp;                                   /* temporary tail mass probability fit to an exponential */
  double  *dnull = NULL;                                /* double version of cm->null, for generating random seqs */
  int      need_seq;                                    /* TRUE if we are ready to generate a new seq */
  int      z;                                           /* counter */
  search_results_t *worker_results;                     /* results for seq si_recv we've just rec'd from a worker, we copy it to results_slist[si_recv] */
  /* *_slist variables: lists of data that are specific to each sequence 0..si..expN-1 */
  ESL_DSQ          **dsq_slist = NULL;                  /* [0..si..expN-1], the digitized sequences to search, when finished, they're freed */
  search_results_t **results_slist = NULL;              /* [0..si..expN-1], the compiled results from searching each seq si, when finished, copied to exp_scA and freed */
  int               *chunks_slist = NULL;               /* [0..si..expN-1], number of chunks of seq si currently being searched by workers */
  int               *sent_slist = NULL;                 /* [0..si..expN-1], TRUE if all chunks of seq si have been sent to workers, FALSE if not */
  /* *_wlist variables: lists of data that are specific to each worker 1..wi..nproc-1 */
  int *si_wlist = NULL;                                 /* [0..wi..nproc-1], the sequence index worker wi is working on */
  int *seqpos_wlist= NULL;                              /* [0..wi..nproc-1] the first position of the sequence that worker wi is searching */
  int *len_wlist = NULL;                                /* [0..wi..nproc-1] length of chunk worker wi is searching */

  /* filter threshold related vars */
  int      filN = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences to search for filter threshold calculation */
  int      fthr_mode = 0;         /* CM mode for filter threshold calculation, FTHR_CM_GC, FTHR_CM_GI, FTHR_CM_LC, FTHR_CM_LI */
  int      fil_cm_cyk_mode;       /* CYK    fthr mode CM is in FTHR_CM_LC or FTHR_CM_GC */
  int      fil_cm_ins_mode;       /* Inside fthr mode CM is in FTHR_CM_LI or FTHR_CM_GI */
  int      fil_nseq_per_worker  = (filN / (cfg->nproc-1)); /* when calcing filters, number of seqs to tell each worker to work on */
  long    *seed_wlist = NULL;      /* [0..wi..nproc-1] seeds for worker's RNGs, we send these to workers */
  /* full arrays of CYK, Inside, Fwd scores, [0..filN-1] */
  float   *fil_cyk_scA = NULL;    /* [0..filN-1] best cm cyk score for each emitted seq */
  float   *fil_ins_scA = NULL;    /* [0..filN-1] best cm insidei score for each emitted seq */
  float   *fil_fwd_scA = NULL;    /* [0..filN-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA   = NULL;    /* [0..filN-1] partition of CM emitted seq */
  /* worker's arrays of CYK, Inside, Fwd scores, [0..nseq_per_worker-1], rec'd from workers, copied to full arrays (ex: fil_cyk_scA) */
  float   *wkr_fil_cyk_scA = NULL;/* rec'd from worker: best cm cyk score for each emitted seq */
  float   *wkr_fil_ins_scA = NULL;/* rec'd from worker: best cm insidei score for each emitted seq */
  float   *wkr_fil_fwd_scA = NULL;/* rec'd from worker: best cp9 Forward score for each emitted seq */
  int     *wkr_fil_partA = NULL;  /* rec'd from worker: partition for seq i */
  
  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { if ((status = print_run_info (go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL)         { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((si_wlist       = malloc(sizeof(int)  * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqpos_wlist   = malloc(sizeof(int)  * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((len_wlist      = malloc(sizeof(int)  * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seed_wlist     = malloc(sizeof(long) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* errbuf was filled above */
  ESL_DPRINTF1(("MPI master is initialized\n"));

  for (wi = 0; wi < cfg->nproc; wi++) {
    si_wlist[wi] = seqpos_wlist[wi] = len_wlist[wi] = -1;
    seed_wlist[wi] = esl_rnd_Roll(cfg->r, 1000000000); /* not sure what to use as max for seed */
    ESL_DPRINTF1(("wi %d seed: %ld\n", wi, seed_wlist[wi]));
  }
  
  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  for (wi = 1; wi < cfg->nproc; wi++)
    MPI_Send(&(seed_wlist[wi]), 1, MPI_LONG, wi, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) cm_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));
  
  /* 3 special (annoying) case:
   * case 1: if we've used the --exp-gc option, we read in a seq file 
   * to fill cfg->gc_freq, and we need to broadcast that info to workers
   *
   * case 2: if we are calculating stats for more than 1 partition, 
   * (--exp-pfile invoked), we need to broadcast that information to 
   * the workers. 
   */
  if(! (esl_opt_IsDefault(go, "--exp-gc"))) { /* receive gc_freq info from master */
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  if(! (esl_opt_IsDefault(go, "--exp-pfile"))) { /* broadcast partition info to workers */
    MPI_Bcast(&(cfg->np),  1,       MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pstart != NULL));
    MPI_Bcast(cfg->pstart, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
  }

  /* Main loop: combining load workers, send/receive, clear workers loops;
   * also, catch error states and die later, after clean shutdown of workers.
   * 
   * When a recoverable error occurs, have_work = FALSE, xstatus !=
   * eslOK, and errbuf is set to an informative message. No more
   * errbuf's can be received after the first one. We wait for all the
   * workers to clear their work units, then send them shutdown signals,
   * then finally print our errbuf and exit.
   * 
   * Unrecoverable errors just crash us out with cm_Fail().
   */
  
  while ((xstatus == eslOK) && ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK)) 
    {
      cfg->ncm++;  
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->cmstatsA, tmp, sizeof(CMStats_t *) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;

      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm     (go, cfg, errbuf, cm))          != eslOK) cm_Fail(errbuf);
      if((status = initialize_cmstats(go, cfg, errbuf, cm))          != eslOK) cm_Fail(errbuf);
      if((status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len))) != eslOK) cm_Fail(errbuf);
      if((status = print_per_cm_column_headings(go, cfg, errbuf, cm))!= eslOK) cm_Fail(errbuf);
      if((status = update_hmm_exp_length(go, cfg, errbuf, cm))       != eslOK) cm_Fail(errbuf);
      if(esl_opt_IsDefault(go, "--exp-gc")) { /* only setup dnull if --exp-gc NOT enabled */
	if((status = set_dnull(cm, errbuf, &dnull))                  != eslOK) cm_Fail(errbuf); 
      }
      
      /* allocate the exp_{a,p}secAA and fil_{a,p}secA arrays that hold {actual,predicted} times */
      ESL_ALLOC(exp_asecAA, sizeof(double *) * (EXP_NMODES));
      ESL_ALLOC(fil_asecA,  sizeof(double) *   (EXP_NMODES));
      ESL_ALLOC(exp_psecAA, sizeof(double *) * (EXP_NMODES));
      ESL_ALLOC(fil_psecA,  sizeof(double) *   (EXP_NMODES));
      esl_vec_DSet(fil_asecA, EXP_NMODES, 0.);
      esl_vec_DSet(fil_psecA, EXP_NMODES, 0.);
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) { 
	ESL_ALLOC(exp_asecAA[exp_mode], sizeof(double) * (cfg->np));
	ESL_ALLOC(exp_psecAA[exp_mode], sizeof(double) * (cfg->np));
	esl_vec_DSet(exp_asecAA[exp_mode], cfg->np, 0.);
	esl_vec_DSet(exp_psecAA[exp_mode], cfg->np, 0.);
      }
      cm_psec = cm_asec = 0.;

      ESL_ALLOC(fil_cyk_scA, sizeof(float) * filN);
      ESL_ALLOC(fil_ins_scA, sizeof(float) * filN);
      ESL_ALLOC(fil_fwd_scA, sizeof(float) * filN);
      ESL_ALLOC(fil_partA,   sizeof(int) *   filN);
  
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
	if(ExpModeIsLocal(exp_mode)) { expN = ExpModeIsForCM(exp_mode) ? cfg->exp_cmN_loc : cfg->exp_hmmN_loc; }
	else                         { expN = ExpModeIsForCM(exp_mode) ? cfg->exp_cmN_glc : cfg->exp_hmmN_glc; }
    
	/* allocate and initialize sequence lists */
	ESL_ALLOC(dsq_slist,     sizeof(ESL_DSQ *) * expN);
	ESL_ALLOC(results_slist, sizeof(search_results_t *) * expN);
	ESL_ALLOC(chunks_slist,  sizeof(int) * expN);
	ESL_ALLOC(sent_slist,    sizeof(int) * expN);
	for(z = 0; z < expN; z++) { 
	  dsq_slist[z]    = NULL; 
	  results_slist[z]= NULL;
	  chunks_slist[z] = 0;
	  sent_slist[z]   = FALSE;
	}
    
	/* do we need to switch from glocal configuration to local? */
	if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf))      != eslOK) cm_Fail(errbuf);
	  if((status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len))) != eslOK) cm_Fail(errbuf);
	}
	chunksize = DetermineSeqChunksize(cfg->nproc, cfg->expL, cm->W);
    
	/* update search info for round 0 (final round) for exp tail mode */
	UpdateSearchInfoForExpMode(cm, 0, exp_mode);
    
	/* We want to use the same seqs for exp tail fittings of all CM modes and HMM modes, 
	 * so we free RNG, then create a new one and reseed it with the initial seed,
	 * The following pairs of modes will have identical sequences used for each member of the pair:
	 * 1. EXP_CP9_GV and EXP_CP9_GF
	 * 2. EXP_CM_GC  and EXP_CM_GI
	 * 3. EXP_CP9_LV and EXP_CP9_LF
	 * 4. EXP_CM_LC  and EXP_CM_LI
	 * Also the first min(--exp-cmN-{loc,glc} <n>, --exp-hmmN-{loc-glc} <n>) sequences between 1 and 2, and between 3 and 4,
	 * will also be identical.
	 */
	seed = esl_randomness_GetSeed(cfg->r);
	esl_randomness_Destroy(cfg->r);
	cfg->r = esl_randomness_Create(seed);
    
	/************************************/
	/* exponential tail fitting section */
	/************************************/
	/* fit exponential tails for this exp mode */
	for(p = 0; p < cfg->np; p++) 
	  {
	    if(cfg->gc_freq != NULL) set_partition_gc_freq(cfg, p);
	    ESL_DPRINTF1(("MPI master: CM: %d exp tail mode: %d partition: %d\n", cfg->ncm, exp_mode, p));
	    /* estimate time for this round, assuming all workers have same processor speed as master */
	    if(p == 0) { 
	      if((status = estimate_time_for_exp_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf); 
	      psec *= expN * cfg->expL; /* psec was per residue */
	      if(cfg->nproc > 1) psec /= (cfg->nproc-1); /* parallelization will speed us up */
	    }
	    else psec = exp_psecAA[exp_mode][0];
	    exp_psecAA[exp_mode][p] = psec;
	    cm_psec    += psec;
	    total_psec += psec;
	    print_exp_line(go, cfg, errbuf, exp_mode, expN, cfg->expL, p, psec);

	    esl_stopwatch_Start(cfg->w_stage);
	    exp_scN  = 0;
	
	    if(xstatus == eslOK) have_work     = TRUE;	/* TRUE while work remains  */
	
	    wi = 1;
	    si   = -1;
	    need_seq = TRUE;
	    have_work = TRUE;	/* TRUE while work remains  */
	    seqpos = 1;
	
	    while (have_work || nproc_working)
	      {
		if (need_seq) 
		  {
		    need_seq = FALSE;
		    /* generate a new seq */
		    si++;
		    if(si < expN)
		      {
			/* generate sequence, either randomly from background null or from hard-wired 5 state HMM that emits genome like sequence */
			if(esl_opt_GetBoolean(go, "--exp-random")) { 
			  if((status = get_random_dsq(cfg, errbuf, cm, dnull, cfg->expL, &(dsq_slist[si]))) != eslOK) goto ERROR; 
			}
			else { 
			  if((status = get_genomic_sequence_from_hmm(cfg, errbuf, cm, cfg->expL, &(dsq_slist[si]))) != eslOK) goto ERROR;
			}
			results_slist[si] = CreateResults(INIT_RESULTS);
			sent_slist[si]    = FALSE;
			chunks_slist[si]  = 0;
			seqpos = 1;
			have_work = TRUE;
		      }
		    else if(si == expN) have_work = FALSE; 
		    else goto ERROR;
		  }
		if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0))
		  {
		    /* we're waiting to receive */
		    if (MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi probe failed");
		    if (MPI_Get_count(&mpistatus, MPI_PACKED, &n)                != 0) cm_Fail("mpi get count failed");
		    wi = mpistatus.MPI_SOURCE;
		    ESL_DPRINTF1(("MPI master sees a result of %d bytes from worker %d\n", n, wi));
		
		    if (n > bn) {
		      if ((buf = realloc(buf, sizeof(char) * n)) == NULL) cm_Fail("reallocation failed");
		      bn = n; 
		    }
		    if (MPI_Recv(buf, bn, MPI_PACKED, wi, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi recv failed");
		    ESL_DPRINTF1(("MPI master has received the buffer\n"));
		
		    /* If we're in a recoverable error state, we're only clearing worker results;
		     * just receive them, don't unpack them or print them.
		     * But if our xstatus is OK, go ahead and process the result buffer.
		     */
		    if (xstatus == eslOK) /* worker reported success. Get the result. */
		      {
			pos = 0;
			if (MPI_Unpack(buf, bn, &pos, &xstatus, 1, MPI_INT, MPI_COMM_WORLD)     != 0)     cm_Fail("mpi unpack failed");
			if (xstatus == eslOK) /* worker reported success. Get the results. */
			  {
			    si_recv = si_wlist[wi];
			    ESL_DPRINTF1(("MPI master sees that the result buffer contains search results\n"));
			    if ((status = cm_search_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_results)) != eslOK) cm_Fail("search results unpack failed");
			    ESL_DPRINTF1(("MPI master has unpacked search results\n"));
			
			    if(worker_results != NULL) { 
			      /* add results to seqlist[si_recv]->results[rclist[wi]] */
			      AppendResults(worker_results, results_slist[si_recv], seqpos_wlist[wi]);
			      /* careful, dbseqlist[si_recv]->results[rclist[wi]] now points to the traces and postal codes in worker_results->data,
			       * don't free those (don't use FreeResults(worker_results)) */
			      free(worker_results->data);
			      free(worker_results);
			      worker_results = NULL;
			    }
			    chunks_slist[si_recv]--;
			    if(sent_slist[si_recv] && chunks_slist[si_recv] == 0) 
			      { /* we're done with sequence si_recv; remove overlapping hits, copy scores of remaining hits, then free data */
				if(results_slist[si_recv]->num_results > 0) 
				  { 
				    RemoveOverlappingHits(results_slist[si_recv], 1, cfg->expL);
				    if(exp_scA == NULL) ESL_ALLOC (exp_scA,      sizeof(float) * (exp_scN + results_slist[si_recv]->num_results));
				    else                ESL_RALLOC(exp_scA, tmp, sizeof(float) * (exp_scN + results_slist[si_recv]->num_results));
				    for(h = 0; h < results_slist[si_recv]->num_results; h++) exp_scA[(exp_scN+h)] = results_slist[si_recv]->data[h].score;
				    exp_scN += results_slist[si_recv]->num_results;
				  }
				free(dsq_slist[si_recv]);
				FreeResults(results_slist[si_recv]);
				dsq_slist[si_recv] = NULL;
				results_slist[si_recv] = NULL;
			      }
			  }
			else	/* worker reported an error. Get the errbuf. */
			  {
			    if (MPI_Unpack(buf, bn, &pos, errbuf, cmERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack of errbuf failed");
			    ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
			    have_work = FALSE;
			    wi_error  = wi;
			  }
		      }
		    nproc_working--;
		  }
		if (have_work)
		  {   
		    /* send new search job */
		    len = (chunksize < (cfg->expL - seqpos + 1)) ? chunksize : (cfg->expL - seqpos + 1);
		    ESL_DPRINTF1(("MPI master is sending sequence i0..j0 %d..%d to search to worker %d\n", seqpos, seqpos+len-1, wi));
		    assert(seqpos > 0);
		
		    if ((status = cm_dsq_MPISend(dsq_slist[si]+seqpos-1, len, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
		
		    si_wlist[wi]     = si;
		    seqpos_wlist[wi] = seqpos;
		    len_wlist[wi]    = len;
		    chunks_slist[si]++;
		
		    wi++;
		    nproc_working++;
		
		    if(len == chunksize) seqpos += len - cm->W + 1;
		    else {
		      need_seq       = TRUE;
		      sent_slist[si] = TRUE; /* we've sent all chunks from this seq */
		    }
		  }
	      } 
	    ESL_DPRINTF1(("MPI master: done with this partition for this exp tail mode. Telling all workers\n"));
	    for (wi = 1; wi < cfg->nproc; wi++) 
	      if ((status = cm_dsq_MPISend(NULL, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");
	
	    /* fill a histogram with the exp_scN scores in exp_scA and fit it to an exponential tail */
	    if(cfg->exptfitfp != NULL) { 
	      fprintf(cfg->exptfitfp, "# CM: %s\n", cm->name);
	      fprintf(cfg->exptfitfp, "# mode: %12s\n", DescribeExpMode(exp_mode));
	    }
	    if((status = fit_histogram(go, cfg, errbuf, exp_scA, exp_scN, exp_mode, &tmp_mu, &tmp_lambda, &tmp_nrandhits, &tmp_tailp)) != eslOK) cm_Fail(errbuf);
	    SetExpInfo(cfg->cmstatsA[cmi]->expAA[exp_mode][p], tmp_lambda, tmp_mu, (long) (cfg->expL * expN), tmp_nrandhits, tmp_tailp);
	
	    for(si = 0; si < expN; si++) {
	      ESL_DASSERT1((dsq_slist[si] == NULL));
	      ESL_DASSERT1((results_slist[si] == NULL));
	      ESL_DASSERT1((chunks_slist[si] == 0));
	      ESL_DASSERT1((sent_slist[si] == TRUE));
	    }
	    esl_stopwatch_Stop(cfg->w_stage);
	    exp_asecAA[exp_mode][p] = cfg->w_stage->elapsed;
	    cm_asec += cfg->w_stage->elapsed;
	    total_asec += cfg->w_stage->elapsed;
	    FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
	    printf("  %10s\n", time_buf);

	    free(exp_scA); 
	    exp_scA = NULL;
	  } /* end of for(p = 0; p < cfg->np; p++) */
	free(dsq_slist);
	free(results_slist);
	free(chunks_slist);
	free(sent_slist);

	/****************************/
	/* filter threshold section */
	/****************************/
	if(exp_mode == EXP_CM_GI || exp_mode == EXP_CM_LI) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  if((status = estimate_time_for_fil_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf);
	  psec *= filN;
	  if(cfg->nproc > 1) psec /= (cfg->nproc-1); /* parallelization will speed us up */
	  fil_psecA[exp_mode] = psec;
	  cm_psec    += psec;
	  total_psec += psec;
	  print_fil_line(go, cfg, errbuf, exp_mode, psec);

	  esl_stopwatch_Start(cfg->w_stage);
	  fthr_mode = ExpModeToFthrMode(exp_mode);
	  ESL_DPRINTF1(("MPI master: CM: %d fthr mode: %d\n", cfg->ncm, fthr_mode));
      
	  if(xstatus == eslOK) have_work = TRUE;  /* TRUE while work remains  */
	  else                 have_work = FALSE; /* we've seen an error and are trying to finish cleanly */
      
	  wi = 1;
      
	  nseq_sent = 0;
	  nseq_recv = 0;
	  while (have_work || nproc_working)
	    {
	      if(have_work) { 
		if(nseq_sent < filN) {
		  nseq_this_worker = ((nseq_sent + fil_nseq_per_worker) <= filN) ? 
		    fil_nseq_per_worker : (filN - nseq_sent);
		}
		else { 
		  have_work = FALSE;
		  ESL_DPRINTF1(("MPI master has run out of numbers of sequences to dole out (%d doled)\n", nseq_sent));
		}
	      }
	      if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0)) {
		/* we're waiting to receive */
		if (MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi probe failed");
		if (MPI_Get_count(&mpistatus, MPI_PACKED, &n)                != 0) cm_Fail("mpi get count failed");
		wi = mpistatus.MPI_SOURCE;
		ESL_DPRINTF1(("MPI master sees a result of %d bytes from worker %d\n", n, wi));
	      
		if (n > bn) {
		  if ((buf = realloc(buf, sizeof(char) * n)) == NULL) cm_Fail("reallocation failed");
		  bn = n; 
		}
		if (MPI_Recv(buf, bn, MPI_PACKED, wi, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi recv failed");
		ESL_DPRINTF1(("MPI master has received the buffer\n"));
		
		/* If we're in a recoverable error state, we're only clearing worker results;
		 * just receive them, don't unpack them or print them.
		 * But if our xstatus is OK, go ahead and process the result buffer.
		 */
		if (xstatus == eslOK) /* worker reported success. Get the result. */
		  {
		    pos = 0;
		    if (MPI_Unpack(buf, bn, &pos, &xstatus, 1, MPI_INT, MPI_COMM_WORLD)     != 0)     cm_Fail("mpi unpack failed");
		    if (xstatus == eslOK) /* worker reported success. Get the results. */
		      {
			ESL_DPRINTF1(("MPI master sees that the result buffer contains HMM filter results\n"));
			if ((status = cmcalibrate_filter_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &wkr_fil_cyk_scA, &wkr_fil_ins_scA, &wkr_fil_fwd_scA, &wkr_fil_partA, &nseq_just_recv)) != eslOK) cm_Fail("cmcalibrate results unpack failed");
			ESL_DPRINTF1(("MPI master has unpacked HMM filter results\n"));
			for(i = 0; i < nseq_just_recv; i++) {
			  fil_cyk_scA[nseq_recv+i] = wkr_fil_cyk_scA[i];
			  fil_ins_scA[nseq_recv+i] = wkr_fil_ins_scA[i];
			  fil_fwd_scA[nseq_recv+i] = wkr_fil_fwd_scA[i];
			  fil_partA[nseq_recv+i]   = wkr_fil_partA[i];
			  ESL_DASSERT1((fil_partA[nseq_recv+i] < cfg->np));
			}
			free(wkr_fil_cyk_scA);
			free(wkr_fil_ins_scA);
			free(wkr_fil_fwd_scA);
			free(wkr_fil_partA);
			nseq_recv += nseq_just_recv;
		      }
		    else	/* worker reported an error. Get the errbuf. */
		      {
			if (MPI_Unpack(buf, bn, &pos, errbuf, cmERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack of errbuf failed");
			ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
			have_work = FALSE;
			wi_error  = wi;
		      }
		  }
		nproc_working--;
	      }

	      if (have_work)
		{   
		  /* send new search job */
		  ESL_DPRINTF1(("MPI master is sending HMM filter nseq %d to worker %d\n", nseq_this_worker, wi));
		  MPI_Send(&(nseq_this_worker), 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	      
		  wi++;
		  nproc_working++;
		  nseq_sent += nseq_this_worker;
		}
	    }

	  if(xstatus == eslOK) { 
	    exp_cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? EXP_CM_LC  : EXP_CM_GC;
	    exp_cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? EXP_CM_LI  : EXP_CM_GI;
	    fil_cm_cyk_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? FTHR_CM_LC : FTHR_CM_GC;
	    fil_cm_ins_mode  = (cm->flags & CMH_LOCAL_BEGIN) ? FTHR_CM_LI : FTHR_CM_GI;
	    /* set cutoffs for forward HMM filters, first for CYK, then for Inside */
	    if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_cyk_scA, fil_fwd_scA, fil_partA, exp_cm_cyk_mode, cfg->cmstatsA[cmi]->hfiA[fil_cm_cyk_mode])) != eslOK) cm_Fail(errbuf);
	    if((status = get_hmm_filter_cutoffs(go, cfg, errbuf, cm, fil_ins_scA, fil_fwd_scA, fil_partA, exp_cm_ins_mode, cfg->cmstatsA[cmi]->hfiA[fil_cm_ins_mode])) != eslOK) cm_Fail(errbuf);
	  }
	  ESL_DPRINTF1(("MPI master: done with HMM filter calc for fthr mode %d for this CM.\n", fthr_mode));
	  
	  for (wi = 1; wi < cfg->nproc; wi++) { 
	    msg = MPI_FINISHED_FILTER;
	    MPI_Send(&msg, 1, MPI_INT, wi, 0, MPI_COMM_WORLD);
	  }
	  esl_stopwatch_Stop(cfg->w_stage);
	  fil_asecA[exp_mode] = cfg->w_stage->elapsed;
	  cm_asec +=  cfg->w_stage->elapsed;
	  total_asec += cfg->w_stage->elapsed;
	  FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
	  printf("  %10s\n", time_buf);
	}
	ESL_DPRINTF1(("MPI master: done with exp tail mode %d for this CM.\n", exp_mode));
      }
      ESL_DPRINTF1(("MPI master: done with this CM.\n"));
      if(xstatus == eslOK) if(cfg->be_verbose) if((status = debug_print_cmstats(cm, errbuf, cfg->cmstatsA[cmi], TRUE)) != eslOK) cm_Fail(errbuf);
      print_per_cm_summary(go, cfg, errbuf, cm, cm_psec, cm_asec);
      if(esl_opt_IsDefault(go, "--forecast")) { if((status = print_post_calibration_info(go, cfg, errbuf, stdout, cm, exp_psecAA, fil_psecA, exp_asecAA, fil_asecA)) != eslOK) cm_Fail(errbuf); }
      free(fil_cyk_scA);
      free(fil_ins_scA);
      free(fil_fwd_scA);
      free(fil_partA);
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) { 
	free(exp_asecAA[exp_mode]); 
	free(exp_psecAA[exp_mode]); 
      }
      free(exp_asecAA);
      free(exp_psecAA);
      free(fil_asecA); 
      free(fil_psecA); 
      FreeCM(cm);
      printf("//\n");
      fflush(stdout);
    }
  free(seed_wlist);
  free(si_wlist);
  free(seqpos_wlist);
  free(len_wlist);
  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if     (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); return xstatus; }
  else if(status != eslEOF) return status;
  else                      return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "mpi_master() memory allocation error.");
  return eslOK; /* NOTREACHED */
}


static int
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;                /* Easel status */
  int      xstatus = eslOK;       /* changes from OK on recoverable error */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm = NULL;             /* the CM */
  int      cmi;                   /* CM index, which model we're working on */
  int      p;                     /* partition index */
  char    *wbuf = NULL;	          /* packed send/recv buffer  */
  int      wn   = 0;	          /* allocation size for wbuf */
  int      sz, n;		  /* size of a packed message */
  int      pos;                   /* posn in wbuf */
  long     seed;                  /* seed for RNG, rec'd from master */
  int      nseq;                  /* number of seqs to emit/search for current job */
  void    *tmp;                   /* ptr for ESL_RALLOC */ 
  MPI_Status mpistatus;           /* MPI status... */

  /* exponential tail related vars */
  ESL_DSQ          *exp_dsq = NULL;     /* dsq chunk received from master */
  int               expL;               /* length of exp_dsq received from master */
  search_results_t *exp_results = NULL; /* hits found in exp_dsq, sent back to master */
  int               exp_mode  = 0;

  /* filter threshold related vars */
  int      fthr_mode = 0;         /* CM mode for filter threshold calculation, FTHR_CM_GC, FTHR_CM_GI, FTHR_CM_LC, FTHR_CM_LI */
  float   *fil_cyk_scA = NULL;    /* [0..nseq-1] best cm cyk score for each emitted seq */
  float   *fil_ins_scA = NULL;    /* [0..nseq-1] best cm insidei score for each emitted seq */
  float   *fil_fwd_scA = NULL;    /* [0..nseq-1] best cp9 Forward score for each emitted seq */
  int     *fil_partA   = NULL;    /* [0..nseq-1] partition of CM emitted seq */
  int      in_fil_section_flag = FALSE; /* set to TRUE while we're in the filter threshold calculation
					 * section, we need to know this when we goto ERROR, b/c we have
					 * to know how many more MPI_Recv() calls to make to match up
					 * with the Master's sends before we can shut down.
					 */

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
	   
  /* Master now sends worker initialization information (RNG seed) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  if (MPI_Recv(&seed, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  if (xstatus == eslOK) { if((cfg->r = esl_randomness_Create(seed)) == NULL)          xstatus = eslEMEM; }
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return xstatus; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized seed: %ld\n", cfg->my_rank, seed));

  /* 2 special (annoying) cases: 
   * case 1: if we've used the --exp-gc option, we read in a seq file to fill
   * cfg->gc_freq, and we need that info here for the worker, so we receive
   * it's broadcast from the master
   * 
   * case 2: if we are calculating stats for more than 1 
   * partition, (--exp-pfile invoked), we need to receive that information 
   * via broadcast from master. Otherwise we need to setup the default partition info
   * (single partition, 0..100 GC content)
   */
  if(! (esl_opt_IsDefault(go, "--exp-gc"))) { /* receive gc_freq info from master */
    ESL_DASSERT1((cfg->gc_freq == NULL));
    ESL_ALLOC(cfg->gc_freq,  sizeof(double) * GC_SEGMENTS);
    ESL_ALLOC(cfg->pgc_freq, sizeof(double) * GC_SEGMENTS);
    MPI_Bcast(cfg->gc_freq, GC_SEGMENTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else cfg->gc_freq = NULL; /* default */
  if(! (esl_opt_IsDefault(go, "--exp-pfile"))) { /* receive partition info from master */
    MPI_Bcast(&(cfg->np),     1, MPI_INT, 0, MPI_COMM_WORLD);
    ESL_DASSERT1((cfg->pstart == NULL));
    ESL_ALLOC(cfg->pstart, sizeof(int) * cfg->np);
    MPI_Bcast(cfg->pstart, cfg->np, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else { /* no --exp-pfile, set up default partition info */  
    cfg->np     = 1;
    ESL_ALLOC(cfg->pstart, sizeof(int) * cfg->np);
    cfg->pstart[0] = 0;
  }
  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      cfg->ncm++;  
      if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
	cfg->cmalloc  += 128;
	ESL_RALLOC(cfg->cmstatsA, tmp, sizeof(CMStats_t *) * cfg->cmalloc);
      }
      cmi = cfg->ncm-1;
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status = initialize_cm(go, cfg, errbuf, cm))               != eslOK) goto ERROR;
      if((status = initialize_cmstats(go, cfg, errbuf, cm))          != eslOK) goto ERROR;
      if((status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len))) != eslOK) goto ERROR;
      
      for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {

	/* do we need to switch from glocal configuration to local? */
	if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	  if((status = switch_global_to_local(go, cfg, cm, errbuf))      != eslOK) goto ERROR;
	  if((status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len))) != eslOK) goto ERROR;
	}
	/* update search info for round 0 (final round) for exp tail mode */
	UpdateSearchInfoForExpMode(cm, 0, exp_mode);

	/************************************/
	/* exponential tail fitting section */
	/************************************/
	for (p = 0; p < cfg->np; p++) { /* for each partition */
	  ESL_DPRINTF1(("worker %d exp_mode: %d partition: %d\n", cfg->my_rank, exp_mode, p));

	  while((status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &exp_dsq, &expL)) == eslOK)
	    {
	      ESL_DPRINTF1(("worker %d: has received dsq chunk of length L: %d\n", cfg->my_rank, expL));
	      if ((status = ProcessSearchWorkunit(cm, errbuf, exp_dsq, expL, &exp_results, esl_opt_GetReal(go, "--mxsize"), cfg->my_rank, NULL, NULL)) != eslOK) goto ERROR;
	      RemoveOverlappingHits(exp_results, 1, expL);
	      ESL_DPRINTF1(("worker %d: has gathered search results\n", cfg->my_rank));

	      n = 0;
	      if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
		ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	      n += sz;
	      if (cm_search_results_MPIPackSize(exp_results, MPI_COMM_WORLD, &sz) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cm_search_results_MPIPackSize() call failed"); 
	      n += sz;  
		
	      if (n > wn) {
		void *tmp;
		ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
		wn = n;
	      }
	      ESL_DPRINTF1(("worker %d: has calculated the search results will pack into %d bytes\n", cfg->my_rank, n));
	      status = eslOK;
		
	      pos = 0;
	      if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
		ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	      if (cm_search_results_MPIPack(exp_results, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cm_search_results_MPIPack() call failed"); 
	      MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	      ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));
		
	      FreeResults(exp_results);
	      free(exp_dsq);
	    }
	  ESL_DPRINTF1(("worker %d exp_mode: %d finished partition: %d\n", cfg->my_rank, exp_mode, p));
	}
	ESL_DPRINTF1(("worker %d finished all partitions for exp_mode: %d\n", cfg->my_rank, exp_mode));

	
	/****************************/
	/* filter threshold section */
	/****************************/
	if(exp_mode == EXP_CM_GI || exp_mode == EXP_CM_LI) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
	  in_fil_section_flag = TRUE;
	  fthr_mode = ExpModeToFthrMode(exp_mode);
	  ESL_DPRINTF1(("worker %d fthr_mode: %d\n", cfg->my_rank, fthr_mode));

	  if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	  while(nseq != MPI_FINISHED_FILTER) {
	    ESL_DPRINTF1(("worker %d: has received hmm filter nseq: %d\n", cfg->my_rank, nseq));
	    
	    if((status = process_filter_workunit (go, cfg, errbuf, cm, nseq, &fil_cyk_scA, &fil_ins_scA, &fil_fwd_scA, &fil_partA)) != eslOK) cm_Fail(errbuf);
	    ESL_DPRINTF1(("worker %d: has gathered HMM filter results\n", cfg->my_rank));
	    n = 0;

	    if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	      ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	    n += sz;
	    if(cmcalibrate_filter_results_MPIPackSize(nseq, MPI_COMM_WORLD, &sz) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPackSize() call failed"); 
	    n += sz;  
	    if (n > wn) {
	      void *tmp;
	      ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	      wn = n;
	    }
	    ESL_DPRINTF1(("worker %d: has calculated the HMM filter results will pack into %d bytes\n", cfg->my_rank, n));
	    status = eslOK;
	    pos = 0;

	    if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	      ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	      if (cmcalibrate_filter_results_MPIPack(fil_cyk_scA, fil_ins_scA, fil_fwd_scA, fil_partA, nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
		ESL_XFAIL(eslFAIL, errbuf, "cmcalibrate_cp9_filter_results_MPIPack() call failed"); 
	    free(fil_cyk_scA);
	    free(fil_ins_scA);
	    free(fil_fwd_scA);
	    free(fil_partA);

	    MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	    ESL_DPRINTF1(("worker %d: has sent CP9 filter results to master in message of %d bytes\n", cfg->my_rank, pos));
	    /* receive next number of sequences, if MPI_FINISHED_EXPTAIL, we'll stop */
	    if(MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XFAIL(eslESYS, errbuf, "mpi recv failed");
	  }
	  in_fil_section_flag = FALSE;
	}
      } /* end of for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) */

      FreeCM(cm);
      cm = NULL;
      ESL_DPRINTF1(("worker %d finished all exp_modes for this cm.\n", cfg->my_rank));
    }
  if (status == eslEOD) ESL_DPRINTF1(("Worker %d told CMs are done.\n", cfg->my_rank));
  else goto ERROR;
  
  if (wbuf != NULL) free(wbuf);
  return eslOK;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  if(status == eslEMEM) sprintf(errbuf, "Memory allocation error.");
  MPI_Pack(&status, 1,               MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

  /* if we get here this worker failed and sent an error message, now the master knows a worker
   * failed but it has to continue through the mpi_master() code, sending the messages that
   * the workers expect, telling them to continue to move through the loops in those functions.
   * Minimal work will be done, but this is necessary so that we shut down cleanly. 
   * Because the master is sending messages to us still, we have to receive them. We can't
   * check that they're the expected messages though (codes telling us to keep moving through
   * the loops) because even if they were the wrong messages we couldn't do anything about it,
   * we've already entered error mode.
   */
  if(in_fil_section_flag) MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  for(; exp_mode < EXP_NMODES; exp_mode++) {
    MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
    if(ExpModeIsForCM(exp_mode)) {
      MPI_Recv(&nseq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
    }
  }
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);

  return eslFAIL; /* recoverable error, master has error message and will print it */
}
#endif /*HAVE_MPI*/

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;
  int nstarts, nexits, nd;
  float exp_cutoff;

  /* config QDB? yes, unless --exp-no-qdb enabled */
  if(esl_opt_GetBoolean(go, "--exp-no-qdb")) { 
    cm->search_opts |= CM_SEARCH_NOQDB; /* don't use QDB to search */
    /* cm->beta_qdb == cm->beta_W, both will be set as beta_W read from cmfile */
  }
  else {
    cm->config_opts |= CM_CONFIG_QDB;   /* configure QDB */
    cm->beta_qdb = esl_opt_GetReal(go, "--exp-beta"); 
  }

  /* set aggregate local begin/end probs, set with --pbegin, --pend, defaults are DEFAULT_PBEGIN, DEFAULT_PEND */
  cm->pbegin = esl_opt_GetReal(go, "--pbegin");
  cm->pend   = esl_opt_GetReal(go, "--pend");
  /* possibly overwrite local begin probs such that all begin points are equiprobable (--pebegin) */
  if(esl_opt_GetBoolean(go, "--pebegin")) {
    nstarts = 0;
    for (nd = 2; nd < cm->nodes; nd++) 
      if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
	nstarts++;
    cm->pbegin = 1.- (1./(1+nstarts));
  }
  /* possibly overwrite cm->pend so that local end prob from all legal states is fixed,
   * this is strange in that cm->pend may be placed as a number greater than 1., this number
   * is then divided by nexits in ConfigLocalEnds() to get the prob for each v --> EL transition,
   * this is guaranteed by the way we calculate it to be < 1.,  it's the argument from --pfend */
  if(! esl_opt_IsDefault(go, "--pfend")) {
    nexits = 0;
    for (nd = 1; nd < cm->nodes; nd++) {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	nexits++;
    }
    cm->pend = nexits * esl_opt_GetReal(go, "--pfend");
  }
  /* process the --fil-gemit option, this option forces all emitted parsetrees to be 'global'
   * in that they'll never contain a local begin or local end. */
  if(esl_opt_GetBoolean(go, "--fil-gemit")) { 
    cm->flags |= CM_EMIT_NO_LOCAL_BEGINS; 
    cm->flags |= CM_EMIT_NO_LOCAL_ENDS;
  }
  cm->search_opts |= CM_SEARCH_NOALIGN;

  if(! esl_opt_GetBoolean(go, "--no-null3")) cm->search_opts |= CM_SEARCH_NULL3;

  /* ALWAYS use the greedy overlap resolution algorithm to return hits for exp calculation
   * it's irrelevant for filter threshold stats, we return best score per seq for that */
  cm->search_opts |= CM_SEARCH_CMGREEDY;
  cm->search_opts |= CM_SEARCH_HMMGREEDY;

  if((status = ConfigCM(cm, errbuf, FALSE, NULL, NULL)) != eslOK) return status; /* FALSE says do not calculate W unless nec b/c we're using QDBs */
  
  /* create and initialize scan info for CYK/Inside scanning functions */
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("initialize_cm(), CreateScanMatrixForCM() call failed.");

  /* create the search info, which holds the thresholds for final round */
  if(esl_opt_IsDefault(go, "--exp-T")) exp_cutoff = -eslINFINITY;
  else exp_cutoff = esl_opt_GetReal(go, "--exp-T");
  CreateSearchInfo(cm, SCORE_CUTOFF, exp_cutoff, -1.);
  ValidateSearchInfo(cm, cm->si);
  
  if((status = update_dp_calcs(go, cfg, errbuf, cm)) != eslOK) return status;
  return eslOK;
}

/* initialize_cmstats()
 * Allocate and initialize a cmstats object in the cfg->cmstatsA array. 
 */
static int
initialize_cmstats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int i;
  int p;
  int cmi = cfg->ncm-1;

  ESL_DPRINTF1(("initializing cmstats for %d partitions\n", cfg->np));

  cfg->cmstatsA[cmi] = AllocCMStats(cfg->np);
  
  ESL_DASSERT1((cfg->pstart[0] == 0));
  for(p = 0; p < cfg->np;     p++) cfg->cmstatsA[cmi]->ps[p] = cfg->pstart[p];
  for(p = 0; p < (cfg->np-1); p++) cfg->cmstatsA[cmi]->pe[p] = cfg->pstart[p+1]-1;
  cfg->cmstatsA[cmi]->pe[(cfg->np-1)] = GC_SEGMENTS-1; /* this is 100 */
  
  for(p = 0; p < cfg->np; p++)
    for(i = cfg->cmstatsA[cmi]->ps[p]; i <= cfg->cmstatsA[cmi]->pe[p]; i++)
      cfg->cmstatsA[cmi]->gc2p[i] = p; 
  return eslOK;
}

/* Function: set_partition_gc_freq()
 * Date:     EPN, Mon Sep 10 08:00:27 2007
 *
 * Purpose:  Set up the GC freq to sample from for the current partition. 
 *           Only used if --exp-gc used to read in dbseq from which to derive
 *           GC distributions for >= 1 partition.
 *
 * Returns:  eslOK on success;
 */
int
set_partition_gc_freq(struct cfg_s *cfg, int p)
{
  int i, begin, end;
  ESL_DASSERT1((cfg->pgc_freq != NULL));
  ESL_DASSERT1((cfg->gc_freq != NULL));

  esl_vec_DSet(cfg->pgc_freq, GC_SEGMENTS, 0.);
  begin = cfg->pstart[p];
  if(p == (cfg->np-1)) end = (GC_SEGMENTS-1); /* this is 100 */
  else end = cfg->pstart[p+1] - 1;
  for (i = begin; i <= end; i++) 
    cfg->pgc_freq[i] = cfg->gc_freq[i];
  esl_vec_DNorm(cfg->pgc_freq, GC_SEGMENTS);

  return eslOK;
}


/* fit_histogram()
 * Create, fill and fit the tail of a histogram to an exponential tail. Data to fill the histogram
 * is given as <scores>.
 */
static int
fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, int exp_mode, double *ret_mu, double *ret_lambda, int *ret_nrandhits, float *ret_tailp)
{
  int status;
  double mu;
  double lambda;
  int i;
  double *xv;         /* raw data from histogram */
  int     n,z;  
  float tailp;
  double  params[2];
  int     nrandhits; 
  float   a;
  float   nhits_to_fit;

  ESL_HISTOGRAM *h = NULL;       /* histogram of scores */

  /* Initialize histogram; these numbers are guesses */
  if((h = esl_histogram_CreateFull(-100., 100., .1)) == NULL) return eslEMEM;    

  /* fill histogram */
  for(i = 0; i < nscores; i++) {
    if((status = esl_histogram_Add(h, scores[i])) != eslOK) ESL_FAIL(status, errbuf, "fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);
    /* printf("%4d %.3f\n", i, scores[i]); */
  }

  /* fit scores to an exponential tail */
  if(cfg->exptfitfp != NULL) { 
    /* fit to 41 different tailp values and print lambda, mu to a save file*/
    fprintf(cfg->exptfitfp, "# %11s  %10s  %10s  %12s\n", "tail pmass",  "lambda",     "mu",         "nhits");
    fprintf(cfg->exptfitfp, "# %11s  %10s  %10s  %12s\n", "-----------", "----------", "----------", "------------");
    for(a = 0.; a >= -4.; a -= 0.1) { 
      tailp = pow(10., a);
      esl_histogram_GetTailByMass(h, tailp, &xv, &n, &z); 
      if(n > 1) { 
	esl_exp_FitComplete(xv, n, &(params[0]), &(params[1]));
	esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params);
	fprintf(cfg->exptfitfp, "  %.9f  %10.6f  %10.4f  %12d\n", tailp, params[1], params[0], n);
      }
      else { 
	fprintf(cfg->exptfitfp, "  %.9f  %10s  %10s  %12d\n", tailp, "N/A", "N/A", n);
      }
    }
    fprintf(cfg->exptfitfp, "//\n");
  }
  /* end of if cfg->exptfitfp != NULL) */

  /* determine the fraction of the tail to fit, if --exp-tail-p, it's easy */
  if(!esl_opt_IsDefault(go, "--exp-tailp")) { 
    tailp = esl_opt_GetReal(go, "--exp-tailp");
    tailp = ESL_MIN(tailp, ((float) esl_opt_GetInteger(go, "--exp-tailxn") / (float) h->n)); /* ensure we don't exceed our max nhits in tail */
  }
  else { /* number of hits is per Mb and specific to local or glocal, CM or HMM fits */
    if(ExpModeIsLocal(exp_mode)) { 
      if(ExpModeIsForCM(exp_mode)) { /* local CM mode */
	nhits_to_fit = (float) esl_opt_GetInteger(go, "--exp-tailn-cloc") * ((cfg->exp_cmN_loc * cfg->expL) / 1000000.);
	tailp = nhits_to_fit / (float) h->n;
	if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--exp-tailn-cloc <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --exp-tailp.", esl_opt_GetInteger(go, "--exp-tailn-cloc"), (h->n / ((float) cfg->exp_cmN_loc * ((float) cfg->expL) / 1000000.)));
      }
      else { /* local HMM mode */
	nhits_to_fit = (float) esl_opt_GetInteger(go, "--exp-tailn-hloc") * ((cfg->exp_hmmN_loc * cfg->expL) / 1000000.);
	tailp = nhits_to_fit / (float) h->n;
	if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--exp-tailn-hloc <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --exp-tailp.", esl_opt_GetInteger(go, "--exp-tailn-hloc"), (h->n / ((float) cfg->exp_hmmN_loc * ((float) cfg->expL) / 1000000.)));
      }
    }
    else { 
      if(ExpModeIsForCM(exp_mode)) { /* glocal CM mode */
	nhits_to_fit = (float) esl_opt_GetInteger(go, "--exp-tailn-cglc") * ((cfg->exp_cmN_glc * cfg->expL) / 1000000.);
	tailp = nhits_to_fit / (float) h->n;
	if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--exp-tailn-cglc <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --exp-tailp.", esl_opt_GetInteger(go, "--exp-tailn-cglc"), (h->n / ((float) cfg->exp_cmN_glc * ((float) cfg->expL) / 1000000.)));
      }
      else { /* glocal HMM mode */
	nhits_to_fit = (float) esl_opt_GetInteger(go, "--exp-tailn-hglc") * ((cfg->exp_hmmN_glc * cfg->expL) / 1000000.);
	tailp = nhits_to_fit / (float) h->n;
	if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--exp-tailn-hglc <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --exp-tailp.", esl_opt_GetInteger(go, "--exp-tailn-hglc"), (h->n / ((float) cfg->exp_hmmN_glc * ((float) cfg->expL) / 1000000.)));
      }
    }
  }

  esl_histogram_GetTailByMass(h, tailp, &xv, &n, &z); /* fit to right 'tailfit' fraction, 0.01 by default */
  if(n <= 1) { 
    if(ExpModeIsLocal(exp_mode)) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), too few points in right tailfit: %f fraction of histogram. Increase --exp-cmL-loc or --exp-hmmLn-loc.", tailp);
    else                         ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), too few points in right tailfit: %f fraction of histogram. Increase --exp-cmL-glc or --exp-hmmLn-glc.", tailp);
  }
  esl_exp_FitComplete(xv, n, &(params[0]), &(params[1]));
  esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params);

  /* printf("# Exponential fit to %.7f%% tail: lambda = %f\n", tailp*100.0, params[1]); */
  mu = params[0];
  lambda = params[1];
  if(isnan(lambda)) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), exp tail fit lambda is NaN, too few hits in histogram. Increase --exp-cmL or --exp-hmmLn.");
  if(isinf(lambda)) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), exp tail fit lambda is inf, too few hits in histogram. Increase --exp-cmL or --exp-hmmLn.");
  nrandhits = h->n; /* total number of hits in the histogram */

  /* print to output files if nec */
  if(cfg->exphfp != NULL)
    esl_histogram_Plot(cfg->exphfp, h);
  if(cfg->expqfp != NULL) {
      esl_histogram_PlotQQ(cfg->expqfp, h, &esl_exp_generic_invcdf, params);
  }

  if (cfg->expsfp != NULL) {
    esl_histogram_PlotSurvival(cfg->expsfp, h);
    esl_exp_Plot(cfg->expsfp, (params[0] - log(1./tailp) / params[1]), 0.693147, esl_exp_surv, h->xmin - 5., h->xmax + 5., 0.1); /* extrapolate mu */
  }

  esl_histogram_Destroy(h);

  *ret_mu     = mu;
  *ret_lambda = lambda;
  *ret_nrandhits = nrandhits;
  *ret_tailp = tailp;
  return eslOK;
}

/* Function: get_random_dsq()
 * Date:     EPN, Tue Sep 11 08:31:47 2007
 * 
 * Purpose:  Generate a random digitized seq and return it.
 *           Two possible modes:
 *           1. if(cfg->pgc_freq == NULL && dnull != NULL) 
 *              use dnull disto (a double version of cm->null) to generate
 *           2. if(cfg->pgc_freq != NULL && dnull == NULL) 
 *              use choose a GC frequency from cfg->pgc_freq
 *              and generate with that
 *
 * Returns:  eslOK on success, ret_dsq filled with newly alloc'ed ESL_DSQ *,
 *           some other status code on failure.
 */
int
get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, double *dnull, int L, ESL_DSQ **ret_dsq)
{
  int status;
  double  gc_comp;
  double *distro = NULL;
  int do_free_distro = FALSE;
  ESL_DSQ *dsq = NULL;

  /* contract check, make sure we're in a valid mode */
  if(cfg->pgc_freq == NULL && dnull == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->pgc_freq == NULL and dnull == NULL");
  if(cfg->pgc_freq != NULL && dnull != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->pgc_freq != NULL and dnull != NULL");

  /* determine mode */ /* generate sequence */
  if      (cfg->pgc_freq == NULL && dnull != NULL) distro = dnull;
  else if (cfg->pgc_freq != NULL && dnull == NULL) {
    assert(cm->abc->K == 4);
    ESL_ALLOC(distro, sizeof(double) * cm->abc->K);
    do_free_distro = TRUE;
    gc_comp = 0.01 * esl_rnd_DChoose(cfg->r, cfg->pgc_freq, GC_SEGMENTS);
    distro[1] = distro[2] = 0.5 * gc_comp;
    distro[0] = distro[3] = 0.5 * (1. - gc_comp);
  }
  /* generate sequence */
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if ((status = esl_rsq_xIID(cfg->r, distro, cm->abc->K, L, dsq) != eslOK)) return status;

  if (do_free_distro) free(distro);
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  return status;
}

/* Function: get_cmemit_dsq()
 * Date:     EPN, Tue Sep 11 08:51:33 2007
 * 
 * Purpose:  Generate a dsq from a CM and return it.
 *
 * Returns:  eslOK on success, ESL_DSQ is filled with newly alloc'ed dsq; some other status code on an error, 
 */
int
get_cmemit_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int *ret_L, int *ret_p, ESL_DSQ **ret_dsq)
{
  int status;
  int p;
  int L;
  ESL_SQ *sq;
  ESL_DSQ *dsq;

  if((status = EmitParsetree(cm, errbuf, cfg->r, "irrelevant", TRUE, NULL, &sq, &L)) != eslOK) return status;
  while(L == 0) { 
    esl_sq_Destroy(sq); 
    if((status = EmitParsetree(cm, errbuf, cfg->r, "irrelevant", TRUE, NULL, &sq, &L)) != eslOK) return status;
  }

  /* determine the partition */
  p = cfg->cmstatsA[cfg->ncm-1]->gc2p[(get_gc_comp(cm->abc, sq->dsq, 1, L))]; /* this is slightly wrong, 1,L for get_gc_comp() should be i and j of best hit */
  assert(p < cfg->np);
  ESL_DASSERT1((p < cfg->np));

  /* free everything allocated by a esl_sqio.c:esl_sq_CreateFrom() call, but the dsq */
  dsq = sq->dsq;
  free(sq->name);
  free(sq->acc);
  free(sq->desc);
  free(sq);

  *ret_L  = L;
  *ret_p  = p;
  *ret_dsq = dsq;
  return eslOK;
}

/* Function: read_partition_file
 * Date:     EPN, Fri Dec  7 08:38:41 2007
 * 
 * Called when --exp-pfile is invoked. 
 * Opens and reads a partition file of 
 * with 2 * <npartitions> tokens, every odd token is
 * a partition start <pstart>, and every even token is 
 * a parititon end <pend>. First <pstart> must be 0,
 * other <pstart>s must be 1 more than previous
 * <pend>. The last <pend> must be 100, other <pends>
 * must be 1 less than following <pstart>.
 *
 * Example of file that implies 3 partitions: 
 * 0..39, 40..60, and 61.100
 * 
 * ~~~~~~~~~~~~~~~~
 * 0 39
 * 40 60
 * 61 100
 * ~~~~~~~~~~~~~~~~
 * 
 * After reading the file and checking it's legit,
 * set up the cfg->np and cfg->pstart data.
 *
 * Returns:  eslOK on success, eslEINVAL if file is 
 *           in wrong format, or doesn't follow rules described above.
 */
int
read_partition_file(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             toklen;
  int            *begin;
  int             end=0;
  int             nread=0;
  int             p;

  ESL_DASSERT1((MAX_PARTITIONS < GC_SEGMENTS));
  if(esl_opt_IsDefault(go, "--exp-pfile")) ESL_FAIL(eslEINVAL, errbuf, "read_partition_file, but --exp-pfile not invoked!\n");

  if (esl_fileparser_Open(esl_opt_GetString(go, "--exp-pfile"), NULL, &efp) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "failed to open %s in read_mask_file\n", esl_opt_GetString(go, "--exp-pfile"));
  esl_fileparser_SetCommentChar(efp, '#');
  
  ESL_ALLOC(begin, sizeof(int) * GC_SEGMENTS);
  begin[0] = 0;

  while((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslEOF) {
    begin[nread] = atoi(tok);
    if(nread == 0) {
      if(atoi(tok) != 0) ESL_FAIL(eslEINVAL, errbuf, "first partition begin must be 0 in %s\n", esl_opt_GetString(go, "--exp-pfile"));
    }
    else if (begin[nread] != (end+1)) {
      if(atoi(tok) != 0) ESL_FAIL(eslEINVAL, errbuf, "partition %d begin point (%d) is not exactly 1 more than prev partition end pt %d in %s\n", (nread+1), begin[nread], end, esl_opt_GetString(go, "--exp-pfile"));
    }      
    if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "no end point for each partition %d's begin (%d) in partition file %s\n", (nread+1), begin[nread], esl_opt_GetString(go, "--exp-pfile"));
    end = atoi(tok);
    if(end < begin[nread]) ESL_FAIL(eslEINVAL, errbuf, "partition %d end point (%d) < begin point (%d) in %s\n", (nread+1), end, begin[nread], esl_opt_GetString(go, "--exp-pfile"));
    nread++;
    if(nread > MAX_PARTITIONS) ESL_FAIL(eslEINVAL, errbuf, "partition file %s has at least %d partitions, but max num partitions is %d\n", esl_opt_GetString(go, "--exp-pfile"), nread, MAX_PARTITIONS);
  }
  if(nread == 0) ESL_FAIL(eslEINVAL, errbuf, "failed to read a single token from %s\n", esl_opt_GetString(go, "--exp-pfile"));
  if(end != 100) ESL_FAIL(eslEINVAL, errbuf, "final partitions end point must be 100, but it's %d in %s\n", end, esl_opt_GetString(go, "--exp-pfile"));

  /* create cfg->pstart */
  ESL_DASSERT1((cfg->pstart == NULL));
  ESL_ALLOC(cfg->pstart, sizeof(int) * nread);
  for(p = 0; p < nread; p++) cfg->pstart[p] = begin[p];
  free(begin);
  cfg->np = nread;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return status;
}

/* Function: switch_global_to_local()
 * Incept:   EPN, Mon Dec 10 08:43:32 2007
 * 
 * Purpose:  Switch a CM and it's CP9 HMM from global configuration
 *           to local configuration. Purposefully a local static function 
 *           in cmcalibrate.c, b/c we don't check if CM is in rsearch mode
 *           or any other jazz that'll never happen in cmcalibrate.
 *
 * Args:      go     - get opts
 *            cfg    - cmcalibrate's cfg
 *            cm     - the model
 *            errbuf - for printing errors
 *
 * Returns:   eslOK on succes, othewise some other easel status code and
 *            errbuf is filled with error message.
 */
int 
switch_global_to_local(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  int status;

  if(cm->flags & CMH_LOCAL_BEGIN) ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CMH_LOCAL_BEGIN flag already raised.\n");
  if(cm->flags & CMH_LOCAL_END)   ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CMH_LOCAL_END flag already raised.\n");
  if(! (cm->flags & CMH_CP9))     ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CMH_CP9 flag down.\n");
  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CPLAN9_LOCAL_BEGIN flag already raised.\n");
  if(cm->cp9->flags & CPLAN9_LOCAL_END)   ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CPLAN9_LOCAL_END flag already raised.\n");
  if(cm->cp9->flags & CPLAN9_EL)          ESL_FAIL(eslEINCOMPAT, errbuf, "switch_global_to_local(), CPLAN9_EL flag already raised.\n");

  /* ConfigLocal() puts CM in local mode, recalcs QDBs (if they exist), remakes cm's scan matrix, 
   * logoddsifies CM, and makes inserts equiprobable (if nec) */
  ConfigLocal(cm, cm->pbegin, cm->pend); 
  /* CPlan9SWConfig() configures CP9 for local alignment, then logoddisfies CP9 (wastefully in this case) */
  CPlan9SWConfig(cm->cp9, cm->pbegin, cm->pbegin, TRUE, cm->ndtype[1]);  /* TRUE means do make I_0, D_1, I_M unreachable to match the CM */
  /* CPlan9ELConfig() configures CP9 for CM EL local ends, then logoddisfies CP9 */
  CPlan9ELConfig(cm);

  /* recalculate cm->W and recalculate QDBs (if the CM has them) */
  if(cm->flags & CMH_QDB) { 
    free(cm->dmin); 
    free(cm->dmax); 
    cm->dmin = cm->dmax = NULL;
    cm->flags &= ~CMH_QDB;
    ConfigQDBAndW(cm, TRUE); /* TRUE says: calculate QDBs */
  }
  else ConfigQDBAndW(cm, FALSE); /* FALSE says: don't calc QDBs */
  /* this will create a new scan matrix, we need to update searchinfo so it points to this scan matrix */

  /* update cfg->fil_cm_ncalcs and cfg->cp9_ncalcs */
  if((status = update_dp_calcs(go, cfg, errbuf, cm)) != eslOK) return status;

  return eslOK;
}

/* Function: update_dp_calcs()
 * Incept:   EPN, Tue Jan 15 15:40:40 2008
 * 
 * Purpose:  Update cfg->fil_ncalcs and cfg->cp9_ncalcs
 *           based on configuration of CM.
 *
 * Args:      go     - get opts
 *            cfg    - cmcalibrate's cfg
 *            errbuf - for printing errors
 *            cm     - the model
 *
 * Returns:   eslOK on succes, othewise some other easel status code and
 *            errbuf is filled with error message.
 */
int 
update_dp_calcs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int  status;
  int *dmin; /* potentially used to calculate temporary d bands */
  int *dmax; /* potentially used to calculate temporary d bands */
  int  safe_windowlen; /* potentially used to calculate temporary d bands */

  /* Count number of DP calcs, these counts are used to determine target, minimum and maximum survival
   * fractions for HMM filter thresholds in get_hmm_filter_cutoffs() (see that code for details).
   * Our predicted running times for searches in get_hmm_filter_cutoffs() are calculated based
   * on cfg->fil_cm_calcs and cfg->fil_cp9_calcs. cfg->fil_cm_calcs is calculated assuming we'll
   * use a QDB filter with beta == cm->beta_qdb. 
   *
   * We want to set for global mode:
   * cfg->fil_cm_ncalcs:  millions of calcs for full CM scan of 1 residue, with QDBs using beta = cm->beta_qdb from cmfile
   * cfg->cp9_ncalcs:     millions of calcs for CP9 HMM scan of 1 residue
   *
   * This function is called twice. Once for global mode and then for local mode when models get localized.
   */

  /* get cfg->fil_cm_ncalcs, we ALWAYS assume that QDB CYK will be used as a second filter after HMM filtering
   * in cmsearch, so we determine target survival fractions based on DP counts for QDB searches, not non-QDB searches */
  safe_windowlen = cm->clen * 2;
  while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta_qdb, FALSE, &(dmin), &(dmax), NULL, NULL))) {
    free(dmin);
    free(dmax);
    safe_windowlen *= 2;
    if(safe_windowlen > (cm->clen * 1000))
      cm_Fail("initialize_cm(), safe_windowlen big: %d\n", safe_windowlen);
  }
  if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*cm->W, dmin, dmax, cm->W, TRUE,  NULL, &(cfg->fil_cm_ncalcs))) != eslOK) return status;
  free(dmin);
  free(dmax);
    
  /* get cfg->cp9_ncalcs, used to determine efficiency of CP9 filters, at first it's global mode, then
   * when switch_global_to_local() is called, cfg->full_cp9_ncalcs is updated to ncalcs in local mode */
  int cp9_ntrans = NHMMSTATETYPES * NHMMSTATETYPES; /* 3*3 = 9 transitions in global mode */
  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) cp9_ntrans++; 
  if(cm->cp9->flags & CPLAN9_LOCAL_END)   cp9_ntrans++; 
  if(cm->cp9->flags & CPLAN9_EL)          cp9_ntrans++; 
  cfg->cp9_ncalcs = (cp9_ntrans * cm->cp9->M) / 1000000.; /* convert to millions of calcs per residue */

  return eslOK;
}

/* Function: get_cmcalibrate_comlog_info
 * Date:     EPN, Mon Dec 31 14:59:52 2007
 *
 * Purpose:  Create the cmcalibrate command info and creation date info 
 *           to eventually be set in the CM's ComLog_t data structure.
 *
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation.
 */
static int
get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int i;
  long seed;
  long temp;
  int  seedlen;
  char *seedstr;

  if(cfg->ccom  != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->ccom  is non-NULL.");
  if(cfg->cdate != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->cdate is non-NULL.");
  
  
  /* Set the cmcalibrate command info, the cfg->ccom string */
  for (i = 0; i < go->optind; i++) { /* copy all command line options, but not the command line args yet, we may need to append '-s ' before the args */
    esl_strcat(&(cfg->ccom),  -1, go->argv[i], -1);
    esl_strcat(&(cfg->ccom),  -1, " ", 1);
  }
  /* if -s NOT enabled, we need to append the seed info also */
  seed = esl_randomness_GetSeed(cfg->r);
  if(esl_opt_IsDefault(go, "-s")) {
    temp = seed; 
    seedlen = 1; 
    while(temp > 0) { temp/=10; seedlen++; } /* determine length of stringized version of seed */
    seedlen += 4; /* strlen(' -s ') */
    ESL_ALLOC(seedstr, sizeof(char) * (seedlen+1));
    sprintf(seedstr, " -s %ld ", seed);
    esl_strcat((&cfg->ccom), -1, seedstr, seedlen);
    free(seedstr);
  }
  else { /* -s was enabled, we'll do a sanity check */
    if(seed != (long) esl_opt_GetInteger(go, "-s")) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "get_cmcalibrate_comlog_info(), cfg->r's seed is %ld, but -s was enabled with argument: %ld!, this shouldn't happen.", seed, (long) esl_opt_GetInteger(go, "-s"));
  }

  for (i = go->optind; i < go->argc; i++) { /* copy command line args yet */
    esl_strcat(&(cfg->ccom), -1, go->argv[i], -1);
    if(i < (go->argc-1)) esl_strcat(&(cfg->ccom), -1, " ", 1);
  }
  
  /* Set the cmcalibrate call date, the cfg->cdate string */
  time_t date = time(NULL);
  if((status = esl_strdup(ctime(&date), -1, &(cfg->cdate))) != eslOK) goto ERROR;
  esl_strchop(cfg->cdate, -1); /* doesn't return anything but eslOK */

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_cmcalibrate_comlog_info() error status: %d, probably out of memory.", status);
  return status; 
}


/* Function: update_comlog
 * Date:     EPN, Mon Dec 31 15:14:26 2007
 *
 * Purpose:  Update the CM's comlog info to reflect the current
 *           cmcalibrate call.
 *
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation.
 */
static int
update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm)
{
  int status;
  if(ccom  == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), ccom  is non-NULL.");
  if(cdate == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), cdate is non-NULL.");

  /* free all cmcalibrate comlog info, we're about to overwrite any information that any previous cmcalibrate
   * call could have written to the cm file.
   */
  if(cm->comlog->ccom  != NULL)  { free(cm->comlog->ccom);  cm->comlog->ccom = NULL;  }
  if(cm->comlog->cdate != NULL)  { free(cm->comlog->cdate); cm->comlog->cdate = NULL; }
  
  if((status = esl_strdup(ccom, -1, &(cm->comlog->ccom)))  != eslOK) goto ERROR; 
  if((status = esl_strdup(cdate,-1, &(cm->comlog->cdate))) != eslOK) goto ERROR; 
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "update_comlog() error status: %d, probably out of memory.", status);
  return status; 
}


/* Function: set_dnull
 * Date:     EPN, Thu Jan 24 09:48:54 2008
 *
 * Purpose:  Allocate, fill and return dnull, a double version of cm->null used
 *           for generating random seqs.
 *
 * Returns:  eslOK on success
 */
static int
set_dnull(CM_t *cm, char *errbuf, double **ret_dnull)
{
  int status;
  double *dnull;
  int i;

  if(ret_dnull == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "set_dnull(), ret_dnull is NULL.");
  ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
  for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
  esl_vec_DNorm(dnull, cm->abc->K);    
  *ret_dnull = dnull;

  return eslOK;

 ERROR:
  ESL_FAIL(eslEINCOMPAT, errbuf, "set_dnull(), memory allocation error.");
}


/* Function: process_filter_workunit()
 * Date:     EPN, Fri Jan 11 11:34:46 2008
 *
 * Purpose:  A filter work unit consists of a CM, an int specifying a 
 *           number of sequences <nseq>. The job is to generate <nseq> sequences 
 *           from the CM and search them first with the CM, both CYK and Inside
 *           and then with the HMM (Forward).
 *           Scores will eventually be used for calc'ing HMM filter thresolds.
 *
 * Args:     go             - getopts
 *           cfg            - cmcalibrate's configuration
 *           errbuf         - for writing out error messages
 *           cm             - the CM (already configured as we want it)
 *           nseq           - number of seqs to generate
 *           ret_cyk_scA    - RETURN: [0..nseq-1] best CM CYK score for each seq
 *           ret_ins_scA    - RETURN: [0..nseq-1] best CM Inside score for each seq
 *           ret_fwd_scA    - RETURN: [0..nseq-1] best CP9 Forward score for each seq
 *           ret_partA      - RETURN: [0..nseq-1] partition of each seq 
 *
 * Returns:  eslOK on success; dies immediately if some error occurs.
 */
static int
process_filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq,
			float **ret_cyk_scA, float **ret_ins_scA, float **ret_fwd_scA, int **ret_partA)
{
  int            status;
  float         *cyk_scA  = NULL;  /* [0..i..nseq-1] best CM CYK score for each state, each seq */
  float         *ins_scA  = NULL;  /* [0..i..nseq-1] best CM Inside score for each state, each seq */
  float         *fwd_scA = NULL;   /* [0..i..nseq-1] best CP9 Viterbi score for each seq */
  int           *partA  = NULL;    /* [0..i..nseq-1] partitions of each seq */
  int            p;                /* what partition we're in */
  int            i;
  int            L;
  ESL_DSQ       *dsq;
  int            orig_search_opts; /* we modify cm->search_opts in this function, then reset it at the end */
  int            do_null3;         /* TRUE to do NULL3 score corrections, FALSE not to */

  if(ret_cyk_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_cyk_scA == NULL.");
  if(ret_ins_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_ins_scA == NULL.");
  if(ret_fwd_scA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_fwd_scA == NULL.");
  if(ret_partA   == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "process_filter_workunit(), ret_partA == NULL.");

  ESL_DPRINTF1(("in process_filter_workunit nseq: %d\n", nseq));

  /* determine algs we'll use and allocate the score arrays we'll pass back */
  ESL_ALLOC(partA, sizeof(int) * nseq); /* will hold partitions */

  ESL_ALLOC(cyk_scA, sizeof(float) * nseq); /* will hold CM CYK scores */
  ESL_ALLOC(ins_scA, sizeof(float) * nseq); /* will hold CM Inside scores */
  ESL_ALLOC(fwd_scA, sizeof(float) * nseq);  /* will hold HMM Forward scores */

  orig_search_opts = cm->search_opts;
  do_null3 = (cm->search_opts & CM_SEARCH_NULL3) ? TRUE : FALSE;

  /* generate dsqs one at a time and collect optimal CM CYK/Inside scores and/or best CP9 Forward score */
  for(i = 0; i < nseq; i++) {
    if((status = get_cmemit_dsq(cfg, errbuf, cm, &L, &p, &dsq)) != eslOK) return status;
    partA[i] = p;
    /*to print seqs to stdout uncomment this block */
    /*ESL_SQ *tmp;
    tmp = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
    esl_sq_Textize(tmp);
    printf(">seq%d\n%s\n", i, tmp->seq);
    esl_sq_Destroy(tmp);
    fflush(stdout);
    */

    /* search dsq thrice, cyk, inside, fwd */
    /* for cyk and inside either use HMM bands, or don't */
    if(esl_opt_GetBoolean(go, "--fil-nonbanded")) { 
      cm->search_opts &= ~CM_SEARCH_INSIDE;
      if((status = FastCYKScan    (cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, do_null3, NULL, &(cyk_scA[i]))) != eslOK) return status; 
      
      cm->search_opts |= CM_SEARCH_INSIDE; 
      if((status = FastIInsideScan(cm, errbuf, cm->smx, dsq, 1, L, 0., NULL, do_null3, NULL, &(ins_scA[i]))) != eslOK) return status; 
    }
    else { /* search with HMM bands */
      cm->search_opts &= ~CM_SEARCH_INSIDE;
      cm->search_opts |= CM_SEARCH_HBANDED;
      cm->tau = esl_opt_GetReal(go, "--fil-tau");
      if(esl_opt_GetBoolean(go, "--fil-aln2bands")) cm->search_opts |= CM_SEARCH_HMMALNBANDS;
      if((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, TRUE, 0)) != eslOK) return status; 
      if((status = FastCYKScanHB(cm, errbuf, dsq, 1, L, 0., NULL, do_null3, cm->hbmx, esl_opt_GetReal(go, "--mxsize"), &(cyk_scA[i]))) != eslOK) return status; 

      cm->search_opts |= CM_SEARCH_INSIDE; 
      if((status = FastFInsideScanHB(cm, errbuf, dsq, 1, L, 0., NULL, do_null3, cm->hbmx, esl_opt_GetReal(go, "--mxsize"), &(ins_scA[i]))) != eslOK) return status; 
    }
    if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL, 
			     TRUE,   /* yes, we are scanning */
			     FALSE,  /* no, we are not aligning */
			     FALSE,  /* don't be memory efficient */
			     do_null3, 
			     NULL,   /* don't want best score at each posn back */
			     NULL,   /* don't want the max scoring posn back */
			     &(fwd_scA[i]))) != eslOK) return status;
    free(dsq);
  }
  /* contract enforced these are all non-NULL */
  *ret_cyk_scA = cyk_scA;
  *ret_ins_scA = ins_scA;
  *ret_fwd_scA = fwd_scA;
  *ret_partA   = partA;

  cm->search_opts = orig_search_opts;

  return eslOK;

 ERROR:
  return status;
}

typedef struct fseq_Eval_s {
  int i;
  float cm_E;
  float fwd_E;
} fseq_Eval_t;

/*
 * Function: compare_fseq_by_{cm,fwd}_Eval()
 * Date:     EPN, Fri Jan 11 12:30:59 2008
 * Purpose:  Compares two fseq_Eval_t's based on CM or HMM Forward E-value
 *           and returns -1 if first is higher E-value than second, 0 if equal, 
 *           1 if first E-value is lower.  This results in sorting by E-value, 
 *           highest first.
 */
int compare_fseq_by_cm_Eval(const void *a_void, const void *b_void) {
  fseq_Eval_t *a, *b;
  a = (fseq_Eval_t *) a_void;
  b = (fseq_Eval_t *) b_void;
  if      (a->cm_E > b->cm_E) return -1;
  else if (a->cm_E < b->cm_E) return  1;
  else                        return  0;
}

int compare_fseq_by_fwd_Eval(const void *a_void, const void *b_void) {
  fseq_Eval_t *a, *b;
  a = (fseq_Eval_t *) a_void;
  b = (fseq_Eval_t *) b_void;
  if      (a->fwd_E > b->fwd_E) return -1;
  else if (a->fwd_E < b->fwd_E) return  1;
  else                          return  0;
}

/* Function: get_hmm_filter_cutoffs()
 * Date:     EPN, Fri Jan 11 12:01:07 2008
 *
 * Purpose:  Given a CM and scores for a CM and HMM Forward scan of
 *           filN target seqs predict the HMM filter threshold for
 *           each possible CM threshold. The CM scores are either CYK
 *           or Inside; so this function is called four times, Once
 *           each for glocal CYK, glocal Inside, local CYK and glocal
 *           Inside.
 *           
 *           For the CM scores, the possible CM thresholds are the
 *           E-values for the first 90% (worst scoring 90%) observed
 *           CYK/Inside scores in a ranked list of such E-values,
 *           stored in sorted order in by_cmA[], E-value from
 *           i=0..filN-1. The first 90% are the elements i=0..imax,
 *           with imax = 0.90 * filN.
 *
 *           The HMM threshold fwd_E_cut[i] for each i=0..imax is the
 *           HMM Forward E-value that recognizes F fraction of the
 *           (filN-i+1) sequences that have a CM E-value better than
 *           by_cmA[i].  Note, that when i == imax, 25% of of the
 *           sequences (0.25 * filN) have an E-value of by_cmA[i].cm_E
 *           or lower, this means that the HMM filter for i == imax
 *           will recognize (F * filN * 0.25) sequences, so filN must
 *           be appropriately large so 0.25 * filN is a reasonable
 *           sample size.
 *
 *           To achieve this, we maintain two separately sorted lists,
 *           one by CM E-values (by_cmA) and one by Forward E-values
 *           (by_fwdA), with the cmi2fwdi array providing a map
 *           between the two lists. These two lists contain the same
 *           data (the CM and Forward scores for the filN CM sampled
 *           sequences), but in a different order. If cmi2fwdi[j] ==
 *           k, this means element with rank j in by_cmA has rank k
 *           in by_fwdA, specifically by_cmA[j].i == by_fwdA[k].i.
 *
 *           With these two sorted lists we can step through the CM
 *           list and at each different threshold point (i=0..imax),
 *           update the Forward list by setting specific elements
 *           (specifically those elements whose CM score is greater
 *           (worse) than by_cmA[i].cm_E) as 'no longer used'.  Thus
 *           at each threshold point i, we have the same curN ==
 *           filN-i+1 elements being 'used' in the by_cmA and by_fwdA
 *           lists. Concurrently, we keep track of what HMM threshold
 *           is necessary to recognize F fraction of those curN hits,
 *           specifically this is the Fidx'th ranking hit OF USED
 *           ELEMENTS in the by_fwdA array. This step-through loop is
 *           implemented in the section of code below marked "main
 *           loop"
 *           
 *           It is somewhat tricky to follow the implementation of
 *           this below, at least the way I've done it, and I've found
 *           it equally tricky to comment it in a clear way, but I
 *           tried to make it clear (without spending an unreasonable
 *           amount of time on it).
 *            
 * Args:     go -      command line options
 *           cfg -     cmcalibrate's cfg object
 *           errbuf -  for printing error messages
 *           cm -      the model
 *           cm_scA -  [0..i..filN-1] best CM score (CYK or Inside) in sequence i 
 *           fwd_scA - [0..i..filN-1] best Foward score in sequence i 
 *           partA   - [0..i..filN-1] partition of sequence i 
 *           cm_mode - CM mode that explain the configuration the CM was in when scores in cm_scA were collected
 *                     either EXP_CM_LC (local CYK), EXP_CM_LI (local Inside), EXP_CM_GC (glocal CYK), or EXP_CM_GI (glocal Inside)
 *           bf      - BestFilterInfo_t object, we'll update this to hold info on Forward filter cutoffs
 *
 * Returns:  Updates BestFilterInfo_t object <bf> to hold HMM filter info
 *           eslOK on success;
 *           Other easel status code on an error with errbuf filled with error message.
 */

int
get_hmm_filter_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_scA, float *fwd_scA, int *partA, int cm_mode, HMMFilterInfo_t *hfi)
{
  int    status;
  float  fil_ncalcs;              /* number of million dp calcs predicted for the HMM filter scan */
  float  cm_ncalcs;               /* number of million dp calcs predicted for a full (next-step-after-filter) CM scan */
  int    i, p;                    /* counters */
  int    cmi   = cfg->ncm-1;      /* CM index we're on */
  float  F     = esl_opt_GetReal   (go, "--fil-F");    /* fraction of CM seqs we require filter to let pass */
  int    filN  = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences we emitted from CM for filter test */
  float  xhmm  = esl_opt_GetReal   (go, "--fil-xhmm"); /* number of sequences we emitted from CM for filter test */
  int    Fidx;                    /* index in sorted scores that threshold will be set at (1-F) * N */
  float  surv_res_per_hit;        /* expected number of residues to survive filter from DB for each hit 2*W-avg_hit_len, twice W minus the average lenght of a hit from QDB calc */
  float  Smax;                    /* maximally useful survival fraction, close to 0.5 */
  float  Starg;                   /* our target survival fraction, if we achieve this fraction the search should take xhmm times longer than a HMM only search */
  float  Smin;                    /* minimally useful survival fraction, any less than this and our filter would be (predicted to be) doing more than 10X the work of the CM */
  float  fwd_Emax;                /* fwd E-value that gives Smax survival fraction */
  float  fwd_Etarg;               /* fwd E-value that gives Starg survival fraction */
  float  fwd_Emin;                /* fwd E-value that gives Smin survival fraction */
  double dbsize;                  /* database size the E-values refer to, just a scaling factor */
  int   *cmi2fwdi;                /* [0..i..filN] map between by_cmA and by_fwdA elements, cmi2fwdi[j] == k ==> by_cmA[j].i == by_fwdA[k].i */
  int    j;                       /* counter */
  int    always_better_than_Smax; /* TRUE if worst observed cm E value (by_cmA[0].cm_E) still gives us a survival fraction better (less) than Smax */
  
  fseq_Eval_t *by_cmA;            /* [0..i..filN] list of fseq_Eval_t for all filN observed seqs, sorted by decreasing CM      E-value */
  fseq_Eval_t *by_fwdA;           /* [0..i..filN] list of fseq_Eval_t for all filN observed seqs, sorted by decreasing Forward E-value */
  float        cm_E, fwd_E;       /* a CM E-value and Forward E-value */
  int          hmm_fwd_mode = (cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? EXP_CP9_LF : EXP_CP9_GF;

  /* contract checks */
  if(! (cfg->cmstatsA[cmi]->expAA[cm_mode][0]->is_valid))      ESL_FAIL(eslEINCOMPAT, errbuf, "get_hmm_filter_cutoffs(), exp tail stats for CM mode: %d are not valid.\n", cm_mode);
  if(! (cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][0]->is_valid)) ESL_FAIL(eslEINCOMPAT, errbuf, "get_hmm_filter_cutoffs(), exp tail stats for HMM fwd mode: %d are not valid.\n", hmm_fwd_mode);

  /* Determine our target, min and max E-value cutoffs for the HMM forward scan
   * these are independent of the observed scores. Each of these is calculated
   * based on how it affects <total_calcs>. <total_calcs> is the sum of the
   * number of DP calculations required for the HMM filter <fil_ncalcs>, plus the predicted
   * number of DP calculations required for the CM to search the survival fraction
   * <cm_ncalcs> * survival fract. So if <total_calcs> is <xhmm> * <fil_ncalcs>, the 
   * required time we predict to do a HMM filter plus CM scan of survivors is <xhmm> times the
   * time it would take to do an HMM only scan. 
   * 
   * fwd_Etarg:  the target  forward E-value cutoff. If used <total_calcs> = <xhmm> * <fil_ncalcs>
   *             <xhmm> is obtained from the --fil-xhmm, by default it is 2.0.
   * fwd_Emin:   the minimal forward E-value cutoff, anything less than this is overkill. 
   *             If used <total_calcs> = 1.1 * fil_ncalcs.
   * fwd_Emax:   the maximum forward E-value cutoff we'll accept as useful. If used 
   *             <total_calcs> = 0.5 * <cm_ncalcs> (we expect HMM filter to give us 
   *                                                   2X speedup versus full CM)
   *
   * forward filter   predicted 
   * E-value cutoff   survival fraction
   * --------------   -----------------
   * fwd_Emax         Smax  = 0.5 -  (fil_ncalcs / cm_ncalcs)
   * fwd_Etarg        Starg = xhmm * (fil_ncalcs / cm_ncalcs)
   * fwd_Emin         Smin  = 0.1  * (fil_ncalcs / cm_ncalcs)
   */

  dbsize           = cfg->cmstatsA[cmi]->expAA[cm_mode][0]->dbsize; 
  /* Update hmm E-value's effective database size as if we're searching a database of size dbsize */
  for(p = 0; p < cfg->cmstatsA[cmi]->np; p++) {
    cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][p]->cur_eff_dbsize = (long) ((((double) dbsize / (double) cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][p]->dbsize) * 
									  ((double) cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][p]->nrandhits)) + 0.5);
  }
  fil_ncalcs       = cfg->cp9_ncalcs;                      /* fil_ncalcs is millions of DP calcs for HMM Forward scan of 1 residue */
  fil_ncalcs      *= dbsize;                               /* now fil_ncalcs is millions of DP calcs for HMM Forward scan of length 1 Mb */
  cm_ncalcs        = cfg->fil_cm_ncalcs;                   /* total number of millions of DP calculations for full CM scan of 1 residue */
  cm_ncalcs       *= dbsize;                               /* now cm_ncalcs corresponds to dbsize (1 Mb)*/
  surv_res_per_hit = (2. * cm->W - (cfg->avg_hit_len));    /* avg length of surviving fraction of db from a single hit (cfg->avg_hit_len is avg subseq len in subtree rooted at v==0, from QDB calculation) */
  Smin             =  0.1      * (fil_ncalcs / cm_ncalcs); /* if survival fraction == Smin, total number of DP calcs (filter + survivors) == (1.1) * fil_ncalcs, so a HMM + CM search will take only 10% longer than HMM only */
  Starg            = (xhmm-1.) * (fil_ncalcs / cm_ncalcs); /* if survival fraction == Starg, the total number DP calcs (filter + survivors) == (1 + xhmm) * fil_ncalcs,
							    * so if xhmm = 2 (which is default), our target is for the HMM filter + CM search of survivors (with QDB possibly) 
							    * to take 2 times long as only the HMM filter would take. */ 
  Smax             = 0.5 -  (fil_ncalcs / cm_ncalcs);      /* if survival fraction == Smax, total number of DP calcs (filter+survivors) == 1/2 * cm_ncalcs, so our predicted speedup
							    * by using the filter is 2 fold */
  if(Smin > Smax) { /* rare case, happens only if fil_ncalcs >= 5/11 * cm_ncalcs, in this case we don't filter */
    /* this will probably never happen, but if it does, we set 1 cut point, set the CM E-value to 0.0 and set
     * always_better_than_Smax to FALSE, then in cmsearch no matter what CM E-value cutoff <E> is used, a filter will
     * not be used, because <E> > 0.0 for all valid <E>. 
     */
    float *cm_E_cut;
    float *fwd_E_cut;
    ESL_ALLOC(fwd_E_cut, sizeof(float) * 1);
    ESL_ALLOC(cm_E_cut,  sizeof(float) * 1);
    cm_E_cut[0]  = 0.0;         
    fwd_E_cut[0] = 1E20; /* this is irrelevant actually */
    always_better_than_Smax = FALSE;
    if((status = SetHMMFilterInfoHMM(hfi, errbuf, F, filN, dbsize, 1, cm_E_cut, fwd_E_cut, always_better_than_Smax)) != eslOK) return status;
    free(cm_E_cut);
    free(fwd_E_cut);
    return eslOK; 
  }
  if(Starg < Smin) { /* another rare case min value for --fil-xhmm is 0.1, so if Starg < Smin it's just because of precision issues, nevertheless we have to deal */
    Starg = Smin;
  }
  if(Starg > Smax) { /* yet another rare case, happens only if fil_ncalcs >= (0.5 * cm_ncalcs) / (xhmm+1), when xhmm = 1, this is when fil_ncalcs <= cm_ncalcs / 4 */
    Starg = Smax;
  }
  fwd_Emax         = (Smax  * (float) dbsize) / surv_res_per_hit;
  fwd_Etarg        = (Starg * (float) dbsize) / surv_res_per_hit;
  fwd_Emin         = (Smin  * (float) dbsize) / surv_res_per_hit;

  assert(fwd_Emax >  -0.000001);
  assert(fwd_Etarg > -0.000001);
  assert(fwd_Emin >  -0.000001);

  /* Copy bit scores to three separate quicksort-able data structures, 
   * Sort one by CM E-value, one by Inside E-value and one by Forward E-value.
   * The seq idx is kept in the data structure, providing a link between 
   * the three different lists after they're sorted. 
   */
  ESL_ALLOC(by_cmA, sizeof(fseq_Eval_t) * filN);
  ESL_ALLOC(by_fwdA, sizeof(fseq_Eval_t) * filN);
  /* convert bit scores to E-values and copy them to the qsortable structures */
  for(i = 0; i < filN; i++) { 
    p     = partA[i];
    cm_E  = Score2E(cm_scA[i],  cfg->cmstatsA[cmi]->expAA[cm_mode][p]->mu_extrap,      cfg->cmstatsA[cmi]->expAA[cm_mode][p]->lambda, cfg->cmstatsA[cmi]->expAA[cm_mode][p]->cur_eff_dbsize);
    fwd_E = Score2E(fwd_scA[i], cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][p]->mu_extrap, cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][p]->lambda, cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][p]->cur_eff_dbsize);
    /* copy E-values to qsortable data structures */
    by_cmA[i].i     = by_fwdA[i].i     = i;
    by_cmA[i].cm_E  = by_fwdA[i].cm_E  = cm_E;
    by_cmA[i].fwd_E = by_fwdA[i].fwd_E = fwd_E;
    /*printf("TEMP i: %5d CM sc: %.3f (E: %20.10f) HMM sc: %.3f (E: %20.10f)\n", i, cm_scA[i], cm_E, fwd_scA[i], fwd_E);*/
  }

  /* qsort */
  qsort(by_cmA,  filN, sizeof(fseq_Eval_t), compare_fseq_by_cm_Eval);
  qsort(by_fwdA, filN, sizeof(fseq_Eval_t), compare_fseq_by_fwd_Eval);

  /* determine the mapping between the sorted CM list and sorted Fwd list.
   * This is N^2 which might be of concern with large filN's 
   * empirically for filN = 10,000 this block takes 1.25s, 
   * so for filN == 100,000 it would take about 2 minutes. But
   * for filN == 100,000 the time search 100,000 seqs will drown
   * this out. Nevertheless, max <n> for --fil-N option is 100,000.
   */
  ESL_ALLOC(cmi2fwdi, sizeof(int) * filN);
  esl_vec_ISet(cmi2fwdi, filN, -1);
  for(i = 0; i < filN; i++) {
    for(j = 0; j < filN; j++) {
      if(by_cmA[i].i == by_fwdA[j].i) cmi2fwdi[i] = j;
    }
  }
  for(i = 0; i < filN; i++) assert(cmi2fwdi[i] != -1);

  /* A few more preparations before the main loop:
   * we want Fidx such that in the sorted arrays, by_cmA and by_fwdA, F fraction of the elements 
   * have E value <= score of element [Fidx], including element Fidx. Some examples:
   * F      filN    Fidx
   * 0.99    100       1  we'll miss 0..0   ==   1 hit  at this threshold, but recognize   99 (0.99)
   * 0.93   2000     140  we'll miss 0..139 == 140 hits at this threshold, but recognize 1860 (0.93)
   * 0.81    243      46  we'll miss 0..45  ==  46 hits at this threshold, but recognize  197 (0.8106)
   */
  Fidx = (int) ((1. - F) * (float) filN);
  /* deal with a precision issue */
  if(Fidx < ((1. - F) * (float) filN)) { /* we rounded down */
    if (((Fidx - ((1. - F) * (float) filN)) - 1.) < eslSMALLX1) Fidx++; /* due to precision issues we unnecessarily rounded down */
  }

  /* Setup fwd_useme array: fwd_useme[i] is TRUE if element with rank i in by_fwdA array is 'in use', 
   * meaning it corresponds to a sequence with CM score better than current threshold */
  int *fwd_useme;      /* [0..i..filN] TRUE if element with rank i in by_fwd array is 'in use', meaning it corresponds to a sequence with CM score better than current threshold */
  ESL_ALLOC(fwd_useme, sizeof(int) * filN);
  esl_vec_ISet(fwd_useme, filN, TRUE);
  int fwd_all = 0;     /* the index in by_fwdA of the worst ranked  (highest) forward E-value of all elements in use (fwd_useme[i] == TRUE) */
  int fwd_F   = Fidx;  /* the index in by_fwdA of the sequence with forward E-value worse than F fraction of all elements in use (fwd_useme[i] == TRUE) */

  /* Precalculate the i values for which we have to change fwd_F, as step through the i=0..imax loop.
   * when i == 0, curN == filN-i+1 == filN, and fwd_F == Fidx because we want to find F fraction of the filN hits
   * but as i increases, curN decreases (as we step through the ranked CM E-value threshold) and the integer 
   * corresponding  to F fraction of (filN-i+1) also changes, but due to precision difference of ints and floats,
   * it's non-trivial, and worth precalculating to avoid putting the messy code that does
   * it in the loop below. 
   */
  int curN = filN;
  int *change_fwd_F;
  int prv_Fidx;
  int cur_Fidx;
  ESL_ALLOC(change_fwd_F, sizeof(int) * (filN+1));
  prv_Fidx = Fidx;
  for(i = filN-1; i >= 0; i--) { 
    cur_Fidx = (int) ((1. - F) * (float) i);
    if(cur_Fidx < ((1. - F) * (float) i)) { /* we rounded down */
      if (((cur_Fidx - ((1. - F) * (float) i)) - 1.) < eslSMALLX1) cur_Fidx++; /* due to precision issues we unnecessarily rounded down */
    }
    if(prv_Fidx == cur_Fidx) change_fwd_F[i+1] = TRUE;
    else                     change_fwd_F[i+1] = FALSE;
    prv_Fidx = cur_Fidx;
  }
  change_fwd_F[0] = FALSE; /* we'll never access this element */

  /* The main loop:
   * step through 9/10 of the CM E-values starting with worst (so we'll stop at the hit with rank imax = (0.90 * filN)),
   * so we step through i=0..imax, and for each i, determine fwd_E_F[i], fwd_E_all[i], and fwd_E_cut[i]
   * fwd_E_F[i]:   the forward E-value cutoff that will recognize F fraction of CM hits with E-values better than or equal to by_cmA[i].cm_E
   * fwd_E_all[i]: the forward E-value cutoff that will recognize all           CM hits with E-values better than or equal to by_cmA[i].cm_E
   * fwd_E_cut[i]: the forward E-value cutoff we would report for this i, may be different than fwd_E_F[i] and fwd_E_all[i]
   *               because we consider our fwd_Emax, fwd_Etarg, and fwd_Emin values.
   * fwd_E_S[i]:   predicted survival fraction using fwd_E_cut[i], equals (fwd_E_cut[i] * surv_res_per_hit / dbsize).
   */
  int imax = (int) ((0.900001) * (float) filN); /* .900001 is to avoid unnec rounding down (for ex, if filN == 100, and we used 0.90, we may get 89.99999999 and end up rounding down) */
  float *fwd_E_F; 
  float *fwd_E_all;
  float *fwd_E_cut;
  float *fwd_E_S;  
  ESL_ALLOC(fwd_E_F,    sizeof(float) * (imax+1));
  ESL_ALLOC(fwd_E_all,  sizeof(float) * (imax+1));
  ESL_ALLOC(fwd_E_cut,  sizeof(float) * (imax+1));
  ESL_ALLOC(fwd_E_S,    sizeof(float) * (imax+1));
  int imax_above_Emax  = -1; /* at end of loop, imax_above_Emax is max i for which fwd_E_cut[i] >  fwd_Emax,  if -1, fwd_E_cut[i] < fwd_Emax  for all i=0..imax */
  int imin_at_Etarg    = -1; /* at end of loop, imin_at_Etarg   is min i for which fwd_E_cut[i] == fwd_Etarg, if -1, fwd_E_cut[i] > fwd_Etarg for all i=0..imax */
  int imax_at_Etarg    = -1; /* at end of loop, imax_at_Etarg   is max i for which fwd_E_cut[i] == fwd_Etarg, if -1, fwd_E_cut[i] > fwd_Etarg for all i=0..imax */
  int imin_at_Emin     = -1; /* at end of loop, imin_at_Emin    is min i for which fwd_E_cut[i] == fwd_Emin,  if -1, fwd_E_cut[i] > fwd_Emin  for all i=0..imax */

  ESL_DPRINTF1(("fwd_Emax:  %f\n", fwd_Emax));
  ESL_DPRINTF1(("fwd_Etarg: %f\n", fwd_Etarg));
  ESL_DPRINTF1(("fwd_Emin: %f\n", fwd_Emin));

  for(i = 0; i <= imax; i++) { /* the main loop */
    fwd_E_all[i] = by_fwdA[fwd_all].fwd_E; /* by_fwdA[fwd_all] is HMM threshold that recognizes ALL curN CM hits w/E value <= by_cmA[i].cm_E */
    fwd_E_F[i]   = by_fwdA[fwd_F].fwd_E;   /* by_fwdA[fwd_F]   is HMM threshold that recognizes F fraction of curN CM hits w/E value <= by_cmA[i].cm_E */

    /* based on fwd_E_all[i] and fwd_E_F[i], determine what our cutoff should be */
    if(fwd_E_F[i] < fwd_Etarg) { /* if TRUE, using E-value that achieves the target survival fraction as cutoff will recognize F fraction of CM hits */
      fwd_E_cut[i]    = ESL_MIN(fwd_Etarg,    fwd_E_all[i]); /* take minimum of: E value that exactly satisfies target survival fraction, and E value that recognizes all CM hits */
      fwd_E_cut[i]    = ESL_MAX(fwd_E_cut[i], fwd_Emin);     /* never go less than fwd_Emin */
    }
    else fwd_E_cut[i] = fwd_E_F[i];       /* we didn't achieve our target survival fraction */
    fwd_E_S[i] = (fwd_E_cut[i] * (float) surv_res_per_hit) / (float) dbsize;

    ESL_DPRINTF1(("i: %3d N: %5d cm: %g fwd_E_all: %g fwd_E_F: %g cut: %g S: %g\n", i, curN, by_cmA[i].cm_E, fwd_E_all[i], fwd_E_F[i], fwd_E_cut[i], fwd_E_S[i])); 

    if(fwd_E_cut[i]  > fwd_Emax) imax_above_Emax++; /* we didn't even achieve our maximum allowed survival fraction */
    if(imin_at_Emin  == -1 && (fabs(fwd_E_cut[i] - fwd_Emin)  < eslSMALLX1)) { ESL_DPRINTF1(("\tAchieved Smin\n"));   imin_at_Emin  = i;   }
    if(imin_at_Etarg == -1 && (fabs(fwd_E_cut[i] - fwd_Etarg) < eslSMALLX1)) { ESL_DPRINTF1(("\tAchieved Starg\n"));  imin_at_Etarg = i;   }
    if(             imin_at_Etarg != -1 && imax_at_Etarg == -1 && (fabs(fwd_E_cut[i] - fwd_Etarg) > eslSMALLX1)) { ESL_DPRINTF1(("\tDone with Starg\n")); imax_at_Etarg = i-1; }
    if(i == imax && imin_at_Etarg != -1 && imax_at_Etarg == -1)                                                  { ESL_DPRINTF1(("\tBoundary case, final point is at Starg\n")); imax_at_Etarg = imax; }

    /* Now manipulate our lists for next step */
    /* mark hit cm rank i as 'no longer used' */
    fwd_useme[cmi2fwdi[i]] = FALSE;
    /* if hit cm rank i was the worst scoring fwd hit, update fwd_all to the new worst fwd hit (now that prev worst is no longer used) */
    if(cmi2fwdi[i] == fwd_all) {
      fwd_all++;
      while(!(fwd_useme[fwd_all])) fwd_all++; /* skip all 'no longer used' hits */
    }
    /* if hit cm rank i was within the (1.-F) fraction of worst hits, update fwd_F to new idx in sorted fwd hits that will recognize F fraction of CM hits */
    if(cmi2fwdi[i] <= fwd_F) { /* should we increment fwd_F? */
      if((change_fwd_F[curN]) && ((fwd_F+1) < filN)) { 
	fwd_F++;
	while((!(fwd_useme[fwd_F])) && ((fwd_F+1) < filN)) fwd_F++; /* skip all 'no longer used' hits */
      } 
    }
    curN--;
  }
  ESL_DPRINTF1(("imax_above_Emax: %d\n", imax_above_Emax));
  ESL_DPRINTF1(("imin_at_Etarg:   %d\n", imin_at_Etarg));
  ESL_DPRINTF1(("imax_at_Etarg:   %d\n", imax_at_Etarg));
  ESL_DPRINTF1(("imin_at_Emin:    %d\n", imin_at_Emin));
  ESL_DPRINTF1(("\n"));

#if eslDEBUGLEVEL >= 1
  /* paranoid, expensive check */
  int error_flag = FALSE;
  int nmissed = 0;
  for(i = 0; i <= imax; i++) {
    nmissed = 0;
    for(j = i; j < filN; j++) { 
      if(by_cmA[j].fwd_E > fwd_E_all[i]) { 
	error_flag = 1; 
	ESL_DPRINTF1(("ERROR: i: %d j: %d fwd_E_all[i]: %g < fwd_E[j]: %g\n", i, j, fwd_E_all[i], by_cmA[i].fwd_E));
      }
      if(by_cmA[j].fwd_E > fwd_E_F[i]) { 
	nmissed++;
      }
    }
    if(((float) nmissed/((float) filN-i)) > F) { 
      error_flag = 1;
      ESL_DPRINTF1(("ERROR: i: %d nmissed: %d Fmissed: %f F: %f\n", i, nmissed, ((float) nmissed/((float)filN-i)), F));
    }
  }
  if(error_flag) cm_Fail("Implementation error dude.");
  ESL_DPRINTF1(("Passed expensive paranoia check\n"));
#endif

  /* step through all cutoffs, determining and keeping a representative set */
  int ip_min = imax_above_Emax == -1 ? 0    : imax_above_Emax+1;
  int ip_max = imin_at_Emin    == -1 ? imax : imin_at_Emin;
  int ip = ip_min;
  float max_next_E_cut;
  float prev_E_cut;
  int *saveme;
  int keep_going;
  ESL_ALLOC(saveme, sizeof(int) * (ip_max + 1));
  esl_vec_ISet(saveme, (ip_max+1), FALSE);
  i = 0;
  ESL_DPRINTF1(("ip_min: %d\nip_max: %d\n", ip_min, ip_max));
  while(ip <= ip_max) {
    saveme[ip] = TRUE;
    max_next_E_cut    = fwd_E_cut[ip] * 0.9; /* we'll skip all subsequenct E-value cutoffs that are within 10% of the one we've just added */
    prev_E_cut        = fwd_E_cut[ip];

    ESL_DPRINTF1(("i: %5d fwd_E_cut[ip: %5d]: %12g CM_cut: %12g S: %12g max_next: %12g ", i, ip, fwd_E_cut[ip], by_cmA[ip].cm_E, fwd_E_S[ip], max_next_E_cut)); 
    if     (fabs(fwd_E_cut[ip] - fwd_Emin)  < eslSMALLX1) ESL_DPRINTF1((" [Emin]\n"));
    else if(fabs(fwd_E_cut[ip] - fwd_Etarg) < eslSMALLX1) ESL_DPRINTF1((" [Etarg]\n"));
    else if(fabs(fwd_E_cut[ip] - fwd_Emax)  < eslSMALLX1) ESL_DPRINTF1((" [Emax]\n"));
    else ESL_DPRINTF1(("\n"));

    i++;
    ip++;

    keep_going = TRUE;
    while(keep_going) { 
      if(ip == imin_at_Emin)                  keep_going = FALSE; /* min i for which we acheive Smin,  we want to add this point no matter what */
      else if(ip == imin_at_Etarg)            keep_going = FALSE; /* min i for which we acheive Starg, we want to add this point no matter what */
      else if(ip >  imax)                     keep_going = FALSE; /* we're done with the list */
      else if(fwd_E_cut[ip] < max_next_E_cut) keep_going = FALSE; /* the difference between this E-value cut, and the prev one we added is > 10% */ 
      if(keep_going) ip++; /* skip all points that fail 4 ifs() above */
    }
  }

  float *fwd_E_cut2save;
  float *cm_E_cut2save;
  int n2save = i;
  ESL_ALLOC(fwd_E_cut2save, sizeof(float) * n2save);
  ESL_ALLOC(cm_E_cut2save,  sizeof(float) * n2save);
  i = 0;
  for(ip = ip_min; ip <= ip_max; ip++) {
    if(saveme[ip]) { 
      fwd_E_cut2save[i] = fwd_E_cut[ip];
      cm_E_cut2save[i]  = by_cmA[ip].cm_E;
      i++;
    }
  }
  assert(i == n2save);
  always_better_than_Smax = (imax_above_Emax == -1) ? TRUE : FALSE;
  if((status = SetHMMFilterInfoHMM(hfi, errbuf, F, filN, dbsize, n2save, cm_E_cut2save, fwd_E_cut2save, always_better_than_Smax)) != eslOK) return status;

  /* if --fil-dfile option enabled, print bit scores and cutoffs to the file */
  if(cfg->fildfp != NULL) { 
    assert(cfg->cmstatsA[cmi]->np == 1);
    int ncut = i;
    fprintf(cfg->fildfp, "# Printing cmcalibrate HMM filter threshold determination data for CM: %s\n", cm->name);
    float fwd_bitmax, fwd_bittarg, fwd_bitmin;
    if((status = E2ScoreGivenExpInfo(cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][0], errbuf, fwd_Emax,  &fwd_bitmax))  != eslOK)  return status;
    if((status = E2ScoreGivenExpInfo(cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][0], errbuf, fwd_Etarg, &fwd_bittarg)) != eslOK)  return status;
    if((status = E2ScoreGivenExpInfo(cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][0], errbuf, fwd_Emin,  &fwd_bitmin))  != eslOK)  return status;
    fprintf(cfg->fildfp, "# Max    bit score: %f\n", fwd_bitmax);
    fprintf(cfg->fildfp, "# Target bit score: %f\n", fwd_bittarg);
    fprintf(cfg->fildfp, "# Min    bit score: %f\n", fwd_bitmin);
    fprintf(cfg->fildfp, "# Number of scores: %d\n", filN);
    fprintf(cfg->fildfp, "# HMM/CM bit scores for each sampled sequence are listed next\n");
    fprintf(cfg->fildfp, "# Format of following %d lines: <x> <y>\n", filN);
    fprintf(cfg->fildfp, "# Note: points at the beginning may appear unsorted because they all have\n#       the worst possible E-value, which was used to sort\n");
    if     (hmm_fwd_mode == EXP_CP9_LF) fprintf(cfg->fildfp, "# <x>: HMM local Forward bit scores\n");
    else if(hmm_fwd_mode == EXP_CP9_GF) fprintf(cfg->fildfp, "# <x>: HMM glocal Forward bit scores\n");
    if     (cm_mode == EXP_CM_LC)       fprintf(cfg->fildfp, "# <y>: CM local CYK bit scores\n");
    else if(cm_mode == EXP_CM_LI)       fprintf(cfg->fildfp, "# <y>: CM local Inside bit scores\n");
    else if(cm_mode == EXP_CM_GC)       fprintf(cfg->fildfp, "# <y>: CM glocal CYK bit scores\n");
    else if(cm_mode == EXP_CM_GI)       fprintf(cfg->fildfp, "# <y>: CM glocal Inside bit scores\n");
    for(i = 0; i < filN; i++) fprintf(cfg->fildfp, "%f\t%f\n", fwd_scA[by_fwdA[i].i], cm_scA[by_fwdA[i].i]);
    fprintf(cfg->fildfp, "&\n");

    fprintf(cfg->fildfp, "# HMM/CM cut points saved to CM file are listed next\n");
    fprintf(cfg->fildfp, "# Format of following %d lines: <x> <y>\n", ncut);
    float fwd_bitcut, cm_bitcut;
    if     (hmm_fwd_mode == EXP_CP9_LF) fprintf(cfg->fildfp, "# <x>: Cut point HMM local Forward bit score cutoff\n");
    else if(hmm_fwd_mode == EXP_CP9_GF) fprintf(cfg->fildfp, "# <x>: Cut point HMM glocal Forward bit score cutoff\n");
    if     (cm_mode == EXP_CM_LC)       fprintf(cfg->fildfp, "# <y>: CM local CYK bit score cutoff\n");
    else if(cm_mode == EXP_CM_LI)       fprintf(cfg->fildfp, "# <y>: CM local Inside bit score cutoff\n");
    else if(cm_mode == EXP_CM_GC)       fprintf(cfg->fildfp, "# <y>: CM glocal CYK bit score cutoff\n");
    else if(cm_mode == EXP_CM_GI)       fprintf(cfg->fildfp, "# <y>: CM glocal Inside bit score cutoff\n");
    for(i = 0; i < ncut; i++) { 
      if((status = E2ScoreGivenExpInfo(cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][0], errbuf, fwd_E_cut2save[i], &fwd_bitcut)) != eslOK)  return status;
      if((status = E2ScoreGivenExpInfo(cfg->cmstatsA[cmi]->expAA[cm_mode][0],      errbuf, cm_E_cut2save[i],  &cm_bitcut))  != eslOK)  return status;
      fprintf(cfg->fildfp, "%f\t%f\n", fwd_bitcut, cm_bitcut);
    }
    fprintf(cfg->fildfp, "&\n");

    fprintf(cfg->fildfp, "# All HMM/CM bit score cutoffs listed next\n");
    fprintf(cfg->fildfp, "# Format of following %d lines: <x> <y>\n", imax+1);
    if     (hmm_fwd_mode == EXP_CP9_LF) fprintf(cfg->fildfp, "# <x>: HMM local Forward bit score cutoff\n");
    else if(hmm_fwd_mode == EXP_CP9_GF) fprintf(cfg->fildfp, "# <x>: HMM glocal Forward bit score cutoff\n");
    if     (cm_mode == EXP_CM_LC)       fprintf(cfg->fildfp, "# <y>: CM local CYK bit score cutoff\n");
    else if(cm_mode == EXP_CM_LI)       fprintf(cfg->fildfp, "# <y>: CM local Inside bit score cutoff\n");
    else if(cm_mode == EXP_CM_GC)       fprintf(cfg->fildfp, "# <y>: CM glocal CYK bit score cutoff\n");
    else if(cm_mode == EXP_CM_GI)       fprintf(cfg->fildfp, "# <y>: CM glocal Inside bit score cutoff\n");
    for(i = 0; i <= imax; i++) { 
      if((status = E2ScoreGivenExpInfo(cfg->cmstatsA[cmi]->expAA[hmm_fwd_mode][0], errbuf, fwd_E_cut[i],  &fwd_bitcut))  != eslOK)  return status;
      fprintf(cfg->fildfp, "%f\t%f\n", fwd_bitcut, cm_scA[by_cmA[i].i]);
    }
    fprintf(cfg->fildfp, "&\n");
  }

  free(by_cmA);
  free(by_fwdA);
  free(cmi2fwdi);
  free(fwd_useme);
  free(change_fwd_F);
  free(fwd_E_F);
  free(fwd_E_all);
  free(fwd_E_cut);
  free(fwd_E_S);
  free(saveme);
  free(fwd_E_cut2save);
  free(cm_E_cut2save);
  return eslOK;

 ERROR:
  return status; 
}


/* Function: print_run_info
 * Date:     EPN, Sun Mar  2 16:57:25 2008
 *
 * Purpose:  Print information on this run of cmcalibrate.
 *           Command used to run it, and execution date.
 *
 * Returns:  eslOK on success
 */
static int
print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  int status;
  char *command;
  char *date;

  if((status = get_command(go, errbuf, &command)) != eslOK) return status;
  if((status = GetDate    (errbuf, &date))        != eslOK) return status;

  fprintf(stdout, "%-10s %s\n",  "# command:", command);
  fprintf(stdout, "%-10s %s\n",  "# date:",    date);
  fprintf(stdout, "%-10s %ld\n", "# seed:", esl_randomness_GetSeed(cfg->r));
  if(cfg->nproc > 1) fprintf(stdout, "# %-8s %d\n", "nproc:", cfg->nproc);

  free(command);
  free(date);
  return eslOK;
}

/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call cmscore
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

/* Function: print_post_calibration_info
 * Date:     EPN, Wed Mar  5 05:25:02 2008
 * Purpose:  Print info about calibration for a CM we just calibrated including
 *           timings. 
 */
int print_post_calibration_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, double **exp_psecAA, double *fil_psecA, double **exp_asecAA, double *fil_asecA)
{
  char  time_buf[128];	      /* for printing run time */
  int   exp_mode;             /* counter over exp tail modes */
  double total_psec = 0.;     /* predicted number of seconds for cm, all stages */
  double total_asec = 0.;     /* actual number of seconds for cm, all stages */
  int   expN;                 /* nseq we'll calibrate on */
  int   p, ps, pe;            /* partition vars */
  float L_Mb;                 /* total seq length we'll calibrate exp tails on in Mb */
  ExpInfo_t *exp;             /* pointer to current exp tail info, for convenience */

  if(exp_psecAA == NULL || fil_psecA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "print_post_calibration_info, exp_psecAA or fil_psecA is NULL");
  if(exp_asecAA == NULL || fil_asecA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "print_post_calibration_info, exp_asecAA or fil_asecA is NULL");

  fprintf(fp, "#\n");
  fprintf(fp, "# Post-calibration info for CM %d: %s\n", cfg->ncm, cm->name);
  fprintf(fp, "#\n");
  fprintf(fp, "# Exponential tail fitting:\n");
  fprintf(fp, "#\n");
  if(cfg->np != 1) { /* --exp-pfile invoked */
    fprintf(fp, "# %-3s  %3s  %3s %4s %3s %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",    "",     "",   "",        "",   "",       "",        "",       "     running time      ");
    fprintf(fp, "# %-3s  %3s  %3s %4s %3s %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",    "",     "",   "",        "",   "",       "",        "",       "-----------------------");
    fprintf(fp, "# %-3s  %3s  %3s %4s %3s %3s %7s %6s %6s %7s %10s %10s\n","mod", "cfg", "alg", "part", "ps", "pe",   "L (Mb)",  "mu",     "lambda",   "nhits", "predicted",     "actual");
    fprintf(fp, "# %3s  %3s  %3s %4s %3s %3s %7s %6s %6s %7s %10s %10s\n", "---", "---", "---", "----", "---", "---", "-------", "------", "------", "-------", "----------", "----------");
  }
  else { 
    fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",  "",   "",       "",        "",       "     running time      ");
    fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",  "",   "",       "",        "",       "-----------------------");
    fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %10s %10s\n","mod", "cfg", "alg", "L (Mb)",  "mu",     "lambda",   "nhits", "predicted",     "actual");
    fprintf(fp, "# %3s  %3s  %3s %7s %6s %6s %7s %10s %10s\n", "---", "---", "---", "-------", "------", "------", "-------", "----------", "----------");
  }
  for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
    if(ExpModeIsLocal(exp_mode)) { expN = ExpModeIsForCM(exp_mode) ? cfg->exp_cmN_loc : cfg->exp_hmmN_loc; }
    else                         { expN = ExpModeIsForCM(exp_mode) ? cfg->exp_cmN_glc : cfg->exp_hmmN_glc; }
    L_Mb = ((float) expN * (float) cfg->expL) / 1000000.;
    for (p = 0; p < cfg->np; p++) {
      total_psec += exp_psecAA[exp_mode][p];
      total_asec += exp_asecAA[exp_mode][p];
      ps = cfg->pstart[p];
      pe = (p == (cfg->np-1)) ? 100 : cfg->pstart[p+1]-1;
      FormatTimeString(time_buf, exp_psecAA[exp_mode][p], FALSE);
      exp = cfg->cmstatsA[cfg->ncm-1]->expAA[exp_mode][p];
      if(cfg->np != 1) fprintf(fp, "  %-12s %4d %3d %3d %7.2f %6.2f %6.3f %7d %10s", DescribeExpMode(exp_mode), p+1, ps, pe, L_Mb, exp->mu_orig, exp->lambda, exp->nrandhits, time_buf);
      else             fprintf(fp, "  %-12s %7.2f %6.2f %6.3f %7d %10s", DescribeExpMode(exp_mode), L_Mb, exp->mu_orig, exp->lambda, exp->nrandhits, time_buf);
      FormatTimeString(time_buf, exp_asecAA[exp_mode][p], FALSE);
      fprintf(fp, " %10s\n", time_buf);
    }
  }
  /* print filter threshold stats */
  fprintf(fp, "#\n");
  fprintf(fp, "# HMM filter threshold determination:\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# %3s  %6s  %22s\n", "",    "",             "     running time     ");
  fprintf(fp, "# %3s  %6s  %22s\n", "",    "",             "----------------------");
  fprintf(fp, "# %3s  %6s  %10s  %10s\n", "cfg",   "nseq", "predicted",  "actual");
  fprintf(fp, "# %3s  %6s  %10s  %10s\n", "---", "------", "----------", "----------");

  for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
    if(exp_mode == EXP_CM_GI || exp_mode == EXP_CM_LI) { /* CM Inside mode, only time we do filter threshold calculations, we'll fill in CYK AND Inside thresholds */
      total_psec += fil_psecA[exp_mode];
      total_asec += fil_asecA[exp_mode];
      FormatTimeString(time_buf, fil_psecA[exp_mode], FALSE);
      fprintf(fp, "  %3s  %6d  %10s", ((exp_mode == EXP_CM_GI) ? "glc" : "loc"), esl_opt_GetInteger(go, "--fil-N"), time_buf);
      FormatTimeString(time_buf, fil_asecA[exp_mode], FALSE);
      fprintf(fp, "  %10s\n", time_buf);
    }
  }
  fprintf(fp, "#\n");
  FormatTimeString(time_buf, total_psec, FALSE);
  fprintf(fp, "# total predicted time for CM: %s\n", time_buf);
  FormatTimeString(time_buf, total_asec, FALSE);
  fprintf(fp, "# total actual    time for CM: %s\n", time_buf);
  
  fflush(fp);
  return eslOK;
}
  
/* Function: estimate_time_for_exp_round
 * Date:     EPN, Wed Mar  5 05:46:45 2008
 * Purpose:  Estimate search time for round of exp tail fitting
 *           of exp mode <exp_mode>. This is done by actually 
 *           searching a sequence with the appropriate algorithm. 
 *           The length of the sequence to search is set such 
 *           that it should take about <targ_sec> seconds.
 */
int estimate_time_for_exp_round(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_res)
{
  int    status;
  int    L;                /* length of sequence we'll generate and search to get time estimate */
  double psec_per_Mc;      /* rough prediction at seconds per Mc based on empirical run times I've witnessed, 
                            * doesn't need to be very accurate as we just use it to set length of seq to search to get real prediction */
  float  Mc;               /* millions of DP calculations we're going to do */
  float  Mc_per_res;       /* millions of dp calcs per residue, if searching with CM, corrects for first W residues requiring less dp calcs */
  int    irrelevant_W;     /* temporary W */
  int    orig_search_opts; /* cm->search_opts when function was entered */
  float  sec_per_res;      /* seconds per residue */
  float  targ_sec = 0.1;   /* target number of seconds our timing expt will take */
  int    Lmin = 100;       /* minimum number of residues to search to get timing */
  int    Lmax = 10000;     /* maximum number of residues to search to get timing */
  int    use_qdb;          /* TRUE if we're using QDB, FALSE if not */
  double *dnull = NULL;    /* background distro for generating random seqs */
  int     i;               /* counter */
  ESL_DSQ *dsq;            /* the random seq we'll create and search to get predicted time */
  ESL_STOPWATCH *w  = esl_stopwatch_Create(); /* for timings */
  
  if(w == NULL)               ESL_FAIL(status, errbuf, "estimate_time_for_exp_round(): memory error, stopwatch not created.\n");
  if(ret_sec_per_res == NULL) ESL_FAIL(status, errbuf, "estimate_time_for_exp_round(): ret_sec_per_res is NULL");

  /* update search info for round 0 (final round) for exp tail mode */
  UpdateSearchInfoForExpMode(cm, 0, exp_mode);
  orig_search_opts = cm->search_opts;
  cm->search_opts  = cm->si->search_opts[0]; /* we'll restore cm->search_opts to orig_search_opts at end of the function */

  if(ExpModeIsForCM(exp_mode)) { 
    if(cm->smx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "estimate_time_for_exp_round(), cm->smx is NULL");
    use_qdb  = (cm->smx->dmin == NULL && cm->smx->dmax == NULL) ? FALSE : TRUE;
    if(use_qdb) { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, cm->smx->beta_qdb, &Mc_per_res, &irrelevant_W))  != eslOK) return status; }
    else        { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, TRUE,  cm->beta_W,        &Mc_per_res, &irrelevant_W))  != eslOK) return status; }
    psec_per_Mc = (cm->search_opts & CM_SEARCH_INSIDE) ? (1. /  75.) : (1. / 275.);  /*  75 Mc/S inside;  275 Mc/S CYK */
    /* determine L that will take about <targ_sec> seconds */
    L = targ_sec / (psec_per_Mc * Mc_per_res);
    L = ESL_MAX(L, Lmin); /* we have to search at least <Lmin> residues */
    L = ESL_MIN(L, Lmax); /* we want to search at most  <Lmax> residues */
    /* now determine exactly how many dp calculations we'd do if we search L residues, 
     * this won't be the same as Mc_per_res * L b/c Mc_per_res from cm_GetNCalcsPerResidueGivenBeta
     * b/c that was calculated after correcting for the fact that the first W residues have fewer
     * DP calcs than all other residues, b/c d < W for first W residues.
     */  
    if((status = cm_CountSearchDPCalcs(cm, errbuf, L, cm->smx->dmin, cm->smx->dmax, cm->smx->W, FALSE,  NULL, &Mc)) != eslOK) return status;
    /* FALSE says don't correct for fewer dp calcs for first W residues, we want to know how many total DP calcs
     * there will be in L residues */
    Mc *= L; /* Mc was for 1 residue, multiply by L to get Mc for L residues */
  }
  else { /* HMM mode */
    if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &Mc_per_res)) != eslOK) return status;
    psec_per_Mc = (cm->search_opts & CM_SEARCH_HMMFORWARD) ? (1. / 175.) : (1. / 380.);  /* 175 Mc/S forward; 380 Mc/S viterbi */
    /* determine L that will take about <targ_sec. seconds */
    L  = targ_sec / (psec_per_Mc * Mc_per_res);
    L  = ESL_MAX(L, Lmin); /* we have to search at least <Lmin> residues */
    L  = ESL_MIN(L, Lmax); /* we want to search at most  <Lmax> residues */
    /* how many millions of DP cells will it be? */
    Mc = Mc_per_res * L;
  }

  /* create dnull for generating seqs */
  if(esl_opt_IsDefault(go, "--exp-gc")) { /* only setup dnull if --exp-gc NOT enabled */
    ESL_ALLOC(dnull, sizeof(double) * cm->abc->K);
    for(i = 0; i < cm->abc->K; i++) dnull[i] = (double) cm->null[i];
    esl_vec_DNorm(dnull, cm->abc->K);    
  }

  search_results_t *results;
  esl_stopwatch_Start(w);
  /* simulate a workunit, generate a sequence, search it, and remove overlaps */
  /*printf("exptL: %d\n", L);*/
  if(esl_opt_GetBoolean(go, "--exp-random")) { 
    if((status = get_random_dsq(cfg, errbuf, cm, dnull, L, &dsq)) != eslOK) return status;
  }
  else { 
    if((status = get_genomic_sequence_from_hmm(cfg, errbuf, cm, L, &dsq)) != eslOK) return status;
  }
  if((status = ProcessSearchWorkunit (cm,  errbuf, dsq, L, &results, esl_opt_GetReal(go, "--mxsize"), 0, NULL, NULL)) != eslOK) return status;
  RemoveOverlappingHits(results, 1, L);

  esl_stopwatch_Stop(w);
  if(w != NULL) esl_stopwatch_Destroy(w);
  FreeResults(results);
  if(dnull != NULL) free(dnull);
  free(dsq);

  cm->search_opts = orig_search_opts;
  sec_per_res = w->user * (Mc_per_res / Mc);
  ESL_DPRINTF1(("L: %d\n", L));
  ESL_DPRINTF1(("w->user: %f\n", w->user));
  ESL_DPRINTF1(("sec_per_res: %f\n", sec_per_res));
  ESL_DPRINTF1(("Mc_per_res: %f\n", Mc_per_res));
  ESL_DPRINTF1(("Mc: %f\n", Mc));

  /*printf("L: %d\n", L);
    printf("w->user: %f\n", w->user);
    printf("sec_per_res: %f\n", sec_per_res);
    printf("Mc_per_res: %f\n", Mc_per_res);
    printf("Mc: %f\n", Mc);*/
  
  *ret_sec_per_res = sec_per_res;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "estimate_time_for_exp_round(): memory error.\n");
  return status; /* NEVERREACHED */
}

/* Function: estimate_time_for_fil_round
 * Date:     EPN, Wed Mar  5 05:46:45 2008
 * Purpose:  Estimate search time for round of filter threshold calculation
 *           with CM in exp mode <exp_mode>. 
 */
int estimate_time_for_fil_round(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_seq)
{
  int    status;
  int    orig_search_opts; /* cm->search_opts when function was entered */
  float  sec_per_seq;      /* seconds per residue */
  int    targ_len = 1000;
  int    nseq;             /* number of sequences to sample and search for our expt */
  /* score arrays for our nseq sample seqs, irrelevant really, only nec to avoid contract violations to process_filter_workunit() */
  float   *tmp_cyk_scA = NULL;    
  float   *tmp_ins_scA = NULL;    
  float   *tmp_fwd_scA = NULL;    
  int     *tmp_partA   = NULL;    
  int      min_nseq_global = 10;
  int      min_nseq_local  = 25;

  nseq = (int) (((float) targ_len / (float) cm->clen) + 0.5); 
  /* we have to sample at least <min_nseq_{global,local} sequences */
  if(ExpModeIsLocal(exp_mode)) nseq = ESL_MAX(nseq, min_nseq_local); 
  else                         nseq = ESL_MAX(nseq, min_nseq_global); 

  ESL_STOPWATCH *w  = esl_stopwatch_Create();

  /* update search info for round 0 (final round) for current mode */
  UpdateSearchInfoForExpMode(cm, 0, exp_mode);
  orig_search_opts = cm->search_opts;
  cm->search_opts  = cm->si->search_opts[0]; /* we'll restore cm->search_opts to orig_search_opts at end of the function */

  if(w == NULL)               ESL_FAIL(status, errbuf, "estimate_time_for_fil_round(): memory error, stopwatch not created.\n");
  if(ret_sec_per_seq == NULL) ESL_FAIL(status, errbuf, "estimate_time_for_fil_round(): ret_sec_per_res is NULL");

  esl_stopwatch_Start(w);
  if((status = process_filter_workunit (go, cfg, errbuf, cm, nseq, &tmp_cyk_scA, &tmp_ins_scA, &tmp_fwd_scA, &tmp_partA)) != eslOK) return status;
  esl_stopwatch_Stop(w);

  free(tmp_cyk_scA);
  free(tmp_ins_scA);
  free(tmp_fwd_scA);
  free(tmp_partA);

  cm->search_opts = orig_search_opts;
  sec_per_seq = w->user / (float) nseq;

  ESL_DPRINTF1(("nseq: %d\n", nseq));
  ESL_DPRINTF1(("w->user: %f\n", w->user));
  ESL_DPRINTF1(("sec_per_seq: %f\n", sec_per_seq));

  esl_stopwatch_Destroy(w);
  *ret_sec_per_seq = sec_per_seq;
  return eslOK;
}


/* Function: print_per_cm_column_headings
 * Date:     EPN, Tue Jan  8 05:51:47 2008
 *
 * Purpose:  Print per-CM info to stdout. 
 *
 * Returns:  eslOK on success
 */
static int
print_per_cm_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  printf("#\n");
  if(esl_opt_IsDefault(go, "--forecast")) { 
    printf("# Calibrating CM %d: %s\n", cfg->ncm, cm->name);
    printf("#\n");
    if(cfg->np != 1) { /* --exp-pfile invoked */
      printf("# %8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %22s\n",           "",      "",    "",    "",      "",   "",    "",           "",       "", "     running time    ");
      printf("# %8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %22s\n",           "",      "",    "",    "",      "",   "",    "",           "",       "", "----------------------");
      printf("# %-8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %10s  %10s\n", "stage",    "mod", "cfg", "alg", "part", "ps",  "pe",  "expL (Mb)", "filN",   "predicted", "actual");
      printf("# %8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %10s  %10s\n",  "--------", "---", "---", "---", "----", "---", "---", "---------", "------", "----------", "----------");
    }
    else { /* no partitions */
      printf("# %8s  %3s  %3s  %3s  %9s %6s %22s\n",           "",      "",    "",    "",    "",       "", "     running time    ");
      printf("# %8s  %3s  %3s  %3s  %9s %6s %22s\n",           "",      "",    "",    "",    "",       "", "----------------------");
      printf("# %-8s  %3s  %3s  %3s  %9s %6s %10s  %10s\n", "stage",    "mod", "cfg", "alg", "expL (Mb)", "filN",   "predicted", "actual");
      printf("# %8s  %3s  %3s  %3s  %9s %6s %10s  %10s\n",  "--------", "---", "---", "---", "---------", "------", "----------", "----------");
    }
  }
  else { 
    if(cfg->np != 1) { /* --exp-pfile invoked */
      printf("# Forecasting time for %d processor(s) to calibrate CM %d: %s\n", esl_opt_GetInteger(go, "--forecast"), cfg->ncm, cm->name);
      printf("#\n");
      printf("# %-8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %14s\n", "stage",    "mod", "cfg", "alg", "part", "ps",  "pe",  "expL (Mb)", "filN",   "predicted time");
      printf("# %8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %14s\n",  "--------", "---", "---", "---", "----", "---", "---", "---------", "------", "--------------");
    }
    else { /* no partitions */
      printf("# Forecasting time for %d processor(s) to calibrate CM %d: %s\n", esl_opt_GetInteger(go, "--forecast"), cfg->ncm, cm->name);
      printf("#\n");
      printf("# %-8s  %3s  %3s  %3s  %9s %6s %14s\n", "stage",    "mod", "cfg", "alg", "expL (Mb)", "filN",   "predicted time");
      printf("# %8s  %3s  %3s  %3s  %9s %6s %14s\n",  "--------", "---", "---", "---", "---------", "------", "--------------");
    }
  }
  return eslOK;
}


/* Function: print_per_cm_summary
 * Date:     EPN, Thu Mar  6 10:35:49 2008
 *
 * Purpose:  Print per-CM summary to stdout. 
 *
 * Returns:  eslOK on success
 */
static int
print_per_cm_summary(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, double psec, double asec)
{
  char  time_buf[128];	      /* for printing run time */
  if(esl_opt_IsDefault(go, "--forecast")) { 
    if(cfg->np != 1) { /* --exp-pfile invoked */
      FormatTimeString(time_buf, psec, FALSE);
      printf("# %8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %10s  %10s\n", "--------", "---", "---", "---", "----", "---", "---", "---------", "-----", "----------", "----------");
      printf("# %-8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %10s",       "all",        "-",   "-",   "-",    "-",   "-",   "-",         "-",     "-",   time_buf);
      FormatTimeString(time_buf, asec, FALSE);
      printf("  %10s\n", time_buf);
    }
    else { /* no partitions */
      FormatTimeString(time_buf, psec, FALSE);
      printf("# %8s  %3s  %3s  %3s  %9s %6s %10s  %10s\n", "--------", "---", "---", "---", "---------", "-----", "----------", "----------");
      printf("# %-8s  %3s  %3s  %3s  %9s %6s %10s",       "all",        "-",   "-",   "-",         "-",     "-",   time_buf);
      FormatTimeString(time_buf, asec, FALSE);
      printf("  %10s\n", time_buf);
    }
  }
  else { 
    if(cfg->np != 1) { /* --exp-pfile invoked */
      FormatTimeString(time_buf, psec, FALSE);
      printf("# %8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %14s\n", "--------", "---", "---", "---",  "----", "---", "---", "---------", "-----", "--------------");
      printf("# %-8s  %3s  %3s  %3s  %4s %3s %3s %9s %6s %14s\n", "all",       "-",   "-",   "-",     "-",   "-",   "-",         "-",     "-",         time_buf);
    }
    else { /* no partitions */
      FormatTimeString(time_buf, psec, FALSE);
      printf("# %8s  %3s  %3s  %3s  %9s %6s %14s\n", "--------", "---", "---", "---",  "---------", "-----", "--------------");
      printf("# %-8s  %3s  %3s  %3s  %9s %6s %14s\n", "all",        "-",   "-",   "-",         "-",     "-",         time_buf);
    }
  }
  return eslOK;
}

/* Function: print_exp_line
 * Date:     EPN, Thu Mar  6 13:38:01 2008
 *
 * Purpose:  Print a line describing exp tail fitting for a given mode.
 */
static int
print_exp_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, int expN, int expL, int p, double psec)
{
  char  time_buf[128];	      /* for printing run time */
  float expL_Mb;              /* exp seq length in Mb */

  int ps;
  int pe;

  ps = cfg->pstart[p];
  pe = (p == (cfg->np-1)) ? 100 : cfg->pstart[p+1]-1;

  FormatTimeString(time_buf, psec, FALSE);

  expL_Mb =  (float) expN * (float) expL; 
  expL_Mb /= 1000000.;

  if(!esl_opt_IsDefault(go, "--forecast")) { 
    if(cfg->np != 1) printf("  %-8s  %-12s  %4d %3d %3d %9.2f %6s %14s\n",  "exp tail", DescribeExpMode(exp_mode), p+1, ps, pe, expL_Mb, "-", time_buf);
    else             printf("  %-8s  %-12s  %9.2f %6s %14s\n",              "exp tail", DescribeExpMode(exp_mode), expL_Mb, "-", time_buf);
  }
  else { 
    if(cfg->np != 1) printf("  %-8s  %-12s  %4d %3d %3d %9.2f %6s %10s",  "exp tail", DescribeExpMode(exp_mode), p+1, ps, pe, expL_Mb, "-", time_buf);
    else             printf("  %-8s  %-12s  %9.2f %6s %10s",              "exp tail", DescribeExpMode(exp_mode), expL_Mb, "-", time_buf);
  }
  fflush(stdout);
  return eslOK;
}

/* Function: print_fil_line
 * Date:     EPN, Thu Mar  6 13:38:01 2008
 *
 * Purpose:  Print a line describing exp tail fitting for a given mode.
 */
static int
print_fil_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, double psec)
{
  char  time_buf[128];	      /* for printing run time */
  int   filN = esl_opt_GetInteger(go, "--fil-N"); /* number of sequences to search for filter threshold calculation */

  FormatTimeString(time_buf, psec, FALSE);
  if(!esl_opt_IsDefault(go, "--forecast")) { 
    if(cfg->np != 1) printf("  %-8s  %3s  %3s  %3s  %4s %3s %3s %9s %6d %14s\n",  "filter", "-", ((exp_mode == EXP_CM_GI) ? "glc" : "loc"), "-", "-", "-", "-", "-", filN, time_buf);
    else             printf("  %-8s  %3s  %3s  %3s  %9s %6d %14s\n",              "filter", "-", ((exp_mode == EXP_CM_GI) ? "glc" : "loc"), "-", "-", filN, time_buf);
  }
  else { 
    if(cfg->np != 1) printf("  %-8s  %3s  %3s  %3s  %4s %3s %3s %9s %6d %10s",  "filter", "-", ((exp_mode == EXP_CM_GI) ? "glc" : "loc"), "-", "-", "-", "-", "-", filN, time_buf);
    else             printf("  %-8s  %3s  %3s  %3s  %9s %6d %10s",              "filter", "-", ((exp_mode == EXP_CM_GI) ? "glc" : "loc"), "-", "-", filN, time_buf);
  }
  fflush(stdout);
  return eslOK;
}
  
/* Function: update_hmm_exp_length
 * Date:     EPN, Fri Mar  7 05:12:05 2008
 * Purpose:  Potentially reset the sequence length used for fitting
 *           HMM exp tails based on the ratio of HMM to CM DP calculations
 *           such that the number of calculations with the HMM is
 *           at least <esl_opt_GetReal(go, "--exp-fract")> 
 *           the fraction of CM calculations.
 *           Based on the idea that we're willing to spend at least
 *           <esl_opt_GetReal(go, "--exp-fract")> the time we spend 
 *           fitting the CM exp tails on the HMM. <min_frac> is usually 
 *           small, it's 0.1 by default.
 */
int update_hmm_exp_length(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int    status;
  float  cm_Mc_per_res;    /* millions of dp calcs per residue searching with CM */
  float  hmm_Mc_per_res;   /* millions of dp calcs per residue searching with HMM */
  int    irrelevant_W;
  float  hmmL_loc, hmmL_glc;
  int    hmmL_loc_int, hmmL_glc_int;
  /* update search info for round 0 (final round) for exp tail mode */

  if(cm->smx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_hmm_exp_length(), cm->smx is NULL");
  /* estimate millions of CM calcs per residue */
  int use_qdb  = (cm->smx->dmin == NULL && cm->smx->dmax == NULL) ? FALSE : TRUE;
  if(use_qdb) { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, cm->smx->beta_qdb, &cm_Mc_per_res, &irrelevant_W))  != eslOK) return status; }
  else        { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, TRUE,  cm->beta_W,        &cm_Mc_per_res, &irrelevant_W))  != eslOK) return status; }
  /* estimate millions of HMM calcs per residue */
  if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &hmm_Mc_per_res)) != eslOK) return status;

  /* when calibrating the HMM, the calibration step will do a cp9_ViterbiBackward or cp9_Backward()
   * run for ALL hits the HMM finds, this is to get actual score of a single i..j parse, 
   * theoretically this could double number of calcs for HMM calibrations.
   * so we multiply number of calcs by 2 here.
   */
  hmm_Mc_per_res *= 2.;

  /* set seq length to search for local HMM */
  hmmL_loc = (cm_Mc_per_res / hmm_Mc_per_res) * esl_opt_GetReal(go, "--exp-fract") * esl_opt_GetReal(go, "--exp-cmL-loc");
  hmmL_loc = ESL_MAX(hmmL_loc, esl_opt_GetReal(go, "--exp-hmmLn-loc"));
  hmmL_loc = ESL_MIN(hmmL_loc, esl_opt_GetReal(go, "--exp-hmmLx"));
  hmmL_loc *= 1000000.; /* convert to Mb */
  hmmL_loc_int = (int) hmmL_loc; 
  cfg->exp_hmmN_loc = (int) (((float) hmmL_loc_int / (float) cfg->expL) + 0.999999); 

  /* set seq length to search for gglcal HMM */
  hmmL_glc = (cm_Mc_per_res / hmm_Mc_per_res) * esl_opt_GetReal(go, "--exp-fract") * esl_opt_GetReal(go, "--exp-cmL-glc");
  hmmL_glc = ESL_MAX(hmmL_glc, esl_opt_GetReal(go, "--exp-hmmLn-glc"));
  hmmL_glc = ESL_MIN(hmmL_glc, esl_opt_GetReal(go, "--exp-hmmLx"));
  hmmL_glc *= 1000000.; /* convert to Mb */
  hmmL_glc_int = (int) hmmL_glc; 
  cfg->exp_hmmN_glc = (int) (((float) hmmL_glc_int / (float) cfg->expL) + 0.999999); 

  ESL_DPRINTF1(("cm  ncalcs: %f\n", cm_Mc_per_res));
  ESL_DPRINTF1(("hmm ncalcs: %f\n", hmm_Mc_per_res));
  ESL_DPRINTF1(("cm/hmm ratio: %f\n", cm_Mc_per_res / hmm_Mc_per_res));
  ESL_DPRINTF1(("min hmmL_loc: %f\n", esl_opt_GetReal(go, "--exp-hmmLn-loc")));
  ESL_DPRINTF1(("min hmmL_glc: %f\n", esl_opt_GetReal(go, "--exp-hmmLn-glc")));
  ESL_DPRINTF1(("hmmL_loc: %f\n", hmmL_loc));
  ESL_DPRINTF1(("hmmL_loc: %f\n", hmmL_glc));
  ESL_DPRINTF1(("cfg->exp_hmmN_loc: %d\n", cfg->exp_hmmN_loc));
  ESL_DPRINTF1(("cfg->exp_hmmN_glc: %d\n", cfg->exp_hmmN_glc));
  ESL_DPRINTF1(("cfg->exp_hmmL_loc * cfg->exp_hmmN in Mb: %f\n", (float) (cfg->expL * cfg->exp_hmmN_loc) / 1000000.));
  ESL_DPRINTF1(("cfg->exp_hmmL_glc * cfg->exp_hmmN in Mb: %f\n", (float) (cfg->expL * cfg->exp_hmmN_glc) / 1000000.));

  return eslOK;
}

/* Function: get_genomic_sequence_from_hmm
 * Date:     EPN, Tue May 20 17:40:54 2008
 * 
 * Purpose:  Emit sequence from a fully connected 
 *           5 state HMM that was trained by EM 
 *           from 30 Mb of 100 Kb chunks of real 
 *           genomes of hand selected GC contents
 *           (10 Mb each from Archaea, Bacteria,
 *            Eukarya genomes). See 
 *           ~nawrockie/notebook/8_0326_inf_default_gc/ 
 *            for more info. 
 *           There were larger HMMs that 'performed'
 *           better, but this 5 state guy was a good
 *           balance b/t performance and number of 
 *           parameters. Performance was judged by
 *           how similar the generated sequence was
 *           to the training 30 Mb genomic sequence.
 */
int
get_genomic_sequence_from_hmm(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int L, ESL_DSQ **ret_dsq)
{
  int      status;
  ESL_DSQ *dsq = NULL;
  int      nstates = 5;
  int      i, si, x;

  /* contract check, make sure we're in a valid mode */
  if(cm->abc->type != eslRNA && cm->abc->type != eslDNA) ESL_FAIL(eslEINCOMPAT, errbuf, "get_genomic_sequence_from_hmm(), cm->abc is not eslRNA nor eslDNA");

  /* start probabilities */
  double *sA;
  ESL_ALLOC(sA, sizeof(double) * nstates);

  sA[0] = 0.157377049180328;
  sA[1] = 0.39344262295082;
  sA[2] = 0.265573770491803; 
  sA[3] = 0.00327868852459016; 
  sA[4] = 0.180327868852459;
  esl_vec_DNorm(sA, nstates);

  /* transition probabilities */
  double **tAA;
  ESL_ALLOC(tAA, sizeof(double *) * nstates);
  for(i = 0; i < nstates; i ++) ESL_ALLOC(tAA[i], sizeof(double) * nstates);

  tAA[0][0] = 0.999483637183643;
  tAA[0][1] = 0.000317942006440604; 
  tAA[0][2] = 0.000185401071732768; 
  tAA[0][3] = 2.60394763669618e-07; 
  tAA[0][4] = 1.27593434198113e-05;
  esl_vec_DNorm(tAA[0], nstates);

  tAA[1][0] = 9.76333640771184e-05; 
  tAA[1][1] = 0.99980020511745; 
  tAA[1][2] = 9.191359010352e-05; 
  tAA[1][3] = 7.94413051888677e-08; 
  tAA[1][4] = 1.01684870641751e-05;
  esl_vec_DNorm(tAA[1], nstates);

  tAA[2][0] = 1.3223694798182e-07; 
  tAA[2][1] = 0.000155642887774602; 
  tAA[2][2] = 0.999700615549769; 
  tAA[2][3] = 9.15079680034191e-05; 
  tAA[2][4] = 5.21013575048369e-05;
  esl_vec_DNorm(tAA[2], nstates);

  tAA[3][0] = 0.994252873563218; 
  tAA[3][1] = 0.0014367816091954; 
  tAA[3][2] = 0.0014367816091954; 
  tAA[3][3] = 0.0014367816091954; 
  tAA[3][4] = 0.0014367816091954;
  esl_vec_DNorm(tAA[3], nstates);

  tAA[4][0] = 8.32138798088677e-06; 
  tAA[4][1] = 2.16356087503056e-05; 
  tAA[4][2] = 6.42411152124459e-05; 
  tAA[4][3] = 1.66427759617735e-07; 
  tAA[4][4] = 0.999905635460297;
  esl_vec_DNorm(tAA[4], nstates);

  /* emission probabilities */
  double **eAA;
  ESL_ALLOC(eAA, sizeof(double *) * nstates);
  for(i = 0; i < nstates; i ++) ESL_ALLOC(eAA[i], sizeof(double) * cm->abc->K);

  eAA[0][0] = 0.370906566523225;
  eAA[0][1] = 0.129213995153577;
  eAA[0][2] = 0.130511270043053;
  eAA[0][3] = 0.369368168280145;
  esl_vec_DNorm(eAA[0], cm->abc->K);

  eAA[1][0] = 0.305194882571888;
  eAA[1][1] = 0.194580936415687;
  eAA[1][2] = 0.192343972160245;
  eAA[1][3] = 0.307880208852179;
  esl_vec_DNorm(eAA[1], cm->abc->K);

  eAA[2][0] = 0.238484980800698;
  eAA[2][1] = 0.261262845707113;
  eAA[2][2] = 0.261810301531792;
  eAA[2][3] = 0.238441871960397;
  esl_vec_DNorm(eAA[2], cm->abc->K);

  eAA[3][0] = 0.699280575539568;
  eAA[3][1] = 0.00143884892086331;
  eAA[3][2] = 0.00143884892086331;
  eAA[3][3] = 0.297841726618705;
  esl_vec_DNorm(eAA[3], cm->abc->K);

  eAA[4][0] = 0.169064007664923;
  eAA[4][1] = 0.331718611320207;
  eAA[4][2] = 0.33045427183482;
  eAA[4][3] = 0.16876310918005;
  esl_vec_DNorm(eAA[4], cm->abc->K);

  /* generate sequence */
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;

  /* pick initial state to emit from */
  si = esl_rnd_DChoose(cfg->r, sA, nstates);
  for (x = 1; x <= L; x++) {
    dsq[x] = esl_rnd_DChoose(cfg->r, eAA[si], cm->abc->K); /* emit residue */
    si = esl_rnd_DChoose(cfg->r, tAA[si], nstates);        /* make transition */
  }
  dsq[x] = '\0';

  for(i = 0; i < nstates; i++) { 
    free(eAA[i]); 
    free(tAA[i]); 
  }
  free(eAA);
  free(tAA);
  free(sA);

  /* TEMPORARY! */
  /*FILE *fp;
  fp = fopen("hmm.fa", "a");
  ESL_SQ *sq;
  sq = esl_sq_CreateDigitalFrom(cm->abc, "irrelevant", dsq, L, NULL, NULL, NULL);
  esl_sq_Textize(sq);
  esl_sqio_Write(fp, sq, eslSQFILE_FASTA);
  fclose(fp);*/

  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  return status;
}  

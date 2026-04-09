/* cmsim: scoring CMs against simulated sequences.
 * [INCOMPLETE: NOT CURRENTLY COMPILED EPN, Tue Mar 20 05:45:57 2012]
 *
 * Main testbed for exploring the statistical behavior of Infernal
 * scores on random sequences, and importance sampling.
 *
 * EPN, Fri Apr 29 14:00:58 2011
 */
#include "config.h"
#include <esl_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_exponential.h>
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_randomseq.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "hmmer.h"

#include "infernal.h"

#define ALPHOPTS "--rna,--dna"           /* Exclusive options for alphabet choice */
#define OUTOPTS "-u,-c,-a,--ahmm,--shmm" /* Exclusive options for output */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles  reqs  incomp  help  docgroup*/
  { "-h", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",
    1 },
  { "-v", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "verbose: print scores", 1 },
  { "-i", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "do importance sampling", 1 },
  { "-g", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "put CM in glocal mode", 1 },
  { "--rL", eslARG_INT, "10000", NULL, "n>0", NULL, NULL, NULL, "length of random target seqs", 1 },
  { "--rN", eslARG_INT, "10", NULL, "n>0", NULL, NULL, NULL, "number of random target seqs", 1 },
  { "--rhmm", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "generate random sequences from realistic HMM", 1 },
  { "--rtailp", eslARG_REAL, NULL, NULL, "0.0<x<0.6", NULL, NULL, NULL,
    "override: tail fraction to fit for random seqs (default: use hits/Mb)", 1 },
  { "--ltailn", eslARG_INT, "750", NULL, "n>=1", NULL, NULL, NULL,
    "hits/Mb to fit for local mode tail [df: 750, matches cmcalibrate]", 1 },
  { "--gtailn", eslARG_INT, "250", NULL, "n>=1", NULL, NULL, NULL,
    "hits/Mb to fit for glocal mode tail [df: 250, matches cmcalibrate]", 1 },
  { "--iN", eslARG_INT, "1000", NULL, "n>0", NULL, NULL, NULL, "number of sampled target seqs", 1 },
  { "--iT", eslARG_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
    "set min bit sc for hits in sampled seqs to <x>", 1 },
  { "--ilocal", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "allow local begins/ends in sampled target seqs", 1 },
  { "--itailp", eslARG_REAL, NULL, NULL, "0.0<x<=1.0", NULL, NULL, NULL,
    "override: tail fraction to fit for IS seqs (default: use hits/Mb)", 1 },
  { "--inonbanded", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "do not use HMM bands to score sampled sequences", 1 },
  { "--null3", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "use NULL3 correction for random/sampled seqs", 1 },
  { "--beta", eslARG_REAL, "1e-7", NULL, "0<x<1", NULL, NULL, NULL,
    "set tail loss prob for QDB calculation to <x>", 1 },
  { "--noqdb", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, "--beta", "do not use QDBs", 1 },
  { "--ilo", eslARG_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
    "set min parsetree score for accepted samples", 1 },
  { "--ihi", eslARG_REAL, NULL, NULL, NULL, NULL, NULL, NULL,
    "set max parsetree score for accepted samples", 1 },
  { "--imix", eslARG_REAL, NULL, NULL, NULL, NULL, NULL, "--exp",
    "set target avg parsetree score via null mixing", 1 },
  { "--isubtr", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, "--ilo,--ihi",
    "IS mode: sub-parsetree splice into null flanks", 1 },
  { "--imu", eslARG_REAL, NULL, NULL, NULL, NULL, "--isubtr", NULL,
    "target local Inside score (mu) for sub-parsetree IS", 1 },
  { "--imutol", eslARG_REAL, "2.0", NULL, "x>0", NULL, "--isubtr", NULL,
    "tolerance (bits) for sub-parsetree score match to --imu", 1 },
  { "--iflank",   eslARG_NONE, FALSE, NULL, NULL, NULL, "--isubtr", NULL,
    "IS mode: scan flanks only, add qc_sc directly as CM-region hit", 1 },
  { "--iflank-W", eslARG_INT,  "0",   NULL, "n>=0", NULL, "--iflank", NULL,
    "IS flank scan: cap window size W at <n> nt (0: use model W)", 1 },
  { "--ipaint",   eslARG_NONE, FALSE, NULL, NULL, NULL, "--isubtr", "--iflank",
    "painting mode: pack sub-scanned seqs into ~10Kb chunks", 1 },
  { "--ipaint-L", eslARG_INT,  "10000", NULL, "n>0", NULL, "--ipaint", NULL,
    "target mega-sequence length per chunk for --ipaint", 1 },
  { "--ipaint-lo", eslARG_REAL, NULL, NULL, NULL, NULL, "--ipaint", NULL,
    "only keep seqs with qc_sc >= <x> (Inside score filter)", 1 },
  { "--ipaint-hi", eslARG_REAL, NULL, NULL, NULL, NULL, "--ipaint", NULL,
    "only keep seqs with qc_sc < <x> (Inside score filter)", 1 },
  { "--ipaint-noqcsc", eslARG_NONE, FALSE, NULL, NULL, NULL, "--ipaint", NULL,
    "discard the v*-rooted [il..ir] hit; keep other overlapping hits", 1 },
  { "--ipaint-flankonly", eslARG_NONE, FALSE, NULL, NULL, NULL, "--ipaint", NULL,
    "discard hits overlapping CM region [il..ir]; keep flank hits only", 1 },
  { "--ipaint-allrand", eslARG_NONE, FALSE, NULL, NULL, NULL, "--ipaint", NULL,
    "also randomize CM region [il..ir] (pure random mega-seq)", 1 },
  { "--iinside-wt", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, "--isubtr",
    "use best Inside hit score for IS weight (not parsetree Viterbi)", 1 },
  { "--ipaint-qcsconly", eslARG_NONE, FALSE, NULL, NULL, NULL, "--ipaint", NULL,
    "collect only v*-rooted qc_sc hits (skip mega-seq scan)", 1 },
  { "--ipaint-sumv", eslARG_NONE, FALSE, NULL, NULL, NULL, "--ipaint", NULL,
    "use sum-over-all-v Inside score for IS weight (not just v*)", 1 },
  { "--ibest1", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "keep only the best (highest-scoring) hit per sequence", 1 },
  { "--ihbanded", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "use HMM-banded Inside scan instead of non-banded", 1 },
  { "--ifloat", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, "--ihbanded",
    "use float (not integer) unbanded Inside scan", 1 },
  { "--emit-cmfile", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
    "save emit_cm (alpha-mixed CM) to file <f>", 1 },
  { "--glocal", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, "--ilocal",
    "use glocal Inside mode (not local)", 1 },
  { "--ewt-lo", eslARG_REAL, NULL, NULL, NULL, NULL, "--iinside-wt", NULL,
    "reject if emit_cm Inside score < <x> (IS weight floor)", 1 },
  { "--ewt-hi", eslARG_REAL, NULL, NULL, NULL, NULL, "--iinside-wt", NULL,
    "reject if emit_cm Inside score > <x> (IS weight ceiling)", 1 },
  { "--tau", eslARG_REAL, "5e-6", NULL, "0<x<0.5", NULL, "--ihbanded", NULL,
    "set HMM band tail loss probability to <x>", 1 },
  { "--no-weight", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "diagnostic: use weight=1 for all IS seqs", 1 },
  { "--exp", eslARG_REAL, NULL, NULL, "x>0", NULL, NULL, "--imix",
    "exponentiate CM probabilities by <x> before sampling", 1 },
  { "--seed", eslARG_INT, "181", NULL, "n>=0", NULL, NULL, NULL,
    "set RNG seed to <n> (if 0: one-time arbitrary seed)", 1 },
  { "--mxsize", eslARG_REAL, "2048.0", NULL, "x>0.", NULL, NULL, NULL,
    "set max HMM banded DP mx size to <x> Mb", 1 },
  { "--ifile", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
    "save impt sample exp tail fits to <f>", 2 },
  { "--infit", eslARG_INT, "100", NULL, NULL, NULL, "--ifile", NULL,
    "with --ifile, do tail fits to <n> equally spaced tail probs", 2 },
  { "--imax", eslARG_REAL, "1.00", NULL, NULL, NULL, "--rfile", NULL,
    "with --ifile, max tail prob to fit is <x>", 2 },
  { "--imin", eslARG_REAL, "0.01", NULL, NULL, NULL, "--rfile", NULL,
    "with --ifile, min tail prob to fit is <x>", 2 },
  { "--no-rand", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "skip the expensive random sequence simulation", 1 },
  { "--isscfile", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
    "save IS scores/weights to <f> (score weight per line)", 1 },
  { "--refN", eslARG_INT, NULL, NULL, "n>0", NULL, NULL, NULL,
    "reference: search <n> random seqs of length clen", 1 },
  { "--refscfile", eslARG_OUTFILE, NULL, NULL, NULL, NULL, "--refN", NULL,
    "with --refN, save all ref scores to <f>", 1 },
  { "--reftailp", eslARG_REAL, "0.02", NULL, "0.0<x<0.6", NULL, "--refN", NULL,
    "with --refN, tail fraction to fit to exp", 1 },
  { "--rfile", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,
    "save random sample exp tail fits to <f>", 2 },
  { "--rnfit", eslARG_INT, "100", NULL, NULL, NULL, "--rfile", NULL,
    "with --rfile, do tail fits to <n> equally spaced tail probs", 2 },
  { "--rmax", eslARG_REAL, "0.10", NULL, NULL, NULL, "--rfile", NULL,
    "with --rfile, max tail prob to fit is <x>", 2 },
  { "--rmin", eslARG_REAL, "0.002", NULL, NULL, NULL, "--rfile", NULL,
    "with --rfile, min tail prob to fit is <x>", 2 },
  { "--idesign", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "design sequences with target Inside scores using Outside-guided mutations", 1 },
  { "--itarget", eslARG_REAL, "-5.0", NULL, NULL, NULL, "--idesign", NULL,
    "target Inside score for --idesign", 1 },
  { "--itarget-tol", eslARG_REAL, "0.5", NULL, "x>0", NULL, "--idesign", NULL,
    "tolerance in bits for --idesign convergence", 1 },
  { "--idesign-n", eslARG_INT, "10", NULL, "n>0", NULL, "--idesign", NULL,
    "number of sequences to design with --idesign", 1 },
  { "--istep", eslARG_REAL, "99.0", NULL, "x>0", NULL, "--idesign", NULL,
    "max score change per mutation step in bits (smaller = more intermediates)", 1 },
  { "--iwt-cap", eslARG_REAL, "10.0", NULL, "x>0", NULL, "--idesign", NULL,
    "stop collecting from a design run when weight exceeds this", 1 },
  { "--imcmc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
    "MCMC sampling from the tail of the null distribution", 1 },
  { "--imcmc-mu", eslARG_REAL, "-10.0", NULL, NULL, NULL, "--imcmc", NULL,
    "mu threshold: only sample sequences with Inside >= mu", 1 },
  { "--imcmc-chains", eslARG_INT, "10", NULL, "n>0", NULL, "--imcmc", NULL,
    "number of independent MCMC chains", 1 },
  { "--imcmc-steps", eslARG_INT, "100", NULL, "n>0", NULL, "--imcmc", NULL,
    "number of accepted steps per chain (after burn-in)", 1 },
  { "--imcmc-burnin", eslARG_INT, "20", NULL, "n>=0", NULL, "--imcmc", NULL,
    "number of burn-in steps to discard per chain", 1 },
  { "--imcmc-maxv", eslARG_NONE, FALSE, NULL, NULL, NULL, "--imcmc", NULL,
    "use max_v [beginsc[v]+alpha[v][L][L]] instead of alpha[0][L][L]", 1 },

  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 *
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 * This strategy is used despite the fact that a MPI version of cmemit does not
 * yet exist!
 */
struct cfg_s {
  char *cmfile;      /* name of input CM file  */
  CM_FILE *cmfp;     /* open input CM file stream       */
  ESL_ALPHABET *abc; /* digital alphabet for CM */
  ESL_RANDOMNESS *r; /* source of randomness */
  int ncm;           /* number CM we're at in file */
  int rN;            /* number of random sequences to search */
  int rL;            /* length of random sequences to search */
  int sN;            /* number of sampled sequences to search */
  float sT;          /* bit sc threshold for collecting hits from sampled seqs */
  int my_rank;

  /* optional output files */
  FILE *ifp;     /* output file for impt sample fits */
  FILE *rfp;     /* output file for random sample fits */
  FILE *refscfp; /* output file for reference scores (--refscfile) */
  FILE *isscfp;  /* output file for IS score+weight pairs (--isscfile) */
};

static char usage[] = "[-options] <cmfile>";
static char banner[] = "score random sequences with a covariance model";

static int init_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void master (const ESL_GETOPTS *go, struct cfg_s *cfg);

static int initialize_cm (const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, int do_local,
                          char *errbuf);
static int print_run_info (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int get_command (const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int collect_scores (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm,
                           CM_t *emit_cm, int N, int L, int *ret_scN, float **ret_scA,
                           float **ret_wtA, double *ret_dbsize);
static int fit_histogram (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float tailp,
                          double dbsize_nt, int do_impt, float *scores, float *weights, int nscores,
                          int exp_mode, double *ret_mu, double *ret_lambda, double *ret_nrandhits,
                          float *ret_tailp);
static int sample_sequence_from_cm (struct cfg_s *cfg, char *errbuf, CM_t *emit_cm, CM_t *score_cm,
                                    int *ret_L, ESL_DSQ **ret_dsq, Parsetree_t **ret_tr,
                                    float *ret_parsetree_sc);
static int impt_exp_FitComplete (double *x, double *w, int n, double *ret_mu, double *ret_lambda,
                                 double *ret_scaled_nhits);
static double cm_ExpectedParsetreeScore (CM_t *cm, double alpha);
static void cm_MixWithNull (CM_t *cm, double alpha);
static int cp9_EnforceQDBBands (CM_t *cm, CP9Bands_t *cp9b, CM_SCAN_MX *smx, int qdbidx, int L, char *errbuf);
static float cm_InsideScoreAfterMutation (CM_t *cm, CM_MX *ins_mx, CM_MX *out_mx, ESL_DSQ *dsq, int L, int pos, int new_res);
static int cm_DesignSequence (CM_t *cm, CM_t *emit_cm, struct cfg_s *cfg, const ESL_GETOPTS *go,
                              char *errbuf, float target_sc, float tol, int max_iter, int verbose,
                              ESL_DSQ **ret_dsq, int *ret_L, float *ret_sc, int *ret_niter,
                              double *ret_log_q);
static float cm_MaxLocalBeginScore (CM_t *cm, CM_MX *ins_mx, int L);
static float cm_BestLocalHitScore (CM_t *cm, CM_MX *ins_mx, int L);
static float cm_BestLocalHitScoreVJD (CM_t *cm, CM_MX *ins_mx, int L, int *ret_v, int *ret_j, int *ret_d);
static int cm_InsideAlign_partial (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L,
                                   CM_MX *mx, int p, float *ret_sc);
static int cm_MCMC_tail (CM_t *cm, CM_t *emit_cm, struct cfg_s *cfg, const ESL_GETOPTS *go,
                         char *errbuf, float mu, int n_chains, int n_steps, int n_burnin,
                         int use_maxv, int verbose, float **ret_scores, int *ret_N);

int
main (int argc, char **argv) {
  ESL_GETOPTS *go = NULL; /* command line processing                     */
  ESL_STOPWATCH *w = esl_stopwatch_Create ();
  if (w == NULL)
    cm_Fail ("Memory allocation error, stopwatch could not be created.");
  esl_stopwatch_Start (w);
  struct cfg_s cfg;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum ();
  FLogsumInit ();

  /***********************************************
   * Parse command line
   ***********************************************/

  /* Process command line options.
   */
  go = esl_getopts_Create (options);
  if (esl_opt_ProcessCmdline (go, argc, argv) != eslOK || esl_opt_VerifyConfig (go) != eslOK) {
    printf ("Failed to parse command line: %s\n", go->errbuf);
    esl_usage (stdout, argv[0], usage);
    printf ("\nTo see more help on available options, do %s -h\n\n", argv[0]);
    exit (1);
  }
  if (esl_opt_GetBoolean (go, "-h") == TRUE) {
    cm_banner (stdout, argv[0], banner);
    esl_usage (stdout, argv[0], usage);
    puts ("\nwhere general options are:");
    esl_opt_DisplayHelp (stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
    puts ("\nmiscellaneous output options are:");
    esl_opt_DisplayHelp (stdout, go, 2, 2, 80);
    exit (0);
  }
  if (esl_opt_ArgNumber (go) != 1) {
    puts ("Incorrect number of command line arguments.");
    esl_usage (stdout, argv[0], usage);
    puts ("\n  where basic options are:");
    esl_opt_DisplayHelp (stdout, go, 1, 2, 80);
    printf ("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
    exit (1);
  }
  /* Initialize what we can in the config structure (without knowing the alphabet yet).
   * We could assume RNA, but this HMMER3 based approach is more general.
   */
  cfg.cmfile = esl_opt_GetArg (go, 1);
  cfg.cmfp = NULL; /* opened in init_cfg() */
  cfg.abc = NULL;  /* created in init_cfg() */
  cfg.r = NULL;    /* created in init_cfg() */

  cfg.ifp = NULL;
  cfg.rfp = NULL;
  cfg.refscfp = NULL;
  cfg.isscfp = NULL;

  cm_banner (stdout, argv[0], banner);

  /* do work */
  master (go, &cfg);

  /* Clean up the cfg.
   */
  if (cfg.abc != NULL) {
    esl_alphabet_Destroy (cfg.abc);
    cfg.abc = NULL;
  }
  if (cfg.cmfp != NULL)
    cm_file_Close (cfg.cmfp);
  if (cfg.r != NULL)
    esl_randomness_Destroy (cfg.r);

  /* master specific cleaning */
  if (cfg.ifp != NULL) {
    fclose (cfg.ifp);
    printf ("# Important sampling fits to various tail masses saved to file %s.\n",
            esl_opt_GetString (go, "--ifile"));
  }
  if (cfg.rfp != NULL) {
    fclose (cfg.rfp);
    printf ("# Random sequence histogram fits to various tail masses saved to file %s.\n",
            esl_opt_GetString (go, "--rfile"));
  }

  esl_getopts_Destroy (go);
  esl_stopwatch_Stop (w);
  printf ("#\n");
  esl_stopwatch_Display (stdout, w, "# CPU time: ");
  esl_stopwatch_Destroy (w);
  return 0;
}

/* init_cfg()
 * Already set:
 *    cfg->cmfile  - command line arg 1
 * Sets:
 *    cfg->cmfp    - open CM file
 *    cfg->ifp     - optional output file (--ifile)
 *    cfg->rfp     - optional output file (--rfile)
 *    cfg->r       - source of randomness
 */
static int
init_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf) {
  int status;

  /* open CM file for reading */
  if ((status = cm_file_Open (cfg->cmfile, NULL, FALSE, &(cfg->cmfp), errbuf)) != eslOK)
    return status;

  /* open optional output files, if nec */
  if (esl_opt_GetString (go, "--ifile") != NULL) {
    if ((cfg->ifp = fopen (esl_opt_GetString (go, "--ifile"), "w")) == NULL)
      ESL_FAIL (eslFAIL, errbuf, "Failed to open important sampling fit save file %s for writing\n",
                esl_opt_GetString (go, "--ifile"));
  }
  if (esl_opt_GetString (go, "--rfile") != NULL) {
    if ((cfg->rfp = fopen (esl_opt_GetString (go, "--rfile"), "w")) == NULL)
      ESL_FAIL (eslFAIL, errbuf, "Failed to open important sampling fit save file %s for writing\n",
                esl_opt_GetString (go, "--rfile"));
  }
  if (esl_opt_GetString (go, "--refscfile") != NULL) {
    if ((cfg->refscfp = fopen (esl_opt_GetString (go, "--refscfile"), "w")) == NULL)
      ESL_FAIL (eslFAIL, errbuf, "Failed to open reference scores file %s for writing\n",
                esl_opt_GetString (go, "--refscfile"));
  }
  if (esl_opt_GetString (go, "--isscfile") != NULL) {
    if ((cfg->isscfp = fopen (esl_opt_GetString (go, "--isscfile"), "w")) == NULL)
      ESL_FAIL (eslFAIL, errbuf, "Failed to open IS scores file %s for writing\n",
                esl_opt_GetString (go, "--isscfile"));
  }

  /* create RNG */
  cfg->r = esl_randomness_CreateFast (esl_opt_GetInteger (go, "--seed"));
  if (cfg->r == NULL)
    ESL_FAIL (eslEINVAL, errbuf,
              "Failed to create random number generator: probably out of memory");

  cfg->rN = esl_opt_GetInteger (go, "--rN");
  cfg->rL = esl_opt_GetInteger (go, "--rL");
  cfg->sN = esl_opt_GetInteger (go, "--iN");

  cfg->my_rank = 0;

  return eslOK;
}

/* master()
 * The serial version of cmsim. (There is no parallel version yet).
 * For each CM, generate random sequences and search them.
 *
 * We only return if successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
master (const ESL_GETOPTS *go, struct cfg_s *cfg) {
  int status;
  char errbuf[eslERRBUFSIZE];
  CM_t *cm = NULL;
  CM_t *emit_cm = NULL; /* exponentiated CM for emission (NULL if no --exp) */
  int exp_mode;         /* exp tail mode */
  int rscN = 0;         /* number of hits in random seqs reported thus far, for all seqs */
  float *rscA = NULL;   /* [0..rscN-1] hit scores for all random seqs */
  float *rwtA = NULL; /* [0..rscN-1] importance weights for random seqs (not used, will be NULL) */
  int sscN = 0;       /* number of hits in CM-sampled seqs reported thus far, for all seqs */
  float *sscA = NULL; /* [0..sscN-1] hit scores for all CM-sampled seqs */
  float *swtA = NULL; /* [0..sscN-1] importance weights for CM-sampled seqs */
  int mscN = 0;       /* number of hits in matched-length random seqs */
  float *mscA = NULL; /* [0..mscN-1] hit scores for matched-length random seqs */
  float *mwtA = NULL; /* not used, will be NULL */
  float min_wt = 1.;
  float max_wt = 0.;
  float sum_wt = 0.;
  float sum_wt_sq = 0.;
  double ess; /* effective sample size */
  ExpInfo_t *impt_expinfo;
  ExpInfo_t *rand_expinfo;
  ExpInfo_t *match_expinfo; /* matched-length random control */

  double mu, lambda;   /* temporary mu and lambda used for setting exp tails */
  double nrandhits;    /* temporary number of rand hits found */
  double nsamphits;    /* temporary number of rand hits found */
  double nmatchhits;   /* number of hits in matched-length random seqs */
  float tailp;         /* temporary tail mass probability fit to an exponential */
  float sc_tailp;      /* scaled tailp */
  double rand_dbsize;  /* effective nt searched in random seq run */
  double is_dbsize;    /* IS-weighted effective nt searched in IS run */
  double match_dbsize; /* effective nt searched in matched-length random run */
  float avg_hitlen;
  float tailp_step; /* size of change in tailp parameter for --ifile, --rfile */
  int nfits;        /* number of fits for --ifile, --rfile */
  int i;            /* counter over fits */

  if ((status = init_cfg (go, cfg, errbuf)) != eslOK)
    cm_Fail (errbuf);
  if ((status = print_run_info (go, cfg, errbuf)) != eslOK)
    cm_Fail (errbuf);

  cfg->ncm = 0;

  while ((status = cm_file_Read (cfg->cmfp, TRUE, &(cfg->abc), &cm)) == eslOK) {
    if (cm == NULL)
      cm_Fail ("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
    cfg->ncm++;

    ESL_ALLOC (impt_expinfo, sizeof (ExpInfo_t *));
    impt_expinfo = CreateExpInfo ();

    ESL_ALLOC (rand_expinfo, sizeof (ExpInfo_t *));
    rand_expinfo = CreateExpInfo ();

    ESL_ALLOC (match_expinfo, sizeof (ExpInfo_t *));
    match_expinfo = CreateExpInfo ();

    /* Create proposal distribution for importance sampling.
     * --imix <target_sc>: mix CM emissions with null to achieve
     *   target average parsetree score (binary search for alpha).
     * --exp <x>: exponentiate CM probabilities by x (legacy method).
     * The original CM is always used for searching. */
    emit_cm = NULL;
    if (esl_opt_IsOn (go, "--imix")) {
      /* Need tsc/beginsc configured for cm_ExpectedParsetreeScore().
       * Configure a temporary clone to get transition scores, then use it
       * for the expected score calculation and free it afterward. */
      CM_t *tmp_cm = NULL;
      if ((status = cm_Clone (cm, errbuf, &tmp_cm)) != eslOK)
        cm_Fail (errbuf);
      if ((status = initialize_cm (go, cfg, tmp_cm, !esl_opt_GetBoolean(go, "--glocal"), errbuf)) != eslOK)
        cm_Fail (errbuf);

      double target_sc = esl_opt_GetReal (go, "--imix");
      double orig_sc = cm_ExpectedParsetreeScore (tmp_cm, 0.0);
      double lo = 0.0, hi = 1.0, mid, mid_sc;
      int iter;

      if (target_sc >= orig_sc) {
        printf ("Warning: --imix target %.1f >= expected score %.1f, using unmodified CM\n",
                target_sc, orig_sc);
      } else {
        if ((status = cm_Clone (cm, errbuf, &emit_cm)) != eslOK)
          cm_Fail (errbuf);

        /* Binary search for alpha that gives target expected score */
        for (iter = 0; iter < 100; iter++) {
          mid = (lo + hi) / 2.0;
          mid_sc = cm_ExpectedParsetreeScore (tmp_cm, mid);
          if (fabs (mid_sc - target_sc) < 0.01)
            break; /* close enough */
          if (mid_sc > target_sc)
            lo = mid; /* score too high, mix more null */
          else
            hi = mid; /* score too low, mix less null */
        }
        printf ("--imix: target_sc=%.1f orig_sc=%.1f alpha=%.6f achieved_sc=%.2f (%d "
                "iterations)\n",
                target_sc, orig_sc, mid, mid_sc, iter);
        cm_MixWithNull (emit_cm, mid);

        /* Save emit_cm BEFORE configure (cm_file_WriteASCII needs unconfigured CM) */
        if (esl_opt_IsOn (go, "--emit-cmfile")) {
          FILE *emit_fp = fopen (esl_opt_GetString (go, "--emit-cmfile"), "w");
          if (emit_fp == NULL)
            cm_Fail ("Failed to open emit_cm output file %s", esl_opt_GetString (go, "--emit-cmfile"));
          if ((status = cm_file_WriteASCII (emit_fp, -1, emit_cm)) != eslOK)
            cm_Fail ("cm_file_WriteASCII failed for emit_cm (status=%d)", status);
          fclose (emit_fp);
          printf ("Saved emit_cm (alpha=%.6f) to %s\n", mid, esl_opt_GetString (go, "--emit-cmfile"));
        }

        if ((status = initialize_cm (go, cfg, emit_cm, !esl_opt_GetBoolean(go, "--glocal"), errbuf)) != eslOK)
          cm_Fail (errbuf);
      }
      FreeCM (tmp_cm);
    } else if (esl_opt_IsOn (go, "--exp")) {
      if ((status = cm_Clone (cm, errbuf, &emit_cm)) != eslOK)
        cm_Fail (errbuf);
      cm_Exponentiate (emit_cm, esl_opt_GetReal (go, "--exp"));
      if ((status = initialize_cm (go, cfg, emit_cm, !esl_opt_GetBoolean(go, "--glocal"), errbuf)) != eslOK)
        cm_Fail (errbuf);
    }

    /* (emit_cm saved above, before initialize_cm, if --emit-cmfile was set) */

    { int dbg_k; long dbg_sum = 0;
      if (cm->cp9 != NULL && (cm->cp9->flags & CPLAN9_HASBITS)) {
        for (dbg_k = 0; dbg_k <= cm->cp9->M; dbg_k++) {
          dbg_sum += cm->cp9->msc[0][dbg_k] + cm->cp9->isc[0][dbg_k];
          dbg_sum += cm->cp9->tsc[CTMM][dbg_k] + cm->cp9->tsc[CTMI][dbg_k];
          dbg_sum += cm->cp9->bsc[dbg_k] + cm->cp9->esc[dbg_k];
        }
      }
      { double t_sum = 0., e_sum = 0.; int dbg_v2;
        for (dbg_v2 = 0; dbg_v2 < cm->M; dbg_v2++) {
          t_sum += cm->t[dbg_v2][0];
          e_sum += cm->e[dbg_v2][0];
        }
        /* Full CM dump: checksum all float arrays */
        { double all_t = 0., all_e = 0., all_tsc = 0., all_esc = 0.;
          for (dbg_v2 = 0; dbg_v2 < cm->M; dbg_v2++) {
            int dbg_j;
            for (dbg_j = 0; dbg_j < MAXCONNECT; dbg_j++) all_t += cm->t[dbg_v2][dbg_j];
            for (dbg_j = 0; dbg_j < cm->abc->K * cm->abc->K; dbg_j++) all_e += cm->e[dbg_v2][dbg_j];
            for (dbg_j = 0; dbg_j < MAXCONNECT; dbg_j++) all_tsc += cm->tsc[dbg_v2][dbg_j];
            for (dbg_j = 0; dbg_j < cm->abc->K * cm->abc->K; dbg_j++) all_esc += cm->esc[dbg_v2][dbg_j];
          }
          printf ("BEFORE init_cm: t=%.10f e=%.10f tsc=%.10f esc=%.10f config=0x%x search=0x%x flags=0x%x\n",
                  all_t, all_e, all_tsc, all_esc, cm->config_opts, cm->search_opts, cm->flags);
        }
      }
    }

    if ((status = initialize_cm (go, cfg, cm, !esl_opt_GetBoolean(go, "--glocal"), errbuf)) != eslOK)
      cm_Fail (errbuf);

    { int dbg_k; long dbg_sum = 0;
      for (dbg_k = 0; dbg_k <= cm->cp9->M; dbg_k++) {
        dbg_sum += cm->cp9->msc[0][dbg_k] + cm->cp9->isc[0][dbg_k];
        dbg_sum += cm->cp9->tsc[CTMM][dbg_k] + cm->cp9->tsc[CTMI][dbg_k];
        dbg_sum += cm->cp9->bsc[dbg_k] + cm->cp9->esc[dbg_k];
      }
      printf ("AFTER init cm: cp9_chksum=%ld flags=0x%x\n", dbg_sum, cm->cp9->flags);
    }

    printf ("CM %d: %s\n", cfg->ncm, cm->name);

    /* --refN mode: search N random sequences of length clen, output all scores,
     * fit exp tail, then skip IS and long random seq work for this CM. */
    if (esl_opt_IsOn (go, "--refN")) {
      int refN = esl_opt_GetInteger (go, "--refN");
      float reftailp = esl_opt_GetReal (go, "--reftailp");
      int refscN = 0;
      float *refscA = NULL;
      float *refwtA = NULL;
      double ref_mu, ref_lambda, ref_nhits;
      int j;
      cm->search_opts |= CM_SEARCH_INSIDE; /* local Inside (do_local=TRUE above) */
      if ((status = collect_scores (go, cfg, errbuf, cm, NULL, refN, cm->clen, &refscN, &refscA,
                                    &refwtA, NULL))
          != eslOK)
        cm_Fail (errbuf);
      if ((status = fit_histogram (go, cfg, errbuf, reftailp, 0., FALSE, refscA, NULL, refscN,
                                   EXP_CM_LI, &ref_mu, &ref_lambda, &ref_nhits, &reftailp))
          != eslOK)
        cm_Fail (errbuf);
      printf ("Reference (L=%d, N=%d) fit:\n\t%12s: %9.5f\n\t%12s: %9.5f\n\t%12s: "
              "%9.0f\n\t%12s: %9.5f\n\n",
              cm->clen, refN, "mu", ref_mu, "lambda", ref_lambda, "nhits", ref_nhits, "tailp",
              reftailp);
      if (cfg->refscfp != NULL) {
        for (j = 0; j < refscN; j++)
          fprintf (cfg->refscfp, "%.4f\n", refscA[j]);
      }
      free (refscA);
      FreeCM (cm);
      continue; /* skip IS and long random seq for this CM */
    }

    /* Set exp_mode: local or glocal Inside */
    exp_mode = esl_opt_GetBoolean (go, "--glocal") ? EXP_CM_GI : EXP_CM_LI;
    /* set CM_SEARCH_INSIDE flag for Inside mode (CYK is the default) */
    if (ExpModeIsInside (exp_mode)) {
      cm->search_opts |= CM_SEARCH_INSIDE;
      if (emit_cm != NULL)
        emit_cm->search_opts |= CM_SEARCH_INSIDE;
    } else {
      cm->search_opts &= ~CM_SEARCH_INSIDE;
      if (emit_cm != NULL)
        emit_cm->search_opts &= ~CM_SEARCH_INSIDE;
    }

    /* Search random sequences and collect score histograms.
     * Diagnostic: if --no-weight, use emit_cm (the modified proposal) for
     * random sequence searching so random and IS use the same model.
     * Skip this expensive block if --no-rand is set. */
    nrandhits = 0.;
    if (!esl_opt_GetBoolean (go, "--no-rand")) {
      {
        CM_t *rand_cm = (esl_opt_GetBoolean (go, "--no-weight") && emit_cm != NULL) ? emit_cm : cm;
        if ((status = collect_scores (go, cfg, errbuf, rand_cm, NULL, cfg->rN, cfg->rL, &rscN,
                                      &rscA, &rwtA, &rand_dbsize))
            != eslOK)
          cm_Fail (errbuf);
      }
      tailp = esl_opt_IsOn (go, "--rtailp") ? esl_opt_GetReal (go, "--rtailp") : 0.;
      if ((status = fit_histogram (go, cfg, errbuf, tailp, rand_dbsize, FALSE, rscA, NULL, rscN,
                                   exp_mode, &mu, &lambda, &nrandhits, &tailp))
          != eslOK)
        cm_Fail (errbuf);
      avg_hitlen = rand_dbsize / (double)nrandhits;
      printf ("Random  seq fit histogram:\n\t%12s: %9.5f\n\t%12s: %9.5f\n\t%12s: "
              "%9.5f\n\t%12s: %9.5f\n\t%12s: %9.5f\n\n",
              "mu", mu, "lambda", lambda, "nrandhits", nrandhits, "tailp", tailp, "avg_len",
              avg_hitlen);
      SetExpInfo (rand_expinfo, lambda, mu, (double)(cfg->rL * cfg->rN), (int)nrandhits, tailp);
      debug_print_expinfo (rand_expinfo);

      /* output to --rfile, if nec */
      if (cfg->rfp != NULL) {
        tailp = esl_opt_GetReal (go, "--rmax");
        nfits = esl_opt_GetInteger (go, "--rnfit");
        tailp_step = (tailp - esl_opt_GetReal (go, "--rmin")) / (float)nfits;
        for (i = 0; i < nfits; i++) {
          if ((status = fit_histogram (go, cfg, errbuf, tailp, 0., FALSE, rscA, NULL, rscN,
                                       exp_mode, &mu, &lambda, &nrandhits, NULL))
              != eslOK)
            cm_Fail (errbuf);
          fprintf (cfg->rfp, "%g  %g  %g  %g  %g\n", tailp, lambda,
                   (mu - log (1. / tailp) / lambda), mu, nrandhits);
          tailp -= tailp_step;
        }
      }
    }

    /* --idesign mode: Outside-guided stochastic mutations to collect IS
     * samples in the tail. Runs N trajectories; at each iteration of each
     * trajectory, collects (score, weight) if score >= itarget AND weight
     * <= iwt_cap. Stops a trajectory when weight exceeds the cap.
     * Tests three mutation strategies to find best lambda match. */
    if (esl_opt_GetBoolean (go, "--idesign")) {
      float itarget   = (float) esl_opt_GetReal (go, "--itarget");
      float wt_cap    = (float) esl_opt_GetReal (go, "--iwt-cap");
      int   idesign_n = esl_opt_GetInteger (go, "--idesign-n");
      int   verbose   = esl_opt_GetBoolean (go, "-v");
      float max_step  = (float) esl_opt_GetReal (go, "--istep");
      CM_t *ecm       = (emit_cm != NULL) ? emit_cm : cm;
      int   K         = cm->abc->K;
      float beta_softmax = 2.0;

      printf ("Design-collect mode: mu=%.3f wt_cap=%.1f N=%d max_step=%.1f\n",
              itarget, wt_cap, idesign_n, max_step);

      int total_collected = 0;
      int total_trajectories = 0;

      for (i = 0; i < idesign_n; i++) {
        ESL_DSQ *dsq = NULL;
        Parsetree_t *tr = NULL;
        ESL_SQ *sq = NULL;
        int L;
        float sc_ptree, sc_inside, emit_cm_inside;
        double log_q;
        CM_MX *ins_mx = NULL, *out_mx = NULL;

        /* Generate starting sequence: from emit_cm if available, else random.
         * With emit_cm: log_q = emit_cm_inside (corrects for CM bias).
         * Without (random): log_q = 0 (weight starts at 1). */
        if (ecm != cm) {
          /* Emit from emit_cm (alpha-mixed proposal) */
          if ((status = EmitParsetree (ecm, errbuf, cfg->r, "design", TRUE, &tr, &sq, &L)) != eslOK)
            cm_Fail (errbuf);
          while (L == 0) {
            esl_sq_Destroy (sq);
            if ((status = EmitParsetree (ecm, errbuf, cfg->r, "design", TRUE, &tr, &sq, &L)) != eslOK)
              cm_Fail (errbuf);
          }
          ESL_ALLOC (dsq, sizeof (ESL_DSQ) * (sq->n + 2));
          memcpy (dsq, sq->dsq, sizeof (ESL_DSQ) * (sq->n + 2));
          esl_sq_Destroy (sq);
          if (tr != NULL) { FreeParsetree (tr); tr = NULL; }

          /* Compute initial log_q from emit_cm Inside */
          CM_MX *emit_mx = cm_mx_Create (ecm->M);
          if ((status = cm_InsideAlign (ecm, errbuf, dsq, L, 512.0, emit_mx, &emit_cm_inside)) != eslOK)
            cm_Fail (errbuf);
          cm_mx_Destroy (emit_mx);
          log_q = (double) emit_cm_inside;
        } else {
          /* No emit_cm: generate random i.i.d. null sequence */
          L = cm->clen;
          ESL_ALLOC (dsq, sizeof (ESL_DSQ) * (L + 2));
          if ((status = esl_rsq_xfIID (cfg->r, cm->null, cm->abc->K, L, dsq)) != eslOK)
            cm_Fail ("ERROR generating random sequence");
          log_q = 0.;
        }
        total_trajectories++;

        ins_mx = cm_mx_Create (cm->M);
        out_mx = cm_mx_Create (cm->M);

        int max_cand = 3 * (L + 1);
        float *cand_sc;
        int   *cand_pos, *cand_res;
        ESL_ALLOC (cand_sc,  sizeof (float) * max_cand);
        ESL_ALLOC (cand_pos, sizeof (int)   * max_cand);
        ESL_ALLOC (cand_res, sizeof (int)   * max_cand);

        int iter;
        for (iter = 0; iter < 200; iter++) {
          /* Run Inside */
          if ((status = cm_InsideAlign (cm, errbuf, dsq, L, 512.0, ins_mx, &sc_inside)) != eslOK)
            cm_Fail (errbuf);

          double weight = pow (2.0, -log_q);

          /* Collect this intermediate if in tail AND weight is acceptable */
          if (sc_inside >= itarget && weight <= wt_cap) {
            if (cfg->isscfp != NULL)
              fprintf (cfg->isscfp, "%.4f\t%.6g\n", sc_inside, weight);
            total_collected++;
            if (verbose)
              printf ("  COLLECT traj=%d iter=%d: Inside=%.3f weight=%.4g log_q=%.3f\n",
                      i, iter, sc_inside, weight, log_q);
          }

          /* Stop if weight already exceeds cap */
          if (weight > wt_cap) {
            if (verbose)
              printf ("  STOP traj=%d iter=%d: weight=%.4g > cap=%.1f\n",
                      i, iter, weight, wt_cap);
            break;
          }

          /* Run Outside */
          if ((status = cm_OutsideAlign (cm, errbuf, dsq, L, 512.0, FALSE, out_mx, ins_mx, NULL)) != eslOK)
            cm_Fail (errbuf);

          /* Enumerate candidate mutations */
          int n_cand = 0;
          float cur_dist = fabs (sc_inside - itarget);
          for (int p = 1; p <= L; p++) {
            for (int r = 0; r < K; r++) {
              if (r == dsq[p]) continue;
              float new_sc = cm_InsideScoreAfterMutation (cm, ins_mx, out_mx, dsq, L, p, r);
              float new_dist = fabs (new_sc - itarget);
              float step_size = fabs (new_sc - sc_inside);
              if (new_dist < cur_dist && step_size <= max_step) {
                cand_sc[n_cand]  = new_sc;
                cand_pos[n_cand] = p;
                cand_res[n_cand] = r;
                n_cand++;
              }
            }
          }
          /* Relax step limit if no candidates found */
          if (n_cand == 0) {
            for (int p = 1; p <= L; p++) {
              for (int r = 0; r < K; r++) {
                if (r == dsq[p]) continue;
                float new_sc = cm_InsideScoreAfterMutation (cm, ins_mx, out_mx, dsq, L, p, r);
                float new_dist = fabs (new_sc - itarget);
                if (new_dist < cur_dist) {
                  cand_sc[n_cand]  = new_sc;
                  cand_pos[n_cand] = p;
                  cand_res[n_cand] = r;
                  n_cand++;
                }
              }
            }
          }
          if (n_cand == 0) break;

          /* Softmax selection */
          double *prob;
          ESL_ALLOC (prob, sizeof (double) * n_cand);
          double max_logp = -eslINFINITY;
          for (int c = 0; c < n_cand; c++) {
            prob[c] = -beta_softmax * fabs (cand_sc[c] - itarget);
            if (prob[c] > max_logp) max_logp = prob[c];
          }
          double Z = 0.;
          for (int c = 0; c < n_cand; c++) {
            prob[c] = exp (prob[c] - max_logp);
            Z += prob[c];
          }
          for (int c = 0; c < n_cand; c++) prob[c] /= Z;

          /* Sample */
          double u = esl_random (cfg->r);
          double cum = 0.;
          int chosen = n_cand - 1;
          for (int c = 0; c < n_cand; c++) {
            cum += prob[c];
            if (u <= cum) { chosen = c; break; }
          }

          log_q += log2 (prob[chosen]);
          dsq[cand_pos[chosen]] = cand_res[chosen];
          free (prob);
        } /* end iterations */

        free (cand_sc); free (cand_pos); free (cand_res);
        cm_mx_Destroy (ins_mx);
        cm_mx_Destroy (out_mx);
        free (dsq);
      } /* end trajectories */

      printf ("Collected %d data points from %d trajectories\n\n",
              total_collected, total_trajectories);
    }

    /* Diagnostic: test cm_InsideAlign_partial for correctness.
     * Procedure:
     *  1. Generate random sequence x_A
     *  2. Compute alpha_A using full cm_InsideAlign
     *  3. Mutate one residue to get x_B
     *  4. Compute alpha_B_full using full cm_InsideAlign on x_B
     *  5. Compute alpha_B_partial using partial DP starting from alpha_A
     *  6. Compare alpha_B_full vs alpha_B_partial — should be identical
     */
    if (esl_opt_GetBoolean (go, "--imcmc")) {
      printf ("DIAG: testing cm_InsideAlign_partial correctness on 20 sequences\n");
      int test_L = cm->clen;
      CM_MX *test_mx_a = cm_mx_Create (cm->M);
      CM_MX *test_mx_b = cm_mx_Create (cm->M);
      int n_match = 0, n_close = 0, n_diff = 0;
      float max_diff = 0.;

      for (int ti = 0; ti < 20; ti++) {
        ESL_DSQ *test_dsq;
        ESL_ALLOC (test_dsq, sizeof (ESL_DSQ) * (test_L + 2));
        esl_rsq_xfIID (cfg->r, cm->null, cm->abc->K, test_L, test_dsq);

        /* alpha_A: full Inside on original sequence */
        float sc_a;
        cm_InsideAlign (cm, errbuf, test_dsq, test_L, 512.0, test_mx_a, &sc_a);

        /* Pick a random mutation */
        int mut_p = 1 + esl_rnd_Roll (cfg->r, test_L);
        int mut_r = esl_rnd_Roll (cfg->r, cm->abc->K - 1);
        if (mut_r >= test_dsq[mut_p]) mut_r++;
        int old_res = test_dsq[mut_p];
        test_dsq[mut_p] = mut_r;

        /* alpha_B_full: full Inside on mutated sequence */
        float sc_b_full;
        cm_InsideAlign (cm, errbuf, test_dsq, test_L, 512.0, test_mx_b, &sc_b_full);

        /* alpha_B_partial: partial Inside starting from alpha_A
         * (test_mx_a was filled for the original sequence; we mutated dsq) */
        float sc_b_partial;
        cm_InsideAlign_partial (cm, errbuf, test_dsq, test_L, test_mx_a, mut_p, &sc_b_partial);

        float diff = fabs (sc_b_full - sc_b_partial);

        /* Also check the full matrix - cm_BestLocalHitScore should match */
        float best_full = cm_BestLocalHitScore (cm, test_mx_b, test_L);
        float best_partial = cm_BestLocalHitScore (cm, test_mx_a, test_L);
        float best_diff = fabs (best_full - best_partial);

        /* Worst-case cell-by-cell diff over the matrix */
        float worst_cell_diff = 0.;
        for (int vv = 0; vv < cm->M; vv++) {
          for (int jj = 0; jj <= test_L; jj++) {
            for (int dd = 0; dd <= jj; dd++) {
              if (NOT_IMPOSSIBLE(test_mx_a->dp[vv][jj][dd]) &&
                  NOT_IMPOSSIBLE(test_mx_b->dp[vv][jj][dd])) {
                float cd = fabs(test_mx_a->dp[vv][jj][dd] - test_mx_b->dp[vv][jj][dd]);
                if (cd > worst_cell_diff) worst_cell_diff = cd;
              }
            }
          }
        }

        if (diff < 1e-4 && best_diff < 1e-4 && worst_cell_diff < 1e-3) n_match++;
        else if (diff < 0.01) n_close++;
        else n_diff++;
        if (diff > max_diff) max_diff = diff;

        printf ("  test %2d: pos=%d %c->%c  root_full=%.4f root_part=%.4f  best_full=%.4f best_part=%.4f  worst_cell=%.4f\n",
                ti, mut_p, "ACGU"[old_res], "ACGU"[mut_r],
                sc_b_full, sc_b_partial, best_full, best_partial, worst_cell_diff);

        free (test_dsq);
      }
      printf ("  Summary: %d exact match, %d close (<0.01), %d different (max diff=%.4f)\n\n",
              n_match, n_close, n_diff, max_diff);

      /* ITERATIVE TEST: do several mutations in sequence and check that the
       * partial DP matches full DP after each one */
      printf ("DIAG: testing iterative partial DP (10 sequential mutations)\n");
      ESL_DSQ *iter_dsq;
      ESL_ALLOC (iter_dsq, sizeof (ESL_DSQ) * (test_L + 2));
      esl_rsq_xfIID (cfg->r, cm->null, cm->abc->K, test_L, iter_dsq);

      /* Initial: full Inside on starting sequence into test_mx_a */
      float iter_sc;
      cm_InsideAlign (cm, errbuf, iter_dsq, test_L, 512.0, test_mx_a, &iter_sc);
      printf ("  iter 0: full Inside sc = %.4f\n", iter_sc);

      for (int iter = 1; iter <= 10; iter++) {
        /* Random mutation */
        int mp = 1 + esl_rnd_Roll (cfg->r, test_L);
        int mr = esl_rnd_Roll (cfg->r, cm->abc->K - 1);
        if (mr >= iter_dsq[mp]) mr++;
        int old_r = iter_dsq[mp];
        iter_dsq[mp] = mr;

        /* Apply partial DP to test_mx_a (the running matrix) */
        float partial_sc;
        cm_InsideAlign_partial (cm, errbuf, iter_dsq, test_L, test_mx_a, mp, &partial_sc);

        /* Compute full DP into test_mx_b for ground truth */
        float full_sc;
        cm_InsideAlign (cm, errbuf, iter_dsq, test_L, 512.0, test_mx_b, &full_sc);

        /* Compare scores */
        float root_diff = fabs (full_sc - partial_sc);

        /* Compare full matrix */
        float worst_cell = 0.;
        int n_bad = 0;
        for (int vv = 0; vv < cm->M; vv++) {
          for (int jj = 0; jj <= test_L; jj++) {
            for (int dd = 0; dd <= jj; dd++) {
              if (NOT_IMPOSSIBLE(test_mx_a->dp[vv][jj][dd]) &&
                  NOT_IMPOSSIBLE(test_mx_b->dp[vv][jj][dd])) {
                float cd = fabs(test_mx_a->dp[vv][jj][dd] - test_mx_b->dp[vv][jj][dd]);
                if (cd > worst_cell) worst_cell = cd;
                if (cd > 0.001) n_bad++;
              }
            }
          }
        }

        float best_full = cm_BestLocalHitScore (cm, test_mx_b, test_L);
        float best_partial = cm_BestLocalHitScore (cm, test_mx_a, test_L);

        printf ("  iter %2d: pos=%d %c->%c  full=%.4f partial=%.4f diff=%.4f  best_full=%.4f best_partial=%.4f  worst_cell=%.4f n_bad=%d\n",
                iter, mp, "ACGU"[old_r], "ACGU"[mr],
                full_sc, partial_sc, root_diff,
                best_full, best_partial, worst_cell, n_bad);
      }
      free (iter_dsq);

      cm_mx_Destroy (test_mx_a);
      cm_mx_Destroy (test_mx_b);
    }

    /* Diagnostic: compare 4 scoring methods on the same random sequences:
     *  1. align_dL    = cm_InsideAlign alpha[0][L][L]            (sum at d=L)
     *  2. align_maxLL = max_v [beginsc[v] + alpha[v][L][L]]      (max at d=L)
     *  3. align_max   = max over all (v,j,d) [bsc + alpha]       (max over all)
     *  4. scan_best   = best scan hit from FastIInsideScan
     */
    if (0 && esl_opt_GetBoolean (go, "--imcmc")) {
      int diag_n = 1000;
      printf ("DIAG: 4 scoring methods on %d random seqs (L=%d):\n", diag_n, cm->clen);
      int diag_L = cm->clen;
      CM_MX *diag_mx = cm_mx_Create (cm->M);

      /* Open output file for all 4 columns + best (v,j,d) location */
      FILE *diag_fp = fopen ("/tmp/diag_scores.txt", "w");
      if (diag_fp == NULL) cm_Fail ("could not open /tmp/diag_scores.txt");
      fprintf (diag_fp, "# align_dL align_maxLL align_max scan_best best_v best_j best_d\n");

      for (int di = 0; di < diag_n; di++) {
        ESL_DSQ *diag_dsq;
        ESL_ALLOC (diag_dsq, sizeof (ESL_DSQ) * (diag_L + 2));
        esl_rsq_xfIID (cfg->r, cm->null, cm->abc->K, diag_L, diag_dsq);

        /* Method 1: cm_InsideAlign alpha[0][L][L] */
        float align_dL;
        cm_InsideAlign (cm, errbuf, diag_dsq, diag_L, 512.0, diag_mx, &align_dL);

        /* Method 2: max_v [beginsc[v] + alpha[v][L][L]] */
        float align_maxLL = cm_MaxLocalBeginScore (cm, diag_mx, diag_L);

        /* Method 3: max over all (v,j,d) [beginsc[v] + alpha[v][j][d]]
         * Also record which (v,j,d) is the best */
        float align_max = IMPOSSIBLE;
        int best_v = -1, best_j = -1, best_d = -1;
        for (int v = 0; v < cm->M; v++) {
          float bsc = (v == 0) ? 0.0f : (NOT_IMPOSSIBLE(cm->beginsc[v]) ? cm->beginsc[v] : IMPOSSIBLE);
          if (! NOT_IMPOSSIBLE(bsc)) continue;
          for (int j = 0; j <= diag_L; j++) {
            for (int d = 0; d <= j; d++) {
              if (NOT_IMPOSSIBLE(diag_mx->dp[v][j][d])) {
                float sc = bsc + diag_mx->dp[v][j][d];
                if (sc > align_max) {
                  align_max = sc;
                  best_v = v;
                  best_j = j;
                  best_d = d;
                }
              }
            }
          }
        }

        /* Method 4: best scan hit from FastIInsideScan */
        CM_TOPHITS *diag_th = cm_tophits_Create ();
        FastIInsideScan (cm, errbuf, cm->smx, SMX_NOQDB,
                         diag_dsq, 1, diag_L, -1000.0, diag_th,
                         FALSE, 0., NULL, NULL, NULL, NULL, NULL,
                         -1, -1, -1, NULL, NULL);
        float scan_best = IMPOSSIBLE;
        for (int h = 0; h < (int) diag_th->N; h++)
          if (diag_th->unsrt[h].score > scan_best) scan_best = diag_th->unsrt[h].score;
        cm_tophits_Destroy (diag_th);

        fprintf (diag_fp, "%.4f %.4f %.4f %.4f %d %d %d\n",
                 align_dL, align_maxLL, align_max, scan_best, best_v, best_j, best_d);
        if (di < 10)
          printf ("  %4d  align_dL=%8.3f  align_maxLL=%8.3f  align_max=%8.3f  scan_best=%8.3f\n",
                  di, align_dL, align_maxLL, align_max, scan_best);
        free (diag_dsq);
      }
      fclose (diag_fp);
      cm_mx_Destroy (diag_mx);
      printf ("  (full data in /tmp/diag_scores.txt)\n\n");
    }

    /* --imcmc mode: MCMC sampling from the tail of the null distribution.
     * No IS weights — samples are from the correct distribution by construction. */
    if (esl_opt_GetBoolean (go, "--imcmc")) {
      float mcmc_mu     = (float) esl_opt_GetReal (go, "--imcmc-mu");
      int   n_chains    = esl_opt_GetInteger (go, "--imcmc-chains");
      int   n_steps     = esl_opt_GetInteger (go, "--imcmc-steps");
      int   n_burnin    = esl_opt_GetInteger (go, "--imcmc-burnin");
      int   verbose     = esl_opt_GetBoolean (go, "-v");
      float *mcmc_scores = NULL;
      int    mcmc_N = 0;

      int use_maxv = esl_opt_GetBoolean (go, "--imcmc-maxv");
      if ((status = cm_MCMC_tail (cm, emit_cm, cfg, go, errbuf,
                                   mcmc_mu, n_chains, n_steps, n_burnin,
                                   use_maxv, verbose,
                                   &mcmc_scores, &mcmc_N)) != eslOK)
        cm_Fail (errbuf);

      printf ("MCMC: collected %d scores (mu=%.3f, %d chains x %d steps, %d burnin)\n",
              mcmc_N, mcmc_mu, n_chains, n_steps, n_burnin);

      /* Fit lambda — unweighted MLE */
      if (mcmc_N >= 10) {
        double sum_diff = 0.;
        for (i = 0; i < mcmc_N; i++) sum_diff += mcmc_scores[i] - mcmc_mu;
        double mcmc_lambda = mcmc_N / sum_diff;
        printf ("  lambda = %.4f  (N=%d)\n", mcmc_lambda, mcmc_N);
      }

      /* Write to score file if requested */
      if (cfg->isscfp != NULL && mcmc_scores != NULL) {
        for (i = 0; i < mcmc_N; i++)
          fprintf (cfg->isscfp, "%.4f\t1.000000\n", mcmc_scores[i]);
      }

      if (mcmc_scores != NULL) free (mcmc_scores);
      printf ("\n");
    }

    /* Search CM-sampled sequences and collect score histograms.
     * emit_cm (if non-NULL) is the proposal distribution (exponentiated CM);
     * cm is the original CM used for searching.
     * Diagnostic: if --no-weight, use emit_cm for searching too, so both
     * random and IS sequences use the same model. */
    {
      CM_t *is_search_cm
          = (esl_opt_GetBoolean (go, "--no-weight") && emit_cm != NULL) ? emit_cm : cm;
      if ((status = collect_scores (go, cfg, errbuf, is_search_cm, emit_cm, cfg->sN, -1, &sscN,
                                    &sscA, &swtA, &is_dbsize))
          != eslOK)
        cm_Fail (errbuf); /* the -1 passed as L tells collect_scores to sample from the CM */
    }

    /* Display weight and ESS statistics */
    if (swtA != NULL) {
      min_wt = max_wt = swtA[0];
      sum_wt = sum_wt_sq = 0.;
      for (i = 0; i < sscN; i++) {
        min_wt = ESL_MIN (min_wt, swtA[i]);
        max_wt = ESL_MAX (max_wt, swtA[i]);
        sum_wt += swtA[i];
        sum_wt_sq += swtA[i] * swtA[i];
      }
      ess = (sum_wt_sq > 0.) ? (sum_wt * sum_wt) / sum_wt_sq : 0.;
      printf ("Importance sampling weight statistics (N=%d hits from %d seqs):\n", sscN, cfg->sN);
      printf ("\t%12s: %12.6f\n", "min_weight", min_wt);
      printf ("\t%12s: %12.6f\n", "max_weight", max_wt);
      printf ("\t%12s: %12.6f\n", "sum_weight", sum_wt);
      printf ("\t%12s: %12.1f\n", "ESS", ess);
      printf ("\t%12s: %12.4f\n", "ESS/N", (sscN > 0) ? ess / (double)sscN : 0.);
      printf ("\n");
    }

    /* Write IS scores and weights to --isscfile, if requested */
    if (cfg->isscfp != NULL && sscA != NULL && swtA != NULL) {
      for (i = 0; i < sscN; i++)
        fprintf (cfg->isscfp, "%.4f\t%.6f\n", sscA[i], swtA[i]);
    }

    tailp = esl_opt_IsOn (go, "--itailp") ? esl_opt_GetReal (go, "--itailp") : 0.;
    status = fit_histogram (go, cfg, errbuf, tailp, is_dbsize, TRUE, sscA, swtA, sscN, exp_mode,
                            &mu, &lambda, &nsamphits, &tailp);
    if (status == eslOK) {
      printf ("Impt sampled seq fit histogram:\n\t%12s: %9.5f\n\t%12s: %9.5f\n\t%12s: "
              "%9.5f\n\t%12s: %9.5f\n\n",
              "mu", mu, "lambda", lambda, "nsamphits", nsamphits, "tailp", tailp);
    } else {
      printf ("Impt sampled seq fit histogram: SKIPPED (too few points; N=%d)\n\n", sscN);
      status = eslOK; /* non-fatal when --isscfile is being used for post-processing */
    }

    /* Diagnostic: fraction of hits in 1-bit bins [ilo, ilo+10).
     * Helps determine if ilo is low enough to cover the tail without truncation.
     * A healthy distribution has few hits near ilo; pile-up near ilo signals truncation. */
    if (esl_opt_IsOn (go, "--ilo") && sscA != NULL && sscN > 0) {
      float ilo_val = esl_opt_GetReal (go, "--ilo");
      int nbins = 10;
      int counts[10] = { 0 };
      int b;
      for (i = 0; i < sscN; i++) {
        b = (int)floor (sscA[i] - ilo_val);
        if (b >= 0 && b < nbins)
          counts[b]++;
      }
      printf ("Search score proximity to ilo (%.1f): fraction of hits in 1-bit bins:\n", ilo_val);
      for (b = 0; b < nbins; b++) {
        printf ("\t[%5.1f,%5.1f): %.4f\n", ilo_val + b, ilo_val + b + 1,
                (float)counts[b] / (float)sscN);
      }
      printf ("\n");
    }
    if (nrandhits > 0.) {
      sc_tailp = ((float)nsamphits / (float)nrandhits);
      SetExpInfo (impt_expinfo, lambda, mu, (double)(cfg->rL * cfg->rN),
                  (int)nrandhits, /* actually this is scaled_nhits */
                  sc_tailp);
      debug_print_expinfo (impt_expinfo);
    }

    /* output to --ifile, if nec */
    if (cfg->ifp != NULL && nrandhits > 0.) {
      tailp = esl_opt_GetReal (go, "--imax");
      nfits = esl_opt_GetInteger (go, "--infit");
      tailp_step = (tailp - esl_opt_GetReal (go, "--imin")) / (float)nfits;
      for (i = 0; i < nfits; i++) {
        if ((status = fit_histogram (go, cfg, errbuf, tailp, 0., TRUE, sscA, swtA, sscN, exp_mode,
                                     &mu, &lambda, &nsamphits, NULL))
            != eslOK)
          cm_Fail (errbuf);
        sc_tailp = ((float)nsamphits / (float)nrandhits);
        fprintf (cfg->ifp, "%g  %g  %g  %g  %g  %g  %g\n", tailp, sc_tailp, lambda,
                 (mu - log (1. / sc_tailp) / lambda), mu, nsamphits, nrandhits);
        tailp -= tailp_step;
      }
    }

    /* Matched-length random sequence control: generate sN random
     * sequences of length clen (same as typical emitted parsetree
     * length), search and fit — apples-to-apples comparison with
     * importance sampling results.
     * Skip when --no-rand is set (reference distribution already known from --refN). */
    if (!esl_opt_GetBoolean (go, "--no-rand")) {
      if ((status = collect_scores (go, cfg, errbuf, cm, NULL, cfg->sN, cm->clen, &mscN, &mscA,
                                    &mwtA, &match_dbsize))
          != eslOK)
        cm_Fail (errbuf);
      tailp = esl_opt_IsOn (go, "--rtailp") ? esl_opt_GetReal (go, "--rtailp") : 0.;
      if ((status = fit_histogram (go, cfg, errbuf, tailp, match_dbsize, FALSE, mscA, NULL, mscN,
                                   exp_mode, &mu, &lambda, &nmatchhits, &tailp))
          != eslOK)
        cm_Fail (errbuf);
      printf ("Matched-length (L=%d) random seq fit histogram:\n\t%12s: %9.5f\n\t%12s: "
              "%9.5f\n\t%12s: %9.5f\n\t%12s: %9.5f\n\n",
              cm->clen, "mu", mu, "lambda", lambda, "nhits", nmatchhits, "tailp", tailp);
      SetExpInfo (match_expinfo, lambda, mu, (double)(cm->clen * cfg->sN), (int)nmatchhits, tailp);
      debug_print_expinfo (match_expinfo);
    } /* end if(!--no-rand) */

    if (mscA != NULL)
      free (mscA);
    if (rscA != NULL)
      free (rscA);
    if (sscA != NULL)
      free (sscA);
    if (swtA != NULL)
      free (swtA);
    mscA = rscA = sscA = NULL;
    swtA = NULL;
    mscN = rscN = sscN = 0;

    if (emit_cm != NULL)
      FreeCM (emit_cm);
    FreeCM (cm);
  }

  if (status != eslEOF)
    cm_Fail (cfg->cmfp->errbuf);
  return;

ERROR:
  cm_Fail ("Out of memory.");
  return;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults.
 * Follows cmcalibrate.c:initialize_cm() pattern.
 */
static int
initialize_cm (const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, int do_local,
               char *errbuf) {
  int status;

  /* config QDB? yes unless --noqdb enabled */
  if (esl_opt_GetBoolean (go, "--noqdb")) {
    cm->search_opts |= CM_SEARCH_NONBANDED; /* don't use QDB to search */
  } else {
    cm->search_opts |= CM_SEARCH_QDB; /* use QDB to search */
    if (CheckCMQDBInfo (cm->qdbinfo, 0., FALSE, esl_opt_GetReal (go, "--beta"), TRUE) != eslOK) {
      cm->config_opts |= CM_CONFIG_QDB; /* configure QDB */
      cm->qdbinfo->beta1 = esl_opt_GetReal (go, "--beta");
      cm->qdbinfo->beta2 = esl_opt_GetReal (go, "--beta");
    }
  }

  cm->search_opts |= CM_SEARCH_NOALIGN;

  if (esl_opt_GetBoolean (go, "--null3"))
    cm->search_opts |= CM_SEARCH_NULL3;

  /* ALWAYS use the greedy overlap resolution algorithm to return hits for exp calculation
   * it's irrelevant for filter threshold stats, we return best score per seq for that */
  /* don't turn on CM_SEARCH_CMNOTGREEDY */

  if (do_local) {
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }

  /* Note: --exp exponentiation is handled in master() on a separate
   * emit_cm clone, not here. This CM is the original for searching. */

  /* we'll need a scan matrix */
  cm->config_opts |= CM_CONFIG_SCANMX;

  /* configure */
  if ((status = cm_Configure (cm, errbuf, -1)) != eslOK)
    return status;

  if (cm->smx == NULL)
    ESL_FAIL (eslEINVAL, errbuf, "unable to create scan matrix for CM");

  /* Set emit flags AFTER cm_Configure() (cm_nonconfigured_Verify()
   * requires these to be down on an unconfigured CM) */
  if (!esl_opt_GetBoolean (go, "--ilocal")) {
    cm->flags |= CM_EMIT_NO_LOCAL_BEGINS;
    cm->flags |= CM_EMIT_NO_LOCAL_ENDS;
  }

  return eslOK;
}

/* Function: print_run_info
 * Date:     EPN, Mon Mar  3 06:01:13 2008
 *
 * Purpose:  Print information on this run of cmsim.
 *           Command used to run it, and execution date.
 *
 * Returns:  eslOK on success
 */
static int
print_run_info (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf) {
  int status;
  char *command;
  char *date;

  if ((status = get_command (go, errbuf, &command)) != eslOK)
    return status;
  if ((status = GetDate (errbuf, &date)) != eslOK)
    return status;

  fprintf (stdout, "%-10s %s\n", "# command:", command);
  fprintf (stdout, "%-10s %s\n", "# date:", date);
  fprintf (stdout, "%-10s %" PRIu32 "\n", "# seed:", esl_randomness_GetSeed (cfg->r));

  fprintf (stdout, "#\n");
  free (command);
  free (date);
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
get_command (const ESL_GETOPTS *go, char *errbuf, char **ret_command) {
  int status;
  int i;
  char *command = NULL;

  for (i = 0; i < go->argc; i++) { /* copy all command line options and args */
    if ((status = esl_strcat (&(command), -1, go->argv[i], -1)) != eslOK)
      goto ERROR;
    if (i < (go->argc - 1))
      if ((status = esl_strcat (&(command), -1, " ", 1)) != eslOK)
        goto ERROR;
  }
  *ret_command = command;

  return eslOK;

ERROR:
  ESL_FAIL (status, errbuf, "get_command(): memory allocation error.");
  return status;
}

/* collect_scores()
 *
 * Generate and score sequences with a CM.
 *
 * Two different modes:
 * 1. generate from the CM.
 * 2. generated random sequences as either
 *    25% ACGU (default) or from a hard-coded HMM
 *    that generates genome-like sequences (if
 *    --rhmm).
 *
 * Return scores in <ret_scA> and number of scores in
 * <ret_scN>. If sampling from CM, also return importance
 * weights in <ret_wtA>.
 */
static int
collect_scores (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, CM_t *emit_cm,
                int N, int L, int *ret_scN, float **ret_scA, float **ret_wtA, double *ret_dbsize) {
  int status;
  int scN = 0;           /* number of hits reported thus far, for all seqs */
  float *scA = NULL;     /* [0..rscN-1] hit scores for all seqs */
  float *wtA = NULL;     /* [0..rscN-1] importance weights for all seqs (if sampling from CM) */
  ESL_DSQ *dsq = NULL;   /* digitized sequence to search */
  int i, h;              /* counters */
  int do_sample;         /* TRUE to sample from the CM, FALSE to sample random seqs */
  void *tmp;             /* ptr for ESL_RALLOC */
  CM_TOPHITS *th = NULL; /* hit list from search */
  int use_qdbs;          /* are we using QDBs? */
  float cutoff;          /* minimum score to report from scan functions */
  Parsetree_t *tr = NULL;
  float parsetree_sc; /* parsetree score from CM (for importance weight) */
  float search_cm_parsetree_sc; /* parsetree score under search CM (same parsetree, different scores) */
  float cp9_trace_sc;           /* CP9 trace score (parsetree mapped to HMM path) */
  float cp9_fwd_sc;             /* CP9 Forward score (sum over all HMM paths) */
  float weight;       /* importance weight = 2^(-parsetree_sc) */
  float ilo, ihi;     /* parsetree score range for rejection sampling */
  int do_filter;      /* TRUE if --ilo or --ihi is set */
  int do_isubtr;              /* TRUE if --isubtr: sub-parsetree IS mode */
  int do_iflank;              /* TRUE if --iflank: scan flanks only, add qc_sc as direct hit */
  int do_ipaint;              /* TRUE if --ipaint: painting mode */
  CM_SCAN_MX *smx_flank = NULL; /* reduced-W scan matrix for flank scans (--iflank-W) */
  float imu;                  /* minimum candidate score (--imu) for sub-parsetree IS */
  int   isubtr_best_v;        /* v* state index for current sequence (do_isubtr) */
  int   isubtr_il;            /* emitl of chosen sub-parsetree */
  int   isubtr_ir;            /* emitr of chosen sub-parsetree */
  float isubtr_best_candidate;/* candidate_sc of chosen sub-parsetree (fallback weight) */
  int   n_isubtr_hit_found;   /* sequences where th had a hit with root==v* */
  int   n_isubtr_hit_missing; /* sequences where th had NO hit with root==v* */
  int n_emitted;      /* total number of parsetrees emitted (including rejected) */
  int n_rejected;     /* number of parsetrees rejected by score filter */
  double dbsize = 0.; /* effective searched nt: sum(w_i*L_i) for IS, N*L for random */

  /* the HMM that generates sequences for exponential tail fitting */
  int ghmm_nstates = 0;     /* number of states in the HMM */
  double *ghmm_sA = NULL;   /* start probabilities [0..ghmm_nstates-1] */
  double **ghmm_tAA = NULL; /* transition probabilities [0..nstates-1][0..nstates-1] */
  double **ghmm_eAA = NULL; /* emission probabilities   [0..nstates-1][0..abc->K-1] */

  do_sample = (L == -1) ? TRUE : FALSE;

  /* Set up parsetree score filtering for rejection sampling */
  do_filter  = (esl_opt_IsOn (go, "--ilo") || esl_opt_IsOn (go, "--ihi")) ? TRUE : FALSE;
  ilo        = esl_opt_IsOn (go, "--ilo") ? esl_opt_GetReal (go, "--ilo") : -eslINFINITY;
  ihi        = esl_opt_IsOn (go, "--ihi") ? esl_opt_GetReal (go, "--ihi") : eslINFINITY;
  do_isubtr  = esl_opt_GetBoolean (go, "--isubtr");
  do_iflank  = esl_opt_GetBoolean (go, "--iflank");
  do_ipaint  = esl_opt_GetBoolean (go, "--ipaint");
  imu        = esl_opt_IsOn (go, "--imu") ? (float) esl_opt_GetReal (go, "--imu") : 0.;
  n_emitted             = 0;
  n_rejected            = 0;
  n_isubtr_hit_found    = 0;
  n_isubtr_hit_missing  = 0;
  isubtr_best_v         = -1;
  isubtr_il             = -1;
  isubtr_ir             = -1;
  isubtr_best_candidate = 0.;

  use_qdbs = (cm->search_opts & CM_SEARCH_QDB) ? TRUE : FALSE;
  cutoff = -eslINFINITY; /* collect all hits */

  /* Create reduced-W scan matrix for flank scans if --iflank-W is set.
   * FastIInsideScan uses integer matrices (ialpha/ialpha_begl), so do_int=TRUE.
   * We temporarily reduce cm->W to W_cap, create the smx, then restore cm->W. */
  if (do_iflank) {
    int iflank_W = esl_opt_GetInteger(go, "--iflank-W");
    if (iflank_W > 0 && iflank_W < cm->smx->W) {
      int saved_W = cm->W;
      cm->W = iflank_W;
      if ((status = cm_scan_mx_Create(cm, errbuf, FALSE, TRUE, &smx_flank)) != eslOK)
        cm_Fail(errbuf);
      cm->W = saved_W;
      if (esl_opt_GetBoolean(go, "-v"))
        printf("# --iflank-W: created reduced-W scan matrix W=%d (model W=%d, %.2f Mb)\n",
               iflank_W, saved_W, smx_flank->size_Mb);
    }
  }

  /* get HMM for generating random seqs, if nec */
  if (esl_opt_GetBoolean (go, "--rhmm")) {
    if ((status = CreateGenomicHMM (cm->abc, errbuf, &ghmm_sA, &ghmm_tAA, &ghmm_eAA, &ghmm_nstates))
        != eslOK)
      cm_Fail ("ERROR unable to make HMM for generating random seqs");
  }

  /* Search sequences and collect score histograms */

  scN = 0;

  /* ================================================================
   * --ipaint: painting mode.
   * Pack multiple sub-scanned IS sequences into ~10Kb mega-sequences,
   * then scan each mega-sequence once.  Each kept sequence contributes
   * a guaranteed hit at its qc_sc score; flank regions (randomized)
   * contribute additional hits.
   * ================================================================ */
  if (do_ipaint && do_isubtr && do_sample) {
    int   ipaint_target_L = esl_opt_GetInteger (go, "--ipaint-L");
    float imutol          = (float) esl_opt_GetReal (go, "--imutol");
    float ipaint_lo       = esl_opt_IsOn (go, "--ipaint-lo") ? (float) esl_opt_GetReal (go, "--ipaint-lo") : -eslINFINITY;
    float ipaint_hi       = esl_opt_IsOn (go, "--ipaint-hi") ? (float) esl_opt_GetReal (go, "--ipaint-hi") :  eslINFINITY;
    int   do_ipaint_filter = (esl_opt_IsOn (go, "--ipaint-lo") || esl_opt_IsOn (go, "--ipaint-hi"));
    int   n_chunks        = 0;
    int   n_subscan       = 0;   /* total sub-scans performed */
    int   n_paint_kept    = 0;   /* total sequences kept across all chunks */
    int   n_paint_reject  = 0;   /* total sequences rejected (Viterbi or Inside filter) */
    int   n_qcsc_reject   = 0;   /* rejected by qc_sc bin filter (after sub-scan) */

    /* Per-chunk buffers (reused each chunk) */
    int    seg_alloc = (ipaint_target_L / 50) + 16;   /* generous initial alloc */
    int    n_seg     = 0;
    int    mega_L    = 0;
    ESL_DSQ **seg_dsq  = NULL;
    int      *seg_L    = NULL;
    int      *seg_il   = NULL;
    int      *seg_ir   = NULL;
    int      *seg_v    = NULL;
    float    *seg_qcsc = NULL;
    float    *seg_wt   = NULL;

    ESL_ALLOC (seg_dsq,  sizeof (ESL_DSQ *) * seg_alloc);
    ESL_ALLOC (seg_L,    sizeof (int)       * seg_alloc);
    ESL_ALLOC (seg_il,   sizeof (int)       * seg_alloc);
    ESL_ALLOC (seg_ir,   sizeof (int)       * seg_alloc);
    ESL_ALLOC (seg_v,    sizeof (int)       * seg_alloc);
    ESL_ALLOC (seg_qcsc, sizeof (float)     * seg_alloc);
    ESL_ALLOC (seg_wt,   sizeof (float)     * seg_alloc);

    for (i = 0; i < N; i++) {
      /* --- Phase 1: emit, find v*, sub-scan, accept/reject --- */
      if ((status = sample_sequence_from_cm (cfg, errbuf, (emit_cm != NULL) ? emit_cm : cm,
                                             (emit_cm != NULL) ? emit_cm : cm, &L, &dsq, &tr,
                                             &parsetree_sc))
          != eslOK)
        cm_Fail (errbuf);
      n_emitted++;

      /* Find v* via Viterbi sub-parsetree scores (same as --isubtr) */
      {
        float *subtree_sc = NULL;
        int    best_tidx  = -1;
        float  best_above = eslINFINITY;
        float  best_subtree_sc_val = 0.;
        float  best_cand  = 0.;
        int    tidx, il, ir;

        if ((status = ParsetreeSubtreeScores (cm, errbuf, tr, dsq, &subtree_sc)) != eslOK)
          cm_Fail (errbuf);

        for (tidx = 0; tidx < tr->n; tidx++) {
          int v_t = tr->state[tidx];
          if (v_t == cm->M) continue;
          if (cm->sttype[v_t] == E_st) continue;
          if (cm->sttype[v_t] == B_st) continue;
          if (NOT_IMPOSSIBLE (cm->beginsc[v_t])) {
            float candidate = cm->beginsc[v_t] + subtree_sc[tidx];
            float above     = candidate - imu;
            if (above >= 0. && above < imutol && above < best_above) {
              best_above          = above;
              best_tidx           = tidx;
              best_subtree_sc_val = subtree_sc[tidx];
              best_cand           = candidate;
            }
          }
        }
        free (subtree_sc);

        if (best_tidx == -1) {
          /* No qualifying sub-parsetree; reject */
          n_rejected++;
          n_paint_reject++;
          free (dsq); dsq = NULL;
          FreeParsetree (tr); tr = NULL;
          i--;
          continue;
        }

        isubtr_best_v         = tr->state[best_tidx];
        isubtr_il             = tr->emitl[best_tidx];
        isubtr_ir             = tr->emitr[best_tidx];
        isubtr_best_candidate = best_cand;
      }

      /* Sub-scan [il..ir] to get true Inside score (qc_sc) */
      {
        int L_v = isubtr_ir - isubtr_il + 1;
        float qc_sc = IMPOSSIBLE;
        float qc_sc_sumv = IMPOSSIBLE;
        int   do_sumv = esl_opt_GetBoolean (go, "--ipaint-sumv");

        /* --ipaint-allrand: randomize [il..ir] before sub-scan */
        if (esl_opt_GetBoolean (go, "--ipaint-allrand")) {
          int j;
          for (j = isubtr_il; j <= isubtr_ir; j++)
            dsq[j] = esl_rnd_FChoose (cfg->r, cm->null, cm->abc->K);
        }

        CM_TOPHITS *th_sub = cm_tophits_Create ();
        if (th_sub == NULL) ESL_FAIL (eslEMEM, errbuf, "out of memory");
        if ((status = FastIInsideScan (cm, errbuf, cm->smx, SMX_NOQDB,
                                       dsq, (int64_t) isubtr_il, (int64_t) isubtr_ir, cutoff,
                                       th_sub, cm->search_opts & CM_SEARCH_NULL3, 0.,
                                       NULL, NULL, NULL, NULL, NULL,
                                       isubtr_best_v, (int64_t) isubtr_ir, L_v, &qc_sc,
                                       do_sumv ? &qc_sc_sumv : NULL))
            != eslOK)
          cm_Fail (errbuf);
        cm_tophits_Destroy (th_sub);
        n_subscan++;

        if (qc_sc == IMPOSSIBLE)
          qc_sc = isubtr_best_candidate;

        /* --ipaint-sumv: use sum-over-all-v score instead of single-v* score */
        if (do_sumv && qc_sc_sumv != IMPOSSIBLE)
          qc_sc = qc_sc_sumv;

        weight = esl_opt_GetBoolean (go, "--no-weight") ? 1.0 : (float) pow (2.0, -qc_sc);

        if (esl_opt_GetBoolean (go, "-v"))
          printf ("  PAINT sub-scan %5d: v*=%d [%d..%d] L_v=%d qc_sc=%.3f%s cand=%.3f wt=%.4g\n",
                  i, isubtr_best_v, isubtr_il, isubtr_ir, L_v, qc_sc,
                  do_sumv ? "(sumv)" : "", isubtr_best_candidate, weight);

        /* Filter by target qc_sc bin if --ipaint-lo/hi are set */
        if (do_ipaint_filter && (qc_sc < ipaint_lo || qc_sc >= ipaint_hi)) {
          n_qcsc_reject++;
          free (dsq); dsq = NULL;
          FreeParsetree (tr); tr = NULL;
          i--;  /* retry this slot */
          continue;
        }

        n_paint_kept++;

        /* --ipaint-qcsconly: just collect the qc_sc hit directly, skip mega-seq */
        if (esl_opt_GetBoolean (go, "--ipaint-qcsconly")) {
          if (scN == 0) {
            ESL_ALLOC  (scA, sizeof (float) * (scN + 1));
            ESL_ALLOC  (wtA, sizeof (float) * (scN + 1));
          } else {
            ESL_RALLOC (scA, tmp, sizeof (float) * (scN + 1));
            ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + 1));
          }
          scA[scN] = qc_sc;
          wtA[scN] = weight;
          scN++;
          dbsize += (double) (isubtr_ir - isubtr_il + 1);
          free (dsq); dsq = NULL;
          FreeParsetree (tr); tr = NULL;
          continue;
        }

        /* Grow buffers if needed */
        if (n_seg >= seg_alloc) {
          seg_alloc *= 2;
          ESL_RALLOC (seg_dsq,  tmp, sizeof (ESL_DSQ *) * seg_alloc);
          ESL_RALLOC (seg_L,    tmp, sizeof (int)       * seg_alloc);
          ESL_RALLOC (seg_il,   tmp, sizeof (int)       * seg_alloc);
          ESL_RALLOC (seg_ir,   tmp, sizeof (int)       * seg_alloc);
          ESL_RALLOC (seg_v,    tmp, sizeof (int)       * seg_alloc);
          ESL_RALLOC (seg_qcsc, tmp, sizeof (float)     * seg_alloc);
          ESL_RALLOC (seg_wt,   tmp, sizeof (float)     * seg_alloc);
        }
        seg_dsq[n_seg]  = dsq;   dsq = NULL;   /* transfer ownership */
        seg_L[n_seg]    = L;
        seg_il[n_seg]   = isubtr_il;
        seg_ir[n_seg]   = isubtr_ir;
        seg_v[n_seg]    = isubtr_best_v;
        seg_qcsc[n_seg] = qc_sc;
        seg_wt[n_seg]   = weight;
        n_seg++;
        mega_L += L;

      }

      FreeParsetree (tr); tr = NULL;

      /* --- Phase 2: when mega-chunk is full, build + scan --- */
      if (mega_L >= ipaint_target_L || i == N - 1) {
        int   k, pos;
        int  *seg_start = NULL;  /* start position of each segment in mega-dsq (1-based) */
        ESL_DSQ *mega_dsq = NULL;

        ESL_ALLOC (seg_start, sizeof (int) * n_seg);
        ESL_ALLOC (mega_dsq,  sizeof (ESL_DSQ) * (mega_L + 2));
        mega_dsq[0] = eslDSQ_SENTINEL;

        /* Build the mega-sequence: randomize flanks, concatenate */
        pos = 1;
        for (k = 0; k < n_seg; k++) {
          int j;
          seg_start[k] = pos;
          /* Randomize flanks [1..il-1] and [ir+1..L] */
          for (j = 1; j < seg_il[k]; j++)
            seg_dsq[k][j] = esl_rnd_FChoose (cfg->r, cm->null, cm->abc->K);
          for (j = seg_ir[k] + 1; j <= seg_L[k]; j++)
            seg_dsq[k][j] = esl_rnd_FChoose (cfg->r, cm->null, cm->abc->K);
          /* --ipaint-allrand: also randomize CM region [il..ir] → pure random mega-seq */
          if (esl_opt_GetBoolean (go, "--ipaint-allrand")) {
            for (j = seg_il[k]; j <= seg_ir[k]; j++)
              seg_dsq[k][j] = esl_rnd_FChoose (cfg->r, cm->null, cm->abc->K);
          }
          /* Copy segment to mega-dsq (dsq is 1-based; copy residues 1..L) */
          memcpy (mega_dsq + pos, seg_dsq[k] + 1, sizeof (ESL_DSQ) * seg_L[k]);
          pos += seg_L[k];
        }
        mega_dsq[pos] = eslDSQ_SENTINEL;

        /* Scan the mega-sequence */
        th = cm_tophits_Create ();
        if (th == NULL) ESL_FAIL (eslEMEM, errbuf, "out of memory");
        if ((status = FastIInsideScan (cm, errbuf, cm->smx, use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB,
                                       mega_dsq, 1, mega_L, cutoff, th,
                                       cm->search_opts & CM_SEARCH_NULL3, 0.,
                                       NULL, NULL, NULL, NULL, NULL,
                                       -1, -1, -1, NULL, NULL))
            != eslOK)
          cm_Fail (errbuf);

        /* Collect hits: assign per-segment IS weights */
        if (th->N > 0) {
          if (scN == 0) {
            ESL_ALLOC  (scA, sizeof (float) * (scN + th->N));
            ESL_ALLOC  (wtA, sizeof (float) * (scN + th->N));
          } else {
            ESL_RALLOC (scA, tmp, sizeof (float) * (scN + th->N));
            ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + th->N));
          }
          for (h = 0; h < (int) th->N; h++) {
            /* Determine which segment this hit belongs to.
             * th->unsrt[h].start and .stop are 1-based positions in mega-dsq. */
            int hit_start = th->unsrt[h].start;
            int hit_stop  = th->unsrt[h].stop;
            int seg_idx   = -1;
            for (k = 0; k < n_seg; k++) {
              int seg_end = seg_start[k] + seg_L[k] - 1;
              if (hit_start >= seg_start[k] && hit_stop <= seg_end) {
                seg_idx = k;
                break;
              }
            }
            if (seg_idx == -1) continue;  /* hit spans segment boundary, discard */

            /* Weight depends on overlap with the CM region [il..ir].
             * CM region in mega-dsq coords:
             *   cm_start = seg_start[k] + (il - 1)
             *   cm_end   = seg_start[k] + (ir - 1)
             * Fractional IS weight: w = 2^(-qc_sc * overlap_frac / L_v*)
             * where overlap_frac = # positions in hit that are in [il..ir],
             * L_v* = ir - il + 1.  Gives w=1 for pure flank, w=2^(-qc_sc)
             * for full [il..ir] overlap, smooth interpolation between. */
            {
              int cm_start = seg_start[seg_idx] + (seg_il[seg_idx] - 1);
              int cm_end   = seg_start[seg_idx] + (seg_ir[seg_idx] - 1);
              int overlap  = 0;  /* # positions of hit overlapping CM region */
              float hit_wt;

              if (hit_stop >= cm_start && hit_start <= cm_end) {
                int ov_start = (hit_start > cm_start) ? hit_start : cm_start;
                int ov_end   = (hit_stop  < cm_end)   ? hit_stop  : cm_end;
                overlap = ov_end - ov_start + 1;
              }

              /* --ipaint-noqcsc: discard the specific v*-rooted hit at [il..ir] */
              if (esl_opt_GetBoolean (go, "--ipaint-noqcsc") &&
                  hit_start == cm_start && hit_stop == cm_end)
                continue;

              /* --ipaint-flankonly: discard any hit overlapping CM region */
              if (overlap > 0 && esl_opt_GetBoolean (go, "--ipaint-flankonly"))
                continue;

              if (esl_opt_GetBoolean (go, "--no-weight")) {
                hit_wt = 1.0;
              } else if (overlap > 0) {
                int L_v = seg_ir[seg_idx] - seg_il[seg_idx] + 1;
                double frac = (double) overlap / (double) L_v;
                hit_wt = (float) pow (2.0, -seg_qcsc[seg_idx] * frac);
              } else {
                hit_wt = 1.0;
              }
              scA[scN] = th->unsrt[h].score;
              wtA[scN] = hit_wt;
            }
            scN++;
          }
        }
        dbsize += (double) mega_L;

        n_chunks++;
        if (esl_opt_GetBoolean (go, "-v"))
          printf ("  PAINT chunk %d: %d seqs, mega_L=%d, %d hits collected\n",
                  n_chunks, n_seg, mega_L, (int) th->N);

        /* Free chunk resources */
        cm_tophits_Destroy (th); th = NULL;
        free (mega_dsq);
        free (seg_start);
        for (k = 0; k < n_seg; k++) free (seg_dsq[k]);
        n_seg  = 0;
        mega_L = 0;
      }
    } /* end for (i = 0; i < N; ...) */

    /* Report painting statistics */
    printf ("Painting: %d emitted, %d kept, %d viterbi-rejected, %d qcsc-rejected, %d sub-scans, %d chunks\n",
            n_emitted, n_paint_kept, n_paint_reject, n_qcsc_reject, n_subscan, n_chunks);
    if (do_ipaint_filter)
      printf ("  qc_sc filter: [%.1f, %.1f)  acceptance rate: %.2f%% (%d/%d sub-scans)\n",
              ipaint_lo, ipaint_hi,
              n_subscan > 0 ? 100.0 * n_paint_kept / (double)(n_paint_kept + n_qcsc_reject) : 0.,
              n_paint_kept, n_paint_kept + n_qcsc_reject);

    /* Free painting buffers */
    free (seg_dsq);
    free (seg_L);
    free (seg_il);
    free (seg_ir);
    free (seg_v);
    free (seg_qcsc);
    free (seg_wt);

  } else
  /* ================================================================ */

  for (i = 0; i < N; i++) {
    /* generate sequence */
    if (do_sample) {
      /* Emit from emit_cm (exponentiated proposal) if available, else from cm.
       * ParsetreeScore is computed by sample_sequence_from_cm using the
       * emission model, giving the correct importance weight. */
      if ((status = sample_sequence_from_cm (cfg, errbuf, (emit_cm != NULL) ? emit_cm : cm,
                                             (emit_cm != NULL) ? emit_cm : cm, &L, &dsq, &tr,
                                             &parsetree_sc))
          != eslOK)
        cm_Fail (errbuf);
      n_emitted++;

      /* Score the same parsetree under search_cm to get a cheap predictor
       * of search_cm Inside score. O(L), essentially free. */
      search_cm_parsetree_sc = 0.;
      if (emit_cm != NULL) {
        if ((status = ParsetreeScore (cm, NULL, errbuf, tr, dsq, FALSE, &search_cm_parsetree_sc,
                                      NULL, NULL, NULL, NULL)) != eslOK)
          cm_Fail (errbuf);
      } else {
        search_cm_parsetree_sc = parsetree_sc; /* no emit_cm, they're the same */
      }

      /* CP9 diagnostics: map parsetree to CP9 trace and score it, then
       * run CP9 Forward to get sum-over-all-paths score.
       * cp9_fwd - cp9_trace estimates the Inside-parsetree gap at HMM level.
       * Adding basepair_bonus (= search_cm_ptree - cp9_trace) back gives:
       *   cm_inside_est = cp9_fwd + basepair_bonus = cp9_fwd + search_cm_ptree - cp9_trace */
      cp9_trace_sc = 0.;
      cp9_fwd_sc = 0.;
      if (cm->cp9 != NULL && cm->cp9map != NULL) {
        CP9trace_t *cp9_tr = NULL;
        if (Parsetree2CP9trace (cm, tr, &cp9_tr) == eslOK) {
          cp9_trace_sc = CP9TraceScore (cm->cp9, dsq, cp9_tr);
          CP9FreeTrace (cp9_tr);
        }
        /* CP9 Forward: do_scan=FALSE, doing_align=TRUE, be_efficient=TRUE */
        if ((status = cp9_Forward (cm->cp9, errbuf, cm->cp9_mx, dsq, 1, L,
                                    FALSE, TRUE, TRUE, NULL, NULL, &cp9_fwd_sc)) != eslOK)
          cm_Fail (errbuf);
      }

      if (do_isubtr) {
        /* Sub-parsetree IS mode: find the subtree rooted at v* whose
         * candidate_sc = cm->beginsc[v*] + subtree_sc[v*] is in [imu, imu+imutol].
         * Among qualifying subtrees pick the one with the smallest candidate_sc
         * (closest to imu from above).  If none qualify, reject and re-emit.
         * Splice null random sequence into the flanking positions.
         * IS weight = 2^(-subtree_sc[v*]). */
        float *subtree_sc = NULL;
        int   best_tidx   = -1;
        float best_above  = eslINFINITY; /* candidate_sc - imu for best so far */
        float best_subtree_sc = 0.;
        float best_candidate  = 0.;
        float imutol = (float) esl_opt_GetReal (go, "--imutol");
        int   tidx, il, ir, j;

        if ((status = ParsetreeSubtreeScores (cm, errbuf, tr, dsq, &subtree_sc)) != eslOK)
          cm_Fail (errbuf);

        /* Find v* with minimum candidate_sc in [imu, imu + imutol]. */
        for (tidx = 0; tidx < tr->n; tidx++) {
          int v_t = tr->state[tidx];
          if (v_t == cm->M) continue;                   /* EL state */
          if (cm->sttype[v_t] == E_st) continue;
          if (cm->sttype[v_t] == B_st) continue;
          if (NOT_IMPOSSIBLE (cm->beginsc[v_t])) {
            float candidate = cm->beginsc[v_t] + subtree_sc[tidx];
            float above     = candidate - imu;
            if (above >= 0. && above < imutol && above < best_above) {
              best_above      = above;
              best_tidx       = tidx;
              best_subtree_sc = subtree_sc[tidx];
              best_candidate  = candidate;
            }
          }
        }
        free (subtree_sc);
        subtree_sc = NULL;

        if (best_tidx == -1) {
          /* No sub-parsetree in [imu, imu+imutol]; reject and re-emit. */
          n_rejected++;
          free (dsq);
          FreeParsetree (tr);
          tr  = NULL;
          dsq = NULL;
          i--;
          continue;
        }

        /* Splice: overwrite flanking positions with null random residues */
        il = tr->emitl[best_tidx];
        ir = tr->emitr[best_tidx];
        for (j = 1; j < il; j++)
          dsq[j] = esl_rnd_FChoose (cfg->r, cm->null, cm->abc->K);
        for (j = ir + 1; j <= L; j++)
          dsq[j] = esl_rnd_FChoose (cfg->r, cm->null, cm->abc->K);

        /* Save v*, [il..ir], candidate_sc for th lookup after FastIInsideScan */
        isubtr_best_v         = tr->state[best_tidx];
        isubtr_il             = il;
        isubtr_ir             = ir;
        isubtr_best_candidate = best_candidate;

        weight = esl_opt_GetBoolean (go, "--no-weight") ? 1.0 : (float) pow (2.0, -best_candidate);
        if (esl_opt_GetBoolean (go, "-v")) {
          printf ("SEQ %5d  L: %4d  il: %4d  ir: %4d  candidate_sc: %8.3f  "
                  "subtree_sc: %8.3f  weight: %12.6g  above_imu: %8.3f  (emitted: %d  rejected: %d)\n",
                  i, L, il, ir, best_candidate, best_subtree_sc, weight, best_above, n_emitted,
                  n_rejected);
        }
      } else {
        /* Standard IS: rejection sampling on full parsetree score, then weight. */

        /* Reject sequences with -inf parsetree score (impossible transitions) */
        if (parsetree_sc <= -eslINFINITY || !isfinite(parsetree_sc)) {
          n_rejected++;
          free (dsq);
          FreeParsetree (tr);
          i--;
          continue;
        }

        /* Rejection sampling: if parsetree score is outside [ilo, ihi], reject
         * and re-emit. Emission + scoring is fast (O(clen)); the expensive DP
         * search only runs on accepted sequences. */
        if (do_filter && (parsetree_sc < ilo || parsetree_sc > ihi)) {
          n_rejected++;
          free (dsq);
          FreeParsetree (tr);
          i--; /* retry this slot */
          continue;
        }

        /* compute importance weight: w = 2^(-parsetree_sc), or 1.0 if --no-weight */
        weight = esl_opt_GetBoolean (go, "--no-weight") ? 1.0 : (float) pow (2.0, -parsetree_sc);
        if (esl_opt_GetBoolean (go, "-v")) {
          printf ("SEQ %5d  L: %4d  parsetree_sc: %8.3f  search_cm_ptree_sc: %8.3f  weight: %12.6g  (emitted: %d "
                  "rejected: %d)\n",
                  i, L, parsetree_sc, search_cm_parsetree_sc, weight, n_emitted, n_rejected);
        }
      }
    } else { /* generate random sequence, either iid (25% ACGU) or from a 'genome-like' HMM */
      /* L was passed in as the desired sequence length (e.g. cm->clen for --refN, cfg->rL for
       * normal random seq mode). Use L directly; do NOT overwrite with cfg->rL. */
      if (esl_opt_GetBoolean (go, "--rhmm")) {
        if ((status = SampleGenomicSequenceFromHMM (cfg->r, cm->abc, errbuf, ghmm_sA, ghmm_tAA,
                                                    ghmm_eAA, ghmm_nstates, L, &dsq))
            != eslOK)
          cm_Fail (errbuf);
      } else {
        ESL_ALLOC (dsq, sizeof (ESL_DSQ) * (L + 2));
        if ((status = esl_rsq_xfIID (cfg->r, cm->null, cm->abc->K, L, dsq) != eslOK))
          cm_Fail ("ERROR, couldn't generate random sequence");
      }
    }

    /* Search the sequence; collect hits and IS weight.
     * --iflank (Idea 1): scan [il..ir] to extract qc_sc, add qc_sc as direct
     * CM-region hit, scan flanks [1..il-1] and [ir+1..L] for null hits.
     * Default: full scan [1..L] with CYK or Inside. */
    float isubtr_qc_sc = IMPOSSIBLE;

    if (do_iflank && do_isubtr && do_sample && isubtr_best_v != -1) {
      /* --iflank: two-flank approach.
       * Scan [il..ir] with FastIInsideScan to get true qc_sc = beginsc[v*] + Inside([il..ir], v*).
       * The query cell (v*, ir, L_cm) extracts this directly from the DP matrix.
       * Use qc_sc as CM-region hit score and to update IS weight. */
      {
        int L_cm = isubtr_ir - isubtr_il + 1;
        CM_TOPHITS *th_cm = cm_tophits_Create ();
        if (th_cm == NULL) ESL_FAIL (eslEMEM, errbuf, "out of memory");
        if ((status = FastIInsideScan (cm, errbuf, cm->smx, SMX_NOQDB,
                                       dsq, (int64_t) isubtr_il, (int64_t) isubtr_ir, cutoff,
                                       th_cm, cm->search_opts & CM_SEARCH_NULL3, 0.,
                                       NULL, NULL, NULL, NULL, NULL,
                                       isubtr_best_v, (int64_t) isubtr_ir, L_cm, &isubtr_qc_sc, NULL))
            != eslOK)
          cm_Fail (errbuf);
        cm_tophits_Destroy (th_cm);
      }
      if (isubtr_qc_sc == IMPOSSIBLE) /* fallback: shouldn't happen for valid subtree */
        isubtr_qc_sc = isubtr_best_candidate;

      /* Update IS weight from true qc_sc */
      if (!esl_opt_GetBoolean (go, "--no-weight"))
        weight = (float) pow (2.0, -isubtr_qc_sc);

      if (esl_opt_GetBoolean (go, "-v"))
        printf ("  INSIDE qc (iflank): v*=%d [%d..%d] qc_sc=%.3f  candidate_sc=%.3f  weight=%.6g\n",
                isubtr_best_v, isubtr_il, isubtr_ir,
                isubtr_qc_sc, isubtr_best_candidate, weight);
      n_isubtr_hit_found++;
      dbsize += (double) L;

      /* Add CM-region hit using true qc_sc as the hit score. */
      if (scN == 0) {
        ESL_ALLOC  (scA, sizeof (float) * (scN + 1));
        if (do_sample) ESL_ALLOC  (wtA, sizeof (float) * (scN + 1));
      } else {
        ESL_RALLOC (scA, tmp, sizeof (float) * (scN + 1));
        if (do_sample) ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + 1));
      }
      scA[scN] = isubtr_qc_sc;
      if (do_sample) wtA[scN] = weight;
      scN++;

      /* (d) Scan left flank [1..il-1] for null-background hits. */
      if (isubtr_il > 1) {
        CM_TOPHITS *th_flank = cm_tophits_Create ();
        if (th_flank == NULL) ESL_FAIL (eslEMEM, errbuf, "out of memory");
        if (cm->search_opts & CM_SEARCH_INSIDE) {
          CM_SCAN_MX *smx_use = (smx_flank != NULL) ? smx_flank : cm->smx;
          if ((status = FastIInsideScan (cm, errbuf, smx_use, SMX_NOQDB,
                                         dsq, 1, (int64_t) (isubtr_il - 1), cutoff,
                                         th_flank, cm->search_opts & CM_SEARCH_NULL3, 0.,
                                         NULL, NULL, NULL, NULL, NULL,
                                         -1, -1, -1, NULL, NULL))
              != eslOK)
            cm_Fail (errbuf);
        }
        if (th_flank->N > 0) {
          if (scN == 0) {
            ESL_ALLOC  (scA, sizeof (float) * (scN + th_flank->N));
            if (do_sample) ESL_ALLOC  (wtA, sizeof (float) * (scN + th_flank->N));
          } else {
            ESL_RALLOC (scA, tmp, sizeof (float) * (scN + th_flank->N));
            if (do_sample) ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + th_flank->N));
          }
          for (h = 0; h < (int) th_flank->N; h++) {
            scA[scN + h] = th_flank->unsrt[h].score;
            if (do_sample) wtA[scN + h] = weight;
          }
          scN += th_flank->N;
        }
        cm_tophits_Destroy (th_flank);
      }

      /* (e) Scan right flank [ir+1..L] for null-background hits. */
      if (isubtr_ir < L) {
        CM_TOPHITS *th_flank = cm_tophits_Create ();
        if (th_flank == NULL) ESL_FAIL (eslEMEM, errbuf, "out of memory");
        if (cm->search_opts & CM_SEARCH_INSIDE) {
          CM_SCAN_MX *smx_use = (smx_flank != NULL) ? smx_flank : cm->smx;
          if ((status = FastIInsideScan (cm, errbuf, smx_use, SMX_NOQDB,
                                         dsq, (int64_t) (isubtr_ir + 1), L, cutoff,
                                         th_flank, cm->search_opts & CM_SEARCH_NULL3, 0.,
                                         NULL, NULL, NULL, NULL, NULL,
                                         -1, -1, -1, NULL, NULL))
              != eslOK)
            cm_Fail (errbuf);
        }
        if (th_flank->N > 0) {
          if (scN == 0) {
            ESL_ALLOC  (scA, sizeof (float) * (scN + th_flank->N));
            if (do_sample) ESL_ALLOC  (wtA, sizeof (float) * (scN + th_flank->N));
          } else {
            ESL_RALLOC (scA, tmp, sizeof (float) * (scN + th_flank->N));
            if (do_sample) ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + th_flank->N));
          }
          for (h = 0; h < (int) th_flank->N; h++) {
            scA[scN + h] = th_flank->unsrt[h].score;
            if (do_sample) wtA[scN + h] = weight;
          }
          scN += th_flank->N;
        }
        cm_tophits_Destroy (th_flank);
      }

    } else {
      /* Default: full-sequence scan with CYK or Inside. */
      th = cm_tophits_Create ();
      if (th == NULL)
        ESL_FAIL (eslEMEM, errbuf, "out of memory");

      /* Debug: print sequence before scan if verbose */
      if (esl_opt_GetBoolean (go, "-v") && L <= 200) {
        printf ("  DSQ[1..%d]: ", L);
        int p;
        for (p = 1; p <= L; p++) printf ("%c", cm->abc->sym[dsq[p]]);
        printf ("\n");
      }

      /* --iinside-wt with --ewt-lo/--ewt-hi: run emit_cm Inside FIRST to get
       * the IS weight score. Reject sequences outside [ewt-lo, ewt-hi] BEFORE
       * the expensive search_cm Inside scan. This saves one full Inside scan
       * for each rejected sequence. */
      float iinside_wt_sc = IMPOSSIBLE;  /* emit_cm Inside score for IS weight */
      if (esl_opt_GetBoolean (go, "--iinside-wt") && do_sample &&
          ! esl_opt_GetBoolean (go, "--no-weight") && emit_cm != NULL) {
        CM_t *wt_cm = emit_cm;
        CM_TOPHITS *th_wt = cm_tophits_Create ();
        if (th_wt == NULL) ESL_FAIL (eslEMEM, errbuf, "out of memory");

        /* Run emit_cm Inside scan to get proposal probability.
         * Use HMM banding if --ihbanded (emit_cm seqs are CM-like → tight bands).
         * Otherwise use unbanded FastIInsideScan with v=0 query cell for glocal. */
        float wt_qc_sc = IMPOSSIBLE;
        if (esl_opt_GetBoolean (go, "--ihbanded")) {
          float hb_mxsize = esl_opt_GetReal (go, "--mxsize");
          float hb_Mb;
          double save_tau = wt_cm->tau;
          wt_cm->tau = esl_opt_GetReal (go, "--tau");
          if ((status = cp9_Seq2Bands (wt_cm, errbuf, wt_cm->cp9_mx, wt_cm->cp9_bmx, wt_cm->cp9_bmx,
                                        dsq, 1, L, wt_cm->cp9b, TRUE, PLI_PASS_STD_ANY, 0))
              != eslOK)
            cm_Fail (errbuf);
          if ((status = FastFInsideScanHB (wt_cm, errbuf, wt_cm->hb_mx, hb_mxsize,
                                            dsq, 1, L, cutoff, th_wt,
                                            wt_cm->search_opts & CM_SEARCH_NULL3, 0.,
                                            NULL, NULL, NULL))
              != eslOK)
            cm_Fail (errbuf);
          wt_cm->tau = save_tau;
          /* Use best hit from banded scan */
          for (h = 0; h < (int) th_wt->N; h++)
            if (th_wt->unsrt[h].score > iinside_wt_sc)
              iinside_wt_sc = th_wt->unsrt[h].score;
        } else {
          if ((status = FastIInsideScan (wt_cm, errbuf, wt_cm->smx,
                                         use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB,
                                         dsq, 1, L, cutoff, th_wt,
                                         wt_cm->search_opts & CM_SEARCH_NULL3, 0.,
                                         NULL, NULL, NULL, NULL, NULL,
                                         esl_opt_GetBoolean (go, "--glocal") ? 0         : -1,
                                         esl_opt_GetBoolean (go, "--glocal") ? (int64_t)L : -1,
                                         esl_opt_GetBoolean (go, "--glocal") ? L          : -1,
                                         esl_opt_GetBoolean (go, "--glocal") ? &wt_qc_sc  : NULL,
                                         NULL))
              != eslOK)
            cm_Fail (errbuf);

          /* Get emit_cm score: v=0 query cell for glocal, best hit for local */
          if (esl_opt_GetBoolean (go, "--glocal") && wt_qc_sc != IMPOSSIBLE) {
            iinside_wt_sc = wt_qc_sc;
          } else {
            for (h = 0; h < (int) th_wt->N; h++)
              if (th_wt->unsrt[h].score > iinside_wt_sc)
                iinside_wt_sc = th_wt->unsrt[h].score;
          }
        }
        cm_tophits_Destroy (th_wt);

        /* Reject based on emit_cm Inside score floor/ceiling */
        {
          float ewt_lo = esl_opt_IsOn (go, "--ewt-lo") ? (float) esl_opt_GetReal (go, "--ewt-lo") : -eslINFINITY;
          float ewt_hi = esl_opt_IsOn (go, "--ewt-hi") ? (float) esl_opt_GetReal (go, "--ewt-hi") :  eslINFINITY;
          if (iinside_wt_sc < ewt_lo || iinside_wt_sc > ewt_hi) {
            n_rejected++;
            free (dsq); dsq = NULL;
            if (tr != NULL) { FreeParsetree (tr); tr = NULL; }
            i--;
            continue;
          }
        }
        weight = (float) pow (2.0, -iinside_wt_sc);
      }

      /* Query cell: get beginsc[v*] + Inside(x[il..ir], v*) directly from the DP.
       * qc_v/qc_j/qc_d are set above in the do_isubtr block (or -1 if not applicable). */
      if (esl_opt_GetBoolean (go, "--ihbanded") && (cm->search_opts & CM_SEARCH_INSIDE)) {
        /* HMM-banded Inside scan */
        float hb_mxsize = esl_opt_GetReal (go, "--mxsize");
        float hb_Mb;
        int64_t hb_ncells;
        double save_tau = cm->tau;
        cm->tau = esl_opt_GetReal (go, "--tau");

        if (esl_opt_GetBoolean (go, "-v")) {
          int dbg_k; long dbg_sum = 0;
          for (dbg_k = 0; dbg_k <= cm->cp9->M; dbg_k++) {
            dbg_sum += cm->cp9->msc[0][dbg_k] + cm->cp9->isc[0][dbg_k];
            dbg_sum += cm->cp9->tsc[CTMM][dbg_k] + cm->cp9->tsc[CTMI][dbg_k];
            dbg_sum += cm->cp9->bsc[dbg_k] + cm->cp9->esc[dbg_k];
          }
          printf ("  PRE-BAND cm: tau=%g cp9_chksum=%ld el_selfsc=%d flags=0x%x\n",
                  cm->tau, dbg_sum, cm->cp9->el_selfsc, cm->cp9->flags);
        }

        if (esl_opt_GetBoolean (go, "-v")) {
          /* Checksum of CP9 matrix state before banding */
          long mx_sum = 0;
          if (cm->cp9_mx->ncells_valid > 0) {
            int dbg_c;
            for (dbg_c = 0; dbg_c < ESL_MIN(100, (int)cm->cp9_mx->ncells_valid); dbg_c++)
              mx_sum += cm->cp9_mx->mmx_mem[dbg_c] + cm->cp9_mx->imx_mem[dbg_c];
          }
          long bmx_sum = 0;
          if (cm->cp9_bmx->ncells_valid > 0) {
            int dbg_c2;
            for (dbg_c2 = 0; dbg_c2 < ESL_MIN(100, (int)cm->cp9_bmx->ncells_valid); dbg_c2++)
              bmx_sum += cm->cp9_bmx->mmx_mem[dbg_c2] + cm->cp9_bmx->imx_mem[dbg_c2];
          }
          /* Checksum of CM begin/end/beginsc/endsc and cp9 begin/end/bsc/esc */
          double be_sum = 0.; long cp9be_sum = 0;
          { int dbg_k2;
            for (dbg_k2 = 0; dbg_k2 < cm->M; dbg_k2++) {
              be_sum += cm->begin[dbg_k2] + cm->end[dbg_k2];
              be_sum += cm->beginsc[dbg_k2] + cm->endsc[dbg_k2];
            }
            for (dbg_k2 = 0; dbg_k2 <= cm->cp9->M; dbg_k2++) {
              cp9be_sum += cm->cp9->bsc[dbg_k2] + cm->cp9->esc[dbg_k2];
            }
          }
          printf ("  PRE-SEQ2BANDS: L=%d fmx=%ld bmx=%ld cm_be=%.6f cp9_be=%ld cp9_flags=0x%x\n",
                  L, mx_sum, bmx_sum, be_sum, cp9be_sum, cm->cp9->flags);
        }

        if ((status = cp9_Seq2Bands (cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                      dsq, 1, L, cm->cp9b, TRUE, PLI_PASS_STD_ANY, 0))
            != eslOK)
          cm_Fail (errbuf);

        /* TODO: Enforce QDB bounds on HMM bands. Currently disabled because
         * post-hoc clamping of hdmin/hdmax breaks the band memory layout
         * used by cm_hb_mx_GrowTo. Needs deeper integration with the band
         * allocation infrastructure (cp9_GrowHDBands). For now, use --noqdb
         * with --ihbanded to avoid QDB/HMM band conflicts. */
#if 0
        if (use_qdbs) {
          if ((status = cp9_EnforceQDBBands (cm, cm->cp9b, cm->smx,
                                             SMX_QDB2_LOOSE, L, errbuf)) != eslOK)
            cm_Fail (errbuf);
        }
#endif

        /* Report banding statistics */
        if (esl_opt_GetBoolean (go, "-v")) {
          int64_t unbanded_cells = (int64_t) cm->M * L * ESL_MIN(L, cm->W);
          cm_hb_mx_SizeNeeded (cm, errbuf, cm->cp9b, L, &hb_ncells, NULL);
          printf ("  HB search_cm: cells=%lld/%lld (%.1f%%)\n",
                  (long long)hb_ncells, (long long)unbanded_cells,
                  100.0 * hb_ncells / (double) unbanded_cells);
          printf ("  === search_cm bands ===\n");
          debug_print_ij_bands (cm);
        }

        if ((status = FastFInsideScanHB (cm, errbuf, cm->hb_mx, hb_mxsize,
                                          dsq, 1, L, cutoff, th,
                                          cm->search_opts & CM_SEARCH_NULL3, 0.,
                                          NULL, NULL, NULL))
            != eslOK)
          cm_Fail (errbuf);
        cm->tau = save_tau;
      } else if (cm->search_opts & CM_SEARCH_INSIDE) {
        if (esl_opt_GetBoolean (go, "--ifloat")) {
          /* Float unbanded Inside (same numerical method as banded) */
          if ((status = FastFInsideScan (cm, errbuf, cm->smx, use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB,
                                         dsq, 1, L, cutoff, th, cm->search_opts & CM_SEARCH_NULL3, 0.,
                                         NULL, NULL, NULL, NULL))
              != eslOK)
            cm_Fail (errbuf);
        } else {
          if ((status = FastIInsideScan (cm, errbuf, cm->smx, use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB,
                                         dsq, 1, L, cutoff, th, cm->search_opts & CM_SEARCH_NULL3, 0.,
                                         NULL, NULL, NULL, NULL, NULL,
                                         (do_isubtr && do_sample) ? isubtr_best_v         : -1,
                                         (do_isubtr && do_sample) ? (int64_t) isubtr_ir   : -1,
                                         (do_isubtr && do_sample) ? isubtr_ir - isubtr_il + 1 : -1,
                                         (do_isubtr && do_sample) ? &isubtr_qc_sc         : NULL,
                                         NULL))
              != eslOK)
            cm_Fail (errbuf);
        }
      } else {
        if ((status = FastCYKScan (cm, errbuf, cm->smx, use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB, dsq, 1,
                                   L, cutoff, th, cm->search_opts & CM_SEARCH_NULL3, 0., NULL, NULL,
                                   NULL, NULL))
            != eslOK)
          cm_Fail (errbuf);
      }
      /* overlaps already removed inside FastCYKScan/FastIInsideScan/FastFInsideScanHB */

      /* --iinside-wt: IS weight from emit_cm Inside score.
       * If --ewt-lo/--ewt-hi are set, the emit_cm scan + rejection was already
       * done above (before the search_cm scan). Otherwise do it here.
       * Also reject sequences whose search CM best hit is outside [ilo, ihi]. */
      if (esl_opt_GetBoolean (go, "--iinside-wt") && do_sample && th->N > 0) {
        /* Rejection based on search CM's best hit */
        float best_inside = th->unsrt[0].score;
        for (h = 1; h < (int) th->N; h++)
          if (th->unsrt[h].score > best_inside)
            best_inside = th->unsrt[h].score;

        if (do_filter && (best_inside < ilo || best_inside > ihi)) {
          n_rejected++;
          cm_tophits_Destroy (th);
          free (dsq); dsq = NULL;
          if (tr != NULL) { FreeParsetree (tr); tr = NULL; }
          i--;
          continue;
        }

        /* Compute IS weight from emit_cm (proposal) Inside score.
         * If iinside_wt_sc was already computed (emit_cm pre-scan above),
         * skip the emit_cm scan here — just use the precomputed weight. */
        if (! esl_opt_GetBoolean (go, "--no-weight") && iinside_wt_sc == IMPOSSIBLE) {
          CM_t *wt_cm = (emit_cm != NULL) ? emit_cm : cm;
          CM_TOPHITS *th_wt = cm_tophits_Create ();
          float wt_best = IMPOSSIBLE;
          if (th_wt == NULL) ESL_FAIL (eslEMEM, errbuf, "out of memory");

          if (esl_opt_GetBoolean (go, "--ihbanded")) {
            /* HMM-banded Inside on emit_cm for IS weight */
            float hb_mxsize = esl_opt_GetReal (go, "--mxsize");
            float hb_Mb;
            double save_tau = wt_cm->tau;
            wt_cm->tau = esl_opt_GetReal (go, "--tau");
            if ((status = cp9_Seq2Bands (wt_cm, errbuf, wt_cm->cp9_mx, wt_cm->cp9_bmx, wt_cm->cp9_bmx,
                                          dsq, 1, L, wt_cm->cp9b, TRUE, PLI_PASS_STD_ANY, 0))
                != eslOK)
              cm_Fail (errbuf);
            /* TODO: Enforce QDB bounds on emit_cm HMM bands (disabled, see #if 0 above) */
#if 0
            if (use_qdbs) {
              if ((status = cp9_EnforceQDBBands (wt_cm, wt_cm->cp9b, wt_cm->smx,
                                                  SMX_QDB2_LOOSE, L, errbuf)) != eslOK)
                cm_Fail (errbuf);
            }
#endif
            if (esl_opt_GetBoolean (go, "-v")) {
              int64_t hb_ncells_wt, hb_ncells_cm;
              cm_hb_mx_SizeNeeded (wt_cm, errbuf, wt_cm->cp9b, L, &hb_ncells_wt, NULL);
              cm_hb_mx_SizeNeeded (cm, errbuf, cm->cp9b, L, &hb_ncells_cm, NULL);
              printf ("  HB emit_cm: cells=%lld  search_cm: cells=%lld  (cp9b: cm=%p wt=%p  cp9_mx: cm=%p wt=%p)\n",
                      (long long)hb_ncells_wt, (long long)hb_ncells_cm,
                      (void*)cm->cp9b, (void*)wt_cm->cp9b,
                      (void*)cm->cp9_mx, (void*)wt_cm->cp9_mx);
            }
            if ((status = FastFInsideScanHB (wt_cm, errbuf, wt_cm->hb_mx, hb_mxsize,
                                              dsq, 1, L, cutoff, th_wt,
                                              wt_cm->search_opts & CM_SEARCH_NULL3, 0.,
                                              NULL, NULL, NULL))
                != eslOK)
              cm_Fail (errbuf);
            wt_cm->tau = save_tau;
          } else if (esl_opt_GetBoolean (go, "--ifloat")) {
            if ((status = FastFInsideScan (wt_cm, errbuf, wt_cm->smx,
                                           use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB,
                                           dsq, 1, L, cutoff, th_wt,
                                           wt_cm->search_opts & CM_SEARCH_NULL3, 0.,
                                           NULL, NULL, NULL, NULL))
                != eslOK)
              cm_Fail (errbuf);
          } else {
            /* In glocal mode, use query cell (v=0, j=L, d=L) to get the
             * full-sequence Inside score directly from the DP, bypassing
             * the greedy non-overlapping hit resolution which may suppress
             * the d=L hit in favor of shorter overlapping hits. */
            float wt_qc_sc = IMPOSSIBLE;
            if ((status = FastIInsideScan (wt_cm, errbuf, wt_cm->smx,
                                           use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB,
                                           dsq, 1, L, cutoff, th_wt,
                                           wt_cm->search_opts & CM_SEARCH_NULL3, 0.,
                                           NULL, NULL, NULL, NULL, NULL,
                                           esl_opt_GetBoolean (go, "--glocal") ? 0         : -1,
                                           esl_opt_GetBoolean (go, "--glocal") ? (int64_t)L : -1,
                                           esl_opt_GetBoolean (go, "--glocal") ? L          : -1,
                                           esl_opt_GetBoolean (go, "--glocal") ? &wt_qc_sc  : NULL,
                                           NULL))
                != eslOK)
              cm_Fail (errbuf);

            if (esl_opt_GetBoolean (go, "--glocal") && wt_qc_sc != IMPOSSIBLE) {
              wt_best = wt_qc_sc;
            } else {
              for (h = 0; h < (int) th_wt->N; h++)
                if (th_wt->unsrt[h].score > wt_best)
                  wt_best = th_wt->unsrt[h].score;
            }
          }
          cm_tophits_Destroy (th_wt);

          if (wt_best != IMPOSSIBLE)
            weight = (float) pow (2.0, -wt_best);
          if (esl_opt_GetBoolean (go, "-v"))
            printf ("  INSIDE wt: emit_cm_inside=%.3f  search_cm_inside=%.3f  parsetree_sc=%.3f  search_cm_ptree=%.3f  cp9_trace=%.3f  cp9_fwd=%.3f  weight=%.6g\n",
                    wt_best, best_inside, parsetree_sc, search_cm_parsetree_sc, cp9_trace_sc, cp9_fwd_sc, weight);
        } else if (iinside_wt_sc != IMPOSSIBLE) {
          /* Weight was already computed in emit_cm pre-scan above */
          if (esl_opt_GetBoolean (go, "-v"))
            printf ("  INSIDE wt (pre): emit_cm_inside=%.3f  search_cm_inside=%.3f  parsetree_sc=%.3f  search_cm_ptree=%.3f  cp9_trace=%.3f  cp9_fwd=%.3f  weight=%.6g\n",
                    iinside_wt_sc, best_inside, parsetree_sc, search_cm_parsetree_sc, cp9_trace_sc, cp9_fwd_sc, weight);
        }
      }

      /* do_isubtr: update IS weight from query-cell Inside score if available.
       * The initial weight (2^(-candidate_sc), set above) is the fallback.
       * isubtr_qc_sc = beginsc[v*] + Inside(x[il..ir], v*) from the DP matrix.
       * This marginalizes over all parsetrees generating x[il..ir] from v*,
       * which is theoretically the correct IS weight for w(x) = P_null(x)/q(x).
       * For short sub-trees it agrees with candidate_sc to within rounding. */
      if (do_isubtr && do_sample && isubtr_best_v != -1) {
        if (isubtr_qc_sc != IMPOSSIBLE) {
          if (! esl_opt_GetBoolean (go, "--no-weight"))
            weight = (float) pow (2.0, -isubtr_qc_sc);
          if (esl_opt_GetBoolean (go, "-v"))
            printf ("  INSIDE qc: v*=%d [%d..%d] inside_sc=%.3f  candidate_sc=%.3f  weight=%.6g\n",
                    isubtr_best_v, isubtr_il, isubtr_ir,
                    isubtr_qc_sc, isubtr_best_candidate, weight);
          n_isubtr_hit_found++;
        } else {
          n_isubtr_hit_missing++;
        }
      }

      /* accumulate dbsize: actual nt searched (unweighted for both IS and random).
       * The hits/Mb criterion counts data points, not IS-equivalent null sequence. */
      dbsize += (double) L;

      if (th->N > 0) {
        if (esl_opt_GetBoolean (go, "--ibest1")) {
          /* --ibest1: keep only the best hit from this sequence */
          float best_sc = th->unsrt[0].score;
          for (h = 1; h < (int) th->N; h++)
            if (th->unsrt[h].score > best_sc)
              best_sc = th->unsrt[h].score;
          if (scN == 0) {
            ESL_ALLOC  (scA, sizeof (float) * (scN + 1));
            if (do_sample) ESL_ALLOC  (wtA, sizeof (float) * (scN + 1));
          } else {
            ESL_RALLOC (scA, tmp, sizeof (float) * (scN + 1));
            if (do_sample) ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + 1));
          }
          scA[scN] = best_sc;
          if (do_sample) wtA[scN] = weight;
          scN++;
        } else {
          /* collect all hits */
          if (scN == 0) {
            ESL_ALLOC  (scA, sizeof (float) * (scN + th->N));
            if (do_sample) ESL_ALLOC  (wtA, sizeof (float) * (scN + th->N));
          } else {
            ESL_RALLOC (scA, tmp, sizeof (float) * (scN + th->N));
            if (do_sample) ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + th->N));
          }
          for (h = 0; h < (int) th->N; h++) {
            scA[(scN + h)] = th->unsrt[h].score;
            if (do_sample) wtA[(scN + h)] = weight;
          }
          scN += th->N;
        }
      }

      cm_tophits_Destroy (th);
    }

    free (dsq);
    if (tr != NULL) {
      FreeParsetree (tr);
      tr = NULL;
    }
  }
  /* free HMM if nec */
  if (esl_opt_GetBoolean (go, "--rhmm")) {
    for (i = 0; i < ghmm_nstates; i++) {
      free (ghmm_eAA[i]);
      free (ghmm_tAA[i]);
    }
    free (ghmm_eAA);
    free (ghmm_tAA);
    free (ghmm_sA);
  }

  /* Report rejection sampling statistics */
  if (do_sample && (do_filter || do_isubtr)) {
    printf ("Rejection sampling: %d emitted, %d accepted, %d rejected (%.4f%% acceptance rate)\n",
            n_emitted, n_emitted - n_rejected, n_rejected,
            (n_emitted > 0) ? 100.0 * (n_emitted - n_rejected) / (double)n_emitted : 0.);
  }
  /* Report IS weight lookup statistics (do_isubtr only) */
  if (do_isubtr && do_sample) {
    int total = n_isubtr_hit_found + n_isubtr_hit_missing;
    printf ("IS weight lookup: %d/%d (%.1f%%) used Inside score; %d (%.1f%%) fell back to parsetree\n",
            n_isubtr_hit_found,  total,
            (total > 0) ? 100.0 * n_isubtr_hit_found  / (double) total : 0.,
            n_isubtr_hit_missing,
            (total > 0) ? 100.0 * n_isubtr_hit_missing / (double) total : 0.);
  }

  if (smx_flank != NULL) cm_scan_mx_Destroy(cm, smx_flank);


  *ret_scN = scN;
  *ret_scA = scA;
  if (ret_wtA != NULL)
    *ret_wtA = wtA;
  if (ret_dbsize != NULL)
    *ret_dbsize = dbsize;

  return eslOK;

ERROR:
  cm_Fail ("Out of memory.");
  return eslEMEM;
}

/* ScoreWeight_t and compare_sw_asc: helper for sorting (score, weight) pairs
 * together by score, used in fit_histogram() to keep weights aligned with scores.
 */
typedef struct {
  double sc;
  double wt;
} ScoreWeight_t;
static int
compare_sw_asc (const void *a, const void *b) {
  const ScoreWeight_t *sa = (const ScoreWeight_t *)a;
  const ScoreWeight_t *sb = (const ScoreWeight_t *)b;
  if (sa->sc < sb->sc)
    return -1;
  if (sa->sc > sb->sc)
    return 1;
  return 0;
}

/* fit_histogram()
 * Create, fill and fit the tail of a histogram to an exponential tail. Data to fill the histogram
 * is given as <scores>. If do_impt is TRUE, <weights> provides importance sampling weights.
 */
static int
fit_histogram (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float tailp,
               double dbsize_nt, int do_impt, float *scores, float *weights, int nscores,
               int exp_mode, double *ret_mu, double *ret_lambda, double *ret_nrandhits,
               float *ret_tailp) {
  int status;
  double mu;
  double lambda;
  int i;
  double *xv;        /* raw data from histogram */
  double *wv = NULL; /* weights corresponding to xv (if do_impt) */
  int n, z;
  double params[2];
  double nrandhits;
  double scaled_nhits_total;
  double scaled_nhits_tail;

  ESL_HISTOGRAM *h = NULL; /* histogram of scores */

  /* Initialize histogram; these numbers are guesses */
  if ((h = esl_histogram_CreateFull (-100., 100., .1)) == NULL)
    return eslEMEM;

  /* fill histogram */
  for (i = 0; i < nscores; i++) {
    if ((status = esl_histogram_Add (h, scores[i])) != eslOK)
      ESL_FAIL (status, errbuf,
                "fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);
    /* printf("%4d %.3f\n", i, scores[i]); */
  }

  /* fit scores to an exponential tail */
#if 0
  if(cfg->rtfitfp != NULL) { 
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
  /* end of if cfg->rtfitfp != NULL) */
#endif

  /* Determine tailp: use hits/Mb criterion (matching cmcalibrate) unless overridden.
   * For IS (do_impt): use WEIGHTED sum to determine where to set the tail threshold,
   * so that the weighted tail contains ~nhits_to_fit effective hits. */
  if (tailp <= 0. && dbsize_nt > 0.) {
    int tailn = ExpModeIsLocal (exp_mode) ? esl_opt_GetInteger (go, "--ltailn")
                                          : esl_opt_GetInteger (go, "--gtailn");
    float nhits_to_fit = (float)tailn * (dbsize_nt / 1e6);

    if (do_impt && weights != NULL) {
      /* Sort scores descending, walk from top accumulating weights until
       * weighted sum >= nhits_to_fit. Set tailp = raw_count / h->n. */
      ScoreWeight_t *sw_tmp = NULL;
      ESL_ALLOC (sw_tmp, sizeof (ScoreWeight_t) * nscores);
      for (i = 0; i < nscores; i++) {
        sw_tmp[i].sc = (double)scores[i];
        sw_tmp[i].wt = (double)weights[i];
      }
      qsort (sw_tmp, nscores, sizeof (ScoreWeight_t), compare_sw_asc);
      /* walk from high scores (end of sorted array) down */
      double wt_cum = 0.;
      int    raw_tail_n = 0;
      for (i = nscores - 1; i >= 0; i--) {
        wt_cum += sw_tmp[i].wt;
        raw_tail_n++;
        if (wt_cum >= nhits_to_fit) break;
      }
      free (sw_tmp);
      tailp = (float) raw_tail_n / (float) h->n;
      if (tailp > 1.) tailp = 1.;
    } else {
      tailp = nhits_to_fit / (float)h->n;
    }
    if (tailp > 1.)
      ESL_FAIL (
          eslERANGE, errbuf,
          "fit_histogram(): only %.1f hits/Mb but need %d; increase iN or lower --ltailn/--gtailn.",
          (float)h->n / (dbsize_nt / 1e6), tailn);
  }

  esl_histogram_GetTailByMass (h, tailp, &xv, &n, &z); /* fit to right 'tailp' fraction */
  if (n <= 1)
    ESL_FAIL (eslERANGE, errbuf,
              "fit_histogram(), too few points in right tailfit: %f fraction of histogram.", tailp);

  if (do_impt) {
    /* Sort (score, weight) pairs together by score so tail weights align with
     * tail scores. weights[] is in insertion order; xv[] from GetTailByMass is
     * in sorted order — they don't correspond by index. */
    ScoreWeight_t *sw = NULL;
    double *xv_sorted = NULL;
    ESL_ALLOC (sw, sizeof (ScoreWeight_t) * nscores);
    for (i = 0; i < nscores; i++) {
      sw[i].sc = (double)scores[i];
      sw[i].wt = (double)weights[i];
    }
    qsort (sw, nscores, sizeof (ScoreWeight_t), compare_sw_asc);

    /* tail starts at index z in the sorted array (same z from GetTailByMass) */
    ESL_ALLOC (xv_sorted, sizeof (double) * n);
    ESL_ALLOC (wv, sizeof (double) * n);
    for (i = 0; i < n; i++) {
      xv_sorted[i] = sw[z + i].sc;
      wv[i] = sw[z + i].wt;
    }

    /* total weight across all scores */
    scaled_nhits_total = 0.;
    for (i = 0; i < nscores; i++)
      scaled_nhits_total += sw[i].wt;

    free (sw);

    /* fit exponential tail using correctly-matched (score, weight) pairs */
    impt_exp_FitComplete (xv_sorted, wv, n, &(params[0]), &(params[1]), &scaled_nhits_tail);
    free (xv_sorted);
  } else {
    esl_exp_FitComplete (xv, n, &(params[0]), &(params[1]));
  }
  esl_histogram_SetExpectedTail (h, params[0], tailp, &esl_exp_generic_cdf, &params);

  mu = params[0];
  lambda = params[1];
  if (isnan (lambda))
    ESL_FAIL (
        eslERANGE, errbuf,
        "fit_histogram(), exp tail fit lambda is NaN, too few hits in histogram. Increase --rL");
  if (isinf (lambda))
    ESL_FAIL (
        eslERANGE, errbuf,
        "fit_histogram(), exp tail fit lambda is inf, too few hits in histogram. Increase --rL");

  if (do_impt) {
    nrandhits = scaled_nhits_tail;
  } else {
    nrandhits = h->n; /* total number of hits in the histogram */
  }
  /* print to output files if nec */
  // if(cfg->exphfp != NULL)
  // esl_histogram_Plot(cfg->exphfp, h);
  // if(cfg->expqfp != NULL) {
  // esl_histogram_PlotQQ(cfg->expqfp, h, &esl_exp_generic_invcdf, params);
  // }

  // if (cfg->expsfp != NULL) {
  // esl_histogram_PlotSurvival(cfg->expsfp, h);
  // esl_exp_Plot(cfg->expsfp, (params[0] - log(1./tailp) / params[1]), 0.693147, esl_exp_surv,
  // h->xmin - 5., h->xmax + 5., 0.1); /* extrapolate mu */
  // }

  esl_histogram_Destroy (h);
  if (wv != NULL)
    free (wv);

  *ret_mu = mu;
  *ret_lambda = lambda;
  *ret_nrandhits = nrandhits;
  if (ret_tailp != NULL)
    *ret_tailp = tailp;
  return eslOK;

ERROR:
  ESL_FAIL (eslEMEM, errbuf, "fit_histogram(): memory allocation error.");
  return eslEMEM;
}

/* Function: sample_sequence_from_cm()
 * Date:     EPN, Mon May  2 08:44:53 2011
 *
 * Purpose:  Generate a dsq from a CM and return it, along with the parsetree score.
 *
 * Returns:  eslOK on success, ESL_DSQ is filled with newly alloc'ed dsq; some other status code on
 * an error,
 */
int
sample_sequence_from_cm (struct cfg_s *cfg, char *errbuf, CM_t *emit_cm, CM_t *score_cm, int *ret_L,
                         ESL_DSQ **ret_dsq, Parsetree_t **ret_tr, float *ret_parsetree_sc) {
  int status;
  int L;
  ESL_SQ *sq;
  ESL_DSQ *dsq;
  Parsetree_t *tr = NULL;
  float parsetree_sc;

  /* Emit from emit_cm (the proposal distribution) */
  if ((status = EmitParsetree (emit_cm, errbuf, cfg->r, "irrelevant", TRUE, &tr, &sq, &L)) != eslOK)
    return status;
  while (L == 0) {
    esl_sq_Destroy (sq);
    if ((status = EmitParsetree (emit_cm, errbuf, cfg->r, "irrelevant", TRUE, &tr, &sq, &L))
        != eslOK)
      return status;
  }

  ESL_ALLOC (dsq, sizeof (ESL_DSQ) * (sq->n + 2));
  memcpy (dsq, sq->dsq, sizeof (ESL_DSQ) * (sq->n + 2));

  /* Calculate the parsetree score under score_cm (the proposal distribution)
   * for the importance sampling weight w = 1/2^parsetree_sc.
   * score_cm should be the same model used for emission (emit_cm),
   * so the weight reflects P(seq|null)/P(seq|proposal). */
  if ((status = ParsetreeScore (score_cm, NULL, errbuf, tr, dsq, FALSE, &parsetree_sc, NULL, NULL,
                                NULL, NULL))
      != eslOK) {
    esl_sq_Destroy (sq);
    return status;
  }

  esl_sq_Destroy (sq);

  *ret_L = L;
  *ret_dsq = dsq;
  *ret_tr = tr;
  if (ret_parsetree_sc != NULL)
    *ret_parsetree_sc = parsetree_sc;

  return eslOK;

ERROR:
  cm_Fail ("Out of memory.");
  return eslEMEM;
}

/* Function:  impt_exp_FitComplete()
 * Incept:    SRE, Wed Aug 10 10:53:47 2005 [St. Louis]
 *            Modified EPN for importance sampling weights
 *
 * Purpose:   Given an array of <n> samples <x[0]..x[n-1]> and
 *            corresponding importance weights <w[0]..w[n-1]>, fit
 *            them to an exponential distribution using weighted MLE.
 *            Return maximum likelihood parameters <ret_mu> and <ret_lambda>.
 *
 * Args:      x          - complete exponentially-distributed data [0..n-1]
 *            w          - importance sampling weights [0..n-1]
 *            n          - number of samples in <x>
 *            ret_mu     - RETURN: lower bound of the distribution (all x_i >= mu)
 *            ret_lambda - RETURN: maximum likelihood estimate of lambda
 *            ret_scaled_nhits - RETURN: sum of weights
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      STL9/138.
 */
int
impt_exp_FitComplete (double *x, double *w, int n, double *ret_mu, double *ret_lambda,
                      double *ret_scaled_nhits) {
  double mu, mean;
  int i;

  double weighted_sum = 0;
  double weight_total = 0;
  double diff = 0;

  /* ML mu is the lowest score. mu=x is ok in the exponential.
   */
  mu = x[0];
  for (i = 1; i < n; i++)
    if (x[i] < mu)
      mu = x[i];

  mean = 0.;
  for (i = 0; i < n; i++) {
    diff = x[i] - mu;
    weighted_sum += diff * w[i];
    weight_total += w[i];

    /*printf("\t\ti: %4d  x[i]: %12.10f  diff: %12.10f  w[i]: %12.10f  prod: %12.10f weight_total:
      %12.10f\n", i, x[i], diff, w[i], diff*w[i], weight_total);*/
  }
  mean = weighted_sum / weight_total;

  // printf("impt n:            %d\n", n);
  // printf("impt weight_total: %.3f\n", weight_total);
  // printf("impt mu:           %.3f\n", mu);
  // printf("impt lambda:       %.3f\n", 1./mean);

  *ret_mu = mu;
  *ret_lambda = 1. / mean;          /* ML estimate trivial & analytic */
  *ret_scaled_nhits = weight_total; /* total weight */
  return eslOK;
}

/* Function: cm_ExpectedParsetreeScore()
 *
 * Purpose:  Compute the expected parsetree score (in bits) for a CM
 *           whose emission probabilities are mixed with the null model
 *           at level <alpha>:
 *              q_k = (1-alpha)*e[v][k] + alpha*null[k]
 *           The expected score is:
 *              E = sum_v psi[v] * sum_k q_k * log2(q_k / null_k)
 *           where psi[v] is the expected occupancy of state v
 *           (from cm_ExpectedStateOccupancy).
 *
 *           alpha=0 gives the expected score for the original CM.
 *           alpha=1 gives 0 (all emissions match null).
 *
 * Args:     cm    - the covariance model (must be in probability form)
 *           alpha - mixing parameter, 0..1
 *
 * Returns:  expected parsetree score in bits
 */
double
cm_ExpectedParsetreeScore (CM_t *cm, double alpha) {
  double *psi = NULL;
  double E = 0.;
  double q; /* mixed emission probability */
  int v, k, l, c;
  int K = cm->abc->K;

  psi = cm_ExpectedStateOccupancy (cm);

  for (v = 0; v < cm->M; v++) {
    if (psi[v] == 0.)
      continue;

    /* Emission scores (alpha-dependent) */
    if (cm->sttype[v] == MP_st) {
      /* pair emitter: q_{k,l} = (1-alpha)*e[v][k*K+l] + alpha*null[k]*null[l] */
      for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
          q = (1.0 - alpha) * cm->e[v][k * K + l] + alpha * cm->null[k] * cm->null[l];
          if (q > 0.)
            E += psi[v] * q * log2 (q / (cm->null[k] * cm->null[l]));
        }
      }
    } else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st
               || cm->sttype[v] == IR_st) {
      /* singlet emitter: q_k = (1-alpha)*e[v][k] + alpha*null[k] */
      for (k = 0; k < K; k++) {
        q = (1.0 - alpha) * cm->e[v][k] + alpha * cm->null[k];
        if (q > 0.)
          E += psi[v] * q * log2 (q / cm->null[k]);
      }
    }

    /* Transition scores (alpha-independent).
     * tsc[v][c] = log2(t[v][c]); always <= 0.
     * B_st and E_st have no transitions. */
    if (cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
      if (v == 0 && (cm->flags & CMH_LOCAL_BEGIN)) {
        /* Local begin: root uses beginsc[y] instead of tsc.
         * cm->begin[y] is the begin transition probability to state y. */
        int y;
        for (y = 0; y < cm->M; y++) {
          if (cm->begin[y] > 0.)
            E += psi[v] * cm->begin[y] * cm->beginsc[y];
        }
      } else {
        /* Normal transitions to children */
        for (c = 0; c < cm->cnum[v]; c++) {
          if (cm->t[v][c] > 0.)
            E += psi[v] * cm->t[v][c] * cm->tsc[v][c];
        }
      }
      /* Local end contribution */
      if ((cm->flags & CMH_LOCAL_END) && cm->end[v] > 0.)
        E += psi[v] * cm->end[v] * cm->endsc[v];
    }
  }

  free (psi);
  return E;
}

/* Function: cm_MixWithNull()
 *
 * Purpose:  Modify emission probabilities of a CM by mixing with the
 *           null model:
 *              e[v][k] = (1-alpha)*e[v][k] + alpha*null[k]
 *           for singlet emitters, and
 *              e[v][k*K+l] = (1-alpha)*e[v][k*K+l] + alpha*null[k]*null[l]
 *           for pair emitters.
 *
 *           Transitions are left untouched.
 *
 *           The CMH_BITS flag is cleared because log-odds scores are
 *           now invalid and must be recalculated (cm_Configure will
 *           do this).
 *
 * Args:     cm    - the covariance model
 *           alpha - mixing parameter, 0..1
 *
 * Returns:  void
 */
void
cm_MixWithNull (CM_t *cm, double alpha) {
  int v, k, l;
  int K = cm->abc->K;

  for (v = 0; v < cm->M; v++) {
    if (cm->sttype[v] == MP_st) {
      for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
          cm->e[v][k * K + l]
              = (1.0 - alpha) * cm->e[v][k * K + l] + alpha * cm->null[k] * cm->null[l];
        }
      }
    } else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st
               || cm->sttype[v] == IR_st) {
      for (k = 0; k < K; k++) {
        cm->e[v][k] = (1.0 - alpha) * cm->e[v][k] + alpha * cm->null[k];
      }
    }
  }
  cm->flags &= ~CMH_BITS; /* log-odds scores are now invalid */
}

/* Function: cp9_EnforceQDBBands()
 *
 * Purpose:  After HMM bands have been computed by cp9_Seq2Bands(),
 *           clamp the d-bands (hdmin/hdmax) to respect QDB bounds
 *           from the scan matrix (smx->dnAAA/dxAAA).  These are
 *           per-position per-state d bounds, matching what
 *           FastIInsideScan uses.
 *
 * Args:     cm      - the CM
 *           cp9b    - the HMM bands to clamp
 *           smx     - scan matrix with QDB d bounds (dnAAA/dxAAA)
 *           qdbidx  - which QDB set: SMX_QDB1_TIGHT or SMX_QDB2_LOOSE
 *           L       - sequence length
 *           errbuf  - for error messages
 *
 * Returns:  eslOK on success
 */
static int
cp9_EnforceQDBBands (CM_t *cm, CP9Bands_t *cp9b, CM_SCAN_MX *smx, int qdbidx, int L, char *errbuf) {
  int v, jp, j;
  int **dnAA, **dxAA;
  int W;

  if (smx == NULL) return eslOK;
  if (qdbidx == SMX_NOQDB) return eslOK;

  dnAA = smx->dnAAA[qdbidx];
  dxAA = smx->dxAAA[qdbidx];
  W    = smx->W;

  for (v = 0; v < cm->M; v++) {
    if (cp9b->jmin[v] > cp9b->jmax[v]) continue;

    for (jp = 0; jp <= cp9b->jmax[v] - cp9b->jmin[v]; jp++) {
      j = jp + cp9b->jmin[v];

      /* Get position-specific QDB d bounds for state v at position j.
       * dnAA/dxAA are indexed by jp_g = j - i0 + 1 (1-based position in seq).
       * For jp_g >= W, use the W entry (steady state). */
      int jp_g = j;  /* j is already 1-based position */
      int qdb_dn, qdb_dx;
      if (jp_g >= W) { qdb_dn = dnAA[W][v]; qdb_dx = dxAA[W][v]; }
      else           { qdb_dn = dnAA[jp_g][v]; qdb_dx = dxAA[jp_g][v]; }

      /* Clamp HMM band d range to QDB bounds, but only if the
       * intersection is non-empty.  If QDB and HMM bands don't
       * overlap in d, keep the original HMM band (conservative). */
      {
        int new_dn = ESL_MAX(cp9b->hdmin[v][jp], qdb_dn);
        int new_dx = ESL_MIN(cp9b->hdmax[v][jp], qdb_dx);
        if (new_dn <= new_dx) {
          cp9b->hdmin[v][jp] = new_dn;
          cp9b->hdmax[v][jp] = new_dx;
        }
      }
    }
  }

  /* Recalculate hd_needed and reallocate hdmin/hdmax memory to match
   * the clamped bands. cp9_GrowHDBands reallocates the flat arrays
   * (hdmin_mem/hdmax_mem) and resets the 2D pointers. */
  cp9b->hd_needed = 0;
  for (v = 0; v < cm->M; v++) {
    if (cp9b->jmin[v] > cp9b->jmax[v]) continue;
    for (jp = 0; jp <= cp9b->jmax[v] - cp9b->jmin[v]; jp++)
      cp9b->hd_needed += cp9b->hdmax[v][jp] - cp9b->hdmin[v][jp] + 1;
  }

  /* Save clamped hdmin/hdmax values, regrow, then restore */
  {
    int *saved_hdmin = NULL, *saved_hdmax = NULL;
    int total = cp9b->hd_needed;
    int status;
    int idx = 0;

    ESL_ALLOC (saved_hdmin, sizeof(int) * ESL_MAX(total, 1));
    ESL_ALLOC (saved_hdmax, sizeof(int) * ESL_MAX(total, 1));

    for (v = 0; v < cm->M; v++) {
      if (cp9b->jmin[v] > cp9b->jmax[v]) continue;
      for (jp = 0; jp <= cp9b->jmax[v] - cp9b->jmin[v]; jp++) {
        saved_hdmin[idx] = cp9b->hdmin[v][jp];
        saved_hdmax[idx] = cp9b->hdmax[v][jp];
        idx++;
      }
    }

    if ((status = cp9_GrowHDBands (cp9b, errbuf)) != eslOK) { free(saved_hdmin); free(saved_hdmax); return status; }

    idx = 0;
    for (v = 0; v < cm->M; v++) {
      if (cp9b->jmin[v] > cp9b->jmax[v]) continue;
      for (jp = 0; jp <= cp9b->jmax[v] - cp9b->jmin[v]; jp++) {
        cp9b->hdmin[v][jp] = saved_hdmin[idx];
        cp9b->hdmax[v][jp] = saved_hdmax[idx];
        idx++;
      }
    }
    free (saved_hdmin);
    free (saved_hdmax);
  }

  return eslOK;

ERROR:
  return eslEMEM;
}


/* Function: cm_InsideScoreAfterMutation()
 *
 * Purpose:  Given Inside matrix alpha and Outside matrix beta for sequence
 *           dsq of length L, compute the EXACT Inside score that would
 *           result from changing dsq[pos] to new_res.
 *
 *           Uses the identity: the (v,j,d) cells that emit a given
 *           position partition the parsetree space, so:
 *             new_sc = FLogsum over {(v,j,d) emitting pos}:
 *                      alpha[v][j][d] + beta[v][j][d] + delta_esc(v,pos,new_res)
 *           where delta_esc = esc(v, new_residue) - esc(v, old_residue).
 *
 *           This is exact because every parsetree emits position pos
 *           through exactly one state, so the sum covers all parsetrees
 *           exactly once.
 *
 * Args:     cm      - the covariance model (configured, with oesc scores)
 *           ins_mx  - Inside DP matrix (from cm_InsideAlign)
 *           out_mx  - Outside DP matrix (from cm_OutsideAlign)
 *           dsq     - digital sequence [1..L]
 *           L       - sequence length
 *           pos     - position to mutate (1..L)
 *           new_res - new residue (0=A, 1=C, 2=G, 3=U)
 *
 * Returns:  the exact Inside score with the mutation applied
 */
static float
cm_InsideScoreAfterMutation (CM_t *cm, CM_MX *ins_mx, CM_MX *out_mx,
                             ESL_DSQ *dsq, int L, int pos, int new_res)
{
  float ***alpha = ins_mx->dp;
  float ***beta  = out_mx->dp;
  float new_sc = IMPOSSIBLE;
  int v, d, j;
  int Kp = cm->abc->Kp;
  int old_res = dsq[pos];

  for (v = 0; v < cm->M; v++) {
    float *esc_v = cm->oesc[v];
    int sd = StateDelta (cm->sttype[v]);

    /* LEFT-emitting states: ML, IL, or MP emitting left residue.
     * v emits position i = j-d+1 on the left. If i = pos, then j = pos+d-1. */
    if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st || cm->sttype[v] == MP_st) {
      for (d = sd; d <= L; d++) {
        j = pos + d - 1;
        if (j > L) break;
        if (! NOT_IMPOSSIBLE(alpha[v][j][d])) continue;
        if (! NOT_IMPOSSIBLE(beta[v][j][d]))  continue;

        float delta_esc;
        if (cm->sttype[v] == MP_st) {
          /* pair emission: changing left residue at pos, right at j is fixed */
          delta_esc = esc_v[new_res * Kp + dsq[j]] - esc_v[old_res * Kp + dsq[j]];
        } else {
          delta_esc = esc_v[new_res] - esc_v[old_res];
        }
        new_sc = FLogsum (new_sc, alpha[v][j][d] + beta[v][j][d] + delta_esc);
      }
    }

    /* RIGHT-emitting states: MR, IR, or MP emitting right residue.
     * v emits position j on the right. If j = pos. */
    if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st || cm->sttype[v] == MP_st) {
      j = pos;
      for (d = sd; d <= j; d++) {
        if (! NOT_IMPOSSIBLE(alpha[v][j][d])) continue;
        if (! NOT_IMPOSSIBLE(beta[v][j][d]))  continue;
        int i = j - d + 1;

        float delta_esc;
        if (cm->sttype[v] == MP_st) {
          /* pair emission: left at i is fixed, changing right residue at pos.
           * For MP, sd=2, so d >= 2, which means i = pos-d+1 <= pos-1.
           * Left branch had j = pos+d-1 >= pos+1, so no overlap. */
          delta_esc = esc_v[dsq[i] * Kp + new_res] - esc_v[dsq[i] * Kp + old_res];
        } else {
          delta_esc = esc_v[new_res] - esc_v[old_res];
        }
        new_sc = FLogsum (new_sc, alpha[v][j][d] + beta[v][j][d] + delta_esc);
      }
    }
  }

  return new_sc;
}


/* Function: cm_DesignSequence()
 *
 * Purpose:  Design a sequence whose glocal Inside score on CM <cm> is
 *           approximately <target_sc>, using stochastic Outside-guided
 *           mutations with trackable proposal probability q(x).
 *
 *           Algorithm:
 *           1. Emit x_0 from <emit_cm> (or <cm>), pre-filter on parsetree score
 *           2. Run cm_InsideAlign to get current Inside score
 *           3. Run cm_OutsideAlign
 *           4. For each (position, residue), compute exact new Inside score
 *           5. STOCHASTIC selection: sample mutation from softmax over
 *              candidates that improve toward target, with known probability
 *           6. Accumulate log q(x) for IS weight
 *           7. Repeat from step 2 until |Inside - target| < tol
 *
 *           The proposal probability is:
 *             q(x_final) = P_emit(x_0) * prod_i P(mut_i | x_{i-1})
 *           where each P(mut_i) is a softmax over improving mutations.
 *
 * Args:     cm         - search CM (configured for glocal Inside)
 *           emit_cm    - proposal CM for emission (alpha-mixed), or NULL for cm
 *           cfg        - config with RNG etc.
 *           errbuf     - for error messages
 *           target_sc  - target Inside score in bits
 *           tol        - tolerance in bits (stop when |Inside - target| < tol)
 *           max_iter   - maximum mutation iterations
 *           verbose    - if TRUE, print progress
 *           ret_dsq    - RETURN: designed sequence (caller frees)
 *           ret_L      - RETURN: sequence length
 *           ret_sc     - RETURN: final Inside score
 *           ret_niter  - RETURN: number of iterations used
 *           ret_log_q  - RETURN: log2(q(x_final)), the proposal log-probability
 *
 * Returns:  eslOK on success
 */
static int
cm_DesignSequence (CM_t *cm, CM_t *emit_cm, struct cfg_s *cfg, const ESL_GETOPTS *go,
                   char *errbuf, float target_sc, float tol, int max_iter, int verbose,
                   ESL_DSQ **ret_dsq, int *ret_L, float *ret_sc, int *ret_niter,
                   double *ret_log_q)
{
  int status;
  ESL_DSQ *dsq = NULL;
  Parsetree_t *tr = NULL;
  int L;
  float sc_inside, sc_ptree;
  float emit_cm_inside;       /* emit_cm Inside score of x_0 for P_emit(x_0) */
  CM_MX *ins_mx = NULL;
  CM_MX *out_mx = NULL;
  int iter;
  int K = cm->abc->K;
  CM_t *ecm = (emit_cm != NULL) ? emit_cm : cm;
  double log_q = 0.;          /* accumulates log2(q(x)) */
  float beta_softmax = 2.0;   /* softmax temperature: higher = more greedy */
  float max_step = (float) esl_opt_GetReal (go, "--istep"); /* max score change per step */

  /* Step 1: Emit and pre-filter on search_cm parsetree score. */
  int n_emitted = 0;
  float prefilter_hi = target_sc + 20.0;
  while (1) {
    ESL_SQ *sq = NULL;
    if (tr != NULL) { FreeParsetree (tr); tr = NULL; }
    if (dsq != NULL) { free (dsq); dsq = NULL; }

    if ((status = EmitParsetree (ecm, errbuf, cfg->r, "design", TRUE, &tr, &sq, &L)) != eslOK)
      return status;
    while (L == 0) {
      esl_sq_Destroy (sq);
      if ((status = EmitParsetree (ecm, errbuf, cfg->r, "design", TRUE, &tr, &sq, &L)) != eslOK)
        return status;
    }
    ESL_ALLOC (dsq, sizeof (ESL_DSQ) * (sq->n + 2));
    memcpy (dsq, sq->dsq, sizeof (ESL_DSQ) * (sq->n + 2));
    esl_sq_Destroy (sq);
    n_emitted++;

    if ((status = ParsetreeScore (cm, NULL, errbuf, tr, dsq, FALSE, &sc_ptree,
                                  NULL, NULL, NULL, NULL)) != eslOK)
      return status;

    if (sc_ptree <= prefilter_hi) break;
    if (n_emitted > 10000) {
      ESL_FAIL (eslERANGE, errbuf, "cm_DesignSequence: 10000 emissions without finding sc_ptree <= %.1f", prefilter_hi);
    }
  }
  FreeParsetree (tr); tr = NULL;

  /* Compute P_emit(x_0) = emit_cm Inside score of the starting sequence.
   * This is the first factor of q(x): q = P_emit(x_0) * prod(P_mut_i).
   * log_q starts as log2(P_emit(x_0) / P_null(x_0)) = emit_cm_inside(x_0).
   * (We track log-odds ratios; the P_null factors cancel in the final IS weight.) */
  {
    CM_MX *emit_mx = cm_mx_Create (ecm->M);
    if ((status = cm_InsideAlign (ecm, errbuf, dsq, L, 512.0, emit_mx, &emit_cm_inside)) != eslOK)
      cm_Fail (errbuf);
    cm_mx_Destroy (emit_mx);
  }
  log_q = (double) emit_cm_inside;  /* log2(P_emit(x_0) / P_null(x_0)) */

  if (verbose)
    printf ("DESIGN: emitted %d seqs, sc_ptree=%.3f emit_inside=%.3f log_q=%.3f (target=%.3f)\n",
            n_emitted, sc_ptree, emit_cm_inside, log_q, target_sc);

  /* Allocate Inside and Outside matrices */
  ins_mx = cm_mx_Create (cm->M);
  out_mx = cm_mx_Create (cm->M);

  /* Candidate mutations array: at most 3*L candidates (3 alternative residues per position) */
  int max_cand = 3 * (L + 1);
  float *cand_sc;       /* predicted new Inside score for each candidate */
  int   *cand_pos;      /* position of each candidate */
  int   *cand_res;      /* new residue of each candidate */
  ESL_ALLOC (cand_sc,  sizeof (float) * max_cand);
  ESL_ALLOC (cand_pos, sizeof (int)   * max_cand);
  ESL_ALLOC (cand_res, sizeof (int)   * max_cand);

  /* Step 2-7: Iterative stochastic mutation */
  for (iter = 0; iter < max_iter; iter++) {

    /* Step 2: Run Inside */
    if ((status = cm_InsideAlign (cm, errbuf, dsq, L, 512.0, ins_mx, &sc_inside)) != eslOK)
      ESL_FAIL (status, errbuf, "cm_DesignSequence: cm_InsideAlign failed");

    if (verbose)
      printf ("  iter %3d: Inside=%.3f  target=%.3f  diff=%+.3f  log_q=%.3f\n",
              iter, sc_inside, target_sc, sc_inside - target_sc, log_q);

    /* Check convergence */
    if (fabs (sc_inside - target_sc) < tol) {
      if (verbose) printf ("  CONVERGED at iter %d\n", iter);
      break;
    }

    /* Step 3: Run Outside */
    if ((status = cm_OutsideAlign (cm, errbuf, dsq, L, 512.0, FALSE, out_mx, ins_mx, NULL)) != eslOK)
      ESL_FAIL (status, errbuf, "cm_DesignSequence: cm_OutsideAlign failed");

    /* Step 4: Enumerate candidate mutations.
     * Keep candidates that improve toward target AND don't jump more than max_step bits.
     * Small max_step = more intermediate data points with better weights. */
    int n_cand = 0;
    float cur_dist = fabs (sc_inside - target_sc);

    for (int p = 1; p <= L; p++) {
      for (int r = 0; r < K; r++) {
        if (r == dsq[p]) continue;
        float new_sc = cm_InsideScoreAfterMutation (cm, ins_mx, out_mx, dsq, L, p, r);
        float new_dist = fabs (new_sc - target_sc);
        float step_size = fabs (new_sc - sc_inside);
        if (new_dist < cur_dist && step_size <= max_step) {
          cand_sc[n_cand]  = new_sc;
          cand_pos[n_cand] = p;
          cand_res[n_cand] = r;
          n_cand++;
        }
      }
    }

    /* If no candidates within max_step, relax the step limit for this iteration */
    if (n_cand == 0) {
      for (int p = 1; p <= L; p++) {
        for (int r = 0; r < K; r++) {
          if (r == dsq[p]) continue;
          float new_sc = cm_InsideScoreAfterMutation (cm, ins_mx, out_mx, dsq, L, p, r);
          float new_dist = fabs (new_sc - target_sc);
          if (new_dist < cur_dist) {
            cand_sc[n_cand]  = new_sc;
            cand_pos[n_cand] = p;
            cand_res[n_cand] = r;
            n_cand++;
          }
        }
      }
    }

    if (n_cand == 0) {
      if (verbose) printf ("  NO IMPROVING MUTATION found at iter %d\n", iter);
      break;
    }

    /* Step 5: Stochastic selection via softmax.
     * P(candidate i) = exp(-beta * |new_sc_i - target|) / Z
     * Higher beta = more greedy (concentrates on best candidates).
     * We use the negative distance as the "energy" so closer-to-target
     * candidates have higher probability. */
    double *log_prob;
    ESL_ALLOC (log_prob, sizeof (double) * n_cand);
    double max_logp = -eslINFINITY;
    for (int c = 0; c < n_cand; c++) {
      log_prob[c] = -beta_softmax * fabs (cand_sc[c] - target_sc);
      if (log_prob[c] > max_logp) max_logp = log_prob[c];
    }
    /* Convert to probabilities (numerically stable softmax) */
    double Z = 0.;
    for (int c = 0; c < n_cand; c++) {
      log_prob[c] = exp (log_prob[c] - max_logp);
      Z += log_prob[c];
    }
    for (int c = 0; c < n_cand; c++)
      log_prob[c] /= Z;  /* now log_prob[c] is actually prob[c] */

    /* Sample from the distribution */
    double u = esl_random (cfg->r);
    double cum = 0.;
    int chosen = n_cand - 1;  /* default to last */
    for (int c = 0; c < n_cand; c++) {
      cum += log_prob[c];
      if (u <= cum) { chosen = c; break; }
    }

    /* Step 6: Accumulate log q(x).
     * The mutation probability is log_prob[chosen] (which is P(mut | x_{iter})).
     * We add log2(P_mut) to log_q. */
    log_q += log2 (log_prob[chosen]);

    if (verbose)
      printf ("  MUTATE pos=%d %c->%c  predicted_sc=%.3f  P(mut)=%.6f  log_q=%.3f  (n_cand=%d)\n",
              cand_pos[chosen], "ACGU"[dsq[cand_pos[chosen]]], "ACGU"[cand_res[chosen]],
              cand_sc[chosen], log_prob[chosen], log_q, n_cand);

    /* Apply mutation */
    dsq[cand_pos[chosen]] = cand_res[chosen];
    free (log_prob);
  }

  /* Final Inside score verification */
  if ((status = cm_InsideAlign (cm, errbuf, dsq, L, 512.0, ins_mx, &sc_inside)) != eslOK)
    ESL_FAIL (status, errbuf, "cm_DesignSequence: final cm_InsideAlign failed");

  if (verbose)
    printf ("  FINAL: Inside=%.3f  target=%.3f  diff=%+.3f  niter=%d  log_q=%.3f\n",
            sc_inside, target_sc, sc_inside - target_sc, iter, log_q);

  cm_mx_Destroy (ins_mx);
  cm_mx_Destroy (out_mx);
  free (cand_sc);
  free (cand_pos);
  free (cand_res);

  *ret_dsq   = dsq;
  *ret_L     = L;
  *ret_sc    = sc_inside;
  if (ret_niter != NULL) *ret_niter = iter;
  if (ret_log_q != NULL) *ret_log_q = log_q;
  return eslOK;

ERROR:
  if (ins_mx != NULL) cm_mx_Destroy (ins_mx);
  if (out_mx != NULL) cm_mx_Destroy (out_mx);
  if (dsq != NULL) free (dsq);
  ESL_FAIL (eslEMEM, errbuf, "cm_DesignSequence: memory allocation error");
}


/* Function: cm_MaxLocalBeginScore()
 *
 * Purpose:  Given a filled Inside alignment matrix, return
 *           max_v [beginsc[v] + alpha[v][L][L]] for all states v
 *           with valid local begins, restricted to d=L. This is the
 *           best FULL-SEQUENCE local Inside score.
 *           For glocal mode (no local begins), returns alpha[0][L][L].
 */
static float
cm_MaxLocalBeginScore (CM_t *cm, CM_MX *ins_mx, int L)
{
  float ***alpha = ins_mx->dp;
  float best = IMPOSSIBLE;

  if (cm->flags & CMH_LOCAL_BEGIN) {
    for (int v = 1; v < cm->M; v++) {
      if (NOT_IMPOSSIBLE(cm->beginsc[v]) && NOT_IMPOSSIBLE(alpha[v][L][L])) {
        float sc = cm->beginsc[v] + alpha[v][L][L];
        if (sc > best) best = sc;
      }
    }
  }
  /* Also include alpha[0][L][L] (the root, which is always valid) */
  if (NOT_IMPOSSIBLE(alpha[0][L][L]) && alpha[0][L][L] > best)
    best = alpha[0][L][L];

  return best;
}


/* Function: cm_BestLocalHitScore()
 *
 * Purpose:  Given a filled Inside alignment matrix, return the best
 *           hit score over ALL (v, j, d) — i.e., the highest-scoring
 *           local alignment to any subsequence. This matches what
 *           FastIInsideScan reports as the best hit.
 *
 *           Score = beginsc[v] + alpha[v][j][d] for v with local begin
 *                   (or v=0, the root, with implicit beginsc=0)
 *
 *           Cost: O(M * L²) — fast (just scans the matrix).
 */
static float
cm_BestLocalHitScore (CM_t *cm, CM_MX *ins_mx, int L)
{
  float ***alpha = ins_mx->dp;
  float best = IMPOSSIBLE;
  int v, j, d;

  for (v = 0; v < cm->M; v++) {
    float bsc;
    if (v == 0) {
      bsc = 0.;  /* root, always reachable */
    } else if (cm->flags & CMH_LOCAL_BEGIN && NOT_IMPOSSIBLE(cm->beginsc[v])) {
      bsc = cm->beginsc[v];
    } else {
      continue;  /* no local begin to this state */
    }
    for (j = 0; j <= L; j++) {
      for (d = 0; d <= j; d++) {
        if (NOT_IMPOSSIBLE(alpha[v][j][d])) {
          float sc = bsc + alpha[v][j][d];
          if (sc > best) best = sc;
        }
      }
    }
  }
  return best;
}


/* Same as cm_BestLocalHitScore but also returns the (v, j, d) of the max */
static float
cm_BestLocalHitScoreVJD (CM_t *cm, CM_MX *ins_mx, int L, int *ret_v, int *ret_j, int *ret_d)
{
  float ***alpha = ins_mx->dp;
  float best = IMPOSSIBLE;
  int v, j, d;
  int best_v = -1, best_j = -1, best_d = -1;

  for (v = 0; v < cm->M; v++) {
    float bsc;
    if (v == 0) {
      bsc = 0.;
    } else if (cm->flags & CMH_LOCAL_BEGIN && NOT_IMPOSSIBLE(cm->beginsc[v])) {
      bsc = cm->beginsc[v];
    } else {
      continue;
    }
    for (j = 0; j <= L; j++) {
      for (d = 0; d <= j; d++) {
        if (NOT_IMPOSSIBLE(alpha[v][j][d])) {
          float sc = bsc + alpha[v][j][d];
          if (sc > best) { best = sc; best_v = v; best_j = j; best_d = d; }
        }
      }
    }
  }
  if (ret_v) *ret_v = best_v;
  if (ret_j) *ret_j = best_j;
  if (ret_d) *ret_d = best_d;
  return best;
}


/* Function: cm_InsideAlign_partial()
 *
 * Purpose:  Incremental Inside DP update after a single-residue mutation.
 *           Given an alpha matrix that represents the Inside DP for the
 *           sequence BEFORE the mutation at position p, update only the
 *           cells whose subsequence contains position p (the "affected"
 *           cells). Cells whose subsequence does NOT contain p have
 *           unchanged alpha values and are skipped.
 *
 *           Cells are processed in the same topological order as
 *           cm_InsideAlign: states v from M-1 down to 0, then for each
 *           v, j and d in increasing order. Within affected cells,
 *           the standard Inside recurrence is applied.
 *
 *           For unaffected child cells, we read the cached value from
 *           the matrix (which is the correct old value, since the cell
 *           was not affected by the mutation).
 *
 *           ASSUMPTIONS for this implementation (subset of full Inside):
 *             - non-banded glocal mode (no QDB, no HMM bands)
 *             - no local ends (CMH_LOCAL_END not set)
 *             - dsq has been MUTATED (position p has the new residue)
 *             - mx already contains alpha values for the OLD sequence
 *
 * Args:     cm     - the CM (configured for glocal Inside)
 *           errbuf - error buffer
 *           dsq    - mutated digital sequence (position p has new residue)
 *           L      - sequence length
 *           mx     - DP matrix, contains OLD alpha; will be updated in place
 *           p      - mutation position (1..L)
 *           ret_sc - RETURN: new alpha[0][L][L]
 *
 * Returns:  eslOK on success.
 */
static int
cm_InsideAlign_partial (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L,
                        CM_MX *mx, int p, float *ret_sc)
{
  float ***alpha = mx->dp;
  int v, j, d, i, k;
  int yoffset;
  float tsc;
  int Kp = cm->abc->Kp;
  int status;

  /* Macro to check if cell (j, d) is affected by mutation at p:
   * cell's subsequence is x[j-d+1..j], affected iff j-d+1 <= p <= j */
  #define IS_AFFECTED(jj, dd) (((jj) - (dd) + 1) <= p && p <= (jj))

  /* Precompute EL self-loop scores: el_scA[d] = cm->el_selfsc * d
   * This is the score for emitting d residues via the EL (local end) loop. */
  float *el_scA = NULL;
  if (cm->flags & CMH_LOCAL_END) {
    ESL_ALLOC (el_scA, sizeof (float) * (L + 1));
    for (d = 0; d <= L; d++) el_scA[d] = cm->el_selfsc * d;
  }

  /* Process states v from M-1 down to 0 (topological order) */
  for (v = cm->M - 1; v >= 0; v--) {
    float const *esc_v = cm->oesc[v];
    float const *tsc_v = cm->tsc[v];
    int sd  = StateDelta (cm->sttype[v]);
    int sdl = StateLeftDelta (cm->sttype[v]);
    int sdr = StateRightDelta (cm->sttype[v]);

    if (cm->sttype[v] == E_st) {
      /* E_st: alpha[v][j][0] = 0, all other d are IMPOSSIBLE.
       * No emission, no transitions, not affected by any mutation.
       * Skip entirely - cached values are correct. */
      continue;
    }

    if (cm->sttype[v] == B_st) {
      /* B_st: bifurcation. alpha[v][j][d] = FLogsum_k alpha[w][j-k][d-k] + alpha[z][j][k]
       * where w = cfirst[v], z = cnum[v]. */
      int w = cm->cfirst[v];
      int z = cm->cnum[v];
      for (j = 0; j <= L; j++) {
        for (d = 0; d <= j; d++) {
          if (! IS_AFFECTED(j, d)) continue;
          /* Reset to initial value: EL contribution if local end, else IMPOSSIBLE.
           * Note: B_st has sd=0, so INIT_CELL_VAL uses el_scA[d]. */
          alpha[v][j][d] = (cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v]) ?
                           (el_scA[d] + cm->endsc[v]) : IMPOSSIBLE;
          for (k = 0; k <= d; k++) {
            alpha[v][j][d] = FLogsum (alpha[v][j][d],
                                       alpha[w][j-k][d-k] + alpha[z][j][k]);
          }
        }
      }
      continue;
    }

    /* Non-E, non-B state. Includes IL, IR (self-transit) and ML, MR, MP, D, S
     * (no self-transit). For self-transit states, the loop order matters:
     * j outer, d inner (so when computing (j,d), (j,d-1) is already done).
     * For non-self-transit, the order is flexible.
     * We use the same order for all (j outer, d inner) for simplicity. */
    for (j = sdr; j <= L; j++) {
      int j_sdr = j - sdr;
      for (d = sd; d <= j; d++) {
        if (! IS_AFFECTED(j, d)) continue;
        int d_sd = d - sd;
        i = j - d + 1;

        /* Reset to initial value: EL contribution if local end, else IMPOSSIBLE */
        alpha[v][j][d] = (cm->flags & CMH_LOCAL_END) && NOT_IMPOSSIBLE(cm->endsc[v]) ?
                         (el_scA[d - sd] + cm->endsc[v]) : IMPOSSIBLE;

        /* Sum transitions from children */
        for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++) {
          int y = cm->cfirst[v] + yoffset;
          tsc = tsc_v[yoffset];
          alpha[v][j][d] = FLogsum (alpha[v][j][d],
                                     alpha[y][j_sdr][d_sd] + tsc);
        }

        /* Add emission score, if any */
        switch (cm->sttype[v]) {
          case ML_st:
          case IL_st:
            alpha[v][j][d] += esc_v[dsq[i]];
            break;
          case MR_st:
          case IR_st:
            alpha[v][j][d] += esc_v[dsq[j]];
            break;
          case MP_st:
            alpha[v][j][d] += esc_v[dsq[i] * Kp + dsq[j]];
            break;
          default:
            break; /* D, S: no emission */
        }

        /* Clamp to IMPOSSIBLE floor */
        if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
      }
    }
  }

  #undef IS_AFFECTED

  /* Local begin handling: in local mode, alpha[0][L][L] gets contributions
   * from local begins to other states. cm_InsideAlign does:
   *   bsc = FLogsum over v of (alpha[v][L][L] + cm->beginsc[v])
   *   alpha[0][L][L] = FLogsum(alpha[0][L][L], bsc)
   *
   * The mutation might have changed alpha[v][L][L] for v's whose subtree
   * subsequence (always [1..L]) contains p — which is always true since p
   * is in [1..L]. So we need to recompute the root cell's contribution
   * from local begins after the partial update.
   *
   * To do this correctly, we need to know what alpha[0][L][L] WOULD be
   * without local begins (just the standard recurrence). We re-derive it:
   *   stripped_root = standard_alpha[0][L][L]
   *   bsc = FLogsum over v with valid begin of (alpha[v][L][L] + beginsc[v])
   *   new alpha[0][L][L] = FLogsum(stripped_root, bsc)
   *
   * Problem: we don't know stripped_root from the matrix alone (alpha[0][L][L]
   * already has bsc added in from previous full Inside calls, and our partial
   * DP just recomputed it WITHOUT the bsc, since the recurrence for v=0 doesn't
   * include local begin contributions).
   *
   * So our partial DP's alpha[0][L][L] IS the stripped_root (without bsc).
   * We need to add bsc to it. */
  if (cm->flags & CMH_LOCAL_BEGIN) {
    float bsc = IMPOSSIBLE;
    for (v = 1; v < cm->M; v++) {
      if (NOT_IMPOSSIBLE(cm->beginsc[v]) && NOT_IMPOSSIBLE(alpha[v][L][L])) {
        bsc = FLogsum (bsc, alpha[v][L][L] + cm->beginsc[v]);
      }
    }
    if (NOT_IMPOSSIBLE(bsc))
      alpha[0][L][L] = FLogsum (alpha[0][L][L], bsc);
  }

  if (el_scA != NULL) free (el_scA);
  if (ret_sc != NULL) *ret_sc = alpha[0][L][L];
  return eslOK;

ERROR:
  if (el_scA != NULL) free (el_scA);
  ESL_FAIL (eslEMEM, errbuf, "cm_InsideAlign_partial: memory allocation error");
}


/* Function: cm_MCMC_tail()
 *
 * Purpose:  MCMC sampling from the tail of the null distribution.
 *           Generates sequences from the uniform distribution over
 *           {x : Inside(x, CM) >= mu}, using Metropolis-Hastings with
 *           single-residue mutation proposals.
 *
 *           No IS weights needed — samples are from the correct distribution
 *           by construction.
 *
 *           Algorithm per chain:
 *           1. Design a starting sequence with Inside ≈ mu (using cm_DesignSequence)
 *           2. For each step:
 *              a. Pick random position p, random alternative residue r
 *              b. Compute Inside(x') using cm_InsideScoreAfterMutation (O(M))
 *              c. If Inside(x') >= mu: accept, recompute Inside+Outside (O(L²M))
 *              d. If Inside(x') < mu: reject, try another proposal
 *              e. After acceptance, record score (if past burn-in)
 *           3. Repeat for n_steps accepted steps
 *
 * Args:     cm         - search CM (configured for glocal Inside)
 *           emit_cm    - proposal CM for initial sequence design, or NULL
 *           cfg        - config with RNG etc.
 *           go         - command line options
 *           errbuf     - for error messages
 *           mu         - score threshold for tail
 *           n_chains   - number of independent chains
 *           n_steps    - accepted steps per chain (after burn-in)
 *           n_burnin   - accepted steps to discard at start of each chain
 *           verbose    - if TRUE, print progress
 *           ret_scores - RETURN: array of scores [0..ret_N-1] (caller frees)
 *           ret_N      - RETURN: number of scores collected
 *
 * Returns:  eslOK on success
 */
static int
cm_MCMC_tail (CM_t *cm, CM_t *emit_cm, struct cfg_s *cfg, const ESL_GETOPTS *go,
              char *errbuf, float mu, int n_chains, int n_steps, int n_burnin,
              int use_maxv, int verbose, float **ret_scores, int *ret_N)
{
  int status;
  int K = cm->abc->K;
  int total_N = n_chains * n_steps;
  float *all_scores = NULL;
  int sc_idx = 0;

  ESL_ALLOC (all_scores, sizeof (float) * total_N);

  for (int chain = 0; chain < n_chains; chain++) {
    ESL_DSQ *dsq = NULL;
    int L;
    float sc_inside;
    int design_niter;
    double design_log_q;
    CM_MX *ins_mx = NULL;

    /* Step 1: Get a starting sequence with score >= mu.
     * For sum mode: design a sequence targeting alpha[0][L][L] ≈ mu.
     * For max-v mode: cm_DesignSequence targets the wrong quantity, so
     * instead generate random sequences until we find one with best
     * local hit >= mu. */
    L = cm->clen;
    ins_mx = cm_mx_Create (cm->M);

    if (use_maxv) {
      /* Random sampling: find a starting sequence with best hit >= mu */
      int n_tries = 0;
      ESL_ALLOC (dsq, sizeof (ESL_DSQ) * (L + 2));
      while (1) {
        if ((status = esl_rsq_xfIID (cfg->r, cm->null, cm->abc->K, L, dsq)) != eslOK)
          cm_Fail ("ERROR generating random sequence");
        n_tries++;
        if ((status = cm_InsideAlign (cm, errbuf, dsq, L, 512.0, ins_mx, &sc_inside)) != eslOK)
          cm_Fail (errbuf);
        sc_inside = cm_BestLocalHitScore (cm, ins_mx, L);
        if (sc_inside >= mu) break;
        if (n_tries > 100000) {
          ESL_FAIL (eslERANGE, errbuf, "cm_MCMC_tail: 100K random tries without finding seq with best hit >= mu=%.3f", mu);
        }
      }
      if (verbose)
        printf ("MCMC chain %d: start L=%d best_hit=%.3f (random, %d tries)\n",
                chain, L, sc_inside, n_tries);
    } else {
      /* Design a starting sequence with alpha[0][L][L] ≈ mu */
      if ((status = cm_DesignSequence (cm, emit_cm, cfg, go, errbuf,
                                        mu, 0.5, 200, FALSE,
                                        &dsq, &L, &sc_inside, &design_niter,
                                        &design_log_q)) != eslOK)
        ESL_FAIL (status, errbuf, "cm_MCMC_tail: failed to design starting sequence for chain %d", chain);

      if (verbose)
        printf ("MCMC chain %d: start L=%d Inside=%.3f (designed in %d iters)\n",
                chain, L, sc_inside, design_niter);

      /* Compute initial Inside */
      if ((status = cm_InsideAlign (cm, errbuf, dsq, L, 512.0, ins_mx, &sc_inside)) != eslOK)
        cm_Fail (errbuf);
    }

    /* Step 2-3: Simple Metropolis MCMC.
     * Propose random single-residue mutation, accept if score stays >= mu.
     * For sum mode (alpha[0][L][L]): symmetric proposal, no MH correction needed.
     * For max-v mode: also symmetric (uniform proposal over positions × residues).
     * Both use cm_InsideAlign per step (O(L²M)). */
    int accepted = 0;
    int collected = 0;
    int total_proposals = 0;

    while (collected < n_steps) {
      /* Propose: random position, random alternative residue */
      int p = 1 + esl_rnd_Roll (cfg->r, L);
      int r = esl_rnd_Roll (cfg->r, K - 1);
      if (r >= dsq[p]) r++;
      total_proposals++;

      /* Apply mutation tentatively */
      int old_res = dsq[p];
      dsq[p] = r;

      /* Capture state of OLD best cell (from previous accept) BEFORE mutation.
       * old_best_v/j/d are persistent across steps within a chain. */
      static int old_best_v = -1, old_best_j = -1, old_best_d = -1;
      float old_score = sc_inside;  /* score at the old best cell, pre-mutation */

      /* For TopK analysis: collect the top-K cells from the OLD matrix
       * (before mutation). We use these to check if the new best cell came from
       * the top-K of the old matrix. */
      #define TOPK_TRACK 50
      int topk_v[TOPK_TRACK], topk_j[TOPK_TRACK], topk_d[TOPK_TRACK];
      float topk_sc[TOPK_TRACK];
      int topk_n = 0;
      static int do_topk_analysis = -1;
      if (do_topk_analysis == -1)
        do_topk_analysis = use_maxv ? 1 : 0;
      if (do_topk_analysis) {
        /* Collect top-K cells from the current alpha matrix (before mutation) */
        for (int vv = 0; vv < cm->M; vv++) {
          float bsc = (vv == 0) ? 0.0f :
            ((cm->flags & CMH_LOCAL_BEGIN) && NOT_IMPOSSIBLE(cm->beginsc[vv])) ?
            cm->beginsc[vv] : -INFINITY;
          if (bsc < -1e30) continue;
          for (int jj = 0; jj <= L; jj++) {
            for (int dd = 0; dd <= jj; dd++) {
              if (! NOT_IMPOSSIBLE(ins_mx->dp[vv][jj][dd])) continue;
              float sc = bsc + ins_mx->dp[vv][jj][dd];
              /* Insert into top-K if better than current min */
              if (topk_n < TOPK_TRACK) {
                /* Find insertion point */
                int pos = topk_n;
                while (pos > 0 && topk_sc[pos-1] < sc) {
                  topk_sc[pos] = topk_sc[pos-1];
                  topk_v[pos] = topk_v[pos-1];
                  topk_j[pos] = topk_j[pos-1];
                  topk_d[pos] = topk_d[pos-1];
                  pos--;
                }
                topk_sc[pos] = sc;
                topk_v[pos] = vv;
                topk_j[pos] = jj;
                topk_d[pos] = dd;
                topk_n++;
              } else if (sc > topk_sc[TOPK_TRACK-1]) {
                /* Better than worst in topk; insert */
                int pos = TOPK_TRACK - 1;
                while (pos > 0 && topk_sc[pos-1] < sc) {
                  topk_sc[pos] = topk_sc[pos-1];
                  topk_v[pos] = topk_v[pos-1];
                  topk_j[pos] = topk_j[pos-1];
                  topk_d[pos] = topk_d[pos-1];
                  pos--;
                }
                topk_sc[pos] = sc;
                topk_v[pos] = vv;
                topk_j[pos] = jj;
                topk_d[pos] = dd;
              }
            }
          }
        }
      }

      /* Recompute Inside for proposed sequence using partial DP.
       * The matrix ins_mx contains the Inside DP for the unmutated sequence.
       * cm_InsideAlign_partial updates only cells affected by mutation at p. */
      float new_sc;
      if ((status = cm_InsideAlign_partial (cm, errbuf, dsq, L, ins_mx, p, &new_sc)) != eslOK)
        cm_Fail (errbuf);
      int new_v = -1, new_j = -1, new_d = -1;
      if (use_maxv)
        new_sc = cm_BestLocalHitScoreVJD (cm, ins_mx, L, &new_v, &new_j, &new_d);

      if (new_sc >= mu) {
        /* Accept */
        sc_inside = new_sc;
        accepted++;

        /* TRACKING: log enriched info per accept */
        if (use_maxv) {
          static FILE *track_fp = NULL;
          if (track_fp == NULL) {
            track_fp = fopen ("/tmp/mcmc_track.txt", "w");
            if (track_fp) fprintf (track_fp, "# accepted mut_pos old_res new_res old_j old_d new_j new_d old_score new_score score_at_old rank_in_topk topk_vs_true topk_best_rank\n");
          }
          if (track_fp) {
            /* Compute score at the OLD best cell after the mutation */
            float score_at_old = -999.0;
            if (old_best_v >= 0 && old_best_j >= 0 && old_best_d >= 0) {
              float bsc = (old_best_v == 0) ? 0.0 :
                          ((cm->flags & CMH_LOCAL_BEGIN) && NOT_IMPOSSIBLE(cm->beginsc[old_best_v])) ?
                          cm->beginsc[old_best_v] : -999.0;
              if (bsc > -998.0)
                score_at_old = bsc + ins_mx->dp[old_best_v][old_best_j][old_best_d];
            }

            /* Find rank of new best cell in OLD top-K. -1 if not in top-K. */
            int rank_in_topk = -1;
            for (int kk = 0; kk < topk_n; kk++) {
              if (topk_v[kk] == new_v && topk_j[kk] == new_j && topk_d[kk] == new_d) {
                rank_in_topk = kk;
                break;
              }
            }
            /* Also: would tracking top-K alone (with their new scores after mutation)
             * find the new best? For each cell in top-K, compute its NEW score
             * (the value is now in ins_mx after partial DP) and find the max. */
            float topk_best = -INFINITY;
            int topk_best_rank = -1;
            for (int kk = 0; kk < topk_n; kk++) {
              float bsc = (topk_v[kk] == 0) ? 0.0f :
                ((cm->flags & CMH_LOCAL_BEGIN) && NOT_IMPOSSIBLE(cm->beginsc[topk_v[kk]])) ?
                cm->beginsc[topk_v[kk]] : -1e30f;
              float sc = bsc + ins_mx->dp[topk_v[kk]][topk_j[kk]][topk_d[kk]];
              if (sc > topk_best) { topk_best = sc; topk_best_rank = kk; }
            }
            float topk_vs_true = topk_best - new_sc;  /* difference from true best */

            fprintf (track_fp, "%d %d %d %d %d %d %d %d %.4f %.4f %.4f %d %.4f %d\n",
                     accepted, p, old_res, r,
                     old_best_j, old_best_d, new_j, new_d,
                     old_score, new_sc, score_at_old,
                     rank_in_topk, topk_vs_true, topk_best_rank);
            fflush (track_fp);
          }
          /* Update tracked old best for next step */
          old_best_v = new_v;
          old_best_j = new_j;
          old_best_d = new_d;
        }

        /* Collect score if past burn-in */
        if (accepted > n_burnin) {
          all_scores[sc_idx++] = sc_inside;
          collected++;
          if (verbose && (collected % 20 == 0 || collected == n_steps))
            printf ("  chain %d: collected %d/%d  score=%.3f  accept_rate=%.1f%%\n",
                    chain, collected, n_steps, sc_inside,
                    100.0 * accepted / total_proposals);
        }
      } else {
        /* Reject: revert mutation. Use partial DP to undo the change in ins_mx
         * (revert is also a single mutation: change dsq[p] back from r to old_res) */
        dsq[p] = old_res;
        if ((status = cm_InsideAlign_partial (cm, errbuf, dsq, L, ins_mx, p, &sc_inside)) != eslOK)
          cm_Fail (errbuf);
        if (use_maxv)
          sc_inside = cm_BestLocalHitScore (cm, ins_mx, L);
      }

      /* Safety valve */
      if (total_proposals > 100 * (n_steps + n_burnin)) {
        printf ("WARNING: chain %d: %d proposals for %d accepts, stopping early\n",
                chain, total_proposals, accepted);
        break;
      }
    }

    if (verbose)
      printf ("  chain %d: done. %d collected, %d accepted, %d proposals (%.1f%% accept rate)\n",
              chain, collected, accepted, total_proposals,
              100.0 * accepted / total_proposals);

    cm_mx_Destroy (ins_mx);
    free (dsq);
  }

  *ret_scores = all_scores;
  *ret_N      = sc_idx;
  return eslOK;

ERROR:
  if (all_scores != NULL) free (all_scores);
  ESL_FAIL (eslEMEM, errbuf, "cm_MCMC_tail: memory allocation error");
}

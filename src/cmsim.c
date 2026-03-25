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
      double target_sc = esl_opt_GetReal (go, "--imix");
      double orig_sc = cm_ExpectedParsetreeScore (cm, 0.0);
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
          mid_sc = cm_ExpectedParsetreeScore (cm, mid);
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
        if ((status = initialize_cm (go, cfg, emit_cm, TRUE, errbuf)) != eslOK)
          cm_Fail (errbuf);
      }
    } else if (esl_opt_IsOn (go, "--exp")) {
      if ((status = cm_Clone (cm, errbuf, &emit_cm)) != eslOK)
        cm_Fail (errbuf);
      cm_Exponentiate (emit_cm, esl_opt_GetReal (go, "--exp"));
      if ((status = initialize_cm (go, cfg, emit_cm, TRUE, errbuf)) != eslOK)
        cm_Fail (errbuf);
    }
    if ((status = initialize_cm (go, cfg, cm, TRUE, errbuf)) != eslOK)
      cm_Fail (errbuf);

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

    /* For now, search only with local inside */
    exp_mode = EXP_CM_LI;
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
  float weight;       /* importance weight = 2^(-parsetree_sc) */
  float ilo, ihi;     /* parsetree score range for rejection sampling */
  int do_filter;      /* TRUE if --ilo or --ihi is set */
  int do_isubtr;              /* TRUE if --isubtr: sub-parsetree IS mode */
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

  /* get HMM for generating random seqs, if nec */
  if (esl_opt_GetBoolean (go, "--rhmm")) {
    if ((status = CreateGenomicHMM (cm->abc, errbuf, &ghmm_sA, &ghmm_tAA, &ghmm_eAA, &ghmm_nstates))
        != eslOK)
      cm_Fail ("ERROR unable to make HMM for generating random seqs");
  }

  /* Search sequences and collect score histograms */

  scN = 0;
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
          printf ("SEQ %5d  L: %4d  parsetree_sc: %8.3f  weight: %12.6g  (emitted: %d "
                  "rejected: %d)\n",
                  i, L, parsetree_sc, weight, n_emitted, n_rejected);
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

    /* Search the sequence with CYK or Inside (follows cmcalibrate.c:process_search_workunit()
     * pattern) */
    th = cm_tophits_Create ();
    if (th == NULL)
      ESL_FAIL (eslEMEM, errbuf, "out of memory");

    if (cm->search_opts & CM_SEARCH_INSIDE) {
      if ((status = FastIInsideScan (cm, errbuf, cm->smx, use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB,
                                     dsq, 1, L, cutoff, th, cm->search_opts & CM_SEARCH_NULL3, 0.,
                                     NULL, NULL, NULL, NULL))
          != eslOK)
        cm_Fail (errbuf);
    } else {
      if ((status = FastCYKScan (cm, errbuf, cm->smx, use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB, dsq, 1,
                                 L, cutoff, th, cm->search_opts & CM_SEARCH_NULL3, 0., NULL, NULL,
                                 NULL, NULL))
          != eslOK)
        cm_Fail (errbuf);
    }
    /* overlaps already removed inside FastCYKScan/FastIInsideScan */

    /* do_isubtr: refine IS weight using the Inside score of the sub-region.
     * Parsetree score underestimates proposal density (Inside >= parsetree);
     * using the actual Inside score of the hit rooted at v* corrects this.
     * th->unsrt[h].root == v* and hit overlaps [il..ir] identifies the target hit.
     * If not found (shadowed by a higher-scoring null-flank hit after overlap removal),
     * we keep the parsetree-based weight (best_candidate) as fallback. */
    if (do_isubtr && do_sample && !esl_opt_GetBoolean (go, "--no-weight") && isubtr_best_v != -1) {
      int   found   = FALSE;
      float best_sc = isubtr_best_candidate;  /* fallback: parsetree-based candidate_sc */
      for (h = 0; h < (int) th->N; h++) {
        if (th->unsrt[h].root  == isubtr_best_v         &&
            th->unsrt[h].start <= (int64_t) isubtr_ir   &&
            th->unsrt[h].stop  >= (int64_t) isubtr_il) {
          if (!found || th->unsrt[h].score > best_sc) {
            best_sc = th->unsrt[h].score;
            found   = TRUE;
          }
        }
      }
      weight = (float) pow (2.0, -best_sc);
      if (found) n_isubtr_hit_found++;
      else        n_isubtr_hit_missing++;
    }

    /* accumulate dbsize: actual nt searched (unweighted for both IS and random).
     * The hits/Mb criterion counts data points, not IS-equivalent null sequence. */
    dbsize += (double)L;

    if (th->N > 0) {
      /* collect all hits */
      if (scN == 0) {
        ESL_ALLOC (scA, sizeof (float) * (scN + th->N));
        if (do_sample)
          ESL_ALLOC (wtA, sizeof (float) * (scN + th->N));
      } else {
        ESL_RALLOC (scA, tmp, sizeof (float) * (scN + th->N));
        if (do_sample)
          ESL_RALLOC (wtA, tmp, sizeof (float) * (scN + th->N));
      }
      for (h = 0; h < (int)th->N; h++) {
        scA[(scN + h)] = th->unsrt[h].score;
        if (do_sample)
          wtA[(scN + h)] = weight;
      }
      scN += th->N;
    }

    cm_tophits_Destroy (th);
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

  /* Determine tailp: use hits/Mb criterion (matching cmcalibrate) unless overridden */
  if (tailp <= 0. && dbsize_nt > 0.) {
    int tailn = ExpModeIsLocal (exp_mode) ? esl_opt_GetInteger (go, "--ltailn")
                                          : esl_opt_GetInteger (go, "--gtailn");
    float nhits_to_fit = (float)tailn * (dbsize_nt / 1e6);
    tailp = nhits_to_fit / (float)h->n;
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
  int v, k, l;
  int K = cm->abc->K;

  psi = cm_ExpectedStateOccupancy (cm);

  for (v = 0; v < cm->M; v++) {
    if (psi[v] == 0.)
      continue;

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

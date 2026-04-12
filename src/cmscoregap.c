/* cmscoregap: Compute parsetree / CYK / Inside score distributions under
 *             three different sequence-generation schemes:
 *
 *   1) default (CM-emit): emit from CM, HMM-banded ALIGNMENT for cyk/inside
 *                         (alpha[0][L][L], full-sequence scores)
 *   2) --random:          generate random sequences from GENOMIC HMM (matching
 *                         cmcalibrate), QDB-banded SCANNING for cyk/inside
 *                         (best hit scores, with null3 correction)
 *   3) --imix <target_sc>: emit from null-mixed emit_cm (binary search alpha
 *                          for expected parsetree score = target_sc), score
 *                          under ORIGINAL cm with HMM-banded ALIGNMENT
 *
 * Usage:   cmscoregap [options] <cmfile>
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "infernal.h"

#define CMSCOREGAP_DEFAULT_N       200
#define CMSCOREGAP_DEFAULT_MXSIZE  2048.0

static ESL_OPTIONS options[] = {
  { "-h",          eslARG_NONE,  FALSE,         NULL, NULL,   NULL,   NULL, NULL,        "show brief help on version and usage",                    1 },
  { "-N",          eslARG_INT,   "200",         NULL, "n>0",  NULL,   NULL, NULL,        "number of sequences per (CM, mode)",                      1 },
  { "-o",          eslARG_OUTFILE,NULL,         NULL, NULL,   NULL,   NULL, NULL,        "write TSV output to file <f> instead of stdout",          1 },
  { "--seed",      eslARG_INT,   "42",          NULL, "n>=0", NULL,   NULL, NULL,        "set random number seed to <n>",                           1 },
  { "--local",     eslARG_NONE,  FALSE,         NULL, NULL,   NULL,   NULL, "--glocal",  "only local mode",                                         1 },
  { "--glocal",    eslARG_NONE,  FALSE,         NULL, NULL,   NULL,   NULL, "--local",   "only glocal mode",                                        1 },
  { "--random",    eslARG_NONE,  FALSE,         NULL, NULL,   NULL,   NULL, "--imix",    "random seqs from genomic HMM (QDB scan, like cmcalibrate)", 1 },
  { "--imix",      eslARG_REAL,  NULL,          NULL, NULL,   NULL,   NULL, "--random",  "target expected parsetree score (bits) for null-mixing",   1 },
  { "--rnd-L",     eslARG_INT,   NULL,          NULL, "n>0",  NULL, "--random", NULL,    "length for --random sequences (default: 10*cm->W)",       1 },
  { "--beta",      eslARG_REAL,  "1e-15",       NULL, "x>0",  NULL,   NULL, NULL,        "QDB beta for --random mode scanning",                     1 },
  { "--tau",       eslARG_REAL,  "1e-7",        NULL, "x>0",  NULL,   NULL, NULL,        "HMM band tau for --imix and default modes",               1 },
  { "--mxsize",    eslARG_REAL,  "2048.0",      NULL, "x>0",  NULL,   NULL, NULL,        "DP matrix size limit in Mb",                              1 },
  { "--maxtau",    eslARG_REAL,  "1e-3",        NULL, "x>0",  NULL,   NULL, NULL,        "max tau for band iteration",                              1 },
  { "--nonull3",   eslARG_NONE,  FALSE,         NULL, NULL,   NULL,   NULL, NULL,        "turn OFF null3 score correction for --random scanning",    1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

static char usage[]  = "[options] <cmfile>";
static char banner[] = "score gap distributions under CM-emit / random / null-mix proposals";

/* ---- Copied from cmsim.c (static there) ---- */
static double
cm_ExpectedParsetreeScore(CM_t *cm, double alpha) {
  double *psi = NULL;
  double E = 0.;
  double q;
  int v, k, l, c;
  int K = cm->abc->K;

  psi = cm_ExpectedStateOccupancy(cm);

  for (v = 0; v < cm->M; v++) {
    if (psi[v] == 0.) continue;
    if (cm->sttype[v] == MP_st) {
      for (k = 0; k < K; k++)
        for (l = 0; l < K; l++) {
          q = (1.0 - alpha) * cm->e[v][k * K + l] + alpha * cm->null[k] * cm->null[l];
          if (q > 0.) E += psi[v] * q * log2(q / (cm->null[k] * cm->null[l]));
        }
    } else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
               cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      for (k = 0; k < K; k++) {
        q = (1.0 - alpha) * cm->e[v][k] + alpha * cm->null[k];
        if (q > 0.) E += psi[v] * q * log2(q / cm->null[k]);
      }
    }
    if (cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
      if (v == 0 && (cm->flags & CMH_LOCAL_BEGIN)) {
        int y;
        for (y = 0; y < cm->M; y++)
          if (cm->begin[y] > 0.) E += psi[v] * cm->begin[y] * cm->beginsc[y];
      } else {
        for (c = 0; c < cm->cnum[v]; c++)
          if (cm->t[v][c] > 0.) E += psi[v] * cm->t[v][c] * cm->tsc[v][c];
      }
      if ((cm->flags & CMH_LOCAL_END) && cm->end[v] > 0.)
        E += psi[v] * cm->end[v] * cm->endsc[v];
    }
  }
  free(psi);
  return E;
}

static void
cm_MixWithNull(CM_t *cm, double alpha) {
  int v, k, l;
  int K = cm->abc->K;
  for (v = 0; v < cm->M; v++) {
    if (cm->sttype[v] == MP_st) {
      for (k = 0; k < K; k++)
        for (l = 0; l < K; l++)
          cm->e[v][k * K + l] = (1.0 - alpha) * cm->e[v][k * K + l] + alpha * cm->null[k] * cm->null[l];
    } else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
               cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      for (k = 0; k < K; k++)
        cm->e[v][k] = (1.0 - alpha) * cm->e[v][k] + alpha * cm->null[k];
    }
  }
  cm->flags &= ~CMH_BITS;
}
/* ---- end copied from cmsim.c ---- */


/* Configure a CM for HMM-banded alignment (for cm-emit and imix modes). */
static int
configure_cm_align(CM_t *cm, int do_local, char *errbuf)
{
  if (do_local) {
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  return cm_Configure(cm, errbuf, -1);
}


/* Configure a CM for QDB-banded scanning (for --random mode, matching
 * cmcalibrate's setup). */
static int
configure_cm_scan(CM_t *cm, int do_local, double beta, int do_null3, char *errbuf)
{
  int status;

  if (do_local) {
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  cm->search_opts |= CM_SEARCH_QDB;
  cm->search_opts |= CM_SEARCH_NOALIGN;
  if (do_null3) cm->search_opts |= CM_SEARCH_NULL3;

  if (CheckCMQDBInfo(cm->qdbinfo, 0., FALSE, beta, TRUE) != eslOK) {
    cm->config_opts |= CM_CONFIG_QDB;
    cm->qdbinfo->beta1 = beta;
    cm->qdbinfo->beta2 = beta;
  }
  cm->config_opts |= CM_CONFIG_SCANMX;

  if ((status = cm_Configure(cm, errbuf, -1)) != eslOK) return status;
  return eslOK;
}


/* Score a sequence with HMM-banded CYK and Inside alignment.
 * Returns alpha[0][L][L] full-sequence scores. */
static int
score_align_hb(CM_t *cm, ESL_DSQ *dsq, int L, float mxsize, double maxtau,
               float *ret_cyk_sc, float *ret_ins_sc, char *errbuf)
{
  int   status;
  float cyk_sc, ins_sc;

  status = cp9_IterateSeq2Bands(cm, errbuf, dsq, 1, L,
                                PLI_PASS_STD_ANY, mxsize,
                                FALSE, FALSE, FALSE, TRUE, maxtau, NULL);
  if (status != eslOK) return status;

  status = cm_CYKInsideAlignHB(cm, errbuf, dsq, L, mxsize,
                               cm->hb_mx, cm->hb_shmx, NULL, &cyk_sc);
  if (status != eslOK) return status;

  status = cm_InsideAlignHB(cm, errbuf, dsq, L, mxsize, cm->hb_mx, &ins_sc);
  if (status != eslOK) return status;

  *ret_cyk_sc = cyk_sc;
  *ret_ins_sc = ins_sc;
  return eslOK;
}


/* Score a random sequence with QDB-banded CYK and Inside scanning.
 * Returns best-hit scores from the scan, matching cmcalibrate's scoring.
 * Uses float CYK (FastCYKScan) and float Inside (FastFInsideScan). */
static int
score_scan_qdb(CM_t *cm, ESL_DSQ *dsq, int L, int do_null3,
               float *ret_cyk_sc, float *ret_ins_sc, char *errbuf)
{
  int        status;
  float      cyk_sc = -eslINFINITY;
  float      ins_sc = -eslINFINITY;
  CM_TOPHITS *th    = NULL;
  int        qdbidx = (cm->search_opts & CM_SEARCH_QDB) ? SMX_QDB2_LOOSE : SMX_NOQDB;

  /* CYK scan (CM_SEARCH_INSIDE must be OFF) */
  cm->search_opts &= ~CM_SEARCH_INSIDE;
  th = cm_tophits_Create();
  status = FastCYKScan(cm, errbuf, cm->smx, qdbidx,
                       dsq, 1, L, -eslINFINITY, th,
                       do_null3, 0., NULL, NULL, NULL, &cyk_sc);
  cm_tophits_Destroy(th);
  if (status != eslOK) return status;

  /* Float Inside scan (CM_SEARCH_INSIDE must be ON) */
  cm->search_opts |= CM_SEARCH_INSIDE;
  th = cm_tophits_Create();
  status = FastFInsideScan(cm, errbuf, cm->smx, qdbidx,
                           dsq, 1, L, -eslINFINITY, th,
                           do_null3, 0., NULL, NULL, NULL, &ins_sc);
  cm_tophits_Destroy(th);
  cm->search_opts &= ~CM_SEARCH_INSIDE;  /* restore */
  if (status != eslOK) return status;

  *ret_cyk_sc = cyk_sc;
  *ret_ins_sc = ins_sc;
  return eslOK;
}


/* Build a null-mixed emit_cm for a target expected parsetree score. */
static int
build_imix_emit_cm(CM_t *cm, int do_local, double target_sc,
                   CM_t **ret_emit_cm, double *ret_alpha, double *ret_orig_sc,
                   char *errbuf)
{
  int    status;
  CM_t  *tmp_cm   = NULL;
  CM_t  *emit_cm  = NULL;
  double lo = 0., hi = 1., mid = 0.5, mid_sc = 0., orig_sc;
  int    iter;

  if ((status = cm_Clone(cm, errbuf, &tmp_cm)) != eslOK) goto ERROR;
  if ((status = configure_cm_align(tmp_cm, do_local, errbuf)) != eslOK) goto ERROR;

  orig_sc = cm_ExpectedParsetreeScore(tmp_cm, 0.0);
  if (target_sc >= orig_sc) {
    ESL_FAIL(eslEINVAL, errbuf,
             "--imix target %.2f >= original expected score %.2f",
             target_sc, orig_sc);
  }

  for (iter = 0; iter < 100; iter++) {
    mid = (lo + hi) / 2.0;
    mid_sc = cm_ExpectedParsetreeScore(tmp_cm, mid);
    if (fabs(mid_sc - target_sc) < 0.01) break;
    if (mid_sc > target_sc) lo = mid;
    else                    hi = mid;
  }

  if ((status = cm_Clone(cm, errbuf, &emit_cm)) != eslOK) goto ERROR;
  cm_MixWithNull(emit_cm, mid);
  if ((status = configure_cm_align(emit_cm, do_local, errbuf)) != eslOK) goto ERROR;

  FreeCM(tmp_cm);
  *ret_emit_cm = emit_cm;
  *ret_alpha   = mid;
  *ret_orig_sc = orig_sc;
  return eslOK;

 ERROR:
  if (tmp_cm  != NULL) FreeCM(tmp_cm);
  if (emit_cm != NULL) FreeCM(emit_cm);
  return status;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = NULL;
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_STOPWATCH  *w2      = esl_stopwatch_Create();
  char           *cmfile  = NULL;
  FILE           *ofp     = stdout;
  char            errbuf[eslERRBUFSIZE];
  CM_FILE        *cmfp    = NULL;
  ESL_ALPHABET   *abc     = NULL;
  ESL_RANDOMNESS *rng     = NULL;
  CM_t           *cm      = NULL;
  int             N, seed, do_local_only, do_glocal_only, do_random, do_imix, do_null3;
  int             rnd_L_override = -1;
  float           mxsize;
  double          tau, maxtau, beta, imix_target = 0.0;
  int             status;

  /* Genomic HMM for --random mode (matching cmcalibrate) */
  int     ghmm_nstates = 0;
  double *ghmm_sA      = NULL;
  double **ghmm_tAA    = NULL;
  double **ghmm_eAA    = NULL;

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessEnvironment(go)         != eslOK) cm_Fail("env error");
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cm_Fail("cmdline error: %s", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cm_Fail("config error: %s", go->errbuf);

  if (esl_opt_GetBoolean(go, "-h")) {
    cm_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);
    puts("\nOptions:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
    exit(0);
  }
  if (esl_opt_ArgNumber(go) != 1) {
    printf("Incorrect number of command line arguments.\n");
    esl_usage(stdout, argv[0], usage);
    exit(1);
  }

  cmfile         = esl_opt_GetArg(go, 1);
  N              = esl_opt_GetInteger(go, "-N");
  seed           = esl_opt_GetInteger(go, "--seed");
  do_local_only  = esl_opt_GetBoolean(go, "--local");
  do_glocal_only = esl_opt_GetBoolean(go, "--glocal");
  do_random      = esl_opt_GetBoolean(go, "--random");
  do_imix        = esl_opt_IsOn(go, "--imix");
  do_null3       = !esl_opt_GetBoolean(go, "--nonull3");
  if (do_imix) imix_target = esl_opt_GetReal(go, "--imix");
  if (esl_opt_IsOn(go, "--rnd-L")) rnd_L_override = esl_opt_GetInteger(go, "--rnd-L");
  mxsize         = (float) esl_opt_GetReal(go, "--mxsize");
  tau            = esl_opt_GetReal(go, "--tau");
  maxtau         = esl_opt_GetReal(go, "--maxtau");
  beta           = esl_opt_GetReal(go, "--beta");

  if (esl_opt_IsOn(go, "-o")) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
      cm_Fail("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  }

  rng = esl_randomness_Create(seed);
  FLogsumInit();  /* REQUIRED before any Inside call */

  /* Header */
  fprintf(ofp, "# cmscoregap output\n");
  fprintf(ofp, "# cmfile=%s N=%d seed=%d tau=%g beta=%g mxsize=%.0f null3=%s\n",
          cmfile, N, seed, tau, beta, mxsize, do_null3 ? "on" : "off");
  if (do_random) fprintf(ofp, "# mode=random (genomic HMM sequences, QDB-banded scan)\n");
  else if (do_imix) fprintf(ofp, "# mode=imix target_sc=%.2f\n", imix_target);
  else fprintf(ofp, "# mode=cm-emit\n");
  fprintf(ofp, "cm_name\tacc\tproposal\tmode\talpha\tseq_idx\tL\t"
               "parsetree_sc\tcyk_sc\tinside_sc\ttime_ms\n");

  if ((status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf)) != eslOK)
    cm_Fail("Failed to open CM file %s\n%s", cmfile, errbuf);

  esl_stopwatch_Start(w);

  /* Create genomic HMM once (used for all CMs in --random mode).
   * Need abc from first CM read, so we read it, create the HMM, then rewind.
   * Alternatively, create it after the first cm_file_Read. */
  int ghmm_created = 0;

  while ((status = cm_file_Read(cmfp, TRUE, &abc, &cm)) == eslOK) {
    const char *cm_name = cm->name ? cm->name : "unnamed";
    const char *cm_acc  = cm->acc  ? cm->acc  : "-";

    /* Create genomic HMM on first CM (needs abc) */
    if (do_random && !ghmm_created) {
      if ((status = CreateGenomicHMM(abc, errbuf, &ghmm_sA, &ghmm_tAA, &ghmm_eAA, &ghmm_nstates)) != eslOK)
        cm_Fail("CreateGenomicHMM failed: %s", errbuf);
      ghmm_created = 1;
    }

    int modes[2]; int n_modes = 0;
    if      (do_local_only)  { modes[0] = 1; n_modes = 1; }
    else if (do_glocal_only) { modes[0] = 0; n_modes = 1; }
    else                     { modes[0] = 0; modes[1] = 1; n_modes = 2; }

    for (int mi = 0; mi < n_modes; mi++) {
      int   do_local     = modes[mi];
      const char *mode_label = do_local ? "local" : "glocal";
      const char *proposal   = do_random ? "random" : (do_imix ? "imix" : "cm-emit");
      double alpha_used = 0.0;
      double orig_sc    = 0.0;

      /* Fresh clone for this (CM, mode) */
      CM_t *cm_search = NULL;
      if ((status = cm_Clone(cm, errbuf, &cm_search)) != eslOK)
        cm_Fail("cm_Clone: %s", errbuf);
      cm_search->tau = tau;

      /* Configure: QDB scan for --random, HMM-banded align for emit/imix */
      if (do_random) {
        if ((status = configure_cm_scan(cm_search, do_local, beta, do_null3, errbuf)) != eslOK)
          cm_Fail("configure_cm_scan failed: %s", errbuf);
      } else {
        if ((status = configure_cm_align(cm_search, do_local, errbuf)) != eslOK)
          cm_Fail("configure_cm_align failed: %s", errbuf);
      }

      /* For --imix, build the null-mixed emit_cm */
      CM_t *emit_cm = NULL;
      if (do_imix) {
        if ((status = build_imix_emit_cm(cm, do_local, imix_target,
                                         &emit_cm, &alpha_used, &orig_sc, errbuf)) != eslOK) {
          fprintf(stderr, "# WARNING %s %s: imix setup failed: %s -- skipping\n",
                  cm_name, mode_label, errbuf);
          FreeCM(cm_search);
          continue;
        }
        fprintf(stderr, "# %s %s imix: target=%.2f orig_sc=%.2f alpha=%.4f\n",
                cm_name, mode_label, imix_target, orig_sc, alpha_used);
      }

      int n_scored = 0, n_failed = 0, attempts = 0;
      int max_attempts = N * 20;

      while (n_scored < N && attempts < max_attempts) {
        attempts++;
        Parsetree_t *tr   = NULL;
        ESL_SQ      *sq   = NULL;
        ESL_DSQ     *dsq  = NULL;
        int          L    = 0;
        float        parsetree_sc = 0., cyk_sc = 0., ins_sc = 0.;
        int          have_parsetree = 0;
        double       elapsed_ms;

        esl_stopwatch_Start(w2);

        if (do_random) {
          /* Generate random sequence from genomic HMM (matching cmcalibrate) */
          L = (rnd_L_override > 0) ? rnd_L_override : 10 * cm_search->W;
          if ((status = SampleGenomicSequenceFromHMM(rng, abc, errbuf, ghmm_sA, ghmm_tAA,
                                                     ghmm_eAA, ghmm_nstates, L, &dsq)) != eslOK) {
            fprintf(stderr, "# WARNING SampleGenomicSequenceFromHMM failed: %s\n", errbuf);
            n_failed++;
            continue;
          }
        } else {
          CM_t *emit_from = do_imix ? emit_cm : cm_search;
          if ((status = EmitParsetree(emit_from, errbuf, rng, "seq", TRUE,
                                       &tr, &sq, &L)) != eslOK) {
            fprintf(stderr, "# WARNING EmitParsetree failed: %s\n", errbuf);
            n_failed++;
            continue;
          }
          if (L == 0) {
            esl_sq_Destroy(sq); FreeParsetree(tr);
            continue;
          }
          dsq = sq->dsq;  /* alias */
          /* Parsetree score under the ORIGINAL (search) cm */
          if ((status = ParsetreeScore(cm_search, NULL, errbuf, tr, dsq, FALSE,
                                        &parsetree_sc, NULL, NULL, NULL, NULL)) != eslOK) {
            fprintf(stderr, "# WARNING ParsetreeScore failed seq %d: %s\n", n_scored, errbuf);
            n_failed++;
            esl_sq_Destroy(sq); FreeParsetree(tr);
            continue;
          }
          if (!isfinite(parsetree_sc)) {
            n_failed++;
            esl_sq_Destroy(sq); FreeParsetree(tr);
            continue;
          }
          have_parsetree = 1;
        }

        /* Score CYK and Inside */
        if (do_random) {
          status = score_scan_qdb(cm_search, dsq, L, do_null3,
                                  &cyk_sc, &ins_sc, errbuf);
        } else {
          status = score_align_hb(cm_search, dsq, L, mxsize, maxtau,
                                  &cyk_sc, &ins_sc, errbuf);
        }
        esl_stopwatch_Stop(w2);
        elapsed_ms = w2->elapsed * 1000.0;

        if (status != eslOK) {
          fprintf(stderr, "# WARNING %s %s seq %d: scoring failed: %s\n",
                  cm_name, mode_label, n_scored, errbuf);
          n_failed++;
          if (do_random) { if (dsq) free(dsq); }
          else           { esl_sq_Destroy(sq); FreeParsetree(tr); }
          continue;
        }

        fprintf(ofp, "%s\t%s\t%s\t%s\t%.4f\t%d\t%d\t",
                cm_name, cm_acc, proposal, mode_label, alpha_used, n_scored, L);
        if (have_parsetree) fprintf(ofp, "%.6f\t", parsetree_sc);
        else                fprintf(ofp, "NA\t");
        fprintf(ofp, "%.6f\t%.6f\t%.2f\n", cyk_sc, ins_sc, elapsed_ms);

        n_scored++;
        if (do_random) { free(dsq); }
        else           { esl_sq_Destroy(sq); FreeParsetree(tr); }
        continue;

      ERROR:
        if (do_random && dsq) free(dsq);
        if (sq) esl_sq_Destroy(sq);
        if (tr) FreeParsetree(tr);
        n_failed++;
      }

      fprintf(stderr, "# %s %s %s%s: scored %d/%d (%d failed, %d attempts)\n",
              cm_name, mode_label, proposal,
              do_imix ? " (imix)" : "",
              n_scored, N, n_failed, attempts);

      if (emit_cm) FreeCM(emit_cm);
      FreeCM(cm_search);
    }

    FreeCM(cm);
    cm = NULL;
  }
  if (status != eslEOF) cm_Fail("cm_file_Read failed: %s", cmfp->errbuf);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stderr, w, "# Total CPU: ");

  /* Clean up genomic HMM */
  if (ghmm_created) {
    int i;
    for (i = 0; i < ghmm_nstates; i++) {
      free(ghmm_eAA[i]);
      free(ghmm_tAA[i]);
    }
    free(ghmm_eAA);
    free(ghmm_tAA);
    free(ghmm_sA);
  }

  cm_file_Close(cmfp);
  if (abc != NULL) esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_stopwatch_Destroy(w);
  esl_stopwatch_Destroy(w2);
  if (ofp != stdout) fclose(ofp);
  esl_getopts_Destroy(go);
  return 0;
}

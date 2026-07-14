/* cm_fast_calibrate_regression.c
 * In-tree regression test for cm_FastCalibrate().
 *
 * Usage: ./cm_fast_calibrate_regression <testdir>
 *
 * <testdir> must contain:
 *   - 30 CM files (*.cm) — 15 STR and 15 NOSS
 *   - cm_fast_calibrate_expected.tsv — expected predictions from Python reference
 *
 * For each CM, calls cm_FastCalibrate() and compares predictions to expected.
 * Tolerance: 1e-6 for STR CMs, 1e-3 for NOSS CMs.
 *
 * The larger NOSS tolerance accounts for amplified floating-point rounding:
 * noend_p_full_length has tiny feature_std (~7e-5), which amplifies ~7e-10
 * differences between C sequential summation and Python numpy pairwise
 * summation by ~45000x, yielding max observed error ~2e-4.
 *
 * Returns 0 if all comparisons pass, 1 if any fail.
 *
 * To regenerate expected values after model changes:
 *   python3 scripts/gen_expected_tsv.py <testdir>/cm_fast_calibrate_expected.tsv
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "infernal.h"
#include "cm_fast_calibrate.h"

#define MAX_CMS    64
#define MAX_PATH   512
#define MAX_NAME   256

/* One expected-value record */
typedef struct {
  char   cm_basename[MAX_NAME];
  char   mode[8];
  double lambda;
  double mu_extrap;
  double mu_orig;
} ExpRow_t;

/* Mode name → expA index */
static int
mode_index(const char *mode)
{
  if (strcmp(mode, "ECMGC") == 0) return EXP_CM_GC;
  if (strcmp(mode, "ECMGI") == 0) return EXP_CM_GI;
  if (strcmp(mode, "ECMLC") == 0) return EXP_CM_LC;
  if (strcmp(mode, "ECMLI") == 0) return EXP_CM_LI;
  return -1;
}

int
main(int argc, char **argv)
{
  char          testdir[MAX_PATH];
  char          tsv_path[MAX_PATH];
  char          cm_path[MAX_PATH];
  char          errbuf[eslERRBUFSIZE];
  FILE         *tsvfp  = NULL;
  CM_FILE      *cmfp   = NULL;
  CM_t         *cm     = NULL;
  ESL_ALPHABET *abc    = NULL;
  ExpRow_t     *rows   = NULL;
  int           nrows  = 0;
  int           nrows_alloc = 256;
  int           npass  = 0;
  int           nfail  = 0;
  int           status, hstatus;
  char          line[1024];
  int           i;

  if (argc != 2) {
    fprintf(stderr, "usage: cm_fast_calibrate_regression <testdir>\n");
    fprintf(stderr, "  <testdir> must contain *.cm files and cm_fast_calibrate_expected.tsv\n");
    return 1;
  }
  strncpy(testdir, argv[1], MAX_PATH - 1);
  testdir[MAX_PATH - 1] = '\0';

  /* ---- Load expected TSV ---- */
  snprintf(tsv_path, MAX_PATH, "%s/cm_fast_calibrate_expected.tsv", testdir);
  if ((tsvfp = fopen(tsv_path, "r")) == NULL) {
    fprintf(stderr, "Error: cannot open %s\n", tsv_path);
    return 1;
  }
  rows = malloc(nrows_alloc * sizeof(ExpRow_t));
  if (!rows) { fprintf(stderr, "malloc failed\n"); return 1; }

  while (fgets(line, sizeof(line), tsvfp)) {
    if (line[0] == '#' || line[0] == '\n') continue;
    if (nrows >= nrows_alloc) {
      nrows_alloc *= 2;
      rows = realloc(rows, nrows_alloc * sizeof(ExpRow_t));
      if (!rows) { fprintf(stderr, "realloc failed\n"); return 1; }
    }
    if (sscanf(line, "%255s %7s %lf %lf %lf",
               rows[nrows].cm_basename,
               rows[nrows].mode,
               &rows[nrows].lambda,
               &rows[nrows].mu_extrap,
               &rows[nrows].mu_orig) == 5) {
      nrows++;
    }
  }
  fclose(tsvfp);

  if (nrows == 0) {
    fprintf(stderr, "Error: no data rows read from %s\n", tsv_path);
    free(rows);
    return 1;
  }

  /* ---- Collect CM files from testdir ---- */
  DIR           *dp  = NULL;
  struct dirent *ent = NULL;
  char          cm_names[MAX_CMS][MAX_NAME];
  int           ncm_files = 0;

  dp = opendir(testdir);
  if (!dp) { fprintf(stderr, "Error: cannot opendir %s\n", testdir); free(rows); return 1; }
  while ((ent = readdir(dp)) != NULL && ncm_files < MAX_CMS) {
    int len = strlen(ent->d_name);
    if (len > 3 && strcmp(ent->d_name + len - 3, ".cm") == 0) {
      strncpy(cm_names[ncm_files], ent->d_name, MAX_NAME - 1);
      cm_names[ncm_files][MAX_NAME - 1] = '\0';
      ncm_files++;
    }
  }
  closedir(dp);

  /* Sort CM names for reproducible output */
  for (i = 0; i < ncm_files - 1; i++) {
    int j;
    for (j = i + 1; j < ncm_files; j++) {
      if (strcmp(cm_names[i], cm_names[j]) > 0) {
        char tmp[MAX_NAME];
        strncpy(tmp,        cm_names[i], MAX_NAME);
        strncpy(cm_names[i], cm_names[j], MAX_NAME);
        strncpy(cm_names[j], tmp,         MAX_NAME);
      }
    }
  }

  /* ---- Run calibration and compare ---- */
  for (i = 0; i < ncm_files; i++) {
    const char *basename = cm_names[i];
    snprintf(cm_path, MAX_PATH, "%s/%s", testdir, basename);

    /* Determine tolerance: NOSS CMs have "_noss" suffix */
    int is_noss = (strstr(basename, "_noss") != NULL);
    double tol  = is_noss ? 1e-3 : 1e-6;

    /* Open + read CM */
    if ((status = cm_file_Open(cm_path, NULL, FALSE, &cmfp, errbuf)) != eslOK) {
      fprintf(stderr, "FAIL (open): %s: %s\n", basename, errbuf);
      nfail++; continue;
    }
    if ((hstatus = cm_file_Read(cmfp, TRUE, &abc, &cm)) != eslOK) {
      fprintf(stderr, "FAIL (read): %s: hstatus=%d\n", basename, hstatus);
      cm_file_Close(cmfp); nfail++; continue;
    }
    cm_file_Close(cmfp);

    /* Calibrate */
    if ((status = cm_FastCalibrate(cm)) != eslOK) {
      fprintf(stderr, "FAIL (calibrate): %s: status=%d\n", basename, status);
      FreeCM(cm); cm = NULL; nfail++; continue;
    }

    /* Compare all 4 modes */
    int cm_pass = 1;
    double cm_max_diff = 0.0;
    char   cm_worst[32] = "";
    int    j;
    for (j = 0; j < nrows; j++) {
      if (strcmp(rows[j].cm_basename, basename) != 0) continue;
      int mi = mode_index(rows[j].mode);
      if (mi < 0) continue;
      double c_lam = cm->expA[mi]->lambda;
      double c_mue = cm->expA[mi]->mu_extrap;
      double c_muo = cm->expA[mi]->mu_orig;
      double d_lam = fabs(c_lam - rows[j].lambda);
      double d_mue = fabs(c_mue - rows[j].mu_extrap);
      double d_muo = fabs(c_muo - rows[j].mu_orig);
      if (d_lam > tol) {
        fprintf(stderr, "  MISMATCH %s %s lambda: C=%.12g Py=%.12g diff=%.2e tol=%.0e\n",
                basename, rows[j].mode, c_lam, rows[j].lambda, d_lam, tol);
        cm_pass = 0;
      }
      if (d_mue > tol) {
        fprintf(stderr, "  MISMATCH %s %s mu_extrap: C=%.12g Py=%.12g diff=%.2e tol=%.0e\n",
                basename, rows[j].mode, c_mue, rows[j].mu_extrap, d_mue, tol);
        cm_pass = 0;
      }
      if (d_muo > tol) {
        fprintf(stderr, "  MISMATCH %s %s mu_orig: C=%.12g Py=%.12g diff=%.2e tol=%.0e\n",
                basename, rows[j].mode, c_muo, rows[j].mu_orig, d_muo, tol);
        cm_pass = 0;
      }
      double max3 = d_lam > d_mue ? d_lam : d_mue;
      if (d_muo > max3) max3 = d_muo;
      if (max3 > cm_max_diff) {
        cm_max_diff = max3;
        snprintf(cm_worst, sizeof(cm_worst), "%s.all", rows[j].mode);
      }
    }

    if (cm_pass) {
      printf("PASS %-30s  max_abs=%.2e\n", basename, cm_max_diff);
      npass++;
    } else {
      printf("FAIL %-30s  max_abs=%.2e (tol=%.0e)\n", basename, cm_max_diff, tol);
      nfail++;
    }

    FreeCM(cm);
    cm = NULL;
  }

  cm_FastCalibrateCleanup();
  esl_alphabet_Destroy(abc);
  free(rows);

  printf("\n=== SUMMARY: %d passed, %d failed (of %d CMs) ===\n", npass, nfail, ncm_files);
  if (nfail > 0) {
    fprintf(stderr, "REGRESSION FAILED: %d/%d CMs mismatched beyond tolerance.\n", nfail, ncm_files);
    fprintf(stderr, "To regenerate expected values after model changes:\n");
    fprintf(stderr, "  python3 <proj>/scripts/gen_expected_tsv.py %s/cm_fast_calibrate_expected.tsv\n", testdir);
  }
  return nfail > 0 ? 1 : 0;
}

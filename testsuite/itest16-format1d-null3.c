/* itest16-format1d-null3.c
 *
 * Round-trip harness for the unified INFERNAL1/d CM file format
 * (null3-OFF E-value block + consensus pknot annotation union).
 * Companion to itest13-pknot.pl (pknot side) and the null3 store-both
 * format port (brief 080).
 *
 * Since the null3 *producer* feature code (cm_FastCalibrate off-set,
 * cmcalibrate --nonull3) is NOT part of the format port, this harness
 * synthesizes a cm->expA_nonull3[] with known values, exercises the
 * WriteASCII/WriteBinary + Read paths, and asserts the block + flag +
 * format survive round-trips. It covers all four content combinations
 * (pknot-only, null3-only, both, neither) through BOTH ASCII and binary,
 * plus 1/c rejection and 1/a,1/b backward read.
 *
 * Compile (against the already-built worktree libs; run from worktree root):
 *   gcc -std=gnu99 -O2 -o itest16 \
 *       -I src -I easel -I hmmer/src \
 *       testsuite/itest16-format1d-null3.c \
 *       src/libinfernal.a hmmer/src/libhmmer.a easel/libeasel.a -lm -lpthread
 *
 * Usage:    ./itest16 <pknot CM file (CMH_PKNOT set)> <tmp prefix>
 * Prints per-check "ok N: ..." lines; prints "ok\n" and exits 0 on
 * success; prints "FAIL: ..." and exits nonzero otherwise.
 */
#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "hmmer.h"
#include "infernal.h"

static int nfail = 0;

#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s\n", (msg)); nfail++; } } while (0)

/* known, "clean" per-mode null3 values (chosen to survive %.5f/%.6f ASCII
 * formatting exactly for the integer fields, within 1e-9 for the doubles). */
static void
fill_nonull3(CM_t *cm)
{
  int z;
  int status;
  if (cm->expA_nonull3 == NULL) {
    ESL_ALLOC(cm->expA_nonull3, sizeof(ExpInfo_t *) * EXP_NMODES);
    for (z = 0; z < EXP_NMODES; z++) cm->expA_nonull3[z] = CreateExpInfo();
  }
  for (z = 0; z < EXP_NMODES; z++) {
    cm->expA_nonull3[z]->lambda    = 0.40000 + 0.01 * z;
    cm->expA_nonull3[z]->mu_extrap = 10.00000 + (double) z;
    cm->expA_nonull3[z]->mu_orig   = 20.00000 + (double) z;
    cm->expA_nonull3[z]->dbsize    = 1000000.0 * (z + 1);
    cm->expA_nonull3[z]->nrandhits = 1000 * (z + 1);
    cm->expA_nonull3[z]->tailp     = 0.010000 * (z + 1);
    cm->expA_nonull3[z]->is_valid  = TRUE;
  }
  cm->flags |= CMH_EXPTAIL_NONULL3_STATS;
  return;
 ERROR:
  cm_Fail("alloc failure in fill_nonull3");
}

/* assert a read-back CM's expA_nonull3 matches fill_nonull3, given a
 * per-double tolerance (0.0 == exact, for binary). */
static void
check_nonull3(CM_t *cm, double tol, const char *where)
{
  int z;
  char buf[256];
  if (! (cm->flags & CMH_EXPTAIL_NONULL3_STATS)) { snprintf(buf,sizeof(buf),"%s: CMH_EXPTAIL_NONULL3_STATS not set on read-back", where); CHECK(0, buf); return; }
  if (cm->expA_nonull3 == NULL)                  { snprintf(buf,sizeof(buf),"%s: expA_nonull3 NULL on read-back", where);              CHECK(0, buf); return; }
  for (z = 0; z < EXP_NMODES; z++) {
    double e_lambda = 0.40000 + 0.01 * z, e_mue = 10.0 + z, e_muo = 20.0 + z, e_tailp = 0.010000 * (z + 1);
    long   e_db = 1000000L * (z + 1); int e_nr = 1000 * (z + 1);
    if (fabs(cm->expA_nonull3[z]->lambda    - e_lambda) > (tol > 0 ? tol : 0.0)) { snprintf(buf,sizeof(buf),"%s: mode %d lambda %.10f != %.10f", where, z, cm->expA_nonull3[z]->lambda, e_lambda); CHECK(0, buf); }
    if (fabs(cm->expA_nonull3[z]->mu_extrap - e_mue)    > (tol > 0 ? tol : 0.0)) { snprintf(buf,sizeof(buf),"%s: mode %d mu_extrap mismatch", where, z); CHECK(0, buf); }
    if (fabs(cm->expA_nonull3[z]->mu_orig   - e_muo)    > (tol > 0 ? tol : 0.0)) { snprintf(buf,sizeof(buf),"%s: mode %d mu_orig mismatch", where, z); CHECK(0, buf); }
    if (fabs(cm->expA_nonull3[z]->tailp     - e_tailp)  > (tol > 0 ? tol : 1e-9)) { snprintf(buf,sizeof(buf),"%s: mode %d tailp mismatch", where, z); CHECK(0, buf); }
    if ((long) (cm->expA_nonull3[z]->dbsize + 0.5) != e_db) { snprintf(buf,sizeof(buf),"%s: mode %d dbsize %ld != %ld", where, z, (long)(cm->expA_nonull3[z]->dbsize+0.5), e_db); CHECK(0, buf); }
    if (cm->expA_nonull3[z]->nrandhits != e_nr)            { snprintf(buf,sizeof(buf),"%s: mode %d nrandhits %d != %d", where, z, cm->expA_nonull3[z]->nrandhits, e_nr); CHECK(0, buf); }
  }
}

/* write CM via ASCII (fmt=-1 default select), return detected format on
 * re-read; opt sets *ret_cm to the re-read CM (caller frees). */
static int
roundtrip_ascii(CM_t *cm, const char *path, char *errbuf, ESL_ALPHABET **abc, CM_t **ret_cm, int *ret_fmt)
{
  FILE *fp; CM_FILE *cmfp; int status;
  if ((fp = fopen(path, "w")) == NULL) return eslFAIL;
  cm_file_WriteASCII(fp, -1, cm);
  fclose(fp);
  if ((status = cm_file_Open((char*)path, NULL, FALSE, &cmfp, errbuf)) != eslOK) return status;
  *ret_fmt = cmfp->format;
  status = cm_file_Read(cmfp, TRUE, abc, ret_cm);
  if (status != eslOK) strcpy(errbuf, cmfp->errbuf);
  cm_file_Close(cmfp);
  return status;
}

static int
roundtrip_bin(CM_t *cm, const char *path, char *errbuf, ESL_ALPHABET **abc, CM_t **ret_cm, int *ret_fmt)
{
  FILE *fp; CM_FILE *cmfp; int status;
  if ((fp = fopen(path, "wb")) == NULL) return eslFAIL;
  cm_file_WriteBinary(fp, -1, cm, NULL);
  fclose(fp);
  if ((status = cm_file_Open((char*)path, NULL, FALSE, &cmfp, errbuf)) != eslOK) return status;
  *ret_fmt = cmfp->format;
  status = cm_file_Read(cmfp, TRUE, abc, ret_cm);
  if (status != eslOK) strcpy(errbuf, cmfp->errbuf);
  cm_file_Close(cmfp);
  return status;
}

static char *first_line(const char *path, char *buf, int n)
{
  FILE *fp = fopen(path, "r");
  if (fp == NULL) { buf[0] = '\0'; return buf; }
  if (fgets(buf, n, fp) == NULL) buf[0] = '\0';
  fclose(fp);
  return buf;
}

int
main(int argc, char **argv)
{
  char          *cmfile, *pfx;
  char           path[1024], errbuf[eslERRBUFSIZE], line[256];
  ESL_ALPHABET  *abc = NULL;
  CM_FILE       *cmfp = NULL;
  CM_t          *cm = NULL, *cm2 = NULL;
  int            status, fmt;
  int            had_pknot;

  if (argc != 3) { fprintf(stderr, "Usage: %s <pknot CM file> <tmp prefix>\n", argv[0]); return 1; }
  cmfile = argv[1]; pfx = argv[2];

  /* load the pknot CM (CMH_PKNOT set) */
  if ((status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf)) != eslOK) { fprintf(stderr, "FAIL: open %s: %s\n", cmfile, errbuf); return 1; }
  if ((status = cm_file_Read(cmfp, TRUE, &abc, &cm)) != eslOK) { fprintf(stderr, "FAIL: read: %s\n", cmfp->errbuf); cm_file_Close(cmfp); return 1; }
  cm_file_Close(cmfp);
  had_pknot = (cm->flags & CMH_PKNOT) ? 1 : 0;
  CHECK(had_pknot, "input CM must carry CMH_PKNOT");

  /* ---------- CASE: BOTH (pknot + null3), ASCII ---------- */
  fill_nonull3(cm);   /* cm now has pknot + null3 */
  snprintf(path, sizeof(path), "%s.both.a.cm", pfx);
  status = roundtrip_ascii(cm, path, errbuf, &abc, &cm2, &fmt);
  CHECK(status == eslOK, "both/ascii: round-trip read failed");
  if (status == eslOK) {
    CHECK(fmt == CM_FILE_1d, "both/ascii: detected format != 1d");
    CHECK(strncmp(first_line(path,line,sizeof(line)), "INFERNAL1/d", 11) == 0, "both/ascii: header not INFERNAL1/d");
    CHECK(cm2->flags & CMH_PKNOT, "both/ascii: pknot flag lost");
    CHECK(cm2->pknot != NULL, "both/ascii: pknot string lost");
    check_nonull3(cm2, 1e-9, "both/ascii");
    printf("ok 1: both (pknot+null3) round-trips through ASCII as 1/d, both blocks survive\n");
    FreeCM(cm2); cm2 = NULL;
  }

  /* ---------- CASE: BOTH, binary (exact) ---------- */
  snprintf(path, sizeof(path), "%s.both.b.cm", pfx);
  status = roundtrip_bin(cm, path, errbuf, &abc, &cm2, &fmt);
  CHECK(status == eslOK, "both/bin: round-trip read failed");
  if (status == eslOK) {
    CHECK(fmt == CM_FILE_1d, "both/bin: detected format != 1d");
    CHECK(cm2->flags & CMH_PKNOT, "both/bin: pknot flag lost");
    CHECK(cm2->pknot != NULL, "both/bin: pknot string lost");
    check_nonull3(cm2, 0.0, "both/bin");  /* binary is bit-exact */
    printf("ok 2: both (pknot+null3) round-trips through binary as 1/d, values bit-exact\n");
    FreeCM(cm2); cm2 = NULL;
  }

  /* ---------- CASE: NULL3-ONLY (drop pknot) ---------- */
  cm->flags &= ~CMH_PKNOT;
  if (cm->pknot != NULL) { free(cm->pknot); cm->pknot = NULL; }
  /* ASCII */
  snprintf(path, sizeof(path), "%s.n3.a.cm", pfx);
  status = roundtrip_ascii(cm, path, errbuf, &abc, &cm2, &fmt);
  CHECK(status == eslOK, "null3/ascii: round-trip read failed");
  if (status == eslOK) {
    CHECK(fmt == CM_FILE_1d, "null3/ascii: detected format != 1d");
    CHECK(strncmp(first_line(path,line,sizeof(line)), "INFERNAL1/d", 11) == 0, "null3/ascii: header not INFERNAL1/d");
    CHECK(!(cm2->flags & CMH_PKNOT), "null3/ascii: pknot flag unexpectedly set");
    check_nonull3(cm2, 1e-9, "null3/ascii");
    printf("ok 3: null3-only round-trips through ASCII as 1/d (placeholder pknot fields, no PKNOT flag)\n");
    FreeCM(cm2); cm2 = NULL;
  }
  /* binary */
  snprintf(path, sizeof(path), "%s.n3.b.cm", pfx);
  status = roundtrip_bin(cm, path, errbuf, &abc, &cm2, &fmt);
  CHECK(status == eslOK, "null3/bin: round-trip read failed");
  if (status == eslOK) {
    CHECK(fmt == CM_FILE_1d, "null3/bin: detected format != 1d");
    CHECK(!(cm2->flags & CMH_PKNOT), "null3/bin: pknot flag unexpectedly set");
    check_nonull3(cm2, 0.0, "null3/bin");
    printf("ok 4: null3-only round-trips through binary as 1/d, values bit-exact\n");
    FreeCM(cm2); cm2 = NULL;
  }

  /* ---------- CASE: NEITHER (drop null3 too) ---------- */
  cm->flags &= ~CMH_EXPTAIL_NONULL3_STATS;
  /* ASCII: must fall back to 1/b */
  snprintf(path, sizeof(path), "%s.none.a.cm", pfx);
  status = roundtrip_ascii(cm, path, errbuf, &abc, &cm2, &fmt);
  CHECK(status == eslOK, "neither/ascii: round-trip read failed");
  if (status == eslOK) {
    CHECK(fmt == CM_FILE_1b, "neither/ascii: detected format != 1b (should NOT bump to 1d)");
    CHECK(strncmp(first_line(path,line,sizeof(line)), "INFERNAL1/b", 11) == 0, "neither/ascii: header not INFERNAL1/b");
    CHECK(!(cm2->flags & CMH_EXPTAIL_NONULL3_STATS), "neither/ascii: null3 flag unexpectedly set");
    CHECK(!(cm2->flags & CMH_PKNOT), "neither/ascii: pknot flag unexpectedly set");
    printf("ok 5: plain CM (neither block) writes as 1/b, not 1/d\n");
    FreeCM(cm2); cm2 = NULL;
  }
  /* binary: must fall back to 1/b */
  snprintf(path, sizeof(path), "%s.none.b.cm", pfx);
  status = roundtrip_bin(cm, path, errbuf, &abc, &cm2, &fmt);
  CHECK(status == eslOK, "neither/bin: round-trip read failed");
  if (status == eslOK) {
    CHECK(fmt == CM_FILE_1b, "neither/bin: detected format != 1b");
    CHECK(!(cm2->flags & CMH_EXPTAIL_NONULL3_STATS), "neither/bin: null3 flag unexpectedly set");
    printf("ok 6: plain CM (neither block) writes binary as 1/b, not 1/d\n");
    FreeCM(cm2); cm2 = NULL;
  }

  /* ---------- CASE: 1/a and 1/b backward read (explicit format write) ---------- */
  {
    FILE *fp;
    snprintf(path, sizeof(path), "%s.explicit1b.cm", pfx);
    fp = fopen(path, "w"); cm_file_WriteASCII(fp, CM_FILE_1b, cm); fclose(fp);
    status = cm_file_Open(path, NULL, FALSE, &cmfp, errbuf);
    CHECK(status == eslOK, "backward/1b: open failed");
    if (status == eslOK) { status = cm_file_Read(cmfp, TRUE, &abc, &cm2); CHECK(status==eslOK && cmfp->format==CM_FILE_1b, "backward/1b: read failed or wrong format"); cm_file_Close(cmfp); if(cm2){FreeCM(cm2);cm2=NULL;} }

    snprintf(path, sizeof(path), "%s.explicit1a.cm", pfx);
    fp = fopen(path, "w"); cm_file_WriteASCII(fp, CM_FILE_1a, cm); fclose(fp);
    status = cm_file_Open(path, NULL, FALSE, &cmfp, errbuf);
    CHECK(status == eslOK, "backward/1a: open failed");
    if (status == eslOK) { status = cm_file_Read(cmfp, TRUE, &abc, &cm2); CHECK(status==eslOK && cmfp->format==CM_FILE_1a, "backward/1a: read failed or wrong format"); cm_file_Close(cmfp); if(cm2){FreeCM(cm2);cm2=NULL;} }
    printf("ok 7: explicit 1/a and 1/b writes still read back correctly (backward compat)\n");
  }

  /* ---------- CASE: 1/c rejection (ASCII) ---------- */
  {
    FILE *fp;
    snprintf(path, sizeof(path), "%s.reject1c.cm", pfx);
    /* write a real 1/b file then rewrite its header line to INFERNAL1/c */
    fp = fopen(path, "w"); cm_file_WriteASCII(fp, CM_FILE_1b, cm); fclose(fp);
    {
      /* slurp, replace the first "INFERNAL1/b" with "INFERNAL1/c" */
      FILE *r = fopen(path, "r"); char *all; long sz;
      fseek(r, 0, SEEK_END); sz = ftell(r); fseek(r, 0, SEEK_SET);
      all = malloc(sz + 1); if (fread(all, 1, sz, r) != (size_t)sz) { /* ignore */ } all[sz] = '\0'; fclose(r);
      char *p = strstr(all, "INFERNAL1/b"); if (p) p[10] = 'c';
      fp = fopen(path, "w"); fwrite(all, 1, sz, fp); fclose(fp); free(all);
    }
    status = cm_file_Open(path, NULL, FALSE, &cmfp, errbuf);
    /* Open may succeed (format detect) but flag as rejected; or fail at open. */
    if (status == eslOK) {
      status = cm_file_Read(cmfp, TRUE, &abc, &cm2);
      cm_file_Close(cmfp);
      CHECK(status != eslOK, "reject/1c-ascii: an INFERNAL1/c ASCII file was NOT rejected");
      if (cm2) { FreeCM(cm2); cm2 = NULL; }
    } else {
      CHECK(status == eslEFORMAT, "reject/1c-ascii: open failed with unexpected status");
      if (cmfp) cm_file_Close(cmfp);
    }
    printf("ok 8: INFERNAL1/c ASCII header is rejected (retired ambiguous dev format)\n");
  }

  /* ---------- CASE: 1/c rejection (binary magic 0xe3edb0b4) ---------- */
  {
    FILE *fp; uint32_t bad = 0xe3edb0b4;
    snprintf(path, sizeof(path), "%s.reject1c.bcm", pfx);
    fp = fopen(path, "wb"); fwrite(&bad, sizeof(uint32_t), 1, fp);
    /* pad with some bytes so it's not empty */
    { int z=0; for (z=0; z<64; z++) fputc(0, fp); }
    fclose(fp);
    status = cm_file_Open(path, NULL, FALSE, &cmfp, errbuf);
    CHECK(status != eslOK, "reject/1c-bin: a 0xe3edb0b4 binary magic was NOT rejected at open");
    if (status == eslOK && cmfp) cm_file_Close(cmfp);
    printf("ok 9: binary magic 0xe3edb0b4 (INFERNAL1/c) is rejected at open\n");
  }

  if (cm  != NULL) FreeCM(cm);
  if (abc != NULL) esl_alphabet_Destroy(abc);

  if (nfail == 0) { printf("ok\n"); return 0; }
  printf("FAIL: %d check(s) failed\n", nfail);
  return 1;
}

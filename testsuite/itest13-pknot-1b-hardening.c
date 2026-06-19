/* itest13-pknot-1b-hardening.c
 *
 * Regression harness for the binary pknot-block format gate
 * (cm_file.c:879 write, :2421 read). Companion to itest13-pknot.pl.
 *
 * Guards the [CORRECTNESS-RISK] finding from review summary 007: the
 * binary pknot block must be FORMAT-AND-flag gated (>= CM_FILE_1c), not
 * flag-gated alone. If it were flag-gated only, a CM that still carries
 * CMH_PKNOT but is written at an older format (1a/1b) would emit a
 * clen+2-byte pknot block that an older reader misreads as cm->map ->
 * silent corruption. cmconvert clears the flag before writing 1b, but
 * this test proves the GATE -- not the caller's flag-clear -- is what
 * protects older formats: we write a flag-carrying CM directly via
 * cm_file_WriteBinary() at 1b, WITHOUT any flag-clear.
 *
 * Compile (against the already-built worktree libs; run from worktree root):
 *   gcc -std=gnu99 -O2 -o itest13-1b-hardening \
 *       -I src -I easel -I hmmer/src \
 *       testsuite/itest13-pknot-1b-hardening.c \
 *       src/libinfernal.a hmmer/src/libhmmer.a easel/libeasel.a -lm -lpthread
 *
 * Usage:    ./itest13-1b-hardening <1c CM file with CMH_PKNOT set> <tmp prefix>
 * Example:  ./itest13-1b-hardening pkhav.1c.cm tmp13
 *
 * Prints "ok\n" and exits 0 on success; prints "FAIL: ..." and exits
 * nonzero otherwise.
 */
#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "hmmer.h"
#include "infernal.h"

static long
fsize(const char *path)
{
  FILE *fp = fopen(path, "rb");
  long  n;
  if (fp == NULL) return -1;
  if (fseek(fp, 0L, SEEK_END) != 0) { fclose(fp); return -1; }
  n = ftell(fp);
  fclose(fp);
  return n;
}

static int
write_bin(CM_t *cm, int format, const char *path)
{
  FILE *fp = fopen(path, "wb");
  int   status;
  if (fp == NULL) return eslFAIL;
  status = cm_file_WriteBinary(fp, format, cm, NULL);
  fclose(fp);
  return status;
}

int
main(int argc, char **argv)
{
  char          *cmfile;
  char          *pfx;
  char           f_1b_set[1024], f_1b_clr[1024], f_1c_set[1024];
  char           errbuf[eslERRBUFSIZE];
  ESL_ALPHABET  *abc  = NULL;
  CM_FILE       *cmfp = NULL;
  CM_t          *cm   = NULL;
  long           s_1b_set, s_1b_clr, s_1c_set, pkbytes;
  int            status;
  int            nfail = 0;

  if (argc != 3) { fprintf(stderr, "Usage: %s <1c CM file> <tmp prefix>\n", argv[0]); return 1; }
  cmfile = argv[1];
  pfx    = argv[2];
  snprintf(f_1b_set, sizeof(f_1b_set), "%s.1b_flagset.cm",   pfx);
  snprintf(f_1b_clr, sizeof(f_1b_clr), "%s.1b_flagclear.cm", pfx);
  snprintf(f_1c_set, sizeof(f_1c_set), "%s.1c_flagset.cm",   pfx);

  /* Read the input 1c CM (must carry CMH_PKNOT). */
  status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf);
  if (status != eslOK) { fprintf(stderr, "FAIL: cannot open %s: %s\n", cmfile, errbuf); return 1; }
  status = cm_file_Read(cmfp, TRUE, &abc, &cm);
  if (status != eslOK) { fprintf(stderr, "FAIL: cannot read CM: %s\n", cmfp->errbuf); cm_file_Close(cmfp); return 1; }
  cm_file_Close(cmfp);

  if (! (cm->flags & CMH_PKNOT)) { fprintf(stderr, "FAIL: input CM lacks CMH_PKNOT; need a 1c CM with pseudoknots\n"); goto DONE; }
  if (cm->pknot == NULL)         { fprintf(stderr, "FAIL: CMH_PKNOT set but cm->pknot is NULL\n");                     goto DONE; }
  pkbytes = (long) cm->clen + 2;

  /* (a) write 1b binary WITH the flag still set (no cmconvert flag-clear) */
  if (write_bin(cm, CM_FILE_1b, f_1b_set) != eslOK) { fprintf(stderr, "FAIL: WriteBinary 1b (flag set)\n");  goto DONE; }
  /* (b) write 1c binary WITH the flag set (the block SHOULD be present here) */
  if (write_bin(cm, CM_FILE_1c, f_1c_set) != eslOK) { fprintf(stderr, "FAIL: WriteBinary 1c (flag set)\n");  goto DONE; }
  /* (c) clear the flag, write 1b again (reference: definitely no block) */
  cm->flags &= ~CMH_PKNOT;
  if (write_bin(cm, CM_FILE_1b, f_1b_clr) != eslOK) { fprintf(stderr, "FAIL: WriteBinary 1b (flag clear)\n"); goto DONE; }

  s_1b_set = fsize(f_1b_set);
  s_1b_clr = fsize(f_1b_clr);
  s_1c_set = fsize(f_1c_set);
  if (s_1b_set < 0 || s_1b_clr < 0 || s_1c_set < 0) { fprintf(stderr, "FAIL: could not stat output files\n"); goto DONE; }

  printf("# clen=%d  pknot_block_bytes(clen+2)=%ld\n", cm->clen, pkbytes);
  printf("# 1b_flagset=%ld  1b_flagclear=%ld  1c_flagset=%ld\n", s_1b_set, s_1b_clr, s_1c_set);

  /* CHECK 1 (the Fix-1 guarantee): a 1b write with the flag still set emits
   * NO pknot block, so its size equals the flag-cleared 1b write. (The flags
   * word is the same 4 bytes either way.) If the block leaked into 1b, the
   * flag-set file would be clen+2 bytes larger. */
  if (s_1b_set != s_1b_clr) {
    printf("FAIL: 1b output with flag SET (%ld) != with flag CLEAR (%ld); a %ld-byte pknot block leaked into 1b\n",
           s_1b_set, s_1b_clr, s_1b_set - s_1b_clr);
    nfail++;
  } else {
    printf("ok 1: 1b output size is independent of CMH_PKNOT -> no pknot block at 1b (format gate is load-bearing)\n");
  }

  /* CHECK 2: the block IS emitted at 1c -- the 1c file is exactly clen+2 bytes
   * larger than the 1b file (binary format has no other >=1c-gated field). */
  if (s_1c_set - s_1b_set != pkbytes) {
    printf("FAIL: 1c-minus-1b size delta = %ld, expected %ld (the pknot block)\n", s_1c_set - s_1b_set, pkbytes);
    nfail++;
  } else {
    printf("ok 2: 1c output is exactly clen+2 (%ld) bytes larger than 1b -> pknot block present only at 1c\n", pkbytes);
  }

  /* CHECK 3: the flag-set 1b file re-reads cleanly as 1b (read-side gate skips
   * the absent block; map/W stay aligned -> no corruption, no read error). */
  {
    ESL_ALPHABET *abc2  = NULL;
    CM_FILE      *cmfp2 = NULL;
    CM_t         *cm2   = NULL;
    int           clen0 = cm->clen;
    status = cm_file_Open(f_1b_set, NULL, FALSE, &cmfp2, errbuf);
    if (status != eslOK) { printf("FAIL: cannot reopen %s: %s\n", f_1b_set, errbuf); nfail++; }
    else {
      status = cm_file_Read(cmfp2, TRUE, &abc2, &cm2);
      if (status != eslOK)            { printf("FAIL: 1b flag-set file did not re-read cleanly: %s\n", cmfp2->errbuf); nfail++; }
      else if (cm2->clen != clen0)    { printf("FAIL: re-read clen %d != original %d (misalignment)\n", cm2->clen, clen0); nfail++; }
      else                            { printf("ok 3: flag-set 1b file re-reads cleanly as 1b, clen intact (no misalignment)\n"); }
      if (cm2  != NULL) FreeCM(cm2);
      if (abc2 != NULL) esl_alphabet_Destroy(abc2);
    }
    if (cmfp2 != NULL) cm_file_Close(cmfp2);
  }

  remove(f_1b_set); remove(f_1b_clr); remove(f_1c_set);

 DONE:
  if (cm  != NULL) FreeCM(cm);
  if (abc != NULL) esl_alphabet_Destroy(abc);
  if (nfail == 0) { printf("ok\n"); return 0; }
  return 1;
}
